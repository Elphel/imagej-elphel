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
//import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicInteger;


public class TileProcessor {
	public ArrayList <CLTPass3d> clt_3d_passes = null;
	public int       clt_3d_passes_size = 0; //clt_3d_passes size after initial processing
	private int      tilesX;
	private int      tilesY;
	private double   corr_magic_scale =    0.85;  // reported correlation offset vs. actual one (not yet understood)
	private double   trustedCorrelation =  4.0;   // trusted measured disparity difference (before scaling)
	private double   maxOverexposure =     0.5;
	private int      tileSize =            8; // number of linear pixels in a tile (tile is  square tileSize*tileSize)
	int               superTileSize =      8; // number of linear tiles in a super-tile (supertile is  square superTileSize*superTileSize tiles
	                                        // or (superTileSize*tileSize) * (superTileSize*tileSize) pixels, currently 64x64 pixels)
	public int       threadsMax =       100; // maximal number of frames to run
	public int       globalDebugLevel = 0;

	public double [][] dbg_filtered_disp_strength;

	// All parameters are set only once, during instantiation

	public TileProcessor(
			int tilesX,
			int tilesY,
			int tileSize,
			int superTileSize,
			double scale,
			double trustedCorrelation,
			double maxOverexposure,
			int threadsMax)
	{
		this.tilesX = tilesX;
		this.tilesY = tilesY;
		this.tileSize = tileSize;
		this.superTileSize = superTileSize;
		this.corr_magic_scale = scale;
		this.trustedCorrelation = trustedCorrelation;
		this.maxOverexposure =      maxOverexposure;
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
	public void setMaxOverexposure(double maxOverexposure)
	{
		this.maxOverexposure = maxOverexposure;
	}
	public double getMaxOverexposure()
	{
		return this.maxOverexposure;
	}
	public double getMagicScale()
	{
		return this.corr_magic_scale;
	}

	public int getThreadsMax(){
		return this.threadsMax;
	}
//	public void setMagicScale (double scale)
//	{
//		this.corr_magic_scale = scale;
//	}

	public void resetCLTPasses(){
		clt_3d_passes = new ArrayList<CLTPass3d>();
		clt_3d_passes_size = 0;
		Runtime runtime = Runtime.getRuntime();
	    runtime.gc();
		System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
	}

	public void saveCLTPasses(){
		clt_3d_passes_size = clt_3d_passes.size();
	}

	public void trimCLTPasses(int keep){
		clt_3d_passes_size = keep;
		trimCLTPasses();
	}


	public void trimCLTPasses(){
		while (clt_3d_passes.size() > clt_3d_passes_size){
			clt_3d_passes.remove(clt_3d_passes_size);
		}
		Runtime runtime = Runtime.getRuntime();
	    runtime.gc();
		System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");
	}

	public void removeNonMeasurement(){
		for (int i = clt_3d_passes.size()-1; i > 0; i--) if (!clt_3d_passes.get(i).isMeasured()){
			clt_3d_passes.remove(i);
		}
	}





	/**
	 * Basic combining: find smallest residual disparity and use the tile data from it
	 * Copy link to texture tile from the same pass, "forced" bit in tile_op is copied too
	 * Even when this method compares calculated values, it still only copies raw ones, all derivatives should
	 * be re-calculated for the new combined pass
	 *
	 * Calculates max_tried_disparity that shows maximal tried disparity for each tile, regardless of the resulsts/strength
	 * @param passes list of passes to merge
	 * @param firstPass first index in the list to use
	 * @param lastPass last index in the list to use
	 * @param debugLevel debug level
	 * @return combined pass, contains same data as after the measurement of the actual one
	 */
	public CLTPass3d combinePasses_old(
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
	 * When expanding over previously identified background (may be in error) remove tiles from the
	 * composite scan (made by compositeScan()) that have small disparity (those should be identified
	 * on the first passes during background scan) or are not extra strong
	 * @param pass
	 * @param minStrength     full correlation strength to consider data to be extra reliable (after conditioning)
	 * @param minStrengthHor  horizontal (for vertical features) correlation strength to consider data to be extra reliable (after conditioning)
	 * @param minStrengthVert vertical (for horizontal features) correlation strength to consider data to be extra reliable (after conditioning)
	 * @param bg_tiles
	 * @param ex_min_over
	 * @param filtered_disp_strength disparity/strength determined by a floating plate (now 5x5 tiles). [1]Strength < 0 means that last measurement
	 *        does not include this tile, 0.0 - to weak/too bad data, >0.0 - strength. If the whole
	 */
	public void filterOverBackground(
			final CLTPass3d              pass,
			 final double                minStrength,
			 final double                minStrengthHor,
			 final double                minStrengthVert,
			 final boolean []            bg_tiles,          // get from selected in clt_3d_passes.get(0);
			 final double                ex_min_over,        // when expanding over previously detected (by error) background, disregard far tiles
			 final double [][]           filtered_disp_strength
			){
		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++) if (pass.tile_op[ty][tx] != 0) {
				int nt = ty * tilesX + tx;
				if (bg_tiles[nt]){
					boolean [] good_comp = {true,true,true};
					if ((filtered_disp_strength != null) && (filtered_disp_strength[1][nt] == 0.0)){ // -1 - not measured last time, 0.0 - too weak last time
						good_comp = null;
					} else {
						if ((pass.calc_disparity[nt]      < ex_min_over) || (pass.strength[nt])      <= minStrength) good_comp[0] = false;
						if ((pass.calc_disparity_hor[nt]  < ex_min_over) || (pass.strength_hor[nt])  <= minStrengthHor) good_comp[1] = false;
						if ((pass.calc_disparity_vert[nt] < ex_min_over) || (pass.strength_vert[nt]) <= minStrengthVert) good_comp[2] = false;
					}
					if ((good_comp == null) || (!good_comp[0] && !good_comp[1] && !good_comp[2])){
						pass.tile_op[ty][tx] =    0;
						pass.texture_tiles = null;
						for (int i = 0; i< ImageDtt.QUAD; i++)  pass.disparity_map[ImageDtt.IMG_DIFF0_INDEX + i][nt] = 0.0;
					}
					if ((good_comp == null) || !good_comp[0]) {
						pass.calc_disparity[nt] = Double.NaN;
						pass.strength[nt] =       0.0;
						if (pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM]!=null )       pass.disparity_map[ImageDtt.DISPARITY_INDEX_CM][nt] = Double.NaN;
						if (pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX]!=null ) pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt] =      0.0;
					}
					if ((good_comp == null) || !good_comp[1]) {
						pass.calc_disparity[nt] = Double.NaN;
						pass.strength[nt] =       0.0;
						if (pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR]!=null )          pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt] = Double.NaN;
						if (pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH]!=null ) pass.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt] =      0.0;
					}
					if ((good_comp == null) || !good_comp[2]) {
						pass.calc_disparity[nt] = Double.NaN;
						pass.strength[nt] =       0.0;
						if (pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT]!=null )       pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt] = Double.NaN;
						if (pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH]!=null ) pass.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt] =      0.0;
					}
				}
			}
		}
	}

	public void filterOverBackground(
			final double [][]           ds,
			final double []             bg_strength,
			final boolean []            bg_tiles,          // get from selected in clt_3d_passes.get(0);
			final double                minStrength,
			final double                ex_min_over        // when expanding over previously detected (by error) background, disregard far tiles
			){
		for (int nt = 0; nt < bg_tiles.length; nt++)  if (bg_tiles[nt]){
			if ((ds[1][nt] <= minStrength ) || (ds[1][nt] < bg_strength[nt]) ||(ds[0][nt] < ex_min_over )){
				ds[1][nt] =0.0;
				ds[0][nt] =0.0;
			}
		}
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
	 * @param max_overexposure maximal fraction of (near) overexposed pixels in a tile to discard it
	 * @param disp_far        lowest disparity value to consider (does not apply to max_tried_disparity)
	 * @param disp_near       highest disparity value to consider (does not apply to max_tried_disparity)
	 * @param minStrength     full correlation strength to consider data to be reliable
	 * @param minStrengthHor  horizontal (for vertical features) correlation strength to consider data to be reliable. NaN - do not use
	 * @param minStrengthVert vertical (for horizontal features) correlation strength to consider data to be reliable. NaN - do not use
	 * @param no_weak         Do not use weak at all
	 * @param use_last        use last scan data if nothing strong enough (false - use the strongest)
	 * @param usePoly         use polynomial method to find max for full correlation, false - use center of mass
	 * @param copyDebug       copy data that is only needed for debug purposes
	 * @return new composite scan pass (not added to the list
	 */
	public CLTPass3d compositeScan(
			 final ArrayList <CLTPass3d> passes,
			 final int                   firstPass,
			 final int                   lastPassPlus1,
			 final double                trustedCorrelation,
			 final double                max_overexposure,
			 final double                disp_far,   // limit results to the disparity range
			 final double                disp_near,
			 final double                minStrength,
			 final double                minStrengthHor,
			 final double                minStrengthVert,
			 final boolean               no_weak,
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

						boolean overexposed =  (max_overexposure > 0.0) &&
								(pass.disparity_map[ImageDtt.OVEREXPOSED] != null) &&
								(pass.disparity_map[ImageDtt.OVEREXPOSED][nt] > max_overexposure);
						if (!overexposed) {
							double adiff =      Math.abs(mdisp);
							double adiff_hor =  Math.abs(mdisp_hor);
							double adiff_vert = Math.abs(mdisp_vert);

							if (adiff <= trustedCorrelation){
								double disp = mdisp/corr_magic_scale +  pass.disparity[ty][tx];
								// do not consider tiles over background if they are far and initially identified as background
								//							if ((bg_tiles == null) || !bg_tiles[nt] || (disp >= ex_min_over)) {
								if ((disp >= disp_far) && (disp <= disp_near) && !Double.isNaN(adiff)){
									if (strength >= minStrength) {
										if (!(adiff >= adiff_best)){ // adiff_best == Double.NaN works too
											adiff_best = adiff;
											best_index = ipass;
										}
									} else {
										if ((last && (strength > 0.0)) || (!no_weak && (strength > strongest_weak))){
											strongest_weak = strength;
											best_weak_index = ipass;
										}
									}
								}
								//							}
							}

							if (!Double.isNaN(minStrengthHor) && (adiff_hor <= trustedCorrelation)){
								double disp_hor = mdisp_hor/corr_magic_scale +  pass.disparity[ty][tx];
								if ((disp_hor >= disp_far) && (disp_hor <= disp_near) && !Double.isNaN(adiff_hor)){
									if (strength_hor >= minStrengthHor) {
										if (!(adiff_hor >= adiff_best_hor)){ // adiff_best == Double.NaN works too
											adiff_best_hor = adiff_hor;
											best_index_hor = ipass;
										}
									} else {
										if ((last && (strength_hor > 0.0))  || (!no_weak && (strength_hor > strongest_weak_hor))){
											strongest_weak_hor = strength_hor;
											best_weak_index_hor = ipass;
										}
									}
								}
							}

							if (!Double.isNaN(minStrengthVert) && (adiff_vert <= trustedCorrelation)){
								double disp_vert = mdisp_vert/corr_magic_scale +  pass.disparity[ty][tx];
								if ((disp_vert >= disp_far) && (disp_vert <= disp_near) && !Double.isNaN(adiff_vert)){
									if (strength_vert >= minStrengthVert) {
										if (!(adiff_vert >= adiff_best_vert)){ // adiff_best == Double.NaN works too
											adiff_best_vert = adiff_vert;
											best_index_vert = ipass;
										}
									} else {
										if ((last && (strength_vert > 0.0))  || (!no_weak && (strength_vert > strongest_weak_vert))){
											strongest_weak_vert = strength_vert;
											best_weak_index_vert = ipass;
										}
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
					// Only copy for full disparity
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
		combo_pass.getDisparity(); // See if it does not break anything - triggers calculation if not done yet
		combo_pass.fixNaNDisparity(); // mostly for debug, measured disparity should be already fixed from NaN
		return combo_pass;
	}

	/**
	 * Create next measurement scan that can handle multiple disparities, using quad correlations only
	 * @param passes
	 * @param firstPass
	 * @param lastPassPlus1
	 * @param trustedCorrelation
	 * @param disp_far
	 * @param disp_near
	 * @param minStrength
	 * @param unique_tolerance
	 * @param usePoly
	 * @param debugLevel
	 * @return
	 */

	public CLTPass3d RefineQuadMulti(
			 final ArrayList <CLTPass3d> passes,
			 final int                   firstPass,
			 final int                   lastPassPlus1,
			 final double                trustedCorrelation,
			 final double                disp_far,   // limit results to the disparity range
			 final double                disp_near,
			 final double                minStrength,
			 final double                unique_tolerance,

			 // TODO: when useCombo - pay attention to borders (disregard)
			 final boolean               usePoly,  // use polynomial method to find max), valid if useCombo == false
			 final int                   debugLevel)
	{
		final int dbg_tile = (debugLevel > 0)? -839: -1; // x = 122, y= 108; -1; // 27669;
		CLTPass3d combo_pass =new CLTPass3d(this);
		final int disparity_index = usePoly ? ImageDtt.DISPARITY_INDEX_POLY : ImageDtt.DISPARITY_INDEX_CM;
		combo_pass.tile_op =              new int [tilesY][tilesX];
		combo_pass.disparity =            new double [tilesY][tilesX];
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);
		int num_new = 0;
		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++){
				int nt = ty * tilesX + tx;
				if (nt == dbg_tile) {
					System.out.println("RefineQuadMulti(): nt = "+nt);
				}
				for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
					CLTPass3d pass = passes.get(ipass);
					if (nt == dbg_tile) {
						System.out.println("RefineQuadMulti(): ipass = "+ipass+" nt = "+nt+" pass.tile_op["+ty+"]["+tx+"]="+pass.tile_op[ty][tx]+
								" pass.isCombo()="+(pass.isCombo())+" pass.isProcessed()="+(pass.isProcessed()));
					}
					if ( pass.isMeasured() && (pass.tile_op[ty][tx] != 0 )) { // current tile has valid data


//						if (	(Double.isNaN(combo_pass.max_tried_disparity[ty][tx]) ||
//								(pass.disparity[ty][tx] > combo_pass.max_tried_disparity[ty][tx]))){
//							combo_pass.max_tried_disparity[ty][tx] = pass.disparity[ty][tx];
//						}
						double mdisp =         pass.disparity_map[disparity_index][nt];
						double strength =      pass.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][nt];
						double adiff =      Math.abs(mdisp);

						if ((strength >= minStrength) && (adiff <= trustedCorrelation)){
							double disp = mdisp/corr_magic_scale +  pass.disparity[ty][tx];
//							double disp_low =  (disp > pass.disparity[ty][tx]) ? pass.disparity[ty][tx] : (2 * disp - pass.disparity[ty][tx]);
//							double disp_high = (disp > pass.disparity[ty][tx]) ? (2 * disp - pass.disparity[ty][tx]) : pass.disparity[ty][tx];
							double disp_low =  Math.min(disp, pass.disparity[ty][tx]);
							double disp_high = Math.max(disp, pass.disparity[ty][tx]);
//							if ((disp_high - disp_low) > 2 * unique_tolerance) { // suggested correction is not too small
							if ((disp_high - disp_low) >= unique_tolerance) { // suggested correction is not too small
								boolean duplicate = false;
								for (int iother = firstPass;  iother <lastPassPlus1; iother++ ) {
									CLTPass3d other = passes.get(iother);
									if ( other.isMeasured() && (other.tile_op[ty][tx] != 0 )){ //  && (iother != ipass)){
										if (    (Math.abs(other.disparity[ty][tx] - disp) < unique_tolerance) ||
												((other.disparity[ty][tx] - disp) * (other.disparity[ty][tx] - pass.disparity[ty][tx]) <0) // between
												){
											duplicate = true;
											break;
										}
									}
								}
								if (!duplicate){
									combo_pass.tile_op[ty][tx] = op;
									combo_pass.disparity[ty][tx] = disp;
									num_new++;
									break; // for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
								}
							}
						}
					}
				}
			}
		}
//		combo_pass.getDisparity(); // See if it does not break anything - triggers calculation if not done yet
//		combo_pass.fixNaNDisparity(); // mostly for debug, measured disparity should be already fixed from NaN
		if (debugLevel > -1){
			System.out.println("RefineQuadMulti(): prepared "+num_new+" tiles to be measured");
		}
		return (num_new > 0) ? combo_pass: null;
	}



	/**
	 * Verify that selected points are not all on the same line
	 * @param sel 2-d sample selection in linescan order
	 * @param side square samples side
	 * @return true if there are enough samples for plane extraction, false otherwise
	 */
	public boolean notColinear (
			boolean [] sel,
			int side)
	{
		int indx0, indx1;
		for (indx0 = 0; indx0 < sel.length; indx0++){
			if (sel[indx0]) break;
		}
		for (indx1 = indx0+1; indx1 < sel.length; indx1++){
			if (sel[indx1]) break;
		}
		if (indx1 >= sel.length) return false; // too few points;
		int sx0 = indx0 % side;
		int sy0 = indx0 / side;
		int sx1 = indx1 % side;
		int sy1 = indx1 / side;
		for (int indx = indx1 +1; indx < sel.length; indx++){
			int sx = indx % side;
			int sy = indx / side;
			if ((sx - sx0) * (sy - sy1) != (sx - sx1) * (sy - sy0)){
				return true;
			}
		}
		return false;
	}

	/**
	 * Verify that selected points are not all on the same line, even if the specified one is removed
	 * @param indx index of the point to be removed
	 * @param sel 2-d sample selection in linescan order
	 * @param side square samples side
	 * @return true if there are enough samples for plane extraction, false otherwise
	 */

	public boolean notColinearWithout (
			int        indx,
			boolean [] sel,
			int side)
	{
		if (!sel[indx]){
			throw new IllegalArgumentException ("notCoplanarWithout(): specified is the non existing index");
		}
		sel[indx] = false;
		boolean rslt = notColinear ( sel, side);
		sel[indx] = true;
		return rslt;
	}

	/**
	 * Get window function for tile samples (currently just 3x3) in a line-scan order
	 * @param smplSide square sample side
	 * @return [smplSide * smplSide]array of weights, 1.0 in the center
	 */

	public double [] getSampleWindow(int smplSide, boolean all1){
		double [] weights = new double [smplSide * smplSide];
		for (int sy = 0; sy < smplSide; sy++){
			for (int sx = 0; sx < smplSide; sx++){
				weights[sy * smplSide + sx] = all1? 1.0 : Math.sin(Math.PI*(sy+0.5)/smplSide) * Math.sin(Math.PI*(sx+0.5)/smplSide);
			}
		}
		return weights;
	}



	/**
	 * Get disparity/strength for last measured tiles (measured_scan), using neighbors measured in the same or earlier
	 * scans, assuming the sample square can (almost) freely tilt for the best fit
	 * @param measured_scan_index index in the clt_3d_passes list of the latest measured scan
	 * @param start_scan_index lowest scan to use data from
	 * @param bg_tiles background tiles selection, may be null
	 * @param ex_min_over only consider tiles that are nearer than this, if they are previously identified (by error?) as background
	 * @param disp_index disparity index in disparity_map (normally ImageDtt.DISPARITY_INDEX_CM)
	 * @param str_index strength index in disparity_map (normally[ImageDtt.DISPARITY_STRENGTH_INDEX)
	 * @param tiltXY null (free floating) or fixed {tilt_x, tilt_y) of the sample
	 * @param trustedCorrelation maximal absolute value of measured correlation (no scaling) to trust (may use global trustedCorrelation)
	 * @param strength_floor subtract from the strength, discard negative
	 * @param strength_pow raise strength to this power
	 * @param smplSide side of tyhe sample square
	 * @param smplNum number of the best fit neighbors to keep
	 * @param smplRms maximal disparity RMS of the best remaining tiles to keep sample
	 * @param smplRelRms relative RMS to add to  smplRms: effective RMS = smplRms + disparity * smplRelRms
	 * @param smplWnd use sample window function (false - all set to 1.0)
	 * @param max_abs_tilt maximal disparity tilt in disparity difference (pix) per tile;
	 * @param max_rel_tilt same as max_abs_tilt, but divided by disparity
	 * @param dbg_x debug tile X
	 * @param dbg_y debug tile Y
	 * @param debugLevel debug level
	 * @return {disparity,strength};
	 */
//					combo_pass.calc_disparity[nt] =         pass.disparity_map[disparity_index][nt]/corr_magic_scale +  pass.disparity[ty][tx];

	public double [][] getFilteredDisparityStrength(
			final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
			final int        measured_scan_index, // will not look at higher scans (OK to be non-measured, last measured will be used then)
			final int        start_scan_index,
			final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
			final double     ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
			final int        disp_index,
			final int        str_index,
			final double []  tiltXY, // null - free with limit on both absolute (2.0?) and relative (0.2) values
		    final double     trustedCorrelation,
			final double     strength_floor,
			final double     strength_pow,
			final int        smplSide, //        = 2;      // Sample size (side of a square)
			final int        smplNum, //         = 3;      // Number after removing worst (should be >1)
			final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			final double     smplRelRms, //      = 0.005;  // Maximal RMS/disparity in addition to smplRms
			final boolean    smplWnd, //
			final double     max_abs_tilt, //  = 2.0; // pix per tile
			final double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
			final int        dbg_x,
			final int        dbg_y,
			final int        debugLevel)

	{
		final double [] smpl_weights = getSampleWindow(smplSide, !smplWnd);
		final int smplLen = smplSide*smplSide;
		final int smpl_center = smplSide / 2;
		final int center_index = smpl_center * ( 1 + smplSide);
		final double smpl_dcenter = (smplSide -1.0) /2;
//		final TileNeibs tnSurface = new TileNeibs(tilesX, tilesY);
		final int tiles = tilesX * tilesY;
		final boolean [] measured = new boolean [tiles];
//		final double smplVar = smplRms * smplRms; // maximal variance (weighted average of the squared difference from the mean)
		final double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		final double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
		final double [][] disp_strength = new double [2][tiles]; // strength < 0 - not measured!



		final int dbg_tile = (debugLevel > 0) ?(dbg_x + dbg_y * tilesX) : -1;
		// Create a list of only measurement scans
		final ArrayList<CLTPass3d> measured_list = new ArrayList<CLTPass3d> ();
		for (int ipass = start_scan_index; (ipass < passes.size()) && (ipass <= measured_scan_index) ; ipass++){
			CLTPass3d pass = passes.get(ipass);
			if ((pass != null) && pass.isMeasured()){
				measured_list.add(pass);
			}
		}
		final CLTPass3d measured_scan =measured_list.get(measured_list.size() - 1);// use last measured scan

		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
//						int dl = (nTile == dbg_tile) ? debugLevel : -1;
						int tileX = nTile % tilesX;
						int tileY = nTile / tilesX;
						measured[nTile] = measured_scan.tile_op[tileY][tileX] != 0;
						if (!measured[nTile]) disp_strength[1][nTile] = -1.0;
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);

		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					PolynomialApproximation pa = new PolynomialApproximation();
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (measured[nTile]){
						int dl = (nTile == dbg_tile) ? debugLevel : -1;
						int tileX = nTile % tilesX;
						int tileY = nTile / tilesX;
						// check that the tile is strong enough and fractional disparity is in the trusted range
						double w = measured_scan.disparity_map[str_index][nTile] - strength_floor;
//						if ((dl > 0) && (measured_scan.disparity_map[str_index][nTile] != 0)){
						if ((dl > 0) && (measured_scan.tile_op[tileY][tileX] != 0)){
							System.out.println("getFilteredDisparityStrength(): tileX="+tileX+", tileY="+tileY+" w="+w);
						}
						if (w > 0) {
							if (strength_pow != 1.0){
								if (strength_pow == 0.0) {
									w = 1.0;
								} else  {
									w= Math.pow(w, strength_pow);
								}
							}
							double disp = measured_scan.disparity_map[disp_index][nTile];
							if (Math.abs(disp) <= trustedCorrelation){
								disp = disp/corr_magic_scale + measured_scan.disparity[tileY][tileX];
								if ((bg_tiles == null) || !bg_tiles[nTile] || (disp >= ex_min_over)) {
									// Center tile is valid, makes sense to investigate neighbors from this or earlier tiles
									int num_in_sample = 1;
									double sum_wnd = 1.0;
									boolean [] smpl_sel = new boolean [smplLen];
									double [] smpl_d =    new double [smplLen];
									double [] smpl_p =    new double [smplLen];
									double [] smpl_w =    new double [smplLen];
									// add center tile
									smpl_sel[center_index] = true;
									smpl_d[center_index] = disp;
									smpl_w[center_index] = w;
									// process tiles around, looking for the best fit among other tiles, skipping the center
									for (int sy = 0; sy < smplSide; sy++){
										int y = tileY + sy - smpl_center;
										if ((y >= 0) && (y < tilesY)) {
											for (int sx = 0; sx < smplSide; sx++){
												int x = tileX + sx - smpl_center;
												if ((x >= 0) && (x < tilesX)) {
													int indxs = sy*smplSide + sx;
													if (indxs != center_index) { // already assigned
														int nTile1 = y * tilesX + x;
														int best_pass_index = -1;
														double best_adiff = Double.NaN;
														for (int ipass = 0; ipass < measured_list.size(); ipass++){
															CLTPass3d scan = measured_list.get(ipass);
															double w1 = scan.disparity_map[str_index][nTile1] - strength_floor;
															if (w1 > 0) {
																if (strength_pow != 1.0){
																	if (strength_pow == 0.0) {
																		w1 = 1.0;
																	} else  {
																		w1= Math.pow(w, strength_pow);
																	}
																}
																double disp1 = scan.disparity_map[disp_index][nTile1];
																if (Math.abs(disp1) <= trustedCorrelation){
																	disp1 = disp1/corr_magic_scale + scan.disparity[y][x];
																	double adiff = Math.abs(disp1 - disp);
																	if ((best_pass_index < 0) || (adiff < best_adiff)){
																		best_pass_index = ipass;
																		best_adiff = adiff;
																		smpl_d[indxs] = disp1;
																		smpl_w[indxs] = w1;
																		// not summed yet!
																	}
																}
															}
														}
														if (best_pass_index >= 0){
															smpl_w  [indxs] *= smpl_weights[indxs];
															smpl_sel[indxs] = true;
															sum_wnd += smpl_weights[indxs];
															num_in_sample ++;
														}

													}
												}
											}
										}
									}

									// filled all neighbor tiles data. See if sufficient to proceed, then remove extra tiles
									// and process the remaining ones
									if ((dl > 0) && (measured_scan.tile_op[tileY][tileX] != 0)){
										System.out.println("getFilteredDisparityStrength() 1: tileX="+tileX+", tileY="+tileY+" num_in_sample="+num_in_sample);
									}
									if (num_in_sample >= smplNum){ // try, remove worst
										sample_loop:
										{
										boolean en_tilt = (tiltXY == null);
										if (en_tilt) { // make sure there are enough samples and not all of them are on the same line
											if (!notColinear(smpl_sel,smplSide)){
												en_tilt = false;
											}
										}
										if (en_tilt) { // enable floating tilt
											double sd2 = 0.0, d_center = 0.0,  sw = 0.0;

											// TODO: making simple - recalculate after removing. Can be done more efficient.
											while (num_in_sample >= smplNum) { // try, remove worst
												sd2 = 0.0;
												d_center = 0.0;
												sw = 0.0;

												double [][][] mdata = new double [num_in_sample][3][];
												int mindx = 0;
												for (int sy = 0; sy < smplSide; sy++){
													for (int sx = 0; sx < smplSide; sx++){
														int indxs = sy * smplSide + sx;
														if (smpl_sel[indxs]) {
															mdata[mindx][0] = new double [2];
															mdata[mindx][0][0] =  sx - smpl_dcenter;
															mdata[mindx][0][1] =  sy - smpl_dcenter;
															mdata[mindx][1] = new double [1];
															mdata[mindx][1][0] = smpl_d[indxs];
															mdata[mindx][2] = new double [1];
															mdata[mindx][2][0] =  smpl_w[indxs];
															mindx ++;
														}
													}
												}
												double[][] approx2d = pa.quadraticApproximation(
														mdata,
														true,          // boolean forceLinear,  // use linear approximation
														thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
														thresholdQuad, // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
														debugLevel);
												if (approx2d == null){
													if (debugLevel > -1){
														System.out.println("getFilteredDisparityStrength(): can not find linear approximation");
													}
													break sample_loop;
												}
												// limit tilt to be within range
												//											double     max_abs_tilt, //  = 2.0; // pix per tile
												//											double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
												double max_tilt = Math.min(max_abs_tilt, max_rel_tilt * approx2d[0][2]);
												boolean overlimit = (Math.abs(approx2d[0][0]) > max_tilt) || (Math.abs(approx2d[0][1]) > max_tilt);
												if (overlimit) {
													approx2d[0][0] = Math.min(approx2d[0][0],  max_tilt);
													approx2d[0][1] = Math.min(approx2d[0][1],  max_tilt);
													approx2d[0][0] = Math.max(approx2d[0][0], -max_tilt);
													approx2d[0][1] = Math.max(approx2d[0][1], -max_tilt);
												}
												// subtract tilt from disparity
												for (int sy = 0; sy < smplSide; sy++){
													for (int sx = 0; sx < smplSide; sx++){
														int indxs = sy * smplSide + sx;
														if (smpl_sel[indxs]) {
															smpl_p[indxs] = approx2d[0][0] * (sx - smpl_dcenter) + approx2d[0][1] * (sy - smpl_dcenter) + approx2d[0][2];
														}
													}
												}

												if (overlimit){ // re-calculate disparity average (in the center)
													double sd=0.0;
													for (int indxs = 0; indxs < smplLen;indxs++) if (smpl_sel[indxs]) {
														double d = smpl_d[indxs] - smpl_p[indxs];
														double dw = d * smpl_w[indxs];
														sd += dw;
														sw += smpl_w[indxs];
													}
													sd /= sw;
													for (int indxs = 0; indxs < smplLen;indxs++) if (smpl_sel[indxs]) {
														smpl_p[indxs] += sd;
													}
													approx2d[0][2] += sd;

												}
												d_center = approx2d[0][2];
												sw = 0.0;
												for (int indxs = 0; indxs < smplLen;indxs++) if (smpl_sel[indxs]) {
													double d = smpl_d[indxs] - smpl_p[indxs];
													double dw = d * smpl_w[indxs];
													//									sd += dw;
													sd2 += dw * smpl_d[indxs];
													sw +=       smpl_w[indxs];
												}


												// remove worst - it should not make remaining set
												if (num_in_sample > smplNum) { // remove worst if it is not the last run where only calculations are needed
													//												double d_mean = sd/sw;
													int iworst = -1;
													double dworst2 = 0.0;
													for (int indxs = 0; indxs < smplLen; indxs++) if (smpl_sel[indxs]) {
														//													double d2 = (smpl_d[i] - d_mean);
														double d2 = smpl_d[indxs] - smpl_p[indxs];
														d2 *=d2;
														if (d2 > dworst2) {
															if (notColinearWithout (
																	indxs, // int        indx,
																	smpl_sel, // boolean [] sel,
																	smplSide)) { // int side))
																iworst = indxs;
																dworst2 = d2;
															}
														}
													}
													if (iworst < 0){
														System.out.println("**** this is a BUG in getFilteredDisparityStrength() can not find the worst sample ****");
														break;
													}
													// remove worst sample
													smpl_sel[iworst] = false;
													//												double dw = smpl_d[iworst] * smpl_w[iworst];
													//												sd -= dw;
													//												sd2 -= dw * smpl_d[iworst];
													//												sw -=       smpl_w[iworst];
													sum_wnd -= smpl_weights[iworst];
													num_in_sample --;
												} else {
													break;
												}

											} // removing worst tiles, all done,
											// calculate variance of the remaining set
											if (sw > 0.0) {
												//											sd /= sw;
												//											sd2 /= sw;
												double var = sd2/sw; //   - sd * sd;
												double smplVar = smplRms + d_center * smplRelRms;
												smplVar *= smplVar;
												if (var < smplVar) { // good, save in the result array
													disp_strength[0][nTile] = d_center;
													//								ds[1][indx] = sw * lapWeight[dy][dx] /num_in_sample; // average weights, multiply by window //** TODO: change
													disp_strength[1][nTile] = sw /sum_wnd; // average weights, multiply by window //** TODO: change
												}
												if (dl > 0){
													System.out.println("getFilteredDisparityStrength() 2: tileX="+tileX+", tileY="+tileY+" num_in_sample="+num_in_sample+" var="+var);
												}
											}

										} else { // fixed tilt
											// tilt around center
											if ((tiltXY != null) && ((tiltXY[0] != 0.0) || (tiltXY[1] != 0.0))){
												for (int sy = 0; sy < smplSide; sy++){
													for (int sx = 0; sx < smplSide; sx++){
														int indxs = sy * smplSide + sx;
														if (smpl_w[indxs] > 0.0) {
															smpl_d[indxs] -= tiltXY[0]* (sx - smpl_dcenter) + tiltXY[1]* (sy - smpl_dcenter);
														}
													}
												}
											}


											// calculate
											double sd=0.0, sd2 = 0.0, sw = 0.0;
											for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
												double dw = smpl_d[i] * smpl_w[i];
												sd += dw;
												sd2 += dw * smpl_d[i];
												sw +=       smpl_w[i];
											}
											// remove worst, update sd2, sd and sw
											while ((num_in_sample > smplNum) && (sw > 0)){ // try, remove worst
												double d_mean = sd/sw;
												int iworst = -1;
												double dworst2 = 0.0;
												for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
													double d2 = (smpl_d[i] - d_mean);
													d2 *=d2;
													if (d2 > dworst2) {
														iworst = i;
														dworst2 = d2;
													}
												}
												if (iworst < 0){
													System.out.println("**** this is a BUG in getFilteredDisparityStrength() ****");
													break;
												}
												// remove worst sample
												smpl_sel[iworst] = false;
												double dw = smpl_d[iworst] * smpl_w[iworst];
												sd -= dw;
												sd2 -= dw * smpl_d[iworst];
												sw -=       smpl_w[iworst];
												sum_wnd -= smpl_weights[iworst];
												num_in_sample --;
											}
											// calculate variance of the remaining set
											if (sw > 0.0) {
												sd /= sw;
												sd2 /= sw;
												double var = sd2 - sd * sd;
												double smplVar = smplRms + sd * smplRelRms;
												smplVar *= smplVar;
												if (var < smplVar) { // good, save in the result array
													disp_strength[0][nTile] = sd;
													//								ds[1][indx] = sw * lapWeight[dy][dx] /num_in_sample; // average weights, multiply by window //** TODO: change
													disp_strength[1][nTile] = sw  /sum_wnd; // average weights, multiply by window //** TODO: change
												}
											} else {
												num_in_sample = 0;
												System.out.println("**** this is a BUG in getFilteredDisparityStrength(), shoud not happen ? ****");
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
		return disp_strength;
	}




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
			 final double                grow_disp_max,
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
		int dbg_remain = 0;
		for (int ty = 0; ty < tilesY; ty ++) {
			for (int tx = 0; tx < tilesX; tx ++){
				if (new_scan.tile_op[ty][tx] != 0){
					total ++;
					for (int ipass = firstPass;  ipass <lastPassPlus1; ipass++ ){
						CLTPass3d pass = passes.get(ipass);
						if ((pass.tile_op[ty][tx] != 0) && pass.isMeasured()){
							// Double.isNaN should not happen
							if (     Double.isNaN(pass.disparity[ty][tx]) ||
									(Math.abs(new_scan.disparity[ty][tx] - pass.disparity[ty][tx]) < unique_tolerance) ||
									(new_scan.disparity[ty][tx] > grow_disp_max) ||
									(new_scan.disparity[ty][tx] <= 0.0)){
								new_scan.tile_op[ty][tx] = 0;
								removed++;
								break;
							}
						}
					}
					if (new_scan.tile_op[ty][tx] != 0){
						dbg_remain++;
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

	/**
	 * Get selection of the tiles that were tried before among selected
	 * @param passes ArrayList of scan passes
	 * @param firstPass index of the first pass to compare against
	 * @param lastPassPlus1 index of the last pass + 1to compare against
	 * @param selected tiles to process
	 * @param prohibited tiles already marked as "bad" (already tried or out of range), may be null. If not null - will be modified
	 * @param disparity array of disparity values to compare with already tried. Double.NaN values will be skipped (not marked)
	 * @param mod_selected if true, selected array will be modified to exclude tried before cells
	 * @param mod_disparity if true, disparity array will be modified to have Double.NaN for the tried before cells
	 * @param en_near permit increasing disparity if that one was already tried
	 * @param en_far permit decreasing disparity if that one was already tried (if both en_near and en_far - try near first)
	 * @param grow_disp_min maximal valid disparity (all below will be marked as "bad" (as if they were already tried)
	 * @param grow_disp_max maximal valid disparity (all above  will be marked as "bad" (as if they were already tried)
	 * @param grow_disp_step disparity increase/decrease step if that one was used and en_near or en_far are true
	 * @param unique_pre_tolerance mark tiles that have disparity within +/- this values from already tried. Set to 0 to disable (only remove invalid)
	 * @return array of the tiles that should not be used
	 */

	public boolean [] getTriedBefore(
			 final ArrayList <CLTPass3d> scan_list,
			 final int                   start_scan_index, // firstPass,
			 final int                   finish_scan_index, // lastPassPlus1,
			 final boolean []            selected,
			 final boolean []            prohibited,
			 final double []             disparity,
			 final boolean               mod_selected,
			 final boolean               mod_disparity,
			 final boolean               en_near,
			 final boolean               en_far,
			 final double                grow_disp_min,
			 final double                grow_disp_max,
			 final double                grow_disp_step,
			 final double                unique_pre_tolerance)
	{
		final int tilesX = getTilesX();
		final int tilesY = getTilesY();
		final int tlen = tilesX * tilesY;
		final boolean [] tried_before = (prohibited != null)? prohibited : (new boolean [tlen]);
		final ArrayList<CLTPass3d> measured_list = new ArrayList<CLTPass3d> ();
		for (int ipass = start_scan_index; ipass < finish_scan_index; ipass++){
			CLTPass3d pass = clt_3d_passes.get(ipass);
			if ((pass != null) && pass.isMeasured()){
				measured_list.add(pass);
			}
		}
		final int firstPass = 0;
		final int lastPassPlus1 = measured_list.size();

		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tlen; nTile = ai.getAndIncrement()) if (selected[nTile]){
						if (!Double.isNaN(disparity[nTile])) {
							if ((disparity[nTile] <= grow_disp_min) || (disparity[nTile] > grow_disp_max)){
								tried_before[nTile] = true;
								if (mod_selected) selected[nTile] = false;
								if (mod_disparity) disparity[nTile] = Double.NaN;
							} else {
								int ty = nTile / tilesX;
								int tx = nTile % tilesX;
								boolean duplicate = false;
								for (int i =  firstPass; i < lastPassPlus1; i++) {
									if ((unique_pre_tolerance>0) &&
											(Math.abs(disparity[nTile] - measured_list.get(i).disparity[ty][tx]) < unique_pre_tolerance)){
										duplicate = true;
										break;
									}
								}
								// Try to replace a duplicate disparity to try with a larger one
								if (duplicate && en_near) {
									for (double disp = disparity[nTile] + grow_disp_step; disp <= grow_disp_max;  disp += grow_disp_step){
										boolean step_duplicate = false;
										for (int i =  firstPass; i < lastPassPlus1; i++) {
											if (Math.abs(disp - measured_list.get(i).disparity[ty][tx]) < unique_pre_tolerance){
												step_duplicate = true;
												break;
											}
										}
										if (!step_duplicate){
											duplicate = false;
											if (mod_disparity) disparity[nTile] = disp;
											break;
										}
									}
								}
								// Try to replace a duplicate disparity to try with a smaller one (including failed larger one)
								if (duplicate && en_far) {
									for (double disp = disparity[nTile] - grow_disp_step; disp > grow_disp_min;  disp -= grow_disp_step){
										boolean step_duplicate = false;
										for (int i =  firstPass; i < lastPassPlus1; i++) {
											if (Math.abs(disp - measured_list.get(i).disparity[ty][tx]) < unique_pre_tolerance){
												step_duplicate = true;
												break;
											}
										}
										if (!step_duplicate){
											duplicate = false;
											if (mod_disparity) disparity[nTile] = disp;
											break;
										}
									}
								}

								if (duplicate) {
									tried_before[nTile] = true;
									if (mod_selected) selected[nTile] = false;
									if (mod_disparity) disparity[nTile] = Double.NaN;
								}

							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return tried_before;
	}



	public int [] prepareExpandVariant(
			final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
			final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
			final CLTPass3d   bg_scan,         // background scan data, null - ignore background
			final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
			final int         firstPass,
			final int         lastPassPlus1,
			final int         min_new,         // discard variant if there are less new tiles
			final boolean []  variants_flags,
			final int         num_steps,  // how far to extend
			final int         num_steps_disc,  // how far to extend
			final int         smpl_size, // == 5
			final int         min_points, // == 3
			final boolean     use_wnd,   // use window function fro the neighbors
			final double      tilt_cost,
			final double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double      same_range,      // modify
			final double      diff_continue,   // maximal difference from the old value (for previously defined tiles
			final double      max_abs_tilt, //   = 2.0; // pix per tile
			final double      max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
			final int         max_tries, // maximal number of smoothing steps
			final double      final_diff, // maximal change to finish iterations
		    final double      grow_disp_min, // minimal disparity to use (0.0?)
			final double      grow_disp_max, // maximal disparity to try
			final double      grow_disp_step, // disparity step (increment/decrement if that one was already tried)
			final double      unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
			final double      unique_tolerance,
			final boolean     show_expand,
			final boolean     show_unique,
			final int         dbg_x,
			final int         dbg_y,
			final int         debugLevel)
	{
		class Variant{
			boolean     expand_discontinuity; // true - expand discontinuity over known tiles
			boolean     smooth_only;   // Do not expand ambiguous cells (with discontinuity) (only valid for expand_discontinuity = false)
			boolean     extend_far;    // extend far (false - extend near)
			boolean     try_next;      // try next closer disparity if requested was tried before
			int         num_steps;  // how far to extend
			double      split_threshold; // if full range of the values around the cell higher, need separate fg, bg
			double      same_range;      // modify
			double      diff_continue;   // maximal difference from the old value (for previously defined tiles

			Variant (
					boolean     expand_discontinuity, // true - expand discontinuity over known tiles
					boolean     smooth_only,   // Do not expand ambiguous cells (with discontinuity) (only valid for expand_discontinuity = false)
					boolean     extend_far,    // extend far (false - extend near)
					boolean     try_next,
					int         num_steps,  // how far to extend
					double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
					double      same_range,      // modify
					double      diff_continue)
			{
				this.expand_discontinuity = expand_discontinuity; // true - expand discontinuity over known tiles
				this.smooth_only =          smooth_only;   // Do not expand ambiguous cells (with discontinuity) (only valid for expand_discontinuity = false)
				this.extend_far =           extend_far;    // extend far (false - extend near)
				this.try_next =             try_next;
				this.num_steps =            num_steps;  // how far to extend
				this.split_threshold =      split_threshold; // if full range of the values around the cell higher, need separate fg, bg
				this.same_range =           same_range;
				this.diff_continue =        diff_continue;
			}
		}
		Variant [] variants = {
				//          disc  smooth   far   next
				new Variant(false, true,  false, false, num_steps,      split_threshold, same_range, diff_continue),
				new Variant(false, false, false, false, num_steps,      split_threshold, same_range, diff_continue),
				new Variant(true,  false, false, false, num_steps_disc, split_threshold, same_range, diff_continue),
				new Variant(false, false, true,  false, num_steps,      split_threshold, same_range, diff_continue),
				new Variant(true,  false, true,  false, num_steps_disc, split_threshold, same_range, diff_continue),
				new Variant(false, false, false, true,  num_steps_disc, split_threshold,  same_range, diff_continue),
		};
		for (int num_var = 0; num_var < variants.length; num_var++){
			if (variants_flags[num_var]) { // (variants_mask & (1 << num_var)) != 0) {
//				scan.saveTileOpDisparity(); //  it is not used from this scan, but from last_scan
				int num_op = setupExpand(
						scan,                                   // final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
						last_scan,                              // final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
						bg_scan,                                // final CLTPass3d   bg_scan,         // background scan data, null - ignore background
						passes,                                 //  final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
						firstPass,                              // final int         firstPass,
						lastPassPlus1,                          // final int         lastPassPlus1,
						variants[num_var].expand_discontinuity, // final boolean     expand_discontinuity, // true - expand discontinuity over known tiles
						variants[num_var].smooth_only,          // final boolean     smooth_only,   // Do not expand ambiguous cells (with discontinuity) (only valid for expand_discontinuity = false)
						variants[num_var].extend_far,           // final boolean     extend_far,    // extend far (false - extend near)
						variants[num_var].try_next,             // extend far (false - extend near)
						variants[num_var].num_steps,            // final int         num_steps,  // how far to extend
						smpl_size,                              // final int         smpl_size, // == 5
						min_points,                             // final int        min_points, // == 3
						use_wnd,                                // final boolean     use_wnd,   // use window function fro the neighbors
						tilt_cost,                              // final double      tilt_cost,
						variants[num_var].split_threshold,      //  final double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
						variants[num_var].same_range,           //  final double      same_range,      // modify
						variants[num_var].diff_continue,        // final double      diff_continue,   // maximal difference from the old value (for previously defined tiles
						max_abs_tilt,                           // final double     max_abs_tilt, //   = 2.0; // pix per tile
						max_rel_tilt,                           // final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
						max_tries,                              // final int         max_tries, // maximal number of smoothing steps
						final_diff,                             // final double      final_diff, // maximal change to finish iterations
						grow_disp_min,                          // final double     grow_disp_min,
						grow_disp_max,                          // final double     grow_disp_max,
						grow_disp_step,                         // final double     grow_disp_step,
						unique_pre_tolerance,                   // final double     unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
						show_expand,                            // final boolean     show_expand,
						dbg_x,                                  // final int         dbg_x,
						dbg_y,                                  // final int         dbg_y,
						debugLevel);                            // final int         debugLevel)
				if (num_op < min_new) {
//					scan.restoreKeepTileOpDisparity(); // it is not used from this scan, but from last_scan
					continue;
				}

				// Remove tiles that were already tried, see if anything is left
				int [] unique_stats = makeUnique(
						passes, // clt_3d_passes,       // final ArrayList <CLTPass3d> passes,
						firstPass,           //  final int                   firstPass,
						lastPassPlus1, // clt_3d_passes.size(),// last_scan,           // - 1,                   //  final int                   lastPassPlus1,
						scan,                //clt_3d_passes.get(refine_pass), //  final CLTPass3d             new_scan,
						grow_disp_max,       // final double                grow_disp_max,
						unique_tolerance,    //  final double                unique_tolerance,
						show_unique);        // final boolean               show_unique)
				int num_left = unique_stats[0];
				if (num_left < min_new) {
					if (debugLevel > 0){
						System.out.println("prepareExpandVariants(): discarding variant "+num_var+ " as remaining "+num_left+" < "+min_new);
					}
//					scan.restoreKeepTileOpDisparity();//  it is not used from this scan, but from last_scan
					continue;
				}
				if (debugLevel > -10){
					System.out.println("prepareExpandVariants(): remaining "+num_left+" tiles to be processed, used variant "+num_var+" ("+num_left+")");
				}
				int [] rslt = {num_left, num_var};
				return rslt;
			}
		}
		return new int [2]; // nothing found
	}


	public int setupExpand(
			final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
			final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
			final CLTPass3d   bg_scan,         // background scan data, null - ignore background
			final ArrayList <CLTPass3d> passes,// List, first, last - to search for the already tried disparity
			final int         firstPass,
			final int         lastPassPlus1,
			final boolean     expand_discontinuity, // true - expand discontinuity over known tiles
			final boolean     smooth_only,   // Do not expand ambiguous cells (with discontinuity) (only valid for expand_discontinuity = false)
			final boolean     extend_far,    // extend far (false - extend near)
			final boolean     try_next,
			final int         num_steps,  // how far to extend
			final int         smpl_size, // == 5
			final int         min_points, // == 3
			final boolean     use_wnd,   // use window function fro the neighbors
			final double      tilt_cost,
			final double      split_threshold, // if full range of the values around the cell higher, need separate fg, bg
			final double      same_range,      // modify
			final double      diff_continue,   // maximal difference from the old value (for previously defined tiles
			final double      max_abs_tilt, //   = 2.0; // pix per tile
			final double      max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
			final int         max_tries, // maximal number of smoothing steps
			final double      final_diff, // maximal change to finish iterations
			final double      grow_disp_min,
			final double      grow_disp_max,
			final double      grow_disp_step,
			final double      unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
			final boolean     show_expand,
			final int         dbg_x,
			final int         dbg_y,
			final int         debugLevel)
	{

		final int tlen = tilesY * tilesX;
		double [] disparity =     scan.getDisparity().clone(); // to modify in-place **** remove clone after debugged ****
//		double [] dbg_disparity = scan.getDisparity().clone(); // to modify in-place **** remove clone after debugged ****
		boolean [] these_no_border = new boolean [tlen];
		for (int i = 0; i < these_no_border.length; i++) {
			these_no_border[i] = last_scan.selected[i] && !last_scan.border_tiles[i];
		}

		boolean [] dbg_last_selected = last_scan.selected.clone();
		boolean [] dbg_last_border =   last_scan.border_tiles.clone();
		boolean [] dbg_no_border =     these_no_border.clone();

		boolean [] known_tiles = these_no_border.clone();
		// known are background or these tiles
		if (bg_scan != null) {
			for (int i = 0; i < known_tiles.length; i++) {
				known_tiles[i] |= bg_scan.selected[i];
			}
		}
		// set combo disparity from  last prepared
		for (int nt = 0; nt < known_tiles.length; nt++){
			int ty = nt / tilesX;
			int tx = nt % tilesX;
			disparity[nt] = 0.0;
			if (these_no_border[nt]) disparity[nt] = last_scan.disparity[ty][tx];
		}
		double [] dbg_disparity = disparity.clone();

		boolean [] expanded;
		ExtendSurfaces es = new ExtendSurfaces(this);
		if (expand_discontinuity) {
			expanded = es.expandDiscontinuity( // return maximal correction value
					passes,          // final ArrayList <CLTPass3d> passes,
					firstPass,       // final int                   firstPass,
					lastPassPlus1,   // final int                   lastPassPlus1,
					disparity,       // final double []  values,     // will be modified, Double.NaN is not yet assigned
					known_tiles,     // final boolean [] known,      // cells with known values (will not be modified)
					num_steps,       // final int        num_steps,  // how far to extend
					smpl_size,       // final int        smpl_size, // == 5
					min_points,      // final int        min_points, // == 3
					use_wnd,         // final boolean    use_wnd,   // use window function fro the neighbors
					tilt_cost,       // final double     tilt_cost,
					split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
					same_range,      // final double     same_range,      // modify
					diff_continue,   // final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
					max_abs_tilt,    // final double     max_abs_tilt, //   = 2.0; // pix per tile
					max_rel_tilt,    // final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
					extend_far,      // final boolean    extend_far,    // extend far (false - extend near)
					max_tries,       // final int        max_tries, // maximal number of smoothing steps
					final_diff,      // final double     final_diff, // maximal change to finish iterations
					grow_disp_min,   // final double     grow_disp_min,
					grow_disp_max,   // final double     grow_disp_max,
					grow_disp_step,  // final double     grow_disp_step,
					unique_pre_tolerance, // final double     unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
					try_next,        // final boolean    try_next, // try next disparity range if the current one was already tried
					dbg_x,           // final int        dbg_x,
					dbg_y,           // final int        dbg_y,
					debugLevel);     // final int        debugLevel)

		} else {
			expanded = es.expandKnown( // return maximal correction value
					passes,          // final ArrayList <CLTPass3d> passes,
					firstPass,       // final int                   firstPass,
					lastPassPlus1,   // final int                   lastPassPlus1,
					disparity,       // final double []  values,     // will be modified, Double.NaN is not yet assigned
					known_tiles,     // final boolean [] known,      // cells with known values (will not be modified)
					num_steps,       // final int        num_steps,  // how far to extend
					smpl_size,       // final int        smpl_size, // == 5
					min_points,      // final int        min_points, // == 3
					use_wnd,         // final boolean    use_wnd,   // use window function fro the neighbors
					tilt_cost,       // final double     tilt_cost,
					split_threshold, // final double     split_threshold, // if full range of the values around the cell higher, need separate fg, bg
					same_range,      // final double     same_range,      // modify
					diff_continue,   // final double     diff_continue,   // maximal difference from the old value (for previously defined tiles
					max_abs_tilt,    // final double     max_abs_tilt, //   = 2.0; // pix per tile
					max_rel_tilt,    // final double     max_rel_tilt, //   = 0.2; // (pix / disparity) per tile
					smooth_only,     // final boolean    smooth_only,   // Do not expand ambiguous cells (with discontinuity)
					extend_far,      // final boolean    extend_far,    // extend far (false - extend near)
					max_tries,       // final int        max_tries, // maximal number of smoothing steps
					final_diff,      // final double     final_diff, // maximal change to finish iterations
					grow_disp_min,   // final double     grow_disp_min,
					grow_disp_max,   // final double     grow_disp_max,
					grow_disp_step,  // final double     grow_disp_step,
					unique_pre_tolerance, // final double     unique_pre_tolerance, // usually larger than clt_parameters.unique_tolerance
					try_next,        // final boolean    try_next, // try next disparity range if the current one was already tried
					dbg_x,           // final int        dbg_x,
					dbg_y,           // final int        dbg_y,
					debugLevel);     // final int        debugLevel)
		}
		// combine what was prepared by refine and here
		if (expanded == null) return 0;
		boolean [] selected = new boolean[tlen];
		boolean [] border =   new boolean[tlen];
		for (int nTile = 0; nTile < tlen; nTile++){
			//			selected[nTile] = known_tiles[nTile] || last_scan.selected[nTile] || expanded[nTile];
//			if ((last_scan.selected == null) || (expanded == null) || (disparity == null)) {
//				System.out.println( " BUG in setupExpand()");
//			}  else {
				selected[nTile] = (last_scan.selected[nTile] || expanded[nTile]) && !Double.isNaN(disparity[nTile]); // selected by refine and new java.lang.NullPointerException
				border[nTile] = selected[nTile] && !known_tiles[nTile];
//			}
		}
		scan.setSelected(selected);
		scan.setBorderTiles(border);
		int num_op = scan.setTileOpDisparity(
				selected,   // boolean [] selection,
				disparity); // double []  disparity)
		if (debugLevel > 0){
			System.out.println("setupExpand(): set "+num_op+" tiles to be processed (some will be removed)");
		}
		// show debug images
		if (show_expand || try_next){
			String [] titles = {"disp_prev", "disparity", "masked", "sel_border","expanded","tile_op","pre-disparity", "this_tile_op", "this_pre-disparity","filt_disp","filt_str"};
			double [][] dbg_img = new double [titles.length][];
			dbg_img[0] = dbg_disparity;
			dbg_img[1] = disparity;
			dbg_img[2] = disparity.clone();
			dbg_img[3] = new double[tlen];
			dbg_img[4] = new double[tlen];
			dbg_img[5] = new double[tlen];
			dbg_img[6] = new double[tlen];
			dbg_img[7] = new double[tlen];
			dbg_img[8] = new double[tlen];
			if (this.dbg_filtered_disp_strength != null) {
				dbg_img[ 9] = this.dbg_filtered_disp_strength[0];
				dbg_img[10] = this.dbg_filtered_disp_strength[1];
			}


			for (int i = 0; i < tlen; i++) {
				dbg_img[3][i] = (dbg_last_selected[i]? 1.0:0.0)+ (dbg_last_border[i]? 2.0:0.0);
				if (expanded != null) {
					dbg_img[4][i] = (expanded[i]? 1.0:0.0);
					if (!expanded[i]) dbg_img[2][i] = Double.NaN;
				}
			}
			if (last_scan.tile_op!= null){
				for (int ty = 0; ty < tilesY; ty++){
					for (int tx = 0; tx < tilesX; tx++){
						if (last_scan.tile_op[ty] != null) dbg_img[5][ty*tilesX + tx] = last_scan.tile_op[ty][tx];
					}
				}
			}
			if (last_scan.disparity!= null){
				for (int ty = 0; ty < tilesY; ty++){
					for (int tx = 0; tx < tilesX; tx++){
						if (last_scan.disparity[ty] != null) dbg_img[6][ty*tilesX + tx] = last_scan.disparity[ty][tx];
					}
				}
			}
			if (scan.tile_op!= null){
				for (int ty = 0; ty < tilesY; ty++){
					for (int tx = 0; tx < tilesX; tx++){
						if (scan.tile_op[ty] != null) dbg_img[7][ty*tilesX + tx] = scan.tile_op[ty][tx];
					}
				}
			}
			if (scan.disparity!= null){
				for (int ty = 0; ty < tilesY; ty++){
					for (int tx = 0; tx < tilesX; tx++){
						if (scan.disparity[ty] != null) dbg_img[8][ty*tilesX + tx] = scan.disparity[ty][tx];
					}
				}
			}

			(new showDoubleFloatArrays()).showArrays(dbg_img,  tilesX, tilesY, true, "setupExpand",titles);
		}

		return num_op;
	}


	public int setupExtendDisparity( // returns number of new tiles to try
			final CLTPass3d   scan,            // combined scan with max_tried_disparity, will be modified to re-scan
			final CLTPass3d   last_scan,       // last prepared tile - can use last_scan.disparity, .border_tiles and .selected
			final CLTPass3d   bg_scan,         // background scan data, null - ignore background
			final int         grow_sweep,      // 8; // Try these number of tiles around known ones
			final double      grow_disp_max,   //  =   50.0; // Maximal disparity to try
			final double      tried_margin,    //  =  4.0; // consider already tried if within this margin from already tried
			final double      grow_disp_step,  //  =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?
			// not used
			final double      grow_min_diff,   // =    0.5; // Grow more only if at least one channel has higher variance from others for the tile
			final boolean     grow_retry_far,  // Retry tiles around known foreground that have low max_tried_disparity
			final boolean     grow_pedantic,   // Scan full range between max_tried_disparity of the background and known foreground
			final boolean     grow_retry_inf,  // Retry border tiles that were identified as infinity earlier

			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			GeometryCorrection geometryCorrection,
			final boolean     show_debug,
			final int         threadsMax,      // maximal number of threads to launch
			final boolean     updateStatus,
			final int         debugLevel)
	{
		final int dbg_tile = -1; // 39379; // 54627; // ty=159, tx = 249
		final int tlen = tilesY * tilesX;
		final boolean retryTwoSteps = true; // false; //  true;

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
		if (bg_scan != null) {
			for (int i = 0; i < known_tiles.length; i++) {
				known_tiles[i] |= bg_scan.selected[i];
			}
		}
		// set combo disparity from  last prepared
		for (int nt = 0; nt < known_tiles.length; nt++){
			int ty = nt / tilesX;
			int tx = nt % tilesX;
			disparity[nt] = 0.0;
			if (these_no_border[nt]) disparity[nt] = last_scan.disparity[ty][tx];
		}

		boolean [] untested_bgnd = null;
		if (grow_retry_far){
			if (retryTwoSteps) {
			untested_bgnd = scan.getUntestedBackgroundBorder2 (
					known_tiles,    // final boolean    [] known,
					disparity,      // final double     [] disparity,
					grow_disp_step, // final double        grow_disp_step,
					debugLevel);    // final int     debugLevel)
			} else {
				untested_bgnd = scan.getUntestedBackgroundBorder (
						known_tiles,    // final boolean    [] known,
						disparity,      // final double     [] disparity,
						grow_disp_step, // final double        grow_disp_step,
						debugLevel);    // final int     debugLevel)

			}
			if (!grow_retry_inf && (bg_scan == null)){
				for (int i = 0; i < untested_bgnd.length; i++) if (bg_scan.selected[i]) untested_bgnd[i] = false;
			}
			// turn these tiles as if they are not known (really known or belong to bgnd)
			for (int i = 0; i < known_tiles.length; i++) known_tiles[i] &= !untested_bgnd[i];
			if (debugLevel > -1){
				int num_untested_far = 0;
				for (int i = 0; i < untested_bgnd.length; i++) if (untested_bgnd[i]) num_untested_far++;
				System.out.println("setupExtendDisparity(): Discovered "+num_untested_far+" tiles that may have disparity between bg and (yet untested) fg");
			}
		}

		// now known_tiles have untested background removed. most likely will be added when growing known

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

		double [] dbg_pre_disp = null;
			if (show_debug) {
				dbg_pre_disp = scan.getDisparity().clone();
			}
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

		// replace  	untested_bgnd[] tiles with max_tried_disparity + grow_disp_step if it is below disparity[...]
		if (grow_retry_far && grow_pedantic) {
			int num_pedantic = 0;
//			disparity = scan.getDisparity();
			for (int nTile = 0; nTile < untested_bgnd.length; nTile++) if  (untested_bgnd[nTile]) {
				int tX = nTile% tilesX;
				int tY = nTile% tilesY;
				double disp = scan.max_tried_disparity[tY][tX] + grow_disp_step;
				if (disp < disparity[nTile]){
					disparity[nTile] = disp;
					num_pedantic++;
				}
				if (debugLevel > -1){
					System.out.println("setupExtendDisparity(): replaced "+num_pedantic+" tiles to discover possible planes between fg and bg");
				}

			}


		}

		// increase disparity if new suggested value was already tried
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
// Wrong calculation of num_extended?

		int num_op_tiles = 0;

		for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
			int indx =  tilesX * ty + tx;
			if (scan.selected[indx]) {
				scan.disparity[ty][tx] = disparity[indx];
				scan.tile_op[ty][tx] = op;
				num_op_tiles ++;
			} else {
				scan.disparity[ty][tx] = 0.0;
				scan.tile_op[ty][tx] = 0;
			}
		}


		if (show_debug){
			String [] dbg_titles0 = {
					"tried",          // 0
					"pre_disp",       // 1
					"disparity",      // 2
					"strength",       // 3
					"bgnd",           // 4
					"these",          // 5
					"known_in",       // 6
					"known",          // 7
					"last_sel_border",// 8
					"no_border",      // 9
					"tried_before",   // 10
					"retest",         // 11
					"diff2max",       // 12
					"diff2maxAvg",    // 13
					"normStrength"};  //14

			dbg_titles = dbg_titles0;
			dbg_img = new double[dbg_titles.length][tilesY * tilesX];
			//dbg_pre_disp
			double [] strength = scan.getStrength();
			double boolean_val = 50.0;
			double [] diff2max = scan.getSecondMaxDiff(false);
			double [] diff2maxAvg = scan.getSecondMaxDiff(true);

			dbg_img[1] = dbg_pre_disp;
			for (int i = 0; i <tlen; i++){
				int ty = i / tilesX;
				int tx = i % tilesX;
				dbg_img[0][i] = scan.max_tried_disparity[ty][tx];
				dbg_img[2][i] = scan.disparity[ty][tx];
				dbg_img[3][i] = strength[i];
				if (bg_scan != null) dbg_img[4][i] = ((bg_scan.selected != null) &&    (bg_scan.selected[i]))? boolean_val:0.0 ;
				dbg_img[5][i] = ((last_scan.selected != null) &&  (last_scan.selected[i]))? boolean_val:0.0 ;
				dbg_img[6][i] = 0.05 * (dbg_img[4][i] + dbg_img[5][i]);
				dbg_img[7][i] = (scan.selected[i]?4:0)+ (scan.border_tiles[i]?8:0);

				dbg_img[8][i] = (dbg_last_selected[i]?4:0)+ (dbg_last_border[i]?8:0);
				dbg_img[9][i] = dbg_no_border[i]?boolean_val:0;
				dbg_img[10][i] = tried_before[i]?boolean_val:0;
				if (untested_bgnd != null){
					dbg_img[11][i] = untested_bgnd[i]? scan.disparity[ty][tx]: Double.NaN;
				}
				if (diff2maxAvg[i] > 0.0){
					dbg_img[14][i] = strength[i]/ diff2maxAvg[i];
				}

			}
			dbg_img[12] = diff2max;
			dbg_img[13] = diff2maxAvg;

			sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "extend_disparity-"+clt_3d_passes.size(),dbg_titles);
		}
		return num_op_tiles; // num_extended;
	}

	public void showScan(
			CLTPass3d   scan,
			String title)
	{
		showDoubleFloatArrays sdfa_instance = null;
		sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		String [] titles = {
				"tile_op",       //  0
				"final",         //  1 - calculated, filtered, combined disparity
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
				"max_tried",     // 12
				"diff0",         // 13
				"diff1",         // 14
				"diff2",         // 15
				"diff3",         // 16
				"diff2max",      // 17
				"diff2maxAvg",   // 18
				"normStrength",  // 19
				"overexp"};      // 20


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
		double [] strength4 = null;
		if (scan.disparity_map != null){
			strength4 = scan.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX];
			dbg_img[ 3] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM];
			dbg_img[ 4] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_HOR];
			dbg_img[ 5] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_VERT];
			dbg_img[ 7] = scan.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX];
			dbg_img[ 8] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH];
			dbg_img[ 9] = scan.disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH];
			dbg_img[20] = scan.disparity_map[ImageDtt.OVEREXPOSED];
		}
		dbg_img[1] = scan.calc_disparity_combo;
		dbg_img[6] = scan.strength;
		double  [][] these_diffs =     scan.getDiffs();
		double [] diff2max = scan.getSecondMaxDiff(false);
		double [] diff2maxAvg = scan.getSecondMaxDiff(true);
		if (these_diffs != null) {
			dbg_img[13] = these_diffs[0];
			dbg_img[14] = these_diffs[1];
			dbg_img[15] = these_diffs[2];
			dbg_img[16] = these_diffs[3];
			dbg_img[17] = diff2max;
			dbg_img[18] = diff2maxAvg;
			dbg_img[19] = new double[tlen];
			if ((diff2maxAvg != null) && (strength4 != null)) {
				for (int i = 0; i < tlen; i++){
					if (diff2maxAvg[i] > 0.0){
						dbg_img[19][i] = strength4[i]/ diff2maxAvg[i];
					}
				}
			}
		}

		title += "-";
		if (scan.isMeasured())  title += "M";
		if (scan.isProcessed()) title += "P";
		if (scan.isCombo())     title += "C";
		sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, title,titles);
		System.out.println("showScan("+title+"): isMeasured()="+scan.isMeasured()+", isProcessed()="+scan.isProcessed()+", isCombo()="+scan.isCombo());

	}



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
		if (debugLevel > -1)
			sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

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
				sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, "tiles_getBackgroundMask");
			}

		}
		return bgnd_tiles;
	}


	//sure_smth
	public boolean [] getBackgroundMask_new( // which tiles do belong to the background
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
		if (debugLevel > -1)
			sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

		//		  boolean [] new_tiles =        new boolean [tilesY * tilesX]; // grow selection by 1 tile over non-background?
		CLTPass3d bgnd_data = clt_3d_passes.get(0);
		//		  double [][][][] texture_tiles = bgnd_data.texture_tiles;
		double [][]     disparity_map=  bgnd_data.disparity_map;
		double [] dbg_worst2 = new double [disparity_map[ImageDtt.IMG_DIFF0_INDEX].length];
		final TileNeibs tnSurface = new TileNeibs(tilesX, tilesY);

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
				double worst1 = 0.0;
				// modified to look around for the 1-st maximum
				for (int i = 1; i< quad; i++){
//					if (disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][tindx] > disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax1][tindx]){
					if (disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][tindx] > worst1){
						imax1 = i;
						worst1 = disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][tindx];
					}
					for (int dir =0; dir < 8; dir++){
						int dtindx = tnSurface.getNeibIndex(tindx, dir);
						if (dtindx >= 0){
							if (disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][dtindx] > worst1){
								imax1 = i;
								worst1 = disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][dtindx];
							}
						}
					}
				}
				int imax2 = (imax1 == 0)? 1 : 0;
				for (int i = 0; i< quad; i++) if (i != imax1) {
					if (disparity_map[ImageDtt.IMG_DIFF0_INDEX+i][tindx] > disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax2][tindx]) imax2 = i;
				}
				dbg_worst2[tindx] = disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax2][tindx];

				block_propagate[tindx] = (disparity_map[ImageDtt.IMG_DIFF0_INDEX + imax2][tindx] > sure_smth);
			}
		}
		// grow block by 1 ?
		tnSurface.growSelection(
				2, // grow,
				block_propagate, // tiles,
				null); // prohibit);
		if ((debugLevel > -1) && show_bgnd_nonbgnd){
			 new showDoubleFloatArrays().showArrays(bgnd_data.disparity_map,  tilesX, tilesY, true, "bgnd_map",ImageDtt.DISPARITY_TITLES);
			 new showDoubleFloatArrays().showArrays(dbg_worst2,  tilesX, tilesY, "worst2");
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
				sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, "tiles_getBackgroundMask_new");
			}

		}
		// was not in tghe original
		tnSurface.growSelection(
				4, // grow,
				bgnd_tiles, // tiles,
				null); // prohibit);
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
//		showDoubleFloatArrays sdfa_instance = null;
//		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?

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
			// disabled below
			final double      ex_min_over,       // when expanding over previously detected (by error) background, disregard far tiles
//			final boolean  [] border_tiles,      // last measured border tiles
			final CLTPass3d   last_meas,         // last measured scan (with border_tiles and getSecondMaxDiff
			final int      [] horVertMod, // +1 - modified by hor correlation, +2 - modified by vert correlation (or null)
			final double      disparity_far,    //
			final double      disparity_near,   //
			final double      ex_strength,      // minimal 4-corr strength to trust tile
			final double      ex_nstrength,     // minimal 4-corr strength divided by channel diff for new (border) tiles
			final double      this_maybe,       // maximal strength to ignore as non-background
			final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
		// using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
			final double [][] plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
			final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
			final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
			final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)

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
//		final double super_trust = 1.6; // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength

		final int dbg_tile = -41030; // x = 206, y = 126 (324*126+206)
		final boolean show_scan = show_filter_scan || (debugLevel > 1);
		showDoubleFloatArrays sdfa_instance = null;
		if ((debugLevel > -2) && ((debugLevel > -1) || show_scan)) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		if (debugLevel > -2){ //-2
			if (debugLevel > -1){ //-2
				System.out.print("FilterScan(,,"+disparity_far+", " +disparity_near+", "+ sure_smth+"... ");
			}
			System.out.print(".");
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
		double [] second_max_diff = null;
		final boolean  [] border_tiles = (last_meas != null) ? last_meas.getBorderTiles() : null;
		if ((ex_nstrength > 0.0) && (border_tiles != null)){
			second_max_diff = last_meas.getSecondMaxDiff(true); // false - single tile, true - 3x3 weighted average
		}

		for (int i = 0; i <tlen; i++) if (!Double.isNaN(this_disparity[i])){
//			if ((debugLevel > 0) && (i==dbg_tile)){
			if ((debugLevel > -2) && (i==dbg_tile)){
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
				// Higher difference, higher the correlation strength is needed
				// assuming this_strength is combine of the strength4 and hor/vert, apply to any
//				if (!bg_tiles[i] || (this_disparity[i] >= ex_min_over)) {
					if ((this_strength[i] > ex_strength) || (horVertMod != null) && (horVertMod[i] != 0)){
						these_tiles[i] = true;
						if (this_strength[i] < ex_strength* super_trust) {
							if (second_max_diff != null) { // only apply to the border tiles, so by next pass they () weak can be used?
								if (border_tiles[i] && (this_strength[i] < second_max_diff[i] * ex_nstrength)) {
									these_tiles[i] = false;
								}
							}
						}
					}
//				}
				// this one is not (yet) compared to normalized difference between the channels
//				if ((horVertMod != null) && (horVertMod[i] != 0)){
//					these_tiles[i] = true;
//				}
			}
			// combine with plates if available
			if ((plate_ds != null) && (plate_ds[1][i] >=0) && (this_strength[i] < ex_strength* super_trust)) { // plate_ds[1][i] == -1 - no data, was not measured
//				if (!bg_tiles[i] || (plate_ds[0][i] >= ex_min_over)) { // prevent expanding over correct background far tiles
					// expanding over (wrong) background is only OK for near tiles
					if (plate_ds[1][i] > 0){
						if ( // !refined_selected[i] ||
								!these_tiles[i] ||
								(this_disparity[i] <= plate_ds[0][i]) ||  // keep old way detected (including poles) is closer, replace if raw is farther
								(this_strength[i] < plate_ds[1][i] * scale_filtered_strength_pre) // keep old way detected (including poles) is closer, replace if raw is farther
								) {
							these_tiles[i] = true;
							this_disparity[i] = plate_ds[0][i];
							this_strength[i] = plate_ds[1][i];
						}
					} else if (these_tiles[i]){ // no plates and raw is not reliably
						these_tiles[i] = false;
					}
//				}
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

//ex_nstrength


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

			String [] titles = {"masked","map","orig_map","hor_map","vert_map","bg_sel","far", // 0-6
					"before_small","before_lone","before_gaps","these","near","block", // 7-12
					"strength","hor-strength","vert-strength", //13-15
					"diff0","diff1","diff2","diff3", "enum_clusters", "disp_cm", "disp_poly", "disp_hor", "disp_vert", "poison","plate_d","plate_s"};
			double [][] dbg_img = new double[titles.length][tilesY * tilesX];
			for (int i = 0; i<dbg_img[0].length;i++){
				dbg_img[ 0][i] =   this_disparity_masked[i];
				dbg_img[ 1][i] =   this_disparity[i];
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
			if (plate_ds != null) {
				dbg_img[26] =  plate_ds[0]; // ? same
				dbg_img[26] =  plate_ds[1]; // ? same
			}

			sdfa_instance.showArrays(dbg_img,  tilesX, tilesY, true, "FilterScan"+clt_3d_passes.size(),titles);
			System.out.println("FilterScan"+clt_3d_passes.size());
		}

		return these_tiles;
	}

// Next runs before filter, maybe should be after or combined otherwise (i.e filter strengths, but keep original disparities. Or use disparity
// to some power when averaging here)
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
	  		final double      or_maxDisp,   // 5.0;   // Maximal disparity to apply ortho correction
	  		final boolean     show_ortho,   // show replacement of disparity/strength by those of hor/vert
			final int         debugLevel)
	{
		double [] disparity =      scan.getDisparity(0);  // calculated, to be modified
		double [] strength =       scan.getStrength();    // combo, to be modified
		double [] disparity_hor =  scan.getDisparity(2); // .clone(); // calculated
		double [] disparity_vert = scan.getDisparity(3); // .clone(); // calculated
		double [] strength_hor =   scan.getHorStrength(); // .clone();
		double [] strength_vert =  scan.getVertStrength(); // .clone();
		int [] replaced = new int[strength.length];
		double [][] dbg_img = null;
		String []  dbg_titles = {"disparity",         // 0
								 "orig_disparity",    // 1
								 "hor_disparity",     // 2
								 "hor_disparity_mod", // 3
								 "vert_disparity",    // 4
								 "vert_disparity_mod",// 5
								 "strength",          // 6
								 "orig_strength",     // 7
								 "hor_strength",      // 8
								 "hor_strength_mod",  // 9
								 "vert_strength",     // 10
								 "vert_strength_mod", // 11
								 "replaced"};         // 12

		if (show_ortho) {
			dbg_img=new double [dbg_titles.length][];
			dbg_img[1] =  disparity.clone();
			dbg_img[2] =  disparity_hor.clone();
			dbg_img[4] =  disparity_vert.clone();
			dbg_img[7] =  strength.clone();
			dbg_img[8] =  strength_hor.clone();
			dbg_img[10] = strength_vert.clone();

		}
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
			if (show_ortho) {
				dbg_img[3] =  disparity_hor.clone();
				dbg_img[9] =  strength_hor.clone();
			}

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
			if (show_ortho) {
				dbg_img[5] =  disparity_vert.clone();
				dbg_img[11] =  strength_vert.clone();
			}
		}
		double ko = (or_threshold - 1) * or_offset;
		double ao = (or_asym - 1) * or_offset;
		if (or_hor){
			for (int i = 0; i < strength_hor.length; i++){
				if (    (disparity_hor[i] < or_maxDisp) &&
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
				if (    (disparity_vert[i] < or_maxDisp) &&
						(strength_vert[i] > or_absVert) &&
						(strength_vert[i] > or_threshold * strength[i] - ko) &&
						(!or_hor || (strength_vert[i] > or_asym * strength_hor[i] - ao))){
					strength[i] =  strength_vert[i];
					disparity[i] = disparity_vert[i];
					replaced[i] |= 2;
				}
			}
		}
		if (show_ortho) {
			dbg_img[0] =  disparity.clone();
			dbg_img[6] =  strength.clone();
			dbg_img[12] = new double [replaced.length];
			for (int i = 0; i < replaced.length; i++) {
				dbg_img[12][i] = replaced[i];
			}
			(new showDoubleFloatArrays()).showArrays(dbg_img, tilesX, tilesY, true, "ortho_internal",dbg_titles);
		}
		return replaced;
	}

	public boolean [] combineHorVertDisparity_old(
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

		boolean use_ortho = true; // false;
		if (use_ortho) {
			for (int i = 0; i < tilesY; i++){
				for (int j = 0; j < tilesX; j++){
					int indx = i * tilesX + j;
					if (!Double.isNaN(this_disparity[i])){
						double disp = this_disparity[indx];
						// Enhance foreground detection: compare horizontal-only (for vertical features) and vertical (for horizontal) and replace 4-pair disparity
						// with that value
						if ((hor_strength_conv[indx] >= clt_parameters.ortho_min_hor) &&
								(hor_strength_conv[indx] / vert_strength_conv[indx]>= clt_parameters.ortho_asym) &&
								((strength[indx] ==0.0) || (hor_strength_conv[indx] / strength[indx] >= clt_parameters.ortho_over4)) &&
								(this_hor_disparity[indx] > disp)) {
							disp = this_hor_disparity[indx];
							used_hor[indx] = true;
						}
						if ((vert_strength_conv[indx] >= clt_parameters.ortho_min_vert) &&
								(vert_strength_conv[indx] / hor_strength_conv[indx]>= clt_parameters.ortho_asym) &&
								((strength[indx] ==0.0) || (vert_strength_conv[indx] / strength[indx] >= clt_parameters.ortho_over4)) &&
								(this_vert_disparity[indx] > disp)) {
							disp = this_vert_disparity[indx];
							used_vert[indx] = true;
						}
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

			// TODO: check - now no limit on the strength of the offending selections, only on these ones

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
	public CLTPass3d getLastMeasured(
			int indx){
		if (indx < 0) {
			indx = clt_3d_passes.size() -1;
		}
		for (int i = indx; i >= 0; i--){
			if (clt_3d_passes.get(i).isMeasured()){
				return clt_3d_passes.get(i);
			}
		}
		return null;
	}


	public CLTPass3d refinePassSetup( // prepare tile tasks for the second pass based on the previous one(s)
			//			  final double [][][]       image_data, // first index - number of image in a quad
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			// disparity range - differences from
			boolean           use_supertiles,   // false
			int               bg_scan_index,    // 0
			double            disparity_far,    // 0.3
			double            disparity_near,   // 150
			double            ex_strength,      // 0.18  minimal 4-corr strength to trust raw tile - above will prefer it to multi-tile plates
	  		double            ex_nstrength,     // 0.4 minimal 4-corr strength divided by channel diff for new (border) tiles
			double            this_maybe,       // 0.1 maximal strength to ignore as non-background
			double            sure_smth,        // 2.0 if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd

			double            super_trust,      // 1.6 If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
			double [][]       plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
			boolean           keep_raw_fg,      // true do not replace raw tiles by the plates, if raw is closer (like poles)
			double            scale_filtered_strength_pre, // 1.5  scale plate_ds[1] before comparing to raw strength
			double            scale_filtered_strength_post,// 2.5 scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)

			int               disparity_index,  // index of disparity value in disparity_map == 2 (0,2 or 4)
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1);
		CLTPass3d scan_bg =   clt_3d_passes.get(bg_scan_index); //
		CLTPass3d scan_lm =   getLastMeasured(-1);
//		boolean [] border_tiles = (scan_lm != null) ? scan_lm.getBorderTiles() : null;
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		//TODO: for next passes - combine all selected for previous passes (all passes with smaller disparity)
		int [] replaced =  null; // +1 - hor, +2 - vert
		int [] replaced0 = null; // +1 - hor, +2 - vert
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
				clt_parameters.or_maxDisp,   // 5.0;   // Maximal disparity to apply ortho correction
				show_ortho,                  //  show replacement of disparity/strength by those of hor/vert
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
			if (debugLevel > 0){
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

		if (debugLevel > 10){
			showScan( scan_prev, "preFilter_last");
			showScan( scan_lm, "preFilter_meas");
		}



		boolean [] these_tiles = FilterScan(
				scan_prev,        // final CLTPass3d   scan,
				scan_bg.selected, // get from selected in clt_3d_passes.get(0);
				clt_parameters.ex_min_over,// when expanding over previously detected (by error) background, disregard far tiles
//				border_tiles,      // final boolean  [] border_tiles,      // last measured boirder tiles
				scan_lm,          // final CLTPass3d   last_meas,         // last measured scan (with border_tiles and getSecondMaxDiff
				replaced,         // final int      [] horVertMod, // +1 - modified by hor correlation, +2 - modified by vert correlation (or null)
				disparity_far,    // final double      disparity_far,    //
				disparity_near,   // final double      disparity_near,   //
				ex_strength,      // final double      ex_strength,   // minimal 4-corr strength to trust tile
				ex_nstrength,     // final double      ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
				this_maybe,       // final double      this_maybe,       // maximal strength to ignore as non-background
				sure_smth,        // final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				super_trust,      // final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
				// using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
				plate_ds,        // final double [][] plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
				keep_raw_fg,  // final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
				scale_filtered_strength_pre, // 					final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
				scale_filtered_strength_post,// final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)
				clt_parameters.show_filter_scan,
				clt_parameters.min_clstr_seed, // clt_parameters.min_clstr_seed
				clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
				clt_parameters.min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
				clt_parameters.min_clstr_lone,  // int        min_area,  // minimal number of pixels
				clt_parameters.fill_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
				2, // 0, // final int     poison_gaps,      // Do not fill gaps that have even single "poisoned" tile
				false,             // final boolean     zero_gap_strength, // set strength to zero when covering gaps
				debugLevel);


		//************************************************
		// Show supertiles histograms

		//		if (clt_parameters.stShow){

		// try renovated supertiles. Do twice to show both original and blured histograms
		double [] dbg_orig_disparity =  null;
		double [] dbg_with_super_disp = null;
		double [] dbg_outlayers =       null;

		if (use_supertiles) { // *** currently is not used (obsolete ***
			String [] dbg_st_titles = {"raw", "blurred"+clt_parameters.stSigma,"max-min-max"};
			double [][] dbg_hist = new double[dbg_st_titles.length][];

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					0.0,// NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,        // filter parameters
//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples
//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.getSuperTiles().showDisparityHistogram();

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					clt_parameters.stSigma, // with blur double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,        // filter parameters
//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples
//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.getSuperTiles().showDisparityHistogram();

			dbg_hist[2] = scan_prev.getSuperTiles().showMaxMinMax();

			int hist_width0 =  scan_prev.getSuperTiles().showDisparityHistogramWidth();
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
		if (show_super && false){ // something is broken, java.lang.NullPointerException at TileProcessor.refinePassSetup(TileProcessor.java:4488)
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
//		clt_3d_passes.add(scan_next); // was here
		//		}
		return scan_next;
	}


	public boolean assignTilesToSurfaces(
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch
			final boolean     updateStatus,
//			final boolean     batch_mode,
			final int         debugLevel)
	{
		final boolean    batch_mode = clt_parameters.batch_run;
		final int debugLevelInner =  batch_mode ? -5: debugLevel;
		trimCLTPasses(); // make possible to run this method multiple times - remove extra passes added by it last time
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
//		boolean show_st =    clt_parameters.stShow || (debugLevel > 1);
		SuperTiles st = scan_prev.getSuperTiles();
		TileSurface tileSurface = st.getTileSurface();
		if (tileSurface == null){
			return false;
		}
		// show testure_tiles

		double [][][][] texture_tiles = scan_prev.getTextureTiles();
		ImageDtt image_dtt = new ImageDtt();

		double [][][]  dispStrength = st.getDisparityStrengths(
				clt_parameters.stMeasSel); // int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)

		boolean [][] tileSel =  st.getMeasurementSelections(
				clt_parameters.stMeasSel); // int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)


		if (texture_tiles != null){
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			double [][] texture_nonoverlap = null;
			double [][] texture_overlap = null;
			String [] rgba_titles = {"red","blue","green","alpha"};
			String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
			String name = scan_prev.getTextureName();

			boolean show_nonoverlap = false; // true; // clt_parameters.show_nonoverlap
			boolean show_overlap =    false; //true; // clt_parameters.show_overlap
			boolean show_rgba_color = false; //true; // clt_parameters.show_rgba_color

			if (!batch_mode && show_nonoverlap){
				texture_nonoverlap = image_dtt.combineRGBATiles(
						texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						clt_parameters.transform_size,
						false,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				sdfa_instance.showArrays(
						texture_nonoverlap,
						tilesX * (2 * clt_parameters.transform_size),
						tilesY * (2 * clt_parameters.transform_size),
						true,
						name + "-TXTNOL-D",
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

			}

			if (!batch_mode && (show_overlap || show_rgba_color)){
				int alpha_index = 3;
				texture_overlap = image_dtt.combineRGBATiles(
						texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
						clt_parameters.transform_size,
						true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
						clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						threadsMax,                    // maximal number of threads to launch
						debugLevel);
				if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
					double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
					for (int i = 0; i < texture_overlap[alpha_index].length; i++){
						double d = texture_overlap[alpha_index][i];
						if      (d >=clt_parameters.alpha1) d = 1.0;
						else if (d <=clt_parameters.alpha0) d = 0.0;
						else d = scale * (d- clt_parameters.alpha0);
						texture_overlap[alpha_index][i] = d;
					}
				}

				if (show_overlap) {
					sdfa_instance.showArrays(
							texture_overlap,
							tilesX * clt_parameters.transform_size,
							tilesY * clt_parameters.transform_size,
							true,
							name + "-TXTOL-D",
							(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
				}
			}
			if (!batch_mode && (debugLevel > -1)) {
				double [][] tiles_tone = scan_prev.getTileRBGA(
						12); // int num_layers);
				sdfa_instance.showArrays(
						tiles_tone,
						tilesX,
						tilesY,
						true,
						name + "tiles_tone",
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));
			}
		}
		TileAssignment ta = new TileAssignment(
				clt_parameters, // EyesisCorrectionParameters.CLTParameters clt_parameters,
				tileSurface,    // TileSurface ts,
				scan_prev,      // CLTPass3d p3d,
				//					tileSel,        // boolean[][] tile_sel,
				dispStrength, // double [][][] dispStrength,
				Double.NaN,     //  double kR,
				Double.NaN,     // double kB,
				Double.NaN);    // double fatZero)
		if (!batch_mode && (debugLevel > -1)) {
			ta.showToneDiffWeights3("Raw_");
			ta.showToneDiffWeights1("Raw_");
		}
		ta.blurMixTones(
				false, // final boolean use_sqrt,
				true); //final boolean weighted) // blur weighted
		if (!batch_mode && (debugLevel >-1)) {
			ta.showToneDiffWeights3("Mixed_");
			ta.showToneDiffWeights1("Mixed_");
		}
		//		}






		// end of show testure_tiles
/*
		tileSurface.testSimpleConnected(
				230, // clt_parameters.tileX,
				131);//clt_parameters.tileY);
*/



		// Reset/initialize assignments - if not done so yet or specifically requested
		boolean first_run = !tileSurface.isInit() || clt_parameters.tsReset;
		tileSurface.InitTilesAssignment(
				clt_parameters.tsReset,
				dispStrength, // final double [][][]                            dispStrength,
				tileSel,      // final boolean [][]                             tileSel,
                debugLevelInner);    // final int                                    debugLevel,
        // assign tiles that do not depend on other assigned tiles - single pass
		int [][] tile_layers = 	             tileSurface.getTileLayersCopy();
		int [][] tile_layers_single_surf =   tileSurface.newTileLayers(tileSel); //  final boolean [][]                             tileSel,)
		int [][] tile_layers_planes =        tileSurface.newTileLayers(tileSel); //  final boolean [][]                             tileSel,)

		int [] stats_single_surf=  tileSurface.assignTilesToSingleSurface(
				tile_layers_single_surf,                      //final int [][]      tileLayers,
				clt_parameters.tsNoEdge ,       // final boolean       noEdge,
				batch_mode ? -5: 0, // -1,                       // debugLevel,                  // final int        debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
		System.out.print("Assign the only surface :");
		tileSurface.printStats(stats_single_surf);
		// grow the only surface assignments
		if (clt_parameters.tsEnGrow){
			double [] growStrengths =     new double [tile_layers_single_surf.length];
			double [] contMinStrengths =  new double [tile_layers_single_surf.length];
			double [] growMaxDiffFar =    new double [tile_layers_single_surf.length];
			double [] growMaxDiffNear =   new double [tile_layers_single_surf.length];

			int [][] assign_conflicts = null;
			for (int i = 0; i < growStrengths.length; i++){
				growStrengths[i] =    clt_parameters.tsGrowStrength; // save all the same, later can use different for different types of correlation
				contMinStrengths[i] = clt_parameters.tsContStrength;
				growMaxDiffFar[i]=    clt_parameters.tsContDiff;
				growMaxDiffNear[i]=   clt_parameters.tsContDiff;
				if ((tile_layers_single_surf[i] != null) && (assign_conflicts == null)){
					assign_conflicts = new int[tile_layers_single_surf[i].length][];
				}
			}
			stats_single_surf = tileSurface.growWeakAssigned(
					tile_layers_single_surf,                                           //final int [][]      tileLayers,
					assign_conflicts,                                      // final int [][]      conflicts,
					clt_parameters.tsNoEdge ,                              // final boolean       noEdge,
					growStrengths,                                         // final double []     maxStrength,
					(clt_parameters.tsGrowStrong? contMinStrengths: null), //final double []     minStrengthContinue,
					(clt_parameters.tsEnGrow?     growMaxDiffFar: null),   // final double []     maxDiffFar, // null
					(clt_parameters.tsEnGrow?     growMaxDiffNear: null),  // final double []     maxDiffNear, // null
					clt_parameters.plDispNorm,                             // final double        dispNorm, // disparity normalize (proportionally scale down disparity difference if above
					dispStrength,                                          // final double [][][] dispStrength,
					debugLevelInner, // 2, // -1,                                              // debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			System.out.print("growWeakAssigned():");
			tileSurface.printStats(stats_single_surf);
		}

		int [] stats_planes = tileSurface.assignPlanesTiles(
				true,                // final boolean                  force,
				tile_layers_planes,         //final int [][]      tileLayers,
				st.planes_mod,       // final TilePlanes.PlaneData[][] planes
				clt_parameters.tsNoEdge ,                              // final boolean       noEdge,
				debugLevelInner, //  2, // -1,            // debugLevel,                  // final int        debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
		System.out.print("assignPlanesTiles():");
		tileSurface.printStats(stats_planes);


		if (clt_parameters.tsEnOnly){
			tileSurface.combineTileLayers(
				true, // final boolean overwrite,
				tile_layers, // final int [][] dst,
				tile_layers_single_surf); // final int [][] src);
		}

		if (clt_parameters.tsEnPlaneSeed) {
				tileSurface.combineTileLayers(
					true, // final boolean overwrite,
					tile_layers, // final int [][] dst,
					tile_layers_planes); // final int [][] src);
		}


		if (clt_parameters.tsEnSingle) {
			int [] stats=  tileSurface.assignTilesToSurfaces(
					tile_layers,                    //final int [][]      tileLayers,
					clt_parameters.tsNoEdge ,       // final boolean       noEdge,
					clt_parameters.tsUseCenter,     // final boolean       useCenter,
					clt_parameters.tsMaxDiff,       //final double        maxDiff,
					clt_parameters.tsMinDiffOther,  //final double        minDiffOther, // should be >= maxDiff
					clt_parameters.tsMinStrength,   //final double        minStrength,
					clt_parameters.tsMaxStrength,   //final double        maxStrength,
					clt_parameters.tsMinSurface,    //final double        minSurfStrength, // minimal surface strength at the tile location
					clt_parameters.tsMoveDirs,      //final int           moveDirs, // 1 increase disparity, 2 - decrease disparity, 3 - both directions
					false,                          // clt_parameters.tsEnMulti,       //final boolean       enMulti,
					clt_parameters.tsSurfStrPow,    //final double        surfStrPow, // surface strength power
					clt_parameters.tsAddStrength,   //final double        addStrength,
					clt_parameters.tsSigma,         //final double        sigma,
					clt_parameters.tsNSigma,        //final double        nSigma,
					clt_parameters.tsMinPull,       //final double        minPull,
					clt_parameters.tsMinAdvantage,  //final double        minAdvantage,
					clt_parameters.plDispNorm,      // final double        dispNorm, // disparity normalize (proportionally scale down disparity difference if above
					dispStrength,                   // final double [][][] dispStrength,
					batch_mode ? -5: 0,   // -1,    // debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			System.out.print("Assign to nearest with single candidate :");
			tileSurface.printStats(stats);
		}

		// assign tiles that depend on other assigned tiles (more neighbors on the same surface - more attractive it is)
		// Single pass or until more tiles are assigned
		if (clt_parameters.tsEnMulti) {
			for (int nTry = 0; nTry < 100; nTry++) {
				int [] stats=  tileSurface.assignTilesToSurfaces(
						tile_layers,                      //final int [][]      tileLayers,
						clt_parameters.tsNoEdge ,       // final boolean       noEdge,
						clt_parameters.tsUseCenter,     // final boolean       useCenter,
						clt_parameters.tsMaxDiff,       //final double        maxDiff,
						clt_parameters.tsMinDiffOther,  //final double        minDiffOther, // should be >= maxDiff
						clt_parameters.tsMinStrength,   //final double        minStrength,
						clt_parameters.tsMaxStrength,   //final double        maxStrength,
						clt_parameters.tsMinSurface,    //final double        minSurfStrength, // minimal surface strength at the tile location
						clt_parameters.tsMoveDirs,      //final int           moveDirs, // 1 increase disparity, 2 - decrease disparity, 3 - both directions
						true, // clt_parameters.tsEnMulti,       //final boolean       enMulti,
						clt_parameters.tsSurfStrPow,    //final double        surfStrPow, // surface strength power
						clt_parameters.tsAddStrength,   //final double        addStrength,
						clt_parameters.tsSigma,         //final double        sigma,
						clt_parameters.tsNSigma,        //final double        nSigma,
						clt_parameters.tsMinPull,       //final double        minPull,
						clt_parameters.tsMinAdvantage,  //final double        minAdvantage,
						clt_parameters.plDispNorm,      // final double        dispNorm, // disparity normalize (proportionally scale down disparity difference if above
						dispStrength,                   // final double [][][] dispStrength,
						batch_mode ? -5: 0,             // 0, // -1,  final int        debugLevel)
						clt_parameters.tileX,
						clt_parameters.tileY);
				System.out.print("Assign to best fit surface :");
				tileSurface.printStats(stats);
				if (!clt_parameters.tsLoopMulti || !tileSurface.makesSensToTry(stats)){
					break;
				}
			}
		}

		// Remove weak clusters that have too few connected tiles or the total weight of the cluster is insufficient.
		// Single pass, before growing assigned tiles
		if (clt_parameters.tsRemoveWeak1) {
			int [] stats=  tileSurface.removeSmallClusters(
					tile_layers,                      //final int [][]      tileLayers,
					clt_parameters.tsClustSize , // final int           minSize,      // **
					clt_parameters. tsClustWeight, // final double        minStrength,  // **
					dispStrength,                   // final double [][][] dispStrength,
					batch_mode ? -5: 0,             // 0, // -1,  final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			tileSurface.printStats(stats);
		}

		// assign tiles that have sufficient neighbors around, run single or multiple times (TODO: Split parameters?)
		// new tile is assigned one of the neighbor surfaces, that has the lowest disparity (farthest from the camera)
		// in the current location
		if (clt_parameters.tsGrowSurround) {
			for (int nTry = 0; nTry < 100; nTry++) {
				int [] stats=  tileSurface.assignFromFarthest(
						tile_layers,                      //final int [][]      tileLayers,
						clt_parameters.tsNoEdge ,         // final boolean       noEdge,
						clt_parameters.tsMinNeib ,        // final int           minNeib,           // **
						clt_parameters.tsMaxSurStrength , // final double        maxStrength,       // **
						clt_parameters.tsCountDis,        // final boolean       includeImpossible, // ** // count prohibited neighbors as assigned
						dispStrength,                     // final double [][][] dispStrength,
						batch_mode ? -5: 0,             // 0, // -1,  final int        debugLevel)
						clt_parameters.tileX,
						clt_parameters.tileY);
				System.out.print("Assign from neighbors :");
				tileSurface.printStats(stats);
				if (!clt_parameters.tsLoopMulti || (tileSurface.newAssigned(stats) ==0)){
					break;
				}
			}
		}




		// Remove weak clusters that have too few connected tiles or the total weight of the cluster is insufficient.
		// Single pass, after growing assigned tiles
		if (clt_parameters.tsRemoveWeak2) {
			int [] stats=  tileSurface.removeSmallClusters(
					tile_layers,                      //final int [][]      tileLayers,
					clt_parameters.tsClustSize , // final int           minSize,      // **
					clt_parameters. tsClustWeight, // final double        minStrength,  // **
					dispStrength,                   // final double [][][] dispStrength,
					batch_mode ? -5: 0,             // 0, // -1,  final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			tileSurface.printStats(stats);
		}

		// Just show current state
		tileSurface.statTileLayers(
				tile_layers, // final int [][]      tileLayers,
				tileSel,      // final boolean [][]       tileSel,
				1); // final int                debugLevel)

// removed tile_layers from globals, so multiple parallel assignments can be generated and then combined

//		int [][] tile_layers = 	             tileSurface.getTileLayersCopy();
//		int [][] tile_layers_single_surf =   tileSurface.newTileLayers(tileSel); //  final boolean [][]                             tileSel,)
//		int [][] tile_layers_planes =        tileSurface.newTileLayers(tileSel); //  final boolean [][]                             tileSel,)
		int [][][] tile_assignments = {tile_layers, tile_layers_single_surf, tile_layers_planes};
		for (int i = 0; i < tile_assignments.length; i++){
			if ((clt_parameters.tsConsensMode & (1 << i)) == 0) {
				tile_assignments[i] = null;
			}
		}
		tileSurface.compareAssignments(
				tile_assignments, // final int [][][] tileAssignments)
				debugLevelInner); //

		int [][][] opinions = new int [tile_layers.length][][];
		tile_layers = tileSurface.getConsensusAssignment(
				clt_parameters.tsConsensAgree, //  final int        min_agree,
				opinions, // int [][][]       opinions_in,
				tile_assignments); // final int [][][] tileAssignments)

		int [][] tile_layers_surf = ta.imgToSurf(tile_layers);
		if (!batch_mode && (debugLevel > -1)) {
			double [][] dbg_tls = new double [tile_layers_surf.length][];
			for (int ml = 0; ml < tile_layers_surf.length; ml++) if (tile_layers_surf[ml] != null){
				dbg_tls[ml] = new double [tile_layers_surf[ml].length];
				for (int i = 0; i < tile_layers_surf[ml].length; i++){
					dbg_tls[ml][i] = tile_layers_surf[ml][i];
				}
			}
			(new showDoubleFloatArrays()).showArrays(dbg_tls, ta.getSurfTilesX(), ta.getSurfTilesY(), true, "tile_layers_surf");
		}
		if (!batch_mode && (debugLevel > -1)) {
			ta.showTileCost("before_",tile_layers_surf);
			ta.showTileCosts("before_",tile_layers_surf);
		}
		TileAssignment.TACosts [] ta_stats = ta.statTileCosts(tile_layers_surf);

		for (int i = 0; i < ta_stats.length; i++){
			System.out.println(ta_stats[i].toString());
		}

		ta.optimizeAssignment25(
				clt_parameters.tsNoEdge ,         // final boolean       noEdge,
				tile_layers_surf, // final int [][]  tileLayers,
				(batch_mode ? -5 : 2),  //	2, // final int       debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);


		if (!batch_mode && (debugLevel > -1)) {
			double [][] dbg_tls = new double [tile_layers_surf.length][];
			for (int ml = 0; ml < tile_layers_surf.length; ml++) if (tile_layers_surf[ml] != null){
				dbg_tls[ml] = new double [tile_layers_surf[ml].length];
				for (int i = 0; i < tile_layers_surf[ml].length; i++){
					dbg_tls[ml][i] = tile_layers_surf[ml][i];
				}
			}
			(new showDoubleFloatArrays()).showArrays(dbg_tls, ta.getSurfTilesX(), ta.getSurfTilesY(), true, "optimized_tile_layers_surf");
		}
		if (!batch_mode && (debugLevel > -2)) {
			ta.showTileCost("after_",tile_layers_surf);
		}
		if (!batch_mode && (debugLevel > -1)) {
			ta.showTileCosts("after_",tile_layers_surf);
		}
		ta_stats = ta.statTileCosts(tile_layers_surf);
		System.out.println("Optimized:");
		for (int i = 0; i < ta_stats.length; i++){
			System.out.println(ta_stats[i].toString());
		}

		tile_layers = ta.surfToImg(tile_layers_surf);

//==============
		tileSurface.setTileLayers(tile_layers);



		// Just show current state
		tileSurface.InitTilesAssignment(
				false,
				dispStrength, // final double [][][]                            dispStrength,
				tileSel,      // final boolean [][]                             tileSel,
                debugLevel);    // final int                                    debugLevel,

		if (!batch_mode && clt_parameters.tsShow && (debugLevel > -2)){
			tileSurface.showAssignment(
					"assignments", // String title,
					dispStrength); // final double [][][] dispStrength)
		}
		boolean [][] assigned_sel = tileSurface.extractSelection(
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);


		int [][] clusters1 = tileSurface.enumerateClusters(
				assigned_sel, //final boolean [][] selection,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		int [][] cluster_stats1 = tileSurface.clusterStats(
				clusters1, // int [][] clusters,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		tileSurface.showClusterStats(
				cluster_stats1, // int [][] cluster_stats,
				-1); //clt_parameters.tsNumClust); // int max_clusters){

		if (!batch_mode && clt_parameters.tsShow && (debugLevel > -1)){
			tileSurface.showClusters(
					"clusters_individual", // String title,
					cluster_stats1, // int [][] cluster_stats,
					clt_parameters.tsNumClust, // int max_clusters
					clusters1); // int [][] clusters); // final double [][][] dispStrength)
		}

		// Try splitting (currently no initial conflicts):
		int [][] clusters1a = tileSurface.spitConflictClusters(
				clusters1, // final int [][] clusters,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);
		int [][] cluster_stats1a = tileSurface.clusterStats(
				clusters1a, // int [][] clusters,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		tileSurface.showClusterStats(
				cluster_stats1a, // int [][] cluster_stats,
				-1); //clt_parameters.tsNumClust); // int max_clusters){

		if (!batch_mode && clt_parameters.tsShow && (debugLevel > -1)){
			tileSurface.showClusters(
					"clusters_individual_fixed", // String title,
					cluster_stats1a, // int [][] cluster_stats,
					clt_parameters.tsNumClust, // int max_clusters
					clusters1a); // int [][] clusters); // final double [][][] dispStrength)
		}



		boolean [][] grown_sel = tileSurface.growSelection(
				2,              // int grow,
				assigned_sel,   // final boolean [][] sel_in,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		int [][] clusters2 = tileSurface.enumerateClusters(
				grown_sel,      // final boolean [][] selection,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		int [][] cluster_stats2 = tileSurface.clusterStats(
				clusters2, // int [][] clusters,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		tileSurface.showClusterStats(
				cluster_stats2, // int [][] cluster_stats,
				-1); //clt_parameters.tsNumClust); // int max_clusters){

		if (!batch_mode && clt_parameters.tsShow && (debugLevel > -1)){
			tileSurface.showClusters(
					"clusters_merged", // String title,
					cluster_stats2, // int [][] cluster_stats,
					clt_parameters.tsNumClust, // int max_clusters
					clusters2); // int [][] clusters); // final double [][][] dispStrength)
		}

		// Just for testing: splitting combined clusters
		int [][] clusters2a = tileSurface.spitConflictClusters(
				clusters2, // final int [][] clusters,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		int [][] cluster_stats2a = tileSurface.clusterStats(
				clusters2a, // int [][] clusters,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		tileSurface.showClusterStats(
				cluster_stats2a, // int [][] cluster_stats,
				-1); //clt_parameters.tsNumClust); // int max_clusters){

		if (!batch_mode && clt_parameters.tsShow && (debugLevel > -1)){
			tileSurface.showClusters(
					"clusters_merged_split", // String title,
					cluster_stats2a, // int [][] cluster_stats,
					clt_parameters.tsNumClust, // int max_clusters
					clusters2a); // int [][] clusters); // final double [][][] dispStrength)
		}

		tileSurface.mergeAndGrow( // T ODO: add result
				assigned_sel, // final boolean [][] sel_in,
				batch_mode ? -5: 0, // final int           debugLevel,
				clt_parameters.tileX,
				clt_parameters.tileY);

		if (!batch_mode && clt_parameters.tsShow && (debugLevel > -1)){
			int numToShow = clt_parameters.tsNumClust;
			int num_surf = tileSurface.getSurfaceDataLength();
			if (numToShow > num_surf) numToShow = num_surf;
			String [] titles = new String[2*numToShow];
			double [][] img_data = new double [2*numToShow][];
			for (int i = 0; i < numToShow; i++){
				titles[2 * i    ] = "disp_"+i;
				titles[2 * i + 1] = "pure_disp_"+i;
				img_data[2 * i    ] = tileSurface.getSurfaceData(i).getDisparityNaN(true);
				img_data[2 * i + 1] = tileSurface.getSurfaceData(i).getDisparityNaN(false);

			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(img_data, tilesX, tilesY, true, "final_surfaces",titles);
		}
		return true;
	}

//======================
	public void showPlanes(
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch
			final boolean     updateStatus,
			final int         debugLevel)
	{
		final boolean                       batch_mode = clt_parameters.batch_run; //disable any debug images
		trimCLTPasses(); // make possible to run this method multiple time - remove extra passes added by it last time
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one

		boolean show_st =    clt_parameters.stShow || (debugLevel > 1);
		// recalculate supertiles (may be removed later)
		//		if (use_supertiles || show_st) {
		String [] dbg_st_titles = {"raw", "sampled", "blurred"+clt_parameters.stSigma,"max-min-max"};
		double [][] dbg_hist = new double[dbg_st_titles.length][];
		if (show_st) { // otherwise only blurred version is needed
			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters
//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples
//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.getSuperTiles().showDisparityHistogram();
			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters
//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples
//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.getSuperTiles().showDisparityHistogram();
		}

		//		SuperTiles st =
		scan_prev.setSuperTiles(
				clt_parameters.stStepNear,       // double     step_disparity,
				clt_parameters.stStepFar,        // double     step_near,
				clt_parameters.stStepThreshold,  // double     step_threshold,
				clt_parameters.stMinDisparity,   // double     min_disparity,
				clt_parameters.grow_disp_max,   // double     max_disparity,
//				clt_parameters.stFloor,          // double     strength_floor,
//				clt_parameters.stPow,            // double     strength_pow,
				clt_parameters.stSigma,          // with blur double     stBlurSigma)
				false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
				clt_parameters.mlfp,         // Filter parameters
//				clt_parameters.stSmplSide,  // Sample size (side of a square)
//				clt_parameters.stSmplNum,   // Number after removing worst
//				clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//				clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples
//				clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//				clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//				clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//				clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//				clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//				clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
				clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
		if (show_st) { // otherwise only blured version is needed
			dbg_hist[2] = scan_prev.getSuperTiles().showDisparityHistogram();
			dbg_hist[3] = scan_prev.getSuperTiles().showMaxMinMax();
		}


		if (show_st){
			int hist_width0 =  scan_prev.getSuperTiles().showDisparityHistogramWidth();
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
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters
//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples
//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
		}

		SuperTiles st = scan_prev.getSuperTiles();
		st.setTileSurface(geometryCorrection); // tileSurface);

		double []  world_hor = {0.0, 1.0, 0.0};
		st.processPlanes5(
				clt_parameters.stGrowSel,       // = 2; // = -1;  //Grow initial selection before processing supertiles, odd - ortho. <0 - use all tiles
				clt_parameters.stMeasSel,       //      =     1   //Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				clt_parameters.plDispNorm,      //      =   2.0;  // Normalize disparities to the average if above
				clt_parameters.plMinPoints,     //      =     5;  // Minimal number of points for plane detection
				clt_parameters.plTargetEigen,   //      =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
				clt_parameters.plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
				clt_parameters.plMaxOutliers,   //      =    20;  // Maximal number of outliers to remove\
				clt_parameters.plPreferDisparity,
				geometryCorrection,
				clt_parameters.correct_distortions,

				clt_parameters.stSmplMode ,     // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				clt_parameters.mlfp,         // Filter parameters
//				clt_parameters.stSmplSide ,     // final int                        smplSide, //        = 2;      // Sample size (side of a square)
//				clt_parameters.stSmplNum ,      // final int                        smplNum, //         = 3;      // Number after removing worst
//				clt_parameters.stSmplRms ,      // final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples

//				clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//				clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//				clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//				clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//				clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//				clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

				clt_parameters.plBlurBinHor,    // final double     bin_blur_hor,   // Blur disparity histograms for horizontal clusters by this sigma (in bins)
				clt_parameters.plBlurBinVert,   // final double     bin_blur_vert,  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
				clt_parameters.plMaxDiffHor,    // final double     max_diff_hor,   // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for horizontal planes
				clt_parameters.plMaxDiffVert,   // final double     max_diff_vert,  // maximal disparity difference (to assign to a cluster (of Double.NaN) at first run for vertical plane
				clt_parameters.plInitPasses,    // final int        max_tries,       // on last run - assign all remaining pixels to some cluster (disregard max_diff)

				clt_parameters.msUseSel,        // final boolean    msUseSel,        // final boolean                   use_sel,
				clt_parameters.msDivideByArea,  // final boolean    msDivideByArea,  // final boolean                   divide_by_area,
				clt_parameters.msScaleProj,     //final double     msScaleProj,     // final double                    scale_projection,


				clt_parameters.stSmallDiff,     //       = 0.4;   // Consider merging initial planes if disparity difference below
				clt_parameters.stHighMix,       // stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
				world_hor,                      // final double []  world_hor, // horizontal plane normal (default [0.0, 1.0, 0.0])
				clt_parameters.show_histograms, // final boolean    show_histograms,
				clt_parameters.batch_run?-1:1, // -1,                       // debugLevel,                  // final int        debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
//		showDoubleFloatArrays sdfa_instance = null;
//		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?


		//==========================

		// Trying new class
		LinkPlanes lp = new LinkPlanes (clt_parameters, st);
		if (!batch_mode && clt_parameters.show_planes && (debugLevel > -1)){
			showPlaneData(
					"initial",
					clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					st, // SuperTiles st,
					null); // TilePlanes.PlaneData[][][] split_planes
		}

		// condition supertiles (create and manage links, merge)
		lp.conditionSuperTiles(
				st.planes,              // final TilePlanes.PlaneData [][] planes,
				10,                 // final int   max_num_merge_try,
				clt_parameters.batch_run?-2:0); // 1); // debugLevel);        // final int         debugLevel);
// Used only by conflicts (not processed currently)
		lp.calcStarValueStrength(
				true, // boolean set_start_planes,
				clt_parameters.plStarOrtho,    // orthoWeight,      // final double         orthoWeight,
				clt_parameters.plStarDiag,     // diagonalWeight,   // final double         diagonalWeight,
				clt_parameters.plStarPwr,      // starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
				clt_parameters.plStarWeightPwr,// starWeightPwr,    // final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
				clt_parameters.plWeightToDens, // weightToDens,     // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
				clt_parameters.plStarValPwr,   // starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
				2, // starSteps,        // final int            steps,
				st.planes,      // final TilePlanes.PlaneData [][] planes,
				clt_parameters.plPreferDisparity, // preferDisparity,  // final boolean        preferDisparity)
				debugLevel-2);

		// re-generate planes in the supertiles using previously calculated planes (for tghe tiles and their neighbors)
		// as hints, new planes will be assumed parallel to the known and possibly slightly offset in disparity
		if (clt_parameters.plDiscrEn) {
			st.regeneratePlanes(
					st.planes,                        // final TilePlanes.PlaneData [][] planes,
					clt_parameters.stMeasSel,         //      =     1   //Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					clt_parameters.plDispNorm,        //      =   2.0;  // Normalize disparities to the average if above
					clt_parameters.plMinPoints,       //      =     5;  // Minimal number of points for plane detection
					clt_parameters.plTargetEigen,     //      =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					clt_parameters.plFractOutliers,   //      =   0.3;  // Maximal fraction of outliers to remove
					clt_parameters.plMaxOutliers,     //      =    20;  // Maximal number of outliers to remove\
					clt_parameters.plPreferDisparity,
					geometryCorrection,
					clt_parameters.correct_distortions,

					clt_parameters.stSmplMode,        // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters
//					clt_parameters.stSmplSide,        // final int                        smplSide, //        = 2;      // Sample size (side of a square)
//					clt_parameters.stSmplNum,         // final int                        smplNum, //         = 3;      // Number after removing worst
//					clt_parameters.stSmplRms,         // final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,         // final boolean                 smplWnd,  // use window functions for the samples

//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

					clt_parameters.plDiscrTolerance,  // final double     plDiscrTolerance,     //     =   0.4;  // Maximal disparity difference from the plane to consider tile
					clt_parameters.plDiscrDispRange,  // final double     plDiscrDispRange,     //     =   0.6;  // Parallel move known planes around original know value for the best overall fit
					clt_parameters.plDiscrSteps,      // final int        plDiscrSteps,         //         =   3;    // Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
//					clt_parameters.plDiscrVariants,   // final int        plDiscrVariants,      //      =   100;  // Total number of variants to try (protect from too many planes)
					clt_parameters.plDiscrMode,       // final int        plDiscrMode,          //          =   3;    // What plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined

					clt_parameters.plDiscrVarFloor,   // final double     plDiscrVarFloor, //       =   0.03;  // Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
					clt_parameters.plDiscrSigma,      // final double     plDiscrSigma, //          =   0.05;  // Gaussian sigma to compare how measured data is attracted to planes
					clt_parameters.plDiscrBlur,       // final double     plDiscrBlur, //           =   0.1;   // Sigma to blur histograms while re-discriminating
					clt_parameters.plDiscrExclusivity,// final double     plDiscrExclusivity, //    =   0.5;   // Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
					clt_parameters.plDiscrExclus2,    // final double     plDiscrExclus2, //    =   0.5;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
					clt_parameters.plDiscrStrict,     // final boolean    plDiscrStrict, //         = true;   // When growing selection do not allow any offenders around (false - more these than others)
					clt_parameters.plDiscrCorrMax,    // final double     plDiscrCorrMax, //         = 0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
					clt_parameters.plDiscrCorrMerge,  // final double     plDiscrCorrMerge,  //     = 0.85;  // Attraction to different planes correlation that is high enough to merge planes
					clt_parameters.plDiscrSteal,      // final int        plDiscrSteal,      //         =   4;     // If offender has this number of tiles (including center) the cell can not be used
					clt_parameters.plDiscrGrown,      // final int        plDiscrGrown,      //         =   0;     // Only use tiles within this range from original selection

					clt_parameters.plDiscrXMedian,    // final double     plDiscrXMedian, //         = 1.5;   // Remove outliers from the final selection that have distance more than scaled median
					debugLevel, // -1,                         // debugLevel,                  // final int        debugLevel)
					clt_parameters.batch_run ? -1 : clt_parameters.tileX, // clt_parameters.tileX,
					clt_parameters.tileY);

			// condition the redefined planes
			lp.conditionSuperTiles(
					st.planes,          // final TilePlanes.PlaneData [][] planes,
					10,                 // final int   max_num_merge_try,
					debugLevel);        // final int         debugLevel);
			lp.calcStarValueStrength(
					true, // boolean set_start_planes,
					clt_parameters.plStarOrtho,    // orthoWeight,      // final double         orthoWeight,
					clt_parameters.plStarDiag,     // diagonalWeight,   // final double         diagonalWeight,
					clt_parameters.plStarPwr,      // starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
					clt_parameters.plStarWeightPwr,// starWeightPwr,    // final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
					clt_parameters.plWeightToDens, // weightToDens,     // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
					clt_parameters.plStarValPwr,   // starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
					2, // starSteps,        // final int            steps,
					st.planes,      // final TilePlanes.PlaneData [][] planes,
					clt_parameters.plPreferDisparity, // preferDisparity,  // final boolean        preferDisparity)
					debugLevel); //  - 2);

		}

		lp.setExclusiveLinks(
				st.planes, // final TilePlanes.PlaneData [][] planes,
				lp.getExNeibCost(), // final double                    max_cost,
				debugLevel,                  // final int        debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
/*
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
*/
		double [][][]  dispStrength = st.getDisparityStrengths(
				clt_parameters.stMeasSel); // int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)
		boolean [][] tileSel =  st.getMeasurementSelections(
				clt_parameters.stMeasSel); // int        stMeasSel) //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert)


		TilePlanes.PlaneData[][][]       split_planes =   // use original (measured planes. See if smoothed are needed here)
				st.breakPlanesToPairs(
						st.getPlanes(), // Mod(),             // final TilePlanes.PlaneData[][] center_planes, // measured_planes,
						st.getPlanes(), // Mod(),             // final TilePlanes.PlaneData[][] neib_planes,   //mod_planes,
						clt_parameters.plSplitPull ,          // final double                   center_pull,
						clt_parameters.plSplitMinNeib ,       // min_neibs, // 2
						clt_parameters.plSplitMinWeight,      // final double splitMinWeight,  //     =  2.0;  // Minimal weight of split plains to show
						clt_parameters.plSplitMinQuality,     // final double splitMinQuality, //    =  1.1;  // Minimal split quality to show

						clt_parameters.plPreferDisparity,
						debugLevel,
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
//					clt_parameters.stFloor,          // final double                     strength_floor,
//					clt_parameters.stPow,            // final double                     strength_pow,
					//					clt_parameters.dbg_migrate && clt_parameters.stSmplMode , // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					clt_parameters.stSmplMode , // final boolean                    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters
//					clt_parameters.stSmplSide , // final int                        smplSide, //        = 2;      // Sample size (side of a square)
//					clt_parameters.stSmplNum , // final int                        smplNum, //         = 3;      // Number after removing worst
//					clt_parameters.stSmplRms , // final double                     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples

//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

					debugLevel, // 1,                               // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			if (debugLevel > -1){
				System.out.println("replaceBrokenPlanes(): Replaced " + numSplitPlanes + " planes by pairs.");
			}
			// now re-create connections and best neighbors
			lp = new LinkPlanes (clt_parameters, st);
			lp.matchPlanes(
					st.planes, // final TilePlanes.PlaneData [][] planes,
					debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);

			lp.interPlaneCosts( // not used yet, just for testing
					true, // final boolean                   en_sticks, // allow merging with bad plates
					st.planes, // final TilePlanes.PlaneData [][] planes,
					debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);

			lp.filterNeighborPlanes(
					st.planes, // final TilePlanes.PlaneData [][] planes,
					true, // final boolean merge_low_eigen,
					debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);

// Try replacing lp.selectNeighborPlanesMutual with
//			lp.setExclusiveLinks()
			double [][] quality_stats2 = lp.selectNeighborPlanesMutual(
//					false, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
					true, // final boolean                   en_sticks, // allow merging with bad plates
					st.planes,              // final TilePlanes.PlaneData [][] planes,
					debugLevel, // 0, // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
			if (debugLevel>100) System.out.println(quality_stats2.length);

			System.out.println("Testing - overwriting selectNeighborPlanesMutual() results with setExclusiveLinks()");

// Just overwrite results of the previous method
			lp.setExclusiveLinks(
					st.planes, // final TilePlanes.PlaneData [][] planes,
					lp.getExNeibCost(), // final double                    max_cost,
					debugLevel,                  // final int        debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);

			st.resolveConflicts(
					lp, // LinkPlanes lp,
					clt_parameters.plMaxEigen,
					clt_parameters.plConflDualTri,  // boolean    conflDualTri,  // Resolve dual triangles conflict (odoodo)
					clt_parameters.plConflMulti,    // boolean    conflMulti,    // Resolve multiple odo triangles conflicts
					clt_parameters.plConflDiag,	    // boolean    conflDiag,     // Resolve diagonal (ood) conflicts
					clt_parameters.plConflStar,     // boolean    conflStar,     // Resolve all conflicts around a supertile
					clt_parameters.plStarSteps,     // int starSteps, // How far to look around when calculationg connection cost
					clt_parameters.plStarOrtho,     // double     orthoWeight,
					clt_parameters.plStarDiag,      // double     diagonalWeight,
					clt_parameters.plStarPwr,       // double     starPwr, // Divide cost by number of connections to this power
					clt_parameters.plStarWeightPwr, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
					clt_parameters.plWeightToDens,  // double     weightToDens,    // // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
					clt_parameters.plStarValPwr,    // double     starValPwr, //  Raise value of each tile before averaging
					clt_parameters.plDblTriLoss,    // double     diagonalWeight,
					clt_parameters.plNewConfl,      // boolean    preferDisparity, // Allow more conflicts if overall cost is reduced
					clt_parameters.plMaxChanges,    // int        maxChanges,  // Maximal number of simultaneous connection changes around one tile (0 - any)
					clt_parameters.plPreferDisparity,
					debugLevel, // 1, // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);
		} // if (clt_parameters.plSplitApply)

/*
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

*/

		int max_num_tries = 20;
		if (clt_parameters.plIterations > 0) {
			st.planes_mod = lp.planesSmoothAndMerge(
					st.planes, // final TilePlanes.PlaneData[][] planes, // planes will be modified
					max_num_tries, // final int                      max_num_tries,
					debugLevel); // final int                      debugLevel)

			// currently results of below calcStarValueStrength() are not used, just to fill instances fields
			lp.calcStarValueStrength(
					true, // boolean set_start_planes,
					clt_parameters.plStarOrtho,    // orthoWeight,      // final double         orthoWeight,
					clt_parameters.plStarDiag,     // diagonalWeight,   // final double         diagonalWeight,
					clt_parameters.plStarPwr,      // starPwr,          // final double         starPwr,    // Divide cost by number of connections to this power
					clt_parameters.plStarWeightPwr,// starWeightPwr,    // final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
					clt_parameters.plWeightToDens, // weightToDens,     // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
					clt_parameters.plStarValPwr,   // starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
					2, // starSteps,        // final int            steps,
					st.planes,      // final TilePlanes.PlaneData [][] planes,
					clt_parameters.plPreferDisparity, // preferDisparity,  // final boolean        preferDisparity)
					0); // debugLevel);
		} else { //if (clt_parameters.plIterations > 0)
			st.planes_mod = st.planes; // just use the measured ones
		}

		// Add "missing" links with relaxed threshold
		if (lp.getRelaxComplete() > 0.0) {
			int total = 0;
			int num_added = lp.addFinalLinks (
					st.planes_mod, // final TilePlanes.PlaneData[][] planes,
					lp.getRelaxComplete() * lp.getExNeibCost(), // final double                   threshold_sq,
					lp.getRelaxComplete() * lp.getExNeibCost(), // final double                   threshold_corn,
					lp.getRelaxComplete() * lp.getExNeibCost(), //final double                   threshold_hyp,
					clt_parameters.plFillSquares, // final boolean                  en_sq,
					clt_parameters.plCutCorners, // final boolean                  en_corn,
					clt_parameters.plHypotenuse, // final boolean                  en_hyp,
					1); // final int                      debugLevel)
			total += num_added;
			if (lp.getRelaxComplete2() > 0.0) {
				num_added = lp.addFinalLinks (
						st.planes_mod, // final TilePlanes.PlaneData[][] planes,
						lp.getRelaxComplete2() * lp.getExNeibCost(), // final double                   threshold_sq,
						lp.getRelaxComplete2() * lp.getExNeibCost(), // final double                   threshold_corn,
						lp.getRelaxComplete2() * lp.getExNeibCost(), //final double                   threshold_hyp,
						clt_parameters.plFillSquares, // final boolean                  en_sq,
						clt_parameters.plCutCorners, // final boolean                  en_corn,
						clt_parameters.plHypotenuse, // final boolean                  en_hyp,
						1); // final int                      debugLevel)
				total += num_added;
				System.out.println("Total: added "+total+ " links between supertiles");
			}

/*
			while (true) {
				int num_added = 0;
				if (clt_parameters.plFillSquares){
					num_added += lp.fillSquares(
							st.planes_mod,
							lp.getRelaxComplete() * lp.getExNeibCost(),
							2, // debugLevel, // 0, // final int debugLevel)
							clt_parameters.tileX,
							clt_parameters.tileY);

				}
				if (debugLevel > -1) {
					System.out.println("after fillSquares() added "+num_added);
				}
				if (clt_parameters.plCutCorners){
					num_added += lp.cutCorners(
							st.planes_mod,
							lp.getRelaxComplete() * lp.getExNeibCost(),
							2, // debugLevel, // 0, // final int debugLevel)
							clt_parameters.tileX,
							clt_parameters.tileY);
				}
				if (debugLevel > -1) {
					System.out.println("after cutCorners() added (cumulative) "+num_added);
				}
				if (num_added == 0) break;
			}
			*/
		}


		// filter out weak planes, create boolean array [per-supertile][per disparity plane]
		boolean [][] selected_planes =	st.selectPlanes(
				clt_parameters.plDispNorm,
				clt_parameters.plMaxEigen, //  maxEigen,
				clt_parameters.plMinStrength, //  minWeight,
				true, // final boolean weak_connected, // select connected planes even if they are weak/thick
				st.getPlanesMod());

		if (!batch_mode && clt_parameters.show_planes && (debugLevel > -2)){
			showPlaneData(
					"final",
					clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					st, // SuperTiles st,
					split_planes); // TilePlanes.PlaneData[][][] split_planes
		}

		st.tileSurface.createTileShells (
				clt_parameters.msUseSel,            // final boolean                   use_sel,
				clt_parameters.msDivideByArea,      // final boolean                   divide_by_area,
				clt_parameters.msScaleProj,         // final double                    scale_projection,
				clt_parameters.msFractUni,          // final double                    fraction_uni,
				st.planes_mod, // st.planes,        // final TilePlanes.PlaneData [][] planes,
				debugLevel,                         // final int        debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);

		st.tileSurface.InitTilesAssignment(
				true,
				dispStrength, // final double [][][]                            dispStrength,
				tileSel,      // final boolean [][]                             tileSel,
				debugLevel);    // final int                                      debugLevel,


		if (debugLevel > -10){
			return; // just cut off the rest
		}

		createShells_old (
				clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				st, // SuperTiles st,
				selected_planes, // boolean [][] selected_planes,
				geometryCorrection, // GeometryCorrection geometryCorrection,
				debugLevel); // final int         debugLevel)
	}


	public void showPlaneData(
			String suffix,
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			SuperTiles st,
			TilePlanes.PlaneData[][][] split_planes
			)
	{
			double [] split_lines = (split_planes==null) ? null: st.showSplitLines(
					split_planes,
					1,                              // final int debugLevel)
					clt_parameters.tileX,
					clt_parameters.tileY);

			int [] wh = st.getShowPlanesWidthHeight();
			double [][] plane_data_nonan_meas = st.getShowPlanes(
					st.getPlanes(),
					clt_parameters.plMinStrength, //  minWeight,
					clt_parameters.plMaxEigen, //  maxEigen,
					clt_parameters.plDispNorm,
					false, //boolean use_NaN)
					0.0,
					10.0);
			double [][] plane_data_nan_meas = st.getShowPlanes(
					st.getPlanes(),
					clt_parameters.plMinStrength, //  minWeight,
					clt_parameters.plMaxEigen, //  maxEigen,
					clt_parameters.plDispNorm,
					true, //boolean use_NaN)
					0.0,
					10.0);
			double [][] plane_data_nonan = st.getShowPlanes(
					(st.planes_mod != null) ? st.getPlanesMod():st.getPlanes(),
							clt_parameters.plMinStrength, //  minWeight,
							clt_parameters.plMaxEigen, //  maxEigen,
							clt_parameters.plDispNorm,
							false, //boolean use_NaN)
							0.0,
							10.0);
			double [][] plane_data_nan = st.getShowPlanes(
					(st.planes_mod != null) ? st.getPlanesMod():st.getPlanes(),
							clt_parameters.plMinStrength, //  minWeight,
							clt_parameters.plMaxEigen, //  maxEigen,
							clt_parameters.plDispNorm,
							true, //boolean use_NaN)
							0.0,
							10.0);
			double [][] plane_data = new double [plane_data_nonan_meas.length + plane_data_nan_meas.length+ plane_data_nonan.length + plane_data_nan.length + 3][];
			int indx = 0;
			for (int i = 0; i < plane_data_nonan_meas.length; i++){
				plane_data[indx++] = plane_data_nonan_meas[i];
			}
			for (int i = 0; i < plane_data_nan_meas.length; i++){
				plane_data[indx++] = plane_data_nan_meas[i];
			}
			for (int i = 0; i < plane_data_nonan.length; i++){
				plane_data[indx++] = plane_data_nonan[i];
			}
			for (int i = 0; i < plane_data_nan.length; i++){
				plane_data[indx++] = plane_data_nan[i];
			}
			plane_data[indx++] = split_lines;
			if (indx <2) { // out of bound
				System.out.println("BUG: insufficient data");
				return;
			}
			plane_data[indx] = plane_data[indx-2].clone(); // java.lang.ArrayIndexOutOfBoundsException: -1
			for (int i = 0; i < plane_data[indx].length;i++){
				if (Double.isNaN(plane_data[indx][i])) plane_data[indx][i] = 0.0;
				if ((plane_data[indx-1] != null) && (plane_data[indx] != null)) {
					if (plane_data[indx-1][i] > 0) plane_data[indx][i] = Double.NaN;
				}
			}
			indx++;
			plane_data[indx] = new double [wh[0]*wh[1]];
			for (int i = 0; i < plane_data[indx].length; i++){
				int dbg_stx = (i % wh[0]) /superTileSize;
				int dbg_sty = (i / wh[0]) /superTileSize;
				int dbg_nsTile = dbg_stx + dbg_sty * (wh[0]/superTileSize);
				if (st.planes_mod != null) {
					if (st.planes_mod[dbg_nsTile] == null){
						plane_data[indx][i] = Double.NaN;
					}else {
						int dbg_nl = 0;
						for (int j = 0; j < st.planes_mod[dbg_nsTile].length; j++){
							if (st.planes_mod[dbg_nsTile][j] != null) dbg_nl++;
						}
						plane_data[indx][i] = dbg_nl;
					}
				} else if (st.planes != null) {
					if (st.planes[dbg_nsTile] == null){
						plane_data[indx][i] = Double.NaN;
					}else {
						int dbg_nl = 0;
						for (int j = 0; j < st.planes[dbg_nsTile].length; j++){
							if (st.planes[dbg_nsTile][j] != null) dbg_nl++;
						}
						plane_data[indx][i] = dbg_nl;
					}

				}
			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "plane_data-"+suffix);
	}


	public void createShells_old (
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			SuperTiles st,
			boolean [][] selected_planes,
			GeometryCorrection geometryCorrection,
			final int         debugLevel)
	{

		//*******************************************************************************
		// detect connected "shells" by running wave algorithm over multi-plane supertiles, create integer array (same structure as planes and selected_planes
		// Each element [per supertile][per disparity plane] is either 0 (not used) or unique shell index (started from 1)
		int [][] shells = st.createSupertileShells(
				selected_planes, // boolean[][]               selection, // may be null
				false,             // boolean                   use_all, // use plane 0 even if there are more than 1
				clt_parameters.plKeepOrphans,             // true, // boolean                   keep_orphans, // single-cell shells
				clt_parameters.plMinOrphan, // orphan_strength // add separate parameter
				st.getPlanesMod(), // TilePlanes.PlaneData [][] planes,
				1, // final int debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);
		// save shell indices with SuperTiles instance
		st.setShellMap(shells); // persistent
		// Create array [per shell][per tile] of shell disparities. tiles of unused supertiles are Double.NaN
		double [][] surfaces = st.getShowShells(
				0, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
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



		if (clt_parameters.show_planes){
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			int [] wh = st.getShowPlanesWidthHeight();
			double [][] plane_data = st.getShowShells(
					0, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
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
					1, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					false,               // boolean show_connections,
					false,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)

			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells2_all");

			plane_data = st.getShowShells(
					0, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
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
					1, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					true,                // boolean show_connections,
					false,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells2_all_arrows");


			plane_data = st.getShowShells(
					0, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					false,               // boolean show_connections,
					true,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells");

			plane_data = st.getShowShells(
					1, // int  nlayer, // over multi-layer - do not render more than nlayer on top of each other
					st.getPlanesMod(),   // TilePlanes.PlaneData [][] planes,
					st.getShellMap(), // shells,              // int [][] shells,
					100,                 // int max_shells,
					clt_parameters.plFuse,// boolean fuse,
					false,               // boolean show_connections,
					true,               // boolean use_NaN,
					10.0,                 // double arrow_dark,
					10.0);               // double arrow_white)
			sdfa_instance.showArrays(plane_data, wh[0], wh[1], true, "shells2");
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

			CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
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
//			double            this_sure,        // minimal strength to be considered definitely background
			double            ex_strength,   // minimal 4-corr strength to trust tile
			double            ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles
			double            this_maybe,       // maximal strength to ignore as non-background
			double            sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
			double            super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds

			int               disparity_index,  // index of disparity value in disparity_map == 2 (0,2 or 4)
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch
			final boolean     updateStatus,
			final int         debugLevel)
	{
		final boolean                       batch_mode = clt_parameters.batch_run; //disable any debug images
		int debugLevelInner = batch_mode ? -5: debugLevel;
		CLTPass3d scan_bg =   clt_3d_passes.get(bg_scan_index); //
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
		CLTPass3d scan_lm =   getLastMeasured(-1);
//		boolean [] border_tiles = (scan_lm != null) ? scan_lm.getBorderTiles() : null;

		showDoubleFloatArrays sdfa_instance = null;
		if (!batch_mode && (debugLevel > -1)) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		//TODO: for next passes - combine all selected for previous passes (all passes with smaller disparity)
		boolean show_st =   !batch_mode && ( clt_parameters.stShow || (debugLevel > 1));

		boolean [] these_tiles;
		int [] replaced =  null; // +1 - hor, +2 - vert
//		int [] replaced0 = null; // +1 - hor, +2 - vert

		if (!clt_parameters.ortho_old) {
			replaced = combineOrthoDisparity(
					scan_prev, // final CLTPass3d   scan,        // scan data
					clt_parameters.or_hor,       // true;  // Apply ortho correction to horizontal correlation (vertical features)
					clt_parameters.or_vert,      // true;  // Apply ortho correction to vertical correlation (horizontal features)
					clt_parameters.or_sigma,     // 2.0;   // Blur sigma: verically for horizontal correlation, horizontally - for vertically
					clt_parameters.or_sharp,     // 0.5;   // 3-point sharpening (-k, +2k+1, -k)
					clt_parameters.or_scale,     // 2.0;   // Scale ortho correlation strength relative to 4-directional one
					clt_parameters.or_offset,    // 0.1;   // Subtract from scaled correlation strength, limit by 0
					clt_parameters.or_asym ,     // 1.5;   // Minimal ratio of orthogonal strengths required for dis[parity replacement
					clt_parameters.or_threshold, // 1.5;   // Minimal scaled offsetg ortho strength to normal strength needed for replacement
					clt_parameters.or_absHor,    // 0.15;  // Minimal horizontal absolute scaled offset ortho strength needed for replacement
					clt_parameters.or_absVert,   // 0.19;  // Minimal vertical absolute scaled offset ortho strength needed for replacement
					clt_parameters.or_maxDisp,   // 5.0;   // Maximal disparity to apply ortho correction
					clt_parameters.show_ortho_combine, //  show replacement of disparity/strength by those of hor/vert

					debugLevelInner);

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
				if (debugLevel > 0){
					System.out.println("fixVerticalPoles() replaced "+ numFixed+ " tiles.");
				}
//				replaced0 = replaced.clone();
				for (int i = 0; i < replaced.length; i++){
					if (tilesHor[i]) replaced[i] |= 1;
				}
			}
			these_tiles = FilterScan(
					scan_prev,        // final CLTPass3d   scan,
					scan_bg.selected, // get from selected in clt_3d_passes.get(0);
					clt_parameters.ex_min_over,// when expanding over previously detected (by error) background, disregard far tiles
//					border_tiles,      // final boolean  [] border_tiles,      // last measured boirder tiles
					scan_lm,          // final CLTPass3d   last_meas,         // last measured scan (with border_tiles and getSecondMaxDiff
					replaced,         // final int      [] horVertMod, // +1 - modified by hor correlation, +2 - modified by vert correlation (or null)
					disparity_far,    // final double      disparity_far,    //
					disparity_near,   // final double      disparity_near,   //
//					this_sure,        // final double      this_sure,        // minimal strength to be considered definitely background
					ex_strength,      // final double      ex_strength,   // minimal 4-corr strength to trust tile
					ex_nstrength,     // final double      ex_nstrength, // minimal 4-corr strength divided by channel diff for new (border) tiles

					this_maybe,       // final double      this_maybe,       // maximal strength to ignore as non-background
					sure_smth,        // final double      sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
					super_trust,      // final double      super_trust,      // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
					// using plates disparity/strength - averaged for small square sets of tiles. If null - just use raw tiles
					null,        // final double [][] plate_ds,  // disparity/strength last time measured for the multi-tile squares. Strength =-1 - not measured. May be null
					true,  // final boolean     keep_raw_fg,  // do not replace raw tiles by the plates, if raw is closer (like poles)
					0.0, // 					final double      scale_filtered_strength_pre, // scale plate_ds[1] before comparing to raw strength
					0.0,// final double      scale_filtered_strength_post,// scale plate_ds[1] when replacing raw (generally plate_ds is more reliable if it exists)
					batch_mode, // true, // clt_parameters.show_filter_scan,
					clt_parameters.min_clstr_seed, // clt_parameters.min_clstr_seed
					clt_parameters.min_clstr_weight,  // double     min_weight // minimal total weight of the cluster
					clt_parameters.min_clstr_max,    // double     min_max_weight // minimal value of the maximal strengh in the cluster
					clt_parameters.min_clstr_lone,  // int        min_area,  // minimal number of pixels
					clt_parameters.fill_gaps, // fill_gaps, // int depth, // same as grow - odd - 4 directions, even - 8
					clt_parameters.fill_final,             // final boolean     poison_gaps,      // Do not fill gaps that have even single "poisoned" tile
					true,             // final boolean     zero_gap_strength, // set strength to zero when covering gaps
					//					final int         threadsMax,  // maximal number of threads to launch
					//					final boolean     updateStatus,
					batch_mode ? -5: 2); //debugLevel);
		} else { //!ortho_old
			these_tiles= combineHorVertDisparity_old(
				scan_prev,                     // final CLTPass3d   scan,
				scan_bg.selected, // clt_3d_passes.get(0).selected, // final boolean [] bg_tiles,          // get from selected in clt_3d_passes.get(0);
				disparity_far,    //
				disparity_near,   //
				ex_strength,      // this_sure,        // minimal strength to be considered definitely background
				this_maybe,       // maximal strength to ignore as non-background
				sure_smth,        // if 2-nd worst image difference (noise-normalized) exceeds this - do not propagate bgnd
				clt_parameters,
				debugLevelInner);

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

		// try renovated supertiles. Do twice to show both original and blurred histograms
		double [] dbg_orig_disparity =  null;
		double [] dbg_with_super_disp = null;
		double [] dbg_outlayers =       null;
//		boolean [] grown = these_tiles.clone();

		if (use_supertiles || show_st) {// currently not use (only to show)
			String [] dbg_st_titles = {"raw", "blurred"+clt_parameters.stSigma,"max-min-max"};
			double [][] dbg_hist = new double[dbg_st_titles.length][];

			scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					0.0, // NO BLUR double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters

//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // Use window functions for the samples

//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[0] = scan_prev.getSuperTiles().showDisparityHistogram();

			SuperTiles st = scan_prev.setSuperTiles(
					clt_parameters.stStepNear,       // double     step_disparity,
					clt_parameters.stStepFar,        // double     step_near,
					clt_parameters.stStepThreshold,  // double     step_threshold,
					clt_parameters.stMinDisparity,   // double     min_disparity,
					clt_parameters.grow_disp_max,   // double     max_disparity,
//					clt_parameters.stFloor,          // double     strength_floor,
//					clt_parameters.stPow,            // double     strength_pow,
					clt_parameters.stSigma, // with blur double     stBlurSigma)
					false, //clt_parameters.stSmplMode,  // Use sample mode (false - regular tile mode)
					clt_parameters.mlfp,         // Filter parameters
//					clt_parameters.stSmplSide,  // Sample size (side of a square)
//					clt_parameters.stSmplNum,   // Number after removing worst
//					clt_parameters.stSmplRms,   // Maximal RMS of the remaining tiles in a sample
//					clt_parameters.stSmplWnd,   // boolean                 smplWnd,  // use window functions for the samples

//					clt_parameters.fs_max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//					clt_parameters.fs_max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//					clt_parameters.fs_damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//					clt_parameters.fs_min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//					clt_parameters.fs_transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//					clt_parameters.fs_far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//					clt_parameters.fs_far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

					clt_parameters.stMeasSel); // bitmask of the selected measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			dbg_hist[1] = scan_prev.getSuperTiles().showDisparityHistogram();

			dbg_hist[2] = scan_prev.getSuperTiles().showMaxMinMax();

			int hist_width0 =  scan_prev.getSuperTiles().showDisparityHistogramWidth();
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
	public void thirdPassSetup_old( // prepare tile tasks for the second pass based on the previous one(s)
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double            disparity_far,    //
			double            disparity_near,   //
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		boolean [] these_tiles = scan_prev.getSelected();
		double [] this_disparity     = scan_prev.getDisparity().clone(); // returns a copy of the FPGA-generated disparity combined with the target one
		double [] this_strength =      scan_prev.getStrength().clone(); // cloned, can be modified/ read back

		scan_prev.fixNaNDisparity(
				null, // border, // boolean [] select,   // which tiles to correct (null - all)
				this_disparity, // scan_prev.getDisparity(), // double [] disparity,
				this_strength); // scan_prev.getStrength()); // double [] strength)


		SuperTiles st = scan_prev.getSuperTiles();

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



	public void thirdPassSetupSurf( // prepare tile tasks for the second pass based on the previous one(s)
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			double            disparity_far,    //
			double            disparity_near,   //
			GeometryCorrection geometryCorrection,
			final int         threadsMax,  // maximal number of threads to launch
			final boolean     updateStatus,
			final int         debugLevel)
	{
		CLTPass3d scan_prev = clt_3d_passes.get(clt_3d_passes.size() -1); // get last one
		SuperTiles st = scan_prev.getSuperTiles();
		TileSurface tileSurface = st.getTileSurface();
		if (tileSurface == null){
			System.out.println("tileSurface does not exist, aborting");
			return;
		}

		int num_surf_proc = clt_parameters.tsNumClust;
		int num_surf = tileSurface.getSurfaceDataLength();
		if (num_surf_proc > num_surf) num_surf_proc = num_surf;
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);
		int numClusters = 0; // needed?
		for (int ns = 0 ; ns < num_surf_proc; ns ++){
			if (ns == 73){
				System.out.println("thirdPassSetupSurf() ns = "+ns);
			}
			// return to original dimensions
			double [][] disparityTask = new double [tilesY][tilesX];
			int [][]    tile_op =       new int [tilesY][tilesX];
			boolean [] borderTiles =    new boolean[tilesY*tilesX]; // to zero alpha in the images
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				int indx =  tilesX * ty + tx;
				disparityTask[ty][tx] = 0.0;
				tile_op[ty][tx] =       0;
				borderTiles[indx] =     false;
			}
			SurfaceData sd = tileSurface.getSurfaceData(ns);
			boolean [] sel = sd.getSelected();
			for (int i = 0; i < sel.length; i++) if (sel[i]){
				int indx = sd.getImageIndex(i);
				int tx = indx % tilesX;
				int ty = indx / tilesX;
				tile_op[ty][tx] = op;
				disparityTask[ty][tx] = sd.getDisparity(i);
				borderTiles[indx] = sd.getBorder(i);
			}
			// Create FPGA task for this cluster

			CLTPass3d scan_next =new CLTPass3d(this);
			scan_next.disparity =     disparityTask;
			scan_next.tile_op =       tile_op;
			scan_next.border_tiles =  borderTiles;
			clt_3d_passes.add(scan_next);
			numClusters++;
		}
		if (debugLevel > -1){
			System.out.println("thirdPassSetupSurf():  created "+ numClusters+ " FPGA passes.");
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
//					double disp = disparity[(bounds.y + y) * tilesX + (bounds.x + x)];
					double disp = (disparity == null)? min_disparity:( disparity[(bounds.y + y) * tilesX + (bounds.x + x)]);
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
			double [] disparity, // null - use min_disparity
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
					double disp = (disparity == null)? min_disparity:( disparity[(bounds.y + y) * tilesX + (bounds.x + x)]);
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
			double d = selected[i]? ((disparity.length >1) ? disparity[i] : disparity[0]):Double.NaN;
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
