/**
 **
 ** CLTPass3d - A single processing "pass" over the image set. May be both actual
 ** FPGA operation or result of merging data from multiple passes
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  CLTPass3d.java is free software: you can redistribute it and/or modify
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
import java.util.concurrent.atomic.AtomicInteger;

public class CLTPass3d{
		public  double [][]     disparity; // per-tile disparity set for the pass[tileY][tileX] 
		public  int    [][]     tile_op;   // what was done in the current pass
		public  double [][]     disparity_map =  null; // add 4 layers - worst difference for the port
		double []               calc_disparity = null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		                                       // using horizontal features and corr_magic_scale
		// used directly in TileProcessor.compositeScan()
		double []               calc_disparity_hor =   null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		double []               calc_disparity_vert =  null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		double []               calc_disparity_combo = null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		double []               strength =             null; // composite strength, initially uses a copy of raw 4-sensor correleation strength
		double []               strength_hor =         null; // updated hor strength, initially uses a copy of raw measured
		double []               strength_vert =        null; // updated hor strength, initially uses a copy of raw measured
		// Bg disparity & strength is calculated from the supertiles and used instead of the tile disparity if it is too weak. Assuming, that
		// foreground features should have good correlation details, and if the tile does not nhave them it likely belongs to the background.
		// calculate disparity and strength from the (lapped) supertiles, using lowest allowed (>= minBgDisparity) disparity histogram maximums
		// of the supertiles this tile belongs to  
		private double          minBgDisparity =        0.0;
		private double          minBgFract =            0.0; // Use the lowest maximum if the strength strength (of all maximus >= minBgDisparity)
		                                                     // exceeds minBgFract, otherwise proceed to the next one (and accumulate strength)
		private double []       bgTileDisparity =      null;
		private double []       bgTileStrength =       null;
		public  boolean []      border_tiles =         null;      // these are border tiles, zero out alpha
		public  boolean []      selected =             null;          // which tiles are selected for this layer
		public  double [][][][] texture_tiles;
		public double [][]      max_tried_disparity =  null; //[ty][tx] used for combined passes, shows maximal disparity wor this tile, regardless of results
		public  boolean         is_combo =             false;
		public  boolean         is_measured =          false;
		public  String          texture = null; // relative (to x3d) path
		public  Rectangle       bounds;
		public  int             dbg_index;
		public  int             disparity_index = ImageDtt.DISPARITY_INDEX_CM; // may also be ImageDtt.DISPARITY_INDEX_POLY
		
		SuperTiles              superTiles = null;
		TileProcessor           tileProcessor;
		public CLTPass3d (TileProcessor tileProcessor)
		{
			this.tileProcessor = tileProcessor;
		}
		public TileProcessor getTileProcessor()
		{
			return this.tileProcessor;
		}
		
		public  void            updateSelection(){ // add updating border tiles?
			int tilesX = tileProcessor.getTilesX();
			int tilesY = tileProcessor.getTilesY();
			selected = new boolean[tilesY*tilesX];
			int minX = tilesX, minY = tilesY, maxX = -1, maxY = -1;
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx < tilesX; tx++){
				if (texture_tiles[ty][tx] != null) {
					selected[ty * tilesX + tx] = true;
					if (maxX < tx) maxX  = tx;
					if (minX > tx) minX  = tx;
					if (maxY < ty) maxY  = ty;
					if (minY > ty) minY  = ty;
				} else {
					selected[ty * tilesX + tx] = false; // may be omitted 
				}
			}
			bounds = new Rectangle(minX, minY, maxX - minX +1, maxY - minY +1 );
		}
		
		public boolean isProcessed(){
			return calc_disparity != null;	
		}
		
		public boolean isMeasured(){
			return is_measured;
//			return (disparity_map != null) && (disparity != null); // 	disparity == null for composite scans
		}
		
		public boolean isCombo(){
			return is_combo;	
		}
		
		/**
		 * Called after each measurement
		 */
		public void resetProcessed(){ 
			
			fixNaNDisparity();
			
			calc_disparity =       null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			calc_disparity_hor =   null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			calc_disparity_vert =  null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			calc_disparity_combo = null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			strength =             null; // composite strength, initially uses a copy of raw 4-sensor correleation strength
			strength_hor =         null; // updated hor strength, initially uses a copy of raw measured
			strength_vert =        null; // updated hor strength, initially uses a copy of raw measured
			bgTileDisparity =      null;
			bgTileStrength =       null;
//			border_tiles =         null;      // these are border tiles, zero out alpha
//			selected =             null;          // which tiles are selected for this layer
			superTiles =           null;
			
			
		}
		/**
		 * Get FPGA-calculated per-tile maximal differences between the particular image and the average one.
		 * @return per-camera sesnor array of line-scan differences
		 */

		public double [][] getDiffs (){
			double  [][] these_diffs =    new double[ImageDtt.QUAD][];
			for (int i = 0; i< ImageDtt.QUAD; i++) these_diffs[i] = disparity_map[ImageDtt.IMG_DIFF0_INDEX + i];
			return these_diffs;
		}
		
		public void resetCalc(){ // only needed if the same task was reused
			calc_disparity = null;
			strength =       null;
			strength_hor =   null;
			strength_vert =  null;
			superTiles =     null;
			
		}
		public boolean [] getSelected(){
			return selected;
		}

		public void fixNaNDisparity()
		{
			fixNaNDisparity(
					null,
					disparity_map[disparity_index],
					disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX]);
			fixNaNDisparity(
					null,
					disparity_map[ImageDtt.DISPARITY_INDEX_HOR],
					disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH]);
			fixNaNDisparity(
					null,
					disparity_map[ImageDtt.DISPARITY_INDEX_VERT],
					disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH]);
		}
		
		
		public void fixNaNDisparity(
				boolean [] select,   // which tiles to correct (null - all)
				double [] disparity,
				double [] strength)
		{
			// depends on direction, but that is OK - just converge faster when smoothing
			int tilesX = tileProcessor.getTilesX();
			int tilesY = tileProcessor.getTilesY();
			int [] dirs8 = {-tilesX,  -tilesX + 1, 1, tilesX +1, tilesX, tilesX - 1, -1, -tilesX - 1};
			for (int ty = 1; ty < (tilesY -1); ty ++) for (int tx = 1; tx < (tilesX -1); tx++){
				int nt = ty * tilesX + tx;
				if (Double.isNaN(disparity[nt]) && ((select == null) || select[nt])) {
					if (strength != null) strength[nt] = 0.0;
					double sd = 0.0, sw = 0.0;
					for (int dir=0; dir < dirs8.length; dir++){
						int nt1 = nt + dirs8[dir];
//						if (!Double.isNaN(disparity[nt1]) && ((select == null) || !select[nt1])) {
						if (!Double.isNaN(disparity[nt1])) { // for wide borders - use neighbors already defined too
							double w = (strength == null) ? 1.0 : strength[nt1];
							sd += w * disparity[nt1];
							sw += w;
						}
					}
					if (sw > 0.0) sd /= sw;
					disparity[nt] = sd;
				}
			}
			// on top/bottom/right/left rows replace NaN disparity with 0.0;
			for (int ty = 0; ty < tilesY; ty ++) {
				int nt = ty * tilesX + 0;
				if (Double.isNaN(disparity[nt]) && ((select == null) || select[nt])) {
					if (strength != null) strength[nt] = 0.0;
					disparity[nt] = 0.0;
				}
				nt = ty * tilesX + tilesX -1;
				if (Double.isNaN(disparity[nt]) && ((select == null) || select[nt])) {
					if (strength != null) strength[nt] = 0.0;
					disparity[nt] = 0.0;
				}
			}
			for (int tx = 0; tx < tilesX; tx ++) {
				int nt = 0 * tilesX + tx;
				if (Double.isNaN(disparity[nt]) && ((select == null) || select[nt])) {
					if (strength != null) strength[nt] = 0.0;
					disparity[nt] = 0.0;
				}
				nt = (tilesY -1) * tilesX + tx;
				if (Double.isNaN(disparity[nt]) && ((select == null) || select[nt])) {
					if (strength != null) strength[nt] = 0.0;
					disparity[nt] = 0.0;
				}
			}
		}
		
		public double [] combineHorVertStrength(
				boolean combineHor,
				boolean combineVert)
		{
			getStrength();     // clone if not done yet
			if (combineHor){
				double [] hstrength = getHorStrength();
				for (int i = 0; i < strength.length; i++) {
					if (strength[i] < hstrength[i]) strength[i] = hstrength[i];
				}
			}
			if (combineVert){
				double [] vstrength = getVertStrength();
				for (int i = 0; i < strength.length; i++) {
					if (strength[i] < vstrength[i]) strength[i] = vstrength[i];
				}
			}
			return strength;
		}
		
		public double [] combineSuper(
				boolean updateStrength, // use ST strength if true, keep original (update disparity only) if false
				double  stStrengthScale,
				double useSuper){
			if (bgTileDisparity == null) { // no supertile disparity is available
				return null;
			}

			double [] strength = getStrength();
			double [] disparity = getDisparity(0);
			
			for (int i = 0; i < disparity.length; i++){
				if (strength[i] < useSuper)  {
					disparity[i] = bgTileDisparity[i];
					if (updateStrength) strength[i] =  stStrengthScale*bgTileStrength[i];
				}
			}
			return disparity;
		}
		
		/**
		 * Returns per-tile correlation "strength". Initially - copy of the FPGA-generated data, b ut later may be replaced by a combination
		 * of the combined data from 4-sensor (4-pair) correlation and horizontal/vertical pairs only to improve detection of vertical/
		 * horizontal features  
		 * @return line-scan array of per-tile correlation strength by reference (not a copy), so it can be modified
		 */
		public double [] getStrength(){
			double trustedCorrelation = tileProcessor.getTrustedCorrelation();
			if (strength == null){
				strength =  disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
				if (trustedCorrelation > 0.0){
					for (int i = 0; i < strength.length; i++){
						if (Math.abs(disparity_map[disparity_index][i]) > trustedCorrelation) strength[i] = 0.0; // too far 
					}
				}
			}
			return strength;
		}
		/**
		 * Get four pairs (original) correlation strength. Not a copy
		 * @return line-scan array of per-tile horizontal pairs correlation strength by reference (not a copy)
		 */
		public double [] getOriginalStrength(){
			return disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX];
		}
		/**
		 * Get horizontal pairs correlation strength for vertical features. Not a copy
		 * @return line-scan array of per-tile horizontal pairs correlation strength by reference (not a copy)
		 */
		public double [] getHorStrength(){
			double trustedCorrelation = tileProcessor.getTrustedCorrelation();
			if (strength_hor == null) {
				strength_hor = disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH].clone();
				if (trustedCorrelation > 0.0){
					for (int i = 0; i < strength_hor.length; i++){
						if (Math.abs(disparity_map[ImageDtt.DISPARITY_INDEX_HOR][i]) > trustedCorrelation) strength_hor[i] = 0.0; // too far 
					}
				}
			}
			return strength_hor;
		}
		/**
		 * Get veriical pairs correlation strength for horizontal features. Not a copy
		 * @return line-scan array of per-tile horizontal pairs correlation strength by reference (not a copy)
		 */
		public double [] getVertStrength(){
			double trustedCorrelation = tileProcessor.getTrustedCorrelation();
			if (strength_vert == null) {
				strength_vert = disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH].clone();
				if (trustedCorrelation > 0.0){
					for (int i = 0; i < strength_vert.length; i++){
						if (Math.abs(disparity_map[ImageDtt.DISPARITY_INDEX_VERT][i]) > trustedCorrelation) strength_vert[i] = 0.0; // too far 
					}
				}
				
			}
			return strength_vert;
		}
		
		/**
		 * Get Get combine per-tile disparity values from correlation combined with pre-programmed initial disparity shift.
		 * @return line-scan array of per-tile disparity by reference (not a copy), so it can be modified
		 */
		
		public double [] getDisparity() // get calculated combo disparity
		{
			return getDisparity(0);
		}
		/**
		 * Get one of the line-scan per-tile correlation data.
		 * @param mode 0 - final data (initially copy FPGA generated 4-pair correation)
		 *             1 - original FPGA generated 4-sensor correlation
		 *             2 - 2 - horizontal pairs correlation, detecting vertical features
		 *             3 - 2 - vertical pairs correlation, detecting horizontal features
		 * @return line-scan array of per-tile disparity by reference (not a copy), so it can be modified
		 */

		public double [] getDisparity(int mode) // mode = 0 - normal disparity, 1 - hor, 2 - vert
		{
			if (calc_disparity == null) conditionDisparity();
			switch (mode) {
			case 1: return calc_disparity;
			case 2: return calc_disparity_hor;
			case 3: return calc_disparity_vert;
			default: if (calc_disparity_combo == null) calc_disparity_combo = calc_disparity.clone();
				return calc_disparity_combo;
			}
		}
		
		// methods to "condition" measured disparity values
		public void conditionDisparity()
		{
			conditionDisparity(disparity_index);
		}

		public void conditionDisparity(int disparity_index)
		{
			int tilesX = tileProcessor.getTilesX();
			int tilesY = tileProcessor.getTilesY();
			double corr_magic_scale = tileProcessor.getMagicScale();


			this.disparity_index = disparity_index;
			calc_disparity =      new double[tilesY*tilesX];
			calc_disparity_hor =  new double[tilesY*tilesX];
			calc_disparity_vert = new double[tilesY*tilesX];
			for (int i = 0; i < tilesY; i++){
				for (int j = 0; j < tilesX; j++){
					int indx = i * tilesX + j;
					calc_disparity[indx] =      disparity_map[disparity_index][indx]/corr_magic_scale +               this.disparity[i][j];
					calc_disparity_hor[indx] =  disparity_map[ImageDtt.DISPARITY_INDEX_HOR][indx]/corr_magic_scale +  this.disparity[i][j];
					calc_disparity_vert[indx] = disparity_map[ImageDtt.DISPARITY_INDEX_VERT][indx]/corr_magic_scale + this.disparity[i][j];
				}
			}
			calc_disparity_combo = calc_disparity.clone(); // for now - just clone, can be modified separately and combined with hor/vert
		}
		
		/**
		 * Replaces current combo disparity for tiles that are weak and do not have any neighbor within disparity range from this one
		 * @param selection optional boolean mask of tiles to use/update
		 * @param weakStrength maximal strength of the tile to be considered weak one
		 * @param maxDiff maximal difference from the most similar neighbor to be considered an outlayer
		 * @param disparityFar minimal acceptable disparity for weak tiles
		 * @param disparityNear maximal acceptable disparity for weak tiles
		 * @return mask of weak (replaced) tiles
		 * 
		 * Replace weak by a weighted average of non-weak. If there are none - use weak ones, including this one too. 
		 */
		public boolean[] replaceWeakOutlayers(
				final boolean [] selection,
				final double weakStrength,    // strength to be considered weak, subject to this replacement
				final double maxDiff,
				final double maxDiffPos,      // Replace weak outlayer tiles that have higher disparity than weighted average
				final double maxDiffNeg,      // Replace weak outlayer tiles that have lower disparity than weighted average
				final double disparityFar,
				final double disparityNear,
				final int debugLevel)
		{
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			
			final int nTiles = tilesX*tilesY;
			final boolean [] weakOutlayers = new boolean [nTiles];
			int [] dirs8 = {-tilesX,  -tilesX + 1, 1, tilesX +1, tilesX, tilesX - 1, -1, -tilesX - 1};
			final int [] dirs = dirs8;
			final double [] disparity = getDisparity(0);
			final double [] strength =  getStrength();
			final double absMinDisparity = 0.5 * disparityFar; // adjust? below this is definitely wrong (weak) 
			final double absMaxDisparity = 1.5 * disparityNear; // change?
			final int dbg_nTile = (debugLevel > 0) ? 43493: -1; // x=77,y=134; // 42228; // x = 108, y = 130 46462; // 41545;
			final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
			// first pass = find outlayers
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (((strength[nTile] < weakStrength) ||
									(disparity[nTile] < absMinDisparity) ||
									(disparity[nTile] > absMaxDisparity))&& ((selection == null) || selection[nTile])) {
								if (nTile == dbg_nTile){
									System.out.println("replaceWeakOutlayers():1 nTile="+nTile);
								}
								double [] dbg_disparity = disparity;
								double dbg_disparity_nTile = disparity[nTile];
								double dbg_disparityFar = disparityFar;
								double dbg_disparityNear = disparityNear;
								boolean [] dbg_weakOutlayers = weakOutlayers;
								int tileY = nTile / tilesX;
								int tileX = nTile % tilesX;
								if ((tileY > 0) && (tileY < (tilesY -1)) &&(tileX > 0) && (tileX < (tilesX -1))){ // disregard outer row/cols
									weakOutlayers[nTile] = true;
									boolean hasNeighbors = false;
									double sd = 0.0, sw = 0.0;
									for (int dir = 0; dir< dirs.length; dir++){
										int nTile1 = nTile + dirs[dir];
										double dbg_disparity_nTile1 = disparity[nTile1];
										if (((selection == null) || selection[nTile1]) &&
												 (disparity[nTile1] >= disparityFar) && // don't count on too near/too far for averaging
												 (disparity[nTile1] <= disparityNear)){
											double w = strength[nTile1];
											sw += w;
											sd += w * disparity[nTile1];
											hasNeighbors = true;
											if (Math.abs(disparity[nTile]-disparity[nTile1]) <= maxDiff){ // any outlayer - will be false
												weakOutlayers[nTile] = false;
//												break;
											}
										}
									}
									if (sw >= 0.0) {
										sd /= sw;
										if      (disparity[nTile] < (sd - maxDiffNeg)) weakOutlayers[nTile] = true;
										else if (disparity[nTile] > (sd + maxDiffPos)) weakOutlayers[nTile] = true;
									}
									if (disparity[nTile] < disparityFar)  weakOutlayers[nTile] = true;
									if (disparity[nTile] > disparityNear) weakOutlayers[nTile] = true;
									if (!hasNeighbors) {
										weakOutlayers[nTile] = false; // lone tile or NaN among NaNs
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			
			// second pass - replace outlayers
			final double [] src_disparity = disparity.clone();
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (nTile == dbg_nTile){
								System.out.println("replaceWeakOutlayers():2 nTile="+nTile);
							}
							if (weakOutlayers[nTile]) {
								double sw = 0.0, sd = 0.0;
								for (int dir = 0; dir< dirs.length; dir++){
									int nTile1 = nTile + dirs[dir];
									if (!weakOutlayers[nTile1] && ((selection == null) || selection[nTile1 ]) ) {
										double w = strength[nTile1];
										sw += w;
										sd += w * src_disparity[nTile1]; 
									}
								}
								if (sw == 0) { // Nothing strong around - repeat with weak and this one too.
									double w = strength[nTile];
									if (!Double.isNaN( src_disparity[nTile])) {
										sw += w;
										sd += w * src_disparity[nTile];
									}
									for (int dir = 0; dir< dirs.length; dir++){
										int nTile1 = nTile + dirs[dir];
										if ((selection == null) || selection[nTile1 ]) {
											w = strength[nTile1];
											if (!Double.isNaN( src_disparity[nTile1])) {
												sw += w;
												sd += w * src_disparity[nTile1]; 
											}
										}
									}
								}
								if (sw > 0) { // should be, do nothing if not
									disparity[nTile] = sd/sw;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return weakOutlayers;
		}
		
		
		public SuperTiles setSuperTiles(
				double     step_near,
				double     step_far,
				double     step_threshold,
				double     min_disparity,
				double     max_disparity,
				double     strength_floor,
				double     strength_pow,
				double     stBlurSigma)
		{
			this.superTiles = new SuperTiles(
					this,
					step_near,
					step_far,
					step_threshold,
					min_disparity,
					max_disparity,
					strength_floor,
					strength_pow,
					stBlurSigma);
			return this.superTiles;
		}
		public double [] showDisparityHistogram()
		{
			if (this.superTiles == null){
				return null;
			}
			return this.superTiles.showDisparityHistogram();
		}
		
		public double [] showDisparityHistogram(double [][] dispHist)
		{
			if (this.superTiles == null){
				return null;
			}
			return this.superTiles.showDisparityHistogram(dispHist);
		}
		
		
		public int showDisparityHistogramWidth()
		{
			return this.superTiles.showDisparityHistogramWidth();
		}
		
		public double [][][] getMaxMinMax(){
			if (this.superTiles == null){
				return null;
			}
			return superTiles.getMaxMinMax();
		}
		
		
		public double [] showMaxMinMax(){
			if (this.superTiles == null){
				return null;
			}
			return this.superTiles.showMaxMinMax();
		}
		public int getNumBins(){
			if (this.superTiles == null){
				return 0;
			}
			return superTiles.numBins;
		}
		
		public double[] getSuperTileStrength()
		{
			if (this.superTiles == null){
				return null;
			}
			return superTiles.stStrength;
		}

		public double [][] getBgDispStrength()
		{
			if ((bgTileDisparity == null) || (bgTileStrength == null)){
				double [][] rslt = {bgTileDisparity,bgTileStrength};
				return rslt;
			}
			return getBgDispStrength(
					this.minBgDisparity, 
					this.minBgFract);
		}
		
		
		public double [][] getBgDispStrength(
				final double minBgDisparity, 
				final double minBgFract)
		{
			if (superTiles == null){
				return null;
			}
			if ((minBgDisparity != this.minBgDisparity) || (minBgFract != this.minBgFract)){
				this.minBgDisparity = minBgDisparity;
				this.minBgFract = minBgFract; 
				superTiles.bgDisparity = null; // per super-tile
				superTiles.bgStrength = null; // per super-tile
				bgTileDisparity = null; // per tile
				bgTileStrength = null; // per tile
			}
			if ((superTiles.bgDisparity == null) || (superTiles.bgStrength == null)){
				if (superTiles.getBgDispStrength(
						minBgDisparity, 
						minBgFract) == null) {
					superTiles.bgDisparity = null; // per super-tile
					superTiles.bgStrength = null; // per super-tile
					bgTileDisparity = null; // per tile
					bgTileStrength = null; // per tile
					return null; // failed
				}
				// now lap-combine supertiles, get this.* from superTiles.*
				
				double [][] bgTileDispStrength = superTiles.getBgTileDispStrength();
				bgTileDisparity = bgTileDispStrength[0];
				bgTileStrength =  bgTileDispStrength[1];
			}
			double [][] rslt = {bgTileDisparity,bgTileStrength};
			return rslt;
		}
		
		public double [] getBgDisparity(){
			return bgTileDisparity;
		}
		public double [] getBgStrength(){
			return bgTileStrength;
		}

	} // end of class CLTPass3d
