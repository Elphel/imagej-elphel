package com.elphel.imagej.tileprocessor;
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
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

public class CLTPass3d{
//	static double  max_overexposed = 0.8; // TODO: make parameter
		public   double [][]    disparity; // per-tile disparity set for the pass[tileY][tileX]
		public   int    [][]    tile_op;   // what was done in the current pass
		private  double [][]    disparity_sav; // saved disparity
		private  int    [][]    tile_op_sav;   // saved tile_op
		public   double [][]    disparity_map =            null; // add 4 layers - worst difference for the port
		public   double [][]    lazy_eye_data =            null;
		public   int            lma_cluster_size =         -1;
		public   boolean []     lazy_eye_force_disparity = null;
		public   int []         num_tile_max =             null; // number of maximums in a tile

		double []               calc_disparity = null; // composite disparity, calculated from "disparity", and "disparity_map" fields

		                                       // using horizontal features and corr_magic_scale
		// used directly in TileProcessor.compositeScan()
		double []               calc_disparity_hor =   null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		double []               calc_disparity_vert =  null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		double []               calc_disparity_combo = null; // composite disparity, calculated from "disparity", and "disparity_map" fields
		double []               strength =             null; // composite strength, initially uses a copy of raw 4-sensor correlation strength
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
		private boolean []      border_tiles =         null; // these are border tiles, zero out alpha
		private boolean []      selected =             null; // which tiles are selected for this layer
		boolean []              has_lma =              null; // which tiles are measured with LMA (reliable strength)
		public  double [][][][] texture_tiles;
		// texture_selection is only used for the GPU and if not null means it is for the GPU
		public boolean []       texture_selection =    null; // use by the GPU to set texture to generate
				
		public  double [][]     max_tried_disparity =  null; //[ty][tx] used for combined passes, shows maximal disparity for this tile, regardless of results
		public  int             clust_radius =         0;    // 0 - no clusters, 1 = 1x1 (same as 0, just different program), 2 - 3x3, 3 - 5x5
		public  boolean         is_combo =             false;
		public  boolean         is_measured =          false;
		public  String          texture = null; // relative (to x3d) path
		public  Rectangle       texture_bounds; // in tiles, not pixels !
		public  int             dbg_index;
		public  int             disparity_index =      ImageDtt.DISPARITY_INDEX_CM; // may also be ImageDtt.DISPARITY_INDEX_POLY
		public  int             lma_disparity_index =  ImageDtt.DISPARITY_INDEX_POLY; // set to -1 to ignore and always use just CM (also applies to lma_strength - next)
		public double [][]      tiles_RBGA =           null;

		SuperTiles              superTiles = null;
		final TileProcessor     tileProcessor;
		public int getNumSensors() {
			return tileProcessor.getNumSensors();
		}
		
		public CLTPass3d (TileProcessor tileProcessor)
		{
			this.tileProcessor = tileProcessor;
		}

		public CLTPass3d (TileProcessor tileProcessor, int mode)
		{
			this.tileProcessor = tileProcessor;
			switch (mode){
			case 0:
				tile_op =   new int [tileProcessor.getTilesY()][tileProcessor.getTilesX()];
				disparity = new double [tileProcessor.getTilesY()][tileProcessor.getTilesX()];
				break;

			}
		}
		public int getClustRadius() {
			return clust_radius;
		}
		
		public void setClustRadius (int clust_radius) {
			this.clust_radius = clust_radius;
		}
		
		public TileProcessor getTileProcessor()
		{
			return this.tileProcessor;
		}

		public double [][][][] getTextureTiles()
		{
			return 	texture_tiles;
		}
		// texture_selection is only used for the GPU and if not null means it is for the GPU
		public boolean [] getTextureSelection()
		{
			return texture_selection;
		}
		
		public void setTextureSelection(boolean [] selection) {
			texture_selection = selection;
		}

		public void setTextureSelection(int indx, boolean sel) {
			texture_selection[indx] = sel;
		}

		public double [][] getTilesRBGA()
		{
			return tiles_RBGA;
		}

		public void setTilesRBGA(double [][] rgba)
		{
			tiles_RBGA = rgba;
		}
		
		/*
		public float [][] getTextureImages()
		{
			return 	texture_img;
		}

		public Rectangle getTextureWoi()
		{
			return 	texture_woi;
		}
		*/
		
		public double [][] getMaxTriedDisparity()
		{
			return max_tried_disparity;
		}
		
		
		public double [][] getTileRBGA(
				int num_layers) // 4 or 12
		{
//			if (texture_img != null) {
//				System.out.println("FIXME: implement replacement for the GPU-generated textures (using macro mode?)");
//			}
			double [][] tones = getTilesRBGA();
			if (tones != null) {
				int nl = tones.length;
				if (nl > num_layers) {
					double [][] tones1 = new double [num_layers][tones[0].length];
					for (int n = 0; n < num_layers; n++) {
						tones1[n] = tones[n];
					}
					return tones1;
				} else {
					return tones;
				}
			}
			System.out.println("FIXME: should not get here, tones should be calculated earlier");
			
			if (texture_tiles == null) return null;
			int tilesY = texture_tiles.length;
			int tilesX = 0;
			int nl = 0;
			for (int ty = 0; ty < tilesY; ty++){
				if (texture_tiles[ty] != null){
					tilesX = texture_tiles[ty].length;
					for (int tx = 0; tx < tilesX; tx++){
						if (texture_tiles[ty][tx] != null){
							nl = texture_tiles[ty][tx].length;
							break;
						}
					}
					if (nl > 0) break;
				}
				if (nl > 0) break;
			}
			if (num_layers > nl) num_layers = nl;
			int numTiles = tilesX * tilesY;
			double [] scales = new double [num_layers];
			for (int n = 0; n < num_layers; n++){
				if       (n < 3)  scales[n] = 1.0/255.0; // R,B,G
				else if  (n == 3) scales[n] = 1.0; //alpha
				else if  (n < 8)  scales[n] = 1.0; // ports 0..3
				else              scales[n] = 1.0/255.0; // RBG rms, in 1/255 units, but small
			}
			double [][] tileTones = new double [num_layers][numTiles];
			for (int ty = 0; ty < tilesY; ty++ ) if (texture_tiles[ty] != null){
				for (int tx = 0; tx < tilesX; tx++ ) if (texture_tiles[ty][tx] != null) {
					int indx = ty * tilesX + tx;
					for (int n = 0; n < num_layers; n++) if (texture_tiles[ty][tx][n] != null){
						double s = 0.0;
						for (int i = 0; i < texture_tiles[ty][tx][n].length; i++){
							s += texture_tiles[ty][tx][n][i];
						}
						s /= (texture_tiles[ty][tx][n].length/4); // overlapping tiles
						s *= scales[n];
						tileTones[n][indx] = s;
					}
				}
			}
			return tileTones;
		}

		public String getTextureName()
		{
			if (texture != null) {
				return texture;
			} else {
				return "null-texture-name";
			}
		}


		// Will not work if texture is disabled
		public  void            updateSelection(){ // add updating border tiles? -- only for CPU!
			int tilesX = tileProcessor.getTilesX();
			int tilesY = tileProcessor.getTilesY();
			selected = new boolean[tilesY*tilesX];
			int minX = tilesX, minY = tilesY, maxX = -1, maxY = -1;
			/*
			if (texture_img != null) { // using GPU output
//tileProcessor.getTileSize()
				if (texture_woi != null) {
					int tile_size = tileProcessor.getTileSize();
					texture_bounds = new Rectangle(
							texture_woi.x/tile_size, texture_woi.y/tile_size, texture_woi.width/tile_size, texture_woi.height/tile_size);
					// setting full rectangle as selected, not just textures? Use some other method?
					for (int ty = texture_bounds.y; ty < (texture_bounds.y + texture_bounds.height); ty++) {
						for (int tx = texture_bounds.x; tx < (texture_bounds.x + texture_bounds.width); tx++) {
							selected[ty*tilesX+tx] = true;
						}
					}
					return;
				}
			}
			*/
			if (texture_tiles != null) {
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
			}
			if (maxX < 0) {
				texture_bounds = null;
			} else {
				texture_bounds = new Rectangle(minX, minY, maxX - minX +1, maxY - minY +1 );
			}
		}
		public  Rectangle  getTextureBounds(){
			return texture_bounds;
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

			if (disparity_map != null) fixNaNDisparity();
			calc_disparity =       null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			calc_disparity_hor =   null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			calc_disparity_vert =  null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			calc_disparity_combo = null; // composite disparity, calculated from "disparity", and "disparity_map" fields
			strength =             null; // composite strength, initially uses a copy of raw 4-sensor correleation strength
			strength_hor =         null; // updated hor strength, initially uses a copy of raw measured
			strength_vert =        null; // updated hor strength, initially uses a copy of raw measured
			bgTileDisparity =      null;
			bgTileStrength =       null;
//			border_tiles =         null; // these are border tiles, zero out alpha
//			selected =             null; // which tiles are selected for this layer
			superTiles =           null;
			setTilesRBGA(null);

		}
		/**
		 * Get FPGA-calculated per-tile maximal differences between the particular image and the average one.
		 * @return per-camera sesnor array of line-scan differences
		 */

		public double [][] getDiffs (){
			if (disparity_map == null) return null;
			double  [][] these_diffs =    new double[getNumSensors()][];
			for (int i = 0; i< these_diffs.length; i++) these_diffs[i] = disparity_map[ImageDtt.IMG_DIFF0_INDEX + i]; // IMG_DIFF0_INDEX does not depend on num sensors
			return these_diffs;
		}

		public void resetCalc(){ // only needed if the same task was reused
			calc_disparity = null;
			strength =       null;
			strength_hor =   null;
			strength_vert =  null;
			superTiles =     null;

		}

		public boolean [] getLazyEyeForceDisparity() {
			return lazy_eye_force_disparity;
		}

		public void setLazyEyeForceDisparity(boolean [] lazy_eye_force_disparity) {
			this.lazy_eye_force_disparity = lazy_eye_force_disparity;
		}

		public double [][] getLazyEyeData() {
			return lazy_eye_data;
		}

		public void setLazyEyeData(double [][] lazy_eye_data) {
			this.lazy_eye_data = lazy_eye_data;
		}

		public int getLazyEyeClusterSize() {
			return lma_cluster_size;
		}
		public void setLazyEyeClusterSize(int lma_cluster_size) {
			this.lma_cluster_size = lma_cluster_size;
		}



		public boolean [] getSelected(){
			return selected;
		}
		
		public boolean [] getSelectedOrTileOp(){ // to use instead of getSelected if it is null
			if (selected != null) return selected;
			if (tile_op== null) return null;
			boolean [] to_selected = new boolean [tile_op.length * tile_op[0].length];
			int indx = 0;
			for (int ty = 0; ty < tile_op.length; ty++) {
				for (int tx = 0; tx < tile_op[ty].length; tx++) {
					to_selected[indx++] = tile_op[ty][tx] != 0;
				}
			}
			return to_selected;
		}
		

		public boolean [] getBorderTiles(){
			return this.border_tiles;
		}

		public void setSelected (boolean [] selected) {
			this.selected = selected;
		}
		public void setBorderTiles (boolean [] border_tiles) {
			this.border_tiles = border_tiles;
		}
		
		public boolean [] getLMA() {
			if (has_lma == null) {
				has_lma = hasLMADefined();
			}
			return has_lma;
		}

		public void setLMA(boolean [] has_lma) {// use for combo tiles
			this.has_lma = has_lma;
		}

		public void resetLMA() {
			this.has_lma = null;
		}
		public void setLMA(double [] disparity_lma) {// use for combo tiles
			this.has_lma = new boolean [disparity_lma.length];
			for (int i = 0; i < disparity_lma.length; i++) {
				this.has_lma[i] = !Double.isNaN(disparity_lma[i]);
			}
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
			if ((disparity == null) || (strength == null)) {
				return;
			}
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

		public double [] getOverexposedFraction(){
			return (disparity_map != null)? disparity_map[ImageDtt.OVEREXPOSED] : null;
		}


		/**
		 * Returns per-tile correlation "strength". Initially - copy of the FPGA-generated data, but later may be replaced by a combination
		 * of the combined data from 4-sensor (4-pair) correlation and horizontal/vertical pairs only to improve detection of vertical/
		 * horizontal features
		 * @return line-scan array of per-tile correlation strength by reference (not a copy), so it can be modified
		 */
		public double [] getStrength(){
			double trustedCorrelation = tileProcessor.getTrustedCorrelation();
			double max_overexposure = tileProcessor.getMaxOverexposure();
			if (strength == null){
				strength =  disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
				if (trustedCorrelation > 0.0){
					for (int i = 0; i < strength.length; i++){
						if (Math.abs(disparity_map[disparity_index][i]) > trustedCorrelation) strength[i] = 0.0; // too far
					}
				}
				double [] overexposed = disparity_map[ImageDtt.OVEREXPOSED];
				if ((max_overexposure > 0.0) && (overexposed != null)){
					for (int i = 0; i < strength.length; i++){
						if (overexposed[i] > max_overexposure) strength[i] = 0.0; // too overexposed
					}
				}
			}
			return strength;
		}
		
		public void setStrength(double [] strength) {
			this.strength = strength;
		}
		
		/**
		 * Get four pairs (original) correlation strength. Not a copy
		 * @return line-scan array of per-tile horizontal pairs correlation strength by reference (not a copy)
		 */
		public double [] getOriginalStrength(){
			if (disparity_map != null) return disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX];
			else return getStrength(); // after replacing with rig data
		}
		/**
		 * Get horizontal pairs correlation strength for vertical features. Not a copy
		 * @return line-scan array of per-tile horizontal pairs correlation strength by reference (not a copy)
		 */
		public double [] getHorStrength(){
			double trustedCorrelation = tileProcessor.getTrustedCorrelation();
			double max_overexposure = tileProcessor.getMaxOverexposure();
			if (strength_hor == null) {
				if (disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH] == null) {
					return null;
				}
				strength_hor = disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH].clone();
				if (trustedCorrelation > 0.0){
					for (int i = 0; i < strength_hor.length; i++){
						if (Math.abs(disparity_map[ImageDtt.DISPARITY_INDEX_HOR][i]) > trustedCorrelation) strength_hor[i] = 0.0; // too far
					}
				}
				double [] overexposed = disparity_map[ImageDtt.OVEREXPOSED];
				if ((max_overexposure > 0.0) && (overexposed != null)){
					for (int i = 0; i < strength_hor.length; i++){
						if (overexposed[i] > max_overexposure) strength_hor[i] = 0.0; // too overexposed
					}
				}

			}
			return strength_hor;
		}
		/**
		 * Get vertical pairs correlation strength for horizontal features. Not a copy
		 * @return line-scan array of per-tile horizontal pairs correlation strength by reference (not a copy)
		 */
		public double [] getVertStrength(){
			double trustedCorrelation = tileProcessor.getTrustedCorrelation();
			double max_overexposure = tileProcessor.getMaxOverexposure();
			if (strength_vert == null) {
				if (disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH] == null) {
					return null;
				}
				strength_vert = disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH].clone();
				if (trustedCorrelation > 0.0){
					for (int i = 0; i < strength_vert.length; i++){
						if (Math.abs(disparity_map[ImageDtt.DISPARITY_INDEX_VERT][i]) > trustedCorrelation) strength_vert[i] = 0.0; // too far
					}
				}
				double [] overexposed = disparity_map[ImageDtt.OVEREXPOSED];
				if ((max_overexposure > 0.0) && (overexposed != null)){
					for (int i = 0; i < strength_hor.length; i++){
						if (overexposed[i] > max_overexposure) strength_vert[i] = 0.0; // too overexposed
					}
				}

			}
			return strength_vert;
		}

		/**
		 * Combine per-tile disparity values from correlation combined with pre-programmed initial disparity shift.
		 * @return line-scan array of per-tile disparity by reference (not a copy), so it can be modified
		 */

		public double [] getDisparity() // get calculated combo disparity
		{
			return getDisparity(0);
		}
		/**
		 * Get one of the line-scan per-tile correlation data.
		 * @param mode 0 - final data (initially copy FPGA generated 4-pair correlation)
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
		public double [] getDisparityLMA() {
			if (calc_disparity == null) conditionDisparity(); // sets has_lma;
			double [] disparityLMA = calc_disparity.clone();
			for (int i = 0; i < disparityLMA.length; i++) {
				if (!has_lma[i]) {
					disparityLMA[i] = Double.NaN;
				}
			}
			return disparityLMA;
		}
		public int [] getNumTileMax() {
			if (num_tile_max == null) {
				if (disparity_map[ImageDtt.DISPARITY_VARIATIONS_INDEX] == null) {
					return null;
				}
				num_tile_max = new int [disparity_map[ImageDtt.DISPARITY_VARIATIONS_INDEX].length];
				for (int i = 0; i < num_tile_max.length; i++) {
					///https://stackoverflow.com/questions/6341896/cast-int-to-double-then-back-to-int-in-java
					num_tile_max[i] = (int)	disparity_map[ImageDtt.DISPARITY_VARIATIONS_INDEX][i];
				}
			}
			return num_tile_max;
		}
		
		public double getDisprityHor(int nt) {
			return (disparity_map[ImageDtt.DISPARITY_INDEX_HOR] == null) ? Double.NaN : disparity_map[ImageDtt.DISPARITY_INDEX_HOR][nt]; 
		}
		public double getDisprityVert(int nt) {
			return (disparity_map[ImageDtt.DISPARITY_INDEX_VERT] == null) ? Double.NaN : disparity_map[ImageDtt.DISPARITY_INDEX_VERT][nt]; 
		}
		public double getStrengthHor(int nt) {
			return (disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH] == null) ? 0.0 : disparity_map[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH][nt]; 
		}
		public double getStrengthVert(int nt) {
			return (disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH] == null) ? 0.0 : disparity_map[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH][nt]; 
		}
		
		
		
		// methods to "condition" measured disparity values
		public void conditionDisparity()
		{
			int tilesX = tileProcessor.getTilesX();
			int tilesY = tileProcessor.getTilesY();
			double corr_magic_scale =     tileProcessor.getMagicScale();
			double corr_magic_scale_LMA = 1.0;
			calc_disparity =      new double[tilesY*tilesX];
			calc_disparity_hor =  new double[tilesY*tilesX];
			calc_disparity_vert = new double[tilesY*tilesX];
			Arrays.fill(calc_disparity,      Double.NaN);
			Arrays.fill(calc_disparity_hor,  Double.NaN);
			Arrays.fill(calc_disparity_vert, Double.NaN);
			
			double [] lma_disparity = (lma_disparity_index >= 0) ? disparity_map[lma_disparity_index] : null;
			double [] lma_strength =  (lma_disparity_index >= 0) ? disparity_map[lma_disparity_index+1] : null;
			has_lma = (lma_disparity != null) ? new boolean[tilesX * tilesY] : null;

			for (int i = 0; i < tilesY; i++){
				for (int j = 0; j < tilesX; j++){
					int indx = i * tilesX + j;
					if ((lma_disparity != null) && !Double.isNaN(lma_disparity[indx]) && (lma_strength[indx] > 0.0)) {
						calc_disparity[indx] =  lma_disparity[indx]/corr_magic_scale_LMA +                            this.disparity[i][j];
						has_lma[indx] = true;
					} else {
						calc_disparity[indx] =  disparity_map[disparity_index][indx]/corr_magic_scale +               this.disparity[i][j];
					}
				}
			}
			calc_disparity_combo = calc_disparity.clone(); // for now - just clone, can be modified separately and combined with hor/vert
			if (disparity_map[ImageDtt.DISPARITY_INDEX_HOR] != null) {
				for (int i = 0; i < tilesY; i++){
					for (int j = 0; j < tilesX; j++){
						int indx = i * tilesX + j;
						calc_disparity_hor[indx] =  disparity_map[ImageDtt.DISPARITY_INDEX_HOR][indx]/corr_magic_scale +  this.disparity[i][j];
					}
				}
			}
			if (disparity_map[ImageDtt.DISPARITY_INDEX_VERT] != null) {
				for (int i = 0; i < tilesY; i++){
					for (int j = 0; j < tilesX; j++){
						int indx = i * tilesX + j;
						calc_disparity_vert[indx] = disparity_map[ImageDtt.DISPARITY_INDEX_VERT][indx]/corr_magic_scale + this.disparity[i][j];
					}
				}
			}
		}
		
		public boolean [] hasLMADefined(){ // will try not to create this.has_lma
			if (disparity_map == null) {
				if (has_lma==null) {
					return null;
				} else {
					return has_lma;
				}
				
			}
			int tilesX = tileProcessor.getTilesX();
			int tilesY = tileProcessor.getTilesY();
			double [] lma_disparity = (lma_disparity_index >= 0) ? disparity_map[lma_disparity_index] : null;
			double [] lma_strength =  (lma_disparity_index >= 0) ? disparity_map[lma_disparity_index+1] : null;
			boolean [] lma_defined = new boolean[tilesX * tilesY];
			if (lma_disparity != null) {
				for (int i = 0; i < lma_disparity.length; i++) {
					lma_defined[i] = !Double.isNaN(lma_disparity[i]) && (lma_strength[i] > 0);
				}
			}
			return lma_defined;
		}
		
		
		// bypassing calculations
		public void setCalcDisparityStrength(
				double [] disparity,
				double [] strength) {
			this.calc_disparity = disparity;
			this.strength =       strength;
			calc_disparity_combo = calc_disparity.clone();
		}

		/**
		 * Replaces current combo disparity for tiles that are weak and do not have any neighbor within disparity range from this one
		 * @param selection optional boolean mask of tiles to use/update
		 * @param weakStrength maximal strength of the tile to be considered weak one
		 * @param maxDiff maximal difference from the most similar neighbor to be considered an outlier
		 * @param disparityFar minimal acceptable disparity for weak tiles
		 * @param disparityNear maximal acceptable disparity for weak tiles
		 * @return mask of weak (replaced) tiles
		 *
		 * Replace weak by a weighted average of non-weak. If there are none - use weak ones, including this one too.
		 */
		public boolean[] replaceWeakOutliers( // does not replace tiles with LMA available, busts LMA-defined strengths when averaging
				final boolean [] selection,
				final double weakStrength,    // strength to be considered weak, subject to this replacement
				final double maxDiff,
				final double maxDiffPos,      // Replace weak outlier tiles that have higher disparity than weighted average
				final double maxDiffNeg,      // Replace weak outlier tiles that have lower disparity than weighted average
				final double disparityFar,
				final double disparityNear,
				final int debugLevel)
		{
			final double  scale_strength_lma = 5.0; // increase LMA-defined strength during averaging
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();

			final int nTiles = tilesX*tilesY;
			final boolean [] weakOutliers = new boolean [nTiles];
			int [] dirs8 = {-tilesX,  -tilesX + 1, 1, tilesX +1, tilesX, tilesX - 1, -1, -tilesX - 1};
			final int [] dirs = dirs8;
			final double [] disparity = getDisparity(0);
			final double [] strength =  getStrength();
			final boolean [] has_lma =  getLMA();

			final double absMinDisparity = 0.5 * disparityFar; // adjust? below this is definitely wrong (weak)
			final double absMaxDisparity = 1.5 * disparityNear; // change?
			final int dbg_nTile = (debugLevel > 0) ? 43493: -1; // x=77,y=134; // 42228; // x = 108, y = 130 46462; // 41545;
			final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
			// first pass = find outliers
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (((!has_lma[nTile] && (strength[nTile] < weakStrength)) ||
									(disparity[nTile] < absMinDisparity) ||
									(disparity[nTile] > absMaxDisparity))&& ((selection == null) || selection[nTile])) {
								if (nTile == dbg_nTile){
									System.out.println("replaceWeakOutliers():1 nTile="+nTile);
								}
								double [] dbg_disparity = disparity;
								double dbg_disparity_nTile = disparity[nTile];
								double dbg_disparityFar = disparityFar;
								double dbg_disparityNear = disparityNear;
								boolean [] dbg_weakOutliers = weakOutliers;
								int tileY = nTile / tilesX;
								int tileX = nTile % tilesX;
								if ((tileY > 0) && (tileY < (tilesY -1)) &&(tileX > 0) && (tileX < (tilesX -1))){ // disregard outer row/cols
									weakOutliers[nTile] = true;
									boolean hasNeighbors = false;
									double sd = 0.0, sw = 0.0;
									for (int dir = 0; dir< dirs.length; dir++){
										int nTile1 = nTile + dirs[dir];
										double dbg_disparity_nTile1 = disparity[nTile1];
										if (((selection == null) || selection[nTile1]) &&
												 (disparity[nTile1] >= disparityFar) && // don't count on too near/too far for averaging
												 (disparity[nTile1] <= disparityNear)){
											double w = strength[nTile1];
											if (has_lma[nTile1]) {
												w *= scale_strength_lma;
											}
											sw += w;
											sd += w * disparity[nTile1];
											hasNeighbors = true;
											if (Math.abs(disparity[nTile]-disparity[nTile1]) <= maxDiff){ // any outlier - will be false
												weakOutliers[nTile] = false;
//												break;
											}
										}
									}
									if (sw >= 0.0) {
										sd /= sw;
										if      (disparity[nTile] < (sd - maxDiffNeg)) weakOutliers[nTile] = true;
										else if (disparity[nTile] > (sd + maxDiffPos)) weakOutliers[nTile] = true;
									}
									if (disparity[nTile] < disparityFar)  weakOutliers[nTile] = true;
									if (disparity[nTile] > disparityNear) weakOutliers[nTile] = true;
									if (!hasNeighbors) {
										weakOutliers[nTile] = false; // lone tile or NaN among NaNs
									}
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);

			// second pass - replace outliers
			final double [] src_disparity = disparity.clone();
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (nTile == dbg_nTile){
								System.out.println("replaceWeakOutliers():2 nTile="+nTile);
							}
							if (weakOutliers[nTile]) {
								double sw = 0.0, sd = 0.0;
								for (int dir = 0; dir< dirs.length; dir++){
									int nTile1 = nTile + dirs[dir];
									if (!weakOutliers[nTile1] && ((selection == null) || selection[nTile1 ]) ) {
										double w = strength[nTile1];
										if (has_lma[nTile1]) {
											w *= scale_strength_lma;
										}
										sw += w;
										sd += w * src_disparity[nTile1];
									}
								}
								if (sw == 0) { // Nothing strong around - repeat with weak and this one too.
									double w = strength[nTile];
									if (has_lma[nTile]) {
										w *= scale_strength_lma;
									}
									if (!Double.isNaN( src_disparity[nTile])) {
										sw += w;
										sd += w * src_disparity[nTile];
									}
									for (int dir = 0; dir< dirs.length; dir++){
										int nTile1 = nTile + dirs[dir];
										if ((selection == null) || selection[nTile1 ]) {
											w = strength[nTile1];
											if (has_lma[nTile1]) {
												w *= scale_strength_lma;
											}
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
			return weakOutliers;
		}

		public boolean [] getUntestedBackgroundBorder (
				final boolean    [] known,
				final double     [] disparity,
				final double        grow_disp_step,
				final int     debugLevel)
		{
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			final int num_tiles = tilesX * tilesY;
			final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
			final boolean [] untested_bgnd = new boolean [num_tiles];
			final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
							if (known[nTile]){
								int tX = nTile % tilesX;
								int tY = nTile / tilesX;
								double max_disp = max_tried_disparity[tY][tX] + grow_disp_step;
								for (int dir = 0; dir < 8; dir++){
									int nTile1 = tnImage.getNeibIndex(nTile, dir);
									if ((nTile1 >=0) && known[nTile1] && (disparity[nTile1] > max_disp)) {
										untested_bgnd[nTile] = true;
										break;
									}
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			return untested_bgnd;
		}

		public boolean [] measuredTiles ()
		{
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			final int num_tiles = tilesX * tilesY;
			final boolean [] measured = new boolean [num_tiles];
			final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
								int tX = nTile % tilesX;
								int tY = nTile / tilesX;
								measured[nTile] = tile_op[tY][tX] != 0;
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			return measured;
		}

		public void saveTileOpDisparity()
		{
			disparity_sav =disparity.clone();
			for (int i = 0; i < disparity.length; i++ ) if (disparity[i]!= null) disparity_sav[i] = disparity[i].clone();
			tile_op_sav =tile_op.clone();
			for (int i = 0; i < tile_op.length; i++ ) if (tile_op[i]!= null) tile_op_sav[i] = tile_op[i].clone();
		}
		public void restoreTileOpDisparity()
		{
			disparity = disparity_sav;
			tile_op = tile_op_sav;
		}
		public void restoreKeepTileOpDisparity()
		{
			restoreTileOpDisparity();
			saveTileOpDisparity();
		}
		public int setTileOpDisparity(
				double []  disparity) {
			boolean [] selection = new boolean [disparity.length];
			for (int i = 0; i < disparity.length; i++) {
				selection[i] = !Double.isNaN(disparity[i]);
			}
			return setTileOpDisparity(
					selection, // boolean [] selection,
					disparity); // double []  disparity)
		}

		public int setTileOpDisparity(
				boolean [] selection,
				double []  disparity)
		{
			int op = ImageDtt.setImgMask(0, 0xf);
			op =     ImageDtt.setPairMask(op,0xf);
			op =     ImageDtt.setForcedDisparity(op,true);
			return setTileOpDisparity(
					op,         // int        tile_op,
					selection,  // boolean [] selection,
					disparity); // double []  disparity)
		}

		public int setTileOpDisparity(
				int        tile_op,
				boolean [] selection,
				double []  disparity)
		{
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			this.disparity =   new double [tilesY][tilesX];
			if (selection != null) {
				this.tile_op =     new int [tilesY][tilesX];
			}
			int num_op_tiles = 0;

			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				int indx =  tilesX * ty + tx;
				if (selection == null) {
					this.disparity[ty][tx] = (disparity == null)? 0.0: disparity[indx];
					if (this.tile_op[ty][tx] != 0){
						num_op_tiles ++;
					}
				} else {
					if (selection[indx]) {
//						this.disparity[ty][tx] = (disparity == null)? 0.0: disparity[indx];
						this.tile_op[ty][tx] = tile_op;
						num_op_tiles ++;
					} else {
//						this.disparity[ty][tx] = 0.0;
						this.tile_op[ty][tx] = 0;
					}
					this.disparity[ty][tx] = (disparity == null)? 0.0: disparity[indx];
					
				}
			}
			return num_op_tiles;
		}
		
		public void setTileOp(
				int  tile_op)
		{
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			this.tile_op =     new int [tilesY][tilesX];

			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				this.tile_op[ty][tx] = tile_op;
			}
		}

		public double [] getDA() {
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			double [] da = new double [tilesX*tilesY];
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				da[tilesX * ty + tx] = this.disparity[ty][tx];
			}
			return da;
		}

		public void setDA(double [] da) {
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			for (int ty = 0; ty < tilesY; ty++) for (int tx = 0; tx <tilesX; tx++){
				this.disparity[ty][tx] = da[tilesX * ty + tx];
			}
		}
		
		/**
		 * Set next measurement disparity from last calculated
		 */
		public void updateDisparity()
		{
			if (disparity_map != null) {
				setTileOpDisparity(null, getDisparity(0));
			}
		}


		public double [] getSecondMaxDiff (
				final boolean averaged)
		{
			final double [][] diffs = getDiffs();
			if (diffs == null) return null;
			for (int i = 0; i < diffs.length; i++) {
				if (diffs[i] == null) return null;
			}
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			final int num_tiles = tilesX * tilesY;
			final double [] second_max = new double [num_tiles];
			final boolean [] measured =  measuredTiles ();
			final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
							int imax1 = 0;
							for (int ip = 1; ip < diffs.length; ip++){
								if (diffs[ip][nTile] > diffs[imax1][nTile]) imax1 = ip;
							}
							int imax2 = (imax1 == 0)? 1 : 0;
							for (int ip = 0; ip < diffs.length; ip++) if (ip != imax1) {
								if (diffs[ip][nTile] > diffs[imax2][nTile]) imax2 = ip;
							}
							second_max[nTile] = diffs[imax2][nTile];
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			if (!averaged) return second_max;
			final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
			final double [] second_max_averaged = new double [num_tiles];
			final double [] dir_weights = {1.0/16, 1.0/8, 1.0/16, 1.0/8, 1.0/16, 1.0/8, 1.0/16, 1.0/8, 1.0/4};
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (measured[nTile]) {
							double sw = 0.0;
							double swd = 0.0;
							for (int dir = 0; dir < 9; dir++){ // including 8 - center
								int nTile1 = tnImage.getNeibIndex(nTile, dir);
								if ((nTile1 >=0) && measured[nTile1]) {
									sw +=  dir_weights[dir];
									swd += dir_weights[dir] * second_max[nTile1] ;
								}
							}
							second_max_averaged[nTile] = swd/sw;
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			return second_max_averaged;
		}



		// same, but 2 steps around
		public boolean [] getUntestedBackgroundBorder2 (
				final boolean    [] known,
				final double     [] disparity,
				final double        grow_disp_step,
				final int     debugLevel)
		{
			final int tilesX = tileProcessor.getTilesX();
			final int tilesY = tileProcessor.getTilesY();
			final int num_tiles = tilesX * tilesY;
			final TileNeibs tnImage = new TileNeibs(tilesX, tilesY); // num_tiles/tilesX);
			final boolean [] untested_bgnd = new boolean [num_tiles];
			final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
							if (known[nTile]){
								int tX = nTile % tilesX;
								int tY = nTile / tilesX;
								double max_disp = max_tried_disparity[tY][tX] + grow_disp_step;
								for (int dir = 0; dir < 24; dir++){
									int nTile1 = tnImage.getNeibIndex2(nTile, dir);
									if ((nTile1 >=0) && known[nTile1] && (disparity[nTile1] > max_disp)) {
										untested_bgnd[nTile] = true;
										break;
									}
								}
							}
						}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			return untested_bgnd;
		}






		public SuperTiles getSuperTiles()
		{
			return this.superTiles;
		}

		public SuperTiles setSuperTiles(
				double     step_near,
				double     step_far,
				double     step_threshold,
				double     min_disparity,
				double     max_disparity,
				double     stBlurSigma,
				boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				MeasuredLayersFilterParameters mlfp,
				int        measSel)
		{
			this.superTiles = new SuperTiles(
					this,
					step_near,
					step_far,
					step_threshold,
					min_disparity,
					max_disparity,
					stBlurSigma,
					smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					mlfp,
					measSel);
			return this.superTiles;
		}
		public double [] showDisparityHistogram(
				double [][][][] disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
				boolean [][]    tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all

				boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				MeasuredLayersFilterParameters mlfp,
				int        measSel)
		{
			if (this.superTiles == null){
				return null;
			}
			return this.superTiles.showDisparityHistogram(
					disparity_strength, // pre-calculated disparity/strength [per super-tile][per-measurement layer][2][tiles] or null
					tile_sel, // null  or per-measurement layer, per-tile selection. For each layer null - do not use, {} - use all

					smplMode, //        = true;   // Use sample mode (false - regular tile mode)
					mlfp,
					measSel);
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
