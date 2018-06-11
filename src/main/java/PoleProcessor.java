/**
 ** PoleProcessor - handling poles like street lights and other vertical objects
 ** common in artificial environments
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  PoleProcessor.java is free software: you can redistribute it and/or modify
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
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;


public class PoleProcessor {
	final static double THRESHOLD_LIN = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
	final static double THRESHOLD_QUAD = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
	BiCamDSI biCamDSI;
	ArrayList<PoleCluster> pole_clusters = new ArrayList<PoleCluster>();
	int tilesX;
	int tilesY;

	class PoleCluster{
		double mean_disparity =   0.0;
		double target_disparity = 0.0;
		double target_previous =  0.0;
		double mean_x =    0.0;
		double mean_y =    0.0;
		double strength =  0.0;
		ArrayList <Integer> tiles = new ArrayList<Integer>();
		Rectangle bBox = new Rectangle();
		Rectangle eBox; // not calculated, assigned by caller
		ArrayList<DsPair> ds_list =     new ArrayList<DsPair>();
		ArrayList<DsPair> ds_selected = new ArrayList<DsPair>();
		int    layer =     -1;
		double  xc_center = Double.NaN;
		double  slope =  Double.NaN; // dx/dy
		boolean disabled = false;
//		int     split = -1;
		boolean [] pole_selection; // tiles in eBox that belong to the pole
		double [][] bg_ds; // background disparity/strength for eBox dimensions;
		class DsPair{
			int    nTile;
			double disparity;
			double strength;
			DsPair(
					int nTile,
					double disparity,
					double strength)
			{
				this.nTile = nTile;
				this.disparity = disparity;
				this.strength = strength;
			}
		}

		double  getMeanDisparity()   {return mean_disparity;}
		double  getTargetDisparity() {return target_disparity;}
		void    setTargetDisparity(double disparity) {target_disparity = disparity;}

		double  getStrength() {return strength;} // to get average strength - need to divide by getNumTiles()
		double  getMeanX()    {return mean_x;}
		double  getMeanY()    {return mean_y;}
		int     getNumTiles() {return tiles.size();}
		ArrayList<Integer>  getTiles() {return tiles;}
		public Rectangle getBox(boolean extended) {return extended? eBox: bBox;}
		public Rectangle getBBox() {return bBox;}
		public Rectangle getEBox() {return eBox;}
		public void trimEBoxBottom(int height) {eBox.height = height;}
		public boolean [] getMask () {return pole_selection;}
		public boolean    getMask (int dx, int dy) {return pole_selection[dx + dy * eBox.width];}

		@Override
		public String toString() {
			return String.format("{%3d, %3d, %2d, %2d}: %7.5f",	eBox.x, eBox.y, eBox.width,eBox.height, target_disparity);
		}
		//-1 - definitely not,+1 - definitely yes, 0 - uncertain
		public int poleFilter() {


//			int min_mask = 10; // good - min 12
			int min_height = 10; // good - min 12
			boolean headlessOK = false;
//			double min_fraction = 0.3;    // worst 0.3478
			double min_fraction = 0.5;    // worst 0.3478 , next - 0.57
			double min_frac_height = 0.5;    // worst 0.25 , next - 0.556
			double min_frac_box =    0.7;    // worst 0.75 , next - 0.556

/*
				boolean unselected_only = true;
				double bg_strength = cluster.getAverageBgStrength(unselected_only);
				double meas_strength = cluster.getAverageStrength(true); // divided by measured tiles only
				double mask_strength = cluster.getAverageStrength(false); // divided by a total number of pole tiles (by mask)
				double rstrength = mask_strength/bg_strength;
 */



//			if (getMaskSize() < min_mask) return false;

// return 0;
//			if (!disabled) return 0;
			if (disabled)                            return -1;
			if (eBox.height < min_height)            return -1;
			if (!headlessOK && (getNumTiles() == 0)) return -1; // skip headless clusters
			if (getFractSelected() < min_fraction)   return -1;
			if (getFractHeight() < min_frac_height)  return -1;
			if (getFractBoxHeight() < min_frac_box)  return -1;
			return 0;
		}

// get cluster stats for filtering
		public int getMaskSize() {
			int n =0;
			for (int i = 0; i< pole_selection.length; i++) if (pole_selection[i]) n++;
			return n;
		}
		public int getSelectionSize() {
			return ds_selected.size();
		}

		public double getFractSelected() {
			return 1.0* ds_selected.size()/getMaskSize();
		}

		public double getAverageStrength(boolean measured_only) { // denominator - measured - only average among measured tiles, false - all pole mask
			double sw = 0.0;
			for (DsPair dsp: ds_selected) {
				sw += dsp.strength;
			}
			return sw / (measured_only? ds_selected.size() : getMaskSize());
		}

		public double getAverageBgStrength( // empty counts as zero strengths
				boolean unselected)
		{
			double sw = 0.0;
			int ns = 0;
			for (int indx = 0; indx < bg_ds.length; indx++) {
				if (!unselected || !pole_selection[indx]) {
					ns++;
					double [] ds = bg_ds[indx];
					if (ds != null){
						sw += bg_ds[indx][1];
					}
				}
			}
			return sw/ns;
		}

		public boolean [][] getVertMasks(){
			boolean [] y_mask =     new boolean[eBox.height];
			boolean [] y_selection = new boolean[eBox.height];
			for (int i = 0; i < pole_selection.length; i++) {
				if (pole_selection[i]) {
					y_mask[i/eBox.width] = true;
				}
			}
			for (int dy = 0; dy < (bBox.y-eBox.y); dy++){ // don't count all what is above the seed
				y_mask[dy] = false;
			}
			// remove seed (original head) area
			for (int nt:tiles) {
				int dy = nt / tilesX - eBox.y;
				y_mask[dy] = false;
			}

			for (DsPair dsp:ds_selected) {
				y_selection[dsp.nTile / tilesX - eBox.y] = true;
			}
			boolean [][] masks = {y_selection, y_mask}; //{ num, denom}
			return masks;
		}

		public double getFractHeight() {
			boolean [][] masks = getVertMasks(); // {y_selection, y_mask};
			int num_y_mask = 0, num_y_sel = 0;
			for (int i = 0; i < eBox.height; i++) {
				if (masks[1][i]) {
					num_y_mask++;
					if (masks[0][i]) num_y_sel++;
				}
			}
			return (num_y_mask > 0) ? (1.0 * num_y_sel / num_y_mask) : 0.0;
		}
		public double getFractBoxHeight() {
			boolean [][] masks = getVertMasks(); // {y_selection, y_mask};
			int num_y_sel = 0;
			int y_top = bBox.x-eBox.x;
			int y_bot = 0;
			for (int i = y_top; i < eBox.height; i++) {
				if (masks[1][i]) {
					y_bot = i;
				}
				if (masks[0][i]) num_y_sel++;
			}
//			return 1.0 * num_y_sel / (eBox.height - y_top);
			return 1.0 * num_y_sel / (y_bot- y_top + 1);
		}

		// get longest missing stem over the "visible" stem mask
		public int getFreeStemGap() {
			boolean [][] masks = getVertMasks(); // {y_selection, y_mask};
			int gs = 0;
			int ln = masks[1].length;
			int max_gap = 0;
			int gap_len  = 0;
			for (; gs < ln -1;) {
				for (; gs < ln-1; gs++) {
					if (masks[1][gs] && !masks[0][gs]) {
						gap_len = 1;
						break; // gs now points at the beginning of the gap
					}
				}
				if (gs >= ln -1) {
					break; // too late
				}
				int ge;
				for (ge = gs+1; ge < ln; ge++) {
					if (masks[1][ge] && !masks[0][ge]) {
						gap_len++;
					} else if (masks[0][ge]){
						break; // ge now points at the detected tile (first after the gap)
					}
				}
				if (gap_len > max_gap) {
					max_gap = gap_len;
				}
				gs = ge;
			}
			return max_gap;
		}

		public double getMDisparityRMS() {
			double sw = 0.0, swd = 0.0, swd2 = 0.0;
			for (DsPair dsp:ds_selected) {
				double w = dsp.strength;
				sw += w;
				double d = dsp.disparity;
				swd+=w*d;
				swd2 += w*d*d;
			}
			return Math.sqrt((swd2 * sw + swd*swd)/(sw*sw));
		}

		public boolean checkRightLeft(
				double min_strength,
				double disp_adiff,
				double disp_rdiff,
				double [][] norn_ds)
		{
			int extra_margin = 2;
			int debugLevel = -1;
			int dbg_tileX = 128; // 221;// 128;
			int dbg_tileY = 87; // 130; // 87;
			if (debugLevel > -2) {
				if (eBox.contains(new Point(dbg_tileX, dbg_tileY))) {
					System.out.println("checkRightLeft(): eBox={"+eBox.x+", "+eBox.y+", "+eBox.width+", "+eBox.height+"}");
					System.out.println("checkRightLeft(): eBox={"+eBox.x+", "+eBox.y+", "+eBox.width+", "+eBox.height+"}");
				}
			}

			double disp_diff = disp_adiff + target_disparity*disp_rdiff;
			double disp_margin = mean_disparity - disp_diff; // need to have tiles both to right and left from the head with disparity below this margin (or be too weak)
			double [] hor_profile = new double[eBox.width+2*extra_margin];
			TileNeibs  tnImage = biCamDSI.tnImage;
			int nTile = eBox.x + eBox.y*tilesX;
			for (int dx = -extra_margin; dx<eBox.width+extra_margin; dx++) {
				for (int dy = (bBox.y-eBox.y); dy < (bBox.y-eBox.y+bBox.height); dy++) {
					int nTile1 = tnImage.getNeibIndex(nTile, dx, dy);
					if (nTile1 >= 0) {
						if ((norn_ds[1][nTile1] > min_strength) && (norn_ds[1][nTile1] > hor_profile[dx+extra_margin])) {
							hor_profile[dx+extra_margin] = norn_ds[0][nTile1];
						}
					}
				}

			}
			double y = bBox.y-eBox.y + 0.5*bBox.height - 0.5*eBox.height;
			double x = xc_center + y* slope + 0.5*eBox.width;
			int x0 = (int) Math.round(x);
			boolean sep = false;
			if ((x0 < 0) || (x0 > hor_profile.length )) {
				System.out.println(String.format("Disabling cluster {%d, %d, %d, %d} disp= %f as it is not separated on right or left. x0= %d",
						eBox.x, eBox.y, eBox.width, eBox.height, target_disparity, x0));
				return false;
			}
			for (int i = x0+extra_margin; i <  hor_profile.length; i++) {
				if  (hor_profile[i] <  disp_margin) {
					sep = true;
					break;
				}
			}
			if (!sep) {
				System.out.println(String.format("Disabling cluster {%d, %d, %d, %d} disp= %f as it is not separated on right. x0= %d",
						eBox.x, eBox.y, eBox.width, eBox.height, target_disparity, x0));
				return false;
			}
			sep = false;
			for (int i = x0+extra_margin; i >= 0; i--) {
				if  (hor_profile[i] <  disp_margin) {
					sep = true;
					break;
				}
			}
			if (!sep) {
				System.out.println(String.format("Disabling cluster {%d, %d, %d, %d} disp= %f as it is not separated on left. x0= %d",
						eBox.x, eBox.y, eBox.width, eBox.height, target_disparity, x0));

			}

			return sep;
		}






		public double [] getPoleLine() {
			if (Double.isNaN(xc_center)) {
				return null;
			}
			double [] pole_line = {xc_center,slope};
			return pole_line;
		}

		public void resetDsPairs() {
			ds_list =     new ArrayList<DsPair>();
			ds_selected = new ArrayList<DsPair>();
		}
		public void addDs(
				int nTile,
				double disparity,
				double strength)
		{
//			ds_list.add(new DsPair(nTile,disparity,strength));
			addDs(new DsPair(nTile,disparity,strength));
		}

		public void addDs(
				DsPair dsp)
		{
			ds_list.add(dsp);
			if (pole_selection != null) {
				int dx = dsp.nTile % tilesX - eBox.x;
				int dy = dsp.nTile / tilesX - eBox.y;
				int indx = dy*eBox.width + dx;
				if ((indx < pole_selection.length) && pole_selection[indx]) {
					ds_selected.add(dsp);
				}
			}
		}


		public void removeDs(
				int nTile)
		{
			// may have multiple for the same tile - remove all
			for (int nt = 0; nt < ds_list.size(); nt++) {
				if (ds_list.get(nt).nTile == nTile) {
					ds_list.remove(nt);
					nt--; // compensate for loop increment
				}
			}
		}


		public double [][] getDsPairs(boolean selected){
			double [][] ds = new double [2][selected?ds_selected.size():ds_list.size()];
			for (int i = 0; i < ds[0].length; i++) {
				DsPair ds_pair = selected? (ds_selected.get(i)):(ds_list.get(i));
				ds[0][i] = ds_pair.disparity;
				ds[1][i] = ds_pair.strength;
			}
			return ds;
		}

		public int [] getDsTiles(boolean selected){
			int[] ds_tiles = new int[selected?ds_selected.size():ds_list.size()];
			for (int i = 0; i < ds_tiles.length; i++) {
				DsPair ds_pair = selected? (ds_selected.get(i)):(ds_list.get(i));
				ds_tiles[i] = ds_pair.nTile;
			}
			return ds_tiles;
		}

		public void setEBox(
				int ext_left,
				int ext_right,
				int ext_up,
				int ext_down)
		{
			this.eBox = new Rectangle(
					bBox.x - ext_left,
					bBox.y - ext_up,
					bBox.width +  ext_right + ext_left,
					bBox.height + ext_up +   ext_down);
			if (eBox.x < 0) {
				eBox.width+=eBox.x;
				eBox.x = 0;
			}
			if (eBox.y < 0) {
				eBox.height +=eBox.y;
				eBox.y = 0;
			}
			if ((eBox.x + eBox.width) >= tilesX) {
				eBox.width = tilesX - eBox.x;
			}
			if ((eBox.y + eBox.height) >= tilesY) {
				eBox.height = tilesY - eBox.y;
			}
		}
		public boolean intersects(PoleCluster cluster) {
			return this.eBox.intersects(cluster.eBox);
		}

		public void setLayer(int layer) {this.layer = layer;}

		public int getLayer() {return layer;}

		void addCluster(
				double [][] norm_ds,
				PoleCluster cluster)
		{
			for (int otherTile:cluster.getTiles()) {
				addTile(norm_ds, otherTile);
			}
			// merge eBox-es if any
			if ((eBox != null) && (cluster.eBox !=null)) {
				eBox.add(cluster.eBox);
			}
		}


		void addTile(
				double [][] norm_ds,
				int nTile)
		{
			if (tiles.indexOf(nTile) >=0) {
				System.out.println("*** BUG in PoleCluster.addTile("+nTile+") - tile already exists here ***");
				return; //
			}
			if (norm_ds[1][nTile] <= 0) {
				System.out.println("*** BUG in PoleCluster.addTile("+nTile+") - strength= "+norm_ds[1][nTile]+" ***");
				return; //
			}
			TileNeibs  tnImage = biCamDSI.tnImage;
			int tileX = nTile % tnImage.sizeX;
			int tileY = nTile / tnImage.sizeX;
			mean_disparity = (mean_disparity * strength + norm_ds[0][nTile] * norm_ds[1][nTile])/(strength + norm_ds[1][nTile]);
			mean_x =    (mean_x *    strength + tileX *             norm_ds[1][nTile])/(strength + norm_ds[1][nTile]);
			mean_y =    (mean_y *    strength + tileY *             norm_ds[1][nTile])/(strength + norm_ds[1][nTile]);
			strength += norm_ds[1][nTile];
			if (tiles.isEmpty() || ( tileX < bBox.x))  bBox.x = tileX;
			if (tiles.isEmpty() || ((tileX - bBox.x) > (bBox.width -1)))  bBox.width = tileX - bBox.x + 1;
			if (tiles.isEmpty() || ( tileY < bBox.y))  bBox.y = tileY;
			if (tiles.isEmpty() || ((tileY - bBox.y) > (bBox.height -1)))  bBox.height = tileY - bBox.y + 1;
			tiles.add(nTile);
			target_disparity = mean_disparity;
		}
		// remove would be trickier - will need to recalculate bBox
		void removeTile(
				double [][] norm_ds,
				int nTile)
		{
			int indx = tiles.indexOf(nTile);
			if (indx < 0) {
				System.out.println("*** BUG in PoleCluster.removeTile("+nTile+") - tile does not belong here ***");
				return; //
			}
			TileNeibs  tnImage = biCamDSI.tnImage;
			int tileX = nTile % tnImage.sizeX;
			int tileY = nTile / tnImage.sizeX;
			mean_disparity = (mean_disparity * strength - norm_ds[0][nTile] * norm_ds[1][nTile])/(strength - norm_ds[1][nTile]);
			mean_x =         (mean_x *         strength - tileX *             norm_ds[1][nTile])/(strength - norm_ds[1][nTile]);
			mean_y =         (mean_y *         strength - tileY *             norm_ds[1][nTile])/(strength - norm_ds[1][nTile]);
			strength -= norm_ds[1][nTile];
			tiles.remove(indx);
			// re-calculate bounding box
			for (int i = 0; i < tiles.size(); i++ ) {
				int nt = tiles.get(i);
				tileX = nt % tnImage.sizeX;
				tileY = nt / tnImage.sizeX;
				if ((i==0) || ( tileX < bBox.x))  bBox.x = tileX;
				if ((i==0) || ((tileX - bBox.x) > (bBox.width -1)))  bBox.width = tileX - bBox.x + 1;
				if ((i==0) || ( tileY < bBox.y))  bBox.y = tileX;
				if ((i==0) || ((tileY - bBox.y) > (bBox.height -1)))  bBox.height = tileY - bBox.y + 1;
			}
			target_disparity = mean_disparity;
		}
		PoleCluster [] two_clusters = null;
		// see if the cluster needs split. returns -1 if it does not, > 0 if all tiles with offset >= go to the second cluster
		// returns null or a pair of new clusters instead of the current (does not trry to split in more than 2
		public PoleCluster[]  checkAndSplit(
				double [][] norm_ds,
				int min_dist,          // minimal distance
				boolean must_zero) {
			if ((eBox.x == 61) && (eBox.y== 128)) {
				System.out.println("Debugging "+this.toString());
			}
			if (eBox.width < 2) {
				return null;
			}
			double [] histogram = new double [eBox.width];
			for (DsPair dsp:ds_list) {
				histogram [dsp.nTile % tilesX - eBox.x] += dsp.strength;
			}
			boolean [] is_max = new boolean[histogram.length];
			boolean [] is_zero = new boolean[histogram.length];
			int num_lmax = 0;
			for (int i = 0; i < histogram.length;i++) {
				is_zero[i] = histogram[i] == 0.0;
				if (!is_zero[i]) {
					is_max[i] = ((i == 0) || (histogram[i] >= histogram[i-1])) && ((i == (histogram.length - 1)) || (histogram[i] >= histogram[i+1]));
					if (is_max[i]) num_lmax ++;
				}
			}
			if (num_lmax < 2) {
				return null;
			}
			int gs = 0, ge = 0;
			boolean is_gap = false;
			int cut = -1;
			// find first (leftmost maximum) - skip all non-maximums
			for (;gs < is_max.length-1;gs++) { //find start of the non-maximum run
				if (is_max[gs]) break;
			}

			for (; gs < is_max.length-1;) {
				for (;gs < is_max.length-1;gs++) { //find start of the non-maximum run
					if (!is_max[gs]) break;
				}
				if (gs >= is_max.length-1) { // too late to have another maximum
					break;
				}
				boolean has_zero = false;
				ge = gs + 1;
				for (;ge < is_max.length; ge++) { //find end of the non-maximum run (next maximum)
					has_zero |= is_zero[ge];
					is_gap = true;
					if (is_max[ge]) break;
				}
				if ((ge < is_max.length) && is_gap && ((ge - gs) >= min_dist) && (!must_zero || has_zero)) { // gap valid for split - find the cut place
					cut = gs;
					for(int i = gs+1; i < ge; i++) {
						if (histogram[i] < histogram[cut]) {
							cut = i;
						}
					}
					// if a long rung of 0 - find center
					int cut_end = cut+1;
					for (;cut_end < ge;cut_end++) {
						if (histogram[cut_end] > histogram[cut]) break;
					}
					cut = (cut+cut_end+1)/2;
					break;
				}
				is_gap = false; // that was a bad gap
				gs = ge;
			}
			if (cut > 0) {
				// split in 2 clusters, adjust eBox without normal gaps not to change layers
				// will have to use wrong (old) mean_disparity as it was used for measurements
				two_clusters = new PoleCluster [2];
				two_clusters[0]= new PoleCluster();
				two_clusters[1]= new PoleCluster();
				for (int nTile: tiles) {
					int dx = nTile % tilesX - eBox.x;
					two_clusters[(dx < cut)? 0 : 1].addTile(norm_ds, nTile);
				}
				// overwrite mean disparity to that of the original
				two_clusters[0].setTargetDisparity(getTargetDisparity());
				two_clusters[1].setTargetDisparity(getTargetDisparity());
				// manually set eBox-es, so they do not overlap on the same layer. OK, as for vertical objects there is no need for horizontal padding for measurements
				two_clusters[0].eBox = new Rectangle(eBox.x,       eBox.y, cut,              eBox.height);
				two_clusters[1].eBox = new Rectangle(eBox.x + cut, eBox.y, eBox.width - cut, eBox.height);
				// set layer to that of the original one
				two_clusters[0].layer = layer;
				two_clusters[1].layer = layer;
				// split ds_list (but not  ds_selected - it should be calculated after that)
				for (DsPair dsp:ds_list) {
					int dx = dsp.nTile % tilesX - eBox.x;
					two_clusters[(dx < cut)? 0 : 1].addDs(dsp);
				}
				System.out.println("Splitting cluster "+this.toString()+" into 2: "+two_clusters[0]+toString()+" and "+two_clusters[1].toString());

			}
			return two_clusters;
		}
		// approximate pole centerline as offset from eBox.x at eBox.y and slope (dx/dy)
		public void calcPoleLine(
				double max_diff,
				double max_tilt,
				double max_rms,
				int    min_tiles,
				double damp_tilt,
				int debugLevel)
		{
			double x0 = eBox.x + 0.5*eBox.width;
			double y0 = eBox.y + 0.5*eBox.height;
			{
				double sw = 0.0, swx = 0.0;
				for (DsPair dsp:ds_list) {
					int x = dsp.nTile % tilesX;
					double w = dsp.strength;
					sw += w;
					swx += w * x;
				}
				if (sw <=0) {
					return;
				}
				xc_center = swx/sw - x0;
				slope = 0.0;
			}
			// preliminary sort:
			ds_selected = new ArrayList<DsPair>();
			for (DsPair dsp: ds_list) {
				double x = dsp.nTile % tilesX - x0; // relative to the eBox center
				double y = dsp.nTile / tilesX - y0;
				if (Math.abs(x - xc_center) <= (max_diff + Math.abs(y * max_tilt))) {
					ds_selected.add(dsp);
				}
			}
			if (ds_selected.size() < min_tiles) {
				return; // use initial vertical line
			}
			for (int nstep = 0; nstep < ds_selected.size(); nstep++) { // should normally break;
				// calculate rms, and worst tile
				double worst_e2 = 0.0;
				int worst_index = -1;
				double sw = 0.0, swe2 = 0.0;
				for (int indx = 0; indx < ds_selected.size(); indx++) {
					DsPair dsp  = ds_selected.get(indx);
					double x = dsp.nTile % tilesX - x0; // relative to the eBox center
					double y = dsp.nTile / tilesX - y0;
					double w = dsp.strength;
					double e = x - (xc_center + y* slope);
					double e2 = e*e;
					if (e2 > worst_e2) {
						worst_e2 =    e2;
						worst_index = indx;
					}
					sw += w;
					swe2 += w * e2;
				}
				double rms = Math.sqrt(swe2/sw);
				if (nstep > 0){
					if ((rms < max_rms) || (ds_selected.size() < min_tiles)){
						return;
					}
					// remove worst
					ds_selected.remove(worst_index);
				}
				// perform linear approximation, at least once
				double [] approx_lin = fitLine(
						damp_tilt,    // double damp_tilt,
						ds_selected,        // ArrayList<DsPair> plist,
						debugLevel);  // int debugLevel);
				if (approx_lin == null) {
					if (debugLevel > -1){
						System.out.println("fitLine(): can not find linear approximation (should not happen)");
					}
					return;
				}
				xc_center = approx_lin[0];
				slope =     approx_lin[1];
			}
		}

		public double [] fitLine(
				double damp_tilt,
				ArrayList<DsPair> plist,
				int debugLevel) {
			PolynomialApproximation pa = new PolynomialApproximation();
			// using 2d approximation and duplicate samples
			final double [] damping = {damp_tilt, damp_tilt, 0.0}; // 0.0 will be applied to average value, tilt_cost - to both tilts
			double [][][] mdata = new double [2 * plist.size()][3][];
			int mindx = 0;
			double x0 = eBox.x + 0.5*eBox.width;
			double y0 = eBox.y + 0.5*eBox.height;
			for (DsPair dsp: plist) {
				int tileX = dsp.nTile % tilesX;
				int tileY = dsp.nTile / tilesX;
				double v = tileX - x0;
				double w = dsp.strength;
				double x = tileY - y0;

				for (int y = 0; y < 2; y++) { //
					mdata[mindx][0] = new double [2];
					mdata[mindx][0][0] =  x; // sx - smpl_radius;
					mdata[mindx][0][1] =  y; // sy - smpl_radius;
					mdata[mindx][1] = new double [1];
					mdata[mindx][1][0] =  v; // smpl_d[indxs];
					mdata[mindx][2] = new double [1];
					mdata[mindx][2][0] =  w; // smpl_w[indxs];
					mindx ++;

				}
			}
			double[][] approx2d = pa.quadraticApproximation(
					mdata,
					true,          // boolean forceLinear,  // use linear approximation
					damping,       // double [] damping,
					THRESHOLD_LIN,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
					THRESHOLD_QUAD, // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
					debugLevel);
			if (approx2d == null){
				if (debugLevel > -1){
					System.out.println("fitLine(): can not find linear approximation");
				}
				return null;
			}
			double [] rslt = {approx2d[0][2], approx2d[0][0]};
			return rslt;
		}

		public void copyNormDS(
				double [][] norm_ds)
		{
			TileNeibs  tnImage = biCamDSI.tnImage;
			bg_ds = new double [eBox.width*eBox.height][];
			int nTile = eBox.x + eBox.y * tilesX;
			for (int dy = 0; dy< eBox.height; dy++) {
				for (int dx = 0; dx < eBox.width; dx++) {
					int nTile1 = tnImage.getNeibIndex(nTile, dx, dy);
					if ((nTile1 >=0) && !Double.isNaN(norm_ds[0][nTile1]) && (norm_ds[1][nTile1] > 0.0)){
						double [] ds = {norm_ds[0][nTile1], norm_ds[1][nTile1]};
						bg_ds[dx + dy * eBox.width] = ds;

					}
				}
			}
		}

		public void createPoleMask( // shoud run copyNormDS before
				//				double [][] norm_ds,
				int         min_neibs,
				boolean     use_seed,
				double      width,
				double      disp_aover,
				double      disp_rover,
				int         debugLevel)
		{
			boolean square_head = true; // only allow seed in the top square of eBox
			int dbg_tileX = 99; // 221;// 128;
			int dbg_tileY = 112; // 130; // 87;
			if (debugLevel > -2) {
				if (eBox.contains(new Point(dbg_tileX, dbg_tileY))) {
					System.out.println("createPoleMask(): eBox={"+eBox.x+", "+eBox.y+", "+eBox.width+", "+eBox.height+"}");
				}
			}
			double head_edge_strength = 0.1;
			pole_selection = new boolean [eBox.width*eBox.height];
			double disp_near = target_disparity + disp_aover + target_disparity * disp_rover;
			if (bBox.width > 0) {
				for (int dy = bBox.y - eBox.y; dy < eBox.height; dy++) {
					double x = 0.5 * eBox.width +xc_center + slope * (dy - 0.5*eBox.height);
					for (int dx = (int) Math.round(x - width); dx <= (int) Math.round(x + width); dx++) {
						if ((dx >= 0) && (dx < eBox.width)) {
							int indx = dx + eBox.width*dy;
							if ((bg_ds[indx] == null) || !(bg_ds[indx][0] > disp_near)) {
								pole_selection[indx] = true;
							}
						}
					}
				}

			} else {
				System.out.println("bBox is not initialized,  eBox={"+eBox.x+","+eBox.y+","+eBox.width+","+eBox.height+"}");
				System.out.println("bBox is not initialized");
			}
			TileNeibs tnEBox  = new TileNeibs(eBox.width, eBox.height);
			if (use_seed) { // will not be masked by nearer objects?
				boolean [] seed_mask = new boolean [eBox.width*eBox.height];
				for (int nTile:tiles){
					int dy = nTile / tilesX - eBox.y;
					if (!square_head || (dy < eBox.width)) {
						int dx = nTile % tilesX - eBox.x;
						int indx = dx + eBox.width*dy;
						pole_selection[indx] = true;
						seed_mask[indx] =      true;
					}
				}
				// add tiles 1 tile around seed if it was measured and passed filter
				tnEBox.growSelection(
						4,         // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
						seed_mask, // boolean [] tiles,
						null);    // boolean [] prohibit)
				for (DsPair dsp:ds_list) {
					if (dsp.strength > head_edge_strength) {
						int dx = dsp.nTile % tilesX - eBox.x;
						int dy = dsp.nTile / tilesX - eBox.y;
						int indx = dx + dy*eBox.width;
						if (seed_mask[indx]) {
							pole_selection[indx] = true;
						}
					}
				}
			}
			if (min_neibs > 0) {
				int num_removed = tnEBox.removeFewNeibs(
						pole_selection, // boolean [] selection, // should be the same size
						min_neibs); // int min_neibs)
				if ((debugLevel > -2) && (num_removed > 0)) {
					System.out.println("createPoleMask(): removed "+num_removed+" loosely connected/isolated tile (num_neibs < "+min_neibs+")");
				}
			}
		}


		public void filterByMask(
				boolean trim_bottom)
		{

			if (trim_bottom) {
				int indx;
				for (indx = pole_selection.length-1; indx >=0; indx--) {
					if (pole_selection[indx]) break;
				}
				int height = indx/eBox.width + 1;
				if (height < eBox.height) {
					System.out.println(String.format("filterByMask() Trimming {%3d, %3d, %2d, %2d}:%7.5f to height %2d",
							eBox.x, eBox.y, eBox.width,eBox.height, target_disparity, height));
					xc_center -= 0.5 * slope * (eBox.height - height);
					eBox.height = height;
					boolean []  mask =      new boolean[eBox.width * eBox.height];
					double [][] trimmed_ds = new double[eBox.width * eBox.height][];
					for (int i = 0; i < mask.length; i++) {
						mask[i] =        pole_selection[i];
						trimmed_ds[i] =  bg_ds[i];
					}
					pole_selection = mask;
					bg_ds =    trimmed_ds;
				}
			}
			for (int np = 0; np< ds_selected.size(); np++) {
				DsPair dsp = ds_selected.get(np);
				int dx = dsp.nTile % tilesX - eBox.x;
				int dy = dsp.nTile / tilesX - eBox.y;
				int indx = dx + dy*eBox.width;
				if ((indx >= pole_selection.length) || !pole_selection[indx]) {
					ds_selected.remove(np);
					np--; // to compensate for loop
				}
			}
		}

		public double getTargetDiff() {
			return target_disparity - target_previous;
		}
		public double applyMeasuredDisparity(
				double disparity_scale,  // target disparity to differential disparity scale (baseline ratio) (~0.2)
				double diff_power,       // bias towards higher disparities - (disparity+offset) is raised to this power and applied to weight
				                         // if 0.0 - do not apply value to weight
				double diff_offset,      // add to measured differential disaprity before raising to specified power
				int    cut_bottom,       // cut few tile rows from the very bottom - they may be influenced by ground objects
				double keep_bottom,      // do not cut more that this fraction of the bounding bow height
				int debugLevel)          // debug level
		{
			int cut = (int) Math.round(Math.min(eBox.height * keep_bottom, cut_bottom));
			int height = eBox.height - cut;
			double sw=0.0, swd = 0.0;
			for (DsPair dsp:ds_selected) {
				int dx = dsp.nTile % tilesX - eBox.x;
				int dy = dsp.nTile / tilesX - eBox.y;
				int indx = dx + dy*eBox.width;
				if (pole_selection[indx] && (dy < height)) {
					double w = dsp.strength;
					double d = dsp.disparity;
					if (diff_power > 0.0) {
						double v=d + diff_offset;
						if (v < 0.0) v = 0;
						w *= Math.pow(v, diff_power);
					}
					sw += w;
					swd += w * d;
				}
			}
			swd /= sw;
			double diff = disparity_scale * swd;
			target_previous = target_disparity;
			target_disparity += diff;
			if (debugLevel > 0) {
				System.out.println("applyMeasuredDisparity() for cluster "+this.toString()+" target disparity change: "+diff);
			}
			return diff;
		}




	} // end of PoleCluster class



	public PoleProcessor (
			BiCamDSI   biCamDSI,
			int tilesX,
			int tilesY)
	{
		this.biCamDSI = biCamDSI;
		this.tilesX = tilesX;
		this.tilesY = tilesY;
	}

	public double [][] conditionDisparityStrength(
			final BiScan biScan,
			final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,
			final double     strength_pow
			) {
		final double     strength_floor = trusted_strength * strength_rfloor;
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] cond_ds = new double [2] [num_tiles];
		final double [][] disparity_strength = biScan.getDisparityStrength(
				false, // only_strong,
				false, // only_trusted,
				true); // only_enabled);
//		final AtomicInteger ai_num_seeds = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
						if (!Double.isNaN(disparity_strength[0][nTile]) &&(disparity_strength[1][nTile] > strength_floor)){
							double w = disparity_strength[1][nTile] - strength_floor;
							if (strength_pow != 1.0) {
								w = Math.pow(w, strength_pow);
							}
							cond_ds[0][nTile] = disparity_strength[0][nTile];
							cond_ds[1][nTile] = w;

						} else {
							cond_ds[0][nTile] = Double.NaN;
							cond_ds[1][nTile] = 0.0;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return cond_ds;
	}




	public boolean [] findSeeds(
			final double [][] norm_ds,
			final double      min_strength,     // after subtracting floor and possible applying pow
			final int         seed_down,
			final double      seed_aover,
			final double      seed_rover,
			final double      max_disparity
			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final boolean [] seeds = new boolean [num_tiles];
//		final AtomicInteger ai_num_seeds = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (norm_ds[1][nTile] > min_strength){
						double d_head = norm_ds[0][nTile];
						if ((d_head <= max_disparity) && (d_head > seed_aover)) {
							double sw = 0.0;
							double swd = 0.0;
							int nTile1 = nTile; // tnImage.getNeibIndex(indx, TileNeibs.DIR_DOWN);
							for (int i = 0; i < seed_down; i++) {
								nTile1 = tnImage.getNeibIndex(nTile1, TileNeibs.DIR_DOWN);
								if ((nTile1 >= 0) && (norm_ds[1][nTile1] > 0.0) && !Double.isNaN(norm_ds[0][nTile1])) {
									double w = norm_ds[1][nTile1];
									sw += w;
									swd += w * norm_ds[0][nTile1];
								}
							}
							if (sw > 0.0){
								swd /= sw; // weighted average of the disparities right below
								double threshold = swd + seed_aover + seed_rover * d_head;
								if (d_head >= threshold) {
									seeds[nTile] = true;
//									ai_num_seeds.getAndIncrement();
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return seeds;
	}

	public ArrayList<PoleCluster> initPoleClusters(
			final double max_dx,
			final double max_dy,
			final double max_dd,
			boolean []   bseeds,
			double [][]  norm_ds,
			final int    debugLevel)
	{
		final TileNeibs  tnImage = biCamDSI.tnImage;
		ArrayList<PoleCluster> clusters = new ArrayList<PoleCluster>();
		for (int nTile = 0; nTile < bseeds.length; nTile++) if (bseeds[nTile]){
			double this_disp = norm_ds[0][nTile];
			int tileX = nTile % tnImage.sizeX;
			int tileY = nTile / tnImage.sizeX;
			label_gotit: {
				for (PoleCluster cluster: clusters) {
					if (
							!(Math.abs(cluster.getMeanDisparity() - this_disp) > max_dd) && // to allow Double.NaN
							!(Math.abs(cluster.getMeanX() - tileX) > max_dx) &&
							!(Math.abs(cluster.getMeanY() - tileY) > max_dy)) {
						cluster.addTile(norm_ds, nTile);
						break label_gotit;
					}
				}
				// did not find, need to create a new cluster
				PoleCluster cluster = new PoleCluster();
				cluster.addTile(norm_ds, nTile);
				clusters.add(cluster);
			}
		}
		sortClustersByDisparity(clusters);
		return clusters;
	}

	public void extendClusters(
			int ext_left,
			int ext_right,
			int ext_up,
			int ext_down,
			ArrayList<PoleCluster> clusters)
	{
		for (PoleCluster cluster: clusters) {
			// Assign extended bounding box to each cluster
			cluster.setEBox(
					ext_left,
					ext_right,
					ext_up,
					ext_down);
		}
	}


	public void cutClusterBottoms(
			final double [][] norm_ds,
			final double      pedestal_strength,     // bottom tiles should be at least this strong (normalized strength)
			final double      pedestal_disp_over,         // bottom tiles should have disparity at least by this above cluster
			final ArrayList<PoleCluster> clusters)
	{

		// can be multithreaded
		final TileNeibs  tnImage = biCamDSI.tnImage;
//		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						Rectangle bBox = cluster.getBBox(); // never extends beyond ?
						Rectangle eBox = cluster.getEBox();
						int nTile = eBox.x + tilesX * eBox.y;
						double disp_marg = cluster.getMeanDisparity() + pedestal_disp_over;
						for (int dy = bBox.y - eBox.y + bBox.height; dy < eBox.height; dy++) {
							// label
							label_still_some_far_weak:{
							for ( int dx = 0; dx < eBox.width;  dx++) {
								int nTile1 = tnImage.getNeibIndex(nTile, dx, dy); // should always be inside image
								if (    (norm_ds[1][nTile1] < pedestal_strength) ||
										(norm_ds[0][nTile1] < disp_marg)) {
									break label_still_some_far_weak;
								}
							}
							// found pedestal - trim eBox
							cluster.trimEBoxBottom(dy);
							break;
						}

						}


					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}


	public int removeFilteredClusters(
			ArrayList<PoleCluster> clusters,
			int debugLevel)
	{
		int original_size = clusters.size();
		for (int nClust = 0; nClust < clusters.size(); nClust++) {
			PoleCluster cluster = clusters.get(nClust);
			if (cluster.poleFilter() < 0) {
				if (debugLevel > -2) {
					System.out.println("Removing filtered out cluster: "+cluster.toString());
				}
				clusters.remove(nClust);
				nClust--;
			}
		}
		return original_size - clusters.size();
	}



	public int assignClustersToLayers(
			ArrayList<PoleCluster> clusters)
	{

		ArrayList<ArrayList<PoleCluster>> layer_lists = new ArrayList<ArrayList<PoleCluster>>();

		for (PoleCluster cluster: clusters) {
			// Assign extended bounding box to each cluster
			// see if this cluster fits to any of the existing layers
			label_layers:{
				for (int nlayer = 0; nlayer < layer_lists.size(); nlayer++) {
					ArrayList<PoleCluster> layer = layer_lists.get(nlayer);
					// see if current cluster intersects with any on this layer
					label_clusters: {
						for (PoleCluster other_cluster: layer){
							if (cluster.intersects(other_cluster)) {
								break label_clusters; // this layer won't work
							}
						}
						// Layer OK
						cluster.setLayer(nlayer);
						layer.add(cluster);
						break label_layers;
					}
				}
				// need new layer
				cluster.setLayer(layer_lists.size());
				ArrayList<PoleCluster> new_layer = new ArrayList<PoleCluster>();
				new_layer.add(cluster);
				layer_lists.add(new_layer);
			}
		}
		return layer_lists.size();
	}

	public int mergeOverlappingClusters(
			boolean     extended,
			double      disp_tolerance,
			double [][] norm_ds,
			ArrayList<PoleCluster> clusters,
			int debugLevel)
	{
		int num_merge = 0;
		for (int nCluster = 0; nCluster < clusters.size(); nCluster++) {
			PoleCluster cluster = clusters.get(nCluster);
			double this_disp = cluster.getMeanDisparity();
			Rectangle this_box = extended?cluster.getEBox():cluster.getBBox();
			for (int nOther = 0; nOther < nCluster; nOther++) {
				PoleCluster other_cluster = clusters.get(nOther);
				double other_disp = other_cluster.getMeanDisparity();
				Rectangle other_box = extended?other_cluster.getEBox():other_cluster.getBBox();
				if (this_box.intersects(other_box) && (Math.abs(other_disp - this_disp) < disp_tolerance)) {
					if (debugLevel > -2) {
						System.out.println("Merging cluster  " +nCluster+ " to cluster " +nOther);
						System.out.println("Merging cluster   ("+this_box.x+","+this_box.y+","+this_box.width+","+this_box.height+") disp="+this_disp+
								" to cluster ("+other_box.x+","+other_box.y+","+other_box.width+","+other_box.height+") disp="+other_disp);
					}
					other_cluster.addCluster(norm_ds, cluster);
					clusters.remove(nCluster);
					nCluster--; // so when incremented by for() it will stay the same
					num_merge++;
					if (debugLevel > -2) {
						other_disp = other_cluster.getMeanDisparity();
						other_box =  extended ? other_cluster.getEBox(): other_cluster.getBBox();
						System.out.println("Combined cluster: ("+other_box.x+","+other_box.y+","+other_box.width+","+other_box.height+") disp="+other_disp);
					}
					break;
				}
			}
		}
		return num_merge;
	}

	public int  splitClusters(
			double [][] norm_ds,
			int min_dist,          // minimal distance between local maximums
			boolean must_zero,     // zero histogram should exist between local maximum to split
			ArrayList<PoleCluster> clusters,
			int debugLevel)
	{
		int num_split = 0;
		for (int nClust = 0; nClust < clusters.size(); nClust++){
			PoleCluster cluster = clusters.get(nClust);
			PoleCluster [] two_clusters = cluster.checkAndSplit(
					norm_ds,
					min_dist,
					must_zero);
			if (two_clusters != null) {
				num_split++;
				clusters.remove(nClust);
				// do not add "headless" clusters, at least one of them will remain
				if (two_clusters[1].getNumTiles() > 0) {
					clusters.add(nClust,two_clusters[1]);
				} else {
					System.out.println("Cluster 1 "+two_clusters[1]+ " has no tiles - discarding");
				}
				if (two_clusters[0].getNumTiles() > 0) {
					clusters.add(nClust,two_clusters[0]);
				} else {
					System.out.println("Cluster 0 "+two_clusters[0]+ " has no tiles - discarding");
				}
				System.out.println("Number of clusters: "+clusters.size());
			}


		}


		return num_split;
	}




	double [][] dbgClusterLayers( // layer and eBox should be set
			boolean show_bbox,
			boolean show_ebox,
			boolean show_lines,
			boolean show_masks,
			boolean selected_only,
			boolean show_rdisparity,
			boolean show_strength,
			boolean filter,
			ArrayList<PoleCluster> clusters)
	{
		final boolean headlessOK = false;
		final TileNeibs  tnImage = biCamDSI.tnImage;
		int tilesX = tnImage.sizeX;
		int tilesY = tnImage.sizeY;
		int num_tiles =  tilesY*tilesX;
		int num_layers = 0;
		for (PoleCluster cluster: clusters){
			int cl = cluster.getLayer();
			if (cl > num_layers) num_layers = cl;
		}
		num_layers++;
		double [][] dbg_layers = new double [num_layers][num_tiles];
		for (int i = 0; i < num_layers; i++) {
			for (int j=0; j < num_tiles; j++) {
				dbg_layers[i][j] = Double.NaN;
			}
		}
		for (PoleCluster cluster: clusters) if (headlessOK ||(cluster.getNumTiles() > 0)){ // skip headless clusters
			int layer = cluster.getLayer();
			if (filter && (cluster.poleFilter() <0)) continue;
			double disp = cluster.getTargetDisparity(); // getMeanDisparity();
			Rectangle [] rectangles = {show_bbox?cluster.getBBox():null,show_ebox?cluster.getEBox():null};
			for (int nbox = 0; nbox< rectangles.length; nbox++) if (rectangles[nbox] != null) {
				Rectangle box = new Rectangle(rectangles[nbox]);
				int nTile = box.y * tilesX + box.x;
				for (int dx = 0; dx < box.width;dx++) {
					if (tnImage.getNeibIndex(nTile, dx, 0) < 0 ) {
						System.out.println("dbgClusterLayers() bug 1: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dx="+dx);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, 0)] = disp; // should never be out of bounds
					if (tnImage.getNeibIndex(nTile, dx, box.height -1) < 0 ) {
						System.out.println("dbgClusterLayers() bug 2: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dx="+dx);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, box.height -1)] = disp; // should never be out of bounds
				}
				for (int dy = 1; dy < box.height - 1; dy++) {
					if (tnImage.getNeibIndex(nTile, 0, dy) < 0 ) {
						System.out.println("dbgClusterLayers() bug 3: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dy="+dy);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, 0,            dy)] = disp; // should never be out of bounds
					if (tnImage.getNeibIndex(nTile, box.width -1, dy) < 0 ) {
						System.out.println("dbgClusterLayers() bug 4: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dy="+dy);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, box.width -1, dy)] = disp; // should never be out of bounds
				}
			}
		}

		if (show_lines) {
			int ends = 1;
			for (PoleCluster cluster: clusters) if (headlessOK ||(cluster.getNumTiles() > 0)){
				if (filter && (cluster.poleFilter() <0)) continue;
				double [] pole_line = cluster.getPoleLine();
				if (pole_line != null) {
					int layer = cluster.getLayer();
					double disp = cluster.getTargetDisparity(); // getMeanDisp();
					Rectangle box = cluster.getEBox();
					int nTile = box.y * tilesX + box.x;
					double rx0 = 0.5*box.width;
					double ry0 = 0.5*box.height;

					for (int dy = -ends; dy < box.height + ends; dy++) {
						int dx = (int) Math.round(rx0 + pole_line[0] +  (dy - ry0)*pole_line[1]);
						int nTile1 = tnImage.getNeibIndex(nTile, dx, dy);
						if (nTile1 >= 0) {
							dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, dy)] = disp;
						}
					}
				}
			}
		}
		if (show_masks) {
			for (PoleCluster cluster: clusters) if (headlessOK ||(cluster.getNumTiles() > 0)){
				if (filter && (cluster.poleFilter() <0)) continue;
				double [] pole_line = cluster.getPoleLine();
				if (pole_line != null) {
					int layer = cluster.getLayer();
					double disp = cluster.getTargetDisparity(); // getMeanDisp();
					Rectangle box = cluster.getEBox();
					int nTile = box.y * tilesX + box.x;
					for (int dy = 0; dy < box.height; dy++) {
						for (int dx = 0; dx < box.width; dx++) {
							if (cluster.getMask (dx, dy)) {
								int nTile1 = tnImage.getNeibIndex(nTile, dx, dy);
								if (nTile1 >= 0) {
									dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, dy)] = disp;
								}
							}
						}
					}
				}
			}
		}

		double [][] dbg_layers_2 = dbg_layers;


		if (show_rdisparity && show_strength) {
			dbg_layers_2= new double [2*num_layers][];
			for (int i = 0; i < num_layers; i++) {
				dbg_layers_2[i] =            dbg_layers[i];
				dbg_layers_2[i+num_layers] = dbg_layers[i].clone();
			}
		}
		if (show_rdisparity || show_strength) {
			int offs_strength =  show_rdisparity ? num_layers:0;
			for (PoleCluster cluster: clusters) if (headlessOK ||(cluster.getNumTiles() > 0)){
				if (filter && (cluster.poleFilter() <0)) continue;
				int layer = cluster.getLayer();
				double [][] ds = cluster.getDsPairs(selected_only);
				int [] nTiles =  cluster.getDsTiles(selected_only);
				for (int i = 0; i < nTiles.length; i++) {
					int nTile = nTiles[i];
					if (show_rdisparity) {
						dbg_layers_2[layer][nTile] =  ds[0][i];
					}
					if (show_strength) {
						dbg_layers_2[layer+offs_strength][nTile] =  ds[1][i];
					}
				}
			}
		}
		return dbg_layers_2;
	}


	double [][] getClusterLayers( // layer and eBox should be set
			ArrayList<PoleCluster> clusters)
	{
		final TileNeibs  tnImage = biCamDSI.tnImage;
		int tilesX = tnImage.sizeX;
		int tilesY = tnImage.sizeY;
		int num_tiles =  tilesY*tilesX;
		int num_layers = 0;
		for (PoleCluster cluster: clusters){
			int cl = cluster.getLayer();
			if (cl > num_layers) num_layers = cl;
		}
		num_layers++;
		double [][] dbg_layers = new double [num_layers][num_tiles];
		for (int i = 0; i < num_layers; i++) {
			for (int j=0; j < num_tiles; j++) {
				dbg_layers[i][j] = Double.NaN;
			}
		}
		for (PoleCluster cluster: clusters){
			int layer = cluster.getLayer();
			double disp = cluster.getTargetDisparity(); // getMeanDisp();
			Rectangle box = new Rectangle(cluster.getEBox());
			int nTile = box.y * tilesX + box.x;
			for (int dy = 0; dy < box.height; dy++) {
				for (int dx = 0; dx < box.width;dx++) {
					dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, dy)] = disp; // should never be out of bounds
				}
			}
		}
		return dbg_layers;
	}


	public void addMeasuredTiles(
			final boolean       reset_data,
			final double [][][] measured_layers,
			final ArrayList<PoleCluster> clusters
			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
//		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final int tilesX = tnImage.sizeX;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						Rectangle box = cluster.getEBox();
						int layer = cluster.getLayer();
						if (reset_data) {
							cluster.resetDsPairs();
						}
						int nTile = box.x + tilesX * box.y;
						for (int dy = 0; dy < box.height; dy++) {
							for (int dx = 0; dx < box.width; dx++) {
								int nTile1 = tnImage.getNeibIndex(nTile, dx, dy); // should always be inside image
								if (!Double.isNaN(measured_layers[layer][0][nTile1])) {
									cluster.addDs(nTile1, measured_layers[layer][0][nTile1], measured_layers[layer][1][nTile1]);
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}


	// returns number of clusters that need to be split
	public int calcPoleLines(
			final double max_diff,
			final double max_tilt,
			final double max_rms,
			final int    min_tiles,
			final double damp_tilt,
			final ArrayList<PoleCluster> clusters,
			final int debugLevel)
	{
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						cluster.calcPoleLine(
								max_diff,    // double max_diff,
								max_tilt,    // double max_tilt,
								max_rms,     // double max_rms,
								min_tiles,   // int    min_tiles,
								damp_tilt,   // double damp_tilt,
								debugLevel); // int debugLevel);
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return 0;
	}


	public void createPoleMasks(
			final double [][] norm_ds,
			final int         min_neibs,
			final boolean     use_seed,
			final double      width,
			final double      disp_aover,
			final double      disp_rover,
			final ArrayList<PoleCluster> clusters,
			final int debugLevel)
	{
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						cluster.copyNormDS(norm_ds);
						cluster.createPoleMask(
//								norm_ds,     // double [][] norm_ds,
								min_neibs,   // int         min_neibs,
								use_seed,    // boolean     use_seed,
								width,       // double      width,
								disp_aover,  // double      disp_aover,
								disp_rover,  // double      disp_rover,
								debugLevel); // int         debugLevel)
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}


	public void filterByMask(
			final boolean trim_bottoms,
			final ArrayList<PoleCluster> clusters)
	{
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						cluster.filterByMask(trim_bottoms);
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}

	public void filterByRightLeftSeparation(
			final double min_strength,
			final double disp_adiff,
			final double disp_rdiff,
			final double [][] norn_ds,
			final ArrayList<PoleCluster> clusters)
	{
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						boolean sep = cluster.checkRightLeft(
								min_strength, // double min_strength,
								disp_adiff, // double disp_adiff,
								disp_rdiff, // double disp_rdiff,
								norn_ds); //  [][] norn_ds);
						if (!sep) {
							System.out.println(String.format("Disabling cluster {%d, %d, %d, %d} disp= %f as it is not separated on right or left",
									cluster.eBox.x, cluster.eBox.y, cluster.eBox.width, cluster.eBox.height, cluster.target_disparity));
							cluster.disabled = true;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}

	public double applyMeasuredDisparity(
			final double disparity_scale,  // target disparity to differential disparity scale (baseline ratio)
			final double diff_power,       // bias towards higher disparities - (disparity+offset) is raised to this power and applied to weight
			                         // if 0.0 - do not apply value to weight
			final double diff_offset,      // add to measured differential disaprity before raising to specified power
			final int    cut_bottom,       // cut few tile rows from the very bottom - they may be influenced by ground objects
			final double keep_bottom,      // do not cut more that this fraction of the bounding bow height
			final ArrayList<PoleCluster> clusters,
			final int debugLevel)          // debug level
	{
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						double diff = cluster.applyMeasuredDisparity(
								disparity_scale,  // double disparity_scale,  // target disparity to differential disparity scale (baseline ratio)
								diff_power,       // double diff_power,       // bias towards higher disparities - (disparity+offset) is raised to this power and applied to weight
								//                                               if 0.0 - do not apply value to weight
								diff_offset,      // double diff_offset,      // add to measured differential disaprity before raising to specified power
								cut_bottom,       // int    cut_bottom,       // cut few tile rows from the very bottom - they may be influenced by ground objects
								keep_bottom,      // double keep_bottom,      // do not cut more that this fraction of the bounding bow height
								debugLevel);      // int    debugLevel);      // debug level
						if (debugLevel > -2) {
							System.out.println("applyMeasuredDisparity() for cluster "+cluster.toString()+" -> target disparity change: "+diff);
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		double max_diff = 0.0;
		for (PoleCluster cluster: clusters) {
			if (Math.abs(cluster.getTargetDiff()) > max_diff) {
				max_diff = Math.abs(cluster.getTargetDiff());
			}
		}
		return max_diff;
	}

	public void printClusterStats(
			final ArrayList<PoleCluster> clusters)
	{
		int [][] real_poles = {
				{232,109},
				{233,109},
				{226,109},
				{219,108},
				{208,106},
				{191,101},
				{ 99,112}, //{101,112},
				{ 35,112}, //{37,112},
				{123,112}, //{124,112},
				{ 84,110},
				{128,87}};


		for (int nClust = 0; nClust < clusters.size(); nClust++) {
			PoleCluster cluster = clusters.get(nClust);
			if (cluster.poleFilter() >= 0) {
				boolean real_pole = false;
				for (int i = 0; i < real_poles.length; i++) {
					if ((real_poles[i][0] == cluster.eBox.x) && (real_poles[i][1] == cluster.eBox.y)) {
						real_pole = true;
						break;
					}
				}
				Rectangle box = cluster.eBox;
				boolean unselected_only = true;
				double bg_strength = cluster.getAverageBgStrength(unselected_only);
				double meas_strength = cluster.getAverageStrength(true); // divided by measured tiles only
				double mask_strength = cluster.getAverageStrength(false); // divided by a total number of pole tiles (by mask)
				double rstrength = mask_strength/bg_strength;


				System.out.println(String.format("%1s  %1s: %3d {%3d,%3d,%3d,%3d} %6.4f %3d (%3d) frac=%6.4f str=%7.5f str_mask=%7.5f frac_height=%5.3f fract_box =%5.3f gap=%2d x0=%6.3f/%6.3f rms=%8.6f, rstrength=%7.4f",
						((cluster.poleFilter() > 0)?"*":" "),
						(real_pole?"#":"-"),
						nClust, box.x, box.y, box.width, box.height, cluster.getTargetDisparity(), cluster.getSelectionSize(), cluster.getMaskSize(),cluster.getFractSelected(),
						meas_strength, mask_strength,cluster.getFractHeight() , cluster.getFractBoxHeight(), cluster.getFreeStemGap(),
						cluster.xc_center, cluster.slope, cluster.getMDisparityRMS(),
						rstrength
						));
			}
		}
	}








	public void sortClustersByDisparity(
			ArrayList<PoleCluster> clusters)
	{
		Collections.sort(clusters, new Comparator<PoleCluster>() {
			@Override
			public int compare(PoleCluster lhs, PoleCluster rhs) {
				// -1 - less than, 1 - greater than, 0 - equal
				return (rhs.target_disparity > lhs.target_disparity) ? -1 : (rhs.target_disparity < lhs.target_disparity) ? 1 : 0;
			}
		});

	}

}
