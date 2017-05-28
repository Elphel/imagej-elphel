/**
 **
 ** TilePlanes - detect planes in tile clusters
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  TilePlanes.java is free software: you can redistribute it and/or modify
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
import Jama.EigenvalueDecomposition;
import Jama.Matrix;


public class TilePlanes {
	private int tileSize = 0; // 8;
	private int stSize =   0; // 8;
	GeometryCorrection   geometryCorrection = null;
	public TilePlanes(
			int tileSize,
			int stSize,
			GeometryCorrection geometryCorrection){
		this.tileSize = tileSize;
		this.stSize = stSize;
		this.geometryCorrection =geometryCorrection;
	}
	public TilePlanes(
			int tileSize,
			int stSize){
		this.tileSize = tileSize;
		this.geometryCorrection = null;
	}
	public class PlaneData{
		GeometryCorrection   geometryCorrection = null;
		boolean     correctDistortions = false;
		// just for visualization - no there can be several measured layers and same tile can be used multiple times. Will just logical or 
		boolean []  plane_sel =  null; // tile selection - has twice supertile size in each direction
		double []   zxy =        null; // [3] - plane center point {disparity, x, y), x=0, y=0 is a 4,4 point of an 8x8 supertile (in pixels, relative to this supertile center)
		double [][] vectors =    null; // [3][3] - re-ordered/re-directed eigenvectors(transposed): [0] - plane normal, most Z-like, towards camera, [1] - X-like, [2] - Y-like
		double []   values =     null; // [3] -eigenvalues
		int         num_points = 0;
		double      weight =     0.0;
		double []   center_xyz = null; // center of this this "plane" (ellipsoid) center in world coordinates
		double []   world_xyz =  null; // world coordinates of the nearest point of the plane, in meters
		double []   world_v1 =   null; // world in-plane vector, corresponding to vectors[1] 
		double []   world_v2 =   null; // world in-plane vector, corresponding to vectors[1]
//		double []   daxy      =  null; // disparity and 2 relative angles (ax and ay) corresponding to fisheye view, near (0,0) scale is pixel size
		// for now keeping both weighted and equal weight merged value - later remove less useful
		double [][] merged_eig_val = null; // for each of the directions (N, NE, .. NW) quality match for each layer 
		double [][] merged_eig_eq =  null; // for each of the directions (N, NE, .. NW) quality match for each layer - ignoring weights 
		boolean [][] merged_valid = null; // for each of the directions (N, NE, .. NW) if it is possible to connect with link swaps
		
		int    []   neib_best =  null; // new int [8]; // for each of the directions (N, NE, .. NW) index of best match, -1 if none 
// stores "worsening" of merging 2 	planes. if L1,L2,L = values[0] of plane1, plane2 plane composite: w1, w2 - weights for plane1, plane2
//		Lav = Math.sqrt((L1 * L1 * w1 + L2 * L2 * w2)/(w1 + w2))
// worsening_12 = (L - Lav) * (w1 + w2) * (w1 + w2) / (Lav * x1 * w2)		
		
		int         tileSize;
		int         superTileSize;
		int []      sTileXY =    null; // X and Y indices of this superTile in the image
		
		MeasuredLayers measuredLayers     =  null;
		boolean [][] measuredSelection =     null; // [number of layers in measuredLayers][2*superTileSize * 2*superTileSize]
		boolean   [] sel_mask          =     null; // selection mask - may be used for splitting plane along a line - each half can have mask
		double       measured_strength_pow = 1.0;
		double       strength_floor =        0.0;
		double       min_weight =            0.0;  // minimal weight of the ellipsoid
		int          min_tiles =             10;
		double       dispNorm =              5.0;  //  Normalize disparities to the average if above
		
		boolean      smplMode = true;   // Use sample mode (false - regular tile mode)
		int          smplSide = 2;      // Sample size (side of a square)
		int          smplNum  = 3;      // Number after removing worst
		double       smplRms  = 0.1;    // Maximal RMS of the remaining tiles in a sample
		
		double [] starValueWeight = null; 
		double    conn_density = Double.NaN; // 
		
		boolean      preferDisparity = false;
		
		public PlaneData clone(){
			PlaneData pd = new PlaneData(
					this.sTileXY,
					this.tileSize,
					this.superTileSize,
					this.geometryCorrection,
					this.correctDistortions);
//			pd.correctDistortions = this.correctDistortions;
			pd.num_points = this.num_points; 
			pd.weight =     this.weight; 
			if (this.plane_sel != null)  pd.plane_sel =  this.plane_sel.clone();
			if (this.values != null)     pd.values =     this.values.clone();
			
			if (this.zxy != null)        pd.zxy =        this.zxy.clone();
			// World calculations should be invalidated during cloning?
/*			
			if (this.center_xyz != null) pd.center_xyz = this.center_xyz.clone(); 
			if (this.world_xyz != null)  pd.world_xyz =  this.world_xyz.clone(); 
			if (this.world_v1 != null)   pd.world_v1 =   this.world_v1.clone(); 
			if (this.world_v2 != null)   pd.world_v2 =   this.world_v2.clone();
*/
			if (this.vectors != null) {
				pd.vectors = new double[3][];
				pd.vectors[0] = this.vectors[0].clone();
				pd.vectors[1] = this.vectors[1].clone();
				pd.vectors[2] = this.vectors[2].clone();
			}
			if (this.measuredLayers != null) pd.measuredLayers = this.measuredLayers;
			
			pd.setMeasSelection(this.measuredSelection);
			
			if (this.sel_mask != null) pd.sel_mask = this.sel_mask.clone();
			
			pd.measured_strength_pow = this.measured_strength_pow;
			pd.strength_floor =        this.strength_floor;

			pd.min_weight =            this.min_weight;
			pd.min_tiles =             this.min_tiles;
			pd.dispNorm =              this.dispNorm;
			
			pd.smplMode =              this.smplMode;
			pd.smplSide =              this.smplSide;
			pd.smplNum =               this.smplNum;
			pd.smplRms =               this.smplRms;
			
			pd.preferDisparity =       this.preferDisparity;
			
			copyNeib(this,pd);
			
			if (starValueWeight != null){
				pd.starValueWeight = starValueWeight.clone();
			}
			pd.conn_density =      this.conn_density;
			return pd;
		}
//		public void setConnectionDensity(double density){
//			conn_density = density;
//		}
		
		public double getConnectionDensity(){
			return conn_density;
		}
		
//		public void setStarValueWeight(double value, double weight){
//			this.starValueWeight = new double[2];
//			this.starValueWeight[0] = value;
//			this.starValueWeight[1] = weight;
//			System.out.println("setStarValueWeight(): conn_density is not set");
//		}

		public void setStarValueWeight(double[] val_weight){
			this.starValueWeight = new double[2];
			this.starValueWeight[0] = val_weight[0];
			this.starValueWeight[1] = val_weight[1];
			this.conn_density = 0.0;
//			if (val_weight.length > 2){
			this.conn_density = val_weight[2];
//			}
		}

		public double [] getStarValueWeight()
		{
			return starValueWeight;
		}
		
		
		public void setSelMask (boolean []sel_mask)
		{
			this.sel_mask = sel_mask;
		}
		public boolean [] getSelMask ()
		{
			return this.sel_mask;
		}
		
		public void setMeasSelection(boolean [][] meas_sel)
		{
			if (meas_sel == null) 
				this.measuredSelection = null;
			else {
				this.measuredSelection = meas_sel.clone();
				for (int i = 0; i < meas_sel.length; i++){
					if (meas_sel[i] != null) {
						this.measuredSelection[i] = meas_sel[i].clone();
					}
				}
			}
		}

		public void orMeasSelection(boolean [][] meas_sel)
		{
			if (meas_sel == null) 
				this.measuredSelection = null;
			else {
				if (this.measuredSelection == null) {
					this.measuredSelection = meas_sel.clone();
				}
				for (int i = 0; i < meas_sel.length; i++){
					if (meas_sel[i] != null) {
						if (this.measuredSelection[i] == null) {
							this.measuredSelection[i] = meas_sel[i].clone();
						} else {
							for (int j = 0; j < meas_sel[i].length; j++){
								this.measuredSelection[i][j] |= meas_sel[i][j];
							}
						}
					} 
				}
			}
		}
		
		
		
		public boolean [] getMeasSelection(int nl){
			if (this.measuredSelection == null) {
				return null;
			}
			return 	this.measuredSelection[nl];
		}

		public boolean [][] getMeasSelection(){
			return 	this.measuredSelection;
		}
		
		
		public MeasuredLayers getMeasuredLayers()
		{
			return this.measuredLayers;
		}
		
		public void copyNeib(
				PlaneData src,
				PlaneData dst)
		{
			if (src.merged_eig_val != null){
				dst.merged_eig_val = src.merged_eig_val.clone();
				for (int i = 0; i < src.merged_eig_val.length; i++){
					if (src.merged_eig_val[i] != null){
						dst.merged_eig_val[i] = src.merged_eig_val[i].clone();
					}
				}
			}
			
			if (src.merged_eig_eq != null){
				dst.merged_eig_eq = src.merged_eig_eq.clone();
				for (int i = 0; i < src.merged_eig_eq.length; i++){
					if (src.merged_eig_eq[i] != null){
						dst.merged_eig_eq[i] = src.merged_eig_eq[i].clone();
					}
				}
			}
			
			
			if (src.merged_valid != null){
				dst.merged_valid = src.merged_valid.clone();
				for (int i = 0; i < src.merged_valid.length; i++){
					if (src.merged_valid[i] != null){
						dst.merged_valid[i] = src.merged_valid[i].clone();
					}
				}
			}
			
			if (src.neib_best != null) dst.neib_best = src.neib_best.clone();
			
			// also copy original plane parameters - tile selection and number of points
			
			dst.num_points = src.num_points; 
			if (src.plane_sel != null)  dst.plane_sel =  src.plane_sel.clone(); 
		}
		
		public void invalidateCalculated()
		{
			this.center_xyz = null; // center of this supertile this plane center in world coordinates
			this.world_xyz =  null; // world coordinates of the nearest point of the plane, in meters
			this.world_v1 =   null; // world in-plane vector, corresponding to vectors[1] 
			this.world_v2 =   null; // world in-plane vector, corresponding to vectors[1]
		}
		
		
		public PlaneData (
				int [] sTileXY, 
				int tileSize,
				int superTileSize,
				GeometryCorrection   geometryCorrection,
				boolean              correctDistortions)
		{
			this.geometryCorrection = geometryCorrection;
			this.correctDistortions = correctDistortions;
			this.tileSize = tileSize;
			this.superTileSize = superTileSize;
			this.sTileXY = sTileXY.clone();
		}

		public PlaneData (
				int [] sTileXY, 
				int tileSize,
				GeometryCorrection   geometryCorrection,
				boolean              correctDistortions,
				MeasuredLayers measuredLayers,
				boolean preferDisparity)
		{
			this.geometryCorrection = geometryCorrection;
			this.correctDistortions = correctDistortions;
			this.tileSize = tileSize;
			this.superTileSize = measuredLayers.getSuperTileSize();
			this.sTileXY = sTileXY.clone();
			this.measuredLayers = measuredLayers;
			this.preferDisparity = preferDisparity;
		}

		/**
		 * Create separation masks for two crossing planes. Uses neighbors connections to distinguish between concave and convex 
		 * @param pd1 first plane data instance to create sel_mak
		 * @param pd2 second plane data instance to create sel_mak
		 * @param tolerance disparity tolerance for separation (will be normalized by dispNorm for large disparities)
		 * @param min_tiles minimal number of tiles in each half
		 * @param debugLevel debug level
		 * @return true if OK, false if the planes are not crossing (masks not created)
		 */
		public boolean calcSelMasks(
				PlaneData pd1,
				PlaneData pd2,
				double    tolerance,
				int       min_tiles,
				int       debugLevel
				)
		{
			int st2 = 2 * superTileSize;
			int [] dirs = {-st2, -st2 + 1, 1, st2 + 1, st2, st2 - 1, -1, -st2 - 1};
			double [][] planes = {
					pd1.getDoublePlaneDisparity(false),
					pd2.getDoublePlaneDisparity(false)};
			PlaneData [] pair = {pd1,pd2};
			int len2 = planes[0].length;
			boolean [][] sel_masks = new boolean [2][len2];
			double concave = 0; // positive - concave, negative - convex
			int cent_indx = (st2/2 -1) * (st2+1);
			for (int np = 0; np < 2; np++){
				int [] neibs = pair[np].getNeibBest();
				for (int dir = 0; dir < dirs.length; dir++) {
					if (neibs[dir] >= 0) { // neighbor connected
						int indx = cent_indx + dirs[dir];
						concave += (planes[np][indx] - planes[np][cent_indx]) - (planes[1 - np][indx] - planes[1 - np][cent_indx]);
					}
				}
			}
			// swap plane selections for convex plane intersection (edge is the closest)
			int first = (concave < 0) ? 1 : 0;
			int [] nums = {0, 0};
			
			for (int i = 0; i <len2; i++){
				double d_av = 0.5 * (planes[0][i] + planes[1][i]);
				double diff = planes[0][i] - planes[1][i];
				if ((dispNorm > 0.0) && (d_av > dispNorm)){
					diff *= dispNorm/d_av;
				}
				if (diff > -tolerance) {
					sel_masks[first][i] = true;
					nums[first] ++;
				}
				if (-diff > -tolerance) {
					sel_masks[1-first][i] = true;
					nums[1-first] ++;
				}
			}
			if ((nums[0] < min_tiles) || (nums[1] < min_tiles)){
				return false;
			}
			// apply selections
			pd1.setSelMask(sel_masks[0]);
			pd2.setSelMask(sel_masks[1]);
			return true;
		}
		
		/**
		 * Spit tiles belonging to this between multiple PlaneData instances
		 * @param pd_set array of plane data instances
		 * @param single_plane it is a single plane to split - use all tiles, not just previously
		 *  selected
		 * @param max_diff maximal normalized disparity difference from the plane to consider
		 * @param other_diff maximal difference of the added tile ratio to the average
		 *  disparity difference of the exclusively selected tiles
		 * 
		 * @param non_exclusive allow the same tile data to belong to multiple PD instances
		 * @param use_other_planes allow the same tile not included in this PD to be used
		 * @param smplMode use square sample mode, false - single-tile samples 
		 * @param smplSide size of the square sample side 
		 * @param smplNum number of averaged samples (should be <= smplSide * smplSide and > 1) 
		 * @param smplRms maximal square root of variance (in disparity pixels) to accept the result
		 * @param measSel (with use_other_planes) select measurements for supertiles :
		 *  +1 - combo, +2 - quad +4 - hor +8 - vert
		 * @param allow_parallel allow parallel shift of each plane before adding more data
		 * @param debugLevel debug level
		 * @return true if OK
		 */
		
		public boolean splitPlaneTiles (
				PlaneData [] pd_set,
				boolean      single_plane,
				double       max_diff, // maximal disparity difference (0 - any), will be normalized by dispNorm
				double       other_diff, 
				boolean      non_exclusive,
				boolean      use_other_planes,
				boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				int          smplSide, //        = 2;      // Sample size (side of a square)
				int          smplNum, //         = 3;      // Number after removing worst
				double       smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				int          measSel, // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				boolean      allow_parallel,
				int          debugLevel)
		{
			if (debugLevel>0){
				System.out.println("splitPlaneTiles");
			}
			
			double [][] planes = new double [pd_set.length][];
			boolean [][][] tsel = new boolean [pd_set.length][measuredLayers.getNumLayers()][];
			boolean [][] sel_masks = new boolean [pd_set.length][];
			for (int np = 0; np < pd_set.length; np ++) if (pd_set[np] != null){
				planes[np] = pd_set[np].getDoublePlaneDisparity(false);
				sel_masks[np] = pd_set[np].getSelMask();
				for (int nl = 0; nl < measuredLayers.getNumLayers(); nl++){
					if (measuredSelection[nl] != null){
						tsel[np][nl] = new boolean [measuredSelection[nl].length];
					}
				}
//				pd_set[np].setMeasSelection(measuredSelection); // copy same selection to all planes
			}
			double max_diff2 = max_diff * max_diff;
			double other_diff2 =other_diff * other_diff; 
			double [] sd2 = new double[pd_set.length];
			double [] sd2_av = new double[pd_set.length];
			double [] sw =  new double[pd_set.length];
			double [][][] disp_strength = new double[measuredLayers.getNumLayers()][][];
			// split exclusively, calculate rms for each, then add others if RMS is not increased
			for (int nl = 0; nl < measuredLayers.getNumLayers(); nl ++){
				if ((measuredSelection[nl] != null) &&  ((measSel & (1 << nl)) !=0)) {
					if (smplMode) {
						disp_strength[nl] = measuredLayers.getDisparityStrength( // expensive to calculate (improve removing outlayers
								nl, // int num_layer,
								getSTileXY()[0],        // int stX,
								getSTileXY()[1],        // int stY,
								(single_plane ? null : measuredSelection[nl]),  // boolean [] sel_in,
								strength_floor,
								measured_strength_pow, //
								smplSide, //        = 2;      // Sample size (side of a square)
								smplNum, //         = 3;      // Number after removing worst
								smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
								true); // boolean null_if_none)
					} else {
						disp_strength[nl] = measuredLayers.getDisparityStrength(
								nl,                     // int num_layer,
								getSTileXY()[0],        // int stX,
								getSTileXY()[1],        // int stY,
								(single_plane ? null : measuredSelection[nl]),  // boolean [] sel_in,
								strength_floor,         //  double strength_floor,
								measured_strength_pow,  // double strength_pow,
								true);                  // boolean null_if_none);
					}
					if (disp_strength[nl] != null) {
						for (int indx = 0; indx < disp_strength[nl][1].length; indx++){
							double w = disp_strength[nl][1][indx];
							if (w > 0.0){
								double d = disp_strength[nl][0][indx];
								double d2_best = Double.NaN;
								int np_best = 0;
								for (int np = 0; np < planes.length; np++) if (planes[np] != null){
									if ((sel_masks[np] == null) || sel_masks[np][indx]) {
										double d_av = 0.5 * (d + planes[np][indx]);
										double d2 = (d - planes[np][indx]);
										if ((dispNorm > 0.0) && (d_av > dispNorm)){
											d2 *= dispNorm/d_av;
										}
										d2 *= d2;
										if (!(d2 >= d2_best)){ // so d2_best NaN is OK
											d2_best = d2;
											np_best = np;
										}
									}
								}
								if ((max_diff2 == 0.0) || (d2_best < max_diff2)) {
									tsel[np_best][nl][indx] = true;
									sd2[np_best] += w * d2_best;
									sw[np_best] += w;
								}
							}
						}
					}
				}
			}
			for (int np = 0; np < planes.length; np++) if (planes[np] != null){
				if (sw[np] > 0.0){
					sd2_av[np] = sd2[np] / sw[np];
				}
			}
			double [] sd = new double[pd_set.length]; // will be all 0.0;

			if (allow_parallel && (non_exclusive || use_other_planes)){
				for (int nl = 0; nl < measuredLayers.getNumLayers(); nl ++){
					if (disp_strength[nl] != null) {
						for (int indx = 0; indx < disp_strength[nl][1].length; indx++){
							double w = disp_strength[nl][1][indx];
							if (w > 0.0){
								double d = disp_strength[nl][0][indx];
								for (int np = 0; np < planes.length; np++) if (planes[np] != null){
									if (!tsel[np][nl][indx]) { // not already in that plane
										sd[np] += w * d;
									}
								}
							}
						}
					}
				}
				for (int np = 0; np < planes.length; np++) if (planes[np] != null){
					if (sw[np] > 0.0){
						sd[np] /= sw[np];
						sd2_av[np] -= sd[np] * sd[np];
						sd2[np] = sd2_av[np] * sw[np]; // to add more 
					}
				}
				
			}
			
			
			if (non_exclusive) {
				for (int nl = 0; nl < measuredLayers.getNumLayers(); nl ++){
					if ((measuredSelection[nl] != null) &&  ((measSel & (1 << nl)) !=0)) {
						if (disp_strength[nl] != null) {
							for (int indx = 0; indx < disp_strength[nl][1].length; indx++){
								double w = disp_strength[nl][1][indx];
								if (w > 0.0){
									double d = disp_strength[nl][0][indx];
									for (int np = 0; np < planes.length; np++) if (planes[np] != null){
										if ((sel_masks[np] == null) || sel_masks[np][indx]) {
											if (!tsel[np][nl][indx]) { // not already in that plane
												double d_av = 0.5 * (d - planes[np][indx]);
												double d2 = (d - planes[np][indx]);
												d2 -= sd[np]; // subtract parallel shift
												if ((dispNorm > 0.0) && (d_av > dispNorm)){
													d2 *= dispNorm/d_av;
												}
												d2 *= d2;
												if (d2 <= sd2_av[np] * other_diff2) { // not more than exclusive tile variance 
													tsel[np][nl][indx] = true;
													sd2[np] += w * d2;
													sw[np] += w;
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
			
			if (use_other_planes && !single_plane) { // no need if single_plane - it already got all planes it could
				for (int nl = 0; nl < measuredLayers.getNumLayers(); nl ++){
					if ((measuredSelection[nl] != null) &&  ((measSel & (1 << nl)) !=0)) {
						// recalculate for all measure tiles, not just selected in the original PD
						if (smplMode) {
							disp_strength[nl] = measuredLayers.getDisparityStrength( // expensive to calculate (improve removing outlayers
									nl, // int num_layer,
									getSTileXY()[0],        // int stX,
									getSTileXY()[1],        // int stY,
									null,  // boolean [] sel_in,
									strength_floor,
									measured_strength_pow, //
									smplSide, //        = 2;      // Sample size (side of a square)
									smplNum, //         = 3;      // Number after removing worst
									smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
									true); // boolean null_if_none)
						} else {
							disp_strength[nl] = measuredLayers.getDisparityStrength(
									nl,                     // int num_layer,
									getSTileXY()[0],        // int stX,
									getSTileXY()[1],        // int stY,
									null,  // boolean [] sel_in,
									strength_floor,         //  double strength_floor,
									measured_strength_pow,  // double strength_pow,
									true);                  // boolean null_if_none);
						}
						//disp_strength[nl] = measuredLayers.getDisparityStrength(
						for (int indx = 0; indx < disp_strength[nl][1].length; indx++){
							double w = disp_strength[nl][1][indx];
							if (w > 0.0){
								double d = disp_strength[nl][0][indx];
								for (int np = 0; np < planes.length; np++) if (planes[np] != null){
									if ((sel_masks[np] == null) || sel_masks[np][indx]) {
										if (!tsel[np][nl][indx]) { // not already in that plane
											double d_av = 0.5 * (d - planes[np][indx]);
											double d2 = (d - planes[np][indx]);
											d2 -= sd[np]; // subtract parallel shift
											if ((dispNorm >0.0) && (d_av > dispNorm)){
												d2 *= dispNorm/d_av;
											}
											d2 *= d2;
											if (d2 <= sd2_av[np]  * other_diff2) { // not more than exclusive tile variance 
												tsel[np][nl][indx] = true;
												sd2[np] += w * d2;
												sw[np] += w;
											}
										}
									}
								}
							}
						}					
					}
				}
			}
			
			// re-calculate variance (just for debug, not needed
			for (int np = 0; np < planes.length; np++)  if (planes[np] != null){
				if (sw[np] > 0.0){
					sd2_av[np] = sd2[np] / sw[np];
				}
			}
			// apply selections to each PD
			for (int np = 0; np < planes.length; np++)  if (planes[np] != null){
				pd_set[np].setMeasSelection(tsel[np]);
			}
			// need to re-calculate new planes, remove outliers
			return true;
		}
		
		/**
		 * Tilt disparity values around the supertile center (this.zxy) so constant disparity in the output
		 * corresponds to the real world plane parallel to the provided one. Used to discriminate tiles by
		 * the effective disparity value (disparity in the center of the supertile of the parallel plane)
		 * Adding protection from behind the horizon - areas where disparity is negative zero the
		 * result strengths
		 * @param world_normal_xyz real world 3d vector of the plane normal (0.0, 1.0, 0.0 - horizontal)
		 * @param disp_center dispariy in at the center of the supertile (to rotate around)
		 * @param tile_sel multi-layer tile selection (or null to use all available tiles)
		 * @param disp_str multi-layer disparity/strength array
		 * @param debugLevel
		 * @return same
		 */
		
		public double [][][] getDisparityToPlane(
				double []     world_normal_xyz,
				double        disp_center,
				boolean [][]  tile_sel, // null - do not use, {} use all (will be modified)
				double [][][] disp_str, // calculate just once if null 
				int           debugLevel)
		{
			final int stSize2 = 2* stSize;
			invalidateCalculated();			
			double [] px_py = getCenterPxPy();
			// find world coordinates of the center of tile intersection with the plane
			Matrix st_xyz = new Matrix(geometryCorrection.getWorldCoordinates(
						px_py[0],
						px_py[1],
						disp_center,
						this.correctDistortions),3);
			Matrix normal_row = new Matrix(world_normal_xyz,1); // 1x3
			double [][][] eff_disp_str = disp_str.clone();
			for (int ml = 0; ml < disp_str.length; ml++) {
				if (disp_str[ml] != null){
					eff_disp_str[ml] =    new double [2][];
					eff_disp_str[ml][0] = new double [disp_str[ml][0].length];
					eff_disp_str[ml][1] = disp_str[ml][1]; // keep same strengths
				}
			}
			double n_by_w =  normal_row.times(st_xyz).get(0, 0);
			if (debugLevel > 1) {
				System.out.println("st_xyz = {"+st_xyz.get(0, 0)+","+st_xyz.get(1, 0)+","+st_xyz.get(2, 0)+"}"+" ="+n_by_w);
			}
			for (int ml = 0; ml < disp_str.length; ml++) if (disp_str[ml] != null){
				for (int dy = 0; dy < stSize2; dy ++ ){
					double y = (dy -  stSize + 0.5) * tileSize;
					for (int dx = 0; dx < stSize2; dx ++ ){
						double x = (dx -  stSize + 0.5) * tileSize;
						int indx = dy * stSize2 + dx;
						if ((disp_str[ml][1][indx] > 0) && ((tile_sel == null) || ((tile_sel[ml] != null) && tile_sel[ml][indx]))){ // do not bother with zero-strength
							// Find world coordinates of a measured tile
							Matrix w_xyz = new Matrix(geometryCorrection.getWorldCoordinates(
									px_py[0] + x,
									px_py[1] + y,
									disp_str[ml][0][indx],
									this.correctDistortions),3);
							// now find intersection of the view line (0,0,0) to world_xyz with a plane through wxyz perpendicular to world_normal_xyz
							// then calculate disparity from z of that point
							// inner product of transposed
							double n_by_p = normal_row.times(w_xyz).get(0, 0);
							if (disp_str[ml][1][indx] > 0){ // do not bother with zero-strength
								double z;
								if ((n_by_p * n_by_w) > 0.0) { 
									z = st_xyz.get(2, 0)*n_by_p / n_by_w;
									// convert z to disparity
									eff_disp_str[ml][0][indx] = geometryCorrection.getDisparityFromZ (-z);
								} else {
									z = 0.0;
									eff_disp_str[ml][1][indx] = 0.0; // behind the horizon
								}
								if (debugLevel > 1) {
									System.out.println("dy = "+dy+", dx=" + dx+ " {"+w_xyz.get(0, 0)+","+w_xyz.get(1, 0)+","+w_xyz.get(2, 0)+"}"+" z="+z+" n_by_p = "+n_by_p
											+" disp = "+disp_str[ml][0][indx]+" px = "+(px_py[0] + x)+" py = "+(px_py[1] + y));
								}
							}
						}
					}
				}
			}
			return eff_disp_str;
		}
		
		/**
		 * Remove outliers from the set of tiles contributing to a single plane ellipsoid
		 * Should run after getPlaneFromMeas as some parameter4s will be copied from that run
		 * @param disp_str - pre-calculated array or null (will be calculated). disp_str
		 * has the same format as the output of getPlaneFromMeas - [measurement layer][2][tile index],
		 * so it can be used for input.
		 * @param targetEigen Target value for the lowest eigenvalue (thickness of the ellipsoid)
		 * @param maxRemoved maximal number of tiles to be removed
		 * @param debugLevel debug level
		 * @return true if OK (currently false is when a program bug happens), TODO: change to "goal reached" 
		 */
		
		public boolean removeOutliers( // getPlaneFromMeas should already have run
				double [][][] disp_str, // calculate just once when removing outliers (null - OK, will generate it) 
				double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
				int        maxRemoved,  // maximal number of tiles to remove (not a constant)
				int        debugLevel)
		{
			int stSize2 = 2 * stSize;
			if (maxRemoved > (getNumPoints() - this.min_tiles)) maxRemoved = getNumPoints() - this.min_tiles;
			boolean need_disp_str = false;
			if (disp_str == null) {
				disp_str =      new double [measuredSelection.length][][];
				for (int nl = 0; nl < measuredSelection.length; nl++){
					if (measuredSelection[nl] != null){
						if (smplMode) {
							if (need_disp_str) {
								disp_str[nl] = measuredLayers.getDisparityStrength( // expensive to calculate (improve removing outlayers
										nl, // int num_layer,
										sTileXY[0], // int stX,
										sTileXY[1], // int stY,
										null, // measuredSelection[nl], // boolean [] sel_in,
										strength_floor,
										measured_strength_pow, //
										smplSide, //        = 2;      // Sample size (side of a square)
										smplNum, //         = 3;      // Number after removing worst
										smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
										true); // boolean null_if_none)
							}
						} else {						
							disp_str[nl] = measuredLayers.getDisparityStrength(
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									null, // measuredSelection[nl], // boolean [] sel_in,
									strength_floor, // 
									measured_strength_pow, //
									true); // boolean null_if_none)
						}
					}
				}
			}
			int numRemoved = 0;
			boolean no_bugs = true;
			for (; (getValue() > targetEigen) && (numRemoved < maxRemoved); numRemoved++){
				if (debugLevel > 2){
					System.out.println("removePlaneOutliers("+sTileXY[0]+":"+sTileXY[1]+"): numRemoved = "+numRemoved+" eigenValue = " + getValue()+" target = "+targetEigen);
				}
				// make a plane and find the worst (largest disparity difference) tile
				// z = -(x*Vx + y*Vy)/Vz
				
				double worst_d2 = 0.0;
				int [] worst_layer_index = {-1,-1};
				double [] v = getVector();
				double [] zxy0 = getZxy();
				for (int nl = 0; nl < measuredSelection.length; nl++){
					if (measuredSelection[nl] != null){
						// already calculated, but not masked by selection!
						if (disp_str[nl] != null) {
							for (int indx = 0; indx < disp_str[nl][0].length; indx++){
								double w = disp_str[nl][1][indx];
								if (measuredSelection[nl][indx] && (w > 0.0)){
									double x = ((indx % stSize2) - stSize) - zxy0[1];
									double y = ((indx / stSize2) - stSize) - zxy0[2];
									double d = disp_str[nl][0][indx];
									d -= zxy0[0];
									d += (x * v[1]+y*v[2])/v[0];
									double d2 = d*d;
									if (d2 > worst_d2){
										worst_d2 = d2;
										worst_layer_index[0] = nl;
										worst_layer_index[1] = indx;
									}
								}
							}
						}
					}
				}
				if (worst_layer_index[0] < 0) {
					System.out.println("This is a BUG in removePlaneOutliers()");
					no_bugs = false;
					break;
				}
				measuredSelection[worst_layer_index[0]][worst_layer_index[1]] = false;
				if (debugLevel > 2){
					System.out.println("removePlaneOutliers() worst_layer="+worst_layer_index[0]+", worst_index="+worst_layer_index[1]);
				}
				
				boolean OK = (getPlaneFromMeas(
						measuredSelection, // null,            // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
						disp_str,						
						Double.NaN,      // double       disp_far, // minimal disparity to select (or NaN)
						Double.NaN,      // double       disp_near, // maximal disparity to select (or NaN)
						this.dispNorm,   // double       dispNorm,   //  Normalize disparities to the average if above
						this.min_weight, // double       min_weight,
						this.min_tiles,  // int          min_tiles,
						strength_floor,
						measured_strength_pow, // double       strength_pow,
						smplMode,
						smplSide,
						smplNum,
						smplRms,
						debugLevel-1) != null);
				
				if (!OK){ // restore last selection, re-run getPlaneFromMeas
					measuredSelection[worst_layer_index[0]][worst_layer_index[1]] = true;
					OK = (getPlaneFromMeas(
							measuredSelection, // null,            // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
							disp_str,						
							Double.NaN,      // double       disp_far, // minimal disparity to select (or NaN)
							Double.NaN,      // double       disp_near, // maximal disparity to select (or NaN)
							this.dispNorm,   // double       dispNorm,   //  Normalize disparities to the average if above
							this.min_weight, // double       min_weight,
							this.min_tiles,  // int          min_tiles,
							strength_floor,
							measured_strength_pow, // double       strength_pow,
							smplMode,
							smplSide,
							smplNum,
							smplRms,
							debugLevel-1) != null);
					if (!OK) {
						System.out.println("This is a BUG in removePlaneOutliers() - run with previous selection and failed");
						no_bugs = false;
					}
					break;
				}
			}			
			return no_bugs;
		}
		
		/**
		 * Get "plane" - ellipsoid from covariance matrix of measured data
		 * @param tile_sel multi-layer selection of the tiles to use (first dimension should
		 *                 match number of measurement layers. each element can be either:
		 *                 a) boolean array 4* superTileSize * superTileSize length
		 *                 b) zero-length array - it will be calculated from all measurement
		 *                    data on that layer matching optional disparity limits
		 *                 c) null - this measuremen5t layer will not be used   
		 * @param disp_str - pre-calculated array or null (will be calculated). disp_str
		 * has the same format as the output [measurement layer][2][tile index]
		 * @param disp_far optional low limit for tile disparity (NaN - do not check)
		 * @param disp_near optional high limit for tile disparity (NaN - do not check)
		 * @param dispNorm reduce scale of the disparity differences for ellipsoids with 
		 *                 disparity center above that value
		 * @param min_weight minimal total weight of the ellipsoid do process
		 * @param min_tiles minimal number of tiles used for calculation
		 * @param strength_pow raise correlation strength to this power
		 * @param smplMode use square sample mode, false - single-tile samples 
		 * @param smplSide size of the square sample side 
		 * @param smplNum number of averaged samples (should be <= smplSide * smplSide and > 1) 
		 * @param smplRms maximal square root of variance (in disparity pixels) to accept the result
		 * 
		 * @param debugLevel debug level
		 * @return per measurement layer disparity/strengths, or null if failed
		 */
		public double [][][] getPlaneFromMeas(
				boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
				double [][][] disp_str, // calculate just once when removing outlayers 
				double       disp_far, // minimal disparity to select (or NaN)
				double       disp_near, // maximal disparity to select (or NaN)
				double       dispNorm,   //  Normalize disparities to the average if above
				double       min_weight,
				int          min_tiles,
				double       strength_floor,
				double       strength_pow,

				boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				int          smplSide, //        = 2;      // Sample size (side of a square)
				int          smplNum, //         = 3;      // Number after removing worst
				double       smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
				int          debugLevel)
		{
			double mindet = 1E-15;
			int stSize2 = 2 * stSize;
			int num_tiles = 0;
			double sw = 0.0;
			if (tile_sel != null) {
				this.measuredSelection =      tile_sel;
			} else {
				tile_sel = this.measuredSelection;
			}
			
			this.strength_floor =         strength_floor;
			this.measured_strength_pow =  strength_pow;
			this.min_weight =             min_weight;
			this.min_tiles =              min_tiles;
			this.dispNorm =               dispNorm;
			this.smplMode =               smplMode; //        = true;   // Use sample mode (false - regular tile mode)
			this.smplSide =               smplSide; //        = 2;      // Sample size (side of a square)
			this.smplNum =                smplNum;    //         = 3;      // Number after removing worst
			this.smplRms =                smplRms;    //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			
			if (debugLevel > 2){
				System.out.println("getPlaneFromMeas()");
			}
			boolean need_disp_str = false;
			if (disp_str == null) {
				disp_str =      new double [tile_sel.length][][];
				need_disp_str = true;
			}
			for (int nl = 0; nl < tile_sel.length; nl++){
				if (tile_sel[nl] != null){
					if (smplMode) {
						if (need_disp_str) {
							disp_str[nl] = measuredLayers.getDisparityStrength( // expensive to calculate (improve removing outlayers
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									//tile_sel[nl], // boolean [] sel_in,
									strength_floor,
									strength_pow, //
									smplSide, //        = 2;      // Sample size (side of a square)
									smplNum, //         = 3;      // Number after removing worst
									smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
									true); // boolean null_if_none)
						}
						if (disp_str[nl] == null)	continue;
						if (Double.isNaN(disp_far) && Double.isNaN(disp_near)){
							tile_sel[nl] =  measuredLayers.getSupertileSelection(
									disp_str[nl],
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									true); // boolean null_if_none)
						} else {
							tile_sel[nl] =  measuredLayers.getSupertileSelection(
									disp_str[nl],
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									disp_far,      // 		double     disp_far,
									disp_near,     // double     disp_near,
									true);         // boolean null_if_none)
						}
						sw += MeasuredLayers.getSumStrength(disp_str[nl],tile_sel[nl]);
						num_tiles += MeasuredLayers.getNumSelected(tile_sel[nl]);
						
					} else {
						if (Double.isNaN(disp_far) && Double.isNaN(disp_near)){
							tile_sel[nl] =  measuredLayers.getSupertileSelection(
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									strength_floor,
									true); // boolean null_if_none)
						} else {
							tile_sel[nl] =  measuredLayers.getSupertileSelection(
									nl,            // int num_layer,
									sTileXY[0],    // int stX,
									sTileXY[1],    // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									disp_far,      // 		double     disp_far,
									disp_near,     // double     disp_near,
									strength_floor,
									true);         // boolean null_if_none)

						}
						num_tiles += MeasuredLayers.getNumSelected(tile_sel[nl]);
						if (tile_sel[nl] != null){
							disp_str[nl] = measuredLayers.getDisparityStrength(
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									tile_sel[nl], // boolean [] sel_in,
									strength_floor,
									strength_pow, //
									true); // boolean null_if_none)
							sw += MeasuredLayers.getSumStrength(disp_str[nl]);
						}
					}
					
					
					if ((debugLevel > 3) && (disp_str[nl] != null)){
//					if ((debugLevel > 1) && (disp_str[nl] != null)){
						  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
						  double [][] dbg_img = new double [3][];
						  dbg_img[0] = disp_str[nl][0];
						  dbg_img[1] = disp_str[nl][1];
						  dbg_img[2] = new double [stSize2*stSize2];
						  for (int i = 0; i < dbg_img[2].length; i++){
							  dbg_img[2][i] = tile_sel[nl][i]?1.0:0.0;
						  }
						  sdfa_instance.showArrays(dbg_img,  stSize2, stSize2, true, "disp_str_x"+sTileXY[0]+"_y"+sTileXY[1]+"_"+nl);
					}
				}
			}
			this.measuredSelection =      tile_sel; // it may be modified			
			
			if ((sw < min_weight) || (num_tiles < min_tiles)) {
				if (debugLevel > 1){
					System.out.println("getPlaneFromMeas():return false");
				}
				return null; // too weak plane or too few tiles selected
			}
			
			
			double [][] acovar = new double [3][3];
			double swz = 0.0, swx = 0.0, swy = 0.0;
			sw =0.0;
			for (int nl = 0; nl < tile_sel.length; nl++){
				if (disp_str[nl] != null) {
					for (int indx = 0; indx < disp_str[nl][0].length; indx++){
						if (tile_sel[nl][indx]) {
							double w = disp_str[nl][1][indx];
							if (w > 0.0){
								double d = disp_str[nl][0][indx];
								// referencing samples to centers of pixels
								double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5; // in pixels, not in tiles
								double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5;
								sw  += w;
								swz += w * d;
								swx += w * x;
								swy += w * y;
							}
						}
					}
				}
			}
			if (sw == 0.0) {
				return null; //
			}
			swz /= sw;
			swx /= sw;
			swy /= sw;

			double kz = ((dispNorm > 0.0) && (swz > dispNorm)) ? (dispNorm / swz) : 1.0; 
			
			if (debugLevel > 0){
				System.out.println("getPlaneFromMeas(): num_tiles="+num_tiles+", sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy+", kz="+kz);
			}

			// TODO: scale disparity to make same scale for 3 axes?
			
			for (int nl = 0; nl < tile_sel.length; nl++){
				if (disp_str[nl] != null) {
					for (int indx = 0; indx < disp_str[nl][0].length; indx++){
						if (tile_sel[nl][indx]) {
							double w = disp_str[nl][1][indx] / sw;
							if (w > 0.0){
								double d =  kz * (disp_str[nl][0][indx] - swz);
								double wd = w*d;
								double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5 - swx;
								double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5 - swy;
								acovar [0][0] += wd * d;
								acovar [0][1] += wd * x;
								acovar [0][2] += wd * y;
								acovar [1][1] += w * x * x;
								acovar [1][2] += w * x * y;
								acovar [2][2] += w * y * y;
							}
						}
					}
				}
			}
			acovar [1][0] = acovar [0][1]; 
			acovar [2][0] = acovar [0][2]; 
			acovar [2][1] = acovar [1][2];
			Matrix covar = new Matrix(acovar);

			EigenvalueDecomposition eig = covar.eig();
			if (Double.isNaN(eig.getV().get(0, 0))){
				System.out.println("getCovar(): Double.isNaN(eig.getV().get(0, 0))");
				debugLevel = 20;
			}
			
			if (debugLevel > 3){
//				if (debugLevel > 0){
				System.out.println("getCovar(): sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy +", covar.det() = "+covar.det());
				System.out.println("getCovar(): covarianvce matrix, number of used points:"+num_tiles);
				covar.print(10, 6); // w,d
				System.out.println("getCovar(): eigenvalues");
				eig.getD().print(10, 6); // w,d
				System.out.println("getCovar(): eigenvectors");
				eig.getV().print(10, 6); // w,d
			}
			if ((eig.getD().get(0, 0) == 0.0) || (Math.abs(covar.det()) < mindet)) {
				return null; // testing with zero eigenvalue
				// Problem with zero eigenvalue is with derivatives and coordinate conversion
			}

			double [][] eig_val =  eig.getD().getArray(); // rslt[0];
			double [][] eig_vect = eig.getV().getArray(); // rslt[1];
			// find vector most orthogonal to view // (anyway it all works with that assumption), make it first
			// TODO normalize to local linear scales
			int oindx = 0;
			if (preferDisparity) {
				for (int i = 1; i <3; i++){
					if (Math.abs(eig_vect[0][i]) > Math.abs(eig_vect[0][oindx])){
						oindx = i;
					}
				}
			} else {
				for (int i = 1; i < 3 ; i++){
					if (eig_val[i][i] < eig_val[oindx][oindx]){
						oindx = i;
					}
				}
			}
			if (eig_val[oindx][oindx] == 0.0){
				System.out.println("getPlane(): zero eigenvalue!!");
			}
			
			// select 2 other axes for increasing eigenvalues (so v is short axis, h  is the long one)
			int vindx = (oindx == 0)? 1 : 0;
			int hindx = (oindx == 0)? 2 : ((oindx == 1) ? 2 : 1);
			if (eig_val[vindx][vindx] > eig_val[hindx][hindx]){
				int tmp = vindx;
				vindx = hindx;
				hindx = tmp;
			}

			double [][] plane = {
					{eig_vect[0][oindx],eig_vect[1][oindx],eig_vect[2][oindx]},  // plane normal to camera
					{eig_vect[0][vindx],eig_vect[1][vindx],eig_vect[2][vindx]},  // "horizontal" axis // to detect skinny planes and poles
					{eig_vect[0][hindx],eig_vect[1][hindx],eig_vect[2][hindx]}}; //  "vertical"   axis  // to detect skinny planes and poles

			// Make normal be towards camera (positive disparity), next vector - positive in X direction (right)
			for (int v = 0; v < 2; v++) {
				if (plane[v][v] < 0.0) for (int i = 0; i < 3; i ++) plane[v][i] = -plane[v][i];
			}
			
			// make  direction last vector so px (x) py (.) disp < 0 (left-hand coordinate system) 
			if (new Matrix(plane).det() > 0){
				for (int i = 0; i < 3; i ++) plane[2][i] = -plane[2][i];
			}
			
			setZxy(swz, swx, swy);
			setWeight(sw);
			setValues(eig_val[oindx][oindx],eig_val[vindx][vindx],eig_val[hindx][hindx]); // eigenvalues [0] - thickness, 2 other to detect skinny (poles)
			setVectors   (plane);
			setNumPoints (num_tiles);
			boolean [] plane_sel = null;
			boolean need_clone = true;
			for (int nl = 0; nl < tile_sel.length; nl++){
				if (tile_sel[nl] != null) {
					if (plane_sel == null) {
						plane_sel = tile_sel[nl];
					} else {
						if (need_clone) {
							plane_sel = plane_sel.clone();
							need_clone = false;
						}
						for (int i = 0; i < plane_sel.length; i++){
							plane_sel[i] |=  tile_sel[nl][i];
						}
					} 
				}
			}
			
			setPlaneSelection(plane_sel);
			return disp_str;
		}
		
		public double [][] initMergedValue()
		{
			this.merged_eig_val = new double[8][];
			this.merged_eig_eq = new double[8][];
			this.merged_valid = new boolean[8][];
			return this.merged_eig_val;
		}
		public double [][] getMergedValue()
		{
			return this.merged_eig_val;
		}

		public double [][] getMergedValueEq()
		{
			return this.merged_eig_eq;
		}

		public double [] initMergedValue(int dir, int leng)
		{
			this.merged_eig_val[dir] = new double[leng];
			this.merged_eig_eq[dir] =  new double[leng];
			this.merged_valid[dir] =   new boolean[leng];
			for (int i = 0; i < leng; i++) {
				this.merged_eig_val[dir][i] = Double.NaN;
				this.merged_eig_eq[dir][i] = Double.NaN;
			}
			return getMergedValue(dir);
		}

		public double [] getMergedValue(int dir)
		{
			if (this.merged_eig_val == null) {
				return null;
			}
			return this.merged_eig_val[dir];
		}

		public double [] getMergedValueEq(int dir)
		{
			if (this.merged_eig_eq == null) {
				return null;
			}
			return this.merged_eig_eq[dir];
		}
		
		public double getMergedValue(int dir, int plane)
		{
			if ((this.merged_eig_val == null) ||(this.merged_eig_val[dir] == null)){
				return Double.NaN;
			}
			return this.merged_eig_val[dir][plane];
		}

		public double getMergedValueEq(int dir, int plane)
		{
			if ((this.merged_eig_eq == null) ||(this.merged_eig_eq[dir] == null)){
				return Double.NaN;
			}
			return this.merged_eig_eq[dir][plane];
		}
		
		public void setNeibMatch(int dir, int plane, double value)
		{
			this.merged_eig_val[dir][plane] = value;
		}
		public void setNeibMatchEq(int dir, int plane, double value)
		{
			this.merged_eig_eq[dir][plane] = value;
		}
		
		
		public boolean [][] getMergedValid()
		{
			return this.merged_valid;
		}
		
		public boolean [] getMergedValid(int dir)
		{
			if (this.merged_valid == null) {
				return null;
			}
			return this.merged_valid[dir];
		}
		
		public boolean isMergedValid(int dir, int plane)
		{
			if ((this.merged_valid == null) || (this.merged_valid[dir] == null)){
				return false;
			}
			return this.merged_valid[dir][plane];
		}

		public void setMergedValid(int dir, int plane, boolean valid)
		{
			this.merged_valid[dir][plane] = valid;
		}

		public void setMergedValid(int dir, int plane, boolean valid, int leng)
		{
			if (this.merged_valid == null){
				this.merged_valid = new boolean[8][];
			}
			if (this.merged_valid[dir] == null){
				this.merged_valid[dir] = new boolean[leng];
			}
			this.merged_valid[dir][plane] = valid;
		}
		
		

		public int [] initNeibBest()
		{
			this.neib_best = new int[8];
			for (int i = 0; i < 8; i++) this.neib_best[i] = -1; 
			return this.neib_best;
		}

		public int [] getNeibBest()
		{
			return this.neib_best;
		}

		public int getNumNeibBest()
		{
			if (this.neib_best == null) {
				return 0;
			}
			int num = 0;
			for (int i = 0; i < this.neib_best.length; i++) if (this.neib_best[i] >= 0) num++;
			return num;
		}
		
		public int getNeibBest(int dir)
		{
			if (this.neib_best == null) {
				return -1;
			}
			return this.neib_best[dir];
		}

		public void setNeibBest(int [] vals)
		{
			this.neib_best = vals;
		}

		public void setNeibBest(int dir, int val)
		{
			this.neib_best[dir] = val;
		}

		
		
		
		public void setCorrectDistortions(
				boolean correctDistortions)
		{
			this.correctDistortions = correctDistortions;
		}
		public boolean getCorrectDistortions() {
			return correctDistortions;
		}
		
		public boolean [] getPlaneSelection(){
			return plane_sel;
		}
		public void setPlaneSelection(boolean [] plane_sel){
			this.plane_sel = plane_sel;
		}
		
		public int[] getSTileXY() {
			return sTileXY;
		}
		public double[] getZxy() {
			return zxy;
		}
		
		
		public void setZxy(double[] zxy) {
			this.zxy = zxy;
		}
		public void setZxy(
				double z,
				double x,
				double y) {
			this.zxy = new double [3];
			this.zxy[0] = z;
			this.zxy[1] = x;
			this.zxy[2] = y;
			
		}
		
		public double[][] getVectors() {
			return vectors;
		}
		public double[] getVector() {
			return vectors[0];
		}
		public void setVectors(double[][] vectors) {
			this.vectors = vectors;
		}
		public double[] getValues() {
			return values;
		}
		public double getValue() {
			return values[0];
		}
		public void setValues(double[] values) {
			this.values = values;
		}
		public void setValues(double v1, double v2, double v3) {
			this.values = new double[3];
			this.values[0] = v1;
			this.values[1] = v2;
			this.values[2] = v3;
		}
		public int getNumPoints() {
			return num_points;
		}
		public void setNumPoints(int num_points) {
			this.num_points = num_points;
		}
		public double getWeight() {
			return weight;
		}
		public void setWeight(double weight) {
			this.weight = weight;
		}
		public void scaleWeight(double scale) {
			this.weight *= scale;
		}
		
		public double [] getWorldXYZ(boolean correct_distortions)
		{
			return getWorldXYZ(correct_distortions,0);
		}
		public double [] getCenterPxPy() // supertile center (not plane center)
		{
			double [] px_py = {
					tileSize*(superTileSize * sTileXY[0] + superTileSize/2), // + zxy[1],
					tileSize*(superTileSize * sTileXY[1] + superTileSize/2)}; //  + zxy[2]};
			return px_py;
		}
		public GeometryCorrection getGeometryCorrection() {
			return geometryCorrection;
		}
		public void setGeometryCorrection(GeometryCorrection geometryCorrection) {
			this.geometryCorrection = geometryCorrection;
		}
		
		/**
		 * Get disparity values for the tiles of this supertile as [superTileSize*superTileSize] array
		 * @param useNaN replace unselected tiles with Double.NaN
		 * @return array of disparity values for the plane (not including overlapped areas)
		 */
		public double[] getSinglePlaneDisparity( // getPlaneDisparity(
				boolean useNaN)
		{
			double [] disparities = new double[superTileSize*superTileSize];
			int indx = 0;
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			for (int sy = -superTileSize/2; sy < superTileSize/2; sy++){
				int indx_sel = (2*sy + superTileSize) * superTileSize + superTileSize/2;
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				for (int sx = -superTileSize/2; sx < superTileSize/2; sx++){
					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					if (plane_sel[indx_sel] || !useNaN ||  (plane_sel==null)){
						disparities[indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
					} else {
						disparities[indx] = Double.NaN;
					}
					indx++;
					indx_sel++;
				}
			}
			return disparities;
		}

		public double[] getDoublePlaneDisparity(
				boolean useNaN)
		{
			return getDoublePlaneDisparity(
					true,
					useNaN);
		}
		
		/**
		 * Get disparity values for the tiles of this overlapping supertile as [2*superTileSize * 2*superTileSize] array
		 * @param useWorld calculate disparity in the real world (false - just px, py, disparity plane)
		 * @param useNaN replace unselected tiles with Double.NaN
		 * @return array of disparity values for the plane (not including overlapped areas)
		 */
		public double[] getDoublePlaneDisparity(
				boolean useWorld,
				boolean useNaN)
		{
			double [] disparities = new double[4*superTileSize*superTileSize];
			int indx = 0;
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			double [] pxyc = getCenterPxPy(); // center of this supertile, not plane center 
			for (int sy = -superTileSize; sy < superTileSize; sy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				for (int sx = -superTileSize; sx < superTileSize; sx++){
					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					if (plane_sel[indx] || !useNaN ||  (plane_sel==null)){
						if (useWorld) {
							disparities[indx] = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
									getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
									x + pxyc[0] + zxy[1],
									y + pxyc[1] + zxy[2],
									this.correctDistortions);
						} else {
							disparities[indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
						}
						
					} else {
						disparities[indx] = Double.NaN;
					}
					indx++;
				}
			}
			return disparities;
		}
		
		public double[] getDoublePlaneDisparity(
				int     dir,
				boolean useNaN)
		{
			return getDoublePlaneDisparity(
					true,
					dir,
					useNaN);
		}

		public double[] getDoublePlaneDisparity(
				boolean useWorld,
				int     dir,
				boolean useNaN)
		{
			double [] plane_all = getDoublePlaneDisparity(useWorld,useNaN);
			double [] plane = new double [superTileSize * superTileSize];
			int [] start_index = {
					superTileSize/2,                            // N    4
					superTileSize,                              // NE   8
					(superTileSize + 1) * superTileSize,        // E   72
					(2* superTileSize + 1) * superTileSize,     // SE 136
					(4* superTileSize + 1) * superTileSize / 2, // S  132
					2* superTileSize * superTileSize,           // SW 128
					superTileSize * superTileSize,              // W   64
					0};                                         // NW   0
			for (int y = 0; y < superTileSize; y++) {
				System.arraycopy(plane_all, start_index[dir] + 2 * superTileSize * y, plane,  superTileSize * y, superTileSize);
			}
			return plane;
		}

		public double[] getTriplePlaneDisparity()
		{
			return getTriplePlaneDisparity(true);
		}

		public double[] getTriplePlaneDisparity(
				boolean   useWorld)
		{
			double [] disparities = new double[9*superTileSize*superTileSize];
			int indx = 0;
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			double [] pxyc = getCenterPxPy(); // center of this supertile, not plane center 
			for (int sy = -3 * superTileSize / 2; sy < 3* superTileSize / 2; sy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				for (int sx = -3 * superTileSize/2; sx < 3 * superTileSize / 2; sx++){
					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					if (useWorld) {
						disparities[indx] = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
								getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
								x + pxyc[0] + zxy[1],
								y + pxyc[1] + zxy[2],
								this.correctDistortions);
					} else {
						disparities[indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
					}
					indx++;
				}
			}
			return disparities;
		}

		public double[] getTriplePlaneDisparity(
				int     dir)
		{
			return getTriplePlaneDisparity(true, dir);
		}

		public double[] getTriplePlaneDisparity(
				boolean   useWorld,
				int     dir)
		{
			double [] plane_all = getTriplePlaneDisparity(useWorld);
			double [] plane = new double [superTileSize * superTileSize];
			int [] start_index = {
					superTileSize,                               // N    8
					2 * superTileSize,                           // NE  16
					(3 * superTileSize + 2) * superTileSize,     // E  208
					(3 * superTileSize + 1) * superTileSize * 2, // SE 400
					(6 * superTileSize + 1) * superTileSize ,    // S  392
					6 * superTileSize * superTileSize,           // SW 384
					3 * superTileSize * superTileSize,           // W  192
					0};                                          // NW   0
			for (int y = 0; y < superTileSize; y++) {
				System.arraycopy(plane_all, start_index[dir] + 3 * superTileSize * y, plane,  superTileSize * y, superTileSize);
			}
			return plane;
		}

		public double[][] getDoublePlaneDisparityStrength(
				double [] window,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				int       debugLevel)
		{		return getDoublePlaneDisparityStrength(
				true,
				window,
				use_sel,
				divide_by_area,
				scale_projection,
				debugLevel);
		}
		
		public EigenvalueDecomposition get2dDecomposition()
		{
			double [] vals3d =      getValues();
			double [][] vectors3d = getVectors();
			double [][] acovar = new double [2][2];
			for (int i = 0; i < 2; i++){
				for (int j = i; j < 2; j++){
					acovar[i][j] = 0.0;
					for (int k = 0; k < 3; k++){
						acovar[i][j] += vals3d[k] * vectors3d[k][i+1] * vectors3d[k][j+1]; // 0 - z, disparity == 0
					}
					if (i != j) {
						acovar[j][i] =acovar[i][j];
					}
				}
			}
			Matrix covar = new Matrix(acovar); // 2d, x y only
			return covar.eig();
		}
		
		/**
		 * Get disparity values for the tiles of this overlapping supertile as [2*superTileSize * 2*superTileSize] array
		 * and weights combined from provided window function, optional selection and using ellipsoid projection on the
		 * px, py plane (constant disparity
		 * Sharp weights - when selecting the best match - use exponent of (delta_disp) ^2 ?
		 * Or divide weight by ellipse area?
		 * @param useWorld calculate disparity in the real world (false - just px, py, disparity plane)
		 * @param window null or window function as [2*superTileSize * 2*superTileSize] array
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid 
		 * @return a pair of arrays {disparity, strength}, each [2*superTileSize * 2*superTileSize]
		 */
		// obsolete, convert to another version ?
		public double[][] getDoublePlaneDisparityStrength(
				boolean   useWorld,
				double [] window,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				int       debugLevel)
		{
			double [][] disp_strength = new double[2][4*superTileSize*superTileSize];
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			double    weight = getWeight();
			double k_gauss = 0;
			Matrix val2d = null, vect2d = null;
			if (scale_projection > 0.0){
				EigenvalueDecomposition eig = get2dDecomposition();
				val2d = eig.getD();
				vect2d = eig.getV().transpose();
				k_gauss = 0.5/(scale_projection*scale_projection);
				if (divide_by_area) {
					double area = Math.sqrt(val2d.get(0, 0)*val2d.get(1, 1));
					if (area > 0){
						weight /= area;
					}
				}
			}
			double [] pxyc = getCenterPxPy(); // center of this supertile, not plane center
			if (debugLevel > 0) {
				double world_disp = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
						getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
						zxy[1] + pxyc[0],
						zxy[2] + pxyc[1],
						this.correctDistortions);
				System.out.println("getDoublePlaneDisparityStrength(): zxy = ["+zxy[0]+", "+zxy[1]+", "+zxy[2]+"], disp = "+world_disp);
			}
			int indx = 0;
			for (int sy = -superTileSize; sy < superTileSize; sy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (sy + 0.5) + 0.5 - zxy[2];
				for (int sx = -superTileSize; sx < superTileSize; sx++){
					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					if (useWorld) {
						disp_strength[0][indx] = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
								getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
								x + pxyc[0] + zxy[1],
								y + pxyc[1] + zxy[2],
								this.correctDistortions);
					} else {
						disp_strength[0][indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
					}
					
					double w = weight;
					if (window != null) w *= window[indx];
					if (use_sel && (sel_mask != null) && !(sel_mask[indx])) w = 0.0;
					if ((w > 0.0) && (scale_projection > 0.0)){
						double [] xy = {x,y};
						Matrix vxy = vect2d.times(new Matrix(xy,2)); // verify if it is correct
						double r2 = 0; 
						for (int i = 0; i <2; i++){
							double d = vxy.get(i,0);
							r2 += d * d / val2d.get(i, i);
						}
						w *= Math.exp(-k_gauss*r2);
					}					
					disp_strength[1][indx] = w;
					indx++;
				}
			}
			return disp_strength;
		}

		public double[][] getDoublePlaneDisparityStrength(
				double [] window,
				int       dir,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				double    fraction_uni,
				int       debugLevel)
		{		
			return getDoublePlaneDisparityStrength(
					true,
					window,
					dir,
					use_sel,
					divide_by_area,
					scale_projection,
					fraction_uni,
					debugLevel);
		}
		
		/**
		 * Get disparity values for the tiles of this overlapping supertile as [2*superTileSize * 2*superTileSize] array
		 * and weights combined from provided window function, optional selection and using ellipsoid projection on the
		 * px, py plane (constant disparity
		 * Sharp weights - when selecting the best match - use exponent of (delta_disp) ^2 ?
		 * Or divide weight by ellipse area?
		 * @param useWorld calculate disparity in the real world (false - just px, py, disparity plane)
		 * @param window null or window function as [2*superTileSize * 2*superTileSize] array
		 * @param dir - source tile shift from the target: -1 center, 0 - N, 1 - NE
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid 
		 * @param fraction_uni add fraction of the total weight to each tile
		 * @return a pair of arrays {disparity, strength}, each [2 * superTileSize * 2 * superTileSize], only 1/2 or 1/4 used for offset tiles\
		 * TODO: add a combination of the ellipses and infinite planes?
		 * 
		 */
		public double[][] getDoublePlaneDisparityStrength(
				boolean   useWorld,
				double [] window,
				int       dir,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				double    fraction_uni,
				int       debugLevel)
		{
			double [][] disp_strength = new double[2][4*superTileSize*superTileSize];
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			double    weight = getWeight();
			double k_gauss = 0;
			Matrix val2d = null, vect2d = null;
			if (scale_projection > 0.0){
				EigenvalueDecomposition eig = get2dDecomposition();
				val2d = eig.getD();
				vect2d = eig.getV().transpose();
				k_gauss = 0.5/(scale_projection*scale_projection);
				if (divide_by_area) {
					double area = Math.sqrt(val2d.get(0, 0)*val2d.get(1, 1));
					if (area > 0){
						weight /= area;
					}
				}
			}
			int ss2 = superTileSize;
			int ss4 = 2 * superTileSize;
			
			
			int [][] offsets = {
					// ymin, ymax, xmin,xmax, offsy, offsx
					{  0, ss4,   0, ss4,    0,    0 },  // center	
					{ss2, ss4,   0, ss4, -ss2,    0 },  // N	
					{ss2, ss4,   0, ss2, -ss2,  ss2 },  // NE	
					{  0, ss4,   0, ss2,    0,  ss2 },  // E	
					{  0, ss2,   0, ss2,  ss2,  ss2 },  // SE	
					{  0, ss2,   0, ss4,  ss2,    0 },  // S	
					{  0, ss2, ss2, ss4,  ss2, -ss2 },  // SW	
					{  0, ss4, ss2, ss4,    0, -ss2 },  // W	
					{ss2, ss4, ss2, ss4, -ss2, -ss2 }}; // NW
			int dir1 = dir + 1;
			double [] pxyc = getCenterPxPy(); // center of this supertile, not plane center 
			for (int iy = offsets[dir1][0]; iy < offsets[dir1][1]; iy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (iy - ss2 + 0.5) + 0.5  - zxy[2];
				int oy = iy + offsets[dir1][4]; //vert index in the result tile
				for (int ix = offsets[dir1][2]; ix < offsets[dir1][3]; ix++){
					double x = tileSize * (ix - ss2 + 0.5) + 0.5 - zxy[1];
					int indx = ss4 * oy + ix + offsets[dir1][5];
					int indx_i = iy * ss4 + ix; // input index
					if (useWorld) {
						disp_strength[0][indx] = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
								getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
								x + pxyc[0] + zxy[1],
								y + pxyc[1] + zxy[2],
								this.correctDistortions);
					} else {
						disp_strength[0][indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
					}
					double w = weight;
					if ((w > 0.0) && (scale_projection > 0.0)){
						double [] xy = {x,y};
						Matrix vxy = vect2d.times(new Matrix(xy,2)); // verify if it is correct
						double r2 = 0; 
						for (int i = 0; i <2; i++){
							double d = vxy.get(i,0);
							r2 += d * d / val2d.get(i, i);
						}
						w *= ((1.0 - fraction_uni) * Math.exp(-k_gauss*r2) + fraction_uni);
						
						if (window != null) w *= window[indx_i];
						if (use_sel && (sel_mask != null) && !(sel_mask[indx_i])) w = 0.0;
					}					
					disp_strength[1][indx] = w;
				}
			}
			return disp_strength;
		}

		public double[][] getSinglePlaneDisparityStrength(
				double [] window,
				int       dir,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				int       debugLevel)
		{
			return getSinglePlaneDisparityStrength(
					true,
					window,
					dir,
					use_sel,
					divide_by_area,
					scale_projection,
					debugLevel);
		}
		/**
		 * Get disparity values for the tiles of this overlapping supertile as [superTileSize * superTileSize] array
		 * and weights combined from provided window function, optional selection and using ellipsoid projection on the
		 * px, py plane (constant disparity
		 * Sharp weights - when selecting the best match - use exponent of (delta_disp) ^2 ?
		 * Or divide weight by ellipse area?
		 * @param useWorld calculate disparity in the real world (false - just px, py, disparity plane)
		 * @param window null or window function as [2*superTileSize * 2*superTileSize] array
		 * @param dir - source tile shift from the target: -1 center, 0 - N, 1 - NE
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid 
		 * @return a pair of arrays {disparity, strength}, each [superTileSize * superTileSize], only 1/2 or 1/4 used for offset tiles\
		 * TODO: add a combination of the ellipses and infinite planes?
		 * 
		 */
		public double[][] getSinglePlaneDisparityStrength(
				boolean   useWorld,
				double [] window,
				int       dir,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				int       debugLevel)
		{
			double [][] disp_strength = new double[2][superTileSize*superTileSize];
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			double    weight = getWeight();
			double k_gauss = 0;
			Matrix val2d = null, vect2d = null;
			if (scale_projection > 0.0){
				EigenvalueDecomposition eig = get2dDecomposition();
				val2d = eig.getD();
				vect2d = eig.getV().transpose();
				k_gauss = 0.5/(scale_projection*scale_projection);
				if (divide_by_area) {
					double area = Math.sqrt(val2d.get(0, 0)*val2d.get(1, 1));
					if (area > 0){
						weight /= area;
					}
				}
			}
			int ss1 = superTileSize / 2;
			int ss2 = superTileSize;
			int ss3 = 3 *ss1;
			int ss4 = 2 * superTileSize;
			
			
			int [][] offsets = {
					// ymin, ymax, xmin,xmax, offsy, offsx
					{ss1, ss3, ss1, ss3, -ss1, -ss1 },  // center	
					{ss3, ss4, ss1, ss3, -ss3, -ss1 },  // N	
					{ss3, ss4,   0, ss1, -ss3,  ss1 },  // NE	
					{ss1, ss3,   0, ss1, -ss1,  ss1 },  // E	
					{  0, ss1,   0, ss1,  ss1,  ss1 },  // SE	
					{  0, ss1, ss1, ss3,  ss1, -ss1 },  // S	
					{  0, ss1, ss3, ss4,  ss1, -ss3 },  // SW	
					{ss1, ss3, ss3, ss4, -ss1, -ss3 },  // W	
					{ss3, ss4, ss3, ss4, -ss3, -ss3 }}; // NW
			int dir1 = dir + 1;
			double [] pxyc = getCenterPxPy(); // center of this supertile, not plane center 
			for (int iy = offsets[dir1][0]; iy < offsets[dir1][1]; iy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (iy - ss2 + 0.5) + 0.5  - zxy[2];
				int oy = iy + offsets[dir1][4]; //vert index in the result tile
				for (int ix = offsets[dir1][2]; ix < offsets[dir1][3]; ix++){
					double x = tileSize * (ix - ss2 + 0.5) + 0.5 - zxy[1];
					int indx = ss2 * oy + ix + offsets[dir1][5];
					int indx_i = iy * ss4 + ix; // ss2;
					if (useWorld) {
						disp_strength[0][indx] = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
								getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
								x + pxyc[0] + zxy[1],
								y + pxyc[1] + zxy[2],
								this.correctDistortions);
					} else {
						disp_strength[0][indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
					}
					
					double w = weight;
					if (window != null) w *= window[indx_i];
					if (use_sel && (sel_mask != null) && !(sel_mask[indx_i])) w = 0.0;
					if ((w > 0.0) && (scale_projection > 0.0)){
						double [] xy = {x,y};
						Matrix vxy = vect2d.times(new Matrix(xy,2)); // verify if it is correct
						double r2 = 0; 
						for (int i = 0; i <2; i++){
							double d = vxy.get(i,0);
							r2 += d * d / val2d.get(i, i);
						}
						w *= Math.exp(-k_gauss*r2);
					}					
					disp_strength[1][indx] = w;
				}
			}
			return disp_strength;
		}

		
		
		
		/**
		 * Cross product of 2 3-d vectors as column matrices
		 * @param v1
		 * @param v2
		 * @return v1 x v2 as a column (3x1) matrix
		 */
		public Matrix cross3d (Matrix v1, Matrix v2)
		{
			double [][]av1 = v1.getArray();
			double [][]av2 = v2.getArray();
			double [][] ar = {
					{av1[1][0]*av2[2][0] - av1[2][0]*av2[1][0]},
					{av1[2][0]*av2[0][0] - av1[0][0]*av2[2][0]},
					{av1[0][0]*av2[1][0] - av1[1][0]*av2[0][0]}};
			return new Matrix (ar);
		}

		/** 
		 * Get sin squared of the angle between planes in the real world
		 * @param otherPd other plane data
		 * @param correct_distortions true if the lens distortions should be corrected 
		 * @return sine squared of the angle between normals to the planes
		 */
		public double getWorldSin2(
				PlaneData otherPd,
				boolean correct_distortions)
		{
			Matrix this_wv =  new Matrix(this.getWorldXYZ(correct_distortions, 0),3);
			Matrix other_wv = new Matrix(otherPd.getWorldXYZ(correct_distortions, 0),3);
			Matrix cp = cross3d(this_wv, other_wv);
			double cp2 = cp.transpose().times(cp).get(0, 0);
			double this_wv2 = this_wv.transpose().times(this_wv).get(0, 0);
			double other_wv2 = other_wv.transpose().times(other_wv).get(0, 0);
			return cp2/(this_wv2 * other_wv2);
		}
		public double getWorldSin2(
				PlaneData otherPd)
		{
			return getWorldSin2(otherPd, this.correctDistortions);
		}

		/**
		 * Get distance squared from the other plane and the center of the current "plane" (ellipsoid) 
		 * @param otherPd other plane data
		 * @param correct_distortions true if the lens distortions should be corrected 
		 * @return distance squared from other plane to the (ellipsoid) center of this one (in meters)
		 */
		public double getWorldPlaneDist2(
				PlaneData otherPd,
				boolean correct_distortions)
		{
			Matrix this_center =  new Matrix(this.getCenterXYZ(correct_distortions, 0),3);
			Matrix other_wv =     new Matrix(otherPd.getWorldXYZ(correct_distortions, 0),3);
			double w_dot_w = other_wv.transpose().times(other_wv).get(0, 0);
			double w_dot_p = other_wv.transpose().times(this_center).get(0, 0);
			return (w_dot_w - w_dot_p)*(w_dot_w - w_dot_p)/w_dot_w;
		}

		public double getWorldPlaneDist2(
				PlaneData otherPd)
		{
			return getWorldPlaneDist2(otherPd, this.correctDistortions);
		}
		
		/**
		 * Get squared relative (to the z of the center) distance from the other plane to the center of the current "plane" (ellipsoid) 
		 * @param otherPd other plane data
		 * @return squared ratio of the distance the other plane to the (ellipsoid) center of this one over Z-distance
		 */
		public double getWorldPlaneRDist2(
				PlaneData otherPd)
		{
			double dist2 = getWorldPlaneDist2(otherPd, this.correctDistortions);
			double z =getCenterXYZ(this.correctDistortions, 0)[2];
			return dist2/(z*z);
		}

		
		
		/**
		 * Combine 2 Plane instances using centers, eigenvalues eihenvectors and total weights of this and other PlaneData objects
		 * other plane should already be transformed to the same supertile coordinate system (with getPlaneToThis() method)
		 * @param otherPd PlaneData object to merge with. Should be transformed to the same supertile ( same sTileXY values)
		 * @param scale_other scale total weight of the other PlaneData object
		 * @param ignore_weights assume original weights of the planes where equal, only use scale_other
		 * @param debugLevel debug level
		 * @return PlaneData object representing merged planes with combined weight (scale_other*otherPd.weight + this.weight),
		 * recalculated center, eigenvalues and eigenvectors
		 */
		public PlaneData mergePlaneToThis1(
				PlaneData otherPd,
				double    scale_other,
				boolean   ignore_weights,
				boolean   sum_weights,
				boolean   preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				int       debugLevel)
		{
			return mergePlaneToThis(
					otherPd,
					scale_other,
					1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
					ignore_weights,
					sum_weights,
					preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
					debugLevel);
		}
		public PlaneData mergePlaneToThis(
				PlaneData otherPd,
				double    scale_other,
				double    starWeightPwr,    // Use this power of tile weight when calculating connection cost
				boolean   ignore_weights,
				boolean   sum_weights,
				boolean   preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				int       debugLevel)
		{
			if (debugLevel > 0) {
				System.out.println("mergePlaneToThis()");
			}
			double [][] this_eig_avals = {
					{values[0], 0.0,       0.0},
					{0.0,       values[1], 0.0},
					{0.0,       0.0,       values[2]}};
			double [][] other_eig_avals = {
					{otherPd.values[0], 0.0,               0.0},
					{0.0,               otherPd.values[1], 0.0},
					{0.0,               0.0,               otherPd.values[2]}};
			Matrix this_eig_vals =     new Matrix(this_eig_avals);
			Matrix other_eig_vals =    new Matrix(other_eig_avals);
			Matrix other_eig_vectors = new Matrix(otherPd.vectors).transpose(); // vectors are saved as rows
			Matrix this_eig_vectors =  new Matrix(this.vectors).transpose();    // vectors are saved as rows 
			Matrix this_center    =    new Matrix(this.getZxy(),3);
			Matrix other_center    =   new Matrix(otherPd.getZxy(),3); // should already be relative to this supertile center
			double this_weight = this.weight;
			double other_weight = otherPd.weight;
			if (starWeightPwr == 0){
				ignore_weights = true;
			} else if (starWeightPwr != 1.0){
				this_weight =  Math.pow(this_weight, starWeightPwr);
				other_weight = Math.pow(other_weight,starWeightPwr);
			}
			
			double sum_weight =      scale_other *  other_weight + this_weight;// should be the same for
			double other_fraction = ignore_weights? (scale_other/(scale_other + 1.0)): ((scale_other *  other_weight) / sum_weight); 
			Matrix common_center =  this_center.times(1.0 - other_fraction).plus(other_center.times(other_fraction));
			Matrix other_offset =   other_center.minus(this_center); // other center from this center
			if ((this.values[0] == 0.0) || (otherPd.values[0] == 0.0)) {
				System.out.println("Zero eigenvalue");
				debugLevel = 10;
			}
			if (debugLevel > 0) {
				System.out.println("other_eig_vals");
				other_eig_vals.print(8, 6);
				System.out.println("this_eig_vals");
				this_eig_vals.print(8, 6);
				System.out.println("other_eig_vectors");
				other_eig_vectors.print(8, 6);
				System.out.println("this_eig_vectors");
				this_eig_vectors.print(8, 6);
				System.out.println("other_center");
				other_center.print(8, 6);
				System.out.println("this_center");
				this_center.print(8, 6);
				System.out.println("common_center");
				common_center.print(8, 6);
				System.out.println("other_offset");
				other_offset.print(8, 6);
				System.out.println("other_fraction="+other_fraction);
			}

			double [][] acovar = { // covariance matrix of center masses (not yet scaled by weight)
					{other_offset.get(0,0)*other_offset.get(0,0), other_offset.get(0,0)*other_offset.get(1,0), other_offset.get(0,0)*other_offset.get(2,0)},
					{other_offset.get(1,0)*other_offset.get(0,0), other_offset.get(1,0)*other_offset.get(1,0), other_offset.get(1,0)*other_offset.get(2,0)},
					{other_offset.get(2,0)*other_offset.get(0,0), other_offset.get(2,0)*other_offset.get(1,0), other_offset.get(2,0)*other_offset.get(2,0)}};

			Matrix other_covar = other_eig_vectors.times(other_eig_vals.times(other_eig_vectors.transpose())); 
			Matrix this_covar =  this_eig_vectors.times(this_eig_vals.times(this_eig_vectors.transpose()));
			Matrix covar = (new Matrix(acovar)).times(other_fraction*(1.0-other_fraction)); // only centers with all masses
			if (debugLevel > 0) {
				System.out.println("other_covar");
				other_covar.print(8, 6);
				System.out.println("this_covar");
				this_covar.print(8, 6);
				System.out.println("covar");
				covar.print(8, 6);
			}			

			covar.plusEquals(other_covar.times(other_fraction));
			if (debugLevel > 0) {
				System.out.println("covar with other_covar");
				covar.print(8, 6);
			}
			covar.plusEquals(this_covar.times(1.0 - other_fraction));
			if (debugLevel > 0) {
				System.out.println("covar with other_covar and this_covar");
				covar.print(8, 6);
			}
			if (Double.isNaN(covar.get(0, 0))){
				System.out.println("covar is NaN !");
				covar.print(8, 6);
			}
			// extract new eigenvalues, eigenvectors
			EigenvalueDecomposition eig = covar.eig(); // verify NaN - it gets stuck
			//			eig.getD().getArray(),
			//			eig.getV().getArray(),
			if (debugLevel > 0) {
				System.out.println("eig.getV()");
				eig.getV().print(8, 6);
				System.out.println("eig.getD()");
				eig.getD().print(8, 6);
			}



			double [][] eig_vect = eig.getV().getArray();
			double [][] eig_val = eig.getD().getArray();
			// make towards camera, left coordinate system
			/*			
			int oindx = 0;
			for (int i = 1; i <3; i++){
				if (Math.abs(eig_vect[0][i]) > Math.abs(eig_vect[0][oindx])){
					oindx = i;
				}
			}
			 */			
			int oindx = 0;
			if (preferDisparity) {
				for (int i = 1; i <3; i++){
					if (Math.abs(eig_vect[0][i]) > Math.abs(eig_vect[0][oindx])){
						oindx = i;
					}
				}
			} else {
				for (int i = 1; i <3; i++){
					if (eig_val[i][i] < eig_val[oindx][oindx]){
						oindx = i;
					}
				}

			}

			// select 2 other axes for increasing eigenvalues (so v is short axis, h  is the long one)
			int vindx = (oindx == 0)? 1 : 0;
			int hindx = (oindx == 0)? 2 : ((oindx == 1) ? 2 : 1);
			if (eig_val[vindx][vindx] > eig_val[hindx][hindx]){
				int tmp = vindx;
				vindx = hindx;
				hindx = tmp;
			}

			double [][] plane = {
					{eig_vect[0][oindx],eig_vect[1][oindx],eig_vect[2][oindx]},  // plane normal to camera
					{eig_vect[0][vindx],eig_vect[1][vindx],eig_vect[2][vindx]},  // "horizontal" axis // to detect skinny planes and poles
					{eig_vect[0][hindx],eig_vect[1][hindx],eig_vect[2][hindx]}}; //  "vertical"   axis  // to detect skinny planes and poles
			// Make normal be towards camera (positive disparity), next vector - positive in X direction (right)
			for (int v = 0; v < 2; v++) {
				if (plane[v][v] < 0.0) for (int i = 0; i < 3; i ++) plane[v][i] = -plane[v][i];
			}

			// make  direction last vector so px (x) py (.) disp < 0 (left-hand coordinate system) 
			if (new Matrix(plane).det() > 0){
				for (int i = 0; i < 3; i ++) plane[2][i] = -plane[2][i];
			}


			PlaneData pd = this.clone(); // will copy selections too
			pd.invalidateCalculated();   // real world vectors
			pd.setValues(eig_val[oindx][oindx],eig_val[vindx][vindx],eig_val[hindx][hindx]); // eigenvalues [0] - thickness, 2 other to detect skinny (poles)
			pd.setVectors(plane);

			pd.setZxy(common_center.getColumnPackedCopy()); // set new center
			// what weight to use? cloned is original weight for this supertile
			// or use weighted average like below?
			double new_weight;
			if (sum_weights) {
				new_weight = sum_weight; // normalize while averaging by the caller	
			} else { // how it was before
				new_weight = other_fraction * other_weight + (1.0 - other_fraction) * this_weight;
			}
			if (!ignore_weights && ((starWeightPwr != 1.0))){
				new_weight = Math.pow(new_weight,1.0/starWeightPwr);
			}
			pd.setWeight(new_weight);
			/*
			if (sum_weights) {
				pd.setWeight(sum_weight); // normalize while averaging by the caller	
			} else { // how it was before
				pd.setWeight(other_fraction * other_weight + (1.0 - other_fraction) * this_weight);
			}
			*/
			return pd;
		}

		
		/**
		 * Convert plane data from other supertile to this one (disparity, px, py) for the center of this supertile
		 * through the plane in world coordinates. orientation of the "main" eigenvector (most disparity-like) is
		 * calculated, two other axes preserve just px/py ratio.
		 * Transformed center is calculated for original px, py of the center and recalculated plane to get disparity value
		 * The converted plane data can be combined with the current one to get plane data for combined tiles.
		 * Converted disparity is not limited to positive values. 
		 * @param otherPd other supertile plane data (will not be modified)
		 * @return plane data converted to local supertile disparity, px, py (center and eigenvectors)  
		 */
		public PlaneData getPlaneToThis(
				PlaneData otherPd,
				int       debugLevel)
		{
			PlaneData pd = otherPd.clone(); // TODO: use clone of this, copy only needed info from otherPD
			// keep world vectors from otherPd
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis()");
			}
			
			double [] px_py = getCenterPxPy(); // this supertile center
			double [] px_py_other = pd.getCenterPxPy(); // other supertile center
			double [] wv1 = otherPd.getWorldV12( // 1-st in-plane vector will calculate, but that will not be used below?
					false, // false - first, true - second vector
					this.correctDistortions);
			
			double [] wv2 = otherPd.getWorldV12( // 2-nd in-plane vector
					true, // false - first, true - second vector
					this.correctDistortions);
			double disp = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
					pd.getWorldXYZ(this.correctDistortions), // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
					px_py[0],
					px_py[1],
					this.correctDistortions);
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), px_py = {"+px_py[0]+", "+px_py[1]+"}, px_py_other = {"+px_py_other[0]+", "+px_py_other[1]+"}");
				System.out.println("getPlaneToThis(), disp = "+disp);
				System.out.println("getPlaneToThis(), wv1 = {"+ wv1[0]+", "+ wv1[1]+", "+ wv1[2]+"}");
				System.out.println("getPlaneToThis(), wv2 = {"+ wv2[0]+", "+ wv2[1]+", "+ wv2[2]+"}");
			}
			
			// find world coodinates of the center of tile intersection with the plane
			Matrix xyz = new Matrix(geometryCorrection.getWorldCoordinates(
						px_py[0],
						px_py[1],
						disp,
						this.correctDistortions),3);
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), Matrix xyz=");
				xyz.print(10,6);
			}
			
			// Jacobian to convert world vectors to 
			Matrix img_jacobian =  new Matrix(geometryCorrection.getImageJacobian(
					xyz.getColumnPackedCopy(),
					this.correctDistortions,
					1)); // debug_level
			
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), Matrix img_jacobian=");
				img_jacobian.print(10,6);
			}
			
			
			// now get both in-plane vectors transformed to this supertile disparity, px, py
			Matrix v1 = img_jacobian.times(new Matrix(wv1,3)); // 3 rows, 1 column
			Matrix v2 = img_jacobian.times(new Matrix(wv2,3)); // 3 rows, 1 column
//			Matrix v0 = cross3d(v1,v2); // orthogonal to both - this is a plane normal vector in local supertile disparity, px. py
			Matrix v0 = cross3d(v2,v1); // orthogonal to both - this is a plane normal vector in local supertile disparity, px. py

			
			
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), Matrix v0 =");
				v0.print(10,6);
				System.out.println("getPlaneToThis(), Matrix v1 =");
				v1.print(10,6);
				System.out.println("getPlaneToThis(), Matrix v2 =");
				v2.print(10,6);
			}
			
			
			// normalize v0, update v1, v2 to be orthonormal with v0 
			v0.timesEquals(1.0/v0.normF()); // unity vector;
			v1.minusEquals(v0.times(v0.transpose().times(v1)));
			v1.timesEquals(1.0/v1.normF()); // unity vector;
			v2.minusEquals(v0.times(v0.transpose().times(v2)));
			v2.minusEquals(v1.times(v1.transpose().times(v2)));
			v2.timesEquals(1.0/v2.normF()); // unity vector;
			
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), ortho-normalized Matrix v0 =");
				v0.print(10,6);
				System.out.println("getPlaneToThis(), ortho-normalized Matrix v1 =");
				v1.print(10,6);
				System.out.println("getPlaneToThis(), ortho-normalized Matrix v2 =");
				v2.print(10,6);
			}
			
			
			pd.vectors[0] = v0.getColumnPackedCopy();
			pd.vectors[1] = v1.getColumnPackedCopy();
			pd.vectors[2] = v2.getColumnPackedCopy();
			
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), pd.vectors[0] = {"+ pd.vectors[0][0]+", "+ pd.vectors[0][1]+", "+ pd.vectors[0][2]+"}");
				System.out.println("getPlaneToThis(), pd.vectors[1] = {"+ pd.vectors[1][0]+", "+ pd.vectors[1][1]+", "+ pd.vectors[1][2]+"}");
				System.out.println("getPlaneToThis(), pd.vectors[2] = {"+ pd.vectors[2][0]+", "+ pd.vectors[2][1]+", "+ pd.vectors[2][2]+"}");
			}
			
			
//			double [] dxy = otherPd.getZxy(); // {disparity, px, py}  px, py in pixels, relative to the other supertile center
			double [] dxy = pd.getZxy(); // {disparity, px, py}  px, py in pixels, relative to the other supertile center
			if (debugLevel > 0) {
				System.out.println("getPlaneToThis(), otherPd.getZxy() = {"+ dxy[0]+", "+ dxy[1]+", "+ dxy[2]+"}");
			}
			
			
			
			dxy[1] += px_py_other[0] - px_py[0]; // relative to this tile
			dxy[2] += px_py_other[1] - px_py[1]; // relative to this tile
// find center of plane segment - get intersection of the plane orthogonal to v0 through the point px_py to the point pd.
			dxy[0] = disp - (pd.vectors[0][1]*dxy[1] + pd.vectors[0][2]*dxy[2]) / pd.vectors[0][0];
			if (debugLevel > 0) { 
				System.out.println("getPlaneToThis(), dxy(modified) = {"+ dxy[0]+", "+ dxy[1]+", "+ dxy[2]+"}");
			}
			pd.sTileXY = this.sTileXY; 
			if (pd.vectors[0][0] > 0) {
//				System.out.println("getPlaneToThis(): "+pd.sTileXY[0]+":"+pd.sTileXY[1]+" -> "+pd.vectors[0][0]+", disp = "+disp+
//						", other_det = "+((new Matrix(otherPd.vectors).det()) +", pdr_det = "+((new Matrix(pd.vectors).det()))));
			} else {
//				System.out.println("getPlaneToThis(): "+pd.sTileXY[0]+":"+pd.sTileXY[1]+" -> "+pd.vectors[0][0]+", disp = "+disp+
//						", other_det = "+((new Matrix(otherPd.vectors).det()) +", pdr_det = "+((new Matrix(pd.vectors).det()))));
			}
			copyNeib(this, pd);
			return pd; // make sure pd are updated // "this" is not used. Should it be used instead of pd?  
		}
		
		
		public double [] getWorldV12(
				boolean v2, // false - first, true - second vector
				boolean correct_distortions)
		{
			return getWorldV12(v2, correct_distortions,0);
		}
		
		public double [] getWorldV12(
				boolean v2, // false - first, true - second vector
				boolean correct_distortions,
				int debugLevel)
		{
			if (v2) {
				if (world_v2 == null){
					getWorldXYZ(
							correct_distortions,
							debugLevel);
				}
				return world_v2;
			} else {
				if (world_v1 == null){
					getWorldXYZ(
							correct_distortions,
							debugLevel);
				}
				return world_v1;
			}
		}
		public double [] getCenterXYZ(
				boolean correct_distortions,
				int debugLevel)
		{
			double delta = 0.0001;
			if (center_xyz != null) return center_xyz;
			setCorrectDistortions(correct_distortions);
			// get pixel coordinates of the plane origin point
			double [] px_py = getCenterPxPy();
			
			double px = px_py[0] + zxy[1];
			double py = px_py[1] + zxy[2];
			double disp =  zxy[0];
			center_xyz = geometryCorrection.getWorldCoordinates(
					px,
					py,
					disp,
					this.correctDistortions);
			return center_xyz;
		}

		public double [] getWorldXYZ(
				boolean correct_distortions,
				int debugLevel)
		{
			double delta = 0.0001;
			if (world_xyz != null) return world_xyz;
			setCorrectDistortions(correct_distortions);
			// get pixel coordinates of the plane origin point
			double [] px_py = getCenterPxPy();
			double px = px_py[0] + zxy[1];
			double py = px_py[1] + zxy[2];
			double disp =  zxy[0];
			center_xyz = geometryCorrection.getWorldCoordinates(
					px,
					py,
					disp,
					this.correctDistortions);
			Matrix xyz = new Matrix(center_xyz, 3); // column matrix
			Matrix dpxpy = new Matrix(vectors[0],3); // 3 rows, 1 column
			if (debugLevel > 0){
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"), correctDistortions="+correctDistortions+", xyz= {"+
						xyz.get(0, 0)+","+xyz.get(1, 0)+","+xyz.get(2, 0)+"}, weight = "+getWeight());
				
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): disp, px, py="+disp+","+px+","+py); // + ", reversed:"+dpxpy[0]+","+dpxpy[1]+","+dpxpy[2]);
				
				Matrix xyz1 = new Matrix(geometryCorrection.getWorldCoordinates(
						px + delta * dpxpy.get(1, 0),
						py + delta * dpxpy.get(2, 0),
						disp + delta * dpxpy.get(0, 0),
						this.correctDistortions),3); // column matrix
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+") delta diff: {"+
						((xyz1.get(0, 0)- xyz.get(0, 0))/delta)+","+
						((xyz1.get(1, 0)- xyz.get(1, 0))/delta)+","+
						((xyz1.get(2, 0)- xyz.get(2, 0))/delta)+"}");
			}
			Matrix jacobian =  new Matrix(geometryCorrection.getWorldJacobian(
					px,
					py,
					disp,
					this.correctDistortions,
					(debugLevel > 2)
					));
			
			// convert both orthogonal axes, normalize their cross product
			Matrix v1 = jacobian.times(new Matrix(vectors[1],3)); // 3 rows, 1 column
			Matrix v2 = jacobian.times(new Matrix(vectors[2],3)); // 3 rows, 1 column
			
			world_v1 = v1.getColumnPackedCopy();
			world_v2 = v2.getColumnPackedCopy();
			
			Matrix norm_xyz =  cross3d(v1,v2);
			if (debugLevel > 0){
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): norm_xyz={"+
						norm_xyz.get(0, 0)+", "+norm_xyz.get(1, 0)+", "+norm_xyz.get(2, 0)+"}, (dpxpy={"+
						dpxpy.get(0, 0)+   ", "+dpxpy.get(1, 0)+   ", "+dpxpy.get(2, 0)+"})");
			}
			if (debugLevel > 2){
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): dpxpy=");
				dpxpy.print(10, 6); // w,d
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): jacobian=");
				jacobian.print(10, 6); // w,d
			}
			if (debugLevel > 2){
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): norm_xyz=");
				norm_xyz.print(10, 6); // w,d
			}
			if (debugLevel > 2){
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): image jacobian=");
				Matrix img_jacobian =  new Matrix(geometryCorrection.getImageJacobian(
						xyz.getColumnPackedCopy(),
						this.correctDistortions,
						1));
				img_jacobian.print(10, 6); // w,d
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): jacobian.times(image_jacobian)=");
				jacobian.times(img_jacobian).print(10, 6); // w,d
			}
			
			norm_xyz.timesEquals(1.0/norm_xyz.normF()); // unity normal vector;
			if (debugLevel > 1){
				System.out.println("+getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): unit plane normal={"+
						norm_xyz.get(0, 0)+", "+norm_xyz.get(1, 0)+", "+norm_xyz.get(2, 0)+"})");
				double dotprod = xyz.transpose().times(norm_xyz).get(0,0);
				Matrix wn = norm_xyz.times(dotprod);
				System.out.println(":getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): xyz.transpose().times(norm_xyz).get(0,0) ="+dotprod);
				System.out.println("?getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"):  plane normal={"+
						wn.get(0, 0)+", "+wn.get(1, 0)+", "+wn.get(2, 0)+"})");
			}
			
			
			if (debugLevel > 2){
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): norm_xyz (normalized) =");
				norm_xyz.print(10, 6); // w,d
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"): xyz.times(norm_xyz.transpose()).get(0,0) ="+xyz.times(norm_xyz.transpose()).get(0,0));
			}
			
			// convert plane normal vector to world coordinates
			//world_xyz
			world_xyz = norm_xyz.times((xyz.transpose().times(norm_xyz).get(0,0))).getColumnPackedCopy();
			return world_xyz;
		}
	}
	
	//TODO: Remove below methods and promote PlaneData (no TilePlanes) after tested
	
	/**
	 * Calculate covariance matrix for a subset of tile data (disparities)
	 * Subtract weight floor from weight
	 * @param data    data array - square (2*stSize) * (2*stSize)  
	 * @param weight  per sample weight (should have floor already subtracted)
	 * @param select  sample selection
	 * @return covariance (diagonal) matrix: [0]: disparity, [1]: d<disparity>/dx, [2]: d<disparity>/dy,
	 */
	public double [][][] getCovar(
			double []  data,
			double []  weight,
			boolean [] select,
			double     plDispNorm, //  Normalize disparities to the average if above
			int        debugLevel){
		double mindet = 1E-15;
		int stSize2 = 2 * stSize;
//		Matrix covar = new Matrix(3,3);
		double [][] acovar = new double [3][3];
		int numPoints = 0;
		double sw =0.0, swz = 0.0, swx = 0.0, swy = 0.0;

		for (int indx = 0; indx < data.length; indx++){
			if (select[indx] && (weight[indx] > 0)){
				numPoints++;
				double w = weight[indx];
				double d = data[indx];
				// referencing samples to centers of pixels
				double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5; // in pixels, not in tiles
				double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5;
				sw  += w;
				swz += w * d;
				swx += w * x;
				swy += w * y;
			}
		}
		if (sw == 0.0) {
			return null;
		}
		swz /= sw;
		swx /= sw;
		swy /= sw;
		
		
		// TODO: scale disparity to make same scale for 3 axes?
		
		double kz = ((plDispNorm > 0.0) && (swz > plDispNorm)) ? (plDispNorm / swz) : 1.0; 
		if (debugLevel > 0){
			System.out.println("getCovar(): sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy+" kz = "+kz);
		}
		for (int indx = 0; indx < data.length; indx++){
			if (select[indx] && (weight[indx] > 0)){
				double w = weight[indx] / sw;
				double d = kz * (data[indx] - swz);
				double wd = w*d;
				double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5 - swx;
				double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5 - swy;
				acovar [0][0] += wd * d;
				acovar [0][1] += wd * x;
				acovar [0][2] += wd * y;
				acovar [1][1] += w * x * x;
				acovar [1][2] += w * x * y;
				acovar [2][2] += w * y * y;
			}
		}
		acovar [1][0] = acovar [0][1]; 
		acovar [2][0] = acovar [0][2]; 
		acovar [2][1] = acovar [1][2];
		Matrix covar = new Matrix(acovar);
//		if (Math.abs(covar.det()) < mindet){
//			debugLevel = 5;
//		}

		EigenvalueDecomposition eig = covar.eig();
		if (Double.isNaN(eig.getV().get(0, 0))){
			System.out.println("getCovar(): Double.isNaN(eig.getV().get(0, 0))");
			debugLevel = 20;
		}
		if (debugLevel > 3){
			System.out.println("getCovar(): sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy +", covar.det() = "+covar.det());
			System.out.println("getCovar(): covarianvce matrix, number of used points:"+numPoints);
			covar.print(10, 6); // w,d
			System.out.println("getCovar(): eigenvalues");
			eig.getD().print(10, 6); // w,d
			System.out.println("getCovar(): eigenvectors");
			eig.getV().print(10, 6); // w,d
		}
		if ((eig.getD().get(0, 0) == 0.0) || (Math.abs(covar.det()) < mindet)) {
			return null; // testing with zero eigenvalue
			// Problem with zero eigenvalue is with derivatives and coordinate conversion
		}
		double [][][] rslt = {
				eig.getD().getArray(),
				eig.getV().getArray(),
				{
					{sw,kz,numPoints},
					{swz, swx, swy}}};
		return rslt;
	}
	// TODO: obsolete - remove
	public PlaneData getPlane(
			int [] sTileXY,
			double []  data,
			double []  weight,
			boolean [] select, // null OK, will enable all tiles
			boolean    correctDistortions,
			boolean    preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			int        debugLevel){
		if (select == null) {
			select = new boolean [4*stSize];
			for (int i = 0; i < select.length; i++) select[i] = true;
		}
		if (debugLevel > 0){
			System.out.println("getPlane()");
		}
		double [][][] rslt = getCovar(
				data,
				weight,
				select,
				0.0,
				debugLevel); // debugLevel1); //0); // debugLevel);
		if (rslt == null) return null;
		int       numPoints =  (int) rslt[2][0][2];
		double    swc =  rslt[2][0][0];
		double [] szxy = rslt[2][1];
		double [][] eig_val =  rslt[0];
		double [][] eig_vect = rslt[1];
		// find vector most orthogonal to view // (anyway it all works with that assumption), make it first
		// TODO normalize to local linear scales
		int oindx = 0;
		if (preferDisparity) {
			for (int i = 1; i <3; i++){
				if (Math.abs(eig_vect[0][i]) > Math.abs(eig_vect[0][oindx])){
					oindx = i;
				}
			}
		} else {
			for (int i = 1; i < 3 ; i++){
				if (eig_val[i][i] < eig_val[oindx][oindx]){
					oindx = i;
				}
			}
		}
		if (eig_val[oindx][oindx] == 0.0){
			System.out.println("getPlane(): zero eigenvalue!!");
		}
		
		/*
		// Find two other axis - "mostly X" (horizontal) and "mostly Y" (vertical) 
		int vindx = (oindx == 0)? 1 : 0;
		int hindx = (oindx == 0)? 2 : ((oindx == 1) ? 2 : 1);
		if (Math.abs(eig_vect[2][vindx]) < Math.abs(Math.abs(eig_vect[2][hindx]))){
			int tmp = vindx;
			vindx = hindx;
			hindx = tmp;
		}
		*/
		// select 2 other axes for increasing eigenvalues (so v is short axis, h  is the long one)
		int vindx = (oindx == 0)? 1 : 0;
		int hindx = (oindx == 0)? 2 : ((oindx == 1) ? 2 : 1);
		if (eig_val[vindx][vindx] > eig_val[hindx][hindx]){
			int tmp = vindx;
			vindx = hindx;
			hindx = tmp;
		}
		
		PlaneData pd = new PlaneData(
				sTileXY,
				this.tileSize,
				this.stSize,
				this.geometryCorrection,
				correctDistortions);
		pd.setZxy(szxy);
		
//		pd.setValues(eig_val[oindx][oindx],eig_val[hindx][hindx],eig_val[vindx][vindx]); // eigenvalues [0] - thickness, 2 other to detect skinny (poles)
		pd.setValues(eig_val[oindx][oindx],eig_val[vindx][vindx],eig_val[hindx][hindx]); // eigenvalues [0] - thickness, 2 other to detect skinny (poles)

		double [][] plane = {
				{eig_vect[0][oindx],eig_vect[1][oindx],eig_vect[2][oindx]},  // plane normal to camera
				{eig_vect[0][vindx],eig_vect[1][vindx],eig_vect[2][vindx]},  // "horizontal" axis // to detect skinny planes and poles
				{eig_vect[0][hindx],eig_vect[1][hindx],eig_vect[2][hindx]}}; //  "vertical"   axis  // to detect skinny planes and poles
		/*
		// Make normal be towards camera (positive disparity), next vector - positive in X direction (right), last one - in positive Y (down)
		for (int v = 0; v <3; v++) {
			if (plane[v][v] < 0.0) for (int i = 0; i < 3; i ++) plane[v][i] = -plane[v][i];
		}
		*/
		// Make normal be towards camera (positive disparity), next vector - positive in X direction (right)
		for (int v = 0; v < 2; v++) {
			if (plane[v][v] < 0.0) for (int i = 0; i < 3; i ++) plane[v][i] = -plane[v][i];
		}
		
		// make  direction last vector so px (x) py (.) disp < 0 (left-hand coordinate system) 
		if (new Matrix(plane).det() > 0){
			for (int i = 0; i < 3; i ++) plane[2][i] = -plane[2][i];
		}
		
		pd.setVectors   (plane);
		pd.setNumPoints (numPoints);
		pd.setWeight    (swc);
		pd.setPlaneSelection(select);
		return pd;
	}
	
	public PlaneData removePlaneOutliers(
			PlaneData  pd,     // already found or null 
			int [] sTileXY,    // may be null if pd is not null 
			double []  data,
			double []  weight,
			boolean [] select, // will be modified
			boolean    correctDistortions,
			double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
			int        maxRemoved,  // maximal number of tiles to remove (not a constant)
			int        minLeft,     // minimal number of tiles to keep
			boolean    preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
			int        debugLevel){
		int stSize2 = 2 * stSize;
		if (pd == null) {
			pd = getPlane(
					sTileXY, 
					data,
					weight,
					select, // null OK
					correctDistortions,
					preferDisparity,
					debugLevel);
		} else if (select != null){
			pd.setPlaneSelection(select);
		}
		if (pd == null){
			return null; // zero eigenvalues
		}
		if (maxRemoved > (pd.getNumPoints() - minLeft)) maxRemoved = pd.getNumPoints() - minLeft;
		int numRemoved = 0;
		for (; (pd.getValue() > targetEigen) && (numRemoved < maxRemoved); numRemoved++){
			if (debugLevel > 2){
				System.out.println("removePlaneOutliers("+sTileXY[0]+":"+sTileXY[1]+"): numRemoved = "+numRemoved+" eigenValue = "+pd.getValue()+" target = "+targetEigen);
			}
			// make a plane and find the worst (largest disparity difference) tile
			// z = -(x*Vx + y*Vy)/Vz
			double worst_d2 = 0.0;
			int worst_index =  -1;
			double [] v = pd.getVector();
			double [] zxy0 = pd.getZxy();
			for (int indx = 0; indx < data.length; indx++){
				if (select[indx] && (weight[indx] > 0)){
//					double w = weight[indx];
					double x = ((indx % stSize2) - stSize) - zxy0[1];
					double y = ((indx / stSize2) - stSize) - zxy0[2];
					double d = data[indx];
					d -= zxy0[0];
					d += (x * v[1]+y*v[2])/v[0];
					double d2 = d*d;
					if (d2 > worst_d2){
						worst_d2 = d2;
						worst_index = indx;
					}
				}
			}
			if (worst_index < 0) {
				System.out.println("This is a BUG in removePlaneOutliers()");
				break;
			}
			select[worst_index] = false;
			if (debugLevel > 2){
				System.out.println("removePlaneOutliers() worst_index = " + worst_index);
			}
			
			pd = getPlane(    // re-calculate with point removed
					pd.getSTileXY(),
					data,
					weight,
					select,
					correctDistortions,
					preferDisparity,
					debugLevel);
			if (pd == null) {
				return null;
			}
		}
		return pd;
	}
	
}
