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
import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;

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
		double []   world_v2 =   null; // world in-plane vector, corresponding to vectors[2]
//		double []   daxy      =  null; // disparity and 2 relative angles (ax and ay) corresponding to fisheye view, near (0,0) scale is pixel size
		// for now keeping both weighted and equal weight merged value - later remove less useful
		double [][] merged_eig_val = null; // for each of the directions (N, NE, .. NW) quality match for each layer
		double [][] merged_eig_eq =  null; // for each of the directions (N, NE, .. NW) quality match for each layer - ignoring weights
		boolean [][] merged_valid = null; // for each of the directions (N, NE, .. NW) if it is possible to connect with link swaps
		boolean [][] merged_strong_valid = null; // for each of the directions (N, NE, .. NW) if it is possible to connect with link swaps (no "discounts"
		// for low weight(s)

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
		double       min_weight =            0.0;  // minimal weight of the ellipsoid
		int          min_tiles =             10;
		double       dispNorm =              5.0;  //  Normalize disparities to the average if above
		boolean      smplMode =              true;   // Use sample mode (false - regular tile mode)

		MeasuredLayersFilterParameters mlfp = new MeasuredLayersFilterParameters();     // filter parameters

//		double       measured_strength_pow = 1.0;
//		double       strength_floor =        0.0;

//		int          smplSide = 2;      // Sample size (side of a square)
//		int          smplNum  = 3;      // Number after removing worst
//		double       smplRms  = 0.1;    // Maximal RMS of the remaining tiles in a sample
//		boolean      smplWnd =  false;   // Use sample mode (false - regular tile mode)

//		double       max_abs_tilt  = 2.0; // Maximal absolute tilt in pixels/tile
//		double       max_rel_tilt  = 0.2; // Maximal relative tilt in pixels/tile/disparity
//		double       damp_tilt  =    0.001; // Damp tilt to handle insufficient  (co-linear)data
//		double       min_tilt_disp = 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//		double       transition    = 1.0; // Mode transition range (between tilted and maximal disparity)
//		int          far_mode  =     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//		double       far_power =     1.0; // Raise disparity to this power before averaging for far objects



		double [] starValueWeight = null;
		PlaneData starPlane = null;
		PlaneData nonexclusiveStar =   null;
		PlaneData nonexclusiveStarEq = null;

		double    conn_density = Double.NaN; //

		boolean      preferDisparity = false;

		// alternative "plane" calculations in the world coordinates
		double []    wxyz =        null; // [3] - plane center point when calculated in world coordinates (x, y , z)
		double [][]  wvectors =    null; // [3][3] - eigenvectors calculated in the real world
		double []    wvalues =     null; // [3] -eigenvalues calculated in the real world
		double [][]  merged_weig_val = null; // for each of the directions (N, NE, .. NW) quality match for each layer
		double [][]  merged_weig_eq =  null; // for each of the directions (N, NE, .. NW) quality match for each layer - ignoring weights
		double [][]  link_costs; // for each of the directions (N, NE, .. NW) composite cost of connection

		int mark0 = -1; // just for temporary labeling the plane

		@Override
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
			if (this.vectors != null) {
				pd.vectors = new double[3][];
				pd.vectors[0] = this.vectors[0].clone();
				pd.vectors[1] = this.vectors[1].clone();
				pd.vectors[2] = this.vectors[2].clone();
			}
			// Adding cloning of the calculated center_xyz and world_xyz (normal). Check that it did not break anything
			if (this.center_xyz != null) pd.center_xyz = this.center_xyz.clone();
			if (this.world_xyz != null) pd.world_xyz = this.world_xyz.clone();
			if (this.world_v1 != null) pd.world_v1 = this.world_v1.clone();
			if (this.world_v2 != null) pd.world_v2 = this.world_v2.clone();

			if (this.measuredLayers != null) pd.measuredLayers = this.measuredLayers;

			pd.setMeasSelection(this.measuredSelection);

			if (this.sel_mask != null) pd.sel_mask = this.sel_mask.clone();

			pd.min_weight =            this.min_weight;
			pd.min_tiles =             this.min_tiles;
			pd.dispNorm =              this.dispNorm;

			pd.smplMode =              this.smplMode;

			pd.mlfp = this.mlfp.clone();
//			pd.measured_strength_pow = this.measured_strength_pow;
//			pd.strength_floor =        this.strength_floor;

//			pd.smplSide =              this.smplSide;
//			pd.smplNum =               this.smplNum;
//			pd.smplRms =               this.smplRms;

//			pd.max_abs_tilt =          this.max_abs_tilt;
//			pd.max_rel_tilt =          this.max_rel_tilt;
//			pd.damp_tilt =             this.damp_tilt;
//			pd.min_tilt_disp =         this.min_tilt_disp;
//			pd.transition =            this.transition;
//			pd.far_mode =              this.far_mode;
//			pd.far_power =             this.far_power;

			pd.preferDisparity =       this.preferDisparity;

			copyNeib(this,pd);
			copyStar(this,pd);
/*
			if (starValueWeight != null){
				pd.starValueWeight = starValueWeight.clone();
			}
			if (this.starPlane !=          null) pd.starPlane =          this.starPlane;
			if (this.nonexclusiveStar !=   null) pd.nonexclusiveStar =   this.nonexclusiveStar;
			if (this.nonexclusiveStarEq != null) pd.nonexclusiveStarEq = this.nonexclusiveStarEq;
			pd.conn_density =          this.conn_density;
*/
			pd.num_points =            this.num_points; // restore, maybe remove from copy_neib?
//
			if (this.wxyz != null)        pd.wxyz =        this.wxyz.clone();

			if (this.wvalues != null)     pd.wvalues =     this.wvalues.clone();

			if (this.wvectors != null) {
				pd.wvectors = new double[3][];
				pd.wvectors[0] = this.wvectors[0].clone();
				pd.wvectors[1] = this.wvectors[1].clone();
				pd.wvectors[2] = this.wvectors[2].clone();
			}
			pd.mark0 = this.mark0;
			return pd;
		}

		public int getMark0(){
			return this.mark0;
		}

		public void setMark0(int mark0){
			this.mark0 = mark0;
		}

		public boolean getPreferDisparity(){
			return this.preferDisparity;
		}

		public String getNeibString()
		{
			if (neib_best == null) {
				return "[      undefined       ] ";
			}
			String s = "[";
			for (int dir = 0; dir < 8; dir++){
				s += (neib_best[dir]>=0) ? neib_best[dir]:"x";
				if (dir < 7) s += ", ";
			}
			s+= "] ";
			return s;

		}
/*
		public boolean isHorizontalW(){
			if (wvectors != null){
				return (Math.abs(wvectors[0][1]) > 0.99) && (checkBadPlate(false) < 0.2);
			}
			return false;
		}
*/
		public boolean isHorizontal(){
			if (world_xyz != null){
				double norm = world_xyz[1] / Math.sqrt(world_xyz[0]*world_xyz[0] + world_xyz[1]*world_xyz[1] + world_xyz[2]*world_xyz[2]);
				return (Math.abs(norm) > 0.99) && (checkBadPlate(false) < 0.2);
			}
			return false;
		}

		public double get2dRatio(){
			if (wvalues != null) return Math.sqrt(wvalues[1]/wvalues[0]);
			return Double.NaN;
		}
		/**
		 * Verify plane normal in real world and calculated from pixels/disparity
		 * Use to check consistency
		 * @return sin squared of the angle between 2 normals, or NaN if there is no data to calculate it
		 */
		public double checkBadPlate(boolean force)
		{
			if (!force && (world_xyz == null)) return Double.NaN;
			if (wvectors == null) return Double.NaN;
			Matrix norm_disp =  new Matrix(this.getWorldXYZ(this.correctDistortions, 0),3); // normal to plane from disparity space
			Matrix norm_world =  new Matrix(wvectors[0],3); // normal to plane from disparity space
			Matrix cp = cross3d(norm_disp, norm_world);
			double cp2 = cp.transpose().times(cp).get(0, 0);
			double this_wv2 = norm_disp.transpose().times(norm_disp).get(0, 0);
			double other_wv2 = norm_world.transpose().times(norm_world).get(0, 0);
			return cp2/(this_wv2 * other_wv2);
		}

		public double checkBadStick(boolean force)
		{
			if (!force && (world_xyz == null)) return Double.NaN;
			if (wvectors == null) return Double.NaN;
			Matrix long_disp =  new Matrix(this.getWorldV12(true,this.correctDistortions, 0),3); // normal to plane from disparity space
			Matrix long_world =  new Matrix(wvectors[2],3); // normal to plane from disparity space
			Matrix cp = cross3d(long_disp, long_world);
			double cp2 = cp.transpose().times(cp).get(0, 0);
			double this_wv2 = long_disp.transpose().times(long_disp).get(0, 0);
			double other_wv2 = long_world.transpose().times(long_world).get(0, 0);
			return cp2/(this_wv2 * other_wv2);
		}

		@Override
		public String toString()
		{

			String s = "          ";
			s += getNeibString();
			if (isHorizontal()){
				s+= "HORIZONTAL ";
			}
			s += String.format("2d=%4.1f ", get2dRatio());
			s += String.format( "np=%3d weight= %8.5f", num_points, weight);
			if (starValueWeight != null) s += String.format(" star=[%8.5f, %8.5f]", starValueWeight[0], starValueWeight[1]);
			else                         s +=               " star=  null";
			s += String.format(" dens=%8.5f", conn_density);
			double [] px_py = getCenterPxPy();
			if (zxy != null) s += String.format("\nzxy =     [%8.3f, %8.3f, %8.3f] (pix)",zxy[0],zxy[1]+px_py[0],zxy[2]+px_py[1]);
			else  s +=                          "\nzxy =     null";
			if (values != null)	s += String.format(", values = [%8.5f, %8.4f, %8.3f] (%8.3f) pix^2",values[0],values[1],values[2], getNormValue());
			else  s +=                             " values = null";
			if (vectors != null) s += String.format("\nvectors = [%8.5f, %8.5f, %8.5f], [%8.5f, %8.5f, %8.5f], [%8.5f, %8.5f, %8.5f]",
					vectors[0][0],vectors[0][1],vectors[0][2], vectors[1][0],vectors[1][1],vectors[1][2], vectors[2][0],vectors[2][1],vectors[2][2]);
			if (center_xyz != null) s += String.format("\ncenter =  [%8.2f, %8.2f, %8.2f]",center_xyz[0],center_xyz[1],center_xyz[2]);
			else  s +=                                 "\ncenter =   null";

			if (world_xyz != null)  s += String.format(" normal = [%8.2f, %8.2f, %8.2f] (m)",world_xyz[0],world_xyz[1],world_xyz[2]);
			else s +=                                  " normal =   null";

			double bad_plate = checkBadPlate(false);
			if (!Double.isNaN(bad_plate)){
				s+=String.format(" bad_plate=%6.4f", bad_plate);
			}
			double bad_stick = checkBadStick(false);
			if (!Double.isNaN(bad_plate)){
				s+=String.format(" bad_stick=%6.4f", bad_stick);
			}
			double [] world_xyz_w = getWorldXYZFromWorld();
			if (world_xyz_w != null) s += String.format("\n%34s world normal = [%8.2f, %8.2f, %8.2f] (m)","",world_xyz_w[0],world_xyz_w[1],world_xyz_w[2]);

			if (wxyz != null)       s += String.format("\nwxyz =    [%8.2f, %8.2f, %8.2f] (m)",wxyz[0],wxyz[1],wxyz[2]);
			else s +=                                  "\nwxyz =  null";
			if (wvalues != null)    s += String.format(" wvals = [%8.4f, %8.3f, %8.2f] (m^2)",wvalues[0],wvalues[1],wvalues[2]);
			else  s +=                                 " wvals =  null";
			if (wvectors != null) s += String.format("\nwvect =   [%8.5f, %8.5f, %8.5f], [%8.5f, %8.5f, %8.5f], [%8.5f, %8.5f, %8.5f]",
					wvectors[0][0],wvectors[0][1],wvectors[0][2], wvectors[1][0],wvectors[1][1],wvectors[1][2], wvectors[2][0],wvectors[2][1],wvectors[2][2]);
			if (nonexclusiveStar != null){
				s+= "\nweighted: ";
				s+= nonexclusiveStar.getNeibString();
				if (nonexclusiveStar.isHorizontal()){
					s+= "HORIZONTAL ";
				}
				s += String.format("2d=%4.1f ", nonexclusiveStar.get2dRatio());
				s += String.format( "np=%3d weight= %8.5f", nonexclusiveStar.num_points, nonexclusiveStar.weight);
				double [] ne_px_py = getCenterPxPy();
				if (nonexclusiveStar.center_xyz != null) s += String.format("\n--center =[%8.2f, %8.2f, %8.2f]",
						nonexclusiveStar.center_xyz[0],
						nonexclusiveStar.center_xyz[1], //  + ne_px_py[0],
						nonexclusiveStar.center_xyz[2]); //  + ne_px_py[1]);
				else  s +=                                 "\n--ncenter =   null";
				if (nonexclusiveStar.world_xyz != null)  s += String.format(" normal = [%8.2f, %8.2f, %8.2f] (m)",
						nonexclusiveStar.world_xyz[0],nonexclusiveStar.world_xyz[1],nonexclusiveStar.world_xyz[2]);
				else s +=                                  " normal =   null";
				bad_plate = nonexclusiveStar.checkBadPlate(false);
				if (!Double.isNaN(bad_plate)){
					s+=String.format(" bad_plate=%6.4f", bad_plate);
				}
				bad_stick = nonexclusiveStar.checkBadStick(false);
				if (!Double.isNaN(bad_plate)){
					s+=String.format(" bad_stick=%6.4f", bad_stick);
				}

			}
			if (nonexclusiveStarEq != null){
				s+= "\nequalized:";
				s+= nonexclusiveStarEq.getNeibString();
				if (nonexclusiveStarEq.isHorizontal()){
					s+= "HORIZONTAL ";
				}
				s += String.format("2d=%4.1f ", nonexclusiveStarEq.get2dRatio());
				s += String.format( "np=%3d weight= %8.5f", nonexclusiveStarEq.num_points, nonexclusiveStarEq.weight);
				double [] ne_px_py = getCenterPxPy();
				if (nonexclusiveStarEq.center_xyz != null) s += String.format("\n--center =[%8.2f, %8.2f, %8.2f]",
						nonexclusiveStarEq.center_xyz[0],
						nonexclusiveStarEq.center_xyz[1], //  + ne_px_py[0],
						nonexclusiveStarEq.center_xyz[2]); //  + ne_px_py[1]);
				else  s +=                                 "\n--ncenter =   null";
				if (nonexclusiveStarEq.world_xyz != null)  s += String.format(" normal = [%8.2f, %8.2f, %8.2f] (m)",
						nonexclusiveStarEq.world_xyz[0],nonexclusiveStarEq.world_xyz[1],nonexclusiveStarEq.world_xyz[2]);
				else s +=                                  " normal =   null";
				bad_plate = nonexclusiveStarEq.checkBadPlate(false);
				if (!Double.isNaN(bad_plate)){
					s+=String.format(" bad_plate=%6.4f", bad_plate);
				}
				bad_stick = nonexclusiveStarEq.checkBadStick(false);
				if (!Double.isNaN(bad_plate)){
					s+=String.format(" bad_stick=%6.4f", bad_stick);
				}
			}
			s+="\n";
			if (link_costs != null){
				for (int dir = 0; dir < link_costs.length; dir++){
					s+=String.format("dir=%d: ", dir);
					if (link_costs[dir] != null) {
						int best_np = -1;
						for (int np = 0; np < link_costs[dir].length; np++){
							if (!Double.isNaN(link_costs[dir][np]) && ((best_np < 0) || (link_costs[dir][np] < link_costs[dir][best_np]))){
								best_np = np;
							}
						}
						for (int np = 0; np < link_costs[dir].length; np++){
							if (np == best_np){
								s+=String.format("%7.3f[%d] ", link_costs[dir][np],np);

							} else {
								s+=String.format("%7.3f    ", link_costs[dir][np]);
							}
						}
					}
					s+="\n";
				}
			}
//			s+="\n";
			return s;
		}

		public PlaneData getNonexclusiveStar()
		{
			return this.nonexclusiveStar;
		}
		public PlaneData getNonexclusiveStarFb() // fallback to this plane if nonexclusiveStar is not available
		{
			if (this.nonexclusiveStar != null) return this.nonexclusiveStar;
			return this;
		}

		public void setNonexclusiveStar( PlaneData pd)
		{
			this.nonexclusiveStar = pd;
		}

		public PlaneData getNonexclusiveStarEq()
		{
			return this.nonexclusiveStarEq;
		}

		public PlaneData getNonexclusiveStarEqFb() // fallback to this plane if nonexclusiveStarEq is not available
		{
			if (this.nonexclusiveStarEq != null) return this.nonexclusiveStarEq;
			return this;
		}


		public void setNonexclusiveStarEq( PlaneData pd)
		{
			this.nonexclusiveStarEq = pd;
		}

		public PlaneData getStarPlane()
		{
			return this.starPlane;
		}
		public void setStarPlane( PlaneData pd)
		{
			this.starPlane = pd;
		}


		public double getConnectionDensity(){
			return conn_density;
		}


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
		public double [] getStarValueWeightDensity()
		{
			double [] vwd = {starValueWeight[0], starValueWeight[1], conn_density};
			return vwd;
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

		public void copyStar(
				PlaneData src,
				PlaneData dst)


		{
			if (src.starValueWeight != null){
				dst.starValueWeight = src.starValueWeight.clone();
			}
			if (src.starPlane !=          null) dst.starPlane =          src.starPlane;
			if (src.nonexclusiveStar !=   null) dst.nonexclusiveStar =   src.nonexclusiveStar;
			if (src.nonexclusiveStarEq != null) dst.nonexclusiveStarEq = src.nonexclusiveStarEq;
			dst.conn_density =          src.conn_density;
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

			if (src.merged_weig_val != null){
				dst.merged_weig_val = src.merged_weig_val.clone();
				for (int i = 0; i < src.merged_weig_val.length; i++){
					if (src.merged_weig_val[i] != null){
						dst.merged_weig_val[i] = src.merged_weig_val[i].clone();
					}
				}
			}

			if (src.merged_weig_eq != null){
				dst.merged_weig_eq = src.merged_weig_eq.clone();
				for (int i = 0; i < src.merged_weig_eq.length; i++){
					if (src.merged_weig_eq[i] != null){
						dst.merged_weig_eq[i] = src.merged_weig_eq[i].clone();
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

			if (src.merged_strong_valid != null){
				dst.merged_strong_valid = src.merged_strong_valid.clone();
				for (int i = 0; i < src.merged_strong_valid.length; i++){
					if (src.merged_strong_valid[i] != null){
						dst.merged_strong_valid[i] = src.merged_strong_valid[i].clone();
					}
				}
			}

			if (src.link_costs != null){
				dst.link_costs = src.link_costs.clone();
				for (int i = 0; i < src.link_costs.length; i++){
					if (src.link_costs[i] != null){
						dst.link_costs[i] = src.link_costs[i].clone();
					}
				}
			}


			if (src.neib_best != null) dst.neib_best = src.neib_best.clone();

			// also copy original plane parameters - tile selection and number of points

//			dst.num_points = src.num_points;
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

				MeasuredLayersFilterParameters mlfp,
//				int          smplSide, //        = 2;      // Sample size (side of a square)
//				int          smplNum, //         = 3;      // Number after removing worst
//				double       smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				boolean      smplWnd,  // use window functions for the samples

//	  			double       max_abs_tilt,  //  2.0;   // pix per tile
//				double       max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//				double       damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				double       min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				double       transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				int          far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				double       far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

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
						disp_strength[nl] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
								nl, // int num_layer,
								getSTileXY()[0],        // int stX,
								getSTileXY()[1],        // int stY,
								(single_plane ? null : measuredSelection[nl]),  // boolean [] sel_in,
								mlfp,
//								strength_floor,
//								measured_strength_pow, //
//								smplSide, //        = 2;      // Sample size (side of a square)
//								smplNum, //         = 3;      // Number after removing worst
//								smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//								smplWnd,  // use window functions for the samples
//								max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//								max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//								damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//								min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//								transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//								far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//								far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
								true, // boolean null_if_none)
								debugLevel);
					} else {
						disp_strength[nl] = measuredLayers.getDisparityStrengthML(
								nl,                     // int num_layer,
								getSTileXY()[0],        // int stX,
								getSTileXY()[1],        // int stY,
								(single_plane ? null : measuredSelection[nl]),  // boolean [] sel_in,
								mlfp,
//								strength_floor,         //  double strength_floor,
//								measured_strength_pow,  // double strength_pow,
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
							disp_strength[nl] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
									nl, // int num_layer,
									getSTileXY()[0],        // int stX,
									getSTileXY()[1],        // int stY,
									null,  // boolean [] sel_in,
									mlfp,
//									strength_floor,
//									measured_strength_pow, //
//									smplSide, //        = 2;      // Sample size (side of a square)
//									smplNum, //         = 3;      // Number after removing worst
//									smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//									smplWnd,  // use window functions for the samples
//									max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//									max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//									damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//									min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//									transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//									far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//									far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
									true, // boolean null_if_none)
									debugLevel);
						} else {
							disp_strength[nl] = measuredLayers.getDisparityStrengthML(
									nl,                     // int num_layer,
									getSTileXY()[0],        // int stX,
									getSTileXY()[1],        // int stY,
									null,  // boolean [] sel_in,
									mlfp,
//									strength_floor,         //  double strength_floor,
//									measured_strength_pow,  // double strength_pow,
									true);                  // boolean null_if_none);
						}
						//disp_strength[nl] = measuredLayers.getDisparityStrengthML(
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
		 * This method is used to filter sample plates by averaging subset of the square
		 * sample and removing outliers. Currently only constant disparity and horizontal
		 * surfaces are used, this method is used for horizontal ones to find tilts
		 * d_disp/d_tx, d_disp/d_ty measured in pixels per tile
		 * @param world_normal_xyz world space normal vector, currently only "up" - (0,1,0) is used
		 * @param tile_disp_strengths [0] - unfiltered disparities for the supertile,
		 *                            [1] - unfiltered correlation strengths for the
		 *                                  supertile (just 0/non-0)
		 * @param debugLevel
		 * @return per tile arrays of either nulls or tilt-x, tilt-y pairs
		 */

		public double [][] getDisparityTilts(
				double []     world_normal_xyz,
				double [][]   tile_disp_strengths,
				int           debugLevel)
		{
			final int stSize2 = 2* stSize;
			final int stCenter = (stSize2 + 1) * superTileSize / 2;
			final double [][] tile_tilts = new double [stSize2*stSize2][];
			final double [] apy_vector = {0.0, 0.0, 1.0};
			final Matrix py_vector = new Matrix(apy_vector,3);
			invalidateCalculated();
			double [] px_py = getCenterPxPy();
			Matrix normal_row = new Matrix(world_normal_xyz,1); // 1x3
			Matrix normal_col = new Matrix(world_normal_xyz,3); // 3x1

			// find world coordinates of the center of tile intersection with the plane
			// disparity 0.0 is invalid, but it somehow got here
			for (int lTile = 0; lTile < tile_disp_strengths[1].length; lTile++) if ((tile_disp_strengths[1][lTile] > 0.0) && (tile_disp_strengths[0][lTile] > 0.0)){
				int tY = lTile / stSize2;
				int tX = lTile % stSize2;
				double px = px_py[0] + tX - stCenter;
				double py = px_py[1] + tY - stCenter;
				double disp = tile_disp_strengths[0][lTile];
				// get world coordinates for each tile as determined individually
				Matrix t_xyz = new Matrix(geometryCorrection.getWorldCoordinates(
						px,
						px,
						disp,
						this.correctDistortions),3);
//				double n_by_w =  normal_row.times(t_xyz).get(0, 0);
				// take pixel vector parallel to py and convert it to the world space
				Matrix jacobian =  new Matrix(geometryCorrection.getWorldJacobian(
						px,
						py,
						disp,
						this.correctDistortions,
						(debugLevel > 2)
						));
				Matrix inv_jacobian = jacobian.inverse();
				Matrix wpy_vector = jacobian.times(py_vector); // 3 rows, 1 column py vector in the real world
				// get 2 in-plane vectors in the real world
				Matrix  wv1 = cross3d(normal_col, wpy_vector);
				Matrix  wv2 = cross3d(normal_col, wv1);

				// convert then to pixel space
				Matrix pv1 = inv_jacobian.times(wv1);
				Matrix pv2 = inv_jacobian.times(wv2);

				// Get plane normal in the pixel space
				Matrix pn = cross3d(pv1, pv2);

				// convert to the two tilts (d/dpx, d/dpy). We will need in pix/tile, not pix/pix
				double [] txty = {
						-pn.get(1, 0)/pn.get(0, 0)*stSize,
						-pn.get(2, 0)/pn.get(0, 0)*stSize};

				if (Double.isNaN(txty[0]) || Double.isNaN(txty[1])){
					System.out.println("**** this is a BUG in getDisparityTilts() ****");
					System.out.println("txty= {"+txty[0]+","+txty[1]+"}");
					jacobian.print(10, 5);
					txty = null;
				}
				tile_tilts[lTile] = txty; // some are nulls?

			}

//			if (debugLevel > 1) {
//				System.out.println("st_xyz = {"+st_xyz.get(0, 0)+","+st_xyz.get(1, 0)+","+st_xyz.get(2, 0)+"}"+" ="+n_by_w);
//			}
			return tile_tilts;

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
			double [] px_py = getCenterPxPy(); // tile center
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
								disp_str[nl] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
										nl, // int num_layer,
										sTileXY[0], // int stX,
										sTileXY[1], // int stY,
										null, // measuredSelection[nl], // boolean [] sel_in,
										this.mlfp,
//										strength_floor,
//										measured_strength_pow, //
//										smplSide, //        = 2;      // Sample size (side of a square)
//										smplNum,  //         = 3;      // Number after removing worst
//										smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//										smplWnd,  // use window functions for the samples

//										max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//										max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//										damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//										min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//										transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//										far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//										far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

										true,     // boolean null_if_none)
										debugLevel);
							}
						} else {
							disp_str[nl] = measuredLayers.getDisparityStrengthML(
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									null, // measuredSelection[nl], // boolean [] sel_in,
									this.mlfp,
//									strength_floor, //
//									measured_strength_pow, //
									true); // boolean null_if_none)
						}
					}
				}
			}
			int numRemoved = 0;
			boolean no_bugs = true;
//			for (; (getValue() > targetEigen) && (numRemoved < maxRemoved); numRemoved++){
			for (; (getNormValue() > targetEigen) && (numRemoved < maxRemoved); numRemoved++){
				if (debugLevel > 2){
					System.out.println("removePlaneOutliers("+sTileXY[0]+":"+sTileXY[1]+"): numRemoved = "+numRemoved+
							" eigenValue = " + getValue()+" norm eigenValue = " + getNormValue()+" target = "+targetEigen);
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
//						strength_floor,
//						measured_strength_pow, // double       strength_pow,
						this.smplMode,
						this.mlfp,

//						smplSide,
//						smplNum,
//						smplRms,
//						smplWnd,  // use window functions for the samples
//						max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//						max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//						damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//						min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//						transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//						far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//						far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
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
//							strength_floor,
//							measured_strength_pow, // double       strength_pow,
							this.smplMode,
							this.mlfp,
//							smplSide,
//							smplNum,
//							smplRms,
//							smplWnd,  // use window functions for the samples
//							max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//							max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//							damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//							min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//							transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//							far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//							far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
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
		 * @param smplWnd  use window functions for the samples
		 *
		 * @param debugLevel debug level
		 * @return per measurement layer : x,y,z, weight, or null if failed. This
		 * value may be re-used in subsequent refinements (as removing outliers)
		 */
		public double [][][] getPlaneFromMeas(
				boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
				double [][][] disp_str, // calculate just once when removing outliers
				double       disp_far, // minimal disparity to select (or NaN)
				double       disp_near, // maximal disparity to select (or NaN)
				double       dispNorm,   //  Normalize disparities to the average if above
				double       min_weight,
				int          min_tiles,
//				double       strength_floor,
//				double       strength_pow,
				boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				MeasuredLayersFilterParameters mlfp,

//				int          smplSide, //        = 2;      // Sample size (side of a square)
//				int          smplNum, //         = 3;      // Number after removing worst
//				double       smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				boolean      smplWnd,        // use window functions for the samples

//	  			double       max_abs_tilt,  //  2.0;   // pix per tile
//				double       max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//				double       damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				double       min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				double       transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				int          far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				double       far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

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

//			this.strength_floor =         strength_floor;
//			this.measured_strength_pow =  strength_pow;
			this.min_weight =             min_weight;
			this.min_tiles =              min_tiles;
			this.dispNorm =               dispNorm;
//			this.smplWnd =                smplWnd;        // use window functions for the samples
			this.smplMode =               smplMode; //        = true;   // Use sample mode (false - regular tile mode)
//			this.smplSide =               smplSide; //        = 2;      // Sample size (side of a square)
//			this.smplNum =                smplNum;    //         = 3;      // Number after removing worst
//			this.smplRms =                smplRms;    //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

//			this.max_abs_tilt =  max_abs_tilt;
//			this.max_rel_tilt =  max_rel_tilt;
//			this.damp_tilt =     damp_tilt;
//			this.min_tilt_disp = min_tilt_disp;
//			this.transition =    transition;
//			this.far_mode =      far_mode;
//			this.far_power =     far_power;

			this.mlfp =                   mlfp.clone();


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
							disp_str[nl] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									//tile_sel[nl], // boolean [] sel_in,
									mlfp,
//									strength_floor,
//									strength_pow, //
//									smplSide, //        = 2;      // Sample size (side of a square)
//									smplNum, //         = 3;      // Number after removing worst
//									smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//									smplWnd, // final boolean    smplWnd,        // use window functions fro the samples
//									max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//									max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//									damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//									min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//									transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//									far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//									far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
									true, // boolean null_if_none)
									debugLevel);
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
									mlfp.strength_floor,
									true); // boolean null_if_none)
						} else {
							tile_sel[nl] =  measuredLayers.getSupertileSelection(
									nl,            // int num_layer,
									sTileXY[0],    // int stX,
									sTileXY[1],    // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									disp_far,      // 		double     disp_far,
									disp_near,     // double     disp_near,
									mlfp.strength_floor,
									true);         // boolean null_if_none)

						}
						num_tiles += MeasuredLayers.getNumSelected(tile_sel[nl]);
						if (tile_sel[nl] != null){
							disp_str[nl] = measuredLayers.getDisparityStrengthML(
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									tile_sel[nl], // boolean [] sel_in,
									mlfp,
//									strength_floor,
//									strength_pow, //
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
						  if (tile_sel[nl] != null) {
							  for (int i = 0; i < dbg_img[2].length; i++){
								  dbg_img[2][i] = tile_sel[nl][i]?1.0:0.0; // exception here?
							  }
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
//								double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5; // in pixels, not in tiles
//								double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5;
								double x = ((indx % stSize2) - stSize + 0.5) * tileSize; // in pixels, not in tiles
								double y = ((indx / stSize2) - stSize + 0.5) * tileSize;
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

//			double kz = ((dispNorm > 0.0) && (swz > dispNorm)) ? (dispNorm / swz) : 1.0;

			if (debugLevel > 0){
				System.out.println("getPlaneFromMeas(): num_tiles="+num_tiles+", sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy); // +", kz="+kz);
			}

			// TODO: scale disparity to make same scale for 3 axes?

			for (int nl = 0; nl < tile_sel.length; nl++){
				if (disp_str[nl] != null) {
					for (int indx = 0; indx < disp_str[nl][0].length; indx++){
						if (tile_sel[nl][indx]) {
							double w = disp_str[nl][1][indx] / sw;
							if (w > 0.0){
//								double d =  kz * (disp_str[nl][0][indx] - swz); // Not here!
								double d =  disp_str[nl][0][indx] - swz;
								double wd = w*d;
//								double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5 - swx;
//								double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5 - swy;
								double x = ((indx % stSize2) - stSize + 0.5) * tileSize  - swx;
								double y = ((indx / stSize2) - stSize + 0.5) * tileSize  - swy;
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


		// similar to getPlaneFromMeas, but building ellipsoids in the real world space
		public double [][][] getWorldPlaneFromMeas(
				boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
				// TODO: mame it accept tiles_xyzw (same as output)
				double [][][] disp_str, // calculate just once when removing outliers
				double       disp_far, // minimal disparity to select (or NaN)
				double       disp_near, // maximal disparity to select (or NaN)
				double       dispNorm,   //  Normalize disparities to the average if above
				double       min_weight,
				int          min_tiles,

//				double       strength_floor,
//				double       strength_pow,

				boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				MeasuredLayersFilterParameters mlfp,

//				int          smplSide, //        = 2;      // Sample size (side of a square)
//				int          smplNum, //         = 3;      // Number after removing worst
//				double       smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

///* was not here */boolean      smplWnd,        // use window functions for the samples

//	  			double       max_abs_tilt,  //  2.0;   // pix per tile
//				double       max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//				double       damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				double       min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				double       transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				int          far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				double       far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

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

//			this.strength_floor =         strength_floor;
//			this.measured_strength_pow =  strength_pow;
			this.min_weight =             min_weight;
			this.min_tiles =              min_tiles;
			this.dispNorm =               dispNorm;
			this.smplMode =               smplMode; //        = true;   // Use sample mode (false - regular tile mode)
//			this.smplSide =               smplSide; //        = 2;      // Sample size (side of a square)
//			this.smplNum =                smplNum;    //         = 3;      // Number after removing worst
//			this.smplRms =                smplRms;    //         = 0.1;    // Maximal RMS of the remaining tiles in a sample

//			this.smplWnd =                smplWnd;    //        was not here !

//			this.max_abs_tilt =           max_abs_tilt;
//			this.max_rel_tilt =           max_rel_tilt;
//			this.damp_tilt =              damp_tilt;
//			this.min_tilt_disp =          min_tilt_disp;
//			this.transition =             transition;
//			this.far_mode =               far_mode;
//			this.far_power =              far_power;
			this.mlfp = mlfp.clone();

			if (debugLevel > 2){
				System.out.println("getWorldPlaneFromMeas()");
			}
			boolean need_disp_str = false;
			if (disp_str == null) {
				disp_str =      new double [tile_sel.length][][];
				need_disp_str = true;
			}
			// TODO: Code duplication with getWorldPlaneFromMeas() - extract common part
			for (int nl = 0; nl < tile_sel.length; nl++){
				if (tile_sel[nl] != null){
					if (smplMode) {
						if (need_disp_str) {
							disp_str[nl] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									//tile_sel[nl], // boolean [] sel_in,
									mlfp,
//									strength_floor,
//									strength_pow, //
//									smplSide, //        = 2;      // Sample size (side of a square)
//									smplNum,  //         = 3;      // Number after removing worst
//									smplRms,  //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
///*was using this. */				smplWnd,  // use window functions for the samples
//									max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//									max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//									damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//									min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//									transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//									far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//									far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
									true,     // boolean null_if_none)
									debugLevel);
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
									mlfp.strength_floor,
									true); // boolean null_if_none)
						} else {
							tile_sel[nl] =  measuredLayers.getSupertileSelection(
									nl,            // int num_layer,
									sTileXY[0],    // int stX,
									sTileXY[1],    // int stY,
									((tile_sel[nl].length == 0)? null:tile_sel[nl]), // boolean [] sel_in,
									disp_far,      // 		double     disp_far,
									disp_near,     // double     disp_near,
									mlfp.strength_floor,
									true);         // boolean null_if_none)

						}
						num_tiles += MeasuredLayers.getNumSelected(tile_sel[nl]);
						if (tile_sel[nl] != null){
							disp_str[nl] = measuredLayers.getDisparityStrengthML(
									nl, // int num_layer,
									sTileXY[0], // int stX,
									sTileXY[1], // int stY,
									tile_sel[nl], // boolean [] sel_in,
									mlfp,
//									strength_floor,
//									strength_pow, //
									true); // boolean null_if_none)
							sw += MeasuredLayers.getSumStrength(disp_str[nl]);
						}
					}


					if ((debugLevel > 2) && (disp_str[nl] != null)){
//					if ((debugLevel > 1) && (disp_str[nl] != null)){
						  showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
						  double [][] dbg_img = new double [3][];
						  dbg_img[0] = disp_str[nl][0];
						  dbg_img[1] = disp_str[nl][1];
						  dbg_img[2] = new double [stSize2*stSize2];
						  if (tile_sel[nl] != null) {
							  for (int i = 0; i < dbg_img[2].length; i++){
								  dbg_img[2][i] = tile_sel[nl][i]?1.0:0.0;
							  }
						  }
						  sdfa_instance.showArrays(dbg_img,  stSize2, stSize2, true, "disp_str_x"+sTileXY[0]+"_y"+sTileXY[1]+"_"+nl);
					}
				}
			}
			this.measuredSelection =      tile_sel; // it may be modified


			if ((sw < min_weight) || (num_tiles < min_tiles)) {
				if (debugLevel > 1){
					System.out.println("getWorldPlaneFromMeas():return false");
				}
				return null; // too weak plane or too few tiles selected
			}

			double [][][] tiles_xyzw = new double [disp_str.length][][];
			for (int nl = 0; nl < tile_sel.length; nl++) if (disp_str[nl] != null) {
				tiles_xyzw[nl] = new double [disp_str[nl][0].length][];
			}
			double [][] acovar = new double [3][3];
			double swz = 0.0, swx = 0.0, swy = 0.0;
			sw =0.0;
			double [] px_py = getCenterPxPy();

			for (int nl = 0; nl < tile_sel.length; nl++){
				if (disp_str[nl] != null) {
					for (int indx = 0; indx < disp_str[nl][0].length; indx++){
						if (tile_sel[nl][indx]) {
							double w = disp_str[nl][1][indx];
							if (w > 0.0){
								tiles_xyzw[nl][indx] = new double [4];
								double d = disp_str[nl][0][indx];
								// referencing samples to centers of pixels
//								double x = ((indx % stSize2) - stSize + 0.5) * tileSize + 0.5 + px_py[0]; // in pixels, not in tiles
//								double y = ((indx / stSize2) - stSize + 0.5) * tileSize + 0.5 + px_py[1];
								double x = ((indx % stSize2) - stSize + 0.5) * tileSize + px_py[0]; // in pixels, not in tiles
								double y = ((indx / stSize2) - stSize + 0.5) * tileSize + px_py[1];
								// difference from getPlaneFromMeas
								double [] wxyz = geometryCorrection.getWorldCoordinates(
										x,
										y,
										d,
										this.correctDistortions);

								tiles_xyzw[nl][indx][0] = wxyz[0];
								tiles_xyzw[nl][indx][1] = wxyz[1];
								tiles_xyzw[nl][indx][2] = wxyz[2];
								tiles_xyzw[nl][indx][3] = w;
								sw  += w;
								swz += w * wxyz[2];
								swx += w * wxyz[0];
								swy += w * wxyz[1];
								if (Double.isNaN(tiles_xyzw[nl][indx][0])) {
									System.out.println("--*--BUG! tiles_xyzw[nl][indx][0] is NaN");
								}
								if (Double.isInfinite(swx) || Double.isInfinite(swz)  || Double.isInfinite(w)) {
									System.out.println("BUG!!!: getPlaneFromMeas(): num_tiles="+num_tiles+", sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy);
								}


								// end of difference from getPlaneFromMeas
							}
						}
					}
					if ((debugLevel > 3) && (disp_str[nl] != null)){
						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
						double [][] dbg_img = new double [4][stSize2*stSize2];
						for (int indx = 0; indx < dbg_img[0].length; indx++){
							if (tiles_xyzw[nl][indx] != null) {
							dbg_img[0][indx] = tiles_xyzw[nl][indx][0];
							dbg_img[1][indx] = tiles_xyzw[nl][indx][1];
							dbg_img[2][indx] = tiles_xyzw[nl][indx][2];
							dbg_img[3][indx] = tiles_xyzw[nl][indx][3];
							} else {
								dbg_img[0][indx] = Double.NaN;
								dbg_img[1][indx] = Double.NaN;
								dbg_img[2][indx] = Double.NaN;
								dbg_img[3][indx] = Double.NaN;
							}
						}
						sdfa_instance.showArrays(dbg_img,  stSize2, stSize2, true, "world_x"+sTileXY[0]+"_y"+sTileXY[1]+"_"+nl);
					}
				}
			}
			if (sw == 0.0) {
				return null; //
			}
			if (Double.isInfinite(swx)) {
				System.out.println("BUG!!!: getPlaneFromMeas(): num_tiles="+num_tiles+", sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy);
			}
			swz /= sw;
			swx /= sw;
			swy /= sw;
			setWxyz(swx, swy, swz);
			if (Double.isInfinite(swx)) {
				System.out.println("BUG!!!: getPlaneFromMeas(): num_tiles="+num_tiles+", sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy);
			}

//			double kz = ((dispNorm > 0.0) && (swz > dispNorm)) ? (dispNorm / swz) : 1.0;

			if (debugLevel > 0){
				System.out.println("getPlaneFromMeas(): num_tiles="+num_tiles+", sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy);
			}

			// TODO: scale disparity to make same scale for 3 axes?

			for (int nl = 0; nl < tile_sel.length; nl++){
				if (disp_str[nl] != null) {
					for (int indx = 0; indx < disp_str[nl][0].length; indx++){
						if (tiles_xyzw[nl][indx] != null) {
							double w = tiles_xyzw[nl][indx][3] / sw;
							if (w > 0.0){
								double x = tiles_xyzw[nl][indx][0] - swx;
								double y = tiles_xyzw[nl][indx][1] - swy;
								double z = tiles_xyzw[nl][indx][2] - swz;
								acovar [0][0] += w * x * x;
								acovar [0][1] += w * x * y;
								acovar [0][2] += w * x * z;
								acovar [1][1] += w * y * y;
								acovar [1][2] += w * y * z;
								acovar [2][2] += w * z * z;
								if (Double.isNaN(acovar [0][0])) {
									System.out.println("--*--BUG! acovar[0][0] is NaN");
								}
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


			// probably this reordering is not needed as they are already ordered
			// for world coordinates the sequence is normal x,y,z (not d,x,y)

			int oindx = 0;
			for (int i = 1; i < 3 ; i++){
				if (eig_val[i][i] < eig_val[oindx][oindx]){
					oindx = i;
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
			double from = 0.0; // see if the first vector ("plane" normal) is away from the camera
			for (int i = 0; i < 3; i++) {
				from += this.wxyz[i]*plane[0][i];
			}
			if (from > 0.0){ // reverse the first vector to be towards the camera
				for (int i = 0; i < 3; i++) {
					plane[0][i] = -plane[0][i];
				}
			}

			// make  direction last vector so px (x) py (.) disp < 0 (left-hand coordinate system)
			if (new Matrix(plane).det() > 0){
				for (int i = 0; i < 3; i ++) plane[2][i] = -plane[2][i];
			}

//			setZxy(swz, swx, swy);
			setWeight(sw); // should be the same

			setWValues(eig_val[oindx][oindx],eig_val[vindx][vindx],eig_val[hindx][hindx]); // eigenvalues [0] - thickness, 2 other to detect skinny (poles)
			setWVectors   (plane);
			setNumPoints (num_tiles);  // should be the same
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
//			return disp_str;
			return disp_str;

		}



		public double [][] initMergedValue()
		{
			this.merged_eig_val =      new double[8][];
			this.merged_eig_eq =       new double[8][];
			this.merged_weig_val =     new double[8][];
			this.merged_weig_eq =      new double[8][];
			this.merged_valid =        new boolean[8][];
			this.merged_strong_valid = new boolean[8][];
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

		public double [][] getMergedWValue()
		{
			return this.merged_weig_val;
		}

		public double [][] getMergedWValueEq()
		{
			return this.merged_weig_eq;
		}
		public double [][] getLinkCosts()
		{
			return this.link_costs;
		}
		public double [] getLinkCosts(int dir)
		{
			if ((link_costs==null) || (link_costs[dir] == null)){
				return null;
			}
			return this.link_costs[dir];
		}
		public double getLinkCosts(int dir, int np)
		{
			if ((link_costs==null) || (link_costs[dir] == null)){
				return Double.MAX_VALUE;
			}
			return this.link_costs[dir][np];
		}

		public void setLinkCosts(int dir, int plane, double value)
		{
			this.link_costs[dir][plane] = value;
		}

		public double [] initMergedValue(int dir, int leng)
		{
			this.merged_eig_val[dir] =        new double[leng];
			this.merged_eig_eq[dir] =         new double[leng];
			this.merged_weig_val[dir] =       new double[leng];
			this.merged_weig_eq[dir] =        new double[leng];
			this.merged_valid[dir] =          new boolean[leng];
			this.merged_strong_valid[dir] =   new boolean[leng];
			for (int i = 0; i < leng; i++) {
				this.merged_eig_val[dir][i] =  Double.NaN;
				this.merged_eig_eq[dir][i] =   Double.NaN;
				this.merged_weig_val[dir][i] = Double.NaN;
				this.merged_weig_eq[dir][i] =  Double.NaN;
			}
			return getMergedValue(dir);
		}

		public double [][] initLinkCosts(){
			this.link_costs =          new double[8][];
			return link_costs;
		}
		public double [] initLinkCosts(int dir, int leng)
		{
			this.link_costs[dir] =            new double[leng];
			for (int i = 0; i < leng; i++) {
				this.link_costs[dir][i] =      Double.NaN;
			}
			return getLinkCosts(dir);
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

		public double [] getMergedWValue(int dir)
		{
			if (this.merged_weig_val == null) {
				return null;
			}
			return this.merged_weig_val[dir];
		}

		public double [] getMergedWValueEq(int dir)
		{
			if (this.merged_weig_eq == null) {
				return null;
			}
			return this.merged_weig_eq[dir];
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

		public double getMergedWValue(int dir, int plane)
		{
			if ((this.merged_weig_val == null) || (this.merged_weig_val[dir] == null)){
				return Double.NaN;
			}
			return this.merged_weig_val[dir][plane];
		}

		public double getMergedWValueEq(int dir, int plane)
		{
			if ((this.merged_weig_eq == null) ||(this.merged_weig_eq[dir] == null)){
				return Double.NaN;
			}
			return this.merged_weig_eq[dir][plane];
		}

		public void setNeibMatch(int dir, int plane, double value)
		{
			this.merged_eig_val[dir][plane] = value;
		}
		public void setNeibMatchEq(int dir, int plane, double value)
		{
			this.merged_eig_eq[dir][plane] = value;
		}

		public void setNeibWMatch(int dir, int plane, double value)
		{
			this.merged_weig_val[dir][plane] = value;
		}
		public void setNeibWMatchEq(int dir, int plane, double value)
		{
			this.merged_weig_eq[dir][plane] = value;
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

		public boolean hasMergedValid(int dir){
			if ((this.merged_valid == null) || (this.merged_valid[dir] == null)){
				return false;
			}
			for (int np = 0; np < this.merged_valid[dir].length; np++){
				if (this.merged_valid[dir][np]) return true;
			}
			return false;

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


		public boolean [][] getMergedStrongValid()
		{
			return this.merged_strong_valid;
		}

		public boolean [] getMergedStrongValid(int dir)
		{
			if (this.merged_strong_valid == null) {
				return null;
			}
			return this.merged_strong_valid[dir];
		}

		public boolean hasMergedStrongValid(int dir){
			if ((this.merged_strong_valid == null) || (this.merged_strong_valid[dir] == null)){
				return false;
			}
			for (int np = 0; np < this.merged_strong_valid[dir].length; np++){
				if (this.merged_strong_valid[dir][np]) return true;
			}
			return false;

		}
		public boolean isMergedStrongValid(int dir, int plane)
		{
			if ((this.merged_strong_valid == null) || (this.merged_strong_valid[dir] == null)){
				return false;
			}
			return this.merged_strong_valid[dir][plane];
		}

		public void setMergedStrongValid(int dir, int plane, boolean valid)
		{
			this.merged_strong_valid[dir][plane] = valid;
		}

		public void setMergedStrongValid(int dir, int plane, boolean valid, int leng)
		{
			if (this.merged_strong_valid == null){
				this.merged_strong_valid = new boolean[8][];
			}
			if (this.merged_strong_valid[dir] == null){
				this.merged_strong_valid[dir] = new boolean[leng];
			}
			this.merged_strong_valid[dir][plane] = valid;
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

		/**
		 * Return "normalized" main eigenvalue - it is reduced for high disparities
		 * @return normalized main eigenvalue
		 */
		public double getNormValue(){
			double val = values[0];
			if ((dispNorm > 0.0) && (zxy[0] > dispNorm)){
				double k = dispNorm / zxy[0];
//				val *= k*k;
				val *= k ; // reducing correction for the large disparities (making it sqrt() )
			}
			return val;
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


		public double[] getWxyz() {
			return wxyz;
		}


		public void setWxyz(double[] wxyz) {
			this.wxyz = wxyz;
		}
		public void setWxyz(
				double x,
				double y,
				double z) {
			this.wxyz = new double [3];
			this.wxyz[0] = x;
			this.wxyz[1] = y;
			this.wxyz[2] = z;
		}




		public double[][] getWVectors() {
			return wvectors;
		}
		public double[] getWVector() {
			return wvectors[0];
		}
		public void setWVectors(double[][] wvectors) {
			this.wvectors = wvectors;
		}
		public double[] getWValues() {
			return wvalues;
		}
		public double getWValue() {
			return wvalues[0];
		}
		public void setWValues(double[] wvalues) {
			this.wvalues = wvalues;
		}
		public void setWValues(double wv1, double wv2, double wv3) {
			this.wvalues = new double[3];
			this.wvalues[0] = wv1;
			this.wvalues[1] = wv2;
			this.wvalues[2] = wv3;
		}

		public int getNumPoints() {
			return num_points;
		}
		public void setNumPoints(int num_points) {
			if (num_points < 0){
				System.out.println("setNumPoints(): Setting negative number of tiles in a plane: "+num_points);
			}
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
//				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				double y = tileSize * (sy + 0.5)  - zxy[2];
				for (int sx = -superTileSize/2; sx < superTileSize/2; sx++){
//					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					double x = tileSize * (sx + 0.5) - zxy[1];
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
//				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				double y = tileSize * (sy + 0.5)   - zxy[2];
				for (int sx = -superTileSize; sx < superTileSize; sx++){
//					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					double x = tileSize * (sx + 0.5) - zxy[1];
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
//				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				double y = tileSize * (sy + 0.5)  - zxy[2];
				for (int sx = -3 * superTileSize/2; sx < 3 * superTileSize / 2; sx++){
//					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					double x = tileSize * (sx + 0.5) - zxy[1];
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

		public double[][] getDoublePlaneDisparityStrengthOld(
				double [] window,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				int       debugLevel)
		{		return getDoublePlaneDisparityStrengthOld(
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
		public double[][] getDoublePlaneDisparityStrengthOld(
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
//				double y = tileSize * (sy + 0.5) + 0.5 - zxy[2];
				double y = tileSize * (sy + 0.5) - zxy[2];
				for (int sx = -superTileSize; sx < superTileSize; sx++){
//					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					double x = tileSize * (sx + 0.5) - zxy[1];
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
			return  getDoublePlaneDisparityStrength(
					useWorld ? getWorldXYZ(this.correctDistortions) : null, // will calculate if not yet done so. Should it use otherPd, not pd? and then clone later?
					window,
					dir,
					use_sel,
					divide_by_area,
					scale_projection,
					fraction_uni,
					debugLevel);
		}

		public double [] getWorldXYZFromWorld()
		{
			double []   wxyz = getWxyz();
			double [][] wvectors = getWVectors();
			if ((wxyz == null) || (wvectors == null)  || (wvectors[0] == null)) return null;
			double normal_dot_normal = 0.0; // supposed to be unit length, but we'll still re-normalize it
			double center_dot_normal = 0.0;
			for (int i = 0; i < 3; i++){
				normal_dot_normal += wvectors[0][i] *  wvectors[0][i];
				center_dot_normal += wxyz[i] *  wvectors[0][i];
			}
			double scale = center_dot_normal / normal_dot_normal;
			double [] rslt = wvectors[0].clone();
			for (int i = 0; i < 3; i++){
				rslt[i] *= scale;
			}
			return rslt;
		}

		public double[][] getDoublePlaneWorldDisparityStrength(
				double [] window,
				int       dir,
				boolean   use_sel,
				boolean   divide_by_area,
				double    scale_projection,
				double    fraction_uni,
				int       debugLevel)
		{

			return  getDoublePlaneDisparityStrength(
					getWorldXYZFromWorld(),
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
		 * @param world_normal Noral vector to the plane from the origin, or null to stay in pixel/disparity space.
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
				double [] world_normal, // either real world normal vector from (0,0,0) to the plane or null, in that case will stay n disparity space
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
//				double y = tileSize * (iy - ss2 + 0.5) + 0.5  - zxy[2];
				double y = tileSize * (iy - ss2 + 0.5)  - zxy[2];
				int oy = iy + offsets[dir1][4]; //vert index in the result tile
				for (int ix = offsets[dir1][2]; ix < offsets[dir1][3]; ix++){
//					double x = tileSize * (ix - ss2 + 0.5) + 0.5 - zxy[1];
					double x = tileSize * (ix - ss2 + 0.5) - zxy[1];
					int indx = ss4 * oy + ix + offsets[dir1][5];
					int indx_i = iy * ss4 + ix; // input index
					if (world_normal != null) {
						disp_strength[0][indx] = geometryCorrection.getPlaneDisparity( // disparity (at this center) for crossing other supertile plane
								world_normal,
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
//				double y = tileSize * (iy - ss2 + 0.5) + 0.5  - zxy[2];
				double y = tileSize * (iy - ss2 + 0.5)  - zxy[2];
				int oy = iy + offsets[dir1][4]; //vert index in the result tile
				for (int ix = offsets[dir1][2]; ix < offsets[dir1][3]; ix++){
//					double x = tileSize * (ix - ss2 + 0.5) + 0.5 - zxy[1];
					double x = tileSize * (ix - ss2 + 0.5) - zxy[1];
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
			if (values == null ) {
				System.out.println("mergePlaneToThis(): values=null:\n"+toString());
				return null;

			}
			if (otherPd.values == null ) {
				System.out.println("mergePlaneToThis(): otherPd.values=null:\n"+otherPd.toString());
				return null;
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
			if (debugLevel > 1) {
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
			if (debugLevel > 1) {
				System.out.println("other_covar");
				other_covar.print(8, 6);
				System.out.println("this_covar");
				this_covar.print(8, 6);
				System.out.println("covar");
				covar.print(8, 6);
			}

			covar.plusEquals(other_covar.times(other_fraction));
			if (debugLevel > 1) {
				System.out.println("covar with other_covar");
				covar.print(8, 6);
			}
			covar.plusEquals(this_covar.times(1.0 - other_fraction));
			if (debugLevel > 1) {
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
			if (debugLevel > 1) {
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
			pd.setNumPoints(otherPd.getNumPoints()+this.getNumPoints());
			// Repeat merging for world-based planes
			return mergePlaneToThisWorld(
					otherPd, // PlaneData otherPd,
					pd,      // PlaneData pd_partial, // disparity-based data is already merged
					scale_other,
					starWeightPwr,    // Use this power of tile weight when calculating connection cost
					ignore_weights,
					sum_weights,
					preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
					debugLevel);
//			return pd;
		}


		// only handles wxyz, wvalues,wvectors - the rest should be done with mergePlaneToThis()
		private PlaneData mergePlaneToThisWorld(
				PlaneData otherPd,
				PlaneData pd_partial, // disparity-based data is already merged
				double    scale_other,
				double    starWeightPwr,    // Use this power of tile weight when calculating connection cost
				boolean   ignore_weights,
				boolean   sum_weights,
				boolean   preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
				int       debugLevel)
		{
			if (debugLevel > 0) {
				System.out.println("mergePlaneToThisWorld()");
			}
			if (wvalues == null ) {
				if (debugLevel > -2) {
					System.out.println("mergePlaneToThisWorld(): wvalues=null:\n"+toString());
				}
				return null;

			}
			if (otherPd.wvalues == null ) {
				if (debugLevel > -1) System.out.println("mergePlaneToThisWorld(): otherPd.wvalues=null:\n"+otherPd.toString());
				return null;
			}
			double [][] this_eig_avals = {
					{wvalues[0], 0.0,        0.0},
					{0.0,        wvalues[1], 0.0},
					{0.0,        0.0,        wvalues[2]}};
			double [][] other_eig_avals = {
					{otherPd.wvalues[0], 0.0,                0.0},
					{0.0,                otherPd.wvalues[1], 0.0},
					{0.0,                0.0,                otherPd.wvalues[2]}};
			Matrix this_eig_vals =     new Matrix(this_eig_avals);
			Matrix other_eig_vals =    new Matrix(other_eig_avals);
			Matrix other_eig_vectors = new Matrix(otherPd.wvectors).transpose(); // vectors are saved as rows
			Matrix this_eig_vectors =  new Matrix(this.wvectors).transpose();    // vectors are saved as rows
			Matrix this_center    =    new Matrix(this.getWxyz(),3);
			Matrix other_center    =   new Matrix(otherPd.getWxyz(),3); // should already be relative to this supertile center
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
			if (debugLevel > 1) {
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
			if (debugLevel > 1) {
				System.out.println("other_covar");
				other_covar.print(8, 6);
				System.out.println("this_covar");
				this_covar.print(8, 6);
				System.out.println("covar");
				covar.print(8, 6);
			}

			covar.plusEquals(other_covar.times(other_fraction));
			if (debugLevel > 1) {
				System.out.println("covar with other_covar");
				covar.print(8, 6);
			}
			covar.plusEquals(this_covar.times(1.0 - other_fraction));
			if (debugLevel > 1) {
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
			if (debugLevel > 1) {
				System.out.println("eig.getV()");
				eig.getV().print(8, 6);
				System.out.println("eig.getD()");
				eig.getD().print(8, 6);
			}



			double [][] eig_vect = eig.getV().getArray();
			double [][] eig_val = eig.getD().getArray();

			// probably this reordering is not needed as they are already ordered
			// for world coordinates the sequence is normal x,y,z (not d,x,y)
			int oindx = 0;
			for (int i = 1; i <3; i++){
				if (eig_val[i][i] < eig_val[oindx][oindx]){
					oindx = i;
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

			// make towards camera, left coordinate system


			// Make normal be towards camera (positive disparity), next vector - positive in X direction (right)
			double from = 0.0; // see if the first vector ("plane" normal) is away from the camera
			for (int i = 0; i < 3; i++) {
				from += this.wxyz[i]*plane[0][i];
			}
			if (from > 0.0){ // reverse the first vector to be towards the camera
				for (int i = 0; i < 3; i++) {
					plane[0][i] = -plane[0][i];
				}
			}

			// make  direction last vector so px (x) py (.) disp < 0 (left-hand coordinate system)
			if (new Matrix(plane).det() > 0){
				for (int i = 0; i < 3; i ++) plane[2][i] = -plane[2][i];
			}

			pd_partial.setWValues(eig_val[oindx][oindx],eig_val[vindx][vindx],eig_val[hindx][hindx]); // eigenvalues [0] - thickness, 2 other to detect skinny (poles)
			pd_partial.setWVectors(plane);

			pd_partial.setWxyz(common_center.getColumnPackedCopy()); // set new center
			// what weight to use? cloned is original weight for this supertile
			// or use weighted average like below?
			if (debugLevel < -1000) { // already done in mergePlaneToThis() that call this method
				double new_weight;
				if (sum_weights) {
					new_weight = sum_weight; // normalize while averaging by the caller
				} else { // how it was before
					new_weight = other_fraction * other_weight + (1.0 - other_fraction) * this_weight;
				}
				if (!ignore_weights && ((starWeightPwr != 1.0))){
					new_weight = Math.pow(new_weight,1.0/starWeightPwr);
				}
				pd_partial.setWeight(new_weight);
			}
			if (debugLevel > 2) {
				double L1 = getWValue();
				double L2 = otherPd.getWValue();
				double W1 = 1.0 - other_fraction;
				double W2 = other_fraction;
				double Lav = (L1*W1 + L2*W2)/(W1+W2);
				double L = pd_partial.getWValue();
				if ((L*1.000001 < Lav) || (debugLevel > 0)){
					System.out.println("========== mergePlaneToThisWorld(): L1="+L1+", L2 = "+L2+", W1="+W1+", W2="+W2+", Lav = "+Lav+", L="+L+
							", L-Lav="+(L-Lav)+", scale_other="+scale_other+", other_fraction="+other_fraction);
					System.out.println();
				}
			}

			return pd_partial;
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
			pd.num_points = otherPd.num_points; // restore, maybe remove from copy_neib?
// copy other data from this tile
			pd.setMeasSelection(this.measuredSelection);
			if (this.measuredLayers != null) pd.measuredLayers = this.measuredLayers;
			if (this.sel_mask != null)       pd.sel_mask = this.sel_mask.clone();
			copyStar(this,pd);
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

		//this.correctDistortions
		public double [] getCenterXYZ(
				int debugLevel){
			return getCenterXYZ(this.correctDistortions, debugLevel);
		}

		public double [] getCenterXYZ(
				boolean correct_distortions,
				int debugLevel)
		{
//			double delta = 0.0001;
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
				int debugLevel){
			return getWorldXYZ(this.correctDistortions, debugLevel);
		}
		/**
		 * Get vector from the camera perpendicular to the plane converted to the real world
		 * @param correct_distortions
		 * @param debugLevel
		 * @return
		 */

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

		public ArrayList<PlaneData> createTilePlanesFromSelections(
				String        suffix,
				boolean [][][] plane_selections, //  = new boolean [nStiles][][][]; // num_tiles
				double  [][][] disp_strength,
				double       dispNorm,   //  Normalize disparities to the average if above
				int          min_tiles,
				double       plTargetEigen, //        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
				double       plFractOutliers, //      =   0.3;  // Maximal fraction of outliers to remove
				int          plMaxOutliers, //        =    20;  // Maximal number of outliers to remove
//				double       strength_floor,
//				double       strength_pow,
				boolean      correct_distortions,
				boolean      smplMode, //        = true;   // Use sample mode (false - regular tile mode)

				MeasuredLayersFilterParameters mlfp,

//				int          smplSide, //        = 2;      // Sample size (side of a square)
//				int          smplNum, //         = 3;      // Number after removing worst
//				double       smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				boolean      smplWnd,  // use window functions for the samples

//	  			double       max_abs_tilt,  //  2.0;   // pix per tile
//				double       max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//				double       damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				double       min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				double       transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				int          far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				double       far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

				int          debugLevel)
		{
			if (debugLevel > 2) {
				debugLevel += 0; // +=1 // no show all eigen stuff (debugLevel > 3)
			}
			if (debugLevel > 0) {
				System.out.println("Debug debugLevel"); // +=1 // no show all eigen stuff (debugLevel > 3)
			}

			// first make a plane from all tiles
			ArrayList<PlaneData> st_planes = new ArrayList<PlaneData>();

			// iterate through all plane selections
			for (int ps = 0; ps < plane_selections.length; ps++) {
				PlaneData pd = this.clone();
				boolean OK = (pd.getPlaneFromMeas(
						plane_selections[ps], // tile_sel,       // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
						disp_strength,
						Double.NaN,    // double       disp_far, // minimal disparity to select (or NaN)
						Double.NaN,    // double       disp_near, // maximal disparity to select (or NaN)
						dispNorm,    // 0.0,            // plDispNorm,  // double       dispNorm,   //  Normalize disparities to the average if above
						0.0,            // double       min_weight,
						min_tiles,    // int          min_tiles,

//						strength_floor, //
//						strength_pow,   // double       strength_pow,

						// update !
						smplMode,
						mlfp,
//						smplSide,
//						smplNum,
//						smplRms,
//						smplWnd,  // use window functions for the samples

//						max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//						max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//						damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//						min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//						transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//						far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//						far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

						debugLevel) != null);            // int          debugLevel)
				if (OK) {
					if (debugLevel > 0) {
						if (pd.getWeight() > 1.0) {
							System.out.println("Processing subplane "+ suffix+
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

					/* Does it needs to be twice?
					double z0 = pd.getZxy()[0];
					if ((dispNorm > 0.0) && (z0 > dispNorm)) { // not needed ?
						double dd = (dispNorm + z0)/ dispNorm; // > 1
						targetV *= dd * dd; // > original
					}
					*/
//					if (pd.getValues()[0] > targetV) {
					if (pd.getNormValue() > targetV) {
						OK = pd.removeOutliers( // getPlaneFromMeas should already have run
								disp_strength,
								targetV,      // double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
								max_outliers, // int        maxRemoved,  // maximal number of tiles to remove (not a constant)
								debugLevel); // int        debugLevel)
						if (!OK) {
							continue;
						}
						if (debugLevel > 0) {
							if (pd.getWeight() > 0.0) { // 1.0) {
								System.out.println("Removed outliers "+ suffix +
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
					if (debugLevel > 0) {
						System.out.println("World normal " + suffix + " = {"+
								norm_xyz[0]+", "+norm_xyz[1]+", "+norm_xyz[2]+"}");

					}
					// calculate the world planes too
					pd.getWorldPlaneFromMeas(
							// use current selection, possibly reduced after removeOutliers()
							null, // plane_selections[ps], // tile_sel,       // boolean [][] tile_sel, // null - do not use, {} use all (will be modified)
							disp_strength,
							Double.NaN,    // double       disp_far, // minimal disparity to select (or NaN)
							Double.NaN,    // double       disp_near, // maximal disparity to select (or NaN)
							dispNorm,      // 0.0,            // plDispNorm,  // double       dispNorm,   //  Normalize disparities to the average if above
							0.0,            // double       min_weight,
							min_tiles,      // int          min_tiles,
//							strength_floor, //
//							strength_pow,   // double       strength_pow,
							// update !
							smplMode,
							mlfp,

//							smplSide,
//							smplNum,
//							smplRms,

//							smplWnd,  // use window functions for the samples

//							max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//							max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//							damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//							min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//							transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//							far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//							far_power,     //    1.0; // Raise disparity to this power before averaging for far objects


							debugLevel);
				}
			}
			if (st_planes.size() > 0){
				// sort planes by increasing disparity (tile center or plane center ? ) Using plane center
				Collections.sort(st_planes, new Comparator<PlaneData>() {
					@Override
					public int compare(PlaneData lhs, PlaneData rhs) {
						// -1 - less than, 1 - greater than, 0 - equal
						return (rhs.getZxy()[0] > lhs.getZxy()[0]) ? -1 : (rhs.getZxy()[0] < lhs.getZxy()[0] ) ? 1 : 0;
					}
				});
				return st_planes;
			}
			return null;
		}

		public boolean [][][]  reDiscriminateTiles_0(
				String           prefix,
				final PlaneData  [] planes,
				final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

				final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				final MeasuredLayersFilterParameters mlfp,

//				final int        smplSide, //        = 2;      // Sample size (side of a square)
//				final int        smplNum, //         = 3;      // Number after removing worst
//				final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				final boolean    smplWnd,  // use window functions for the samples

//				final double     max_abs_tilt,  //  2.0;   // pix per tile
//				final double     max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//				final double     damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				final double     min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				final double     transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				final int        far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				final double     far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

				final double     max_disp_diff,    // maximal disparity difference from the plane to consider tile
				final double     disp_range,       // parallel move known planes around original know value for the best overall fit
				final int        amplitude_steps,  // number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
				final int        num_variants,     // total number of variants to try (protect from too many planes)
				final int        mode, // 0 - weighted, 1 - equalized, 2 - best, 3 - combined
				final int        debugLevel)
		{
			final int size2 = 4 * superTileSize*superTileSize;
			final double max_disp_diff2 = max_disp_diff*max_disp_diff;
			TileNeibs tileNeibs = new TileNeibs(2 * stSize, 2 * stSize);
			if (planes == null) return null;
			// create a list of usable planes according to the mode
			ArrayList<PlaneData> tilePlanes = new ArrayList<PlaneData>();
			for (int np = 0; np < planes.length; np++) if (planes[np] != null){
				PlaneData weighted_pd = planes[np].getNonexclusiveStar();
				PlaneData equal_pd =    planes[np].getNonexclusiveStarEq();
				switch (mode){
				case 0:
					if (weighted_pd == null){
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, trying getNonexclusiveStarEq");
						if (equal_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(equal_pd);
						}
//						throw new IllegalArgumentException ("refineDiscriminateTiles(): getNonexclusiveStar() returned null");
					} else  {
						tilePlanes.add(weighted_pd);
					}
					break;
				case 1:
					if (equal_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles(): getNonexclusiveStarEq() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, trying getNonexclusiveStar");
						if (weighted_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(weighted_pd);
						}

					} else {
						tilePlanes.add(equal_pd);
					}
					break;
				case 2:
					if (weighted_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, trying getNonexclusiveStarEq");
						if (equal_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(equal_pd);
						}
						break;
					}
					if (equal_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, trying getNonexclusiveStar");
						tilePlanes.add(weighted_pd);
						break;
					}
					if (weighted_pd.getWeight() > equal_pd.getWeight()) tilePlanes.add(weighted_pd);
					else                                                tilePlanes.add(equal_pd);
					break;
				case 3:
					if (weighted_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, trying getNonexclusiveStarEq");
						if (equal_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(equal_pd);
						}
						break;
					}
					if (equal_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, trying getNonexclusiveStar");
						tilePlanes.add(weighted_pd);
						break;
					}
					PlaneData combo_pd =  weighted_pd.mergePlaneToThis(
							equal_pd, // PlaneData otherPd,
							1.0,      // double    scale_other,
							1.0,      // double    starWeightPwr,    // Use this power of tile weight when calculating connection cost
							false,    // boolean   ignore_weights,
							true,     // boolean   sum_weights,
							preferDisparity, // boolean   preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
							debugLevel - 3); // int       debugLevel)
					 tilePlanes.add(combo_pd);
					break;
				default:
		    		throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": invalid mode="+mode);
				}
			}
			if (tilePlanes.isEmpty()){
				return null;
			}
			if (debugLevel > 1)	{
				System.out.println ("refineDiscriminateTiles() "+prefix+" - step 1");
			}

           // get measured disparity/strength data, filtered, not tilted
			double [][][] disp_strength = new double[measuredLayers.getNumLayers()][][];
			for (int ml = 0; ml < disp_strength.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
				if (smplMode) {
					disp_strength[ml] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
							ml, // int num_layer,
							getSTileXY()[0],        // int stX,
							getSTileXY()[1],        // int stY,
							null,                   // boolean [] sel_in, - use all
							mlfp,
//							strength_floor,
//							measured_strength_pow,  //
//							smplSide,               // = 2;      // Sample size (side of a square)
//							smplNum,                // = 3;      // Number after removing worst
//							smplRms,                // = 0.1;    // Maximal RMS of the remaining tiles in a sample
//							smplWnd,                // use window functions for the samples

//							max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//							max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//							damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//							min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//							transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//							far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//							far_power,     //    1.0; // Raise disparity to this power before averaging for far objects

							true,                   // boolean null_if_none)
							debugLevel);
				} else {
					disp_strength[ml] = measuredLayers.getDisparityStrengthML(
							ml,                     // int num_layer,
							getSTileXY()[0],        // int stX,
							getSTileXY()[1],        // int stY,
							null,                   // boolean [] sel_in, - use all
							mlfp,
//							strength_floor,         //  double strength_floor,
//							measured_strength_pow,  // double strength_pow,
							true);                  // boolean null_if_none);
				}
			}
//			double [] window =	getWindow (2 * superTileSize);
			if (debugLevel > 1)	{
				System.out.println ("refineDiscriminateTiles() "+prefix+" - step 2");
			}
			// get all original planes
			int num_planes = tilePlanes.size();
			double [][][] pds = new double [num_planes][][];
			for (int np = 0; np < pds.length; np++){
				PlaneData pd = tilePlanes.get(np);
				pds[np] = pd.getDoublePlaneDisparityStrength(
						false, // boolean   useWorld,
						null, // double [] window,
						-1, // int       dir,
						false, // boolean   use_sel,
						false, // boolean   divide_by_area,
						1.0, // double    scale_projection,
						0.0, // double    fraction_uni,
						0); // int       debugLevel)
			}
			double [][][] flatness = new double [num_planes][disp_strength.length][];
			int [][][] num_cells = new int [num_planes][disp_strength.length][];
			for (int np = 0; np < pds.length; np++){
				for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
					flatness[np][ml] = new double[size2];
					num_cells[np][ml] = new int[size2];
					for (int indx = 0; indx < size2; indx++) if (disp_strength[ml][1][indx] > 0.0){
						double d0 = pds[np][0][indx]; // disp_strength[ml][0][indx];
						double sw = disp_strength[ml][1][indx], sd = 0.0, sd2 = 0.0;
						num_cells[np][ml][indx]++;
						for (int dir = 0; dir<8; dir++){
							int indx1 = tileNeibs.getNeibIndex(indx, dir);
							if (indx1 >= 0){
								double w =  disp_strength[ml][1][indx1];
								if (w > 0.0){
									double d = disp_strength[ml][0][indx1] - d0;
									sw += w;
									sd += w * d;
									sd2 += w * d * d;
									num_cells[np][ml][indx]++;
								}
							}
						}
						if (sw > 0.0) {
							sd /= sw;
							sd2 /= sw;
							sd2 -= sd * sd;
						}
						flatness[np][ml][indx] = sd2;
					}
				}
			}
			if (debugLevel > 2){
				int dbg_ml = 0; // to protect from different layer configuration
				for (; dbg_ml < disp_strength.length;  dbg_ml++){
					if (disp_strength[dbg_ml] != null) break;
				}
				//disp_strength[dbg_ml][0], disp_strength[dbg_ml][1] * 10, disp_strength[dbg_ml][0] - pds[np][0], pds[np][0]
				String [] dbg_titles = new String[2 + 5 * num_planes];
				double [][] dbg_img = new double [dbg_titles.length][];
				dbg_titles[0] = "disp_"+dbg_ml;
				dbg_titles[1] = "str_"+dbg_ml;
				dbg_img[0] = disp_strength[dbg_ml][0];
				dbg_img[1] = disp_strength[dbg_ml][1];
				for (int np = 0; np < num_planes; np++){
					dbg_titles[2 + 0 * num_planes + np] = "pln_"+np;
					dbg_titles[2 + 1 * num_planes + np] = "diff_"+np;
					dbg_titles[2 + 2 * num_planes + np] = "flat_"+np;
					dbg_titles[2 + 3 * num_planes + np] = "rflat_"+np;
					dbg_titles[2 + 3 * num_planes + np] = "nrflat_"+np;
					dbg_img[2 + 0 * num_planes + np] = pds[np][0];
					dbg_img[2 + 1 * num_planes + np] = disp_strength[dbg_ml][0].clone();
					dbg_img[2 + 2 * num_planes + np] = flatness[np][dbg_ml];
					dbg_img[2 + 3 * num_planes + np] = new double[size2];
					dbg_img[2 + 4 * num_planes + np] = new double[size2];
					for (int i = 0; i < size2; i++) {
						dbg_img[2 + 1 * num_planes + np][i] -= pds[np][0][i];
						dbg_img[2 + 3 * num_planes + np][i] = 1.0 / (flatness[np][dbg_ml][i] + 0.001);
						dbg_img[2 + 4 * num_planes + np][i] = num_cells[np][dbg_ml][i] * 1.0 / (flatness[np][dbg_ml][i] + 0.001);
						if (disp_strength[dbg_ml][1][i] == 0.0){
							dbg_img[2 + 1 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 2 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 3 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 4 * num_planes + np][i] = Double.NaN;
						}
					}
				}
				showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
				sdfa_instance.showArrays(dbg_img,     2 * superTileSize, 2 * superTileSize, true, "refine-"+prefix,dbg_titles);
			}




			// calculate number of the variants for each plane
			int  extra_vars = (((int) Math.pow(num_variants, 1.0/num_planes)) - 1) / 2; // 0 - single, 1 - 3 (1 each direction), ...
			if (extra_vars > amplitude_steps) extra_vars = amplitude_steps; //
			int steps = 2 * extra_vars + 1;
			double [] disps = new double [steps];
			disps[0] = 0.0; // center
			for (int i = 0; i < extra_vars; i++){
				disps[2 * i + 1] =  (disp_range * (i + 1)) / extra_vars;
				disps[2 * i + 2] = -disps[2 * i + 1];
			}
			int [] state = new int [num_planes]; // variants counter
			boolean [][][] best_selections = null;
//			int best_num_tiles = 0;
			double best_cost = Double.NaN;
			if (debugLevel > 1)	{
				System.out.println ("refineDiscriminateTiles() "+prefix+" - step 3");
			}
			while (true){
				if (debugLevel > 1)	{
					System.out.print ("refineDiscriminateTiles() "+prefix+" - state:");
					for (int np = 0; np < num_planes; np++){
						System.out.print(" "+state[np]);
					}
					System.out.println();
				}
				// evaluate current variant
				double [] weights = new double [num_planes];
				double [] diffs0  = new double [num_planes];
				boolean [][][] var_sels = new boolean [num_planes][disp_strength.length][];
				for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
					for (int nTile = 0; nTile < size2; nTile++) if (disp_strength[ml][1][nTile] != 0){
						// find best fit
						int best_np = -1;
						double best_err2 = max_disp_diff2;
						for (int np = 0; np < num_planes; np++){
							double e2 = disp_strength[ml][0][nTile] - (pds[np][0][nTile] + disps[state[np]]);
							e2 *= e2;
							if (e2 < best_err2){
								best_np = np;
								best_err2 = e2;
							}
						}
						if (best_np >= 0) { // found acceptable assignment
							double w = disp_strength[ml][1][nTile]; //  * window[nTile];
							double d = disp_strength[ml][0][nTile] - (pds[best_np][0][nTile] + disps[state[best_np]]);
							weights[best_np] +=  w;
							diffs0[best_np] +=  w*d;
						}
					}
				}
				if (debugLevel > 1)	{
					System.out.println ("refineDiscriminateTiles() "+prefix+" - pass1 DONE");
				}
					// second pass - use new averages disparity offsets, that may change assignments
				for (int np = 0; np < num_planes; np++){
					if (weights[np] > 0.0) diffs0[np] /= weights[np];
				}
				weights =           new double [num_planes];
				double [] diffs  =  new double [num_planes];
				double [] diffs2  = new double [num_planes];
				int    [] ntiles =  new int [num_planes];
				int num_tiles = 0;
				for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
					for (int np = 0; np < num_planes; np++){
						var_sels[np][ml] = new boolean [size2];
					}
					for (int nTile = 0; nTile < size2; nTile++) if (disp_strength[ml][1][nTile] != 0){
						// find best fit
						int best_np = -1;
						double best_err2 = max_disp_diff2;
						for (int np = 0; np < num_planes; np++){
							double e2 = disp_strength[ml][0][nTile] - (pds[np][0][nTile] + disps[state[np]]+ diffs0[np]);
							e2 *= e2;
							if (e2 < best_err2){
								best_np = np;
								best_err2 = e2;
							}
						}
						if (best_np >= 0) { // found acceptable assignment
							double w = disp_strength[ml][1][nTile]; //  * window[nTile];
							double d = disp_strength[ml][0][nTile] - (pds[best_np][0][nTile] + disps[state[best_np]]);
							weights[best_np] +=  w;
							diffs[best_np] +=    w*d;
							diffs2[best_np] +=   w*d*d;
							ntiles[best_np]++;
							var_sels[best_np][ml][nTile] = true;
							num_tiles ++;
						}
					}
				}
				for (int np = 0; np < num_planes; np++) if (weights[np] > 0.0){
					diffs[np] /= weights[np];
					diffs2[np] /= weights[np];
					diffs2[np] -= diffs[np]*diffs[np];
				}
				// How to decide if this variant is better? Have number of tiles that fit, per-plane weights and per-plane rms=sqrt(diffs2)
				// start with just number of tiles fit
				double cost = 0;
				double weight = 0.0;
				for (int np = 0; np < num_planes; np++) {
					weight += weights[np];
					cost += weights[np]*diffs2[np];
				}

				if ((weight > 0.0) && (num_tiles > 0)){
					double corr_cost = cost/num_tiles;
					if (debugLevel > 1)	{
						System.out.print ("refineDiscriminateTiles() "+prefix+": [");
						for (int np = 0; np < num_planes; np++){
							System.out.print (" "+state[np]);
						}
						System.out.println ("] num_tiles = "+num_tiles+" weight="+weight+" corr_cost="+corr_cost+" cost="+cost);
					}

					if (Double.isNaN(best_cost) || (best_cost > corr_cost)) {
						best_cost = corr_cost;
						best_selections = var_sels;

					}
				}

//				if (num_tiles > best_num_tiles) {
//					best_num_tiles = num_tiles;
//					best_selections = var_sels;
//				}
				if (debugLevel > 1)	{
					System.out.println ("refineDiscriminateTiles() "+prefix+" - pass2 DONE");
				}
				//  calculate next variant
				boolean all_done = true;
				for (int i = 0; i < num_planes; i ++){
					if (state[i] < (steps - 1)){
						state[i] ++;
						for (int j = 0; j < i; j++){
							state[j] = 0;
						}
						all_done = false;
						break;
					}
				}
				if (all_done){
					break;
				}

			}
			return best_selections;
		}

		/**
		 * re-discriminate tiles between supertiles in the list using hints from the already calculated planes and their neighbors
		 * @param prefix text to be used in debug output (such as supertile number)
		 * @param planes array of PlaneData instances of the current supertile to be used as hints
		 * @param merge_planes should be initilaized to int[2], will return indices of planes suggested to be merged or [-1, -1]
		 * @param stMeasSel select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
		 * @param dispNorm normalize disparities to the average if above
		 * @param smplMode use sample mode (false - each tile independent)
		 * @param smplSide sample size (side of a square)
		 * @param smplNum number after removing worst
		 * @param smplRms maximal RMS of the remaining tiles in a sample
		 * @param smplWnd use window functions for the samples
		 *
		 * @param max_abs_tilt pix per tile
		 * @param max_rel_tilt pix / disparity) per tile
		 * @param damp_tilt Damp tilt to handle insufficient  (co-linear)data
		 * @param min_tilt_disp Disparity switch between filtering modes - near objects use tilts, far - use max disparity
		 * @param transition Mode transition range (between tilted and maximal disparity)
		 * @param far_mode Far objects filtering mode (0 - off, 1 - power of disparity)
		 * @param far_power Raise disparity to this power before averaging for far objects
		 *
		 * @param disp_tolerance maximal disparity difference from the plane to consider tile
		 * @param disp_var_floor squared add to variance to calculate reverse flatness (used mostly for single-cell clusters) (reuse disp_sigma)?
		 * @param disp_sigma Gaussian sigma to compare how measured data is attracted to planes
		 * @param disp_range full range of the parallel plane shifts to evaluate (histograms full range)
		 * @param amplitude_steps numer of histogram steps in each direction from the center
		 * @param hist_blur histogram LPF sigma
		 * @param exclusivity 1.0 - tile belongs to one plane only, 0.0 - regardless of others
		 * @param exclusivity2 when exclusivity > 1.0, add tiles to the existing clusters if attraction > this relaxed level
		 * @param exclusivity_strict do not add to the clusters if there is a single offending neighbor (false just more of the these neighbors than offenders)
		 * @param attractionCorrMax do not discriminate if at least one pair has attraction correlation above this level
         * @param attractionCorrMerge attraction to different planes correlation that is high enough to merge planes
         * @param plDiscrSteal if offender has this number of tiles (including center) the cell can not be used
         * @param plDiscrGrown only use tiles within this range from original selection, < 0 - disable this filter
         * @param outliersXMedian remove outliers from the final selection that have distance more than scaled median
		 * @param mode what neighbor-dependent pre-calculated plane to use as hints: 0 - weighted, 1 - equalized, 2 - best, 3 - combined
		 * @param debugLevel debug level
		 * @return per-plane, per-measurement layer, per tile index - use this tile. Return null if could not discriminate,
		 *         in that case merge_planes may contain non-negative indices to be merged and the whole method re-ran
		 */
		public boolean [][][]  reDiscriminateTiles(
				String           prefix,
				final PlaneData  [] planes,
				final int []     merge_planes, // indices of planes suggested to be merged
				final int        stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				final double     dispNorm,   //  Normalize disparities to the average if above

				final boolean    smplMode, //        = true;   // Use sample mode (false - regular tile mode)
				final MeasuredLayersFilterParameters mlfp,
//				final int        smplSide, //        = 2;      // Sample size (side of a square)
//				final int        smplNum, //         = 3;      // Number after removing worst
//				final double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				final boolean    smplWnd,  // use window functions for the samples

//				final double     max_abs_tilt,  //  2.0;   // pix per tile
//				final double     max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//				final double     damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				final double     min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				final double     transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				final int        far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				final double     far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

				final double     disp_tolerance,   // maximal disparity difference from the plane to consider tile
				final double     disp_var_floor,   // squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
				final double     disp_sigma,       // G.sigma to compare how measured data is attracted to planes
				final double     disp_range,       // parallel move known planes around original know value for the best overall fit
				final int        amplitude_steps,  // number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
				final double     hist_blur,        // Sigma to blur histogram
				final double     exclusivity,      // 1.0 - tile belongs to one plane only, 0.0 - regardless of others
			    final double     exclusivity2,        //        =   0.8;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
			    final boolean    exclusivity_strict,  //         = true;   // When growing selection do not allow any offenders around (false - more these than others)
				final double     attractionCorrMax,   //         = 0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
				final double     attractionCorrMerge, //         = 0.85;  // Attraction to different planes correlation that is high enough to merge planes
				final int        plDiscrSteal,        //         =   4;     // If offender has this number of tiles (including center) the cell can not be used
				final int        plDiscrGrown,        //         =   0;     // Only use tiles within this range from original selection
				final double     outliersXMedian,     //         = 1.5;   // Remove outliers from the final selection that have distance more than scaled median

				final int        mode, // 0 - weighted, 1 - equalized, 2 - best, 3 - combined
				final int        debugLevel)
		{
			final int size2 = 4 * superTileSize*superTileSize;
			final double disp_tolerance2 = disp_tolerance * disp_tolerance;
			final double floor2 = disp_var_floor * disp_var_floor;
			final double ksigma = 0.5/(disp_sigma * disp_sigma);

			TileNeibs tileNeibs = new TileNeibs(2 * stSize, 2 * stSize);
			if (planes == null) return null;
			// create a list of usable planes according to the mode
			ArrayList<PlaneData> tilePlanes = new ArrayList<PlaneData>();
			for (int np = 0; np < planes.length; np++) if (planes[np] != null){
				PlaneData weighted_pd = planes[np].getNonexclusiveStar();
				PlaneData equal_pd =    planes[np].getNonexclusiveStarEq();
				switch (mode){
				case 0:
					if (weighted_pd == null){
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, trying getNonexclusiveStarEq");
						if (equal_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(equal_pd);
						}
//						throw new IllegalArgumentException ("refineDiscriminateTiles(): getNonexclusiveStar() returned null");
					} else  {
						tilePlanes.add(weighted_pd);
					}
					break;
				case 1:
					if (equal_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles(): getNonexclusiveStarEq() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, trying getNonexclusiveStar");
						if (weighted_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(weighted_pd);
						}

					} else {
						tilePlanes.add(equal_pd);
					}
					break;
				case 2:
					if (weighted_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, trying getNonexclusiveStarEq");
						if (equal_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(equal_pd);
						}
						break;
					}
					if (equal_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, trying getNonexclusiveStar");
						tilePlanes.add(weighted_pd);
						break;
					}
					if (weighted_pd.getWeight() > equal_pd.getWeight()) tilePlanes.add(weighted_pd);
					else                                                tilePlanes.add(equal_pd);
					break;
				case 3:
					if (weighted_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStar() returned null, trying getNonexclusiveStarEq");
						if (equal_pd == null) {
							if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, using plane itself");
							tilePlanes.add(planes[np]);
						} else {
							tilePlanes.add(equal_pd);
						}
						break;
					}
					if (equal_pd == null) { //throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null");
						if (debugLevel > 0)	System.out.println ("refineDiscriminateTiles() "+prefix+":"+np+": getNonexclusiveStarEq() returned null, trying getNonexclusiveStar");
						tilePlanes.add(weighted_pd);
						break;
					}
					PlaneData combo_pd =  weighted_pd.mergePlaneToThis(
							equal_pd, // PlaneData otherPd,
							1.0,      // double    scale_other,
							1.0,      // double    starWeightPwr,    // Use this power of tile weight when calculating connection cost
							false,    // boolean   ignore_weights,
							true,     // boolean   sum_weights,
							preferDisparity, // boolean   preferDisparity, // Always start with disparity-most axis (false - lowest eigenvalue)
							debugLevel - 3); // int       debugLevel)
					 tilePlanes.add(combo_pd);
					break;
				default:
		    		throw new IllegalArgumentException ("refineDiscriminateTiles() "+prefix+":"+np+": invalid mode="+mode);
				}
				tilePlanes.get(tilePlanes.size()-1).setMark0(np);
			}
			if (tilePlanes.isEmpty()){
				return null;
			}
			if (debugLevel > 1)	{
				System.out.println ("refineDiscriminateTiles() "+prefix+" - step 1");
			}

           // get measured disparity/strength data, filtered, not tilted
			double [][][] disp_strength = new double[measuredLayers.getNumLayers()][][];
			for (int ml = 0; ml < disp_strength.length; ml++) if ((stMeasSel & ( 1 << ml)) != 0){
				if (smplMode) {
					disp_strength[ml] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
							ml, // int num_layer,
							getSTileXY()[0],        // int stX,
							getSTileXY()[1],        // int stY,
							null,                   // boolean [] sel_in, - use all
							mlfp,
//							strength_floor,
//							measured_strength_pow,  //
//							smplSide,               // = 2;      // Sample size (side of a square)
//							smplNum,                // = 3;      // Number after removing worst
//							smplRms,                // = 0.1;    // Maximal RMS of the remaining tiles in a sample
//							smplWnd,                // use window functions for the samples
//							max_abs_tilt,  // 2.0; // Maximal absolute tilt in pixels/tile
//							max_rel_tilt,  // 0.2; // Maximal relative tilt in pixels/tile/disparity
//							damp_tilt,     //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//							min_tilt_disp, // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//							transition,    // 1.0; // Mode transition range (between tilted and maximal disparity)
//							far_mode,      //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//							far_power,     //    1.0; // Raise disparity to this power before averaging for far objects
							true,                  // boolean null_if_none)
							debugLevel);
				} else {
					disp_strength[ml] = measuredLayers.getDisparityStrengthML(
							ml,                     // int num_layer,
							getSTileXY()[0],        // int stX,
							getSTileXY()[1],        // int stY,
							null,                   // boolean [] sel_in, - use all
							mlfp,
//							strength_floor,         //  double strength_floor,
//							measured_strength_pow,  // double strength_pow,
							true);                  // boolean null_if_none);
				}
			}
			double [] window =	getWindow (2 * superTileSize);
			if (debugLevel > 1)	{
				System.out.println ("refineDiscriminateTiles() "+prefix+" - step 2");
			}
			// get all original planes
			int num_planes = tilePlanes.size();
			double [][][] pds = new double [num_planes][][];
			for (int np = 0; np < pds.length; np++){
				PlaneData pd = tilePlanes.get(np);
				pds[np] = pd.getDoublePlaneDisparityStrength(
						false, // boolean   useWorld,
						null, // double [] window,
						-1, // int       dir,
						false, // boolean   use_sel,
						false, // boolean   divide_by_area,
						1.0, // double    scale_projection,
						0.0, // double    fraction_uni,
						0); // int       debugLevel)
			}
			// Create enable and disable masks

			boolean [][] prev_used = new boolean [num_planes][size2];
			boolean [][] mask =      new boolean [num_planes][];
			boolean [][] used_strong = new boolean [num_planes][size2];
			for (int np = 0; np < num_planes; np++){
				boolean [][] meas_sel = tilePlanes.get(np).getMeasSelection();
				for (int ml = 0; ml < meas_sel.length; ml++) if (meas_sel[ml] != null){
					for (int indx = 0; indx < size2; indx++){
						prev_used[np][indx] |= meas_sel[ml][indx];
					}
				}
				mask[np] = prev_used[np].clone();
				if (plDiscrGrown > 100) { // >0
					tileNeibs.growSelection(
							plDiscrGrown, // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
							mask[np], // boolean [] tiles,
							null); // boolean [] prohibit)
				} else if (plDiscrGrown < 100){ //< 0
					for (int indx = 0; indx < size2; indx++) mask[np][indx] = true;

				}
				if (plDiscrSteal > 0){
					for (int indx = 0; indx < size2; indx++) if (prev_used[np][indx]){
						int num_neibs = 1;
						for (int dir = 0; dir < 8; dir++){
							int indx1 = tileNeibs.getNeibIndex(indx, dir);
							if ((indx1 >=0) && prev_used[np][indx1]){
								num_neibs++;
							}
						}
						used_strong[np][indx] = (num_neibs >= plDiscrSteal);
					}
				}
			}
			if (0 * plDiscrSteal > 0){ // 0*
				for (int np = 0; np < num_planes; np++){
					for (int np1 = 0; np1 < num_planes; np1++) if (np1 != np){
						for (int indx = 0; indx < size2; indx++){
							mask[np][indx] &= !used_strong[np1][indx];
						}
					}
				}
			}

			double [][][] flatness = new double [num_planes][disp_strength.length][];
			double [][][] norm_flatness = new double [num_planes][disp_strength.length][];
			int [][][] num_cells = new int [num_planes][disp_strength.length][];
			for (int np = 0; np < pds.length; np++){
				for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
					flatness[np][ml] = new double[size2];
					norm_flatness[np][ml] = new double[size2];
					num_cells[np][ml] = new int[size2];
					for (int indx = 0; indx < size2; indx++) if (disp_strength[ml][1][indx] > 0.0){
						double d0 = pds[np][0][indx]; // disp_strength[ml][0][indx];
//						double sw = disp_strength[ml][1][indx], sd = 0.0, sd2 = 0.0;
						double sw = disp_strength[ml][1][indx]/window[indx], sd = 0.0, sd2 = 0.0;
						num_cells[np][ml][indx]++;
						for (int dir = 0; dir<8; dir++){
							int indx1 = tileNeibs.getNeibIndex(indx, dir);
							if (indx1 >= 0){
//								double w =  disp_strength[ml][1][indx1];
								double w =  disp_strength[ml][1][indx1]/window[indx1];
								if (w > 0.0){
									double d = disp_strength[ml][0][indx1] - d0;
									sw += w;
									sd += w * d;
									sd2 += w * d * d;
									num_cells[np][ml][indx]++;
								}
							}
						}
						if (sw > 0.0) {
							sd /= sw;
							sd2 /= sw;
							sd2 -= sd * sd;
						}
						double rms = Math.sqrt(sd2);
						flatness[np][ml][indx] = rms; // sd2; // will not be used
//						norm_flatness[np][ml][indx] = num_cells[np][ml][indx]/(sd2 + floor2);
//						norm_flatness[np][ml][indx] = (num_cells[np][ml][indx] > 1) ? (1.0/sd2): (1.0/floor2);
						if (rms < disp_var_floor){
							rms = disp_var_floor;
						}
						norm_flatness[np][ml][indx] = (num_cells[np][ml][indx] > 1) ? (1.0/rms): (1.0/disp_var_floor);
					}
				}
			}
			// calculate histograms for each plane, combining weight with flatness
			final int steps = amplitude_steps * 2 + 1;
			double bin_size = disp_range/(2 * amplitude_steps);
			double k_bin = 2 * amplitude_steps /  disp_range; // 1/bin_size
			double [][] histograms = new double [num_planes][steps];

			for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
				for (int indx = 0; indx < size2; indx++){
					double weight = disp_strength[ml][1][indx];
					if (weight > 0.0) {
						double disp =   disp_strength[ml][0][indx];
						for (int np = 0; np < num_planes; np++) if (mask[np][indx]){
							double d = disp - pds[np][0][indx];
							double db = d * k_bin + amplitude_steps;
							int bin = (int) Math.round (db);
							if ((bin >= 0) && (bin < steps)) {
								// combine correlation weight, variance between valid neighbors and number of valid neighbors
//								double w = weight /(flatness[np][ml][indx] + floor2) * num_cells[np][ml][indx]; // TODO: tweak?
								double w = weight * norm_flatness[np][ml][indx]; // TODO: tweak?
								histograms[np][bin] += w;
							}
						}
					}
				}
			}
			// Raise weight of the center part of the histogram - there can be originally distinct parallel planes
			// alternatively - use planes_mod or star* and no histograms at all?
			final double initial_trust_sigma = 0.3 * disp_tolerance; // TODO: use dispNorm !
			for (int np = 0; np < num_planes; np++) {
				double k_sigma = 0.5/(initial_trust_sigma*initial_trust_sigma);
				for (int nb = 0; nb < steps; nb++) {
					double offs = bin_size * (nb - amplitude_steps);
					histograms[np][nb] *= Math.exp(-k_sigma*offs*offs);
				}
			}

			// normalize histograms
			for (int np = 0; np < num_planes; np++) {
				double s = 0.0;
				for (int nb = 0; nb < steps; nb++)
					s += histograms[np][nb];
				if (s > 0.0) {
					for (int nb = 0; nb < steps; nb++)
						histograms[np][nb] /= s;
				}
			}
			if (debugLevel > 2) {
				System.out.println("refineDiscriminateTiles() histograms raw:");
				System.out.print("disparity: ");
				for (int nb = 0; nb < steps; nb++) {
					System.out.print(bin_size * (nb - amplitude_steps));
					if (nb < (steps - 1))
						System.out.print(", ");
				}
				System.out.println();
				for (int np = 0; np < num_planes; np++) {
					System.out.print(np + ": ");
					for (int nb = 0; nb < steps; nb++) {
						System.out.print(histograms[np][nb]);
						if (nb < (steps - 1))
							System.out.print(", ");
					}
					System.out.println();
				}
				System.out.println();
				System.out.println();
			}

			// Add LPF here?
			if (hist_blur > 0.0) {
				DoubleGaussianBlur gb=new DoubleGaussianBlur();
				for (int np = 0; np < num_planes; np++) {
					gb.blur1Direction(histograms[np], steps, 1, hist_blur/bin_size, 0.01,true);
				}
				if (debugLevel > 2) {
					System.out.println("refineDiscriminateTiles() histograms blured:");
					System.out.print("disparity: ");
					for (int nb = 0; nb < steps; nb++) {
						System.out.print(bin_size * (nb - amplitude_steps));
						if (nb < (steps - 1))
							System.out.print(", ");
					}
					System.out.println();
					for (int np = 0; np < num_planes; np++) {
						System.out.print(np + ": ");
						for (int nb = 0; nb < steps; nb++) {
							System.out.print(histograms[np][nb]);
							if (nb < (steps - 1))
								System.out.print(", ");
						}
						System.out.println();
					}
					System.out.println();
				}

			}


			// For each plane, regardless of others
			double [] offsets = new double [num_planes];
			for (int np = 0; np < num_planes; np++) {
				int max_bin = 0;
				for (int nb = 1; nb < steps; nb++){
					if (histograms[np][nb] > histograms[np][max_bin]){
						max_bin = nb;
					}
				}
				// improve max by second degree polynomial
				 offsets[np] = bin_size * (max_bin - amplitude_steps);
				if ((max_bin > 0) && (max_bin < (steps-1))){
					offsets[np] += bin_size * 0.5 *(histograms[np][max_bin + 1] - histograms[np][max_bin - 1]) /
							(histograms[np][max_bin + 1] + histograms[np][max_bin - 1] - 2 * histograms[np][max_bin]);
				}
			}

			// for each tile, each plane calculate "attraction" of the tile to each plane, where attraction is defined by Gaussian
			// of normalized disparity difference and "flatness" for each plane. Then, if the attraction exceeds "exclusivity" times
			// best competitor - include this tile for that plane
			double  [][][] attractions = new double[num_planes][disp_strength.length][];
			for (int np = 0; np < num_planes; np++) {
				for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
					attractions[np][ml] = new double[size2];
				}
			}
			for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
				for (int indx = 0; indx < size2; indx++) if (disp_strength[ml][1][indx] > 0){
					for (int np = 0; np < num_planes; np++) if (mask[np][indx]){
						double d1 = pds[np][0][indx] + offsets[np]; // shifted plane
						double d2 = disp_strength[ml][0][indx];
						double dav = 0.5 * (d1 + d2);
						double diff = d2 - d1;
						if (dav > dispNorm) diff *= dispNorm/dav;
						double diff2 = diff * diff;
						if (diff2 <= disp_tolerance2) {
							attractions[np][ml][indx] = Math.exp(-ksigma*diff2) * norm_flatness[np][ml][indx];

						}
					}
				}
			}
			double [][] attr_corr = new double [num_planes][num_planes];
			for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
				for (int indx = 0; indx < size2; indx++) if (disp_strength[ml][1][indx] > 0){
					for (int np = 0; np < num_planes; np++) {
						for (int np1 = np; np1 < num_planes; np1++) {
							attr_corr[np][np1] += attractions[np][ml][indx] * attractions[np1][ml][indx];
							if (np1 > np){
								attr_corr[np1][np] = attr_corr[np][np1];
							}
						}
					}
				}
			}
			for (int np = 0; np < num_planes; np++) {
				for (int np1 = np+1; np1 < num_planes; np1++) {
					if ((attr_corr[np][np] > 0.0) && (attr_corr[np1][np1] > 0.0)){
						attr_corr[np][np1] /= Math.sqrt(attr_corr[np][np] * attr_corr[np1][np1]);
						attr_corr[np1][np] = attr_corr[np][np1];
					}
				}
				attr_corr[np][np] = 1.0;

			}
			double max_attr_corr = 0.0;
			int [] merge_pair = {-1,-1};
			for (int np = 0; np < num_planes; np++) {
				for (int np1 = np + 1; np1 < num_planes; np1++) {
					if (attr_corr[np][np1] > max_attr_corr) {
						max_attr_corr = attr_corr[np][np1];
						merge_pair[0] = np;
						merge_pair[1] = np1;
					}
				}
			}


			// discriminate
			boolean [][][] best_selections = new boolean[num_planes][disp_strength.length][];
			for (int np = 0; np < num_planes; np++) {
				for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
					best_selections[np][ml] = new boolean[size2];
				}
			}
			for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
				for (int indx = 0; indx < size2; indx++) if (disp_strength[ml][1][indx] > 0){
					if (num_planes == 1){
						best_selections[0][ml][indx] = attractions[0][ml][indx] > 0.0;
					} else { // find to best candidates
						int best_np = 0;
						for (int np = 1; np < num_planes; np++) if (attractions[np][ml][indx] > attractions[best_np][ml][indx]){
							best_np = np;
						}
						int best_np2 = (best_np ==0) ? 1:0;
						for (int np = 1; np < num_planes; np++) if ((np != best_np ) && (attractions[np][ml][indx] > attractions[best_np2][ml][indx])){
							best_np2 = np;
						}
						for (int np = 0; np < num_planes; np++){
							int np_other = (np == best_np) ? best_np2 : best_np;
							best_selections[np][ml][indx] = attractions[np][ml][indx] > exclusivity * attractions[np_other][ml][indx];
						}
					}
				}
			}
			if ((exclusivity > 1.0) && (exclusivity2 > 0.0) && (num_planes > 1)) { // second pass - lower threshold but depend on neighbors
				boolean changed = true;
				while (changed) {
					changed = false;
					boolean [][][] best_selections_prev = best_selections.clone();
					for (int np = 0; np < num_planes; np++) {
						best_selections_prev[np] = best_selections[np].clone();
						for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
							best_selections_prev[np][ml] = best_selections[np][ml].clone();
						}
					}

					for (int ml = 0; ml < disp_strength.length; ml++) if (disp_strength[ml] != null){
						for (int indx = 0; indx < size2; indx++) if (disp_strength[ml][1][indx] > 0){
							int best_np = 0;
							for (int np = 1; np < num_planes; np++) if (attractions[np][ml][indx] > attractions[best_np][ml][indx]){
								best_np = np;
							}
							int best_np2 = (best_np ==0) ? 1:0;
							for (int np = 1; np < num_planes; np++) if ((np != best_np ) && (attractions[np][ml][indx] > attractions[best_np2][ml][indx])){
								best_np2 = np;
							}
							for (int np = 0; np < num_planes; np++) if (!best_selections_prev[np][ml][indx]){ // only new
								int np_other = (np == best_np) ? best_np2 : best_np;
								if (attractions[np][ml][indx] > exclusivity2 * attractions[np_other][ml][indx]) {
									int num_this = 0;
									int num_other = 0;
									boolean used_by_other = false;
									for (int np1 = 0; np1 < num_planes; np1++) if (np1 != np){
										if (best_selections_prev[np1][ml][indx]) {
											used_by_other = true;
											break;
										}
									}
									if (! used_by_other) {
										for (int dir = 0; dir < 8; dir++){
											int indx1 = tileNeibs.getNeibIndex(indx, dir);
											if (indx1 >= 0){
												if (best_selections_prev[np][ml][indx1]) num_this ++;
												for (int np1 = 0; np1 < num_planes; np1++) if (np1 != np){
													if (best_selections_prev[np1][ml][indx1]) {
														num_other ++;
														break;
													}
												}
											}
										}
										if ((num_this > num_other) && (!exclusivity_strict || (num_other == 0))){ // make it stricter and require num_other==0 ?
											best_selections[np][ml][indx] = true;
											changed = true;
										}
									}
								}
							}
						}
					}
				}

			}
			boolean [][][] all_best_selections = null;
			if (outliersXMedian > 0.0) {
				all_best_selections = best_selections.clone();
				for (int np = 0; np < num_planes; np++) if (best_selections[np] != null){
					all_best_selections[np] = best_selections[np].clone();
					for (int ml = 0; ml < best_selections[np].length; ml++) if (best_selections[np][ml] != null){
						all_best_selections[np][ml] = best_selections[np][ml].clone();
					}
				}
				// from each plane remove outliers (without re-calculating mean value) that are farther from the mean than scaled median distance
				class ValIndex {
					int ml;
					int indx;
					double val;
					ValIndex(int ml, int indx, double val){
						this.ml = ml;
						this.indx = indx;
						this.val =  val;
					}
				}
				for (int np = 0; np < num_planes; np++) if (best_selections[np] != null){
					ArrayList<ValIndex> tile_list = new ArrayList<ValIndex>();
					for (int ml = 0; ml < best_selections[np].length; ml++) if (best_selections[np][ml] != null){
						for (int indx = 0; indx < best_selections[np][ml].length; indx++) if (best_selections[np][ml][indx]){
							double d1 = pds[np][0][indx] + offsets[np]; // shifted plane
							double d2 = disp_strength[ml][0][indx];
							double diff = d2 - d1;
//							double dav = 0.5 * (d1 + d2);
//							if (dav > dispNorm) diff *= dispNorm/dav;
							double diff2 = diff * diff;
							tile_list.add(new ValIndex(ml,indx,diff2));
						}
					}
					if (!tile_list.isEmpty()){
						Collections.sort(tile_list, new Comparator<ValIndex>() {
							@Override
							public int compare(ValIndex lhs, ValIndex rhs) {
								// -1 - less than, 1 - greater than, 0 - equal
								return (rhs.val > lhs.val) ? -1 : (rhs.val < lhs.val ) ? 1 : 0;
							}
						});
						int size = tile_list.size();
						double threshold =tile_list.get(size/2).val * outliersXMedian * outliersXMedian *1.0;
						for (int i = (outliersXMedian > 1.0) ? (size / 2 + 1) : 0; i < size; i++){
							if (tile_list.get(i).val > threshold) {
								best_selections[np][tile_list.get(i).ml][tile_list.get(i).indx] = false;
							}
						}
					}
					// now restore outliers if they had non-outliers neighbors
					ArrayList<Point> to_resore_list = new ArrayList<Point>();
					for (int ml = 0; ml < best_selections[np].length; ml++) if (best_selections[np][ml] != null){
						for (int indx = 0; indx < best_selections[np][ml].length; indx++) {
							if (all_best_selections[np][ml][indx] && !best_selections[np][ml][indx]){ // removed as outliers
								for (int dir = 0; dir<8; dir++){
									int indx1 = tileNeibs.getNeibIndex(indx, dir);
									if ((indx1 >= 0) && best_selections[np][ml][indx1]){
										to_resore_list.add(new Point(ml, indx));
										break;
									}
								}
							}
						}
					}
					for (Point p:to_resore_list) {
						best_selections[np][p.x][p.y] = true;
					}
				}
			}
			if (debugLevel > 2){
				int dbg_ml = 0; // to protect from different layer configuration
				for (; dbg_ml < disp_strength.length;  dbg_ml++){
					if (disp_strength[dbg_ml] != null) break;
				}
				//disp_strength[dbg_ml][0], disp_strength[dbg_ml][1] * 10, disp_strength[dbg_ml][0] - pds[np][0], pds[np][0]
				String [] dbg_titles = new String[2 + 10 * num_planes];
				double [][] dbg_img = new double [dbg_titles.length][];
				dbg_titles[0] = "disp_"+dbg_ml;
				dbg_titles[1] = "str_"+dbg_ml;
				dbg_img[0] = disp_strength[dbg_ml][0];
				dbg_img[1] = disp_strength[dbg_ml][1];
				for (int np = 0; np < num_planes; np++){
					dbg_titles[2 + 0 * num_planes + np] = "pln_"+np;
					dbg_titles[2 + 1 * num_planes + np] = "diff_"+np;
					dbg_titles[2 + 2 * num_planes + np] = "flat_"+np;
					dbg_titles[2 + 3 * num_planes + np] = "rflat_"+np;
					dbg_titles[2 + 4 * num_planes + np] = "nrflat_"+np;
					dbg_titles[2 + 5 * num_planes + np] = "mask_"+np;
					dbg_titles[2 + 6 * num_planes + np] = "attr_"+np;
					dbg_titles[2 + 7 * num_planes + np] = "pre_sel_"+np; // before outliers
					dbg_titles[2 + 8 * num_planes + np] = "sel_"+np; // add also old selection?
					dbg_titles[2 + 9 * num_planes + np] = "oldsel_"+np; // add also old selection?
					dbg_img[2 + 0 * num_planes + np] = pds[np][0];
					dbg_img[2 + 1 * num_planes + np] = disp_strength[dbg_ml][0].clone();
					dbg_img[2 + 2 * num_planes + np] = flatness[np][dbg_ml];
					dbg_img[2 + 3 * num_planes + np] = new double[size2];
					dbg_img[2 + 4 * num_planes + np] = norm_flatness[np][dbg_ml];
					dbg_img[2 + 5 * num_planes + np] = new double[size2];
					dbg_img[2 + 6 * num_planes + np] = attractions[np][dbg_ml];
					dbg_img[2 + 7 * num_planes + np] = new double[size2];
					dbg_img[2 + 8 * num_planes + np] = new double[size2];
					dbg_img[2 + 9 * num_planes + np] = new double[size2];
					for (int i = 0; i < size2; i++) {
						dbg_img[2 + 1 * num_planes + np][i] -= pds[np][0][i];
						dbg_img[2 + 3 * num_planes + np][i] = (num_cells[np][dbg_ml][i] > 1) ? (1.0/flatness[np][dbg_ml][i]): (1.0/disp_var_floor); // floor2);
						dbg_img[2 + 5 * num_planes + np][i] =  mask[np][i] ? 1.0 : 0.0;
						if (all_best_selections != null) {
							dbg_img[2 + 7 * num_planes + np][i] = all_best_selections[np][dbg_ml][i] ? 1.0 : 0.0;
						}
						dbg_img[2 + 8 * num_planes + np][i] = best_selections[np][dbg_ml][i] ? 1.0 : 0.0;
						dbg_img[2 + 9 * num_planes + np][i] = planes[np+1].measuredSelection[dbg_ml][i]? 1.0:0.0; // temporarily

						if (disp_strength[dbg_ml][1][i] == 0.0){
							dbg_img[2 + 1 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 2 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 3 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 5 * num_planes + np][i] = Double.NaN;
							dbg_img[2 + 6 * num_planes + np][i] = Double.NaN;
						}
					}
				}
				showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
				sdfa_instance.showArrays(dbg_img,     2 * superTileSize, 2 * superTileSize, true, "refine-"+prefix,dbg_titles);
			}
			if (debugLevel > 1)	{
				System.out.println ("refineDiscriminateTiles() "+prefix+" - step 3");
			}
			// moved here to debug
			if  (merge_planes != null) {
				if ((max_attr_corr > attractionCorrMax) && (max_attr_corr > attractionCorrMerge)) { // attractionCorrMerge can be ==0, then all >..Max will be merged
					tilePlanes.get(merge_pair[0]).getMark0();

					merge_planes[0] = tilePlanes.get(merge_pair[0]).getMark0(); // merge_pair[0];
					merge_planes[1] = tilePlanes.get(merge_pair[1]).getMark0(); // merge_pair[1];
				} else {
					merge_planes[0] = -1;
					merge_planes[1] = -1;
				}
			}

			if (((debugLevel > 0) && ((debugLevel > 1) || (max_attr_corr > attractionCorrMax)) )&& (num_planes > 1)){
				String dbg_s = "refineDiscriminateTiles() plane attraction correlation for "+prefix+": maximal="+max_attr_corr;
				for (int np = 0; np < num_planes; np++) {
					for (int np1 = np + 1; np1 < num_planes; np1++) {
						dbg_s += String.format(" %d-%d:%6.3f",np,np1,attr_corr[np][np1]);
					}
				}
				if (merge_planes[0] >= 0) {
					dbg_s += " will merge pair: ["+merge_pair[0]+", "+merge_pair[1]+"],"+
							" planes["+tilePlanes.get(merge_pair[0]).getMark0()+"] and"+
							" planes["+tilePlanes.get(merge_pair[1]).getMark0()+"] ";
				} else if (max_attr_corr > attractionCorrMax){
					dbg_s += " KEEPENG original tile selections";
				}
				System.out.println(dbg_s);
			}

			if (max_attr_corr > attractionCorrMax) {
				return null;
			}
			return best_selections;
		}


		// detect and split/eliminate planes that are use disconnected tiles
		// and are not separated by higher disparity planes
		// return planes/selections (no remove outliers!)
		boolean[][][] filterBridges(
				ArrayList<PlaneData> tilePlanes,
				int max_grow_these,
				int max_grow_far,
				int debugLevel)
		{
			double [][] pds = new double [tilePlanes.size()][];
			double [] max_disp = null;
			boolean [][] selections = new boolean[pds.length][];
			final int tsize = 4 * stSize * stSize;
			for (int np = 0; np < pds.length; np++){
				PlaneData pd = tilePlanes.get(np);
				if (pd != null){
					pds[np] = pd.getDoublePlaneDisparityStrength(
							false, // boolean   useWorld,
							null, // double [] window,
							-1, // int       dir,
							false, // boolean   use_sel,
							false, // boolean   divide_by_area,
							1.0, // double    scale_projection,
							0.0, // double    fraction_uni,
							0)[0]; // int       debugLevel)
					if (max_disp == null){
						max_disp = pds[np].clone();
					} else {
						for (int i = 0; i < max_disp.length; i++){
							max_disp[i] = Math.max(max_disp[i], pds[np][i]);
						}
					}
					boolean [][] ms = pd.getMeasSelection(); // should not be null
					// combine selections from all measurement layers
					for (int ml = 0; ml < ms.length; ml++) if (ms[ml] != null){
						if (selections[np] == null) {
							selections[np] = ms[ml].clone();
						} else {
							for (int i = 0; i < ms[ml].length; i++) {
								selections[np][i] |= ms[ml][i];
							}
						}
					}
				}
			}
			class TileSelections{
				int        np;   // number of original plane
				boolean [] mask; // selection mask to apply to measurement layers
				TileSelections(int nl, boolean [] mask){
					this.np = nl;
					this.mask = mask;
				}
			}
			ArrayList<TileSelections> split_planes = new ArrayList<TileSelections>();
			TileNeibs tileNeibs = new TileNeibs(2 * stSize, 2 * stSize);
			HashSet<Integer> old_planes = new HashSet<Integer>();
			for (int np = 0; np < pds.length; np++){
				PlaneData pd = tilePlanes.get(np);
				if (pd != null){
					old_planes.add(np);
					// see if all the plane tiles have other planes with higher disparity
					boolean all_far = true;
					for (int i = 0; i < max_disp.length; i++)if (selections[np][i]){
						if (pds[np][i] == max_disp[i]){
							all_far = false;
							break;
						}
					}
					if (!all_far){
						boolean [] grown_sel = selections[np].clone();
						tileNeibs.growSelection( //
								max_grow_these,  // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
								grown_sel, // boolean [] tiles,
								null);      // boolean [] prohibit,
						int [] clusters = tileNeibs.enumerateClusters(
								grown_sel, // boolean [] tiles,
								false); // boolean ordered)

						int num_clusters = tileNeibs.getMax(
								clusters); // int [] data)
						if (num_clusters > 1){
							boolean [] dbg_grown_sel =grown_sel.clone();
							// determine area that should have connections between these clusters, not blocked by known
							// tiles with lower disparity
							boolean [] bound_octo= tileNeibs.boundShape(
									selections[np], // boolean [] selection,
									true); // boolean octo)
							// now combine all "offending" selections - known tiles belonging to farther planes
							boolean [] sel_offend = new boolean[tsize];
							for (int np_other = 0; np_other < pds.length; np_other++) if ((selections[np_other] != null) && (np_other != np)){
								for (int i = 0; i < tsize; i++) if (selections[np_other][i]){
									sel_offend [i] |= (pds[np][i] > pds[np_other][i]);
								}
							}
							// now grow offending selections to fill gaps
							boolean [] dbg_sel_offend = sel_offend.clone();
							tileNeibs.growSelection( //
									max_grow_far,  // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
									sel_offend, // boolean [] tiles,
									null);      // boolean [] prohibit,
							// create prohibited mask - all not in grown selection that are outside bounding octagon or belong to grown offenders

							boolean [] prohibit = new boolean [tsize];
							for (int i = 0; i < prohibit.length; i++) {
								prohibit[i] = !grown_sel[i] && (!bound_octo[i] || sel_offend[i] ); // (pds[np][i] > max_disp[i]);
							}
							tileNeibs.growSelection( //
									stSize,     // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
									grown_sel,  // boolean [] tiles,
									prohibit);      // boolean [] prohibit,
							// re-clusterize what is remaining
							int [] dbg_clusters = clusters.clone();
							clusters = tileNeibs.enumerateClusters(
									grown_sel, // boolean [] tiles,
									false); // boolean ordered)
							num_clusters = tileNeibs.getMax(
									clusters);
							if (num_clusters > 1){
								boolean [][] cluster_sels = new boolean [num_clusters][selections[np].length];
								for (int i = 0; i < selections[np].length; i++) if (selections[np][i]){
									cluster_sels[clusters[i]-1][i] = true;
								}
								for (int nc = 0; nc < num_clusters; nc++){
									split_planes.add(new TileSelections(np, cluster_sels[nc]));
								}
								old_planes.remove(np);
							}
							if (debugLevel > 2) {
								String [] dbg_titles = {"sel","grown", "clust0","bound","offend","off_grown", "prohibit", "clusters","clus_masked"};
								double [][] dbg_img =  new double [dbg_titles.length][tsize];
								double lon = 3.0;
								for (int i = 0; i < tsize; i++) {
									dbg_img[0][i] = selections[np][i] ? lon: 0.0;
									dbg_img[1][i] = dbg_grown_sel[i] ? lon: 0.0;
									dbg_img[2][i] = dbg_clusters[i];
									dbg_img[3][i] = bound_octo[i] ? lon: 0.0;
									dbg_img[4][i] = dbg_sel_offend[i] ? lon: 0.0;
									dbg_img[5][i] = sel_offend[i] ? lon: 0.0;
									dbg_img[6][i] = prohibit[i] ? lon: 0.0;
									dbg_img[7][i] = clusters[i];
									dbg_img[8][i] = selections[np][i] ? clusters[i] : 0.0;
								}

								showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
								sdfa_instance.showArrays(dbg_img, 2 * superTileSize, 2* superTileSize, true, "bridges-"+np+"-"+debugLevel,dbg_titles);
								sdfa_instance.showArrays(pds,     2 * superTileSize, 2* superTileSize, true, "pds-bridges-"+np+"-"+debugLevel,dbg_titles);
							}
						}
					}
				}
			}
			if (!split_planes.isEmpty()){
				boolean [][][] split_selections = new boolean [old_planes.size() + split_planes.size()][][];
				int ns = 0;
				for (Integer np: old_planes){
					split_selections[ns++] = tilePlanes.get(np).getMeasSelection();
				}
				for (TileSelections ts: split_planes){
//					PlaneData pd = tilePlanes.get(ts.np);
					boolean [][] ms = tilePlanes.get(ts.np).getMeasSelection().clone();
					for (int ml = 0; ml < ms.length; ml++) if (ms[ml] != null){
						ms[ml] = ms[ml].clone();
						for (int i = 0; i < ms[ml].length; i++){
							ms[ml][i] &= ts.mask[i];
						}
					}
					split_selections[ns++] = ms;
				}
				return split_selections;

			}
			return null; // no changes
		}

		public double [] getWindow (
				int size)
		{
			double [] wnd1d = new double [size];
			for (int i = 0; i < size/2; i++){
				wnd1d[i] = 0.5 * (1.0 - Math.cos(2*Math.PI*(i+0.5)/size));
				wnd1d[size - i -1] = wnd1d[i];
			}
			double [] wnd = new double [size * size];
			int indx = 0;
			for (int i = 0; i < size; i++){
				for (int j = 0; j < size; j++){
					wnd[indx++] = wnd1d[i]*wnd1d[j];
				}
			}
			return wnd;
		}

		public boolean isWeakForeground(
				final PlaneData fg_plane,
//				final int stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
				final int outliers,  // remove any glitches?
				double min_fg_strength,
				final int debugLevel)
		{
			double [][] bg_ds = getDoublePlaneDisparityStrength(
					false, // boolean   useWorld,
					null, // double [] window,
					-1, // int       dir,
					false, // boolean   use_sel,
					false, // boolean   divide_by_area,
					1.0, // double    scale_projection,
					0.0, // double    fraction_uni,
					0); // int       debugLevel)
			double [][] fg_ds = fg_plane.getDoublePlaneDisparityStrength(
					false, // boolean   useWorld,
					null, // double [] window,
					-1, // int       dir,
					false, // boolean   use_sel,
					false, // boolean   divide_by_area,
					1.0, // double    scale_projection,
					0.0, // double    fraction_uni,
					0); // int       debugLevel)
			boolean [][] fg_sel = fg_plane.getMeasSelection();
			for (int ml = 0; ml < fg_sel.length; ml++) if (fg_sel[ml] != null){ // if ((stMeasSel & ( 1 << ml)) != 0){
				for (int indx = 0; indx < bg_ds[0].length; indx++) if (fg_sel[ml][indx]) {
					if (fg_ds[0][indx] < bg_ds[0][indx]) {
						return false;  // not a foreground
					}
				}
			}

			// it is foreground, now get measured data and find maximal strength (remove outlayers?
			double [] lap_weights = measuredLayers.getLapWeights1d();
			double [][][] disp_strength = new double[measuredLayers.getNumLayers()][][];
			for (int ml = 0; ml < disp_strength.length; ml++) if (fg_sel[ml] != null){ //  if ((stMeasSel & ( 1 << ml)) != 0) {
				if (this.smplMode) {
					disp_strength[ml] = measuredLayers.getDisparityStrengthMLTilted( // expensive to calculate (improve removing outlayers
							ml, // int num_layer,
							getSTileXY()[0],        // int stX,
							getSTileXY()[1],        // int stY,
							null,                   // boolean [] sel_in, - use all
							this.mlfp,
//							strength_floor,
//							measured_strength_pow,  //
//							smplSide,               // = 2;      // Sample size (side of a square)
//							smplNum,                // = 3;      // Number after removing worst
//							smplRms,                // = 0.1;    // Maximal RMS of the remaining tiles in a sample
//							smplWnd,                // use window functions for the samples
//							max_abs_tilt,           // 2.0; // Maximal absolute tilt in pixels/tile
//							max_rel_tilt,           // 0.2; // Maximal relative tilt in pixels/tile/disparity
//							damp_tilt,              //    0.001; // Damp tilt to handle insufficient  (co-linear)data
//							min_tilt_disp,          // 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//							transition,             // 1.0; // Mode transition range (between tilted and maximal disparity)
//							far_mode,               //     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
//							far_power,              //    1.0; // Raise disparity to this power before averaging for far objects
							true,                   // boolean null_if_none)
							debugLevel);
				} else {
					disp_strength[ml] = measuredLayers.getDisparityStrengthML(
							ml,                     // int num_layer,
							getSTileXY()[0],        // int stX,
							getSTileXY()[1],        // int stY,
							null,                   // boolean [] sel_in, - use all
							this.mlfp,
//							strength_floor,         //  double strength_floor,
//							measured_strength_pow,  // double strength_pow,
							true);                  // boolean null_if_none);
				}
			}
			ArrayList<Double> fg_strengths = new ArrayList<Double>();
			for (int ml = 0; ml < disp_strength.length; ml++)  if (fg_sel[ml] != null){ //  if ((stMeasSel & ( 1 << ml)) != 0) {
				double [] strength = disp_strength[ml][1];
				for (int indx = 0; indx < strength.length; indx ++)  if (fg_sel[ml][indx]) {
					fg_strengths.add(strength[indx]/lap_weights[indx]); // strengths were already masked by a window
				}
			}
			if (fg_strengths.isEmpty()){
				return true; // should not happen, but it is OK to merge empty plane
			}
			Collections.sort(fg_strengths);
			int indx = fg_strengths.size() - outliers -1;
			if (indx < 0) indx = 0;
			if (debugLevel > 0) {
				if (fg_strengths.get(indx) >= min_fg_strength) {
					System.out.println("strong pair: "+fg_strengths.get(indx)+" >= "+min_fg_strength);
				} else {
					System.out.println("weak pair: "+fg_strengths.get(indx)+" < "+min_fg_strength);
				}
			}
			return (fg_strengths.get(indx) < min_fg_strength);
		}
	}
}
