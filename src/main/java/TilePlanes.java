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
		boolean []  plane_sel =  null; // tile selection - has twice supertile size in each direction
		double []   zxy =        null; // [3] - plane center point {disparity, x, y), x=0, y=0 is a 4,4 point of an 8x8 supertile (in pixels, relative to this supertile center)
		double [][] vectors =    null; // [3][3] - re-ordered/re-directed eigenvectors(transposed): [0] - plane normal, most Z-like, towards camera, [1] - X-like, [2] - Y-like
		double []   values =     null; // [3] -eigenvalues
		int         num_points = 0;
		double      weight =     0.0;
		double []   center_xyz = null; // center of this supertile this plane center in world coordinates
		double []   world_xyz =  null; // world coordinates of the nearest point of the plane, in meters
		double []   world_v1 =   null; // world in-plane vector, corresponding to vectors[1] 
		double []   world_v2 =   null; // world in-plane vector, corresponding to vectors[1]
//		double []   daxy      =  null; // disparity and 2 relative angles (ax and ay) corresponding to fisheye view, near (0,0) scale is pixel size
		double [][] merged_eig_val = null; // for each of the directions (N, NE, .. NW) quality match for each layer 
		int    []   neib_best =  null; // new int [8]; // for each of the directions (N, NE, .. NW) index of best match, -1 if none 
// stores "worsening" of merging 2 	planes. if L1,L2,L = values[0] of plane1, plane2 plane composite: w1, w2 - weights for plane1, plane2
//		Lav = Math.sqrt((L1 * L1 * w1 + L2 * L2 * w2)/(w1 + w2))
// worsening_12 = (L - Lav) * (w1 + w2) * (w1 + w2) / (Lav * x1 * w2)		
		
		
		
		
		int         tileSize;
		int         superTileSize;
		int []      sTileXY =    null; // X and Y indices of this superTile in the image 
		
		public PlaneData clone(){
			PlaneData pd = new PlaneData(
					this.sTileXY,
					this.tileSize,
					this.superTileSize,
					this.geometryCorrection);
			pd.num_points = this.num_points; 
			pd.weight =     this.weight; 
			if (this.plane_sel != null)  pd.plane_sel =  this.plane_sel.clone(); 
			if (this.zxy != null)        pd.zxy =        this.zxy.clone(); 
			if (this.values != null)     pd.values =     this.values.clone(); 
			if (this.center_xyz != null) pd.center_xyz = this.center_xyz.clone(); 
			if (this.world_xyz != null)  pd.world_xyz =  this.world_xyz.clone(); 
			if (this.world_v1 != null)   pd.world_v1 =   this.world_v1.clone(); 
			if (this.world_v2 != null)   pd.world_v2 =   this.world_v2.clone(); 
//			if (this.daxy != null)       pd.daxy =       this.daxy.clone(); 
			if (this.vectors != null) {
				pd.vectors = new double[3][];
				pd.vectors[0] = this.vectors[0].clone();
				pd.vectors[1] = this.vectors[1].clone();
				pd.vectors[2] = this.vectors[2].clone();
			}
			copyNeib(this,pd);
			/*
			if (this.merged_eig_val != null){
				pd.merged_eig_val = this.merged_eig_val.clone();
				for (int i = 0; i<this.merged_eig_val.length; i++){
					if (this.merged_eig_val[i] != null){
						pd.merged_eig_val[i] = this.merged_eig_val[i].clone();
					}
				}
			}
			if (this.neib_best != null) pd.neib_best = this.neib_best.clone();
			*/
			return pd;
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
				GeometryCorrection   geometryCorrection)
		{
			this.geometryCorrection = geometryCorrection;
			this.tileSize = tileSize;
			this.superTileSize = superTileSize;
			this.sTileXY = sTileXY.clone();
		}
		
		public double [][] initMergedValue()
		{
			this.merged_eig_val = new double[8][];
			return this.merged_eig_val;
		}
		public double [][] getMergedValue()
		{
			return this.merged_eig_val;
		}
		public double [] initMergedValue(int dir, int leng)
		{
			this.merged_eig_val[dir] = new double[leng];
			for (int i = 0; i < leng; i++) this.merged_eig_val[dir][i] = Double.NaN;
			return getMergedValue(dir);
		}

		public double [] getMergedValue(int dir)
		{
			if (this.merged_eig_val == null) {
				return null;
			}
			return this.merged_eig_val[dir];
		}

		public double getMergedValue(int dir, int plane)
		{
			if ((this.merged_eig_val == null) ||(this.merged_eig_val[dir] == null)){
				return Double.NaN;
			}
			return this.merged_eig_val[dir][plane];
		}
		public void setNeibMatch(int dir, int plane, double value)
		{
			this.merged_eig_val[dir][plane] = value;
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

		public int getNeibBest(int dir)
		{
			if (this.neib_best == null) {
				return -1;
			}
			return this.neib_best[dir];
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
		public double[] getPlaneDisparity(
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

		/**
		 * Get disparity values for the tiles of this overlapping supertile as [2*superTileSize * 2*superTileSize] array
		 * @param useNaN replace unselected tiles with Double.NaN
		 * @return array of disparity values for the plane (not including overlapped areas)
		 */
		public double[] getDoublePlaneDisparity(
				boolean useNaN)
		{
			double [] disparities = new double[4*superTileSize*superTileSize];
			int indx = 0;
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			for (int sy = -superTileSize; sy < superTileSize; sy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				for (int sx = -superTileSize; sx < superTileSize; sx++){
					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
					if (plane_sel[indx] || !useNaN ||  (plane_sel==null)){
						disparities[indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
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
			double [] plane_all = getDoublePlaneDisparity(useNaN);
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
			double [] disparities = new double[9*superTileSize*superTileSize];
			int indx = 0;
			double [] normal = getVector();
			double [] zxy =    getZxy(); // {disparity, x center in pixels, y center in pixels (relative to a supertile center)
			for (int sy = -3 * superTileSize / 2; sy < 3* superTileSize / 2; sy++){
				// adding half-tile and half-pixel to match the center of the pixel. Supertile center is between
				// pixel 31 and pixel 32 (counting from 0) in both directions
				double y = tileSize * (sy + 0.5) + 0.5  - zxy[2];
				for (int sx = -3 * superTileSize/2; sx < 3 * superTileSize / 2; sx++){
					double x = tileSize * (sx + 0.5) + 0.5 - zxy[1];
						disparities[indx] = zxy[0] - (normal[1] * x + normal[2] * y)/normal[0];
					indx++;
				}
			}
			return disparities;
		}
		
		public double[] getTriplePlaneDisparity(
				int     dir)
		{
			double [] plane_all = getTriplePlaneDisparity();
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
		
		
//		double px = tileSize*(superTileSize * sTileXY[0] + superTileSize/2) + zxy[1];  // [3] - plane point {disparity, x, y), x=0, y=0 is a 4,4 point of an 8x8 supertile
//		double py = tileSize*(superTileSize * sTileXY[1] + superTileSize/2) + zxy[2];
		
		
		
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
				boolean   ignore_weights,
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
			double sum_weight =      scale_other *  otherPd.weight + this.weight;
			double other_fraction = ignore_weights? (scale_other/(scale_other + 1.0)): ((scale_other *  otherPd.weight) / sum_weight); 
			Matrix common_center =  this_center.times(1.0 - other_fraction).plus(other_center.times(other_fraction));
			Matrix other_offset =   other_center.minus(this_center); // other center from this center 
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
			
			// extract new eigenvalues, eigenvectors
			EigenvalueDecomposition eig = covar.eig();
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
			int oindx = 0;
			for (int i = 1; i <3; i++){
				if (Math.abs(eig_vect[0][i]) > Math.abs(eig_vect[0][oindx])){
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
			pd.setWeight(other_fraction * otherPd.weight + (1.0 - other_fraction) * this.weight);
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

		public double [] getWorldXYZ(
				boolean correct_distortions,
				int debugLevel)
		{
			double delta = 0.0001;
			if (world_xyz != null) return world_xyz;
			setCorrectDistortions(correct_distortions);
			// get pixel coordinates of the plane origin point
//			double px = tileSize*(superTileSize * sTileXY[0] + superTileSize/2) + zxy[1];  // [3] - plane point {disparity, x, y), x=0, y=0 is a 4,4 point of an 8x8 supertile
//			double py = tileSize*(superTileSize * sTileXY[1] + superTileSize/2) + zxy[2];
			double [] px_py = getCenterPxPy();
//			if ((px_py[0] == 1760) && (px_py[1] == 1056)){ // 27, 15
//				System.out.println("getWorldXYZ, px_py = {"+px_py[0]+","+px_py[1]+"}");
//				debugLevel = 2;
//			}
			
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
//				if (sTileXY[0] == 27 ){
//					System.out.println("STOP");
//					
//				}
				System.out.println("getWorldXYZ("+sTileXY[0]+","+sTileXY[1]+"), correctDistortions="+correctDistortions+", xyz= {"+
						xyz.get(0, 0)+","+xyz.get(1, 0)+","+xyz.get(2, 0)+"}, weight = "+getWeight());
//				xyz.print(10, 6); // w,d
				
//				double [] dpxpy = geometryCorrection.getImageCoordinates(xyz.getColumnPackedCopy(),this.correctDistortions);
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
//			norm_xyz = norm_xyz.times(1.0/norm_xyz.normF()); // unity normal vector;
//			norm_xyz = jacobian.times(dpxpy); // plane normal vector in world xyz
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
//				Matrix wn = norm_xyz.times(-dotprod);
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
//			world_xyz = norm_xyz.times(-(xyz.transpose().times(norm_xyz).get(0,0))).getColumnPackedCopy();
			world_xyz = norm_xyz.times((xyz.transpose().times(norm_xyz).get(0,0))).getColumnPackedCopy();
			return world_xyz;
		}
	}
	
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
		
		if (debugLevel > 0){
			System.out.println("getCovar(): sw = "+sw +", swz = "+swz +", swx = "+swx +", swy = "+swy);
		}
		
		// TODO: scale disparity to make same scale for 3 axes?
		
		double kz = ((plDispNorm > 0.0) && (swz > plDispNorm)) ? (plDispNorm / swz) : 1.0; 
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
		EigenvalueDecomposition eig = covar.eig();
		if (debugLevel > 0){
			System.out.println("getCovar(): covarianvce matrix, number of used points:"+numPoints);
			covar.print(10, 6); // w,d
			System.out.println("getCovar(): eigenvalues");
			eig.getD().print(10, 6); // w,d
			System.out.println("getCovar(): eigenvectors");
			eig.getV().print(10, 6); // w,d
		}
		double [][][] rslt = {
				eig.getD().getArray(),
				eig.getV().getArray(),
				{
					{sw,kz,numPoints},
					{swz, swx, swy}}};
		return rslt;
	}
	
	public PlaneData getPlane(
			int [] sTileXY,
			double []  data,
			double []  weight,
			boolean [] select, // null OK, will enable all tiles
			int        debugLevel){
		if (select == null) {
			select = new boolean [4*stSize];
			for (int i = 0; i < select.length; i++) select[i] = true;
		}
//		int debugLevel1 = ((sTileXY[0] == 27) && (sTileXY[1] == 20))? 1: 0; 
//		int debugLevel1 = ((sTileXY[0] == 27) && (sTileXY[1] == 17))? 1: 0; 
		int debugLevel1 = ((sTileXY[0] == 27) && (sTileXY[1] == 16))? 1: 0; // check why v[0][0] <0  
		
		
		double [][][] rslt = getCovar(
				data,
				weight,
				select,
				0.0,
				debugLevel1); //0); // debugLevel);
		if (rslt == null) return null;
		int       numPoints =  (int) rslt[2][0][2];
//		double    kz =   rslt[2][0][1]; // == 1.0
		double    swc =  rslt[2][0][0];
		double [] szxy = rslt[2][1];
		double [][] eig_val =  rslt[0];
		double [][] eig_vect = rslt[1];
		// find vector most orthogonal to view // (anyway it all works with that assumption), make it first
		// TODO normalize to local linear scales
		int oindx = 0;
		for (int i = 1; i <3; i++){
			if (Math.abs(eig_vect[0][i]) > Math.abs(eig_vect[0][oindx])){
				oindx = i;
			}
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
				this.geometryCorrection);
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
			double     targetEigen, // target eigenvalue for primary axis (is disparity-dependent, so is non-constant)
			int        maxRemoved,  // maximal number of tiles to remove (not a constant)
			int        minLeft,     // minimal number of tiles to keep
			int        debugLevel){
		int stSize2 = 2 * stSize;
		if (pd == null) {
			pd = getPlane(
					sTileXY, 
					data,
					weight,
					select, // null OK
					debugLevel);
		} else if (select != null){
			pd.setPlaneSelection(select);
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
			pd = getPlane(    // re-calculate with point removed
					pd.getSTileXY(),
					data,
					weight,
					select,
					debugLevel);
		}
		return pd;
	}
	
}
