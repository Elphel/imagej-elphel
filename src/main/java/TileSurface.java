/**
 **
 ** TileSurface - handle tile surfaces
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  TileSurface.java is free software: you can redistribute it and/or modify
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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;


public class TileSurface {
//	public 
//		private int tileSize;
		private int superTileSize;
		private int imageTilesX;
		private int imageTilesY;
		private int stilesX;
		private int stilesY;
		private int [] st_dirs8;
		private int [] t_dirs8;
		private double [] window;
		private int threadsMax = 100;
		private int [][] tileLayers = null;
		private TileData [][] tileData = null;
		private SurfaceData [] surfaceData = null; // result surfaces
		
		static int STAT_UNASSIGNED = 0; // index of number of unassigned tiles
		static int STAT_ASSIGNED =   1; // index of number of assigned tiles
		static int STAT_PROHIBITED = 2; // index of number of initially prohibited tiles
		static int STAT_IMPOSSIBLE = 3; // index of number of impossible (like no surfaces at that location) tiles
		static int STAT_NUM_ML =     4; // index of number of measurement layers used
		static int STAT_LEN =        5; // number of stat entries
		
		static int UNASSIGNED =   0; //tile marked as invalid
		static int PROHOBITED =  -1; //tile marked as invalid
		static int IMPOSSIBLE =  -2; // give up on this tile (as no surface for it) 
		
		static int NEW_ASSIGNED =     0; // successfully assigned to a surface
		static int NO_SURF =          1; // no surfaces for this tile cell
		static int TOO_WEAK =         2; // tile strength is too low
		static int TOO_STRONG =       3; // tile strength is too high ( for that disparity difference)
		static int TOO_FAR =          4; // no surface candidates within the allowed disparity range 
		static int NOT_UNIQUE =       5; // multiple surfaces are within range
		static int REMOVED_TILES =    6; // number of removed tiles in weak clusters
		static int REMOVED_CLUSTERS = 7; // number of removed weak clusters
		static int NUM_STATS =        8;

		static int CLUST_NUM_INDEX =       0;
		static int CLUST_NUM_TILES =       1;
		static int CLUST_NUM_LAYERS =      2;
		static int CLUST_NUM_OVERLAPS =    3;
		static int CLUST_NUM_MULTI =       4;
		static int CLUST_NUM_CONFLICTS_A = 5;
		static int CLUST_NUM_CONFLICTS_B = 6;
		static int CLUST_NUM_STATS =       7;
		
		
		
//		private int nsTilesstSize =   0; // 8;
		GeometryCorrection   geometryCorrection = null;
		public TileSurface(
				int tileSize,
				int superTileSize,
				int tilesX,
				int tilesY,
				GeometryCorrection geometryCorrection,
				int threadsMax){
//			this.tileSize = tileSize;
			this.superTileSize = superTileSize;
			this.geometryCorrection =geometryCorrection;
			this.imageTilesX =  tilesX;
			this.imageTilesY =  tilesY;
			this.window = getWindow(2*superTileSize);
			this.threadsMax = threadsMax;
			stilesX = (tilesX + superTileSize -1)/superTileSize;
			stilesY = (tilesY + superTileSize -1)/superTileSize;
//			int [] dirs =  {-tilesX, -tilesX + 1, 1, tilesX + 1, tilesX, tilesX - 1, -1, -tilesX - 1};
			int [] dirs =  {-stilesX, -stilesX + 1, 1, stilesX + 1, stilesX, stilesX - 1, -1, -stilesX - 1};
			this.st_dirs8 = dirs;

			int tx = superTileSize * stilesX; 
			int [] tdirs =  {-tx, -tx + 1, 1, tx + 1, tx, tx - 1, -1, -tx - 1};
			this.t_dirs8 = tdirs;
			
		}

		public class TileData{
			double [] disp_strength;
			int indx =         0;
			int new_index =    0;
			boolean enable =   true;
			int [] neighbors = {-1,-1,-1,-1,-1,-1,-1,-1};
			int dbg_nsTile;
			public TileData (
					double disparity,
					double strength)
			{
				setDisparityStrength(disparity,strength);
			}
			
			public void setDbgNsTile(int dbg_nsTile)
			{
				this.dbg_nsTile = dbg_nsTile;
			}

			public int getDbgNsTile()
			{
				return this.dbg_nsTile;
			}

			public void setIndex(int indx)
			{
				this.indx = indx;
				if (indx < 0){
					System.out.println("setIndex("+indx+")");
				}
			}
			public int  getIndex()
			{
				return this.indx; 
			}


			public void setNewIndex(int indx)
			{
				this.new_index = indx; 
			}
			public int  getNewIndex()
			{
				return this.new_index; 
			}
			
			
			public void setNeighbors(int [] neighbors)
			{
				this.neighbors = neighbors; 
			}
			public int [] getNeighbors()
			{
				return this.neighbors; 
			}
			public void setNeighbor(int dir,int neib)
			{
//				if (this.neighbors == null) this.neighbors = new int[8];
				this.neighbors[dir] = neib; 
			}
			public int getNeighbor(int dir)
			{
//				if (this.neighbors == null) this.neighbors = new int[8];
				return this.neighbors[dir]; 
			}
			public void setEnable(boolean enable)
			{
				this.enable = enable;
			}
			public boolean getEnable()
			{
				return this.enable;
			}
			public void setDisparityStrength(
					double disparity,
					double strength)

			{
				this.disp_strength = new double[2];
				this.disp_strength[0] = disparity;
				this.disp_strength[1] =  strength;
			}
			
			public void setDisparity(double disparity)
			{
				if (this.disp_strength == null){
					this.disp_strength = new double[2];
				}
				this.disp_strength[0] = disparity;
			}
			public double getDisparity()
			{
				if (this.disp_strength == null){
					this.disp_strength = new double[2];
				}
				return this.disp_strength[0];
			}
			public double getDisparity(boolean useNaN)
			{
				if (useNaN && (this.disp_strength == null)) return Double.NaN;
				if (this.disp_strength == null){
					this.disp_strength = new double[2];
				}
				if (useNaN && (this.disp_strength[1] == 0.0))  return Double.NaN;
				return this.disp_strength[0];
			}

			public double getDisparityNaN()
			{
				return getDisparity(true);
			}
			
			public void setStrength(double strength)
			{
				if (this.disp_strength == null){
					this.disp_strength = new double[2];
				}
				this.disp_strength[1] = strength;
			}
			public double getStrength()
			{
				if (this.disp_strength == null){
					this.disp_strength = new double[2];
				}
				return this.disp_strength[1];
			}
			
		}

		public int getNStileDir(
				int nsTile,
				int dir)
		{
			if (dir < 0) return nsTile;
			int sty = nsTile / stilesX;
			int stx = nsTile % stilesX;
			if ((stx > 0) && (sty > 0) && (sty < (stilesY - 1)) && (stx < (stilesX - 1))) return nsTile + st_dirs8[dir]; // most likely case
			if ((sty == 0)             && ((dir < 2) || (dir == 7))) return -1; 
			if ((sty == (stilesY - 1)) &&  (dir > 2) && (dir < 6))   return -1; 
			if ((stx == 0)            &&  (dir > 4))                 return -1; 
			if ((stx == (stilesX - 1)) &&   (dir > 0) && (dir < 4))  return -1;
			return nsTile + st_dirs8[dir];
		}

		public int getNtileDir(
				int nTile,
				int dir)
		{
			if (dir < 0) return nTile;
			int tilesX = stilesX * superTileSize;
			int tilesY = stilesY * superTileSize;
			int ty = nTile / tilesX;
			int tx = nTile % tilesX;
			if ((tx > 0) && (ty > 0) && (ty < (tilesY - 1)) && (tx < (tilesX - 1))) return nTile + t_dirs8[dir]; // most likely case
			if ((ty == 0)            && ((dir < 2) || (dir == 7))) return -1; 
			if ((ty == (tilesY - 1)) &&  (dir > 2) && (dir < 6))   return -1; 
			if ((tx == 0)            &&  (dir > 4))                 return -1; 
			if ((tx == (tilesX - 1)) &&  (dir > 0) && (dir < 4))  return -1;
			return nTile + t_dirs8[dir];
		}
		
		
		public int getDirToStile(
				int nsTile,
				int nsTile1)
		{
			
			int sty =  nsTile / stilesX;
			int stx =  nsTile % stilesX;
			int sty1 = nsTile1 / stilesX;
			int stx1 = nsTile1 % stilesX;
			int dx = stx1 - stx; 
			int dy = sty1 - sty;
//			int sdx = (dx > 0) ? 1: ( (dx < 0) ? -1 : 0);
//			int sdy = (dy > 0) ? 1: ( (dy < 0) ? -1 : 0);
			if ((dy == 0 ) && (dx == 0)) return -1; // same tile
			if (dy < 0) {
				if (dx < 0) return 7;
				if (dx > 0) return 1;
				return 0;
			}
			if (dy > 0) {
				if (dx < 0) return 5;
				if (dx > 0) return 3;
				return 4;
			}
			if (dx < 0) return 6;
			if (dx > 0) return 2;
			return -1;
		}

		
		/**
		 * Get tile surface number from supertile number, direction (-1 same) and the supertile plane index
		 * @param nsTile number of the supertile
		 * @param dir direction -1 - same supertile, 0 - N (up), 1 - NE, .. 7 - NW
		 * @param np number of the supertile plane
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @return unique tile surface index, if ((dir == 8) && (np == 0)) returns total number of tile surfaces
		 */
		public int getTileSurfaceNumber (
				int nsTile,
				int dir,              // direction, or -1 (same)
				int np,
				TilePlanes.PlaneData [][] planes)
		{
			if (nsTile < 0) {
				return -1;
			}
			int tsn = (planes[nsTile] == null) ? 0 : planes[nsTile].length; // nsTile = -1 when mesh goes out of the image area
			if (dir < 0) {
				if (np >= tsn){
					return -1;
				}
				return np; 
			}
			int nsTile1 = -1;
			for (int d = 0; d < dir; d ++){
				nsTile1 = getNStileDir(nsTile, d);
				if ((nsTile1 >=0) && (planes[nsTile1] != null)){
					tsn += planes[nsTile1].length;
				}
			}
			if (dir < 8) {
				nsTile1 = getNStileDir(nsTile, dir);
				int last_Length = ((nsTile1< 0) || (planes[nsTile1] == null)) ? 0: planes[nsTile1].length;
				if (np >= last_Length) {
					return -1;
				}
			}
			return tsn + np; 
		}
		/**
		 * Get supertile direction and the plane number that contributeted to a specific tile surface 
		 * @param nsTile supertile index
		 * @param tp tile surface index (generated by getTileSurfaceNumber)
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @return a pair of {dir, plane index}. dir is -1 for the plane in the same supertile, 0..7 for neighbors
		 */

		public int []  getSuperTileDirPlane (
				int nsTile,
				int tp,
				TilePlanes.PlaneData [][] planes)
		{
			int num_planes = (planes[nsTile] == null)? 0: planes[nsTile].length;
			int [] rslt = {-1, tp};
			if (tp < num_planes) return rslt;
			tp -= num_planes;
			for (int d = 0; d < st_dirs8.length; d ++){
				int nsTile1 = getNStileDir(nsTile, d);
				num_planes = ((nsTile1 >=0) && (planes[nsTile1] != null))? planes[nsTile1].length : 0;
				if (tp < num_planes){
					rslt[0] = d;
					rslt[1] = tp;
					return rslt;
				}
				tp -= num_planes; 
			}
			return null; // error - invalid input
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
		
		public double [] getWindow()
		{
			return window;
		}

		/**
		 *  TODO: replace misplaced description Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid
		 * @param fraction_uni add fraction of the total weight to each tile
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @param debugLevel debug level
		 * @param dbg_X debug supertile X coordinate
		 * @param dbg_Y debug supertile Y coordinate
		 * @return per-tile (rounded up to contain whole supertiles) sparse array of TileData instances
		 */
		public double [][][][] fuseSupertilePlanes (
				final boolean                   use_sel,
				final boolean                   divide_by_area,
				final double                    scale_projection,
				final double                    fraction_uni,
				final TilePlanes.PlaneData [][] planes, // why plane here has weighth == 0.0 ?
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final double [][][][] fused_data = new double [nStiles][][][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int dbg_tile = dbg_Y * stilesX + dbg_X;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								int dl = ((debugLevel > -1) && (nsTile == dbg_tile)) ? 3:0;
								if (dl > 0){
									System.out.println("fuseSupertilePlanes(), nsTile = "+nsTile);
								}
								double [][][] disp_strength = new double [planes[nsTile].length][][];
								for (int np = 0; np < disp_strength.length; np++){
									if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
										disp_strength[np] = planes[nsTile][np].getDoublePlaneDisparityStrength(
												getWindow(),      // double [] window,
												-1,              // int dir (-1 - center, 0- N, 1 - NE, .. 7 - NW
												false,           // use_sel,          // boolean   use_sel, apply selection to the result
												divide_by_area,   //boolean   divide_by_area,
												scale_projection, // double    scale_projection,
												fraction_uni, //				double    fraction_uni,
												debugLevel-1);   // int       debugLevel)
										// multiply disparities by strengths to calculate weighted averages
										double [] disp0 = disp_strength[np][0].clone(); // to use where strength == 0 
										for (int i = 0; i < disp_strength[np][1].length; i++){
											disp_strength[np][0][i] *= disp_strength[np][1][i]; 
										}
										//									}
										for (int dir = 0; dir < st_dirs8.length; dir++){
											int sNeib = planes[nsTile][np].getNeibBest(dir);
											if (sNeib >= 0){
												int nsTile1 = getNStileDir(nsTile, dir); //   nsTile + ((dir < 0) ? 0: st_dirs8[dir]);
												if ((nsTile1 >= 0) && (planes[nsTile1] != null)){
//													double [][] ds = planes[nsTile1][np].getSinglePlaneDisparityStrength(
													double [][] ds = planes[nsTile1][sNeib].getDoublePlaneDisparityStrength(
															getWindow(),      // double [] window,
															dir,              // int dir (-1 - center, 0- N, 1 - NE, .. 7 - NW
															false,           // use_sel,          // boolean   use_sel, apply selection to the result
															divide_by_area,   //boolean   divide_by_area,
															scale_projection, // double    scale_projection,
															fraction_uni, //				double    fraction_uni,
															debugLevel-1);   // int       debugLevel)
													for (int i = 0; i < disp_strength[np][1].length; i++){
														if (ds[1][i] > 0.0){
															disp_strength[np][1][i] += ds[1][i]; 
															disp_strength[np][0][i] += ds[1][i] * ds[0][i]; 
														}
													}
												}												
											}
										}
										// calculate weighted average for each tile
										for (int i = 0; i < disp_strength[np][1].length; i++){ 
											if (disp_strength[np][1][i] > 0.0){
												disp_strength[np][0][i] /= disp_strength[np][1][i]; 
											} else {
												disp_strength[np][0][i] = disp0[i]; // may have discontinuity
											}
										}
										if (use_sel){ // zero out selection after averaging, apply to this tile
											boolean [] sel = planes[nsTile][np].getSelMask();
											if (sel != null){
												for (int i = 0; i < disp_strength[np][1].length; i++){
													if (!sel[i]) disp_strength[np][1][i] = 0.0;
												}												
											}
										}
										
										if ((debugLevel > -1) && (dl>0)){
											String str_neib =  "fuseSupertilePlanes_"+nsTile+":"+np;
											for (int dir = 0; dir < 8; dir++){
												str_neib += " " + planes[nsTile][np].getNeibBest(dir);
											}
											System.out.println(str_neib);
										}
									}
								}
								fused_data[nsTile] = disp_strength;
								if ((debugLevel > -1) && (dl>0)){
									String[] titles = new String [3 * disp_strength.length];
									double [][] dbg_img = new double [titles.length][];
									for (int i = 0; i < disp_strength.length; i++) {
										titles [i + 0 * disp_strength.length] = "disp_" + i;
										titles [i + 1 * disp_strength.length] = "mdisp_" + i;
										titles [i + 2 * disp_strength.length] = "str_" + i;
										if (disp_strength[i] != null) {
											dbg_img[i + 0 * disp_strength.length] = disp_strength[i][0]; 
											dbg_img[i + 2 * disp_strength.length] = disp_strength[i][1]; 
											dbg_img[i + 1 * disp_strength.length] = disp_strength[i][0].clone();
											for (int j = 0; j < disp_strength[i][0].length; j++){
												if (disp_strength[i][1][j] == 0.0){
													dbg_img[i + 1 * disp_strength.length][j] = Double.NaN;
												}
											}
										}
									}
									showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
									sdfa_instance.showArrays(dbg_img,  2 * superTileSize, 2 * superTileSize, true, "surf_ds_"+nsTile, titles);
								}
								
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return fused_data;
		}


		/**
		 * Prepare topology of the supertiles connections. For each of the 4 quadrants (0, 1 / 2, 3) of each
		 * used supertile plane, get 4 plane indices that contribute to it (also in linescan (0, 1/ 2,3) order
		 * That will tell which of the overlapping 2x supertile planes can be merged
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @param debugLevel debug level
		 * @param dbg_X debug supertile X coordinate
		 * @param dbg_Y debug supertile Y coordinate
		 * @return per-tile , per plane, per quadrant, per quadrant 4 corners - index of contributing plane (or -1 if none)
		 */
		public int [][][][] getSupertilesTopology (
				final TilePlanes.PlaneData [][] planes,
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
//			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int [][] dir_corn = {
					{ 7,  0,  6, -1},  // 0 (top left)
					{ 0,  1, -1,  2},  // 1 (top right)
					{ 6, -1,  5,  4},  // 2 (bottom left)
					{-1,  2,  4,  3}}; // 3 (bottom right)
			
			final int [][][][] corners = new int [nStiles][][][];
			final int dbg_tile = dbg_Y * stilesX + dbg_X;

			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								int dl = ((debugLevel > -1) && (nsTile == dbg_tile)) ? 3:0;
								if (dl > 0){
									System.out.println("getSupertilesTopology(), nsTile = "+nsTile);
								}
								corners[nsTile] = new int [planes[nsTile].length][][];
								for (int np = 0; np < planes[nsTile].length; np++){
//									if (planes[nsTile][np] != null){
									if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
										int [] neibs = planes[nsTile][np].getNeibBest();
										if (neibs == null) {
											System.out.println("getSupertilesTopology(), nsTile = "+nsTile+" neibs= null");
											
										} else {
											corners[nsTile][np]= new int [4][4];
											for (int i= 0; i < 4; i++){
												for (int j = 0; j < 4; j++){
													if (dir_corn[i][j] < 0){
														corners[nsTile][np][i][j] = np;
													} else {
														corners[nsTile][np][i][j] = neibs[dir_corn[i][j]]; 
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
			return corners;
		}
		
		/**
		 * Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @param corners - topology data generated by getSupertilesTopology() method
		 * @param debugLevel debug level
		 * @param dbg_X debug supertile X coordinate
		 * @param dbg_Y debug supertile Y coordinate
		 * @return per-tile (rounded up to contain whole supertiles) sparse array of TileData instances
		 */
		public int [][][][][] generateOverlappingMeshes (
				final TilePlanes.PlaneData [][] planes,
				final int [][][][]              corners, 
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
//			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final int [][][][][] meshes = new int [nStiles][][][][];
			
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int ss2 = 2 * superTileSize;
			final int ss1 =  superTileSize;
			final int sh  =  superTileSize/2;
			final int len_st2 = ss2  * ss2 ;
			
			final int [][][] quad_check = { // [quadrant 01/23][dir: left, right, diagonal]{dir, quadrant}
					{ // top left quadrant
						{6, 1},  //left
						{0, 2},  //right
						{7, 3}   //diagonal
					},
					{ // top right quadrant
						{0, 3},  //left
						{2, 0},  //right
						{1, 2}   //diagonal
					},
					{ // bottom left quadrant
						{4, 0},  //left
						{6, 3},  //right
						{5, 1}   //diagonal
					},
					{ // bottom right quadrant
						{2, 2},  //left
						{4, 1},  //right
						{3, 0}}};//diagonal
			final int [][][] cut_ortho = { // [quadrant][left, right][index, width, height}
					{ // quadrant 0 - top left
						{0, sh,  ss1}, // left path
						{0, ss1, sh }  // right path
					},
					{ // quadrant 1 - top right
						{ss1,    ss1, sh},
						{3 * sh, sh,  ss1}
					},
					{ // quadrant 2 - bottom left
						{3 * ss2 * sh, ss1, sh},
						{ss1 * ss2, sh, ss1}
					},
					{ // quadrant 3 - bottom right
						{ss1 * ss2 + 3 * sh, sh, ss1},
						{3 * sh * ss2 + ss1, ss1, sh}
					},
			};
			
			final TileNeibs tileNeibs = new TileNeibs(2*superTileSize);

			final int dbg_tile = dbg_Y * stilesX + dbg_X;
			
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								int dl = ((debugLevel > -1) && (nsTile == dbg_tile)) ? 3:0;
								if (dl > 0){
									System.out.println("generateOverlappingMeshes(), nsTile = "+nsTile);
								}
								meshes[nsTile] = new int [planes[nsTile].length][][][];
								for (int np = 0; np < planes[nsTile].length; np++){
//									if (planes[nsTile][np] != null){
									if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
										int [][] pre_mesh = new int [len_st2][2];
										for (int i = 0; i < len_st2; i ++){
											pre_mesh[i][0] = nsTile;
											pre_mesh[i][1] = np;
										}
										int [] neibs = planes[nsTile][np].getNeibBest();
										for (int quadrant = 0; quadrant <4; quadrant ++) {
											int [] these_corner_planes = corners[nsTile][np][quadrant];
//											int [][] neib_id = new int[3][2];
											int [][] neib_id = new int[3][];
											for (int arr = 0; arr < 3; arr++){
												int dir = quad_check[quadrant][arr][0];
												if (neibs[dir] >= 0) {
													neib_id[arr] = new int [2];
													int nsTile1 = nsTile + st_dirs8[dir];
													int [] other_corner_planes = corners[nsTile1][neibs[dir]][quad_check[quadrant][arr][1]];
													neib_id[arr][0] = nsTile1;
													neib_id[arr][1] = neibs[dir];
													for (int i = 0; i < these_corner_planes.length; i++){
														if (other_corner_planes[i] != these_corner_planes[i]){
															neib_id[arr] = null;
															break;
														}
													}
												}
											}
											// depending on match values, cut and join mesh with the neighbor
											// change diagonal first (add corner square later again
											if (neib_id[2] != null){
												switch (quadrant){
												case 0: // top left 
													for (int j = 0; j < (ss1 - 1); j++){
														for (int i = ss1 - 2 - j; i>=0; i--){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 1: // top right
													for (int j = ss1; j < ss2; j++){
														for (int i = j - ss1; i >= 0; i--){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 2: // bottom left
													for (int j = 0; j < (ss1 - 1); j++){
														for (int i = ss1 + 1 + j; i < ss2; i++){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 3: // bottom right
													for (int j = ss1; j < ss2; j++){
														for (int i = ss2 + ss1 - 1 - j; i < ss2; i++){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												}
											}
											// change ortho - on top of diagonal
											for (int arr = 0; arr < 2; arr++) if (neib_id[arr] != null){
												for (int y = 0; y < cut_ortho[quadrant][arr][2]; y++){
													for (int x = 0; x < cut_ortho[quadrant][arr][1]; x++){
														int indx =  cut_ortho[quadrant][arr][0] + y * ss2 + x;
														pre_mesh[indx] = neib_id[arr]; 
													}
												}
											}
											// add corner square corner on top of possible ortho
											if (neib_id[2] != null){
												switch (quadrant){
												case 0: // top left 
													for (int j = 0; j < sh; j++){
														for (int i =  0 ; i < sh; i++){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 1: // top right
													for (int j = ss1 + sh; j < ss2; j++){
														for (int i = 0; i < sh; i++){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 2: // bottom left
													for (int j = 0; j < sh; j++){
														for (int i = ss1 + sh; i < ss2; i++){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 3: // bottom right
													for (int j = ss1 + sh; j < ss2; j++){
														for (int i = ss1 + sh; i < ss2; i++){
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												}
											}
										}
										// build mesh , then add cuts if needed
										meshes[nsTile][np] = new int [len_st2][][];
										int [][][] dbg_meshes = meshes[nsTile][np]; 
										if (dl > 0){
											System.out.println("generateOverlappingMeshes(), dbg_meshes.length = "+dbg_meshes.length);
										}
										
										for (int i = 0; i < len_st2; i ++){
											if ((pre_mesh[i] != null) && (pre_mesh[i][0] == nsTile)){
												meshes[nsTile][np][i] = new int [8][];
												for (int dir = 0; dir < 8; dir++) {
													int ineib = tileNeibs.getNeibIndex(i, dir);
													if (ineib >= 0) meshes[nsTile][np][i][dir] = pre_mesh[ineib];
												}
											}
										}
										// add cuts
										// up
										for (int ncut = 0; ncut <8; ncut++){
											int indx, dir_go = -1, dir_start = -1;
											boolean cut_right = false;
											switch (ncut){
											case 0:	dir_go = 0; dir_start =  7;	cut_right = true;  break;
											case 1:	dir_go = 0; dir_start =  0;	cut_right = false; break;
											case 2:	dir_go = 2; dir_start =  0;	cut_right = true;  break;
											case 3:	dir_go = 2; dir_start = -1;	cut_right = false; break;
											case 4:	dir_go = 4; dir_start = -1;	cut_right = true;  break;
											case 5:	dir_go = 4; dir_start =  6;	cut_right = false; break;
											case 6:	dir_go = 6; dir_start =  6;	cut_right = true;  break;
											case 7:	dir_go = 6; dir_start =  7;	cut_right = false; break;
											}
											int dir_go45 =   (dir_go + (cut_right ? 1:7)) % 8; 
											int dir_go90 =   (dir_go + (cut_right ? 2:6)) % 8; 
											int dir_go135 =  (dir_go + (cut_right ? 3:5)) % 8; 
											int dir_go180 =  (dir_go + 4) % 8; 

											indx = ss1 * (ss2 + 1); // center point
											
											for (int i = 0; i < sh; i++) indx = tileNeibs.getNeibIndex(indx, dir_go);
											if (dir_start >= 0) indx = tileNeibs.getNeibIndex(indx, dir_start);
											
											int indx1 = tileNeibs.getNeibIndex(indx, dir_go90);
//											if ((pre_mesh[indx] != null) && (pre_mesh[indx1] == null)){ // there is a cut
//											if ((pre_mesh[indx][0] == nsTile) && (pre_mesh[indx1][0] != nsTile)){ // there is a cut
											if ((meshes[nsTile][np][indx] != null) && (meshes[nsTile][np][indx1] == null)){ // there is a cut
												for (int i = 0; i < sh; i++){
													int indx_back = tileNeibs.getNeibIndex(indx, dir_go180);
													if (meshes[nsTile][np][indx_back] != null) meshes[nsTile][np][indx_back][dir_go45] = null; // NE for N, right
													if (meshes[nsTile][np][indx] != null) {
														meshes[nsTile][np][indx][dir_go90] = null; // E for N, right
														if (i > 0){
															meshes[nsTile][np][indx][dir_go135] = null; // SE for N, right
														}
													}
													indx = tileNeibs.getNeibIndex(indx, dir_go);
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
			return meshes;
		}

		/**
		 * Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @param fusedSupertilePlanes disparity/strength data generated by fuseSupertilePlanes() method
		 * @param lappingMeshes per super-tile overlapping surface meshes, generateOverlappingMeshes
		 * @param debugLevel debug level
		 * @param dbg_X debug supertile X coordinate
		 * @param dbg_Y debug supertile Y coordinate
		 * @return per-tile (rounded up to contain whole supertiles) sparse array of TileData instances
		 */
		
		public TileData [][] createTileShells (
				final TilePlanes.PlaneData [][] planes,
				final double [][][][]           fusedSupertilePlanes,
				final int [][][][][]            lappingMeshes,				
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int tilesX = stilesX * superTileSize;
			final int tilesY = stilesY * superTileSize;
			final int nTiles =  nStiles * superTileSize * superTileSize;
			final TileData [][] tile_data = new TileData [nTiles][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int ss2 = 2 * superTileSize;
			final int sh =  superTileSize/2;
			final int len2 = ss2  * ss2 ;
			final TileNeibs tileNeibs = new TileNeibs(2 * superTileSize);
			
			final int dbg_tile = dbg_Y * stilesX + dbg_X;
			
			// initialize result structure
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							int dl = ((debugLevel > -1) && (nsTile == dbg_tile)) ? 3:0;
							if (dl > 0){
								System.out.println("createTileShells():1, nsTile = "+nsTile);
							}
							int num_surf = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
									nsTile, // int nsTile,
									8, // int dir,              // direction, or -1 (same)
									0, // int np,
									planes); // TilePlanes.PlaneData [][] planes) 
							if (num_surf > 0) { // 0 - nothing in this supertile, none around - remove 
								if (num_surf > 0) { // 0 - nothing in this supertile, none around 
									int stileY = nsTile / stilesX;  
									int stileX = nsTile % stilesX;
									for (int ty = 0; ty < superTileSize; ty++){
										for (int tx = 0; tx < superTileSize; tx++){
											int indx = ((stileY * superTileSize) + ty) * tilesX + ((stileX * superTileSize) + tx); 
											tile_data[indx] = new TileData[num_surf];
										}
										
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
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								int dl = ((debugLevel > -1) && (nsTile == dbg_tile)) ? 3:0;
								if (dl > 0){
									System.out.println("createTileShells():2, nsTile = "+nsTile);
								}
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;
//								for (int np = 0; np < planes[nsTile].length; np++) if (planes[nsTile][np] != null){
								for (int np = 0; np < planes[nsTile].length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){
									
									int [][][] src_mesh = lappingMeshes[nsTile][np];
									double [][] disp_strength = fusedSupertilePlanes[nsTile][np];
									TileData [] dual_mesh = new TileData [len2]; // full overlapping dual-sized mesh
									if ((planes == null) || (planes[nsTile] == null) || (planes[nsTile][np] == null)){
										System.out.println("createTileShells():2, *** NULL here***  nsTile = "+nsTile+" np="+np);
										continue;
									}
									for (int indx = 0 ; indx < len2; indx++){
										 // src_mesh non-null elements can only be generated by this supertile, while
										 // neighbor links can point to others (connected).
										 // others (connected) should have unique surface index based on their own planes ?
										if (src_mesh[indx] != null){
											if ((dl > 0) && ((indx & 15) == 0)){
												System.out.println("createTileShells():3, nsTile = "+nsTile+", indx="+indx);
											}
											int [][] src_neibs = src_mesh[indx];
											if (src_neibs != null){
												int tsegm = tileNeibs.getSegment(indx);
												int nsTile0 = getNStileDir(nsTile, tsegm);  // supertile over which this tile is
												if (    (tsegm < 0)    || // own tile (center square)
														(nsTile0 >= 0)) { // <0 - out of picture area
													dual_mesh[indx] = new TileData( // now sets neighbors to -1
															disp_strength[0][indx],  // disparity
															disp_strength[1][indx]); // strength
													dual_mesh[indx].setDbgNsTile(nsTile);
													int dirThisfrom0 = getDirToStile(nsTile0, nsTile); // can be -1;
													int surf0 = getTileSurfaceNumber ( // Number of the surface for the tile itself
															nsTile0,      // int nsTile,
															dirThisfrom0, // int dir,              // direction, or -1 (same)
															np,           // int np,
															planes); 
													dual_mesh[indx].setIndex(surf0);
													
													for (int dir = 0; dir < 8; dir++) {
														if (src_neibs[dir] != null){
															int nsTile1 = src_neibs[dir][0];
															int np1 = src_neibs[dir][1];
															int indx1 = tileNeibs.getNeibIndex(indx,dir); // index of the destination tile
															// now find tile location - it may be outside of both nsTile and nsTile_1
															int segm1 = tileNeibs.getSegment(indx1);
															int nsTile2 = getNStileDir(nsTile, segm1);  // negative segm 1 is OK ?
															// now: nsTile - this supertile,
															// nsTile1 - to which supertile surface we switch
															// nsTile2 - non-overlapping supertile where the destination tile belongs
															// Unique surface number should be determined for nsTile2, generated by nsTile1,
															// for direction how nsTile1 is visible from the nsTile2
															int dir1from2 = getDirToStile(nsTile2, nsTile1); // can be -1;
															int surf = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
																	nsTile2,  // int nsTile,
																	dir1from2, // int dir,              // direction, or -1 (same)
																	np1, // int np,
																	planes); 
															dual_mesh[indx].setNeighbor(dir, surf);
														}
													}
												}
											}
										}
									}
									// Now we have a double-sized surface with all tiles set with correct absolute indices, now just split it
									//surf_number =
									for (int ty = 0; ty < ss2; ty++ ){
										for (int tx = 0; tx < ss2; tx++ ){
											int indx = ty * ss2 + tx;
											if (dual_mesh[indx] != null) { // some cells may be missing after merge
												int ix, iy;
												ix = (stileX * superTileSize) + tx -sh; iy = (stileY * superTileSize) + ty -sh;

												if ((ix >= 0) && (ix < tilesX) && (iy >= 0) && (iy < tilesX)) {
													int tindx = iy * tilesX + ix;
													tile_data[tindx][dual_mesh[indx].getIndex()] = dual_mesh[indx]; //oob
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
			return tile_data;
		}		

		public TileData [][]  compactSortShells (
				final TileData [][]     tileData_src,
				final int               debugLevel,
				final int               dbg_X,
				final int               dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final TileData [][] tile_data = new TileData [nTiles][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int dbg_stile = (dbg_Y * superTileSize) * (stilesX * superTileSize) + (dbg_X * superTileSize);
			final int dbg_tile = dbg_Y * stilesX + dbg_X;
			final int dbg_tilesX = stilesX * superTileSize;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							int dbg_stX = nTile % dbg_tilesX;   
							int dbg_stY = nTile / dbg_tilesX;
							int dbg_st = (dbg_stY / superTileSize) * stilesX + (dbg_stX / superTileSize);
							int dl = ((debugLevel > -1) && (dbg_st == dbg_stile)) ? 3:0;
							if (dl > 0){
								System.out.println("compactSortShells():1, nTile = "+nTile+ ", nsTile = "+dbg_st);
							}
							int dls = ((debugLevel > -1) && (nTile == dbg_tile)) ? 3:0;
							if (tileData_src[nTile] != null){
								ArrayList<TileData> tdList = new ArrayList<TileData>();
								for (int nl = 0; nl < tileData_src[nTile].length; nl++){
									if (tileData_src[nTile][nl] != null){
										tdList.add(tileData_src[nTile][nl]);
									}
								}
								Collections.sort(tdList, new Comparator<TileData>() {
									@Override
									public int compare(TileData lhs, TileData rhs) {
										double lhs_d =lhs.getDisparity();
										double rhs_d =rhs.getDisparity();
										if (Double.isNaN(lhs_d) && Double.isNaN(rhs_d)) return 0;
										if (Double.isNaN(lhs_d)) return 1;
										if (Double.isNaN(rhs_d)) return -1;
										int sgn = (lhs.getDisparity() > rhs.getDisparity()) ? 1 : (lhs.getDisparity() < rhs.getDisparity() ) ? -1 : 0;
										return sgn;
									}
								});
								// increasing disparity
								for (int i = 0; i < tdList.size(); i++){
									tdList.get(i).setNewIndex(i);
								}
								if (tdList.size() > 0) {
									tile_data[nTile] = tdList.toArray(new TileData[0] );
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
						TileData [][]  tileData_src_dbg= tileData_src;
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							int dbg_stX = nTile % dbg_tilesX;   
							int dbg_stY = nTile / dbg_tilesX;
							int dbg_st = (dbg_stY / superTileSize) * stilesX + (dbg_stX / superTileSize);
							int dl = ((debugLevel > -1) && (dbg_st == dbg_stile)) ? 3:0;
							if (dl > 0){
								System.out.println("compactSortShells():2, nTile = "+nTile+ ", nsTile = "+dbg_st);
							}
							if (tile_data[nTile] != null){
								for (int i = 0; i < tile_data[nTile].length; i++){
									int [] neibs = tile_data[nTile][i].getNeighbors();
									for (int dir = 0; dir < neibs.length; dir++){
										if (neibs[dir] >= 0){
											int nTile1 = getNtileDir(nTile, dir);
											if (nTile1 >= 0) {
												if ((tile_data[nTile1] == null) || (tileData_src[nTile1][neibs[dir]] == null)){
													int dbg_sstile = tile_data[nTile][i].getDbgNsTile();
													int dbg_stileX = dbg_sstile % stilesX; 
													int dbg_stileY = dbg_sstile / stilesX;
													int dbg_tx = nTile % dbg_tilesX;
													int dbg_ty = nTile / dbg_tilesX;
													int dbg_dx = dbg_tx - (superTileSize * dbg_stileX + superTileSize/2); 
													int dbg_dy = dbg_ty - (superTileSize * dbg_stileY + superTileSize/2); 
													
													System.out.println("Null tile: "+nTile1+ " from "+nTile+", i="+i+", dir = "+dir+
															", dbg_stX="+dbg_stX+", dbg_stY="+dbg_stY+", dbg_st="+dbg_st+", neibs[dir]="+neibs[dir]+
															" dbg_nsTile = "+dbg_sstile +" ("+dbg_stileX+":"+dbg_stileY+")"+
															" nTile="+nTile+" ("+dbg_tx+":"+dbg_ty+")"+
															", deltas from src center: "+dbg_dx+":"+dbg_dy);
													neibs[dir] = -1;
												} else {
//													neibs[dir] = tile_data[nTile1][neibs[dir]].getNewIndex();
													neibs[dir] = tileData_src[nTile1][neibs[dir]].getNewIndex();
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
			return tile_data;
		}

	
		public void  checkShellsConnections (
				final TileData [][]     tileData,
				final int               debugLevel,
				final int               dbg_X,
				final int               dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
//			final TileData [][] tile_data = new TileData [nTiles][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int dbg_stile = (dbg_Y * superTileSize) * (stilesX * superTileSize) + (dbg_X * superTileSize);
//			final int dbg_tile = dbg_Y * stilesX + dbg_X;
			final int dbg_tilesX = stilesX * superTileSize;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							
							int dbg_stX = nTile % dbg_tilesX;   
							int dbg_stY = nTile / dbg_tilesX;
							int dbg_st = (dbg_stY / superTileSize) * stilesX + (dbg_stX / superTileSize);
							int dl = ((debugLevel > -1) && (dbg_st == dbg_stile)) ? 3:0;
							if (dl > 0){
								System.out.println("checkShellsConnections(), nTile = "+nTile+ ", nsTile = "+dbg_st);
							}
							if (tileData[nTile] != null){
								for (int nl = 0; nl < tileData[nTile].length; nl++){
									if (tileData[nTile][nl] != null){
										int  [] neibs = tileData[nTile][nl].getNeighbors();
										for (int dir = 0; dir < neibs.length; dir++){
											if (neibs[dir] >= 0){
												int nTile1 = getNtileDir(nTile, dir);
												if (nTile1 >= 0) {
													if ((tileData[nTile1] == null) || (tileData[nTile1][neibs[dir]] == null)){
														if (debugLevel > -1) {
															int dbg_sstile = tileData[nTile][nl].getDbgNsTile();
															int dbg_stileX = dbg_sstile % stilesX; 
															int dbg_stileY = dbg_sstile / stilesX;
															int dbg_tx = nTile % dbg_tilesX;
															int dbg_ty = nTile / dbg_tilesX;
															int dbg_dx = dbg_tx - (superTileSize * dbg_stileX + superTileSize/2); 
															int dbg_dy = dbg_ty - (superTileSize * dbg_stileY + superTileSize/2); 

															System.out.println("Broken link: "+nTile1+ " from "+nTile+", nl="+nl+", dir = "+dir+
																	", dbg_stX="+dbg_stX+", dbg_stY="+dbg_stY+", dbg_st="+dbg_st+", neibs[dir]="+neibs[dir]+
																	" dbg_nsTile = "+dbg_sstile +" ("+dbg_stileX+":"+dbg_stileY+")"+
																	" nTile="+nTile+" ("+dbg_tx+":"+dbg_ty+")"+
																	", deltas from src center: "+dbg_dx+":"+dbg_dy);
														}
														neibs[dir] = -2; // was broken link
													} else { // check if link is mutual
														int [] neibs_other = tileData[nTile1][neibs[dir]].getNeighbors();
														if (neibs_other[(dir + 4) % 8] != nl){
															if (debugLevel > -1) {
																int dbg_sstile = tileData[nTile][nl].getDbgNsTile();
																int dbg_stileX = dbg_sstile % stilesX; 
																int dbg_stileY = dbg_sstile / stilesX;
																int dbg_tx = nTile % dbg_tilesX;
																int dbg_ty = nTile / dbg_tilesX;
																int dbg_dx = dbg_tx - (superTileSize * dbg_stileX + superTileSize/2); 
																int dbg_dy = dbg_ty - (superTileSize * dbg_stileY + superTileSize/2); 

																System.out.println("Link not mutual: "+nTile1+ " from "+nTile+", nl="+nl+", dir = "+dir+
																		", dbg_stX="+dbg_stX+", dbg_stY="+dbg_stY+", dbg_st="+dbg_st+", neibs[dir]="+neibs[dir]+
																		" dbg_nsTile = "+dbg_sstile +" ("+dbg_stileX+":"+dbg_stileY+")"+
																		" nTile="+nTile+" ("+dbg_tx+":"+dbg_ty+")"+
																		", deltas from src center: "+dbg_dx+":"+dbg_dy+
																		", neibs_other["+((dir + 4) % 8)+"]="+neibs_other[(dir + 4) % 8]+
																		", dbg_nsTile other="+tileData[nTile1][neibs[dir]].getDbgNsTile());
															}
															neibs[dir] = -3; // not a mutual link (break only this side here)
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
		
		public void  addBackShellsConnections (
				final TileData [][]     tileData,
				final int               debugLevel,
				final int               dbg_X,
				final int               dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
//			final TileData [][] tile_data = new TileData [nTiles][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int dbg_stile = (dbg_Y * superTileSize) * (stilesX * superTileSize) + (dbg_X * superTileSize);
//			final int dbg_tile = dbg_Y * stilesX + dbg_X;
			final int dbg_tilesX = stilesX * superTileSize;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							
							int dbg_stX = nTile % dbg_tilesX;   
							int dbg_stY = nTile / dbg_tilesX;
							int dbg_st = (dbg_stY / superTileSize) * stilesX + (dbg_stX / superTileSize);
							int dl = ((debugLevel > -1) && (dbg_st == dbg_stile)) ? 3:0;
							if (dl > 0){
								System.out.println("checkShellsConnections(), nTile = "+nTile+ ", nsTile = "+dbg_st);
							}
							if (tileData[nTile] != null){
								for (int nl = 0; nl < tileData[nTile].length; nl++){
									if (tileData[nTile][nl] != null){
										int  [] neibs = tileData[nTile][nl].getNeighbors();
										for (int dir = 0; dir < neibs.length; dir++){
											if (neibs[dir] >= 0){
												int nTile1 = getNtileDir(nTile, dir);
												if (nTile1 >= 0) {
													if ((tileData[nTile1] == null) || (tileData[nTile1][neibs[dir]] == null)){
														if (debugLevel > 0) {
															int dbg_sstile = tileData[nTile][nl].getDbgNsTile();
															int dbg_stileX = dbg_sstile % stilesX; 
															int dbg_stileY = dbg_sstile / stilesX;
															int dbg_tx = nTile % dbg_tilesX;
															int dbg_ty = nTile / dbg_tilesX;
															int dbg_dx = dbg_tx - (superTileSize * dbg_stileX + superTileSize/2); 
															int dbg_dy = dbg_ty - (superTileSize * dbg_stileY + superTileSize/2); 

															System.out.println("Broken link: "+nTile1+ " from "+nTile+", nl="+nl+", dir = "+dir+
																	", dbg_stX="+dbg_stX+", dbg_stY="+dbg_stY+", dbg_st="+dbg_st+", neibs[dir]="+neibs[dir]+
																	" dbg_nsTile = "+dbg_sstile +" ("+dbg_stileX+":"+dbg_stileY+")"+
																	" nTile="+nTile+" ("+dbg_tx+":"+dbg_ty+")"+
																	", deltas from src center: "+dbg_dx+":"+dbg_dy);
														}
//														neibs[dir] = -2; // was broken link
													} else { // check if link is mutual
														int [] neibs_other = tileData[nTile1][neibs[dir]].getNeighbors();
														if (neibs_other[(dir + 4) % 8] != nl){
															if (debugLevel > 0) {
																int dbg_sstile = tileData[nTile][nl].getDbgNsTile();
																int dbg_stileX = dbg_sstile % stilesX; 
																int dbg_stileY = dbg_sstile / stilesX;
																int dbg_tx = nTile % dbg_tilesX;
																int dbg_ty = nTile / dbg_tilesX;
																int dbg_dx = dbg_tx - (superTileSize * dbg_stileX + superTileSize/2); 
																int dbg_dy = dbg_ty - (superTileSize * dbg_stileY + superTileSize/2); 

																System.out.println("Link not mutual: "+nTile1+ " from "+nTile+", nl="+nl+", dir = "+dir+
																		", dbg_stX="+dbg_stX+", dbg_stY="+dbg_stY+", dbg_st="+dbg_st+", neibs[dir]="+neibs[dir]+
																		" dbg_nsTile = "+dbg_sstile +" ("+dbg_stileX+":"+dbg_stileY+")"+
																		" nTile="+nTile+" ("+dbg_tx+":"+dbg_ty+")"+
																		", deltas from src center: "+dbg_dx+":"+dbg_dy+
																		", neibs_other["+((dir + 4) % 8)+"]="+neibs_other[(dir + 4) % 8]+
																		", dbg_nsTile other="+tileData[nTile1][neibs[dir]].getDbgNsTile());
															}
															if (neibs_other[(dir + 4) % 8] < 0 ){
																neibs_other[(dir + 4) % 8] = nl; // adding back link instead of missing one
															} else {
																int nTile2 = getNtileDir(nTile1, (dir + 4) % 8);
																if (    (nTile2 < 0) ||
																		(tileData[nTile2] == null) ||
																		(tileData[nTile2][neibs_other[(dir + 4) % 8]] == null)) {
																	neibs_other[(dir + 4) % 8] = nl; // adding back link instead of broken one
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
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		

		
		
		
		
		
		public int getTileLayersNumber (
				final TileData [][] tileData)
		{
			int num = 0;
			for (int i = 0; i < tileData.length; i++){
				if ((tileData[i] != null) && (tileData[i].length > num )){
					num  = tileData[i].length; 
				}
			}
			return num;
		}

		public int [] getTilesWH()
		{
			int [] wh = {stilesX*superTileSize, stilesY*superTileSize};
			return wh;
		}
		
		public double [][][] getTileDisparityStrengths (
				final TileData [][] tileData,
				final boolean useNaN)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int numLayers = getTileLayersNumber(tileData);
			final double [][][] disp_strength = new double [numLayers][2][tileData.length];
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (tileData[nTile] != null){
								for (int nl = 0; nl < tileData[nTile].length; nl++) if (tileData[nTile][nl]!=null){ // only beforfe compacting
									disp_strength[nl][0][nTile] = tileData[nTile][nl].getDisparity(useNaN);
									disp_strength[nl][1][nTile] = tileData[nTile][nl].getStrength();
								}
								if (useNaN){
									for (int nl = tileData[nTile].length; nl < numLayers; nl++){
										disp_strength[nl][0][nTile] = Double.NaN;
									}
								}
							} else if (useNaN){
								for (int nl = 0; nl < numLayers; nl++){
									disp_strength[nl][0][nTile] = Double.NaN;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return disp_strength;
		}
		public int [][][] getTileConnections (
				final TileData [][] tileData)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int numLayers = getTileLayersNumber(tileData);
			final int [][][] connections = new int [numLayers][tileData.length][8];
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (tileData[nTile] != null) {
								for (int nl = 0; nl < tileData[nTile].length; nl++) if (tileData[nTile][nl] != null) {
									for (int indx = 0 ; indx < tileData.length; indx++){
										for (int dir = 0; dir < 8; dir ++){
											if (tileData[nTile][nl].getNeighbor(dir) >= 0){
												connections[nl][nTile][dir] = tileData[nTile][nl].getNeighbor(dir)+1;
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
			return connections;
		}

		public int [][] getTileGenerator (
				final TileData [][] tileData)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int numLayers = getTileLayersNumber(tileData);
			final int [][] generators = new int [numLayers][tileData.length];
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (tileData[nTile] != null) {
								for (int nl = 0; nl < tileData[nTile].length; nl++) if (tileData[nTile][nl] != null) {
									for (int indx = 0 ; indx < tileData.length; indx++){
										generators[nl][nTile] = tileData[nTile][nl].getDbgNsTile();
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return generators;
		}
		public int [] getNumSurfaces (
				final TileData [][] tileData)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int [] surfaces = new int [tileData.length];
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (tileData[nTile] != null) {
								for (int nl = 0; nl < tileData[nTile].length; nl++) if (tileData[nTile][nl] != null) {
									surfaces[nTile] ++;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return surfaces;
		}
		
		
		public void showSurfaceDS (
				TileData [][] tileData,
				String title)
		{
			int [] wh = getTilesWH();
			double [][][] tds =  getTileDisparityStrengths (
					tileData,
					false); // useNaN);
			double [][][] tds_nan =  getTileDisparityStrengths (
					tileData,
					true); // useNaN);
			int [][] generators = getTileGenerator(tileData);
			int [] surfaces = getNumSurfaces (tileData);


			String [] titles = new String [5 * tds.length + 1];
			double [][] img_data = new double [titles.length][];
			for (int i = 0; i <tds.length; i++){
				titles[i + 0 * tds.length] = "disp_"+i;
				titles[i + 1 * tds.length] = "str_"+i;
				titles[i + 2 * tds.length] = "mdisp_"+i;
				titles[i + 3 * tds.length] = "mstr_"+i;
				titles[i + 4 * tds.length] = "gen_"+i;
				img_data[i + 0 * tds.length] = tds[i][0];
				img_data[i + 1 * tds.length] = tds[i][1];
				img_data[i + 2 * tds.length] = tds_nan[i][0];
				img_data[i + 3 * tds.length] = tds_nan[i][1];
				img_data[i + 4 * tds.length] = new double [generators[i].length];
				for (int j = 0; j < generators[i].length; j++){
					img_data[i + 4 * tds.length][j] = 0.01*generators[i][j];
				}
			}
			titles[5 * tds.length] = "layers";
			img_data[5 * tds.length] = new double [surfaces.length];
			for (int j = 0; j < surfaces.length; j++){
				img_data[5 * tds.length][j] =surfaces[j];
			}
			
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(img_data,  wh[0], wh[1], true, title, titles);
		}
		
		/**
		 * Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly
		 *                         scale ellipsoid (enlarge) 
		 * @param fraction_uni add fraction of the total weight to each tile
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @param debugLevel debug level
		 * @param dbg_X debug supertile X coordinate
		 * @param dbg_Y debug supertile Y coordinate
		 * @return per-tile (rounded up to contain whole supertiles) array of TileData instances
		 */
		public TileData [][] createTileShells (
				final boolean                   use_sel,
				final boolean                   divide_by_area,
				final double                    scale_projection,
				final double                    fraction_uni,
				final TilePlanes.PlaneData [][] planes,
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{

			double [][][][] fused_planes = fuseSupertilePlanes ( // some planes have zero weight that still participate - should they be eliminated?
					use_sel,           // final boolean                   use_sel,
					divide_by_area,    // final boolean                   divide_by_area,
					scale_projection,  // final double                    scale_projection,
					fraction_uni,      // final double                    fraction_uni,
					planes,            // final TilePlanes.PlaneData [][] planes,
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			
			int [][][][] surf_topology = getSupertilesTopology (
					planes,            // final TilePlanes.PlaneData [][] planes,
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			int [][][][][] overlapped_meshes = generateOverlappingMeshes (
					planes,            // final TilePlanes.PlaneData [][] planes,
					surf_topology ,    // final int [][][][]              corners, 
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			TileData [][] tileData = createTileShells (
					planes,            // final TilePlanes.PlaneData [][] planes,
					fused_planes,      // final double [][][][]           fusedSupertilePlanes,
					overlapped_meshes, // final int [][][][][]            lappingMeshes,				
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			System.out.println("addBackShellsConnections()");
			addBackShellsConnections (
					tileData,          // final TileData [][]     tileData_src,
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
//			
			System.out.println("checkShellsConnections()");
			checkShellsConnections (
					tileData,          // final TileData [][]     tileData_src,
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			tileData = compactSortShells (
					tileData,          // final TileData [][]     tileData_src,
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			
			showSurfaceDS (tileData, "tileData");			
			this.tileData = tileData;
			return tileData;
		}
		
		public int [] getTilesAssignStats(
				final int [][] tileLayers)
		{
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final int numThreads = threads.length;
			int [] stats = new int [STAT_LEN];
			final int [][] stats_all = new int [numThreads][STAT_LEN];
			final AtomicInteger ai_numThread = new AtomicInteger(0);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ml = 0; ml < tileLayers.length; ml++) if (tileLayers[ml] != null){
				final int fml = ml;
				ai_numThread.set(0);
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
							for (int nTile = ai.getAndIncrement(); nTile < tileLayers[fml].length; nTile = ai.getAndIncrement()) {
								if (tileLayers[fml][nTile] > 0){ // index + 1
									stats_all[numThread][STAT_ASSIGNED] ++; 		
								} else if (tileLayers[fml][nTile] == UNASSIGNED) {
									stats_all[numThread][STAT_UNASSIGNED] ++;		
								} else if (tileLayers[fml][nTile] == PROHOBITED) {
									stats_all[numThread][STAT_PROHIBITED] ++;		
								} else if (tileLayers[fml][nTile] == IMPOSSIBLE) {
									stats_all[numThread][STAT_IMPOSSIBLE] ++;		
								} else {
									System.out.println("Bug in getTilesAssignStats(): tileLayers["+fml+"]["+nTile+"]="+tileLayers[fml][nTile]);
									stats_all[numThread][0] ++; // prohibited		
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
				stats[STAT_NUM_ML]++; // number of non-null measurement layers
			}
			for (int nt = 0; nt < numThreads; nt ++){
				for (int i = 0 ; i < stats.length; i++ ){
					stats[i] += stats_all[nt][i];
				}
			}
			return stats;
		}
		
		public double getNormDispFromSurface(
				double disp_tile,
				double disp_surf,
				double disp_norm)
		{
			double disp_avg = 0.5 * (disp_tile + disp_surf);
			if (disp_avg <= disp_norm){
				return disp_tile - disp_surf;
			} else {
				return (disp_tile - disp_surf) * disp_norm / disp_avg;
			}
		}
		
		/**
		 * Convert from image tile index to the surface tile index (surface tiles are all
		 * full superTileSize),
		 * TODO: update/remove if surface grid will be trimmed to fit image
		 * Currently there are 324 tiles horizontally in the image and 328 in the surfaces
		 * @param nTile image tile index in scan order
		 * @return surface tile index in scan order
		 */
		public int getSurfaceTileIndex(
				int nTile)
		{
			// calculate index in tileData (has different dimensions - TODO: trim?
			int surfaceTilesX = stilesX * superTileSize;
			return surfaceTilesX * (nTile /  imageTilesX) + (nTile %  imageTilesX);
		}

		/**
		 * Convert from image tile index to the surface supertile index (surface tiles are all
		 * full superTileSize),
		 * TODO: update/remove if surface grid will be trimmed to fit image
		 * Currently there are 324 tiles horizontally in the image and 328 in the surfaces
		 * @param nTile image tile index in scan order
		 * @return surface tile index in scan order
		 */
		public int getSurfaceSuperTileIndex(
				int nTile)
		{
			// calculate index in tileData (has different dimensions - TODO: trim?
//			int surfaceTilesX = stilesX * superTileSize;
//			return surfaceTilesX * (nTile /  imageTilesX) + (nTile %  imageTilesX);
			return stilesX * (nTile /  imageTilesX / superTileSize) + (nTile %  imageTilesX)/superTileSize;
		}


		
		
		/**
		 * Convert from  surface tile index (surface tiles are all full superTileSize) to
		 * the image tile index
		 * TODO: update/remove if surface grid will be trimmed to fit image
		 * Currently there are 324 tiles horizontally in the image and 328 in the surfaces
		 * @param nSurfTile surface tile index in scan order
		 * @return image tile index in scan order or -1 if outside of the image tiles
		 */
		public int getImageTileIndex(
				int nSurfTile)
		{
			// calculate index in tileData (has different dimensions - TODO: trim?
			int surfaceTilesX = stilesX * superTileSize;
			int tx = nSurfTile %  surfaceTilesX;
			int ty = nSurfTile /  surfaceTilesX;
			if ((tx >= imageTilesX) || (ty >= imageTilesY)){
				return -1; // out of image
			}
			return imageTilesX * ty + tx;
		}
		
		/**
		 * Grow around surface tile, without going back. May still produce multi-layer result
		 * that will need to be filtered. Used to fill large flat gaps
		 * @param nsTile0
		 * @param nl0
		 * @return
		 */
		public boolean [][] growSimpleConnected(
				int nsTile0,
				int nl0)
		{
			int tilesX = stilesX * superTileSize;
			if ((tileData[nsTile0] == null) || (tileData[nsTile0][nl0] == null)){
				return null; // error
			}
			TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);
			boolean [][] simple_surf = new boolean [tileData.length][];
			int x0 = nsTile0 % tilesX;
			int y0 = nsTile0 / tilesX;
			ArrayList<Point> wave_list = new ArrayList<Point>();
			simple_surf[nsTile0] = new boolean [tileData[nsTile0].length];
			simple_surf[nsTile0][nl0] = true;
			wave_list.add(new Point(nsTile0, nl0));
			while (!wave_list.isEmpty()){
				Point p = wave_list.remove(0);
				int nsTile = p.x;
				int nl = p.y;
				int x = nsTile % tilesX;
				int y = nsTile / tilesX;
				int dirs = 0xff;
				if      (x > x0) dirs &= 0x1f;
				else if (x < x0) dirs &= 0xf1;
				if      (y > y0) dirs &= 0x7c;
				else if (y < y0) dirs &= 0xc7;
				int [] neibs = tileData[nsTile][nl].getNeighbors();
				for (int dir = 0; dir < neibs.length; dir ++) {
					int nl1 = neibs[dir];
					if ((nl1 >= 0) && ((dirs & (1 << dir)) != 0)){
						int nsTile1 = tnSurface.getNeibIndex(nsTile, dir);
						if ((nsTile1 >= 0) && ((simple_surf[nsTile1] == null) || !simple_surf[nsTile1][nl1])){
							if (simple_surf[nsTile1] == null) {
								simple_surf[nsTile1] = new boolean [tileData[nsTile1].length];
							}
							simple_surf[nsTile1][nl1] = true;
							wave_list.add(new Point(nsTile1, nl1));
						}
					}
				}
			}
			return simple_surf;
		}

		public int [][] growSimpleIConnected(
				int nsTile0,
				int nl0)
		{
			int tilesX = stilesX * superTileSize;
			if ((tileData[nsTile0] == null) || (tileData[nsTile0][nl0] == null)){
				return null; // error
			}
			TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);
			int [][] simple_surf = new int [tileData.length][];
			int x0 = nsTile0 % tilesX;
			int y0 = nsTile0 / tilesX;
			ArrayList<Point> wave_list = new ArrayList<Point>();
			simple_surf[nsTile0] = new int [tileData[nsTile0].length];
			simple_surf[nsTile0][nl0] = 1;
			wave_list.add(new Point(nsTile0, nl0));
			while (!wave_list.isEmpty()){
				Point p = wave_list.remove(0);
				int nsTile = p.x;
				int nl = p.y;
				int x = nsTile % tilesX;
				int y = nsTile / tilesX;
				int dirs = 0xff;
				if      (x > x0) dirs &= 0x1f;
				else if (x < x0) dirs &= 0xf1;
				if      (y > y0) dirs &= 0x7c;
				else if (y < y0) dirs &= 0xc7;
				int [] neibs = tileData[nsTile][nl].getNeighbors();
				for (int dir = 0; dir < neibs.length; dir ++) {
					int nl1 = neibs[dir];
					if ((nl1 >= 0) && ((dirs & (1 << dir)) != 0)){
						int nsTile1 = tnSurface.getNeibIndex(nsTile, dir);
						if ((nsTile1 >= 0) && ((simple_surf[nsTile1] == null) || (simple_surf[nsTile1][nl1] == 0))){
							if (simple_surf[nsTile1] == null) {
								simple_surf[nsTile1] = new int [tileData[nsTile1].length];
							}
							simple_surf[nsTile1][nl1] = simple_surf[nsTile][nl]+1;
							wave_list.add(new Point(nsTile1, nl1));
						}
					}
				}
			}
			return simple_surf;
		}
		
		
		
		public double [][] getShowSimpleConnected(
				int          nsTile0,
				boolean [][] simple_surf)
		{
			int numSurfaces = 1;
			for (int nsTile = 0; nsTile < simple_surf.length; nsTile++){
				if (simple_surf[nsTile] != null){
					int num_filled = 0;
					for (int i = 0; i < simple_surf[nsTile].length; i++) if (simple_surf[nsTile][i]) num_filled++;
					if (num_filled > numSurfaces) numSurfaces = num_filled; 
				}
			}
			double [][] img_data = new double [numSurfaces][tileData.length];
			for (int nsTile = 0; nsTile < tileData.length; nsTile++){
				if (simple_surf[nsTile] != null){
					int ns = 0;
					for (int nl = 0; nl < simple_surf[nsTile].length; nl++) if (simple_surf[nsTile][nl]){
						img_data[ns][nsTile] = tileData[nsTile][nl].getDisparity();
						ns++;
					}
					for (; ns < numSurfaces; ns++) img_data[ns][nsTile] = Double.NaN;
				} else {
					for (int ns = 0; ns < numSurfaces; ns++) img_data[ns][nsTile] = Double.NaN;
				}
			}
			// mark start point:
			for (int ns = 0; ns < numSurfaces; ns++) img_data[ns][nsTile0] = Double.NaN;
			return img_data;
		}
		public double [][] getShowSimpleConnected(
				int          nsTile0,
				int [][] simple_surf)
		{
			int numSurfaces = 1;
			for (int nsTile = 0; nsTile < simple_surf.length; nsTile++){
				if (simple_surf[nsTile] != null){
					int num_filled = 0;
					for (int i = 0; i < simple_surf[nsTile].length; i++) if (simple_surf[nsTile][i] !=0) num_filled++;
					if (num_filled > numSurfaces) numSurfaces = num_filled; 
				}
			}
			double [][] img_data = new double [numSurfaces][tileData.length];
			for (int nsTile = 0; nsTile < tileData.length; nsTile++){
				if (simple_surf[nsTile] != null){
					int ns = 0;
					for (int nl = 0; nl < simple_surf[nsTile].length; nl++) if (simple_surf[nsTile][nl] != 0){
						img_data[ns][nsTile] = tileData[nsTile][nl].getDisparity();
						ns++;
					}
					for (; ns < numSurfaces; ns++) img_data[ns][nsTile] = Double.NaN;
				} else {
					for (int ns = 0; ns < numSurfaces; ns++) img_data[ns][nsTile] = Double.NaN;
				}
			}
			// mark start point:
			for (int ns = 0; ns < numSurfaces; ns++) img_data[ns][nsTile0] = Double.NaN;
			return img_data;
		}

		public double [][] getShowSimpleConnectedDistanceLayer(
				boolean show_distance, // false - layer number
				int [][] simple_surf)
		{
			int numSurfaces = 1;
			for (int nsTile = 0; nsTile < simple_surf.length; nsTile++){
				if (simple_surf[nsTile] != null){
					int num_filled = 0;
					for (int i = 0; i < simple_surf[nsTile].length; i++) if (simple_surf[nsTile][i] !=0) num_filled++;
					if (num_filled > numSurfaces) numSurfaces = num_filled; 
				}
			}
			double [][] img_data = new double [numSurfaces][tileData.length];
			for (int nsTile = 0; nsTile < simple_surf.length; nsTile++){
				if (simple_surf[nsTile] != null){
					int ns = 0;
					for (int nl = 0; nl < simple_surf[nsTile].length; nl++) if (simple_surf[nsTile][nl] != 0){
						
						img_data[ns][nsTile] = show_distance? simple_surf[nsTile][nl]: nl;
						ns++;
					}
					for (; ns < numSurfaces; ns++) img_data[ns][nsTile] = Double.NaN;
				} else {
					for (int ns = 0; ns < numSurfaces; ns++) img_data[ns][nsTile] = Double.NaN;
				}
			}
			return img_data;
		}
		

		
		public void testSimpleConnected(
				int tileX,
				int tileY)
		{
			testSimpleConnected(stilesX * superTileSize * tileY+ tileX);
			printSurfaceConnections(
					tileX - 7, // int left,
					tileY - 7, // int top,
					15,        // int width,
					15,        // int height,
					1,         // int extraH,
					1);        // int extraV)
			
		}
		
		
		public void testSimpleConnected(
				int nsTile)
		{
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			double [][][] data =  new double [tileData[nsTile].length][][];
			double [][][] dist =  new double [tileData[nsTile].length][][];
			double [][][] layer = new double [tileData[nsTile].length][][];
			int [][][] simple_surf = new int [tileData[nsTile].length][][];
			int num_layers = 0;
			for (int nl = 0; nl < tileData[nsTile].length; nl++) if (tileData[nsTile][nl] != null) {
				simple_surf[nl] = growSimpleIConnected(nsTile, nl);
				data[nl] =  getShowSimpleConnected(nsTile, simple_surf[nl]);
				dist[nl] =  getShowSimpleConnectedDistanceLayer(true,simple_surf[nl]);
				layer[nl] = getShowSimpleConnectedDistanceLayer(false,simple_surf[nl]);
				num_layers += data[nl].length; 
			}
			double [][] img_data = new double [num_layers * 3][];
			String [] titles = new String [num_layers * 3];
			int indx = 0;
			for (int nl = 0; nl < data.length; nl++) if (data[nl] != null) {
				for (int i = 0; i < data[nl].length; i++){
					titles[indx + 0 * num_layers] = ""+nl+((data[nl].length>1)?(":"+i):"");
					titles[indx + 1 * num_layers] = "dist"+nl+((data[nl].length>1)?(":"+i):"");
					titles[indx + 2 * num_layers] = "layer"+nl+((data[nl].length>1)?(":"+i):"");
					img_data[indx + 0 * num_layers]= data[nl][i];
					img_data[indx + 1 * num_layers]= dist[nl][i];
					img_data[indx + 2 * num_layers]= layer[nl][i];
					indx++;
				}
			}
			sdfa_instance.showArrays(img_data,   stilesX * superTileSize,  stilesY * superTileSize, true, "simple_"+nsTile, titles);
		}
		
		public void printSurfaceConnections(
				int left,
				int top,
				int width,
				int height,
				int extraH,
				int extraV)
		{
			int tilesX = stilesX * superTileSize;
			int numSurfaces = 0;
	 			for (int tileY = top; tileY < (top + height); tileY++) {
				for (int tileX = left; tileX < (left + width); tileX++) {
					int nsTile = tileY * tilesX + tileX;
					if (tileData[nsTile] != null){
						if (tileData[nsTile].length > numSurfaces) numSurfaces = tileData[nsTile].length; 
					}
				}
			}
 			String hor_gap = "";
 			for (int i = 0; i < extraH; i++) hor_gap += " ";
 			String vert_gap = "";
 			for (int i = 0; i < extraV; i++) vert_gap += "\n";
 			for (int ns = 0; ns < numSurfaces; ns++){
 				System.out.println("\n ===== Surface number:"+ns+" =====");
 	 			for (int tileY = top; tileY < (top + height); tileY++) {
 	 				for (int l = 0; l < 3; l++){
 	 					for (int tileX = left; tileX < (left + width); tileX++) {
 	 						int nsTile = tileY * tilesX + tileX;
 	 						if ((tileData[nsTile]!= null) && (tileData[nsTile].length > ns) && (tileData[nsTile][ns] != null)){
 	 							int [] neibs = tileData[nsTile][ns].getNeighbors();
 	 							int master_st = tileData[nsTile][ns].getDbgNsTile();
 	 							switch (l){
 	 							case 0: System.out.print(" "+
 	 									((neibs[7] >= 0)?neibs[7]:" ")+
 	 									((neibs[0] >= 0)?neibs[0]:" ")+
 	 									((neibs[1] >= 0)?neibs[1]:" ")+" ");
 	 							break;
 	 							case 1: System.out.print(""+
 	 									((neibs[6] >= 0)?neibs[6]:" ")+
 	 									((master_st % stilesX) % 10)+
 	 									("*")+
 	 									((master_st / stilesX) % 10)+
 	 									((neibs[2] >= 0)?neibs[2]:" "));
 	 							break;
 	 							case 2: System.out.print(" "+
 	 									((neibs[5] >= 0)?neibs[5]:" ")+
 	 									((neibs[4] >= 0)?neibs[4]:" ")+
 	 									((neibs[3] >= 0)?neibs[3]:" ")+" ");
 	 							break;
 	 							}
 	 						} else {
 	 							switch (l){
 	 							case 0:
 	 							case 2: System.out.print("     "); break;
 	 							case 1: System.out.print("  .  "); break;
 	 							}
 	 						}
 	 						if (tileX < (left + width -1)){
 	 							System.out.print(hor_gap);
 	 						} else {
 	 							System.out.println();
 	 						}
 	 					} 	 					
 	 				}
 	 				if (tileY < (top + height -1)){
 	 					System.out.print(vert_gap);
 	 				}
 	 			}
 			}
		}
		
		/**
		 * Assign tiles to a certain disparity surface if there is only one surface candidate
		 * @param noEdge do not assign tiles to the surface edges (can not add border later)
		 * @param maxDiff maximal (normalized) disparity difference
		 * @param minDiffOther minimal disparity difference to closest 2-nd place candidate
		 * @param minStrength minimal processed (floor subtracted) correlation strength of the candidate 
		 * @param maxStrength maximal processed (floor subtracted) correlation strength of the candidate
		 * @param moveDirs +1 - allow moving tile closer to the camera (increase disparity, +2 - allow moving away
		 * @param dispNorm disparity normalization - disparity difference with average above it will be scaled down
		 * @param tileLayers measured tiles assignment (will be modified): -1 - prohibited, 0 - unassigned,
		 * >0 - number of surface where this tile is assigned plus 1.
		 * @param tileData per-tile, per layer array of TileData objects specifying surfaces to snap to
		 * @param dispStrength per measurement layer, combined disparity and strength array ([num_ml [2][])
		 * @param debugLevel debug level
		 * @param dbg_X debug tile X coordinate
		 * @param dbg_Y debug tile Y coordinate
		 * @return 
		 */
		public int [] assignTilesToSingleCandidate( // not used
				final boolean       noEdge,
				final double        maxDiff,
				final double        minDiffOther,
				final double        minStrength,
				final double        maxStrength,
				final int           moveDirs, // 1 increase disparity, 2 - decrease disparity, 3 - both directions
				final double        dispNorm, // disparity normalize (proportionally scale down disparity difference if above
//				final int [][]      tileLayers,
//				final TileData [][] tileData,
				final double [][][] dispStrength,
                final int           debugLevel,
				final int           dbg_X,
				final int           dbg_Y)
		{
			
			int [] stats_new = new int [NUM_STATS];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final int numThreads = threads.length;
			final int [][] stats_all = new int [numThreads][stats_new.length];
			final AtomicInteger ai_numThread = new AtomicInteger(0);
			final AtomicInteger ai = new AtomicInteger(0);
			final boolean en_lower =  (moveDirs & 1) != 0;
			final boolean en_higher = (moveDirs & 2) != 0;
			for (int ml = 0; ml < tileLayers.length; ml++) if (tileLayers[ml] != null){
				final int fml = ml;
				ai_numThread.set(0);
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
							for (int nTile = ai.getAndIncrement(); nTile < tileLayers[fml].length; nTile = ai.getAndIncrement()) {
								if (tileLayers[fml][nTile] == 0){ // unassigned only
									if (dispStrength[fml][1][nTile] < minStrength){
										stats_all[numThread][TOO_WEAK] ++;
									} else	if (dispStrength[fml][1][nTile] > maxStrength){
										stats_all[numThread][TOO_STRONG] ++;
									} else {
										// calculate index in tileData (has different dimensions - TODO: trim?
										int nSurfTile = getSurfaceTileIndex(nTile);
										if ((tileData[nSurfTile] == null) || (tileData[nSurfTile].length == 0)){
											stats_all[numThread][NO_SURF] ++;
											tileLayers[fml][nTile] = IMPOSSIBLE;
										} else {
//											double [] surf_disp_diff = new double [tileData[nSurfTile].length];
											int num_fit = 0;
											int num_fit_other = 0;
											int fit = -1;
											// Eliminate surface edges if prohibited
											boolean [] bad_surface = new boolean [tileData[nSurfTile].length]; // has less than 8 neighbors if controlled
											for (int ns = 0; ns < tileData[nSurfTile].length; ns++){
												double surf_disp_diff = getNormDispFromSurface (
														dispStrength[fml][0][nTile], // double disp_tile,
														tileData[nSurfTile][ns].getDisparity(), // double disp_surf,
														dispNorm); //double disp_norm)
												if (((surf_disp_diff >= 0) && en_higher) || ((surf_disp_diff <= 0) && en_lower)){
													// Eliminate surface edges if prohibited
													if (noEdge){
														int []neibs = tileData[nSurfTile][ns].getNeighbors();
														for (int i = 0; i < neibs.length; i++){
															if (neibs[i] < 0){
																bad_surface[ns] = true;
															}
														}
													}
													if (!bad_surface[ns]) {
														if (Math.abs(surf_disp_diff) <= maxDiff){
															fit = ns;   // no rating for fit "quality" here
															num_fit ++;
														}
														if (Math.abs(surf_disp_diff) <= minDiffOther){
															num_fit_other ++;
														}
													}
												}
											}
											if (num_fit < 1){
												stats_all[numThread][TOO_FAR] ++;
											} else if ((num_fit == 1) && (num_fit_other <= 1)){ // assign
												tileLayers[fml][nTile] = fit + 1;
												stats_all[numThread][NEW_ASSIGNED] ++;
											} else {
												stats_all[numThread][NOT_UNIQUE] ++;
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
			for (int nt = 0; nt < numThreads; nt ++){
				for (int i = 0 ; i < stats_new.length; i++ ){
					stats_new[i] += stats_all[nt][i];
				}
			}
			return stats_new;
		}
		
		public void printStats(int []stats)
		{
			boolean nothing = true;
			for (int i = 0; nothing && (i < stats.length); i++){
				nothing &= stats[i] == 0;
			}
			if (nothing) {
				System.out.println(" -- no changes --");
			} else {
				if (stats[NEW_ASSIGNED]     > 0) System.out.print(" NEW_ASSIGNED = " +     stats[NEW_ASSIGNED]);
				if (stats[NO_SURF]          > 0) System.out.print(" NO_SURF = " +          stats[NO_SURF]);
				if (stats[TOO_WEAK]         > 0) System.out.print(" TOO_WEAK = " +         stats[TOO_WEAK]);
				if (stats[TOO_STRONG]       > 0) System.out.print(" TOO_STRONG = " +       stats[TOO_STRONG]);
				if (stats[TOO_FAR]          > 0) System.out.print(" TOO_FAR = " +          stats[TOO_FAR]);
				if (stats[NOT_UNIQUE]       > 0) System.out.print(" NOT_UNIQUE = " +       stats[NOT_UNIQUE]);
				if (stats[REMOVED_TILES]    > 0) System.out.print(" REMOVED_TILES = " +    stats[REMOVED_TILES]);
				if (stats[REMOVED_CLUSTERS] > 0) System.out.print(" REMOVED_CLUSTERS = " + stats[REMOVED_CLUSTERS]);
			}
			System.out.println();
		}
		
		public boolean makesSensToTry(int [] stats)
		{
			return ((stats[NEW_ASSIGNED] > 0) && (stats[NOT_UNIQUE] > 0));
		}
		public int newAssigned(int [] stats)
		{
			return stats[NEW_ASSIGNED];
		}
		
		public void showAssignment(
				String title,
				final double [][][] dispStrength)
		{
			int layer_disp =     0;
			int layer_a_disp =   1;
			int layer_a_nan =    2;
			int layer_index =    3;
			int layer_strength = 4;
			int ng = 5;
			String [] titles = new String[ng * tileLayers.length];
			double [][] img_data = new double [titles.length][];
			for (int ml = 0; ml < tileLayers.length; ml ++){
				titles[ng * ml + layer_disp] =     "disp_"+ml;
				titles[ng * ml + layer_a_disp] =   "a_disp_"+ml;
				titles[ng * ml + layer_a_nan] =    "a_nan_"+ml;
				titles[ng * ml + layer_index] =    "index_"+ml;
				titles[ng * ml + layer_strength] = "strength_"+ml;
			}
			for (int ml = 0; ml < tileLayers.length; ml ++){
				if (dispStrength[ml] != null) {
					img_data[ng * ml + layer_disp] =     dispStrength[ml][0];
					img_data[ng * ml + layer_strength] = dispStrength[ml][1];
					img_data[ng * ml + layer_a_disp] =   new double [dispStrength[ml][0].length];
					img_data[ng * ml + layer_a_nan] =    new double [dispStrength[ml][0].length];
					img_data[ng * ml + layer_index] =    new double [dispStrength[ml][0].length];
					for (int nTile = 0;  nTile < dispStrength[ml][0].length; nTile++){
						int nSurfTile = getSurfaceTileIndex(nTile);
						if (tileLayers[ml][nTile] > 0){
							img_data[ng * ml + layer_a_disp][nTile] = tileData[nSurfTile][tileLayers[ml][nTile]-1].getDisparity();
							img_data[ng * ml + layer_a_nan][nTile] =  tileData[nSurfTile][tileLayers[ml][nTile]-1].getDisparity();
						} else {
							img_data[ng * ml + layer_a_disp][nTile] = dispStrength[ml][0][nTile];
							img_data[ng * ml + layer_a_nan][nTile] =  Double.NaN;
						}
						img_data[ng * ml + layer_index][nTile] = tileLayers[ml][nTile];
					}
				}
			}			
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(img_data,  imageTilesX, imageTilesY, true, title, titles);
		}
		
		/**
		 * Unassign tiles that have too few connected other tiles (or total weight of the cluster is too small)
		 * This is a single-threaded method
		 * @param minSize minimal tiles in the cluster
		 * @param minStrength minimal total strength of the cluster 
		 * @param dispStrength per measurement layer, combined disparity and strength array ([num_ml][2][])
		 * @param debugLevel debug level
		 * @param dbg_X debug tile X coordinate
		 * @param dbg_Y debug tile Y coordinate
		 * @return {number of tiles, number of clusters} removed  
		 */
		public int [] removeSmallClusters(
				final int           minSize,
				final double        minStrength,
				final double [][][] dispStrength,
                final int           debugLevel,
				final int           dbg_X,
				final int           dbg_Y)
		{
			boolean [][] wave_conf = new boolean [tileLayers.length][]; // true when wave or if confirmed
			int [] stats_new = new int [NUM_STATS];
			for (int ml = 0; ml < tileLayers.length; ml ++) if (tileLayers[ml] != null){
				wave_conf[ml] = new boolean[tileLayers[ml].length];
			}
			TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
//			TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);

			for (int ml = 0; ml < tileLayers.length; ml ++) if (tileLayers[ml] != null){
				for (int nTile0 = 0; nTile0 < tileLayers[ml].length; nTile0++) if ((tileLayers[ml][nTile0] > 0) && !wave_conf[ml][nTile0]){
					ArrayList<Point> wave_list = new ArrayList<Point>();
					double sum_weight = 0.0;
					int tailp = 0; // do not remove elements from the list while building the cluster, just advance tail pointer 
					Point p = new Point(nTile0, ml);
					sum_weight += dispStrength[p.y][1][p.x];
					wave_conf[p.y][p.x] = true;
					wave_list.add(p);
					while (tailp < wave_list.size()){
						Point pt = wave_list.get(tailp++);
						int nSurfTile1 = getSurfaceTileIndex(pt.x);
						int nl1 = tileLayers[pt.y][pt.x] - 1; // zero-based from 1-based
						int [] neibs = tileData[nSurfTile1][nl1].getNeighbors();
						for (int dir  = 0; dir < tnImage.dirs; dir++) if (neibs[dir] >= 0){
							int nTile1 = tnImage.getNeibIndex(pt.x, dir);
							if (nTile1 >= 0) {
								for (int ml1 = 0; ml1 < tileLayers.length; ml1++) {
									// may be several ml tiles on the same surface - count them all
									if ((tileLayers[ml1] != null) && (tileLayers[ml1][nTile1] == (neibs[dir] +1)) && !wave_conf[ml1][nTile1]){
										Point p1 = new Point(nTile1, ml1);
										sum_weight += dispStrength[p1.y][1][p1.x];
										wave_conf[p1.y][p1.x] = true;
										wave_list.add(p1);
									}
								}
							}
						}
					}
					// See if it is a good cluster
					if ((wave_list.size() < minSize) || (sum_weight < minStrength)){
						while (wave_list.size() > 0){
							Point pt = wave_list.remove(0);
							tileLayers[pt.y][pt.x] = 0;
							wave_conf [pt.y][pt.x] = false; // not necessary
							
							stats_new[REMOVED_TILES]++;
						}
						stats_new[REMOVED_CLUSTERS]++;
					} else { // it is a strong cluster, nothing to do here (it is already marked in wave_conf[][]
						
					}
				}
			}			
			return stats_new;
		}

		/**
		 * Assign (weak) tile surrounded by assigned one to the disparity of the farthest tile (lowest disparity).
		 * This is a single-threaded method
		 * @param noEdge do not assign tiles to the surface edges (can not add border later)
		 * @param minNeib minimal number of occupied directions (of 8), several occupied levels count as one
		 * @param maxStrength maximal strength of the tile to assign (strong one may make trust its disparity after all)
		 * @param includeImpossible count impossible (blocked, on the image edge,...) tiles as if assigned towards
		 * the number of occupied directions
		 * @param dispStrength per measurement layer, combined disparity and strength array ([num_ml][2][])
		 * @param debugLevel debug level
		 * @param dbg_X debug tile X coordinate
		 * @param dbg_Y debug tile Y coordinate
		 * @return {number of tiles, number of clusters} removed  
		 */
		
		public int [] assignFromFarthest(
				final boolean       noEdge,
				final int           minNeib,
				final double        maxStrength,
				final boolean       includeImpossible, // count prohibited neighbors as assigned
				final double [][][] dispStrength,
                final int           debugLevel,
				final int           dbg_X,
				final int           dbg_Y)
		{
			
			final int [][] tileLayers_src = tileLayers.clone();
			for (int i = 0; i < tileLayers_src.length; i++){
				if (tileLayers_src[i] != null){
					tileLayers_src[i] = tileLayers[i].clone();
				}
			}
			int [] stats_new = new int [NUM_STATS];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final int numThreads = threads.length;
			final int [][] stats_all = new int [numThreads][stats_new.length];
			final AtomicInteger ai_numThread = new AtomicInteger(0);
			final AtomicInteger ai = new AtomicInteger(0);
			final TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
			final TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);
			
			for (int ml = 0; ml < tileLayers.length; ml++) if (tileLayers[ml] != null){
				final int fml = ml;
				ai_numThread.set(0);
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
							for (int nTile = ai.getAndIncrement(); nTile < tileLayers_src[fml].length; nTile = ai.getAndIncrement()) {
								//nTile is in image, not surface coordinates 
								int dbg_tileX = nTile % imageTilesX;
								int dbg_tileY = nTile / imageTilesX;
								int dl = ((debugLevel > -1) && (dbg_tileX == dbg_X ) && (dbg_tileY == dbg_Y ))?3:0;
								if (dl > 0){
									System.out.println("assignFromFarthest, nTile = " + nTile);
								}
								if (tileLayers_src[fml][nTile] == 0){ // unassigned only
									if (dispStrength[fml][1][nTile] > maxStrength){
										stats_all[numThread][TOO_STRONG] ++;
									} else {
										// find number of tiles around (x,y) that have surface connection to this one
										// (multiple ml count as one), and which one has the lowest disparity
										int nSurfTile = getSurfaceTileIndex(nTile);
										double min_disp = Double.NaN;
										int best_nSurf = -1;
										int numNeibs = 0;
										for (int dir  = 0; dir < tnImage.dirs; dir++) {
											int nTile1 = tnImage.getNeibIndex(nTile,dir);
											boolean neib_exists = false;
											if (nTile1 >= 0){
												for (int ml_other = 0; ml_other < tileLayers_src.length; ml_other++) if (tileLayers_src[ml_other] != null){
													if (tileLayers_src[ml_other][nTile1] < 0 ) { //
														neib_exists |= includeImpossible;
													} else if (tileLayers_src[ml_other][nTile1] > 0 ){
														int nSurfTile1 = tnSurface.getNeibIndex(nSurfTile,dir);
														int nSurf = tileData[nSurfTile1][tileLayers_src[ml_other][nTile1] - 1].getNeighbor(tnSurface.opposite(dir));
														if (nSurf >= 0){
															neib_exists = true;
															double disp = tileData[nSurfTile][nSurf].getDisparity();
															boolean bad_surface = false;
															if (noEdge) {
																int [] neibs =  tileData[nSurfTile][nSurf].getNeighbors();
																for (int i = 0; i < neibs.length; i++) if (neibs[i] <0) {
																	bad_surface = true;
																	break;
																}
															}
															
															if (!(disp >= min_disp) && !bad_surface) {
																best_nSurf = nSurf;
																min_disp = disp;
															}
														}
													}
												}
											} else {
												neib_exists = includeImpossible; 
											}
											if (neib_exists){
												numNeibs++;
											}
										}
										if ((numNeibs >= minNeib) && (best_nSurf >= 0)){
											tileLayers[fml][nTile] = best_nSurf + 1;
											stats_all[numThread][NEW_ASSIGNED] ++;
										}
									}
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}
			for (int nt = 0; nt < numThreads; nt ++){
				for (int i = 0 ; i < stats_new.length; i++ ){
					stats_new[i] += stats_all[nt][i];
				}
			}
			return stats_new;
		}
		
		
		/**
		 * Assign tiles to a certain disparity surface if there is only one surface candidate
		 * @param noEdge do not assign tiles to the surface edges (can not add border later)
		 * @param useCenter only assign outside of 8x8 center if no suitable alternative
		 * @param maxDiff maximal (normalized) disparity difference
		 * @param minDiffOther minimal disparity difference to closest 2-nd place candidate
		 * @param minStrength minimal processed (floor subtracted) correlation strength of the candidate 
		 * @param maxStrength maximal processed (floor subtracted) correlation strength of the candidate
		 * @param minSurfStrength minimal surface strength at the tile location
		 * @param moveDirs +1 - allow moving tile closer to the camera (increase disparity, +2 - allow moving away
		 * @param enMulti allow assignment when several surfaces fit
		 * @param surfStrPow raise surface strengths ratio to this power when comparing candidates
		 * @param addStrength  add to strengths when calculating pull of assigned tiles
		 * @param sigma radius of influence (in tiles) of the previously assigned tiles
		 * @param nSigma maximal relative to radius distance to calculate influence
		 * @param minPull additional pull for no-tile surfaces (to avoid division by zero)
		 * @param minAdvantage minimal ratio of the best surface candidate to the next one to make selection
		 * @param dispNorm disparity normalization - disparity difference with average above it will be scaled down
		 * @param tileLayers measured tiles assignment (will be modified): -1 - prohibited, 0 - unassigned,
		 * >0 - number of surface where this tile is assigned plus 1.
		 * @param tileData per-tile, per layer array of TileData objects specifying surfaces to snap to
		 * @param dispStrength per measurement layer, combined disparity and strength array ([num_ml][2][])
		 * @param debugLevel debug level
		 * @param dbg_X debug tile X coordinate
		 * @param dbg_Y debug tile Y coordinate
		 * @return statistics array
		 */
		
		public int [] assignTilesToSurfaces(
				final boolean       noEdge,
				final boolean       useCenter,
				final double        maxDiff,
				final double        minDiffOther, // should be >= maxDiff
				final double        minStrength,
				final double        maxStrength,
				final double        minSurfStrength, // minimal surface strength at the tile location
				final int           moveDirs, // 1 increase disparity, 2 - decrease disparity, 3 - both directions
				final boolean       enMulti,
				final double        surfStrPow, // surface strength power
				final double        addStrength, //
				final double        sigma,
				final double        nSigma,
				final double        minPull,
				final double        minAdvantage,
				final double        dispNorm, // disparity normalize (proportionally scale down disparity difference if above
				final double [][][] dispStrength,
                final int           debugLevel,
				final int           dbg_X,
				final int           dbg_Y)
		{
			
			final int [][] tileLayers_src = tileLayers.clone();
			for (int i = 0; i < tileLayers_src.length; i++){
				if (tileLayers_src[i] != null){
					tileLayers_src[i] = tileLayers[i].clone();
				}
			}
			int [] stats_new = new int [NUM_STATS];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final int numThreads = threads.length;
			final int [][] stats_all = new int [numThreads][stats_new.length];
			final AtomicInteger ai_numThread = new AtomicInteger(0);
			final AtomicInteger ai = new AtomicInteger(0);
			final boolean en_lower =  (moveDirs & 1) != 0;
			final boolean en_higher = (moveDirs & 2) != 0;
			final double radius = sigma * nSigma;
			final double rsigma2 = 1.0 / ( 2.0 * sigma * sigma);
			final int iradius = (int) Math.round(radius + 0.001);  
			final int field_size = 2 * iradius + 1;
			final int center_index = iradius * (field_size + 1);
			final double cost_start = 1.0;
			final double cost_ortho = 1.0;
			final double cost_diag  = 1.5; // Math.sqrt(2.0);
			final int surfTilesX = stilesX * superTileSize;
			final int [] ldirs8 = {
					-field_size,
					-field_size + 1,
					1,
					field_size + 1,
					field_size,
					field_size - 1,
					-1,
					-field_size - 1}; 
			
			for (int ml = 0; ml < tileLayers.length; ml++) if (tileLayers[ml] != null){
				final int fml = ml;
				ai_numThread.set(0);
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
							for (int nTile = ai.getAndIncrement(); nTile < tileLayers_src[fml].length; nTile = ai.getAndIncrement()) {
								//nTile is in image, not surface coordinates 
								int dbg_tileX = nTile % imageTilesX;
								int dbg_tileY = nTile / imageTilesX;
								int dl = ((debugLevel > -1) && (dbg_tileX == dbg_X ) && (dbg_tileY == dbg_Y ))?3:0;
								
								if (tileLayers_src[fml][nTile] == 0){ // unassigned only
									if (dispStrength[fml][1][nTile] < minStrength){
										stats_all[numThread][TOO_WEAK] ++;
									} else	if (dispStrength[fml][1][nTile] > maxStrength){
										stats_all[numThread][TOO_STRONG] ++;
									} else {
										// calculate index in tileData (has different dimensions - TODO: trim?
										int nSurfTile = getSurfaceTileIndex(nTile);
										int nSurfSuperTile = getSurfaceSuperTileIndex(nTile);
										if ((tileData[nSurfTile] == null) || (tileData[nSurfTile].length == 0)){
											stats_all[numThread][NO_SURF] ++;
											tileLayers[fml][nTile] = IMPOSSIBLE;
										} else {
//											double [] surf_disp_diff = new double [tileData[nSurfTile].length];
											int num_fit = 0;
											int num_fit_other = 0;
											int fit = -1;
											int num_fit_center = 0;
											int num_fit_other_center = 0;
											int fit_center = -1;
											boolean [] bad_surface = new boolean [tileData[nSurfTile].length]; // has less than 8 neighbors if controlled
											boolean [] center_surface = new boolean [tileData[nSurfTile].length]; // is inside the center part
											for (int ns = 0; ns < tileData[nSurfTile].length; ns++){
												double surf_disp_diff = getNormDispFromSurface (
														dispStrength[fml][0][nTile], // double disp_tile,
														tileData[nSurfTile][ns].getDisparity(), // double disp_surf,
														dispNorm); //double disp_norm)
												if (((surf_disp_diff >= 0) && en_higher) || ((surf_disp_diff <= 0) && en_lower)){
													// Eliminate surface edges if prohibited
													if (noEdge){
														int []neibs = tileData[nSurfTile][ns].getNeighbors();
														for (int i = 0; i < neibs.length; i++) if (neibs[i] < 0) {
															bad_surface[ns] = true;
															break;
														}
													}
													if (!bad_surface[ns]) {
														if (Math.abs(surf_disp_diff) <= maxDiff){
															fit = ns;   // no rating for fit "quality" here
															num_fit ++;
														}
														if (Math.abs(surf_disp_diff) <= minDiffOther){
															num_fit_other ++;
														}
														// separately calculate for center of the tiles getDbgNsTile
														TileData td = tileData[nSurfTile][ns];
														if (tileData[nSurfTile][ns].getDbgNsTile() == nSurfSuperTile){
															center_surface[ns] = true;
															if (Math.abs(surf_disp_diff) <= maxDiff){
																fit_center = ns;   // no rating for fit "quality" here
																num_fit_center ++;
															}
															if (Math.abs(surf_disp_diff) <= minDiffOther){
																num_fit_other_center ++;
															}
														}
													}
												}
											}
											if (useCenter && (num_fit_center > 0)) { // use only center if possible
												fit =           fit_center;
												num_fit =       num_fit_center;
												num_fit_other = num_fit_other_center;
												for (int i = 0; i < center_surface.length; i++){
													if (!center_surface[i]) bad_surface[i] = true;
												}
											}
											if (num_fit < 1){
												stats_all[numThread][TOO_FAR] ++;
											} else if ((num_fit == 1) && (num_fit_other <= 1)){ // assign
												tileLayers[fml][nTile] = fit + 1;
												stats_all[numThread][NEW_ASSIGNED] ++;
											} else if (!enMulti) {
												stats_all[numThread][NOT_UNIQUE] ++;
											} else { // multi, enabled
												int [] candidates =       new int     [num_fit_other];
												boolean [] close_enough = new boolean [num_fit_other];
												num_fit_other = 0;
												for (int ns = 0; ns < tileData[nSurfTile].length; ns++) if (!bad_surface[ns]){ // not all 8 neighbors
													double surf_disp_diff = getNormDispFromSurface (
															dispStrength[fml][0][nTile], // double disp_tile,
															tileData[nSurfTile][ns].getDisparity(), // double disp_surf,
															dispNorm); //double disp_norm)
													if (((surf_disp_diff >= 0) && en_higher) || ((surf_disp_diff <= 0) && en_lower)){
														if (Math.abs(surf_disp_diff) <= minDiffOther){
															close_enough[num_fit_other] = (Math.abs(surf_disp_diff) <= maxDiff);
															candidates[num_fit_other++] = ns;
														}
													}
												}
												if (dl > 0) {
													System.out.print("assignTilesToSurfaces(): nTile="+nTile+", candidates=");
													for (int ii = 0; ii < candidates.length; ii++){
														System.out.print(" "+candidates[ii]);	
													}
													System.out.println();	
												}
												double [][][] distances = new double [num_fit_other][field_size  * field_size ][];
												// for each local index get surface tile index
												int [] surfIndices =  new int [field_size  * field_size];
												int [] imageIndices = new int [field_size  * field_size];
												int stx0 = (nTile % imageTilesX) - iradius; // 
												int sty0 = (nTile / imageTilesX) - iradius;
												for (int iy = 0; iy < field_size; iy++){
													for (int ix = 0; ix < field_size; ix++){
														int indx = iy * field_size + ix;
														surfIndices[indx] = (sty0 + iy)*  surfTilesX + (stx0 + ix);
														imageIndices[indx] = getImageTileIndex(surfIndices[indx]);
													}
												}
												// calculate distances from the center point for each surface with wave
												// algorithm, first regardless of who is closer.
												// later, when comparing pairs only use the same side
												for (int isurf = 0; isurf < num_fit_other; isurf++){
													ArrayList<Point> lwave = new ArrayList<Point>();
													Point p0 = new Point(center_index, candidates[isurf]);
													distances[isurf][p0.x] = new double [tileData[nSurfTile].length];
													distances[isurf][p0.x][p0.y] = cost_start;
													if (dl > 0) {
														System.out.println("Add: p0.x="+p0.x+", p0.y="+p0.y);
													}
													lwave.add(p0);
													// run wave build radius (plus 1.0) along each surface connections,
													// until next radius >= radius
													while (!lwave.isEmpty()){
														p0 = lwave.remove(0);
														TileData [] dbg_tileData = tileData[surfIndices[p0.x]];
														int [] neibs = tileData[surfIndices[p0.x]][p0.y].getNeighbors();
														if (dl > 0) {
															System.out.println("Remove: p0.x="+p0.x+", p0.y="+p0.y+" surfIndices[p0.x]="+surfIndices[p0.x]+
																	" neibs:"+
																	" [ "+((neibs[0] >= 0)? neibs[0]:"-")+
																	" | "+((neibs[1] >= 0)? neibs[1]:"-")+
																	" | "+((neibs[2] >= 0)? neibs[2]:"-")+
																	" | "+((neibs[3] >= 0)? neibs[3]:"-")+
																	" | "+((neibs[4] >= 0)? neibs[4]:"-")+
																	" | "+((neibs[5] >= 0)? neibs[5]:"-")+
																	" | "+((neibs[6] >= 0)? neibs[6]:"-")+
																	" | "+((neibs[7] >= 0)? neibs[7]:"-")+
																	" ]");
														}
														// try ortho directions first
														double new_dist =  distances[isurf][p0.x][p0.y] + cost_ortho;
														if (new_dist <= (radius + cost_start)) {
															for (int dir = 0; dir < 8; dir +=2) if ( neibs[dir] >= 0){
																Point pn = new Point (p0.x + ldirs8[dir], neibs[dir]);
																if (distances[isurf][pn.x] == null){
																	distances[isurf][pn.x] = new double [tileData[surfIndices[pn.x]].length];
																}
																if ((distances[isurf][pn.x][pn.y] == 0) || (distances[isurf][pn.x][pn.y] > new_dist)){
																	if (dl > 0) {
																		System.out.println("Add ortho: p0.x="+p0.x+", p0.y="+p0.y+
																				" distances["+isurf+"]["+pn.x+"]["+pn.y+"]="+distances[isurf][pn.x][pn.y]+
																				", new_dist="+new_dist);
																	}
																	distances[isurf][pn.x][pn.y] = new_dist;
																	lwave.add(pn);
																}
															}
														}
														// try diagonal directions second
														new_dist =  distances[isurf][p0.x][p0.y] + cost_diag;
														if (new_dist <= (radius + cost_start)) {
															for (int dir = 1; dir < 8; dir +=2) if ( neibs[dir] >= 0){
																Point pn = new Point (p0.x + ldirs8[dir], neibs[dir]);
																if (distances[isurf][pn.x] == null){
																	distances[isurf][pn.x] = new double [tileData[surfIndices[pn.x]].length];
																}
																if ((distances[isurf][pn.x][pn.y] == 0) || (distances[isurf][pn.x][pn.y] > new_dist)){
																	if (dl > 0) {
																		System.out.println("Add diag: p0.x="+p0.x+", p0.y="+p0.y+
																				" distances["+isurf+"]["+pn.x+"]["+pn.y+"]="+distances[isurf][pn.x][pn.y]+
																				", new_dist="+new_dist);
																	}
																	distances[isurf][pn.x][pn.y] = new_dist;
																	lwave.add(pn);
																}
															}
														}
													}
												}											
												if (dl > 0) {
													for (int cand = 0; cand < distances.length; cand ++){
														int num_dist_layers = 0;
														for (int i = 0; i < distances[cand].length; i++){
															if ((distances[cand][i] != null) && (distances[cand][i].length > num_dist_layers)){
																num_dist_layers = distances[cand][i].length;
															}
														}
														for (int dist_l = 0; dist_l < num_dist_layers; dist_l++){
															System.out.println("Candidate #"+cand+", layer "+dist_l);
															for (int ddy = 0; ddy < field_size; ddy ++){
																for (int ddx = 0; ddx < field_size; ddx ++){
																	if ((distances[cand][ddy * field_size + ddx] == null) ||
																			(distances[cand][ddy * field_size + ddx].length <= dist_l) ||
																			(distances[cand][ddy * field_size + ddx][dist_l] == 0)){
																		System.out.print("--- ");
																	} else {
																		System.out.print(distances[cand][ddy * field_size + ddx][dist_l]+" ");
																	}
																}
																System.out.println();
															}
														}
													}
												}
												
												
												// pulls belong to pairs, not individual surfaces (difference when they cross)
												double [][] surface_pulls = new double [num_fit_other][num_fit_other];
												// now calculate advantage of each one surface (close_enough) to each other (as a ratio)
												// and then see if it is above minAdvantage
												for (int other_ml = 0; other_ml < tileLayers.length; other_ml++) if (tileLayers[other_ml] != null){
													for (int lindx = 0; lindx < imageIndices.length; lindx ++){
														if (imageIndices[lindx] >= 0) {
															int nsurf = tileLayers_src[other_ml][imageIndices[lindx]] - 1; // assigned surface number (>=0)
															if ( nsurf >= 0){
																double strength = dispStrength[other_ml][1][imageIndices[lindx]] + addStrength; // add strength so very weak count
																// see if this tile belongs to any of the considered surfaces
																int num_found = 0;
																boolean [] on_surface = new boolean [num_fit_other];
																for (int i = 0; i < num_fit_other; i++) {
																	if ((distances[i][lindx] != null) && (distances[i][lindx][nsurf] > 0.0)){
																		num_found++;
																		on_surface[i] = true;
																	}
																}
																if (num_found > 0) { // tile lies on at least one of the considered surfaces
																	for (int is1 = 0; is1 < num_fit_other; is1++) if (on_surface[is1]) {
																		// is2 can be any other candidate, just check it is on the same side
																		// of is1 as in the center (new tile)
																		for (int is2 = 0; is2 < num_fit_other; is2++) if (is2 != is1) {
																			boolean good_pair = true;
																			if (distances[is2][lindx] != null) { // otherwise OK
																				for (int i = 0; i < distances[is2][lindx].length; i++){
																					if (distances[is2][lindx][i] >= 0){
																						if (	((is2 > is1) && (i < nsurf)) ||
																								((is2 < is1) && (i > nsurf))) {
																							good_pair = false; // surfaces cross between
																							break;
																						}
																					}
																				}
																			}
																			if (good_pair){
																				double r = distances[is1][lindx][nsurf] - cost_start;
																				// pull to is1 when in pair with is2
																				surface_pulls[is1][is2] += Math.exp(- r * r * rsigma2) * strength ; 
																			}
																		}
																	}
																}
															}
														}
													}
												}
												
												double [][] advantages = new double [num_fit_other][num_fit_other];
												for (int is1 = 0; is1 < num_fit_other; is1++){
													for (int is2 = is1 + 1; is2 < num_fit_other; is2++){
														double ad1 = surface_pulls[is1][is2] + minPull;
														double ad2 = surface_pulls[is2][is1] + minPull;
														// normally minPull >0.0, if not - prevent div by zero
														if ((ad1 == 0) || (ad2 == 0)){
															if ((ad1 == 0) && (ad2 == 0)){
																ad1 = 1.0;
																ad2 = 1.0;
															} else if (ad1 == 0) {
																ad2 = 2.0 * minAdvantage;
																ad1 = 1.0;
															} else {
																ad1 = 2.0 * minAdvantage;
																ad2 = 1.0;
															}
														}
														advantages[is1][is2] = ad1/ad2;
														advantages[is2][is1] = ad2/ad1;
														if (surfStrPow != 0.0){ // consider surface strength also 
															double str1 = tileData[nSurfTile][candidates[is1]].getStrength();
															double str2 = tileData[nSurfTile][candidates[is1]].getStrength();
															if ((str1 > 0.0) && (str2 > 0.0)){
																advantages[is1][is2] *= Math.pow(str1/str2, surfStrPow);
																advantages[is2][is1] = 1.0/advantages[is1][is2];
															} else if (str1 > 0.0) {
																advantages[is1][is2] = 2.0 * minAdvantage; // sure will win
																advantages[is2][is1] = (minAdvantage > 0.0) ? (1.0/advantages[is1][is2]) : 0.0; 
																//minAdvantage
															} else if (str2 > 0.0) {
																advantages[is2][is1] = 2.0 * minAdvantage; // sure will win
																advantages[is1][is2] = (minAdvantage > 0.0) ? (1.0/advantages[is2][is1]) : 0.0; 
															} else { // both zero - do nothing about surface strengths
																
															}
														}
													}
												}
												// Now see if we have a winner that is (2 could not satisfy, look for the first:
												// a) close enough, and
												// b) sufficient advantage over all other candidates
												fit = -1;
												for (int is1 = 0; is1 < num_fit_other; is1++){
													if (close_enough[is1]){ //
														boolean is_a_winner = true; 
														for (int is2 = is1 + 1; is2 < num_fit_other; is2++){
															if (advantages[is1][is2] < minAdvantage){
																if (dl > 0) {
																	System.out.println("assignTilesToSurfaces() advantages["+is1+"]["+is2+"]="+advantages[is1][is2]);
																}
																is_a_winner = false;
																if (dl > 0) {
																	System.out.println("assignTilesToSurfaces(): Not a winner, advantages < "+minAdvantage);
																}
																break;
															}
														}														
														if (is_a_winner){
															fit = is1;
															if (dl > 0) {
																System.out.println("assignTilesToSurfaces(): "+is1+" is a winner!");
															}
															break;
														}
													}
												}
												if (fit >= 0) {
													tileLayers[fml][nTile] = candidates[fit] + 1;
													stats_all[numThread][NEW_ASSIGNED] ++;
												} else {
													stats_all[numThread][NOT_UNIQUE] ++;
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
			for (int nt = 0; nt < numThreads; nt ++){
				for (int i = 0 ; i < stats_new.length; i++ ){
					stats_new[i] += stats_all[nt][i];
				}
			}
			return stats_new;
		}

		
		public int [] InitTilesAssignment(
				final boolean            force,
				final double [][][]      dispStrength,
				final boolean [][]       tileSel,
                final int                debugLevel)
		{
			if (force || (this.tileLayers == null)) {
				int [][] tileLayers = new int [tileSel.length][];
				for (int ml = 0; ml < tileSel.length; ml++){
					if (tileSel[ml] != null){
						tileLayers[ml] = new int [tileSel[ml].length];
						for (int i = 0; i < tileSel[ml].length; i++){
							tileLayers[ml][i] = tileSel[ml][i] ? 0: -1; // 0 - unassigned, -1 - prohibited	
						}
					}
				}
				this.tileLayers = tileLayers;
			}
			int []stats = getTilesAssignStats(tileLayers);
			if (debugLevel >= -1) {
				System.out.println("sortTilesToSurfaces(): using "+stats[STAT_NUM_ML] +" measurement layers"+
						", number of assigned tiles: "+stats[STAT_ASSIGNED]+
						", number of unassigned tiles: "+stats[STAT_UNASSIGNED]+
						", number of prohibited tiles: "+stats[STAT_PROHIBITED]+
						", number of impossible tiles: "+stats[STAT_IMPOSSIBLE]);
			}
			return stats;
		}

		
		public boolean [][] extractSelection(
				final int         debugLevel,
				final int         dbg_X,
				final int         dbg_Y)
		{
			//			boolean [][] wave_conf = new boolean [tileLayers.length][]; // true when wave or if confirmed
			int num_tiles = 0;
			for (int ml = 0; ml < tileLayers.length; ml++) if ((tileLayers[ml] != null) && (tileLayers[ml].length > 0)){
				num_tiles = tileLayers[ml].length;
				break;
			}
			boolean [][] assigned_tiles = new boolean [num_tiles][];
			for (int ml = 0; ml < tileLayers.length; ml ++) if (tileLayers[ml] != null){
				for (int nTile0 = 0; nTile0 < tileLayers[ml].length; nTile0++) if (tileLayers[ml][nTile0] > 0) {
					int nSurfTile = getSurfaceTileIndex(nTile0);
					if (tileData[nSurfTile] != null) {
						int l = tileData[nSurfTile].length;
						assigned_tiles[nTile0] = new boolean [l];
					}
					assigned_tiles[nTile0][tileLayers[ml][nTile0]-1] = true;
				}
			}
			return assigned_tiles;
		}		
		
		
		public int [][] enumerateClusters(
				final boolean [][] selection,
				final int          debugLevel,
				final int          dbg_X,
				final int          dbg_Y)
		{
			int [][] enum_clust = new int [selection.length][];
			for (int i = 0; i < enum_clust.length; i++){
				if (selection[i] != null) enum_clust[i] = new int [selection[i].length];
			}
			TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
			//			TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);
			int num_cluster = 1;
			for (int nTile0 = 0; nTile0 <selection.length; nTile0++) if (selection[nTile0] != null) {
				for (int nl0 = 0; nl0 < selection[nTile0].length; nl0++) if (selection[nTile0][nl0] && ((enum_clust[nTile0][nl0] == 0))) {
					ArrayList<Point> wave_list = new ArrayList<Point>();
					Point p = new Point(nTile0, nl0); // image tile number, surface index (0-based)
					enum_clust[p.x][p.y] = num_cluster;
					wave_list.add(p);
					while (wave_list.size() > 0){
						Point pt = wave_list.remove(0);
						int [] neibs = tileData[getSurfaceTileIndex(pt.x)][pt.y].getNeighbors();
						for (int dir  = 0; dir < tnImage.dirs; dir++) if (neibs[dir] >= 0){
							int nTile1 = tnImage.getNeibIndex(pt.x, dir);
							if (nTile1 >= 0) {
								int nl1 = neibs[dir];
								if ((selection[nTile1] != null) && selection[nTile1][nl1] && (enum_clust[nTile1][nl1] == 0) ){ // exists and not yet counted
									Point p1 = new Point(nTile1, nl1);
									enum_clust[nTile1][nl1] = num_cluster;
									wave_list.add(p1);
								}
							}
						}
					}
					num_cluster ++; 
				}
			}
			num_cluster --;
			if (debugLevel > -1 ) {
				System.out.println("enumerateClusters(): found " + num_cluster + " clusters");
			}
			return enum_clust;
		}
		
		public int getNumClusters(
				int [][] clusters)
		{
			int max_number = 0;
			for (int nTile = 0; nTile <clusters.length; nTile++) if (clusters[nTile] != null) {
				for (int nl = 0; nl < clusters[nTile].length; nl ++){
					if (clusters[nTile][nl] > max_number){
						max_number = clusters[nTile][nl];
					}
				}
			}
			return max_number;
			
		}
		
		public int [][] clusterStats(
				int [][]           clusters,
				final int          debugLevel,
				final int          dbg_X,
				final int          dbg_Y)
		{
			TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
			int max_number = getNumClusters(clusters);
			if (debugLevel > -1){
				System.out.println("clusterStats(): number of clusters: "+max_number);
			}
			int [][] stats0 = new int[max_number][CLUST_NUM_STATS]; // number of tiles, overlaps (each pair), conflicts-a (in dir, but connected to dif.), conflicts-b (no link)
			for (int nTile = 0; nTile <clusters.length; nTile++) if (clusters[nTile] != null) {
				for (int nl = 0; nl < clusters[nTile].length; nl ++){
					int nclust = clusters[nTile][nl];
					if (nclust > 0) {
						stats0[nclust - 1][CLUST_NUM_TILES]++;
						int num_same = 0;
						for (int nl1 = 0; nl1 < clusters[nTile].length; nl1 ++){
							if ((nl1 != nl) && (clusters[nTile][nl1] == nclust)){
								num_same ++;
							}
						}
						if (stats0[nclust - 1][CLUST_NUM_LAYERS] < (num_same + 1)){
							stats0[nclust - 1][CLUST_NUM_LAYERS] = (num_same + 1);
						}
						if (num_same > 0) {
							stats0[nclust - 1][CLUST_NUM_OVERLAPS]++;
							if (num_same > 1) {
								stats0[nclust - 1][CLUST_NUM_MULTI]++;
							}
						}
						if (num_same == 0) { // do not bother with overlapping tiles
							if (stats0[nclust - 1][CLUST_NUM_LAYERS] == 0){
								stats0[nclust - 1][CLUST_NUM_LAYERS] = 1;
							}
							int [] neibs = tileData[getSurfaceTileIndex(nTile)][nl].getNeighbors();
							for (int dir  = 0; dir < tnImage.dirs; dir++) if (neibs[dir] >= 0){
								int nTile1 = tnImage.getNeibIndex(nTile, dir);
								if ((nTile1 >= 0) && (clusters[nTile1] != null)) {
									int num_same1 = 0;
									for (int nl1 = 0; nl1 < clusters[nTile1].length; nl1 ++){
										if (clusters[nTile1][nl1] == nclust){
											num_same1 ++;
										}
									}
									if (num_same1 > 0){
										if (neibs[dir] <0 ){
											stats0[nclust - 1][CLUST_NUM_CONFLICTS_B]++;
										} else if (clusters[nTile1][neibs[dir]] != nclust){
											stats0[nclust - 1][CLUST_NUM_CONFLICTS_A]++;
										}
									}
								}
							}
						}
					}
				}			
			}
			// now sort tiles by the number of tiles in it
			ArrayList<Point> sort_list = new ArrayList<Point>(max_number);
			for (int nclust = 0 ; nclust < max_number; nclust++) {
				sort_list.add(new Point(stats0[nclust][CLUST_NUM_TILES], nclust));
			}
			Collections.sort(sort_list, new Comparator<Point>() {
				@Override
				public int compare(Point lhs, Point rhs) {
					// -1 - less than, 1 - greater than, 0 - equal, all inverted for descending
					return (lhs.x > rhs.x) ? -1 : (lhs.x < rhs.x) ? 1 : 0;
				}
			});
			int [][] stats = new int[max_number][CLUST_NUM_STATS]; // number of tiles, overlaps (each pair), conflicts-a (in dir, but connected to dif.), conflicts-b (no link)
			for (int i = 0; i < max_number; i++){
				stats[i] = stats0[sort_list.get(i).y];
				stats[i][CLUST_NUM_INDEX] = sort_list.get(i).y + 1; // zero-based -> 1-based
			}
			return stats;
		}

		public void showClusterStats(
				int [][] cluster_stats,
				int max_clusters){
			System.out.println("Detected "+cluster_stats.length+" clusters");
			if (max_clusters == 0){
				max_clusters = cluster_stats.length;
			}
			for (int i = 0; (i < cluster_stats.length) && (i < max_clusters); i++){
				System.out.print(i+": #"+ cluster_stats[i][CLUST_NUM_INDEX]);
				System.out.print(  " tiles: "+          cluster_stats[i][CLUST_NUM_TILES]);
				System.out.print(  " layers: "+         cluster_stats[i][CLUST_NUM_LAYERS]);
				
				System.out.print(  " overlaps: "+       cluster_stats[i][CLUST_NUM_OVERLAPS]);
				System.out.print(  " multi-overlaps: "+ cluster_stats[i][CLUST_NUM_MULTI]);
				System.out.print(  " coflicts A: "+ cluster_stats[i][CLUST_NUM_CONFLICTS_A]);
				System.out.print(  " coflicts B: "+ cluster_stats[i][CLUST_NUM_CONFLICTS_B]);
				System.out.println();
			}
		}

		public void showClusters(
				String title,
				int [][] cluster_stats,
				int max_clusters,
				int [][] clusters)
		{
			System.out.println("showClusters("+title+"," +max_clusters+"), numClusters = "+cluster_stats.length);

//			TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
			int num_slices = 0;
			int num_tiles = clusters.length;
			if (max_clusters == 0){
				max_clusters = cluster_stats.length;
			}
			if (max_clusters > cluster_stats.length){
				max_clusters = cluster_stats.length;
			}
			for (int i = 0; i < max_clusters; i++){
				num_slices += cluster_stats[i][CLUST_NUM_LAYERS];
			}
			double [][] img_clust = new double [num_slices][num_tiles];
			for (int i = 0; i < num_slices; i ++){
				for (int j = 0; j < num_tiles; j ++){
					img_clust[i][j] = Double.NaN;
				}
			}
			int [] rev_clust = new int [cluster_stats.length];
			for (int i = 0; i < cluster_stats.length; i ++){
				rev_clust[cluster_stats[i][CLUST_NUM_INDEX] - 1] = i;
			}
			
			String [] titles = new String [num_slices];
			int [] slice_first = new int [num_slices];
			int nslice = 0;
			for (int i = 0; i < max_clusters; i++){
				slice_first[i]=nslice;
				for (int j = 0; j < cluster_stats[i][CLUST_NUM_LAYERS]; j++){
					titles[nslice++] = "cluster-"+i+((cluster_stats[i][CLUST_NUM_LAYERS]>1)?("_"+j):"");
				}
			}
			
			for (int nTile = 0; nTile < num_tiles; nTile++) if (clusters[nTile] != null){
				int nSurfTile = getSurfaceTileIndex(nTile);
				boolean [] rend = new boolean [clusters[nTile].length];
				for (int i = 0; i < rend.length; i++) if (!rend[i] && (clusters[nTile][ i] > 0)) { // if (!rend[i] && (rev_clust[clusters[nTile][ i]] < max_clusters)){
					if (((clusters[nTile][ i] -1) >= rev_clust.length) || ((clusters[nTile][ i] -1) < 0)){
						System.out.println("showClusters() BUG:, clusters["+nTile+"]["+ i+"]-1 ("+(clusters[nTile][ i] -1)+") >= "+rev_clust.length);

					} else {
						if (rev_clust[clusters[nTile][ i] - 1] < max_clusters){
							int clust_index = clusters[nTile][ i] -1;
							int clust_layer = 0;
							for (int j = i; j < rend.length; j++) if ((clusters[nTile][ j] -1) == clust_index) {
								rend[j] = true;
								img_clust[slice_first[rev_clust[clust_index]] + clust_layer][nTile] = tileData[nSurfTile][j].getDisparity();
								clust_layer++;
							}
						}
					}
				}
			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays();
			sdfa_instance.showArrays(img_clust,  imageTilesX, imageTilesY, true, title, titles);
		}

		public boolean [][] growClusterOnce( // used from inside threads
				final boolean [][] cluster_local_pure,
				final int []       window,
				final int []       img_indices,				
				final int []       surf_indices,				
				final int          debugLevel,
				final int          dbg_X,
				final int          dbg_Y)
				
		{
			int num_tiles = img_indices.length;
			int width =  window[2];
			int height = window[3];
			boolean [][] grown_cluster = new boolean [num_tiles][];
			for (int neTile = 0; neTile < num_tiles; neTile++){
				if (cluster_local_pure[neTile] != null){
					grown_cluster[neTile] = cluster_local_pure[neTile].clone();
				}
			}
			int dbg_ntiles1 = 0;

			if (debugLevel > 0) {
				for (int neTile = 0; neTile < num_tiles; neTile++){
					if (grown_cluster[neTile] != null){
						for (int nl = 0; nl < grown_cluster[neTile].length; nl++){
							if (grown_cluster[neTile][nl]) dbg_ntiles1 ++;
						}
					}
				}
				System.out.println("growClusterOnce(): number of tiles = "+dbg_ntiles1);
			}
			final TileNeibs tnWindow =   new TileNeibs(width,height);

			for (int neTile0 = 0; neTile0 < num_tiles; neTile0++) if (cluster_local_pure[neTile0] != null){
				for (int nl0 = 0; nl0 < cluster_local_pure[neTile0].length; nl0++){
					if (cluster_local_pure[neTile0][nl0]){ // source should be single-layer, but ...
						int nSurfTile0 = surf_indices[neTile0];
						int [] neibs = tileData[nSurfTile0][nl0].getNeighbors(); 

						for (int dir = 0; dir <  tnWindow.dirs; dir++) {
							int nl1 = neibs[dir]; 
							if (nl1 < 0){
								if (debugLevel >-1) {
									System.out.println("growClusterOnce(): Expected 8 neighbors for tile nSurfTile0="+
											nSurfTile0+" neibs["+dir+"] = "+nl1);
								}
							} else {
								int neTile1 = tnWindow.getNeibIndex(neTile0, dir);
								if (neTile1 >= 0) {
									if (	((cluster_local_pure[neTile1] == null) || !cluster_local_pure[neTile1][nl1]) &&
											((grown_cluster[neTile1] == null) || !grown_cluster[neTile1][nl1])) {
										if (grown_cluster[neTile1] == null) {
											int nSurfTile1 = surf_indices[neTile1];
											// should never be null as it is connected from nSurfTile0
											grown_cluster[neTile1] = new boolean [tileData[nSurfTile1].length];
										}
										grown_cluster[neTile1][nl1] = true;
									}
								}
							}
						}
					}
				}
			}
			if (debugLevel > 0) {
				int dbg_ntiles2 = 0;
				for (int neTile = 0; neTile < num_tiles; neTile++){
					if (grown_cluster[neTile] != null){
						for (int nl = 0; nl < grown_cluster[neTile].length; nl++){
							if (grown_cluster[neTile][nl]) dbg_ntiles2 ++;
						}
					}
				}
				System.out.println("growClusterOnce(): new number of tiles = "+dbg_ntiles2+" (was "+dbg_ntiles1+")");
			}
			
			return grown_cluster;
		}
		
		public void setSurfaceData(
				final int [][] cluster_groups,
				final int [][] clusters_pure, // each tile/layer belongs to a single pure cluster
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{
			if (debugLevel >-1){
				System.out.println("setSurfaceData(): groups: "+cluster_groups.length);
			}
			final SurfaceData [] sdata = new SurfaceData[cluster_groups.length];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int [] tilesWH = {imageTilesX, imageTilesY};
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iGClust = ai.getAndIncrement(); iGClust < sdata.length; iGClust = ai.getAndIncrement()) {
							int dl = ((debugLevel > -1) && (iGClust == dbg_X)) ? 1:0; 
						
							int [] window = getClusterBBox(
									cluster_groups[iGClust], // final int []   nClust, // 1-based
									1,                       // final int      border,
									clusters_pure);          // final int [][] clusters)
							
							int [][] tile_indices =  getClusterBBoxIndices(
									window,                  // final int [] window,
									0);                      // border); // final int    border) // maybe 0, actual value just saves time
							int [] img_tile_indices = tile_indices[0]; 
							int [] surf_tile_indices = tile_indices[1];
							boolean [][] combined = new boolean [img_tile_indices.length][];
							if (dl > 0){
								String dbg_s = "";
								for (int indxClust = 0; indxClust < cluster_groups[iGClust].length; indxClust++){
									dbg_s += " " +cluster_groups[iGClust][indxClust];
								}
								System.out.println("setSurfaceData(): iGClust = "+iGClust+" window = [ "+window[0]+", "+window[1]+", "+window[2]+", "+window[3]+" ],  clusters: ["+dbg_s+"]");
							}
							for (int indxClust = 0; indxClust< cluster_groups[iGClust].length; indxClust++){
								if (dl > 0){
									System.out.println("setSurfaceData():  cluster_groups["+iGClust+"]["+indxClust+"] = "+cluster_groups[iGClust][indxClust]);
								}
								boolean [][] extracted_cluster = extractCluster(
										cluster_groups[iGClust][indxClust],   // final int      numToSplit, // 1-based
										clusters_pure,       // final int [][] clusters,
										window,              // final int []   window,
										img_tile_indices,    // final int []   bbox_indices,				
										debugLevel,          // final int      debugLevel,
										dbg_X,               // final int      dbg_X,
										dbg_Y);              // final int      dbg_Y)
								for (int neTile = 0; neTile < combined.length; neTile++){
									if (extracted_cluster[neTile] != null){
										if (combined[neTile] == null){
											combined[neTile] = extracted_cluster[neTile].clone();
										} else {
											for (int nl = 0; nl < extracted_cluster[neTile].length; nl++){
												combined[neTile][nl] |= extracted_cluster[neTile][nl];
											}
										}
									}
								}
							}
							boolean [] pure_sel = new boolean [combined.length];
							for (int neTile = 0; neTile < combined.length; neTile++){
								if (combined[neTile] != null){
									for (int nl = 0; nl < combined[neTile].length; nl++){
										if (combined[neTile][nl]) {
											pure_sel[neTile] = true;
											break;
										}
									}									
								}
							}
							boolean [][] grown_combined = growClusterOnce( // used from inside threads
									combined,   // final boolean [][] cluster_local_pure,
									window,            // final int []   window,
									img_tile_indices,  // final int []   bbox_indices,				
									surf_tile_indices, // final int []       surf_indices,				
									0, // dl, // debugLevel,       // final int      debugLevel,
									dbg_X,            // final int      dbg_X,
									dbg_Y);           // final int      dbg_Y)
							boolean [] selected = new boolean [combined.length];
							boolean [] border  =  new boolean [combined.length];
							for (int neTile = 0; neTile < grown_combined.length; neTile++){
								if (grown_combined[neTile] != null){
									for (int nl = 0; nl < grown_combined[neTile].length; nl++){
										if (grown_combined[neTile][nl]) {
											selected[neTile] = true;
											if (!pure_sel[neTile]){
												border[neTile] = true;
											}
											break;
										}
									}									
								}
							}
							sdata[iGClust] = new 	SurfaceData(
									tilesWH,
									window);
							sdata[iGClust].setSelected(selected);
							sdata[iGClust].setBorder(border);
							sdata[iGClust].setCluster(cluster_groups[iGClust]);
							for (int neTile = 0; neTile < grown_combined.length; neTile++){
								if (selected[neTile]){
									for (int nl = 0; nl < grown_combined[neTile].length; nl++){
										if (grown_combined[neTile][nl]){
											sdata[iGClust].setDisparity(
													neTile,
													tileData[surf_tile_indices[neTile]][nl].getDisparity());
										}
									}
								}
							}
							if (dl > 0){
								System.out.println("setSurfaceData(): sdata ["+iGClust+"] prepared");
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ArrayList<SurfaceData> sdata_list = new ArrayList<SurfaceData>();
			for (int i = 0; i < sdata.length; i++){
				sdata_list.add(sdata[i]);
			}
			Collections.sort(sdata_list, new Comparator<SurfaceData>() {
				@Override
				public int compare(SurfaceData lhs, SurfaceData rhs) {
					// -1 - less than, 1 - greater than, 0 - equal, all inverted for descending
					return (lhs.getSize(true) > rhs.getSize(true)) ? -1 : (lhs.getSize(true) < rhs.getSize(true)) ? 1 : 0;
				}
			});
			this.surfaceData = 	sdata_list.toArray(new SurfaceData[0] );
			if (debugLevel >-1){
				if (debugLevel > 0){
					for (int i = 0; i < this.surfaceData.length; i++) {
						System.out.print(i+": "+this.surfaceData[i].getSize(true)+" ("+this.surfaceData[i].getSize(false)+") [");
						for (int j = 0; j < this.surfaceData[i].getCluster().length; j++){
							System.out.print(this.surfaceData[i].getCluster()[j] +", ");
						}
						System.out.println("]");
					}
				}
				System.out.println("setSurfaceData(): DONE");
			}
 		}
		
		public SurfaceData getSurfaceData(int indx){
			return this.surfaceData[indx];
		}

		public int getSurfaceDataLength(){
			return this.surfaceData.length;
		}
		
		public int [][] mergeNoConflict(
				final int [][] matchedGrown,
				final int [][] clusters_grown,
				final int [][] clusters_pure,
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{
			if (debugLevel > 0){
				System.out.println("mergeNoConflict(): groups: "+matchedGrown.length);
			}
			final int num_grown = matchedGrown.length;
			final int [][][] pre_merged_subs = new int [num_grown][][];
			
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iGClust = ai.getAndIncrement(); iGClust < num_grown; iGClust = ai.getAndIncrement()) {
							int num_subs = matchedGrown[iGClust].length;
							if (num_subs < 2){ // nothing to split
								pre_merged_subs[iGClust] = new int [num_subs][];
								if (num_subs > 0) { // should always be so?
									pre_merged_subs[iGClust][0] = matchedGrown[iGClust].clone();
								}
							} else { //>= 2 clusters
								int nGClust = iGClust + 1;
								int dl = ((debugLevel > -1) && (num_subs > 3)) ? 3: 0;
								if (dl > 3){
									System.out.println("mergeNoConflict(): nGClust: "+nGClust+" num_subs= "+num_subs);
								}
								
								// create window for the grown cluster, it will include all subs
								int [] window = getClusterBBox(
										nGClust, // 1-based
										0, // border,
										clusters_grown);

								int [][] tile_indices =  getClusterBBoxIndices(
										window,  // final int [] window,
										0);      // border); // final int    border) // maybe 0, actual value just saves time
								int [] img_tile_indices = tile_indices[0]; 
								int [] surf_tile_indices = tile_indices[1];
								int num_tiles = img_tile_indices.length;
								int [] clust_sizes = new int [num_subs];
								boolean [][][] subs_pure =  new boolean [num_subs][][];
								boolean [][][] subs_grown = new boolean [num_subs][][];
								for (int nSub = 0; nSub < num_subs; nSub++){
									subs_pure[nSub] = extractCluster(
											matchedGrown[iGClust][nSub],   // final int      numToSplit, // 1-based
											clusters_pure,    // final int [][] clusters,
											window,           // final int []   window,
											img_tile_indices, // final int []   bbox_indices,				
											debugLevel,       // final int      debugLevel,
											dbg_X,            // final int      dbg_X,
											dbg_Y);           // final int      dbg_Y)
									
									subs_grown[nSub] = growClusterOnce( // used from inside threads
											subs_pure[nSub],   // final boolean [][] cluster_local_pure,
											window,            // final int []   window,
											img_tile_indices,  // final int []   bbox_indices,				
											surf_tile_indices, // final int []       surf_indices,				
											0, // dl, // debugLevel,       // final int      debugLevel,
											dbg_X,            // final int      dbg_X,
											dbg_Y);           // final int      dbg_Y)
								}
								// now build sets of clusters
								// a) connected to current (and having higher number)
								// b) conflicting with current (and having higher number)
								ArrayList<HashSet<Integer>> touching_list = new ArrayList<HashSet<Integer>>(); 
								ArrayList<HashSet<Integer>> conflict_list = new ArrayList<HashSet<Integer>>();
								for (int nSub = 0; nSub < num_subs; nSub++){
									touching_list.add(new HashSet<Integer>());
									conflict_list.add(new HashSet<Integer>());
								}
								for (int neTile = 0; neTile < num_tiles; neTile++){
									for (int nSub = 0; nSub < num_subs; nSub++){
										if (subs_grown[nSub][neTile] != null){
											for (int nl1 = 0; nl1 < subs_grown[nSub][neTile].length; nl1++){
												if (subs_grown[nSub][neTile][nl1]){
													clust_sizes[nSub] ++;
													for (int nSub2 = nSub+1; nSub2 < num_subs; nSub2++) if (subs_grown[nSub2][neTile] != null){
														for (int nl2 = 0; nl2 < subs_grown[nSub2][neTile].length; nl2++){ //*
															if (subs_grown[nSub2][neTile][nl2]){
																if (nl1 == nl2) {
																	touching_list.get(nSub).add(new Integer(nSub2));
																	touching_list.get(nSub2).add(new Integer(nSub)); // both ways
																} else {
																	conflict_list.get(nSub).add(new Integer(nSub2));
																	conflict_list.get(nSub2).add(new Integer(nSub)); // both ways
																}
															}
														}
													}
												}
											}
										}
									}
								}
								
								ArrayList<HashSet<Integer>> groups_list = new ArrayList<HashSet<Integer>>();
								HashSet<Integer> clusters_left = new HashSet<Integer>();
								for (int nSub = 0; nSub < num_subs; nSub++){
									clusters_left.add(new Integer(nSub));
								}
								while (!clusters_left.isEmpty()){
									HashSet<Integer> new_group = new HashSet<Integer>();
									HashSet<Integer> candidates = new HashSet<Integer>();
									HashSet<Integer> conflicts = new HashSet<Integer>();
//									HashSet<Integer> touching = new HashSet<Integer>();
									
									// start with the largest of the remaining clusters
									Integer best_sub = -1;
									for (Integer sc:clusters_left){
										if ((best_sub < 0) || (clust_sizes[best_sub] < clust_sizes[sc])){
											best_sub = sc;
										}
									}
									if ( !clusters_left.remove(best_sub)) {
										System.out.println("mergeNoConflict() bug: can not remove"+best_sub);
									}
									new_group.add(best_sub);
									conflicts.addAll(conflict_list.get(best_sub));    // keep track of all accumulated conflicts
									candidates.addAll(touching_list.get(best_sub));   // add all clusters that are touching the current selection
									candidates.removeAll(conflicts);                  // remove all conflicts
									
									while (!candidates.isEmpty()) { // add more clusters if possible
										// Find the largest one
										best_sub = -1;
										for (Integer sc:candidates){
											if ((best_sub < 0) || (clust_sizes[best_sub] < clust_sizes[sc])){
												best_sub = sc;
											}
										}
										// add that new cluster
										clusters_left.remove(best_sub);
										new_group.add(best_sub);
										conflicts.addAll(conflict_list.get(best_sub));    // keep track of all accumulated conflicts
										candidates.addAll(touching_list.get(best_sub));   // add all clusters that are touching the current selection
										candidates.removeAll(conflicts);                  // remove all conflicts
										candidates.removeAll(new_group);                  // remove all what is already included
									}
									groups_list.add(new_group);
								}
								pre_merged_subs[iGClust] = new int [groups_list.size()][];
								for (int ng = 0; ng < groups_list.size(); ng++) {
									pre_merged_subs[iGClust][ng] = new int [groups_list.get(ng).size()];
									int nc = 0;
									for (Integer cl: groups_list.get(ng)){
										pre_merged_subs[iGClust][ng][nc++] = matchedGrown[iGClust][cl]; // cl;
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			// "flatten" pre_merged_subs
			int num_new_groups = 0;
			for (int ng = 0; ng < pre_merged_subs.length; ng ++) {
				num_new_groups += pre_merged_subs[ng].length;
			}
			if (debugLevel >-1){
				System.out.println("mergeNoConflict(): groups: " + num_new_groups + " (was : "+matchedGrown.length+")");
			}
			final int [][] merged_subs = new int [num_new_groups][];
			int indx = 0;
			for (int ng = 0; ng < pre_merged_subs.length; ng ++) {
				for (int nc = 0; nc < pre_merged_subs[ng].length; nc ++){
					merged_subs[indx++] = pre_merged_subs[ng][nc];
				}
			}
			
			if (debugLevel > 0) {
				for (int i = 0; i < num_new_groups; i++){
					System.out.print("mergeNoConflict(): "+ (i+1)+" [");
					for (int j = 0; j < merged_subs[i].length; j++) {
						if (j > 0) System.out.print(", ");
						System.out.print(merged_subs[i][j]);
					}
					System.out.println(" ]");
				}
			}
			
			return merged_subs;
			
		}
		
		
		/**
		 * Grow each of the clusters (encoded as positive cluster numbers per tile per layer) by 1 in each
		 * of 8 directions. As they now may overlap they are encoded in boolean array [cluster][tile][layer]
		 * @param clusters_pure cluster numbers (>0) for each tile, each layer
		 * @param debugLevel
		 * @param dbg_X
		 * @param dbg_Y
		 * @return boolean array [cluster-1][tile][layer]
		 */
		public boolean [][][] growEachCluster (
				final int [][] clusters_pure,
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{
			final int num_clusters = getNumClusters(clusters_pure);
			final boolean [][][] clusters = new boolean [num_clusters][clusters_pure.length][];
			final TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);

			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iClust = ai.getAndIncrement(); iClust < num_clusters; iClust = ai.getAndIncrement()) {
							int nClust = iClust + 1; // 1-based
							for (int nTile0 = 0; nTile0 < clusters_pure.length; nTile0++) if (clusters_pure[nTile0] != null){
								for (int nl0 = 0; nl0 < clusters_pure[nTile0].length; nl0++){
									if (clusters_pure[nTile0][nl0] == nClust){ // source should be single-layer, but ...
										if (clusters[iClust][nTile0]==null) {
											clusters[iClust][nTile0] = new boolean [clusters_pure[nTile0].length];
										}
										clusters[iClust][nTile0][nl0] = true;
										int nSurfTile0 = getSurfaceTileIndex(nTile0);
										int [] neibs = tileData[nSurfTile0][nl0].getNeighbors(); 

										for (int dir = 0; dir <  tnImage.dirs; dir++) {
											int nl1 = neibs[dir]; 
											if (nl1 < 0){
												if (debugLevel >-1) {
													System.out.println("growEachCluster(): Expected 8 neighbors for tile nSurfTile0="+
															nSurfTile0+" neibs["+dir+"] = "+nl1);
												}
												int nTile1 = tnImage.getNeibIndex(nTile0, dir);
												if (nTile1 >= 0) {
													if (	((clusters_pure[nTile1] == null) || (clusters_pure[nTile1][nl1] != nClust)) &&
															((clusters[iClust][nTile1] == null) || !clusters[iClust][nTile1][nl1])) {
														if (clusters[iClust][nTile1] == null) {
															int nSurfTile1 = getSurfaceTileIndex(nTile1);
															clusters[iClust][nTile1] =  new boolean [tileData[nSurfTile1].length];
														}
														clusters[iClust][nTile1][nl1] = true;
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
			
			
			return clusters;
		}
		
		
		/**
		 * Match grown/merged clusters to the ones they were made of
		 * @param clusters_grown per-tile, per layer array of positive merged cluster numbers (0 - empty)
		 * @param clusters_pure per-tile, per layer array of positive pure cluster numbers (0 - empty)
		 * @param debugLevel
		 * @param dbg_X
		 * @param dbg_Y
		 * @return array for each of the merged clusters array of pure cluster numbers it is made of
		 */
		public int [][] matchPureGrown (
				final int [][] clusters_grown,
				final int [][] clusters_pure,
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{
			final int num_grown = getNumClusters(clusters_grown);
			final int [][] grown_subs = new int [num_grown][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iparent = ai.getAndIncrement(); iparent < num_grown; iparent = ai.getAndIncrement()) {
							int num_parent = iparent + 1; // 1-based
							HashSet<Integer> subs_set = new HashSet<Integer>();
							for (int nTile = 0; nTile < clusters_grown.length; nTile++) if (clusters_grown[nTile] != null) {
								for (int nl = 0; nl < clusters_grown[nTile].length; nl++) if (clusters_grown[nTile][nl] == num_parent){
									if (clusters_pure[nTile] != null){
										int iclust = clusters_pure[nTile][nl];
										if (iclust > 0){
											subs_set.add(new Integer(iclust));
										}
									}
								}
							}
							grown_subs[iparent] = new int [subs_set.size()];
							int indx = 0;
							for (Integer s:subs_set){
								grown_subs[iparent][indx++] = s;
							}
							Arrays.sort(grown_subs[iparent]);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (debugLevel > 0) {
				for (int i = 0; i < num_grown; i++){
					System.out.print("matchPureGrown(): "+ (i+1)+" [");
					for (int j = 0; j < grown_subs[i].length; j++) {
						if (j > 0) System.out.print(", ");
						System.out.print(grown_subs[i][j]);
					}
					System.out.println(" ]");
				}
			}
			return grown_subs;
		}
		
		
		public int [] getClusterBBox(
				final int      nClust, // 1-based
				final int      border,
				final int [][] clusters)
		{
//			if (numToSplit == 74) { //177){
//				System.out.println("getClusterBBox() numToSplit="+numToSplit);
//			}
			int x_max = -1, y_max = -1, x_min = imageTilesX, y_min = imageTilesY;
			for (int nTile = 0; nTile < clusters.length; nTile++) if (clusters[nTile] != null){
				for (int nl = 0; nl < clusters[nTile].length; nl++){
					if (clusters[nTile][nl] == nClust){
						int tX = nTile % imageTilesX;
						
						int tY = nTile / imageTilesX;
						if (tX < x_min) x_min = tX;
						if (tX > x_max) x_max = tX;
						if (tY < y_min) y_min = tY;
						if (tY > y_max) y_max = tY;
					}
				}
			}
			int x0 = x_min - border; // may be negative
			int y0 = y_min - border; // may be negative
			int width = x_max +1 + border - x0; 
			int height = y_max +1 + border - y0;
			int [] bbox = {x0, y0, width, height};
			return bbox;
		}		
		public int [] getClusterBBox(
				final int []   nClust, // 1-based
				final int      border,
				final int [][] clusters)
		{
			int x_max = -1, y_max = -1, x_min = imageTilesX, y_min = imageTilesY;
			for (int nTile = 0; nTile < clusters.length; nTile++) if (clusters[nTile] != null){
				for (int nl = 0; nl < clusters[nTile].length; nl++){
					boolean sel = false;
					for (int i = 0; i < nClust.length; i++){
						if (clusters[nTile][nl] == nClust[i]){
							sel = true;
							break;
						}
					}
					if (sel){
						int tX = nTile % imageTilesX;
						int tY = nTile / imageTilesX;
						if (tX < x_min) x_min = tX;
						if (tX > x_max) x_max = tX;
						if (tY < y_min) y_min = tY;
						if (tY > y_max) y_max = tY;
					}
				}
			}
			int x0 = x_min - border; // may be negative
			int y0 = y_min - border; // may be negative
			int width = x_max +1 + border - x0; 
			int height = y_max +1 + border - y0;
			int [] bbox = {x0, y0, width, height};
			return bbox;
		}		

		
		
		
		public int [][] getClusterBBoxIndices(
				final int [] window,
				final int    border) // maybe 0, actual value just saves time
		{
			int x0 =     window[0];
			int y0 =     window[1];
			int width =  window[2];
			int height = window[3];
			int [][] bbox_indices = new int [2][width*height]; // [0] - image index, [1] - surface index
			int wb = width - border;
			int hb = height - border;
			
			
			for (int neTile = 0; neTile < bbox_indices[0].length; neTile++){
				bbox_indices[0][neTile] = -1;
				bbox_indices[1][neTile] = -1;
				int bx = neTile % width;
				int by = neTile / width;
				if ((by >= border) && (bx >= border) && (by < hb) && (bx < wb)){
					int tx = bx + x0;
					int ty = by + y0;
					if ((tx >= 0) && (ty >= 0) && (tx < imageTilesX) && (ty < imageTilesY)){
						bbox_indices[0][neTile] = ty * imageTilesX + tx;
						bbox_indices[1][neTile] = ty * stilesX *superTileSize + tx;
					}				
				}
			}
			return bbox_indices;
		}
		
		
		
		public boolean [][] extractCluster(
				final int      numToSplit, // 1-based
				final int [][] clusters,
				final int []   window,
				final int []   bbox_indices,				
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{
			boolean [][] extracted_cluster = new boolean [bbox_indices.length][];
			if (numToSplit == 74) { //177){
				System.out.println("extractCluster() numToSplit="+numToSplit);
			}
			
			for (int neTile = 0; neTile < extracted_cluster.length; neTile++){
				int nTile = bbox_indices[neTile];
				if (nTile >= 0){
					if (clusters[nTile] != null) {
						boolean has_cluster = false;
						for (int nl = 0; nl< clusters[nTile].length; nl++){
							if (clusters[nTile][nl] == numToSplit){
								has_cluster =true;
								break;
							}
						}
						if (has_cluster){
							extracted_cluster[neTile] = new boolean[clusters[nTile].length];
							for (int nl = 0; nl< clusters[nTile].length; nl++){
								extracted_cluster[neTile][nl] = (clusters[nTile][nl] == numToSplit); 
							}
						}
					}
					
				}
			}
			return extracted_cluster;
		}
		

		public int [][] spitConflictClusters(
				final int [][] clusters,
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{
			final int num_clusters = getNumClusters(clusters);
			if (debugLevel > -1){
				System.out.println("spitConflictClusters(): number of original clusters = " + num_clusters);
			}
//			final int [][]   split_windows =  new int [num_clusters][];   // x,y,w,h for the cluster that is split 
			final int [][][]   split_indices =  new int [num_clusters][2][];   // for each tile in the window - number of the corresponding image tile or -1
			final int [][][] split_clusters = new int [num_clusters][][]; // for each tile in the window, a stack of layers, each number of a subcluster or 0
			final int []     split_number =   new int [num_clusters];     // number of sublusters the parent was split to
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int numToSplit = ai.getAndIncrement(); numToSplit < num_clusters; numToSplit = ai.getAndIncrement()) {
							if ((numToSplit + 1) == 74) { //177){
								System.out.println("extractCluster() numToSplit="+numToSplit);
							}

							int [] window = getClusterBBox(
									numToSplit + 1, // 1-based
									0, // border,
									clusters);

							split_indices[numToSplit] =  getClusterBBoxIndices(
									window,  // final int [] window,
									0);      // border); // final int    border) // maybe 0, actual value just saves time

							boolean [][] extracted_cluster = extractCluster(
									numToSplit + 1,   // final int      numToSplit, // 1-based
									clusters,     // final int [][] clusters,
									window,       // final int []   window,
									split_indices[numToSplit][0], // final int []   bbox_indices,				
									debugLevel,   // final int      debugLevel,
									dbg_X,        // final int      dbg_X,
									dbg_Y);       // final int      dbg_Y)

							split_clusters[numToSplit] =  splitCluster(
									numToSplit + 1,        // final int      numToSplit, // 1-based
									clusters,          // final int [][] clusters,
									window,            // final int []   window,
									split_indices[numToSplit][0],      // final int []   bbox_indices,				
									split_indices[numToSplit][1],      // final int []   bbox_surf_indices,		
									extracted_cluster, // final boolean [][] extracted_cluster,
									debugLevel,        // final int      debugLevel,
									dbg_X,             // final int      dbg_X,
									dbg_Y);            // final int      dbg_Y)
							
							split_number[numToSplit] = 	getNumClusters(split_clusters[numToSplit]);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			final int [] new_cluster_index = new int [num_clusters];
			for (int i = 0; i < num_clusters; i++) {
				if (i == 0) new_cluster_index[i] = 0;
				else new_cluster_index[i] = new_cluster_index[i-1] + split_number[i-1];
			}
			int num_new = new_cluster_index[num_clusters - 1] + split_number[num_clusters - 1];
			if (debugLevel > -1){
				System.out.println("spitConflictClusters(): number of new clusters = " + num_new + " ( was " + num_clusters + ")");
			}
			final int [][] new_clusters = new int [clusters.length][];
			for (int i = 0; i < new_clusters.length; i++) if (clusters[i] != null) new_clusters[i] = new int [clusters[i].length];
			
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int numToSplit = ai.getAndIncrement(); numToSplit < num_clusters; numToSplit = ai.getAndIncrement()) {
							int numToSplit1 = numToSplit+1;
							if (numToSplit1 == 74) { //177){
								System.out.println("extractCluster() numToSplit1="+numToSplit1);
							}
							
							for (int neTile = 0; neTile < split_indices[numToSplit][0].length; neTile++) {
								int nTile = split_indices[numToSplit][0][neTile] ;
								if ((nTile >= 0) && (clusters[nTile] != null)){
									for (int nl = 0; nl < clusters[nTile].length; nl++){
										if (clusters[nTile][nl] == numToSplit1) {
											new_clusters[nTile][nl] = new_cluster_index[numToSplit] + split_clusters[numToSplit][neTile][nl];
										}
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (debugLevel > -1){
				System.out.println("spitConflictClusters(): DONE");
			}
			return new_clusters;
		}
		
		/**
		 * 
		 * @param numToSplit
		 * @param clusters
		 * @param window
		 * @param bbox_indices
		 * @param bbox_surf_indices
		 * @param extracted_cluster
		 * @param debugLevel
		 * @param dbg_X
		 * @param dbg_Y
		 * @return for each tile in the window,for each layer in that location number of a subcluster (1-based) or 0
		 */
		public int [][] splitCluster(
				final int      numToSplit, // 1-based
				final int [][] clusters,
				final int []   window,
				final int []   bbox_indices,
				final int []   bbox_surf_indices,
				final boolean [][] extracted_cluster,
				final int      debugLevel,
				final int      dbg_X,
				final int      dbg_Y)
		{

			int width =  window[2];
			int height = window[3];
			int [][] split_cluster = new int [extracted_cluster.length][];
			for (int neTile = 0; neTile < extracted_cluster.length; neTile++){
				if (extracted_cluster[neTile] != null){
					split_cluster[neTile] = new int [extracted_cluster[neTile].length];
				}
			}
			final TileNeibs tnCluster =   new TileNeibs(width, height);
			if (numToSplit == 74) { //177){
				int num_tiles = 0;
				for (int neTile = 0; neTile < extracted_cluster.length; neTile++){
					if (extracted_cluster[neTile] != null){
						for (int i = 0; i < extracted_cluster[neTile].length; i++) {
							if (extracted_cluster[neTile][i]) num_tiles++;
						}
					}
				}
				
				System.out.println("splitCluster() numToSplit="+numToSplit+" num tiles = "+ num_tiles);
				
			}
			
			int num_cluster = 1;
			while (true) {
				
				int neTile0 = 0;
				int nl0 = 0;
				label_loop : {
					for (; neTile0 < extracted_cluster.length; neTile0++){
						if (extracted_cluster[neTile0] != null){
							for (nl0 = 0; nl0 < extracted_cluster[neTile0].length; nl0++) if (extracted_cluster[neTile0][nl0]){
								break label_loop;
							}
						}
					}
				}
				if (neTile0 >= extracted_cluster.length){
					break; // nothing left
				}
				// run wave from neTile0, nl0 checking for conflicts, filling
				ArrayList<Point> clust_list = new ArrayList<Point>();
				{
					Point p = new Point(neTile0, nl0);
					split_cluster[p.x][p.y] = num_cluster;
					extracted_cluster[p.x][p.y] = false;
					clust_list.add(p);
				}
				int dbg_size = 1;
				while (!clust_list.isEmpty()){
					Point p = clust_list.remove(0);
					int neTile = p.x;
					int nl = p.y;
					int nSurfTile = bbox_surf_indices[neTile]; // WRONG - fixed
					int [] neibs = tileData[nSurfTile][nl].getNeighbors(); 
					for (int dir = 0; dir < tnCluster.numNeibs(); dir++) {
						int neTile1 = tnCluster.getNeibIndex(neTile,dir);
						if ((neTile1 >= 0) && (extracted_cluster[neTile1] != null)) { // it should have original tile there, not yet assigned to a new
							// check it is empty on all levels
							boolean tileEmpty = true;
							if (split_cluster[neTile1] != null) {
								for (int i = 0; i < split_cluster[neTile1].length; i++ ){
									if (split_cluster[neTile1][i] == num_cluster){
										tileEmpty = false;
										break;
									}
								}
							}
							if (tileEmpty) {
								if (neibs[dir] < 0){ // should not happen
									if (debugLevel > 0) {
										System.out.println("splitCluster(): existing nSurfTile="+nSurfTile+" surface is not connected in direction "+dir);
									}
								} else {
									int nl1 = neibs[dir];
									if (extracted_cluster[neTile1][nl1]) {
										// verify it does not have connections from assigned tiles to different levels
										label_no_conflict: {
											for (int dir1 = 0; dir1 < tnCluster.numNeibs(); dir1++) {
												int neTile2 = tnCluster.getNeibIndex(neTile1,dir1);
												if ((neTile2 >= 0) && (split_cluster[neTile2] != null)) {
													for (int nl2 = 0; nl2 < split_cluster[neTile2].length; nl2++){
														if (split_cluster[neTile2][nl2] == num_cluster){
															int nl_this = tileData[bbox_surf_indices[neTile2]][nl2].getNeighbor(tnCluster.opposite(dir1));
															if ((nl_this >= 0) && (nl_this != nl1)){
																break label_no_conflict;
															}
														}
													}
												}
											}
											// no conflicts found, add this tile
											extracted_cluster[neTile1][nl1] = false;
											split_cluster[neTile1][nl1] = num_cluster;
											clust_list.add(new Point(neTile1, nl1));
											dbg_size ++;
										}
									}
								}
							}
						}
					}
					
				}
				if (numToSplit == 74) { //177){
					System.out.println("splitCluster() numToSplit="+numToSplit+" num_cluster="+num_cluster+" num tiles="+dbg_size);
				}

				num_cluster++;
			}
			if (numToSplit == 74) { //177){
				int num_tiles = 0, num_tiles1= 0;
				for (int neTile = 0; neTile < extracted_cluster.length; neTile++){
					if (extracted_cluster[neTile] != null){
						for (int i = 0; i < extracted_cluster[neTile].length; i++) {
							if (extracted_cluster[neTile][i]) num_tiles++;
							if (split_cluster[neTile][i] > 0) num_tiles1++;
							if (extracted_cluster[neTile][i]) {
								System.out.println("neTile="+neTile+", nl = "+i);
							}
						}
					}
				}
				
				System.out.println("splitCluster() return  numToSplit="+numToSplit + " num tiles = "+ num_tiles+ " num tiles1 = "+ num_tiles1);
			}
			return split_cluster;
		}
		
		public boolean [][] growSelectionOnce(
//				int grow,
				final boolean [][] sel_in,
				final int          debugLevel,
				final int          dbg_X,
				final int          dbg_Y)
		{
			final boolean [][] selection =  new boolean [sel_in.length][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
			//			final TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);
			ai.set(0);
			// make sure even empty tiles have surface selection arrays defined
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < sel_in.length; nTile = ai.getAndIncrement()) {
							boolean not_empty = false;
							boolean [] stack = null;
							if (sel_in[nTile] != null) {
								stack = sel_in[nTile].clone(); 
							}
							int nSurfTile = getSurfaceTileIndex(nTile);
							if (tileData[nSurfTile] != null) {
								if (stack == null) {
									stack = new boolean[tileData[nSurfTile].length];
								}
								for (int nl = 0; nl < stack.length; nl++){
									if (stack[nl]) {
										not_empty = true;
									} else {
										int [] neibs = tileData[nSurfTile][nl].getNeighbors();
										for (int dir = 0; dir < neibs.length; dir++){
											if (neibs[dir] >= 0){
												int nTile1 = tnImage.getNeibIndex(nTile, dir);
												if ((nTile1 >= 0) && (sel_in[nTile1] != null) && sel_in[nTile1][neibs[dir]]){
													stack[nl] = true;
													not_empty = true;
												}
											}
										}
									}
								}								
							}
							if (not_empty){
								selection[nTile] = stack;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			return selection;
		}

		
		
		public boolean [][] growSelection(
				int grow,
				final boolean [][] sel_in,
				final int          debugLevel,
				final int          dbg_X,
				final int          dbg_Y)
		{
			final int DIR_RIGHT = 2;
			final int DIR_LEFT =  6;
			final int DIR_UP =    0;
			final int DIR_DOWN =  4;
			final int [] DIRS = {DIR_LEFT, DIR_RIGHT, DIR_UP, DIR_DOWN};  
			final boolean [][] selection =   sel_in.clone();
			final boolean [][] selection_tmp = sel_in.clone();
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final TileNeibs tnImage =   new TileNeibs(imageTilesX, imageTilesY);
			//			final TileNeibs tnSurface = new TileNeibs(stilesX * superTileSize, stilesY * superTileSize);
			ai.set(0);
			// make sure even empty tiles have surface selection arrays defined
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < sel_in.length; nTile = ai.getAndIncrement()) {
							int nSurfTile = getSurfaceTileIndex(nTile);
							if (tileData[nSurfTile] != null) {
								if (sel_in[nTile] != null){
									selection[nTile] = sel_in[nTile].clone();
								} else {
									selection[nTile] = new boolean[tileData[nSurfTile].length];
								}
								//								selection_tmp[nTile] = new boolean[tileData[nSurfTile].length];
							}								
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			for (; grow > 0; grow -=2){
				for (int dri = 0; dri < DIRS.length; dri++ ){
					final int fdri = dri;
					ai.set(0);
					// extend to the right
					for (int ithread = 0; ithread < threads.length; ithread++) {
						threads[ithread] = new Thread() {
							public void run() {
								for (int nTile = ai.getAndIncrement(); nTile < sel_in.length; nTile = ai.getAndIncrement()) {
									if (selection[nTile] != null) {
										if (fdri == 0) selection_tmp[nTile] = new boolean[selection[nTile].length]; // new, empty
										int nTile1 = tnImage.getNeibIndex(nTile, DIRS[fdri]);
										if (nTile1 >= 0){
											int nSurfTile = getSurfaceTileIndex(nTile);
											for (int nl = 0; nl < selection_tmp[nTile].length; nl ++) {
												if (tileData[nSurfTile][nl] != null){
													int nl1 = tileData[nSurfTile][nl].getNeighbor(DIRS[fdri]);
													if ((nl1 >=0 ) && selection[nTile1][nl1]){
														selection_tmp[nTile][nl] = true;
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

					if (((dri == 1) &&  (grow > 1)) || (dri == 3)){ // always after both directions, after horizontal if diagonal is enabled
						// copy results back
						ai.set(0);
						// extend to the left
						for (int ithread = 0; ithread < threads.length; ithread++) {
							threads[ithread] = new Thread() {
								public void run() {
									for (int nTile = ai.getAndIncrement(); nTile < sel_in.length; nTile = ai.getAndIncrement()) {
										if (selection_tmp[nTile] != null){
											for (int nl = 0; nl < selection_tmp[nTile].length; nl ++) {
												if (selection_tmp[nTile][nl]) {
													selection[nTile][nl] = true;
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
			}
			return selection;
		}
		
		public void mergeAndGrow( // TODO: add result
				final boolean [][] sel_in,
				final int          debugLevel,
				final int          dbg_X,
				final int          dbg_Y)
		{
			System.out.println("mergeAndGrow() start");
			boolean [][] assigned_sel = extractSelection(
					debugLevel,
					dbg_X,
					dbg_Y);
			
			int [][] clusters_pure = enumerateClusters(
					assigned_sel, //final boolean [][] selection,
					debugLevel,
					dbg_X,
					dbg_Y);
			clusters_pure = spitConflictClusters(
					clusters_pure, // final int [][] clusters,
					debugLevel,
					dbg_X,
					dbg_Y);

			boolean [][] grown_sel = growSelectionOnce(
//					2,              // int grow,
					assigned_sel,   // final boolean [][] sel_in,
					debugLevel,
					dbg_X,
					dbg_Y);
			int [][] clusters_grown = enumerateClusters(
					grown_sel, //final boolean [][] selection,
					debugLevel,
					dbg_X,
					dbg_Y);
			int [][] grown_sub_clusters = matchPureGrown (
					clusters_grown,
					clusters_pure,
					debugLevel,
					dbg_X,
					dbg_Y);
			
			int [][] merged_no_conflict = mergeNoConflict(
					grown_sub_clusters,
					clusters_grown,
					clusters_pure,
					debugLevel,
					dbg_X,
					dbg_Y);

			setSurfaceData(
					merged_no_conflict, // final int [][] cluster_groups,
					clusters_pure, // final int [][] clusters_pure, // each tile/layer belongs to a single pure cluster
					0, // final int      debugLevel,
					dbg_X,
					dbg_Y);
			
			System.out.println("mergeAndGrow() done");
			

		}
		
		
}
