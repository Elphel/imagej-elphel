import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;

/**
 **
 ** TileSurface - hadle tile surfaces
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

public class TileSurface {
//	public 
		private int tileSize;
		private int superTileSize;
		private int tilesX;
		private int tilesY;
		private int stilesX;
		private int stilesY;
		private int [] st_dirs8;
		private int [] t_dirs8;
		private int [] ss_dirs8;
		private double [] window;
		private int threadsMax = 100;
		
//		private int nsTilesstSize =   0; // 8;
		GeometryCorrection   geometryCorrection = null;
		public TileSurface(
				int tileSize,
				int superTileSize,
				int tilesX,
				int tilesY,
				GeometryCorrection geometryCorrection,
				int threadsMax){
			this.tileSize = tileSize;
			this.superTileSize = superTileSize;
			this.geometryCorrection =geometryCorrection;
			this.tilesX =  tilesX;
			this.tilesY =  tilesY;
//			int [] dirs =  {-tilesX, -tilesX + 1, 1, tilesX + 1, tilesX, tilesX - 1, -1, -tilesX - 1};
			int [] dirs =  {-stilesX, -stilesX + 1, 1, stilesX + 1, stilesX, stilesX - 1, -1, -stilesX - 1};
			this.st_dirs8 = dirs;

			int tx = superTileSize * stilesX; 
			int [] tdirs =  {-tx, -tx + 1, 1, tx + 1, tx, tx - 1, -1, -tx - 1};
			this.t_dirs8 = tdirs;
			
			
			int [] dirs_ss =  {-superTileSize, -superTileSize + 1, 1, superTileSize + 1, superTileSize, superTileSize - 1, -1, -superTileSize - 1};
			this.ss_dirs8 = dirs_ss;

			this.window = getWindow(2*superTileSize);
			this.threadsMax = threadsMax;
			stilesX = (tilesX + superTileSize -1)/superTileSize;
			stilesY = (tilesY + superTileSize -1)/superTileSize;
			
		}

		public class TileData{
			double [] disp_strength;
			int indx =         0;
			int new_index =    0;
			boolean enable =   true;
			int [] neighbors = null;
			public TileData (
					double disparity,
					double strength)
			{
				setDisparityStrength(disparity,strength);
			}

			public void setIndex(int indx)
			{
				this.indx = indx; 
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
			public void setNeighbor(int neib, int dir)
			{
				if (this.neighbors == null) this.neighbors = new int[8];
				this.neighbors[dir] = neib; 
			}
			public int getNeighbor(int dir)
			{
				if (this.neighbors == null) this.neighbors = new int[8];
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

		public class TileNeibs{
			int size;
			TileNeibs(int size){
				this.size = size;
			}
			/**
			 * Get 2d element index after step N, NE, ... NW. Returns -1 if leaving array   
			 * @param indx start index
			 * @param dir step direction (CW from up)
			 * @return new index or -1 if leaving 
			 */
			int getNeibIndex(int indx, int dir)
			{
				int y = indx / size;
				int x = indx % size;
				if (dir < 0) return indx;
				switch (dir % 8){
				case 0: return (y == 0) ?                               -1 : (indx - size); 
				case 1: return ((y == 0)         || ( x == (size-1))) ? -1 : (indx - size + 1); 
				case 2: return (                    ( x == (size-1))) ? -1 : (indx        + 1); 
				case 3: return ((y == (size -1)) || ( x == (size-1))) ? -1 : (indx + size + 1); 
				case 4: return ((y == (size -1))                    ) ? -1 : (indx + size); 
				case 5: return ((y == (size -1)) || ( x == 0))        ? -1 : (indx + size - 1); 
				case 6: return (                    ( x == 0))        ? -1 : (indx        - 1); 
				case 7: return ((y == 0)         || ( x == 0))        ? -1 : (indx - size - 1); 
				default: return indx;
				}
			}
			/**
			 * Return tile segment for 50% overlap. -1 - center, 0 N, 1 - NE,... 7 - NW
			 * @param indx element index
			 * @return which of the 9 areas this element belongs
			 */
			int getSegment(int indx)
			{
				int s1 = size / 4;
				int s2 = size /2;
				int s3 = 3 * size / 4;
				int x = indx % size;
				int y = indx / size;
				boolean up = y < s1;
				boolean down = y >= s3;
				boolean left = x < s1;
				boolean right = x >= s3;
				if (up){
					if (left) return 7;
					if (right) return 1;
					return 0;
				}
				if (down){
					if (left) return 5;
					if (right) return 3;
					return 4;
				}
				if (left) return 6;
				if (right) return 2;
				return -1;
			}

			/**
			 * Find if the step leaves the center half of all area
			 * @param indx start point
			 * @param dir direction
			 * @return direction to the new tile (assuming 50% overlap) or -1 if did not cross the border
			 */
			int leaveOvderlapedCenter(int indx, int dir)
			{
				int segm = getSegment(indx);
				int indx1 = getNeibIndex(indx, dir);
				if (indx1 < 0 ) {
					return -1; // should not happen
				}
				int segm1 = getSegment(indx1);
				if (segm1 == segm) return -1;
				if (segm == -1) return segm1;
				int [][] dxy = {
						{ 0, 0},
						{-1, 0},
						{-1, 1},
						{ 0, 1},
						{ 1, 1},
						{ 1, 0},
						{ 1,-1},
						{ 0,-1},
						{-1,-1}};
				int dx = dxy[segm1 + 1][1] - dxy[segm + 1][1]; 
				int dy = dxy[segm1 + 1][0] - dxy[segm + 1][0];
				for (int dp1 = 0; dp1 <=8; dp1++) {
					int sdx = (dx > 0) ? 1: ( (dx < 0) ? -1 : 0);
					int sdy = (dy > 0) ? 1: ( (dy < 0) ? -1 : 0);
					if ((dxy[dp1][0] == sdy) && (dxy[dp1][1] == sdx)){
						return dp1 -1;
					}
				}
				return -1; // shpuld not happen
			}
		}

		public int getNStileDir(
				int nsTile,
				int dir)
		{
			if (dir < 0) return nsTile;
			int sty = nsTile / stilesX;
			int stx = nsTile % stilesX;
			if ((stx > 0) && (sty > 0) && (sty == (stilesY - 1)) && (stx == (stilesX - 1))) return nsTile + st_dirs8[dir]; // most likely case
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
			if ((tx > 0) && (ty > 0) && (ty == (tilesY - 1)) && (tx == (tilesX - 1))) return nTile + t_dirs8[dir]; // most likely case
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
			int sdx = (dx > 0) ? 1: ( (dx < 0) ? -1 : 0);
			int sdy = (dy > 0) ? 1: ( (dy < 0) ? -1 : 0);
			if ((dy ==0 ) && (dx == 0)) return -1; // same tile
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
			if (dx < 0) return 2;
			if (dx > 0) return 6;
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
//			if (planes[nsTile] == null){
//				return -1; // empty supertile or supertile plane
//			}
			int tsn = (planes[nsTile] == null) ? 0 : planes[nsTile].length;
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
				int last_Length = (planes[nsTile1] == null) ? 0: planes[nsTile1].length;
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
				for (int j = 0; i < size; i++){
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
		 * Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid 
		 * @param planes array of the per-supertile, per plane plane data (each level can be null)
		 * @param debugLevel debug level
		 * @param dbg_X debug supertile X coordinate
		 * @param dbg_Y debug supertile Y coordinate
		 * @return per-tile (rounded up to contain whole supertiles) sparse array of TileData instances
		 */
		public TileData [][][] createTileShells0 (
				final boolean                   use_sel,
				final boolean                   divide_by_area,
				final double                    scale_projection,
				final TilePlanes.PlaneData [][] planes,
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final TileData [][][] tile_data = new TileData [nTiles][][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int len_st =  superTileSize * superTileSize;
//			final int ss2 = 2 * superTileSize;
//			final int sh =  superTileSize/2;
//			final int len2 = ss2  * ss2 ;
			final int [][] neib_of_neib_dir = {
				  //  N  NE   E  SE   S  SW   W  NW 
					{-1, -1,  1,  2, -1,  6,  7, -1}, // N  then E  = NE, N  then SE = E,  N then SW = W,  N then W = NW
					{-1, -1, -1, -1,  2, -1,  0, -1}, // NE then S  = E,  NE then W  = N
					{ 1, -1, -1, -1,  3,  4, -1,  0}, // E  then N  = NE, E  then S  = SE, E then SW = S,  E then NW = N
					{ 2, -1,  4, -1, -1, -1, -1, -1}, // SE then N  = E,  SE then W  = S
					{-1,  2,  3, -1, -1, -1,  5,  6}, // S  then NE = E,  S  then E  = SE, S then W  = SW, S then NW = W 
					{ 6, -1,  4, -1, -1, -1, -1, -1}, // SW then N  = W,  SW then E  = S  
					{ 7,  0, -1,  4,  5, -1, -1, -1}, // W  then N  = NW, W  then NE = N,  W then SE = S,  W then S  = SW 
					{-1, -1,  0, -1,  6, -1, -1, -1}, // NW then E  = N,  NW then S  = W
			};
			
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							double [][][][] disp_strengths = new double [9][][][];
							for (int dir = -1; dir < st_dirs8.length; dir++){
								int nsTile1 = getNStileDir(nsTile, dir); //   nsTile + ((dir < 0) ? 0: st_dirs8[dir]);
								if ((nsTile1 >= 0) && (planes[nsTile1] != null)){
									disp_strengths[dir] = new double [planes[nsTile1].length][][];
									for (int np = 0; np < planes[nsTile1].length; np++) {
										disp_strengths[dir + 1][np] = planes[nsTile1][np].getSinglePlaneDisparityStrength(
												getWindow(),      // double [] window,
												dir,              // int dir (-1 - center, 0- N, 1 - NE, .. 7 - NW
												use_sel,          // boolean   use_sel,
												divide_by_area,   //boolean   divide_by_area,
												scale_projection, // double    scale_projection,
												debugLevel-1);   // int       debugLevel)
									}
								}
							}
							// GET shifted/center value
							int num_surf = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
									nsTile, // int nsTile,
									8, // int dir,              // direction, or -1 (same)
									0, // int np,
									planes); // TilePlanes.PlaneData [][] planes) 
							if (num_surf > 0) {
//								int stileY = nsTile / stilesX;  
//								int stileX = nsTile % stilesX;
								tile_data[nsTile] = new TileData[superTileSize * superTileSize][num_surf];
								double [][][] all_disp_strengths = new double [num_surf][][];
								// First - process all surfaces for the existing supertile planes in this supertile
								if (planes[nsTile] != null) {
									for (int np = 0; np < planes[nsTile].length; np++) if (planes[nsTile][np] != null) {
										double [] strength =  disp_strengths[0][np][1].clone();
										double [] disparity = new double [len_st];
										for (int i = 0; i < len_st; i++){
											disparity[i] = disp_strengths[0][np][0][i] * disp_strengths[0][np][1][i]; 
										}
										for (int dir = 0; dir < st_dirs8.length; dir++){
											int sNeib = planes[nsTile][np].getNeibBest(dir);
											// add certain already shifted data from other planes around this one 
											if (sNeib >= 0){
												double [][] ds = disp_strengths[dir + 1][sNeib];
												for (int i = 0; i < len_st; i++){
													if (ds[1][i] > 0.0){
														strength[i] +=  ds[1][i]; 
														disparity[i] += ds[1][i] * ds[0][i]; 
													}
												}												
											}
											
										}
										int ns = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
												nsTile, // int nsTile,
												-1, // int dir,              // direction, or -1 (same)
												np, // int np,
												planes); // TilePlanes.PlaneData [][] planes) 
										all_disp_strengths[ns] = new double [2][];
										all_disp_strengths[ns][0] = disparity;
										all_disp_strengths[ns][1] = strength;
									}
								}
								// now process all tiles that are not connected to this one and let their planes "leak" here
								for (int dir = 0; dir < st_dirs8.length; dir++){
									int nsTile1 = getNStileDir(nsTile, dir); //   nsTile + ((dir < 0) ? 0: st_dirs8[dir]);
									if ((nsTile1 >= 0) && (planes[nsTile1] != null)){
										for (int np = 0; np < planes[nsTile1].length; np++) if (planes[nsTile1][np] != null) {
											// make sure it is not connected to the nsTile;
											int sNeib = planes[nsTile1][np].getNeibBest((dir + st_dirs8.length/2) % st_dirs8.length);
											if (sNeib < 0) { // surfaces connected to the current tile are already processed
												double [] strength =  disp_strengths[dir + 1][np][1].clone();
												double [] disparity = new double [len_st];
												for (int i = 0; i < len_st; i++){
													disparity[i] = disp_strengths[dir + 1][np][0][i] * disp_strengths[dir + 1][np][1][i]; 
												}
												for (int dir1 = 0; dir1 < st_dirs8.length; dir1++){
													int dir2 = neib_of_neib_dir[dir][dir1];
													if (dir2 >= 0) {
														sNeib = planes[nsTile1][np].getNeibBest(dir1);
														// add certain already shifted data from other planes around this one 
														if (sNeib >= 0){
															double [][] ds = disp_strengths[dir2 + 1][sNeib];
															for (int i = 0; i < len_st; i++){
																if (ds[1][i] > 0.0){
																	strength[i] +=  ds[1][i]; 
																	disparity[i] += ds[1][i] * ds[0][i]; 
																}
															}												
														}
													}
												}
												int ns = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
														nsTile, // int nsTile,
														dir, // int dir,              // direction, or -1 (same)
														np, // int np,
														planes); // TilePlanes.PlaneData [][] planes) 
												all_disp_strengths[ns] = new double [2][];
												all_disp_strengths[ns][0] = disparity;
												all_disp_strengths[ns][1] = strength;
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


		/**
		 * Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid 
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
				final TilePlanes.PlaneData [][] planes,
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{
			final int nStiles = stilesX * stilesY; 
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final double [][][][] fused_data = new double [nTiles][][][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								double [][][] disp_strength = new double [planes[nsTile].length][][];
								for (int np = 0; np < disp_strength.length; np++){
									if (planes[nsTile][np] != null){
										disp_strength[np] = planes[nsTile][np].getSinglePlaneDisparityStrength(
												getWindow(),      // double [] window,
												-1,              // int dir (-1 - center, 0- N, 1 - NE, .. 7 - NW
												use_sel,          // boolean   use_sel,
												divide_by_area,   //boolean   divide_by_area,
												scale_projection, // double    scale_projection,
												debugLevel-1);   // int       debugLevel)
										// multiply disparities by strengths to calculate weighted averages
										for (int i = 0; i < disp_strength[np][1].length; i++){
											disp_strength[np][0][i] *= disp_strength[np][1][i]; 
										}
									}
									for (int dir = 0; dir < st_dirs8.length; dir++){
										int sNeib = planes[nsTile][np].getNeibBest(dir);
										if (sNeib >= 0){
											int nsTile1 = getNStileDir(nsTile, dir); //   nsTile + ((dir < 0) ? 0: st_dirs8[dir]);
											if ((nsTile1 >= 0) && (planes[nsTile1] != null)){
												double [][] ds = planes[nsTile1][np].getSinglePlaneDisparityStrength(
														getWindow(),      // double [] window,
														dir,              // int dir (-1 - center, 0- N, 1 - NE, .. 7 - NW
														use_sel,          // boolean   use_sel,
														divide_by_area,   //boolean   divide_by_area,
														scale_projection, // double    scale_projection,
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
										}
									}
								}
								fused_data[nsTile] = disp_strength;
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
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int [][] dir_corn = {
					{ 7,  0,  6, -1},  // 0 (top left)
					{ 0,  1, -1,  2},  // 1 (top right)
					{ 6, -1,  5,  4},  // 2 (bottom left)
					{-1,  2,  4,  3}}; // 3 (bottom right)
			
			final int [][][][] corners = new int [nTiles][][][];
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								corners[nsTile] = new int [planes[nsTile].length][][];
								for (int np = 0; np < planes[nsTile].length; np++){
									if (planes[nsTile][np] != null){
										int [] neibs = planes[nsTile][np].getNeibBest();
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
			final int nTiles =  nStiles * superTileSize * superTileSize; 
			final int [][][][][] meshes = new int [nTiles][][][][];
			
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
						{ss2 * ss1 * 3 * sh, sh, ss1},
						{3 * sh * ss2 + ss1, ss1, sh}
					},
					{ // quadrant 3 - bottom right
						{3 * sh * ss2, ss1, sh},
						{ss1 * ss2,    sh, ss1}
					},
			};
			
			final TileNeibs tileNeibs = new TileNeibs(2*superTileSize);
			
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
							if (planes[nsTile] != null) {
								meshes[nsTile] = new int [planes[nsTile].length][][][];
								for (int np = 0; np < planes[nsTile].length; np++){
									if (planes[nsTile][np] != null){
										int [][] pre_mesh = new int [len_st2][2];
										for (int i = 0; i < len_st2; i ++){
											pre_mesh[i][0] = nsTile;
											pre_mesh[i][1] = np;
										}
										int [] neibs = planes[nsTile][np].getNeibBest();
										for (int quadrant = 0; quadrant <4; quadrant ++) {
											int [] these_corner_planes = corners[nsTile][np][quadrant];
											int [][] neib_id = new int[3][2];
											for (int arr = 0; arr < 3; arr++){
												int dir = quad_check[quadrant][arr][0];
												if (neibs[dir] >= 0) {
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
											// erase ortho
											for (int arr = 0; arr < 2; arr++) if (neib_id[arr] != null){
												for (int y = 0; y < cut_ortho[quadrant][arr][2]; y++){
													for (int x = 0; x < cut_ortho[quadrant][arr][1]; x++){
														int indx =  cut_ortho[quadrant][arr][0] + y * ss2 + x;
														pre_mesh[indx] = neib_id[arr]; 
//														meshes[nsTile][np][indx] = null;
													}
												}
											}
											// erase diagonal
											if (neib_id[2] != null){
												switch (quadrant){
												case 0: // top left 
													for (int j = 0; j < (ss1 - 1); j++){
														for (int i = ss1 - 1 - j; i>=0; i--){
//															meshes[nsTile][np][i * ss2 + j] = null;
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 1: // top right
													for (int j = ss1; j < ss2; j++){
														for (int i = j - ss1; i >= 0; i--){
//															meshes[nsTile][np][i * ss2 + j] = null;
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 2: // bottom left
													for (int j = 0; j < (ss1 - 1); j++){
														for (int i = ss1 + 1 + j; i < ss2; i++){
//															meshes[nsTile][np][i * ss2 + j] = null;
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												case 3: // bottom right
													for (int j = ss1; j < ss2; j++){
														for (int i = ss2 + sh - 1 - j; i < ss2; i++){
//															meshes[nsTile][np][i * ss2 + j] = null;
															pre_mesh[i * ss2 + j] = neib_id[2];
														}
													}
													break;
												}
											}
										}
										// build mesh , then add cuts if needed
										meshes[nsTile][np] = new int [len_st2][][];
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
											case 0:	dir_go = 0; dir_start =  6;	cut_right = true;  break;
											case 1:	dir_go = 0; dir_start = -1;	cut_right = false; break;
											case 2:	dir_go = 2; dir_start =  0;	cut_right = true;  break;
											case 3:	dir_go = 2; dir_start = -1;	cut_right = false; break;
											case 4:	dir_go = 4; dir_start = -1;	cut_right = true;  break;
											case 5:	dir_go = 4; dir_start =  6;	cut_right = false; break;
											case 6:	dir_go = 6; dir_start = -1;	cut_right = true;  break;
											case 7:	dir_go = 6; dir_start =  0;	cut_right = false; break;
											}
											int dir_go45 =   (dir_go + (cut_right ? 1:7)) % 8; 
											int dir_go90 =   (dir_go + (cut_right ? 2:6)) % 8; 
											int dir_go135 =  (dir_go + (cut_right ? 3:5)) % 8; 
											int dir_go180 =  (dir_go + 4) % 8; 

											indx = ss1 * (ss2 + 1); // center point
											
											for (int i = 0; i < sh; i++) indx = tileNeibs.getNeibIndex(indx, dir_go);
											if (dir_start >= 0) indx = tileNeibs.getNeibIndex(indx, dir_start);
											
											int indx1 = tileNeibs.getNeibIndex(indx, dir_go90);
											if ((pre_mesh[indx] != null) && (pre_mesh[indx1] == null)){ // there is a cut
												for (int i = 0; i < sh; i++){
													meshes[nsTile][np][tileNeibs.getNeibIndex(indx, dir_go180)][dir_go45] = null; // NE for N, right
													meshes[nsTile][np][indx][dir_go90] = null; // E for N, right
													if (i > 0){
														meshes[nsTile][np][indx][dir_go135] = null; // SE for N, right
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
			final int nTiles =  nStiles * superTileSize * superTileSize;
			final TileData [][] tile_data = new TileData [nTiles][];
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			final int ss2 = 2 * superTileSize;
			final int sh =  superTileSize/2;
			final int len2 = ss2  * ss2 ;
			final TileNeibs tileNeibs = new TileNeibs(2 * superTileSize);
			// initialize result structure
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
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
								int stileY = nsTile / stilesX;  
								int stileX = nsTile % stilesX;
								for (int np = 0; np < planes[nsTile].length; np++){
									int surf_number =  getTileSurfaceNumber ( // maximal number of surfaces in this supertile
											nsTile, // int nsTile,
											-1, // int dir,              // direction, or -1 (same)
											0, // int np,
											planes); // TilePlanes.PlaneData [][] planes)

									int [][][] src_mesh = lappingMeshes[nsTile][np];
									double [][] disp_strength = fusedSupertilePlanes[nsTile][np];
									TileData [] dual_mesh = new TileData [len2]; // full overlapping dual-sized mesh
									int [] sNeibs = planes[nsTile][np].getNeibBest();
									int [] surface_numbers = new int [8];
									for (int dir = 0; dir < 8; dir++){
										int nsTile1 = getNStileDir(nsTile, dir);
										if (sNeibs[dir] >= 0){
											surface_numbers[dir] = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
													nsTile1,  // int nsTile,
													-1, // int dir,              // direction, or -1 (same)
													sNeibs[dir], // int np,
													planes);
										} else if (nsTile1 >= 0){
											surface_numbers[dir] = getTileSurfaceNumber ( // maximal number of surfaces in this supertile
													nsTile1, // int nsTile,
													((dir + 4) % 4), // int dir,              // direction, or -1 (same)
													0, // int np,
													planes);
										} else { // out of the picture
											surface_numbers[dir] = -1; // out of the picture
										}
									}
									for (int indx = 0 ; indx < len2; indx++){
										if (src_mesh[indx] != null){
											int [][] src_neibs = src_mesh[indx];
											if (src_neibs != null){
												dual_mesh[indx] = new TileData(
														disp_strength[0][indx],  // disparity
														disp_strength[1][indx]); // strength

												int tsegm = tileNeibs.getSegment(indx);
												if (tsegm < 0) {
													dual_mesh[indx].setIndex(surf_number);
												} else {
													dual_mesh[indx].setIndex(surface_numbers[tsegm]);
												}

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
									// Now we have a double-sized surface with all tiles set with correct absolute indices, now just split it
									//surf_number =
									int sh3 = 3 * sh;
									for (int ty = 0; ty < ss2; ty++ ){
										for (int tx = 0; tx < ss2; tx++ ){
											int indx = ty * ss2 + tx;
											int tsegm = tileNeibs.getSegment(indx);
											int nsTile1 = getNStileDir(nsTile,tsegm);
											int ix, iy;
											switch (tsegm){
											case -1 : ix = tx - sh ; iy = ty - sh ; break;
											case  0 : ix = tx - sh ; iy = ty + sh ; break;
											case  1 : ix = tx - sh3; iy = ty + sh ; break;
											case  2 : ix = tx - sh3; iy = ty - sh ; break;
											case  3 : ix = tx - sh3; iy = ty + sh3; break;
											case  4 : ix = tx - sh ; iy = ty + sh3; break;
											case  5 : ix = tx + sh ; iy = ty + sh3; break;
											case  6 : ix = tx + sh ; iy = ty - sh ; break;
											case  7 : ix = tx + sh ; iy = ty + sh ; break;
											default:
												ix = tx -sh; iy = ty -sh;
											}
											if ((ix >= 0) && (ix < superTileSize) && (iy >= 0) && (iy < superTileSize)) {
												int tindx = ((stileY * superTileSize) + iy) * tilesX + ((stileX * superTileSize) + ix); 
												tile_data[tindx][dual_mesh[indx].getIndex()] = dual_mesh[indx];
//												tile_data[nsTile1][iy * superTileSize + ix][dual_mesh[indx].getIndex()] = dual_mesh[indx]; 
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
			final int ss2 = 2 * superTileSize;
			final int sh =  superTileSize/2;
			final int len2 = ss2  * ss2 ;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
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
									tdList.get(0).setNewIndex(i);
								}
								tile_data[nTile] = tdList.toArray(new TileData[0] );
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
						for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
							if (tile_data[nTile] != null){
								for (int i = 0; i < tile_data[nTile].length; i++){
									int [] neibs = tile_data[nTile][i].getNeighbors();
									for (int dir = 0; dir < neibs.length; dir++){
										if (neibs[dir] >= 0){
											int nTile1 = getNtileDir(nTile, dir);
											if (nTile1 >= 0){
												neibs[dir] = tile_data[nTile1][neibs[dir]].getNewIndex();
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
		
		public int getTileLayersNumber (
				final TileData [][] tileData)
		{
			int num = 0;
			for (int i = 0; i < tileData.length; i++){
				if (tileData[i].length > num ){
					num  = tileData[i].length; 
				}
			}
			return num;
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
			final double [][][] disp_strength = new double [numLayers][2][tileData.length];
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
		/**
		 * Calculate per-tile surface data (TileData) including disparity, strength, and 8 neighbors indices
		 * @param use_sel use plane selection (this.sel_mask) to select only some part of the plane
		 * @param divide_by_area divide weights by ellipsoid area
		 * @param scale_projection use plane ellipsoid projection for weight: 0 - do not use, > 0 linearly scale ellipsoid 
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
				final TilePlanes.PlaneData [][] planes,
				final int                       debugLevel,
				final int                       dbg_X,
				final int                       dbg_Y)
		{

			double [][][][] fused_planes = fuseSupertilePlanes (
					use_sel,           // final boolean                   use_sel,
					divide_by_area,    // final boolean                   divide_by_area,
					scale_projection,  // final double                    scale_projection,
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

			tileData = compactSortShells (
					tileData,          // final TileData [][]     tileData_src,
					debugLevel,        // final int                       debugLevel,
					dbg_X,             // final int                       dbg_X,
					dbg_Y);            // final int                       dbg_Y);
			return tileData;
		}
//getNtileDir		
}
