package com.elphel.imagej.tileprocessor;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;

/**
 ** TileNeibs - handles walking inside rectangular area
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TileNeibs.java is free software: you can redistribute it and/or modify
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

public class TileNeibs{
	final public static int DIR_N =       0; // UP
	final public static int DIR_NE =      1;
	final public static int DIR_E =       2; // Right
	final public static int DIR_SE =      3;
	final public static int DIR_S =       4; // Down
	final public static int DIR_SW =      5;
	final public static int DIR_W =       6; // Left
	final public static int DIR_NW =      7;
	final public static int DIR_CENTER = -1;
	final public static int DIR_UP =      0; // UP
	final public static int DIR_LEFT =    2; // Right
	final public static int DIR_DOWN =    4; // Down
	final public static int DIR_RIGHT =   6; // Left
	final public static int [][] DIR_XY = {{0,-1}, {1,-1}, {1,0}, {1,1}, {0,1}, {-1,1}, {-1,0}, {-1,-1}};
	final public static int DIRS =        DIR_XY.length; //8; // total dirs
	
	
	public static int reverseDir(int dir) {
		if ((dir < 0) || (dir >= DIRS)) {
			return dir;
		}
		return (dir+DIRS/2) % DIRS;
	}
	int sizeX;
	int sizeY;
	public int dirs = DIRS;
	public TileNeibs(int size){
		this.sizeX = size;
		this.sizeY = size;
	}
	public TileNeibs(int sizeX, int sizeY){
		this.sizeX = sizeX;
		this.sizeY = sizeY;
	}
	public int numNeibs() // TODO: make configurable to
	{
		return dirs;
	}
	public int opposite(int dir){
		return (dir + dirs / 2) % dirs;
	}

	int getLength(){
		return sizeX * sizeY;
	}

	boolean isBorder(int indx) {
		int [] xy = {indx % sizeX ,indx / sizeX};
		return (xy[0]==0) || (xy[1]==0) || (xy[0]== (sizeX - 1)) || (xy[1]==(sizeY-1)); 
	}
	
	
	/**
	 * Get x,y pair from index
	 * @param indx element index
	 * @return array of {x,y}
	 */

	int [] getXY(int indx)
	{
		int [] xy = {indx % sizeX ,indx / sizeX};
		return xy;
	}
	
	int getX(int indx) {return indx % sizeX;};
	
	int getY(int indx) {return indx / sizeX;};
	
	int getSizeX() {
		return sizeX;
	}
	int getSizeY() {
		return sizeY;
	}

	/**
	 * Get element index from x and y
	 * @param x horizontal position
	 * @param y vertical position
	 * @return element linescan index
	 */

	int getIndex(int x, int y){
		if ((x < 0) || (y < 0) || (x >= sizeX) || (y >= sizeY)) return -1;
		return y * sizeX + x;
	}

	int getIndex(int [] xy){
		if ((xy[0] < 0) || (xy[1] < 0) || (xy[0] >= sizeX) || (xy[1] >= sizeY)) return -1;
		return xy[1] * sizeX + xy[0];
	}

	/**
	 * Get element index offset by dx and dy from the indx
	 * @param indx start index
	 * @param dx offset in x direction
	 * @param dy offset in y direction
	 * @return new index or -1 if leaving
	 */
	int getNeibIndex(int indx, int dx, int dy) {
		int y = indx / sizeX + dy;
		int x = indx % sizeX + dx;
		if ((x < 0) || (y < 0 ) || (x >= sizeX) || (y >= sizeY)) return -1;
		return y * sizeX + x;
	}

	public boolean isInside (int indx, Rectangle roi) {
		if (indx < 0) return false;
		int y = indx / sizeX;
		int x = indx % sizeX;
		return (y >= roi.y) && (x >= roi.x)  || (y < roi.y + roi.height) || (x < roi.x + roi.width);
	}

	/**
	 * Get 2d element index after step N, NE, ... NW. Returns -1 if leaving array
	 * @param indx start index
	 * @param dir step direction (CW from up)
	 * @return new index or -1 if leaving
	 */
	public int getNeibIndex(int indx, int dir)
	{
		int y = indx / sizeX;
		int x = indx % sizeX;
		if (dir < 0) return indx;
		if (dir > 8) {
			System.out.println("getNeibIndex(): indx="+indx+", dir="+dir);
		}
		switch (dir){
		case 0: return (y == 0) ?                                    -1 : (indx - sizeX);
		case 1: return ((y == 0)           || ( x == (sizeX - 1))) ? -1 : (indx - sizeX + 1);
		case 2: return (                      ( x == (sizeX - 1))) ? -1 : (indx        + 1);
		case 3: return ((y == (sizeY - 1)) || ( x == (sizeX - 1))) ? -1 : (indx + sizeX + 1);
		case 4: return ((y == (sizeY - 1))                       ) ? -1 : (indx + sizeX);
		case 5: return ((y == (sizeY - 1)) || ( x == 0))           ? -1 : (indx + sizeX - 1);
		case 6: return (                      ( x == 0))           ? -1 : (indx        - 1);
		case 7: return ((y == 0)           || ( x == 0))           ? -1 : (indx - sizeX - 1);
		default: return indx;
		}
	}

	public static int getDX(int dir)
	{
		if (dir > 8) {
			System.out.println("getDX(): dir="+dir);
		}
		switch (dir){
		case 1:
		case 2:
		case 3:
			return 1;
		case 5:
		case 6:
		case 7:
			return -1;
		default:
			return 0;
		}
	}
	
	public static int getDY(int dir)
	{
		if (dir > 8) {
			System.out.println("getDY(): dir="+dir);
		}
		switch (dir){
		case 7:
		case 0:
		case 1:
			return -1;
		case 3:
		case 4:
		case 5:
			return 1;
		default:
			return 0;
		}
	}
	
	

	public static int getNumDirs(int radius) {
		if (radius < 0) {
			return 0;
		} else if (radius == 0) {
			return 1;
		} else {
			return 8 * radius;
		}
	}
	
	/**
	 * Get 2d element index after step of variable radius:
	 * radius==1 - same as getNeibIndex(int indx, int dir), 8 directions
	 * radius==2 - 16 directions (5x5 square), 0 - still up, north
	 * radius==3 - 24 directions (7x7 square)
	 * ...
	 * @param indx start index
	 * @param dir step direction (CW from up)
	 * @param radius - "distance" from the start point
	 * @return new index or -1 if leaving array in any direction
	 */

	public int getNeibIndexRadius(int indx, int dir, int radius) {
		if (radius < 2) {
			return  getNeibIndex(indx, dir);
		}
		int y = indx / sizeX;
		int x = indx % sizeX;
		if (dir > (8 * radius)) {
			System.out.println("getNeibIndex(): indx="+indx+", dir="+dir+", radius="+radius);
		}
		int dr = (dir + radius) % (8 * radius);
		int quad = dr / (2 * radius);
		int side = dr %  (2 * radius);
		switch (quad) {
		case 0:
			x = x - radius + side;
			y = y - radius;
			break;
		case 1:
			x = x + radius; 
			y = y - radius + side;
			break;
		case 2:
			x = x + radius - side; 
			y = y + radius;
			break;
		case 3:
			x = x - radius; 
			y = y + radius - side;
			break;
		}
		if ((x >= 0) && (y >= 0) && (x < sizeX) && (y < sizeY)) {
			return x + sizeX*y;
		} else {
			return -1;
		}
	}
	
	
	
	
	/**
	 * Get 2d element index after step N, NE, ... NW. Returns -1 if leaving array
	 * And 2 steps for dir = 8(N), 9(NNE),..23(NNW)
	 * @param indx start index
	 * @param dir step direction (CW from up)
	 * @return new index or -1 if leaving
	 */
	int getNeibIndex2(int indx, int dir)
	{
		int y = indx / sizeX;
		int x = indx % sizeX;
		if (dir < 0) return indx;
		if (dir > 24) {
			System.out.println("getNeibIndex(): indx="+indx+", dir="+dir);
		}
//		switch (dir % dirs){
		switch (dir){
		case  0: return (y == 0) ?                                    -1 : (indx - sizeX);
		case  1: return ((y == 0)           || ( x == (sizeX - 1))) ? -1 : (indx - sizeX + 1);
		case  2: return (                      ( x == (sizeX - 1))) ? -1 : (indx        + 1);
		case  3: return ((y == (sizeY - 1)) || ( x == (sizeX - 1))) ? -1 : (indx + sizeX + 1);
		case  4: return ((y == (sizeY - 1))                       ) ? -1 : (indx + sizeX);
		case  5: return ((y == (sizeY - 1)) || ( x == 0))           ? -1 : (indx + sizeX - 1);
		case  6: return (                      ( x == 0))           ? -1 : (indx        - 1);
		case  7: return ((y == 0)           || ( x == 0))           ? -1 : (indx - sizeX - 1);

		case  8: return ( y < 2) ?                                   -1 : (indx - 2 * sizeX);
		case  9: return ((y < 2)            || ( x > (sizeX - 2))) ? -1 : (indx - 2 * sizeX + 1);
		case 10: return ((y < 2)            || ( x > (sizeX - 3))) ? -1 : (indx - 2 * sizeX + 2);
		case 11: return ((y < 1)            || ( x > (sizeX - 3))) ? -1 : (indx - 1 * sizeX + 2);
		case 12: return (                      ( x > (sizeX - 3))) ? -1 : (indx             + 2);
		case 13: return ((y > (sizeY - 2))  || ( x > (sizeX - 3))) ? -1 : (indx + 1 * sizeX + 2);
		case 14: return ((y > (sizeY - 3))  || ( x > (sizeX - 3))) ? -1 : (indx + 2 * sizeX + 2);
		case 15: return ((y > (sizeY - 3))  || ( x > (sizeX - 2))) ? -1 : (indx + 2 * sizeX + 1);
		case 16: return ((y > (sizeY - 3))                       ) ? -1 : (indx + 2 * sizeX);
		case 17: return ((y > (sizeY - 3))  || ( x < 1))           ? -1 : (indx + 2 * sizeX - 1);
		case 18: return ((y > (sizeY - 3))  || ( x < 2))           ? -1 : (indx + 2 * sizeX - 2);
		case 19: return ((y > (sizeY - 2))  || ( x < 2))           ? -1 : (indx + 1 * sizeX - 2);
		case 20: return (                      ( x < 2))           ? -1 : (indx             - 2);
		case 21: return ((y < 1)            || ( x < 2))           ? -1 : (indx - 1 * sizeX - 2);
		case 22: return ((y < 2)            || ( x < 2))           ? -1 : (indx - 2 * sizeX - 2);
		case 23: return ((y < 2)            || ( x < 1))           ? -1 : (indx - 2 * sizeX - 1);
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
		int s1x = sizeX / 4;
		int s3x = 3 * sizeX / 4;
		int s1y = sizeY / 4;
		int s3y = 3 * sizeY / 4;
		int x = indx % sizeX;
		int y = indx / sizeX;
		boolean up = y < s1y;
		boolean down = y >= s3y;
		boolean left = x < s1x;
		boolean right = x >= s3x;
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
		return -1; // should not happen
	}
	public void shrinkSelection(
			int        shrink,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			boolean [] tiles,
			boolean [] prohibit)
	{
		boolean [] itiles = new boolean [tiles.length];
		for (int i = 0; i < tiles.length; i++) itiles[i] = !tiles[i];
		growSelection(
				shrink,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				itiles,
				prohibit);
		for (int i = 0; i < tiles.length; i++) tiles[i] = !itiles[i];
	}

	
	public boolean [] getEdgeSelection(
			int        shrink,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			boolean [] tiles,
			boolean [] prohibit)
	{
		boolean [] etiles = new boolean [tiles.length];
		for (int i = 0; i < tiles.length; i++) etiles[i] = !tiles[i];
		growSelection(
				shrink,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				etiles,
				prohibit);
		for (int i = 0; i < tiles.length; i++) etiles[i] &= tiles[i];
		return etiles;
	}
	
	
	

	public void growSelection(
			int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			boolean [] tiles,
			boolean [] prohibit)
	{
		boolean [] src_tiles = tiles.clone(); // just in case
		// grow
		boolean hor = true;
		int num_prev = 0;
		for (; grow > 0; grow--){
			boolean single = (grow ==1) && hor;
			src_tiles = tiles.clone();
			int num_new = 1; // as if previous pass was successful
			if (hor){
				for (int tileY = 0; tileY < sizeY; tileY++){
					for (int tileX = 0; tileX < (sizeX - 1); tileX++){
						int tindx = tileY * sizeX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx + 1])) {
							if (!tiles[tindx + 1] && src_tiles[tindx]){
								num_new++;
							}
							tiles[tindx + 1] |= src_tiles[tindx];
						}

					}
					for (int tileX = 1; tileX < sizeX; tileX++){
						int tindx = tileY * sizeX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx - 1])) {
							if (!tiles[tindx - 1] && src_tiles[tindx]){
								num_new++;
							}
							tiles[tindx - 1] |= src_tiles[tindx];
						}
					}
				}
			}
			if (!hor || single){ // do vertically, but from previous state
				for (int tileX = 0; tileX < sizeX; tileX++){
					for (int tileY = 0; tileY < (sizeY - 1); tileY++){
						int tindx = tileY * sizeX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx + sizeX])) {
							if (!tiles[tindx + sizeX] && src_tiles[tindx]){
								num_new++;
							}
							tiles[tindx + sizeX] |= src_tiles[tindx];
						}

					}
					for (int tileY = 1; tileY < sizeY; tileY++){
						int tindx = tileY * sizeX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx - sizeX])) {
							if (!tiles[tindx - sizeX] && src_tiles[tindx]){
								num_new++;
							}
							tiles[tindx - sizeX] |= src_tiles[tindx];
						}
					}
				}
			}
			hor = !hor;
			if ((num_new == 0) && (num_prev == 0)){
				break;
			}
			num_prev = num_new;
		}
	}

	public boolean [] boundShape(
			boolean [] selection,
			boolean octo)
	{
		boolean [] bound_shape = new boolean [selection.length];
		int min_x=-1, max_x=-1, min_y=-1, max_y=-1;
		int min_s=-1, max_s=-1, min_d=-1, max_d=-1;
		boolean is_set = false;
		for (int i = 0; i < selection.length; i++) if (selection[i]){
			int [] xy = getXY(i);
			int [] sd = {xy[0]+xy[1],xy[0]-xy[1]};
			if (!is_set) {
				min_x = xy[0];
				max_x = xy[0];
				min_y = xy[1];
				max_y = xy[1];
				if (octo) {
				min_s = sd[0];
				max_s = sd[0];
				min_d = sd[1];
				max_d = sd[1];
				}
				is_set = true;
			} else {
				if      (xy[0] < min_x) min_x = xy[0];
				else if (xy[0] > max_x) max_x = xy[0];
				if      (xy[1] < min_y) min_y = xy[1];
				else if (xy[1] > max_y) max_y = xy[1];
				if (octo) {
					if      (sd[0] < min_s) min_s = sd[0];
					else if (sd[0] > max_s) max_s = sd[0];
					if      (sd[1] < min_d) min_d = sd[1];
					else if (sd[1] > max_d) max_d = sd[1];
				}

			}
		}
		for (int y = min_y; y <= max_y; y++){
			for (int x = min_x; x <= max_x; x++){
				if (!octo ||
				(((x + y) >= min_s) && ((x + y) <= max_s) && ((x - y) >= min_d) && ((x - y) <= max_d))) {
					bound_shape[getIndex(x,y)] = true;

				}
			}
		}
		return bound_shape;
	}

	public static double[] getDoubleWindow(
			Rectangle window,
			double [] data,
			int data_width) {
		double [] window_data = new double [window.width * window.height];
		for (int row = 0; row < window.height; row++) {
			System.arraycopy(
					data,
					(window.y + row) * data_width + window.x,
					window_data,
					row * window.width,
					window.width);
		}
		return window_data;
	}
	
	public static void setDoubleWindow(
			Rectangle window,
			double [] window_data,
			double [] data,
			int       data_width) {
		for (int row = 0; row < window.height; row++) {
			System.arraycopy(
					window_data,
					row * window.width,
					data,
					(window.y + row) * data_width + window.x,
					window.width);
		}
	}
	

	public int [] distanceFromEdge(
			boolean [] tiles) {
		int [] dfe = new int [tiles.length];
		for (int i = 0; i < tiles.length; i++) dfe[i] = tiles[i] ? -1 : 0;
		ArrayList<Integer> front = new ArrayList<Integer>();
		int dist = 0;
		for (boolean has_empty = true; has_empty;) {
			has_empty = false;
			for (int start_indx = 0; start_indx < tiles.length; start_indx++) if (dfe[start_indx] < 0){
				has_empty = true;
				// see if it has wave front around
				boolean has_front_near = false;
				for (int d = 0; (d < dirs) && !has_front_near; d++){
					int ipx1 = 	getNeibIndex(start_indx, d);
					if ((ipx1 >= 0) && (dfe[ipx1] == dist)) {
						has_front_near = true;
					}
				}
				if (has_front_near) {
					dfe[start_indx] = dist+1;
					front.add(start_indx);
					// build wave front (may be several)
					while (!front.isEmpty()) {
						int ipx = front.remove(0);// get oldest element
						for (int d = 0; (d < dirs) && !has_front_near; d++){
							int ipx1 = 	getNeibIndex(ipx, d);
							if ((ipx1 >= 0) && (dfe[ipx1] < 0)) { // fresh cell
								// does it has front neighbor?

								for (int d1 = 0; (d1 < dirs) && !has_front_near; d1++){
									int ipx2 = 	getNeibIndex(ipx1, d1);
									if ((ipx2 >= 0) && (dfe[ipx2] == dist)) { // old front cell
										dfe[ipx1] = dist+1;
										front.add(ipx1);
										break;
									}
								}
							}
						}
					}		

				}
			}
			dist++; 
		}
		return dfe;
	}
	
	/**
	 * Enumerate clusters on rectangular area
	 * @param tiles        selected tiles, size should be sizeX * sizeY
	 * @param num_clusters if non null, will return number of clusters
	 * @param ordered      if true, order tiles from largest to smallest5
	 * @return integer array, where 0 is unused, 1+ cluster it belongs to
	 */
	public int [] enumerateClusters(
			boolean [] tiles,
			int []     num_clusters,
			boolean ordered)
	{
		int [] waves = new int [tiles.length];
		for (int i = 0; i < tiles.length; i++) waves[i] = tiles[i] ? 0: -1;
		ArrayList<Integer> front = new ArrayList<Integer>();
		int [] enum_clust = new int[tiles.length];
		int numClust = 0;
		for (int start_indx = 0; start_indx < tiles.length; start_indx++) if (waves[start_indx] == 0){ // found first pixel of a new cluster
			numClust ++;
			Integer ipx = start_indx;
			Integer ipx1;
			front.clear();
			int area = 1;
			waves[ipx] = area;
			enum_clust[ipx] = numClust;
			front.add(ipx);
			while (!front.isEmpty()) {
				ipx = front.remove(0);// get oldest element
				for (int d = 0; d < dirs; d++){
					ipx1 = 	getNeibIndex(ipx, d);
					if (ipx1 >= 0){
						if (waves[ipx1] == 0) {
							area++;
							waves[ipx1] = area;
							enum_clust[ipx1] = numClust;
							front.add(ipx1);
						}
					}
				}
			}
		}
		if (num_clusters != null) {
			num_clusters[0] = numClust;
		}
		if (!ordered) {
			return enum_clust;
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
		int [] enum_clust_ordered = new int[tiles.length];
		for (int i=0; i < enum_clust_ordered.length; i++){
			enum_clust_ordered[i] = (enum_clust[i] > 0) ? revIndex[enum_clust[i] - 1] : 0;
		}
		return enum_clust_ordered;
	}
	public static int getMax(
			int [] data)
	{
		int mx = data[0];
		for (int i = 1; i < data.length; i++){
			if (data[i] > mx){
				mx= data[i];
			}
		}
		return mx;
	}

	public int removeFewNeibs(
			boolean [] selection, // should be the same size
			int min_neibs)
	{
		int num_total_removed=0;
		int l = getLength();
		while (true) {
			int num_removed = 0;
			boolean [] to_remove = new boolean[l];
			for (int indx = 0; indx < l; indx++) if (selection[indx]){
				int num_neibs = 0;
				label_neibs: {
					for (int dir = 0; dir < 8; dir++) {
						int indx1 = getNeibIndex(indx, dir);
						if ((indx1 >= 0) && selection[indx1]) {
							num_neibs++;
							if (num_neibs >= min_neibs) {
								break label_neibs;
							}

						}
					}
					to_remove[indx] = true;
					num_removed++;
				}
			}
			if (num_removed > 0) {
				for (int i = 0; i < l; i++) if (to_remove[i]) {
					selection[i] = false;
				}
				num_total_removed += num_removed;
			} else {
				break;
			}
		}


		return num_total_removed;
	}
	
	public boolean [] fillConcave(
			boolean [] tiles // should be the same size
			) {
		boolean [] filled = tiles.clone();
		int min_neibs = 4;
		ArrayList<Integer> front = new ArrayList<Integer>();
		for (int start_indx = 0; start_indx < tiles.length; start_indx++) if (!filled[start_indx]) {
			int nneib = 0;
			for (int dir = 0; dir < dirs; dir++) {
				int id = getNeibIndex(start_indx, dir);
				if ((id >= 0) && filled[id]) nneib++;
			}
			if (nneib >= min_neibs) {
				front.add(start_indx); // do not mark yet
				filled[start_indx] = true;
				while (!front.isEmpty()) {
					int seed = front.remove(0);// get oldest element
					for (int d = 0; d < dirs; d++){
						int indx = 	getNeibIndex(seed, d);
						if ((indx >=0) && !filled[indx]) {
							nneib = 0;
							for (int dir = 0; dir < dirs; dir++) {
								int id = getNeibIndex(indx, dir);
								if ((id >= 0) && filled[id]) nneib++;
							}
							if (nneib >= min_neibs) {
								front.add(indx); // do not mark yet
								filled[indx] = true;
							}
						}
					}
				}				
			}
		}
		return filled;
	}
}
