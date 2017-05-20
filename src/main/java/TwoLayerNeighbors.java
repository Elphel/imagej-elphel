import java.awt.Point;
import java.util.ArrayList;

/**
 **
 ** TwoLayerNeighbors - Handle connection swapping for resolving conflicts
 ** between two layer connections
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  TwoLayerNeighbors.java is free software: you can redistribute it and/or modify
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

public class TwoLayerNeighbors {
	/**
	 * connection variants (excluding from the center) as
	 * {start_dir, end_dir, relative_dir}
	 */
	static int [][] PAIRS = {
			{0, 1, 2},
			{1, 2, 4},
			{2, 3, 4},
			{3, 4, 6},
			{4, 5, 6},
			{5, 6, 0},
			{6, 7, 0},
			{7, 0, 2},
			{0, 2, 3},
			{2, 4, 5},
			{4, 6, 7},
			{6, 0, 1}};
	/**
	 *  Direction from dir1 (from center) to dir2 (from center)
	 *  -1 - same (dir1 == dir2)
	 *  -2 - impossible (no direct connection
	 */
	static int [][] rel_dirs = {
			//0   1   2   3   4   5   6   7	// dir2 values			
			{-1,  2,  3, -2, -2, -2,  5,  6}, // dir1 = 0;
			{ 6, -1,  4, -2, -2, -2, -2, -2}, // dir1 = 1;
			{ 7,  0, -1,  4,  5, -2, -2, -2}, // dir1 = 2;
			{-2, -2,  0, -1,  6, -2, -2, -2}, // dir1 = 3;
			{-2, -2,  1,  2, -1,  6,  7, -2}, // dir1 = 4;
			{-2, -2, -2, -2,  2, -1,  0, -2}, // dir1 = 5;
			{ 1, -2, -2, -2,  3,  4, -1,  0}, // dir1 = 6;
			{ 2, -2, -2, -2, -2, -2,  4, -1}};// dir1 = 7
	int nl1;
	int nl2;
	NeibVariant neibs_init = new NeibVariant();
	int [][] layers_around = new int [8][];
	int []   options_around = new int [8];
	int [][] num_se = new int [PAIRS.length][];
	int [][][] conns = new int [PAIRS.length][][];
	int [] selection_star = null;
	int [] selection_conns = null;
	
	class NeibVariant {
		
				
		int [][][] neighbors = new int[9][][];
		public int [][][] toArray()
		{
			return neighbors;
		}
		public void setNeighbors (int [][] neibs, int dir)
		{
			if (dir < 0) neighbors[8] =     neibs;
			else         neighbors[dir] =  neibs;
		}
		public int [][] getNeighbors (int dir)
		{
			if (dir < 0) return neighbors[8];
			else         return neighbors[dir];
		}
		
		public NeibVariant clone(){
			NeibVariant variant = new NeibVariant();
			variant.neighbors = neighbors.clone();
			for (int dir = 0; dir < neighbors.length; dir++){
				if (neighbors[dir] != null) {
					variant.neighbors[dir] = neighbors[dir].clone();
					for (int i = 0; i < neighbors[dir].length; i++){
						if (neighbors[dir][i] != null){
							variant.neighbors[dir][i] = neighbors[dir][i].clone();
						}
					}
				}
			}
			return variant; 
		}
		public int getDir2From1 (
				int dir1,
				int dir2)
		{
			if (dir1 < 0) return dir2;
			if (dir2 < 0) return ((dir1 + 4) % 8);
			return rel_dirs[dir1][dir2]; 
		}
		/**
		 * Connect tile at dir1 (-1 - center), layer nl1 to dir2, layer nl2
		 * Create connect in both directions, reconnect other ends of the broken links or plug with -1
		 * if there was none
		 * @param dir1 direction from the center to the start of the connection (-1 - center)
		 * @param nl1 start layer to connect (-1 - just disconnect the end)
		 * @param dir2 direction from the center to the end of the connection (-1 - center)
		 * @param nl2 end layer to connect (-1 - just disconnect the start)
		 */
		public void connect(
				int dir1,
				int nl1,
				int dir2,
				int nl2){
			int dir12 = getDir2From1(dir1, dir2);
			if (dir12 <0){
				throw new IllegalArgumentException ("Invalid connection from "+dir1+" to "+dir2+": resulted in direction 1->2 = "+dir12);
			}
			int dir21 = (dir12 + 4) % 8;
			int [][] neibs_start = getNeighbors(dir1);
			int [][] neibs_end =   getNeighbors(dir2);
			int old_nl2 = -1, old_nl1 = -1;
			if (nl1 >= 0){
				old_nl2 = neibs_start[nl1][dir12]; // where it was connected before, may be -1
				if (old_nl2 != nl2) {
					neibs_start[nl1][dir12] = nl2;
				}
			}
			if (nl2 >= 0){
				old_nl1 = neibs_end[nl2][dir21];
				neibs_end[nl2][dir21] = nl1;
			}
			// reconnect or plug broken links
			if (old_nl2 >= 0){
				neibs_end[old_nl2][dir21] = old_nl1; // (old_nl1 may be -1 here)
			}
			if (old_nl1 >= 0){
				neibs_start[old_nl1][dir12] = old_nl2; // (old_nl2 may be -1 here)
			}
			
		}
		
	}
	
	/**
	 * Initialize or advance variand selection. Return false if nothing left
	 * @return new selection available 
	 */
	public boolean nextSelection(){
		if (selection_star == null){
			selection_star = new int [options_around.length]; // 8
			selection_conns = new int [PAIRS.length]; // 12
			return true;
		} else {
			// increment connection variant if possible
			for (int np = 0; np < PAIRS.length; np++){
//				if (
//						(num_se[np] == null) || 
//						(conns[np] == null)) {
//					System.out.println("nextSelection BUG");
//					return false;
//				}
				if ((num_se[np] != null) && (num_se[np][0] == 2) && (num_se[np][1] == 2) && (conns[np] != null) && (conns[np].length == 1)){
					if (selection_conns[np] == 0){
						selection_conns[np] = 1;
						for (int i = 0; i < np; i ++){
							selection_conns[i] = 0;
						}
						return true;
					}
				}
			}
			// increment neighbor option, reset connection options;
			for (int i = 0; i < PAIRS.length; i ++){
				selection_conns[i] = 0;
			}
			for (int dir = 0; dir < options_around.length; dir++){
				if ((options_around[dir]>1) && (selection_star[dir] == 0)){
					selection_star[dir] = 1;
					for (int dir1 = 0; dir1 < dir; dir1++){
						selection_star[dir1] = 0;
					}
					return true;
				}
			}
			return false;
		}
	}
	
	/**
	 * Generate variant for the current selection (if consistent)
	 * @return neibVariant instance fro the current selection or null if the
	 * selection leads to conflicts  
	 */
	public NeibVariant generateVariant()
	{
		// verify all connetions are possible
		for (int np = 0; np < PAIRS.length; np++) if (conns[np] != null){
			// single connection for a single variant for start and end - either match or not
			if ((num_se[np] != null) && (conns[np].length == 1) && (num_se[np][0] == 1) && (num_se[np][1] == 1)){
				// Start and end of the connection belong to different groups - they can not be connected 
				if (selection_star[PAIRS[np][0]] != selection_star[PAIRS[np][1]]){
					return null;
				}
			}
		}
		
		// current selection is consistent, generate it
		NeibVariant variant = neibs_init.clone();
		
		// set connections for the center
		for (int dir = 0; dir < 8; dir++) if (options_around[dir] > 0){
			// make a first connection, if there are two - other will be created simultaneously
			variant.connect(
					-1, // 	int dir1,
					((selection_star[dir] > 0) ? nl2 : nl1), // int nl1,
					dir, // int dir2,
					layers_around[dir][0]); // int nl2);
		}
		
		// set all other connections
		for (int np = 0; np < PAIRS.length; np++) if (conns[np] != null){
			int start_dir = PAIRS[np][0];
			int end_dir =   PAIRS[np][1];
			boolean swap = (selection_star[start_dir] != selection_star[end_dir]) ^ (selection_conns[np] > 0);
			int [] opts = {0,0};
			if (swap){
				if (num_se[np][0] > 1){
					opts[0] = 1;
				} else {
					opts[1] = 1; // assuming there are two variants for the connection end as it should be
				}
			}
			variant.connect(
					start_dir, // 	int dir1,
					layers_around[start_dir][opts[0]], // int nl1,
					end_dir, // int dir2,
					layers_around[end_dir][opts[1]]); // int nl2);
			return variant;
		}
		return null;
	}
	
	public int [][][][] getNeighborVariants()
	{
		ArrayList<NeibVariant> variant_list = new ArrayList<NeibVariant>();
		while (nextSelection()){
			NeibVariant variant = generateVariant();
			if (variant != null){
				variant_list.add(variant);
			}
			
		}
		int [][][][] variants = new int [variant_list.size()][][][];
		int indx = 0;
		for(NeibVariant variant : variant_list){
			variants[indx++] = variant.toArray();
		}
		return variants;
	}
	
	
	public void setNeighbors (int [][] neibs, int dir)
	{
		neibs_init.setNeighbors(neibs, dir);
	}

	public int [][] getInitNeighbors (int dir)
	{
		return neibs_init.getNeighbors(dir);
	}

	public int [] getLayers()
	{
		int [] layers = {nl1,nl2};
		return layers; 
	}
	
	/**
	 * Return array pairs of start/end options or null if there are no connections for this pair
	 * @param np
	 * @return
	 */
	int [][] getConnections(int np)
	{
		int start_dir = PAIRS[np][0];
		int end_dir =   PAIRS[np][1];
		int dir =       PAIRS[np][2];
		if ((options_around[start_dir] > 0) && (options_around[end_dir] > 0)){
			ArrayList<Point> conn_list = new ArrayList<Point>();

			for (int opt1 = 0; opt1 < options_around[start_dir]; opt1++){
				int start_layer = layers_around[start_dir][opt1];
				for (int opt2 = 0; opt2 < options_around[end_dir]; opt2++){
					int end_layer = layers_around[end_dir][opt2];
					if (neibs_init.neighbors[start_dir][start_layer][dir] == end_layer){
						conn_list.add(new Point(opt1,opt2));
					}
				}
			}
			if (!conn_list.isEmpty()){
				int [][] pconns = new int [conn_list.size()][2];
				int indx = 0;
				for (Point p :conn_list){
					pconns[indx]  [0] = p.x;
					pconns[indx++][1] = p.y;
				}
				return pconns;
			}
		}
		return null;
	}
	
	public void setLayers(int nl1, int nl2)
	{
		this.nl1 = nl1;
		this.nl2 = nl1;
		layers_around = new int [8][];
		options_around = new int [8];
		int [][] neighbors_center = neibs_init.getNeighbors(-1);
		for (int dir = 0; dir < 8; dir++){
			if ((neighbors_center[nl1][dir] >= 0) || (neighbors_center[nl2][dir] >= 0)){
				options_around[dir] = 1;
				layers_around[dir] = new int [2];
				if (neighbors_center[nl1][dir] >= 0){
					layers_around[dir][0] = neighbors_center[nl1][dir];
					if (neighbors_center[nl2][dir] >= 0){
						layers_around[dir][1] = neighbors_center[nl2][dir];
						options_around[dir] ++;
					} else {
						layers_around[dir][1] = neighbors_center[nl1][dir];
					}
				} else {
					layers_around[dir][0] = neighbors_center[nl2][dir];
					layers_around[dir][1] = neighbors_center[nl2][dir];
				}
			}
		}
		for (int np = 0; np < PAIRS.length; np++){
			if ((options_around[PAIRS[np][0]] > 0) && (options_around[PAIRS[np][1]] > 0)){
				num_se[np] = new int [2]; // {options_around[PAIRS[np][0]],options_around[PAIRS[np][0]]};
				num_se[np][0] = options_around[PAIRS[np][0]];
				num_se[np][1] = options_around[PAIRS[np][1]];
				conns[np] = getConnections(np);
			}
		}
	}
	
}
