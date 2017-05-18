import java.awt.Point;
import java.util.HashMap;
import java.util.HashSet;

/**
 **
 ** ConnectionCosts - calculate and incrementally update cost of supertile connections
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  ConnectionCosts.java is free software: you can redistribute it and/or modify
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

public class ConnectionCosts {
	TilePlanes.PlaneData [][] planes =  null;
	boolean preferDisparity = false;
	TileSurface.TileNeibs tnSurface;
	double         orthoWeight;
	double         diagonalWeight;
	double         starPwr;    // Divide cost by number of connections to this power
	int            steps;
	int [][][]     neibs_init;
	int []         mod_tiles;
	int []         all_tiles;
	double [][][]  val_weights;
	double         init_val;
	double         init_weight;
	HashMap<Integer,Integer> tile_map; // map from tile full index to index in neibs[] and 

	public ConnectionCosts(
			double         orthoWeight,
			double         diagonalWeight,
			double         starPwr,    // Divide cost by number of connections to this power
			int            steps,
			TilePlanes.PlaneData [][] planes,
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity)
	{
		this.planes =          planes;
		this.preferDisparity = preferDisparity;
		this.tnSurface =       tnSurface;
		this.orthoWeight =     orthoWeight;
		this.diagonalWeight =  diagonalWeight;
		this.starPwr =         starPwr;         // Divide cost by number of connections to this power
		this.steps =           steps;
	}
	
	public int [][][] initConnectionCosts(
			int []         nsTiles)
	{
		int [] exp_tiles = nsTiles;
		for (int i = 1; i < steps; i++) exp_tiles = getInvolvedSupertiles(exp_tiles);
		mod_tiles = nsTiles;
		all_tiles = exp_tiles;
		val_weights = new double [all_tiles.length][][];
		neibs_init = new int [mod_tiles.length][][];
		for (int isTile = 0; isTile < mod_tiles.length; isTile++){
			int nsTile = mod_tiles[isTile];
			if (planes[nsTile] != null){
				neibs_init[isTile] = new int [planes[nsTile].length][];
				for (int nl = 0; nl < planes[nsTile].length; nl++) if ( planes[nsTile][nl] != null){
					neibs_init[isTile][nl] =  planes[nsTile][nl].getNeibBest();
				}
			}
		}

		for (int isTile = 0; isTile < all_tiles.length; isTile++){
			int nsTile = all_tiles[isTile];
			if (planes[nsTile] != null){
				val_weights[isTile] = new double [planes[nsTile].length][];
				for (int nl = 0; nl < planes[nsTile].length; nl++) if ( planes[nsTile][nl] != null){
					val_weights[isTile][nl] = new double[2];
				}
			}
		}

		tile_map = new HashMap<Integer,Integer>();
		for (int i = 0; i < mod_tiles.length; i++){
			tile_map.put(mod_tiles[i], i);
		}
		
		switch (steps){
		case 1:
			val_weights = getConnectionsCostSingleStep (
					null,	
					-1); // int        debugLevel)
			break;
		case 2:
			val_weights = getConnectionsCostDualStep (
					null,	
					-1); // int        debugLevel)
			break;
		default:
			val_weights = getConnectionsCostSingleStep (
					null,	
					-1); // int        debugLevel)			
		}
		
		init_val = 0.0;
		init_weight = 0.0; // should not change during update
		for (int isTile = 0; isTile < all_tiles.length; isTile++){
			if (val_weights[isTile] != null){
				for (int nl = 0; nl < val_weights[isTile].length; nl++) if ( val_weights[isTile][nl] != null){
					init_val +=    val_weights[isTile][nl][0] * val_weights[isTile][nl][1];
					init_weight += val_weights[isTile][nl][1];
				}
			}
		}
		// Likely weight will never change except first run, but we will still normalize by weight
		if (init_weight != 0.0) init_val /= init_weight;
		
		return neibs_init; // neighbors to clone
	}
	
	public double [][][] getConnectionsCostSingleStep (
			int [][][] neibs,	
			int        debugLevel)
	{
		boolean force = (neibs == null);
		if (force) neibs = neibs_init;
		double [][][] vw = new double [all_tiles.length][][];
		// now re-calculate val_weights where neibs are different from neibs_prev
		for (int isTile = 0; isTile < all_tiles.length; isTile++){
			int nsTile = all_tiles[isTile];
			if (planes[nsTile] != null) {
				int ineib = (tile_map.containsKey(nsTile))? tile_map.get(nsTile) : -1; 
				if (neibs[isTile] != null){
					vw[isTile] = new double[ planes[nsTile].length][];
					for (int nl = 0; nl < planes[nsTile].length; nl++) if (planes[nsTile][nl] != null){
						if  ((ineib >=0) && (neibs[ineib][nl] != null)) { // here just ignoring all tiles outside core ones (should be none)
							boolean neibs_changed = false;
							if (force){
								neibs_changed = true;
							} else {
								for (int dir = 0; dir < 8; dir++) if (neibs[isTile][nl][dir] != neibs_init[isTile][nl][dir]){
									neibs_changed = true;
									break;
								}
							}
							if (neibs_changed){
								vw[isTile][nl] = getStarValueWeight(
										nsTile,
										nl,
										neibs[isTile][nl],
										orthoWeight,
										diagonalWeight,
										starPwr, // double         starPwr,    // Divide cost by number of connections to this power
										tnSurface,
										preferDisparity,
										-1); // debugLevel);
							} else {
								vw[isTile][nl] = val_weights[isTile][nl];
							}
						} else {
							vw[isTile][nl] = null;
						}
					}				
				} else {
					vw[isTile] = null;
				}
			}
		}
		return vw; // negative - improvement
	}

	public double [][][] getConnectionsCostDualStep (
			int [][][] neibs,	
			int        debugLevel)
	{
		boolean force = (neibs == null);
		if (force) neibs = neibs_init;
		double [][][] vw = new double [all_tiles.length][][];
		// now re-calculate val_weights where neibs are different from neibs_prev
		for (int isTile = 0; isTile < all_tiles.length; isTile++){
			int nsTile = all_tiles[isTile];
			if (planes[nsTile] != null) {
				int ineib = (tile_map.containsKey(nsTile))? tile_map.get(nsTile) : -1; 
				if ((planes[nsTile] != null ) && ((ineib < 0) || (neibs[ineib] != null))){
					vw[isTile] = new double[ planes[nsTile].length][];
					for (int nl = 0; nl < planes[nsTile].length; nl++) if (planes[nsTile][nl] != null){
						int [] neibs0 = (ineib >= 0) ? neibs[ineib][nl] : planes[nsTile][nl].getNeibBest();
						if  (neibs0 != null) {
							boolean neibs_changed = false;
							if (force){
								neibs_changed = true;
							}
							int [][] neibs2 = new int [8][];
							for (int dir = 0; dir < 8; dir++ ){
								if (!neibs_changed && (ineib >= 0) &&
										((neibs_init[ineib] == null) || (neibs_init[ineib][nl] == null) || (neibs0[dir] != neibs_init[ineib][nl][dir]))) {
									neibs_changed = true;
								}
								if (neibs0[dir] >= 0) {
									int nl1 = neibs0[dir];
									int nsTile1 = tnSurface.getNeibIndex(nsTile, dir);
									int ineib1 = (tile_map.containsKey(nsTile1))? tile_map.get(nsTile1) : -1;
									neibs2[dir] = (ineib1 >=0) ? neibs[tile_map.get(nsTile1)][nl1] : planes[nsTile1][nl1].getNeibBest();
									if (!neibs_changed && (ineib1 >= 0)) {
										if ((neibs_init[ineib1] == null) || (neibs_init[ineib1][nl1] == null)) {
											neibs_changed = true;
										} else {
											for (int dir1 = 0; dir1 < 8; dir1++) if (neibs[ineib1][nl1][dir1] != neibs_init[ineib1][nl1][dir1]){
												neibs_changed = true;
												break;
											}
										}
									}
								}
							}
							
							if (neibs_changed){
								vw[isTile][nl] = getStarValueWeight(
										nsTile,
										nl,
										neibs0,
										neibs2,
										orthoWeight,
										diagonalWeight,
										starPwr, // double         starPwr,    // Divide cost by number of connections to this power
										tnSurface,
										preferDisparity,
										-1); // debugLevel);
							} else {
								vw[isTile][nl] = val_weights[isTile][nl];
							}
						} else {
							vw[isTile][nl] = null;
						}
					}				
				} else {
					vw[isTile] = null;
				}
			}
		}
		return vw; // negative - improvement
	}
	
	
	
	public double getConnectionsCostDiff (
			int [][][]     neibs,		// should be initialized at top dimension if neibs_prev==null	
			int    debugLevel)
	{
		double [][][] vw;
		switch (steps){
		case 1:
			vw = getConnectionsCostSingleStep (
					neibs,	
					-1); // int        debugLevel)
			break;
		case 2:
			vw = getConnectionsCostDualStep (
					neibs,	
					-1); // int        debugLevel)
			break;
		default:
			vw = getConnectionsCostSingleStep (
					neibs,	
					-1); // int        debugLevel)			
		}

		// calculate new cost
		double new_value = 0.0;
		double new_weight = 0.0; // should not change during update
		for (int isTile = 0; isTile < all_tiles.length; isTile++){
			if (vw[isTile] != null){
				for (int nl = 0; nl < vw[isTile].length; nl++) if ( vw[isTile][nl] != null){
					new_value +=   vw[isTile][nl][0] * vw[isTile][nl][1];
					new_weight +=  vw[isTile][nl][1];
				}
			}
		}
		if (new_weight != 0.0) new_value /= new_weight;
		return new_value - init_val; // negative - improvement
	}
	
	/**
	 * Calculate main eigenvalue of the current plane and all connected ones - used to estimate advantage of connection swap
	 * @param nsTile supertile index
	 * @param nl surface layer
	 * @param neibs array of 8 neighbors layers (N,NE,...NW), -1 - not connected
	 * @param orthoWeight multiply contribution of ortho neighbors
	 * @param diagonalWeight  multiply contribution of diagonal neighbors
	 * @param diagonalWeight  divide value by number of connections to this power (if !=0)
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return a pair of eigenvalue of the combine plane and its weight
	 */
	public double [] getStarValueWeight(
			int    nsTile,
			int    nl,
			int [] neibs,
			double orthoWeight,
			double diagonalWeight,
			double starPwr,    // Divide cost by number of connections to this power
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity,
			int    debugLevel)
	{
		TilePlanes.PlaneData merged_plane = planes[nsTile][nl];
		for (int dir = 0; dir < 8; dir++){
			if (neibs[dir] >= 0){
				double other_weight = ((dir & 1) != 0) ? diagonalWeight : orthoWeight;
				TilePlanes.PlaneData other_plane = merged_plane.getPlaneToThis(  // layer here does not matter
						planes[tnSurface.getNeibIndex(nsTile, dir)][neibs[dir]],
						debugLevel - 1); // debugLevel);
				merged_plane = merged_plane.mergePlaneToThis(
						other_plane,     // PlaneData otherPd,
						other_weight,    // double    scale_other,
						false,           // boolean   ignore_weights,
						true,            // boolean   sum_weights,
						preferDisparity, 
						debugLevel - 1); // int       debugLevel)
			}
		}
		double [] value_weight = {merged_plane.getValue(),merged_plane.getWeight()};
		if (starPwr != 0){
			value_weight[0] /= (Math.pow((planes[nsTile][nl].getNumNeibBest() + 1.0), starPwr));
		}
		return value_weight;
	}

	/**
	 * Calculate main eigenvalue of the current plane and all connected ones - used to estimate advantage of connection swap
	 * This version uses two steps - not only directly connected, but neighbors' neighbors also, multiple paths to the same
	 * tile add together.
	 * @param nsTile supertile index
	 * @param nl surface layer
	 * @param neibs array of 8 neighbors layers (N,NE,...NW), -1 - not connected
	 * @param neibs2 2-d array of 8 neighbors' neighbors layers (N,NE,...NW), -1 - not connected
	 * @param orthoWeight multiply contribution of ortho neighbors
	 * @param diagonalWeight  multiply contribution of diagonal neighbors
	 * @param diagonalWeight  divide value by number of connections to this power (if !=0)
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @param preferDisparity - the first eigenvalue/vector is the most disparity-like
	 *                          (false - smallest eigenvalue)
	 * @param debugLevel
	 * @return a pair of eigenvalue of the combine plane and its weight
	 */
	
	
	public double [] getStarValueWeight(
			int    nsTile,
			int    nl,
			int [] neibs,
			int [][] neibs2, // neighbors' neighbors
			double orthoWeight,
			double diagonalWeight,
			double starPwr,    // Divide cost by number of connections to this power
			TileSurface.TileNeibs tnSurface,
			boolean preferDisparity,
			int    debugLevel)
	{
		double [] dir_weight = {orthoWeight,  diagonalWeight, orthoWeight,  diagonalWeight, orthoWeight,  diagonalWeight, orthoWeight,  diagonalWeight};
		HashMap<Point, Double> tile_weights = new HashMap<Point, Double>();
		for (int dir = 0; dir < 8; dir++){
			int nl1 =  neibs[dir];
			if (nl1 >= 0){
				int nsTile1 = tnSurface.getNeibIndex(nsTile, dir);
				double weight1 = dir_weight[dir];
				tile_weights.put(new Point(nsTile1, nl1), new Double(weight1)); // no need to check for existence here
				for (int dir1 = 0; dir1 < 8; dir1++){
					if ((dir1 != dir) && (neibs2[dir]!= null)){
						 int nl2 =neibs2[dir][dir1];
						 if (nl2 >= 0){
							 Point p = new Point (tnSurface.getNeibIndex(nsTile1, dir1), nl2);
							 Double w0 = tile_weights.get(p);
							 double weight2 = dir_weight[dir1]*weight1;
							 if (w0 != null) weight2 += w0;
							 tile_weights.put(p, new Double(weight2));
						 }
					}
				}
			}
		}
		TilePlanes.PlaneData merged_plane =  planes[nsTile][nl]; // center point
		for (HashMap.Entry<Point, Double> entry : tile_weights.entrySet()){
			TilePlanes.PlaneData other_plane = merged_plane.getPlaneToThis(  // layer here does not matter
					planes[entry.getKey().x][entry.getKey().y],
					debugLevel - 1); // debugLevel);
			merged_plane = merged_plane.mergePlaneToThis(
					other_plane,      // PlaneData otherPd,
					entry.getValue(), // double    scale_other,
					false,           // boolean   ignore_weights,
					true,            // boolean   sum_weights,
					preferDisparity, 
					debugLevel - 1); // int       debugLevel)
		}
		double [] value_weight = {merged_plane.getValue(),merged_plane.getWeight()};
		if (starPwr != 0){
			value_weight[0] /= (Math.pow(tile_weights.size() + 1.0, starPwr));
		}
		return value_weight;
	}

	/**
	 * Calculate array of supertile indices that need to have connection cost recalculated when they are updated
	 * first entries of the result will be the same in input array
	 * @param mod_supertiles array of supertile indices that will be modified
	 * @param tnSurface TileNeibs instance to navigate tile index and control array borders
	 * @return array of supertile indices to watch connection cost
	 */
	public int [] getInvolvedSupertiles(
			int [] mod_supertiles)
	{
		HashSet<Integer> stiles_set = new HashSet<Integer>();
		for (int i = 0; i < mod_supertiles.length; i++){
			stiles_set.add(new Integer(mod_supertiles[i]));
			for (int dir = 0; dir < 8; dir++){
				Integer nsTile = tnSurface.getNeibIndex(mod_supertiles[i], dir);
				if (nsTile >= 0) stiles_set.add (nsTile);
			}
		}
		int [] stiles = new int [stiles_set.size()];
		int indx = 0;
		for (; indx < mod_supertiles.length; indx++){
			stiles_set.remove(new Integer(mod_supertiles[indx]));
			stiles[indx] = mod_supertiles[indx];
		}
		for (Integer nsTile: stiles_set){
			stiles[indx++] = nsTile;
		}
		return stiles;
	}
	
	
}
