package com.elphel.imagej.tileprocessor;

import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.ShowDoubleFloatArrays;

/**
 ** TexturedModel - Generate 3D mesh with textures
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TexturedModel.java is free software: you can redistribute it and/or modify
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

public class TexturedModel {
	public static final int          threadsMax = 100;  // maximal number of threads to launch
	public static final int TILE_EMPTY =          0; 
	public static final int TILE_BORDER =         1; 
	public static final int TILE_BORDER_FLOAT =   2; 
	public static final int TILE_CONFIRMED =      3; 
	public static final int TILE_CANDIDATE =      4; // not used 
	public static final int CLUSTER_NAN =        -2; // disparity is NaN 
	public static final int CLUSTER_UNASSIGNED =  -1; // not yet assinged (>=0 - cluster number)
	
	public static boolean isBorder(int d) {
		return (d==TILE_BORDER) || (d==TILE_BORDER_FLOAT);
	}
	
	public static TileCluster [] clusterizeFgBg( //
			final int          tilesX,
			final double [][]  disparities, // may have more layers
			final boolean []   selected, // to remove sky (pre-filter by caller, like for ML?)
			final double       disp_adiffo,
			final double       disp_rdiffo,
			final double       disp_adiffd,
			final double       disp_rdiffd,
			final double       disp_fof,    // enable higher difference (scale) for fried of a friend 
			final int          debugLevel) {
		final double disp_border = disp_fof;
		final int tiles = disparities[0].length;
		final int tilesY = tiles/tilesX;
		final int layers = disparities.length;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
//		final boolean [][][] connections = new boolean [tiles][][];
		final double [][][][] connections = new double [tiles][][][];
//		final boolean [][][][] bconn = new double [boolean][][][];
//		final int [][] num_neibs =     new int[tiles][layers];
		final int [][] num_neibs_dir = new int[tiles][layers]; // -1 - none, otherwise - bitmask
		final int [][] ncluster =  new int [tiles][layers]; // -2 - NaN, -1 - not yet assigned, 0+ - cluster number
		for (int tile = 0; tile < tiles; tile++) {
			Arrays.fill(ncluster[tile], CLUSTER_NAN);
		}
		final int [] num_bits = new int [512];
		for (int i = 0; i < num_bits.length; i++) {
			for (int d = i; d != 0; d>>=1) {
				num_bits[i]+=(d & 1);
			}
		}

		final int dbg_tile = (debugLevel>0)? 2021:-1; // 977 : -1;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if ((selected == null) || selected[tile]) {
						if (tile==dbg_tile) {
							System.out.println("clusterizeFgBg().1: tile="+tile);
						}
						for (int layer = 0; layer < layers; layer++) {
							if (!Double.isNaN(disparities[layer][tile])) {
								ncluster[tile][layer] = CLUSTER_UNASSIGNED; // not yet assigned
								if (connections[tile] == null) {
									connections[tile] = new double[layers][][];
								}
//								connections[tile][layer] = new double [TileNeibs.DIRS/2][];
								connections[tile][layer] = new double [TileNeibs.DIRS][]; // leave room for future symmetry
								for (int dir = 0; dir < TileNeibs.DIRS/2; dir++) {
									int tile1 = tn.getNeibIndex(tile, dir);
									if (tile1 >= 0) {
										for (int layer1 = 0; layer1 < layers; layer1++) {
											if (!Double.isNaN(disparities[layer1][tile1])) {
												if (connections[tile][layer][dir] == null) {
													connections[tile][layer][dir] = new double[layers];
													Arrays.fill(connections[tile][layer][dir],Double.NaN);
												}
												
												double mid_disp = Math.max(0.0, 0.5*(disparities[layer][tile] + disparities[layer1][tile1]));
											    double max_disp_diff = ((dir & 1) == 0) ?
											    		(disp_adiffo + mid_disp * disp_rdiffo) :
											    			(disp_adiffd + mid_disp * disp_rdiffd);
											    connections[tile][layer][dir][layer1] = Math.abs(disparities[layer][tile] - disparities[layer1][tile1])/max_disp_diff;
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

		ai.set(0);
		// calculate total number of connections (w/o fof) by combining opposite directions
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (connections[tile] != null) {
						if (tile==dbg_tile) {
							System.out.println("clusterizeFgBg().2: tile="+tile);
						}
						for (int layer = 0; layer < layers; layer++) if (connections[tile][layer] != null){
							for (int dir0 = 0; dir0 < TileNeibs.DIRS/2; dir0++) {
								int dir = TileNeibs.reverseDir(dir0);
								int tile1 = tn.getNeibIndex(tile, dir);
								if ((tile1 >= 0) && (connections[tile1] != null)) {
									if (connections[tile][layer][dir] == null) {
										connections[tile][layer][dir] = new double[layers];
										Arrays.fill(connections[tile][layer][dir],Double.NaN);
									}
									for (int layer1 = 0; layer1 < layers; layer1++) {
										if (    (connections[tile1][layer1] != null) &&
												(connections[tile1][layer1][dir0] != null)) {
											connections[tile][layer][dir][layer1] = connections[tile1][layer1][dir0][layer];
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
		
		ai.set(0);
		// calculate total number of connections (w/o fof) by combining opposite directions
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (connections[tile] != null) {
						if (tile==dbg_tile) {
							System.out.println("clusterizeFgBg().3: tile="+tile);
						}
						for (int layer = 0; layer < layers; layer++) if (connections[tile][layer] != null){
//							num_neibs[tile][layer]++; // has center tile, maximal will be 9
							num_neibs_dir[tile][layer] = 1; // center tile
							for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
								if (connections[tile][layer][dir] != null) {
									for (int layer1 = 0; layer1 < layers; layer1++) {
										if (connections[tile][layer][dir][layer1] <= 1.0) { // Double.NaN - OK
//											num_neibs[tile][layer]++;
											num_neibs_dir[tile][layer] |= 2 << dir; // 9 bits
											break; // increment once per dir
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
		if (debugLevel > 0) {
			String [] dbg_titles = {"FG","BG"};
			double [][] dbg_img = new double[layers][tiles];
			for (int i = 0; i < tiles;i++) {
				for (int j = 0; j < dbg_img.length; j++) {
					dbg_img[j][i] = num_bits[num_neibs_dir[i][j]];
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"num_neibs",
					dbg_titles);
			(new ShowDoubleFloatArrays()).showArrays(
					disparities,
					tilesX,
					tilesY,
					true,
					"disparities",
					dbg_titles);
		}
		final ArrayList <TileCluster> cluster_list = new ArrayList<TileCluster>();
		// build all clusters
		int tile_start = 0;
		double disparity[] = new double[tiles]; // current cluster disparities
		int [] tile_stat = new int [tiles];
//		int [] tile_layer = new int [tiles]; // just to know which layer was used for assigned tiles
//		int current_cluster = 0;
		final boolean debug_index = debugLevel > 0;

		while (true) {
			// find remaining tile with maximal number of neighbors (usually 8) - this will require num_neibs update
			int best_tile = -1, best_layer = -1, nn=0;
			find_start:
			{
				int tile_end = tile_start + tiles;
				for (int tile1=tile_start; tile1 < tile_end; tile1++) {
					int tile = (tile1 >= tiles) ? (tile1 - tiles) : tile1;
					if (!tn.isBorder(tile)) { // do not start on the border
						for (int layer = 0; layer < layers; layer++) {
							int n_neibs = num_bits[num_neibs_dir[tile][layer]];
							if ((ncluster[tile][layer] == CLUSTER_UNASSIGNED) && (n_neibs > nn)) { // not yet assigned and >=0 neibs
								nn = n_neibs;
								best_tile = tile;
								best_layer = layer;
								if (nn == (TileNeibs.DIRS + 1)) {
									break find_start;
								}
							}
						}
					}
				}
			}
			if (best_tile == dbg_tile) {
				System.out.println("clusterizeFgBg().4: best_tile="+best_tile);
			}
			tile_start = best_tile; // will start from this
			if (nn == 0) { // no even single-tile clusters are left
				break; // no more seeds for clusters
			}
			Arrays.fill(tile_stat,0);
			ArrayList<Point> tile_layer_list = new ArrayList<Point>(); // pair of int x tile, int y - layer
			Point p0 = new Point(best_tile, best_layer); // First point will always be good
			int num_tiles = 0;
			tile_layer_list.add(p0);
			while (!tile_layer_list.isEmpty()) {
				Point p = tile_layer_list.remove(0);
				// p is a candidate, initially will always be approved
				int tile = p.x;
				int layer = p.y;
				// see if the current tile candidate is a good one
				boolean confirm = false;
				disparity[tile] = disparities[layer][tile];
				if (tn.isBorder(tile)) {
					tile_stat[tile] = TILE_BORDER_FLOAT; // TILE_BORDER;
					if (Double.isNaN(disparity[tile])) {
						System.out.println("clusterizeFgBg(): 4. Disparity is set NaN for tile "+tile);
					}
				} else {
					confirm = true;
					for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
						int tile1 = tn.getNeibIndex(tile, dir); // should always be > 0 here
						//is it a border tile or already confirmed one ?
//						if ((tile_stat[tile1] == TILE_BORDER) || (tile_stat[tile1] == TILE_CONFIRMED)){
						if (isBorder(tile_stat[tile1]) || (tile_stat[tile1] == TILE_CONFIRMED)){
							double mid_disp = Math.max(0.0, 0.5*(disparities[layer][tile] + disparity[tile1]));
							double max_disp_diff = ((dir & 1) == 0) ?
									(disp_adiffo + mid_disp * disp_rdiffo) :
										(disp_adiffd + mid_disp * disp_rdiffd);
//							max_disp_diff *= (tile_stat[tile1] == TILE_BORDER)? disp_border : disp_fof;
							max_disp_diff *= isBorder(tile_stat[tile1])? disp_border : disp_fof;
							if ((Math.abs(disparities[layer][tile] - disparity[tile1])/max_disp_diff) > 1.0){
								confirm = false;// too large diff
								// make it a border tile, but keep disparity
								//								disparity[tile] = disparities[layer][tile];
								tile_stat[tile] = TILE_BORDER_FLOAT; // TILE_BORDER;
								break;
							}
						}
					}
				}
				if (confirm) {
					if (debugLevel > 1) {
						System.out.println("Confirmed tile "+tile+ " (x="+(tile%tilesX)+", y="+(tile/tilesX)+") -> "+tile_stat[tile]);
					}
					tile_stat[tile] = TILE_CONFIRMED;
					if (Double.isNaN(disparity[tile])) {
						System.out.println("clusterizeFgBg(): 5. Disparity is set NaN for tile "+tile);
					}
					ncluster[tile][layer] = cluster_list.size(); // current cluster number - Mark as assigned
					//						tile_layer[tile] = layer;
					// right here - update number of neighbors for used tile
					num_neibs_dir[tile][layer] = 0; // bitmap, 0- -center, 1..8 - up, ...
					// update number of neighbors pointed here

					num_tiles++; // counting only core tiles, not the border one
					// look around (including layers - select best), add candidates (if fit) and mark borders (if does not fit)
					for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
						int tile1 = tn.getNeibIndex(tile, dir); // should always be > 0 here
						if (tile_stat[tile1] == TILE_EMPTY) { // Do not consider already assigned tiles for this cluster
							// update neighbors bits (for number of neighbors calculation)
							int rdir = TileNeibs.reverseDir(dir);
							num_neibs_dir[tile][layer] &= ~(2 << rdir); // zero out pointed to this tile
							for (int layer1 = 0; layer1 < layers; layer1++) if (ncluster[tile1][layer1] == CLUSTER_UNASSIGNED) { // from where connection
								// ncluster[tile][layer] is already assigned
								for (int layer2 = 0; layer2 < layers; layer2++) if (ncluster[tile][layer2] == CLUSTER_UNASSIGNED){ // to where connection
									if (connections[tile1][layer1][rdir][layer2] <= 1.0) { // Double.NaN - OK
										num_neibs_dir[tile][layer] |= 2 << rdir; // restore as there is another layer that fits
										break;
									}
								}								
							}							
							// does it have any layers in this direction?
							int layer1 = -1;
							int layer_assigned = -1; //earlier assigned, can not continue but use disparity
							if ((connections[tile][layer] != null) && (connections[tile][layer][dir] != null)) {
								for (int l = 0; l < layers; l++) if (connections[tile][layer][dir][l] <= 1.0){
									// verify it was not already used. If used - use used disparity for the border 									
									// **** (select correct layer for that border disparity
									if (ncluster[tile1][l] == CLUSTER_UNASSIGNED) {
										if ((layer1 < 0) || (connections[tile][layer][dir][l] < connections[tile][layer][dir][layer1])) {
											layer1 = l;
										}
									} else { // already assigned to some cluster that was split
										if ((layer_assigned < 0) || (connections[tile][layer][dir][l] < connections[tile][layer][dir][layer_assigned])) {
											layer_assigned = l;
										}
									}
								}
							}
							if (layer_assigned >= 0) { // already assigned to some cluster that was split
								tile_stat[tile1] = TILE_BORDER; // here - fixed disparity, not float
								disparity[tile1] = disparities[layer_assigned][tile1]; // [tile]; // Use interrupted cluster disparity
								if (Double.isNaN(disparity[tile1])) {
									System.out.println("clusterizeFgBg(): *1. Disparity is set NaN for tile "+tile1);
								}
							} else 	if ((layer1 < 0) || (connections[tile][layer][dir][layer1] > 1.0)) { // no connections in this direction - add border using same disparity as this tile
								if (tile_stat[tile1] == TILE_EMPTY) { // not yet set (otherwise keep)
									tile_stat[tile1] = TILE_BORDER_FLOAT; //TILE_BORDER;
									disparity[tile1] = disparity[tile];
									if (Double.isNaN(disparity[tile1])) {
										System.out.println("clusterizeFgBg(): 2. Disparity is set NaN for tile "+tile1);
									}
								}
							} else { // connection exists, add tile1/layer1 as a candidate to the list
								tile_stat[tile1] = TILE_CANDIDATE; // to avoid suggesting the same tile multiple times
								if (Double.isNaN(disparity[tile1])) {
									System.out.println("clusterizeFgBg(): *3. Disparity is set NaN for tile "+tile1);
								}
								tile_layer_list.add(new Point(tile1, layer1));
							}
						}
					}						
				}
				//				}
			} // while (!tile_layer_list.isEmpty())
			if (num_tiles > 0) { // generate/add new cluster, update
				// Create and crop cluster,(update neighbors?)
				int min_x = tilesX, min_y = tilesY, max_x = 0, max_y=0;
				for (int ty = 0; ty < tilesY; ty++) {
					for (int tx = 0; tx < tilesX; tx++) {
						int tile = ty*tilesX + tx;
						if (tile_stat[tile] != TILE_EMPTY) {
							if (ty > max_y) max_y = ty;
							if (ty < min_y) min_y = ty;
							if (tx > max_x) max_x = tx;
							if (tx < min_x) min_x = tx;
						}
					}
				}
				int width =  max_x - min_x + 1;
				int height = max_y - min_y + 1;
				
				Rectangle bounds = new Rectangle(min_x, min_y, width, height);
				double  [] disparity_crop = new double [width * height]; 
				boolean [] border_crop = new boolean   [width * height];
				double [] wdir = {1.0,0.7,1.0,0.7,1.0,0.7,1.0,0.7}; // weights for directions
				for (int dty = 0; dty < height; dty++) {
					int ty = min_y + dty;
					for (int dtx = 0; dtx < width; dtx++) {
						int tdest = dty*width + dtx;
						int tsrc = ty * tilesX + min_x + dtx;
						// replace disparity for the floating border tiles with the weighted average of non-floating 
						if (tile_stat[tsrc] == TILE_BORDER_FLOAT) { // average disparity from non-float defined tiles
							// TODO: use best plane fit (for gradients)
							// TODO: remove some "inner" border tiles?
							double sw = 0.0;
							double swd = 0.0;
							for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
								int tile1 = tn.getNeibIndex(tsrc, dir); // should always be > 0 here
								if (tile1 >= 0) {
									if ((tile_stat[tile1] == TILE_BORDER) || (tile_stat[tile1] == TILE_CONFIRMED)) {
										sw += wdir[dir];
										swd += wdir[dir] * disparity[tile1];
									}
								}
							}
							if (sw > 0.0) {
								disparity[tsrc] = swd/sw;
							}
						}
						border_crop[tdest] = isBorder(tile_stat[tsrc]); //  == TILE_BORDER;
						disparity_crop[tdest] = (tile_stat[tsrc] == TILE_EMPTY) ? Double.NaN : disparity[tsrc];
					}
				}
				TileCluster tileCluster = (new TileCluster(
						bounds,
						(debug_index? cluster_list.size(): -1),
						border_crop,
						disparity_crop));
				cluster_list.add(tileCluster);
				// update

				
			} else {
				if (debugLevel >-2) {
					System.out.println("clusterizeFgBg(): Empty cluster started from tile="+p0.x+", layer="+p0.y);
				}
			}
		} // while (true) { - generating clusters
		
		// consolidate clusters "good enough"
		int [] comb_clusters = new int [cluster_list.size()];
		Arrays.fill(comb_clusters,-1);
		int this_combo = 0;
		for (; ; this_combo++) {
			// find first unassigned cluster
			int index_first = -1;
			for (int i = 0; i < comb_clusters.length; i++) {
				if (comb_clusters[i] < 0) {
					index_first = i;
					break;
				}
			}
			if (index_first < 0) {
				break; // all clusters assigned
			}
			comb_clusters[index_first] = this_combo;
			for (int index_other = index_first; index_other < comb_clusters.length; index_other++) if (comb_clusters[index_other] < 0) {
				// check to intersection with all prior clusters in this combo
				candidate_cluster:
				{
					Rectangle new_bounds = cluster_list.get(index_other).getBounds();
					for (int index_already = index_first; index_already < index_other; index_already++) if (comb_clusters[index_already] == this_combo) {
						if (cluster_list.get(index_already).getBounds().intersects(new_bounds)) {
							break candidate_cluster; // intersects - skip it
						}
					}
					comb_clusters[index_other] = this_combo;
				}
			}
		}
		TileCluster [] consolidated_clusters = new TileCluster[this_combo];
		Rectangle full_tiles = new Rectangle(0, 0, tilesX, tilesY);
		for (int i = 0; i < this_combo; i++) {
			consolidated_clusters[i] = new TileCluster(
					full_tiles,
					(debug_index? 0:-1),
					null,
					null);
		}
		for (int i = 0; i < comb_clusters.length; i++) {
			consolidated_clusters[comb_clusters[i]].add(cluster_list.get(i));
		}

		if (debugLevel > 0) {
			double [][] dbg_img =     new double[this_combo][tiles];
			double [][] dbg_borders = new double[this_combo][tiles];
			double [][] dbg_index = null;
			if (debug_index) {
				dbg_index = new double[this_combo][tiles];
			}
			for (int n = 0; n < dbg_img.length; n++) {
				for (int i = 0; i < tiles;i++) {
					dbg_img[n][i] =     consolidated_clusters[n].getDisparity()[i];
					dbg_borders[n][i] = consolidated_clusters[n].getBorder()[i]? 1.0:0.0;
					if (dbg_index != null) {
						double d = consolidated_clusters[n].getClusterIndex()[i];
						dbg_index[n][i] = (d >=0)? d : Double.NaN;
					}
				}				
			}			
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"cluster_disparity");
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_borders,
					tilesX,
					tilesY,
					true,
					"cluster_borders");
			if (dbg_index != null) {
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_index,
						tilesX,
						tilesY,
						true,
						"cluster_indices");
			}
			
		}
		return consolidated_clusters;
	}

}
