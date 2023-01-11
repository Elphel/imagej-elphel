/**
 ** TexturedModel - Generate 3D mesh with textures
 **
 ** Copyright (C) 2022 Elphel, Inc.
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
package com.elphel.imagej.tileprocessor;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;

import org.json.JSONException;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.gpu.TpTask;
import com.elphel.imagej.x3d.export.GlTfExport;
import com.elphel.imagej.x3d.export.TriMesh;
import com.elphel.imagej.x3d.export.WavefrontExport;
import com.elphel.imagej.x3d.export.X3dOutput;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;

public class TexturedModel {
	public static final int THREADS_MAX = 100;  // maximal number of threads to launch
	// Some of below may be obsolete
	public static final int TILE_EMPTY =          0; 
	public static final int TILE_BORDER =         1; // tile shared between meshes, border with known (not float) disparity
	public static final int TILE_BORDER_FLOAT =   2; // border tile with disparity calculated from neighbors
	public static final int TILE_CONFIRMED =      3; // internal (not border) mesh tile
	public static final int TILE_CANDIDATE =      4; // not used 
	public static final int CLUSTER_NAN =        -2; // disparity is NaN 
	public static final int CLUSTER_UNASSIGNED =  -1; // not yet assinged (>=0 - cluster number)
	public static final int []    NUM_NEIBS_FROM_BITS = new int [512]; 
	
	public static final int TILE_IS_FG_WEAK =    0; 
	public static final int TILE_IS_FG_STRONG =  1; 
	public static final int TILE_HAS_BG_WEAK =   2; 
	public static final int TILE_HAS_BG_STRONG = 3;
	public static final int TILE_KEEP =          4;
	public static final int TILE_STITCH =        5; 
	public static final int TILE_STITCHED =      6; 
	public static final int TILE_BOOLEANS =      TILE_STITCHED + 1;

	public static final int PIX_TRIM_FG =        0; 
	public static final int PIX_HAS_BG =         1;
	public static final int PIX_EDGE_FG =        2; 
	
	public static final int PIX_BOOLEANS =       PIX_EDGE_FG + 1;
	

	
	
	public static boolean isBorder(int d) {
		return (d==TILE_BORDER) || (d==TILE_BORDER_FLOAT);
	}
	
	public static TileCluster [] clusterizeFgBgOld( //
			final int          tilesX,
			final double [][]  disparities, // may have more layers
			final boolean []   blue_sky, // use to expand background by blurring available data?
			final int          blue_sky_layer,
			final int          blue_sky_below,
			final boolean []   selected, // to remove sky (pre-filter by caller, like for ML?)
			final double       disp_adiffo,
			final double       disp_rdiffo,
			final double       disp_adiffd,
			final double       disp_rdiffd,
			final double       disp_fof,    // enable higher difference (scale) for fried of a friend 
			final int          cluster_gap, // gap between clusters 
			final int          debugLevel) {
		final double disp_border = disp_fof;
		final int tiles = disparities[0].length;
		final int tilesY = tiles/tilesX;
		final int layers = disparities.length;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		final double [][][][] connections = new double [tiles][][][];
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
		// calculate "connections - per tile, per layer, per direction (1 of the first 4), per target layer - normalized difference difference 
		final int dbg_tile = (debugLevel>0)? 1090:-1; // 977 : -1;
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
								boolean is_bs = (layer == blue_sky_layer) && blue_sky[tile];
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
												boolean is_bs1 = (layer1 == blue_sky_layer) && blue_sky[tile1];
											    if (is_bs1 == is_bs) { // do not fix bs/no bs
											    	connections[tile][layer][dir][layer1] = Math.abs(disparities[layer][tile] - disparities[layer1][tile1])/max_disp_diff;
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

		ai.set(0);
		// Fill in opposite connections by combining opposite directions
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
		// calculate total number of connections (w/o fof) with value < 1.0, increment once
		// per direction even if there are multiple connected layers
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (connections[tile] != null) {
						if (tile==dbg_tile) {
							System.out.println("clusterizeFgBg().3: tile="+tile);
						}
						for (int layer = 0; layer < layers; layer++) if (connections[tile][layer] != null){
							num_neibs_dir[tile][layer] = 1; // center tile
							for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
								if (connections[tile][layer][dir] != null) {
									for (int layer1 = 0; layer1 < layers; layer1++) {
										if (connections[tile][layer][dir][layer1] <= 1.0) { // Double.NaN - OK
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
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"num_neibs",
					dbg_titles);
			ShowDoubleFloatArrays.showArrays(
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
								if (nn == (TileNeibs.DIRS + 1)) { // No sense to look more - it can not be > 9
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
			boolean sky_cluster = (best_layer == blue_sky_layer) && blue_sky[best_tile];
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
					// see if it has confirmed (or margin) same-layer neighbor 1.5 times more disparity difference
					// in such case - mark this tile TILE_BORDER_FLOAT (to split mesh)  and unconfirm
					for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
						int tile1 = tn.getNeibIndex(tile, dir); // should always be > 0 here
						//is it a border tile or already confirmed one ?
						if (isBorder(tile_stat[tile1]) || (tile_stat[tile1] == TILE_CONFIRMED)){
							double mid_disp = Math.max(0.0, 0.5*(disparities[layer][tile] + disparity[tile1]));
							double max_disp_diff = ((dir & 1) == 0) ?
									(disp_adiffo + mid_disp * disp_rdiffo) :
										(disp_adiffd + mid_disp * disp_rdiffd);
							max_disp_diff *= isBorder(tile_stat[tile1])? disp_border : disp_fof;
							if ((Math.abs(disparities[layer][tile] - disparity[tile1])/max_disp_diff) > 1.0){
								confirm = false;// too large diff
								// make it a border tile, but keep disparity
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
							// check - maybe there is another layer (but why both layer1 and layer2?)
							// do we need to hadle sky here too?
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
								for (int l = 0; l < layers; l++) if (connections[tile][layer][dir][l] <= 1.0){ // connections between sky and not sky will be NaN here
									// verify it was not already used. If used - use used disparity for the border 									
									// **** (select correct layer for that border disparity
									if (ncluster[tile1][l] == CLUSTER_UNASSIGNED) {
										if ((layer1 < 0) || (connections[tile][layer][dir][l] < connections[tile][layer][dir][layer1])) {
											layer1 = l; // best neighbor
										}
									} else { // already assigned to some cluster that was split
										if ((layer_assigned < 0) || (connections[tile][layer][dir][l] < connections[tile][layer][dir][layer_assigned])) {
											layer_assigned = l; // best of assigned layers
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
							} else 	if ((layer1 < 0) || !(connections[tile][layer][dir][layer1] <= 1.0)) { // no connections in this direction - add border using same disparity as this tile
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
			// finished current cluster, see if it has any tiles
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
				if (sky_cluster) { // increase bounding box for sky cluster
					min_y = 0;
					max_y+= blue_sky_below;
					min_x = 0;
					max_x = tilesX -1;
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
						cluster_list.size(), // (debug_index? cluster_list.size(): -1),
						border_crop,
						null,     // int []     border_int,           // will replace border? Provide on-the-fly? 
						0,        // int        border_int_max,       // outer border value
						disparity_crop,
						sky_cluster));     // boolean is_sky));
				cluster_list.add(tileCluster);
				// update
			} else {
				if (debugLevel >-2) {
					System.out.println("clusterizeFgBg(): Empty cluster started from tile="+p0.x+", layer="+p0.y);
				}
			}
		} // while (true) { - generating clusters
		// consolidate clusters "good enough", use bounding box intersections, add cluster_gap to grow extra tiles by Gaussian
		// cluster_gap
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
					Rectangle new_bounds = cluster_list.get(index_other).getBounds(cluster_gap); // cluster_gap should be 2x
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
					null,     // int []     border_int,           // will replace border? Provide on-the-fly? 
					0,        // int        border_int_max,       // outer border value
					null,
					false); // boolean is_sky));
		}
		for (int i = 0; i < comb_clusters.length; i++) {
			consolidated_clusters[comb_clusters[i]].add(cluster_list.get(i));
		}
		// incorrectly combined - first combo has only one cluster and NaN in the areas where other clusters could fit
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
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"cluster_disparity");
			ShowDoubleFloatArrays.showArrays(
					dbg_borders,
					tilesX,
					tilesY,
					true,
					"cluster_borders");
			if (dbg_index != null) {
				ShowDoubleFloatArrays.showArrays(
						dbg_index,
						tilesX,
						tilesY,
						true,
						"cluster_indices");
			}
			
		}
		return consolidated_clusters;
	}

	// generate/update number of neighbors to select clusters' seeds
	/**
	 * Generate and later update array on number of neighbors to find candidates for the next
	 * texture meshes. Initially calculated for all tiles, then updated by after each cluster
	 * extraction in the area possibly affected by the cluster (cluster bounds extended by one
	 * in each direction).
	 * 
	 * @param num_neibs_dir    [layers][tiles] - array of connections (defined by close disparities)
	 *                         per layer, per tile. Undefined tile has 0, isolated tile - 1, maximal
	 *                         (with all 8 neighbors present) - 9. Should be initialized to
	 *                         [layers][tiles].
	 * @param bounds           Rectangle instance specifying the new cluster bounds. If null (used
	 *                         during first run), all tiles are processed
	 * @param disparity_layers [layers][tiles] multi-layer disparity array with NaN
	 *                         for non-existing tiles 
	 * @param blue_sky         per-tile array of the "blue sky" boolean array 
	 * @param blue_sky_layer   layer for the "blue sky" data (should be on a single layer)
	 * @param disp_adiffo      absolute disparity difference for connecting tiles in ortho
	 *                         directions (used for initial connections estimations)
	 * @param disp_rdiffo      relative disparity difference (added to absolute being multiplied
	 *                         by the tile disparity)
	 * @param disp_adiffd      absolute disparity difference for diagonal tiles.
	 * @param disp_rdiffd      relative disparity difference for diagonal tiles.
	 * @param disp_fof         >=1.0 - increased inter-tile tolerance for a friend-of-a-friend.
	 *                         In current code just scales calculated (absolute+relative) tolerance
	 *                         for all steps after initial connections.
	 * @param tilesX           horizontal array dimension 
	 * @param debugLevel       debug level - controls generation of images
	 */
	public static void updateSeeds( // and update num_neibs_dir
			final int [][]     num_neibs_dir, // [tile][layer]
			final Rectangle    bounds,    // null - all
			final double [][]  disparity_layers, // [layer][tile]should not have same tile disparity on multiple layers
			final boolean []   blue_sky, // use to expand background by blurring available data?
			final int          blue_sky_layer,
			final double       disp_adiffo,
			final double       disp_rdiffo,
			final double       disp_adiffd,
			final double       disp_rdiffd,
			final double       disp_fof,    // enable higher difference (scale) for fried of a friend 
			final int          tilesX,
			final int          debugLevel) {
		final int tiles = disparity_layers[0].length;
		final int tiles_wnd = (bounds == null) ? tiles : (bounds.width * bounds.height);
		final int tilesY = tiles/tilesX;
		final int layers = disparity_layers.length;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		final double [][][][] connections = new double [tiles][][][];
		if (NUM_NEIBS_FROM_BITS[511] == 0) {
			for (int i = 0; i < NUM_NEIBS_FROM_BITS.length; i++) {
				for (int d = i; d != 0; d>>=1) {
					NUM_NEIBS_FROM_BITS[i]+=(d & 1);
				}
			}
		};
		final Rectangle bounds_ext = (bounds != null) ?((new Rectangle(bounds.x-1, bounds.y-1, bounds.width+2, bounds.height + 2)).
				intersection(new Rectangle(tilesX, tilesY))) : null;
		final Rectangle bounds_ext2 = (bounds != null) ?((new Rectangle(bounds.x-2, bounds.y-2, bounds.width+4, bounds.height + 4)).
				intersection(new Rectangle(tilesX, tilesY))) : null;
		// calculate "connections - per tile, per layer, per direction (1 of the first 4), per target layer - normalized difference difference 
		final int dbg_tile = (debugLevel>0)? 1090:-1; // 977 : -1;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile_wnd = ai.getAndIncrement(); tile_wnd < tiles_wnd; tile_wnd = ai.getAndIncrement()) {
						int tile = tile_wnd;
						if (bounds_ext2 != null) {
							int tileX = bounds_ext2.x + tile_wnd % bounds_ext2.width;
							int tileY = bounds_ext2.y + tile_wnd / bounds_ext2.width;
							tile = tileY * tilesX + tileX;
						}
						if (tile==dbg_tile) {
							System.out.println("updateSeeds().1: tile="+tile);
						}
						for (int layer = 0; layer < layers; layer++) {
							if (!Double.isNaN(disparity_layers[layer][tile])) {
								if (connections[tile] == null) {
									connections[tile] = new double[layers][][];
								}
								boolean is_bs = (layer == blue_sky_layer) && blue_sky[tile];
								connections[tile][layer] = new double [TileNeibs.DIRS][]; // leave room for future symmetry
								for (int dir = 0; dir < TileNeibs.DIRS/2; dir++) {
									int tile1 = tn.getNeibIndex(tile, dir);
									if (tile1 >= 0) {
										for (int layer1 = 0; layer1 < layers; layer1++) {
											if (!Double.isNaN(disparity_layers[layer1][tile1])) {
												if (connections[tile][layer][dir] == null) {
													connections[tile][layer][dir] = new double[layers];
													Arrays.fill(connections[tile][layer][dir],Double.NaN);
												}
												double mid_disp = Math.max(0.0, 0.5*(disparity_layers[layer][tile] + disparity_layers[layer1][tile1]));
											    double max_disp_diff = ((dir & 1) == 0) ?
											    		(disp_adiffo + mid_disp * disp_rdiffo) :
											    			(disp_adiffd + mid_disp * disp_rdiffd);
												boolean is_bs1 = (layer1 == blue_sky_layer) && blue_sky[tile1];
											    if (is_bs1 == is_bs) { // do not mix bs/no bs
											    	connections[tile][layer][dir][layer1] = Math.abs(disparity_layers[layer][tile] - disparity_layers[layer1][tile1])/max_disp_diff;
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

		ai.set(0);
		// Fill in opposite connections by combining opposite directions
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile_wnd = ai.getAndIncrement(); tile_wnd < tiles_wnd; tile_wnd = ai.getAndIncrement()) {
						int tile = tile_wnd;
						if (bounds != null) {
							int tileX = bounds_ext.x + tile_wnd % bounds_ext.width;
							int tileY = bounds_ext.y + tile_wnd / bounds_ext.width;
							tile = tileY * tilesX + tileX;
						}
						if (tile==dbg_tile) {
							System.out.println("clusterizeFgBg().2: tile="+tile);
						}
						for (int layer = 0; layer < layers; layer++) if (!Double.isNaN(disparity_layers[layer][tile])) {
							if ((connections[tile] != null) && (connections[tile][layer] != null)) {
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
						} // for (int layer = 0; layer < layers; layer++) if (!Double.isNaN(disparity_layers[layer][tile])) {
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		// extend bounds by 1 each side
		ai.set(0);
		// calculate total number of connections (w/o fof) with value < 1.0, increment once
		// per direction even if there are multiple connected layers
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile_wnd = ai.getAndIncrement(); tile_wnd < tiles_wnd; tile_wnd = ai.getAndIncrement()) {
						int tile = tile_wnd;
						if (bounds_ext != null) {
							int tileX = bounds_ext.x + tile_wnd % bounds_ext.width;
							int tileY = bounds_ext.y + tile_wnd / bounds_ext.width;
							tile = tileY * tilesX + tileX;
						}
						Arrays.fill(num_neibs_dir[tile], 0);
						if (connections[tile] != null) {
							if (tile==dbg_tile) {
								System.out.println("updateSeeds().3: tile="+tile);
							}
							for (int layer = 0; layer < layers; layer++) if (connections[tile][layer] != null){
								num_neibs_dir[tile][layer] = 1; // center tile
								for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
									if (connections[tile][layer][dir] != null) {
										for (int layer1 = 0; layer1 < layers; layer1++) {
											if (connections[tile][layer][dir][layer1] <= 1.0) { // Double.NaN - OK
												num_neibs_dir[tile][layer] |= 2 << dir; // 9 bits
												break; // increment once per dir
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
		return;
	}

	/**
	 * Find the next seed tile to build a new mesh (tile cluster). Search concludes with the maximal
	 * Number of neighbors tile or the first encountered 9 neighbors (self plus 8 around) as it is the 
	 * maximal possible number.
	 *                         
	 * @param disparity_layers [layers][tiles] multi-layer disparity array with NaN
	 *                         for non-existing tiles 
	 * @param num_neibs_dir    [layers][tiles] - array of connections (defined by close disparities)
	 *                         per layer, per tile. Undefined tile has 0, isolated tile - 1, maximal
	 *                         (with all 8 neighbors present) - 9. Should be initialized to
	 *                         [layers][tiles].
	 * @param tile_start       tile index to begin search. Provide tile where the last search ended,
	 *                         the search will wrap over the full array size and start from 0.
	 * @param tilesX           horizontal array dimension. 
	 * @return                 a pair {tile, layer} or null if there are no tiles left
	 */
	public static int [] getNextSeed(
		final double [][] disparity_layers, //
		final int [][]    num_neibs_dir, // [tile][layer]
		final int         tile_start,
		final int         tilesX)
	{
		final int tiles = num_neibs_dir.length;
		final int tilesY = tiles/tilesX;
		final int layers =disparity_layers.length;
		final TileNeibs tn =    new TileNeibs(tilesX, tilesY);
		int best_tile = -1, best_layer = -1, nn=0;
		find_start:
		{
			int tile_end = tile_start + tiles;
			for (int tile1=tile_start; tile1 < tile_end; tile1++) {
				int tile = (tile1 >= tiles) ? (tile1 - tiles) : tile1;
				if (!tn.isBorder(tile)) { // do not start on the border
					for (int layer = 0; layer < layers; layer++) {
						int n_neibs = NUM_NEIBS_FROM_BITS[num_neibs_dir[tile][layer]];
//						if ((ncluster[tile][layer] == CLUSTER_UNASSIGNED) && (n_neibs > nn)) { // not yet assigned and >=0 neibs
						if (!Double.isNaN(disparity_layers[layer][tile]) && (n_neibs > nn)) { // not yet assigned and >=0 neibs
							nn = n_neibs;
							best_tile = tile;
							best_layer = layer;
							if (nn == (TileNeibs.DIRS + 1)) { // No sense to look more - it can not be > 9
								break find_start;
							}
						}
					}
				}
			}
		}
		if (best_tile < 0) {
			return null;
		}
		return new int [] {best_tile, best_layer};
	}
	
	/**
	 * Create initial cluster of connected tiles without provisions for overlap resolution
	 * @param disparity_layers [layers][tiles] multi-layer disparity array with NaN
	 *                         for non-existing tiles 
	 * @param seams_layers     [layers][tiles] marked seams from previous clusters
	 * @param seams            [tiles] seams corresponding to result disparity
	 * @param start_layer      seed layer to start growing a cluster
	 * @param start_tile       seed tile to start growing a cluster
	 * @param blue_sky         per-tile array of the "blue sky" boolean array 
	 * @param blue_sky_layer   layer for the "blue sky" data (should be on a single layer)
	 * @param disp_adiffo      absolute disparity difference for connecting tiles in ortho
	 *                         directions (used for initial connections estimations)
	 * @param disp_rdiffo      relative disparity difference (added to absolute being multiplied
	 *                         by the tile disparity)
	 * @param disp_adiffd      absolute disparity difference for diagonal tiles.
	 * @param disp_rdiffd      relative disparity difference for diagonal tiles.
	 * @param disp_fof         >=1.0 - increased inter-tile tolerance for a friend-of-a-friend.
	 *                         In current code just scales calculated (absolute+relative) tolerance
	 *                         for all steps after initial connections.
	 * @param jump_r           "jump" over small gaps when building initial clusters. jump_r == 2 
	 *                         allows jumping to other tiles in 5x5 square around the defined tile,
	 *                         jump_r == 3 - inside 7x7 square. The jump destination should not have
	 *                         any already defined tiles among 8 neighbors (in a 3x3 square). 
	 * @param disp_adiffj      maximal absolute disparity difference for the "jumps".
	 * @param disp_rdiffj      maximal relative disparity difference for the "jumps".       
	 * @param tilesX           horizontal array dimension 
	 * @param debugLevel       debug level - controls generation of images
	 * @return                 [tiles] disparity array of the selected tiles. NaN for unselected
	 */
	public static double[] buildInitialCluster(
		final double [][]     disparity_layers, // should not have same tile disparity on multiple layers
		final int [][]        seams_layers,
		final int []          seams,
		final int             start_layer,
		final int             start_tile,
		final boolean []      blue_sky, // Do not mix bs/no_bs in the same cluster 
		final int             blue_sky_layer,
		final double          disp_adiffo,
		final double          disp_rdiffo,
		final double          disp_adiffd,
		final double          disp_rdiffd,
		final double          disp_fof,    // enable higher difference (scale) for friend of a friend 
		final int             jump_r,      // jump over alien/NaN disparities
		final double          disp_adiffj,
		final double          disp_rdiffj,
		final int             tilesX,
		final int             debugLevel)
	{
		final boolean        is_sky_cluster = (start_layer == blue_sky_layer) &&  blue_sky[start_tile];
		final int  num_layers = disparity_layers.length;
		final int  tiles =      disparity_layers[0].length;
		final int  tilesY =     tiles/tilesX;         
		double disparity[] =    new double[tiles]; // current cluster disparities
		Arrays.fill(disparity, Double.NaN);
		final TileNeibs tn =    new TileNeibs(tilesX, tilesY);
		ArrayList<Integer> tile_layer_list = new ArrayList<Integer>(); // pair of int x tile, int y - layer
		tile_layer_list.add(start_tile);
		disparity[start_tile] = disparity_layers[start_layer][start_tile];
		seams[start_tile] =     seams_layers[start_layer][start_tile];		
		final Thread[] threads =    ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai =    new AtomicInteger(0);
		final AtomicInteger alayer_tile = new AtomicInteger(-1);
		while (true) {
			while (!tile_layer_list.isEmpty()) {
				int tile = tile_layer_list.remove(0);
				double disp = disparity[tile]; 
				double delta_disp_ortho = (disp_adiffo + disp * disp_rdiffo) * disp_fof;
				double delta_disp_diag =  (disp_adiffd + disp * disp_rdiffd) * disp_fof;
				for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
					int tile1 = tn.getNeibIndex(tile, dir); // should always be > 0 here
//					if (tile1 == 2849) {
//						System.out.println("buildInitialCluster(): tile1 = "+tile1);
//					}
					if (tile1 >= 0) {
						double delta_disp = ((dir & 1) == 0)? delta_disp_ortho : delta_disp_diag;
						// see if it already has a tile of the same cluster in this direction
						if (!Double.isNaN(disparity[tile1])) { // already assigned to this cluster
							if (Math.abs(disparity[tile1] - disp) < (disp_fof *  delta_disp)) {
								continue; // many neighbors fall here - already assigned at fit
							}
						}
						// find best fit (then reconsider previous assignment)
						int blayer = -1;
						double bdisp = Double.NaN;
						for (int layer1 = 0; layer1 < num_layers; layer1++) {
							double disp1 = disparity_layers[layer1][tile1];
							boolean is_bs1 =  (layer1 == blue_sky_layer) && blue_sky[tile1];
							if (!Double.isNaN(disp1) && (is_bs1 == is_sky_cluster)) {
								if ((blayer < 0) || ((Math.abs(disp1 - disp) < Math.abs(bdisp - disp)))) {
									blayer = layer1;
									bdisp = disp1;
								}
							}
						}
						if (blayer >= 0) {
							double mid_disp = Math.max(0.0, 0.5*(disp + bdisp));
							double max_disp_diff = ((dir & 1) == 0) ?
									(disp_adiffo + mid_disp * disp_rdiffo) :
										(disp_adiffd + mid_disp * disp_rdiffd);
							if ((Math.abs(disp - bdisp)/max_disp_diff) <= 1.0){ // fits
								if (!Double.isNaN(disparity[tile1])) { // already assigned to this cluster
									if (bdisp > disparity[tile1]) { // new found is FG (higher disparity than the old one) -> replace old
										disparity[tile1] = bdisp;
										seams[tile1] = seams_layers[blayer][tile1];
									}
									// replaced assignment - do not increase number of tiles
								} else {
									disparity[tile1] = bdisp;
									seams[tile1] = seams_layers[blayer][tile1];
								}
								tile_layer_list.add(tile1);
							}
						}
					}
				}
			} // while (!tile_layer_list.isEmpty()) {
			// Try jumping over;
	        ai.set(0);
	        alayer_tile.set(-1);
	        // calculate total number of connections (w/o fof) by combining opposite directions
	        for (int ithread = 0; ithread < threads.length; ithread++) {
	            threads[ithread] = new Thread() {
	                public void run() {
	                    for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
	                    	if (alayer_tile.get() >= 0) {
	                    		break;
	                    	}
	                    	double disp = disparity[tile];
	                    	if (!Double.isNaN(disp)) {
	                    		double delta_disp = disp_adiffj + Math.max(0, disp) * disp_rdiffj;
	                    		int best_r2 =  0, best_tile = -1, best_layer = -1;
	                    		for (int dy = -jump_r; dy <= jump_r; dy++) {
	                    			boolean good_col = Math.abs(dy) > 1;
	                    			for (int dx = -jump_r; dx <= jump_r; dx++) {
	                    				// do not use 8 neighbors - they should be used during wave
	                    				if (good_col || (Math.abs(dx) > 1)) {
	                    					int tile1 = tn.getNeibIndex(tile, dx, dy);
	                    					if ((tile1 >=0) && Double.isNaN(disparity[tile1])) {
	                    						double disp_min = disp - delta_disp; 
	                    						double disp_max = disp + delta_disp; 
	                    						for (int layer = 0; layer < num_layers; layer++) {
	                    							boolean is_bs1 =  (layer == blue_sky_layer) && blue_sky[tile1];
	                    							if (is_bs1 == is_sky_cluster) {
	                    								double disp1 = disparity_layers[layer][tile1];
	                    								if (!Double.isNaN(disp1) && (disp1 >= disp_min) && (disp1 <= disp_max)) {
	                    									int r2 = dy*dy+dx*dx;
	                    									if ((best_tile < 0) || (r2 < best_r2)) {
	                    										// check that new tile does not have selected neighbors already
	                    										boolean no_sel_neibs = true;
	                    										for (int dir1 = 0; dir1 < 8; dir1++) {
	                    											int tile2 = tn.getNeibIndex(tile1, dir1);
	                    											if ((tile2 >= 0) && !Double.isNaN(disparity[tile2])) {
	                    												no_sel_neibs=false;
	                    												break;
	                    											}
	                    										}
	                    										if (no_sel_neibs) {
	                    											best_r2 = r2;
	                    											best_tile = tile1;
	                    											best_layer = layer;
	                    										}
	                    									}
	                    								}
	                    							}
	                    						}
	                    					}
	                    				}
	                    			}
	                    		}
	                    		if (best_tile >= 0) {
	                    			alayer_tile.getAndSet(best_layer * tiles + best_tile);
	                    			break;
	                    		}
	                    	}
	                    }
	                }
	            };
	        }		      
	        ImageDtt.startAndJoin(threads);
	        int slt = alayer_tile.get();
	        if (slt < 0) {
	        	break;
	        }
	        int sl = slt / tiles;
	        int st = slt % tiles;
			//alayer_tile.getAndSet(best_layer * tiles + best_tile);
			tile_layer_list.add(st);
     		disparity[st] =  disparity_layers[sl][st];
     		seams[st] = seams_layers[sl][st];
		} // while (true) {
		return disparity;
	}

	/**
	 * Use initial cluster (source_disparity) to output cluster with appropriate border tiles (currently
	 * trying with two border layers around the selected disparity tiles). Foreground has higher priority,
	 * so if higher disparity conflicts with lower disparity it pushes lower disparity tiles until higher
	 * disparity with specified (now 2) layers of extrapolated (not from disparity_layers array /
	 * source_disparity). After "pushing away" conflicting tiles, additional layers (same number) over
	 * defined (in source_disparity) tiles are removed from source_disparity. Only remaining (not
	 * marked as borders) tiles from the original source_disparity are removed from disparity_layers
	 * used to generate next clusters. Each extracted tileCluster preserves Rectangle bounds later used
	 * to combine multiple clusters in the same full frame array, disparity and bounds arrays are defined
	 * in rectangular arrays corresponding to bounds.   
	 *  
	 * @param cluster_list     ArrayList of previously defined tileCluster instances. Used to get this
	 *                         tileCluster index, the new tileCluster is added to the list. 
	 * @param is_sky_cluster   if this cluster is a special type cluster - "blue sky". This property is
	 *                         saved in the tileCluster instance. If, additionally, the next parameter
	 *                         is >=0 the cluster bounds are extended: left, right and top bounds - to the
	 *                         frame edges, and the bottom one - blue_sky_below below the highest tile Y.
	 * @param blue_sky_below   (only if is_sky_cluster) if >=0 - extend bounds, if <0 treat blue_sky bounds
	 *                         same as regular ones.  
	 * @param disparity_layers [layers][tiles] multi-layer disparity array with NaN for non-existing tiles
	 * @param seams_layers     [layers][tiles] multi-layer border levels to communicate BG stitches
	 * @param seams            [tiles] seams corresponding to source_disparity
	 * @param source_disparity [tiles] continuous tiles used to generate this tileCluster that uses only
	 *                         of source_disparity and extrapolated borders.
	 * @param max_neib_lev     maximal neighbor level of borders used to generate internal int [] neib_lev
	 *                         array. Neighbor level 0 corresponds to internal tiles that uses original
	 *                         disparity, level 1 - neighbors of level 0 on the border or in conflict with
	 *                         level 0. Level 2 - neighbors of 1 or in conflict with level 0, and so all.   
	 * @param disp_adiff       absolute disparity difference for connecting tiles in any direction, 
	 *                         should include disp_fof used for initial cluster generation
	 * @param disp_rdiff       relative disparity difference for connecting tiles in any direction.
	 * @param tilesX           horizontal array dimension 
	 * @param debugLevel       debug level - controls generation of images
	 * @return                 TileCluster instance generated from the initial source_disparity[] array.
	 */
	public static TileCluster [] buildTileCluster(
			// used disparity_layers will be set to Double.NaN
			// make it in a separate method?
			final ArrayList <TileCluster> cluster_list,
			final boolean         is_sky_cluster, // this is a blue sky cluster, mark as such and extend bounds
			final int             blue_sky_below, // extend bounds down from the blue sky lower 
			final double [][]     disparity_layers, // should not have same tile disparity on multiple layers
			final int [][]        seams_layers,
	        final int []          seams,
			final double []       source_disparity, // should not have same tile disparity on multiple layers
			final int             max_neib_lev,
			final double          disp_adiff, // should already include disp_fof,
			final double          disp_rdiff,
			final int             tilesX,
			int             debugLevel_in)
	{
		final int dbg_tile = 2410; // 858; // -3868; // 28+48*80;
		//		final int       num_layers = disparity_layers.length;
		final int        tiles =        source_disparity.length;
		final int        tilesY =       tiles/tilesX;
		final int []     neib_lev =     new int   [tiles];
		final boolean [] initial_seam = new boolean [tiles]; // neib_lev was set from initial seam
		final boolean [] new_seam =     new boolean [tiles];
		final double []  disparity =    new double[tiles]; // current cluster disparities
		final double []  max_neib =     new double[tiles]; // maximal disparity of neibs
		if (!Double.isNaN(source_disparity[dbg_tile])) {
			System.out.println("buildTileCluster(): source_disparity["+dbg_tile+"] is not NaN");
			debugLevel_in = 2;
		}
		final int debugLevel = debugLevel_in; 
		if (debugLevel > 1) {
			String [] dbg_titles = {"Disparity","Seams"};
			double [] dbg_seams = new double [tiles];
			
			for (int i = 0; i < tiles; i++) {
				dbg_seams[i] = 10*seams[i];
			}
			
			double [][] dbg_img = {
					source_disparity,
					dbg_seams};
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"disparity-seams-"+String.format("%02d", cluster_list.size()),
					dbg_titles);
		}		
		
		
//		final boolean [] disp_mod = new boolean[tiles]; // disparity modified from surce_disparity
		Arrays.fill(neib_lev, -1);
		System.arraycopy(source_disparity, 0, disparity, 0, tiles);
		final TileNeibs tn =    new TileNeibs(tilesX, tilesY);
		ArrayList<Integer> loc_list = new ArrayList<Integer>();
		ArrayList<Integer> lor_list = new ArrayList<Integer>();
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ati = new AtomicInteger(0);
		ai.set(0);
		ati.set(0);
		// create list of conflicts and 1 tile around defined, mark known disparity in neib_lev[]
		final ArrayList<ArrayList<Integer>> loc_multi = new ArrayList<ArrayList<Integer>>(threads.length);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			loc_multi.add(new ArrayList<Integer>());
			threads[ithread] = new Thread() {
				public void run() {
					ArrayList<Integer> loc_this = loc_multi.get(ati.getAndIncrement());
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
						if (tile == dbg_tile) {
							System.out.println("buildTileCluster().01: tile="+tile);
						}
						double max_n = Double.NaN;
						for (int dir = 0; dir < 8; dir++) {
							int tile1 = tn.getNeibIndex(tile, dir);
							if (tile1 >= 0) { // has defined neighbor
								double disp1 =  disparity[tile1];
								if (!Double.isNaN(disp1)) {
									if (!(max_n >= disp1)) { //  handles initial max_n==NaN too 
										max_n = disp1;
									}
								}
							}
						}
						double disp = disparity[tile];
						if (!Double.isNaN(disp)) {
// 							neib_lev[tile] = (seams[tile] > 0) ? (max_neib_lev - seams[tile] + 1) : 0; // 2-> 1; 1 -> 2; 0-> 0;
							if (seams[tile] > 0) {
								neib_lev[tile] = max_neib_lev - seams[tile] + 1;
								initial_seam[tile] = true; 
							} else {
								neib_lev[tile] = 0;
							}
						}
						if (!Double.isNaN(max_n)) { // got at least 1 neighbor
							if (Double.isNaN(disp)) {
								max_neib[tile] = max_n; // is it needed? Yes, for ordering
								loc_this.add(tile);
							} else { // disparity defined, is it a conflict?
								// is it a conflict?
								if (disp < max_n) {
									double max_diff = disp_adiff + disp_rdiff * Math.max(0.0, max_n);
									if (disp < (max_n - max_diff)) {
										loc_this.add(tile);
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		// Combine lists from multithreaded output to a common one
		int loc_len = 0;
		for (ArrayList<Integer> part_loc: loc_multi) {
			loc_len+= part_loc.size();
		}
		loc_list.clear();
		loc_list.ensureCapacity(loc_len);
		for (ArrayList<Integer> part_loc: loc_multi) {
			loc_list.addAll(part_loc);
		}
		
		while (!loc_list.isEmpty()) {
			// Sort list by decreasing max_neib;
			Collections.sort(loc_list, new Comparator<Integer>() {
				@Override
				public int compare(Integer lhs, Integer rhs) { // descending
					return (max_neib[rhs] > max_neib[lhs]) ? 1 : ((max_neib[rhs] < max_neib[lhs]) ? -1 : 0) ; // lhs.compareTo(rhs);
				}
			});
			// go through list, if still has conflict - replace disparity and neib_lev[], put into lor_list 
			// max_neib_lev
			lor_list.clear();
			lor_list.ensureCapacity(loc_list.size());
			while (!loc_list.isEmpty()) {
				int tile = loc_list.remove(0);
//				if (((tile >= 4028) && (tile <= 4032)) || ((tile >= 4108) && (tile <= 4112))) {
//					System.out.println("buildTileCluster().11: tile="+tile);
//					System.out.println();
//				}
				// find highest neighbor of neib_lev < max_neib_lev
//				double max_n = Double.NaN;
//				int source_neib_level = 0; // maybe find max separately for each neib_level, and assign
				// lower disparity but lower neib_level if both conflict?
				double [] max_n = new double [max_neib_lev];
				Arrays.fill(max_n, Double.NaN);
				for (int dir = 0; dir < 8; dir++) {
					int tile1 = tn.getNeibIndex(tile, dir);
					if (tile1 >= 0) { // && (neib_lev[tile1] >= 0) && (neib_lev[tile1] < max_neib_lev)){ // has defined neighbor
						int nlev1 = neib_lev[tile1];
						double disp1 =  disparity[tile1];
						if (!Double.isNaN(disp1) && (nlev1 >=0) && (nlev1 < max_neib_lev)) {
							if (!Double.isNaN(disp1)) {
								if (!(max_n[nlev1] >= disp1)) { //  handles initial max_n==NaN too 
									max_n[nlev1] = disp1;
								}
							}
						}
					}
				}
				if (Double.isNaN(disparity[tile])) { // add previously undefined
					for (int i = 0; i < max_n.length; i ++)  if (!Double.isNaN(max_n[i])) {
						disparity[tile] = max_n[i];
						neib_lev[tile] =  i + 1; // was -1
						lor_list.add(tile);
						break;
					}
				} else { // old one, find the lowest neighbor conflict
					for (int i = 0; i < max_n.length; i ++)  if (!Double.isNaN(max_n[i]) && (disparity[tile] < max_n[i])) {
						double max_diff = disp_adiff + disp_rdiff * Math.max(0.0, max_n[i]);
						if (disparity[tile] < (max_n[i] - max_diff)) { // it is a conflict
							disparity[tile] = max_n[i];
							neib_lev[tile] =  i + 1; // was -1
							lor_list.add(tile);
							break;
						}
					}					
				}
			} // while (!loc_list.isEmpty()) { finished with loc_list, created lor_list
			loc_list.clear(); // restarting building new list of conflicts
			while (!lor_list.isEmpty()) { // may be already empty
				int tile0 = lor_list.remove(0);
				// look around, for conflicts, add if was not already there (consider using additional array?)
				for (int dir0 = 0; dir0 < 8; dir0++) {
					int tile = tn.getNeibIndex(tile0, dir0);
//					if (((tile >= 4028) && (tile <= 4032)) || ((tile >= 4108) && (tile <= 4112))) {
//						System.out.println("buildTileCluster().12: tile="+tile+", tile0="+tile0);
//						System.out.println();
//					}
					// tries many times as it does not qualify to be added
					if ((tile >= 0) && !loc_list.contains(tile)) { // Do not check+add same tile
						double disp = disparity[tile];
						// See if there is a conflict
						double max_n = Double.NaN;
						for (int dir = 0; dir < 8; dir++) {
							int tile1 = tn.getNeibIndex(tile, dir);
							if (tile1 >= 0) { // has defined neighbor
								double disp1 =  disparity[tile1];
								int nlev1 = neib_lev[tile1];
								if (!Double.isNaN(disp1) && (nlev1 >= 0) && (nlev1 < max_neib_lev)) {
									if (!(max_n >= disp1)) { //  handles initial max_n==NaN too 
										max_n = disp1;
									}
								}
							}
						}
						//
						if (!Double.isNaN(max_n)) { // got at least 1 neighbor
//							if (((tile >= 4028) && (tile <= 4032)) || ((tile >= 4108) && (tile <= 4112))) {
//								System.out.println("buildTileCluster().13: tile="+tile+", tile0="+tile0);
//								System.out.println();
//							}
							if (Double.isNaN(disparity[tile])) {
								max_neib[tile] = max_n; // is it needed? Yes, for ordering
								loc_list.add(tile);
							} else { // disparity defined, is it a conflict?
								// is it a conflict?
								if (disparity[tile] < max_n) {
									double max_diff = disp_adiff + disp_rdiff * Math.max(0.0, max_n);
									if (disparity[tile] < (max_n - max_diff)) {
										loc_list.add(tile);
									}
								}
							}
						}
					}
				}
			} // while (!lor_list.isEmpty()) { // may be already empty
		} // while (!loc_list.isEmpty()) { - no conflicts left, finalize
		final int [] dbg_neib_lev_preorph = (debugLevel > 0)? neib_lev.clone() : null;
		final double [] dbg_disparity1 = (debugLevel > 0)? disparity.clone() : null;
		// if (debugLevel > 0)
		
		// mark selected tiles that conflict with max_neib_lev
		ai.set(0);
		ati.set(0);
		// re-create list of conflicts of defined tiles with neib_lev[] < max_neib_lev
		loc_multi.clear();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			loc_multi.add(new ArrayList<Integer>());
			threads[ithread] = new Thread() {
				public void run() {
					ArrayList<Integer> loc_this = loc_multi.get(ati.getAndIncrement());
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement())
						if ((neib_lev[tile] >= 0) && (neib_lev[tile] < max_neib_lev)) {
						double max_n = Double.NaN;
						for (int dir = 0; dir < 8; dir++) {
							int tile1 = tn.getNeibIndex(tile, dir);
							if ((tile1 >= 0) && (neib_lev[tile1] == max_neib_lev)) { // only conflicts with max_neib_lev
								double disp1 =  disparity[tile1];
								if (!Double.isNaN(disp1)) {
									if (!(max_n >= disp1)) { //  handles initial max_n==NaN too 
										max_n = disp1;
									}
								}
							}
						}
						if (disparity[tile] < max_n) { // works with Double.isNaN(max_n)
							double max_diff = disp_adiff + disp_rdiff * Math.max(0.0, max_n);
							if (disparity[tile] < (max_n - max_diff)) {
								loc_this.add(tile);
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		// Combine lists from multithreaded output to a common one
		loc_list.clear();
		for (ArrayList<Integer> part_loc: loc_multi) {
			loc_list.addAll(part_loc);
		}
		boolean [] discontinued = new boolean [neib_lev.length];
		// Temporarily mark loc_list with max_neib_lev+1 to remove them from averaging
		for (int tile:loc_list) {
			neib_lev[tile] = max_neib_lev+1;
			discontinued[tile] = true;
		}
		
		// there may be some orphans left neib_lev >0 that do not have neighbors with neib_lev one less
		// Recalculate replaced disparities. For now - just averaging, maybe use 5x5 plane best fit?
		// Some of the actual tiles (that conflict with max_neib_lev) are temporarily marked with max_neib_lev+1
		// to prevent them from being averaged
		final double [] wdir = {1.0, 0.7, 1.0, 0.7, 1.0, 0.7, 1.0, 0.7};
		for (int nlev = 1; nlev <= max_neib_lev; nlev++) {
			final int fnlev = nlev;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							// do not modify border tiles that have original disparity (stitched to previous cluster)
							if ((neib_lev[tile] == fnlev) && (disparity[tile] != source_disparity[tile])){
								double swd = 0, sw = 0;
								for (int dir = 0; dir < 8; dir++) {
									int tile1 = tn.getNeibIndex(tile, dir);
									if ((tile1 >= 0) && (neib_lev[tile1] >= 0) && (neib_lev[tile1] < fnlev)) {
										double w = wdir[dir];
										sw += w;
										swd += w * disparity[tile1];
									}
								}
								if (sw > 0.0) {
									disparity[tile] = swd/sw;
								} else {
									neib_lev[tile] = (fnlev < max_neib_lev) ? (fnlev + 1) : -1;
									if (debugLevel > 0) {
										System.out.println("buildTileCluster() removed orphan tile "+tile+", fnlev="+fnlev+ " for cluster "+cluster_list.size());
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		final int [] dbg_neib_lev_predefined = (debugLevel > 0)? neib_lev.clone() : null;
		final AtomicInteger num_removed = new AtomicInteger(0); 		

		// Recreate border tiles by selecting existing one touching last detected conflicts
		// current loc_list is marked with max_neib_lev+1, will change to max_neib_lev.
		// As the new border may overlap old ones, 
		
		
		if (!loc_list.isEmpty()) {
			for (int nlev = max_neib_lev; nlev > 0; nlev--) {
				if (loc_list.isEmpty()) {
					break;
				}
				for (int tile:loc_list) {
//					if (tile == 2608) {
//						System.out.println("buildTileCluster().31: tile="+tile);
//						System.out.println();
//					}
					neib_lev[tile] = nlev;
					// mark new seams:
					for (int layer = 0; layer < disparity_layers.length; layer++) {
						if (disparity_layers[layer][tile] == disparity[tile]) {
							seams_layers[layer][tile] = nlev;
							break;
						}
					}
					new_seam [tile] = true;
				}
				if (nlev == 1) { // just mark, do not create a new list
					break;
				}
				int llen = loc_list.size();
				for (int i = 0; i < llen; i++) {
					int tile = loc_list.remove(0);
					for (int dir = 0; dir < 8; dir++) {
						int tile1 = tn.getNeibIndex(tile, dir);
						if ((tile1 >= 0) &&
								(neib_lev[tile1] >=0 ) &&
								((neib_lev[tile1] < (nlev -1)) || initial_seam[tile1]) 
								) {
							if (!loc_list.contains(tile1)) {
								loc_list.add(tile1);
							}
							initial_seam[tile1] = false;
						}
					}
				}
			}
		}
		
		// Remove selected inner tiles from disparity_layers
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (!Double.isNaN(disparity[tile])){
						if ((neib_lev[tile] == 0) || // is it an inner tile, or
								// an old seam area now resolved
								(       !new_seam[tile] && // resolved but made a seam again
										(seams[tile] > 0) &&
										((max_neib_lev - seams[tile] - neib_lev[tile] + 1) == 0))) {
							for (int layer = 0; layer < disparity_layers.length; layer++) {
								if (disparity_layers[layer][tile] == disparity[tile]) {
									if (tile == dbg_tile ) {
										System.out.println("buildTileCluster() tile="+tile);
									}
									disparity_layers[layer][tile] = Double.NaN;
									seams_layers[layer][tile] = 0;
									num_removed.getAndIncrement();
									break;
								}
							}
					    }
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (num_removed.get() == 0) {
			System.out.println("buildTileCluster() BUG - no tiles removed from disparity_layers[] for cluster "+cluster_list.size());
		}
		
// Marking	discontinued max_neib_lev with max_neib_lev+1 for building meshes they should not be connected to max_neib_lev	
// Maybe there is a more elegant way to do this.		
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (discontinued[tile]){
						neib_lev[tile] = max_neib_lev+1;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debugLevel > 1) {
			String [] dbg_titles = {"Source","Intermediate","Final", "neib_lev0", "neib_lev1", "neib_lev2",
					"seams", "seams_layers_0", "seams_layers_1", "newseam_discontinued",
					"disparity_layers_0", "disparity_layers_1"};
			double [][] dbg_neib_lev = new double [7][tiles];
			
			for (int i = 0; i < tiles; i++) {
				if (dbg_neib_lev_preorph != null) {
					dbg_neib_lev[0][i] = 10*dbg_neib_lev_preorph[i];
				}
				if (dbg_neib_lev_predefined != null) {
					dbg_neib_lev[1][i] = 10*dbg_neib_lev_predefined[i];
				}
				dbg_neib_lev[2][i] = 10*neib_lev[i];
				dbg_neib_lev[3][i] = 10*seams[i];
				dbg_neib_lev[4][i] = 10*seams_layers[0][i];
				dbg_neib_lev[5][i] = 10*seams_layers[1][i];
				dbg_neib_lev[6][i] = 10*(new_seam[i]?1:0) + 20*(discontinued[i]? 1:0);
			}
			
			double [][] dbg_img = {
					source_disparity,
					dbg_disparity1,
					disparity,
					dbg_neib_lev[0],
					dbg_neib_lev[1],
					dbg_neib_lev[2],
					dbg_neib_lev[3],
					dbg_neib_lev[4],
					dbg_neib_lev[5],
					dbg_neib_lev[6],
					disparity_layers[0],
					disparity_layers[1]};
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"source_final_disparity-"+String.format("%02d", cluster_list.size()),
					dbg_titles);
		}		
		
		
		
		
		// Split result into connected clusters
		int [] pnum_clust = new int[1];
		boolean [] sel_tiles = new boolean[tiles];
		for (int tile = 0; tile < tiles; tile++) {
			sel_tiles[tile] = neib_lev[tile] >= 0;
		}
		int [] enum_clusters = tn.enumerateClusters(
				sel_tiles,
				pnum_clust,
				false);
		int num_sub = pnum_clust[0];
		if ((debugLevel > -2) && (num_sub > 1)) {
			System.out.println("buildTileCluster(): splitting textureCluster into "+
					num_sub + " connected ones, current cluster before them is "+cluster_list.size());
		}
// neib_lev[tile]
// disparity
// border_int
		TileCluster [] tileClusters = new TileCluster[pnum_clust[0]]; 
		for (int nsub = 0; nsub < num_sub; nsub++) {
			int []     neib_lev_sub =  (num_sub > 1) ? (new int [tiles]):    neib_lev;
			double  [] disparity_sub = (num_sub > 1) ? (new double [tiles]): disparity;
			if (num_sub > 1) {
				Arrays.fill(neib_lev_sub,  -1);
				Arrays.fill(disparity_sub, Double.NaN);
				int indx = nsub+1;
				for (int tile = 0; tile < tiles; tile++) if (enum_clusters[tile] == indx){
					neib_lev_sub[tile] =  neib_lev[tile];
					disparity_sub[tile] = disparity[tile];
				}
			}
			// find bounds
			AtomicInteger min_y = new AtomicInteger(tilesY);
			AtomicInteger max_y = new AtomicInteger(0);
			AtomicInteger min_x = new AtomicInteger(tilesX);
			AtomicInteger max_x = new AtomicInteger(0);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (neib_lev_sub[tile] >= 0) {
							int tileY = tile / tilesX;
							int tileX = tile % tilesX;
							min_y.getAndAccumulate(tileY, Math::min);
							max_y.getAndAccumulate(tileY, Math::max);
							min_x.getAndAccumulate(tileX, Math::min);
							max_x.getAndAccumulate(tileX, Math::max);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			//		final boolean sky_cluster = blue_sky_below >=0;
			if (is_sky_cluster && (blue_sky_below >= 0)) { // increase bounding box for sky cluster
				min_y.set(0);
				max_y.addAndGet(blue_sky_below);
				if (max_y.get() >= tilesY) {
					max_y.set(tilesY-1);
				}
				min_x.set(0);
				max_x.set(tilesX -1);
			}
			final int width =  max_x.get() - min_x.get() + 1;
			final int height = max_y.get() - min_y.get() + 1;
			final Rectangle bounds = new Rectangle(min_x.get(), min_y.get(), width, height);
			final double  [] disparity_crop =   new double [width * height]; 
			//		final boolean [] border_crop =      new boolean   [disparity_crop.length];
			final int []     border_int_crop =  new int   [disparity_crop.length];
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile_crop = ai.getAndIncrement(); tile_crop < disparity_crop.length; tile_crop = ai.getAndIncrement()) {
							int tileY = tile_crop / width + bounds.y ;
							int tileX = tile_crop % width + bounds.x;
							int tile = tileX + tileY * tilesX;
							disparity_crop[tile_crop] = disparity_sub[tile];
							border_int_crop[tile_crop] = neib_lev_sub[tile];
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			// Create new TileCluster

			tileClusters[nsub] = (new TileCluster(
					bounds,
					cluster_list.size(), // (debug_index? cluster_list.size(): -1),
					null, // border_crop, // will create from border_int_crop
					border_int_crop,     // int []     border_int,           // will replace border? Provide on-the-fly? 
					max_neib_lev,        // int        border_int_max,       // outer border value
					disparity_crop,
					is_sky_cluster));       // boolean is_sky));
			cluster_list.add(tileClusters[nsub]);
		}
		return tileClusters;
	}
	
	/**
	 * Create texture clusters collectively covering the depth map from multi-layer
	 * disparity arrays and optional "blue sky" binary array
	 * @param tilesX               horizontal array dimension 
	 * @param disparity_layers_src [layers][tiles] multi-layer disparity array with NaN
	 *                             for non-existing tiles 
	 * @param blue_sky             per-tile array of the "blue sky" boolean array 
	 * @param blue_sky_layer       layer for the "blue sky" data (should be on a single layer)
	 * @param blue_sky_below       if >=0, specifies how low to extend blue sky backdrop below
	 *                             the lowest defined blue sky tile. In that case bounds are
	 *                             also extended to the top of frame, right and left margins.
	 *                             If (blue_sky_below <0) cluster its bounds are treated the same
	 *                             as regular ones. 
	 * @param max_neib_lev         Maximal level of additional (to defined in disparity_layers)
	 *                             tiles. Currently max_neib_lev=2 adding 2 border layers. So
	 *                             a single isolated tile results in total 5x5 area with level 2,
	 *                             center 3x3 - level 1, and the center tile(defined) - level 0.
	 *                             Other (undefined) tiles are assigned level of -1.  
	 * @param disp_adiffo          absolute disparity difference for connecting tiles in ortho
	 *                             directions (used for initial connections estimations)
	 * @param disp_rdiffo          relative disparity difference (added to absolute being multiplied
	 *                             by the tile disparity)
	 * @param disp_adiffd          absolute disparity difference for diagonal tiles.
	 * @param disp_rdiffd          relative disparity difference for diagonal tiles.
	 * @param disp_fof             >=1.0 - increased inter-tile tolerance for a friend-of-a-friend.
	 *                             In current code just scales calculated (absolute+relative) tolerance
	 *                             for all steps after initial connections.
	 * @param jump_r               "jump" over small gaps when building initial clusters. jump_r == 2 
	 *                             allows jumping to other tiles in 5x5 square around the defined tile,
	 *                             jump_r == 3 - inside 7x7 square. The jump destination should not have
	 *                             any already defined tiles among 8 neighbors (in a 3x3 square). 
	 * @param disp_adiffj          maximal absolute disparity difference for the "jumps".
	 * @param disp_rdiffj          maximal relative disparity difference for the "jumps".       
	 * @param debugLevel           debug level - controls generation of images
	 * @return                     array of consolidated (multiple non-overlapping clusters) clusters,
	 *                             each having the full frame bounds.
	 */
	public static TileCluster [] clusterizeFgBg( //
			final int          tilesX,
			final double [][]  disparity_layers_src, // may have more layers
			final boolean []   blue_sky, // use to expand background by blurring available data?
			final int          blue_sky_layer,
			final int          blue_sky_below,
			final int          max_neib_lev,
			final double       disp_adiffo,
			final double       disp_rdiffo,
			final double       disp_adiffd,
			final double       disp_rdiffd,
			final double       disp_fof,    // enable higher difference (scale) for friend of a friend
			final int          jump_r,
			final double       disp_adiffj,
			final double       disp_rdiffj,
			final int          debugLevel) {
		final int tiles = disparity_layers_src[0].length;
		final int tilesY = tiles/tilesX;
		final int layers = disparity_layers_src.length;
		final int [][] num_neibs_dir = new int[tiles][layers]; // -1 - none, otherwise - bitmask
		// Needed to tell next clusters of the BG seam tiles (shared by >=2 clusters)
		final int [][] seams_layers = new int [layers][tiles];
		final int []   seams =        new int [tiles]; // seams matching disparity[]
		//copy original disparity_layers_src to disparity_layers - they will be modified
		final double [][] disparity_layers = new double [disparity_layers_src.length][];
		for (int i = 0 ; i < disparity_layers.length; i++) {
			disparity_layers[i] = disparity_layers_src[i].clone();
		}
		// maybe ncluster[][] will not be used at all - disparity_layers will be modified to NaN used tiles
		
		// calculate initial num_neibs_dir
		updateSeeds( // and update num_neibs_dir
				num_neibs_dir,    // final int [][]     num_neibs_dir, // [tile][layer]
				null,             // final Rectangle    bounds,    // null - all
				disparity_layers, // final double [][]  disparity_layers, // [layer][tile]should not have same tile disparity on multiple layers
				blue_sky,         // final boolean []   blue_sky, // use to expand background by blurring available data?
				blue_sky_layer,   // final int          blue_sky_layer,
				disp_adiffo,      // final double       disp_adiffo,
				disp_rdiffo,      // final double       disp_rdiffo,
				disp_adiffd,      // final double       disp_adiffd,
				disp_rdiffd,      // final double       disp_rdiffd,
				disp_fof,         // final double       disp_fof,    // enable higher difference (scale) for fried of a friend 
				tilesX,           // final int          tilesX,
				debugLevel);      // final int          debugLevel)
				
		if (debugLevel > -1) { // was > 0
			String [] dbg_titles = {"FG","BG"};
			double [][] dbg_img = new double[layers][tiles];
			for (int i = 0; i < tiles;i++) {
				for (int j = 0; j < dbg_img.length; j++) {
					dbg_img[j][i] = NUM_NEIBS_FROM_BITS[num_neibs_dir[i][j]];
				}
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"num_neibs",
					dbg_titles);
			ShowDoubleFloatArrays.showArrays(
					disparity_layers,
					tilesX,
					tilesY,
					true,
					"disparity_layers",
					dbg_titles);
		}
		
		final ArrayList <TileCluster> cluster_list = new ArrayList<TileCluster>();
		// build all clusters
		int tile_start = 2820; // 0; // change to debug tile to start with the largest one
		while (true) {
			int [] next_seed_tile_layer = getNextSeed(
					disparity_layers, // final double [][] disparity_layers, //
					num_neibs_dir,    // final int [][]    num_neibs_dir, // [tile][layer]
					tile_start,       // final int         tile_start,
					tilesX) ;         //  final int         tilesX)
			if (next_seed_tile_layer == null) {
				break;
			}
			// next_seed_tile_layer is now {tile, layer}
			final boolean is_sky_cluster = (next_seed_tile_layer[1] == blue_sky_layer) &&  blue_sky[next_seed_tile_layer[0]];
			double [] cluster_initial_disparity = buildInitialCluster(
					disparity_layers,        // final double [][]     disparity_layers, // should not have same tile disparity on multiple layers
			        seams_layers,            // final int [][]        seams_layers,
			        seams,                   // final int []          seams,
					next_seed_tile_layer[1], // final int             start_layer,
					next_seed_tile_layer[0], // final int             start_tile,
					blue_sky,                // final boolean []   blue_sky, // use to expand background by blurring available data?
					blue_sky_layer,          // final int          blue_sky_layer,
					disp_adiffo,             // final double          disp_adiffo,
					disp_rdiffo,             // final double          disp_rdiffo,
					disp_adiffd,             // final double          disp_adiffd,
					disp_rdiffd,             // final double          disp_rdiffd,
					disp_fof,                // final double          disp_fof,    // enable higher difference (scale) for friend of a friend
					jump_r,                  // final int             jump_r,
					disp_adiffj,             // final double          disp_adiffj,
					disp_rdiffj,             // final double          disp_rdiffj,
					tilesX,                  // final int             tilesX,
					debugLevel);             // final int             debugLevel)
			
			final double          disp_adiff = disp_fof * disp_adiffd; // should already include disp_fof,
			final double          disp_rdiff = disp_fof * disp_rdiffd; // should already include disp_fof,
			TileCluster [] tileClusters = buildTileCluster(
					// used disparity_layers will be set to Double.NaN
					// make it in a separate method?
					cluster_list,              // final ArrayList <TileCluster> cluster_list,
					is_sky_cluster,            // is_sky_cluster, // this is a blue sky cluster, mark as such and extend bounds
					blue_sky_below,            // final int             blue_sky_below, // >=0 this is a blue sky cluster, mark as such and extend bounds
					disparity_layers,          // final double [][]     disparity_layers, // should not have same tile disparity on multiple layers
					seams_layers,              // final int [][]        seams_layers,
			        seams,                     // final int []          seams,
					cluster_initial_disparity, // final double []       source_disparity, // should not have same tile disparity on multiple layers
					max_neib_lev,              // final int             max_neib_lev,
					disp_adiff,                // final double          disp_adiff, // should already include disp_fof,
					disp_rdiff,                // final double          disp_rdiff,
					tilesX,                    // final int             tilesX,
					debugLevel);               // final int             debugLevel)
//			if (debugLevel > -1000) {
//				return null;
//			}
			for (int nsub = 0; nsub < tileClusters.length; nsub++) {
				updateSeeds( // and update num_neibs_dir
						num_neibs_dir,           // final int [][]     num_neibs_dir, // [tile][layer]
						tileClusters[nsub].getBounds(), // final Rectangle    bounds,    // null - all
						disparity_layers,        // final double [][]  disparity_layers, // [layer][tile]should not have same tile disparity on multiple layers
						blue_sky,                // final boolean []   blue_sky, // use to expand background by blurring available data?
						blue_sky_layer,          // final int          blue_sky_layer,
						disp_adiffo,             // final double       disp_adiffo,
						disp_rdiffo,             // final double       disp_rdiffo,
						disp_adiffd,             // final double       disp_adiffd,
						disp_rdiffd,             // final double       disp_rdiffd,
						disp_fof,                // final double       disp_fof,    // enable higher difference (scale) for fried of a friend 
						tilesX,                  // final int          tilesX,
						debugLevel);             // final int          debugLevel)
			}
			if (debugLevel > 1) {
				String [] dbg_titles = {"FG","BG"};
				double [][] dbg_img = new double[layers][tiles];
				for (int i = 0; i < tiles;i++) {
					for (int j = 0; j < dbg_img.length; j++) {
						dbg_img[j][i] = NUM_NEIBS_FROM_BITS[num_neibs_dir[i][j]];
					}
				}
				ShowDoubleFloatArrays.showArrays(
						dbg_img,
						tilesX,
						tilesY,
						true,
						"num_neibs-"+String.format("%02d", cluster_list.size()),
						dbg_titles);
				ShowDoubleFloatArrays.showArrays(
						disparity_layers,
						tilesX,
						tilesY,
						true,
						"disparity_layers-"+String.format("%02d", cluster_list.size()),
						dbg_titles);
			}
			
			
			
			tile_start = next_seed_tile_layer[0];
		} // while (true) {		
//		int [] tile_stat = new int [tiles];
//		int [] tile_layer = new int [tiles]; // just to know which layer was used for assigned tiles
		// consolidate clusters "good enough", use bounding box intersections, add cluster_gap to grow extra tiles by Gaussian
		// cluster_gap
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
//					Rectangle new_bounds = cluster_list.get(index_other).getBounds(cluster_gap); // cluster_gap should be 2x
					Rectangle new_bounds = cluster_list.get(index_other).getBounds(); // cluster_gap should be 2x
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
		final boolean debug_index = debugLevel > -2; // 0;
		for (int i = 0; i < this_combo; i++) {
			consolidated_clusters[i] = new TileCluster(
					full_tiles,
					(debug_index? 0:-1),
					null,
					null,     // int []     border_int,           // will replace border? Provide on-the-fly? 
					0,        // int        border_int_max,       // outer border value
					null,
					false); // boolean is_sky));
		}
		for (int i = 0; i < comb_clusters.length; i++) {
			consolidated_clusters[comb_clusters[i]].add(cluster_list.get(i));
		}

		if (debugLevel > 0) {
			double [][] dbg_img =         new double[this_combo][tiles];
			double [][] dbg_borders =     new double[this_combo][tiles];
			double [][] dbg_borders_int = new double[this_combo][tiles];
			double [][] dbg_index = null;
			if (debug_index) {
				dbg_index = new double[this_combo][tiles];
			}
			for (int n = 0; n < dbg_img.length; n++) {
				for (int i = 0; i < tiles;i++) {
					dbg_img[n][i] =     consolidated_clusters[n].getDisparity()[i];
					dbg_borders[n][i] = consolidated_clusters[n].getBorder()[i]? 1.0:0.0;
					dbg_borders_int[n][i] = consolidated_clusters[n].getBorderInt()[i];
					if (dbg_index != null) {
						double d = consolidated_clusters[n].getClusterIndex()[i];
						dbg_index[n][i] = (d >=0)? d : Double.NaN;
					}
				}				
			}			
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					"cluster_disparity");
			ShowDoubleFloatArrays.showArrays(
					dbg_borders,
					tilesX,
					tilesY,
					true,
					"cluster_borders");
			ShowDoubleFloatArrays.showArrays(
					dbg_borders_int,
					tilesX,
					tilesY,
					true,
					"cluster_borders_int");
			if (dbg_index != null) {
				ShowDoubleFloatArrays.showArrays(
						dbg_index,
						tilesX,
						tilesY,
						true,
						"cluster_indices");
			}
		}
		return consolidated_clusters;
	}

	public static boolean output3d( // USED in lwir
			CLTParameters                            clt_parameters,
			ColorProcParameters                      colorProcParameters,
			EyesisCorrectionParameters.RGBParameters rgbParameters,
			final QuadCLT                            parameter_scene, // to use for rendering parameters in multi-series sequences
            // if null - use reference scene 
			QuadCLT []                               scenes,
			double [][]                              combo_dsn_final, // null OK, will read file
			final boolean                            updateStatus,
			final int                                debugLevel)
	{
		final boolean batch_mode =        clt_parameters.multiseq_run; // batch_run;
		final boolean gltf_emissive =     clt_parameters.gltf_emissive;
		final boolean use_alpha_blend =   clt_parameters.gltf_alpha_blend;
		final double  tex_disp_adiffo =   clt_parameters.tex_disp_adiffo;       // 0.35; // 0.3;  disparity absolute tolerance to connect in ortho directions 
		final double  tex_disp_rdiffo =   clt_parameters.tex_disp_rdiffo;       // 0.12; // 0.1;  disparity relative tolerance to connect in ortho directions
		final double  tex_disp_adiffd =   clt_parameters.tex_disp_adiffd;       // 0.6;  // 0.4;  disparity absolute tolerance to connect in diagonal directions
		final double  tex_disp_rdiffd =   clt_parameters.tex_disp_rdiffd;       // 0.18; // 0.12; disparity relative tolerance to connect in diagonal directions
		final double  tex_disp_fof =      clt_parameters.tex_disp_fof;          // 1.5;  // Increase tolerance for friend of a friend
		final int     jump_r =            clt_parameters.tex_jump_tiles;        // 2;
		final double  disp_adiffj =       clt_parameters.tex_disp_adiffj;
		final double  disp_rdiffj =       clt_parameters.tex_disp_rdiffj;
		final double  tex_fg_bg =         clt_parameters.tex_fg_bg;             // 0.1;  // Minimal FG/BG disparity difference (NaN bg if difference from FG < this)
		final double  max_disparity_lim = clt_parameters.tex_max_disparity_lim; // 100.0;  // do not allow stray disparities above this
		final double  min_trim_disparity =clt_parameters.tex_min_trim_disparity;//   2.0;  // do not try to trim texture outlines with lower disparities
		final int     max_neib_lev =      clt_parameters.tex_max_neib_lev;      //   2; // 1 - single tiles layer around, 2 - two layers
		final boolean split_textures =    clt_parameters.tex_split_textures;    // (debugLevel > 1000); // false;
		final int     subdiv_tiles =      clt_parameters.tex_subdiv_tiles;      //   4; // subdivide tiles to smaller triangles
		final int     sky_below =         clt_parameters.tex_sky_below;         //  10; // extend sky these tile rows below lowest
		final boolean debug_disp_tri =    clt_parameters.tex_debug_disp_tri;    //   !batch_mode && (debugLevel > 0); // TODO: use clt_parameters
		
		// If there is gap between clusters, add extra row of background tiles
		int     add_bg_tiles =           clt_parameters.tex_add_bg_tiles;           // 0; // 1;
		final boolean save_full_textures=clt_parameters.tex_save_full_textures  || !clt_parameters.tex_split_textures; //true; // false; // true;
		final double  alpha_threshold =  clt_parameters.tex_alpha_threshold;    // 0.5;
		final boolean renormalize =      clt_parameters.tex_renormalize;        // true; // false - use normalizations from previous scenes to keep consistent colors
		final boolean no_alpha =        !clt_parameters.tex_alpha;              // true; // also - use jpeg?
		final int     jpeg_quality =     clt_parameters.tex_jpeg_quality; //95;   // JPEG quality for textures
		
		final boolean showTri =          !batch_mode && (debugLevel > -1) && (clt_parameters.show_triangles);
		final boolean disp_hires_tri =    clt_parameters.tex_disp_hires_tri;
		final int dbg_scale_mesh =        clt_parameters.tex_dbg_scale_mesh; //  4; // <=0 - do not show
		
		
		
		final int sky_layer = 0; // source disparity layer that contains "blue sky" 
		
		final int           ref_index =  scenes.length - 1;
		final QuadCLT       ref_scene =  scenes[ref_index];
		final TileProcessor tp =         ref_scene.getTileProcessor();

		if (ref_scene.image_data == null){
			return false; // not used in lwir
		}
		double infinity_disparity = 	ref_scene.getGeometryCorrection().getDisparityFromZ(clt_parameters.infinityDistance);
		X3dOutput x3dOutput = null;
		WavefrontExport    wfOutput = null; 
		ArrayList<TriMesh> tri_meshes = null; 
		long startStepTime=System.nanoTime();
		final int    tilesX = tp.getTilesX();
		final int    transform_size = tp.getTileSize();
		
		final int width = tilesX * transform_size;
		final int height = tp.getTilesY() * transform_size;
		// get multi-scene disparity map for FG and BG and filter it
		if (combo_dsn_final == null) {
			combo_dsn_final =scenes[ref_index].readDoubleArrayFromModelDirectory(
					"-INTER-INTRA-LMA", // String      suffix,
					0,                  // int         num_slices, // (0 - all)
					null);              // int []      wh);
		}
		boolean [] sky_tiles =  new boolean[combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_BLUE_SKY].length];
		boolean [] sky_invert = new boolean[combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_BLUE_SKY].length];
		for (int i = 0; i < sky_tiles.length; i++) {
			sky_tiles[i] = combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_BLUE_SKY][i] > 0.0;
			sky_invert[i] =  !sky_tiles[i]; // not used
		}
		// re-load , should create quadCLTs[ref_index].dsi
		double [][] dls_fg = {
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_DISP],
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_LMA],
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_STRENGTH]
		};
		
		// currently conditionInitialDS() zeroes disparity for blue_sky. TODO: allow some FG over blue_sky?
		double [][] ds_fg = OpticalFlow.conditionInitialDS(
				true, // boolean        use_conf,       // use configuration parameters, false - use following  
				clt_parameters,      // CLTParameters  clt_parameters,
				dls_fg,                 // double [][]    dls
				scenes[ref_index], // QuadCLT        scene,
				debugLevel);         // int debug_level)
		double [][] dls_bg = {
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_DISP_BG].clone(),
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_LMA_BG].clone(),
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_STRENGTH_BG].clone()
		};
		for (int i = 0; i < sky_tiles.length; i++) if (Double.isNaN(dls_bg[0][i])){
			dls_bg[0][i] = dls_fg[0][i];
			dls_bg[1][i] = dls_fg[1][i];
			dls_bg[2][i] = dls_fg[2][i];
		}
		double [][] ds_bg = OpticalFlow.conditionInitialDS(
				true, // boolean        use_conf,       // use configuration parameters, false - use following  
				clt_parameters,      // CLTParameters  clt_parameters,
				dls_bg,                 // double [][]    dls
				scenes[ref_index], // QuadCLT        scene,
				debugLevel);         // int debug_level)
		double[][] ds_fg_bg = {ds_fg[0], ds_bg[0].clone()}; 
		for (int i = 0; i < sky_tiles.length; i++) {
			if (Math.abs(ds_fg_bg[1][i]-ds_fg_bg[0][i]) < tex_fg_bg) {
				ds_fg_bg[1][i] = Double.NaN;
			}
		}
		// Create data for consolidated textures (multiple texture segments combined in same "passes" 
		TileCluster [] tileClusters = clusterizeFgBg( // wrong result type, not decided
				tilesX,            // final int          tilesX,
				ds_fg_bg,          // final double [][]  disparities, // may have more layers
				sky_tiles,         // final boolean      blue_sky, // use to expand background by blurring available data?
				sky_layer,         // final int          sky_layer,
				sky_below,         // final int          blue_sky_below,
				max_neib_lev,      // final int          max_neib_lev,
				tex_disp_adiffo,   // final double       disp_adiffo,
				tex_disp_rdiffo,   // final double       disp_rdiffo,
				tex_disp_adiffd,   // final double       disp_adiffd,
				tex_disp_rdiffd,   // final double       disp_rdiffd,
				tex_disp_fof,      // final double       disp_fof,    // enable higher difference (scale) for friend of a friend
				jump_r,            // final int          jump_r,
				disp_adiffj,       // final double       disp_adiffj,
				disp_rdiffj,       // final double       disp_rdiffj,
				debugLevel); //1); //  2); // final int          debugLevel)
/*
// Debugging up to here:
//		if (debugLevel > -1000) {
//			return false;
//		}
		
		if (tileClusters == null) {
			System.out.println("Temporary exit after clusterizeFgBg()");
			return false;
		}
*/
		boolean [] scenes_sel = new boolean[scenes.length];
		//		for (int i = scenes.length - 10; i <  scenes.length; i++) { // start with just one (reference) scene
		for (int i = 0; i <  scenes.length; i++) { // start with just one (reference) scene
			scenes_sel[i] = true;
		}
		// If there is gap between clusters, add extra row of background tiles
		while ((add_bg_tiles--) > 0) { // obsolete
			extendClustersBackground(
					tileClusters, // final TileCluster[] tileClusters,
					max_disparity_lim,   // final double   	    max_disparity_lim,
					min_trim_disparity, // final double        min_trim_disparity,
					transform_size,     // final int           transform_size,
					tilesX);            // final int           tilesX)
		}
		
		
		double[][][] faded_textures = getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
				clt_parameters,      // final CLTParameters  clt_parameters,
				colorProcParameters, // ColorProcParameters  colorProcParameters,
				parameter_scene,     // final QuadCLT        parameter_scene, // to use for rendering parameters in multi-series sequences
				// if null - use reference scene
				ref_index,           // final int            ref_index,
				scenes,              // final QuadCLT []     scenes,
				scenes_sel,          // final boolean []     scenes_sel, // null or which scenes to process
				tileClusters,        // final TileCluster [] tileClusters, // disparities, borders, selections for texture passes
				renormalize,         // final boolean        re-normalize,  // false - use normalizations from previous scenes to keep consistent colors
				max_disparity_lim,   // final double         max_disparity_lim,   //  100.0;  // do not allow stray disparities above this
				min_trim_disparity,  // final double         min_trim_disparity,  //  2.0;  // do not try to trim texture outlines with lower disparities
				debugLevel);         // final int            debug_level)
		ImagePlus[] combined_textures = getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
				clt_parameters,      // final CLTParameters  clt_parameters,
				no_alpha,            // final boolean        no_alpha,
				scenes[ref_index],   // QuadCLT              ref_scene,
				parameter_scene,     // final QuadCLT        parameter_scene, // to use for rendering parameters in multi-series sequences
				faded_textures,      // double [][][]        faded_textures,
				tilesX,              // int                  tilesX,
				ref_scene.getTileProcessor().getTilesY(),              // int                  tilesY,
				transform_size,      // int                  transform_size,
				debugLevel);         // final int            debug_level)
		
				
		EyesisCorrectionParameters.CorrectionParameters correctionsParameters = ref_scene.correctionsParameters;
		String x3d_dir = ref_scene.getX3dDirectory();
		boolean use_png = (jpeg_quality <= 0) || !no_alpha; 
		if (save_full_textures) { //  || !split_textures) {
			for (int nslice = 0; nslice < combined_textures.length; nslice++) {
				EyesisCorrections.saveAndShow(
						combined_textures[nslice], // imp_texture_cluster,
						x3d_dir,
						use_png,      // correctionsParameters.png,
						false,        // (nslice < 4), // clt_parameters.show_textures,
						jpeg_quality, // -1, // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						1);           //
			}
		}
		boolean [][] combined_alphas = null;
		if (subdiv_tiles > 0) { // Use subdiv_tiles==1 to keep alpha
			combined_alphas = new boolean [faded_textures.length][faded_textures[0][0].length];
			for (int i = 0; i < faded_textures.length; i++) { // TODO: accelerate with multi
				for (int j = 0; j < faded_textures[i][1].length; j++) {
					combined_alphas[i][j] = faded_textures[i][1][j] > alpha_threshold;
				}
			}
		}
		ImagePlus [] imp_textures = null;
		if (split_textures) {
			// Maybe will switch to combined textures (less files)
			imp_textures = splitCombinedTextures(
					tileClusters,       // TileCluster [] tileClusters, //should have name <timestamp>-*
					transform_size,     // int            transform_size,
					combined_textures); // ImagePlus []   combo_textures )
			for (int i = 0; i < imp_textures.length; i++) if (imp_textures[i] != null) { // should not be
				EyesisCorrections.saveAndShow(
						imp_textures[i], // imp_texture_cluster,
						x3d_dir,
						use_png, // correctionsParameters.png,
						false, // (nslice < 4), // clt_parameters.show_textures,
						jpeg_quality, // -1, // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						1); //
			}
		}
		// create x3d file
		if (clt_parameters.output_x3d) {
			x3dOutput = new X3dOutput(
					clt_parameters,
					ref_scene.correctionsParameters,
					ref_scene.getGeometryCorrection(),
					null);// tp.clt_3d_passes);
		}
		if (clt_parameters.output_obj && (x3d_dir != null)) {
			try {
				wfOutput = new WavefrontExport(
						x3d_dir,
						ref_scene.correctionsParameters.getModelName(ref_scene.getImageName()),
						clt_parameters,
						ref_scene.correctionsParameters,
						ref_scene.getGeometryCorrection(),
						null); // tp.clt_3d_passes);
			} catch (IOException e) {
				System.out.println("Failed to open Wavefront files for writing");
				// TODO Auto-generated catch block
				e.printStackTrace();
				// do nothing, just keep
			}
		}
		if ((clt_parameters.output_x3d || clt_parameters.output_glTF) && (x3d_dir != null)) {
			tri_meshes = new ArrayList<TriMesh>();
		}		
		
		if (x3dOutput != null) { // 09.18.2022 For now - skipping background
			x3dOutput.generateBackground(clt_parameters.infinityDistance <= 0.0); // needs just first (background) scan
		}

		// 09.18.2022 - skipping background generation
		
		int num_clusters = -1;
		for (int nslice=0; nslice < tileClusters.length; nslice++) {
			for (int indx: tileClusters[nslice].getSubIndices()) {
				if (indx > num_clusters) num_clusters= indx;
			}
		}
		num_clusters++; // was cluster with largest index, became number of clusters
		final double [][] dbg_tri_disp = debug_disp_tri? (new double [tileClusters.length][width*height]): null;
		final double [][] dbg_tri_tri =  debug_disp_tri? (new double [tileClusters.length][width*height]): null;
		
		int dbg_scaled_width =  tp.getTilesX() * transform_size * dbg_scale_mesh; 
		int dbg_scaled_height = tp.getTilesY() * transform_size * dbg_scale_mesh;
		boolean debug_alpha = false;
		double [][] dbg_mesh_imgs = null;
		if ((dbg_scale_mesh > 0) && disp_hires_tri) {
			debug_alpha =   true;
			dbg_mesh_imgs = new double[tileClusters.length][dbg_scaled_width * dbg_scaled_height];
			// maybe fill with NaN?
		}
		for (int nslice = 0; nslice < tileClusters.length; nslice++){
			if (dbg_tri_disp != null) {
				Arrays.fill(dbg_tri_disp[nslice], Double.NaN);
			}

			final double [][] dbg_disp_tri_slice = (dbg_tri_tri != null) ?
					((dbg_tri_disp != null)? (new double[][] {dbg_tri_disp[nslice], dbg_tri_tri[nslice]}):
						(new double[][] {dbg_tri_tri[nslice]})	): null; 
			
			if (debugLevel > -2){
				System.out.println("Generating cluster images from texture slice "+nslice);
			}
			int [] indices = tileClusters[nslice].getSubIndices();
			Rectangle [] bounds = tileClusters[nslice].getSubBounds(); // tiles?
			Rectangle    texture_bounds = null;  // if not null - allows trimmed combo textures tiles?
			int dbg_tri_indx = 3; // showing triangles for cluster 3 
			for (int sub_i = 0; sub_i < indices.length; sub_i++) {
				Rectangle roi = bounds[sub_i];
				int cluster_index = indices[sub_i];
				ImagePlus imp_texture_cluster = combined_textures[nslice];				
				boolean [] alpha = (combined_alphas != null) ? combined_alphas[nslice] : null;
				double [] dbg_alpha = debug_alpha? faded_textures[nslice][1] : null;
				if (imp_textures != null) { // wrong? roi is in tiles or pixels?
					//transform_size
					imp_texture_cluster = imp_textures[cluster_index];
					if (combined_alphas != null) {
						int alpha_width =  roi.width *  transform_size;
						int alpha_height = roi.height * transform_size;
						int alpha_x0 =     roi.x * transform_size;
						int alpha_y0 =     roi.y * transform_size;
						alpha = new boolean[alpha_width * alpha_height];
						for (int row = 0; row < alpha_height; row++) {
							System.arraycopy(
									combined_alphas[nslice],
									(alpha_y0 + row) * width + alpha_x0,
									alpha,
									row * alpha_width,
									alpha_width);
						}
						if (dbg_alpha != null) {
							dbg_alpha = new double[alpha_width * alpha_height];
							for (int row = 0; row < alpha_height; row++) {
								System.arraycopy(
										faded_textures[nslice][1],
										(alpha_y0 + row) * width + alpha_x0,
										dbg_alpha,
										row * alpha_width,
										alpha_width);
							}
						}

					}
				}
				if (imp_texture_cluster == null) {
					if (debugLevel > -1){
						System.out.println("Empty cluster #"+cluster_index);
					}
					continue;
				}
//				String texturePath = imp_texture_cluster.getTitle()+".png";
				String texturePath = imp_texture_cluster.getOriginalFileInfo().fileName;
				
				double [] scan_disparity = tileClusters[nslice].getSubDisparity(sub_i); // limited to cluster bounds
				boolean [] scan_selected = tileClusters[nslice].getSubSelected(sub_i); // limited to cluster bounds
				int [] scan_border_int =   tileClusters[nslice].getSubBorderInt(sub_i); // limited to cluster bounds
				int    max_border =        tileClusters[nslice].getBorderIntMax();

				// skipping averaging disparity for a whole cluster (needs strength and does not seem to be useful)
				try {
					if (alpha == null) {
					TriMesh.generateClusterX3d( // old version also generates wavefront obj
							(imp_textures == null), //   boolean         full_texture, // true - full size image, false - bounds only
							0, //   int             subdivide_mesh, // 0,1 - full tiles only, 2 - 2x2 pixels, 4 - 2x2 pixels
							x3dOutput,
							wfOutput,  // output WSavefront if not null
							tri_meshes, // ArrayList<TriMesh> tri_meshes,
							texturePath,
							"shape_id-"+cluster_index, // id
							null, // class
							roi, // scan.getTextureBounds(),
							scan_selected, // scan.getSelected(),
							scan_disparity, // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
							clt_parameters.transform_size,
							tp.getTilesX(), // int             tilesX,
							tp.getTilesY(), // int             tilesY,
							ref_scene.getGeometryCorrection(), // GeometryCorrection geometryCorrection,
							clt_parameters.correct_distortions, // requires backdrop image to be corrected also
							showTri && (cluster_index == dbg_tri_indx), // (scanIndex < next_pass + 1) && clt_parameters.show_triangles,
							// FIXME: make a separate parameter:
							infinity_disparity, //  0.25 * clt_parameters.bgnd_range,  // 0.3
							clt_parameters.grow_disp_max, // other_range, // 2.0 'other_range - difference from the specified (*_CM)
							clt_parameters.maxDispTriangle,
						    clt_parameters.maxZtoXY,      // double          maxZtoXY,       // 10.0. <=0 - do not use
						    clt_parameters.maxZ,
						    clt_parameters.limitZ,
						    dbg_disp_tri_slice, //   double [][]     dbg_disp_tri_slice,
							debugLevel + 1); //   int             debug_level) > 0
					} else {
						TriMesh.generateClusterX3d( // new version with small triangles for alpha also generates wavefront obj
								(imp_textures == null), // boolean         full_texture, // true - full size image, false - bounds only
								subdiv_tiles,           // int             subdivide_mesh, // 0,1 - full tiles only, 2 - 2x2 pixels, 4 - 2x2 pixels
								alpha,                  // boolean []      alpha,     // boolean alpha - true - opaque, false - transparent. Full/bounds
								dbg_alpha,              // double []       dalpha,             // before boolean
								x3dOutput,
								wfOutput,               // output WSavefront if not null
								tri_meshes,             // ArrayList<TriMesh> tri_meshes,
								texturePath,
								"shape_id-"+cluster_index, // id
								null,                   // class
								roi,                    // scan.getTextureBounds(),
								texture_bounds,         // Rectangle       texture_bounds,     // if not null - allows trimmed combo textures
								scan_selected,          // scan.getSelected(),
								scan_disparity,         // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
								scan_border_int,        //   int     []      border_int,
								max_border, // int             max_border,
								clt_parameters.transform_size,
								tp.getTilesX(),         // int             tilesX,
								tp.getTilesY(),         // int             tilesY,
								ref_scene.getGeometryCorrection(), // GeometryCorrection geometryCorrection,
								clt_parameters.correct_distortions, // requires backdrop image to be corrected also
								dbg_mesh_imgs[nslice],  //   double []       tri_img,   //
								dbg_scaled_width,       // int             tri_img_width,
//								showTri,                // (scanIndex < next_pass + 1) && clt_parameters.show_triangles,
								// FIXME: make a separate parameter:
								infinity_disparity,            //  0.25 * clt_parameters.bgnd_range,  // 0.3
								clt_parameters.grow_disp_max,  // other_range, // 2.0 'other_range - difference from the specified (*_CM)
								clt_parameters.maxDispTriangle,
							    clt_parameters.maxZtoXY,       // double          maxZtoXY,       // 10.0. <=0 - do not use
							    clt_parameters.maxZ,
							    clt_parameters.limitZ,
//							    dbg_disp_tri_slice,           //   double [][]     dbg_disp_tri_slice,
								debugLevel + 1,               //   int             debug_level) > 0
								clt_parameters.tex_dbg_plot_center,              //   boolean         dbg_plot_center, //  = true;
								clt_parameters.tex_dbg_line_color,               //   double          dbg_line_color, //  =  1.0;
								clt_parameters.tex_dbg_center_color);    //  double          dbg_center_color// = 3.0;
								
					}
					//dbg_disp_tri_slice
				} catch (IOException e) {
					e.printStackTrace();
					return false;
				}
			}
			// if (imp_textures[nslice] != null)
		} // for (int nslice = 0; nslice < tileClusters.length; nslice++){
		
		if (dbg_mesh_imgs != null) {
			ShowDoubleFloatArrays.showArrays(
					dbg_mesh_imgs,
					dbg_scaled_width,
					dbg_scaled_height,
					true,
					ref_scene.getImageName()+"-tri-meshes");
		}
		
		if (dbg_tri_disp != null) {
			ShowDoubleFloatArrays.showArrays(
					dbg_tri_disp,
					width,
					height,
					true,
					ref_scene.getImageName()+"-mesh_disparities");
		}
		if (dbg_tri_tri != null) {
			ShowDoubleFloatArrays.showArrays(
					dbg_tri_tri,
					width,
					height,
					true,
					ref_scene.getImageName()+"-mesh_triangles");
		}
		if ((dbg_tri_disp != null) && (dbg_tri_tri != null)) {
			double [][] dbg_tri = new double [2*dbg_tri_tri.length][];
			String [] dbg_titles = new String [dbg_tri.length];
			for (int i = 0; i < dbg_tri_tri.length; i++) {
				dbg_titles[2*i + 0] = "disparity-"+String.format("%02d", i);
				dbg_titles[2*i + 1] = "triangles-"+String.format("%02d", i);
				dbg_tri[2*i + 0] = dbg_tri_disp[i];
				dbg_tri[2*i + 1] = dbg_tri_tri[i];
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_tri,
					width,
					height,
					true,
					ref_scene.getImageName()+"-mesh_disparity_triangles",
					dbg_titles);
		}
		boolean exit_now = (debugLevel > 1000);
		if (exit_now) {
			return false;
		}
		
		
		
		if ((x3d_dir != null) && (x3dOutput != null)){
			x3dOutput.generateX3D(x3d_dir+Prefs.getFileSeparator() + ref_scene.correctionsParameters.getModelName(ref_scene.getImageName())+".x3d");
		}
		
		if (wfOutput != null){
			wfOutput.close();
			System.out.println("Wavefront object file saved to "+wfOutput.obj_path);
			System.out.println("Wavefront material file saved to "+wfOutput.mtl_path);
		}
		if (clt_parameters.output_glTF && (tri_meshes != null)) {
			try {
				GlTfExport.glTFExport(
						x3d_dir,         // String x3d_dir,
						ref_scene.correctionsParameters.getModelName(ref_scene.getImageName()), // String model_name,
						tri_meshes,      // ArrayList<TriMesh> tri_meshes,
						gltf_emissive,   // boolean gltf_emissive,
						use_alpha_blend, // boolean            use_alpha_blend,
						1);
			} catch (JSONException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} // int debugLevel
		}
		
		// Save KML and ratings files if they do not exist (so not to overwrite edited ones), make them world-writable
		  if (!correctionsParameters.use_set_dirs) {
			  System.out.println("**** output3d(): likely a bug (not copied?), temporary fix ***");
			  correctionsParameters.use_set_dirs = true; 
		  }

		ref_scene.writeKml        (null, debugLevel);
		ref_scene.writeRatingFile (debugLevel);
		Runtime.getRuntime().gc();
		System.out.println("output3d(): generating 3d output files  finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-startStepTime),3)+" sec, --- Free memory25="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		return true;
	}

	public static double [][] getComboTexture (
			final double [][][] sensor_texture) {
		final int num_slices =  sensor_texture.length;
		final int num_sensors = sensor_texture[0].length;
		final int img_size =    sensor_texture[0][0].length;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final double [][] combo_texture = new double [num_slices][img_size];
		final AtomicInteger ai = new AtomicInteger(0);
			for (int nslice = 0; nslice < num_slices; nslice++) {
				final int fnslice = nslice;
				Arrays.fill(combo_texture[fnslice], Double.NaN);
				ai.set(0);
				// calculate total number of connections (w/o fof) by combining opposite directions
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int indx = ai.getAndIncrement(); indx < img_size; indx = ai.getAndIncrement()) if (!Double.isNaN(sensor_texture[fnslice][0][indx])) {
								double d = 0;
								for (int nsens = 0; nsens < num_sensors; nsens++) {
									d += sensor_texture[fnslice][nsens][indx];
								}
								combo_texture[fnslice][indx] = d/num_sensors;
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}
		return combo_texture;
	}
	
	public static double getMaxDisparity (
			double [][] slice_disparities,
			double   	max_disparity_lim) {
		final int num_slices = slice_disparities.length;
		final int tiles = slice_disparities[0].length;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ati = new AtomicInteger(0);
		final double [] th_max_disp = new double [threads.length];
		final int tile_slices = tiles*num_slices;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int ti = ati.getAndIncrement();
					for (int indx = ai.getAndIncrement(); indx < tile_slices; indx = ai.getAndIncrement()) {
						int nslice = indx / tiles;
						int tile =   indx % tiles;
						if (!Double.isNaN(slice_disparities[nslice][tile])) {
							if ((slice_disparities[nslice][tile] > th_max_disp[ti]) &&
									(slice_disparities[nslice][tile] < max_disparity_lim)){
								th_max_disp[ti] = slice_disparities[nslice][tile];
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		double disparity_max = 0.0;
		for (int i = 0; i < th_max_disp.length; i++) {
			disparity_max = Math.max(disparity_max, th_max_disp[i]);
		}
		return disparity_max;
	}
	
	/**
	 * Extend background edges of clusters - those that have FG tiles in front of them
	 * There should be gaps between clusters to expand
	 */
	public static void extendClustersBackground(
				final TileCluster[] tileClusters,
				final double   	    max_disparity_lim,
				final double        min_trim_disparity,
				final int           transform_size,
				final int           tilesX){
		final double weight_diag = 0.7;
		final int num_slices = tileClusters.length;
		final double [][] slice_disparities = new double [num_slices][];
	    final int [][]    slice_border_int = new int [num_slices][];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			slice_disparities[nslice] = tileClusters[nslice].getDisparity();  // disparity in the reference view tiles (Double.NaN - invalid)
	        slice_border_int[nslice] =  tileClusters[nslice].getBorderInt(); 
		}
//	    int border_int_max = (num_slices>0) ?tileClusters[0].getBorderIntMax(): 0;
		
		final int tiles = slice_disparities[0].length;
		final int tilesY = tiles / tilesX;
		final boolean [][][] fg_has_bg = get_fg_has_bg_any( // {is_fg, has_bg, has_tile}
				slice_disparities,  // final double [][] slice_disparities,
				slice_border_int,   // final int    [][] slice_border_int,
//				border_int_max,     // final int         border_int_max,
				max_disparity_lim,  // final double   	max_disparity_lim,
				min_trim_disparity, // final double min_trim_disparity,
				transform_size,     // final int   transform_size,
				tilesX);            // final int   tilesX)
		
		showDebugDisparities( // nop if dbg_prefix== null
				slice_disparities, // final double [][] slice_disparities,
				tilesX,            // final int   tilesX,
				"before-extending-disparities");       // String      prefix);
		
		showDebugFgBg( // nop if dbg_prefix== null
				fg_has_bg,         // boolean [][][] fg_has_bg,
				tilesX,            // final int   tilesX,
				"before-extending-fgbg");       // String      prefix);
		
		
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		// find all border tiles that are not FG
		final boolean [][] seed_tiles = new boolean[num_slices][];
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < num_slices; nslice = ai.getAndIncrement()) {
						seed_tiles[nslice] = tn.getEdgeSelection(
								2,          // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
								fg_has_bg[2][nslice], // boolean [] tiles,
								null);     // boolean [] prohibit);
						for (int i = 0; i < seed_tiles[nslice].length; i++) {
							seed_tiles[nslice][i] &= !fg_has_bg[0][nslice][i]; // keep only those that are not FG 
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		showDebugBoolean(
				seed_tiles, // boolean [][] selected,
				tilesX, // final int   tilesX,
				"seed_tiles"); // String      prefix)
		
		final double [][] extended_slice_disparities = new double [num_slices][];
		final boolean [][] new_tiles = new boolean[num_slices][tiles];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (seed_tiles[fnslice][tile]) {
							for (int dir = 0; dir < 8; dir++) {
								int tile1 = tn.getNeibIndex(tile, dir);
								if (tile1 >= 0) {
									if (Double.isNaN(slice_disparities[fnslice][tile1])) {
										new_tiles[fnslice][tile1] = true; // may be concurrent, but that's OK
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);

			extended_slice_disparities[fnslice] = slice_disparities[fnslice].clone();
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (new_tiles[fnslice][tile]) {
							double sw = 0, swd = 0.0;
							for (int dir = 0; dir < 8; dir++) {
								int tile1 = tn.getNeibIndex(tile, dir);
								if ((tile1 >= 0) && !Double.isNaN(slice_disparities[fnslice][tile1])) { // only average old tiles
									double w = ((dir & 1) != 0)? weight_diag : 1.0;
									sw += w;
									swd += slice_disparities[fnslice][tile1] * w;
								}
							}
							extended_slice_disparities[fnslice][tile] = swd/sw; // weighted average, sw should not be 0 here
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		showDebugBoolean(
				new_tiles, // boolean [][] selected,
				tilesX, // final int   tilesX,
				"new_tiles"); // String      prefix)

		showDebugDisparities( // nop if dbg_prefix== null
				extended_slice_disparities, // final double [][] slice_disparities,
				tilesX,            // final int   tilesX,
				"all-extended-disparities");       // String      prefix);
		
		final boolean [][][] fg_has_bg_ext = get_fg_has_bg_any( // {is_fg, has_bg, has_tile}
				extended_slice_disparities,  // final double [][] slice_disparities,
	            slice_border_int,            // final int    [][] slice_border_int,
//	            border_int_max,              // final int         border_int_max,
				max_disparity_lim,           // final double   	max_disparity_lim,
				min_trim_disparity,          // final double min_trim_disparity,
				transform_size,              // final int   transform_size,
				tilesX);                     // final int   tilesX)
		
		showDebugFgBg( // nop if dbg_prefix== null
				fg_has_bg_ext,         // boolean [][][] fg_has_bg,
				tilesX,            // final int   tilesX,
				"all-extended-fgbg");       // String      prefix);

		
		
		// Now find which tiles are FG, but only compare with old tiles
		// need to modify get_fg_has_bg_any() and make sure that if two disparities are equal, both tiles are FG, DONE
		// but none is count as BG for the other.
		// Now update TileCluster[] tileClusters to include new tiles 
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (new_tiles[fnslice][tile]) {
							if (fg_has_bg_ext[0][fnslice][tile]) { // new is foreground - remove
								extended_slice_disparities[fnslice][tile] = Double.NaN;
								new_tiles[fnslice][tile] = false;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		
		showDebugDisparities( // nop if dbg_prefix== null
				extended_slice_disparities, // final double [][] slice_disparities,
				tilesX,            // final int   tilesX,
				"trimmed-extended-disparities");       // String      prefix);
		
		// update tileClusters[fnslice]

		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < num_slices; nslice = ai.getAndIncrement()) {
						tileClusters[nslice].setDisparity(extended_slice_disparities[nslice]);
						tileClusters[nslice].increaseBounds(); // increase bounds by 1 tile each side
						tileClusters[nslice].resetClusterIndex();
						// will rebuild next time tileClusters[nslice].getClusterIndex() will be called 
 //						tileClusters[nslice].getClusterIndex(); 
						
						/*
						Rectangle [] sub_bounds = tileClusters[nslice].getSubBounds();
						for (int nsub = 0; nsub < sub_bounds.length; nsub++) {
							double [] cluster_disparity = TileNeibs.getDoubleWindow(
									sub_bounds[nsub],                   // Rectangle window,
									extended_slice_disparities[nslice], // double [] data,
									tilesX);                            // int data_width)
							tileClusters[nslice].setSubDisparity(
									nsub,               // int indx,
									cluster_disparity); // double [] sub_disparity)
						}
						*/
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return;
	}
	
	/**
	 * Determine which tiles are foreground ones (not obscured by others) and which tiles have other ones
	 * (not transparent) behind them. Only foreground tiles that have background can be trimmed.
	 * Now trying:
	 *   potentially foreground - neb_lev == 0,1 
	 *   can obscure (deny FG) -  neb_lev == 0,1
	 *   are opaque (can be considered as BG behind) - neb_lev == 1 // 0
	 *    
	 * @param slice_disparities   per-slice, per-tile disparities
	 * @param slice_border_int    border tile status: <0 - undefined, 0 - defined disparities, >0 - border_lev 
	 * @param max_disparity_lim   do not consider tiles with higher disparities (stray ones?)
	 * @param min_trim_disparity  do not consider FG with lower disparities (anyway can not obscure much)
	 * @param transform_size      CLT transform size - now 8
	 * @param tilesX              number of tiles in a full row (80)
	 * @return                    boolean [3][tiles]. [0] which tiles are FG, [1] - which tiles have BG [2] - has_tile 
	 */
	public static boolean [][][] get_fg_has_bg_any(
			final double [][] slice_disparities,
			final int    [][] slice_border_int, // not extended
			final double      max_disparity_lim,
			final double      min_trim_disparity,
			final int         transform_size,
			final int         tilesX) {
		// removing 2 last layers - has_tile and border - replace with border_int
//		final double [][] slice_disparities_real = null; // not extended
		final int num_slices = slice_disparities.length;
		final int tiles = slice_disparities[0].length;
		final int tilesY = tiles / tilesX;
		final boolean [][] is_fg_tile =   new boolean [num_slices][tiles];
		final boolean [][] has_bg_tile =  new boolean [num_slices][tiles];
		final boolean [][] has_tile =     new boolean [num_slices][tiles];
//		final boolean [][] border_tiles = new boolean [num_slices][tiles];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final double disparity_max = getMaxDisparity (
				slice_disparities,    // double [][] slice_disparities,
				max_disparity_lim);   // double   	max_disparity_lim)
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		final int max_bg_lev =       1; // 0; // all bg tiles with slice_border_int> max_bg_lev are considered semi-transparent
		final int max_fg_lev =       1; // maximal border_level that can become foreground
		final int max_obscuring_lev =1; // maximal border level that can obscure other tiles
		final int max_wbg_keep = 1;     // keep tiles that have opaque BG ones with level not greater than 
	
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice; 
			ai.set(0);
			// calculate total number of connections (w/o fof) by combining opposite directions
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (slice_border_int[fnslice][tile] >= 0) {
							// }(has_tile[fnslice][tile]){
							// may be a FG tile that needs trimming (not considering yet tiles that both can be obscured and obscure).
							if ((fnslice == -6) && (tile==2333)) {
								System.out.println("fnslice="+fnslice+", tile="+tile);
								System.out.println("fnslice="+fnslice+", tile="+tile);
							}
							if (slice_disparities[fnslice][tile] > min_trim_disparity)	{						
								is_fg_tile[fnslice][tile] =  slice_border_int[fnslice][tile] <= max_fg_lev;// already tested for >=0// true;
								for (int ns = 0; ns < num_slices; ns++)
									if (    (ns != fnslice) &&
											(slice_border_int[ns][tile] <= max_bg_lev) && // these tiles may be semi-transparent - do not count them
											(slice_border_int[ns][tile] >= 0) &&
											(slice_disparities[ns][tile] < slice_disparities[fnslice][tile])) {
										has_bg_tile[fnslice][tile] = true;
										break;
									}
							}
							search_obscuring:
							{
								int tile_range = (int) Math.ceil((disparity_max - slice_disparities[fnslice][tile])/Math.sqrt(2)/transform_size);
								for (int dty = -tile_range; dty <= tile_range; dty++) {
									for (int dtx = -tile_range; dtx <= tile_range; dtx++) {
										int tile1 = tn.getNeibIndex(tile, dtx, dty);
										if (tile1 >= 0) {
											double dd = (Math.sqrt(dty*dty + dtx*dtx) + 0.0)*transform_size + slice_disparities[fnslice][tile]; // is it correct?
											for (int ns = 0; ns < num_slices; ns++) if ((ns != fnslice) || (dty !=0 ) || (dtx != 0)) {
												if ((slice_disparities[ns][tile1] > dd) &&
														(slice_border_int[fnslice][tile] <= max_obscuring_lev)){ // ignore transparent tiles
													is_fg_tile[fnslice][tile] = false;
													break search_obscuring;
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
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (slice_border_int[fnslice][tile] >= 0){
							has_tile[fnslice][tile] = !has_bg_tile[fnslice][tile] || (slice_border_int[fnslice][tile] <=max_wbg_keep);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return new boolean [][][] {is_fg_tile, has_bg_tile, has_tile}; // , border_tiles};
	}	
	
	public static boolean [][][] getTileBooleans(
			final double [][] slice_disparities,
			final int    [][] slice_border_int, // not extended
			final double      max_disparity_lim,
			final double      min_trim_disparity,
			final int         max_neib_lev,
			final int         transform_size,
			final int         tilesX) {
		final int num_slices = slice_disparities.length;
		final int tiles = slice_disparities[0].length;
		final int tilesY = tiles / tilesX;
		final boolean [][] is_fg =         new boolean [num_slices][tiles]; // any strength 
		final boolean [][] is_fg_weak =    new boolean [num_slices][tiles];
		final boolean [][] is_fg_strong =  new boolean [num_slices][tiles];
		final boolean [][] has_bg_weak =   new boolean [num_slices][tiles];
		final boolean [][] has_bg_strong = new boolean [num_slices][tiles];
		final boolean [][] has_tile =      new boolean [num_slices][tiles];
		final boolean [][] stitch_tile =   new boolean [num_slices][tiles];
		final boolean [][] stitched_tile = new boolean [num_slices][tiles]; // stitch for other
		
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final double disparity_max = getMaxDisparity (
				slice_disparities,    // double [][] slice_disparities,
				max_disparity_lim);   // double   	max_disparity_lim)
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		final int max_bg_lev =         max_neib_lev + 1; // 0; // all bg tiles with slice_border_int> max_bg_lev are considered semi-transparent
		final int max_bg_lev_strong =  max_neib_lev - 1; // 0; // all bg tiles with slice_border_int> max_bg_lev are considered semi-transparent
		final int max_fg_lev =         1; // maximal border_level that can become foreground
		final int max_fg_lev_strong =  0; // maximal border_level that can become foreground
		final int max_obscuring_lev =  1; // maximal border level that can obscure other tiles
		final int max_wbg_keep =       1; // keep tiles that have opaque BG ones with level not greater than ??
		final boolean [] stitch_this =   new boolean [tiles];
//		final int dbg_tile = 3948; // 28+49*80;
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							stitch_this[tile] = slice_border_int[fnslice][tile] > max_neib_lev;
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			for (int nlev = max_neib_lev - 1; nlev > 0; nlev--) {
				int fnlev = nlev;
				if (nlev < (max_neib_lev - 1)) {
					System.arraycopy(stitch_tile[fnslice], 0, stitch_this, 0, tiles);
				}
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
								if (stitch_this[tile]) {
									stitch_tile[fnslice][tile] = true;
								} else {
									if (slice_border_int[fnslice][tile] == fnlev) {
										for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
											int tile1 = tn.getNeibIndex(tile, dir);
											if ((tile1 >= 0) && stitch_this[tile1]) {
												stitch_tile[fnslice][tile] = true;
												break;
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
		}
		
		// mark stitched to other
		
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
//							if (tile == dbg_tile ) {
//								System.out.println("getTileBooleans().2 tile="+tile);
//							}
							if (stitch_tile[fnslice][tile]) {
								for (int ns = 0; ns < num_slices; ns++)
									if ((ns != fnslice) &&
											!stitch_tile[ns][tile] && 
											(slice_disparities[ns][tile] == slice_disparities[fnslice][tile])) {
										stitched_tile[ns][tile] = true;
									}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		

		
		
		// restart outer loop to be available for all slices later
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			// calculate total number of connections (w/o fof) by combining opposite directions
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (slice_border_int[fnslice][tile] >= 0) {
							// may be a FG tile that needs trimming (not considering yet tiles that both can be obscured and obscure).
							if ((fnslice == -6) && (tile==2333)) {
								System.out.println("fnslice="+fnslice+", tile="+tile);
								System.out.println("fnslice="+fnslice+", tile="+tile);
							}
							if (slice_disparities[fnslice][tile] > min_trim_disparity)	{
								is_fg       [fnslice][tile] = true;  
								is_fg_weak  [fnslice][tile] =  (slice_border_int[fnslice][tile] <= max_fg_lev) || stitch_tile[fnslice][tile];
								is_fg_strong[fnslice][tile] =  (slice_border_int[fnslice][tile] <= max_fg_lev_strong) || stitch_tile[fnslice][tile];
								for (int ns = 0; ns < num_slices; ns++)
									if (    (ns != fnslice) &&
											((slice_border_int[ns][tile] <= max_bg_lev) || stitch_tile[ns][tile] )&& // these tiles may be semi-transparent - do not count them
											(slice_border_int[ns][tile] >= 0) &&
											(slice_disparities[ns][tile] < slice_disparities[fnslice][tile])) {
										has_bg_weak[fnslice][tile] = true;
										if ((slice_border_int[ns][tile] <= max_bg_lev_strong)  || stitch_tile[ns][tile]) {
											has_bg_strong[fnslice][tile] = true;
											break;
										}
									}
							}
//							if (is_fg_weak  [fnslice][tile]) {
							if (is_fg  [fnslice][tile]) { // check for weakest
								search_obscuring:
								{
									int tile_range = (int) Math.ceil((disparity_max - slice_disparities[fnslice][tile])/Math.sqrt(2)/transform_size);
									for (int dty = -tile_range; dty <= tile_range; dty++) {
										for (int dtx = -tile_range; dtx <= tile_range; dtx++) {
											int tile1 = tn.getNeibIndex(tile, dtx, dty);
											if (tile1 >= 0) {
												double dd = (Math.sqrt(dty*dty + dtx*dtx) + 0.0)*transform_size + slice_disparities[fnslice][tile]; // is it correct?
												for (int ns = 0; ns < num_slices; ns++) if ((ns != fnslice) || (dty !=0 ) || (dtx != 0)) {
													if ((slice_disparities[ns][tile1] > dd) &&
															((slice_border_int[ns][tile] <= max_obscuring_lev) || stitch_tile[ns][tile] )){ // ignore transparent tiles
														is_fg       [fnslice][tile] = false;
														is_fg_weak  [fnslice][tile] = false;
														is_fg_strong[fnslice][tile] = false;
														break search_obscuring;
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
		// if stitch and FG - make others with the same disparity FG strong
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement())
							if (stitch_tile[fnslice][tile] && is_fg[fnslice][tile]){ 
							for (int ns = 0; ns < num_slices; ns++)
								if ((ns != fnslice) && (slice_disparities[fnslice][tile] == slice_disparities[ns][tile])){
									is_fg       [ns][tile] = true;
									is_fg_weak  [ns][tile] = true;
									is_fg_strong[ns][tile] = true;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		

		// Trim keep after processing stitch
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (slice_border_int[fnslice][tile] >= 0){
							has_tile[fnslice][tile] =
									// Even no-bg should not save peripheral (outside of weak) FG tiles
									(!has_bg_strong[fnslice][tile] && (!is_fg[fnslice][tile] || is_fg_weak  [fnslice][tile])) ||
									(slice_border_int[fnslice][tile] <= max_wbg_keep) ||
									stitch_tile[fnslice][tile] ||
									stitched_tile[fnslice][tile] ||
									is_fg_strong [fnslice][tile]; // other's stitch area and fg makes this strong fg
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		
		// Remove outer stitch tiles (max_neib_lev +1) if there is other (non-stitch) with the same disparity
		// (level==l = does it need to be tested) 
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
//							if (tile == dbg_tile ) {
//								System.out.println("getTileBooleans().6 tile="+tile);
//							}
							if (stitch_tile[fnslice][tile] && (slice_border_int[fnslice][tile] == (max_neib_lev +1))){
								for (int ns = 0; ns < num_slices; ns++)
									if (!stitch_tile[ns][tile] &&
											(slice_disparities[ns][tile] == slice_disparities[fnslice][tile])) {
										has_tile[fnslice][tile] = false;
										break;	
									}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		

		boolean[][][] rslt = new boolean [TILE_BOOLEANS][][];
		rslt[TILE_IS_FG_WEAK] =     is_fg_weak;
		rslt[TILE_IS_FG_STRONG] =   is_fg_strong;
		rslt[TILE_HAS_BG_WEAK] =    has_bg_weak;
		rslt[TILE_HAS_BG_STRONG] =  has_bg_strong;
		rslt[TILE_KEEP] =           has_tile;
		rslt[TILE_STITCH] =         stitch_tile;
		rslt[TILE_STITCHED] =       stitched_tile;
		return rslt;
	}
	
	public static double [][][][] getPixelOffsets(
		final TpTask[][][]   tp_tasks_ref, //
		final boolean [][][] tile_booleans, // to filter?
		final int         tilesX)
	{
		final int num_slices = tile_booleans[0].length;
		final int num_tiles = tile_booleans[0][0].length;
		double [][][][] pix_offsets = new double [num_slices][num_tiles][][];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < num_slices; nslice = ai.getAndIncrement()) {
						if ((tp_tasks_ref[nslice]!= null) && (tp_tasks_ref[nslice].length>0) &&  (tp_tasks_ref[nslice][0]!= null)) {
							for (int ntile = 0; ntile <  tp_tasks_ref[nslice][0].length; ntile++) {
								TpTask task = tp_tasks_ref[nslice][0][ntile];
								int tile = task.getTileX()+task.getTileY()*tilesX;
								pix_offsets[nslice][tile] = task.getDoubleXY();
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);

		return pix_offsets;
	}
	
	
	/**
	 * Select pixels between weak tiles and strong tiles for both has_bg (edge where
	 * triangular mesh will end) and is_fg tiles extending 4 pixels over weak foreground
	 * tiles from strong FG tiles
	 * @param weak_tiles     [tilesX*tilesY] selected weak tiles (include strong tiles)
	 * @param strong_tiles   [tilesX*tilesY] selected stgrong tiles (should be true for weak tiles)
	 * @param grow_tiles     grow strong into weak. == transform_size (8) - grow half tile
	 * @param transform_size CLT transform size (==8)
	 * @param tilesX         number of tiles in a row    
	 * @return               [(tilesX*transform_size) * (tilesY*transform_size)] array that includes
	 *                       pixels corresponding to strong tiles and extended by transform_size/2 over
	 *                       weak tiles.
	 */
	
	public static boolean [][] halfStrong( // select pixels between weak and strong 
			final boolean [][]  weak_tiles,
			final boolean [][]  strong_tiles,
			final int           grow_tiles,
			final int           transform_size,
			final int           tilesX) {
		final int num_slices = weak_tiles.length;
		final int tiles = weak_tiles[0].length;
		final int tilesY = tiles / tilesX;
		final int width =  tilesX * transform_size;
		final int height = tilesY * transform_size;
		final int pixels = width * height;
		final boolean [][] half_strong_pix = new boolean [num_slices][pixels];
		final boolean [] prohibit = new boolean [pixels];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs pn =     new TileNeibs(width, height);
		final boolean dbg = width < 0; // never
		final double [][] dbg_img = dbg? new double [3*num_slices] [pixels]:null;
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			Arrays.fill(prohibit, true);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (weak_tiles[fnslice][tile]) {
							boolean is_strong = strong_tiles[fnslice][tile];
							int tileX = tile % tilesX;
							int tileY = tile / tilesX;
							int indx0 = (tileY * width + tileX) * transform_size;
							for (int y = 0; y < transform_size; y++) {
								int indx1 = indx0 + y * width;  
								for (int x = 0; x < transform_size; x++) {
									int indx2 = indx1 + x;
									half_strong_pix[fnslice][indx2] = is_strong;
									prohibit[indx2] = false;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (dbg_img != null) {
				for (int i = 0; i < pixels; i++) {
					dbg_img[3*fnslice+0][i] = prohibit[i]? 0.0 : 1.0;
					dbg_img[3*fnslice+1][i] = half_strong_pix[fnslice][i]? 1.0 : 0.0;
				}
			}
			
			pn.growSelection( // will fill half of "weak" has_bg_tiles - it is where triangular mesh will reach
					grow_tiles,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					half_strong_pix[fnslice],
					prohibit);
			if (dbg_img != null) {
				for (int i = 0; i < pixels; i++) {
					dbg_img[3*fnslice+2][i] = half_strong_pix[fnslice][i]? 1.0 : 0.0;
				}
			}
		}
		if (dbg_img != null) {
			String [] dbg_titles = new String [3*num_slices];
			for (int nslice = 0; nslice < num_slices; nslice++) {
				dbg_titles[3*nslice+0] = "weak-"+nslice;
				dbg_titles[3*nslice+1] = "strong-"+nslice;
				dbg_titles[3*nslice+2] = "result-"+nslice;
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					width,
					height,
					true,
					"halfStrong",
					dbg_titles);
		}
		// 				String [] dbg_titles

		return half_strong_pix;
	}
	
	
	public static boolean [][] tileToPix( // expand tile selection to pixel selection 
			final boolean [][]  sel_tiles,
			final int           transform_size,
			final int           tilesX) {
		final int num_slices = sel_tiles.length;
		final int tiles = sel_tiles[0].length;
		final int tilesY = tiles / tilesX;
		final int width =  tilesX * transform_size;
		final int height = tilesY * transform_size;
		final int pixels = width * height;
		final boolean [][] sel_pixels = new boolean [num_slices][pixels];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs pn =     new TileNeibs(width, height);
		final boolean dbg = width < 0; // never
		final double [][] dbg_img = dbg? new double [3*num_slices] [pixels]:null;
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (sel_tiles[fnslice][tile]) {
							int tileX = tile % tilesX;
							int tileY = tile / tilesX;
							int indx0 = (tileY * width + tileX) * transform_size;
							for (int y = 0; y < transform_size; y++) {
								int indx1 = indx0 + y * width;  
								for (int x = 0; x < transform_size; x++) {
									int indx2 = indx1 + x;
									sel_pixels[fnslice][indx2] = true;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return sel_pixels;
	}
	
	public static boolean [][] getTrimPixels(
			final boolean [][] is_fg_pix,
			final boolean [][] has_bg_pix,
			final boolean [][] is_stitch_tile,
			final int          transform_size,
			final int          tilesX) {
		final boolean [][] trim_pixels = new boolean [is_fg_pix.length][];
		boolean [][] inv_stitch_pixels = null;
		if (is_stitch_tile != null) {
			inv_stitch_pixels = tileToPix(    // expand tile selection to pixel selection 
					is_stitch_tile,        // final boolean [][]  sel_tiles,
					transform_size,                    // final int           transform_size,
					tilesX);                           // final int           tilesX)
		}
		for (int nslice = 0; nslice < trim_pixels.length; nslice++) {
			trim_pixels[nslice] = has_bg_pix[nslice].clone();
			TileNeibs.andSelection(
					is_fg_pix[nslice], // final boolean [] src_tiles,
					trim_pixels[nslice] // final boolean [] dst_tiles
					);//
			if (is_stitch_tile != null) {
				TileNeibs.invertSelection(
						inv_stitch_pixels[nslice],
						inv_stitch_pixels[nslice]);
				TileNeibs.andSelection(
						inv_stitch_pixels[nslice], // final boolean [] src_tiles,
						trim_pixels[nslice]        // final boolean [] dst_tiles
						);
			}
		}
		return trim_pixels;
	}
	
	public static boolean [][]  getFgEdge( 
			final boolean [][]  fg_weak_tiles,
			final boolean [][]  fg_strong_tiles,
			final boolean [][]  stitch_tiles,
			final boolean [][]  fg_pix,
			final int           transform_size,
			final int           tilesX) {
		final int num_slices = fg_weak_tiles.length;
		final int tiles = fg_weak_tiles[0].length;
		final int tilesY = tiles / tilesX;
		final int width =  tilesX * transform_size;
		final int height = tilesY * transform_size;
		final int pixels = width * height;
		final boolean [][] fg_edge_pix = new boolean [num_slices][pixels];
		final boolean [] prohibit = new boolean [pixels];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs pn =     new TileNeibs(width, height);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			Arrays.fill(prohibit, true);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (fg_weak_tiles[fnslice][tile]) {
							boolean enable = 
									fg_weak_tiles[fnslice][tile] &&
									!fg_strong_tiles[fnslice][tile] &&
									!stitch_tiles[fnslice][tile] // no edge on the stitch
									;
							if (enable) {
								int tileX = tile % tilesX;
								int tileY = tile / tilesX;
								int indx0 = (tileY * width + tileX) * transform_size;
								for (int y = 0; y < transform_size; y++) {
									int indx1 = indx0 + y * width;  
									for (int x = 0; x < transform_size; x++) {
										int indx2 = indx1 + x;
										prohibit[indx2] = false;
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			System.arraycopy(fg_pix[fnslice], 0, fg_edge_pix[fnslice], 0, pixels);
			pn.growSelection( // will fill half of "weak" has_bg_tiles - it is where triangular mesh will reach
					2,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
					fg_edge_pix[fnslice],
					prohibit);
			boolean [] inv = TileNeibs.invertSelection(fg_pix[fnslice]);
			TileNeibs.andSelection(
					inv,                   // final boolean [] src_tiles,
					fg_edge_pix[fnslice]); // final boolean [] dst_tiles,
		}
		return fg_edge_pix;
	}
	
	public static boolean [][]  getTrimSeeds(
			final boolean [][]  trim_pix,  // pixels that may be trimmed
			final boolean [][]  seed_pix_in,  // FG edge, just outside of trim_pix. Will be modified. Or null
			final double  [][]  vars_same,
			final double  [][]  vars_inter,
			final double        seed_same_fz, // add to var_same in denominator
			final double        seed_fom,   // minimal value of vars_inter/sqrt(vars_same)
			final double        seed_inter, //  =   150; works for high-contrast over sky
			final int           width) {
		final int num_slices = trim_pix.length;
		final int num_pix =    trim_pix[0].length;
		final boolean [][]  seed_pix = (seed_pix_in != null) ? seed_pix_in : new boolean [num_slices][num_pix];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) if (trim_pix[fnslice][pix]) {
							if (vars_inter[fnslice][pix] > seed_inter) {  // works for high-contrast over sky
								seed_pix[fnslice][pix] = true;
							} else { // for foreground with lower contrast
								double fom = vars_inter[fnslice][pix]/(vars_same[fnslice][pix] + seed_same_fz);
								if (fom >= seed_fom) {
									seed_pix[fnslice][pix] = true;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return seed_pix;
	}

	
	
	
	public static double [][] getTrimFom(
			final boolean [][]  trim_pix,  // pixels that may be trimmed
			final double  [][]  vars_same,
			final double  [][]  vars_inter,
			final double        trim_inter_fz,  // minimal value of vars_same to block propagation
			final double        trim_fom_threshold, // = 120.0; // count only pixels with VAR_SAME > this value 
			final double        trim_fom_boost, // = 5;   // boost high-varinace values that exceed threshold  
			final double        trim_fom_blur,
			final int           width,
			final double [][][] fom_dbg) {
		final int num_slices = trim_pix.length;
		final int num_pix =    trim_pix[0].length;
		final int height = num_pix/width;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] trim_foms = new double [num_slices][num_pix];
		final double []   fom_blurred = (trim_fom_blur > 0) ? new double [num_pix] : null;
		final double boost_scale = trim_fom_boost / trim_fom_threshold;
		final double boost_subtract =  trim_fom_boost - 1.0;
		if (fom_dbg != null) {
			for (int i = 0; i <fom_dbg.length; i++) {
				fom_dbg[i] = new double[num_slices][];
			}
		}
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			Arrays.fill(trim_foms[fnslice], Double.NaN);
			if (trim_fom_blur > 0) {
				Arrays.fill(fom_blurred, 1.0); // trim_fom_threshold);
			}
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) if (!Double.isNaN(vars_same[fnslice][pix])){
							if (trim_pix[fnslice][pix]) { // should be NaN outside of trim_pix
								trim_foms[fnslice][pix] = vars_same[fnslice][pix]/(vars_inter[fnslice][pix] + trim_inter_fz);
							}
							if ((trim_fom_blur > 0) && !Double.isNaN(trim_foms[fnslice][pix])) {
								fom_blurred[pix] = Math.max(vars_same[fnslice][pix],trim_fom_threshold) * boost_scale - boost_subtract;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (fom_dbg != null) {
				fom_dbg[0][nslice] = trim_foms[fnslice].clone();
				fom_dbg[2][nslice] = fom_blurred.clone();
			}
			if (trim_fom_blur > 0) {
				(new DoubleGaussianBlur()).blurDouble(
						fom_blurred,
						width,
						height,
						trim_fom_blur,
						trim_fom_blur,
						0.01);
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) if (!Double.isNaN(trim_foms[fnslice][pix])){
								trim_foms[fnslice][pix] /= fom_blurred[pix];
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
				if (fom_dbg != null) {
					fom_dbg[3][nslice] = fom_blurred.clone();
					fom_dbg[1][nslice] = trim_foms[fnslice].clone();
				}
			}
		}
		return trim_foms;
	}
	
	public static boolean [][]  thresholdAnalog(
			final double  [][]  data,
			final double        threshold,
			final boolean       greater) {
		final int num_slices = data.length;
		final int num_pix =    data[0].length;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final boolean [][] binaries = new boolean [num_slices][num_pix];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			if (greater) {
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) {
								binaries[fnslice][pix] = data[fnslice][pix] >= threshold; 
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			} else {
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) {
								binaries[fnslice][pix] = data[fnslice][pix] <= threshold; 
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}
		}
		return binaries;
	}
		
	
	
	public static void getTrimAlpha(
			final double  [][]  fom_pix,   // should be NaN outside of trim_pix
			final boolean [][]  trim_pix,  // pixels that may be trimmed
			final boolean [][]  seed_pix,  // FG edge (just outside of trim_pix) and seeds from vars_inter mismatch
			final double        trim_fom, // minimal value of vars_same/vars_inter  to block propagation
			final int           trim_grow,      // 3*transform_size?
			final int           width) {
		final int num_slices = trim_pix.length;
		final int num_pix =    trim_pix[0].length;
		final int height = num_pix/width;
		final boolean [] prohibit = new boolean[num_pix];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs pn =     new TileNeibs(width, height);
		final double [][] foms = new double [2][num_pix];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			Arrays.fill(prohibit, true);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) {
							 // NaN should keep prohibit (fom_pix should be NaN outside of trim_pix
							if (((fom_pix[fnslice][pix] <= trim_fom) && trim_pix[fnslice][pix])
									|| seed_pix[fnslice][pix]) {
								prohibit[pix] = false;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			// grow seed
		    pn.growSelection(
		    		trim_grow,
		    		seed_pix[fnslice],
		            prohibit);
		    TileNeibs.andSelection(
		    		trim_pix[fnslice],  // final boolean [] src_tiles,
		    		seed_pix[fnslice]); // final boolean [] dst_tiles,
		    TileNeibs.invertSelection( // invert in place
		    		seed_pix[fnslice],
		    		seed_pix[fnslice]);
		}
		return; //  seed_pix; // extends outside selected tiles, but that's OK
	}
	
	public static void expandTrimAlpha(
			final boolean [][]  trim_pix,   // pixels that may be trimmed
			final boolean [][]  alpha_pix,  //
			final double  [][]  value,      // will grow only in increasing
			final double        min_incr,
			final boolean       dual_pass,
			final int           trim_grow,      // 3*transform_size?
			final int           width) {
		final int num_slices = trim_pix.length;
		final int num_pix =    trim_pix[0].length;
		final int height = num_pix/width;
		final boolean [] removed =  new boolean[num_pix];
		final boolean [] prohibit = new boolean[num_pix];
		final TileNeibs pn =     new TileNeibs(width, height);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
		    TileNeibs.invertSelection(
		    		alpha_pix[fnslice],
		    		removed); // 
		    TileNeibs.invertSelection(
		    		trim_pix[fnslice],
		    		prohibit); // prohibit all that is not trim
		    pn.growSelectionGradient( //
		    		trim_grow,          // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
		    		removed,            // final boolean [] tiles,
		    		prohibit,           // final boolean [] prohibit,
		    		value[fnslice],     // final double []  value,
		    		min_incr,           // final double     min_incr,
					true,               // final boolean    keep_top,     // do not grow over local max
					false);             // final boolean    keep_thin_top)// do not grow over local max only if next after that is already selected
			if (dual_pass) {
			    pn.growSelectionGradient( //
			    		2, // trim_grow,    // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			    		removed,            // final boolean [] tiles,
			    		prohibit,           // final boolean [] prohibit,
			    		value[fnslice],     // final double []  value,
			    		min_incr,           // final double     min_incr,
						false,              // final boolean    keep_top,     // do not grow over local max
						true);              // final boolean    keep_thin_top)// do not grow over local max only if next after that is already selected
			}
		    TileNeibs.invertSelection(
		    		removed,
		    		alpha_pix[fnslice]); // 
		}
		return;
	}

	
	
	
	public static void filterAlpha(
			final boolean [][]  alpha_pix,  // per-pixel alpha
			final boolean [][]  trim_pix,   // pixels that may be trimmed
			final int           min_neibs,  // minimal neighbors to keep alpha
			final int           grow_alpha, // grow alpha selection
			final int           width) {
		final int num_slices = alpha_pix.length;
		final int num_pix =    alpha_pix[0].length;
		final int height = num_pix/width;
		final boolean [] btmp = (min_neibs > 0) ? new boolean[num_pix] : null;
		final boolean [] prohibit = new boolean[num_pix];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs pn =     new TileNeibs(width, height);
		final int dbg_tile = -82228; // 71992; //312/112 or 61800 for 360/96
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			if (btmp != null) {
				Arrays.fill(btmp, false);
				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int pix = ai.getAndIncrement(); pix < num_pix; pix = ai.getAndIncrement()) if (alpha_pix[fnslice][pix]){
								if ((pix == dbg_tile) || (pix == (dbg_tile-1))){
								    System.out.println("filterAlpha(): pix="+pix);
								}
								int nneibs = 0;
								for (int dir = 0; dir < TileNeibs.DIR_S; dir++) {
									int pix1 = pn.getNeibIndex(pix, dir);
									if ((pix1 >= 0) && alpha_pix[fnslice][pix1]) {
										if (++nneibs >= min_neibs) {
											btmp[pix] = true;
											break;
										}
									}
								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
				System.arraycopy(btmp, 0, alpha_pix[fnslice], 0, num_pix);
			} // if (btmp != null) {
			
		    TileNeibs.invertSelection( // invert in place
		    		trim_pix[fnslice],
		    		prohibit);
			if (grow_alpha > 0) {
				pn.growSelection(
						grow_alpha,
						alpha_pix[fnslice],
						prohibit);
			}
		}
	}

		
	
	
	public static void  filterWeakFG(
			final boolean [][]  alpha_pix,
			final boolean [][]  trim_pix,			
			final boolean [][]  has_bg_pix,			
			final boolean [][]  fg_weak_tiles,
			final boolean [][]  fg_strong_tiles,
			final int           transform_size,
			final int           tilesX) {
		final int num_slices = fg_weak_tiles.length;
		final int tiles = fg_weak_tiles[0].length;
		final int tilesY = tiles / tilesX;
		final int width =  tilesX * transform_size;
		final int height = tilesY * transform_size;
		final int pixels = width * height;
		final boolean [][] alpha = new boolean [num_slices][pixels];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice =0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement())
							 // tiles that may be half trim_pix, other half should be removed
							if (fg_weak_tiles[fnslice][tile] && !fg_strong_tiles[fnslice][tile]) {
							int tileX = tile % tilesX;
							int tileY = tile / tilesX;
							int indx0 = (tileY * width + tileX) * transform_size;
							for (int y = 0; y < transform_size; y++) {
								int indx1 = indx0 + y * width;
								System.arraycopy(alpha_pix[fnslice], indx1, alpha[fnslice], indx1, transform_size);
								for (int x = 0; x < transform_size; x++) {
									int indx2 = indx1 + x;
									if (has_bg_pix[fnslice][indx2] && !trim_pix[fnslice][indx2]) {
										alpha_pix[fnslice][indx2] = false;
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

	public static boolean [][]  trimAlphaToTiles(
			final boolean [][]  alpha_pix,
			final boolean [][]  selected_tiles,
			final int           transform_size,
			final int           tilesX) {
		final int num_slices = selected_tiles.length;
		final int tiles = selected_tiles[0].length;
		final int tilesY = tiles / tilesX;
		final int width =  tilesX * transform_size;
		final int height = tilesY * transform_size;
		final int pixels = width * height;
		final boolean [][] alpha = new boolean [num_slices][pixels];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice =0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (selected_tiles[fnslice][tile]) {
							int tileX = tile % tilesX;
							int tileY = tile / tilesX;
							int indx0 = (tileY * width + tileX) * transform_size;
							for (int y = 0; y < transform_size; y++) {
								int indx1 = indx0 + y * width;
								System.arraycopy(alpha_pix[fnslice], indx1, alpha[fnslice], indx1, transform_size);
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return alpha;
	}

	/**
	 * Fix gap (copy alpha row/col) between meshes were strong has_bg is on the edge.
	 * Alpha should already be set everywhere needed.
	 *   
	 * @param alpha_pix
	 * @param selected_tiles
	 * @param has_bg_strong_tiles 
	 * @param slice_border_int    to prevent bridging disparity gap 
	 * @param neib_max maximal neib_lev, now 2 (there can be neib_max+1 - outmost BG) 
	 * @param transform_size
	 * @param tilesX
	 * @return
	 */
	public static void  fix_bg_overlap(
			final boolean [][]  alpha_pix,
			final boolean [][]  selected_tiles,
			final boolean [][]  has_bg_strong_tiles,
			final int      [][] slice_border_int,
			final int           neib_max, // now 2
			final int           transform_size,
			final int           tilesX) {
		final int num_slices = selected_tiles.length;
		final int tiles = selected_tiles[0].length;
		final int tilesY = tiles / tilesX;
		final int width =  tilesX * transform_size;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		final boolean [] new_sel = new boolean [tiles];
		final int [] corner =   {0, transform_size - 1, (transform_size - 1) * (width + 1), (transform_size - 1) * width };
		final int [] src_offs =  {-width, 1,      width, -1};
		final int [] step_offs = {      1, width, -1,     -width};
		final int outmost_fg = neib_max;
		final int outmost_bg = neib_max+1;
		final boolean [][] compat_outmost = new boolean[neib_max+2][neib_max+2];
		for (int i = 0; i < compat_outmost.length; i++) {
			for (int j = i; j < compat_outmost.length; j++) {
				compat_outmost[i][j] = ! (((i == outmost_fg) && (j == outmost_bg)) || ((j == outmost_fg) && (i == outmost_bg)));
				compat_outmost[j][i] = compat_outmost[i][j]; 
			}
		}
		for (int nslice =0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			Arrays.fill(new_sel,false);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement())
							if (has_bg_strong_tiles[fnslice][tile] && !selected_tiles[fnslice][tile]) {
								int px = (tile % tilesX) * transform_size;
								int py = (tile / tilesX) * transform_size;
								int indx0 = px + py * width; 
								for (int dir2 = 0; dir2 < TileNeibs.DIRS; dir2 += 2) {
									int tile1 = tn.getNeibIndex(tile, dir2);
									if ((tile1 >= 0) &&
											selected_tiles[fnslice][tile1] &&
											compat_outmost[slice_border_int[fnslice][tile]][slice_border_int[fnslice][tile1]])
									{
										new_sel[tile] = true;
										int dir = dir2 / 2;
										// copy existing row from that direction
										int dindx = indx0 + corner[dir];
										int sindx = dindx + src_offs[dir];
										for (int i = 0; i < transform_size; i++) {
											alpha_pix[fnslice][dindx] |= alpha_pix[fnslice][sindx];
											sindx += step_offs[dir];
											dindx += step_offs[dir];
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
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement())
							if (new_sel[tile]) {
								selected_tiles[fnslice][tile] = true;
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
	}
	
	
	

	
	//final double [][] slice_disparities,	
	public static void showDebugDisparities(
			final double [][] slice_disparities,
			final int   tilesX,
			String      prefix) {

		if (prefix != null) {
			final int num_slices = slice_disparities.length;
			final int tiles = slice_disparities[0].length;
			final int tilesY = tiles / tilesX;
			String [] dbg_titles = new String[num_slices];
			for (int i = 0; i < dbg_titles.length; i++) {
				dbg_titles[i] = String.format("slice-%02d",i);
			}
			ShowDoubleFloatArrays.showArrays(
					slice_disparities,
					tilesX,
					tilesY,
					true,
					prefix+"-slice_disparities",
					dbg_titles);

		}	
	}
		
	public static void showDebugFgBg(
			boolean [][][] fg_has_bg,
			final int   tilesX,
			String      prefix) {
		if (prefix != null) {
			final int num_slices = fg_has_bg[0].length;
			final int tiles = fg_has_bg[0][0].length;
			final int tilesY = tiles / tilesX;

			String [] dbg_titles = new String[num_slices];
			for (int i = 0; i < dbg_titles.length; i++) {
				dbg_titles[i] = String.format("slice-%02d",i);
			}
			double [][] dbg_fgbg= new double [num_slices][tiles];
			for (int ns = 0; ns < num_slices; ns++) {
				Arrays.fill(dbg_fgbg[ns],Double.NaN);
				for (int tile = 0; tile< tiles; tile++) {
					dbg_fgbg[ns][tile] = (fg_has_bg[2][ns][tile]?0.2:0) + (fg_has_bg[0][ns][tile]?2:0) + (fg_has_bg[1][ns][tile]?1:0);
				}
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_fgbg,
					tilesX,
					tilesY,
					true,
					prefix+"-fgbg",
					dbg_titles);

		}		
	}		

	public static void showDebugBordersInt(
			final int [][] border_int,
			final int      tilesX,
			String         prefix) {
		if (prefix != null) {
			final int num_slices = border_int.length;
			final int tiles = border_int[0].length;
			final int tilesY = tiles / tilesX;

			String [] dbg_titles = new String[num_slices];
			for (int i = 0; i < dbg_titles.length; i++) {
				dbg_titles[i] = String.format("slice-%02d",i);
			}
			double [][] dbg_img= new double [num_slices][tiles];
			for (int ns = 0; ns < num_slices; ns++) {
				Arrays.fill(dbg_img[ns],Double.NaN);
				for (int tile = 0; tile< tiles; tile++)  if (border_int[ns][tile] >= 0){
					dbg_img[ns][tile] = border_int[ns][tile];
				}
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					prefix+"-border_int",
					dbg_titles);

		}		
	}		
	
	
	public static void showDebugBoolean(
			boolean [][] selected,
			final int   tilesX,
			String      prefix) {
		if (prefix != null) {
			final int num_slices = selected.length;
			final int tiles = selected[0].length;
			final int tilesY = tiles / tilesX;

			String [] dbg_titles = new String[num_slices];
			for (int i = 0; i < dbg_titles.length; i++) {
				dbg_titles[i] = String.format("slice-%02d",i);
			}
			double [][] dbg_selected= new double [num_slices][tiles];
			for (int ns = 0; ns < num_slices; ns++) {
				Arrays.fill(dbg_selected[ns],Double.NaN);
				for (int tile = 0; tile< tiles; tile++) {
					dbg_selected[ns][tile] = (selected[ns][tile]?1:0);
				}

			}
			ShowDoubleFloatArrays.showArrays(
					dbg_selected,
					tilesX,
					tilesY,
					true,
					prefix+"-boolean",
					dbg_titles);
		}		
	}		
	

	public static double[][][] getVariances (
			final double [][][] sensor_texture,
			final double [][]   combo_texture,
			final double        var_radius,
			final int           width) { // last slice - ratio
		final int num_slices =  sensor_texture.length;
		final int num_sensors = sensor_texture[0].length;
		final int img_size =    sensor_texture[0][0].length;
		final int height =      img_size/width;
		final int ivar_radius = (int) Math.floor(var_radius);
		final double [][] var_weights = new double [ivar_radius+1][ivar_radius+1];
		for (int i = 0; i < var_weights.length; i++) {
			for (int j = 0; j < var_weights[i].length; j++) {
				var_weights[i][j] = Math.cos(0.5*Math.PI*i/var_radius) * Math.cos(0.5*Math.PI*j/var_radius);
			}
		}
		final double [][][] vars = new double [5][num_slices][img_size]; // same, inter, d/dx, d/dy,sqrt ((d/dx)^2+(d/dy)^2)
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
 		final TileNeibs pn =     new TileNeibs(width, height);

		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			for (int i = 0; i < vars.length; i++) {
				Arrays.fill(vars[i][nslice], Double.NaN);
			}
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int cindx = ai.getAndIncrement(); cindx < img_size; cindx = ai.getAndIncrement())
							if (!Double.isNaN(combo_texture[fnslice][cindx])) {
								int y0 = cindx / width;
								int x0 = cindx % width;
								double var_same = Double.NaN;
								double var_inter = Double.NaN;
								// calculate weighted variance
								double sw = 0.0, swd=0.0, swd2 = 0.0, swxd = 0.0, swyd = 0.0;
								for (int dvy = -ivar_radius; dvy <= ivar_radius; dvy++) {
									for (int dvx = -ivar_radius; dvx <= ivar_radius; dvx++) {
										int indx = pn.getIndex(x0+dvx, y0+dvy);
										if ((indx >= 0) && !Double.isNaN(combo_texture[fnslice][indx])) {
											double w = var_weights[Math.abs(dvy)][Math.abs(dvx)]; // 1.0;
											double d = combo_texture[fnslice][indx];
											double wd = w * d;
											sw += w;
											swd += wd;
											swd2 += wd*d;
											swxd += dvx * wd;
											swyd += dvy * wd;
										}
									}
								}
								if (sw > 0.0) { // always
									double avg =  swd/sw;
									double avg2 = swd2/sw;
									swxd /= sw;
									swyd /= sw;
									var_same = Math.sqrt(avg2-avg*avg);
									vars[0][fnslice][cindx] = var_same;
									vars[2][fnslice][cindx] = swxd;
									vars[3][fnslice][cindx] = swyd;
									vars[4][fnslice][cindx] = Math.sqrt(swxd*swxd + swyd*swyd);
								}
								// calculate inter-sensor variance (add local normalization?)
								sw = 0.0; swd=0.0; swd2 = 0.0;
								for (int nsens = 0; nsens< num_sensors; nsens++) {
									double w = 1.0;
									double d = sensor_texture[fnslice][nsens][cindx];
									sw += w;
									swd += w * d;
									swd2 += w * d*d;
								}
								if (sw > 0.0) { // always
									double avg =  swd/sw;
									double avg2 = swd2/sw;
									var_inter = Math.sqrt(avg2-avg*avg);
									vars[1][fnslice][cindx] = var_inter;
								}
							}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		
		return vars;	
	}
	
	public static boolean [][] getAllTexturePixels (
			final double [][]   combo_texture,
			final boolean [][]  tile_on,
			final int           width,
			final int           transform_size) {
		final int num_slices =  combo_texture.length;
		final int img_size =    combo_texture[0].length;
		final int tilesX = width / transform_size;
		final boolean [][] texture_on = new boolean[num_slices][img_size];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int cindx = ai.getAndIncrement(); cindx < img_size; cindx = ai.getAndIncrement())
							if (!Double.isNaN(combo_texture[fnslice][cindx])) {
								texture_on[fnslice][cindx] = true;
								if (tile_on != null) {
									int tileY = (cindx / width) / transform_size;
									int tileX = (cindx % width) / transform_size;
									if (!tile_on[fnslice][tileX + tilesX * tileY]) {
										texture_on[fnslice][cindx] = false;
									}
								}
							}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}		
		return texture_on;	
	}	
	
	public static boolean [][] getEdgeTexturePixels(
			final boolean [][] texture_on,
			final int          width,
			final int          trim_edge) {
		final int num_slices =  texture_on.length;
		final int img_size =    texture_on[0].length;
		final boolean [][] texture_edge = new boolean[num_slices][];
		final TileNeibs pn =     new TileNeibs(width, img_size/width);
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0); // may remove multithreaded
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < num_slices; nslice = ai.getAndIncrement()) {
						texture_edge[nslice]=pn.getEdgeSelection(
								trim_edge,          // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
								texture_on[nslice], // boolean [] tiles,
								null);     // boolean [] prohibit);
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return texture_edge;
	}
	
	public static boolean [][] filterFgByVariances(
			final double  [][] combo_texture,
			final boolean [][] texture_on, // if null will rely on NaN in combo_texture
			final double  [][] vars_same,
			final double  [][] vars_inter,
			final boolean [][] is_fg_tile,
			final boolean [][] has_bg_tile,
			final double       fg_max_inter,
			final double       fg_max_rel,
			final int          width,
			final int          transform_size) {
		final int num_slices =  texture_on.length;
		final int img_size =    texture_on[0].length;
		final int tilesX = width / transform_size;
		final boolean [][] texture_filt = new boolean[num_slices][img_size];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int cindx = ai.getAndIncrement(); cindx < img_size; cindx = ai.getAndIncrement())
							// is_fg_tile[fnslice][tile]
							if (!Double.isNaN(combo_texture[fnslice][cindx]) && ((texture_on==null)|| texture_on[fnslice][cindx]))  {
								texture_filt[fnslice][cindx] = true;
								int y0 = cindx / width;
								int x0 = cindx % width;
								int tileX = x0 / transform_size;
								int tileY = y0 / transform_size;
								int tile = tileX + tilesX * tileY;
								// only trim if nothing obscures this and has some BG 
								if (is_fg_tile[fnslice][tile] && has_bg_tile[fnslice][tile]) {
									if (vars_inter[fnslice][cindx] > fg_max_inter) {
										texture_filt[fnslice][cindx] = false;
									}
									if (vars_inter[fnslice][cindx]/vars_same[fnslice][cindx] > fg_max_rel) {
										texture_filt[fnslice][cindx] = false;
									}
								}
							}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
		}		
		return texture_filt;
	}

	public static void trimFgEdgePixels(
			final boolean [][] texture_filt,
			final boolean [][] texture_edge,
			final boolean [][] is_fg_tile,
			final boolean [][] has_bg_tile,
			final int          trim_edge_center,			
			final int          width,
			final int          transform_size) {
		final int num_slices =  texture_filt.length;
		final int img_size =    texture_filt[0].length;
		final int tilesX = width / transform_size;
		final int tiles = img_size / transform_size / transform_size;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						int tile_center_offs = ((transform_size / 2) - 1)* (width + 1); // 1923 
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							if (is_fg_tile[fnslice][tile] && has_bg_tile[fnslice][tile]) {
								int tileX = tile % tilesX;
								int tileY = tile / tilesX;
								int indx0 = (tileY * transform_size) * width + (tileX * transform_size);
								int indx1 = indx0  + tile_center_offs;
								int num_cent = 0;
								for (int dy = 0; dy < 2; dy++) {
									for (int dx = 0; dx < 2; dx++) {
										if (texture_filt[fnslice][indx1 + width * dy + dx]) {
											num_cent++;
										}
									}
								}
								if (num_cent < trim_edge_center) {
									for (int dy = 0; dy < transform_size; dy++) {
										for (int dx = 0; dx < transform_size; dx++) {
											int indx = indx0+ width*dy + dx;
											if (texture_edge[fnslice][indx]) {
												texture_filt[fnslice][indx] = false;
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
	}
	
	/**
	 * Trim low-variance FG tiles near the edge
	 * @param texture_en
	 * @param texture_edge
	 * @param is_fg_tile
	 * @param has_bg_tile
	 * @param variance
	 * @param min_variance
	 * @param width
	 * @param transform_size
	 */
	
	public static void trimFgEdgeVariancePixels(
			final boolean [][] texture_en,
			final boolean [][] texture_edge,
			final boolean [][] is_fg_tile,
			final boolean [][] has_bg_tile,
			final double  [][] variance,
			final double       min_edge_variance,			
			final int          width,
			final int          transform_size) {
		boolean dbg = false; // true;
		final int num_slices =      texture_en.length;
		final int img_size =        texture_en[0].length;
		final int tilesX = width/transform_size;
		final int tilesY = img_size/width/transform_size;
		final Thread[] threads =    ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai =    new AtomicInteger(0);
		final TileNeibs pn =        new TileNeibs(width, img_size/width);
		final int grow = 2 * width;
		final String [] dbg_titles = {"sel_in", "edge", "seed", "prohibit", "erase", "sel_out"}; 
		for (int ithread = 0; ithread < threads.length; ithread++) {
			final boolean [] erase =     new boolean [img_size];
			final boolean [] prohibit = new boolean [img_size];
			final double [][] dbg_img = dbg? (new double [dbg_titles.length][img_size]): null;
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < num_slices; nslice = ai.getAndIncrement()) {
						Arrays.fill(prohibit, true);
						Arrays.fill(erase,false);
						if (dbg_img != null) {
							for (int i = 0; i < dbg_img.length; i++) {
								Arrays.fill(dbg_img[i], Double.NaN);
							}
							for (int i = 0; i < texture_en[nslice].length; i++) {
								dbg_img[0][i] = texture_en[nslice][i]? 1.0 : 0.0;
								dbg_img[1][i] = texture_edge[nslice][i]? 1.0 : 0.0;
							}
						}
						for (int tileY=0; tileY <tilesY; tileY++) {
							for (int tileX=0; tileX <tilesX; tileX++) {
								int tile = tileY*tilesX+tileX;
								if (is_fg_tile[nslice][tile] && has_bg_tile[nslice][tile]) {
									int indx0 = (tileY * transform_size) * width + (tileX * transform_size);
									for (int dy = 0; dy < transform_size; dy++) {
										for (int dx = 0; dx < transform_size; dx++) {
											int indx = indx0+ width*dy + dx;
											prohibit[indx] = (variance[nslice][indx] >= min_edge_variance) ||
													!texture_en[nslice][indx];
											erase[indx] = !prohibit[indx] && texture_edge[nslice][indx];
										}								
									}
								}
							}
						}
						if (dbg_img != null) {
							for (int i = 0; i < texture_en[nslice].length; i++) {
								dbg_img[2][i] = erase[i]? 1.0 : 0.0;
								dbg_img[3][i] = prohibit[i]? 1.0 : 0.0;
							}
						}
						
						pn.growSelection(
								grow,       // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
								erase,       // boolean [] tiles,
								prohibit);  // boolean [] prohibit);
						for (int i = 0; i < erase.length; i++) {
							texture_en[nslice][i] &= !erase[i];
						}
						if (dbg_img != null) {
							for (int i = 0; i < texture_en[nslice].length; i++) {
								dbg_img[4][i] = erase[i]? 1.0 : 0.0;
								dbg_img[5][i] = texture_en[nslice][i]? 1.0 : 0.0;
							}

							ShowDoubleFloatArrays.showArrays(
									dbg_img,
									width,
									img_size/width,
									true,
									"VAR_EDGE_FILTER-"+String.format("%02d",nslice),
									dbg_titles);
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	public static void filterSelections(
			final boolean [][] selections, // ** will be modified
			final int          min_neibs,
			final int          grow,
			final int          shrink,
			final int          width) {
		final int num_slices =  selections.length;
		final TileNeibs pn =     new TileNeibs(width, selections[0].length/width);
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < num_slices; nslice = ai.getAndIncrement()) {
						pn.removeFewNeibs(
								selections[nslice], // boolean [] selection, // should be the same size
								min_neibs); // int min_neibs)
						if (grow > 0) {
							pn.growSelection(
									grow,          // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
									selections[nslice], // boolean [] tiles,
									null);     // boolean [] prohibit);
						}
						if (shrink > 0) {
							pn.shrinkSelection(
									shrink,          // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
									selections[nslice], // boolean [] tiles,
									null);     // boolean [] prohibit);
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	public static void applyTextureSelection(
			final boolean [][]   selections,
			final double  [][]   combo_texture // will be modified - NaN where not selected
			) {
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < selections.length; nslice = ai.getAndIncrement()) {
						for (int i = 0; i < selections[nslice].length; i++) {
							if (!selections[nslice][i]) {
								combo_texture[nslice][i] = Double.NaN;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public static double [][] combineTextureAlpha(
			final double         alpha_threshold,
			final double  [][]   textures,
			final double  [][]   alphas
			) {
		final double [][] masked_textures = new double [textures.length][];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < textures.length; nslice = ai.getAndIncrement()) {
						masked_textures[nslice] = textures[nslice].clone();						
						for (int i = 0; i < textures[nslice].length; i++) {
							if (alphas[nslice][i] < alpha_threshold) {
								masked_textures[nslice][i] = Double.NaN;
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return masked_textures;
	}
	
	
	
	public static void extendBlueSKy(
			final TileCluster [] tileClusters,
			final double [][][]  faded_textures,
			final int            shrink_sky_tiles,
			final int            width,
			final int            transform_size) {
//		final int img_size =        faded_textures[0].length;
//		final int tilesX = width/transform_size;
//		final int tilesY = img_size/width/transform_size;
//		final int tiles =  tilesX * tilesY;
		
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		int sky_slice = -1;
		int sky_subindex = -1;
		for (int i = 0; i < tileClusters.length; i++) {
			sky_subindex = tileClusters[i].getSkyClusterIndex();
			if (sky_subindex >=0) {
				sky_slice = i;
				break;
			}
		}
		if (sky_slice >= 0) {
			Rectangle sky_tiles_bounds = tileClusters[sky_slice].getSubBounds(sky_subindex);
			Rectangle sky_pixels_bounds = new Rectangle(
					sky_tiles_bounds.x *      transform_size,
					sky_tiles_bounds.y *      transform_size,
					sky_tiles_bounds.width *  transform_size,
					sky_tiles_bounds.height * transform_size
					);
			final double [] sky_disparity = tileClusters[sky_slice].getSubDisparity(sky_subindex);
			final int num_sky_tiles =    sky_disparity.length;
			final TileNeibs tn_sky_tiles =      new TileNeibs(sky_tiles_bounds.width, sky_tiles_bounds.height);
//			final TileNeibs tn_sky_pixels =     new TileNeibs(sky_pixels_bounds.width, sky_pixels_bounds.height);
			final boolean [] sky_sel = new boolean [sky_disparity.length];
			for (int i = 0; i < sky_sel.length; i++) {
				sky_sel[i] = !Double.isNaN(sky_disparity[i]);
			}
			final boolean [] reliable_sky = sky_sel.clone();
			final double [] sky_pixels = TileNeibs.getDoubleWindow(
					sky_pixels_bounds,            // Rectangle window,
					faded_textures[sky_slice][0], // double [] data,
					width);                       // int data_width)
			if (shrink_sky_tiles > 0) {
				tn_sky_tiles.shrinkSelection(
						shrink_sky_tiles,
						reliable_sky,
						null);
				for (int i = 0; i < reliable_sky.length; i++) if (!reliable_sky[i]){
					sky_disparity[i] = Double.NaN;
				}
				// erase corresponding texture pixels

				ai.set(0);
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int sky_tile = ai.getAndIncrement(); sky_tile < num_sky_tiles; sky_tile = ai.getAndIncrement()) {
								if (sky_sel[sky_tile] && !reliable_sky[sky_tile]) {
									int sky_tileX = sky_tile % sky_tiles_bounds.width; 
									int sky_tileY = sky_tile / sky_tiles_bounds.width;
									int pixX = sky_tileX * transform_size;
									int pixY = sky_tileY * transform_size;
									int pix_indx0 = pixY * sky_pixels_bounds.width + pixX;
									for (int row = 0; row < transform_size; row++) {
										int pix_indx = pix_indx0+ row * sky_pixels_bounds.width;
										Arrays.fill(sky_pixels,
												pix_indx, //  int fromIndex,
												pix_indx + transform_size, //  int toIndex,
												Double.NaN);
									}

								}
							}
						}
					};
				}		      
				ImageDtt.startAndJoin(threads);
			}
			// now fill gaps in disparity and pixels
			double []sky_disparity_filled = TileProcessor.fillNaNs( // multithreaded
					sky_disparity,            // final double [] data,
					null,                     // final boolean [] prohibit,
					sky_tiles_bounds.width,   // int       width,
					2 * Math.max(sky_tiles_bounds.width, sky_tiles_bounds.height), // 16,           // final int grow,
					0.7,                      // double    diagonal_weight, // relative to ortho
					100,                      // int       num_passes,
					0.01,                     // final double     max_rchange, //  = 0.01
					THREADS_MAX);             // final int threadsMax)      // maximal number of threads to launch 
			tileClusters[sky_slice].setSubDisparity(sky_subindex,sky_disparity_filled);

			double[] sky_pixels_filled = TileProcessor.fillNaNs( // multithreaded
					sky_pixels,               // final double [] data,
					null,                     // final boolean [] prohibit,
					sky_pixels_bounds.width,  // int       width,
					3 * Math.min(sky_pixels_bounds.width,sky_pixels_bounds.height) / 2, // 16,           // final int grow,
					0.7,                      // double    diagonal_weight, // relative to ortho
					100,                      // int       num_passes,
					0.01,                     // final double     max_rchange, //  = 0.01
					THREADS_MAX);             // final int threadsMax)      // maximal number of threads to launch 
			
			TileNeibs.setDoubleWindow(
					sky_pixels_bounds,            // Rectangle window,
					sky_pixels_filled,            // double [] window_data,
					faded_textures[sky_slice][0], // double [] data,
					width);                       // int data_width)
			double [] sky_alpha = new double [sky_pixels.length];
			Arrays.fill(sky_alpha, 1.0);
			TileNeibs.setDoubleWindow(
					sky_pixels_bounds,            // Rectangle window,
					sky_alpha,                   // double [] window_data,
					faded_textures[sky_slice][1], // double [] data,
					width);                       // int data_width)
		}
		return;
	}
	
	// Old version with analog alpha
	public static void fixAlphaSameDisparity(
			final TileCluster [] tileClusters,
			final double [][][]  faded_textures,
			final double         alphaOverlapTolerance, // 0 - require exact match
			final int            width,
			final int            transform_size) {
		final int num_slices = tileClusters.length;
		final int img_size =        faded_textures[0][0].length;
		final int tilesX = width/transform_size;
		final int tilesY = img_size/width/transform_size;
		final int tiles =  tilesX * tilesY;
		
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] tile_disparity = new double [num_slices][];
		for (int nslice = 0; nslice < tileClusters.length; nslice++) {
			tile_disparity[nslice] = tileClusters[nslice].getDisparity(); 
		}
		final double min_disparity = 2.0;
		
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [] disparities = new double [num_slices];
					int group [] = new int [num_slices];
					int [][] members = new int [num_slices][];
					int [] group_members = new int [num_slices];
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
						int num_tiles = 0;
						for (int nslice = 0; nslice < num_slices; nslice++) {
							disparities[nslice] = tile_disparity[nslice][tile];
							if (!Double.isNaN(disparities[nslice])) {
								num_tiles++;
							}
						}
						if (num_tiles > 1) {
							Arrays.fill(group, -1);
							int ngroup = 0;
							for (int i = 0; i < (num_slices-1); i++) if ((group[i] < 0) && !Double.isNaN(disparities[i])){
								int nsame = 0;
								group_members[nsame++] = i;
								double max_diff = Math.max(disparities[i], min_disparity) * alphaOverlapTolerance; 
								for (int j = i+1; j < num_slices; j++) if ((group[j] < 0) && !Double.isNaN(disparities[j])){
									boolean same = (alphaOverlapTolerance == 0) ?
											(disparities[j] == disparities[i]) :
												(Math.abs(disparities[j] - disparities[i]) < max_diff);
									if (same) {
										group[j] = ngroup;
										group_members[nsame++] = j;
									}
								}
								if (nsame > 1) {
									members[ngroup]=new int[nsame];
									group[i] = ngroup;
									for (int j = 0; j < nsame; j++) {
										members[ngroup][j] = group_members[j];
									}
									ngroup++;
								}
							}
							if (ngroup > 0) {
								int y0 = (tile / tilesX) * transform_size;
								int x0 = (tile % tilesX) * transform_size;
								int indx0 = y0 * width + x0; 
								for (int ng =0; ng < ngroup; ng++) {
									for (int dy = 0; dy < transform_size; dy++) {
										int indx1 = indx0 + dy * width;
										for (int dx = 0; dx < transform_size; dx++) {
											int indx = indx1 + dx;
											double a = faded_textures[members[ng][0]][1][indx];
											for (int j = 1; j < members[ng].length; j++) {
												a = Math.max(a, faded_textures[members[ng][j]][1][indx]);
											}
											for (int j = 0; j < members[ng].length; j++) {
												faded_textures[members[ng][j]][1][indx] = a;
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
		return;
	}	
	

	public static void setMeshTileSelection(
			final double  [][] slice_disparities,
			final boolean [][] keep_tiles) {
		final int num_slices = slice_disparities.length;
		final int num_tiles =  slice_disparities[0].length;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++){
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < num_tiles; tile = ai.getAndIncrement()) {
							if (!keep_tiles[fnslice][tile]) {
								slice_disparities[fnslice][tile] = Double.NaN;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
	}
	
	
	// new version with binary alpha. I
	public static void fixAlphaSameDisparity(
			final TileCluster [] tileClusters,
			final boolean [][]   keep_tiles, 
			final boolean [][]   alpha_pix,
			final boolean        use_or, // (maximal alpha), false - and (minimal alpha)
			final double         alphaOverlapTolerance, // 0 - require exact match
			final int            width,
			final int            transform_size) {
		final int num_slices = tileClusters.length;
		final int img_size =        alpha_pix[0].length;
		final int tilesX = width/transform_size;
		final int tilesY = img_size/width/transform_size;
		final int tiles =  tilesX * tilesY;
		
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] tile_disparity = new double [num_slices][];
		for (int nslice = 0; nslice < tileClusters.length; nslice++) {
			tile_disparity[nslice] = tileClusters[nslice].getDisparity(); 
		}
		final double min_disparity = 2.0;
		final int dbg_tile=-2122; // : 42/26 // 1872: 32/23
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [] disparities = new double [num_slices];
					int group [] = new int [num_slices];
					int [][] members = new int [num_slices][];
					int [] group_members = new int [num_slices];
					for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
						if (tile==dbg_tile) {
							System.out.println("fixAlphaSameDisparity(): tile="+tile);
						}
						int num_tiles = 0;
						for (int nslice = 0; nslice < num_slices; nslice++) {
							disparities[nslice] = tile_disparity[nslice][tile];
							if (!Double.isNaN(disparities[nslice]) && keep_tiles[nslice][tile]) {
								num_tiles++;
							}
						}
						if (num_tiles > 1) {
							Arrays.fill(group, -1);
							int ngroup = 0;
							for (int i = 0; i < (num_slices-1); i++) if ((group[i] < 0) && !Double.isNaN(disparities[i]) && keep_tiles[i][tile]){
								int nsame = 0;
								group_members[nsame++] = i;
								double max_diff = Math.max(disparities[i], min_disparity) * alphaOverlapTolerance; 
								for (int j = i+1; j < num_slices; j++)
									if ((group[j] < 0) && !Double.isNaN(disparities[j]) && keep_tiles[j][tile]){
										boolean same = (alphaOverlapTolerance == 0) ?
												(disparities[j] == disparities[i]) :
													(Math.abs(disparities[j] - disparities[i]) < max_diff);
										if (same) {
											group[j] = ngroup;
											group_members[nsame++] = j;
										}
									}
								if (nsame > 1) {
									members[ngroup]=new int[nsame];
									group[i] = ngroup;
									for (int j = 0; j < nsame; j++) {
										members[ngroup][j] = group_members[j];
									}
									ngroup++;
								}
							}
							if (ngroup > 0) {
								int y0 = (tile / tilesX) * transform_size;
								int x0 = (tile % tilesX) * transform_size;
								int indx0 = y0 * width + x0; 
								for (int ng =0; ng < ngroup; ng++) {
									for (int dy = 0; dy < transform_size; dy++) {
										int indx1 = indx0 + dy * width;
										for (int dx = 0; dx < transform_size; dx++) {
											int indx = indx1 + dx;
											boolean a = alpha_pix[members[ng][0]][indx];
											if (use_or) {
												for (int j = 1; j < members[ng].length; j++) {
													a |= alpha_pix[members[ng][j]][indx];
												}
											} else {
												for (int j = 1; j < members[ng].length; j++) {
													a &= alpha_pix[members[ng][j]][indx];
												}
											}
											for (int j = 0; j < members[ng].length; j++) {
												alpha_pix[members[ng][j]][indx] = a;
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
		return;
	}
	
	/**
	 * Generate bitmask of sensors that should be removed from the composite texture. Uses
	 * image offsets from TileTask array to get shift between textures rendered for different
	 * disparities. The source (pre-aberration) offsets directly are not used, just difference
	 * for the same sensors.
	 * Considering for being occluded all but strong FG tiles
	 *  
	 * @param channel_pixel_offsets per-slice, per tile (linescan order), per-sensor x,y offsets.
	 * @param alpha_pix         boolean "alpha" - true - opaque, false - transparent
	 * @param slice_disparities per-tile disparities ([slice][tile]).
	 * @param tile_keep         boolean map of kept tiles (tile_booleans[TILE_KEEP])
	 * @param tile_fg_strong    boolean map of strong FG tiles (tile_booleans[TILE_IS_FG_STRONG])
	 * @param tile_stitch       boolean map of stitch tiles (tile_booleans[TILE_STITCH]) - they have duplicates
	 * @param occlusion_frac    thershold for interpolating occlusion - fraction of BG tile being occluded
	 *                          to actually occlude
	 * @param occlusion_min_disp minimal FG BG separation to calculate occlusions on BG                         
	 * @param width             image width in pixels 
	 * @param transform_size    CLT conversion size. Always 8
	 * @return                  [nslice][pix] bit map of occluded sensors to be removed from sources of the
	 *                          combined textures. 
	 */
	public static int [][] getOccludedMap(
			final double [][][][] channel_pixel_offsets,
			final boolean [][]    alpha_pix,
			final double  [][]    slice_disparities,
			final boolean [][]    tile_keep,      // do not check occluded strong foreground
			final boolean [][]    tile_fg_strong, // do not check occluded strong foreground
			final boolean [][]    tile_stitch,    // do not process these - there are duplicates
			final double          occlusion_frac, // ratio of opaque pixel overlap to consider occlusion
			final double          occlusion_min_disp, 
			final int             width,
			final int             transform_size){
		final int num_slices = alpha_pix.length;
		final int img_size =   alpha_pix[0].length;
		final int height = img_size/width;
		final int tilesX = width/transform_size;
		final int tilesY = img_size/width/transform_size;
		final int tiles =  tilesX * tilesY;
		final int dbg_tile = -4123;
		final int dbg_slice = 0;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final int [][] occluded = new int [num_slices][img_size];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							if ((fnslice == dbg_slice) && (tile == dbg_tile )) {
								System.out.println("getNonOccludedMap().1 nslice="+fnslice+", tile="+tile);
							}
							if (    tile_keep[fnslice][tile] &&
									!tile_fg_strong[fnslice][tile] &&
									!tile_stitch[fnslice][tile]) { // ***
								
								double [][] offs_bg = channel_pixel_offsets[fnslice][tile];
								for (int ns = 0; ns < num_slices; ns++) {
									if ((ns != fnslice) &&
											tile_keep[ns][tile] &&
											!tile_stitch[fnslice][tile] && // ***
											(slice_disparities[ns][tile] > slice_disparities[fnslice][tile]) &&
											((slice_disparities[ns][tile] - slice_disparities[fnslice][tile]) > occlusion_min_disp )) {
										double [][] offs_fg = channel_pixel_offsets[ns][tile];
										double [][] pixel_offs = new double [offs_bg.length][2];
										for (int nsens = 0; nsens < pixel_offs.length; nsens++) {
											if (offs_bg[nsens] != null) { // to implement sensor mask later
												pixel_offs[nsens][0] = offs_bg[nsens][0] - offs_fg[nsens][0];
												pixel_offs[nsens][1] = offs_bg[nsens][1] - offs_fg[nsens][1];
											}
										}
										int tileX = tile % tilesX; 
										int tileY = tile / tilesX;
										for (int dy = 0; dy < transform_size; dy++) {
											int py0 = tileY * transform_size + dy;
											for (int dx = 0; dx < transform_size; dx++) {
												int px0 = tileX * transform_size + dx;
												int occluded_mask = 0;
												for (int nsens = 0; nsens < pixel_offs.length; nsens++) if (offs_bg[nsens] != null) {
													double px = px0 + pixel_offs[nsens][0]; 
													double py = py0 + pixel_offs[nsens][1];
													if ((px >= 0) && (px < (width - 1)) && (py >= 0) && (py < (height - 1))) {
														int ipx = (int) Math.floor(px);
														int ipy = (int) Math.floor(py);
														int indx_fg = ipx + ipy*width;
														boolean occl_any = 
																alpha_pix[ns][indx_fg] ||
																alpha_pix[ns][indx_fg + 1] ||
																alpha_pix[ns][indx_fg + width] ||
																alpha_pix[ns][indx_fg + width + 1];
														boolean occl_all = 
																alpha_pix[ns][indx_fg] &&
																alpha_pix[ns][indx_fg + 1] &&
																alpha_pix[ns][indx_fg + width] &&
																alpha_pix[ns][indx_fg + width + 1];
														if (occl_all) {
															occluded_mask |= (1 << nsens);
														} else {
															if (occl_any) {
																double fx =  px - ipx;
																double fy =  py - ipy;
																double d = 0;
																if (alpha_pix[ns][indx_fg])             d += (1 - fx) * (1- fy);
																if (alpha_pix[ns][indx_fg + 1])         d += (    fx) * (1- fy);
																if (alpha_pix[ns][indx_fg + width])     d += (1 - fx) * (   fy);
																if (alpha_pix[ns][indx_fg + width + 1]) d += (    fx) * (   fy);
																if (d >= occlusion_frac) {
																	occluded_mask |= (1 << nsens);
																}
															}
														}
													}
												} // for (int nsens = 0; nsens < pixel_offs.length; nsens++) {
												int indx = (((tileY * width + tileX) * transform_size) + dy * width) + dx;
												occluded[fnslice][indx] |= occluded_mask;
											}
										}
									} // if ((ns != fnslice) && ...
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		// duplicate for stitch tiles
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							if (tile_keep[fnslice][tile] && !tile_fg_strong[fnslice][tile] && tile_stitch[fnslice][tile]) {
								for (int ns = 0; ns < num_slices; ns++) {
									if ((ns != fnslice) && // other layer with same disparity and non-stitch (probably stitched) 
											tile_keep[ns][tile] &&
											!tile_stitch[fnslice][tile] &&
											(slice_disparities[ns][tile] == slice_disparities[fnslice][tile])) {
										int tileX = tile % tilesX; 
										int tileY = tile / tilesX;
										for (int dy = 0; dy < transform_size; dy++) {
											int indx0 = (tileY * transform_size + dy) * width + tileX * transform_size;
											// use or? both ways
											System.arraycopy(
													occluded[ns],
													indx0,
													occluded[fnslice],
													indx0,
													transform_size);
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
		return occluded;
		// for debug - display number of bits from bit_mask
	}
	
	public static boolean [][] getTrimTiles(
			boolean [][] trim_pix,
			final int    width,
			final int    transform_size){
		final int num_slices =    trim_pix.length;
		final int img_size =      trim_pix[0].length;
		final int tilesX =        width/transform_size;
		final int tilesY =        img_size/width/transform_size;
		final int tiles =         tilesX * tilesY;
		final boolean [][] trim_tiles = new boolean [num_slices][tiles];
		final Thread[] threads =  ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai =  new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			ai.set(0);
			final boolean [] trim_this = trim_pix[nslice];
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							int tileY = tile / tilesX;
							int tileX = tile % tilesX;
							int indx = (tileY * width + tileX) * transform_size;
							search_pix:
							for (int dy = 0; dy < transform_size; dy++) {
								for (int dx = 0; dx < transform_size; dx++) {
									if (trim_this[indx++]) {
										trim_tiles[fnslice][tile] = true;
										break search_pix;
									}
								}
								indx += width - transform_size;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return trim_tiles;
	}
	
	
	public static int updateFgAlpha(
			final double [][][][] channel_pixel_offsets,
			final double  [][]    textures,
			final boolean [][]    alpha_pix,      // will be updated
//			final double  [][]    combo_texture,
			final double  [][][]  sensor_texture,
			final int     [][]    occluded_map,   // bitmap of blocked by FG sensors
			final int             min_sensors,    // minimal number of sensors visible from the FG pixel
			final double  [][]    slice_disparities,
			final boolean [][]    tile_keep,      // tiles that have at least one pixel 
			final boolean [][]    tile_stitch,      // tiles that have at least one pixel 
			final boolean [][]    trim_tiles,     // tiles that have at least one pixel 
			final boolean [][]    trim_pix,       // pixels that may be trimmed
			final boolean [][]    transparent,    // definitely transparent
			final boolean [][]    opaque,         // definitely opaque
			final boolean         en_cut,         // enable change FG pixel to transparent from opaque
			final boolean         en_patch,       // enable change FG pixel to opaque from transparent
			final double          min_disp_diff,  // do not consider obscuring too close BG (1 pix or more?)
			// other parameters
			final double          weight_neib,    // weight of same neighbors
			final double          weight_bg,      // weight of BG cost relative to the FG one
//			final double          weight_bg2,     // fraction of BG variance cost (1-weight_bg2) - the BG true one
			final double          best_dir_frac,  // for BG - use this fraction of all sensors in the best direction
			final double          cost_min,       // minimal absolute value of the total cost to make changes
			// debug arrays
			final double [][][]   debug_costs,    // if not null, should be double [nslices][] - will return costs/NaN
			final int    [][]     debug_stats,    // if not null, should be int [nslices][] - will return number of added/removed per slice
			final int             width,
			final int             transform_size){
		final int min_sensors_bg = min_sensors; // maybe reduce? *=best_dir_frac?
		final int num_slices =    alpha_pix.length;
		final int img_size =      alpha_pix[0].length;
		final int height =        img_size/width;
		final int tilesX =        width/transform_size;
		final int tilesY =        img_size/width/transform_size;
		final int tiles =         tilesX * tilesY;
		final int dbg_tile =      -1; // 4123;
		final int dbg_slice =     0;
		final Thread[] threads =  ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai =  new AtomicInteger(0);
		final AtomicInteger aplus =   new AtomicInteger(0); // number of added opaque pixels
		final AtomicInteger aminus =  new AtomicInteger(0); // number of removed opaque pixels
		final TileNeibs pn = new TileNeibs(width,height);
		int num_modified_pixels = 0;
		final int dbg_pix = -168170;
		final boolean [][] new_alpha = new boolean[num_slices][img_size];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			int fnslice = nslice;
			System.arraycopy(alpha_pix[fnslice], 0, new_alpha[fnslice], 0, img_size);
			if (debug_costs != null) {
				debug_costs[fnslice] = new double [3][img_size]; // {cost, cost_fg, cost_bg}
				for (int i = 0; i < debug_costs[fnslice].length; i++) {
					Arrays.fill(debug_costs[fnslice][i], Double.NaN);
				}
			}
			ai.set(0);
			aplus.set(0);
			aminus.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [][][] bg_value = new double [transform_size][transform_size][]; 
						double [][][] bg_disparity = new double [transform_size][transform_size][]; 
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) {
							if ((fnslice == dbg_slice) && (tile == dbg_tile )) {
								System.out.println("updateFgAlpha().1 nslice="+fnslice+", tile="+tile);
							}
							if (     trim_tiles [fnslice][tile] &&
									!tile_stitch[fnslice][tile]) {
								int tileX = tile % tilesX; 
								int tileY = tile / tilesX;
								int pix0 = (tileY  * width + tileX) * transform_size; 
								double [][] offs_fg = channel_pixel_offsets[fnslice][tile];
								int num_sens = offs_fg.length;
								int best_dir_number = (int) Math.round (best_dir_frac * num_sens);
								boolean valid_bg = false;
								for (int ns = 0; ns < num_slices; ns++) {
									if ((ns != fnslice) &&
											tile_keep[ns][tile] &&
											!tile_stitch[fnslice][tile] && // ***
											(slice_disparities[ns][tile] < slice_disparities[fnslice][tile]) &&
											((slice_disparities[fnslice][tile] - slice_disparities[ns][tile]) > min_disp_diff )) {
										double [][] offs_bg = channel_pixel_offsets[ns][tile];
										double [][] pixel_offs = new double [num_sens][2];
										for (int nsens = 0; nsens < num_sens; nsens++) {
											if (offs_bg[nsens] != null) { // to implement sensor mask later
												pixel_offs[nsens][0] = offs_fg[nsens][0] - offs_bg[nsens][0];
												pixel_offs[nsens][1] = offs_fg[nsens][1] - offs_bg[nsens][1];
											}
										}
										for (int dy = 0; dy < transform_size; dy++) {
											int py0 = tileY * transform_size + dy;
											int pix1 = pix0 + dy * width;
											for (int dx = 0; dx < transform_size; dx++) {
												int px0 = tileX * transform_size + dx;
												int pix = pix1 + dx;
												if (pix==dbg_pix) {
													System.out.println("updateFgAlpha().1 pix="+pix+", ns="+ns+", dx="+dx+", dy="+dy+
															", disp_fg="+slice_disparities[fnslice][tile]+", disp_bg="+slice_disparities[ns][tile]);
												}
												if (trim_pix[fnslice][pix]) { // assign for all trim_pix 
													if (!transparent[fnslice][pix] && !opaque[fnslice][pix]) {
														if ((alpha_pix[fnslice][pix] && en_patch) || (!alpha_pix[fnslice][pix] && en_cut)) {
															// calculate costs
															// maybe multiple backgrounds? Then combine them all
															// each sensor - single BG - common array of 16?

															// cost for FG - average w/o center, possibly tilt
															// consider spread normalize to sigma?
															// how to normalize BG error
															// Or do not normalize at all - compare absolute values?
															if (!valid_bg) { // lazy initialization
																for (int i = 0; i < transform_size; i++) {
																	Arrays.fill(bg_disparity[i], null);
																}
																valid_bg = true;
															}
															if (bg_disparity[dy][dx] == null) {
																bg_disparity[dy][dx] = new double [num_sens];
																Arrays.fill(bg_disparity[dy][dx], Double.NaN);
															}
															if (bg_value[dy][dx] == null) { // do not need to initialize
																bg_value[dy][dx] = new double [num_sens];
															}
															for (int nsens = 0; nsens < pixel_offs.length; nsens++) if (offs_bg[nsens] != null) {
																// corresponding BG pixels
																double bgx = px0 + pixel_offs[nsens][0]; 
																double bgy = py0 + pixel_offs[nsens][1];
																int ibgx = (int) Math.round(bgx); // here just center
																int ibgy = (int) Math.round(bgy);
																int bg_pix = ibgx + ibgy * width;
																if (alpha_pix[ns][bg_pix] && !Double.isNaN(textures[ns][bg_pix])) {
																	if (!(bg_disparity[dy][dx][nsens] >= slice_disparities[ns][tile])) { // was NaN -> true
																		bg_disparity[dy][dx][nsens] = slice_disparities[ns][tile];
																		bg_value[dy][dx][nsens] = textures[ns][bg_pix];
																		if (pix==dbg_pix) {
																			System.out.println(String.format(
																					"%2d: bg_pix=%6d ibgx=%3d ibgy=%3d bg_value=%8.2f",
																					nsens, bg_pix, ibgx, ibgy, bg_value[dy][dx][nsens]));
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
								} // for (int ns = 0; ns < num_slices; ns++) {
								// now consider if (!transparent[fnslice][indx] && !opaque[fnslice][indx]),
								// use bg_value[][][], bg_disparity[][][] to calculate bg weighths,
								// calculate FG weights (same as VAR_INTER)
								// add num neibs - 4 weight and make decisions
								for (int dy = 0; dy < transform_size; dy++) {
									int pix1 = pix0 + dy * width;
									for (int dx = 0; dx < transform_size; dx++) {
										int pix = pix1 + dx;
										if (pix==dbg_pix) {
											System.out.println("updateFgAlpha().2 pix="+pix);
										}
										if (trim_pix[fnslice][pix]) { // assign for all trim_pix
											boolean new_transparent = false;
											boolean new_opaque = false;
											if (!transparent[fnslice][pix] && !opaque[fnslice][pix]) {
												if ((alpha_pix[fnslice][pix] && en_patch) || (!alpha_pix[fnslice][pix] && en_cut)) {
													// calculate number of sensors, visible from this FG pixel
													int smask = occluded_map[fnslice][pix];
													int num_fg = 0;
													double s_fg=0, s2_fg = 0;
													for (int nsens = 0; nsens < num_sens; nsens++ ) {
														if ((smask & (1 << nsens)) == 0) {
															double d = sensor_texture[fnslice][nsens][pix];
															s_fg += d;
															s2_fg += d * d;
															num_fg++;
														}
													}
													
													int num_bg=0;
													double s2_bg = 0;
													for (int nsens = 0; nsens < (best_dir_number-1); nsens++ ) {
														if (pix==dbg_pix) {
															System.out.println(String.format(
																	"%2d: fg_value= %8.2f bg_value=%8.2f diff=%8.2f",
																	nsens,
																	sensor_texture[fnslice][nsens][pix],
																	bg_value[dy][dx][nsens],
																	sensor_texture[fnslice][nsens][pix] - bg_value[dy][dx][nsens]));
														}
														if (((smask & (1 << nsens)) == 0) && (bg_disparity[dy][dx] != null) && !Double.isNaN(bg_disparity[dy][dx][nsens])) {
															double db = sensor_texture[fnslice][nsens][pix] - bg_value[dy][dx][nsens];
															s2_bg += db * db;
															num_bg ++;
														}
													}													
													double best_cost_bg = Double.NaN;
													for (int i = 0; i < num_sens; i++ ) {
														int nsens_plus =  (best_dir_number + i - 1) % num_sens;
														int nsens_minus = i;
														if (pix==dbg_pix) {
															System.out.println(String.format(
																	"%2d: fg_value= %8.2f bg_value=%8.2f diff=%8.2f",
																	nsens_plus,
																	sensor_texture[fnslice][nsens_plus][pix],
																	bg_value[dy][dx][nsens_plus],
																	sensor_texture[fnslice][nsens_plus][pix] - bg_value[dy][dx][nsens_plus]));
														}
														if (((smask & (1 << nsens_plus)) == 0) &&
																(bg_disparity[dy][dx] != null) &&
																!Double.isNaN(bg_disparity[dy][dx][nsens_plus])) {
															double db = sensor_texture[fnslice][nsens_plus][pix] - bg_value[dy][dx][nsens_plus];
															s2_bg += db * db;
															num_bg ++;
														}
														if (num_bg >= min_sensors_bg) {
															double avg2_bg = s2_bg/num_bg;
															double cost_bg= Math.sqrt(avg2_bg);
															if (!(cost_bg > best_cost_bg)) {
																best_cost_bg = cost_bg;
															}
															if (pix==dbg_pix) {
																System.out.print(String.format(
																		"avg2_bg= %8.2f cost_bg=%8.2f -> ",
																		avg2_bg, cost_bg));
															}
														}
														if (i < (num_sens -1)) {
															if (((smask & (1 << nsens_minus)) == 0) &&
																	(bg_disparity[dy][dx] != null) &&
																	!Double.isNaN(bg_disparity[dy][dx][nsens_minus])) {
																double db = sensor_texture[fnslice][nsens_minus][pix] - bg_value[dy][dx][nsens_minus];
																s2_bg -= db * db;
																num_bg --;
															}
														}
													}
													// calculate costs
													// maybe multiple backgrounds? Then combine them all
													// each sensor - single BG - common array of 16?
													// cost for FG - average w/o center, possibly tilt
													// consider spread normalize to sigma?
													// how to normalize BG error
													// Or do not normalize at all - compare absolute values?
													if (num_fg >= min_sensors) {//  (do not touch if less)
														double avg_fg = s_fg/num_fg;
														double avg2_fg = s2_fg/num_fg;
														double cost_fg = Math.sqrt(avg2_fg-avg_fg*avg_fg);
														double cost_bg = best_cost_bg;
														// calculate number of opaque neighbors
														int n_opaque = 0, n_neibs=0;
														for (int dir = 0; dir < TileNeibs.DIRS; dir++) {
															int pix_n = pn.getNeibIndex(pix, dir);
															if (pix_n >=0) {
																if (alpha_pix[fnslice][pix_n]) {
																	n_opaque++;
																}
															}
															n_neibs++;
														}
														// positive for more opaque, negative - for more transparent
														double cost_neibs = weight_neib * (n_opaque - 0.5* n_neibs);
														double cost = Double.NaN;
														if (Double.isNaN(cost_bg)) {
															new_opaque = true;
														} else {
															// cost > 0 -> opaque, cost < 0 -> transparent 
															cost = weight_bg*cost_bg - cost_fg + cost_neibs;
															if (Math.abs(cost) > cost_min) {
																new_opaque =      cost > 0;
																new_transparent = cost < 0;
															}
														}
														if (debug_costs != null) {
															debug_costs[fnslice][0][pix] = cost;
															debug_costs[fnslice][1][pix] = cost_fg;
															debug_costs[fnslice][2][pix] = cost_bg;
															if (pix==dbg_pix) {
																System.out.println(String.format(
																		"cost= %8.2f cost_neibs=%8.2f cost_fg=%8.2f cost_bg=%8.2f",
																		cost, cost_neibs, cost_fg, cost_bg));
															}
														}
														
													} // if (num_vis >= min_sensors) {
												}
											} else { // if (!transparent[fnslice][indx] && !opaque[fnslice][indx]) {
												if (transparent[fnslice][pix]) {
													new_transparent = true;
												} else if (opaque[fnslice][pix]) {
													new_opaque = true;
												}
											} // if (!transparent[fnslice][pix] && !opaque[fnslice][pix])
											if (new_opaque) {
												if (!alpha_pix[fnslice][pix] && en_patch) {
													new_alpha[fnslice][pix] = true;
													aplus.getAndIncrement();
												}
											} else if (new_transparent){
												if (alpha_pix[fnslice][pix] && en_cut) {
													new_alpha[fnslice][pix] = false;
													aminus.getAndIncrement();
												}
											}
										}
									} // for (int dx = 0; dx < transform_size; dx++) {
								} // for (int dy = 0; dy < transform_size; dy++)
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (debug_stats != null) {
				debug_stats[nslice] = new int[] {aplus.get(), aminus.get()};
			}
			num_modified_pixels += aplus.get() + aminus.get();
		}
		// replace boolean alphas with the new ones.
		for (int nslice = 0; nslice < num_slices; nslice++) {
			alpha_pix[nslice] = new_alpha[nslice];
		}
		// consider using such method without preliminary methods with using
		// analog (semi-transparent) alpha that finally stick to 0/1
		return num_modified_pixels;
	}

	public static int [][] enhanceOccludedMap(
			final int [][] occluded_map_in,
			final double   keep_frac,
			final int      min_sensors,
			final int      num_sens,
			final int      dbg_indx){
		// TODO: for 16 sensors it will be faster to create int[65536] and recode by table
		boolean create_table = occluded_map_in == null;
		final int [][] occluded_map = create_table? new int [1][1 << num_sens] : occluded_map_in;
		if (create_table) {
			for (int i = 0; i < occluded_map[0].length; i++) {
				occluded_map[0][i] = i;
			}
		}
		final int num_slices = occluded_map.length;
		final int img_size = occluded_map[0].length;
		final int keep_sens = (int) Math.round(num_sens * keep_frac);
		final int [][] map_enh = new int [num_slices][img_size];
		final int mask = (1 << num_sens) -1; // will not work for >=31
		final int dbg_val = dbg_indx; //1018;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < img_size; pix = ai.getAndIncrement()) {
							if (pix ==  dbg_val) {
								System.out.println("enhanceOccludedMap(): pix="+pix);
							}
							if ((occluded_map[fnslice][pix] == 0) || (occluded_map[fnslice][pix] & mask) == mask ){
								map_enh[fnslice][pix] = occluded_map[fnslice][pix];		
							} else {
								int msk = occluded_map[fnslice][pix]; // bit = 0 - exists
								int []  sp0 = new int[num_sens],
										sm0 = new int[num_sens],
										spx = new int[num_sens],
										spx2= new int[num_sens],
										smx = new int[num_sens],
										smx2= new int[num_sens];
								for (int i = 0; i < num_sens; i++) {
									boolean sens_on =  (msk & (1 << i)) == 0;
									if (i < keep_sens) {
										if (sens_on) {
											sp0[0]++;
											spx[0] += i;
											spx2[0] +=i * i;
										}
									} else {
										if (!sens_on) {
											sm0[0]++;
											smx[0] +=  i;
											smx2[0] += i * i;
										}
									}
								}
								int best_num_used = sp0[0];
								for (int i = 1; i < num_sens; i++) {
									sp0[i] =  sp0[i-1];
									spx[i] =  spx[i-1];
									spx2[i] = spx2[i-1];
									sm0[i] =  sm0[i-1];
									smx[i] =  smx[i-1];
									smx2[i] = smx2[i-1];
									if ((msk & (1 << (i - 1))) == 0) {
										sp0[i]--;
										int i1 = i -1;
										spx[i] -= i1;
										spx2[i] -= i1 * i1;
									} else { // for minus count missing sensors
										sm0[i]++;
										int i1 = i + num_sens -1;
										smx[i] += i1;
										smx2[i] += i1*i1;
									}
									if ((msk & (1 << ((i + keep_sens- 1) % num_sens))) == 0) {
										int i1 = i + keep_sens- 1 ;
										sp0[i]++;
										spx[i] += i1;
										spx2[i] += i1*i1;
									} else {
										int i1 = i + keep_sens- 1 ;
										sm0[i]--;
										smx[i] -= i1;
										smx2[i] -= i1*i1;
									}
									if (sp0[i] > best_num_used) {
										best_num_used = sp0[i]; 
									}
								}
								if (best_num_used <  min_sensors) {
									map_enh[fnslice][pix] = mask; // all sensors disabled
								} else {
									boolean [] candidates = new boolean [num_sens];
									double [] pointed2 = new double [num_sens];
									int num_prev_candidates = 0;
									double best_pointed2 = Double.NaN;
									int best_index = -1;
									// find (may be multiple) best number of used sensors
									for (int i = 0; i < num_sens; i++) {
										if (sp0[i] == best_num_used) {
											candidates[i] = true;
											num_prev_candidates ++;
											double avg =  (1.0 * spx[i]) / sp0[i];
											double avg2 = (1.0 * spx2[i]) / sp0[i];
											double center = i + 0.5 * (keep_sens- 1);
											double offset = (avg - center);
											double l2 = avg2 - avg*avg + offset * offset;
											if (!(l2 > best_pointed2)) {
												best_pointed2 = l2;
												best_index = i;
											}
											pointed2[i] = l2;
										}
									}
									double best_pointed2m = Double.NaN;
									if (num_prev_candidates > 1) {
										for (int i = 0; i < num_sens; i++) if (candidates[i]) {
											if (pointed2[i] > best_pointed2) {
												candidates[i] = false;
												num_prev_candidates --;
											} else {
												double avg =  (1.0 *smx[i]) / sm0[i]; // can not be zero
												double avg2 = (1.0 * smx2[i]) / sm0[i];
												double center = i + keep_sens + 0.5 * (num_sens- keep_sens- 1); // remaining
												double offset = (avg - center);
												double l2 = avg2 - avg*avg + offset * offset;
												if (!(l2 > best_pointed2m)) {
													best_pointed2m = l2;
													best_index = i;
												}
												pointed2[i] = l2;
											}
										}
									}
									if (num_prev_candidates > 1) {
										for (int i = 0; i < num_sens; i++) if (candidates[i]) {
											if (pointed2[i] > best_pointed2m) {
												candidates[i] = false;
												num_prev_candidates --;
											} else {
												best_index = i;
											}
										}									
									}
									for (int i = 0; i < (num_sens- keep_sens); i++) {
										int i1 = (best_index + i + keep_sens) % num_sens;
										msk |= 1 << i1;
									}
									map_enh[fnslice][pix] = msk;
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		
		String sin = ""; // ++++++++--------";
		while(testDirectionMap (
				sin, // String sin,
				map_enh[0])); //int [] map)
		
		return map_enh;
	}
	
	public static int [][] enhanceOccludedMap(
			final int [][] occluded_map,
			final int []   map){
		final int num_slices = occluded_map.length;
		final int img_size = occluded_map[0].length;
		final int [][] map_enh = new int [num_slices][img_size];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < img_size; pix = ai.getAndIncrement()) {
						map_enh[fnslice][pix] = map[occluded_map[fnslice][pix]];
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return map_enh;
	}
	
	public static boolean testDirectionMap (
			String sin,
			int [] map) {
		if (sin.isEmpty()) {
			return false;
		}
		int [] d = new int[2];
		for (int i = sin.length()-1; i >= 0; i--) {
			d[0] = (d[0] << 1) + ((sin.charAt(i) == '+')? 0:1); // 1 - occluded
		}
		d[1] = map[d[0]];
		for (int n = 0; n < d.length; n++) {
			for (int i = 15; i >=0; i--) {
				System.out.print((((d[n] >> i) & 1) == 0) ? " +":" -" );
			}
			System.out.println();
		}
		return true;
	}
	
	public static double [][] debugOccludedMap(
			final int     [][]   occluded_map){
		final int num_slices =  occluded_map.length;
		final int img_size =    occluded_map[0].length;
		double [][] dbg_map = new double [num_slices][img_size];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			for (int pix = 0; pix < img_size; pix++) {
				if (occluded_map[nslice][pix] != 0) {
					int n = 0;
					for (int d = occluded_map[nslice][pix]; d != 0;  d >>= 1) {
						if ((d & 1) != 0) {
							n++;
						}
					}
					dbg_map[nslice][pix] = n;
				}
			}
		}		
		return dbg_map;
	}
	

	
	
	public static double [][] combineTexturesWithOcclusions(
			final double  [][][] sensor_texture,
			final double  [][]   combo_texture,
			final int     [][]   occluded_map){
		final int num_slices =  sensor_texture.length;
		final int img_size =    combo_texture[0].length;
		final int num_sensors = sensor_texture[0].length;
		final double [][] occluded_texture = new double [num_slices][img_size];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < img_size; pix = ai.getAndIncrement()) {
							if (occluded_map[fnslice][pix] == 0) {
								occluded_texture[fnslice][pix] = combo_texture[fnslice][pix];
							} else {
								int num_used_sensors = 0;
								int msk = occluded_map[fnslice][pix];
								double s = 0.0;
								for (int nsens = 0; nsens < num_sensors; nsens++ ) {
									if ((msk & (1 << nsens)) == 0) {
										s += sensor_texture[fnslice][nsens][pix];
										num_used_sensors++;
									}
								}
								if (num_used_sensors > 0) {
									s /= num_used_sensors;
								} else {
									s = Double.NaN;
								}
								occluded_texture[fnslice][pix] = s;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return occluded_texture;
	}
	
	public static double [][] fillOcclusionsNaN(
			final double  [][]   combo_texture,
			final double  [][]   combo_occluded_texture,
			final int            grow,
			final int            num_passes,
			final double         max_change,
			final int            width){
		final int num_slices =        combo_texture.length;
		final int img_size =          combo_texture[0].length;
		final double diagonal_weight = 0.7;
		final boolean [] prohibit =   new boolean [img_size];
		final double [][] filled_occluded = new double[num_slices][];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int pix = ai.getAndIncrement(); pix < img_size; pix = ai.getAndIncrement()) {
							prohibit[pix] = Double.isNaN(combo_texture[fnslice][pix]);
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			filled_occluded[fnslice] = TileProcessor.fillNaNs( // multithreaded
					combo_occluded_texture[fnslice],   // final double [] data,
					prohibit,                          // final boolean [] prohibit,
					width,                             // int       width,
					grow,                              // 16,           // final int grow,
					diagonal_weight,                   // double    diagonal_weight, // relative to ortho
					num_passes,                        // int       num_passes,
					max_change,                        // final double     max_rchange, //  = 0.01
					THREADS_MAX);                      // final int threadsMax)      // maximal number of threads to 
		}
		return filled_occluded;
	}
	
	
	
	
	public static double  [][] processBgOcclusions(
			final double  [][][] sensor_texture,
			final double  [][]   combo_texture,
			final boolean [][]   selections,
			final boolean [][]   is_fg_tile,
			final boolean [][]   has_bg_tile,
			final double  [][]   vars_inter,
			double               dir_radius,
			final int            width,
			final int            transform_size,
			final double         try_dir_var,     // 20.0; // try directional if the intersensor variance exceeds this value
			final int            dir_num_start,   //  7;   // start with this number of consecutive sensors			
			final int            dir_num_restart, //  5;   // restart (from best direction) with this number of consecutive sensors
			final double         dir_var_max,     // 15.0; // do not add more sensors if the variance would exceed this
			final double         dir_worsen_rel,  //  0.15;// add more sensors until variance grows by this relative
			final double [][][]  dbg_out // [5][][]
			) {
		final int num_slices =  selections.length;
		final int img_size =    selections[0].length;
		final int num_sensors = sensor_texture[0].length;
		final int tilesX =      width / transform_size;
		final double [][] out_textures = new double [num_slices][];
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final TileNeibs pn =     new TileNeibs(width, img_size/width);
		final int idir_radius = (int) Math.floor(dir_radius);
		final double [][] dir_weights = new double [idir_radius+1][idir_radius+1];
		for (int i = 0; i < dir_weights.length; i++) {
			for (int j = 0; j < dir_weights[i].length; j++) {
				dir_weights[i][j] = Math.cos(0.5*Math.PI*i/dir_radius) * Math.cos(0.5*Math.PI*j/dir_radius);
			}
		}
		if (dbg_out != null) {
			for (int i = 0; i < dbg_out.length; i++) {
				dbg_out[i] = new double [num_slices][];
			}
		}

		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			out_textures[fnslice] = combo_texture[fnslice].clone();
			final double [] vars_dir_initial = new double [img_size];
			final double [] dirs_initial =     new double [img_size];
			final double [] vars_dir_final =   new double [img_size];
			final double [] dirs_final =       new double [img_size];
			final double [] dirs_len =         new double [img_size];
			System.arraycopy(vars_inter[fnslice], 0, vars_dir_initial, 0, img_size);
			Arrays.fill(dirs_initial,     Double.NaN);
			Arrays.fill(dirs_final,       Double.NaN);
			Arrays.fill(dirs_len,         Double.NaN);
			ai.set(0);
			System.arraycopy(vars_inter[fnslice], 0, vars_dir_final, 0, img_size);
			final int dir_samples = (2 * idir_radius + 1) * (2 * idir_radius + 1);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] dirvar = new double [num_sensors];
						double [] sw =     new double[dir_samples];
						double [] swd =    new double[dir_samples];
						double [] swd2 =   new double[dir_samples];
						for (int cindx = ai.getAndIncrement(); cindx < img_size; cindx = ai.getAndIncrement())
							// is_fg_tile[fnslice][tile]
							if (selections[fnslice][cindx])  {
								int y0 = cindx / width;
								int x0 = cindx % width;
								int tileX = x0 / transform_size;
								int tileY = y0 / transform_size;
								int tile = tileX + tilesX * tileY;
								// only trim if nothing obscures this and has some BG 
								if (!(is_fg_tile[fnslice][tile] && has_bg_tile[fnslice][tile]) &&
										(vars_inter[fnslice][cindx] > try_dir_var)) { // try obscuring by others
									Arrays.fill(sw,  0.0);
									Arrays.fill(swd, 0.0);
									Arrays.fill(swd2,0.0);
									for (int dvy = -idir_radius; dvy <= idir_radius; dvy++) {
										for (int dvx = -idir_radius; dvx <= idir_radius; dvx++) {
											int pindx = pn.getIndex(x0+dvx, y0+dvy);
											if ((pindx >= 0) && selections[fnslice][cindx]) {
												int vindx = (dvy + idir_radius) * (2 * idir_radius + 1) + (dvx + idir_radius);
												for (int nsens = 0; nsens < dir_num_start; nsens++) {
													double w = 1.0;
													double d = sensor_texture[fnslice][nsens][pindx];
													sw[vindx] +=   w;
													swd[vindx] +=  w * d;
													swd2[vindx] += w * d * d;
												}
											}
										}
									}
									double svw = 0.0, svwd = 0.0;
									for (int start_dir = 0; start_dir < num_sensors; start_dir++) {
										for (int dvy = -idir_radius; dvy <= idir_radius; dvy++) {
											for (int dvx = -idir_radius; dvx <= idir_radius; dvx++) {
												int vindx = (dvy + idir_radius) * (2 * idir_radius + 1) + (dvx + idir_radius);
												if (sw[vindx] > 0.0) {
													double avg =  swd[vindx]/sw[vindx];
													double avg2 = swd2[vindx]/sw[vindx];
													double v = Math.sqrt(avg2-avg*avg);
													double w = dir_weights[Math.abs(dvy)][Math.abs(dvx)];
													svw +=  w;
													svwd += w * v;
												}
											}
										}
										dirvar[start_dir] = svwd/svw; // weighted
										if (start_dir == (num_sensors - 1)) {
											break;
										}
										// remove previous values pointed by start_dir
										// and add new values pointed by start_dir+dir_num_start
										for (int dvy = -idir_radius; dvy <= idir_radius; dvy++) {
											for (int dvx = -idir_radius; dvx <= idir_radius; dvx++) {
												int pindx = pn.getIndex(x0+dvx, y0+dvy);
												if ((pindx >= 0) && selections[fnslice][cindx]) {
													int vindx = (dvy + idir_radius) * (2 * idir_radius + 1) + (dvx + idir_radius);
													double w0 = 1.0;
													double d0 = sensor_texture[fnslice][start_dir][pindx];
													sw[vindx] -=   w0;
													swd[vindx] -=  w0 * d0;
													swd2[vindx] -= w0 * d0 * d0;
													int dir1 = (start_dir + dir_num_start) % num_sensors;
													double w1 = 1.0;
													sw[vindx] +=   w1;
													double d1 = sensor_texture[fnslice][dir1][pindx];
													swd[vindx] +=  w1 * d1;
													swd2[vindx] += w1 * d1 * d1;
												}
											}
										}
									}
									// find best direction
									int best_dir = 0;
									for (int dir = 1; dir < num_sensors; dir++) {
										if (dirvar[dir] < dirvar[best_dir]) {
											best_dir = dir;
										}
									}
									double center_dir = best_dir + 0.5 * (dir_num_start - 1);
									if (center_dir >= num_sensors) {
										center_dir -= num_sensors;
									}
									vars_dir_initial[cindx] = dirvar[best_dir];
									dirs_initial[cindx] =     center_dir;

									// recalculate variance for the same direction but maybe different number of sensors
									// dir_num_restart
									// Use only center pixel, do not average variance
									int dir_start = (int) Math.round(center_dir - 0.5 * (dir_num_restart - 1));
									if (dir_start < 0) {
										dir_start += num_sensors;
									}
									double scw = 0.0, scwd = 0.0, scwd2 = 0.0;
									for (int i = 0; i < dir_num_restart; i++) {
										int dir = dir_start +i;
										if (dir >= num_sensors) {
											dir -= num_sensors;
										}
										double w = 1.0;
										double d = sensor_texture[fnslice][dir][cindx];
										scw += w;
										scwd += w * d;
										scwd2 += w * d*d;
									}
									double avg =     scwd/scw;
									double avg2 =    scwd2/scw;
									double dir_var = Math.sqrt(avg2-avg*avg);
									double max_var = Math.min(dir_var_max, dir_var*(1.0+dir_worsen_rel));
									// add sensors from both ends (which is better) while variance stays below max_var
									int dir_len ;
									for (dir_len = dir_num_start; dir_len < num_sensors; dir_len++) {
										// try CCW
										int nsens0 = dir_start - 1; //) % num_sensors;
										if (nsens0 < 0) {
											nsens0 += num_sensors;
										}
										double w = 1.0;
										double d = sensor_texture[fnslice][nsens0][cindx];
										double sw_0 =   scw + w;
										double swd_0 =  scwd + w*d;
										double swd2_0 = scwd2 + w*d*d;
										avg =  swd_0 / sw_0;
										avg2 = swd2_0 / sw_0;
										double var_0 = Math.sqrt(avg2 - avg * avg);
										// try CW
										int nsens1 = dir_start + dir_len;
										if (nsens1 >=num_sensors) {
											nsens1 -= num_sensors ;
										}
										w = 1.0;
										d = sensor_texture[fnslice][nsens1][cindx];
										double sw_1 =   scw + w;
										double swd_1 =  scwd + w*d;
										double swd2_1 = scwd2 + w*d*d;
										avg =  swd_1 / sw_1;
										avg2 = swd2_1 / sw_1;
										double var_1 = Math.sqrt(avg2 - avg * avg);
										if (Math.min(var_0, var_1) > max_var) {
											break;
										}
										if (var_0 < var_1) {
											dir_start = nsens0;
											scw = sw_0;
											scwd = swd_0;
											scwd2 = swd2_0;
											dir_var = var_0;
										} else {
											// dir_start stays the same
											scw = sw_1;
											scwd = swd_1;
											scwd2 = swd2_1;
											dir_var = var_1;
										}
									}
									double dir_avg = scwd/scw;
									double ddir = dir_start + 0.5 * (dir_len - 1);
									if (ddir >= num_sensors) {
										ddir -= num_sensors;
									}
									// fill arrays
									dirs_final[cindx] =            ddir;
									dirs_len[cindx] =              dir_len;
									vars_dir_final[cindx] =        dir_var; 
									out_textures[fnslice][cindx] =       dir_avg; 
								} // } else if (vars_inter[cindx] > try_dir_var) { // try obscuring by others
							}
					}
				};
			}
			ImageDtt.startAndJoin(threads);
			if (dbg_out != null) {
				dbg_out[0][fnslice] = vars_dir_initial; 
				dbg_out[1][fnslice] =   vars_dir_final; 
				dbg_out[2][fnslice] =     dirs_initial; 
				dbg_out[3][fnslice] =       dirs_final; 
				dbg_out[4][fnslice] =         dirs_len;
			}
		}
		return out_textures;
	}
	
	public static double [][] copyTexture (double [][] data){
		double [][] result = new double[data.length][];
		for (int i = 0; i < data.length; i++) {
			result[i] = data[i].clone();
		}
		return result;
	}
	public static boolean [][] copyTexture (boolean [][] data){
		boolean [][] result = new boolean[data.length][];
		for (int i = 0; i < data.length; i++) {
			result[i] = data[i].clone();
		}
		return result;
	}
	public static double [][] dbgBooleanTexture (boolean [][] data, double f, double t){
		double [][] result = new double[data.length][];
		for (int i = 0; i < data.length; i++) {
			result[i] = new double [data[i].length];
			for (int j=0; j < result[i].length; j++) {
				result[i][j] = data[i][j]? t : f; 
			}
		}
		return result;
	}

	public static double [][] dbgBooleanTexture (boolean [][] data0, boolean [][] data1, double v0, double v1, double v2, double v3){
		double [][] result = new double[data0.length][];
		for (int i = 0; i < data0.length; i++) {
			result[i] = new double [data0[i].length];
			for (int j=0; j < result[i].length; j++) {
				result[i][j] = data1[i][j]? (data0[i][j]?v3:v2): (data0[i][j]?v1:v0); 
			}
		}
		return result;
	}
	
	public static double [][] generateAlphaTemplates(
			final int transform_size,
			boolean debug){
		double [] fade1d = new double [transform_size];     // 8
		double [][] alpha8 = new double [8][transform_size * transform_size];
		double [][] alpha = new double [256][transform_size * transform_size];
		for (int i = 0; i < fade1d.length; i++) {
			fade1d[i] = 0.5 * (1 + Math.cos((0.0 + i) *Math.PI /transform_size)); //0.5
		}
		int tsm1 = transform_size -1;
		for (int idir = 0; idir < 8; idir++) {
			for (int y = 0; y < transform_size; y++) {
				for (int x = 0; x < transform_size; x++) {
					int indx = -1;
					switch (idir ) {
					case 0: indx =  tsm1 - y;     break;
					case 1: indx =  x - y;        break;
					case 2: indx =  x;            break;
					case 3: indx =  x + y - tsm1; break;
					case 4: indx =  y;            break;
					case 5: indx = -x + y;        break;
					case 6: indx =  tsm1 - x;     break;
					case 7: indx = -x - y + tsm1; break;
					}
					if (indx >= 0) {
						alpha8[idir][y * transform_size + x] = fade1d[indx];
					} else {
						alpha8[idir][y * transform_size + x] = 1.0;
					}
				}
			}
		}

		for (int imode = 0; imode < alpha.length; imode++) {
			int mode = (imode | imode << 8) & 0x1ff;
			Arrays.fill(alpha[imode], 1.0);
			for (int idir = 0; idir < alpha8.length; idir++) {
				if ((mode & (1 << idir)) == 0) {
					if (((idir & 1) != 0) && ((mode & (1 << (idir -1))) == 0)  && ((mode & (1 << (idir +1))) == 0)) {
						for (int i = 0; i < alpha8[idir].length; i++) {
							int i1 = alpha8[0].length - 1 -i;
							alpha[imode][i] = Math.min(alpha[imode][i],(1.0 - alpha8[idir][i1]));
						}
					} else {
						for (int i = 0; i < alpha8[idir].length; i++) {
							alpha[imode][i] = Math.min(alpha[imode][i], alpha8[idir][i]);
						}
					}
				}
			}
		}
		
		
		
		if (debug) {
			ShowDoubleFloatArrays.showArrays(
					alpha8,
					transform_size,
					transform_size,
					true,
					"alpha8");
			ShowDoubleFloatArrays.showArrays(
					alpha,
					transform_size,
					transform_size,
					true,
					"alpha");
		}
		return alpha;
	}
	
	public static double [][] generateTileAlphas(
			boolean [][] texture_tiles,
			int          width,
			int          transform_size){
		final double [][] alpha_templates =  generateAlphaTemplates(
				transform_size, // final int transform_size,
				false); // boolean debug)

		final int num_slices =  texture_tiles.length;
//		final int img_size = texture_tiles.length/transform_size/transform_size;
		final int tilesX = width /transform_size;
		final int tilesY = texture_tiles[0].length / tilesX;
		final int tiles = tilesX * tilesY;
		final int height = tilesY * transform_size; 
		final TileNeibs tn =     new TileNeibs(tilesX, tilesY);
		final double [][] alpha = new double [num_slices][width * height];
		
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int tile = ai.getAndIncrement(); tile < tiles; tile = ai.getAndIncrement()) if (texture_tiles[fnslice][tile]) {
							int dir_mode = 0;
							for (int dir = 0; dir < 8; dir++) {
								int tile1 = tn.getNeibIndex(tile, dir);
								if ((tile1 >= 0) && texture_tiles[fnslice][tile1]) {
									dir_mode |= (1 << dir);
								}
								int x0 = (tile % tilesX) * transform_size; 
								int y0 = (tile / tilesX) * transform_size;
								int indx0 = x0 + width * y0;
								for (int row = 0; row < transform_size; row++) {
									System.arraycopy(
											alpha_templates[dir_mode],
											row * transform_size,
											alpha[fnslice],
											indx0+ row * width,
											transform_size);
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return alpha;
	}
	public static void generateFgAlphas(
			final boolean [][] texture_en,   // non-transparent texture pixels     
			final double  [][] alphas,       // or null will be updated
			final int          edge,         // now should be 2 or 1
			final double       border_alpha, // alpha along the border 
			final int          width,
			final int          transform_size){
		final int num_slices =  texture_en.length;
		final int img_size =    texture_en[0].length;
		final TileNeibs pn =     new TileNeibs(width, img_size/width);
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nslice = ai.getAndIncrement(); nslice < texture_en.length; nslice = ai.getAndIncrement()) {
						boolean [] texture_edge=pn.getEdgeSelection(
								edge,          // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
								texture_en[nslice], // boolean [] tiles,
								null);     // boolean [] prohibit);
						for (int i = 0; i < texture_en[nslice].length; i++) {
							if (texture_en[nslice][i]){
								alphas[nslice][i] *= texture_edge[i] ? border_alpha : 1.0;
							} else {
								alphas[nslice][i] = 0.0;
							}
						}						
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	public static double [][][] processTexture(
			final CLTParameters  clt_parameters,
			final int            tilesX,
			final double [][]    slice_disparities,
			final int    [][]    slice_border_int,
			final int            border_int_max,
			final double [][][]  sensor_texture,    // per-sensor texture value
			final double [][]    combo_texture_in,  // average texture value
			final TileCluster[]  tileClusters, // to process blue_sky?
			final double         max_disparity_lim, //  = 100.0;  // do not allow stray disparities above this
			final double         min_trim_disparity, //  =  2.0;  // do not try to trim texture outlines with lower disparities
			final TpTask[][][]   tp_tasks_ref,       // reference tasks for each slice to get offsets			
			final String         dbg_prefix_in) {
//		clt_parameters.batch_run		
		final int    max_neib_lev =     clt_parameters.tex_max_neib_lev;      //   2; // 1 - single tiles layer around, 2 - two layers
		/*
		final double var_radius =       1.5; // 3.5;   // for variance filter of the combo disparity
		final double        seed_inter =   150; // 120; // 150;
		final double        seed_same_fz =  6.5; // 13; // seed_inter = 50.0;
		final double        seed_fom =      2.0; // 1.9; // 1.2;
		final double        trim_inter_fz = 5.0; // 13.0;
		final double        trim_fom =      0.5; // 0.8; // 1.3; // 1.8; // 1.2; // 0.8; // 0.4;  // 0.7;
		// scale down fom for pixels near high-variance VAR_SAME
		final double        trim_fom_threshold = 120.0; // count only pixels with VAR_SAME > this value 
		final double        trim_fom_boost = 5;   // boost high-varinace values that exceed threshold  
		final double        trim_fom_blur = 10.0; // divide trim_fom array by blurred version to reduce over sky sharp edge
		// Sure values to set unconditionally transparent and unconditionally opaque FG
		final double        seed_fom_sure =     5.0;
		final double        seed_inter_sure = 150.0; // 13.0;
		final double        trim_fom_sure =    10; //  2.0; temporary disabling it
		final double        min_incr =      100; // temporary disable // 5; // 20.0; // 0.5; // only for sky?
		
		final int           min_neibs_alpha = 1;  // minimal neighbors to keep alpha
		final int           grow_alpha = 0; // 2; // grow alpha selection
		
		final double        alphaOverlapTolerance = 0.0; // exact match only
		final int           reduce_has_bg_grow = 2; // 0 - exactly half tile (between strong and weak)
		final double        occlusion_frac = 0.9;
		final double        occlusion_min_disp = 0.3; // do not calculate occlusions for smaller disparity difference
		final boolean         enhance_map =      false; // debugged, but seems worse - disabling
		final int             min_map_sensors =      3;
		final double          keep_map_frac =        0.6;  // for BG - use this fraction of all sensors in the best direction
		// for fillOcclusionsNaN:
		final int             num_fill_passes =    100;
		final double          max_fill_change =      0.1;
		
		
		// Processing BG and FG trim
		final boolean         en_cut =               true; // enable change FG pixel to transparent from opaque
		final boolean         en_patch =             true; // enable change FG pixel to opaque from transparent
		final double          fg_disp_diff =         1.0;  // do not consider obscuring too close BG (1 pix or more?)
		final int             min_sensors =          4;    // minimal number of sensors visible from the FG pixel
		final double          weight_neib =          3.0; // 2.0; // 1.0;  // weight of same neighbors - add to cost multiplied by num_neib-4
		final double          weight_bg =            0.9; // 0.8; // 1.0; // 15.0/16; // 1.0;  // weight of BG cost relative to the FG one
		final double          best_dir_frac =        0.6;  // for BG - use this fraction of all sensors in the best direction
		final double          cost_min =             1.0;  // minimal absolute value of the total cost to make changes
		final int             max_trim_iterations = 10;
		// debug images
		final boolean         show_debug =              true;
		final boolean         show_update_alpha_slice = false; //true;
		final boolean         show_update_alpha_combo = true;
		final boolean         show_textures_slice =     false; //true;
		final boolean         show_textures_combo =     false; //true;
		final boolean         show_textures_tiles =     false; //true;

		 */
		
		final double        var_radius =      clt_parameters.lre_var_radius;       //    1.5; // 3.5;   // for variance filter of the combo disparity
		final double        seed_inter =      clt_parameters.lre_seed_inter;       //  150; // 120; // 150;
		final double        seed_same_fz =    clt_parameters.lre_seed_same_fz;     //  6.5; // 13; // seed_inter = 50.0;
		final double        seed_fom =        clt_parameters.lre_seed_fom;         //  2.0; // 1.9; // 1.2;
		final double        trim_inter_fz =   clt_parameters.lre_trim_inter_fz;    //   5.0; // 13.0;
		final double        trim_fom =        clt_parameters.lre_trim_fom;         //  0.5; // 0.8; // 1.3; // 1.8; // 1.2; // 0.8; // 0.4;  // 0.7;
		// scale down fom for pixels near high-variance VAR_SAME
		final double        trim_fom_threshold = clt_parameters.lre_trim_fom_threshold; //  120.0; // count only pixels with VAR_SAME > this value 
		final double        trim_fom_boost =  clt_parameters.lre_trim_fom_boost;   //  5;   // boost high-varinace values that exceed threshold  
		final double        trim_fom_blur =   clt_parameters.lre_trim_fom_blur;    //  10.0; // divide trim_fom array by blurred version to reduce over sky sharp edge
		// Sure values to set unconditionally transparent and unconditionally opaque FG
		final double        seed_fom_sure =   clt_parameters.lre_seed_fom_sure;    //     5.0;
		final double        seed_inter_sure = clt_parameters.lre_seed_inter_sure;  //  150.0; // 13.0;
		final double        trim_fom_sure =   clt_parameters.lre_trim_fom_sure;    //  10; //  2.0; temporary disabling it
		final double        min_incr =        clt_parameters.lre_min_incr;         //  100; // temporary disable // 5; // 20.0; // 0.5; // only for sky?
		final int           min_neibs_alpha = clt_parameters.lre_min_neibs_alpha;  //  1;  // minimal neighbors to keep alpha
		final int           grow_alpha =      clt_parameters.lre_grow_alpha;       //  0; // 2; // grow alpha selection
		final double        alphaOverlapTolerance = clt_parameters.lre_alphaOverlapTolerance; //   0.0; // exact match only
		final int           reduce_has_bg_grow =    clt_parameters.lre_reduce_has_bg_grow;    //   2; // 0 - exactly half tile (between strong and weak)
		final double        occlusion_frac =        clt_parameters.lre_occlusion_frac;        //   0.9;
		final double        occlusion_min_disp =    clt_parameters.lre_occlusion_min_disp;    //   0.3; // do not calculate occlusions for smaller disparity difference
		final boolean       enhance_map =           clt_parameters.lre_enhance_map;           //   false; // debugged, but seems worse - disabling
		final int           min_map_sensors =       clt_parameters.lre_min_map_sensors;       //        3;
		final double        keep_map_frac =         clt_parameters.lre_keep_map_frac;         //      0.6;  // for BG - use this fraction of all sensors in the best direction
		// for fillOcclusionsNaN:
		final int           num_fill_passes =       clt_parameters.lre_num_fill_passes;       //    100;
		final double        max_fill_change =       clt_parameters.lre_max_fill_change;       //      0.1;
		// Processing BG and FG trim
		final boolean         en_cut =        clt_parameters.lre_en_cut;        //  true; // enable change FG pixel to transparent from opaque
		final boolean         en_patch =      clt_parameters.lre_en_patch;      //  true; // enable change FG pixel to opaque from transparent
		final double          fg_disp_diff =  clt_parameters.lre_fg_disp_diff;  // 1.0;  // do not consider obscuring too close BG (1 pix or more?)
		final int             min_sensors =   clt_parameters.lre_min_sensors;   // 4;    // minimal number of sensors visible from the FG pixel
		final double          weight_neib =   clt_parameters.lre_weight_neib;   // 3.0; // 2.0; // 1.0;  // weight of same neighbors - add to cost multiplied by num_neib-4
		final double          weight_bg =     clt_parameters.lre_weight_bg;     // 0.9; // 0.8; // 1.0; // 15.0/16; // 1.0;  // weight of BG cost relative to the FG one
		final double          best_dir_frac = clt_parameters.lre_best_dir_frac; // 0.6;  // for BG - use this fraction of all sensors in the best direction
		final double          cost_min =      clt_parameters.lre_cost_min;      // 1.0;  // minimal absolute value of the total cost to make changes
		final int             max_trim_iterations = clt_parameters.lre_max_trim_iterations; //  10;
		// debug images
		final boolean         show_debug =              clt_parameters.lre_show_debug && !clt_parameters.multiseq_run;              //  true;
		final boolean         show_update_alpha_slice = clt_parameters.lre_show_update_alpha_slice; //  false; //true;
		final boolean         show_update_alpha_combo = clt_parameters.lre_show_update_alpha_combo; //  true;
		final boolean         show_textures_slice =     clt_parameters.lre_show_textures_slice;     //  false; //true;
		final boolean         show_textures_combo =     clt_parameters.lre_show_textures_combo;     //  false; //true;
		final boolean         show_textures_tiles =     clt_parameters.lre_show_textures_tiles;     //  false; //true;
		
		final String         dbg_prefix = show_debug ? dbg_prefix_in : null;
		
		final int    num_slices =       sensor_texture.length;
		final int    transform_size =   clt_parameters.transform_size;
		final int    width = tilesX * transform_size;
		final int    img_size = sensor_texture[0][0].length;
		final int    height =   img_size/width;
		
		final int    trim_grow_pix = transform_size * 3;      // 3*transform_size?
		final int    fill_grow =     6*transform_size;
		
		
		final double [][] gcombo_texture = // now always calculate as it has lower noise
				(combo_texture_in != null) ?
						combo_texture_in :
							getComboTexture (sensor_texture);

		
		boolean [][][] tile_booleans = getTileBooleans(
				slice_disparities,                 // final double [][] slice_disparities,
				slice_border_int,                  // final int    [][] slice_border_int, // not extended
				max_disparity_lim,                 // final double      max_disparity_lim,
				min_trim_disparity,                // final double      min_trim_disparity,
				max_neib_lev,                      // final int         max_neib_lev,
				transform_size,                    // final int         transform_size,
				tilesX);                           // final int         tilesX)

		double [][][][] channel_pixel_offsets = getPixelOffsets(
				tp_tasks_ref,  //final TpTask[][][]   tp_tasks_ref, //
				tile_booleans, //final boolean [][][] tile_booleans, // to filter?
				tilesX);       // final int         tilesX)

		if ((dbg_prefix != null) && show_textures_tiles) {
			double [][] dbg_img = new double [tile_booleans[0].length * 5][tile_booleans[0][0].length];
			String[] dbg_titles = new String [tile_booleans[0].length * 5];
			for (int nslice = 0; nslice < tile_booleans[0].length; nslice++) {
				dbg_titles[5*nslice + 0] = "BORDER-"+nslice;
				dbg_titles[5*nslice + 1] = "HAS_BG-"+nslice;
				dbg_titles[5*nslice + 2] = "IS_FG-"+nslice;
				dbg_titles[5*nslice + 3] = "STITCH-STITCHED"+nslice;
				dbg_titles[5*nslice + 4] = "ALL_KEEP"+nslice;
				for (int i = 0; i < dbg_img[0].length; i++) {
					dbg_img[5*nslice + 0][i] = slice_border_int[nslice][i];
					dbg_img[5*nslice + 1][i] = 
							(tile_booleans[TILE_HAS_BG_WEAK][nslice][i]?   1.0 : 0.0) +
							(tile_booleans[TILE_HAS_BG_STRONG][nslice][i]? 2.0 : 0.0) ;
					dbg_img[5*nslice + 2][i] = 
							(tile_booleans[TILE_IS_FG_WEAK][nslice][i]?   1.0 : 0.0) +
							(tile_booleans[TILE_IS_FG_STRONG][nslice][i]? 2.0 : 0.0) ;
					dbg_img[5*nslice + 3][i] = 
							(tile_booleans[TILE_STITCH][nslice][i]? 1.0 : 0.0) +
							(tile_booleans[TILE_STITCHED][nslice][i]? 2.0 : 0.0);
					dbg_img[5*nslice + 4][i] = 
							((slice_border_int[nslice][i] >= 0)? 1.0 : 0.0) +
							(tile_booleans[TILE_KEEP][nslice][i]? 2.0 : 0.0);
				}
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX,
					dbg_img[0].length/tilesX,
					true,
					dbg_prefix+"-tile_booleans",
					dbg_titles);
		}
		
		boolean [][] has_bg_pix = halfStrong(      // select pixels between weak and strong 
				tile_booleans[TILE_HAS_BG_WEAK],   // final boolean [][]  weak_tiles,
				tile_booleans[TILE_HAS_BG_STRONG], // final boolean [][]  strong_tiles,
				transform_size-reduce_has_bg_grow, // final int           grow_tiles,
				transform_size,                    // final int           transform_size,
				tilesX);                           // final int           tilesX)

		boolean [][] is_fg_pix = halfStrong(       // select pixels between weak and strong 
				tile_booleans[TILE_IS_FG_WEAK],    // final boolean [][]  weak_tiles,
				tile_booleans[TILE_IS_FG_STRONG],  // final boolean [][]  strong_tiles,
				transform_size,                    // final int           grow_tiles,
				transform_size,                    // final int           transform_size,
				tilesX);                           // final int           tilesX)
		boolean [][] stitch_pixels = null;
		if (dbg_prefix != null) {
			stitch_pixels = tileToPix(    // expand tile selection to pixel selection 
					tile_booleans[TILE_STITCH],        // final boolean [][]  sel_tiles,
					transform_size,                    // final int           transform_size,
					tilesX);                           // final int           tilesX)
		}
		boolean [][] trim_pixels = getTrimPixels(
				is_fg_pix,   // final boolean [][] is_fg_pix,
				has_bg_pix, // final boolean [][] has_bg_pix)
				tile_booleans[TILE_STITCH],        // final boolean [][] is_stitch_tile,
				transform_size, // final int          transform_size,
				tilesX); // final int          tilesX)
		
		boolean [][] unbound_alpha =  getFgEdge( 
				tile_booleans[TILE_IS_FG_WEAK],    // final boolean [][]  fg_weak_tiles,
				tile_booleans[TILE_IS_FG_STRONG],  // final boolean [][]  fg_strong_tiles,
				tile_booleans[TILE_STITCH],        // final boolean [][]  stitch_tiles,
				is_fg_pix,                         // final boolean [][]  fg_pix,
				transform_size,                    // final int         transform_size,
				tilesX);                           // final int         tilesX)
		// Get vars_same, vars_inter and in debug mode - also 
		final double[][][] vars = getVariances (
				sensor_texture, //				final double [][][] sensor_texture,
				gcombo_texture, // 				final double [][]   combo_texture,
				var_radius, // 				final double        var_radius,
				width); // 				final int           width,
		getTrimSeeds(
				trim_pixels,    // final boolean [][]  trim_pix,  // pixels that may be trimmed
				unbound_alpha,  // final boolean [][]  seed_pix,  // FG edge, just outside of trim_pix. Will be modified
				vars[0],        // final double  [][]  vars_same,
				vars[1],        // 	final double  [][]  vars_inter,
				seed_same_fz,   // final double        seed_same_fz, // add to var_same in denominator
				seed_fom,       // final double        seed_fom,   // minimal value of vars_inter/sqrt(vars_same)
				seed_inter,     // final double        seed_inter, //  =   150;
				width);         // final int           width)		
// copy unbound_alpha here for debug		
		
		final boolean [][] trim_seeds = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				trim_seeds[i] = unbound_alpha[i].clone();
			}
		}
		final double [][][] fom_dbg = (dbg_prefix != null)? new double [4][][] : null;
		double [][] trim_fom_pix = getTrimFom(
				trim_pixels,        // final boolean [][]  trim_pix,  // pixels that may be trimmed
				vars[0],            // final double  [][]  vars_same,
				vars[1],            // final double  [][]  vars_inter,
				trim_inter_fz,      // final double        trim_inter_fz,  // minimal value of vars_same to block propagation
				trim_fom_threshold, //final double        trim_fom_threshold, // = 120.0; // count only pixels with VAR_SAME > this value 
				trim_fom_boost,     //final double        trim_fom_boost, // = 5;   // boost high-varinace values that exceed threshold  
				trim_fom_blur,      // final double        trim_fom_blur,
				width,            // final int           width)
				fom_dbg); // final double [][][] fom_dbg)
		
		getTrimAlpha(
				trim_fom_pix,   // final double  [][]  fom_pix,   // should be NaN outside of trim_pix
				trim_pixels,    // final boolean [][]  trim_pix,  // pixels that may be trimmed
				unbound_alpha,  // final boolean [][]  seed_pix,  // FG edge (just outside of trim_pix) and seeds from vars_inter mismatch
				trim_fom,       // final double        trim_fom, // minimal value of vars_same/vars_inter  to block propagation
				trim_grow_pix,  // final int           trim_grow,      // 3*transform_size?
				width);         // final int           width)

		final boolean [][] first_trimmed_alpha = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				first_trimmed_alpha[i] = unbound_alpha[i].clone();
			}
		}
		// not used:
		final boolean       dual_pass = false; // true;
		expandTrimAlpha(
				trim_pixels,    // final boolean [][]  trim_pix,  // pixels that may be trimmed
				unbound_alpha,  // final boolean [][]  alpha_pix,  //
				vars[0],        // final double  [][]  value,      // will grow only in increasing
				min_incr,       // final double        min_incr,
				dual_pass,      // final boolean       dual_pass,
				trim_grow_pix,  // final int           trim_grow,      // 3*transform_size?
			    width);         // final int           width) {

		final boolean [][] unfiltered_alpha = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				unfiltered_alpha[i] = unbound_alpha[i].clone();
			}
		}
		
		filterAlpha(
				unbound_alpha,   // final boolean [][]  alpha_pix,  // pixels that may be trimmed
				trim_pixels,     // final boolean [][]  trim_pix,   // pixels that may be trimmed
				min_neibs_alpha, // final int           min_neibs,  // minimal neighbors to keep alpha
				grow_alpha,      // final int           grow_alpha, // grow alpha selection
				width);          // final int           width) {

		final boolean [][] filtered_alpha = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				filtered_alpha[i] = unbound_alpha[i].clone();
			}
		}
		
		// remove remaining border FG half-tiles 
		filterWeakFG(
				unbound_alpha,                    // final boolean [][]  alpha_pix,
				trim_pixels,                      // final boolean [][]  trim_pix,			
				has_bg_pix,                       // final boolean [][]  has_bg_pix,			
				tile_booleans[TILE_IS_FG_WEAK],   // final boolean [][]  fg_weak_tiles,
				tile_booleans[TILE_IS_FG_STRONG], // final boolean [][]  fg_strong_tiles,
				transform_size,                   // final int           transform_size,
				tilesX);                          // final int           tilesX);
		
		final boolean [][] weak_fg_alpha = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				weak_fg_alpha[i] = unbound_alpha[i].clone();
			}
		}

		unbound_alpha =  trimAlphaToTiles( // reuse same array
				unbound_alpha,             // final boolean [][]  alpha_pix,
				tile_booleans[TILE_KEEP],  // final boolean [][]  selected_tiles,
				transform_size,            // final int           transform_size,
				tilesX);                   // final int           tilesX)
		
		final boolean [][] before_fix_bg_overlap = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				before_fix_bg_overlap[i] = unbound_alpha[i].clone();
			}
		}

		fix_bg_overlap(
				unbound_alpha,                     // final boolean [][]  alpha_pix,
				tile_booleans[TILE_KEEP],          // final boolean [][]  selected_tiles,
				tile_booleans[TILE_HAS_BG_STRONG], // final boolean [][]  has_bg_strong_tiles,
				slice_border_int,                  // final int      [][] slice_border_int,
				border_int_max,                    // final int           neib_max, // now 2
				transform_size,                    // final int           transform_size,
				tilesX);                           // final int           tilesX)
		
		
		final boolean [][] before_fix_same = (dbg_prefix != null)? new boolean [num_slices][] : null;
		if (dbg_prefix != null) {
			for (int i = 0; i < num_slices; i++) {
				before_fix_same[i] = unbound_alpha[i].clone();
			}
		}
		
		// use minimal alpha if disparity is exactly the same (stitch area)
		fixAlphaSameDisparity( // uses disparities from tileCluster, not slice_disparities ?
				tileClusters,             // final TileCluster [] tileClusters,
				tile_booleans[TILE_KEEP], // final boolean [][]   keep_tiles,
				unbound_alpha,            // final boolean [][]   alpha_pix,
				false,                    // final boolean        use_or, // (maximal alpha), false - and (minimal alpha)
				alphaOverlapTolerance,    // final double         alphaOverlapTolerance, // 0 - require exact match
				width,                    // final int            width,
				transform_size);          // final int            transform_size)

// Processing BG and FG trim
		int [][] occluded_map =                  null;
		int [][] occluded_map_enh =              null;
		double  [][]  dbg_occluded_map =         null;
		double  [][]  occluded_textures =        null;
		double  [][]  occluded_filled_textures = null;
		boolean [][]  sure_transparent =         null;
		boolean [][]  sure_opaque =              null;
		double [][][] debug_costs = ((dbg_prefix != null) && show_update_alpha_slice)? new double [trim_pixels.length][][] : null; 
		int [][]      debug_stats =  (dbg_prefix != null)? new int [trim_pixels.length][] : null;
		boolean [][]  debug_alpha = ((dbg_prefix != null) && show_update_alpha_slice )? new boolean [trim_pixels.length][] : null;
		boolean [][][] dbg_alpha_mods = ((dbg_prefix != null) && show_update_alpha_combo)? new boolean [max_trim_iterations][num_slices][] : null;
		
		if (debug_alpha != null) {
			for (int i = 0; i < unbound_alpha.length; i++) {
				debug_alpha[i] = unbound_alpha[i].clone();
			}
		}

		if (dbg_alpha_mods != null) {
			for (int i = 0; i < unbound_alpha.length; i++) {
				dbg_alpha_mods[0][i] = unbound_alpha[i].clone();
			}
		}
		
		boolean [][] trim_tiles = getTrimTiles(
				trim_pixels,       // boolean [][] trim_pix,
				width,             // final int             width,
				transform_size);   // final int             transform_size);
		int [] occlusionLookUp = null;
		if (enhance_map) {
			final int      dbg_indx = 256;
			occlusionLookUp =  enhanceOccludedMap(
					null,              // final int [][] occluded_map_in,
					keep_map_frac,     // final double   keep_frac,
					min_map_sensors,   // final int      min_sensors,
					sensor_texture[0].length,  // final int      num_sens)
					dbg_indx)[0]; // final int      dbg_indx
		}
		int updated_tiles = 0;
		for (int niter = 0; niter < max_trim_iterations; niter++) {
			occluded_map = getOccludedMap(
					channel_pixel_offsets,            // final double [][][][] channel_pixel_offsets,
					unbound_alpha,                    // final boolean [][]    alpha_pix,
					slice_disparities,                // final double  [][]    slice_disparities,
					tile_booleans[TILE_KEEP],         // final boolean [][]    tile_keep,      // do not check occluded strong foreground
					tile_booleans[TILE_IS_FG_STRONG], // final boolean [][]    tile_fg_strong, // do not check occluded strong foreground
					tile_booleans[TILE_STITCH],       // final boolean [][]    tile_stitch,    // do not process these - there are duplicates
					occlusion_frac,                   // final double          occlusion_frac, // ratio of opaque pixel overlap to consider occlusion
					occlusion_min_disp,               // final double          occlusion_min_disp,
					width,                            // final int             width,
					transform_size);                  // final int             transform_size);
			occluded_map_enh = occluded_map;
			if (occlusionLookUp != null) {
				occluded_map_enh = enhanceOccludedMap(
						occluded_map,     // final int [][] occluded_map,
						occlusionLookUp); //final int []   map)
			}
			dbg_occluded_map = ((dbg_prefix == null) && show_update_alpha_slice)? null:debugOccludedMap(occluded_map_enh); // occluded_map);
			occluded_textures = combineTexturesWithOcclusions(
					sensor_texture,    // final double  [][][] sensor_texture,
					gcombo_texture,    // final double  [][]   combo_texture,
					occluded_map_enh); // occluded_map);     // final int     [][]   occluded_map);
			occluded_filled_textures = fillOcclusionsNaN(
					gcombo_texture,    // final double  [][]   combo_texture,
					occluded_textures, // final double  [][]   combo_occluded_texture,
					fill_grow,         // final int            grow,
					num_fill_passes,   // final int            num_passes,
					max_fill_change,   // final double         max_change,
					width);            // final int            width)
			// Occluded textures should be calculated after updateFgAlpha(), so skip updateFgAlpha() during last iteration 
			if (niter < (max_trim_iterations-1)) {
				sure_transparent = getTrimSeeds(
						trim_pixels,     // final boolean [][]  trim_pix,  // pixels that may be trimmed
						null,            // final boolean [][]  seed_pix_in,  // FG edge, just outside of trim_pix. Will be modified. Or null
						vars[0],         // final double  [][]  vars_same,
						vars[1],         // 	final double  [][]  vars_inter,
						seed_same_fz,    // final double        seed_same_fz, // add to var_same in denominator
						seed_fom_sure,   // final double        seed_fom,   // minimal value of vars_inter/sqrt(vars_same)
						seed_inter_sure, // final double        seed_inter, //  =   150;
						width);          // final int           width)		

				sure_opaque = thresholdAnalog(
						trim_fom_pix,    // final double  [][]  data,
						trim_fom_sure,   // final double        threshold,
						true);           // final boolean       greater)
				if (debug_alpha != null) {
					for (int i = 0; i < unbound_alpha.length; i++) {
						debug_alpha[i] = unbound_alpha[i].clone();
					}
				}

				updated_tiles = updateFgAlpha(
						channel_pixel_offsets,      // final double [][][][] channel_pixel_offsets,
						occluded_filled_textures,   // final double  [][]    textures,
						unbound_alpha,              // final boolean [][]    alpha_pix,
						sensor_texture,             // final double  [][][]  sensor_texture,
						occluded_map,// occluded_map_enh  // final int     [][]    occluded_map,   // bitmap of blocked by FG sensors
						min_sensors,                // final int             min_sensors,    // minimal number of sensors visible from the FG pixel
						slice_disparities,          // final double  [][]    slice_disparities,
						tile_booleans[TILE_KEEP],   // final boolean [][]    tile_keep,      // tiles that have at least one pixel 
						tile_booleans[TILE_STITCH], // final boolean [][]    tile_stitch,      // tiles that have at least one pixel 
						trim_tiles,                 // final boolean [][]    trim_tiles,     // tiles that have at least one pixel 
						trim_pixels,                // final boolean [][]    trim_pix,       // pixels that may be trimmed
						sure_transparent,           // final boolean [][]    transparent,    // definitely transparent
						sure_opaque,                // final boolean [][]    opaque,         // definitely opaque
						en_cut,                     // final boolean         en_cut,         // enable change FG pixel to transparent from opaque
						en_patch,                   // final boolean         en_patch,       // enable change FG pixel to opaque from transparent
						fg_disp_diff,               // final double          min_disp_diff,  // do not consider obscuring too close BG (1 pix or more?)
						// other parameters
						weight_neib,                // final double          weight_neib,    // weight of same neighbors
						weight_bg,                  // final double          weight_bg,      // weight of BG cost relative to the FG one
						best_dir_frac,              // final double          best_dir_frac,  // for BG - use this fraction of all sensors in the best direction
						cost_min,                   // final double          cost_min,       // minimal absolute value of the total cost to make changes
						debug_costs,                // final double [][]     debug_cost,     // if not null, should be double [nslices][] - will return costs/NaN
						debug_stats,                // final int    [][]     debug_stats,    // if not null, should be int [nslices][] - will return number of added/removed per slice
						width,                      // final int             width,
						transform_size);            // final int             transform_size){
				if (dbg_alpha_mods != null) {
					for (int i = 0; i < unbound_alpha.length; i++) {
						dbg_alpha_mods[niter+1][i] = unbound_alpha[i].clone();
					}
				}
			}

			if ((dbg_prefix != null) && show_update_alpha_slice) {
				// TODO:
				// 1. Display alpha mod sequence
				// 2. occluded_map improvements (similar as in updateFgAlpha) - 
				// remove some "unreliable" sensors  
				for (int nslice = 0; nslice < debug_stats.length; nslice++) {
					System.out.println (String.format("#%02d: %5d added, %5d removed (total %5d) opaque FG pixels",
							nslice, debug_stats[nslice][0], debug_stats[nslice][1], debug_stats[nslice][0]+debug_stats[nslice][1]));
				}
				
				String [] dbg_titles0 = {"sure","before","after", "cost", "cost_fg",
						"cost_bg", "combo", "occluded-filed", "occluded"};
				int dbg_len = width * height;
				int sublen = dbg_titles0.length;
				String [] dbg_titles = new String [sublen * num_slices];
				double [][] dbg_img = new double [dbg_titles.length][];
				for (int nslice = 0; nslice< num_slices; nslice++) {
					for (int i = 0; i < dbg_titles0.length; i++) {
						dbg_titles[nslice * sublen + i] = dbg_titles0[i]+"-"+nslice;
					}
					dbg_img[nslice * sublen + 0] = new double [dbg_len];
					dbg_img[nslice * sublen + 1] = new double [dbg_len];
					dbg_img[nslice * sublen + 2] = new double [dbg_len];
					for (int i = 0; i < dbg_len; i++) {
						dbg_img[nslice * sublen + 0][i] =
								(sure_transparent[nslice][i]? 0 : 1) + (sure_opaque[nslice][i]? 2 : 0); 
						dbg_img[nslice * sublen + 1][i] = (debug_alpha[nslice][i]? 3 : 0); 
						dbg_img[nslice * sublen + 2][i] = (unbound_alpha[nslice][i]? 3 : 0); 
					}
					dbg_img[nslice * sublen + 3] = debug_costs[nslice][0];
					dbg_img[nslice * sublen + 4] = debug_costs[nslice][1];
					dbg_img[nslice * sublen + 5] = debug_costs[nslice][2];
					dbg_img[nslice * sublen + 6] = gcombo_texture[nslice];
					dbg_img[nslice * sublen + 7] = occluded_filled_textures[nslice];
					dbg_img[nslice * sublen + 8] = occluded_textures[nslice];
				}
				ShowDoubleFloatArrays.showArrays(
						dbg_img,
						width,
						height,
						true,
						dbg_prefix+"-update_fg-"+niter, // +nslice,
						dbg_titles);
				System.out.println("updateFgAlpha() -> "+updated_tiles);
			}
		}
		
		if (dbg_alpha_mods != null) {
			// TODO:
			// 1. Display alpha mod sequence
			String [] dbg_titles = new String [num_slices * max_trim_iterations];
			int dbg_len = width * height;
			double [][] dbg_img = new double [dbg_titles.length][dbg_len];
			for (int nslice = 0; nslice < num_slices; nslice++) {
				for (int niter = 0; niter < max_trim_iterations; niter++) {
					int slice_indx = nslice * max_trim_iterations + niter; 
					dbg_titles[slice_indx] = "alpha-"+nslice+":"+niter;
					for (int i = 0; i < dbg_len; i++) {
						dbg_img[slice_indx][i] = (dbg_alpha_mods[niter][nslice][i]? 1 : 0); 
					}
				}
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					width,
					height,
					true,
					dbg_prefix+"-update_fg-alpha", // +nslice,
					dbg_titles);
		}
		
		double [][] alphas = new double[num_slices][];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			// replace old alpha with the new binary one
			alphas[nslice] = new double [unbound_alpha[nslice].length];
			for (int i = 0; i < unbound_alpha[nslice].length; i++) {
				alphas[nslice][i] = unbound_alpha[nslice][i] ? 1.0 : 0.0;
			}
		}
		
		if ((dbg_prefix != null) && show_textures_slice) {
			for (int nslice = 0; nslice < num_slices; nslice++) {
				final double [] vars_ratio =               new double [img_size];
				final double [] vars_fom =                 new double [img_size]; // inter/sqrt(same)
				final double [] half_pix =                 new double [img_size];
				final double [] stitch_trim_pix =          new double [img_size];
				final double [] trim_seed_pix =            new double [img_size];
				final double [] seed_trim_grow_pix =       new double [img_size];
				final double [] unfilt_filt_pix =          new double [img_size];
				final double [] weak_fg_pix =              new double [img_size];
				final double [] trim_tiles_pix =           new double [img_size];
				final double [] fix_bg_pix =               new double [img_size];
				final double [] fix_same_pix =             new double [img_size];
				final double [] trim_alpha_pix =           new double [img_size];
				for (int i = 0; i <img_size; i++) {
					vars_fom[i] =       vars[1][nslice][i]/(vars[0][nslice][i]+seed_same_fz);
					if (Double.isNaN(vars_ratio[i])) vars_ratio[i] = 0;
					if (Double.isNaN(vars_fom[i]))   vars_fom[i] = 0;
					half_pix[i] = 
							(has_bg_pix [nslice][i]? 1.0:0.0) +
							(is_fg_pix  [nslice][i]? 2.0:0.0);
					stitch_trim_pix[i] = 
							(stitch_pixels [nslice][i]? 1.0:0.0) +
							(trim_pixels    [nslice][i]? 2.0:0.0);
					trim_seed_pix[i] = 
							(trim_pixels [nslice][i]? 1.0:0.0) +
							(trim_seeds  [nslice][i]? 2.0:0.0);
					
					seed_trim_grow_pix[i] = 
							(trim_seeds         [nslice][i]? 1.0:0.0) +
							(first_trimmed_alpha[nslice][i]? 2.0:0.0) + // after first trimming
							(unfiltered_alpha   [nslice][i]? 4.0:0.0);  // grown by var_same increase 
					
					unfilt_filt_pix[i] = 
							(unfiltered_alpha  [nslice][i]? 1.0:0.0) +
							(filtered_alpha    [nslice][i]?  2.0:0.0);
					weak_fg_pix[i] = 
							(filtered_alpha    [nslice][i]? 1.0:0.0) +
							(weak_fg_alpha     [nslice][i]?  2.0:0.0);
					trim_tiles_pix[i] = 
							(weak_fg_alpha     [nslice][i]? 1.0:0.0) +
							(before_fix_same   [nslice][i]?  2.0:0.0);
					fix_bg_pix[i] = 
							(before_fix_same [nslice][i]? 1.0:0.0) +
							((before_fix_same [nslice][i] ^ before_fix_bg_overlap[nslice][i])? 2.0:0.0);
					fix_same_pix[i] = 
							(unbound_alpha  [nslice][i]? 1.0:0.0) +
							((unbound_alpha  [nslice][i] ^ before_fix_same  [nslice][i])?  2.0:0.0);
					
					trim_alpha_pix[i] = 
							(trim_pixels  [nslice][i]? 1.0:0.0) +
							(unbound_alpha[nslice][i]?  2.0:0.0);
				}
				double [][] dbg_img = {
						vars[0][nslice],
						vars[1][nslice],
						vars[2][nslice],
						vars[3][nslice],
						vars[4][nslice],
						trim_fom_pix[nslice], // normalized by blurred
						
						fom_dbg[0][nslice],
						fom_dbg[1][nslice],
						fom_dbg[2][nslice],
						fom_dbg[3][nslice],
						vars_fom,
						half_pix,
						stitch_trim_pix,
						trim_seed_pix,
						seed_trim_grow_pix,
						unfilt_filt_pix,
						weak_fg_pix,
						trim_tiles_pix,
						fix_bg_pix,
						fix_same_pix,
						trim_alpha_pix,
						dbg_occluded_map[nslice],
						gcombo_texture[nslice],
						occluded_filled_textures[nslice], // put before occluded_textures to compare with gcombo_texture
						occluded_textures[nslice],
						sensor_texture[nslice][ 0],
						sensor_texture[nslice][ 1],
						sensor_texture[nslice][ 2],
						sensor_texture[nslice][ 3],
						sensor_texture[nslice][ 4],
						sensor_texture[nslice][ 5],
						sensor_texture[nslice][ 6],
						sensor_texture[nslice][ 7],
						sensor_texture[nslice][ 8],
						sensor_texture[nslice][ 9],
						sensor_texture[nslice][10],
						sensor_texture[nslice][11],
						sensor_texture[nslice][12],
						sensor_texture[nslice][13],
						sensor_texture[nslice][14],
						sensor_texture[nslice][15]
				};
				String [] dbg_titles = {
						"VAR_SAME",
						"VAR_INTER",
						"GRAD_X",
						"GRAD_Y",
						"GRAD_ABS",
						"TRIM_FOM_NORM", // same/(inter+trim_inter_fz) normalized by blurred version
						"TRIM_FOM_INI", // initial fom
						"TRIM_FOM_FIN", // final fom
						"TRIM_FOM_THRESH", // var_same_thresholded
						"TRIM_FOM_THRESH_BLUR", // var_same_thresholded_blured
						"SEED_FOM", // inter/(same+seed_same_fz)
						"HALF_BG_FG",
						"STITCH_TRIM",
						"TRIM_SEED",
						"SEED_TRIMMED_MORE",
						"UNFILT_FILT",
						"WEAK_FG",
						"TRIM_TILES",
						"FIX_HAS_BG",
						"FIX_SAME",
						"TRIM_ALPHA",
						"OCCLUSIONS_MAP",
						"COMBO_TEXTURE",
						"OCCLUDED_FILLED_TEXTURES",
						"OCCLUDED_TEXTURES",
						"T00",
						"T01",
						"T02",
						"T03",
						"T04",
						"T05",
						"T06",
						"T07",
						"T08",
						"T09",
						"T10",
						"T11",
						"T12",
						"T13",
						"T14",
						"T15"
				};
				ShowDoubleFloatArrays.showArrays(
						dbg_img,
						width,
						height,
						true,
						dbg_prefix+"-textures-"+nslice,
						dbg_titles);
				assert true;
			}
		}
		if ((dbg_prefix != null) && show_textures_combo) {
			ShowDoubleFloatArrays.showArrays(
					occluded_filled_textures, // out_textures,
					width,
					height,
					true,
					dbg_prefix+"-out_textures");
			ShowDoubleFloatArrays.showArrays(
					alphas,
					width,
					height,
					true,
					dbg_prefix+"-alphas");
			
			double [][] masked_textures = 	combineTextureAlpha(
					0.5,                      // final double         alpha_threshold,
					occluded_filled_textures, // out_textures, // final double  [][]   textures,
					alphas);                  // final double  [][]   alphas
			ShowDoubleFloatArrays.showArrays(
					masked_textures,
					width,
					height,
					true,
					dbg_prefix+"-masked_textures");
		}

		double [][][] textures_alphas = new double [num_slices][][];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			alphas[nslice] = new double [unbound_alpha[nslice].length];
			for (int i = 0; i < unbound_alpha[nslice].length; i++) {
				alphas[nslice][i] = unbound_alpha[nslice][i] ? 1.0 : 0.0;
			}
			textures_alphas[nslice] = new double [][] {occluded_filled_textures[nslice], alphas[nslice]};
		}
		// set slice_disparities to NaN for unselected tiles - it will update tileClusters
		setMeshTileSelection(
				slice_disparities,         // final double  [][] slice_disparities,
				tile_booleans[TILE_KEEP]); //final boolean [][] keep_tiles)
		String dbg_prefix1 = (dbg_prefix==null)? null: (dbg_prefix+"-masked");
		if ((dbg_prefix1!=null) && show_textures_tiles) {
		showDebugDisparities( // nop if dbg_prefix== null
				slice_disparities, // final double [][] slice_disparities,
				tilesX,            // final int   tilesX,
				dbg_prefix1);       // String      prefix);
		}
		return textures_alphas; // What about colors? 
	}
	
	/**
	 * Generate combined full-size image textures, currently they are split, maybe will be used as is
	 * to reduce total number of images.
	 *  
	 * @param clt_parameters      processing parameters.
	 * @param colorProcParameters (older) parameters related to color image representation.
	 * @param rgbParameters       (older) another color-processing parameters.
	 * @param parameter_scene     scene (QuadCLT instance) to use for rendering parameters in multi-series sequences
	 *                            if null - use reference scene instead.
	 * @param ref_index           index of the reference scene (currently the last one in the scenes sequence).
	 * @param scenes              array of scenes (QuadCLT instances).
	 * @param scenes_sel          null ar binary array to select which scenes to process (now all true).
	 * @param tileClusters        combined tile clusters (each containing info about containing sub-clusters).
	 * @param renormalize         false - use normalizations from previous scenes to keep consistent colors.
	 *                            true - re-normalize rendered textures.
	 * @param max_disparity_lim   do not allow stray disparities above this (now 100.0).
	 * @param min_trim_disparity  do not try to trim texture outlines with lower disparities (now 2.0).
	 * @param debugLevel          debug level - controls generation of images.
	 * @return                    array of ImagePlus instances corresponding to tileClusters array
	 */
//	public static ImagePlus[] getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
	public static double[][][] getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
			final CLTParameters  clt_parameters,
			ColorProcParameters  colorProcParameters,
//			EyesisCorrectionParameters.RGBParameters rgbParameters,
//			final boolean        no_alpha,
			QuadCLT              parameter_scene,
			final int            ref_index,
			final QuadCLT []     scenes,
			final boolean []     scenes_sel,
			final TileCluster [] tileClusters,
			final boolean        renormalize,
			final double         max_disparity_lim,
			final double         min_trim_disparity,
			final int            debugLevel)
	{
		// TODO: ***** scenes with high motion blur also have high ERS to be corrected ! *****
		final QuadCLT ref_scene = scenes[ref_index];
		if (parameter_scene == null) {
			parameter_scene = ref_scene;
		}
		final int earliestScene = ref_scene.getEarliestScene(scenes);
		final ErsCorrection ers_reference = ref_scene.getErsCorrection();
		final int tilesX =                  ref_scene.getTileProcessor().getTilesX();
		final int tilesY =                  ref_scene.getTileProcessor().getTilesY();
		final int tiles =                   tilesX * tilesY;
		final int transform_size=           ref_scene.getTileProcessor().getTileSize();
		final int tile_len =                transform_size * transform_size;
		final boolean filter_bg =           true; // make a clt parameter?
		final boolean mb_en =               clt_parameters.imp.mb_en;
		final double  mb_tau =              clt_parameters.imp.mb_tau;      // 0.008;// time constant, sec
		final double  mb_max_gain =         clt_parameters.imp.mb_max_gain; // 5.0;  // motion blur maximal gain (if more - move second point more than a pixel

		final double  max_distortion =      clt_parameters.tex_distort;   // 0.5;  // Maximal texture distortion to accumulate multiple scenes (0 - any)
		final double  tex_mb =              clt_parameters.tex_mb;        // 1.0;  // Reduce texture weight if motion blur exceeds this (as square of MB length)
//		final boolean sharp_alpha =         clt_parameters.sharp_alpha;
		final boolean is_lwir =             ref_scene.isLwir();
		final boolean tex_um =              clt_parameters.tex_um;        // imp.um_mono; // TODO: add own parameter
		final double  tex_um_sigma =        clt_parameters.tex_um_sigma;  // imp.um_sigma;
		final double  tex_um_weight =       clt_parameters.tex_um_weight; // imp.um_weight;
		// TODO: - make texture variants, tex_um_fixed/tex_um_range apply only to unsharp mask, regardless of colors
		final boolean lwir_autorange =      is_lwir && clt_parameters.tex_lwir_autorange; // colorProcParameters.lwir_autorange;
		final boolean tex_um_fixed =        clt_parameters.tex_um_fixed;  // imp.mono_fixed; //  true; // normalize to fixed range when converting to 8 bits 
		final double  tex_um_range =        clt_parameters.tex_um_range;  // imp.mono_range; // 500.0;  // monochrome full-scale range (+/- half)
		final boolean tex_hist_norm =       clt_parameters.tex_hist_norm; //  true;  
		final double  tex_hist_amount =     clt_parameters.tex_hist_amount; // clt_parameters. 0.7;  
		final int     tex_hist_bins =       clt_parameters.tex_hist_bins;   //  1024 ;   
		final int     tex_hist_segments =   clt_parameters.tex_hist_segments; // 32 ;   
///		final boolean tex_color =           clt_parameters.tex_color;     //  true;  
///		final int     tex_palette =         clt_parameters.tex_palette;     // 2 ;   
//		final boolean extend_sky =          true;
		final int     shrink_sky_tiles =    4; // 2; sum of 2 +bg extend 
		final boolean grow_sky =            true;
		final boolean alphaOverlapFix =     true; // if multiple tiles have the same (+/-?) disparity, make alpha max of them
		final double  alphaOverlapTolerance = 0.0; // compare same disparity with tolerance (relative to disparity? make absolute meters?)
		ImageDtt image_dtt;
		image_dtt = new ImageDtt(
				ref_scene.getNumSensors(), // numSens,
				transform_size,
				clt_parameters.img_dtt,
				ref_scene.isAux(),
				ref_scene.isMonochrome(),
				ref_scene.isLwir(),
				clt_parameters.getScaleStrength(ref_scene.isAux()),
				ref_scene.getGPU());
		if (ref_scene.getGPU() != null) {
			ref_scene.getGPU().setGpu_debug_level(debugLevel);
		}
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		final int num_slices =               tileClusters.length;
		double [][][]     inter_weights =    new double [num_slices][tilesY][tilesX]; // per-tile texture weights for inter-scene accumulation;
		// weighted sum
		double [][][][][] inter_textures_wd= new double [num_slices][tilesY][tilesX][][]; // [channel][64] - overlapping textures
		// weighted sum of squares
		double [][][] ref_pXpYDs =           new double [num_slices][][]; // individual for each slice
		int    [][] cluster_indices =        (max_distortion > 0.0) ? (new int [num_slices][]): null;
		boolean [][] borders =               new boolean [num_slices][];
		for (int nslice = 0; nslice < num_slices; nslice++) { // prepare and measure textures for each combo textures
			ref_pXpYDs[nslice] = OpticalFlow.transformToScenePxPyD( // now should work with offset ref_scene
					null, // fov_tiles,                  // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
					tileClusters[nslice].getDisparity(), // final double []   disparity_ref, // invalid tiles - NaN in disparity
					OpticalFlow.ZERO3,                   // final double []   scene_xyz, // camera center in world coordinates
					OpticalFlow.ZERO3,                   // final double []   scene_atr, // camera orientation relative to world frame
					scenes[ref_index],                   // final QuadCLT     scene_QuadClt,
					scenes[ref_index],                   // final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
					THREADS_MAX);                        // int               threadsMax)
			borders[nslice] = tileClusters[nslice].getBorder();
			if (max_distortion > 0.0) {
				cluster_indices[nslice] = tileClusters[nslice].getClusterIndex();
			}
		}		
		
		final int num_sensors = parameter_scene.getNumSensors();
		final int num_colors =  parameter_scene.isMonochrome()?1:3;
		
		final double [][][] sensor_textures = new double [num_slices][num_sensors][];
		final double [][] combo_textures = new double [num_slices][];
		final TpTask[][][] tp_tasks_ref = new TpTask [num_slices][][];
		for (int nscene = earliestScene; nscene < scenes.length; nscene++) if ((scenes_sel == null) || scenes_sel[nscene]){
			String ts = scenes[nscene].getImageName();
			double []   scene_xyz = OpticalFlow.ZERO3;
			double []   scene_atr = OpticalFlow.ZERO3;
			if (nscene != ref_index) {
				scene_xyz = ers_reference.getSceneXYZ(ts);
				scene_atr = ers_reference.getSceneATR(ts);
				if ((scene_xyz == null) || (scene_atr == null)){
					continue; // scene is not matched
				}
				double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
				scenes[nscene].getErsCorrection().setErsDt(
						scene_ers_xyz_dt, // double []    ers_xyz_dt,
						scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
			}
			double [][] dxyzatr_dt = null;
			// should get velocities from HashMap at reference scene from timestamp , not re-calculate.
			if (mb_en) { // all scenes have the same name/path
				dxyzatr_dt = new double[][] { // for all, including ref
					scenes[nscene].getErsCorrection().getErsXYZ_dt(),
					scenes[nscene].getErsCorrection().getErsATR_dt()};				
			}
			scenes[nscene].saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU)
			//parameter_scene
			for (int nslice = 0; nslice < num_slices; nslice++) { // prepare and measure textures for each combo textures
				final double [] disparity_ref = tileClusters[nslice].getDisparity();  // disparity in the reference view tiles (Double.NaN - invalid)
				// Motion blur vectors are individual per-slice
				// Calculate motion blur vectors - may be used to modify weights of averaged textures
				final double [][] motion_blur = (mb_en && (dxyzatr_dt != null))? OpticalFlow.getMotionBlur(
						scenes[ref_index],   // QuadCLT        ref_scene,
						scenes[nscene],      // QuadCLT        scene,         // can be the same as ref_scene
						ref_pXpYDs[nslice],  // double [][]    ref_pXpYD,     // here it is scene, not reference!
						scene_xyz,           // double []      camera_xyz,
						scene_atr,           // double []      camera_atr,
						dxyzatr_dt[0],       // double []      camera_xyz_dt,
						dxyzatr_dt[1],       // double []      camera_atr_dt,
						0,                   // int            shrink_gaps,  // will gaps, but not more that grow by this
						debugLevel) : null;        // int            debug_level)
				if (debugLevel > 0) {
					System.out.println("nscene="+nscene+", nslice="+nslice+" will run texturesGPUFromDSI() that needs debug >2");
					System.out.print("");
				}
				if ((debugLevel > -1) && (nscene == ref_index)) { // change to "-2" to activate
					System.out.println("Processing reference scene");
					System.out.print("");
				}
				final TpTask[][][] tp_tasks_ret = ((nscene == ref_index) && (tp_tasks_ref != null))?
								new TpTask[1][][] : null;
				double [][][][] slice_texture88 = QuadCLT.texturesNoOverlapGPUFromDSI(
						clt_parameters,          // CLTParameters     clt_parameters,
						disparity_ref,           // double []         disparity_ref,
						// motion blur compensation 
						mb_tau,                  // double            mb_tau,      // 0.008; // time constant, sec
						mb_max_gain,             // double            mb_max_gain, // 5.0;   // motion blur maximal gain (if more - move second point more than a pixel
						motion_blur,             // double [][]       mb_vectors,  // now [2][ntiles];
						scene_xyz,               // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,               // final double []   scene_atr, // camera orientation relative to world frame
						scenes[nscene],          // final QuadCLT     scene,
						scenes[ref_index],       // final QuadCLT     ref_scene, // now - may be null - for testing if scene is rotated ref
						filter_bg && (nscene != ref_index),    // final boolean filter_bg, // remove bg tiles (possibly occluded)
						max_distortion,          // final double      max_distortion, // maximal neighbor tiles offset as a fraction of tile size (8)
						cluster_indices[nslice], // final int []      cluster_index,  //
						borders[nslice],         // final boolean []  border, // border tiles
						// without the following uniform sky develops horizontal lines caused by image edge tiles on scenes where
						// they are over clear part of the reference scene window
						10,                      // final int         discard_frame_edges, // do not use tiles that have pixels closer to the frame margins 
						1,                       // final int         keep_frame_tiles, // do not discard pixels for border tiles in reference frame 
						true,                    // keep_channels,           // final boolean     keep_channels,
						tp_tasks_ret,            // final TpTask[][][] tp_tasks_ret, // if not null, should be [1] - will return tp_tasks_ret[0] = tp_tasks
						debugLevel);             // final int         debugLevel);
				if (tp_tasks_ret != null) {
					tp_tasks_ref[nslice] = 	tp_tasks_ret[0];
				}
				if (slice_texture88 != null) { // will just accumulate
					// Use MB vectors for texture weights				
					final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
					final AtomicInteger ai = new AtomicInteger(0);
					final int fnslice = nslice;
					final double mb_tau2 = mb_tau * mb_tau / tex_mb / tex_mb;
					for (int ithread = 0; ithread < threads.length; ithread++) {
						threads[ithread] = new Thread() {
							public void run() {
								for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
									int tileX = nTile % tilesX;
									int tileY = nTile / tilesX;
									if (slice_texture88[tileY][tileX] != null) {
										double w = 1.0;
										if (tex_mb > 0.0) {
											double mb_l2 = mb_tau2 * ( motion_blur[0][nTile]*motion_blur[0][nTile] + // motion_blur == null;
													motion_blur[1][nTile]*motion_blur[1][nTile]);
											if (mb_l2 > 1.0) {
												w /= mb_l2; // 1/(squared mb)
											}
										}
										if (w > 0) {
											inter_weights[fnslice][tileY][tileX] +=w;
											if (inter_textures_wd[fnslice][tileY][tileX] == null) { // create if it did not exist
												inter_textures_wd[fnslice][tileY][tileX] =  new double [slice_texture88[tileY][tileX].length + num_colors][slice_texture88[tileY][tileX][0].length];
											}
											for (int nchn = 0; nchn < slice_texture88[tileY][tileX].length; nchn++) {
												for (int i = 0; i < slice_texture88[tileY][tileX][nchn].length; i++) {
													double d = slice_texture88[tileY][tileX][nchn][i];
													inter_textures_wd [fnslice][tileY][tileX][nchn][i] += w * d; 	
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
				if (debugLevel > -1) { // -2
					if (nscene == ref_index) {
						System.out.println("Textures from the reference scene, nslice = " + nslice +((slice_texture88 == null)? " - EMPTY":""));
					} else {
						System.out.println("Textures from scene "+nscene+", slice="+nslice +((slice_texture88 == null)? " - EMPTY":""));
					}
				}
			} // for (int nslice = 0; nslice < num_slices; nslice++) {
		} // for (int nscene = 0; nscene < num_scenes; nscene++) {
		
		// Divide accumulated data by weights
		final double [][] dbg_weights = (debugLevel > 0 )?(new double [num_slices][tiles]) : null;
		final int width = tilesX * transform_size;
		final int height = tilesY * transform_size;
		final int y_color = num_colors-1;
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int nslice = 0; nslice < num_slices; nslice++) {
			for (int nsens = 0; nsens < num_sensors; nsens++) {
				sensor_textures[nslice][nsens] = new double [width*height];
				Arrays.fill(sensor_textures[nslice][nsens], Double.NaN);
			}
			combo_textures[nslice] = new double [width*height];
			Arrays.fill(combo_textures[nslice], Double.NaN);
			final int fnslice = nslice;
			if (dbg_weights != null) {
				Arrays.fill(dbg_weights[nslice], Double.NaN);
			}
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							int tileX = nTile % tilesX;
							int tileY = nTile / tilesX;
							if (inter_weights[fnslice][tileY][tileX] > 0.0) {
								if (dbg_weights != null) {
									dbg_weights[fnslice][nTile] = inter_weights[fnslice][tileY][tileX];		
								}
								double w = 1.0/ inter_weights[fnslice][tileY][tileX];
								for (int nchn = 0; nchn < inter_textures_wd[fnslice][tileY][tileX].length - num_colors; nchn++) {
									for (int i = 0; i < inter_textures_wd[fnslice][tileY][tileX][nchn].length; i++) {
										double d = inter_textures_wd[fnslice][tileY][tileX][nchn][i] * w; // average
										inter_textures_wd[fnslice][tileY][tileX][nchn][i] = d;
									}
								}
								for (int ncol = 0; ncol < num_colors; ncol++) {
									int navg = inter_textures_wd[fnslice][tileY][tileX].length - num_colors + ncol;
									for (int i = 0; i < tile_len; i++) {
										inter_textures_wd[fnslice][tileY][tileX][navg][i] = 0;
										for (int nsens = 0; nsens < num_sensors; nsens++) {
											inter_textures_wd[fnslice][tileY][tileX][navg][i] += 
													inter_textures_wd[fnslice][tileY][tileX][1 + (nsens + 1) * num_colors][i]/num_sensors;
										}
									}
								}
								// important that texture is not overlapped here!
								for (int row = 0; row < transform_size; row++) {
									for (int nsens = 0; nsens < num_sensors; nsens++) {
										System.arraycopy(
												inter_textures_wd[fnslice][tileY][tileX][num_colors* (nsens + 1) + 1 + y_color],
												row*transform_size,
												sensor_textures[fnslice][nsens] ,// dbg_textures[n],
												(tileY * transform_size + row) * width + (tileX * transform_size),
												transform_size);
									}
									int navg = y_color; // inter_textures_wd[fnslice][tileY][tileX].length - num_colors + y_color;
									System.arraycopy(
											inter_textures_wd[fnslice][tileY][tileX][navg],
											row*transform_size,
											combo_textures[fnslice] ,// dbg_textures[n],
											(tileY * transform_size + row) * width + (tileX * transform_size),
											transform_size);
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
		} // for (int nslice = 0; nslice < num_slices; nslice++) {	
		
		final double [][] slice_disparities = new double [num_slices][];
		final int [][]    slice_border_int = new int [num_slices][];
		for (int nslice = 0; nslice < num_slices; nslice++) {
			slice_disparities[nslice] = tileClusters[nslice].getDisparity();  // disparity in the reference view tiles (Double.NaN - invalid)
			slice_border_int[nslice] =  tileClusters[nslice].getBorderInt(); 
		}
		int border_int_max = tileClusters[0].getBorderIntMax();
		double [][][] faded_textures = processTexture( //[slice]{texture, alpha}
				clt_parameters,            // final CLTParameters  clt_parameters,
				tilesX,                    // final int            tilesX,
				slice_disparities,         // final double [][]    slice_disparities,
				slice_border_int,          // final int    [][]    slice_border_int,
				border_int_max,            // final int            border_int_max,
				sensor_textures,           // final double [][]    sensor_texture,    // per-sensor texture value
				null, // combo_textures,   // null, // final double []      combo_texture_in,  // average texture value
				tileClusters,              // final TileCluster[]  tileClusters, // to process blue_sky?
				max_disparity_lim,         // final double max_disparity_lim,  // do not allow stray disparities above this
				min_trim_disparity,        // final double min_trim_disparity, // do not try to trim texture outlines with lower disparities
				tp_tasks_ref,              // final TpTask[][][]   tp_tasks_ref,       // reference tasks for each slice to get offsets			
				ref_scene.getImageName()); // null); // ref_scene.getImageName()); // final String         dbg_prefix);
		if (debugLevel > -1) {
			double [][] dbg_textures = new double [faded_textures.length * faded_textures[0].length][faded_textures[0][0].length];
			String [] dbg_titles = new String[dbg_textures.length];
			String [] dbg_subtitles = new String [faded_textures[0].length];
			for (int i = 0; i < dbg_subtitles.length; i++) {
				dbg_subtitles[i] = (i <  (dbg_subtitles.length -1)) ? ("Y"+i):"alpha";
			}
			for (int i = 0; i < dbg_textures.length; i++) {
				dbg_textures[i] = faded_textures[i / faded_textures[0].length][i % faded_textures[0].length];
				dbg_titles[i] = dbg_subtitles[i % dbg_subtitles.length] + "-" + (i / dbg_subtitles.length);
				
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_textures,
					tilesX * transform_size,
					tilesY * transform_size,
					true,
					ref_scene.getImageName()+"-combined_textures-prenorm-pre_UM",
					dbg_titles);
			
		}
		// Grow sky
		if (grow_sky) {
			extendBlueSKy(
					tileClusters,     // final TileCluster [] tileClusters,
					faded_textures,   // final double [][][]  faded_textures,
					shrink_sky_tiles, // final int            shrink_sky_tiles,
					width,            // final int            width,
					transform_size);  // final int            transform_size);
		}
		// fix alpha
		if (alphaOverlapFix) {
			fixAlphaSameDisparity(
					tileClusters,          // final TileCluster [] tileClusters,
					faded_textures,        // final double [][][]  faded_textures,
					alphaOverlapTolerance, // final int            alphaOverlapTolerance, // 0 - require exact match
					width,                 // final int            width,
					transform_size);       // final int            transform_size)			
		}
		
		// Is it needed here? Or move to processTexture()? - Slow
		for (int nslice = 0; nslice < faded_textures.length; nslice++) {
			faded_textures[nslice][0] = TileProcessor.fillNaNs(
					faded_textures[nslice][0], // final double [] data,
					null,                     // final boolean [] prohibit,
					width,        // int       width,
					16,           // final int grow,
					0.7,          // double    diagonal_weight, // relative to ortho
					100,          // int       num_passes,
					0.01,                     // final double     max_rchange, //  = 0.01
					THREADS_MAX); // final int threadsMax)      // maximal number of threads to launch 
		}
		
		if (debugLevel > -1) {
			double [][] dbg_textures = new double [faded_textures.length * faded_textures[0].length][faded_textures[0][0].length];
			String [] dbg_titles = new String[dbg_textures.length];
			String [] dbg_subtitles = new String [faded_textures[0].length];
			for (int i = 0; i < dbg_subtitles.length; i++) {
				dbg_subtitles[i] = (i <  (dbg_subtitles.length -1)) ? ("Y"+i):"alpha";
			}
			
			for (int i = 0; i < dbg_textures.length; i++) {
				dbg_textures[i] = faded_textures[i / faded_textures[0].length][i % faded_textures[0].length];
				dbg_titles[i] = dbg_subtitles[i % dbg_subtitles.length] + "-" + (i / dbg_subtitles.length);
				
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_textures,
					tilesX * transform_size,
					tilesY * transform_size,
					true,
					ref_scene.getImageName()+"-combined_textures-filled-NaNs",
					dbg_titles);
		}		
		
		// Optionally apply UM (before auto/manual range)
		if (tex_um) {
			QuadCLTCPU.umTextures(
					faded_textures, // final double [][][] textures, //  [nslices][nchn][i]
					tilesX * transform_size, // final int    width,
					tex_um_sigma, // final double um_sigma,
					tex_um_weight); // final double um_weight)
		}
		
		if (debugLevel > -1) {
			double [][] dbg_textures = new double [faded_textures.length * faded_textures[0].length][faded_textures[0][0].length];
			String [] dbg_titles = new String[dbg_textures.length];
			String [] dbg_subtitles = new String [faded_textures[0].length];
			for (int i = 0; i < dbg_subtitles.length; i++) {
				dbg_subtitles[i] = (i <  (dbg_subtitles.length -1)) ? ("Y"+i):"alpha";
			}
			
			for (int i = 0; i < dbg_textures.length; i++) {
				dbg_textures[i] = faded_textures[i / faded_textures[0].length][i % faded_textures[0].length];
				dbg_titles[i] = dbg_subtitles[i % dbg_subtitles.length] + "-" + (i / dbg_subtitles.length);
				
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_textures,
					tilesX * transform_size,
					tilesY * transform_size,
					true,
					ref_scene.getImageName()+"-combined_textures-prenorm",
					dbg_titles);
			if (dbg_weights != null) {
				ShowDoubleFloatArrays.showArrays(
						dbg_weights,
						tilesX,
						tilesY,
						true,
						ref_scene.getImageName()+"-texture_weights-prenorm");
			}
		}
		
		
		//renormalize
		// normalize all slices together if LWIR
		// FIXME: Should it be here? Will setColdHot() change photometric calibration ? Or should it be disabled?
		double [] norm_table = null; // first try, then make save to properties with cold/hot
		if (renormalize) {
			if (lwir_autorange) {
				double rel_low;
				double rel_high;
				boolean force_min_max = true;
				if (!tex_um && !force_min_max) { // for UM will use min/max
					rel_low =  colorProcParameters.lwir_low;
					rel_high = colorProcParameters.lwir_high;
					if (!Double.isNaN(parameter_scene.getLwirOffset())) { // ref_scene or parameter_scene? Or both?
						rel_low -=  parameter_scene.getLwirOffset();
						rel_high -= parameter_scene.getLwirOffset();
					}
				} else { // for UM need to calculate min and max (probably OK for non-UM too !)
					double [] minmax = QuadCLTCPU.getMinMaxTextures(
							faded_textures ); //double [][][] textures //  [slices][nchn][i]
					rel_low =  minmax[0]; // absolute min
					rel_high = minmax[1]; // absolute max
				}
				double [] cold_hot =  QuadCLTCPU.autorangeTextures(
						faded_textures,                    // double [][][] textures, //  [nslices][nchn][i]
						rel_low,                           // double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
						rel_high,                          // double hard_hot,   // matches data, DC (this.lwir_offset)  subtracted
						colorProcParameters.lwir_too_cold, // double too_cold, // pixels per image
						colorProcParameters.lwir_too_hot,  // double too_hot,  // pixels per image
						tex_hist_bins); // int num_bins)
				if ((cold_hot != null) && !tex_um && !force_min_max) {
					if (!Double.isNaN(parameter_scene.getLwirOffset())) {
						cold_hot[0] += parameter_scene.getLwirOffset();
						cold_hot[1] += parameter_scene.getLwirOffset();
					}
				}
				parameter_scene.setColdHot(cold_hot); // will be used for shifted images and for texture tiles
			} else if (tex_um && tex_um_fixed) { // apply fixed range, but for UM only (what about RGB?)
				parameter_scene.setColdHot(-0.5*tex_um_range, 0.5*tex_um_range);
			}
			
			if (tex_hist_norm) { // will normalize (0..1) keeping cold_hot to apply during rendering
				// last norm_table element is <=1.0, first >=0;
				 norm_table = QuadCLTCPU.getHistogramNormalization(
						  faded_textures,               // double [][][] textures, //  [nslices][nchn][i]
						  parameter_scene.getColdHot(), // double [] minmax,
						  tex_hist_bins,                    // int       num_bins,
						  tex_hist_segments,               //int       num_nodes
						  tex_hist_amount);  //double    hist_normalize_amount // 1.0 - full
			}
		}
		if (tex_hist_norm && (norm_table != null)) {
			// apply histogram normalization
			double [] cold_hot = parameter_scene.getColdHot(); // used in linearStackToColor
			double [] inverted_table = QuadCLTCPU.invertHistogramNormalization(
					norm_table, // double [] direct_table, // last is <1.0, first > 0
					tex_hist_bins); //  int       num_bins)
			QuadCLTCPU.applyTexturesNormHist(
					faded_textures,  // final double [][][] textures, //  [nslices][nchn][i]
					cold_hot,        // final double []     min_max,
					inverted_table); //  final double []     inv_table)
		}
		
		
		if (debugLevel > -1) {
			double [][] dbg_textures = new double [faded_textures.length * faded_textures[0].length][faded_textures[0][0].length];
			String [] dbg_titles = new String[dbg_textures.length];
			String [] dbg_subtitles = new String [faded_textures[0].length];
			for (int i = 0; i < dbg_subtitles.length; i++) {
				dbg_subtitles[i] = (i <  (dbg_subtitles.length -1)) ? ("Y"+i):"alpha";
			}
			
			for (int i = 0; i < dbg_textures.length; i++) {
				dbg_textures[i] = faded_textures[i / faded_textures[0].length][i % faded_textures[0].length];
				dbg_titles[i] = dbg_subtitles[i % dbg_subtitles.length] + "-" + (i / dbg_subtitles.length);
				
			}
			ShowDoubleFloatArrays.showArrays(
					dbg_textures,
					tilesX * transform_size,
					tilesY * transform_size,
					true,
					ref_scene.getImageName()+"-combined_textures",
					dbg_titles);
			if (dbg_weights != null) {
				ShowDoubleFloatArrays.showArrays(
						dbg_weights,
						tilesX,
						tilesY,
						true,
						ref_scene.getImageName()+"-texture_weights");
			}
		}
		return faded_textures;
	}
	
	public static ImagePlus[] getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
			final CLTParameters  clt_parameters,
			boolean              no_alpha,
			QuadCLT              ref_scene,
			QuadCLT              parameter_scene,
			double [][][]        faded_textures,
			int                  tilesX,
			int                  tilesY,
			int                  transform_size,
			int                  debugLevel)
	{
		final boolean tex_color =           clt_parameters.tex_color;     //  true;  
		final int     tex_palette =         clt_parameters.tex_palette;     // 2 ;   
		int num_slices = faded_textures.length;
		if (parameter_scene == null) {
			parameter_scene = ref_scene;
		}
		double [] minmax = parameter_scene.getColdHot(); // used in linearStackToColor
		ImagePlus [] imp_tex = new ImagePlus[num_slices];
		for (int nslice = 0; nslice < num_slices; nslice++) {
            String title=String.format("%s-combo%03d-texture",ref_scene.getImageName(), nslice);
            double [][] rendered_texture = faded_textures[nslice].clone(); // shallow !
            if (no_alpha) {
            	rendered_texture[1] = new double [rendered_texture[0].length];
            	for (int i = 0; i < rendered_texture[0].length; i++) {
            		rendered_texture[1][i] = Double.isNaN(rendered_texture[0][i])? 0.0: 1.0;
            	}
            }
			imp_tex[nslice] =  	  QuadCLTCPU.linearStackToColorLWIR(
					clt_parameters, // CLTParameters  clt_parameters,
					tex_palette, // int            lwir_palette, // <0 - do not convert
					  minmax, // double []      minmax,
					  title, // String         name,
					  "", // String         suffix, // such as disparity=...
					  tex_color, // boolean        toRGB,
					  rendered_texture, // faded_textures[nslice], // double [][]    texture_data,
					  tilesX * transform_size, // int            width, // int tilesX,
					  tilesY * transform_size, // int            height, // int tilesY,
					  debugLevel); // int            debugLevel )
			// Add synthetic mesh only with higher resolution? or just any by a specified period?what king of mesh - vertical random, ...
			// Split and save as png
		}
		// Process accumulated textures: average, apply borders, convert to color or apply UM, add synthetic mesh, ... 
		return imp_tex; // ImagePlus[] ? with alpha, to be split into png and saved with alpha.
	}

	
	
	
	
	
	public static ImagePlus [] splitCombinedTextures(
			TileCluster [] tileClusters, //should have name <timestamp>-*
			int            transform_size,
			ImagePlus []   combo_textures ) {
		int max_cluster = -1;
		for (int nslice=0; nslice < tileClusters.length; nslice++) {
			for (int indx: tileClusters[nslice].getSubIndices()) {
				if (indx > max_cluster) max_cluster= indx;
			}
		}
		ImagePlus [] tex_clusters = new ImagePlus[max_cluster + 1];
		for (int nslice=0; nslice < tileClusters.length; nslice++) {
			String basename = combo_textures[nslice].getTitle();
			basename = basename.substring(0,basename.indexOf('-'));
			ImageStack combo_stack = combo_textures[nslice].getStack();
			int nSlices = combo_stack.getSize();
			int [] indices = tileClusters[nslice].getSubIndices();
			Rectangle [] bounds = tileClusters[nslice].getSubBounds();
			// try to deal with a single-slice stack?
			for (int i = 0; i < indices.length; i++) {
				Rectangle roi = bounds[i];
				ImageStack sub_stack = combo_stack.crop(
						roi.x * transform_size,
						roi.y * transform_size,
						0, // z - start slice
						roi.width * transform_size,
						roi.height* transform_size,
						nSlices);
				int cluster_index = indices[i];
	            String title=String.format("%s-img%04d-texture",basename, cluster_index);
				tex_clusters[cluster_index] = new ImagePlus(title, sub_stack);
//				tex_clusters[cluster_index].getProcessor().resetMinAndMax(); // probably not needed for png
			}
		}
		return tex_clusters;
	}
	
	
	
	// modifying so mono has 2 layers (Y + alpha) from IDC6585: public double [][] combineRBGATiles(
	/**
	 * Convert texture tiles  [tilesY][tilesX][2..4][4*transform_size] to image-like [2..4][width*height] 
	 *  
	 * @param texture_tiles [tilesY][tilesX][2..4][4*transform_size] {Y,alpha} or {r,b,g, alpha}
	 * @param overlap   false join full 16x16 texture tiles into twice larger image in each direction, true combine overlapping tiles
	 * @param sharp_alpha false: treat alpha same as Y/RBG, true - use alpha from the center 8x8
	 * @param transform_size now 8 everywhere
	 * @param num_channels 4 for RGBA, 2 for mono
	 * @param debugLevel
	 * @return [num_channels][width*height], where width is tilesX * transform_size for overlap=true and 2 * tilesX * transform_size otherwise 
	 */
	public static double [][] combineYRBGATiles(
			final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
			final boolean         overlap,    // when false - output each tile as 16x16, true - overlap to make 8x8
			final boolean         sharp_alpha, // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
			final int             transform_size, // 
			final int             num_channels, // 4 for RGBA, 2 for Y (should match textures)
			final int             debugLevel)
	{
		final int tilesY=texture_tiles.length;
		final int tilesX=texture_tiles[0].length;
		final int width=  (overlap?1:2)*tilesX * transform_size;
		final int height=  (overlap?1:2)*tilesY * transform_size;
		if (debugLevel > 1) {
			System.out.println("iclt_2d():tilesX=        "+tilesX);
			System.out.println("iclt_2d():tilesY=        "+tilesY);
			System.out.println("iclt_2d():width=         "+width);
			System.out.println("iclt_2d():height=        "+height);
			System.out.println("iclt_2d():overlap=       "+overlap);
			System.out.println("iclt_2d():sharp_alpha=   "+sharp_alpha);
		}
		final boolean is_mono = num_channels < 3;
		final int     alpha_chn = is_mono ? 1: 3;
///		final boolean has_alpha = num_channels > alpha_chn;
		final double [][] dpixels = new double[num_channels][width*height]; // assuming java initializes them to 0
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int n=0; n<4; n++){
			final int fn = n;
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						int tileY,tileX;
						int n2 = transform_size * 2;
						int n_half = transform_size / 2;
						int lastY = tilesY-1;
						int lastX = tilesX-1;
						int offset = n_half * (transform_size * tilesX) + n_half; // 4 pixels left and down (right/up when subtracted below)
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[fn].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[fn][nTile][0];
							tileY = tiles_list[fn][nTile][1];
							double [][] texture_tile =texture_tiles[tileY][tileX];
							if (texture_tile != null) {
								if (overlap) {
									if ((tileY >0) && (tileX > 0) && (tileY < lastY) && (tileX < lastX)) { // fast, no extra checks - ignore first/last rows and columns
										for (int i = 0; i < n2; i++){
											int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size - offset;
											for (int chn = 0; chn < texture_tile.length; chn++) {
												if (texture_tile[chn] == null) { // should never happen
													dpixels[chn] = new double [n2*n2];
												} else {
													// should it be better to multiply each color by alpha before accumulating? No, it is already windowed!
													if ((chn != alpha_chn) || !sharp_alpha) {
														for (int j = 0; j<n2;j++) {
															dpixels[chn][start_line + j] += texture_tile[chn][n2 * i + j];
														}
													} else if ((i >= n_half) && (i < (n2-n_half))) {
														for (int j = n_half; j < (n2 - n_half); j++) { // not used in lwir
															dpixels[chn][start_line + j] += texture_tile[chn][n2 * i + j];
														}
													}
												}
											}
										}
									} else { // be careful with margins
										for (int i = 0; i < n2;i++){
											if (	((tileY > 0) && (tileY < lastY)) ||
													((tileY == 0) && (i >= n_half)) ||
													((tileY == lastY) && (i < (n2 - n_half)))) {
												int start_line = ((tileY*transform_size + i) * tilesX + tileX)*transform_size  - offset;
												for (int chn = 0; chn < texture_tile.length; chn++) {
													if (texture_tile[chn] == null) {  // should never happen
														dpixels[chn] = new double [n2*n2];
													} else {
														if ((chn != alpha_chn) || !sharp_alpha) {
															for (int j = 0; j<n2;j++) {
																if (	((tileX > 0) && (tileX < lastX)) ||
																		((tileX == 0) && (j >= n_half)) ||
																		((tileX == lastX) && (j < (n2 - n_half)))) {
																	dpixels[chn][start_line + j] += texture_tile[chn][n2 * i + j];
																}
															}
														} else if ((i >= n_half) && (i < (n2-n_half))) {
															for (int j = n_half; j < (n2 - n_half); j++) { // not used in lwir
																if (	((tileX > 0) && (tileX < lastX)) ||
																		((tileX == 0) && (j >= n_half)) ||
																		((tileX == lastX) && (j < (n2 - n_half)))) {
																	dpixels[chn][start_line + j] += texture_tile[chn][n2 * i + j];
																}
															}
														}
													}
												}
											}
										}

									}
								} else { //if (overlap) - just copy tiles w/o overlapping
									for (int i = 0; i < n2;i++){ // not used in lwir
										for (int chn = 0; chn < texture_tile.length; chn++) {
											if (texture_tile[chn] == null) {  // should never happen
												dpixels[chn] = new double [n2*n2];
											} else {
												System.arraycopy(
														texture_tile[chn],
														i * n2,
														dpixels[chn],
														(tileY * n2 + i)* width + tileX*n2,
														n2);
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
		return dpixels;
	}
	
	/**
	 * Convert texture tiles, each having multiple non-overlapping clusters into YA or RBGA images,
	 * to be converted to pseudocolors or subject to unsharp mask together, then cut into individual
	 * rectangular areas using data in tileCluster
	 * 
	 * @param clt_parameters processing parameters
	 * @param texture_tiles  texture tiles as [tilesY][tilesX][2..4][256] or [tilesY][tilesX][][] 
	 * @param tileCluster    instance of TileCluster with individual bounds , disparities and borders
	 * @param sharp_alpha    combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
	 * @param transform_size always 8 now
	 * @param num_channels   2 for YA, 4 for RBGA
	 * @param debugLevel     debug level
	 * @return [2..4][width*height]
	 */
	public static double [][] getFadedTextures( // get image from a single pass, return relative path for x3d // USED in lwir
			CLTParameters         clt_parameters,
			final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
			final TileCluster     tileCluster,    // disparities, borders, selections for texture passes
			final boolean         sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
			final int             transform_size, // 
			final int             num_channels,   // 4 for RGBA, 2 for Y (should match textures)
			final int             debugLevel)
	{
		final int        tilesY =        texture_tiles.length;
		final int        tilesX =        texture_tiles[0].length;
		final boolean    is_mono =       num_channels < 3;
		final int        alpha_chn =     is_mono ? 1: 3;
		final int []     cluster_index = tileCluster.getClusterIndex();
		final boolean [] borderTiles =   tileCluster.getBorder(); // .getBorderTiles();
//		final boolean [] selected =      tileCluster.getSelected(); // maybe not needed
		final double [][]alphaFade =     TileProcessor.getAlphaFade(transform_size);
		final double []  alpha_zero =    new double [4*transform_size*transform_size];
		final double [][][][] texture_tiles_cluster = new double[tilesY][tilesX][][];

		for (int i = 0; i < alpha_zero.length; i++) alpha_zero[i]=0.0;
		// border tiles are copied, alpha from alphaFade (not multiplied?)
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				texture_tiles_cluster[tileY][tileX]= null;
				if (texture_tiles[tileY][tileX] != null) {
					int tile = tileY * tilesX + tileX;
					if (borderTiles[tile]) {
						int ci = cluster_index[tile];
						texture_tiles_cluster[tileY][tileX] = texture_tiles[tileY][tileX].clone();
						if (clt_parameters.shAggrFade) { // not used in lwir
							texture_tiles_cluster[tileY][tileX][alpha_chn] = alpha_zero;
						} else {
//							if ((debugLevel > -1) && (scanIndex == 1)) {
//								System.out.println("getPassImage(): tileY="+tileY+", tileX = "+tileX+", tileY="+tileY);
//							}
							int fade_mode = 0;
							int index_u = tile - tilesX;
							int index_r = tile + 1;
							int index_d = tile + tilesX;
							int index_l = tile - 1;
							if ((tileY > 0) &&           (texture_tiles[tileY - 1][tileX] != null) && !borderTiles[index_u] && (cluster_index[index_u] == ci)) fade_mode |= 1;
							if ((tileX < (tilesX -1)) && (texture_tiles[tileY][tileX + 1] != null) && !borderTiles[index_r] && (cluster_index[index_r] == ci)) fade_mode |= 2;
							if ((tileY < (tilesY -1)) && (texture_tiles[tileY + 1][tileX] != null) && !borderTiles[index_d] && (cluster_index[index_d] == ci)) fade_mode |= 4;
							if ((tileX > 0) &&           (texture_tiles[tileY][tileX - 1] != null) && !borderTiles[index_l] && (cluster_index[index_l] == ci)) fade_mode |= 8;
//							texture_tiles_cluster[tileY][tileX][alpha_chn] = alphaFade[fade_mode]; // alpha_zero;
							texture_tiles_cluster[tileY][tileX][alpha_chn] = texture_tiles_cluster[tileY][tileX][alpha_chn].clone(); // is it needed?
							for (int i = 0; i < texture_tiles_cluster[tileY][tileX][alpha_chn].length; i++) {
								texture_tiles_cluster[tileY][tileX][alpha_chn][i] *= alphaFade[fade_mode][i];
							}
						}
					}else{
						texture_tiles_cluster[tileY][tileX]= texture_tiles[tileY][tileX];
					}
				}
			}
		}
		
		double [][] texture_overlap = combineYRBGATiles(
				texture_tiles_cluster, // final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
				true,                  // final boolean         overlap,    // when false - output each tile as 16x16, true - overlap to make 8x8
				sharp_alpha,           // final boolean         sharp_alpha, // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
				transform_size,        // final int             transform_size, // 
				num_channels,          // final int             num_channels, // 4 for RGBA, 2 for Y (should match textures)
				debugLevel);           // final int             debugLevel)
		
		// Update alpha to sharpen "tree branches"
		if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
			double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
			for (int i = 0; i < texture_overlap[alpha_chn].length; i++){
				double d = texture_overlap[alpha_chn][i];
				if      (d >=clt_parameters.alpha1) d = 1.0;
				else if (d <=clt_parameters.alpha0) d = 0.0;
				else d = scale * (d- clt_parameters.alpha0);
				texture_overlap[alpha_chn][i] = d;
			}
		}
		return texture_overlap;
		/*
		// for now - use just RGB. Later add option for RGBA (?)
		double [][] texture_rgb = {texture_overlap[0],texture_overlap[1],texture_overlap[2]};
		double [][] texture_rgba = {texture_overlap[0],texture_overlap[1],texture_overlap[2],texture_overlap[3]};
		double [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);
		
		return  texture_rgbx;
		*/

		/*
		// Resize was extracting rectangle from the full size texture. With consolidated texture (9.12.2022) multiple
		// rectangles can be extracted from a single texture array
		boolean resize = true;
		if (resize) {
			texture_rgbx = resizeGridTexture(
					texture_rgbx,
					image_dtt.transform_size,
					tilesX,
					tilesY,
					scan.getTextureBounds());
		}

		int width = resize ? (image_dtt.transform_size * scan.getTextureBounds().width): (image_dtt.transform_size * tilesX);
		int height = resize ? (image_dtt.transform_size * scan.getTextureBounds().height): (image_dtt.transform_size * tilesY);
		if ((width <= 0) || (height <= 0)) {
			System.out.println("***** BUG in getPassImage(): width="+width+", height="+height+", resize="+resize+" ****"); // not used in lwir
		}
		// 9.12.2022: should be done for each cluster
		ImagePlus imp_texture_cluster = linearStackToColor(
				clt_parameters,
				colorProcParameters,
				rgbParameters,
				name+"-texture", // String name,
				"", //String suffix, // such as disparity=...
				true, // toRGB,
				!this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				true, // boolean saveShowIntermediate, // save/show if set globally
				false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				texture_rgbx,
				width, //tp.tilesX *  image_dtt.transform_size,
				height, //tp.tilesY *  image_dtt.transform_size,
				1.0,         // double scaleExposure, // is it needed?
				debugLevel);


		String path= correctionsParameters.selectX3dDirectory(
				//TODO: Which one to use - name or this.image_name ?
				correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				correctionsParameters.x3dModelVersion,
				true,  // smart,
				true);  //newAllowed, // save
		// only show/save original size if debug or debug_filters)
		eyesisCorrections.saveAndShow(
				imp_texture_cluster,
				path,
				correctionsParameters.png,
				clt_parameters.show_textures,
				-1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		return imp_texture_cluster.getTitle()+".png"; // imp_texture_cluster;
		*/
	}
	
	


}
