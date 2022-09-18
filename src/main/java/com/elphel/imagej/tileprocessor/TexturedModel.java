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
package com.elphel.imagej.tileprocessor;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.gpu.GpuQuad;
import com.elphel.imagej.gpu.TpTask;
import com.elphel.imagej.x3d.export.WavefrontExport;
import com.elphel.imagej.x3d.export.X3dOutput;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;

public class TexturedModel {
	public static final int THREADS_MAX = 100;  // maximal number of threads to launch
	public static final int TILE_EMPTY =          0; 
	public static final int TILE_BORDER =         1; 
	public static final int TILE_BORDER_FLOAT =   2; 
	public static final int TILE_CONFIRMED =      3; 
	public static final int TILE_CANDIDATE =      4; // not used 
	public static final int CLUSTER_NAN =        -2; // disparity is NaN 
	public static final int CLUSTER_UNASSIGNED =  -1; // not yet assinged (>=0 - cluster number)
//	
//	long                    startStepTime;
	
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
		final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
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
						cluster_list.size(), // (debug_index? cluster_list.size(): -1),
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
	
// Starting redoing
	public static boolean output3d( // USED in lwir
			CLTParameters                            clt_parameters,
			ColorProcParameters                      colorProcParameters,
			EyesisCorrectionParameters.RGBParameters rgbParameters,
			final QuadCLT                            parameter_scene, // to use for rendering parameters in multi-series sequences
            // if null - use reference scene 
			QuadCLT []                               scenes,
			double [][]                              combo_dsn_final, // null OK, will read file
//			final int                                threadsMax,  // maximal number of threads to launch
			final boolean                            updateStatus,
			final int                                debugLevel)
	{
		final boolean       batch_mode = clt_parameters.batch_run;
		final int           ref_index =  scenes.length - 1;
		final QuadCLT       ref_scene =  scenes[ref_index];
		final TileProcessor tp =         ref_scene.getTileProcessor();
		if (ref_scene.image_data == null){
			return false; // not used in lwir
		}
		double infinity_disparity = 	ref_scene.getGeometryCorrection().getDisparityFromZ(clt_parameters.infinityDistance);
		X3dOutput x3dOutput = null;
		WavefrontExport wfOutput = null; 
		long startStepTime=System.nanoTime();
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int transform_size = tp.getTileSize();
		final double     tex_disp_adiffo = clt_parameters.tex_disp_adiffo; // 0.35; // 0.3;  disparity absolute tolerance to connect in ortho directions 
		final double     tex_disp_rdiffo = clt_parameters.tex_disp_rdiffo; // 0.12; // 0.1;  disparity relative tolerance to connect in ortho directions
		final double     tex_disp_adiffd = clt_parameters.tex_disp_adiffd; // 0.6;  // 0.4;  disparity absolute tolerance to connect in diagonal directions
		final double     tex_disp_rdiffd = clt_parameters.tex_disp_rdiffd; // 0.18; // 0.12; disparity relative tolerance to connect in diagonal directions
		final double     tex_disp_fof =    clt_parameters.tex_disp_fof;    // 1.5;  // Increase tolerance for friend of a friend
		final double     tex_fg_bg =       clt_parameters.tex_fg_bg;       // 0.1;  // Minimal FG/BG disparity difference (NaN bg if difference from FG < this)
		final double     tex_distort =     clt_parameters.tex_distort;     // 0.1;  // Maximal texture distortion to accumulate multiple scenes
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
			sky_invert[i] =  !sky_tiles[i];
		}
		// re-load , should create quadCLTs[ref_index].dsi
		double [][] dls_fg = {
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_DISP],
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_LMA],
				combo_dsn_final[OpticalFlow.COMBO_DSN_INDX_STRENGTH]
		};
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
				tilesX, // final int          tilesX,
				ds_fg_bg, // final double [][]  disparities, // may have more layers
				sky_invert, // final boolean []   selected, // to remove sky (pre-filter by caller, like for ML?)
				tex_disp_adiffo, // final double       disp_adiffo,
				tex_disp_rdiffo, // final double       disp_rdiffo,
				tex_disp_adiffd, // final double       disp_adiffd,
				tex_disp_rdiffd, // final double       disp_rdiffd,
				tex_disp_fof, // final double       disp_fof,    // enable higher difference (scale) for fried of a friend 
				debugLevel); //1); //  2); // final int          debugLevel)
		double [][]     inter_weights = new double [tilesY][tilesX]; // per-tile texture weights for inter-scene accumulation;
		double [][][][] inter_textures= new double [tilesY][tilesX][][]; // [channel][256] - non-overlapping textures
		boolean [] scenes_sel = new boolean[scenes.length];
		//		for (int i = scenes.length - 10; i <  scenes.length; i++) { // start with just one (reference) scene
		for (int i = 0; i <  scenes.length; i++) { // start with just one (reference) scene
			scenes_sel[i] = true;
		}
		boolean        renormalize = true;// false - use normalizations from previous scenes to keep consistent colors
		ImagePlus[] combined_textures = getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
				clt_parameters,      // final CLTParameters  clt_parameters,
				colorProcParameters, // ColorProcParameters  colorProcParameters,
				rgbParameters,       // EyesisCorrectionParameters.RGBParameters rgbParameters,
				parameter_scene,     // final QuadCLT        parameter_scene, // to use for rendering parameters in multi-series sequences
				// if null - use reference scene 
				scenes,              // final QuadCLT []     scenes,
				scenes_sel,          // final boolean []     scenes_sel, // null or which scenes to process
				null,                // final boolean []     selection, // may be null, if not null do not  process unselected tiles
				tileClusters,        // final TileCluster [] tileClusters, // disparities, borders, selections for texture passes
				//				final int            margin,
				renormalize,  // final boolean        renormalize,  // false - use normalizations from previous scenes to keep consistent colors
				debugLevel);         // final int            debug_level)

		boolean save_full_textures = false; // true;
		EyesisCorrectionParameters.CorrectionParameters correctionsParameters = ref_scene.correctionsParameters;
		String x3d_dir= correctionsParameters.selectX3dDirectory(
				//TODO: Which one to use - name or this.image_name ?
				correctionsParameters.getModelName(ref_scene.getImageName()), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				correctionsParameters.x3dModelVersion,
				true,  // smart,
				true);  //newAllowed, // save
		
		if (save_full_textures) {
			for (int nslice = 0; nslice < combined_textures.length; nslice++) {
				EyesisCorrections.saveAndShow(
						combined_textures[nslice], // imp_texture_cluster,
						x3d_dir,
						correctionsParameters.png,
						false, // (nslice < 4), // clt_parameters.show_textures,
						-1, // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						1); //
			}
		}
		boolean save_individual_textures = true;
		if (save_individual_textures) {
			ImagePlus [] imp_textures = splitCombinedTextures(
					tileClusters,       // TileCluster [] tileClusters, //should have name <timestamp>-*
					transform_size,     // int            transform_size,
					combined_textures); // ImagePlus []   combo_textures )
			for (int i = 0; i < imp_textures.length; i++) if (imp_textures[i] != null) { // should not be
				EyesisCorrections.saveAndShow(
						imp_textures[i], // imp_texture_cluster,
						x3d_dir,
						correctionsParameters.png,
						false, // (nslice < 4), // clt_parameters.show_textures,
						-1, // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						1); //
			}
		}

		if (debugLevel > -100) {
			return true;
		}

		ArrayList<CLTPass3d> clt_3d_passes = new ArrayList<CLTPass3d>(); // will use its own passes
		int next_pass = 1; // temporarily
		/*

		  if (clt_parameters.remove_scans){
			  System.out.println("Removing all scans but the first(background) and the last to save memory");
			  System.out.println("Will need to re-start the program to be able to output differently");
			  CLTPass3d   latest_scan = tp.clt_3d_passes.get(tp.clt_3d_passes.size() - 1);
			  tp.trimCLTPasses(1);
			  tp.clt_3d_passes.add(latest_scan); // put it back
		  }
		  int next_pass = tp.clt_3d_passes.size(); //
		  // Create tasks to scan, have tasks, disparity and border tiles in tp.clt_3d_passes
		  tp.thirdPassSetupSurf( // prepare tile tasks for the second pass based on the previous one(s) // needs last scan
				  clt_parameters,
				  //FIXME: make a special parameter?
				  infinity_disparity, //0.25 * clt_parameters.bgnd_range, // double            disparity_far,
				  clt_parameters.grow_disp_max, // other_range, //double            disparity_near,   //
				  ref_scene.getGeometryCorrection(),
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  0); // final int         debugLevel)
		 */



		if (!batch_mode && (debugLevel > -1)) {
			tp.showScan(
					tp.clt_3d_passes.get(next_pass-1),   // CLTPass3d   scan,
					"after_pass3-"+(next_pass-1)); //String title)
		}
		String x3d_path = ref_scene.getX3dDirectory();
		// create x3d file
		if (clt_parameters.output_x3d) {
			x3dOutput = new X3dOutput(
					clt_parameters,
					ref_scene.correctionsParameters,
					ref_scene.getGeometryCorrection(),
					tp.clt_3d_passes);
		}
		if (clt_parameters.output_obj && (x3d_path != null)) {
			try {
				wfOutput = new WavefrontExport(
						x3d_path,
						ref_scene.correctionsParameters.getModelName(ref_scene.getImageName()),
						clt_parameters,
						ref_scene.correctionsParameters,
						ref_scene.getGeometryCorrection(),
						tp.clt_3d_passes);
			} catch (IOException e) {
				System.out.println("Failed to open Wavefront files for writing");
				// TODO Auto-generated catch block
				e.printStackTrace();
				// do nothing, just keep
			}
		}
		if (x3dOutput != null) {
			x3dOutput.generateBackground(clt_parameters.infinityDistance <= 0.0); // needs just first (background) scan
		}

		int bgndIndex = 0; // it already exists?
		CLTPass3d bgndScan = tp.clt_3d_passes.get(bgndIndex);
		//TODO make it w/o need for  bgndScan.texture as GPU will calculate texture right before output
		//use selection? or texture_selection instead?		  
		if (bgndScan.texture != null) { // TODO: same for the backdrop too
			if (clt_parameters.infinityDistance > 0.0){ // generate background as a billboard
				// grow selection, then grow once more and create border_tiles
				// create/rstore, probably not needed
				boolean [] bg_sel_backup = bgndScan.getSelected().clone();
				boolean [] bg_border_backup = (bgndScan.getBorderTiles() == null) ? null: bgndScan.getBorderTiles().clone();
				boolean [] bg_selected = bgndScan.getSelected();
				boolean [] border_tiles = bg_selected.clone();
				tp.growTiles(
						2,                   // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
						bg_selected,
						null); // prohibit

				for (int i = 0; i < border_tiles.length; i++){
					border_tiles[i] = !border_tiles[i] && bg_selected[i];
				}
				// update texture_tiles (remove what is known not to be background
				if (bgndScan.texture_tiles != null) { // for CPU
					for (int ty = 0; ty < tilesY; ty++){
						for (int tx = 0; tx < tilesX; tx++){
							if (!bg_selected[tx + tilesX*ty]){
								bgndScan.texture_tiles[ty][tx] = null; //
							}

						}
					}
				} else {//  for GPU
					for (int i = 0; i < bg_selected.length; i++) {
						if (!bg_selected[i]) {
							bgndScan.setTextureSelection(i,false);
						}
					}
				}

				//TODO2020: set texture_selection
				bgndScan.setBorderTiles(border_tiles);
				// limit tile_op to selection?
				// updates selection from non-null texture tiles
				String texturePath = ref_scene.getPassImage( // get image from a single pass - both CPU and GPU
						clt_parameters,
						colorProcParameters,
						rgbParameters,
						ref_scene.correctionsParameters.getModelName(ref_scene.getImageName())+"-img_infinity", // +scanIndex,
						bgndIndex,
						THREADS_MAX,  // maximal number of threads to launch
						updateStatus,
						batch_mode ? -5: debugLevel);
				if (texturePath != null) { // null if empty image
					double [] scan_disparity = new double [tilesX * tilesY];
					int indx = 0;
					//		  boolean [] scan_selected = scan.getSelected();
					for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
						scan_disparity[indx++] = infinity_disparity;
					}
					//			  tp.showScan(
					//    				  scan, // CLTPass3d   scan,
					//    				  "infinityDistance");

					boolean showTri = false; // ((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||((scanIndex - next_pass) == 73);
					try {
						ref_scene.generateClusterX3d(
								x3dOutput,
								wfOutput,  // output WSavefront if not null
								texturePath,
								"INFINITY", // id (scanIndex - next_pass), // id
								"INFINITY", // class
								bgndScan.getTextureBounds(),
								bgndScan.getSelected(), // selected,
								scan_disparity, // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
								clt_parameters.transform_size,
								clt_parameters.correct_distortions, // requires backdrop image to be corrected also
								showTri, // (scanIndex < next_pass + 1) && clt_parameters.show_triangles,
								infinity_disparity,  // 0.3
								clt_parameters.grow_disp_max, // other_range, // 2.0 'other_range - difference from the specified (*_CM)
								clt_parameters.maxDispTriangle);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
						return false;
					}
					// maybe not needed
					bgndScan.setBorderTiles(bg_border_backup);
					bgndScan.setSelected(bg_sel_backup);
				}
			}
		}

		// With GPU - do nothing here or copy selected -> texture_selection?
		for (int scanIndex = next_pass; scanIndex < tp.clt_3d_passes.size(); scanIndex++){
			if (debugLevel > 0){
				System.out.println("FPGA processing scan #"+scanIndex);
			}
			ref_scene.CLTMeasureTextures( // has GPU version - will just copy selection
					clt_parameters,
					scanIndex,
					THREADS_MAX,  // maximal number of threads to launch
					updateStatus,
					batch_mode ? -5: debugLevel);
		}

		for (int scanIndex = next_pass; (scanIndex < tp.clt_3d_passes.size()) && (scanIndex < clt_parameters.max_clusters); scanIndex++){ // just temporary limiting
			if (debugLevel > -1){
				System.out.println("Generating cluster images (limit is set to "+clt_parameters.max_clusters+") largest, scan #"+scanIndex);
			}
			String texturePath = ref_scene.getPassImage( // get image from a single pass
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					ref_scene.correctionsParameters.getModelName(ref_scene.getImageName())+"-img"+scanIndex,
					scanIndex,
					THREADS_MAX,  // maximal number of threads to launch
					updateStatus,
					batch_mode ? -5: debugLevel);
			if (texturePath == null) { // not used in lwir
				continue; // empty image
			}
			CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);

			// TODO: use new updated disparity, for now just what was forced for the picture
			double [] scan_disparity = new double [tilesX * tilesY];
			int indx = 0;
			for (int ty = 0; ty < tilesY; ty ++) for (int tx = 0; tx < tilesX; tx ++){
				scan_disparity[indx++] = scan.disparity[ty][tx];
			}
			if (clt_parameters.avg_cluster_disp){
				double sw = 0.0, sdw = 0.0;
				for (int i = 0; i< scan_disparity.length; i++){
					//					  if (scan.selected[i] && !scan.border_tiles[i]){
					if (scan.getSelected()[i] && !scan.getBorderTiles()[i]){
						double w = scan.disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX][i];
						sw +=w;
						sdw += scan_disparity[i]*w;
					}
				}
				sdw/=sw;
				for (int i = 0; i< scan_disparity.length; i++){
					scan_disparity[i] = sdw;
				}
			}
			boolean showTri = !batch_mode && (debugLevel > -1) && (((scanIndex < next_pass + 1) && clt_parameters.show_triangles) ||((scanIndex - next_pass) == 73));

			try {
				ref_scene.generateClusterX3d(
						x3dOutput,
						wfOutput,  // output WSavefront if not null
						texturePath,
						"shape_id-"+(scanIndex - next_pass), // id
						null, // class
						scan.getTextureBounds(),
						scan.getSelected(),
						scan_disparity, // scan.disparity_map[ImageDtt.DISPARITY_INDEX_CM],
						clt_parameters.transform_size,
						clt_parameters.correct_distortions, // requires backdrop image to be corrected also
						showTri, // (scanIndex < next_pass + 1) && clt_parameters.show_triangles,
						// FIXME: make a separate parameter:
						infinity_disparity, //  0.25 * clt_parameters.bgnd_range,  // 0.3
						clt_parameters.grow_disp_max, // other_range, // 2.0 'other_range - difference from the specified (*_CM)
						clt_parameters.maxDispTriangle);
			} catch (IOException e) {
				e.printStackTrace();
				return false;
			}
		}

		if ((x3d_path != null) && (x3dOutput != null)){
			x3dOutput.generateX3D(x3d_path+Prefs.getFileSeparator() + ref_scene.correctionsParameters.getModelName(ref_scene.getImageName())+".x3d");
		}
		if (wfOutput != null){
			wfOutput.close();
			System.out.println("Wavefront object file saved to "+wfOutput.obj_path);
			System.out.println("Wavefront material file saved to "+wfOutput.mtl_path);
		}

		// Save KML and ratings files if they do not exist (so not to overwrite edited ones), make them world-writable
		ref_scene.writeKml        (null, debugLevel);
		ref_scene.writeRatingFile (debugLevel);
		Runtime.getRuntime().gc();
		System.out.println("output3d(): generating 3d output files  finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-startStepTime),3)+" sec, --- Free memory25="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
		return true;
	}

	// get multi-scene textures (GPU-only)
	public static ImagePlus[] getInterCombinedTextures( // return ImagePlus[] matching tileClusters[], with alpha
			final CLTParameters  clt_parameters,
			ColorProcParameters  colorProcParameters,
			EyesisCorrectionParameters.RGBParameters rgbParameters,
			QuadCLT              parameter_scene, // to use for rendering parameters in multi-series sequences
			                                      // if null - use reference scene 
			final QuadCLT []     scenes,
			final boolean []     scenes_sel, // null or which scenes to process
			final boolean []     selection, // may be null, if not null do not  process unselected tiles
			final TileCluster [] tileClusters, // disparities, borders, selections for texture passes
			final boolean        renormalize,  // false - use normalizations from previous scenes to keep consistent colors
			final int            debug_level)
	{
		// TODO: ***** scenes with high motion blur also have high ERS to be corrected ! *****
		final int               ref_index = scenes.length -1;
		final QuadCLT ref_scene = scenes[ref_index];
		if (parameter_scene == null) {
			parameter_scene = ref_scene;
		}
		final int earliestScene = ref_scene.getEarliestScene(scenes);
		final ErsCorrection ers_reference = ref_scene.getErsCorrection();
		final int tilesX =        ref_scene.getTileProcessor().getTilesX();
		final int tilesY =        ref_scene.getTileProcessor().getTilesY();
		final int tiles =         tilesX * tilesY;
		final int transform_size= ref_scene.getTileProcessor().getTileSize();
		final int num_channels =  ref_scene.isMonochrome()?2:4;  //
		
		final boolean filter_bg =      true; // make a clt parameter?
		final boolean mb_en =          clt_parameters.imp.mb_en;
		final double  mb_tau =         clt_parameters.imp.mb_tau;      // 0.008;// time constant, sec
		final double  mb_max_gain =    clt_parameters.imp.mb_max_gain; // 5.0;  // motion blur maximal gain (if more - move second point more than a pixel

		final double  max_distortion = clt_parameters.tex_distort;   // 0.5;  // Maximal texture distortion to accumulate multiple scenes (0 - any)
		final double  tex_mb =         clt_parameters.tex_mb;        // 1.0;  // Reduce texture weight if motion blur exceeds this (as square of MB length)
		final boolean sharp_alpha =    clt_parameters.sharp_alpha;
		final boolean is_lwir =        ref_scene.isLwir();
		
		final boolean tex_um =           clt_parameters.tex_um;        // imp.um_mono; // TODO: add own parameter
		final double  tex_um_sigma =     clt_parameters.tex_um_sigma;  // imp.um_sigma;
		final double  tex_um_weight =    clt_parameters.tex_um_weight; // imp.um_weight;
		// TODO: - make texture variants, tex_um_fixed/tex_um_range apply only to unsharp mask, regardless of colors
		
		final boolean lwir_autorange =   is_lwir && clt_parameters.tex_lwir_autorange; // colorProcParameters.lwir_autorange;
		final boolean tex_um_fixed =     clt_parameters.tex_um_fixed;  // imp.mono_fixed; //  true; // normalize to fixed range when converting to 8 bits 
		final double  tex_um_range =     clt_parameters.tex_um_range;  // imp.mono_range; // 500.0;  // monochrome full-scale range (+/- half)
		final boolean tex_hist_norm =    clt_parameters.tex_hist_norm; //  true;  
		final double  tex_hist_amount =  clt_parameters.tex_hist_amount; // clt_parameters. 0.7;  
		final int     tex_hist_bins =    clt_parameters.tex_hist_bins;   //  1024 ;   
		final int     tex_hist_segments =clt_parameters.tex_hist_segments; // 32 ;   

		final boolean tex_color =        clt_parameters.tex_color;     //  true;  
		final int     tex_palette =      clt_parameters.tex_palette;     // 2 ;   
		
		
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
			ref_scene.getGPU().setGpu_debug_level(debug_level);
		}
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		final int num_slices = tileClusters.length;
		double [][][]     inter_weights = new double [num_slices][tilesY][tilesX]; // per-tile texture weights for inter-scene accumulation;
		double [][][][][] inter_textures= new double [num_slices][tilesY][tilesX][][]; // [channel][256] - non-overlapping textures
///		double [][] scene_pXpYD;
///		final double disparity_corr = 0.00; // (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
///		TpTask[] tp_tasks_ref = null;
		double [][][] ref_pXpYDs = new double [num_slices][][]; // individual for each slice
		int    [][] cluster_indices = (max_distortion > 0.0) ? (new int [num_slices][]): null;
		boolean [][] borders = new boolean [num_slices][];
		for (int nslice = 0; nslice < num_slices; nslice++) { // prepare and measure textures for each combo textures
			ref_pXpYDs[nslice] = OpticalFlow.transformToScenePxPyD( // now should work with offset ref_scene
					null, // fov_tiles,            // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
					tileClusters[nslice].getDisparity(), // final double []   disparity_ref, // invalid tiles - NaN in disparity
					OpticalFlow.ZERO3,                // final double []   scene_xyz, // camera center in world coordinates
					OpticalFlow.ZERO3,                // final double []   scene_atr, // camera orientation relative to world frame
					scenes[ref_index],  // final QuadCLT     scene_QuadClt,
					scenes[ref_index],  // final QuadCLT     reference_QuadClt, // now - may be null - for testing if scene is rotated ref
					THREADS_MAX);          // int               threadsMax)
			borders[nslice] = tileClusters[nslice].getBorder();
			if (max_distortion > 0.0) {
				cluster_indices[nslice] = tileClusters[nslice].getClusterIndex();
			}
		}		
		
		
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
			if (mb_en) {
				dxyzatr_dt = OpticalFlow.getVelocities( // looks at previous/next scene poses
						scenes,   // QuadCLT []     quadCLTs,
						nscene);  // int            nscene)
			}
			scenes[nscene].saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU)
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
						debug_level) : null;        // int            debug_level)
				if (debug_level > 0) {
					System.out.println("nscene="+nscene+", nslice="+nslice+" will run texturesGPUFromDSI() that needs debug >2");
					System.out.print("");
				}
				double [][][][] slice_texture = QuadCLT.texturesGPUFromDSI(
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
						debug_level);            // final int         debugLevel);
				if (slice_texture != null) {
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
									if (slice_texture[tileY][tileX] != null) {
										double w = 1.0;
										if (tex_mb > 0.0) {
											double mb_l2 = mb_tau2 * ( motion_blur[0][nTile]*motion_blur[0][nTile] +
													motion_blur[1][nTile]*motion_blur[1][nTile]);
											if (mb_l2 > 1.0) {
												w /= mb_l2; // 1/(squared mb)
											}
										}
										if (w > 0) {
											inter_weights[fnslice][tileY][tileX] +=w;
											if (inter_textures[fnslice][tileY][tileX] == null) { // create if it did not exist
												inter_textures[fnslice][tileY][tileX] = new double [slice_texture[tileY][tileX].length][slice_texture[tileY][tileX][0].length];
											}
											for (int nchn = 0; nchn < inter_textures[fnslice][tileY][tileX].length; nchn++) {
												for (int i = 0; i < inter_textures[fnslice][tileY][tileX][nchn].length; i++) {
													inter_textures[fnslice][tileY][tileX][nchn][i] += w *slice_texture[tileY][tileX][nchn][i]; 	
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

				if (debug_level > -1) {
					if (nscene == ref_index) {
						System.out.println("Textures from the reference scene, nslice = " + nslice +((slice_texture == null)? " - EMPTY":""));
					} else {
						System.out.println("Textures from scene "+nscene+", slice="+nslice +((slice_texture == null)? " - EMPTY":""));
					}
				}
			} // for (int nslice = 0; nslice < num_slices; nslice++) {
		} // for (int nscene = 0; nscene < num_scenes; nscene++) {
		
		double [][][] faded_textures = new double [num_slices][][];
		final double [][] dbg_weights = (debug_level > 0 )?(new double [num_slices][tiles]) : null;
		final double [][] dbg_overlap = (debug_level > 0 )?(new double [num_slices*num_channels][]) : null;
		for (int nslice = 0; nslice < num_slices; nslice++) {
			final int fnslice = nslice;
			if (dbg_weights != null) {
				Arrays.fill(dbg_weights[nslice], Double.NaN);
			}
			final Thread[] threads = ImageDtt.newThreadArray(THREADS_MAX);
			final AtomicInteger ai = new AtomicInteger(0);
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
								for (int nchn = 0; nchn < inter_textures[fnslice][tileY][tileX].length; nchn++) {
									for (int i = 0; i < inter_textures[fnslice][tileY][tileX][nchn].length; i++) {
										inter_textures[fnslice][tileY][tileX][nchn][i] *= w; 	
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
    // Process slice of textures: apply borders, convert to color or apply UM, add synthetic mesh, ...			
	// 2 layers for mono, 4 layers - for color		
	// First - merge overlapped and apply borders, alpha is the last slice
			faded_textures[nslice] = getFadedTextures( // get image from a single pass, return relative path for x3d // USED in lwir
					clt_parameters,          // CLTParameters         clt_parameters,
					inter_textures[fnslice], // final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
					tileClusters[fnslice],   // final TileCluster     tileCluster,    // disparities, borders, selections for texture passes
					sharp_alpha,             // final boolean         sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
					transform_size,          // final int             transform_size, // 
					num_channels,            // final int             num_channels,   // 4 for RGBA, 2 for Y (should match textures)
					debug_level);            // final int             debugLevel)
			
			if (dbg_overlap != null) {
				double [][] non_overlap =  combineYRBGATiles(
						inter_textures[fnslice], // final double [][][][] texture_tiles,  // array [tilesY][tilesX][4][4*transform_size] or [tilesY][tilesX]{null}
						false, // final boolean         overlap,    // when false - output each tile as 16x16, true - overlap to make 8x8
						sharp_alpha, // final boolean         sharp_alpha, // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
						transform_size, // final int             transform_size, // 
						num_channels,            //final int             num_channels, // 4 for RGBA, 2 for Y (should match textures)
						debug_level); // final int             debugLevel)
				for (int i = 0; i < num_channels; i++) {
				  dbg_overlap[fnslice*num_channels+i] = non_overlap[i];
				}
			}
		}
		// Optionally apply UM (before auto/manual range)
		if (tex_um) {
			QuadCLTCPU.umTextures(
					faded_textures, // final double [][][] textures, //  [nslices][nchn][i]
					tilesX * transform_size, // final int    width,
					tex_um_sigma, // final double um_sigma,
					tex_um_weight); // final double um_weight)
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
							faded_textures ); //double [][][] textures //  [nslices][nchn][i]
					rel_low =  minmax[0];
					rel_high = minmax[1];
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
		
		
		if (debug_level > -1) {
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
			
			if (dbg_overlap != null) {
				ShowDoubleFloatArrays.showArrays(
						dbg_overlap,
						2 * tilesX * transform_size,
						2 * tilesY * transform_size,
						true,
						ref_scene.getImageName()+"-non-overlap_textures",
						dbg_titles);
			}
			if (dbg_weights != null) {
				ShowDoubleFloatArrays.showArrays(
						dbg_weights,
						tilesX,
						tilesY,
						true,
						ref_scene.getImageName()+"-texture_weights");
			}
		}

		double [] minmax = parameter_scene.getColdHot(); // used in linearStackToColor
		ImagePlus [] imp_tex = new ImagePlus[num_slices];
		for (int nslice = 0; nslice < num_slices; nslice++) {
            String title=String.format("%s-combo%03d-texture",ref_scene.getImageName(), nslice);
			imp_tex[nslice] =  	  QuadCLTCPU.linearStackToColorLWIR(
					clt_parameters, // CLTParameters  clt_parameters,
					tex_palette, // int            lwir_palette, // <0 - do not convert
					  minmax, // double []      minmax,
					  title, // String         name,
					  "", // String         suffix, // such as disparity=...
					  tex_color, // boolean        toRGB,
					  faded_textures[nslice], // double [][]    texture_data,
					  tilesX * transform_size, // int            width, // int tilesX,
					  tilesY * transform_size, // int            height, // int tilesY,
					  debug_level); // int            debugLevel )
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
		for (int nscene=0; nscene < tileClusters.length; nscene++) {
			for (int indx: tileClusters[nscene].getSubIndices()) {
				if (indx > max_cluster) max_cluster= indx;
			}
		}
		ImagePlus [] tex_clusters = new ImagePlus[max_cluster + 1];
		for (int nscene=0; nscene < tileClusters.length; nscene++) {
			String basename = combo_textures[nscene].getTitle();
			basename = basename.substring(0,basename.indexOf('-'));
			ImageStack combo_stack = combo_textures[nscene].getStack();
			int nSlices = combo_stack.getSize();
			int [] indices = tileClusters[nscene].getSubIndices();
			Rectangle [] bounds = tileClusters[nscene].getSubBounds();
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
		final boolean [] selected =      tileCluster.getSelected(); // maybe not needed
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
