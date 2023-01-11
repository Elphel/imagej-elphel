package com.elphel.imagej.x3d.export;

import java.awt.Rectangle;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.tileprocessor.GeometryCorrection;
import com.elphel.imagej.tileprocessor.ImageDtt;
import com.elphel.imagej.tileprocessor.TexturedModel;
import com.elphel.imagej.tileprocessor.TileNeibs;

/**
 **
 ** TriMesh - triangular mesh representation
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TriMesh.java is free software: you can redistribute it and/or modify
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

public class TriMesh {
	public static final int       TRI_DOWN_LEFT =      0; // down, then left (independent)
	public static final int       TRI_RIGHT_DOWNLEFT = 1; // right, then down-left (one of this and two next)
	public static final int       TRI_DOWNRIGHT_LEFT = 2; // down-right, then left
	public static final int       TRI_RIGHT_DOWN =     3; // right, then down
	public static final int[]     TRI_NONE =   {-1,-1,-1,-1}; // four possible triangles
	public static final int[][]   TRI_OFFS_XY = {{0,0}, {1,0}, {1,1}, {0,1}}; // X and Y relative to top-left
	
	public static final int[][][] TRI_SET_xy = {
			{{0,0},{0,1},{-1,1}}, // TRI_DOWN_LEFT
			{{0,0},{1,0},{ 0,1}}, // TRI_RIGHT_DOWNLEFT
			{{0,0},{1,1},{ 0,1}}, // TRI_DOWNRIGHT_LEFT
			{{0,0},{1,0},{ 1,1}}};// TRI_RIGHT_DOWN 
	
	public String texture_image;
	public double [][] worldXYZ;
	public double [][] texCoord;
	public int [][]    triangles;
	public TriMesh (
			String texture_image,
			double [][] worldXYZ,
			double [][] texCoord,
			int [][] triangles) {
		this.texture_image = texture_image;
		this.worldXYZ = worldXYZ;
		this.texCoord = texCoord;
		this.triangles = triangles;
	}
	public String getImage() {return texture_image;}
	int [][] getTriangles()  {return triangles;}
	double [][] getTexCoord() {
		return texCoord;
	}
	public double [][] getTexCoord(boolean inv_x, boolean inv_y, boolean swap_xy) {
		if (!inv_x && !inv_y && !swap_xy) {
			return texCoord;
		}
		double [][] inv_tex_coord = new double [texCoord.length][2];
		double scale_x = inv_x ? -1.0 : 1.0;
		double scale_y = inv_y ? -1.0 : 1.0;
		if (swap_xy) {
			for (int i = 0; i <texCoord.length; i++) {
				inv_tex_coord[i][0] = scale_y * texCoord[i][1];
				inv_tex_coord[i][1] = scale_x * texCoord[i][0];
			}
		} else {
			for (int i = 0; i <texCoord.length; i++) {
				inv_tex_coord[i][0] = scale_x * texCoord[i][0];
				inv_tex_coord[i][1] = scale_y * texCoord[i][1];
			}
		}
		return inv_tex_coord;	
	}
	public double [][] getCoordinates(){
		return worldXYZ;
	}
	/**
	 * 0: XYZ -> XYZ
	 * 1: XYZ -> YZX
	 * 2: XYZ -> ZXY
	 * 3: XYZ -> XZY
	 * 4: XYZ -> ZYX
	 * 5: XYZ -> YXZ
	 * @param inv_x
	 * @param inv_y
	 * @param inv_z
	 * @param swap3
	 * @return
	 */
	public double [][] getCoordinates(boolean inv_x, boolean inv_y, boolean inv_z, int swap3) {
		if (!inv_x && !inv_y && !inv_y && (swap3 == 0)) {
			return worldXYZ;
		}
		double [][] inv_worldXYZ = new double [worldXYZ.length][3];
		double scale_x = inv_x ? -1.0 : 1.0;
		double scale_y = inv_y ? -1.0 : 1.0;
		double scale_z = inv_z ? -1.0 : 1.0;
		switch (swap3) {
		case 0: // XYZ -> XYZ
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][0] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][1] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][2] = scale_z * worldXYZ[i][2];
			}
			break;
		case 1: // XYZ -> YZX
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][1] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][2] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][0] = scale_z * worldXYZ[i][2];
			}
			break;
		case 2: // XYZ -> ZXY
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][2] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][0] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][1] = scale_z * worldXYZ[i][2];
			}
			break;
		case 3: // XYZ -> XZY
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][0] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][2] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][1] = scale_z * worldXYZ[i][2];
			}
			break;
		case 4: // XYZ -> ZYX
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][2] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][1] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][0] = scale_z * worldXYZ[i][2];
			}
			break;
		case 5: // XYZ -> YXZ
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][1] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][0] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][2] = scale_z * worldXYZ[i][2];
			}
			break;
		default: return null;
		}
		return inv_worldXYZ;
	}
	
	
// moved here from TileProcessor
	/**
	 * Enumerate selected tiles. Input boolean array maybe either full size (tilesY*tilesX)
	 * or bounds.width*bounds.height
	 * @param bounds Rectangle specifing selection area 
	 * @param selected boolean array of populated tiles
	 * @param tilesX full width of image array
	 * @return [y][x] array of incremental indices, -1 for unselected tiles.
	 */
	public static int [][] getCoordIndices( // starting with 0, -1 - not selected
			Rectangle bounds,
			boolean [] selected,
			int        tilesX)
	{
		int [][] indices = new int [bounds.height][bounds.width];
		int indx = 0;
		if (selected.length > (bounds.height * bounds.width)) { // old version - selected is full size
			for (int y = 0; y < bounds.height; y++) {
				for (int x = 0; x < bounds.width; x++){
					if (selected[tilesX * (bounds.y + y) + (bounds.x + x)]){
						indices[y][x] = indx++;
					} else {
						indices[y][x] = -1;
					}
				}
			}
		} else { // 09.18.2022
			for (int y = 0; y < bounds.height; y++) {
				for (int x = 0; x < bounds.width; x++){
					if (selected[bounds.width * y + x]){
						indices[y][x] = indx++;
					} else {
						indices[y][x] = -1;
					}
				}
			}
		}
		return indices;
	}
	
	/**
	 * Enumerate "large" and "small" tiles, where "large" are actual tiles and "small" are
	 * subdivided (by subdiv in each direction) ones to increase lateral mesh resolution.
	 * sub-tiles are populated if at least one pixel in it is opaque. Input selected_tiles
	 * and alpha arrays may correspond to either full image or rectangular bounds, output
	 * array always corresponds to bounds.
	 * 
	 * @param bounds         Rectangle ROI for output and optionally input data
	 * @param selected_tiles selected tiles, all unselected are ignored
	 * @param tilesX         full image width in tiles
	 * @param tile_size      tile size (8)
	 * @param alpha          pixel array, where true means "opaque". may be either full image
	 *                       or be bound to bounds (scaled to pixels from tiles)
	 * @param subdiv         subdivide tiles. Best if is equal to 1,2,4 and 8 that results in
	 *                       uniform tiles 
	 * @param num_indices    return parameter if int[1] - will provide total number of selected
	 *                       tiles - large and small.
	 * @return               array with tile indices (all different) int [bounds.height][bounds.width][][]
	 *                       Each "large" tile may be null if empty, contain a single int [][] {{index}}
	 *                       for "large" tiles or [subdiv][subdiv] array of the sub-tile indices with
	 *                       -1 for empty subtile.  
	 */
	

	public static int [][][][] getCoordIndices( // starting with 0, -1 - not selected
			final Rectangle  bounds,
			final boolean [] selected_tiles,
			final int        tilesX,
			final int        tile_size,		
			final boolean [] alpha,
			final int        subdiv,
			final int []     num_indices
			) {
		//TODO: add optimization (merging some in-tile triangles)
		final boolean full_selection = selected_tiles.length > (bounds.height * bounds.width); // applies to selected_tiles 
		final boolean full_alpha = alpha.length > (bounds.height * bounds.width * tile_size * tile_size);
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final int source_tile_width =  full_selection? tilesX : bounds.width;
		final int source_tile_offsx =  full_selection? bounds.x : 0;
		final int source_tile_offsy =  full_selection? bounds.y : 0;
		final int source_pix_offsx = (full_alpha? bounds.x : 0) * tile_size;
		final int source_pix_offsy = (full_alpha? bounds.y : 0) * tile_size;
		final int source_pix_width = (full_alpha? tilesX : bounds.width) * tile_size;
		final boolean [][][][] nodes = new boolean [bounds.height][bounds.width][][]; // selected nodes to be enumerated
		final boolean [][] emty_alpha = new boolean [bounds.height][bounds.width]; // selected tile with no opaque subtiles
		
		final int btiles = bounds.width * bounds.height;
		
		ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                public void run() {
                	boolean [][] pix_sel = new boolean[tile_size][tile_size];
					for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
						int btilex = btile % bounds.width;
						int btiley = btile / bounds.width;
						int stile = (btilex + source_tile_offsx) + (btiley + source_tile_offsy) * source_tile_width;
						if (selected_tiles[stile]) {
							boolean all_sel = true, some_sel = false;
							int pindx0 = (source_pix_offsx + btilex * tile_size) + (source_pix_offsy + btiley * tile_size) * source_pix_width;
							for (int y = 0; y < tile_size; y++) {
								int pindx1 = pindx0 + y*source_pix_width;
								for (int x = 0; x < tile_size; x++) {
									int pindx = pindx1 + x;
									boolean psel = alpha[pindx];
									pix_sel[y][x] = psel;
									all_sel &= psel;
									some_sel |= psel;
								}								
							}
							if (some_sel) {
								if (all_sel) {
//									indices[btiley][btilex] = new int [][] {{aindx.getAndIncrement()}}; // single center index
									nodes[btiley][btilex] = new boolean [][] {{true}}; // single center index
								} else {
//									indices[btiley][btilex] = new int [subdiv][subdiv];
									nodes[btiley][btilex] = new boolean [subdiv][subdiv];
									for (int ny = 0; ny < subdiv; ny++) {
										int y0 = (ny * tile_size)/subdiv;
										int y1 = Math.min(((ny+1) * tile_size)/subdiv, tile_size);
										for (int nx = 0; nx < subdiv; nx++) {
											int x0 = (nx * tile_size)/subdiv;
											int x1 = Math.min(((nx+1) * tile_size)/subdiv, tile_size);
											boolean has_any = false;
											lhas_any:
												for (int y = y0; y < y1; y++) {
													for (int x = x0; x < x1; x++) {
														if (pix_sel[y][x]) {
															has_any = true;
															break lhas_any;
														}
													}
												}
											nodes[btiley][btilex][ny][nx] = has_any;
										}										
									}
								}
							} else {
								emty_alpha[btiley][btilex] = true;
							}
						}
					}
                }
            };
        }		      
        ImageDtt.startAndJoin(threads);
        // split neighbors of emty_alpha. Direction in outer cycle to prevent thread conflicts
        for (int dir = 0; dir < TileNeibs.DIR_S; dir++) {
        	final int fdir = dir;
    		ai.set(0);
            for (int ithread = 0; ithread < threads.length; ithread++) {
                threads[ithread] = new Thread() {
                    public void run() {
    					for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
    						int btilex = btile % bounds.width;
    						int btiley = btile / bounds.width;
    						if (emty_alpha[btiley][btilex]) { // empty_alpha a rare, this is why WE iterate over them 
    							int btilex1 = btilex +  TileNeibs.DIR_XY[fdir][0];
    							int btiley1 = btiley +  TileNeibs.DIR_XY[fdir][1];
    							if ((btilex1 >= 0) && (btiley1 >= 0) &&
    									(btilex1 < bounds.width) && (btiley1 < bounds.height) &&
    									(nodes[btiley1][btilex1] != null) && (nodes[btiley1][btilex1].length == 1)) {
    								// split single-tile
									nodes[btiley1][btilex1] = new boolean [subdiv][subdiv];
									for (int ny = 0; ny < subdiv; ny++) {
										Arrays.fill(nodes[btiley1][btilex1][ny], true);
									}
    							}
    						}
    					}
                    }
                };
            }		      
            ImageDtt.startAndJoin(threads);
        }
        
        // Assign indices to selected tiles/subtiles
		final AtomicInteger aindx = new AtomicInteger(0);
		final int [][][][] indices = new int [bounds.height][bounds.width][][];
		ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                public void run() {
					for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
						int btilex = btile % bounds.width;
						int btiley = btile / bounds.width;
						if (nodes[btiley][btilex] != null) {
							if (nodes[btiley][btilex].length == 1) {
								indices[btiley][btilex] = new int [][] {{aindx.getAndIncrement()}}; // single center index
							} else {
								indices[btiley][btilex] = new int [subdiv][subdiv];
								for (int ny = 0; ny < subdiv; ny++) {
									for (int nx = 0; nx < subdiv; nx++) {
										if (nodes[btiley][btilex][ny][nx]) {
											indices[btiley][btilex][ny][nx] = aindx.getAndIncrement();
										} else {
											indices[btiley][btilex][ny][nx] = -1;
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
        if (num_indices != null) {
        	num_indices[0] = aindx.get();
        }
        return indices;
	}

	/**
	 * Designed to prevent bridging between neighbor tiles over disparity discontinuity gap. 
	 * neib_lev==2 is the outermost FG tile, while neib_lev3 - the outermost of the BG. For the FG
	 * levels go from 0 (inner tiles) 0-1-2, for the BG - 0-1-3, so tiles 3 should not be connected to 2.
	 * Mark border tiles that should not be connected by triangles. Tiles with border_int[] of (max_border)
	 * should not be connected to (max_border+1). THey are marked with 1 and 2 (all others are 0). 
	 * @param bounds      Rectangle ROI for output and optionally input data
	 * @param indices     array of [height][width]{{index}} for large tiles and [height][width][py][px]
	 *                    for small ones. This array will be modified and re-indexed if needed.
	 * @param border_int  border values array (same dimensions as disparity and selected) -1 - unassigned,
	 *                    0 - not a border, 1 inner border, 2 and 3 (for max_border==2) are both outer borders,
	 *                    but they are for different "leaves" and should not be meshed 
	 * @param max_border  maximal border_int value - now 2
	 * @param tilesX      full image width in tiles
	 * @param tile_size   tile size (8)
	 * @return            2D array corresponding to top indices. Value 0 - nop, 1 and 2 should not be connected
	 *                    by any triangle.
	 */
	public static int [][] getNoConnect( // 0 - neutral 1 - max_neib_lev, 2 - max_neib_lev+1
			final Rectangle  bounds,
			final int [][][][] indices,
			final int []     border_int,
			final int        max_border,
			final int        tilesX,
			final int        tile_size		
			) {
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final boolean full_selection = border_int.length > (bounds.height * bounds.width); // applies to selected_tiles 
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final int source_tile_width =  full_selection? tilesX : bounds.width;
		final int source_tile_offsx =  full_selection? bounds.x : 0;
		final int source_tile_offsy =  full_selection? bounds.y : 0;
		final int btiles = bounds.width * bounds.height;
		final int [][] no_connect = new int [bheight][bwidth];
        ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				if (indices[btiley][btilex] != null) {
    						int stile = (btilex + source_tile_offsx) + (btiley + source_tile_offsy) * source_tile_width;
    						if (border_int[stile] >= max_border) {
    							no_connect[btiley][btilex] = border_int[stile] - max_border + 1;
    						}
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
		return no_connect;
	}	
	
	/**
	 * Convert "large" tiles to arrays of small ones if it has a small-tile neighbor with
	 * gaps along the border with this one
	 * @param indices -   array of [height][width]{{index}} for large tiles and [height][width][py][px]
	 *                    for small ones. This array will be modified and re-indexed if needed.
	 * @param subdiv      subdivide tiles. Best if is equal to 1,2,4 and 8 that results in
	 *                    uniform tiles 
	 * @return            number of indices after re-indexing.                   
	 */
	public static int splitLargeTileIndices(
			int [][][][] indices,
			final int    subdiv) {
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final boolean [][] need_split = new boolean[bheight][bwidth];
		final int btiles = bwidth * bheight;
		final int [][] dxy = {{0, -1},{1,0},{0,1},{-1,0}};
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				if ((indices[btiley][btilex] != null) && (indices[btiley][btilex].length == 1)) { // full tile, not split
        					l_split:
        						for (int dir = 0; dir < dxy.length; dir++) {
        							int btilex1 = btilex+dxy[dir][0];
        							int btiley1 = btiley+dxy[dir][1];
        							if ((btilex1 >= 0) && (btiley1 >= 0) && (btilex1 < bwidth) && (btiley1 < bheight)) {
        								int [][] index_xy= indices[btiley1][btilex1];
        								if ((index_xy != null) && (index_xy.length >1)) { // split tile
        									int x0 = (dxy[dir][0] > 0) ? 0: (index_xy[0].length - 1);
        									int y0 = (dxy[dir][1] > 0) ? 0: (index_xy.length - 1);
        									if (dxy[dir][0] == 0) { // top or bottom row
        										for (int x = 0; x < index_xy[0].length; x++) {
        											if (index_xy[y0][x] < 0) {
        												need_split[btiley][btilex] = true;
        												break l_split;
        											}
        										}
        									} else {//  right or left column 
        										for (int y = 0; y < index_xy.length; y++) {
        											if (index_xy[y][x0] < 0) {
        												need_split[btiley][btilex] = true;
        												break l_split;
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
        // re-index
		final AtomicInteger aindx = new AtomicInteger(0);
        ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
						int [][] index_xy= indices[btiley][btilex];
						if (index_xy != null) {
							if (need_split[btiley][btilex]) {
								indices[btiley][btilex] = new int [subdiv][subdiv];
								for (int y = 0; y < subdiv; y++) {
									for (int x = 0; x < subdiv; x++) {
										indices[btiley][btilex][y][x] = aindx.getAndIncrement();
									}								
								}
							} else {
								for (int y = 0; y < index_xy.length; y++) {
									for (int x = 0; x < index_xy[y].length; x++) {
										if (index_xy[y][x] >= 0) {
											index_xy[y][x] = aindx.getAndIncrement(); // re-index any non-negative
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
		return aindx.get();
	}

	/**
	 * Assuming only neighbor tiles
	 * @param x
	 * @param y
	 * @param size
	 * @return
	 */
	static int [] getNeibNode(int x, int y, int size) {
		int dir = 8;
		if (x < 0) {
			if (y < 0) {
				x += size;
				y += size;
				dir = 7;
			} else if (y >= size) {
				x += size;
				y -= size;
				dir = 5;
			} else {
				x += size;
				dir = 6;
			}
		} else if (x >= size) {
			if (y < 0) {
				x -= size;
				y += size;
				dir = 1;
			} else if (y >= size) {
				x -= size;
				y -= size;
				dir = 3;
			} else {
				x -= size;
				dir = 2;
			}
		} else {
			if (y < 0) {
				y += size;
				dir = 0;
			} else if (y >= size) {
				y -= size;
				dir = 4;
			} else {
				dir=8;
			}
		}
		return new int[] {x, y, dir};
		 
	}

	/**
	 * get edge indices, CCW: 0 - top edge, 1 right edge, 2 - bottom edge, 3 - left edge 
	 * @param tile square tile [y][x]
	 * @param dir direction/edge 0..3
	 * @return edge values of the tile in counter-clockwise order
	 */
	public static int [] getEdgeIndices(int [][] tile, int dir) {
		if (tile == null) return null;
		int subdiv = tile.length;
		int [] edge = new int[subdiv];
		switch (dir) {
		case 0: for (int i = 0; i < subdiv; i++) {edge[i] = tile[0][subdiv - i - 1];}        break;
		case 1: for (int i = 0; i < subdiv; i++) {edge[i] = tile[subdiv - i - 1][subdiv-1];} break;
		case 2: for (int i = 0; i < subdiv; i++) {edge[i] = tile[subdiv-1][i];}              break;
		case 3: for (int i = 0; i < subdiv; i++) {edge[i] = tile[i][0];}                     break;
		}
		return edge;
	}
	
	/**
	 * Check if triangle vertices do not cross mesh split
	 * @param no_connect per-tile 2D array [y][x] corresponding to indices, indicating tile border status:
	 *                   0 - not a border,  1 and 2 - incompatible outer border tiles  
	 * @param x0         x index of the first vertex
	 * @param y0         y index of the first vertex
	 * @param x1         x index of the second vertex
	 * @param y1         y index of the second vertex
	 * @param x2         x index of the third vertex
	 * @param y2         y index of the third vertex
	 * @return           true if such triangle is possible, false - if not
	 */
	public static boolean sameLeafTri(
			int [][] no_connect,
			int x0, int y0,
			int x1, int y1,
			int x2, int y2) {
		return sameLeafTri (new int[] {no_connect[y0][x0], no_connect[y1][x1], no_connect[y2][x2]});
	}
	
	public static boolean sameLeafTri(
			int [][] no_connect,
			int x0, int y0,
			int x1, int y1) {
		return sameLeafTri (new int[] {no_connect[y0][x0],no_connect[y1][x1]});
	}
	
	public static boolean sameLeafTri(
			int [] samples) {
		for (int i = 0; i < samples.length; i ++) {
			if (samples[i] != 0) {
				for (int j = i+1; j < samples.length; j ++) {
					if ((samples[j] != 0) && (samples[j] != samples[i])) {
						return false;
					}
				}
				return true;
			}
		}
		return true;
	}
	
	/**
	 * Triangulate large and small equilateral 45-degree triangles 
	 * @param indices - array of [height][width]{{index}} for large tiles and [heigh][width][py][px]
	 *                  for small ones. This array will be modified and re-indexed if needed.
	 * @return          int [][3] - array of triangles 3 vertex indices, clockwise                 
	 */
	public static int [][] triangulateSameSize(
			int [][][][] indices,
			int [][]     no_connect)
	{
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final int btiles = bwidth * bheight;
		final int [][][][][] tris = new int [bheight][bwidth][][][];

		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger atri = new AtomicInteger(0);
		// initialize triangles array
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				if (indices[btiley][btilex] != null) {
        					tris[btiley][btilex] = new int [indices[btiley][btilex].length][indices[btiley][btilex][0].length][];
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
        
        // triangulate macro (large triangles) - similar as tit was implemented before
        ai.set(0);
        final int btiles_wo_last_row = bwidth * (bheight - 1);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles_wo_last_row; btile = ai.getAndIncrement()) {
        				int x = btile % bwidth;
        				int y = btile / bwidth;
        				if ((indices[y][x] != null) && (indices[y][x].length == 1)){
        					int tris_en = 0;
        					if ((x > 0) &&
        							(indices[y + 1][x - 1] !=null) && (indices[y + 1][x - 1].length == 1) &&
        							(indices[y + 1][x] != null) && (indices[y + 1][x].length == 1)){
        						if (sameLeafTri(no_connect,
        								x,y,
        								x - 1, y + 1,
        								x, y+1)) {
        							tris_en |= (1 << TRI_DOWN_LEFT);
        						}
        					}
        					if (x < (bwidth - 1)) {
        						if ((indices[y + 1][x] != null) && (indices[y + 1][x].length == 1)
        								&& sameLeafTri(no_connect, x, y, x, y + 1)){
        							if ((indices[y][x + 1] != null) && (indices[y][x + 1].length == 1)){
        								if (sameLeafTri(no_connect, x, y, x + 1, y, x, y + 1)) {
        									tris_en |= (1 << TRI_RIGHT_DOWNLEFT);
        								}
        							} else if ((indices[y + 1][x + 1] != null) && (indices[y + 1][x + 1].length ==1)){
        								if (sameLeafTri(no_connect, x, y, x + 1, y + 1, x, y + 1)) {
        									tris_en |= (1 << TRI_DOWNRIGHT_LEFT);
        								}
        							}
        						} else if ((indices[y][x + 1] != null) && (indices[y][x + 1].length == 1) &&
        								(indices[y + 1][x + 1] != null) && (indices[y + 1][x + 1].length ==1)) {
        							if (sameLeafTri(no_connect, x, y, x + 1, y, x + 1, y + 1)) {
        								tris_en |= (1 << TRI_RIGHT_DOWN);
        							}
        						}
        					}
        					if (tris_en != 0) {
        						tris[y][x][0][0]=TRI_NONE.clone();
        						for (int i = 0; i < TRI_NONE.length; i++) {
        							if ((tris_en & (1 << i)) != 0) {
        								tris[y][x][0][0][i] = atri.getAndIncrement();
        							}
        						}
        					}
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);

        // triangulate micro (small triangles) - without crossing tiles borders
        ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				if ((indices[btiley][btilex] != null) && (indices[btiley][btilex].length > 1)){ // subdivided
        					int [][] tindices = indices[btiley][btilex];
        					for (int y = 0;  y < (tindices.length - 1); y++){
        						for (int x = 0; x < tindices[y].length; x++){
        							if (tindices[y][x] >= 0){
        	        					int tris_en = 0;
        								if ((x > 0) && (tindices[y + 1][x - 1] >= 0) && (tindices[y + 1][x] >= 0)){
        									tris_en |= (1 << TRI_DOWN_LEFT);
        								}
        								if (x < (tindices[y].length - 1)) {
        									if (tindices[y + 1][x] >= 0){
        										if (tindices[y][x + 1] >= 0){
        	                						tris_en |= (1 << TRI_RIGHT_DOWNLEFT);
        										} else if (tindices[y + 1][x + 1] >= 0){
        	                						tris_en |= (1 << TRI_DOWNRIGHT_LEFT);
        										}
        									} else if ((tindices[y][x + 1] >= 0) && (tindices[y + 1][x + 1] >= 0)) {
        	            						tris_en |= (1 << TRI_RIGHT_DOWN);
        									}
        								}
        	        					if (tris_en != 0) {
        	        						tris[btiley][btilex][y][x]=TRI_NONE.clone();
        	        						for (int i = 0; i < TRI_NONE.length; i++) {
        	        							if ((tris_en & (1 << i)) != 0) {
        	        								tris[btiley][btilex][y][x][i] = atri.getAndIncrement();
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

        // now triangulate between subdivided tiles - first right and down
        ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			boolean [] quad_corners =    new boolean [TRI_NONE.length]; // 4 corners: top-left, top=right, bottom-right and borrom-left
        			int []     quad_no_connect = new int     [TRI_NONE.length]; // 4 corners: top-left, top=right, bottom-right and borrom-left
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				if ((indices[btiley][btilex] != null) && (indices[btiley][btilex].length > 1)){ // subdivided
        					int [][][] tneib_indices = new int [TileNeibs.DIR_XY.length + 1][][];
        					int [] no_connect_local = new int [TileNeibs.DIR_XY.length + 1];
        					tneib_indices[8] = indices[btiley][btilex];
        					no_connect_local[8] = no_connect[btiley][btilex];
        					for (int dir = 2; dir <7; dir++) { // not all directions needed
        						int btx = btilex + TileNeibs.DIR_XY[dir][0];
        						int bty = btiley + TileNeibs.DIR_XY[dir][1];
        						if ((btx >= 0) && (btx < bwidth) && (bty >= 0) && (bty < bheight) &&
        								(indices[bty][btx] != null) && (indices[bty][btx].length > 1)) {
        							tneib_indices[dir] =    indices[bty][btx];
                					no_connect_local[dir] = no_connect[bty][btx];
        						}
        					}
        					
        					int subdiv = tneib_indices[8].length; // assuming square
        					// First trying TRI_DOWN_LEFT going along left and bottom edge. This point is top-right, not top-left as for others
        					for (int i = 0; i < (2 * subdiv -1); i++) {
        						int x = i;
        						int y = subdiv - 1;
        						if (i >= subdiv) {
        							x = 0;
        							y = i - subdiv;
        						}
        						int num_corn = 0;
        						for (int dir = 1; dir < 4; dir++) { //  skipping top-left corner
        							int x1 = x + TRI_OFFS_XY[dir][0] - 1; // -1 as this point is top-right, not top-left
        							int y1 = y + TRI_OFFS_XY[dir][1];
        							int [] xyd = getNeibNode(x1, y1, subdiv);
        							boolean exists = false;
        							if (tneib_indices[xyd[2]] != null) {
        								exists = tneib_indices[xyd[2]][xyd[1]][xyd[0]] >= 0; // is populated
        								quad_corners[dir] =    exists;
        			        			quad_no_connect[dir] = no_connect_local[xyd[2]];
        							}
        							if (exists) {
        								num_corn ++;
        							} else {
        								break; // here all 3 tested corners should be good
        							}
        						}
        						if (num_corn >= 3) { // all 3 corners exist
        							if (sameLeafTri(new int [] {quad_no_connect[1], quad_no_connect[2], quad_no_connect[3]})){
        								if (tris[btiley][btilex][y][x] == null) {
        									tris[btiley][btilex][y][x] = TRI_NONE.clone();
        								}
        								tris[btiley][btilex][y][x][TRI_DOWN_LEFT] = atri.getAndIncrement();
        							}
        						}
        					}
        					Arrays.fill(quad_corners, false);
        					// now try 3 other triangles
        					for (int i = 0; i < (2 * subdiv -1); i++) {
        						int x = i;
        						int y = subdiv - 1;
        						if (i >= subdiv) {
        							x = subdiv - 1;
        							y =  i - subdiv;
        						}
        						quad_corners[0] = tneib_indices[8][y][x] >= 0; // this subtile
			        			quad_no_connect[0] = no_connect_local[8];
        						if (quad_corners[0]) { // all following triangles assume that top-left corner exists 
        							int num_corn = 1;
        							for (int dir = 1; dir < 4; dir++) { //  skipping top-left corner
        								int x1 = x + TRI_OFFS_XY[dir][0];
        								int y1 = y + TRI_OFFS_XY[dir][1];
        								int [] xyd = getNeibNode(x1, y1, subdiv);
        								boolean exists = false;
        								exists = (tneib_indices[xyd[2]] != null) && (tneib_indices[xyd[2]][xyd[1]][xyd[0]] >= 0); // is populated
        								quad_corners[dir] = exists;
        			        			quad_no_connect[dir] = no_connect_local[xyd[2]];
        								if (exists) {
        									num_corn ++;
        								}
        							}
        							if (num_corn >= 3) {
        								boolean was_null = false; // to remove after
        								int ntri = 0;
        								if (tris[btiley][btilex][y][x] == null) {
        									tris[btiley][btilex][y][x] = TRI_NONE.clone();
        									was_null = true;
        								}
        								if (quad_corners[3]) {
        									if (quad_corners[1]) {
        	        							if (sameLeafTri(new int [] 
        	        									{quad_no_connect[0], quad_no_connect[1], quad_no_connect[3]})){
        	        								tris[btiley][btilex][y][x][TRI_RIGHT_DOWNLEFT] = atri.getAndIncrement();
        	        								ntri++;
        	        							}
        									} else {
        										if (sameLeafTri(new int [] 
        												{quad_no_connect[0], quad_no_connect[2], quad_no_connect[3]})){
        											tris[btiley][btilex][y][x][TRI_DOWNRIGHT_LEFT] = atri.getAndIncrement();
        											ntri++;
        										}
        									}
        								} else {
        									if (sameLeafTri(new int [] 
        											{quad_no_connect[0], quad_no_connect[1], quad_no_connect[2]})){
        										tris[btiley][btilex][y][x][TRI_RIGHT_DOWN] = atri.getAndIncrement();
        										ntri++;
        									}
        								}
        								if ((ntri == 0) && was_null) { // only remove if it was null before
        									tris[btiley][btilex][y][x] = null;
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
        int num_tris =  atri.get();

        // initialize triangles array
        final int [][] tri_indices = new int [num_tris][3];
        // Collect equilateral large and small triangles, generate triangles array
		ai.set(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				int [][][] tile_tris = tris[btiley][btilex];
        				int subdiv = 1;
        				if (tile_tris != null) {
        					int [][][] tneib_indices = null;
        					subdiv = tile_tris.length; 
        					if (subdiv > 1) {
            					tneib_indices = new int [TileNeibs.DIR_XY.length + 1][][];
            					tneib_indices[8] = indices[btiley][btilex];
            					for (int dir = 2; dir <7; dir++) { // not all directions needed
            						int btx = btilex + TileNeibs.DIR_XY[dir][0];
            						int bty = btiley + TileNeibs.DIR_XY[dir][1];
            						if ((btx >= 0) && (btx < bwidth) && (bty >= 0) && (bty < bheight) &&
            								(indices[bty][btx] != null) && (indices[bty][btx].length > 1)) {
            							tneib_indices[dir] = indices[bty][btx];
            						}
            					}
        					}
        					
        					for (int y = 0; y < tile_tris.length; y++) {
            					for (int x = 0; x < tile_tris[y].length; x++) {
            						if (tile_tris[y][x] != null) {
            							for (int nt = 0; nt < tile_tris[y][x].length; nt++) {
            								int tri_indx = tile_tris[y][x][nt];
            								if (tri_indx >= 0) {
            									if (tile_tris.length > 1) { // mini
            										for (int v = 0; v < 3; v++) {
            											int x1 = x + TRI_SET_xy[nt][v][0];
            											int y1 = y + TRI_SET_xy[nt][v][1];
            											int [] xyd = getNeibNode(x1, y1, subdiv);
            											tri_indices[tri_indx][v] =
            													tneib_indices[xyd[2]][xyd[1]][xyd[0]];
            										}            										
            									} else { // macro
            										for (int v = 0; v < 3; v++) {
            											tri_indices[tri_indx][v] =
            													indices [btiley + TRI_SET_xy[nt][v][1]]
            															[btilex + TRI_SET_xy[nt][v][0]][0][0];
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
		return tri_indices;
	}
	
	
	
	/**
	 * Triangulate connections between large (tile centers) and small triangles (subdivided tiles).
	 * All subdivided tiles have full row/column facing non-divided tiles. 
	 * @param indices - array of [height][width]{{index}} for large tiles and [heigh][width][py][px]
	 *                  for small ones. This array will be modified and re-indexed if needed.
	 * @return          int [][3] - array of triangles 3 vertex indices, clockwise                 
	 */
	public static int [][] connectLargeSmallTriangles(
			int [][][][] indices,
			int [][]     no_connect)
	{
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final int btiles = bwidth * bheight;
		final int [][][][] tris = new int [bheight][bwidth][][];
		// subdivided tiles that have non-subdivided neighbors have or may have different triangle configurations:
		// Below subdivided is 2, full tile - 1, empty 0
		// type 0 (4 dirs, subdiv-1 triangles):
		// 2->1 in any of 4 directions produce subdiv - 1 triangles (subdiv-1) connected to tile center
		// type 1 (4 dirs, 1 triangle):
		// two 1 in ortho directions with 0 or 1 between them - produces 1 triangle between
		// centers of 1  and the corner of center 2.
		// type 2 (4 dirs, 1 triangle):
		// 1 in ortho, 2 CCW from it, 0 or 2 (not 1) CCW from 2
		// type 3 (4 dirs, 1 triangle) - mirror of type2:
		// 1 in ortho, 2 CW from it, 0 or 2 (not 1) CW from 2
		// type 4 (4 dirs, 2 triangles):
		// 1 in ortho, 1 CCW from it, 2 CCW from 1. Does not need mirror, as the mirror will be if looking
		// from the last 2 (90 degrees CCW from the first 1)
		// type 5 (4 dirs, 1 triangle): 
		// 1 in ortho, 2 CCW from it, 1 CCW from 1. Does not need mirror
		// type 6 (4 dirs, 1 triangle): 
		// 1 in ortho, 1 CCW from it, 0 CCW from 1.
		// type 7 (4 dirs, 1 triangle) (mirror of 6) 
		// 1 in ortho, 1 CW from it, 0 CW from 1.
	
		// First pass - reserving triangles indices
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger atri = new AtomicInteger(0);
//		final int dbg_x =  32;
//		final int dbg_y =  22;
		// initialize triangles array
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
//        				if ((btiley==dbg_y) && (btilex==dbg_x)) {
//        					System.out.println("connectLargeSmallTriangles().1: btilex="+btilex+", btiley="+btiley);
//       				}
        				if ((indices[btiley][btilex] != null) && (indices[btiley][btilex].length >1)) { // only for subdivided
        					int [] tneib_types = new int [TileNeibs.DIR_XY.length + 1];
        					tneib_types[8] = 2;
        					int [] no_connect_local = new int [TileNeibs.DIR_XY.length + 1];
        					no_connect_local[8] =     no_connect[btiley][btilex];
        					int subdiv = indices[btiley][btilex].length;
        					boolean has_full_neib = false;
        					for (int dir = 0; dir < 8; dir++) { // not all directions needed
        						int btx = btilex + TileNeibs.DIR_XY[dir][0];
        						int bty = btiley + TileNeibs.DIR_XY[dir][1];
        						if ((btx >= 0) && (btx < bwidth) && (bty >= 0) && (bty < bheight) &&
        								(indices[bty][btx] != null)) {
        							tneib_types[dir] =      (indices[bty][btx].length > 1) ? 2 : 1;
        			                no_connect_local[dir] = no_connect[bty][btx];
        							// all modes 0..4 require 2-1 in ortho direction
        							has_full_neib |= (tneib_types[dir] == 1) && ((dir & 1) == 0);
        						}
        					}
        					if (has_full_neib) {
        						tris[btiley][btilex] = new int [8][4]; // [types][directions
        						for (int i = 0; i < tris[btiley][btilex].length; i++) {
        							Arrays.fill(tris[btiley][btilex][i], -1);
        						}
        						// reserve indices for type0:
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];
        							if ((tneib == 1) && sameLeafTri(new int []
        									{no_connect_local[8], no_connect_local[2*dir]})) {
        								tris[btiley][btilex][0][dir] = atri.getAndAdd(subdiv - 1);
        							}
        						}
        						// reserve indices for type1:
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+7) % 8]; // CCW 1 from pointed
        							int tneib2= tneib_types[(2*dir+6) % 8]; // CCW 2 from pointed
        							if ((tneib == 1) && (tneib1 != 2) && (tneib2 == 1) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+6) % 8]})) {
        								tris[btiley][btilex][1][dir] = atri.getAndAdd(1);
        							}
        						}
        						// reserve indices for type2:
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+7) % 8]; // CCW 1 from pointed
        							int tneib2= tneib_types[(2*dir+6) % 8]; // CCW 2 from pointed
        							if ((tneib == 1) && (tneib1 == 2) && (tneib2 != 1) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+7) % 8]})) {
        								tris[btiley][btilex][2][dir] = atri.getAndAdd(1);
        							}
        						}
        						// reserve indices for type3:
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+1) % 8]; // CW 1 from pointed
        							int tneib2= tneib_types[(2*dir+2) % 8]; // CW 2 from pointed
        							if ((tneib == 1) && (tneib1 == 2) && (tneib2 != 1) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+1) % 8]})) {
        								tris[btiley][btilex][3][dir] = atri.getAndAdd(1);
        							}
        						}
        						// reserve indices for type4:
        						// simplified for discontinuities/borders: turn on/off both triangles
        						// at once.
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+7) % 8]; // CCW 1 from pointed
        							int tneib2= tneib_types[(2*dir+6) % 8]; // CCW 2 from pointed
        							if ((tneib == 1) && (tneib1 == 1) && (tneib2 == 2) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+7) % 8],
        											no_connect_local[(2*dir+6) % 8]})) {
        								tris[btiley][btilex][4][dir] = atri.getAndAdd(2);
        							}
        						}
        						// reserve indices for type5:
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+7) % 8]; // CCW 1 from pointed
        							int tneib2= tneib_types[(2*dir+6) % 8]; // CCW 2 from pointed
        							if ((tneib == 1) && (tneib1 == 2) && (tneib2 == 1) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+7) % 8]})) {
        								tris[btiley][btilex][5][dir] = atri.getAndAdd(1);
        							}
        						}
        						// reserve indices for type6:
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+7) % 8]; // CCW 1 from pointed
        							int tneib2= tneib_types[(2*dir+6) % 8]; // CCW 2 from pointed
        							if ((tneib == 1) && (tneib1 == 1) && (tneib2 == 0) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+7) % 8]})) {
        								tris[btiley][btilex][6][dir] = atri.getAndAdd(1);
        							}
        						}
        						// reserve indices for type7: (mirror of 6)
        						for (int dir = 0; dir < 4; dir++) {
        							int tneib = tneib_types[2*dir];         // pointed
        							int tneib1= tneib_types[(2*dir+1) % 8]; // CW 1 from pointed
        							int tneib2= tneib_types[(2*dir+2) % 8]; // CW 2 from pointed
        							if ((tneib == 1) && (tneib1 == 1) && (tneib2 == 0) && 
        									sameLeafTri(new int [] {
        											no_connect_local[8],
        											no_connect_local[2*dir],
        											no_connect_local[(2*dir+1) % 8]})) {
        								tris[btiley][btilex][7][dir] = atri.getAndAdd(1);
        							}
        						}
        						
        					}
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
        int num_tris =  atri.get();

        // initialize triangles array
        final int [][] tri_indices = new int [num_tris][3];
        // Collect seam triangles, generate triangles array
		ai.set(0);
		
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
//        				if ((btiley==dbg_y) && (btilex==dbg_x)) {
//        					System.out.println("connectLargeSmallTriangles().2: btilex="+btilex+", btiley="+btiley);
//        				}
        				if (tris[btiley][btilex] != null) { // only for subdivided
        					int [][][] tneib_indices = new int [TileNeibs.DIR_XY.length + 1][][];
        					tneib_indices [8] = indices[btiley][btilex];
        					int subdiv_m1 = indices[btiley][btilex].length - 1;
        					for (int dir = 0; dir < 8; dir++) { // not all directions needed
        						int btx = btilex + TileNeibs.DIR_XY[dir][0];
        						int bty = btiley + TileNeibs.DIR_XY[dir][1];
        						if ((btx >= 0) && (btx < bwidth) && (bty >= 0) && (bty < bheight) &&
        								(indices[bty][btx] != null)) {
        							tneib_indices[dir] = indices[bty][btx]; 
        						}
        					}
    						// build triangles for type0:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][0][dir4]; // type0
        						if (tri_index >= 0) {
        							int [] edge = getEdgeIndices(tneib_indices [8], dir4);
        							int indx_1 = tneib_indices[2 * dir4][0][0]; // null pointer
        							for (int i = 0; i < subdiv_m1; i++) {
        								tri_indices[tri_index + i][0] = indx_1;
        								tri_indices[tri_index + i][1] = edge[i];
        								tri_indices[tri_index + i][2] = edge[i + 1];
        							}
        						}
        					}
    						// build triangles for type1:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][1][dir4]; // type1
        						if (tri_index >= 0) {
        							int [] edge = getEdgeIndices(tneib_indices [8], dir4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							int indx_2 = tneib_indices[(2 * dir4 + 6) % 8][0][0];
        							tri_indices[tri_index][0] = indx_1;
        							tri_indices[tri_index][1] = edge[subdiv_m1];
        							tri_indices[tri_index][2] = indx_2;
        						}
        					}
    						// build triangles for type2:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][2][dir4]; // type2
        						if (tri_index >= 0) {
        							int [] edge =  getEdgeIndices(tneib_indices [8], dir4);
        							int [] edge1 = getEdgeIndices(tneib_indices [(2 * dir4 + 7) % 8], (dir4 + 1) % 4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							tri_indices[tri_index][0] = indx_1;
        							tri_indices[tri_index][1] = edge[subdiv_m1];
        							tri_indices[tri_index][2] = edge1[0];
        						}
        					}
    						// build triangles for type3:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][3][dir4]; // type3
        						if (tri_index >= 0) {
        							int [] edge =  getEdgeIndices(tneib_indices [8], dir4);
        							int [] edge1 = getEdgeIndices(tneib_indices [(2 * dir4 + 1) % 8], (dir4 + 3) % 4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							tri_indices[tri_index][0] = indx_1;
        							tri_indices[tri_index][1] = edge1[subdiv_m1];
        							tri_indices[tri_index][2] = edge[0];
        						}
        					}
    						// build triangles for type4:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][4][dir4]; // type4
        						if (tri_index >= 0) {
        							int [] edge =  getEdgeIndices(tneib_indices [8], dir4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							int indx_2 = tneib_indices[(2 * dir4 + 7) % 8][0][0];
        							int [] edge1 = getEdgeIndices(tneib_indices [(2 * dir4 + 6) % 8], dir4);
        							
        							tri_indices[tri_index][0] = indx_1;
        							tri_indices[tri_index][1] = edge[subdiv_m1];
        							tri_indices[tri_index][2] = indx_2;
        							
        							tri_indices[tri_index + 1][0] = indx_2;
        							tri_indices[tri_index + 1][1] = edge[subdiv_m1];
        							tri_indices[tri_index + 1][2] = edge1[0];
        						}
        					}
    						// build triangles for type5:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][5][dir4]; // type5
        						if (tri_index >= 0) {
        							int [] edge =  getEdgeIndices(tneib_indices [8], dir4);
        							int [] edge1 = getEdgeIndices(tneib_indices [(2 * dir4 + 7) % 8], (dir4 + 2) % 4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							tri_indices[tri_index][0] = edge [subdiv_m1];
        							tri_indices[tri_index][1] = edge1[subdiv_m1];
        							tri_indices[tri_index][2] = indx_1;
        						}
        					}
    						// build triangles for type6:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][6][dir4]; // type6
        						if (tri_index >= 0) {
        							int [] edge =  getEdgeIndices(tneib_indices [8], dir4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							int indx_2 = tneib_indices[(2 * dir4 + 7) % 8][0][0];
        							tri_indices[tri_index][0] = indx_1;
        							tri_indices[tri_index][1] = edge[subdiv_m1];
        							tri_indices[tri_index][2] = indx_2;
        						}
        					}        					
    						// build triangles for type7:
        					for (int dir4 = 0; dir4 < 4; dir4++) {
        						int tri_index = tris[btiley][btilex][7][dir4]; // type7
        						if (tri_index >= 0) {
        							int [] edge =  getEdgeIndices(tneib_indices [8], dir4);
        							int indx_1 = tneib_indices[2 * dir4][0][0];
        							int indx_2 = tneib_indices[(2 * dir4 + 1) % 8][0][0];
        							tri_indices[tri_index][0] = indx_2;
        							tri_indices[tri_index][1] = edge[0];
        							tri_indices[tri_index][2] = indx_1;
        						}
        					}        					
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
		return tri_indices;
	}
	
	  /**
	   * Triangulate all vertice indices - combine triangulation of same-size equailateral 45-degree
	   * large (tile size) and small (tile subdivisions) and add connections between large and small ones 
	   * @param indices - array of [height][width]{{index}} for large tiles and [heigh][width][py][px]
	   *                  for small ones. This array will be modified and re-indexed if needed.
	   * @return          int [][3] - array of triangles 3 vertex indices, clockwise                 
	   */
	  public static int [][] triangulateAll(
			  int [][][][] indices,
			  int [][]     no_connect)
	  {
		  int [][] tri_same =  triangulateSameSize        (indices, no_connect);
		  int [][] tri_inter = connectLargeSmallTriangles (indices, no_connect);
		  int [][] triangles = new int [tri_same.length + tri_inter.length][];
		  System.arraycopy(tri_same,  0, triangles, 0,               tri_same.length);
		  System.arraycopy(tri_inter, 0, triangles, tri_same.length, tri_inter.length);
		  return triangles;
	  }
	
	
	/**
	 * Get texture coordinates (0..1) for horizontal (positive - to the right) and vertical (positive - up)
	 * @param wh      null if texture image is cropped one, or {full_width, full_height} otherwise 
	 * @param indices pairs of [x,y] in integer tiles
	 * @return texture x,y for the tile centers
	 */
	
	public static double [][] getTexCoords( // get texture coordinates for indices
			int []   wh,     // 0 or full width, full height of the image in tiles
			int [][] indices)
	{
		int maxIndex = -1;
		int height = indices.length;
		int width = indices[0].length;
		int tex_width =  (wh != null)? wh[0]: width;
		int tex_height = (wh != null)? wh[1]: height;
		outer_label:{
			for (int y = height - 1 ; y >= 0; y--) {
				for (int x = width - 1; x >= 0; x--){
					if (indices[y][x] >=0){
						maxIndex = indices[y][x];
						break outer_label;
					}
				}
			}
		}
		double [][] textureCoordinate = new double [maxIndex+1][2];
		int indx = 0;
		for (int y = 0;  indx <= maxIndex; y++) {
			for (int x = 0; (x < width) && (indx <= maxIndex); x++){
				if (indices[y][x] >=0){
					textureCoordinate[indx][0] = (x + 0.5)/tex_width;
					textureCoordinate[indx][1] = (tex_height - y - 0.5) / tex_height; // y is up
					indx ++;
				}
			}
		}
		return textureCoordinate;
	}
	
	/**
	 * Create texture coordinates for 2 levels of triangles - tile centers and subdivided tiles 
	 * @param rect        null if texture image matches indices[][][][] array, otherwise
	 *                    x,y - top left corner corresponding indices, width, height - full image
	 *                    size in tiles
	 * @param num_indices total number of vertices in indices[][][][] array
	 * @param indices     two level vertices indices - [tilesy][tilesx][1][1] for "large" tiles
	 *                    (single vertice in the tile centre) or [tilesy][tilesx][subdiv][subdiv]
	 *                    for subdivided indices[tilesy][tilesx] may be null if the tile is empty,
	 *                    [tilesy][tilesx][y][x] is either unique index or -1 for missing subtile.
	 * @return   [num_indices][2] texture coordinates in 0..1.0 range, positive texture y is up
	 */
	public static double [][] getTexCoords( // get texture coordinates for indices
			final Rectangle    rect,     // 0 or full width, full height of the image in tiles
			final int          num_indices,
			final int [][][][] indices)
	{
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final int btiles = bwidth * bheight;
		final int tex_width =  (rect != null)? rect.width: bwidth;
		final int tex_height = (rect != null)? rect.height: bheight;
		final int tex_x =      (rect != null)? rect.x: 0;
		final int tex_y =      (rect != null)? rect.y: 0;
		final double [][] textureCoordinate = new double [num_indices][2];
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				if (indices[btiley][btilex] != null) {
        					int subdiv = indices[btiley][btilex].length;
        					for (int y = 0; y < indices[btiley][btilex].length; y++) {
            					for (int x = 0; x < indices[btiley][btilex][y].length; x++) {
            						int indx = indices[btiley][btilex][y][x];
            						if (indx >= 0) {
            							textureCoordinate[indx][0] = (tex_x + btilex + (x + 0.5) / subdiv) / tex_width;
            							textureCoordinate[indx][1] = (tex_height - tex_y - btiley - (y + 0.5) / subdiv) / tex_height;
            						}
            					}
        					}
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
		return textureCoordinate;
	}
	
	public static double [] getIndexedDisparities( // get disparity for each index
			final double []     disparity,
			final double        min_disparity,
			final double        max_disparity,
			final Rectangle     bounds,
			final int           num_indices,
			final int [][][][]  indices,
			final int           tilesX)
	{
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final int btiles = bwidth * bheight;
		final boolean full_disparity = disparity.length > (bounds.height * bounds.width);
		final double [] indexedDisparity = new double [num_indices];

		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				int tile = btile;
        				if (full_disparity) {
        					tile = (btiley + bounds.y) * tilesX + (btilex + bounds.x);  
        				}
        				if (indices[btiley][btilex] != null) {
        					int subdiv = indices[btiley][btilex].length;
        					double disp = disparity[tile];
    						if      (disp < min_disparity) disp = min_disparity;
    						else if (disp > max_disparity) disp = max_disparity;
        					for (int y = 0; y < indices[btiley][btilex].length; y++) {
            					for (int x = 0; x < indices[btiley][btilex][y].length; x++) {
            						int indx = indices[btiley][btilex][y][x];
            						if (indx >= 0) {
            							indexedDisparity[indx] = disp; // currently - same disparity for subdivisions
            							// TODO - linear interpolate
            						}
            					}
        					}
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
		return indexedDisparity;
	}
	
	public static double [][] getCoords( // get disparity for each index
			final double []          disparity,
			final double             min_disparity,
			final double             max_disparity,
			final Rectangle          bounds,
			final int                num_indices,
			final int [][][][]       indices,
			final int                tilesX,
			final int                tile_size,
			final boolean            correctDistortions, // requires backdrop image to be corrected also
			final GeometryCorrection geometryCorrection)
			
	{
		final int bwidth=indices[0].length;
		final int bheight=indices.length;
		final int btiles = bwidth * bheight;
		final boolean full_disparity = disparity.length > (bounds.height * bounds.width);
		final double [][] coordinate = new double [num_indices][];

		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
        				int btilex = btile % bwidth;
        				int btiley = btile / bwidth;
        				int tile = btile;
        				if (full_disparity) {
        					tile = (btiley + bounds.y) * tilesX + (btilex + bounds.x);  
        				}
        				if (indices[btiley][btilex] != null) {
        					int subdiv = indices[btiley][btilex].length;
        					double disp = disparity[tile];
    						if      (disp < min_disparity) disp = min_disparity;
    						else if (disp > max_disparity) disp = max_disparity;
        					for (int y = 0; y < indices[btiley][btilex].length; y++) {
            					for (int x = 0; x < indices[btiley][btilex][y].length; x++) {
            						// TODO - linear interpolate disparity here
            						int indx = indices[btiley][btilex][y][x];
            						if (indx >= 0) {
                						double px = (bounds.x + btilex + (x + 0.5) / subdiv) * tile_size - 0.5;
                						double py = (bounds.y + btiley + (y + 0.5) / subdiv) * tile_size - 0.5;
            							coordinate[indx] = geometryCorrection.getWorldCoordinates(
            									px,
            									py,
            									disp, // currently - same for all tile
            									correctDistortions);
            						}
            					}
        					}
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
		return coordinate;
	}
	
	
	
	/** Plot generated triangles on a provided canvas. Can combine multiple meshes on the same canvas
	 * 
	 * @param canvas         2D canvas in a line-scan order
	 * @param width          canvas width in pixels
	 * @param tex_coord      texture coordinates as [vertice_index]{x,y}, where 0 <=x <1.0, 0 <=y <1.0, y - upward 
	 * @param triangles      [tri_index]{p0_index, p1_index, p2_index}, where p*_index is vertice_index in tex_coord 
	 * @param plot_center    plot triangle center
	 * @param line_color     pixel to use for lines
	 * @param center_color   pixel to use for triangle centers
	 */
	public static void plotMesh(
			final double []   canvas,
			final int         width,
			final double [][] tex_coord,
			final int [][]    triangles,
			final boolean     plot_center,
			final double      line_color,
			final double      center_color)
	{
		final int height = canvas.length / width;
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
        	threads[ithread] = new Thread() {
        		public void run() {
        			int width_m1 = width-1;
        			int height_m1 = height-1;
        			for (int ntri  = ai.getAndIncrement(); ntri < triangles.length; ntri = ai.getAndIncrement()) {
        				for (int i0 = 0; i0 < 3; i0++) {
        					int i1 = (i0+1) % 3;
        					double [] pxy0 = {
        							width*tex_coord[triangles[ntri][i0]][0],
        							height*(1.0 - tex_coord[triangles[ntri][i0]][1])};
        					double [] pxy1 = {
        							width*tex_coord[triangles[ntri][i1]][0],
        							height*(1.0 - tex_coord[triangles[ntri][i1]][1])};
        					double dx = pxy1[0] - pxy0[0];
        					double dy = pxy1[1] - pxy0[1];
        					double l = Math.sqrt(dx*dx+dy*dy);
        					for (int j = 0; j < l; j++) {
        						int px = (int) Math.round(pxy0[0]+j*dx/l);
        						int py = (int) Math.round(pxy0[1]+j*dy/l);
        						px = Math.min(Math.max(0, px), width_m1);
        						py = Math.min(Math.max(0, py), height_m1);
        						canvas[py*width+px] = line_color;
        					}
        				}
        				if (plot_center) {
        					double sx = 0.0, sy = 0.0;
            				for (int i = 0; i < 3; i++) {
            					sx+= width*tex_coord[triangles[ntri][i]][0];
            					sy+= height*(1.0 - tex_coord[triangles[ntri][i]][1]);
            				}
            				sx /= 3;
            				sy /= 3;
            				int px = (int) Math.round(sx);
            				int py = (int) Math.round(sy);
            				px = Math.min(Math.max(0, px), width_m1);
            				py = Math.min(Math.max(0, py), height_m1);
            				canvas[py*width+px] = center_color;
        				}
        			}                	
        		}
        	};
        }		      
        ImageDtt.startAndJoin(threads);
	}

	final static double [][] plotForMesh(
			final boolean    full_selection, // false
			final boolean    full_alpha,     // true
			final Rectangle  bounds,
			final int        out_width,
			final int        tilesX,
			final int        tilesY,
			final int        tile_size,		
			final int        subdiv,
			final double []  disparity,
			final boolean [] selection,
			final int []     border_int,
			final boolean [] alpha,
			final double []  dalpha){
		// scale pixels to match mesh triangles
		final double scale_disparity = 0.05;
		final double scale_selection = 1.0;
		final double scale_border =    0.3;
		final double scale_alpha =     1.0;
		final double scale_dalpha =    1.0;
		final int out_height = out_width * tilesY / tilesX;
		final int out_scale = out_width / (tilesX * tile_size);
		final int out_tile = out_scale * tile_size; 
		final double [][] rslt_img = new double [5][];
		if (disparity != null)  rslt_img[0] = new double [out_height * out_width];
		if (selection != null)  rslt_img[1] = new double [out_height * out_width];
		if (border_int != null) rslt_img[2] = new double [out_height * out_width];
		if (alpha != null)      rslt_img[3] = new double [out_height * out_width];
		if (dalpha != null)     rslt_img[4] = new double [out_height * out_width];
		final int source_tile_width =  full_selection? tilesX : bounds.width;
		final int source_tile_offsx =  full_selection? bounds.x : 0;
		final int source_tile_offsy =  full_selection? bounds.y : 0;
		final int source_pix_offsx = (full_alpha? bounds.x : 0) * tile_size;
		final int source_pix_offsy = (full_alpha? bounds.y : 0) * tile_size;
		final int source_pix_width = (full_alpha? tilesX : bounds.width) * tile_size;
		final int btiles = bounds.width * bounds.height;
		final Thread[] threads = ImageDtt.newThreadArray(TexturedModel.THREADS_MAX);
		final AtomicInteger ai = new AtomicInteger(0);
		// mark disparity, selection, borders
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                public void run() {
					for (int btile = ai.getAndIncrement(); btile < btiles; btile = ai.getAndIncrement()) {
						int btilex = btile % bounds.width;
						int btiley = btile / bounds.width;
						int stile = (btilex + source_tile_offsx) + (btiley + source_tile_offsy) * source_tile_width;
						int ftilex = bounds.x + btilex; // relative to full tile frame tilesX*tilesY
						int ftiley = bounds.y + btiley;
						// mark disparity, selection, borders
						// out_scale * tile_size
						int out_indx0 = (ftilex * out_tile) + (ftiley * out_tile) * out_width;
						double d_disp = (disparity != null) ?                      scale_disparity * disparity[stile] : Double.NaN;
						double d_sel =  ((selection != null) && selection[stile])? scale_selection : Double.NaN;
						double d_bord = ((border_int != null) && (border_int[stile] >= 0))? scale_border * border_int[stile] : Double.NaN;
						for (int oy = 0; oy < out_tile; oy++) {
							int out_indx1 = out_indx0 + oy*out_width;
							if (disparity != null)  Arrays.fill(rslt_img[0], out_indx1, out_indx1 + out_tile, d_disp);
							if (selection != null)  Arrays.fill(rslt_img[1], out_indx1, out_indx1 + out_tile, d_sel);
							if (border_int != null) Arrays.fill(rslt_img[2], out_indx1, out_indx1 + out_tile, d_bord);
						}						
						// index for alpha - top left of as tile
						if (alpha != null) {
							int aindx0 = (source_pix_offsx + btilex * tile_size) + (source_pix_offsy + btiley * tile_size) * source_pix_width;
							for (int ay = 0; ay < tile_size; ay++) {
								int aindx1 = aindx0 + ay*source_pix_width;
								int out_indx1 = out_indx0 + ay * out_scale * out_width;
								for (int ax = 0; ax < tile_size; ax++) {
									int aindx = aindx1 + ax;
									int out_indx2 = out_indx1 + ax * out_scale;
									double d_alpha = alpha[aindx]? scale_alpha : 0.0;
									for (int row = 0; row < out_scale; row++) {
										int out_indx = out_indx2 + row * out_width;
										Arrays.fill(rslt_img[3], out_indx, out_indx + out_scale, d_alpha);
									}
								}
							}
						}
						// index for dalpha - top left of as tile
						if (dalpha != null) {
							int aindx0 = (source_pix_offsx + btilex * tile_size) + (source_pix_offsy + btiley * tile_size) * source_pix_width;
							for (int ay = 0; ay < tile_size; ay++) {
								int aindx1 = aindx0 + ay*source_pix_width;
								int out_indx1 = out_indx0 + ay * out_scale * out_width;
								for (int ax = 0; ax < tile_size; ax++) {
									int aindx = aindx1 + ax;
									int out_indx2 = out_indx1 + ax * out_scale;
									double d_alpha = dalpha[aindx] * scale_dalpha;
									for (int row = 0; row < out_scale; row++) {
										int out_indx = out_indx2 + row * out_width;
										Arrays.fill(rslt_img[4], out_indx, out_indx + out_scale, d_alpha);
									}
								}
							}
						}
					}
                }
            };
        }		      
        ImageDtt.startAndJoin(threads);
		return rslt_img;
	}
	
	
	
	
	
	
	public static double [] getIndexedDisparities( // get disparity for each index
			double [] disparity,
			double min_disparity,
			double max_disparity,
			Rectangle bounds,
			int [][]  indices,
			int       tile_size,
			int       tilesX)
	{
		//		  int height = indices.length;
		int width = indices[0].length;
		int maxIndex = getMaxIndex(indices);
		double [] indexedDisparity = new double [maxIndex+1];
		int indx = 0;
		if (disparity.length > (bounds.height * bounds.width)) { // old version - selected is full size
			for (int y = 0;  indx <= maxIndex; y++) {
				for (int x = 0; (x < width) && (indx <= maxIndex); x++){
					if (indices[y][x] >=0){
						// center coordinates for 8*8 tile is [3.5,3.5]
						double disp = (disparity == null)? min_disparity:( disparity[(bounds.y + y) * tilesX + (bounds.x + x)]);
						if      (disp < min_disparity) disp = min_disparity;
						else if (disp > max_disparity) disp = max_disparity;
						indexedDisparity[indx] =disp;
						indx ++;
					}
				}
			} 
		} else { // 09.18.2022
			for (int y = 0;  indx <= maxIndex; y++) {
				for (int x = 0; (x < width) && (indx <= maxIndex); x++){
					if (indices[y][x] >=0){
						// center coordinates for 8*8 tile is [3.5,3.5]
						double disp = (disparity == null)? min_disparity:( disparity[bounds.width *y + x]);
						if      (disp < min_disparity) disp = min_disparity;
						else if (disp > max_disparity) disp = max_disparity;
						indexedDisparity[indx] =disp;
						indx ++;
					}
				}
			} 
		}
		return indexedDisparity;
	}

	public static double [][] getCoords( // get world XYZ in meters for indices
			double [] disparity, // null - use min_disparity
			double min_disparity,
			double max_disparity,
			Rectangle bounds,
			int [][]  indices,
			int       tile_size,
			int       tilesX,
			boolean   correctDistortions, // requires backdrop image to be corrected also
			GeometryCorrection geometryCorrection)
	{
		//		  int height = indices.length;
		int width = indices[0].length;
		int maxIndex = getMaxIndex(indices);
		double [][] coordinate = new double [maxIndex+1][];
		int indx = 0;
		if (disparity.length > (bounds.height * bounds.width)) { // old version - selected is full size
			for (int y = 0;  indx <= maxIndex; y++) {
				for (int x = 0; (x < width) && (indx <= maxIndex); x++){
					if (indices[y][x] >=0){
						// center coordinates for 8*8 tile is [3.5,3.5]
						double px = (bounds.x + x + 0.5) * tile_size - 0.5;
						double py = (bounds.y + y + 0.5) * tile_size - 0.5;
						double disp = (disparity == null)? min_disparity:( disparity[(bounds.y + y) * tilesX + (bounds.x + x)]);
						if      (disp < min_disparity) disp = min_disparity;
						else if (disp > max_disparity) disp = max_disparity;
						coordinate[indx] = geometryCorrection.getWorldCoordinates(
								px,
								py,
								disp,
								correctDistortions);
						indx ++;
					}
				}
			}
		} else { // 09.18.2022
			for (int y = 0;  indx <= maxIndex; y++) {
				for (int x = 0; (x < width) && (indx <= maxIndex); x++){
					if (indices[y][x] >=0){
						// center coordinates for 8*8 tile is [3.5,3.5]
						double px = (bounds.x + x + 0.5) * tile_size - 0.5;
						double py = (bounds.y + y + 0.5) * tile_size - 0.5;
						double disp = (disparity == null)? min_disparity:( disparity[bounds.width * y + x]);
						if      (disp < min_disparity) disp = min_disparity;
						else if (disp > max_disparity) disp = max_disparity;
						coordinate[indx] = geometryCorrection.getWorldCoordinates(
								px,
								py,
								disp,
								correctDistortions);
						indx ++;
					}
				}
			}

		}
		return coordinate;
	}

	public static int [][] filterTriangles(
			int  [][] triangles,
			double [] disparity, // disparities per vertex index
			double    maxDispDiff, // maximal relative disparity difference in a triangle
			int       debug_level)
	{
		final double min_avg = 3.0; // 0.5; // minimal average disparity to normalize triangle
		class Triangle {
			int [] points = new int [3];
			Triangle (int i1, int i2, int i3){
				points[0] = i1;
				points[1] = i2;
				points[2] = i3;
			}
		}
		ArrayList<Triangle> triList = new ArrayList<Triangle>();
		for (int i = 0; i < triangles.length; i++){
			double disp_avg = (disparity[triangles[i][0]] + disparity[triangles[i][1]]+ disparity[triangles[i][2]])/3.0; // fixed 09.18.2022!
			if (disp_avg < min_avg) disp_avg = min_avg;
			loop:{
				for (int j = 0; j < 3; j++){
					int j1 = (j + 1) % 3;
					if (Math.abs(disparity[triangles[i][j]] - disparity[triangles[i][j1]]) > (disp_avg* maxDispDiff)) {
						if (debug_level > 1) {
							System.out.println("removed triangle "+i+": "+
									disparity[triangles[i][0]]+". "+disparity[triangles[i][1]]+". "+disparity[triangles[i][2]]+
									". Avg = "+disp_avg);
						}
						break loop;
					}
				}
				triList.add(new Triangle(
						triangles[i][0],
						triangles[i][1],
						triangles[i][2]));
			}
		}
		int [][] filteredTriangles = new int [triList.size()][3];
		for (int i = 0; i < filteredTriangles.length; i++){
			filteredTriangles[i] = triList.get(i).points;
		}
		return filteredTriangles;
	}

	public static int [][] filterTrianglesWorld(
			int  [][] triangles,
			double [][] worldXYZ,  // world per vertex index
			double      maxZtoXY,
			double maxZ) 
	{
		final double maxZtoXY2 = maxZtoXY * maxZtoXY;
		class Triangle {
			int [] points = new int [3];
			Triangle (int i1, int i2, int i3){
				points[0] = i1;
				points[1] = i2;
				points[2] = i3;
			}
		}
		ArrayList<Triangle> triList = new ArrayList<Triangle>();
		for (int i = 0; i < triangles.length; i++){
			double [][] min_max = new double[3][2];
			boolean not_too_far = true;
			for (int di = 0; di < 3; di++) {
				min_max[di][0] = worldXYZ[triangles[i][0]][di];
				min_max[di][1] = min_max[di][0]; // both min and max to the same vertex 0
			}
			if (maxZ != 0) {
				not_too_far &=  worldXYZ[triangles[i][0]][2] > -maxZ; 
			}
			for (int vi = 1; vi < 3; vi++) {
				for (int di = 0; di < 3; di++) {
					min_max[di][0] = Math.min(min_max[di][0], worldXYZ[triangles[i][vi]][di]);
					min_max[di][1] = Math.max(min_max[di][1], worldXYZ[triangles[i][vi]][di]);
				}
			}
			double dx = min_max[0][1]-min_max[0][0];
			double dy = min_max[1][1]-min_max[1][0];
			double dz = min_max[2][1]-min_max[2][0];
			double ratio2 = dz*dz/(dx*dx+dy*dy + 0.001);
			if (not_too_far && ((maxZtoXY == 0) || (ratio2 < maxZtoXY2))) {
				triList.add(new Triangle(
						triangles[i][0],
						triangles[i][1],
						triangles[i][2]));
			}
		}
		int [][] filteredTriangles = new int [triList.size()][3];
		for (int i = 0; i < filteredTriangles.length; i++){
			filteredTriangles[i] = triList.get(i).points;
		}
		return filteredTriangles;
	}
	
	public static int [] reIndex(
			int [][]  indices,
			int [][] triangles) {
		int last_index = -1;
		for (int i = 0; i < indices.length; i++) {
			for (int j = 0; j < indices[i].length; j++) {
				if (indices[i][j] > last_index) {
					last_index = indices[i][j]; 
				}
			}
		}
		boolean [] used_indices = new boolean[last_index+1];
		for (int i = 0; i < triangles.length; i++) {
			for (int j = 0; j < triangles[i].length; j++) { // always 3
				used_indices[triangles[i][j]] = true;
			}
		}
		int new_len = 0;
		for (int i = 0; i < used_indices.length; i++) if (used_indices[i]) {
			new_len++;
		}
		if (new_len == used_indices.length) {
			return null; // no re-indexing is needed
		}
		int [] re_index = new int [new_len];
		int indx = 0;
		for (int i = 0; i < indices.length; i++) {
			for (int j = 0; j < indices[i].length; j++) {
				int old_index=indices[i][j];
				if (old_index >= 0) {
					if (used_indices[old_index]) { // keep
						re_index[indx] = old_index;
						indices[i][j] = indx++;
					} else {
						indices[i][j] = -1;
					}
				}
			}
		}
		return re_index;
	}
	
	public static int [][] triangulateCluster(
			int [][]  indices)
	{
		int height = indices.length;
		int width = indices[0].length;
		class Triangle {
			int [] points = new int [3];
			Triangle (int i1, int i2, int i3){
				points[0] = i1;
				points[1] = i2;
				points[2] = i3;
			}
		}
		ArrayList<Triangle> triList = new ArrayList<Triangle>();
		for (int y = 0;  y < (height - 1); y++){
			for (int x = 0; x < width; x++){
				if (indices[y][x] >= 0){
					if ((x > 0) && (indices[y + 1][x - 1] >= 0) && (indices[y + 1][x] >= 0)){
						triList.add(new Triangle(
								indices[y][x],
								indices[y + 1][x],
								indices[y + 1][x - 1]));
					}
					if (x < (width - 1)) {
						if (indices[y + 1][x] >= 0){
							if (indices[y][x + 1] >= 0){
								triList.add(new Triangle(
										indices[y][x],
										indices[y][x + 1],
										indices[y + 1][x]));

							} else if (indices[y + 1][x + 1] >= 0){
								triList.add(new Triangle(
										indices[y][x],
										indices[y + 1][x + 1],
										indices[y + 1][x]));

							}
						} else if ((indices[y][x + 1] >= 0) && (indices[y + 1][x + 1] >= 0)) {
							triList.add(new Triangle(
									indices[y][x],
									indices[y][x + 1],
									indices[y + 1][x + 1]));
						}
					}
				}
			}
		}
		int [][] triangles = new int [triList.size()][3];
		for (int i = 0; i < triangles.length; i++){
			triangles[i] = triList.get(i).points;
		}
		return triangles;
	}
	
	
	
	public static void testTriangles(
			String     texturePath, // if not null - will show
			Rectangle  bounds,
			boolean [] selected,
			double []  disparity,
			int        tile_size,
			int        tilesX,
			int        tilesY,
			int [][]   indices,
			int [][]   triangles,
			double [][] debug_triangles) // if not null - should be [2][width* height], will mark disparity and triangles
	{
		String [] titles = {"disparity","triangles"};
		double [][] dbg_img = new double [titles.length][tilesX*tilesY*tile_size*tile_size];
		Arrays.fill(dbg_img[0], Double.NaN);
		if (selected.length > (bounds.height * bounds.width)) { // old version - selected is full size
			for (int i = 0; i < selected.length; i++ ){
				double d = selected[i]? ((disparity.length >1) ? disparity[i] : disparity[0]):Double.NaN;
				int y = i / tilesX;
				int x = i % tilesX;
				for (int dy = 0; dy <tile_size; dy ++){
					for (int dx = 0; dx <tile_size; dx ++){
						dbg_img[0][(y * tile_size + dy)*(tile_size*tilesX) + (x * tile_size + dx)] = d;
					}
				}
			}
		} else { // 09.18.2022
			for (int i = 0; i < selected.length; i++ ){
				double d = selected[i]? ((disparity.length > 1) ? disparity[i] : disparity[0]):Double.NaN;
				int y = i / bounds.width + bounds.y;
				int x = i % bounds.width + bounds.x;
				for (int dy = 0; dy <tile_size; dy ++){
					for (int dx = 0; dx <tile_size; dx ++){
						dbg_img[0][(y * tile_size + dy)*(tile_size*tilesX) + (x * tile_size + dx)] = d;
					}
				}
			}
		}
		int maxIndex = getMaxIndex(indices);
		int [][] pxy = new int [maxIndex+1][2];
		int height = indices.length;
		int width = indices[0].length;

		for (int y = 0; y < height; y++){
			for (int x = 0; x < width; x++){
				if (indices[y][x] >= 0){
					pxy[indices[y][x]][0] =  (bounds.x + x)*tile_size + (tile_size/2);
					pxy[indices[y][x]][1] =  (bounds.y + y)*tile_size + (tile_size/2);
				}
			}
		}

		for (int i = 0; i < triangles.length; i++ ){
			for (int side = 0; side < triangles[i].length; side++){
				int [] pntIndx = {
						triangles[i][side],
						(side == (triangles[i].length -1)? triangles[i][0]:triangles[i][side+1])};
				int dx = iSign(pxy[pntIndx[1]][0] - pxy[pntIndx[0]][0]);
				int dy = iSign(pxy[pntIndx[1]][1] - pxy[pntIndx[0]][1]);
				for (int j = 0; j < tile_size; j++){
					int x = pxy[pntIndx[0]][0] + dx*j;
					int y = pxy[pntIndx[0]][1] + dy*j;
					dbg_img[1][y * tile_size * tilesX + x] = 10.0; //1711748
				}
			}
		}
		if (texturePath != null) {
			ShowDoubleFloatArrays.showArrays(
					dbg_img,
					tilesX * tile_size,
					tilesY * tile_size,
					true,
					"triangles-"+texturePath,
					titles);
		}
		if (debug_triangles != null) {
			int indx_tri = (debug_triangles.length>1) ? 1 : 0;
			for (int i = 0; i < debug_triangles[indx_tri].length; i++) {
				if (dbg_img[1][i] > 0) {
					debug_triangles[indx_tri][i] = dbg_img[1][i]; // 10.0 to have the same scale as disparity
				}
			}
			if (indx_tri > 0) {
				for (int i = 0; i < debug_triangles[indx_tri].length; i++) {
					if (!Double.isNaN(dbg_img[0][i])) {
						debug_triangles[0][i] = dbg_img[0][i]; // disparity if not NaN
					}
				}
			}
		}
	}
	
	static int iSign (int a) {return (a > 0) ? 1 : ((a < 0)? -1 : 0);}
	static int getMaxIndex(int [][] indices)
	{
		int height = indices.length;
		int width = indices[0].length;
		for (int y = height - 1 ; y >= 0; y--) {
			for (int x = width - 1; x >= 0; x--){
				if (indices[y][x] >= 0){
					return indices[y][x];
				}
			}
		}
		return -1;
	}

	  public static void generateClusterX3d( 
			  boolean         full_texture, // true - full size image, false - bounds only
			  int             subdivide_mesh, // 0,1 - full tiles only, 2 - 2x2 pixels, 4 - 2x2 pixels
//			  boolean []      alpha,     // boolean alpha - true - opaque, false - transparent. Full/bounds
			                             // matching selection
			  X3dOutput       x3dOutput, // output x3d if not null
			  WavefrontExport wfOutput,  // output WSavefront if not null
			  ArrayList<TriMesh> tri_meshes,
			  String          texturePath,
			  String          id,
			  String          class_name,
			  Rectangle       bounds,
			  boolean []      selected, // may be either tilesX * tilesY or bounds.width*bounds.height
			  double []       disparity, // if null, will use min_disparity
			  int             tile_size,
			  int             tilesX,
			  int             tilesY,
			  GeometryCorrection geometryCorrection,
			  boolean         correctDistortions, // requires backdrop image to be corrected also
			  boolean         show_triangles,
			  double          min_disparity,
			  double          max_disparity,
			  double          maxDispTriangle, // relative <=0 - do not use
			  double          maxZtoXY,        // 10.0. <=0 - do not use
			  double          maxZ,            // far clip (0 - do not clip). Negative - limit by max
			  boolean         limitZ,
			  double [][]     dbg_disp_tri_slice,
			  int             debug_level
			  ) throws IOException
	  {
//		  int debug_level = 1;
		  if (bounds == null) {
			  return; // not used in lwir
		  }
		  int [][] indices =  getCoordIndices( // starting with 0, -1 - not selected // updated 09.18.2022
				  bounds, 
				  selected,
				  tilesX); 
		  double [][] texCoord = getTexCoords( // get texture coordinates for indices
				  full_texture ? (new int[] {tilesX, tilesY}): null,
				  indices);

		  double [][] worldXYZ = getCoords( // get world XYZ in meters for indices // updated 09.18.2022
				  disparity,
				  min_disparity,
				  max_disparity,
				  bounds,
				  indices,
				  tile_size,
				  tilesX,
				  correctDistortions, // requires backdrop image to be corrected also
				  geometryCorrection);

          double [] indexedDisparity = getIndexedDisparities( // get disparity for each index // updated 09.18.2022
							disparity,
							min_disparity,
							max_disparity,
							bounds,
							indices,
							tile_size,
							tilesX);

		  int [][] triangles = 	triangulateCluster(
				  indices);

		  int num_removed = 0;
		  if (maxDispTriangle > 0.0) {
			  int pre_num = triangles.length;
			  triangles = 	filterTriangles( // remove crazy triangles with large disparity difference
					  triangles,
					  indexedDisparity, // disparities per vertex index
					  maxDispTriangle,  // maximal disparity difference in a triangle
					  debug_level + 0); // 	int       debug_level);
			  if (triangles.length < pre_num) {
				  num_removed += pre_num - triangles.length;
				  if (debug_level > 0) {
					  System.out.println("filterTriangles() removed "+ (pre_num - triangles.length)+" triangles");
				  }
			  }
			  
		  }
		  if ((maxZ != 0.0) && limitZ) {
			  for (int i = 0; i < worldXYZ.length; i++) {
				  if (worldXYZ[i][2] < -maxZ) {
					  double k = -maxZ/worldXYZ[i][2];
					  worldXYZ[i][0] *= k;
					  worldXYZ[i][1] *= k;
					  worldXYZ[i][2] *= k;
				  }
			  }
		  }
		  if ((maxZtoXY > 0.0) || ((maxZ != 0) && !limitZ) ) {
			  int pre_num = triangles.length;
			  triangles = 	filterTrianglesWorld(
					  triangles,
					  worldXYZ,  // world per vertex index
					  maxZtoXY,
					  maxZ);
			  if (triangles.length < pre_num) {
				  num_removed += pre_num - triangles.length;
				  if (debug_level > 0) {
					  System.out.println("filterTrianglesWorld() removed "+ (pre_num - triangles.length)+" triangles");
				  }
			  }
		  }
		  if (triangles.length == 0) {
			  if (debug_level > 0) {
				  System.out.println("generateClusterX3d() no triangles left in a cluster");
			  }
			  return; // all triangles removed
			  
		  }
		  if (num_removed > 0) { 
			  int [] re_index = reIndex( // Move to TriMesh?
					  indices, // will be modified if needed (if some indices are removed
					  triangles);
			  if (re_index != null) {// need to update other arrays: texCoord, worldXYZ. indexedDisparity[] will not be used
				  int num_indices_old = worldXYZ.length;
				  int [] inv_index = new int [num_indices_old];
				  Arrays.fill(inv_index,-1); // just to get an error
				  for (int i = 0; i < re_index.length; i++) {
					  inv_index[re_index[i]] = i;
				  }
				  double [][]  texCoord_new = new double [re_index.length][];
				  double [][]  worldXYZ_new = new double [re_index.length][];
				  for (int i = 0; i < re_index.length; i++) {
					  texCoord_new[i] = texCoord[re_index[i]];
					  worldXYZ_new[i] = worldXYZ[re_index[i]];
				  }
				  texCoord = texCoord_new;
				  worldXYZ = worldXYZ_new;
				  for (int i = 0; i < triangles.length; i++) {
					  for (int j=0; j < triangles[i].length; j++) {
						  triangles[i][j] = inv_index[triangles[i][j]];
					  }
				  }
			  }
			  if (debug_level > 0) {
				  show_triangles = true; // show after removed
			  }
		  }

		  if (show_triangles || (dbg_disp_tri_slice != null)) {
			  double [] ddisp = (disparity == null)?(new double[1]):disparity;
			  if (disparity == null) {
				  ddisp[0] = min_disparity;
			  }
			  testTriangles(
					  (show_triangles? texturePath: null),
					  bounds,
					  selected,
					  ddisp, // disparity, // if disparity.length == 1 - use for all
					  tile_size,
					  tilesX, // int        tilesX,
					  tilesY, // 						int        tilesY,
					  indices,
					  triangles,
					  dbg_disp_tri_slice); // double [][] debug_triangles);
		  }
		  if (x3dOutput != null) {
		  x3dOutput.addCluster(
				  texturePath,
				  id,
				  class_name,
				  texCoord,
				  worldXYZ,
				  triangles);
		  }
		  if (wfOutput != null) {
			  wfOutput.addCluster(
				  texturePath,
				  id,
//				  class_name,
				  texCoord,
				  worldXYZ,
				  triangles);
		  }
		  if (tri_meshes != null) {
			  tri_meshes.add(new TriMesh(
					  texturePath, // String texture_image,
					  worldXYZ,    // double [][] worldXYZ,
					  texCoord,    // double [][] texCoord,
					  triangles)); // int [][] triangles					  
		  }
	  }

	  // New version with subdivision
	  public static void generateClusterX3d(      // New version with alpha
			  boolean         full_texture,       // true - full size image, false - bounds only
			  int             subdivide_mesh,     // 0,1 - full tiles only, 2 - 2x2 pixels, 4 - 2x2 pixels
			  boolean []      alpha,              // boolean alpha - true - opaque, false - transparent. Full/bounds
			                                      // matching selection
			  double []       dalpha,             // before boolean
			  X3dOutput       x3dOutput,          // output x3d if not null
			  WavefrontExport wfOutput,           // output WSavefront if not null
			  ArrayList<TriMesh> tri_meshes,
			  String          texturePath,
			  String          id,
			  String          class_name,
			  Rectangle       bounds,
			  Rectangle       texture_bounds,     // if not null - allows trimmed combo textures
			  // Below selected and disparity are bounds.width*bounds.height
			  boolean []      selected,           // may be either tilesX * tilesY or bounds.width*bounds.height
			  double []       disparity,          // if null, will use min_disparity
			  int     []      border_int,
			  int             max_border,
			  int             tile_size,
			  int             tilesX,
			  int             tilesY,
			  GeometryCorrection geometryCorrection,
			  boolean         correctDistortions, // requires backdrop image to be corrected also
			  double []       tri_img,   //
			  int             tri_img_width,
			  double          min_disparity,
			  double          max_disparity,
			  double          maxDispTriangle,    // relative <=0 - do not use
			  double          maxZtoXY,           // 10.0. <=0 - do not use
			  double          maxZ,               // far clip (0 - do not clip). Negative - limit by max
			  boolean         limitZ,
//			  double [][]     dbg_disp_tri_slice,
			  int             debug_level,
			  boolean         dbg_plot_center, //  = true;
			  double          dbg_line_color, //  =  1.0;
			  double          dbg_center_color// = 3.0;
			  ) throws IOException
	  {
//		  boolean         show_triangles = tri_img != null; 
		  if (bounds == null) {
			  return; // not used in lwir
		  }
		  boolean display_triangles = debug_level > 0;
		  boolean display_src = debug_level > 1;
		  boolean display_for_mesh = debug_level > 1;
		  
		  
		  
		  if (display_src) {
			  double [][] dbg_img = new double [3][selected.length];
			  for (int i = 0; i < dbg_img[0].length; i++) {
				  dbg_img[0][i] = disparity[i];
				  dbg_img[1][i] = selected[i]? 20.0 : 0.0;
				  dbg_img[2][i] = border_int[i] * 10;
			  }
				ShowDoubleFloatArrays.showArrays(
						dbg_img,
						bounds.width,
						bounds.height,
						true,
						"src_for_triangles",
						new String[] {"disparity", "selected","borders"});
		  }
		  
		  double [][] d_for_mesh = null;
		  if (display_for_mesh) {
				final boolean full_selection = selected.length > (bounds.height * bounds.width); // applies to selected_tiles 
				final boolean full_alpha = alpha.length > (bounds.height * bounds.width * tile_size * tile_size);

			  d_for_mesh =plotForMesh(
					  full_selection, // final boolean    full_selection, // false
					  full_alpha,     // final boolean    full_alpha,     // true
					  bounds,         // final Rectangle  bounds,
					  tri_img_width,  // final int        out_width,
					  tilesX,         // final int        tilesX,
					  tilesY,         // final int        tilesY,
					  tile_size,      // final int        tile_size,		
					  subdivide_mesh, // final int        subdiv,
					  disparity,      // final double []  disparity,
					  selected,       // final boolean [] selection,
					  border_int,     // final int []     border_int,
					  alpha,          // final boolean [] alpha)
					  dalpha);        // final double []  dalpha){

			  if (display_src) {
					ShowDoubleFloatArrays.showArrays(
							d_for_mesh,
							tri_img_width,
							d_for_mesh[0].length / tri_img_width,
							true,
							"fullsize_src_for_triangles",
							new String[] {"disparity", "selected","borders","alpha", "dalpha"});
			  }
		  }
		  
		  /*
		  int [][] indices =  getCoordIndices( // starting with 0, -1 - not selected // updated 09.18.2022
				  bounds, 
				  selected,
				  tilesX); 
		  */
		  int [] pnum_indices = new int[1];
		  /*
		   * Enumerate "large" and "small" tiles, where "large" are actual tiles and "small" are
		   * subdivided (by subdiv in each direction) ones to increase lateral mesh resolution.
		   * sub-tiles are populated if at least one pixel in it is opaque. Input selected_tiles
		   * and alpha arrays may correspond to either full image or rectangular bounds, output
		   * array always corresponds to bounds.
		   */
		  int [][][][] indices =   getCoordIndices( // starting with 0, -1 - not selected
				  bounds,         // final Rectangle  bounds,
				  selected,       // final boolean [] selected_tiles, can not be null
				  tilesX,         // final int        tilesX,
				  tile_size,      // final int        tile_size,		
				  alpha,          // final boolean [] alpha,
				  subdivide_mesh, // final int        subdiv,
				  pnum_indices);  // final int []     num_indices
		  int [][] no_connect = getNoConnect( // 0 - neutral 1 - max_neib_lev, 2 - max_neib_lev+1
				  bounds,         // final Rectangle  bounds,
				  indices,        // final int [][][][] indices,
				  border_int,     // final int []     border_int,
				  max_border,     // final int        max_border,
				  tilesX,         // final int        tilesX,
				  tile_size);     // final int        tile_size,		
		  /*
		   * Convert "large" tiles to arrays of small ones if it has a small-tile neighbor with
		   * gaps along the border with this one
		   */
		  pnum_indices[0] = splitLargeTileIndices(
				  indices,          // int [][][][] indices)
				  subdivide_mesh);  // final int [] num_indices) {
		  
		  Rectangle tex_rect = full_texture ? (
				  (texture_bounds != null) ?
						  new Rectangle( // when combo texture is trimmed (alpha should still be for the full image)
								  bounds.x - texture_bounds.x,
								  bounds.y - texture_bounds.y,
								  texture_bounds.width,
								  texture_bounds.height):
							  new Rectangle( // when combo texture is full tilesX * tilesY
									  bounds.x,
									  bounds.y,
									  tilesX,
									  tilesY)) :
								  null; // when texture is for cluster only
		  /*
		   * Get texture coordinates (0..1) for horizontal (positive - to the right) and vertical (positive - up)
		   */
		  double [][] texCoord = getTexCoords( // get texture coordinates for indices
				  tex_rect, // final Rectangle    rect,     // 0 or full width, full height of the image in tiles
				  pnum_indices[0], // final int          num_indices,
				  indices); // final int [][][][] indices)

		  double [][] worldXYZ = getCoords( // get world XYZ in meters for indices // updated 09.18.2022
				  disparity,          // final double []          disparity,
				  min_disparity,      // final double             min_disparity,
				  max_disparity,      // final double             max_disparity,
				  bounds,             // final Rectangle          bounds,
				  pnum_indices[0],    // final int                num_indices, 
				  indices,            // final int [][][][]       indices,
				  tilesX,             // final int                tilesX,
				  tile_size,          // final int                tile_size,
				  correctDistortions, // requires backdrop image to be corrected also
				  geometryCorrection);// final GeometryCorrection geometryCorrection)
		  /*
		   * Triangulate all vertice indices - combine triangulation of same-size equailateral 45-degree
		   * large (tile size) and small (tile subdivisions) and add connections between large and small ones 
		   */
		  int [][] triangles = 	triangulateAll(
				  indices,           // int [][][][] indices,
				  no_connect);       // int [][]     no_connect)

		  if (triangles.length == 0) {
			  System.out.println("generateClusterX3d(): got NO triangles, do not output 3D mesh");
			  return;
		  } else {
			  System.out.println("generateClusterX3d(): got "+triangles.length+" triangles");

		  }
		  
		  if (tri_img != null) {
			  plotMesh(
					  tri_img,       // final double []   canvas,
					  tri_img_width, // final int         width,
					  texCoord,      // final double [][] tex_coord,
					  triangles,     // final int [][]    triangles,
					  dbg_plot_center,   // final boolean     plot_center,
					  dbg_line_color,    // final double      line_color,
					  dbg_center_color); // final double      center_color)
			  if (display_triangles) {
				  if (d_for_mesh != null) {
					  ShowDoubleFloatArrays.showArrays(
							  new double[][] {d_for_mesh[0], d_for_mesh[1], d_for_mesh[2], d_for_mesh[3], d_for_mesh[4], tri_img},
							  tri_img_width,
							  tri_img.length / tri_img_width,
							  true,
							  "full-src-triangles",
							  new String[] {"disparity", "selected", "borders", "alpha", "dalpha", "mesh"});

				  } else {
					  ShowDoubleFloatArrays.showArrays(
							  tri_img,
							  tri_img_width,
							  tri_img.length / tri_img_width,
							  "this-triangles");
				  }
			  }
		  }
		  
		  if (x3dOutput != null) {
			  x3dOutput.addCluster(
					  texturePath,
					  id,
					  class_name,
					  texCoord,
					  worldXYZ,
					  triangles);
		  }
		  if (wfOutput != null) {
			  wfOutput.addCluster(
				  texturePath,
				  id,
				  texCoord,
				  worldXYZ,
				  triangles);
		  }
		  if (tri_meshes != null) {
			  tri_meshes.add(new TriMesh(
					  texturePath, // String texture_image,
					  worldXYZ,    // double [][] worldXYZ,
					  texCoord,    // double [][] texCoord,
					  triangles)); // int [][] triangles					  
		  }
	  }
	  
	  
}
