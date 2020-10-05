/**
 **
 ** OpticalFlow - Process scene pairs
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  OpticalFlow.java is free software: you can redistribute it and/or modify
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

import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

public class OpticalFlow {
	public static double [] ZERO3 = {0.0,0.0,0.0};
	public static double  LINE_ERR = 0.1;
	public int            threadsMax = 100;  // maximal number of threads to launch
	public boolean        updateStatus = true;

	public OpticalFlow (
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus) {
		this.threadsMax = threadsMax;
		this.updateStatus = updateStatus;
		
	}

	public void fillTilesNans(
			final double [][][] nan_tiles,
			final QuadCLT     qthis,
			final int         num_passes,
			final double      max_change,
			final int         debug_level)
	{
		final double   diagonal_weight = 0.5 * Math.sqrt(2.0); // relative to ortho
		final int margin =               0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         qthis.getTileProcessor();
		final int transform_size =       tp.getTileSize();
		final int fullTileSize =         2 * (transform_size + margin);
///		final int fullTileLen =          fullTileSize * fullTileSize;
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int dbg_mtile = 0; // 203;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iMTile = ai.getAndIncrement(); iMTile < nan_tiles.length; iMTile = ai.getAndIncrement()) {
						if (iMTile == dbg_mtile) {
							System.out.println("iMTile = "+iMTile);
						}
						if (nan_tiles[iMTile] != null) {
							tilesFillNaN(
									neibw, // final double []   neibw,
									nan_tiles[iMTile], // final double [][] slices,
									num_passes, // final int         num_passes,
									max_change, // final double      max_change,
									fullTileSize); // final int         width
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		if (debug_level > 0) {
			// show debug image
			String title = qthis.getImageName()+"-NO-NaN";
			showMacroTiles(
					title,        // String title,
					nan_tiles, // double [][][] source_tiles,
					qthis,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
		}
		System.out.println("fillTilesNans() DONE.");
	}


	public double [][][] prepareSceneTiles(// to match to reference
			// null for {scene,reference}{xyz,atr} uses instances globals 
			final double []   scene_xyz,     // camera center in world coordinates
			final double []   scene_atr,     // camera orientation relative to world frame
//			final double []   reference_xyz, // camera center in world coordinates
//			final double []   reference_atr, // camera orientation relative to world frame
			final QuadCLT     scene_QuadClt,
			final QuadCLT     reference_QuadClt,
			final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
			final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
			final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional tile shifts [-0.5, 0.5)
			final double      tolerance_absolute, // absolute disparity half-range in each tile
			final double      tolerance_relative, // relative disparity half-range in each tile
			final double      occupancy,          // fraction of remaining  tiles (<1.0)
			final int         num_passes,
			final double      max_change,
			final int         debug_level)
	{
		final double   diagonal_weight = 0.5 * Math.sqrt(2.0); // relative to ortho
		final int         margin = 0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         reference_QuadClt.getTileProcessor();
		final double [][] dsrbg_scene =  scene_QuadClt.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int macroX0 =              (tilesX - macroTilesX * transform_size)/2; // distribute extra tiles symmetrically ==0 for 324
		final int macroY0 =              (tilesY - macroTilesY * transform_size)/2; // distribute extra tiles symmetrically (242 - 1 tile above and below)
		final double [][][] scene_tiles = new double [macroTilesX*macroTilesY][][];
		final int fullTileSize =         2 * (transform_size + margin);
		final int fullTileLen =          fullTileSize * fullTileSize; 
//		final int min_remain_center = (int) Math.round(center_occupancy * transform_size * transform_size);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final double [][] pXpYD = new double [tilesX * tilesY][]; // to hold scene pX, pY, disparity, where pX, pY include flowXY (image pre-shift)
		final double []  hist_weights = new double [fullTileSize];
		for (int i = 0; i < transform_size; i++) {
			hist_weights[margin + transform_size/2+ i] = Math.sin(Math.PI * (i +0.5) / transform_size);
		}
		final ErsCorrection ers_reference = reference_QuadClt.getErsCorrection();
		final ErsCorrection ers_scene =     scene_QuadClt.getErsCorrection();
		ers_reference.setupERS(); // just in case - setUP using instance paRAMETERS
		ers_scene.setupERS();     // just in case - setUP using instance paRAMETERS
		final int hist_len = 8;
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 
		final int dbg_mtile = 453; // 500;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					double [][] pXpYD = new double [fullTileLen][];
					double [][] tXtYD = new double [fullTileLen][]; // measured in tiles, not pixels (disparity - still pixels)
					for (int iMTile = ai.getAndIncrement(); iMTile < reference_tiles.length; iMTile = ai.getAndIncrement()) if (reference_tiles[iMTile] != null){
						if (iMTile == dbg_mtile) {
							System.out.println("iMTile = "+iMTile);
						}
						int mtileY = iMTile / macroTilesX; 
						int mtileX = iMTile % macroTilesX;
						int tY0 = mtileY * transform_size + macroY0 -transform_size/2 - margin;
						int tX0 = mtileX * transform_size + macroX0 -transform_size/2 - margin;
//						Arrays.fill(pXpYD, null);
						Arrays.fill(tXtYD, null);
						for (int iY = 0; iY < fullTileSize; iY++) {
							int tileY = tY0 + iY;
							if ((tileY >= 0) && (tileY < tilesY)) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int tileX = tX0 + iX;
									if ((tileX >= 0) && (tileX < tilesX)) {
//										int nTile = tileX + tileY * tilesX;
										int iTile = iX + iY * fullTileSize;
										
										double disparity = reference_tiles[iMTile][QuadCLT.DSRBG_DISPARITY][iTile];
										if (disparity < 0) {
											disparity = 0.0; // is it needed or should it work with negative too?
										}
										double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
										double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
										if (!Double.isNaN(disparity)) {
											double [] xyd = ers_reference.getImageCoordinatesERS(
													scene_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
													centerX,        // double px,                // pixel coordinate X in this camera view
													centerY,        //double py,                // pixel coordinate Y in this camera view
													disparity,      // double disparity,         // this view disparity 
													true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
													ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
													ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
													true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
													scene_xyz,     // double [] camera_xyz,     // camera center in world coordinates
													scene_atr,     // double [] camera_atr,     // camera orientation relative to world frame
													LINE_ERR);       // double    LINE_ERR)       // threshold error in scan lines (1.0)
											if (xyd != null) {
												tXtYD[iTile] = new double [] {
														((xyd[0] + flowXY[iMTile][0])/transform_size),
														((xyd[1] + flowXY[iMTile][0])/transform_size),
														xyd[2]};
											} else {
												tXtYD[iTile] = null;
											}
										} else {
											tXtYD[iTile] = null;
										}
									}
								}
							}
						}
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							String [] dbg_titles= {"tX","tY","disp"};
							String dbg_title= "tXtYD-MX"+mtileX+"_MY"+mtileY;
							double [][] dbg_img = new double [3][fullTileLen];
							for (int nt =0; nt < fullTileLen; nt++) {
								if(tXtYD[nt] != null) {
									for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = tXtYD[nt][i];
								} else {
									for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = Double.NaN;
								}
							}
							(new ShowDoubleFloatArrays()).showArrays(
									dbg_img,
									fullTileSize,
									fullTileSize,
									true,
									dbg_title,
									dbg_titles);
						}
						
						// Find best fractional pixel offset
						double [] hist_fx = new double [hist_len];
						double [] hist_fy = new double [hist_len];
						for (int iY = 0; iY < fullTileSize; iY++) if (hist_weights[iY] > 0.0){
							for (int iX = 0; iX < fullTileSize; iX++) if (hist_weights[iX] > 0.0){
								int iTile = iX + iY * fullTileSize;
								if (tXtYD[iTile] != null) {
									int hidx = (int) Math.floor((tXtYD[iTile][0] - Math.floor(tXtYD[iTile][0])) * hist_len);
									if (hidx < 0) hidx = 0; 
									else if (hidx >= hist_len) hidx = hist_len - 1;
									hist_fx[hidx] += hist_weights[iX];
									int hidy = (int) Math.floor((tXtYD[iTile][1] - Math.floor(tXtYD[iTile][1])) * hist_len);
									if (hidy < 0) hidy = 0; 
									else if (hidy >= hist_len) hidy = hist_len - 1;
									hist_fy[hidy] += hist_weights[iY];
								}
							}
						}
						int hist_fx_mx = 0;
						int hist_fy_mx = 0;
						for (int i = 1; i < hist_len; i++) {
							if (hist_fx[i] > hist_fx[hist_fx_mx]) hist_fx_mx = i;
							if (hist_fy[i] > hist_fy[hist_fy_mx]) hist_fy_mx = i;
						}
						double offsX = (0.5 + hist_fx_mx) / hist_len;
						double offsY = (0.5 + hist_fy_mx) / hist_len;
						double swx = 0.0, swy = 0.0, sx =  0.0, sy = 0.0;
						for (int iY = 0; iY < fullTileSize; iY++) {
							for (int iX = 0; iX < fullTileSize; iX++) {
								int iTile = iX + iY * fullTileSize;
								if (tXtYD[iTile] != null) {
									double wx = hist_weights[iX];
									double wy = hist_weights[iX];
									double dx = tXtYD[iTile][0] - offsX;
									double dy = tXtYD[iTile][1] - offsY;
									dx = dx - Math.round(dx);
									dy = dy - Math.round(dy);
									swx += wx;
									swy += wy;
									sx += wx * dx;
									sy += wy * dy;
								}
							}
						}
						offsX += sx/swx;
						offsY += sy/swy;
						// use offsX, offsY as fractional shift and for data interpolation
						if (offsX >= .5) offsX -= 1.0;
						if (offsY >= .5) offsY -= 1.0;
						flowXY_frac[iMTile] = new double [] {offsX, offsY};
						double min_tX = Double.NaN, max_tX = Double.NaN, min_tY = Double.NaN, max_tY = Double.NaN;
						for (int iY = 0; iY < fullTileSize; iY++) {
							for (int iX = 0; iX < fullTileSize; iX++) {
								int iTile = iX + iY * fullTileSize;
								if (tXtYD[iTile] != null) {
									tXtYD[iTile][0]-=offsX;
									tXtYD[iTile][1]-=offsY;
									if (Double.isNaN(min_tX)) {
										min_tX =tXtYD[iTile][0]; 
										min_tY =tXtYD[iTile][1];
										max_tX = min_tX;
										max_tY = min_tY;
									}
									if (min_tX > tXtYD[iTile][0]) min_tX = tXtYD[iTile][0];
									if (min_tY > tXtYD[iTile][1]) min_tY = tXtYD[iTile][1];
									if (max_tX < tXtYD[iTile][0]) max_tX = tXtYD[iTile][0];
									if (max_tY < tXtYD[iTile][1]) max_tY = tXtYD[iTile][1];
								}
							}
						}
						int imin_tX = (int) Math.floor(min_tX);
						int imin_tY = (int) Math.floor(min_tY);
						int imax_tX = (int) Math.ceil (max_tX);
						int imax_tY = (int) Math.ceil (max_tY);
						// See if at least some of fits into the frame
						if ((imin_tX >= tilesX) || (imin_tY >= tilesY) || (imax_tX < 0)  || (imax_tX < 0)) {
							continue; // no overlap at all
						}
						int iwidth =  imax_tX - imin_tX + 1;
						int iheight = imax_tY - imin_tY + 1;
						double [][] scene_slices = new double [dsrbg_scene.length][iwidth*iheight];
						for (int iY = 0; iY < iheight; iY++) {
							int tY = imin_tY + iY;
							if ((tY >= 0) && (tY < tilesY)) {
								for (int iX = 0; iX < iwidth; iX++) {
									int tX = imin_tX + iX;
									if ((tX >= 0) && (tX < tilesX)) {
										int iTile = iX + iY * iwidth;
										int tile = tX + tilesX * tY;
										if (dsrbg_scene[QuadCLT.DSRBG_STRENGTH][tile] > 0.0) {
											double d = dsrbg_scene[QuadCLT.DSRBG_DISPARITY][tile];
											scene_slices[QuadCLT.DSRBG_DISPARITY][iTile] = d;
											if (!Double.isNaN(d)) {
												scene_slices[QuadCLT.DSRBG_STRENGTH][iTile] = dsrbg_scene[QuadCLT.DSRBG_STRENGTH][tile];
											}
										}
									}
								}
							}
						}
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							String [] dbg_titles= {"d","s","r","b","g"};
							String dbg_title= "before_disparity-MX"+mtileX+"_MY"+mtileY;
							(new ShowDoubleFloatArrays()).showArrays(
									scene_slices,
									iwidth,
									iheight,
									true,
									dbg_title,
									dbg_titles);
						}
						// filter by disparity - 1 tile around rounded
						final TileNeibs tn =  new TileNeibs(iwidth, iheight);
						for (int iY = 0; iY < fullTileSize; iY++) {
							for (int iX = 0; iX < fullTileSize; iX++) {
								int iTile = iX + iY * fullTileSize;
								if (tXtYD[iTile] != null) {
									double disp_tolerance = tolerance_absolute + tXtYD[iTile][2] * tolerance_relative;
									double disp_min = tXtYD[iTile][2] - disp_tolerance;
									double disp_max = tXtYD[iTile][2] + disp_tolerance;
									int nt = tn.getIndex(
											(int) Math.round(tXtYD[iTile][0]) - imin_tX,
											(int) Math.round(tXtYD[iTile][1]) - imin_tY);
									if (nt >= 0) { // should always be
										for (int dir = 0; dir <9; dir++) {
											int nt1 = tn.getNeibIndex(nt, dir);
											if ((nt1 >= 0) && (scene_slices[QuadCLT.DSRBG_STRENGTH][nt1] > 0.0)) {
												if (    (scene_slices[QuadCLT.DSRBG_DISPARITY][nt1] < disp_min) ||
														(scene_slices[QuadCLT.DSRBG_DISPARITY][nt1] > disp_max)) {
													scene_slices[QuadCLT.DSRBG_STRENGTH][nt1] = 0.0; // disable tile
												}
											}
										}
									}
								}
							}
						}
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							String [] dbg_titles= {"d","s","r","b","g"};
							String dbg_title= "after_disparity-MX"+mtileX+"_MY"+mtileY;
							(new ShowDoubleFloatArrays()).showArrays(
									scene_slices,
									iwidth,
									iheight,
									true,
									dbg_title,
									dbg_titles);
						}
						
						// copy rest of the data (all but strength) where strength > 0
						for (int iY = 0; iY < iheight; iY++) {
							int tY = imin_tY + iY;
							if ((tY >= 0) && (tY < tilesY)) {
								for (int iX = 0; iX < iwidth; iX++) {
									int tX = imin_tX + iX;
									if ((tX >= 0) && (tX < tilesX)) {
										int iTile = iX + iY * iwidth;
										int tile = tX + tilesX * tY;
										if (scene_slices[QuadCLT.DSRBG_STRENGTH][iTile] > 0.0) {
											for (int i = 0; i < scene_slices.length; i++) if (i != QuadCLT.DSRBG_STRENGTH) {
												scene_slices[i][iTile] = dsrbg_scene[i][tile];
											}
										}
									}
								}
							}
						}
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							String [] dbg_titles= {"d","s","r","b","g"};
							String dbg_title= "scene1-MX"+mtileX+"_MY"+mtileY;
							(new ShowDoubleFloatArrays()).showArrays(
									scene_slices,
									iwidth,
									iheight,
									true,
									dbg_title,
									dbg_titles);
						}
						
						// set NAN everywhere where strength is 0 (including border tiles
						int num_dead = 0;
						for (int iTile = 0; iTile < scene_slices[0].length; iTile++) {
							if (scene_slices[QuadCLT.DSRBG_STRENGTH][iTile] <=0.0) {
								for (int i = 0; i < scene_slices.length; i++) {
									scene_slices[i][iTile] = Double.NaN;
								}
								num_dead++;
							}
						}
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							System.out.println("scene2-MX"+mtileX+"_MY"+mtileY+" num_dead="+num_dead);
							String [] dbg_titles= {"d","s","r","b","g"};
							String dbg_title= "scene2-MX"+mtileX+"_MY"+mtileY;
							(new ShowDoubleFloatArrays()).showArrays(
									scene_slices,
									iwidth,
									iheight,
									true,
									dbg_title,
									dbg_titles);
						}
						
						
						//center_occupancy -> all occupancy						
						double fract_active = 1.0*(scene_slices[0].length - num_dead)/scene_slices[0].length; // all tiles, not center
						if (fract_active < occupancy) {
							flowXY_frac[iMTile] = null;
							continue;
						}
						// Here need to fill NaNs, then
						
						tilesFillNaN(
								neibw,        // final double []   neibw,
								scene_slices, // final double [][] slices,
								num_passes,   // final int         num_passes,
								max_change,   // final double      max_change,
								iwidth);      // final int         width
						
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							String [] dbg_titles= {"d","s","r","b","g"};
							String dbg_title= "scene3NaN-MX"+mtileX+"_MY"+mtileY;
							(new ShowDoubleFloatArrays()).showArrays(
									scene_slices,
									iwidth,
									iheight,
									true,
									dbg_title,
									dbg_titles);
						}
						
// bi-linear interpolate						
						double [][] scene_mapped = new double [dsrbg_scene.length][fullTileLen];
						boolean need_nan_filter = false;;
						for (int iY = 0; iY < fullTileSize; iY++) {
							for (int iX = 0; iX < fullTileSize; iX++) {
								int iTile = iX + iY * fullTileSize;
								if (tXtYD[iTile] != null) {
									int itX = (int) Math.floor(tXtYD[iTile][0]);
									int itY = (int) Math.floor(tXtYD[iTile][1]);
									double kX = tXtYD[iTile][0]-itX;
									double kY = tXtYD[iTile][1]-itY;
									itX -= imin_tX; // relative to scene_slices = new double [dsrbg_scene.length][iwidth*iheight];
									itY -= imin_tY; // relative to scene_slices = new double [dsrbg_scene.length][iwidth*iheight];
									int indx = itX + itY * iwidth; 
									for (int i = 0; i < scene_mapped.length; i++) {
										scene_mapped[i][iTile] = 
												(1.0 - kY) * ((1.0-kX) * scene_slices[i][indx] +          kX * scene_slices[i][indx +          1])+
												(      kY) * ((1.0-kX) * scene_slices[i][indx + iwidth] + kX * scene_slices[i][indx + iwidth + 1]);
									}									
								} else {
									for (int i = 0; i < scene_mapped.length; i++) {
										scene_mapped[i][iTile] = Double.NaN; // will need to filter?
									}
									need_nan_filter=true;
								}
							}
						}
						if ((debug_level>1) && (iMTile == dbg_mtile)) {
							String [] dbg_titles= {"d","s","r","b","g"};
							String dbg_title= "mapped-MX"+mtileX+"_MY"+mtileY;
							(new ShowDoubleFloatArrays()).showArrays(
									scene_mapped,
									fullTileSize,
									fullTileSize,
									true,
									dbg_title,
									dbg_titles);
						}
						
						if (need_nan_filter) {
							tilesFillNaN(
									neibw,         // final double []   neibw,
									scene_mapped,  // final double [][] slices,
									num_passes,    // final int         num_passes,
									max_change,    // final double      max_change,
									fullTileSize); // final int         width
							if ((debug_level>1) && (iMTile == dbg_mtile)) {
								String [] dbg_titles= {"d","s","r","b","g"};
								String dbg_title= "mappedNaN-MX"+mtileX+"_MY"+mtileY;
								(new ShowDoubleFloatArrays()).showArrays(
										scene_mapped,
										fullTileSize,
										fullTileSize,
										true,
										dbg_title,
										dbg_titles);
							}
							
						}
						scene_tiles[iMTile] = scene_mapped;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		
		if (debug_level > 0) {
			
			
			// show debug image
			String title = scene_QuadClt.getImageName()+"-scene_tiles";
			showMacroTiles(
					title,        // String title,
					scene_tiles, // double [][][] source_tiles,
					scene_QuadClt,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
			
			
			String [] dbg_titles= {"dX","dY"};
			String dbg_title= "flowXY_frac";
			double [][] dbg_img = new double [2][macroTilesX*macroTilesY];
			for (int nt =0; nt < flowXY_frac.length; nt++) {
				if(flowXY_frac[nt] != null) {
					for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = flowXY_frac[nt][i];
				} else {
					for (int i=0; i < dbg_img.length; i++) dbg_img[i][nt] = Double.NaN;
				}
			}
			if (debug_level > 1) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					macroTilesX,
					macroTilesY,
					true,
					dbg_title,
					dbg_titles);
			}
			
			
		}
		System.out.println("fillTilesNans() DONE.");
		return scene_tiles;
		
	}
	
	
	// helper to be called from thread
	private void tilesFillNaN(
			final double []   neibw,
			final double [][] slices,
			final int         num_passes,
			final double      max_change,
			final int         width
			) {
		final int tiles = slices[0].length;
		final int height = tiles/width;
		double [] strength = slices[QuadCLT.DSRBG_STRENGTH]; 
		final TileNeibs tn =  new TileNeibs(width, height);
		double [] slice_in =  new double [tiles];
		double [] slice_out = new double [tiles];
		
		//first - process strength, then use calculated strength to fill other slices
		boolean [] fixed = new boolean [tiles]; 
		for (int i = 0; i < tiles; i++) {
//			if (!Double.isNaN(slices[QuadCLT.DSRBG_DISPARITY][i]) &&  (strength[i] > 0.0)) {
			if (strength[i] > 0.0) {
				fixed[i] = true;
			} else {
				strength[i] = 0.0; // get rid of NaN; for non-strength will use average
			}
		}
		System.arraycopy(strength, 0, slice_in, 0, tiles);
		for (int pass = 0; pass < num_passes; pass ++) {
			double pass_change = 0.0;
			for (int nt = 0; nt < tiles; nt++) if (!fixed[nt]){
				double s = 0.0;
				double sw = 0.0;
				double d;
				for (int dir = 0; dir < 8; dir++) {
					int nt1 = tn.getNeibIndex(nt, dir);
					if (nt1 >=0) {
						if (fixed[nt1]) {
							d = strength[nt1];
						}else {
							d = slice_in[nt1];
						}
						s += d * neibw[dir];
						sw += neibw[dir];
					}
				}
				// sw should never be 0;
				s /= sw;
				pass_change = Math.max(pass_change, Math.abs(slice_out[nt] - s));
				slice_out[nt] = s;
			}
			if (pass_change < max_change) {
				break;
			}
			System.arraycopy(slice_out, 0, slice_in, 0, tiles);
		}
		for (int i = 0; i < fixed.length; i++) if (!fixed[i]){
			strength[i] = slice_out[i];
		}
		//non-strength							
		for (int iSlice = 0; iSlice < slices.length; iSlice++) if (iSlice != QuadCLT.DSRBG_STRENGTH){
			double [] slice =    slices[iSlice];
			System.arraycopy(slice, 0, slice_in, 0, tiles);
			double fs =0.0;
			double fsw = 0.0;
			for (int i = 0; i < fixed.length; i++) {
				if (!Double.isNaN(slice[i])  &&  (strength[i] > 0.0)) { //  - now already non-null
					fixed[i] = true;
					fs +=  slice[i] * strength[i];
					fsw += strength[i];
				}
			}
			if (fsw <= 0.0) {
				continue; // should not happen
			}
			fs /= fsw; // average value
			for (int i = 0; i < fixed.length; i++) if (! fixed[i]){
				slice_in[i] = fs;
			}								
			for (int pass = 0; pass < num_passes; pass ++) {
				double pass_change = 0.0;
				for (int nt = 0; nt < tiles; nt++) if (!fixed[nt]){
					double s = 0.0;
					double sw = 0.0;
					double d;
					for (int dir = 0; dir < 8; dir++) {
						int nt1 = tn.getNeibIndex(nt, dir);
						if (nt1 >=0) {
							if (fixed[nt1]) {
								d = slice[nt1];
							}else {
								d = slice_in[nt1];
							}
							double w = neibw[dir]; //  * strength[nt1];
							s += d * w ;
							sw += w;
						}
					}
					if (sw > 0) {
						s /= sw;
					}
					pass_change = Math.max(pass_change, Math.abs(slice_out[nt] - s));
					slice_out[nt] = s;
				}
				if (pass_change < max_change) {
					break;
				}
				System.arraycopy(slice_out, 0, slice_in, 0, tiles);
			}
			for (int i = 0; i < fixed.length; i++) if (!fixed[i]){
				slice[i] = slice_out[i];
			}
		}
	}
	
	
/*
	public double [] getImageCoordinatesERS(
			QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
			double px,                // pixel coordinate X in this camera view
			double py,                // pixel coordinate Y in this camera view
			double disparity,         // this view disparity 
			boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
			double [] reference_xyz,  // this view position in world coordinates (typically zero3)
			double [] reference_atr,  // this view orientation relative to world frame  (typically zero3)
			boolean distortedCamera,  // camera view is distorted (false - rectilinear)
			double [] camera_xyz,     // camera center in world coordinates
			double [] camera_atr,     // camera orientation relative to world frame
			double    line_err)       // threshold error in scan lines (1.0)
	
 */
	
	public double [][][] prepareReferenceTiles(
			final QuadCLT     qthis,
//			final int         margin, // extra margins over 16x16 tiles to accommodate distorted destination tiles
			final double      tolerance_absolute, // absolute disparity half-range in each tile
			final double      tolerance_relative, // relative disparity half-range in each tile
			final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
			final int         debug_level)
	{
		final int margin =               0; // 1; // extra margins over 16x16 tiles to accommodate distorted destination tiles
		final TileProcessor tp =         qthis.getTileProcessor();
		final double [][] dsrbg =        qthis.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int macroX0 =              (tilesX - macroTilesX * transform_size)/2; // distribute extra tiles symmetrically ==0 for 324
		final int macroY0 =              (tilesY - macroTilesY * transform_size)/2; // distribute extra tiles symmetrically (242 - 1 tile above and below)
		final double [][][] source_tiles = new double [macroTilesX*macroTilesY][][];
		final int fullTileSize =         2 * (transform_size + margin);
		final int fullTileLen =          fullTileSize * fullTileSize; 
		final int min_remain_center = (int) Math.round(center_occupancy * transform_size * transform_size);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		// indices of the center 8x8 of the full 20x20 tile (assuming margin = 4) 
		final Integer [] order_indices = new Integer [transform_size*transform_size];
		for (int i = 0; i <transform_size; i++) {
			int i1 = i + margin + transform_size/2;
			for (int j = 0; j <transform_size; j++) {
				int j1 = j + margin + transform_size/2;
				order_indices[i * transform_size + j] = i1 * fullTileSize + j1;
			}
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					final double [] disparity = new double [fullTileLen];
					final double [] strength =  new double [fullTileLen];
					Integer [] local_indices = new Integer [transform_size*transform_size];
					for (int iMTile = ai.getAndIncrement(); iMTile < source_tiles.length; iMTile = ai.getAndIncrement()) {
						int mtileY = iMTile / macroTilesX; 
						int mtileX = iMTile % macroTilesX;
						int tY0 = mtileY * transform_size + macroY0 -transform_size/2 - margin;
						int tX0 = mtileX * transform_size + macroX0 -transform_size/2 - margin;
						Arrays.fill(strength,  0.0);
						Arrays.fill(disparity, 0.0); // Double.NaN);
						System.arraycopy(order_indices,0,local_indices,0,local_indices.length);
						
						for (int iY = 0; iY < fullTileSize; iY++) {
							int tileY = tY0 + iY;
							if ((tileY >= 0) && (tileY < tilesY)) {
								for (int iX = 0; iX < fullTileSize; iX++) {
									int tileX = tX0 + iX;
									if ((tileX >= 0) && (tileX < tilesX)) {
										int nTile = tileX + tileY * tilesX;
										int iTile = iX + iY * fullTileSize;
										double d = dsrbg[QuadCLT.DSRBG_DISPARITY][nTile];
										double s = dsrbg[QuadCLT.DSRBG_STRENGTH][nTile];
										if (!Double.isNaN(d) && (s > 0.0)){
											disparity[iTile] = d;
											strength[iTile]  = s;
										}
									}
								}
							}
						}
						double sw =   0.0;
						double swd =  0.0;
						int num_remain = 0;
						for (int i = 0; i < local_indices.length; i++) {
							double d = disparity[local_indices[i]];
							double s = strength[local_indices[i]];
							sw += s;
							swd += s*d;
							num_remain++;
						}
						if (num_remain < min_remain_center) {
							continue; // already too few tiles
						}
						Arrays.sort(local_indices, new Comparator<Integer>() {
						    @Override
						    public int compare(Integer lhs, Integer rhs) {
						        // -1 - less than, 1 - greater than, 0 - equal, not inverted for ascending disparity
						        //										return lhs.disparity > rhs.disparity ? -1 : (lhs.disparity < rhs.disparity ) ? 1 : 0;
//						        return disparity[lhs] < disparity[rhs] ? -1 : (disparity[lhs] > disparity[rhs] ) ? 1 : 0;
						    	// modifying to make  NaN greater than all non-NaN and equal to each other
						        return (!(disparity[lhs] >= disparity[rhs])) ? -1 : (!(disparity[lhs] <= disparity[rhs] )) ? 1 : 0;
						    }
						});
						int indx_min = 0;
						int indx_max = num_remain - 1;

						double d_low = Double.NaN;
						double d_high = Double.NaN; 

						while (true) {
							double disp_avg = swd / sw;
							double d_tol = tolerance_absolute + Math.max(disp_avg * tolerance_relative, 0.0);
							d_low =  disp_avg - d_tol;
							d_high = disp_avg + d_tol;
							// see if both min and max are within tolerance
							double d_min = disp_avg - disparity[local_indices[indx_min]];
							double d_max = disparity[local_indices[indx_max]] - disp_avg; 
							
							if ((d_min <= d_tol) && (d_max <= d_tol)) {
								break; // already OK
							}
							num_remain --;
							if (num_remain < min_remain_center) {
								break;
							}
							int indx_gone = -1;
							if (d_min > d_max) {
								indx_gone = indx_min;
								indx_min++;
							} else {
								indx_gone = indx_max;
								indx_max--;
							}
							double d = disparity[local_indices[indx_gone]];
							double s = strength[local_indices[indx_gone]];
							sw  -= s;
							swd -= d * s;
						}
						if (num_remain < min_remain_center) {
							continue; // too few remains in this tile
						}
						// now put all tiles in fullTileSize * fullTileSize that fit in disparity tolerance and positive strength
						source_tiles[iMTile] = new double[dsrbg.length][fullTileLen];
						for (int l = 0; l < dsrbg.length; l++) {
							Arrays.fill(source_tiles[iMTile][l], Double.NaN);
						}
						for (int i = 0; i < fullTileLen; i++) {
							if (!Double.isNaN(disparity[i]) && (strength[i] > 0.0) && (disparity[i] >= d_low) && (disparity[i] <= d_high)) {
								int tileY = tY0 + i / fullTileSize;
								int tileX = tX0 + i % fullTileSize;
								if ((tileY >= 0) && (tileY < tilesY) && (tileX >= 0) && (tileX < tilesX)) {
									int tile = tileY * tilesX + tileX; 
									for (int l = 0; l < dsrbg.length; l++) {
										source_tiles[iMTile][l][i] = dsrbg[l][tile];
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debug_level > 0) {
			// show debug image
			String title = qthis.getImageName()+"-F"+center_occupancy+"-A"+tolerance_absolute+"-R"+tolerance_relative;
			showMacroTiles(
					title,        // String title,
					source_tiles, // double [][][] source_tiles,
					qthis,        // final QuadCLT qthis,
					margin);      // final int     margin); // extra margins over 16x16 tiles to accommodate distorted destination tiles
			
		}

		return source_tiles;
	}
	
	public void showMacroTiles(
			String title,
			double [][][] source_tiles,
			final QuadCLT qthis,
			final int     margin) // extra margins over 16x16 tiles to accommodate distorted destination tiles
	{
		final TileProcessor tp =         qthis.getTileProcessor();
		final double [][] dsrbg =        qthis.getDSRBG();
		final int tilesX =               tp.getTilesX();
		final int tilesY =               tp.getTilesY();
		final int transform_size =       tp.getTileSize();
		final int macroTilesX =          tilesX/transform_size;
		final int macroTilesY =          tilesY/transform_size;
		final int fullTileSize =         2 * (transform_size + margin);
		
		// show debug image
		final int dbg_with =   macroTilesX * (fullTileSize +1) - 1;
		final int dbg_height = macroTilesY * (fullTileSize +1) - 1;
		final double [][] dbg_img = new double [dsrbg.length][dbg_with * dbg_height];
		for (int l = 0; l < dbg_img.length; l++) {
			Arrays.fill(dbg_img[l],  Double.NaN);
		}
		for (int mtile = 0; mtile < source_tiles.length; mtile++) if (source_tiles[mtile] != null){
			int mTileY = mtile / macroTilesX;
			int mTileX = mtile % macroTilesX;
			for (int iY = 0; iY < fullTileSize; iY++) {
				int tileY = (fullTileSize +1) * mTileY + iY;
				for (int iX = 0; iX < fullTileSize; iX++) {
					int tileX = (fullTileSize +1) * mTileX + iX;
					for (int l = 0; l < dbg_img.length; l++) {
						dbg_img[l][tileY * dbg_with + tileX] = source_tiles[mtile][l][iY * fullTileSize + iX];
					}							
				}
			}
		}
		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				dbg_with,
				dbg_height,
				true,
				title,
				dsrbg_titles);
		
	}
	
	
	
	
	public double[][][]  get_pair(
			double k_prev,
			QuadCLT qthis,
			QuadCLT qprev,
			double corr_scale, //  = 0.75
			int debug_level)
	{
		TileProcessor tp = qthis.getTileProcessor();
		final int iscale = 8;
		double ts =        qthis.getTimeStamp();
		double ts_prev =   ts;
		double [] camera_xyz0 = ZERO3.clone();
		double [] camera_atr0 = ZERO3.clone();
		
		ErsCorrection ersCorrection = qthis.getErsCorrection();
		String this_image_name = qthis.getImageName();
		
		System.out.println("\n"+this_image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", this_image_name,
				ersCorrection.ers_wxyz_center[0], ersCorrection.ers_wxyz_center[1],ersCorrection.ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f",	this_image_name,
				ersCorrection.ers_wxyz_center_dt[0], ersCorrection.ers_wxyz_center_dt[1],ersCorrection.ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", this_image_name,
				ersCorrection.ers_wxyz_center_d2t[0], ersCorrection.ers_wxyz_center_d2t[1],ersCorrection.ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", this_image_name,
				ersCorrection.ers_watr_center_dt[0], ersCorrection.ers_watr_center_dt[1],ersCorrection.ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", this_image_name,
				ersCorrection.ers_watr_center_d2t[0], ersCorrection.ers_watr_center_d2t[1],ersCorrection.ers_watr_center_d2t[2] ));
		
		double dt = 0.0;
		if (qprev == null) {
			qprev = qthis;
		}
		if (qprev != null) {
			ts_prev = qprev.getTimeStamp();
			dt = ts-ts_prev;
			if (dt < 0) {
				k_prev = (1.0-k_prev);
			}
			if (Math.abs(dt) > 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
				k_prev = 0.5;
				System.out.println("Non-consecutive frames, dt = "+dt);
			}
			ErsCorrection ersCorrectionPrev = (ErsCorrection) (qprev.geometryCorrection);
			double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
			double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt;
			double [] wxyz_delta = new double[3];
			double [] watr_delta = new double[3];
			for (int i = 0; i <3; i++) {
				wxyz_delta[i] = corr_scale * dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_wxyz_center_dt[i]);
				watr_delta[i] = corr_scale * dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_watr_center_dt[i]);
			}
			camera_xyz0 = wxyz_delta;
			camera_atr0 = watr_delta;
		}
		
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		
		double [][] dsrbg = transformCameraVew( // shifts previous image correctly (right)
				qthis,       // reference
				qprev,       // QuadCLT   camera_QuadClt,
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr, // camera orientation relative to world frame
				iscale);
		double [][][] pair = {qthis.getDSRBG(),dsrbg};
		
		// combine this scene with warped previous one
		if (debug_level > 0) {
			String [] rtitles = new String[2* dsrbg_titles.length];
			double [][] dbg_rslt = new double [rtitles.length][];
			for (int i = 0; i < dsrbg_titles.length; i++) {
				rtitles[2*i] =    dsrbg_titles[i]+"0";
				rtitles[2*i+1] =  dsrbg_titles[i];
				dbg_rslt[2*i] =   pair[0][i];
				dbg_rslt[2*i+1] = pair[1][i];
			}
			String title = this_image_name+"-"+qprev.image_name+"-dt"+dt;
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_rslt,
					tilesX,
					tilesY,
					true,
					title,
					rtitles);
		}
		/* */
		double      tolerance_absolute = 0.25; // absolute disparity half-range in each tile
		double      tolerance_relative = 0.2; // relative disparity half-range in each tile
		double      center_occupancy =   0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)
		int         num_passes = 100;
		double      max_change = 0.005 ;

		double      tolerance_absolute_inter = 0.25; // absolute disparity half-range in each tile
		double      tolerance_relative_inter = 0.2; // relative disparity half-range in each tile
		double      occupancy_inter =         0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)

		
		double [][][] reference_tiles = prepareReferenceTiles(
				qthis,        // final QuadCLT     qthis,
				tolerance_absolute, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				tolerance_relative, // final double      tolerance_relative, // relative disparity half-range in each tile
				center_occupancy,   // final double      center_occupancy,   // fraction of remaining  tiles in the center 8x8 area (<1.0)
				2); // final int         debug_level)
		
		fillTilesNans(
				reference_tiles,          // final double [][][] nan_tiles,
				qthis,                 // final QuadCLT     qthis,
				num_passes,            // final int         num_passes,
				max_change,            // final double      max_change,
				2);                    // final int         debug_level)
		
		double [][] flowXY = new double [reference_tiles.length][2]; // zero pre-shifts
		double [][] flowXY_frac = new double [reference_tiles.length][]; // Will contain fractional X/Y shift for CLT
		
		double [][][] scene_tiles = prepareSceneTiles(// to match to reference
				// null for {scene,reference}{xyz,atr} uses instances globals 
				camera_xyz0,              // final double []   scene_xyz,     // camera center in world coordinates
				camera_atr0,              // final double []   scene_atr,     // camera orientation relative to world frame
				qprev,                    // final QuadCLT     scene_QuadClt,
				qthis,                    // final QuadCLT     reference_QuadClt,
				reference_tiles,          // final double [][][] reference_tiles, // prepared with prepareReferenceTiles() + fillTilesNans();
				flowXY,                   // final double [][] flowXY, // per macro tile {mismatch in image pixels in X and Y directions
				flowXY_frac,              // final double [][] flowXY_frac, // should be initialized as [number of macro tiles][] - returns fractional shifts [-0.5, 0.5)
				tolerance_absolute_inter, // final double      tolerance_absolute, // absolute disparity half-range in each tile
				tolerance_relative_inter, // final double      tolerance_relative, // relative disparity half-range in each tile
				occupancy_inter,          // final double      occupancy,          // fraction of remaining  tiles (<1.0)
				num_passes,               // final int         num_passes,
				max_change,               // final double      max_change,
				2);                       // final int         debug_level)
		/* */
		return pair;
	}

	
	
	public double [][] transformCameraVew(
			QuadCLT   reference_QuadClt,
			QuadCLT   camera_QuadClt,
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr, // camera orientation relative to world frame
			int       iscale)
	{
//		double    line_err = 0.1; // 10.0; // 0.1; // BUG
		TileProcessor tp = reference_QuadClt.getTileProcessor();
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		int rel_num_passes = 10;
		int num_passes =    transform_size; // * 2;

//		double [] zero3 =               {0.0,0.0,0.0};	
		int stilesX = iscale*tilesX; 
		int stilesY = iscale*tilesY;
		int stiles = stilesX*stilesY;
		double sigma = 0.5 * iscale;
		double scale =  1.0 * iscale/transform_size;
		double [][] dsrbg_camera =    camera_QuadClt.getDSRBG();
		double [][] dsrbg_reference = reference_QuadClt.getDSRBG();
		double [][] ds =        new double [dsrbg_camera.length][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
			}
		}
		
		ErsCorrection ersReferenceCorrection = reference_QuadClt.getErsCorrection();
		ersReferenceCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
		double [] zbuffer = new double [tiles];
		for (int tileY = 0; tileY < tilesY; tileY++) {
//			int stileY = iscale * tileY + iscale/2;
			for (int tileX = 0; tileX < tilesX; tileX++) {
//				int stileX = iscale * tileX + iscale/2;
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_camera[QuadCLT.DSRBG_DISPARITY][nTile];
//				double disparity = dsrbg_reference[QuadCLT.DSRBG_DISPARITY][nTile];
				if (disparity < 0) {
					disparity = 0.0;
				}
				// found that there are tiles with strength == 0.0, while disparity is not NaN
				if (!Double.isNaN(disparity) && (dsrbg_camera[QuadCLT.DSRBG_STRENGTH][nTile] > 0.0)) {
//				if (!Double.isNaN(disparity) && (dsrbg_reference[QuadCLT.DSRBG_STRENGTH][nTile] > 0.0)) {
					/*
					double [] pXpYD = ersReferenceCorrection.getImageCoordinatesERS( // ersCorrection - reference
					camera_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
					centerX,        // double px,                // pixel coordinate X in the reference view
					centerY,        // double py,                // pixel coordinate Y in the reference view
					disparity,      // double disparity,         // reference disparity 
					true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
					ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
					ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
					true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
					camera_xyz,     // double [] camera_xyz,     // camera center in world coordinates
					camera_atr,     // double [] camera_atr,     // camera orientation relative to world frame
					LINE_ERR);       // double    line_err)       // threshold error in scan lines (1.0)
*/					
					double [] pXpYD = ersReferenceCorrection.getImageCoordinatesReferenceERS( // ersCorrection - reference
							
							camera_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
							centerX,        // double px,                // pixel coordinate X in the reference view
							centerY,        // double py,                // pixel coordinate Y in the reference view
							disparity,      // double disparity,         // reference disparity 
							true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
							ZERO3,          // double [] reference_xyz,  // this view position in world coordinates (typically ZERO3)
							ZERO3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically ZERO3)
							true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
							camera_xyz,     // double [] camera_xyz,     // camera center in world coordinates
							camera_atr,     // double [] camera_atr,     // camera orientation relative to world frame
							LINE_ERR);       // double    line_err)       // threshold error in scan lines (1.0)
					if (pXpYD != null) {
						int px = (int) Math.round(pXpYD[0]/transform_size);
						int py = (int) Math.round(pXpYD[1]/transform_size);
						int spx = (int) Math.round(pXpYD[0]*scale);
						int spy = (int) Math.round(pXpYD[1]*scale);
						if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
							//Z-buffer
							if (!(pXpYD[2] < zbuffer[px + py* tilesX])) {
								zbuffer[px + py* tilesX] = pXpYD[2];
								if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
									int sTile = spx + spy* stilesX;
									ds[QuadCLT.DSRBG_DISPARITY][sTile] = pXpYD[2]; //reduce*
									for (int i = QuadCLT.DSRBG_STRENGTH; i < dsrbg_camera.length; i++) {
										ds[i][sTile] = dsrbg_camera[i][nTile]; // reduce * 
									}
								}								
							}
						}
					}
				}
			}
		}
		
		//dsrbg_out[DSRBG_DISPARITY]
		for (int i = 0; i < ds.length; i++) {
			ds[i] = (new DoubleGaussianBlur()).blurWithNaN(
					ds[i], // double[] pixels,
					null,  // double [] in_weight, // or null
					stilesX, // int width,
					stilesY, // int height,
					sigma, // double sigmaX,
					sigma, // double sigmaY,
					0.01); // double accuracy);
		}
		double [][] dsrbg_out = new double [dsrbg_camera.length][tiles];
		int [][] num_non_nan = new int [dsrbg_out.length] [tiles];
		
		for (int stileY = 0; stileY < stilesY; stileY++) {
			int tileY = stileY / iscale; 
			for (int stileX = 0; stileX < stilesX; stileX++) {
				int tileX = stileX / iscale;
				int stile = stileX + stileY * stilesX;
				int tile =  tileX +  tileY *  tilesX;
				for (int i = 0; i < dsrbg_out.length; i++) {
					double d = ds[i][stile];
					if (!Double.isNaN(d)) {
						num_non_nan[i][tile] ++;
						dsrbg_out[i][tile] += d;
					}	
				}
			}
		}
		for (int i = 0; i < dsrbg_out.length; i++) {
			for (int j = 0; j < tiles; j++) {
				if (num_non_nan[i][j] == 0) {
					dsrbg_out[i][j] = Double.NaN;
				} else {
					dsrbg_out[i][j]/=num_non_nan[i][j];
				}
			}
		}

		if (num_passes > 0) {
			for (int i = 0; i < dsrbg_out.length; i++) {
//				dsrbg_out[i] = reference_QuadClt.fillNaNGaps(dsrbg_out[i], num_passes, rel_num_passes, 100); // threadsMax);
				dsrbg_out[i] = tp.fillNaNs(
						dsrbg_out[i],                  // double [] data,
						tilesX,                        //int       width, 
						2 * num_passes,                // int       grow,
						0.5 * Math.sqrt(2.0),          // double    diagonal_weight, // relative to ortho
						num_passes * rel_num_passes,   // int       num_passes,
						threadsMax);                   // final int threadsMax) // maximal number of threads to launch                         
			}
		}
		return dsrbg_out;
	}
	

}
