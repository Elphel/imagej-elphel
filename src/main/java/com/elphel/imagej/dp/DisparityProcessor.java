package com.elphel.imagej.dp;
/**
 **
 ** DisparityProcessor
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  DisparityProcessor.java is free software: you can redistribute it and/or modify
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

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.showDoubleFloatArrays;

public class DisparityProcessor {
	static int [] corn_side_neib = { // of +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			0b00011100,    // top left corner
			0b01111100,    // top middle
			0b01110000,    // top right
			0b00011111,    // middle left
			0b11111111,    // middle
			0b11110001,    // middle right
			0b00000111,    // bottom left
			0b11000111,    // bottom middle
			0b11000001};   // bottom right

	TileProcessor tp;
	//   disparity*scale_dz_dx - disparity difference between neighbor tiles to have 45 degree XZ
	double scale_dz_dx;  // == tile_size* ( 0.001 * this.pixelSize) / this.focalLength
	public DisparityProcessor (
			TileProcessor tp,
			double scale_dz_dx)
	{
		this.tp = tp;
		this.scale_dz_dx = scale_dz_dx;
	}

	public int [] getNeighbors( // creates neighbors mask from bitmask
			int [] surf_indices,
			int    indx,
			int tilesX)
	{
		boolean [] selected = new boolean [surf_indices.length];
		for (int i = 0; i < surf_indices.length; i++){
			selected[i] = surf_indices[i] == indx;
		}
		return getNeighbors(selected, tilesX);
	}
	
	
	public int [] getNeighbors( // creates neighbors mask from bitmask
			boolean [] selected,
			int tilesX
			)
	{
		int [] neibs = new int [selected.length];
		int tilesY = selected.length/tilesX;
		final int [] dirs8 = {-tilesX,  -tilesX + 1, 1, tilesX +1, tilesX, tilesX - 1, -1, -tilesX - 1};

		for (int nTile = 0; nTile < selected.length; nTile++) {
			int tileY = nTile/tilesX; 
			int tileX = nTile - (tileY * tilesX);
			int tileType = 4; 
			if (tileY == 0){
				if (tileX == 0){
					tileType = 0;  
				} else if (tileX == (tilesX - 1)) {
					tileType = 2;  
				} else {
					tileType = 1;  
				}
			} else if (tileY == (tilesY - 1)) {
				if (tileX == 0){
					tileType = 6;  
				} else if (tileX == (tilesX - 1)) {
					tileType = 8;  
				} else {
					tileType = 7;  
				}
			} else {
				if (tileX == 0){
					tileType = 3;  
				} else if (tileX == (tilesX - 1)) {
					tileType = 5;  
				}
			}
			if (selected[nTile]){
				neibs[nTile] = corn_side_neib[tileType];
				for (int i = 0; i < 8; i++) {
					int b = 1 << i;
					if (((neibs[nTile] & b) != 0) && !selected [nTile + dirs8[i]]){
						neibs[nTile] &= ~b;
					}
				}
			}
			else {
				neibs[nTile] = 0; // change to  -1?
			}
		}
		return neibs;
	}

	public double [] dbgShowNeighbors(
			boolean [] selected,
			int [] neighbors,
			int    tile_size,
			double bgnd,
			double fgnd)
	{
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int [][] dirXY8 = {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}};
		int width = tilesX * tile_size;
		int height = tilesY * tile_size;
		double [] rslt = new double [width*height];
		for (int i = 0; i < rslt.length; i++) rslt[i] = bgnd;
		for (int nTile = 0; nTile < neighbors.length; nTile++) {
			if ((neighbors [nTile] >=0) && ((selected == null) || selected[nTile])) { //  if (neighbors[nTile] != 0) {
				int tileY = nTile / tilesX;
				int tileX = nTile % tilesX;
				for (int i = -1; i<= 1; i++){
					for (int j = -1; j<= 1; j++){
						rslt[(tileY*tile_size + tile_size/2 + i) * width + (tileX*tile_size + tile_size/2 + j)] = fgnd;
					}
				}
				for (int ib = 0; ib < dirXY8.length; ib++) if ((neighbors[nTile] & (1 << ib)) != 0){
					for (int i = 0; i <= tile_size/2; i++){
						int x = tileX*tile_size + tile_size/2 + i * dirXY8[ib][0];
						int y = tileY*tile_size + tile_size/2 + i * dirXY8[ib][1];
						rslt[y*width+x] = fgnd;
					}
				}
			}
		}
		return rslt;
	}

	public double [][] dbgShowOverlaps(
//			boolean [] selected,
			int [][] flaps,
			int    tile_size,
			double bgnd,
			double fgnd)
	{
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int [][] dirXY8 = {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}};
		int width = tilesX * tile_size;
		int height = tilesY * tile_size;
		double [][] rslt = new double [8][width*height];
		for (int i = 0; i < rslt[0].length; i++) {
			for (int l=0; l<rslt.length; l++)	rslt[l][i] = bgnd;
		}
		for (int nTile = 0; nTile < flaps.length; nTile++) if (flaps[nTile]!=null) {
			int tileY = nTile / tilesX;
			int tileX = nTile % tilesX;
			for (int l = 0; l < 8; l++) {
				for (int i = -1; i<= 1; i++){
					for (int j = -1; j<= 1; j++){
						rslt[l][(tileY*tile_size + tile_size/2 + i) * width + (tileX*tile_size + tile_size/2 + j)] = fgnd;
					}
				}
				for (int ibr = 0; ibr < dirXY8.length; ibr++){
					if ((flaps[nTile][l] & (1 << ibr)) != 0){
						int ib = (ibr + 4) % 8;
						for (int i = 0; i <= tile_size/2; i++){
							int x = tileX*tile_size + tile_size/2 + i * dirXY8[ib][0];
							int y = tileY*tile_size + tile_size/2 + i * dirXY8[ib][1];
							rslt[l][y*width+x] = fgnd;
						}
					}
				}
			}
		}
		return rslt;
	}
	
	
	
	public double [] dbgShowStress(
			double [][] stress,
			int    tile_size)
	{
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
//		int [][] dirXY8 = {{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1}};
		int width = tilesX * tile_size;
		int height = tilesY * tile_size;
		double [] rslt = new double [width*height];
		for (int nTile = 0; nTile < stress[0].length; nTile++) {
			int tileY = nTile / tilesX;
			int tileX = nTile % tilesX;
			if ((tileY < (tilesY -1)) && (tileX < (tilesX -1))) {
				for (int i = 0; i <= tile_size/2; i ++){
					int x0 = tileX*tile_size +     tile_size/2 + i;
					int x1 = (tileX+1)*tile_size + tile_size/2 - i;
					int y0 = tileY*tile_size +     tile_size/2 + i;
					int y1 = (tileY+1)*tile_size + tile_size/2 - i;
					for (int j = 0; (j <= i) && (j < 2); j++){
						int yy0 = tileY*tile_size +  tile_size/2 - j;
						int yy1 = tileY*tile_size +  tile_size/2 + j;
						int xx0 = tileX*tile_size +  tile_size/2 - j;
						int xx1 = tileX*tile_size +  tile_size/2 + j;
						rslt[width * yy0 + x0] = stress[0][nTile];
						rslt[width * yy1 + x0] = stress[0][nTile];
						rslt[width * yy0 + x1] = stress[0][nTile];
						rslt[width * yy1 + x1] = stress[0][nTile];

						rslt[width * y0 + xx0] = stress[1][nTile];
						rslt[width * y0 + xx1] = stress[1][nTile];
						rslt[width * y1 + xx0] = stress[1][nTile];
						rslt[width * y1 + xx1] = stress[1][nTile];
					}
				}
			}
		}
		return rslt;
	}
	
	public double [] dbgRescaleToPixels(
			double [] data,
			int       tile_size)
	{
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int width =  tilesX * tile_size;
		int height = tilesY * tile_size;
		double [] rslt = new double [width*height];
		for (int nTile = 0; nTile < data.length; nTile++) {
			int tileY = nTile / tilesX;
			int tileX = nTile % tilesX;
			for (int y = tileY*tile_size; y < (tileY + 1)*tile_size; y++){
				for (int x = tileX*tile_size; x < (tileX + 1)*tile_size; x++){
					rslt[width * y + x] = data[nTile];
				}
			}
		}
		return rslt;
	}

	
	
	
	public void smoothDisparity(
//			final double      break3, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg)
			final double      dispPull, // clt_parameters.tiDispPull or 0.0
			final int         mask,     // 1 - work on internal elements, 2 - on border elements, 3 - both (internal first);
			final int         num_passes,
			final double      maxDiff, // maximal change in any of the disparity values
			final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			final double  []  disparity,          // current disparity value
			final double  []  measured_disparity, // measured disparity
			final double  []  strength,
			final double  []  hor_disparity, // not yet used
			final double  []  hor_strength, // not yet used
			final boolean []  selected,
			final boolean []  border,       // may be null
			final  EyesisCorrectionParameters.CLTParameters  clt_parameters,
			final int         threadsMax,      // maximal number of threads to launch                         
			final int         debugLevel)
	{
		final int dbg_tile = -1; // 28643; // x=131, y=88
		showDoubleFloatArrays sdfa_instance = null;
		if (debugLevel > 0) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final int numThreads = threads.length;
		if (debugLevel > 0) System.out.println("smoothDisparity(): using "+numThreads+" threads");
		final int len = disparity.length;
		int numBorder = 0, numInternal = 0;
		for (int i = 0; i < len; i++){
			if ((border != null) && border[i]) numBorder++;
			else if (selected[i]  && (neighbors[i] >= 0)) numInternal++; // only if not border
		}
		final int numTiles = (((mask & 1) != 0)? numInternal : 0) + (((mask & 2) != 0)? numBorder : 0); 
		final int [] indices =   new int [numBorder + numInternal]; // internal excludes border
		int indx = 0;
		for (int i = 0; i < len; i++){
			if ((border != null) && border[i]){
				if ((mask & 2) != 0) indices[indx++] = i;
			} else if (selected[i]  && (neighbors[i] >= 0)) {
				if ((mask & 1) != 0) indices[indx++] = i;
			}
		}
		final int [] dirs8 =     {-tp.getTilesX(),  -tp.getTilesX() + 1, 1, tp.getTilesX() +1, tp.getTilesX(), tp.getTilesX() - 1, -1, -tp.getTilesX() - 1};
		final double [] rigid8 = {
				clt_parameters.tiRigidVertical,
				clt_parameters.tiRigidDiagonal,
				clt_parameters.tiRigidHorizontal,
				clt_parameters.tiRigidDiagonal,
				clt_parameters.tiRigidVertical,
				clt_parameters.tiRigidDiagonal,
				clt_parameters.tiRigidHorizontal,
				clt_parameters.tiRigidDiagonal};

		final double [] rslt_diffs = new double [numThreads];
		final double [] zero_diffs = new double [numThreads]; // all 0;
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
//		final double [][] disp_data = {disparity, new double [len]};
		final double [][] disp_data = {disparity, disparity.clone()};
		// neighbors
		final double [][] dbg_pull = new double [3][len];
		if (debugLevel > 0) {
			System.out.println("smoothDisparity()");
		}

		for (int pass = 0; (pass < num_passes) || (num_passes ==0); pass++) {
			final int dbg_pass = pass;
			System.arraycopy(zero_diffs, 0, rslt_diffs, 0, numThreads); // set all to 0
			ai.set(0);
			ai_numThread.set(0);
			final int dbg_fpass = pass;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double diff;
						int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
						for (int iTile = ai.getAndIncrement(); iTile < numTiles; iTile = ai.getAndIncrement()) {
							int nTile = indices[iTile];
							int tileY = nTile/tilesX; 
							int tileX = nTile - (tileY * tilesX);
							// calculate pull by neighbors - does not depend on strength of disparity
							// to access from the debugger
							double [] dbg_disp_data = disp_data[0];
							boolean [] dbg_selected = selected;
							boolean [] dbg_border = border;
							int [] dbg_neighbors = neighbors; 
							double neib_avg = 0.0, neib_weight = 0.0;
							if ((debugLevel > 0) &&(nTile == dbg_tile)) {
								System.out.println("smoothDisparity() nTile = "+nTile+" tileX="+tileX+" tileY="+tileY);
							}
							for (int i = 0; i < dirs8.length; i++) if ((neighbors[nTile] >= 0) && ((neighbors[nTile] & (1 << i)) != 0)) {
								double w = rigid8[i]; // no strength here  strength[nTile + dirs8[i]]
								int dbg_dirs8 = nTile + dirs8[i];
								double dbg_d = disp_data[0][nTile + dirs8[i]];
								if ((debugLevel > 0) &&(nTile == dbg_tile)) {
									System.out.println("smoothDisparity(), neib_avg = "+neib_avg+", neib_weight="+neib_weight);
								}
								neib_weight += w;
								neib_avg += w * disp_data[0][nTile + dirs8[i]];
								int dbg_pass = dbg_fpass;
//								if (((debugLevel > 0) &&(nTile == dbg_tile)) ||Double.isNaN(disp_data[0][nTile + dirs8[i]])) {
								if (Double.isNaN(disp_data[0][nTile + dirs8[i]])) {
									System.out.println("nTile="+nTile+": smoothDisparity(),dirs8["+i+"]="+ dirs8[i]+", disp_data[0]["+(nTile + dirs8[i])+"]="+disp_data[0][nTile + dirs8[i]]+", pass="+dbg_pass);
								}
								
//								if (nTile == 28967){//  (tileY == 89) { // (nTile == 28964){
//									System.out.println("smoothDisparity() c1: i = "+i+" w="+w+" dbg_dirs8="+dbg_dirs8+" dbg_d="+dbg_d);
//								}
							}
							double new_disp = disp_data[0][nTile];
							
							if ((debugLevel > 0) &&(nTile == dbg_tile)) {
								System.out.println("smoothDisparity(), neib_avg = "+neib_avg+", neib_weight="+neib_weight);
							}
							
							
							
							if (neib_weight > 0.0){
								neib_avg/= neib_weight;
//								if (nTile == 28967){//  (tileY == 89) { // (nTile == 28964){
//									System.out.println("smoothDisparity() c2: neib_weight = "+neib_weight+" neib_avg="+neib_avg);
//								}
								// calculate pull by the measured disparity
								double disparity_diff = (measured_disparity[nTile] - disp_data[0][nTile]);
								double eff_strength = ((border != null) && border[nTile])? 0.0: (strength[nTile] - clt_parameters.tiStrengthOffset);
								if (eff_strength < 0) eff_strength = 0;
								double disparity_pull = eff_strength;
//								if (tileY == 89){
//									System.out.println("smoothDisparity() d: tileY="+tileY+", tileX="+tileX+" nTile="+nTile+" neib_avg="+neib_avg);
//								}

								if (clt_parameters.tiDispPow > 0.0){
									disparity_pull *= Math.abs(disparity_diff)/clt_parameters.tiDispScale;
									if (clt_parameters.tiDispPow > 1.0){
										disparity_pull = Math.pow(disparity_pull, clt_parameters.tiDispPow);
									}
								}
								//									disparity_pull *= clt_parameters.tiDispPull;
								disparity_pull *= dispPull;
								dbg_pull[0][nTile] = neib_avg -                   disp_data[0][nTile];
								dbg_pull[1][nTile] = (measured_disparity[nTile] -  disp_data[0][nTile])*disparity_pull;
								dbg_pull[2][nTile] = disparity_pull;


								new_disp = (neib_avg * neib_weight +  measured_disparity[nTile] * disparity_pull) / (neib_weight + disparity_pull);
								//									System.out.println("neib_avg = "+neib_avg+" disp_data[0]["+nTile+"]="+disp_data[0][nTile]+" new_disp="+new_disp);

							}
							disp_data[1][nTile] = new_disp; // update with new value that will be used at new iteration

							diff = Math.abs(disp_data[1][nTile] - disp_data[0][nTile]);
							if (diff > rslt_diffs[numThread]) rslt_diffs[numThread] = diff;
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			double [][] dbg_data = {measured_disparity,disp_data[0],disp_data[1],dbg_pull[0],dbg_pull[1],dbg_pull[2],strength};
			String [] titles = {"measured","[0]","[1]","avg","meas","pull","strength"};
			if ((debugLevel> 2) && ((pass ==0) || (pass >= (num_passes-2)))){
				sdfa_instance.showArrays(dbg_data,tilesX, tilesY, true, "disp_smoothed",titles);
			}

			double [] tmp= disp_data[0]; // swap, new data will be in disp_data[0], disp_data[1] to be written to by threads in next run 
			disp_data[0] = disp_data[1];
			disp_data[1] = tmp;
//			if ((debugLevel > 0) && (dbg_tile >= 0)) {
//				System.out.println("smoothDisparity() pass="+pass+", disp_data[0]["+dbg_tile+"] = "+disp_data[0][dbg_tile]+" disp_data[1]["+dbg_tile+"] = "+disp_data[1][dbg_tile]);
//			}

			if (maxDiff > 0){
				double diff = 0.0;
				for (int i = 0; (i < numThreads) && (diff <= maxDiff) ; i++) if (diff < rslt_diffs[i]) diff = rslt_diffs[i];
				if (diff <= maxDiff) {
					if (debugLevel > 0) System.out.println("smoothDisparity(): pass = "+pass+", diff = "+diff+" <= "+maxDiff);
					break;
				}
			}
		}
		// dbgDeriv
		System.arraycopy(disp_data[0], 0, disp_data[1], 0, len); // set all to 0 (odd/even, disparity may be now any of 2)

		//		double [][] dbg_data = {measured_disparity,disp_data[0],disp_data[1],strength, disparity, dbgDeriv[0], dbgDeriv[1]};
		//		String [] titles = {"measured", "[0]", "[1]", "strength", "disp", "deriv0", "deriv1"};
		//		sdfa_instance.showArrays(dbg_data,tilesX, tilesY, true, "disp_smoothed", titles);
	}

	/**
	 * Select tiles where measured disparity is closer/farther than currently approximated by a sufficient margin. Return null
	 * if none are
	 */
	public boolean [] findNearFar(
			boolean findFar,
			double  threshold, // select tiles with non-zero strength that are far/near
			final double  []  disparity,          // current (approximated) disparity value
			final double  []  measured_disparity, // measured disparity
			final double  []  strength,           // masked by previously removed tiles as far/near
			final boolean []  selected)
	{
		int leng = disparity.length;
		boolean [] found = new boolean [leng];
		int numFound = 0;
		for (int i = 0; i<leng; i++)
			if ((strength[i] > 0.0) && ((selected == null) || selected[i]) &&
					((findFar? (disparity[i] - measured_disparity[i]) : (measured_disparity[i]-disparity[i])) > threshold)){
				found[i] = true;
				numFound++;
		} else {
			found[i] = false;
		}
		return (numFound > 0)? found : null;		
	}
	
	public void breakDisparity( // break using 3-th derivative
			final double      break3, // clt_parameters.tiBreak/0 allow disconnecting from neighbors (fg/bg). if 0.0 - do not break, just calculate stresses
			final int         mode, // 0: 3-rd derivative, 1: - third * first (compare positive threshold only), 2: second by abs(first) (compare positive threshold only) 
			final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			final double  []  disparity,          // current disparity value
			final boolean []  selected,
			final boolean     extend_flat, // if the tile is on the hor/vert edge, assume same disparity on the other side
			final double      k_same,
			final double      k_turn,
			final  EyesisCorrectionParameters.CLTParameters  clt_parameters,
			final double [][] dbgDeriv, //double [2][len] or null;
			final int         threadsMax,      // maximal number of threads to launch                         
			final int         debugLevel)
	{
//		showDoubleFloatArrays sdfa_instance = null;
//		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		final int dbg_tile = (debugLevel > -1)  ? 27991 : -1; // x = 127, y = 86
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final int numThreads = threads.length;
		if (debugLevel > -1) System.out.println("breakDisparity3(): using "+numThreads+" threads");
		final int len = disparity.length;
		int numTiles = 0;
		for (int i = 0; i < len; i++){
			if (selected[i] && (neighbors[i] >= 0)) numTiles++; // only if not border
		}
		final int [] indices =   new int [numTiles]; // internal excludes border
		int indx = 0;
		for (int i = 0; i < len; i++){
			if (selected[i] && (neighbors[i] >= 0)) {
				indices[indx++] = i;
			}
		}
		final AtomicInteger ai = new AtomicInteger(0);
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final double [][] deriv3HV = new double [2][len]; // 4-th derivative [0] - horizontal, 1 - vertical 
		final double [] k4_3 =  {-1.0,  3.0, -3.0,  1.0}; // 3-rd
		final double [] k4_2 =  {-1.0,  1.0, 1.0, -1.0};  // - 2-nd (average for 2 sides)
		final double [] k4 = (mode == 2)? k4_2 : k4_3;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iTile = ai.getAndIncrement(); iTile < indices.length; iTile = ai.getAndIncrement()) {
						int nTile = indices[iTile];
//						int tileY = nTile/tilesX; 
//						int tileX = nTile - (tileY * tilesX);
//						double deriv4;
						int neib = neighbors[nTile];
						if (neib < 0) neib = 0; // should not get here 
						// try break horizontally (always - to the right)
						int b = (1 << 2); // "E" (right)
						int rb = (1 << 6); // "W" (left)
						if (nTile == dbg_tile) {
							System.out.println("breakDisparity{}: nTile = "+nTile);
						}
						if ((neib & b) != 0){
							hor_label:
							{
//							double deriv = k_inner * (disparity[nTile + 1]- disparity[nTile]); // 3 - for 3-rd derivative, 1 - for second (won't work with extension)
							double deriv = k4[2] * disparity[nTile + 1] +k4[1]* disparity[nTile]; // 3 - for 3-rd derivative, 1 - for second (won't work with extension)
							if ((neib & rb) != 0) deriv += k4[0]*disparity[nTile - 1];
							else if (extend_flat) deriv += k4[0]*disparity[nTile];
							else break hor_label;
							if ((neighbors[nTile + 1] >=0) && ((neighbors[nTile + 1] & b) != 0)) deriv += k4[3]* disparity[nTile + 2];
							else if (extend_flat)                                                deriv += k4[3]* disparity[nTile + 1];
							else break hor_label;
							switch (mode){
							case 1:
								deriv *= disparity[nTile]- disparity[nTile + 1]; //  3-rd derivative * 1-st (reversed order to get -1)
								break;
							case 2:
								deriv *= Math.abs(disparity[nTile + 1]- disparity[nTile]); // 2-nd derivative * abs(1-st derivative)
							default: // do nothing
							}
							double disp_avg = 0.5*(disparity[nTile + 1] + disparity[nTile]);
							if ((clt_parameters.tiBreakNorm > 0.0) && (disp_avg > clt_parameters.tiBreakNorm )){
								deriv *= clt_parameters.tiBreakNorm / disp_avg;
							}
							deriv3HV[0][nTile] = deriv;
							}
						}
						b = (1 << 4); // "S" (down)
						rb = (1 << 0); // "N" (up)
						// try break vertical (always - down)
						if ((neib & b) != 0){
							vert_label:
							{
//							double deriv = k_inner * (disparity[nTile + tilesX]- disparity[nTile]);
							double deriv = k4[2] * disparity[nTile + tilesX] +k4[1] * disparity[nTile]; // 3 - for 3-rd derivative, 1 - for second (won't work with extension)
							if ((neib & rb) != 0) deriv += k4[0]*disparity[nTile - tilesX];
							else if (extend_flat) deriv += k4[0]*disparity[nTile];
							else break vert_label;
							if ((neighbors[nTile +  + tilesX] >=0) && ((neighbors[nTile + tilesX] & b) != 0)) deriv += k4[3] * disparity[nTile + 2 * tilesX];
							else if (extend_flat)                                                             deriv += k4[3] * disparity[nTile + tilesX];
							else break vert_label;
							switch (mode){
							case 1:
								deriv *= disparity[nTile]- disparity[nTile + tilesX]; //  3-rd derivative * - 1-st (reversed order to get -1)
								break;
							case 2:
								deriv *= Math.abs(disparity[nTile + tilesX]- disparity[nTile]); // 2-nd derivative * abs(1-st derivative)
							default: // do nothing
							}
							double disp_avg = 0.5*(disparity[nTile + tilesX] + disparity[nTile]);
							if ((clt_parameters.tiBreakNorm > 0.0) && (disp_avg > clt_parameters.tiBreakNorm )){
								deriv *= clt_parameters.tiBreakNorm / disp_avg;
							}
							
							deriv3HV[1][nTile] = deriv; // vertical 4-th derivative
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		// Amplify same-direction (and turn) stress
		if ((k_same != 0.0) || (k_turn != 0.0)){
			final double [][] deriv3HVCopy = {deriv3HV[0].clone(),deriv3HV[1].clone()};
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int iTile = ai.getAndIncrement(); iTile < indices.length; iTile = ai.getAndIncrement()) {
							int nTile = indices[iTile];
							int tileY = nTile/tilesX; 
							int tileX = nTile - (tileY * tilesX);
							if (nTile == dbg_tile) {
								System.out.println("breakDisparity{}: nTile = "+nTile);
							}
							if (k_same !=0 ) {
								// vertical "gradient" in the same horizontal row
								if ((tileX > 0) &&           selected[nTile - 1])      deriv3HV[1][nTile] += k_same * deriv3HVCopy[1][nTile - 1];
								if ((tileX < (tilesX -1)) && selected[nTile + 1])      deriv3HV[1][nTile] += k_same * deriv3HVCopy[1][nTile + 1];
								
								// horizontal  "gradient" in the same vertical column
								if ((tileY > 0) &&           selected[nTile - tilesX]) deriv3HV[0][nTile] += k_same * deriv3HVCopy[0][nTile - tilesX];
								if ((tileY < (tilesY -1)) && selected[nTile + tilesX]) deriv3HV[0][nTile] += k_same * deriv3HVCopy[0][nTile + tilesX];
							}
							if (k_turn !=0 ) {
								if ((tileY > 0) &&                                    selected[nTile - tilesX])     deriv3HV[0][nTile] += k_turn * deriv3HVCopy[1][nTile - tilesX    ];
								if ((tileY > 0) && (tileX < (tilesX -1)) &&           selected[nTile - tilesX + 1]) deriv3HV[0][nTile] += k_turn * deriv3HVCopy[1][nTile - tilesX + 1];
								if (true)                                                                           deriv3HV[0][nTile] += k_turn * deriv3HVCopy[1][nTile];
								if (               (tileX < (tilesX -1)) &&           selected[nTile          + 1]) deriv3HV[0][nTile] += k_turn * deriv3HVCopy[1][nTile           +1];
								
								if (                         (tileX > 0) &&           selected[nTile          - 1]) deriv3HV[1][nTile] += k_turn * deriv3HVCopy[0][nTile          - 1];
								if ((tileY < (tilesY -1)) && (tileX > 0) &&           selected[nTile + tilesX - 1]) deriv3HV[1][nTile] += k_turn * deriv3HVCopy[0][nTile + tilesX - 1];
								if ((tileY < (tilesY -1)) &&                          selected[nTile + tilesX    ]) deriv3HV[1][nTile] += k_turn * deriv3HVCopy[0][nTile + tilesX   ];
								if (true)                                                                           deriv3HV[1][nTile] += k_turn * deriv3HVCopy[0][nTile];
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		// compare to threshold and update neighbors (does twice, but multi-threading compatible)
		if (break3 >0.0) {
			final int [] neighbors_in = neighbors.clone();
			ai.set(0);
			if (mode == 0) {
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int iTile = ai.getAndIncrement(); iTile < indices.length; iTile = ai.getAndIncrement()) {
								int nTile = indices[iTile];
								int tileY = nTile/tilesX; 
								int tileX = nTile - (tileY * tilesX);
								int neib = neighbors_in[nTile];
								if (neib < 0) neib = 0;
								if (Math.abs(deriv3HV[0][nTile]) >= break3) neib &= ~0b00001110; // break all right side
								if (Math.abs(deriv3HV[1][nTile]) >= break3) neib &= ~0b00111000; // break all down  side
								if ((tileX > 0) && selected[nTile - 1] &&      (Math.abs(deriv3HV[0][nTile - 1]) >= break3))      neib &= ~0b11100000; // break all left side
								if ((tileY > 0) && selected[nTile - tilesX] && (Math.abs(deriv3HV[1][nTile - tilesX]) >= break3)) neib &= ~0b10000011; // break all up side
								// corners
								if ((tileY > 0) && (tileX <(tilesX-1)) && selected[nTile - tilesX +1] && (Math.abs(deriv3HV[1][nTile - tilesX + 1]) >= break3)) neib &= ~0b00000010; // break NE
								if ((tileY > 0) &&                        selected[nTile - tilesX   ] && (Math.abs(deriv3HV[0][nTile - tilesX    ]) >= break3)) neib &= ~0b00000010; // break NE

								if (               (tileX <(tilesX-1)) && selected[nTile          +1] && (Math.abs(deriv3HV[1][nTile          + 1]) >= break3)) neib &= ~0b00001000; // break SE
								if ((tileY <(tilesY-1))                && selected[nTile+ tilesX    ] && (Math.abs(deriv3HV[0][nTile + tilesX    ]) >= break3)) neib &= ~0b00001000; // break SE

								if ((tileX > 0)                        && selected[nTile          -1] && (Math.abs(deriv3HV[1][nTile          - 1]) >= break3)) neib &= ~0b00100000; // break SW
								if ((tileX > 0) && (tileY <(tilesY-1)) && selected[nTile + tilesX -1] && (Math.abs(deriv3HV[0][nTile + tilesX - 1]) >= break3)) neib &= ~0b00100000; // break SW

								if ((tileX > 0) && (tileY > 0)         && selected[nTile - tilesX -1] && (Math.abs(deriv3HV[1][nTile - tilesX - 1]) >= break3)) neib &= ~0b10000000; // break NW
								if ((tileX > 0) && (tileY > 0)         && selected[nTile - tilesX -1] && (Math.abs(deriv3HV[0][nTile - tilesX - 1]) >= break3)) neib &= ~0b10000000; // break NW

								neighbors[nTile] = neib;
							}
						}
					};
				}
			} else { // modes 1 & 2 - use positive values comparison, not Math.abs()
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						public void run() {
							for (int iTile = ai.getAndIncrement(); iTile < indices.length; iTile = ai.getAndIncrement()) {
								int nTile = indices[iTile];
								int tileY = nTile/tilesX; 
								int tileX = nTile - (tileY * tilesX);
								int neib = neighbors_in[nTile];
								if (neib < 0) neib = 0;
								if (deriv3HV[0][nTile] >= break3) neib &= ~0b00001110; // break all right side
								if (deriv3HV[1][nTile] >= break3) neib &= ~0b00111000; // break all down  side
								if ((tileX > 0) && selected[nTile - 1] &&      (deriv3HV[0][nTile - 1] >= break3))      neib &= ~0b11100000; // break all left side
								if ((tileY > 0) && selected[nTile - tilesX] && (deriv3HV[1][nTile - tilesX] >= break3)) neib &= ~0b10000011; // break all up side
								// corners
								if ((tileY > 0) && (tileX <(tilesX-1)) && selected[nTile - tilesX +1] && (deriv3HV[1][nTile - tilesX + 1] >= break3)) neib &= ~0b00000010; // break NE
								if ((tileY > 0) &&                        selected[nTile - tilesX   ] && (deriv3HV[0][nTile - tilesX    ] >= break3)) neib &= ~0b00000010; // break NE

								if (               (tileX <(tilesX-1)) && selected[nTile          +1] && (deriv3HV[1][nTile          + 1] >= break3)) neib &= ~0b00001000; // break SE
								if ((tileY <(tilesY-1))                && selected[nTile+ tilesX    ] && (deriv3HV[0][nTile + tilesX    ] >= break3)) neib &= ~0b00001000; // break SE

								if ((tileX > 0)                        && selected[nTile          -1] && (deriv3HV[1][nTile          - 1] >= break3)) neib &= ~0b00100000; // break SW
								if ((tileX > 0) && (tileY <(tilesY-1)) && selected[nTile + tilesX -1] && (deriv3HV[0][nTile + tilesX - 1] >= break3)) neib &= ~0b00100000; // break SW

								if ((tileX > 0) && (tileY > 0)         && selected[nTile - tilesX -1] && (deriv3HV[1][nTile - tilesX - 1] >= break3)) neib &= ~0b10000000; // break NW
								if ((tileX > 0) && (tileY > 0)         && selected[nTile - tilesX -1] && (deriv3HV[0][nTile - tilesX - 1] >= break3)) neib &= ~0b10000000; // break NW

								neighbors[nTile] = neib;
							}
						}
					};
				}

			}

			ImageDtt.startAndJoin(threads);
		}
		// return derivatives for debug purposes if array is provided (not null)

		if (dbgDeriv != null){
			dbgDeriv[0] = deriv3HV[0];
			dbgDeriv[1] = deriv3HV[1];
		}
	}

	
	
	public void reconnectDisparity( // connect disconnected tiles if they have close approximated disparity
			final double      maxDiffOrto,
			final double      maxDiffDiagonal,
			final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			final double  []  disparity,          // current disparity value
			final boolean []  selected,
			final int         threadsMax,      // maximal number of threads to launch                         
			final int         debugLevel)
	{
//		showDoubleFloatArrays sdfa_instance = null;
//		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final int numThreads = threads.length;
		if (debugLevel > -1) System.out.println("reconnectDisparity(): using "+numThreads+" threads");
		final int len = disparity.length;
		int numTiles = 0;
		for (int i = 0; i < len; i++){
			if (selected[i] && (neighbors[i] >= 0)) numTiles++; // only if not border
		}
		final int [] indices =   new int [numTiles]; // internal excludes border
		int indx = 0;
		for (int i = 0; i < len; i++){
			if (selected[i] && (neighbors[i] >= 0)) {
				indices[indx++] = i;
			}
		}
		final AtomicInteger ai = new AtomicInteger(0);
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iTile = ai.getAndIncrement(); iTile < indices.length; iTile = ai.getAndIncrement()) {
						int nTile = indices[iTile];
						int tileY = nTile/tilesX; 
						int tileX = nTile - (tileY * tilesX);
						int neib = neighbors[nTile]; 
						if (neib < 0) neib = 0;
						int b = (1 << 2); // "E" (right)
						if (((neib & b) == 0) && (tileX < (tilesX - 1) && selected[nTile + 1])) {
							if (Math.abs(disparity[nTile] - disparity[nTile + 1]) < maxDiffOrto){
								neib |= b;
							}
						}
						b = (1 << 4); // "S" (down)
						if (((neib & b) == 0) && (tileY < (tilesY - 1)) && selected[nTile + tilesX]) {
							if (Math.abs(disparity[nTile] - disparity[nTile + tilesX]) < maxDiffOrto){
								neib |= b;
							}
						}
						b = (1 << 6); // "W" (left)
						if (((neib & b) == 0) && (tileX > 0) && selected[nTile - 1]) {
							if (Math.abs(disparity[nTile] - disparity[nTile - 1]) < maxDiffOrto){
								neib |= b;
							}
						}
						b = (1 << 0); // "N" (up)
						if (((neib & b) == 0) && (tileY > 0) && selected[nTile - tilesX]) { // was  not connected
							if (Math.abs(disparity[nTile] - disparity[nTile - tilesX]) < maxDiffOrto){
								neib |= b;
							}
						}
						
						b = (1 << 1); // "NE"
						if (((neib & b) == 0) && (tileX < (tilesX - 1)) && (tileY > 0) && selected[nTile - tilesX + 1]) {
							if (Math.abs(disparity[nTile] - disparity[nTile - tilesX + 1]) < maxDiffDiagonal){
								neib |= b;
							}
						}
						b = (1 << 3); // "SE"
						if (((neib & b) == 0) && (tileX < (tilesX - 1)) && (tileY < (tilesY - 1)) && selected[nTile + tilesX + 1]) {
							if (Math.abs(disparity[nTile] - disparity[nTile + tilesX + 1]) < maxDiffDiagonal){
								neib |= b;
							}
						}
						b = (1 << 5); // "SW"
						if (((neib & b) == 0) && (tileX > 0) && (tileY < (tilesY - 1)) && selected[nTile + tilesX - 1]) {
							if (Math.abs(disparity[nTile] - disparity[nTile + tilesX - 1]) < maxDiffDiagonal){
								neib |= b;
							}
						}
						b = (1 << 7); // "NW"
						if (((neib & b) == 0) && (tileX > 0) && (tileY > 0) && selected[nTile - tilesX - 1]) {
							if (Math.abs(disparity[nTile] - disparity[nTile - tilesX - 1]) < maxDiffDiagonal){
								neib |= b;
							}
						}
						neighbors[nTile] = neib;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		reconnectDiagonals( // connect diagonals if 2 sides are connected
				neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				selected,
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		
	}

	public void reconnectDiagonals( // connect diagonals if 2 sides are connected
			final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			final boolean []  selected,
			final int         threadsMax,      // maximal number of threads to launch                         
			final int         debugLevel)
	{
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final int numThreads = threads.length;
		if (debugLevel > -1) System.out.println("reconnectDiagonals(): using "+numThreads+" threads");
		final int len = neighbors.length;
		int numTiles = 0;
		for (int i = 0; i < len; i++){
			if (selected[i] && (neighbors[i] >= 0)) numTiles++; // only if not border
		}
		final int [] indices =   new int [numTiles]; // internal excludes border
		int indx = 0;
		for (int i = 0; i < len; i++){
			if (selected[i] && (neighbors[i] >= 0)) {
				indices[indx++] = i;
			}
		}
		final AtomicInteger ai = new AtomicInteger(0);
		final int tilesX = tp.getTilesX();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int iTile = ai.getAndIncrement(); iTile < indices.length; iTile = ai.getAndIncrement()) {
						int nTile = indices[iTile];
						int neib = neighbors[nTile];
						if (neib < 0) neib = 0;
						if (
								((neib & (1 << 0)) != 0) &&
								((neighbors[nTile - tilesX] & (1 << 2)) != 0)){
							neib |= (1 <<1);
						} else if (
								((neib & (1 << 2)) != 0) &&
								((neighbors[nTile + 1] &      (1 << 0)) != 0)){
							neib |= (1 <<1);
						}
						
						if (
								((neib & (1 << 4)) != 0) &&
								((neighbors[nTile + tilesX] & (1 << 2)) != 0)){
							neib |= (1 << 3);
						} else if (
								((neib & (1 << 2)) != 0) &&
								((neighbors[nTile + 1] &      (1 << 4)) != 0)){
							neib |= (1 << 3);
						}

						if (
								((neib & (1 << 4)) != 0) &&
								((neighbors[nTile + tilesX] & (1 << 6)) != 0)){
							neib |= (1 << 5);
						} else if (
								((neib & (1 << 6)) != 0) &&
								((neighbors[nTile - 1] &      (1 << 4)) != 0)){
							neib |= (1 << 5);
						}
						if (
								((neib & (1 << 0)) != 0) &&
								((neighbors[nTile - tilesX] & (1 << 6)) != 0)){
							neib |= (1 << 7);
						} else if (
								((neib & (1 << 6)) != 0) &&
								((neighbors[nTile - 1] &      (1 << 0)) != 0)){
							neib |= (1 << 7);
						}
						neighbors[nTile] = neib;
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	/**
	 * "Heals" small gaps in connections between the members of the same cluster, if the missing connections
	 * All start/ end in the same cluster (so no borders with others), if the total length does not exceed the
	 * threshold. 
	 * @param neighbors array of connections of the current tile with 8 others - should be consistent
	 * @param maxlen maximal length of the break to heal 
	 * @param selected   selected tiles (only used for filling diagonal connections)
	 * @param threadsMax maximal number of CPU therads to use
	 * @param debugLevel debug level
	 * @return number of the new orthogonal connections
	 */
	
	public int healSame( // returns number of new ortho connections
			int [] neighbors,
			int maxlen,
			// just to fill in diagonals
			final boolean []  selected,        
			final int         threadsMax,                         
			final int         debugLevel)
	{
//		int leng = neighbors.length;
		int NEIB_RIGHT = 4;
		int NEIB_LEFT = 64;
		int NEIB_UP =    1;
		int NEIB_DOWN = 16;
		int FLD_PROHIB = -2;
		int FLD_CONN =   -1;
		int FLD_EMPTY =   0;
		int FLD_WALKED =  1;
		int [][][] TRANSITIONS = {
				{// direction was 0 (right, going up
					{ 1, -1, 1},  // x += 1; y -= 1; newdir = 1 (down, going right)  first  (rightmost choice)
					{ 0, -1, 0},  // x += 0; y -= 1; newdir = 0 (right, going up)    second (straight choice) 
					{ 0, -1, 3},  // x += 0; y -= 1; newdir = 3 (down, going left)   third  (left choice)      
					{ 0,  0, 2}}, // x += 0; y += 0; newdir = 2 (right, going down)  last   (back choice)     
				{// direction was 1 (down, going right)
					{ 0,  1, 2},  // x += 0; y += 1; newdir = 2 (right, going down)
					{ 1,  0, 1},  // x += 1; y += 0; newdir = 1 (down, going right)
					{ 0,  0, 0},  // x += 0; y += 0; newdir = 0 (right, going up)
					{ 0,  0, 3}}, // x += 0; y += 0; newdir = 3 (down, going left)
				{// direction was 2 (right, going down)
					{ 0,  0, 3},  // x += 0; y += 0; newdir = 3 (down, going left)
					{ 0,  1, 2},  // x += 0; y += 1; newdir = 2 (right, going down)
					{ 1,  0, 1},  // x += 1; y += 0; newdir = 1 (down, going right)
					{ 0,  0, 0}}, // x += 0; y += 0; newdir = 0 (right, going up)
				{// direction was 3 (down, going left)
					{ -1, 0, 0},  // x -= 1; y += 0; newdir = 0 (right, going up)
					{ -1, 0, 3},  // x -= 1; y += 0; newdir = 3 (down, going left)
					{ -1, 1, 2},  // x -= 1; y += 1; newdir = 2 (right, going down)
					{  0, 0, 1}}};// x += 0; y += 0; newdir = 1 (down, going right)
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		class XYDir {
			int x;
			int y;
			int dir;
			XYDir (int x, int y, int dir){
				this.x = x; this.y = y; this.dir = dir;
			}
		}
		ArrayList<XYDir> xYDirList = new ArrayList<XYDir>();
		// adding top row and left column
		int [][][] field = new int [tilesY + 1][tilesX + 1][4]; //  right, going up; down, going right; right going down; down, going left;
		for (int i = 0; i < tilesY + 1; i++){
			for (int j = 0; j < tilesX + 1; j++){
				for (int k = 0; k < 4; k++){
					field[i][j][k] = FLD_PROHIB;
				}
			}
		}
		for (int ty = 0; ty < tilesY; ty++){
			for (int tx = 0; tx < tilesX; tx++){
				int neib = neighbors[ty * tilesX + tx];
				if (neib < 0) neib = 0;
				int []cell = field[ty + 1][tx + 1];
				if ( neib!= 0) {
					cell[0] = ((neib & NEIB_RIGHT) == 0)? FLD_EMPTY:FLD_CONN;
					cell[2] = ((neib & NEIB_RIGHT) == 0)? FLD_EMPTY:FLD_CONN;
					cell[1] = ((neib & NEIB_DOWN) == 0)?  FLD_EMPTY:FLD_CONN;
					cell[3] = ((neib & NEIB_DOWN) == 0)?  FLD_EMPTY:FLD_CONN;
				}
			}				
		}
		int updated = 0;
		for (int startY = 1; startY < tilesY; startY++){
			for (int startX = 1; startX < tilesX; startX++){				
				for (int startDir = 0; startDir < 2; startDir++){
					if (field[startY][startX][startDir] == FLD_EMPTY){
						boolean prohib = false; // met at least one prohibited cell
						int tx =   startX;
						int ty =   startY;
						int tDir = startDir;
						
						int gapLen = 0; // number of cells traversed
//						boolean fwd = true; // up for right, right for down
						xYDirList.clear();
						gapLen++;
						field[ty][tx][tDir] = FLD_WALKED;
						xYDirList.add(new XYDir(tx,ty,tDir));
						// walk right hand, mark path with FLD_WALKED, look for prohibited. Verify each cell walked twice
						// when done, mark with either prohibited (too long or met prohibited, or not twice) or as connected (and update neib)
						int dirChoice = -1;
						walking:{
							while (true) {
								for (dirChoice =0; dirChoice < 4; dirChoice++){
									int tx1 = tx + TRANSITIONS[tDir][dirChoice][0];
									int ty1 = ty + TRANSITIONS[tDir][dirChoice][1];
									int tDir1 =    TRANSITIONS[tDir][dirChoice][2];
									// Was already here (in the same direction?
									if ((tx1 < 0) || (ty1 < 0) || (tx1 > tilesX) || (ty1 > tilesY)){
										prohib = true;	
									} else if (field[ty1][tx1][tDir1] == FLD_WALKED) { // 325
										break walking;
									} else if (field[ty1][tx1][tDir1] == FLD_EMPTY) {
										tx = tx1;
										ty = ty1;
										tDir = tDir1;
										gapLen++;
										field[ty][tx][tDir] = FLD_WALKED;
										xYDirList.add(new XYDir(tx,ty,tDir));
										break; // just for loop
									} else if (field[ty1][tx1][tDir1] == FLD_PROHIB) {
										prohib = true;
									}
								}
								// for loop (dirChoice) ended . Nothing is needed here, just debug if endless loop 
							}
						}
						// finished walking. Was it good? TODO: more tests are needed - each has forward+back
						boolean islands = false; // islands or outer border
						for (XYDir xYdir: xYDirList) {
							if (field[xYdir.y][xYdir.x][xYdir.dir ^ 2] != FLD_WALKED) { // "^2" - opposite direction
								islands = true;
								break;
							}
						}						
						if (!prohib && !islands && (gapLen <= 2 * maxlen)) { // add connections to neighbors
							for (XYDir xYdir: xYDirList) {
								switch (xYdir.dir){
								case 0: // (right, going up)
									neighbors[(xYdir.y -1) * tilesX + (xYdir.x -1)] |= NEIB_RIGHT;
									// update right neighbor too
									neighbors[(xYdir.y -1) * tilesX + (xYdir.x)] |=    NEIB_LEFT;
									updated++;
									break;
								case 1: // (down, going right)
									neighbors[(xYdir.y -1) * tilesX + (xYdir.x -1)] |= NEIB_DOWN;
									// update down neighbor too
									neighbors[(xYdir.y   ) * tilesX + (xYdir.x -1)] |= NEIB_UP;
									updated++;
									break;
								// do nothing for other (back) cases	
								}
								field[xYdir.y][xYdir.x][xYdir.dir] = FLD_CONN; // mark as connection
								//								field[xYdir.y][xYdir.x][xYdir.dir] = FLD_PROHIB;
							}
						} else { // mark all the path as prohibited
							for (XYDir xYdir: xYDirList) {
								field[xYdir.y][xYdir.x][xYdir.dir] = FLD_PROHIB;
							}
						}
					}
				}
			}
		}
		if (updated > 0) {
		reconnectDiagonals( // connect diagonals if 2 sides are connected
				neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
				selected,
				threadsMax,      // maximal number of threads to launch                         
				debugLevel);
		}
		return updated;
	}
	
	
	/**
	*  For each selected neighbors tile create array[8] of the shared incoming directions
	*  each element index is an incoming direction, value is a bitmask of this and other directions sharing the same tile.
	*  Only non-connected directions from selected tiles are covered (for connected (nighbors[i] +4) % 8 is non-zero)
	*  Not only selected, but one away from selected (with border)
	*  0 - coming in N direction (from S), 1 - coming in NE (from SW), etc. 
	*/  
	public int [][] createOverlapGeometry(
			final int     []  neighbors, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			final boolean []  selected, // only inner?
			final boolean []  border,
			final int         threadsMax,      // maximal number of threads to launch                         
			final int         debugLevel)
	{
		final int debugTile = 41631; 
//		showDoubleFloatArrays sdfa_instance = null;
//		if (debugLevel > -1) sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final int numThreads = threads.length;
		if (debugLevel > -1) System.out.println("createOverlapGeometry(): using "+numThreads+" threads");
		final int len = selected.length;
		int numTiles = 0;
		int numTilesAll = 0;
//		final boolean [] selectedAll = new boolean[len];
		for (int i = 0; i < len; i++){
			if ((selected[i] && (neighbors[i] >= 0)) || border[i]) {
//				selectedAll[i] = true;
				numTilesAll++;
				if (!border[i])	numTiles++;
			}
		}
		final int [] indices =      new int [numTiles];
		final int [] indicesAll =   new int [numTilesAll];
		int indx = 0;
		int indxAll = 0;
		for (int i = 0; i < len; i++){
			if ((selected[i] && (neighbors[i] >= 0)) || border[i]) {
				indicesAll[indxAll++] = i;
				if (!border[i])	indices[indx++] = i;
			}
		}
		final AtomicInteger ai = new AtomicInteger(0);
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int [][] geom = new int [len][];
		

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					int [] b;
					for (int iTile = ai.getAndIncrement(); iTile < indicesAll.length; iTile = ai.getAndIncrement()) {
						int nTile = indicesAll[iTile];
						boolean debug = (debugLevel > -1) && (nTile == debugTile);
						if (debug){
							System.out.println("createOverlapGeometry(): nTile="+nTile);
						}
	
						geom[nTile] = new int [8];
						int [] dbg_geom = geom[nTile];
//						int [] dbg_neib = {neighbors[nTile - tilesX],neighbors[nTile + 1] , neighbors[nTile + tilesX],neighbors[nTile - 1]}; // out of bounds
						int neib = neighbors[nTile];
						if (neib != 0b11111111){ // do nothing for internal tiles
							int tileY = nTile/tilesX; 
							int tileX = nTile - (tileY * tilesX);
							if ((tileY > 0 ) &&                                      selected[nTile - tilesX]     && ((neib & (1 << 0)) == 0)) geom[nTile][4] = 1 << 4;
							if ((tileY > 0 ) &&           (tileX < (tilesX - 1))  && selected[nTile - tilesX + 1] && ((neib & (1 << 1)) == 0)) geom[nTile][5] = 1 << 5;
							if (                          (tileX < (tilesX - 1))  && selected[nTile          + 1] && ((neib & (1 << 2)) == 0)) geom[nTile][6] = 1 << 6;
							if ((tileY < (tilesY - 1)) && (tileX < (tilesX - 1))  && selected[nTile + tilesX + 1] && ((neib & (1 << 3)) == 0)) geom[nTile][7] = 1 << 7;
							if ((tileY < (tilesY - 1)) &&                            selected[nTile + tilesX] &&     ((neib & (1 << 4)) == 0)) geom[nTile][0] = 1 << 0;
							if ((tileY < (tilesY - 1)) && (tileX > 0) &&             selected[nTile + tilesX - 1] && ((neib & (1 << 5)) == 0)) geom[nTile][1] = 1 << 1;
							if (                          (tileX > 0) &&             selected[nTile          - 1] && ((neib & (1 << 6)) == 0)) geom[nTile][2] = 1 << 2;
							if ((tileY > 0) &&            (tileX > 0) &&             selected[nTile - tilesX - 1] && ((neib & (1 << 7)) == 0)) geom[nTile][3] = 1 << 3;
							
							if (tileY > 0) {
								if (((neighbors[nTile - tilesX] &  4) != 0) && ((neighbors[nTile] & 1) == 0)  && ((neighbors[nTile] &   2) == 0)) { // E from N
									geom[nTile][4] |= geom[nTile][5];
									geom[nTile][5]  = geom[nTile][4];
								}
								if (((neighbors[nTile - tilesX] & 64) != 0) && ((neighbors[nTile] & 1) == 0)  && ((neighbors[nTile] & 128) == 0)) { // W from N
									geom[nTile][4] |= geom[nTile][3];
									geom[nTile][3]  = geom[nTile][4];
									if (((neighbors[nTile - tilesX] & 8) != 0) && ((neighbors[nTile] &   4) == 0)) { // SE from N
										geom[nTile][4] |= geom[nTile][6];
										geom[nTile][6]  = geom[nTile][4];
									}
								}
							}
							if (tileX < (tilesX-1)) {
								if (((neighbors[nTile + 1] &  1) != 0) && ((neighbors[nTile] &  4) == 0)  && ((neighbors[nTile] &  2) == 0)) { // N from E
									geom[nTile][6] |= geom[nTile][5];
									geom[nTile][5]  = geom[nTile][6];
								}
								if (((neighbors[nTile + 1] & 16) != 0) && ((neighbors[nTile] &  4) == 0)  && ((neighbors[nTile] &  8) == 0)){ // S from E
									geom[nTile][6] |= geom[nTile][7];
									geom[nTile][7]  = geom[nTile][6];
									if (((neighbors[nTile - tilesX] & 32) != 0) && ((neighbors[nTile] & 16) == 0)) { // SW from E
										geom[nTile][6] |= geom[nTile][0];
										geom[nTile][0]  = geom[nTile][6];
									}
								}
							}
							if (tileY < (tilesY-1)) {
								if (((neighbors[nTile + tilesX] &  4) != 0) &&  ((neighbors[nTile] & 16) == 0)  && ((neighbors[nTile] &  8) == 0)) { // E from S
									geom[nTile][0] |= geom[nTile][7];
									geom[nTile][7]  = geom[nTile][0];
								}
								if (((neighbors[nTile + tilesX] & 64) != 0) &&  ((neighbors[nTile] & 16) == 0)  && ((neighbors[nTile] & 32) == 0)) { // W from S
									geom[nTile][0] |= geom[nTile][1];
									geom[nTile][1]  = geom[nTile][0];
									if (((neighbors[nTile - tilesX] & 128) != 0) && ((neighbors[nTile] & 64) == 0)) { // NW from S
										geom[nTile][0] |= geom[nTile][2];
										geom[nTile][2]  = geom[nTile][0];
									}
								}
							}
							if (tileX > 0) {
								if (((neighbors[nTile - 1] & 16) != 0) && ((neighbors[nTile] & 64) == 0)  && ((neighbors[nTile] & 32) == 0)) { // S from W
									geom[nTile][2] |= geom[nTile][1];
									geom[nTile][1]  = geom[nTile][2];
								}
								if (((neighbors[nTile - 1] & 1) != 0) &&  ((neighbors[nTile] & 64) == 0)  && ((neighbors[nTile] & 128) == 0)) { // N from W
									geom[nTile][2] |= geom[nTile][3];
									geom[nTile][3]  = geom[nTile][2];
									if (((neighbors[nTile - tilesX] & 2) != 0) && ((neighbors[nTile] & 1) == 0)) { // NE from W
										geom[nTile][2] |= geom[nTile][4];
										geom[nTile][4]  = geom[nTile][2];
									}
								}
							}
							boolean propagated = false;
							while (!propagated){
								propagated = true;
								for (int i = 0; i < 8; i++){
									for (int j = 0; j < 8; j++){
										if ((geom[nTile][i] & (1 <<j)) != 0) {
											if (geom[nTile][i] != geom[nTile][j]) {
												propagated = false;
												geom[nTile][i] |= geom[nTile][j];
												geom[nTile][j] = geom[nTile][i];
											}
										}
									}
								}
							}
						}
						// try break horizontally (always - to the right)
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return geom;
	}
	
	/**
	 *  For each cluster returns 3 arrays: indiced of internal cells, indices of fixed border cells (alpha = 0, disparity fixed)
	 *  and floating border cells (alpha = 0, disparity - from neighbors)
	 * @param diag_en true - 8 directions, false - only 4 orthogonal
	 * @param neighbors array of per-cell bitmaps - which neighbors current cell is connected to (+1 - to N, +2 - to NE,...+128 = to NW)
	 * @param selected boolean array of celected (internal) cells
	 * @param border boolean array of border cells (added for alpha)
	 * @param threadsMax maximal number of threads for multi-threaded application
	 * @param debugLevel  
	 * @return array [cluster_number][list_type: 0..2][element index]
	 */
	public int [][][] extractNonOlerlap(
			final boolean     diag_en,
			final int     []  neighbors0, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
			final boolean []  selected0, // only inner?
			final boolean []  border0,  
			final int         threadsMax,      // maximal number of threads to launch                         
			final int         debugLevel)
			{
		
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		
		int [][] flaps = createOverlapGeometry(
						neighbors0, // +1 - up (N), +2 - up-right - NE, ... +0x80 - NW
						selected0, // only inner?
						border0,
						threadsMax,      // maximal number of threads to launch                         
						debugLevel);
		
		// adding 1-tile frame around to avoid checking for the borders,
		final int tilesX2 = tilesX + 2;
		final int tilesY2 = tilesY + 2;
		int len2 = tilesX2*tilesY2;
		int [] dirs4 = {-tilesX2, 1, tilesX2,-1};
		int [] dirs8 = {-tilesX2,  -tilesX2 + 1, 1, tilesX2 +1, tilesX2, tilesX2 - 1, -1, -tilesX2 - 1};
		int [] dirs = diag_en? dirs8 : dirs4;
		int [] bits4 = {0, 2, 4, 6};
		int [] bits8 = {0, 1, 2, 3, 4, 5, 6, 7};
		int [] neib_bits = diag_en? bits8 : bits4;

		boolean [] selected2 =    new boolean [len2];
		boolean [] border2 =      new boolean [len2];
		int []     neighbors2 =   new int [len2];
		int [][]   flaps2 =        new int [len2][];
		// not needed in java, but for porting
		int [][] startLenInc = {{0,tilesX2-1, 1},{tilesX2-1, len2, tilesX2},{0, len2, tilesX2},{(tilesY2-1)*tilesX2,len2,1}};
		for (int n = 0; n < startLenInc.length; n++) {
			for (int i = startLenInc[n][0]; i <startLenInc[n][1]; i+= startLenInc[n][2]){
				selected2[i] = false;
				border2[i] =   false;
				neighbors2[i] = 0;
				flaps2[i] =     null;
			}
		}
		
		for (int i = 0; i < tilesY; i++) {
			for (int j = 0; j < tilesX; j++) {
				int indx_src = i * tilesX + j;
				int indx_dst = (i+1)*tilesX2 + j + 1;
				selected2[indx_dst] =  selected0[indx_src] && !border0[indx_src]; // so selected2 does not include border (it should not)
				border2[indx_dst] =    border0[indx_src];
				neighbors2[indx_dst] = neighbors0[indx_src];
				flaps2[indx_dst] =     flaps[indx_src];
			}
		}
		boolean [] clustersLeft2 = selected2.clone();
		class ClusterData{
			int [] internal;       // normal cell, obey calculated disparity
			int [] border_fixed;   // border (fade alpha), obey calculated disparity
			int [] border_float;   // border (fade alpha), reconstruct disparity from neighbors
			int indexFromExtended(int ext_index)
			{
				int tY = ext_index / tilesX2;
				int tX = ext_index % tilesX2;
				return (tY-1) * tilesX + (tX-1);
			}
		}
		final int INNER_USED = 0; // used fro normal cells (not border 
		final int UNUSED = -1;
		final int FIXED_BORDER = -2;
		int [] wave = new int [len2];
		
//		final int dbg_tile2 = 35990; // extended tileX = 130, tileY = 110
//		final int dbg_tile2 = 28168; // extended tileX = 132, tileY = 86 (131/85)
//		final int dbg_tile2 = 27842; // extended tileX = 132, tileY = 85 (131/84)
//		final int dbg_tile2 = 43359; // extended tileX = 132, tileY = 85 (131/84)
//		final int dbg_tile2 = 36601; // extended tileX = 132, tileY = 85 (131/84)
//		final int dbg_tile2 = 39276; // extended tileX = 156, tileY = 120 (155/119)
//		final int dbg_tile2 = 39277; // extended tileX = 157, tileY = 120 (155/119)
		final int dbg_tile2 = 39602; // extended tileX = 156, tileY = 121 (155/120)
		
		for (int i = 0; i < len2; i++) wave[i] = UNUSED; // 0 - used normally, >0 (from flaps) - added flaps/border,-2 - already found impossible
		ArrayList<ClusterData> clusterList = new ArrayList<ClusterData>(); 
		ArrayList<Integer> waveList =   new ArrayList<Integer>(); // builds list for clusterList.internal
		ArrayList<Integer> borderFloatList = new ArrayList<Integer>(); // builds list for clusterList.border_float
		ArrayList<Integer> borderFixedList = new ArrayList<Integer>(); // builds list for clusterList.border_fixed
		for (int start_indx = 0; start_indx < len2; start_indx++) if (clustersLeft2[start_indx]){ // found first pixel of a new cluster
			waveList.clear();
			borderFloatList.clear();
			borderFixedList.clear();
			// first internal tile is always OK to use
			wave[start_indx] = INNER_USED; // normal cell
			waveList.add(new Integer(start_indx));
			int frontTail = 0; // advance frontTail, keep list
			// grow wave list (not just wave), add to impossibleList when needed
			while (frontTail < waveList.size()) {
				int indx0 = waveList.get(frontTail++); // just advance pointer, do not remove
				if (indx0 == dbg_tile2){
					System.out.println("extractNonOlerlap(): indx0 = "+indx0);
				}
				// go around and if not yet assigned, add border (should be possible) or new cell (if it is possible), otherwise add it
				// as a fixed border
				for (int dir0 = 0; dir0 < dirs.length; dir0++){
					int indx1 = indx0 + dirs[dir0]; // guaranteed to be inside bounds as indx0 is at least 1 tile away.
					if (indx1 == dbg_tile2){
						System.out.println("extractNonOlerlap(): indx1 = "+indx1);
					}
					if (wave[indx1] == UNUSED){
						if ((neighbors2[indx0] & (1 << neib_bits[dir0])) != 0){ // actual cell - check if it is possible
							if (!clustersLeft2[indx1]){ // connected but not needed - already encoded, used as fixed border
								borderFixedList.add(new Integer(indx1));
								wave[indx1] = FIXED_BORDER;
							} else  {
								possible_label:{
								for (int dir1 = 0; dir1 < dirs.length; dir1++){
									int dir1_back = (dir1 +4) % 8;
									int indx2 = indx1 + dirs[dir1]; // not guaranteed to be inside bounds, check neighbors and flaps first
									if ((neighbors2[indx1] & (1 << neib_bits[dir1])) != 0) { // connected to other in this direction
										if (wave[indx2] > 0){
											// TODO:add as a fixed border
											borderFixedList.add(new Integer(indx1));
											wave[indx1] = FIXED_BORDER;
											if (indx1 == dbg_tile2){
												System.out.println("extractNonOlerlap() A: indx2 = "+indx2);
											}
											
											break possible_label; // impossible
										}
									} else if (flaps2[indx1][neib_bits[dir1_back]] != 0){ // that neighbor is not connected, but here can be reached from there
										///											if ((wave[indx2] != UNUSED) && (wave[indx2] != flaps2[indx1][neib_bits[dir1]])) {
										///											if ((wave[indx2] != UNUSED) && (wave[indx2] != flaps2[indx2][neib_bits[dir1]])) { // there when reached from here
										// got here bumping into FIXED_BORDER from previous shell. Should be treated as unused?
										if ((wave[indx2] != UNUSED) && (wave[indx2] != FIXED_BORDER) && (wave[indx2] != flaps2[indx2][neib_bits[dir1]])) { // there when reached from here

											if (indx1 == dbg_tile2){
												System.out.println("extractNonOlerlap() B: indx2 = "+indx2);
											}
											// TODO:add as a fixed border ????
											borderFixedList.add(new Integer(indx1));
											wave[indx1] = FIXED_BORDER;
											break possible_label; // impossible
										}
										// should be marked here - no, only afer all tested!
									}
								}
								// all conditions are met, add this cell and flaps around it
								waveList.add(new Integer(indx1));
								wave[indx1] = INNER_USED;
								// add all flaps around
								for (int dir1 = 0; dir1 < dirs.length; dir1++){
									int dir1_back = (dir1 +4) % 8;
									// not connected, but reachable
									if (indx1 == dbg_tile2){
										System.out.println("extractNonOlerlap(): indx1 = "+indx1);
										System.out.println("extractNonOlerlap(): dir1 = "+dir1);
										System.out.println("extractNonOlerlap(): neib_bits[dir1] = "+neib_bits[dir1]);
										System.out.println("extractNonOlerlap(): neighbors2[indx1] = "+neighbors2[indx1]);
										System.out.println("extractNonOlerlap(): flaps2[indx1][neib_bits["+dir1_back+"]] = "+flaps2[indx1][neib_bits[dir1_back]]);
									}

									///										if (((neighbors2[indx1] & (1 << neib_bits[dir1])) == 0) &&  (flaps2[indx1][neib_bits[dir1]] != 0))  {
									// that neighbor is not connected, but here can be reached from there
									if (((neighbors2[indx1] & (1 << neib_bits[dir1])) == 0) &&  (flaps2[indx1][neib_bits[dir1_back]] != 0))  {
										int indx2 = indx1 + dirs[dir1]; // not guaranteed to be inside bounds, check neighbors and flaps first
										///											if (wave[indx2] == UNUSED){
										if ((wave[indx2] == UNUSED) || (wave[indx2] == FIXED_BORDER)){
											///												wave[indx2] = flaps2[indx1][neib_bits[dir1]];
											wave[indx2] = flaps2[indx2][neib_bits[dir1]];
											borderFloatList.add(new Integer(indx2));
										}
									}
								}								
							}
							/// impossible to add - should it be a fixed border??
//							borderFixedList.add(new Integer(indx1));
//							wave[indx1] = FIXED_BORDER;

							}

						} else { // unused but not connected (only around initial cell - no need to check possibility for flaps
							// add flaps, normally the should be added when a new possible cell is added, only for the first one needs here
							// got here - not after first. Should it be marked when considering 27842?
							//							wave[indx1] = flaps2[indx0][neib_bits[dir0]];
							if (flaps2[indx1] == null){
								// null may be on the border (not reachable tiles even if they have neighbors)
								int tileX2 = indx1 % tilesX2;
								int tileY2 = indx1 / tilesX2;
								if ((tileX2 == 0) || (tileX2 == tilesX2 - 1) || (tileY2 == 0) || (tileY2 == tilesY2 - 1)) {
									if (debugLevel >-1){
										System.out.println("Unreachable tile on the border, OK: tilesX2="+tileX2+", tilesY2="+tileY2);
									}
								} else {
									System.out.println("BUG1061!! indx0="+indx0+ ", indx1="+indx1+", tilesX2="+tileX2+", tilesY2="+tileY2);
								}
							} else {
								wave[indx1] = flaps2[indx1][neib_bits[dir0]]; 
								borderFloatList.add(new Integer(indx1));
							}
						}
					}
				}
			}
			// wave died, encode what it got and clean up
			ClusterData clusterData =  new ClusterData();
			clusterData.internal =     new int [waveList.size()];
			clusterData.border_fixed = new int [borderFixedList.size()];
			clusterData.border_float = new int [borderFloatList.size()];
			int indx = 0;
			for (Integer I : waveList){
				clusterData.internal[indx++] = clusterData.indexFromExtended(I);
				// cleanup;
				clustersLeft2[I] = false; // used in this cluster
				wave[I] = UNUSED;
			}
			indx = 0;
			for (Integer I : borderFixedList){
				clusterData.border_fixed[indx++] = clusterData.indexFromExtended(I); 	
				wave[I] = UNUSED; // TODO: handle when encoding around !!
			}
			indx = 0;
			for (Integer I : borderFloatList){
				clusterData.border_float[indx++] = clusterData.indexFromExtended(I); 	
				wave[I] = UNUSED; // TODO: handle when encoding around !!
			}
			clusterList.add(clusterData);
			if (debugLevel > -1){
				System.out.println("extractNonOlerlap(): added #"+clusterList.size()+" internal:"+clusterData.internal.length+
						" fixed:"+clusterData.border_fixed.length+" float:"+clusterData.border_float.length);
			}
		}
		// convert list to array
		

		int [][][] clusters = new int [clusterList.size()][][];
		int indx = 0;
		for (ClusterData cd : clusterList){
			int [][] cda = {cd.internal, cd.border_fixed, cd.border_float};
			clusters[indx++] = cda;
		}
		return clusters;
	}

	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	static Thread[] newThreadArray(int maxCPUs) {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
	/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	public static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try
		{   
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}

}
