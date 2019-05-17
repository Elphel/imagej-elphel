package com.elphel.imagej.tileprocessor;
/**
 **
 ** Conflicts - Represent "conflicts" (instances of Conflict) between connected supertiles
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  Conflicts.java is free software: you can redistribute it and/or modify
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
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;

public class Conflicts {
	private int [] num_ortho_diag_ortho = new int [8];
	private int [] num_ortho_ortho_diag = new int [16];
	private  int [] num_all_conflicts =        new int [24];
	private  int    num_ortho_incompat = 0;
	private  int    num_ortho_dual = 0;
	private  int    num_conflicts = 0;
	private  SuperTiles st = null;
	
	public Conflicts(
			SuperTiles st){
		this.st = st;
	}
	
	public Conflicts(
			int [][][] conflicts,
			SuperTiles st,
			int        debugLevel)
	{
		this.st = st;
		addConflicts( conflicts, debugLevel);
	}
	public Conflicts(
			int [][][] conflicts,
			SuperTiles st)
	{
		this.st = st;
		addConflicts( conflicts, -1);
	}
	
	public Conflicts(
			Conflict [][] conflicts,
			SuperTiles st,
			int        debugLevel)
	{
		this.st = st;
		addConflicts( conflicts, debugLevel);
	}
	public Conflicts(
			Conflict [][] conflicts,
			SuperTiles st)
	{
		this.st = st;
		addConflicts( conflicts, -1);
	}
	public void addConflicts(Conflicts conflicts){
		addsubConflicts(conflicts, false);
	}
	public void subConflicts(Conflicts conflicts){
		addsubConflicts(conflicts, true);
	}
	
	public void addsubConflicts(Conflicts conflicts, boolean sub)
	{
		int s = sub? -1: 1;
		for (int i = 0; i < num_ortho_diag_ortho.length; i++) num_ortho_diag_ortho[i] += s * conflicts.num_ortho_diag_ortho[i];  
		for (int i = 0; i < num_ortho_ortho_diag.length; i++) num_ortho_ortho_diag[i] += s * conflicts.num_ortho_ortho_diag[i];  
		for (int i = 0; i < num_all_conflicts.length; i++)    num_all_conflicts[i] += s * conflicts.num_all_conflicts[i];  
		num_ortho_incompat += s * conflicts.num_ortho_incompat;
		num_ortho_dual     += s * conflicts.num_ortho_dual;
		num_conflicts      += s * conflicts.num_conflicts;
		
	}

	public void resetConflicts()
	{
		num_ortho_diag_ortho = new int [8];
		num_ortho_ortho_diag = new int [16];
		num_all_conflicts =    new int [24];
		num_ortho_incompat =   0;
		num_ortho_dual =       0;
		num_conflicts =        0;
	}
	
	
	
	public int addConflicts(
			int [][][] conflicts,
			int        debugLevel)
	{
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null){
			for (int nc = 0; nc < conflicts[nsTile].length; nc++){
				Conflict conf = new Conflict(nsTile, conflicts[nsTile][nc]);
				if (conf.getNumOrthoDiagOrthoConflicts() > 0)  num_ortho_diag_ortho[conf.getNumOrthoDiagOrthoConflicts() - 1]++;
				if (conf.getNumOrthoOrthoDiagConflicts() > 0)  num_ortho_ortho_diag[conf.getNumOrthoOrthoDiagConflicts() - 1]++;
				if (conf.getNumConflicts() > 0)                num_all_conflicts[conf.getNumConflicts() - 1]++;
				num_ortho_incompat += conf.getIncompatibleOrthoDiagOrthoConflicts();
				num_ortho_dual += conf.getDualTriOrthoDiagOrthoConflicts();
				num_conflicts += conf.getNumConflicts();
				if (debugLevel > 0){
					int tilesX =        st.tileProcessor.getTilesX();
					int superTileSize = st.tileProcessor.getSuperTileSize();
					int stilesX = (tilesX + superTileSize -1)/superTileSize;  
					int ty = nsTile / stilesX;
					int tx = nsTile % stilesX;
					printConflict("addConflicts() nsTile = "+nsTile+" ["+tx+":"+ty+"] ", conf);
				}
			}
		}
		return num_conflicts;
	}

	public int addConflicts(
			Conflict [][] conflicts,
			int           debugLevel)
	{
		for (int nsTile = 0; nsTile < conflicts.length; nsTile++) if (conflicts[nsTile] != null){
			for (int nc = 0; nc < conflicts[nsTile].length; nc++){
				Conflict conf = conflicts[nsTile][nc];
				if (conf.getNumOrthoDiagOrthoConflicts() > 0)  num_ortho_diag_ortho[conf.getNumOrthoDiagOrthoConflicts() - 1]++;
				if (conf.getNumOrthoOrthoDiagConflicts() > 0)  num_ortho_ortho_diag[conf.getNumOrthoOrthoDiagConflicts() - 1]++;
				if (conf.getNumConflicts() > 0)                num_all_conflicts[conf.getNumConflicts() - 1]++;
				num_ortho_incompat += conf.getIncompatibleOrthoDiagOrthoConflicts();
				num_ortho_dual += conf.getDualTriOrthoDiagOrthoConflicts();
				num_conflicts += conf.getNumConflicts();
				if (debugLevel > 0){
					int tilesX =        st.tileProcessor.getTilesX();
					int superTileSize = st.tileProcessor.getSuperTileSize();
					int stilesX = (tilesX + superTileSize -1)/superTileSize;  
					int ty = nsTile / stilesX;
					int tx = nsTile % stilesX;
					printConflict("addConflicts() nsTile = "+nsTile+" ["+tx+":"+ty+"] ", conf);
				}
			}
		}
		return num_conflicts;
	}
	
	public int [] getNumOrthoDiagOrtho(){
		return num_ortho_diag_ortho;
	}
	public int [] getNumOrthoOrthoDiag(){
		return num_ortho_ortho_diag;
	}
	public int [] getNumAllConflicts(){
		return num_all_conflicts;
	}
	public int getNumOrthoIncompat(){
		return num_ortho_incompat;
	}
	public int getNumOrthoDual(){
		return num_ortho_dual;
	}
	public int getNumConflicts(){
		return num_conflicts;
	}
	
	public void printConflict(String prefix, Conflict conf)
	{
		System.out.println(prefix+conf.toString());
	}
	
	public int numBetterWorse(
			boolean better,
			boolean use_all,
			boolean use_odo,
			boolean use_ood)
	{
		int num = 0;
		if (use_all) {
			for (int i = 0; i < num_all_conflicts.length; i++){
				if (better?(num_all_conflicts[i] < 0):(num_all_conflicts[i] > 0)) num++;
			}
		}
		if (use_odo) {
			for (int i = 0; i < num_ortho_diag_ortho.length; i++){
				if (better?(num_ortho_diag_ortho[i] < 0):(num_ortho_diag_ortho[i] > 0)) num++;
			}
			if (better? (num_ortho_incompat < 0) : (num_ortho_incompat > 0)) num++;
			if (better? (num_ortho_dual < 0) : (num_ortho_dual > 0)) num++;

		}
		if (use_ood) {
			for (int i = 0; i < num_ortho_ortho_diag.length; i++){
				if (better?(num_ortho_ortho_diag[i] < 0):(num_ortho_ortho_diag[i] > 0)) num++;
			}
		}
		return num;
	}
	public int sumConflicts(
			boolean use_all,
			boolean use_odo,
			boolean use_ood)
	{
		int num = 0;
		if (use_all) {
			for (int i = 0; i < num_all_conflicts.length; i++){
				num += num_all_conflicts[i];
			}
		}
		if (use_odo) {
			for (int i = 0; i < num_ortho_diag_ortho.length; i++){
				num +=num_ortho_diag_ortho[i];
			}
		}
		if (use_ood) {
			for (int i = 0; i < num_ortho_ortho_diag.length; i++){
				num +=num_ortho_ortho_diag[i];
			}
		}
		return num;
	}
	
	
	public void printConflictSummary(
			String prefix,
			boolean use_all,
			boolean use_odo,
			boolean use_ood)
	{
		System.out.print(prefix);
		if (use_all) {
			for (int i = 0; i < num_all_conflicts.length; i++){
				if (num_all_conflicts[i] != 0) System.out.print(" all_"+(i + 1)+": "+num_all_conflicts[i]);
			}
		}

		if (use_odo) {
			for (int i = 0; i < num_ortho_diag_ortho.length; i++){
				if (num_ortho_diag_ortho[i] != 0) System.out.print(" odo_"+(i + 1)+": "+num_ortho_diag_ortho[i]);
			}
			if (num_ortho_incompat != 0) {
				System.out.print(" number of incompatible odo triangles = " + num_ortho_incompat);
			}
			if (num_ortho_dual != 0) {
				System.out.print(" number of dual odo triangles = " + num_ortho_dual);
			}
		}
		if (use_ood) {
			for (int i = 0; i < num_ortho_ortho_diag.length; i++){
				if (num_ortho_ortho_diag[i] != 0) System.out.print(" ood_"+(i + 1)+": "+num_ortho_ortho_diag[i]);
			}
		}		
		System.out.println();		
	}
//	TilePlanes.PlaneData [][] planes =  null;

	/**
	 * Find "triangular" conflicts after running selectNeighborPlanesMutual.
	 * Such conflicts happen is when starting N (or other ortho direction) from some node,
	 * then turning 135 degrees right, and then 135 right again it will get to the different
	 * layer than started (and all 3 connections exist. 
	 * @param debugLevel
	 * @return 3-d array, first index - tile number, 2-nd - tile conflict number. Each conflict
	 * is stored as {start layer, end layer, bitmask of start directions} (+1 - N, +2 - E, +4 - S, +8 - W)
	 */

	
	public int [][][] detectTriangularConflicts(
			final int debugLevel)
	{
		final int tilesX =        st.getTileProcessor().getTilesX();
		final int tilesY =        st.getTileProcessor().getTilesY();
		final int superTileSize = st.getTileProcessor().getSuperTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final int [][][] conflicts = new int [st.getPlanes().length][][];
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(st.getTileProcessor().threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						if ( st.getPlanes()[nsTile0] != null) {
							conflicts[nsTile0] =  detectTriangularTileConflicts(
									nsTile0,
									null,    // HashMap<Integer,Integer> replacement_tiles, //
									null,    // int [][][] replacement_neibs,
									tnSurface);
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		if (debugLevel > -1){
			addConflicts(
					conflicts,
					debugLevel);
			printConflictSummary("Detected conflicts:",true,false,false);
			printConflictSummary("Detected ortho-diagonal-ortho conflicts:",false, true, false);
			printConflictSummary("Detected ortho-ortho-diagonal conflicts:",false, false, true);
		}
		return conflicts;
	}
	
	/**
	 * Calculate cost of all conflicts around supertile nsTile0 by adding "star" weight of each tile involved.
	 * Tile weight can be either from the planes[][] array or from the replacement values in replacement_val_weights
	 * (when the tiles configuration is not yet committed)
	 * @param nsTile0 supertile index to process
	 * @param scaleStartEnd for each conflict triangle add start and end tiles (at the center nsTile0) scaled by this value
	 * @param conflicts array of calculated conflicts (each is {start_layer, end_layer, direction/type bitmask} or null,
	 * in that case it will be calculated
	 * @param replacement_tiles a map of supertile indices to replacement indices (for replacement_neibs and
	 *        replacement_val_weights) or null. If not null, and entry for the supertile full index exists, 
	 *        replacement_neibs and replacement_val_weights will be used, otherwise planes[][] data.
	 *        replacement_tiles may be null, in that case planes[][] data will be used unconditionally.
	 * @param replacement_neibs array of neighbors to use instead of the planes data. First index is an index
	 *        of the replacement supertile (values in the  replacement_tiles map), second - layer number and
	 *        the 3-rd one - direction index for 8 connections: N, NE, E...NW. Value is the destination layer
	 *        number. Can be null if not used (when replacement_tiles is null)
	 * @param replacement_val_weights similar array for tile value/weight data (per-tile index, per layer).
	 *  The innermost data is a tuple {value, weight}
	 * @param tnSurface TileNeibs instance to navigate through the 2-d array encoded in linescan order
	 * @return sum ao the weights of all conflicts fro this tile
	 */
	
	public double getConflictsCost(
			int           nsTile0,
			double        scaleStartEnd, // include start and and layer tiles in the center in overall cost for each triangle (1.0)
			int [][]      conflicts, //
			HashMap<Integer,Integer> replacement_tiles, // null is OK
			int [][][]    replacement_neibs,               // null OK if  replacement_tiles == null
			double [][][] replacement_val_weights,
			TileNeibs tnSurface)
	{
		TilePlanes.PlaneData [][] planes = st.getPlanes();
		// generate conflicts if not provided
		if (conflicts == null) {
			conflicts = detectTriangularTileConflicts(
					nsTile0,
					replacement_tiles, //
					replacement_neibs,
					tnSurface);
		}
		double cost = 0.0;
		if (conflicts != null){
			Integer isTile = (replacement_tiles != null) ? replacement_tiles.get(nsTile0) : null;
			for (int nConfl = 0; nConfl < conflicts.length; nConfl++){
				Conflict conflict = new Conflict(nsTile0, conflicts[nConfl]);
				int start_layer = conflict.getStartLayer();
				int end_layer =   conflict.getEndLayer();
				int [] neibs1 = ((isTile == null) ?
						planes[nsTile0][start_layer].getNeibBest() :
							replacement_neibs[isTile][start_layer]);
				int [] neibs2 = ((isTile == null) ?
						planes[nsTile0][end_layer].getNeibBest() :
							replacement_neibs[isTile][end_layer]);

				int [][] involved = conflict.getInvolvedTiles();
				for (int nTri = 0; nTri < involved.length; nTri++){
					int [] nsTiles = {
							tnSurface.getNeibIndex(nsTile0, involved[nTri][0]),
							tnSurface.getNeibIndex(nsTile0, involved[nTri][1])};
					int [] layers = {neibs1[involved[nTri][0]], neibs2[involved[nTri][1]]};
					for (int it = 0; it < nsTiles.length; it++){
						Integer isTile1 = (replacement_tiles != null) ? replacement_tiles.get(nsTiles[it]) : null;
						if (isTile1 != null){
							cost += replacement_val_weights[isTile1][layers[it]][1];
						} else {
							cost += planes[nsTiles[it]][layers[it]].getStarValueWeight()[1];
						}
					}
				}
				if (scaleStartEnd != 0.0) {
					if (isTile != null){
						cost += scaleStartEnd * (replacement_val_weights[isTile][start_layer][1]+
								replacement_val_weights[isTile][end_layer][1]);
					} else {
						cost += scaleStartEnd * (planes[nsTile0][start_layer].getStarValueWeight()[1]+
								planes[nsTile0][end_layer].getStarValueWeight()[1]);
					}
				}
			}
		}
		return cost;
	}


	public double getConflictsCost(
			int           nsTile0,
			double        scaleStartEnd, // include start and and layer tiles in the center in overall cost for each triangle (1.0)
			int [][]      conflicts, //
			HashMap<Integer,Integer> replacement_tiles, // null is OK
			int [][][]    replacement_neibs,               // null OK if  replacement_tiles == null
			ConnectionCosts connectionCosts,
			TileNeibs tnSurface)
	{
		TilePlanes.PlaneData [][] planes = st.getPlanes();
		// generate conflicts if not provided
		if (conflicts == null) {
			conflicts = detectTriangularTileConflicts(
					nsTile0,
					replacement_tiles, //
					replacement_neibs,
					tnSurface);
		}
		double cost = 0.0;
		if (conflicts != null){
			Integer isTile = (replacement_tiles != null) ? replacement_tiles.get(nsTile0) : null;
			for (int nConfl = 0; nConfl < conflicts.length; nConfl++){
				Conflict conflict = new Conflict(nsTile0, conflicts[nConfl]);
				int start_layer = conflict.getStartLayer();
				int end_layer =   conflict.getEndLayer();
				int [] neibs1 = ((isTile == null) ?
						planes[nsTile0][start_layer].getNeibBest() :
							replacement_neibs[isTile][start_layer]);
				int [] neibs2 = ((isTile == null) ?
						planes[nsTile0][end_layer].getNeibBest() :
							replacement_neibs[isTile][end_layer]);

				int [][] involved = conflict.getInvolvedTiles();
				for (int nTri = 0; nTri < involved.length; nTri++){
					int [] nsTiles = {
							tnSurface.getNeibIndex(nsTile0, involved[nTri][0]),
							tnSurface.getNeibIndex(nsTile0, involved[nTri][1])};
					int [] layers = {neibs1[involved[nTri][0]], neibs2[involved[nTri][1]]};
					for (int it = 0; it < nsTiles.length; it++){
						cost += connectionCosts. getValWeightLast(
								nsTiles[it], // int nsTile,
								layers[it], // int nl,
								false)[1]; // boolean initialValue)
					}
				}
				if (scaleStartEnd != 0.0) {
					cost += scaleStartEnd * connectionCosts. getValWeightLast(
							nsTile0,     // int nsTile,
							start_layer, // int nl,
							false)[1];   // boolean initialValue)
					cost += scaleStartEnd * connectionCosts. getValWeightLast(
							nsTile0,     // int nsTile,
							end_layer, // int nl,
							false)[1];   // boolean initialValue)
				}
			}
		}
		return cost;
	}


	
	
	public int [][] detectTriangularTileConflicts(
			int nsTile0,
			HashMap<Integer,Integer> replacement_tiles, // null is OK - will use only planes data
			int [][][] replacement_neibs,               // null OK if  replacement_tiles == null
			TileNeibs tnSurface)
	{
		TilePlanes.PlaneData [][] planes = st.getPlanes();
		ArrayList<Conflict> conflicts_list= new ArrayList<Conflict>();
		if ( planes[nsTile0] != null) {
			Integer repl_indx = (replacement_tiles != null) ? replacement_tiles.get(new Integer(nsTile0)): null;
			for (int np0 = 0; np0 < planes[nsTile0].length; np0++){
				if (planes[nsTile0][np0] != null) {
					int [] neibs0 = ((repl_indx == null) || (replacement_neibs[repl_indx][np0] == null)) ?
							planes[nsTile0][np0].getNeibBest() :
								replacement_neibs[repl_indx][np0];
					for (int dir = 0; dir < 8; dir +=2){
						int np1;
						np1= neibs0[dir]; // planes[nsTile0][np0].getNeibBest(dir);
						if (np1 >=0) {
							int nsTile1 = tnSurface.getNeibIndex(nsTile0, dir);
							if (nsTile1 >= 0){
								Integer repl_indx1 = (replacement_tiles != null) ? replacement_tiles.get(new Integer(nsTile1)): null;
								int [] neibs1 = (((repl_indx1 == null) || (replacement_neibs[repl_indx1][np1] == null) )?
										planes[nsTile1][np1].getNeibBest() :
											replacement_neibs[repl_indx1][np1]);
								int dir1 = (dir + 3) % 8;
								int np2 = neibs1[dir1]; // planes[nsTile1][np1].getNeibBest(dir1);
								if (np2 >= 0){
									int nsTile2 = tnSurface.getNeibIndex(nsTile1, dir1);
									if (nsTile2 >= 0){
										Integer repl_indx2 = (replacement_tiles != null) ? replacement_tiles.get(new Integer(nsTile2)): null;
										int [] neibs2 = (((repl_indx2 == null) || (replacement_neibs[repl_indx2][np2] == null)) ?
												planes[nsTile2][np2].getNeibBest() :
													replacement_neibs[repl_indx2][np2]);
										int dir2 = (dir1 + 3) % 8;
										int np3 = neibs2[dir2]; // planes[nsTile2][np2].getNeibBest(dir2);
										if ((np3 >= 0) && (np3 != np0)){
											Conflict conflict = new Conflict(nsTile0, np0, np3, dir/2);
											label_apply:
											{
												for (Conflict conf_old:conflicts_list){
													if (conf_old.combine(conflict)){
														break label_apply;
													}
												}
												conflicts_list.add(conflict);
											}
										}
									}
								}
								// ortho-ortho-diagonal, right hand from np0
								dir1 = (dir + 2) % 8;
								np2 = neibs1[dir1]; // planes[nsTile1][np1].getNeibBest(dir1);
								if (np2 >= 0){
									int nsTile2 = tnSurface.getNeibIndex(nsTile1, dir1);
									if (nsTile2 >= 0){
										Integer repl_indx2 = (replacement_tiles != null) ? replacement_tiles.get(new Integer(nsTile2)): null;
										int [] neibs2 = (((repl_indx2 == null) || (replacement_neibs[repl_indx2][np2] == null)) ?
												planes[nsTile2][np2].getNeibBest() :
													replacement_neibs[repl_indx2][np2]);
										int dir2 = (dir1 + 3) % 8;
										int np3 = neibs2[dir2]; // planes[nsTile2][np2].getNeibBest(dir2);
										if ((np3 >= 0) && (np3 != np0)){
											Conflict conflict = new Conflict(nsTile0, np0, np3, dir/2, true); // ood, right
											label_apply:
											{
												for (Conflict conf_old:conflicts_list){
													if (conf_old.combine(conflict)){
														break label_apply;
													}
												}
												conflicts_list.add(conflict);
											}
										}
									}
								}
								// ortho-ortho-diagonal, left hand from np0
								dir1 = (dir + 6) % 8;
								np2 = neibs1[dir1]; // planes[nsTile1][np1].getNeibBest(dir1);
								if (np2 >= 0){
									int nsTile2 = tnSurface.getNeibIndex(nsTile1, dir1);
									if (nsTile2 >= 0){
										Integer repl_indx2 = (replacement_tiles != null) ? replacement_tiles.get(new Integer(nsTile2)): null;
										int [] neibs2 = (((repl_indx2 == null) || (replacement_neibs[repl_indx2][np2] == null)) ?
												planes[nsTile2][np2].getNeibBest() :
													replacement_neibs[repl_indx2][np2]);
										int dir2 = (dir1 + 5) % 8;
										int np3 = neibs2[dir2]; // planes[nsTile2][np2].getNeibBest(dir2);
										if ((np3 >= 0) && (np3 != np0)){
											Conflict conflict = new Conflict(nsTile0, np0, np3, dir/2, false); // ood, left
											label_apply:
											{
												for (Conflict conf_old:conflicts_list){
													if (conf_old.combine(conflict)){
														break label_apply;
													}
												}
												conflicts_list.add(conflict);
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
		if(conflicts_list.isEmpty()){
			return null;
		}
		int [][] conflicts = new int [conflicts_list.size()][];
		int indx=0;
		for (Conflict conflict:conflicts_list){
			conflicts[indx++] = conflict.toArray();
		}
		return conflicts;
	}
	
}
