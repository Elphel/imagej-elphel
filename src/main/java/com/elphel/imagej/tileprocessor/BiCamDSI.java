package com.elphel.imagej.tileprocessor;
/**
 ** BiCamDSI - Building DSI using correlation between two quad cameras
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  BiCamDSI.java is free software: you can redistribute it and/or modify
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

import com.elphel.imagej.common.ShowDoubleFloatArrays;


public class BiCamDSI {
	TileNeibs         tnImage; //  = new TileNeibs(tilesX, tilesY)
	int               threadsMax;
	public ArrayList<BiScan> biScans;



	public int addBiScan(
			double [] disparity,
			double [] strength,
			boolean [] trusted,
			boolean [] disabled,
			int        scan_type) {

		if (biScans == null) biScans = new ArrayList<BiScan>();
		biScans.add(new BiScan(this, biScans.size(), disparity, strength, trusted, disabled, scan_type));
		return biScans.size()-1;

	}

	public int addBiScan(
			double [][] disparity_bimap,
			int        scan_type
			) {
		return addBiScan(
				disparity_bimap[ImageDtt.BI_TARGET_INDEX], // double [] disparity,
				disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX], // double [] strength,
				null, // boolean [] trusted,
				null, // boolean [] disabled)
				scan_type);
	}

	/**
	 * Get last BiScan of specific type
	 * @param scan_type -1 - any (returns last scan), >=0 must match type. If none of the type exists, returns null
	 * @return
	 */
	public BiScan getLastBiScan( int scan_type) {
		if (scan_type < 0) 	return getBiScan(-1);
		if (biScans.isEmpty()) return null;
		int last_index = -1;
		for (int i = 0; i < biScans.size(); i++) {
			if (biScans.get(i).scan_type == scan_type) {
				last_index = i;
			}
		}
		if (last_index >= 0) {
			return biScans.get(last_index);
		} else {
			return null;
		}
	}
	public BiScan getBiScan(int indx) {
		if (biScans.isEmpty()) return null;
		return biScans.get((indx>=0)? indx: (biScans.size() - 1));
	}

	public double [] getTargetDisparity(int indx) {
		BiScan biScan = getBiScan(indx);
		if (biScan == null) return null;
		return biScan.target_disparity;
	}

	public BiCamDSI(
			int tilesX,
			int tilesY,
			int threadsMax) {
		this.threadsMax = threadsMax;
		tnImage  = new TileNeibs(tilesX, tilesY);
	}

	private double [][] get9weightMasks(
			int        neib_dist, // >=1
			double     rsigma) {
		int sample_size = 2*neib_dist + 1;
		int center = neib_dist; //+1;
		int last = 2*neib_dist;

		double [][] weight_mask = new double [9][sample_size * sample_size];
		for (int row = 0; row <= neib_dist; row++) {
			int row1 = (last - row);
			for (int col = 0; col < sample_size; col++) {
				weight_mask[0][row *sample_size + col ] = 1.0;
				weight_mask[2][col *sample_size + row1] = 1.0;
				weight_mask[4][row1*sample_size + col ] = 1.0;
				weight_mask[6][col *sample_size + row ] = 1.0;
			}
		}
		for (int i = 0; i < sample_size; i++) {
			for (int j = 0; j <= i; j++) {
				int row0 =  j;
				int col0 = last - i + j;
				int row1 = (last - row0);
				int col1 = (last - col0);
				weight_mask[1][row0*sample_size + col0 ] = 1.0;
				weight_mask[3][col0*sample_size + row1 ] = 1.0;
				weight_mask[5][row1*sample_size + col1 ] = 1.0;
				weight_mask[7][col1*sample_size + row0 ] = 1.0;
			}
		}
		// Fill the center-symmetrical weights
		int k =0;
		double w4 = 0.0;
		double w8 = 0.0;
		int n8 = 0;
		int left_cells = sample_size * (neib_dist + 1); //number of cells in "halves"
		for (; ((2*k + 1)*(2*k + 1)) <= left_cells; k++); // seems that it never equals
		k--; // k is the largest square with odd side and area less than for 8 halves above
		left_cells -= (2*k + 1)*(2*k + 1);
		if (left_cells < 4) {
			w4 = left_cells/4.0;
		} else {
			w4 = 1.0;
			left_cells -=4;
			n8 = left_cells/8;
			left_cells -= 8*n8;
			if (n8 == k) {
				w8 = left_cells/4.0; // using very corners, so there are only 4, not 8 different cells for w8
			} else {
				w8 = left_cells/8.0;
			}


		}
		// Fill center square
		for (int row = neib_dist-k; row <= (neib_dist+k); row++) {
			for (int col = neib_dist-k; col <= (neib_dist+k); col++) {
				weight_mask[8][row*sample_size + col ] = 1.0;
			}
		}
		int low =  center - k - 1;
		int high = center + k + 1;
		weight_mask[8][center*sample_size + low   ] = w4;
		weight_mask[8][low   *sample_size + center] = w4;
		weight_mask[8][center*sample_size + high  ] = w4;
		weight_mask[8][high  *sample_size + center] = w4;
		if (w8 > 0.0) {
			n8++;
		} else {
			w8 = 1.0;
		}
		if (n8 > 0){
			for (int in8 = 1; in8 <= n8; in8++) {
				double w = (in8 == n8)? w8: 1.0;
				weight_mask[8][(center - in8)*sample_size +           low ] = w;
				weight_mask[8][(center + in8)*sample_size +           low ] = w;
				weight_mask[8][low           *sample_size + (center - in8)] = w;
				weight_mask[8][low           *sample_size + (center + in8)] = w;
				weight_mask[8][(center - in8)*sample_size +           high] = w;
				weight_mask[8][(center + in8)*sample_size +           high] = w;
				weight_mask[8][high          *sample_size + (center - in8)] = w;
				weight_mask[8][high          *sample_size + (center + in8)] = w;
			}
		}
		if (rsigma > 0) {
			double sigma = rsigma * center;
			double [] w1d = new double [center+1];
			for (int i = 0; i < w1d.length; i++) w1d[i] = Math.exp(-i*i/(2*sigma*sigma));
			for (int row = 0; row < sample_size; row++) {
				int ady = (row > center)? (row - center) : (center - row);
				for (int col = 0; col < sample_size; col++) {
					int adx = (col > center)? (col - center) : (center - col);
					for (int n = 0; n < weight_mask.length; n++) {
						weight_mask[n][row*sample_size + col] *= w1d[ady]*w1d[adx];
					}
				}
			}
		}
		return weight_mask;
	}

	/**
	 * Select near tiles (those that are not infinity)
	 * @param min_strength minimal tile strength (may be 0 as the tiles are already filtered)
	 * @param infinity_select tile, previously assigned to infinity (or null)
	 * @param selection existing selection (to add to) or null
	 * @param disparity array of tile disparity values (may have NaN-s)
	 * @param strength array of tile strength values
	 * @return new selection (new or added to existing)
	 */

	public boolean [] selectNearObjects(
			double     min_strength,
			boolean [] infinity_select,
			boolean [] selection,
			double []  disparity,
			double []  strength) {
		int num_tiles = disparity.length;
		if (infinity_select == null) {
			infinity_select = new boolean [num_tiles];
			for (int nTile = 0; nTile < num_tiles; nTile++) infinity_select[nTile] = true;
		}
		if (selection == null) selection = new boolean [num_tiles];
		for (int nTile =0; nTile <num_tiles; nTile++) {
			if (!infinity_select[nTile] && (strength[nTile] > 0.0) && (strength[nTile] >= min_strength)) {
				selection[nTile] = true;
			}
		}
		return selection;
	}

	/**
	 * Select small/thin near objects over infinity
	 * @param min_strength minimal strength of the tile and neighbors
	 * @param min_neibs minimal number of neighbors with close disparity
	 * @param min_disparity minimal disaparity
	 * @param disp_atolerance disparity absolute tolerance (to qualify for a similar neighbor)
	 * @param disp_rtolerance disparity relative tolerance (adds to absolute)
	 * @param infinity_select tile, previously assigned to infinity (or null)
	 * @param selection existing selection (to add to) or null
	 * @param disparity array of tile disparity values (may have NaN-s)
	 * @param strength array of tile strength values
	 * @return new selection (new or added to existing)
	 */
	public boolean [] selectNearOverInfinity(
			double     min_strength,
			int        min_neibs,
			double     min_disparity,
			double     disp_atolerance,
			double     disp_rtolerance,
			boolean [] infinity_select,
			boolean [] selection,
			double []  disparity,
			double []  strength) {
		int num_tiles = disparity.length;
		if (infinity_select == null) {
			infinity_select = new boolean [num_tiles];
			for (int nTile = 0; nTile < num_tiles; nTile++) infinity_select[nTile] = true;
		}
		if (selection == null) selection = new boolean [num_tiles];
		boolean [] new_sel = selection.clone();
		//tnImage
		for (int nTile =0; nTile <num_tiles; nTile++)
			if (     infinity_select[nTile] &&
					!selection[nTile] &&
					(strength[nTile] >=  min_strength ) &&
					(disparity[nTile] >= min_disparity)){
			int nn = 0;
			double tolerance = disp_atolerance + disparity[nTile] * disp_rtolerance;
			double min_d = Math.max(disparity[nTile] - tolerance, min_disparity);
			double max_d = disparity[nTile] + tolerance;
			for (int dir = 0; dir < 8; dir++) {
				int nTile1 = tnImage.getNeibIndex(nTile, dir);
				if (    (nTile1 >=0) &&
						((strength [nTile1] >=  min_strength ) || selection[nTile1]) && // strong or already selected
						(disparity[nTile1] >= min_d) &&
						(disparity[nTile1] <= max_d)){
					nn++;
					/// if (selection[nTile]) nn++; // count selected as 2 ?
				}
			}
			if (nn >= min_neibs) {
				new_sel[nTile] = true;
			}
		}
		for (int nTile =0; nTile <num_tiles; nTile++) {
			selection[nTile] |= new_sel[nTile];
		}
		return selection;
	}

	/**
	 * Select far tiles (even lone) over infinity areas that have sufficient disparity and strength
	 * @param min_far_strength minimal tile strength (may be 0 as the tiles are already filtered)
	 * @param min_far_disparity minimal tile disparity to accept lone strong tiles for far objects
	 * @param infinity_select tile, previously assigned to infinity (or null)
	 * @param selection existing selection (to add to) or null
	 * @param disparity array of tile disparity values (may have NaN-s)
	 * @param strength array of tile strength values
	 * @return new selection (new or added to existing)
	 */

	public boolean [] selectLoneFar(
			double     min_far_strength,
			double     min_far_disparity,
			boolean [] infinity_select,
			boolean [] selection,
			double []  disparity,
			double []  strength) {
		int num_tiles = disparity.length;
		if (infinity_select == null) {
			infinity_select = new boolean [num_tiles];
			for (int nTile = 0; nTile < num_tiles; nTile++) infinity_select[nTile] = true;
		}
		if (selection == null) selection = new boolean [num_tiles];
		for (int nTile =0; nTile <num_tiles; nTile++) {
			if (infinity_select[nTile] && (strength[nTile] > 0.0) && (strength[nTile] >= min_far_strength) && (disparity[nTile] >= min_far_disparity)) {
				selection[nTile] = true;
			}
		}
		return selection;
	}

	public double [][] filterDisparityStrength(
			double  min_disparity,
			double [] disparity,
			double [] strength,
			boolean [] selected){
		int num_tiles = disparity.length;
		double [][] fds = new double [2][num_tiles];
		for (int nTile = 0; nTile < num_tiles; nTile++) {
			fds[0][nTile] = Double.NaN;
			if ((strength[nTile] > 0.0) && (disparity[nTile] >= min_disparity) && ((selected == null) || selected[nTile])) {
				fds[0][nTile] = disparity[nTile];
				fds[1][nTile] = strength[nTile];
			}
		}
		return fds;
	}


	/**
	 * Filter selection by expanding, then shrinking (fills small gaps) shrinking (combined with
	 * previous step) and expanding again (removes small clusters)
	 * @param pre_expand number of steps for initial expand (odd last is - only hor/vert, even - last step includes diagonals)
	 * @param shrink number of shrink steps normally equals to pre_expand+post_expand
	 * @param post_expand numaber of final expansion steps
	 * @param selection selection to be modified
	 * @param prohibit tiles to keep from expansion/shrinking (such as infinity selection to modify only near regions (may be null)
	 */
	public void expandShrinkExpandFilter(
			int        pre_expand,
			int        shrink,
			int        post_expand,
			boolean [] selection,
			boolean [] prohibit)
	{
		tnImage.growSelection(
				pre_expand,  // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				selection,   // boolean [] tiles,
				prohibit);   // boolean [] prohibit)
		tnImage.shrinkSelection(
				shrink,      // int        shrink,
				selection,   // boolean [] tiles,
				prohibit);   // boolean [] prohibit)
		tnImage.growSelection(
				post_expand,  // int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				selection,   // boolean [] tiles,
				prohibit);   // boolean [] prohibit)
	}
	/**
	 * Combine (OR) two selections
	 * @param selection1 first selection
	 * @param selection2 second selection
	 * @param enable only use second selection if corresponding element in enable is true. May be null
	 * @param invert_enable make 'enable' 'disable'
	 * @return OR-ed selections
	 */
	public boolean [] combineSelections(
			boolean [] selection1,
			boolean [] selection2,
			boolean [] enable,
			boolean    invert_enable) {
		boolean [] selection = selection1.clone();
		if (enable == null) {
			for (int i = 0; i < selection.length; i++) selection[i] |= selection2[i];
		} else if (invert_enable) {
			for (int i = 0; i < selection.length; i++) selection[i] |= (selection2[i] && !enable[i]);
		} else {
			for (int i = 0; i < selection.length; i++) selection[i] |= (selection2[i] && enable[i]);
		}
		return selection;
	}



	public boolean [] selectFarObjects(
			double     strength_floor,
			double     min_disparity,
			double     min_mean,
			double     max_disparity,
			double     max_mean,
			double     min_disp_to_rms,
			double     min_strength,
			int        neib_dist, // >=1
			double     rsigma,
			double     tile_frac,
			boolean [] infinity_select,
			double []  disparity,
			double []  strength,
			int        tilesX,
			int        debugLevel) {
		int num_tiles = disparity.length;
		if (infinity_select == null) {
			infinity_select = new boolean [num_tiles];
			for (int nTile = 0; nTile < num_tiles; nTile++) infinity_select[nTile] = true;
		}
		if (Double.isNaN(strength_floor)) {
			strength_floor = Double.NaN;
			for (int nTile = 0; nTile < num_tiles; nTile++) if (
					infinity_select[nTile] &&
					!Double.isNaN(disparity[nTile])&&
					(strength[nTile] > 0) &&
					(Double.isNaN(strength_floor) || (strength[nTile] < strength_floor))){
				strength_floor = strength[nTile];
			}
		}
		if (debugLevel >-2) {
			System.out.println("selectFarObjects(): strength_floor = "+strength_floor);
		}
		if (debugLevel > 0) {
			double [] rsigmas = {0.0, 1.0, 0.5};
			for (double rs:rsigmas) {
				for (int nd = 1; nd < 7; nd++) {
					double [][] weight_mask =  get9weightMasks(
							nd, // >=1int        neib_dist, // >=1
							rs);
					(new ShowDoubleFloatArrays()).showArrays(
							weight_mask,
							2*nd+1,
							2*nd+1,
							true,
							"w9m_"+rs);

				}
			}
		}

		double [][] weight_mask =  get9weightMasks(
				neib_dist, // >=1int        neib_dist, // >=1
				rsigma); // int debugLevel
		int [] min_tiles = new int [weight_mask.length];
		for (int n = 0; n< min_tiles.length; n++) {
			int numnz = 0;
			for (int i = 0; i < weight_mask[n].length; i++) {
				if (weight_mask[n][i] >0.0) numnz++;

			}
			min_tiles[n] = (int) Math.round(numnz * tile_frac);

		}
		int sample_size = 2*neib_dist + 1;
//		int last = 2*neib_dist; tile_frac

		double [][] mean_to_rms = (debugLevel >-2) ?(new double [4][num_tiles]): null;
		double [][] dbg_dirs =    (debugLevel > 0) ?(new double [3*weight_mask.length][num_tiles]): null;
		if (dbg_dirs != null) {
			for (int n = 0; n < dbg_dirs.length; n++) {
				for (int nTile =0; nTile <num_tiles; nTile++) {
					dbg_dirs[n][nTile] = Double.NaN;
				}
			}
		}
		boolean [] selection = new boolean [num_tiles];
		if (mean_to_rms != null) for (int nTile =0; nTile <num_tiles; nTile++) for (int i = 0; i < mean_to_rms.length; i++) mean_to_rms[i][nTile] = Double.NaN; // only needed for debug?
		int    [] num_best = new int [weight_mask.length]; //debug feature

		for (int nTile =0; nTile <num_tiles; nTile++) {
			// disparity[nTile] may be Double.NaN
			if (infinity_select[nTile] && (disparity[nTile] >= min_disparity) && (disparity[nTile] <= max_disparity) && (strength[nTile] >= min_strength)) {
				double d0 = disparity[nTile];
				double [] s0 =  new double[weight_mask.length];
				double [] s0a = new double[weight_mask.length];
				double [] s1 =  new double[weight_mask.length];
				double [] s2 =  new double[weight_mask.length];
				int [] num_nz = new int[weight_mask.length];
				// calculate weights, means, and rms for each of 9 windows
				for (int dy = -neib_dist; dy <= neib_dist; dy++) {
					for (int dx = -neib_dist; dx <= neib_dist; dx++) {
						int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
						if (    (nTile1 >=0) &&
								!Double.isNaN(disparity[nTile1])) {
							double w = strength[nTile1] - strength_floor;
							if (w > 0.0) {
								int tindx = (dy + neib_dist) * sample_size + (dx + neib_dist);
								double d = disparity[nTile1];
								for (int n = 0; n < weight_mask.length; n++) {
									double ww = w * weight_mask[n][tindx];
									s0[n] += ww;
									s1[n] += ww* d;
									if ((dy!=0) || (dx!=0)){
										double dd = (d - d0);
										s0a[n] +=ww;
										s2[n] += ww*dd*dd;
									}

									if (ww >0.0) {
										num_nz[n]++;
									}
								}
							}
						}
					}
				}
				//min_tiles[n]
				double rm = 0.0;
				int best_dir = -1;
				for (int n = 0; n < weight_mask.length; n++) {
					if (s0[n] > 0) {
						s1[n] /= s0[n]; // mean
					}
					if (num_nz[n] < min_tiles[n]) {
						s1[n] =0.0;
					}
					if (s0a[n] > 0) {
						s2[n] = Math.sqrt(s2[n]/s0a[n]);
						double r = s1[n]/s2[n]; // denominator - rms (not including itself
						if (r > rm) {
							rm = r;
							best_dir = n;
						}
					}
				}
				if (dbg_dirs != null) {
					for (int n = 0; n < weight_mask.length; n++) {
						dbg_dirs[n + 0 * weight_mask.length][nTile] = s1[n]/s2[n];
						dbg_dirs[n + 1 * weight_mask.length][nTile] = s1[n];
						dbg_dirs[n + 2 * weight_mask.length][nTile] = s2[n];
					}
				}

//		double [][] dbg_dirs =    (debugLevel >-2) ?(new double [3*weight_mask.length][num_tiles]): null;

				if ((rm > min_disp_to_rms) && (s1[best_dir] >= min_mean ) && (s1[best_dir] <= max_mean )){
					num_best[best_dir]++;
					if (mean_to_rms != null) {
						mean_to_rms[0][nTile] = rm;
						if (best_dir >= 0) {
							mean_to_rms[1][nTile] = s1[best_dir];
							mean_to_rms[2][nTile] = s2[best_dir];
							mean_to_rms[3][nTile] = 0.1 * best_dir;
						}
					}
					selection[nTile] = true;
				}
			}
		}
//			double     min_mean,
//		double     max_disparity,
//		double     max_mean,
		if (debugLevel>-2) {
			System.out.println("selectFarObjects(): number of wins:");
			for (int i = 0; i < num_best.length; i++) {
				System.out.println(i+": "+num_best[i]);
			}

		}
		if (debugLevel > 0) {
			(new ShowDoubleFloatArrays()).showArrays(
					weight_mask,
					2*neib_dist+1,
					2*neib_dist+1,
					true,
					"weight_mask");
		}

		if (mean_to_rms != null) {
			double [][] dbg_img = {disparity, strength, mean_to_rms[0], mean_to_rms[1], mean_to_rms[2], mean_to_rms[3]};
			String [] titles = {"disparity", "strength", "mean_to_rms","mean","rms", "dir/10"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					mean_to_rms[0].length/tilesX,
					true,
					"mean_to_rms_"+neib_dist,
					titles);
		}
		if (dbg_dirs != null) {
			String [] titles = {
					"dr0","dr1","dr2","dr3","dr4","dr5","dr6","dr7","dr8",
					"m0","m1","m2","m3","m4","m5","m6","m7","m8",
					"r0","r1","r2","r3","r4","r5","r6","r7","r8"};
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_dirs,
					tilesX,
					mean_to_rms[0].length/tilesX,
					true,
					"mean_and_rms_"+neib_dist,
					titles);
		}
		return selection;
	}



	public int removeFalseMatches(
			double [][] disparity_bimap,
			double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
			double      stray_rstrength, // = 2.0;    // Relative to trusted strength - trust above that
			double      strength_rfloor, // = 0.28;    // Relative to trusted strength
			int         stray_dist,      // = 2;      // How far to look (and erase) around a potentially falsely matched tile
			double      stray_over,      // = 2.0;    // Stray tile should be this stronger than the strongest neighbor to be recognized
			int         debugLevel)
	{
			double max_stray_strength = trusted_strength * stray_rstrength;
			double strength_floor =     trusted_strength * strength_rfloor;
			double max_under = 1.0/stray_over;
			ArrayList<Integer> stray_list = new ArrayList<Integer>();
			final double [] disparity = disparity_bimap[ImageDtt.BI_TARGET_INDEX];
			final double [] strength =  disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX];
			for (int nTile=0; nTile < strength.length; nTile++) {
				if (!Double.isNaN(disparity[nTile]) && (strength[nTile] >= trusted_strength) && (strength[nTile] < max_stray_strength)){
					double max_stength = 0.0;
					for (int dy = -stray_dist; dy <= stray_dist; dy++) {
						for (int dx = -stray_dist; dx <= stray_dist; dx++)  if ((dy!=0) || (dx!=0)){
							int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
							if (    (nTile1 >=0) &&
									!Double.isNaN(disparity[nTile1]) &&
									(strength[nTile1] > max_stength)) {
								max_stength = strength[nTile1];
							}
						}
					}
					if (((max_stength - strength_floor)/(strength[nTile] - strength_floor) <= max_under) &&
							(max_stength < trusted_strength)){
						stray_list.add(nTile);
						if (debugLevel > -2) {
							System.out.println(String.format("removeFalseMatches(): adding tile %d (%d/%d) ",nTile, nTile%tnImage.sizeX, nTile/tnImage.sizeX));
						}
					}
				}
			}
			// now erase around tiles in stray_list
			for (int nTile: stray_list) {
				for (int dy = -stray_dist; dy <= stray_dist; dy++) {
					for (int dx = -stray_dist; dx <= stray_dist; dx++) {
						int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
						disparity[nTile1] = Double.NaN;
						strength[nTile1] = 0.0;
					}
				}
			}
			return stray_list.size();

	}

	public boolean [] getLTTrusted(
			double [][] disparity_bimap,
			double      min_disparity, //  =         0.0;    // apply low texture to near objects
			double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
			double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
			double      friends_diff, // =           0.15;   // pix difference to neighbors to be considered a match (TODO: use tilted)
			double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
			int         min_friends_any, // =        2;      // minimal number of even weak friends
			int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
			int         friends_dist, // =           3;      // how far to look for friends
			int         debugLevel
			) {
		final double cond_strength = trusted_strength * need_friends;
		final double [] disparity = disparity_bimap[ImageDtt.BI_TARGET_INDEX];
		final double [] strength =  disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX];
		final boolean [] trusted = new boolean[strength.length];
		final boolean [] cond_trusted = trusted.clone();
		for (int nTile = 0; nTile < trusted.length; nTile ++) if (strength[nTile] >= cond_strength){
			cond_trusted[nTile] = true;
			trusted[nTile] = strength[nTile] >= trusted_strength;
		}
		for (int new_tiles = 0;  ; new_tiles = 0) {
			for (int nTile = 0; nTile < trusted.length; nTile ++) {
				if (cond_trusted[nTile] && !trusted[nTile]) {
					int num_trusted = 0;
					int num_friends = 0;
					double low_friend =  disparity[nTile] - friends_diff -  friends_rdiff*disparity[nTile];
					double high_friend = disparity[nTile] + friends_diff +  friends_rdiff*disparity[nTile];
					label_tile:
					{
						for (int dy = -friends_dist; dy <= friends_dist; dy++) {
							for (int dx = -friends_dist; dx <= friends_dist; dx++) if ((dy!=0) || (dx!=0)){
								int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
								if ((nTile1 >= 0) && (disparity[nTile1] >= low_friend)  && (disparity[nTile1] <= high_friend)){ // disparity[nTile1] may be NaN!
									if (cond_trusted[nTile1]) {
										num_friends++;
										if (num_friends >= min_friends_any){
											trusted[nTile] = true;
											new_tiles++;
											break label_tile;
										} else if (trusted[nTile1]) {
											num_trusted++;
											if (num_trusted >= min_friends_trusted){
												trusted[nTile] = true;
												new_tiles++;
												break label_tile;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			if (debugLevel > -2) System.out.print ("new tiles = "+new_tiles); // find out why second pass always returns 0
			if (new_tiles == 0) break;
		}
		if (debugLevel > -2) System.out.println();
		return trusted;
	}

	public int removeLTUntrusted(
			double [][] disparity_bimap,
			double      min_disparity, //  =         0.0;    // apply low texture to near objects
			double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
			double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
			double      friends_diff, // =           0.15;   // pix difference to neighbors to be considered a match (TODO: use tilted)
			double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
			int         min_friends_any, // =        2;      // minimal number of even weak friends
			int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
			int         friends_dist, // =           3;      // how far to look for friends
			int         debugLevel
			) {
		boolean [] trusted = 	getLTTrusted(
				disparity_bimap,     // double [][] disparity_bimap,
				min_disparity,       //double      min_disparity, //  =         0.0;    // apply low texture to near objects
				trusted_strength,    //double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
				need_friends,        //double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
				friends_diff,        //double      friends_diff, // =           0.15;   // pix difference to neighbors to be considered a match (TODO: use tilted)
				friends_rdiff,       //double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
				min_friends_any,     //int         min_friends_any, // =        2;      // minimal number of even weak friends
				min_friends_trusted, //int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
				friends_dist,        //int         friends_dist, // =           3;      // how far to look for friends
				debugLevel);         //int         debugLevel
		int num_trusted = 0;
		for (int nTile = 0; nTile < trusted.length; nTile++) {
			if (trusted[nTile]) {
				num_trusted++;
			} else {
				disparity_bimap[ImageDtt.BI_TARGET_INDEX][nTile] = Double.NaN;
				disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX][nTile] = 0.0;
			}
		}
		return num_trusted;
	}




	public double [] suggestLTTiles(
			double [][] disparity_bimap,
			boolean []  trusted,       // may be null if disparity is alreasdy NaN-ed
			double      min_disparity, //  =         0.0;    // apply low texture to near objects
			double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
			double      strength_rfloor,  // =       0.28;   // strength floor relative to trusted_strength
			double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
			int         extend_dist, // =            3;      // how far to extend around known tiles (probably should increase this value up to?
			// dealing with neighbors variance
			double      wsigma,     //  = 1.0; // influence of far neighbors diminish as a Gaussian with this sigma
			double      max_asigma, // =             .15;     // Maximal acceptable standard deviation of the neighbors (remove, then add)
			double      max_rsigma, // =             .05;     // Maximal acceptable standard deviation of the neighbors (remove, then add)
			int         debugLevel
			) {
		final double rsigma = max_rsigma; //pix/pix
		final double asigma = max_asigma; // 1.0;

		final double [][] weights = new double [extend_dist+1][extend_dist+1];
		//final double cond_strength = trusted_strength * need_friends;
		final double strength_floor = strength_rfloor * trusted_strength;
		for (int i = 0; i <weights.length; i++) {
			for (int j = i; j <weights[i].length; j++) {
				weights[i][j]=Math.exp(-(i*i+j*j)/(2*wsigma*wsigma));
				weights[j][i] = weights[i][j];
			}
		}
		final double [] disparity = disparity_bimap[ImageDtt.BI_TARGET_INDEX];
		final double [] strength =  disparity_bimap[ImageDtt.BI_STR_CROSS_INDEX];
		if (trusted == null) {
			trusted = new boolean[disparity.length];
			for (int nTile = 0; nTile < trusted.length; nTile++) trusted[nTile] =  !Double.isNaN(disparity[nTile]);
		}

		double sigma = asigma;
		double sigma2 = sigma*sigma;
		final double [] to_measure = new double [disparity.length];
		for (int nTile = 0; nTile < to_measure.length; nTile++) {
			to_measure[nTile] = Double.NaN;
		}
		final boolean [] candidates = new boolean[strength.length];
		// can be done faster if jump on empty by square side
		for (int nTile = 0; nTile < candidates.length; nTile++) if (!trusted[nTile]){
			label_tile1:
			for (int dy = -extend_dist; dy <= extend_dist; dy++) {
				for (int dx = -extend_dist; dx <= extend_dist; dx++) if ((dy!=0) ||(dx !=0)){
					int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
					if ((nTile1 >= 0) && trusted[nTile1]){
						candidates[nTile] = true;
						break label_tile1;
					}
				}
			}
		}

		int num_sigma = 0;
		for (int nTile = 0; nTile < candidates.length; nTile++) if (candidates[nTile]){
			ArrayList<Integer> sample_list = new ArrayList<Integer>();
			double s0 = 0.0, s1 = 0.0, s2 = 0.0;
			for (int dy = -extend_dist; dy <= extend_dist; dy++) {
				for (int dx = -extend_dist; dx <= extend_dist; dx++) if ((dy!=0) ||(dx !=0)){
					int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
					if ((nTile1 >= 0) && trusted[nTile1]){
						double w = (strength[nTile1] - strength_floor) * weights[(dy>0)?dy:-dy][(dx>0)?dx:-dx];
						s0 += w;
						s1 += w * disparity[nTile1];
						s2 += w * disparity[nTile1] * disparity[nTile1];
						sample_list.add(nTile1);
					}
				}
			}
			if (s0 <= 0) {
				System.out.println("suggestLTTiles() BUG? nTile = "+nTile+" s0 ="+s0);
				continue;
			}
			double s_mean = s1/s0;
			double smpl_var = s2/s0 -s_mean*s_mean;
			sigma = asigma + rsigma * s_mean;
			sigma2 = sigma*sigma;
// FIXME: use tilted planes
			if (smpl_var > sigma2) {
				if (debugLevel > -2) {
					System.out.print (String.format("%3d/%3d mean=%8f sigma2=%f var=%8f tiles=%3d ",nTile%tnImage.sizeX, nTile/tnImage.sizeX, s_mean, sigma2, smpl_var, sample_list.size()));
				}
				num_sigma++;
				ArrayList<Integer> rejected_list = new ArrayList<Integer>();
				int [] xy0 = tnImage.getXY(nTile);
				while (smpl_var > sigma2) {
					// find worst sample
					int worst_indx = -1;
					double worst_diff = 0;
					double d;
					for (int i = 0; i < sample_list.size(); i++) {
						int nTile1 = sample_list.get(i);
						d = Math.abs(disparity[nTile1] - s_mean);
						if (d > worst_diff) {
							worst_diff = d;
							worst_indx = i;
						}
					}
					// remove worst sample, add to reject list
					int nTile1 = sample_list.get(worst_indx);
					rejected_list.add(nTile1);
					sample_list.remove(worst_indx);
					// recalculate statistics
					int [] xy = tnImage.getXY(nTile1);
					int dx =xy[0] - xy0[0];
					int dy =xy[1] - xy0[1];
					double w = (strength[nTile1] - strength_floor) * weights[(dy>0)?dy:-dy][(dx>0)?dx:-dx];
					s0 -= w;
					s1 -= w * disparity[nTile1];
					s2 -= w * disparity[nTile1] * disparity[nTile1];
					s_mean = s1/s0;
					smpl_var = s2/s0 -s_mean*s_mean;
				}
				if (debugLevel > -2) {
//					System.out.print (" -> s_mean = "+s_mean+", smpl_var = "+smpl_var+" ... "+ " ntiles="+(sample_list.size()));
					System.out.print (String.format(" (-)-> mean=%8f var=%8f tiles=%3d ",s_mean, smpl_var, sample_list.size()));

				}

				// Try to add best of the rejected back (trying to deal with 2-maximums histogram)
				while (smpl_var < sigma2) { // then remove last added
					// find best rejected sample
					int best_indx = -1;
					double best_diff = 0;
					double d;
					for (int i = 0; i < rejected_list.size(); i++) {
						int nTile1 = rejected_list.get(i);
						d = Math.abs(disparity[nTile1] - s_mean);
						if ((best_indx <0) || (d < best_diff)) {
							best_diff = d;
							best_indx = i;
						}
					}
					// restore best rejected sample
					int nTile1 = rejected_list.remove(best_indx);
					sample_list.add(nTile1); // best_indx); // will be last - easy to remove
					// recalculate statistics
					int [] xy = tnImage.getXY(nTile1);
					int dx =xy[0] - xy0[0];
					int dy =xy[1] - xy0[1];
					double w = (strength[nTile1] - strength_floor) * weights[(dy>0)?dy:-dy][(dx>0)?dx:-dx];
					s0 += w;
					s1 += w * disparity[nTile1];
					s2 += w * disparity[nTile1] * disparity[nTile1];
					s_mean = s1/s0;
					smpl_var = s2/s0 -s_mean*s_mean;
				}
				if (debugLevel > -2) {
					System.out.print (String.format(" (+)-> mean=%8f var=%8f tiles=%3d ",s_mean, smpl_var, sample_list.size()));
				}
				// remove last added sample
				// remove worst sample, add to reject list
				int nTile1 = sample_list.get(sample_list.size()-1); // last added, no need to actually remove
				// recalculate statistics
				int [] xy = tnImage.getXY(nTile1);
				int dx =xy[0] - xy0[0];
				int dy =xy[1] - xy0[1];
				double w = (strength[nTile1] - strength_floor) * weights[(dy>0)?dy:-dy][(dx>0)?dx:-dx];
				s0 -= w;
				s1 -= w * disparity[nTile1];
				s2 -= w * disparity[nTile1] * disparity[nTile1];
				s_mean = s1/s0;
				smpl_var = s2/s0 -s_mean*s_mean;

				if (debugLevel > -2) {
//					System.out.println (" s_mean = "+s_mean+", smpl_var = "+smpl_var+ " ntiles="+(sample_list.size()-1) );
					System.out.println (String.format(" => mean=%8f var=%8f tiles=%3d ", s_mean, smpl_var, sample_list.size()-1));
				}

			} // if (smpl_var > sigma2) {
			to_measure[nTile] = s_mean;
			if (debugLevel > -2) {
				disparity_bimap[ImageDtt.BI_DBG2_INDEX][nTile] = s_mean;
				disparity_bimap[ImageDtt.BI_DBG3_INDEX][nTile] = Math.sqrt(smpl_var);
			}
		}

		if (debugLevel > -2) {
			int num_to_measure =0;
			for (int nTile = 0; nTile < to_measure.length; nTile++) {
				if (!Double.isNaN(to_measure[nTile]))num_to_measure++;
			}

			System.out.println ("Updated sigma in tiles:"+num_sigma+" (sigma = "+sigma+", sigma2 = "+sigma2);
			System.out.println ("Tiles to meaure:"+num_to_measure);
			disparity_bimap[ImageDtt.BI_DBG3_INDEX] = to_measure; // overwrites old debug data
			(new ShowDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tnImage.sizeX,
					tnImage.sizeY,
					true,
					"BiCamDSI-2",
					ImageDtt.BIDISPARITY_TITLES);
		}
		return to_measure;
	}
/*
 * sigma2
	  static int  BI_DISP_FULL_INDEX =            0;  // 0 - disparity for all directions of the main camera
	  static int  BI_DISP_HOR_INDEX =             1;  // 1 - disparity for 2 horizontal pairs of the main camera
	  static int  BI_DISP_VERT_INDEX =            2;  // 2 - disparity for 2 vertical pairs of the main camera
	  static int  BI_DISP_DIAGM_INDEX =           3;  // 3 - disparity for main diagonal pair of the main camera
	  static int  BI_DISP_DIAGO_INDEX =           4;  // 4 - disparity for main diagonal pair of the main camera
	  static int  BI_ADISP_FULL_INDEX =           5;  // 5 - disparity for all directions of the aux camera
	  static int  BI_ADISP_HOR_INDEX =            6;  // 6 - disparity for 2 horizontal pairs of the aux camera
	  static int  BI_ADISP_VERT_INDEX =           7;  // 7 - disparity for 2 vertical pairs of the aux camera
	  static int  BI_ADISP_DIAGM_INDEX =          8;  // 8 - disparity for main diagonal pair of the aux camera
	  static int  BI_ADISP_DIAGO_INDEX =          9;  // 9 - disparity for main diagonal pair of the aux camera
	  static int  BI_DISP_CROSS_INDEX =          10;  //10 - disparity between the main the aux camera
	  static int  BI_DISP_CROSS_DX_INDEX =       11;  //11 - delta disparity between the main the aux camera (horizontal)
	  static int  BI_DISP_CROSS_DY_INDEX =       12;  //12 - delta disparity between the main the aux camera (vertical)
	  static int  BI_STR_FULL_INDEX =            13;  //13 - strength for all directions of the main camera
	  static int  BI_STR_HOR_INDEX =             14;  //14 - strength for 2 horizontal pairs of the main camera
	  static int  BI_STR_VERT_INDEX =            15;  //15 - strength for 2 vertical pairs of the main camera
	  static int  BI_STR_DIAGM_INDEX =           16;  //16 - strength for main diagonal pair of the main camera
	  static int  BI_STR_DIAGO_INDEX =           17;  //17 - strength for main diagonal pair of the main camera
	  static int  BI_ASTR_FULL_INDEX =           18;  //18 - strength for all directions of the aux camera
	  static int  BI_ASTR_HOR_INDEX =            19;  //19 - strength for 2 horizontal pairs of the aux camera
	  static int  BI_ASTR_VERT_INDEX =           20;  //20 - strength for 2 vertical pairs of the aux camera
	  static int  BI_ASTR_DIAGM_INDEX =          21;  //21 - strength for main diagonal pair of the aux camera
	  static int  BI_ASTR_DIAGO_INDEX =          22;  //22 - strength for main diagonal pair of the aux camera
	  static int  BI_STR_CROSS_INDEX =           23;  //23 - strength between the main the aux camera
	  static int  BI_STR_ALL_INDEX =             24;  //23 - average strength (product of strengths to 1/3 power), TODO: strength at cross disparity
	  static int  BI_TARGET_INDEX =              25;  //24 - target disparity
	  static int  BI_DBG1_INDEX =                26;  //26 - debug layer 1
	  static int  BI_DBG2_INDEX =                27;  //27 - debug layer 2

 */

}
