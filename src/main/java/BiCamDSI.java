/**
 ** BiCamDSI - Building DSI using correlation between two quad cameras
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TwoQuadCLT.java is free software: you can redistribute it and/or modify
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


public class BiCamDSI {
//	public int tilesX;
//	public int tilesY;
	TileNeibs tnImage; //  = new TileNeibs(tilesX, tilesY)

	public BiCamDSI(
			int tilesX,
			int tilesY) {
		tnImage  = new TileNeibs(tilesX, tilesY);
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
							for (int dx = -friends_dist; dx <= friends_dist; dx++) if ((dy!=0) ||(dx !=0)){
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
			System.out.println("new tiles = "+new_tiles); // find out why second pass always returns 0
			if (new_tiles == 0) break;
		}
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
			}
		}
		return num_trusted;
	}


	public double [] suggestLTTiles(
			double [][] disparity_bimap,
			boolean []  trusted,       // may be null if disparity is alreasdy NaN-ed
			double      min_disparity, //  =         0.0;    // apply low texture to near objects
			double      trusted_strength, // =       0.2;    // strength sufficient without neighbors
			double      need_friends, // =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
//			double      friends_diff, // =           0.15;   // pix difference to neighbors to be considered a match (TODO: use tilted)
//			double      friends_rdiff, // =          0.04;   // additional relative pix per pixel of disparity
//			int         min_friends_any, // =        2;      // minimal number of even weak friends
//			int         min_friends_trusted, // =    2;      // minimal number of trusted (strong or already confirmed)
//			int         friends_dist, // =           3;      // how far to look for friends
//			boolean     replace_lone, // =           true;   // try to overwrite lone weak
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
		final double cond_strength = trusted_strength * need_friends;
		final double strength_floor = 0.7*cond_strength;
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
//		final boolean [] trusted = new boolean[strength.length];
//		final boolean [] cond_trusted = trusted.clone();
//		for (int nTile = 0; nTile < trusted.length; nTile ++) if (strength[nTile] >= cond_strength){
//			cond_trusted[nTile] = true;
//			trusted[nTile] = strength[nTile] >= trusted_strength;
//		}
		double sigma = asigma;
		double sigma2 = sigma*sigma;
		final double [] to_measure = new double [disparity.length];
		for (int nTile = 0; nTile < to_measure.length; nTile++) {
			to_measure[nTile] = Double.NaN;
		}
/*
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
							for (int dx = -friends_dist; dx <= friends_dist; dx++) if ((dy!=0) ||(dx !=0)){
								int nTile1 =  tnImage.getNeibIndex(nTile, dx, dy);
								if ((nTile1 >= 0) && (disparity[nTile1] >= low_friend)  && (disparity[nTile1] <= high_friend)){
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
			System.out.println("new tiles = "+new_tiles);
			if (new_tiles == 0) break;
		}

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

		final boolean [] dbg_trusted = trusted.clone();
*/
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
/*
		if (debugLevel > -2) {
			for (int nTile = 0; nTile < trusted.length; nTile ++) {
				disparity_bimap[ImageDtt.BI_DBG1_INDEX][nTile] = Double.NaN;
				disparity_bimap[ImageDtt.BI_DBG2_INDEX][nTile] = Double.NaN;
				disparity_bimap[ImageDtt.BI_DBG3_INDEX][nTile] = Double.NaN;
				disparity_bimap[ImageDtt.BI_DBG4_INDEX][nTile] = Double.NaN;
				if (trusted[nTile]) {
					disparity_bimap[ImageDtt.BI_DBG1_INDEX][nTile] = disparity[nTile];
					disparity_bimap[ImageDtt.BI_DBG2_INDEX][nTile] = disparity[nTile];
					disparity_bimap[ImageDtt.BI_DBG4_INDEX][nTile] = 2.0;
//					if (dbg_trusted[nTile])disparity_bimap[ImageDtt.BI_DBG4_INDEX][nTile] = 3.0;
//				} else if (cond_trusted[nTile]) {
//					disparity_bimap[ImageDtt.BI_DBG4_INDEX][nTile] = 1.0;
				}
			}
			for (int nTile = 0; nTile < trusted.length; nTile ++) {
				if (candidates[nTile]) 	disparity_bimap[ImageDtt.BI_DBG4_INDEX][nTile] = 0.0;

			}
			(new showDoubleFloatArrays()).showArrays(
					disparity_bimap,
					tnImage.sizeX,
					tnImage.sizeY,
					true,
					"BiCamDSI",
					ImageDtt.BIDISPARITY_TITLES);
		}
*/
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
//			if (debugLevel > -2) {
//				disparity_bimap[ImageDtt.BI_DBG2_INDEX][nTile] = s_mean;
//				disparity_bimap[ImageDtt.BI_DBG3_INDEX][nTile] = Math.sqrt(smpl_var);
//			}
//			final double rsigma = 0.05; //pix/pix
//			final double asigma = max_sigma; // 1.0;
			sigma = asigma + rsigma * s_mean;
			sigma2 = sigma*sigma;
// FIXME: use tilted planes
			if (smpl_var > sigma2) {
				if (debugLevel > -2) {
//					System.out.print ((nTile%tnImage.sizeX)+"/"+(nTile/tnImage.sizeX)+": s_mean = "+s_mean+", smpl_var = "+smpl_var+" ... "+ " ntiles="+(sample_list.size()));
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
			(new showDoubleFloatArrays()).showArrays(
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
