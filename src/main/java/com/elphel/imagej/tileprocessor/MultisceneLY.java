package com.elphel.imagej.tileprocessor;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.gpu.GpuQuad;
import com.elphel.imagej.gpu.TpTask;

public class MultisceneLY {
	public static double [] ZERO3 = {0.0,0.0,0.0};
	public static double  LINE_ERR = 0.1;
	public int            threadsMax = 100;  // maximal number of threads to launch
	public boolean        updateStatus = true;
	public int            numSens;
    public static enum    MSLY_MODE {INF_ONLY, NOINF_ONLY, INF_NOINF};
    public static String [] SINF_NOINF= {"inf","noinf"};  	
	public MultisceneLY (
			int            numSens,
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus) {
			this.numSens =                numSens;
			this.threadsMax =             threadsMax;
			this.updateStatus =           updateStatus;
	}
	/**
	 * Remeasure scene (from data in the model directory) to find number of correlation maximums
	 * in each tile. Needed to keep only single-maximum tiles during Lazy Eye adjustment
	 * @param clt_parameters parameters
	 * @param lma_only keep only tiles that have LMA results
	 * @param scene scene data 
	 * @param debug_level debug level
	 * @return number of correlation maximums per tile (now 0,1, or 2) in linescan order
	 */
	public int [] getNumCorrMax(
			CLTParameters  clt_parameters,
			boolean lma_only,
			QuadCLT scene, // ordered by increasing timestamps
			int debug_level)
	{
		scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU)
		CLTPass3d scan = new CLTPass3d(scene.tp);
		double [] target_disparity = scene.getDSRBG()[0];
		boolean run_lma = true; // lma_only; // or always?
		scan.setTileOpDisparity(target_disparity); // measure w/o LMA, use just strength
		scene.CLTMeas( // perform single pass according to prepared tiles operations and disparity
				clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
				scan,           // final CLTPass3d   scan,
				false,          // final boolean     save_textures,
				false,          // final boolean     need_diffs,     // calculate diffs even if textures are not needed 
				0,              // final int         clust_radius,
				true,           // final boolean     save_corr,
				run_lma,          // final boolean     run_lma, // =    true;
				threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
				updateStatus,   // final boolean     updateStatus,
				debug_level);    // final int         debugLevel);
		int [] num_corr_max = scan.getNumTileMax();
		if (lma_only) {
			boolean [] has_lma = scan.hasLMADefined();
			for (int i = 0; i < has_lma.length; i++) {
				if (!has_lma[i]) {
					num_corr_max[i] = -1;
				}
			}
		}
		return num_corr_max;
	}
	
	/**
	 * Process multiscene combined depth map to extract infinity tiles
	 * by analyzing depth map histogram
	 * @param tp tile processor instance
	 * @param composite_ds [tile]{disparity, strength} - per-tile disparity/strength pair
	 * @param far_inf far limit for a tile considered to be at infinity (typical -0.5)
	 * @param near_inf near limit for infinity tiles (typical +0.5)
	 * @param far_fract mode - fraction of all pixels in a  histogram, typical 0.05.
	 *                  0.0 corresponds to exact minimum of all tiles far_inf<=disparity<=near_inf,
	 *                  1.0 - exact maximum. 0.5 - median.
	 * @param inf_range_offs - offset (add to) the mode value by this parameter to get average
	 *                  expected disparity of objects at infinity. The result is returned in
	 *                  inf_avg[0] that should be initialized to int[1] by the caller.
	 * @param inf_range full range (symmetrical around offset mode value) of potentially
	 *                  infinity tiles
	 * @param min_inf_str minimal strength of infinity tiles
	 * @param min_fg_str minimal strength of non-infinity tiles
	 * @param inf_avg (should be initialized as double[1] - return parameter of average disparity
	 *                  at infinity 
	 * @param debug     generate debug images if true.
	 * @return          per-tile array, where true is for infinity tiles
	 */
	public static boolean [] getComboInfinity(
			TileProcessor tp,
			double [][]   composite_ds,
			double        far_inf,
			double        near_inf,
			double        far_fract,
			double        inf_range_offs,
			double        inf_range,
			double        min_inf_str,
			double        min_fg_str,
			double []     inf_avg,
			boolean       debug
			) 	{
		int num_bins = 256;
	    int [] hist = new int [num_bins];
	    double [][] dbg_img = debug? (new double[2][]): null;
		double [] disparity_lma = composite_ds[OpticalFlow.COMBO_DSN_INDX_LMA]; //COMBO_DSN_INDX_DISP_BG?
		double [] strength = composite_ds[OpticalFlow.COMBO_DSN_INDX_STRENGTH];
		int num_good = 0;
		for (int nTile = 0; nTile< disparity_lma.length; nTile++) {
			if ((disparity_lma[nTile] >= far_inf) && (disparity_lma[nTile] < near_inf) && (strength[nTile] >= min_inf_str)) {
				int ibin = (int) ((disparity_lma[nTile] - far_inf)/(near_inf - far_inf) * num_bins);
				if ((ibin >= 0) && (ibin < num_bins)) {
					hist[ibin]++;
					num_good++;
				}
			}
		}
		int ifar_cumul = (int) (far_fract * num_good);
		double cumul = 0.0;
		double bin = num_bins;
		double cumul_prev = 0.0;
		for (int i = 0; i < num_bins; i++) {
			cumul += hist[i];
			if (cumul >= ifar_cumul) {
				bin = i + (ifar_cumul - cumul_prev) / (cumul - cumul_prev);
				break;
			}
			cumul_prev = cumul;
		}
		double far_lim = far_inf + bin * (near_inf - far_inf) / num_bins - inf_range_offs;
		double near_lim = far_lim + inf_range;
		boolean [] is_inf = new boolean[disparity_lma.length];
		for (int nTile = 0; nTile< is_inf.length; nTile++) {
			is_inf[nTile] = ((disparity_lma[nTile] >= far_lim) &&
					(disparity_lma[nTile] < near_lim) &&
					(strength[nTile] >= min_inf_str));
		}		
		boolean [] is_fg = new boolean[disparity_lma.length];
		for (int nTile = 0; nTile< is_fg.length; nTile++) {
			is_fg[nTile] = ((disparity_lma[nTile] >= near_lim) &&
					(strength[nTile] >= min_fg_str));
		}	
		
		if (inf_avg != null) {
			double s = 0.0;
			int num = 0;
			for (int nTile = 0; nTile < is_inf.length; nTile++) {
				if (is_inf[nTile]) {
					s+= disparity_lma[nTile];
					num++;
				}
			}
			inf_avg[0] = s/num;
		}
		
		if (dbg_img != null) {
			for (int i = 0; i < dbg_img.length; i++) {
				dbg_img[i] = new double [is_inf.length];
			}
			for (int nTile = 0; nTile < is_inf.length; nTile++) {
				if (is_inf[nTile]) dbg_img[0][nTile] += 1.0;
				if (is_fg[nTile])  dbg_img[1][nTile] += 1.0;
			}
		}
		// grow sure fg
		tp.growTiles(
				2,      // int        grow,   // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				is_fg,  // boolean [] tiles,
				null);  // boolean [] prohibit)
		// Remove even core infinity if it is expanden FG
		for (int nTile = 0; nTile < is_inf.length; nTile++) {
			is_inf[nTile] &= !is_fg[nTile];
		}

		tp.growTiles(
				2,      // int        grow,   // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				is_inf, // boolean [] tiles,
				is_fg); // boolean [] prohibit)
		
		tp.fillGaps( // grows, then shrinks
				2, // int depth, // same as grow - odd - 4 directions, even - 8
				true, // false, // boolean    poison, // do not fill gaps that even touch prohibited
				is_inf, // boolean [] tiles,
				is_fg); //boolean [] prohibit)
		
		if (dbg_img != null) {
			for (int nTile = 0; nTile < is_inf.length; nTile++) {
				if (is_inf[nTile]) dbg_img[0][nTile] += 2.0;
				if (is_fg[nTile])  dbg_img[1][nTile] += 2.0;
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tp.getTilesX(),
					tp.getTilesY(),
					true,
					"inf_fg",
					new String[] {"infinity", "foreground"});
		}
		return is_inf;
	}
	
	/**
	 * Map reference scene (the last one in a sequence) disparity map and map of
	 * the infinity tiles to all the scenes.  
	 * @param scenes sequence of scenes
	 * @param inf_disp_ref expected disparity at infinity
	 * @param infinity_ref map of infinity tiles for the reference (last) scene
	 * @param debug show debug images
	 * @param threadsMax maximal number of threads to use
	 * @return boolean array per tile - which tiles are infinity - per scene, per tile
	 */
	public static boolean [][] infinityPerScene(
			QuadCLT [] scenes,
			double     inf_disp_ref, // average disparity at infinity for ref scene
			boolean [] infinity_ref,
			boolean    debug,
			int        threadsMax
			){
		final double [] disparity_ref = new double [infinity_ref.length];
		final int num_scenes = scenes.length;
		int last_scene_index = scenes.length-1;
		QuadCLT last_scene = scenes[last_scene_index];
		final ErsCorrection ers_reference = last_scene.getErsCorrection();
		final int tilesX = last_scene.getTileProcessor().getTilesX();
		final int tilesY = last_scene.getTileProcessor().getTilesY();
		final int tileSize = last_scene.getTileProcessor().getTileSize();
		boolean [][] inf_scenes = new boolean [scenes.length][tilesX*tilesY];
		for (int nTile = 0; nTile < disparity_ref.length; nTile++) {
			disparity_ref[nTile] = infinity_ref[nTile]? inf_disp_ref : (Double.NaN);
		}
		TileNeibs tn = new TileNeibs(tilesX,tilesY);
		final double [][][] scenes_pXpYD = new double [scenes.length][][];
		final int dbg_tile = debug? 320 : -1;
		for (int nscene = 0; nscene < num_scenes; nscene++) {
			String ts = scenes[nscene].getImageName();
			if (nscene == last_scene_index) {
				scenes_pXpYD[nscene] = OpticalFlow.transformToScenePxPyD(
						null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
						disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
						OpticalFlow.ZERO3,              // final double []   scene_xyz, // camera center in world coordinates
						OpticalFlow.ZERO3,              // final double []   scene_atr, // camera orientation relative to world frame
						last_scene,          // final QuadCLT     scene_QuadClt,
						last_scene,          // final QuadCLT     reference_QuadClt)
						threadsMax);
			} else {
				double []   scene_xyz = ers_reference.getSceneXYZ(ts);
				double []   scene_atr = ers_reference.getSceneATR(ts);
				double []   scene_ers_xyz_dt = ers_reference.getSceneErsXYZ_dt(ts);
				double []   scene_ers_atr_dt = ers_reference.getSceneErsATR_dt(ts);
				scenes[nscene].getErsCorrection().setErsDt(
						scene_ers_xyz_dt, // double []    ers_xyz_dt,
						scene_ers_atr_dt); // double []    ers_atr_dt)(ers_scene_original_xyz_dt);
				//setupERS() will be inside transformToScenePxPyD()
				scenes_pXpYD[nscene] = OpticalFlow.transformToScenePxPyD( // will be null for disparity == NaN, total size - tilesX*tilesY
						null, // final Rectangle [] extra_woi,    // show larger than sensor WOI (or null)
						disparity_ref,      // final double []   disparity_ref, // invalid tiles - NaN in disparity (maybe it should not be masked by margins?)
						scene_xyz,          // final double []   scene_xyz, // camera center in world coordinates
						scene_atr,          // final double []   scene_atr, // camera orientation relative to world frame
						scenes[nscene],      // final QuadCLT     scene_QuadClt,
						last_scene,          // final QuadCLT     reference_QuadClt)
						threadsMax);
			}
		}
		for (int nTile = 0; nTile < infinity_ref.length; nTile++) if (infinity_ref[nTile]){
			if (nTile == dbg_tile) {
				System.out.println("infinityPerScene(): nTile = "+nTile);
			}
			// mark tile itself
			for (int nscene = 0; nscene < num_scenes; nscene++) {
				if (scenes_pXpYD[nscene][nTile] != null) {
					int [] pXY0 = {
							(int) Math.floor(scenes_pXpYD[nscene][nTile][0]/tileSize),
							(int) Math.floor(scenes_pXpYD[nscene][nTile][1]/tileSize)};
					int indx0 = tn.getIndex(pXY0);
					if (indx0 >= 0) {
						inf_scenes[nscene][indx0] = true;
					}
				}
			}
			for (int dir = 0; dir < 4; dir++) {// Only half of all directions
				int nTile1 = tn.getNeibIndex(nTile, dir);
				if ((nTile1 >=0) && infinity_ref[nTile1]) {
					for (int nscene = 0; nscene < num_scenes; nscene++) {
						if ((scenes_pXpYD[nscene][nTile] != null) && (scenes_pXpYD[nscene][nTile1] != null)) {
							int [] pXY0 = {
									(int) Math.floor(scenes_pXpYD[nscene][nTile][0]/tileSize),
									(int) Math.floor(scenes_pXpYD[nscene][nTile][1]/tileSize)};
							int [] pXY1 = {
									(int) Math.floor(scenes_pXpYD[nscene][nTile1][0]/tileSize),
									(int) Math.floor(scenes_pXpYD[nscene][nTile1][1]/tileSize)};
							int indx0 = tn.getIndex(pXY0);
							int indx1 = tn.getIndex(pXY1);
							if ((indx0 >= 0) && (indx1 >= 0)) {
								int [] pXY2 = { // half-way
									(int) Math.round(0.5 * (pXY0[0] + pXY1[0])),
									(int) Math.round(0.5 * (pXY0[1] + pXY1[1]))};
								int indx2 = tn.getIndex(pXY2);
								inf_scenes[nscene][indx0] = true;
								inf_scenes[nscene][indx1] = true;
								inf_scenes[nscene][indx2] = true;
							}
						}
					}
				}
			}
		}
		if (debug) {
			double [][] dbg_inf_scenes = new double [num_scenes][inf_scenes[0].length];
			for (int nscene = 0; nscene < num_scenes; nscene++) {
				for (int nTile = 0; nTile < dbg_inf_scenes[nscene].length; nTile++){
					dbg_inf_scenes[nscene][nTile] = inf_scenes[nscene][nTile]? 1.0 : 0.0;
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_inf_scenes,
					tilesX,
					tilesY,
					true,
					"inf_scenes");
		}
		return inf_scenes;
	}
	
	/**
	 * Filter tiles and generate target disparities array with NaN for the removed tiles.
	 * Originally planned to generate both infinity tiles with target disparity equal to
	 * the common expected infinity value and non-infinity ones with their own target
	 * disparities, but later decided to do that in 2 passes - first processing only
	 * infinity tiles (much smaller number), second pass - treat all as non-infinity
	 * (own disparity) to get more data for LY offsets.   
	 * @param clust_size cluster size
	 * @param inf_range full range for the acceptable infinity tiles centered at inf_disp_ref
	 * @param scene_range disparity range for non-infinity in the same cluster in the same scene
	 * @param min_num_inf minimal number of tiles (in all scenes total) in an infinity cluster
	 * @param scenes sequence of scenes ordered by increasing timestamps
	 * @param valid_tile [scene_index][tile_index] array indicating tiles with exactly one
	 *                   correlation maximum and valid LMA
	 * @param inf_disp_ref average disparity at infinity for the reference (last) scene
	 * @param is_infinity which tile of which scene is infinity. May be null (all non-infinity),
	 *                    if not - may be infinity from the composite depth map
	 * @param only_infinity - only generate tiles for infinity, skip non-infinity (set to NaN)
	 * @param in_num_tiles int array of total number of clusters or null. If not null, will return
	 *                     number of used tiles in each cluster
	 * @param in_inf_cluster boolean array of total number of clusters or null. If not null will
	 *                       return if the cluster is infinity (all non-infinity tiles will be
	 *                       removed from it).
	 * @param threadsMax maximal number of threads to use
	 * @param debug show debug images
	 * @return target disparity per scene per tile. Uses inf_disp_ref for infinity tiles and disparity
	 *                maps for non-ifinity tiles (restored from the model directory to the scenes array).
	 *                All unused (filtered out) tiles have NaN for disparity.
	 */
	public static double [][] useTilesLY(
			final int          clust_size,
			final double       inf_range,     // full range centered at inf_disp_ref to be used as infinity
			final double       scene_range,   // disparity range for non-infinity in the same cluster
			final int          min_num_inf,   // Minimal number of tiles (in all scenes total) in an infinity cluster
			final QuadCLT []   scenes,        // ordered by increasing timestamps
			final boolean [][] valid_tile,    // tile with lma and single correlation maximum
			final double       inf_disp_ref,  // average disparity at infinity for ref scene
			final boolean [][] is_infinity,   // may be null (all non-infinity), if not - may be infinity from the composite depth map
			final boolean      only_infinity, // do not process non-infinity
			int []             in_num_tiles,     // null or array of number of clusters
			boolean []         in_inf_cluster,   // null or array of number of clusters, will return if cluster uses only infinite tiles
			final int          threadsMax,
			final boolean      debug){
		int last_scene_index = scenes.length-1;
		QuadCLT last_scene =   scenes[last_scene_index];
		final int num_scenes = scenes.length;
		final double inf_hrange = inf_range/2;
		final int tilesX = last_scene.tp.getTilesX();
		final int tilesY = last_scene.tp.getTilesY();
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clust_size);
		final int clusters = clustersX * clustersY;
		final double [][] target_disparities = new double [num_scenes][];
		final double [][] rslt_disparities =   new double [num_scenes][tilesX*tilesY];
		
		final int []     num_tiles =   (in_num_tiles != null)?   in_num_tiles :   (new int [clusters]);
		final boolean [] inf_cluster = (in_inf_cluster != null)? in_inf_cluster : (new boolean [clusters]);
		Arrays.fill(num_tiles, 0);
		Arrays.fill(inf_cluster, false);
		
		final int [] num_inf_cluster = new int[clusters];
		for (int nscene = 0; nscene < num_scenes; nscene++) {
			target_disparities[nscene] = scenes[nscene].getDSRBG()[0];
			Arrays.fill(rslt_disparities[nscene], Double.NaN);
			
		}
 		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		if (is_infinity != null) {
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) {
							int clustX = nClust % clustersX;
							int clustY = nClust / clustersX;
							for (int ctY = 0; ctY < clust_size; ctY++) {
								int tileY = clustY * clust_size + ctY;
								if (tileY < tilesY) {
									for (int ctX = 0; ctX < clust_size; ctX++) {
										int tileX = clustX * clust_size + ctX;
										if (tileX < tilesX) {
											int nTile = tileY*tilesX+tileX;
											for (int nscene = 0; nscene < num_scenes; nscene++) {
												if (is_infinity[nscene][nTile] &&
														valid_tile[nscene][nTile] && // next may be NaN
														(Math.abs(target_disparities[nscene][nTile] - inf_disp_ref) <= inf_hrange)) {
													num_inf_cluster[nClust]++;
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
			// mark clusters as infinity if they have enough infinity tiles
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) {
							if (num_inf_cluster[nClust] >= min_num_inf) {
								inf_cluster[nClust] = true;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			
			// Mark suitable infinity tiles (they are NaN now)
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) if (inf_cluster[nClust]) {
							int clustX = nClust % clustersX;
							int clustY = nClust / clustersX;
							num_tiles[nClust] = 0;
							for (int ctY = 0; ctY < clust_size; ctY++) {
								int tileY = clustY * clust_size + ctY;
								if (tileY < tilesY) {
									for (int ctX = 0; ctX < clust_size; ctX++) {
										int tileX = clustX * clust_size + ctX;
										if (tileX < tilesX) {
											int nTile = tileY*tilesX+tileX;
											for (int nscene = 0; nscene < num_scenes; nscene++) {
												if (is_infinity[nscene][nTile] &&
														valid_tile[nscene][nTile] && // next may be NaN
														(Math.abs(target_disparities[nscene][nTile] - inf_disp_ref) <= inf_hrange)) {
													rslt_disparities[nscene][nTile] = inf_disp_ref;
													num_tiles[nClust]++;
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
		
		if (!only_infinity) {
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] clust_disp = new double[clusters];
						// only use clusters that are not infinity
						for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) if (!inf_cluster[nClust]){
							int clustX = nClust % clustersX;
							int clustY = nClust / clustersX;
							//						int num_tot = 0;
							for (int nscene = 0; nscene < num_scenes; nscene++) {
								Arrays.fill(clust_disp,Double.NaN);
								double sum_disp = 0;
								int num_disp = 0;

								for (int ctY = 0; ctY < clust_size; ctY++) {
									int tileY = clustY * clust_size + ctY;
									if (tileY < tilesY) {
										for (int ctX = 0; ctX < clust_size; ctX++) {
											int tileX = clustX * clust_size + ctX;
											if (tileX < tilesX) {
												int ct = ctY * clust_size + ctX;
												int nTile = tileY*tilesX+tileX;
												if (valid_tile[nscene][nTile]){
													clust_disp[ct] = target_disparities[nscene][nTile];
													sum_disp += clust_disp[ct];
													num_disp ++;
												}
											}
										}
									}
								}
								while (num_disp > 0) {
									int imin = 0;
									int imax = 0;
									for (int ct = 0; ct < clusters; ct++) if (!Double.isNaN(clust_disp[ct])){
										if (!(clust_disp[ct] >= clust_disp[imin] )) imin = ct; 
										if (!(clust_disp[ct] <= clust_disp[imax] )) imax = ct; 
									}
									if ((clust_disp[imax] - clust_disp[imin]) <= scene_range) {
										break;
									}
									double disp_avg = sum_disp/num_disp;
									if ((clust_disp[imax] - disp_avg) > (disp_avg - clust_disp[imin])){
										sum_disp -= clust_disp[imax]; 
										clust_disp[imax] = Double.NaN;

									} else {
										sum_disp -= clust_disp[imin]; 
										clust_disp[imin] = Double.NaN;
									}
									num_disp --;
								}
								for (int ctY = 0; ctY < clust_size; ctY++) {
									int tileY = clustY * clust_size + ctY;
									if (tileY < tilesY) {
										for (int ctX = 0; ctX < clust_size; ctX++) {
											int tileX = clustX * clust_size + ctX;
											if (tileX < tilesX) {
												int ct = ctY * clust_size + ctX;
												int nTile = tileY*tilesX+tileX;
												if (!Double.isNaN(clust_disp[ct])) {
													rslt_disparities[nscene][nTile] = clust_disp[ct];
												}
											}
										}
									}
								}
								num_tiles[nClust] += num_disp;
								//								num_tot += num_disp;
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		}
		return rslt_disparities;
	}

	/**
	 * Calculate Lazy Eye data in 2 passes - one for infinity adjustment, another - for LY.
	 * @param clt_parameters parameters of calculation.
	 * @param scenes sequence of scenes ordered by increasing timestamps.
	 * @param valid_tile [scene_index][tile_index] array indicating tiles with exactly one
	 *                   correlation maximum and valid LMA.
	 * @param inf_disp_ref average disparity at infinity for the reference (last) scene
	 * @param is_scene_infinity which tile of which scene is infinity. 
	 * @param target_disparities if not null, should be initialized to double[2][][],
	 *        will return target disparities array used for infinity [0][][] and non-infinity [1][][].
	 * @param dbg_disparity_offset - add to all disparities (for testing).
	 * @param in_num_tiles null or int[2][] - will return number of tiles per cluster for
	 *                     infinity [0][] and non-infinity - [1][].
	 * @param corr_vector_delta null or extrinsic vector offset, applied to all scenes.
	 *        Only ERS parameters will be different, all the rest are the same for all scenes.
	 * @param nrefine number of disparity refines for non-infinity.
	 * @param threadsMax maximal number of threads to use.
	 * @param debug_level debug level, will show debug images if >-2.
	 * @return [2][clusters][LY_sclices] LY data for infinity [0][][] and non-infinity[1][][].
	 */
	public static double [][][] getLYDataInfNoinf(
			final CLTParameters  clt_parameters,
			final QuadCLT []     scenes,            // ordered by increasing timestamps
			final boolean [][]   valid_tile,        // tile with lma and single correlation maximum
			final double         inf_disp_ref,      // average disparity at infinity for ref scene // is_scene_infinity
			final boolean [][]   is_scene_infinity, // may be null, if not - may be infinity from the composite depth map
			final double  [][][] target_disparities,
			final double         dbg_disparity_offset,
			final int [][]       in_num_tiles,      // null or number of tiles per cluster to multiply strength
			final CorrVector     corr_vector_delta, // null or extrinsic vector offset, applied to all scenes
			final int            nrefine,           // number of disparity refines for non-inf
			final int            threadsMax,
			final int            debug_level){
		final double [][][] inf_noinf_lazy_eye_data = new double [2][][];
		int last_scene_index = scenes.length-1;
		QuadCLT last_scene = scenes[last_scene_index];
		int tilesX = last_scene.tp.getTilesX();
		int tilesY = last_scene.tp.getTilesY();
		int clustersX = (int) Math.ceil(1.0 * tilesX / clt_parameters.lyms_clust_size);
		int clustersY = (int) Math.ceil(1.0 * tilesY / clt_parameters.lyms_clust_size);
		int clusters = clustersX * clustersY;
		final int [][]     num_tiles =   (in_num_tiles != null)?   in_num_tiles :   (new int [2][]);
		for (int i = 0; i < num_tiles.length; i++) {
			num_tiles[i] = new int [clusters];
		}
		boolean [] inf_cluster  = new boolean [clusters]; // null;		
		boolean       debug      = debug_level > -2;
		// get for infinity only
		double [][] target_disparities_inf = MultisceneLY.useTilesLY(
				clt_parameters.lyms_clust_size,    // final int          clust_size,
				clt_parameters.lyms_inf_range,     // final double       inf_range,     // full range centered at inf_disp_ref to be used as infinity
				clt_parameters.lyms_scene_range,   // final double       scene_range,   // disparity range for non-infinity in the same cluster
				clt_parameters.lyms_min_num_inf,   // final int          min_num_inf,   // Minimal number of tiles (in all scenes total) in an infinity cluster
				scenes,        // final QuadCLT []   scenes,        // ordered by increasing timestamps
				valid_tile,    // final boolean [][] valid_tile,    // tile with lma and single correlation maximum
				inf_disp_ref,  // final double       inf_disp_ref,  // average disparity at infinity for ref scene // is_scene_infinity
				is_scene_infinity, // final boolean [][] is_infinity,   // may be null, if not - may be infinity from the composite depth map
				true,          // final boolean      only_infinity, // do not process non-infinity
				num_tiles[0],  // int []             in_num_tiles,     // null or array of number of clusters
				inf_cluster,   // boolean []         in_inf_cluster,   // null or array of number of clusters, will return if cluster uses only infinite tiles
				threadsMax,    // final int          threadsMax,
				debug);        // final boolean      debug)
		
		if (dbg_disparity_offset != 0.00) {
			for (int nscene = 0; nscene < target_disparities_inf.length; nscene++) {
				for (int i = 0; i < target_disparities_inf[nscene].length; i++) {
					if (!Double.isNaN(target_disparities_inf[nscene][i])) {
						target_disparities_inf[nscene][i] += dbg_disparity_offset;
					}
				}
			}
		}
		
		inf_noinf_lazy_eye_data[0] = 	MultisceneLY.getLYData( // TODO: show lazy_eye_data[][]
				clt_parameters,         // final CLTParameters  clt_parameters,
				clt_parameters.lyms_clust_size,             // final int            clust_size,
				scenes,                 // final QuadCLT []     scenes,        // ordered by increasing timestamps
				target_disparities_inf, // final double[][]     target_disparities,
				num_tiles[0],           // final int []         num_tiles,
				corr_vector_delta,      // final CorrVector     corr_vector_delta, // null or extrinsic vector offset, applied to all scenes
				threadsMax,             // final int            threadsMax,
				debug_level,            // final int            debug_level);
				((debug_level >-2)?"INF":null)); // String debug_suffix);
		// now for non-infinity:
		double [][] target_disparities_noinf = MultisceneLY.useTilesLY(
				clt_parameters.lyms_clust_size,    // final int          clust_size,
				clt_parameters.lyms_inf_range,     // final double       inf_range,     // full range centered at inf_disp_ref to be used as infinity
				clt_parameters.lyms_scene_range,   // final double       scene_range,   // disparity range for non-infinity in the same cluster
				clt_parameters.lyms_min_num_inf,   // final int          min_num_inf,   // Minimal number of tiles (in all scenes total) in an infinity cluster
				scenes,        // final QuadCLT []   scenes,        // ordered by increasing timestamps
				valid_tile,    // final boolean [][] valid_tile,    // tile with lma and single correlation maximum
				inf_disp_ref,  // final double       inf_disp_ref,  // average disparity at infinity for ref scene // is_scene_infinity
				null,          // final boolean [][] is_infinity,   // may be null, if not - may be infinity from the composite depth map
				false,         // final boolean      only_infinity, // do not process non-infinity
				num_tiles[1],  // int []             in_num_tiles,     // null or array of number of clusters
				inf_cluster,   // boolean []         in_inf_cluster,   // null or array of number of clusters, will return if cluster uses only infinite tiles
				threadsMax,    // final int          threadsMax,
				debug);        // final boolean      debug)
		if (dbg_disparity_offset != 0.00) {
			for (int nscene = 0; nscene < target_disparities_noinf.length; nscene++) {
				for (int i = 0; i < target_disparities_noinf[nscene].length; i++) {
					if (!Double.isNaN(target_disparities_noinf[nscene][i])) {
						target_disparities_noinf[nscene][i] += dbg_disparity_offset;
					}
				}
			}
		}
		inf_noinf_lazy_eye_data[1] = 	MultisceneLY.getLYData( // TODO: show lazy_eye_data[][]
				clt_parameters,           // final CLTParameters  clt_parameters,
				clt_parameters.lyms_clust_size,               // final int            clust_size,
				scenes,                   // final QuadCLT []     scenes,        // ordered by increasing timestamps
				target_disparities_noinf, // final double[][]     target_disparities,
				num_tiles[1],             // final int []         num_tiles,
				corr_vector_delta,        // final CorrVector     corr_vector_delta, // null or extrinsic vector offset, applied to all scenes
				threadsMax,               // final int            threadsMax,
				debug_level,              // final int            debug_level);
				((debug_level >-2)?"NOINF":null)); // String debug_suffix);

		if (debug) {
			showLY(
					clt_parameters, // CLTParameters clt_parameters,
					last_scene.getGeometryCorrection(), // GeometryCorrection gc,
					clustersX, // int           clustersX,
					clustersY, // int           clustersY,
					inf_noinf_lazy_eye_data[1], // double [][]   ly_data,
					"LY_noinf-initial"); // String        title);					
		}
		
		for (int n = 0; n < nrefine; n++) {
			// modify target disparity with disparity difference 
			updateTargetDisparities(
					clt_parameters,             // final CLTParameters clt_parameters,
					last_scene.tp,              // final TileProcessor tp,
					target_disparities_noinf,   // double [][] target_disparities,
					inf_noinf_lazy_eye_data[1], // double [][] ly_data)
					threadsMax);                // final int         threadsMax);
			inf_noinf_lazy_eye_data[1] = 	MultisceneLY.getLYData( // TODO: show lazy_eye_data[][]
					clt_parameters,           // final CLTParameters  clt_parameters,
					clt_parameters.lyms_clust_size,               // final int            clust_size,
					scenes,                   // final QuadCLT []     scenes,        // ordered by increasing timestamps
					target_disparities_noinf, // final double[][]     target_disparities,
					num_tiles[1],             // final int []         num_tiles,
					corr_vector_delta,        // final CorrVector     corr_vector_delta, // null or extrinsic vector offset, applied to all scenes
					threadsMax,               // final int            threadsMax,
					debug_level,              // final int            debug_level);
					((debug_level >-2)?("NOINF"+n):null)); // String debug_suffix);
			
			if (debug) {
				showLY(
						clt_parameters, // CLTParameters clt_parameters,
						last_scene.getGeometryCorrection(), // GeometryCorrection gc,
						clustersX, // int           clustersX,
						clustersY, // int           clustersY,
						inf_noinf_lazy_eye_data[1], // double [][]   ly_data,
						"LY_noinf-"+n); // String        title);					
			}
			
		}
				
		if (target_disparities != null) {
			target_disparities[0] = target_disparities_inf;
			target_disparities[1] = target_disparities_noinf;
		}
		return inf_noinf_lazy_eye_data;
	}	
	
	/**
	 * Display measure Lazy Eye data
	 * 
	 * @param clt_parameters configured parameters
	 * @param gc GeometryCorrection instance (e.g from one scene), needed to create
	 *           and initialize ExtrinsicAdjustment instance      
	 * @param clustersX number of clusters in a row
	 * @param clustersY number of clusters rows
	 * @param ly_data Lazy Eye data - array in cluster line-scan order, each may be null
	 *        or an array of LY data
	 * @param title Image title to use
	 */
	public static void showLY(
			CLTParameters clt_parameters,
			GeometryCorrection gc,
			int           clustersX,
			int           clustersY,
			double [][]   ly_data,
			String        title
			) {
		ExtrinsicAdjustment ea = new ExtrinsicAdjustment (
				gc,                            // GeometryCorrection gc,
				clt_parameters.lyms_clust_size, // int         clusterSize,
				clustersX,                      // int         clustersX,
				clustersY);                     // int         clustersY)
		int numSens = gc.getNumSensors();
		double [][] dbg_cluster = new double [ExtrinsicAdjustment.get_INDX_LENGTH(numSens)][clustersY * clustersX];
		for (int i = 0; i < dbg_cluster.length; i++) {
			Arrays.fill(dbg_cluster[i], Double.NaN);
		}
		for (int n = 0; n < ly_data.length; n++) {
			if (ly_data[n] != null) {
				for (int i = 0; i < ExtrinsicAdjustment.get_INDX_LENGTH(numSens); i++ ) {
					dbg_cluster[i][n] = ly_data[n][i];
				}
			}
		}
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_cluster,
				clustersX,
				clustersY,
				true,
				title,
				ea.data_titles); //  ExtrinsicAdjustment.DATA_TITLES);
	}
	
	/**
	 * Update per-tile target_disparity values for non-infinity tiles using LY data
	 * (LMA disparity difference) adding that values to all tiles of the corresponding 
	 * clusters in each scene where they are defined (non-NAN).
	 *  
	 * @param clt_parameters configuration parameters
	 * @param tp TileProcessor instance
	 * @param target_disparities per-scene, per-tile array of target disparities
	 *        (NaN for undefined tiles)
	 * @param ly_data Lazy Eye data (only differential disparity used)
	 * @param threadsMax maximal number of threads
	 */
	private static void updateTargetDisparities(
		final CLTParameters clt_parameters,
		final TileProcessor tp,
		final double [][]   target_disparities,
		final double [][]   ly_data,
		final int           threadsMax)
	{
		final int clust_size = clt_parameters.lyms_clust_size;
		final int tilesX = tp.getTilesX();
		final int tilesY = tp.getTilesY();
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clust_size);
		final int clusters = clustersX * clustersY;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
        for (int ithread = 0; ithread < threads.length; ithread++) {
            threads[ithread] = new Thread() {
                public void run() {
                    for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) if (ly_data[nClust] != null) {
                        int clustX = nClust % clustersX;
                        int clustY = nClust / clustersX;
                        double disp_diff = ly_data[nClust][ExtrinsicAdjustment.INDX_DIFF];
                        for (int ctY = 0; ctY < clust_size; ctY++) {
                            int tileY = clustY * clust_size + ctY;
                            if (tileY < tilesY) {
                                for (int ctX = 0; ctX < clust_size; ctX++) {
                                    int tileX = clustX * clust_size + ctX;
                                    if (tileX < tilesX) {
                                    	int nTile = tileY*tilesX+tileX;
                                        for (int nscene = 0; nscene < target_disparities.length; nscene++){
                                        	if ((target_disparities[nscene] != null) && !Double.isNaN(target_disparities[nscene][nTile])) {
                                        		target_disparities[nscene][nTile] += disp_diff;
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
	
	
	/**
	 * Calculate Lazy Eye data (now is called twice - separately for infinity and non-infinity)
	 * @param clt_parameters parameters of calculation.
	 * @param clust_size cluster size.
	 * @param scenes sequence of scenes ordered by increasing timestamps.
	 * @param target_disparities - per scene, per tile. NaN for unused tiles
	 * @param num_tiles null or number of tiles per cluster, used to multiply strength in the LY results
	 * @param corr_vector_delta null or extrinsic vector offset, applied to all scenes.
	 *        Only ERS parameters will be different, all the rest are the same for all scenes.
	 * @param threadsMax maximal number of threads
	 * @param debug_level debug level. Generates debug images if >-2.
	 * @param debug_suffix - add to image name, null - no images
	 * @return Lazy Eye data [clusters][LY_sclices]
	 */
	public static double [][] getLYData(
			final CLTParameters  clt_parameters,
			final int            clust_size,
			final QuadCLT []     scenes,        // ordered by increasing timestamps
			final double[][]     target_disparities,
			final int []         num_tiles, // null or number of tiles per cluster to multiply strength
			final CorrVector     corr_vector_delta, // null or extrinsic vector offset, applied to all scenes
			final int            threadsMax,
			final int            debug_level,
			final String         debug_suffix){
		int last_scene_index = scenes.length-1;
		QuadCLT last_scene = scenes[last_scene_index];
		int numSens = last_scene.getNumSensors();
		final int num_scenes = scenes.length;
		final int tilesX = last_scene.tp.getTilesX();
		final int tilesY = last_scene.tp.getTilesY();
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clust_size);
		final int clusters = clustersX * clustersY;
		final int num_pairs = Correlation2d.getNumPairs(last_scene.getNumSensors());
		
		boolean show_corr = (debug_suffix != null) ; //debug_level > -2;
		final float  [][][]     fclt_corr = show_corr ? (new float [tilesX * tilesY][][]) : null;

		final double [][][][][] dcorr_td_acc  = new double[num_pairs][][][][];
		final float  [][][][]   fcorr_td_acc  = new float [tilesY][tilesX][][];
		final float  [][][]     num_acc = new float [tilesY][tilesX][num_pairs];
	    
	    final double gpu_sigma_corr =     clt_parameters.getGpuCorrSigma(last_scene.isMonochrome());
	    final double gpu_sigma_rb_corr =  last_scene.isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr;
	    final double gpu_sigma_log_corr = clt_parameters.getGpuCorrLoGSigma(last_scene.isMonochrome());
	    double              disparity_corr = 0.0;
		int debugLevel = debug_level;
		final double [][] lazy_eye_data = new double [clustersY*clustersX][];
		ImageDtt image_dtt= new ImageDtt(
				numSens,
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				last_scene.isAux(),
				last_scene.isMonochrome(),
				last_scene.isLwir(),
				clt_parameters.getScaleStrength(last_scene.isAux()),
				last_scene.getGPU());
		if (last_scene.getGPU() != null) {
			last_scene.getGPU().setGpu_debug_level(debug_level);
		}
		final TpTask[][] tp_tasks_scenes = new TpTask[num_scenes][];
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt, last_scene.getNumSensors());

		for (int nscene = 0; nscene < num_scenes; nscene++) {
			QuadCLT scene = scenes[nscene];
			GeometryCorrection scene_gc = scene.getGeometryCorrection();
		    GeometryCorrection cond_gc = scene.hasGPU() ? null: scene_gc; // to skip calculating left for GPU
		    CorrVector corr_vector =  null; // backup copy of original, will be used to restore
		    CorrVector corr_vectorp = null; // copy to be modified and applied
		    if (corr_vector_delta != null) {
		    	corr_vector = scene_gc.getCorrVector().clone();
		    	corr_vectorp = corr_vector.clone();
		    	corr_vectorp.incrementVector(corr_vector_delta,  1.0); // 0.5 for p/m
		    	scene_gc.setCorrVector(corr_vectorp);
		    }
			scene.saveQuadClt(); // to re-load new set of Bayer images to the GPU (do nothing for CPU)
			scene. gpuResetCorrVector(); // is it needed - no, it is included in above
			if (debug_level > -2) {
				if (nscene == last_scene_index) {
					System.out.println("Correlating last scene");
				} else {
					System.out.println("Correlating scene "+nscene);
				}
			}
			double [][] disparity_array = new double[tilesY][tilesX];
			for (int tileY = 0; tileY < tilesY; tileY++) {
				System.arraycopy(target_disparities[nscene], tilesX*tileY, disparity_array[tileY], 0, tilesX);
			}
			tp_tasks_scenes[nscene] = GpuQuad.setTasks(
					numSens, // num_sensors,                  // final int                      num_cams,
					disparity_array,              // final double [][]	           disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
					disparity_corr,               // final double                   disparity_corr,
					null, // tile_op,             // final int [][]                 tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile
					cond_gc,                      // final GeometryCorrection       geometryCorrection,
					threadsMax);                  // final int                      threadsMax)       // maximal number of threads to launch

			if (scene.hasGPU()) {
				float  [][][][]     fcorr_td =       new float[tilesY][tilesX][][];
				image_dtt.quadCorrTD( // maybe remove "imageDtt."
						clt_parameters.img_dtt,            // final ImageDttParameters imgdtt_params,    // Now just extra correlation parameters, later will include, most others
						tp_tasks_scenes[nscene],           // *** will be updated inside from GPU-calculated geometry
						fcorr_td,                          // fcorrs_td[nscene],                 // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
//						scene_gc,     //
						clt_parameters.gpu_sigma_r,        // 0.9, 1.1
						clt_parameters.gpu_sigma_b,        // 0.9, 1.1
						clt_parameters.gpu_sigma_g,        // 0.6, 0.7
						clt_parameters.gpu_sigma_m,        //  =       0.4; // 0.7;
						gpu_sigma_rb_corr,                 // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
						gpu_sigma_corr,                    //  =    0.9;gpu_sigma_corr_m
						gpu_sigma_log_corr,                // final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
						clt_parameters.corr_red,           // +used
						clt_parameters.corr_blue,          // +used
						mcorr_sel,                         // final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
						threadsMax,       // maximal number of threads to launch
						debugLevel);
				if (image_dtt.getGPU().getGpu_debug_level() > -1) {
					System.out.println("==ooo=after image_dtt.quadCorrTD()");
				}
// Verify tasks are now updated
				OpticalFlow.accumulateCorrelations(
						num_acc,       // final int [][][]     num_acc,          // number of accumulated tiles [tilesY][tilesX][pair]
						fcorr_td,      // final float [][][][] fcorr_td,         // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs 
						fcorr_td_acc, // final float [][][][] fcorr_td_acc      // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
						threadsMax);
				if (image_dtt.getGPU().getGpu_debug_level() > -1) {
					System.out.println("==ooo=accumulateCorrelations()");
				}
			} else {
				image_dtt.getCorrelation2d();
				double [][][][]     dcorr_td = new double[tp_tasks_scenes[nscene].length][][][];
				image_dtt.quadCorrTD( // clt_data [task][sensor][color][][];
						scene.getImageData(),                           // final double [][][]       image_data,      // first index - number of image in a quad
						scene_gc.getSensorWH()[0], // final int                 width,
						tp_tasks_scenes[nscene],                        // tp_tasks,                            // final TpTask []           tp_tasks,
						clt_parameters.img_dtt,              // final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
						// dcorr_td should be either null, or double [tp_tasks.length][][];
						dcorr_td,                            // final double [][][][]     dcorr_td,        // [tile][pair][4][64] sparse by pair transform domain representation of corr pairs
						// no combo here - rotate, combine in pixel domain after interframe
						scene.getCltKernels(), //   clt_kernels,                         // final double [][][][][][] clt_kernels,     // [sensor][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
//						clt_parameters.kernel_step,          // final int                 kernel_step,
						clt_parameters.clt_window,           // final int                 window_type,
						clt_parameters.corr_red,             // final double              corr_red,
						clt_parameters.corr_blue,            // final double              corr_blue,
						mcorr_sel,                           // final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert
						clt_parameters.tileX,                // final int                 debug_tileX,
						clt_parameters.tileY,                // final int                 debug_tileY,
						threadsMax,                          // final int                 threadsMax,       // maximal number of threads to launch
						debugLevel);                         // final int                 globalDebugLevel);
				OpticalFlow.accumulateCorrelations(
						tp_tasks_scenes[nscene],      // final TpTask []         tp_tasks,
						num_acc,                      // final int [][][]        num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
						dcorr_td,                     // final double [][][][][] dcorr_td,    // [tile][pair][4][64] sparse transform domain representation of corr pairs 
						dcorr_td_acc,                 // final double [][][][][] dcorr_td_acc // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
						threadsMax);
			}
			// restore scene's geometry correction corr_vector
			if (corr_vector != null) {
				scene_gc.setCorrVector(corr_vector) ; // restore
				scene.gpuResetCorrVector();
			}
		} // for (int nscene = 0; nscene < num_scenes; nscene++)
		// Normalize accumulated correlations
		double [][]    tile_clust_weights = new double[tilesY][tilesX];
		if (last_scene.hasGPU()) {
			OpticalFlow.accumulateCorrelationsAcOnly( // normalize
					clust_size,         // final int            clust_size,
					num_acc,            // final float [][][]   num_acc,        // number of accumulated tiles [tilesY][tilesX][pair]
					fcorr_td_acc,       // final float [][][][] fcorr_td_acc,   // [tilesY][tilesX][pair][256] sparse transform domain representation of corr pairs
					tile_clust_weights, // final double [][]    tile_clust_weights,  // null or  [tilesY][tilesX]
					threadsMax);        // final int            threadsMax
			final double [][] combo_pXpYD= getAveragePxPyD(
					clust_size,                  // final int            clust_size, 
					tp_tasks_scenes, // [last_scene_index], // final TpTask []      tp_tasks,
					tile_clust_weights,          // final double [][]    tile_clust_weights,
					threadsMax);                 //final int            threadsMax)
			TpTask[] tp_tasks_combo =  GpuQuad.setInterTasks(
					last_scene.getNumSensors(),
					last_scene.getErsCorrection().getSensorWH()[0],
					!last_scene.hasGPU(),          // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
					combo_pXpYD,                   // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
					null,                          // final boolean []          selection, // may be null, if not null do not  process unselected tiles
					last_scene.getErsCorrection(), // final GeometryCorrection  geometryCorrection,
					disparity_corr,                // final double              disparity_corr,
					0, // margin,                  // final int                 margin,      // do not use tiles if their centers are closer to the edges
					null,                          // final boolean []          valid_tiles,            
					threadsMax);                   // final int                 threadsMax)  // maximal number of threads to launch
			// Need to initialize and update tp_tasks_combo from GPU
			image_dtt.setUpdateTasksGPU(
					tp_tasks_combo); // final TpTask[]            tp_tasks
			
			double [][][]     dcorr_tiles = show_corr? (new double [tp_tasks_combo.length][][]):null;
			final double[][] disparity_map = new double [image_dtt.getDisparityTitles().length][];
			final double [][][][]     ddnd = new double [tilesY][tilesX][][];
			image_dtt.clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					clt_parameters.img_dtt,		   // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					fcorr_td_acc,		 	     // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
					num_acc,                       // float [][][]                num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null)       
					null, // dcorr_weight,                  // double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
					clt_parameters.gpu_corr_scale, //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					clt_parameters.getGpuFatZero(last_scene.isMonochrome()),   // final double     gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
					image_dtt.transform_size - 1,  // final int                 gpu_corr_rad,    // = transform_size - 1 ?
			        // The tp_tasks data should be decoded from GPU to get coordinates
					tp_tasks_combo,                  // final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					last_scene.getErsCorrection().getRXY(false), // final double [][]         rXY,             // from geometryCorrection
					// next both can be nulls
					null,                          // final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
				    // combo will be added as extra pair if mcorr_comb_width > 0 and clt_corr_out has a slot for it
					// to be converted to float
					dcorr_tiles,                   // final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// When clt_mismatch is non-zero, no far objects extraction will be attempted
					//optional, may be null
					disparity_map,                 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					ddnd,                          // [tilesY][tilesX][num_sensors][2] data for LY. Should be either null or [tilesY][tilesX][][]. disparity_map should be non-null
					clt_parameters.correlate_lma,  // final boolean             run_lma,         // calculate LMA, false - CM only
		  		    // define combining of all 2D correlation pairs for CM (LMA does not use them)
					clt_parameters.img_dtt.mcorr_comb_width, //final int                 mcorr_comb_width,  // combined correlation tile width (set <=0 to skip combined correlations)
					clt_parameters.img_dtt.mcorr_comb_height,//final int                 mcorr_comb_height, // combined correlation tile full height
					clt_parameters.img_dtt.mcorr_comb_offset,//final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					clt_parameters.img_dtt.mcorr_comb_disp,	 //final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
					clt_parameters.clt_window,     // final int                 window_type,     // GPU: will not be used
					clt_parameters.tileX,          // final int                 debug_tileX,
					clt_parameters.tileY,          // final int                 debug_tileY,
					threadsMax,                    // final int                 threadsMax,      // maximal number of threads to launch
					null,                          // final String              debug_suffix,
					debug_level + 2); // -1 );              // final int                 globalDebugLevel)
			ImageDtt.convertFcltCorr(
					dcorr_tiles, // double [][][] dcorr_tiles,// [tile][sparse, correlation pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					fclt_corr);  // float  [][][] fclt_corr) //  new float [tilesX * tilesY][][] or null

			if (show_corr){ // -1
				float [][] accum_2d_img = ImageDtt.corr_partial_dbg( // not used in lwir
						fclt_corr,      // final float  [][][]     fcorr_data,       // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
						tp_tasks_combo, // final TpTask []         tp_tasks,        //
						tilesX,         //final int               tilesX,
						tilesY,         //final int               tilesX,
						2*image_dtt.transform_size - 1,	// final int               corr_size,
						1000,           // will be limited by available layersfinal int               layers0,
						clt_parameters.corr_border_contrast, // final double            border_contrast,
						threadsMax,     // final int               threadsMax,     // maximal number of threads to launch
						debug_level);   // final int               globalDebugLevel)
				int [] wh = new int[2];
				float [][] accum_2d_decimated = ImageDtt.corr2d_decimate( // not used in lwir
						accum_2d_img,  // final float [][]        corr2d_img,
						tilesX,        // final int               tilesX,
						tilesY,        // final int               tilesY,
						clust_size,    // final int               clust_size,
						wh,            // final int []            wh, 
						threadsMax,    // final int               threadsMax,     // maximal number of threads to launch
						debug_level);  // final int               globalDebugLevel)
				if (accum_2d_img != null) {
					String [] titles = new String [accum_2d_img.length]; // dcorr_tiles[0].length];
					int ind_length = image_dtt.getCorrelation2d().getCorrTitles().length;

					System.arraycopy(image_dtt.getCorrelation2d().getCorrTitles(), 0, titles, 0, ind_length);
					for (int i = ind_length; i < titles.length; i++) {
						titles[i] = "combo-"+(i - ind_length);
					}
					(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
							accum_2d_decimated,
							wh[0],
							wh[1],
							true,
							last_scene.getImageName()+"-CORR-DECIMATED"+clust_size+"-"+debug_suffix,
							titles);
				}
				double [][] disparity_map_decimated = ImageDtt.corr2d_decimate( // not used in lwir
						disparity_map,  // final float [][]        corr2d_img,
						tilesX,        // final int               tilesX,
						tilesY,        // final int               tilesY,
						clust_size,    // final int               clust_size,
						wh,            // final int []            wh, 
						threadsMax,    // final int               threadsMax,     // maximal number of threads to launch
						debug_level);  // final int               globalDebugLevel)
				
				(new ShowDoubleFloatArrays()).showArrays(
						disparity_map_decimated,
						wh[0],
						wh[1],
						true,
						"disparity_map_decimated"+"-"+debug_suffix,
						ImageDtt.getDisparityTitles(last_scene.getNumSensors(),last_scene.isMonochrome()) // ImageDtt.DISPARITY_TITLES
						);
				
			}
			
			final double [][] rXY = last_scene.getErsCorrection().getRXY(false);
			final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) {
							int clustX = nClust % clustersX;
							int clustY = nClust / clustersX;
							int tileX0 = clustX * clust_size;
							int tileY0 = clustY * clust_size;
							int nTile0 = tileY0 * tilesX + tileX0;
							//Maybe single test is enough
							if ((ddnd[tileY0][tileX0] != null) && !Double.isNaN(disparity_map[ImageDtt.DISPARITY_INDEX_POLY][nTile0])) {
								lazy_eye_data[nClust] = new double [ExtrinsicAdjustment.get_INDX_LENGTH(numSens)];
								// Number of tiles, strength or a product
								double strength = disparity_map[ImageDtt.DISPARITY_INDEX_POLY + 1][nTile0];
								if (num_tiles != null) {
									strength *= num_tiles[nClust];
								}
								lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_STRENGTH] = strength;
								lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_TARGET] = combo_pXpYD[nTile0][2];
								lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_DIFF] =   disparity_map[ImageDtt.DISPARITY_INDEX_POLY][nTile0];
								lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_DISP] = lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_TARGET]+
										lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_DIFF];
								lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_PX + 0] = combo_pXpYD[nTile0][0];
								lazy_eye_data[nClust][ExtrinsicAdjustment.INDX_PX + 1] = combo_pXpYD[nTile0][1];
								for (int cam = 0; cam < ddnd[tileY0][tileX0].length; cam++) {
									lazy_eye_data[nClust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] =
											ddnd[tileY0][tileX0][cam][0] * rXY[cam][0] - ddnd[tileY0][tileX0][cam][1] * rXY[cam][1];
									lazy_eye_data[nClust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] =
											ddnd[tileY0][tileX0][cam][0] * rXY[cam][1] + ddnd[tileY0][tileX0][cam][1] * rXY[cam][0];
									lazy_eye_data[nClust][ExtrinsicAdjustment.get_INDX_DD0(numSens) + cam] = ddnd[tileY0][tileX0][cam][0];
									lazy_eye_data[nClust][ExtrinsicAdjustment.get_INDX_ND0(numSens) + cam] = ddnd[tileY0][tileX0][cam][1];
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			// fill what is needed from tasks
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTask = ai.getAndIncrement(); nTask < tp_tasks_combo.length; nTask = ai.getAndIncrement()) {
							int tileY =tp_tasks_combo[nTask].getTileY();
							int tileX =tp_tasks_combo[nTask].getTileX();
							int clustX = tileX / clust_size; // should be no remainder !
							int clustY = tileY / clust_size; // should be no remainder !
							if (((tileX % clust_size) != 0) || ((tileY % clust_size) != 0)) {
								System.out.println("fcltCorrTD(): tileX="+tileX+",tileY="+tileY+" are not in the top left corner of a cluster");
							}
							double [][] centersXY = tp_tasks_combo[nTask].getDoubleXY();
							double [][] disp_dist = tp_tasks_combo[nTask].getDoubleDispDist();
							int nClust = clustY * clustersX + clustX;
							if (lazy_eye_data[nClust] == null) {
								System.out.println("getLYData(): lazy_eye_data["+nClust+"] == null");
								continue;
							}
							for (int cam = 0; cam < numSens; cam++) {
								lazy_eye_data[nClust][ExtrinsicAdjustment.get_INDX_DYDDISP0(numSens) + cam] += disp_dist[cam][2]; // for ERS only
								lazy_eye_data[nClust][ExtrinsicAdjustment.get_INDX_PYDIST(numSens) + cam] += centersXY[cam][1]; // for ERS only
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
		} else {
			OpticalFlow.accumulateCorrelations(
					num_acc,       // final int [][][]        num_acc,     // number of accumulated tiles [tilesY][tilesX][pair]
					dcorr_td_acc, // final double [][][][][] dcorr_td_acc // [pair][tilesY][tilesX][4][64] sparse transform domain representation of corr pairs
					threadsMax);					
			throw new IllegalArgumentException ("CPU version not yet supported");
		}
		
		return lazy_eye_data;
	}
	
	/**
	 * Calculate average pX, pY, Disparity (pXpYD) for each cluster to be used in LY data.
	 * @param clust_size cluster size
	 * @param tp_tasks TpTask array to read center positions and disparity
	 * @param tile_clust_weights [tilesY][tilesX] - tile weights within each cluster
	 *                           (sum of the values for each cluster should be 1.0) 
	 * @param threadsMax maximal number of threads to use
	 * @return pXpYD array [tiles][] with only top-left corner tile of each cluster is
	 *         non-null and contain {pX, pY, Disparity}
	 */
	public static double [][] getAveragePxPyD(
			final int            clust_size, 
			final TpTask [][]    tp_tasks,
			final double [][]    tile_clust_weights,
			final int            threadsMax
			){
		final int tilesY =    tile_clust_weights.length;
		final int tilesX =    tile_clust_weights[0].length;
		final int tiles = tilesY * tilesX;
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clust_size);
		final int clusters =  clustersX * clustersY;
		final double [][] combo_pXpYD = new double[tiles][]; // will only have non-zero for top left corners of each cluster
		
		// average disparity, set centerXY. Using equal weights for each scene
		final double []   scene_tile_weight = new double [tiles];
		final double [][] sum_pXpYD = new double [tiles][];
		
		
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int iscene = 0; iscene < tp_tasks.length; iscene++) {
			final int nscene = iscene;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int nTask = ai.getAndIncrement(); nTask < tp_tasks[nscene].length; nTask = ai.getAndIncrement()) {
							int tileY =tp_tasks[nscene][nTask].getTileY();
							int tileX =tp_tasks[nscene][nTask].getTileX();
							int tile = tileY * tilesX + tileX;
							
							scene_tile_weight[tile] += 1.0;
							if (sum_pXpYD[tile] == null) {
								sum_pXpYD[tile] = new double[3];
							}
							sum_pXpYD[tile][0] += tp_tasks[nscene][nTask].getDoubleCenterXY()[0];
							sum_pXpYD[tile][1] += tp_tasks[nscene][nTask].getDoubleCenterXY()[1];
							sum_pXpYD[tile][2] += tp_tasks[nscene][nTask].getTargetDisparity();
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
		}
		//iterate by clusters, tiles in each cluster, ...
		for (int ithread = 00; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < clusters; nClust = ai.getAndIncrement()) {
						int clustX = nClust % clustersX;
						int clustY = nClust / clustersX;
						int tileX0 = clust_size * clustX;  
						int tileY0 = clust_size * clustY;
						int tile0 = tileY0 * tilesX + tileX0;
						for (int ctY = 0; ctY < clust_size; ctY++) {
							int tileY = clustY * clust_size + ctY;
							if (tileY < tilesY) {
								for (int ctX = 0; ctX < clust_size; ctX++) {
									int tileX = clustX * clust_size + ctX;
									if (tileX < tilesX) {
										if (tile_clust_weights[tileY][tileX] >  0) {
											if (combo_pXpYD[tile0] == null) {
												combo_pXpYD[tile0] = new double [3];
											}
											int tile = tileY * tilesX + tileX;
											if ( scene_tile_weight[tile] <= 0) {
												System.out.println("getAveragePxPyD(): unexpected zero weight for tile="+tile);
											}
											double w = tile_clust_weights[tileY][tileX]/ scene_tile_weight[tile]; // sum of tile_clust_weights[tileY][tileX] == 1.0
											for (int i = 0; i < 3; i++) {
												combo_pXpYD[tile0][i] += w * sum_pXpYD[tile][i];
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
		return combo_pXpYD;
	}
	
	/**
	 * Merge Lazy Eye data from two separate ones - for infinity (measuring disparity offset
	 * with a goal to have disparity be exactly zero for objects at infinity (such as clouds)
	 * and "lazy eye" measurements. It is needed as there are too few measurements for infinity
	 * objects, so lazy eye data is measured separately for the same tiles and combined with
	 * disparity measurements for infinity only   
	 * @param adjust_mode uses enum MSLY_MODE: INF_ONLY, NOINF_ONLY, and INF_NOINF (typical)
	 * @param lazy_eye_data2 [2][clusters][LY parameters] [0][][] - infinity measurements
	 *        (disparity at infinity is used) ,[1][][] - non-infinity measurements (LY data
	 *        is used). 
	 * @param force_disparity if non null, should be initialized as boolean[clusters], will be
	 *        updated to have true for tiles that have disparity data (offset for infinity)
	 * @return combined LY array, compatible with (older) single-scene LY adjustment.
	 */
	public static double [][] mergeLY(
			MSLY_MODE      adjust_mode,	
			double [][][]  lazy_eye_data2,
			boolean []     force_disparity // null or [clusters]
			)
			{
		double [][] lazy_eye_data = null;
		int clusters = 0;
		for (int i = 0; i < lazy_eye_data2.length; i++) {
			if (lazy_eye_data2[i] != null) {
				clusters = lazy_eye_data2[i].length;
				break;
			}
		}
		if (adjust_mode == MSLY_MODE.INF_NOINF) {
			if ((lazy_eye_data2.length < 2) || (lazy_eye_data2[1] == null)) {
				adjust_mode = MSLY_MODE.INF_ONLY;
			} else if (lazy_eye_data2[0] == null) {
				adjust_mode = MSLY_MODE.NOINF_ONLY;
			}
		}
		switch (adjust_mode) {
		case INF_ONLY: 
			if (lazy_eye_data2[0] == null) {
				return null;
			}
			lazy_eye_data = lazy_eye_data2[0];
			if (force_disparity != null) {
				for (int i = 0; i < clusters; i++) {
					force_disparity[i] = lazy_eye_data[i] != null;
				}
			}
		break;
		case NOINF_ONLY:
			if ((lazy_eye_data2.length > 1) && (lazy_eye_data2[1] == null)) {
				return null;
			}
			lazy_eye_data = lazy_eye_data2[1];
			if (force_disparity != null) {
				for (int i = 0; i < clusters; i++) {
					force_disparity[i] = false;
				}
			}
			
		break;
		case INF_NOINF: // both lazy_eye_data2[0] and lazy_eye_data2[1] are non-nulls
			if (force_disparity != null) {
				for (int i = 0; i < clusters; i++) {
					force_disparity[i] = (lazy_eye_data2[0][i] != null) && (lazy_eye_data2[1][i] != null); // do not process w/o noinf
				}
			}
			lazy_eye_data = new double [clusters][];
			for (int i = 0; i < clusters; i++) if (lazy_eye_data2[1][i] != null){
				lazy_eye_data[i] = lazy_eye_data2[1][i].clone();
				if (lazy_eye_data2[0][i] != null) {
					for (int n = 0;n < lazy_eye_data2[0][i].length; n++) {
						switch (n){
						case ExtrinsicAdjustment.INDX_DISP:
						case ExtrinsicAdjustment.INDX_STRENGTH:
						case ExtrinsicAdjustment.INDX_TARGET:
						case ExtrinsicAdjustment.INDX_DIFF:
							lazy_eye_data[i][n] =  lazy_eye_data2[0][i][n];
						}
					}
				}
			}			
			break;
		}
		return lazy_eye_data;
	}
	
	/**
	 * A placeholder to move from TwoQuadCLT
	 * @param clt_parameters
	 * @param adjust_mode
	 * @param scenes
	 * @param lazy_eye_data2
	 * @param valid_tile
	 * @param inf_disp_ref
	 * @param is_scene_infinity
	 * @param update_disparity
	 * @param threadsMax
	 * @param updateStatus
	 * @param debugLevel
	 * @return
	 */
	public static boolean processLYdata(
			final CLTParameters  clt_parameters,
			MSLY_MODE            adjust_mode,			
			final QuadCLT []     scenes,        // ordered by increasing timestamps
			final double [][][]  lazy_eye_data2,
			final boolean [][]   valid_tile,        // tile with lma and single correlation maximum
			final double         inf_disp_ref,      // average disparity at infinity for ref scene // is_scene_infinity
			final boolean [][]   is_scene_infinity, // may be null, if not - may be infinity from the composite depth map
			boolean              update_disparity, // re-measure disparity before measuring LY
			final int            threadsMax,  // maximal number of threads to launch
			final boolean        updateStatus,
			final int            debugLevel) {
		int last_scene_index = scenes.length-1;
		QuadCLT last_scene = scenes[last_scene_index];
		int numSens = last_scene.getNumSensors();
		final int num_scenes = scenes.length;
		final int tilesX = last_scene.tp.getTilesX();
		final int tilesY = last_scene.tp.getTilesY();
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clt_parameters.lyms_clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clt_parameters.lyms_clust_size);
		final int clusters = clustersX * clustersY;
		final int num_pairs = Correlation2d.getNumPairs(last_scene.getNumSensors());
	
		ExtrinsicAdjustment ea = new ExtrinsicAdjustment (
				last_scene.getErsCorrection(),  // GeometryCorrection gc,
				clt_parameters.lyms_clust_size, // int         clusterSize,
				clustersX,                      // int         clustersX,
				clustersY);                     // int         clustersY)
		double [][][] dbg_cluster = new double [lazy_eye_data2.length][ExtrinsicAdjustment.get_INDX_LENGTH(last_scene.getNumSensors())][clustersY * clustersX];
		for (int m = 0; m < lazy_eye_data2.length; m++) {
			for (int i = 0; i < dbg_cluster.length; i++) {
				Arrays.fill(dbg_cluster[m][i], Double.NaN);
			}
			for (int n = 0; n < lazy_eye_data2[m].length; n++) {
				if (lazy_eye_data2[m][n] != null) {
					for (int i = 0; i < ExtrinsicAdjustment.get_INDX_LENGTH(last_scene.getNumSensors()); i++ ) {
						dbg_cluster[m][i][n] = lazy_eye_data2[m][n][i];
					}
				}
			}
		}
		for (int m = 0; m < lazy_eye_data2.length; m++) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_cluster[m],
					clustersX,
					clustersY,
					true,
					last_scene.getImageName()+"-lazy_eye_data-"+m,
					ea.data_titles); //  ExtrinsicAdjustment.DATA_TITLES);
		}
		ea.showInput(lazy_eye_data2[0],"first_data-inf");
		ea.showInput(lazy_eye_data2[1],"first_data-noinf");
		return true;
	}
	
	// apply delta to each parameter, perform LY measurement and calculate difference
	
	/**
	 * Measure approximate derivatives of the LY data the subcameras poses
	 * @param clt_parameters configuration parameters
	 * @param scenes array of consecutive camera scenes, ordered by increasing timestamps
	 * @param lazy_eye_data2 a pair of {infinity LY, non-infinity LY}
	 * @param valid_tile per-scene, per tile array of "valid" tiles - ones that have exactly
	 *        one obtained with LMA correlation maximum
	 * @param inf_disp_ref avarage disparity of infinity tiles
	 * @param is_scene_infinity (may be null) - per scene, per tile - tile predicted to be
	 *        of infinity object, determined by applying rotation+movement of the composite
	 *        depth map 
	 * @param threadsMax maximal number of threads
	 * @param delta delta for measuring derivatives, scaled for different parameters internally.
	 *        Typical value 0.001 (0.01 is too high for some tiles)
	 * @param use_tarz if true - use sensor pose angles (tilt, azimuth, roll, zoom), false -
	 *        use symmetrical parameters (linear combination of T,A,R,Z)
	 * @param debugLevel debug level (>-2 - threshold)
	 */
	public static void debugLYDerivatives(
			final CLTParameters     clt_parameters,
			final QuadCLT []        scenes,            // ordered by increasing timestamps
			double [][][]           lazy_eye_data2, // inf, no_inf
			final boolean [][]      valid_tile,        // tile with lma and single correlation maximum
			final double            inf_disp_ref,      // average disparity at infinity for ref scene // is_scene_infinity
			final boolean [][]      is_scene_infinity, // may be null, if not - may be infinity from the composite depth map
			final int               threadsMax,  // maximal number of threads to launch
			double                  delta,
			boolean                 use_tarz,  // derivatives by tarz, not symmetrical vectors
			final int               debugLevel)
	{
		double scale_tl = 0.3;
		double scale_az = 0.3;
		double scale_rl0 = 2.0;
		double scale_rl = 0.3;
		double scale_zoom = 0.3;
		double [] scales_imu = {
				150, // 8.814346720176489E-1,  // 5,  // 13 10000x omega-tilt
				150, // 7.071297501674136E-1,  // 5,  // 14 10000x omega az
				150, // 1.306306793587865E-0,  // 4,  // 15 10000x omega roll
				300, // 2.8929916645453735E-0, // 4,  // 16 10000x vx
				300, // 2.943408022525927E-0,  // 4,  // 17 10000x vy
				500.0}; // 390.6185365641268};    //4};  // 18 100000x vz
		//		  delta = 0.001; // should be 0.001
		final boolean debug_img = false; // true; // false;
		final int last_scene_index = scenes.length-1;
		final QuadCLT last_scene = scenes[last_scene_index];
		final int numSens = last_scene.getNumSensors();
		final int tilesX = last_scene.tp.getTilesX();
		final int tilesY = last_scene.tp.getTilesY();
		final int clustersX = (int) Math.ceil(1.0 * tilesX / clt_parameters.lyms_clust_size);
		final int clustersY = (int) Math.ceil(1.0 * tilesY / clt_parameters.lyms_clust_size);
		final int clusters = clustersX * clustersY;
		ExtrinsicAdjustment ea = new ExtrinsicAdjustment (
				last_scene.getErsCorrection(),  // GeometryCorrection gc,
				clt_parameters.lyms_clust_size, // int         clusterSize,
				clustersX,                      // int         clustersX,
				clustersY);                     // int         clustersY)

		double [] parameter_scales = new double [CorrVector.getLength(numSens)];	
		for (int i = 0; i < numSens; i++) {
			parameter_scales    [CorrVector.getRollIndex(numSens)+   i] = ((i > 0) || use_tarz)? scale_rl : scale_rl0;
			if (i < numSens - 1) {
				parameter_scales[CorrVector.getTiltIndex(numSens)+   i]=scale_tl;
				parameter_scales[CorrVector.getAzimuthIndex(numSens)+i]=scale_az;
				parameter_scales[CorrVector.getZoomIndex(numSens)+   i]=scale_zoom;
			}
		}

		for (int i = 0; i < scales_imu.length; i++) {
			parameter_scales[CorrVector.getIMUIndex(numSens)+   i] = scales_imu[i];
		}
		
		int num_ly = ExtrinsicAdjustment.get_INDX_LENGTH(numSens); // scan.getLazyEyeData().length;
		int num_pars = CorrVector.getLength(numSens); // curr_corr_arr.length;
		double [][][][] ly_diff2 = new double [lazy_eye_data2.length][][][]; //[num_pars][num_ly][clusters];
		for (int nly = 0; nly < ly_diff2.length; nly++) if (lazy_eye_data2[nly] != null){
			ly_diff2[nly] = new double [num_pars][num_ly][clusters];  
			for (int np = 0; np < num_pars; np++) {
				for (int nl = 0; nl < num_ly; nl++) {
					Arrays.fill(ly_diff2[nly][np][nl], Double.NaN);  
				}
			}
		}
		// save initial ly data
		double [][][] ly_initial2 = new double [lazy_eye_data2.length][][]; // [clusters][];
		double [][][] ly2 = lazy_eye_data2; // scan.getLazyEyeData();
		for (int nly = 0; nly < ly2.length; nly++) if (ly2[nly] != null) {
			ly_initial2[nly] = new double[clusters][];
			for (int cluster = 0; cluster < clusters; cluster++) if (ly2[nly][cluster]!=null){
				ly_initial2[nly][cluster] = ly2[nly][cluster].clone();
			}
		}
		System.out.println(last_scene.getGeometryCorrection().getCorrVector().toString());
		String [] titles;
		if (use_tarz) {
			titles = CorrVector.getCorrNames(numSens);
		} else {
			titles = new String [num_pars]; //ea.getSymNames(); // why "S" here, while it is tarz???
			for (int i = 0; i < num_pars; i++) {
				titles[i] = "S"+i;
			}
		}
	
		if (debug_img) {
			for (int nly = 0; nly < ly_initial2.length; nly++) if (ly_initial2[nly] != null ) {
				ea.showInput(
						ly_initial2[nly], // double[][] data,
						"drv_reference-"+SINF_NOINF[nly]);// String title);
			}
			System.out.println("Initial:\n"+last_scene.getGeometryCorrection().getCorrVector().toString(true));  // true - short out
			double min_strength = 0.1; // 0.23		  
			int [] pfmt = {8,3};
			if (debugLevel > -3) {
				for (int nly = 0; nly < lazy_eye_data2.length; nly++) if (lazy_eye_data2[nly] != null) {
					System.out.println(ea.stringWeightedLY(
							lazy_eye_data2[nly],   // double [][] data,
							null,                  // double [][] ref_data,
							min_strength,          // double min_strength,
							pfmt,                  // int [] format,
							"_00"));               // String suffix))
				}
			}
		}
		int            nrefine = 4;
		for (int npar = 0; npar < num_pars; npar++) {
			// perform asymmetric delta
			double [] par_inc = new double [num_pars];
			par_inc[npar] = delta * parameter_scales[npar];
			CorrVector corr_delta;			  
			if (use_tarz) {
				corr_delta = new CorrVector (last_scene.getGeometryCorrection(),par_inc);
			} else {
				corr_delta = last_scene.getGeometryCorrection().getCorrVector(par_inc,null); // , par_mask); all parameters			  
			}

			double rdelta = 1.0/ par_inc[npar];
			System.out.println(npar+": "+ titles[npar]+", scale="+rdelta); // +"\n"+(geometryCorrection.getCorrVector().toString()));
			System.out.println("delta:\n"+corr_delta.toString(true));  // true - short out
			// Measure new LY
			ly2 = getLYDataInfNoinf(
					clt_parameters,    // final CLTParameters  clt_parameters,
					scenes,            // final QuadCLT []     scenes,            // ordered by increasing timestamps
					valid_tile,        // final boolean [][]   valid_tile,        // tile with lma and single correlation maximum
					inf_disp_ref,      // final double         inf_disp_ref,      // average disparity at infinity for ref scene // is_scene_infinity
					is_scene_infinity, // final boolean [][]   is_scene_infinity, // may be null, if not - may be infinity from the composite depth map
					null,              // final double  [][][] target_disparities,
					0.0,               // final double         dbg_disparity_offset,
					null,              // final int [][]       in_num_tiles,      // null or number of tiles per cluster to multiply strength
					corr_delta,        // final CorrVector     corr_vector_delta, // null or extrinsic vector offset, applied to all scenes
					nrefine,              // final int            nrefine,           // number of disparity refines for non-inf
					threadsMax,        // final int            threadsMax,
					debugLevel);       //  final int            debug_level);			

			if (debug_img) {
				for (int nly = 0; nly < ly2.length; nly++) if (ly2[nly] != null ) {
					ea.showInput(
							ly2[nly], // double[][] data,
							"drv_par-"+SINF_NOINF[nly]);// String title);
				}
			}

			for (int nly = 0; nly < ly_diff2.length; nly++) if (ly_diff2[nly] != null ) {
				for (int cluster = 0; cluster < clusters; cluster++) if ((ly_initial2[nly][cluster] != null) && (ly2[nly][cluster]!=null)){
					for (int nl = 0; nl < ly_initial2[nly][cluster].length; nl++) {
						ly_diff2[nly][npar][nl][cluster] =  rdelta * (ly2[nly][cluster][nl] - ly_initial2[nly][cluster][nl]); 
					}
				}
			}
			if (debugLevel > -3) {
				double [][][] ly_diff12 = new double [ly_initial2.length][][];
				for (int nly = 0; nly < ly_initial2.length; nly++) if (ly_initial2[nly] != null ) {
					ly_diff12[nly] = new double [ly_initial2[nly].length][];
					for (int cluster = 0; cluster < clusters; cluster++) if ((ly_initial2[nly][cluster] != null) && (ly2[nly][cluster]!=null)){
						ly_diff12[nly][cluster] = new double [ly_initial2[nly][cluster].length]; // nullp
						for (int nl = 0; nl < ly_initial2[nly][cluster].length; nl++) {
							ly_diff12[nly][cluster][nl] =  rdelta * (ly2[nly][cluster][nl] - ly_initial2[nly][cluster][nl]); 
						}
					}
					double min_strength = 0.10; // 0.23		  
					int [] pfmt = {8,3};
					System.out.println(ea.stringWeightedLY( // nullp
							ly_diff12[nly],         // double [][] data,
							ly_initial2[nly],            // double [][] ref_data,
							min_strength,          // double min_strength,
							pfmt,                  // int [] format,
							"_d"+titles[npar]+"-"+SINF_NOINF[nly]));    // String suffix))
					
					if (debug_img) {
						// check nd bug
						String [] sens_titles = new String [numSens+1];
						double [][] dbg_nd = new double [numSens+1][clusters];
						for (int n = 0; n < dbg_nd.length; n++) {
							Arrays.fill(dbg_nd[n], Double.NaN);
							sens_titles[n] = "sens-"+n;
						}
						sens_titles[numSens] = "sum";
						for (int cluster = 0; cluster < clusters; cluster++) if (ly_diff12[nly][cluster] != null){
							dbg_nd[dbg_nd.length - 1][cluster] = 0.0;
							for (int n = 0; n < numSens; n++) {
								dbg_nd[n][cluster] =                  ly_diff12[nly][cluster][n + ExtrinsicAdjustment.get_INDX_ND0(numSens)] / rdelta;
								dbg_nd[dbg_nd.length - 1][cluster] += ly_diff12[nly][cluster][n + ExtrinsicAdjustment.get_INDX_ND0(numSens)] / rdelta / numSens;
							}
						}
						(new ShowDoubleFloatArrays()).showArrays(
								dbg_nd,
								clustersX,
								clustersY,
								true,
								"ND_test_npar-_"+npar+"-"+SINF_NOINF[nly],
								sens_titles);
					}
				}
			}
		}
		int gap = 10;
		int rows = ea.getRowsCols()[0];
		int cols = ea.getRowsCols()[1];
		int width  = cols * (ea.clustersX + gap) - gap;
		int height = rows * (ea.clustersY + gap) - gap;
		double [][][]dbg_img2 = new double [ly_diff2.length][num_pars][width*height];
		for (int nly = 0; nly < ly_initial2.length; nly++) if (ly_initial2[nly] != null ) {
			for (int par = 0; par < num_pars; par++) {
				for (int mode = 0; mode < ExtrinsicAdjustment.get_POINTS_SAMPLE(numSens); mode++) {
					int x0 = (mode % cols) * (ea.clustersX + gap);
					int y0 = (mode / cols) * (ea.clustersY + gap);
					for (int cluster = 0; cluster < clusters;  cluster++) {
						int x = x0 + (cluster % ea.clustersX);
						int y = y0 + (cluster / ea.clustersX);
						int pix = x + y * width;
						int indx = (mode == 0) ? ExtrinsicAdjustment.INDX_DIFF : (ExtrinsicAdjustment.get_INDX_DD0(numSens) + mode - 1);
						if (mode == 0) {
							dbg_img2[nly][par][pix] = -ly_diff2[nly][par][indx][cluster];
						} else {
							dbg_img2[nly][par][pix] =  ly_diff2[nly][par][indx][cluster];
						}
					}
				}
			}
		}
		for (int nly = 0; nly < dbg_img2.length; nly++) if (dbg_img2[nly] != null ) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img2[nly],
					width,
					height,
					true,
					"dLY_dpar_"+delta+"DINV"+"-"+SINF_NOINF[nly],
					titles);
		}
		dbg_img2 = new double [ly_diff2.length][num_pars][width*height];
		for (int nly = 0; nly < ly_initial2.length; nly++) if (ly_initial2[nly] != null ) {
			for (int par = 0; par < num_pars; par++) {
				for (int mode = 0; mode < ExtrinsicAdjustment.get_POINTS_SAMPLE(numSens); mode++) {
					int x0 = (mode % cols) * (ea.clustersX + gap);
					int y0 = (mode / cols) * (ea.clustersY + gap);
					for (int cluster = 0; cluster < clusters;  cluster++) {
						int x = x0 + (cluster % ea.clustersX);
						int y = y0 + (cluster / ea.clustersX);
						int pix = x + y * width;
						int indx = (mode == 0) ? ExtrinsicAdjustment.INDX_DIFF : (ExtrinsicAdjustment.INDX_X0 + mode - 1);
						if (mode == 0) {
							dbg_img2[nly][par][pix] = -ly_diff2[nly][par][indx][cluster];
						} else {
							dbg_img2[nly][par][pix] = ly_diff2[nly][par][indx][cluster];
						}
					}
				}
			}
		}
		for (int nly = 0; nly < dbg_img2.length; nly++) if (dbg_img2[nly] != null ) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img2[nly],
					width,
					height,
					true,
					"dLY_dpar_"+delta+"DINV"+"-XY"+"-"+SINF_NOINF[nly],
					titles);
		}
		return;
	}


}
