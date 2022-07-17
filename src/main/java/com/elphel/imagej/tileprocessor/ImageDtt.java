package com.elphel.imagej.tileprocessor;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.gpu.GPUTileProcessor;
import com.elphel.imagej.gpu.GpuQuad;
import com.elphel.imagej.gpu.TpTask;

import ij.ImagePlus;

//import Jama.Matrix;

public class ImageDtt extends ImageDttCPU {
	public boolean debug_strengths = false; // true;
	private final GpuQuad gpuQuad;

	public ImageDtt(
			int numSensors,
			int transform_size,
			ImageDttParameters imgdtt_params,
			boolean aux,
			boolean mono,
			boolean lwir,
			double scale_strengths,
			GpuQuad gpuQuadIn){
		super ( numSensors,
				transform_size,
				imgdtt_params,
				aux,
				mono,
				lwir,
				scale_strengths);
		gpuQuad = gpuQuadIn;
	}

	public ImageDtt(
			int numSensors,
			int transform_size,
			ImageDttParameters imgdtt_params,
			boolean aux,
			boolean mono,
			boolean lwir,
			double scale_strengths){
		super ( numSensors,
				transform_size,
				imgdtt_params,
				aux,
				mono,
				lwir,
				scale_strengths);
		gpuQuad = null;
	}
	
	public GpuQuad getGPU() {
		return this.gpuQuad;
	}
	
	public boolean [] getCorrMask() {
		if (gpuQuad != null) {
			return getCorrMask();
		}
		return super.getCorrMask();
	}
	
//	public double [][][][][][]
	
	//FIXME: pair combining that will not work with non-quad !!!
	public void clt_aberrations_quad_corr_GPU(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][4*64] TD of combo corrs: qud, cross, hor,vert
			                                           // each of the top elements may be null to skip particular combo type
			final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
			final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
                                                       // [tilesY][tilesX] should be set by caller
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
			final double [][][][]     texture_tiles,   // compatible with the CPU ones      
			final float [][]          texture_img,     // null or [3][] (RGB) or [4][] RGBA
			final Rectangle           texture_woi_pix, // null or generated texture location/size in pixels
			final float [][][]        iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
			// new parameters, will replace some other?
			final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
			final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
			final double              gpu_sigma_corr,   //  =    0.9;gpu_sigma_corr_m
			final int                 gpu_corr_rad,     // = transform_size - 1 ?
			final double              corr_red, // +used
			final double              corr_blue,// +used
			final double              max_corr_radius, // 3.9;
			final int                 corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
			
			final double              min_shot,        // +10.0;  // Do not adjust for shot noise if lower than
			final double              scale_shot,      // +3.0;   // scale when dividing by sqrt ( <0 - disable correction)
			final double              diff_sigma,      // +5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
			final double              diff_threshold,  // +5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
			final boolean             diff_gauss,      // true;  // when averaging images, use Gaussian around average as weight (false - sharp all/nothing)
			final double              min_agree,       // +3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			final boolean             dust_remove,     // +Do not reduce average weight when only one image differs much from the average
			final GeometryCorrection  geometryCorrection, // for GPU TODO: combine geometry corrections if geometryCorrection_main is not null
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final int                 window_type,    // GPU: will not be used
			final double              disparity_corr, // disparity at infinity
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return;
		}
		final boolean [][] saturation_imp = gpuQuad.quadCLT.saturation_imp;               // boolean [][] saturation_imp, // (near) saturated pixels or null
//gpuQuad
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;

//		final double [][] debug_offsets = new double[imgdtt_params.lma_dbg_offset.length][2];
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}


		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int numcol = isMonochrome()?1:3;
		// FIXME maybe something else is needed.
		// When switching from larger images to smaller ones requested texture was smaller than
		// still increased GPU window size
		if (texture_tiles != null) { // maybe something else is needed
			// GPUTileProcessor.DTT_SIZE
			if ((texture_tiles.length != gpuQuad.getTilesY()) ||
					(texture_tiles[0].length != gpuQuad.getTilesX())) {
				gpuQuad.deallocate4Images();
			}
		}
		

		final int width =  gpuQuad.getImageWidth();
		final int height = gpuQuad.getImageHeight();
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size; // still old - before
		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*getNumSensors()][tilesX*tilesY]) : null;
		// keep for now for mono, find out  what do they mean for macro mode
		
		if (isMonochrome()) {
			col_weights[0] = 1.0; // 2] = 1.0;// green color/mono
		} else {
			if (macro_mode) { // all the same as they now mean different
				//compensating Bayer correction
				col_weights[0] = 0.25; //  1.0/3;
				col_weights[1] = 0.25; //  1.0/3;
				col_weights[2] = 0.5; // 1.0/3;
			} else {
				col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
				col_weights[0] = corr_red *  col_weights[2];
				col_weights[1] = corr_blue * col_weights[2];
			}
		}
		final int corr_size = transform_size * 2 - 1;
		
		if ((globalDebugLevel > -10) && (disparity_corr != 0.0)){
			System.out.println(String.format("Using manual infinity disparity correction of %8.5f pixels",disparity_corr));
		}

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 1){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}


		if (globalDebugLevel > 1) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		
	
		// add optional initialization of debug layers here
		boolean need_macro = false;
		boolean need_corr = (clt_mismatch != null) || (fcorr_combo_td !=null) || (fcorr_td !=null) ; // (not the only reason)
		// skipping DISPARITY_VARIATIONS_INDEX - it was not used
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++) {
				if (isSliceBit(i) && ((disparity_modes & (1 << i)) != 0)) { 
					if ((i == OVEREXPOSED) && (saturation_imp == null)) {
						continue;
					}
					disparity_map[i] = new double [tilesY*tilesX];
					if (isCorrBit (i)) {
						need_corr = true;
					}
				} else if (isDiffIndex(i) && needImgDiffs(disparity_modes)){
					disparity_map[i] = new double [tilesY*tilesX];
					need_macro = true;
				} else if (isToneRGBIndex(i) && needTonesRGB(disparity_modes)){
					disparity_map[i] = new double [tilesY*tilesX];
					need_macro = true;
				}
			}
		}

		
		
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}

		/*
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256] - never used
		final double [] lt_window2 = new double [lt_window.length]; // squared - never used

		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];
		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
*/
		if (globalDebugLevel > 0) {
			System.out.println("macro_mode="+macro_mode);
		}

		final boolean use_main = geometryCorrection_main != null;
		boolean [] used_corrs = new boolean[1];
	    final int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		final TpTask[] tp_tasks =  gpuQuad.setTpTask( // tile-op is 80x64
				disparity_array, // final double [][]  disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
				disparity_corr,  // final double       disparity_corr,
				used_corrs,      // final boolean []   need_corrs,       // should be initialized to boolean[1] or null
				tile_op,         // final int [][]     tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				all_pairs,       // final int                      corr_mask,        // <0 - use corr mask from the tile tile_op, >=0 - overwrite all with non-zero corr_mask_tp
				threadsMax);     // final int          threadsMax,       // maximal number of threads to launch
		
		if (tp_tasks.length == 0) {
			System.out.println("Empty tasks - nothing to do");
			return;
		}
		//texture_tiles
		final boolean fneed_macro = need_macro;
		final boolean fneed_corr =  need_corr && used_corrs[0]; // *** tasks should include correlation

		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > 2);

		final float [] lpf_flat = floatGetCltLpfFd(gpu_sigma_corr);

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				globalDebugLevel > 2);

		final float [] lpf_rb_flat = floatGetCltLpfFd(gpu_sigma_rb_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				globalDebugLevel > 2);

		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				use_main,                  // use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify

		gpuQuad.execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers

		if ((fdisp_dist != null) || (fpxpy != null)) { // skip
			final TpTask[] tp_tasks_full = gpuQuad.getTasks(use_main);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int indx_tile = ai.getAndIncrement(); indx_tile < tp_tasks_full.length; indx_tile = ai.getAndIncrement()) {
							TpTask task = tp_tasks_full[indx_tile];
							if (fdisp_dist != null) {
								fdisp_dist[task.getTileY()][task.getTileX()] = task.getDispDist();
							}
							if (fpxpy != null) {
								fpxpy[task.getTileY()][task.getTileX()] = task.getXY(); // use_main); // boolean use_aux);
							}
						} // end of tile
					}
				};
			}
			startAndJoin(threads);
			ai.set(0);
		}
		
		
		gpuQuad.execConvertDirect(-1); // boolean erase_clt
		if (iclt_fimg != null) {
			gpuQuad.execImcltRbgAll(isMonochrome());  // execute GPU kernel
			for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
				iclt_fimg[ncam] = gpuQuad.getRBG(ncam); // retrieve data from GPU (not used !)
			}
		} else {gpuQuad.execImcltRbgAll(isMonochrome());} // just for testing
		// does it need texture tiles to be output?
		if (texture_img != null) {
			Rectangle woi = new Rectangle(); // will be filled out to match actual available image
			gpuQuad.execRBGA(
					col_weights,                   // double [] color_weights,
					isLwir(),         // boolean   is_lwir,
					min_shot,       // double    min_shot,           // 10.0
					scale_shot,     // double    scale_shot,         // 3.0
					diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
					diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove);   // boolean   dust_remove,
			float [][] rbga = gpuQuad.getRBGA(
					(isMonochrome() ? 1 : 3), // int     num_colors,
					(texture_woi_pix != null)? texture_woi_pix : woi);
			for (int ncol = 0; ncol < texture_img.length; ncol++) if (ncol < rbga.length) {
				texture_img[ncol] = rbga[ncol];
			}
		}
		// does it need macro data?
		if (fneed_macro) {
			//Generate non-overlapping (16x16) texture tiles, prepare 
			gpuQuad.execTextures(
				col_weights,                   // double [] color_weights,
				isLwir(),         // boolean   is_lwir,
				min_shot,       // double    min_shot,           // 10.0
				scale_shot,     // double    scale_shot,         // 3.0
				diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				diff_threshold, // double    diff_threshold,     // pixel value/pixel change - never used in GPU ?
				min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				dust_remove,    // boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
				false,          // boolean   calc_textures,
				true,           // boolean   calc_extra
				false);         // boolean   linescan_order) // TODO: use true to avoid reordering of the low-res output 
			float [][] extra = gpuQuad.getExtra(); // now 4*numSensors
//			int num_cams = gpuQuad.getNumCams();
			int num_cams = getNumSensors();
			for (int ncam = 0; ncam < num_cams; ncam++) {
				int indx = ncam + IMG_DIFF0_INDEX;
//				if ((disparity_modes & (1 << indx)) != 0){
				if (needImgDiffs(disparity_modes)){
					disparity_map[indx] = new double [extra[ncam].length];
					for (int i = 0; i < extra[ncam].length; i++) {
						disparity_map[indx][i] = extra[ncam][i];
					}
				}
			}
			for (int nc = 0; nc < (extra.length - num_cams); nc++) {
				int sindx = nc + num_cams;
				/*
				int indx = nc + IMG_TONE_RGB;
				if ((disparity_modes & (1 << indx)) != 0){
					disparity_map[indx] = new double [extra[sindx].length];
					for (int i = 0; i < extra[sindx].length; i++) {
						disparity_map[indx][i] = extra[sindx][i];
					}
				}
	            */
				int indx = nc + getImgToneRGB(); // IMG_TONE_RGB;
//				if ((disparity_modes & (1 << indx)) != 0){
				if (needTonesRGB(disparity_modes)){
					disparity_map[indx] = new double [extra[sindx].length];
					for (int i = 0; i < extra[sindx].length; i++) {
						disparity_map[indx][i] = extra[sindx][i];
					}
				}

				
			}			
		}
		// does it need non-overlapping texture tiles
		if (texture_tiles != null) {   // compatible with the CPU ones
			//Generate non-overlapping (16x16) texture tiles, prepare 
			gpuQuad.execTextures(
				col_weights,                   // double [] color_weights,
				isLwir(),         // boolean   is_lwir,
				min_shot,       // double    min_shot,           // 10.0
				scale_shot,     // double    scale_shot,         // 3.0
				diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				diff_threshold, // double    diff_threshold,     // pixel value/pixel change - never used in GPU ?
				min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				dust_remove,    // boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
				true,           // boolean   calc_textures,
				false,          // boolean   calc_extra
				false);         // boolean   linescan_order) 

			int [] texture_indices =  gpuQuad.getTextureIndices();
			int    num_src_slices =   numcol + 1; //  + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
			float [] flat_textures =  gpuQuad.getFlatTextures( // fatal error has been detected by the Java Runtime Environment:
					texture_indices.length,
					numcol,    // int     num_colors,
		    		false);    // clt_parameters.keep_weights); // boolean keep_weights);

			gpuQuad.doubleTextures(
		    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
		    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
		    		texture_indices,                     // int []       indices,
		    		flat_textures,                       // float [][][] ftextures,
		    		tilesX,                              // int          full_width,
		    		isMonochrome()? 2: 4,                // rbga only /int          num_slices
		    		num_src_slices                       // int          num_src_slices
		    		);
		}
		
		
		// does it need correlations?
		if (fneed_corr) {
			int mcorr_sel = Correlation2d.corrSelEncode(imgdtt_params,numSensors);
			int [] i_mcorr_sel = Correlation2d.intCorrPairs(
					mcorr_sel,
					numSensors,
					4); // int num_out); should be 4 int
			//Generate 2D phase correlations from the CLT representation
			gpuQuad.execCorr2D_TD(col_weights,i_mcorr_sel); // Get TD version of correlations (may be read out and saved) 
			final int [] corr_indices = gpuQuad.getCorrIndices();
			if (fcorr_td != null) {
				gpuQuad.getCorrTilesTd(fcorr_td); // generate transform domain correlation pairs
			}
			
			gpuQuad.execCorr2D_normalize(
	        		false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
			
			// Combine 4 ortho pairs
			final int num_pairs = Correlation2d.getNumPairs(numSensors);
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x0f); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 0) && (fcorr_combo_td[0] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[0]); // generate transform domain correlation pairs for quad ortho combination
			}
			// normalize and convert to pixel domain
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final int [] corr_quad_indices = gpuQuad.getCorrComboIndices(); // get quad
			final float [][] fcorr2D_quad =   gpuQuad.getCorr2DCombo(gpu_corr_rad);

			// Combine 2 diagonal pairs			
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x30); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 1) && (fcorr_combo_td[1] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[1]); // generate transform domain correlation pairs for cross diagonal combination
			}
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_cross =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			
			// Combine 2 horizontal pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x03); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 2) && (fcorr_combo_td[2] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[2]); // generate transform domain correlation pairs for horizontal combination
			}
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_hor =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			// Combine 2 vertical pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
					0x0c, // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
					true); // boolean       no_transpose_vertical
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 3) && (fcorr_combo_td[3] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[3]); // generate transform domain correlation pairs for vertical combination
			}
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_vert =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			
			
			
			if (corr_indices.length > 0) {
				final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size
				// assuming that the correlation pairs sets are the same for each tile that has correlations
				// find this number
				int nt0 = (corr_indices[0] >> GPUTileProcessor.CORR_NTILE_SHIFT);
				int nc0 = 1;
				for (int i = 1; (i < corr_indices.length) && ((corr_indices[i] >> GPUTileProcessor.CORR_NTILE_SHIFT) == nt0) ; i++) {
					nc0++;
				}
				final int num_tile_corr = nc0; // normally 6
				final int num_tiles = corr_indices.length / num_tile_corr; 
				

				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						@Override
						public void run() {
							//							int tileY,tileX,tIndex;
//							double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
							
							Correlation2d corr2d = new Correlation2d(
									numSensors,
									imgdtt_params,              // ImageDttParameters  imgdtt_params,
									transform_size,             // int transform_size,
									2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
									isMonochrome(), // boolean monochrome,
									(globalDebugLevel > -1));   //   boolean debug)
							corr2d.createOrtoNotch(
									imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
									imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
									(imgdtt_params.lma_debug_level > 1)); // boolean debug);

							for (int indx_tile = ai.getAndIncrement(); indx_tile < num_tiles; indx_tile = ai.getAndIncrement()) {
								// double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
								// added quad and cross combos
								double [][]  corrs = new double [num_pairs + 4][corr_length]; // 225-long (15x15)
								int indx_corr = indx_tile * num_tile_corr;
								int nt = (corr_indices[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
								int tileX = nt % tilesX;
								int tileY = nt / tilesX;
								int tIndex = tileY * tilesX + tileX;
								
								// Prepare the same (currently 10-layer) corrs as double [][], as in CPU version
								int pair_mask = 0;
								for (int indx_pair = 0; indx_pair < num_tile_corr; indx_pair++) {
									int pair = corr_indices[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
									assert pair < num_pairs : "invalid correllation pair";
									pair_mask |= (1 << pair);
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D[indx_corr][i]; // from float to double
									}
									indx_corr++; 
								}
								// add 4 combo layers : quad, cross, hor, vert
								int pair = num_pairs; // 6
								nt = (corr_quad_indices[indx_tile] >> GPUTileProcessor.CORR_NTILE_SHIFT); // corr_quad_indices - different sequence
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_quad[indx_tile][i]; // from float to double
								}
								// indices for cross are the same as for quad
								pair++;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_cross[indx_tile][i]; // from float to double
								}
								
								// indices for hor are the same as for quad
								pair++;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_hor[indx_tile][i]; // from float to double
								}

								// indices for vert are the same as for quad
								pair++;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_vert[indx_tile][i]; // from float to double
								}
								
								int used_pairs = pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
								boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
								
								corr_common_GPU(
										imgdtt_params,        // final ImageDttParameters  imgdtt_params,
										clt_corr_partial,     // final double [][][][][]   clt_corr_partial,			
										used_pairs,           // final int           used_pairs,
										disparity_map,        // final double [][]   disparity_map,
										clt_mismatch,         // final double [][]   clt_mismatch,
										saturation_imp,       // final boolean [][]  saturation_imp,
										fneed_macro,          // final boolean       fneed_macro,
										corr2d,               // final Correlation2d corr2d,
										corrs,                // final double [][]   corrs,
										tileX,                // final int           tileX,
										tileY,                // final int           tileY,
										max_corr_radius,      // final double        max_corr_radius, // 3.9;
										tile_lma_debug_level, // int                 tile_lma_debug_level,
										debugTile,            // boolean             debugTile,
										globalDebugLevel);    // final int           globalDebugLevel)							
								if ((disparity_map != null) && Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
									System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
								}
								// only debug is left
								// old (per-color correlation)
								// removed
								
							} // end of tile
						}
					};
				}
				startAndJoin(threads);
			} else {
				// no correlation tiles to process
			}
		}
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
	}

	// calculate FD correlations for specified pXpYD array to use for interscene accumulation
	@Deprecated
	public void quadCorrTD(
			final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
   			final double [][]         pXpYD,            // per-tile array of pX,pY,disparity triplets (or nulls)
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: quad, cross, hor,vert
			final GeometryCorrection  geometryCorrection,
			final double              disparity_corr,   // disparity offset at infinity
			final int                 margin,           // do not use tiles if their centers are closer to the edges
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final double              gpu_sigma_rb_corr,    //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 :
			final double              gpu_sigma_corr,       //  =    0.9;gpu_sigma_corr_m
			final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
			final double              corr_red, // +used
			final double              corr_blue,// +used
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		// prepare tasks
       	TpTask[] tp_tasks =  gpuQuad.setInterTasks(
       			false, // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
       			pXpYD,              // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				null,               // final boolean []          selection, // may be null, if not null do not  process unselected tiles
       			geometryCorrection, // final GeometryCorrection  geometryCorrection,
       			disparity_corr,     // final double              disparity_corr,
       			margin,             // final int                 margin,      // do not use tiles if their centers are closer to the edges
    			null,               // final boolean []          valid_tiles,            
    			threadsMax);        // final int                 threadsMax)  // maximal number of threads to launch
       	quadCorrTD(
    			imgdtt_params,    // Now just extra correlation parameters, later will include, most others
    			tp_tasks,
    			fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
    			fcorr_combo_td,  // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: quad, cross, hor,vert
    			geometryCorrection,
    			margin,           // do not use tiles if their centers are closer to the edges
    			gpu_sigma_r,     // 0.9, 1.1
    			gpu_sigma_b,     // 0.9, 1.1
    			gpu_sigma_g,     // 0.6, 0.7
    			gpu_sigma_m,     //  =       0.4; // 0.7;
    			gpu_sigma_rb_corr,    //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 :
    			gpu_sigma_corr,       //  =    0.9;gpu_sigma_corr_m
    			gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
    			corr_red, // +used
    			corr_blue,// +used
    			threadsMax,       // maximal number of threads to launch
    			globalDebugLevel);
	}

	// Run after quadCorrTD when clt is in memory
	public double [][][][] process_texture_tiles(
				double corr_red,
				double corr_blue,
				double    min_shot,           // 10.0
				double    scale_shot,         // 3.0
				double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				double    diff_threshold,     // pixel value/pixel change
				double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				boolean   dust_remove
			){
		int numcol = isMonochrome()? 1 : 3;
		double [] col_weights = new double[numcol];
		if (isMonochrome()) {
			col_weights[0] = 1.0;
		} else {
			col_weights[2] = 1.0/(1.0 +  corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}
		
		gpuQuad.execTextures(
				col_weights,    // double [] color_weights,
				isLwir(),       // boolean   is_lwir,
				min_shot,       // double    min_shot,           // 10.0
				scale_shot,     // double    scale_shot,         // 3.0
				diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				dust_remove,    // boolean   dust_remove,
				true,           // boolean   calc_textures,
				false,          // boolean   calc_extra
				false);         // boolean   linescan_order)
		int    num_src_slices = numcol + 1; //  + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
		int [] texture_indices = gpuQuad.getTextureIndices();
		float [] flat_textures =  gpuQuad.getFlatTextures(
				texture_indices.length,
				numcol, // int     num_colors,
	    		false); // clt_parameters.keep_weights); // boolean keep_weights);
		int tilesX = gpuQuad.img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY = gpuQuad.img_height / GPUTileProcessor.DTT_SIZE;
		double [][][][] texture_tiles = new double [tilesY][tilesX][][];
		gpuQuad.doubleTextures(
	    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
	    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
	    		texture_indices,                     // int []       indices,
	    		flat_textures,                       // float [][][] ftextures,
	    		tilesX,                              // int          full_width,
	    		isMonochrome()? 2: 4,                // rbga only /int          num_slices Same number
	    		num_src_slices                       // int          num_src_slices
	    		);
		return texture_tiles;
	}
	
	public float [][] get_diffs_lowres(
			double      corr_red,
			double      corr_blue,
			double      min_shot,           // 10.0
			double      scale_shot,         // 3.0
			double      diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
			double      diff_threshold,     // pixel value/pixel change
			double      min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean     dust_remove,
			boolean     save_diff,
			boolean     save_lowres,
			double [][] disparity_map   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate			
			){
		if (gpuQuad.num_task_tiles == 0) {
			System.out.println("get_diffs_lowres(): num_task_tiles=0!");
			return null;
		}

		int numcol = isMonochrome()? 1 : 3;
		double [] col_weights = new double[numcol];
		if (isMonochrome()) {
			col_weights[0] = 1.0;
		} else {
			col_weights[2] = 1.0/(1.0 +  corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}
		gpuQuad.execTextures( // CUDA_ERROR_INVALID_VALUE
				col_weights,    // double [] color_weights,
				isLwir(),       // boolean   is_lwir,
				min_shot,       // double    min_shot,           // 10.0
				scale_shot,     // double    scale_shot,         // 3.0
				diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				dust_remove,    // boolean   dust_remove,
				false,          // boolean   calc_textures,
				true,           // boolean   calc_extra
				false);         // boolean   linescan_order)
		float [][] extra = gpuQuad.getExtra(); // now 4*numSensors
		int ports_rgb_len = numSensors*numcol;  // 12 / 16

		if (disparity_map != null) { // convert to double and copy to the correct layers of disparity_map
//			int num_cams = getNumSensors();
			if (disparity_map != null){ // Only init layers that are going to be used, keep others
				if (save_diff && (disparity_map.length >= (IMG_DIFF0_INDEX + numSensors))) {
					for (int i = 0; i < numSensors; i++) {
						int i_dest = IMG_DIFF0_INDEX + i;
						disparity_map[i_dest] = new double[extra[i].length];
						for (int j = 0; j < extra[i].length; j++) {
							disparity_map[i_dest][j] = extra[i][j];
						}
					}
				}
				if (save_lowres && (disparity_map.length >= (getImgToneRGB() + ports_rgb_len))) {
		            for (int i = 0; i < ports_rgb_len; i++){
		            	int i_src = numSensors + i;
						int i_dest = getImgToneRGB() + i;
		                disparity_map[i_dest] = new double[extra[i_src].length];
						for (int j = 0; j < extra[i].length; j++) {
							disparity_map[i_dest][j] = extra[i_src][j];
						}
		            }
				}
			}
		}
		return extra;
	}
	
	public void setUpdateTasksGPU(
			final TpTask[]            tp_tasks
			) {
		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				false,                     // use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify
		// FIXME: change back to false !!!!
		// Testing, remove when done
//		gpuQuad.resetGeometryCorrection();
//		gpuQuad.setConvolutionKernels(true); // set kernels if they are not set already
//		gpuQuad.setBayerImages(true);     // set Bayer images if this.quadCLT instance has new ones
		// Why always NON-UNIFORM grid? Already set in tp_tasks
		gpuQuad.execSetTilesOffsets(false); // false); // prepare tiles offsets in GPU memory, using NON-UNIFORM grid (pre-calculated)
		// update tp_tasks
		gpuQuad.updateTasks(
				tp_tasks,
				false); // boolean use_aux    // while is it in class member? - just to be able to free
	
		
	}
	
	public void quadCorrTD(
			final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
			final TpTask[]            tp_tasks,
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
//			final GeometryCorrection  geometryCorrection, // not used !!!
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final double              gpu_sigma_rb_corr,    //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 :
			final double              gpu_sigma_corr,       //  =    0.9;gpu_sigma_corr_m
			final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
			final double              corr_red, // +used
			final double              corr_blue,// +used
			final int                 mcorr_sel,    // Which pairs to correlate // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert. 0 - no correlations
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int numcol = isMonochrome()?1:3;
		final double [] col_weights= new double [numcol]; // colors are RBG
		if (isMonochrome()) {
			col_weights[0] = 1.00; // Was 0 ! 0;
		} else {
			col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}

		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > 2);

		final float [] lpf_flat = floatGetCltLpfFd(gpu_sigma_corr);

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				globalDebugLevel > 2);

		final float [] lpf_rb_flat = floatGetCltLpfFd(gpu_sigma_rb_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				globalDebugLevel > 2);
		
		final float [] log_flat = floatGetCltHpfFd(gpu_sigma_log_corr);
		if (globalDebugLevel < -100) {
			double dbg_sum = 0.0;
			for (int i = 0; i < log_flat.length; i++) dbg_sum +=log_flat[i];
			System.out.println("dbg_sum("+gpu_sigma_log_corr+")="+dbg_sum);
			(new ShowDoubleFloatArrays()).showArrays(
					log_flat,
					8,
					8,
					"hpf_"+gpu_sigma_log_corr);
			final float [] log_flat0 = floatGetCltHpfFd(4.0);
			dbg_sum = 0.0;
			for (int i = 0; i < log_flat.length; i++) dbg_sum +=log_flat0[i];
			System.out.println("dbg_sum("+4.0+")="+dbg_sum);
			(new ShowDoubleFloatArrays()).showArrays(
					log_flat0,
					8,
					8,
					"hpf_"+4.0);
			final float [] log_flat1 = floatGetCltHpfFd(1.0);
			dbg_sum = 0.0;
			for (int i = 0; i < log_flat.length; i++) dbg_sum +=log_flat1[i];
			System.out.println("dbg_sum("+1.0+")="+dbg_sum);
			(new ShowDoubleFloatArrays()).showArrays(
					log_flat1,
					8,
					8,
					"hpf_"+1.0);
			System.out.println("dbg_sum("+1.0+")="+dbg_sum);
		}
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"LoG_corr", // String const_name, // "lpf_corr"
				log_flat,
				globalDebugLevel > -1);

		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				false,                     // use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify
		// used alternative method to prepare tasks, not centered in the tile centers
		// FIXME: change back to false !!!!
		// Testing, remove when done
		// gpuQuad.resetGeometryCorrection();
		// gpuQuad.setConvolutionKernels(true); // set kernels if they are not set already
		// gpuQuad.setBayerImages(true);     // set Bayer images if this.quadCLT instance has new ones
		
		// Why always NON-UNIFORM grid? Already set in tp_tasks
		gpuQuad.execSetTilesOffsets(false); // false); // prepare tiles offsets in GPU memory, using NON-UNIFORM grid (pre-calculated)
		// update tp_tasks
		gpuQuad.updateTasks(
				tp_tasks,
				false); // boolean use_aux    // while is it in class member? - just to be able to free
		
       	// Skipping if ((fdisp_dist != null) || (fpxpy != null)) {...

		gpuQuad.execConvertDirect(-1); // boolean erase_clt
		if (mcorr_sel == 0) { // no correlation at all
			return;
		}

		int [] i_mcorr_sel = Correlation2d.intCorrPairs(
				mcorr_sel,
				numSensors,
				4); // int num_out); should be 4 int
		gpuQuad.setCorrMask(i_mcorr_sel);
		
		boolean test_execCorr2D = false;
		if (test_execCorr2D) {
			double fat_zero = 2000.0; // 30.0; // 2000.0; // 30.00;
			int gpu_corr_rad =  7;
			int corr_size = 2* gpu_corr_rad + 1;
			//Generate 2D phase correlations from the CLT representation
			
    // Try to zero out memory before calculating?	
//			gpuQuad.eraseGpuCorrs();
			gpuQuad.execCorr2D(
					gpuQuad.getCorrMask(),    // boolean [] pair_select,
					col_weights,              // double [] scales,
					fat_zero,                 // double fat_zero);
					gpu_corr_rad);            // int corr_radius
			// Should be done before execCorr2D_TD as corr_indices are shared to save memory
			int [] corr_indices = gpuQuad.getCorrIndices();
			// the following is not yet shared
			float [][] corr2D =  gpuQuad.getCorr2D(
					gpu_corr_rad); //  int corr_rad);
			final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
			final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
			int num_tiles = tilesX * tilesY;
//			int with =   tilesX * GPUTileProcessor.DTT_SIZE;
//			int height = tilesY * GPUTileProcessor.DTT_SIZE;
			int sq = 16;
			int num_pairs = gpuQuad.getNumUsedPairs();
			float [][] corr_img = new float [num_pairs][tilesY * sq * tilesX * sq];
			for (int pair = 0; pair < num_pairs; pair++) {
				Arrays.fill(corr_img[pair], Float.NaN);
			}
		    for (int ict = 0; ict < corr_indices.length; ict++){
		    	//    	int ct = cpu_corr_indices[ict];
		    	int ctt = ( corr_indices[ict] >>  GPUTileProcessor.CORR_NTILE_SHIFT);
		    	int cpair = corr_indices[ict] & ((1 << GPUTileProcessor.CORR_NTILE_SHIFT) - 1);
		    	int ty = ctt / tilesX;
		    	int tx = ctt % tilesX;
		    	int dst_offs0 = (ty * sq * tilesX * sq) + (tx * sq);
		    	for (int iy = 0; iy < corr_size; iy++){
		    		int dst_offs = dst_offs0 + iy * (tilesX * sq);
		    		System.arraycopy(corr2D[ict], iy * corr_size, corr_img[cpair], dst_offs, corr_size);
		    	}
		    }
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					corr_img,
					tilesX * sq,
					tilesY * sq,
					true,
					"test-corr");
		}

		//Generate 2D phase correlations from the CLT representation
		gpuQuad.execCorr2D_TD(col_weights, i_mcorr_sel); // Get TD version of correlations 
//		final int [] corr_indices = gpuQuad.getCorrIndices();
		if (fcorr_td != null) {
			gpuQuad.getCorrTilesTd(fcorr_td); // generate transform domain correlation pairs
		}
//		// FIXME: will not work with combining pairs !!!
//		final int num_pairs = Correlation2d.getNumPairs(numSensors);
	}
	
    /**
     * Correlate two scenes - reference (should be set with setReferenceTD() ) and this one, keep results in TD
     * results include selected sensors and the sum of them
     * @param imgdtt_params
     * @param tp_tasks (tasks should not include the tiles that are missing from the reference scene)
     * @param fcorr_td null or float [tilesY][tilesX][][] - will return [number_of_selected_sensors + 1][256] for non-empty
     * @param gpu_sigma_r
     * @param gpu_sigma_b
     * @param gpu_sigma_g
     * @param gpu_sigma_m
     * @param gpu_sigma_rb_corr
     * @param gpu_sigma_corr
     * @param gpu_sigma_log_corr
     * @param corr_red
     * @param corr_blue
     * @param sensor_mask_inter The bitmask - which sensors to correlate, -1 - all.
     * @param threadsMax
     * @param globalDebugLevel
     */
	public void interCorrTD(
			final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
			final TpTask[]            tp_tasks,
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final double              gpu_sigma_rb_corr,    //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 :
			final double              gpu_sigma_corr,       //  =    0.9;gpu_sigma_corr_m
			final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
			final double              corr_red, // +used
			final double              corr_blue,// +used
			final int                 sensor_mask_inter, // The bitmask - which sensors to correlate, -1 - all.
			final int                 threadsMax,        // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int numcol = isMonochrome()?1:3;
		final double [] col_weights= new double [numcol]; // colors are RBG
		if (isMonochrome()) {
			col_weights[0] = 1.00; // Was 0 ! 0;
		} else {
			col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}

		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > 2);

		final float [] lpf_flat = floatGetCltLpfFd(gpu_sigma_corr);

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				globalDebugLevel > 2);

		final float [] lpf_rb_flat = floatGetCltLpfFd(gpu_sigma_rb_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				globalDebugLevel > 2);
		
		final float [] log_flat = floatGetCltHpfFd(gpu_sigma_log_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"LoG_corr", // String const_name, // "lpf_corr"
				log_flat,
				globalDebugLevel > 2);

		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				false,                     // use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify
		// used alternative method to prepare tasks, not centered in the tile centers
		// FIXME: change back to false !!!!
		// Testing, remove when done
		// gpuQuad.resetGeometryCorrection();
		// gpuQuad.setConvolutionKernels(true); // set kernels if they are not set already
		// gpuQuad.setBayerImages(true);     // set Bayer images if this.quadCLT instance has new ones
		
		// Why always NON-UNIFORM grid? Already set in tp_tasks
		gpuQuad.execSetTilesOffsets(false); // false); // prepare tiles offsets in GPU memory, using NON-UNIFORM grid (pre-calculated)
		// update tp_tasks
		gpuQuad.updateTasks(
				tp_tasks,
				false); // boolean use_aux    // while is it in class member? - just to be able to free
		
       	// Skipping if ((fdisp_dist != null) || (fpxpy != null)) {...

		gpuQuad.execConvertDirect(-1); // boolean erase_clt
		if (sensor_mask_inter == 0) { // no correlation at all
			return;
		}
		gpuQuad.setSensorMaskInter(sensor_mask_inter);
		//Generate 2D phase correlations from the CLT representation
		gpuQuad.execCorr2D_inter_TD(
				col_weights); // double [] scales,
		if (fcorr_td != null) {
			gpuQuad.getCorrTilesTd(
					true, //boolean        inter,
					fcorr_td); // generate transform domain correlation pairs
		}
		return;
	}
	
	/**
	 * Convert reference scene to FD and save result in extra GPU array for the future interscene correlation
	 * Geometry correction and images will come from gpuQuad instance - 
	 * @param erase_clt erase CLT (<0 - do not erase, 0 - erase to 0.0, >0 - erase to NaN). Needed only for later IMCLT
	 *                  end rendering images. NaN produces sharp, distinct borders; 0f - blended
	 * @param wh if null, will uses sensor dimensions. Otherwise {width, height} in pixels
	 * @param imgdtt_params
	 * @param use_reference_buffer true - use extra GPU array, false - use main one
	 * @param tp_tasks
	 * @param gpu_sigma_r
	 * @param gpu_sigma_b
	 * @param gpu_sigma_g
	 * @param gpu_sigma_m
	 * @param threadsMax
	 * @param globalDebugLevel
	 */
	public void setReferenceTD(
			final int                 erase_clt,
			final int []              wh,               // null (use sensor dimensions) or pair {width, height} in pixels
			final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
			final boolean             use_reference_buffer,
			final TpTask[]            tp_tasks,
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > 2);
		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				false,                     // use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify
		// Why always NON-UNIFORM grid? Already set in tp_tasks
		gpuQuad.execSetTilesOffsets(false); // false); // prepare tiles offsets in GPU memory, using NON-UNIFORM grid (pre-calculated)
		// update tp_tasks
		gpuQuad.updateTasks(
				tp_tasks,
				false); // boolean use_aux    // while is it in class member? - just to be able to free
		// Skipping if ((fdisp_dist != null) || (fpxpy != null)) {...
//		int [] wh = null;
//		int  erase_clt = 1; // NaN;
		gpuQuad.execConvertDirect(use_reference_buffer, wh, erase_clt); // put results into a "reference" buffer
	}



	
	
	
	@Deprecated
	public void quadCorrTD(
			final ImageDttParameters  imgdtt_params,    // Now just extra correlation parameters, later will include, most others
			final TpTask[] tp_tasks,
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: quad, cross, hor,vert
			final GeometryCorrection  geometryCorrection,
//			final double              disparity_corr,   // disparity offset at infinity
			final int                 margin,           // do not use tiles if their centers are closer to the edges
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final double              gpu_sigma_rb_corr,    //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 :
			final double              gpu_sigma_corr,       //  =    0.9;gpu_sigma_corr_m
			final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
			final double              corr_red, // +used
			final double              corr_blue,// +used
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int numcol = isMonochrome()?1:3;
		final double [] col_weights= new double [numcol]; // colors are RBG
		if (isMonochrome()) {
			col_weights[0] = 0;
//			col_weights[2] = 1.0;// green color/mono
//			col_weights[0] = 0;
//			col_weights[1] = 0;
		} else {
			col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
			col_weights[0] = corr_red *  col_weights[2];
			col_weights[1] = corr_blue * col_weights[2];
		}
		
		//texture_tiles

		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > 2);

		final float [] lpf_flat = floatGetCltLpfFd(gpu_sigma_corr);

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				globalDebugLevel > 2);

		final float [] lpf_rb_flat = floatGetCltLpfFd(gpu_sigma_rb_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				globalDebugLevel > 2);
		
		final float [] log_flat = floatGetCltHpfFd(gpu_sigma_log_corr);
		if (globalDebugLevel < -100) {
			double dbg_sum = 0.0;
			for (int i = 0; i < log_flat.length; i++) dbg_sum +=log_flat[i];
			System.out.println("dbg_sum("+gpu_sigma_log_corr+")="+dbg_sum);
			(new ShowDoubleFloatArrays()).showArrays(
					log_flat,
					8,
					8,
					"hpf_"+gpu_sigma_log_corr);
			final float [] log_flat0 = floatGetCltHpfFd(4.0);
			dbg_sum = 0.0;
			for (int i = 0; i < log_flat.length; i++) dbg_sum +=log_flat0[i];
			System.out.println("dbg_sum("+4.0+")="+dbg_sum);
			(new ShowDoubleFloatArrays()).showArrays(
					log_flat0,
					8,
					8,
					"hpf_"+4.0);
			final float [] log_flat1 = floatGetCltHpfFd(1.0);
			dbg_sum = 0.0;
			for (int i = 0; i < log_flat.length; i++) dbg_sum +=log_flat1[i];
			System.out.println("dbg_sum("+1.0+")="+dbg_sum);
			(new ShowDoubleFloatArrays()).showArrays(
					log_flat1,
					8,
					8,
					"hpf_"+1.0);
			System.out.println("dbg_sum("+1.0+")="+dbg_sum);
		}
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"LoG_corr", // String const_name, // "lpf_corr"
				log_flat,
				globalDebugLevel > -1);
		

		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				false,                     //  use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify
		// used alternative method to prepare tasks, not centered in the tile centers
///		gpuQuad.execSetTilesOffsets(); // prepare tiles offsets in GPU memory
       	// Skipping if ((fdisp_dist != null) || (fpxpy != null)) {...
       	
		gpuQuad.execConvertDirect(-1); // boolean erase_clt

		//Generate 2D phase correlations from the CLT representation
		int mcorr_sel = Correlation2d.corrSelEncode(imgdtt_params,numSensors);
		int [] i_mcorr_sel = Correlation2d.intCorrPairs(
				mcorr_sel,
				numSensors,
				4); // int num_out); should be 4 int
		gpuQuad.execCorr2D_TD(col_weights,i_mcorr_sel); // Get TD version of correlations 
		final int [] corr_indices = gpuQuad.getCorrIndices();
		if (fcorr_td != null) {
			gpuQuad.getCorrTilesTd(fcorr_td); // generate transform domain correlation pairs
		}
		// FIXME: will not work with combining pairs !!!
		final int num_pairs = Correlation2d.getNumPairs(numSensors);
		// Combine 4 ortho pairs
		if (fcorr_combo_td != null) {
			if ((fcorr_combo_td.length >= 0) && (fcorr_combo_td[0] != null)) {
				gpuQuad.execCorr2D_combine( // calculate cross pairs
						true, // boolean init_corr,    // initialize output tiles (false - add to current)
						num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
						0x0f); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[0]); // generate transform domain correlation pairs for quad ortho combination
			}
			// Combine 2 diagonal pairs			
			if ((fcorr_combo_td.length >= 1) && (fcorr_combo_td[1] != null)) {
				gpuQuad.execCorr2D_combine( // calculate cross pairs
						true, // boolean init_corr,    // initialize output tiles (false - add to current)
						num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
						0x30); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[1]); // generate transform domain correlation pairs for cross diagonal combination
			}
			// Combine 2 horizontal pairs
			if ((fcorr_combo_td.length >= 2) && (fcorr_combo_td[2] != null)) {
				gpuQuad.execCorr2D_combine( // calculate cross pairs
						true, // boolean init_corr,    // initialize output tiles (false - add to current)
						num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
						0x03); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[2]); // generate transform domain correlation pairs for horizontal combination
			}
			// Combine 2 vertical pairs
			if ((fcorr_combo_td.length >= 3) && (fcorr_combo_td[3] != null)) {
				gpuQuad.execCorr2D_combine( // calculate cross pairs
						true, // boolean init_corr,    // initialize output tiles (false - add to current)
						num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
						0x0c, // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
						true); // boolean       no_transpose_vertical
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[3]); // generate transform domain correlation pairs for vertical combination
			}
		}
	}
	
	
	@Deprecated
	public TpTask[][] clt_aberrations_quad_corr_GPU_test( // Not used
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
			                                           // each of the top elements may be null to skip particular combo type
			final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
			final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
			final float  [][][][]     fpxpy_test,      // [tilesY][tilesX][cams][2], tile {pX,pY}
			final double [][][][]     clt_corr_combo,  // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			                                           // [type][tilesY][tilesX] should be set by caller
													   // types: 0 - selected correlation (product+offset), 1 - sum
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
                                                       // [tilesY][tilesX] should be set by caller
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
			final double [][][][]     texture_tiles,   // compatible with the CPU ones      
			final float [][]          texture_img,     // null or [3][] (RGB) or [4][] RGBA
			final Rectangle           texture_woi_pix, // null or generated texture location/size in pixels
			final float [][][]        iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
			// new parameters, will replace some other?
			final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
			final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
			final double              gpu_sigma_r,     // 0.9, 1.1
			final double              gpu_sigma_b,     // 0.9, 1.1
			final double              gpu_sigma_g,     // 0.6, 0.7
			final double              gpu_sigma_m,     //  =       0.4; // 0.7;
			final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
			final double              gpu_sigma_corr,   //  =    0.9;gpu_sigma_corr_m
			final int                 gpu_corr_rad,     // = transform_size - 1 ?
			final double              corr_red, // +used
			final double              corr_blue,// +used
			final double              max_corr_radius, // 3.9;
			final int                 corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
			
			final double              min_shot,        // +10.0;  // Do not adjust for shot noise if lower than
			final double              scale_shot,      // +3.0;   // scale when dividing by sqrt ( <0 - disable correction)
			final double              diff_sigma,      // +5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
			final double              diff_threshold,  // +5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
			final boolean             diff_gauss,      // true;  // when averaging images, use Gaussian around average as weight (false - sharp all/nothing)
			final double              min_agree,       // +3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			final boolean             dust_remove,     // +Do not reduce average weight when only one image differs much from the average
			final GeometryCorrection  geometryCorrection, // for GPU TODO: combine geometry corrections if geometryCorrection_main is not null
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final int                 window_type,    // GPU: will not be used
			final double              disparity_corr, // disparity at infinity
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,       // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return null;
		}
		final boolean [][] saturation_imp = gpuQuad.quadCLT.saturation_imp;               // boolean [][] saturation_imp, // (near) saturated pixels or null
//gpuQuad
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;

//		final double [][] debug_offsets = new double[imgdtt_params.lma_dbg_offset.length][2];
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}


		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int quad = 4;   // number of subcameras
//		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green
		final int numcol = isMonochrome()?1:3;

		final int width =  gpuQuad.getImageWidth();
		final int height = gpuQuad.getImageHeight();
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
////		final int nTilesInChn=tilesX*tilesY;
////		final double [][][][][][] clt_data = new double[quad][numcol][tilesY][tilesX][][];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*quad][tilesX*tilesY]) : null;
		// not yet used with GPU
/**		
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}
*/
		// keep for now for mono, find out  what do they mean for macro mode
		
		if (macro_mode) { // all the same as they now mean different
			//compensating Bayer correction
			col_weights[0] = 0.25; //  1.0/3;
			col_weights[1] = 0.25; //  1.0/3;
			col_weights[2] = 0.5; // 1.0/3;
		} else {
			if (isMonochrome()) {
				col_weights[2] = 1.0;// green color/mono
				col_weights[0] = 0;
				col_weights[1] = 0;
			} else {
				col_weights[2] = 1.0/(1.0 + corr_red + corr_blue);    // green color
				col_weights[0] = corr_red *  col_weights[2];
				col_weights[1] = corr_blue * col_weights[2];
			}
		}
		final int corr_size = transform_size * 2 - 1;
		
		if ((globalDebugLevel > -10) && (disparity_corr != 0.0)){
			System.out.println(String.format("Using manual infinity disparity correction of %8.5f pixels",disparity_corr));
		}

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 1){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}


		if (globalDebugLevel > 1) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		
	
		// add optional initialization of debug layers here
		boolean need_macro = false;
		boolean need_corr = (clt_mismatch != null) || (fcorr_combo_td !=null) || (fcorr_td !=null) ; // (not the only reason)
		// skipping DISPARITY_VARIATIONS_INDEX - it was not used
		/*
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++) if ((disparity_modes & (1 << i)) != 0){
				if ((i == OVEREXPOSED) && (saturation_imp == null)) {
					continue;
				}
				disparity_map[i] = new double [tilesY*tilesX];
				if ((i >= IMG_TONE_RGB) || ((i >= IMG_DIFF0_INDEX) && (i < (IMG_DIFF0_INDEX + 4)))) {
					need_macro = true;
				}
				if (i <=DISPARITY_STRENGTH_INDEX) {
					need_corr = true;
				}
			}
		}
		*/
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++) {
				if (isSliceBit(i) && ((disparity_modes & (1 << i)) != 0)) { 
					if ((i == OVEREXPOSED) && (saturation_imp == null)) {
						continue;
					}
					disparity_map[i] = new double [tilesY*tilesX];
					if (isCorrBit (i)) {
						need_corr = true;
					}
				} else if (isDiffIndex(i) && needImgDiffs(disparity_modes)){
					disparity_map[i] = new double [tilesY*tilesX];
					need_macro = true;
				} else if (isToneRGBIndex(i) && needTonesRGB(disparity_modes)){
					disparity_map[i] = new double [tilesY*tilesX];
					need_macro = true;
				}
			}
		}

		
		
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}

/*		
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256] - never used
		final double [] lt_window2 = new double [lt_window.length]; // squared - never used

		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];

		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
*/		
		if (globalDebugLevel > 0) {
			System.out.println("macro_mode="+macro_mode);
		}

		final boolean use_main = geometryCorrection_main != null;
		boolean [] used_corrs = new boolean[1];
	    final int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		final TpTask[] tp_tasks =  gpuQuad.setTpTask(
				disparity_array, // final double [][]  disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
				disparity_corr,  // final double       disparity_corr,
				used_corrs,      // final boolean []   need_corrs,       // should be initialized to boolean[1] or null
				tile_op,         // final int [][]     tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				all_pairs,       // final int                      corr_mask,        // <0 - use corr mask from the tile tile_op, >=0 - overwrite all with non-zero corr_mask_tp
				threadsMax);     // final int          threadsMax,       // maximal number of threads to launch
		
		if (tp_tasks.length == 0) {
			System.out.println("Empty tasks - nothing to do");
			return null;
		}
		//texture_tiles
		final boolean fneed_macro = need_macro;
		final boolean fneed_corr =  need_corr && used_corrs[0]; // *** tasks should include correlation

		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > 2);

		final float [] lpf_flat = floatGetCltLpfFd(gpu_sigma_corr);

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				globalDebugLevel > 2);

		final float [] lpf_rb_flat = floatGetCltLpfFd(gpu_sigma_rb_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				globalDebugLevel > 2);

		gpuQuad.setTasks(                  // copy tp_tasks to the GPU memory
				tp_tasks,                  // TpTask [] tile_tasks,
				use_main,                  // use_aux); // boolean use_aux)
				imgdtt_params.gpu_verify); // boolean verify

		gpuQuad.execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers

		TpTask[][] test_tasks = new TpTask[3][];
		test_tasks[2] = tp_tasks;
		if ((fdisp_dist != null) || (fpxpy != null)) {
			final TpTask[] tp_tasks_full = gpuQuad.getTasks(use_main);
			test_tasks[0] = tp_tasks_full;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int indx_tile = ai.getAndIncrement(); indx_tile < tp_tasks_full.length; indx_tile = ai.getAndIncrement()) {
							TpTask task = tp_tasks_full[indx_tile];
							if (fdisp_dist != null) {
								fdisp_dist[task.getTileY()][task.getTileX()] = task.getDispDist();
							}
							if (fpxpy != null) {
								fpxpy[task.getTileY()][task.getTileX()] = task.getXY(); // use_main); // boolean use_aux);
							}
						} // end of tile
					}
				};
			}
			startAndJoin(threads);
			ai.set(0);
		}
		if (fpxpy_test != null) {

			gpuQuad.execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers
			
			final TpTask[] tp_tasks_full = gpuQuad.getTasks(use_main); // reads the same
			test_tasks[1] = tp_tasks_full;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						for (int indx_tile = ai.getAndIncrement(); indx_tile < tp_tasks_full.length; indx_tile = ai.getAndIncrement()) {
							TpTask task = tp_tasks_full[indx_tile];
							if (fpxpy != null) {
								fpxpy_test[task.getTileY()][task.getTileX()] = task.getXY(); // use_main); // boolean use_aux);
							}
						} // end of tile
					}
				};
			}
			startAndJoin(threads);
			ai.set(0);
		}

//			GPUTileProcessor.TpTask[] tp_tasks_full1,
//		GPUTileProcessor.TpTask[] tp_tasks_full2)
		
		
		gpuQuad.execConvertDirect(-1); // boolean erase_clt
		if (iclt_fimg != null) {
			gpuQuad.execImcltRbgAll(isMonochrome());  // execute GPU kernel
			for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
				iclt_fimg[ncam] = gpuQuad.getRBG(ncam); // retrieve data from GPU not used, but now width/height are nominal, not increased 
			}
		} else {gpuQuad.execImcltRbgAll(isMonochrome());} // just for testing
		// does it need texture tiles to be output?
		if (texture_img != null) {
			Rectangle woi = new Rectangle(); // will be filled out to match actual available image
			gpuQuad.execRBGA(
					col_weights,                   // double [] color_weights,
					isLwir(),         // boolean   is_lwir,
					min_shot,       // double    min_shot,           // 10.0
					scale_shot,     // double    scale_shot,         // 3.0
					diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
					diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove);   // boolean   dust_remove,
			float [][] rbga = gpuQuad.getRBGA(
					(isMonochrome() ? 1 : 3), // int     num_colors,
					(texture_woi_pix != null)? texture_woi_pix : woi);
			for (int ncol = 0; ncol < texture_img.length; ncol++) if (ncol < rbga.length) {
				texture_img[ncol] = rbga[ncol];
			}
		}
		// does it need macro data?
		if (fneed_macro) {
			//Generate non-overlapping (16x16) texture tiles, prepare 
			gpuQuad.execTextures(
				col_weights,                   // double [] color_weights,
				isLwir(),         // boolean   is_lwir,
				min_shot,       // double    min_shot,           // 10.0
				scale_shot,     // double    scale_shot,         // 3.0
				diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				diff_threshold, // double    diff_threshold,     // pixel value/pixel change - never used in GPU ?
				min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				dust_remove,    // boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
				false,          // boolean   calc_textures,
				true,           // boolean   calc_extra
				false);         // boolean   linescan_order) // TODO: use true to avoid reordering of the low-res output 
			float [][] extra = gpuQuad.getExtra();
//			int num_cams = gpuQuad.getNumCams();
			int num_cams = getNumSensors();
			for (int ncam = 0; ncam < num_cams; ncam++) {
				int indx = ncam + IMG_DIFF0_INDEX;
//				if ((disparity_modes & (1 << indx)) != 0){
				if (needImgDiffs(disparity_modes)){
					disparity_map[indx] = new double [extra[ncam].length];
					for (int i = 0; i < extra[ncam].length; i++) {
						disparity_map[indx][i] = extra[ncam][i];
					}
				}
			}
			for (int nc = 0; nc < (extra.length - num_cams); nc++) {
				int sindx = nc + num_cams;
				int indx = nc + getImgToneRGB(); // IMG_TONE_RGB;
//				if ((disparity_modes & (1 << indx)) != 0){
				if (needTonesRGB(disparity_modes)){
					disparity_map[indx] = new double [extra[sindx].length];
					for (int i = 0; i < extra[sindx].length; i++) {
						disparity_map[indx][i] = extra[sindx][i];
					}
				}
			}			
		}
		// does it need non-overlapping texture tiles
		if (texture_tiles != null) {   // compatible with the CPU ones
			//Generate non-overlapping (16x16) texture tiles, prepare 
			gpuQuad.execTextures(
				col_weights,                   // double [] color_weights,
				isLwir(),         // boolean   is_lwir,
				min_shot,       // double    min_shot,           // 10.0
				scale_shot,     // double    scale_shot,         // 3.0
				diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				diff_threshold, // double    diff_threshold,     // pixel value/pixel change - never used in GPU ?
				min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				dust_remove,    // boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
				true,           // boolean   calc_textures,
				false,          // boolean   calc_extra
				false);         // boolean   linescan_order) 

			int [] texture_indices = gpuQuad.getTextureIndices();
			int          num_src_slices = numcol + 1; //  + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
			float [] flat_textures =  gpuQuad.getFlatTextures( // fatal error has been detected by the Java Runtime Environment:
					texture_indices.length,
					numcol, // int     num_colors,
		    		false); // clt_parameters.keep_weights); // boolean keep_weights);

			gpuQuad.doubleTextures(
		    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
		    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
		    		texture_indices,                     // int []       indices,
		    		flat_textures,                       // float [][][] ftextures,
		    		tilesX,                              // int          full_width,
		    		4,                                   // rbga only /int          num_slices
		    		num_src_slices                       // int          num_src_slices
		    		);
		}
		
		
		// does it need correlations?
		// FIXME: will not work with combining pairs !!!
		final int num_pairs = Correlation2d.getNumPairs(numSensors);
		if (fneed_corr) {
			//Generate 2D phase correlations from the CLT representation
			int mcorr_sel = Correlation2d.corrSelEncode(imgdtt_params,numSensors);
			int [] i_mcorr_sel = Correlation2d.intCorrPairs(
					mcorr_sel,
					numSensors,
					4); // int num_out); should be 4 int
			gpuQuad.execCorr2D_TD(col_weights,i_mcorr_sel); // Get TD version of correlations (may be read out and saved) 
			final int [] corr_indices = gpuQuad.getCorrIndices();
			if (fcorr_td != null) {
				gpuQuad.getCorrTilesTd(fcorr_td); // generate transform domain correlation pairs
			}
			
			gpuQuad.execCorr2D_normalize(
	        		false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
			
			// Combine 4 ortho pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x0f); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 0) && (fcorr_combo_td[0] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[0]); // generate transform domain correlation pairs for quad ortho combination
			}
			// normalize and convert to pixel domain
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final int [] corr_quad_indices = gpuQuad.getCorrComboIndices(); // get quad
			final float [][] fcorr2D_quad =   gpuQuad.getCorr2DCombo(gpu_corr_rad);

			// Combine 2 diagonal pairs			
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x30); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 1) && (fcorr_combo_td[1] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[1]); // generate transform domain correlation pairs for cross diagonal combination
			}
			gpuQuad.execCorr2D_normalize(
	        		true,          // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero,  // double fat_zero);
					null,          // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_cross =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			
			// Combine 2 horizontal pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x03); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 2) && (fcorr_combo_td[2] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[2]); // generate transform domain correlation pairs for horizontal combination
			}
			gpuQuad.execCorr2D_normalize(
	        		true,          // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero,  // double fat_zero);
					null,          // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_hor =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			// Combine 2 vertical pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        num_pairs,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
					0x0c, // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
					true); // boolean       no_transpose_vertical
			if ((fcorr_combo_td != null) && (fcorr_combo_td.length >= 3) && (fcorr_combo_td[3] != null)) {
				gpuQuad.getCorrTilesComboTd(fcorr_combo_td[3]); // generate transform domain correlation pairs for vertical combination
			}
			gpuQuad.execCorr2D_normalize(
	        		true,          // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero,  // double fat_zero);
					null,          // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_vert =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			
			
			
			if (corr_indices.length > 0) {
				/*				
				if (true) { // debugging only
					int [] wh = new int[2];
					double [][] dbg_corr = GPUTileProcessor.getCorr2DView(
							tilesX,
							tilesY,
							corr_indices,
							fcorr2D,
							wh);
					(new ShowDoubleFloatArrays()).showArrays(
							dbg_corr,
							wh[0],
							wh[1],
							true,
							"dbg-corr2D-initial", // name+"-CORR2D-D"+clt_parameters.disparity,
							GPUTileProcessor.getCorrTitles());
				}
				*/
				final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size
				// assuming that the correlation pairs sets are the same for each tile that has correlations
				// find this number
				int nt0 = (corr_indices[0] >> GPUTileProcessor.CORR_NTILE_SHIFT);
				int nc0 = 1;
				for (int i = 1; (i < corr_indices.length) && ((corr_indices[i] >> GPUTileProcessor.CORR_NTILE_SHIFT) == nt0) ; i++) {
					nc0++;
				}
				final int num_tile_corr = nc0; // normally 6
				final int num_tiles = corr_indices.length / num_tile_corr; 
				
				for (int ithread = 0; ithread < threads.length; ithread++) {
					threads[ithread] = new Thread() {
						@Override
						public void run() {
							//							int tileY,tileX,tIndex;
//							double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
							
							Correlation2d corr2d = new Correlation2d(
									numSensors,
									imgdtt_params,              // ImageDttParameters  imgdtt_params,
									transform_size,             // int transform_size,
									2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
									isMonochrome(), // boolean monochrome,
									(globalDebugLevel > -1));   //   boolean debug)
							corr2d.createOrtoNotch(
									imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
									imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
									(imgdtt_params.lma_debug_level > 1)); // boolean debug);

							for (int indx_tile = ai.getAndIncrement(); indx_tile < num_tiles; indx_tile = ai.getAndIncrement()) {
								// double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
								// added quad and cross combos
								double [][]  corrs = new double [num_pairs + 4][corr_length]; // 225-long (15x15)
								int indx_corr = indx_tile * num_tile_corr;
								int nt = (corr_indices[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
								int tileX = nt % tilesX;
								int tileY = nt / tilesX;
								int tIndex = tileY * tilesX + tileX;
								
								// Prepare the same (currently 10-layer) corrs as double [][], as in CPU version
								int pair_mask = 0;
								for (int indx_pair = 0; indx_pair < num_tile_corr; indx_pair++) {
									int pair = corr_indices[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
									assert pair < num_pairs : "invalid correllation pair";
									pair_mask |= (1 << pair);
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D[indx_corr][i]; // from float to double
									}
									indx_corr++; 
								}
								// add 4 combo layers : quad, cross, hor, vert
								int pair = num_pairs; // 6
								nt = (corr_quad_indices[indx_tile] >> GPUTileProcessor.CORR_NTILE_SHIFT); // corr_quad_indices - different sequence
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_quad[indx_tile][i]; // from float to double
								}
								// indices for cross are the same as for quad
								pair++;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_cross[indx_tile][i]; // from float to double
								}
								
								// indices for hor are the same as for quad
								pair++;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_hor[indx_tile][i]; // from float to double
								}

								// indices for vert are the same as for quad
								pair++;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D_vert[indx_tile][i]; // from float to double
								}
								
								int used_pairs = pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
								boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
								
								corr_common_GPU(
										imgdtt_params,        // final ImageDttParameters  imgdtt_params,
										clt_corr_partial,     // final double [][][][][]   clt_corr_partial,			
										used_pairs,           // final int           used_pairs,
										disparity_map,        // final double [][]   disparity_map,
										clt_mismatch,         // final double [][]   clt_mismatch,
										saturation_imp,       // final boolean [][]  saturation_imp,
										fneed_macro,          // final boolean       fneed_macro,
										corr2d,               // final Correlation2d corr2d,
										corrs,                // final double [][]   corrs,
										tileX,                // final int           tileX,
										tileY,                // final int           tileY,
										max_corr_radius,      // final double        max_corr_radius, // 3.9;
										tile_lma_debug_level, // int                 tile_lma_debug_level,
										debugTile,            // boolean             debugTile,
										globalDebugLevel);    // final int           globalDebugLevel)							
								
								// double extra_disparity = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
/*
// Disabled for GPU
								if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
								else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
								else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
								else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];  // not used in lwir
								else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];  // not used in lwir
								if (Double.isNaN(extra_disparity)) extra_disparity = 0;  // used in lwir
*/
								if ((disparity_map != null) && Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
									System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
								}

								// only debug is left
								// old (per-color correlation)
								// removed
								
							} // end of tile
						}
					};
				}
				startAndJoin(threads);
			} else {
				// no correlation tiles to process
			}
		}
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
		return test_tasks;
		
//			GPUTileProcessor.TpTask[] tp_tasks_full1,
//		GPUTileProcessor.TpTask[] tp_tasks_full2)
		
	}
	
	@Deprecated
	public void clt_process_tl_correlations_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			// both arrays should have same non-null tiles
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][4*64] TD of combo corrs: qud, cross, hor,vert
			                                           // each of the top elements may be null to skip particular combo type
			final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
            // 											  [tilesY][tilesX] should be set by caller
			final float  [][][]       fcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
			final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
			final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
			final int                 gpu_corr_rad,    // = transform_size - 1 ?
			final double              max_corr_radius, // 3.9;
			final int                 window_type,     // GPU: will not be used
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		clt_process_tl_correlations_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
				imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				fcorr_combo_td,  // [4][tilesY][tilesX][4*64] TD of combo corrs: qud, cross, hor,vert
				                                           // each of the top elements may be null to skip particular combo type
				corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
				clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
	            // 											  [tilesY][tilesX] should be set by caller
				fcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				// When clt_mismatch is non-zero, no far objects extraction will be attempted
				clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                           // last 2 - contrast, avg/ "geometric average)
				disparity_modes, // bit mask of disparity_map slices to calculate/return
				imgdtt_params.dbg_pair_mask, // //corr_pairs_mask, // which pairs to use in LMA
				gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
				gpu_corr_rad,    // = transform_size - 1 ?
				max_corr_radius, // 3.9;
				window_type,     // GPU: will not be used
				debug_tileX,
				debug_tileY,
				threadsMax,      // maximal number of threads to launch
				globalDebugLevel);
	}
	
	@Deprecated
	public void clt_process_tl_correlations_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			// both arrays should have same non-null tiles
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][4*64] TD of combo corrs: qud, cross, hor,vert
			                                           // each of the top elements may be null to skip particular combo type
			final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
            // 											  [tilesY][tilesX] should be set by caller
			final float  [][][]       fcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
			final int                 corr_pairs_mask, // which pairs to use in LMA
			final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
			final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
			final int                 gpu_corr_rad,    // = transform_size - 1 ?
			final double              max_corr_radius, // 3.9;
			final int                 window_type,     // GPU: will not be used
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int [] dbg_tiles = (globalDebugLevel > -1)?
				new int []{37721, 39022, 40975, 38042, 40963, 42253, 38694, 39343, 36443} : null;
		final float gpu_fcorr_scale = (float) gpu_corr_scale;
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return;
		}
		final boolean [][] saturation_imp = gpuQuad.quadCLT.saturation_imp;               // boolean [][] saturation_imp, // (near) saturated pixels or null
//gpuQuad
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;

//		final double [][] debug_offsets = new double[imgdtt_params.lma_dbg_offset.length][2];
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}
//dbg_pair_mask
		final int quad = 4;   // number of subcameras
//		final int numcol = isMonochrome()?1:3;

		final int width =  gpuQuad.getImageWidth();
		final int height = gpuQuad.getImageHeight();
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*quad][tilesX*tilesY]) : null;
		// not yet used with GPU
		// keep for now for mono, find out  what do they mean for macro mode
		
		final int corr_size = transform_size * 2 - 1;
		

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 1){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}


		if (globalDebugLevel > 1) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++) if ((disparity_modes & (1 << i)) != 0){
				if ((i == OVEREXPOSED) && (saturation_imp == null)) {
					continue;
				}
				disparity_map[i] = new double [tilesY*tilesX];
			}
		}
		
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}

		/*
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256] - never used
		final double [] lt_window2 = new double [lt_window.length]; // squared - never used
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];
		*/
///		final boolean use_main = geometryCorrection_main != null;
///	    final int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		
		int [] corr_indices_ = null; 
		float [][] fcorr2D_ = null;
		
		if (fcorr_td != null) {
			
///			int [][] pairs_map = {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5}};
///			corr_indices_ = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
///					fcorr_td, // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
///					pairs_map); // int [][] pairs) // typically {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5} [0] - 3rd index in corr_tiles, [1] -
			
			corr_indices_ = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
					fcorr_td, // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
		            null, // int [][][]           num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] (or null)
		            null); // float [][]           pfcorr_weights); // null or one per correlation tile (num_corr_tiles) to divide fat zero2
			
			gpuQuad.execCorr2D_normalize(
					false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
					gpu_corr_rad); // int corr_radius
			fcorr2D_ = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
		}
		
		final int [] corr_indices = corr_indices_; 
		final float [][] fcorr2D = fcorr2D_; // len = 49311

		
		
		
		
		
		final int num_combo= (fcorr_combo_td != null)? fcorr_combo_td.length : 0;
		final int [][] corr_combo_indices = (fcorr_combo_td != null)? new int [num_combo][] : null;
		final float [][][] fcorr2D_combo = (fcorr_combo_td != null)? new float[num_combo][][] : null;
		if (num_combo > 0) {
			int [] ipairs = {0xf, 0x30, 0x3, 0xc}; // does not matter. using quad, cross, hor, vert
			for (int i = 0; i < num_combo; i++) {
				corr_combo_indices[i] = gpuQuad.setCorrTilesComboTd(
						fcorr_combo_td[i], // final float [][][] corr_tiles, // [tileY][tileX][4*64]
						ipairs[i]); // change to smth? int ipair) // just to set in the index low bits
				// normalize and convert to pixel domain
				gpuQuad.execCorr2D_normalize(
						true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
						gpu_fat_zero, // double fat_zero);
						null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
						gpu_corr_rad); // int corr_radius
				fcorr2D_combo[i] =  gpuQuad.getCorr2DCombo(gpu_corr_rad);
			}
		}
		final boolean fneed_macro = false;
		
		if (corr_indices.length > 0) {
			/*
			if (true) { // debugging only
				int [] wh = new int[2];
				double [][] dbg_corr = GPUTileProcessor.getCorr2DView(
						tilesX,
						tilesY,
						corr_indices,
						fcorr2D,
						wh);
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_corr,
						wh[0],
						wh[1],
						true,
						"dbg-corr2D-fromTD", // name+"-CORR2D-D"+clt_parameters.disparity,
						GPUTileProcessor.getCorrTitles());
			}
			 */
			final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size
			// assuming that the correlation pairs sets are the same for each tile that has correlations
			// find this number
			int nt0 = (corr_indices[0] >> GPUTileProcessor.CORR_NTILE_SHIFT);
			int nc0 = 1;
			for (int i = 1; (i < corr_indices.length) && ((corr_indices[i] >> GPUTileProcessor.CORR_NTILE_SHIFT) == nt0) ; i++) {
				nc0++;
			}
			final int num_tile_corr = nc0; // normally 6
			final int num_tiles = corr_indices.length / num_tile_corr;
			// FIXME: will not work with combining pairs !!!
			final int num_pairs = Correlation2d.getNumPairs(numSensors);


			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {

						Correlation2d corr2d = new Correlation2d(
								numSensors,
								imgdtt_params,              // ImageDttParameters  imgdtt_params,
								transform_size,             // int transform_size,
								2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
								isMonochrome(), // boolean monochrome,
								(globalDebugLevel > -1));   //   boolean debug)
						corr2d.createOrtoNotch(
								imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
								imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
								(imgdtt_params.lma_debug_level > 1)); // boolean debug);

						for (int indx_tile = ai.getAndIncrement(); indx_tile < num_tiles; indx_tile = ai.getAndIncrement()) {
							// double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
							// added quad and cross combos
							double [][]  corrs = new double [num_pairs + num_combo][corr_length]; // 225-long (15x15)
							float  [][] fcorrs = (fcorr_tiles == null) ? null : new float  [num_pairs + num_combo][corr_length]; // 225-long (15x15)
							int indx_corr = indx_tile * num_tile_corr;
							int nt = (corr_indices[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
							int tileX = nt % tilesX;
							int tileY = nt / tilesX;
							int tIndex = tileY * tilesX + tileX;

							// Prepare the same (currently 10-layer) corrs as double [][], as in CPU version
							int pair_mask = 0;
							if (fcorr_td != null) {
								for (int indx_pair = 0; indx_pair < num_tile_corr; indx_pair++) {
									int pair = corr_indices[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
									assert pair < num_pairs : "invalid correllation pair";
									pair_mask |= (1 << pair);
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D[indx_corr][i]; // from float to double
									}
									if (fcorrs != null) for (int i = 0; i < corr_length; i++) {
										fcorrs[pair][i] = gpu_fcorr_scale * fcorr2D[indx_corr][i];
									}
									indx_corr++; 
								}
							}
							// add 4 combo layers : quad, cross, hor, vert
							if (num_combo > 0) {
								for (int ncm = 0; ncm < num_combo; ncm++) if (corr_combo_indices[ncm]!=null){
									nt = (corr_combo_indices[ncm][indx_tile] >> GPUTileProcessor.CORR_NTILE_SHIFT); // corr_quad_indices - different sequence
									int pair = num_pairs + ncm; // 6+
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D_combo[ncm][indx_tile][i]; // from float to double
									}
									if (fcorrs != null) for (int i = 0; i < corr_length; i++) {
										fcorrs[pair][i] = gpu_fcorr_scale * fcorr2D_combo[ncm][indx_tile][i];
									}
								}
							}
							if (corr_tiles != null) {
								corr_tiles[tileY][tileX] = corrs; 
							}
							if (fcorr_tiles != null) {
								fcorr_tiles[tileY * tilesX + tileX] = fcorrs; // does not require corr_common_GPU()
							}
							if ((disparity_map != null) || (clt_corr_partial != null) || (clt_mismatch != null)) {
//								int used_pairs = pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int used_pairs = pair_mask & corr_pairs_mask; // imgdtt_params.dbg_pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
								boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
								if (dbg_tiles != null) for (int ii : dbg_tiles){
									if (ii == tIndex) {
										debugTile=true;
										System.out.println("Debugging tile "+tIndex+" ( "+tileX+" / "+tileY+" )");
										break;
									}
								}
								corr_common_GPU(
										imgdtt_params,        // final ImageDttParameters  imgdtt_params,
										clt_corr_partial,     // final double [][][][][]   clt_corr_partial,
										used_pairs,           // final int           used_pairs,
										disparity_map,        // final double [][]   disparity_map,
										clt_mismatch,         // final double [][]   clt_mismatch,
										saturation_imp,       // final boolean [][]  saturation_imp,
										fneed_macro,          // final boolean       fneed_macro,
										corr2d,               // final Correlation2d corr2d,
										corrs,                // final double [][]   corrs,
										tileX,                // final int           tileX,
										tileY,                // final int           tileY,
										max_corr_radius,      // final double        max_corr_radius, // 3.9;
										tile_lma_debug_level, // int                 tile_lma_debug_level,
										debugTile,            // boolean             debugTile,
										globalDebugLevel);    // final int           globalDebugLevel)							
							}
							// double extra_disparity = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
							/*
// Disabled for GPU
								if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
								else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
								else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
								else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];  // not used in lwir
								else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];  // not used in lwir
								if (Double.isNaN(extra_disparity)) extra_disparity = 0;  // used in lwir
							 */
							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}

							// only debug is left
							// old (per-color correlation)
							// removed

						} // end of tile
					}
				};
			}
			startAndJoin(threads);
		} else {
			// no correlation tiles to process
		}
		if ((dbg_distort != null) &&(globalDebugLevel >= 0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
/*
//		final double [][] dbg_distort = debug_distort? (new double [4*quad][tilesX*tilesY]) : null;
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}

		if (dbg_ports_coords != null) {
			(new showDoubleFloatArrays()).showArrays(dbg_ports_coords,  tilesX, tilesY, true, "ports_coordinates", dbg_titles);
		}
*/
//		return clt_data;
	}
	
	public double [][][] clt_process_tl_interscene( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
	        final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
	        float [][][]              num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null). Can be inner null if not used in tp_tasks
	        double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
	        final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
	        final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
	        final int                 gpu_corr_rad,    // = transform_size - 1 ?
	        // The tp_tasks data should be decoded from GPU to get coordinates
			final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMA for the integrated correlations
			// to be converted to float (may be null)
			final double  [][][]      dcorr_tiles,     // [tile][pair absolute, sparse][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			final double [][]         pXpYD,           // pXpYD for the reference scene 
			final double [][]         fpn_offsets,     // null, or per-tile X,Y offset to be blanked
			final double              fpn_radius,      // radius to be blanked around FPN offset center
			final boolean             fpn_ignore_border, // only if fpn_mask != null - ignore tile if maximum touches fpn_mask			
			final double[][][]        motion_vectors,  // [tilesY*tilesX][][] -> [][num_sel_sensors+1 or 2][3]
			final boolean             run_poly,        // polynomial max, if false - centroid
			final boolean             use_partial,     // find motion vectors for individual pairs, false - for sum only
			final double              centroid_radius, // 0 - use all tile, >0 - cosine window around local max
			final int                 n_recenter,      // when cosine window, re-center window this many times
			final double              td_weight,       // mix correlations accumulated in TD with 
			final double              pd_weight,       // correlations (post) accumulated in PD
			final boolean             td_nopd_only,    // only use TD accumulated data if no safe PD is available for the tile.
			final double              min_str_nofpn,    //  = 0.25;
			final double              min_str_sum_nofpn,// = 0.8; // 5;
			final double              min_str_fpn,      //  = 0.25;
			final double              min_str_sum_fpn,  // = 0.8; // 5;
			final int                 min_neibs,       //   2;	   // minimal number of strong neighbors (> min_str)
			final double              weight_zero_neibs,//  0.2;   // Reduce weight for no-neib (1.0 for all 8)
			final double              half_disparity,  //   5.0;   // Reduce weight twice for this disparity
			final double              half_avg_diff,   //   0.2;   // when L2 of x,y difference from average of neibs - reduce twice
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
				
		if (this.gpuQuad == null) {
			System.out.println("clt_process_tl_interscene(): this.gpuQuad is null, bailing out");
			return null;
		}
//		final int min_neibs = clt_parameters.imp.min_neibs;
		final boolean extra_sum = true; // use sum of pixel-domain correlations (TD have artifacts for low contrast
		// - maybe  -related to float vs. double - not tested yet 
//		final int width =  gpuQuad.getImageWidth();
//		final int height = gpuQuad.getImageHeight();
		final int tilesX=  gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=  gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final double [][][] coord_motion = new double [(pXpYD != null)?2:1][tilesX * tilesY][];
		// not yet used with GPU
		// keep for now for mono, find out  what do they mean for macro mode
		final int corr_size = transform_size * 2 - 1;
		final double [][]         debug_lma = imgdtt_params.lmamask_dbg? (new double [6][tilesX*tilesY]):null;
		if (debug_lma != null) {
			for (int i = 0; i < debug_lma.length; i++) {
				Arrays.fill(debug_lma[i], Double.NaN);
			}
		}
		final float [][]           pfcorr_weights = ((num_acc != null) || (dcorr_weight != null))? new float[1][] : null;
		// This version obeys tp_task order and fills fcorr_td gaps (should be none_) with zeros.
		int [] corr_indices ;
		if (fcorr_td == null) { // used with no accumulation, assume TD correlation data is still in GPU
			corr_indices = gpuQuad.getCorrIndices(); // also sets num_corr_tiles
		} else {
			if (num_acc != null) { // version with float [][][]           num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] (or null)
				corr_indices = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
						tp_tasks,        // final TpTask []      tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
						true,            // final boolean        inter_mode, // interscene mode
						fcorr_td,        // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
						num_acc,         // float [][][]           num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] (or null)
						pfcorr_weights); // float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2
			} else { // version with // double []            dcorr_weight, // [ntile] (or null), compatible with the CPU version 
				corr_indices = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
						tp_tasks,        // final TpTask []      tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
						true, // final boolean        inter_mode, // interscene mode
						fcorr_td,        // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
						dcorr_weight,    // double []            dcorr_weight, // [ntile] (or null)
						pfcorr_weights); // float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2
			}
		}
		int dbg_imax = 0;
		for (int ii = 1; ii < corr_indices.length; ii++) {
			if (corr_indices[ii] > corr_indices[dbg_imax]) {
				dbg_imax=ii;
			}
		}
//		System.out.println("dbg_imax = "+dbg_imax+", corr_indices[dbg_imax]="+corr_indices[dbg_imax]+" tile="+((dbg_imax-16)/17));
		if (corr_indices.length == 0) {
			return null;
		}
		float [] fcorr_weights = ((num_acc != null) || (dcorr_weight != null))? pfcorr_weights[0] : null;
		gpuQuad.execCorr2D_normalize(
				false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
				gpu_fat_zero,            // double fat_zero);
				fcorr_weights,           // fcorr_weights,           // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
				gpu_corr_rad);           // int corr_radius
		
		final float [][] fcorr2D = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);

		final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size
		
		final int [] used_sensors_list = gpuQuad.getSensInter();
		
		final int extra_len = extra_sum? 1 : 0;
		final int corrs_len = (use_partial?used_sensors_list.length:1); // without optional extra_len but including GPU sum

		
		final int num_tiles = corr_indices.length / used_sensors_list.length; // number of correlated tiles (not in tp_tasks) 

// Add (and init by caller) if needed, so far static is enough
//		if (correlation2d == null) {
//			throw new IllegalArgumentException ("clt_process_tl_correlations(): correlation2d == null!");
//		}

		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}

		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX,nTile; // , chn;
					for (int iCorrTile = ai.getAndIncrement(); iCorrTile < num_tiles; iCorrTile = ai.getAndIncrement()) {
						nTile = (corr_indices[iCorrTile* used_sensors_list.length] >> GPUTileProcessor.CORR_NTILE_SHIFT);
						tileY = nTile / tilesX;
						tileX = nTile % tilesX;
						boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > 2); // 0);
						if (debugTile0) {
							System.out.println("clt_process_tl_correlations(): tileX="+tileX+", tileY="+tileY+", nTile="+nTile+", nTile="+nTile);
						}
						// zero FPN right in the fcorr2D so it will go to both processing and display
						boolean [] fpn_mask = null;
						double min_str = min_str_nofpn;         // higher threshold when FPN is possible
						double min_str_sum = min_str_sum_nofpn; // higher threshold when FPN is possible
						if ((fpn_offsets != null) && (fpn_offsets[nTile] != null)) {
							double fpn_x = transform_size - 1 - fpn_offsets[nTile][0]; //  0 -> 7.0 
							double fpn_y = transform_size - 1 - fpn_offsets[nTile][1]; //  0 -> 7.0
							int min_x = (int) Math.max(Math.round(fpn_x - fpn_radius),0); 
							int max_x = (int) Math.min(Math.round(fpn_x + fpn_radius), corr_size-1); 
							int min_y = (int) Math.max(Math.round(fpn_y - fpn_radius),0); 
							int max_y = (int) Math.min(Math.round(fpn_y + fpn_radius), corr_size-1); 
							int fcorr2D_indx = (iCorrTile + 1)*  used_sensors_list.length -1; // last in each group - sum in TD
							fpn_mask = new boolean[fcorr2D[fcorr2D_indx].length];
							for (int iy = min_y; iy <= max_y; iy++) {
								for (int ix = min_x; ix <= max_x; ix++) {
									int indx = iy * corr_size + ix;
									fcorr2D[fcorr2D_indx][indx] = 0;
									fpn_mask[indx] = true;
								}
							}
							min_str = min_str_fpn;
							min_str_sum = min_str_sum_fpn;
						}
						double [][] corrs = new double [corrs_len + extra_len][];
						// copy correlation tiles from the GPU's floating point arrays
						double scale = 1.0/getNumSensors();
						if (extra_sum) {
							corrs[corrs_len] = new double [corr_length];
						}
						for (int isens = corrs_len - 1; isens >= 0; isens--) {
							int nsens = used_sensors_list.length - corrs_len + isens;
							corrs[isens] = new double[corr_length];
							int fcorr2D_indx = iCorrTile *  used_sensors_list.length + nsens;
							for (int i = 0; i < corr_length; i++) {
								corrs[isens][i] = gpu_corr_scale * fcorr2D[fcorr2D_indx][i]; // copy one-by-one converting from floats to doubles
							}
							if (use_partial && (isens < (corrs_len - 1))) { // not including sum
								for (int i = 0; i < corr_length; i++) {
									corrs[corrs_len][i] += scale*corrs[isens][i];
								}
							}
						}
						if (!use_partial && extra_sum) {
							scale *= gpu_corr_scale;
							for (int nsens = 0; nsens < (used_sensors_list.length - 1); nsens++) {
								int fcorr2D_indx = iCorrTile *  used_sensors_list.length + nsens;
								for (int i = 0; i < corr_length; i++) {
									corrs[corrs_len][i] += scale * fcorr2D[fcorr2D_indx][i]; // copy one-by-one converting from floats to doubles
								}
							}
						}

						if (dcorr_tiles != null) { // This will be visualized
							int index_es = getNumSensors() + extra_len;
							dcorr_tiles[iCorrTile] = new double[getNumSensors()+1 + extra_len][];
							if (extra_sum) {
								dcorr_tiles[iCorrTile][index_es] = new double[corr_length];
							}
							for (int nsens = 0; nsens < used_sensors_list.length; nsens++) {
								int abs_sens = used_sensors_list[nsens];
								if (abs_sens >= getNumSensors()) {
									abs_sens = getNumSensors(); // last - sum of all sensors
								}
								if ((abs_sens < getNumSensors()) && extra_sum) {
									int fcorr2D_indx = iCorrTile *  used_sensors_list.length + nsens;
									for (int i = 0; i < corr_length; i++) {
										dcorr_tiles[iCorrTile][index_es][i] += scale * gpu_corr_scale * fcorr2D[fcorr2D_indx][i]; // copy one-by-one converting from floats to doubles 	
									}
								}
								dcorr_tiles[iCorrTile][abs_sens] = new double[corr_length];
								int fcorr2D_indx = iCorrTile *  used_sensors_list.length + nsens;
								for (int i = 0; i < corr_length; i++) {
									dcorr_tiles[iCorrTile][abs_sens][i] = gpu_corr_scale * fcorr2D[fcorr2D_indx][i]; // copy one-by-one converting from floats to doubles 	
								}
							}
						}
						//			final double [][][][]     motion_vectors,  // [tilesY][tilesX][][] -> [][][num_sel_sensors+1][2]
						if (motion_vectors != null) { // TODO: now used only as debug, may be removed later
							motion_vectors[nTile] = new double [corrs.length][];
							for (int nsens = 0; nsens < corrs.length; nsens++) { // all 18
								double min_vstr =  (nsens == (corrs_len - 1))? min_str_sum : min_str;
								motion_vectors[nTile][nsens] = Correlation2d.getMaxXYCm(
										corrs[nsens],    // double [] data,
										corr_size,       // int       data_width,      //  = 2 * transform_size - 1;
										centroid_radius, // double    radius, // 0 - all same weight, > 0 cosine(PI/2*sqrt(dx^2+dy^2)/rad)
										n_recenter,      // int       refine, //  re-center window around new maximum. 0 -no refines (single-pass)
										null,            // boolean [] fpn_mask,
										false,           // boolean    ignore_border, // only if fpn_mask != null - ignore tile if maximum touches fpn_mask
										false);          // boolean   debug)
								if (motion_vectors[nTile][nsens] != null) {
									if (motion_vectors[nTile][nsens][2] < min_vstr) {
										motion_vectors[nTile][nsens] = null;
									}
								}
							}
						}
						// now calculate only for composite
						double [] mv_pd = new double [3];
						double [] mv_td = new double [3];
						if (pd_weight > 0.0) {
							mv_pd = Correlation2d.getMaxXYCm(
								corrs[corrs.length-1], // double [] data,
								corr_size,             // int       data_width,      //  = 2 * transform_size - 1;
								centroid_radius,       // double    radius, // 0 - all same weight, > 0 cosine(PI/2*sqrt(dx^2+dy^2)/rad)
								n_recenter,            // int       refine, //  re-center window around new maximum. 0 -no refines (single-pass)
								null,                  // boolean [] fpn_mask,
								false,                 // boolean    ignore_border, // only if fpn_mask != null - ignore tile if maximum touches fpn_mask
								false);                // boolean   debug)
							if (mv_pd != null) {
								if (mv_pd[2] < min_str) {
									mv_pd = null;
								} else {
									mv_pd[2] -= min_str;
								}
							}
						}
						if (td_weight > 0.0) {
							mv_td = Correlation2d.getMaxXYCm(
								corrs[corrs.length-2], // double [] data,
								corr_size,             // int       data_width,      //  = 2 * transform_size - 1;
								centroid_radius,       // double    radius, // 0 - all same weight, > 0 cosine(PI/2*sqrt(dx^2+dy^2)/rad)
								n_recenter,            // int       refine, //  re-center window around new maximum. 0 -no refines (single-pass)
								fpn_mask,              // boolean [] fpn_mask,
								false,                 // boolean    ignore_border, // only if fpn_mask != null - ignore tile if maximum touches fpn_mask
								false);                // boolean   debug)
							if (mv_td != null) {
								if (mv_td[2] < min_str_sum) {
									mv_td = null;
								} else {
									mv_td[2] -= min_str_sum;
								}
							}
						}
						if ((mv_td != null) || (mv_pd != null)) {
							double [] mv = new double[3];
							if (mv_pd != null) {
								mv = mv_pd;
								mv[2] *= pd_weight;
								if ((mv_td != null) && !td_nopd_only) {  // mix
									mv[0] = (mv[0] * pd_weight + mv_td[0] * td_weight)/ (pd_weight + td_weight);
									mv[1] = (mv[1] * pd_weight + mv_td[1] * td_weight)/ (pd_weight + td_weight);
									mv[2] += mv_td[2] * td_weight;
								} // mix
							} else { // (mv_pd == null) &&  (mv_td != null)  below
								mv = mv_td;
								mv[2] *= td_weight;
							}
							if (mv != null) {
								if (pXpYD == null) {
									coord_motion[0][nTile] = mv;
								} else {
									if (pXpYD[nTile] != null) { // seems always
										coord_motion[0][nTile] = pXpYD[nTile].clone();
										coord_motion[1][nTile] = mv;
									}
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		
		ai.set(0);
		final int tiles = tilesX * tilesY;
		final double [][] mv = coord_motion[coord_motion.length - 1];
		final double [][] pxd = (coord_motion.length>1) ? coord_motion[0] : null;
		final double scale_num_neib = ((weight_zero_neibs >= 0) && (weight_zero_neibs < 1.0)) ? (weight_zero_neibs * 8/(1.0 - weight_zero_neibs)): 0.0;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					double l2;
					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
//						if (nTile==162) {
//							System.out.println("nTile="+nTile);
//						}
						if ((mv[nTile] != null) && (pXpYD[nTile] != null)) { 
							int num_neibs=0;
							double sx=0.0, sy=0.0;
							for (int dir = 0; dir < 8; dir++) {
								int nTile1 = tn.getNeibIndex(nTile, dir);
								if ((nTile1 >= 0) &&
										(mv[nTile1] != null) &&
										(pXpYD[nTile1] != null) &&
										!Double.isNaN(mv[nTile1][2]) &&
										!Double.isNaN(mv[nTile1][0]) &&
										!Double.isNaN(mv[nTile1][1])){
									num_neibs++;
									sx += mv[nTile1][0];
									sy += mv[nTile1][1];
								}
							}
							if (num_neibs < min_neibs) { // filter by minimal neighbors
								mv[nTile][2] = 0;
								/*
								mv[nTile] = null;
								if (pxd != null) {
									pxd[nTile] = null;
								}
								*/
								continue;
							}
							if ((weight_zero_neibs >= 0) && (weight_zero_neibs < 1.0)) { // scale weight by number of neighbors
								mv[nTile][2] *= (num_neibs + scale_num_neib)/(8.0 + scale_num_neib);
							}
							if (half_disparity > 0.0) { // scale by disparity
								mv[nTile][2] *= half_disparity/(half_disparity + pxd[nTile][2]);
							}
							if ((half_avg_diff > 0.0) &&(num_neibs > 0)) {
								double dx = mv[nTile][0] - sx / num_neibs;
								double dy = mv[nTile][1] - sy / num_neibs;
								l2 = Math.sqrt(dx * dx + dy * dy);
								mv[nTile][2] *= half_avg_diff/(half_avg_diff + l2);
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);

		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
//					double l2;
//					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
//						if (nTile==162) {
//							System.out.println("nTile="+nTile);
//						}
						if ((mv[nTile] != null) && (pXpYD[nTile] != null)) {
							if (mv[nTile][2] <= 0) {
								mv[nTile] = null;
								if (pxd != null) {
									pxd[nTile] = null;
								}
							}
						}
					}
				}
			};
		}
		startAndJoin(threads);
		return coord_motion;
	}

	
	
	// using most of the ImageDttCPU.clt_process_tl_correlations
	public void clt_process_tl_correlations( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
	        final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of all selected corr pairs
	        float [][][]              num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null). Can be inner null if not used in tp_tasks
	        double []                 dcorr_weight,    // alternative to num_acc, compatible with CPU processing (only one non-zero enough)
	        final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
	        final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
	        final int                 gpu_corr_rad,    // = transform_size - 1 ?
	        // The tp_tasks data should be decoded from GPU to get coordinates
			final TpTask []           tp_tasks,        // data from the reference frame - will be applied to LMA for the integrated correlations
			final double [][]         rXY,             // from geometryCorrection
			// next both can be nulls
		    final double [][][][]     clt_corr_out,   // sparse (by the first index) [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] or null
		    // combo will be added as extra pair if mcorr_comb_width > 0 and clt_corr_out has a slot for it
			// to be converted to float
			final double  [][][]      dcorr_tiles,     // [tile][pair][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			//optional, may be null
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			final double [][][][]     ddnd,            // [tilesY][tilesX][num_sensors][2] data for LY. Should be either null or [tilesY][tilesX][][]. disparity_map should be non-null
			final boolean             run_lma,         // calculate LMA, false - CM only
  		    // define combining of all 2D correlation pairs for CM (LMA does not use them)
			final int                 mcorr_comb_width,  // combined correlation tile width (set <=0 to skip combined correlations)
			final int                 mcorr_comb_height, // combined correlation tile full height
			final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
			final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
			final int                 window_type,     // GPU: will not be used
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final double disparity_scale = 1.0/Math.sqrt(2); // combo pixels -> disparity pixels
		final boolean diameters_combo = (imgdtt_params.mcorr_dual_fract > 0.0); // add diameters-only combo after all-combo
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return;
		}

		
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}

		final int width =  gpuQuad.getImageWidth();
		final int height = gpuQuad.getImageHeight();
		final int tilesX=  gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=  gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		// not yet used with GPU
		// keep for now for mono, find out  what do they mean for macro mode
		
		final int corr_size = transform_size * 2 - 1;
		//FIXME: lmaDisparityStrengths expects length = 14;
		final String[] debug_lma_titles_nobi = {"disp_samples","num_cnvx_samples","num_comb_samples",
				"num_lmas","num_iters","rms"};
		final String[] debug_lma_titles_bi  = {"disparity","strength_mod","strength", "area","ac","min(a,c)",
				"max(a,c)","a","c","b","str1","rrms","rms","fail_reason"};
		final String[] debug_lma_titles= imgdtt_params.bimax_dual_LMA? debug_lma_titles_bi:debug_lma_titles_nobi;
//		final double [][]         debug_lma = imgdtt_params.lmamask_dbg? (new double [6][tilesX*tilesY]):null;
		final double [][]         debug_lma = imgdtt_params.lmamask_dbg? (new double [debug_lma_titles.length][tilesX*tilesY]):null;

		if (debug_lma != null) {
			for (int i = 0; i < debug_lma.length; i++) {
				Arrays.fill(debug_lma[i], Double.NaN);
			}
		}
		

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 1){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}


		if (globalDebugLevel > 1) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		
		final float [][]           pfcorr_weights = ((num_acc != null) || (dcorr_weight != null))? new float[1][] : null;
// This version obeys tp_task order and fills fcorr_td gaps (should be none_) with zeros.
		
		int [] corr_indices;
		if (num_acc != null) { // version with float [][][]           num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] (or null)
			corr_indices = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
					tp_tasks,        // final TpTask []      tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					false, // final boolean        inter_mode, // interscene mode
					fcorr_td,        // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
					num_acc,         // float [][][]           num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] (or null)
					pfcorr_weights); // float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		} else { // version with // double []            dcorr_weight, // [ntile] (or null), compatible with the CPU version 
			corr_indices = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
					tp_tasks,        // final TpTask []      tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
					false, // final boolean        inter_mode, // interscene mode
					fcorr_td,        // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
					dcorr_weight,    // double []            dcorr_weight, // [ntile] (or null)
					pfcorr_weights); // float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		}
		float [] fcorr_weights = ((num_acc != null) || (dcorr_weight != null))? pfcorr_weights[0] : null;
		gpuQuad.execCorr2D_normalize(
				false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
				gpu_fat_zero,            // double fat_zero);
				fcorr_weights,           // fcorr_weights,           // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
				gpu_corr_rad);           // int corr_radius
		final float [][] fcorr2D = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
		
		if (corr_indices.length > 0) {
			final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size

			// CPU version below
			// CPU version uses sparse 2D correlation arrays (length 120 for 16 cameras, unused are null. GPU version compacts these arrays
			int [] used_pairs_list = gpuQuad.getUsedPairsList(); // GPU index -> CPU index


			if (correlation2d == null) {
				throw new IllegalArgumentException ("clt_process_tl_correlations(): correlation2d == null!");
			}
			correlation2d.setCorrPairs(gpuQuad.getCorrMask()); // copy correlation selection from GPU to Correlation2d instance 
			final int num_pairs = correlation2d.getCorrPairs().length; // will throw if null

			// which to use - num_pairs or num_used_pairs? Or set correlation2d to match num_used_pairs
			
			final boolean combine_corrs = (mcorr_comb_width > 0); // use 0 to disable combining (use LMA only)
			if (!combine_corrs) {
				System.out.println("**** Warning: clt_process_tl_correlations() does not use combine_corrs, new LMA is disbled ****");
			}

			if (clt_corr_out != null) {  // not used so far? REMOVE?
				boolean [] calc_corr_pairs = correlation2d.getCorrPairs();
				for (int i = 0; i < calc_corr_pairs.length; i++) if (calc_corr_pairs[i]){
					clt_corr_out[i] = new double[tilesY][tilesX][]; 
				}
				if (clt_corr_out.length > num_pairs) {
					clt_corr_out[num_pairs] = new double[tilesY][tilesX][]; // for combo
				}
			}

			if (combine_corrs) {
				correlation2d.generateResample( // should be called before *** This can be done in CPU, table(s) copied to GPU 
						mcorr_comb_width,  // combined correlation tile width
						mcorr_comb_height, // combined correlation tile full height
						mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
						mcorr_comb_disp);
			}
			final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
					transform_size,
					imgdtt_params.lma_wnd);
			final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
			if (corr_wnd_inv_limited != null) {
				double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
				for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
					for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
						corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
					}
				}
			}

			if (disparity_map != null){ // only slices that are needed, keep others
				boolean use_bimax = imgdtt_params.bimax_dual_LMA;
				ArrayList<Integer> used_slices = new ArrayList<Integer>();
				for (int i : new int[] {DISPARITY_INDEX_CM, DISPARITY_INDEX_CM+1,DISPARITY_STRENGTH_INDEX})	used_slices.add(i);
				if (run_lma) for (int i : new int[] {DISPARITY_INDEX_POLY,DISPARITY_INDEX_POLY+1})	        used_slices.add(i);
				if (use_bimax) used_slices.add(DISPARITY_VARIATIONS_INDEX);
				ArrayList<Integer> nan_slices = new ArrayList<Integer>();
				nan_slices.add(DISPARITY_INDEX_CM);
				if (run_lma) nan_slices.add(DISPARITY_INDEX_POLY);
				for (int indx:used_slices) {
					disparity_map[indx] = new double[tilesY*tilesX];
				}
				for (int indx:nan_slices) {
					Arrays.fill(disparity_map[indx], Double.NaN);
				}
			}

			final Thread[] threads = newThreadArray(threadsMax);
			final AtomicInteger ai = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						DttRad2 dtt = new DttRad2(transform_size);
						dtt.set_window(window_type);
						int tileY,tileX,nTile; // , chn;
						for (int iTile = ai.getAndIncrement(); iTile < tp_tasks.length; iTile = ai.getAndIncrement()) {
							tileY = tp_tasks[iTile].getTileY(); //  /tilesX;
							tileX = tp_tasks[iTile].getTileX(); //nTile % tilesX;
							if (fcorr_td[tileY][tileX] == null) {
								continue; // nothing accumulated for this tile 
							}
							nTile = tileY * tilesX + tileX;
							if (tp_tasks[iTile].getTask() == 0) continue; // nothing to do for this tile
							boolean debugTile0 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > 0); // 1);
							debugTile0 |=(tileX == debug_tileX-1) && (tileY == debug_tileY) && (globalDebugLevel > 0); // 1);
							boolean debugTile1 =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -10);
							if (debugTile0) {
								System.out.println("clt_process_tl_correlations(): tileX="+tileX+", tileY="+tileY+", iTile="+iTile+", nTile="+nTile);
							}
							// TODO: move port coordinates out of color channel loop
							// double [][] centersXY = tp_tasks[iTile].getDoubleXY(); // isAux());
							
							
							// Generate +1/+2
							double [][] disp_dist = tp_tasks[iTile].getDoubleDispDist();
							int num_combo = combine_corrs? (diameters_combo? 2 : 1) : 0;
							double [][] corrs =      new double [correlation2d.getNumPairs() + num_combo][]; // extra for combo "all"
							// copy correlation tiles from the GPU's floating point arrays 
							for (int ipair = 0; ipair < used_pairs_list.length; ipair++) {
								int pair = used_pairs_list[ipair];
								corrs[pair] = new double[corr_length];
								int fcorr2D_indx = iTile *  used_pairs_list.length + ipair;
								for (int i = 0; i < corr_length; i++) {
									corrs[pair][i] = gpu_corr_scale * fcorr2D[fcorr2D_indx][i]; // copy one-by-one converting from floats to doubles 	
								}
							}
							if (dcorr_tiles != null) {
								dcorr_tiles[iTile] =  corrs;
							}
							if (clt_corr_out != null) { // now not used
								for (int num_pair = 0; num_pair < corrs.length; num_pair++) if (corrs[num_pair] != null){
									clt_corr_out[num_pair][tileY][tileX] = corrs[num_pair];
								}
							}
							// get CM disparity/strength
							double [] disp_str = {0.0, 0.0}; // disparity = 0 will be initial approximation for LMA if no averaging
							if (combine_corrs) {
								double [] corr_combo_tile = correlation2d.accumulateInit(); // combine all available pairs
								double [] corr_dia_tile = diameters_combo ? correlation2d.accumulateInit(): null; // combine diameters
								double sumw = 0.0, sumw_dia = 0.0;
								if (imgdtt_params.mcorr_static_weights || imgdtt_params.mcorr_dynamic_weights) {
									double [] weights = new double [correlation2d.getNumPairs()];
									if (imgdtt_params.mcorr_static_weights) {
										int [] pair_lengths = correlation2d.getPairLengths();
										for (int npair = 0; npair< weights.length; npair++) if (corrs[npair] != null) {
											weights[npair] = imgdtt_params.mcorr_weights[pair_lengths[npair] - 1];
										}
									} else {
										for (int npair = 0; npair< weights.length; npair++) if (corrs[npair] != null) {
											weights[npair] = 1.0;
										}
										//									Arrays.fill(weights, 1.0);
									}
									if (imgdtt_params.mcorr_dynamic_weights) {
										for (int npair = 0; npair< weights.length; npair++) if (corrs[npair] != null) {
											double pair_width = correlation2d. getPairWidth(
													corrs[npair] , //double [] corr_tile,
													npair);        // int       num_pair)
											if (pair_width < 0.1) {
												if (globalDebugLevel > 1) {
													System.out.println("pair_width["+npair+"]="+pair_width);
												}
											} else {
												weights[npair] /= Math.pow(pair_width, imgdtt_params.mcorr_weights_power);
											}
										}

									}
									if (debugTile0) {
										System.out.println("clt_process_tl_correlations(): per-pair weights:");
										for (int npair = 0; npair< weights.length; npair++) if (weights[npair] > 0.0) {
											int [] pair = correlation2d.getPair(npair);
											int pair_length = correlation2d.getPairLength(npair);
											System.out.println(String.format("%3d: %2d -> %2d [%2d]: %8.5f",npair,pair[0],pair[1],pair_length,weights[npair]));
										}
									}
									sumw = correlation2d.accummulatePairs(
											corr_combo_tile,      // double []   accum_tile,
											corrs,                // double [][] corr_tiles, may be longer than selection, extra will be ignored 
											weights);             // double []     weights);
									if (corr_dia_tile != null) {
										double [] weights_dia = weights.clone();
										boolean [] sel_dia = correlation2d.selectDiameters(null);
										for (int i= 0; i < sel_dia.length; i++) {
											if (!sel_dia[i]) {
												weights_dia[i] = 0.0;
											}
										}
										sumw_dia = correlation2d.accummulatePairs(
												corr_dia_tile,        // double []   accum_tile,
												corrs,                // double [][] corr_tiles, may be longer than selection, extra will be ignored 
												weights_dia);         // double []     weights);
									}
								} else { // old way, same weight for all pairs
									sumw = correlation2d.accummulatePairs(
											corr_combo_tile,           // double []   accum_tile,
											corrs,                     // double [][] corr_tiles, may be longer than selection, extra will be ignored 
											correlation2d.selectAll(), // boolean []  selection,
											1.0);             // double      weight);
									if (corr_dia_tile != null) {
										sumw_dia = correlation2d.accummulatePairs(
												corr_dia_tile,           // double []   accum_tile,
												corrs,                     // double [][] corr_tiles, may be longer than selection, extra will be ignored 
												correlation2d.selectDiameters(null), // boolean []  selection,
												1.0);             // double      weight);
									}
								}

								correlation2d.normalizeAccumulatedPairs(
										corr_combo_tile,
										sumw);
								corrs[correlation2d.getNumPairs()] = corr_combo_tile;
								if (corr_dia_tile != null) {
									correlation2d.normalizeAccumulatedPairs(
											corr_dia_tile,
											sumw_dia);
									corrs[correlation2d.getNumPairs()+1] = corr_dia_tile;
								}
								
								// copy to output for monitoring if non-null;
								if ((clt_corr_out != null) && (clt_corr_out.length > num_pairs)) {
									clt_corr_out[num_pairs][tileY][tileX] = corr_combo_tile;
									if ((clt_corr_out.length > (num_pairs+1)) && (corr_dia_tile != null)) {
										clt_corr_out[num_pairs + 1][tileY][tileX] = corr_dia_tile;
									}
								}
								// calculate 0,1, or 2 maximums
								if (debugTile1) {
									System.out.println("clt_process_tl_correlations(): debugTile1, tp_task["+iTile+"]target_disparity="+tp_tasks[iTile].getTargetDisparity());
//									debugTile0 = true;
								}
								double [][] disp_str_sel = null;
								double [][] disp_str_lma = null;
								int [] sel_fg_bg = null;
								// now use older LMA (single-max) when LY is needed
								if (imgdtt_params.bimax_dual_LMA && (ddnd == null)) {
									double [][] maxes = correlation2d.getDoublePoly(
											disparity_scale, // double    disparity_scale,
											((corr_dia_tile != null) ? corr_dia_tile : corr_combo_tile), // double [] combo_corrs,
											imgdtt_params.mcorr_dual_fract,   //double    min_fraction
											imgdtt_params.mcorr_dual_min_max,  // double    min_max, // =      0.2;  // Minimal absolute strength of the strongest in a dual-max to consider second one
											imgdtt_params.mcorr_dual_min_min); // double    min_min);  // =      0.08; // Minimal absolute strength of a weakest in a dual-max to consider second one
									if ((disparity_map != null) && (disparity_map[DISPARITY_VARIATIONS_INDEX] != null)) {
										disparity_map[DISPARITY_VARIATIONS_INDEX    ][nTile] = maxes.length;
									}
									if (maxes.length > 0) {
										// TODO: add corr layer - copy of combo with singles as nulls
										// just for debugging to indicate tiles with dual maximums
										if ((maxes.length < 2) && (corr_dia_tile!=null)) { //FIXME: Debug
											//									corrs[correlation2d.getNumPairs()] = null; // temporarily keep only with pairs
											Arrays.fill(corrs[correlation2d.getNumPairs()+1], Double.NaN);
										}
										if ((maxes.length >= 2) || ! imgdtt_params.bimax_dual_only) {
											/// create disparity/strength results w/o LMA. Use it if LMA fails
											sel_fg_bg = new int[] { // first - index FG, second - index BG
													Correlation2d.selDispStrIndex( // FG
															Correlation2d.CAMEL_FG, // int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
															maxes),                           // double[][]          disp_str_dual)  // -preliminary center x in pixels for largest baseline
													Correlation2d.selDispStrIndex( // BG
															Correlation2d.CAMEL_BG, // int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
															maxes)                           // double[][]          disp_str_dual)  // -preliminary center x in pixels for largest baseline
											};

											int sel_max = Correlation2d.selDispStrIndex( // single tile
													imgdtt_params.bimax_combine_mode, // int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
													maxes);                           // double[][]          disp_str_dual)  // -preliminary center x in pixels for largest baseline

											disp_str_sel = Correlation2d.selDispStr( // single tile
													imgdtt_params.bimax_combine_mode, // int                 combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
													maxes);                           // double[][]          disp_str_dual)  // -preliminary center x in pixels for largest baseline
											int combine_mode_corrected = imgdtt_params.bimax_combine_mode;
											if (    (maxes[sel_fg_bg[0]][1] < imgdtt_params.mcorr_fb_fract * maxes[sel_fg_bg[1]][1]) ||
													(maxes[sel_fg_bg[1]][1] < imgdtt_params.mcorr_bf_fract * maxes[sel_fg_bg[0]][1])) {
												//maxes = new double [][] {maxes[0]}; // strongest
												disp_str_sel = new double [][] {maxes[0]}; // strongest
												sel_max = 0;
												sel_fg_bg = new int[] {0,0};
												combine_mode_corrected = Correlation2d.CAMEL_STRONGEST;
											}

											if (run_lma) {
												if (debugTile1) {
													System.out.println("clt_process_tl_correlations() maxes, tp_task["+iTile+"]target_disparity="+tp_tasks[iTile].getTargetDisparity());
													for (int i = 0; i < maxes.length; i++) {
														System.out.println(String.format("maxes[%d][0]=%f (quadcam disparity pixels, not combo pixels), maxes[%d][1]=%f", i, maxes[i][0], i, maxes[i][1]));
													}
												}
												int combine_mode_pre = imgdtt_params.bimax_post_LMA ? Correlation2d.CAMEL_BOTH : combine_mode_corrected;
												Corr2dLMA lma_dual = correlation2d.corrLMA2DualMax( // null pointer
														imgdtt_params,                // ImageDttParameters  imgdtt_params,
														combine_mode_pre, // 3, // 0, // 1,// int     combine_mode,   // 0 - both,  1 - strongest, 2 - nearest to zero, 3 - FG, 4 - BG
														corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
														corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
														corrs,                        // corrs,          // double [][]         corrs,
														disp_dist,
														rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
														// all that are not null in corr_tiles
														correlation2d.selectAll(),    // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
														maxes,                        //double[][]          disp_str_dual,          // -preliminary center x in pixels for largest baseline
														null, // debug_lma_tile,               // double []           debug_lma_tile,
														(debugTile0 ? 1: -2),         // int                 debug_level,
														tileX,                        // int                 tileX, // just for debug output
														tileY );
												if (debugTile1 && (lma_dual!= null)) { // addad (lma_dual != null)
													System.out.println("clt_process_tl_correlations() corrLMA2DualMax() done, lma_dual="+
															((lma_dual== null)? "null": " not null"));
												}
												if (lma_dual != null) {
													boolean dbg_dispStrs = (debug_lma != null);
													double [][][] dispStrs = lma_dual.lmaDisparityStrengths( //TODO: add parameter to filter out negative minimums ?
															dbg_dispStrs, // false, // boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
															imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
															imgdtt_params.lmas_min_amp_bg,   //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
															imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
															imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
															imgdtt_params.lmas_min_max_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
															imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
															imgdtt_params.lmas_max_area,     //double  lma_max_area,     // maximal half-area (if > 0.0)
															imgdtt_params.lma_str_scale,     // convert lma-generated strength to match previous ones - scale
															imgdtt_params.lma_str_offset,    // convert lma-generated strength to match previous ones - add to result
															dbg_dispStrs, // false // boolean dbg_mode
										    				imgdtt_params.lma_ac_offset     // Add to A, C coefficients for near-lines where A,C could become negative because of window
															);
													disp_str_lma = new double [dispStrs.length][]; // order matching input ones
													for (int nmax = 0;nmax < dispStrs.length; nmax++) {
														if ((dispStrs[nmax] != null) && (dispStrs[nmax].length >0)) {
															disp_str_lma[nmax] = dispStrs[nmax][0];
															if (dbg_dispStrs) {
																for (int i = 0; i < dispStrs[nmax][0].length; i++) {
																	debug_lma[i][nTile] = dispStrs[nmax][0][i];
																}
															}
														}
													}
													if (imgdtt_params.bimax_post_LMA && (sel_max >=0)) {
														disp_str_lma = new double [][] {disp_str_lma[sel_max]}; // was already selected before LMA
													}
												} // if (lma_dual != null) {
											}
										} // no positive maximums on y==0 - nothing in this tile
									}
									if ((disp_str_sel != null) && (disp_str_sel.length > 0)) {
										if (debugTile1) { // FIXME: remove debugTile1!
											System.out.println("clt_process_tl_correlations() disp_str_sel=");
											for (int nmax = 0; nmax < disp_str_sel.length; nmax++) {
												System.out.println(String.format("disp_str_sel[%d][0]=%f, disp_str_sel[%d][1]=%f LMA=%s",
														nmax, disp_str_sel[nmax][0], nmax, disp_str_sel[nmax][1], (((disp_str_lma!=null) && (disp_str_lma[nmax]!=null))?"true":"false")));
											}
										}
										if (disparity_map != null) {
											disparity_map[DISPARITY_INDEX_CM      ][nTile] = disp_str_sel[0][0]; // disparity non-LMA
											disparity_map[DISPARITY_INDEX_CM + 1  ][nTile] = disp_str_sel[0][1]; // strength non-LMA;
											disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = disp_str_sel[0][1]; // strength non-LMA;
											if ((disp_str_lma!=null) && (disp_str_lma.length>0)) {
												if (disp_str_lma.length == 1) {
													if (!Double.isNaN(disp_str_lma[0][0])) {
														disparity_map[DISPARITY_INDEX_POLY    ][nTile] = disp_str_lma[0][0]; // disparity LMA
														disparity_map[DISPARITY_INDEX_POLY + 1][nTile] = disp_str_lma[0][2]; // strength LMA
														// keep strength from CM strength (not LMA)
//														disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = disp_str_lma[0][1]; // overwrite with LMA strength
													}
												} else {
													int indx = sel_fg_bg[0];
													if (!Double.isNaN(disp_str_lma[indx][0])) {
														disparity_map[DISPARITY_INDEX_POLY    ][nTile] = disp_str_lma[indx][0]; // disparity LMA
														disparity_map[DISPARITY_INDEX_POLY + 1][nTile] = disp_str_lma[indx][2]; // strength LMA
														// keep strength from CM strength (not LMA)
//														disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = disp_str_lma[indx][1]; // overwrite with LMA strength
													}
												}
											}
										}
									}
								} else {  // if (imgdtt_params.bimax_dual_LMA && (ddnd == null)) else
									if (disparity_map != null) {
										int [] ixy =   correlation2d.getMaxXYInt( // find integer pair or null if below threshold // USED in lwir
												corr_combo_tile, // double [] data,      // [data_size * data_size]
												null, // disp_str_combo,
												correlation2d.getCombWidth(), //       data_width,
												correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(), // int       center_row, ??????????????
												true, // boolean   axis_only,
												-1.0, // imgdtt_params.min_corr,  // ???? double    minMax,    // minimal value to consider (at integer location, not interpolated)
												false); // debugCluster); // tile_lma_debug_level > 0); // boolean   debug);
										double [] corr_stat = correlation2d.getMaxXCm(         // get fractional center as a "center of mass" inside circle/square from the integer max
												corr_combo_tile,                      // double [] data,      // [data_size * data_size]
												correlation2d.getCombWidth(),        // int       data_width,      //  = 2 * transform_size - 1;
												correlation2d.getCombHeight()/2 - correlation2d.getCombOffset(),// int       center_row,
												ixy[0],                              // int       ixcenter,  // integer center x
												false); // debugCluster); // (tile_lma_debug_level > 0)); // boolean   debug);
										if (corr_stat != null) { // almost always
											disp_str = new double [] {-corr_stat[0], corr_stat[1]};
											disparity_map[DISPARITY_INDEX_CM      ][nTile] = disp_str[0];
											disparity_map[DISPARITY_INDEX_CM + 1  ][nTile] = disp_str[1];
											disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = disp_str[1];
										}
									}
								} // if (imgdtt_params.bimax_dual_LMA) else
							} // if (combine_corrs) {
							
							// 	New method depends on combine_corrs, so if it is disabled - use old way
							// ddnd != null - use corrLMA2Single()
							if ((!imgdtt_params.bimax_dual_LMA || !combine_corrs || (ddnd != null)) && run_lma && (disparity_map != null)) { // old way
								// debug the LMA correlations
								if (debugTile0) { // should be debugTile
									System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
								}
								double [] poly_disp = {Double.NaN, 0.0};
								double [] debug_lma_tile = (debug_lma != null) ? (new double [debug_lma.length]):null;
								if (debug_lma_tile != null) {
									for (int i = 0; i < debug_lma.length; i++) {
										debug_lma_tile[i] = debug_lma[i][nTile];
									}
								}
								boolean use_LY = imgdtt_params.lmas_LY_single || (ddnd != null);
								Corr2dLMA lma2 = correlation2d.corrLMA2Single( // null pointer
										imgdtt_params,                // ImageDttParameters  imgdtt_params,
										use_LY, // imgdtt_params.lmas_LY_single, // false,    // boolean             adjust_ly, // adjust Lazy Eye
										corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
										corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
										corrs,                        // corrs,          // double [][]         corrs,
										disp_dist,
										rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
										// all that are not null in corr_tiles
										correlation2d.selectAll(),    // longToArray(imgdtt_params.dbg_pair_mask),  // int                 pair_mask, // which pairs to process
										disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
										poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
										imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
										debug_lma_tile,               // double []           debug_lma_tile,
										(debugTile0 ? 1: -2),         // int                 debug_level,
//										-2, //0,                            // tile_lma_debug_level, // +2,         // int                 debug_level,
										tileX,                        // int                 tileX, // just for debug output
										tileY ); 
								// int                 tileY
								if (debugTile0) { // should be debugTile
									System.out.println("Ran LMA for tileX="+tileX+", tileY="+tileY);
								}
								double [][] ds = null;
								if (lma2 != null) {
									ds = lma2.lmaDisparityStrength(
						    				false, // boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
											imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
											imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
											imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
											imgdtt_params.lmas_min_max_ac,       // maximal of A and C coefficients minimum (measures sharpest point/line)
											imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
											imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
											imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
											imgdtt_params.lma_str_offset,    // convert lma-generated strength to match previous ones - add to result
						    				imgdtt_params.lma_ac_offset     // Add to A, C coefficients for near-lines where A,C could become negative because of window
											);
									if (ds != null) { // always true
//										if (disparity_map!=null) {
											disparity_map[DISPARITY_INDEX_POLY    ][nTile] = ds[0][0];
											disparity_map[DISPARITY_INDEX_POLY + 1][nTile] = ds[0][2]; // LMA strength as is
											disparity_map[DISPARITY_STRENGTH_INDEX][nTile] = ds[0][1]; // overwrite with LMA strength
//										}
										if (ddnd != null) {
											ddnd[tileY][tileX] = lma2.getDdNd();
										}
										if (debugTile0) {
											lma2.printStats(ds,1);
											if (use_LY) { // imgdtt_params.lmas_LY_single) {
//												double [][] ddnd = lma2.getDdNd(); // will not be used here
												if (ddnd != null) {
													double [][] dxy= new double [ddnd.length][2];
													for (int i = 0; i < dxy.length; i++) {
														dxy[i][0] = ddnd[tileY][tileX][i][0] * rXY[i][0] - ddnd[tileY][tileX][i][1] * rXY[i][1];
														dxy[i][1] = ddnd[tileY][tileX][i][0] * rXY[i][1] + ddnd[tileY][tileX][i][1] * rXY[i][0];
													}
													System.out.print("       Port:  ");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format("   %2d   ", i)); System.out.println();
													System.out.print("Radial_in =  [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[tileY][tileX][i][0])); System.out.println("]");
													System.out.print("Tangent_CW = [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", ddnd[tileY][tileX][i][1])); System.out.println("]");
													System.out.print("X =          [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][0])); System.out.println("]");
													System.out.print("Y =          [");
													for (int i = 0; i < dxy.length; i++) System.out.print(String.format(" %6.3f,", dxy[i][1])); System.out.println("]");
													System.out.println();
												} else {
													System.out.println("No dd/nd and x/y offsets data is available ");
												}
											} else {
												System.out.println("LY offsets are not measured");
											}
										}
									}
								}

								if (debug_lma_tile != null) {
									for (int i = 0; i < debug_lma.length; i++) {
										debug_lma[i][nTile] = debug_lma_tile[i];
									}
								}
							} // if (run_lma && (disparity_map != null))
						}
					}
				};
			}
			startAndJoin(threads);
			if (debug_lma != null) {
				if (imgdtt_params.bimax_dual_LMA) {
					(new ShowDoubleFloatArrays()).showArrays(
							debug_lma,
							tilesX,
							tilesY,
							true,
							"lma_debug_dual_LMA",
							debug_lma_titles
							);
				} else {
					(new ShowDoubleFloatArrays()).showArrays(
							debug_lma,
							tilesX,
							tilesY,
							true,
							"lma_debug",
							debug_lma_titles
							);
				}

			}			
		}
		return;
	}	
	
		
	
	
	
	
	public void clt_process_tl_correlations_GPU_DBG( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			// both arrays should have same non-null tiles
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
			final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][4*64] TD of combo corrs: qud, cross, hor,vert
			                                           // each of the top elements may be null to skip particular combo type
			final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
			final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
            // 											  [tilesY][tilesX] should be set by caller
			final float  [][][]       fcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
			// When clt_mismatch is non-zero, no far objects extraction will be attempted
			final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
			                                           // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
			final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
			                                           // last 2 - contrast, avg/ "geometric average)
			final int                 corr_pairs_mask, // which pairs to use in LMA
			final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
			final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
			final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0
			final int                 gpu_corr_rad,    // = transform_size - 1 ?
			final double              max_corr_radius, // 3.9;
			final int                 window_type,     // GPU: will not be used
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int [] dbg_tiles = (globalDebugLevel > -1)?
				new int []{37721, 39022, 40975, 38042, 40963, 42253, 38694, 39343, 36443} : null;
		final float gpu_fcorr_scale = (float) gpu_corr_scale;
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return;
		}
		final boolean [][] saturation_imp = gpuQuad.quadCLT.saturation_imp;               // boolean [][] saturation_imp, // (near) saturated pixels or null
//gpuQuad
		final boolean debug_distort= globalDebugLevel > 0; ///false; // true;

//		final double [][] debug_offsets = new double[imgdtt_params.lma_dbg_offset.length][2];
		final double [][] debug_offsets = new double[getNumSensors()][2];  
		for (int i = 0; i < imgdtt_params.lma_dbg_offset.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}
//dbg_pair_mask
		final int quad = 4;   // number of subcameras
//		final int numcol = isMonochrome()?1:3;

		final int width =  gpuQuad.getImageWidth();
		final int height = gpuQuad.getImageHeight();
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
//		final double [] col_weights= new double [numcol]; // colors are RBG
		final double [][] dbg_distort = debug_distort? (new double [4*quad][tilesX*tilesY]) : null;
		// not yet used with GPU
/**		
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],inv_pwr), imgdtt_params.lma_min_wnd);
				}
			}
		}
*/
		// keep for now for mono, find out  what do they mean for macro mode
		
		final int corr_size = transform_size * 2 - 1;
		

		// reducing weight of on-axis correlation values to enhance detection of vertical/horizontal lines
		// multiply correlation results inside the horizontal center strip  2*enhortho_width - 1 wide by enhortho_scale

		final double [] enh_ortho_scale = new double [corr_size];
		for (int i = 0; i < corr_size; i++){
			if ((i < (transform_size - imgdtt_params.getEnhOrthoWidth(isAux()))) || (i > (transform_size - 2 + imgdtt_params.getEnhOrthoWidth(isAux())))) {
				enh_ortho_scale[i] = 1.0;
			} else {
				enh_ortho_scale[i] = imgdtt_params.getEnhOrthoScale(isAux());
			}
			if (i == (transform_size-1)) enh_ortho_scale[i] = 0.0 ; // hardwired 0 in the center
			enh_ortho_scale[i] *= Math.sin(Math.PI*(i+1.0)/(2*transform_size));
		}
		if (globalDebugLevel > 1){
			System.out.println("getEnhOrthoWidth(isAux())="+ imgdtt_params.getEnhOrthoWidth(isAux())+" getEnhOrthoScale(isAux())="+ imgdtt_params.getEnhOrthoScale(isAux()));
			for (int i = 0; i < corr_size; i++){
				System.out.println(" enh_ortho_scale["+i+"]="+ enh_ortho_scale[i]);

			}
		}

		// Create window  to select center correlation strip using
		// ortho_height - full width of non-zero elements
		// ortho_eff_height - effective height (ration of the weighted column sum to the center value)
		int wcenter = transform_size - 1;
		final double [] ortho_weights = new double [corr_size]; // [15]
		for (int i = 0; i < corr_size; i++){
			if ((i >= wcenter - imgdtt_params.ortho_height/2) && (i <= wcenter + imgdtt_params.ortho_height/2)) {
				double dx = 1.0*(i-wcenter)/(imgdtt_params.ortho_height/2 + 1);
				ortho_weights[i] = 0.5*(1.0+Math.cos(Math.PI*dx))/imgdtt_params.ortho_eff_height;
			}
		}
		if (globalDebugLevel > 1){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}


		if (globalDebugLevel > 1) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		
		if (disparity_map != null){
			for (int i = 0; i<disparity_map.length;i++) if ((disparity_modes & (1 << i)) != 0){
				if ((i == OVEREXPOSED) && (saturation_imp == null)) {
					continue;
				}
				disparity_map[i] = new double [tilesY*tilesX];
			}
		}
		
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}

		/*
		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256] - never used
		final double [] lt_window2 = new double [lt_window.length]; // squared - never used
		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];
		*/
///		final boolean use_main = geometryCorrection_main != null;
///	    final int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		
		int [] corr_indices_ = null; 
		float [][] fcorr2D_ = null;
		
		if (fcorr_td != null) {
			int [][] pairs_map = {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5}};
			corr_indices_ = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
					fcorr_td, // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
					pairs_map); // int [][] pairs) // typically {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5} [0] - 3rd index in corr_tiles, [1] -
			gpuQuad.execCorr2D_normalize(
					false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
					gpu_corr_rad); // int corr_radius
			fcorr2D_ = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
		}
		
		final int [] corr_indices = corr_indices_; 
		final float [][] fcorr2D = fcorr2D_; // len = 49311

		final int num_combo= (fcorr_combo_td != null)? fcorr_combo_td.length : 0;
		final int [][] corr_combo_indices = (fcorr_combo_td != null)? new int [num_combo][] : null;
		final float [][][] fcorr2D_combo = (fcorr_combo_td != null)? new float[num_combo][][] : null;
		if (num_combo > 0) {
			int [] ipairs = {0xf, 0x30, 0x3, 0xc}; // does not matter. using quad, cross, hor, vert
			for (int i = 0; i < num_combo; i++) {
				corr_combo_indices[i] = gpuQuad.setCorrTilesComboTd(
						fcorr_combo_td[i], // final float [][][] corr_tiles, // [tileY][tileX][4*64]
						ipairs[i]); // change to smth? int ipair) // just to set in the index low bits
				// normalize and convert to pixel domain
				gpuQuad.execCorr2D_normalize(
						true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
						gpu_fat_zero, // double fat_zero);
						null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
						gpu_corr_rad); // int corr_radius
				fcorr2D_combo[i] =  gpuQuad.getCorr2DCombo(gpu_corr_rad);
			}
		}
		final boolean fneed_macro = false;
		if (corr_indices.length > 0) {
			/*
			if (true) { // debugging only
				int [] wh = new int[2];
				double [][] dbg_corr = GPUTileProcessor.getCorr2DView(
						tilesX,
						tilesY,
						corr_indices,
						fcorr2D,
						wh);
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_corr,
						wh[0],
						wh[1],
						true,
						"dbg-corr2D-fromTD", // name+"-CORR2D-D"+clt_parameters.disparity,
						GPUTileProcessor.getCorrTitles());
			}
			 */
			final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size
			// assuming that the correlation pairs sets are the same for each tile that has correlations
			// find this number
			int nt0 = (corr_indices[0] >> GPUTileProcessor.CORR_NTILE_SHIFT);
			int nc0 = 1;
			for (int i = 1; (i < corr_indices.length) && ((corr_indices[i] >> GPUTileProcessor.CORR_NTILE_SHIFT) == nt0) ; i++) {
				nc0++;
			}
			final int num_tile_corr = nc0; // normally 6
			final int num_tiles = corr_indices.length / num_tile_corr;
			// FIXME: will not work with combining pairs !!!
			final int num_pairs = Correlation2d.getNumPairs(numSensors);
//			Runtime.getRuntime().gc();
//			System.out.println("--- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {

						Correlation2d corr2d = new Correlation2d(
								numSensors,
								imgdtt_params,              // ImageDttParameters  imgdtt_params,
								transform_size,             // int transform_size,
								2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
								isMonochrome(), // boolean monochrome,
								(globalDebugLevel > -1));   //   boolean debug)
						corr2d.createOrtoNotch(
								imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
								imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
								(imgdtt_params.lma_debug_level > 1)); // boolean debug);

						for (int indx_tile = ai.getAndIncrement(); indx_tile < num_tiles; indx_tile = ai.getAndIncrement()) {
							// double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
							// added quad and cross combos
							double [][]  corrs = new double [num_pairs + num_combo][corr_length]; // 225-long (15x15)
							float  [][] fcorrs = (fcorr_tiles == null) ? null : new float  [num_pairs + num_combo][corr_length]; // 225-long (15x15)
							int indx_corr = indx_tile * num_tile_corr;
							int nt = (corr_indices[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
							int tileX = nt % tilesX;
							int tileY = nt / tilesX;
							int tIndex = tileY * tilesX + tileX;

							// Prepare the same (currently 10-layer) corrs as double [][], as in CPU version
							int pair_mask = 0;
							if (fcorr_td != null) {
								for (int indx_pair = 0; indx_pair < num_tile_corr; indx_pair++) {
									int pair = corr_indices[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
									assert pair < num_pairs : "invalid correllation pair";
									pair_mask |= (1 << pair);
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D[indx_corr][i]; // from float to double
									}
									if (fcorrs != null) for (int i = 0; i < corr_length; i++) {
										fcorrs[pair][i] = gpu_fcorr_scale * fcorr2D[indx_corr][i];
									}
									indx_corr++; 
								}
							}
							// add 4 combo layers : quad, cross, hor, vert
							if (num_combo > 0) {
								for (int ncm = 0; ncm < num_combo; ncm++) if (corr_combo_indices[ncm]!=null){
									nt = (corr_combo_indices[ncm][indx_tile] >> GPUTileProcessor.CORR_NTILE_SHIFT); // corr_quad_indices - different sequence
									int pair = num_pairs + ncm; // 6+
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D_combo[ncm][indx_tile][i]; // from float to double
									}
									if (fcorrs != null) for (int i = 0; i < corr_length; i++) {
										fcorrs[pair][i] = gpu_fcorr_scale * fcorr2D_combo[ncm][indx_tile][i];
									}
								}
							}
							if (corr_tiles != null) {
								corr_tiles[tileY][tileX] = corrs; 
							}
							if (fcorr_tiles != null) {
								fcorr_tiles[tileY * tilesX + tileX] = fcorrs; // does not require corr_common_GPU()
							}
							if ((disparity_map != null) || (clt_corr_partial != null) || (clt_mismatch != null)) {
//								int used_pairs = pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int used_pairs = pair_mask & corr_pairs_mask; // imgdtt_params.dbg_pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
								boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
								if (dbg_tiles != null) for (int ii : dbg_tiles){
									if (ii == tIndex) {
										debugTile=true;
										System.out.println("Debugging tile "+tIndex+" ( "+tileX+" / "+tileY+" )");
										break;
									}
								}
								corr_common_GPU(
										imgdtt_params,        // final ImageDttParameters  imgdtt_params,
										clt_corr_partial,     // final double [][][][][]   clt_corr_partial,
										used_pairs,           // final int           used_pairs,
										disparity_map,        // final double [][]   disparity_map,
										clt_mismatch,         // final double [][]   clt_mismatch,
										saturation_imp,       // final boolean [][]  saturation_imp,
										fneed_macro,          // final boolean       fneed_macro,
										corr2d,               // final Correlation2d corr2d,
										corrs,                // final double [][]   corrs,
										tileX,                // final int           tileX,
										tileY,                // final int           tileY,
										max_corr_radius,      // final double        max_corr_radius, // 3.9;
										tile_lma_debug_level, // int                 tile_lma_debug_level,
										debugTile,            // boolean             debugTile,
										globalDebugLevel);    // final int           globalDebugLevel)							
							}
							// double extra_disparity = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
							/*
// Disabled for GPU
								if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
								else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
								else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
								else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];  // not used in lwir
								else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];  // not used in lwir
								if (Double.isNaN(extra_disparity)) extra_disparity = 0;  // used in lwir
							 */
							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 3. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}

							// only debug is left
							// old (per-color correlation)
							// removed

						} // end of tile
					}
				};
			}
			startAndJoin(threads);
		} else {
			// no correlation tiles to process
		}
		if ((dbg_distort != null) &&(globalDebugLevel >= 0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}
/*
//		final double [][] dbg_distort = debug_distort? (new double [4*quad][tilesX*tilesY]) : null;
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
			(new ShowDoubleFloatArrays()).showArrays(dbg_distort,  tilesX, tilesY, true, "disparity_distortions"); // , dbg_titles);
		}

		if (dbg_ports_coords != null) {
			(new showDoubleFloatArrays()).showArrays(dbg_ports_coords,  tilesX, tilesY, true, "ports_coordinates", dbg_titles);
		}
*/
//		return clt_data;
	}
	
	
	
	public void corr_common_GPU(
			final ImageDttParameters  imgdtt_params,
			final double [][][][][]   clt_corr_partial,
			final int           used_pairs,
			final double [][]   disparity_map,
			final double [][]   clt_mismatch,
			final boolean [][]  saturation_imp,
			final boolean       fneed_macro,
			final Correlation2d corr2d,
			final double [][]   corrs,
			final int           tileX,
			final int           tileY,
			final double        max_corr_radius, // 3.9; only used with clt_mismatch
			int                 tile_lma_debug_level,
			boolean             debugTile,
			final int           globalDebugLevel)
	{
		final int quad = 4;   // number of subcameras
		final int numcol = isMonochrome()?1:3;
		// does not include combo
//		int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
//		boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
		// non-GPU initialization of the data structures
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
//		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		int tIndex = tileY * tilesX + tileX;
		
		final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;
		if (disparity_map != null) {
			for (int i = 0; i < disparity_map.length; i++) {
				if (disparity_map[i] != null) {
					if ((((1 << i) & BITS_FROM_GPU) == 0) || !fneed_macro) { // do not touch data already set up
						disparity_map[i][tIndex] = (
								(i == DISPARITY_STRENGTH_INDEX) ||
								(i == DISPARITY_INDEX_HOR_STRENGTH) ||
								(i == DISPARITY_INDEX_VERT_STRENGTH)) ? 0.0 : Double.NaN; // once and for all
					}
				}
			}
			// calculate overexposed fraction
			if (saturation_imp != null){
				// not yet implemented in GPU
				disparity_map[OVEREXPOSED][tIndex] = 0.0; // (1.0 * overexp_all[0]) / overexp_all[1];
			}
			//clt_mismatch should only be used with disparity_map != null;
			if (clt_mismatch != null) {
				for (int np = 0; np < clt_mismatch.length/3; np++) {
					clt_mismatch[3 * np + 0 ][tIndex] = Double.NaN;
					clt_mismatch[3 * np + 1 ][tIndex] = Double.NaN;
					clt_mismatch[3 * np + 2 ][tIndex] = 0;
				}
			}
		}


		// calculate all selected pairs correlations
		//int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		// Code that was after correlations calculation

		double [][] strips = corr2d.scaleRotateInterpoateCorrelations(
				corrs,                          // double [][] correlations,
				used_pairs,                      // int         pairs_mask,
				imgdtt_params.corr_strip_hight, //);    // int         hwidth);
				(tile_lma_debug_level > 0) ? used_pairs:0); // debugMax);

		// Combine strips for selected pairs. Now using only for all available pairs.
		// Other combinations are used only if requested (clt_corr_partial != null)

		double [] strip_combo = corr2d.combineInterpolatedCorrelations(
				strips,                        // double [][] strips,
				used_pairs,                     // int         pairs_mask,
				imgdtt_params.corr_offset,     // double      offset);
				imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)

		double [][] strips_intra = null;
		double [] strip_combo_intra = null;
		if (corrs.length > 6) {
			strips_intra = new double[2][];
			strips_intra[0] = corr2d.scaleRotateInterpoateSingleCorrelation(
					corrs[6], // double []   corr,  // quad 
					imgdtt_params.corr_strip_hight, // int         hwidth,
					Correlation2d.PAIR_HORIZONTAL,  // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
					1,                              // int         ss, // 1
					false // boolean     debug
					); 
			strips_intra[1] = corr2d.scaleRotateInterpoateSingleCorrelation(
					corrs[7], // double []   corr,  // quad 
					imgdtt_params.corr_strip_hight, // int         hwidth,
					Correlation2d.PAIR_DIAGONAL_MAIN,  // int         dir, // 0 - hor, 1 - vert, 2 - parallel to row = col (main) diagonal (0->3), 3 -2->1
					1,                              // int         ss, // 1
					false // boolean     debug
					); 
			strip_combo_intra = corr2d.combineInterpolatedCorrelations(
					strips_intra,                  // double [][] strips,
					3,                     // int         pairs_mask,
					imgdtt_params.corr_offset,     // double      offset);
					imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)
		}
		// Debug feature - only calculated if requested
		if ((clt_corr_partial != null) && (imgdtt_params.corr_mode_debug || imgdtt_params.gpu_mode_debug)) {
			@SuppressWarnings("unused")
			double [] strip_ortho = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x0f,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)
			@SuppressWarnings("unused")
			double [] strip_diag = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x30,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)
			@SuppressWarnings("unused")
			double [] strip_all = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x3f,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

			double [] strip_hor = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x03,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

			double [] strip_vert = corr2d.combineInterpolatedCorrelations(
					strips,                         // double [][] strips,
					0x0c,                           // int         pairs_mask,
					imgdtt_params.corr_offset,      // double      offset);
					imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

			// re-using arrays that were made for color channels
			clt_corr_partial[tileY][tileX] = new double[quad][numcol+1][];
			clt_corr_partial[tileY][tileX][0][0] = corrs[0];                        // 1
			clt_corr_partial[tileY][tileX][0][1] = corrs[1];                        // 2
			clt_corr_partial[tileY][tileX][0][2] = corrs[2];                        // 3
			clt_corr_partial[tileY][tileX][0][3] = corrs[3];                        // 4
			clt_corr_partial[tileY][tileX][1][0] = corrs[4];                        // 5
			clt_corr_partial[tileY][tileX][1][1] = corrs[5];                        // 6
			clt_corr_partial[tileY][tileX][1][2] = corrs[6];                        // 7
			clt_corr_partial[tileY][tileX][1][3] = corrs[7];                        // 8
			//												    	clt_corr_partial[tileY][tileX][1][2] = corrs_ortho;                     // 7
			//												    	clt_corr_partial[tileY][tileX][1][3] = corrs_cross;                     // 8
			//												    	clt_corr_partial[tileY][tileX][1][2] = corr2d.debugStrip(strip_hor);    // 7
			//												    	clt_corr_partial[tileY][tileX][1][3] = corr2d.debugStrip(strip_vert);   // 8
			//strip_combo_intra						    	
			clt_corr_partial[tileY][tileX][2][0] = corrs[8];                        // 9
			clt_corr_partial[tileY][tileX][2][1] = corrs[9];                        // 10
			//												    	clt_corr_partial[tileY][tileX][2][0] = corr2d.debugStrip(strips[4]);    // 9
			//												    	clt_corr_partial[tileY][tileX][2][1] = corr2d.debugStrip(strips[5]);    // 10
			clt_corr_partial[tileY][tileX][2][2] = corr2d.debugStrip2(strip_hor);   // 11
			clt_corr_partial[tileY][tileX][2][3] = corr2d.debugStrip2(strip_vert);  // 12
			if (strips_intra != null) {
				clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strips_intra[0]); // 13
				clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strips_intra[1]);  // 14
			}
			//												    	clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strip_ortho); // 13
			//												    	clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strip_diag);  // 14
			if (strip_combo_intra != null) {
				clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip2(strip_combo_intra);    // 15
			}
			//												    	clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip(strip_all);    // 15
			clt_corr_partial[tileY][tileX][3][3] = corr2d.debugStrip2(strip_combo); // 16
		}
		if (imgdtt_params.pcorr_use && (strip_combo_intra != null)) {
			strip_combo = strip_combo_intra;
		}
		// calculate CM maximums for all mixed channels
		// First get integer correlation center, relative to the center
		for (int i = 0; i < strip_combo.length; i++) if (Double.isNaN(strip_combo[i])){
			strip_combo[i] = 0.0; // ????
		}
		int [] ixy =  corr2d.getMaxXYInt( // find integer pair or null if below threshold
				strip_combo,              // double [] data,
				true,                     // boolean   axis_only,
				imgdtt_params.min_corr,   //  double    minMax,    // minimal value to consider (at integer location, not interpolated)
				tile_lma_debug_level > 0); // boolean   debug);

		double [] corr_stat = null;

		// if integer argmax was strong enough, calculate CM argmax
		// will not fill out DISPARITY_INDEX_INT+1, DISPARITY_INDEX_CM+1, DISPARITY_INDEX_POLY+1
		// use clt_mismatch for that
		//							double strength = 0.0;
		//							double disparity = 0.0;
		double [] disp_str = new double[2];
		if (ixy != null) {
			disp_str[1] = strip_combo[ixy[0]+transform_size-1]; // strength at integer max on axis
			if (disparity_map != null) {
				disparity_map[DISPARITY_INDEX_INT][tIndex] =      -ixy[0];
				disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = disp_str[1];
				if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
					System.out.println("BUG: 1. disparity_map[DISPARITY_STRENGTH_INDEX]["+tIndex+"] should not be NaN");
				}
			}
			corr_stat = corr2d.getMaxXCm(   // get fractional center as a "center of mass" inside circle/square from the integer max
					strip_combo,                      // double [] data,      // [data_size * data_size]
					ixy[0],                           // int       ixcenter,  // integer center x
					// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
					// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
					(tile_lma_debug_level > 0)); // boolean   debug);
		}
		if (disparity_map != null) {
			if (imgdtt_params.pcorr_use_hv) {
				// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
				double [] hor_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						0x100, //  corrs[8] Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
						imgdtt_params.corr_offset,          // double      corr_offset,
						true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						false,                              // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0);          // boolean   debug);
				if (hor_pair1 != null){
					disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
					disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
				}

				double [] vert_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						0x200, // corrs[9] Correlation2d.getMaskVertical(1), // int         pairs_mask,
						imgdtt_params.corr_offset,        // double      corr_offset,
						true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						true, // not anymore transposed false, // already transposed  // true,                             // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0); // boolean   debug);
				if (vert_pair1 != null) {
					disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
					disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
				}
			} else  {								
				// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
				double [] hor_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
						imgdtt_params.corr_offset,          // double      corr_offset,
						true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						false,                              // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0);          // boolean   debug);
				if ((hor_pair1 != null)) {
					disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
					disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
				}

				double [] vert_pair1 = corr2d.getMaxXSOrtho(
						corrs,                              // double [][] correlations,
						Correlation2d.getMaskVertical(1), // int         pairs_mask,
						imgdtt_params.corr_offset,        // double      corr_offset,
						true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						true,                             // boolean     is_vert,      // transpose X/Y
						debugTile); // tile_lma_debug_level > 0); // boolean   debug);
				if (vert_pair1 != null) {
					disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
					disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
				}
			}
		}
		// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
		if (corr_stat != null) {
			// skipping DISPARITY_VARIATIONS_INDEX - it was not used
			disp_str[0] = -corr_stat[0];
			if (disparity_map != null) {
				disparity_map[DISPARITY_INDEX_CM][tIndex] = disp_str[0]; // disparity is negative X
			}
			if (tile_lma_debug_level > 0) {
				System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
			}

			// debug new LMA correlations
			if (debugTile) {
				System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY+"\n\n\n*** Disabled for the GPU ***\n\n\n");
				// temporarily disabling for the GPU (can be restored as disp_dist is available)
				/*
				double [] poly_disp = {Double.NaN, 0.0};
				Corr2dLMA lma2 = corr2d.corrLMA2(
						imgdtt_params,                // ImageDttParameters  imgdtt_params,
						corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
						corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
						corrs,                        // double [][]         corrs,
						disp_dist,
						rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
						imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
						disp_str,  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
						poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
						imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
						tile_lma_debug_level+2,         // int                 debug_level,
						tileX,                        // int                 tileX, // just for debug output
						tileY );                      // int                 tileY
				double [][] ds = null;
				if (lma2 != null) {
					ds = lma2.lmaDisparityStrength(
							imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
							imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
							imgdtt_params.lmas_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
							imgdtt_params.lmas_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
							imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
							imgdtt_params.lma_str_offset    // convert lma-generated strength to match previous ones - add to result
							);
					lma2.printStats(ds,1);
				}
				 */
			}

			//								disparity_map[DISPARITY_INDEX_CM + 1][tIndex] = // y not available here
			// calculate/fill out hor and vert
			// convert to multi-baseline combining results from several integer scales

			// see if strength is enough to proceed with LMA/poly (otherwise keep disp/strength
			if (disp_str[1] > imgdtt_params.min_poly_strength) {
				// create LMA instance, calculate LMA composite argmax
				// Create 2 groups: ortho & diag
				Correlations2dLMA lma;
				if (imgdtt_params.pcorr_use) { // new group phase correlation
					double [][] fake_corrs = {corrs[6],null,null,null,corrs[7],null};
					lma = corr2d.corrLMA(
							imgdtt_params,                // ImageDttParameters  imgdtt_params,
							fake_corrs,                   // double [][]         corrs,
							corr2d.longToArray(0x11), // 0x11, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
							false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
							corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
							imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
							tile_lma_debug_level,         // int                 debug_level,
							tileX,                        // int                 tileX, // just for debug output
							tileY );                      // int                 tileY
				} else {
					lma = corr2d.corrLMA(
							imgdtt_params,                // ImageDttParameters  imgdtt_params,
							corrs,                        // double [][]         corrs,
							corr2d.longToArray(used_pairs), // used_pairs, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
							false,                        // boolean             run_poly_instead, // true - run LMA, false - run 2d polynomial approximation
							corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
							imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
							tile_lma_debug_level,         // int                 debug_level,
							tileX,                        // int                 tileX, // just for debug output
							tileY );                      // int                 tileY
				}
				double [] lma_disparity_strength = null;
				double max_disp_diff_lma = 3.0;
				if (lma != null) {
					lma_disparity_strength = lma.getDisparityStrength();
					if (Math.abs(lma_disparity_strength[0] - disp_str[0] ) > max_disp_diff_lma) {
						if (globalDebugLevel > -1) {
							System.out.println("Crazy LMA for tileX="+tileX+", tileY="+tileY+": disparity="+disp_str[0]+",lma_disparity_strength[0]="+lma_disparity_strength[0]);
						}
						lma = null;
					}
				}
				if (lma != null) {
					double []   mod_disparity_diff = null;
					double [][] dir_corr_strength =  null;
					lma_disparity_strength = lma.getDisparityStrength();
					if (tile_lma_debug_level > 0){
						System.out.println(String.format("Tile X/Y = %d/%d LMA disparity = %7.4f, strength = %7.4f",
								tileX, tileY,
								lma_disparity_strength[0],lma_disparity_strength[1]));
					}
					if (disparity_map != null) {
						disparity_map[DISPARITY_INDEX_POLY]         [tIndex] = lma_disparity_strength[0];
					}

					// if enabled overwrite - replace  DISPARITY_INDEX_CM and DISPARITY_STRENGTH_INDEX
					if (imgdtt_params.mix_corr_poly) { //true
						disp_str[0] =  lma_disparity_strength[0];
						disp_str[1] =  lma_disparity_strength[1];
						if (disparity_map != null) {
							disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
							disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
							if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
								System.out.println("BUG: 2. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
							}
						}
					}
					// store debug data
					// if strong enough and enabled - try to improve far objects
					// temporarily removed strength requirement for debugging, restore later to make faster
					///										if ((imgdtt_params.fo_correct && (strength > imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {

					if ((imgdtt_params.fo_correct && (disp_str[1] > 0 * imgdtt_params.fo_min_strength)) || (clt_mismatch != null)) {
						// try all dirs:
						dir_corr_strength = corr2d.corr4dirsLMA(
								imgdtt_params,                // ImageDttParameters  imgdtt_params,
								corrs,                        // double [][]         corrs,
								used_pairs, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								-disp_str[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								tile_lma_debug_level,         // int                 debug_level,
								tileX,                        // int                 tileX, // just for debug output
								tileY );                      // int                 tileY
						if ((tile_lma_debug_level > 0) && (dir_corr_strength != null)) {
							double [] nan2 = {Double.NaN, Double.NaN, Double.NaN, Double.NaN};
							for (int ii = 0; ii < dir_corr_strength.length; ii++) {
								if (dir_corr_strength[ii] == null) dir_corr_strength[ii] = nan2;
							}
							System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⇖:%7.4f (%7.4f) ⇗:%7.4f (%7.4f)",
									dir_corr_strength[0][0],dir_corr_strength[0][1],dir_corr_strength[1][0],dir_corr_strength[1][1],
									dir_corr_strength[2][0],dir_corr_strength[2][1],dir_corr_strength[3][0],dir_corr_strength[3][1]));
						}

						mod_disparity_diff =     corr2d.foregroundCorrect(
								imgdtt_params.fo_far,            // boolean   bg,
								imgdtt_params.fo_ortho,          // boolean   ortho,
								dir_corr_strength,               // double [] dir_disp,
								disp_str[0],                     // double    full_disp,
								imgdtt_params.fo_min_strength,   // double      min_strength,
								imgdtt_params.fo_min_eff,        // double      min_eff,
								imgdtt_params.fo_min_eff_ratio,  // double      min_eff_ratio,
								imgdtt_params.fo_max_hwidth,    // double      max_hwidth, //  =          3.0;  // maximal half-width disparity  direction to try to correct far objects
								imgdtt_params.fo_min_diff,       // double    fo_min_diff,
								imgdtt_params.fo_overcorrection, // double    fo_overcorrection,
								imgdtt_params.fo_lim_overcorr,   // double    fo_lim_overcorr,
								(tile_lma_debug_level > 0) );    // boolean debug);

						// Do not use modified far object distance when mismatch is measured
						if ((mod_disparity_diff[0] != disp_str[0]) && (clt_mismatch == null)){ // if it changed
							if (imgdtt_params.fo_correct && (disp_str[1] > imgdtt_params.fo_min_strength)) { // always
								disp_str[0] = mod_disparity_diff[0];
								if (disparity_map != null) {
									disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
								}
							}
						}
					}
					if (tile_lma_debug_level > -1) {
						System.out.println("debug12348973591");
					}
					if (clt_mismatch != null) { // mod_disparity_diff should be calculated
						// bypass difference or zero strength if disparity difference is too high (will influence mismatch correction)
						// but setting it too low will make it impossible to correct larger mismatches. Maybe multi-pass?
						if (mod_disparity_diff[2] <= imgdtt_params.mismatch_max_diff) { // may be NaN, will fail test as intended
							if (tile_lma_debug_level > -1) {
								System.out.println("debug12348973590");
							}
							double [] mismatch_result = null;
							boolean need_CM = true;
							if (imgdtt_params.ly_poly) {  // not used in lwir
								mismatch_result = corr2d.mismatchPairs( // returns x-xcenter, y, strength (sign same as disparity)
										imgdtt_params,                // ImageDttParameters  imgdtt_params,
										corrs,                        // double [][]         corrs,
										used_pairs,                    // int                 pair_mask, // which pairs to process
										-disp_str[0],                   // double    xcenter,   // preliminary center x in pixels for largest baseline
										max_corr_radius,              // double    vasw_pwr,  // value as weight to this power,
										tile_lma_debug_level,// int                 debug_level,
										tileX,         // int                 tileX, // just for debug output
										tileY );       // int                 tileY
								// check if got crazy poly, then retry with CM
								boolean has_NaN = false;
								boolean need_dbg = false;
								for (int dir = 0; dir < (mismatch_result.length/3); dir ++) {
									if (Double.isNaN(mismatch_result[3*dir + 0]) || Double.isNaN(mismatch_result[3*dir + 1])) {
										has_NaN = true;
									} else if ((mismatch_result[3*dir + 2] != 0.0) &&
											((Math.abs(mismatch_result[3*dir + 0]) > imgdtt_params.ly_crazy_poly) ||
													(Math.abs(mismatch_result[3*dir + 1]) > imgdtt_params.ly_crazy_poly))) {
										mismatch_result[3*dir + 2] = 0;
										has_NaN = true;
										need_dbg = true;
									}
								}
								need_CM = imgdtt_params.ly_poly_backup && has_NaN;
								if (need_dbg && (imgdtt_params.lma_debug_level > 0)) {
									System.out.println("Failed polynomial mismatch for tileX="+tileX+", tileY="+tileY);
									for (int dir = 0; dir < (mismatch_result.length/3); dir ++) {
										System.out.println(String.format("%d: dxy[%d]=%f8.5, dxy[%d]=%f8.5 strength=%7.5f",
												dir, (dir*2+0), mismatch_result[dir*3 + 0], (dir*2 + 1), mismatch_result[dir*3 + 1], mismatch_result[dir*3 + 2]));
									}
								}
							}
							// TODO: use magic_scale for CM?
							if (need_CM) { // if poly was off or gave crazy poly
								mismatch_result = corr2d.mismatchPairsCM( // returns x-xcenter, y, strength (sign same as disparity)
										imgdtt_params,                // ImageDttParameters  imgdtt_params,
										corrs,                        // double [][]         corrs,
										used_pairs,                    // int                 pair_mask, // which pairs to process
										-disp_str[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
										max_corr_radius, // imgdtt_params.ortho_vasw_pwr, // radius,    // positive - within that distance, negative - within 2*(-radius)+1 square
										tile_lma_debug_level,// int                 debug_level,
										tileX,         // int                 tileX, // just for debug output
										tileY );       // int                 tileY

								if (imgdtt_params.ly_poly && (imgdtt_params.lma_debug_level > 0)) {
									System.out.println("Corrected by CM failed polynomial mismatch for tileX="+tileX+", tileY="+tileY);
									for (int dir = 0; dir < (mismatch_result.length/3); dir ++) {
										System.out.println(String.format("%d: dxy[%d]=%f8.5, dxy[%d]=%f8.5 strength=%7.5f",
												dir, (dir*2+0), mismatch_result[dir*3 + 0], (dir*2 + 1), mismatch_result[dir*3 + 1], mismatch_result[dir*3 + 2]));
									}
								}
							}
							if (tile_lma_debug_level > 0) {
								System.out.println("Lazy eye mismatch:");
								for (int np = 0; np < mismatch_result.length/3; np++) {
									System.out.println(String.format("%2d: dx = %7.4f, dy = %7.4f, strength = %7.4f,",
											np, mismatch_result[3 * np + 0], mismatch_result[3 * np + 1], mismatch_result[3 * np + 2]));
								}
							}

							for (int np = 0; np < clt_mismatch.length/3; np++) if (np < mismatch_result.length/3){
								clt_mismatch[3 * np + 0 ][tIndex] = mismatch_result[3 * np + 0 ];
								clt_mismatch[3 * np + 1 ][tIndex] = mismatch_result[3 * np + 1 ];
								clt_mismatch[3 * np + 2 ][tIndex] = mismatch_result[3 * np + 2 ];
							}
						}
					}
				} else { // if (lma != null)
					if (imgdtt_params.corr_poly_only) { // discard tile if LMA failed
						disp_str[0] =  Double.NaN;
						disp_str[1] =  0.0;
						if (disparity_map != null) {
							disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
							disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
						}
					}
				}
			}
		} // end of if (corr_stat != null)
		// double extra_disparity = 0.0; // used for textures:  if allowed, shift images extra before trying to combine
		/*
//Disabled for GPU
			if      (corr_mode == 0) extra_disparity = disparity_map[DISPARITY_INDEX_INT][tIndex];
			else if (corr_mode == 1) extra_disparity = disparity_map[DISPARITY_INDEX_CM][tIndex];
			else if (corr_mode == 2) extra_disparity = disparity_map[DISPARITY_INDEX_POLY][tIndex];
			else if (corr_mode == 3) extra_disparity = disparity_map[DISPARITY_INDEX_HOR][tIndex];  // not used in lwir
			else if (corr_mode == 4) extra_disparity = disparity_map[DISPARITY_INDEX_VERT][tIndex];  // not used in lwir
			if (Double.isNaN(extra_disparity)) extra_disparity = 0;  // used in lwir
		 */
	} 
	
	
	
	
	public float [][][][]  blur_corr_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
//			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
//			final int                 debug_tileX,
//			final int                 debug_tileY,
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return null;
		}
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final int numTiles = tilesX*tilesY;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		
		final float [][][][] fcorr_td_out = new float [tilesY][tilesX][][];
		final float [] fweights = {0.125f, 0.0625f, 0.125f, 0.0625f,0.125f, 0.0625f,0.125f, 0.0625f, 0.25f};  // sum = 1.0

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					for (int indx_tile = ai.getAndIncrement(); indx_tile < numTiles; indx_tile = ai.getAndIncrement()) {
						float s = 0;
						float [] weights = new float [9];
						int corr_len = 0;
						int num_pairs = 0;
						for (int dir = 0; dir < 9; dir++) { // 8 - center
							int indx = tn.getNeibIndex(indx_tile, dir);
							if (indx >=0) {
								int [] xy = tn.getXY(indx);
								if (fcorr_td[xy[1]][xy[0]] != null) {
									float fw = fweights[dir];
									s += fw;
									weights[dir] = fw;
									if (num_pairs == 0) {
										num_pairs = fcorr_td[xy[1]][xy[0]].length;
										corr_len = fcorr_td[xy[1]][xy[0]][0].length;
									}
								}
							}
						}
						if (s > 0.0) {
							for (int i = 0; i < weights.length; i++) {
								weights[i] /= s;
							}
							float [][] tile_corrs = new float [num_pairs][corr_len];
							for (int dir = 0; dir < 9; dir++) if (weights[dir] > 0.0f){ // 8 - center
								int [] xy = tn.getXY(tn.getNeibIndex(indx_tile, dir));
								for (int np = 0; np < num_pairs; np++) {
									float [] ft = fcorr_td[xy[1]][xy[0]][np];
									for (int i = 0; i < corr_len; i++) {
										tile_corrs[np][i] += weights[dir] * ft[i];
									}
								}
							}
							int [] xy = tn.getXY(indx_tile);
							fcorr_td_out[xy[1]][xy[0]] = tile_corrs;
						}
					} // end of tile
				}
			};
		}
		startAndJoin(threads);
		return fcorr_td_out;
	}
	
	
	public float [][][][]  blur_corr_combo_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
			final float  [][][][]     fcorr_combo_td,  // [n][tilesY][tilesX][4*64] transform domain representation of combined corr pairs
			final int                 threadsMax,      // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		if (this.gpuQuad == null) {
			System.out.println("clt_aberrations_quad_corr_GPU(): this.gpuQuad is null, bailing out");
			return null;
		}
		final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
		final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final int numTiles = tilesX*tilesY;
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		
		final float [][][][] fcorr_td_combo_out = new float [fcorr_combo_td.length][][][];//[tilesY][tilesX][][];
		for (int nl = 0; nl < fcorr_combo_td.length; nl++) {
			if (fcorr_combo_td[nl] != null) fcorr_td_combo_out[nl] = new float [tilesY][tilesX][];
		}
		final float [] fweights = {0.125f, 0.0625f, 0.125f, 0.0625f,0.125f, 0.0625f,0.125f, 0.0625f, 0.25f};  // sum = 1.0

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					TileNeibs tn = new TileNeibs(tilesX,tilesY);
					for (int indx_tile = ai.getAndIncrement(); indx_tile < numTiles; indx_tile = ai.getAndIncrement()) {
						for (int nl = 0; nl < fcorr_combo_td.length; nl++) if (fcorr_combo_td[nl] != null){
							float [][][] fcorr_combo = fcorr_combo_td[nl];
							float s = 0;
							float [] weights = new float [9];
							int corr_len = 0;
							for (int dir = 0; dir < 9; dir++) { // 8 - center
								int indx = tn.getNeibIndex(indx_tile, dir);
								if (indx >=0) {
									int [] xy = tn.getXY(indx);
									if (fcorr_combo[xy[1]][xy[0]] != null) {
										float fw = fweights[dir];
										s += fw;
										weights[dir] = fw;
										if (corr_len == 0) {
											corr_len = fcorr_combo[xy[1]][xy[0]].length;
										}
									}
								}
							}
							if (s > 0.0) {
								for (int i = 0; i < weights.length; i++) {
									weights[i] /= s;
								}
								float [] tile_corrs = new float [corr_len];
								for (int dir = 0; dir < 9; dir++) if (weights[dir] > 0.0f){ // 8 - center
									int [] xy = tn.getXY(tn.getNeibIndex(indx_tile, dir));
									float [] ft = fcorr_combo[xy[1]][xy[0]];
									for (int i = 0; i < corr_len; i++) {
										tile_corrs[i] += weights[dir] * ft[i];
									}
								}
								int [] xy = tn.getXY(indx_tile);
								fcorr_td_combo_out[nl][xy[1]][xy[0]] = tile_corrs;
							}
						}
					} // end of tile
				}
			};
		}
		startAndJoin(threads);
		return fcorr_td_combo_out;
	}
	
	public double [][] cltMeasureLazyEyeGPU ( // returns d,s lazy eye parameters
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
			                                           // if exactly zero - considered to be infinity capturing
			final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr
			final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
			final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
			
			final double              gpu_corr_scale,  // 1.0 now! 0.75; // reduce GPU-generated correlation values
			final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
			
			final GeometryCorrection  geometryCorrection,
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
			final double              disparity_corr, // disparity at infinity
			final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
			final int                 super_radius, // 0 - none, 1 - 3x3, 2 - 5x5, (2)
			final int                 debug_tileX,
			final int                 debug_tileY,
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 globalDebugLevel)
	{
		final int quad =      4; // number of subcameras
		final int num_pairs = 6;
		final boolean use_main = geometryCorrection_main != null;
		final double [][] rXY = geometryCorrection.getRXY(use_main);
		final int         tilesX=gpuQuad.getTilesX(); // width/transform_size;
		final int         tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
		final int         clustersX= (tilesX + tileStep - 1) / tileStep;
		final int         clustersY= (tilesY + tileStep - 1) / tileStep;
		final double [][] lazy_eye_data = new double [clustersY*clustersX][];
		final int         gpu_corr_rad = transform_size - 1;
		final int nClustersInChn=clustersX * clustersY;

		final int debug_clustX = debug_tileX / tileStep;
		final int debug_clustY = debug_tileY / tileStep;
		
		// calculate which tiles to use for each cluster
		// will generate sparse array for cluster central tiles to match CPU software
		final float  [][][][]     fcorr_td_centers = new float [tilesY][tilesX][][]; // sparse, only in cluster centers
		final int    [][]         num_in_cluster = new int [clustersY][clustersX];  // only in cluster centers
		final double [][][][]     disp_dist = new double [clustersY][clustersX][][];
		final double [][][]       clust_pY =  new double  [clustersY][clustersX][];
		final double [][][]       pxpy = new double [clustersY][clustersX][];
		final double [][]         disparity_array_center = new double [clustersY][clustersX];
		final boolean [][]        bg_cluster = new boolean [clustersY][clustersX];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double shiftX = 0.0;
		final double shiftY = 0.0;
		// TODO: Maybe calculate full energy in each TD tile for normalization
		// First pass merge correlation result for each cluster
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nCluster = ai.getAndIncrement(); nCluster < nClustersInChn; nCluster = ai.getAndIncrement()) {
						int clustY = nCluster / clustersX;
						int clustX = nCluster % clustersX;
						boolean debugCluster =  (clustX == debug_clustX) && (clustY == debug_clustY);
						///						boolean debugCluster1 = (Math.abs(clustX - debug_clustX) < 10) && (Math.abs(clustY - debug_clustY) < 10);
						if (debugCluster) {
							System.out.println("debugCluster1");
						}
// filter only tiles with similar disparity to enable lazy eye for the ERS.
						int num_good_tiles = 0;
						double avg= 0.0;
						boolean has_bg = false;
						for (int cTileY = 0; (cTileY < tileStep)  && !has_bg; cTileY++) {
							int tileY = clustY * tileStep + cTileY ;
							if (tileY < tilesY) {
								for (int cTileX = 0; (cTileX < tileStep) && !has_bg; cTileX++) {
									int tileX = clustX * tileStep + cTileX ;
									if ((tileX < tilesX) && (fcorr_td[tileY][tileX] != null)) {
										double d = disparity_array [tileY][tileX];
										if (d == 0) {
											has_bg = true;
											break;
										}
									}
								}
							}
						}
						bg_cluster[clustY][clustX] = has_bg;
						float [][] ftd = new float [num_pairs][4* transform_size * transform_size];

						// select only tiles with close enough disparity values (now 5 pix)
						while (true) {
							int mnTx = -1, mnTy = -1, mxTx = -1, mxTy = -1;
							double mn = Double.NaN;
							double mx = Double.NaN;
							num_good_tiles = 0;
							avg= 0.0;
							for (int cTileY = 0; cTileY < tileStep; cTileY++) {
								int tileY = clustY * tileStep + cTileY ;
								if (tileY < tilesY) {
									for (int cTileX = 0; cTileX < tileStep; cTileX++) {
										int tileX = clustX * tileStep + cTileX ;
										if ((tileX < tilesX) && (fcorr_td[tileY][tileX] != null)) {
											double d = disparity_array [tileY][tileX];
											if (has_bg && (d != 0)) {
												continue; // for bg tiles do not mix with non-bg
											}
											avg += d;
											if (!(d <= mx)) {
												mx = d;
												mxTx = tileX;
												mxTy = tileY;
											}
											if (!(d >= mn)) {
												mn = d;
												mnTx = tileX;
												mnTy = tileY;
											}
											num_good_tiles++;
										}
									}
								}
							}
							avg /= num_good_tiles;
							if (num_good_tiles ==0) {
								break;
							}
							if ((mx-mn) <= imgdtt_params.lma_disp_range ) {
								break;
							}
							if ((mx-avg) > (avg-mn)) {
								fcorr_td[mxTy][mxTx] = null;
							} else {
								fcorr_td[mnTy][mnTx] = null;
							}
							//imgdtt_params.lma_disp_range
						}
						// for now - num_good_tiles - the only strength measure
						// should it be calculated from the correlation normalization?
						// will try to calculate through the sum of squared normalizations
						if (num_good_tiles > 0) {
							double dscale = 1.0/num_good_tiles;

							// average fdisp_dist and fpxpy over remaining tiles 
							//fpxpy
							double []   avg_py = new double [quad]; 
							double [][] avg_disp_dist = new double [quad][4]; 
							//							double [][] avg_pxpy = new double [quad][2];
							double [] avg_pxpy = new double [2];
							for (int cTileY = 0; cTileY < tileStep; cTileY++) {
								int tileY = clustY * tileStep + cTileY ;
								if (tileY < tilesY) {
									for (int cTileX = 0; cTileX < tileStep; cTileX++) {
										int tileX = clustX * tileStep + cTileX ;
										if ((tileX < tilesX) && (fcorr_td[tileY][tileX] != null)) {
											for (int nc = 0; nc < quad; nc++) {
												for (int i = 0; i < 4; i++) {
													avg_disp_dist[nc][i] += fdisp_dist[tileY][tileX][nc][i];
												}
												avg_py[nc] += fpxpy[tileY][tileX][nc][1]; // averaging some tiles, but it will be the same for all channels
											}
											double centerX = tileX * transform_size + transform_size/2 - shiftX;
											double centerY = tileY * transform_size + transform_size/2 - shiftY;
											avg_pxpy[0] += centerX;
											avg_pxpy[1] += centerY;
										}
									}
								}
							}

							num_in_cluster        [clustY][clustX] = num_good_tiles;
							disparity_array_center[clustY][clustX] = avg;

							for (int nc = 0; nc < quad; nc++) {
								for (int i = 0; i < 4; i++) {
									avg_disp_dist[nc][i] *= dscale;
								}
								avg_py[nc] *= dscale;
							}
							for (int i = 0; i < 2; i++) {
								avg_pxpy[i] *= dscale;
							}

							disp_dist[clustY][clustX] = avg_disp_dist;
							pxpy     [clustY][clustX] = avg_pxpy;
							clust_pY [clustY][clustX] = avg_py; // to use for ERS
							// accumulate TD tiles
							// If needed - save average
							for (int cTileY = 0; cTileY < tileStep; cTileY++) {
								int tileY = clustY * tileStep + cTileY ;
								if (tileY < tilesY) {
									for (int cTileX = 0; cTileX < tileStep; cTileX++) {
										int tileX = clustX * tileStep + cTileX ;
										if ((tileX < tilesX) && (fcorr_td[tileY][tileX] != null)) {
											for (int np = 0; np <num_pairs; np++) {
												for (int i = 0; i < ftd[np].length; i++) {
													ftd[np][i]+=fcorr_td[tileY][tileX][np][i];
												}
											}
										}
									}
								}
							}
						}
						if (num_good_tiles > 0) {
							int centerY = clustY * tileStep + tileStep/2; // integer was 242!
							int centerX = clustX * tileStep + tileStep/2; // integer
							if (centerY >= tilesY) {
								centerY = tilesY - 1;
							}
							if (centerX >= tilesX) {
								centerX = tilesX - 1;
							}
							float fscale = 1.0f/num_good_tiles;
							for (int np = 0; np <num_pairs; np++) {
								for (int i = 0; i < ftd[np].length; i++) {
									ftd[np][i] *= fscale;
								}
							}
							fcorr_td_centers[centerY][centerX]= ftd;				
						}
					} // end of cluster
				}
			};
		}
		startAndJoin(threads);
		
				int [][] pairs_map = {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5}};
		final int [] corr_indices = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
				fcorr_td_centers, // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
				pairs_map); // int [][] pairs) // typically {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5} [0] - 3rd index in corr_tiles, [1] -
		gpuQuad.execCorr2D_normalize(
				false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
				gpu_fat_zero, // double fat_zero);
				null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
				gpu_corr_rad); // int corr_radius
		final float [][] fcorr2D = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
		
		final int corr_length = fcorr2D[0].length;// all correlation tiles have the same size
		final int num_tiles = corr_indices.length / num_pairs; 
		final double [][] corr_wnd = Corr2dLMA.getCorrWnd(
				transform_size,
				imgdtt_params.lma_wnd);
		final double [] corr_wnd_inv_limited = (imgdtt_params.lma_min_wnd <= 1.0)?  new double [corr_wnd.length * corr_wnd[0].length]: null;
		if (corr_wnd_inv_limited != null) {
			double inv_pwr = imgdtt_params.lma_wnd_pwr - (imgdtt_params.lma_wnd - 1.0); // compensate for lma_wnd
			for (int i = imgdtt_params.lma_hard_marg; i < (corr_wnd.length - imgdtt_params.lma_hard_marg); i++) {
				for (int j = imgdtt_params.lma_hard_marg; j < (corr_wnd.length - imgdtt_params.lma_hard_marg); j++) {
					corr_wnd_inv_limited[i * (corr_wnd.length) + j] = 1.0/Math.max(Math.pow(corr_wnd[i][j],
							inv_pwr),
							imgdtt_params.lma_min_wnd);
				}
			}
		}
		final double [][] dbg_img = debug_strengths? (new double[19][clustersX*clustersY]):null;
		// Second pass - process normalized per-cluster correlations
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					Correlation2d corr2d = new Correlation2d(
							numSensors,
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)
					corr2d.createOrtoNotch(
							imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
							imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
							(imgdtt_params.lma_debug_level > 1)); // boolean debug);
					for (int indx_tile = ai.getAndIncrement(); indx_tile < num_tiles; indx_tile = ai.getAndIncrement()) {
						double [][]  corrs = new double [num_pairs][corr_length]; // 225-long (15x15)
						int indx_corr = indx_tile * num_pairs;
						int nt = (corr_indices[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
						int tileX = nt % tilesX;
						int tileY = nt / tilesX;
						int clustX = tileX/tileStep;
						int clustY = tileY/tileStep;
						int nclust = clustX + clustY * clustersX;
						if (dbg_img != null) dbg_img[0][nclust] = 1.0;
						double []  disp_str =  null;

						for (int indx_pair = 0; indx_pair < num_pairs; indx_pair++) {
							int pair = corr_indices[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
							assert pair < num_pairs : "invalid correllation pair";
							for (int i = 0; i < corr_length; i++) {
								corrs[pair][i] = gpu_corr_scale * fcorr2D[indx_corr][i]; // from float to double
							}
							indx_corr++; 
						}
						//			final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
						double [][] tile_disp_dist = disp_dist[clustY][clustX];
						
						// debug new LMA correlations
						boolean debugCluster =  (clustX == debug_clustX) && (clustY == debug_clustY);
						if (debugCluster) {
							System.out.println("debugCluster2");
						}
						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? imgdtt_params.lma_debug_level : -1;
						int tdl = debugCluster ? tile_lma_debug_level : -3;
						if (debugCluster && (globalDebugLevel > -1)) { // -2)) {
							System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
						}
						double [] poly_disp = {Double.NaN, 0.0};
						Corr2dLMA lma2 = corr2d.corrLMA2Single( // multitile num_tiles == 1
								imgdtt_params,                // ImageDttParameters  imgdtt_params,
					    		true,                         // boolean             adjust_ly, // adjust Lazy Eye
								corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
								corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
								corrs,                        // double [][]         corrs,
								tile_disp_dist,
								rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
								corr2d.longToArray(imgdtt_params.dbg_pair_mask), // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								null,                         // disp_str[cTile],  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
								imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								tdl, // tile_lma_debug_level, //+2,         // int                 debug_level,
								tileX,                        // int                 tileX, // just for debug output
								tileY);                       // int                 tileY
						
						
						/*
						double [][] poly_disp2 = {{Double.NaN, 0.0}};
						double [][][] corrs2 = {corrs};
						double [][][] tile_disp_dist2 = {tile_disp_dist};
						// TODO: maybe use corrLMA2Single again, but take care of initVector!
						Corr2dLMA lma2 = corr2d.corrLMA2Multi( // multitile num_tiles == 1
								imgdtt_params,                // ImageDttParameters  imgdtt_params,
								1,                            // int                 clust_width,
								corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
								corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
								corrs2, // corrs,                        // double [][]         corrs,
								tile_disp_dist2,
								rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
								imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
//								null,                         // disp_str[cTile],  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								poly_disp2,                    // double[]            poly_ds,    // null or pair of disparity/strength
								imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								tdl, // tile_lma_debug_level, //+2,         // int                 debug_level,
								tileX,                        // int                 tileX, // just for debug output
								tileY);                       // int                 tileY
						*/
						if (lma2 != null) {
							if (dbg_img != null) dbg_img[1][nclust] = 1.0;
							// was for single tile
							disp_str = lma2.lmaDisparityStrength(
				    				false, // boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
				    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
									imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
									imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
									imgdtt_params.lmas_min_max_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
									imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
									imgdtt_params.lmas_max_area,     // double  lma_max_area,     // maximal half-area (if > 0.0)
									imgdtt_params.lma_str_scale,     // convert lma-generated strength to match previous ones - scale
									imgdtt_params.lma_str_offset,     // convert lma-generated strength to match previous ones - add to result
				    				imgdtt_params.lma_ac_offset     // Add to A, C coefficients for near-lines where A,C could become negative because of window
									)[0];
							if (tile_lma_debug_level > 0) {
								double [][] ds_dbg = {disp_str};
								lma2.printStats(ds_dbg,1);
							}
							double [][] ddnd = lma2.getDdNd();
							double [] stats  = lma2.getStats (num_in_cluster[clustY][clustX]);
							if (dbg_img != null) {
								dbg_img[2][nclust] = stats[0];
								dbg_img[3][nclust] = stats[1];
								dbg_img[4][nclust] = stats[2];
							}

		    				double [][] lma_ds = lma2.lmaDisparityStrengthLY( // [1][2]
		    						imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
		    						imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
		    						imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
									imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
									imgdtt_params.lma_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
				    				1.0, // imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
				    				0.0); // imgdtt_params.lma_str_offset);  // convert lma-generated strength to match previous ones - add to result
							
							double strengh_k = 1.0; // 0.2*Math.sqrt(num_in_cluster[clustY][clustX]); // approximately matching old/multitile
							if (dbg_img != null) {
								double [][] dbg_ext_stat = lma2.lmaGetExtendedStats(
										imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
										imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
										imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
										imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
										imgdtt_params.lma_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
										1.0, // imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
										0.0); // imgdtt_params.lma_str_offset);  // convert lma-generated strength to match previous ones - add to result
								for (int ii = 0; ii < dbg_ext_stat[0].length; ii++) {
									dbg_img[5+ii][nclust] = dbg_ext_stat[0][ii];
								}
								dbg_img[16][nclust] = num_in_cluster[clustY][clustX];
								dbg_img[17][nclust] = strengh_k * lma_ds[0][1]  * num_in_cluster[clustY][clustX]; 
								dbg_img[18][nclust] = lma_ds[0][0] + disparity_array_center[clustY][clustX] + disparity_corr; 
							}
							if ((lma_ds[0] != null) && (lma_ds[0][1]> 0.0)) {
								lazy_eye_data[nclust] = new double [ExtrinsicAdjustment.get_INDX_LENGTH(numSensors)];
								lazy_eye_data[nclust][ExtrinsicAdjustment.INDX_STRENGTH] =           strengh_k * lma_ds[0][1]  * num_in_cluster[clustY][clustX]; 
								lazy_eye_data[nclust][ExtrinsicAdjustment.INDX_DISP] =               lma_ds[0][0] + disparity_array_center[clustY][clustX] + disparity_corr;
								lazy_eye_data[nclust][ExtrinsicAdjustment.INDX_TARGET] =             (disparity_array_center[clustY][clustX] + disparity_corr);
								lazy_eye_data[nclust][ExtrinsicAdjustment.INDX_DIFF]  =              lma_ds[0][0];
								lazy_eye_data[nclust][ExtrinsicAdjustment.INDX_PX + 0] =             pxpy[clustY][clustX][0];
								lazy_eye_data[nclust][ExtrinsicAdjustment.INDX_PX + 1] =             pxpy[clustY][clustX][1];
								for (int cam = 0; cam < quad; cam++) {
									lazy_eye_data[nclust][ExtrinsicAdjustment.get_INDX_DYDDISP0(numSensors) + cam] = tile_disp_dist[cam][2];
									lazy_eye_data[nclust][ExtrinsicAdjustment.get_INDX_PYDIST(numSensors) + cam] =   clust_pY [clustY][clustX][cam];
									
								}
								for (int cam = 0; cam < ddnd.length; cam++) {
									if (ddnd[cam] != null) { //convert to x,y from dd/nd
										lazy_eye_data[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] = ddnd[cam][0] * rXY[cam][0] - ddnd[cam][1] * rXY[cam][1];
										lazy_eye_data[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] = ddnd[cam][0] * rXY[cam][1] + ddnd[cam][1] * rXY[cam][0];
										lazy_eye_data[nclust][ExtrinsicAdjustment.get_INDX_DD0(numSensors) + cam] =        ddnd[cam][0];
										lazy_eye_data[nclust][ExtrinsicAdjustment.get_INDX_ND0(numSensors) + cam] =        ddnd[cam][1];
									} else {
										lazy_eye_data[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] = Double.NaN;
										lazy_eye_data[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] = Double.NaN;
										lazy_eye_data[nclust][ExtrinsicAdjustment.get_INDX_DD0(numSensors) + cam] =        Double.NaN;
										lazy_eye_data[nclust][ExtrinsicAdjustment.get_INDX_ND0(numSensors) + cam] =        Double.NaN;
									}
								}
							}
						}
					} // end of tile
				}
			};
		}
		startAndJoin(threads);
		
		if (dbg_img != null) {
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					clustersX,
					clustersY,
					true,
					"ly_dbg"); // name+"-CORR2D-D"+clt_parameters.disparity,
		}
		if (super_radius == 0) {
			return lazy_eye_data; // no processing of clouds in the sky
		}
		// save a copy of 		
		final double [][]         lazy_eye_data_final = new double [clustersY*clustersX][];
		final int    [][]         num_in_cluster_final = new int [clustersY][clustersX];  // only in cluster centers
		final float  [][][][]     fcorr_td_super = new float [tilesY][tilesX][][]; // sparse, only in cluster centers
		final double [][][][]     disp_dist_super = new double [clustersY][clustersX][][];
		final double [][][]       clust_pY_super =  new double  [clustersY][clustersX][];
		final double [][][]       pxpy_super = new double [clustersY][clustersX][];
		
		// pass 3 - prepare superclusters for weak infinity (clouds)
		ai.set(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nCluster = ai.getAndIncrement(); nCluster < nClustersInChn; nCluster = ai.getAndIncrement()) {
						int clustY = nCluster / clustersX;
						int clustX = nCluster % clustersX;
						boolean debugCluster =  (clustX == debug_clustX) && (clustY == debug_clustY);
						///						boolean debugCluster1 = (Math.abs(clustX - debug_clustX) < 10) && (Math.abs(clustY - debug_clustY) < 10);
						if (debugCluster) {
							System.out.println("debugCluster3");
						}
						if (lazy_eye_data[nCluster] != null) {
							lazy_eye_data_final[nCluster] =  lazy_eye_data[nCluster];
							num_in_cluster_final[clustY][clustX] = num_in_cluster[clustY][clustX]; 
						} else if (bg_cluster[clustY][clustX]) { // only for empty (too weak) infinity
							// only process empty infinity clusters
							// filter only tiles with similar disparity to enable lazy eye for the ERS.
///							int num_good_super = 0; // number of clusters - remove?
							int num_good_tiles = 0;
//							boolean [][] usable_clust = new boolean [2 * super_radius + 1][2 * super_radius + 1]; 
//							int [][] centersX = new int [2 * super_radius + 1][2 * super_radius + 1]; 
//							int [][] centersY = new int [2 * super_radius + 1][2 * super_radius + 1];
							// accumulate averages
							float [][]  ftd =           new float [num_pairs][4* transform_size * transform_size];
							double []   avg_py =        new double [quad]; 
							double [][] avg_disp_dist = new double [quad][4]; 
							double []   avg_pxpy =      new double [2];
							
							for (int dClustY = -super_radius; dClustY <= super_radius; dClustY++) {
								int iClustY = clustY + dClustY;
								if ((iClustY >= 0) && (iClustY< clustersY)) {
									for (int dClustX = -super_radius; dClustX <= super_radius; dClustX++) {
										int iClustX = clustX + dClustX;
										if ((iClustX >= 0) && (iClustX < clustersX)) {
											int iCluster = iClustY * clustersX + iClustX;
											int iCenterY = iClustY * tileStep + tileStep/2; // integer was 242!
											int iCenterX = iClustX * tileStep + tileStep/2; // integer
											if (iCenterY >= tilesY) {
												iCenterY = tilesY - 1;
											}
											if (iCenterX >= tilesX) {
												iCenterX = tilesX - 1;
											}
											if (bg_cluster[clustY][clustX] &&
													(lazy_eye_data[iCluster] == null) &&
													(fcorr_td_centers[iCenterY][iCenterX] != null)) {
//												usable_clust[dClustY + super_radius][dClustX + super_radius] = true;
///												num_good_super++;
												num_good_tiles += num_in_cluster[iClustY][iClustX];
												double dscale = num_in_cluster[iClustY][iClustX];
												float fscale = (float) dscale;
												for (int np = 0; np <num_pairs; np++) {
													for (int i = 0; i < ftd[np].length; i++) {
														ftd[np][i]+= fscale * fcorr_td_centers[iCenterY][iCenterX][np][i];
													}
												}
												for (int nc = 0; nc < quad; nc++) {
													for (int i = 0; i < 4; i++) {
														avg_disp_dist[nc][i] += dscale * disp_dist[iClustY][iClustX][nc][i];
													}
													avg_py[nc] += dscale * clust_pY[clustY][clustX][nc];
												}
												for (int i = 0; i < 2; i++) {
													avg_pxpy[i] += dscale * pxpy [iClustY][iClustX][i] ;
												}
												
											}
										}
									}
								}
							}
							num_in_cluster_final[clustY][clustX] = num_good_tiles;
				            double dscale = 1.0/num_good_tiles;
				            float fscale = (float) dscale;
							for (int np = 0; np <num_pairs; np++) {
								for (int i = 0; i < ftd[np].length; i++) {
									ftd[np][i] *= fscale;
								}
							}
							for (int nc = 0; nc < quad; nc++) {
								for (int i = 0; i < 4; i++) {
									avg_disp_dist[nc][i] *= dscale;
								}
								avg_py[nc] *= dscale;
							}
							for (int i = 0; i < 2; i++) {
								avg_pxpy[i] *= dscale;
							}
							disp_dist_super[clustY][clustX] = avg_disp_dist;
							pxpy_super     [clustY][clustX] = avg_pxpy;
							clust_pY_super [clustY][clustX] = avg_py; // to use for ERS
							int centerY = clustY * tileStep + tileStep/2; // integer was 242!
							int centerX = clustX * tileStep + tileStep/2; // integer
							if (centerY >= tilesY) {
								centerY = tilesY - 1;
							}
							if (centerX >= tilesX) {
								centerX = tilesX - 1;
							}
							fcorr_td_super[centerY][centerX]= ftd;				
						} // else if (bg_cluster[clustY][clustX]) { // only for empty (too weak) infinity
					} // end of cluster
				}
			};
		}
		startAndJoin(threads);
//		1) process  fcorr_td_super with GPU
		final int [] corr_indices_super = gpuQuad.setCorrTilesTd( // .length = 295866 should set num_corr_tiles!
				fcorr_td_super, // final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
				pairs_map); // int [][] pairs) // typically {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5} [0] - 3rd index in corr_tiles, [1] -
		gpuQuad.execCorr2D_normalize(
				false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
				gpu_fat_zero, // double fat_zero);
				null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
				gpu_corr_rad); // int corr_radius
		final float [][] fcorr2D_super = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
		final int corr_length_super = fcorr2D_super[0].length + 0;// all correlation tiles have the same size
		final int num_tiles_super = corr_indices_super.length / num_pairs; 
		final double [][] dbg_img2 = debug_strengths? (new double[19][clustersX*clustersY]):null;

		// Fourth pass - process normalized per-cluster correlations for low-contrast infinity tiles (clouds in the sky)
		ai.set(0);
// TODO: 	2) calculate lazy_eye_data_final 

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					Correlation2d corr2d = new Correlation2d(
							numSensors,
							imgdtt_params,              // ImageDttParameters  imgdtt_params,
							transform_size,             // int transform_size,
							2.0,                        //  double wndx_scale, // (wndy scale is always 1.0)
							isMonochrome(), // boolean monochrome,
							(globalDebugLevel > -1));   //   boolean debug)
					corr2d.createOrtoNotch(
							imgdtt_params.getEnhOrthoWidth(isAux()), // double getEnhOrthoWidth(isAux()),
							imgdtt_params.getEnhOrthoScale(isAux()), //double getEnhOrthoScale(isAux()),
							(imgdtt_params.lma_debug_level > 1)); // boolean debug);
					for (int indx_tile = ai.getAndIncrement(); indx_tile < num_tiles_super; indx_tile = ai.getAndIncrement()) {
						double [][]  corrs = new double [num_pairs][corr_length_super]; // 225-long (15x15)
						int indx_corr = indx_tile * num_pairs;
						int nt = (corr_indices_super[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
						int tileX = nt % tilesX;
						int tileY = nt / tilesX;
						int clustX = tileX/tileStep;
						int clustY = tileY/tileStep;
						int nclust = clustX + clustY * clustersX;
						if (dbg_img2 != null) {
							dbg_img2[0][nclust] = 1.0;
						}
						double []  disp_str =  null;

						for (int indx_pair = 0; indx_pair < num_pairs; indx_pair++) {
							int pair = corr_indices_super[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
							assert pair < num_pairs : "invalid correllation pair";
							for (int i = 0; i < corr_length_super; i++) {
								corrs[pair][i] = gpu_corr_scale * fcorr2D_super[indx_corr][i]; // from float to double
							}
							indx_corr++; 
						}
						//			final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
						double [][] tile_disp_dist = disp_dist_super[clustY][clustX];
						
						// debug new LMA correlations
						boolean debugCluster =  (clustX == debug_clustX) && (clustY == debug_clustY);
						if (debugCluster) {
							System.out.println("debugCluster2");
						}
						int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? imgdtt_params.lma_debug_level : -1;
						int tdl = debugCluster ? tile_lma_debug_level : -3;
						if (debugCluster && (globalDebugLevel > -1)) { // -2)) {
							System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY);
						}
						double [] poly_disp = {Double.NaN, 0.0};
						Corr2dLMA lma2 = corr2d.corrLMA2Single( // multitile num_tiles_super == 1
								imgdtt_params,                // ImageDttParameters  imgdtt_params,
					    		true,                         // boolean             adjust_ly, // adjust Lazy Eye
								corr_wnd,                     // double [][]         corr_wnd, // correlation window to save on re-calculation of the window
								corr_wnd_inv_limited,         // corr_wnd_limited, // correlation window, limited not to be smaller than threshold - used for finding max/convex areas (or null)
								corrs,                        // double [][]         corrs,
								tile_disp_dist,
								rXY,                          // double [][]         rXY, // non-distorted X,Y offset per nominal pixel of disparity
								corr2d.longToArray(imgdtt_params.dbg_pair_mask), // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
								null,                         // disp_str[cTile],  //corr_stat[0],                 // double    xcenter,   // preliminary center x in pixels for largest baseline
								poly_disp,                    // double[]            poly_ds,    // null or pair of disparity/strength
								imgdtt_params.ortho_vasw_pwr, // double    vasw_pwr,  // value as weight to this power,
								tdl, // tile_lma_debug_level, //+2,         // int                 debug_level,
								tileX,                        // int                 tileX, // just for debug output
								tileY);                       // int                 tileY
						
						if (lma2 != null) {
							// was for single tile
							disp_str = lma2.lmaDisparityStrength(
				    				false, // boolean bypass_tests,     // keep even weak for later analysis. Normally - only in test mode
				    				imgdtt_params.lmas_min_amp,      //  minimal ratio of minimal pair correlation amplitude to maximal pair correlation amplitude
									imgdtt_params.lmas_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
									imgdtt_params.lmas_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
									imgdtt_params.lmas_min_max_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
									imgdtt_params.lmas_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
									imgdtt_params.lmas_max_area,     // double  lma_max_area,     // maximal half-area (if > 0.0)
									imgdtt_params.lma_str_scale,     // convert lma-generated strength to match previous ones - scale
									imgdtt_params.lma_str_offset,     // convert lma-generated strength to match previous ones - add to result
				    				imgdtt_params.lma_ac_offset     // Add to A, C coefficients for near-lines where A,C could become negative because of window
									)[0];
							if (tile_lma_debug_level > 0) {
								double [][] ds_dbg = {disp_str};
								lma2.printStats(ds_dbg,1);
							}
							double [][] ddnd = lma2.getDdNd();
							double [] stats  = lma2.getStats (num_in_cluster_final[clustY][clustX]);
							double k = 2.5;
		    				double [][] lma_ds = lma2.lmaDisparityStrengthLY( // [1][2]
		    						imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
		    						imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
		    						(1/k)*imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
		    						(1/k)*imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
									k* imgdtt_params.lma_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
				    				1.0, // imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
				    				0.0); // imgdtt_params.lma_str_offset);  // convert lma-generated strength to match previous ones - add to result
							double strengh_k = 1.0; // 0.2*Math.sqrt(num_in_cluster[clustY][clustX]); // approximately matching old/multitile
							strengh_k /= (2 * super_radius + 1)*(2 * super_radius + 1);
		    				if (dbg_img2 != null) {
		    					dbg_img2[1][nclust] = 1.0;
		    					dbg_img2[2][nclust] = stats[0];
		    					dbg_img2[3][nclust] = stats[1];
		    					dbg_img2[4][nclust] = stats[2];

		    					double [][] dbg_ext_stat = lma2.lmaGetExtendedStats(
		    							imgdtt_params.lma_max_rel_rms,  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
		    							imgdtt_params.lma_min_strength, // minimal composite strength (sqrt(average amp squared over absolute RMS)
		    							(1/k)*imgdtt_params.lma_min_ac,       // minimal of A and C coefficients maximum (measures sharpest point/line)
		    							(1/k)*imgdtt_params.lma_min_min_ac,   // minimal of A and C coefficients minimum (measures sharpest point)
		    							k* imgdtt_params.lma_max_area,      //double  lma_max_area,     // maximal half-area (if > 0.0)
		    							1.0, // imgdtt_params.lma_str_scale,    // convert lma-generated strength to match previous ones - scale
		    							0.0); // imgdtt_params.lma_str_offset);  // convert lma-generated strength to match previous ones - add to result
		    					for (int ii = 0; ii < dbg_ext_stat[0].length; ii++) {
		    						dbg_img2[5+ii][nclust] = dbg_ext_stat[0][ii];
		    					}
		    					dbg_img2[16][nclust] = num_in_cluster_final[clustY][clustX];
		    					dbg_img2[17][nclust] = strengh_k * lma_ds[0][1]  * num_in_cluster_final[clustY][clustX]; 
		    					dbg_img2[18][nclust] = lma_ds[0][0] + disparity_array_center[clustY][clustX] + disparity_corr;
		    				}

							if ((lma_ds[0] != null) && (lma_ds[0][1]> 0.0)) {
								lazy_eye_data_final[nclust] = new double [ExtrinsicAdjustment.get_INDX_LENGTH(numSensors)];
								lazy_eye_data_final[nclust][ExtrinsicAdjustment.INDX_STRENGTH] =           strengh_k * lma_ds[0][1]  * num_in_cluster_final[clustY][clustX]; 
								lazy_eye_data_final[nclust][ExtrinsicAdjustment.INDX_DISP] =               lma_ds[0][0] + disparity_array_center[clustY][clustX] + disparity_corr;
								lazy_eye_data_final[nclust][ExtrinsicAdjustment.INDX_TARGET] =             disparity_array_center[clustY][clustX] + disparity_corr;
								lazy_eye_data_final[nclust][ExtrinsicAdjustment.INDX_DIFF] =               lma_ds[0][0];
								lazy_eye_data_final[nclust][ExtrinsicAdjustment.INDX_PX + 0] =             pxpy_super[clustY][clustX][0];
								lazy_eye_data_final[nclust][ExtrinsicAdjustment.INDX_PX + 1] =             pxpy_super[clustY][clustX][1];
								for (int cam = 0; cam < quad; cam++) {
									lazy_eye_data_final[nclust][ExtrinsicAdjustment.get_INDX_DYDDISP0(numSensors) + cam] = tile_disp_dist[cam][2];
									lazy_eye_data_final[nclust][ExtrinsicAdjustment.get_INDX_PYDIST(numSensors) + cam] =   clust_pY_super [clustY][clustX][cam];
								}
								for (int cam = 0; cam < ddnd.length; cam++) {
									if (ddnd[cam] != null) { //convert to x,y from dd/nd
										lazy_eye_data_final[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] = ddnd[cam][0] * rXY[cam][0] - ddnd[cam][1] * rXY[cam][1];
										lazy_eye_data_final[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] = ddnd[cam][0] * rXY[cam][1] + ddnd[cam][1] * rXY[cam][0];
										lazy_eye_data_final[nclust][ExtrinsicAdjustment.get_INDX_DD0(numSensors) + cam] =        ddnd[cam][0];
										lazy_eye_data_final[nclust][ExtrinsicAdjustment.get_INDX_ND0(numSensors) + cam] =        ddnd[cam][1];
									} else {
										lazy_eye_data_final[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 0] = Double.NaN;
										lazy_eye_data_final[nclust][2 * cam + ExtrinsicAdjustment.INDX_X0 + 1] = Double.NaN;
										lazy_eye_data_final[nclust][ExtrinsicAdjustment.get_INDX_DD0(numSensors) + cam] =        Double.NaN;
										lazy_eye_data_final[nclust][ExtrinsicAdjustment.get_INDX_ND0(numSensors) + cam] =        Double.NaN;
									}
								}
							}
						}
					} // end of tile
				}
			};
		}
		startAndJoin(threads);
		if (dbg_img2 != null) {									
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img2,
					clustersX,
					clustersY,
					true,
					"ly_dbg_clouds"); // name+"-CORR2D-D"+clt_parameters.disparity,

			double[][] dbg_img_combo = new double [dbg_img.length][clustersX*clustersY];
			int dbg_w_indx = 16;
			for (int i = 0; i <  dbg_img_combo.length; i++) {
				dbg_img_combo[i] = dbg_img[i].clone();
				for (int j = 0; j < dbg_img[i].length; j++) {
					if (dbg_img2[dbg_w_indx][j] > 0.0) {
						dbg_img_combo[i][j] = dbg_img2[i][j]; 					
					}
				}
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img_combo,
					clustersX,
					clustersY,
					true,
					"ly_dbg_combo");
		}
		return lazy_eye_data_final;
	}	

	
	
	
}
	
