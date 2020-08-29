package com.elphel.imagej.tileprocessor;

import java.awt.Rectangle;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.gpu.GPUTileProcessor;

import Jama.Matrix;

public class ImageDtt extends ImageDttCPU {
	
	private final GPUTileProcessor.GpuQuad gpuQuad;

	public ImageDtt(
			int transform_size,
			boolean mono,
			boolean lwir,
			double scale_strengths,
			GPUTileProcessor.GpuQuad gpuQuadIn){
		super (
				transform_size,
				mono,
				lwir,
				scale_strengths);
		gpuQuad = gpuQuadIn;
	}

	public ImageDtt(
			int transform_size,
			boolean mono,
			boolean lwir,
			double scale_strengths){
		super (
				transform_size,
				mono,
				lwir,
				scale_strengths);
		gpuQuad = null;

	}
	
	public GPUTileProcessor.GpuQuad getGPU() {
		return this.gpuQuad;
	}
	
//	public double [][][][][][]
	public void clt_aberrations_quad_corr_GPU(
			final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
			final int                 macro_scale,     // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
			final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
//// Uses quadCLT from gpuQuad			
////		final double [][][]       image_data,      // first index - number of image in a quad
////		    final boolean [][]        saturation_imp,  // (near) saturated pixels or null
			 // correlation results - final and partial
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
			
			final float [][]          texture_img,     // null or [3][] (RGB) or [4][] RGBA
			final Rectangle           texture_woi,     // null or generated texture location/size
////		final double [][][][]     texture_tiles,   // [tilesY][tilesX]["RGBA".length()][];  null - will skip images combining
			final float [][][]        iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
////		final int                 width,
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
////			final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum (use for absolute)
////			final boolean             corr_sym,
////			final double              corr_offset,
			final double              corr_red, // +used
			final double              corr_blue,// +used
////			final double              corr_sigma,
////			final boolean             corr_normalize,  // normalize correlation results by rms
////	  		final double              min_corr,        // 0.02; // minimal correlation value to consider valid
////			final double              max_corr_sigma,  // 1.2;  // weights of points around global max to find fractional
			final double              max_corr_radius, // 3.9;
////			final boolean 			  max_corr_double, //"Double pass when masking center of mass to reduce preference for integer values
			final int                 corr_mode, // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
			
			final double              min_shot,        // +10.0;  // Do not adjust for shot noise if lower than
			final double              scale_shot,      // +3.0;   // scale when dividing by sqrt ( <0 - disable correction)
			final double              diff_sigma,      // +5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
			final double              diff_threshold,  // +5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
			final boolean             diff_gauss,      // true;  // when averaging images, use Gaussian around average as weight (false - sharp all/nothing)
			final double              min_agree,       // +3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			final boolean             dust_remove,     // +Do not reduce average weight when only one image differs much from the average
////			final boolean             keep_weights,    // Add port weights to RGBA stack (debug feature)
			final GeometryCorrection  geometryCorrection, // for GPU TODO: combine geometry corrections if geometryCorrection_main is not null
			final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
////		using quadCLT instance
////			final double [][][][][][] clt_kernels,    // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
////			final int                 kernel_step,
			final int                 window_type,    // GPU: will not be used
////			final double [][]         shiftXY,        // [port]{shiftX,shiftY} // GPU: will not be used
			final double              disparity_corr, // disparity at infinity
////			final double [][][]       fine_corr, // quadratic coefficients for fine correction (or null) // GPU: will not be used
////			final double              corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
////			final double              shiftX, // shift image horizontally (positive - right) - just for testing // GPU: will not be used
////			final double              shiftY, // shift image vertically (positive - down)                       // GPU: will not be used
			final int                 debug_tileX,
			final int                 debug_tileY,
//			final boolean             no_fract_shift,   // GPU: will not be used
//			final boolean             no_deconvolution, // GPU: will not be used
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

		final double [][] debug_offsets = new double[imgdtt_params.lma_dbg_offset.length][2];
		for (int i = 0; i < debug_offsets.length; i++) for (int j = 0; j < debug_offsets[i].length; j++) {
			debug_offsets[i][j] = imgdtt_params.lma_dbg_offset[i][j]*imgdtt_params.lma_dbg_scale;
		}


		final boolean macro_mode = macro_scale != 1;      // correlate tile data instead of the pixel data
		final int quad = 4;   // number of subcameras
		final int numcol = 3; // number of colors // keep the same, just do not use [0] and [1], [2] - green
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
/*
		double [] scales = isMonochrome() ?
				(new double [] {1.0}) :
					(macro_mode?
							(new double [] {0.25,0.25,0.5}) :
								(new double [] {
										corr_red, // 0.25
										corr_blue, // 0.25
										1.0 - corr_red - corr_blue})); // 0.5
		
*/		
		final int corr_size = transform_size * 2 - 1;
//		final int [][] transpose_indices = new int [corr_size*(corr_size-1)/2][2];
		
		if ((globalDebugLevel > -10) && (disparity_corr != 0.0)){
			System.out.println(String.format("Using manual infinity disparity correction of %8.5f pixels",disparity_corr));
		}
/*		
		{ int indx = 0;
		for (int i =0; i < corr_size-1; i++){
			for (int j = i+1; j < corr_size; j++){
				transpose_indices[indx  ][0] = i * corr_size + j;
				transpose_indices[indx++][1] = j * corr_size + i;
			}
		}
		}
*/

////		final int first_color = isMonochrome()? MONO_CHN : 0; // color that is non-zero

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
		if (globalDebugLevel > 0){
			System.out.println("ortho_height="+ imgdtt_params.ortho_height+" ortho_eff_height="+ imgdtt_params.ortho_eff_height);
			for (int i = 0; i < corr_size; i++){
				System.out.println(" ortho_weights["+i+"]="+ ortho_weights[i]);
			}
		}




		if (globalDebugLevel > 0) {
			System.out.println("clt_aberrations_quad_corr(): width="+width+" height="+height+" transform_size="+transform_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		
	
		// add optional initialization of debug layers here
		boolean need_macro = false;
		boolean need_corr = (clt_mismatch != null);
		// skipping DISPARITY_VARIATIONS_INDEX - it was not used
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
		
		if (clt_mismatch != null){
			for (int i = 0; i<clt_mismatch.length;i++){
				clt_mismatch[i] = new double [tilesY*tilesX]; // will use only "center of mass" centers
			}
		}

		DttRad2 dtt = new DttRad2(transform_size);
		dtt.set_window(window_type);
		final double [] lt_window = dtt.getWin2d();	// [256] - never used
		final double [] lt_window2 = new double [lt_window.length]; // squared - never used

		for (int i = 0; i < lt_window.length; i++) lt_window2[i] = lt_window[i] * lt_window[i];

		if (globalDebugLevel > 1) {
			ShowDoubleFloatArrays sdfa_instance = new ShowDoubleFloatArrays(); // just for debugging?
			sdfa_instance.showArrays(lt_window,  2*transform_size, 2*transform_size, "lt_window");
		}
		if (globalDebugLevel > 0) {
			System.out.println("macro_mode="+macro_mode);
		}

		Matrix [] corr_rots_aux = null;
		if (geometryCorrection_main != null) {
			corr_rots_aux = geometryCorrection.getCorrVector().getRotMatrices(geometryCorrection.getRotMatrix(true));
		}
		
		final boolean use_main = corr_rots_aux != null;
////		final Matrix [] corr_rots = use_main ? corr_rots_aux : geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		boolean [] used_corrs = new boolean[1];
	    final int all_pairs = imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
		final GPUTileProcessor.TpTask[] tp_tasks =  gpuQuad.setTpTask(
				disparity_array, // final double [][]  disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
				disparity_corr,  // final double       disparity_corr,
				used_corrs,      // final boolean []   need_corrs,       // should be initialized to boolean[1] or null
				tile_op,         // final int [][]     tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				all_pairs,       // final int                      corr_mask,        // <0 - use corr mask from the tile tile_op, >=0 - overwrite all with non-zero corr_mask_tp
				threadsMax);     // final int          threadsMax,       // maximal number of threads to launch
		final boolean fneed_macro = need_macro;
		final boolean fneed_corr =  need_corr && used_corrs[0];

		final float [][] lpf_rgb = new float[][] {
			floatGetCltLpfFd(gpu_sigma_r),
			floatGetCltLpfFd(gpu_sigma_b),
			floatGetCltLpfFd(gpu_sigma_g),
			floatGetCltLpfFd(gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				globalDebugLevel > -1);

		final float [] lpf_flat = floatGetCltLpfFd(gpu_sigma_corr);

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				globalDebugLevel > -1);

		final float [] lpf_rb_flat = floatGetCltLpfFd(gpu_sigma_rb_corr);
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				globalDebugLevel > -1);

		gpuQuad.setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_main); // use_aux); // boolean use_aux)

		gpuQuad.execSetTilesOffsets(); // prepare tiles offsets in GPU memory
		gpuQuad.execConvertDirect();
		if (iclt_fimg != null) {
			gpuQuad.execImcltRbgAll(isMonochrome());  // execute GPU kernel
			for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
				iclt_fimg[ncam] = gpuQuad.getRBG(ncam); // retrieve data from GPU
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
					(texture_woi != null)? texture_woi : woi);
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
				false,                 // boolean   calc_textures,
				true);                   // boolean   calc_extra)
			float [][] extra = gpuQuad.getExtra();
			int num_cams = gpuQuad.getNumCams();
			for (int ncam = 0; ncam < num_cams; ncam++) {
				int indx = ncam + IMG_DIFF0_INDEX;
				if ((disparity_modes & (1 << indx)) != 0){
					disparity_map[indx] = new double [extra[ncam].length];
					for (int i = 0; i < extra[ncam].length; i++) {
						disparity_map[indx][i] = extra[ncam][i];
					}
				}
			}
			for (int nc = 0; nc < (extra.length - num_cams); nc++) {
				int sindx = nc + num_cams;
				int indx = nc + IMG_TONE_RGB;
				if ((disparity_modes & (1 << indx)) != 0){
					disparity_map[indx] = new double [extra[sindx].length];
					for (int i = 0; i < extra[sindx].length; i++) {
						disparity_map[indx][i] = extra[sindx][i];
					}
				}
			}			
			
		}
		// does it need correlations?
		if (fneed_corr) {
			//Generate 2D phase correlations from the CLT representation
			gpuQuad.execCorr2D_TD(col_weights); // Get TD version of correlations (may be read out and saved) 
			final int [] corr_indices = gpuQuad.getCorrIndices();
			gpuQuad.execCorr2D_normalize(
	        		false, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D = gpuQuad.getCorr2D(gpu_corr_rad); //  int corr_rad);
			
			// Combine 4 ortho pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        GPUTileProcessor.NUM_PAIRS,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x0f); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			// normalize and convert to pixel domain
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
		    		gpu_corr_rad); // int corr_radius
			final int [] corr_quad_indices = gpuQuad.getCorrComboIndices(); // get quad
			final float [][] fcorr2D_quad =   gpuQuad.getCorr2DCombo(gpu_corr_rad);

			// Combine 2 diagonal pairs			
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        GPUTileProcessor.NUM_PAIRS,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x30); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_cross =   gpuQuad.getCorr2DCombo(gpu_corr_rad);
			
			// Combine 2 horizontal pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        GPUTileProcessor.NUM_PAIRS,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x03); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
		    		gpu_corr_rad); // int corr_radius
			final float [][] fcorr2D_hor =   gpuQuad.getCorr2DCombo(gpu_corr_rad);

			// Combine 2 vertical pairs
			gpuQuad.execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        GPUTileProcessor.NUM_PAIRS,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x0c); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			gpuQuad.execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
					gpu_fat_zero, // double fat_zero);
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
							"dbg-corr2D", // name+"-CORR2D-D"+clt_parameters.disparity,
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
//							double [][]  corrs = new double [GPUTileProcessor.NUM_PAIRS][corr_length]; // 225-long (15x15)
							
							Correlation2d corr2d = new Correlation2d(
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
								// double [][]  corrs = new double [GPUTileProcessor.NUM_PAIRS][corr_length]; // 225-long (15x15)
								// added quad and cross combos
								double [][]  corrs = new double [GPUTileProcessor.NUM_PAIRS + 4][corr_length]; // 225-long (15x15)
								int indx_corr = indx_tile * num_tile_corr;
								int nt = (corr_indices[indx_corr] >> GPUTileProcessor.CORR_NTILE_SHIFT);
								int tileX = nt % tilesX;
								int tileY = nt / tilesX;
								int tIndex = tileY * tilesX + tileX;
								
								// Prepare the same (currently 10-layer) corrs as double [][], as in CPU version
								int pair_mask = 0;
								for (int indx_pair = 0; indx_pair < num_tile_corr; indx_pair++) {
									int pair = corr_indices[indx_corr] & GPUTileProcessor.CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
									assert pair < GPUTileProcessor.NUM_PAIRS : "invalid correllation pair";
									pair_mask |= (1 << pair);
									for (int i = 0; i < corr_length; i++) {
										corrs[pair][i] = gpu_corr_scale * fcorr2D[indx_corr][i]; // from float to double
									}
									indx_corr++; 
								}
								// add 4 combo layers : quad, cross, hor, vert
								int pair = GPUTileProcessor.NUM_PAIRS; // 6
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
								
								// does not include combo
								int used_pairs = pair_mask; // imgdtt_params.dbg_pair_mask; //TODO: use tile tasks
								
								int tile_lma_debug_level =  ((tileX == debug_tileX) && (tileY == debug_tileY))? (imgdtt_params.lma_debug_level-1) : -2;
								boolean debugTile =(tileX == debug_tileX) && (tileY == debug_tileY) && (globalDebugLevel > -1);
								
								// non-GPU initializaqtion of the data structures
								final int [] overexp_all = (saturation_imp != null) ? ( new int [2]): null;
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
								//clt_mismatch should only be used with disparity_map != null;
								if (clt_mismatch != null) {
									for (int np = 0; np < clt_mismatch.length/3; np++) {
										clt_mismatch[3 * np + 0 ][tIndex] = Double.NaN;
										clt_mismatch[3 * np + 1 ][tIndex] = Double.NaN;
										clt_mismatch[3 * np + 2 ][tIndex] = 0;
									}
								}

								// calculate overexposed fraction
								if (saturation_imp != null){
									disparity_map[OVEREXPOSED][tIndex] = (1.0 * overexp_all[0]) / overexp_all[1];
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

						        double [][] strips_intra = new double[2][];
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
						        double [] strip_combo_intra = corr2d.combineInterpolatedCorrelations(
						                strips_intra,                  // double [][] strips,
						                3,                     // int         pairs_mask,
						                imgdtt_params.corr_offset,     // double      offset);
						                imgdtt_params.twice_diagonal); //    		boolean     twice_diagonal)
						        
						        
								// Debug feature - only calculated if requested
								if ((clt_corr_partial != null) && (imgdtt_params.corr_mode_debug || imgdtt_params.gpu_mode_debug)) {
									double [] strip_ortho = corr2d.combineInterpolatedCorrelations(
											strips,                         // double [][] strips,
											0x0f,                           // int         pairs_mask,
											imgdtt_params.corr_offset,      // double      offset);
											imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

									double [] strip_diag = corr2d.combineInterpolatedCorrelations(
											strips,                         // double [][] strips,
											0x30,                           // int         pairs_mask,
											imgdtt_params.corr_offset,      // double      offset);
											imgdtt_params.twice_diagonal);  //    		boolean     twice_diagonal)

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
						            clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strips_intra[0]); // 13
						            clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strips_intra[1]);  // 14
//												    	clt_corr_partial[tileY][tileX][3][0] = corr2d.debugStrip2(strip_ortho); // 13
//												    	clt_corr_partial[tileY][tileX][3][1] = corr2d.debugStrip2(strip_diag);  // 14
						            clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip2(strip_combo_intra);    // 15
//												    	clt_corr_partial[tileY][tileX][3][2] = corr2d.debugStrip(strip_all);    // 15
						            clt_corr_partial[tileY][tileX][3][3] = corr2d.debugStrip2(strip_combo); // 16
								}
						        if (imgdtt_params.pcorr_use) {
						            strip_combo = strip_combo_intra;
						        }
								if ((clt_corr_combo != null) && imgdtt_params.corr_mode_debug) {
									// reuse it too?

								}
								// calculate CM maximums for all mixed channels
								// First get integer correlation center, relative to the center
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
									disparity_map[DISPARITY_INDEX_INT][tIndex] =      -ixy[0];
									//								disparity_map[DISPARITY_INDEX_INT + 1][tIndex] =
									disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] = disp_str[1];
									if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
										System.out.println("BUG: 1. disparity_map[DISPARITY_STRENGTH_INDEX]["+tIndex+"] should not be NaN");
									}
									corr_stat = corr2d.getMaxXCm(   // get fractional center as a "center of mass" inside circle/square from the integer max
											strip_combo,                      // double [] data,      // [data_size * data_size]
											ixy[0],                           // int       ixcenter,  // integer center x
											// corr_wndy,                        // double [] window_y,  // (half) window function in y-direction(perpendicular to disparity: for row0  ==1
											// corr_wndx,                        // double [] window_x,  // half of a window function in x (disparity) direction
											(tile_lma_debug_level > 0)); // boolean   debug);
								}
								if (imgdtt_params.pcorr_use) {
						            // for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
						            double [] hor_pair1 = corr2d.getMaxXSOrtho(
						                    corrs,                              // double [][] correlations,
						                    0x100, //  corrs[8] Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
						                    imgdtt_params.corr_offset,          // double      corr_offset,
						                    true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						                    false,                              // boolean     is_vert,      // transpose X/Y
						                    tile_lma_debug_level > 0);          // boolean   debug);
						            if (hor_pair1 != null) {
						                disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
						                disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
						            }

						            double [] vert_pair1 = corr2d.getMaxXSOrtho(
						                    corrs,                              // double [][] correlations,
						                    0x200, // corrs[9] Correlation2d.getMaskVertical(1), // int         pairs_mask,
						                    imgdtt_params.corr_offset,        // double      corr_offset,
						                    true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
						                    false, // already transposed  // true,                             // boolean     is_vert,      // transpose X/Y
						                    tile_lma_debug_level > 0); // boolean   debug);
						            if (vert_pair1 != null) {
						                disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
						                disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
						            }
								} else {								
									// for compatibility with old code executed unconditionally. TODO: Move to if (corr_stat != null) ... condition below
									double [] hor_pair1 = corr2d.getMaxXSOrtho(
											corrs,                              // double [][] correlations,
											Correlation2d.getMaskHorizontal(1), // int         pairs_mask,
											imgdtt_params.corr_offset,          // double      corr_offset,
											true,                               // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
											false,                              // boolean     is_vert,      // transpose X/Y
											tile_lma_debug_level > 0);          // boolean   debug);
									if (hor_pair1 != null) {
										disparity_map[DISPARITY_INDEX_HOR][tIndex] =          -hor_pair1[0];
										disparity_map[DISPARITY_INDEX_HOR_STRENGTH][tIndex] =  hor_pair1[1];
									}

									double [] vert_pair1 = corr2d.getMaxXSOrtho(
											corrs,                              // double [][] correlations,
											Correlation2d.getMaskVertical(1), // int         pairs_mask,
											imgdtt_params.corr_offset,        // double      corr_offset,
											true,                             // boolean     symmetric,   // for comparing with old implementation average with symmetrical before multiplication
											true,                             // boolean     is_vert,      // transpose X/Y
											tile_lma_debug_level > 0); // boolean   debug);
									if (vert_pair1 != null) {
										disparity_map[DISPARITY_INDEX_VERT][tIndex] =         -vert_pair1[0];
										disparity_map[DISPARITY_INDEX_VERT_STRENGTH][tIndex] = vert_pair1[1];
									}
								}
								// proceed only if CM correlation result is non-null // for compatibility with old code we need it to run regardless of the strength of the normal correlation
								if (corr_stat != null) {
									// skipping DISPARITY_VARIATIONS_INDEX - it was not used
									disp_str[0] = -corr_stat[0];
									disparity_map[DISPARITY_INDEX_CM][tIndex] = disp_str[0]; // disparity is negative X
									if (tile_lma_debug_level > 0) {
										System.out.println("Will run getMaxXSOrtho( ) for tileX="+tileX+", tileY="+tileY);
									}

									// debug new LMA correlations
									if (debugTile) {
										System.out.println("Will run new LMA for tileX="+tileX+", tileY="+tileY+"\n\n\n*** Disabled fro the GPU ***\n\n\n");
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
										if (imgdtt_params.pcorr_use) {
											double [][] fake_corrs = {corrs[6],null,null,null,corrs[7],null};
						                    lma = corr2d.corrLMA(
						                            imgdtt_params,                // ImageDttParameters  imgdtt_params,
						                            fake_corrs,                   // double [][]         corrs,
						                            0x11, // imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
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
													imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
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
											disparity_map[DISPARITY_INDEX_POLY]         [tIndex] = lma_disparity_strength[0];

											// if enabled overwrite - replace  DISPARITY_INDEX_CM and DISPARITY_STRENGTH_INDEX
											if (imgdtt_params.mix_corr_poly) { //true
												disp_str[0] =  lma_disparity_strength[0];
												disp_str[1] =  lma_disparity_strength[1];
												disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
												disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
												if (Double.isNaN(disparity_map[DISPARITY_STRENGTH_INDEX][tIndex])) {
													System.out.println("BUG: 2. disparity_map[DISPARITY_STRENGTH_INDEX][tIndex] should not be NaN");
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
														imgdtt_params.dbg_pair_mask,  // int                 pair_mask, // which pairs to process
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
													System.out.println(String.format("corr4dirsLMA -> ↔: %7.4f (%7.4f) ↕:%7.4f (%7.4f) ⤡:%7.4f (%7.4f) ⤢:%7.4f (%7.4f)",
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
														disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
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
											    disparity_map[DISPARITY_INDEX_CM]       [tIndex] = disp_str[0];
											    disparity_map[DISPARITY_STRENGTH_INDEX] [tIndex] = disp_str[1];
												
											}
										}
									}
								} // end of if (corr_stat != null)
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
		}
		if ((dbg_distort != null) &&(globalDebugLevel >=0)) {
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
	
}
