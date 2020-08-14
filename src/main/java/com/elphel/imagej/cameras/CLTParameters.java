package com.elphel.imagej.cameras;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Properties;
import java.util.Set;

import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.lwir.LwirReaderParameters;
import com.elphel.imagej.tileprocessor.BiQuadParameters;
import com.elphel.imagej.tileprocessor.ImageDtt;
import com.elphel.imagej.tileprocessor.ImageDttParameters;
import com.elphel.imagej.tileprocessor.MeasuredLayersFilterParameters;
import com.elphel.imagej.tileprocessor.PoleProcessorParameters;

import ij.gui.GenericDialog;

public class CLTParameters {
	public int        transform_size =      8; //
	public int        clt_window =          1; // currently only 3 types of windows - 0 (none), 1 and 2
	public double     shift_x =           0.0;
	public double     shift_y =           0.0;
	public int        tileStep =            4;  // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
	public int        iclt_mask =          15;  // which transforms to combine
	public int        tileX =            -258;  // number of kernel tile (0..163)
	public int        tileY =             133;  // number of kernel tile (0..122)
	public int        dbg_mode =            0;  // 0 - normal, +1 - no DCT/IDCT
	public int        ishift_x =            0;  // debug feature - shift source image by this pixels left
	public int        ishift_y =            0;  // debug feature - shift source image by this pixels down

	private double    fat_zero =          0.05; // modify phase correlation to prevent division by very small numbers
	private double    fat_zero_mono =     0.1;  // modify phase correlation to prevent division by very small numbers
	private double    corr_sigma =        0.8;  // LPF correlation sigma
	private double    corr_sigma_mono =   0.1;  // LPF correlation sigma for monochrome images
	private double    scale_strength_main = 1.0; // leave as is
	private double    scale_strength_aux =  0.3; // reduce to match lower sigma
	public boolean    norm_kern =         true; // normalize kernels
	public boolean    gain_equalize =     false;// equalize green channel gain (bug fix for wrong exposure in Exif ?)
	public boolean    colors_equalize =   true; // equalize R/G, B/G of the individual channels
	public boolean    nosat_equalize =    true; // Skip saturated when adjusting gains
	public double     sat_level =         0.95; // Saturation level of the most saturated color channel
	public double     max_overexposure =  0.6;  // Do not use tiles with higher fraction of (near) saturated tiles
	public double     novignetting_r    = 0.2644; // reg gain in the center of sensor calibration R (instead of vignetting)
	public double     novignetting_g    = 0.3733; // green gain in the center of sensor calibration G
	public double     novignetting_b    = 0.2034; // blue gain in the center of sensor calibration B
	public double     scale_r =           1.0; // extra gain correction after vignetting or non-vignetting, before other processing
	public double     scale_g =           1.0;
	public double     scale_b =           1.0;
	public double     vignetting_max    = 0.4; // value in vignetting data to correspond to 1x in the kernel
	public double     vignetting_range  = 5.0; // do not try to correct vignetting less than vignetting_max/vignetting_range
	public int        kernel_step =       16;  // source kernels step in pixels (have 1 kernel margin on each side)
	public double     disparity  =        0.0; // nominal disparity between side of square cameras (pix)
	public double     z_correction  =     0.0; // Inverse distance to infinity (misalignment correction)
	public boolean    correlate =         true; // calculate correlation
	public int        corr_mask =         15;  // bitmask of pairs to combine in the composite
	public boolean    corr_sym =          false; // combine correlation with mirrored around disparity direction
	public boolean    corr_keep =         true;  // keep all partial correlations (otherwise - only combined one)
	public boolean    corr_show =         false; // Show combined correlations
	public boolean    corr_mismatch=      false; // calculate per-pair X/Y variations of measured correlations
	// TODO: what to do if some occlusion is present (only some channels correlate)
	public double     corr_offset =       0.1; //0.1;  // add to pair correlation before multiplying by other pairs (between sum and product)
	// negative - add, not mpy
	public double     corr_red =          0.5;  // Red to green correlation weight
	public double     corr_blue =         0.2;  // Blue to green correlation weight
	public boolean    corr_normalize =    false; // normalize each correlation tile by rms
	public double     min_corr =          0.02; // minimal correlation value to consider valid
	public double     min_corr_normalized =  2.0; // minimal correlation value to consider valid when normalizing correlation results
	public double     max_corr_sigma =    1.2;  // weights of points around global max to find fractional
	// pixel location by quadratic approximation
	public double     corr_r2_offset =    0.0;  // 0.0 - for compatibility with old, normally should be 2.0 (or 1.0)

	public double     max_corr_radius =   3.9;  // maximal distance from int max to consider

	//  		public int        enhortho_width =    2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
	//  		public double     enhortho_scale =    0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)

	public boolean    max_corr_double =   false; // NOT USED double pass when masking center of mass to reduce preference for integer values
	public int        corr_mode =         2;    // which correlation mode to use: 0 - integer max, 1 - center of mass, 2 - polynomial

	// pixel location by quadratic approximation
	public double     corr_border_contrast = 0.01; // contrast of dotted border on correlation results

	public int        tile_task_op =      0xff;   // bitmask of operation modes applied to tiles (0 - nothing), bits TBD later
	// +(0..f) - images, +(00.f0) - process pairs + 256 - force disparity when combining images
	// window to process tiles (later arbitrary masks will be generated to follow particular stages);
	public int        tile_task_wl =      0;   //
	public int        tile_task_wt =      0;   //
	public int        tile_task_ww =      324; //
	public int        tile_task_wh =      242; //
	public double     min_shot =          10.0;  // Do not adjust for shot noise if lower than
	public double     scale_shot =        3.0;   // scale when dividing by sqrt

	public double     diff_sigma =        5.0;   // RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
	public double     diff_threshold =    1.5;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
	public boolean    diff_gauss =        true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
	public double     min_agree =         3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
	public boolean    dust_remove =       true;  // Do not reduce average weight when only one image differes much from the average

	public boolean    black_back =        true;  // use Black for backdrop outside of the FOV
	public boolean    keep_weights =      true;  // add port weights to RGBA stack (debug feature)
	public boolean    sharp_alpha =       false; // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
	public double     alpha0 = 	          0.6; // > .525 Alpha channel 0.0 thereshold (lower - transparent) (watch for far objects)
	public double     alpha1 = 	          0.8; // Alpha channel 1.0 threshold (higher - opaque) (watch for window dust)

	public boolean    gen_chn_stacks =    false; // generate shifted channel rgb stacks
	public boolean    gen_chn_img =       true;  // generate shifted channel images
	public boolean    gen_4_img =         true;  // Generate shifted channel images and save with the model
	public boolean    show_nonoverlap =   true;  // show result RGBA before overlap combined (first channels, then RGBA combined?)
	public boolean    show_overlap =      true;  // show result RGBA (first channels, then RGBA combined?)
	public boolean    show_rgba_color =   true;  // show combined color image
	public boolean    show_map =          true;  // show disparity maps
	public boolean    show_corr =         true;  // show GPU correlation

	public double     disp_scan_start =   0.0;   // disparity scan start value
	public double     disp_scan_step =    1.0;   // disparity scan step
	public int        disp_scan_count =   10;    // disparity scan number of measurements

	public boolean    fine_dbg =          false; // Debug infinity/lazy eye correction
	public double     fine_corr_x_0 =     0.0;   // additionally shift image in port 0 in x direction
	public double     fine_corr_y_0 =     0.0;   // additionally shift image in port 0 in y direction
	public double     fine_corr_x_1 =     0.0;   // additionally shift image in port 1 in x direction
	public double     fine_corr_y_1 =     0.0;   // additionally shift image in port 1 in y direction
	public double     fine_corr_x_2 =     0.0;   // additionally shift image in port 2 in x direction
	public double     fine_corr_y_2 =     0.0;   // additionally shift image in port 2 in y direction
	public double     fine_corr_x_3 =     0.0;   // additionally shift image in port 3 in x direction
	public double     fine_corr_y_3 =     0.0;   // additionally shift image in port 3 in y direction
	public boolean    fine_corr_ignore =  false; // Ignore manual pixel correction
	public boolean    fine_corr_apply =   true;  // Apply and set to ignore manual pixel correction after extrinsics correction

	public double     fcorr_radius =       0.75 ; // Do not try to correct outside this fraction of width/hight
	public double     fcorr_min_strength = 0.15 ; // 0.005 minimal correlation strength to apply fine correction
	public double     fcorr_disp_diff =   1.5;   // consider only tiles with absolute residual disparity lower than
	public boolean    fcorr_quadratic =   true;  // Use quadratic polynomial for fine correction (false - only linear)
	public boolean    fcorr_ignore =      false; // Ignore currently calculated fine correction
	public double     fcorr_inf_strength = 0.20 ; // Minimal correlation strength to use for infinity correction
	public double     fcorr_inf_diff =    0.2;   // Disparity half-range for infinity
	public boolean    fcorr_inf_quad =    true;  // Use quadratic polynomial for infinity correction (false - only linear)
	public boolean    fcorr_inf_vert =    false; // Correct infinity in vertical direction (false - only horizontal)

	//--
	public boolean    inf_disp_apply =   true;   // Apply disparity correction to zero at infinity
	public int        inf_repeat =       5;      // Re run disparity correction at infinity multiple times
	//  		public boolean    inf_mism_apply =   true;   // Apply lazy eye correction at infinity

	public int        inf_iters =        20;     // Infinity extraction - maximum iterations
	public double     inf_final_diff =   0.0001; // Coefficients maximal increment to exit iterations
	public double     inf_far_pull =     0.0;    // include farther tiles than tolerance, but scale their weights

	// infinity filter
	public double     inf_str_pow =      1.0;    // Strength power for infinity filtering
	public int        inf_smpl_side =    3;      // Sample size (side of a square) for infinity filtering
	public int        inf_smpl_num =     5;      // Number after removing worst (should be >1) for infinity filtering
	public double     inf_smpl_rms =     0.1;    // Maximal RMS of the remaining tiles in a sample for infinity filtering

	//Histogram infinity filter
	public int        ih_smpl_step =     8;      // Square sample step (50% overlap)
	public double     ih_disp_min =     -1.0;    // Minimal disparity
	public double     ih_disp_step =     0.05;   // Disparity step
	public int        ih_num_bins =     40;      // Number of bins
	public double     ih_sigma =         0.1;    // Gaussian sigma (in disparity pixels)
	public double     ih_max_diff =      0.1;    // Keep samples within this difference from farthest maximum
	public int        ih_min_samples =  10;      // Minimal number of remaining samples
	public boolean    ih_norm_center =  true;    // Replace samples with a single average with equal weight
	public boolean    inf_restore_disp = true;   // Add disparity back to d{x,y}[i] (debug feature)
	// Lazy eye parameters

	public boolean    ly_lma_ers =      true;    // Use 2020 LMA-based measurement of mismatch
	public double     ly_gt_strength =  0.18;    // use some configurable parameters
	public boolean    ly_gt_use_wnd =   true;    //
	public double     ly_gt_rms =       0.2;     // split small source samples to FG/BG if all aux tile RMS exceeds this value

// Use for LWIR
	public boolean    lylw_inf_en =       true;    // Simultaneously correct disparity at infinity (both poly and extrinsic)
	public boolean    lylw_aztilt_en =    true;    // Adjust azimuths and tilts
	public boolean    lylw_diff_roll_en = true;    // Adjust differential rolls (3 of 4 angles)
	public boolean    lylw_focalLength=   true;    // Correct scales (focal length temperature? variations)
	public boolean    lylw_com_roll=      false;   // Enable common roll (valid for high disparity range only)
	public int        lylw_par_sel   =    0;       // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use checkbox selections above)

	public double     ly_marg_fract =   0.2;     // part of half-width, and half-height to reduce weights
	public boolean    ly_on_scan =      true;    // Calculate and apply lazy eye correction after disparity scan (poly or extrinsic)
	public boolean    ly_inf_en =       true;    // Simultaneously correct disparity at infinity (both poly and extrinsic)
	public int        ly_min_forced =     20;    // Minimal number of clusters with forced disparity to use it
	public boolean    ly_aztilt_en =    true;    // Adjust azimuths and tilts
	public boolean    ly_diff_roll_en = true;    // Adjust differential rolls (3 of 4 angles)
	public boolean    ly_focalLength=   true;    // Correct scales (focal length temperature? variations)
	public boolean    ly_com_roll=      false;   // Enable common roll (valid for high disparity range only)
	public boolean    ly_ers_rot=       true;    // Enable ERS correction of the camera rotation
	public boolean    ly_ers_forw=      true;    // Enable ERS correction of the camera linear movement in z direction
	public boolean    ly_ers_side=      false;   // true;    // Enable ERS correction of the camera linear movement in x direction
	public boolean    ly_ers_vert=      false;   // true;    // Enable ERS correction of the camera linear movement in y direction

	public int        ly_par_sel   =    0;       // Manually select the parameter mask bit 0 - sym0, bit1 - sym1, ... (0 - use checkbox selections above)

	public int        ly_debug_level =  0;       // LY debug level

	public boolean    ly_right_left=    true;    // equalize weights of right/left FoV (use with horizon in both halves and gross infinity correction)

	public int        ly_per_quad =     10;      // minimal tiles per quadrant (not counting the worst) tp proceed
	public double     ly_per_quad_r =   0.003;    // minimal relative tiles per quadrant (not counting the worst) tp proceed
	public int        ly_inf =          10;      // minimal number of tiles at infinity to proceed
	public double     ly_inf_r =        0.0; //0.01;    // minimal relative number of tiles at infinity to proceed
	public int        ly_inf_scale =    20;      // minimal number of tiles at infinity to apply weight scaling
	public double     ly_inf_scale_r =  0.02;    // minimal relative number of tiles at infinity to apply weight scaling

	public double     ly_inf_frac =     0.7;     // Relative weight of infinity calibration data
	public double     ly_inf_max_disparity = 0.2;   // Maximal disparity to be treated as infinity when adjusting with rig data
	public boolean    ly_inf_disp=      false;   // Correct disparity for infinity tiles
	public boolean    ly_inf_force=     false;   // Force convergence correction during extrinsic, even with no infinity data
	public boolean    ly_poly =         false;   // Use polynomial correction, false - correct tilt/azimuth/roll of each sensor

	public int        ly_smpl_side =    3;       // Sample size (side of a square) disp/strength filter
	public int        ly_smpl_num =     5;       // Number after removing worst (should be >1)
	// 		public double     ly_meas_disp =    1.5;     // Maximal measured relative disparity - using  (0.8*disp_scan_step)
	public double     ly_smpl_rms =     0.05;     // 1;     // Maximal RMS of the remaining tiles in a sample

	public double     ly_disp_var =     0.1;     // Maximal full disparity difference to 8 neighbors
	public double     ly_disp_rvar =    0.01;    // Maximal relative full disparity difference to 8 neighbors

	public double     ly_disp_var_gt =  0.5;     // Maximal full disparity difference to 8 neighbors when  GT is available
	public double     ly_disp_rvar_gt = 0.05;    // Maximal relative full disparity difference to 8 neighbors when  GT is available


	public double     ly_norm_disp =    5.0;     // Reduce weight of higher disparity tiles

	// Lazy eye multi-step fitting
	public double     lym_overexp     = 0.0001;  // Any (near) saturated pixels - discard tile (see sat_level also)
	public boolean    lym_update_disp = true;    // Update target disparity after each step
	public int        lym_iter =        25;      // Maximal number of iterations
	private double    lym_change =      0.5e-5;    // Parameter vector difference to exit 4e-6 - OK
	private double    lym_change_aux =  1e-4;    // same for aux camera (currntly)lwir
	public double     lym_poly_change = 0.002;   // Parameter vector difference to exit from polynomial correction

	public boolean    lyf_filter =      false;   // Filter lazy eye pairs by their values
	public int        lyf_smpl_side =   3;       // 8 x8 masked, 16x16 sampled
	public double     lyf_rms_max =     0.1;     // Maximal RMS (all components to components average)
	public double     lyf_frac_keep =   0.5;     // Keep best fit samples, discard worst
	public int        lyf_min_samples = 5;       // Minimal number of tiles remaining in the sample
	public boolean    lyf_norm_center = true;    // Replace samples with a single average with equal weight
	public double     ly_corr_scale =   1.0;     // Scale calculated correction vector
	public boolean    lyr_filter_ds =   false;    // true;
	public boolean    lyr_filter_lyf =  false;   // ~clt_parameters.lyf_filter, but may be different, now off for a single cameras


	// old fcorr parameters, reuse?
	// 		public int        fcorr_sample_size = 32;    // Use square this size side to detect outliers
	// 		public int        fcorr_mintiles =    8;     // Keep tiles only if there are more in each square
	// 		public double     fcorr_reloutliers = 0.5;   // Remove this fraction of tiles from each sample
	// 		public double     fcorr_sigma =       20.0;  // Gaussian blur channel mismatch data

	public double     corr_magic_scale =  0.85;  // reported correlation offset vs. actual one (not yet understood)
	public int        corr_select =       0;     // Select correlation type for all-pair supertiles and assignments (0 - CM, 1 - polynomial)

	// 3d reconstruction
	public boolean    show_textures    = true;  // show generated textures
	public boolean    debug_filters    = false;// show intermediate results of filtering
	// not used anywhere so far
	public double     min_smth         = 0.25;  // 0.25 minimal noise-normalized pixel difference in a channel to suspect something
	public double     sure_smth        = 2.0;   // reliable noise-normalized pixel difference in a channel to have something
	public double     bgnd_range       = 0.3;   // disparity range to be considered background
	public double     other_range      = 2.0;   // disparity difference from center (provided) disparity to trust

	public double     ex_strength      = 0.18;  // minimal 4-corr strength to trust tile
	public double     ex_nstrength     = 0.4;   // minimal 4-corr strength divided by channel diff for new (border) tiles

	public boolean    ex_over_bgnd     = false; // Allow expansion over previously identified background (infinity)
	public double     ex_min_over      = 1.0;   // When expanding over background, disregard lower disparity

	public double     pt_super_trust   = 1.6;   // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
	public boolean    pt_keep_raw_fg   = true;  // Do not replace raw tiles by the plates, if raw is closer (like poles)
	public double     pt_scale_pre     = 1.5;   // Scale plates strength before comparing to raw strength
	public double     pt_scale_post    = 2.5;   // Scale plates strength when replacing raw (plates d/s data is more reliable if it exists)

	public double     bgnd_sure        = 0.18;  // minimal strength to be considered definitely background
	public double     bgnd_maybe       = 0.1; // maximal strength to ignore as non-background
	//  		public double     bgnd_2diff       = 0.005; // maximal strength to ignore as non-background
	public int        min_clstr_seed   = 4; //2;     // number of tiles in a cluster to seed (just background?)
	public int        min_clstr_lone   = 4;     // number of tiles in a cluster not close to other clusters (more than 2 tiles apart)
	public double     min_clstr_weight = 0.0;   // Minimal total strength of the cluster
	public double     min_clstr_max    = 0.25;  // Minimal maximal strength of the cluster

	public int        fill_gaps        = 4;     // same as in grow - 1:  4 directions by 1 step, 2: 8 directions by 1 step. +2*n - alternating hor/vert
	public int        fill_final       = 50;    // same as fill_gaps, on the final pass
	public int        min_clstr_block  = 3;     // number of tiles in a cluster to block (just non-background?)
	public int        bgnd_grow        = 2;     // number of tiles to grow (1 - hor/vert, 2 - hor/vert/diagonal)

	//  		public double     ortho_min        = 0.09;  // minimal strength of hor/vert correlation to be used instead of full 4-pair correlation
	public boolean    ortho_old        = false; // use old ortho features processing (orth0_* parameters, false - use or_*)
	public double     ortho_min_hor    = 0.07;  // minimal strength of hor correlation to be used instead of full 4-pair correlation -
	public double     ortho_min_vert   = 0.15;  // minimal strength of vert correlation to be used instead of full 4-pair correlation
	public double     ortho_asym       = 1.2;   // vert/hor (or hor/vert) strength to be used instead of the full correlation
	public double     ortho_over4      = 0.8;   // vert/hor (or hor/vert) strength exceeding scaled 4-pair strength
	public double     ortho_sustain    = 0.05;  // minimal strength of hor/vert to bridge over
	public int        ortho_run        = 3;     // minimal run of hor/vert tiles to be considered (at least from one side)
	public double     ortho_minmax     = 0.09;  // minimal maximal strength in an ortho run
	public int        ortho_bridge     = 10;    // number of tiles to bridge over hor/vert gaps
	public double     ortho_rms        = 0.3;   // maximal disparity RMS in a run to replace by average
	public int        ortho_half_length = 4;    // convolve hor/vert strength by 3*(2*l+1) kernels to detect multi-tile features
	public double     ortho_mix        = 0.5;   // Fraction ovf convolved ortho in a mix with raw

	// Alternative mixing of ortho disparity/strength
	public boolean    or_hor           = true;  // Apply ortho correction to horizontal correlation (vertical features)
	public boolean    or_vert          = true;  // Apply ortho correction to vertical correlation (horizontal features)
	public double     or_sigma         = 2.0;   // Blur sigma: verically for horizontal correlation, horizontally - for vertically
	public double     or_sharp         = 0.0;   // 0.5;   // 3-point sharpening (-k, +2k+1, -k)
	public double     or_scale         = 2.5;   // Scale ortho correletion strength relative to 4-directional one
	public double     or_offset        = 0.1;   // Subtract from scaled correlation strength, limit by 0
	public double     or_asym          = 1.5;   // Minimal ratio of orthogonal strengths required for dis[parity replacement
	public double     or_threshold     = 0.3;   // 1.5;   // Minimal scaled offset ortho strength to normal strength needed for replacement
	public double     or_absHor        = 0.15;  // Minimal horizontal absolute scaled offset ortho strength needed for replacement
	public double     or_absVert       = 0.19;  // Minimal vertical absolute scaled offset ortho strength needed for replacement
	public double     or_maxDisp       = 5.0;   // Maximal disparity to apply ortho correction

	public boolean    poles_fix        = true;  // Continue vertical structures to the ground
	public int        poles_len        = 25;    // Number of tiles to extend over the poles bottoms
	public double     poles_ratio      = 1.0;   // Maximal ratio of invisible to visible pole length
	public double     poles_min_strength = 0.1; // Set new pole segment strength to max of horizontal correlation and this value
	public boolean    poles_force_disp = true;  // Set disparity to that of the bottom of existing segment (false - use hor. disparity)

	public int        max_clusters     = 500;   // Maximal number of clusters to generate for one run
	public boolean    remove_scans     = true;  // Remove all unneeded scans when generating x3d output to save memory
	public boolean    output_x3d       = true;  // Generate x3d output
	public boolean    output_obj       = true;  // Generate Wavefront obj output


	public boolean    correct_distortions = false; // Correct lens geometric distortions in a model (will need backdrop to be corrected too)
	public boolean    show_triangles =    true;  // Show generated triangles
	public boolean    avg_cluster_disp =  false;  // Weight-average disparity for the whole cluster
	public double     maxDispTriangle   = 0.2;    // Maximal relative disparity difference in a triangle face
	public double     infinityDistance  = 10000;  // Distance to generate backdrop (0 - use regular backdrop)
	public int        min_bgnd_tiles    = 10;     // Minimal number of background tiles to generate background
	public boolean    shUseFlaps        = true;  // Split into shells with flaps
	public boolean    shAggrFade        = false; // true;  // Aggressive fade alpha (whole boundary)
	public int        shMinArea         = 1;     // Minimal shell area (not counting flaps
	public double     shMinStrength    = 0.2;   // Minimal value of the shell maximum strength

	// Thin ice parameters
	public double     tiRigidVertical   = 3.0;   // 2.0 relative disparity rigidity in vertical direction
	public double     tiRigidHorizontal = 1.0;   // 2.0 relative disparity rigidity in horizontal direction
	public double     tiRigidDiagonal   = 0.5;   // 0.5 relative disparity rigidity in diagonal direction
	public double     tiStrengthOffset  = 0.1;   // 0.1 strength "floor" - subtract (limit 0) before applying
	public double     tiDispScale       = 0.5;   // 0.5 divide actual (disparity*eff_Strength) by this  before Math.pow and applying to T.I.
	public double     tiDispPow         = 0.0;   // 0.0 apply pow to disparity (restore sign) for disparity difference pressure on ice
	public double     tiDispPull        =  .1;   // 10.0 tiDispPull: multiply strength*disparity difference to pull force
	public double     tiDispPullPreFinal=  .1;   // 5.0 Scale tiDispPull for pre-final pass
	public double     tiDispPullFinal   =  .01;  // 2.0 Scale tiDispPull for final pass

	public double     tiBreakNorm       =   .5;  // Normalize stresses to average disparity if it is above threshold
	public double     tiBreak3          = 0.6;   // 0.5 TI break value of abs(d0-3d1+3d2-d3)
	public double     tiBreak31         = 0.1;   // 0.04 TI break value of (d0-3d1+3d2-d3) * (d1 - d2)
	public double     tiBreak21         = 0.1;   // 0.1 TI break value of (-d0+d1+d2-d3) * abs(d1 - d2)
	public double     tiBreakFar        = 0.3;   // 0.3  TI disparity threshold to remove as too far tiles
	public double     tiBreakNear       = 0.3;   // 0.3 TI disparity threshold to remove as too near tiles
	public int        tiBreakMode       = 0;     // 1 TI break mode: +1: abs(3-rd derivative), +2: -(3-rd * 1-st), +4: -(2-nd * abs(1-st)) , +8 - remove far, +16 - remove near
	public double     tiBreakSame       = 0.5;   // 0.75 Amplify colinear breaks in neighbor tiles
	public double     tiBreakTurn       = 0.125; // 0.125 Amplify 90-degree turnintg breaks in neighbor tiles

	public double     tiHealPreLast     = 0.1;    // 0.1 Heal disparity gap before pre-last smooth
	public double     tiHealLast        = 0.05;   // 0.05 Heal disparity gap before last smooth
	public int        tiHealSame        = 5;      // 10 Maximal length of an internal break in the cluster to heal

	public int        tiIterations      = 1000;  // 300 maximal number of TI iterations for each step
	public int        tiPrecision       = 6;     // 6 iteration maximal error (1/power of 10)
	public int        tiNumCycles       = 5;     // 5 Number of cycles break-smooth (after the first smooth)

	// FG/BG separation
	public boolean    stUseRefine =       false; // Apply super-tiles during refine passes
	public boolean    stUsePass2 =        true;  // Apply super-tiles during pass2
	public boolean    stUseRender =       true;  // Apply super-tiles during render

	public boolean    stShow =            false; // Calculate and show supertiles histograms
	public int        stSize            = 8;     // Super tile size (square, in tiles)
	public double     stStepFar         = 0.1;   // Disparity histogram step for far objects
	public double     stStepNear        = 0.5;   // Disparity histogram step for near objects
	public double     stStepThreshold   = 1.0;   // Disparity threshold to switch from linear to logarithmic steps
	public double     stMinDisparity    = 0.0;   // Minimal disparity (center of a bin)

	// moved to MeasuredLayersFilterParameters
	//  		public double     stFloor           = 0.1;  // Subtract from strength, discard negative
	//  		public double     stPow             = 1.0;   // raise strength to this power
	public double     stSigma           = 1.5;   // Blur disparity histogram (sigma in bins)
	public double     stMinBgDisparity  = 0.0;   // Minimal backgroubnd disparity to extract as a maximum from the supertiles
	public double     stMinBgFract      = 0.1;   // Minimal fraction of the disparity histogram to use as background
	public double     stUseDisp         = 0.12;  // Use background disparity from supertiles if tile strength is less
	public double     stStrengthScale   = 50.0;  // Multiply st strength if used instead of regular strength
	public boolean    stSmplMode        = true;   // Use sample mode (false - regular tile mode)

	//  		public int        stSmplSide        = 3;      // Sample size (side of a square)
	//  		public int        stSmplNum         = 5;      // Number after removing worst
	//  		public double     stSmplRms         = 0.3;    // Maximal RMS of the remaining tiles in a sample
	//		public boolean    stSmplWnd         = false;  // Use window function for the samples (TODO: change default to true after testing)

	// Additional parameters to prevent increase of the  distance to far (small) objects. It happens because of averaging bg and fg disparities
	// These parameters are also used when expanding tiles from the known ones.

	//		public double     fs_max_abs_tilt  = 2.0; // Maximal absolute tilt in pixels/tile
	//		public double     fs_max_rel_tilt  = 0.2; // Maximal relative tilt in pixels/tile/disparity
	//		public double     fs_damp_tilt  =    0.001; // Tilt cost for damping insufficient plane data
	//		public double     fs_min_tilt_disp = 4.0; // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
	//		public double     fs_transition    = 1.0; // Mode transition range (between tilted and maximal disparity)
	//		public int        fs_far_mode  =     1;   // Far objects filtering mode (0 - off, 1 - power of disparity)
	//		public double     fs_far_power =     3.0; // Raise disparity to this power before averaging for far objects


	public int        stGrowSel         = 2;     // Grow initial selection before processing supertiles, odd - ortho. <0 - use all tiles
	public int        stMeasSel         = 1;     // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert

	public double     stSmallDiff       = 0.4;   // Consider merging initial planes if disparity difference below
	public double     stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above


	public double     outlierStrength  =  0.3;   // Outlier tiles weaker than this may be replaced from neighbors
	public double     outlierDiff      =  0.4;   // Replace weak outlier tiles that do not have neighbors within this disparity difference
	public double     outlierDiffPos   =  1.0;   // Replace weak outlier tiles that have higher disparity than weighted average
	public double     outlierDiffNeg   =  0.4;   // Replace weak outlier tiles that have lower disparity than weighted average

	// TODO: Make refine skip if already good?
	public boolean    combine_refine    = true; // combine with all previous after refine pass
	public double     combine_min_strength = 0.12; // Disregard weaker tiles when combining scans
	public double     combine_min_hor =      0.12; // Disregard weaker tiles when combining scans for horizontal correlation
	public double     combine_min_vert =     0.12; // Disregard weaker tiles when combining scans for vertical correlation

	// TODO: Move together with similar parameters
	//  		public double     unique_tolerance = 0.1; // Do not re-measure correlation if target disparity differs from some previous by this

	// Multi-pass growing disparity
	public int        grow_sweep         = 8; // Try these number of tiles around known ones
	public double     grow_disp_max =   160.0;// Maximal disparity to try
	public double     grow_disp_trust =  4.0; // Trust measured disparity within +/- this value
	public double     grow_disp_step =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?
	public double     grow_min_diff =    0.5; // Grow more only if at least one channel has higher variance from others for the tile

	public boolean    grow_retry_far =  false; // Retry tiles around known foreground that have low max_tried_disparity
	public boolean    grow_pedantic =   false; // Scan full range between max_tried_disparity of the background and known foreground
	public boolean    grow_retry_inf =  false;  // Retry border tiles that were identified as infinity earlier

	// New for initial growing
	public boolean    gr_new_expand        =   true;
	//		public int        gr_max_expand        = 500; // 150; // 30;

	public double     fds_str_floor        =   0.09; // Should be normally less than combine_min_strength
	public double     gr_ovrbg_cmb         =   0.3; // 0.3;
	public double     gr_ovrbg_cmb_hor     =   0.3; // 0.3;
	public double     gr_ovrbg_cmb_vert    =   0.3; // 0.3;
	public double     gr_ovrbg_filtered    =   0.3; // 0.3;
	public double     fds_str_pow          =   1.0;
	public int        fds_smpl_side        =   5; // 3;      // Sample size (side of a square)
	public int        fds_smpl_num         =  13; // 13; // 5;      // Number after removing worst (should be >1)
	public double     fds_smpl_rms         =   0.15; // Maximal RMS of the remaining tiles in a sample
	public double     fds_smpl_rel_rms     =   0.01; // 05;  // Maximal RMS/disparity in addition to smplRms
	public boolean    fds_smpl_wnd         =   true; //
	public double     fds_abs_tilt         =   2.0; // pix per tile
	public double     fds_rel_tilt         =   0.2; // (pix / disparity) per tile

	public boolean    per_filter =             true; // detect and filter periodic structures
	public double     per_trustedCorrelation = 2.0;
	public double     per_initial_diff =       0.8;   // initial disparity difference to merge to maximum
	public double     per_strength_floor =     0.1;   //
	public double     per_strength_max_over =  0.03;  //
	public double     per_min_period =         4.0;   //
	public int        per_min_num_periods =    3; // minimal number of periods
	public double     per_disp_tolerance =     1.0;   // maximal difference between the average of fundamental and
	public double     per_disp_match =         1.5;   // disparity difference to match neighbors
	public double     per_strong_match_inc =   0.05;  // extra strength to treat match as strong (for hysteresis)

	// Macro disparity scanning parameters
	public double     mc_disp8_step        =   2.0;   // Macro disparity scan step (actual disparity step is 8x)
	public double     mc_disp8_trust       =   2.0;   //Trust measured disparity within +/- this value
	public double     mc_strength          =   0.2;   // Minimal composite correlation to process (0.2..0.3)
	public double     mc_unique_tol        =   0.05;  // Do not re-measure macro correlation if target disparity differs from some previous by this

	public double     mc_trust_fin         =   0.3;   // When consolidating macro results, exclude high residual disparity
	public double     mc_trust_sigma       =   0.2;   // Gaussian sigma to reduce weight of large residual disparity
	public double     mc_ortho_weight      =   0.5;   // Weight from ortho neighbor supertiles
	public double     mc_diag_weight       =   0.25;  // Weight from diagonal neighbor supertiles
	public double     mc_gap               =   0.4;   // Do not remove measurements farther from the kept ones


	public double     mc_weight_var        =   1.0;   // weight of variance data (old, detects thin wires?)
	public double     mc_weight_Y          =   1.0;   // weight of average intensity
	public double     mc_weight_RBmG       =   5.0;   // weight of average color difference (0.5*(R+B)-G), shoukld be ~5*weight_Y

	//		  0x1e, // 0x1f, // final int         variants_mask,
	public int        gr_min_new           =  20;    // Discard variant if it requests too few tiles
	public boolean    gr_var_new_sngl      =   false;// Expand only unambiguous tiles over previously undefined
	public boolean    gr_var_new_fg        =   true; // Expand unambiguous and foreground tiles over previously undefined
	public boolean    gr_var_all_fg        =   true;
	public boolean    gr_var_new_bg        =   false;
	public boolean    gr_var_all_bg        =   false;
	public boolean    gr_var_next          =   false; // try next disparity range  TODO: add related statements
	public int        gr_num_steps         =   8;    // How far to extend over previously undefined disparity tiles
	public int        gr_steps_over        =   4;    // How far to extend over previously determined disparity tiles
	public int        gr_smpl_size         =   5;    // Extend sample square side
	public int        gr_min_pnts          =   3;    // Extend at least this number of the seed tiles
	public boolean    gr_use_wnd           =   true; // Use window function for square sample
	public double     gr_tilt_damp         =   0.001;// Tilt cost for damping insufficient plane data
	public double     gr_split_rng         =   5.0;  // When growing, range of disparities to be extended without far/near division
	public double     gr_same_rng          =   3.0;  // consider far/near tiles within that range from the farthest/closest
	public double     gr_diff_cont         =   2.0;  // Maximal difference from the old value when smoothing
	public double     gr_abs_tilt          =   2.0;  // Maximal filter disparity absolute tilt (pix per tile)
	public double     gr_rel_tilt          =   0.2;  // Maximal filter disparity tilt (pix / disparity) per tile
	public int        gr_smooth            =   50;   // Maximal number of smoothing steps (reduce if long?)
	public double     gr_fin_diff          =   0.01; // Maximal change to finish smoothing iterations
	public double     gr_unique_tol        =   0.15; // Do not re-measure correlation if target disparity differs from some previous by this
	public double     gr_unique_pretol     =   0.5;  // Larger tolerance for expanding (not refining)


	public boolean    ft_mod_strength =          true;    // When set, multiply each tile strength by the number of selected neighbors
	public boolean    ft_clusterize_by_highest = true;    // Clusterize using disparity horizontal maximums for fronto planes and minimums - for horizontal. False - use histograms
	public double     ft_clust_sigma =             0.7;   // Blur disparity before argmax/argmin for initial clusterization
	public double     ft_disp_arange_vert =        0.07;  // Absolute disparity range for fronto clusters
	public double     ft_disp_rrange_vert =        0.01;  // Relative disparity range for fronto clusters
	public double     ft_disp_arange_hor =         0.035; // Absolute disparity range for horizontal clusters
	public double     ft_disp_rrange_hor =          0.005; // Relative disparity range for horizontal clusters
	public double     ft_tolerance_above_near =  100.0;   // Actual disparity positive tolerance over blurred disparity argmax range
	public double     ft_tolerance_below_near =   -0.01;  // Actual disparity negative tolerance under blurred disparity argmax range
	public double     ft_tolerance_above_far =     0.07;  // Actual disparity positive tolerance over blurred disparity argmin range
	public double     ft_tolerance_below_far =     0.1;   // Actual disparity negative tolerance under blurred disparity argmin range
	public int        ft_hor_vert_overlap =        2;     // Allow clusters tile sharing between fronto and horizontal. 2 - 1 tile in 8 directions, 1 - 1 tile in 4 directions
	public int        ft_used_companions =         5;     // Cell that has this many new used companions is considered used (borders and already use3d are considered used too)
	public int        ft_used_true_companions =    1;     // There should be at least this many new selected tiles among neighbors.,

	public boolean    plPreferDisparity    =   false;// Always start with disparity-most axis (false - lowest eigenvalue)
	public double     plDispNorm           =   5.0;  // Normalize disparities to the average if above (now only for eigenvalue comparison)
	public double     plFrontoTol          =   0.0;  // for compatibility with old //0.1;  // Fronto tolerance (pix) - treat almost fronto as fronto (constant disparity). <= 0 - disable
	public double     plFrontoRms          =   0.05; // Target rms for the fronto planes - same as sqrt(plMaxEigen) for other planes
	public double     plFrontoOffs         =   0.2;  // increasing weight of the near tiles by using difference between the reduced average as weight. <= 0 - disable
	public double     PlFrontoPow          =   1.0;  // increase weight even more

	public double     plBlurBinVert        =   1.2;  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
	public double     plBlurBinHor         =   0.8;  // Blur disparity histograms for horizontal clusters by this sigma (in bins)
	public double     plMaxDiffVert        =   0.4;  // Maximal normalized disparity difference when initially assigning to vertical plane
	public double     plMaxDiffHor         =   0.2;  // Maximal normalized disparity difference when initially assigning to horizontal plane
	public int        plInitPasses         =     3;  // Number of initial passes to assign tiles to vert (const disparity) and hor planes

	public int        plMinPoints          =     5;  // Minimal number of points for plane detection
	public double     plTargetEigen        =   0.02; // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
	public double     plFractOutliers      =   0.3;  // Maximal fraction of outliers to remove
	public int        plMaxOutliers        =   200;  // Maximal number of outliers to remove
	public double     plMinStrength        =   0.01; // Minimal total strength of a plane
	public double     plMaxEigen           =   0.06; // Maximal eigenvalue of a plane
	public double     plEigenFloor         =   0.005;// Add to eigenvalues of each participating plane and result to validate connections
	public double     plEigenStick         =   25.0; // Consider plane to be a "stick" if second eigenvalue is below
	public double     plBadPlate           =   0.2;  // Not a plate if sin^2 between normals from disparity and world exceeds this
	public boolean    plDbgMerge           =   true; // Combine 'other' plane with current
	public double     plWorstWorsening     =   2.0;  // Worst case worsening after merge
	public double     plWorstWorsening2    =   5.0;  // Worst case worsening for thin planes
	public double     plWorstEq            =   1.0;  // Worst case worsening after merge with equal weights
	public double     plWorstEq2           =   2.0;  // Worst case worsening for thin planes with equal weights
	public double     plOKMergeEigen       =   0.03; // If result of the merged planes is below, OK to use thin planes (higher) threshold
	public double     plMaxWorldSin2       =   0.1;  // Maximal sine squared of the world angle between planes to merge. Set to >= 1.0 to disable
	public double     pl2dForSin           =   2.5;  // Do not compare sin() between planes, if at least one has too small axis ratio
	public double     plWeakWorsening      =   1.0;  // Relax merge requirements for weaker planes
	public double     plMaxOverlap         =   0.1;  // Maximal overlap between the same supertile planes to merge
	// Merge same supetile planes if at least one is weak and they do not differ much
	public double     plWeakWeight         =   0.2 ; // Maximal weight of the weak plane to merge
	public double     plWeakEigen          =   0.1;  // Maximal eigenvalue of the result of non-weighted merge
	public double     plWeakWeight2        =  10.0 ; // Maximal weight of the weak plane to merge (second variant)
	public double     plWeakEigen2         =   0.05; // Maximal eigenvalue of the result of non-weighted merge  (second variant)
	public double     plSumThick           =   1.6;  // Do not merge if any sqrt of merged eigenvalue exceeds scaled sum of components
	public double     plNeNeibCost         =   5.0;  // When calculating non-exclusive planes, do not use neighbors with high cost
	public double     plNeOwn              =   5.0;  // When calculating non-exclusive planes, use center plane relative weight

	public double     plExNeibCost         =   5.0;  // When calculating exclusive planes links, do not use neighbors with high cost
	//  	    public double     plExNeibCostSngl     =  10.0;  // When calculating exclusive planes links, do not use no-link neighbors with high cost
	public double     plExNeibSmooth       =   0.5;  // Scale down maximal costs for smoothed planes (tighter requirements)
	public double     plMergeCostStar      =   5.0;  // Cost threshold for merging same tile planes if the plane has connected neighbors
	public double     plMergeCost          =  10.0;  // Cost threshold for merging same tile planes if not connected

	public boolean    plConflMerge         =  true;  // Try to merge conflicting planes
	public double     plConflRelax         =   1.5;  // Scale parameters to relax planes fit for merging conflicting planes
	public boolean    plConflSngl          =  true;  // Only merge conflicting planes if this is the only conflicting pair in the supertile
	public boolean    plConflSnglPair      =  true;  // Only merge conflicting planes only if there are just two planes in the supertile

	public double     plWeakFgStrength     =   0.15; // Consider merging plane if it is foreground and maximal strength below this
	public int        plWeakFgOutliers     =   1;    // Remove these strongest from foreground when determining the maximal strength
	public double     plWeakFgRelax        =   2.0;  // Relax cost requirements when merging with weak foreground


	public double     plThickWorld         =   0.2;  // Maximal real-world thickness of merged overlapping planes (meters)
	public double     plThickWorldConfl    =   0.4;  // Maximal real-world merged thickness for conflicting planes
	public double     plRelaxComplete      =   1.5;  // Relax cost requirements when adding exclusive links to complete squares and triangles
	public double     plRelaxComplete2     =   2.0;  // Relax cost requirements during the second pass


	public double     plMaxZRatio          =   2.0;  // Maximal ratio of Z to allow plane merging
	public double     plMaxDisp            =   0.6;  // Maximal disparity of one of the planes to apply  maximal ratio
	public double     plCutTail            =   1.4;  // When merging with neighbors cut the tail that is worse than scaled best
	public double     plMinTail            =   0.015;// Set cutoff value level not less than


	// parameters to recreate planes from tiles disparity/strengths using determined plane connections to neighbors
	public boolean    plDiscrEn            =   true; // Enable planes tiles selection regeneration hinted by supertile neighbors
	public double     plDiscrTolerance     =   0.4;  // Maximal disparity difference from the plane to consider tile
	public double     plDiscrDispRange     =   1.0;  // Parallel move known planes around original know value for the best overall fit
	public int        plDiscrSteps         =   10;   // Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
	//  		public int        plDiscrVariants      =   500;  // total number of variants to try (protect from too many planes)
	public int        plDiscrMode          =   3;    // 0 - weighted, 1 - equalized, 2 - best, 3 - combined

	public double     plDiscrVarFloor      =   0.03;  // Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
	public double     plDiscrSigma         =   0.05;  // Gaussian sigma to compare how measured data is attracted to planes
	public double     plDiscrBlur          =   0.05;  // Sigma to blur histograms while re-discriminating
	public double     plDiscrExclusivity   =   1.5;   // Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
	public double     plDiscrExclus2       =   0.8;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors
	public boolean    plDiscrStrict        =   false; // When growing selection do not allow any offenders around (false - more these than others)
	public double     plDiscrCorrMax       =   0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
	public double     plDiscrCorrMerge     =   0.85;  // Attraction to different planes correlation that is high enough to merge planes
	public int        plDiscrSteal         =   4;     // If offender has this number of tiles (including center) the cell can not be used
	public int        plDiscrGrown         =   4;     // Only use tiles within this range from original selection

	public double     plDiscrXMedian       =   1.5;   // Remove outliers from the final selection that have distance more than scaled median

	// comparing merge quality for plane pairs
	public double     plCostDist           =   4.0;  // Disparity (pix) - closer cost will use more of the real world, farther - disparity
	public double     plCostKrq            =   0.8;  // Cost of merge quality sqrt(weighted*equal) in disparity space
	public double     plCostKrqEq          =   0.2;  // Cost of merge quality average of weighted and equal weight in disparity space
	public double     plCostWrq            =   0.8;  // Cost of merge quality sqrt(weighted*equal) in world space
	public double     plCostWrqEq          =   0.2;  // Cost of merge quality average of weighted and equal weight in world space
	public double     plCostSin2           =  10.0;  // Cost of sin squared between normals
	public double     plCostRdist2         =1000.0;  // Cost of squared relative distances


	public boolean    plConflDualTri       =   false; // Resolve dual triangles conflict (odoodo)
	public boolean    plConflMulti         =   false; // Resolve multiple odo triangles conflicts
	public boolean    plConflDiag          =   false; // Resolve diagonal (ood) conflicts
	public boolean    plConflStar          =   true;  // Resolve all conflicts around a supertile

	public int        plStarSteps          =   2;    // How far to look around when calculating connection cost
	public double     plStarOrtho          =   0.5;  // When calculating cost for the connections scale 4 ortho neighbors
	public double     plStarDiag           =   0.25; // When calculating cost for the connections scale 4 diagonal neighbors
	public double     plStarPwr            =   0.5;  // Divide cost by number of connections to this power
	public double     plStarWeightPwr      =   0.5;  // use this power of tile weight when calculating connection cost
	public double     plWeightToDens       =   0.3;  // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
	public double     plStarValPwr         =   1.0;  // Raise value of each tile before averaging
	public double     plDblTriLoss         =   0.0001; // When resolving double triangles allow minor degradation (0.0 - strict)
	public boolean    plNewConfl           =   false; // Allow more conflicts if overall cost is reduced
	public int        plMaxChanges         =   0;     // Maximal number of simultaneous connection changes around one tile (0 - any)

	public boolean    plMutualOnly         =   true; // keep only mutual links, remove weakest if conflict
	public boolean    plFillSquares        =   true; // Add diagonals to full squares
	public boolean    plCutCorners         =   true; // Add ortho to 45-degree corners
	public boolean    plHypotenuse         =   true; // Add hypotenuse connection if both legs exist

	public double     plPull               =  5.0; // .3;   // Relative weight of original (measured) plane compared to average neighbor pull
	// when combing with neighbors
	public int        plIterations         =  10;   // Maximal number of smoothing iterations for each step
	public boolean    plStopBad            =  true; // Do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
	public int        plPrecision          =  6;    // Maximal step difference (1/power of 10)

	public double     plSplitPull          =  .5;   // Relative weight of center plane when splitting into pairs
	public double     plNormPow            =  .5;   // 0.0: 8 neighbors pull 8 times as 1, 1.0 - same as 1
	public int        plSplitMinNeib       =  2;    // Minimal number of neighbors to split plane in pairs
	public double     plSplitMinWeight     =  2.0;  // Minimal weight of split plains to show
	public double     plSplitMinQuality    =  1.1;  // Maximal normalized disparity difference from the plane to consider

	public boolean    plSplitApply         = true;  // Apply plane split to pairs
	public boolean    plNonExclusive       = true;  // Allow tiles to belong to both planes of the pair
	public boolean    plUseOtherPlanes     = false; // Allow other tiles from the same supertile
	public boolean    plAllowParallel      = true;  // Allow parallel shift of the specified planes before adding
	public double     plMaxDiff            = 0.3;   // Maximal normalized tile disparity difference from the plane to consider
	public double     plOtherDiff          = 1.4;   // Maximal difference of the added tile ratio to the average  disparity difference


	public boolean    plSplitXY            = true;  // Separate tiles for split planes by X, Y
	public double     plSplitXYTolerance   = 0.2;   // Disparity tolerance when separating by X, Y


	public boolean    plFuse               =  true; // Fuse planes together (off for debug only)
	public boolean    plKeepOrphans        =  true; // Keep unconnected supertiles
	public double     plMinOrphan          =  2.0;  // Minimal strength unconnected supertiles to keep

	public double     plSnapDispAny        =  .2;   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
	public double     plSnapStrengthAny    =  .2;   // Maximal strength to fit any distance (if does not fit otherwise - treat as zero str4ength
	public double     plSnapNegAny         =  .5;   // Maximal negative disparity difference from the best match
	public double     plSnapDispMax        =  .5;   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
	public double     plSnapDispWeight     =  .5;   // Maximal disparity diff. by weight product to snap to plane
	public int        plSnapZeroMode       =  1;    // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest

	public boolean    msUseSel             =   true; // Use planes selection masks (generated when splitting to intersecting pairs
	public boolean    msDivideByArea       =   true; // Divide plane strengths by ellipsoid area
	public double     msScaleProj          =   1.5;  // Scale projection of the plane ellipsoid
	public double     msFractUni           =   0.3;  // Spread this fraction of the ellipsoid weight among extended (double) supertile

	public boolean    tsNoEdge             = true;   // Do not assign tiles to the surface edges (not having all 8 neighbors)
	public boolean    tsUseCenter          = true;   // Only assign outside of 8x8 center if no suitable alternative
	public double     tsMaxDiff            = 0.3;    // Maximal disparity difference when assigning tiles
	public double     tsMinDiffOther       = 0.35;   // Minimal disparity difference to be considered as a competitor surface
	public double     tsMinStrength        = 0.05;   // Minimal tile correlation strength to be assigned
	public double     tsMaxStrength        = 10.0;   // Maximal tile correlation strength to be assigned
	public double     tsMinSurface         = 0.001;  // Minimal surface strength at the tile location
	public int        tsMoveDirs           = 3;      // Allowed tile disparity correction: 1 increase, 2 - decrease, 3 - both directions
	public double     tsSurfStrPow         = 0.0;    // Raise surface strengths ratio to this power when comparing candidates
	public double     tsAddStrength        = 0.01;   // Add to strengths when calculating pull of assigned tiles
	public double     tsSigma              = 2.0;    // Radius of influence (in tiles) of the previously assigned tiles
	public double     tsNSigma             = 2.0;    // Maximal relative to radius distance to calculate influence
	public double     tsMinPull            = 0.001;  // Additional pull of each surface
	public double     tsMinAdvantage       = 3.0;    // Minimal ratio of the best surface candidate to the next one to make selection

	public int        tsClustSize          = 1;      // Minimal size of a cluster to keep
	public double     tsClustWeight        = 0.2;    // Minimal total weight of a cluster to keep

	public int        tsMinNeib            = 6;      // Minimal number of neighbors of unassigned tile to join (the farthest)
	public double     tsMaxSurStrength     = 0.05;   // Maximal strength of the surrounded unassigned tile to join
	public boolean    tsCountDis           = true;   // Include disabled tiles/borders when counting assigned neighbors

	public boolean    tsEnPlaneSeed        = true;   // Assign tiles that were used to generate planes
	public boolean    tsEnOnly             = true;   // Allow assignment only surface
	public boolean    tsEnGrow             = true;   // Grow the only surface assignments
	public double     tsGrowStrength       = 0.01;   // Maximal strength when growing the only surfaces
	public boolean    tsGrowStrong         = true;   // Grow over strong if disparity matches
	public double     tsContStrength       = 0.1;    // Minimal strength to continue grow with disparity match
	public double     tsContDiff           = 0.1;    // Maximal normalized disparity error to grow over strong tiles


	public boolean    tsEnSingle           = true;   // Allow assignment to the nearest surface with no competitors
	public boolean    tsEnMulti            = true;   // Allow assignment when several surfaces fit

	public boolean    tsRemoveWeak1        = false;  // Remove weak clusters before growing
	public boolean    tsGrowSurround       = true;   // Assign tiles that have neighbors to the lowest disparity
	public boolean    tsRemoveWeak2        = true;   // Remove weak clusters after growing

	public boolean    tsLoopMulti          = true;   // Repeat multi-choice assignment while succeeding
	public boolean    tsReset              = false;  // Reset tiles to surfaces assignment
	public boolean    tsShow               = false;  // Show results of tiles to surfaces assignment
	public int        tsNumClust           = 500;     // Number of clusters to keep

	public int        tsConsensMode        = 7;      // Which assignments to match +1 - combo, +2 grown single, +4 plane seeds
	public int        tsConsensAgree       = 1;      // Minimal number of assignments to agree

	// Tile assignment parameters
	public double     taMinFgBg            = 0.1;    // Minimal foreground/ background separation to look for weak FG edge
	public double     taMinFgEdge          = 0.2;    // Minimal foreground edge strength (stronger edges will have proportionally smaller costs)
	public double     taMinColSep          = 0.05;   // Minimal surface separation that requires color change
	public double     taMinColDiff         = 0.01;   // Minimal color variation (larger proportionally reduces cost)
	public double     taOutlier            = 1.0;    // Disparity difference limit
	public double     taDiffPwr            = 0.25;   // Strength power when calculating disparity error
	public double     taBestPwr            = 0.0;    // Strength power when calculating disparity error over best
	public double     taDiff9Pwr           = 0.5;    // Strength power when calculating disparity error for group of 9
	public double     taColSigma           = 1.5;    // Gaussian sigma to blur color difference between tiles along each direction
	public double     taColFraction        = 0.3;    // Relative amount of the blurred color difference in the mixture


	public double     taCostEmpty          = 1.0;    // Cost of a tile that is not assigned
	public double     taCostNoLink         = 1.0;    // Cost of a tile not having any neighbor in particular direction
	public double     taCostSwitch         = 1.0;    // Cost of a tile switching to a neighbor that does not have a link
	public double     taCostColor          = 1.0;    // Cost of a tile switching to a disconnected neighbor divided by a color mismatch
	public double     taCostDiff           = 1.0;    // Cost of a weighted normalized tile disparity error
	public double     taCostDiffBest       = 1.0;    // Cost of a weighted normalized tile disparity error above best surface
	public double     taCostDiff9          = 1.0;    // Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)
	public double     taCostWeakFgnd       = 1.0;    // Cost of a weak foreground edge
	public double     taCostFlaps          = 1.0;    // Cost of using supertile "flaps" (not in the center 8x8 tiles area)
	public double     taCostMismatch       = 1.0;    // Cost of a measurement layer not having same layer in the same location or near

	public boolean    taEnEmpty            = true;   // Enable cost of a tile that is not assigned
	public boolean    taEnNoLink           = true;   // Enable cost of a tile not having any neighbor in particular direction
	public boolean    taEnSwitch           = true;   // Enable cost of a tile switching to a neighbor that does not have a link
	public boolean    taEnColor            = true;   // Enable cost of a tile switching to a disconnected neighbor divided by a color mismatch
	public boolean    taEnDiff             = true;   // Enable cost of a weighted normalized tile disparity error
	public boolean    taEnDiffBest         = true;   // Enable cost of a weighted normalized tile disparity error above best surface
	public boolean    taEnDiff9            = true;   // Enable cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)
	public boolean    taEnWeakFgnd         = true;   // Enable cost of a weak foreground edge
	public boolean    taEnFlaps            = true;   // Enable cost of using supertile "flaps" (not in the center 8x8 tiles area)
	public boolean    taEnMismatch         = false;  // Enable cost of a measurement layer not having same layer in the same location or near

// gpu processing parameters
	public int        gpu_corr_rad =        7;  // size of the correlation to save - initially only 15x15
	public double     gpu_weight_r =      0.25;
	public double     gpu_weight_b =      0.25; // weight g = 1.0 - gpu_weight_r - gpu_weight_b
	public double     gpu_sigma_r =       0.9; // 1.1;
	public double     gpu_sigma_b =       0.9; // 1.1;
	public double     gpu_sigma_g =       0.6; // 0.7;
	public double     gpu_sigma_m =       0.4; // 0.7;
	public double     gpu_sigma_rb_corr = 0.5; // apply LPF after accumulating R and B correlation before G,
	public double     gpu_sigma_corr =    0.9;
	public double     gpu_sigma_corr_m =  0.15;
	public double     gpu_fatz =          30.0;
	public double     gpu_fatz_m =        30.0;

	public boolean    gpu_woi =             false; // if true - use gpu_woi_tx, ...
	public int        gpu_woi_tx =              0;
	public int        gpu_woi_ty =              0;
	public int        gpu_woi_twidth =        324;
	public int        gpu_woi_theight =       242;
	public boolean    gpu_woi_round =       false;
	public boolean    gpu_save_ports_xy =   false; // debug feature - save calculated ports X,Y to compare with Java-generated
	public boolean    gpu_show_jtextures =  true;  // debug feature - show Java-generated textures from non-overlapping in GPU (will not generate if false)
	public boolean    gpu_show_extra =      true;  // show low-res data for macro
	
	public boolean    gpu_use_main =        false; // accelerate tile processor for the main quad camera
	public boolean    gpu_use_main_macro =  false; // accelerate tile processor for the main quad camera in macro mode
	public boolean    gpu_use_main_adjust = false; // accelerate tile processor for the main quad camera for field calibration
	public boolean    gpu_use_aux =         false; // accelerate tile processor for the aux (lwir) quad camera
	public boolean    gpu_use_aux_macro =   false; // accelerate tile processor for the aux (lwir) quad camera in macro mode
	public boolean    gpu_use_aux_adjust =  false; // accelerate tile processor for the aux (lwir) quad camera for field calibration

	public boolean    replaceWeakOutliers =   true; // false;

	public boolean    debug_initial_discriminate = false;

	public boolean    dbg_migrate =            true;

	// other debug images
	public int        dbg_early_exit =         0;     // Exit from early stages
	public boolean    show_first_bg =          false; // show extrinsic adjustment differences

	public boolean    show_extrinsic =         false; // show extrinsic adjustment differences
	public boolean    show_ortho_combine =     false; // Show 'ortho_combine'
	public boolean    show_refine_supertiles = false; // show 'refine_disparity_supertiles'
	public boolean    show_bgnd_nonbgnd =      false; // show 'bgnd_nonbgnd'
	public boolean    show_filter_scan =       false; // show 'FilterScan'
	public boolean    show_combined =          false; // show 'combo_scan' (combined multiple scans)
	public boolean    show_unique =            false; // show 'unique_scan' (removed already measured tiles with the same disparity)
	public boolean    show_histograms =        false; // show supertile disparity histograms
	public boolean    show_init_refine =       false; // show debug images during initial refinement
	public boolean    show_expand =            false; // show debug images during disparity expansion
	public boolean    show_variant =           false; // show prepareExpandVariant when elevating variant number
	public boolean    show_retry_far =         false; // show debug images related to retrying far tiles near foreground
	public boolean    show_macro =             false; // show debug images related to macro correlation

	public boolean    show_shells =            false; // show 'shells'
	public boolean    show_neighbors =         false; // show 'neighbors'
	public boolean    show_flaps_dirs =        false; // show 'flaps-dirs'
	public boolean    show_first_clusters =    false; // show 'first_N_clusters'
	public boolean    show_planes =            false; // show planes
	public double []  vertical_xyz =           {0.0,1.0,0.0}; // real world up unit vector in camera CS (x - right, y - up, z - to camera};

	public ImageDttParameters             img_dtt = new ImageDttParameters();
	public BiQuadParameters               rig =     new BiQuadParameters();
	public PoleProcessorParameters        poles =   new PoleProcessorParameters();
	public MeasuredLayersFilterParameters mlfp =    new MeasuredLayersFilterParameters();
	public LwirReaderParameters           lwir =    new LwirReaderParameters();


	public HashMap<String,Double> z_corr_map = new HashMap<String,Double>(); //old one
	public HashMap<String,Double> infinity_distace_map = new HashMap<String,Double>(); //new one
	public static String Z_CORR_PREFIX = "z_corr.";
	public static String INFINITY_DISTANCE_PREFIX = "infinity_distance.";

	public boolean   batch_run =               false; // turned on only while running in batch mode

	public boolean useGPU() {
		return  useGPU(false) || useGPU(true);
	
	}
	public boolean useGPU(boolean aux) {
		return  aux? (gpu_use_aux  || gpu_use_aux_macro || gpu_use_aux_adjust)
				: (gpu_use_main || gpu_use_main_macro || gpu_use_main_adjust);
	}

	
	
	public double getCorrSigma(boolean monochrome) {
		return monochrome ? corr_sigma_mono : corr_sigma;
	}
	public double getFatZero(boolean monochrome) {
		return monochrome ? fat_zero_mono : fat_zero;
	}

	public double getGpuFatZero(boolean monochrome) {
		return monochrome ? gpu_fatz_m : gpu_fatz;
	}

	public double getGpuCorrSigma(boolean monochrome) {
		return monochrome ? gpu_sigma_corr_m : gpu_sigma_corr;
	}

	public double getGpuCorrRBSigma(boolean monochrome) {
		return monochrome ? 1.0 : gpu_sigma_rb_corr;
	}

	public double getScaleStrength(boolean aux) {
		return aux ? scale_strength_aux : scale_strength_main;
	}

	public double getLymChange(boolean aux) {
		return aux ? lym_change_aux : lym_change;
	}

	public int getLyPerQuad (int num_tiles) {
		return (int) Math.max(ly_per_quad_r * num_tiles , ly_per_quad );
	}
	public int getLyInf (int num_tiles) {
		return (int) Math.max(ly_inf_r * num_tiles , ly_inf );
	}
	public int getLyInfScale (int num_tiles) {
		return (int) Math.max(ly_inf_scale_r * num_tiles , ly_inf_scale );
	}

	public CLTParameters(){}
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"transform_size",             this.transform_size+"");
		properties.setProperty(prefix+"clt_window",                 this.clt_window+"");
		properties.setProperty(prefix+"shift_x",                    this.shift_x+"");
		properties.setProperty(prefix+"shift_y",                    this.shift_y+"");
		properties.setProperty(prefix+"tileStep",                   this.tileStep+"");
		properties.setProperty(prefix+"iclt_mask",                  this.iclt_mask+"");
		properties.setProperty(prefix+"tileX",                      this.tileX+"");
		properties.setProperty(prefix+"tileY",                      this.tileY+"");
		properties.setProperty(prefix+"dbg_mode",                   this.dbg_mode+"");
		properties.setProperty(prefix+"ishift_x",                   this.ishift_x+"");
		properties.setProperty(prefix+"ishift_y",                   this.ishift_y+"");
		properties.setProperty(prefix+"fat_zero",                   this.fat_zero+"");
		properties.setProperty(prefix+"fat_zero_mono",              this.fat_zero_mono+"");
		properties.setProperty(prefix+"corr_sigma",                 this.corr_sigma+"");
		properties.setProperty(prefix+"corr_sigma_mono",            this.corr_sigma_mono+"");
		properties.setProperty(prefix+"scale_strength_main",        this.scale_strength_main+"");
		properties.setProperty(prefix+"scale_strength_aux",         this.scale_strength_aux+"");

		properties.setProperty(prefix+"norm_kern",                  this.norm_kern+"");
		properties.setProperty(prefix+"gain_equalize",              this.gain_equalize+"");
		properties.setProperty(prefix+"colors_equalize",            this.colors_equalize+"");
		properties.setProperty(prefix+"nosat_equalize",             this.nosat_equalize+"");
		properties.setProperty(prefix+"sat_level",                  this.sat_level+"");
		properties.setProperty(prefix+"max_overexposure",           this.max_overexposure+"");

		properties.setProperty(prefix+"novignetting_r",             this.novignetting_r+"");
		properties.setProperty(prefix+"novignetting_g",             this.novignetting_g+"");
		properties.setProperty(prefix+"novignetting_b",             this.novignetting_b+"");
		properties.setProperty(prefix+"scale_r",                    this.scale_r+"");
		properties.setProperty(prefix+"scale_g",                    this.scale_g+"");
		properties.setProperty(prefix+"scale_b",                    this.scale_b+"");
		properties.setProperty(prefix+"vignetting_max",             this.vignetting_max+"");
		properties.setProperty(prefix+"vignetting_range",           this.vignetting_range+"");
		properties.setProperty(prefix+"kernel_step",                this.kernel_step+"");
		properties.setProperty(prefix+"disparity",                  this.disparity +"");
		properties.setProperty(prefix+"z_correction",               this.z_correction +"");
		properties.setProperty(prefix+"correlate",                  this.correlate+"");
		properties.setProperty(prefix+"corr_mask",                  this.corr_mask+"");
		properties.setProperty(prefix+"corr_sym",                   this.corr_sym+"");
		properties.setProperty(prefix+"corr_keep",                  this.corr_keep+"");
		properties.setProperty(prefix+"corr_show",                  this.corr_show+"");
		properties.setProperty(prefix+"corr_mismatch",              this.corr_mismatch+"");
		properties.setProperty(prefix+"corr_offset",                this.corr_offset +"");
		properties.setProperty(prefix+"corr_red",                   this.corr_red +"");
		properties.setProperty(prefix+"corr_blue",                  this.corr_blue +"");
		properties.setProperty(prefix+"corr_normalize",             this.corr_normalize+"");
		properties.setProperty(prefix+"min_corr",                   this.min_corr +"");
		properties.setProperty(prefix+"min_corr_normalized",        this.min_corr_normalized +"");
		properties.setProperty(prefix+"max_corr_sigma",             this.max_corr_sigma +"");
		properties.setProperty(prefix+"max_corr_radius",            this.max_corr_radius +"");

		//  			properties.setProperty(prefix+"enhortho_width",             this.enhortho_width +"");
		//  			properties.setProperty(prefix+"enhortho_scale",             this.enhortho_scale +"");


		properties.setProperty(prefix+"max_corr_double",            this.max_corr_double+"");
		properties.setProperty(prefix+"corr_mode",                  this.corr_mode+"");
		properties.setProperty(prefix+"corr_border_contrast",       this.corr_border_contrast +"");
		properties.setProperty(prefix+"tile_task_op",               this.tile_task_op+"");
		properties.setProperty(prefix+"tile_task_wl",               this.tile_task_wl+"");
		properties.setProperty(prefix+"tile_task_wt",               this.tile_task_wt+"");
		properties.setProperty(prefix+"tile_task_ww",               this.tile_task_ww+"");
		properties.setProperty(prefix+"tile_task_wh",               this.tile_task_wh+"");
		properties.setProperty(prefix+"min_shot",                   this.min_shot +"");
		properties.setProperty(prefix+"scale_shot",                 this.scale_shot +"");
		properties.setProperty(prefix+"diff_sigma",                 this.diff_sigma +"");
		properties.setProperty(prefix+"diff_threshold",             this.diff_threshold +"");
		properties.setProperty(prefix+"diff_gauss",                 this.diff_gauss+"");
		properties.setProperty(prefix+"min_agree",                  this.min_agree +"");
		properties.setProperty(prefix+"dust_remove",                this.dust_remove+"");
		properties.setProperty(prefix+"black_back",                 this.black_back+"");
		properties.setProperty(prefix+"keep_weights",               this.keep_weights+"");
		properties.setProperty(prefix+"sharp_alpha",                this.sharp_alpha+"");

		properties.setProperty(prefix+"alpha0",                     this.alpha0 +"");
		properties.setProperty(prefix+"alpha1",                     this.alpha1 +"");

		properties.setProperty(prefix+"gen_chn_stacks",             this.gen_chn_stacks+"");
		properties.setProperty(prefix+"gen_chn_img",                this.gen_chn_img+"");
		properties.setProperty(prefix+"gen_4_img",                  this.gen_4_img+"");
		properties.setProperty(prefix+"show_nonoverlap",            this.show_nonoverlap+"");
		properties.setProperty(prefix+"show_overlap",               this.show_overlap+"");
		properties.setProperty(prefix+"show_rgba_color",            this.show_rgba_color+"");
		properties.setProperty(prefix+"show_map",                   this.show_map+"");
		properties.setProperty(prefix+"show_corr",                  this.show_corr+"");
		properties.setProperty(prefix+"disp_scan_start",            this.disp_scan_start +"");
		properties.setProperty(prefix+"disp_scan_step",             this.disp_scan_step +"");
		properties.setProperty(prefix+"disp_scan_count",            this.disp_scan_count+"");

		properties.setProperty(prefix+"fine_dbg",                   this.fine_dbg+"");
		properties.setProperty(prefix+"fine_corr_x_0",              this.fine_corr_x_0 +"");
		properties.setProperty(prefix+"fine_corr_y_0",              this.fine_corr_y_0 +"");
		properties.setProperty(prefix+"fine_corr_x_1",              this.fine_corr_x_1 +"");
		properties.setProperty(prefix+"fine_corr_y_1",              this.fine_corr_y_1 +"");
		properties.setProperty(prefix+"fine_corr_x_2",              this.fine_corr_x_2 +"");
		properties.setProperty(prefix+"fine_corr_y_2",              this.fine_corr_y_2 +"");
		properties.setProperty(prefix+"fine_corr_x_3",              this.fine_corr_x_3 +"");
		properties.setProperty(prefix+"fine_corr_y_3",              this.fine_corr_y_3 +"");
		properties.setProperty(prefix+"fine_corr_ignore",           this.fine_corr_ignore+"");
		properties.setProperty(prefix+"fine_corr_apply",            this.fine_corr_apply+"");

		properties.setProperty(prefix+"fcorr_radius",               this.fcorr_radius +"");
		properties.setProperty(prefix+"fcorr_min_strength",         this.fcorr_min_strength +"");
		properties.setProperty(prefix+"fcorr_disp_diff",            this.fcorr_disp_diff +"");
		properties.setProperty(prefix+"fcorr_quadratic",            this.fcorr_quadratic+"");
		properties.setProperty(prefix+"fcorr_ignore",               this.fcorr_ignore+"");

		properties.setProperty(prefix+"fcorr_inf_strength",         this.fcorr_inf_strength +"");
		properties.setProperty(prefix+"fcorr_inf_diff",             this.fcorr_inf_diff +"");
		properties.setProperty(prefix+"fcorr_inf_quad",             this.fcorr_inf_quad+"");
		properties.setProperty(prefix+"fcorr_inf_vert",             this.fcorr_inf_vert+"");

		properties.setProperty(prefix+"inf_disp_apply",             this.inf_disp_apply+"");
		properties.setProperty(prefix+"inf_repeat",                 this.inf_repeat+"");

		//			properties.setProperty(prefix+"inf_mism_apply",             this.inf_mism_apply+"");
		properties.setProperty(prefix+"inf_iters",                  this.inf_iters+"");
		properties.setProperty(prefix+"inf_final_diff",             this.inf_final_diff +"");
		properties.setProperty(prefix+"inf_far_pull",               this.inf_far_pull +"");
		properties.setProperty(prefix+"inf_str_pow",                this.inf_str_pow +"");
		properties.setProperty(prefix+"inf_smpl_side",              this.inf_smpl_side+"");
		properties.setProperty(prefix+"inf_smpl_num",               this.inf_smpl_num+"");
		properties.setProperty(prefix+"inf_smpl_rms",               this.inf_smpl_rms +"");

		properties.setProperty(prefix+"ih_smpl_step",               this.ih_smpl_step+"");
		properties.setProperty(prefix+"ih_disp_min",                this.ih_disp_min +"");
		properties.setProperty(prefix+"ih_disp_step",               this.ih_disp_step +"");
		properties.setProperty(prefix+"ih_num_bins",                this.ih_num_bins+"");
		properties.setProperty(prefix+"ih_sigma",                   this.ih_sigma +"");
		properties.setProperty(prefix+"ih_max_diff",                this.ih_max_diff +"");
		properties.setProperty(prefix+"ih_min_samples",             this.ih_min_samples+"");
		properties.setProperty(prefix+"ih_norm_center",             this.ih_norm_center+"");
		properties.setProperty(prefix+"inf_restore_disp",           this.inf_restore_disp+"");

		properties.setProperty(prefix+"ly_lma_ers",                 this.ly_lma_ers+"");
		properties.setProperty(prefix+"ly_gt_strength",             this.ly_gt_strength+"");
		properties.setProperty(prefix+"ly_gt_use_wnd",              this.ly_gt_use_wnd+"");
		properties.setProperty(prefix+"ly_gt_rms",                  this.ly_gt_rms+"");

		properties.setProperty(prefix+"lylw_inf_en",                  this.lylw_inf_en+"");
		properties.setProperty(prefix+"lylw_aztilt_en",               this.lylw_aztilt_en+"");
		properties.setProperty(prefix+"lylw_diff_roll_en",            this.lylw_diff_roll_en+"");
		properties.setProperty(prefix+"lylw_focalLength",             this.lylw_focalLength+"");
		properties.setProperty(prefix+"lylw_com_roll",                this.lylw_com_roll+"");
		properties.setProperty(prefix+"lylw_par_sel",                 this.lylw_par_sel+"");

		properties.setProperty(prefix+"ly_marg_fract",              this.ly_marg_fract+"");
		properties.setProperty(prefix+"ly_on_scan",                 this.ly_on_scan+"");
		properties.setProperty(prefix+"ly_inf_en",                  this.ly_inf_en+"");
		properties.setProperty(prefix+"ly_min_forced",              this.ly_min_forced+"");
		properties.setProperty(prefix+"ly_aztilt_en",               this.ly_aztilt_en+"");
		properties.setProperty(prefix+"ly_diff_roll_en",            this.ly_diff_roll_en+"");
		properties.setProperty(prefix+"ly_focalLength",             this.ly_focalLength+"");
		properties.setProperty(prefix+"ly_com_roll",                this.ly_com_roll+"");
		properties.setProperty(prefix+"ly_ers_rot",                 this.ly_ers_rot+"");
		properties.setProperty(prefix+"ly_ers_forw",                this.ly_ers_forw+"");
		properties.setProperty(prefix+"ly_ers_side",                this.ly_ers_side+"");
		properties.setProperty(prefix+"ly_ers_vert",                this.ly_ers_vert+"");



		properties.setProperty(prefix+"ly_par_sel",                 this.ly_par_sel+"");
		properties.setProperty(prefix+"ly_debug_level",             this.ly_debug_level+"");


		properties.setProperty(prefix+"ly_right_left",              this.ly_right_left+"");

		properties.setProperty(prefix+"ly_per_quad",                this.ly_per_quad +"");
		properties.setProperty(prefix+"ly_per_quad_r",              this.ly_per_quad_r +"");
		properties.setProperty(prefix+"ly_inf",                     this.ly_inf +"");
		properties.setProperty(prefix+"ly_inf_r",                   this.ly_inf_r +"");
		properties.setProperty(prefix+"ly_inf_scale",               this.ly_inf_scale +"");
		properties.setProperty(prefix+"ly_inf_scale_r",             this.ly_inf_scale_r +"");

		properties.setProperty(prefix+"ly_inf_frac",                this.ly_inf_frac +"");

		properties.setProperty(prefix+"ly_inf_max_disparity",       this.ly_inf_max_disparity +"");
		properties.setProperty(prefix+"ly_inf_disp",                this.ly_inf_disp+"");
		properties.setProperty(prefix+"ly_inf_force",               this.ly_inf_force+"");
		properties.setProperty(prefix+"ly_poly",                    this.ly_poly+"");

		properties.setProperty(prefix+"ly_smpl_side",               this.ly_smpl_side+"");
		properties.setProperty(prefix+"ly_smpl_num",                this.ly_smpl_num+"");
		properties.setProperty(prefix+"ly_smpl_rms",                this.ly_smpl_rms +"");
		properties.setProperty(prefix+"ly_disp_var",                this.ly_disp_var +"");
		properties.setProperty(prefix+"ly_disp_rvar",               this.ly_disp_rvar +"");
		properties.setProperty(prefix+"ly_disp_var_gt",             this.ly_disp_var_gt +"");
		properties.setProperty(prefix+"ly_disp_rvar_gt",            this.ly_disp_rvar_gt +"");
		properties.setProperty(prefix+"ly_norm_disp",               this.ly_norm_disp +"");
		properties.setProperty(prefix+"lym_overexp",                this.lym_overexp +"");
		properties.setProperty(prefix+"lym_update_disp",            this.lym_update_disp+"");
		properties.setProperty(prefix+"lym_iter",                   this.lym_iter+"");
		properties.setProperty(prefix+"lym_change",                 this.lym_change +"");
		properties.setProperty(prefix+"lym_change_aux",             this.lym_change_aux +"");
		properties.setProperty(prefix+"lym_poly_change",            this.lym_poly_change +"");

		properties.setProperty(prefix+"lyf_filter",                 this.lyf_filter+"");
		properties.setProperty(prefix+"lyf_smpl_side",              this.lyf_smpl_side+"");
		properties.setProperty(prefix+"lyf_rms_max",                this.lyf_rms_max +"");
		properties.setProperty(prefix+"lyf_frac_keep",              this.lyf_frac_keep +"");
		properties.setProperty(prefix+"lyf_min_samples",            this.lyf_min_samples+"");
		properties.setProperty(prefix+"lyf_norm_center",            this.lyf_norm_center+"");
		properties.setProperty(prefix+"ly_corr_scale",              this.ly_corr_scale +"");

		properties.setProperty(prefix+"lyr_filter_ds",              this.lyr_filter_ds +"");
		properties.setProperty(prefix+"lyr_filter_lyf",             this.lyr_filter_lyf +"");

		properties.setProperty(prefix+"corr_magic_scale",           this.corr_magic_scale +"");
		properties.setProperty(prefix+"corr_select",                this.corr_select +"");

		properties.setProperty(prefix+"show_textures",              this.show_textures+"");
		properties.setProperty(prefix+"debug_filters",              this.debug_filters+"");

		properties.setProperty(prefix+"min_smth",                   this.min_smth +"");
		properties.setProperty(prefix+"sure_smth",                  this.sure_smth +"");
		properties.setProperty(prefix+"bgnd_range",                 this.bgnd_range +"");
		properties.setProperty(prefix+"other_range",                this.other_range +"");

		properties.setProperty(prefix+"ex_strength",                this.ex_strength +"");
		properties.setProperty(prefix+"ex_nstrength",               this.ex_nstrength +"");

		properties.setProperty(prefix+"ex_over_bgnd",               this.ex_over_bgnd+"");
		properties.setProperty(prefix+"ex_min_over",                this.ex_min_over +"");

		properties.setProperty(prefix+"pt_super_trust",             this.pt_super_trust +"");
		properties.setProperty(prefix+"pt_keep_raw_fg",             this.pt_keep_raw_fg+"");
		properties.setProperty(prefix+"pt_scale_pre",               this.pt_scale_pre +"");
		properties.setProperty(prefix+"pt_scale_post",              this.pt_scale_post +"");

		properties.setProperty(prefix+"bgnd_sure",                  this.bgnd_sure +"");
		properties.setProperty(prefix+"bgnd_maybe",                 this.bgnd_maybe +"");
		properties.setProperty(prefix+"min_clstr_seed",             this.min_clstr_seed+"");
		properties.setProperty(prefix+"min_clstr_lone",             this.min_clstr_lone+"");
		properties.setProperty(prefix+"min_clstr_weight",           this.min_clstr_weight +"");
		properties.setProperty(prefix+"min_clstr_max",              this.min_clstr_max +"");

		properties.setProperty(prefix+"fill_gaps",                  this.fill_gaps+"");
		properties.setProperty(prefix+"fill_final",                 this.fill_final+"");
		properties.setProperty(prefix+"min_clstr_block",            this.min_clstr_block+"");
		properties.setProperty(prefix+"bgnd_grow",                  this.bgnd_grow+"");

		properties.setProperty(prefix+"ortho_old",                  this.ortho_old+"");
		properties.setProperty(prefix+"ortho_min_hor",              this.ortho_min_hor +"");
		properties.setProperty(prefix+"ortho_min_vert",             this.ortho_min_vert +"");
		properties.setProperty(prefix+"ortho_asym",                 this.ortho_asym +"");
		properties.setProperty(prefix+"ortho_over4",                this.ortho_over4 +"");
		properties.setProperty(prefix+"ortho_sustain",              this.ortho_sustain +"");
		properties.setProperty(prefix+"ortho_run",                  this.ortho_run+"");
		properties.setProperty(prefix+"ortho_minmax",               this.ortho_minmax +"");
		properties.setProperty(prefix+"ortho_bridge",               this.ortho_bridge+"");
		properties.setProperty(prefix+"ortho_rms",                  this.ortho_rms +"");
		properties.setProperty(prefix+"ortho_half_length",          this.ortho_half_length+"");
		properties.setProperty(prefix+"ortho_mix",                  this.ortho_mix +"");


		properties.setProperty(prefix+"or_hor",                     this.or_hor+"");
		properties.setProperty(prefix+"or_vert",                    this.or_vert+"");
		properties.setProperty(prefix+"or_sigma",                   this.or_sigma +"");
		properties.setProperty(prefix+"or_sharp",                   this.or_sharp +"");
		properties.setProperty(prefix+"or_scale",                   this.or_scale +"");
		properties.setProperty(prefix+"or_offset",                  this.or_offset +"");
		properties.setProperty(prefix+"or_asym",                    this.or_asym +"");
		properties.setProperty(prefix+"or_threshold",               this.or_threshold +"");
		properties.setProperty(prefix+"or_absHor",                  this.or_absHor +"");
		properties.setProperty(prefix+"or_absVert",                 this.or_absVert +"");
		properties.setProperty(prefix+"or_maxDisp",                 this.or_maxDisp +"");

		properties.setProperty(prefix+"poles_fix",                  this.poles_fix+"");
		properties.setProperty(prefix+"poles_len",                  this.poles_len+"");
		properties.setProperty(prefix+"poles_ratio",                this.poles_ratio +"");
		properties.setProperty(prefix+"poles_min_strength",         this.poles_min_strength +"");
		properties.setProperty(prefix+"or_maxDisp",                 this.or_maxDisp+"");



		properties.setProperty(prefix+"max_clusters",               this.max_clusters+"");
		properties.setProperty(prefix+"remove_scans",               this.remove_scans+"");
		properties.setProperty(prefix+"output_x3d",                 this.output_x3d+"");
		properties.setProperty(prefix+"output_obj",                 this.output_obj+"");
		properties.setProperty(prefix+"correct_distortions",        this.correct_distortions+"");
		properties.setProperty(prefix+"show_triangles",             this.show_triangles+"");
		properties.setProperty(prefix+"avg_cluster_disp",           this.avg_cluster_disp+"");
		properties.setProperty(prefix+"maxDispTriangle",            this.maxDispTriangle +"");
		properties.setProperty(prefix+"infinityDistance",           this.infinityDistance +"");
		properties.setProperty(prefix+"min_bgnd_tiles",             this.min_bgnd_tiles+"");
		properties.setProperty(prefix+"shUseFlaps",                 this.shUseFlaps+"");
		properties.setProperty(prefix+"shAggrFade",                 this.shAggrFade+"");
		properties.setProperty(prefix+"shMinArea",                  this.shMinArea+"");
		properties.setProperty(prefix+"shMinStrength",              this.shMinStrength +"");
		properties.setProperty(prefix+"tiRigidVertical",            this.tiRigidVertical +"");
		properties.setProperty(prefix+"tiRigidHorizontal",          this.tiRigidHorizontal +"");
		properties.setProperty(prefix+"tiRigidDiagonal",            this.tiRigidDiagonal +"");
		properties.setProperty(prefix+"tiStrengthOffset",           this.tiStrengthOffset +"");
		properties.setProperty(prefix+"tiDispScale",                this.tiDispScale +"");
		properties.setProperty(prefix+"tiDispPow",                  this.tiDispPow +"");
		properties.setProperty(prefix+"tiDispPull",                 this.tiDispPull +"");
		properties.setProperty(prefix+"tiDispPullPreFinal",          this.tiDispPullPreFinal +"");
		properties.setProperty(prefix+"tiDispPullFinal",            this.tiDispPullFinal +"");
		properties.setProperty(prefix+"tiBreakNorm",                this.tiBreakNorm +"");
		properties.setProperty(prefix+"tiBreak3",                   this.tiBreak3 +"");
		properties.setProperty(prefix+"tiBreak31",                  this.tiBreak31 +"");
		properties.setProperty(prefix+"tiBreak21",                  this.tiBreak21 +"");
		properties.setProperty(prefix+"tiBreakFar",                 this.tiBreakFar +"");
		properties.setProperty(prefix+"tiBreakNear",                this.tiBreakNear +"");
		properties.setProperty(prefix+"tiBreakMode",                this.tiBreakMode +"");
		properties.setProperty(prefix+"tiBreakSame",                this.tiBreakSame +"");
		properties.setProperty(prefix+"tiBreakTurn",                this.tiBreakTurn +"");

		properties.setProperty(prefix+"tiHealPreLast",              this.tiHealPreLast +"");
		properties.setProperty(prefix+"tiHealLast",                 this.tiHealLast +"");
		properties.setProperty(prefix+"tiHealSame",                 this.tiHealSame+"");

		properties.setProperty(prefix+"tiIterations",               this.tiIterations+"");
		properties.setProperty(prefix+"tiPrecision",                this.tiPrecision+"");
		properties.setProperty(prefix+"tiNumCycles",                this.tiNumCycles+"");

		properties.setProperty(prefix+"stUseRefine",                this.stUseRefine+"");
		properties.setProperty(prefix+"stUsePass2",                 this.stUsePass2+"");
		properties.setProperty(prefix+"stUseRender",                this.stUseRender+"");

		properties.setProperty(prefix+"stShow",                     this.stShow+"");
		properties.setProperty(prefix+"stSize",                     this.stSize+"");
		properties.setProperty(prefix+"stStepFar",                  this.stStepFar +"");
		properties.setProperty(prefix+"stStepNear",                 this.stStepNear +"");
		properties.setProperty(prefix+"stStepThreshold",            this.stStepThreshold +"");
		properties.setProperty(prefix+"stMinDisparity",             this.stMinDisparity +"");

		///			properties.setProperty(prefix+"stFloor",                    this.stFloor +"");   // mlfp.strength_floor
		///			properties.setProperty(prefix+"stPow",                      this.stPow +"");     // mlfp.strength_pow
		properties.setProperty(prefix+"stSigma",                    this.stSigma +"");          // keep
		properties.setProperty(prefix+"stMinBgDisparity",           this.stMinBgDisparity +""); // keep
		properties.setProperty(prefix+"stMinBgFract",               this.stMinBgFract +"");     // keep
		properties.setProperty(prefix+"stUseDisp",                  this.stUseDisp +"");        // keep
		properties.setProperty(prefix+"stStrengthScale",            this.stStrengthScale +"");  // keep

		properties.setProperty(prefix+"stSmplMode",                 this.stSmplMode+"");        // keep

		///			properties.setProperty(prefix+"stSmplSide",                 this.stSmplSide+"");        // mlfp.smplSide
		///			properties.setProperty(prefix+"stSmplNum",                  this.stSmplNum+"");         // mlfp.smplNum
		///			properties.setProperty(prefix+"stSmplRms",                  this.stSmplRms +"");        // mlfp.smplRms
		///			properties.setProperty(prefix+"stSmplWnd",                  this.stSmplWnd+"");         // mlfp.smplWnd

		///			properties.setProperty(prefix+"fs_max_abs_tilt",            this.fs_max_abs_tilt+"");   // mlfp.max_abs_tilt
		///			properties.setProperty(prefix+"fs_max_rel_tilt",            this.fs_max_rel_tilt+"");   // mlfp.max_rel_tilt
		///			properties.setProperty(prefix+"fs_damp_tilt",               this.fs_damp_tilt+"");      // mlfp.damp_tilt
		///			properties.setProperty(prefix+"fs_min_tilt_disp",           this.fs_min_tilt_disp+"");  // mlfp.min_tilt_disp
		///			properties.setProperty(prefix+"fs_transition",              this.fs_transition+"");     // mlfp.transition
		///			properties.setProperty(prefix+"fs_far_mode",                this.fs_far_mode+"");       // mlfp.far_mode
		///			properties.setProperty(prefix+"fs_far_power",               this.fs_far_power+"");      // mlfp.far_power

		properties.setProperty(prefix+"stGrowSel",                  this.stGrowSel+"");
		properties.setProperty(prefix+"stMeasSel",                  this.stMeasSel+"");
		properties.setProperty(prefix+"stSmallDiff",                this.stSmallDiff +"");
		properties.setProperty(prefix+"stHighMix",                  this.stHighMix +"");


		properties.setProperty(prefix+"outlierStrength",            this.outlierStrength +"");
		properties.setProperty(prefix+"outlierDiff",                this.outlierDiff +"");
		properties.setProperty(prefix+"outlierDiffPos",             this.outlierDiffPos +"");
		properties.setProperty(prefix+"outlierDiffNeg",             this.outlierDiffNeg +"");

		properties.setProperty(prefix+"combine_refine",             this.combine_refine+"");

		properties.setProperty(prefix+"combine_min_strength",       this.combine_min_strength +"");
		properties.setProperty(prefix+"combine_min_hor",            this.combine_min_hor +"");
		properties.setProperty(prefix+"combine_min_vert",           this.combine_min_vert +"");
		//			properties.setProperty(prefix+"unique_tolerance",           this.unique_tolerance +"");
		properties.setProperty(prefix+"grow_sweep",                 this.grow_sweep+"");
		properties.setProperty(prefix+"grow_disp_max",              this.grow_disp_max +"");
		properties.setProperty(prefix+"grow_disp_trust",            this.grow_disp_trust +"");
		properties.setProperty(prefix+"grow_disp_step",             this.grow_disp_step +"");
		properties.setProperty(prefix+"grow_min_diff",              this.grow_min_diff +"");
		properties.setProperty(prefix+"grow_retry_far",             this.grow_retry_far+"");
		properties.setProperty(prefix+"grow_pedantic",              this.grow_pedantic+"");
		properties.setProperty(prefix+"grow_retry_inf",             this.grow_retry_inf+"");

		properties.setProperty(prefix+"gr_new_expand",              this.gr_new_expand+"");
		//  			properties.setProperty(prefix+"gr_max_expand",              this.gr_max_expand+"");
		properties.setProperty(prefix+"fds_str_floor",              this.fds_str_floor +"");
		properties.setProperty(prefix+"gr_ovrbg_cmb",               this.gr_ovrbg_cmb +"");
		properties.setProperty(prefix+"gr_ovrbg_cmb_hor",           this.gr_ovrbg_cmb_hor +"");
		properties.setProperty(prefix+"gr_ovrbg_cmb_vert",          this.gr_ovrbg_cmb_vert +"");
		properties.setProperty(prefix+"gr_ovrbg_filtered",          this.gr_ovrbg_filtered +"");
		properties.setProperty(prefix+"fds_str_pow",                this.fds_str_pow +"");
		properties.setProperty(prefix+"fds_smpl_side",              this.fds_smpl_side+"");
		properties.setProperty(prefix+"fds_smpl_num",               this.fds_smpl_num+"");
		properties.setProperty(prefix+"fds_smpl_rms",               this.fds_smpl_rms +"");
		properties.setProperty(prefix+"fds_smpl_rel_rms",           this.fds_smpl_rel_rms +"");
		properties.setProperty(prefix+"fds_smpl_wnd",               this.fds_smpl_wnd+"");
		properties.setProperty(prefix+"fds_abs_tilt",               this.fds_abs_tilt +"");
		properties.setProperty(prefix+"fds_rel_tilt",               this.fds_rel_tilt +"");

		properties.setProperty(prefix+"per_filter",                 this.per_filter +"");
		properties.setProperty(prefix+"per_trustedCorrelation",     this.per_trustedCorrelation +"");
		properties.setProperty(prefix+"per_initial_diff",           this.per_initial_diff +"");
		properties.setProperty(prefix+"per_strength_floor",         this.per_strength_floor +"");
		properties.setProperty(prefix+"per_strength_max_over",      this.per_strength_max_over +"");
		properties.setProperty(prefix+"per_min_period",             this.per_min_period +"");
		properties.setProperty(prefix+"per_min_num_periods",        this.per_min_num_periods +"");
		properties.setProperty(prefix+"per_disp_tolerance",         this.per_disp_tolerance +"");
		properties.setProperty(prefix+"per_disp_match",             this.per_disp_match +"");
		properties.setProperty(prefix+"per_strong_match_inc",       this.per_strong_match_inc +"");

		properties.setProperty(prefix+"mc_disp8_step",              this.mc_disp8_step +"");
		properties.setProperty(prefix+"mc_disp8_trust",             this.mc_disp8_trust +"");
		properties.setProperty(prefix+"mc_strength",                this.mc_strength +"");
		properties.setProperty(prefix+"mc_unique_tol",              this.mc_unique_tol +"");

		properties.setProperty(prefix+"mc_trust_fin",               this.mc_trust_fin +"");
		properties.setProperty(prefix+"mc_trust_sigma",             this.mc_trust_sigma +"");
		properties.setProperty(prefix+"mc_ortho_weight",            this.mc_ortho_weight +"");
		properties.setProperty(prefix+"mc_diag_weight",             this.mc_diag_weight +"");
		properties.setProperty(prefix+"mc_gap",                     this.mc_gap +"");

		properties.setProperty(prefix+"mc_weight_var",              this.mc_weight_var +"");
		properties.setProperty(prefix+"mc_weight_Y",                this.mc_weight_Y +"");
		properties.setProperty(prefix+"mc_weight_RBmG",             this.mc_weight_RBmG +"");

		properties.setProperty(prefix+"gr_min_new",                 this.gr_min_new+"");
		properties.setProperty(prefix+"gr_var_new_sngl",            this.gr_var_new_sngl+"");
		properties.setProperty(prefix+"gr_var_new_fg",              this.gr_var_new_fg+"");
		properties.setProperty(prefix+"gr_var_all_fg",              this.gr_var_all_fg+"");
		properties.setProperty(prefix+"gr_var_new_bg",              this.gr_var_new_bg+"");
		properties.setProperty(prefix+"gr_var_all_bg",              this.gr_var_all_bg+"");
		properties.setProperty(prefix+"gr_var_next",                this.gr_var_next+"");
		properties.setProperty(prefix+"gr_num_steps",               this.gr_num_steps+"");
		properties.setProperty(prefix+"gr_steps_over",              this.gr_steps_over+"");
		properties.setProperty(prefix+"gr_smpl_size",               this.gr_smpl_size+"");
		properties.setProperty(prefix+"gr_min_pnts",                this.gr_min_pnts+"");
		properties.setProperty(prefix+"gr_use_wnd",                 this.gr_use_wnd+"");
		properties.setProperty(prefix+"gr_tilt_damp",               this.gr_tilt_damp +"");
		properties.setProperty(prefix+"gr_split_rng",               this.gr_split_rng +"");
		properties.setProperty(prefix+"gr_same_rng",                this.gr_same_rng +"");
		properties.setProperty(prefix+"gr_diff_cont",               this.gr_diff_cont +"");
		properties.setProperty(prefix+"gr_abs_tilt",                this.gr_abs_tilt +"");
		properties.setProperty(prefix+"gr_rel_tilt",                this.gr_rel_tilt +"");
		properties.setProperty(prefix+"gr_smooth",                  this.gr_smooth+"");
		properties.setProperty(prefix+"gr_fin_diff",                this.gr_fin_diff +"");
		properties.setProperty(prefix+"gr_unique_tol",              this.gr_unique_tol +"");
		properties.setProperty(prefix+"gr_unique_pretol",           this.gr_unique_pretol +"");

		properties.setProperty(prefix+"ft_mod_strength",                    this.ft_mod_strength +"");
		properties.setProperty(prefix+"ft_clusterize_by_highest",           this.ft_clusterize_by_highest +"");
		properties.setProperty(prefix+"ft_clust_sigma",                     this.ft_clust_sigma +"");
		properties.setProperty(prefix+"ft_disp_arange_vert",                this.ft_disp_arange_vert +"");
		properties.setProperty(prefix+"ft_disp_rrange_vert",                this.ft_disp_rrange_vert +"");
		properties.setProperty(prefix+"ft_disp_arange_hor",                 this.ft_disp_arange_hor +"");
		properties.setProperty(prefix+"ft_disp_rrange_hor",                  this.ft_disp_rrange_hor +"");
		properties.setProperty(prefix+"ft_tolerance_above_near",            this.ft_tolerance_above_near +"");
		properties.setProperty(prefix+"ft_tolerance_below_near",            this.ft_tolerance_below_near +"");
		properties.setProperty(prefix+"ft_tolerance_above_far",             this.ft_tolerance_above_far +"");
		properties.setProperty(prefix+"ft_tolerance_below_far",             this.ft_tolerance_below_far +"");
		properties.setProperty(prefix+"ft_hor_vert_overlap",                this.ft_hor_vert_overlap +"");
		properties.setProperty(prefix+"ft_used_companions",                 this.ft_used_companions +"");
		properties.setProperty(prefix+"ft_used_true_companions",            this.ft_used_true_companions +"");

		properties.setProperty(prefix+"plPreferDisparity",          this.plPreferDisparity+"");
		properties.setProperty(prefix+"plDispNorm",                 this.plDispNorm +"");
		properties.setProperty(prefix+"plFrontoTol",                this.plFrontoTol +"");
		properties.setProperty(prefix+"plFrontoRms",                this.plFrontoRms +"");
		properties.setProperty(prefix+"plFrontoOffs",               this.plFrontoOffs +"");
		properties.setProperty(prefix+"PlFrontoPow",                this.PlFrontoPow +"");

		properties.setProperty(prefix+"plBlurBinVert",              this.plBlurBinVert +"");
		properties.setProperty(prefix+"plBlurBinHor",               this.plBlurBinHor +"");
		properties.setProperty(prefix+"plMaxDiffVert",              this.plMaxDiffVert +"");
		properties.setProperty(prefix+"plMaxDiffHor",               this.plMaxDiffHor +"");
		properties.setProperty(prefix+"plInitPasses",               this.plInitPasses+"");

		properties.setProperty(prefix+"plMinPoints",                this.plMinPoints+"");
		properties.setProperty(prefix+"plTargetEigen",              this.plTargetEigen +"");
		properties.setProperty(prefix+"plFractOutliers",            this.plFractOutliers +"");
		properties.setProperty(prefix+"plMaxOutliers",              this.plMaxOutliers+"");
		properties.setProperty(prefix+"plMinStrength",              this.plMinStrength +"");
		properties.setProperty(prefix+"plMaxEigen",                 this.plMaxEigen +"");
		properties.setProperty(prefix+"plEigenFloor",               this.plEigenFloor +"");
		properties.setProperty(prefix+"plEigenStick",               this.plEigenStick +"");
		properties.setProperty(prefix+"plBadPlate",                 this.plBadPlate +"");
		properties.setProperty(prefix+"plDbgMerge",                 this.plDbgMerge+"");
		properties.setProperty(prefix+"plWorstWorsening",           this.plWorstWorsening +"");
		properties.setProperty(prefix+"plWorstWorsening2",          this.plWorstWorsening2 +"");
		properties.setProperty(prefix+"plWorstEq",                  this.plWorstEq +"");
		properties.setProperty(prefix+"plWorstEq2",                 this.plWorstEq2 +"");
		properties.setProperty(prefix+"plOKMergeEigen",             this.plOKMergeEigen +"");
		properties.setProperty(prefix+"plMaxWorldSin2",             this.plMaxWorldSin2 +"");
		properties.setProperty(prefix+"pl2dForSin",                 this.pl2dForSin +"");
		properties.setProperty(prefix+"plWeakWorsening",            this.plWeakWorsening +"");
		properties.setProperty(prefix+"plMaxOverlap",               this.plMaxOverlap +"");

		properties.setProperty(prefix+"plWeakWeight",               this.plWeakWeight +"");
		properties.setProperty(prefix+"plWeakEigen",                this.plWeakEigen +"");
		properties.setProperty(prefix+"plWeakWeight2",              this.plWeakWeight2 +"");
		properties.setProperty(prefix+"plWeakEigen2",               this.plWeakEigen2 +"");
		properties.setProperty(prefix+"plSumThick",                 this.plSumThick +"");
		properties.setProperty(prefix+"plNeNeibCost",               this.plNeNeibCost +"");
		properties.setProperty(prefix+"plNeOwn",                    this.plNeOwn +"");

		properties.setProperty(prefix+"plExNeibCost",               this.plExNeibCost +"");
		properties.setProperty(prefix+"plExNeibSmooth",             this.plExNeibSmooth +"");
		properties.setProperty(prefix+"plMergeCostStar",            this.plMergeCostStar +"");
		properties.setProperty(prefix+"plMergeCost",                this.plMergeCost +"");

		properties.setProperty(prefix+"plConflMerge",               this.plConflMerge+"");
		properties.setProperty(prefix+"plConflRelax",               this.plConflRelax +"");
		properties.setProperty(prefix+"plConflSngl",                this.plConflSngl+"");
		properties.setProperty(prefix+"plConflSnglPair",            this.plConflSnglPair+"");

		properties.setProperty(prefix+"plWeakFgStrength",           this.plWeakFgStrength +"");
		properties.setProperty(prefix+"plWeakFgOutliers",           this.plWeakFgOutliers+"");
		properties.setProperty(prefix+"plWeakFgRelax",              this.plWeakFgRelax +"");

		properties.setProperty(prefix+"plThickWorld",               this.plThickWorld +"");
		properties.setProperty(prefix+"plThickWorldConfl",          this.plThickWorldConfl +"");
		properties.setProperty(prefix+"plRelaxComplete",            this.plRelaxComplete +"");
		properties.setProperty(prefix+"plRelaxComplete2",           this.plRelaxComplete2 +"");

		properties.setProperty(prefix+"plMaxZRatio",                this.plMaxZRatio +"");
		properties.setProperty(prefix+"plMaxDisp",                  this.plMaxDisp +"");
		properties.setProperty(prefix+"plCutTail",                  this.plCutTail +"");
		properties.setProperty(prefix+"plMinTail",                  this.plMinTail +"");

		properties.setProperty(prefix+"plDiscrEn",                  this.plDiscrEn+"");
		properties.setProperty(prefix+"plDiscrTolerance",           this.plDiscrTolerance +"");
		properties.setProperty(prefix+"plDiscrDispRange",           this.plDiscrDispRange +"");
		properties.setProperty(prefix+"plDiscrSteps",               this.plDiscrSteps+"");
		//  			properties.setProperty(prefix+"plDiscrVariants",            this.plDiscrVariants+"");
		properties.setProperty(prefix+"plDiscrMode",                this.plDiscrMode+"");
		properties.setProperty(prefix+"plDiscrVarFloor",            this.plDiscrVarFloor +"");
		properties.setProperty(prefix+"plDiscrSigma",               this.plDiscrSigma +"");
		properties.setProperty(prefix+"plDiscrBlur",                this.plDiscrBlur +"");
		properties.setProperty(prefix+"plDiscrExclusivity",          this.plDiscrExclusivity +"");
		properties.setProperty(prefix+"plDiscrExclus2",             this.plDiscrExclus2 +"");
		properties.setProperty(prefix+"plDiscrStrict",              this.plDiscrStrict+"");
		properties.setProperty(prefix+"plDiscrCorrMax",             this.plDiscrCorrMax +"");
		properties.setProperty(prefix+"plDiscrCorrMerge",           this.plDiscrCorrMerge +"");
		properties.setProperty(prefix+"plDiscrSteal",               this.plDiscrSteal+"");
		properties.setProperty(prefix+"plDiscrGrown",               this.plDiscrGrown+"");
		properties.setProperty(prefix+"plDiscrXMedian",             this.plDiscrXMedian +"");

		properties.setProperty(prefix+"plCostDist",                 this.plCostDist +"");
		properties.setProperty(prefix+"plCostKrq",                  this.plCostKrq +"");
		properties.setProperty(prefix+"plCostKrqEq",                this.plCostKrqEq +"");
		properties.setProperty(prefix+"plCostWrq",                  this.plCostWrq +"");
		properties.setProperty(prefix+"plCostWrqEq",                this.plCostWrqEq +"");
		properties.setProperty(prefix+"plCostSin2",                 this.plCostSin2 +"");
		properties.setProperty(prefix+"plCostRdist2",               this.plCostRdist2 +"");

		properties.setProperty(prefix+"plConflDualTri",             this.plConflDualTri+"");
		properties.setProperty(prefix+"plConflMulti",               this.plConflMulti+"");
		properties.setProperty(prefix+"plConflDiag",                this.plConflDiag+"");
		properties.setProperty(prefix+"plConflStar",                this.plConflStar+"");
		properties.setProperty(prefix+"plStarSteps",                this.plStarSteps+"");
		properties.setProperty(prefix+"plStarOrtho",                this.plStarOrtho +"");
		properties.setProperty(prefix+"plStarDiag",                 this.plStarDiag +"");
		properties.setProperty(prefix+"plStarPwr",                  this.plStarPwr +"");
		properties.setProperty(prefix+"plStarWeightPwr",            this.plStarWeightPwr +"");
		properties.setProperty(prefix+"plWeightToDens",             this.plWeightToDens +"");
		properties.setProperty(prefix+"plStarValPwr",               this.plStarValPwr +"");
		properties.setProperty(prefix+"plDblTriLoss",               this.plDblTriLoss +"");
		properties.setProperty(prefix+"plNewConfl",                 this.plNewConfl+"");
		properties.setProperty(prefix+"plMaxChanges",               this.plMaxChanges+"");

		properties.setProperty(prefix+"plMutualOnly",               this.plMutualOnly+"");
		properties.setProperty(prefix+"plFillSquares",              this.plFillSquares+"");
		properties.setProperty(prefix+"plCutCorners",               this.plCutCorners+"");
		properties.setProperty(prefix+"plHypotenuse",               this.plHypotenuse+"");
		properties.setProperty(prefix+"plPull",                     this.plPull +"");
		properties.setProperty(prefix+"plNormPow",                  this.plNormPow +"");
		properties.setProperty(prefix+"plIterations",               this.plIterations+"");
		properties.setProperty(prefix+"plStopBad",                  this.plStopBad+"");
		properties.setProperty(prefix+"plPrecision",                this.plPrecision+"");

		properties.setProperty(prefix+"plSplitPull",                this.plSplitPull +"");
		properties.setProperty(prefix+"plSplitMinNeib",             this.plSplitMinNeib+"");
		properties.setProperty(prefix+"plSplitMinWeight",           this.plSplitMinWeight +"");
		properties.setProperty(prefix+"plSplitMinQuality",          this.plSplitMinQuality +"");
		properties.setProperty(prefix+"plSplitApply",               this.plSplitApply+"");
		properties.setProperty(prefix+"plNonExclusive",             this.plNonExclusive+"");
		properties.setProperty(prefix+"plUseOtherPlanes",           this.plUseOtherPlanes+"");
		properties.setProperty(prefix+"plAllowParallel",            this.plAllowParallel+"");
		properties.setProperty(prefix+"plMaxDiff",                  this.plMaxDiff +"");
		properties.setProperty(prefix+"plOtherDiff",                this.plOtherDiff +"");
		properties.setProperty(prefix+"plSplitXY",                  this.plSplitXY+"");
		properties.setProperty(prefix+"plSplitXYTolerance",          this.plSplitXYTolerance +"");

		properties.setProperty(prefix+"plFuse",                     this.plFuse+"");
		properties.setProperty(prefix+"plKeepOrphans",              this.plKeepOrphans+"");
		properties.setProperty(prefix+"plMinOrphan",                this.plMinOrphan +"");

		properties.setProperty(prefix+"plSnapDispAny",              this.plSnapDispAny +"");
		properties.setProperty(prefix+"plSnapStrengthAny",          this.plSnapStrengthAny +"");
		properties.setProperty(prefix+"plSnapNegAny",               this.plSnapNegAny +"");
		properties.setProperty(prefix+"plSnapDispMax",              this.plSnapDispMax +"");
		properties.setProperty(prefix+"plSnapDispWeight",           this.plSnapDispWeight +"");
		properties.setProperty(prefix+"plSnapZeroMode",             this.plSnapZeroMode+"");

		properties.setProperty(prefix+"msUseSel",                   this.msUseSel+"");
		properties.setProperty(prefix+"msDivideByArea",             this.msDivideByArea+"");
		properties.setProperty(prefix+"msScaleProj",                this.msScaleProj +"");
		properties.setProperty(prefix+"msFractUni",                 this.msFractUni +"");

		properties.setProperty(prefix+"tsNoEdge",                   this.tsNoEdge+"");
		properties.setProperty(prefix+"tsUseCenter",                this.tsUseCenter+"");
		properties.setProperty(prefix+"tsMaxDiff",                  this.tsMaxDiff +"");
		properties.setProperty(prefix+"tsMinDiffOther",             this.tsMinDiffOther +"");
		properties.setProperty(prefix+"tsMinStrength",              this.tsMinStrength +"");
		properties.setProperty(prefix+"tsMaxStrength",              this.tsMaxStrength +"");
		properties.setProperty(prefix+"tsMinSurface",               this.tsMinSurface +"");
		properties.setProperty(prefix+"tsMoveDirs",                 this.tsMoveDirs+"");
		properties.setProperty(prefix+"tsSurfStrPow",               this.tsSurfStrPow +"");
		properties.setProperty(prefix+"tsAddStrength",              this.tsAddStrength +"");
		properties.setProperty(prefix+"tsSigma",                    this.tsSigma +"");
		properties.setProperty(prefix+"tsNSigma",                   this.tsNSigma +"");
		properties.setProperty(prefix+"tsMinPull",                  this.tsMinPull +"");
		properties.setProperty(prefix+"tsMinAdvantage",             this.tsMinAdvantage +"");

		properties.setProperty(prefix+"tsClustSize",                this.tsClustSize +"");
		properties.setProperty(prefix+"tsClustWeight",              this.tsClustWeight +"");
		properties.setProperty(prefix+"tsMinNeib",                  this.tsMinNeib +"");
		properties.setProperty(prefix+"tsMaxSurStrength",           this.tsMaxSurStrength +"");
		properties.setProperty(prefix+"tsCountDis",                 this.tsCountDis +"");

		properties.setProperty(prefix+"tsEnPlaneSeed",              this.tsEnPlaneSeed+"");
		properties.setProperty(prefix+"tsEnOnly",                   this.tsEnOnly+"");
		properties.setProperty(prefix+"tsEnGrow",                   this.tsEnGrow+"");
		properties.setProperty(prefix+"tsGrowStrength",             this.tsGrowStrength +"");
		properties.setProperty(prefix+"tsGrowStrong",               this.tsGrowStrong+"");
		properties.setProperty(prefix+"tsContStrength",             this.tsContStrength +"");
		properties.setProperty(prefix+"tsContDiff",                 this.tsContDiff +"");

		properties.setProperty(prefix+"tsEnSingle",                 this.tsEnSingle+"");
		properties.setProperty(prefix+"tsEnMulti",                  this.tsEnMulti+"");
		properties.setProperty(prefix+"tsRemoveWeak1",              this.tsRemoveWeak1+"");
		properties.setProperty(prefix+"tsGrowSurround",             this.tsGrowSurround+"");
		properties.setProperty(prefix+"tsRemoveWeak2",              this.tsRemoveWeak2+"");
		properties.setProperty(prefix+"tsLoopMulti",                this.tsLoopMulti+"");
		properties.setProperty(prefix+"tsShow",                     this.tsShow+"");
		properties.setProperty(prefix+"tsNumClust",                 this.tsNumClust +"");

		properties.setProperty(prefix+"tsConsensMode",              this.tsConsensMode +"");
		properties.setProperty(prefix+"tsConsensAgree",             this.tsConsensAgree +"");

		properties.setProperty(prefix+"taMinFgBg",                  this.taMinFgBg +"");
		properties.setProperty(prefix+"taMinFgEdge",                this.taMinFgEdge +"");
		properties.setProperty(prefix+"taMinColSep",                this.taMinColSep +"");
		properties.setProperty(prefix+"taMinColDiff",               this.taMinColDiff +"");
		properties.setProperty(prefix+"taOutlier",                  this.taOutlier +"");
		properties.setProperty(prefix+"taDiffPwr",                  this.taDiffPwr +"");
		properties.setProperty(prefix+"taBestPwr",                  this.taBestPwr +"");
		properties.setProperty(prefix+"taDiff9Pwr",                 this.taDiff9Pwr +"");
		properties.setProperty(prefix+"taColSigma",                 this.taColSigma +"");
		properties.setProperty(prefix+"taColFraction",              this.taColFraction +"");

		properties.setProperty(prefix+"taCostEmpty",                this.taCostEmpty +"");
		properties.setProperty(prefix+"taCostNoLink",               this.taCostNoLink +"");
		properties.setProperty(prefix+"taCostSwitch",               this.taCostSwitch +"");
		properties.setProperty(prefix+"taCostColor",                this.taCostColor +"");
		properties.setProperty(prefix+"taCostDiff",                 this.taCostDiff +"");
		properties.setProperty(prefix+"taCostDiffBest",             this.taCostDiffBest +"");
		properties.setProperty(prefix+"taCostDiff9",                this.taCostDiff9 +"");
		properties.setProperty(prefix+"taCostWeakFgnd",             this.taCostWeakFgnd +"");
		properties.setProperty(prefix+"taCostFlaps",                this.taCostFlaps +"");
		properties.setProperty(prefix+"taCostMismatch",             this.taCostMismatch +"");

		properties.setProperty(prefix+"taEnEmpty",                  this.taEnEmpty +"");
		properties.setProperty(prefix+"taEnNoLink",                 this.taEnNoLink +"");
		properties.setProperty(prefix+"taEnSwitch",                 this.taEnSwitch +"");
		properties.setProperty(prefix+"taEnColor",                  this.taEnColor +"");
		properties.setProperty(prefix+"taEnDiff",                   this.taEnDiff +"");
		properties.setProperty(prefix+"taEnDiffBest",               this.taEnDiffBest +"");
		properties.setProperty(prefix+"taEnDiff9",                  this.taEnDiff9 +"");
		properties.setProperty(prefix+"taEnWeakFgnd",               this.taEnWeakFgnd +"");
		properties.setProperty(prefix+"taEnFlaps",                  this.taEnFlaps +"");
		properties.setProperty(prefix+"taEnMismatch",               this.taEnMismatch +"");


		properties.setProperty(prefix+"gpu_corr_rad",               this.gpu_corr_rad +"");
		properties.setProperty(prefix+"gpu_weight_r",               this.gpu_weight_r +"");
		properties.setProperty(prefix+"gpu_weight_b",               this.gpu_weight_b +"");
		properties.setProperty(prefix+"gpu_sigma_r",                this.gpu_sigma_r +"");
		properties.setProperty(prefix+"gpu_sigma_b",                this.gpu_sigma_b +"");
		properties.setProperty(prefix+"gpu_sigma_g",                this.gpu_sigma_g +"");
		properties.setProperty(prefix+"gpu_sigma_m",                this.gpu_sigma_m +"");
		properties.setProperty(prefix+"gpu_sigma_rb_corr",          this.gpu_sigma_rb_corr +"");
		properties.setProperty(prefix+"gpu_sigma_corr",             this.gpu_sigma_corr +"");
		properties.setProperty(prefix+"gpu_sigma_corr_m",           this.gpu_sigma_corr_m +"");
		properties.setProperty(prefix+"gpu_fatz",                   this.gpu_fatz +"");
		properties.setProperty(prefix+"gpu_fatz_m",                 this.gpu_fatz_m +"");

		properties.setProperty(prefix+"gpu_woi",                    this.gpu_woi +"");
		properties.setProperty(prefix+"gpu_woi_tx",                 this.gpu_woi_tx +"");
		properties.setProperty(prefix+"gpu_woi_ty",                 this.gpu_woi_ty +"");
		properties.setProperty(prefix+"gpu_woi_twidth",             this.gpu_woi_twidth +"");
		properties.setProperty(prefix+"gpu_woi_theight",            this.gpu_woi_theight +"");
		properties.setProperty(prefix+"gpu_woi_round",              this.gpu_woi_round +"");
		properties.setProperty(prefix+"gpu_save_ports_xy",          this.gpu_save_ports_xy +"");
		properties.setProperty(prefix+"gpu_show_jtextures",         this.gpu_show_jtextures +"");
		properties.setProperty(prefix+"gpu_show_extra",             this.gpu_show_extra +"");

		properties.setProperty(prefix+"gpu_use_main",               this.gpu_use_main +"");
		properties.setProperty(prefix+"gpu_use_main_macro",         this.gpu_use_main_macro +"");
		properties.setProperty(prefix+"gpu_use_main_adjust",        this.gpu_use_main_adjust +"");
		properties.setProperty(prefix+"gpu_use_aux",                this.gpu_use_aux +"");
		properties.setProperty(prefix+"gpu_use_aux_macro",          this.gpu_use_aux_macro +"");
		properties.setProperty(prefix+"gpu_use_aux_adjust",         this.gpu_use_aux_adjust +"");
		
		properties.setProperty(prefix+"debug_initial_discriminate",       this.debug_initial_discriminate+"");
		properties.setProperty(prefix+"dbg_migrate",                      this.dbg_migrate+"");

		properties.setProperty(prefix+"dbg_early_exit",                   this.dbg_early_exit+"");
		properties.setProperty(prefix+"show_first_bg",                    this.show_first_bg+"");

		properties.setProperty(prefix+"show_extrinsic",                   this.show_extrinsic+"");
		properties.setProperty(prefix+"show_ortho_combine",               this.show_ortho_combine+"");
		properties.setProperty(prefix+"show_refine_supertiles",           this.show_refine_supertiles+"");
		properties.setProperty(prefix+"show_bgnd_nonbgnd",                this.show_bgnd_nonbgnd+"");
		properties.setProperty(prefix+"show_filter_scan",                 this.show_filter_scan+"");
		properties.setProperty(prefix+"show_combined",                    this.show_combined+"");
		properties.setProperty(prefix+"show_unique",                      this.show_unique+"");
		properties.setProperty(prefix+"show_histograms",                  this.show_histograms+"");
		properties.setProperty(prefix+"show_init_refine",                 this.show_init_refine+"");
		properties.setProperty(prefix+"show_expand",                      this.show_expand+"");
		properties.setProperty(prefix+"show_variant",                     this.show_variant+"");
		properties.setProperty(prefix+"show_retry_far",                   this.show_retry_far+"");
		properties.setProperty(prefix+"show_macro",                       this.show_macro+"");
		properties.setProperty(prefix+"show_shells",                      this.show_shells+"");
		properties.setProperty(prefix+"show_neighbors",                   this.show_neighbors+"");
		properties.setProperty(prefix+"show_flaps_dirs",                  this.show_flaps_dirs+"");
		properties.setProperty(prefix+"show_first_clusters",              this.show_first_clusters+"");
		properties.setProperty(prefix+"show_planes",                      this.show_planes+"");

		properties.setProperty(prefix+"vertical_xyz.x",                   this.vertical_xyz[0]+"");
		properties.setProperty(prefix+"vertical_xyz.y",                   this.vertical_xyz[1]+"");
		properties.setProperty(prefix+"vertical_xyz.z",                   this.vertical_xyz[2]+"");

		if (z_corr_map != null) {
			for (HashMap.Entry<String,Double> entry : z_corr_map.entrySet()){
				properties.setProperty(prefix+Z_CORR_PREFIX+entry.getKey(),   entry.getValue().toString());
			}
		}
		setPropertiesInfinityDistance(prefix, properties);
		/*
			if (infinity_distace_map != null) {
				for (HashMap.Entry<String,Double> entry : infinity_distace_map.entrySet()){
					properties.setProperty(prefix+INFINITY_DISTANCE_PREFIX+entry.getKey(),   entry.getValue().toString());
				}
			}
		 */
		img_dtt.setProperties (prefix+"_img_dtt", properties);
		mlfp.setProperties    (prefix+"_mlfp",    properties);
		rig.setProperties     (prefix+"_rig",     properties);
		poles.setProperties   (prefix+"_poles",   properties);
		lwir.setProperties    (prefix+"_lwir",   properties);
	}

	public void setPropertiesInfinityDistance(String prefix,Properties properties){
		if (infinity_distace_map != null) {
			for (HashMap.Entry<String,Double> entry : infinity_distace_map.entrySet()){
				properties.setProperty(prefix+INFINITY_DISTANCE_PREFIX+entry.getKey(),   entry.getValue().toString());
			}
		}
	}


	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"transform_size")!=null)                this.transform_size=Integer.parseInt(properties.getProperty(prefix+"transform_size"));
		if (properties.getProperty(prefix+"clt_window")!=null)                    this.clt_window=Integer.parseInt(properties.getProperty(prefix+"clt_window"));
		if (properties.getProperty(prefix+"shift_x")!=null)                       this.shift_x=Double.parseDouble(properties.getProperty(prefix+"shift_x"));
		if (properties.getProperty(prefix+"shift_y")!=null)                       this.shift_y=Double.parseDouble(properties.getProperty(prefix+"shift_y"));

		if (properties.getProperty(prefix+"tileStep")!=null)                      this.tileStep=Integer.parseInt(properties.getProperty(prefix+"tileStep"));
		if (properties.getProperty(prefix+"iclt_mask")!=null)                     this.iclt_mask=Integer.parseInt(properties.getProperty(prefix+"iclt_mask"));
		if (properties.getProperty(prefix+"tileX")!=null)                         this.tileX=Integer.parseInt(properties.getProperty(prefix+"tileX"));
		if (properties.getProperty(prefix+"tileY")!=null)                         this.tileY=Integer.parseInt(properties.getProperty(prefix+"tileY"));
		if (properties.getProperty(prefix+"dbg_mode")!=null)                      this.dbg_mode=Integer.parseInt(properties.getProperty(prefix+"dbg_mode"));
		if (properties.getProperty(prefix+"ishift_x")!=null)                      this.ishift_x=Integer.parseInt(properties.getProperty(prefix+"ishift_x"));
		if (properties.getProperty(prefix+"ishift_y")!=null)                      this.ishift_y=Integer.parseInt(properties.getProperty(prefix+"ishift_y"));
		if (properties.getProperty(prefix+"fat_zero")!=null)                      this.fat_zero=Double.parseDouble(properties.getProperty(prefix+"fat_zero"));
		if (properties.getProperty(prefix+"fat_zero_mono")!=null)                 this.fat_zero_mono=Double.parseDouble(properties.getProperty(prefix+"fat_zero_mono"));
		if (properties.getProperty(prefix+"corr_sigma")!=null)                    this.corr_sigma=Double.parseDouble(properties.getProperty(prefix+"corr_sigma"));
		if (properties.getProperty(prefix+"corr_sigma_mono")!=null)               this.corr_sigma_mono=Double.parseDouble(properties.getProperty(prefix+"corr_sigma_mono"));
		if (properties.getProperty(prefix+"scale_strength_main")!=null)           this.scale_strength_main=Double.parseDouble(properties.getProperty(prefix+"scale_strength_main"));
		if (properties.getProperty(prefix+"scale_strength_aux")!=null)            this.scale_strength_aux=Double.parseDouble(properties.getProperty(prefix+"scale_strength_aux"));

		if (properties.getProperty(prefix+"norm_kern")!=null)                     this.norm_kern=Boolean.parseBoolean(properties.getProperty(prefix+"norm_kern"));
		if (properties.getProperty(prefix+"gain_equalize")!=null)                 this.gain_equalize=Boolean.parseBoolean(properties.getProperty(prefix+"gain_equalize"));
		if (properties.getProperty(prefix+"colors_equalize")!=null)               this.colors_equalize=Boolean.parseBoolean(properties.getProperty(prefix+"colors_equalize"));
		if (properties.getProperty(prefix+"nosat_equalize")!=null)                this.nosat_equalize=Boolean.parseBoolean(properties.getProperty(prefix+"nosat_equalize"));
		if (properties.getProperty(prefix+"sat_level")!=null)                     this.sat_level=Double.parseDouble(properties.getProperty(prefix+"sat_level"));
		if (properties.getProperty(prefix+"max_overexposure")!=null)              this.max_overexposure=Double.parseDouble(properties.getProperty(prefix+"max_overexposure"));

		if (properties.getProperty(prefix+"novignetting_r")!=null)                this.novignetting_r=Double.parseDouble(properties.getProperty(prefix+"novignetting_r"));
		if (properties.getProperty(prefix+"novignetting_g")!=null)                this.novignetting_g=Double.parseDouble(properties.getProperty(prefix+"novignetting_g"));
		if (properties.getProperty(prefix+"novignetting_b")!=null)                this.novignetting_b=Double.parseDouble(properties.getProperty(prefix+"novignetting_b"));
		if (properties.getProperty(prefix+"scale_r")!=null)                       this.scale_r=Double.parseDouble(properties.getProperty(prefix+"scale_r"));
		if (properties.getProperty(prefix+"scale_g")!=null)                       this.scale_g=Double.parseDouble(properties.getProperty(prefix+"scale_g"));
		if (properties.getProperty(prefix+"scale_b")!=null)                       this.scale_b=Double.parseDouble(properties.getProperty(prefix+"scale_b"));
		if (properties.getProperty(prefix+"vignetting_max")!=null)                this.vignetting_max=Double.parseDouble(properties.getProperty(prefix+"vignetting_max"));
		if (properties.getProperty(prefix+"vignetting_range")!=null)              this.vignetting_range=Double.parseDouble(properties.getProperty(prefix+"vignetting_range"));
		if (properties.getProperty(prefix+"kernel_step")!=null)                   this.kernel_step=Integer.parseInt(properties.getProperty(prefix+"kernel_step"));
		if (properties.getProperty(prefix+"disparity")!=null)                     this.disparity=Double.parseDouble(properties.getProperty(prefix+"disparity"));
		if (properties.getProperty(prefix+"z_correction")!=null)                  this.z_correction=Double.parseDouble(properties.getProperty(prefix+"z_correction"));
		if (properties.getProperty(prefix+"correlate")!=null)                     this.correlate=Boolean.parseBoolean(properties.getProperty(prefix+"correlate"));
		if (properties.getProperty(prefix+"corr_mask")!=null)                     this.corr_mask=Integer.parseInt(properties.getProperty(prefix+"corr_mask"));
		if (properties.getProperty(prefix+"corr_sym")!=null)                      this.corr_sym=Boolean.parseBoolean(properties.getProperty(prefix+"corr_sym"));
		if (properties.getProperty(prefix+"corr_keep")!=null)                     this.corr_keep=Boolean.parseBoolean(properties.getProperty(prefix+"corr_keep"));
		if (properties.getProperty(prefix+"corr_show")!=null)                     this.corr_show=Boolean.parseBoolean(properties.getProperty(prefix+"corr_show"));
		if (properties.getProperty(prefix+"corr_mismatch")!=null)                 this.corr_mismatch=Boolean.parseBoolean(properties.getProperty(prefix+"corr_mismatch"));
		if (properties.getProperty(prefix+"corr_offset")!=null)                   this.corr_offset=Double.parseDouble(properties.getProperty(prefix+"corr_offset"));
		if (properties.getProperty(prefix+"corr_red")!=null)                      this.corr_red=Double.parseDouble(properties.getProperty(prefix+"corr_red"));
		if (properties.getProperty(prefix+"corr_blue")!=null)                     this.corr_blue=Double.parseDouble(properties.getProperty(prefix+"corr_blue"));
		if (properties.getProperty(prefix+"corr_normalize")!=null)                this.corr_normalize=Boolean.parseBoolean(properties.getProperty(prefix+"corr_normalize"));
		if (properties.getProperty(prefix+"min_corr")!=null)                      this.min_corr=Double.parseDouble(properties.getProperty(prefix+"min_corr"));
		if (properties.getProperty(prefix+"min_corr_normalized")!=null)           this.min_corr_normalized=Double.parseDouble(properties.getProperty(prefix+"min_corr_normalized"));
		if (properties.getProperty(prefix+"max_corr_sigma")!=null)                this.max_corr_sigma=Double.parseDouble(properties.getProperty(prefix+"max_corr_sigma"));
		if (properties.getProperty(prefix+"max_corr_radius")!=null)               this.max_corr_radius=Double.parseDouble(properties.getProperty(prefix+"max_corr_radius"));

		// for compatibility with old settings
		if (properties.getProperty(prefix+"enhortho_width")!=null)                this.img_dtt.setEnhOrthoWidth(Integer.parseInt(properties.getProperty(prefix+"enhortho_width")));
		if (properties.getProperty(prefix+"enhortho_scale")!=null)                this.img_dtt.setEnhOrthoScale(Double.parseDouble(properties.getProperty(prefix+"enhortho_scale")));
		// for compatibility with old settings

		if (properties.getProperty(prefix+"max_corr_double")!=null)               this.max_corr_double=Boolean.parseBoolean(properties.getProperty(prefix+"max_corr_double"));
		if (properties.getProperty(prefix+"corr_mode")!=null)                     this.corr_mode=Integer.parseInt(properties.getProperty(prefix+"corr_mode"));
		if (properties.getProperty(prefix+"corr_border_contrast")!=null)          this.corr_border_contrast=Double.parseDouble(properties.getProperty(prefix+"corr_border_contrast"));
		if (properties.getProperty(prefix+"tile_task_op")!=null)                  this.tile_task_op=Integer.parseInt(properties.getProperty(prefix+"tile_task_op"));
		if (properties.getProperty(prefix+"tile_task_wl")!=null)                  this.tile_task_wl=Integer.parseInt(properties.getProperty(prefix+"tile_task_wl"));
		if (properties.getProperty(prefix+"tile_task_wt")!=null)                  this.tile_task_wt=Integer.parseInt(properties.getProperty(prefix+"tile_task_wt"));
		if (properties.getProperty(prefix+"tile_task_ww")!=null)                  this.tile_task_ww=Integer.parseInt(properties.getProperty(prefix+"tile_task_ww"));
		if (properties.getProperty(prefix+"tile_task_wh")!=null)                  this.tile_task_wh=Integer.parseInt(properties.getProperty(prefix+"tile_task_wh"));
		if (properties.getProperty(prefix+"min_shot")!=null)                      this.min_shot=Double.parseDouble(properties.getProperty(prefix+"min_shot"));
		if (properties.getProperty(prefix+"scale_shot")!=null)                    this.scale_shot=Double.parseDouble(properties.getProperty(prefix+"scale_shot"));
		if (properties.getProperty(prefix+"diff_sigma")!=null)                    this.diff_sigma=Double.parseDouble(properties.getProperty(prefix+"diff_sigma"));
		if (properties.getProperty(prefix+"diff_threshold")!=null)                this.diff_threshold=Double.parseDouble(properties.getProperty(prefix+"diff_threshold"));
		if (properties.getProperty(prefix+"diff_gauss")!=null)                    this.diff_gauss=Boolean.parseBoolean(properties.getProperty(prefix+"diff_gauss"));
		if (properties.getProperty(prefix+"min_agree")!=null)                     this.min_agree=Double.parseDouble(properties.getProperty(prefix+"min_agree"));
		if (properties.getProperty(prefix+"dust_remove")!=null)                   this.dust_remove=Boolean.parseBoolean(properties.getProperty(prefix+"dust_remove"));
		if (properties.getProperty(prefix+"black_back")!=null)                    this.black_back=Boolean.parseBoolean(properties.getProperty(prefix+"black_back"));
		if (properties.getProperty(prefix+"keep_weights")!=null)                  this.keep_weights=Boolean.parseBoolean(properties.getProperty(prefix+"keep_weights"));
		if (properties.getProperty(prefix+"sharp_alpha")!=null)                   this.sharp_alpha=Boolean.parseBoolean(properties.getProperty(prefix+"sharp_alpha"));
		if (properties.getProperty(prefix+"alpha0")!=null)                        this.alpha0=Double.parseDouble(properties.getProperty(prefix+"alpha0"));
		if (properties.getProperty(prefix+"alpha1")!=null)                        this.alpha1=Double.parseDouble(properties.getProperty(prefix+"alpha1"));
		if (properties.getProperty(prefix+"gen_chn_stacks")!=null)                this.gen_chn_stacks=Boolean.parseBoolean(properties.getProperty(prefix+"gen_chn_stacks"));
		if (properties.getProperty(prefix+"gen_chn_img")!=null)                   this.gen_chn_img=Boolean.parseBoolean(properties.getProperty(prefix+"gen_chn_img"));
		if (properties.getProperty(prefix+"gen_4_img")!=null)                     this.gen_4_img=Boolean.parseBoolean(properties.getProperty(prefix+"gen_4_img"));
		if (properties.getProperty(prefix+"show_nonoverlap")!=null)               this.show_nonoverlap=Boolean.parseBoolean(properties.getProperty(prefix+"show_nonoverlap"));
		if (properties.getProperty(prefix+"show_overlap")!=null)                  this.show_overlap=Boolean.parseBoolean(properties.getProperty(prefix+"show_overlap"));
		if (properties.getProperty(prefix+"show_rgba_color")!=null)               this.show_rgba_color=Boolean.parseBoolean(properties.getProperty(prefix+"show_rgba_color"));
		if (properties.getProperty(prefix+"show_map")!=null)                      this.show_map=Boolean.parseBoolean(properties.getProperty(prefix+"show_map"));
		if (properties.getProperty(prefix+"show_corr")!=null)                     this.show_corr=Boolean.parseBoolean(properties.getProperty(prefix+"show_corr"));

		if (properties.getProperty(prefix+"disp_scan_start")!=null)               this.disp_scan_start=Double.parseDouble(properties.getProperty(prefix+"disp_scan_start"));
		if (properties.getProperty(prefix+"disp_scan_step")!=null)                this.disp_scan_step=Double.parseDouble(properties.getProperty(prefix+"disp_scan_step"));
		if (properties.getProperty(prefix+"disp_scan_count")!=null)               this.disp_scan_count=Integer.parseInt(properties.getProperty(prefix+"disp_scan_count"));

		if (properties.getProperty(prefix+"fine_dbg")!=null)                      this.fine_dbg=Boolean.parseBoolean(properties.getProperty(prefix+"fine_dbg"));
		if (properties.getProperty(prefix+"fine_corr_x_0")!=null)                 this.fine_corr_x_0=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_0"));
		if (properties.getProperty(prefix+"fine_corr_y_0")!=null)                 this.fine_corr_y_0=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_0"));
		if (properties.getProperty(prefix+"fine_corr_x_1")!=null)                 this.fine_corr_x_1=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_1"));
		if (properties.getProperty(prefix+"fine_corr_y_1")!=null)                 this.fine_corr_y_1=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_1"));
		if (properties.getProperty(prefix+"fine_corr_x_2")!=null)                 this.fine_corr_x_2=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_2"));
		if (properties.getProperty(prefix+"fine_corr_y_2")!=null)                 this.fine_corr_y_2=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_2"));
		if (properties.getProperty(prefix+"fine_corr_x_3")!=null)                 this.fine_corr_x_3=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_3"));
		if (properties.getProperty(prefix+"fine_corr_y_3")!=null)                 this.fine_corr_y_3=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_3"));
		if (properties.getProperty(prefix+"fine_corr_ignore")!=null)              this.fine_corr_ignore=Boolean.parseBoolean(properties.getProperty(prefix+"fine_corr_ignore"));
		if (properties.getProperty(prefix+"fine_corr_apply")!=null)               this.fine_corr_apply=Boolean.parseBoolean(properties.getProperty(prefix+"fine_corr_apply"));

		if (properties.getProperty(prefix+"fcorr_radius")!=null)                  this.fcorr_radius=Double.parseDouble(properties.getProperty(prefix+"fcorr_radius"));
		if (properties.getProperty(prefix+"fcorr_min_strength")!=null)            this.fcorr_min_strength=Double.parseDouble(properties.getProperty(prefix+"fcorr_min_strength"));
		if (properties.getProperty(prefix+"fcorr_disp_diff")!=null)               this.fcorr_disp_diff=Double.parseDouble(properties.getProperty(prefix+"fcorr_disp_diff"));
		if (properties.getProperty(prefix+"fcorr_quadratic")!=null)               this.fcorr_quadratic=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_quadratic"));
		if (properties.getProperty(prefix+"fcorr_ignore")!=null)                  this.fcorr_ignore=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_ignore"));

		if (properties.getProperty(prefix+"fcorr_inf_strength")!=null)            this.fcorr_inf_strength=Double.parseDouble(properties.getProperty(prefix+"fcorr_inf_strength"));
		if (properties.getProperty(prefix+"fcorr_inf_diff")!=null)                this.fcorr_inf_diff=Double.parseDouble(properties.getProperty(prefix+"fcorr_inf_diff"));
		if (properties.getProperty(prefix+"fcorr_inf_quad")!=null)                this.fcorr_inf_quad=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_inf_quad"));
		if (properties.getProperty(prefix+"fcorr_inf_vert")!=null)                this.fcorr_inf_vert=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_inf_vert"));



		if (properties.getProperty(prefix+"inf_disp_apply")!=null)                this.inf_disp_apply=Boolean.parseBoolean(properties.getProperty(prefix+"inf_disp_apply"));
		if (properties.getProperty(prefix+"inf_repeat")!=null)                    this.inf_repeat=Integer.parseInt(properties.getProperty(prefix+"inf_repeat"));
		//  			if (properties.getProperty(prefix+"inf_mism_apply")!=null)              this.inf_mism_apply=Boolean.parseBoolean(properties.getProperty(prefix+"inf_mism_apply"));

		if (properties.getProperty(prefix+"inf_iters")!=null)                     this.inf_iters=Integer.parseInt(properties.getProperty(prefix+"inf_iters"));
		if (properties.getProperty(prefix+"inf_final_diff")!=null)                this.inf_final_diff=Double.parseDouble(properties.getProperty(prefix+"inf_final_diff"));
		if (properties.getProperty(prefix+"inf_far_pull")!=null)                  this.inf_far_pull=Double.parseDouble(properties.getProperty(prefix+"inf_far_pull"));

		if (properties.getProperty(prefix+"inf_str_pow")!=null)                   this.inf_str_pow=Double.parseDouble(properties.getProperty(prefix+"inf_str_pow"));
		if (properties.getProperty(prefix+"inf_smpl_side")!=null)                 this.inf_smpl_side=Integer.parseInt(properties.getProperty(prefix+"inf_smpl_side"));
		if (properties.getProperty(prefix+"inf_smpl_num")!=null)                  this.inf_smpl_num=Integer.parseInt(properties.getProperty(prefix+"inf_smpl_num"));
		if (properties.getProperty(prefix+"inf_smpl_rms")!=null)                  this.inf_smpl_rms=Double.parseDouble(properties.getProperty(prefix+"inf_smpl_rms"));

		if (properties.getProperty(prefix+"ih_smpl_step")!=null)                  this.ih_smpl_step=Integer.parseInt(properties.getProperty(prefix+"ih_smpl_step"));
		if (properties.getProperty(prefix+"ih_disp_min")!=null)                   this.ih_disp_min=Double.parseDouble(properties.getProperty(prefix+"ih_disp_min"));
		if (properties.getProperty(prefix+"ih_disp_step")!=null)                  this.ih_disp_step=Double.parseDouble(properties.getProperty(prefix+"ih_disp_step"));
		if (properties.getProperty(prefix+"ih_num_bins")!=null)                   this.ih_num_bins=Integer.parseInt(properties.getProperty(prefix+"ih_num_bins"));
		if (properties.getProperty(prefix+"ih_sigma")!=null)                      this.ih_sigma=Double.parseDouble(properties.getProperty(prefix+"ih_sigma"));
		if (properties.getProperty(prefix+"ih_max_diff")!=null)                   this.ih_max_diff=Double.parseDouble(properties.getProperty(prefix+"ih_max_diff"));
		if (properties.getProperty(prefix+"ih_min_samples")!=null)                this.ih_min_samples=Integer.parseInt(properties.getProperty(prefix+"ih_min_samples"));
		if (properties.getProperty(prefix+"ih_norm_center")!=null)                this.ih_norm_center=Boolean.parseBoolean(properties.getProperty(prefix+"ih_norm_center"));
		if (properties.getProperty(prefix+"inf_restore_disp")!=null)              this.inf_restore_disp=Boolean.parseBoolean(properties.getProperty(prefix+"inf_restore_disp"));

		if (properties.getProperty(prefix+"ly_lma_ers")!=null)                    this.ly_lma_ers=Boolean.parseBoolean(properties.getProperty(prefix+"ly_lma_ers"));
		if (properties.getProperty(prefix+"ly_gt_strength")!=null)                this.ly_gt_strength=Double.parseDouble(properties.getProperty(prefix+"ly_gt_strength"));
		if (properties.getProperty(prefix+"ly_gt_use_wnd")!=null)                 this.ly_gt_use_wnd=Boolean.parseBoolean(properties.getProperty(prefix+"ly_gt_use_wnd"));
		if (properties.getProperty(prefix+"ly_gt_rms")!=null)                     this.ly_gt_rms=Double.parseDouble(properties.getProperty(prefix+"ly_gt_rms"));

		if (properties.getProperty(prefix+"lylw_inf_en")!=null)                     this.lylw_inf_en=Boolean.parseBoolean(properties.getProperty(prefix+"lylw_inf_en"));
		if (properties.getProperty(prefix+"lylw_aztilt_en")!=null)                  this.lylw_aztilt_en=Boolean.parseBoolean(properties.getProperty(prefix+"lylw_aztilt_en"));
		if (properties.getProperty(prefix+"lylw_diff_roll_en")!=null)               this.lylw_diff_roll_en=Boolean.parseBoolean(properties.getProperty(prefix+"lylw_diff_roll_en"));
		if (properties.getProperty(prefix+"lylw_focalLength")!=null)                this.lylw_focalLength=Boolean.parseBoolean(properties.getProperty(prefix+"lylw_focalLength"));
		if (properties.getProperty(prefix+"lylw_com_roll")!=null)                   this.lylw_com_roll=Boolean.parseBoolean(properties.getProperty(prefix+"lylw_com_roll"));
		if (properties.getProperty(prefix+"lylw_par_sel")!=null)                    this.lylw_par_sel=Integer.parseInt(properties.getProperty(prefix+"lylw_par_sel"));

		if (properties.getProperty(prefix+"ly_marg_fract")!=null)                 this.ly_marg_fract=Double.parseDouble(properties.getProperty(prefix+"ly_marg_fract"));
		if (properties.getProperty(prefix+"ly_on_scan")!=null)                    this.ly_on_scan=Boolean.parseBoolean(properties.getProperty(prefix+"ly_on_scan"));
		if (properties.getProperty(prefix+"ly_inf_en")!=null)                     this.ly_inf_en=Boolean.parseBoolean(properties.getProperty(prefix+"ly_inf_en"));
		if (properties.getProperty(prefix+"ly_min_forced")!=null)                 this.ly_min_forced=Integer.parseInt(properties.getProperty(prefix+"ly_min_forced"));
		if (properties.getProperty(prefix+"ly_aztilt_en")!=null)                  this.ly_aztilt_en=Boolean.parseBoolean(properties.getProperty(prefix+"ly_aztilt_en"));
		if (properties.getProperty(prefix+"ly_diff_roll_en")!=null)               this.ly_diff_roll_en=Boolean.parseBoolean(properties.getProperty(prefix+"ly_diff_roll_en"));
		if (properties.getProperty(prefix+"ly_focalLength")!=null)                this.ly_focalLength=Boolean.parseBoolean(properties.getProperty(prefix+"ly_focalLength"));
		if (properties.getProperty(prefix+"ly_com_roll")!=null)                   this.ly_com_roll=Boolean.parseBoolean(properties.getProperty(prefix+"ly_com_roll"));
		if (properties.getProperty(prefix+"ly_ers_rot")!=null)                    this.ly_ers_rot=Boolean.parseBoolean(properties.getProperty(prefix+"ly_ers_rot"));
		if (properties.getProperty(prefix+"ly_ers_forw")!=null)                   this.ly_ers_forw=Boolean.parseBoolean(properties.getProperty(prefix+"ly_ers_forw"));
		if (properties.getProperty(prefix+"ly_ers_side")!=null)                   this.ly_ers_side=Boolean.parseBoolean(properties.getProperty(prefix+"ly_ers_side"));
		if (properties.getProperty(prefix+"ly_ers_vert")!=null)                   this.ly_ers_vert=Boolean.parseBoolean(properties.getProperty(prefix+"ly_ers_vert"));



		if (properties.getProperty(prefix+"ly_par_sel")!=null)                    this.ly_par_sel=Integer.parseInt(properties.getProperty(prefix+"ly_par_sel"));
		if (properties.getProperty(prefix+"ly_debug_level")!=null)                this.ly_debug_level=Integer.parseInt(properties.getProperty(prefix+"ly_debug_level"));

		if (properties.getProperty(prefix+"ly_right_left")!=null)                 this.ly_right_left=Boolean.parseBoolean(properties.getProperty(prefix+"ly_right_left"));

		if (properties.getProperty(prefix+"ly_per_quad")!=null)                   this.ly_per_quad=Integer.parseInt(properties.getProperty(prefix+"ly_per_quad"));
		if (properties.getProperty(prefix+"ly_per_quad_r")!=null)                 this.ly_per_quad_r=Double.parseDouble(properties.getProperty(prefix+"ly_per_quad_r"));
		if (properties.getProperty(prefix+"ly_inf")!=null)                        this.ly_inf=Integer.parseInt(properties.getProperty(prefix+"ly_inf"));
		if (properties.getProperty(prefix+"ly_inf_r")!=null)                      this.ly_inf_r=Double.parseDouble(properties.getProperty(prefix+"ly_inf_r"));
		if (properties.getProperty(prefix+"ly_inf_scale")!=null)                  this.ly_inf_scale=Integer.parseInt(properties.getProperty(prefix+"ly_inf_scale"));
		if (properties.getProperty(prefix+"ly_inf_scale_r")!=null)                this.ly_inf_scale_r=Double.parseDouble(properties.getProperty(prefix+"ly_inf_scale_r"));

		if (properties.getProperty(prefix+"ly_inf_frac")!=null)                   this.ly_inf_frac=Double.parseDouble(properties.getProperty(prefix+"ly_inf_frac"));
		if (properties.getProperty(prefix+"ly_inf_max_disparity")!=null)          this.ly_inf_max_disparity=Double.parseDouble(properties.getProperty(prefix+"ly_inf_max_disparity"));

		if (properties.getProperty(prefix+"ly_inf_disp")!=null)                   this.ly_inf_disp=Boolean.parseBoolean(properties.getProperty(prefix+"ly_inf_disp"));

		if (properties.getProperty(prefix+"ly_inf_force")!=null)                  this.ly_inf_force=Boolean.parseBoolean(properties.getProperty(prefix+"ly_inf_force"));


		if (properties.getProperty(prefix+"ly_poly")!=null)                       this.ly_poly=Boolean.parseBoolean(properties.getProperty(prefix+"ly_poly"));

		if (properties.getProperty(prefix+"ly_smpl_side")!=null)                  this.ly_smpl_side=Integer.parseInt(properties.getProperty(prefix+"ly_smpl_side"));
		if (properties.getProperty(prefix+"ly_smpl_num")!=null)                   this.ly_smpl_num=Integer.parseInt(properties.getProperty(prefix+"ly_smpl_num"));
		if (properties.getProperty(prefix+"ly_smpl_rms")!=null)                   this.ly_smpl_rms=Double.parseDouble(properties.getProperty(prefix+"ly_smpl_rms"));
		if (properties.getProperty(prefix+"ly_disp_var")!=null)                   this.ly_disp_var=Double.parseDouble(properties.getProperty(prefix+"ly_disp_var"));
		if (properties.getProperty(prefix+"ly_disp_rvar")!=null)                  this.ly_disp_rvar=Double.parseDouble(properties.getProperty(prefix+"ly_disp_rvar"));
		if (properties.getProperty(prefix+"ly_disp_var_gt")!=null)                this.ly_disp_var_gt=Double.parseDouble(properties.getProperty(prefix+"ly_disp_var_gt"));
		if (properties.getProperty(prefix+"ly_disp_rvar_gt")!=null)               this.ly_disp_rvar_gt=Double.parseDouble(properties.getProperty(prefix+"ly_disp_rvar_gt"));
		if (properties.getProperty(prefix+"ly_norm_disp")!=null)                  this.ly_norm_disp=Double.parseDouble(properties.getProperty(prefix+"ly_norm_disp"));
		if (properties.getProperty(prefix+"lym_overexp")!=null)                   this.lym_overexp=Double.parseDouble(properties.getProperty(prefix+"lym_overexp"));
		if (properties.getProperty(prefix+"lym_update_disp")!=null)               this.lym_update_disp=Boolean.parseBoolean(properties.getProperty(prefix+"lym_update_disp"));
		if (properties.getProperty(prefix+"lym_iter")!=null)                      this.lym_iter=Integer.parseInt(properties.getProperty(prefix+"lym_iter"));
		if (properties.getProperty(prefix+"lym_change")!=null)                    this.lym_change=Double.parseDouble(properties.getProperty(prefix+"lym_change"));
		if (properties.getProperty(prefix+"lym_change_aux")!=null)                this.lym_change_aux=Double.parseDouble(properties.getProperty(prefix+"lym_change_aux"));

		if (properties.getProperty(prefix+"lym_poly_change")!=null)               this.lym_poly_change=Double.parseDouble(properties.getProperty(prefix+"lym_poly_change"));

		if (properties.getProperty(prefix+"lyf_filter")!=null)                    this.lyf_filter=Boolean.parseBoolean(properties.getProperty(prefix+"lyf_filter"));
		if (properties.getProperty(prefix+"lyf_smpl_side")!=null)                 this.lyf_smpl_side=Integer.parseInt(properties.getProperty(prefix+"lyf_smpl_side"));
		if (properties.getProperty(prefix+"lyf_rms_max")!=null)                   this.lyf_rms_max=Double.parseDouble(properties.getProperty(prefix+"lyf_rms_max"));
		if (properties.getProperty(prefix+"lyf_frac_keep")!=null)                 this.lyf_frac_keep=Double.parseDouble(properties.getProperty(prefix+"lyf_frac_keep"));
		if (properties.getProperty(prefix+"lyf_min_samples")!=null)               this.lyf_min_samples=Integer.parseInt(properties.getProperty(prefix+"lyf_min_samples"));
		if (properties.getProperty(prefix+"lyf_norm_center")!=null)               this.lyf_norm_center=Boolean.parseBoolean(properties.getProperty(prefix+"lyf_norm_center"));
		if (properties.getProperty(prefix+"ly_corr_scale")!=null)                 this.ly_corr_scale=Double.parseDouble(properties.getProperty(prefix+"ly_corr_scale"));

		if (properties.getProperty(prefix+"lyr_filter_ds")!=null)                 this.lyr_filter_ds=Boolean.parseBoolean(properties.getProperty(prefix+"lyr_filter_ds"));
		if (properties.getProperty(prefix+"lyr_filter_lyf")!=null)                this.lyr_filter_lyf=Boolean.parseBoolean(properties.getProperty(prefix+"lyr_filter_lyf"));

		if (properties.getProperty(prefix+"corr_magic_scale")!=null)              this.corr_magic_scale=Double.parseDouble(properties.getProperty(prefix+"corr_magic_scale"));
		if (properties.getProperty(prefix+"corr_select")!=null)                   this.corr_select=Integer.parseInt(properties.getProperty(prefix+"corr_select"));

		if (properties.getProperty(prefix+"show_textures")!=null)                 this.show_textures=Boolean.parseBoolean(properties.getProperty(prefix+"show_textures"));
		if (properties.getProperty(prefix+"debug_filters")!=null)                 this.debug_filters=Boolean.parseBoolean(properties.getProperty(prefix+"debug_filters"));

		if (properties.getProperty(prefix+"min_smth")!=null)                      this.min_smth=Double.parseDouble(properties.getProperty(prefix+"min_smth"));
		if (properties.getProperty(prefix+"sure_smth")!=null)                     this.sure_smth=Double.parseDouble(properties.getProperty(prefix+"sure_smth"));
		if (properties.getProperty(prefix+"bgnd_range")!=null)                    this.bgnd_range=Double.parseDouble(properties.getProperty(prefix+"bgnd_range"));
		if (properties.getProperty(prefix+"other_range")!=null)                   this.other_range=Double.parseDouble(properties.getProperty(prefix+"other_range"));

		if (properties.getProperty(prefix+"ex_strength")!=null)                   this.ex_strength=Double.parseDouble(properties.getProperty(prefix+"ex_strength"));
		if (properties.getProperty(prefix+"ex_nstrength")!=null)                  this.ex_nstrength=Double.parseDouble(properties.getProperty(prefix+"ex_nstrength"));

		if (properties.getProperty(prefix+"ex_over_bgnd")!=null)                  this.ex_over_bgnd=Boolean.parseBoolean(properties.getProperty(prefix+"ex_over_bgnd"));
		if (properties.getProperty(prefix+"ex_min_over")!=null)                   this.ex_min_over=Double.parseDouble(properties.getProperty(prefix+"ex_min_over"));

		if (properties.getProperty(prefix+"pt_super_trust")!=null)                this.pt_super_trust=Double.parseDouble(properties.getProperty(prefix+"pt_super_trust"));
		if (properties.getProperty(prefix+"pt_keep_raw_fg")!=null)                this.pt_keep_raw_fg=Boolean.parseBoolean(properties.getProperty(prefix+"pt_keep_raw_fg"));
		if (properties.getProperty(prefix+"pt_scale_pre")!=null)                  this.pt_scale_pre=Double.parseDouble(properties.getProperty(prefix+"pt_scale_pre"));
		if (properties.getProperty(prefix+"pt_scale_post")!=null)                 this.pt_scale_post=Double.parseDouble(properties.getProperty(prefix+"pt_scale_post"));

		if (properties.getProperty(prefix+"bgnd_sure")!=null)                     this.bgnd_sure=Double.parseDouble(properties.getProperty(prefix+"bgnd_sure"));
		if (properties.getProperty(prefix+"bgnd_maybe")!=null)                    this.bgnd_maybe=Double.parseDouble(properties.getProperty(prefix+"bgnd_maybe"));
		if (properties.getProperty(prefix+"min_clstr_seed")!=null)                this.min_clstr_seed=Integer.parseInt(properties.getProperty(prefix+"min_clstr_seed"));
		if (properties.getProperty(prefix+"min_clstr_lone")!=null)                this.min_clstr_lone=Integer.parseInt(properties.getProperty(prefix+"min_clstr_lone"));
		if (properties.getProperty(prefix+"min_clstr_weight")!=null)              this.min_clstr_weight=Double.parseDouble(properties.getProperty(prefix+"min_clstr_weight"));
		if (properties.getProperty(prefix+"min_clstr_max")!=null)                 this.min_clstr_max=Double.parseDouble(properties.getProperty(prefix+"min_clstr_max"));


		if (properties.getProperty(prefix+"fill_gaps")!=null)                     this.fill_gaps=Integer.parseInt(properties.getProperty(prefix+"fill_gaps"));
		if (properties.getProperty(prefix+"fill_final")!=null)                    this.fill_final=Integer.parseInt(properties.getProperty(prefix+"fill_final"));
		if (properties.getProperty(prefix+"min_clstr_block")!=null)               this.min_clstr_block=Integer.parseInt(properties.getProperty(prefix+"min_clstr_block"));
		if (properties.getProperty(prefix+"bgnd_grow")!=null)                     this.bgnd_grow=Integer.parseInt(properties.getProperty(prefix+"bgnd_grow"));

		if (properties.getProperty(prefix+"ortho_old")!=null)                     this.ortho_old=Boolean.parseBoolean(properties.getProperty(prefix+"ortho_old"));
		if (properties.getProperty(prefix+"ortho_min_hor")!=null)                 this.ortho_min_hor=Double.parseDouble(properties.getProperty(prefix+"ortho_min_hor"));
		if (properties.getProperty(prefix+"ortho_min_vert")!=null)                this.ortho_min_vert=Double.parseDouble(properties.getProperty(prefix+"ortho_min_vert"));
		if (properties.getProperty(prefix+"ortho_asym")!=null)                    this.ortho_asym=Double.parseDouble(properties.getProperty(prefix+"ortho_asym"));
		if (properties.getProperty(prefix+"ortho_over4")!=null)                   this.ortho_over4=Double.parseDouble(properties.getProperty(prefix+"ortho_over4"));
		if (properties.getProperty(prefix+"ortho_sustain")!=null)                 this.ortho_sustain=Double.parseDouble(properties.getProperty(prefix+"ortho_sustain"));
		if (properties.getProperty(prefix+"ortho_run")!=null)                     this.ortho_run=Integer.parseInt(properties.getProperty(prefix+"ortho_run"));
		if (properties.getProperty(prefix+"ortho_minmax")!=null)                  this.ortho_minmax=Double.parseDouble(properties.getProperty(prefix+"ortho_minmax"));
		if (properties.getProperty(prefix+"ortho_bridge")!=null)                  this.ortho_bridge=Integer.parseInt(properties.getProperty(prefix+"ortho_bridge"));
		if (properties.getProperty(prefix+"ortho_rms")!=null)                     this.ortho_rms=Double.parseDouble(properties.getProperty(prefix+"ortho_rms"));
		if (properties.getProperty(prefix+"ortho_half_length")!=null)             this.ortho_half_length=Integer.parseInt(properties.getProperty(prefix+"ortho_half_length"));
		if (properties.getProperty(prefix+"ortho_mix")!=null)                     this.ortho_mix=Double.parseDouble(properties.getProperty(prefix+"ortho_mix"));

		if (properties.getProperty(prefix+"or_hor")!=null)                        this.or_hor=Boolean.parseBoolean(properties.getProperty(prefix+"or_hor"));
		if (properties.getProperty(prefix+"or_vert")!=null)                       this.or_vert=Boolean.parseBoolean(properties.getProperty(prefix+"or_vert"));
		if (properties.getProperty(prefix+"or_sigma")!=null)                      this.or_sigma=Double.parseDouble(properties.getProperty(prefix+"or_sigma"));
		if (properties.getProperty(prefix+"or_sharp")!=null)                      this.or_sharp=Double.parseDouble(properties.getProperty(prefix+"or_sharp"));
		if (properties.getProperty(prefix+"or_scale")!=null)                      this.or_scale=Double.parseDouble(properties.getProperty(prefix+"or_scale"));
		if (properties.getProperty(prefix+"or_offset")!=null)                     this.or_offset=Double.parseDouble(properties.getProperty(prefix+"or_offset"));
		if (properties.getProperty(prefix+"or_asym")!=null)                       this.or_asym=Double.parseDouble(properties.getProperty(prefix+"or_asym"));
		if (properties.getProperty(prefix+"or_threshold")!=null)                  this.or_threshold=Double.parseDouble(properties.getProperty(prefix+"or_threshold"));
		if (properties.getProperty(prefix+"or_absHor")!=null)                     this.or_absHor=Double.parseDouble(properties.getProperty(prefix+"or_absHor"));
		if (properties.getProperty(prefix+"or_absVert")!=null)                    this.or_absVert=Double.parseDouble(properties.getProperty(prefix+"or_absVert"));
		if (properties.getProperty(prefix+"or_maxDisp")!=null)                    this.or_maxDisp=Double.parseDouble(properties.getProperty(prefix+"or_maxDisp"));

		if (properties.getProperty(prefix+"poles_fix")!=null)                     this.poles_fix=Boolean.parseBoolean(properties.getProperty(prefix+"poles_fix"));
		if (properties.getProperty(prefix+"poles_len")!=null)                     this.poles_len=Integer.parseInt(properties.getProperty(prefix+"poles_len"));
		if (properties.getProperty(prefix+"poles_ratio")!=null)                   this.poles_ratio=Double.parseDouble(properties.getProperty(prefix+"poles_ratio"));
		if (properties.getProperty(prefix+"poles_min_strength")!=null)            this.poles_min_strength=Double.parseDouble(properties.getProperty(prefix+"poles_min_strength"));
		if (properties.getProperty(prefix+"poles_force_disp")!=null)              this.poles_force_disp=Boolean.parseBoolean(properties.getProperty(prefix+"poles_force_disp"));

		if (properties.getProperty(prefix+"max_clusters")!=null)                  this.max_clusters=Integer.parseInt(properties.getProperty(prefix+"max_clusters"));
		if (properties.getProperty(prefix+"remove_scans")!=null)                  this.remove_scans=Boolean.parseBoolean(properties.getProperty(prefix+"remove_scans"));
		if (properties.getProperty(prefix+"output_x3d")!=null)                    this.output_x3d=Boolean.parseBoolean(properties.getProperty(prefix+"output_x3d"));
		if (properties.getProperty(prefix+"output_obj")!=null)                    this.output_obj=Boolean.parseBoolean(properties.getProperty(prefix+"output_obj"));
		if (properties.getProperty(prefix+"correct_distortions")!=null)           this.correct_distortions=Boolean.parseBoolean(properties.getProperty(prefix+"correct_distortions"));
		if (properties.getProperty(prefix+"show_triangles")!=null)                this.show_triangles=Boolean.parseBoolean(properties.getProperty(prefix+"show_triangles"));
		if (properties.getProperty(prefix+"avg_cluster_disp")!=null)              this.avg_cluster_disp=Boolean.parseBoolean(properties.getProperty(prefix+"avg_cluster_disp"));
		if (properties.getProperty(prefix+"maxDispTriangle")!=null)               this.maxDispTriangle=Double.parseDouble(properties.getProperty(prefix+"maxDispTriangle"));
		if (properties.getProperty(prefix+"infinityDistance")!=null)              this.infinityDistance=Double.parseDouble(properties.getProperty(prefix+"infinityDistance"));
		if (properties.getProperty(prefix+"min_bgnd_tiles")!=null)                this.min_bgnd_tiles=Integer.parseInt(properties.getProperty(prefix+"min_bgnd_tiles"));
		if (properties.getProperty(prefix+"shUseFlaps")!=null)                    this.shUseFlaps=Boolean.parseBoolean(properties.getProperty(prefix+"shUseFlaps"));
		if (properties.getProperty(prefix+"shAggrFade")!=null)                    this.shAggrFade=Boolean.parseBoolean(properties.getProperty(prefix+"shAggrFade"));
		if (properties.getProperty(prefix+"shMinArea")!=null)                     this.shMinArea=Integer.parseInt(properties.getProperty(prefix+"shMinArea"));
		if (properties.getProperty(prefix+"shMinStrength")!=null)                 this.shMinStrength=Double.parseDouble(properties.getProperty(prefix+"shMinStrength"));
		if (properties.getProperty(prefix+"tiRigidVertical")!=null)               this.tiRigidVertical=Double.parseDouble(properties.getProperty(prefix+"tiRigidVertical"));
		if (properties.getProperty(prefix+"tiRigidHorizontal")!=null)             this.tiRigidHorizontal=Double.parseDouble(properties.getProperty(prefix+"tiRigidHorizontal"));
		if (properties.getProperty(prefix+"tiRigidDiagonal")!=null)               this.tiRigidDiagonal=Double.parseDouble(properties.getProperty(prefix+"tiRigidDiagonal"));
		if (properties.getProperty(prefix+"tiStrengthOffset")!=null)              this.tiStrengthOffset=Double.parseDouble(properties.getProperty(prefix+"tiStrengthOffset"));
		if (properties.getProperty(prefix+"tiDispScale")!=null)                   this.tiDispScale=Double.parseDouble(properties.getProperty(prefix+"tiDispScale"));
		if (properties.getProperty(prefix+"tiDispPow")!=null)                     this.tiDispPow=Double.parseDouble(properties.getProperty(prefix+"tiDispPow"));
		if (properties.getProperty(prefix+"tiDispPull")!=null)                    this.tiDispPull=Double.parseDouble(properties.getProperty(prefix+"tiDispPull"));
		if (properties.getProperty(prefix+"tiDispPullPreFinal")!=null)            this.tiDispPullPreFinal=Double.parseDouble(properties.getProperty(prefix+"tiDispPullPreFinal"));
		if (properties.getProperty(prefix+"tiDispPullFinal")!=null)               this.tiDispPullFinal=Double.parseDouble(properties.getProperty(prefix+"tiDispPullFinal"));
		if (properties.getProperty(prefix+"tiBreakNorm")!=null)                   this.tiBreakNorm=Double.parseDouble(properties.getProperty(prefix+"tiBreakNorm"));
		if (properties.getProperty(prefix+"tiBreak3")!=null)                      this.tiBreak3=Double.parseDouble(properties.getProperty(prefix+"tiBreak3"));
		if (properties.getProperty(prefix+"tiBreak31")!=null)                     this.tiBreak31=Double.parseDouble(properties.getProperty(prefix+"tiBreak31"));
		if (properties.getProperty(prefix+"tiBreak21")!=null)                     this.tiBreak21=Double.parseDouble(properties.getProperty(prefix+"tiBreak21"));
		if (properties.getProperty(prefix+"tiBreakFar")!=null)                    this.tiBreakFar=Double.parseDouble(properties.getProperty(prefix+"tiBreakFar"));
		if (properties.getProperty(prefix+"tiBreakNear")!=null)                   this.tiBreakNear=Double.parseDouble(properties.getProperty(prefix+"tiBreakNear"));
		if (properties.getProperty(prefix+"tiBreakMode")!=null)                   this.tiBreakMode=Integer.parseInt(properties.getProperty(prefix+"tiBreakMode"));
		if (properties.getProperty(prefix+"tiBreakSame")!=null)                   this.tiBreakSame=Double.parseDouble(properties.getProperty(prefix+"tiBreakSame"));
		if (properties.getProperty(prefix+"tiBreakTurn")!=null)                   this.tiBreakTurn=Double.parseDouble(properties.getProperty(prefix+"tiBreakTurn"));

		if (properties.getProperty(prefix+"tiHealPreLast")!=null)                 this.tiHealPreLast=Double.parseDouble(properties.getProperty(prefix+"tiHealPreLast"));
		if (properties.getProperty(prefix+"tiHealLast")!=null)                    this.tiHealLast=Double.parseDouble(properties.getProperty(prefix+"tiHealLast"));
		if (properties.getProperty(prefix+"tiHealSame")!=null)                    this.tiHealSame=Integer.parseInt(properties.getProperty(prefix+"tiHealSame"));

		if (properties.getProperty(prefix+"tiIterations")!=null)                  this.tiIterations=Integer.parseInt(properties.getProperty(prefix+"tiIterations"));
		if (properties.getProperty(prefix+"tiPrecision")!=null)                   this.tiPrecision=Integer.parseInt(properties.getProperty(prefix+"tiPrecision"));
		if (properties.getProperty(prefix+"tiNumCycles")!=null)                   this.tiNumCycles=Integer.parseInt(properties.getProperty(prefix+"tiNumCycles"));

		if (properties.getProperty(prefix+"stUseRefine")!=null)                   this.stUseRefine=Boolean.parseBoolean(properties.getProperty(prefix+"stUseRefine"));
		if (properties.getProperty(prefix+"stUsePass2")!=null)                    this.stUsePass2=Boolean.parseBoolean(properties.getProperty(prefix+"stUsePass2"));
		if (properties.getProperty(prefix+"stUseRender")!=null)                   this.stUseRender=Boolean.parseBoolean(properties.getProperty(prefix+"stUseRender"));

		if (properties.getProperty(prefix+"stShow")!=null)                        this.stShow=Boolean.parseBoolean(properties.getProperty(prefix+"stShow"));
		if (properties.getProperty(prefix+"stSize")!=null)                        this.stSize=Integer.parseInt(properties.getProperty(prefix+"stSize"));
		if (properties.getProperty(prefix+"stStepFar")!=null)                     this.stStepFar=Double.parseDouble(properties.getProperty(prefix+"stStepFar"));
		if (properties.getProperty(prefix+"stStepNear")!=null)                    this.stStepNear=Double.parseDouble(properties.getProperty(prefix+"stStepNear"));
		if (properties.getProperty(prefix+"stStepThreshold")!=null)               this.stStepThreshold=Double.parseDouble(properties.getProperty(prefix+"stStepThreshold"));
		if (properties.getProperty(prefix+"stMinDisparity")!=null)                this.stMinDisparity=Double.parseDouble(properties.getProperty(prefix+"stMinDisparity"));
		//  			if (properties.getProperty(prefix+"stMaxDisparity")!=null)                this.stMaxDisparity=Double.parseDouble(properties.getProperty(prefix+"stMaxDisparity"));

		if (properties.getProperty(prefix+"stSigma")!=null)                       this.stSigma=Double.parseDouble(properties.getProperty(prefix+"stSigma"));
		if (properties.getProperty(prefix+"stMinBgDisparity")!=null)              this.stMinBgDisparity=Double.parseDouble(properties.getProperty(prefix+"stMinBgDisparity"));
		if (properties.getProperty(prefix+"stMinBgFract")!=null)                  this.stMinBgFract=Double.parseDouble(properties.getProperty(prefix+"stMinBgFract"));
		if (properties.getProperty(prefix+"stUseDisp")!=null)                     this.stUseDisp=Double.parseDouble(properties.getProperty(prefix+"stUseDisp"));
		if (properties.getProperty(prefix+"stStrengthScale")!=null)               this.stStrengthScale=Double.parseDouble(properties.getProperty(prefix+"stStrengthScale"));

		if (properties.getProperty(prefix+"stSmplMode")!=null)                    this.stSmplMode=Boolean.parseBoolean(properties.getProperty(prefix+"stSmplMode"));


		// next kept for compatibility with old configs
		if (properties.getProperty(prefix+"stFloor")!=null)                       this.mlfp.strength_floor=Double.parseDouble(properties.getProperty(prefix+"stFloor"));
		if (properties.getProperty(prefix+"stPow")!=null)                         this.mlfp.strength_pow=Double.parseDouble(properties.getProperty(prefix+"stPow"));
		if (properties.getProperty(prefix+"stSmplSide")!=null)                    this.mlfp.smplSide=Integer.parseInt(properties.getProperty(prefix+"stSmplSide"));
		if (properties.getProperty(prefix+"stSmplNum")!=null)                     this.mlfp.smplNum=Integer.parseInt(properties.getProperty(prefix+"stSmplNum"));
		if (properties.getProperty(prefix+"stSmplRms")!=null)                     this.mlfp.smplRms=Double.parseDouble(properties.getProperty(prefix+"stSmplRms"));
		if (properties.getProperty(prefix+"stSmplWnd")!=null)                     this.mlfp.smplWnd=Boolean.parseBoolean(properties.getProperty(prefix+"stSmplWnd"));
		if (properties.getProperty(prefix+"fs_max_abs_tilt")!=null)               this.mlfp.max_abs_tilt=Double.parseDouble(properties.getProperty(prefix+"fs_max_abs_tilt"));
		if (properties.getProperty(prefix+"fs_max_rel_tilt")!=null)               this.mlfp.max_rel_tilt=Double.parseDouble(properties.getProperty(prefix+"fs_max_rel_tilt"));
		if (properties.getProperty(prefix+"fs_damp_tilt")!=null)                  this.mlfp.damp_tilt=Double.parseDouble(properties.getProperty(prefix+"fs_damp_tilt"));
		if (properties.getProperty(prefix+"fs_min_tilt_disp")!=null)              this.mlfp.min_tilt_disp=Double.parseDouble(properties.getProperty(prefix+"fs_min_tilt_disp"));
		if (properties.getProperty(prefix+"fs_transition")!=null)                 this.mlfp.transition=Double.parseDouble(properties.getProperty(prefix+"fs_transition"));
		if (properties.getProperty(prefix+"fs_far_mode")!=null)                   this.mlfp.far_mode=Integer.parseInt(properties.getProperty(prefix+"fs_far_mode"));
		if (properties.getProperty(prefix+"fs_far_power")!=null)                  this.mlfp.far_power=Double.parseDouble(properties.getProperty(prefix+"fs_far_power"));
		// end of next kept for compatibility with old configs


		if (properties.getProperty(prefix+"stGrowSel")!=null)                     this.stGrowSel=Integer.parseInt(properties.getProperty(prefix+"stGrowSel"));
		if (properties.getProperty(prefix+"stMeasSel")!=null)                     this.stMeasSel=Integer.parseInt(properties.getProperty(prefix+"stMeasSel"));
		if (properties.getProperty(prefix+"stSmallDiff")!=null)                   this.stSmallDiff=Double.parseDouble(properties.getProperty(prefix+"stSmallDiff"));
		if (properties.getProperty(prefix+"stHighMix")!=null)                     this.stHighMix=Double.parseDouble(properties.getProperty(prefix+"stHighMix"));

		// old wrong spelling
		if (properties.getProperty(prefix+"outlayerStrength")!=null)              this.outlierStrength=Double.parseDouble(properties.getProperty(prefix+"outlayerStrength"));
		if (properties.getProperty(prefix+"outlayerDiff")!=null)                  this.outlierDiff=Double.parseDouble(properties.getProperty(prefix+"outlayerDiff"));
		if (properties.getProperty(prefix+"outlayerDiffPos")!=null)               this.outlierDiffPos=Double.parseDouble(properties.getProperty(prefix+"outlayerDiffPos"));
		if (properties.getProperty(prefix+"outlayerDiffNeg")!=null)               this.outlierDiffNeg=Double.parseDouble(properties.getProperty(prefix+"outlayerDiffNeg"));


		if (properties.getProperty(prefix+"outlierStrength")!=null)               this.outlierStrength=Double.parseDouble(properties.getProperty(prefix+"outlierStrength"));
		if (properties.getProperty(prefix+"outlierDiff")!=null)                   this.outlierDiff=Double.parseDouble(properties.getProperty(prefix+"outlierDiff"));
		if (properties.getProperty(prefix+"outlierDiffPos")!=null)                this.outlierDiffPos=Double.parseDouble(properties.getProperty(prefix+"outlierDiffPos"));
		if (properties.getProperty(prefix+"outlierDiffNeg")!=null)                this.outlierDiffNeg=Double.parseDouble(properties.getProperty(prefix+"outlierDiffNeg"));

		if (properties.getProperty(prefix+"combine_refine")!=null)                this.combine_refine=Boolean.parseBoolean(properties.getProperty(prefix+"combine_refine"));

		if (properties.getProperty(prefix+"combine_min_strength")!=null)          this.combine_min_strength=Double.parseDouble(properties.getProperty(prefix+"combine_min_strength"));
		if (properties.getProperty(prefix+"combine_min_hor")!=null)               this.combine_min_hor=Double.parseDouble(properties.getProperty(prefix+"combine_min_hor"));
		if (properties.getProperty(prefix+"combine_min_vert")!=null)              this.combine_min_vert=Double.parseDouble(properties.getProperty(prefix+"combine_min_vert"));
		//  			if (properties.getProperty(prefix+"unique_tolerance")!=null)              this.unique_tolerance=Double.parseDouble(properties.getProperty(prefix+"unique_tolerance"));
		if (properties.getProperty(prefix+"grow_sweep")!=null)                    this.grow_sweep=Integer.parseInt(properties.getProperty(prefix+"grow_sweep"));
		if (properties.getProperty(prefix+"grow_disp_max")!=null)                 this.grow_disp_max=Double.parseDouble(properties.getProperty(prefix+"grow_disp_max"));
		if (properties.getProperty(prefix+"grow_disp_trust")!=null)               this.grow_disp_trust=Double.parseDouble(properties.getProperty(prefix+"grow_disp_trust"));
		if (properties.getProperty(prefix+"grow_disp_step")!=null)                this.grow_disp_step=Double.parseDouble(properties.getProperty(prefix+"grow_disp_step"));
		if (properties.getProperty(prefix+"grow_min_diff")!=null)                 this.grow_min_diff=Double.parseDouble(properties.getProperty(prefix+"grow_min_diff"));
		if (properties.getProperty(prefix+"grow_retry_far")!=null)                this.grow_retry_far=Boolean.parseBoolean(properties.getProperty(prefix+"grow_retry_far"));
		if (properties.getProperty(prefix+"grow_pedantic")!=null)                 this.grow_pedantic=Boolean.parseBoolean(properties.getProperty(prefix+"grow_pedantic"));
		if (properties.getProperty(prefix+"grow_retry_inf")!=null)                this.grow_retry_inf=Boolean.parseBoolean(properties.getProperty(prefix+"grow_retry_inf"));

		if (properties.getProperty(prefix+"gr_new_expand")!=null)                 this.gr_new_expand=Boolean.parseBoolean(properties.getProperty(prefix+"gr_new_expand"));
		//  			if (properties.getProperty(prefix+"gr_max_expand")!=null)                 this.gr_max_expand=Integer.parseInt(properties.getProperty(prefix+"gr_max_expand"));
		if (properties.getProperty(prefix+"fds_str_floor")!=null)                 this.fds_str_floor=Double.parseDouble(properties.getProperty(prefix+"fds_str_floor"));
		if (properties.getProperty(prefix+"gr_ovrbg_cmb")!=null)                  this.gr_ovrbg_cmb=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_cmb"));
		if (properties.getProperty(prefix+"gr_ovrbg_cmb_hor")!=null)              this.gr_ovrbg_cmb_hor=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_cmb_hor"));
		if (properties.getProperty(prefix+"gr_ovrbg_cmb_vert")!=null)             this.gr_ovrbg_cmb_vert=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_cmb_vert"));
		if (properties.getProperty(prefix+"gr_ovrbg_filtered")!=null)             this.gr_ovrbg_filtered=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_filtered"));
		if (properties.getProperty(prefix+"fds_str_pow")!=null)                   this.fds_str_pow=Double.parseDouble(properties.getProperty(prefix+"fds_str_pow"));
		if (properties.getProperty(prefix+"fds_smpl_side")!=null)                 this.fds_smpl_side=Integer.parseInt(properties.getProperty(prefix+"fds_smpl_side"));
		if (properties.getProperty(prefix+"fds_smpl_num")!=null)                  this.fds_smpl_num=Integer.parseInt(properties.getProperty(prefix+"fds_smpl_num"));
		if (properties.getProperty(prefix+"fds_smpl_rms")!=null)                  this.fds_smpl_rms=Double.parseDouble(properties.getProperty(prefix+"fds_smpl_rms"));
		if (properties.getProperty(prefix+"fds_smpl_rel_rms")!=null)              this.fds_smpl_rel_rms=Double.parseDouble(properties.getProperty(prefix+"fds_smpl_rel_rms"));
		if (properties.getProperty(prefix+"fds_smpl_wnd")!=null)                  this.fds_smpl_wnd=Boolean.parseBoolean(properties.getProperty(prefix+"fds_smpl_wnd"));
		if (properties.getProperty(prefix+"fds_abs_tilt")!=null)                  this.fds_abs_tilt=Double.parseDouble(properties.getProperty(prefix+"fds_abs_tilt"));
		if (properties.getProperty(prefix+"fds_rel_tilt")!=null)                  this.fds_rel_tilt=Double.parseDouble(properties.getProperty(prefix+"fds_rel_tilt"));

		if (properties.getProperty(prefix+"per_filter")!=null)                    this.per_filter=Boolean.parseBoolean(properties.getProperty(prefix+"per_filter"));
		if (properties.getProperty(prefix+"per_trustedCorrelation")!=null)        this.per_trustedCorrelation=Double.parseDouble(properties.getProperty(prefix+"per_trustedCorrelation"));
		if (properties.getProperty(prefix+"per_initial_diff")!=null)              this.per_initial_diff=Double.parseDouble(properties.getProperty(prefix+"per_initial_diff"));
		if (properties.getProperty(prefix+"per_strength_floor")!=null)            this.per_strength_floor=Double.parseDouble(properties.getProperty(prefix+"per_strength_floor"));
		if (properties.getProperty(prefix+"per_strength_max_over")!=null)         this.per_strength_max_over=Double.parseDouble(properties.getProperty(prefix+"per_strength_max_over"));
		if (properties.getProperty(prefix+"per_min_period")!=null)                this.per_min_period=Double.parseDouble(properties.getProperty(prefix+"per_min_period"));
		if (properties.getProperty(prefix+"per_min_num_periods")!=null)           this.per_min_num_periods=Integer.parseInt(properties.getProperty(prefix+"per_min_num_periods"));
		if (properties.getProperty(prefix+"per_disp_tolerance")!=null)            this.per_disp_tolerance=Double.parseDouble(properties.getProperty(prefix+"per_disp_tolerance"));
		if (properties.getProperty(prefix+"per_disp_match")!=null)                this.per_disp_match=Double.parseDouble(properties.getProperty(prefix+"per_disp_match"));
		if (properties.getProperty(prefix+"per_strong_match_inc")!=null)          this.per_strong_match_inc=Double.parseDouble(properties.getProperty(prefix+"per_strong_match_inc"));


		if (properties.getProperty(prefix+"mc_disp8_step")!=null)                 this.mc_disp8_step=Double.parseDouble(properties.getProperty(prefix+"mc_disp8_step"));
		if (properties.getProperty(prefix+"mc_disp8_trust")!=null)                this.mc_disp8_trust=Double.parseDouble(properties.getProperty(prefix+"mc_disp8_trust"));
		if (properties.getProperty(prefix+"mc_strength")!=null)                   this.mc_strength=Double.parseDouble(properties.getProperty(prefix+"mc_strength"));
		if (properties.getProperty(prefix+"mc_unique_tol")!=null)                 this.mc_unique_tol=Double.parseDouble(properties.getProperty(prefix+"mc_unique_tol"));
		if (properties.getProperty(prefix+"mc_trust_fin")!=null)                  this.mc_trust_fin=Double.parseDouble(properties.getProperty(prefix+"mc_trust_fin"));
		if (properties.getProperty(prefix+"mc_trust_sigma")!=null)                this.mc_trust_sigma=Double.parseDouble(properties.getProperty(prefix+"mc_trust_sigma"));
		if (properties.getProperty(prefix+"mc_ortho_weight")!=null)               this.mc_ortho_weight=Double.parseDouble(properties.getProperty(prefix+"mc_ortho_weight"));
		if (properties.getProperty(prefix+"mc_diag_weight")!=null)                this.mc_diag_weight=Double.parseDouble(properties.getProperty(prefix+"mc_diag_weight"));
		if (properties.getProperty(prefix+"mc_gap")!=null)                        this.mc_gap=Double.parseDouble(properties.getProperty(prefix+"mc_gap"));

		if (properties.getProperty(prefix+"mc_weight_var")!=null)                 this.mc_weight_var=Double.parseDouble(properties.getProperty(prefix+"mc_weight_var"));
		if (properties.getProperty(prefix+"mc_weight_Y")!=null)                   this.mc_weight_Y=Double.parseDouble(properties.getProperty(prefix+"mc_weight_Y"));
		if (properties.getProperty(prefix+"mc_weight_RBmG")!=null)                this.mc_weight_RBmG=Double.parseDouble(properties.getProperty(prefix+"mc_weight_RBmG"));

		if (properties.getProperty(prefix+"gr_min_new")!=null)                    this.gr_min_new=Integer.parseInt(properties.getProperty(prefix+"gr_min_new"));
		if (properties.getProperty(prefix+"gr_var_new_sngl")!=null)               this.gr_var_new_sngl=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_new_sngl"));
		if (properties.getProperty(prefix+"gr_var_new_fg")!=null)                 this.gr_var_new_fg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_new_fg"));
		if (properties.getProperty(prefix+"gr_var_all_fg")!=null)                 this.gr_var_all_fg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_all_fg"));
		if (properties.getProperty(prefix+"gr_var_new_bg")!=null)                 this.gr_var_new_bg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_new_bg"));
		if (properties.getProperty(prefix+"gr_var_all_bg")!=null)                 this.gr_var_all_bg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_all_bg"));
		if (properties.getProperty(prefix+"gr_var_next")!=null)                   this.gr_var_next=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_next"));
		if (properties.getProperty(prefix+"gr_num_steps")!=null)                  this.gr_num_steps=Integer.parseInt(properties.getProperty(prefix+"gr_num_steps"));
		if (properties.getProperty(prefix+"gr_steps_over")!=null)                 this.gr_steps_over=Integer.parseInt(properties.getProperty(prefix+"gr_steps_over"));
		if (properties.getProperty(prefix+"gr_smpl_size")!=null)                  this.gr_smpl_size=Integer.parseInt(properties.getProperty(prefix+"gr_smpl_size"));
		if (properties.getProperty(prefix+"gr_min_pnts")!=null)                   this.gr_min_pnts=Integer.parseInt(properties.getProperty(prefix+"gr_min_pnts"));
		if (properties.getProperty(prefix+"gr_use_wnd")!=null)                    this.gr_use_wnd=Boolean.parseBoolean(properties.getProperty(prefix+"gr_use_wnd"));
		if (properties.getProperty(prefix+"gr_tilt_damp")!=null)                  this.gr_tilt_damp=Double.parseDouble(properties.getProperty(prefix+"gr_tilt_damp"));
		if (properties.getProperty(prefix+"gr_split_rng")!=null)                  this.gr_split_rng=Double.parseDouble(properties.getProperty(prefix+"gr_split_rng"));
		if (properties.getProperty(prefix+"gr_same_rng")!=null)                   this.gr_same_rng=Double.parseDouble(properties.getProperty(prefix+"gr_same_rng"));
		if (properties.getProperty(prefix+"gr_diff_cont")!=null)                  this.gr_diff_cont=Double.parseDouble(properties.getProperty(prefix+"gr_diff_cont"));
		if (properties.getProperty(prefix+"gr_abs_tilt")!=null)                   this.gr_abs_tilt=Double.parseDouble(properties.getProperty(prefix+"gr_abs_tilt"));
		if (properties.getProperty(prefix+"gr_rel_tilt")!=null)                   this.gr_rel_tilt=Double.parseDouble(properties.getProperty(prefix+"gr_rel_tilt"));
		if (properties.getProperty(prefix+"gr_smooth")!=null)                     this.gr_smooth=Integer.parseInt(properties.getProperty(prefix+"gr_smooth"));
		if (properties.getProperty(prefix+"gr_fin_diff")!=null)                   this.gr_fin_diff=Double.parseDouble(properties.getProperty(prefix+"gr_fin_diff"));
		if (properties.getProperty(prefix+"gr_unique_tol")!=null)                 this.gr_unique_tol=Double.parseDouble(properties.getProperty(prefix+"gr_unique_tol"));
		if (properties.getProperty(prefix+"gr_unique_pretol")!=null)              this.gr_unique_pretol=Double.parseDouble(properties.getProperty(prefix+"gr_unique_pretol"));


		if (properties.getProperty(prefix+"ft_mod_strength")!=null)                    this.ft_mod_strength=Boolean.parseBoolean(properties.getProperty(prefix+"ft_mod_strength"));
		if (properties.getProperty(prefix+"ft_clusterize_by_highest")!=null)           this.ft_clusterize_by_highest=Boolean.parseBoolean(properties.getProperty(prefix+"ft_clusterize_by_highest"));
		if (properties.getProperty(prefix+"ft_clust_sigma")!=null)                     this.ft_clust_sigma=Double.parseDouble(properties.getProperty(prefix+"ft_clust_sigma"));
		if (properties.getProperty(prefix+"ft_disp_arange_vert")!=null)                this.ft_disp_arange_vert=Double.parseDouble(properties.getProperty(prefix+"ft_disp_arange_vert"));
		if (properties.getProperty(prefix+"ft_disp_rrange_vert")!=null)                this.ft_disp_rrange_vert=Double.parseDouble(properties.getProperty(prefix+"ft_disp_rrange_vert"));
		if (properties.getProperty(prefix+"ft_disp_arange_hor")!=null)                 this.ft_disp_arange_hor=Double.parseDouble(properties.getProperty(prefix+"ft_disp_arange_hor"));
		if (properties.getProperty(prefix+"ft_disp_rrange_hor")!=null)                 this.ft_disp_rrange_hor=Double.parseDouble(properties.getProperty(prefix+"ft_disp_rrange_hor"));
		if (properties.getProperty(prefix+"ft_tolerance_above_near")!=null)            this.ft_tolerance_above_near=Double.parseDouble(properties.getProperty(prefix+"ft_tolerance_above_near"));
		if (properties.getProperty(prefix+"ft_tolerance_below_near")!=null)            this.ft_tolerance_below_near=Double.parseDouble(properties.getProperty(prefix+"ft_tolerance_below_near"));
		if (properties.getProperty(prefix+"ft_tolerance_above_far")!=null)             this.ft_tolerance_above_far=Double.parseDouble(properties.getProperty(prefix+"ft_tolerance_above_far"));
		if (properties.getProperty(prefix+"ft_tolerance_below_far")!=null)             this.ft_tolerance_below_far=Double.parseDouble(properties.getProperty(prefix+"ft_tolerance_below_far"));
		if (properties.getProperty(prefix+"ft_hor_vert_overlap")!=null)                this.ft_hor_vert_overlap=Integer.parseInt(properties.getProperty(prefix+"ft_hor_vert_overlap"));
		if (properties.getProperty(prefix+"ft_used_companions")!=null)                 this.ft_used_companions=Integer.parseInt(properties.getProperty(prefix+"ft_used_companions"));
		if (properties.getProperty(prefix+"ft_used_true_companions")!=null)            this.ft_used_true_companions=Integer.parseInt(properties.getProperty(prefix+"ft_used_true_companions"));



		if (properties.getProperty(prefix+"plPreferDisparity")!=null)             this.plPreferDisparity=Boolean.parseBoolean(properties.getProperty(prefix+"plPreferDisparity"));
		if (properties.getProperty(prefix+"plDispNorm")!=null)                    this.plDispNorm=Double.parseDouble(properties.getProperty(prefix+"plDispNorm"));
		if (properties.getProperty(prefix+"plFrontoTol")!=null)                   this.plFrontoTol=Double.parseDouble(properties.getProperty(prefix+"plFrontoTol"));
		if (properties.getProperty(prefix+"plFrontoRms")!=null)                   this.plFrontoRms=Double.parseDouble(properties.getProperty(prefix+"plFrontoRms"));
		if (properties.getProperty(prefix+"plFrontoOffs")!=null)                  this.plFrontoOffs=Double.parseDouble(properties.getProperty(prefix+"plFrontoOffs"));
		if (properties.getProperty(prefix+"PlFrontoPow")!=null)                   this.PlFrontoPow=Double.parseDouble(properties.getProperty(prefix+"PlFrontoPow"));

		if (properties.getProperty(prefix+"plBlurBinVert")!=null)                 this.plBlurBinVert=Double.parseDouble(properties.getProperty(prefix+"plBlurBinVert"));
		if (properties.getProperty(prefix+"plBlurBinHor")!=null)                  this.plBlurBinHor=Double.parseDouble(properties.getProperty(prefix+"plBlurBinHor"));
		if (properties.getProperty(prefix+"plMaxDiffVert")!=null)                 this.plMaxDiffVert=Double.parseDouble(properties.getProperty(prefix+"plMaxDiffVert"));
		if (properties.getProperty(prefix+"plMaxDiffHor")!=null)                  this.plMaxDiffHor=Double.parseDouble(properties.getProperty(prefix+"plMaxDiffHor"));
		if (properties.getProperty(prefix+"plInitPasses")!=null)                  this.plInitPasses=Integer.parseInt(properties.getProperty(prefix+"plInitPasses"));

		if (properties.getProperty(prefix+"plMinPoints")!=null)                   this.plMinPoints=Integer.parseInt(properties.getProperty(prefix+"plMinPoints"));
		if (properties.getProperty(prefix+"plTargetEigen")!=null)                 this.plTargetEigen=Double.parseDouble(properties.getProperty(prefix+"plTargetEigen"));
		if (properties.getProperty(prefix+"plFractOutliers")!=null)               this.plFractOutliers=Double.parseDouble(properties.getProperty(prefix+"plFractOutliers"));
		if (properties.getProperty(prefix+"plMaxOutliers")!=null)                 this.plMaxOutliers=Integer.parseInt(properties.getProperty(prefix+"plMaxOutliers"));
		if (properties.getProperty(prefix+"plMinStrength")!=null)                 this.plMinStrength=Double.parseDouble(properties.getProperty(prefix+"plMinStrength"));
		if (properties.getProperty(prefix+"plMaxEigen")!=null)                    this.plMaxEigen=Double.parseDouble(properties.getProperty(prefix+"plMaxEigen"));
		if (properties.getProperty(prefix+"plEigenFloor")!=null)                  this.plEigenFloor=Double.parseDouble(properties.getProperty(prefix+"plEigenFloor"));
		if (properties.getProperty(prefix+"plEigenStick")!=null)                  this.plEigenStick=Double.parseDouble(properties.getProperty(prefix+"plEigenStick"));
		if (properties.getProperty(prefix+"plBadPlate")!=null)                    this.plBadPlate=Double.parseDouble(properties.getProperty(prefix+"plBadPlate"));
		if (properties.getProperty(prefix+"plDbgMerge")!=null)                    this.plDbgMerge=Boolean.parseBoolean(properties.getProperty(prefix+"plDbgMerge"));
		if (properties.getProperty(prefix+"plWorstWorsening")!=null)              this.plWorstWorsening=Double.parseDouble(properties.getProperty(prefix+"plWorstWorsening"));
		if (properties.getProperty(prefix+"plWorstWorsening2")!=null)             this.plWorstWorsening2=Double.parseDouble(properties.getProperty(prefix+"plWorstWorsening2"));
		if (properties.getProperty(prefix+"plWorstEq")!=null)                     this.plWorstEq=Double.parseDouble(properties.getProperty(prefix+"plWorstEq"));
		if (properties.getProperty(prefix+"plWorstEq2")!=null)                    this.plWorstEq2=Double.parseDouble(properties.getProperty(prefix+"plWorstEq2"));
		if (properties.getProperty(prefix+"plOKMergeEigen")!=null)                this.plOKMergeEigen=Double.parseDouble(properties.getProperty(prefix+"plOKMergeEigen"));
		if (properties.getProperty(prefix+"plMaxWorldSin2")!=null)                this.plMaxWorldSin2=Double.parseDouble(properties.getProperty(prefix+"plMaxWorldSin2"));
		if (properties.getProperty(prefix+"pl2dForSin")!=null)                    this.pl2dForSin=Double.parseDouble(properties.getProperty(prefix+"pl2dForSin"));
		if (properties.getProperty(prefix+"plWeakWorsening")!=null)               this.plWeakWorsening=Double.parseDouble(properties.getProperty(prefix+"plWeakWorsening"));
		if (properties.getProperty(prefix+"plMaxOverlap")!=null)                  this.plMaxOverlap=Double.parseDouble(properties.getProperty(prefix+"plMaxOverlap"));

		if (properties.getProperty(prefix+"plWeakWeight")!=null)                  this.plWeakWeight=Double.parseDouble(properties.getProperty(prefix+"plWeakWeight"));
		if (properties.getProperty(prefix+"plWeakEigen")!=null)                   this.plWeakEigen=Double.parseDouble(properties.getProperty(prefix+"plWeakEigen"));
		if (properties.getProperty(prefix+"plWeakWeight2")!=null)                 this.plWeakWeight2=Double.parseDouble(properties.getProperty(prefix+"plWeakWeight2"));
		if (properties.getProperty(prefix+"plWeakEigen2")!=null)                  this.plWeakEigen2=Double.parseDouble(properties.getProperty(prefix+"plWeakEigen2"));
		if (properties.getProperty(prefix+"plSumThick")!=null)                    this.plSumThick=Double.parseDouble(properties.getProperty(prefix+"plSumThick"));
		if (properties.getProperty(prefix+"plNeNeibCost")!=null)                  this.plNeNeibCost=Double.parseDouble(properties.getProperty(prefix+"plNeNeibCost"));
		if (properties.getProperty(prefix+"plNeOwn")!=null)                       this.plNeOwn=Double.parseDouble(properties.getProperty(prefix+"plNeOwn"));

		if (properties.getProperty(prefix+"plExNeibCost")!=null)                  this.plExNeibCost=Double.parseDouble(properties.getProperty(prefix+"plExNeibCost"));
		if (properties.getProperty(prefix+"plExNeibSmooth")!=null)                this.plExNeibSmooth=Double.parseDouble(properties.getProperty(prefix+"plExNeibSmooth"));
		if (properties.getProperty(prefix+"plMergeCostStar")!=null)               this.plMergeCostStar=Double.parseDouble(properties.getProperty(prefix+"plMergeCostStar"));
		if (properties.getProperty(prefix+"plMergeCost")!=null)                   this.plMergeCost=Double.parseDouble(properties.getProperty(prefix+"plMergeCost"));

		if (properties.getProperty(prefix+"plConflMerge")!=null)                  this.plConflMerge=Boolean.parseBoolean(properties.getProperty(prefix+"plConflMerge"));
		if (properties.getProperty(prefix+"plConflRelax")!=null)                  this.plConflRelax=Double.parseDouble(properties.getProperty(prefix+"plConflRelax"));
		if (properties.getProperty(prefix+"plConflSngl")!=null)                   this.plConflSngl=Boolean.parseBoolean(properties.getProperty(prefix+"plConflSngl"));
		if (properties.getProperty(prefix+"plConflSnglPair")!=null)               this.plConflSnglPair=Boolean.parseBoolean(properties.getProperty(prefix+"plConflSnglPair"));

		if (properties.getProperty(prefix+"plWeakFgStrength")!=null)              this.plWeakFgStrength=Double.parseDouble(properties.getProperty(prefix+"plWeakFgStrength"));
		if (properties.getProperty(prefix+"plWeakFgOutliers")!=null)              this.plWeakFgOutliers=Integer.parseInt(properties.getProperty(prefix+"plWeakFgOutliers"));
		if (properties.getProperty(prefix+"plWeakFgRelax")!=null)                 this.plWeakFgRelax=Double.parseDouble(properties.getProperty(prefix+"plWeakFgRelax"));

		if (properties.getProperty(prefix+"plThickWorld")!=null)                  this.plThickWorld=Double.parseDouble(properties.getProperty(prefix+"plThickWorld"));
		if (properties.getProperty(prefix+"plThickWorldConfl")!=null)             this.plThickWorldConfl=Double.parseDouble(properties.getProperty(prefix+"plThickWorldConfl"));
		if (properties.getProperty(prefix+"plRelaxComplete")!=null)               this.plRelaxComplete=Double.parseDouble(properties.getProperty(prefix+"plRelaxComplete"));
		if (properties.getProperty(prefix+"plRelaxComplete2")!=null)              this.plRelaxComplete2=Double.parseDouble(properties.getProperty(prefix+"plRelaxComplete2"));

		if (properties.getProperty(prefix+"plMaxZRatio")!=null)                   this.plMaxZRatio=Double.parseDouble(properties.getProperty(prefix+"plMaxZRatio"));
		if (properties.getProperty(prefix+"plMaxDisp")!=null)                     this.plMaxDisp=Double.parseDouble(properties.getProperty(prefix+"plMaxDisp"));
		if (properties.getProperty(prefix+"plCutTail")!=null)                     this.plCutTail=Double.parseDouble(properties.getProperty(prefix+"plCutTail"));
		if (properties.getProperty(prefix+"plMinTail")!=null)                     this.plMinTail=Double.parseDouble(properties.getProperty(prefix+"plMinTail"));

		if (properties.getProperty(prefix+"plDiscrEn")!=null)                     this.plDiscrEn=Boolean.parseBoolean(properties.getProperty(prefix+"plDiscrEn"));
		if (properties.getProperty(prefix+"plDiscrTolerance")!=null)              this.plDiscrTolerance=Double.parseDouble(properties.getProperty(prefix+"plDiscrTolerance"));
		if (properties.getProperty(prefix+"plDiscrDispRange")!=null)              this.plDiscrDispRange=Double.parseDouble(properties.getProperty(prefix+"plDiscrDispRange"));
		if (properties.getProperty(prefix+"plDiscrSteps")!=null)                  this.plDiscrSteps=Integer.parseInt(properties.getProperty(prefix+"plDiscrSteps"));
		//  			if (properties.getProperty(prefix+"plDiscrVariants")!=null)               this.plDiscrVariants=Integer.parseInt(properties.getProperty(prefix+"plDiscrVariants"));
		if (properties.getProperty(prefix+"plDiscrMode")!=null)                   this.plDiscrMode=Integer.parseInt(properties.getProperty(prefix+"plDiscrMode"));
		if (properties.getProperty(prefix+"plDiscrVarFloor")!=null)               this.plDiscrVarFloor=Double.parseDouble(properties.getProperty(prefix+"plDiscrVarFloor"));
		if (properties.getProperty(prefix+"plDiscrSigma")!=null)                  this.plDiscrSigma=Double.parseDouble(properties.getProperty(prefix+"plDiscrSigma"));
		if (properties.getProperty(prefix+"plDiscrBlur")!=null)                   this.plDiscrBlur=Double.parseDouble(properties.getProperty(prefix+"plDiscrBlur"));
		if (properties.getProperty(prefix+"plDiscrExclusivity")!=null)            this.plDiscrExclusivity=Double.parseDouble(properties.getProperty(prefix+"plDiscrExclusivity"));
		if (properties.getProperty(prefix+"plDiscrExclus2")!=null)                this.plDiscrExclus2=Double.parseDouble(properties.getProperty(prefix+"plDiscrExclus2"));
		if (properties.getProperty(prefix+"plDiscrStrict")!=null)                 this.plDiscrStrict=Boolean.parseBoolean(properties.getProperty(prefix+"plDiscrStrict"));
		if (properties.getProperty(prefix+"plDiscrCorrMax")!=null)                this.plDiscrCorrMax=Double.parseDouble(properties.getProperty(prefix+"plDiscrCorrMax"));
		if (properties.getProperty(prefix+"plDiscrCorrMerge")!=null)              this.plDiscrCorrMerge=Double.parseDouble(properties.getProperty(prefix+"plDiscrCorrMerge"));
		if (properties.getProperty(prefix+"plDiscrSteal")!=null)                  this.plDiscrSteal=Integer.parseInt(properties.getProperty(prefix+"plDiscrSteal"));
		if (properties.getProperty(prefix+"plDiscrGrown")!=null)                  this.plDiscrGrown=Integer.parseInt(properties.getProperty(prefix+"plDiscrGrown"));
		if (properties.getProperty(prefix+"plDiscrXMedian")!=null)                this.plDiscrXMedian=Double.parseDouble(properties.getProperty(prefix+"plDiscrXMedian"));

		if (properties.getProperty(prefix+"plCostDist")!=null)                    this.plCostDist=Double.parseDouble(properties.getProperty(prefix+"plCostDist"));
		if (properties.getProperty(prefix+"plCostKrq")!=null)                     this.plCostKrq=Double.parseDouble(properties.getProperty(prefix+"plCostKrq"));
		if (properties.getProperty(prefix+"plCostKrqEq")!=null)                   this.plCostKrqEq=Double.parseDouble(properties.getProperty(prefix+"plCostKrqEq"));
		if (properties.getProperty(prefix+"plCostWrq")!=null)                     this.plCostWrq=Double.parseDouble(properties.getProperty(prefix+"plCostWrq"));
		if (properties.getProperty(prefix+"plCostWrqEq")!=null)                   this.plCostWrqEq=Double.parseDouble(properties.getProperty(prefix+"plCostWrqEq"));
		if (properties.getProperty(prefix+"plCostSin2")!=null)                    this.plCostSin2=Double.parseDouble(properties.getProperty(prefix+"plCostSin2"));
		if (properties.getProperty(prefix+"plCostRdist2")!=null)                  this.plCostRdist2=Double.parseDouble(properties.getProperty(prefix+"plCostRdist2"));

		if (properties.getProperty(prefix+"plConflDualTri")!=null)                this.plConflDualTri=Boolean.parseBoolean(properties.getProperty(prefix+"plConflDualTri"));
		if (properties.getProperty(prefix+"plConflMulti")!=null)                  this.plConflMulti=Boolean.parseBoolean(properties.getProperty(prefix+"plConflMulti"));
		if (properties.getProperty(prefix+"plConflDiag")!=null)                   this.plConflDiag=Boolean.parseBoolean(properties.getProperty(prefix+"plConflDiag"));
		if (properties.getProperty(prefix+"plConflStar")!=null)                   this.plConflStar=Boolean.parseBoolean(properties.getProperty(prefix+"plConflStar"));
		if (properties.getProperty(prefix+"plStarSteps")!=null)                   this.plStarSteps=Integer.parseInt(properties.getProperty(prefix+"plStarSteps"));
		if (properties.getProperty(prefix+"plStarOrtho")!=null)                   this.plStarOrtho=Double.parseDouble(properties.getProperty(prefix+"plStarOrtho"));
		if (properties.getProperty(prefix+"plStarDiag")!=null)                    this.plStarDiag=Double.parseDouble(properties.getProperty(prefix+"plStarDiag"));
		if (properties.getProperty(prefix+"plStarPwr")!=null)                     this.plStarPwr=Double.parseDouble(properties.getProperty(prefix+"plStarPwr"));
		if (properties.getProperty(prefix+"plStarWeightPwr")!=null)               this.plStarWeightPwr=Double.parseDouble(properties.getProperty(prefix+"plStarWeightPwr"));
		if (properties.getProperty(prefix+"plWeightToDens")!=null)                this.plWeightToDens=Double.parseDouble(properties.getProperty(prefix+"plWeightToDens"));
		if (properties.getProperty(prefix+"plStarValPwr")!=null)                  this.plStarValPwr=Double.parseDouble(properties.getProperty(prefix+"plStarValPwr"));
		if (properties.getProperty(prefix+"plDblTriLoss")!=null)                  this.plDblTriLoss=Double.parseDouble(properties.getProperty(prefix+"plDblTriLoss"));
		if (properties.getProperty(prefix+"plNewConfl")!=null)                    this.plNewConfl=Boolean.parseBoolean(properties.getProperty(prefix+"plNewConfl"));
		if (properties.getProperty(prefix+"plMaxChanges")!=null)                  this.plMaxChanges=Integer.parseInt(properties.getProperty(prefix+"plMaxChanges"));

		if (properties.getProperty(prefix+"plMutualOnly")!=null)                  this.plMutualOnly=Boolean.parseBoolean(properties.getProperty(prefix+"plMutualOnly"));
		if (properties.getProperty(prefix+"plFillSquares")!=null)                 this.plFillSquares=Boolean.parseBoolean(properties.getProperty(prefix+"plFillSquares"));
		if (properties.getProperty(prefix+"plCutCorners")!=null)                  this.plCutCorners=Boolean.parseBoolean(properties.getProperty(prefix+"plCutCorners"));
		if (properties.getProperty(prefix+"plHypotenuse")!=null)                  this.plHypotenuse=Boolean.parseBoolean(properties.getProperty(prefix+"plHypotenuse"));

		if (properties.getProperty(prefix+"plPull")!=null)                        this.plPull=Double.parseDouble(properties.getProperty(prefix+"plPull"));
		if (properties.getProperty(prefix+"plNormPow")!=null)                     this.plNormPow=Double.parseDouble(properties.getProperty(prefix+"plNormPow"));
		if (properties.getProperty(prefix+"plIterations")!=null)                  this.plIterations=Integer.parseInt(properties.getProperty(prefix+"plIterations"));
		if (properties.getProperty(prefix+"plStopBad")!=null)                     this.plStopBad=Boolean.parseBoolean(properties.getProperty(prefix+"plStopBad"));
		if (properties.getProperty(prefix+"plPrecision")!=null)                   this.plPrecision=Integer.parseInt(properties.getProperty(prefix+"plPrecision"));

		if (properties.getProperty(prefix+"plSplitPull")!=null)                   this.plSplitPull=Double.parseDouble(properties.getProperty(prefix+"plSplitPull"));
		if (properties.getProperty(prefix+"plSplitMinNeib")!=null)                this.plSplitMinNeib=Integer.parseInt(properties.getProperty(prefix+"plSplitMinNeib"));
		if (properties.getProperty(prefix+"plSplitMinWeight")!=null)              this.plSplitMinWeight=Double.parseDouble(properties.getProperty(prefix+"plSplitMinWeight"));
		if (properties.getProperty(prefix+"plSplitMinQuality")!=null)             this.plSplitMinQuality=Double.parseDouble(properties.getProperty(prefix+"plSplitMinQuality"));
		if (properties.getProperty(prefix+"plSplitApply")!=null)                  this.plSplitApply=Boolean.parseBoolean(properties.getProperty(prefix+"plSplitApply"));
		if (properties.getProperty(prefix+"plNonExclusive")!=null)                this.plNonExclusive=Boolean.parseBoolean(properties.getProperty(prefix+"plNonExclusive"));
		if (properties.getProperty(prefix+"plUseOtherPlanes")!=null)            this.plUseOtherPlanes=Boolean.parseBoolean(properties.getProperty(prefix+"plUseOtherPlanes"));
		if (properties.getProperty(prefix+"plAllowParallel")!=null)             this.plAllowParallel=Boolean.parseBoolean(properties.getProperty(prefix+"plAllowParallel"));
		if (properties.getProperty(prefix+"plMaxDiff")!=null)                   this.plMaxDiff=Double.parseDouble(properties.getProperty(prefix+"plMaxDiff"));
		if (properties.getProperty(prefix+"plOtherDiff")!=null)                 this.plOtherDiff=Double.parseDouble(properties.getProperty(prefix+"plOtherDiff"));
		if (properties.getProperty(prefix+"plSplitXY")!=null)                   this.plSplitXY=Boolean.parseBoolean(properties.getProperty(prefix+"plSplitXY"));
		if (properties.getProperty(prefix+"plSplitXYTolerance")!=null)          this.plSplitXYTolerance=Double.parseDouble(properties.getProperty(prefix+"plSplitXYTolerance"));

		if (properties.getProperty(prefix+"plFuse")!=null)                      this.plFuse=Boolean.parseBoolean(properties.getProperty(prefix+"plFuse"));
		if (properties.getProperty(prefix+"plKeepOrphans")!=null)               this.plKeepOrphans=Boolean.parseBoolean(properties.getProperty(prefix+"plKeepOrphans"));
		if (properties.getProperty(prefix+"plMinOrphan")!=null)                 this.plMinOrphan=Double.parseDouble(properties.getProperty(prefix+"plMinOrphan"));

		if (properties.getProperty(prefix+"plSnapDispAny")!=null)               this.plSnapDispAny=Double.parseDouble(properties.getProperty(prefix+"plSnapDispAny"));
		if (properties.getProperty(prefix+"plSnapStrengthAny")!=null)           this.plSnapStrengthAny=Double.parseDouble(properties.getProperty(prefix+"plSnapStrengthAny"));
		if (properties.getProperty(prefix+"plSnapNegAny")!=null)                this.plSnapNegAny=Double.parseDouble(properties.getProperty(prefix+"plSnapNegAny"));
		if (properties.getProperty(prefix+"plSnapDispMax")!=null)               this.plSnapDispMax=Double.parseDouble(properties.getProperty(prefix+"plSnapDispMax"));
		if (properties.getProperty(prefix+"plSnapDispWeight")!=null)            this.plSnapDispWeight=Double.parseDouble(properties.getProperty(prefix+"plSnapDispWeight"));
		if (properties.getProperty(prefix+"plSnapZeroMode")!=null)              this.plPrecision=Integer.parseInt(properties.getProperty(prefix+"plSnapZeroMode"));

		if (properties.getProperty(prefix+"msUseSel")!=null)                    this.msUseSel=Boolean.parseBoolean(properties.getProperty(prefix+"msUseSel"));
		if (properties.getProperty(prefix+"msDivideByArea")!=null)              this.msDivideByArea=Boolean.parseBoolean(properties.getProperty(prefix+"msDivideByArea"));
		if (properties.getProperty(prefix+"msScaleProj")!=null)                 this.msScaleProj=Double.parseDouble(properties.getProperty(prefix+"msScaleProj"));
		if (properties.getProperty(prefix+"msFractUni")!=null)                  this.msFractUni=Double.parseDouble(properties.getProperty(prefix+"msFractUni"));


		if (properties.getProperty(prefix+"tsNoEdge")!=null)                    this.tsNoEdge=Boolean.parseBoolean(properties.getProperty(prefix+"tsNoEdge"));
		if (properties.getProperty(prefix+"tsUseCenter")!=null)                 this.tsUseCenter=Boolean.parseBoolean(properties.getProperty(prefix+"tsUseCenter"));
		if (properties.getProperty(prefix+"tsMaxDiff")!=null)                   this.tsMaxDiff=Double.parseDouble(properties.getProperty(prefix+"tsMaxDiff"));
		if (properties.getProperty(prefix+"tsMinDiffOther")!=null)              this.tsMinDiffOther=Double.parseDouble(properties.getProperty(prefix+"tsMinDiffOther"));
		if (properties.getProperty(prefix+"tsMinStrength")!=null)               this.tsMinStrength=Double.parseDouble(properties.getProperty(prefix+"tsMinStrength"));
		if (properties.getProperty(prefix+"tsMaxStrength")!=null)               this.tsMaxStrength=Double.parseDouble(properties.getProperty(prefix+"tsMaxStrength"));
		if (properties.getProperty(prefix+"tsMinSurface")!=null)                this.tsMinSurface=Double.parseDouble(properties.getProperty(prefix+"tsMinSurface"));
		if (properties.getProperty(prefix+"tsMoveDirs")!=null)                  this.tsMoveDirs=Integer.parseInt(properties.getProperty(prefix+"tsMoveDirs"));
		if (properties.getProperty(prefix+"tsSurfStrPow")!=null)                this.tsSurfStrPow=Double.parseDouble(properties.getProperty(prefix+"tsSurfStrPow"));
		if (properties.getProperty(prefix+"tsAddStrength")!=null)               this.tsAddStrength=Double.parseDouble(properties.getProperty(prefix+"tsAddStrength"));
		if (properties.getProperty(prefix+"tsSigma")!=null)                     this.tsSigma=Double.parseDouble(properties.getProperty(prefix+"tsSigma"));
		if (properties.getProperty(prefix+"tsNSigma")!=null)                    this.tsNSigma=Double.parseDouble(properties.getProperty(prefix+"tsNSigma"));
		if (properties.getProperty(prefix+"tsMinPull")!=null)                   this.tsMinPull=Double.parseDouble(properties.getProperty(prefix+"tsMinPull"));
		if (properties.getProperty(prefix+"tsMinAdvantage")!=null)              this.tsMinAdvantage=Double.parseDouble(properties.getProperty(prefix+"tsMinAdvantage"));

		if (properties.getProperty(prefix+"tsClustSize")!=null)                 this.tsClustSize=Integer.parseInt(properties.getProperty(prefix+"tsClustSize"));
		if (properties.getProperty(prefix+"tsClustWeight")!=null)               this.tsClustWeight=Double.parseDouble(properties.getProperty(prefix+"tsClustWeight"));
		if (properties.getProperty(prefix+"tsMinNeib")!=null)                   this.tsMinNeib=Integer.parseInt(properties.getProperty(prefix+"tsMinNeib"));
		if (properties.getProperty(prefix+"tsMaxSurStrength")!=null)            this.tsMaxSurStrength=Double.parseDouble(properties.getProperty(prefix+"tsMaxSurStrength"));
		if (properties.getProperty(prefix+"tsCountDis")!=null)                  this.tsCountDis=Boolean.parseBoolean(properties.getProperty(prefix+"tsCountDis"));

		if (properties.getProperty(prefix+"tsEnPlaneSeed")!=null)               this.tsEnPlaneSeed=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnPlaneSeed"));
		if (properties.getProperty(prefix+"tsEnOnly")!=null)                    this.tsEnOnly=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnOnly"));
		if (properties.getProperty(prefix+"tsEnGrow")!=null)                    this.tsEnGrow=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnGrow"));
		if (properties.getProperty(prefix+"tsGrowStrength")!=null)              this.tsGrowStrength=Double.parseDouble(properties.getProperty(prefix+"tsGrowStrength"));
		if (properties.getProperty(prefix+"tsGrowStrong")!=null)                this.tsGrowStrong=Boolean.parseBoolean(properties.getProperty(prefix+"tsGrowStrong"));
		if (properties.getProperty(prefix+"tsContStrength")!=null)              this.tsGrowStrength=Double.parseDouble(properties.getProperty(prefix+"tsContStrength"));
		if (properties.getProperty(prefix+"tsContDiff")!=null)                  this.tsContDiff=Double.parseDouble(properties.getProperty(prefix+"tsContDiff"));

		if (properties.getProperty(prefix+"tsEnSingle")!=null)                  this.tsEnSingle=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnSingle"));
		if (properties.getProperty(prefix+"tsEnMulti")!=null)                   this.tsEnMulti=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnMulti"));
		if (properties.getProperty(prefix+"tsRemoveWeak1")!=null)               this.tsRemoveWeak1=Boolean.parseBoolean(properties.getProperty(prefix+"tsRemoveWeak1"));
		if (properties.getProperty(prefix+"tsGrowSurround")!=null)              this.tsGrowSurround=Boolean.parseBoolean(properties.getProperty(prefix+"tsGrowSurround"));
		if (properties.getProperty(prefix+"tsRemoveWeak2")!=null)               this.tsRemoveWeak2=Boolean.parseBoolean(properties.getProperty(prefix+"tsRemoveWeak2"));

		if (properties.getProperty(prefix+"tsLoopMulti")!=null)                 this.tsLoopMulti=Boolean.parseBoolean(properties.getProperty(prefix+"tsLoopMulti"));
		if (properties.getProperty(prefix+"tsShow")!=null)                      this.tsShow=Boolean.parseBoolean(properties.getProperty(prefix+"tsShow"));
		if (properties.getProperty(prefix+"tsNumClust")!=null)                  this.tsNumClust=Integer.parseInt(properties.getProperty(prefix+"tsNumClust"));

		if (properties.getProperty(prefix+"tsConsensMode")!=null)               this.tsConsensMode=Integer.parseInt(properties.getProperty(prefix+"tsConsensMode"));
		if (properties.getProperty(prefix+"tsConsensAgree")!=null)              this.tsConsensAgree=Integer.parseInt(properties.getProperty(prefix+"tsConsensAgree"));

		if (properties.getProperty(prefix+"taMinFgBg")!=null)                   this.taMinFgBg=Double.parseDouble(properties.getProperty(prefix+"taMinFgBg"));
		if (properties.getProperty(prefix+"taMinFgEdge")!=null)                 this.taMinFgEdge=Double.parseDouble(properties.getProperty(prefix+"taMinFgEdge"));
		if (properties.getProperty(prefix+"taMinColSep")!=null)                 this.taMinColSep=Double.parseDouble(properties.getProperty(prefix+"taMinColSep"));
		if (properties.getProperty(prefix+"taMinColDiff")!=null)                this.taMinColDiff=Double.parseDouble(properties.getProperty(prefix+"taMinColDiff"));
		if (properties.getProperty(prefix+"taOutlier")!=null)                   this.taOutlier=Double.parseDouble(properties.getProperty(prefix+"taOutlier"));
		if (properties.getProperty(prefix+"taDiffPwr")!=null)                   this.taDiffPwr=Double.parseDouble(properties.getProperty(prefix+"taDiffPwr"));
		if (properties.getProperty(prefix+"taBestPwr")!=null)                   this.taBestPwr=Double.parseDouble(properties.getProperty(prefix+"taBestPwr"));
		if (properties.getProperty(prefix+"taDiff9Pwr")!=null)                  this.taDiff9Pwr=Double.parseDouble(properties.getProperty(prefix+"taDiff9Pwr"));
		if (properties.getProperty(prefix+"taColSigma")!=null)                  this.taDiff9Pwr=Double.parseDouble(properties.getProperty(prefix+"taColSigma"));
		if (properties.getProperty(prefix+"taColFraction")!=null)               this.taColFraction=Double.parseDouble(properties.getProperty(prefix+"taColFraction"));

		if (properties.getProperty(prefix+"taCostEmpty")!=null)                 this.taCostEmpty=Double.parseDouble(properties.getProperty(prefix+"taCostEmpty"));
		if (properties.getProperty(prefix+"taCostNoLink")!=null)                this.taCostNoLink=Double.parseDouble(properties.getProperty(prefix+"taCostNoLink"));
		if (properties.getProperty(prefix+"taCostSwitch")!=null)                this.taCostSwitch=Double.parseDouble(properties.getProperty(prefix+"taCostSwitch"));
		if (properties.getProperty(prefix+"taCostColor")!=null)                 this.taCostColor=Double.parseDouble(properties.getProperty(prefix+"taCostColor"));
		if (properties.getProperty(prefix+"taCostDiff")!=null)                  this.taCostDiff=Double.parseDouble(properties.getProperty(prefix+"taCostDiff"));
		if (properties.getProperty(prefix+"taCostDiffBest")!=null)              this.taCostDiffBest=Double.parseDouble(properties.getProperty(prefix+"taCostDiffBest"));
		if (properties.getProperty(prefix+"taCostDiff9")!=null)                 this.taCostDiff9=Double.parseDouble(properties.getProperty(prefix+"taCostDiff9"));
		if (properties.getProperty(prefix+"taCostWeakFgnd")!=null)              this.taCostWeakFgnd=Double.parseDouble(properties.getProperty(prefix+"taCostWeakFgnd"));
		if (properties.getProperty(prefix+"taCostFlaps")!=null)                 this.taCostFlaps=Double.parseDouble(properties.getProperty(prefix+"taCostFlaps"));
		if (properties.getProperty(prefix+"taCostMismatch")!=null)              this.taCostMismatch=Double.parseDouble(properties.getProperty(prefix+"taCostMismatch"));

		if (properties.getProperty(prefix+"taEnEmpty")!=null)                   this.taEnEmpty=Boolean.parseBoolean(properties.getProperty(prefix+"taEnEmpty"));
		if (properties.getProperty(prefix+"taEnNoLink")!=null)                  this.taEnNoLink=Boolean.parseBoolean(properties.getProperty(prefix+"taEnNoLink"));
		if (properties.getProperty(prefix+"taEnSwitch")!=null)                  this.taEnSwitch=Boolean.parseBoolean(properties.getProperty(prefix+"taEnSwitch"));
		if (properties.getProperty(prefix+"taEnColor")!=null)                   this.taEnColor=Boolean.parseBoolean(properties.getProperty(prefix+"taEnColor"));
		if (properties.getProperty(prefix+"taEnDiff")!=null)                    this.taEnDiff=Boolean.parseBoolean(properties.getProperty(prefix+"taEnDiff"));
		if (properties.getProperty(prefix+"taEnDiffBest")!=null)                this.taEnDiffBest=Boolean.parseBoolean(properties.getProperty(prefix+"taEnDiffBest"));
		if (properties.getProperty(prefix+"taEnDiff9")!=null)                   this.taEnDiff9=Boolean.parseBoolean(properties.getProperty(prefix+"taEnDiff9"));
		if (properties.getProperty(prefix+"taEnWeakFgnd")!=null)                this.taEnWeakFgnd=Boolean.parseBoolean(properties.getProperty(prefix+"taEnWeakFgnd"));
		if (properties.getProperty(prefix+"taEnFlaps")!=null)                   this.taEnFlaps=Boolean.parseBoolean(properties.getProperty(prefix+"taEnFlaps"));
		if (properties.getProperty(prefix+"taEnMismatch")!=null)                this.taEnMismatch=Boolean.parseBoolean(properties.getProperty(prefix+"taEnMismatch"));

		if (properties.getProperty(prefix+"gpu_corr_rad")!=null)                this.gpu_corr_rad=Integer.parseInt(properties.getProperty(prefix+"gpu_corr_rad"));
		if (properties.getProperty(prefix+"gpu_weight_r")!=null)                this.gpu_weight_r=Double.parseDouble(properties.getProperty(prefix+"gpu_weight_r"));
		if (properties.getProperty(prefix+"gpu_weight_b")!=null)                this.gpu_weight_b=Double.parseDouble(properties.getProperty(prefix+"gpu_weight_b"));
		if (properties.getProperty(prefix+"gpu_sigma_r")!=null)                 this.gpu_sigma_r=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_r"));
		if (properties.getProperty(prefix+"gpu_sigma_b")!=null)                 this.gpu_sigma_b=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_b"));
		if (properties.getProperty(prefix+"gpu_sigma_g")!=null)                 this.gpu_sigma_g=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_g"));
		if (properties.getProperty(prefix+"gpu_sigma_m")!=null)                 this.gpu_sigma_m=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_m"));
		if (properties.getProperty(prefix+"gpu_sigma_rb_corr")!=null)           this.gpu_sigma_rb_corr=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_rb_corr"));
		if (properties.getProperty(prefix+"gpu_sigma_corr")!=null)              this.gpu_sigma_corr=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_corr"));
		if (properties.getProperty(prefix+"gpu_sigma_corr_m")!=null)            this.gpu_sigma_corr_m=Double.parseDouble(properties.getProperty(prefix+"gpu_sigma_corr_m"));
		if (properties.getProperty(prefix+"gpu_fatz")!=null)                    this.gpu_fatz=Double.parseDouble(properties.getProperty(prefix+"gpu_fatz"));
		if (properties.getProperty(prefix+"gpu_fatz_m")!=null)                  this.gpu_fatz_m=Double.parseDouble(properties.getProperty(prefix+"gpu_fatz_m"));

		if (properties.getProperty(prefix+"gpu_woi")!=null)                     this.gpu_woi=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_woi"));
		if (properties.getProperty(prefix+"gpu_woi_tx")!=null)                  this.gpu_woi_tx=Integer.parseInt(properties.getProperty(prefix+"gpu_woi_tx"));
		if (properties.getProperty(prefix+"gpu_woi_ty")!=null)                  this.gpu_woi_ty=Integer.parseInt(properties.getProperty(prefix+"gpu_woi_ty"));
		if (properties.getProperty(prefix+"gpu_woi_twidth")!=null)              this.gpu_woi_twidth=Integer.parseInt(properties.getProperty(prefix+"gpu_woi_twidth"));
		if (properties.getProperty(prefix+"gpu_woi_theight")!=null)             this.gpu_woi_theight=Integer.parseInt(properties.getProperty(prefix+"gpu_woi_theight"));
		if (properties.getProperty(prefix+"gpu_woi_round")!=null)               this.gpu_woi_round=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_woi_round"));
		if (properties.getProperty(prefix+"gpu_save_ports_xy")!=null)           this.gpu_save_ports_xy=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_save_ports_xy"));
		if (properties.getProperty(prefix+"gpu_show_jtextures")!=null)          this.gpu_show_jtextures=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_show_jtextures"));
		if (properties.getProperty(prefix+"gpu_show_extra")!=null)              this.gpu_show_extra=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_show_extra"));
		
		if (properties.getProperty(prefix+"gpu_use_main")!=null)                this.gpu_use_main=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_use_main"));
		if (properties.getProperty(prefix+"gpu_use_main_macro")!=null)          this.gpu_use_main_macro=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_use_main_macro"));
		if (properties.getProperty(prefix+"gpu_use_main_adjust")!=null)         this.gpu_use_main_adjust=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_use_main_adjust"));
		if (properties.getProperty(prefix+"gpu_use_aux")!=null)                 this.gpu_use_aux=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_use_aux"));
		if (properties.getProperty(prefix+"gpu_use_aux_macro")!=null)           this.gpu_use_aux_macro=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_use_aux_macro"));
		if (properties.getProperty(prefix+"gpu_use_aux_adjust")!=null)          this.gpu_use_aux_adjust=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_use_aux_adjust"));

		if (properties.getProperty(prefix+"debug_initial_discriminate")!=null)       this.debug_initial_discriminate=Boolean.parseBoolean(properties.getProperty(prefix+"debug_initial_discriminate"));
		if (properties.getProperty(prefix+"dbg_migrate")!=null)                      this.dbg_migrate=Boolean.parseBoolean(properties.getProperty(prefix+"dbg_migrate"));

		if (properties.getProperty(prefix+"dbg_early_exit")!=null)                   this.dbg_early_exit=Integer.parseInt(properties.getProperty(prefix+"dbg_early_exit"));
		if (properties.getProperty(prefix+"show_first_bg")!=null)                    this.show_first_bg=Boolean.parseBoolean(properties.getProperty(prefix+"show_first_bg"));

		if (properties.getProperty(prefix+"show_extrinsic")!=null)                   this.show_extrinsic=Boolean.parseBoolean(properties.getProperty(prefix+"show_extrinsic"));
		if (properties.getProperty(prefix+"show_ortho_combine")!=null)               this.show_ortho_combine=Boolean.parseBoolean(properties.getProperty(prefix+"show_ortho_combine"));
		if (properties.getProperty(prefix+"show_refine_supertiles")!=null)           this.show_refine_supertiles=Boolean.parseBoolean(properties.getProperty(prefix+"show_refine_supertiles"));
		if (properties.getProperty(prefix+"show_bgnd_nonbgnd")!=null)                this.show_bgnd_nonbgnd=Boolean.parseBoolean(properties.getProperty(prefix+"show_bgnd_nonbgnd"));
		if (properties.getProperty(prefix+"show_filter_scan")!=null)                 this.show_filter_scan=Boolean.parseBoolean(properties.getProperty(prefix+"show_filter_scan"));
		if (properties.getProperty(prefix+"show_combined")!=null)                    this.show_combined=Boolean.parseBoolean(properties.getProperty(prefix+"show_combined"));
		if (properties.getProperty(prefix+"show_unique")!=null)                      this.show_unique=Boolean.parseBoolean(properties.getProperty(prefix+"show_unique"));
		if (properties.getProperty(prefix+"show_histograms")!=null)                  this.show_histograms=Boolean.parseBoolean(properties.getProperty(prefix+"show_histograms"));
		if (properties.getProperty(prefix+"show_init_refine")!=null)                 this.show_init_refine=Boolean.parseBoolean(properties.getProperty(prefix+"show_init_refine"));
		if (properties.getProperty(prefix+"show_expand")!=null)                      this.show_expand=Boolean.parseBoolean(properties.getProperty(prefix+"show_expand"));
		if (properties.getProperty(prefix+"show_variant")!=null)                     this.show_variant=Boolean.parseBoolean(properties.getProperty(prefix+"show_variant"));
		if (properties.getProperty(prefix+"show_retry_far")!=null)                   this.show_retry_far=Boolean.parseBoolean(properties.getProperty(prefix+"show_retry_far"));
		if (properties.getProperty(prefix+"show_macro")!=null)                       this.show_macro=Boolean.parseBoolean(properties.getProperty(prefix+"show_macro"));
		if (properties.getProperty(prefix+"show_shells")!=null)                      this.show_shells=Boolean.parseBoolean(properties.getProperty(prefix+"show_shells"));
		if (properties.getProperty(prefix+"show_neighbors")!=null)                   this.show_neighbors=Boolean.parseBoolean(properties.getProperty(prefix+"show_neighbors"));
		if (properties.getProperty(prefix+"show_flaps_dirs")!=null)                  this.show_flaps_dirs=Boolean.parseBoolean(properties.getProperty(prefix+"show_flaps_dirs"));
		if (properties.getProperty(prefix+"show_first_clusters")!=null)              this.show_first_clusters=Boolean.parseBoolean(properties.getProperty(prefix+"show_first_clusters"));
		if (properties.getProperty(prefix+"show_planes")!=null)                      this.show_planes=Boolean.parseBoolean(properties.getProperty(prefix+"show_planes"));


		if (properties.getProperty(prefix+"vertical_xyz.x")!=null)               this.vertical_xyz[0]=Double.parseDouble(properties.getProperty(prefix+"vertical_xyz.x"));
		if (properties.getProperty(prefix+"vertical_xyz.y")!=null)               this.vertical_xyz[1]=Double.parseDouble(properties.getProperty(prefix+"vertical_xyz.y"));
		if (properties.getProperty(prefix+"vertical_xyz.z")!=null)               this.vertical_xyz[2]=Double.parseDouble(properties.getProperty(prefix+"vertical_xyz.z"));

		Set<String> ss = properties.stringPropertyNames();
		String full_prefix = prefix+Z_CORR_PREFIX;
		int li = full_prefix.length();
		for (String s:ss){
			if (s.indexOf(full_prefix) == 0){
				z_corr_map.put(s.substring(li), Double.parseDouble(properties.getProperty(s)));
			}
		}
		String full_prefix_infinity_distance = prefix+INFINITY_DISTANCE_PREFIX;
		int lii = full_prefix_infinity_distance.length();
		for (String s:ss){
			if (s.indexOf(full_prefix_infinity_distance) == 0){
				infinity_distace_map.put(s.substring(lii), Double.parseDouble(properties.getProperty(s)));
			}
		}
		img_dtt.getProperties (prefix+"_img_dtt", properties);
		mlfp.getProperties    (prefix+"_mlfp",    properties);
		rig.getProperties     (prefix+"_rig",     properties);
		poles.getProperties   (prefix+"_poles",   properties);
		lwir.getProperties    (prefix+"_lwir",   properties);
	}

	public boolean showJDialog() {
		//  			GenericDialog gd = new GenericDialog("Set CLT parameters");
		GenericJTabbedDialog gd = new GenericJTabbedDialog("Set CLT parameters",800,900);
		gd.addTab         ("General", "General parameters");
		gd.addNumericField("Nominal (rectilinear) disparity between side of square cameras (pix)",              this.disparity,  3,7,"pix",
				"Used when rendering 4 images");
		gd.addNumericField("Transform size (default 8)",                                                        this.transform_size,            0, 6, "pixels","Should always be 8");
		gd.addNumericField("Lapped transform window type (0- rectangular, 1 - sinus)",                          this.clt_window,                0);
		gd.addNumericField("shift_x",                                                                           this.shift_x,                   4);
		gd.addNumericField("shift_y",                                                                           this.shift_y,                   4);

		gd.addNumericField("Lazy eye cluster size",                                                             this.tileStep,                  0, 6, "tiles",
				"Process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters");
		gd.addNumericField("Bit mask - which of 4 transforms to combine after iclt",                            this.iclt_mask,                 0);
		gd.addNumericField("Tile X to extract (0..163)",                                                        this.tileX,                     0);
		gd.addNumericField("Tile Y to extract (0..122)",                                                        this.tileY,                     0);
		gd.addNumericField("dbg_mode: 0 - normal, +1 - no DCT/IDCT, just fold",                                 this.dbg_mode,                  0);
		gd.addNumericField("ishift_x: shift source image by this pixels left",                                  this.ishift_x,                  0);
		gd.addNumericField("ishift_y: shift source image by this pixels down",                                  this.ishift_y,                  0);
		gd.addNumericField("Modify phase correlation to prevent division by very small numbers",                this.fat_zero,                  4);
		gd.addNumericField("Modify phase correlation for monochrome images",                                    this.fat_zero_mono,             4);
		gd.addNumericField("LPF correlarion sigma for Bayer color images",                                      this.corr_sigma,                4);
		gd.addNumericField("LPF correlarion sigma for monochrome images",                                       this.corr_sigma_mono,           4);

		gd.addNumericField("Scale all correlation strengths (to compensate correlation sigma) for main camera", this.scale_strength_main,       4);
		gd.addNumericField("Scale all correlation strengths (to compensate correlation sigma) for aux camera",  this.scale_strength_aux,        4);

		gd.addCheckbox    ("Normalize kernels",                                                                 this.norm_kern);
		gd.addCheckbox    ("Equalize green channel gain of the individual cnannels (bug fix for exposure)",           this.gain_equalize);
		gd.addCheckbox    ("Equalize R/G, B/G balance of the individual channels",                              this.colors_equalize);
		gd.addCheckbox    ("Skip saturated when adjusting gains",                                               this.nosat_equalize);
		gd.addNumericField("Saturation level of the most saturated color channel",                              this.sat_level,   4);
		gd.addNumericField("Do not use tiles with higher fraction of (near) saturated tiles",                   this.max_overexposure,   4);

		gd.addNumericField("Red gain in the center of sensor calibration R (instead of vignetting)",            this.novignetting_r,   4);
		gd.addNumericField("Green gain in the center of sensor calibration G (instead of vignetting)",          this.novignetting_g, 4);
		gd.addNumericField("Blue gain in the center of sensor calibration B (instead of vignetting)",           this.novignetting_b,  4);
		gd.addNumericField("Extra red correction to compensate for light temperature",                          this.scale_r,  4);
		gd.addNumericField("Extra green correction to compensate for light temperature",                        this.scale_g,  4);
		gd.addNumericField("Extra blue correction to compensate for light temperature",                         this.scale_b,  4);
		gd.addNumericField("Value (max) in vignetting data to correspond to 1x in the kernel",                  this.vignetting_max,      3);
		gd.addNumericField("Do not try to correct vignetting smaller than this fraction of max",                this.vignetting_range,  3);
		gd.addNumericField("Kernel step in pixels (has 1 kernel margin on each side)",                          this.kernel_step,            0);
		gd.addNumericField("Inverse distance to infinity (misalignment correction)",                            this.z_correction,  6);
		gd.addCheckbox    ("Perform correlation",                                                               this.correlate);
		gd.addNumericField("Bitmask of pairs to combine in the composite (top, bottom, left,righth)",           this.corr_mask,            0);
		gd.addCheckbox    ("Combine correlation with mirrored around disparity direction",                      this.corr_sym);
		gd.addCheckbox    ("Keep all partial correlations (otherwise - only combined one)",                     this.corr_keep);
		gd.addCheckbox    ("Show combined correlations",                                                        this.corr_show);
		gd.addCheckbox    ("Calculate per-pair X/Y variations of measured correlations ",                       this.corr_mismatch);
		gd.addNumericField("Add to pair correlation before multiplying by other pairs (between sum and product)",              this.corr_offset,  6);
		gd.addNumericField("Red to green correlation weight",                                                   this.corr_red,  4);
		gd.addNumericField("Blue to green correlation weight",                                                  this.corr_blue,  4);
		gd.addCheckbox    ("Normalize each correlation tile by rms",                                            this.corr_normalize);
		gd.addNumericField("Minimal correlation value to consider valid",                                       this.min_corr,  6);
		gd.addNumericField("Minimal correlation value to consider valid when normalizing results",              this.min_corr_normalized,  6);
		gd.addNumericField("Sigma for weights of points around global max to find fractional",                  this.max_corr_sigma,  3);
		gd.addNumericField("Maximal distance from int max to consider",                                         this.max_corr_radius,  3);
		gd.addNumericField("Calculated from correlation offset vs. actual one (not yet understood)",            this.corr_magic_scale,  3,6,"", "CM and poly");
		gd.addNumericField("Select all-pair correlation type to use 0 - CM, 1 - poly",                          this.corr_select,  0);

		gd.addTab         ("imageDtt", "Setup extra ImageDtt parameters - eventually all will be set that way");
		this.img_dtt.dialogQuestions(gd);
		gd.addTab         ("Rig", "Parameters for the wide baseline rig with two quad cameras");
		this.rig.dialogQuestions(gd);
		gd.addTab         ("Poles", "Parameters for processing vertical poles (such as street lights)");
		this.poles.dialogQuestions(gd);

		gd.addTab         ("vert/hor", "Enhance detection of horizontal/vertical features (when enh_ortho is enabled for tile");

		gd.addMessage("--- Enhance detection of horizontal/vertical features (when enh_ortho is enabled for tile ---");


		gd.addCheckbox    ("Double pass when masking center of mass to reduce preference for integer values",           this.max_corr_double);
		gd.addNumericField("Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial",             this.corr_mode,            0,4,"",
				"Only applies to textures");
		gd.addNumericField("Contrast of dotted border on correlation results",                                  this.corr_border_contrast,  6);

		gd.addTab         ("tileTasks", "tiles tasks (current tile_task_op = "+          this.tile_task_op+")");

		gd.addMessage("--- tiles tasks (current tile_task_op = "+          this.tile_task_op+") ---");
		gd.addCheckbox    ("Enhance ortho lines detection (enh_ortho)",                                              ImageDtt.getOrthoLines(          this.tile_task_op));
		gd.addCheckbox    ("Force disparity for image rendering (false - use found from tile correlation)",          ImageDtt.getForcedDisparity(          this.tile_task_op));
		gd.addNumericField("Bitmask of used images (1 - top left, 2 - top right, 4 - bottom left, 8 bottom right)",  ImageDtt.getImgMask(          this.tile_task_op),            0);
		gd.addNumericField("Bitmask of used pairs  (1 - top, 2 - bottom, 4 - left, 8 - right)",                      ImageDtt.getPairMask(          this.tile_task_op),            0);
		gd.addNumericField("Tile operations window left (in 8x8 tiles)",                                        this.tile_task_wl,            0);
		gd.addNumericField("Tile operations window top",                                                        this.tile_task_wt,            0);
		gd.addNumericField("Tile operations window width",                                                      this.tile_task_ww,            0);
		gd.addNumericField("Tile operations window height",                                                     this.tile_task_wh,            0);

		gd.addNumericField("Do not adjust for shot noise (~sqrt) if lower than this",                           this.min_shot,  4);
		gd.addNumericField("Scale when dividing by sqrt for shot noise compensation of pixel differences (<0 - disable)",           this.scale_shot,  4);

		gd.addNumericField("RMS difference from average to reduce weights (255 full scale image)",              this.diff_sigma,  4);
		gd.addNumericField("RMS difference from average in sigmas to discard channel",                          this.diff_threshold,  4);
		gd.addCheckbox    ("Gaussian as weight when averaging images (false - sharp all/nothing)",              this.diff_gauss);
		gd.addNumericField("Minimal number of channels to agree on a point (real number to work with fuzzy averages)",             this.min_agree,  2);
		gd.addCheckbox    ("Do not reduce average weight when only one image differes much from the average",           this.dust_remove);
		gd.addCheckbox    ("Use black for backdrop outside of the FOV",                                         this.black_back);
		gd.addCheckbox    ("Add port weights to RGBA stack (debug feature)",                                    this.keep_weights);
		gd.addCheckbox    ("Alpha channel: use center 8x8 (unchecked - treat same as RGB)",                     this.sharp_alpha);
		gd.addNumericField("Alpha channel 0.0 thereshold (lower - transparent)",                                this.alpha0,   3);
		gd.addNumericField("Alpha channel 1.0 threshold (higher - opaque)",                                     this.alpha1,   3);
		gd.addCheckbox    ("Generate shifted channel linear RGB stacks",                                        this.gen_chn_stacks);
		gd.addCheckbox    ("Generate shifted channel color image stack",                                        this.gen_chn_img);
		gd.addCheckbox    ("Generate shifted channel images and save with the model 'CLT process corr'",        this.gen_4_img);
		gd.addCheckbox    ("Show result RGBA before overlap combined",                                          this.show_nonoverlap);
		gd.addCheckbox    ("Show result RGBA",                                                                  this.show_overlap);
		gd.addCheckbox    ("Show result color",                                                                 this.show_rgba_color);
		gd.addCheckbox    ("Show disparity maps",                                                               this.show_map);
		gd.addCheckbox    ("Show correlation tiles",                                                            this.show_corr);

		gd.addNumericField("Disparity scan start value",                                                        this.disp_scan_start,  2);
		gd.addNumericField("Disparity scan step",                                                               this.disp_scan_step,  2);
		gd.addNumericField("Disparity scan number of disparity values to scan",                                 this.disp_scan_count,            0);

		gd.addTab         ("Manual", "camera fine correction: X/Y for images 0..3");

		gd.addMessage("--- camera fine correction: X/Y for images 0..3  ---");
		gd.addCheckbox    ("Debug infinity/lazy eye correction",                                                this.fine_dbg);
		gd.addNumericField("X 0",                                                                               this.fine_corr_x_0,  3);
		gd.addNumericField("Y 0",                                                                               this.fine_corr_y_0,  3);
		gd.addNumericField("X 1",                                                                               this.fine_corr_x_1,  3);
		gd.addNumericField("Y 1",                                                                               this.fine_corr_y_1,  3);
		gd.addNumericField("X 2",                                                                               this.fine_corr_x_2,  3);
		gd.addNumericField("Y 2",                                                                               this.fine_corr_y_2,  3);
		gd.addNumericField("X 3",                                                                               this.fine_corr_x_3,  3);
		gd.addNumericField("Y 3",                                                                               this.fine_corr_y_3,  3);
		gd.addCheckbox    ("Ignore manual pixel correction (set/reset automatically if enabled below)",             this.fine_corr_ignore);
		gd.addCheckbox    ("Apply and set to ignore manual pixel correction after extrinsics correction",           this.fine_corr_apply);

		gd.addNumericField("fcorr_radius",                                                                      this.fcorr_radius,  3);
		gd.addNumericField("Do not try to correct outside this fraction of width/hight",                        this.fcorr_min_strength,3);
		gd.addNumericField("Consider only tiles with absolute residual disparity lower than",                   this.fcorr_disp_diff,  3);
		gd.addCheckbox    ("Use quadratic polynomial for fine correction (false - only linear)",                this.fcorr_quadratic);
		gd.addCheckbox    ("Ignore current calculated fine correction (use manual only)",                       this.fcorr_ignore);

		gd.addNumericField("Minimal correlation strength to use for infinity correction",                       this.fcorr_inf_strength,3);
		gd.addNumericField("Disparity half-range for infinity",                                                 this.fcorr_inf_diff,  3);
		gd.addCheckbox    ("Use quadratic polynomial for infinity correction (false - only linear)",            this.fcorr_inf_quad);
		gd.addCheckbox    ("Correct infinity in vertical direction (false - only horizontal)",                  this.fcorr_inf_vert);


		gd.addCheckbox    ("Apply disparity correction to zero at infinity",                                    this.inf_disp_apply);
		gd.addNumericField("Re run disparity correction at infinity multiple times",                            this.inf_repeat,      0);

		//  			gd.addCheckbox    ("Apply lazy eye correction at infinity",                                             this.inf_mism_apply);
		gd.addNumericField("Infinity extraction - maximum iterations",                                          this.inf_iters,       0);
		gd.addNumericField("Coefficients maximal increment to exit iterations",                                 this.inf_final_diff,  6);
		gd.addNumericField("Include farther tiles than tolerance, but scale their weights",                     this.inf_far_pull,    3);

		gd.addTab         ("Infinity", "Infinity filter and Infinity histogram filter");

		gd.addMessage     ("--- Infinity filter ---");
		gd.addNumericField("Strength power",                                                                    this.inf_str_pow,     3);
		gd.addNumericField("Sample size (side of a square)",                                                    this.inf_smpl_side,   0);
		gd.addNumericField("Number after removing worst (should be >1)",                                        this.inf_smpl_num,    0);
		gd.addNumericField("Maximal RMS of the remaining tiles in a sample",                                    this.inf_smpl_rms,    3);
		gd.addMessage     ("--- Infinity histogram filter ---");
		gd.addNumericField("Square sample step (50% overlap)",                                                  this.ih_smpl_step,    0);
		gd.addNumericField("Histogram minimal disparity",                                                       this.ih_disp_min,     3);
		gd.addNumericField("Histogram disparity step",                                                          this.ih_disp_step,    3);
		gd.addNumericField("Histogram number of bins",                                                          this.ih_num_bins,     0);
		gd.addNumericField("Histogram Gaussian sigma (in disparity pixels)",                                    this.ih_sigma,        3);
		gd.addNumericField("Keep samples within this difference from farthest maximum",                         this.ih_max_diff,     3);
		gd.addNumericField("Minimal number of remaining samples",                                               this.ih_min_samples,  0);
		gd.addCheckbox    ("Replace samples with a single average with equal weight",                           this.ih_norm_center);
		gd.addCheckbox    ("Add disparity back to d{x,y}[i] (debug feature)",                                   this.inf_restore_disp);


		gd.addTab         ("Lazy eye", "Lazy eye parameters");

		gd.addCheckbox    ("Use 2020 LMA-based measurement of mismatch",                                        this.ly_lma_ers);
		gd.addMessage     ("--- main-to-aux depth map parameters ---");
		gd.addNumericField("Minimal reference (main) channel correlation strength",                             this.ly_gt_strength,        3);
		gd.addCheckbox    ("Use window for AUX tiles to reduce weight of the hi-res tiles near low-res tile boundaries", this.ly_gt_use_wnd);
		gd.addNumericField("Aux disparity thershold to split FG and BG (and disable AUX tile for adjustment)",  this.ly_gt_rms,        3);

		gd.addMessage     ("--- Lazy eye from GT ---");
		gd.addCheckbox    ("Adjust disparity using objects at infinity by changing individual tilt and azimuth ", this.lylw_inf_en," disable if there are no really far objects in the scene");
		gd.addCheckbox    ("Adjust azimuths and tilts",                                                           this.lylw_aztilt_en,"Adjust azimuths and tilts excluding those that change disparity");
		gd.addCheckbox    ("Adjust differential rolls",                                                           this.lylw_diff_roll_en,"Adjust differential rolls (3 of 4 rolls, keeping average roll)");
		gd.addCheckbox    ("Correct scales (focal length temperature? variations)",                               this.lylw_focalLength);
		gd.addCheckbox    ("Enable common roll adjustment (valid for high disparity range scans only)",           this.lylw_com_roll);
		gd.addNumericField("Manual parameter mask selection (0 use checkboxes above)",                            this.lylw_par_sel,  0, 5,"",
				"bit 0 - sym0, bit1 - sym1, ...");

		gd.addMessage     ("--- other LMA parameters ---");

		gd.addNumericField("Relative weight margins (0.0 - all 1.0, 1.0 sin^2",                                 this.ly_marg_fract,  8,3,"",
				"Reduce weigt of peripheral tiles");

		gd.addCheckbox    ("Calculate and apply lazy eye correction after disparity scan (poly or extrinsic), may repeat", this.ly_on_scan);
		gd.addCheckbox    ("Adjust disparity using objects at infinity by changing individual tilt and azimuth ",          this.ly_inf_en," disable if there are no really far objects in the scene");
		gd.addNumericField("Minimal number of clusters with forced disparity to use it (otherwise keep current)",this.ly_min_forced,  0);
		gd.addCheckbox    ("Adjust azimuths and tilts",                                                         this.ly_aztilt_en,"Adjust azimuths and tilts excluding those that change disparity");
		gd.addCheckbox    ("Adjust differential rolls",                                                         this.ly_diff_roll_en,"Adjust differential rolls (3 of 4 rolls, keeping average roll)");
		gd.addCheckbox    ("Correct scales (focal length temperature? variations)",                             this.ly_focalLength);
		gd.addCheckbox    ("Enable common roll adjustment (valid for high disparity range scans only)",         this.ly_com_roll);
		gd.addCheckbox    ("Enable ERS correction of the camera rotation",                                      this.ly_ers_rot);
		gd.addCheckbox    ("Enable ERS correction of the camera forward motion",                                this.ly_ers_forw);
		gd.addCheckbox    ("Enable ERS correction of the camera sideways motion",                               this.ly_ers_side);
		gd.addCheckbox    ("Enable ERS correction of the camera vertical motion",                               this.ly_ers_vert);
		gd.addNumericField("Manual parameter mask selection (0 use checkboxes above)",                          this.ly_par_sel,  0, 5,"",
				"bit 0 - sym0, bit1 - sym1, ...");
		gd.addNumericField("Debug level for lazy eye/ers processing",                                           this.ly_debug_level,  0, 5,"",
				"Active when global debug level > -1, 1 - min, 2 - lma steps, 3 - images");



		gd.addCheckbox    ("Equalize weights of right/left FoV",           this.ly_right_left,
				"Use this mode use with horizon visible in both FoV halves when gross infinity correction is needed");


		gd.addNumericField("Minimal tiles per quadrant (not counting the worst) tp proceed",                         this.ly_per_quad,  0);
		gd.addNumericField("Minimal tiles per quadrant (not counting the worst) tp proceed - fraction of all tiles", this.ly_per_quad_r,  3);
		gd.addNumericField("Minimal number of tiles at infinity to proceed",                                         this.ly_inf,  0);
		gd.addNumericField("Minimal number of tiles at infinity to proceed - fraction of all tiles",                 this.ly_inf_r,  3);
		gd.addNumericField("Minimal number of tiles at infinity to apply weight scaling",                            this.ly_inf_scale,  0);
		gd.addNumericField("Minimal number of tiles at infinity to apply weight scaling - fraction of all tiles",    this.ly_inf_scale_r,  3);

		gd.addNumericField("Relative weight of infinity calibration data",                                           this.ly_inf_frac,  3);

		gd.addNumericField("Maximal disparity to be treated as infinity when adjusting with the rig data",           this.ly_inf_max_disparity,  8,3,"pix",
				"Only used in guided (by rig data) mode");

		gd.addCheckbox    ("Correct disparity for infinity tiles )has to disable until code fixed)",                 this.ly_inf_disp);
		gd.addCheckbox    ("Force convergence correction during extrinsic, even with no infinity data",              this.ly_inf_force);
		gd.addCheckbox    ("*Use polynomial correction, false - correct tilt/azimuth/roll of each sensor)",          this.ly_poly);

		gd.addMessage     ("--- Lazy eye parameters ---");
		gd.addNumericField("Sample size (side of a square)",                                                    this.ly_smpl_side,  0);
		gd.addNumericField("Number after removing worst (should be >1)",                                        this.ly_smpl_num,  0);
		gd.addMessage     ("Maximal measured relative disparity = "+ (0.8*disp_scan_step)+" (0.8 * disp_scan_step)");
		//			gd.addNumericField("Maximal measured relative disparity",                                               this.ly_meas_disp,  3);
		gd.addNumericField("Maximal RMS of the remaining tiles in a sample",                                    this.ly_smpl_rms,  5);
		gd.addNumericField("Maximal full disparity difference to 8 neighbors",                                  this.ly_disp_var, 8,5,"pix",
				"Full allowed mismatch is a sum of absolute and disparity times relative");
		gd.addNumericField("Maximal relative full disparity difference to 8 neighbors",                         this.ly_disp_rvar, 8,5,"",
				"Full allowed mismatch is a sum of absolute and disparity times relative");
		gd.addNumericField("Maximal full disparity difference to 8 neighbors with GT",                          this.ly_disp_var_gt, 8,5,"pix",
				"Full allowed mismatch is a sum of absolute and disparity times relative (relaxed when ground truth is available)");
		gd.addNumericField("Maximal relative full disparity difference to 8 neighbors with GT",                 this.ly_disp_rvar_gt, 8,5,"",
				"Full allowed mismatch is a sum of absolute and disparity times relative (relaxed when ground truth is available)");
		gd.addNumericField("Reduce weight of higher disparity tiles",                                           this.ly_norm_disp, 5);
		gd.addMessage     ("--- Lazy eye multi-step fitting ---");
		gd.addNumericField("Any (near) saturated pixels - discard tile (see sat_level also)",                   this.lym_overexp,  10);
		gd.addCheckbox    ("Update target disparity after each step",                                           this.lym_update_disp);
		gd.addNumericField("Maximal number of iterations",                                                      this.lym_iter,  0);
		gd.addNumericField("Parameter vector difference to exit (main camera)",                                 this.lym_change,  10);
		gd.addNumericField("Parameter vector difference to exit (aux camera)",                                  this.lym_change_aux,  10);

		gd.addNumericField("Parameter vector difference to exit from polynomial correction",                    this.lym_poly_change,  10);

		gd.addMessage     ("--- Lazy eye samples filter ---");
		gd.addCheckbox    ("Filter lazy eye pairs by their values",                                             this.lyf_filter);
		gd.addNumericField("Fileter sample side (if 8, 8 x8 masked, 16x16 sampled)",                            this.lyf_smpl_side,  0);
		gd.addNumericField("Maximal RMS (all components to components average)",                                this.lyf_rms_max,  3);
		gd.addNumericField("Keep best fit samples, discard worst",                                              this.lyf_frac_keep,  3);
		gd.addNumericField("Minimal number of tiles remaining in the sample",                                   this.lyf_min_samples,  0);
		gd.addCheckbox    ("Replace samples with a single average with equal weight",                           this.lyf_norm_center);
		gd.addNumericField("Scale calculated correction vector",                                                this.ly_corr_scale,  3);
		gd.addMessage     ("--- Parameters specific to LY adjustments with dual-camera rig data (as ground truth ---");

		gd.addCheckbox    ("Use samples filter with rig data",                                                  this.lyr_filter_ds,
				"Raw measured data may not need filtering when ground truth disparity is available");
		gd.addCheckbox    ("Filter lazy eye pairs by their values wen GT data is available",                    this.lyr_filter_lyf,
				"Same as \"Filter lazy eye pairs by their values\" above, but for the rig-guided adjustments");


		//  			gd.addNumericField("Use square this size side to detect outliers",                                      this.fcorr_sample_size,  0);
		//  			gd.addNumericField("Keep tiles only if there are more in each square",                                  this.fcorr_mintiles,     0);
		//  			gd.addNumericField("Remove this fraction of tiles from each sample",                                    this.fcorr_reloutliers,  3);
		//  			gd.addNumericField("Gaussian blur channel mismatch data",                                               this.fcorr_sigma,        3);

		///  			gd.addNumericField("Calculated from correlation offset vs. actual one (not yet understood)",            this.corr_magic_scale,  3);

		gd.addTab         ("3D", "3D reconstruction");
		gd.addMessage     ("--- 3D reconstruction ---");
		gd.addCheckbox    ("Show generated textures",                                                                this.show_textures);
		gd.addCheckbox    ("show intermediate results of filtering",                                                 this.debug_filters);

		gd.addNumericField("Minimal noise-normalized pixel difference in a channel to suspect something",            this.min_smth,  3);
		gd.addNumericField("Reliable noise-normalized pixel difference in a channel to have something ",             this.sure_smth,  3);
		gd.addNumericField("Disparity range to be considered background",                                            this.bgnd_range,  3);
		gd.addNumericField("Disparity difference from the center (provided) disparity to trust",                     this.other_range,  3);

		gd.addNumericField("Minimal 4-corr strength to trust tile",                                                  this.ex_strength,  3);
		gd.addNumericField("Minimal 4-corr strength divided by channel diff for new (border) tiles",                 this.ex_nstrength,  3);

		gd.addCheckbox    ("Allow expansion over previously identified background (infinity)",                       this.ex_over_bgnd);
		gd.addNumericField("When expanding over background, disregard lower disparity ",                             this.ex_min_over,  3);

		gd.addTab         ("Plates", "Plates filtering when building initial z-map");
		gd.addMessage     ("********* Plates filtering when building initial z-map *********");
		gd.addNumericField("If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds",           this.pt_super_trust,  3);
		gd.addCheckbox    ("Do not replace raw tiles by the plates, if raw is closer (like poles)",                  this.pt_keep_raw_fg);
		gd.addNumericField("Scale plates strength before comparing to raw strength",                                 this.pt_scale_pre,  3);
		gd.addNumericField("Scale plates strength when replacing raw (plates d/s data is more reliable if it exists)",           this.pt_scale_post,  3);

		gd.addNumericField("Minimal strength to be considered definitely background",                                this.bgnd_sure,  3);
		gd.addNumericField("Maximal strength to ignore as non-background",                                           this.bgnd_maybe,  3);

		gd.addNumericField("Number of tiles in a cluster to seed (just background?)",                                this.min_clstr_seed,   0);
		gd.addNumericField("Number of tiles in a cluster not close to other clusters (more than 2 tiles apart)",           this.min_clstr_lone,   0);
		gd.addNumericField("Minimal total strength of the cluster",                                                  this.min_clstr_weight,  3);
		gd.addNumericField("Minimal maximal strength of the cluster",                                                this.min_clstr_max,  3);

		gd.addNumericField("Fill gaps betsween clusters, see comments for 'grow'",                                   this.fill_gaps,   0);
		gd.addNumericField("Same as fill_gaps above, on the final pass",                                             this.fill_final,   0);
		gd.addNumericField("Number of tiles in a cluster to block (just non-background?)",                           this.min_clstr_block,   0);
		gd.addNumericField("Number of tiles to grow tile selection (1 - hor/vert, 2 - hor/vert/diagonal)",           this.bgnd_grow,   0);

		gd.addCheckbox    ("Use old ortho features processing (ortho_* parameters, false - use or_*)",               this.ortho_old);
		gd.addMessage     ("--- old ones, new are in  \"Ortho+4\" tab---");
		gd.addNumericField("Minimal strength of hor correlation to be used instead of full 4-pair correlation",      this.ortho_min_hor,  3);
		gd.addNumericField("Minimal strength of vert correlation to be used instead of full 4-pair correlation",     this.ortho_min_vert,  3);
		gd.addNumericField("Vert/hor (or hor/vert) strength to be used instead of the full correlation",             this.ortho_asym,  3);
		gd.addNumericField("Vert/hor (or hor/vert) strength exceeding scaled 4-pair strength",                       this.ortho_over4,  3);

		gd.addNumericField("Minimal strength of hor/vert to bridge over",                                            this.ortho_sustain,  3);
		gd.addNumericField("minimal run of hor/vert tiles to be considered (at least from one side)",                this.ortho_run,      0);
		gd.addNumericField("Minimal maximal strength in an ortho run (max. in run >=...)",                           this.ortho_minmax,  3);
		gd.addNumericField("Number of tiles to bridge over hor/vert gaps",                                           this.ortho_bridge,   0);
		gd.addNumericField("Maximal disparity RMS in a run to replace by average)",                                  this.ortho_rms,  3);
		gd.addNumericField("Convolve hor/vert strength by 3*(2*l+1) kernels to detect multi-tile features",          this.ortho_half_length,   0);
		gd.addNumericField("Fraction of convolved ortho in a mix with raw",                                          this.ortho_mix,  3);

		gd.addTab         ("Ortho+4", "Combining of ortho and 4-pair correlations");
		gd.addMessage     ("--- Combining of ortho and 4-pair correlations ---");
		gd.addCheckbox    ("Apply ortho correction to horizontal correlation (vertical features)",                      this.or_hor);
		gd.addCheckbox    ("Apply ortho correction to vertical correlation (horizontal features)",                      this.or_vert);
		gd.addNumericField("Blur sigma: verically for horizontal correlation, horizontally - for vertically",           this.or_sigma,  3);
		gd.addNumericField("3-point sharpening (-k, +2k+1, -k)",                                                        this.or_sharp,  3);
		gd.addNumericField("Scale ortho correletion strength relative to 4-directional one",                            this.or_scale,  3);
		gd.addNumericField("Subtract from scaled correlation strength, limit by 0",                                     this.or_offset,  3);
		gd.addNumericField("Minimal ratio of orthogonal strengths required for disparity replacement",                  this.or_asym,  3);
		gd.addNumericField("Minimal scaled offset ortho strength to normal strength needed for replacement",            this.or_threshold,  3);
		gd.addNumericField("Minimal horizontal absolute scaled offset ortho strength needed for replacement",           this.or_absHor,  3);
		gd.addNumericField("Minimal vertical absolute scaled offset ortho strength needed for replacement",             this.or_absVert,  3);
		gd.addNumericField("Maximal disparity to apply ortho correction",                                               this.or_maxDisp,  3,8,"pix",
				"Maximal disparity to apply ortho correction");

		gd.addTab         ("Vert", "Fix vertical structures, such as street poles");
		gd.addMessage     ("--- Fix vertical structures, such as street poles ---");
		gd.addCheckbox    ("Continue vertical structures to the ground",                                                this.poles_fix);
		gd.addNumericField("Number of tiles to extend over the poles bottoms",                                          this.poles_len,   0);
		gd.addNumericField("Maximal ratio of invisible to visible pole length",                                         this.poles_ratio,  3);
		gd.addNumericField("Set new pole segment strength to max of horizontal correlation and this value",             this.poles_min_strength,  3);
		gd.addCheckbox    ("Set disparity to that of the bottom of existing segment (false - use hor. disparity)",          this.poles_force_disp);

		gd.addNumericField("Maximal number of output meshes to generate",                                            this.max_clusters,   0);
		gd.addCheckbox    ("Remove all unneeded scans when generating x3d output to save memory",                    this.remove_scans);
		gd.addCheckbox    ("Generate x3d output",                                                                    this.output_x3d);
		gd.addCheckbox    ("Generate Wavefront obj output",                                                          this.output_obj);
		gd.addCheckbox    ("Correct lens geometric distortions in a model (will need backdrop to be corrected too)",           this.correct_distortions);
		gd.addCheckbox    ("Show generated triangles",                                                               this.show_triangles);
		gd.addCheckbox    ("Weight-average disparity for the whole cluster ",                                        this.avg_cluster_disp);
		gd.addNumericField("Maximal disparity difference in a triangle face to show",                                this.maxDispTriangle,  6);
		gd.addNumericField("Distance to generate backdrop (0 - use regular backdrop)",                               this.infinityDistance,  8);
		gd.addNumericField(" Minimal number of background tiles to generate background",                             this.min_bgnd_tiles,   0);

		gd.addCheckbox    ("Split into shells with flaps",                                                           this.shUseFlaps);
		gd.addCheckbox    ("Aggressive fade alpha (whole boundary)",                                                 this.shAggrFade);
		gd.addNumericField("Minimal shell area (not counting flaps",                                                 this.shMinArea,   0);
		gd.addNumericField("Minimal value of the shell maximum strength",                                            this.shMinStrength,  6);

		gd.addTab         ("Thin ice", "Thin ice parameters (obsolete)");
		gd.addMessage     ("--- Thin ice parameters (obsolete) ---");
		gd.addNumericField("Relative disparity rigidity in vertical direction",                                      this.tiRigidVertical,  6);
		gd.addNumericField("Relative disparity rigidity in horizontal direction",                                    this.tiRigidHorizontal,  6);
		gd.addNumericField("Relative disparity rigidity in diagonal   direction",                                    this.tiRigidDiagonal,  6);
		gd.addNumericField("Strength floor - subtract (limit with 0) before applying",                               this.tiStrengthOffset,  6);
		gd.addNumericField("Divide actual disparity by this  before Math.pow and applying to TI",                    this.tiDispScale,  6);
		gd.addNumericField("Apply pow to disparity (restore sign) for disparity difference pressure on TI",          this.tiDispPow,  6);
		gd.addNumericField("tiDispPull: multiply strength*disparity difference to pull force",                       this.tiDispPull,  6);
		gd.addNumericField("Scale tiDispPull for pre-final pass",                                                    this.tiDispPullPreFinal,  6);
		gd.addNumericField("Scale tiDispPull for final pass",                                                        this.tiDispPullFinal,  6);
		gd.addNumericField("Normalize stresses to average disparity if it is above threshold",                       this.tiBreakNorm,  6);
		gd.addNumericField("TI break value of abs(d0-3d1+3d2-d3)",                                                   this.tiBreak3,  6);
		gd.addNumericField("TI break value of (d0-3d1+3d2-d3) * (d1 - d2)",                                          this.tiBreak31,  6);
		gd.addNumericField("TI break value of (-d0+d1+d2-d3) * abs(d1 - d2)",                                        this.tiBreak21,  6);
		gd.addNumericField("TI disparity threshold to remove as too far tiles",                                      this.tiBreakFar,  6);
		gd.addNumericField("TI disparity threshold to remove as too near tiles",                                     this.tiBreakNear,  6);
		gd.addNumericField("TI break mode: +1: abs(3-rd derivative), +2: -(3-rd * 1-st), +4: -(2-nd * abs(1-st))",           this.tiBreakMode,  0);
		gd.addNumericField("Amplify colinear breaks in neighbor tiles",                                              this.tiBreakSame,  6);
		gd.addNumericField("Amplify 90-degree turnintg breaks in neighbor tiles",                                    this.tiBreakTurn,  6);

		gd.addNumericField("Heal disparity gap before pre-last smooth",                                              this.tiHealPreLast,  6);
		gd.addNumericField("Heal disparity gap before last smooth",                                                  this.tiHealLast,  6);
		gd.addNumericField("Maximal length of an internal break in the cluster to heal",                             this.tiHealSame, 0);

		gd.addNumericField("Maximal number of TI iterations for each step",                                          this.tiIterations, 0);
		gd.addNumericField("Iteration maximal error (1/power of 10)",                                                this.tiPrecision,  0);
		gd.addNumericField("Number of cycles break-smooth (after the first smooth)",                                 this.tiNumCycles,  0);

		gd.addTab         ("Fg/Bg", "Fg/Bg separation");
		gd.addMessage     ("--- Fg/Bg separation ---");
		gd.addCheckbox    ("Apply super-tiles during refine passes",                                                 this.stUseRefine);
		gd.addCheckbox    ("Apply super-tiles during pass2 ",                                                        this.stUsePass2);
		gd.addCheckbox    ("Apply super-tiles during render",                                                        this.stUseRender);

		gd.addCheckbox    ("Show supertiles histograms",                                                             this.stShow);
		gd.addNumericField("Super tile size (square, in tiles)",                                                     this.stSize,  0);
		gd.addNumericField("Disparity histogram step for far objects",                                               this.stStepFar,  6);
		gd.addNumericField("Disparity histogram step for near objects",                                              this.stStepNear,  6);
		gd.addNumericField("Disparity threshold to switch from linear to logarithmic steps",                         this.stStepThreshold,  6);
		gd.addNumericField("Minimal disparity (center of a bin)",                                                    this.stMinDisparity,  6);

		gd.addTab         ("Supertiles", "Build supertiles (plates)");
		//  			gd.addNumericField("Maximal disparity (center of a bin)",                                                    this.stMaxDisparity,  6);
		gd.addMessage     ("Maximal disparity (center of a bin) - using grow_disp_max="+          this.grow_disp_max);
		//  			gd.addNumericField("Subtract from strength, discard negative",                                               this.stFloor,  6);
		//  			gd.addNumericField("Raise strength to this power ",                                                          this.stPow,  6);
		gd.addNumericField("Blur disparity histogram (sigma in bins)",                                               this.stSigma,  6);
		gd.addNumericField("Minimal backgroubnd disparity to extract as a maximum from the supertiles",              this.stMinBgDisparity,  6);
		gd.addNumericField("Minimal fraction of the disparity histogram to use as background",                       this.stMinBgFract,  6);
		gd.addNumericField("Use background disparity from supertiles if tile strength is less",                      this.stUseDisp,  6);
		gd.addNumericField("Multiply st strength if used instead of regular strength ",                              this.stStrengthScale,  6);

		//  			gd.addMessage     ("Sample filtering mode");

		gd.addMessage     ("Other supertile parameters");
		gd.addNumericField("Grow initial selection before processing supertiles, odd - ortho. <0 - use all tiles",          this.stGrowSel,  0);
		gd.addNumericField("Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert",          this.stMeasSel,  0);
		gd.addNumericField("Consider merging initial planes if disparity difference below",                          this.stSmallDiff,  6);
		gd.addNumericField("Consider merging initial planes if jumps between ratio above",                           this.stHighMix,  6);

		gd.addNumericField("Outlier tiles weaker than this may be replaced from neighbors",                          this.outlierStrength,  6);
		gd.addNumericField("Replace weak outlier tiles that do not have neighbors within this disparity difference",           this.outlierDiff,  6);
		gd.addNumericField("Replace weak outlier tiles that have higher disparity than weighted average",            this.outlierDiffPos,  6);
		gd.addNumericField("Replace weak outlier tiles that have lower disparity than weighted average",             this.outlierDiffNeg,  6);

		gd.addCheckbox    ("Combine with all previous after refine pass",                                            this.combine_refine);
		gd.addNumericField("Disregard weaker tiles when combining scans",                                            this.combine_min_strength,  6);
		gd.addNumericField("Disregard weaker tiles when combining scans  for horizontal correlation",                this.combine_min_hor,  6);
		gd.addNumericField("Disregard weaker tiles when combining scans  for vertical correlation",                  this.combine_min_vert,  6);
		//  			gd.addNumericField("Do not re-measure correlation if target disparity differs from some previous by this",          this.unique_tolerance,  6);

		gd.addTab         ("Tile Filter", "Filter for measured tile data");
		gd.addCheckbox    ("Use sample mode (false - regular tile mode)",                                            this.stSmplMode);
		gd.addNumericField("Sure strength",                                                                          this.mlfp.strength_sure,  4,6,"",
				"Do not filter tile disparities above this strength level (on top of strength floor");
		gd.addNumericField("Strength floor",                                                                         this.mlfp.strength_floor,  4,6,"",
				"Subtract from strength, discard negative");
		gd.addNumericField("Raise strength to this power ",                                                          this.mlfp.strength_pow,  6);
		gd.addNumericField("Sample size (side of a square)",                                                         this.mlfp.smplSide,  0);
		gd.addNumericField("Number after removing worst",                                                            this.mlfp.smplNum,  0);
		gd.addNumericField("Maximal RMS of the remaining tiles in a sample",                                         this.mlfp.smplRms,  6);
		gd.addCheckbox    ("Use window function for the samples",                                                    this.mlfp.smplWnd);

		gd.addMessage     ("Additional parameters to prevent increase of the  distance to far (small) objects",
				"These parameters are also used when expanding tiles from the known ones");
		gd.addNumericField("Maximal absolute plane \"tilt\" in pixels/tile",                                         this.mlfp.max_abs_tilt,  4, 6,
				"pix/tile", "difference between disparity of the neighbor tiles in disparity space image");
		gd.addNumericField("Maximal relative plane \"tilt\"",                                                        this.mlfp.max_rel_tilt,  4, 6,
				"1/tile", "pixel/tile tilte divided by the full disparity of the tle");

		gd.addNumericField("Tilt cost for damping insufficient plane data",                                          this.mlfp.damp_tilt,  6, 8,
				"", "Allows to find plane when there are not enough tiles to process");

		gd.addNumericField("Disparity switch between filtering modes",                                               this.mlfp.min_tilt_disp,  4,6,
				"pix","Objects that are closer (larger disparity) use tilted plane model, far objects use maximal among neighbors disparity");
		gd.addNumericField("Mode transition range (between tilted and maximal disparity)",                           this.mlfp.transition,  4,6,
				"pix","Disparity range to gradually switch between maximal and tilted modes");
		gd.addNumericField("Far objects filtering mode (0 - off, 1,2 - power of disparity)",                         this. mlfp.far_mode,  0,3,
				"","So far only 2 far modes: 1 - remove outliers from both sides, 2 - only from far side, 0 - use tilted mode for all tiles");
		gd.addNumericField("Raise disparity to this power before averaging for far objects",                         this.mlfp.far_power,  4,6,
				"", "Parameter to specify preference of the larger disparity when averaging");


		gd.addTab         ("Grow", "Growing disparity range to scan");
		gd.addMessage     ("========= Growing disparity range to scan ========");
		gd.addNumericField("Try these number of tiles around known ones",                                                   this.grow_sweep,  0);
		gd.addNumericField("Maximal disparity to try",                                                                      this.grow_disp_max,  6);
		gd.addNumericField("Trust measured disparity within +/- this value",                                                this.grow_disp_trust,  6);
		gd.addNumericField("Increase disparity (from maximal tried) if nothing found in that tile",                         this.grow_disp_step,  6);
		gd.addNumericField("Grow more only if at least one channel has higher variance from others for the tile",           this.grow_min_diff,  6);
		gd.addCheckbox    ("Retry tiles around known foreground that have low max_tried_disparity",                         this.grow_retry_far);
		gd.addCheckbox    ("Scan full range between max_tried_disparity of the background and known foreground",            this.grow_pedantic);
		gd.addCheckbox    ("Retry border tiles that were identified as infinity earlier",                                   this.grow_retry_inf);

		gd.addTab         ("Filter", "\"plate\" filtering when growing parameters");
		gd.addMessage     ("--- more growing parameters ---");
		gd.addCheckbox    ("New expansion mode",                                                                            this.gr_new_expand);
		//  			gd.addNumericField("Expansion steps limit",                                                                         this.gr_max_expand,  0);
		gd.addNumericField("Strength floor for multi-tile (now 5x5) samples (normally < combine_min_strength) ",            this.fds_str_floor,  6);
		gd.addNumericField("Over background extra reliable strength",                                                       this.gr_ovrbg_cmb,  6);
		gd.addNumericField("Over background extra reliable strength horizontal",                                            this.gr_ovrbg_cmb_hor,  6);
		gd.addNumericField("Over background extra reliable strength vertical",                                              this.gr_ovrbg_cmb_vert,  6);
		gd.addNumericField("Over background filtered extra reliable strength",                                              this.gr_ovrbg_filtered,  6);

		gd.addMessage     ("--- \"plate\" filtering when growing parameters ---");
		gd.addNumericField("Strength power exponent for tilted plates growing",                                             this.fds_str_pow,  6);
		gd.addNumericField("Sample size (side of a square) for tilted plates growing",                                      this.fds_smpl_side,  0);
		gd.addNumericField("Number of tiles in a square tilted plate (should be >1)",                                       this.fds_smpl_num,  0);
		gd.addNumericField("Maximal RMS for the tiles to the tilted plate",                                                 this.fds_smpl_rms,  6);
		gd.addNumericField("Maximal relative RMS for the tiles to the tilted plate - multiply by disparity and add",           this.fds_smpl_rel_rms,  6);
		gd.addCheckbox    ("Use window function for the square sample plates",                                              this.fds_smpl_wnd);
		gd.addNumericField("Maximal growing plate tilt in disparity pix per tile",                                          this.fds_abs_tilt,  6);
		gd.addNumericField("Maximal relative growing plate tilt in disparity pix per tile per disparity pixel",             this.fds_rel_tilt,  6);

		gd.addTab         ("Periodic", "Detect and filter periodic features");
		gd.addCheckbox    ("Apply filter to remove false correlations of periodic structures",                              this.per_filter,
				"Detect repetitive features (such as windows grid) and remove false matches");
		gd.addNumericField("Maximal residual disparity difference to use for periodic features",                            this.per_trustedCorrelation,  4, 6,"pix",
				"do not trust correlation results with larger absolute value of residual disparity");
		gd.addNumericField("Initial half-width of the disparity range per cluster",                                         this.per_initial_diff,  4, 6,"pix",
				"Select tiles this far from the local maximums");
		gd.addNumericField("Strength floor to detect periodic features",                                                    this.per_strength_floor,  4, 6,"",
				"Only use stronger correlations, subtract this value for weighted averages");
		gd.addNumericField("Maximums seeds should have at least this much over strength floor",                             this.per_strength_max_over,  4, 6,"",
				"When growing clusters, only start with the correlations this much stronger than strength floor ");
		gd.addNumericField("Minimal feature period",                                                                        this.per_min_period,  4, 6,"",
				"Only detect/filter feature that have period above this value");
		gd.addNumericField("Minimal number of periods in periodic structures",                                              this.per_min_num_periods,  0, 2,"",
				"Minimal number of full periods to be detected");
		gd.addNumericField("Maximal period variations",                                                                     this.per_disp_tolerance,  4, 6,"pix",
				"The detected periods should be higher than this value");
		gd.addNumericField("Disparity difference to match neighbors",                                                       this.per_disp_match,  0, 2,"",
				"Neighbor tiles should have smaller disparity difference from gthe center to qualify");
		gd.addNumericField("Extra strength to treat match as strong (for hysteresis)",                                      this.per_strong_match_inc,  0, 2,"",
				"The layer that does not have any match in the direction where some other layer has a strong match will be removed");


		gd.addTab         ("Macro", "Macro tiles correlation parameters");
		gd.addMessage     ("--- Macro correlation parameters ---");
		gd.addNumericField("Macro disparity scan step (actual disparity step is 8x)",                                       this.mc_disp8_step,  6);
		gd.addNumericField("Trust measured macro(8x)  disparity within +/- this value ",                                    this.mc_disp8_trust,  6);
		gd.addNumericField("Minimal macro correlation strength to process",                                                 this.mc_strength,  6);
		gd.addNumericField("Do not re-measure macro correlation if target disparity differs less",                          this.mc_unique_tol,  6);

		gd.addNumericField("When consolidating macro results, exclude high residual disparity",                             this.mc_trust_fin,  6);
		gd.addNumericField("Gaussian sigma to reduce weight of large residual disparity",                                   this.mc_trust_sigma,  6);
		gd.addNumericField("Weight from ortho neighbor supertiles",                                                         this.mc_ortho_weight,  6);
		gd.addNumericField("Weight from diagonal neighbor supertiles",                                                      this.mc_diag_weight,  6);
		gd.addNumericField("Do not remove measurements farther from the kept ones",                                         this.mc_gap,  6);


		gd.addNumericField("weight of variance data ",                                                                      this.mc_weight_var,  4, 6,"",
				"This is what was used before - can detect tin wires, probably");
		gd.addNumericField("Weight of average intensity",                                                                   this.mc_weight_Y,  4, 6,"",
				"Combination of R,G,B fro the tile");
		gd.addNumericField("Weight of average color difference",                                                            this.mc_weight_RBmG,  4, 6,"",
				"0.5*(R+B)-G");


		gd.addTab         ("Grow2", "More disparity range growing parameters");
		gd.addMessage     ("--- more growing parameters ---");
		gd.addNumericField("Discard variant if it requests too few tiles",                                                  this.gr_min_new,  0);
		gd.addCheckbox    ("Expand only unambiguous tiles over previously undefined",                                       this.gr_var_new_sngl);
		gd.addCheckbox    ("Expand unambiguous and FOREGROUND tiles over previously UNDEFINED",                             this.gr_var_new_fg);
		gd.addCheckbox    ("Expand unambiguous and FOREGROUND tiles over already DEFINED",                                  this.gr_var_all_fg);
		gd.addCheckbox    ("Expand unambiguous and BACKGROUND tiles over previously UNDEFINED",                             this.gr_var_new_bg);
		gd.addCheckbox    ("Expand unambiguous and BACKGROUND tiles over already DEFINED",                                  this.gr_var_all_bg);
		gd.addCheckbox    ("Try next disparity range that was not tried before",                                            this.gr_var_next);
		gd.addNumericField("How far to extend over previously undefined disparity tiles",                                   this.gr_num_steps,  0);
		gd.addNumericField("How far to extend over previously determined disparity tiles",                                  this.gr_steps_over,  0);
		gd.addNumericField("Extend sample plate square side",                                                               this.gr_smpl_size,  0);
		gd.addNumericField("Extend at least this number of the seed tiles",                                                 this.gr_min_pnts,  0);
		gd.addCheckbox    ("Use window function for square sample ",                                                        this.gr_use_wnd);
		gd.addNumericField("Tilt cost for damping insufficient plane data",                                                 this.gr_tilt_damp,  6);
		gd.addNumericField("When growing, range of disparities to be extended without far/near subdivision",                this.gr_split_rng,  6);
		gd.addNumericField("Consider far/near tiles within that range from the farthest/closest",                           this.gr_same_rng,  6);
		gd.addNumericField("Maximal difference from the old value when smoothing",                                          this.gr_diff_cont,  6);
		gd.addNumericField("Maximal filter disparity absolute tilt (pix per tile)",                                         this.gr_abs_tilt,  6);
		gd.addNumericField("Maximal filter disparity tilt (pix / disparity) per tile",                                      this.gr_rel_tilt,  6);
		gd.addNumericField("Maximal number of smoothing steps (reduce if long?)",                                           this.gr_smooth,  0);
		gd.addNumericField("Maximal change to finish smoothing iterations",                                                 this.gr_fin_diff,  6);
		gd.addNumericField("Do not re-measure correlation if target disparity differs from some previous less",             this.gr_unique_tol,  6);
		gd.addNumericField("Larger tolerance for expanding (not refining)",                                                 this.gr_unique_pretol,  6);

		gd.addTab         ("Alt CLusterize", "Alternative initial tiles clusterization");

		gd.addCheckbox    ("Modify cluster strengths",                                                                      this.ft_mod_strength,
				"Supplement sum of strengths with other parameters, such as density and height ");
		gd.addCheckbox    ("Enable alternative initial tile clusterization",                                                this.ft_clusterize_by_highest,
				"Clusterize using disparity horizontal maximums for fronto planes and minimums - for horizontal. False - use histograms");
		gd.addNumericField("Disparity blur sigma",                                                                          this.ft_clust_sigma, 4, 6,"pix",
				"Blur disparity before finding each line max (for fronto planes ) or min (for horizontal planes) during initial clusterization");
		gd.addNumericField("Absolute disparity range for fronto clusters",                                                  this.ft_disp_arange_vert, 4, 6,"pix",
				"Disparity range for blurred disparity (down from max disparity) for fronto planes");
		gd.addNumericField("Relative disparity range for fronto clusters",                                                  this.ft_disp_rrange_vert, 4, 6,"pix/pix",
				"Increase disparity range for fronto clusters for each disparity pixel");
		gd.addNumericField("Absolute disparity range for horizontal clusters",                                              this.ft_disp_arange_hor, 4, 6,"pix",
				"Disparity range for blurred disparity (up from min disparity to horizontal difference) for horizontal planes");
		gd.addNumericField("Relative disparity range for horizontal clusters",                                              this.ft_disp_rrange_hor, 4, 6,"pix/pix",
				"Increase disparity range for horizontal clusters for each disparity pixel");
		gd.addNumericField("Actual disparity positive tolerance over blurred disparity max range",                          this.ft_tolerance_above_near, 4, 6,"pix",
				"Allow measured tile disparity above cluster disparity range for fronto clusters");
		gd.addNumericField("Actual disparity negative tolerance over blurred disparity max range",                          this.ft_tolerance_below_near, 4, 6,"pix",
				"Allow measured tile disparity below cluster disparity range for fronto planes");
		gd.addNumericField("Actual disparity positive tolerance over blurred disparity min range",                          this.ft_tolerance_above_far, 4, 6,"pix",
				"Allow measured tile disparity above cluster disparity range for horizontal planes");
		gd.addNumericField("Actual disparity negative tolerance over blurred disparity min range",                          this.ft_tolerance_below_far, 4, 6,"pix",
				"Allow measured tile disparity below cluster disparity range for horizontal planes");
		gd.addNumericField("Fronto/horizontal selections overlap",                                                          this.ft_hor_vert_overlap,  0,6,"",
				"Allow clusters tile sharing between fronto and horizontal. 2 - 1 tile in 8 directions, 1 - 1 tile in 4 directions");
		gd.addNumericField("Mark as used if has used/disabled neighbors ",                                                  this.ft_used_companions,  0,6,"",
				"Cell that has this many new used companions is considered used (borders and already use3d are considered used too)");
		gd.addNumericField("Minimal number of new used cells among new/old used and marginal tiles",                        this.ft_used_true_companions,  0,6,"",
				"There should be at least this many new selected tiles among neighbors");


		gd.addTab         ("Plane Det", "Planes detection");
		gd.addMessage     ("--- Planes detection ---");
		gd.addCheckbox    ("Always start with disparity-most axis (false - lowest eigenvalue)",                      this.plPreferDisparity);
		gd.addNumericField("Normalize disparities to the average if above",                                          this.plDispNorm,  4,6, "pix");
		gd.addNumericField("Fronto tolerance",                                                                       this.plFrontoTol,  4,6,"pix",
				"Fronto tolerance (pix) - treat almost fronto planes as fronto (constant disparity). If <= 0 - disable this feature");
		gd.addNumericField("Fronto RMS",                                                                             this.plFrontoRms,  4,6,"pix",
				"Target half-thikness of the fronto planes. Similar to  sqrt(plMaxEigen) for other planes");
		gd.addNumericField("Fronto offset",                                                                          this.plFrontoOffs,  4,6,"pix",
				"Increasing weights of the near tiles by using difference between tile disparity and reduced by this value average as weight. If <= 0 - disable feature");
		gd.addNumericField("Fronto power",                                                                           this.PlFrontoPow,  4,6,"pix",
				"Increasing weights of the near tiles by even more (see previous parameter) by raising disparity difference to this power");
		gd.addNumericField("Blur disparity histograms for constant disparity clusters by this sigma (in bins)",             this.plBlurBinVert,  6);
		gd.addNumericField("Blur disparity histograms for horizontal clusters by this sigma (in bins)",                     this.plBlurBinHor,  6);
		gd.addNumericField("Maximal normalized disparity difference when initially assigning to vertical plane",            this.plMaxDiffVert,  6);
		gd.addNumericField("Maximal normalized disparity difference when initially assigning to horizontal plane",          this.plMaxDiffHor,  6);
		gd.addNumericField("Number of initial passes to assign tiles to vert (const disparity) and hor planes",             this.plInitPasses,  0);

		gd.addNumericField("Minimal number of points for plane detection",                                           this.plMinPoints,  0);
		gd.addNumericField("Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below",           this.plTargetEigen,  6);
		gd.addNumericField("Maximal fraction of outliers to remove",                                                 this.plFractOutliers,  6);
		gd.addNumericField("Maximal number of outliers to remove",                                                   this.plMaxOutliers,  0);
		gd.addNumericField("Minimal total strength of a plane",                                                      this.plMinStrength,  6);
		gd.addNumericField("Maximal eigenvalue of a plane",                                                          this.plMaxEigen,  6);
		gd.addNumericField("Add to eigenvalues of each participating plane and result to validate connections",          this.plEigenFloor,  6);
		gd.addNumericField("Consider plane to be a \"stick\" if second eigenvalue is below",                         this.plEigenStick,  6);
		gd.addNumericField("Not a plate if sin^2 between normals from disparity and world exceeds this",             this.plBadPlate,  6);
		gd.addCheckbox    ("Combine 'other' plane with the current (unused)",                                        this.plDbgMerge);
		gd.addNumericField("Worst case worsening after merge",                                                       this.plWorstWorsening,  6);
		gd.addNumericField("Worst case worsening for thin planes",                                                   this.plWorstWorsening2,  6);
		gd.addNumericField("Worst case worsening after merge with equal weights",                                    this.plWorstEq,  6);
		gd.addNumericField("Worst case worsening for thin planes with equal weights",                                this.plWorstEq2,  6);
		gd.addNumericField("If result of the merged planes is below, OK to use thin planes (higher) threshold ",          this.plOKMergeEigen,  6);
		gd.addNumericField("Maximal sine squared of the world angle between planes to merge. Set to >= 1.0 to disable",           this.plMaxWorldSin2,  6);
		gd.addNumericField("Do not compare sin() between planes, if at least one has too small axis ratio",          this.pl2dForSin,  6);
		gd.addNumericField("Relax merge requirements for weaker planes",                                             this.plWeakWorsening,  6);
		gd.addNumericField("Maximal overlap between the same supertile planes to merge",                             this.plMaxOverlap,  6);

		gd.addTab         ("Merge", "Merge same supetile planes if at least one is weak and they do not differ much");
		gd.addMessage     ("--- Merge same supetile planes if at least one is weak and they do not differ much ---");
		gd.addNumericField("Maximal weight of the weak plane to merge (first variant)",                              this.plWeakWeight,  6);
		gd.addNumericField("Maximal eigenvalue of the result of non-weighted merge (first variant)",                 this.plWeakEigen,  6);
		gd.addNumericField("Maximal weight of the weak plane to merge (second variant)",                             this.plWeakWeight2,  6);
		gd.addNumericField("Maximal eigenvalue of the result of non-weighted merge (second variant)",                this.plWeakEigen2,  6);
		gd.addNumericField("Do not merge if any sqrt of merged eigenvalue exceeds scaled sum of components",           this.plSumThick,  6);
		gd.addNumericField("When calculating non-exclusive planes, do not use neighbors with high cost",             this.plNeNeibCost,  6);
		gd.addNumericField("When calculating non-exclusive planes, use cenrter plane relative weight",               this.plNeOwn,  6);
		gd.addNumericField("When calculating exclusive planes links, do not use neighbors with high cost",           this.plExNeibCost,  6);
		gd.addNumericField("Scale down maximal costs for smoothed planes links (tighter requirements)",              this.plExNeibSmooth,  6);
		gd.addNumericField("Cost threshold for merging same tile planes if the plane has connected neighbors",           this.plMergeCostStar,  6);
		gd.addNumericField("Cost threshold for merging same tile planes if not connected",                           this.plMergeCost,  6);

		gd.addMessage     ("--- Merging planes with topological conflicts ---");
		gd.addCheckbox    ("Try to merge conflicting planes",                                                               this.plConflMerge);
		gd.addNumericField("Scale parameters to relax planes fit for merging conflicting planes",                           this.plConflRelax,  6);
		gd.addCheckbox    ("Only merge conflicting planes if this is the only conflicting pair in the supertile",           this.plConflSngl);
		gd.addCheckbox    ("Only merge conflicting planes only if there are just two planes in the supertile",              this.plConflSnglPair);

		gd.addNumericField("Consider merging plane if it is foreground and maximal strength below this",                    this.plWeakFgStrength,  6);
		gd.addNumericField("Remove these strongest from foreground when determining the maximal strength",                  this.plWeakFgOutliers,  0);
		gd.addNumericField("Relax cost requirements when merging with weak foreground",                                     this.plWeakFgRelax,  6);

		gd.addNumericField("Maximal real-world thickness of merged overlapping planes (meters)",                            this.plThickWorld,  6);
		gd.addNumericField("Maximal real-world merged thickness for conflicting planes",                                    this.plThickWorldConfl,  6);
		gd.addNumericField("Relax cost requirements when adding exclusive links to complete squares and triangles",          this.plRelaxComplete,  6);
		gd.addNumericField("Relax cost requirements more during the second pass",                                           this.plRelaxComplete2,  6);

		gd.addMessage     ("---  ---");
		gd.addNumericField("Maximal ratio of Z to allow plane merging",                                              this.plMaxZRatio,  6);
		gd.addNumericField("Maximal disparity of one of the planes to apply  maximal ratio",                         this.plMaxDisp,  6);
		gd.addNumericField("When merging with neighbors cut the tail that is worse than scaled best",                this.plCutTail,  6);
		gd.addNumericField("Set cutoff value level not less than this",                                              this.plMinTail,  6);


		gd.addTab         ("Regenerate", "Parameters to regenerate planes using preliminary tile connections");
		gd.addMessage     ("--- Parameters to regenerate planes using preliminary tile connections  ---");
		gd.addCheckbox    ("Enable planes tiles selection regeneration hinted by supertile neighbors",               this.plDiscrEn);
		gd.addNumericField("Maximal disparity difference from the plane to consider tile",                           this.plDiscrTolerance,  6);
		gd.addNumericField("Parallel move known planes around original know value for the best overall fit",          this.plDiscrDispRange,  6);
		gd.addNumericField("Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)",          this.plDiscrSteps,  0);
		gd.addNumericField("What plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined",           this.plDiscrMode,  0);
		gd.addNumericField("Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)",          this.plDiscrVarFloor,  6);
		gd.addNumericField("Gaussian sigma to compare how measured data is attracted to planes",                     this.plDiscrSigma,  6);
		gd.addNumericField("Sigma to blur histograms while re-discriminating",                                       this.plDiscrBlur,  6);
		gd.addNumericField("Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others",          this.plDiscrExclusivity,  6);
		gd.addNumericField("For second pass if exclusivity > 1.0 - will assign only around strong neighbors",          this.plDiscrExclus2,  6);
		gd.addCheckbox    ("When growing selection do not allow any offenders around (false - more these than others)",           this.plDiscrStrict);
		gd.addNumericField("Attraction to different planes correlation that is too high for re-discrimination",          this.plDiscrCorrMax,  6);
		gd.addNumericField("Attraction to different planes correlation that is high enough to merge",                this.plDiscrCorrMerge,  6);
		gd.addNumericField("If offender has this number of tiles (including center) the cell can not be used",           this.plDiscrSteal,  0);
		gd.addNumericField("Only use tiles within this range from original selection",                               this.plDiscrGrown,  0);
		gd.addNumericField("Remove outliers from the final selection that have distance more than scaled median",          this.plDiscrXMedian,  6);

		gd.addTab         ("Merge costs", "Planes merge costs");
		gd.addMessage     ("--- Planes merge costs ---");
		gd.addNumericField(" Disparity (pix) - closer cost will use more of the real world, farther - disparity",          this.plCostDist,  6);
		gd.addNumericField("Cost of merge quality sqrt(weighted*equal) in disparity space",                          this.plCostKrq,     6);
		gd.addNumericField("Cost of merge quality average of weighted and equal weight in disparity space",          this.plCostKrqEq,   6);
		gd.addNumericField("Cost of merge quality sqrt(weighted*equal) in world space",                              this.plCostWrq,     6);
		gd.addNumericField("Cost of merge quality average of weighted and equal weight in world space",              this.plCostWrqEq,   6);
		gd.addNumericField("Cost of sin squared between normals",                                                    this.plCostSin2,    6);
		gd.addNumericField("Cost of squared relative plane-to-other-center distances",                               this.plCostRdist2,  6);

		gd.addCheckbox    ("Resolve dual triangles conflict (odoodo)",                                               this.plConflDualTri);
		gd.addCheckbox    ("Resolve multiple odo triangles conflicts",                                               this.plConflMulti);
		gd.addCheckbox    ("Resolve diagonal (ood) conflicts",                                                       this.plConflDiag);
		gd.addCheckbox    ("Resolve all conflicts around a supertile",                                               this.plConflStar);
		gd.addNumericField("How far to look around when calculationg connection cost",                               this.plStarSteps,  0);
		gd.addNumericField("When calculating cost for the connections scale 4 ortho neighbors",                      this.plStarOrtho,  6);
		gd.addNumericField("When calculating cost for the connections scale 4 diagonal neighbors",                   this.plStarDiag,  6);
		gd.addNumericField("Divide cost by number of connections to this power",                                     this.plStarPwr,  6);
		gd.addNumericField("Use this power of tile weight when calculating connection cost",                         this.plStarWeightPwr,  6);
		gd.addNumericField("Balance weighted density against density. 0.0 - density, 1.0 - weighted density",           this.plWeightToDens,  6);
		gd.addNumericField("Raise value of each tile before averaging",                                              this.plStarValPwr,  6);
		gd.addNumericField("When resolving double triangles allow minor degradation (0.0 - strict)",                 this.plDblTriLoss,  6);
		gd.addCheckbox    ("Allow more conflicts if overall cost is reduced",                                        this.plNewConfl);
		gd.addNumericField("aximal number of simultaneous connection changes around one tile (0 - any)",             this.plMaxChanges,  0);

		gd.addCheckbox    ("Keep only mutual links, remove weakest if conflict",                                     this.plMutualOnly);

		gd.addCheckbox    ("Add diagonals to full squares",                                                          this.plFillSquares);
		gd.addCheckbox    ("Add ortho to 45-degree corners",                                                         this.plCutCorners);
		gd.addCheckbox    ("Add hypotenuse connection if both legs exist",                                           this.plHypotenuse);

		gd.addNumericField("Relative (to average neighbor) weight of the measured plane when combing with neighbors",               this.plPull,  6);
		gd.addNumericField("0.0: 8 neighbors pull 8 times as 1, 1.0 - same as 1",                                    this.plNormPow,  6);
		gd.addNumericField("Maximal number of smoothing iterations for each step",                                   this.plIterations,  0);
		gd.addCheckbox    ("Do not update supertile if any of connected is not good (false: just skip that neighbor)",           this.plStopBad);
		gd.addNumericField("Maximal step difference (1/power of 10)",                                                this.plPrecision,  0);

		gd.addNumericField("Relative weight of center plane when splitting into pairs",                              this.plSplitPull,  6);
		gd.addNumericField("Minimal number of neighbors to split plane in pairs",                                    this.plSplitMinNeib,  0);
		gd.addNumericField("Minimal weight of split plains to show",                                                 this.plSplitMinWeight,  6);
		gd.addNumericField("Minimal split quality to show",                                                          this.plSplitMinQuality,  6);
		gd.addCheckbox    ("Apply plane split to pairs",                                                             this.plSplitApply);
		gd.addCheckbox    ("Allow tiles to belong to both planes of the pair",                                       this.plNonExclusive);
		gd.addCheckbox    ("Allow other tiles from the same supertile",                                              this.plUseOtherPlanes);
		gd.addCheckbox    ("Allow parallel shift of the specified planes before adding",                             this.plAllowParallel);
		gd.addNumericField("Maximal normalized tile disparity difference from the plane to consider",                this.plMaxDiff,  6);
		gd.addNumericField("Maximal difference of the added tile ratio to the average  disparity difference",          this.plOtherDiff,  6);
		gd.addCheckbox    ("Separate tiles for split planes by X, Y",                                                this.plSplitXY);
		gd.addNumericField("Disparity tolerance when separating by X, Y",                                            this.plSplitXYTolerance,  6);

		gd.addCheckbox    ("Fuse planes together (off for debug only)",                                              this.plFuse);
		gd.addCheckbox    ("Keep unconnected supertiles",                                                            this.plKeepOrphans);
		gd.addNumericField("Minimal strength unconnected supertiles to keep",                                        this.plMinOrphan,  6);

		gd.addTab         ("Snap", "Snap to planes");
		gd.addMessage     ("--- Snap to planes ---");
		gd.addNumericField("Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength",               this.plSnapDispAny,  6);
		gd.addNumericField("Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength",           this.plSnapStrengthAny,  6);
		gd.addNumericField("Maximal negative disparity difference from the best match",                              this.plSnapNegAny,  6);
		gd.addNumericField("Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength",               this.plSnapDispMax,  6);
		gd.addNumericField("Maximal disparity diff. by weight product to snap to plane",                             this.plSnapDispWeight,  6);
		gd.addNumericField("Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest",          this.plSnapZeroMode,  0);

		gd.addCheckbox    ("Use planes selection masks (generated when splitting to intersecting pairs",             this.msUseSel);
		gd.addCheckbox    ("Divide plane strengths by ellipsoid area",                                               this.msDivideByArea);
		gd.addNumericField("Scale projection of the plane ellipsoid",                                                this.msScaleProj,  6);
		gd.addNumericField("Spread this fraction of the ellipsoid weight among extended (double) supertile",          this.msFractUni,  6);

		gd.addTab         ("Assign", "Tiles assignment");
		gd.addMessage     ("--- Tiles assignment ---");
		gd.addCheckbox    ("Do not assign tiles to the surface edges (not having all 8 neighbors)",                     this.tsNoEdge);
		gd.addCheckbox    ("Only assign outside of 8x8 center if no suitable alternative",                              this.tsUseCenter);
		gd.addNumericField("Maximal disparity difference when assigning tiles",                                         this.tsMaxDiff,  6);
		gd.addNumericField("Minimal disparity difference to be considered as a competitor surface",                     this.tsMinDiffOther,  6);
		gd.addNumericField("Minimal tile correlation strength to be assigned",                                          this.tsMinStrength,  6);
		gd.addNumericField("Maximal tile correlation strength to be assigned",                                          this.tsMaxStrength,  6);
		gd.addNumericField("Minimal surface strength at the tile location",                                             this.tsMinSurface,  6);
		gd.addNumericField("Allowed tile disparity correction: 1 increase, 2 - decrease, 3 - both directions",          this.tsMoveDirs,  0);
		gd.addNumericField("Raise surface strengths ratio to this power when comparing candidates",                     this.tsSurfStrPow,  6);
		gd.addNumericField("Add to strengths when calculating pull of assigned tiles",                                  this.tsAddStrength,  6);
		gd.addNumericField("Radius of influence (in tiles) of the previously assigned tiles",                           this.tsSigma,  6);
		gd.addNumericField("Maximal relative to radius distance to calculate influence",                                this.tsNSigma,  6);
		gd.addNumericField(" Additional pull of each surface ",                                                         this.tsMinPull,  6);
		gd.addNumericField("Minimal ratio of the best surface candidate to the next one to make selection",             this.tsMinAdvantage,  6);

		gd.addNumericField("Minimal size of a cluster to keep",                                                         this.tsClustSize,  0);
		gd.addNumericField("Minimal total weight of a cluster to keep",                                                 this.tsClustWeight,  6);
		gd.addNumericField("Minimal number of neighbors of unassigned tile to join (the farthest)",                     this.tsMinNeib,  0);
		gd.addNumericField("Maximal strength of the surrounded unassigned tile to join",                                this.tsMaxSurStrength,  6);
		gd.addCheckbox    ("Include disabled tiles/borders when counting assigned neighbors",                           this.tsCountDis);

		gd.addCheckbox    ("Assign tiles that were used to generate planes",                                            this.tsEnPlaneSeed);
		gd.addCheckbox    ("Allow assignment only surface",                                                             this.tsEnOnly);
		gd.addCheckbox    ("Grow the only surface assignments",                                                         this.tsEnGrow);
		gd.addNumericField("Maximal strength when growing the only surfaces",                                           this.tsGrowStrength,  6);
		gd.addCheckbox    ("Grow over strong if disparity matches",                                                     this.tsGrowStrong);
		gd.addNumericField("Minimal strength to continue grow with disparity match",                                    this.tsContStrength,  6);
		gd.addNumericField("Maximal normalized disparity error to grow over strong tiles",                              this.tsContDiff,  6);

		gd.addCheckbox    ("Allow assignment to the nearest surface with no competitors",                               this.tsEnSingle);
		gd.addCheckbox    ("Allow assignment when several surfaces fit",                                                this.tsEnMulti);
		gd.addCheckbox    ("Remove weak clusters before growing",                                                       this.tsRemoveWeak1);
		gd.addCheckbox    ("Assign tiles that have neighbors to the lowest disparity",                                  this.tsGrowSurround);
		gd.addCheckbox    ("Remove weak clusters after growing",                                                        this.tsRemoveWeak2);

		gd.addCheckbox    ("Repeat multi-choice assignment while succeeding",                                           this.tsLoopMulti);
		gd.addCheckbox    ("Show results of tiles to surfaces assignment",                                              this.tsShow);
		gd.addNumericField("Number of clusters to keep",                                                                this.tsNumClust,  0);

		gd.addNumericField("Which assignments to match +1 - combo, +2 grown single, +4 plane seeds",                    this.tsConsensMode,  0);
		gd.addNumericField("Minimal number of assignments to agree",                                                    this.tsConsensAgree,  0);

		gd.addMessage     ("--- Tile assignment parameters ---");
		gd.addNumericField("Minimal foreground/ background separation to look for weak FG edge",                        this.taMinFgBg,    6);
		gd.addNumericField("Minimal foreground edge strength (stronger edges will have proportionally smaller costs)",           this.taMinFgEdge,  6);
		gd.addNumericField("Minimal surface separation that requires color change",                                     this.taMinColSep,  6);
		gd.addNumericField("Minimal color variation (larger proportionally reduces cost)",                              this.taMinColDiff, 6);
		gd.addNumericField("Disparity difference limit (to handle outliers)",                                           this.taOutlier, 6);
		gd.addNumericField("Strength power when calculating disparity error",                                           this.taDiffPwr, 6);
		gd.addNumericField("Strength power when calculating disparity error over best",                                 this.taBestPwr, 6);
		gd.addNumericField("Strength power when calculating disparity error for group of 9",                            this.taDiff9Pwr, 6);
		gd.addNumericField("Gaussian sigma to blur color difference between tiles along each direction",                this.taColSigma, 6);
		gd.addNumericField("Relative amount of the blurred color difference in the mixture",                            this.taColFraction, 6);

		gd.addNumericField("Cost of a tile that is not assigned",                                                       this.taCostEmpty,  6);
		gd.addNumericField("Cost of a tile not having any neighbor in particular direction",                            this.taCostNoLink,  6);
		gd.addNumericField("Cost of a tile switching to a neighbor that does not have a link",                          this.taCostSwitch,  6);
		gd.addNumericField("Cost of a tile switching to a disconnected neighbor divided by a color",                    this.taCostColor,  6);
		gd.addNumericField("Cost of a weighted normalized tile disparity error",                                        this.taCostDiff,  6);
		gd.addNumericField("Cost of a weighted normalized tile disparity error above best surface",                     this.taCostDiffBest,  6);
		gd.addNumericField("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",          this.taCostDiff9,  6);
		gd.addNumericField("Cost of a weak foreground edge",                                                            this.taCostWeakFgnd,  6);
		gd.addNumericField("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",                      this.taCostFlaps,  6);
		gd.addNumericField("Cost of a measurement layer not having same layer in the same location or near",            this.taCostMismatch,  6);

		gd.addCheckbox    ("Cost of a tile that is not assigned",                                                       this.taEnEmpty);
		gd.addCheckbox    ("Cost of a tile not having any neighbor in particular direction",                            this.taEnNoLink);
		gd.addCheckbox    ("Cost of a tile switching to a neighbor that does not have a link",                          this.taEnSwitch);
		gd.addCheckbox    ("Cost of a tile switching to a disconnected neighbor divided by a color",                    this.taEnColor);
		gd.addCheckbox    ("Cost of a weighted normalized tile disparity error",                                        this.taEnDiff);
		gd.addCheckbox    ("Cost of a weighted normalized tile disparity error above best surface",                     this.taEnDiffBest);
		gd.addCheckbox    ("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",          this.taEnDiff9);
		gd.addCheckbox    ("Cost of a weak foreground edge",                                                            this.taEnWeakFgnd);
		gd.addCheckbox    ("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",                      this.taEnFlaps);
		gd.addCheckbox    ("Cost of a measurement layer not having same layer in the same location or near",            this.taEnMismatch);

		gd.addTab         ("GPU", "Parameters for GPU development");
		gd.addMessage     ("--- GPU processing parameters ---");
		gd.addNumericField("Correlation radius",                                                                        this.gpu_corr_rad, 0, 6,"pix",
				"Size of the 2D correlation - maximal radius = 7 corresponds to full 15x15 pixel tile");
		gd.addNumericField("Correlation weight R",                                                                      this.gpu_weight_r, 4, 6,"",
				"Weight of R for composite 2D correlation (green weight is 1.0 -gpu_weight_r - gpu_weight_b");
		gd.addNumericField("Correlation weight B",                                                                      this.gpu_weight_b, 4, 6,"",
				"Weight of R for composite 2D correlation (green weight is 1.0 -gpu_weight_r - gpu_weight_b");
		gd.addMessage     ("--- LPF sigmas used for images/textures ---");
		gd.addNumericField("Color LPF sigma R",                                                                         this.gpu_sigma_r, 4, 6,"pix",
				"LPF sigma to process color components during aberration correction");
		gd.addNumericField("Color LPF sigma B",                                                                         this.gpu_sigma_b, 4, 6,"pix",
				"LPF sigma to process color components during aberration correction");
		gd.addNumericField("Color LPF sigma G",                                                                         this.gpu_sigma_g, 4, 6,"pix",
				"LPF sigma to process color components during aberration correction");
		gd.addNumericField("Monochrome LPF sigma",                                                                      this.gpu_sigma_m, 4, 6,"pix",
				"LPF sigma to process monochrome (e.g.LWIR) during aberration correction");
		gd.addMessage     ("--- LPF sigmas used for 2D phase correlation ---");
		gd.addNumericField("LPF sigma Red/Blue (extra)",                                                                this.gpu_sigma_rb_corr, 4, 6,"pix",
				"LPF sigma to apply to the composite R+B 2D correlation before adding G (next LPF will be applied to R+B+G correlation");
		gd.addNumericField("LPF sigma for correlation, color",                                                          this.gpu_sigma_corr, 4, 6,"pix",
				"LPF sigma to apply to the composite 2D correlation for RGB images");
		gd.addNumericField("LPF sigma for correlation, mono",                                                           this.gpu_sigma_corr_m, 4, 6,"pix",
				"LPF sigma to apply to the composite 2D correlation for monochrome images");
		gd.addNumericField("Fat zero (absolute) for phase correlation of color images",                                 this.gpu_fatz, 4, 6,"",
				"Add squared fat zero to the sum of squared amplitudes, color images");
		gd.addNumericField("Fat zero (absolute) for phase correlation of monochrome images",                            this.gpu_fatz_m, 4, 6,"",
				"Add squared fat zero to the sum of squared amplitudes, monochrome images");
		gd.addMessage     ("--- GPU WOI selection ---");

		gd.addCheckbox    ("Use following WOI for GPU processing (unchecked - ignore WOI dimensions)",                  this.gpu_woi);

		gd.addNumericField("WOI left",                                                                                  this.gpu_woi_tx, 0, 6,"tiles",
				"Left WOI margin, in tiles (0..323");
		gd.addNumericField("WOI top",                                                                                   this.gpu_woi_ty, 0, 6,"tiles",
				"Top WOI margin, in tiles (0..241");
		gd.addNumericField("WOI width",                                                                                 this.gpu_woi_twidth, 0, 6,"tiles",
				"WOI width, in tiles (1..324");
		gd.addNumericField("WOI height",                                                                                this.gpu_woi_theight, 0, 6,"tiles",
				"WOI height, in tiles (1..242");
		gd.addCheckbox    ("Select circle/ellipse within the rectanghular WOI",                                         this.gpu_woi_round);
		gd.addCheckbox    ("Debug feature - save calculated ports X,Y to compare with Java-generated",                  this.gpu_save_ports_xy);
		gd.addCheckbox    ("Show Java-generated textures from non-overlapping in GPU (will not generate if false)",     this.gpu_show_jtextures);
		gd.addCheckbox    ("Show low-res data for macro (will not generate if false)",                                  this.gpu_show_extra);

		gd.addCheckbox    ("Accelerate tile processor for the main quad camera",                                        this.gpu_use_main);
		gd.addCheckbox    ("Accelerate tile processor for the main quad camera in macro mode",                          this.gpu_use_main_macro);
		gd.addCheckbox    ("Accelerate tile processor for the main quad camera for field calibration",                  this.gpu_use_main_adjust);
		gd.addCheckbox    ("Accelerate tile processor for the aux (lwir) quad camera",                                  this.gpu_use_aux);
		gd.addCheckbox    ("Accelerate tile processor for the aux (lwir) quad camera in macro mode",                    this.gpu_use_aux_macro);
		gd.addCheckbox    ("Accelerate tile processor for the aux (lwir) quad camera for field calibration",            this.gpu_use_aux_adjust);
		
		gd.addTab         ("LWIR", "parameters for LWIR/EO 8-camera rig");
		this.lwir.dialogQuestions(gd);
		gd.addTab         ("Debug", "Other debug images");
		gd.addMessage     ("--- Other debug images ---");
		//	clt_parameters.debug_initial_discriminate, // final boolean    debug_initial_discriminate,
		gd.addCheckbox    ("Debug initial clusterization of the supertile tiles",                                    this.debug_initial_discriminate);
		gd.addCheckbox    ("Test new mode after migration",                                                          this.dbg_migrate);

		gd.addNumericField("Temporay exit stage (0- normal execution)",                                              this.dbg_early_exit,  0,6,"","Temporary exit at intermediate stage (0 - normal)");
		gd.addCheckbox    ("Show first infinity scan",                                                               this.show_first_bg, "Show results of the first calculated background scan");

		gd.addCheckbox    ("Show extrinsic adjustment differences",                                                  this.show_extrinsic);
		gd.addCheckbox    ("Show 'ortho_combine'",                                                                   this.show_ortho_combine);
		gd.addCheckbox    ("Show 'refine_disparity_supertiles'",                                                     this.show_refine_supertiles);
		gd.addCheckbox    ("Show 'bgnd_nonbgnd'",                                                                    this.show_bgnd_nonbgnd);
		gd.addCheckbox    ("Show 'FilterScan'",                                                                      this.show_filter_scan);
		gd.addCheckbox    ("Show 'combo_scan' (combined multiple scans)",                                            this.show_combined);
		gd.addCheckbox    ("Show 'unique_scan' (removed already measured tiles with the same disparity)",            this.show_unique);
		gd.addCheckbox    ("Show supertile disparity histograms ",                                                   this.show_histograms);
		gd.addCheckbox    ("Show debug images during initial refinement",                                            this.show_init_refine,"Show debug images during initial refinement");
		gd.addCheckbox    ("Show debug images during disparity expansion",                                           this.show_expand);
		gd.addCheckbox    ("Show prepareExpandVariant when elevating variant",                                       this.show_variant);
		gd.addCheckbox    ("Show debug images related to retrying far tiles near foreground",                        this.show_retry_far);
		gd.addCheckbox    ("Show debug images related to macro correlation",                                         this.show_macro);
		gd.addCheckbox    ("Show 'shells'",                                                                          this.show_shells);
		gd.addCheckbox    ("show 'neighbors'",                                                                       this.show_neighbors);
		gd.addCheckbox    ("Show 'flaps-dirs'",                                                                      this.show_flaps_dirs);
		gd.addCheckbox    ("Show 'first_N_clusters'",                                                                this.show_first_clusters);
		gd.addCheckbox    ("Show planes",                                                                            this.show_planes);

		gd.addMessage     ("Unity up vector in camera coordinate system (x - right, y - up, z - to camera): {"+
				this.vertical_xyz[0]+","+          this.vertical_xyz[1]+","+          this.vertical_xyz[2]+"}");


		//  			gd.buildDialog();
		gd.showDialog();
		//  			System.out.println("gd.wasCanceled="+gd.wasCanceled());
		/*
  			if (true) return false;
  			WindowTools.addScrollBars(gd);
  			gd.showDialog();
		 */
		if (gd.wasCanceled()) return false;
		this.disparity=             gd.getNextNumber();
		this.transform_size=  (int) gd.getNextNumber();
		this.clt_window=      (int) gd.getNextNumber();
		this.shift_x =              gd.getNextNumber();
		this.shift_y =              gd.getNextNumber();

		this.tileStep=        (int) gd.getNextNumber();
		this.iclt_mask=       (int) gd.getNextNumber();
		this.tileX=           (int) gd.getNextNumber();
		this.tileY=           (int) gd.getNextNumber();
		this.dbg_mode=        (int) gd.getNextNumber();
		this.ishift_x=        (int) gd.getNextNumber();
		this.ishift_y=        (int) gd.getNextNumber();
		this.fat_zero =             gd.getNextNumber();
		this.fat_zero_mono =        gd.getNextNumber();
		this.corr_sigma =           gd.getNextNumber();
		this.corr_sigma_mono =      gd.getNextNumber();
		this.scale_strength_main =  gd.getNextNumber();
		this.scale_strength_aux =   gd.getNextNumber();
		this.norm_kern=             gd.getNextBoolean();
		this.gain_equalize=         gd.getNextBoolean();
		this.colors_equalize=       gd.getNextBoolean();
		this.nosat_equalize=        gd.getNextBoolean();
		this.sat_level=             gd.getNextNumber();
		this.max_overexposure=      gd.getNextNumber();

		this.novignetting_r=        gd.getNextNumber();
		this.novignetting_g=        gd.getNextNumber();
		this.novignetting_b=        gd.getNextNumber();
		this.scale_r=               gd.getNextNumber();
		this.scale_g=               gd.getNextNumber();
		this.scale_b=               gd.getNextNumber();
		this.vignetting_max=        gd.getNextNumber();
		this.vignetting_range=      gd.getNextNumber();
		this.kernel_step=     (int) gd.getNextNumber();
		this.z_correction=          gd.getNextNumber();
		this.correlate=             gd.getNextBoolean();
		this.corr_mask=       (int) gd.getNextNumber();
		this.corr_sym=              gd.getNextBoolean();
		this.corr_keep=             gd.getNextBoolean();
		this.corr_show=             gd.getNextBoolean();
		this.corr_mismatch=         gd.getNextBoolean();
		this.corr_offset=           gd.getNextNumber();
		this.corr_red=              gd.getNextNumber();
		this.corr_blue=             gd.getNextNumber();
		this.corr_normalize=        gd.getNextBoolean();
		this.min_corr=              gd.getNextNumber();
		this.min_corr_normalized=   gd.getNextNumber();
		this.max_corr_sigma=        gd.getNextNumber();
		this.max_corr_radius=       gd.getNextNumber();
		this.corr_magic_scale=      gd.getNextNumber();
		this.corr_select =    (int) gd.getNextNumber();


		this.img_dtt.dialogAnswers(gd);
		this.rig.dialogAnswers(gd);
		this.poles.dialogAnswers(gd);

		this.max_corr_double=       gd.getNextBoolean();
		this.corr_mode=       (int) gd.getNextNumber();
		this.corr_border_contrast=  gd.getNextNumber();

		//  			          this.tile_task_op=    (int) gd.getNextNumber();
		this.tile_task_op = ImageDtt.setOrthoLines      (          this.tile_task_op, gd.getNextBoolean());
		this.tile_task_op = ImageDtt.setForcedDisparity (          this.tile_task_op, gd.getNextBoolean());
		this.tile_task_op = ImageDtt.setImgMask         (          this.tile_task_op, (int) gd.getNextNumber());
		this.tile_task_op = ImageDtt.setPairMask        (          this.tile_task_op, (int) gd.getNextNumber());

		this.tile_task_wl=    (int) gd.getNextNumber();
		this.tile_task_wt=    (int) gd.getNextNumber();
		this.tile_task_ww=    (int) gd.getNextNumber();
		this.tile_task_wh=    (int) gd.getNextNumber();
		this.min_shot=              gd.getNextNumber();
		this.scale_shot=            gd.getNextNumber();
		this.diff_sigma=            gd.getNextNumber();
		this.diff_threshold=        gd.getNextNumber();
		this.diff_gauss=            gd.getNextBoolean();
		this.min_agree=             gd.getNextNumber();
		this.dust_remove=           gd.getNextBoolean();
		this.black_back=            gd.getNextBoolean();
		this.keep_weights=          gd.getNextBoolean();
		this.sharp_alpha=           gd.getNextBoolean();
		this.alpha0=                gd.getNextNumber();
		this.alpha1=                gd.getNextNumber();
		this.gen_chn_stacks=        gd.getNextBoolean();
		this.gen_chn_img=           gd.getNextBoolean();
		this.gen_4_img=             gd.getNextBoolean();
		this.show_nonoverlap=       gd.getNextBoolean();
		this.show_overlap=          gd.getNextBoolean();
		this.show_rgba_color=       gd.getNextBoolean();
		this.show_map=              gd.getNextBoolean();
		this.show_corr=             gd.getNextBoolean();
		this.disp_scan_start=       gd.getNextNumber();
		this.disp_scan_step=        gd.getNextNumber();
		this.disp_scan_count= (int) gd.getNextNumber();

		this.fine_dbg=              gd.getNextBoolean();
		this.fine_corr_x_0=         gd.getNextNumber();
		this.fine_corr_y_0=         gd.getNextNumber();
		this.fine_corr_x_1=         gd.getNextNumber();
		this.fine_corr_y_1=         gd.getNextNumber();
		this.fine_corr_x_2=         gd.getNextNumber();
		this.fine_corr_y_2=         gd.getNextNumber();
		this.fine_corr_x_3=         gd.getNextNumber();
		this.fine_corr_y_3=         gd.getNextNumber();
		this.fine_corr_ignore=      gd.getNextBoolean();
		this.fine_corr_apply=       gd.getNextBoolean();

		this.fcorr_radius=          gd.getNextNumber();
		this.fcorr_min_strength=    gd.getNextNumber();
		this.fcorr_disp_diff=       gd.getNextNumber();
		this.fcorr_quadratic=       gd.getNextBoolean();
		this.fcorr_ignore=          gd.getNextBoolean();

		this.fcorr_inf_strength=    gd.getNextNumber();
		this.fcorr_inf_diff=        gd.getNextNumber();
		this.fcorr_inf_quad=        gd.getNextBoolean();
		this.fcorr_inf_vert=        gd.getNextBoolean();

		this.inf_disp_apply=        gd.getNextBoolean();
		this.inf_repeat=      (int) gd.getNextNumber();
		this.inf_iters=       (int) gd.getNextNumber();
		this.inf_final_diff=        gd.getNextNumber();
		this.inf_far_pull=          gd.getNextNumber();

		this.inf_str_pow=           gd.getNextNumber();
		this.inf_smpl_side=   (int) gd.getNextNumber();
		this.inf_smpl_num=    (int) gd.getNextNumber();
		this.inf_smpl_rms=          gd.getNextNumber();

		this.ih_smpl_step=    (int) gd.getNextNumber();
		this.ih_disp_min=           gd.getNextNumber();
		this.ih_disp_step=          gd.getNextNumber();
		this.ih_num_bins=     (int) gd.getNextNumber();
		this.ih_sigma=              gd.getNextNumber();
		this.ih_max_diff=           gd.getNextNumber();
		this.ih_min_samples=  (int) gd.getNextNumber();
		this.ih_norm_center=        gd.getNextBoolean();
		this.inf_restore_disp=      gd.getNextBoolean();

		this.ly_lma_ers =           gd.getNextBoolean();
		this.ly_gt_strength=        gd.getNextNumber();
		this.ly_gt_use_wnd=         gd.getNextBoolean();
		this.ly_gt_rms=             gd.getNextNumber();

		this.lylw_inf_en =          gd.getNextBoolean();
		this.lylw_aztilt_en =       gd.getNextBoolean();
		this.lylw_diff_roll_en =    gd.getNextBoolean();
		this.lylw_focalLength =     gd.getNextBoolean();
		this.lylw_com_roll =        gd.getNextBoolean();
		this.lylw_par_sel=    (int) gd.getNextNumber();

		this.ly_marg_fract=         gd.getNextNumber();
		this.ly_on_scan=            gd.getNextBoolean();
		this.ly_inf_en=             gd.getNextBoolean();
		this.ly_min_forced=   (int) gd.getNextNumber();
		this.ly_aztilt_en=          gd.getNextBoolean();
		this.ly_diff_roll_en=       gd.getNextBoolean();
		this.ly_focalLength=        gd.getNextBoolean();
		this.ly_com_roll=           gd.getNextBoolean();
		this.ly_ers_rot=            gd.getNextBoolean();
		this.ly_ers_forw=           gd.getNextBoolean();
		this.ly_ers_side=           gd.getNextBoolean();
		this.ly_ers_vert=           gd.getNextBoolean();
		this.ly_par_sel=      (int) gd.getNextNumber();
		this.ly_debug_level=  (int) gd.getNextNumber();

		this.ly_right_left=         gd.getNextBoolean();

		this.ly_per_quad=     (int) gd.getNextNumber();
		this.ly_per_quad_r=         gd.getNextNumber();
		this.ly_inf=          (int) gd.getNextNumber();
		this.ly_inf_r=              gd.getNextNumber();
		this.ly_inf_scale=    (int) gd.getNextNumber();
		this.ly_inf_scale_r=        gd.getNextNumber();

		this.ly_inf_frac=           gd.getNextNumber();
		this.ly_inf_max_disparity=  gd.getNextNumber();

		this.ly_inf_disp=           gd.getNextBoolean();
		this.ly_inf_force=          gd.getNextBoolean();
		this.ly_poly=               gd.getNextBoolean();

		this.ly_smpl_side=    (int) gd.getNextNumber();
		this.ly_smpl_num=     (int) gd.getNextNumber();
		this.ly_smpl_rms=           gd.getNextNumber();
		this.ly_disp_var=           gd.getNextNumber();
		this.ly_disp_rvar=          gd.getNextNumber();
		this.ly_disp_var_gt=        gd.getNextNumber();
		this.ly_disp_rvar_gt=       gd.getNextNumber();
		this.ly_norm_disp=          gd.getNextNumber();
		this.lym_overexp=           gd.getNextNumber();
		this.lym_update_disp=       gd.getNextBoolean();
		this.lym_iter=        (int) gd.getNextNumber();
		this.lym_change=            gd.getNextNumber();
		this.lym_change_aux=        gd.getNextNumber();

		this.lym_poly_change=       gd.getNextNumber();

		this.lyf_filter=            gd.getNextBoolean();
		this.lyf_smpl_side=   (int) gd.getNextNumber();
		this.lyf_rms_max=           gd.getNextNumber();
		this.lyf_frac_keep=         gd.getNextNumber();
		this.lyf_min_samples= (int) gd.getNextNumber();
		this.lyf_norm_center=       gd.getNextBoolean();
		this.ly_corr_scale=         gd.getNextNumber();

		this.lyr_filter_ds=         gd.getNextBoolean();
		this.lyr_filter_lyf=        gd.getNextBoolean();

		//  			          this.fcorr_sample_size= (int)gd.getNextNumber();
		//  			          this.fcorr_mintiles= (int)  gd.getNextNumber();
		//  			          this.fcorr_reloutliers=     gd.getNextNumber();
		//  			          this.fcorr_sigma=           gd.getNextNumber();

		///  			          this.corr_magic_scale=      gd.getNextNumber();

		this.show_textures=         gd.getNextBoolean();
		this.debug_filters=         gd.getNextBoolean();
		this.min_smth=              gd.getNextNumber();
		this.sure_smth=             gd.getNextNumber();
		this.bgnd_range=            gd.getNextNumber();
		this.other_range=           gd.getNextNumber();

		this.ex_strength=           gd.getNextNumber();
		this.ex_nstrength=          gd.getNextNumber();

		this.ex_over_bgnd=          gd.getNextBoolean();
		this.ex_min_over=           gd.getNextNumber();

		this.pt_super_trust=        gd.getNextNumber();
		this.pt_keep_raw_fg=        gd.getNextBoolean();
		this.pt_scale_pre=          gd.getNextNumber();
		this.pt_scale_post=         gd.getNextNumber();

		this.bgnd_sure=             gd.getNextNumber();
		this.bgnd_maybe=            gd.getNextNumber();
		this.min_clstr_seed=  (int) gd.getNextNumber();
		this.min_clstr_lone=  (int) gd.getNextNumber();
		this.min_clstr_weight=      gd.getNextNumber();
		this.min_clstr_max=         gd.getNextNumber();

		this.fill_gaps=       (int) gd.getNextNumber();
		this.fill_final=      (int) gd.getNextNumber();
		this.min_clstr_block= (int) gd.getNextNumber();
		this.bgnd_grow=       (int) gd.getNextNumber();

		this.ortho_old=             gd.getNextBoolean();
		this.ortho_min_hor=         gd.getNextNumber();
		this.ortho_min_vert=        gd.getNextNumber();
		this.ortho_asym=            gd.getNextNumber();
		this.ortho_over4=           gd.getNextNumber();

		this.ortho_sustain=         gd.getNextNumber();
		this.ortho_run=       (int) gd.getNextNumber();
		this.ortho_minmax=          gd.getNextNumber();
		this.ortho_bridge=    (int) gd.getNextNumber();
		this.ortho_rms=             gd.getNextNumber();
		this.ortho_half_length=(int)gd.getNextNumber();
		this.ortho_mix=             gd.getNextNumber();

		this.or_hor=                gd.getNextBoolean();
		this.or_vert=               gd.getNextBoolean();
		this.or_sigma=              gd.getNextNumber();
		this.or_sharp=              gd.getNextNumber();
		this.or_scale=              gd.getNextNumber();
		this.or_offset=             gd.getNextNumber();
		this.or_asym=               gd.getNextNumber();
		this.or_threshold=          gd.getNextNumber();
		this.or_absHor=             gd.getNextNumber();
		this.or_absVert=            gd.getNextNumber();
		this.or_maxDisp=            gd.getNextNumber();

		this.poles_fix=             gd.getNextBoolean();
		this.poles_len=        (int)gd.getNextNumber();
		this.poles_ratio=           gd.getNextNumber();
		this.poles_min_strength=    gd.getNextNumber();
		this.poles_force_disp=      gd.getNextBoolean();

		this.max_clusters=    (int) gd.getNextNumber();
		this.remove_scans=          gd.getNextBoolean();
		this.output_x3d=            gd.getNextBoolean();
		this.output_obj=            gd.getNextBoolean();
		this.correct_distortions=   gd.getNextBoolean();
		this.show_triangles=        gd.getNextBoolean();
		this.avg_cluster_disp=      gd.getNextBoolean();
		this.maxDispTriangle=       gd.getNextNumber();
		this.infinityDistance=      gd.getNextNumber();
		this.min_bgnd_tiles=  (int) gd.getNextNumber();
		this.shUseFlaps=            gd.getNextBoolean();
		this.shAggrFade=            gd.getNextBoolean();
		this.shMinArea=       (int) gd.getNextNumber();
		this.shMinStrength=         gd.getNextNumber();
		this.tiRigidVertical=       gd.getNextNumber();
		this.tiRigidHorizontal=     gd.getNextNumber();
		this.tiRigidDiagonal=       gd.getNextNumber();
		this.tiStrengthOffset=      gd.getNextNumber();
		this.tiDispScale=           gd.getNextNumber();
		this.tiDispPow=             gd.getNextNumber();
		this.tiDispPull=            gd.getNextNumber();
		this.tiDispPullPreFinal=    gd.getNextNumber();
		this.tiDispPullFinal=       gd.getNextNumber();
		this.tiBreakNorm=           gd.getNextNumber();
		this.tiBreak3=              gd.getNextNumber();
		this.tiBreak31=             gd.getNextNumber();
		this.tiBreak21=             gd.getNextNumber();
		this.tiBreakFar=            gd.getNextNumber();
		this.tiBreakNear=           gd.getNextNumber();
		this.tiBreakMode=     (int) gd.getNextNumber();
		this.tiBreakSame=           gd.getNextNumber();
		this.tiBreakTurn=           gd.getNextNumber();

		this.tiHealPreLast=         gd.getNextNumber();
		this.tiHealLast=            gd.getNextNumber();
		this.tiHealSame=      (int) gd.getNextNumber();

		this.tiIterations=    (int) gd.getNextNumber();
		this.tiPrecision=     (int) gd.getNextNumber();
		this.tiNumCycles=     (int) gd.getNextNumber();

		this.stUseRefine=           gd.getNextBoolean();
		this.stUsePass2=            gd.getNextBoolean();
		this.stUseRender=           gd.getNextBoolean();

		this.stShow=                gd.getNextBoolean();
		this.stSize=          (int) gd.getNextNumber();
		this.stStepFar=             gd.getNextNumber();
		this.stStepNear=            gd.getNextNumber();
		this.stStepThreshold=       gd.getNextNumber();
		this.stMinDisparity=        gd.getNextNumber();
		//  			          this.stMaxDisparity=        gd.getNextNumber();
		this.stSigma=               gd.getNextNumber();
		this.stMinBgDisparity=      gd.getNextNumber();
		this.stMinBgFract=          gd.getNextNumber();
		this.stUseDisp=             gd.getNextNumber();
		this.stStrengthScale=       gd.getNextNumber();

		this.stGrowSel=       (int) gd.getNextNumber();
		this.stMeasSel=       (int) gd.getNextNumber();
		this.stSmallDiff=           gd.getNextNumber();
		this.stHighMix=             gd.getNextNumber();

		this.outlierStrength=      gd.getNextNumber();
		this.outlierDiff=          gd.getNextNumber();
		this.outlierDiffPos=       gd.getNextNumber();
		this.outlierDiffNeg=       gd.getNextNumber();

		this.combine_refine=        gd.getNextBoolean();

		this.combine_min_strength=  gd.getNextNumber();
		this.combine_min_hor=       gd.getNextNumber();
		this.combine_min_vert=      gd.getNextNumber();
		///--------------------
		this.stSmplMode=            gd.getNextBoolean();
		this.mlfp.strength_sure=    gd.getNextNumber();
		this.mlfp.strength_floor=   gd.getNextNumber();
		this.mlfp.strength_pow=     gd.getNextNumber();

		this.mlfp.smplSide=   (int) gd.getNextNumber();
		this.mlfp.smplNum=    (int) gd.getNextNumber();
		this.mlfp.smplRms=          gd.getNextNumber();
		this.mlfp.smplWnd=          gd.getNextBoolean();

		this.mlfp.max_abs_tilt=     gd.getNextNumber();
		this.mlfp.max_rel_tilt=     gd.getNextNumber();
		this.mlfp.damp_tilt =       gd.getNextNumber();
		this.mlfp.min_tilt_disp=    gd.getNextNumber();
		this.mlfp.transition=       gd.getNextNumber();
		this.mlfp.far_mode=   (int) gd.getNextNumber();
		this.mlfp.far_power=        gd.getNextNumber();

		//  			          this.unique_tolerance=      gd.getNextNumber();
		this.grow_sweep=      (int) gd.getNextNumber();
		this.grow_disp_max=         gd.getNextNumber();
		this.grow_disp_trust=       gd.getNextNumber();
		this.grow_disp_step=        gd.getNextNumber();
		this.grow_min_diff=         gd.getNextNumber();
		this.grow_retry_far=        gd.getNextBoolean();
		this.grow_pedantic=         gd.getNextBoolean();
		this.grow_retry_inf=        gd.getNextBoolean();

		this.gr_new_expand=         gd.getNextBoolean();
		//  			          this.gr_max_expand=   (int) gd.getNextNumber();
		this.fds_str_floor=         gd.getNextNumber();
		this.gr_ovrbg_cmb=          gd.getNextNumber();
		this.gr_ovrbg_cmb_hor=      gd.getNextNumber();
		this.gr_ovrbg_cmb_vert=     gd.getNextNumber();
		this.gr_ovrbg_filtered=     gd.getNextNumber();
		this.fds_str_pow=           gd.getNextNumber();
		this.fds_smpl_side=   (int) gd.getNextNumber();
		this.fds_smpl_num=    (int) gd.getNextNumber();
		this.fds_smpl_rms=          gd.getNextNumber();
		this.fds_smpl_rel_rms=      gd.getNextNumber();
		this.fds_smpl_wnd=          gd.getNextBoolean();
		this.fds_abs_tilt=          gd.getNextNumber();
		this.fds_rel_tilt=          gd.getNextNumber();

		this.per_filter=               gd.getNextBoolean();
		this.per_trustedCorrelation=   gd.getNextNumber();
		this.per_initial_diff=         gd.getNextNumber();
		this.per_strength_floor=       gd.getNextNumber();
		this.per_strength_max_over=    gd.getNextNumber();
		this.per_min_period=           gd.getNextNumber();
		this.per_min_num_periods=(int) gd.getNextNumber();
		this.per_disp_tolerance=       gd.getNextNumber();
		this.per_disp_match=           gd.getNextNumber();
		this.per_strong_match_inc=     gd.getNextNumber();

		this.mc_disp8_step=         gd.getNextNumber();
		this.mc_disp8_trust=        gd.getNextNumber();
		this.mc_strength=           gd.getNextNumber();
		this.mc_unique_tol=         gd.getNextNumber();

		this.mc_trust_fin=          gd.getNextNumber();
		this.mc_trust_sigma=        gd.getNextNumber();
		this.mc_ortho_weight=       gd.getNextNumber();
		this.mc_diag_weight=        gd.getNextNumber();
		this.mc_gap=                gd.getNextNumber();

		this.mc_weight_var=         gd.getNextNumber();
		this.mc_weight_Y=           gd.getNextNumber();
		this.mc_weight_RBmG=        gd.getNextNumber();

		this.gr_min_new=      (int) gd.getNextNumber();
		this.gr_var_new_sngl=       gd.getNextBoolean();
		this.gr_var_new_fg=         gd.getNextBoolean();
		this.gr_var_all_fg=         gd.getNextBoolean();
		this.gr_var_new_bg=         gd.getNextBoolean();
		this.gr_var_all_bg=         gd.getNextBoolean();
		this.gr_var_next=           gd.getNextBoolean();
		this.gr_num_steps=    (int) gd.getNextNumber();
		this.gr_steps_over=   (int) gd.getNextNumber();
		this.gr_smpl_size=    (int) gd.getNextNumber();
		this.gr_min_pnts=     (int) gd.getNextNumber();
		this.gr_use_wnd=            gd.getNextBoolean();
		this.gr_tilt_damp=          gd.getNextNumber();
		this.gr_split_rng=          gd.getNextNumber();
		this.gr_same_rng=           gd.getNextNumber();
		this.gr_diff_cont=          gd.getNextNumber();
		this.gr_abs_tilt=           gd.getNextNumber();
		this.gr_rel_tilt=           gd.getNextNumber();
		this.gr_smooth=       (int) gd.getNextNumber();
		this.gr_fin_diff=           gd.getNextNumber();
		this.gr_unique_tol=         gd.getNextNumber();
		this.gr_unique_pretol=      gd.getNextNumber();

		this.ft_mod_strength=               gd.getNextBoolean();
		this.ft_clusterize_by_highest=      gd.getNextBoolean();
		this.ft_clust_sigma=                gd.getNextNumber();
		this.ft_disp_arange_vert=           gd.getNextNumber();
		this.ft_disp_rrange_vert=           gd.getNextNumber();
		this.ft_disp_arange_hor=            gd.getNextNumber();
		this.ft_disp_rrange_hor=            gd.getNextNumber();
		this.ft_tolerance_above_near=       gd.getNextNumber();
		this.ft_tolerance_below_near=       gd.getNextNumber();
		this.ft_tolerance_above_far=        gd.getNextNumber();
		this.ft_tolerance_below_far=        gd.getNextNumber();
		this.ft_hor_vert_overlap=     (int) gd.getNextNumber();
		this.ft_used_companions=      (int) gd.getNextNumber();
		this.ft_used_true_companions= (int) gd.getNextNumber();

		this.plPreferDisparity=     gd.getNextBoolean();
		this.plDispNorm=            gd.getNextNumber();
		this.plFrontoTol =          gd.getNextNumber();
		this.plFrontoRms =          gd.getNextNumber();
		this.plFrontoOffs =         gd.getNextNumber();
		this.PlFrontoPow =          gd.getNextNumber();

		this.plBlurBinVert=         gd.getNextNumber();
		this.plBlurBinHor=          gd.getNextNumber();
		this.plMaxDiffVert=         gd.getNextNumber();
		this.plMaxDiffHor=          gd.getNextNumber();
		this.plInitPasses=    (int) gd.getNextNumber();

		this.plMinPoints=     (int) gd.getNextNumber();
		this.plTargetEigen=         gd.getNextNumber();
		this.plFractOutliers=       gd.getNextNumber();
		this.plMaxOutliers=   (int) gd.getNextNumber();
		this.plMinStrength=         gd.getNextNumber();
		this.plMaxEigen=            gd.getNextNumber();
		this.plEigenFloor=          gd.getNextNumber();
		this.plEigenStick=          gd.getNextNumber();
		this.plBadPlate=            gd.getNextNumber();
		this.plDbgMerge=            gd.getNextBoolean();
		this.plWorstWorsening=      gd.getNextNumber();
		this.plWorstWorsening2=     gd.getNextNumber();
		this.plWorstEq=             gd.getNextNumber();
		this.plWorstEq2=            gd.getNextNumber();
		this.plOKMergeEigen=        gd.getNextNumber();
		this.plMaxWorldSin2=        gd.getNextNumber();
		this.pl2dForSin=            gd.getNextNumber();
		this.plWeakWorsening=       gd.getNextNumber();
		this.plMaxOverlap=          gd.getNextNumber();

		this.plWeakWeight=          gd.getNextNumber();
		this.plWeakEigen=           gd.getNextNumber();
		this.plWeakWeight2=         gd.getNextNumber();
		this.plWeakEigen2=          gd.getNextNumber();
		this.plSumThick=            gd.getNextNumber();
		this.plNeNeibCost=          gd.getNextNumber();
		this.plNeOwn=               gd.getNextNumber();

		this.plExNeibCost=          gd.getNextNumber();
		this.plExNeibSmooth=        gd.getNextNumber();
		this.plMergeCostStar=       gd.getNextNumber();
		this.plMergeCost=           gd.getNextNumber();

		this.plConflMerge=          gd.getNextBoolean();
		this.plConflRelax=          gd.getNextNumber();
		this.plConflSngl=           gd.getNextBoolean();
		this.plConflSnglPair=       gd.getNextBoolean();

		this.plWeakFgStrength=      gd.getNextNumber();
		this.plWeakFgOutliers=(int) gd.getNextNumber();
		this.plWeakFgRelax=         gd.getNextNumber();

		this.plThickWorld=          gd.getNextNumber();
		this.plThickWorldConfl=     gd.getNextNumber();
		this.plRelaxComplete=       gd.getNextNumber();
		this.plRelaxComplete2=      gd.getNextNumber();

		this.plMaxZRatio=           gd.getNextNumber();
		this.plMaxDisp=             gd.getNextNumber();
		this.plCutTail=             gd.getNextNumber();
		this.plMinTail=             gd.getNextNumber();

		this.plDiscrEn=             gd.getNextBoolean();
		this.plDiscrTolerance=      gd.getNextNumber();
		this.plDiscrDispRange=      gd.getNextNumber();
		this.plDiscrSteps=    (int) gd.getNextNumber();
		//  			          this.plDiscrVariants= (int) gd.getNextNumber();
		this.plDiscrMode=     (int) gd.getNextNumber();
		this.plDiscrVarFloor=       gd.getNextNumber();
		this.plDiscrSigma=          gd.getNextNumber();
		this.plDiscrBlur=           gd.getNextNumber();
		this.plDiscrExclusivity=    gd.getNextNumber();
		this.plDiscrExclus2=        gd.getNextNumber();
		this.plDiscrStrict=         gd.getNextBoolean();
		this.plDiscrCorrMax=        gd.getNextNumber();
		this.plDiscrCorrMerge=      gd.getNextNumber();
		this.plDiscrSteal=    (int) gd.getNextNumber();
		this.plDiscrGrown=    (int) gd.getNextNumber();
		this.plDiscrXMedian=        gd.getNextNumber();

		this.plCostDist=            gd.getNextNumber();
		this.plCostKrq=             gd.getNextNumber();
		this.plCostKrqEq=           gd.getNextNumber();
		this.plCostWrq=             gd.getNextNumber();
		this.plCostWrqEq=           gd.getNextNumber();
		this.plCostSin2=            gd.getNextNumber();
		this.plCostRdist2=          gd.getNextNumber();

		this.plConflDualTri=        gd.getNextBoolean();
		this.plConflMulti=          gd.getNextBoolean();
		this.plConflDiag=           gd.getNextBoolean();
		this.plConflStar=           gd.getNextBoolean();
		this.plStarSteps=     (int) gd.getNextNumber();
		this.plStarOrtho=           gd.getNextNumber();
		this.plStarDiag=            gd.getNextNumber();
		this.plStarPwr=             gd.getNextNumber();
		this.plStarWeightPwr=       gd.getNextNumber();
		this.plWeightToDens=        gd.getNextNumber();
		this.plStarValPwr=          gd.getNextNumber();
		this.plDblTriLoss=          gd.getNextNumber();
		this.plNewConfl=            gd.getNextBoolean();
		this.plMaxChanges=    (int) gd.getNextNumber();

		this.plMutualOnly=          gd.getNextBoolean();

		this.plFillSquares=         gd.getNextBoolean();
		this.plCutCorners=          gd.getNextBoolean();
		this.plHypotenuse=          gd.getNextBoolean();

		this.plPull=                gd.getNextNumber();
		this.plNormPow=             gd.getNextNumber();
		this.plIterations=    (int) gd.getNextNumber();
		this.plStopBad=             gd.getNextBoolean();
		this.plPrecision=     (int) gd.getNextNumber();

		this.plSplitPull=           gd.getNextNumber();
		this.plSplitMinNeib=  (int) gd.getNextNumber();
		this.plSplitMinWeight=      gd.getNextNumber();
		this.plSplitMinQuality=     gd.getNextNumber();
		this.plSplitApply=          gd.getNextBoolean();
		this.plNonExclusive=        gd.getNextBoolean();
		this.plUseOtherPlanes=      gd.getNextBoolean();
		this.plAllowParallel=       gd.getNextBoolean();
		this.plMaxDiff=             gd.getNextNumber();
		this.plOtherDiff=           gd.getNextNumber();
		this.plSplitXY=             gd.getNextBoolean();
		this.plSplitXYTolerance=    gd.getNextNumber();

		this.plFuse=                gd.getNextBoolean();
		this.plKeepOrphans=         gd.getNextBoolean();
		this.plMinOrphan=           gd.getNextNumber();

		this.plSnapDispAny=         gd.getNextNumber();
		this.plSnapStrengthAny=     gd.getNextNumber();
		this.plSnapNegAny=          gd.getNextNumber();
		this.plSnapDispMax=         gd.getNextNumber();
		this.plSnapDispWeight=      gd.getNextNumber();
		this.plSnapZeroMode=  (int) gd.getNextNumber();

		this.msUseSel=              gd.getNextBoolean();
		this.msDivideByArea=        gd.getNextBoolean();
		this.msScaleProj=           gd.getNextNumber();
		this.msFractUni=            gd.getNextNumber();

		this.tsNoEdge=              gd.getNextBoolean();
		this.tsUseCenter=           gd.getNextBoolean();
		this.tsMaxDiff=             gd.getNextNumber();
		this.tsMinDiffOther=        gd.getNextNumber();
		this.tsMinStrength=         gd.getNextNumber();
		this.tsMaxStrength=         gd.getNextNumber();
		this.tsMinSurface=          gd.getNextNumber();
		this.tsMoveDirs=      (int) gd.getNextNumber();
		this.tsSurfStrPow=          gd.getNextNumber();
		this.tsAddStrength=         gd.getNextNumber();
		this.tsSigma=               gd.getNextNumber();
		this.tsNSigma=              gd.getNextNumber();
		this.tsMinPull=             gd.getNextNumber();
		this.tsMinAdvantage=        gd.getNextNumber();

		this.tsClustSize=     (int) gd.getNextNumber();
		this.tsClustWeight=         gd.getNextNumber();
		this.tsMinNeib=       (int) gd.getNextNumber();
		this.tsMaxSurStrength=      gd.getNextNumber();
		this.tsCountDis=            gd.getNextBoolean();
		this.tsReset              = false;  // Reset tiles to surfaces assignment

		this.tsEnPlaneSeed=         gd.getNextBoolean();
		this.tsEnOnly=              gd.getNextBoolean();
		this.tsEnGrow=              gd.getNextBoolean();
		this.tsGrowStrength=        gd.getNextNumber();
		this.tsGrowStrong=          gd.getNextBoolean();
		this.tsContStrength=        gd.getNextNumber();
		this.tsContDiff=            gd.getNextNumber();

		this.tsEnSingle=            gd.getNextBoolean();
		this.tsEnMulti=             gd.getNextBoolean();
		this.tsRemoveWeak1=         gd.getNextBoolean();
		this.tsGrowSurround=        gd.getNextBoolean();
		this.tsRemoveWeak2=         gd.getNextBoolean();

		this.tsLoopMulti=           gd.getNextBoolean();
		this.tsShow               = gd.getNextBoolean();
		this.tsNumClust =     (int) gd.getNextNumber();

		this.tsConsensMode =  (int) gd.getNextNumber();
		this.tsConsensAgree = (int) gd.getNextNumber();

		this.taMinFgBg=             gd.getNextNumber();
		this.taMinFgEdge=           gd.getNextNumber();
		this.taMinColSep=           gd.getNextNumber();
		this.taMinColDiff=          gd.getNextNumber();
		this.taOutlier=             gd.getNextNumber();
		this.taDiffPwr=             gd.getNextNumber();
		this.taBestPwr=             gd.getNextNumber();
		this.taDiff9Pwr=            gd.getNextNumber();
		this.taColSigma=            gd.getNextNumber();
		this.taColFraction=         gd.getNextNumber();


		this.taCostEmpty=           gd.getNextNumber();
		this.taCostNoLink=          gd.getNextNumber();
		this.taCostSwitch=          gd.getNextNumber();
		this.taCostColor=           gd.getNextNumber();
		this.taCostDiff=            gd.getNextNumber();
		this.taCostDiffBest=        gd.getNextNumber();
		this.taCostDiff9=           gd.getNextNumber();
		this.taCostWeakFgnd=        gd.getNextNumber();
		this.taCostFlaps=           gd.getNextNumber();
		this.taCostMismatch=        gd.getNextNumber();

		this.taEnEmpty=             gd.getNextBoolean();
		this.taEnNoLink=            gd.getNextBoolean();
		this.taEnSwitch=            gd.getNextBoolean();
		this.taEnColor=             gd.getNextBoolean();
		this.taEnDiff=              gd.getNextBoolean();
		this.taEnDiffBest=          gd.getNextBoolean();
		this.taEnDiff9=             gd.getNextBoolean();
		this.taEnWeakFgnd=          gd.getNextBoolean();
		this.taEnFlaps=             gd.getNextBoolean();
		this.taEnMismatch=          gd.getNextBoolean();

		this.gpu_corr_rad =   (int) gd.getNextNumber();
		this.gpu_weight_r =         gd.getNextNumber();
		this.gpu_weight_b =         gd.getNextNumber();
		this.gpu_sigma_r =          gd.getNextNumber();
		this.gpu_sigma_b =          gd.getNextNumber();
		this.gpu_sigma_g =          gd.getNextNumber();
		this.gpu_sigma_m =          gd.getNextNumber();
		this.gpu_sigma_rb_corr =    gd.getNextNumber();
		this.gpu_sigma_corr =       gd.getNextNumber();
		this.gpu_sigma_corr_m =     gd.getNextNumber();
		this.gpu_fatz =             gd.getNextNumber();
		this.gpu_fatz_m =           gd.getNextNumber();

		this.gpu_woi=               gd.getNextBoolean();
		this.gpu_woi_tx =     (int) gd.getNextNumber();
		this.gpu_woi_ty =     (int) gd.getNextNumber();
		this.gpu_woi_twidth = (int) gd.getNextNumber();
		this.gpu_woi_theight =(int) gd.getNextNumber();
		this.gpu_woi_round=         gd.getNextBoolean();
		this.gpu_save_ports_xy=     gd.getNextBoolean();
		this.gpu_show_jtextures=    gd.getNextBoolean();
		this.gpu_show_extra=        gd.getNextBoolean();

		this.gpu_use_main=          gd.getNextBoolean();
		this.gpu_use_main_macro=    gd.getNextBoolean();
		this.gpu_use_main_adjust=   gd.getNextBoolean();
		this.gpu_use_aux=           gd.getNextBoolean();
		this.gpu_use_aux_macro=     gd.getNextBoolean();
		this.gpu_use_aux_adjust=    gd.getNextBoolean();
		
		this.lwir.dialogAnswers(gd);

		this.debug_initial_discriminate= gd.getNextBoolean();
		this.dbg_migrate=                gd.getNextBoolean();

		this.dbg_early_exit = (int) gd.getNextNumber();
		this.show_first_bg=         gd.getNextBoolean();

		this.show_extrinsic=        gd.getNextBoolean();
		this.show_ortho_combine=    gd.getNextBoolean();
		this.show_refine_supertiles=gd.getNextBoolean();
		this.show_bgnd_nonbgnd=     gd.getNextBoolean(); // first on second pass
		this.show_filter_scan=      gd.getNextBoolean(); // first on refine
		this.show_combined=         gd.getNextBoolean();
		this.show_unique=           gd.getNextBoolean();
		this.show_histograms=       gd.getNextBoolean();
		this.show_init_refine=      gd.getNextBoolean();
		this.show_expand=           gd.getNextBoolean();
		this.show_variant=          gd.getNextBoolean();
		this.show_retry_far=        gd.getNextBoolean();
		this.show_macro=            gd.getNextBoolean();
		this.show_shells=           gd.getNextBoolean();
		this.show_neighbors=        gd.getNextBoolean();
		this.show_flaps_dirs=       gd.getNextBoolean();
		this.show_first_clusters=   gd.getNextBoolean();
		this.show_planes=           gd.getNextBoolean();

		return true;
	}


	public boolean showTsDialog() {
		GenericDialog gd = new GenericDialog("Set CLT tiles to surfaces assignment parameters");
		gd.addCheckbox    ("Do not assign tiles to the surface edges (not having all 8 neighbors)",                     this.tsNoEdge);
		gd.addCheckbox    ("Only assign outside of 8x8 center if no suitable alternative",                              this.tsUseCenter);
		gd.addNumericField("Maximal disparity difference when assigning tiles",                                         this.tsMaxDiff,  6);
		gd.addNumericField("Minimal disparity difference to be considered as a competitor surface",                     this.tsMinDiffOther,  6);
		gd.addNumericField("Minimal tile correlation strength to be assigned",                                          this.tsMinStrength,  6);
		gd.addNumericField("Maximal tile correlation strength to be assigned",                                          this.tsMaxStrength,  6);
		gd.addNumericField("Minimal surface strength at the tile location",                                             this.tsMinSurface,  6);
		gd.addNumericField("Allowed tile disparity correction: 1 increase, 2 - decrease, 3 - both directions",          this.tsMoveDirs,  0);
		gd.addNumericField("Raise surface strengths ratio to this power when comparing candidates",                     this.tsSurfStrPow,  6);
		gd.addNumericField("Add to strengths when calculating pull of assigned tiles",                                  this.tsAddStrength,  6);
		gd.addNumericField("Radius of influence (in tiles) of the previously assigned tiles",                           this.tsSigma,  6);
		gd.addNumericField("Maximal relative to radius distance to calculate influence",                                this.tsNSigma,  6);
		gd.addNumericField("Additional pull of each surface ",                                                          this.tsMinPull,  6);
		gd.addNumericField("Minimal ratio of the best surface candidate to the next one to make selection",             this.tsMinAdvantage,  6);
		gd.addNumericField("Minimal size of a cluster to keep",                                                         this.tsClustSize,  0);
		gd.addNumericField("Minimal total weight of a cluster to keep",                                                 this.tsClustWeight,  6);
		gd.addNumericField("Minimal number of neighbors of unassigned tile to join (the farthest)",                     this.tsMinNeib,  0);
		gd.addNumericField("Maximal strength of the surrounded unassigned tile to join",                                this.tsMaxSurStrength,  6);
		gd.addCheckbox    ("Include disabled tiles/borders when counting assigned neighbors",                           this.tsCountDis);
		gd.addCheckbox    ("Reset tiles to surfaces assignment",                                              false);

		gd.addCheckbox    ("Assign tiles that were used to generate planes",                                            this.tsEnPlaneSeed);
		gd.addCheckbox    ("Allow assignment only surface",                                                             this.tsEnOnly);
		gd.addCheckbox    ("Grow the only surface assignments",                                                         this.tsEnGrow);
		gd.addNumericField("Maximal strength when growing the only surfaces",                                           this.tsGrowStrength,  6);
		gd.addCheckbox    ("Grow over strong if disparity matches",                                                     this.tsGrowStrong);
		gd.addNumericField("Minimal strength to continue grow with disparity match",                                    this.tsContStrength,  6);
		gd.addNumericField("Maximal normalized disparity error to grow over strong tiles",                              this.tsContDiff,  6);

		gd.addCheckbox    ("Allow assignment to the nearest surface with no competitors",                               this.tsEnSingle);
		gd.addCheckbox    ("Allow assignment when several surfaces fit",                                                this.tsEnMulti);
		gd.addCheckbox    ("Remove weak clusters before growing",                                                       this.tsRemoveWeak1);
		gd.addCheckbox    ("Assign tiles that have neighbors to the lowest disparity",                                  this.tsGrowSurround);
		gd.addCheckbox    ("Remove weak clusters after growing",                                                        this.tsRemoveWeak2);

		gd.addCheckbox    ("Repeat multi-choice assignment while succeeding",                                           this.tsLoopMulti);
		gd.addCheckbox    ("Show results of tiles to surfaces assignment",                                              this.tsShow);
		gd.addNumericField("Number of clusters to keep",                                                                this.tsNumClust,  0);

		gd.addNumericField("Which assignments to match +1 - combo, +2 grown single, +4 plane seeds",                    this.tsConsensMode,  0);
		gd.addNumericField("Minimal number of assignments to agree",                                                    this.tsConsensAgree,  0);

		gd.addMessage     ("--- Tile assignment parameters ---");
		gd.addNumericField("Minimal foreground/ background separation to look for weak FG edge",                        this.taMinFgBg,    6);
		gd.addNumericField("Minimal foreground edge strength (stronger edges will have proportionally smaller costs)",           this.taMinFgEdge,  6);
		gd.addNumericField("Minimal surface separation that requires color change",                                     this.taMinColSep,  6);
		gd.addNumericField("Minimal color variation (larger proportionally reduces cost)",                              this.taMinColDiff, 6);
		gd.addNumericField("Disparity difference limit (to handle outliers)",                                           this.taOutlier, 6);
		gd.addNumericField("Strength power when calculating disparity error",                                           this.taDiffPwr, 6);
		gd.addNumericField("Strength power when calculating disparity error over best",                                 this.taBestPwr, 6);
		gd.addNumericField("Strength power when calculating disparity error for group of 9",                            this.taDiff9Pwr, 6);
		gd.addNumericField("Gaussian sigma to blur color difference between tiles along each direction",                this.taColSigma, 6);
		gd.addNumericField("Relative amount of the blurred color difference in the mixture",                            this.taColFraction, 6);

		gd.addNumericField("Cost of a tile that is not assigned",                                                       this.taCostEmpty,  6);
		gd.addNumericField("Cost of a tile not having any neighbor in particular direction",                            this.taCostNoLink,  6);
		gd.addNumericField("Cost of a tile switching to a neighbor that does not have a link",                          this.taCostSwitch,  6);
		gd.addNumericField("Cost of a tile switching to a disconnected neighbor divided by a color",                    this.taCostColor,  6);
		gd.addNumericField("Cost of a weighted normalized tile disparity error",                                        this.taCostDiff,  6);
		gd.addNumericField("Cost of a weighted normalized tile disparity error above best surface",                     this.taCostDiffBest,  6);
		gd.addNumericField("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",          this.taCostDiff9,  6);
		gd.addNumericField("Cost of a weak foreground edge",                                                            this.taCostWeakFgnd,  6);
		gd.addNumericField("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",                      this.taCostFlaps,  6);
		gd.addNumericField("Cost of a measurement layer not having same layer in the same location or near",            this.taCostMismatch,  6);

		gd.addCheckbox    ("Cost of a tile that is not assigned",                                                       this.taEnEmpty);
		gd.addCheckbox    ("Cost of a tile not having any neighbor in particular direction",                            this.taEnNoLink);
		gd.addCheckbox    ("Cost of a tile switching to a neighbor that does not have a link",                          this.taEnSwitch);
		gd.addCheckbox    ("Cost of a tile switching to a disconnected neighbor divided by a color",                    this.taEnColor);
		gd.addCheckbox    ("Cost of a weighted normalized tile disparity error",                                        this.taEnDiff);
		gd.addCheckbox    ("Cost of a weighted normalized tile disparity error above best surface",                     this.taEnDiffBest);
		gd.addCheckbox    ("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",          this.taEnDiff9);
		gd.addCheckbox    ("Cost of a weak foreground edge",                                                            this.taEnWeakFgnd);
		gd.addCheckbox    ("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",                      this.taEnFlaps);
		gd.addCheckbox    ("Cost of a measurement layer not having same layer in the same location or near",            this.taEnMismatch);

		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		this.tsNoEdge=              gd.getNextBoolean();
		this.tsUseCenter=           gd.getNextBoolean();
		this.tsMaxDiff=             gd.getNextNumber();
		this.tsMinDiffOther=        gd.getNextNumber();
		this.tsMinStrength=         gd.getNextNumber();
		this.tsMaxStrength=         gd.getNextNumber();
		this.tsMinSurface=          gd.getNextNumber();
		this.tsMoveDirs=      (int) gd.getNextNumber();
		this.tsSurfStrPow=          gd.getNextNumber();
		this.tsAddStrength=         gd.getNextNumber();
		this.tsSigma=               gd.getNextNumber();
		this.tsNSigma=              gd.getNextNumber();
		this.tsMinPull=             gd.getNextNumber();
		this.tsMinAdvantage=        gd.getNextNumber();

		this.tsClustSize=     (int) gd.getNextNumber();
		this.tsClustWeight=         gd.getNextNumber();
		this.tsMinNeib=       (int) gd.getNextNumber();
		this.tsMaxSurStrength=      gd.getNextNumber();
		this.tsCountDis=            gd.getNextBoolean();
		this.tsReset=               gd.getNextBoolean();

		this.tsEnPlaneSeed=         gd.getNextBoolean();
		this.tsEnOnly=              gd.getNextBoolean();
		this.tsEnGrow=              gd.getNextBoolean();
		this.tsGrowStrength=        gd.getNextNumber();
		this.tsGrowStrong=          gd.getNextBoolean();
		this.tsContStrength=        gd.getNextNumber();
		this.tsContDiff=            gd.getNextNumber();

		this.tsEnSingle=            gd.getNextBoolean();
		this.tsEnMulti=             gd.getNextBoolean();
		this.tsRemoveWeak1=         gd.getNextBoolean();
		this.tsGrowSurround=        gd.getNextBoolean();
		this.tsRemoveWeak2=         gd.getNextBoolean();

		this.tsLoopMulti=           gd.getNextBoolean();
		this.tsShow =               gd.getNextBoolean();
		this.tsNumClust=      (int) gd.getNextNumber();
		this.tsConsensMode =  (int) gd.getNextNumber();
		this.tsConsensAgree = (int) gd.getNextNumber();

		this.taMinFgBg=             gd.getNextNumber();
		this.taMinFgEdge=           gd.getNextNumber();
		this.taMinColSep=           gd.getNextNumber();
		this.taMinColDiff=          gd.getNextNumber();
		this.taOutlier=             gd.getNextNumber();
		this.taDiffPwr=             gd.getNextNumber();
		this.taBestPwr=             gd.getNextNumber();
		this.taDiff9Pwr=            gd.getNextNumber();
		this.taColSigma=            gd.getNextNumber();
		this.taColFraction=         gd.getNextNumber();

		this.taCostEmpty=           gd.getNextNumber();
		this.taCostNoLink=          gd.getNextNumber();
		this.taCostSwitch=          gd.getNextNumber();
		this.taCostColor=           gd.getNextNumber();
		this.taCostDiff=            gd.getNextNumber();
		this.taCostDiffBest=        gd.getNextNumber();
		this.taCostDiff9=           gd.getNextNumber();
		this.taCostWeakFgnd=        gd.getNextNumber();
		this.taCostFlaps=           gd.getNextNumber();
		this.taCostMismatch=        gd.getNextNumber();

		this.taEnEmpty=             gd.getNextBoolean();
		this.taEnNoLink=            gd.getNextBoolean();
		this.taEnSwitch=            gd.getNextBoolean();
		this.taEnColor=             gd.getNextBoolean();
		this.taEnDiff=              gd.getNextBoolean();
		this.taEnDiffBest=          gd.getNextBoolean();
		this.taEnDiff9=             gd.getNextBoolean();
		this.taEnWeakFgnd=          gd.getNextBoolean();
		this.taEnFlaps=             gd.getNextBoolean();
		this.taEnMismatch=          gd.getNextBoolean();
		return true;
	}
	public boolean modifyZCorr (String title) {
		if (z_corr_map == null){
			z_corr_map = new HashMap<String,Double>();
		}
		GenericDialog gd = new GenericDialog(title);
		if (z_corr_map.size() > 0) {
			gd.addMessage("Edit infinity disparity correction (in 1/m), set >= 1.0 to remove for the following");
			for (HashMap.Entry<String,Double> entry : z_corr_map.entrySet()){
				gd.addNumericField(entry.getKey(),   entry.getValue(), 9,12, "m-1");
			}
		}
		gd.addMessage("Add new infinity correction");
		gd.addStringField ("Timestamp string (seconds_microseconds):",                              "", 40);
		gd.addNumericField("Infinity correction (in 1/m)",                                          0,  9,12,"m-1");
		gd.addCheckbox    ("Clear list",                                                            false);
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		HashMap<String,Double> new_map = new HashMap<String,Double>();
		for (HashMap.Entry<String,Double> entry : z_corr_map.entrySet()){
			double d = gd.getNextNumber();
			if (d < 1.0) {
				new_map.put(entry.getKey(),d);
			}
		}
		String new_ts =  gd.getNextString();
		double d =       gd.getNextNumber();
		if (gd.getNextBoolean()){
			z_corr_map = new HashMap<String,Double>();
		} else {
			z_corr_map = new_map;
		}
		if(new_ts.length() > 0){
			z_corr_map.put(new_ts, d);
		}
		return true;
	}
	public boolean modifyInfCorr (String title) {
		if (infinity_distace_map == null){
			infinity_distace_map = new HashMap<String,Double>();
		}
		class StringDouble{
			String str;
			double dbl;
			StringDouble (String s, double d){
				this.str = s;
				this.dbl = d;
			}
		}
		ArrayList<StringDouble> sd_list = new ArrayList<StringDouble>();
		for (HashMap.Entry<String,Double> entry : infinity_distace_map.entrySet()){
			sd_list.add(new StringDouble(entry.getKey(),   entry.getValue()));
		}
		Collections.sort(sd_list, new Comparator<StringDouble>() {
			@Override
			public int compare(StringDouble lhs, StringDouble rhs) {
				// -1 - less than, 1 - greater than, 0 - equal, not inverted for ascending disparity
				return  lhs.str.compareTo(rhs.str);
			}
		});


		GenericDialog gd = new GenericDialog(title);
		if (infinity_distace_map.size() > 0) {
			gd.addMessage("Edit \"infinity\" (when mountains are too close) distance (in m), set == 0.0 to remove for the following");
			gd.addMessage("Press \"Cancel\" to exit");
			//				for (HashMap.Entry<String,Double> entry : infinity_distace_map.entrySet()){
			for (StringDouble entry : sd_list){
				gd.addNumericField(entry.str,   entry.dbl, 1,6, "m");
			}
		}
		gd.addMessage("Add new infinity correction");
		gd.addStringField ("Timestamp string (seconds_microseconds):",                              "", 40);
		gd.addNumericField("Infinity distance",                                          0,  1,6,"m");
		gd.addCheckbox    ("Clear list",                                                            false);
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		HashMap<String,Double> new_map = new HashMap<String,Double>();
		for (StringDouble entry : sd_list){
			double d = gd.getNextNumber();
			if (d != 0.0) {
				new_map.put(entry.str,d);
			}
		}


		String new_ts =  gd.getNextString();
		double d =       gd.getNextNumber();
		if (gd.getNextBoolean()){
			infinity_distace_map = new HashMap<String,Double>();
		} else {
			infinity_distace_map = new_map;
		}
		if(new_ts.length() > 0){
			infinity_distace_map.put(new_ts, d);
		}
		return true;
	}
}