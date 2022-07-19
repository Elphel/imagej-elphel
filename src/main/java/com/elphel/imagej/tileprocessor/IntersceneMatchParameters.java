package com.elphel.imagej.tileprocessor;
/**
 **
 ** IntersceneMatchParameters - Class for handling multiple configuration parameters
 ** related to the interscene match 
 **
 ** Copyright (C) 202 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  IntersceneMatchParameters.java is free software: you can redistribute it and/or modify
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

import java.awt.Color;
import java.io.IOException;
import java.util.Properties;
import java.util.StringTokenizer;

import com.elphel.imagej.common.GenericJTabbedDialog;

public class IntersceneMatchParameters {
	public static String [] MODES3D = {"RAW", "INF", "FG", "BG"}; // RAW:-1
	public static String [] MODES_AVI = {"RAW", "JPEG", "PNG"};
	// Maybe add parameters to make sure there is enough data? Enough in each zone? Enough spread?
	public  boolean force_ref_dsi =      false; // true;
	public  boolean force_orientations = false;
	public  int min_num_orient =         2; // make from parameters, should be >= 1
	public  int min_num_interscene =     2; // make from parameters, should be >= 1

	public  boolean generate_mapped =    true; // generate any of the sequences - Tiff, video, mono or stereo
	public  boolean reuse_video =        false; // dry-run video generation, just reassemble video fragments (active if !generate_mapped)
	public  boolean save_mapped_color =  false; // may be not needed when using AVI
	public  boolean save_mapped_mono =   true;  // may be not needed when using AVI
	public  boolean generate_raw =       false;
	public  boolean generate_inf =       false;
	public  boolean generate_fg =        true;
	public  boolean generate_bg =        false;

	public  boolean generate_stereo =    false; // applies to all stereobases, including 0
	
	public  boolean export_images =      true;  // pseudo-color 16-slice images (same disparity, COMBO_DSN_INDX_DISP_FG and COMBO_DSN_INDX_DISP_BG_ALL,
	public  boolean show_images =        false; // color, infinity
	public  boolean show_images_bgfg =   false; // bg and fg
	public  boolean show_images_mono =   false; // float, monochrome 16-slice images (same disparity, COMBO_DSN_INDX_DISP_FG and COMBO_DSN_INDX_DISP_BG_ALL,
	public  boolean show_color_nan =     true;  // use NAN background for color images (sharp, but distinct black)
	public  boolean show_mono_nan =      false; // use NAN background for monochrome images (sharp, but distinct black)
	
//	public  double [] stereo_bases =   {0.0, 200.0, 500.0, 1000.0}; 
	public  double [][] stereo_views = {  // base, up, back
			{   0.0,   0.0,    0.0},
			{ 200.0,   0.0,    0.0},
			{ 500.0,   0.0, 2000.0},
			{1000.0, 500.0, 3000.0}};
	
//	public  boolean [] generate_stereo_var = new boolean[stereo_bases.length];
	public  boolean [] generate_stereo_var = new boolean[stereo_views.length];

	// Other parameters
	public  int     min_num_scenes =     10;     // abandon series if there are less than this number of scenes in it 
	public  double  blur_egomotion =     2.0;
	
	
	public  boolean export_ranges =      true;
	public  boolean show_ranges =        true;
	
	public  boolean export_ml_files =    false; // export files for DNN training
	
	public  boolean stereo_merge =       true;
	public  int     stereo_gap =         32;    // pixels between right and left frames
	public  double  stereo_intereye =    63.5;  // mm 
	public  double  stereo_phone_width = 143.0; // mm if 0 - no padding. Will be padded only when merging

	public  int     extra_hor_tile =     15;
	public  int     extra_vert_tile =    10;
	public  boolean crop_3d =            true; // do not show extra of FG/BG views (currently they only ref scene has disparity)
	public  int     sensor_mask =        1; // -1 - all
	public  boolean merge_all =          true; // merge all 16 channels for 3D modes 

	public  boolean show_mapped_color =  true;
	public  boolean show_mapped_mono =   true;

	public  boolean gen_avi_color =      true; // will use save_mapped_color
	public  boolean gen_avi_mono =       true; // will use save_mapped_mono
	public  double  video_fps =          15; // 4x slowmo
	public  double  sensor_fps =         60; // sensor fps
	public  int     mode_avi =           0;  // 0 - raw, 1 - JPEG, 2 - PNG
	public  int     avi_JPEG_quality =  90;
	public  boolean run_ffmpeg =         true; // only after AVI
	public  String  video_ext =          ".webm";
	public  String  video_codec =        "vp8";
	public  int     video_crf =         40;    // lower - better, larger file size
	public  boolean remove_avi =         true; // remove avi after conversion to webm
	public  String  video_codec_combo =  "vp8"; // applies when combining videos
	public  int     video_crf_combo =   40;    // lower - better, larger file size applies when combining videos
	public  boolean um_mono =            true; // applies to both TIFF and AVI 
	public  double  um_sigma =          10;
	public  double  um_weight =          0.97; //
	public  boolean mono_fixed =         true; // normalize to fixed range when converting to 8 bits 
	public  double  mono_range =       500.0;  // monochrome full-scale range (+/- half)
	public  boolean anaglyth_en =        true; // applies to both TIFF and AVI 
	public static Color   anaglyph_left_default = new Color (255, 0, 0); // red
	public static Color   anaglyph_right_default =new Color (0, 255, 255); // cyan
	public Color   anaglyph_left =    anaglyph_left_default;
	public Color   anaglyph_right =   anaglyph_right_default; // cyan
	
	public  boolean annotate_color =     true; // annotate pseudo-color video frames with timestamps 
	public  boolean annotate_mono =      true; // annotate monochrome video frames with timestamps 
	public Color   annotate_color_color =new Color( 255, 255, 255); // greenish over "fire"
	public Color annotate_color_mono =   new Color( 255, 180,  50); // reddish over grey
	public boolean annotate_transparent_mono =  false; // // black if not transparent
	
	
	
	public  double  range_disparity_offset =   -0.08;
	public  double  range_min_strength = 0.5;
	public  double  range_max =       5000.0;
	
	
// Other parameters for filtering depth maps	
	public  int         num_bottom =                      3; // 6; // average this number of lowest disparity neighbors (of 8)
	public  int         num_passes =                    100;
	public  double      max_change =                      1.0e-3;
	public  double      min_disparity =                  -0.15; // 0.2;
	public  double      max_sym_disparity =               0.1; // 0.2;
	// 2 next are wrong currently - minimal strength is ~0.25
	public  double      min_strength_lma =                0.3;  // no real filtering
	public  double      min_strength_replace =            0.05; ///  0.14; /// Before /// - LWIR, after - RGB
	public  double      min_strength_blur =               0.06; ///  0.2;
	public  double      sigma =                           2; /// 5;
	public  int         num_blur =                        1; // 3;
	public  double      disparity_corr =                  0.0;
	public  int         outliers_nth_fromextrem =         1; // second from min/max - removes dual-tile max/mins
	public  double      outliers_tolerance_absolute =     0.2;
	public  double      outliers_tolerance_relative =     0.02;
	public  int         outliers_max_iter =             100;
	public  double      outliers_max_strength2 =          1.0; // 0.5; any
	public  int         outliers_nth_fromextrem2 =        0; // second from min/max - removes dual-tile max/mins
	public  double      outliers_tolerance_absolute2 =    0.5;
	public  double      outliers_tolerance_relative2 =    0.1;
	public  double      outliers_lma_max_strength =       0.4; // 0.5;
	public  double      outliers_max_strength =           0.1; ///  0.25;
	public  double      outliers_from_lma_max_strength =  0.8;
	public  int         search_radius =                   3;  // Search farther if no LMA neighbor is found closer. Original value - 1 (8 neighbors)
	public  boolean     remove_no_lma_neib =              true;
	public  double      diff_from_lma_pos =             100.0;
	public  double      diff_from_lma_neg =               2.0;
	public  int         outliers_lma_nth_fromextrem =     0; // 1; 
	public  int         filter_margin =                  -8; // 8; // pixels
	
	public  double      weak_tolerance_absolute=          0.25;
	public  double      weak_tolerance_relative=          0.025;
	public  int         weak_min_neibs =                  5;
	public  double      strong_strength=                  0.5;
	public  double      weak_strength=                    0.2; // none is below 
	
	

	// Some "AGC" to adjust how much to discard
	public  int     margin =                       1;      // do not use tiles if their centers are closer to the image edge
	public  int     sensor_mask_inter =           -1;      // bitmask of the sensors to use (-1 - all)
	public  boolean use_partial =                   true;  // find motion vectors for individual pairs, false - for sum only
	public  boolean run_poly =                      false; // not yet implemented
	public  double  centroid_radius =               4.0;   // 
	public  int     n_recenter =                    2;     // when cosine window, re-center window this many times
	// filtering motion vectors
	// TD accumulation of the inter-scene correlations  demonstrated artifacts (horizontally offset by 8 pixels
	// false maximum that is sharper than the real one. Still not understood - maybe float precision related. 
	public  double  td_weight =                     0.5;   // mix correlations accumulated in TD with 
	public  double  pd_weight =                     0.5;   // correlations (post) accumulated in PD
	public  boolean td_nopd_only =                  false; // true;  // only use TD accumulated data if no safe PD is available for the tile.

	public  double  min_str_fpn =                   0.2; // 0.25; // minimal correlation strength for all but TD-accumulated layer
	public  double  min_str_sum_fpn =               0.5; // 0.8;  // minimal correlation strength for TD-accumulated layer

	public  double  min_str =                       0.18;  // tiles w/o FPN: minimal correlation strength for all but TD-accumulated layer
	public  double  min_str_sum =                   0.33; // tiles w/o FPN: minimal correlation strength for TD-accumulated layer

	public  int     min_neibs =                     2;	   // minimal number of strong neighbors (> min_str)
	public  double  weight_zero_neibs =             0.2;   // Reduce weight for no-neib (1.0 for all 8)
	public  double  half_disparity =                5.0;   // Reduce weight twice for this disparity
	public  double  half_avg_diff =                 0.2;   // when L2 of x,y difference from average of neibs - reduce twice

	// Photometric calibration (move elsewhere)?
	public boolean  photo_en =                      false; // perform photogrammetric calibration to equalize pixel values
	public int      photo_num_full =                1;     // Number of full recalibrations with re-processing of the images  
	public int      photo_num_refines =             3;     // Calibrate, remove outliers, recalibrate, ... 
	public double   photo_min_strength =            0.0;   // maybe add to filter out weak tiles
	public double   photo_max_diff =                40.0;  // To filter mismatches. Normal (adjusted) have RMSE ~9
	public boolean  photo_debug =                   false; // Generate images and text
	
	
	
	// Detect initial match
	public  double  min_ref_str =                   0.22;  // For orientations: use only tiles of the reference scene DSI_MAIN is stronger  
	public  int     pix_step =                      4;     // Azimuth/tilt search step in pixels
	public  int     search_rad =                   10;     // Search radius in steps
	public  double  maybe_sum =                     1.0;   // minimal sum of strengths (will search for the best)
	public  double  sure_sum =                     50.0;   // definitely good sum of strengths (no more search) too high, will be FOM-defined
	public  double  maybe_avg =                     0.01;  // maybe average strength
	public  double  sure_avg =                      0.15;  // 0.015; // sure average strength disabling (was 0.03)
	public  double  max_search_rms =                1.5;   // good - 0.34, so-so - 0.999
	public  double  maybe_fom =                     1.0;   // good - 38, second good - 4.5
	public  double  sure_fom =                     12.0;   // good - 38, second good - 4.5
	
	
	// Reference scene disparity 
	public  boolean use_combo_dsi =                 true;  // use interscene DSI if available (instead of the single-scene)
	public  boolean use_lma_dsi =                   true;  // only use reference DSI tiles that have LMA (strong) disparity
	
	// Remove correlation caused by FPN
	public  boolean fpn_remove =                    true;  // only use reference DSI tiles that have LMA (strong) disparity
	public  double  fpn_max_offset =                8.0;   // pix - ignore larger FPN offsets
	public  double  fpn_radius =                    0.75;  // pix - zero around center 
	public  boolean fpn_ignore_border =             false; // only if fpn_mask != null - ignore tile if maximum touches fpn_mask
	
	// Remove moving objects (goal is not to detect slightest movement, but to improve pose matching
	public  boolean mov_en =                        true;  // enable detection/removal of the moving objects during pose matching
	public  double  mov_sigma =                     1.5;   // pix - weighted-blur offsets before detection
	// next two to prevent motion detection while errors are too big
	public  double  mov_max_std =                   0.5;   // pix
	public  double  mov_thresh_rel =                3.5;   // .0;   // exceed average error
	public  double  mov_thresh_abs=                 0.5;   // sqrt(dx^2+dy^2) in moving areas 
	public  double  mov_clust_max =                 1.5;   // cluster maximum should exceed threshold this times
	public  int     mov_grow =                      4;     // grow detected moving area
	public  boolean mov_show =                      true; // show debug images for movement detection
	public  int     mov_debug_level =               1;     // >0 verbose
	
	
	//LMA parameters
	
	public  boolean [] adjust_atr = new boolean [] {true,true,true};
	public  boolean [] adjust_xyz = new boolean [] {true,true,true};
	public  double  exit_change_atr =               1.0E-5;//rad,  L2 norm for difference to ext LMA 
	public  double  exit_change_xyz =               1.0E-3;//meters, L2 norm for difference to ext LMA 
	public  int     max_cycles =                   10;     // hard limit on full correlation/LMA cycles
	public  int     max_LMA =                      25;     // hard limit on LMA iterations
	public  double  max_rms =                       2.0;   // maximal RMS to consider adjustment to be a failure
	
	// Debug and visualization
	public  boolean scene_is_ref_test=              false; // correlate ref-2-ref for testing
	private boolean render_ref =                    true;  // render reference scene
	private boolean render_scene =                  true;  // render scene to be correlated with ref
	public  boolean toRGB =                         true;  // render scenes in pseudo-colors DOES NOT WORK - use ColorProcParameters
	private boolean show_2d_correlations =          true;  // show raw 2D correlations (individual and combined)
	private boolean show_motion_vectors =           true;  // show calculated motion vectors
	public  int     debug_level =                     -1;  // all renders are disable for debug_level < 0, scene "renders" for for debug_level < 1
    

	// Pairwise ERS testing
	public boolean  test_ers =                     false;
	public  int     test_ers0 =                       -1; // try adjusting a pair of scenes with ERS. Reference scene index
	public  int     test_ers1 =                       -1; // try adjusting a pair of scenes with ERS. Other scene index
	
	
	public boolean renderRef()            {return renderRef             (debug_level);}
	public boolean renderScene()          {return renderScene           (debug_level);}
	public boolean show2dCorrelations()   {return show2dCorrelations    (debug_level);}
	public boolean showMotionVectors()    {return showMotionVectors     (debug_level);}
	public boolean showCorrMotion()       {return showCorrMotion        (debug_level);}
	public boolean showMovementDetection(){return showMovementDetection (debug_level);}
	public int movDebugLevel()            {return movDebugLevel         (debug_level);}

	public boolean renderRef(int dlev)            {return (dlev>1) && render_ref;}
	public boolean renderScene(int dlev)          {return (dlev>1) && render_scene;}
	public boolean show2dCorrelations(int dlev)   {return (dlev>1) && show_2d_correlations;}
	public boolean showMotionVectors(int dlev)    {return (dlev>1) && show_motion_vectors;}
	public boolean showCorrMotion(int dlev)       {return (dlev>0) && show_motion_vectors;}
	public boolean showMovementDetection(int dlev){return (dlev>0) && mov_show;}
	public int movDebugLevel(int dlev)            {return (dlev > -1) ? mov_debug_level : 0;}
	
	
	public IntersceneMatchParameters() {
		
	}

	public void dialogQuestions(GenericJTabbedDialog gd) {
		//		gd.addMessage  ("Scene parameters selection");

		gd.addMessage  ("Build series options");
		gd.addCheckbox ("Force reference scene DSI calculation",     this.force_ref_dsi,
				"Calculate reference scene DSI even if the file exists.");
		gd.addCheckbox ("Force egomotion calculation",               this.force_orientations,
				"Calculate relative poses of each scene camera relative to the reference scene even if the data exists.");
		gd.addNumericField("Minimal number of egomotion calculations",this.min_num_orient, 0,3,"",
				"Minimal number of fitting scenes cycles, should be >=1. First cycle includes spiral search for the first scene");
		gd.addNumericField("Minimal number of interscene accumulations", this.min_num_interscene, 0,3,"",
				"Minimal required number of re-calculations of the interscene-accumulated DSI.");
		
		gd.addMessage  ("Generate/show scene sequences");
		gd.addCheckbox ("Generate mapped scene sequence",            this.generate_mapped,
				"Generate scene sequence mapped to the reference scene.");
		gd.addCheckbox ("Reuse existing video fragments, disable other options",this.reuse_video,
				"Dry-run video generation, just re-assemble earlier generated video fragments. Disables all other file generation.");
		gd.addCheckbox ("Save scene sequences in (pseudo)colors as TIFF",this.save_mapped_color,
				"Save generated scene sequences in (pseudo)color mode as a multi-slice TIFF.");
		gd.addCheckbox ("Save scene sequences in monochrome as TIFF",    this.save_mapped_mono,
				"Show generated scene sequences in monochrome mode as a multi-slice TIFF. May use Unsharp Mask.");
		gd.addCheckbox ("Generate color video",                      this.gen_avi_color,
				"Generate video for color scene sequences.");
		gd.addCheckbox ("Generate monochrome video",                 this.gen_avi_mono,
				"Generate video for monochrome scene sequences.");
		gd.addCheckbox ("Show scene sequences in (pseudo) colors",    this.show_mapped_color,
				"Show generated scene sequences in (pseudo) color mode. Disabled in batch mode");
		gd.addCheckbox ("Show scene sequences in monochrome",        this.show_mapped_mono,
				"Show generated scene sequences in monochrome mode. May use Unsharp Mask. Disabled in batch mode.");
		
		
		gd.addCheckbox ("Generate RAW images",                       this.generate_raw,
				"Raw images from single (top) camera using original view - just after aberration correction and aligning sensors .");
		gd.addCheckbox ("Generate INF images",                       this.generate_inf,
				"Single-camera images aligned at infinity - 2D correction functionally similar to camera stabilization.");
		gd.addCheckbox ("Generate 3D-corrected FG images",           this.generate_fg,
				"Correct in 3D scene images (from all 16 sensors) matching reference (last in sequence) scene with foreground (FG) priority.");
		gd.addCheckbox ("Generate 3D-corrected BG images",           this.generate_bg,
				"Correct in 3D scene images (from all 16 sensors) matching reference (last in sequence) scene with background (BG) priority.");
		gd.addCheckbox ("Generate binocular stereo pairs",           this.generate_stereo,
				"Generate stereo-pairs for 3D-corrected videos (FG,BG). Ebables specific modes (including 0-baseline / mono).");
		/*
		for (int i = 0; i < stereo_bases.length; i++) {
			double base = stereo_bases[i];
			String title = (base == 0.0)?
					"Generate mono (single camera) scene sequences":
						"Generate "+base+"mm-baseline stereo scene sequences";
			String tooltip = (base == 0.0)?
					"Generate mono (single camera) scene sequences as Tiff and/or video.":
						"Generate "+base+"mm-baseline stereo scene sequences as Tiff and/or video.";
			gd.addCheckbox (title,                                   this.generate_stereo_var[i], tooltip);
		}
		*/

		for (int i = 0; i < stereo_views.length; i++) {
//			String stereo_view = doublesToString(stereo_views[i]);
			double base = stereo_views[i][0]; // stereo_bases[i];
			String ub = String.format("(%.0fmm up, %.0fmm back) ",stereo_views[i][1],stereo_views[i][2]);
			if ((stereo_views[i][1]==0) && (stereo_views[i][2]==0)){
				ub="";
			}
			String title = (base == 0.0)?
					"Generate mono (single camera) scene sequences"+ub:
						"Generate "+base+"mm-baseline stereo scene sequences"+ub;
			String tooltip = (base == 0.0)?
					"Generate mono (single camera) scene sequences "+ub+"as Tiff and/or video.":
						"Generate "+base+"mm-baseline stereo scene sequences "+ub+"as Tiff and/or video.";
			gd.addCheckbox (title,                                   this.generate_stereo_var[i], tooltip);
		}
		
		
		
		
		gd.addMessage  ("Generate/save reference (last) scene images");
		gd.addCheckbox ("Export all-sensor images",                  this.export_images,
				"Export multi-slice images: with constant disparity, with foreground disparity, and with background disparity");
		gd.addCheckbox ("Show exported images (same disparity)",     this.show_images,
				"Display generated/saved image set, pseudocolors");
		gd.addCheckbox ("Show exported FG/BG 3d-corrected",          this.show_images_bgfg,
				"Show foreground and background exported images");
		gd.addCheckbox ("Show floating-point monochrome images",     this.show_images_mono,
				"Display generated/saved monochrome images (in addition to color) ");
		gd.addCheckbox ("Generate metric depth map",                 this.export_ranges,
				"Calculate strength, distance, X, and Y in meters");
		gd.addCheckbox ("Show metric depth map",                     this.show_ranges,
				"Display calculated depth map in meters");
		gd.addCheckbox ("Export files for DNN training/testing",     this.export_ml_files,
				"Display calculated depth map in meters");
		
		gd.addMessage  ("Additional parameters");

		gd.addCheckbox ("Color NaN background",                      this.show_color_nan,
				"Use NaN for undefined tiles (false - 0.0f). NaN produces sharp distinct result, 0.0f - blended");
		gd.addCheckbox ("Mono NaN background",                       this.show_mono_nan,
				"Use NaN for undefined tiles (false - 0.0f). NaN produces sharp distinct result, 0.0f - blended");
		gd.addNumericField("Minimal number of scenes to keep series",this.min_num_scenes, 0,3,"",
				"Scrap all seriest if less numer of scenes can be matched to the reference scene (including reference itself)");
		gd.addNumericField("LPF egomotion sigma",                    this.blur_egomotion, 3,5,"scenes",
				"LPF egomotion components with this sigma before using as ERS (not implemented).");

		gd.addMessage  ("Metric distance map generation");
		gd.addNumericField("Disparity at infinity",                  this.range_disparity_offset, 5,7,"pix",
				"Disparity at infinity - subtract from measured disparity when converting to ranges.");
		gd.addNumericField("Minimal strength for range calculation", this.range_min_strength, 5,7,"",
				"Disregard weaker results when measuring range.");
		gd.addNumericField("Maximal displayed range",                this.range_max, 5,7,"m",
				"Do not display extremely far objects.");
		

		gd.addMessage  ("Depth map filtering parameters");
		gd.addNumericField("Average lowest disparity neighbors",     this.num_bottom, 0,3,"",
				"Average this number of lowest disparity neighbors (of 8)");
		gd.addNumericField("Number of filter iterations",            this.num_passes, 0,3,"",
				"Number of filter iterations");
		gd.addNumericField("Maximal change to exit",                 this.max_change, 5,7,"pix",
				"Maximal change to exit iterations");
		gd.addNumericField("Minimal disparity",                      this.min_disparity, 5,7,"pix",
				"Minimal acceptable disparity");
		gd.addNumericField("Disparity range for symmetrical pull",   this.max_sym_disparity, 5,7,"pix",
				"For larger disparities undefined tiles are pull predominantly down (to average of lowest), "+
		        "but for small disparities (sky) pull is to average of all neighbors.");
		gd.addNumericField("Minimal strength LMA",                   this.min_strength_lma, 5,7,"",
				"Lower strength for LMA tiles - treat as non-LMA");
		gd.addNumericField("Minimal strength to replace",            this.min_strength_replace, 5,7,"",
				"Minimal strength to replace (now not used)");
		gd.addNumericField("Minimal strength to blur",               this.min_strength_blur, 5,7,"",
				"Minimal strength to blur  (now not used)");
		gd.addNumericField("Blur weak sigma",                        this.sigma, 5,7,"",
				"Blur weak sigma");
		gd.addNumericField("Number of blur passes",                  this.num_blur, 0,3,"",
				"Number of blur passes");
		gd.addNumericField("Disparity offset (not used)",            this.disparity_corr, 5,7,"",
				"for setInterTasks() - currently not used ");
		gd.addNumericField("Outlier N-th from extreme",              this.outliers_nth_fromextrem, 0,3,"",
				"0 - use min/max, 1 - use second min/max, ... second from min/max - removes dual-tile clusters max/mins");
		gd.addNumericField("Outliers tolerance absolute",            this.outliers_tolerance_absolute, 5,7,"pix",
				"Absolute difference from extreme neighbor");
		gd.addNumericField("Outliers tolerance realtive",            this.outliers_tolerance_relative, 5,7,"pix/pix",
				"Add to tolerance as a fraction of absolute tile disparity");
		gd.addNumericField("Maximal number of outlier iterations",   this.outliers_max_iter, 0,3,"",
				"Maximal number of iterations for removing outliers");
		gd.addNumericField("Maximal outlier strength 2",             this.outliers_max_strength2, 5,7,"",
				"Maximal outlier strength for second pass filtering");
		gd.addNumericField("Outlier N-th from extreme 2",            this.outliers_nth_fromextrem2, 0,3,"",
				"Second filter: 0 - use min/max, 1 - use second min/max, ... second from min/max - removes dual-tile clusters max/mins");
		gd.addNumericField("Outliers tolerance absolute 2",          this.outliers_tolerance_absolute2, 5,7,"",
				"Absolute difference from extreme neighbor, second filter");
		gd.addNumericField("Outliers tolerance realtive 2",          this.outliers_tolerance_relative2, 5,7,"",
				"Add to tolerance as a fraction of absolute tile disparity, second filter");
		gd.addNumericField("Outliers LMA maximal strength",          this.outliers_lma_max_strength, 5,7,"",
				"Maximal LMA strength for an outlier");
		gd.addNumericField("Outliers non-LMA maximal strength",      this.outliers_max_strength, 5,7,"",
				"Outliers non-LMA (centroid) maximal strength");
		gd.addNumericField("Outliers from LMA max strength",         this.outliers_from_lma_max_strength, 5,7,"",
				"Non-LMA outliers from LMA neighbors maximal strength");
		
		gd.addNumericField("Search for LMA radius",                  this.search_radius, 0,3,"tile",
				"Search farther if no LMA tiles are found around this non-LMA tile");
		gd.addCheckbox ("Remove non-LMA tile if no LMA near",        this.remove_no_lma_neib,
				"Remove non-LMA tile if no LMA one is found within specified radius");
		
		gd.addNumericField("LMA difference positive",                this.diff_from_lma_pos, 5,7,"",
				"Maximal non-LMA tile positive difference from the LMA neighbor");
		gd.addNumericField("LMA difference negative",                this.diff_from_lma_neg, 5,7,"",
				"Maximal non-LMA tile negative difference from the LMA neighbor");
		gd.addNumericField("LMA outliers N-th from extreme",         this.outliers_lma_nth_fromextrem, 0,3,"",
				"0 - use min/max, 1 - use second min/max, ... ");
		gd.addNumericField("Margin from the sensor FOV edge ",       this.filter_margin, 0,3,"pix",
				"Disregard tiles with centers closer to the sensor FoV in pixels");
		
		gd.addMessage  ("Removing weak tiles with insufficient matching neighbors");
		gd.addNumericField("Neighbor disparity tolerance (absolute)",this.weak_tolerance_absolute, 5,7,"pix",
				"Absolute disparity difference to be considered a neighbor");
		gd.addNumericField("Neighbor disparity tolerance (relative)", this.weak_tolerance_relative, 5,7,"",
				"Relative disparity difference to be considered a neighbor");
		gd.addNumericField("Minimal matching neighbors",              this.weak_min_neibs, 0,3,"",
				"Minimal number of weak matching neighbors (at least one is needed if there are strong non-neighbors).");
		gd.addNumericField("Minimal strength to be strong",           this.strong_strength, 5,7,"",
				"Minimal strength to be considered strong(not weak). LMA are also considered strong.");
		gd.addNumericField("Minimal weak strength",                   this.weak_strength, 5,7,"",
				"Weaker are removed unconditionally (not used now).");
		
		gd.addMessage  ("Stereo");
		/*
		if (stereo_bases.length > 0) {
			String [] stereo_choices = new String [stereo_bases.length + 1];
			stereo_choices[0] = "--none--";
			for (int i = 0; i < stereo_bases.length; i++) {
				stereo_choices[i+1] = stereo_bases[i]+"mm";	
			}
			gd. addChoice("Remove stereo-base",                 stereo_choices, stereo_choices[0], 
					"Remove selected stereo-base");
		}
		*/
		if (stereo_views.length > 0) {
			String [] stereo_choices = new String [stereo_views.length + 1];
			stereo_choices[0] = "--none--";
			for (int i = 0; i < stereo_views.length; i++) {
				stereo_choices[i+1] = doublesToString(stereo_views[i],"%.0f")+" mm";	
			}
			gd. addChoice("Remove stereo-view (base, up, back)",      stereo_choices, stereo_choices[0], 
					"Remove selected stereo-view, consisting of streo-base, viewpoint above camera, viewpoint behing camera - all in mm");
		}

		
		gd.addStringField("Add another stereo view (baseline, above, behind)",                 "", 40,
				"Add another stereo view by providing baseline, above camera, behind camera (mm).");
		
//		gd.addNumericField("Stereo baseline",                        this.stereo_baseline, 5,7,"mm",
//				"Synthetic 3D with possibly exagerrated stereo baseline");
		gd.addCheckbox ("Stereo merge right, left",                  this.stereo_merge,
				"Combine stereo pair in a single (wide) frame. Unchecked - generate separate videos.");
		gd.addNumericField("Stereo gap",                             this.stereo_gap, 0,3,"pix",
				"Distance (pixels) between right and left sub-images in a stereo frame.");
		gd.addNumericField("VR headset stereo base",                 this.stereo_intereye, 5,7,"mm",
				"Inter-lens distance of the VR phone adapter (Google Cardboard - 63.5)");
		gd.addNumericField("Phone screen width",                     this.stereo_phone_width, 5,7,"mm",
				"Phone screen width in full screen video mode.");
		
		
		gd.addNumericField("Scene sequence horizontal extra",        this.extra_hor_tile, 0,3,"tiles",
				"Enlarge reference scene window horizontally in each direction to accommodate other scenes in a sequence");
		gd.addNumericField("Scene sequence vertical extra",          this.extra_vert_tile, 0,3,"tiles",
				"Enlarge reference scene window vertically in each direction to accommodate other scenes in a sequence");
		gd.addCheckbox ("Crop 3D",                                   this.crop_3d,
				"Do not enlarge reference scene windows fro 3D views (FG, BG)");
		
		gd.addNumericField("Sensor mask (bitmask, -1 - all sensors)",this.sensor_mask, 0,3,"",
				"Select which sensors to be included in each scene of the sequence");
		gd.addCheckbox ("Merge all channels in 3D modes",                                   this.merge_all,
				"Ignore sensor mask, use all channels and merge them into one in 3D modes (FG and BG)");
		
///		gd. addChoice("3D mode",                 MODES3D, MODES3D[this.mode3d + 1], 
///				"3D mode for rendering scenes in a sequence: RAW - raw images, INF - no 3D, use infinity; FG - Foreground; BG - Background");
		
		
		gd.addMessage  ("Generate video for scene sequences");
		gd.addNumericField("Video output frame rate",                this.video_fps, 5,7,"fps",
				"Video frames per second (original is 60fps).");
		gd.addNumericField("Sensor frame rate (60 for Boson)",       this.sensor_fps, 5,7,"fps",
				"Sensor (real) frame rate to be used in realtime videos after slowmo.");
		gd. addChoice("AVI mode",                                     MODES_AVI, MODES_AVI[this.mode_avi], 
				"AVI generation mode");
		gd.addNumericField("AVI JPEG Quality",                       this.avi_JPEG_quality, 0,3,"",
				"Only used in JPEG compression mode");		
		
		gd.addCheckbox ("Convert AVI to WEBM",                       this.run_ffmpeg,
				"Run ffmpeg on AVI video to convert to smaller modern video.");
		gd.addStringField ("Result video extension (e.g. \".webm\")",this.video_ext, 60,
				"Converted video extension, starting with dot.");
		gd.addStringField ("Video encoder",                          this.video_codec, 60,
				"FFMPEG video encoder, such as \"VP8\" or \"VP9\".");
		gd.addNumericField("Video CRF",                              this.video_crf, 0,3,"",
				"Quality - the lower the better. 40 - OK");		

		gd.addCheckbox ("Remove AVI",                                this.remove_avi,
				"Remove large AVI files after (and only) conversion with ffmpeg.");
		
		gd.addStringField ("Video encoder for combining",            this.video_codec_combo, 60,
				"FFMPEG video encoder, such as \"VP8\" or \"VP9\". Applies when merging segments.");
		gd.addNumericField("Video CRF for combining",                this.video_crf_combo, 0,3,"",
				"Quality - the lower the better. 40 - OK.  Applies when merging segments.");		
		
		gd.addCheckbox ("Apply unsharp mask to mono",                this.um_mono,
				"Apply unsharp mask to monochrome image sequences/video.  Applies to TIFF generatiojn too");
		gd.addNumericField("Unsharp mask sigma (radius)",            this.um_sigma, 5,7,"pix",
				"Unsharp mask Gaussian sigma.");
		gd.addNumericField("Unsharp mask weight",                    this.um_weight, 5,7,"",
				"Unsharp mask weightt (multiply blurred version before subtraction from the original).");
		
		gd.addCheckbox ("Fixed monochrome range",                    this.mono_fixed,
				"Normalize monochrome (after UM) to a fixed range when converting to 8 bit RGB.");
		gd.addNumericField("Monochrome full range",                  this.mono_range, 5,7,"",
				"Monochrome full range to convert to 0..255.");
		

		gd.addCheckbox ("Generate anaglyph stereo",                  this.anaglyth_en,
				"Apply unsharp mask to monochrome image sequences/video.  Applies to TIFF generatiojn too");

		{String scolor = String.format("%08x", getLongColor(this.anaglyph_left));
		gd.addStringField ("Anaglyph color left",scolor, 8, "Any invalid hex number sets default red");}
		
		{String scolor = String.format("%08x", getLongColor(this.anaglyph_right));
		gd.addStringField ("Anaglyph color right",scolor, 8, "Any invalid hex number sets default cyan");}
		
		
		
		
		
		gd.addCheckbox ("Timestamp color videos",                    this.annotate_color,
				"Annotate pseudo-color video frames with timestamps.");
		gd.addCheckbox ("Timestamp monochrome videos",               this.annotate_mono,
				"Annotate monochrome video frames with timestamps.");
		String scolor = (this.annotate_color_color==null)?"none":String.format("%08x", getLongColor(this.annotate_color_color));
		gd.addStringField ("Timestamp color for pseudocolor frames",scolor, 8, "Any invalid hex number disables annotation");
		scolor = (this.annotate_color_mono==null)?"none":String.format("%08x", getLongColor(this.annotate_color_mono));
		
		gd.addStringField ("Timestamp color for monochrome frames",scolor, 8, "Any invalid hex number disables annotation");
		gd.addCheckbox ("Transparent timestamp background (monochrome)",this.annotate_transparent_mono,
				"Put monochrome timestamp over image (unchecked - use black background). Color - always black.");

		
		gd.addMessage  ("Interscene match parameters");
		gd.addNumericField("Image margin",                           this.margin, 0,5,"pix",
				"Do not use tiles if their centers are closer to the virtual image edge");
		gd.addNumericField("Used sensors mask",                      this.sensor_mask_inter, 0,5,"",
				"Bitmask of the sensors to use, -1 - use all");
		gd.addCheckbox ("Preserve partial",                          this.use_partial,
				"Preserve motion vectors for each of the selected sensors, fallse - use only combined (in transform domain and in pixel domain)");
		gd.addCheckbox ("Poly maximums",                             this.run_poly,
				"Use polinomial interpolationn to determin correlation maximums instead of centroids (not yet implemented, ignored)");
		gd.addNumericField("Centroid radius",                        this.centroid_radius, 5,7,"pix",
				"Calculate centroids after multiplication by a half-cosine window. All correlation data farther than this value from the center is ignored");
		gd.addNumericField("Refine centroids",                       this.n_recenter, 0,5,"",
				"Repeat centroids after moving the window center to the new centroid location this many times (0 - calculate once)");

	//	gd.addMessage  ("Filter tiles to use for scene poses");
	//	gd.addNumericField("DSI_MAIN minimal strength",              this.min_ref_str, 5,7,"",
	//			"Match only tiles where DSI_MAIN is stronger than that (and has LMA).");
		
		
		gd.addMessage  ("Mixing TD and PD accumulation of 2d correlations");
		gd.addNumericField("TD-accumulated weight",                  this.td_weight, 5,7,"",
				"Mix argmax from TD-accumulated correlation.");
		gd.addNumericField("PD-accumulated weight",                  this.pd_weight, 5,7,"",
				"Mix argmax from PD-accumulated correlation.");
		gd.addCheckbox ("TD when no PD only",                        this.td_nopd_only,
				"Use argmax from TD only if PD data is not available for this tile.");
		
		gd.addMessage  ("Filtering motion vectors");
		gd.addNumericField("Minimal correlation strength (non-sum)",       this.min_str, 5,7,"",
				"Minimal correlation strength for individual correlation and for pixel-domain averaged one. Weeker tiles results are removed.");
		gd.addNumericField("Minimal correlation strength (non-sum) w/FPN", this.min_str_fpn, 5,7,"",
				"Similar to above, but for small offsets where FPN correlation may be present");
		gd.addNumericField("Minimal correlation strength (sum only)",      this.min_str_sum, 5,7,"",
				"Minimal correlation strength for transform-domain averaging. Weeker tiles results are removed.");
		gd.addNumericField("Minimal correlation strength (sum only) w/FPN",this.min_str_sum_fpn, 5,7,"",
				"Similar to above, but for small offsets where FPN correlation may be present");
		
		gd.addNumericField("Minimal number of neighbors (of 8)",     this.min_neibs, 0,3,"",
				"Remove motion vectors with less than this number of defined (passing min_str) neighbors.");
		gd.addNumericField("No-neighbors weight (<1.0)",             this.weight_zero_neibs, 5,7,"",
				"Punish for reduced neighbors - weigh for no-neighbors), weight of 8 neighbors = 1.0.");
		gd.addNumericField("Disparity to reduce weight twice from infinity", this.half_disparity, 5,7,"",
				"Weight at this disparity is 0.5, at infinity - 1.0.");
		gd.addNumericField("Difference from neighbors average ",     this.half_avg_diff, 5,7,"",
				"Reduce twice for high difference from neighbors average.");

		gd.addMessage  ("Photometric calibration (move elsewhere?)");
		gd.addCheckbox ("Enable photometric calibration",            this.photo_en,
				"Equalize per- sensor gains and offsets. Requires disparity map. Save to reference scene and with current scene (to .corr-zml).");
		gd.addNumericField("Full photometric (re)calibrations",      this.photo_num_full, 0,3,"pix",
				"Full recalibratrions include re-importing raw images with updated offsets/gains");
		gd.addNumericField("Refines",                                this.photo_num_refines, 0,3,"pix",
				"Calculate calibration, remove outliers (e.g. FG/BG) and repeat");
		gd.addNumericField("Minimal DSI strength",                   this.photo_min_strength, 5,7,"",
				"Do not use weak tiles.");
		gd.addNumericField("Maximal channel mismatch",               this.photo_max_diff, 5,7,"",
				"Detect (and remove outliers). Adjusted images have RMSE ~9 counts.");
		gd.addCheckbox ("Debug pphotometric calibration",            this.photo_debug,
				"Generate debug images an text output.");
		
		
		gd.addMessage  ("Initial search for the inter-scene match");
		gd.addNumericField("DSI_MAIN minimal strength",              this.min_ref_str, 5,7,"",
					"Match only tiles where DSI_MAIN is stronger than that (and has LMA).");
		
		gd.addNumericField("Azimuth/tilt step",                      this.pix_step, 0,3,"pix",
				"Search in a spiral starting with no-shift with this step between probes, in approximate pixels");
		gd.addNumericField("Search spiral radius",                   this.search_rad, 0,3,"steps",
				"Maximal search radius in steps");
		gd.addNumericField("\"Maybe\" sum of strengths",             this.maybe_sum, 5,7,"",
				"Minimal acceptable sum of defined tiles strengths (will look for the best among matching)");
		gd.addNumericField("\"Sure\" sum of strengths",              this.sure_sum, 5,7,"",
				"Definitely sufficient sum of defined tiles strengths (will not continue looking for better).");
		gd.addNumericField("\"Maybe\" average of strengths",         this.maybe_avg, 5,7,"",
				"Minimal acceptable average of defined tiles strengths (will look for the best among matching)");
		gd.addNumericField("\"Sure\" average of strengths",          this.sure_avg, 5,7,"",
				"Definitely sufficient average of defined tiles strengths (will not continue looking for better).");
		
		gd.addNumericField("Maximal offset standard deviation",      this.max_search_rms, 5,7,"pix",
				"Standard deviation for X,Y offsets.");
		gd.addNumericField("\"Maybe\" FOM",                          this.maybe_fom, 5,7,"",
				"Minimal acceptable Figure of Merit (semi-total strength divided by standard deviation of offsets), will look for the best among matching.");
		gd.addNumericField("\"Sure\" FOM",                           this.sure_fom, 5,7,"",
				"Definitely sufficient FOM (semi-total strength divided by standard deviation of offsets), will non continue looking for better.");
		
		
		gd.addMessage  ("Reference scene disparity");
		gd.addCheckbox ("Use combo DSI (if available)",              this.use_combo_dsi,
				"Use interscene DSI if available (instead of the single-scene)");
		gd.addCheckbox ("LMA tiles only",                            this.use_lma_dsi,
				"Use only strong correlation tiles (with LMA available) for interscene correlation (pose matching)");

		gd.addMessage  ("Supress FPN in interscene correlations");
		gd.addCheckbox ("Enable FPN suppression",                    this.fpn_remove,
				"Zero-out integrated inter-scene correlation around zero-shift offset");
		gd.addNumericField("Maximal FPN offset to consider",         this.fpn_max_offset, 6,7,"pix",
				"Maximal offset from the zero shift to consider (normally just 8 pix)");
		gd.addNumericField("FPN suppression radius",                 this.fpn_radius, 6,7,"pix",
				"Blank correlation pixels closer than this distance from the FPN offset");
		gd.addCheckbox ("Ignore maximums \"touching\" FPN",          this.fpn_ignore_border,
				"Discard TD integrated tiles where local maximum is a neighbor (including diagonal) to blanked FPN correlation pixels");
		
		gd.addMessage  ("Detect and remove moving objects from pose matching");
		gd.addCheckbox ("Enable movement detection/elimination",     this.mov_en,
				"Detect and mask areas with detected movement to improve pose matching");
		gd.addNumericField("Detection blur sigma",                   this.mov_sigma, 6,7,"pix",
				"Blur squared difference (dx^2+dy^2) befopre thresholding for movement detection");
		gd.addNumericField("Max weighted mismatch std",              this.mov_max_std, 6,7,"pix",
				"Do not try to detect moving objects until approximate pose match is achieved (standard deviation)");
		gd.addNumericField("Relative threshold",                     this.mov_thresh_rel, 6,7,"x",
				"Moving areas over standard deviation. Both relative and absolute should be exceeded.");
		gd.addNumericField("Absolute threshold",                     this.mov_thresh_abs, 6,7,"pix",
				"Moving areas sqrt(dx*dx+dy*dy). Both relative and absolute should be exceeded.");
		gd.addNumericField("Cluster max over threshold",             this.mov_clust_max, 6,7,"",
				"Moving cluster should contain tile with this exceed over thresholds");
		gd.addNumericField("Moving cluster grow",                    this.mov_grow, 0,3,"",
				"Standard grow values - 1 - ortho, 2 - diagonal, 3 - twice orto, 4 - twice diagonal");
		gd.addCheckbox ("Show movement debug images",                this.mov_show,
				"Disabled if 'Debug Level for interscene match' < 1");
		gd.addNumericField("Debug level for movement detection (0/1)", this.mov_debug_level, 0,3,"",
				"Disabled if 'Debug Level for interscene match' < 0");
		
		gd.addMessage  ("LMA parameters");
		gd.addCheckbox ("Azimuth",                                   this.adjust_atr[0],
				"Adjust scene camera azimuth with LMA");
		gd.addCheckbox ("Tilt",                                      this.adjust_atr[1],
				"Adjust scene camera tilt with LMA");
		gd.addCheckbox ("Roll",                                      this.adjust_atr[2],
				"Adjust scene camera roll with LMA");
		gd.addCheckbox ("X",                                         this.adjust_xyz[0],
				"Adjust scene camera X with LMA");
		gd.addCheckbox ("Y",                                         this.adjust_xyz[1],
				"Adjust scene camera Y with LMA");
		gd.addCheckbox ("Z",                                         this.adjust_xyz[2],
				"Adjust scene camera Z with LMA");

		gd.addNumericField("Exit ATR change",                        this.exit_change_atr, 6,7,"pix",
				"L2 of azimuth, tilt, roll change (in pixels at infinity) to exit fitting");
		gd.addNumericField("Exit XYZ change",                        this.exit_change_xyz, 5,7,"m",
				"L2 of linear scene camera position change (in meters) to exit fitting");


		gd.addNumericField("Max full iterations",                    this.max_cycles, 0,3,"",
				"Hard limit on interscene correlations plus LMA cycles");
		gd.addNumericField("Max LMA iterations",                     this.max_LMA, 0,3,"",
				"Hard limit on LMA iterations");
		gd.addNumericField("Maximal RMS to fail",                    this.max_rms, 5,7,"pix",
				"Maximal fitting RMS to consider fitting failure");


		gd.addMessage  ("Debug and visialization parameters");
		gd.addCheckbox ("Replace scene with reference scene",        this.scene_is_ref_test,
				"Correlate reference scene with itself for testing (may want to manually change scene_atr and scene_xyz in debug mode)");
		gd.addCheckbox    ("Render reference scene",                 this.render_ref,
				"Render reference scene as a multi-layer (per-port) image");
		gd.addCheckbox    ("Render compared scene",                  this.render_scene,
				"Render reference compared as a multi-layer (per-port) image");
		gd.addCheckbox    ("Pseudo-color render (does not work)",    this.toRGB,
				"Use pseudo-color palette for scene render  - DOES NOT WORK - use ColorProcParameters");
		gd.addCheckbox    ("Show 2D correlations",                   this.show_2d_correlations,
				"Show raw 2D correlations (per channel and combined)");
		gd.addCheckbox    ("Show motion vectors",                    this.show_motion_vectors,
				"Show motion vectors and strengths per channel and combines, NaN-s are when strengths corresponding are below thresholds");
		gd.addNumericField("Debug Level for interscene match",       this.debug_level, 0,3,"",
				"Debug Level for the above parameters: renders and raw correlations need >1,  motion vectors > 0");

		gd.addMessage  ("Pairwise ERS testing");
		gd.addCheckbox ("Test pairwise matching with ERS",        this.test_ers,
				"Test pairwise dispaly and pose correction to test ERS compensation");
		gd.addNumericField("Test scene reference index",             this.test_ers0, 0,3,"",
				"Reference scene index in a scene sequence");
		gd.addNumericField("Test scene other scene index",           this.test_ers1, 0,3,"",
				"Other scene index in a scene sequence (should have a very different angular/linear velocity component)");
		
		
	}

	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.force_ref_dsi =                  gd.getNextBoolean();
		this.force_orientations =             gd.getNextBoolean();
		this.min_num_orient =           (int) gd.getNextNumber(); if (min_num_orient < 1) min_num_orient = 1;
		this.min_num_interscene =       (int) gd.getNextNumber(); if (min_num_interscene < 1) min_num_interscene = 1;
		this.generate_mapped =                gd.getNextBoolean();
		this.reuse_video =                    gd.getNextBoolean();
		this.save_mapped_color =              gd.getNextBoolean(); // as tiff
		this.save_mapped_mono =               gd.getNextBoolean(); // as tiff
		this.gen_avi_color =                  gd.getNextBoolean(); // as video
		this.gen_avi_mono =                   gd.getNextBoolean(); // as video
		this.show_mapped_color =              gd.getNextBoolean();
		this.show_mapped_mono =               gd.getNextBoolean();
		this.generate_raw =                   gd.getNextBoolean();
		this.generate_inf =                   gd.getNextBoolean();
		this.generate_fg =                    gd.getNextBoolean();
		this.generate_bg =                    gd.getNextBoolean();
		this.generate_stereo =                gd.getNextBoolean();
		/*
		for (int i = 0; i < stereo_bases.length; i++) {
			this.generate_stereo_var[i] =     gd.getNextBoolean();
		}
		*/
		for (int i = 0; i < stereo_views.length; i++) {
			this.generate_stereo_var[i] =     gd.getNextBoolean();
		}
		
		this.export_images =                  gd.getNextBoolean();
		this.show_images =                    gd.getNextBoolean();
		this.show_images_bgfg =               gd.getNextBoolean();
		this.show_images_mono =               gd.getNextBoolean();
		this.export_ranges =                  gd.getNextBoolean();
		this.show_ranges =                    gd.getNextBoolean();
		this.export_ml_files =                gd.getNextBoolean();

		// additional parameters		
		this.show_color_nan =                 gd.getNextBoolean();
		this.show_mono_nan =                  gd.getNextBoolean();
		
		
		this.min_num_scenes =           (int) gd.getNextNumber();
		this.blur_egomotion =                 gd.getNextNumber();
		this.range_disparity_offset =         gd.getNextNumber();
		this.range_min_strength =             gd.getNextNumber();
		this.range_max =                      gd.getNextNumber();

		this.num_bottom =               (int) gd.getNextNumber();
		this.num_passes =               (int) gd.getNextNumber();
		this.max_change =                     gd.getNextNumber();
		this.min_disparity =                  gd.getNextNumber();
		this.max_sym_disparity =              gd.getNextNumber();
		this.min_strength_lma =               gd.getNextNumber();
		this.min_strength_replace =           gd.getNextNumber();
		this.min_strength_blur =              gd.getNextNumber();
		this.sigma =                          gd.getNextNumber();
		this.num_blur =                 (int) gd.getNextNumber();
		this.disparity_corr =                 gd.getNextNumber();
		this.outliers_nth_fromextrem =  (int) gd.getNextNumber();
		this.outliers_tolerance_absolute =    gd.getNextNumber();
		this.outliers_tolerance_relative =    gd.getNextNumber();
		this.outliers_max_iter =        (int) gd.getNextNumber();
		this.outliers_max_strength2 =         gd.getNextNumber();
		this.outliers_nth_fromextrem2 = (int) gd.getNextNumber();
		this.outliers_tolerance_absolute2 =   gd.getNextNumber();
		this.outliers_tolerance_relative2 =   gd.getNextNumber();
		this.outliers_lma_max_strength =      gd.getNextNumber();
		this.outliers_max_strength =          gd.getNextNumber();
		this.outliers_from_lma_max_strength = gd.getNextNumber();
		this.search_radius =            (int) gd.getNextNumber();
		this.remove_no_lma_neib =             gd.getNextBoolean();
		this.diff_from_lma_pos =              gd.getNextNumber();
		this.diff_from_lma_neg =              gd.getNextNumber();
		this.outliers_lma_nth_fromextrem=(int)gd.getNextNumber();
		this.filter_margin =            (int) gd.getNextNumber();
		
		this.weak_tolerance_absolute =        gd.getNextNumber();
		this.weak_tolerance_relative =        gd.getNextNumber();
		this.weak_min_neibs =           (int) gd.getNextNumber();
		this.strong_strength =                gd.getNextNumber();
		this.weak_strength =                  gd.getNextNumber();

		/*
		if (stereo_bases.length > 0) {
			int i =                           gd.getNextChoiceIndex();
			if (i > 0) {
				removeStereo(i-1);
			}
		}
		String s =                            gd.getNextString();
		if (s.trim().length() > 0) {
			try {
				double d = Double.parseDouble(s);
				addStereo(d, true);
//				orderStereo();
			} catch (Exception e) {

			}
		}
		 */
		if (stereo_views.length > 0) {
			int i =                           gd.getNextChoiceIndex();
			if (i > 0) {
				removeStereoView(i-1);
			}
		}
		String s =                            gd.getNextString();
		addStereoView(s, true);
		
//		this.stereo_baseline =          gd.getNextNumber();
		this.stereo_merge =             gd.getNextBoolean();
		this.stereo_gap =         (int) gd.getNextNumber();
		this.stereo_intereye =          gd.getNextNumber();
		this.stereo_phone_width =       gd.getNextNumber();
		
		this.extra_hor_tile =     (int) gd.getNextNumber();
		this.extra_vert_tile =    (int) gd.getNextNumber();
		this.crop_3d =                  gd.getNextBoolean();
		this.sensor_mask =        (int) gd.getNextNumber();
		this.merge_all =                gd.getNextBoolean();
		
//		this.mode3d =                   gd.getNextChoiceIndex() - 1;
		
		this.video_fps =                gd.getNextNumber();
		this.sensor_fps =               gd.getNextNumber();
		this.mode_avi =                 gd.getNextChoiceIndex();
		this.avi_JPEG_quality =   (int) gd.getNextNumber();
		this.run_ffmpeg =               gd.getNextBoolean();
		this.video_ext=                 gd.getNextString();	
		this.video_codec=               gd.getNextString();	
		this.video_crf =          (int) gd.getNextNumber();
		this.remove_avi =               gd.getNextBoolean();
		this.video_codec_combo=         gd.getNextString();	
		this.video_crf_combo =    (int) gd.getNextNumber();
		this.um_mono =                  gd.getNextBoolean();
		this.um_sigma =                 gd.getNextNumber();
		this.um_weight =                gd.getNextNumber();
		this.mono_fixed =               gd.getNextBoolean();
		this.mono_range =               gd.getNextNumber();
		
		this.anaglyth_en =              gd.getNextBoolean();
		{
			String scolor =             gd.getNextString();
			long lcolor = -1;
			try {
				lcolor = Long.parseLong(scolor,16);
				this.anaglyph_left = setLongColor(lcolor);
			} catch(NumberFormatException e){
				this.anaglyph_left = anaglyph_left_default;
			}
		}

		{
			String scolor =             gd.getNextString();
			long lcolor = -1;
			try {
				lcolor = Long.parseLong(scolor,16);
				this.anaglyph_right = setLongColor(lcolor);
			} catch(NumberFormatException e){
				this.anaglyph_right = anaglyph_right_default;
			}
		}
		
		this.annotate_color =           gd.getNextBoolean();
		this.annotate_mono =            gd.getNextBoolean();
		{
			String scolor =             gd.getNextString();
			long lcolor = -1;
			try {
				lcolor = Long.parseLong(scolor,16);
				this.annotate_color_color = setLongColor(lcolor);
			} catch(NumberFormatException e){
				this.annotate_color_color = null;
			}
		}
		{
			String scolor =             gd.getNextString();
			long lcolor = -1;
			try {
				lcolor = Long.parseLong(scolor,16);
				this.annotate_color_mono = setLongColor(lcolor);
			} catch(NumberFormatException e){
				this.annotate_color_mono = null;
			}
		}
		this.annotate_transparent_mono= gd.getNextBoolean();
		
		this.margin =             (int) gd.getNextNumber();
		this.sensor_mask_inter=   (int) gd.getNextNumber();
		this.use_partial =              gd.getNextBoolean();
		this.run_poly =                 gd.getNextBoolean();
		this.centroid_radius =          gd.getNextNumber();
		this.n_recenter =         (int) gd.getNextNumber();

//		this.min_ref_str =              gd.getNextNumber();

		this.td_weight =                gd.getNextNumber();
		this.pd_weight =                gd.getNextNumber();
		this.td_nopd_only =             gd.getNextBoolean();
		
		this.min_str =                  gd.getNextNumber();
		this.min_str_fpn =              gd.getNextNumber();
		this.min_str_sum =              gd.getNextNumber();
		this.min_str_sum_fpn =          gd.getNextNumber();

		this.min_neibs =          (int) gd.getNextNumber();
		this.weight_zero_neibs =        gd.getNextNumber();
		this.half_disparity =           gd.getNextNumber();
		this.half_avg_diff =            gd.getNextNumber();
		
		this.photo_en =                 gd.getNextBoolean();
		this.photo_num_full =     (int) gd.getNextNumber();
		this.photo_num_refines =  (int) gd.getNextNumber();
		this.photo_min_strength =       gd.getNextNumber();
		this.photo_max_diff =           gd.getNextNumber();
		this.photo_debug =              gd.getNextBoolean();
		
		this.min_ref_str =              gd.getNextNumber();
		this.pix_step =           (int) gd.getNextNumber();
		this.search_rad =         (int) gd.getNextNumber();
		this.maybe_sum =                gd.getNextNumber();
		this.sure_sum =                 gd.getNextNumber();
		this.maybe_avg =                gd.getNextNumber();
		this.sure_avg =                 gd.getNextNumber();
		
		this.max_search_rms =           gd.getNextNumber();
		this.maybe_fom =                gd.getNextNumber();
		this.sure_fom =                 gd.getNextNumber();

		this.use_combo_dsi =            gd.getNextBoolean();
		this.use_lma_dsi =              gd.getNextBoolean();
		
		this.fpn_remove =               gd.getNextBoolean();
		this.fpn_max_offset =           gd.getNextNumber();
		this.fpn_radius =               gd.getNextNumber();
		this.fpn_ignore_border =        gd.getNextBoolean();
		
		this.mov_en =                   gd.getNextBoolean();
		this.mov_sigma =                gd.getNextNumber();
		this.mov_max_std =              gd.getNextNumber();
		this.mov_thresh_rel =           gd.getNextNumber();
		this.mov_thresh_abs =           gd.getNextNumber();
		this.mov_clust_max =            gd.getNextNumber();
		this.mov_grow =           (int) gd.getNextNumber();
		this.mov_show =                 gd.getNextBoolean();
		this.mov_debug_level =    (int) gd.getNextNumber();
		
		this.adjust_atr[0] =            gd.getNextBoolean();
		this.adjust_atr[1] =            gd.getNextBoolean();
		this.adjust_atr[2] =            gd.getNextBoolean();
		this.adjust_xyz[0] =            gd.getNextBoolean();
		this.adjust_xyz[1] =            gd.getNextBoolean();
		this.adjust_xyz[2] =            gd.getNextBoolean();
		this.exit_change_atr =          gd.getNextNumber();
		this.exit_change_xyz =          gd.getNextNumber();
		this.max_cycles =         (int) gd.getNextNumber();
		this.max_LMA =            (int) gd.getNextNumber();
		this.max_rms =                  gd.getNextNumber();
		this.scene_is_ref_test =        gd.getNextBoolean();
		this.render_ref =               gd.getNextBoolean();
		this.render_scene =             gd.getNextBoolean();
		this.toRGB =                    gd.getNextBoolean();
		this.show_2d_correlations =     gd.getNextBoolean();
		this.show_motion_vectors =      gd.getNextBoolean();
		this.debug_level =        (int) gd.getNextNumber();

		this.test_ers =                 gd.getNextBoolean();
		this.test_ers0 =          (int) gd.getNextNumber();
		this.test_ers1 =          (int) gd.getNextNumber();
		
		if (this.weight_zero_neibs > 1.0) this.weight_zero_neibs = 1.0;
	}
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"force_ref_dsi",                 this.force_ref_dsi + "");               // boolean
		properties.setProperty(prefix+"force_orientations",            this.force_orientations + "");          // boolean
//		properties.setProperty(prefix+"readjust_orient",               this.readjust_orient + "");             // boolean
//		properties.setProperty(prefix+"force_interscene",              this.force_interscene + "");            // boolean
		properties.setProperty(prefix+"min_num_orient",                this.min_num_orient+"");                // int
		properties.setProperty(prefix+"min_num_interscene",            this.min_num_interscene+"");            // int
		properties.setProperty(prefix+"generate_mapped",               this.generate_mapped+"");               // boolean
		properties.setProperty(prefix+"reuse_video",                   this.reuse_video+"");                   // boolean
		properties.setProperty(prefix+"save_mapped_color",             this.save_mapped_color+"");             // boolean
		properties.setProperty(prefix+"save_mapped_mono",              this.save_mapped_mono+"");              // boolean
		properties.setProperty(prefix+"gen_avi_color",                 this.gen_avi_color+"");       // boolean
		properties.setProperty(prefix+"gen_avi_mono",                  this.gen_avi_mono+"");        // boolean
		properties.setProperty(prefix+"show_mapped_color",    this.show_mapped_color+"");   // boolean
		properties.setProperty(prefix+"show_mapped_mono",     this.show_mapped_mono+"");    // boolean
		
		properties.setProperty(prefix+"generate_raw",                  this.generate_raw+"");                  // boolean
		properties.setProperty(prefix+"generate_inf",                  this.generate_inf+"");                  // boolean
		properties.setProperty(prefix+"generate_fg",                   this.generate_fg+"");                   // boolean
		properties.setProperty(prefix+"generate_bg",                   this.generate_bg+"");                   // boolean
		properties.setProperty(prefix+"generate_stereo",               this.generate_stereo+"");               // boolean
		/*
		properties.setProperty(prefix+"stereo_bases_num",              this.stereo_bases.length+"");           // int
		for (int i = 0; i < this.stereo_bases.length; i++) {
			properties.setProperty(prefix+"stereo_bases_"+i,           this.stereo_bases[i]+"");               // double
			properties.setProperty(prefix+"generate_stereo_var_"+i,    this.generate_stereo_var[i]+"");        // boolean
		}
	   */
		properties.setProperty(prefix+"stereo_views_num",                  this.stereo_views.length+"");           // int
		for (int i = 0; i < this.stereo_views.length; i++) {
			properties.setProperty(prefix+"stereo_views_"+i,           doublesToString(this.stereo_views[i],"%.0f")); // String
			properties.setProperty(prefix+"generate_stereo_var_"+i,    this.generate_stereo_var[i]+"");        // boolean
		}
		
		
		
		properties.setProperty(prefix+"export_images",                 this.export_images + "");               // boolean
		properties.setProperty(prefix+"show_images",                   this.show_images + "");                 // boolean
		properties.setProperty(prefix+"show_images_bgfg",              this.show_images_bgfg + "");            // boolean
		properties.setProperty(prefix+"show_images_mono",              this.show_images_mono + "");            // boolean
		properties.setProperty(prefix+"export_ranges",                 this.export_ranges + "");               // boolean
		properties.setProperty(prefix+"show_ranges",                   this.show_ranges + "");                 // boolean
		properties.setProperty(prefix+"export_ml_files",               this.export_ml_files + "");                 // boolean
		
		
		properties.setProperty(prefix+"show_color_nan",                this.show_color_nan + "");              // boolean
		properties.setProperty(prefix+"show_mono_nan",                 this.show_mono_nan + "");               // boolean
		properties.setProperty(prefix+"min_num_scenes",                this.min_num_scenes+"");                // int
		properties.setProperty(prefix+"blur_egomotion",                this.blur_egomotion+"");                // double
		properties.setProperty(prefix+"range_disparity_offset",        this.range_disparity_offset+"");        // double
		properties.setProperty(prefix+"range_min_strength",            this.range_min_strength+"");            // double
		properties.setProperty(prefix+"range_max",                     this.range_max+"");                     // double
		properties.setProperty(prefix+"num_bottom",                    this.num_bottom+"");                    // int
		properties.setProperty(prefix+"num_passes",                    this.num_passes+"");                    // int
		properties.setProperty(prefix+"max_change",                    this.max_change+"");                    // double
		properties.setProperty(prefix+"min_disparity",                 this.min_disparity+"");                 // double
		properties.setProperty(prefix+"max_sym_disparity",             this.max_sym_disparity+"");             // double
		properties.setProperty(prefix+"min_strength_lma",              this.min_strength_lma+"");          // double
		properties.setProperty(prefix+"min_strength_replace",          this.min_strength_replace+"");          // double
		properties.setProperty(prefix+"min_strength_blur",             this.min_strength_blur+"");             // double
		properties.setProperty(prefix+"sigma",                         this.sigma+"");                         // double
		properties.setProperty(prefix+"num_blur",                      this.num_blur+"");                      // int
		properties.setProperty(prefix+"disparity_corr",                this.disparity_corr+"");                // double
		properties.setProperty(prefix+"outliers_nth_fromextrem",       this.outliers_nth_fromextrem+"");       // int
		properties.setProperty(prefix+"outliers_tolerance_absolute",   this.outliers_tolerance_absolute+"");   // double
		properties.setProperty(prefix+"outliers_tolerance_relative",   this.outliers_tolerance_relative+"");   // double
		properties.setProperty(prefix+"outliers_max_iter",             this.outliers_max_iter+"");             // int
		properties.setProperty(prefix+"outliers_max_strength2",        this.outliers_max_strength2+"");        // double
		properties.setProperty(prefix+"outliers_nth_fromextrem2",      this.outliers_nth_fromextrem2+"");      // int
		properties.setProperty(prefix+"outliers_tolerance_absolute2",  this.outliers_tolerance_absolute2+"");  // double
		properties.setProperty(prefix+"outliers_tolerance_relative2",  this.outliers_tolerance_relative2+"");  // double
		properties.setProperty(prefix+"outliers_lma_max_strength",     this.outliers_lma_max_strength+"");     // double
		properties.setProperty(prefix+"outliers_max_strength",         this.outliers_max_strength+"");         // double
		properties.setProperty(prefix+"outliers_from_lma_max_strength",this.outliers_from_lma_max_strength+"");// double
		properties.setProperty(prefix+"search_radius",                 this.search_radius+"");                 // int
		properties.setProperty(prefix+"remove_no_lma_neib",            this.remove_no_lma_neib+"");            // boolean
		properties.setProperty(prefix+"diff_from_lma_pos",             this.diff_from_lma_pos+"");             // double
		properties.setProperty(prefix+"diff_from_lma_neg",             this.diff_from_lma_neg+"");             // double
		properties.setProperty(prefix+"outliers_lma_nth_fromextrem",   this.outliers_lma_nth_fromextrem+"");   // int
		properties.setProperty(prefix+"filter_margin",                 this.filter_margin+"");                 // int
		
		properties.setProperty(prefix+"weak_tolerance_absolute",       this.weak_tolerance_absolute+"");       // double
		properties.setProperty(prefix+"weak_tolerance_relative",       this.weak_tolerance_relative+"");       // double
		properties.setProperty(prefix+"weak_min_neibs",                this.weak_min_neibs+"");                // int
		properties.setProperty(prefix+"strong_strength",               this.strong_strength+"");               // double
		properties.setProperty(prefix+"weak_strength",                 this.weak_strength+"");                 // double
		
		properties.setProperty(prefix+"stereo_merge",                  this.stereo_merge+"");                  // boolean
		properties.setProperty(prefix+"stereo_gap",                    this.stereo_gap+"");                    // int
		properties.setProperty(prefix+"stereo_intereye",               this.stereo_intereye+"");               // double
		properties.setProperty(prefix+"stereo_phone_width",            this.stereo_phone_width+"");            // double
		
		properties.setProperty(prefix+"extra_hor_tile",       this.extra_hor_tile+"");      // int
		properties.setProperty(prefix+"extra_vert_tile",      this.extra_vert_tile+"");     // int
		properties.setProperty(prefix+"crop_3d",              this.crop_3d+"");             // boolean
		properties.setProperty(prefix+"sensor_mask",          this.sensor_mask+"");         // int
		properties.setProperty(prefix+"merge_all",            this.merge_all+"");             // boolean
		
		properties.setProperty(prefix+"video_fps",            this.video_fps+"");           // double
		properties.setProperty(prefix+"sensor_fps",           this.sensor_fps+"");          // double
		properties.setProperty(prefix+"mode_avi",             this.mode_avi+"");            // int
		properties.setProperty(prefix+"avi_JPEG_quality",     this.avi_JPEG_quality+"");    // int
		properties.setProperty(prefix+"run_ffmpeg",           this.run_ffmpeg+"");          // boolean
		properties.setProperty(prefix+"video_ext",            this.video_ext+"");           // String
		properties.setProperty(prefix+"video_codec",          this.video_codec+"");         // String
		properties.setProperty(prefix+"video_crf",            this.video_crf+"");           // int
		properties.setProperty(prefix+"remove_avi",           this.remove_avi+"");          // boolean

		properties.setProperty(prefix+"video_codec_combo",    this.video_codec_combo+"");         // String
		properties.setProperty(prefix+"video_crf_combo",      this.video_crf_combo+"");           // int
		
		properties.setProperty(prefix+"um_mono",              this.um_mono+"");             // boolean
		properties.setProperty(prefix+"um_sigma",             this.um_sigma+"");            // double
		properties.setProperty(prefix+"um_weight",            this.um_weight+"");           // double
		properties.setProperty(prefix+"mono_fixed",           this.mono_fixed+"");          // boolean
		properties.setProperty(prefix+"mono_range",           this.mono_range+"");          // double

		properties.setProperty(prefix+"anaglyth_en",          this.anaglyth_en+"");         // boolean
		
		properties.setProperty(prefix+"anaglyph_left",        getLongColor(anaglyph_left)+""); // Color
		properties.setProperty(prefix+"anaglyph_right",       getLongColor(anaglyph_right)+""); // Color
		
		
		properties.setProperty(prefix+"annotate_color",       this.annotate_color+"");      // boolean
		properties.setProperty(prefix+"annotate_mono",        this.annotate_mono+"");       // boolean
		{
			long lcolor_annotate =        (annotate_color_color == null) ? -1 : getLongColor(annotate_color_color);
			properties.setProperty(prefix+"annotate_color_color",              lcolor_annotate+"");
		}
		{
			long lcolor_annotate =        (annotate_color_mono == null) ? -1 : getLongColor(annotate_color_mono);
			properties.setProperty(prefix+"annotate_color_mono",              lcolor_annotate+"");
		}
		properties.setProperty(prefix+"annotate_transparent_mono",this.annotate_transparent_mono+"");// boolean
		properties.setProperty(prefix+"margin",               this.margin+"");              // int
		properties.setProperty(prefix+"sensor_mask_inter",    this.sensor_mask_inter+"");   // int
		properties.setProperty(prefix+"use_partial",          this.use_partial+"");         // boolean
		properties.setProperty(prefix+"run_poly",             this.run_poly+"");            // boolean
		properties.setProperty(prefix+"centroid_radius",      this.centroid_radius+"");     // double
		properties.setProperty(prefix+"n_recenter",           this.n_recenter+"");          // int

		properties.setProperty(prefix+"min_ref_str",          this.min_ref_str+"");           // double
		
		properties.setProperty(prefix+"td_weight",            this.td_weight+"");           // double
		properties.setProperty(prefix+"pd_weight",            this.pd_weight+"");           // double
		properties.setProperty(prefix+"td_nopd_only",         this.td_nopd_only+"");        // boolean
		
		properties.setProperty(prefix+"min_str",              this.min_str+"");             // double
		properties.setProperty(prefix+"min_str_fpn",          this.min_str_fpn+"");         // double
		properties.setProperty(prefix+"min_str_sum",          this.min_str_sum+"");         // double
		properties.setProperty(prefix+"min_str_sum_fpn",      this.min_str_sum_fpn+"");     // double
		
		properties.setProperty(prefix+"min_neibs",            this.min_neibs+"");           // int
		properties.setProperty(prefix+"weight_zero_neibs",    this.weight_zero_neibs+"");   // double
		properties.setProperty(prefix+"half_disparity",       this.half_disparity+"");      // double
		properties.setProperty(prefix+"half_avg_diff",        this.half_avg_diff+"");       // double
		
		properties.setProperty(prefix+"photo_en",             this.photo_en+"");            // boolean
		properties.setProperty(prefix+"photo_num_full",       this.photo_num_full+"");      // int
		properties.setProperty(prefix+"photo_num_refines",    this.photo_num_refines+"");   // int
		properties.setProperty(prefix+"photo_min_strength",   this.photo_min_strength+"");  // double
		properties.setProperty(prefix+"photo_max_diff",       this.photo_max_diff+"");      // double
		properties.setProperty(prefix+"photo_debug",          this.photo_debug+"");         // boolean
		
		properties.setProperty(prefix+"pix_step",             this.pix_step+"");            // int
		properties.setProperty(prefix+"search_rad",           this.search_rad+"");          // int
		properties.setProperty(prefix+"maybe_sum",            this.maybe_sum+"");           // double
		properties.setProperty(prefix+"sure_sum",             this.sure_sum+"");            // double
		properties.setProperty(prefix+"maybe_avg",            this.maybe_avg+"");           // double
		properties.setProperty(prefix+"sure_avg",             this.sure_avg+"");            // double
		
		properties.setProperty(prefix+"max_search_rms",       this.max_search_rms+"");      // double
		properties.setProperty(prefix+"maybe_fom",            this.maybe_fom+"");           // double
		properties.setProperty(prefix+"sure_fom",             this.sure_fom+"");            // double
		
		properties.setProperty(prefix+"use_combo_dsi",        this.use_combo_dsi+"");       // boolean
		properties.setProperty(prefix+"use_lma_dsi",          this.use_lma_dsi+"");         // boolean
		
		properties.setProperty(prefix+"fpn_remove",           this.fpn_remove+"");          // boolean
		properties.setProperty(prefix+"fpn_max_offset",       this.fpn_max_offset+"");      // double
		properties.setProperty(prefix+"fpn_radius",           this.fpn_radius+"");          // double
		properties.setProperty(prefix+"fpn_ignore_border",    this.fpn_ignore_border+"");   // boolean
		
		properties.setProperty(prefix+"mov_en",               this.mov_en+"");              // boolean
		properties.setProperty(prefix+"mov_sigma",            this.mov_sigma+"");           // double
		properties.setProperty(prefix+"mov_max_std",          this.mov_max_std+"");         // double
		properties.setProperty(prefix+"mov_thresh_rel",       this.mov_thresh_rel+"");      // double
		properties.setProperty(prefix+"mov_thresh_abs",       this.mov_thresh_abs+"");      // double
		properties.setProperty(prefix+"mov_clust_max",        this.mov_clust_max+"");       // double
		properties.setProperty(prefix+"mov_grow",             this.mov_grow+"");            // int
		properties.setProperty(prefix+"mov_show",             this.mov_show+"");            // boolean
		properties.setProperty(prefix+"mov_debug_level",      this.mov_debug_level+"");     // int
		
		properties.setProperty(prefix+"adjust_atr_0",         this.adjust_atr[0]+"");       // boolean
		properties.setProperty(prefix+"adjust_atr_1",         this.adjust_atr[1]+"");       // boolean
		properties.setProperty(prefix+"adjust_atr_2",         this.adjust_atr[2]+"");       // boolean
		properties.setProperty(prefix+"adjust_xyz_0",         this.adjust_xyz[0]+"");       // boolean
		properties.setProperty(prefix+"adjust_xyz_1",         this.adjust_xyz[1]+"");       // boolean
		properties.setProperty(prefix+"adjust_xyz_2",         this.adjust_xyz[2]+"");       // boolean
		properties.setProperty(prefix+"exit_change_atr",      this.exit_change_atr+"");     // double
		properties.setProperty(prefix+"exit_change_xyz",      this.exit_change_xyz+"");     // double
		properties.setProperty(prefix+"max_cycles",           this.max_cycles+"");          // int
		properties.setProperty(prefix+"max_LMA",              this.max_LMA+"");             // int
		properties.setProperty(prefix+"max_rms",              this.max_rms+"");             // double
		properties.setProperty(prefix+"scene_is_ref_test",    this.scene_is_ref_test+"");   // boolean
		properties.setProperty(prefix+"render_ref",           this.render_ref+"");          // boolean
		properties.setProperty(prefix+"render_scene",         this.render_scene+"");        // boolean
		properties.setProperty(prefix+"toRGB",                this.toRGB+"");               // boolean
		properties.setProperty(prefix+"show_2d_correlations", this.show_2d_correlations+"");// boolean
		properties.setProperty(prefix+"show_motion_vectors",  this.show_motion_vectors+""); // boolean
		properties.setProperty(prefix+"debug_level",          this.debug_level+"");         // int
		properties.setProperty(prefix+"test_ers",             this.test_ers+"");            // boolean
		properties.setProperty(prefix+"test_ers0",            this.test_ers0+"");           // int
		properties.setProperty(prefix+"test_ers1",            this.test_ers1+"");           // int
	
	}
	
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"force_ref_dsi")!=null)                 this.force_ref_dsi=Boolean.parseBoolean(properties.getProperty(prefix+"force_ref_dsi"));		
		if (properties.getProperty(prefix+"force_orientations")!=null)            this.force_orientations=Boolean.parseBoolean(properties.getProperty(prefix+"force_orientations"));		
//		if (properties.getProperty(prefix+"readjust_orient")!=null)               this.readjust_orient=Boolean.parseBoolean(properties.getProperty(prefix+"readjust_orient"));		
//		if (properties.getProperty(prefix+"force_interscene")!=null)              this.force_interscene=Boolean.parseBoolean(properties.getProperty(prefix+"force_interscene"));		
		if (properties.getProperty(prefix+"min_num_orient")!=null)                this.min_num_orient=Integer.parseInt(properties.getProperty(prefix+"min_num_orient"));
		if (properties.getProperty(prefix+"min_num_interscene")!=null)            this.min_num_interscene=Integer.parseInt(properties.getProperty(prefix+"min_num_interscene"));
		if (properties.getProperty(prefix+"generate_mapped")!=null)               this.generate_mapped=Boolean.parseBoolean(properties.getProperty(prefix+"generate_mapped"));
		if (properties.getProperty(prefix+"reuse_video")!=null)                   this.reuse_video=Boolean.parseBoolean(properties.getProperty(prefix+"reuse_video"));
		
		if (properties.getProperty(prefix+"save_mapped_color")!=null)             this.save_mapped_color=Boolean.parseBoolean(properties.getProperty(prefix+"save_mapped_color"));
		if (properties.getProperty(prefix+"save_mapped_mono")!=null)              this.save_mapped_mono=Boolean.parseBoolean(properties.getProperty(prefix+"save_mapped_mono"));
		if (properties.getProperty(prefix+"gen_avi_color")!=null)                 this.gen_avi_color=Boolean.parseBoolean(properties.getProperty(prefix+"gen_avi_color"));
		if (properties.getProperty(prefix+"gen_avi_mono")!=null)                  this.gen_avi_mono=Boolean.parseBoolean(properties.getProperty(prefix+"gen_avi_mono"));
		if (properties.getProperty(prefix+"show_mapped_color")!=null)    this.show_mapped_color=Boolean.parseBoolean(properties.getProperty(prefix+"show_mapped_color"));
		if (properties.getProperty(prefix+"show_mapped_mono")!=null)     this.show_mapped_mono=Boolean.parseBoolean(properties.getProperty(prefix+"show_mapped_mono"));
		
		if (properties.getProperty(prefix+"generate_raw")!=null)                  this.generate_raw=Boolean.parseBoolean(properties.getProperty(prefix+"generate_raw"));
		if (properties.getProperty(prefix+"generate_inf")!=null)                  this.generate_inf=Boolean.parseBoolean(properties.getProperty(prefix+"generate_inf"));
		if (properties.getProperty(prefix+"generate_fg")!=null)                   this.generate_fg=Boolean.parseBoolean(properties.getProperty(prefix+"generate_fg"));
		if (properties.getProperty(prefix+"generate_bg")!=null)                   this.generate_bg=Boolean.parseBoolean(properties.getProperty(prefix+"generate_bg"));
		
		if (properties.getProperty(prefix+"generate_stereo")!=null)               this.generate_stereo=Boolean.parseBoolean(properties.getProperty(prefix+"generate_stereo"));
		/*
		if (properties.getProperty(prefix+"stereo_bases_num")!=null) {
			 int stereo_bases_num=Integer.parseInt(properties.getProperty(prefix+"stereo_bases_num"));
			 this.stereo_bases =        new double[stereo_bases_num];
			 this.generate_stereo_var = new boolean[stereo_bases_num];
			 for (int i = 0; i < stereo_bases_num; i++) {
				 if (properties.getProperty(prefix+"stereo_bases_"+i)!=null)       this.stereo_bases[i]=Double.parseDouble(properties.getProperty(prefix+"stereo_bases_"+i));
				 if (properties.getProperty(prefix+"generate_stereo_var_"+i)!=null)this.generate_stereo_var[i]=Boolean.parseBoolean(properties.getProperty(prefix+"generate_stereo_var_"+i));
			 }
			 orderStereo(); 
		}
		*/
		if (properties.getProperty(prefix+"stereo_views_num")!=null) {
			 int stereo_views_num=Integer.parseInt(properties.getProperty(prefix+"stereo_views_num"));
			 this.stereo_views =        new double[stereo_views_num][];
			 this.generate_stereo_var = new boolean[stereo_views_num];
			 for (int i = 0; i < stereo_views_num; i++) {
				 if (properties.getProperty(prefix+"stereo_views_"+i)!=null) {
					 this.stereo_views[i]=StringToDoubles(properties.getProperty(prefix+"stereo_views_"+i),3);
				 }
				 if (properties.getProperty(prefix+"generate_stereo_var_"+i)!=null) {
					 this.generate_stereo_var[i]=Boolean.parseBoolean(properties.getProperty(prefix+"generate_stereo_var_"+i));
				 }
			 }
			 orderStereoViews(); 
		}
		
		
		
		
		if (properties.getProperty(prefix+"export_images")!=null)                 this.export_images=Boolean.parseBoolean(properties.getProperty(prefix+"export_images"));		
		if (properties.getProperty(prefix+"show_images")!=null)                   this.show_images=Boolean.parseBoolean(properties.getProperty(prefix+"show_images"));		
		if (properties.getProperty(prefix+"show_images_bgfg")!=null)              this.show_images_bgfg=Boolean.parseBoolean(properties.getProperty(prefix+"show_images_bgfg"));		
		if (properties.getProperty(prefix+"show_images_mono")!=null)              this.show_images_mono=Boolean.parseBoolean(properties.getProperty(prefix+"show_images_mono"));		
		if (properties.getProperty(prefix+"export_ranges")!=null)                 this.export_ranges=Boolean.parseBoolean(properties.getProperty(prefix+"export_ranges"));
		if (properties.getProperty(prefix+"show_ranges")!=null)                   this.show_ranges=Boolean.parseBoolean(properties.getProperty(prefix+"show_ranges"));
		if (properties.getProperty(prefix+"export_ml_files")!=null)               this.export_ml_files=Boolean.parseBoolean(properties.getProperty(prefix+"export_ml_files"));

		if (properties.getProperty(prefix+"show_color_nan")!=null)                this.show_color_nan=Boolean.parseBoolean(properties.getProperty(prefix+"show_color_nan"));		
		if (properties.getProperty(prefix+"show_mono_nan")!=null)                 this.show_mono_nan=Boolean.parseBoolean(properties.getProperty(prefix+"show_mono_nan"));
		if (properties.getProperty(prefix+"min_num_scenes")!=null)                this.min_num_scenes=Integer.parseInt(properties.getProperty(prefix+"min_num_scenes"));
		if (properties.getProperty(prefix+"blur_egomotion")!=null)                this.blur_egomotion=Double.parseDouble(properties.getProperty(prefix+"blur_egomotion"));
		if (properties.getProperty(prefix+"range_disparity_offset")!=null)        this.range_disparity_offset=Double.parseDouble(properties.getProperty(prefix+"range_disparity_offset"));
		if (properties.getProperty(prefix+"range_min_strength")!=null)            this.range_min_strength=Double.parseDouble(properties.getProperty(prefix+"range_min_strength"));
		if (properties.getProperty(prefix+"range_max")!=null)                     this.range_max=Double.parseDouble(properties.getProperty(prefix+"range_max"));
		if (properties.getProperty(prefix+"num_bottom")!=null)                    this.num_bottom=Integer.parseInt(properties.getProperty(prefix+"num_bottom"));
		if (properties.getProperty(prefix+"num_passes")!=null)                    this.num_passes=Integer.parseInt(properties.getProperty(prefix+"num_passes"));
		if (properties.getProperty(prefix+"max_change")!=null)                    this.max_change=Double.parseDouble(properties.getProperty(prefix+"max_change"));
		if (properties.getProperty(prefix+"min_disparity")!=null)                 this.min_disparity=Double.parseDouble(properties.getProperty(prefix+"min_disparity"));
		if (properties.getProperty(prefix+"max_sym_disparity")!=null)             this.max_sym_disparity=Double.parseDouble(properties.getProperty(prefix+"max_sym_disparity"));
		if (properties.getProperty(prefix+"min_strength_lma")!=null)              this.min_strength_lma=Double.parseDouble(properties.getProperty(prefix+"min_strength_lma"));
		if (properties.getProperty(prefix+"min_strength_replace")!=null)          this.min_strength_replace=Double.parseDouble(properties.getProperty(prefix+"min_strength_replace"));
		if (properties.getProperty(prefix+"min_strength_blur")!=null)             this.min_strength_blur=Double.parseDouble(properties.getProperty(prefix+"min_strength_blur"));
		if (properties.getProperty(prefix+"sigma")!=null)                         this.sigma=Double.parseDouble(properties.getProperty(prefix+"sigma"));
		if (properties.getProperty(prefix+"num_blur")!=null)                      this.num_blur=Integer.parseInt(properties.getProperty(prefix+"num_blur"));
		if (properties.getProperty(prefix+"disparity_corr")!=null)                this.disparity_corr=Double.parseDouble(properties.getProperty(prefix+"disparity_corr"));
		if (properties.getProperty(prefix+"outliers_nth_fromextrem")!=null)       this.outliers_nth_fromextrem=Integer.parseInt(properties.getProperty(prefix+"outliers_nth_fromextrem"));
		if (properties.getProperty(prefix+"outliers_tolerance_absolute")!=null)   this.outliers_tolerance_absolute=Double.parseDouble(properties.getProperty(prefix+"outliers_tolerance_absolute"));
		if (properties.getProperty(prefix+"outliers_tolerance_relative")!=null)   this.outliers_tolerance_relative=Double.parseDouble(properties.getProperty(prefix+"outliers_tolerance_relative"));
		if (properties.getProperty(prefix+"outliers_max_iter")!=null)             this.outliers_max_iter=Integer.parseInt(properties.getProperty(prefix+"outliers_max_iter"));
		if (properties.getProperty(prefix+"outliers_max_strength2")!=null)        this.outliers_max_strength2=Double.parseDouble(properties.getProperty(prefix+"outliers_max_strength2"));
		if (properties.getProperty(prefix+"outliers_nth_fromextrem2")!=null)      this.outliers_nth_fromextrem2=Integer.parseInt(properties.getProperty(prefix+"outliers_nth_fromextrem2"));
		if (properties.getProperty(prefix+"outliers_tolerance_absolute2")!=null)  this.outliers_tolerance_absolute2=Double.parseDouble(properties.getProperty(prefix+"outliers_tolerance_absolute2"));
		if (properties.getProperty(prefix+"outliers_tolerance_relative2")!=null)  this.outliers_tolerance_relative2=Double.parseDouble(properties.getProperty(prefix+"outliers_tolerance_relative2"));
		if (properties.getProperty(prefix+"outliers_lma_max_strength")!=null)     this.outliers_lma_max_strength=Double.parseDouble(properties.getProperty(prefix+"outliers_lma_max_strength"));
		if (properties.getProperty(prefix+"outliers_max_strength")!=null)         this.outliers_max_strength=Double.parseDouble(properties.getProperty(prefix+"outliers_max_strength"));
		if (properties.getProperty(prefix+"outliers_from_lma_max_strength")!=null)this.outliers_from_lma_max_strength=Double.parseDouble(properties.getProperty(prefix+"outliers_from_lma_max_strength"));
		if (properties.getProperty(prefix+"search_radius")!=null)                 this.search_radius=Integer.parseInt(properties.getProperty(prefix+"search_radius"));
		if (properties.getProperty(prefix+"remove_no_lma_neib")!=null)            this.remove_no_lma_neib=Boolean.parseBoolean(properties.getProperty(prefix+"remove_no_lma_neib"));
		if (properties.getProperty(prefix+"diff_from_lma_pos")!=null)             this.diff_from_lma_pos=Double.parseDouble(properties.getProperty(prefix+"diff_from_lma_pos"));
		if (properties.getProperty(prefix+"diff_from_lma_neg")!=null)             this.diff_from_lma_neg=Double.parseDouble(properties.getProperty(prefix+"diff_from_lma_neg"));
		if (properties.getProperty(prefix+"outliers_lma_nth_fromextrem")!=null)   this.outliers_lma_nth_fromextrem=Integer.parseInt(properties.getProperty(prefix+"outliers_lma_nth_fromextrem"));
		if (properties.getProperty(prefix+"filter_margin")!=null)                 this.filter_margin=Integer.parseInt(properties.getProperty(prefix+"filter_margin"));

		if (properties.getProperty(prefix+"weak_tolerance_absolute")!=null)       this.weak_tolerance_absolute=Double.parseDouble(properties.getProperty(prefix+"weak_tolerance_absolute"));
		if (properties.getProperty(prefix+"weak_tolerance_relative")!=null)       this.weak_tolerance_relative=Double.parseDouble(properties.getProperty(prefix+"weak_tolerance_relative"));
		if (properties.getProperty(prefix+"weak_min_neibs")!=null)                this.weak_min_neibs=Integer.parseInt(properties.getProperty(prefix+"weak_min_neibs"));
		if (properties.getProperty(prefix+"strong_strength")!=null)               this.strong_strength=Double.parseDouble(properties.getProperty(prefix+"strong_strength"));
		if (properties.getProperty(prefix+"weak_strength")!=null)                 this.weak_strength=Double.parseDouble(properties.getProperty(prefix+"weak_strength"));
		
		if (properties.getProperty(prefix+"stereo_merge")!=null)                  this.stereo_merge=Boolean.parseBoolean(properties.getProperty(prefix+"stereo_merge"));
		if (properties.getProperty(prefix+"stereo_gap")!=null)                    this.stereo_gap=Integer.parseInt(properties.getProperty(prefix+"stereo_gap"));
		if (properties.getProperty(prefix+"stereo_intereye")!=null)               this.stereo_intereye=Double.parseDouble(properties.getProperty(prefix+"stereo_intereye"));
		if (properties.getProperty(prefix+"stereo_phone_width")!=null)            this.stereo_phone_width=Double.parseDouble(properties.getProperty(prefix+"stereo_phone_width"));
		
		if (properties.getProperty(prefix+"extra_hor_tile")!=null)       this.extra_hor_tile=Integer.parseInt(properties.getProperty(prefix+"extra_hor_tile"));
		if (properties.getProperty(prefix+"extra_vert_tile")!=null)      this.extra_vert_tile=Integer.parseInt(properties.getProperty(prefix+"extra_vert_tile"));
		if (properties.getProperty(prefix+"crop_3d")!=null)              this.crop_3d=Boolean.parseBoolean(properties.getProperty(prefix+"crop_3d"));		
		if (properties.getProperty(prefix+"sensor_mask")!=null)          this.sensor_mask=Integer.parseInt(properties.getProperty(prefix+"sensor_mask"));
		if (properties.getProperty(prefix+"merge_all")!=null)            this.merge_all=Boolean.parseBoolean(properties.getProperty(prefix+"merge_all"));		
		
		if (properties.getProperty(prefix+"video_fps")!=null)            this.video_fps=Double.parseDouble(properties.getProperty(prefix+"video_fps"));
		if (properties.getProperty(prefix+"sensor_fps")!=null)           this.sensor_fps=Double.parseDouble(properties.getProperty(prefix+"sensor_fps"));
		if (properties.getProperty(prefix+"mode_avi")!=null)             this.mode_avi=Integer.parseInt(properties.getProperty(prefix+"mode_avi"));
		if (properties.getProperty(prefix+"avi_JPEG_quality")!=null)     this.avi_JPEG_quality=Integer.parseInt(properties.getProperty(prefix+"avi_JPEG_quality"));
		if (properties.getProperty(prefix+"run_ffmpeg")!=null)           this.run_ffmpeg=Boolean.parseBoolean(properties.getProperty(prefix+"run_ffmpeg"));
		if (properties.getProperty(prefix+"video_ext")!=null)            this.video_ext=(String) properties.getProperty(prefix+"video_ext");
		if (properties.getProperty(prefix+"video_codec")!=null)          this.video_codec=(String) properties.getProperty(prefix+"video_codec");
		if (properties.getProperty(prefix+"video_crf")!=null)            this.video_crf=Integer.parseInt(properties.getProperty(prefix+"video_crf"));
		if (properties.getProperty(prefix+"remove_avi")!=null)           this.remove_avi=Boolean.parseBoolean(properties.getProperty(prefix+"remove_avi"));
		if (properties.getProperty(prefix+"video_codec_combo")!=null)    this.video_codec_combo=(String) properties.getProperty(prefix+"video_codec_combo");
		if (properties.getProperty(prefix+"video_crf_combo")!=null)      this.video_crf_combo=Integer.parseInt(properties.getProperty(prefix+"video_crf_combo"));
		if (properties.getProperty(prefix+"um_mono")!=null)              this.um_mono=Boolean.parseBoolean(properties.getProperty(prefix+"um_mono"));
		if (properties.getProperty(prefix+"um_sigma")!=null)             this.um_sigma=Double.parseDouble(properties.getProperty(prefix+"um_sigma"));
		if (properties.getProperty(prefix+"um_weight")!=null)            this.um_weight=Double.parseDouble(properties.getProperty(prefix+"um_weight"));
		if (properties.getProperty(prefix+"mono_fixed")!=null)           this.mono_fixed=Boolean.parseBoolean(properties.getProperty(prefix+"mono_fixed"));
		if (properties.getProperty(prefix+"mono_range")!=null)           this.mono_range=Double.parseDouble(properties.getProperty(prefix+"mono_range"));
		
		if (properties.getProperty(prefix+"anaglyth_en")!=null)           this.anaglyth_en=Boolean.parseBoolean(properties.getProperty(prefix+"anaglyth_en"));
		if  (properties.getProperty(prefix+"anaglyph_left") != null) {
			try {
				this.anaglyph_left = setLongColor(Long.parseLong(properties.getProperty(prefix+"anaglyph_left")));
			} catch(NumberFormatException e){
				this.anaglyph_left = anaglyph_left_default;
			}
		}
		if  (properties.getProperty(prefix+"anaglyph_right") != null) {
			try {
				this.anaglyph_right = setLongColor(Long.parseLong(properties.getProperty(prefix+"anaglyph_right")));
			} catch(NumberFormatException e){
				this.anaglyph_right = anaglyph_right_default;
			}
		}
		
		if (properties.getProperty(prefix+"annotate_color")!=null)       this.annotate_color=Boolean.parseBoolean(properties.getProperty(prefix+"annotate_color"));
		if (properties.getProperty(prefix+"annotate_mono")!=null)        this.annotate_mono=Boolean.parseBoolean(properties.getProperty(prefix+"annotate_mono"));

		if  (properties.getProperty(prefix+"annotate_color_color")!=null) {
			long lcolor_annotate = Long.parseLong(properties.getProperty(prefix+"annotate_color_color"));
			if (lcolor_annotate < 0) this.annotate_color_color = null;
			else this.annotate_color_color = setLongColor(lcolor_annotate);
		}
		if  (properties.getProperty(prefix+"annotate_color_mono")!=null) {
			long lcolor_annotate = Long.parseLong(properties.getProperty(prefix+"annotate_color_mono"));
			if (lcolor_annotate < 0) this.annotate_color_mono = null;
			else this.annotate_color_mono = setLongColor(lcolor_annotate);
		}
		if (properties.getProperty(prefix+"annotate_transparent_mono")!=null) this.annotate_transparent_mono=Boolean.parseBoolean(properties.getProperty(prefix+"annotate_transparent_mono"));
		
		if (properties.getProperty(prefix+"margin")!=null)               this.margin=Integer.parseInt(properties.getProperty(prefix+"margin"));
		if (properties.getProperty(prefix+"sensor_mask_inter")!=null)    this.sensor_mask_inter=Integer.parseInt(properties.getProperty(prefix+"sensor_mask_inter"));
		if (properties.getProperty(prefix+"use_partial")!=null)          this.use_partial=Boolean.parseBoolean(properties.getProperty(prefix+"use_partial"));		
		if (properties.getProperty(prefix+"run_poly")!=null)             this.run_poly=Boolean.parseBoolean(properties.getProperty(prefix+"run_poly"));
		if (properties.getProperty(prefix+"centroid_radius")!=null)      this.centroid_radius=Double.parseDouble(properties.getProperty(prefix+"centroid_radius"));
		if (properties.getProperty(prefix+"n_recenter")!=null)           this.n_recenter=Integer.parseInt(properties.getProperty(prefix+"n_recenter"));
		
		if (properties.getProperty(prefix+"min_ref_str")!=null)          this.min_ref_str=Double.parseDouble(properties.getProperty(prefix+"min_ref_str"));

		if (properties.getProperty(prefix+"td_weight")!=null)            this.td_weight=Double.parseDouble(properties.getProperty(prefix+"td_weight"));
		if (properties.getProperty(prefix+"pd_weight")!=null)            this.pd_weight=Double.parseDouble(properties.getProperty(prefix+"pd_weight"));
		if (properties.getProperty(prefix+"td_nopd_only")!=null)         this.td_nopd_only=Boolean.parseBoolean(properties.getProperty(prefix+"td_nopd_only"));
		
		if (properties.getProperty(prefix+"min_str")!=null)              this.min_str=Double.parseDouble(properties.getProperty(prefix+"min_str"));
		if (properties.getProperty(prefix+"min_str_fpn")!=null)          this.min_str_fpn=Double.parseDouble(properties.getProperty(prefix+"min_str_fpn"));
		if (properties.getProperty(prefix+"min_str_sum")!=null)          this.min_str_sum=Double.parseDouble(properties.getProperty(prefix+"min_str_sum"));
		if (properties.getProperty(prefix+"min_str_sum_fpn")!=null)      this.min_str_sum_fpn=Double.parseDouble(properties.getProperty(prefix+"min_str_sum_fpn"));

		if (properties.getProperty(prefix+"min_neibs")!=null)            this.min_neibs=Integer.parseInt(properties.getProperty(prefix+"min_neibs"));
		if (properties.getProperty(prefix+"weight_zero_neibs")!=null)    this.weight_zero_neibs=Double.parseDouble(properties.getProperty(prefix+"weight_zero_neibs"));
		if (properties.getProperty(prefix+"half_disparity")!=null)       this.half_disparity=Double.parseDouble(properties.getProperty(prefix+"half_disparity"));
		if (properties.getProperty(prefix+"half_avg_diff")!=null)        this.half_avg_diff=Double.parseDouble(properties.getProperty(prefix+"half_avg_diff"));
		
		if (properties.getProperty(prefix+"photo_en")!=null)             this.photo_en=Boolean.parseBoolean(properties.getProperty(prefix+"photo_en"));		
		if (properties.getProperty(prefix+"photo_num_full")!=null)       this.photo_num_full=Integer.parseInt(properties.getProperty(prefix+"photo_num_full"));
		if (properties.getProperty(prefix+"photo_num_refines")!=null)    this.photo_num_refines=Integer.parseInt(properties.getProperty(prefix+"photo_num_refines"));
		if (properties.getProperty(prefix+"photo_min_strength")!=null)   this.photo_min_strength=Double.parseDouble(properties.getProperty(prefix+"photo_min_strength"));
		if (properties.getProperty(prefix+"photo_max_diff")!=null)       this.photo_max_diff=Double.parseDouble(properties.getProperty(prefix+"photo_max_diff"));
		if (properties.getProperty(prefix+"photo_debug")!=null)          this.photo_debug=Boolean.parseBoolean(properties.getProperty(prefix+"photo_debug"));		
		
		if (properties.getProperty(prefix+"pix_step")!=null)             this.pix_step=Integer.parseInt(properties.getProperty(prefix+"pix_step"));
		if (properties.getProperty(prefix+"search_rad")!=null)           this.search_rad=Integer.parseInt(properties.getProperty(prefix+"search_rad"));
		if (properties.getProperty(prefix+"maybe_sum")!=null)            this.maybe_sum=Double.parseDouble(properties.getProperty(prefix+"maybe_sum"));
		if (properties.getProperty(prefix+"sure_sum")!=null)             this.sure_sum=Double.parseDouble(properties.getProperty(prefix+"sure_sum"));
		if (properties.getProperty(prefix+"maybe_avg")!=null)            this.maybe_avg=Double.parseDouble(properties.getProperty(prefix+"maybe_avg"));
		if (properties.getProperty(prefix+"sure_avg")!=null)             this.sure_avg=Double.parseDouble(properties.getProperty(prefix+"sure_avg"));

		if (properties.getProperty(prefix+"max_search_rms")!=null)       this.max_search_rms=Double.parseDouble(properties.getProperty(prefix+"max_search_rms"));
		if (properties.getProperty(prefix+"maybe_fom")!=null)            this.maybe_fom=Double.parseDouble(properties.getProperty(prefix+"maybe_fom"));
		if (properties.getProperty(prefix+"sure_fom")!=null)             this.sure_fom=Double.parseDouble(properties.getProperty(prefix+"sure_fom"));
		
		
		
		if (properties.getProperty(prefix+"use_combo_dsi")!=null)        this.use_combo_dsi=Boolean.parseBoolean(properties.getProperty(prefix+"use_combo_dsi"));		
		if (properties.getProperty(prefix+"use_lma_dsi")!=null)          this.use_lma_dsi=Boolean.parseBoolean(properties.getProperty(prefix+"use_lma_dsi"));

		if (properties.getProperty(prefix+"fpn_remove")!=null)           this.fpn_remove=Boolean.parseBoolean(properties.getProperty(prefix+"fpn_remove"));
		if (properties.getProperty(prefix+"fpn_max_offset")!=null)       this.fpn_max_offset=Double.parseDouble(properties.getProperty(prefix+"fpn_max_offset"));
		if (properties.getProperty(prefix+"fpn_radius")!=null)           this.fpn_radius=Double.parseDouble(properties.getProperty(prefix+"fpn_radius"));
		if (properties.getProperty(prefix+"fpn_ignore_border")!=null)    this.fpn_ignore_border=Boolean.parseBoolean(properties.getProperty(prefix+"fpn_ignore_border"));
		
		if (properties.getProperty(prefix+"mov_en")!=null)               this.mov_en=Boolean.parseBoolean(properties.getProperty(prefix+"mov_en"));
		if (properties.getProperty(prefix+"mov_sigma")!=null)            this.mov_sigma=Double.parseDouble(properties.getProperty(prefix+"mov_sigma"));
		if (properties.getProperty(prefix+"mov_max_std")!=null)          this.mov_max_std=Double.parseDouble(properties.getProperty(prefix+"mov_max_std"));
		if (properties.getProperty(prefix+"mov_thresh_rel")!=null)       this.mov_thresh_rel=Double.parseDouble(properties.getProperty(prefix+"mov_thresh_rel"));
		if (properties.getProperty(prefix+"mov_thresh_abs")!=null)       this.mov_thresh_abs=Double.parseDouble(properties.getProperty(prefix+"mov_thresh_abs"));
		if (properties.getProperty(prefix+"mov_clust_max")!=null)        this.mov_clust_max=Double.parseDouble(properties.getProperty(prefix+"mov_clust_max"));
		if (properties.getProperty(prefix+"mov_grow")!=null)             this.mov_grow=Integer.parseInt(properties.getProperty(prefix+"mov_grow"));
		if (properties.getProperty(prefix+"mov_show")!=null)             this.mov_show=Boolean.parseBoolean(properties.getProperty(prefix+"mov_show"));
		if (properties.getProperty(prefix+"mov_debug_level")!=null)      this.mov_debug_level=Integer.parseInt(properties.getProperty(prefix+"mov_debug_level"));
		
		if (properties.getProperty(prefix+"adjust_atr_0")!=null)         this.adjust_atr[0]=Boolean.parseBoolean(properties.getProperty(prefix+"adjust_atr_0"));		
		if (properties.getProperty(prefix+"adjust_atr_1")!=null)         this.adjust_atr[1]=Boolean.parseBoolean(properties.getProperty(prefix+"adjust_atr_1"));		
		if (properties.getProperty(prefix+"adjust_atr_2")!=null)         this.adjust_atr[2]=Boolean.parseBoolean(properties.getProperty(prefix+"adjust_atr_2"));		
		if (properties.getProperty(prefix+"adjust_xyz_0")!=null)         this.adjust_xyz[0]=Boolean.parseBoolean(properties.getProperty(prefix+"adjust_xyz_0"));		
		if (properties.getProperty(prefix+"adjust_xyz_1")!=null)         this.adjust_xyz[1]=Boolean.parseBoolean(properties.getProperty(prefix+"adjust_xyz_1"));		
		if (properties.getProperty(prefix+"adjust_xyz_2")!=null)         this.adjust_xyz[2]=Boolean.parseBoolean(properties.getProperty(prefix+"adjust_xyz_2"));		
		if (properties.getProperty(prefix+"exit_change_atr")!=null)      this.exit_change_atr=Double.parseDouble(properties.getProperty(prefix+"exit_change_atr"));
		if (properties.getProperty(prefix+"exit_change_xyz")!=null)      this.exit_change_xyz=Double.parseDouble(properties.getProperty(prefix+"exit_change_xyz"));
		if (properties.getProperty(prefix+"max_cycles")!=null)           this.max_cycles=Integer.parseInt(properties.getProperty(prefix+"max_cycles"));
		if (properties.getProperty(prefix+"max_LMA")!=null)              this.max_LMA=Integer.parseInt(properties.getProperty(prefix+"max_LMA"));
		if (properties.getProperty(prefix+"max_rms")!=null)              this.max_rms=Double.parseDouble(properties.getProperty(prefix+"max_rms"));
		if (properties.getProperty(prefix+"scene_is_ref_test")!=null)    this.scene_is_ref_test=Boolean.parseBoolean(properties.getProperty(prefix+"scene_is_ref_test"));		
		if (properties.getProperty(prefix+"render_ref")!=null)           this.render_ref=Boolean.parseBoolean(properties.getProperty(prefix+"render_ref"));		
		if (properties.getProperty(prefix+"render_scene")!=null)         this.render_scene=Boolean.parseBoolean(properties.getProperty(prefix+"render_scene"));		
		if (properties.getProperty(prefix+"toRGB")!=null)                this.toRGB=Boolean.parseBoolean(properties.getProperty(prefix+"toRGB"));		
		if (properties.getProperty(prefix+"show_2d_correlations")!=null) this.show_2d_correlations=Boolean.parseBoolean(properties.getProperty(prefix+"show_2d_correlations"));		
		if (properties.getProperty(prefix+"show_motion_vectors")!=null)  this.show_motion_vectors=Boolean.parseBoolean(properties.getProperty(prefix+"show_motion_vectors"));		
		if (properties.getProperty(prefix+"debug_level")!=null)          this.debug_level=Integer.parseInt(properties.getProperty(prefix+"debug_level"));
		if (properties.getProperty(prefix+"test_ers")!=null)             this.test_ers=Boolean.parseBoolean(properties.getProperty(prefix+"test_ers"));		
		if (properties.getProperty(prefix+"test_ers0")!=null)            this.test_ers0=Integer.parseInt(properties.getProperty(prefix+"test_ers0"));
		if (properties.getProperty(prefix+"test_ers1")!=null)            this.test_ers1=Integer.parseInt(properties.getProperty(prefix+"test_ers1"));
	}
	
	@Override
	public IntersceneMatchParameters clone() throws CloneNotSupportedException {
		IntersceneMatchParameters imp =     new IntersceneMatchParameters();
 		imp.force_ref_dsi                 = this.force_ref_dsi;
		imp.force_orientations            = this.force_orientations;
//		imp.readjust_orient               = this.readjust_orient;
//		imp.force_interscene              = this.force_interscene;
		imp.min_num_orient                = this.min_num_orient;
		imp.min_num_interscene            = this.min_num_interscene;
		imp.generate_mapped               = this.generate_mapped;
		imp.reuse_video                   = this.reuse_video;
		imp.save_mapped_color             = this.save_mapped_color;
		imp.save_mapped_mono              = this.save_mapped_mono;
		imp.gen_avi_color                 = this.gen_avi_color;
		imp.gen_avi_mono                  = this.gen_avi_mono;
		imp.show_mapped_color             = this.show_mapped_color;
		imp.show_mapped_mono              = this.show_mapped_mono;
		
		imp.generate_raw                  = this.generate_raw;
		imp.generate_inf                  = this.generate_inf;
		imp.generate_fg                   = this.generate_fg;
		imp.generate_bg                   = this.generate_bg;
		
		imp.generate_stereo               = this.generate_stereo;
		
//		imp.stereo_bases                  = this.stereo_bases.clone();
		imp.stereo_views                  = this.stereo_views.clone();
		for (int i = 0; i < this.stereo_views.length; i++) {
			imp.stereo_views[i]           = this.stereo_views[i].clone();
		}
		imp.generate_stereo_var           = this.generate_stereo_var.clone();
		
		imp.export_images                 = this.export_images;
		imp.show_images                   = this.show_images;
		imp.show_images_bgfg              = this.show_images_bgfg;
		imp.show_images_mono              = this.show_images_mono;
		imp.export_ranges                 = this.export_ranges;
		imp.show_ranges                   = this.show_ranges;
		imp.export_ml_files               = this.export_ml_files;

		imp.show_color_nan                = this.show_color_nan;
		imp.show_mono_nan                 = this.show_mono_nan;
		imp.min_num_scenes                = this.min_num_scenes;
		imp.blur_egomotion                = this.blur_egomotion;
		imp.range_disparity_offset        = this.range_disparity_offset;
		imp.range_min_strength            = this.range_min_strength;
		imp.range_max                     = this.range_max;
		imp.num_bottom                    = this.num_bottom;
		imp.num_passes                    = this.num_passes;
		imp.max_change                    = this.max_change;
		imp.min_disparity                 = this.min_disparity;
		imp.max_sym_disparity             = this.max_sym_disparity;
		imp.min_strength_lma              = this.min_strength_lma;
		imp.min_strength_replace          = this.min_strength_replace;
		imp.min_strength_blur             = this.min_strength_blur;
		imp.sigma                         = this.sigma;
		imp.num_blur                      = this.num_blur;
		imp.disparity_corr                = this.disparity_corr;
		imp.outliers_nth_fromextrem       = this.outliers_nth_fromextrem;
		imp.outliers_tolerance_absolute   = this.outliers_tolerance_absolute;
		imp.outliers_tolerance_relative   = this.outliers_tolerance_relative;
		imp.outliers_max_iter             = this.outliers_max_iter;
		imp.outliers_max_strength2        = this.outliers_max_strength2;
		imp.outliers_nth_fromextrem2      = this.outliers_nth_fromextrem2;
		imp.outliers_tolerance_absolute2  = this.outliers_tolerance_absolute2;
		imp.outliers_tolerance_relative2  = this.outliers_tolerance_relative2;
		imp.outliers_lma_max_strength     = this.outliers_lma_max_strength;
		imp.outliers_max_strength         = this.outliers_max_strength;
		imp.outliers_from_lma_max_strength= this.outliers_from_lma_max_strength;
		imp.search_radius                 = this.search_radius;
		imp.remove_no_lma_neib            = this.remove_no_lma_neib;
		imp.diff_from_lma_pos             = this.diff_from_lma_pos;
		imp.diff_from_lma_neg             = this.diff_from_lma_neg;
		imp.outliers_lma_nth_fromextrem   = this.outliers_lma_nth_fromextrem;
		imp.filter_margin                 = this.filter_margin;

		imp.weak_tolerance_absolute       = this.weak_tolerance_absolute;
		imp.weak_tolerance_relative       = this.weak_tolerance_relative;
		imp.weak_min_neibs                = this.weak_min_neibs;
		imp.strong_strength               = this.strong_strength;
		imp.weak_strength                 = this.weak_strength;
		
		imp.stereo_merge                  = this.stereo_merge;
		imp.stereo_gap                    = this.stereo_gap;
		imp.stereo_intereye               = this. stereo_intereye;
		imp.stereo_phone_width            = this. stereo_phone_width;
		
		imp.extra_hor_tile        = this.extra_hor_tile;
		imp.extra_vert_tile       = this.extra_vert_tile;
		imp.crop_3d               = this.crop_3d;
		imp.sensor_mask           = this.sensor_mask;
		imp.merge_all             = this.merge_all;

		imp.video_fps =             this. video_fps;
		imp.sensor_fps =            this. sensor_fps;
		imp.mode_avi =              this. mode_avi;
		imp.avi_JPEG_quality =      this. avi_JPEG_quality;
		imp.run_ffmpeg =            this. run_ffmpeg;
		imp.video_ext =             this. video_ext;
		imp.video_codec =           this. video_codec;
		imp.video_crf =             this. video_crf;
		imp.remove_avi =            this. remove_avi;
		imp.video_codec_combo =     this. video_codec_combo;
		imp.video_crf_combo =       this. video_crf_combo;
		
		imp.um_mono =               this. um_mono;
		imp.um_sigma =              this. um_sigma;
		imp.um_weight =             this. um_weight;
		imp.mono_fixed =            this. mono_fixed;
		imp.mono_range =            this. mono_range;
		
		imp.anaglyth_en =           this. anaglyth_en;
		imp.anaglyph_left =         this. anaglyph_left;
		imp.anaglyph_right =        this. anaglyph_right;
		
		imp.annotate_color =        this. annotate_color;
		imp.annotate_mono =         this. annotate_mono;
		imp.annotate_color_color =  this. annotate_color_color;
		imp.annotate_color_mono =   this. annotate_color_mono;
		imp.annotate_transparent_mono = this. annotate_transparent_mono;

		imp.margin                = this.margin;
		imp.sensor_mask_inter     = this.sensor_mask_inter;
		imp.use_partial           = this.use_partial;
		imp.run_poly              = this.run_poly;
		imp.centroid_radius       = this.centroid_radius;
		imp.n_recenter            = this.n_recenter;
		
		imp.min_ref_str           = this.min_ref_str;

		imp.td_weight             = this.td_weight;
		imp.pd_weight             = this.pd_weight;
		imp.td_nopd_only          = this.td_nopd_only;

		imp.min_str               = this.min_str;
		imp.min_str_fpn           = this.min_str_fpn;
		imp.min_str_sum           = this.min_str_sum;
		imp.min_str_sum_fpn       = this.min_str_sum_fpn;
		imp.min_neibs             = this.min_neibs;
		imp.weight_zero_neibs     = this.weight_zero_neibs;
		imp.half_disparity        = this.half_disparity;
		imp.half_avg_diff         = this.half_avg_diff;
		
		imp.photo_en              = this.photo_en;
		imp.photo_num_full        = this.photo_num_full;
		imp.photo_num_refines     = this.photo_num_refines;
		imp.photo_min_strength    = this.photo_min_strength;
		imp.photo_max_diff        = this.photo_max_diff;
		imp.photo_debug           = this.photo_debug;
		
		imp.pix_step =              this.pix_step;
		imp.search_rad =            this.search_rad;
		imp.maybe_sum =             this.maybe_sum;
		imp.sure_sum =              this.sure_sum;
		imp.maybe_avg =             this.maybe_avg;
		imp.sure_avg =              this.sure_avg;
		
		imp.max_search_rms =        this.max_search_rms;
		imp.maybe_fom =             this.maybe_fom;
		imp.sure_fom =              this.sure_fom;
		
		imp.use_combo_dsi         = this.use_combo_dsi;
		imp.use_lma_dsi           = this.use_lma_dsi;
		
		imp.fpn_remove            = this.fpn_remove;
		imp.fpn_max_offset        = this.fpn_max_offset;
		imp.fpn_radius            = this.fpn_radius;
		imp.fpn_ignore_border     = this.fpn_ignore_border;
		
		imp.mov_en                = this.mov_en;
		imp.mov_sigma             = this.mov_sigma;
		imp.mov_max_std           = this.mov_max_std;
		imp.mov_thresh_rel        = this.mov_thresh_rel;
		imp.mov_thresh_abs        = this.mov_thresh_abs;
		imp.mov_clust_max         = this.mov_clust_max;
		imp.mov_grow              = this.mov_grow;
		imp.mov_show              = this.mov_show;
		imp.mov_debug_level       = this.mov_debug_level;
		
		imp.adjust_atr[0]         = this.adjust_atr[0];
		imp.adjust_atr[1]         = this.adjust_atr[1];
		imp.adjust_atr[2]         = this.adjust_atr[2];
		imp.adjust_xyz[0]         = this.adjust_xyz[0];
		imp.adjust_xyz[1]         = this.adjust_xyz[1];
		imp.adjust_xyz[2]         = this.adjust_xyz[2];
		imp.exit_change_atr       = this.exit_change_atr;
		imp.exit_change_xyz       = this.exit_change_xyz;
		imp.max_cycles            = this.max_cycles;
		imp.max_LMA               = this.max_LMA;
		imp.max_rms               = this.max_rms;
		imp.scene_is_ref_test     = this.scene_is_ref_test;
		imp.render_ref            = this.render_ref;
		imp.render_scene          = this.render_scene;
		imp.toRGB                 = this.toRGB;
		imp.show_2d_correlations  = this.show_2d_correlations;
		imp.show_motion_vectors   = this.show_motion_vectors ;
		imp.debug_level           = this.debug_level;
		return imp;
	}
	public static long getLongColor(Color color) {
		return ((long) color.getRGB()) & 0xffffffffL;
	}
	public static Color setLongColor(long lcolor) {
		if (lcolor < (1 << 24)) { // no alpha
			return new Color((int) lcolor);
		} else { // has alpha, may or may not fit into int
			if (lcolor > Integer.MAX_VALUE) {
				lcolor -= (1L << 32);
			}
			return new Color((int) lcolor, true);
		}
	}

	/*
	public void orderStereo(){
		boolean ordered;
		do {
			ordered=true;
			for (int i = 0; i < (stereo_bases.length - 1); i++) {
				if (stereo_bases[i+1]<stereo_bases[i]) {
					boolean en = generate_stereo_var[i+1];
					generate_stereo_var[i+1] = generate_stereo_var[i];
					generate_stereo_var[i] = en;
					double base =  stereo_bases[i+1];
					stereo_bases[i+1] = stereo_bases[i];
					stereo_bases[i] = base;
					ordered = false;
				}
			}
			
		} while (!ordered);
	}
	public void addStereo(double base, boolean en) {
		double [] bases = new double [stereo_bases.length + 1];
		boolean [] ens = new boolean [stereo_bases.length + 1];
		bases[0] = base;
		ens[0] = en;
		System.arraycopy(stereo_bases, 0, bases, 1,      stereo_bases.length);
		System.arraycopy(generate_stereo_var, 0, ens, 1, stereo_bases.length);
		stereo_bases = bases;
		generate_stereo_var = ens;
		orderStereo();
	}
	public void removeStereo(int indx) {
		if ((indx >=0) && (indx <stereo_bases.length)) {
			double [] bases = new double [stereo_bases.length - 1];
			boolean [] ens = new boolean [stereo_bases.length - 1];
			if (indx > 0) {
				System.arraycopy(stereo_bases,        0, bases, 0, indx);
				System.arraycopy(generate_stereo_var, 0, ens,   0, indx);
			}
			if (indx < (stereo_bases.length - 1)) {
				System.arraycopy(stereo_bases,        indx+1, bases, indx, stereo_bases.length - indx - 1);
				System.arraycopy(generate_stereo_var, indx+1, ens,   indx, stereo_bases.length - indx - 1);
			}
			stereo_bases = bases;
			generate_stereo_var = ens;
		}
	}
	

	 */
	
	public void orderStereoViews(){
		boolean ordered;
		do {
			ordered=true;
			for (int i = 0; i < (stereo_views.length - 1); i++) {
				if (stereo_views[i+1][0] > stereo_views[i][0]) {
					continue;
				}
				if (    (stereo_views[i+1][0] == stereo_views[i][0]) &&
						(stereo_views[i+1][1] >  stereo_views[i][1])) {
					continue;
				}
				if (    (stereo_views[i+1][0] == stereo_views[i][0]) &&
						(stereo_views[i+1][1] == stereo_views[i][1]) &&
						(stereo_views[i+1][2] >  stereo_views[i][2])) {
					continue;
				}
				if (    (stereo_views[i+1][0] == stereo_views[i][0]) &&
						(stereo_views[i+1][1] == stereo_views[i][1]) &&
						(stereo_views[i+1][2] == stereo_views[i][2])) {
					// all same values - remove extra
					generate_stereo_var[i] |= generate_stereo_var[i+1];
					for (int j = i+1; j < (stereo_views.length - 1); j++) {
						generate_stereo_var[j] = generate_stereo_var[j+1];
						stereo_views[j] = stereo_views[j + 1]; 
					}
					ordered = false;
					break; // next while
				}
				boolean en = generate_stereo_var[i+1];
				generate_stereo_var[i+1] = generate_stereo_var[i];
				generate_stereo_var[i] = en;
				double [] view =  stereo_views[i+1];
				stereo_views[i+1] = stereo_views[i];
				stereo_views[i] = view;
				ordered = false;
			}

		} while (!ordered);
		return;
	}

	public void addStereoView(String stereo_view_string, boolean en) {
		double[] stereo_view = StringToDoubles(stereo_view_string,3);
		if (stereo_view != null) {
			addStereoView(stereo_view, en);
		}
	}
	
	public void addStereoView(double[] stereo_view, boolean en) {
		double [][]  views = new double [stereo_views.length + 1][];
		boolean [] ens = new boolean [stereo_views.length + 1];
		views[0] = stereo_view;
		ens[0] = en;
		System.arraycopy(stereo_views, 0, views, 1,      stereo_views.length);
		System.arraycopy(generate_stereo_var, 0, ens, 1, stereo_views.length);
		stereo_views = views;
		generate_stereo_var = ens;
		orderStereoViews();
	}
	
	public void removeStereoView(int indx) {
		if ((indx >=0) && (indx <stereo_views.length)) {
			double [][] views = new double [stereo_views.length - 1][];
			boolean [] ens = new boolean [stereo_views.length - 1];
			if (indx > 0) {
				System.arraycopy(stereo_views,        0, views, 0, indx);
				System.arraycopy(generate_stereo_var, 0, ens,   0, indx);
			}
			if (indx < (stereo_views.length - 1)) {
				System.arraycopy(stereo_views,        indx+1, views, indx, stereo_views.length - indx - 1);
				System.arraycopy(generate_stereo_var, indx+1, ens,   indx, stereo_views.length - indx - 1);
			}
			stereo_views = views;
			generate_stereo_var = ens;
		}
	}

	public static String doublesToString(double [] data) {
		return doublesToString(data, null);
	}
	public static String doublesToString(double [] data, String fmt) {
//		if ((fmt == null) || (fmt.trim().length()==0)) {
//			fmt = "%.0f";
//		}
		String s = "";
		for (int i = 0; i < data.length; i++) {
			if (fmt==null) {
				s += data[i]; // unformatted
			} else { 
				s+=String.format(fmt,data[i]);
			}
			if (i < (data.length - 1)) {
				s+= ", ";
			}
		}
		return s;
	}
	
	public static double [] StringToDoubles(String s, int len) {
		StringTokenizer st = new StringTokenizer(s, " \t\n\r\f,");
		if (st.countTokens() == 0) {
			return null;
		}
		if (len <= 0) {
			len = st.countTokens();
		}
		double [] data = new double [len]; 
		int i = 0;
		while (st.hasMoreTokens() && (i < len)) {
			double d = 0;
			try {
				d = Double.parseDouble(st.nextToken());
			} catch(NumberFormatException e){
				d = 0;
			}

			data[i++] = d;
		}
		return data;
	}
	
		
	
}
