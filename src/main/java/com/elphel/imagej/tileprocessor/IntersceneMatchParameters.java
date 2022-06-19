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

import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;

public class IntersceneMatchParameters {
	// Maybe add parameters to make sure there is enough data? Enough in each zone? Enough spread?
	public  boolean force_ref_dsi =      false; // true;
	public  boolean force_orientations = false;
	public  boolean force_interscene =   false; // true;
	public  boolean export_images =      true;  // pseudo-color 16-slice images (same disparity, COMBO_DSN_INDX_DISP_FG and COMBO_DSN_INDX_DISP_BG_ALL,
	public  boolean show_images =        false; // color, infinity
	public  boolean show_images_bgfg =   false; // bg and fg
	public  boolean show_images_mono =   false; // float, monochrome 16-slice images (same disparity, COMBO_DSN_INDX_DISP_FG and COMBO_DSN_INDX_DISP_BG_ALL,
	public  boolean show_color_nan =     true;  // use NAN background for color images (sharp, but distinct black)
	public  boolean show_mono_nan =      false; // use NAN background for monochrome images (sharp, but distinct black)
	public  boolean show_ranges =        true;
	
	public  boolean generate_mapped =    true;
	public  int     extra_hor_tile =     15;
	public  int     extra_vert_tile =    10;
	public  int     sensor_mask =        1; // -1 - all
	public  int     mode3d =             1; // -1 - raw, 0 - infinity, 1 - FG, 2 - BG
	public  boolean show_mapped_color =  true;
	public  boolean show_mapped_mono =   false;
	
	public  double  range_disparity_offset =   -0.08;
	public  double  range_min_strength = 0.5;
	public  double  range_max =       5000.0;
	
// Other parameters for filtering depth maps	
	public  int         num_bottom =                      6; // average this number of lowest disparity neighbors (of 8)
	public  int         num_passes =                    100;
	public  double      max_change =                      1.0e-3;
	public  double      min_disparity =                  -0.2;
	public  double      max_sym_disparity =               0.2;		
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
	public  double      outliers_lma_max_strength =       0.5;
	public  double      outliers_max_strength =           0.1; ///  0.25;
	public  double      outliers_from_lma_max_strength =  0.8;
	public  int         search_radius =                   3;  // Search farther if no LMA neighbor is found closer. Original value - 1 (8 neighbors)
	public  boolean     remove_no_lma_neib =              true;
	public  double      diff_from_lma_pos =             100.0;
	public  double      diff_from_lma_neg =               2.0;
	public  int         outliers_lma_nth_fromextrem =     1; 
	public  int         filter_margin =                  -8; // 8; // pixels
	
	

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
	public  boolean td_nopd_only =                  true;  // only use TD accumulated data if no safe PD is available for the tile.

	public  double  min_str =                       0.25;  // minimal correlation strength for all but TD-accumulated layer
	public  double  min_str_sum =                   0.8;	  // minimal correlation strength for TD-accumulated layer
	public  int     min_neibs =                     2;	   // minimal number of strong neighbors (> min_str)
	public  double  weight_zero_neibs =             0.2;   // Reduce weight for no-neib (1.0 for all 8)
	public  double  half_disparity =                5.0;   // Reduce weight twice for this disparity
	public  double  half_avg_diff =                 0.2;   // when L2 of x,y difference from average of neibs - reduce twice

	// Detect initial match
	public  int     pix_step =                      4;     // Azimuth/tilt search step in pixels
	public  int     search_rad =                   10;     // Search radius in steps
	public  double  maybe_sum =                     8.0;   // minimal sum of strengths (will search for the best)
	public  double  shure_sum =                    30.0;   // definitely good sum of strengths (no farther search)
	public  double  maybe_avg =                     0.005; // maybe average strength
	public  double  shure_avg =                     0.015; // sure average strength
	
	
	//LMA parameters
	public  boolean [] adjust_atr = new boolean [] {true,true,true};
	public  boolean [] adjust_xyz = new boolean [] {true,true,true};
	public  double  exit_change_atr =               1.0E-5;//rad,  L2 norm for difference to ext LMA 
	public  double  exit_change_xyz =               1.0E-3;//meters, L2 norm for difference to ext LMA 
	public  int     max_cycles =                   10;     // hard limit on full correlation/LMA cycles
	public  int     max_LMA =                      25;     // hard limit on LMA iterations
	public  double  max_rms =                       2.0;   // maximal RMS to consider adjustment to be a failure
	
	// Debug and visualisation
	public  boolean scene_is_ref_test=              false; // correlate ref-2-ref for testing
	private boolean render_ref =                    true;  // render reference scene
	private boolean render_scene =                  true;  // render scene to be correlated with ref
	public  boolean toRGB =                         true;  // render scenes in pseudo-colors DOES NOT WORK - use ColorProcParameters
	private boolean show_2d_correlations =          true;  // show raw 2D correlations (individual and combined)
	private boolean show_motion_vectors =           true;  // show calculated motion vectors
	public  int     debug_level =                     -1;  // all renders are disable for debug_level < 0, scene "renders" for for debug_level < 1
    
	public boolean renderRef()          {return (debug_level>1) && render_ref;}
	public boolean renderScene()        {return (debug_level>1) && render_scene;}
	public boolean show2dCorrelations() {return (debug_level>1) && show_2d_correlations;}
	public boolean showMotionVectors()  {return (debug_level>1) && show_motion_vectors;}
	public boolean showCorrMotion()     {return (debug_level>0) && show_motion_vectors;}

	public IntersceneMatchParameters() {
		
	}

	public void dialogQuestions(GenericJTabbedDialog gd) {
		//		gd.addMessage  ("Scene parameters selection");

		gd.addMessage  ("Build series options");
		gd.addCheckbox ("Force reference scene DSI calculation",     this.force_ref_dsi,
				"Calculate reference scene DSI even if the file exists.");
		gd.addCheckbox ("Force egomotion calculation",               this.force_orientations,
				"Calculate relative poses of each scene camera relative to the reference scene even if the data exists.");
		gd.addCheckbox ("Force interscene DSI accumulation",         this.force_interscene,
				"Force interscene calculation (+ML export) even if it was performed before.");
		gd.addCheckbox ("Export all-sensor images",                  this.export_images,
				"Export multi-slice images: with constant disparity, with foreground disparity, and with background disparity");
		gd.addCheckbox ("Show exported images (same disparity)",     this.show_images,
				"Display generated/saved image set, pseudocolors");
		gd.addCheckbox ("Show exported FG/BG",                       this.show_images_bgfg,
				"Show foreground and background exported images");
		gd.addCheckbox ("Show floating-point monochrome images",     this.show_images_mono,
				"Display generated/saved monochrome images");
		gd.addCheckbox ("Color NaN background",                      this.show_color_nan,
				"Use NaN for undefined tiles (false - 0.0f). NaN produces sharp distinct result, 0.0f - blended");
		gd.addCheckbox ("Mono NaN background",                       this.show_mono_nan,
				"Use NaN for undefined tiles (false - 0.0f). NaN produces sharp distinct result, 0.0f - blended");

		gd.addMessage  ("Metric distance map generation");
		gd.addCheckbox ("Show distances in meters",                  this.show_ranges,
				"Calculate strength, distance, X, and Y in meters");
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
		gd.addNumericField("Max neib. disparity difference",         this.max_sym_disparity, 5,7,"pix",
				"Difference from the nearest of neighbors should not exceed this");
		gd.addNumericField("Minimal strength to replace",            this.min_strength_replace, 5,7,"",
				"Minimal strength to replace");
		gd.addNumericField("Minimal strength to blur",               this.min_strength_blur, 5,7,"",
				"Minimal strength to blur");
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
		gd.addNumericField("Margin from the sensor FOV edge ",        this.filter_margin, 0,3,"pix",
				"Disregard tiles with centers closer to the sensor FoV in pixels");
		
		
		gd.addMessage  ("Show scene sequence");
		gd.addCheckbox ("Generate mapped scene sequence",            this.generate_mapped,
				"Generate scene sequence mapped to the reference scene");
		gd.addNumericField("Scene sequence horizontal extra",        this.extra_hor_tile, 0,3,"tiles",
				"Enlarge reference scene window horizontally in each direction to accommodate other scenes in a sequence");
		gd.addNumericField("Scene sequence vertical extra",          this.extra_vert_tile, 0,3,"tiles",
				"Enlarge reference scene window vertically in each direction to accommodate other scenes in a sequence");
		gd.addNumericField("Sensor mask (bitmask, -1 - all sensors)",this.sensor_mask, 0,3,"",
				"Select which sensors to be included in each scene of the sequence");
		gd.addNumericField("3D mode (-1 - RAW, 0 - infinity, 1 - FG, 2 BG)",   this.mode3d, 0, 3, "",
				"3D mode for rendering scenes in a sequence: -1 - raw images, 0 - no 3D, use infinity; 1 - Foreground; 2 - Background");
		gd.addCheckbox ("Show scene sequences in (pseudo)colors",    this.show_mapped_color,
				"Show generated scene sequences in (pseudo)color mode");
		gd.addCheckbox ("Show scene sequences in monochrome",        this.show_mapped_mono,
				"Show generated scene sequences in monochrome mode");
		
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

		gd.addMessage  ("Mixing TD and PD accumulation of 2d correlations");
		gd.addNumericField("TD-accumulated weight",                  this.td_weight, 5,7,"",
				"Mix argmax from TD-accumulated correlation.");
		gd.addNumericField("PD-accumulated weight",                  this.pd_weight, 5,7,"",
				"Mix argmax from PD-accumulated correlation.");
		gd.addCheckbox ("TD when no PD only",                        this.td_nopd_only,
				"Use argmax from TD only if PD data is not available for this tile.");
		
		gd.addMessage  ("Filtering motion vectors");
		gd.addNumericField("Minimal correlation strength (non-sum)", this.min_str, 5,7,"",
				"Minimal correlation strength for individual correlation and for pixel-domain averaged one. Weeker tiles results are removed.");
		gd.addNumericField("Minimal correlation strength (sum only)",this.min_str_sum, 5,7,"",
				"Minimal correlation strength for transform-domain averaging. Weeker tiles results are removed.");
		gd.addNumericField("Minimal number of neighbors (of 8)",     this.min_neibs, 0,3,"",
				"Remove motion vectors with less than this number of defined (passing min_str) neighbors.");
		gd.addNumericField("No-neighbors weight (<1.0)",             this.weight_zero_neibs, 5,7,"",
				"Punish for reduced neighbors - weigh for no-neighbors), weight of 8 neighbors = 1.0.");
		gd.addNumericField("Disparity to reduce weight twice from infinity", this.half_disparity, 5,7,"",
				"Weight at this disparity is 0.5, at infinity - 1.0.");
		gd.addNumericField("Difference from neighbors average ",     this.half_avg_diff, 5,7,"",
				"Reduce twice for high difference from neighbors average.");

		gd.addMessage  ("Initial search for the inter-scene match");
		gd.addNumericField("Azimuth/tilt step",                      this.pix_step, 0,3,"pix",
				"Search in a spiral starting with no-shift with this step between probes, in approximate pixels");
		gd.addNumericField("Search spiral radius",                   this.search_rad, 0,3,"steps",
				"Maximal search radius in steps");
		gd.addNumericField("\"Maybe\" sum of strengths",             this.maybe_sum, 5,7,"",
				"Minimal acceptable sum of defined tiles strengths (will look for the best among matching)");
		gd.addNumericField("\"Sure\" sum of strengths",              this.shure_sum, 5,7,"",
				"Definitely sufficient sum of defined tiles strengths (will non continue looking for better).");
		gd.addNumericField("\"Maybe\" average of strengths",         this.maybe_avg, 5,7,"",
				"Minimal acceptable average of defined tiles strengths (will look for the best among matching)");
		gd.addNumericField("\"Sure\" average of strengths",          this.shure_avg, 5,7,"",
				"Definitely sufficient average of defined tiles strengths (will non continue looking for better).");
		
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

	}

	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.force_ref_dsi =            gd.getNextBoolean();
		this.force_orientations =       gd.getNextBoolean();
		this.force_interscene =         gd.getNextBoolean();
		this.export_images =            gd.getNextBoolean();
		this.show_images =              gd.getNextBoolean();
		this.show_images_bgfg =         gd.getNextBoolean();
		this.show_images_mono =         gd.getNextBoolean();
		this.show_color_nan =           gd.getNextBoolean();
		this.show_mono_nan =            gd.getNextBoolean();
		this.show_ranges =              gd.getNextBoolean();
		this.range_disparity_offset =   gd.getNextNumber();
		this.range_min_strength =       gd.getNextNumber();
		this.range_max =                gd.getNextNumber();

		this.num_bottom =               (int) gd.getNextNumber();
		this.num_passes =               (int) gd.getNextNumber();
		this.max_change =                     gd.getNextNumber();
		this.min_disparity =                  gd.getNextNumber();
		this.max_sym_disparity =              gd.getNextNumber();
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
		
		this.generate_mapped =          gd.getNextBoolean();
		this.extra_hor_tile =     (int) gd.getNextNumber();
		this.extra_vert_tile =    (int) gd.getNextNumber();
		this.sensor_mask =        (int) gd.getNextNumber();
		this.mode3d =             (int) gd.getNextNumber();
		this.show_mapped_color =        gd.getNextBoolean();
		this.show_mapped_mono =         gd.getNextBoolean();
		
		this.margin =             (int) gd.getNextNumber();
		this.sensor_mask_inter=   (int) gd.getNextNumber();
		this.use_partial =              gd.getNextBoolean();
		this.run_poly =                 gd.getNextBoolean();
		this.centroid_radius =          gd.getNextNumber();
		this.n_recenter =         (int) gd.getNextNumber();

		this.td_weight =                gd.getNextNumber();
		this.pd_weight =                gd.getNextNumber();
		this.td_nopd_only =             gd.getNextBoolean();
		
		this.min_str =                  gd.getNextNumber();
		this.min_str_sum =              gd.getNextNumber();
		this.min_neibs =          (int) gd.getNextNumber();
		this.weight_zero_neibs =        gd.getNextNumber();
		this.half_disparity =           gd.getNextNumber();
		this.half_avg_diff =            gd.getNextNumber();
		this.pix_step =           (int) gd.getNextNumber();
		this.search_rad =         (int) gd.getNextNumber();
		this.maybe_sum =                gd.getNextNumber();
		this.shure_sum =                gd.getNextNumber();
		this.maybe_avg =                gd.getNextNumber();
		this.shure_avg =                gd.getNextNumber();
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
		
		if (this.weight_zero_neibs > 1.0) this.weight_zero_neibs = 1.0;
	}
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"force_ref_dsi",        this.force_ref_dsi + "");     // boolean
		properties.setProperty(prefix+"force_orientations",   this.force_orientations + "");// boolean
		properties.setProperty(prefix+"force_interscene",     this.force_interscene + "");  // boolean
		properties.setProperty(prefix+"export_images",        this.export_images + "");     // boolean
		properties.setProperty(prefix+"show_images",          this.show_images + "");       // boolean
		properties.setProperty(prefix+"show_images_bgfg",     this.show_images_bgfg + "");  // boolean
		properties.setProperty(prefix+"show_images_mono",     this.show_images_mono + "");  // boolean
		properties.setProperty(prefix+"show_color_nan",       this.show_color_nan + "");    // boolean
		properties.setProperty(prefix+"show_mono_nan",        this.show_mono_nan + "");     // boolean
		properties.setProperty(prefix+"show_ranges",          this.show_ranges + "");       // boolean
		properties.setProperty(prefix+"range_disparity_offset",this.range_disparity_offset+""); // double
		properties.setProperty(prefix+"range_min_strength",   this.range_min_strength+"");  // double
		properties.setProperty(prefix+"range_max",            this.range_max+"");           // double

		properties.setProperty(prefix+"num_bottom",                    this.num_bottom+"");                    // int
		properties.setProperty(prefix+"num_passes",                    this.num_passes+"");                    // int
		properties.setProperty(prefix+"max_change",                    this.max_change+"");                    // double
		properties.setProperty(prefix+"min_disparity",                 this.min_disparity+"");                 // double
		properties.setProperty(prefix+"max_sym_disparity",             this.max_sym_disparity+"");             // double
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
		
		properties.setProperty(prefix+"generate_mapped",      this.generate_mapped+"");     // boolean
		properties.setProperty(prefix+"extra_hor_tile",       this.extra_hor_tile+"");      // int
		properties.setProperty(prefix+"extra_vert_tile",      this.extra_vert_tile+"");     // int
		properties.setProperty(prefix+"sensor_mask",          this.sensor_mask+"");         // int
		properties.setProperty(prefix+"mode3d",                this.mode3d+"");         // int
		properties.setProperty(prefix+"show_mapped_color",    this.show_mapped_color+"");   // boolean
		properties.setProperty(prefix+"show_mapped_mono",     this.show_mapped_mono+"");    // boolean
		
		properties.setProperty(prefix+"margin",               this.margin+"");              // int
		properties.setProperty(prefix+"sensor_mask_inter",    this.sensor_mask_inter+"");   // int
		properties.setProperty(prefix+"use_partial",          this.use_partial+"");         // boolean
		properties.setProperty(prefix+"run_poly",             this.run_poly+"");            // boolean
		properties.setProperty(prefix+"centroid_radius",      this.centroid_radius+"");     // double
		properties.setProperty(prefix+"n_recenter",           this.n_recenter+"");          // int
		
		properties.setProperty(prefix+"td_weight",            this.td_weight+"");           // double
		properties.setProperty(prefix+"pd_weight",            this.pd_weight+"");           // double
		properties.setProperty(prefix+"td_nopd_only",         this.td_nopd_only+"");        // boolean
		
		properties.setProperty(prefix+"min_str",              this.min_str+"");             // double
		properties.setProperty(prefix+"min_str_sum",          this.min_str_sum+"");         // double
		properties.setProperty(prefix+"min_neibs",            this.min_neibs+"");           // int
		properties.setProperty(prefix+"weight_zero_neibs",    this.weight_zero_neibs+"");   // double
		properties.setProperty(prefix+"half_disparity",       this.half_disparity+"");      // double
		properties.setProperty(prefix+"half_avg_diff",        this.half_avg_diff+"");       // double
		properties.setProperty(prefix+"pix_step",             this.pix_step+"");            // int
		properties.setProperty(prefix+"search_rad",           this.search_rad+"");          // int
		properties.setProperty(prefix+"maybe_sum",            this.maybe_sum+"");           // double
		properties.setProperty(prefix+"shure_sum",            this.shure_sum+"");           // double
		properties.setProperty(prefix+"maybe_avg",            this.maybe_avg+"");           // double
		properties.setProperty(prefix+"shure_avg",            this.shure_avg+"");           // double
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
	}
	
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"force_ref_dsi")!=null)        this.force_ref_dsi=Boolean.parseBoolean(properties.getProperty(prefix+"force_ref_dsi"));		
		if (properties.getProperty(prefix+"force_orientations")!=null)   this.force_orientations=Boolean.parseBoolean(properties.getProperty(prefix+"force_orientations"));		
		if (properties.getProperty(prefix+"force_interscene")!=null)     this.force_interscene=Boolean.parseBoolean(properties.getProperty(prefix+"force_interscene"));		
		if (properties.getProperty(prefix+"export_images")!=null)        this.export_images=Boolean.parseBoolean(properties.getProperty(prefix+"export_images"));		
		if (properties.getProperty(prefix+"show_images")!=null)          this.show_images=Boolean.parseBoolean(properties.getProperty(prefix+"show_images"));		
		if (properties.getProperty(prefix+"show_images_bgfg")!=null)     this.show_images_bgfg=Boolean.parseBoolean(properties.getProperty(prefix+"show_images_bgfg"));		
		if (properties.getProperty(prefix+"show_images_mono")!=null)     this.show_images_mono=Boolean.parseBoolean(properties.getProperty(prefix+"show_images_mono"));		
		if (properties.getProperty(prefix+"show_color_nan")!=null)       this.show_color_nan=Boolean.parseBoolean(properties.getProperty(prefix+"show_color_nan"));		
		if (properties.getProperty(prefix+"show_mono_nan")!=null)        this.show_mono_nan=Boolean.parseBoolean(properties.getProperty(prefix+"show_mono_nan"));		
		if (properties.getProperty(prefix+"show_ranges")!=null)          this.show_images=Boolean.parseBoolean(properties.getProperty(prefix+"show_ranges"));
		if (properties.getProperty(prefix+"range_disparity_offset")!=null) this.range_disparity_offset=Double.parseDouble(properties.getProperty(prefix+"range_disparity_offset"));
		if (properties.getProperty(prefix+"range_min_strength")!=null)   this.range_min_strength=Double.parseDouble(properties.getProperty(prefix+"range_min_strength"));
		if (properties.getProperty(prefix+"range_max")!=null)            this.range_max=Double.parseDouble(properties.getProperty(prefix+"range_max"));
		
		if (properties.getProperty(prefix+"num_bottom")!=null)                    this.num_bottom=Integer.parseInt(properties.getProperty(prefix+"num_bottom"));
		if (properties.getProperty(prefix+"num_passes")!=null)                    this.num_passes=Integer.parseInt(properties.getProperty(prefix+"num_passes"));
		if (properties.getProperty(prefix+"max_change")!=null)                    this.max_change=Double.parseDouble(properties.getProperty(prefix+"max_change"));
		if (properties.getProperty(prefix+"min_disparity")!=null)                 this.min_disparity=Double.parseDouble(properties.getProperty(prefix+"min_disparity"));
		if (properties.getProperty(prefix+"max_sym_disparity")!=null)             this.max_sym_disparity=Double.parseDouble(properties.getProperty(prefix+"max_sym_disparity"));
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
		
		if (properties.getProperty(prefix+"generate_mapped")!=null)      this.generate_mapped=Boolean.parseBoolean(properties.getProperty(prefix+"generate_mapped"));
		if (properties.getProperty(prefix+"extra_hor_tile")!=null)       this.extra_hor_tile=Integer.parseInt(properties.getProperty(prefix+"extra_hor_tile"));
		if (properties.getProperty(prefix+"extra_vert_tile")!=null)      this.extra_vert_tile=Integer.parseInt(properties.getProperty(prefix+"extra_vert_tile"));
		if (properties.getProperty(prefix+"sensor_mask")!=null)          this.sensor_mask=Integer.parseInt(properties.getProperty(prefix+"sensor_mask"));
		if (properties.getProperty(prefix+"mode3d")!=null)               this.mode3d=Integer.parseInt(properties.getProperty(prefix+"mode3d"));
		if (properties.getProperty(prefix+"show_mapped_color")!=null)    this.show_mapped_color=Boolean.parseBoolean(properties.getProperty(prefix+"show_mapped_color"));
		if (properties.getProperty(prefix+"show_mapped_mono")!=null)     this.show_mapped_mono=Boolean.parseBoolean(properties.getProperty(prefix+"show_mapped_mono"));
		if (properties.getProperty(prefix+"margin")!=null)               this.margin=Integer.parseInt(properties.getProperty(prefix+"margin"));
		if (properties.getProperty(prefix+"sensor_mask_inter")!=null)    this.sensor_mask_inter=Integer.parseInt(properties.getProperty(prefix+"sensor_mask_inter"));
		if (properties.getProperty(prefix+"use_partial")!=null)          this.use_partial=Boolean.parseBoolean(properties.getProperty(prefix+"use_partial"));		
		if (properties.getProperty(prefix+"run_poly")!=null)             this.run_poly=Boolean.parseBoolean(properties.getProperty(prefix+"run_poly"));
		if (properties.getProperty(prefix+"centroid_radius")!=null)      this.centroid_radius=Double.parseDouble(properties.getProperty(prefix+"centroid_radius"));
		if (properties.getProperty(prefix+"n_recenter")!=null)           this.n_recenter=Integer.parseInt(properties.getProperty(prefix+"n_recenter"));
		
		if (properties.getProperty(prefix+"td_weight")!=null)            this.td_weight=Double.parseDouble(properties.getProperty(prefix+"td_weight"));
		if (properties.getProperty(prefix+"pd_weight")!=null)            this.pd_weight=Double.parseDouble(properties.getProperty(prefix+"pd_weight"));
		if (properties.getProperty(prefix+"td_nopd_only")!=null)         this.td_nopd_only=Boolean.parseBoolean(properties.getProperty(prefix+"td_nopd_only"));
		
		if (properties.getProperty(prefix+"min_str")!=null)              this.min_str=Double.parseDouble(properties.getProperty(prefix+"min_str"));
		if (properties.getProperty(prefix+"min_str_sum")!=null)          this.min_str_sum=Double.parseDouble(properties.getProperty(prefix+"min_str_sum"));
		if (properties.getProperty(prefix+"min_neibs")!=null)            this.min_neibs=Integer.parseInt(properties.getProperty(prefix+"min_neibs"));
		if (properties.getProperty(prefix+"weight_zero_neibs")!=null)    this.weight_zero_neibs=Double.parseDouble(properties.getProperty(prefix+"weight_zero_neibs"));
		if (properties.getProperty(prefix+"half_disparity")!=null)       this.half_disparity=Double.parseDouble(properties.getProperty(prefix+"half_disparity"));
		if (properties.getProperty(prefix+"half_avg_diff")!=null)        this.half_avg_diff=Double.parseDouble(properties.getProperty(prefix+"half_avg_diff"));
		if (properties.getProperty(prefix+"pix_step")!=null)             this.pix_step=Integer.parseInt(properties.getProperty(prefix+"pix_step"));
		if (properties.getProperty(prefix+"search_rad")!=null)           this.search_rad=Integer.parseInt(properties.getProperty(prefix+"search_rad"));
		if (properties.getProperty(prefix+"maybe_sum")!=null)            this.maybe_sum=Double.parseDouble(properties.getProperty(prefix+"maybe_sum"));
		if (properties.getProperty(prefix+"shure_sum")!=null)            this.shure_sum=Double.parseDouble(properties.getProperty(prefix+"shure_sum"));
		if (properties.getProperty(prefix+"maybe_avg")!=null)            this.maybe_avg=Double.parseDouble(properties.getProperty(prefix+"maybe_avg"));
		if (properties.getProperty(prefix+"shure_avg")!=null)            this.shure_avg=Double.parseDouble(properties.getProperty(prefix+"shure_avg"));
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
	}
	
	@Override
	public IntersceneMatchParameters clone() throws CloneNotSupportedException {
		IntersceneMatchParameters imp =     new IntersceneMatchParameters();
 		imp.force_ref_dsi         = this.force_ref_dsi;
		imp.force_orientations    = this.force_orientations;
		imp.force_interscene      = this.force_interscene;
		imp.export_images         = this.export_images;
		imp.show_images           = this.show_images;
		imp.show_images_bgfg      = this.show_images_bgfg;
		imp.show_images_mono      = this.show_images_mono;
		imp.show_color_nan        = this.show_color_nan;
		imp.show_mono_nan         = this.show_mono_nan;
		imp.show_ranges           = this.show_ranges;
		imp.range_disparity_offset= this.range_disparity_offset;
		imp.range_min_strength    = this.range_min_strength;
		imp.range_max             = this.range_max;
		
		imp.num_bottom                    = this.num_bottom;
		imp.num_passes                    = this.num_passes;
		imp.max_change                    = this.max_change;
		imp.min_disparity                 = this.min_disparity;
		imp.max_sym_disparity             = this.max_sym_disparity;
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

		imp.generate_mapped       = this.generate_mapped;
		imp.extra_hor_tile        = this.extra_hor_tile;
		imp.extra_vert_tile       = this.extra_vert_tile;
		imp.sensor_mask           = this.sensor_mask;
		imp.mode3d                = this.mode3d;               
		imp.show_mapped_color     = this.show_mapped_color;
		imp.show_mapped_mono      = this.show_mapped_mono;
		
		imp.margin                = this.margin;
		imp.sensor_mask_inter     = this.sensor_mask_inter;
		imp.use_partial           = this.use_partial;
		imp.run_poly              = this.run_poly;
		imp.centroid_radius       = this.centroid_radius;
		imp.n_recenter            = this.n_recenter;
		
		imp.td_weight             = this.td_weight;
		imp.pd_weight             = this.pd_weight;
		imp.td_nopd_only          = this.td_nopd_only;

		imp.min_str               = this.min_str;
		imp.min_str_sum           = this.min_str_sum;
		imp.min_neibs             = this.min_neibs;
		imp.weight_zero_neibs     = this.weight_zero_neibs;
		imp.half_disparity        = this.half_disparity;
		imp.half_avg_diff         = this.half_avg_diff;
		imp.pix_step =              this.pix_step;
		imp.search_rad =            this.search_rad;
		imp.maybe_sum =             this.maybe_sum;
		imp.shure_sum =             this.shure_sum;
		imp.maybe_avg =             this.maybe_avg;
		imp.shure_avg =             this.shure_avg;
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
}
