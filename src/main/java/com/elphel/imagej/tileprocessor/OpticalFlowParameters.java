package com.elphel.imagej.tileprocessor;
/**
 **
 ** OpticalFlowParameters - parameters defining Optical Flow operations
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  OpticalFlowParameters.java is free software: you can redistribute it and/or modify
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

public class OpticalFlowParameters {
    public double      scale_no_lma_disparity = 1.0; // multiply strength were disparity_lma = NaN; 
	public double      k_prev = 0.75;
	public double      ers_to_pose_scale = 0.75;
	public double      tolerance_absolute_ref = 0.25; // absolute disparity half-range in each tile
	public double      tolerance_relative_ref = 0.2; // relative disparity half-range in each tile
	public double      center_occupancy_ref =   0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)
	public int         num_laplassian = 100;
	public double      change_laplassian = 0.005 ;
	public double      tolerance_absolute_inter = 0.25; // absolute disparity half-range in each tile
	public double      tolerance_relative_inter = 0.2; // relative disparity half-range in each tile
	public double      occupancy_inter =          0.25;   // fraction of remaining  tiles in the center 8x8 area (<1.0)
	public double      nsigma =                   1.5; // Remove outliers by more that this scaled sigma
	public double      nsigma2 =                  2.0; // Second time  outlier filter (<0 - disable)
	public double []   chn_weights = {1.0,1.0,1.0,1.0}; // strength, r,b,g
//	double []   chn_weights = {1.0,0.0,0.0,0.0}; // strength, r,b,g
//	double []   chn_weights = {0.0,1.0,1.0,1.0}; // strength, r,b,g
	public double      corr_sigma =               0.5;
	public double      fat_zero =                 0.05;
	public double      frac_radius =              0.9;  // add to integer radius for window calculation
	public double      tolerance_absolute_macro = 0.25; // absolute disparity half-range to consolidate macro tiles
	public double      tolerance_relative_macro = 0.2;  // relative disparity half-range to consolidate macro tiles
	public int         iradius_cm =               3;      // half-size of the square to process 
	public double      dradius_cm =               1.5;      // weight calculation (1/(r/dradius)^2 + 1)
	public int         refine_num_cm =            5;   // number of iterations to apply weights around new center
	public double      magic_scale = 0.85; // 2.0 * 0.85;
	public int         max_refines =             50;
	public int         num_refine_all =           3;
	public double      min_change =               0.1; // 01;//   sqrt (dx*dx + dy*dy) for correction (int tiles) in pixels
	public int         best_neibs_num =           4; // use  4 best neighbors to calculate std deviation
	public double      ref_stdev =                5.0; // strength 0.5 if standard deviation of best neighbors to tile difference is this.
	public boolean     ignore_ers =               false; // ignore velocities from individual ERS (LWIR - ignore, RGB - do not ignore
	public double      lpf_pairs =                5.0;  // velocities LPF during pairwise fitting
	public double      lpf_series =               5.0;  // velocities LPF during all-to-reference fitting

	public boolean     pattern_mode =             false; // true; // not yet used
	public double      max_rms_maybe =            7.0;
	public double      max_rms_sure =             2.1; // too high - false positives (for selected set >= 2.2), too low - longer testing
	public int         pix_step =                20;
	public int         search_rad =               2; // 0;
	public int         center_index =             -1; // when <0 find center index after pair-wise matching and set (in pattern_mode)
	
	
	
	// "testing and debug" section
	public boolean     combine_empty_only =       true; // false;
	public boolean     late_normalize_iterate =   true;
	public int         test_corr_rad_max =        3;
	public boolean     show_result_images =       true;
	
	// for recalculateFlowXY()
	
	public int         debug_level_optical =  1;
	public int         debug_level_iterate = -1;
	public boolean     enable_debug_images =  true;
	
	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addNumericField("Scale strength for the tiles with no LMA disparity",              this.scale_no_lma_disparity,3,6,"",
				"Reduce strength (but kip less reliable CM-measured disparity for the tiles with no disparity measured with LMA");
		gd.addMessage("Intraframe ERS to pose");
		gd.addNumericField("Previous (in time) frame weight (0.75)",                          this.k_prev,             3,6,"",
				"Earlier frame ERS was sampled closer to interframe");
		gd.addNumericField("Correction coefficient (0.75)",                                   this.ers_to_pose_scale,  3,6,"",
				"Still to fing out the reason why ERS estimates higher value");

		gd.addMessage("Reference and scene macrotiles calculation");
		gd.addNumericField("Absolute disparity tolerance for reference scene",                this.tolerance_absolute_ref,  3,6,"pix",
				"Filter reference scene tiles in a macrotile by disparity (absolute term)");
		gd.addNumericField("Relative disparity tolerance for reference scene",                this.tolerance_relative_ref,  3,6,"pix/pix",
				"Filter reference scene tiles in a macrotile by disparity (relative to disparity)");
		gd.addNumericField("Minimal fraction of disparity-compliant tiles in the reference center", this.center_occupancy_ref,  4,6,"",
				"Discard reference macrotile with less fraction of remaining tiles after filtering by the same disparity in the center 8x8 of the 16x16 macrotile");
		gd.addNumericField("Number of replace-NaN Laplassian passes",                         this.num_laplassian,  0,4,"",
				"Maximal number of repetitions replacing the removed by disparity filtered tiles by weighted average of the 8 neighbors");
		gd.addNumericField("Maximal change threshold for Laplassian NaN replacement",         this.change_laplassian,  5,8,"",
				"Exit replacement passes when maximal change falls below this threshold");
		gd.addNumericField("Absolute disparity interscene tolerance",                         this.tolerance_absolute_inter,  3,6,"pix",
				"Filter scene tiles in a macrotile by interscene disparity comparison (absolute term)");
		gd.addNumericField("Relative disparity interscene tolerance",                         this.tolerance_relative_inter,  3,6,"pix/pix",
				"Filter scene tiles in a macrotile by interscene disparity comparison (relative to disparity term)");
		gd.addNumericField("Occupancy of the scene macrotile",                                this.occupancy_inter,  3,6,"",
				"Fraction of the remaining tiles in a scene macrotile after interframe disparity filtering");

		
		gd.addMessage("Interscene 2D correlation");
		gd.addNumericField("Remove outliers (relative to sigma)",                             this.nsigma,  3,6,"",
				"Remove optical flow macrotiles that differ from weighted average by more that this factor scaled standard deviation");
		gd.addNumericField("Remove outliers (relative to sigma), second pass (-1 disable)",   this.nsigma2,  3,6,"",
				"Second pass of outlier removal (with new sigma). Set to -1 to disable second pass");
		gd.addNumericField("Correlation weight of the intrascene strength channel",           this.chn_weights[0],  3,6,"",
				"Weight of the intrascene correlation strength in interscene macrotile correlation (will be normalized)");
		gd.addNumericField("Correlation weight of the intrascene red channel",                this.chn_weights[1],  3,6,"",
				"Weight of the intrascene average red color in interscene macrotile correlation (will be normalized)");
		gd.addNumericField("Correlation weight of the intrascene blue channel",               this.chn_weights[2],  3,6,"",
				"Weight of the intrascene average blue color in interscene macrotile correlation (will be normalized)");
		gd.addNumericField("Correlation weight of the intrascene green channel",              this.chn_weights[3],  3,6,"",
				"Weight of the intrascene average green color in interscene macrotile correlation (will be normalized)");
		gd.addNumericField("Interscene correlation LPF sigma",                                this.corr_sigma,  3,6,"tiles",
				"LPF 2D interscene correlations before copnverting to pixel domain");
		gd.addNumericField("Fat zero for transform-domain correlation normalization",         this.fat_zero,  3,6,"",
				"Add fat zero to regularize phase correlation");
		gd.addNumericField("Add fraction pixel radius for correlation consolidation",         this.frac_radius,  3,6,"",
				"Just for correlation testing with variable number of consolidated tiles - add this to integer radius");
		gd.addNumericField("Absolute disparity tolerance for macrotile consolidation",        this.tolerance_absolute_macro,  3,6,"pix",
				"Consolidate only those macrotiles correlations that have close average disparitgy");
		gd.addNumericField("Relartive disparity tolerance for macrotile consolidation",       this.tolerance_relative_macro,  3,6,"",
				"Consolidate only those macrotiles correlations that have close average disparitgy");
		
		gd.addMessage("Argmax calculation with center-of-mass");
		gd.addNumericField("Half-size of the square to use for CM argmax",                    this.iradius_cm,    0,4,"tiles",
				"Use this half-size square to calculate center of mass for argmax");
		gd.addNumericField("CM correlation radius (weight = (1/(r/dradius_cm)^2 + 1)",        this.dradius_cm,    3,6,"tiles",
				"Fade-out function for CM argmax() calculation");
		gd.addNumericField("Number of CM iterations",                                         this.refine_num_cm, 0,4,"",
				"Re-calculate weights of the tiles for CM calculation after refining argmax");
		
		gd.addMessage("Optical flow vector refines");
		gd.addNumericField("Magic scale (0.85) after CM argmax",                              this.magic_scale,   3,6,"",
				"Optical vector refinement should be divided by this value before application");
		gd.addNumericField("Maximal number of Optical Flow refine cycles",                    this.max_refines,   0,4,"",
				"Hard limit on th number of 2D correlations/ optical flow vectors refines");
		gd.addNumericField("Number of unconditional refines",                                 this.num_refine_all,0,4,"",
				"Number of refines even if the last vector increment was larger than the previous one");
		gd.addNumericField("Minimal Optical Flow vector change to exit macrotile refinement", this.min_change,    3,6,"pix",
				"Stop re-correlating tile when change is smaller than this threshold (in image pixels)");
		
		gd.addMessage("Confidence calculation");
		gd.addNumericField("Number of the best (closest) Optical Flow neighbors",             this.best_neibs_num,0,4,"",
				"Use this number of the neighbors (of total 8) that are most similar to the current tile to calculate confidence. Zero confidence if there are not enough neighbors.");
		gd.addNumericField("Expected standard deviation of the Optical Flow",                 this.ref_stdev,     3,6,"pix",
				"Calculate for the best neighbors around the current tile: confidence= (ref_stdev ^ 2)/(ref_stdev ^2 + stdev^2)");

		gd.addMessage("Linear and Rotational Velocities");
	    gd.addCheckbox    ("Ignore velocities from individual ERS",                           this.ignore_ers,
			"Ignore relative pose from individa-scene ERS when calculating relative poses");
		gd.addNumericField("Velocities LPF during pairwise fitting",                          this.lpf_pairs,     3,6,"samples",
				"Half of the full width LPF to smooth ERS velocities during consecutive pairs pose fitting");
		gd.addNumericField("Velocities LPF during all-to-reference fitting",                  this.lpf_series,     3,6,"samples",
				"Half of the full width LPF to smooth ERS velocities during all scenes to reference pose fitting");
		
		gd.addMessage("Processing factory calibration (pattern) images, some settings may be useful for field images too");
	    gd.addCheckbox    ("Pattern mode",                                                    this.pattern_mode,
			     "When processing factory calibration images. Not used yet. In pattern mode XYZ adjustemnt (Inter-LMA tab) should be disabled");
		gd.addNumericField("RMS larger than any possible RMS in correct fitting",             this.max_rms_maybe,     3,6,"pix",
				"Now not important, can be set to sufficiently large value (now 7.0)");
		gd.addNumericField("Definitely good RMS for scene pair-wise fitting",                 this.max_rms_sure,     3,6,"pix",
				"When set too high leads to false positive fitting, too low - causes too much testing around.");
		gd.addNumericField("Fitting scan search step",                                        this.pix_step,0,4,"pix",
				"Search around in a spiral with specified pixel step. With the current LWIR resolution/focal length should be 20pix(safe)-50pix(probably OK).");
		gd.addNumericField("Search distance",                                                 this.search_rad,0,4,"steps",
				"Search in a spiral, 0 - just single center, 1 - 3x3 pattern, 2 - 5x5 pattern.");
		gd.addNumericField("Center (reference) index",                                                 this.center_index,0,4,"steps",
				"Will be used as a reference frame for matching all to it. If < 0 it will be calculated and set");

		gd.addMessage("Testing and Debug");
	    gd.addCheckbox    ("Consolidate correlation macrotiles only if the current macrotile is null",  this.combine_empty_only,
			"When false - consolidate all macrotiles, including defined ones");
	    gd.addCheckbox    ("Late normalize (after combining channels in transform domain)",   this.late_normalize_iterate,
			"Assumed true when consolidating macrotiles");
		gd.addNumericField("Maximal consolidation radius to try",                             this.test_corr_rad_max,  0,4,"tiles",
				"Test and dispaly consolidation from from zero to this radius around each macrotile");
	    gd.addCheckbox    ("Show ERS debug images",   this.show_result_images,
			"To be continued - debugging EDRS for LWIR, but ERS was too small, need to debug with faster rotations.");
		
		gd.addNumericField("Debug level for Optical Flow testing",                            this.debug_level_optical,  0,4,"",
				"Apply to series or Optical FLow tests");
		
		gd.addNumericField("Debug level correlation iterations",                              this.debug_level_iterate,  0,4,"",
				"Apply during Optical Flow refinement itgerations");
	    gd.addCheckbox    ("Enable debug images",  this.enable_debug_images,
			"When false - no debug images will be generated regardless of debug level");
	}	
	
	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.scale_no_lma_disparity =       gd.getNextNumber();
		this.k_prev =                       gd.getNextNumber();
		this.ers_to_pose_scale =            gd.getNextNumber();
		this.tolerance_absolute_ref =       gd.getNextNumber();
		this.tolerance_relative_ref =       gd.getNextNumber();
		this.center_occupancy_ref =         gd.getNextNumber();
		this.num_laplassian =         (int) gd.getNextNumber();
		this.change_laplassian =            gd.getNextNumber();
		this.tolerance_absolute_inter =     gd.getNextNumber();
		this.tolerance_relative_inter =     gd.getNextNumber();
		this.occupancy_inter =              gd.getNextNumber();
		this.nsigma =                       gd.getNextNumber();
		this.nsigma2 =                      gd.getNextNumber();
		this.chn_weights[0] =               gd.getNextNumber(); 
		this.chn_weights[1] =               gd.getNextNumber(); 
		this.chn_weights[2] =               gd.getNextNumber(); 
		this.chn_weights[3] =               gd.getNextNumber(); 

		this.corr_sigma =                   gd.getNextNumber();
		this.fat_zero =                     gd.getNextNumber();
		this.frac_radius =                  gd.getNextNumber();
		this.tolerance_absolute_macro =     gd.getNextNumber();
		this.tolerance_relative_macro =     gd.getNextNumber();
		
		this.iradius_cm =             (int) gd.getNextNumber();
		this.dradius_cm =                   gd.getNextNumber();
		this.refine_num_cm =          (int) gd.getNextNumber();
		this.magic_scale =                  gd.getNextNumber();
		this.max_refines =            (int) gd.getNextNumber();
		this.num_refine_all =         (int) gd.getNextNumber();
		this.min_change =                   gd.getNextNumber();

		this.best_neibs_num =         (int) gd.getNextNumber();
		this.ref_stdev =                    gd.getNextNumber();

		this.ignore_ers =                   gd.getNextBoolean();
		this.lpf_pairs =                    gd.getNextNumber();
		this.lpf_series =                   gd.getNextNumber();
		
		this.pattern_mode =                 gd.getNextBoolean();
		this.max_rms_maybe =                gd.getNextNumber();
		this.max_rms_sure =                 gd.getNextNumber();
		this.pix_step =               (int) gd.getNextNumber();
		this.search_rad =             (int) gd.getNextNumber();
		this.center_index =           (int) gd.getNextNumber();
		
		this.combine_empty_only =           gd.getNextBoolean();
		this.late_normalize_iterate =       gd.getNextBoolean();
		this.test_corr_rad_max =      (int) gd.getNextNumber();
		this.show_result_images =           gd.getNextBoolean();
		
		this.debug_level_optical =    (int) gd.getNextNumber();
		this.debug_level_iterate =    (int) gd.getNextNumber();
		this.enable_debug_images =          gd.getNextBoolean();
		
	}	
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"scale_no_lma_disparity",   this.scale_no_lma_disparity+"");
		properties.setProperty(prefix+"k_prev",                   this.k_prev+"");
		properties.setProperty(prefix+"ers_to_pose_scale",        this.ers_to_pose_scale+"");
		properties.setProperty(prefix+"tolerance_absolute_ref",   this.tolerance_absolute_ref+"");
		properties.setProperty(prefix+"tolerance_relative_ref",   this.tolerance_relative_ref+"");
		properties.setProperty(prefix+"center_occupancy_ref",     this.center_occupancy_ref+"");
		properties.setProperty(prefix+"num_laplassian",           this.num_laplassian+"");
		properties.setProperty(prefix+"change_laplassian",        this.change_laplassian+"");
		properties.setProperty(prefix+"tolerance_absolute_inter", this.tolerance_absolute_inter+"");
		properties.setProperty(prefix+"tolerance_relative_inter", this.tolerance_relative_inter+"");
		properties.setProperty(prefix+"occupancy_inter",          this.occupancy_inter+"");
		properties.setProperty(prefix+"nsigma",                   this.nsigma+"");
		properties.setProperty(prefix+"nsigma2",                  this.nsigma2+"");
		for (int i = 0; i < chn_weights.length; i++) {
			properties.setProperty(prefix+"chn_weights_"+i,       this.chn_weights[i]+"");
		}
		properties.setProperty(prefix+"corr_sigma",               this.corr_sigma+"");
		properties.setProperty(prefix+"fat_zero",                 this.fat_zero+"");
		properties.setProperty(prefix+"frac_radius",              this.frac_radius+"");
		properties.setProperty(prefix+"tolerance_absolute_macro", this.tolerance_absolute_macro+"");
		properties.setProperty(prefix+"tolerance_relative_macro", this.tolerance_relative_macro+"");

		properties.setProperty(prefix+"iradius_cm",               this.iradius_cm+"");
		properties.setProperty(prefix+"dradius_cm",               this.dradius_cm+"");
		properties.setProperty(prefix+"refine_num_cm",            this.refine_num_cm+"");
		properties.setProperty(prefix+"magic_scale",              this.magic_scale+"");
		properties.setProperty(prefix+"max_refines",              this.max_refines+"");
		properties.setProperty(prefix+"num_refine_all",           this.num_refine_all+"");
		properties.setProperty(prefix+"min_change",               this.min_change+"");
	
		properties.setProperty(prefix+"best_neibs_num",           this.best_neibs_num+"");
		properties.setProperty(prefix+"ref_stdev",                this.ref_stdev+"");
		
		properties.setProperty(prefix+"ignore_ers",               this.ignore_ers+"");
		properties.setProperty(prefix+"lpf_pairs",                this.lpf_pairs+"");
		properties.setProperty(prefix+"lpf_series",               this.lpf_series+"");
		
		properties.setProperty(prefix+"pattern_mode",             this.pattern_mode+"");
		properties.setProperty(prefix+"max_rms_maybe",            this.max_rms_maybe+"");
		properties.setProperty(prefix+"max_rms_sure",             this.max_rms_sure+"");
		properties.setProperty(prefix+"pix_step",                 this.pix_step+"");
		properties.setProperty(prefix+"search_rad",               this.search_rad+"");
		properties.setProperty(prefix+"center_index",             this.center_index+"");

		properties.setProperty(prefix+"combine_empty_only",       this.combine_empty_only+"");
		properties.setProperty(prefix+"late_normalize_iterate",   this.late_normalize_iterate+"");
		properties.setProperty(prefix+"test_corr_rad_max",        this.test_corr_rad_max+"");
		properties.setProperty(prefix+"show_result_images",       this.show_result_images+"");

		properties.setProperty(prefix+"debug_level_optical",      this.debug_level_optical+"");
		properties.setProperty(prefix+"debug_level_iterate",      this.debug_level_iterate+"");
		properties.setProperty(prefix+"enable_debug_images",      this.enable_debug_images+"");
	}	
	
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"scale_no_lma_disparity")!=null)   this.scale_no_lma_disparity=Double.parseDouble(properties.getProperty(prefix+"scale_no_lma_disparity"));
		if (properties.getProperty(prefix+"k_prev")!=null)                   this.k_prev=Double.parseDouble(properties.getProperty(prefix+"k_prev"));
		if (properties.getProperty(prefix+"ers_to_pose_scale")!=null)        this.ers_to_pose_scale=Double.parseDouble(properties.getProperty(prefix+"ers_to_pose_scale"));
		if (properties.getProperty(prefix+"tolerance_absolute_ref")!=null)   this.tolerance_absolute_ref=Double.parseDouble(properties.getProperty(prefix+"tolerance_absolute_ref"));
		if (properties.getProperty(prefix+"tolerance_relative_ref")!=null)   this.tolerance_relative_ref=Double.parseDouble(properties.getProperty(prefix+"tolerance_relative_ref"));
		if (properties.getProperty(prefix+"center_occupancy_ref")!=null)     this.center_occupancy_ref=Double.parseDouble(properties.getProperty(prefix+"center_occupancy_ref"));
		if (properties.getProperty(prefix+"num_laplassian")!=null)           this.num_laplassian=Integer.parseInt(properties.getProperty(prefix+"num_laplassian"));
		if (properties.getProperty(prefix+"change_laplassian")!=null)        this.change_laplassian=Double.parseDouble(properties.getProperty(prefix+"change_laplassian"));
		if (properties.getProperty(prefix+"tolerance_absolute_inter")!=null) this.tolerance_absolute_inter=Double.parseDouble(properties.getProperty(prefix+"tolerance_absolute_inter"));
		if (properties.getProperty(prefix+"tolerance_relative_inter")!=null) this.tolerance_relative_inter=Double.parseDouble(properties.getProperty(prefix+"tolerance_relative_inter"));
		if (properties.getProperty(prefix+"occupancy_inter")!=null)          this.occupancy_inter=Double.parseDouble(properties.getProperty(prefix+"occupancy_inter"));
		if (properties.getProperty(prefix+"nsigma")!=null)                   this.nsigma=Double.parseDouble(properties.getProperty(prefix+"nsigma"));
		if (properties.getProperty(prefix+"nsigma2")!=null)                  this.nsigma2=Double.parseDouble(properties.getProperty(prefix+"nsigma2"));
		for (int i = 0; i < chn_weights.length; i++) {
			String s_chn_weight = "chn_weights_"+i;
			if (properties.getProperty(prefix+s_chn_weight)!=null)           this.chn_weights[i]=Double.parseDouble(properties.getProperty(prefix+s_chn_weight));
		}
		if (properties.getProperty(prefix+"corr_sigma")!=null)               this.corr_sigma=Double.parseDouble(properties.getProperty(prefix+"corr_sigma"));
		if (properties.getProperty(prefix+"fat_zero")!=null)                 this.fat_zero=Double.parseDouble(properties.getProperty(prefix+"fat_zero"));
		if (properties.getProperty(prefix+"frac_radius")!=null)              this.frac_radius=Double.parseDouble(properties.getProperty(prefix+"frac_radius"));
		if (properties.getProperty(prefix+"tolerance_absolute_macro")!=null) this.tolerance_absolute_macro=Double.parseDouble(properties.getProperty(prefix+"tolerance_absolute_macro"));
		if (properties.getProperty(prefix+"tolerance_relative_macro")!=null) this.tolerance_relative_macro=Double.parseDouble(properties.getProperty(prefix+"tolerance_relative_macro"));

		if (properties.getProperty(prefix+"iradius_cm")!=null)               this.iradius_cm=Integer.parseInt(properties.getProperty(prefix+"iradius_cm"));
		if (properties.getProperty(prefix+"dradius_cm")!=null)               this.dradius_cm=Double.parseDouble(properties.getProperty(prefix+"dradius_cm"));
		if (properties.getProperty(prefix+"refine_num_cm")!=null)            this.refine_num_cm=Integer.parseInt(properties.getProperty(prefix+"refine_num_cm"));
		if (properties.getProperty(prefix+"magic_scale")!=null)              this.magic_scale=Double.parseDouble(properties.getProperty(prefix+"magic_scale"));
		if (properties.getProperty(prefix+"max_refines")!=null)              this.max_refines=Integer.parseInt(properties.getProperty(prefix+"max_refines"));
		if (properties.getProperty(prefix+"num_refine_all")!=null)           this.num_refine_all=Integer.parseInt(properties.getProperty(prefix+"num_refine_all"));
		if (properties.getProperty(prefix+"min_change")!=null)               this.min_change=Double.parseDouble(properties.getProperty(prefix+"min_change"));
		
		if (properties.getProperty(prefix+"best_neibs_num")!=null)           this.best_neibs_num=Integer.parseInt(properties.getProperty(prefix+"best_neibs_num"));
		if (properties.getProperty(prefix+"ref_stdev")!=null)                this.ref_stdev=Double.parseDouble(properties.getProperty(prefix+"ref_stdev"));
		
		if (properties.getProperty(prefix+"ignore_ers")!=null)               this.ignore_ers=Boolean.parseBoolean(properties.getProperty(prefix+"ignore_ers"));
		if (properties.getProperty(prefix+"lpf_pairs")!=null)                this.lpf_pairs=Double.parseDouble(properties.getProperty(prefix+"lpf_pairs"));
		if (properties.getProperty(prefix+"lpf_series")!=null)               this.lpf_series=Double.parseDouble(properties.getProperty(prefix+"lpf_series"));

		if (properties.getProperty(prefix+"pattern_mode")!=null)             this.pattern_mode=Boolean.parseBoolean(properties.getProperty(prefix+"pattern_mode"));
		if (properties.getProperty(prefix+"max_rms_maybe")!=null)            this.max_rms_maybe=Double.parseDouble(properties.getProperty(prefix+"max_rms_maybe"));
		if (properties.getProperty(prefix+"max_rms_sure")!=null)             this.max_rms_sure=Double.parseDouble(properties.getProperty(prefix+"max_rms_sure"));
		if (properties.getProperty(prefix+"pix_step")!=null)                 this.pix_step=Integer.parseInt(properties.getProperty(prefix+"pix_step"));
		if (properties.getProperty(prefix+"search_rad")!=null)               this.search_rad=Integer.parseInt(properties.getProperty(prefix+"search_rad"));
		if (properties.getProperty(prefix+"center_index")!=null)             this.center_index=Integer.parseInt(properties.getProperty(prefix+"center_index"));
		
		if (properties.getProperty(prefix+"combine_empty_only")!=null)       this.combine_empty_only=Boolean.parseBoolean(properties.getProperty(prefix+"combine_empty_only"));
		if (properties.getProperty(prefix+"late_normalize_iterate")!=null)   this.late_normalize_iterate=Boolean.parseBoolean(properties.getProperty(prefix+"late_normalize_iterate"));
		if (properties.getProperty(prefix+"test_corr_rad_max")!=null)        this.test_corr_rad_max=Integer.parseInt(properties.getProperty(prefix+"test_corr_rad_max"));
		if (properties.getProperty(prefix+"show_result_images")!=null)       this.show_result_images=Boolean.parseBoolean(properties.getProperty(prefix+"show_result_images"));
		
		if (properties.getProperty(prefix+"debug_level_optical")!=null)      this.debug_level_optical=Integer.parseInt(properties.getProperty(prefix+"debug_level_optical"));
		if (properties.getProperty(prefix+"debug_level_iterate")!=null)      this.debug_level_iterate=Integer.parseInt(properties.getProperty(prefix+"debug_level_iterate"));
		if (properties.getProperty(prefix+"enable_debug_images")!=null)      this.enable_debug_images=Boolean.parseBoolean(properties.getProperty(prefix+"enable_debug_images"));
		
	}
	@Override
	public OpticalFlowParameters clone() throws CloneNotSupportedException {
		OpticalFlowParameters ofp =     new OpticalFlowParameters();
		ofp.scale_no_lma_disparity =        this.scale_no_lma_disparity;
		ofp.k_prev =                        this.k_prev;
		ofp.ers_to_pose_scale =             this.ers_to_pose_scale;
		ofp.tolerance_absolute_ref =        this.tolerance_absolute_ref;
		ofp.tolerance_relative_ref =        this.tolerance_relative_ref;
		ofp.center_occupancy_ref =          this.center_occupancy_ref;
		ofp.num_laplassian =                this.num_laplassian;
		ofp.change_laplassian =             this.change_laplassian;
		ofp.tolerance_absolute_inter =      this.tolerance_absolute_inter;
		ofp.tolerance_relative_inter =      this.tolerance_relative_inter;
		ofp.occupancy_inter =               this.occupancy_inter;
		ofp.nsigma =                        this.nsigma;
		ofp.nsigma2 =                       this.nsigma2;
		ofp.chn_weights =                   this.chn_weights.clone();
		ofp.corr_sigma =                    this.corr_sigma;
		ofp.fat_zero =                      this.fat_zero;
		ofp.frac_radius =                   this.frac_radius;
		ofp.tolerance_absolute_macro =      this.tolerance_absolute_macro;
		ofp.tolerance_relative_macro =      this.tolerance_relative_macro;

		ofp.iradius_cm =                    this.iradius_cm;
		ofp.dradius_cm =                    this.dradius_cm;
		ofp.refine_num_cm =                 this.refine_num_cm;
		ofp.magic_scale =                   this.magic_scale;
		ofp.max_refines =                   this.max_refines;
		ofp.num_refine_all =                this.num_refine_all;
		ofp.min_change =                    this.min_change;
		
		ofp.best_neibs_num =                this.best_neibs_num;
		ofp.ref_stdev =                     this.ref_stdev;
		
		ofp.ignore_ers =                    this.ignore_ers;
		ofp.lpf_pairs =                     this.lpf_pairs;
		ofp.lpf_series =                    this.lpf_series; 
		
		ofp.pattern_mode =                  this.pattern_mode; 
		ofp.max_rms_maybe =                 this.max_rms_maybe; 
		ofp.max_rms_sure =                  this.max_rms_sure; 
		ofp.pix_step =                      this.pix_step; 
		ofp.search_rad =                    this.search_rad; 
		ofp.center_index =                  this.center_index; 
		
		ofp.combine_empty_only =            this.combine_empty_only;
		ofp.late_normalize_iterate =        this.late_normalize_iterate;
		ofp.test_corr_rad_max =             this.test_corr_rad_max;
		ofp.show_result_images =            this.show_result_images;
		
		ofp.debug_level_optical =           this.debug_level_optical;
		ofp.debug_level_iterate =           this.debug_level_iterate;
		ofp.enable_debug_images =           this.enable_debug_images;
		return ofp;
	}	
}
