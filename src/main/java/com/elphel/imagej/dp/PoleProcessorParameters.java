package com.elphel.imagej.dp;
import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;

/**
 **
 ** PoleProcessorParameters - parameters for PoleProcessor class that identifies
 ** thin vertical object (such as street lights) and measures distances
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  PoleProcessorParameters.java is free software: you can redistribute it and/or modify
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

public class PoleProcessorParameters {
	public boolean poles_en =      true;
	// Creating initial pole seeds (tops need to be detected and available in the available DSI
	public double  seed_rsrength =      0.8;   // pole top seed relative strength (fraction of the strong trusted value that may change with correlation fat zero or other parameters
	public int     seed_down =          5;     // seed disparity should be greater than average of this number of tiles below
	public double  seed_aover =         0.2;   // how much seed disparity should exceed that of the tiles below
	public double  seed_rover =         0.05;  // increase the previous by each 1 pixel of disparity
	public double  max_disparity =      5.0;   // maximal disparity to consider (near objects can be handled well without these tricks
	// Seed (pole heads) tiles initial clustering
	public double  max_dd =             0.1;   // Maximal difference among seed tiles (compares to average of the already included)
	public double  max_dx =             5;     // Maximal horizontal difference between tiles to be included in the same initial cluster
	public double  max_dy =            -1.0;   // Maximal horizontal difference between tiles to be included in the same initial cluster
	public int     ext_side =           2;     // Extend initial cluster bounding box horizontally
	public int     ext_up =             2;     // Extend initial cluster bounding box up
	public int     ext_down =          30;     // Extend initial cluster bounding box down
	public double  pedestal_strength =  0.15;  // Pole pedestal minimal normalized strength
	public double  pedestal_disp_over = 0.3;   // Pole pedestal disparity over the head
	// cluster merge and split parameters
	public boolean merge_extended =     false; // merge if only original clusters intersect
	public double  merge_rtolerance =   1.5;   // max_dd;
	public int     split_min_dist =     3;     // minimal distance between local maximums to split
	public boolean split_must_zero =    false; // true; // false; // there should be complete zero-strength on the profile to enable split
    // Creating pole centerline and selection mask
	public double  max_diff =           2.0;   // maximal horizontal distance between the pole centerline and contributing tiles
	public double  max_tilt =           0.1;   // maximal pole center line tilt (dx/dy) from the vertical
	public double  max_rms =            1.0;   // Maximal horizontal weighted RMS between the contributing tiles and the centerline
	public int     min_tiles =          3;     // Minimal number of the remaining tiles to calculate centerline
	public double  damp_tilt =          0.1;   // "Force" to make centerline vertical
	public boolean use_seed =           true;  // use the original seed (pole head) in the selection mask
	public double  hwidth =             0.75;  // maximal horizontal distance between the centerline and pole tiles in the selection
	public double  disp_aover =         0.2;   // interrupt mask if there is some FG object in front of it, with disparity exceeding by this margin
	public double  disp_rover =         0.05;  // increment above value per each pixel of disparity
	public int     min_neibs =          3;     // minimal number of neighbors to keep tile in the pole mask
	// More tile filtering: trim unoccupied bottoms and check that the cluster head is separated on both right and left sides
	public boolean trim_bottoms  =      true;  // Trim bounding box bottom that does not contain any measured tiles
	public int     extra_margin =       2;     // check for the pole head separation from the same/closer objects on both (right and left) sides, up to this far from the extended bounding box
	public double  sep_min_strength =   0.07;  // require that the pole head is separated from right and left with the weaker tiles or farther ones
	public double  sep_disp_adiff =     0.15;  // separation should have objects with disparity by this lower than the head
	public double  sep_disp_rdiff =     0.05;  // increase required disparity separation margin proportionally to the absolute disparity
	// measurements/re-measurements parameters with notch filtering to emphasize vertical lines
	public int     remeasure_hheight =  4;     // half maximum number of tiles vertical column to combine correlations
	public double  high_disp_tolerance =0.3;   // count only tiles that had measured residual disparity not exceeding this value
	public double  low_disp_tolerance = 0.3;   // count only tiles that had measured residual disparity not less than minus this value
	public double  min_strength =       0.15;  // minimal measured strength to consider
	public boolean use_fittest =        true;  // when selecting measurement among ones with different vertical integration, use the one closest to 0, false - use strongest result
	// Updating disparity using the measurement results
	public int     max_refines =        20;    // maximal number of re-measure/update disparity cycles
	public double  min_refine_diff =    0.001; // minimal disparity difference to keep refining
	public double  disparity_scale =    0.3; // 2;  // target disparity to differential disparity scale (baseline ratio)
	public double  diff_power =         1.0;  // bias towards higher disparities - (disparity+offset) is raised to this power and applied to weight
	//                                      if 0.0 - do not apply value to weight
	public double  diff_offset =        0.5;  // add to measured differential disparity before raising to specified power
	public int     cut_bottom =         2;    // cut few tile rows from the very bottom - they may be influenced by ground objects
	public double  keep_bottom =        0.1;  // do not cut more that this fraction of the bounding box height
	// Filtering clusters (removing those that fail
	public boolean filter_headlessOK =     false;  // allow clusters that do not have heads (generated by splitting)
	public boolean filter_en =             true;   // enable cluster filtering (off - pass all)
	public int     filter_min_height =    10;      // minimal extended bounding box height (good - min 12)
	public double  filter_min_fraction =   0.5;    // minimal fraction of the pole mask that has measured valid tiles (worst 0.3478 , next - 0.57)
	public double  filter_min_frac_height =0.5;    // vertical profile: minimal fraction of the pole selection mask that has valid measurements worst 0.25 , next - 0.556
	public double  filter_min_frac_box =   0.7;    // same as above, but denominator includes full mask height including gaps if any worst 0.75 , next - 0.556
	public double  filter_min_disparity =  0.3;    // minimal pole absolutre disparity
	public double  filter_max_disparity =  3.05;   // maximal mole absolute disparity (2.95 - 3.15)

	public int     poles_debug_level =    -1;      // debug level for the poles



	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addCheckbox    ("Enable poles detection and usage",                                                                   this.poles_en,
				"Run automatically when creating ground truth data. If disable you may run it with \"Poles GT\" command button");
		gd.addMessage("Creating initial pole seeds (tops need to be detected and available in the available DSI");
		gd.addNumericField("Pole relative (to trusted tiles strength) to qualify for the pole seed",                             this.seed_rsrength,  3,6,"",
				"Pole top seed relative strength (fraction of the strong trusted value that may change with correlation fat zero or other parameters");
		gd.addNumericField("Seed disparity should be greater than average of this number of tiles below it",                     this.seed_down,  0,3,"tiles",
				"Seed (pole head) tiles should appear closer (higher disparity) than the (presumably BG) tiles resolved below it. This is how many tiles below the seed should be averaged to compare");
		gd.addNumericField("Seed disparity should exceed BG below by this margin",                                               this.seed_aover,  3,6,"pix",
				"How much seed disparity should exceed that of the tiles below");
		gd.addNumericField("Relative seed disparity margin (per each disparity pixel)",                                          this.seed_rover,  3,6,"pix/pix",
				"Add to the previous parameter after multiplying by the absolute disparity");
		gd.addNumericField("Maximal disparity to consider as a seed",                                                            this.max_disparity,  3,6,"pix",
				"Near objects are handled nicely without such tricks, and they usually have larger horizontal dimension (more than a tile wide)");

		gd.addMessage("Seed (pole heads) tiles initial clustering");
		gd.addNumericField("Maximal disparity difference among seed tiles",                                                      this.max_dd,  3,6,"pix",
				"Disparity of each new tile is compared to average disparity of the already included tiles. Set to -1 if don't care");
		gd.addNumericField("Maximal horizontal difference between tiles to be included in the same initial cluster",             this.max_dx,  3,6,"tiles",
				"Compare added tile horizontal position to that of the average of the already included tiles. Set to -1.0 if don't care");
		gd.addNumericField("Maximal vertical difference between tiles to be included in the same initial cluster",               this.max_dy,  3,6,"tiles",
				"Compare added tile vertical position to that of the average of the already included tiles. Set to -1.0 if don't care");
		gd.addNumericField("Extend initial cluster horizontally (each side)",                                                    this.ext_side,  0,3,"tiles",
				"Extend initial cluster bounding box (and search area) horizontally, each side");
		gd.addNumericField("Extend initial cluster up from the initial seed cluster",                                            this.ext_up,  0,3,"tiles",
				"Extend initial cluster bounding box (and search area) up");
		gd.addNumericField("Extend initial cluster bounding box down",                                                           this.ext_down,  0,3,"tiles",
				"Search for the pole pedestal (tiles definitely closer than the seed tiles) this many tiles down from the seed");
		gd.addNumericField("Pole pedestal minimal normalized strength",                                                          this.pedestal_strength,  3,6,"",
				"Minimal normalized strength (after subtracting strength floor and optionally raising to specified power of teh pole base");
		gd.addNumericField("Pole pedestal disparity over the head",                                                              this.pedestal_disp_over,  3,6,"pix",
				"Pole base should have tiles closer (and so have higher disparity) than the pole seed by this margin");

		gd.addMessage("Cluster merge and split parameters");
		gd.addCheckbox    ("Merge if even extended bounding boxes intersect",                                                    this.merge_extended,
				"Unchecked - merge only if the seed bounding boxes intersect");
		gd.addNumericField("Merge relative (to maximum cluster disparity difference) disparity tolerance",                       this.merge_rtolerance,  3,6,"",
				"Merge if the disparity differs by less than this times scaled seed disparity allowed difference");
		gd.addNumericField("Minimal distance between local maximums to split",                                                   this.split_min_dist,  0,3,"tiles",
				"Calculate horizontal profile for the cluster seed and split if the distance between local maximums does not exceed this value");
		gd.addCheckbox    ("There should be a complete zero (empty) point on the horizontal profile to allow split",             this.split_must_zero,
				"Local maximums for the cluster split should be separated by at least one complete empty point");

		gd.addMessage("Creating pole centerline and selection mask");
 		gd.addNumericField("Distance from the centerline",                                                                        this.max_diff,  3,6,"tiles",
				"Maximal horizontal distance between the pole centerline and contributing tiles");
 		gd.addNumericField("Maximal pole tilt",                                                                                   this.max_tilt,  3,6,"",
				"Maximal pole center line tilt (dx/dy) from the vertical");
 		gd.addNumericField("Maximal horizontal RMS from the centerline",                                                          this.max_rms,  3,6,"tiles",
				"Maximal horizontal weighted RMS between the contributing tiles and the centerline");
		gd.addNumericField("Minimal number of tiles for the centerline calculation",                                              this.min_tiles,  0,3,"tiles",
				"Minimal number of the remaining tiles to calculate centerline");
 		gd.addNumericField("\"Force\" to make centerline vertical",                                                               this.damp_tilt,  3,6,"",
				"Relative penalty for tilted pole centerline");
		gd.addCheckbox    ("Include cluster seed tiles in the cluster tile selection",                                            this.use_seed,
				"Include cluster seed tiles in addition to the centerline in the cluster tile selection");
 		gd.addNumericField("Half-width of the pole selection",                                                                    this.hwidth,  3,6,"tiles",
				"Maximal horizontal distance between the centerline and pole tiles in the selection");
 		gd.addNumericField("Disapritry difference to occlude the pole by the FG objects",                                         this.disp_aover,  3,6,"pix",
				"Interrupt pole mask if there is some FG object in front of it, with disparity exceeding by this margin");
 		gd.addNumericField("Relative disparity difference required for occlusion",                                                this.disp_rover,  3,6,"pix/pix",
				"Increase the required disparity difference proportionally to the absolute disparity");
		gd.addNumericField("Minimal number of neighbors to keep tile in the pole mask",                                           this.min_neibs,  0,3,"tiles",
				"Filtering out lone or poorly connected tiles in the pole selection mask");

		gd.addMessage("More tile filtering: trim unoccupied bottoms and check that the cluster head is separated on both right and left sides");
		gd.addCheckbox    ("Trim bounding box bottom that does not contain any measured tiles",                                   this.trim_bottoms,
				"Remove unused bottom portion of the extended bounding box of each cluster");
		gd.addNumericField("Extra horizontal margins to look for cluster separation",                                             this.extra_margin,  0,3,"tiles",
				"Cluster head should be an \"island\" not the overhanging part of some larger objet. This value specifies how far to look for connections to other objects");
 		gd.addNumericField("Minimal normalized strength of the right/left head connection to be considered",                      this.sep_min_strength,  3,6,"",
				"When looking for the suspected pole head connection to right/left objects consider only tiles stronger thatn that");
 		gd.addNumericField("Disconnection disparity difference",                                                                  this.sep_disp_adiff,  3,6,"pix",
				"To be disconnected, the pole head has to have neighbors with disparity this less than its own disparity");
 		gd.addNumericField("Relative disconnection disparity difference addition",                                                this.sep_disp_rdiff,  3,6,"pix/pix",
				"Increase the previous value proportionally to the absolute disparity");

		gd.addMessage("Measurements/re-measurements parameters with notch filtering to emphasize vertical lines and suppress random texture of the background");
		gd.addNumericField("Half maximum number of tiles vertical column to combine correlations",                                this.remeasure_hheight,  0,3,"tiles",
				"Try to acquire poles using different number of vertically combined tiles, staring from 0 (center only). 1 - up to 3 tiles, 2 - up to 5 tiles, and so on");
 		gd.addNumericField("Maximal measured (residual) disparity",                                                               this.high_disp_tolerance,  3,6,"pix/pix",
				"Count only tiles that had measured residual disparity not exceeding this value");
 		gd.addNumericField("Absolute value of the minimal measured (residual) disparity",                                         this.low_disp_tolerance,  3,6,"pix/pix",
				"Count only tiles that had measured residual disparity not less than negative of this value");
 		gd.addNumericField("Minimal measured strength to consider",                                                               this.min_strength,  3,6,"pix/pix",
				"Increase the previous value proportionally to the absolute disparity");
		gd.addCheckbox    ("Use fittest of the measurements for each tile (unchecked - use strongest)",                            this.use_fittest,
				"when selecting measurement among ones with different vertical integration, use the one closest to 0, false - use strongest result");

		gd.addMessage("Updating disparity using the measurement results");
		gd.addNumericField("Maximal number of re-measure/update disparity cycles",                                                this.max_refines,  0,3,"",
				"Number of update pole disparity - re-measure cycles to perform (if not reached maximal diff limit below");
 		gd.addNumericField("Minimal disparity difference to keep refining",                                                       this.min_refine_diff,  3,6,"",
				"Stop updating disparity for a pole if the last difference is less than this value");
 		gd.addNumericField("Target disparity to differential disparity scale (baseline ratio)",                                   this.disparity_scale,  3,6,"",
 				"Ideally it should be the ratio of the master camera baseline to inter-camera baseline, but because the interaction"+
 				"of the correlation window and correlation maximum width reduces reported disparity, it should be increased");
 		gd.addNumericField("Increase weight of the higher disparity measurements (>=1.0), 0.0 - use normal weighted average ",    this.diff_power,  3,6,"",
				"Bias towards higher disparities - (disparity+offset) is raised to this power and applied to weight if 0.0 - do not apply value to weight");
 		gd.addNumericField("Add to disparity before raising to the power specified above",                                        this.diff_offset,  3,6,"",
				"Disregard samples below negative of this value, multiply weight by the sum of the residual disparity (around 0.0) and this value raised to certain power");
		gd.addNumericField("Do not use this many bottom tiles as their disparity may be influenced by the pedestal",              this.cut_bottom,  0,3,"tiles",
				"These cut-off tiles are still used in the pole mask. The amount of cut tiles may be reduced for the short objects by the next parameter");
 		gd.addNumericField("Do not cut pole bottoms by more than this fraction of its height",                                    this.keep_bottom,  3,6,"",
				"This parameter allows to reduce the number specified in the previous one for the short poles");

 		gd.addMessage("Filtering clusters (removing those that fail)");
		gd.addCheckbox    ("Allow headless clusters ",                                                                            this.filter_headlessOK,
				"Allow clusters that do not have heads (generated by splitting)");
		gd.addCheckbox    ("Enable cluster filtering",                                                                            this.filter_en,
				"Uncheck to process and report all clusters (headless are controlled independently)");
		gd.addNumericField("Minimal cluster height",                                                                              this.filter_min_height,  0,3,"tiles",
				"Minimal cluster extended bounding box height");
 		gd.addNumericField("minimal tile ratio measured/pole mask",                                                               this.filter_min_fraction,  3,6,"",
				"Minimal fraction of the pole mask that has measured valid tiles");
 		gd.addNumericField("Minimal ratio of measured/mask in the vertical profile",                                              this.filter_min_frac_height,  3,6,"",
				"Vertical profile: minimal fraction of the pole selection mask that has valid measurements");
 		gd.addNumericField("Minimal ration of the measured/full height tiles in the vertical profile",                            this.filter_min_frac_box,  3,6,"",
				"Denominator includes full selection height, including gaps if any");
 		gd.addNumericField("Minimal final pole disparity",                                                                        this.filter_min_disparity,  3,6,"pix",
				"Exclude too far poles as determined after disparity refinement");
 		gd.addNumericField("Minimal final pole disparity",                                                                        this.filter_max_disparity,  3,6,"pix",
 				"Exclude too far poles as determined after disparity refinement. Closer poles are usually resolved by normal correlation, pole extraction"+
 				"make the result objects fronto parallel billboards");
		gd.addNumericField("Debug level for the pole processing",                                                                 this.poles_debug_level,  0,3,"",
				"Change independently of the global debug level, that later may still be used to disable ");
/*
	// Filtering clusters (removing those that fail

 		gd.addNumericField("",                       this.merge_rtolerance,  3,6,"",
				"");
		gd.addNumericField("",                       this.split_min_dist,  0,3,"tiles",
				"");
		gd.addCheckbox    ("",                                                    this.merge_extended,
				"");
	// Filtering clusters (removing those that fail
		this.filter_headlessOK,
		this.filter_en,
		this.filter_min_height,  0,3,"tiles",
 		this.filter_min_fraction,  3,6,"",
 		this.filter_min_frac_height,  3,6,"",
 		this.filter_min_frac_box,  3,6,"",
 		this.filter_min_disparity,  3,6,"pix",
 		this.filter_max_disparity,  3,6,"pix",
		this.poles_debug_level,  0,3,"",
 */
	}


	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.poles_en =                     gd.getNextBoolean();
		// Creating initial pole seeds (tops need to be detected and available in the available DSI
		this.seed_rsrength =                gd.getNextNumber();
		this.seed_down =              (int) gd.getNextNumber();
		this.seed_aover =                   gd.getNextNumber();
		this.seed_rover =                   gd.getNextNumber();
		this.max_disparity =                gd.getNextNumber();
		// Seed (pole heads) tiles initial clustering
		this.max_dd =                       gd.getNextNumber();
		this.max_dx =                       gd.getNextNumber();
		this.max_dy =                       gd.getNextNumber();
		this.ext_side =               (int) gd.getNextNumber();
		this.ext_up =                 (int) gd.getNextNumber();
		this.ext_down =               (int) gd.getNextNumber();
		this.pedestal_strength =            gd.getNextNumber();
		this.pedestal_disp_over =           gd.getNextNumber();
		// Cluster merge and split parameters
		this.merge_extended =               gd.getNextBoolean();
		this.merge_rtolerance =             gd.getNextNumber();
		this.split_min_dist =         (int) gd.getNextNumber();
		this.split_must_zero =              gd.getNextBoolean();
	    // Creating pole centerline and selection mask
 		this.max_diff =                     gd.getNextNumber();
 		this.max_tilt =                     gd.getNextNumber();
 		this.max_rms =                      gd.getNextNumber();
		this.min_tiles =              (int) gd.getNextNumber();
 		this.damp_tilt =                    gd.getNextNumber();
		this.use_seed=                      gd.getNextBoolean();
 		this.hwidth =                       gd.getNextNumber();
 		this.disp_aover =                   gd.getNextNumber();
 		this.disp_rover =                   gd.getNextNumber();
		this.min_neibs =              (int) gd.getNextNumber();
		// More tile filtering: trim unoccupied bottoms and check that the cluster head is separated on both right and left sides
		this.trim_bottoms =                 gd.getNextBoolean();
		this.extra_margin =           (int) gd.getNextNumber();
 		this.sep_min_strength =             gd.getNextNumber();
 		this.sep_disp_adiff =               gd.getNextNumber();
 		this.sep_disp_rdiff =               gd.getNextNumber();
 		//Measurements/re-measurements parameters with notch filtering to emphasize vertical lines and suppress random texture of the background
		this.remeasure_hheight =      (int) gd.getNextNumber();
 		this.high_disp_tolerance =          gd.getNextNumber();
 		this.low_disp_tolerance =           gd.getNextNumber();
 		this.min_strength =                 gd.getNextNumber();
		this.use_fittest =                  gd.getNextBoolean();
	 	// Updating disparity using the measurement results
 	 	this.max_refines =            (int) gd.getNextNumber();
 		this.min_refine_diff =              gd.getNextNumber();
 		this.disparity_scale =              gd.getNextNumber();
 		this.diff_power =                   gd.getNextNumber();
 		this.diff_offset =                  gd.getNextNumber();
		this.cut_bottom =             (int) gd.getNextNumber();
 		this.keep_bottom =                  gd.getNextNumber();
 		// Filtering clusters (removing those that fail
		this.filter_headlessOK =            gd.getNextBoolean();
		this.filter_en =                    gd.getNextBoolean();
		this.filter_min_height =      (int) gd.getNextNumber();
 		this.filter_min_fraction =          gd.getNextNumber();
 		this.filter_min_frac_height =       gd.getNextNumber();
 		this.filter_min_frac_box =          gd.getNextNumber();
 		this.filter_min_disparity =         gd.getNextNumber();
 		this.filter_max_disparity =         gd.getNextNumber();
		this.poles_debug_level =      (int) gd.getNextNumber();
	}

	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"poles_en",                      this.poles_en+"");
		// Creating initial pole seeds (tops need to be detected and available in the available DSI
		properties.setProperty(prefix+"seed_rsrength",                 this.seed_rsrength+"");
		properties.setProperty(prefix+"seed_down",                     this.seed_down+"");
		properties.setProperty(prefix+"seed_aover",                    this.seed_aover+"");
		properties.setProperty(prefix+"seed_rover",                    this.seed_rover+"");
		properties.setProperty(prefix+"max_disparity",                 this.max_disparity+"");
		// Seed (pole heads) tiles initial clustering
		properties.setProperty(prefix+"max_dd",                        this.max_dd+"");
		properties.setProperty(prefix+"max_dx",                        this.max_dx+"");
		properties.setProperty(prefix+"max_dy",                        this.max_dy+"");
		properties.setProperty(prefix+"ext_side",                      this.ext_side+"");
		properties.setProperty(prefix+"ext_up",                        this.ext_up+"");
		properties.setProperty(prefix+"ext_down",                      this.ext_down+"");
		properties.setProperty(prefix+"pedestal_strength",             this.pedestal_strength+"");
		properties.setProperty(prefix+"pedestal_disp_over",            this.pedestal_disp_over+"");
		// Cluster merge and split parameters
		properties.setProperty(prefix+"merge_extended",                this.merge_extended+"");
		properties.setProperty(prefix+"merge_rtolerance",              this.merge_rtolerance+"");
		properties.setProperty(prefix+"split_min_dist",                this.split_min_dist+"");
		properties.setProperty(prefix+"split_must_zero",               this.split_must_zero+"");
	    // Creating pole centerline and selection mask
		properties.setProperty(prefix+"max_diff",                      this.max_diff+"");
		properties.setProperty(prefix+"max_tilt",                      this.max_tilt+"");
		properties.setProperty(prefix+"max_rms",                       this.max_rms+"");
		properties.setProperty(prefix+"min_tiles",                     this.min_tiles+"");
		properties.setProperty(prefix+"damp_tilt",                     this.damp_tilt+"");
		properties.setProperty(prefix+"use_seed",                      this.use_seed+"");
		properties.setProperty(prefix+"hwidth",                        this.hwidth+"");
		properties.setProperty(prefix+"disp_aover",                    this.disp_aover+"");
		properties.setProperty(prefix+"disp_rover",                    this.disp_rover+"");
		properties.setProperty(prefix+"min_neibs",                     this.min_neibs+"");
		// More tile filtering: trim unoccupied bottoms and check that the cluster head is separated on both right and left sides
		properties.setProperty(prefix+"trim_bottoms",                  this.trim_bottoms+"");
		properties.setProperty(prefix+"extra_margin",                  this.extra_margin+"");
		properties.setProperty(prefix+"sep_min_strength",              this.sep_min_strength+"");
		properties.setProperty(prefix+"sep_disp_adiff",                this.sep_disp_adiff+"");
		properties.setProperty(prefix+"sep_disp_rdiff",                this.sep_disp_rdiff+"");
 		//Measurements/re-measurements parameters with notch filtering to emphasize vertical lines and suppress random texture of the background
		properties.setProperty(prefix+"remeasure_hheight",             this.remeasure_hheight+"");
		properties.setProperty(prefix+"high_disp_tolerance",           this.high_disp_tolerance+"");
		properties.setProperty(prefix+"low_disp_tolerance",            this.low_disp_tolerance+"");
		properties.setProperty(prefix+"min_strength",                  this.min_strength+"");
		properties.setProperty(prefix+"use_fittest",                   this.use_fittest+"");
	 	// Updating disparity using the measurement results
		properties.setProperty(prefix+"max_refines",                   this.max_refines+"");
		properties.setProperty(prefix+"min_refine_diff",               this.min_refine_diff+"");
		properties.setProperty(prefix+"disparity_scale",               this.disparity_scale+"");
		properties.setProperty(prefix+"diff_power",                    this.diff_power+"");
		properties.setProperty(prefix+"diff_offset",                   this.diff_offset+"");
		properties.setProperty(prefix+"cut_bottom",                    this.cut_bottom+"");
		properties.setProperty(prefix+"keep_bottom",                   this.keep_bottom+"");
		// Filtering clusters (removing those that fail
		properties.setProperty(prefix+"filter_headlessOK",             this.filter_headlessOK+"");
		properties.setProperty(prefix+"filter_en",                     this.filter_en+"");
		properties.setProperty(prefix+"filter_min_height",             this.filter_min_height+"");
		properties.setProperty(prefix+"filter_min_fraction",           this.filter_min_fraction+"");
		properties.setProperty(prefix+"filter_min_frac_height",        this.filter_min_frac_height+"");
		properties.setProperty(prefix+"filter_min_frac_box",           this.filter_min_frac_box+"");
		properties.setProperty(prefix+"filter_min_disparity",          this.filter_min_disparity+"");
		properties.setProperty(prefix+"filter_max_disparity",          this.filter_max_disparity+"");
		properties.setProperty(prefix+"poles_debug_level",             this.poles_debug_level+"");
	}
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"poles_en")!=null)                this.poles_en=Boolean.parseBoolean(properties.getProperty(prefix+"poles_en"));
		// Creating initial pole seeds (tops need to be detected and available in the available DSI
		if (properties.getProperty(prefix+"seed_rsrength")!=null)           this.seed_rsrength=Double.parseDouble(properties.getProperty(prefix+"seed_rsrength"));
		if (properties.getProperty(prefix+"seed_down")!=null)               this.seed_down=Integer.parseInt(properties.getProperty(prefix+"seed_down"));
		if (properties.getProperty(prefix+"seed_aover")!=null)              this.seed_aover=Double.parseDouble(properties.getProperty(prefix+"seed_aover"));
		if (properties.getProperty(prefix+"seed_rover")!=null)              this.seed_rover=Double.parseDouble(properties.getProperty(prefix+"seed_rover"));
		if (properties.getProperty(prefix+"max_disparity")!=null)           this.max_disparity=Double.parseDouble(properties.getProperty(prefix+"max_disparity"));
		// Seed (pole heads) tiles initial clustering
		if (properties.getProperty(prefix+"max_dd")!=null)                  this.max_dd=Double.parseDouble(properties.getProperty(prefix+"max_dd"));
		if (properties.getProperty(prefix+"max_dx")!=null)                  this.max_dx=Double.parseDouble(properties.getProperty(prefix+"max_dx"));
		if (properties.getProperty(prefix+"max_dy")!=null)                  this.max_dy=Double.parseDouble(properties.getProperty(prefix+"max_dy"));
		if (properties.getProperty(prefix+"ext_side")!=null)                this.ext_side=Integer.parseInt(properties.getProperty(prefix+"ext_side"));
		if (properties.getProperty(prefix+"ext_up")!=null)                  this.ext_up=Integer.parseInt(properties.getProperty(prefix+"ext_up"));
		if (properties.getProperty(prefix+"ext_down")!=null)                this.ext_down=Integer.parseInt(properties.getProperty(prefix+"ext_down"));
		if (properties.getProperty(prefix+"pedestal_strength")!=null)       this.pedestal_strength=Double.parseDouble(properties.getProperty(prefix+"pedestal_strength"));
		if (properties.getProperty(prefix+"pedestal_disp_over")!=null)      this.pedestal_disp_over=Double.parseDouble(properties.getProperty(prefix+"pedestal_disp_over"));
		// Cluster merge and split parameters
		if (properties.getProperty(prefix+"merge_extended")!=null)          this.merge_extended=Boolean.parseBoolean(properties.getProperty(prefix+"merge_extended"));
		if (properties.getProperty(prefix+"merge_rtolerance")!=null)        this.merge_rtolerance=Double.parseDouble(properties.getProperty(prefix+"merge_rtolerance"));
		if (properties.getProperty(prefix+"split_min_dist")!=null)          this.split_min_dist=Integer.parseInt(properties.getProperty(prefix+"split_min_dist"));
		if (properties.getProperty(prefix+"split_must_zero")!=null)         this.split_must_zero=Boolean.parseBoolean(properties.getProperty(prefix+"split_must_zero"));
	    // Creating pole centerline and selection mask
		if (properties.getProperty(prefix+"max_diff")!=null)                this.max_diff=Double.parseDouble(properties.getProperty(prefix+"max_diff"));
		if (properties.getProperty(prefix+"max_tilt")!=null)                this.max_tilt=Double.parseDouble(properties.getProperty(prefix+"max_tilt"));
		if (properties.getProperty(prefix+"max_rms")!=null)                 this.max_rms=Double.parseDouble(properties.getProperty(prefix+"max_rms"));
		if (properties.getProperty(prefix+"min_tiles")!=null)               this.min_tiles=Integer.parseInt(properties.getProperty(prefix+"min_tiles"));
		if (properties.getProperty(prefix+"damp_tilt")!=null)               this.damp_tilt=Double.parseDouble(properties.getProperty(prefix+"damp_tilt"));
		if (properties.getProperty(prefix+"use_seed")!=null)                this.use_seed=Boolean.parseBoolean(properties.getProperty(prefix+"use_seed"));
		if (properties.getProperty(prefix+"hwidth")!=null)                  this.hwidth=Double.parseDouble(properties.getProperty(prefix+"hwidth"));
		if (properties.getProperty(prefix+"disp_aover")!=null)              this.disp_aover=Double.parseDouble(properties.getProperty(prefix+"disp_aover"));
		if (properties.getProperty(prefix+"disp_rover")!=null)              this.disp_rover=Double.parseDouble(properties.getProperty(prefix+"disp_rover"));
		if (properties.getProperty(prefix+"min_neibs")!=null)               this.min_neibs=Integer.parseInt(properties.getProperty(prefix+"min_neibs"));
		// More tile filtering: trim unoccupied bottoms and check that the cluster head is separated on both right and left sides
		if (properties.getProperty(prefix+"trim_bottoms")!=null)            this.trim_bottoms=Boolean.parseBoolean(properties.getProperty(prefix+"trim_bottoms"));
		if (properties.getProperty(prefix+"extra_margin")!=null)            this.extra_margin=Integer.parseInt(properties.getProperty(prefix+"extra_margin"));
		if (properties.getProperty(prefix+"sep_min_strength")!=null)        this.sep_min_strength=Double.parseDouble(properties.getProperty(prefix+"sep_min_strength"));
		if (properties.getProperty(prefix+"sep_disp_adiff")!=null)          this.sep_disp_adiff=Double.parseDouble(properties.getProperty(prefix+"sep_disp_adiff"));
		if (properties.getProperty(prefix+"sep_disp_rdiff")!=null)          this.sep_disp_rdiff=Double.parseDouble(properties.getProperty(prefix+"sep_disp_rdiff"));
 		//Measurements/re-measurements parameters with notch filtering to emphasize vertical lines and suppress random texture of the background
		if (properties.getProperty(prefix+"remeasure_hheight")!=null)       this.remeasure_hheight=Integer.parseInt(properties.getProperty(prefix+"remeasure_hheight"));
		if (properties.getProperty(prefix+"high_disp_tolerance")!=null)     this.high_disp_tolerance=Double.parseDouble(properties.getProperty(prefix+"high_disp_tolerance"));
		if (properties.getProperty(prefix+"low_disp_tolerance")!=null)      this.low_disp_tolerance=Double.parseDouble(properties.getProperty(prefix+"low_disp_tolerance"));
		if (properties.getProperty(prefix+"min_strength")!=null)            this.min_strength=Double.parseDouble(properties.getProperty(prefix+"min_strength"));
		if (properties.getProperty(prefix+"use_fittest")!=null)             this.use_fittest=Boolean.parseBoolean(properties.getProperty(prefix+"use_fittest"));
		// Updating disparity using the measurement results
		if (properties.getProperty(prefix+"max_refines")!=null)             this.max_refines=Integer.parseInt(properties.getProperty(prefix+"max_refines"));
		if (properties.getProperty(prefix+"min_refine_diff")!=null)         this.min_refine_diff=Double.parseDouble(properties.getProperty(prefix+"min_refine_diff"));
		if (properties.getProperty(prefix+"disparity_scale")!=null)         this.disparity_scale=Double.parseDouble(properties.getProperty(prefix+"disparity_scale"));
		if (properties.getProperty(prefix+"diff_power")!=null)              this.diff_power=Double.parseDouble(properties.getProperty(prefix+"diff_power"));
		if (properties.getProperty(prefix+"diff_offset")!=null)             this.diff_offset=Double.parseDouble(properties.getProperty(prefix+"diff_offset"));
		if (properties.getProperty(prefix+"cut_bottom")!=null)              this.cut_bottom=Integer.parseInt(properties.getProperty(prefix+"cut_bottom"));
		if (properties.getProperty(prefix+"keep_bottom")!=null)             this.keep_bottom=Double.parseDouble(properties.getProperty(prefix+"keep_bottom"));
		// Filtering clusters (removing those that fail
		if (properties.getProperty(prefix+"filter_headlessOK")!=null)       this.filter_headlessOK=Boolean.parseBoolean(properties.getProperty(prefix+"filter_headlessOK"));
		if (properties.getProperty(prefix+"filter_en")!=null)               this.filter_en=Boolean.parseBoolean(properties.getProperty(prefix+"filter_en"));
		if (properties.getProperty(prefix+"filter_min_height")!=null)       this.filter_min_height=Integer.parseInt(properties.getProperty(prefix+"filter_min_height"));
		if (properties.getProperty(prefix+"filter_min_fraction")!=null)     this.filter_min_fraction=Double.parseDouble(properties.getProperty(prefix+"filter_min_fraction"));
		if (properties.getProperty(prefix+"filter_min_frac_height")!=null)  this.filter_min_frac_height=Double.parseDouble(properties.getProperty(prefix+"filter_min_frac_height"));
		if (properties.getProperty(prefix+"filter_min_frac_box")!=null)     this.filter_min_frac_box=Double.parseDouble(properties.getProperty(prefix+"filter_min_frac_box"));
		if (properties.getProperty(prefix+"filter_min_disparity")!=null)    this.filter_min_disparity=Double.parseDouble(properties.getProperty(prefix+"filter_min_disparity"));
		if (properties.getProperty(prefix+"filter_max_disparity")!=null)    this.filter_max_disparity=Double.parseDouble(properties.getProperty(prefix+"filter_max_disparity"));
		if (properties.getProperty(prefix+"poles_debug_level")!=null)       this.poles_debug_level=Integer.parseInt(properties.getProperty(prefix+"poles_debug_level"));
	}
	@Override
	public PoleProcessorParameters clone() { //  throws CloneNotSupportedException {
		PoleProcessorParameters ppp =       new PoleProcessorParameters();
		ppp.poles_en=                   this.poles_en;
		// Creating initial pole seeds (tops need to be detected and available in the available DSI
		ppp.seed_rsrength=              this.seed_rsrength;
		ppp.seed_down=                  this.seed_down;
		ppp.seed_aover=                 this.seed_aover;
		ppp.seed_rover=                 this.seed_rover;
		ppp.max_disparity=              this.max_disparity;
		// Seed (pole heads) tiles initial clustering
		ppp.max_dd=                     this.max_dd;
		ppp.max_dx=                     this.max_dx;
		ppp.max_dy=                     this.max_dy;
		ppp.ext_side=                   this.ext_side;
		ppp.ext_up=                     this.ext_up;
		ppp.ext_down=                   this.ext_down;
		ppp.pedestal_strength=          this.pedestal_strength;
		ppp.pedestal_disp_over=         this.pedestal_disp_over;
		// Cluster merge and split parameters
		ppp.merge_extended =            this.merge_extended;
		ppp.merge_rtolerance =          this.merge_rtolerance;
		ppp.split_min_dist =            this.split_min_dist;
		ppp.split_must_zero =           this.split_must_zero;
	    // Creating pole centerline and selection mask
		ppp.max_diff =                  this.max_diff;
		ppp.max_tilt =                  this.max_tilt;
		ppp.max_rms =                   this.max_rms;
		ppp.min_tiles =                 this.min_tiles;
		ppp.damp_tilt =                 this.damp_tilt;
		ppp.use_seed =                  this.use_seed;
		ppp.hwidth =                    this.hwidth;
		ppp.disp_aover =                this.disp_aover;
		ppp.disp_rover =                this.disp_rover;
		ppp.min_neibs =                 this.min_neibs;
		// More tile filtering: trim unoccupied bottoms and check that the cluster head is separated on both right and left sides
		ppp.trim_bottoms =              this.trim_bottoms;
		ppp.extra_margin =              this.extra_margin;
 		ppp.sep_min_strength =          this.sep_min_strength;
 		ppp.sep_disp_adiff =            this.sep_disp_adiff;
 		ppp.sep_disp_rdiff =            this.sep_disp_rdiff;
 		//Measurements/re-measurements parameters with notch filtering to emphasize vertical lines and suppress random texture of the background
 		ppp.remeasure_hheight =         this.remeasure_hheight;
 		ppp.high_disp_tolerance =       this.high_disp_tolerance;
 		ppp.low_disp_tolerance =        this.low_disp_tolerance;
 		ppp.min_strength =              this.min_strength;
 		ppp.use_fittest=                this.use_fittest;
		// Updating disparity using the measurement results
 		ppp.max_refines =               this.max_refines;
 		ppp.min_refine_diff =           this.min_refine_diff;
 		ppp.disparity_scale =           this.disparity_scale;
 		ppp.diff_power =                this.diff_power;
 		ppp.diff_offset =               this.diff_offset;
 		ppp.cut_bottom =                this.cut_bottom;
 		ppp.keep_bottom =               this.keep_bottom;
 		// Filtering clusters (removing those that fail
 		ppp.filter_headlessOK =         this.filter_headlessOK;
 		ppp.filter_en =                 this.filter_en;
 		ppp.filter_min_height =         this.filter_min_height;
 		ppp.filter_min_fraction =       this.filter_min_fraction;
 		ppp.filter_min_frac_height =    this.filter_min_frac_height;
 		ppp.filter_min_frac_box =       this.filter_min_frac_box;
 		ppp.filter_min_disparity =      this.filter_min_disparity;
 		ppp.filter_max_disparity =      this.filter_max_disparity;
 		ppp.poles_debug_level =         this.poles_debug_level;
		return ppp;
	}


}
