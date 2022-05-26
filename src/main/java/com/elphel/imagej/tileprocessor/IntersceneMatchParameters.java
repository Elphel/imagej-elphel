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
	public  boolean export_images =      true; // 16-slice images (same disparity, COMBO_DSN_INDX_DISP_FG and COMBO_DSN_INDX_DISP_BG_ALL,
	public  boolean show_images =        false;
	public  boolean show_ranges =        true;
	
	public  double  range_disparity_offset =   -0.08;
	public  double  range_min_strength = 0.5;
	public  double  range_max =       5000.0;
	

	// Some "AGC" to adjust how much to discard
	public  int     margin =                       1;      // do not use tiles if their centers are closer to the image edge
	public  int     sensor_mask_inter =           -1;      // bitmask of the sensors to use (-1 - all)
	public  boolean use_partial =                   true;  // find motion vectors for individual pairs, false - for sum only
	public  boolean run_poly =                      false; // not yet implemented
	public  double  centroid_radius =               4.0;   // 
	public  int     n_recenter =                    2;     // when cosine window, re-center window this many times
	// filtering motion vectors
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
		gd.addCheckbox ("Show exported images",                      this.show_images,
				"Display generated/saved images");
		gd.addCheckbox ("Show distances in meters",                  this.show_ranges,
				"Calculate strength, distance, X, and Y in meters");
		
		gd.addNumericField("Disparity at infinity",                  this.range_disparity_offset, 5,7,"pix",
				"Disparity at infinity - subtract from measured disparity when converting to ranges.");
		gd.addNumericField("Minimal strength for range calculation", this.range_min_strength, 5,7,"",
				"Disregard weaker results when measuring range.");
		gd.addNumericField("Maximal displayed range",                this.range_max, 5,7,"m",
				"Do not display extremely far objects.");
		
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
		this.show_ranges =              gd.getNextBoolean();
		this.range_disparity_offset =   gd.getNextNumber();
		this.range_min_strength =       gd.getNextNumber();
		this.range_max =                gd.getNextNumber();

		this.margin =             (int) gd.getNextNumber();
		this.sensor_mask_inter=   (int) gd.getNextNumber();
		this.use_partial =              gd.getNextBoolean();
		this.run_poly =                 gd.getNextBoolean();
		this.centroid_radius =          gd.getNextNumber();
		this.n_recenter =         (int) gd.getNextNumber();
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
		properties.setProperty(prefix+"show_ranges",          this.show_ranges + "");       // boolean

		properties.setProperty(prefix+"range_disparity_offset",this.range_disparity_offset+""); // double
		properties.setProperty(prefix+"range_min_strength",   this.range_min_strength+"");  // double
		properties.setProperty(prefix+"range_max",            this.range_max+"");           // double
		
		properties.setProperty(prefix+"margin",               this.margin+"");              // int
		properties.setProperty(prefix+"sensor_mask_inter",    this.sensor_mask_inter+"");   // int
		properties.setProperty(prefix+"use_partial",          this.use_partial+"");         // boolean
		properties.setProperty(prefix+"run_poly",             this.run_poly+"");            // boolean
		properties.setProperty(prefix+"centroid_radius",      this.centroid_radius+"");     // double
		properties.setProperty(prefix+"n_recenter",           this.n_recenter+"");          // int
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
		if (properties.getProperty(prefix+"show_ranges")!=null)          this.show_images=Boolean.parseBoolean(properties.getProperty(prefix+"show_ranges"));
		if (properties.getProperty(prefix+"range_disparity_offset")!=null) this.range_disparity_offset=Double.parseDouble(properties.getProperty(prefix+"range_disparity_offset"));
		if (properties.getProperty(prefix+"range_min_strength")!=null)   this.range_min_strength=Double.parseDouble(properties.getProperty(prefix+"range_min_strength"));
		if (properties.getProperty(prefix+"range_max")!=null)            this.range_max=Double.parseDouble(properties.getProperty(prefix+"range_max"));
		if (properties.getProperty(prefix+"margin")!=null)               this.margin=Integer.parseInt(properties.getProperty(prefix+"margin"));
		if (properties.getProperty(prefix+"sensor_mask_inter")!=null)    this.sensor_mask_inter=Integer.parseInt(properties.getProperty(prefix+"sensor_mask_inter"));
		if (properties.getProperty(prefix+"use_partial")!=null)          this.use_partial=Boolean.parseBoolean(properties.getProperty(prefix+"use_partial"));		
		if (properties.getProperty(prefix+"run_poly")!=null)             this.run_poly=Boolean.parseBoolean(properties.getProperty(prefix+"run_poly"));
		if (properties.getProperty(prefix+"centroid_radius")!=null)      this.centroid_radius=Double.parseDouble(properties.getProperty(prefix+"centroid_radius"));
		if (properties.getProperty(prefix+"n_recenter")!=null)           this.n_recenter=Integer.parseInt(properties.getProperty(prefix+"n_recenter"));
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
		imp.show_ranges           = this.show_ranges;
		imp.range_disparity_offset = this.range_disparity_offset;
		imp.range_min_strength    = this.range_min_strength;
		imp.range_max             = this.range_max;
		imp.margin                = this.margin;
		imp.sensor_mask_inter     = this.sensor_mask_inter;
		imp.use_partial           = this.use_partial;
		imp.run_poly              = this.run_poly;
		imp.centroid_radius       = this.centroid_radius;
		imp.n_recenter            = this.n_recenter;
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
