import java.util.Properties;

/**
 **
 ** BiQuadParameters - parameters defining Operation of a two quad camera rig
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  BiQuadParameters.java is free software: you can redistribute it and/or modify
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

public class BiQuadParameters {
	public boolean rig_mode_debug =            true;
	public boolean use_poly =                  true;
	public boolean use_xy_poly =               true;

	public double  min_poly_strength =         0.2;
	public double  min_xy_poly_strength =      1.0; // never
	public double  inf_min_strength_main =     0.12;
	public double  inf_min_strength_aux =      0.12;
	public double  inf_min_strength_rig =      0.25;
	public double  inf_max_disp_main =         0.15;
	public double  inf_max_disp_aux =          0.15;
	public double  inf_max_disp_rig =          0.2;    // maybe even higher (2.0) to lock to initially high mismatch
	public double  inf_neg_tolerance =         2.5;    // increase negative disparity for infinity tolerance
	public double  inf_weight =                0.7;    // weight of infinity measurements of all measurements

	public double  first_max_disp_main =       0.5;    // before refinement
	public double  first_max_disp_aux =        0.5;
	public double  first_max_disp_rig =        4.0;    //

	public int     num_refine_master =         5;      // Number of disparity refinement passes
	public int     num_refine_inter =          2;      // Number of disparity refinement passes
	public double  near_max_disp_main =        0.3;    // ater refinement
	public double  near_max_disp_aux =         0.2;
	public double  near_max_disp_rig =         0.5;    //

	public double  refine_min_strength =       0.1;    //
	public double  refine_tolerance =          0.03;    //

	public double  rig_disp_range =           10.0;    // maximal disparity to try to adjust auxiliary camera pose
	public int     rig_num_disp_steps =        2;      // number of disparity vaues to try - 1 - only infinity, 2 - infinity and near, ...
	public int     rig_adjust_full_cycles=     1;      // number of full rig adjustment steps (with re-creating tiles lists)
	public int     rig_adjust_short_cycles=    15;     // maximal number of full rig adjustment steps (with current tiles lists)
	public double  rig_adjust_short_threshold= 0.00001;// rms relative improvement to exit short cycle
	public boolean rig_adjust_orientation =    true;   //
	public boolean rig_adjust_zoom =           true;   //
	public boolean rig_adjust_angle =          true;   //
	public boolean rig_adjust_distance =       false;  // distance between camera centers
	public boolean rig_adjust_forward =        false;  // aux camera forward from the principal plane (not implemented)
	public double  rig_correction_scale=       1.0;    // scale calculaated correction

	public int     min_new =                   100;    // Minimal number of he new tiles during rig refine
	public int     num_inf_refine =            20;     // Number of infinity refine passes
	public int     num_near_refine =           20;     // Number of non-infinity refine passes
	public double  min_trusted_strength =      0.1;//14// Minimal trusted combo strength;
	public double  trusted_tolerance  =        1.0;    // Trusted tolerance for small baseline camera(s)

	// rig LT (poor textured areas)
	public double  lt_min_disparity =          0.0;    // apply low texture to near objects
	public double  lt_trusted_strength =       0.2;    // strength sufficient without neighbors
	public double  lt_need_friends =           0.4;    // strength sufficient with neighbors support, fraction of lt_trusted_strength
	public double  lt_friends_diff =           0.15;   // pix difference to neighbors to be considered a match (TODO: use tilted)
	public double  lt_friends_rdiff =          0.04;   // Add per each disparity pixel
	public int     lt_min_friends_any =        2;      // minimal number of even weak friends
	public int     lt_min_friends_trusted =    2;      // minimal number of trusted (strong or already confirmed)
	public int     lt_friends_dist =           3;      // how far to look for friends
	public boolean lt_replace_lone =           true;   // try to overwrite lone weak
	public int     lt_extend_dist =            3;      // how far to extend around known tiles (probably should increase this value up to?
	// dealing with neighbors variance
	public double  lt_wsigma =                 1.0;    // Reduce influence of far neighbors with this Gaussian sigma
	public double  lt_max_asigma =             .15;    // Maximal acceptable standard deviation of the neighbors (remove, then add)
	public double  lt_max_rsigma =             .04;    // Additional standard deviation for each pixel of disparity (relative)


	public int     ml_hwidth =                 2;      // Half-width of the ML tiles to export (0-> 1x1, 1->3x3, 2 -> 5x5)
	public double  ml_disparity_sweep  =       2.0;    // Disparity sweep around ground truth, each side
	public int     ml_sweep_steps =            5;      // Number of disparity sweep steps

	public boolean ml_keep_aux =               true; // include auxiliary camera data in the ML output
	public boolean ml_keep_inter =             true; // include inter-camera correlation data in the ML output
	public boolean ml_keep_hor_vert =          true; // include combined horizontal and vertical pairs data in the ML output
	public boolean ml_keep_tbrl =              true; // include individual top, bottom, right, left pairs
	public boolean ml_keep_debug=              true; // include debug layer(s) data in the ML output
	public boolean ml_8bit=                    true; // output in 8-bit format (default - 32-bit TIFF
	public double  ml_limit_extrim =           0.00001; // ignore lowest and highest values when converting to 8 bpp
	public boolean ml_show_ml =                true; // show each generated MLoutput file
	public double  ml_fatzero =                0.05; // Use this value for correlation





	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addCheckbox    ("Debug rig/bi-camera functionality ",                              this.rig_mode_debug,"Enable debugging of the methods related to dual camera rig");
		gd.addCheckbox    ("Use poly for main/aux/rig correlations (false - CM)",             this.use_poly,"Use LMA/polynomial if correlation is strong enough");
		gd.addCheckbox    ("Use poly for rig X/Y mismatch measurements",                      this.use_xy_poly,"Use polynomial forX/Y offset measurements if correlation is strong enough");

		gd.addNumericField("Use poly mode for disparity if strength is greater than",                             this.min_poly_strength,  3,6,"",
				"Minimal correlation strength to use poly mode (below use CM)");
		gd.addNumericField("Use poly mode for X/Y rig correction  if strength is greater than",                   this.min_xy_poly_strength,  3,6,"",
				"Minimal correlation strength to use poly mode (below use CM) for measureking X/Y offset at infinity");
		gd.addNumericField("Minimal strength for main camera correlation to use for infinity rig adjustment",     this.inf_min_strength_main,  3,6,"",
				"Do not use tile for infinity adjustment if main correlation strength is less than");
		gd.addNumericField("Minimal strength for aux camera correlation to use for infinity rig adjustment",      this.inf_min_strength_aux,  3,6,"",
				"Do not use tile for infinity adjustment if aux correlation strength is less than");
		gd.addNumericField("Minimal strength for inter-camera correlation to use for infinity rig adjustment",    this.inf_min_strength_rig,  3,6,"",
				"Do not use tile for infinity adjustment if inter-camera correlation strength is less than");
		gd.addNumericField("Maximal absolute value of main camera disparity to use for infinity rig adjustment",  this.inf_max_disp_main,  3,6,"pix",
				"Do not use tile for infinity adjustment if absolute value of the main camera disparity is too high");
		gd.addNumericField("Maximal absolute value of aux camera disparity to use for infinity rig adjustment",   this.inf_max_disp_aux,  3,6,"pix",
				"Do not use tile for infinity adjustment if absolute value of the main camera disparity is too high");
		gd.addNumericField("Maximal absolute value of inter-camera disparity to use for infinity rig adjustment", this.inf_max_disp_rig,  3,6,"pix",
				"Do not use tile for infinity adjustment if absolute value of the inter-camera disparity is too high");
		gd.addNumericField("Loosen negative disparity tolerance for infinity",                                    this.inf_neg_tolerance,  3,6,"",
				"Allow farther negative than positive disparity tiles for infinity (only for main/rig pair)");
		gd.addNumericField("Weight of infinity measurements in all measurements",                                 this.inf_weight,  3,6,"",
				"Set importance of infinity matching among all measurements");


		gd.addNumericField("Maximal absolute value of main camera disparity difference near objects, before refinement", this.first_max_disp_main,  3,6,"pix",
				"First pass for near objects, before refinement pass");
		gd.addNumericField("Maximal absolute value of aux camera disparity difference near objects, before refinement",  this.first_max_disp_aux,  3,6,"pix",
				"First pass for near objects, before refinement pass");
		gd.addNumericField("Maximal absolute value of inter-camera disparity difference near objects, before refinement",this.first_max_disp_rig,  3,6,"pix",
				"First pass for near objects, before refinement pass");

		gd.addNumericField("Number of disparity refinement passes for the non-infinity measurements, master camera",       this.num_refine_master,  0,3,"",
				"Remeasure after setting per-tile expected disparity (phase shift before correlation)");
		gd.addNumericField("Number of disparity refinement passes for the non-infinity measurements, inter-cameras",       this.num_refine_inter,  0,3,"",
				"Remeasure after setting per-tile expected disparity (phase shift before correlation)");

		gd.addNumericField("Maximal absolute value of main camera disparity difference near objects, after refinement",   this.near_max_disp_main,  3,6,"pix",
				"After refinement pass (target disparity set to the measured one");
		gd.addNumericField("Maximal absolute value of aux camera disparity difference near objects, after refinement",    this.near_max_disp_aux,  3,6,"pix",
				"After refinement pass (target disparity set to the measured one");
		gd.addNumericField("Maximal absolute value of inter-camera disparity difference near objects, after refinement",  this.near_max_disp_rig,  3,6,"pix",
				"After refinement pass (target disparity set to the measured one");

		gd.addNumericField("Minimal correlation strength to try to refine disparity",                                     this.refine_min_strength,  3,6,"",
				"Re-scan tile with updated expected disparity if previously measured strength was above this");
		gd.addNumericField("Refine tolerance - do not try to refine if residual disparity is lower",                      this.refine_tolerance,  3,6,"pix",
				"Re-scan tile only if residual disparity absolute value is lower that this");


		gd.addNumericField("Maximal disparity to try for the rig adjustment",                                     this.rig_disp_range,  3,6,"pix",
				"Maximal disparity of the main camera when adjusting auxiliary camera pose");
		gd.addNumericField("Number of disparity values to try between zero (infinity) and the maximal one",       this.rig_num_disp_steps,  0,3,"",
				"1 - only infinity, 2 - infinity and maximal, 3 - infinuty, half-maximal and maximal");
		gd.addNumericField("Number of full rig adjustment cycles",                                                this.rig_adjust_full_cycles,  0,3,"",
				"Includes full-frame scanning and re-creating the sample tile list");
		gd.addNumericField("Number of short rig adjustment cycles",                                               this.rig_adjust_short_cycles,  0,3,"",
				"Re-scan only previously selected tiles");
		gd.addNumericField("Minimal RMS improvement after short cycle to exit",                                   this.rig_adjust_short_threshold,  3,6,"",
				"Short cycle iteration finish after relative RMS improvement is less than this of maximal number of steps is exceeded");
		gd.addCheckbox    ("Adjust orientation (azimuth, tilt, roll) of the auxiliary camera",                    this.rig_adjust_orientation,"3 angles, most basic adjustment");
		gd.addCheckbox    ("Adjust orientation relative zoom of the auxiliary camera",                            this.rig_adjust_zoom,"Not likely to change, may be frozen");
		gd.addCheckbox    ("Adjust position of the aux camera in the main camera pricipal plane",                 this.rig_adjust_angle,"0 - exactly  to the right, positive - higher than main");
		gd.addCheckbox    ("Adjust aux camera distance from the main",                                            this.rig_adjust_distance,"Normally not practical, rely on direct measurement");
		gd.addCheckbox    ("Adjust aux camera distance from the main principal plane",                            this.rig_adjust_forward,"Not implemented, assumed zero");
		gd.addNumericField("Scale calculated correction before applying",                                         this.rig_correction_scale,  3,6,"",
				"Debug feature");

		gd.addNumericField("Minimal number of he new tiles during rig refine",                                    this.min_new,  0,3,"",
				"Exit from refine infinity cycle if number of new refine candidates is less than this number");
		gd.addNumericField("Number of infinity refine passes",                                                    this.num_inf_refine,  0,3,"",
				"Exit from refine non- infinity cycle if number of new refine candidates is less than this number");
		gd.addNumericField("Number of non-infinity refine passes",                                                this.num_near_refine,  0,3,"",
				"Re-scan only previously selected tiles");
		gd.addNumericField("Minimal trusted combo strength",                                                      this.min_trusted_strength,  3,6,"",
				"Combo strength is cubic root of the product of main, aux and intyer correlation strengths");
		gd.addNumericField("Trusted tolerance for small baseline camera(s)",                                      this.trusted_tolerance,  3,6,"",
				"When downscaling valid residual disparity from the most sensitive inter-camera, do not reduce it to be lower than this");

        gd.addTab("Rig LT","Deal with the low textured areas");

		gd.addNumericField("Apply low texture to near objects (disparity above)",                                 this.lt_min_disparity,  3,6,"pix",
				"Handle low textured objects with disparity (main camera pixels) above this threshold");
		gd.addNumericField("Inter-camera correlation strength sufficient without neighbors",                      this.lt_trusted_strength,  3,6,"",
				"Consider such tiles valid regardless of surrounding tiles");
		gd.addNumericField("Strength sufficient with neighbors support, fraction of the trusted (above)",         this.lt_need_friends,  3,6,"",
				"Tiles above theis inter-camera strength may be valid if they have same disparity neighbors");
		gd.addNumericField("Maximal disparity difference to neighbors to be considered a match",                  this.lt_friends_diff,  3,6,"pix",
				"TODO: use non fronto parallel surfaces (e.g. horizontal) ");
		gd.addNumericField("Add to allowed disparity difference per each pixel of dipsarity",                     this.lt_friends_rdiff,  3,6,"pix/pix",
				"Additional relative disparity difference");

		gd.addNumericField("Minimal number of even weak friends sufficient for support",                          this.lt_min_friends_any,  0,3,"",
				"Tiles having this number of same-disparity neighbors with strength above minimal are considered trusted");
		gd.addNumericField("Minimal number of trusted (strong or already confirmed)",                             this.lt_min_friends_trusted,  0,3,"",
				"Tiles having this number of same-disparity neighbors with trusted strength or already confirmed by other neighbors are trusted too");
		gd.addNumericField("How far to look for friends",                                                         this.lt_friends_dist,  0,3,"",
				"Look for matching tiles inside a square of 2*<this value>+1 around this tile");
		gd.addCheckbox    ("Discard lone weak",                                                                   this.lt_replace_lone,
				"Discard lone (no matching neighbors) tiles that have strength below trusted");
		gd.addNumericField("How far to extend around known tiles",                                                this.lt_extend_dist,  0,3,"",
				"Try new tiles that have neighbors within a square around it");
		gd.addNumericField("Reduce influence of far neighbors with this Gaussian sigma",                          this.lt_wsigma,  4,6,"pix",
				"Neighbor weight is multipled by a Gaussian with this sigma around the new tile");


		gd.addNumericField("Maximal acceptable standard deviation of the neighbors (remove, then add)",           this.lt_max_asigma,  4,6,"pix",
				"When multiple neighbors have different disparity, try to remove outliers, then add some back");
		gd.addNumericField("Additional standard deviation for each pixel of disparity (relative)",                this.lt_max_rsigma,  4,6,"pix/pix",
				"Loosen sigma requirements for high disparity by adding this value for each 1 pixel of disparity");

        gd.addTab("ML","Parameters related to the ML files generation for the dual-quad camera rig");

		gd.addNumericField("Half-width of the ML tiles to export (0-> 1x1, 1->3x3, 2 -> 5x5)",                    this.ml_hwidth,  0,3,"",
				"Amount of data to export to the ML system");
		gd.addNumericField("Disparity sweep around ground truth, each side",                                      this.ml_disparity_sweep,  3,6,"",
				"Sweep symmetrically target disparity around the ground truth disparity, each side");
		gd.addNumericField("Number of target disparity sweep steps",                                                               this.ml_sweep_steps,  0,3,"",
				"Generate this many files for each file set. Each tile results depend on the target disparity and this tile data, do not depend on other tiles target disparity");
		gd.addCheckbox    ("Include auxiliary camera data in the ML output",                                      this.ml_keep_aux,
				"ML output will have the second set of the layers for the auxiliary camera. Disparity values should be scaled for the camera baseline");
		gd.addCheckbox    ("Keep inter-camera correlation data",                                                  this.ml_keep_inter,
				"Inter-camera correlation data has only one layer (and one correlation pair). It is used to generate ground truth data. Usable disparity range (measured in the main camera pixels) is ~1/5 of teh main camera");
		gd.addCheckbox    ("Keep individual top, bottom, right, and left pairs",                                  this.ml_keep_tbrl,
				"Each of these two layers per camera are calculated from a pair of top/bottom and left/right pairs. Can possibly be used instead of originals to reduce amount of input data");
		gd.addCheckbox    ("Keep combined horizonta/vertical pairs",                                              this.ml_keep_hor_vert,
				"Individual horizontal and vertical pairs (4 total). Can be replaced by two combined (horizontal+vertical) ones");

		gd.addCheckbox    ("Keep debug layer(s)",                                                                 this.ml_keep_debug,
				"Keep additional (debug) layers that may change for different file versions");
		gd.addCheckbox    ("Use 8 bpp TIFF (default - 32 bpp)",                                                   this.ml_8bit,
				"Reduce file size by lowering bpp");
		gd.addNumericField("When converting to 8bpp, limit fraction of extreme values",                           1E6 * this.ml_limit_extrim,  1,8,"ppm",
				"Use values histogram to find min/max values, ignoring(limiting) this fraction (parts per million) of pixels at both extremes");
		gd.addCheckbox    ("Show each generated ML file",                                                         this.ml_show_ml,
				"Use only for small number of generated files to reduce memory usage");
		gd.addNumericField("Use this phase correlation fat zero when generating ML files",                        this.ml_fatzero,  5,8,"",
				"Replace normal fat zero value used elsethere");

	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.rig_mode_debug=                gd.getNextBoolean();
		this.use_poly=                      gd.getNextBoolean();
		this.use_xy_poly=                   gd.getNextBoolean();
		this.min_poly_strength=             gd.getNextNumber();
		this.min_xy_poly_strength=          gd.getNextNumber();
		this.inf_min_strength_main=         gd.getNextNumber();
		this.inf_min_strength_aux=          gd.getNextNumber();
		this.inf_min_strength_rig=          gd.getNextNumber();
		this.inf_max_disp_main=             gd.getNextNumber();
		this.inf_max_disp_aux=              gd.getNextNumber();
		this.inf_max_disp_rig=              gd.getNextNumber();
		this.inf_neg_tolerance=             gd.getNextNumber();
		this.inf_weight=                    gd.getNextNumber();
		this.first_max_disp_main=           gd.getNextNumber();
		this.first_max_disp_aux=            gd.getNextNumber();
		this.first_max_disp_rig=            gd.getNextNumber();
		this.num_refine_master=       (int) gd.getNextNumber();
		this.num_refine_inter=        (int) gd.getNextNumber();
		this.near_max_disp_main=            gd.getNextNumber();
		this.near_max_disp_aux=             gd.getNextNumber();
		this.near_max_disp_rig=             gd.getNextNumber();

		this.refine_min_strength=           gd.getNextNumber();
		this.refine_tolerance=              gd.getNextNumber();

		this.rig_disp_range=                gd.getNextNumber();
		this.rig_num_disp_steps=      (int) gd.getNextNumber();
		this.rig_adjust_full_cycles=  (int) gd.getNextNumber();
		this.rig_adjust_short_cycles= (int) gd.getNextNumber();
		this.rig_adjust_short_threshold=    gd.getNextNumber();
		this.rig_adjust_orientation=        gd.getNextBoolean();
		this.rig_adjust_zoom=               gd.getNextBoolean();
		this.rig_adjust_angle=              gd.getNextBoolean();
		this.rig_adjust_distance=           gd.getNextBoolean();
		this.rig_adjust_forward=            gd.getNextBoolean();
		this.rig_correction_scale=          gd.getNextNumber();

		this.min_new=                 (int) gd.getNextNumber();
		this.num_inf_refine=          (int) gd.getNextNumber();
		this.num_near_refine=         (int) gd.getNextNumber();
		this.min_trusted_strength=          gd.getNextNumber();
		this.trusted_tolerance=             gd.getNextNumber();

		this.lt_min_disparity=              gd.getNextNumber();
		this.lt_trusted_strength=           gd.getNextNumber();
		this.lt_need_friends=               gd.getNextNumber();
		this.lt_friends_diff=               gd.getNextNumber();
		this.lt_friends_rdiff=              gd.getNextNumber();
		this.lt_min_friends_any=      (int) gd.getNextNumber();
		this.lt_min_friends_trusted=  (int) gd.getNextNumber();
		this.lt_friends_dist=         (int) gd.getNextNumber();
		this.lt_replace_lone=               gd.getNextBoolean();
		this.lt_extend_dist=          (int) gd.getNextNumber();
		this.lt_wsigma=                     gd.getNextNumber();
		this.lt_max_asigma=                 gd.getNextNumber();
		this.lt_max_rsigma=                 gd.getNextNumber();

		this.ml_hwidth=               (int) gd.getNextNumber();
		this.ml_disparity_sweep=            gd.getNextNumber();
		this.ml_sweep_steps=          (int) gd.getNextNumber();
		this.ml_keep_aux=                   gd.getNextBoolean();
		this.ml_keep_inter=                 gd.getNextBoolean();
		this.ml_keep_tbrl=                  gd.getNextBoolean();
		this.ml_keep_hor_vert=              gd.getNextBoolean();
		this.ml_keep_debug=                 gd.getNextBoolean();
		this.ml_8bit=                       gd.getNextBoolean();
		this.ml_limit_extrim=               gd.getNextNumber() * 1E-6;
		this.ml_show_ml=                    gd.getNextBoolean();
		this.ml_fatzero=                    gd.getNextNumber();
	}

	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"rig_mode_debug",            this.rig_mode_debug+"");
		properties.setProperty(prefix+"use_poly",                  this.use_poly+"");
		properties.setProperty(prefix+"use_xy_poly",               this.use_xy_poly+"");
		properties.setProperty(prefix+"min_poly_strength",         this.min_poly_strength+"");
		properties.setProperty(prefix+"min_xy_poly_strength",      this.min_xy_poly_strength+"");
		properties.setProperty(prefix+"inf_min_strength_main",     this.inf_min_strength_main+"");
		properties.setProperty(prefix+"inf_min_strength_aux",      this.inf_min_strength_aux+"");
		properties.setProperty(prefix+"inf_min_strength_rig",      this.inf_min_strength_rig+"");
		properties.setProperty(prefix+"inf_max_disp_main",         this.inf_max_disp_main+"");
		properties.setProperty(prefix+"inf_max_disp_aux",          this.inf_max_disp_aux+"");
		properties.setProperty(prefix+"inf_max_disp_rig",          this.inf_max_disp_rig+"");
		properties.setProperty(prefix+"inf_neg_tolerance",         this.inf_neg_tolerance+"");
		properties.setProperty(prefix+"inf_weight",                this.inf_weight+"");


		properties.setProperty(prefix+"first_max_disp_main",       this.first_max_disp_main+"");
		properties.setProperty(prefix+"first_max_disp_aux",        this.first_max_disp_aux+"");
		properties.setProperty(prefix+"first_max_disp_rig",        this.first_max_disp_rig+"");
		properties.setProperty(prefix+"num_refine_master",         this.num_refine_master+"");
		properties.setProperty(prefix+"num_refine_inter",          this.num_refine_inter+"");
		properties.setProperty(prefix+"near_max_disp_main",        this.near_max_disp_main+"");
		properties.setProperty(prefix+"near_max_disp_aux",         this.near_max_disp_aux+"");
		properties.setProperty(prefix+"near_max_disp_rig",         this.near_max_disp_rig+"");

		properties.setProperty(prefix+"refine_min_strength",       this.refine_min_strength+"");
		properties.setProperty(prefix+"refine_tolerance",          this.refine_tolerance+"");

		properties.setProperty(prefix+"rig_disp_range",            this.rig_disp_range+"");
		properties.setProperty(prefix+"rig_num_disp_steps",        this.rig_num_disp_steps+"");
		properties.setProperty(prefix+"rig_adjust_full_cycles",    this.rig_adjust_full_cycles+"");
		properties.setProperty(prefix+"rig_adjust_short_cycles",   this.rig_adjust_short_cycles+"");
		properties.setProperty(prefix+"rig_adjust_short_threshold",this.rig_adjust_short_threshold+"");
		properties.setProperty(prefix+"rig_adjust_orientation",    this.rig_adjust_orientation+"");
		properties.setProperty(prefix+"rig_adjust_zoom",           this.rig_adjust_zoom+"");
		properties.setProperty(prefix+"rig_adjust_angle",          this.rig_adjust_angle+"");
		properties.setProperty(prefix+"rig_adjust_distance",       this.rig_adjust_distance+"");
		properties.setProperty(prefix+"rig_adjust_forward",        this.rig_adjust_forward+"");
		properties.setProperty(prefix+"rig_correction_scale",      this.rig_correction_scale+"");

		properties.setProperty(prefix+"min_new",                   this.min_new+"");
		properties.setProperty(prefix+"num_inf_refine",            this.num_inf_refine+"");
		properties.setProperty(prefix+"num_near_refine",           this.num_near_refine+"");
		properties.setProperty(prefix+"min_trusted_strength",      this.min_trusted_strength+"");
		properties.setProperty(prefix+"trusted_tolerance",         this.trusted_tolerance+"");

		properties.setProperty(prefix+"lt_min_disparity",          this.lt_min_disparity+"");
		properties.setProperty(prefix+"lt_trusted_strength",       this.lt_trusted_strength+"");
		properties.setProperty(prefix+"lt_need_friends",           this.lt_need_friends+"");
		properties.setProperty(prefix+"lt_friends_diff",           this.lt_friends_diff+"");
		properties.setProperty(prefix+"lt_friends_rdiff",          this.lt_friends_rdiff+"");
		properties.setProperty(prefix+"lt_min_friends_any",        this.lt_min_friends_any+"");
		properties.setProperty(prefix+"lt_min_friends_trusted",    this.lt_min_friends_trusted+"");
		properties.setProperty(prefix+"lt_friends_dist",           this.lt_friends_dist+"");
		properties.setProperty(prefix+"lt_replace_lone",           this.lt_replace_lone+"");
		properties.setProperty(prefix+"lt_extend_dist",            this.lt_extend_dist+"");
		properties.setProperty(prefix+"lt_wsigma",                 this.lt_wsigma+"");
		properties.setProperty(prefix+"lt_max_asigma",             this.lt_max_asigma+"");
		properties.setProperty(prefix+"lt_max_rsigma",             this.lt_max_rsigma+"");

		properties.setProperty(prefix+"ml_hwidth",                 this.ml_hwidth+"");
		properties.setProperty(prefix+"ml_disparity_sweep",        this.ml_disparity_sweep+"");
		properties.setProperty(prefix+"ml_sweep_steps",            this.ml_sweep_steps+"");
		properties.setProperty(prefix+"ml_keep_aux",               this.ml_keep_aux+"");
		properties.setProperty(prefix+"ml_keep_inter",             this.ml_keep_inter+"");
		properties.setProperty(prefix+"ml_keep_tbrl",              this.ml_keep_tbrl+"");
		properties.setProperty(prefix+"ml_keep_hor_vert",          this.ml_keep_hor_vert+"");
		properties.setProperty(prefix+"ml_keep_debug",             this.ml_keep_debug+"");
		properties.setProperty(prefix+"ml_8bit",                   this.ml_8bit+"");
		properties.setProperty(prefix+"ml_limit_extrim",           this.ml_limit_extrim+"");
		properties.setProperty(prefix+"ml_show_ml",                this.ml_show_ml+"");
		properties.setProperty(prefix+"ml_fatzero",                this.ml_fatzero+"");
	}
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"rig_mode_debug")!=null)        this.rig_mode_debug=Boolean.parseBoolean(properties.getProperty(prefix+"rig_mode_debug"));
		if (properties.getProperty(prefix+"use_poly")!=null)              this.use_poly=Boolean.parseBoolean(properties.getProperty(prefix+"use_poly"));
		if (properties.getProperty(prefix+"use_xy_poly")!=null)           this.use_xy_poly=Boolean.parseBoolean(properties.getProperty(prefix+"use_xy_poly"));
		if (properties.getProperty(prefix+"min_poly_strength")!=null)     this.min_poly_strength=Double.parseDouble(properties.getProperty(prefix+"min_poly_strength"));
		if (properties.getProperty(prefix+"min_xy_poly_strength")!=null)  this.min_xy_poly_strength=Double.parseDouble(properties.getProperty(prefix+"min_xy_poly_strength"));
		if (properties.getProperty(prefix+"inf_min_strength_main")!=null) this.inf_min_strength_main=Double.parseDouble(properties.getProperty(prefix+"inf_min_strength_main"));
		if (properties.getProperty(prefix+"inf_min_strength_aux")!=null)  this.inf_min_strength_aux=Double.parseDouble(properties.getProperty(prefix+"inf_min_strength_aux"));
		if (properties.getProperty(prefix+"inf_min_strength_rig")!=null)  this.inf_min_strength_rig=Double.parseDouble(properties.getProperty(prefix+"inf_min_strength_rig"));
		if (properties.getProperty(prefix+"inf_max_disp_main")!=null)     this.inf_max_disp_main=Double.parseDouble(properties.getProperty(prefix+"inf_max_disp_main"));
		if (properties.getProperty(prefix+"inf_max_disp_aux")!=null)      this.inf_max_disp_aux=Double.parseDouble(properties.getProperty(prefix+"inf_max_disp_aux"));
		if (properties.getProperty(prefix+"inf_max_disp_rig")!=null)      this.inf_max_disp_rig=Double.parseDouble(properties.getProperty(prefix+"inf_max_disp_rig"));
		if (properties.getProperty(prefix+"inf_neg_tolerance")!=null)     this.inf_neg_tolerance=Double.parseDouble(properties.getProperty(prefix+"inf_neg_tolerance"));
		if (properties.getProperty(prefix+"inf_weight")!=null)            this.inf_weight=Double.parseDouble(properties.getProperty(prefix+"inf_weight"));

		if (properties.getProperty(prefix+"first_max_disp_main")!=null)     this.first_max_disp_main=Double.parseDouble(properties.getProperty(prefix+"first_max_disp_main"));
		if (properties.getProperty(prefix+"first_max_disp_aux")!=null)      this.first_max_disp_aux=Double.parseDouble(properties.getProperty(prefix+"first_max_disp_aux"));
		if (properties.getProperty(prefix+"first_max_disp_rig")!=null)      this.first_max_disp_rig=Double.parseDouble(properties.getProperty(prefix+"first_max_disp_rig"));
		if (properties.getProperty(prefix+"num_refine_master")!=null)       this.num_refine_master=Integer.parseInt(properties.getProperty(prefix+"num_refine_master"));
		if (properties.getProperty(prefix+"num_refine_inter")!=null)        this.num_refine_inter=Integer.parseInt(properties.getProperty(prefix+"num_refine_inter"));
		if (properties.getProperty(prefix+"near_max_disp_main")!=null)      this.near_max_disp_main=Double.parseDouble(properties.getProperty(prefix+"near_max_disp_main"));
		if (properties.getProperty(prefix+"near_max_disp_aux")!=null)       this.near_max_disp_aux=Double.parseDouble(properties.getProperty(prefix+"near_max_disp_aux"));
		if (properties.getProperty(prefix+"near_max_disp_rig")!=null)       this.near_max_disp_rig=Double.parseDouble(properties.getProperty(prefix+"near_max_disp_rig"));
		if (properties.getProperty(prefix+"refine_min_strength")!=null)     this.refine_min_strength=Double.parseDouble(properties.getProperty(prefix+"refine_min_strength"));
		if (properties.getProperty(prefix+"refine_tolerance")!=null)        this.refine_tolerance=Double.parseDouble(properties.getProperty(prefix+"refine_tolerance"));
		if (properties.getProperty(prefix+"rig_disp_range")!=null)          this.rig_disp_range=Double.parseDouble(properties.getProperty(prefix+"rig_disp_range"));
		if (properties.getProperty(prefix+"rig_num_disp_steps")!=null)      this.rig_num_disp_steps=Integer.parseInt(properties.getProperty(prefix+"rig_num_disp_steps"));
		if (properties.getProperty(prefix+"rig_adjust_full_cycles")!=null)  this.rig_adjust_full_cycles=Integer.parseInt(properties.getProperty(prefix+"rig_adjust_full_cycles"));
		if (properties.getProperty(prefix+"rig_adjust_short_cycles")!=null) this.rig_adjust_short_cycles=Integer.parseInt(properties.getProperty(prefix+"rig_adjust_short_cycles"));
		if (properties.getProperty(prefix+"rig_adjust_short_threshold")!=null) this.rig_adjust_short_threshold=Double.parseDouble(properties.getProperty(prefix+"rig_adjust_short_threshold"));
		if (properties.getProperty(prefix+"rig_adjust_orientation")!=null)     this.rig_adjust_orientation=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_orientation"));
		if (properties.getProperty(prefix+"rig_adjust_zoom")!=null)         this.rig_adjust_zoom=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_zoom"));
		if (properties.getProperty(prefix+"rig_adjust_angle")!=null)        this.rig_adjust_angle=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_angle"));
		if (properties.getProperty(prefix+"rig_adjust_distance")!=null)     this.rig_adjust_distance=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_distance"));
		if (properties.getProperty(prefix+"rig_adjust_forward")!=null)      this.rig_adjust_forward=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_forward"));
		if (properties.getProperty(prefix+"rig_correction_scale")!=null)    this.rig_correction_scale=Double.parseDouble(properties.getProperty(prefix+"rig_correction_scale"));

		if (properties.getProperty(prefix+"min_new")!=null)                 this.min_new=Integer.parseInt(properties.getProperty(prefix+"min_new"));
		if (properties.getProperty(prefix+"num_inf_refine")!=null)          this.num_inf_refine=Integer.parseInt(properties.getProperty(prefix+"num_inf_refine"));
		if (properties.getProperty(prefix+"num_near_refine")!=null)         this.num_near_refine=Integer.parseInt(properties.getProperty(prefix+"num_near_refine"));
		if (properties.getProperty(prefix+"min_trusted_strength")!=null)    this.min_trusted_strength=Double.parseDouble(properties.getProperty(prefix+"min_trusted_strength"));
		if (properties.getProperty(prefix+"trusted_tolerance")!=null)       this.trusted_tolerance=Double.parseDouble(properties.getProperty(prefix+"trusted_tolerance"));
		if (properties.getProperty(prefix+"ml_hwidth")!=null)               this.ml_hwidth=Integer.parseInt(properties.getProperty(prefix+"ml_hwidth"));

		if (properties.getProperty(prefix+"lt_min_disparity")!=null)        this.lt_min_disparity=Double.parseDouble(properties.getProperty(prefix+"lt_min_disparity"));
		if (properties.getProperty(prefix+"lt_trusted_strength")!=null)     this.lt_trusted_strength=Double.parseDouble(properties.getProperty(prefix+"lt_trusted_strength"));
		if (properties.getProperty(prefix+"lt_need_friends")!=null)         this.lt_need_friends=Double.parseDouble(properties.getProperty(prefix+"lt_need_friends"));
		if (properties.getProperty(prefix+"lt_friends_diff")!=null)         this.lt_friends_diff=Double.parseDouble(properties.getProperty(prefix+"lt_friends_diff"));
		if (properties.getProperty(prefix+"lt_friends_rdiff")!=null)        this.lt_friends_rdiff=Double.parseDouble(properties.getProperty(prefix+"lt_friends_rdiff"));
		if (properties.getProperty(prefix+"lt_min_friends_any")!=null)      this.lt_min_friends_any=Integer.parseInt(properties.getProperty(prefix+"lt_min_friends_any"));
		if (properties.getProperty(prefix+"lt_min_friends_trusted")!=null)  this.lt_min_friends_trusted=Integer.parseInt(properties.getProperty(prefix+"lt_min_friends_trusted"));
		if (properties.getProperty(prefix+"lt_friends_dist")!=null)         this.lt_friends_dist=Integer.parseInt(properties.getProperty(prefix+"lt_friends_dist"));
		if (properties.getProperty(prefix+"lt_replace_lone")!=null)         this.lt_replace_lone=Boolean.parseBoolean(properties.getProperty(prefix+"lt_replace_lone"));
		if (properties.getProperty(prefix+"lt_extend_dist")!=null)          this.lt_extend_dist=Integer.parseInt(properties.getProperty(prefix+"lt_extend_dist"));
		if (properties.getProperty(prefix+"lt_wsigma")!=null)               this.lt_wsigma=Double.parseDouble(properties.getProperty(prefix+"lt_wsigma"));
		if (properties.getProperty(prefix+"lt_max_asigma")!=null)           this.lt_max_asigma=Double.parseDouble(properties.getProperty(prefix+"lt_max_sigma"));
		if (properties.getProperty(prefix+"lt_max_rsigma")!=null)           this.lt_max_rsigma=Double.parseDouble(properties.getProperty(prefix+"lt_max_sigma"));

		if (properties.getProperty(prefix+"ml_disparity_sweep")!=null)      this.ml_disparity_sweep=Double.parseDouble(properties.getProperty(prefix+"ml_disparity_sweep"));
		if (properties.getProperty(prefix+"ml_sweep_steps")!=null)          this.ml_sweep_steps=Integer.parseInt(properties.getProperty(prefix+"ml_sweep_steps"));
		if (properties.getProperty(prefix+"ml_keep_aux")!=null)             this.ml_keep_aux=Boolean.parseBoolean(properties.getProperty(prefix+"ml_keep_aux"));
		if (properties.getProperty(prefix+"ml_keep_inter")!=null)           this.ml_keep_inter=Boolean.parseBoolean(properties.getProperty(prefix+"ml_keep_inter"));
		if (properties.getProperty(prefix+"ml_keep_tbrl")!=null)            this.ml_keep_tbrl=Boolean.parseBoolean(properties.getProperty(prefix+"ml_keep_tbrl"));
		if (properties.getProperty(prefix+"ml_keep_hor_vert")!=null)        this.ml_keep_hor_vert=Boolean.parseBoolean(properties.getProperty(prefix+"ml_keep_hor_vert"));
		if (properties.getProperty(prefix+"ml_keep_debug")!=null)           this.ml_keep_debug=Boolean.parseBoolean(properties.getProperty(prefix+"ml_keep_debug"));
		if (properties.getProperty(prefix+"ml_8bit")!=null)                 this.ml_8bit=Boolean.parseBoolean(properties.getProperty(prefix+"ml_8bit"));
		if (properties.getProperty(prefix+"ml_limit_extrim")!=null)         this.ml_limit_extrim=Double.parseDouble(properties.getProperty(prefix+"ml_limit_extrim"));
		if (properties.getProperty(prefix+"ml_show_ml")!=null)              this.ml_show_ml=Boolean.parseBoolean(properties.getProperty(prefix+"ml_show_ml"));
		if (properties.getProperty(prefix+"ml_fatzero")!=null)              this.ml_fatzero=Double.parseDouble(properties.getProperty(prefix+"ml_fatzero"));
	}
	@Override
	public BiQuadParameters clone() throws CloneNotSupportedException {
		BiQuadParameters bqp =       new BiQuadParameters();
		bqp.rig_mode_debug=             this.rig_mode_debug;
		bqp.use_poly=                   this.use_poly;
		bqp.use_xy_poly=                this.use_xy_poly;
		bqp.min_poly_strength =         this.min_poly_strength;
		bqp.min_xy_poly_strength =      this.min_xy_poly_strength;
		bqp.inf_min_strength_main =     this.inf_min_strength_main;
		bqp.inf_min_strength_aux =      this.inf_min_strength_aux;
		bqp.inf_min_strength_rig =      this.inf_min_strength_rig;
		bqp.inf_max_disp_main =         this.inf_max_disp_main;
		bqp.inf_max_disp_aux =          this.inf_max_disp_aux;
		bqp.inf_max_disp_rig =          this.inf_max_disp_rig;
		bqp.inf_neg_tolerance =         this.inf_neg_tolerance;
		bqp.inf_weight =                this.inf_weight;

		bqp.first_max_disp_main=        this.first_max_disp_main;
		bqp.first_max_disp_aux=         this.first_max_disp_aux;
		bqp.first_max_disp_rig=         this.first_max_disp_rig;
		bqp.num_refine_master=          this.num_refine_master;
		bqp.num_refine_inter=           this.num_refine_inter;
		bqp.near_max_disp_main=         this.near_max_disp_main;
		bqp.near_max_disp_aux=          this.near_max_disp_aux;
		bqp.near_max_disp_rig=          this.near_max_disp_rig;
		bqp.refine_min_strength=        this.refine_min_strength;
		bqp.refine_tolerance=           this.refine_tolerance;
		bqp.rig_disp_range=             this.rig_disp_range;
		bqp.rig_num_disp_steps=         this.rig_num_disp_steps;
		bqp.rig_adjust_full_cycles=     this.rig_adjust_full_cycles;
		bqp.rig_adjust_short_cycles=    this.rig_adjust_short_cycles;
		bqp.rig_adjust_short_threshold= this.rig_adjust_short_threshold;
		bqp.rig_adjust_orientation=     this.rig_adjust_orientation;
		bqp.rig_adjust_zoom=            this.rig_adjust_zoom;
		bqp.rig_adjust_angle=           this.rig_adjust_angle;
		bqp.rig_adjust_distance=        this.rig_adjust_distance;
		bqp.rig_adjust_forward=         this.rig_adjust_forward;
		bqp.rig_correction_scale=       this.rig_correction_scale;

		bqp.min_new=                    this.min_new;
		bqp.num_inf_refine=             this.num_inf_refine;
		bqp.num_near_refine=            this.num_near_refine;
		bqp.min_trusted_strength=       this.min_trusted_strength;
		bqp.trusted_tolerance=          this.trusted_tolerance;

		bqp.lt_min_disparity=           this.lt_min_disparity;
		bqp.lt_trusted_strength=        this.lt_trusted_strength;
		bqp.lt_need_friends=            this.lt_need_friends;
		bqp.lt_friends_diff=            this.lt_friends_diff;
		bqp.lt_friends_rdiff=           this.lt_friends_rdiff;
		bqp.lt_min_friends_any=         this.lt_min_friends_any;
		bqp.lt_min_friends_trusted=     this.lt_min_friends_trusted;
		bqp.lt_friends_dist=            this.lt_friends_dist;
		bqp.lt_replace_lone=            this.lt_replace_lone;
		bqp.lt_extend_dist=             this.lt_extend_dist;
		bqp.lt_wsigma=                  this.lt_wsigma;
		bqp.lt_max_asigma=              this.lt_max_asigma;
		bqp.lt_max_rsigma=              this.lt_max_rsigma;

		bqp.ml_hwidth=                  this.ml_hwidth;
		bqp.ml_disparity_sweep=         this.ml_disparity_sweep;
		bqp.ml_sweep_steps=             this.ml_sweep_steps;
		bqp.ml_keep_aux=                this.ml_keep_aux;
		bqp.ml_keep_inter=              this.ml_keep_inter;
		bqp.ml_keep_tbrl=               this.ml_keep_tbrl;
		bqp.ml_keep_hor_vert=           this.ml_keep_hor_vert;
		bqp.ml_keep_debug=              this.ml_keep_debug;
		bqp.ml_8bit=                    this.ml_8bit;
		bqp.ml_limit_extrim=            this.ml_limit_extrim;
		bqp.ml_show_ml=                 this.ml_show_ml;
		bqp.ml_fatzero =                this.ml_fatzero;
		return bqp;
	}
}