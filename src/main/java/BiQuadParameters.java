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
		return bqp;
	}
}
