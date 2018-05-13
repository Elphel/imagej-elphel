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
	public double  inf_max_disp_rig =          0.5;    // maybe even higher (2.0) to lock to initially high mismatch

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
	public double  rig_correction_scale=       1.0;    // scale calcualated correction



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

		this.first_max_disp_main=           gd.getNextNumber();
		this.first_max_disp_aux=            gd.getNextNumber();
		this.first_max_disp_rig=            gd.getNextNumber();
		this.num_refine_master=       (int) gd.getNextNumber();
		this.num_refine_inter=        (int) gd.getNextNumber();
		this.near_max_disp_main=            gd.getNextNumber();
		this.near_max_disp_aux=             gd.getNextNumber();
		this.near_max_disp_rig=             gd.getNextNumber();

		this.refine_min_strength=             gd.getNextNumber();
		this.refine_tolerance=             gd.getNextNumber();

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

		properties.setProperty(prefix+"first_max_disp_main",       this.first_max_disp_main+"");
		properties.setProperty(prefix+"first_max_disp_aux",        this.first_max_disp_aux+"");
		properties.setProperty(prefix+"first_max_disp_rig",        this.first_max_disp_rig+"");
		properties.setProperty(prefix+"num_refine_master",         this.num_refine_master+"");
		properties.setProperty(prefix+"num_refine_inter",          this.num_refine_inter+"");
		properties.setProperty(prefix+"near_max_disp_main",        this.near_max_disp_main+"");
		properties.setProperty(prefix+"near_max_disp_aux",         this.near_max_disp_aux+"");
		properties.setProperty(prefix+"near_max_disp_rig",         this.near_max_disp_rig+"");

		properties.setProperty(prefix+"refine_min_strength",         this.refine_min_strength+"");
		properties.setProperty(prefix+"refine_tolerance",         this.refine_tolerance+"");

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
		if (properties.getProperty(prefix+"rig_adjust_short_threshold")!=null)   this.rig_adjust_short_threshold=Double.parseDouble(properties.getProperty(prefix+"rig_adjust_short_threshold"));
		if (properties.getProperty(prefix+"rig_adjust_orientation")!=null)       this.rig_adjust_orientation=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_orientation"));
		if (properties.getProperty(prefix+"rig_adjust_zoom")!=null)         this.rig_adjust_zoom=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_zoom"));
		if (properties.getProperty(prefix+"rig_adjust_angle")!=null)        this.rig_adjust_angle=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_angle"));
		if (properties.getProperty(prefix+"rig_adjust_distance")!=null)     this.rig_adjust_distance=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_distance"));
		if (properties.getProperty(prefix+"rig_adjust_forward")!=null)      this.rig_adjust_forward=Boolean.parseBoolean(properties.getProperty(prefix+"rig_adjust_forward"));
		if (properties.getProperty(prefix+"rig_correction_scale")!=null)    this.rig_correction_scale=Double.parseDouble(properties.getProperty(prefix+"rig_correction_scale"));

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

		return bqp;


	}
}
