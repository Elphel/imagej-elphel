package com.elphel.imagej.calibration;

import java.awt.Color;
import java.net.MalformedURLException;
import java.util.Properties;

import com.elphel.imagej.cameras.EyesisCameraParameters;
import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.lwir.LwirReaderParameters;

import ij.gui.GenericDialog;

public class CalibrationIllustrationParameters {
	static double          DFLT_LWIR_LO = 22500.0;
	static double          DFLT_LWIR_HI = 23500.0;
	static double          DFLT_EO_R2G =  1.03;  // gain red relative to green
	static double          DFLT_EO_B2G =  1.03;  // gain blue relative to green
	static double          DFLT_EO_HI =   255.0; // range to map to 255.0
	LwirReaderParameters   lwirReaderParameters;
	EyesisCameraParameters eyesisCameraParameters;
	
	double [][]            lwir_ranges; //  = new double [lwirReaderParameters.getLwirChannels(false).length][2];
	double [][]            lwir_offset_gains; // for balancing channels [0] - offset (subtract) [1] - gain (divide) 
	int                    palette = 0; // 0 - white - hot, 1 - black - hot, 2+ - colored
	String                 src_chn_prefix="src_chn-";
	boolean                save_png = true;
	int                    JPEG_quality = 90;
	String                 channel_dir_prefix = "chn_";
	double [][]            eo_rb2g_hi; // []{r2g,b2g,rgb_hi}
	double                 eo_gamma = 0.57;
	double                 eo_minlin_gamma = 20.0;
	double                 eo_saturation = 1.5; // color saturation
	double                 threshold_contrast = 10.0;
	int                    threshold_number =   10;// grid should have this number of nodes with above-threshold contrast
	int                    station_sel =        0xf; // bitmask: 0 - station 0; 1 - station 1, ...
	boolean                station_in_filenames = true;
    Color                  color_grid =         new Color(250, 0,  0, 255);
    Color                  color_grid_weak =    new Color(200, 100,  0, 128);
	Color                  color_grid_extra =   new Color(150, 150,0, 255); // at least one end points to extra (unreliable) nodes (may be null)
	int                    line_width_eo =      3;
	int                    line_width_lwir =    1;

	// kernel illustration
	int                    kernel_dia_direct_lwir =       7; // partial and direct
	int                    kernel_dia_direct_eo =         7;
	int                    kernel_dia_inverse_lwir =      7;
	int                    kernel_dia_inverse_eo =        7;
	int                    kernel_decimate_direct_lwir =  1;
	int                    kernel_decimate_direct_eo =    1;
	int                    kernel_decimate_inverse_lwir = 4;
	int                    kernel_decimate_inverse_eo =   4;
	
	double                 kernel_direct_lwir_min =       0.0;   // -0.007
	double                 kernel_direct_lwir_max =       0.15;  // 0.176
	double                 kernel_direct_red_min =        0.0;   // -0.001
	double                 kernel_direct_red_max =        0.055; // 0.0448
	double                 kernel_direct_blue_min =       0.0;   // -0.0008
	double                 kernel_direct_blue_max =       0.055; // 0.0344
	double                 kernel_direct_green_min =      0.0;   // -0.00118
	double                 kernel_direct_green_max =      0.055; // 0.0627
	
	double                 kernel_inverse_lwir_min =     -0.5;   // -0.73  (-0.44 most) -.144(-.07)
	double                 kernel_inverse_lwir_max =      1.8;   //  0.2   (1.8 - most) .874 (.7) 
	double                 kernel_inverse_red_min =      -0.8;   // -0.73  (-.3)
	double                 kernel_inverse_red_max =       2.75;  //  2.3   (1.7)
	double                 kernel_inverse_blue_min =     -0.8;   // -0.6   (-0.4)
	double                 kernel_inverse_blue_max =      2.75;  //  2.3   (1.94)
	double                 kernel_inverse_green_min =    -0.8;   // -1.1   (-0.53)
	double                 kernel_inverse_green_max =     2.75;  //  3.6   (2.87)
	int                    kernel_lwir_palette =          0;     // 0 - normal, 2 - heat map
	boolean                kernels_normalize =            false; // make each kernel each component maximum be 1.0 
	
	
	boolean                kernel_process_partial =       true;
	boolean                kernel_process_direct =        true;
	boolean                kernel_process_inverse =       true;
	
	boolean                captures_all_lwir =            true;
	boolean                captures_all_eo =              true;
	boolean                captures_all =                 false;
	double                 percentile_min =               0.0; // absolute min            
	double                 percentile_max =               0.0; // absolute max
	double                 max_range =                 1000.0; // 600-1000
	double                 hot_importance =                1.0; // 0 - use value zero, 1.0 - minimal in frame
	boolean                auto_range =                   true; // false use per-channel lo/high
	boolean                auto_lim_range =               true; // if difference min/max exceeds max_range, find "optimal" position of fixed range 
	boolean                captures_annotate =            true; // Imprint timestamp 
    Color                  color_annotate =               new Color(  0,200, 200, 255);
	int                    captures_palette =             0; // 0 - white - hot, 1 - black - hot, 2+ - colored
	long                   selected_channels =            0; // bitmask of selected channels
	double                 min_ts =                       0.0;                   
	double                 max_ts =              2000000000.0;

	int                    auto_range_wnd_type =          3; // 0 - piramid, 1 half-sin, 2-piramid squared, 3 - sin^2
	boolean                calib_offs_gain =              true; // perform offset/gain calibration
	double                 calib_offs_gain_ts =           0.0;  // timestamp to start gain calibration
	double                 calib_offs_gain_dur = 2000000000.0;  // duration for gain calibration 
	double                 noise_sigma =                  5.0;  // image sigma for no-signal
	double                 min_sigma =                  100.0;  // do not use images with smaller sigmas for gain calibration
	double                 autorange_offs_up =            0.2;  // shift range up (fraction of new value)
	double                 autorange_offs_down =          0.05; // shift range down (fraction of new value)
	double                 autorange_range_up =           0.1;  // increase range (fraction of new value)
	double                 autorange_range_down =         0.02; // decrease range (fraction of new value)
	
	
	
	
	public CalibrationIllustrationParameters (
			LwirReaderParameters lwirReaderParameters,
			EyesisCameraParameters eyesisCameraParameters) {
		this.lwirReaderParameters = lwirReaderParameters;
		this.eyesisCameraParameters = eyesisCameraParameters;
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
	
	
	
	public void setProperties(String prefix,Properties properties){
//		properties.setProperty(prefix+"camera_name",         this.camera_name+"");
		set_parameters();
		for (int i = 0; i < lwir_ranges.length; i++) {
			properties.setProperty(prefix+"lwir_range_lo_"+i,         this.lwir_ranges[i][0]+"");
			properties.setProperty(prefix+"lwir_range_hi_"+i,         this.lwir_ranges[i][1]+"");
		}		
		for (int i = 0; i < lwir_offset_gains.length; i++) {
			properties.setProperty(prefix+"lwir_offset_"+i,         this.lwir_offset_gains[i][0]+"");
			properties.setProperty(prefix+"lwir_gain_"+i,           this.lwir_offset_gains[i][1]+"");
		}		
		properties.setProperty(prefix+"palette",            this.palette+"");
		properties.setProperty(prefix+"save_png",           this.save_png+"");
		properties.setProperty(prefix+"JPEG_quality",       this.JPEG_quality+"");
		properties.setProperty(prefix+"channel_dir_prefix", this.channel_dir_prefix);
		for (int i = 0; i < eo_rb2g_hi.length; i++) {
			properties.setProperty(prefix+"eo_r2g_"+i,      this.eo_rb2g_hi[i][0]+"");
			properties.setProperty(prefix+"eo_b2g_"+i,      this.eo_rb2g_hi[i][1]+"");
			properties.setProperty(prefix+"eo_hi_"+i,       this.eo_rb2g_hi[i][2]+"");
		}		
		properties.setProperty(prefix+"eo_gamma",           this.eo_gamma+"");
		properties.setProperty(prefix+"eo_minlin_gamma",    this.eo_minlin_gamma+"");
		properties.setProperty(prefix+"eo_saturation",      this.eo_saturation+"");
		properties.setProperty(prefix+"threshold_contrast", this.threshold_contrast+"");
		properties.setProperty(prefix+"threshold_number",   this.threshold_number+"");
		properties.setProperty(prefix+"station_sel",        this.station_sel+"");
		properties.setProperty(prefix+"station_in_filenames", this.station_in_filenames+"");
		long lcolor_grid =       (color_grid == null) ? -1 : getLongColor(color_grid);
		long lcolor_grid_weak =  (color_grid_weak == null) ? -1 : getLongColor(color_grid_weak);
		long lcolor_grid_extra = (color_grid_extra == null) ? -1 : getLongColor(color_grid_extra);
		properties.setProperty(prefix+"color_grid",         lcolor_grid+"");
		properties.setProperty(prefix+"color_grid_weak",    lcolor_grid_weak+"");
		properties.setProperty(prefix+"color_grid_extra",   lcolor_grid_extra+"");

		properties.setProperty(prefix+"line_width_eo",      this.line_width_eo+"");
		properties.setProperty(prefix+"line_width_lwir",    this.line_width_lwir+"");

		
		properties.setProperty(prefix+"kernel_dia_direct_lwir",      this.kernel_dia_direct_lwir+"");
		properties.setProperty(prefix+"kernel_dia_direct_eo",        this.kernel_dia_direct_eo+"");
		properties.setProperty(prefix+"kernel_dia_inverse_lwir",     this.kernel_dia_inverse_lwir+"");
		properties.setProperty(prefix+"kernel_dia_inverse_eo",       this.kernel_dia_inverse_eo+"");
		properties.setProperty(prefix+"kernel_decimate_direct_lwir", this.kernel_decimate_direct_lwir+"");
		properties.setProperty(prefix+"kernel_decimate_direct_eo",   this.kernel_decimate_direct_eo+"");
		properties.setProperty(prefix+"kernel_decimate_inverse_lwir",this.kernel_decimate_inverse_lwir+"");
		properties.setProperty(prefix+"kernel_decimate_inverse_eo",  this.kernel_decimate_inverse_eo+"");
		properties.setProperty(prefix+"kernel_process_partial",      this.kernel_process_partial+"");
		properties.setProperty(prefix+"kernel_process_direct",       this.kernel_process_direct+"");
		properties.setProperty(prefix+"kernel_process_inverse",      this.kernel_process_inverse+"");

		properties.setProperty(prefix+"kernel_direct_lwir_min",      this.kernel_direct_lwir_min+"");
		properties.setProperty(prefix+"kernel_direct_lwir_max",      this.kernel_direct_lwir_max+"");
		properties.setProperty(prefix+"kernel_direct_red_min",       this.kernel_direct_red_min+"");
		properties.setProperty(prefix+"kernel_direct_red_max",       this.kernel_direct_red_max+"");
		properties.setProperty(prefix+"kernel_direct_blue_min",      this.kernel_direct_blue_min+"");
		properties.setProperty(prefix+"kernel_direct_blue_max",      this.kernel_direct_blue_max+"");
		properties.setProperty(prefix+"kernel_direct_green_min",     this.kernel_direct_green_min+"");
		properties.setProperty(prefix+"kernel_direct_green_max",     this.kernel_direct_green_max+"");

		properties.setProperty(prefix+"kernel_inverse_lwir_min",     this.kernel_inverse_lwir_min+"");
		properties.setProperty(prefix+"kernel_inverse_lwir_max",     this.kernel_inverse_lwir_max+"");
		properties.setProperty(prefix+"kernel_inverse_red_min",      this.kernel_inverse_red_min+"");
		properties.setProperty(prefix+"kernel_inverse_red_max",      this.kernel_inverse_red_max+"");
		properties.setProperty(prefix+"kernel_inverse_blue_min",     this.kernel_inverse_blue_min+"");
		properties.setProperty(prefix+"kernel_inverse_blue_max",     this.kernel_inverse_blue_max+"");
		properties.setProperty(prefix+"kernel_inverse_green_min",    this.kernel_inverse_green_min+"");
		properties.setProperty(prefix+"kernel_inverse_green_max",    this.kernel_inverse_green_max+"");

		properties.setProperty(prefix+"kernel_lwir_palette",         this.kernel_lwir_palette+"");
		properties.setProperty(prefix+"kernels_normalize",           this.kernels_normalize+"");
		
		properties.setProperty(prefix+"captures_all_lwir",           this.captures_all_lwir+"");
		properties.setProperty(prefix+"captures_all_eo",             this.captures_all_eo+"");
		properties.setProperty(prefix+"captures_all",                this.captures_all+"");
		properties.setProperty(prefix+"percentile_min",              this.percentile_min+"");
		properties.setProperty(prefix+"percentile_max",              this.percentile_max+"");
		properties.setProperty(prefix+"max_range",                   this.max_range+"");
		properties.setProperty(prefix+"hot_importance",               this.hot_importance+"");
		properties.setProperty(prefix+"auto_range",                  this.auto_range+"");
		properties.setProperty(prefix+"auto_lim_range",              this.auto_lim_range+"");
	

		properties.setProperty(prefix+"captures_annotate",           this.captures_annotate+"");
		long lcolor_annotate =       (color_annotate == null) ? -1 : getLongColor(color_annotate);
		properties.setProperty(prefix+"color_annotate",              lcolor_annotate+"");
		properties.setProperty(prefix+"captures_palette",            this.captures_palette+"");
		properties.setProperty(prefix+"selected_channels",           this.selected_channels+"");
		properties.setProperty(prefix+"min_ts",                      this.min_ts+"");
		properties.setProperty(prefix+"max_ts",                      this.max_ts+"");

		properties.setProperty(prefix+"auto_range_wnd_type",         this.auto_range_wnd_type+"");
		properties.setProperty(prefix+"calib_offs_gain",             this.calib_offs_gain+"");
		properties.setProperty(prefix+"calib_offs_gain_ts",          this.calib_offs_gain_ts+"");
		properties.setProperty(prefix+"calib_offs_gain_dur",         this.calib_offs_gain_dur+"");
		properties.setProperty(prefix+"noise_sigma",                 this.noise_sigma+"");
		properties.setProperty(prefix+"min_sigma",                   this.min_sigma+"");
		properties.setProperty(prefix+"autorange_offs_up",           this.autorange_offs_up+"");
		properties.setProperty(prefix+"autorange_offs_down",         this.autorange_offs_down+"");
		properties.setProperty(prefix+"autorange_range_up",          this.autorange_range_up+"");
		properties.setProperty(prefix+"autorange_range_down",        this.autorange_range_down+"");
	}
	
	public void getProperties(String prefix,Properties properties){
		set_parameters();
		for (int i = 0; i < lwir_ranges.length; i++) {
			if (properties.getProperty(prefix+"lwir_range_lo_"+i)!=null) {
				this.lwir_ranges[i][0] = Double.parseDouble(properties.getProperty(prefix+"lwir_range_lo_"+i));
			}
			if (properties.getProperty(prefix+"lwir_range_hi_"+i)!=null) {
				this.lwir_ranges[i][1] = Double.parseDouble(properties.getProperty(prefix+"lwir_range_hi_"+i));
			}
		}
		for (int i = 0; i < lwir_offset_gains.length; i++) {
			if (properties.getProperty(prefix+"lwir_offset_"+i)!=null) {
				this.lwir_offset_gains[i][0] = Double.parseDouble(properties.getProperty(prefix+"lwir_offset_"+i));
			}
			if (properties.getProperty(prefix+"lwir_gain_"+i)!=null) {
				this.lwir_offset_gains[i][1] = Double.parseDouble(properties.getProperty(prefix+"lwir_gain_"+i));
			}
		}
		
		/*
		for (int i = 0; i < lwir_offset_gains.length; i++) {
			properties.setProperty(prefix+"lwir_offset_"+i,         this.lwir_offset_gains[i][0]+"");
			properties.setProperty(prefix+"lwir_gain_"+i,           this.lwir_offset_gains[i][1]+"");
		}		
		 */
		
		
		if (properties.getProperty(prefix+"palette")!=null)      this.palette = Integer.parseInt(properties.getProperty(prefix+"palette"));
		if (properties.getProperty(prefix+"save_png")!=null)     this.save_png = Boolean.parseBoolean(properties.getProperty(prefix+"save_png"));
		if (properties.getProperty(prefix+"JPEG_quality")!=null) this.JPEG_quality = Integer.parseInt(properties.getProperty(prefix+"JPEG_quality"));
		if (properties.getProperty(prefix+"channel_dir_prefix")!=null) this.channel_dir_prefix = (String) properties.getProperty(prefix+"channel_dir_prefix");

		for (int i = 0; i < eo_rb2g_hi.length; i++) {
			if (properties.getProperty(prefix+"eo_r2g_"+i)!=null) {
				this.eo_rb2g_hi[i][0] = Double.parseDouble(properties.getProperty(prefix+"eo_r2g_"+i));
			}
			if (properties.getProperty(prefix+"eo_b2g_"+i)!=null) {
				this.eo_rb2g_hi[i][1] = Double.parseDouble(properties.getProperty(prefix+"eo_b2g_"+i));
			}
			if (properties.getProperty(prefix+"eo_hi_"+i)!=null) {
				this.eo_rb2g_hi[i][2] = Double.parseDouble(properties.getProperty(prefix+"eo_hi_"+i));
			}
		}		
		if (properties.getProperty(prefix+"eo_gamma")!=null)             this.eo_gamma = Double.parseDouble(properties.getProperty(prefix+"eo_gamma"));
		if (properties.getProperty(prefix+"eo_minlin_gamma")!=null)      this.eo_minlin_gamma = Double.parseDouble(properties.getProperty(prefix+"eo_minlin_gamma"));
		if (properties.getProperty(prefix+"eo_saturation")!=null)        this.eo_saturation = Double.parseDouble(properties.getProperty(prefix+"eo_saturation"));
		if (properties.getProperty(prefix+"threshold_contrast")!=null)   this.threshold_contrast = Double.parseDouble(properties.getProperty(prefix+"threshold_contrast"));
		if (properties.getProperty(prefix+"threshold_number")!=null)     this.threshold_number = Integer.parseInt(properties.getProperty(prefix+"threshold_number"));
		if (properties.getProperty(prefix+"station_sel")!=null)          this.station_sel = Integer.parseInt(properties.getProperty(prefix+"station_sel"));
		if (properties.getProperty(prefix+"station_in_filenames")!=null) this.station_in_filenames = Boolean.parseBoolean(properties.getProperty(prefix+"station_in_filenames"));
		if  (properties.getProperty(prefix+"color_grid")!=null) {
			long lcolor_grid = Long.parseLong(properties.getProperty(prefix+"color_grid"));
			if (lcolor_grid < 0) this.color_grid = null;
			else this.color_grid = setLongColor(lcolor_grid);
		}
		if  (properties.getProperty(prefix+"color_grid_weak")!=null) {
			long lcolor_grid_weak = Long.parseLong(properties.getProperty(prefix+"color_grid_weak"));
			if (lcolor_grid_weak < 0) this.color_grid_weak = null;
			else this.color_grid_weak = setLongColor(lcolor_grid_weak);
		}
		if  (properties.getProperty(prefix+"color_grid_extra")!=null) {
			long lcolor_grid_extra = Long.parseLong(properties.getProperty(prefix+"color_grid_extra"));
			if (lcolor_grid_extra < 0) this.color_grid_extra = null;
			else this.color_grid_extra = setLongColor(lcolor_grid_extra);
		}
		if (properties.getProperty(prefix+"line_width_eo")!=null)      this.line_width_eo = Integer.parseInt(properties.getProperty(prefix+"line_width_eo"));
		if (properties.getProperty(prefix+"line_width_lwir")!=null)    this.line_width_lwir = Integer.parseInt(properties.getProperty(prefix+"line_width_lwir"));
		
		if (properties.getProperty(prefix+"kernel_dia_direct_lwir")!=null)       this.kernel_dia_direct_lwir = Integer.parseInt(properties.getProperty(prefix+"kernel_dia_direct_lwir"));
		if (properties.getProperty(prefix+"kernel_dia_direct_eo")!=null)         this.kernel_dia_direct_eo = Integer.parseInt(properties.getProperty(prefix+"kernel_dia_direct_eo"));
		if (properties.getProperty(prefix+"kernel_dia_inverse_lwir")!=null)      this.kernel_dia_inverse_lwir = Integer.parseInt(properties.getProperty(prefix+"kernel_dia_inverse_lwir"));
		if (properties.getProperty(prefix+"kernel_dia_inverse_eo")!=null)        this.kernel_dia_inverse_eo = Integer.parseInt(properties.getProperty(prefix+"kernel_dia_inverse_eo"));
		if (properties.getProperty(prefix+"kernel_decimate_direct_lwir")!=null)  this.kernel_decimate_direct_lwir = Integer.parseInt(properties.getProperty(prefix+"kernel_decimate_direct_lwir"));
		if (properties.getProperty(prefix+"kernel_decimate_direct_eo")!=null)    this.kernel_decimate_direct_eo = Integer.parseInt(properties.getProperty(prefix+"kernel_decimate_direct_eo"));
		if (properties.getProperty(prefix+"kernel_decimate_inverse_lwir")!=null) this.kernel_decimate_inverse_lwir = Integer.parseInt(properties.getProperty(prefix+"kernel_decimate_inverse_lwir"));
		if (properties.getProperty(prefix+"kernel_decimate_inverse_eo")!=null)   this.kernel_decimate_inverse_eo = Integer.parseInt(properties.getProperty(prefix+"kernel_decimate_inverse_eo"));
		if (properties.getProperty(prefix+"kernel_process_partial")!=null)       this.kernel_process_partial = Boolean.parseBoolean(properties.getProperty(prefix+"kernel_process_partial"));
		if (properties.getProperty(prefix+"kernel_process_direct")!=null)        this.kernel_process_direct = Boolean.parseBoolean(properties.getProperty(prefix+"kernel_process_direct"));
		if (properties.getProperty(prefix+"kernel_process_inverse")!=null)       this.kernel_process_inverse = Boolean.parseBoolean(properties.getProperty(prefix+"kernel_process_inverse"));
		
		if (properties.getProperty(prefix+"kernel_direct_lwir_min")!=null)       this.kernel_direct_lwir_min = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_lwir_min"));
		if (properties.getProperty(prefix+"kernel_direct_lwir_max")!=null)       this.kernel_direct_lwir_max = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_lwir_max"));
		if (properties.getProperty(prefix+"kernel_direct_red_min")!=null)        this.kernel_direct_red_min = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_red_min"));
		if (properties.getProperty(prefix+"kernel_direct_red_max")!=null)        this.kernel_direct_red_max = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_red_max"));
		if (properties.getProperty(prefix+"kernel_direct_blue_min")!=null)       this.kernel_direct_blue_min = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_blue_min"));
		if (properties.getProperty(prefix+"kernel_direct_blue_max")!=null)       this.kernel_direct_blue_max = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_blue_max"));
		if (properties.getProperty(prefix+"kernel_direct_green_min")!=null)      this.kernel_direct_green_min = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_green_min"));
		if (properties.getProperty(prefix+"kernel_direct_green_max")!=null)      this.kernel_direct_green_max = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_green_max"));

		if (properties.getProperty(prefix+"kernel_inverse_lwir_min")!=null)      this.kernel_inverse_lwir_min = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_lwir_min"));
		if (properties.getProperty(prefix+"kernel_inverse_lwir_max")!=null)      this.kernel_inverse_lwir_max = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_lwir_max"));
		if (properties.getProperty(prefix+"kernel_inverse_red_min")!=null)       this.kernel_inverse_red_min = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_red_min"));
		if (properties.getProperty(prefix+"kernel_inverse_red_max")!=null)       this.kernel_inverse_red_max = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_red_max"));
		if (properties.getProperty(prefix+"kernel_inverse_blue_min")!=null)      this.kernel_inverse_blue_min = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_blue_min"));
		if (properties.getProperty(prefix+"kernel_inverse_blue_max")!=null)      this.kernel_inverse_blue_max = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_blue_max"));
		if (properties.getProperty(prefix+"kernel_inverse_green_max")!=null)     this.kernel_inverse_green_max = Double.parseDouble(properties.getProperty(prefix+"kernel_inverse_green_max"));
		if (properties.getProperty(prefix+"kernel_direct_lwir_min")!=null)       this.kernel_direct_lwir_min = Double.parseDouble(properties.getProperty(prefix+"kernel_direct_lwir_min"));
		if (properties.getProperty(prefix+"kernel_lwir_palette")!=null)          this.kernel_lwir_palette = Integer.parseInt(properties.getProperty(prefix+"kernel_lwir_palette"));
		if (properties.getProperty(prefix+"kernels_normalize")!=null)            this.kernels_normalize = Boolean.parseBoolean(properties.getProperty(prefix+"kernels_normalize"));
		if (properties.getProperty(prefix+"captures_all_lwir")!=null)            this.captures_all_lwir = Boolean.parseBoolean(properties.getProperty(prefix+"captures_all_lwir"));
		if (properties.getProperty(prefix+"captures_all_eo")!=null)              this.captures_all_eo = Boolean.parseBoolean(properties.getProperty(prefix+"captures_all_eo"));
		if (properties.getProperty(prefix+"captures_all")!=null)                 this.captures_all = Boolean.parseBoolean(properties.getProperty(prefix+"captures_all"));
		if (properties.getProperty(prefix+"percentile_min")!=null)               this.percentile_min = Double.parseDouble(properties.getProperty(prefix+"percentile_min"));
		if (properties.getProperty(prefix+"percentile_max")!=null)               this.percentile_max = Double.parseDouble(properties.getProperty(prefix+"percentile_max"));
		if (properties.getProperty(prefix+"max_range")!=null)                    this.max_range = Double.parseDouble(properties.getProperty(prefix+"max_range"));
		if (properties.getProperty(prefix+"hot_importance")!=null)               this.hot_importance = Double.parseDouble(properties.getProperty(prefix+"hot_importance"));
		if (properties.getProperty(prefix+"auto_range")!=null)                   this.auto_range = Boolean.parseBoolean(properties.getProperty(prefix+"auto_range"));
		if (properties.getProperty(prefix+"auto_lim_range")!=null)               this.auto_lim_range = Boolean.parseBoolean(properties.getProperty(prefix+"auto_lim_range"));
		if (properties.getProperty(prefix+"captures_annotate")!=null)            this.captures_annotate = Boolean.parseBoolean(properties.getProperty(prefix+"captures_annotate"));
		if  (properties.getProperty(prefix+"color_annotate")!=null) {
			long lcolor_annotate = Long.parseLong(properties.getProperty(prefix+"color_annotate"));
			if (lcolor_annotate < 0) this.color_annotate = null;
			else this.color_annotate = setLongColor(lcolor_annotate);
		}
		if (properties.getProperty(prefix+"captures_palette")!=null)             this.captures_palette = Integer.parseInt(properties.getProperty(prefix+"captures_palette"));
		if (properties.getProperty(prefix+"selected_channels")!=null)            this.selected_channels = Long.parseLong(properties.getProperty(prefix+"selected_channels"));
		if (properties.getProperty(prefix+"min_ts")!=null)                       this.min_ts = Double.parseDouble(properties.getProperty(prefix+"min_ts"));
		if (properties.getProperty(prefix+"max_ts")!=null)                       this.max_ts = Double.parseDouble(properties.getProperty(prefix+"max_ts"));

		if (properties.getProperty(prefix+"auto_range_wnd_type")!=null)          this.auto_range_wnd_type = Integer.parseInt(properties.getProperty(prefix+"auto_range_wnd_type"));
		if (properties.getProperty(prefix+"calib_offs_gain")!=null)              this.calib_offs_gain = Boolean.parseBoolean(properties.getProperty(prefix+"calib_offs_gain"));
		if (properties.getProperty(prefix+"calib_offs_gain_ts")!=null)           this.calib_offs_gain_ts = Double.parseDouble(properties.getProperty(prefix+"calib_offs_gain_ts"));
		if (properties.getProperty(prefix+"calib_offs_gain_dur")!=null)          this.calib_offs_gain_dur = Double.parseDouble(properties.getProperty(prefix+"calib_offs_gain_dur"));
		if (properties.getProperty(prefix+"noise_sigma")!=null)                  this.noise_sigma = Double.parseDouble(properties.getProperty(prefix+"noise_sigma"));
		if (properties.getProperty(prefix+"min_sigma")!=null)                    this.min_sigma = Double.parseDouble(properties.getProperty(prefix+"min_sigma"));
		if (properties.getProperty(prefix+"autorange_offs_up")!=null)            this.autorange_offs_up = Double.parseDouble(properties.getProperty(prefix+"autorange_offs_up"));
		if (properties.getProperty(prefix+"autorange_offs_down")!=null)          this.autorange_offs_down = Double.parseDouble(properties.getProperty(prefix+"autorange_offs_down"));
		if (properties.getProperty(prefix+"autorange_range_up")!=null)           this.autorange_range_up = Double.parseDouble(properties.getProperty(prefix+"autorange_range_up"));
		if (properties.getProperty(prefix+"autorange_range_down")!=null)         this.autorange_range_down = Double.parseDouble(properties.getProperty(prefix+"autorange_range_down"));
	}
	
	public void dialogQuestions(GenericJTabbedDialog gd) {
		final int lwir0 = getLwirReaderParameters().getLwirChn0();
		final int eo0 =   getLwirReaderParameters().getEoChn0();
		final int numStations = eyesisCameraParameters.getNumStations();
        gd.addTab("LWIR","LWIR photometric parameters");
		gd.addNumericField("Thermal color palette", this.palette,  0,3,"","0 - white-hot, 1 - black-hot, 2+ - colored");
		for (int i = 0; i < lwir_ranges.length; i++) {
			gd.addNumericField("LWIR chn:"+(i+lwir0)+" center", 0.5*(this.lwir_ranges[i][1]+this.lwir_ranges[i][0]),  1,8,"","LWIR sensor range center");
			gd.addNumericField("LWIR chn:"+(i+lwir0)+" range", this.lwir_ranges[i][1]-this.lwir_ranges[i][0],  1,8,"","LWIR sensor range");
		}
        gd.addTab("EO","EO photometric parameters");
		for (int i = 0; i < eo_rb2g_hi.length; i++) {
			gd.addNumericField("EO chn:"+(i+eo0)+" r2g", this.eo_rb2g_hi[i][0],  4,8,"","gain ratio red to green");
			gd.addNumericField("EO chn:"+(i+eo0)+" b2g", this.eo_rb2g_hi[i][1],  4,8,"","gain ratio red to green");
			gd.addNumericField("EO chn:"+(i+eo0)+" full range", this.eo_rb2g_hi[i][2],  1,8,"","Full intensity range to map to 255 in the output");
		}
		gd.addNumericField("Gamma correction",  this.eo_gamma,        3,6,"","Gamma correction: 1.0 - linear, lower increases dynamic range.");
		gd.addNumericField("Linear range",      this.eo_minlin_gamma, 1,6,"","Apply gamma only to higher input values");
		gd.addNumericField("Color saturation",  this.eo_saturation,   3,6,"","Boost color saturation from the linear representation");
		
        gd.addTab("Grid lines","Grid lines parameters");
		gd.addNumericField("Line width for LWIR images",   this.line_width_lwir,  0,3,"pix","line width for the grid in LWIR images");
		gd.addNumericField("Line width for LWIR images",   this.line_width_eo,  0,3,"pix","line width for the grid in EO images");
		String scolor = (this.color_grid==null)?"none":String.format("%08x", getLongColor(this.color_grid));
		gd.addStringField ("Line color for the grid lines (actual grid)",scolor, 8, "Any invalid hex number disables drawing lines");
		String scolor_weak = (this.color_grid_weak==null)?"none":String.format("%08x", getLongColor(this.color_grid_weak));
		gd.addStringField ("Line color for the weak grid lines (below threshold)",scolor_weak, 8, "Any invalid hex number disables drawing lines");
		String scolor_extra = (this.color_grid_extra==null)?"none":String.format("%08x", getLongColor(this.color_grid_extra));
		gd.addStringField ("Line color for the extended grid lines (extra grid)",scolor_extra, 8, "Any invalid hex number disables drawing lines");
        
        gd.addTab("General","Output format and other parameters");
    	gd.addStringField ("Channel directory prefix",this.channel_dir_prefix, 15,"Prefix to a directory name to save channel annotated files");
        gd.addMessage("Stations to include");
        for (int i = 0; i < numStations;i++) {
        	gd.addCheckbox("Use station "+i, ((this.station_sel >> i) & 1) != 0, "Include Station "+i+" in generated files.");
        }
    	gd.addCheckbox("Include station number in result file names", this.station_in_filenames, "Use station number as a part of the result file names.");
		gd.addCheckbox("Save as PNG instead of JPEG", save_png);
		gd.addNumericField("JPEG quality", this.JPEG_quality,  0,3,"","JPEG quality, 0 - use Tiff");
		gd.addNumericField("Threshold contrast",  this.threshold_contrast,   3,6,"","Consider grid nodes with higher contrast to determine bad grids.");
		gd.addNumericField("Minimal number of high-contrast nodes", this.threshold_number,  0,3,"","Consider a failed grid if the number of strong nodes is below this.");

		gd.addTab("Kernels","'Humanize' partial, direct and inverse kernels");
		gd.addNumericField("Kernel (partial, direct) size to keep, LWIR", this.kernel_dia_direct_lwir,  0,3,"pix","Size=3 preserves 3x3 center from each kernel, even add 1 pixel gap");
		gd.addNumericField("Kernel (partial, direct) size to keep, EO",   this.kernel_dia_direct_eo,  0,3,"pix","Size=3 preserves 3x3 center from each kernel, even add 1 pixel gap");
		gd.addNumericField("Kernel (inverted) size to keep, LWIR",        this.kernel_dia_inverse_lwir,  0,3,"pix","Size=3 preserves 3x3 center from each kernel, even add 1 pixel gap");
		gd.addNumericField("Kernel (inverted) size to keep, EO",          this.kernel_dia_inverse_eo,  0,3,"pix","Size=3 preserves 3x3 center from each kernel, even add 1 pixel gap");
		gd.addNumericField("Decimate direct LWIR kernels",                this.kernel_decimate_direct_lwir,  0,3,"","reduce output kernel density");
		gd.addNumericField("Decimate direct EO kernels",                  this.kernel_decimate_direct_eo,  0,3,"pix","reduce output kernel density");
		gd.addNumericField("Decimate inverted LWIR kernels",              this.kernel_decimate_inverse_lwir,  0,3,"pix","reduce output kernel density");
		gd.addNumericField("Decimate inverted EO kernels",                this.kernel_decimate_inverse_eo,  0,3,"pix","reduce output kernel density");
    	gd.addCheckbox    ("Illustrate partial kernels",                  this.kernel_process_partial, "Process and save partial kernels");
    	gd.addCheckbox    ("Illustrate direct kernels",                   this.kernel_process_direct,  "Process and save combined direct kernels");
    	gd.addCheckbox    ("Illustrate inverted kernels",                 this.kernel_process_inverse, "Process and save inverted kernels");
    	
    	gd.addMessage("Levels for direct LWIR kernels");
		gd.addNumericField("Low level for LWIR direct kernels",           this.kernel_direct_lwir_min,   4,8,"","Applies to partial and direct kernels (from histograms: -0.007)");
		gd.addNumericField("High level for LWIR direct kernels",          this.kernel_direct_lwir_max,   4,8,"","Applies to partial and direct kernels (from histograms:  0.176)");

    	gd.addMessage("Levels for direct RGB kernels (same for all color components to merge into a color image)");
		gd.addNumericField("Low level for RGB/red direct kernels",        this.kernel_direct_red_min,    4,8,"","Applies to partial and direct kernels (from histograms: -0.001)");
		gd.addNumericField("High level for RGB/red direct kernels",       this.kernel_direct_red_max,    4,8,"","Applies to partial and direct kernels (from histograms:  0.0448)");
		gd.addNumericField("Low level for RGB/green direct kernels",      this.kernel_direct_green_min,  4,8,"","Applies to partial and direct kernels (from histograms: -0.00118)");
		gd.addNumericField("High level for RGB/green direct kernels",     this.kernel_direct_green_max,  4,8,"","Applies to partial and direct kernels (from histograms:  0.0627)");
		gd.addNumericField("Low level for RGB/blue direct kernels",       this.kernel_direct_blue_min,   4,8,"","Applies to partial and direct kernels (from histograms: -0.0008)");
		gd.addNumericField("High level for RGB/blue direct kernels",      this.kernel_direct_blue_max,   4,8,"","Applies to partial and direct kernels (from histograms:  0.0344)");
    	
    	gd.addMessage("Levels for inverse LWIR kernels");
		gd.addNumericField("Low level for LWIR inverse kernels",           this.kernel_inverse_lwir_min,   4,8,"","Applies to inverse kernels (from histograms: -0.144  (-0.07 most))");
		gd.addNumericField("High level for LWIR inverse kernels",          this.kernel_inverse_lwir_max,   4,8,"","Applies to inverse kernels (from histograms:  .874   (.7 - most))");

    	gd.addMessage("Levels for inverse RGB kernels");
		gd.addNumericField("Low level for RGB/red inverse kernels",        this.kernel_inverse_red_min,  4,8,"","Applies to inverse kernels (from histograms: -0.73  (-0.3 most))");
		gd.addNumericField("High level for RGB/red inverse kernels",       this.kernel_inverse_red_max,  4,8,"","Applies to inverse kernels (from histograms:  0.23  ( 1.7 most))");
		gd.addNumericField("Low level for RGB/green inverse kernels",      this.kernel_inverse_green_min,4,8,"","Applies to inverse kernels (from histograms: -1.1   (-0.53 most))");
		gd.addNumericField("High level for RGB/green inverse kernels",     this.kernel_inverse_green_max,4,8,"","Applies to inverse kernels (from histograms:  3.6   ( 2.87 most))");
		gd.addNumericField("Low level for RGB/blue inverse kernels",       this.kernel_inverse_blue_min, 4,8,"","Applies to inverse kernels (from histograms: -0.6   (-0.4 most))");
		gd.addNumericField("High level for RGB/blue inverse kernels",      this.kernel_inverse_blue_max, 4,8,"","Applies to inverse kernels (from histograms:  2.3   (1.94 most))");

		gd.addNumericField("Palette to use for monochrome kernels",        this.kernel_lwir_palette,     0,3,"","0 - normal, 1 - inverted, 2+ - heatmap");
    	gd.addCheckbox    ("Normalize each kernel/component",              this.kernels_normalize, "Make each kernel component maximum be 1.0");
    	
        gd.addTab("Field footage","Coverting raw captured images into PNG/JPEG sequences to convert to video");
    	gd.addCheckbox    ("Ignore scenes with partial LWIR channels",    this.captures_all_lwir, "Skip scenes where not all or none LWIR channels are available");
    	gd.addCheckbox    ("Ignore scenes with partial EO channels",      this.captures_all_eo,   "Skip scenes where not all or none EO channels are available");
    	gd.addCheckbox    ("Ignore partial scenes",                       this.captures_all,      "Skip scenes where not all (20) images are present");
		gd.addNumericField("Ignore coldest percentile",               100*this.percentile_min,    4,8,"%","This fraction of all pixels will be considered 'too cold'");
		gd.addNumericField("Ignore hottest percentile",               100*this.percentile_max,    4,8,"%","This fraction of all pixels will be considered 'too hot'");
		gd.addNumericField("Maximal range",                               this.max_range,    4,8,"counts","Maximal range (in 16-bit raw TIFF counts) to use while presenting output data");
		gd.addNumericField("Hot importance",                              this.hot_importance,    4,8,"","Maximal range (in 16-bit raw TIFF counts) to use while presenting output data");
    	gd.addCheckbox    ("Auto range LWIR images",                      this.auto_range, "Auto range images, if unchecked - use per-channel ranges specified in LWIR tab");
    	gd.addCheckbox    ("Limit range",                                 this.auto_lim_range, "When difference between max and min exceeds maximal range, use maximal range and find best offset");
    	gd.addCheckbox    ("Annotate images with timestamps",             this.captures_annotate, "Imprint timestamps into output images");
		scolor = (this.color_annotate==null)?"none":String.format("%08x", getLongColor(this.color_annotate));
		gd.addStringField ("Line color for the grid lines (actual grid)",scolor, 8, "Any invalid hex number disables annotation");
		gd.addNumericField("Thermal color palette", this.captures_palette,  0,3,"","0 - white-hot, 1 - black-hot, 2+ - colored");
		gd.addNumericField("Minimal timestamp to process",                this.min_ts, 4,20,"s","Do not process scenes earlier than that timestamp");
		gd.addNumericField("Maximal timestamp to process",                this.max_ts, 4,20,"s","Do not process scenes later than that timestamp");
    	
		gd.addNumericField("Sensor balancing window type",                this.auto_range_wnd_type,     0,3,"","0 - piramid, 1 half-sin, 2-piramid squared, 3 - sin^2");
    	gd.addCheckbox    ("Perform LWIR sensors offset/gain balancing ", this.calib_offs_gain, "Calculate sensor channels offset and gains from selected imagery");
		gd.addNumericField("Balancing start timestamp",                   this.calib_offs_gain_ts, 4,20,"s","Earliest scene to use for sensor channel balancing");
		gd.addNumericField("Balancing duration",                          this.calib_offs_gain_dur, 4,20,"s","Balancing duration (end timestamp - start timestamp)");
		gd.addNumericField("Noise sigma",                                 this.noise_sigma, 2,8,"counts","Typical sigma of the LWIR images with uniform illumination");
		gd.addNumericField("Minimal LWIR sigma",                          this.min_sigma, 2,8,"counts","Do not use images were sigma is below for gain calibration");
		gd.addNumericField("Shift range up (fraction of new value)",      this.autorange_offs_up, 4,6,"","If system wants to shift range up, use this (<=1.0) fraction of the new value");
		gd.addNumericField("Shift range down (fraction of new value)",    this.autorange_offs_down, 4,6,"","If system wants to shift range up, use this (<=1.0) fraction of the new value");
		gd.addNumericField("Increase range (fraction of new value)",      this.autorange_range_up, 4,6,"","If system wants to increase range, use this (<=1.0) fraction of the new value");
		gd.addNumericField("Decrease range (fraction of new value)",      this.autorange_range_down, 4,6,"","If system wants to decrease range, use this (<=1.0) fraction of the new value");

		gd.addCheckbox    ("Edit selected channels",                      false, "check to edit selected channels");
	}
	
	//kernels_normalize
	public void dialogAnswers(GenericJTabbedDialog gd) {
// --- LWIR ---		
		this.palette =          (int) gd.getNextNumber();
		for (int i = 0; i < lwir_ranges.length; i++) {
			double center = gd.getNextNumber();
			double range =  gd.getNextNumber();
			this.lwir_ranges[i][0] =  center - 0.5*range;
			this.lwir_ranges[i][1] =  center + 0.5*range;
		}
// --- EO ---		
		for (int i = 0; i < eo_rb2g_hi.length; i++) {
			this.eo_rb2g_hi[i][0] =  gd.getNextNumber();
			this.eo_rb2g_hi[i][1] =  gd.getNextNumber();
			this.eo_rb2g_hi[i][2] =  gd.getNextNumber();
		}
		this.eo_gamma =              gd.getNextNumber();
		this.eo_minlin_gamma =       gd.getNextNumber();
		this.eo_saturation =         gd.getNextNumber();		
// --- Grid lines ---
		this.line_width_lwir = (int) gd.getNextNumber();
		this.line_width_eo =   (int) gd.getNextNumber();
		String scolor =              gd.getNextString();
		long lcolor = -1;
		try {
			lcolor = Long.parseLong(scolor,16);
			this.color_grid = setLongColor(lcolor);
		} catch(NumberFormatException e){
			this.color_grid = null;
		}
		scolor =                     gd.getNextString();
		try {
			lcolor = Long.parseLong(scolor,16);
			this.color_grid_weak = setLongColor(lcolor);
		} catch(NumberFormatException e){
			this.color_grid_weak = null;
		}	

		scolor =                     gd.getNextString();
		try {
			lcolor = Long.parseLong(scolor,16);
			this.color_grid_extra = setLongColor(lcolor);
		} catch(NumberFormatException e){
			this.color_grid_extra = null;
		}	
		
// --- General ---		
		this.channel_dir_prefix =    gd.getNextString();
		this.station_sel = 0;
		for (int i = 0; i < eyesisCameraParameters.getNumStations();i++) {
			if (gd.getNextBoolean()) {
				this.station_sel |=  1 << i;
			}
		}
		this.station_in_filenames =               gd.getNextBoolean();
		this.save_png =                           gd.getNextBoolean();
		this.JPEG_quality =                 (int) gd.getNextNumber();
		this.threshold_contrast =                 gd.getNextNumber();
		this.threshold_number=              (int) gd.getNextNumber();
// --- Kernels ---
		this.kernel_dia_direct_lwir =       (int) gd.getNextNumber();
		this.kernel_dia_direct_eo =         (int) gd.getNextNumber();
		this.kernel_dia_inverse_lwir =      (int) gd.getNextNumber();
		this.kernel_dia_inverse_eo =        (int) gd.getNextNumber();
		this.kernel_decimate_direct_lwir =  (int) gd.getNextNumber();
		this.kernel_decimate_direct_eo =    (int) gd.getNextNumber();
		this.kernel_decimate_inverse_lwir = (int) gd.getNextNumber();
		this.kernel_decimate_inverse_eo =   (int) gd.getNextNumber();
		this.kernel_process_partial =             gd.getNextBoolean();
		this.kernel_process_direct =              gd.getNextBoolean();
		this.kernel_process_inverse =             gd.getNextBoolean();
// kernel levels
		this.kernel_direct_lwir_min =             gd.getNextNumber();
		this.kernel_direct_lwir_max =             gd.getNextNumber();
		this.kernel_direct_red_min =              gd.getNextNumber();
		this.kernel_direct_red_max =              gd.getNextNumber();
		this.kernel_direct_green_min =            gd.getNextNumber();
		this.kernel_direct_green_max =            gd.getNextNumber();
		this.kernel_direct_blue_min =             gd.getNextNumber();
		this.kernel_direct_blue_max =             gd.getNextNumber();
		this.kernel_inverse_lwir_min =            gd.getNextNumber();
		this.kernel_inverse_lwir_max =            gd.getNextNumber();
		this.kernel_inverse_red_min =             gd.getNextNumber();
		this.kernel_inverse_red_max =             gd.getNextNumber();
		this.kernel_inverse_green_min =           gd.getNextNumber();
		this.kernel_inverse_green_max =           gd.getNextNumber();
		this.kernel_inverse_blue_min =            gd.getNextNumber();
		this.kernel_inverse_blue_max =            gd.getNextNumber();
		this.kernel_lwir_palette =          (int) gd.getNextNumber();
		this.kernels_normalize =                  gd.getNextBoolean();

		this.captures_all_lwir =                  gd.getNextBoolean();
		this.captures_all_eo =                    gd.getNextBoolean();
		this.captures_all =                       gd.getNextBoolean();
		this.percentile_min =                0.01*gd.getNextNumber();
		this.percentile_max =                0.01*gd.getNextNumber();
		this.max_range =                          gd.getNextNumber();
		this.hot_importance =                     gd.getNextNumber();
		this.auto_range =                         gd.getNextBoolean();
		this.auto_lim_range =                     gd.getNextBoolean();
		this.captures_annotate =                  gd.getNextBoolean();
		scolor =                                  gd.getNextString();
		lcolor = -1;
		try {
			lcolor = Long.parseLong(scolor,16);
			this.color_annotate = setLongColor(lcolor);
		} catch(NumberFormatException e){
			this.color_annotate = null;
		}
		this.captures_palette =             (int) gd.getNextNumber();
		this.min_ts =                             gd.getNextNumber();
		this.max_ts =                             gd.getNextNumber();
		
		this.auto_range_wnd_type =          (int) gd.getNextNumber();
		this.calib_offs_gain =                    gd.getNextBoolean();
		this.calib_offs_gain_ts =                 gd.getNextNumber();
		this.calib_offs_gain_dur =                gd.getNextNumber();
		this.noise_sigma =                        gd.getNextNumber();
		this.min_sigma =                          gd.getNextNumber();
		this.autorange_offs_up =                  gd.getNextNumber();
		this.autorange_offs_down =                gd.getNextNumber();
		this.autorange_range_up =                 gd.getNextNumber();
		this.autorange_range_down =               gd.getNextNumber();
		
		if (gd.getNextBoolean()) {
			selectChannelsToProcess("Select channels to process");
		}
	}
	
	public boolean showJDialog() {
		set_parameters();
		GenericJTabbedDialog gd = new GenericJTabbedDialog("Set illustration parameters",800,900);
		dialogQuestions(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		dialogAnswers(gd);
		return true;
	}
	public void set_parameters () {
		//		this.lwirReaderParameters = lwirReaderParameters;
		if ((lwir_ranges == null) || (lwir_ranges.length != lwirReaderParameters.getLwirChannels(false).length)){
			lwir_ranges = new double [lwirReaderParameters.getLwirChannels(false).length][2];
			for (int i = 0; i < lwir_ranges.length; i++) {
				this.lwir_ranges[i][0] = DFLT_LWIR_LO;
				this.lwir_ranges[i][1] = DFLT_LWIR_HI;
			}
		}
		if ((eo_rb2g_hi == null) || (eo_rb2g_hi.length != lwirReaderParameters.getEoChannels(false).length)) {
			eo_rb2g_hi = new double [lwirReaderParameters.getEoChannels(false).length][3];
			for (int i = 0; i < lwir_ranges.length; i++) {
				eo_rb2g_hi[i][0] = DFLT_EO_R2G;  // gain red relative to green
				eo_rb2g_hi[i][1] = DFLT_EO_B2G;  // gain blue relative to green
				eo_rb2g_hi[i][2] = DFLT_EO_HI;   // range to map to 255.0
			}
		}
		if ((lwir_offset_gains == null) || (lwir_offset_gains.length != lwirReaderParameters.getLwirChannels(false).length)){
			lwir_offset_gains = new double [lwirReaderParameters.getLwirChannels(false).length][2];
			for (int i = 0; i < lwir_offset_gains.length; i++) {
				this.lwir_offset_gains[i][0] = 0.0;
				this.lwir_offset_gains[i][1] = 1.0;
			}
		}
	}
	
	public void setLWIROffsetGains(double [][] offset_gains) {
		for (int i = 0; i < lwir_offset_gains.length; i++) {  // offset_gains may be [20], while lwir_offset_gains - just 16
			if (offset_gains[i] != null) {
				lwir_offset_gains[i] = offset_gains[i].clone(); 
			}
		}
	}

	public double [][] getLWIROffsetGains(){
		return lwir_offset_gains;
	}
	
	public double getGamma() {
		return this.eo_gamma;
	}
	
	public double getMinLin() {
		return this.eo_minlin_gamma;
	}
	public double getSaturation() {
		return this.eo_saturation;
	}

	public int getLineWidthLwir() {
		return this.line_width_lwir;
	}

	public int getLineWidthEo() {
		return this.line_width_eo;
	}
	public Color getGridColor() {
		return this.color_grid;
	}
	public Color getGridWeakColor() {
		return this.color_grid_weak;
	}
	public Color getGridExtraColor() {
		return this.color_grid_extra;
	}
	public String getChannelPrefix() {
		return this.channel_dir_prefix;
	}
	public boolean useStation(int st_num) {
		return (this.station_sel & (1 << st_num)) != 0;
	}
	public boolean useStationInFilenames() {
		return this.station_in_filenames;
	}
	public boolean usePNG() {
		return this.save_png;
	}
	public int getJPEG_quality() {
		return this.JPEG_quality;
	}
	
	public double getThresholdContrast() {
		return this.threshold_contrast;
	}

	public int getThresholdNumber() {
		return this.threshold_number;
	}

	public int getPalette() {
		return this.palette;
	}
	public double [] getLwirRange(int lwir_index) {
		return lwir_ranges[lwir_index];
	}
	public LwirReaderParameters getLwirReaderParameters() {
		return lwirReaderParameters;
	}
	public boolean [] getSelectedChannels() {
		int num_lwir_channels = lwirReaderParameters.getLwirChannels(true).length;
		int num_eo_channels =   lwirReaderParameters.getEoChannels(true).length;
		boolean [] selection = new boolean[num_lwir_channels + num_eo_channels];
		for (int i = 0; i < selection.length; i++) {
			selection[i] = ((this.selected_channels >> i) & 1) != 0;
		}
		return selection;
		
	}
	public boolean selectChannelsToProcess(String title) {
		int num_lwir_channels = lwirReaderParameters.getLwirChannels(true).length;
		int num_eo_channels =   lwirReaderParameters.getEoChannels(true).length;
		boolean [] selection = getSelectedChannels();
		int numChannels=selection.length;
		while (true) {
			GenericJTabbedDialog gd = new GenericJTabbedDialog(title);
//			GenericDialog gd = new GenericDialog(title);
			gd.addMessage("LWIR channels");
			for (int i=0;i<num_lwir_channels;i++) gd.addCheckbox("channel "+i, selection[i]);
			gd.addMessage("EO channels");
			for (int i=num_lwir_channels; i < (num_lwir_channels + num_eo_channels);i++) gd.addCheckbox("channel "+i, selection[i]);
			String [] labels =   {"OK", "All like channel 0/"+num_lwir_channels, "Cancel"};
			String [] actions = {"OK", "_All like channel 0/"+num_lwir_channels, "Cancel"};
			String [] tooltips = labels;
			gd.addButtons(labels, actions, tooltips);
//			gd.enableYesNoCancel("OK", "All like channel 0");
//			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			for (int i=0;i<numChannels;i++) selection[i]=gd.getNextBoolean();
			if (gd.wasOKed()){
				this.selected_channels = 0;
				for (int i=0;i<numChannels;i++) {
					this.selected_channels |= (selection[i] ? 1 : 0) << i;
				}
				return true;
			} else {
				for (int i=1;i<num_lwir_channels;i++) selection[i]=selection[0];
				
				for (int i=1;i<num_eo_channels;i++) selection[num_lwir_channels+i]=selection[num_lwir_channels];
			}
		}
	}
}
