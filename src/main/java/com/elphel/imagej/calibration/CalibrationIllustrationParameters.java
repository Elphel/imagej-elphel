package com.elphel.imagej.calibration;

import java.awt.Color;
import java.net.MalformedURLException;
import java.util.Properties;

import com.elphel.imagej.cameras.EyesisCameraParameters;
import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.lwir.LwirReaderParameters;

public class CalibrationIllustrationParameters {
	static double          DFLT_LWIR_LO = 22500.0;
	static double          DFLT_LWIR_HI = 23500.0;
	static double          DFLT_EO_R2G =  1.03;  // gain red relative to green
	static double          DFLT_EO_B2G =  1.03;  // gain blue relative to green
	static double          DFLT_EO_HI =   255.0; // range to map to 255.0
	LwirReaderParameters   lwirReaderParameters;
	EyesisCameraParameters eyesisCameraParameters;
	
	double [][]            lwir_ranges; //  = new double [lwirReaderParameters.getLwirChannels(false).length][2];
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
		if (properties.getProperty(prefix+"eo_gamma")!=null)           this.eo_gamma = Double.parseDouble(properties.getProperty(prefix+"eo_gamma"));
		if (properties.getProperty(prefix+"eo_minlin_gamma")!=null)    this.eo_minlin_gamma = Double.parseDouble(properties.getProperty(prefix+"eo_minlin_gamma"));
		if (properties.getProperty(prefix+"eo_saturation")!=null)      this.eo_saturation = Double.parseDouble(properties.getProperty(prefix+"eo_saturation"));
		if (properties.getProperty(prefix+"threshold_contrast")!=null) this.threshold_contrast = Double.parseDouble(properties.getProperty(prefix+"threshold_contrast"));
		if (properties.getProperty(prefix+"threshold_number")!=null)   this.threshold_number = Integer.parseInt(properties.getProperty(prefix+"threshold_number"));
		if (properties.getProperty(prefix+"station_sel")!=null)        this.station_sel = Integer.parseInt(properties.getProperty(prefix+"station_sel"));
		if (properties.getProperty(prefix+"station_in_filenames")!=null)        this.station_in_filenames = Boolean.parseBoolean(properties.getProperty(prefix+"station_in_filenames"));
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
			if (lcolor_grid_extra < 0) this.color_grid = null;
			else this.color_grid_extra = setLongColor(lcolor_grid_extra);
		}
		if (properties.getProperty(prefix+"line_width_eo")!=null)      this.line_width_eo = Integer.parseInt(properties.getProperty(prefix+"line_width_eo"));
		if (properties.getProperty(prefix+"line_width_lwir")!=null)    this.line_width_lwir = Integer.parseInt(properties.getProperty(prefix+"line_width_lwir"));
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
		gd.addNumericField("Line width for LWIR images", this.line_width_lwir,  0,3,"pix","line width for the grid in LWIR images");
		gd.addNumericField("Line width for LWIR images", this.line_width_eo,  0,3,"pix","line width for the grid in EO images");
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
	}
	
	
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
		this.station_in_filenames =  gd.getNextBoolean();
		this.save_png =              gd.getNextBoolean();
		this.JPEG_quality =    (int) gd.getNextNumber();
		this.threshold_contrast =    gd.getNextNumber();
		this.threshold_number= (int) gd.getNextNumber();
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
}
