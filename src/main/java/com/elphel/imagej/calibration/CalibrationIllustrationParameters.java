package com.elphel.imagej.calibration;

import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.lwir.LwirReaderParameters;

public class CalibrationIllustrationParameters {
	double               dflt_lwir_lo = 22500.0;
	double               dflt_lwir_hi = 23500.0;
	LwirReaderParameters lwirReaderParameters;
	double [][]          lwir_ranges; //  = new double [lwirReaderParameters.getLwirChannels(false).length][2];
	int                  palette = 0; // 0 - white - hot, 1 - black - hot, 2+ - colored
	String               src_chn_prefix="src_chn-";
	boolean              save_png = true;
	int                  JPEG_quality = 90;
	String               channel_dir_prefix = "chn_";
	
	public CalibrationIllustrationParameters (LwirReaderParameters lwirReaderParameters) {
		this.lwirReaderParameters = lwirReaderParameters;
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

	}
	public void dialogQuestions(GenericJTabbedDialog gd) {
		for (int i = 0; i < lwir_ranges.length; i++) {
			gd.addNumericField("LWIR chn:"+i+" low range", this.lwir_ranges[i][0],  0,8,"","LWIR sensor range low level ");
			gd.addNumericField("LWIR chn:"+i+" high range", this.lwir_ranges[i][1],  0,8,"","LWIR sensor range high level ");
		}
		gd.addNumericField("Thermal color palette", this.palette,  0,3,"","0 - white-hot, 1 - black-hot, 2+ - colored");
		gd.addCheckbox("Save as PNG instead of JPEG", save_png);
		gd.addNumericField("JPEG quality", this.JPEG_quality,  0,3,"","Jpeg quality, 0 - use Tiff");
    	gd.addStringField ("Channel directory prefix",this.channel_dir_prefix, 15,"prefix to a directory name to save channel annotated files");
	}
	
	public void dialogAnswers(GenericJTabbedDialog gd) {
		for (int i = 0; i < lwir_ranges.length; i++) {
			this.lwir_ranges[i][0] = gd.getNextNumber();
			this.lwir_ranges[i][1] = gd.getNextNumber();
		}
		this.palette =      (int) gd.getNextNumber();
		this.save_png =           gd.getNextBoolean();
		this.JPEG_quality = (int) gd.getNextNumber();
		this.channel_dir_prefix = gd.getNextString();
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
				this.lwir_ranges[i][0] = dflt_lwir_lo;
				this.lwir_ranges[i][1] = dflt_lwir_hi;
			}
		}
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
