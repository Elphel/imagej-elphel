/**
 ** -----------------------------------------------------------------------------**
 ** LwirReaderParameters.java
 **
 ** Parameters handling for LwirReader class
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwirReaderParameters.java is free software: you can redistribute it and/or modify
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
package com.elphel.imagej.lwir;

import java.io.File;
import java.io.FilenameFilter;
import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;

import ij.Prefs;

public class LwirReaderParameters {
	private boolean parameters_updated = false;
	protected int     avg_number =            4; // number of measurements to average
	protected boolean lwir_ffc =              true;
	protected boolean avg_all =               true;
	protected String  lwir_ip =               "192.168.0.36";
	protected String  vnir_ip =               "192.168.0.38";
	protected int []  lwir_channels =         {0, 1, 2 ,3};
	protected int []  vnir_channels =         {0, 1, 2 ,3};
	protected boolean lwir_telemetry =        true;
	protected double  vnir_quality =          98.0;
	protected boolean vnir_scale =           false; // restore sensor pixel values, undo camera white balancing
	protected boolean vnir_autoexp =         false;
	protected double  vnir_max_autoexp_ms =  20.0;
	protected double  vnir_exposure_ms =      5.0;
	protected boolean vnir_whitebal =        false;
	protected double  vnir_gain_g =           2.0;
	protected double  vnir_gain_rg =          0.7705; // 1.116; halogen/fluorescent
	protected double  vnir_gain_bg =          2.401;  // 1.476;
	protected boolean [] selected_channels =  {true, true, true, true, true, true, true, true};

/*
	protected double [] vnir_exp_corr = {1.0, 1.0, 1.0, 1.0};
	protected double [] vnir_gcorr_rbgb = {
			1.0,    1.0,    1.0,
			1.0,    1.0,    1.0,
			1.0,    1.0,    1.0,
			1.0,    1.0,    1.0};

 */
	protected double [] vnir_exp_corr =	{
			1.0, 1.0026, 0.9868, 1.0211};
	protected double [] vnir_gcorr_rbgb = { // example after autobalance for 192.168.0.38 with halogen lamps
			1.0,    1.0,    0.9743,
			1.0668, 1.1055, 1.0006,
			1.1533, 1.0780, 1.0015,
			0.9966, 1.0445, 1.0023};
	protected int     lwir_trig_dly   =     9000000; //in 100MHz clock cycle, current FPGA requires >0 (used 1000)
	protected int     vnir_lag   =          2;//1; // frames
	protected double  max_mismatch_ms =     0.05;
	protected int     max_frame_diff =      2;// 1; // 2;
	protected int     debug_level =           0;//-3: OFF, -2:Fatal, -1:ERROR, 0:WARN, 1:INFO,2:DEBUG

	// --- interface methods
	public int getDebugLevel() {
		return this.debug_level;
	}
	public void setDebugLevel(int level) {
		this.debug_level = level;
	}

	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"avg_number",          this.avg_number+"");
		properties.setProperty(prefix+"lwir_ffc",            this.lwir_ffc+"");
		properties.setProperty(prefix+"avg_all",             this.avg_all+"");
		properties.setProperty(prefix+"lwir_ip",             this.lwir_ip+"");
		properties.setProperty(prefix+"vnir_ip",             this.vnir_ip+"");
		properties.setProperty(prefix+"lwir_channels",       arr_to_str(this.lwir_channels));
		properties.setProperty(prefix+"vnir_channels",       arr_to_str(this.vnir_channels));
		properties.setProperty(prefix+"lwir_telemetry",      this.lwir_telemetry+"");
		properties.setProperty(prefix+"vnir_quality",        this.vnir_quality+"");
		properties.setProperty(prefix+"vnir_scale",          this.vnir_scale+"");
		properties.setProperty(prefix+"vnir_autoexp",        this.vnir_autoexp+"");
		properties.setProperty(prefix+"vnir_max_autoexp_ms", this.vnir_max_autoexp_ms+"");
		properties.setProperty(prefix+"vnir_exposure_ms",    this.vnir_exposure_ms+"");
		properties.setProperty(prefix+"vnir_whitebal",       this.vnir_whitebal+"");
		properties.setProperty(prefix+"vnir_gain_g",         this.vnir_gain_g+"");
		properties.setProperty(prefix+"vnir_gain_rg",        this.vnir_gain_rg+"");
		properties.setProperty(prefix+"vnir_gain_bg",        this.vnir_gain_bg+"");
		properties.setProperty(prefix+"vnir_exp_corr",       arr_to_str(this.vnir_exp_corr));
		properties.setProperty(prefix+"vnir_gcorr_rbgb",     arr_to_str(this.vnir_gcorr_rbgb));
		properties.setProperty(prefix+"lwir_trig_dly",       this.lwir_trig_dly+"");
		properties.setProperty(prefix+"vnir_lag",            this.vnir_lag+"");
		properties.setProperty(prefix+"max_mismatch_ms",     this.max_mismatch_ms+"");
		properties.setProperty(prefix+"max_frame_diff",      this.max_frame_diff+"");
		properties.setProperty(prefix+"debug_level",         this.debug_level+"");
		properties.setProperty(prefix+"selected_channels",   arr_to_str(this.selected_channels));

	}

	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"avg_number")!=null)          this.avg_number=Integer.parseInt(properties.getProperty(prefix+"avg_number"));
		if (properties.getProperty(prefix+"lwir_ffc")!=null)            this.lwir_ffc= Boolean.parseBoolean(properties.getProperty(prefix+"lwir_ffc"));
		if (properties.getProperty(prefix+"avg_all")!=null)             this.avg_all= Boolean.parseBoolean(properties.getProperty(prefix+"avg_all"));
		if (properties.getProperty(prefix+"lwir_ip")!=null)             this.lwir_ip=      properties.getProperty(prefix+"lwir_ip");
		if (properties.getProperty(prefix+"vnir_ip")!=null)             this.vnir_ip=      properties.getProperty(prefix+"vnir_ip");
		if (properties.getProperty(prefix+"lwir_channels")!=null)       this.lwir_channels=str_to_iarr(properties.getProperty(prefix+"lwir_channels"));
		if (properties.getProperty(prefix+"vnir_channels")!=null)       this.vnir_channels=str_to_iarr(properties.getProperty(prefix+"vnir_channels"));
		if (properties.getProperty(prefix+"lwir_telemetry")!=null)      this.lwir_telemetry= Boolean.parseBoolean(properties.getProperty(prefix+"lwir_telemetry"));
		if (properties.getProperty(prefix+"vnir_quality")!=null)        this.vnir_quality=Double.parseDouble(properties.getProperty(prefix+"vnir_quality"));
		if (properties.getProperty(prefix+"vnir_scale")!=null)          this.vnir_scale= Boolean.parseBoolean(properties.getProperty(prefix+"vnir_scale"));
		if (properties.getProperty(prefix+"vnir_autoexp")!=null)        this.vnir_autoexp= Boolean.parseBoolean(properties.getProperty(prefix+"vnir_autoexp"));
		if (properties.getProperty(prefix+"vnir_max_autoexp_ms")!=null) this.vnir_max_autoexp_ms=Double.parseDouble(properties.getProperty(prefix+"vnir_max_autoexp_ms"));
		if (properties.getProperty(prefix+"vnir_exposure_ms")!=null)    this.vnir_exposure_ms=Double.parseDouble(properties.getProperty(prefix+"vnir_exposure_ms"));
		if (properties.getProperty(prefix+"vnir_whitebal")!=null)       this.vnir_whitebal= Boolean.parseBoolean(properties.getProperty(prefix+"vnir_whitebal"));
		if (properties.getProperty(prefix+"vnir_gain_g")!=null)         this.vnir_gain_g=Double.parseDouble(properties.getProperty(prefix+"vnir_gain_g"));
		if (properties.getProperty(prefix+"vnir_gain_rg")!=null)        this.vnir_gain_rg=Double.parseDouble(properties.getProperty(prefix+"vnir_gain_rg"));
		if (properties.getProperty(prefix+"vnir_gain_bg")!=null)        this.vnir_gain_bg=Double.parseDouble(properties.getProperty(prefix+"vnir_gain_bg"));
		if (properties.getProperty(prefix+"vnir_exp_corr")!=null)       this.vnir_exp_corr=str_to_darr(properties.getProperty(prefix+"vnir_exp_corr"));
		if (properties.getProperty(prefix+"vnir_gcorr_rbgb")!=null)     this.vnir_gcorr_rbgb=str_to_darr(properties.getProperty(prefix+"vnir_gcorr_rbgb"));
		if (properties.getProperty(prefix+"lwir_trig_dly")!=null)       this.lwir_trig_dly=Integer.parseInt(properties.getProperty(prefix+"lwir_trig_dly"));
		if (properties.getProperty(prefix+"vnir_lag")!=null)            this.vnir_lag=Integer.parseInt(properties.getProperty(prefix+"vnir_lag"));
		if (properties.getProperty(prefix+"max_mismatch_ms")!=null)     this.max_mismatch_ms=Double.parseDouble(properties.getProperty(prefix+"max_mismatch_ms"));
		if (properties.getProperty(prefix+"max_frame_diff")!=null)      this.max_frame_diff=Integer.parseInt(properties.getProperty(prefix+"max_frame_diff"));
		if (properties.getProperty(prefix+"debug_level")!=null)         this.debug_level=Integer.parseInt(properties.getProperty(prefix+"debug_level"));
		if (properties.getProperty(prefix+"selected_channels")!=null)   this.selected_channels=str_to_barr(properties.getProperty(prefix+"selected_channels"));
		parameters_updated = true;
	}
	@Override
	public LwirReaderParameters clone() { //  throws CloneNotSupportedException {
		LwirReaderParameters lrp =       new LwirReaderParameters();
		lrp.avg_number=                this.avg_number;
		lrp.lwir_ffc=                  this.lwir_ffc;
		lrp.avg_all=                   this.avg_all;
		lrp.lwir_ip=                   this.lwir_ip;
		lrp.vnir_ip=                   this.vnir_ip;
		lrp.lwir_channels=             this.lwir_channels.clone();
		lrp.vnir_channels=             this.vnir_channels.clone();
		lrp.lwir_telemetry=            this.lwir_telemetry;
		lrp.vnir_quality=              this.vnir_quality;
		lrp.vnir_scale=                this.vnir_scale;
		lrp.vnir_autoexp=              this.vnir_autoexp;
		lrp.vnir_max_autoexp_ms=       this.vnir_max_autoexp_ms;
		lrp.vnir_exposure_ms=          this.vnir_exposure_ms;
		lrp.vnir_whitebal=             this.vnir_whitebal;
		lrp.vnir_gain_g=               this.vnir_gain_g;
		lrp.vnir_gain_rg=              this.vnir_gain_rg;
		lrp.vnir_gain_bg=              this.vnir_gain_bg;
		lrp.vnir_exp_corr=             this.vnir_exp_corr.clone();
		lrp.vnir_gcorr_rbgb=           this.vnir_gcorr_rbgb.clone();
		lrp.lwir_trig_dly=             this.lwir_trig_dly;
		lrp.vnir_lag=                  this.vnir_lag;
		lrp.max_mismatch_ms=           this.max_mismatch_ms;
		lrp.max_frame_diff=            this.max_frame_diff;
		lrp.debug_level=               this.debug_level;
		lrp.selected_channels =        this.selected_channels.clone();
		return lrp;
	}

	@Override
	public boolean equals(Object o) {
		if (o == this) {
			return true;
		}
		if (!(o instanceof LwirReaderParameters)) {
			return false;
		}
		LwirReaderParameters lrp = (LwirReaderParameters) o;
		return  (lrp.avg_number == this.avg_number) &&
				(lrp.lwir_ffc == this.lwir_ffc) &&
				(lrp.avg_all == this.avg_all) &&
				(lrp.lwir_ip.equals(this.lwir_ip)) &&
				(lrp.vnir_ip.equals(this.vnir_ip)) &&
				(java.util.Arrays.equals(lrp.lwir_channels, this.lwir_channels)) &&
				(java.util.Arrays.equals(lrp.vnir_channels, this.vnir_channels)) &&
				(lrp.lwir_telemetry == this.lwir_telemetry) &&
				(lrp.vnir_quality == this.vnir_quality) &&
				(lrp.vnir_scale == this.vnir_scale) &&
				(lrp.vnir_autoexp == this.vnir_autoexp) &&
				(lrp.vnir_max_autoexp_ms == this.vnir_max_autoexp_ms) &&
				(lrp.vnir_exposure_ms == this.vnir_exposure_ms) &&
				(lrp.vnir_whitebal == this.vnir_whitebal) &&
				(lrp.vnir_gain_g == this.vnir_gain_g) &&
				(lrp.vnir_gain_rg == this.vnir_gain_rg) &&
				(lrp.vnir_gain_bg == this.vnir_gain_bg) &&
				(java.util.Arrays.equals(lrp.vnir_exp_corr, this.vnir_exp_corr)) &&
				(java.util.Arrays.equals(lrp.vnir_gcorr_rbgb, this.vnir_gcorr_rbgb)) &&
				(lrp.lwir_trig_dly == this.lwir_trig_dly) &&
				(lrp.vnir_lag == this.vnir_lag) &&
				(lrp.max_mismatch_ms == this.max_mismatch_ms) &&
				(lrp.max_frame_diff == this.max_frame_diff) &&
				(lrp.debug_level ==    this.debug_level) &&
				(java.util.Arrays.equals(lrp.selected_channels, this.selected_channels));
	}

	@Override
	public int hashCode() {
		int prime = 31;
		int result = 1;
		result = prime * result + (new Integer(avg_number)).hashCode();
		result = prime * result + (lwir_ffc?1:0);
		result = prime * result + (avg_all?1:0);
		result = prime * result + lwir_ip.hashCode();
		result = prime * result + vnir_ip.hashCode();
		result = prime * result + arr_to_str(lwir_channels).hashCode();
		result = prime * result + arr_to_str(vnir_channels).hashCode();
		result = prime * result + (lwir_telemetry?1:0);
		result = prime * result + (new Double(vnir_quality)).hashCode();
		result = prime * result + (vnir_scale?1:0);
		result = prime * result + (vnir_autoexp?1:0);
		result = prime * result + (new Double(vnir_max_autoexp_ms)).hashCode();
		result = prime * result + (new Double(vnir_exposure_ms)).hashCode();
		result = prime * result + (vnir_whitebal?1:0);
		result = prime * result + (new Double(vnir_gain_g)).hashCode();
		result = prime * result + (new Double(vnir_gain_rg)).hashCode();
		result = prime * result + (new Double(vnir_gain_bg)).hashCode();
		result = prime * result + arr_to_str(vnir_exp_corr).hashCode();
		result = prime * result + arr_to_str(vnir_gcorr_rbgb).hashCode();
		result = prime * result + (new Integer(lwir_trig_dly)).hashCode();
		result = prime * result + arr_to_str(selected_channels).hashCode();
		// next are not needed to be programmed to the cameras
//		result = prime * result + (new Integer(vnir_lag)).hashCode();
//		result = prime * result + (new Double(max_mismatch_ms)).hashCode();
//		result = prime * result + (new Integer(max_frame_diff)).hashCode();
		return 0;
	}

	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addNumericField("Number to average",this.avg_number,     0,3,"","Number of acquired consecutive images to average");
		gd.addCheckbox    ("Run FFC",       this.lwir_ffc, "Perform calibration before each measurements to average (takes ~1.6 sec, 15 frames)");
		gd.addCheckbox    ("Average all",   this.avg_all, "Average all simultaneously acquired images (unchecked - only requested number to average)");
		gd.addStringField ("LWIR IP",       this.lwir_ip, 20, "LWIR camera IP address");
		gd.addStringField ("VNIR IP",       this.vnir_ip, 20, "Visible range high resolution camera IP address");
		gd.addStringField ("LWIR channels", arr_to_str(this.lwir_channels), 20, "Space-separated list of used LWIR camera channels, such as '0 1 2 3'");
		gd.addStringField ("VNIR channels", arr_to_str(this.vnir_channels), 20, "Space-separated list of used visible range camera channels, such as '0 1 2 3'");
		gd.addCheckbox    ("LWIR telemetry",    this.lwir_telemetry, "Set LWIR sesnors to provide telemetry data in the last 2 lines (may become mandatory later)");
		gd.addNumericField("VNIR quality",  this.vnir_quality,  3,6,"%", "Visible range camera JPEG compression quality (all channels)");
		gd.addCheckbox    ("VNIR undo white balance",    this.vnir_scale, "Undo in-camera white balancing");
		gd.addCheckbox    ("VNIR autoexposure", this.vnir_autoexp, "Enable autoexposure for the visible range camera");
		gd.addNumericField("VNIR vnir_max_autoexp_ms", this.vnir_max_autoexp_ms,  3,6,"ms", "Visible range camera maximal exposure in autoexposure mode");
		gd.addNumericField("VNIR exposure", this.vnir_exposure_ms,  3,6,"ms", "Visible range camera exposure time (all channels)");
		gd.addCheckbox    ("VNIR white balance", this.vnir_whitebal, "Enable automatic white balancing for the visible range camera");
		gd.addNumericField("VNIR gain G",   this.vnir_gain_g,  3,6,"","Analog gain for green channel for all visible range camera channels (normally 2.0)");
		gd.addNumericField("VNIR gain R/G", this.vnir_gain_rg, 3,6,"","Red to green gain ratio for all visible range camera channels");
		gd.addNumericField("VNIR gain B/G", this.vnir_gain_bg, 3,6,"","Blue to green gain ratio for all visible range camera channels");
		gd.addStringField ("VNIR exposuere corrections", arr_to_str(this.vnir_exp_corr), 50, "Fine corrections of channel exposures (4 channel relative exposures)");
		gd.addStringField ("VNIR gain corrections", arr_to_str(this.vnir_gcorr_rbgb), 100, "Fine corrections to per channel, per color gains:'r0 b0 gb0 r1 b1 gb1 ...'");
		gd.addNumericField("LWIR trig dly", this.lwir_trig_dly,     0,10,"x10ns","Output trigger delay, should eventually match Lepton+FPGA latency to trigger VNIR exactly 1 frame after LWIR. 0 does not work with current FPGA - usec do not match sec in transmitted timestamp");
		gd.addNumericField("VNIR lag",      this.vnir_lag,     0,3,"","Visible camera lag (in frames) relative to LWIR one");
		gd.addNumericField("Max mismatch",  this.max_mismatch_ms,    3,6,"ms","Maximal mismatch between image timestamps. Larger mismatch requires LWIR sinsor reinitialization");
		gd.addNumericField("Max frame diff",this.max_frame_diff,     0,3,"","Maximal difference in frames between simultaneously acquired channels as calculated from the timestamps");
		gd.addNumericField("Debug level",   this.debug_level,        0,3,"","Image acquisition log level: -3: OFF, -2:FATAL, -1:ERROR, 0:WARN, 1:INFO, 2:DEBUG");
		gd.addStringField ("Selected channels", arr_to_str(this.selected_channels), 20, "Space-separated channel selection (1 - selected, 0 - unselected)");
	}

	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.avg_number =             (int) gd.getNextNumber();
		this.lwir_ffc =                     gd.getNextBoolean();
		this.avg_all =                      gd.getNextBoolean();
		this.lwir_ip =                      gd.getNextString();
		this.vnir_ip =                      gd.getNextString();
		this.lwir_channels =                str_to_iarr(gd.getNextString());
		this.vnir_channels =                str_to_iarr(gd.getNextString());
		this.lwir_telemetry =               gd.getNextBoolean();
		this.vnir_quality =                 gd.getNextNumber();
		this.vnir_scale =                   gd.getNextBoolean();
		this.vnir_autoexp =                 gd.getNextBoolean();
		this.vnir_max_autoexp_ms =          gd.getNextNumber();
		this.vnir_exposure_ms =             gd.getNextNumber();
		this.vnir_whitebal =                gd.getNextBoolean();
		this.vnir_gain_g =                  gd.getNextNumber();
		this.vnir_gain_rg =                 gd.getNextNumber();
		this.vnir_gain_bg =                 gd.getNextNumber();
		this.vnir_exp_corr =                str_to_darr(gd.getNextString());
		this.vnir_gcorr_rbgb =              str_to_darr(gd.getNextString());
		this.lwir_trig_dly =          (int) gd.getNextNumber();
		this.vnir_lag =               (int) gd.getNextNumber();
		this.max_mismatch_ms =              gd.getNextNumber();
		this.max_frame_diff =         (int) gd.getNextNumber();
		this.debug_level =            (int) gd.getNextNumber();
		this.selected_channels =            str_to_barr(gd.getNextString());
		parameters_updated = true;
	}


	public boolean showJDialog() {
		GenericJTabbedDialog gd = new GenericJTabbedDialog("Set CLT parameters",800,900);
//		gd.addTab         ("General", "General parameters");
		dialogQuestions(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		dialogAnswers(gd);
		return true;
	}


	public boolean is_updated() {
		return parameters_updated;
	}

	public void reset_updated() {
		parameters_updated = false;
	}

	public boolean [] getSelected(){
		return  selected_channels;
	}

	public boolean [] getSelectedLwir(){
		boolean [] sel = selected_channels.clone();
		for (int i = lwir_channels.length; i < sel.length; i++ ) {
			sel[i] = false;
		}
		return  sel;
	}

	public boolean [] getSelectedVnir(){
		boolean [] sel = selected_channels.clone();
		for (int i = 0; i < lwir_channels.length; i++ ) {
			sel[i] = false;
		}
		return  sel;
	}

	public String [] getSourceFilesFlat(String [] sets, boolean[] channels) {
		int num_sel = 0;
		for (boolean s: channels) if (s) num_sel++;
		String [] files = new String [sets.length * num_sel];
		int indx = 0;
		for (String set:sets) {
			for (int i = 0; i < channels.length; i++) if (channels[i]) {
				String file_base = set;
				if (set.indexOf(Prefs.getFileSeparator()) >=0) {
					file_base = set.substring(set.lastIndexOf(Prefs.getFileSeparator())+1);
				}
			files[indx++] =set + Prefs.getFileSeparator()+ file_base+"_"+i+".tiff";
			}
		}
		return files;
	}


	// Image set may have different timestamps, only lwir0 matches
	public String [][] getSourceFiles(String [] sets, boolean[] channels) {
		String [][] files = new String [sets.length][channels.length];
		for (int nset= 0; nset < sets.length; nset++) {
			// read all files in the directory
			File  set_dir = new File(sets[nset]);
			File [] channel_files = set_dir.listFiles(new FilenameFilter() {
				@Override
				public boolean accept(File current, String name) {
					if (!(new File(current, name).isFile()) || !name.endsWith(".tiff")) return false;
					String base = name.substring(0, name.lastIndexOf(".tiff"));
					int undr = base.lastIndexOf("_");
					if (undr < 0) return false;
					int chn = -1;
					try {
						chn = Integer.parseInt(base.substring(undr + 1));
					} catch (Exception e) {

					}
					return chn >= 0;
				}
			});
			for (File f: channel_files) {
				String base = f.getName().substring(0, f.getName().lastIndexOf(".tiff"));
				int undr = base.lastIndexOf("_");
				int chn =  Integer.parseInt(base.substring(undr + 1));
				if ((chn >= 0) && (chn < channels.length) && channels[chn]) {
					files[nset][chn] = f.getPath();
				}
			}
		}
		return files;
	}



	// --- internal methods

	private String arr_to_str(int [] arr) {
		String s = "";
		for (int c:arr) s+=c+" ";
		return s.trim();
	}

	private String arr_to_str(double [] arr) {
		String s = "";
		for (double c:arr) s+=c+" ";
		return s.trim();
	}


	private String arr_to_str(boolean [] arr) {
		String s = "";
		for (boolean c:arr) s+= (c?1:0)+" ";
		return s.trim();
	}

	private boolean [] str_to_barr(String s) {
		int [] iarr = str_to_iarr(s);
		if (iarr == null) return null;
		boolean [] barr = new boolean [iarr.length];
		for (int i = 0; i < barr.length; i++) {
			barr[i] = iarr[i] != 0;
		}
		return barr;
	}

	private int [] str_to_iarr(String s) {
		String [] sa;
		if (s.indexOf(",") >= 0) {
			sa = s.trim().split(",");
		} else {
			sa = s.trim().split(" ");
		}
		int [] iarr = new int [sa.length];
		for (int i = 0; i < sa.length; i++) {
			iarr[i]=Integer.parseInt(sa[i].trim());
		}
		return iarr;
	}

	private double [] str_to_darr(String s) {
		String [] sa;
		if (s.indexOf(",") >= 0) {
			sa = s.trim().split(",");
		} else {
			sa = s.trim().split(" ");
		}
		double [] darr = new double [sa.length];
		for (int i = 0; i < sa.length; i++) {
			darr[i]=Double.parseDouble(sa[i].trim());
		}
		return darr;
	}


	public boolean [] selectSourceChannels(){
		return selectSourceChannels(selected_channels);
	}

	public boolean [] selectSourceChannels(boolean [] sel) {
		GenericJTabbedDialog gd = new GenericJTabbedDialog("Set CLT parameters",300,500);
		for (int i = 0; i < sel.length; i++) {
			gd.addCheckbox    ("Channel "+i,       sel[i], "Enable processing camera channel "+i);
		}
		gd.showDialog();
		if (gd.wasCanceled()) return null;
		for (int i = 0; i < sel.length; i++) {
			sel[i] = gd.getNextBoolean();
		}
    	return sel;
    }

}
