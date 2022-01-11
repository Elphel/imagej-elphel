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
import com.elphel.imagej.readers.EyesisTiff;

import ij.ImagePlus;
import ij.Prefs;

public class LwirReaderParameters {
	public final static String    NAME_TALON=        "Talon";
	public final static String    NAME_LWIR16=       "LWIR16";
	public final static String [] CAMERA_NAMES= {NAME_TALON,NAME_LWIR16};
	public final static int []    BOSON_FFC_FRAMES= {2,4,8,16};
	public final static int []    FFC_GROUPS=          {1,2,4};
	public static final String [] SENSOR_TYPES = {"EO","LWIR"};
	public static final String    SENSOR_TYPE =  "SENSOR_TYPE";
	public final static int       TYPE_EO= 0;
	public final static int       TYPE_LWIR= 1;
	
	protected static int MAX_LWIR_WIDTH =     1024; //
	
	private boolean parameters_updated =     false;
	protected String  camera_name =          "Talon"; // "LWIR16";
	protected int     avg_number =           4; // number of LWIR measurements to average
	protected int     avg_number_eo =        1; // number of EO measurements to average
	protected int     num_frames =           7; // misses images if high quality large size?
	protected boolean lwir_ffc =             true;
	protected int     ffc_frames =           8;
	protected int     ffc_wait_frames =      10; // extra wait after FFC finished (2xffc_frames)
	protected int     ffc_groups =           2; // 1/2/4 - do not run FFC on all channels simultaneously (43 failed)
	protected double  acquire_delay =        3.0; // seconds. Schedule later from the master channel current time to start all cameras
	protected boolean avg_all =              true;
//	protected String  lwir_ip =              "192.168.0.36";
	protected String[] lwir_ips =             {"192.168.0.41","192.168.0.42","192.168.0.43","192.168.0.44"}; // first - master
	protected String  eo_ip =                "192.168.0.38";
	protected int []  lwir_channels =        {0, 1, 2 ,3};
	protected int []  eo_channels =          {0, 1, 2 ,3};
	protected boolean lwir_telemetry =       true;
	protected boolean eo_full_window =       true;
	protected double  eo_quality =           98.0;
	protected boolean eo_scale =             false; // restore sensor pixel values, undo camera white balancing
	protected boolean eo_autoexp =           false;
	protected double  eo_max_autoexp_ms =    20.0;
	protected double  eo_exposure_ms =       5.0;
	protected boolean eo_whitebal =          false;
	protected double  eo_gain_g =            2.0;
	protected double  eo_gain_rg =           0.7705; // 1.116; halogen/fluorescent
	protected double  eo_gain_bg =           2.401;  // 1.476;
	protected boolean [] selected_channels = {true, true, true, true, true, true, true, true};
/*
	protected double [] eo_exp_corr = {1.0, 1.0, 1.0, 1.0};
	protected double [] eo_gcorr_rbgb = {
			1.0,    1.0,    1.0,
			1.0,    1.0,    1.0,
			1.0,    1.0,    1.0,
			1.0,    1.0,    1.0};
 */
	protected double [] eo_exp_corr =	{
			1.0, 1.0, 1.0, 1.0};
	protected double [] eo_gcorr_rbgb = { // example after autobalance for 192.168.0.38 with halogen lamps
			1.0,    1.0,    0.9743,
			1.0668, 1.1055, 1.0006,
			1.1533, 1.0780, 1.0015,
			0.9966, 1.0445, 1.0023};
	protected int     lwir_trig_dly   =     9000000; //in 100MHz clock cycle, current FPGA requires >0 (used 1000)
	protected int     eo_lag   =          2;//1; // frames
	protected double  max_mismatch_ms =     0.05;
	protected int     max_frame_diff =      2;// 1; // 2;
	protected int     debug_level =           0;//-3: OFF, -2:Fatal, -1:ERROR, 0:WARN, 1:INFO,2:DEBUG
	protected boolean show_images = false;

	protected int     lwir_chn0 =             0; // not configurable
//	protected int     eo_chn0 =               4; // not configurable

	public LwirReaderParameters() {
		
	}
	
	public int [] getLwirChannels(boolean absolote) {
		int [] absolute_chn = new int [lwir_channels.length];
		for (int i = 0; i < absolute_chn.length; i++) {
			absolute_chn[i] = lwir_channels[i] + (absolote? getLwirChn0() : 0);
		}
		return absolute_chn;
	}

	public int [] getEoChannels(boolean absolote) {
		int [] absolute_chn = new int [eo_channels.length];
		for (int i = 0; i < absolute_chn.length; i++) {
			absolute_chn[i] = eo_channels[i] +  (absolote? getEoChn0():0);
		}
		return absolute_chn;
	}
//	protected int []  lwir_channels =        {0, 1, 2 ,3};
//	protected int []  eo_channels =          {0, 1, 2 ,3};
	
	
	public LwirReaderParameters(String name) {
		if (NAME_TALON.equals(name)) camera_name = NAME_TALON;
		else if (NAME_LWIR16.equals(name)) camera_name = NAME_LWIR16;
	}

	public static boolean is_LWIR(String type) {
		return ((type != null) && type.equals(SENSOR_TYPES[1]));
	}
	public static boolean is_EO(String type) {
		return ((type != null) && type.equals(SENSOR_TYPES[1]));
	}
	public static String sensorName(boolean lwir) {
		return SENSOR_TYPES[lwir?1:0];
	}
	
	// --- interface methods
	public int getAvgNumberLwir() {
		return avg_number;
	}
	public int getAvgNumberEO() {
		return avg_number_eo;
	}

	
	public String getCameraName() {
		return this.camera_name;
	}
	
	public int getNumFrames() {
		return num_frames;
	}
	public int getLwirChn0() {
		return lwir_chn0;
	}

	public int getEoChn0() {
		return isLwir16() ? 16:4;// eo_chn0;
	}
	
	public int getNumChannels() {
		return isLwir16()?20:8;
	}
	
	public int [] getTypeMap() { // eo - 0, lwir - 1
		int [] types = new int [getNumChannels()];
		for (int i = 0; i < types.length; i++) {
			types[i] = (i >= getEoChn0())? TYPE_EO : TYPE_LWIR;
		}
		return types;
	}
	
	public static boolean is_LWIR(int width) {
		return width <= MAX_LWIR_WIDTH;
	}

	public static boolean is_LWIR(ImagePlus imp){
		// See if image has LwirReaderParameters.SENSOR_TYPE property, then use is_LWIR(String property_value),
		// if not - use old width property
		if (imp.getProperty("WOI_WIDTH")==null) {
			EyesisTiff.decodeProperiesFromInfo(imp);
		}
		if (imp.getProperty(SENSOR_TYPE)!=null) {
			return is_LWIR((String) imp.getProperty(SENSOR_TYPE));
		}
		return is_LWIR(imp.getWidth());
	}
	public static int sensorType(ImagePlus imp) {
		if (imp.getProperty("WOI_WIDTH")==null) {
			EyesisTiff.decodeProperiesFromInfo(imp);
		}
		if (imp.getProperty(SENSOR_TYPE)!=null) {
			return (is_LWIR((String) imp.getProperty(SENSOR_TYPE)))? 1:0;
		}
		return is_LWIR(imp.getWidth())? TYPE_LWIR : TYPE_EO;
		
	}

	public int getDebugLevel() {
		return this.debug_level;
	}

	public void setDebugLevel(int level) {
		this.debug_level = level;
	}

	public boolean isShowImages() {
		return show_images;
	}
	public String getLwirIP() {
		return lwir_ips[0];
	}
	public String getLwirIP(int n) {
		return lwir_ips[n];
	}
	public String [] getLwirIPs() {
		return lwir_ips;
	}

	public String getEOIP() {
		return eo_ip;
	}
	public String [] getAllIPs() {
		String [] all_ips = new String [lwir_ips.length+1];
		int i = 0;
		for (i = 0; i < lwir_ips.length; i++) {
			all_ips[i] = lwir_ips[i];
		}
		all_ips[i] = eo_ip;
		return all_ips;
	}
	boolean [] getSelectedChannels() {
		return selected_channels;
	}
	
	// call after switch from Lwir16 to Talon
	public boolean isLwir16() {
		return camera_name.equals(NAME_LWIR16) || ((lwir_channels != null) && (lwir_channels.length == 16));
	}
	
	public int getFFCFrames() {
		return ffc_frames;
	}

	public int getFFCGroups() {
		return ffc_groups;
	}

	public int getFFCWaitFrames() {
		return ffc_wait_frames;
	}
	
	public double getAcquireDelay() {
		return acquire_delay;
	}
	
	public boolean runFFCBeforeMeasurments() {
		return lwir_ffc;
	}

	public void setDefaultChannels() {
		boolean[] old_selected = selected_channels;
		if (isLwir16()) {
			lwir_ips =      new String[]{"192.168.0.41","192.168.0.42","192.168.0.43","192.168.0.43"}; 
			lwir_channels = new int[] {0,1,2,3,4,5,6,7,8,9,10,11,12,12,14,15};
			eo_ip =         "192.168.0.45";
		} else {
			lwir_ips =      new String[]{"192.168.0.36"}; 
			lwir_channels = new int[] {0,1,2,3};
			eo_ip =         "192.168.0.38";
		}
		selected_channels = new boolean[lwir_channels.length + eo_channels.length];
		for (int i = 0; i < selected_channels.length; i++) {
			if (i < old_selected.length) {
				selected_channels[i] = old_selected[i];
			} else {
				selected_channels[i] = true;
			}
		}
	}
	
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"camera_name",         this.camera_name+"");
		properties.setProperty(prefix+"avg_number",          this.avg_number+"");
		properties.setProperty(prefix+"avg_number_eo",       this.avg_number_eo+"");
		properties.setProperty(prefix+"num_frames",          this.num_frames+"");
		properties.setProperty(prefix+"lwir_ffc",            this.lwir_ffc+"");
		properties.setProperty(prefix+"ffc_frames",          this.ffc_frames+"");
		properties.setProperty(prefix+"ffc_wait_frames",     this.ffc_wait_frames+"");
		properties.setProperty(prefix+"ffc_groups",          this.ffc_groups+"");
		properties.setProperty(prefix+"acquire_delay",       this.acquire_delay+"");
		properties.setProperty(prefix+"avg_all",             this.avg_all+"");
		properties.setProperty(prefix+"lwir_ips",            arr_to_str(this.lwir_ips)); // this.lwir_ip+"");; this.lwir_ip+"");
		properties.setProperty(prefix+"eo_ip",               this.eo_ip+"");
		properties.setProperty(prefix+"lwir_channels",       arr_to_str(this.lwir_channels));
		properties.setProperty(prefix+"eo_channels",         arr_to_str(this.eo_channels));
		properties.setProperty(prefix+"lwir_telemetry",      this.lwir_telemetry+"");
		properties.setProperty(prefix+"eo_full_window",      this.eo_full_window+"");
		properties.setProperty(prefix+"eo_quality",          this.eo_quality+"");
		properties.setProperty(prefix+"eo_scale",            this.eo_scale+"");
		properties.setProperty(prefix+"eo_autoexp",          this.eo_autoexp+"");
		properties.setProperty(prefix+"eo_max_autoexp_ms",   this.eo_max_autoexp_ms+"");
		properties.setProperty(prefix+"eo_exposure_ms",      this.eo_exposure_ms+"");
		properties.setProperty(prefix+"eo_whitebal",         this.eo_whitebal+"");
		properties.setProperty(prefix+"eo_gain_g",           this.eo_gain_g+"");
		properties.setProperty(prefix+"eo_gain_rg",          this.eo_gain_rg+"");
		properties.setProperty(prefix+"eo_gain_bg",          this.eo_gain_bg+"");
		properties.setProperty(prefix+"eo_exp_corr",         arr_to_str(this.eo_exp_corr));
		properties.setProperty(prefix+"eo_gcorr_rbgb",       arr_to_str(this.eo_gcorr_rbgb));
		properties.setProperty(prefix+"lwir_trig_dly",       this.lwir_trig_dly+"");
		properties.setProperty(prefix+"eo_lag",              this.eo_lag+"");
		properties.setProperty(prefix+"max_mismatch_ms",     this.max_mismatch_ms+"");
		properties.setProperty(prefix+"max_frame_diff",      this.max_frame_diff+"");
		properties.setProperty(prefix+"debug_level",         this.debug_level+"");
		properties.setProperty(prefix+"selected_channels",   arr_to_str(this.selected_channels));
		properties.setProperty(prefix+"show_images",         this.show_images+"");
	}

	public void getProperties(String prefix,Properties properties){
		String old_camera_name =       this.camera_name; 
		if (properties.getProperty(prefix+"camera_name")!=null) {
			this.camera_name = properties.getProperty(prefix+"camera_name");
		} else {
			this.camera_name = NAME_TALON; // old configs
		}
		if (!this.camera_name.equals(old_camera_name)) {
			setDefaultChannels(); // modifies lwir_IPs, eo_ip, channels,
		}
		if (properties.getProperty(prefix+"vnir_ip")!=null)             this.eo_ip=      properties.getProperty(prefix+"vnir_ip");
		if (properties.getProperty(prefix+"vnir_channels")!=null)       this.eo_channels=str_to_iarr(properties.getProperty(prefix+"vnir_channels"));
		if (properties.getProperty(prefix+"vnir_quality")!=null)        this.eo_quality=Double.parseDouble(properties.getProperty(prefix+"vnir_quality"));
		if (properties.getProperty(prefix+"vnir_scale")!=null)          this.eo_scale= Boolean.parseBoolean(properties.getProperty(prefix+"vnir_scale"));
		if (properties.getProperty(prefix+"vnir_autoexp")!=null)        this.eo_autoexp= Boolean.parseBoolean(properties.getProperty(prefix+"vnir_autoexp"));
		if (properties.getProperty(prefix+"vnir_max_autoexp_ms")!=null) this.eo_max_autoexp_ms=Double.parseDouble(properties.getProperty(prefix+"vnir_max_autoexp_ms"));
		if (properties.getProperty(prefix+"vnir_exposure_ms")!=null)    this.eo_exposure_ms=Double.parseDouble(properties.getProperty(prefix+"vnir_exposure_ms"));
		if (properties.getProperty(prefix+"vnir_whitebal")!=null)       this.eo_whitebal= Boolean.parseBoolean(properties.getProperty(prefix+"vnir_whitebal"));
		if (properties.getProperty(prefix+"vnir_gain_g")!=null)         this.eo_gain_g=Double.parseDouble(properties.getProperty(prefix+"vnir_gain_g"));
		if (properties.getProperty(prefix+"vnir_gain_rg")!=null)        this.eo_gain_rg=Double.parseDouble(properties.getProperty(prefix+"vnir_gain_rg"));
		if (properties.getProperty(prefix+"vnir_gain_bg")!=null)        this.eo_gain_bg=Double.parseDouble(properties.getProperty(prefix+"vnir_gain_bg"));
		if (properties.getProperty(prefix+"vnir_exp_corr")!=null)       this.eo_exp_corr=str_to_darr(properties.getProperty(prefix+"vnir_exp_corr"));
		if (properties.getProperty(prefix+"vnir_gcorr_rbgb")!=null)     this.eo_gcorr_rbgb=str_to_darr(properties.getProperty(prefix+"vnir_gcorr_rbgb"));
		if (properties.getProperty(prefix+"vnir_lag")!=null)            this.eo_lag=Integer.parseInt(properties.getProperty(prefix+"vnir_lag")); // old version
		if (properties.getProperty(prefix+"avg_number")!=null)          this.avg_number=Integer.parseInt(properties.getProperty(prefix+"avg_number"));
		if (properties.getProperty(prefix+"avg_number_eo")!=null)       this.avg_number_eo=Integer.parseInt(properties.getProperty(prefix+"avg_number_eo"));
		if (properties.getProperty(prefix+"num_frames")!=null)          this.num_frames=Integer.parseInt(properties.getProperty(prefix+"num_frames"));
		if (properties.getProperty(prefix+"lwir_ffc")!=null)            this.lwir_ffc= Boolean.parseBoolean(properties.getProperty(prefix+"lwir_ffc"));
		if (properties.getProperty(prefix+"ffc_frames")!=null)          this.ffc_frames=Integer.parseInt(properties.getProperty(prefix+"ffc_frames"));
		if (properties.getProperty(prefix+"ffc_wait_frames")!=null)     this.ffc_wait_frames=Integer.parseInt(properties.getProperty(prefix+"ffc_wait_frames"));
		if (properties.getProperty(prefix+"ffc_groups")!=null)          this.ffc_groups=Integer.parseInt(properties.getProperty(prefix+"ffc_groups"));
		if (properties.getProperty(prefix+"acquire_delay")!=null)       this.acquire_delay=Double.parseDouble(properties.getProperty(prefix+"acquire_delay"));
		if (properties.getProperty(prefix+"avg_all")!=null)             this.avg_all= Boolean.parseBoolean(properties.getProperty(prefix+"avg_all"));
		if (properties.getProperty(prefix+"lwir_ip")!=null)             this.lwir_ips[0]=    properties.getProperty(prefix+"lwir_ip");
		if (properties.getProperty(prefix+"lwir_ips")!=null)            this.lwir_ips=  str_to_sarr(properties.getProperty(prefix+"lwir_ips"));
		if (properties.getProperty(prefix+"eo_ip")!=null)               this.eo_ip=      properties.getProperty(prefix+"eo_ip");
		if (properties.getProperty(prefix+"lwir_channels")!=null)       this.lwir_channels=str_to_iarr(properties.getProperty(prefix+"lwir_channels"));
		if (properties.getProperty(prefix+"eo_channels")!=null)         this.eo_channels=str_to_iarr(properties.getProperty(prefix+"eo_channels"));
		if (properties.getProperty(prefix+"lwir_telemetry")!=null)      this.lwir_telemetry= Boolean.parseBoolean(properties.getProperty(prefix+"lwir_telemetry"));
		if (properties.getProperty(prefix+"eo_full_window")!=null)      this.eo_full_window= Boolean.parseBoolean(properties.getProperty(prefix+"eo_full_window"));
		if (properties.getProperty(prefix+"eo_quality")!=null)          this.eo_quality=Double.parseDouble(properties.getProperty(prefix+"eo_quality"));
		if (properties.getProperty(prefix+"eo_scale")!=null)            this.eo_scale= Boolean.parseBoolean(properties.getProperty(prefix+"eo_scale"));
		if (properties.getProperty(prefix+"eo_autoexp")!=null)          this.eo_autoexp= Boolean.parseBoolean(properties.getProperty(prefix+"eo_autoexp"));
		if (properties.getProperty(prefix+"eo_max_autoexp_ms")!=null)   this.eo_max_autoexp_ms=Double.parseDouble(properties.getProperty(prefix+"eo_max_autoexp_ms"));
		if (properties.getProperty(prefix+"eo_exposure_ms")!=null)      this.eo_exposure_ms=Double.parseDouble(properties.getProperty(prefix+"eo_exposure_ms"));
		if (properties.getProperty(prefix+"eo_whitebal")!=null)         this.eo_whitebal= Boolean.parseBoolean(properties.getProperty(prefix+"eo_whitebal"));
		if (properties.getProperty(prefix+"eo_gain_g")!=null)           this.eo_gain_g=Double.parseDouble(properties.getProperty(prefix+"eo_gain_g"));
		if (properties.getProperty(prefix+"eo_gain_rg")!=null)          this.eo_gain_rg=Double.parseDouble(properties.getProperty(prefix+"eo_gain_rg"));
		if (properties.getProperty(prefix+"eo_gain_bg")!=null)          this.eo_gain_bg=Double.parseDouble(properties.getProperty(prefix+"eo_gain_bg"));
		if (properties.getProperty(prefix+"eo_exp_corr")!=null)         this.eo_exp_corr=str_to_darr(properties.getProperty(prefix+"eo_exp_corr"));
		if (properties.getProperty(prefix+"eo_gcorr_rbgb")!=null)       this.eo_gcorr_rbgb=str_to_darr(properties.getProperty(prefix+"eo_gcorr_rbgb"));
		if (properties.getProperty(prefix+"lwir_trig_dly")!=null)       this.lwir_trig_dly=Integer.parseInt(properties.getProperty(prefix+"lwir_trig_dly"));
		if (properties.getProperty(prefix+"eo_lag")!=null)              this.eo_lag=Integer.parseInt(properties.getProperty(prefix+"eo_lag"));
		if (properties.getProperty(prefix+"max_mismatch_ms")!=null)     this.max_mismatch_ms=Double.parseDouble(properties.getProperty(prefix+"max_mismatch_ms"));
		if (properties.getProperty(prefix+"max_frame_diff")!=null)      this.max_frame_diff=Integer.parseInt(properties.getProperty(prefix+"max_frame_diff"));
		if (properties.getProperty(prefix+"debug_level")!=null)         this.debug_level=Integer.parseInt(properties.getProperty(prefix+"debug_level"));
		if (properties.getProperty(prefix+"selected_channels")!=null)   this.selected_channels=str_to_barr(properties.getProperty(prefix+"selected_channels"));
		if (properties.getProperty(prefix+"show_images")!=null)         this.show_images= Boolean.parseBoolean(properties.getProperty(prefix+"show_images"));
		parameters_updated = true;
	}
	@Override
	public LwirReaderParameters clone() { //  throws CloneNotSupportedException {
		LwirReaderParameters lrp =  new LwirReaderParameters();
		lrp.camera_name =           this.camera_name;
		lrp.avg_number=             this.avg_number;
		lrp.avg_number_eo=          this.avg_number_eo;
		lrp.num_frames=             this.num_frames;
		lrp.lwir_ffc=               this.lwir_ffc;
		lrp.ffc_frames=             this.ffc_frames;
		lrp.ffc_wait_frames=        this.ffc_wait_frames;
		lrp.ffc_groups=             this.ffc_groups;
		lrp.acquire_delay=          this.acquire_delay;
		lrp.avg_all=                this.avg_all;
		lrp.lwir_ips=               this.lwir_ips.clone();
		lrp.eo_ip=                  this.eo_ip;
		lrp.lwir_channels=          this.lwir_channels.clone();
		lrp.eo_channels=            this.eo_channels.clone();
		lrp.lwir_telemetry=         this.lwir_telemetry;
		lrp.eo_full_window=         this.eo_full_window;
		lrp.eo_quality=             this.eo_quality;
		lrp.eo_scale=               this.eo_scale;
		lrp.eo_autoexp=             this.eo_autoexp;
		lrp.eo_max_autoexp_ms=      this.eo_max_autoexp_ms;
		lrp.eo_exposure_ms=         this.eo_exposure_ms;
		lrp.eo_whitebal=            this.eo_whitebal;
		lrp.eo_gain_g=              this.eo_gain_g;
		lrp.eo_gain_rg=             this.eo_gain_rg;
		lrp.eo_gain_bg=             this.eo_gain_bg;
		lrp.eo_exp_corr=            this.eo_exp_corr.clone();
		lrp.eo_gcorr_rbgb=          this.eo_gcorr_rbgb.clone();
		lrp.lwir_trig_dly=          this.lwir_trig_dly;
		lrp.eo_lag=                 this.eo_lag;
		lrp.max_mismatch_ms=        this.max_mismatch_ms;
		lrp.max_frame_diff=         this.max_frame_diff;
		lrp.debug_level=            this.debug_level;
		lrp.selected_channels =     this.selected_channels.clone();
		lrp.show_images =           this.show_images;
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
		boolean all_ips_eq = lrp.lwir_ips.length == this.lwir_ips.length;
		if (all_ips_eq) {
			for (int i = 0; (i < lrp.lwir_ips.length) && all_ips_eq; i++) {
				all_ips_eq &= (lrp.lwir_ips[i].equals(this.lwir_ips[i]));
			}
			
		}
		return  (lrp.camera_name.equals(this.camera_name)) &&
				(lrp.avg_number == this.avg_number) &&
				(lrp.avg_number_eo == this.avg_number_eo) &&
				(lrp.lwir_ffc == this.lwir_ffc) &&
				(lrp.ffc_frames == this.ffc_frames) &&
				(lrp.ffc_wait_frames == this.ffc_wait_frames) &&
				(lrp.ffc_groups == this.ffc_groups) &&
				(lrp.acquire_delay == this.acquire_delay) &&
				(lrp.avg_all == this.avg_all) &&
				all_ips_eq && //(lrp.lwir_ip.equals(this.lwir_ip)) &&
				(lrp.eo_ip.equals(this.eo_ip)) &&
				(java.util.Arrays.equals(lrp.lwir_channels, this.lwir_channels)) &&
				(java.util.Arrays.equals(lrp.eo_channels, this.eo_channels)) &&
				(lrp.lwir_telemetry == this.lwir_telemetry) &&
				(lrp.eo_full_window == this.eo_full_window) &&
				(lrp.eo_quality == this.eo_quality) &&
				(lrp.eo_scale == this.eo_scale) &&
				(lrp.eo_autoexp == this.eo_autoexp) &&
				(lrp.eo_max_autoexp_ms == this.eo_max_autoexp_ms) &&
				(lrp.eo_exposure_ms == this.eo_exposure_ms) &&
				(lrp.eo_whitebal == this.eo_whitebal) &&
				(lrp.eo_gain_g == this.eo_gain_g) &&
				(lrp.eo_gain_rg == this.eo_gain_rg) &&
				(lrp.eo_gain_bg == this.eo_gain_bg) &&
				(java.util.Arrays.equals(lrp.eo_exp_corr, this.eo_exp_corr)) &&
				(java.util.Arrays.equals(lrp.eo_gcorr_rbgb, this.eo_gcorr_rbgb)) &&
				(lrp.lwir_trig_dly == this.lwir_trig_dly) &&
				(lrp.eo_lag == this.eo_lag) &&
				(lrp.max_mismatch_ms == this.max_mismatch_ms) &&
				(lrp.max_frame_diff == this.max_frame_diff) &&
//				(lrp.debug_level ==    this.debug_level) &&
				(java.util.Arrays.equals(lrp.selected_channels, this.selected_channels));
//				&& (lrp.show_images ==    this.show_images);
	}

	@Override
	public int hashCode() {
		int prime = 31;
		int result = 1;
		result = prime * result + camera_name.hashCode();
		result = prime * result + ((Integer) avg_number).hashCode();
		result = prime * result + ((Integer) avg_number_eo).hashCode();
		result = prime * result + (lwir_ffc?1:0);
		result = prime * result + ((Integer) ffc_frames).hashCode();
		result = prime * result + ((Integer) ffc_wait_frames).hashCode();
		result = prime * result + ((Integer) ffc_groups).hashCode();
		result = prime * result + ((Double) acquire_delay).hashCode();
		result = prime * result + (avg_all?1:0);
		for (String ip:lwir_ips) {
			result = prime * result + ip.hashCode();
		}
		result = prime * result + eo_ip.hashCode();
		result = prime * result + arr_to_str(lwir_channels).hashCode();
		result = prime * result + arr_to_str(eo_channels).hashCode();
		result = prime * result + (lwir_telemetry?1:0);
		result = prime * result + (eo_full_window?1:0);
		result = prime * result + ((Double) eo_quality).hashCode();
		result = prime * result + (eo_scale?1:0);
		result = prime * result + (eo_autoexp?1:0);
		result = prime * result + ((Double) eo_max_autoexp_ms).hashCode();
		result = prime * result + ((Double) eo_exposure_ms).hashCode();
		result = prime * result + (eo_whitebal?1:0);
		result = prime * result + ((Double) eo_gain_g).hashCode();
		result = prime * result + ((Double) eo_gain_rg).hashCode();
		result = prime * result + ((Double) eo_gain_bg).hashCode();
		result = prime * result + arr_to_str(eo_exp_corr).hashCode();
		result = prime * result + arr_to_str(eo_gcorr_rbgb).hashCode();
		result = prime * result + ((Integer) lwir_trig_dly).hashCode();
		result = prime * result + arr_to_str(selected_channels).hashCode();
		// next are not needed to be programmed to the cameras
//		result = prime * result + (new Integer(eo_lag)).hashCode();
//		result = prime * result + (new Double(max_mismatch_ms)).hashCode();
//		result = prime * result + (new Integer(max_frame_diff)).hashCode();
		return 0;
	}

	public void dialogQuestions(GenericJTabbedDialog gd) {
		// TODO: separate parameters into common and type -specific
		gd. addChoice("Camera name", CAMERA_NAMES, camera_name, "Camera type - now Talon (4*Lepton+4*MT9P006) and LWIR16 (16xBoson+4*MT9P006");

		gd.addNumericField("Number to average LWIR",this.avg_number,      0,3,"","Number of acquired consecutive LWIR images to average");
		gd.addNumericField("Number to average EO",  this.avg_number_eo,   0,3,"","Number of acquired consecutive EO images to average");
		gd.addNumericField("Number of sets to acquire",this.num_frames,   0,3,"","Total number, may be limited by EO quality/frame size (last will be missed/null)");
		gd.addCheckbox    ("Run FFC",       this.lwir_ffc, "Perform calibration before each measurements to average (Lepton takes ~1.6 sec, 15 frames), Boson - twice ffc_frames");

		String [] sffc_frames = new String [BOSON_FFC_FRAMES.length];
		for (int i = 0; i < sffc_frames.length; i++) sffc_frames[i]=""+BOSON_FFC_FRAMES[i]; 
		gd. addChoice     ("FFC frames", sffc_frames, ""+ffc_frames,"Number of frames to integrate during FFC (Full FFC takes 2x this number of frames)" );
		gd.addNumericField("FFC wait frames",this.ffc_wait_frames,        0,3,"","Wait this number of frames after FFC is supposed to finish");
		
		String [] sffc_groups = new String [FFC_GROUPS.length];
		for (int i = 0; i < sffc_groups.length; i++) sffc_groups[i]=""+FFC_GROUPS[i];
		gd. addChoice     ("FFC groups", sffc_groups, ""+ffc_groups,"Split FFC actions in groups to limit power consumption" );
		gd.addNumericField("Acquire delay", this.acquire_delay,  3,6,"s", "Schedule synchronous acquisition ahead of current master camera time");
		gd.addCheckbox    ("Average all",   this.avg_all, "Average all simultaneously acquired images (unchecked - only requested number to average)");
//		gd.addStringField ("LWIR IP",       this.lwir_ip, 20, "LWIR camera IP address");
		for (int i = 0; i < this.lwir_ips.length; i++) {
			gd.addStringField ("LWIR IP"+i,       this.lwir_ips[i], 20, "LWIR camera "+(i+((i==0)?" (master)":""))+" IP address");
		}
		gd.addStringField ("EO IP",         this.eo_ip, 20, "Visible range high resolution camera IP address");
		gd.addStringField ("LWIR channels", arr_to_str(this.lwir_channels), 40, "Space-separated list of used LWIR camera channels, such as '0 1 2 3'");
		gd.addStringField ("EO channels", arr_to_str(this.eo_channels), 20, "Space-separated list of used visible range camera channels, such as '0 1 2 3'");
		gd.addCheckbox    ("LWIR telemetry",this.lwir_telemetry, "Set LWIR sesnors to provide telemetry data in the last 2 lines (may become mandatory later)");
		gd.addCheckbox    ("EO full window",this.eo_full_window, "Set EO to full window, disregarding preset margins (need reset to restore original)");
		gd.addNumericField("EO quality",    this.eo_quality,  3,6,"%", "Visible range camera JPEG compression quality (all channels)");
		gd.addCheckbox    ("EO undo white balance",    this.eo_scale, "Undo in-camera white balancing");
		gd.addCheckbox    ("EO autoexposure",this.eo_autoexp, "Enable autoexposure for the visible range camera");
		gd.addNumericField("EO eo_max_autoexp_ms", this.eo_max_autoexp_ms,  3,6,"ms", "Visible range camera maximal exposure in autoexposure mode");
		gd.addNumericField("EO exposure",   this.eo_exposure_ms,  3,6,"ms", "Visible range camera exposure time (all channels)");
		gd.addCheckbox    ("EO white balance", this.eo_whitebal, "Enable automatic white balancing for the visible range camera");
		gd.addNumericField("EO gain G",     this.eo_gain_g,  3,6,"","Analog gain for green channel for all visible range camera channels (normally 2.0)");
		gd.addNumericField("EO gain R/G",   this.eo_gain_rg, 3,6,"","Red to green gain ratio for all visible range camera channels");
		gd.addNumericField("EO gain B/G",   this.eo_gain_bg, 3,6,"","Blue to green gain ratio for all visible range camera channels");
		gd.addStringField ("EO exposure corrections", arr_to_str(this.eo_exp_corr), 50, "Fine corrections of channel exposures (4 channel relative exposures)");
		gd.addStringField ("EO gain corrections", arr_to_str(this.eo_gcorr_rbgb), 100, "Fine corrections to per channel, per color gains:'r0 b0 gb0 r1 b1 gb1 ...'");
		gd.addNumericField("LWIR trig dly", this.lwir_trig_dly,     0,10,"x10ns","Output trigger delay, should eventually match Lepton+FPGA latency to trigger EO exactly 1 frame after LWIR. 0 does not work with current FPGA - usec do not match sec in transmitted timestamp");
		gd.addNumericField("EO lag",        this.eo_lag,     0,3,"","Visible camera lag (in frames) relative to LWIR one");
		gd.addNumericField("Max mismatch",  this.max_mismatch_ms,    3,6,"ms","Maximal mismatch between image timestamps. Larger mismatch requires LWIR sinsor reinitialization");
		gd.addNumericField("Max frame diff",this.max_frame_diff,     0,3,"","Maximal difference in frames between simultaneously acquired channels as calculated from the timestamps");
		gd.addNumericField("Debug level",   this.debug_level,        0,3,"","Image acquisition log level: -3: OFF, -2:FATAL, -1:ERROR, 0:WARN, 1:INFO, 2:DEBUG");
		gd.addStringField ("Selected channels", arr_to_str(this.selected_channels), 40, "Space-separated channel selection (1 - selected, 0 - unselected)");
		gd.addCheckbox    ("Show images",   this.show_images, "Show acquired images after averaging)");
	}

	public void dialogAnswers(GenericJTabbedDialog gd) {
		String new_camera_name =       CAMERA_NAMES[gd.getNextChoiceIndex()];
		this.avg_number =        (int) gd.getNextNumber();
		this.avg_number_eo =     (int) gd.getNextNumber();
		this.num_frames =        (int) gd.getNextNumber();
		this.lwir_ffc =                gd.getNextBoolean();
		this.ffc_frames = BOSON_FFC_FRAMES[gd.getNextChoiceIndex()];
		
		this.ffc_wait_frames =   (int) gd.getNextNumber();
		this.ffc_groups =   FFC_GROUPS[gd.getNextChoiceIndex()];
		this.acquire_delay =           gd.getNextNumber();		
		this.avg_all =                 gd.getNextBoolean();
//		this.lwir_ip =                 gd.getNextString();
		for (int i = 0; i < this.lwir_ips.length; i++) {
			this.lwir_ips[i] =         gd.getNextString();
		}
		
		this.eo_ip =                   gd.getNextString();
		this.lwir_channels =           str_to_iarr(gd.getNextString());
		this.eo_channels =             str_to_iarr(gd.getNextString());
		this.lwir_telemetry =          gd.getNextBoolean();
		this.eo_full_window =          gd.getNextBoolean();
		this.eo_quality =              gd.getNextNumber();
		this.eo_scale =                gd.getNextBoolean();
		this.eo_autoexp =              gd.getNextBoolean();
		this.eo_max_autoexp_ms =       gd.getNextNumber();
		this.eo_exposure_ms =          gd.getNextNumber();
		this.eo_whitebal =             gd.getNextBoolean();
		this.eo_gain_g =               gd.getNextNumber();
		this.eo_gain_rg =              gd.getNextNumber();
		this.eo_gain_bg =              gd.getNextNumber();
		this.eo_exp_corr =             str_to_darr(gd.getNextString());
		this.eo_gcorr_rbgb =           str_to_darr(gd.getNextString());
		this.lwir_trig_dly =     (int) gd.getNextNumber();
		this.eo_lag =            (int) gd.getNextNumber();
		this.max_mismatch_ms =         gd.getNextNumber();
		this.max_frame_diff =    (int) gd.getNextNumber();
		this.debug_level =       (int) gd.getNextNumber();
		this.selected_channels =       str_to_barr(gd.getNextString());
		this.show_images =             gd.getNextBoolean();
		if (!this.camera_name.equals(new_camera_name)) {
			this.camera_name = new_camera_name; // in the very end to to disturb the number of items  
			setDefaultChannels(); // modifies lwir_IPs, eo_ip, channels,
		}
		parameters_updated = true;
	}


	public boolean showJDialog() {
		GenericJTabbedDialog gd = new GenericJTabbedDialog("Set LWIR/EO parameters",800,900);
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

	public boolean [] getSelectedEo(){
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

	private String arr_to_str(String [] arr) {
		String s = "";
		for (int i = 0; i < arr.length; i++) {
			s += arr[i].trim();
			if (i < (arr.length - 1)){
				s += ",";
			}
		}
		return s;
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

	private String [] str_to_sarr(String s) {
		String [] sa;
		if (s.indexOf(",") >= 0) {
			sa = s.trim().split(",");
		} else {
			sa = s.trim().split(" ");
		}
		String [] sarr = new String [sa.length];
		for (int i = 0; i < sa.length; i++) {
			sarr[i]=sa[i].trim();
		}
		return sarr;
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
