package com.elphel.imagej.cameras;
/**
 **
 ** InterNoiseParameters - Class for handling configuration parameters
 ** to measure disparity map calculation with images degraded by the artificially
 ** introduced noise.
 ** related to the interscene LMA 
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  InterNoiseParameters.java is free software: you can redistribute it and/or modify
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
import com.elphel.imagej.tileprocessor.NoiseParameters;


public class InterNoiseParameters {
	public NoiseParameters noise = new NoiseParameters (0.1, 0.0, 1.5, 1.0, 16);
//	public double noise_sigma =        1.5;
//	public double noise_scale =        0.01;
//	public double initial_offset =     1.0;
	public boolean ref_only =          false; // also see imgdtt_params.dbg_pair_mask to switch between all pairs (63) and binocular only (1)
	public boolean show_final_2d =     true;  // show 2d correlations during last interation step
	public int     noise_debug_level = -1; // Noise testing debug level

	
	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addMessage("Additive noise parameters");
		gd.addMessage("LMA other parameters");
		gd.addNumericField("Noise Gaussian sigma",                            this.noise.sigma,          3,5,"pix",
				"Blur noise with 2D Gaussian");
		gd.addNumericField("Scale noise random (each scene indepemdent)",     this.noise.scale_random,   6,8,"",
				"Scale noise relative to the average value of the color component");
		gd.addNumericField("Scale noise FPN (same for each scene)",           this.noise.scale_fpn,      6,8,"",
				"Scale noise relative to the average value of the color component.");
		gd.addNumericField("Offset target disparity",                         this.noise.initial_offset, 3,5,"pix",
				"Offset target disparity before attempting to correlate and refine.");
		gd.addNumericField("Subset number of sensors (2/4/8/16)",             this.noise.used_sensors,   0,3,"",
				"Performance comparison - use only some of 16 sensors");
		
		gd.addCheckbox    ("Reference scene only",                            this.ref_only,
				"Process only reference scene (intra-scene, no inter-scene accumulation).");
		gd.addCheckbox    ("Show 2d correlations for the last iteration",     this.show_final_2d,
				"Show single-scene and interframe 2d correlations for the last iteration");
		gd.addNumericField("Debug level",                                     this.noise_debug_level,    0,3,"",
				"Debug level of interscene noise testing.");
		
	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.noise.sigma =                 gd.getNextNumber();
		this.noise.scale_random =          gd.getNextNumber();
		this.noise.scale_fpn =             gd.getNextNumber();
		this.noise.initial_offset =        gd.getNextNumber();
		this.noise.used_sensors =    (int) gd.getNextNumber();
		this.ref_only =                    gd.getNextBoolean();
		this.show_final_2d =               gd.getNextBoolean();
		this.noise_debug_level =     (int) gd.getNextNumber();
	}
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"noise.sigma",          this.noise.sigma+"");	
		properties.setProperty(prefix+"noise.scale_random",   this.noise.scale_random+"");	
		properties.setProperty(prefix+"noise.scale_fpn",      this.noise.scale_fpn+"");	
		properties.setProperty(prefix+"noise.initial_offset", this.noise.initial_offset+"");	
		properties.setProperty(prefix+"noise.used_sensors",   this.noise.used_sensors+"");	
		properties.setProperty(prefix+"ref_only",             this.ref_only+"");	
		properties.setProperty(prefix+"show_final_2d",        this.show_final_2d+"");	
		properties.setProperty(prefix+"noise_debug_level",    this.noise_debug_level+"");	
	}
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"noise.sigma")!=null) {
			this.noise.sigma=Double.parseDouble(properties.getProperty(prefix+"noise.sigma"));
		} else if (properties.getProperty(prefix+"noise_sigma")!=null) { // old format
			this.noise.sigma=Double.parseDouble(properties.getProperty(prefix+"noise_sigma"));
		}
		if (properties.getProperty(prefix+"noise.scale_random")!=null) {
			this.noise.scale_random=Double.parseDouble(properties.getProperty(prefix+"noise.scale_random"));
		} else if (properties.getProperty(prefix+"noise_scale")!=null) {
			this.noise.scale_random=Double.parseDouble(properties.getProperty(prefix+"noise_scale"));
		}
		if (properties.getProperty(prefix+"noise.scale_fpn")!=null) this.noise.scale_fpn=Double.parseDouble(properties.getProperty(prefix+"noise.scale_fpn"));
		
		if (properties.getProperty(prefix+"noise.initial_offset")!=null) {
			this.noise.initial_offset=Double.parseDouble(properties.getProperty(prefix+"noise.initial_offset"));
		} else if (properties.getProperty(prefix+"initial_offset")!=null) {
			this.noise.initial_offset=Double.parseDouble(properties.getProperty(prefix+"initial_offset"));
		}
		if (properties.getProperty(prefix+"noise.used_sensors")!=null) this.noise.used_sensors=Integer.parseInt(properties.getProperty(prefix+"noise.used_sensors"));

		if (properties.getProperty(prefix+"ref_only")!=null)          this.ref_only=Boolean.parseBoolean(properties.getProperty(prefix+"ref_only"));
		
		if (properties.getProperty(prefix+"show_final_2d")!=null)     this.show_final_2d=Boolean.parseBoolean(properties.getProperty(prefix+"show_final_2d"));
		
		if (properties.getProperty(prefix+"noise_debug_level")!=null) this.noise_debug_level=Integer.parseInt(properties.getProperty(prefix+"noise_debug_level"));
	}
	@Override
	public InterNoiseParameters clone() throws CloneNotSupportedException {
		InterNoiseParameters inp = new InterNoiseParameters();
		inp.noise =                this.noise.clone();
		inp.ref_only =             this.ref_only;
		inp.show_final_2d =        this.show_final_2d;
		inp.noise_debug_level =    this.noise_debug_level;
		return inp;
	}	
}
