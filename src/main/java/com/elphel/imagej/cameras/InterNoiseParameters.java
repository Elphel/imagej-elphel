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


public class InterNoiseParameters {
	public double noise_sigma =        1.5;
	public double noise_scale =        0.01;
	public double initial_offset =     1.0;
	public boolean ref_only =          false; // also see imgdtt_params.dbg_pair_mask to switch between all pairs (63) and binocular only (1)
	public int     noise_debug_level = -1; // Noise testing debug level

	
	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addMessage("Additive noise parameters");
		gd.addMessage("LMA other parameters");
		gd.addNumericField("Noise Gaussian sigma",                            this.noise_sigma, 3,5,"pix",
				"Blur noise with 2D Gaussian");
		gd.addNumericField("Scale noise",                                     this.noise_scale, 6,8,"",
				"Scale noise relative to the average value of the color component.");
		gd.addNumericField("Offset target disparity",                         this.initial_offset, 3,5,"pix",
				"Offset target disparity before attempting to correlate and refine.");
		gd.addCheckbox    ("Reference scene only",                            this.ref_only,
				"Process only reference scene (intra-scene, no inter-scene accumulation).");
		gd.addNumericField("Debug level",                                     this.noise_debug_level, 0,3,"",
				"Debug level of interscene noise testing.");
		
	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.noise_sigma =                 gd.getNextNumber();
		this.noise_scale =                 gd.getNextNumber();
		this.initial_offset =              gd.getNextNumber();
		this.ref_only =                    gd.getNextBoolean();
		this.noise_debug_level =     (int) gd.getNextNumber();
	}
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"noise_sigma",          this.noise_sigma+"");	
		properties.setProperty(prefix+"noise_scale",          this.noise_scale+"");	
		properties.setProperty(prefix+"initial_offset",       this.initial_offset+"");	
		properties.setProperty(prefix+"ref_only",       this.ref_only+"");	
		properties.setProperty(prefix+"noise_debug_level",       this.noise_debug_level+"");	
	}
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"noise_sigma")!=null)       this.noise_sigma=Double.parseDouble(properties.getProperty(prefix+"noise_sigma"));
		if (properties.getProperty(prefix+"noise_scale")!=null)       this.noise_scale=Double.parseDouble(properties.getProperty(prefix+"noise_scale"));
		if (properties.getProperty(prefix+"initial_offset")!=null)    this.initial_offset=Double.parseDouble(properties.getProperty(prefix+"initial_offset"));
		if (properties.getProperty(prefix+"ref_only")!=null)          this.ref_only=Boolean.parseBoolean(properties.getProperty(prefix+"ref_only"));
		if (properties.getProperty(prefix+"noise_debug_level")!=null) this.noise_debug_level=Integer.parseInt(properties.getProperty(prefix+"noise_debug_level"));
	}
	@Override
	public InterNoiseParameters clone() throws CloneNotSupportedException {
		InterNoiseParameters inp =     new InterNoiseParameters();
		inp.noise_sigma =          this.noise_sigma;
		inp.noise_scale =          this.noise_scale;
		inp.initial_offset =       this.initial_offset;
		inp.ref_only =             this.ref_only;
		inp.noise_debug_level =    this.noise_debug_level;
		return inp;
	}	
	
	
	
}
