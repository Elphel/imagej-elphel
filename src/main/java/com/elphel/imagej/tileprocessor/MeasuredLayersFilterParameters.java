package com.elphel.imagej.tileprocessor;
import java.util.Properties;

/**
 **
 ** MeasuredLayersFilterParameters - parameters defining tile filters
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ImageDtt.java is free software: you can redistribute it and/or modify
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

public class MeasuredLayersFilterParameters {
	public double     strength_sure  = 0.25;  // Do not filter disparity above this strength
	public double     strength_floor = 0.1;   // Subtract from strength, discard negative
	public double     strength_pow =   1.0;   // raise strength to this power
	public int        smplSide =       3;     // Sample size (side of a square)
	public int        smplNum =        5;     // Number after removing worst
	public double     smplRms =        0.3;   // Maximal RMS of the remaining tiles in a sample
	public boolean    smplWnd =        false; // Use window function for the samples (TODO: change default to true after testing)
	public double     max_abs_tilt =   2.0;   // Maximal absolute tilt in pixels/tile
	public double     max_rel_tilt =   0.2;   // Maximal relative tilt in pixels/tile/disparity
	public double     damp_tilt  =     0.001; // Tilt cost for damping insufficient plane data
	public double     min_tilt_disp =  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
	public double     transition =     1.0;   // Mode transition range (between tilted and maximal disparity)
	public int        far_mode  =      1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
	public double     far_power =      3.0;   // Raise disparity to this power before averaging for far objects

	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"strength_sure",  this.strength_sure +"");
		properties.setProperty(prefix+"strength_floor", this.strength_floor +"");
		properties.setProperty(prefix+"strength_pow",   this.strength_pow +"");

		properties.setProperty(prefix+"smplSide",       this.smplSide+"");
		properties.setProperty(prefix+"smplNum",        this.smplNum+"");
		properties.setProperty(prefix+"smplRms",        this.smplRms +"");
		properties.setProperty(prefix+"smplWnd",        this.smplWnd+"");

		properties.setProperty(prefix+"max_abs_tilt",   this.max_abs_tilt+"");
		properties.setProperty(prefix+"max_rel_tilt",   this.max_rel_tilt+"");
		properties.setProperty(prefix+"damp_tilt",      this.damp_tilt+"");
		properties.setProperty(prefix+"min_tilt_disp",  this.min_tilt_disp+"");
		properties.setProperty(prefix+"transition",     this.transition+"");
		properties.setProperty(prefix+"far_mode",       this.far_mode+"");
		properties.setProperty(prefix+"far_power",      this.far_power+"");
	}

	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"strength_sure")!=null)  this.strength_sure=Double.parseDouble(properties.getProperty(prefix+"strength_sure"));
		if (properties.getProperty(prefix+"strength_floor")!=null) this.strength_floor=Double.parseDouble(properties.getProperty(prefix+"strength_floor"));
		if (properties.getProperty(prefix+"strength_pow")!=null)   this.strength_pow=Double.parseDouble(properties.getProperty(prefix+"strength_pow"));

		if (properties.getProperty(prefix+"smplSide")!=null)       this.smplSide=Integer.parseInt(properties.getProperty(prefix+"smplSide"));
		if (properties.getProperty(prefix+"smplNum")!=null)        this.smplNum=Integer.parseInt(properties.getProperty(prefix+"smplNum"));
		if (properties.getProperty(prefix+"smplRms")!=null)        this.smplRms=Double.parseDouble(properties.getProperty(prefix+"smplRms"));
		if (properties.getProperty(prefix+"smplWnd")!=null)        this.smplWnd=Boolean.parseBoolean(properties.getProperty(prefix+"smplWnd"));


		if (properties.getProperty(prefix+"max_abs_tilt")!=null)   this.max_abs_tilt=Double.parseDouble(properties.getProperty(prefix+"max_abs_tilt"));
		if (properties.getProperty(prefix+"max_rel_tilt")!=null)   this.max_rel_tilt=Double.parseDouble(properties.getProperty(prefix+"max_rel_tilt"));
		if (properties.getProperty(prefix+"damp_tilt")!=null)      this.damp_tilt=Double.parseDouble(properties.getProperty(prefix+"damp_tilt"));

		if (properties.getProperty(prefix+"min_tilt_disp")!=null)  this.min_tilt_disp=Double.parseDouble(properties.getProperty(prefix+"min_tilt_disp"));
		if (properties.getProperty(prefix+"transition")!=null)     this.transition=Double.parseDouble(properties.getProperty(prefix+"transition"));
		if (properties.getProperty(prefix+"far_mode")!=null)       this.far_mode=Integer.parseInt(properties.getProperty(prefix+"far_mode"));
		if (properties.getProperty(prefix+"far_power")!=null)      this.far_power=Double.parseDouble(properties.getProperty(prefix+"far_power"));
	}

	@Override
	public MeasuredLayersFilterParameters clone() {
		MeasuredLayersFilterParameters mlfp = new MeasuredLayersFilterParameters();
		mlfp.strength_sure =  this.strength_sure;
		mlfp.strength_floor = this.strength_floor;
		mlfp.strength_pow =   this.strength_pow;
		mlfp.smplSide =       this.smplSide;
		mlfp.smplNum =        this.smplNum;
		mlfp.smplRms =        this.smplRms;
		mlfp.smplWnd =        this.smplWnd;
		mlfp.max_abs_tilt =   this.max_abs_tilt;
		mlfp.max_rel_tilt =   this.max_rel_tilt;
		mlfp.damp_tilt  =     this.damp_tilt;
		mlfp.min_tilt_disp =  this.min_tilt_disp;
		mlfp.transition =     this.transition;
		mlfp.far_mode  =      this.far_mode;
		mlfp.far_power =      this.far_power;
		return mlfp;

	}

	public boolean equals(MeasuredLayersFilterParameters mlfp) {
		return
		(mlfp.strength_sure ==  this.strength_sure ) &&
		(mlfp.strength_floor == this.strength_floor ) &&
		(mlfp.strength_pow ==   this.strength_pow ) &&
		(mlfp.smplSide ==       this.smplSide ) &&
		(mlfp.smplNum ==        this.smplNum ) &&
		(mlfp.smplRms ==        this.smplRms ) &&
		(mlfp.smplWnd ==        this.smplWnd ) &&
		(mlfp.max_abs_tilt ==   this.max_abs_tilt ) &&
		(mlfp.max_rel_tilt ==   this.max_rel_tilt ) &&
		(mlfp.damp_tilt  ==     this.damp_tilt ) &&
		(mlfp.min_tilt_disp ==  this.min_tilt_disp ) &&
		(mlfp.transition ==     this.transition ) &&
		(mlfp.far_mode  ==      this.far_mode ) &&
		(mlfp.far_power ==      this.far_power);
	}

}
