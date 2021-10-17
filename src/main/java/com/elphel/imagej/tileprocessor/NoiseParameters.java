package com.elphel.imagej.tileprocessor;

import com.elphel.imagej.cameras.InterNoiseParameters;

public class NoiseParameters {
	public double scale_random =   0.1;
	public double scale_fpn    =   0.0;
	public double sigma        =   1.5;
	public double initial_offset = 1.0;
	public int    used_sensors =   16; // 2/4/8/16
	
	public NoiseParameters (
			double scale_random,
			double scale_fpn,
			double sigma,
			double initial_offset,
			int used_sensors) {
		this.scale_random = scale_random;
		this.scale_fpn = scale_fpn;
		this.sigma = sigma;
		this.initial_offset= initial_offset;
		this.used_sensors = used_sensors;
	}
	@Override
	public NoiseParameters clone() throws CloneNotSupportedException {
		return  new NoiseParameters(
				scale_random,
				scale_fpn,
				sigma,
				initial_offset,
				used_sensors);
	}	
	
}
