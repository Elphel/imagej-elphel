/**
 **
 ** LWIRWorldParameters - Class for handling multiple configuration parameters
 ** related to the octree-based world(s) parameters  
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LWIRWorldParameters.java is free software: you can redistribute it and/or modify
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

package com.elphel.imagej.tileprocessor;

import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;

public class LWIRWorldParameters {
	public double    world_hsize =    1000.0; // m			
	public double    min_hsize =         0.3; // m
	public double    max_hsize =     10000.0; // m			
	public int       max_mesh_centers = 10; //
	public int       max_cameras      = 10; //
	
	public LWIRWorldParameters() {
		
	}

	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addTab      ("LWIR World", "Building 3D World from multiple scenes/ scenes series");
		gd.addMessage  ("Octree world parameters");
		gd.addNumericField("World half-size",                          this.world_hsize, 1,7,"m",
				"Initial half-size of the world (it may grow if needed).");
		gd.addNumericField("Minimal node half-size",                   this.min_hsize, 1,7,"m",
				"Do not subdivide nodes if its half-size is smaller than this.");
		gd.addNumericField("Maximal node half-size",                   this.max_hsize, 1,7,"m",
				"Do not grow the root node (whole world) if its half-size is >= this.");
		gd.addNumericField("Maximal number of mesh centers",           this.max_mesh_centers, 0,3,"",
				"Subdivide octree node that has more than this number of mesh centers.");
		gd.addNumericField("Maximal number of scenes/cameras",         this.max_cameras, 0,3,"",
				"Subdivide octree node that has more than this number of scenes/cameras.");
	}	
	
	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.world_hsize =                  gd.getNextNumber();
		this.min_hsize =                    gd.getNextNumber();
		this.max_hsize =                    gd.getNextNumber();
		this.max_mesh_centers =       (int) gd.getNextNumber();
		this.max_cameras =            (int) gd.getNextNumber();
	}
	
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"world_hsize",                     this.world_hsize+"");                   // double
		properties.setProperty(prefix+"min_hsize",                       this.min_hsize+"");                     // double
		properties.setProperty(prefix+"max_hsize",                       this.max_hsize+"");                     // double
		properties.setProperty(prefix+"max_mesh_centers",                this.max_mesh_centers+"");              // int
		properties.setProperty(prefix+"max_cameras",                     this.max_cameras+"");                   // int
	}
	
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"world_hsize")!=null)          this.world_hsize=Double.parseDouble(properties.getProperty(prefix+"world_hsize"));
		if (properties.getProperty(prefix+"min_hsize")!=null)            this.min_hsize=Double.parseDouble(properties.getProperty(prefix+"min_hsize"));
		if (properties.getProperty(prefix+"max_hsize")!=null)            this.max_hsize=Double.parseDouble(properties.getProperty(prefix+"max_hsize"));
		if (properties.getProperty(prefix+"max_mesh_centers")!=null)     this.max_mesh_centers=Integer.parseInt(properties.getProperty(prefix+"max_mesh_centers"));
		if (properties.getProperty(prefix+"max_cameras")!=null)          this.max_cameras=Integer.parseInt(properties.getProperty(prefix+"max_cameras"));
	}
	
	@Override
	public LWIRWorldParameters clone() throws CloneNotSupportedException {
		LWIRWorldParameters lwp =     new LWIRWorldParameters();
		lwp.world_hsize =               this.world_hsize;
		lwp.min_hsize =                 this.min_hsize;
		lwp.max_hsize =                 this.max_hsize;
		lwp.max_mesh_centers =          this.max_mesh_centers;
		lwp.max_cameras =               this.max_cameras;
		return lwp;
	}	
}
