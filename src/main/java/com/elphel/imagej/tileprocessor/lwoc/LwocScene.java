/**
 **
 ** LwocScene.java - Scene of the world built from LWIR16 images
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwocScene.java is free software: you can redistribute it and/or modify
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

package com.elphel.imagej.tileprocessor.lwoc;

import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import com.elphel.imagej.tileprocessor.GeometryCorrection;

public class LwocScene {
	private static final long    serialVersionUID = 1L;
	static AtomicInteger SCENE_ID = new AtomicInteger();
	static ArrayList<LwocScene> LWOC_SCENES;
	
	int       id;         // assign unique ID
	String    stimestamp;
	GeometryCorrection geometryCorrection;
	double [] camera_xyz;
	double [] camera_atr;
	double [] camera_xyz_dt;
	double [] camera_atr_dt;
	double [][][] tile_layers_dsm; // per tile, per layer, {disparity, strength, mode}. For mode -
	// interpret as Double doubleToLongBits(), [22] - used LMA, 0..15 - which sensors were used

	public static void resetScenes() {
		SCENE_ID.set(0);
		LWOC_SCENES = new ArrayList<LwocScene>();
	}
	
	public LwocScene (
		String             stimestamp,
		GeometryCorrection geometryCorrection,
		double []          camera_xyz,
		double []          camera_atr,
		double []          camera_xyz_dt,
		double []          camera_atr_dt,
		double [][][]      tile_layers_dsm) {
		
		this.stimestamp =         stimestamp;
		this.geometryCorrection = geometryCorrection;
		this.camera_xyz =         camera_xyz;
		this.camera_atr =         camera_atr;
		this.camera_xyz_dt =      camera_xyz_dt;
		this.camera_atr_dt =      camera_atr_dt;
		this.tile_layers_dsm =    tile_layers_dsm;
		id = SCENE_ID.getAndIncrement();
		LWOC_SCENES.add(this);
	};
	// add functionality to save/restore
	public double [] getCameraXYZ() {
		return camera_xyz;
	}
	public double [] getCameraATR() {
		return camera_atr;
	}
	public double [] getCameraXYZdt() {
		return camera_xyz_dt;
	}
	public double [] getCameraATRdt() {
		return camera_atr_dt;
	}
	public double [][][] getTileLayersDSM() {
		return tile_layers_dsm;
	}

	
	
	public String getTimestamp() {
    	return stimestamp;
    }
	
	public static String getTimestamp(double dts) {
    	return String.format("%.6f", dts).replace('.', '-');
    }

	public double getDoubleTimestamp() {
    	return Double.parseDouble(stimestamp.replace('_', '.'));
    }
	
	public static double getDoubleTimestamp(String sts) {
    	return Double.parseDouble(sts.replace('_', '.'));
    }
	
	public int getIntTimestamp() {
    	return (int) Math.round(1E6*getDoubleTimestamp());
    }
	public void setTimeStamp(String ts) {
    	stimestamp =  ts;
    }
	public void setTimeStamp(double dts) {
    	stimestamp = getTimestamp(dts);
    }
	public void setTimeStamp(int its) {
    	stimestamp =  String.format("%d_%06d", its / 1000000, its % 1000000);
    }
}
