/**
 **
 ** LwocMesh.java - Single mesh of the world built from LWIR16 images
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwocMesh.java is free software: you can redistribute it and/or modify
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

import java.io.Serializable;
import java.util.concurrent.atomic.AtomicInteger;

public class LwocMesh implements Serializable {
	private static final long  serialVersionUID = 1L;
	static AtomicInteger       MESH_ID = new AtomicInteger();
//	static ArrayList<LwocMesh> LWOC_MESHES;
	int                        id;         // assign unique ID
	String                     stimestamp; // mesh is always referenced to a single scene
	double []                  world_xyz;  // world coordinates of the center
	double []                  hdims_xyz;  // world bounding box half-dimensions along x,y,z
	
	public static void resetMeshes() {
		MESH_ID.set(0);
//		LWOC_MESHES = new ArrayList<LwocMesh>();
	}
	// maybe use mesh properties instead of id?
	public boolean equals (LwocMesh other_mesh) {
		return id == other_mesh.id;
	}
	public LwocMesh(
			String    stimestamp
			) {
		this.stimestamp = stimestamp;
		id = MESH_ID.getAndIncrement();
//		LWOC_MESHES.add(this);
	}
	
	public double [] getCenter() {
		return world_xyz;
	}
	
	public double [] getDims() {
		return hdims_xyz;
	}
}
