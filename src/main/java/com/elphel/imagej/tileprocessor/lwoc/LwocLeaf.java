/**
 **
 ** LwocLeaf.java - Octree leaf of the world built from LWIR16 images
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwocLeaf.java is free software: you can redistribute it and/or modify
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

public class LwocLeaf {
	ArrayList <LwocMesh>  meshes;
	ArrayList <LwocMesh>  mesh_centers;
	ArrayList <LwocScene> scenes;
	public LwocLeaf() {
		meshes =       new ArrayList <LwocMesh>();  // all meshes BB intersecting this node
		mesh_centers = new ArrayList <LwocMesh>();  // mesh centers in this node
		scenes =       new ArrayList <LwocScene>(); // cameras located in this node
	}
	
	public void addScene(LwocScene scene,
			boolean check_existed) {
		if (!check_existed || !scenes.contains(scene)) {
			scenes.add(scene);
		}
	}
	
	public void addMeshCenter(LwocMesh mesh,
			boolean check_existed) {
		if (!check_existed || !mesh_centers.contains(mesh)) {
			mesh_centers.add(mesh);
		}
	}

	public void addMesh(LwocMesh mesh,
			boolean check_existed) {
		if (!check_existed || !meshes.contains(mesh)) {
			meshes.add(mesh);
		}
	}
	
	public ArrayList <LwocMesh> getMeshes(){
		return meshes;
	}
	
	public ArrayList <LwocMesh> getMeshCenters(){
		return mesh_centers;
	}
	
	public ArrayList <LwocScene> getScenes(){
		return scenes;
	}
}
