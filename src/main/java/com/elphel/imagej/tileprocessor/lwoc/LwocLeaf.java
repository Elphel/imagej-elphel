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

import java.io.Serializable;
import java.util.ArrayList;

public class LwocLeaf implements Serializable{
	private static final long    serialVersionUID = 1L;
	transient ArrayList <LwocMesh>  meshes; // will be rebuilt with LwocOctree.rebuildMeshLists()
	ArrayList <LwocMesh>  mesh_centers;
	ArrayList <LwocScene> scenes;
	public LwocLeaf() {
		meshes =       new ArrayList <LwocMesh>();  // all meshes BB intersecting this node
		mesh_centers = new ArrayList <LwocMesh>();  // mesh centers in this node
		scenes =       new ArrayList <LwocScene>(); // cameras located in this node
	}
	
	public void initMeshes() {
		meshes =       new ArrayList <LwocMesh>();  // all meshes BB intersecting this node
	}
	public synchronized void addScene( // IS synchronized needed?
			LwocScene scene,
			boolean check_existed) {
		if (!check_existed || !scenes.contains(scene)) {
			scenes.add(scene);
		}
	}
	
	public synchronized void addMeshCenter( // IS synchronized needed?
			LwocMesh mesh,
			boolean check_existed) {
		if (!check_existed || !mesh_centers.contains(mesh)) {
			mesh_centers.add(mesh);
		}
	}

	public synchronized void removeScene(
			LwocScene scene) {
		while (scenes.remove(scene)); // will remove all
	}

	public synchronized void removeMeshCenter(
			LwocMesh mesh) {
		while (mesh_centers.remove(mesh)); // will remove all
	}

	public synchronized void removeMesh(
			LwocMesh mesh) {
		while (meshes.remove(mesh)); // will remove all
	}
	
	
	public synchronized void addMesh( // synchronized IS needed
			LwocMesh mesh,
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
