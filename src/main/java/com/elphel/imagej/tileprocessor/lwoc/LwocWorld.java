/**
 **
 ** LwocWorld.java - Octree world parameters (previous static)
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwocWorld.java is free software: you can redistribute it and/or modify
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

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class LwocWorld  implements Serializable {
	private static final long   serialVersionUID = 1L;
	static List<LwocWorld>      lwoc_worlds;   //  
	
	public double               min_hsize; //  =         0.3; // meter - do not subdivide more
	public double               max_hsize; //  =     10000.0; // meter - do not grow more
	public int                  max_mesh_centers; //  = 10;   // maximal number of meshes in a leaf;
	public int                  max_cameras; //  =      10;   // maximal number of cameras in a leaf;
//	public transient LwocOctree lwoc_root; //  =        null;
	public LwocOctree           lwoc_root; //  =        null;
	public double []            atr; //  =              null; // Azimuth, tilt, roll - may be used to merge
	public double []            xyz; //  =              null; // this world offset (before rotation) or
    // null if not known.
	public LwocWorld (
			double     min_hsize,
			double     max_hsize,
			int        max_mesh_centers,
			int        max_cameras,
			LwocOctree lwoc_root,
			double []  atr,
			double []  xyz) {
		setMinHsize      (min_hsize);
		setMaxHsize      (max_hsize);
		setMaxMeshCenters(max_mesh_centers);
		setMaxCameras    (max_cameras);
		setLwocRoot      (lwoc_root);
		setATR           (atr);
		setXYZ           (xyz);
	}
	public double     getMinHsize()       {return min_hsize;}
	public double     getMaxHsize()       {return max_hsize;}
	public int        getMaxMeshCenters() {return max_mesh_centers;}
	public int        getMaxCameras()     {return max_cameras;}
	public LwocOctree getLwocRoot()       {return lwoc_root;}
	public double []  getATR()            {return atr;}
	public double []  getXYZ()            {return xyz;}
	
	public void       setMinHsize      (double     min_hsize)        {this.min_hsize =        min_hsize;}
	public void       setMaxHsize      (double     max_hsize)        {this.max_hsize =        max_hsize;}
	public void       setMaxMeshCenters(int        max_mesh_centers) {this.max_mesh_centers = max_mesh_centers;}
	public void       setMaxCameras    (int        max_cameras)      {this.max_cameras =      max_cameras;}
	public void       setLwocRoot      (LwocOctree lwoc_root)        {this.lwoc_root =        lwoc_root;}
	public void       setATR           (double []  atr)              {this.atr =              atr;}
	public void       setXYZ           (double []  xyz)              {this.xyz =              xyz;}
	
	
	/**
	 * Serialize and write world to an output stream
	 * @param oos ObjectOutputStream to write to.
	 * @throws IOException
	 */
	private void writeObject(ObjectOutputStream oos) 
			throws IOException {
		oos.defaultWriteObject();
		// as of now - nothing else, custom for reads only
		
//		if (parent == null) { // write world members
			/*
			oos.writeObject(world.getMinHsize());
			oos.writeObject(world.getMaxHsize());
			oos.writeObject(world.getMaxMeshCenters());
			oos.writeObject(world.getMaxCameras());
			oos.writeObject(world.getATR());
			oos.writeObject(world.getXYZ());
			*/
//			oos.writeObject(world);
			
			// Remove those
//			oos.writeObject(OCTREE_ID.get());
//			oos.writeObject(LwocMesh.MESH_ID.get());
			
//		}
	}

	/**
	 * Read input stream and deserialize it to world octree structure.
	 * @param ois ObjectInputStream to read from
	 * @throws ClassNotFoundException
	 * @throws IOException
	 */
	private void readObject(ObjectInputStream ois) 
			throws ClassNotFoundException, IOException {
		ois.defaultReadObject();
		// now rebuild transient .parent and .world field for the whole tree 
		lwoc_root.restoreWorldParent(
	    		this, // LwocWorld world,
	    		null); //LwocOctree parent)
		lwoc_root.rebuildMeshLists(
	    		true); // final boolean check_existed)
	}
	
	/**
	 * Init worlds and their list. Multiple worlds may be needed if some world initially can not
	 * be matched. They will grow independently and possibly merged later
	 */
	public static void initWorlds() {
		LwocMesh.resetMeshes();
		LwocScene.resetScenes();
		LwocOctree.resetPending();
		lwoc_worlds = new ArrayList<LwocWorld>();
	}
	
	/**
	 * Create a new empty world and add it to the list of worlds in the universe.
	 * @param center X,Y,Z of the new world 
	 * @param world_hsize new world half-size
	 * @param min_hsize do not subdivide nodes if half-size is smaller than this
	 * @param max_hsize do not grow world if its half-size it >= this
	 * @param max_mesh_centers maximal number of mesh centers in a node to trigger subdivision
	 * @param max_cameras maximal number of scenes (cameras) in a node to trigger subdivision
	 * @param atr world rotation relative to the universe (non-null when defined)
	 * @param xyz world offset in the universe (non-null when defined)
	 * @return new world instance
	 */
	public static LwocWorld newWorld(
			double [] center,
			double    world_hsize,			
			double    min_hsize,
			double    max_hsize,			
			int       max_mesh_centers,
			int       max_cameras,
			double [] atr,
			double [] xyz) {
		LwocWorld new_world = new  LwocWorld (
				min_hsize,        // double     min_hsize,
				max_hsize,        // double     max_hsize,
				max_mesh_centers, // int        max_mesh_centers,
				max_cameras,      // int        max_cameras,
				null,             // LwocOctree lwoc_root,
				null,             //double []  atr,
				null);            //double []  xyz);
		LwocOctree root_node = new LwocOctree (
				new_world, // LwocWorld  lwocWorld,
				null,      // LwocOctree parent,
				center,    // double []  xyz,
				world_hsize); //double     hsize)
		root_node.leaf = new LwocLeaf();
		new_world.setLwocRoot(root_node);
		lwoc_worlds.add(new_world);
		return new_world;
	}

	/**
	 * Delete world from the universe (if existed)
	 * @param world world instance to remove
	 * @return True if removed, false if did not exist. 
	 */
	public static boolean deleteWorld (LwocWorld world) {
		return lwoc_worlds.remove(world);
	}

	
	
}
