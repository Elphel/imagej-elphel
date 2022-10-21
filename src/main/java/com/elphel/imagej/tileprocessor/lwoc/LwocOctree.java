/**
 **
 ** LwocOctree.java - Octree of the world built from LWIR16 images
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwocOctree.java is free software: you can redistribute it and/or modify
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
import com.elphel.imagej.common.MultiThreading;
public class LwocOctree {
	public static final int      NUM_CHILDREN =      8;
	public static double         MIN_HSIZE =         0.3; // meter - do not subdivide more
	// can not be used as intersecting meshes' bb will not decrease number after splitting
	public static int            MAX_MESH_CENTERS = 10;   // maximal number of meshes in a leaf;
	public static int            MAX_CAMERAS =      10;   // maximal number of cameras in a leaf;
	public static LwocOctree     lwoc_root =        null;
	static ArrayList<LwocOctree> LWOK_OCTREE =      new ArrayList<LwocOctree>();
	static AtomicInteger         OCTREE_ID =        new AtomicInteger();
	// synchronized lists to handle thread-unsafe operations: growing world and splitting leaves
///	static List<LwocScene>       pendingScenes;   // synchronized - add not-yet-added scenes pending growing world 
///	static List<LwocMesh>        pendingMeshes;   // synchronized - add not-yet-added meshes pending growing world
///	static List<LwocOctree>      pendingLeafNodes;// synchronized - leaf nodes needed to be split
	
	static AtomicInteger         numPendingScenes;
	static AtomicInteger         numPendingMeshes;
	static AtomicInteger         numPendingLeafNodes;
	
	
	int                          id;         // assign unique ID
	LwocOctree                   parent;
	LwocOctree []                children;
	LwocLeaf                     leaf;       //
	double []                    center;     // x-center, y-center, z-center
	double                       hsize;
// just for testing	
//	List<List<LwocScene>> group = new ArrayList<List<LwocScene>>(4);	
	
	/**
	 * Initialize Octree data and set initial world node. It is possible to grow
	 * world as needed later.
	 * @param center world center: x,y,z in meters. Usually {0,0,0}
	 * @param world_hsize Half linear dimension of the world
	 */
	public static void initOctree(double [] center, double world_hsize) {
		LwocMesh.resetMeshes();
		LwocScene.resetScenes();
		resetPending();
		LWOK_OCTREE =    new ArrayList<LwocOctree>();
		// Add a single leaf node
		lwoc_root =      new LwocOctree(null, center, world_hsize);
		lwoc_root.leaf = new LwocLeaf();
	}
	
	/**
	 * Synchronized lists of pending scenes, meshes and nodes are used to delay modification of
	 * the world octree and its elements from the multithreaded environment. These lists should
	 * be tested after multithreaded methods that potentially can modify the world and then
	 * processes by a single-thread method. 
	 */
	public static void resetPending() {
		numPendingScenes.set(0);    // = Collections.synchronizedList(new ArrayList<LwocScene>());
		numPendingMeshes.set(0);    // = Collections.synchronizedList(new ArrayList<LwocMesh>());
		numPendingLeafNodes.set(0); // = Collections.synchronizedList(new ArrayList<LwocOctree>());
	}
	// use per-thread lists and then combine them avoiding synchronizations?
	//Pass list to any function that may grow it
	
	/**
	 * Test if it is leaf node
	 * @return True if it is a leaf node
	 */
	public boolean isLeaf() {
		return children == null;
	}
	
	/**
	 * Check if the world needs growing to accommodate new scenes or meshes 
	 * @return if world needs growing.
	 */
	public static boolean pendingGrow() {
//		return !pendingScenes.isEmpty() || !pendingMeshes.isEmpty();
		return (numPendingScenes.get()  + numPendingMeshes.get()) > 0;
	}
	
	/**
	 * Check if any leaf nodes need splitting
	 * @return if any leaf nodes require splitting
	 */
	public static boolean pendingSplit() {
		return numPendingLeafNodes.get() > 0;
	}
	
	/**
	 * LwocOctree constructor
	 * @param parent Parent node
	 * @param xyz position of the node center in meters
	 * @param hsize Node half size - each of the X,Y,Z coordinates of the internal 
	 *              points are limited within +/-hsize from the node center
	 */
    public LwocOctree (LwocOctree parent,
    		           double []  xyz,
    		           double     hsize) {
    	id = OCTREE_ID.getAndIncrement();
    	this.parent = parent;
    	this.hsize =  hsize;
    	LWOK_OCTREE.add(this);
    }
    
    /**
     * Find the leaf node that includes specified point. Thread safe.
     * @param xyz X, Y, and Z coordinates of the point
     * @return Octree node that includes the point or null if the point
     *         is outside the root node (the whole world)
     */
    public static LwocOctree getLeafNode( // if null - needs growing world
    		double [] xyz) { // thread safe
    	if (!lwoc_root.contains(xyz)) return null; // needs growing world
    	LwocOctree node = lwoc_root;
    	while (lwoc_root.children != null) {
    		int indx = 0;
    		for (int dm = 0; dm < 3; dm++) {
    			if (xyz[dm] >= lwoc_root.center[dm] ) {
    				indx += 1 << dm;
    			}
    		}
    		node = node.children[indx];
    	}
		return node;
    }
    
    /**
     * Add scene to the world in a thread-safe mode, or put it to pending list if the world
     * needs growing.
     * @param scene Scene to add.
     * @param pendingScenes per thread list of scenes to be combined after merging threads
     * @param pendingLeafNodes  per thread list of nodes to be combined after merging threads
     * @return An octree node where scene is added or null if the world needs growing
     */
    public static LwocOctree addScene(
    		LwocScene             scene,
    		ArrayList<LwocScene>  pendingScenes,
    		ArrayList<LwocOctree> pendingLeafNodes
    		) {// returns null and adds to pendingScenes if needs growing
    	LwocOctree node = getLeafNode(scene.getCameraXYZ());
    	if (node == null) {
    		if (pendingScenes != null) {
    			pendingScenes.add(scene);
    		}
    	} else {
    		node.leaf.addScene(scene);
    		// Check if leaf node requires splitting
    		if (node.leaf.getScenes().size() > MAX_CAMERAS) {
    			if (node.hsize > MIN_HSIZE) {
    				pendingLeafNodes.add(node);
    			}
    		}
    	}
    	return node;
    }
    
    /**
     * Tend to pending scenes list. Not thread-safe, should be run in a single-thread mode.
     * @param pendingScenes    combined list of scenes to be added growing world (merged from per-thread ones)
     * @param pendingLeafNodes a list of nodes that will need splitting (here in single-thread mode - single one)
     */
    public static void tendPendingScenes(
    		ArrayList<LwocScene>  pendingScenes, 
    		ArrayList<LwocOctree> pendingLeafNodes
    		) {
    	if (pendingScenes.isEmpty()) {
    		return; // nothing to do
      	} else {
    		for (LwocScene scene: pendingScenes) {
    			growXYZ(scene.getCameraXYZ());
    			LwocOctree node= addScene(  // should not be null
    					scene,
    		    		pendingScenes,
    		    		pendingLeafNodes);
    			assert node != null : "addScene() should not fail after growing world";
    		} 
    	}
    }
    
    /**
     * Tend to pending meshes list. Not thread-safe, should be run in a single-thread mode.
     * @param check_existed    do not add duplicate meshes
     * @param pendingMeshes    combined list of meshes to be added (merged from per-thread ones). Will only be read.
     * @param pendingLeafNodes a list of nodes that will need splitting (here in single-thread mode - single one). May grow.
     */
    public static void tendPendingMeshCenters(
    		boolean               check_existed,
    		ArrayList<LwocMesh>   pendingMeshes,   // will read
    		ArrayList<LwocOctree> pendingLeafNodes // may add to
    		) {
    	if (pendingMeshes.isEmpty()) {
    		return; // nothing to do
      	} else {
    		for (LwocMesh mesh: pendingMeshes) {
    	    	double [] center = mesh.getCenter();
    	    	double [] dims =   mesh.getDims();
    			double [] corner_xyz = center.clone();
    	    	for (int dm = 0; dm < center.length; dm++) {
    	    		corner_xyz[dm] += dims[dm];
    	    		growXYZ(corner_xyz);
    	    		corner_xyz[dm] -= 2*dims[dm];
    	    		growXYZ(corner_xyz);
    	    		corner_xyz[dm] += dims[dm];
    	    	}    			
    			boolean ok= addMeshCenter(  // should not be null
    					mesh,
    					check_existed,
    					null, // pendingMeshes,
    		    		pendingLeafNodes);
    			assert ok : "addMeshCenter() should not fail after growing world";
//    			int num_added = 
    			lwoc_root.addMesh(
    					mesh,
    					check_existed);
    		} 
    	}
    }
    
    public static void tendPendingLeafNodes(
    		ArrayList<LwocOctree> pendingLeafNodesIn) {
    	// remove any possible duplicates
    	final ArrayList<LwocOctree> pendingLeafNodes = new ArrayList<LwocOctree>();
    	for (LwocOctree node:pendingLeafNodesIn) {
    		if (!pendingLeafNodes.contains(node)) {
    			pendingLeafNodes.add(node);
    		}
    	}
    	// run splitting multithreaded
		final Thread[] threads = MultiThreading.newThreadArray();
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int indx_node = ai.getAndIncrement(); indx_node < pendingLeafNodes.size(); indx_node = ai.getAndIncrement()) {
						LwocOctree node = pendingLeafNodes.get(indx_node);
					    // Split recursively until satisfied conditions
					} // end of tile
				}
			};
		}
		MultiThreading.startAndJoin(threads);
		ai.set(0);
    }
    
    // Split recursively until satisfied conditions
    public int splitNode(
    		double min_hsize,
    		int    max_mesh_centers,
    		int    max_cameras
    		) {
    	if ((hsize <= min_hsize) || (leaf == null)) {
    		return 0;
    	}
    	if ((leaf.getScenes().size() <= max_cameras) &&
    			(leaf.getMeshCenters().size() <= max_mesh_centers)) {
    		return 0;
    	}
    	// actually need splitting
    	int num_added = NUM_CHILDREN - 1;
    	children = new LwocOctree[NUM_CHILDREN];
    	double new_hsize = hsize/2;
    	for (int nchild = 0; nchild < children.length; nchild++) {
    		double [] xyz = center.clone();
    		for (int dm = 0; dm < center.length; dm++) {
    			xyz[dm] +=new_hsize * ((((nchild >> dm) & 1) > 0)? 1 : -1);
    		}
    		children[nchild] = new LwocOctree(
    				this,
    				xyz,
    				new_hsize);
    		LwocLeaf new_leaf = new LwocLeaf();
    		children[nchild].leaf = new_leaf;
    	}
    	ArrayList <LwocMesh> meshes =       leaf.getMeshes();
    	// distribute meshes between 8 children
    	ArrayList <LwocMesh> mesh_centers = leaf.getMeshCenters();
    	// distribute mesh_centers  between 8 children
    	ArrayList <LwocScene> scenes =      leaf.getScenes();
    	// distribute scnes (cameras) between 8 children
    	// delete old leaf (if it is in a list) and its lists - no need
    	leaf = null;
    	for (int nchild = 0; nchild < children.length; nchild++) {
    		num_added += children[nchild].splitNode(
    				min_hsize,
    				max_mesh_centers,
    				max_cameras);
    	}
    	return num_added;
    }
    
    
    // not thread safe
    // grow world to include xyz (all 8 corners for a bounding box of as mesh
    
    
    /**
     * Grow the world to include specified point Will repeat growing twice (in each direction)
     * until the specified point gets inside.
     * @param xyz The point to be included in the world.
     */
    public static void growXYZ(double [] xyz) {
       	while (!lwoc_root.contains(xyz)) {// already fits, do not grow
    	    // grow once
    		int indx=0; //direction to grow
    		double [] new_center = new double[3];
    		for (int dm = 0; dm < new_center.length; dm++) {
    			if (xyz[dm] >= lwoc_root.center[dm] ) {
    				indx += 1 << dm;
        			new_center[dm]= lwoc_root.center[dm] + lwoc_root.hsize;
    			} else {
        			new_center[dm]= lwoc_root.center[dm] - lwoc_root.hsize;
    			}
    		}
    		LwocOctree new_root = new LwocOctree(null, new_center, lwoc_root.hsize * 2);
    		new_root.children = new LwocOctree[NUM_CHILDREN];
    		int child =  (~indx ) & (NUM_CHILDREN -1); // opposite direction, from new root to old root 
    		for (int ichild = 0; ichild < NUM_CHILDREN; ichild++) {
    			if (ichild == child) {
    				lwoc_root.parent = new_root;
    				new_root.children[ichild] = lwoc_root;
    			} else { // create empty leaf nodes
    				new_root.children[ichild] = new LwocOctree(
							new_root,
    						new double [] {
    								new_center[0] +	lwoc_root.hsize * ((((child >> 0) & 1) > 0)? 1 : -1),
    								new_center[1] +	lwoc_root.hsize * ((((child >> 1) & 1) > 0)? 1 : -1),
    								new_center[2] +	lwoc_root.hsize * ((((child >> 2) & 1) > 0)? 1 : -1)
    						},
    						lwoc_root.hsize);
    				
    				// add leaf 
    				new_root.children[ichild].leaf = new LwocLeaf();
    			}
    		}
    		lwoc_root = new_root;
    	}
    }
    /**
     * Test if the specified box intersects this node
     * @param xyz box center
     * @param half_whd  half width(x), height(y) and depth(z)
     * @return True if it intersects
     */
    public boolean intersects (
    		double [] xyz,
    		double [] half_whd) {
    	for (int dm = 0; dm < center.length; dm++) {
    		if ((xyz[dm] - half_whd[dm]) >= (center[dm] + hsize)) { // semiinterval
    			return false;
    		}
    		if ((xyz[dm] + half_whd[dm]) < (center[dm] - hsize)) {
    			return false;
    		}
    	}
    	return true;
    }
    
    /**
     * Test if the specified point fits in this node
     * @param xyz point
     * @return True if it intersects
     */
    public boolean contains (
    		double [] xyz) {
    	for (int dm = 0; dm < center.length; dm++) {
    		if (xyz[dm] < (center[dm] - hsize)) { // semi-interval
    			return false;
    		}
    		if (xyz[dm] >= (center[dm] + hsize)) {
    			return false;
    		}
    	}
    	return true;
    }
    
    
    
    /**
     * Test if the specified box completely fits in this node
     * @param xyz box center
     * @param half_whd  half width(x), height(y) and depth(z)
     * @return True if it intersects
     */
    public boolean contains (
    		double [] xyz,
    		double [] half_whd) {
    	for (int dm = 0; dm < center.length; dm++) {
    		if ((xyz[dm] - half_whd[dm]) < (center[dm] - hsize)) { // semiinterval
    			return false;
    		}
    		if ((xyz[dm] + half_whd[dm]) >= (center[dm] + hsize)) {
    			return false;
    		}
    	}
    	return true;
    }
    
    /**
     * Add mesh to each node intersecting with the mesh bounding box,
     * add to pending list if does not fit in a root node
     * @param mesh new mesh to add
     * @param check_existed Add only if the same mesh did not exist
     * @param pendingMeshes per thread list of meshes to be combined after merging threads
     * @param pendingLeafNodes  per thread list of nodes to be combined after merging threads
     * @return number of nodes it was added to (intersecting)
     */
    public static boolean addMeshCenter( // returns 0 and adds to pendingMeshes if needs growing
    		LwocMesh              mesh,
    		boolean               check_existed,
    		ArrayList<LwocMesh>   pendingMeshes, 
    		ArrayList<LwocOctree> pendingLeafNodes
    		) {
    	// should in check before adding?
    	// mesh.equals(mesh);
    	double [] center = mesh.getCenter();
    	double [] dims =   mesh.getDims();
    	// First - check that all 8 corners fit in the world
		double [] corner_xyz = center.clone();
    	for (int dm = 0; dm < center.length; dm++) {
    		corner_xyz[dm] += dims[dm];
    		if (!lwoc_root.contains(corner_xyz)){
    			if (pendingMeshes != null) {
    				pendingMeshes.add(mesh);
    			}
    			return false;
    		}
    		corner_xyz[dm] -= 2*dims[dm];
    		if (!lwoc_root.contains(corner_xyz)){
    			if (pendingMeshes != null) {
    				pendingMeshes.add(mesh);
    			}
    			return false;
    		}
    		corner_xyz[dm] += dims[dm];
    	}
    	// all corners fit - add center
    	LwocOctree node = getLeafNode(center);
		node.leaf.addMeshCenter(mesh);
		// Check if leaf node requires splitting
		if (node.leaf.getMeshCenters().size() > MAX_MESH_CENTERS) {
			if (node.hsize > MIN_HSIZE) {
				pendingLeafNodes.add(node);
			}
		}
    	// Now recursively add mesh to all nodes intersected by this mesh bounding box 
    	return true;
    }

    public int addMesh( // returns 0 and adds to pendingMeshes if needs growing
    		LwocMesh mesh,
    		boolean check_existed) {
    	if (!intersects(
    			mesh.getCenter(),
    			mesh.getDims())) {
    		return 0;
    	}
    	if (isLeaf()) {
    		leaf.addMesh(mesh, check_existed);
    		return 1;
    	}
    	int added=0;
    	for (LwocOctree node: children) {
    		added += node.addMesh(mesh, check_existed);
    	}
    	return added;
    }
    
    /**
     * Remove mesh to from each node intersecting with the mesh bounding box.
     * Should fit in a root node (as it was earlier added)
     * @param mesh new mesh to add
     * @param check_existed Add only if the same mesh did not exist 
     * @return The number of nodes it was removed from
     */
    public static int removeMesh(
    		LwocMesh mesh,
    		boolean remove_all) { // check and remove duplicates
    	// should in check before adding?
    	return 0;
    }

}