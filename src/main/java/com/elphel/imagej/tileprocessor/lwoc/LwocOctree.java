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

//import java.io.FileInputStream;
//import java.io.FileOutputStream;
//import java.io.IOException;
//import java.io.ObjectInputStream;
//import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import com.elphel.imagej.common.MultiThreading;
public class LwocOctree implements Serializable {
	private static final long    serialVersionUID = 1L;
	public static final int      NUM_CHILDREN =      8;
	static List<LwocScene>       pendingScenes;   // synchronized - add not-yet-added scenes pending growing world 
	static List<LwocMesh>        pendingMeshes;   // synchronized - add not-yet-added meshes pending growing world
	static List<LwocOctree>      pendingLeafNodes;// synchronized - leaf nodes needed to be split
	
	int                          id;         // assign unique ID
	transient LwocWorld          world; 	
	transient LwocOctree         parent;
	LwocOctree []                children;
	LwocLeaf                     leaf;       //
	double []                    center;     // x-center, y-center, z-center
	double                       hsize;
	int                          infinity = 0; 
	
	public static void addPendingScene(LwocScene scene) {
		synchronized(pendingScenes) {
			pendingScenes.add(scene);
		}
	}
	
	public static void addPendingMesh(LwocMesh mesh) {
		synchronized(pendingMeshes) {
			pendingMeshes.add(mesh);
		}
	}
	
	public static void addPendingLeafNode(LwocOctree node) {
		synchronized(pendingLeafNodes) {
			pendingLeafNodes.add(node);
		}
	}

	/**
	 * Synchronized lists of pending scenes, meshes and nodes are used to delay modification of
	 * the world octree and its elements from the multithreaded environment. These lists should
	 * be tested after multithreaded methods that potentially can modify the world and then
	 * processes by a single-thread method. 
	 */
	public static void resetPending() {
		pendingScenes =    new ArrayList<LwocScene>();
		pendingMeshes =    new ArrayList<LwocMesh>();
		pendingLeafNodes = new ArrayList<LwocOctree>();
	}
	
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
		return (pendingScenes.size()  + pendingMeshes.size()) > 0;
	}
	
	/**
	 * Check if any leaf nodes need splitting
	 * @return if any leaf nodes require splitting
	 */
	public static boolean pendingSplit() {
		return pendingLeafNodes.size() > 0;
	}
	
	/**
	 * LwocOctree constructor. Sets infinity bitmap +1 - no limit for negative X,
	 *  +2 - no limit for positive X, +4 -  no limit for negative Y, ...,
	 *  +32 - no limit for positive Z.
	 * @param lwocWorld World global parameters
	 * @param parent_in Parent node
	 * @param xyz       Position of the node center in meters
	 * @param h_size    Node half size - each of the X,Y,Z coordinates of the internal 
	 *                  points are limited within +/-hsize from the node center
	 */
	public LwocOctree (
			LwocWorld  lwocWorld,
			LwocOctree parent_in,
			double []  xyz,
			double     h_size) {
		world =  lwocWorld;
    	parent = parent_in;
    	center = xyz;
    	hsize =  h_size;
    	// set which directions are to infinity (not limited by hsize)
    	infinity = 0;
    	double max_hsize = world.getMaxHsize();
    	LwocOctree lwoc_root = world.getLwocRoot();
    	for (int dm = 0; dm < center.length;  dm++) {
    		if ((center[dm] - hsize) <= (lwoc_root.center[dm] - max_hsize)) {
    			infinity |= (1 << (2 * dm));
    		}
    		if ((center[dm] + hsize) >= (lwoc_root.center[dm] + max_hsize)) {
    			infinity |= (1 << (2 * dm + 1));
    		}
    	}
    }
    
    /**
     * Find the leaf node that includes specified point. Thread safe.
     * @param xyz X, Y, and Z coordinates of the point
     * @return Octree node that includes the point or null if the point
     *         is outside the root node (the whole world)
     */
    public LwocOctree getLeafNode( // if null - needs growing world
    		double [] xyz) { // thread safe
    	if (!world.getLwocRoot().contains(xyz)) return null; // needs growing world
    	LwocOctree node = world.getLwocRoot();
    	while (world.getLwocRoot().children != null) {
    		int indx = 0;
    		for (int dm = 0; dm < 3; dm++) {
    			if (xyz[dm] >= world.getLwocRoot().center[dm] ) {
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
     * @param check_existed - do not add duplicate scenes and nodes
     * @return An octree node where scene is added or null if the world needs growing
     */
    public LwocOctree addScene(
    		LwocScene             scene,
    		boolean               check_existed
    		) {// returns null and adds to pendingScenes if needs growing
    	LwocOctree node = getLeafNode(scene.getCameraXYZ());
    	if (node == null) {
    		if (pendingScenes != null) {
    			if (!check_existed || !pendingScenes.contains(scene)) {
    				addPendingScene(scene);
    			}
    		}
    	} else {
    		node.leaf.addScene(scene,check_existed);
    		// Check if leaf node requires splitting
    		if (node.leaf.getScenes().size() > world.getMaxCameras()) {
    			if (node.hsize > world.getMinHsize()) {
        			if (!check_existed || !pendingLeafNodes.contains(node)) {
        				addPendingLeafNode(node);
        			}
    			}
    		}
    	}
    	return node;
    }
    
    /**
     * Tend to pending scenes list. Not thread-safe, should be run in a single-thread mode.
     * @param check_existed - do not add duplicate scenes and nodes
     */
    public void tendPendingScenes(
    		boolean               check_existed
    		) {
    	if (pendingScenes.isEmpty()) {
    		return; // nothing to do
      	} else {
    		for (LwocScene scene: pendingScenes) {
    			growXYZ(scene.getCameraXYZ());  // will not grow in infinity direction
    			LwocOctree node= addScene(  // should not be null
    					scene,
    		    		check_existed);
    			assert node != null : "addScene() should not fail after growing world";
    		} 
    	}
    }
    
    /**
     * Tend to pending meshes list. Not thread-safe, should be run in a single-thread mode.
     * Only adds mesh centers and grows the world if needed.
     * @param check_existed    do not add duplicate meshes
     */
    public void tendPendingMeshCentersOnly(
    		boolean               check_existed
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
    	    		growXYZ(corner_xyz); // will not grow in infinity direction
    	    		corner_xyz[dm] -= 2*dims[dm];
    	    		growXYZ(corner_xyz); // will not grow in infinity direction
    	    		corner_xyz[dm] += dims[dm];
    	    	}    			
    			boolean ok= addMeshCenter(  // should not be null
    					mesh,
    					check_existed);
    			assert ok : "addMeshCenter() should not fail after growing world";
    		} 
    	}
    }

    /**
     * Tend to pending meshes list after thge world is grown and related mesh centers are added.
     * Only adds mesh themselves that do not trigger node splits.
     * Should be run after tendPendingMeshCentersOnly().
     * @param check_existed    do not add duplicate meshes
     */
    public void tendPendingMeshes(
    		boolean               check_existed ) {
    	if (pendingMeshes.isEmpty()) {
    		return; // nothing to do
      	} else {
    		for (LwocMesh mesh: pendingMeshes) {
//    			int num_added = 
    			world.lwoc_root.addMesh(
    					mesh,
    					check_existed);
    		} 
    	}
    }
    
    /**
     * Recursively (if needed) splits octree nodes to reduce number of scenes/cameras
     * and mesh centers below specified thresholds (limited by the minimal node size) 
     * @param check_existed Filter identical meshes or scenes from corresponding lists.
     */
    public static void tendPendingLeafNodes(
    		final boolean         check_existed) {
    	// remove any possible duplicates
    	final ArrayList<LwocOctree> pendingLeafNodesFiltered = new ArrayList<LwocOctree>();
    	for (LwocOctree node:pendingLeafNodes) {
    		if (!pendingLeafNodesFiltered.contains(node) && node.isLeaf()) {// filter already split nodes 
    			pendingLeafNodesFiltered.add(node);
    		}
    	}
    	if (!pendingLeafNodesFiltered.isEmpty()) {
    		LwocWorld world = pendingLeafNodesFiltered.get(0).world;
    		final double min_hsize =        world.getMinHsize();
    		final int    max_mesh_centers = world.getMaxMeshCenters();
    		final int    max_cameras =      world.getMaxCameras();
    		// run splitting multithreaded
    		final Thread[] threads = MultiThreading.newThreadArray();
    		final AtomicInteger ai = new AtomicInteger(0);
    		for (int ithread = 0; ithread < threads.length; ithread++) {
    			threads[ithread] = new Thread() {
    				@Override
    				public void run() {
    					for (int indx_node = ai.getAndIncrement(); indx_node < pendingLeafNodesFiltered.size(); indx_node = ai.getAndIncrement()) {
    						LwocOctree node = pendingLeafNodesFiltered.get(indx_node);
    						// Split recursively until satisfied conditions
    						node.splitNode(
    								min_hsize,        // double  min_hsize,
    								max_mesh_centers, // int     max_mesh_centers,
    								max_cameras,      // int     max_cameras,
    								check_existed);   // boolean check_existed)
    					} // end of tile
    				}
    			};
    		}
    		MultiThreading.startAndJoin(threads);
    	}
    }
    
    /**
     * Split recursively the octree node (and its children if needed) to keep number of scenes and
     * mesh centers below thresholds.  
     * @param min_hsize do not split nodes smaller than this (half size in each direction) even if
     *        it has over-threshold number of scenes and/or mesh centers. 
     * @param max_mesh_centers Split node if it has more mesh centers in its list.
     * @param max_cameras Split node if it has more scenes (camera positions) than this.
     * @param check_existed 
     * @return
     */
    public int splitNode(
    		double  min_hsize,
    		int     max_mesh_centers,
    		int     max_cameras,
    		boolean check_existed
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
    				world,
    				this,
    				xyz,
    				new_hsize);
    		LwocLeaf new_leaf = new LwocLeaf();
    		children[nchild].leaf = new_leaf;
    	}
    	ArrayList <LwocMesh> meshes =       leaf.getMeshes();
    	// distribute meshes between 8 children
    	for (LwocMesh mesh: meshes) { // may go to any number >=1 of the children
    		double [] mesh_center = mesh.getCenter();
    		double [] mesh_dims = mesh.getDims();
    		for (LwocOctree node:children) {
    			if (node.intersects(
    					mesh_center,
    					mesh_dims)) {
    				node.leaf.addMesh(
    						mesh,
    						check_existed);
    			}
    		}
    	}
    	ArrayList <LwocMesh> mesh_centers = leaf.getMeshCenters();
    	// distribute mesh_centers  between 8 children
    	for (LwocMesh mesh: mesh_centers) { // goes to a single child only
    		double [] mesh_center = mesh.getCenter();
    		for (LwocOctree node:children) {
    			if (node.contains(mesh_center)) {
    				node.leaf.addMeshCenter(
    						mesh,
    						check_existed);
    				break; // single point - goes to one child only 
    			}
    		}
    	}
    	ArrayList <LwocScene> scenes =      leaf.getScenes();
    	// distribute scenes (cameras) between 8 children
    	for (LwocScene scene: scenes) { // goes to a single child only
    		double [] camera_xyz = scene. getCameraXYZ();
    		for (LwocOctree node:children) {
    			if (node.contains(camera_xyz)) {
    				node.leaf.addScene(
    						scene,
    						check_existed);
    				break; // single point - goes to one child only 
    			}
    		}
    	}
    	// delete old leaf (if it is in a list) and its lists - no need
    	leaf = null;
    	for (int nchild = 0; nchild < children.length; nchild++) {
    		num_added += children[nchild].splitNode(
    				min_hsize,
    				max_mesh_centers,
    				max_cameras,
    				check_existed);
    	}
    	return num_added;
    }
    
    /**
     * Grow the world to include specified point Will repeat growing twice (in each direction)
     * until the specified point gets inside. Not thread safe, should run in single-thread mode
     * Will mark infinity directions (inside LwocOctree constructor), so contains() will return
     * true when the world size is big enough.
     * @param xyz The point to be included in the world.
     */
    public void growXYZ(double [] xyz) {
       	while (!world.lwoc_root.contains(xyz)) {// already fits, do not grow. Will obey infinity
    	    // grow once
    		int indx=0; //direction to grow
    		double [] new_center = new double[3];
    		double [] root_center = world.lwoc_root.center;
			double root_hsize = world.lwoc_root.hsize;
    		for (int dm = 0; dm < new_center.length; dm++) {
    			if (xyz[dm] >= root_center[dm] ) {
    				indx += 1 << dm;
        			new_center[dm]= root_center[dm] + root_hsize;
    			} else {
        			new_center[dm]= root_center[dm] - root_hsize;
    			}
    		}
    		LwocOctree new_root = new LwocOctree( // will mark infinities
    				world,
    				null,
    				new_center,
    				world.lwoc_root.hsize * 2);
    		new_root.children = new LwocOctree[NUM_CHILDREN];
    		int child =  (~indx ) & (NUM_CHILDREN -1); // opposite direction, from new root to old root 
    		for (int ichild = 0; ichild < NUM_CHILDREN; ichild++) {
    			if (ichild == child) {
    				world.lwoc_root.parent = new_root;
    				new_root.children[ichild] = world.lwoc_root;
    			} else { // create empty leaf nodes
    				new_root.children[ichild] = new LwocOctree( // will mark infinities
    						world,
							new_root,
    						new double [] {
    								new_center[0] +	root_hsize * ((((child >> 0) & 1) > 0)? 1 : -1),
    								new_center[1] +	root_hsize * ((((child >> 1) & 1) > 0)? 1 : -1),
    								new_center[2] +	root_hsize * ((((child >> 2) & 1) > 0)? 1 : -1)
    						},
    						root_hsize);
    				
    				// add leaf 
    				new_root.children[ichild].leaf = new LwocLeaf();
    			}
    		}
    		world.lwoc_root = new_root;
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
    	if (infinity == 0) { // this branch is just for performance, it may be removed 
    		for (int dm = 0; dm < center.length; dm++) {
    			if ((xyz[dm] - half_whd[dm]) >= (center[dm] + hsize)) { // semiinterval
    				return false;
    			}
    			if ((xyz[dm] + half_whd[dm]) < (center[dm] - hsize)) {
    				return false;
    			}
    		}
    	} else {
    		for (int dm = 0; dm < center.length; dm++) {
    			if (((xyz[dm] - half_whd[dm]) >= (center[dm] + hsize))  && (((infinity >> (2 * dm + 1)) & 1 ) == 0)) { // semiinterval
    				return false;
    			}
    			if (((xyz[dm] + half_whd[dm]) < (center[dm] - hsize))   && (((infinity >> (2 * dm)) & 1 ) == 0)) {
    				return false;
    			}
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
    	if (infinity == 0) { // this branch is just for performance, it may be removed 
    		for (int dm = 0; dm < center.length; dm++) {
    			if (xyz[dm] < (center[dm] - hsize)) { // semi-interval
    				return false;
    			}
    			if (xyz[dm] >= (center[dm] + hsize)) {
    				return false;
    			}
    		}
    		return true;
    	} else {
    		for (int dm = 0; dm < center.length; dm++) {
    			if ((xyz[dm] < (center[dm] - hsize)) && (((infinity >> (2 * dm)) & 1 ) == 0)) { // semi-interval
    				return false;
    			}
    			if ((xyz[dm] >= (center[dm] + hsize)) && (((infinity >> (2 * dm + 1)) & 1 ) == 0)) {
    				return false;
    			}
    		}
    		return true;
    	}
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
    	if (infinity == 0) { // this branch is just for performance, it may be removed 
    		for (int dm = 0; dm < center.length; dm++) {
    			if ((xyz[dm] - half_whd[dm]) < (center[dm] - hsize)) { // semiinterval
    				return false;
    			}
    			if ((xyz[dm] + half_whd[dm]) >= (center[dm] + hsize)) {
    				return false;
    			}
    		}
    	} else {
    		for (int dm = 0; dm < center.length; dm++) {
    			if (((xyz[dm] - half_whd[dm]) < (center[dm] - hsize)) && (((infinity >> (2 * dm)) & 1 ) == 0)){ // semiinterval
    				return false;
    			}
    			if (((xyz[dm] + half_whd[dm]) >= (center[dm] + hsize)) && (((infinity >> (2 * dm + 1)) & 1 ) == 0)){
    				return false;
    			}
    		}
    	}
    	return true;
    }
    
    /**
     * Add mesh to each node intersecting with the mesh bounding box,
     * add to pending list if does not fit in a root node
     * @param mesh new mesh to add
     * @param check_existed Add only if the same mesh did not exist
     * @return number of nodes it was added to (intersecting)
     */
    public boolean addMeshCenter( // returns 0 and adds to pendingMeshes if needs growing. Obeys infinity
    		LwocMesh              mesh,
    		boolean               check_existed
    		) {
    	double [] center = mesh.getCenter();
    	double [] dims =   mesh.getDims();
    	// First - check that all 8 corners fit in the world
		double [] corner_xyz = center.clone();
    	for (int dm = 0; dm < center.length; dm++) {
    		corner_xyz[dm] += dims[dm];
    		if (!world.lwoc_root.contains(corner_xyz)){
    			if (pendingMeshes != null) {
    				if (!check_existed || !pendingMeshes.contains(mesh)) {
    					addPendingMesh(mesh);
    				}
    			}
    			return false;
    		}
    		corner_xyz[dm] -= 2*dims[dm];
    		if (!world.lwoc_root.contains(corner_xyz)){
    			if (pendingMeshes != null) {
    				if (!check_existed || !pendingMeshes.contains(mesh)) {
    					addPendingMesh(mesh);
    				}
    			}
    			return false;
    		}
    		corner_xyz[dm] += dims[dm];
    	}
    	// all corners fit - add center
    	LwocOctree node = getLeafNode(center);
		node.leaf.addMeshCenter(mesh, check_existed);
		// Check if leaf node requires splitting
		if (node.leaf.getMeshCenters().size() > world.getMaxMeshCenters()) {
			if (node.hsize >  world.getMinHsize()) {
				if (!check_existed || !pendingLeafNodes.contains(node)) {
					addPendingLeafNode(node);
				}
			}
		}
    	// Not here: Now recursively add mesh to all nodes intersected by this mesh bounding box 
    	return true;
    }

    /**
     * Recursively add mesh to all leaf nodes that intersect with its bounding box. 
     * Usually starts with lwoc_root.
     * @param mesh Mesh instance to be added
     * @param check_existed verify node lists does not already have this mesh if true
     * @return number of nodes to which this mesh was added
     */
    public int addMesh(
    		LwocMesh mesh,
    		boolean check_existed) {
    	if (!intersects( // supports infinity
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
     * Remove mesh center from the world.
     * Should fit in a root node (as it was earlier added)
     * @param mesh new mesh to add
     */
    public void removeMeshCenter(
    		LwocMesh mesh) { // check and remove duplicates
    	LwocOctree node = getLeafNode(mesh.getCenter());
		node.leaf.removeMeshCenter(mesh);
    }

    /**
     * Remove mesh from each node intersecting with the mesh bounding box.
     * Should fit in a root node (as it was earlier added)
     * @param mesh mesh to remove
     */
    public int removeMesh(
    		LwocMesh mesh) { // check and remove duplicates
    	if (!intersects(
    			mesh.getCenter(),
    			mesh.getDims())) {
    		return 0;
    	}
    	if (isLeaf()) {
    		leaf.removeMesh(mesh);
    		return 1;
    	}
    	int removed = 0;
    	for (LwocOctree node: children) {
    		removed += node.removeMesh(mesh);
    	}
    	return removed;
    }
    
    /**
     * Restore pointers to .parent and .world after de-serialization
     * @param world
     * @param parent
     */
    public void restoreWorldParent(
    		LwocWorld world,
    		LwocOctree parent) {
    	this.parent = parent;
    	this.world = world;
    	if (!isLeaf()) {
    		for (LwocOctree child:children) {
    			child.restoreWorldParent(
    		    		world, // LwocWorld world,
    		    		this); // LwocOctree parent)
    		}    		
    	}
    }
    
    
    /**
     * Recursively create a list of all leave nodes in the world.
     * @return A list of all leaf nodes
     */
    public List<LwocOctree> getAllLeafNodes(){
    	List<LwocOctree> leaves = new ArrayList<LwocOctree>();
    	if (isLeaf()) {
    		leaves.add(this);
    	} else {
    		for (LwocOctree child:children) {
    			leaves.addAll(child.getAllLeafNodes());
    		}
    	}
    	return leaves;
    }
    
    /**
     * Rebuild meshes lists from mesh_centers for all leaf nodes after restoring the world.
     * Used after restoring world from serialized form
     */
    public void rebuildMeshLists(
    		final boolean check_existed) {
    	final List<LwocOctree> leaves = world.lwoc_root.getAllLeafNodes();
    	final List<LwocMesh> meshes = new ArrayList<LwocMesh>();
    	for (LwocOctree node : leaves) {
    		node.leaf.initMeshes();
    		meshes.addAll(node.leaf.getMeshCenters());
    	}
    	
		final Thread[] threads = MultiThreading.newThreadArray();
		final AtomicInteger ai = new AtomicInteger(0);
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int indx_mesh = ai.getAndIncrement(); indx_mesh < meshes.size(); indx_mesh = ai.getAndIncrement()) {
						LwocMesh mesh = meshes.get(indx_mesh);
						world.lwoc_root.addMesh(mesh, check_existed);
					}
				}
			};
		}
		MultiThreading.startAndJoin(threads);
    }
}
