package com.elphel.imagej.tileprocessor;
/**
 ** TileCluster - clusters from disparity map to generate textured mesh
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TileCluster.java is free software: you can redistribute it and/or modify
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

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;

class TileCluster{
	Rectangle  bounds;
	boolean [] border;
	double []  disparity;            // all and only unused - NaN 
	int []     cluster_index = null; // for debug purposes, index of the source cluster
	int        index = -1;
	ArrayList<IndexedRectanle> clust_list;
	class IndexedRectanle{
		int index;
		Rectangle bounds;
		IndexedRectanle (
				int index,
				Rectangle bounds){
			this.index = index;
			this.bounds = bounds;
		}
	}
	// to use cluster_index - set index >= 0, <0 - do not use.
	public TileCluster (
			Rectangle  bounds,
			int index, // <0 to skip
			boolean [] border,
			double []  disparity){
		this.bounds =    bounds;
		this.index = index; 
		/**
		if (index >= 0) {
			this.cluster_index = new int [bounds.width * bounds.height];
			Arrays.fill(cluster_index, -1);
			if (disparity != null) {
				for (int i = 0; i < cluster_index.length; i++) if (!Double.isNaN(disparity[i])){
					cluster_index[i] = index;
				}
			}
		}
		*/
		if (border == null) {
			border = new boolean[bounds.width * bounds.height];
		}
		this.border =    border;
		if (disparity == null) {
			disparity = new double[bounds.width * bounds.height];
			Arrays.fill(disparity, Double.NaN);
		}
		this.disparity = disparity;
	}
	
	public Rectangle  getBounds() {return bounds;}
	public boolean [] getBorder() {return border;}
	public double []  getDisparity() {return disparity;}
	public double []  getSubDisparity(int indx) { // disparity should be NaN for unused !
		if (clust_list == null) {
			return null;
		}
		Rectangle sub_bounds = clust_list.get(indx).bounds;
		double [] sub_disparity = new double [sub_bounds.width * sub_bounds.height];
		int src_x = sub_bounds.x - bounds.x;
		for (int dst_y = 0; dst_y < sub_bounds.height; dst_y++) {
			int src_y = dst_y + sub_bounds.y - bounds.y;
			System.arraycopy(
					disparity,
					src_y * bounds.width + src_x,
					sub_disparity,
					dst_y * sub_bounds.width,
					sub_bounds.width);
		}		
		return sub_disparity;
	}
	
	public boolean []  getSubBorder(int indx) { // disparity should be NaN for unused !
		if (clust_list == null) {
			return null;
		}
		Rectangle sub_bounds = clust_list.get(indx).bounds;
		boolean [] sub_border = new boolean [sub_bounds.width * sub_bounds.height];
		int src_x = sub_bounds.x - bounds.x;
		for (int dst_y = 0; dst_y < sub_bounds.height; dst_y++) {
			int src_y = dst_y + sub_bounds.y - bounds.y;
			System.arraycopy(
					border,
					src_y * bounds.width + src_x,
					sub_border,
					dst_y * sub_bounds.width,
					sub_bounds.width);
		}		
		return sub_border;
	}
	
	public boolean []  getSubSelected(int indx) { // disparity should be NaN for unused !
		if (clust_list == null) {
			return null;
		}
		double [] sub_disparity = getSubDisparity(indx);
		boolean [] sub_selection = new boolean [sub_disparity.length];
		for (int i = 0; i < sub_disparity.length; i++) {
			sub_selection[i] = !Double.isNaN(sub_disparity[i]);
		}
		return sub_selection;
	}
	
	

	
	
	public boolean [] getSelected() {
		if (disparity == null) {
			return null;
		}
		boolean [] selected = new boolean [disparity.length];
		for (int i = 0; i < selected.length; i++) {
			selected[i] = !Double.isNaN(disparity[i]);
		}
		return selected;
	}
	
	public Rectangle [] getSubBounds() {
		if (clust_list == null) {
			return null;
		} else {
			Rectangle [] sub_bounds = new Rectangle[clust_list.size()];
			for (int i = 0; i <clust_list.size(); i++) {
				sub_bounds[i] = clust_list.get(i).bounds;
			}
			return sub_bounds;
		}
	}
	public int [] getSubIndices() {
		if (clust_list == null) {
			return null;
		} else {
			int [] sub_indices = new int[clust_list.size()];
			for (int i = 0; i <clust_list.size(); i++) {
				sub_indices[i] = clust_list.get(i).index;
			}
			return sub_indices;
		}
	}
	public int [] getClusterIndex() { // (Now not) just a debug feature, no need to optimize?
		if (clust_list == null) {
			return null;
		} else if (cluster_index == null){
			cluster_index = new int [bounds.width * bounds.height];
			Arrays.fill(cluster_index, -1);
			for (IndexedRectanle ir : clust_list) {
				for (int y_src = 0; y_src < ir.bounds.height; y_src++) {
					int y = y_src + ir.bounds.y - bounds.y;
					int offs = y * bounds.width + ir.bounds.x - bounds.x;
					for (int x_src = 0; x_src < ir.bounds.width; x_src++) {
						if (!Double.isNaN(disparity[offs])) {
							cluster_index[offs] = ir.index;
						}
						offs++;
					}
				}
			}
		}
		return cluster_index;

	}

	public void add(TileCluster tileCluster) {
		if (!bounds.contains(tileCluster.bounds)) {
			throw new IllegalArgumentException ("TileCluster.add(): Added cluster should fit into this ");
		}
		cluster_index =   null;
		if (clust_list == null) {
			clust_list =  new ArrayList<IndexedRectanle>();
		}
		clust_list.add(new IndexedRectanle(tileCluster.index, tileCluster.bounds));
		
		int dst_x = tileCluster.bounds.x - bounds.x;
		for (int src_y = 0; src_y < tileCluster.bounds.height; src_y++) {
			int dst_y = src_y + tileCluster.bounds.y - bounds.y;
			System.arraycopy(
					tileCluster.border,
					src_y * tileCluster.bounds.width,
					border,
					dst_y * bounds.width + dst_x,
					tileCluster.bounds.width);
			System.arraycopy(
					tileCluster.disparity,
					src_y * tileCluster.bounds.width,
					disparity,
					dst_y * bounds.width + dst_x,
					tileCluster.bounds.width);
		}
		return;
	}
	
	
	
}