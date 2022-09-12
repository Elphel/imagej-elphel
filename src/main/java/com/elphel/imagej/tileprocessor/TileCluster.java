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
import java.util.Arrays;

class TileCluster{
	Rectangle  bounds;
	boolean [] border;
	double []  disparity;
	public TileCluster (Rectangle  bounds, boolean [] border, double []  disparity){
		this.bounds =    bounds;
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
	public Rectangle getBounds() {return bounds;}
	public boolean [] getBorder() {return border;}
	public double [] getDisparity() {return disparity;}
	/*
	public TileCluster combine (TileCluster tileCluster) {
		TileCluster outer, inner;
		if (bounds.contains(tileCluster.bounds)) {
			outer = this;
			inner = tileCluster;
		} else if (tileCluster.bounds.contains(bounds)) {
			outer = tileCluster;
			inner = this;
		} else {
			Rectangle outer_bounds = bounds.union(tileCluster.bounds);
			outer = new TileCluster(outer_bounds, null, null);
			outer.combine(this); // 
			inner = tileCluster;
		}
		int dst_x = inner.bounds.x - outer.bounds.x;
		for (int src_y = 0; src_y < bounds.height; src_y++) {
			int dst_y = src_y + inner.bounds.y - outer.bounds.y;
			System.arraycopy(
					inner.border,
					src_y * bounds.width,
					outer.border,
					dst_y * outer.bounds.width + dst_x,
					bounds.width);
			System.arraycopy(
					inner.disparity,
					src_y * bounds.width,
					outer.disparity,
					dst_y * outer.bounds.width + dst_x,
					bounds.width);
		}
		return outer;
	}
    */
	public void add(TileCluster tileCluster) {
		if (!bounds.contains(tileCluster.bounds)) {
			throw new IllegalArgumentException ("TileCluster.add(): Added cluster should fit into this ");
		}
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