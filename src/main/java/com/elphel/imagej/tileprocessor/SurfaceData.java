package com.elphel.imagej.tileprocessor;
/**
 **
 ** SurfaceData - represents rectangular area with per-tile disparity
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  SurfaceData.java is free software: you can redistribute it and/or modify
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

public class SurfaceData {
	private int []     window;
	private int        tilesX;
	private int        tilesY;
	private double []  disparity;
	private boolean [] selected; // should include border
	private boolean [] border;
	private int []     sizes;
	private int []     clusters;

	public  SurfaceData(int [] tilesWH, int [] window)
	{
		this.tilesX = tilesWH[0];
		this.tilesY = tilesWH[1];
		this.window = window;
		int len = window[2]*window[3];
		this.disparity = new double [len];
		this.selected =  new boolean [len];
		this.border =    new boolean [len];
	}
	public SurfaceData(int tilesX, int tilesY, int x0, int y0, int width,int height)
	{
		this.tilesX = tilesX;
		this.tilesY = tilesY;
		int [] window = {x0,y0,width,height}; 
		this.window = window;
		int len = window[2]*window[3];
		this.disparity = new double [len];
		this.selected =  new boolean [len];
		this.border =    new boolean [len];
	}
	
	public void setCluster(int [] clusters)
	{
		this.clusters = clusters;
	}
	public int []  getCluster()
	{
		return this.clusters;
	}
	
	public int getWidth()
	{
		return window[2];
	}
	public int getHeight()
	{
		return window[3];
	}
	public int getLeft()
	{
		return window[0];
	}
	public int getTop()
	{
		return window[1];
	}
	public int [] getWindow()
	{
		return window;
	}
	
	
	public void setDisparity(double [] disparity)
	{
		this.disparity = disparity;
	}
	public void setDisparity(int indx, double disparity)
	{
		this.disparity[indx] = disparity;
	}
	public double [] getDisparity()
	{
		return this.disparity;
	}
	public double getDisparity(int indx)
	{
		return this.disparity[indx];
	}

	
	public void setSelected(boolean [] selected)
	{
		this.selected = selected;
		this.sizes = null;
	}
	public void setSelected(int indx, boolean selected)
	{
		this.selected[indx] = selected;
		this.sizes = null;
	}
	public boolean [] getSelected()
	{
		return this.selected;
	}
	public boolean getSelected(int indx)
	{
		return this.selected[indx];
	}

	
	public void setBorder(boolean [] border)
	{
		this.border = border;
		this.sizes = null;
	}
	public void setBorder(int indx, boolean border)
	{
		this.border[indx] = border;
		this.sizes = null;
	}
	public boolean [] getBorder()
	{
		return this.border;
	}
	public boolean getBorder(int indx)
	{
		return this.border[indx];
	}
	
	int getImageIndex(int indx)
	{
		int x = (indx % window[2]) + window[0];
		int y = (indx / window[2]) + window[1];
		if ((y <0) || (x < 0) || (y >= tilesY) || (x >= tilesX)){
			return -1;
		}
		return  y * tilesX + x;
	}
	
	public int getSize(
			boolean include_border)
	{
		if (this.sizes == null){
			int [] sizes = {0,0};
			this.sizes = sizes;
			for (int i = 0; i < selected.length; i++){
				if (selected[i]) {
					sizes[1]++;
					if (!border[i]){
						sizes[0]++;
					}
				}
			}
		}
		return this.sizes[include_border ? 1 : 0];

	}
	public double [] getDisparityNaN(
			boolean include_border)
	{
		double [] data = new double [tilesX*tilesY];
		for (int i = 0; i < data.length; i++){
			data[i] = Double.NaN;
		}
		for (int i = 0; i < disparity.length; i++){
			if (selected[i] && (include_border || !border[i])){
				data[getImageIndex(i)] = disparity[i];
			}
		}
		return data;
	}

}
