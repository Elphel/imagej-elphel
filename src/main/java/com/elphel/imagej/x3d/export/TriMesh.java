package com.elphel.imagej.x3d.export;
/**
 **
 ** TriMesh - triangular mesh representation
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  TriMesh.java is free software: you can redistribute it and/or modify
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

public class TriMesh {
	public String texture_image;
	public double [][] worldXYZ;
	public double [][] texCoord;
	public int [][]    triangles;
	public TriMesh (
			String texture_image,
			double [][] worldXYZ,
			double [][] texCoord,
			int [][] triangles) {
		this.texture_image = texture_image;
		this.worldXYZ = worldXYZ;
		this.texCoord = texCoord;
		this.triangles = triangles;
	}
	public String getImage() {return texture_image;}
	int [][] getTriangles()  {return triangles;}
	double [][] getTexCoord() {
		return texCoord;
	}
	public double [][] getTexCoord(boolean inv_x, boolean inv_y, boolean swap_xy) {
		if (!inv_x && !inv_y && !swap_xy) {
			return texCoord;
		}
		double [][] inv_tex_coord = new double [texCoord.length][2];
		double scale_x = inv_x ? -1.0 : 1.0;
		double scale_y = inv_y ? -1.0 : 1.0;
		if (swap_xy) {
			for (int i = 0; i <texCoord.length; i++) {
				inv_tex_coord[i][0] = scale_y * texCoord[i][1];
				inv_tex_coord[i][1] = scale_x * texCoord[i][0];
			}
		} else {
			for (int i = 0; i <texCoord.length; i++) {
				inv_tex_coord[i][0] = scale_x * texCoord[i][0];
				inv_tex_coord[i][1] = scale_y * texCoord[i][1];
			}
		}
		return inv_tex_coord;	
	}
	public double [][] getCoordinates(){
		return worldXYZ;
	}
	/**
	 * 0: XYZ -> XYZ
	 * 1: XYZ -> YZX
	 * 2: XYZ -> ZXY
	 * 3: XYZ -> XZY
	 * 4: XYZ -> ZYX
	 * 5: XYZ -> YXZ
	 * @param inv_x
	 * @param inv_y
	 * @param inv_z
	 * @param swap3
	 * @return
	 */
	public double [][] getCoordinates(boolean inv_x, boolean inv_y, boolean inv_z, int swap3) {
		if (!inv_x && !inv_y && !inv_y && (swap3 == 0)) {
			return worldXYZ;
		}
		double [][] inv_worldXYZ = new double [worldXYZ.length][3];
		double scale_x = inv_x ? -1.0 : 1.0;
		double scale_y = inv_y ? -1.0 : 1.0;
		double scale_z = inv_z ? -1.0 : 1.0;
		switch (swap3) {
		case 0: // XYZ -> XYZ
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][0] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][1] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][2] = scale_z * worldXYZ[i][2];
			}
			break;
		case 1: // XYZ -> YZX
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][1] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][2] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][0] = scale_z * worldXYZ[i][2];
			}
			break;
		case 2: // XYZ -> ZXY
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][2] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][0] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][1] = scale_z * worldXYZ[i][2];
			}
			break;
		case 3: // XYZ -> XZY
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][0] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][2] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][1] = scale_z * worldXYZ[i][2];
			}
			break;
		case 4: // XYZ -> ZYX
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][2] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][1] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][0] = scale_z * worldXYZ[i][2];
			}
			break;
		case 5: // XYZ -> YXZ
			for (int i = 0; i <texCoord.length; i++) {
				inv_worldXYZ[i][1] = scale_x * worldXYZ[i][0];
				inv_worldXYZ[i][0] = scale_y * worldXYZ[i][1];
				inv_worldXYZ[i][2] = scale_z * worldXYZ[i][2];
			}
			break;
		default: return null;
		}
		return inv_worldXYZ;
	}



}
