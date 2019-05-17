package com.elphel.imagej.x3d.export;
import java.io.FileWriter;
import java.io.IOException;

/**
 **
 ** generate Wavefront representation of the model
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  WavefrontExport.java is free software: you can redistribute it and/or modify
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

import java.util.ArrayList;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters.CLTParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters.CorrectionParameters;
import com.elphel.imagej.dp.CLTPass3d;
import com.elphel.imagej.dp.GeometryCorrection;

import ij.Prefs;

public class WavefrontExport {
	static final String MTL_EXT=".mtl";
	static final String OBJ_EXT=".obj";
		GeometryCorrection                                     geometry_correction;
		public ArrayList <CLTPass3d>             clt_3d_passes;
		public  EyesisCorrectionParameters.CLTParameters       clt_parameters;
		public EyesisCorrectionParameters.CorrectionParameters correctionsParameters;
		public int debugLevel = 1;
		FileWriter obj_writer; // f0 = new FileWriter("output.txt");
		FileWriter mtl_writer; // f0 = new FileWriter("output.txt");
		Document    x3dDoc;
		Element     el_X3d;
		Element     el_Scene;
		Element     el_TopGroup;
		public int  max_line_length = 0; // 100; // 0 - no limit
		public int v_index =  1; // first index is 1, not 0
		public int vt_index = 1; // first index is 1, not 0
		public int f_index =  1; // first index is 1, not 0
		public String obj_path;
		public String mtl_path;

		public WavefrontExport(
				String                                          out_dir,
				String                                          project_name,
				EyesisCorrectionParameters.CLTParameters        clt_parameters,
				EyesisCorrectionParameters.CorrectionParameters correctionsParameters,
				GeometryCorrection                              geometry_correction,
				ArrayList <CLTPass3d>             clt_3d_passes) throws IOException{
			this.clt_parameters =        clt_parameters;
			this.correctionsParameters = correctionsParameters;
			this.geometry_correction =    geometry_correction;
			this.clt_3d_passes = clt_3d_passes;
			mtl_path = out_dir+Prefs.getFileSeparator()+project_name+MTL_EXT;
			obj_path = out_dir+Prefs.getFileSeparator()+project_name+OBJ_EXT;
			mtl_writer = new FileWriter(mtl_path);
			obj_writer = new FileWriter(obj_path); // .close();
			mtl_writer.write("#\n# Wavefront material file\n#\n");
			obj_writer.write("#\n# Wavefront object file\n#\n");
			obj_writer.write("mtllib ./"+project_name+MTL_EXT+"\n\n"); // add "./" to indicate relative path? // ./1488240527_408296.obj.mtl\n");
			
		}
		public void close()
		{
			try {
				mtl_writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
				obj_writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		public void addCluster(
				String      texture_path,
				String      material_id,
				double [][] texCoord,
				double [][] coordinate,
				int [][] triangles) throws IOException
		{
			// Write to material file
/*
newmtl material_0
Ka 0.200000 0.200000 0.200000
Kd 1.000000 1.000000 1.000000
Ks 1.000000 1.000000 1.000000
Tr 1.000000
illum 2
Ns 0.000000
map_Kd 1488240527_408296-img2-texture.png
 */
			mtl_writer.write("\nnewmtl "+material_id+"\n");
			mtl_writer.write("Ka 0.200000 0.200000 0.200000\n");
			mtl_writer.write("Kd 1.000000 1.000000 1.000000\n");
			mtl_writer.write("Ks 1.000000 1.000000 1.000000\n");
			mtl_writer.write("Tr 1.000000\n");
			mtl_writer.write("illum 2\n");
//			mtl_writer.write("illum 0\n"); // should it be 0 orf 2?
			mtl_writer.write("Ns 0.000000\n");
			mtl_writer.write("map_Kd "+texture_path+"\n");
			
			// Write OBJ file
			// start using new material
			
			obj_writer.write("usemtl "+material_id+"\n");
			// output all vertices
			obj_writer.write("# start from vertex :         "+ v_index+ "\n");
			obj_writer.write("# start from texture vertex:  "+ vt_index+"\n");
			obj_writer.write("# start from face (triangle): "+ f_index+ "\n");
			for (int i = 0; i < coordinate.length; i++){
				obj_writer.write(String.format("v %.3f %.3f %.3f\n", coordinate[i][0], coordinate[i][1], coordinate[i][2]));
			}
			// output all texture vertices
			for (int i = 0; i < texCoord.length; i++){
				obj_writer.write(String.format("vt %.5f %.5f\n", texCoord[i][0], texCoord[i][1]));
			}
			// output all triangles (faces)
			for (int i = 0; i < triangles.length; i++){
				obj_writer.write(String.format("f %d/%d %d/%d %d/%d\n",
///						triangles[i][0]+v_index,triangles[i][0]+vt_index, // actually v_index and vt_index should =be the same
///						triangles[i][1]+v_index,triangles[i][1]+vt_index,
///						triangles[i][2]+v_index,triangles[i][2]+vt_index));
						// wrong normals
				triangles[i][2]+v_index,triangles[i][2]+vt_index, // actually v_index and vt_index should =be the same
				triangles[i][1]+v_index,triangles[i][1]+vt_index,
				triangles[i][0]+v_index,triangles[i][0]+vt_index));
			}
			obj_writer.write("# end of material "+   material_id+"\n");
			obj_writer.write("# vertices:          "+         coordinate.length+"\n");
			obj_writer.write("# texture vertices:  "+ texCoord.length+"\n");
			obj_writer.write("# faces (triangles): "+ triangles.length+"\n");
			
			obj_writer.write("\n");
			// increment indices
			v_index +=  coordinate.length;
			vt_index += texCoord.length;
			f_index +=  triangles.length;
//f 27151/27139/27151 27141/27140/27141 27140/27141/27140			 * 
			
		}
		
		
		
}

