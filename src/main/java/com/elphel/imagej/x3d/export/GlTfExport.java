package com.elphel.imagej.x3d.export;
/**
 **
 ** GlTfExport - generate glTF representation of the model
 **
 ** Copyright (C) 2022 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  GlTfExport.java is free software: you can redistribute it and/or modify
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


import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.HashMap;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;

import ij.Prefs;

public class GlTfExport {
	public static int NEAREST =                  9728;
	public static int LINEAR =                   9729;
	public static int NEAREST_MIPMAP_NEAREST =   9984;
	public static int LINEAR_MIPMAP_NEAREST =    9985;
	public static int NEAREST_MIPMAP_LINEAR =    9986;
	public static int LINEAR_MIPMAP_LINEAR =     9986;
	
	public static int CLAMP_TO_EDGE =           33071;
	public static int MIRRORED_REPEAT =         33648;
	public static int REPEAT =                  10497;

	public static int ARRAY_BUFFER =            34962; // for attributes
	public static int ELEMENT_ARRAY_BUFFER =    34963; // for indices
	
	public static int COMPONENT_SIGNED_BYTE =    5120;
	public static int COMPONENT_UNSIGNED_BYTE =  5121;
	public static int COMPONENT_SIGNED_SHORT =   5122;
	public static int COMPONENT_UNSIGNED_SHORT = 5123;
	public static int COMPONENT_UNSIGNED_INT =   5125;
	public static int COMPONENT_FLOAT =          5126;
	
	public static String ALPHAMODE_OPAQUE = "OPAQUE";
	public static String ALPHAMODE_MASK =   "MASK";
	public static String ALPHAMODE_BLEND =  "BLEND";
	
	public static String TYPE_SCALAR =      "SCALAR"; //  1
	public static String TYPE_VEC2 =        "VEC2";   //  2
	public static String TYPE_VEC3 =        "VEC3";   //  3
	public static String TYPE_VEC4 =        "VEC4";   //  4
	public static String TYPE_MAT2 =        "MAT2";   //  4
	public static String TYPE_MAT3 =        "MAT3";   //  9
	public static String TYPE_MAT4 =        "MAT4";   // 16

	public static int MODE_POINTS =         0;
	public static int MODE_LINES =          1;
	public static int MODE_LINE_LOOP =      2;
	public static int MODE_LINE_STRIP =     3;
	public static int MODE_TRIANGLES =      4;
	public static int MODE_TRIANGLE_STRIP = 5;
	public static int MODE_TRIANGLE_FAN =   6;

	
	public static void glTFExport(
			String             x3d_dir,
			String             model_name,
			ArrayList<TriMesh> tri_meshes,
			boolean            gltf_emissive,
			boolean            use_alpha_blend,
			int                debugLevel
			) throws IOException, JSONException {
		boolean invert_faces = true; // false; //  true; // false; // true;
		boolean [] inv_xyz = {false,false,false};
		boolean [] inv_tex = {false,true};
		int [] triangle = {0,1,2};
		if (invert_faces) {
			triangle[0]=2;
			triangle[2]=0;
		}
		String base_name = x3d_dir + Prefs.getFileSeparator() + model_name;
		String bin_path = base_name+".bin";
		String bin_name = model_name+".bin";
		String gltf_path = base_name+".gltf";
		int total_meshes =    tri_meshes.size(); 
		int total_vertices =  0;
//		int total_triangles = 0;
		int [] mesh_vert = new int [total_meshes];
		int [] mesh_tri =  new int [total_meshes];
		int [][] offset_len_indx = new int [total_meshes][2];
		int [][] offset_len_tex =  new int [total_meshes][2];
		int [][] offset_len_coor = new int [total_meshes][2];
		int [][]     minmax_indx = new int [total_meshes][2];
		float [][][] minmax_tex =  new float [total_meshes][2][2]; // [][{x,y}[{min, max}]
		float [][][] minmax_coor = new float [total_meshes][3][2]; // [][{x,y,z}[{min, max}]
		int index_size = 4;// 2; // 4 bytes per index entry - will use int - probably faster render even if more memory
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			mesh_vert[nmesh] = tri_meshes.get(nmesh).getCoordinates().length;
			mesh_tri[nmesh] = tri_meshes.get(nmesh).getTriangles().length;
			total_vertices += mesh_vert[nmesh];
///			total_triangles += mesh_tri[nmesh];
			offset_len_indx[nmesh][1] = mesh_tri[nmesh] * index_size * 3;
			offset_len_tex [nmesh][1] = mesh_vert[nmesh] * 4 * 2; // float, VEC2
			offset_len_coor[nmesh][1] = mesh_vert[nmesh] * 4 * 3; // float, VEC3
			if (nmesh > 0) {
				int index_prev_len4 = 4 * ((offset_len_indx[nmesh - 1][1] + 3)/4);
				offset_len_indx[nmesh][0] = offset_len_indx[nmesh - 1][0] + index_prev_len4; // pad to 4 bytes
				offset_len_tex [nmesh][0] = offset_len_tex [nmesh - 1][0] + offset_len_tex [nmesh - 1][1];
				offset_len_coor[nmesh][0] = offset_len_coor[nmesh - 1][0] + offset_len_coor[nmesh - 1][1]; 
			}
		}
		
		int buffer_index_bytes = offset_len_indx[total_meshes-1][0]+offset_len_indx[total_meshes-1][1];
		buffer_index_bytes =  4 * ((buffer_index_bytes + 3)/4);
		int buffer_texture_bytes =    4 * total_vertices * 2; // VEC2
		int buffer_coordinate_bytes = 4 * total_vertices * 3; // VEC3
		int buffer_total_bytes = buffer_index_bytes+ buffer_texture_bytes+ buffer_coordinate_bytes;
		if (debugLevel > -1) {
			System.out.println("glTFExport(): allocating "+buffer_total_bytes+" bytes buffer");
			System.out.println(buffer_index_bytes+" for indices,\n" +
			buffer_texture_bytes+" for texture coordinates, and \n" +
			buffer_coordinate_bytes+" for vertice coordinates.");
		}
		ByteBuffer bb0 = ByteBuffer.allocate(buffer_total_bytes);
		ByteBuffer bb = bb0.order(ByteOrder.LITTLE_ENDIAN);
		// will be 3 buffer views: 0 for indices, 1 - for texture coordinates and 2 - for coordinates
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			int [][] triangles = tri_meshes.get(nmesh).getTriangles();
			if (index_size == 2) {
				for (int i = 0; i < triangles.length; i++) {
					bb.putShort((short) (triangles[i][triangle[0]]));
					bb.putShort((short) (triangles[i][triangle[1]]));
					bb.putShort((short) (triangles[i][triangle[2]]));
					bb.putShort((short) 0);
				}
			} else {
				for (int i = 0; i < triangles.length; i++) {
					bb.putInt(triangles[i][triangle[0]]);
					bb.putInt(triangles[i][triangle[1]]);
					bb.putInt(triangles[i][triangle[2]]);
				}
			}
			minmax_indx[nmesh][0] = triangles[0][0]; // zero triangles:Index 0 out of bounds for length 0
			minmax_indx[nmesh][1] = minmax_indx[nmesh][0];
			for (int i = 0; i < triangles.length; i++) {
				for (int j = 0; j < triangles[i].length; j++) {
					if (triangles[i][j] < minmax_indx[nmesh][0]) {
						minmax_indx[nmesh][0] = triangles[i][j];
					} else if (triangles[i][j] > minmax_indx[nmesh][1]) {
						minmax_indx[nmesh][1] = triangles[i][j];
					}
				}
			}
		}
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			double [][] tex_coord =   tri_meshes.get(nmesh).getTexCoord();
			for (int i = 0; i < tex_coord.length; i++) {
				for (int j = 0; j <  tex_coord[i].length; j++) {
					if (inv_tex[j]) tex_coord[i][j] = 1.0-tex_coord[i][j];
				}
			}
			for (int j = 0; j < tex_coord[0].length; j++) {
				minmax_tex[nmesh][j][0] = (float) tex_coord[0][j];
				minmax_tex[nmesh][j][1] = minmax_tex[nmesh][j][0];
			}
			for (int i = 0; i < tex_coord.length; i++) {
				bb.putFloat((float) tex_coord[i][0]);
				bb.putFloat((float) tex_coord[i][1]);
				for (int j = 0; j < tex_coord[i].length; j++) {
					float f = (float) tex_coord[i][j];
					if (f < minmax_tex[nmesh][j][0]) {
						minmax_tex[nmesh][j][0] = f;
					} else if (f > minmax_tex[nmesh][j][1]) {
						minmax_tex[nmesh][j][1] = f;
					} 
				}
			}
		}		
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			double [][] coordinates = tri_meshes.get(nmesh).getCoordinates();
			for (int i = 0; i < coordinates.length; i++) {
				for (int j = 0; j <  coordinates[i].length; j++) {
					if (inv_xyz[j]) coordinates[i][j] = -coordinates[i][j];
				}
			}
			for (int j = 0; j < coordinates[0].length; j++) {
				minmax_coor[nmesh][j][0] = (float) coordinates[0][j];
				minmax_coor[nmesh][j][1] = minmax_coor[nmesh][j][0];
			}
			
			for (int i = 0; i < coordinates.length; i++) {
				bb.putFloat((float) coordinates[i][0]);
				bb.putFloat((float) coordinates[i][1]);
				bb.putFloat((float) coordinates[i][2]);
				for (int j = 0; j < coordinates[i].length; j++) {
					float f = (float) coordinates[i][j];
					if (f < minmax_coor[nmesh][j][0]) {
						minmax_coor[nmesh][j][0] = f;
					} else if (f > minmax_coor[nmesh][j][1]) {
						minmax_coor[nmesh][j][1] = f;
					} 
				}
			}		
		}
		bb.flip();
		try (// write bb to file
			@SuppressWarnings("resource")
			FileChannel fc = new FileOutputStream(bin_path).getChannel()) {
			fc.write(bb);
			fc.close();
		}
		if (debugLevel > -1) {
			System.out.println("glTFExport(): wrote binary data to "+bin_path);
		}
		// Now generate glTF json
		// prepare mesh short names such as  "img001"
		String [] short_names = new String[total_meshes];
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			String tex_uri = tri_meshes.get(nmesh).getImage();
//			short_names[nmesh] = tex_uri.substring(tex_uri.indexOf("-")+1,tex_uri.lastIndexOf("-"));
			short_names[nmesh] = tex_uri.substring(tex_uri.indexOf("-")+1,tex_uri.lastIndexOf("-"))+"-"+String.format("%03d", nmesh);
		}		
		
		JSONObject gltf_json =  new JSONObject();
		//create asset
		JSONObject gltf_asset = new JSONObject();
		gltf_asset.put("generator", "imagej-elphel");
		gltf_asset.put("copyright", "2022 (C) Elphel, Inc.");
		gltf_asset.put("version",   "2.0");
		gltf_json.put("asset", gltf_asset);
		// set single scene
		gltf_json.put("scene", 0);
		
		// Create nodes -> change to a single node, single mesh with multiple primitives
		JSONArray gltf_nodes = new JSONArray();
		JSONObject gltf_node0 =  new JSONObject();
		gltf_node0.put("name", "all-meshes-node");
		gltf_node0.put("mesh", 0); // index of a single mesh
		gltf_nodes.put(gltf_node0);
		gltf_json.put("nodes", gltf_nodes);
		
		// Create a single scene0
		JSONArray gltf_scene0_nodes = new JSONArray();
		gltf_scene0_nodes.put(0);
		JSONObject gltf_scene0 =  new JSONObject();
		gltf_scene0.put("name",  model_name);
		gltf_scene0.put("nodes", gltf_scene0_nodes);
		
		// Create scenes
		JSONArray gltf_scenes = new JSONArray();
		gltf_scenes.put(gltf_scene0);
		gltf_json.put("scenes", gltf_scenes);
		
		// Create images
		HashMap<String,Integer> map_texture_url = new HashMap<String,Integer>();
		int [] image_index =      new int [total_meshes];   // nmesh first mentioning the url
		boolean [] is_first_ref = new boolean[total_meshes]; // this nmesh is the first to mention		
		JSONArray gltf_images =   new JSONArray();
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			JSONObject gltf_image =  new JSONObject();
			String uri = tri_meshes.get(nmesh).getImage();
			int sz = map_texture_url.size();
			Integer fm = map_texture_url.putIfAbsent(uri,sz);
			if (fm == null) { // did not yet exist
				gltf_image.put("uri", uri);
				gltf_images.put(gltf_image);
				image_index[nmesh] = sz;
				is_first_ref[nmesh] = true;
			} else {
				image_index[nmesh] = fm;
			}
		}
		gltf_json.put("images", gltf_images);
		
		// Create textures
		JSONArray gltf_textures = new JSONArray();
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			if (is_first_ref[nmesh]) {
				JSONObject gltf_texture =  new JSONObject();
				gltf_texture.put("sampler", 0);
				gltf_texture.put("source", image_index[nmesh]); // nmesh);
				gltf_texture.put("name", short_names[nmesh]+"-tex");
				gltf_textures.put(gltf_texture);
			}
		}
		gltf_json.put("textures", gltf_textures);
		
		// Create sampler0
		JSONObject gltf_sampler0 = new JSONObject();
		// page 46 - otherwise needs power-of-two textures
		/*
3.8.4.5. Non-power-of-two Textures
Client implementations SHOULD resize non-power-of-two textures (so that their horizontal and
vertical sizes are powers of two) when running on platforms that have limited support for such
texture dimensions.

Implementation Note
Specifically, if the sampler the texture references:
- has a wrapping mode (either wrapS or wrapT) equal to repeat or mirrored repeat,
or
- has a minification filter (minFilter) that uses mipmapping.
		 */
		gltf_sampler0.put("magFilter", LINEAR);
		gltf_sampler0.put("minFilter", LINEAR); // NEAREST_MIPMAP_LINEAR);
		gltf_sampler0.put("wrapS", CLAMP_TO_EDGE); // REPEAT); // CLAMP_TO_EDGE
		gltf_sampler0.put("wrapT", CLAMP_TO_EDGE); //REPEAT); // CLAMP_TO_EDGE
		// Create samplers
		JSONArray gltf_samplers = new JSONArray();
		gltf_samplers.put(gltf_sampler0);
		gltf_json.put("samplers", gltf_samplers);
		
		// Create materials (material references texture, so each mesh object - new material
		JSONArray gltf_materials = new JSONArray(); 
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			if (is_first_ref[nmesh]) {
				int indx = image_index[nmesh];
				JSONObject gltf_material =             new JSONObject();
				JSONObject gltf_pbrMetallicRoughness = new JSONObject();
				JSONObject gltf_baseColorTexture =     new JSONObject();
				gltf_baseColorTexture.put("index", indx); // nmesh); // reference corresponding texture
				gltf_pbrMetallicRoughness.put("baseColorTexture", gltf_baseColorTexture);
				gltf_pbrMetallicRoughness.put("metallicFactor", 0.0);
				gltf_material.put("pbrMetallicRoughness",gltf_pbrMetallicRoughness);
				gltf_material.put("name", short_names[nmesh]+"-mat");
				if (use_alpha_blend) {
					gltf_material.put("alphaMode", ALPHAMODE_BLEND);
				} else {
					gltf_material.put("alphaMode", ALPHAMODE_OPAQUE);
				}
				if (gltf_emissive) {
					JSONObject gltf_emissiveTexture =        new JSONObject();
					gltf_emissiveTexture.put("index", indx); //  nmesh);
					gltf_material.put("emissiveTexture", gltf_emissiveTexture);
					JSONArray gltf_material_emissiveFactor = new JSONArray();
					gltf_material_emissiveFactor.put(1.0);
					gltf_material_emissiveFactor.put(1.0);
					gltf_material_emissiveFactor.put(1.0);
					gltf_material.put("emissiveFactor", gltf_material_emissiveFactor);
				}
				gltf_materials.put(gltf_material);
			}
		}
		gltf_json.put("materials", gltf_materials);
		
		// create buffer(s) - a single one
		
		// Create buffer0 (only one)
		JSONObject gltf_buffer0 = new JSONObject();
		gltf_buffer0.put("byteLength", buffer_total_bytes);
		gltf_buffer0.put("uri",        bin_name); // relative
		gltf_buffer0.put("name", "single-buffer");
		// Create buffers
		JSONArray gltf_buffers = new JSONArray();
		gltf_buffers.put(gltf_buffer0);
		gltf_json.put("buffers", gltf_buffers);
		
		// create bufferViews - 3 per mesh: index, tex_coords, coordinates
		JSONArray gltf_bufferViews = new JSONArray();
			// bufferView for triangle indices:
		JSONObject gltf_view_index =  new JSONObject();
		gltf_view_index.put("buffer",     0);
		gltf_view_index.put("byteOffset", 0);
		gltf_view_index.put("byteLength", buffer_index_bytes);
///		gltf_view_index.put("byteStride", 2*index_size); // Probably not needed ? page 86
		gltf_view_index.put("target",     ELEMENT_ARRAY_BUFFER);
		gltf_view_index.put("name",       "buffer-view-index");
		gltf_bufferViews.put(gltf_view_index);

		// bufferView for texture coordinates
		JSONObject gltf_view_tex =    new JSONObject();
		gltf_view_tex.put("buffer",     0);
		gltf_view_tex.put("byteOffset", buffer_index_bytes);
		gltf_view_tex.put("byteLength", buffer_texture_bytes);
		gltf_view_tex.put("byteStride", 2*4);
		gltf_view_tex.put("target",     ARRAY_BUFFER);
		gltf_view_tex.put("name",       "buffer-view-tex-coor");
		gltf_bufferViews.put(gltf_view_tex);

		// bufferView for world coordinates
		JSONObject gltf_view_coor =    new JSONObject();
		gltf_view_coor.put("buffer",     0);
		gltf_view_coor.put("byteOffset", buffer_index_bytes+ buffer_texture_bytes);
		gltf_view_coor.put("byteLength", buffer_coordinate_bytes);
		gltf_view_coor.put("byteStride", 3 * 4);
		gltf_view_coor.put("target",     ARRAY_BUFFER);
		gltf_view_coor.put("name",       "buffer-view-world-coor");
		gltf_bufferViews.put(gltf_view_coor);
		gltf_json.put("bufferViews", gltf_bufferViews);
		// generate accessors - 3 per mesh: index, texture coordinates and world coordinates
		JSONArray gltf_accessors = new JSONArray();
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			// accessor for triangle indices
			JSONObject gltf_accessor_index =  new JSONObject();
			JSONArray gltf_accessor_index_min = new JSONArray();
			JSONArray gltf_accessor_index_max = new JSONArray();
			gltf_accessor_index_min.put(minmax_indx[nmesh][0]);
			gltf_accessor_index_max.put(minmax_indx[nmesh][1]);
			gltf_accessor_index.put("bufferView", 0); // of 3
			gltf_accessor_index.put("byteOffset", offset_len_indx[nmesh][0]);
			// actually it is signed short/signed int
			gltf_accessor_index.put("componentType", (index_size == 2)? COMPONENT_UNSIGNED_SHORT : COMPONENT_UNSIGNED_INT);
			gltf_accessor_index.put("count", 3 * mesh_tri[nmesh]);
			gltf_accessor_index.put("min",  gltf_accessor_index_min);
			gltf_accessor_index.put("max",  gltf_accessor_index_max);
			gltf_accessor_index.put("type", TYPE_SCALAR);
			gltf_accessor_index.put("name", short_names[nmesh]+"-accessor-index-"+(3*nmesh+0));
			gltf_accessors.put(gltf_accessor_index);

			// accessor for texture coordinates
			JSONObject gltf_accessor_texcoor =  new JSONObject();
			JSONArray gltf_accessor_texcoor_min = new JSONArray();
			JSONArray gltf_accessor_texcoor_max = new JSONArray();
			gltf_accessor_texcoor_min.put(minmax_tex[nmesh][0][0]);
			gltf_accessor_texcoor_min.put(minmax_tex[nmesh][1][0]);
			gltf_accessor_texcoor_max.put(minmax_tex[nmesh][0][1]);
			gltf_accessor_texcoor_max.put(minmax_tex[nmesh][1][1]);
			gltf_accessor_texcoor.put("bufferView", 1); // of 3
			gltf_accessor_texcoor.put("byteOffset", offset_len_tex[nmesh][0]);
			// actually it is signed short/signed int
			gltf_accessor_texcoor.put("componentType", COMPONENT_FLOAT);
			gltf_accessor_texcoor.put("count", mesh_vert[nmesh]);
			gltf_accessor_texcoor.put("min",  gltf_accessor_texcoor_min);
			gltf_accessor_texcoor.put("max",  gltf_accessor_texcoor_max);
			gltf_accessor_texcoor.put("type", TYPE_VEC2);
			gltf_accessor_texcoor.put("name", short_names[nmesh]+"-accessor-texcoor-"+(3*nmesh+1));
			gltf_accessors.put(gltf_accessor_texcoor);
			
			// accessor for world coordinates
			JSONObject gltf_accessor_coor =  new JSONObject();
			JSONArray gltf_accessor_coor_min = new JSONArray();
			JSONArray gltf_accessor_coor_max = new JSONArray();
			gltf_accessor_coor_min.put(minmax_coor[nmesh][0][0]);
			gltf_accessor_coor_min.put(minmax_coor[nmesh][1][0]);
			gltf_accessor_coor_min.put(minmax_coor[nmesh][2][0]);
			gltf_accessor_coor_max.put(minmax_coor[nmesh][0][1]);
			gltf_accessor_coor_max.put(minmax_coor[nmesh][1][1]);
			gltf_accessor_coor_max.put(minmax_coor[nmesh][2][1]);
			gltf_accessor_coor.put("bufferView", 2); // of 3
			gltf_accessor_coor.put("byteOffset", offset_len_coor[nmesh][0]);
			// actually it is signed short/signed int
			gltf_accessor_coor.put("componentType", COMPONENT_FLOAT);
			gltf_accessor_coor.put("count", mesh_vert[nmesh]);
			gltf_accessor_coor.put("min",  gltf_accessor_coor_min);
			gltf_accessor_coor.put("max",  gltf_accessor_coor_max);
			gltf_accessor_coor.put("type", TYPE_VEC3);
			gltf_accessor_coor.put("name", short_names[nmesh]+"-accessor-coor-"+(3*nmesh+2));
			gltf_accessors.put(gltf_accessor_coor);
		}
		gltf_json.put("accessors", gltf_accessors);
		
		// Generate meshes - trying a single mesh with multiple primitives
		JSONObject gltf_mesh0 =  new JSONObject();
		JSONArray  gltf_mesh0_primitives = new JSONArray();
		for (int nmesh = 0; nmesh < total_meshes; nmesh++) {
			JSONObject gltf_mesh0_primitive =       new JSONObject(); // does not have "name"
			JSONObject gltf_mesh0_primitive_attr =  new JSONObject();
			gltf_mesh0_primitive_attr.put("POSITION",   3 * nmesh + 2);
			gltf_mesh0_primitive_attr.put("TEXCOORD_0", 3 * nmesh + 1);
			gltf_mesh0_primitive.put("attributes",      gltf_mesh0_primitive_attr);
			gltf_mesh0_primitive.put("indices",         3 * nmesh + 0);
			gltf_mesh0_primitive.put("mode",            MODE_TRIANGLES);
			gltf_mesh0_primitive.put("material",        image_index[nmesh]); // nmesh);
			gltf_mesh0_primitives.put(gltf_mesh0_primitive);
		}
		gltf_mesh0.put("primitives", gltf_mesh0_primitives);
		gltf_mesh0.put("name", "all-model-meshes");
		JSONArray  gltf_meshes = new JSONArray();
		gltf_meshes.put(gltf_mesh0); // just a single mesh with multiple primitives
		gltf_json.put("meshes", gltf_meshes);
 
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		
		JsonElement je = JsonParser.parseString(gltf_json.toString());
		String prettyJsonString = gson.toJson(je);
		FileWriter gltf_writer; // f0 = new FileWriter("output.txt");
		gltf_writer = new FileWriter(gltf_path);
		gltf_writer.write(prettyJsonString);
		gltf_writer.close();
		return;
	}
}
