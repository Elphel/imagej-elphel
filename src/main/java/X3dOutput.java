/**
 **
 ** X3dOutput - generate x3d representation of the model
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  X3dOutput.java is free software: you can redistribute it and/or modify
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

import java.io.File;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Attr;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

// Will use 1m units

public class X3dOutput {
	GeometryCorrection                                     geometry_correction;
	public ArrayList <CLTPass3d>             clt_3d_passes;
	public  EyesisCorrectionParameters.CLTParameters       clt_parameters;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters;
	public int debugLevel = 1;
	Document    x3dDoc;
	Element     el_X3d;
	Element     el_Scene;
	Element     el_TopGroup;
	public int  max_line_length = 0; // 100; // 0 - no limit

	public X3dOutput(
			EyesisCorrectionParameters.CLTParameters        clt_parameters,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters,
			GeometryCorrection                              geometry_correction,
			ArrayList <CLTPass3d>             clt_3d_passes){
		this.clt_parameters =        clt_parameters;
		this.correctionsParameters = correctionsParameters;
		this.geometry_correction =    geometry_correction;
		this.clt_3d_passes = clt_3d_passes;
	}
	// init document, bounding box, backdrop
	public void generateBackground(boolean use_backdrop)
	{
		try {

			DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
			x3dDoc = docBuilder.newDocument();

		} catch (ParserConfigurationException pce) {
			pce.printStackTrace();
		}
		// root elements
		el_X3d = x3dDoc.createElement("x3d");
		el_X3d.setAttribute("profile","Interchange");
		el_X3d.setAttribute("version","3.3");
		x3dDoc.appendChild(el_X3d);
		el_Scene  = x3dDoc.createElement("Scene");
		el_X3d.appendChild(el_Scene);
		
		el_TopGroup = x3dDoc.createElement("Group");
		double [][] bbox = getBBox();
		el_TopGroup.setAttribute("class","GroupTop");
		el_TopGroup.setAttribute("id",   "GroupTop");
		
		el_TopGroup.setAttribute("bboxCenter",  String.format("%.3f %.3f %.3f ",bbox[0][0],bbox[0][1],bbox[0][2]));
		el_TopGroup.setAttribute("bboxSize",    String.format("%.3f %.3f %.3f ",bbox[1][0],bbox[1][1],bbox[1][2]));
		
		el_Scene.appendChild(el_TopGroup);

		CLTPass3d bgnd_pass = clt_3d_passes.get(0);

        Element el_Bgnd = x3dDoc.createElement("Background");
        el_Bgnd.setAttribute("class","Background");
        el_Bgnd.setAttribute("id",   "Background");
        if (use_backdrop) {
        	el_Bgnd.setAttribute("frontUrl",   bgnd_pass.texture);
        	// temporarily - add same picture to all other sides. Actually - any square will work, make some
        	// perspective grids/ colors to simplify orientation when looking wrong way
        	el_Bgnd.setAttribute("backUrl",   bgnd_pass.texture);
        	el_Bgnd.setAttribute("leftUrl",   bgnd_pass.texture);
        	el_Bgnd.setAttribute("rightUrl",   bgnd_pass.texture);
        	el_Bgnd.setAttribute("topUrl",   bgnd_pass.texture);
        	el_Bgnd.setAttribute("bottomUrl",   bgnd_pass.texture);
        }
		el_Scene.appendChild(el_Bgnd);
	}
	
	public void addCluster(
			String      url,
			String      id,
			String      class_name,
			double [][] texCoord,
			double [][] coordinate,
			int [][] triangles)
	{
		//	public int  max_line_length = 100; // 0 - no limit
		int linepos = 0;
		StringBuffer sb_coord_index = new StringBuffer();
		for (int i = 0; i < triangles.length; i++){
			if ((max_line_length > 0) && ((sb_coord_index.length()-linepos) >= max_line_length)){
				sb_coord_index.append("\n");
				linepos = 0;
			} else if ((linepos > 0) || ((max_line_length == 0) && (sb_coord_index.length() > 0))){
				sb_coord_index.append(" ");
			}
//			sb_coord_index.append(triangles[i][0]+" "+triangles[i][1]+" "+triangles[i][2]+" -1");
			sb_coord_index.append(triangles[i][2]+" "+triangles[i][1]+" "+triangles[i][0]+" -1");
		}
		StringBuffer sb_coords = new StringBuffer();
		linepos = 0;
		for (int i = 0; i < coordinate.length; i++){
			if ((max_line_length >0) && ((sb_coords.length()-linepos) >= max_line_length)){
				sb_coords.append("\n");
				linepos = 0;
			} else if ((linepos > 0) || ((max_line_length == 0) && (sb_coords.length() > 0))){
				sb_coords.append(" ");
			}
			sb_coords.append(String.format("%.3f %.3f %.3f", coordinate[i][0], coordinate[i][1], coordinate[i][2]));
		}
		StringBuffer sb_tex_coords = new StringBuffer();
		linepos = 0;
		for (int i = 0; i < texCoord.length; i++){
			if ((max_line_length >0) && ((sb_tex_coords.length()-linepos) >= max_line_length)){
				sb_tex_coords.append("\n");
				linepos = 0;
			} else if ((linepos > 0) || ((max_line_length == 0) && (sb_tex_coords.length() > 0))){
				sb_tex_coords.append(" ");
			}
			sb_tex_coords.append(String.format("%.4f %.4f", texCoord[i][0], texCoord[i][1]));
		}
		String sindex =  sb_coord_index.toString(); // for both coordIndex and texCoordIndex
		String scoord =  sb_coords.toString();
		String stcoord = sb_tex_coords.toString();
		
        Element el_shape =        x3dDoc.createElement("Shape");
		el_Scene.appendChild(el_shape);
		if (id != null) el_shape.setAttribute("id",id);
		if (class_name != null) el_shape.setAttribute("class",class_name);
		
        Element el_appearance =   x3dDoc.createElement("Appearance");
		el_shape.appendChild(el_appearance);

/*
        Element el_material =   x3dDoc.createElement("Material");
		el_appearance.appendChild(el_material);
		el_material.setAttribute("diffuseColor",    "0.376471 0.5 0.376471");
*/		
		
        Element el_imageTexture = x3dDoc.createElement("ImageTexture");
        el_imageTexture.setAttribute("url",url);
        el_appearance.appendChild(el_imageTexture);
        
        Element el_IFC =          x3dDoc.createElement("IndexedFaceSet");
        el_IFC.setAttribute("coordIndex",    sindex); // can it be reused?
        el_IFC.setAttribute("texCoordIndex", sindex);
        el_shape.appendChild(el_IFC);
        
        Element el_coordinate =   x3dDoc.createElement("Coordinate");
        el_coordinate.setAttribute("point",    scoord);
        el_IFC.appendChild(el_coordinate);

        Element el_texCoordinate =   x3dDoc.createElement("TextureCoordinate");
        el_texCoordinate.setAttribute("point",    stcoord);
        el_IFC.appendChild(el_texCoordinate);
        
}
	
	
	
	
	// close document, generate x3d file
	public void generateX3D(String path)
	{
		try {
			// write the content into xml file
			TransformerFactory transformerFactory = TransformerFactory.newInstance();
			Transformer transformer = transformerFactory.newTransformer();
			DOMSource source = new DOMSource(x3dDoc);
			StreamResult result = new StreamResult(new File(path));
			// Output to console for testing
			// StreamResult result = new StreamResult(System.out);
			transformer.setOutputProperty(OutputKeys.INDENT, "yes");
			transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");			
			transformer.transform(source, result);
			System.out.println("x3d file saved to "+path);

		} catch (TransformerException tfe) {
			tfe.printStackTrace();
		}

	}

	public double [][] getBBox() // center: x,y,z, size:x,y,z
	{
		double depth =  geometry_correction.getZFromDisparity(clt_parameters.bgnd_range); 
		double width  = depth * geometry_correction.getFOVWidth();
		double height = depth * geometry_correction.getFOVHeight();
		double [][] bbox = {{0, 0, -depth/2},{width,height,depth}};
		return bbox;
	}
}
