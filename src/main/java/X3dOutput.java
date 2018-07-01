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
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.attribute.PosixFilePermission;
import java.util.ArrayList;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

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

	public X3dOutput() {}

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
        if (use_backdrop && (bgnd_pass.texture != null)) {
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
			System.out.println("x3d file is saved to "+path);

		} catch (TransformerException tfe) {
			tfe.printStackTrace();
		}

	}

	public double [][] getBBox() // center: x,y,z, size:x,y,z
	{
//		double depth =  geometry_correction.getZFromDisparity(clt_parameters.bgnd_range);
		double depth =  clt_parameters.infinityDistance;
		double width  = depth * geometry_correction.getFOVWidth();
		double height = depth * geometry_correction.getFOVHeight();
		double [][] bbox = {{0, 0, -depth/2},{width,height,depth}};
		return bbox;
	}

	public void generateKML(
			String path,
			boolean overwrite,
			String icon_path, //<href>x3d/1487451413_967079.x3d</href> ?
			double timestamp,
			double [] lla)
	{
		if (!overwrite) {
			if (new File(path).exists()) {
				System.out.println("KML file "+path+" exists, skipping KML generation");
				return;
			}
		}
		Document kmlDoc;
		try {

			DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder docBuilder =        docFactory.newDocumentBuilder();
			kmlDoc =                   docBuilder.newDocument();

		} catch (ParserConfigurationException pce) {
			pce.printStackTrace();
			return;
		}
		// root elements
		Element el_kml = kmlDoc.createElement("kml");
		el_kml.setAttribute("xmlns","http://earth.google.com/kml/2.2");
		kmlDoc.appendChild(el_kml);
		Element el_Document = addElement(kmlDoc, el_kml,"Document");
		Element el_PhotoOverlay = addElement(kmlDoc, el_Document, "PhotoOverlay"); // kmlDoc.createElement("PhotoOverlay");
		addTextElement(kmlDoc, el_PhotoOverlay, "name", "test");
		addTextElement(kmlDoc, el_PhotoOverlay, "description", "this kml description");
		addTextElement(kmlDoc, el_PhotoOverlay, "visibility", "1");
		addTextElement(kmlDoc, el_PhotoOverlay, "shape", "rectangle");


		Element el_TimeStamp = addElement(kmlDoc, el_PhotoOverlay, "TimeStamp");
		addTextElement(kmlDoc, el_TimeStamp, "when", String.format("%f",timestamp));

		Element el_Camera = addElement(kmlDoc, el_PhotoOverlay, "Camera");
		addTextElement(kmlDoc, el_Camera, "latitude",  String.format("%f",lla[0]));
		addTextElement(kmlDoc, el_Camera, "longitude", String.format("%f",lla[1]));
		addTextElement(kmlDoc, el_Camera, "altitude",  String.format("%f",lla[2]));
		addTextElement(kmlDoc, el_Camera, "heading",   String.format("%f",0.0));
		addTextElement(kmlDoc, el_Camera, "tilt",      String.format("%f",90.0));
		addTextElement(kmlDoc, el_Camera, "roll",      String.format("%f",0.0));
		if (icon_path != null) {
			Element el_Icon = addElement(kmlDoc, el_PhotoOverlay, "Icon");
			addTextElement(kmlDoc, el_Icon, "href",  icon_path);
		}
		Element el_ExtendedData = addElement(kmlDoc, el_PhotoOverlay, "ExtendedData");
		Element el_OriginalData = addElement(kmlDoc, el_ExtendedData, "OriginalData");
		addTextElement(kmlDoc, el_OriginalData, "latitude",  String.format("%f",lla[0]));
		addTextElement(kmlDoc, el_OriginalData, "longitude", String.format("%f",lla[1]));
		addTextElement(kmlDoc, el_OriginalData, "altitude",  String.format("%f",lla[2]));
		addTextElement(kmlDoc, el_OriginalData, "heading",   String.format("%f",0.0));
		addTextElement(kmlDoc, el_OriginalData, "tilt",      String.format("%f",90.0));
		addTextElement(kmlDoc, el_OriginalData, "roll",      String.format("%f",0.0));

		try {
			// write the content into xml file
			TransformerFactory transformerFactory = TransformerFactory.newInstance();
			Transformer transformer = transformerFactory.newTransformer();
			DOMSource source = new DOMSource(kmlDoc);
			StreamResult result = new StreamResult(new File(path));
			// Output to console for testing
			// StreamResult result = new StreamResult(System.out);
			transformer.setOutputProperty(OutputKeys.INDENT, "yes");
			transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
			transformer.transform(source, result);
			System.out.println("kml file is saved to "+path);

		} catch (TransformerException tfe) {
			tfe.printStackTrace();
		}

		try {
			Path fpath = Paths.get((new File(path)).getCanonicalPath());
			Set<PosixFilePermission> perms =  Files.getPosixFilePermissions(fpath);
			perms.add(PosixFilePermission.OTHERS_WRITE);
			Files.setPosixFilePermissions(fpath, perms);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}



	}


	Element addTextElement(Document doc, Element parent, String name, String text)
	{
		Element child = doc.createElement(name);
		child.setTextContent(text);
		parent.appendChild(child);
		return child;
	}
	Element addElement(Document doc, Element parent, String name)
	{
		Element child = doc.createElement(name);
		parent.appendChild(child);
		return child;

	}


}


/*
 * <?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.2">
<Document>
<PhotoOverlay>
	<name>test</name>
	<visibility>1</visibility>
	<shape>rectangle</shape>
	<TimeStamp>
		<when>1487451413.967079</when>
	</TimeStamp>
	<Camera>
		<longitude>-111.93292339</longitude>
		<latitude>40.72350154</latitude>
		<altitude>1305.1</altitude>
		<heading>68.8186</heading>
		<tilt>90</tilt>
		<roll>0</roll>
	</Camera>
	<Icon>
		<href>x3d/1487451413_967079.x3d</href>
	</Icon>
	<ExtendedData>
		<OriginalData>
			<longitude>-111.9328843</longitude>
			<latitude>40.7233861</latitude>
			<altitude>1305.1</altitude>
			<heading>65</heading>
			<tilt>90</tilt>
			<roll>0</roll>
		</OriginalData>
	</ExtendedData>
<description></description></PhotoOverlay>
</Document>
</kml>

<?xml version="1.0" encoding="UTF-8"?>
<x3d profile="Interchange" version="3.3">
  <Scene>
    <Group bboxCenter="0.000 0.000 -1419.899 " bboxSize="3704.001 2766.569 2839.798 " class="GroupTop" id="GroupTop"/>

 */
