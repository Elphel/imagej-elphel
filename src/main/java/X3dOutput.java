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
// Will use 1m units

public class X3dOutput {
	GeometryCorrection                                     geometry_correction;
	public ArrayList <TileProcessor.CLTPass3d>             clt_3d_passes;
	public  EyesisCorrectionParameters.CLTParameters       clt_parameters;
	public EyesisCorrectionParameters.CorrectionParameters correctionsParameters;
	public int debugLevel = 1;
	Document    x3dDoc;
	Element     el_X3d;
	Element     el_Scene;
	Element     el_TopGroup;

	public X3dOutput(
			EyesisCorrectionParameters.CLTParameters        clt_parameters,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters,
			GeometryCorrection                              geometry_correction,
			ArrayList <TileProcessor.CLTPass3d>             clt_3d_passes){
		this.clt_parameters =        clt_parameters;
		this.correctionsParameters = correctionsParameters;
		this.geometry_correction =    geometry_correction;
		this.clt_3d_passes = clt_3d_passes;
	}
	// init document, bounding box, backdrop
	public void generateBackground()
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

		TileProcessor.CLTPass3d bgnd_pass = clt_3d_passes.get(0);

        Element el_Bgnd = x3dDoc.createElement("Background");
        el_Bgnd.setAttribute("class","Background");
        el_Bgnd.setAttribute("id",   "Background");
        el_Bgnd.setAttribute("frontUrl",   bgnd_pass.texture);
        // temporarily - add same picture to all other sides. Actually - any square will work, make some
        // perspective grids/ colors to simplify orientation when looking wrong way
        el_Bgnd.setAttribute("backUrl",   bgnd_pass.texture);
        el_Bgnd.setAttribute("leftUrl",   bgnd_pass.texture);
        el_Bgnd.setAttribute("rightUrl",   bgnd_pass.texture);
        el_Bgnd.setAttribute("topUrl",   bgnd_pass.texture);
        el_Bgnd.setAttribute("bottomUrl",   bgnd_pass.texture);
		el_Scene.appendChild(el_Bgnd);
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
		double [][] bbox = {{0, 0, depth/2},{width,height,depth}};
		return bbox;
	}
}
