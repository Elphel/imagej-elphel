/**
 ** -----------------------------------------------------------------------------**
 ** ImagejJp4Tiff.java
 **
 ** Uses loci.format compatible readers for Elphel 8/16 bpp monochrome Tiff and
 ** JP4 files to read/parse camera files in a single url read operation by
 ** buffering camera data with Location.mapFile()
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ImagejJp4Tiff.java is free software: you can redistribute it and/or modify
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
package com.elphel.imagej.readers;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

import org.apache.commons.compress.utils.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import loci.common.ByteArrayHandle;
import loci.common.Location;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.ClassList;
import loci.formats.FormatException;
import loci.formats.IFormatReader;
import loci.formats.ImageReader;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;

public class ImagejJp4Tiff {

	// -- Constants --

	private static final Logger LOGGER = LoggerFactory.getLogger(ClassList.class);
	private static final boolean BYPASS_SERVICES = true;

	// -- Fields --

	private ImageReader reader = null;
	private String content_fileName = "undefined"; // from Content-disposition
	private URL url = null; // save here actual URL when reading file to memory
	IMetadata omeMeta = null;

	// -- Constructor --

	public ImagejJp4Tiff() {
		// Bypass readers.txt
		ClassList<IFormatReader> classList =  new ClassList<IFormatReader>(IFormatReader.class);
//		classList.addClass(com.elphel.imagej.readers.ElphelTiffJp4Reader.class);
		classList.addClass(com.elphel.imagej.readers.ElphelJp4Reader.class);
		classList.addClass(com.elphel.imagej.readers.ElphelTiffReader.class);

		if (!BYPASS_SERVICES) {
			ServiceFactory factory = null;

			try {
				factory = new ServiceFactory();
			} catch (DependencyException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			OMEXMLService service = null;
			try {
				service = factory.getInstance(OMEXMLService.class);
			} catch (DependencyException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			try {
				omeMeta = service.createOMEXMLMetadata();
			} catch (ServiceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
	    reader = new ImageReader(classList);
	    if (!BYPASS_SERVICES) {
	    	reader.setMetadataStore(omeMeta);
	    }
	}

	// -- API methods --

	public ImagePlus readTiffJp4(String path) throws IOException, FormatException {
		return readTiffJp4(path, 	"STD_"); // null);
	}

	public ImagePlus readTiffJp4(String path_url, String std ) throws IOException, FormatException { // std - include non-elphel properties with prefix std
		// determine if it is a file or URL and read url to memory
		// If URL, then read to memory, if normal file - use direct access
		url = null;
		//	 		String mime = null; // use to select jp4/tiff later? Or to check it is correct
		content_fileName = null;
		try {
			url = new URL(path_url);
		} catch (MalformedURLException e) {
			LOGGER.warn("Bad URL: " + path_url);
		}
		if (url != null) {
			LOGGER.info("Read "+ path_url +" to memory first");
			URLConnection connection = url.openConnection();

			String content_disposition = connection.getHeaderField("Content-Disposition");
			// raw = "attachment; filename=abc.jpg"
			if(content_disposition != null && content_disposition.indexOf("=") != -1) {
				content_fileName = content_disposition.split("=")[1]; //getting value after '='
				// trim quotes
				content_fileName= content_fileName.substring(1, content_fileName.length()-1);
			} else {
				String mime =  connection.getContentType();
				int slash = mime.lastIndexOf("/");
				String suffix = slash < 0 ? "" : mime.substring(slash+1);
				content_fileName = "unknown." + suffix;
			}
			InputStream is =    url.openStream (); //
			byte[] inBytes = IOUtils.toByteArray(is);
			if (is != null) is.close();
			LOGGER.info("Bytes read: "+ inBytes.length);
			Location.mapFile(content_fileName, new ByteArrayHandle(inBytes));
			//				  HashMap<String,Object> dbg_loc = Location.getIdMap();
			//				super.setId(content_fileName);
		} else { // read file normally
			content_fileName = path_url;
			LOGGER.info("read '"+path_url+"' file directly");
		}
		reader.setId(content_fileName);
		byte [] bytes = null;
		ImagePlus imp=  null;
		bytes = reader.openBytes(0);
		int bpp =   reader.getBitsPerPixel();

		boolean is_le = reader.isLittleEndian();
		int bytes_per_pixel =  (bpp + 7) / 9;
		float [] pixels = new float [bytes.length/bytes_per_pixel];
		if (bytes_per_pixel == 1) {
			for (int i = 0; i < pixels.length; i++) {
				pixels[i] = ((bytes[i])) & 0xff;
			}
		} else {
			ByteBuffer bb = ByteBuffer.wrap(bytes);
			if (is_le) {
				bb.order( ByteOrder.LITTLE_ENDIAN);
			} else {
				bb.order( ByteOrder.BIG_ENDIAN);
			}
			for (int i = 0; i < pixels.length; i++) {
				pixels[i] = ((bb.getShort())) & 0xffff;
			}
		}

		ImageProcessor ip=new FloatProcessor(reader.getSizeX(), reader.getSizeY());
		ip.setPixels(pixels);
		ip.resetMinAndMax();
		Hashtable<String, Object> meta_hash = reader.getGlobalMetadata();
		String prefix = ElphelTiffReader.ELPHEL_PROPERTY_PREFIX;
		String imageName = content_fileName; // path;
		String imageNameKey = prefix+ElphelTiffReader.CONTENT_FILENAME;
		if (meta_hash.containsKey(imageNameKey)) {
			imageName = meta_hash.get(imageNameKey).toString();
		}
		imp =  new ImagePlus(imageName, ip); // original jp46 reader had full path as title

		// first - save all as properties, later - only ELPHEL_*
		for (String key:meta_hash.keySet()) {
			if (key.startsWith(prefix)) {
				imp.setProperty(key.substring(prefix.length()), meta_hash.get(key).toString());
			} else if (std != null) {
				imp.setProperty(std+(key.replace(" ","_").replace("/","_")), meta_hash.get(key).toString());
			}
		}
		encodeProperiesToInfo(imp);
		Location.mapFile(content_fileName, null);
		return imp;
	}

	// -- Helper methods --

	public static ImagePlus encodeProperiesToInfo(ImagePlus imp){
		String info="<?xml version=\"1.0\" encoding=\"UTF-8\"?><properties>";
		Set<Object> jp4_set;
		Properties jp4_prop;
		Iterator<Object> itr;
		String str;
		jp4_prop=imp.getProperties();
		if (jp4_prop!=null) {
			jp4_set=jp4_prop.keySet();
			itr=jp4_set.iterator();
			while(itr.hasNext()) {
				str = (String) itr.next();
				//				if (!str.equals("Info")) info+="<"+str+">\""+jp4_prop.getProperty(str)+"\"</"+str+">";
				if (!str.equals("Info")) info+="<"+str+">"+jp4_prop.getProperty(str)+"</"+str+">";
			}
		}
		info+="</properties>\n";
		imp.setProperty("Info", info);
		return imp;
	}


}
