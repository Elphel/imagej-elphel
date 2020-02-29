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
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

import org.apache.commons.compress.utils.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ij.ImagePlus;
import ij.io.FileInfo;
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
	private static final boolean BYPASS_SERVICES = false; //  true;
	private static final String TELEMETRY_PREFIX = "TLM_";
	private static final String SERVICES_PATH = "services.properties.forelphel";
	private static final boolean KEEP_EXTENSION = false; // remove extension for ImagePlus title
	private static final boolean BOTTOM_TELEMETRY = true; // telemetry (if present) is at the bottom of the image

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
		// Trying ServiceFactory before it is going to be initialized, so static defaultFactory will be initialized
		// with small set of services - only needed for Elphel
		//URL u = this.getClass().getResource(SERVICES_PATH);
//		URL u = this.getClass().getResource("/"+SERVICES_PATH);

		if (!BYPASS_SERVICES) {
			ServiceFactory factory = null;

			try {
//				factory = new ServiceFactory();
				// Still does not work - during initFile->populatePixels it loads all
				// services, including buggy ones

			factory = new ServiceFactory(); // "/"+SERVICES_PATH);

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
		return readTiffJp4(path, true); // null);
	}
	public ImagePlus readTiffJp4(String path, boolean scale) throws IOException, FormatException {
		return readTiffJp4(path, scale,	"STD_"); // null);
	}

	public ImagePlus readTiffJp4(String path_url, boolean scale, String std ) throws IOException, FormatException { // std - include non-elphel properties with prefix std
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
		//https://stackoverflow.com/questions/39086500/read-http-response-header-and-body-from-one-http-request-in-java
		if (url != null) {
			LOGGER.info("Read "+ path_url +" to memory first");
			URLConnection connection = url.openConnection();
// Wrong - waits forever
			String content_disposition = connection.getHeaderField("Content-Disposition"); // reads file
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
//			InputStream is =    url.openStream (); // reads file for the second time
			InputStream is =    connection.getInputStream ();
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
		ByteBuffer bb = null;
		if (bytes_per_pixel == 1) {
			for (int i = 0; i < pixels.length; i++) {
				pixels[i] = ((bytes[i])) & 0xff;
			}
		} else {
			bb = ByteBuffer.wrap(bytes);
			if (is_le) {
				bb.order( ByteOrder.LITTLE_ENDIAN);
			} else {
				bb.order( ByteOrder.BIG_ENDIAN);
			}
			for (int i = 0; i < pixels.length; i++) {
				pixels[i] = ((bb.getShort())) & 0xffff;
			}
		}
		Hashtable<String, Object> meta_hash = reader.getGlobalMetadata();

		boolean degamma = bytes_per_pixel < 2; // both JP4 and 8-bit tiff
		if (deGammaScale(pixels, reader.getSizeX(), meta_hash, degamma, scale ) == null) {
			LOGGER.error("Problem degamma/scaling of "+content_fileName);

		}
		boolean telemetry = (bytes_per_pixel == 2);
		Lepton3Telemetry lepton3Telemetry = null;
		if (telemetry) {
			lepton3Telemetry = new Lepton3Telemetry(bb,pixels, true); // bottom
			telemetry = lepton3Telemetry.hasTelemetry();
		}
		ImageProcessor ip = null;
		if (telemetry) {
			ip=new FloatProcessor(Lepton3Telemetry.LEPTON3_IMAGE_WIDTH, Lepton3Telemetry.LEPTON3_IMAGE_HEIGHT);
			ip.setPixels(lepton3Telemetry.getPixels()); // bo
		} else {
			ip=new FloatProcessor(reader.getSizeX(), reader.getSizeY());
			ip.setPixels(pixels);
		}
		ip.resetMinAndMax();
		String prefix = ElphelTiffReader.ELPHEL_PROPERTY_PREFIX;
		String imageName = content_fileName; // path;
		String imageNameKey = prefix+ElphelTiffReader.CONTENT_FILENAME;
		if (meta_hash.containsKey(imageNameKey)) {
			imageName = meta_hash.get(imageNameKey).toString();
		}
		if (!KEEP_EXTENSION) {
			int dot_indx = imageName.lastIndexOf(".");
			if (dot_indx >= 0) {
				imageName = imageName.substring(0, dot_indx);
			}
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
		if (telemetry) {
			HashMap<String, String> telemetryMap = lepton3Telemetry.parseTelemetry();
			for (String key:telemetryMap.keySet()) {
				imp.setProperty(TELEMETRY_PREFIX+key, telemetryMap.get(key));
			}
		}
		
		FileInfo fi = imp.getFileInfo();
		
		String dir = imageName.substring(0, imageName.lastIndexOf("/")+1);
		
		fi.directory = dir;
		fi.fileName = imageName;
		imp.setFileInfo(fi);
		
		encodeProperiesToInfo(imp);
		
		System.out.println("");
		
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

	public float [] deGammaScale(float [] pixels, int width, Hashtable<String, Object> meta_hash, boolean degamma, boolean scale) {
		int height = pixels.length/width;
		String prefix = ElphelTiffReader.ELPHEL_PROPERTY_PREFIX;
		double [] rgains =        {1.0, 1.0, 1.0, 1.0};
		double [] blacks =        new double[4];
		double[] blacks256=       new double[4];
		double [] gammas =        new double[4];
		long [] gamma_scales =    new long[4];
		double [][] rgammas =     new double[4][];
		double min_gain = 1.0;
		boolean flipv = false, fliph=false;
		if (meta_hash.get(prefix+"FLIPV")!=null) flipv=   Integer.valueOf((String) meta_hash.get(prefix+"FLIPV")).intValue() > 0; // else return null;
		if (meta_hash.get(prefix+"FLIPH")!=null) fliph=    Integer.valueOf((String) meta_hash.get(prefix+"FLIPV")).intValue() > 0;// else return null;
		if (meta_hash.get(prefix+"GAIN")!=null) min_gain=  Double.valueOf((String) meta_hash.get(prefix+"GAIN")).doubleValue();   // else return null;
		for (int i=0;i<4;i++) { /* r,g,gb,b */
			if (scale) {
				if (meta_hash.get(prefix+"gains_"+i)!=null)        rgains[i]=      min_gain/Double.valueOf((String) meta_hash.get(prefix+"gains_"+i)).doubleValue();        else return null;
			}
			if (degamma) {
				if (meta_hash.get(prefix+"blacks_"+i)!=null)       blacks[i]=        Double.valueOf((String) meta_hash.get(prefix+"blacks_"+i)).doubleValue();       else return null;
				if (meta_hash.get(prefix+"gammas_"+i)!=null)       gammas[i]=        Double.valueOf((String) meta_hash.get(prefix+"gammas_"+i)).doubleValue();       else return null;
				if (meta_hash.get(prefix+"gamma_scales_"+i)!=null) gamma_scales[i]=  Integer.valueOf((String) meta_hash.get(prefix+"gamma_scales_"+i)).intValue();   else return null;
			}
		}
		if (flipv) {
			if (scale) {
				ElphelMeta.swapArrayElements (rgains,       1, 3);
				ElphelMeta.swapArrayElements (rgains,       0, 2);
			}
			if (degamma) {
				ElphelMeta.swapArrayElements (blacks,      1, 3);
				ElphelMeta.swapArrayElements (blacks,      0, 2);
				ElphelMeta.swapArrayElements (gammas,      1, 3);
				ElphelMeta.swapArrayElements (gammas,      0, 2);
				ElphelMeta.swapArrayElements (gamma_scales,1, 3);
				ElphelMeta.swapArrayElements (gamma_scales,0, 2);
			}
		}
		if (fliph) {
			if (scale) {
				ElphelMeta.swapArrayElements (rgains,       1, 0);
				ElphelMeta.swapArrayElements (rgains,       3, 2);
			}
			if (degamma) {
				ElphelMeta.swapArrayElements (blacks,      1, 0);
				ElphelMeta.swapArrayElements (blacks,      3, 2);
				ElphelMeta.swapArrayElements (gammas,      1, 0);
				ElphelMeta.swapArrayElements (gammas,      3, 2);
				ElphelMeta.swapArrayElements (gamma_scales,1, 0);
				ElphelMeta.swapArrayElements (gamma_scales,3, 2);
			}
		}
		if (degamma) {
			for (int i=0;i<4;i++) {
				rgammas[i]=ElphelMeta.elphel_gamma_calc (gammas[i], blacks[i], gamma_scales[i]);
				blacks256[i]=256.0*blacks[i];
				for (int j = 0; j < rgammas[i].length; j++) {
					rgammas[i][j] -= blacks256[i];
					if (scale) {
						rgammas[i][j] *= rgains[i];
					}
				}
			}
		}
		if (degamma) {
			for (int y = 0; y < height; y+=2) {
				for (int dy = 0; dy < 2; dy++) {
					int base = (y + dy)*width;
					for (int x = 0; x < width; x+=2) {
						for (int dx = 0; dx<2; dx++) {
							int indx = base + x + dx;
							int iv = (int) pixels[indx];
							if (iv < 0) iv = 0;
							else if (iv > 255) iv = 255;
							pixels[indx] = (float) rgammas[(dy << 1) + 1 - dx][iv];
						}
					}
				}
			}
		} else if (scale) {
			for (int y = 0; y < height; y+=2) {
				for (int dy = 0; dy < 2; dy++) {
					int base = (y + dy)*width;
					for (int x = 0; x < width; x+=2) {
						for (int dx = 0; dx<2; dx++) {
							int indx = base + x + dx;
							pixels[indx] *= rgains[(dy << 1) + 1 - dx];
						}
					}
				}
			}
		}
		//TODO: Add flips here!
		return pixels;
	}

}
