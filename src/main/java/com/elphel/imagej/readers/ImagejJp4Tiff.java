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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.compress.utils.IOUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.tileprocessor.TileNeibs;

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
	private static final boolean BYPASS_SERVICES = false; //  true;
	private static final String TELEMETRY_PREFIX = "TLM_";
	private static final String SERVICES_PATH = "services.properties.forelphel";
	private static final boolean KEEP_EXTENSION = false; // remove extension for ImagePlus title
	private static final boolean BOTTOM_TELEMETRY = true; // telemetry (if present) is at the bottom of the image

	// Fixing bit12 in channel 6 of LWIR 16 camera 00:0E:64:10:C4:35
	private static final String FIXCH6_SERIAL = "00:0E:64:10:C4:35";
	private static final int    FIXCH6_CHANNEL =     2;
	private static final int    FIXCH6_BIT     =    12;
//	private static final int    FIXCH6_MAXVAL  = 23367; // higher - subtract 4096, <19271 -add 4096 
	private static final int    FIXCH6_EXPECTED  = 21319; // expected value 
	private static final String FIXCH6_EARLIEST = "2021-12-01 00:00:00.000";
	private static final String FIXCH6_LATEST =   "2022-04-01 00:00:00.000";
	
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
		if (path_url == null) {
			LOGGER.warn("readTiffJp4(): null path_url");
			return null;
		}
		url = null;
		//	 		String mime = null; // use to select jp4/tiff later? Or to check it is correct
		content_fileName = null;
		boolean bad_file = false;
		try {
			File test_file = new File(path_url); // null
			test_file.getCanonicalPath();
		} catch (IOException e) {
			bad_file = true;
		}
		if (bad_file) {
			try {
				url = new URL(path_url);
			} catch (MalformedURLException e) {
				LOGGER.warn("Bad URL3: " + path_url);
			}
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
		if (bpp == 32) { // already ImageJ file - just read it and decode properties
			imp =  new ImagePlus(content_fileName);
			decodeProperiesFromInfo(imp);
			return imp;
		}

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

		int fixed_center =  fix000E6410C435(
				meta_hash, // Hashtable<String, Object> meta_hash,
				pixels); // short [] pixels)
		
		
		boolean degamma = bytes_per_pixel < 2; // both JP4 and 8-bit tiff
		
		// FIXME: Do not scale telemetry!
		/*
		if (deGammaScale(pixels, reader.getSizeX(), meta_hash, degamma, scale ) == null) {
			LOGGER.error("Problem degamma/scaling of "+content_fileName);

		}
		*/
		boolean telemetry = (bytes_per_pixel == 2);
		boolean lepton_telemetry = telemetry && ((reader.getSizeX() == 160));
		boolean boson_telemetry =  telemetry && ((reader.getSizeX() == 640));
		Lepton3Telemetry lepton3Telemetry = null;
		Boson640Telemetry boson640Telemetry =  null;
		if (telemetry) {
			if (lepton_telemetry) {
				lepton3Telemetry = new Lepton3Telemetry(bb,pixels, true); // bottom
				telemetry = lepton3Telemetry.hasTelemetry();
			} else if (boson_telemetry) {
				boson640Telemetry = new Boson640Telemetry(bb,pixels, false); // top
				telemetry = boson640Telemetry.hasTelemetry();
			}
		}
		ImageProcessor ip = null;
		if (telemetry) {
			if (lepton_telemetry) {
				ip=new FloatProcessor(Lepton3Telemetry.getWidth(), Lepton3Telemetry.getHeight());
				pixels = lepton3Telemetry.getPixels();
//				ip.setPixels(lepton3Telemetry.getPixels()); // bo
			} else if (boson_telemetry) {
				ip=new FloatProcessor(Boson640Telemetry.getWidth(), Boson640Telemetry.getHeight());
				pixels=boson640Telemetry.getPixels(); // bo
//				ip.setPixels(boson640Telemetry.getPixels()); // bo
			}			
		} else {
			ip=new FloatProcessor(reader.getSizeX(), reader.getSizeY());
//			ip.setPixels(pixels);
		}
		
		if (deGammaScale(pixels, reader.getSizeX(), meta_hash, degamma, scale ) == null) {
			LOGGER.error("Problem degamma/scaling of "+content_fileName);

		}
		ip.setPixels(pixels);
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
		if (lepton_telemetry) {
			HashMap<String, String> telemetryMap = lepton3Telemetry.parseTelemetry();
			for (String key:telemetryMap.keySet()) {
				imp.setProperty(TELEMETRY_PREFIX+key, telemetryMap.get(key));
			}
		}
		if (boson_telemetry) {
			HashMap<String, String> telemetryMap = boson640Telemetry.parseTelemetry();
			for (String key:telemetryMap.keySet()) {
				imp.setProperty(TELEMETRY_PREFIX+key, telemetryMap.get(key));
			}
		}
		if (fixed_center > 0) {
			imp.setProperty("FIXED_CHN6", ""+fixed_center);
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
	
	public static boolean decodeProperiesFromInfo(ImagePlus imp){
		if (imp.getProperty("Info")==null) return false;
		String xml= (String) imp.getProperty("Info");

	    DocumentBuilder db=null;
		try {
			db = DocumentBuilderFactory.newInstance().newDocumentBuilder();
		} catch (ParserConfigurationException e) {
			return false;
		}
	    InputSource is = new InputSource();
	    is.setCharacterStream(new StringReader(xml));
    	Document doc = null;
	    try {
	    	doc = db.parse(is);
	    } catch (SAXException e) {
	    	return false;
	    } catch (IOException e) {
	    	return false;
	    }
	    NodeList allNodes=doc.getDocumentElement().getElementsByTagName("*");
	    for (int i=0;i<allNodes.getLength();i++) {
	        String name= allNodes.item(i).getNodeName();
	        String value="";
	        try {
	        	value=allNodes.item(i).getFirstChild().getNodeValue();
	        } catch(Exception e) {

	        }
    		imp.setProperty(name, value);

	    }

		return true;
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
							if (indx >= pixels.length) {
								System.out.println(indx);
								continue;
							}
							pixels[indx] *= rgains[(dy << 1) + 1 - dx];
						}
					}
				}
			}
		}
		//TODO: Add flips here!
		return pixels;
	}
	
	// fixing "00:0E:64:10:C4:35" camera
	public int fix000E6410C435(
			Hashtable<String, Object> meta_hash,
			float [] pixels)
	{
		boolean debug = false;
		int width =  reader.getSizeX();
		int height = reader.getSizeY();
		String serial = (String) meta_hash.get("Serial_Number");
		if ((serial == null) || !serial.equals(FIXCH6_SERIAL)) {
			return -1; // wrong camera 
		}
		String schannel = ((String) meta_hash.get("PageNumber")).substring(0,1);
		if ((schannel == null) || (Integer.parseInt(schannel) != FIXCH6_CHANNEL)) {
			return -2; // wrong channel
		}		
		
		String sfdate = (String) meta_hash.get("DateTime");
		sfdate = sfdate.replaceFirst(":", "-");
		sfdate = sfdate.replaceFirst(":", "-")+".000";
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss.SSS");
		try {
			Date startDate = dateFormat.parse(FIXCH6_EARLIEST);
			Date endDate =   dateFormat.parse(FIXCH6_LATEST);
			Date fileDate =  dateFormat.parse(sfdate);
			if (fileDate.before(startDate) || fileDate.after(endDate)) {
				return -3; // too early or too late 
			}
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int [] hist = new int [4096];
		int max12 = ((int) pixels[640]) & 0xf000;
		int min12 = max12;
		for (int i = 640; i < pixels.length; i++) {
			int sp = (int) pixels[i];
			int sp12 = sp & 0xe000;
			if (sp12 > max12) {
				max12 = sp12;
			} else if (sp12 < min12) {
				min12 = sp12;
			}
			hist[sp & 0xfff]++; 
		}
		int z0 = 0;
		int zbest = 0;
		int lbest = 0;
		while (z0 <= 4095) {
			for (; z0 < (4096-lbest); z0++) {
				if (hist[z0] == 0) {
					break;
				}
			}
			if ( z0 > (4095 - lbest) ) {
				break;
			}
			int l = 1;
			for (; l < 4096; l++) {
				if (hist[(z0 + l) & 0xfff] != 0) {
					break;
				}
			}
			if (l > lbest) {
				zbest = z0;
				lbest = l;
			}
			z0 += l + 1;
		}
		if (lbest == 0) {
			System.out.println ("Failed to fix bit 12 for channel 2 of camera 00:0E:64:10:C4:35");
		}
		int center = 0;
		int split = (zbest+ lbest/2) & 0xfff;
		int vmin = 0;
		if (max12 > min12) {
			vmin = max12 - 4096 + split;
		} else if (lbest < 2048) {
			vmin = min12 + split;
		} else {
			center = min12 + ((split + 2048) & 0xfff);
			boolean use_high = Math.abs(center - FIXCH6_EXPECTED) > Math.abs(center + 4096 - FIXCH6_EXPECTED);
			if (use_high) {
				center += 4096;
			}
			vmin = center - 2048;
		}
		float mn = vmin;
		float mx = mn + 4095;
		float fr = 4096;
		for (int i = 640; i < pixels.length; i++) {
			if (pixels[i] > mx) {
				pixels[i] -= fr; 
			} else if (pixels[i] < mn) {
				pixels[i] += fr;
			}
		}
		center = vmin+2048;
		height = 512; // was 513
		int [] idata = new int [width*height];
		double max_adiff = 4096/2;
		int num_diff = 0;
		double [][] dbg_img = debug? new double [4][width*height] : null;
		for (int iy = 0; iy < height; iy++) {
			for (int ix = 0; ix < width; ix++) {
				int indx = iy*width+ix;
				int sindx = indx + width;
				idata[indx] = (int) pixels[sindx];
				double hdiff=0, vdiff=0, ahdiff = 0, avdiff = 0;
				if (ix < (width - 1)) {
					hdiff = pixels[sindx] - pixels[sindx + 1];
					ahdiff = Math.abs(hdiff);
				}
				if (iy < (height - 1)) {
					vdiff = pixels[sindx] - pixels[sindx + width];
					avdiff = Math.abs(vdiff);
				}
				if (dbg_img != null) {
					dbg_img[0][indx] = pixels[sindx];
					if (ix < (width - 1)) {
						dbg_img[1][indx] = hdiff;
					}					
					if (iy < (height - 1)) {
						dbg_img[2][indx] = vdiff;
					}
					dbg_img[3][indx] = idata[indx] & 0xf000;
				}
				if (Math.max(avdiff, ahdiff) > max_adiff) {
					num_diff++;
				}
			}
		}
		
		if (dbg_img != null) {
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_img,
					width,
					height,
					true,
					"fix000E6410C435");
		}

		if (num_diff == 0) {
			System.out.println ("No discontinuities remain");
			return center;
		}
		System.out.println ("Discontinuities remain, will fix by clusters");
		boolean OK = fixByClusters(
				width,
				height,
				FIXCH6_BIT,         // int    bad_bit,
				idata,              // int [] data)
				(dbg_img != null)); // boolean debug)
		System.out.println ("FixByClusters -> "+OK);
		for (int i = 0; i < idata.length; i++) {
			pixels[i+640] = idata[i];
		}
		return center;
	}

	public static boolean fixByClusters(
			int    width,
			int    height,
			int    bad_bit,
			int [] data,
			boolean debug) {
		TileNeibs tn = new TileNeibs(width,height);
		int [] clusters = new int [width*height];
		int max_diff = (1 << (bad_bit - 1))-1;
		//		  max_diff = (1 << (bad_bit - 2))-1; // verify
		int start_pix =     0;
		int cluster_index = 0;
		ArrayList<Integer> wave_list = new ArrayList<Integer>(1000);
		while (start_pix < clusters.length) {
			// find first empty pixel
			for (; (start_pix < clusters.length) && (clusters[start_pix] != 0); start_pix++);
			if (start_pix >= clusters.length) {
				break;
			}
			cluster_index++;
			wave_list.add(start_pix);
			clusters[start_pix] = cluster_index;
			while (!wave_list.isEmpty()) {
				int indx = wave_list.remove(0);
				int center_value = data[indx];
				int mn = center_value - max_diff;
				int mx = center_value + max_diff;
				for (int dir = 0; dir < 8; dir++) {
					int indx_new = tn.getNeibIndex(indx, dir);
					if ((indx_new >= 0) && (clusters[indx_new] == 0)){
						int val = data[indx_new];
						if ((val > mn) && (val < mx)) {
							wave_list.add(indx_new);
							clusters[indx_new] = cluster_index;
						}
					}
				}
			}
		}
		if (cluster_index < 2) {
			return false;//  nothing to do, probably error - reduce max_diff?
		}
		int [] clust_sizes = new int [cluster_index];
		for (int i = 0; i < clusters.length; i++) {
			clust_sizes[clusters[i]-1]++;
		}
		// Find the largest cluster -assuming it is correct
		int largest = 0;
		int [] cluster_offset = new int [cluster_index];
		boolean [] cluster_offset_set = new boolean[cluster_index];
		for (int i = 1; i < cluster_index; i++) {
			if (clust_sizes[i] > clust_sizes[largest]) {
				largest = i;
			}
		}
		cluster_offset_set[largest] = true; // largest cluster has offset 0
		int num_unset = cluster_index - 1;
		int indx = 0;
		int indx_new = 0;
		int [] dirs = {TileNeibs.DIR_E,TileNeibs.DIR_S};
		int ioffs = 1 << bad_bit; // 4096
		for (;num_unset > 0; num_unset--) {
			int ntry = 0;
			for (; ntry < (2*clusters.length); ntry++) {
				int ndir = ntry & 1;
				int clust1 = clusters[indx]-1;
				indx_new = tn.getNeibIndex(indx, dirs[ndir]);
				if ((indx_new >= 0) && (cluster_offset_set[clust1] != cluster_offset_set[clusters[indx_new]-1])) {
					break;
				}
				if (ndir == 1) {
					indx++;
					if (indx >= clusters.length) { // start from the beginning
						indx -= clusters.length;
					}
				}
			}
			if (ntry >= (2*clusters.length)) {
				System.out.println("fixByClusters() BUG!");
				break;
			}
			int indx_set =  indx;// cluster_offset_set[clusters[indx]] ? indx
			int indx_nset = indx_new;
			if (cluster_offset_set[clusters[indx_new]-1]) {
				indx_set =  indx_new;// cluster_offset_set[clusters[indx]] ? indx
				indx_nset = indx;
			}
			
			if (data[indx_nset] > data[indx_set]) {
				cluster_offset[clusters[indx_nset]-1] = cluster_offset[clusters[indx_set]-1] - ioffs; 
			} else {
				cluster_offset[clusters[indx_nset]-1] = cluster_offset[clusters[indx_set]-1] + ioffs; 
			}
			cluster_offset_set[clusters[indx_nset]-1] = true;
		}
		double [][] dbg_img = null;
		if (debug) {
			dbg_img = new double [4][clusters.length];
			for (int i = 0; i< clusters.length; i++) {
				dbg_img[0][i] = clusters[i];
				dbg_img[1][i] = cluster_offset[clusters[i]-1];
				dbg_img[2][i] = data[i];
			}
		}
		// apply correction
		for (int i = 0; i < clusters.length; i++) {
			data[i] += cluster_offset[clusters[i]-1];
		}
		if (debug) {
			for (int i = 0; i< clusters.length; i++) {
				dbg_img[3][i] = data[i];
			}
		}
		
		if (dbg_img != null) {
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					dbg_img,
					width,
					height,
					true,
					"clusters",
					new String [] {"cluster", "offset","data_in","data_out"});
		}
		

		
		return true;
	}
	
}
