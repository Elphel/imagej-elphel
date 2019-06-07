package com.elphel.imagej.readers;
/**
** -----------------------------------------------------------------------------**
** EyesisTiff.java
**
** Writes Tiff files suitable for Emblend, preserve ImageJ Info data
** Uses bioformat library
**
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  EyesisTiff.java is free software: you can redistribute it and/or modify
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
import java.awt.Frame;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBufferInt;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.lang.reflect.Array;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

//import org.apache.log4j.Logger;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.io.FileInfo;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import loci.common.RandomAccessInputStream;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.ClassList;
import loci.formats.FormatException;
import loci.formats.IFormatReader;
import loci.formats.ImageReader;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.formats.tiff.IFD;
import loci.formats.tiff.IFDList;
import loci.formats.tiff.TiffParser;
import loci.formats.tiff.TiffRational;
import loci.formats.tiff.TiffSaver;

public class EyesisTiff {
	  private static ClassList<IFormatReader> defaultClasses;
	  private static final Logger LOGGER = LoggerFactory.getLogger(ClassList.class);
	  public static ClassList<IFormatReader> getCustomReaderClasses() {
			defaultClasses = null;
		    if (defaultClasses == null) {
		      try {
		        defaultClasses =
		          new ClassList<IFormatReader>(
		        		  EyesisTiff.class.getClassLoader().getResource("readers.txt").getFile(), // @param file Configuration file containing the list of classes.
//		        		  EyesisTiff.class.getClassLoader().getResource("readers1.txt").getFile(), // @param file Configuration file containing the list of classes.
		        		  IFormatReader.class, // @param base Base class to which all classes are assignable.
		        		  null); // @param location Class indicating which package to search for the file.
		      }
		      catch (IOException exc) {
		        defaultClasses = new ClassList<IFormatReader>(IFormatReader.class);
		        LOGGER.info("Could not parse class list; using default classes", exc);
		      }
		    }
//		    ClassList<IFormatReader> dc = defaultClasses;
	        LOGGER.info("Loaded "+defaultClasses.getClasses().length+" classes");
		    return defaultClasses;
		  }



	//	private static org.apache.log4j.Logger log= Logger.getLogger(EyesisTiff.class);

	public EyesisTiff(){
		//	Please initialize the log4j system properly


	}
	public ImagePlus readTiff(String path) {
		return readTiff(path, 	"STD_"); // null);
	}

	public ImagePlus readTiff(String path, String std ) { // std - include non-elphel properties with prefix std

		String inId = null;
		inId = path;

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

	    IMetadata omeMeta = null;
		try {
			omeMeta = service.createOMEXMLMetadata();
		} catch (ServiceException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

//		final Class<? extends IFormatReader>[] defaultClasses =
//				ImageReader.getDefaultReaderClasses().getClasses();
//			final int currentHash = Arrays.hashCode(defaultClasses);
//https://www.javatips.net/api/libbio-formats-java-master/components/scifio/src/loci/formats/in/BaseTiffReader.java
//https://docs.openmicroscopy.org/bio-formats/5.7.2/developers/reader-guide.html
		//https://docs.openmicroscopy.org/bio-formats/5.9.2/developers/java-library.html#file-reading-and-performance

		//defaultClasses = new ClassList<IFormatReader>(IFormatReader.class);
//		ClassList<IFormatReader> cl0 = getCustomReaderClasses(); // new ClassList<IFormatReader>(IFormatReader.class);
		ClassList<IFormatReader> classList =  new ClassList<IFormatReader>(IFormatReader.class);
		classList.addClass(com.elphel.imagej.readers.ElphelTiffReader.class);
	    ImageReader reader = new ImageReader(getCustomReaderClasses());
//	    ImageReader reader = new ImageReader(classList);
	    reader.setMetadataStore(omeMeta);
	    try {
			reader.setId(inId);
		} catch (FormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// a lot of output for each unmatched format
			//e.printStackTrace();
		}
	    /* read-end */



	    int bpp =   reader.getBitsPerPixel();
	    byte [] bytes = null;
	    ImagePlus imp=  null;
	    try {
			bytes = reader.openBytes(0);
		} catch (FormatException e1) {
	        LOGGER.warn("Invalid image format, error "+e1);
		} catch (IOException e1) {
	        LOGGER.warn("I/O error "+e1);
		}
	    if (bytes != null) {
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
		    String imageName = path;
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
					imp.setProperty(std+(key.replace(" ","_")), meta_hash.get(key).toString());
				}
			}
			encodeProperiesToInfo(imp);
	    }

	    try {
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return imp;
	}

/*
 * image width 	getSizeX()
image height 	getSizeY()
number of series per file 	getSeriesCount()
total number of images per series 	getImageCount()
number of slices in the current series 	getSizeZ()
number of timepoints in the current series 	getSizeT()
number of actual channels in the current series 	getSizeC()
number of channels per image 	getRGBChannelCount()
the ordering of the images within the current series 	getDimensionOrder()
whether each image is RGB 	isRGB()
whether the pixel bytes are in little-endian order 	isLittleEndian()
whether the channels in an image are interleaved 	isInterleaved()
the type of pixel data in this file 	getPixelType()


 */

	public void savePNG_ARGB32(
			ImagePlus imp,
			String path
			){
		int width = imp.getWidth();
		int height = imp.getHeight();
		int [] pixels = (int []) imp.getProcessor().getPixels();
		System.out.println("savePNG_ARGB32("+path+"): width="+width+", height="+height+" length="+pixels.length);

		DataBufferInt buffer = new DataBufferInt(pixels, pixels.length);

		int[] bandMasks = {0xFF0000, 0xFF00, 0xFF, 0xFF000000}; // ARGB (yes, ARGB, as the masks are R, G, B, A always) order
		WritableRaster raster = Raster.createPackedRaster(buffer, width, height, width, bandMasks, null);

		ColorModel cm = ColorModel.getRGBdefault();
		BufferedImage bimage = new BufferedImage(cm, raster, cm.isAlphaPremultiplied(), null);
		try{
			(new File(path)).delete(); // Otherwise TiffSaver appends! - is it the same for PNG?
			ImageIO.write(bimage, "PNG", (new File(path)));
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}



	public void saveTiff(
			ImagePlus imp,
			String path,
			int mode,     // 0 - 8-bit, 1 - 16-bit, 2 - 32-bit unsigned, 3 - 32FP (only if source image is a 4-stack of FP)
			double scale, // full scale, absolute
			boolean imageJTags,
			int debugLevel
			) throws IOException, FormatException, ServiceException, DependencyException{
		if (imp.getType()==ImagePlus.COLOR_RGB) { // 4*8 bit - type = 0
			if (debugLevel>1) System.out.println("Saving 8-bit TIFF with alpha-channel: "+path); // is it possible to just add alpha to high bytes?
			saveTiffARGB32(imp, path, imageJTags, debugLevel);
			return;
		} else if (imp.getStackSize()==4) {
			if (debugLevel>1) System.out.println("Saving 32-bit float TIFF with alpha-channel: "+path);
			saveTiffARGB(imp, path, mode, scale, imageJTags,debugLevel);
			return;

		}
		IJ.showMessage("Not yet implemented for this image type");
	}

	public void saveTiffARGB(
			ImagePlus imp,
			String path,
			int mode,     // 0 - 8-bit, 1 - 16-bit, 2 - 32-bit unsigned, 3 - 32FP (only if source image is a 4-stack of FP)
			double scale, // full scale, absolute
			boolean imageJTags,
			int debugLevel
			) throws IOException, FormatException, ServiceException, DependencyException{
		int IFDImageFullWidth= 0x8214; // not defined in loci.formats.tiff.IFD
		int IFDImageFullLength=0x8215; // not defined in loci.formats.tiff.IFD

		int IFDImageJByteCounts= 0xc696; // was array {12( if no slices, roi, etc.), bytes in info}
		int IFDImageJInfo=       0xc697; // ImageJ info, starting with magic IJIJinfo,
		byte [] ImageJInfoMagic={73,74,73,74,105,110,102,111,0,0,0,1};

		int pixelsDenominator=1000;
		String description=(imp.getProperty("description")!=null)?((String) imp.getProperty("description")):"Elphel Eyesis4pi";
		//		int [] iPixels= (int []) imp.getProcessor().getPixels();
		final float [][] imagePixels= new float [4][];
		try {
			for (int c=0;c<imagePixels.length;c++){
				imagePixels[c]=	(float []) imp.getStack().getPixels(c+1); // bytes are little endian
			}
		} catch (ClassCastException e) {
			try {
				for (int c=0;c<imagePixels.length;c++) {
					byte [] bpixels = (byte []) imp.getStack().getPixels(c+1); // bytes are little endian
					imagePixels[c] = new float[bpixels.length];
					for (int i = 0; i <bpixels.length; i++)
						imagePixels[c][i]=	((float) bpixels[i])/255;
				}
			} catch (ClassCastException e1) {
				for (int c=0;c<imagePixels.length;c++) {
					short [] spixels = (short []) imp.getStack().getPixels(c+1); // bytes are little endian
					imagePixels[c] = new float[spixels.length];
					for (int i = 0; i < spixels.length; i++)
						imagePixels[c][i]=	((float) spixels[i])/65535;
				}
			}
		}
		//		for (int i = 0; i < imagePixels[3].length; i++){
		//			imagePixels[3][i] = 1.0f - imagePixels[3][i];
		//		}
		int bw;
		byte [] bytes;
		int pIndex=0;
		//        double s;
		int pixelType=0; //int pixelType:INT8, UINT8, INT16, UINT16, INT32, UINT32, FLOAT, BIT, DOUBLE, COMPLEX, DOUBLECOMPLEX; - used to find bpp
		switch (mode){
		case 0:
			int maxVal= 0xff;
			double [] s0 ={scale*maxVal,scale*maxVal,scale*maxVal,maxVal};
			bw= 1;
			pixelType=0; //1; //0; //3; // 0; // should it be 3?
			bytes=new byte[imagePixels[0].length*4*bw];
			for (int i=0;i<imagePixels[0].length;i++) for (int c=0;c<imagePixels.length;c++){ // r,g,b,alpha
				int d= (int) (imagePixels[c][i]*s0[c]);
				if (d<0) d=0;
				else if (d>maxVal) d= maxVal;
				bytes[pIndex++]=(byte) ((d>> 0 )& 0xff);
			}
			break;
		case 1:
			int iMaxVal= 0xffff;
			double [] s1 ={scale*iMaxVal,scale*iMaxVal,scale*iMaxVal,iMaxVal};
			bw= 2;
			pixelType=2; //3; //4;
			bytes=new byte[imagePixels[0].length*4*bw];
			for (int i=0;i<imagePixels[0].length;i++) for (int c=0;c<imagePixels.length;c++){ // r,g,b,alpha
				int d= (int) (imagePixels[c][i]*s1[c]);
				if (d<0) d=0;
				else if (d>iMaxVal) d= iMaxVal;
				bytes[pIndex++]=(byte) ((d>> 8 )& 0xff);
				bytes[pIndex++]=(byte) ((d>> 0 )& 0xff);
			}
			break;
		case 2:
			long lMaxVal= 0xffffffffL;
			double [] s2 ={scale*lMaxVal,scale*lMaxVal,scale*lMaxVal,lMaxVal};
			bw= 3; //4;
			pixelType=5; // 2; //5;
			bytes=new byte[imagePixels[0].length*4*bw];
			for (int i=0;i<imagePixels[0].length;i++) for (int c=0;c<imagePixels.length;c++){ // r,g,b,alpha
				long d= (long) (imagePixels[c][i]*s2[c]);
				if (d<0) d=0;
				else if (d>lMaxVal) d= lMaxVal;
				bytes[pIndex++]=(byte) ((d>>24 )& 0xff);
				bytes[pIndex++]=(byte) ((d>>16 )& 0xff);
				bytes[pIndex++]=(byte) ((d>> 8 )& 0xff);
				bytes[pIndex++]=(byte) ((d>> 0 )& 0xff);
			}
			break;
		case 3:
			bw= 4;
			pixelType=6;
			bytes=new byte[imagePixels[0].length*4*bw];
			// what scale should be for alpha 0..1 or 0.. scale?
			for (int i=0;i<imagePixels[0].length;i++) for (int c=0;c<imagePixels.length;c++){ // r,g,b,alpha
				int d;
				if (scale==1.0) d=Float.floatToIntBits(imagePixels[c][i]);
				else d=Float.floatToIntBits((float) (imagePixels[c][i]*scale));
				bytes[pIndex++]=(byte) ((d>>24 )& 0xff);
				bytes[pIndex++]=(byte) ((d>>16 )& 0xff);
				bytes[pIndex++]=(byte) ((d>> 8 )& 0xff);
				bytes[pIndex++]=(byte) ((d>> 0 )& 0xff);
			}
			break;
		default:
			IJ.error("saveTiffARGBFloat32", "Unsupported output format mode ="+mode);
			return;
		}
		//        System.out.println("saveTiffARGBFloat32(): mode="+mode+" pixelType="+pixelType+" bw="+bw);
		IFD ifd=new IFD();
		ifd.put(new Integer(IFD.LITTLE_ENDIAN), new Boolean(false));
		//        ifd.put(new Integer(IFD.LITTLE_ENDIAN), new Boolean(true));
		ifd.put(new Integer(IFD.IMAGE_WIDTH), imp.getWidth());
		ifd.put(new Integer(IFD.IMAGE_LENGTH), imp.getHeight());
		ifd.put(new Integer(IFD.SAMPLES_PER_PIXEL), 4);
		ifd.putIFDValue(IFD.SOFTWARE, "Elphel Eyesis");
		ifd.putIFDValue(IFD.IMAGE_DESCRIPTION, description);
		// copy some other data?
		ifd.putIFDValue(IFD.COMPRESSION, 1); //TiffCompression.UNCOMPRESSED);
		ifd.putIFDValue(IFD.PHOTOMETRIC_INTERPRETATION,2); // RGB
		ifd.putIFDValue(IFD.EXTRA_SAMPLES,2); // 0 = Unspecified data  1 = Associated alpha data (with pre-multiplied color) 2 = Unassociated alpha data
		//        int [] bpsArray={8,8,8,8};
		//        ifd.putIFDValue(IFD.BITS_PER_SAMPLE, bpsArray); // will be done automatically
		if (imp.getProperty("XPosition")!=null) {
			ifd.putIFDValue(IFD.X_POSITION,
					new TiffRational((int) Math.round(pixelsDenominator*Double.parseDouble((String) imp.getProperty("XPosition"))) , pixelsDenominator));
		}
		if (imp.getProperty("YPosition")!=null) {
			ifd.putIFDValue(IFD.Y_POSITION,
					new TiffRational((int) Math.round(pixelsDenominator*Double.parseDouble((String) imp.getProperty("YPosition"))) , pixelsDenominator));
		}
		if (imp.getProperty("ImageFullWidth")!=null){
			ifd.putIFDValue(IFDImageFullWidth, 	(long) Integer.parseInt((String) imp.getProperty("ImageFullWidth")));
		}
		if (imp.getProperty("ImageFullLength")!=null){
			ifd.putIFDValue(IFDImageFullLength, (long) Integer.parseInt((String) imp.getProperty("ImageFullLength")));
		}
		//TODO: Seems to match ImageJ Info, but it is not recognized :-(
		if (imageJTags && (imp.getProperty("Info")!=null) && (imp.getProperty("Info") instanceof String)){
			int skipFirstBytes=2;
			String info=(String) imp.getProperty("Info");
			byte [] bInfoBody=info.getBytes("UTF-16");
			int [] bInfo = new int [ImageJInfoMagic.length+bInfoBody.length-skipFirstBytes];
			int index=0;
			for (int i=0;i<ImageJInfoMagic.length;i++) bInfo[index++]=ImageJInfoMagic[i];
			for (int i=skipFirstBytes;i<bInfoBody.length;      i++) bInfo[index++]=bInfoBody[i]; // first 2 bytes {-2, -1} ???
			long [] imageJcounts={12, bInfoBody.length-skipFirstBytes};
			ifd.putIFDValue(IFDImageJByteCounts, imageJcounts);
			ifd.putIFDValue(IFDImageJInfo, bInfo);

		}
		(new File(path)).delete(); // Otherwise TiffSaver appends!
		TiffSaver tiffSaver = new TiffSaver(path);
		tiffSaver.setWritingSequentially(true);
		tiffSaver.setLittleEndian(false);
		tiffSaver.writeHeader();
		//        tiffSaver.writeIFD(ifd,0); //* SHould not write here, some fields are calculated during writeImage, that writes IFD too
		//        System.out.println("bytes.length="+bytes.length);
		tiffSaver.writeImage(bytes,
				ifd,
				0, //int no,
				pixelType, // 0, //int pixelType:INT8, INT16, INT32, UINT8, UINT16, UINT32, FLOAT, BIT, DOUBLE, COMPLEX, DOUBLECOMPLEX; - used to find bpp
				true); // boolean last)
	}

	public void saveTiffARGB32(
			ImagePlus imp,
			String path,
			boolean imageJTags,
			int debugLevel
			) throws IOException, FormatException, ServiceException, DependencyException{
		int IFDImageFullWidth= 0x8214; // not defined in loci.formats.tiff.IFD
		int IFDImageFullLength=0x8215; // not defined in loci.formats.tiff.IFD
		//		public static final int META_DATA_BYTE_COUNTS = 50838; // private tag registered with Adobe
		//		public static final int META_DATA = 50839; // private tag registered with Adobe

		int IFDImageJByteCounts= 0xc696; // was array {12( if no slices, roi, etc.), bytes in info}
		int IFDImageJInfo=       0xc697; // ImageJ info, starting with magic IJIJinfo,
		byte [] ImageJInfoMagic={73,74,73,74,105,110,102,111,0,0,0,1};

		int pixelsDenominator=1000;
		String description=(imp.getProperty("description")!=null)?((String) imp.getProperty("description")):"Elphel Eyesis4pi";
		int [] iPixels= (int []) imp.getProcessor().getPixels();
		byte [] bytes=new byte[iPixels.length*4];
		int pIndex=0;
		for (int i=0;i<iPixels.length;i++){
			bytes[pIndex++]=(byte) ((iPixels[i]>>16)& 0xff); // R
			bytes[pIndex++]=(byte) ((iPixels[i]>> 8)& 0xff); // G
			bytes[pIndex++]=(byte) ((iPixels[i]>> 0)& 0xff); // B
			bytes[pIndex++]=(byte) ((iPixels[i]>>24)& 0xff); // alpha
		}
		IFD ifd=new IFD();
		ifd.put(new Integer(IFD.LITTLE_ENDIAN), new Boolean(false));
		//        ifd.put(new Integer(IFD.LITTLE_ENDIAN), new Boolean(true));
		ifd.put(new Integer(IFD.IMAGE_WIDTH), imp.getWidth());
		ifd.put(new Integer(IFD.IMAGE_LENGTH), imp.getHeight());
		ifd.put(new Integer(IFD.SAMPLES_PER_PIXEL), 4);
		ifd.putIFDValue(IFD.SOFTWARE, "Elphel Eyesis");
		ifd.putIFDValue(IFD.IMAGE_DESCRIPTION, description);
		// copy some other data?
		ifd.putIFDValue(IFD.COMPRESSION, 1); //TiffCompression.UNCOMPRESSED);
		ifd.putIFDValue(IFD.PHOTOMETRIC_INTERPRETATION,2); // RGB
		ifd.putIFDValue(IFD.EXTRA_SAMPLES,2); // extra bytes (over 3) meaning Unassociated alpha data

		//        int [] bpsArray={8,8,8,8};
		//        ifd.putIFDValue(IFD.BITS_PER_SAMPLE, bpsArray); // will be done automatically
		if (imp.getProperty("XPosition")!=null) {
			ifd.putIFDValue(IFD.X_POSITION,
					new TiffRational((int) Math.round(pixelsDenominator*Double.parseDouble((String) imp.getProperty("XPosition"))) , pixelsDenominator));
		}
		if (imp.getProperty("YPosition")!=null) {
			ifd.putIFDValue(IFD.Y_POSITION,
					new TiffRational((int) Math.round(pixelsDenominator*Double.parseDouble((String) imp.getProperty("YPosition"))) , pixelsDenominator));
		}
		if (imp.getProperty("ImageFullWidth")!=null){
			ifd.putIFDValue(IFDImageFullWidth, 	(long) Integer.parseInt((String) imp.getProperty("ImageFullWidth")));
		}
		if (imp.getProperty("ImageFullLength")!=null){
			ifd.putIFDValue(IFDImageFullLength, (long) Integer.parseInt((String) imp.getProperty("ImageFullLength")));
		}
		//TODO: Seems to match ImageJ Info, but it is not recognized :-(
		if (imageJTags &&  (imp.getProperty("Info")!=null) && (imp.getProperty("Info") instanceof String)){
			int skipFirstBytes=2;
			String info=(String) imp.getProperty("Info");
			byte [] bInfoBody=info.getBytes("UTF-16");
			int [] bInfo = new int [ImageJInfoMagic.length+bInfoBody.length-skipFirstBytes];
			int index=0;
			for (int i=0;i<ImageJInfoMagic.length;i++) bInfo[index++]=ImageJInfoMagic[i];
			for (int i=skipFirstBytes;i<bInfoBody.length;      i++) bInfo[index++]=bInfoBody[i]; // first 2 bytes {-2, -1} ???
			/*
        	StringBuffer sb=new StringBuffer("bInfo: ");
        	for (int i=0;i<bInfo.length;i++) sb.append(bInfo[i]+" ");
        	System.out.println(sb.toString());
        	sb=new StringBuffer("ImageJInfoMagic: ");
        	for (int i=0;i<ImageJInfoMagic.length;i++) sb.append(ImageJInfoMagic[i]+" ");
        	System.out.println(sb.toString());
        	sb=new StringBuffer("bInfoBody: ");
        	for (int i=0;i<bInfoBody.length;i++) sb.append(bInfoBody[i]+" ");
        	System.out.println(sb.toString());
        	System.out.println("info[0]="+info.charAt(0));
        	System.out.println("info[1]="+info.charAt(1));
        	System.out.println("info[2]="+info.charAt(2));
			 */
			long [] imageJcounts={12, bInfoBody.length-skipFirstBytes};
			ifd.putIFDValue(IFDImageJByteCounts, imageJcounts);
			ifd.putIFDValue(IFDImageJInfo, bInfo);

		}
		(new File(path)).delete(); // Otherwise TiffSaver appends!
		TiffSaver tiffSaver = new TiffSaver(path);
		tiffSaver.setWritingSequentially(true);
		tiffSaver.setLittleEndian(false);
		tiffSaver.writeHeader();
		//        tiffSaver.writeIFD(ifd,0); //* SHould not write here, some fields are calculated during writeImage, that writes IFD too
		System.out.println("bytes.length="+bytes.length);
		tiffSaver.writeImage(bytes,
				ifd,
				0, //int no,
				0, //int pixelType:INT8, INT16, INT32, UINT8, UINT16, UINT32, FLOAT, BIT, DOUBLE, COMPLEX, DOUBLECOMPLEX; - used to find bpp
				true); // boolean last)
	}




	public void propertiesTiff(ImagePlus imp){
		FileInfo fi = imp.getOriginalFileInfo();
		if ((fi==null) ||(fi.directory==null) ||  (fi.fileFormat!=FileInfo.TIFF)) {
			IJ.error("TIFF Dumper", "File path not available or not TIFF file");
			return;
		}
		String path = fi.directory + fi.fileName;
		IJ.log("\\Clear");
		IJ.log("PATH = "+path);
		try {
			dumpIFDs(path);
		} catch(IOException e) {
			IJ.error("Tiff Dumper", ""+e);
		}
		Frame log = WindowManager.getFrame("Log");
		if (log!=null) log.toFront();

	}

	public static void dumpIFDs(String path) throws IOException {
		IJ.showStatus("Parsing IFDs");
		RandomAccessInputStream in = new RandomAccessInputStream(path);

		//TiffParser parser = new TiffParser(in);
		TiffParser parser = new TiffParser(in);
		IFDList ifdList = parser.getIFDs();
		IJ.showStatus("");
		for (IFD ifd : ifdList) {
			for (Integer key : ifd.keySet()) {
				int k = key.intValue();
				String name = IFD.getIFDTagName(k)+String.format("(%d [0x%x])", k,k);
				String value = prettyValue(ifd.getIFDValue(k), 0);
				IJ.log(name + " = " + value);
			}
		}
		in.close();
	}

	private static String prettyValue(Object value, int indent) {
		if (!value.getClass().isArray()) return value.toString()+" ("+value.getClass().toString()+")";
		char[] spaceChars = new char[indent];
		Arrays.fill(spaceChars, ' ');
		String spaces = new String(spaceChars);
		StringBuilder sb = new StringBuilder();
		sb.append("{\n");
		for (int i=0; i<Array.getLength(value); i++) {
			sb.append(spaces);
			sb.append(" ");
			Object component = Array.get(value, i);
			sb.append(prettyValue(component, indent + 2));
			sb.append("\n");
		}
		sb.append(spaces);
		sb.append("}");
		byte [] bstring=new byte [Array.getLength(value)];

		for (int i=0;i<bstring.length;i++) bstring[i]= (byte) Integer.parseInt(Array.get(value, i).toString());
		//   String astring=new String((byte []) value);
		String astring="";
		try {
			astring = new String(bstring,"UTF-16");
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		sb.append("\n\""+astring+"\"");
		return sb.toString();
	}

// copied from JP46_Reader_camera.java
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

	public static  boolean decodeProperiesFromInfo(ImagePlus imp){
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
//            System.out.print(name+" -> ");
            String value = "";
            try {
            	value=allNodes.item(i).getFirstChild().getNodeValue();
            } catch (Exception e) {
            }
//            System.out.println(value);
            imp.setProperty(name, value);
	    }
		return true;
	}


}
