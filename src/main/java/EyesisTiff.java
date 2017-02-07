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
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.lang.reflect.Array;
import java.util.Arrays;

import org.slf4j.LoggerFactory;
import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;


//import org.apache.log4j.Logger;


import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.io.FileInfo;
import loci.common.RandomAccessInputStream;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.codec.CodecOptions;
import loci.formats.tiff.IFD;
import loci.formats.tiff.IFDList;
import loci.formats.tiff.TiffParser;
import loci.formats.tiff.TiffRational;
import loci.formats.tiff.TiffSaver;
import loci.formats.tiff.TiffCompression;


public class EyesisTiff {
//	private static org.apache.log4j.Logger log= Logger.getLogger(EyesisTiff.class); 
	
private String codec="UNCOMPRESSED";

public EyesisTiff(){
//	Please initialize the log4j system properly


}

public EyesisTiff(String codec){
//	Please initialize the log4j system properly
	this.codec=codec;
}

	public void saveTiff(
			ImagePlus imp,
			String path,
			int mode,     // 0 - 8-bit, 1 - 16-bit, 2 - 32-bit unsigned, 3 - 32FP (only if source image is a 4-stack of FP)
			double scale, // full scale, absolute
			boolean imageJTags,
			int debugLevel
	) throws IOException, FormatException, ServiceException, DependencyException{
		Level [] loggerLevel={Level.OFF,Level.ERROR,Level.WARN,Level.INFO,Level.DEBUG,Level.TRACE};
		Logger TiffSaverLogger = (Logger)LoggerFactory.getLogger(TiffSaver.class);
		TiffSaverLogger.setLevel(debugLevel>5?Level.ALL:loggerLevel[debugLevel]);

		if (imp.getType()==ImagePlus.COLOR_RGB) {
			if (debugLevel>1) System.out.println("Saving 8-bit TIFF with alpha-channel: "+path);
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
		for (int c=0;c<imagePixels.length;c++) imagePixels[c]=	(float []) imp.getStack().getPixels(c+1); // bytes are little endian
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
        ifd.putIFDValue(IFD.COMPRESSION, TiffCompression.valueOf(codec).getCode());
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
        tiffSaver.setCodecOptions(TiffCompression.valueOf(codec).getCompressionCodecOptions(ifd));
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
        ifd.putIFDValue(IFD.COMPRESSION, TiffCompression.valueOf(codec).getCode());
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

}
