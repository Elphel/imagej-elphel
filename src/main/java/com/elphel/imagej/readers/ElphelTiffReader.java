/**
 ** -----------------------------------------------------------------------------**
 ** ElphelTiffReader.java
 **
 ** loci.format compatible reader for Elphel 8/16 bpp monochrome Tiff files with
 ** MakerNote
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ElphelTiffReader.java is free software: you can redistribute it and/or modify
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
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Date;
import java.util.HashMap;
import java.util.Hashtable;

import org.apache.commons.compress.utils.IOUtils;
import org.joda.time.DateTime;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.drew.imaging.ImageMetadataReader;
import com.drew.metadata.Metadata;
import com.drew.metadata.Tag;
import com.drew.metadata.exif.ExifIFD0Directory;
import com.drew.metadata.exif.ExifSubIFDDirectory;
import com.drew.metadata.exif.GpsDirectory;

import loci.common.ByteArrayHandle;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.FormatTools;
import loci.formats.in.MetadataLevel;
import loci.formats.in.TiffReader;
import loci.formats.meta.MetadataStore;
import loci.formats.tiff.IFD;
import loci.formats.tiff.IFDList;
import loci.formats.tiff.TiffRational;
import ome.xml.model.primitives.Timestamp;

/*
  // non-IFD tags (for internal use)
  public static final int LITTLE_ENDIAN = 0;
  public static final int BIG_TIFF = 1;
  public static final int REUSE = 3;
   <MakerNote          tag="0x927c" format="LONG" count="16"  seq="26" dlen="64"/> 37500
package loci.formats.tiff;
IFD.java  public static final int MAKER_NOTE = 37500;

 * IFD data
{       0=false, LITTLE_ENDIAN
        1=false, BIG_TIFF
0x0100 256=160,                                     ImageWidth
0x0101 257=120,                                     ImageLength
0x0102 258= 16,                                     BitsPerSample
0x0106 262=  1,                                     PhotometricInterpretation
0x0129 297=  [I@3d136141,                           PageNumber
0x8769 34665=242,                                   ExifTag (A pointer to the Exif IFD.)
0x010e 270=[Ljava.lang.String;@640a5e9d,            ImageDescription
0x010f 271=Elphel,                                  Make
0xc62f 50735=00:0E:64:10:8F:15,                     CameraSerialNumber
0x0110 272=LEPTON35_15,                             Model
0x0111 273=557,                                     StripOffsets
0x0131 305=https://git.elphel.com/Elphel/elphel393, Software
0x9211 37393=122,                                   ImageNumber
0x0112 274=1,                                       Orientation
0x0132 306=2019:05:06 14:41:51,                     DateTime
0x0115 277=1,                                       SamplesPerPixel
0x0116 278=120,                                     RowsPerStrip
0x0117 279=38400,                                   StripByteCounts
0x00fe 254=0                                        NewSubfileType }
  */
public class ElphelTiffReader extends TiffReader{ // BaseTiffReader {

	// -- Constants --

	public static final String ELPHEL_PROPERTY_PREFIX = "ELPHEL_";
	public static final String CONTENT_FILENAME = "CONTENT_FILENAME";
    public static final int    TAG_IMAGE_NUMBER          = 0x9211; //com.drew.metadata.exif.ExifDirectoryBase
    public static final int    TAG_SERIAL_NUMBER         = 0xc62f;
    public static final int    TAG_IMAGE_DESCRIPTION     = 0x010e;

// added 08/03/2021    
	public static final String MAKER_NOTE =            "Makernote";
	public static final String SUB_SEC_TIME_ORIGINAL = "SUBSEC_TIME_ORIGINAL"; // "Sub-Sec Time Original";
	public static final String EXPOSURE_TIME =         "Exposure Time";
	public static final String DATE_TIME_ORIGINAL =    "DATE_TIME_ORIGINAL"; // "Date/Time Original";
	public static final boolean REORDER = true; // false;
	public static final String[][] REPLACEMENT_TAGS = // to/from! TODO: check/add for GPS?
		   {{"SUBSEC_TIME_ORIGINAL",      "Sub-Sec Time Original"},
			{"DATE_TIME_ORIGINAL",        "Date/Time Original"},
			{"Instrument_Make",           "Make"},
			{"Serial_Number",             "Unknown tag (0xc62f)"},
			{"Instrument_Model",          "Model"}};
    
    
    
    

	  /** Merge SubIFDs into the main IFD list. */
//	  protected transient boolean mergeSubIFDs = true; // false;

	/** Logger for this class. */
	private static final Logger LOGGER =
			LoggerFactory.getLogger(ElphelTiffReader.class);

	//	  public static final String[] ELPHEL_TIFF_SUFFIXES =
	//	    {"tif", "tiff"}; // , "tf2", "tf8", "btf"};

	//	  public static final String[] COMPANION_SUFFIXES = {"xml", "txt"};

	//	  public static final int IMAGEJ_TAG = 50839;

	// -- Fields --
//	private String inId = null; // to close Location.mapFile
	private URL url = null; // save here actual URL when reading file to memory
	private String content_fileName = null; // from Content-disposition
	private boolean mapped_externally = false; // file is read/mapped externally, do not close it here
	private boolean file_initialized = false;
	
	//	  private String companionFile;
	//	  private String description;
	//	  private String calibrationUnit;
	//	  private Double physicalSizeZ;
	//	  private Double timeIncrement;
	//	  private Integer xOrigin, yOrigin;
	private ExifSubIFDDirectory directory;
	private ExifIFD0Directory   directory_ifd0;
	private GpsDirectory        directory_gps;
	private HashMap<String,String> REPLACEMENT_TAG_MAP = null; // per instance
	

	// -- Constructor --

	/** Constructs a new Tiff reader. */
	public ElphelTiffReader() {
		super(); // "Tagged Image File Format", ELPHEL_TIFF_SUFFIXES); // See if we can use TiffReader without its parent
		//		mergeSubIFDs = true; // false;
		// TODO: See if the selection is just between 2 readers (jp4 and tiff - just Elphel cameras),
		// or these readers are combined with all other readers in readers.txt
		suffixNecessary = true; // false
		suffixSufficient = true; // false;
///		mergeSubIFDs = true; // false;
		LOGGER.debug("ElphelTiffReader(), after supper(), mergeSubIFDs = true;");
		if (REPLACEMENT_TAG_MAP == null) {
			REPLACEMENT_TAG_MAP = new HashMap<String,String>();
			for (String [] line: REPLACEMENT_TAGS) {
				REPLACEMENT_TAG_MAP.put(line[1], line[0]);
			}
		}
	}

	// -- IFormatReader API methods --

	/* @see loci.formats.IFormatReader#isThisType(String, boolean) */
	@Override
	public boolean isThisType(String name, boolean open) {
		if (open) {
			return super.isThisType(name, open);
		}
		return checkSuffix(name, getSuffixes());
	}

	/* @see loci.formats.IFormatReader#isThisType(RandomAccessInputStream) */
	/*
	 * Check the first few bytes of a file to determine if the file can be read
	 * by this reader. You can assume that index 0 in the stream corresponds to
	 * the index 0 in the file. Return true if the file can be read; false if not
	 * (or if there is no way of checking).
	 */
	@Override
	public boolean isThisType(RandomAccessInputStream stream) throws IOException
	{
		final int blockLen = 4;
		if (!FormatTools.validStream(stream, blockLen, false)) return false;

		byte[] signature = new byte[blockLen];
		stream.read(signature);

		if (signature[0] != (byte) 0xff || signature[1] != (byte) 0xd8 ||
				signature[2] != (byte) 0xff || (signature[3] & 0xf0) == 0)
		{
			return false;
		}
		return true;
	}

	  @Override
	  public void setId(String id) throws FormatException, IOException {
		  LOGGER.debug("setId("+id+"). before super" );
		  file_initialized = false;
		  mapped_externally = false;
		  if (Location.getIdMap().containsKey(id)) {
			  LOGGER.debug("id '"+id+"' is already mapped" );
			  content_fileName = id; // id; // maybe set to null to handle externally?
			  mapped_externally = true;
			  LOGGER.debug("Starting setId() method, read file directly");
			  super.setId(id);
		  } else {
			  // If URL, then read to memory, if normal file - use direct access
			  //	 		String cameraURL = "http://192.168.0.36:2323/img"; // for testing
			  // just testing - ignore file name and use camera URL
			  //	 		id = cameraURL;
			  url = null;
			  //	 		String mime = null; // use to select jp4/tiff later? Or to check it is correct
			  content_fileName = null;
			  boolean bad_file = false;
			  try {
				  File test_file = new File(id);
				  test_file.getCanonicalPath();
			  } catch (IOException e) {
				  bad_file = true;
			  }
			  if (bad_file) {
				  try {
					  url = new URL(id);
				  } catch (MalformedURLException e) {
					  LOGGER.warn("Bad URL2: " + id);
				  }
			  }

			  if (url != null) {
				  LOGGER.debug("Starting initFile() method, read "+ id +" to memory first");
				  //https://www.rgagnon.com/javadetails/java-0487.html
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
				  //	 			currentId = fileName; //???
				  //	 			LOGGER.debug("Mime type = "+mime);
				  // https://stackoverflow.com/questions/2793150/how-to-use-java-net-urlconnection-to-fire-and-handle-http-requests

				  //https://stackoverflow.com/questions/2295221/java-net-url-read-stream-to-byte
				  InputStream is =    url.openStream (); //
				  byte[] inBytes = IOUtils.toByteArray(is);
				  if (is != null) is.close();
				  LOGGER.debug("Bytes read: "+ inBytes.length);
				  Location.mapFile(content_fileName, new ByteArrayHandle(inBytes));
//				  HashMap<String,Object> dbg_loc = Location.getIdMap();
				  super.setId(content_fileName);
			  } else { // read file normally
				  content_fileName = id;
				  LOGGER.debug("read file directly");
				  super.setId(id);
			  }
		  }

		  //getReader
		  //	    super.setId(id);
		  LOGGER.debug("setId("+id+"). after super" );
		  file_initialized = true;
	  }

//	@Override
	protected void initFileOld(java.lang.String id)
            throws FormatException,
                   java.io.IOException
    {
		// Trying ServiceFactory before it is going to be initialized, so static defaultFactory will be initialized
		// with small set of services - only needed for Elphel
		LOGGER.debug("Starting initFile() method");
		super.initFile(id);
		LOGGER.debug("Ending initFile() method");
    }

	/* @see loci.formats.FormatReader#initFile(String) */
	// copied from ElphelJp4Reader
	@Override
	protected void initFile(String id) throws FormatException, IOException {
		LOGGER.debug("initFile("+id+"), currentId="+currentId+",  before super" );
		try {
			super.initFile(id); // fails class_not_found
		}
		catch (IllegalArgumentException e) {
			throw new FormatException(e);
		}
		LOGGER.debug("initFile("+id+"), currentId="+currentId+",  after super" );
		// Below needs to be modified - EXIFService does not work with mapFile
		MetadataStore store = makeFilterMetadata();
		LOGGER.debug("Parsing TIFF EXIF data");
		HashMap<String, String> tags = null;
		try {
			// Reimplementing ExifServiceImpl as original does not have ExifIFD0Directory
		    try (RandomAccessInputStream jpegFile = new RandomAccessInputStream(id)) {
		        try {
		          Metadata metadata = ImageMetadataReader.readMetadata(jpegFile);
		          directory =      metadata.getFirstDirectoryOfType(ExifSubIFDDirectory.class);
		          directory_ifd0 = metadata.getFirstDirectoryOfType(ExifIFD0Directory.class);
		          directory_gps =  metadata.getFirstDirectoryOfType(GpsDirectory.class);
		        }
		        catch (Throwable e) {
		          throw new ServiceException("Could not read EXIF data", e);
		        }
		      }
		    if (directory==null) { // trying to read ImageJ file 640x512
		    	return;
		    }
		    Date date =  directory.getDate(ExifSubIFDDirectory.TAG_DATETIME_ORIGINAL);

		    if (date != null) {
				Timestamp timestamp = new Timestamp(new DateTime(date));
				store.setImageAcquisitionDate(timestamp, 0);
			}

			tags = new HashMap<String, String>();
			if (directory != null) {
				for (Tag tag : directory.getTags()) {
					String tag_name = tag.getTagName();
					if (REPLACEMENT_TAG_MAP.containsKey(tag_name)) {
						tags.put(REPLACEMENT_TAG_MAP.get(tag_name), tag.getDescription());
					} else {
						tags.put(tag.getTagName(), tag.getDescription());
					}
				}
			}
			if (directory_ifd0 != null) {
				for (Tag tag : directory_ifd0.getTags()) {
					String tag_name = tag.getTagName();
					if (REPLACEMENT_TAG_MAP.containsKey(tag_name)) {
						tags.put(REPLACEMENT_TAG_MAP.get(tag_name), tag.getDescription());
					} else {
						tags.put(tag.getTagName(), tag.getDescription());
					}
				}
			}
			if (directory_gps != null) {
				for (Tag tag : directory_gps.getTags()) {
					String tag_name = tag.getTagName();
					if (REPLACEMENT_TAG_MAP.containsKey(tag_name)) {
						tags.put(REPLACEMENT_TAG_MAP.get(tag_name), tag.getDescription());
					} else {
						tags.put(tag.getTagName(), tag.getDescription());
					}
				}
			}
			// remove "sec" from exposure
			if (tags.containsKey(EXPOSURE_TIME)){
				tags.put(EXPOSURE_TIME, tags.get(EXPOSURE_TIME).split(" ")[0]);
			}

			for (String tagName : tags.keySet()) {
				addGlobalMeta(tagName, tags.get(tagName));
			} //{Makernote=105455 131072 127570 300581 171508736 171508736 171508736 171508736 169869312 124780556 1118544 0 0 327779 648 1296, Sub-Sec Time Original=560439, Exposure Time=11167/500000 sec, Date/Time Original=2019:05:13 04:30:26}

		}
		catch (ServiceException e) {
			LOGGER.debug("Could not parse EXIF data", e);
		}
		long [] maker_note = null;
		double exposure = Double.NaN;
		String date_time = null;
		if (tags.containsKey(MAKER_NOTE)){
			String [] smn = tags.get(MAKER_NOTE).split(" ");
			maker_note = new long[smn.length];
			for (int i = 0; i < maker_note.length; i++) {
				maker_note[i] = Integer.parseInt(smn[i]);
			}
		}
		if (tags.containsKey(EXPOSURE_TIME)){
			if (tags.get(EXPOSURE_TIME).contains("/")) {
				String [] s = tags.get(EXPOSURE_TIME).split("/");
				exposure = 1.0 * Integer.parseInt(s[0]) / Integer.parseInt(s[1].split(" ")[0]);
			} else {
				exposure = Double.parseDouble(tags.get(EXPOSURE_TIME));
			}
		}
		if (tags.containsKey(DATE_TIME_ORIGINAL)){
			date_time = tags.get(DATE_TIME_ORIGINAL);
			if (tags.containsKey(SUB_SEC_TIME_ORIGINAL)){
				date_time += "."+tags.get(SUB_SEC_TIME_ORIGINAL);
			}
		}
//		int bytes_per_pixel =  1;
		int bpp =   getBitsPerPixel();
		int bytes_per_pixel =  (bpp + 7) / 9;
		Hashtable<String, String> property_table = ElphelMeta.getMeta(
				null, maker_note, exposure, date_time, bytes_per_pixel, true );
		LOGGER.debug("Created elphelMeta table, size="+property_table.size());
		for (String key:property_table.keySet()) {
			addGlobalMeta(ELPHEL_PROPERTY_PREFIX+key,property_table.get(key));
		}
		MetadataLevel level = getMetadataOptions().getMetadataLevel();
		if (level != MetadataLevel.MINIMUM) {
			//			Integer[] tags = ifds.get(0).keySet().toArray(new Integer[0]);
			//			LOGGER.debug("initStandardMetadata() - got "+tags.length+" tags");
		}
		addGlobalMeta(ELPHEL_PROPERTY_PREFIX+CONTENT_FILENAME,content_fileName);
	}
	
	
	/* @see loci.formats.IFormatReader#getSeriesUsedFiles(boolean) */
	@Override
	public String[] getSeriesUsedFiles(boolean noPixels) {
		return super.getSeriesUsedFiles(noPixels);
		//	    if (noPixels) {
		//	      return companionFile == null ? null : new String[] {companionFile};
		//	    }
		//	    if (companionFile != null) return new String[] {companionFile, currentId};
		//	    return new String[] {currentId};
	}

	/* @see loci.formats.IFormatReader#close(boolean) */
	@Override
	public void close(boolean fileOnly) throws IOException {
		LOGGER.debug("close("+fileOnly+") before super");
		super.close(fileOnly); // curerent_id == null only during actual close?
		LOGGER.debug("close("+fileOnly+") after super");
//		if ((content_fileName != null) && file_initialized){
		if (!mapped_externally && file_initialized){ // will try to unmap non-mapped file, OK
			Location.mapFile(content_fileName, null);
			file_initialized = false;
		}
//	    HashMap<String,Object> dbg_loc = Location.getIdMap();

		if (!fileOnly) {
			//	      companionFile = null;
			//	      description = null;
			//	      calibrationUnit = null;
			//	      physicalSizeZ = null;
			//	      timeIncrement = null;
			//	      xOrigin = null;
			//	      yOrigin = null;
		}
	}

	// -- Internal BaseTiffReader API methods --

	/* @see BaseTiffReader#initStandardMetadata() */
	/* Removed 08/03/2021 to match IP4
	@Override
	protected void initStandardMetadata() throws FormatException, IOException {
		LOGGER.debug("initStandardMetadata() - before super()");
		super.initStandardMetadata();
		String comment = ifds.get(0).getComment(); // IMAGE_DESCRIPTION
		LOGGER.debug("initStandardMetadata() - after super()");
		long[] maker_note = null;
		double exposure = Double.NaN;
		String date_time = null;
		IFDList exifIFDs = tiffParser.getExifIFDs();
		//ifds.get(0) has 37393=297894 (ImageNumbder), but not exifIFD
		IFD ifd0 = ifds.get(0);
		if (ifd0.containsKey(TAG_IMAGE_NUMBER))  addGlobalMeta("Image_Number",          ifd0.get(TAG_IMAGE_NUMBER));
		if (ifd0.containsKey(TAG_SERIAL_NUMBER)) addGlobalMeta("Serial_Number",         ifd0.get(TAG_SERIAL_NUMBER));
		if (ifd0.containsKey(TAG_IMAGE_DESCRIPTION)) addGlobalMeta("Image_Description", ifd0.get(TAG_IMAGE_DESCRIPTION));

		//TAG_SERIAL_NUMBER

		if (exifIFDs.size() > 0) {
			IFD exifIFD = exifIFDs.get(0);
			tiffParser.fillInIFD(exifIFD);
			if (exifIFD.containsKey(IFD.MAKER_NOTE)) {
				maker_note = (long[]) exifIFD.get(IFD.MAKER_NOTE);
			}
			if (exifIFD.containsKey(IFD.EXPOSURE_TIME)) {
				Object exp = exifIFD.get(IFD.EXPOSURE_TIME);
				if (exp instanceof TiffRational) {
					TiffRational texp = (TiffRational) exp;
					exposure = 1.0*texp.getNumerator()/texp.getDenominator();
				}
			}
			if (exifIFD.containsKey(IFD.DATE_TIME_ORIGINAL)) {
				date_time = exifIFD.get(IFD.DATE_TIME_ORIGINAL).toString();
				if (exifIFD.containsKey(IFD.SUB_SEC_TIME_ORIGINAL)) {
					date_time += "."+exifIFD.get(IFD.SUB_SEC_TIME_ORIGINAL).toString();
				}
			}
		}
		int bpp =   getBitsPerPixel();
		int bytes_per_pixel =  (bpp + 7) / 9;
		Hashtable<String, String> property_table = ElphelMeta.getMeta(
				null, maker_note, exposure, date_time, bytes_per_pixel, true );


		LOGGER.debug("Created elphelMeta table, size="+property_table.size());
		for (String key:property_table.keySet()) {
			addGlobalMeta(ELPHEL_PROPERTY_PREFIX+key,property_table.get(key));
		}
		MetadataLevel level = getMetadataOptions().getMetadataLevel();
		if (level != MetadataLevel.MINIMUM) {
			Integer[] tags = ifds.get(0).keySet().toArray(new Integer[0]);
			LOGGER.debug("initStandardMetadata() - got "+tags.length+" tags");
	    }
		addGlobalMeta(ELPHEL_PROPERTY_PREFIX+CONTENT_FILENAME,content_fileName);
		// convert MAKER_NOTE to the same text format as in com.drew.metadata
		if (maker_note != null) {
			StringBuffer sb_maker_note = new StringBuffer();
			for (int i = 0; i < maker_note.length; i++) {
				if (i > 0) sb_maker_note.append(" ");
				sb_maker_note.append(maker_note[i]);
			}
			addGlobalMeta("MAKER_NOTE", sb_maker_note.toString());
		}

		// check for ImageJ-style TIFF comment
		boolean ij = checkCommentImageJ(comment);
//		if (ij) parseCommentImageJ(comment);
	}
	*/
	  /* @see BaseTiffReader#initMetadataStore() */
	/* Removed 08/03/2021 to match IP4
	  @Override
	protected void initMetadataStore() throws FormatException {
	    super.initMetadataStore();
	    MetadataStore store = makeFilterMetadata();
//	    if (description != null) {
//	      store.setImageDescription(description, 0);
//	    }
	    populateMetadataStoreImageJ(store);
	  }
	 */
	  /**
	   * @see loci.formats.IFormatReader#openBytes(int, byte[], int, int, int, int)
	   */
	  @Override
	  public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h)
	    throws FormatException, IOException
	  {
		LOGGER.debug("openBytes() - before super()");
	    super.openBytes(no, buf, x, y, w, h);
		LOGGER.debug("openBytes() - after super()");
	    return buf;
	  }




	  // -- Helper methods --

	  // Convert to Elphel-specific parameters
	  /**
	   * Checks the original metadata table for ImageJ-specific information
	   * to propagate into the metadata store.
	   */
	  private void populateMetadataStoreImageJ(MetadataStore store) {
	    // TODO: Perhaps we should only populate the physical Z size if the unit is
	    //       a known, physical quantity such as "micron" rather than "pixel".
	    //       e.g.: if (calibrationUnit.equals("micron"))
		/*
	    if (physicalSizeZ != null) {
	      double zDepth = physicalSizeZ.doubleValue();
	      if (zDepth < 0) zDepth = -zDepth;
	      store.setPixelsPhysicalSizeZ(new PositiveFloat(zDepth), 0);
	    }
	    if (timeIncrement != null) {
	      store.setPixelsTimeIncrement(timeIncrement, 0);
	    }
	    */
	  }
	//
	private boolean checkCommentImageJ(String comment) {
		return comment != null && comment.startsWith("ImageJ=");
	}

}
