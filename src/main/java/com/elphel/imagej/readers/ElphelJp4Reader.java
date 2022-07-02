/**
 ** -----------------------------------------------------------------------------**
 ** ElphelJp4Reader.java
 **
 ** loci.format compatible reader for Elphel JP4 files
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ElphelJp4Reader.java is free software: you can redistribute it and/or modify
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
import loci.formats.in.ImageIOReader;
import loci.formats.in.MetadataLevel;
//import loci.formats.services.EXIFService;
import ome.xml.meta.MetadataStore;
import ome.xml.model.primitives.Timestamp;

//ElphelTiffReader
public class ElphelJp4Reader  extends ImageIOReader{
	// -- Constants --
	public static final String MAKER_NOTE =            "Makernote";
	public static final String SUB_SEC_TIME_ORIGINAL = "SUBSEC_TIME_ORIGINAL"; // "Sub-Sec Time Original";
	public static final String EXPOSURE_TIME =         "Exposure Time";
	public static final String DATE_TIME_ORIGINAL =    "DATE_TIME_ORIGINAL"; // "Date/Time Original";

	public static final String ELPHEL_PROPERTY_PREFIX = "ELPHEL_";
	public static final String CONTENT_FILENAME =       "CONTENT_FILENAME";
	public static final boolean REORDER = true; // false;
	public static final String[][] REPLACEMENT_TAGS = // to/from! TODO: check/add for GPS?
		   {{"SUBSEC_TIME_ORIGINAL",      "Sub-Sec Time Original"},
			{"DATE_TIME_ORIGINAL",        "Date/Time Original"},
			{"Instrument_Make",           "Make"},
			{"Serial_Number",             "Unknown tag (0xc62f)"},
			{"Instrument_Model",          "Model"}};





	/** Logger for this class. */
	private static final Logger LOGGER =
			LoggerFactory.getLogger(ElphelTiffReader.class);

	// -- Fields --
	private URL url = null; // save here actual URL when reading file to memory
	private String content_fileName = null; // from Content-disposition
	private boolean mapped_externally = false; // file is read/mapped externally, do not close it here
	private boolean file_initialized = false;
	private byte [] image_bytes = null;
	private ExifSubIFDDirectory directory;
	private ExifIFD0Directory   directory_ifd0;
	private GpsDirectory        directory_gps;
	private HashMap<String,String> REPLACEMENT_TAG_MAP = null; // per instance


	// -- Constructor --

	/** Constructs a new Tiff reader. */
	public ElphelJp4Reader() {
		super("JP4", new String[] {"jp4"});
		//		mergeSubIFDs = true; // false;
		// TODO: See if the selection is just between 2 readers (jp4 and tiff - just Elphel cameras),
		// or these readers are combined with all other readers in readers.txt
		suffixNecessary = true; // false
		suffixSufficient = true; // false;
		LOGGER.info("ElphelTiffReader(), after super()");
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
	/* https://docs.openmicroscopy.org/bio-formats/5.7.2/developers/reader-guide.html
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
	public void setId(String id) throws FormatException, IOException { // same as for tiff?
		image_bytes = null;
//		buffered_data = null;
		LOGGER.info("setId("+id+"). before super" );
		file_initialized = false;
		mapped_externally = false;

		if (Location.getIdMap().containsKey(id)) {
			LOGGER.info("id '"+id+"' is already mapped" );
			content_fileName = id; // id; // maybe set to null to handle externally?
			mapped_externally = true;
			LOGGER.info("Starting initFile() method, read file directly");
			super.setId(id);
		} else {
			// If URL, then read to memory, if normal file - use direct access
			url = null;
			//	 		String mime = null; // use to select jp4/tiff later? Or to check it is correct
			content_fileName = null;
			try {
				url = new URL(id);
			} catch (MalformedURLException e) {
//				LOGGER.warn("Bad URL1: " + id); // will try direct file, not an error
			}
			if (url != null) {
				LOGGER.info("Starting initFile() method, read "+ id +" to memory first");
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
				//	 			LOGGER.info("Mime type = "+mime);
				// https://stackoverflow.com/questions/2793150/how-to-use-java-net-urlconnection-to-fire-and-handle-http-requests

				//https://stackoverflow.com/questions/2295221/java-net-url-read-stream-to-byte
				InputStream is =    url.openStream (); //
				byte[] inBytes = IOUtils.toByteArray(is);
				if (is != null) is.close();
				LOGGER.info("Bytes read: "+ inBytes.length);
				Location.mapFile(content_fileName, new ByteArrayHandle(inBytes));
//				HashMap<String,Object> dbg_loc = Location.getIdMap();
				super.setId(content_fileName);
			} else { // read file normally
				content_fileName = id;
				LOGGER.info("read file directly");
				super.setId(id);
			}
		}

		//getReader
		//	    super.setId(id);
		LOGGER.info("setId("+id+"). after super" );
		file_initialized = true;
	}

	/* @see loci.formats.FormatReader#initFile(String) */
	@Override
	protected void initFile(String id) throws FormatException, IOException {
		LOGGER.info("initFile("+id+"), currentId="+currentId+",  before super" );
		try {
			super.initFile(id); // fails class_not_found
		}
		catch (IllegalArgumentException e) {
			throw new FormatException(e);
		}
		LOGGER.info("initFile("+id+"), currentId="+currentId+",  after super" );
		// Below needs to be modified - EXIFService does not work with mapFile
		MetadataStore store = makeFilterMetadata();
		LOGGER.info("Parsing JPEG EXIF data");
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
			LOGGER.info("Could not parse EXIF data", e);
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
		int bytes_per_pixel =  1;
		Hashtable<String, String> property_table = ElphelMeta.getMeta(
				null, maker_note, exposure, date_time, bytes_per_pixel, true );
		LOGGER.info("Created elphelMeta table, size="+property_table.size());
		for (String key:property_table.keySet()) {
			addGlobalMeta(ELPHEL_PROPERTY_PREFIX+key,property_table.get(key));
		}
		MetadataLevel level = getMetadataOptions().getMetadataLevel();
		if (level != MetadataLevel.MINIMUM) {
			//			Integer[] tags = ifds.get(0).keySet().toArray(new Integer[0]);
			//			LOGGER.info("initStandardMetadata() - got "+tags.length+" tags");
		}
		addGlobalMeta(ELPHEL_PROPERTY_PREFIX+CONTENT_FILENAME,content_fileName);
	}


	/* @see loci.formats.IFormatReader#close(boolean) */
	@Override
	public void close(boolean fileOnly) throws IOException {
//		HashMap<String,Object> dbg_loc = Location.getIdMap();
		String saveCurrentId = currentId;
		currentId = null;
		LOGGER.info("close("+fileOnly+") before super");
		super.close(fileOnly); // curerent_id == null only during actual close?
		LOGGER.info("close("+fileOnly+") after super");
		currentId = saveCurrentId;
//		if ((content_fileName != null) && file_initialized){
		if (!mapped_externally && file_initialized){ // will try to unmap non-mapped file, OK
			Location.mapFile(content_fileName, null);
			file_initialized = false;
		}
//		dbg_loc = Location.getIdMap();

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

	/**
	 * @see loci.formats.IFormatReader#openBytes(int, byte[], int, int, int, int)
	 *
	 * openBytes(int, byte[], int, int, int, int) Returns a byte array containing
	 * the pixel data for a specified subimage from the given file. The dimensions
	 * of the subimage (upper left X coordinate, upper left Y coordinate, width,
	 * and height) are specified in the final four int parameters.
	 * This should throw a FormatException if the image number is invalid (less than
	 * 0 or >= the number of images). The ordering of the array returned by openBytes
	 * should correspond to the values returned by isLittleEndian and isInterleaved.
	 * Also, the length of the byte array should be
	 * [image_width * image height * bytes_per_pixel]. Extra bytes will generally
	 * be truncated. It is recommended that the first line of this method be
	 * FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h) - this ensures
	 * that all of the parameters are valid.
	 */
	@Override
	public byte[] openBytes(int no, byte[] buf, int x, int y, int w, int h)
			throws FormatException, IOException
	{
		LOGGER.info("openBytes() - before super()");
		FormatTools.checkPlaneParameters(this, no, buf.length, x, y, w, h);
		if (image_bytes == null) {
			jp4Decode(no);
		}
		int width =  getSizeX();
		int y1 = y + h;
		int dest=0;
		for (int line = y; line < y1; line++) {
			System.arraycopy(image_bytes,
					         line*width+x,
					         buf,
					         dest,
					         w);
			dest += w;
		}
		LOGGER.info("openBytes() - after super()");
		return buf;
	}
	public void jp4Decode(int no) throws FormatException, IOException {
		int width =  getSizeX();
		int height = getSizeY();
		image_bytes = new byte[width*height];
		byte [] ib =  new byte [width*height];
		byte [][]     macroblock=new byte[16][16];

		super.openBytes(no, ib, 0, 0, width, height);
		if (REORDER) {
			int yb,y,xb,x,offset,nb,xbyr,ybyr; //,i;
			for (yb=0;yb<(height>>4); yb++) for (xb=0;xb<(width>>4); xb++) { /* iterating macroblocks */
				for (nb=0; nb<4;nb++) {
					xbyr=nb & 1;
					ybyr=(nb>>1) & 1;
					for (y=0;y<8;y++) {
						offset=((yb<<4)+y)*width+ (nb<<3) +((xb>=(width>>5))?(((xb<<5)-width)+(width<<3)):(xb<<5));
						for (x=0;x<8;x++) {
							macroblock[(y<<1) | ybyr][(x<<1) | xbyr]=ib[offset+x];
						}
					}
				}
				for (y=0;y<16;y++) {
					System.arraycopy(macroblock[y],
					                 0,
					                 image_bytes,
					                 ((yb<<4)+y)*width + (xb<<4),
					                 16);
				}
			}
		} else {
			image_bytes = ib; // temporary
		}
		LOGGER.info("jp4Decode()");
	}

}
