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

import loci.common.ByteArrayHandle;
import loci.common.Location;
import loci.common.RandomAccessInputStream;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.FormatTools;
import loci.formats.in.ImageIOReader;
import loci.formats.in.MetadataLevel;
import loci.formats.services.EXIFService;
//import loci.formats.services.EXIFService;
import ome.xml.meta.MetadataStore;
import ome.xml.model.primitives.Timestamp;

//ElphelTiffReader
public class ElphelJp4Reader  extends ImageIOReader{
	// -- Constants --
	public static final String MAKER_NOTE =            "Makernote";
	public static final String SUB_SEC_TIME_ORIGINAL = "Sub-Sec Time Original";
	public static final String EXPOSURE_TIME =         "Exposure Time";
	public static final String DATE_TIME_ORIGINAL =    "Date/Time Original";

	public static final String ELPHEL_PROPERTY_PREFIX = "ELPHEL_";
	public static final String CONTENT_FILENAME = "CONTENT_FILENAME";

	/** Logger for this class. */
	private static final Logger LOGGER =
			LoggerFactory.getLogger(ElphelTiffReader.class);

	// -- Fields --
	private URL url = null; // save here actual URL when reading file to memory
	private String content_fileName = null; // from Content-disposition
	private boolean file_initialized = false;
	private ElphelTiffReader elphelTiffReader = null;
//	private ImageReader helperReader;


	// -- Constructor --

	/** Constructs a new Tiff reader. */
	public ElphelJp4Reader() {
		super("JP4", new String[] {"jp4"});
		//		mergeSubIFDs = true; // false;
		suffixNecessary = true; // false
		suffixSufficient = true; // false;

		LOGGER.info("ElphelTiffReader(), after super()");
		elphelTiffReader = new ElphelTiffReader();
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
		LOGGER.debug("setId("+id+"). before super" );
		file_initialized = false;
		if (Location.getIdMap().containsKey(id)) {
			LOGGER.debug("id '"+id+"' is already mapped" );
			content_fileName = null; // id; // maybe set to null to handle externally?
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
				LOGGER.warn("Bad URL: " + id);
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
				HashMap<String,Object> dbg_loc = Location.getIdMap();
				super.setId(content_fileName);
			} else { // read file normally
				content_fileName = id;
				LOGGER.info("read file directly");
				super.setId(id);
			}
		}

		//getReader
		//	    super.setId(id);
		LOGGER.debug("setId("+id+"). after super" );
		file_initialized = true;
	}

	/* @see loci.formats.FormatReader#initFile(String) */
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

//		elphelTiffReader.setId(id); //error here - invalid tiff file


		// Below needs to be modified - EXIFService does not work with mapFile

		MetadataStore store = makeFilterMetadata();
		LOGGER.info("Parsing JPEG EXIF data");
		HashMap<String, String> tags = null;
		try {
	        EXIFService exif = new ServiceFactory().getInstance(EXIFService.class);
	        if (exif == null) {
				return;
			}
			exif.initialize(id);

			// Set the acquisition date
			Date date = exif.getCreationDate();
			if (date != null) {
				Timestamp timestamp = new Timestamp(new DateTime(date));
				store.setImageAcquisitionDate(timestamp, 0);
			}

			tags = exif.getTags();
			for (String tagName : tags.keySet()) {
				addGlobalMeta(tagName, tags.get(tagName));
			} //{Makernote=105455 131072 127570 300581 171508736 171508736 171508736 171508736 169869312 124780556 1118544 0 0 327779 648 1296, Sub-Sec Time Original=560439, Exposure Time=11167/500000 sec, Date/Time Original=2019:05:13 04:30:26}
		}
		catch (ServiceException e) {
			LOGGER.debug("Could not parse EXIF data", e);
		}
		catch (DependencyException e) {
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
			String [] s = tags.get(EXPOSURE_TIME).split("/");
			exposure = 1.0 * Integer.parseInt(s[0]) / Integer.parseInt(s[1].split(" ")[0]);
		}
		if (tags.containsKey(DATE_TIME_ORIGINAL)){
			date_time = tags.get(DATE_TIME_ORIGINAL);
			if (tags.containsKey(SUB_SEC_TIME_ORIGINAL)){
				date_time += tags.get(SUB_SEC_TIME_ORIGINAL);
			}
		}
		Hashtable<String, String> property_table = ElphelMeta.getMeta(
				null, maker_note, exposure, date_time, true );
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
		HashMap<String,Object> dbg_loc = Location.getIdMap();
		String saveCurrentId = currentId;
		currentId = null;
		LOGGER.info("close("+fileOnly+") before super");
		super.close(fileOnly); // curerent_id == null only during actual close?
		LOGGER.info("close("+fileOnly+") after super");
		currentId = saveCurrentId;
		if ((content_fileName != null) && file_initialized){
			Location.mapFile(content_fileName, null);
			file_initialized = false;
		}
		dbg_loc = Location.getIdMap();

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






}
