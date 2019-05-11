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
import java.io.IOException;
import java.util.Hashtable;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import loci.formats.FormatException;
import loci.formats.in.MetadataLevel;
import loci.formats.in.TiffReader;
import loci.formats.meta.MetadataStore;

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

	//	  private String companionFile;
	//	  private String description;
	//	  private String calibrationUnit;
	//	  private Double physicalSizeZ;
	//	  private Double timeIncrement;
	//	  private Integer xOrigin, yOrigin;

	// -- Constructor --

	/** Constructs a new Tiff reader. */
	public ElphelTiffReader() {
		super(); // "Tagged Image File Format", ELPHEL_TIFF_SUFFIXES); // See if we can use TiffReader without its parent
		mergeSubIFDs = true; // false;
		LOGGER.info("ElphelTiffReader(), after supper(), mergeSubIFDs = true;");
	}

	// -- IFormatReader API methods --


	@Override
	protected void initFile(java.lang.String id)
            throws FormatException,
                   java.io.IOException
    {
		LOGGER.info("Starting initFile() method");
		super.initFile(id);
		LOGGER.info("Ending initFile() method");
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
		LOGGER.info("close("+fileOnly+") before super");
		super.close(fileOnly);
		LOGGER.info("close("+fileOnly+") after super");
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
	@Override
	protected void initStandardMetadata() throws FormatException, IOException {
		LOGGER.info("initStandardMetadata() - before super()");
		super.initStandardMetadata();
		String comment = ifds.get(0).getComment(); // IMAGE_DESCRIPTION
		LOGGER.info("initStandardMetadata() - after super()");
		ElphelMeta elphelMeta = new ElphelMeta(tiffParser, true);
		Hashtable<String, String> property_table = elphelMeta.getPropertyTable();
		LOGGER.info("Created elphelMeta table, size="+property_table.size());
		for (String key:property_table.keySet()) {
			addGlobalMeta(ELPHEL_PROPERTY_PREFIX+key,property_table.get(key));
		}
		MetadataLevel level = getMetadataOptions().getMetadataLevel();
		if (level != MetadataLevel.MINIMUM) {
			Integer[] tags = ifds.get(0).keySet().toArray(new Integer[0]);
			LOGGER.info("initStandardMetadata() - got "+tags.length+" tags");
	    }
		// check for ImageJ-style TIFF comment
		boolean ij = checkCommentImageJ(comment);
//		if (ij) parseCommentImageJ(comment);
		/*
		// check for MetaMorph-style TIFF comment
		boolean metamorph = checkCommentMetamorph(comment);
		if (metamorph && level != MetadataLevel.MINIMUM) {
			parseCommentMetamorph(comment);
		}
		put("MetaMorph", metamorph ? "yes" : "no");

		// check for other INI-style comment
		if (!ij && !metamorph && level != MetadataLevel.MINIMUM) {
			parseCommentGeneric(comment);
		}
		 */
		// check for another file with the same name
		/*
		if (isGroupFiles()) {
			Location currentFile = new Location(currentId).getAbsoluteFile();
			String currentName = currentFile.getName();
			Location directory = currentFile.getParentFile();
			String[] files = directory.list(true);
			if (files != null) {
				for (String file : files) {
					String name = file;
					if (name.indexOf(".") != -1) {
						name = name.substring(0, name.indexOf("."));
					}

					if (currentName.startsWith(name) &&
							checkSuffix(name, COMPANION_SUFFIXES))
					{
						companionFile = new Location(directory, file).getAbsolutePath();
						break;
					}
				}
			}
		}
		*/
	}

	  /* @see BaseTiffReader#initMetadataStore() */
	  @Override
	protected void initMetadataStore() throws FormatException {
	    super.initMetadataStore();
	    MetadataStore store = makeFilterMetadata();
//	    if (description != null) {
//	      store.setImageDescription(description, 0);
//	    }
	    populateMetadataStoreImageJ(store);
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
