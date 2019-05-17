/**
 ** -----------------------------------------------------------------------------**
 ** LwirReader.java
 **
 ** Multithreaded reading cameras - visible and LWIR
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  LwirReader.java is free software: you can redistribute it and/or modify
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
package com.elphel.imagej.lwir;

import java.io.IOException;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.elphel.imagej.readers.ImagejJp4TiffMulti;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import loci.formats.FormatException;

public class LwirReader {
	public static String [] BASE_URLS= {
			"http://192.168.0.36:2323/bchn0",
			"http://192.168.0.36:2324/bchn0",
			"http://192.168.0.36:2325/bchn0",
			"http://192.168.0.36:2326/bchn0",
			"http://192.168.0.38:2323/bchn4",
			"http://192.168.0.38:2324/bchn4",
			"http://192.168.0.38:2325/bchn4",
			"http://192.168.0.38:2326/bchn4",
			};

	public static String IMAGE_URL_FIRST="/towp/wait/img/save";     // next image name, same image number/time
	public static String IMAGE_URL_NEXT= "/torp/next/wait/img/save"; // same image name, same image number/time

//	public static String IMAGE_URL_FIRST="/towp/wait/bimg/save";     // next image name, same image number/time
//	public static String IMAGE_URL_NEXT= "/torp/next/wait/bimg/save"; // same image name, same image number/time

	/** Logger for this class. */
	private static final Logger LOGGER =
			LoggerFactory.getLogger(LwirReader.class);

	private ImagejJp4TiffMulti imagejJp4TiffMulti;





	public LwirReader() {
		imagejJp4TiffMulti = null;
	}


	public ImagePlus[][] readAllMultiple(
			final int     num_frames,
			final boolean show,
			final boolean scale) {
		return readAllMultiple(num_frames, show, scale, "STD_");
	}

	public ImagePlus[][] readAllMultiple(
			final int     num_frames,
			final boolean show,
			final boolean scale,
			final String  std)
	{
		String [] urls0 = new String[BASE_URLS.length];
		String [] urls1 = new String[BASE_URLS.length];
		ImagePlus [][] imps = new ImagePlus[num_frames][BASE_URLS.length];
		for (int i = 0; i < BASE_URLS.length; i++) {
			urls0[i] = BASE_URLS[i] + IMAGE_URL_FIRST;
			urls1[i] = BASE_URLS[i] + IMAGE_URL_NEXT;
		}
		if (imagejJp4TiffMulti == null) {
			imagejJp4TiffMulti = new ImagejJp4TiffMulti();
			for (int i = 0; i < 2; i++) {
				try {
					imagejJp4TiffMulti.getMultiImages(urls0, imps[0],  scale,  std);
				} catch (IOException e) {
					LOGGER.error("readAllMultiple0: IOException, priming" );
				} catch (FormatException e) {
					LOGGER.error("readAllMultiple0:FormatException, priming");
				}
				LOGGER.info("priming..."+(i+1));
				try {
					Thread.sleep(1000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		for (int n = 0; n < num_frames; n++) {
			LOGGER.info("---- Acquiring frame set "+n);
			try {
				imagejJp4TiffMulti.getMultiImages( (n==0)? urls0:urls1, imps[n],  scale,  std);
			} catch (IOException e) {
				LOGGER.error("readAllMultiple0: IOException, n = " + n);
			} catch (FormatException e) {
				LOGGER.error("readAllMultiple0:FormatException, n = " + n);
			}
		}
		double [][] img_seconds = new double [num_frames][BASE_URLS.length];
		int [][]    img_numbers = new int [num_frames][BASE_URLS.length];
		String [][] img_names = new String [num_frames][BASE_URLS.length];
		for (int n = 0; n < num_frames; n++) {
			for (int i = 0; i < imps[n].length; i++) {
				String dt = null;
				try {
					dt = (String) imps[n][i].getProperty("DATE_TIME");
				} catch (Exception e) {
					LOGGER.error("Something is wrong!");
				}

				img_seconds[n][i] = Double.parseDouble(dt.substring(dt.lastIndexOf(":")+1));
				try {
					img_numbers[n][i] = Integer.parseInt((String) imps[n][i].getProperty("STD_Image_Number"));
				} catch (Exception e) {
					img_numbers[n][i] = -1;
				}
				img_names[n][i] = (String) imps[n][i].getProperty("CONTENT_FILENAME");
				LOGGER.info("Seconds for" + n+":"+i+" - "+img_seconds[n][i]+", number"+img_numbers[n][i]+", name "+img_names[n][i]);
			}
		}


		for (ImagePlus [] imp_sets :imps) {
			for (ImagePlus imp: imp_sets) {
//				imp.show();
			}
		}
		return imps;
	}

	public ImagePlus [] averageMultiFrames(ImagePlus [][] sets) {
		int num_frames = sets.length;
		int num_channels = sets[0].length;
		ImagePlus [] imps_avg = new ImagePlus [num_channels];

		for (int chn = 0; chn < num_channels; chn++) {
			int width =  sets[0][chn].getWidth();
			int height = sets[0][chn].getHeight();
			String title = sets[0][chn].getTitle()+"_average"+num_frames;
			float [] pixels_avg = (float []) sets[0][chn].getProcessor().getPixels(); //null pointer
			for (int n = 1; n < num_frames; n++) {
				float [] pixels = (float []) sets[n][chn].getProcessor().getPixels();
				for (int i = 0; i < pixels_avg.length; i++) {
					pixels_avg[i] += pixels[i];
				}
			}
			double scale = 1.0/num_frames;
			for (int i = 0; i < pixels_avg.length; i++) {
				pixels_avg[i] *= scale;
			}
			ImageProcessor ip=new FloatProcessor(width,height);
			ip.setPixels(pixels_avg);
			ip.resetMinAndMax();
			imps_avg[chn]=  new ImagePlus(title, ip);
			Properties properties0 = sets[0][chn].getProperties();
			for (String key:properties0.stringPropertyNames()) {
				imps_avg[chn].setProperty(key, properties0.getProperty(key));
			}
			// TODO: Overwrite some properties?
		}

		return imps_avg;
	}


	public ImagePlus [][] matchSets(ImagePlus [][] sets, double max_mismatch, int max_frame_diff){
		int num_frames = sets.length;
		int num_channels = sets[0].length;
		double [][] img_seconds = new double [num_frames][num_channels];
		int [] frame_offsets = new int[num_channels];
		double [] time_offsets = new double [num_channels];
		boolean some_notsynced = false;
		int fr_min = 0, fr_max = 0;
		for (int n = 0; n < num_frames; n++) {
			for (int i = 0; i < num_channels; i++) {
				String dt = (String) sets[n][i].getProperty("DATE_TIME");
				img_seconds[n][i] = Double.parseDouble(dt.substring(dt.lastIndexOf(":")+1));

			}
		}
		for (int i = 1; i < num_channels; i++) { // compare all to the first channel
			frame_offsets[i] = 0;
			time_offsets[i] = secOffs(img_seconds[0][0], img_seconds[0][i]);
			if (num_frames > max_frame_diff) {
				if (time_offsets[i] > max_mismatch) {
					for (int fr_diff = 1; fr_diff < max_frame_diff; fr_diff++) {
						double aoff = secOffs(img_seconds[fr_diff][0], img_seconds[0][i]);
						if (aoff < time_offsets[i]) {
							time_offsets[i] = aoff;
							frame_offsets[i] = fr_diff; // second channel starts later, than [0]
						}
						if (time_offsets[i] <=max_mismatch) break;
						aoff = secOffs(img_seconds[0][0], img_seconds[fr_diff][i]);
						if (aoff < time_offsets[i]) {
							time_offsets[i] = aoff;
							frame_offsets[i] = -fr_diff; // second channel starts earlier, than [0]
						}
						if (time_offsets[i] <=max_mismatch) break;
					}
				}
			}
			if (time_offsets[i] > max_mismatch) {
				some_notsynced = true;
			}
			if (frame_offsets[i] > 0) {
				fr_max = frame_offsets[i];
			}
			if (frame_offsets[i] < 0) {
				fr_min = frame_offsets[i];
			}
		}
		if (some_notsynced) {
			LOGGER.error("*** Some channels are not synchronized, reboot or sensors re-start is needed ***");
			for (int i = 0; i < num_channels; i++) {
				LOGGER.error("Channel "+ i+" frame offset="+frame_offsets[i]+ ", time offset = "+time_offsets[i]+" sec");
			}
			return null;
		}
		if (((fr_max - fr_min) > max_frame_diff) || ((fr_max - fr_min) > (num_frames - 1))) {
			LOGGER.error("*** Earliest/latest channels differ by more than 1 frame, that should not happen! ***");
			for (int i = 0; i < num_channels; i++) {
				LOGGER.error("Channel "+ i+" frame offset="+frame_offsets[i]+ ", time offset = "+time_offsets[i]+" sec");
			}
			return null;
		}
		for (int i = 0; i < num_channels; i++) {
			// change to info later:
			LOGGER.info("Channel "+ i+" frame offset="+frame_offsets[i]+ ", time offset = "+time_offsets[i]+" sec");
		}
		ImagePlus [][] imps_synced = new ImagePlus [num_frames - fr_max + fr_min][num_channels];
		for (int n = 0; n < imps_synced.length; n++) {
			imps_synced[n][0]= sets[n+fr_max][0];
			for (int i = 1; i < num_channels; i++) {
				imps_synced[n][i]= sets[n - fr_min][i];
			}
		}
		return imps_synced;
	}

	private double secOffs(double sec1, double sec2) {
		double aoff = Math.abs(sec2 - sec1);
		if (aoff > 30) {
			aoff = Math.abs(aoff - 30);
		}
		return aoff;
	}

}
