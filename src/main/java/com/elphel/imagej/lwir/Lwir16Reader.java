/**
 ** -----------------------------------------------------------------------------**
 ** Lwir16Reader.java
 **
 ** Multithreaded reading cameras - visible and LWIR
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Lwir16Reader.java is free software: you can redistribute it and/or modify
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

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import com.elphel.imagej.readers.ImagejJp4Tiff;
import com.elphel.imagej.readers.ImagejJp4TiffMulti;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import loci.formats.FormatException;

public class Lwir16Reader {
	public static String [] BASE_URLS= {
			"http://192.168.0.41:2323/bchn0",
			"http://192.168.0.41:2324/bchn0",
			"http://192.168.0.41:2325/bchn0",
			"http://192.168.0.41:2326/bchn0",
			"http://192.168.0.42:2323/bchn4",
			"http://192.168.0.42:2324/bchn4",
			"http://192.168.0.42:2325/bchn4",
			"http://192.168.0.42:2326/bchn4",
			"http://192.168.0.43:2323/bchn8",
			"http://192.168.0.43:2324/bchn8",
			"http://192.168.0.43:2325/bchn8",
			"http://192.168.0.43:2326/bchn8",
			"http://192.168.0.44:2323/bchn12",
			"http://192.168.0.44:2324/bchn12",
			"http://192.168.0.44:2325/bchn12",
			"http://192.168.0.44:2326/bchn12",
			"http://192.168.0.45:2323/bchn16",
			"http://192.168.0.45:2324/bchn16",
			"http://192.168.0.45:2325/bchn16",
			"http://192.168.0.45:2326/bchn16",
			};

	public static final int [] IMGSRV_PORTS = {2323,2324,2325,2326};
	public static final String IMAGE_URL_STATUS="/pointers"; // get current status
	public static final String IMAGE_URL_RESET= "/towp/save/xframe";       // reset read pointer while compressors are stopped
//	public static final String IMAGE_URL_META=  "/torp/next/wait/meta";      // wait for frame at current read poiter and read meta 
//	public static final String IMAGE_URL_NEXT=  "/torp/next/wait/img/save"; // same image name, same image number/time
	public static final String IMAGE_URL_NEXT=  "/torp/wait/img/next/save"; // same image name, same image number/time
	public static final String SCRIPT_CAPTURE=  "capture_range.php"; // capture synchronized images
	public static final String SCRIPT_RESET=    "reset_frames.php";  // Reset frame numbers
	public static final String SCRIPT_WAIT=     "wait_frame.php";    // wait for absolute frame (positive) or skip frames (nefative)
	// will remove later
	public static final String IMAGE_URL_FIRST="/towp/wait/img/save";     // next image name, same image number/time
	public static final String SKIP_FRAME_URL= "/towp/wait/meta";          // Wait for the next frame, return meta data
	
	public static final String REG_FFC_FRAMES= "SENSOR_REGS4";      // Register for the number of FFC frames to integrate
	public static final String REG_FFC_RUN=    "SENSOR_REGS26";     // Register to trigger FFC
	
//	public static final String SKIP_FRAME_URL=  "/towp/wait/meta";          // Wait for the next frame, return meta data
	public static final int COLOR_JP4 = 5;
	public static final int COLOR_RAW = 15;
	public static final int LWIR_HEIGHT = 512;
	public static final int LWIR_TELEMETRY_LINES = 1;
	public static final int FRAMES_SKIP = 4;
	public static final int MAX_THREADS = 100; // combine from all classes?

	public static final int FRAMES_AHEAD = 5;

	/** Logger for this class. */
	private static final Logger LOGGER =
			LoggerFactory.getLogger(LwirReader.class);

	private ImagejJp4TiffMulti imagejJp4TiffMulti;
	private boolean out_of_sync = false;

	/** Configuration parameters */
	private LwirReaderParameters lwirReaderParameters = null;
	private LwirReaderParameters last_programmed = null;


	private int [] motorsPosition = null;
	public boolean reportTiming=false;

	// -- constructors

	public Lwir16Reader() {
		this.lwirReaderParameters = new LwirReaderParameters(LwirReaderParameters.NAME_LWIR16); // default
		imagejJp4TiffMulti = null;
	}

	public Lwir16Reader(LwirReaderParameters lwirReaderParameters) {
		if (lwirReaderParameters == null) {
			this.lwirReaderParameters = new LwirReaderParameters(LwirReaderParameters.NAME_LWIR16); // default
		} else {
			this.lwirReaderParameters = lwirReaderParameters;
		}
		imagejJp4TiffMulti = null;
	}


   	public void setMotorsPosition (int [] position) {
   		this.motorsPosition=position;
   	}

   	public int [] getMotorsPosition() {
   		return this.motorsPosition;  // may be null
   	}


	// TODO: create
	public boolean reSyncNeeded() {
		return out_of_sync;
	}

	public ImagePlus[][] readAllMultiple(
			final int     num_frames,
//			final boolean telemetry,
			final boolean show,
			final boolean scale) {
		return readAllMultiple(num_frames, show, scale, "STD_");
//		return readAllMultiple(num_frames, telemetry, show, scale, "STD_");
	}

	public ImagePlus[][] readAllMultiple(
			final int     num_frames,
//			final boolean telemetry,
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
//					imagejJp4TiffMulti.getMultiImages(urls0, imps[0],  telemetry, scale,  std);
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
//				imagejJp4TiffMulti.getMultiImages( (n==0)? urls0:urls1, imps[n], telemetry,  scale,  std);
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
				if (dt==null) {
					System.out.println("n="+n+", i="+i+", dt is NULL! (may need to reduce number of acquire images to "+(n-1));
				}
//				System.out.println("n="+n+", i="+i+", dt="+dt);
				img_seconds[n][i] = Double.parseDouble(dt.substring(dt.lastIndexOf(":")+1));
				try {
					img_numbers[n][i] = Integer.parseInt((String) imps[n][i].getProperty("STD_Image_Number"));
				} catch (Exception e) {
					img_numbers[n][i] = -1;
				}
				img_names[n][i] = (String) imps[n][i].getProperty("CONTENT_FILENAME");
				LOGGER.warn("Seconds for" + n+":"+i+" - "+img_seconds[n][i]+", number"+img_numbers[n][i]+", name "+img_names[n][i]);
			}
		}


		for (ImagePlus [] imp_sets :imps) {
			for (ImagePlus imp: imp_sets) {
//				imp.show();
			}
		}
		return imps;
	}

	public ImagePlus [] averageMultiFrames(
			final ImagePlus [][] sets) {
		final int num_frames = sets.length;
		final int num_channels = sets[0].length;
		final ImagePlus [] imps_avg = new ImagePlus [num_channels];
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);

   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int chn = indxAtomic.getAndIncrement(); chn < num_channels; chn = indxAtomic.getAndIncrement())
   					{
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
   						imps_avg[chn].setProperty("average", ""+num_frames);
   						if (motorsPosition!=null) for (int m=0;m<motorsPosition.length;m++ ) {
   							imps_avg[chn].setProperty("MOTOR"+(m+1), ""+motorsPosition[m]);
   						}
   						ImagejJp4Tiff.encodeProperiesToInfo(imps_avg[chn]);
   					}
   				}
   			};
   		}
   		startAndJoin(threads);

		return imps_avg;
	}

	public ImagePlus [][] matchSets(ImagePlus [][] sets, double max_mismatch, int max_frame_diff){
		return matchSets(sets, max_mismatch, max_frame_diff, null);
	}
	public ImagePlus [][] matchSets(ImagePlus [][] sets, double max_mismatch, int max_frame_diff, int [] lags){
		int num_frames = sets.length;
		int num_channels = sets[0].length;
		double [][] img_seconds = new double [num_frames][num_channels];
		int [] frame_offsets = new int[num_channels];
		double [] time_offsets = new double [num_channels];
		boolean some_notsynced = false;
		int fr_min = 0, fr_max = 0;
		frame_offsets[0] = 0;
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
					for (int fr_diff = 1; fr_diff <= max_frame_diff; fr_diff++) {
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
			if (frame_offsets[i] > fr_max) 	fr_max = frame_offsets[i];
			if (frame_offsets[i] < fr_min)  fr_min = frame_offsets[i];
		}
		if (some_notsynced) {
			LOGGER.error("*** Some channels are not synchronized, reboot or sensors re-start is needed ***");
			out_of_sync = true;
			for (int i = 0; i < num_channels; i++) {
				LOGGER.error("Channel "+ i+" frame offset="+frame_offsets[i]+ ", time offset = "+time_offsets[i]+" sec");
			}
			return null;
		}

		if (((fr_max - fr_min) > max_frame_diff) || ((fr_max - fr_min) > (num_frames - 1))) {
			out_of_sync = true;
			LOGGER.error("*** Earliest/latest channels differ by more than "+max_frame_diff+" frames, that should not happen! ***");
			for (int i = 0; i < num_channels; i++) {
				LOGGER.error("Channel "+ i+" frame offset="+frame_offsets[i]+ ", time offset = "+time_offsets[i]+" sec");
			}
			return null;
		}
		out_of_sync = false;
		for (int i = 0; i < num_channels; i++) {
			// change to info later:
			LOGGER.warn("Channel "+ i+" frame offset="+frame_offsets[i]+ ", time offset = "+time_offsets[i]+" sec");
		}
		// 		if (lags == null) lags = new int [num_channels];

		// recalculate  frame_offsets, fr_max, fr_min considering provided lags (>0 - channel images are acquired later)
		if (lags != null) {
			fr_max = lags[0];
			fr_min = lags[0];
			for (int i = 0; i < num_channels; i++ ) {
				frame_offsets[i] += lags[i];
				if (frame_offsets[i] > fr_max) 	fr_max = frame_offsets[i];
				if (frame_offsets[i] < fr_min)  fr_min = frame_offsets[i];
			}
		}
		ImagePlus [][] imps_synced = null;

		imps_synced = new ImagePlus [num_frames - fr_max + fr_min][num_channels];
		for (int n = 0; n < imps_synced.length; n++) {
			for (int i = 0; i < num_channels; i++) {
				imps_synced[n][i]= sets[n + fr_max -frame_offsets[i]][i];
			}
		}
		return imps_synced;
	}

	public boolean calibrate(LwirReaderParameters lrp) {
		// FFC takes 15 frames that are not output (so no need to check operation is over)
		if (lrp.lwir_channels.length == 0) {
			LOGGER.error("calibrate(): No LWIR channels are configured");
			return false;
		}
		// Skipping frame to synchronize setting calibration command
		if (!skipFrame(lrp)) {
			LOGGER.error("programLWIRCamera():Failed to skip frame");
		}

		int chn = lrp.lwir_channels[0];
		int channels = 0;
		for (int c : lrp.lwir_channels) {
			channels |= 1 << c;
		}
		String hex_chan = String.format("0x%x", channels);

		String url = "http://"+lrp.lwir_ips[0]+"/parsedit.php?immediate&sensor_port="+chn+
				"&SENSOR_REGS67=0*"+FRAMES_AHEAD+"!"+hex_chan;
//		String url = "http://"+lrp.lwir_ip+"/parsedit.php?immediate&sensor_port="+chn+
//				"&SENSOR_REGS67=0&*SENSOR_REGS67="+hex_chan;
		Document dom=null;
		LOGGER.warn("calibrate(): Perform calibration (instead of 15 frames), url="+url);
		try {
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = dbf.newDocumentBuilder();
			dom = db.parse(url);
			if (!dom.getDocumentElement().getNodeName().equals("parameters")) {
				LOGGER.error("programLWIRCamera() in " + url+
						": Root element: expected 'parameters', got \"" + dom.getDocumentElement().getNodeName()+"\"");
				return false;
			}
		} catch(MalformedURLException e){
			LOGGER.error("programLWIRCamera() in " + url+ ": " + e.toString());
			return false;
		} catch(IOException  e1){
			LOGGER.error("programLWIRCamera() in " + url+ " - camera did not respond: " + e1.toString());
			return false;
		}catch(ParserConfigurationException pce) {
			LOGGER.error("programLWIRCamera() in " + url+ " - PCE error: " + pce.toString());
			return false;
		}catch(SAXException se) {
			LOGGER.error("programLWIRCamera() in " + url+ " - SAX error: " + se.toString());
			return false;
		}
		for (int i = 0; i < (FRAMES_SKIP+FRAMES_AHEAD); i++) {
			if (!skipFrame(lrp)) {
				LOGGER.error("programLWIRCamera():Failed to skip frame");
			}
		}

		return true;
	}


	private double secOffs(double sec1, double sec2) {
		double aoff = Math.abs(sec2 - sec1);
		if (aoff > 30) {
			aoff = Math.abs(aoff - 30);
		}
		return aoff;
	}
	public boolean skipFrame() {
		return skipFrame(lwirReaderParameters);
	}

	public boolean skipFrame(LwirReaderParameters lrp) {
		if (lrp.lwir_channels.length == 0) {
			LOGGER.error("skipFrame(): No LWIR channels are configured");
			return false;
		}
		int chn = lrp.lwir_channels[0];
		String url =  "http://"+lrp.lwir_ips[0]+":"+IMGSRV_PORTS[chn] +SKIP_FRAME_URL;
			Document dom=null;
			try {
				DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				DocumentBuilder db = dbf.newDocumentBuilder();
				dom = db.parse(url);
				if (!dom.getDocumentElement().getNodeName().equals("meta")) {
					LOGGER.error("skipFrame() in " + url+
							": Root element: expected 'me3ta', got \"" + dom.getDocumentElement().getNodeName()+"\"");
						return false;
				}
	   		} catch(MalformedURLException e){
				LOGGER.error("skipFrame() in " + url+ ": " + e.toString());
				return false;
	   		} catch(IOException  e1){
			LOGGER.error("skipFrame() in " + url+ " - camera did not respond: " + e1.toString());
			return false;
	   		}catch(ParserConfigurationException pce) {
			LOGGER.error("skipFrame() in " + url+ " - PCE error: " + pce.toString());
			return false;
	   		}catch(SAXException se) {
			LOGGER.error("skipFrame() in " + url+ " - SAX error: " + se.toString());
			return false;
	   		}
		return true;
	}

	public ImagePlus [] acquire(String dirpath) {
		return acquire(lwirReaderParameters, dirpath);
	}


	public ImagePlus [] acquire_lepton(LwirReaderParameters lrp, String dirpath) {
		if (!condProgramLWIRCamera(lrp)) {
			LOGGER.error("acquire(): failed to program cameras");
			return null;
		}
		if (lrp.lwir_ffc) {
			calibrate(lrp); // seems to work. test calibration duration and if any images are sent during calibration
		}
		int num_frames = lrp.avg_number + lrp.eo_lag + 2 * lrp.max_frame_diff;
		if (num_frames > lrp.num_frames) {
			LOGGER.warn("acquire(): Wanted to acquire "+num_frames+", but hard limit is set to "+lrp.num_frames+
					" TDOD: stop EO compressor after acquisition of the required number of frames") ;
			num_frames = lrp.num_frames;
		}
		ImagePlus [][] imps = readAllMultiple(
				num_frames,
//				lrp.lwir_telemetry,
				false, // final boolean show,
				lrp.eo_scale);
		if (imps == null) {
			LOGGER.error("acquire(): failed to acquire images");
			return null;
		}
		LOGGER.debug("LWIR_ACQUIRE: got "+imps.length+" image sets");
// Verify fresh FFC for all channels
		double [] after_ffc = new double [lrp.lwir_channels.length];
		for (int chn = 0; chn < after_ffc.length; chn++) {
			double uptime = Double.NaN, ffctime = Double.NaN;
			try {
				uptime =  Double.parseDouble(((String) imps[0][chn].getProperty("TLM_UPTIME")));
				ffctime = Double.parseDouble(((String) imps[0][chn].getProperty("TLM_FFC_TIME")));
				after_ffc[chn] = uptime-ffctime;
			} catch (Exception e) {
				LOGGER.warn("acquire(): TLM_UPTIME and/or TLM_FFC_TIME properties are not available for"+imps[0][chn].getTitle());
			}
			LOGGER.warn(String.format("LWIR channel %d time from FFC: %.3f", chn, after_ffc[chn]));
		}

		int [] lags = new int [lrp.lwir_channels.length + lrp.eo_channels.length];
		for (int i = 0; i < lags.length; i++) {
			lags[i] = (i >= lrp.lwir_channels.length) ? lrp.eo_lag : 0;
		}
		ImagePlus [][] imps_sync =  matchSets(
				imps,
				lrp.max_mismatch_ms * 0.001, // 0.001,
				lrp.max_frame_diff,// 3); // double max_mismatch)
				lags);
		if (imps_sync == null) {
			return null;
		}
        int num_keep = imps_sync.length;
        if (!lrp.avg_all && (lrp.avg_number < imps_sync.length)) {
        	num_keep = lrp.avg_number;
        }
		LOGGER.debug("LWIR_ACQUIRE: got "+imps_sync.length+", requested "+lrp.avg_number+", keeping " + num_keep);
		if (num_keep < imps_sync.length) {
			ImagePlus [][] imps_sync0 = imps_sync;
			imps_sync = new ImagePlus[num_keep][];
			for (int i = 0; i < num_keep; i++) {
				imps_sync[i] = imps_sync0[i];
			}
		}
		ImagePlus [] imps_avg = averageMultiFrames(imps_sync);



		if (dirpath != null) {
			File  f_dir = new File(dirpath);
			// get series path from the first channel file title
			String first_name = imps_sync[0][0].getTitle();
			String set_name = first_name.substring(0, first_name.lastIndexOf('_'));
			String set_path = dirpath+Prefs.getFileSeparator()+set_name;
			File  set_dir = new File(set_path);
			set_dir.mkdirs(); // including parent
			LOGGER.warn("Saving image set to: "+set_dir.getAbsolutePath());
			for (ImagePlus imp:imps_avg) {
				String fname = imp.getTitle();
				fname = fname.substring(0, fname.lastIndexOf('_')) + ".tiff"; // remove _average
   				FileSaver fs=new FileSaver(imp);
   				String path=set_path+Prefs.getFileSeparator()+fname;
   				IJ.showStatus("Saving "+path);
   				LOGGER.info("LWIR_ACQUIRE: 'Saving "+path );
   				fs.saveAsTiff(path);
			}
		}

		if (lrp.isShowImages()) {
			if (imps_avg != null) {
				for (ImagePlus imp: imps_avg) {
					imp.show();
				}
			}
		}
		return imps_avg;
	}

	// TODO: Implement LWIR restart

	// Program cameras only if parameters had changed (including those, that do not actually need to be programmed)
	public boolean condProgramLWIRCamera() {
		return condProgramLWIRCamera(lwirReaderParameters);
	}


	public boolean condProgramLWIRCamera(LwirReaderParameters lrp) {
		boolean ok = true;
		if ((last_programmed == null) || !last_programmed.equals(lrp)) {
			ok = programLWIRCamera(lrp);
			if (ok) {
				last_programmed = lrp.clone();
			}
		}
		return ok;
	}

	//actually (unconditionally) program cameras parameters
	public boolean programLWIRCamera() {
		return programLWIRCamera(lwirReaderParameters);
	}
	
	public Document [] collectXmlResponses(final String [] urls) {
		final Document [] documents = new Document[urls.length]; 
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int indx = indxAtomic.getAndIncrement(); indx < urls.length; indx = indxAtomic.getAndIncrement()) {
   						documents[indx] = collectXmlResponse(urls[indx]);
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
		return documents;
	}

	
	public Document collectXmlResponse(String url) {
		Document document = null;
		try {
			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
			DocumentBuilder db = dbf.newDocumentBuilder();
			document = db.parse(url);
		} catch(MalformedURLException e){
			LOGGER.error("collectXmlResponse() MalformedURLException in " + url+ ": " + e.toString());
			return null;
		} catch(IOException  e1){
			LOGGER.error("collectXmlResponse() in " + url+ " - camera did not respond: " + e1.toString());
			return null;
		}catch(ParserConfigurationException pce) {
			LOGGER.error("collectXmlResponse() in " + url+ " - PCE error: " + pce.toString());
			return null;
		}catch(SAXException se) {
			LOGGER.error("collectXmlResponse() in " + url+ " - SAX error: " + se.toString());
			return null;
		}
		return document;
	}
	
	
	
	
	public boolean resetAllFrames(LwirReaderParameters lrp) {
   		String [] all_ips = lrp.getAllIPs() ;
   		String [] urls = new String [all_ips.length];
   		boolean OK = true; 
   		for (int i = 0; i < all_ips.length; i++) {
//   			urls[i] = "http://"+all_ips[i]+"/"+SCRIPT_RESET; // +"?frame=-2";
   			urls[i] = "http://"+all_ips[i]+"/"+SCRIPT_RESET+"?frame=-2";// works even with default -1
   		}
   		Document [] docs = collectXmlResponses(urls);
   		for (int i = 0; i < urls.length; i++) {
   			if (docs[i] == null) {
   				OK=false;
   			} else {
   				if (!docs[i].getDocumentElement().getNodeName().equals("reset_frames")) {
   					LOGGER.error("resetAllFrames() in " + urls[i]+
   							": Root element: expected 'reset_frames', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
   					OK = false;
   				}
   			}
   		}
		return OK;
	}
	
	public boolean skipMasterFrames(String [] ips, int nskip) {
   		boolean OK = true;
   		String [] urls = new String[ips.length];
   		for (int i = 0; i < ips.length; i++) {
   			urls[i] = "http://"+ips[i]+"/"+SCRIPT_WAIT+"?frame="+(-nskip);
   		}
   		Document [] docs = collectXmlResponses(urls);
   		for (int i = 0; i < urls.length; i++) {
   			if (docs[i] == null) {
   				OK=false;
   			} else {
   				if (!docs[i].getDocumentElement().getNodeName().equals("wait_frames")) {
   					LOGGER.error("resetAllFrames() in " + urls[i]+
   							": Root element: expected 'wait_frames', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
   					OK = false;
   				}
   			}
   		}
		return OK;
	}

	public void printXmlDocument(Document doc) {
		if (doc==null) {
			System.out.println ("***document=null***");
		}
	    NodeList allNodes=doc.getDocumentElement().getElementsByTagName("*");
	    for (int i=0;i<allNodes.getLength();i++) {
	        String name= allNodes.item(i).getNodeName();
	        String value="";
	        try {
	        	value=allNodes.item(i).getFirstChild().getNodeValue();
	        } catch(Exception e) {

	        }
	        System.out.println(name+": "+value);
	    }

	}
	
	public double getFutureTimestamnp (String ip, double delay) { // delay in seconds from master
		if (lwirReaderParameters.getDebugLevel() > 0) {
			System.out.println("Increasing delay by 3.0s when there is some printing");
			delay += 3.0;
		}
		String url= String.format("http://%s/%s?sensor_port=0",ip,SCRIPT_CAPTURE);
		if (lwirReaderParameters.getDebugLevel() > 0) System.out.println("getFutureTimestamnp(): "+url);
		Document doc = collectXmlResponse(url);
		if (doc == null) return Double.NaN;
        if (!doc.getDocumentElement().getNodeName().equals("capture_range")) {
            IJ.showMessage("Root element: expected 'capture_range', got'" + doc.getDocumentElement().getNodeName()+"'"); 
            return Double.NaN;
         }
		if (lwirReaderParameters.getDebugLevel() > 0)  printXmlDocument(doc);
		double this_ts = Double.parseDouble(((Node) 
				(((Node) doc.getDocumentElement().getElementsByTagName("this_timestamp").item(0)).getChildNodes().item(0))).getNodeValue());
		if (lwirReaderParameters.getDebugLevel() > 0) System.out.println(String.format("this_ts=%18.6f",this_ts));
		double late_ts = this_ts + delay; // may point to a fraction of frame
		url= String.format("http://%s/%s?sensor_port=0&ts=%f",ip,SCRIPT_CAPTURE, late_ts);
		doc = collectXmlResponse(url);
		if (lwirReaderParameters.getDebugLevel() > 0) {
			System.out.println("getFutureTimestamnp(): "+url);
	        printXmlDocument(doc);
		}
		if (doc == null) return Double.NaN;
        if (!doc.getDocumentElement().getNodeName().equals("capture_range")) {
            IJ.showMessage("Root element: expected 'capture_range', got'" + doc.getDocumentElement().getNodeName()+"'"); 
            return Double.NaN;
        }
		double target_ts = Double.parseDouble(((Node) 
				(((Node) doc.getDocumentElement().getElementsByTagName("timestamp").item(0)).getChildNodes().item(0))).getNodeValue());
		if (lwirReaderParameters.getDebugLevel() > 0) System.out.println(String.format("getFutureTimestamnp(): target_ts=%18.6f",target_ts));
		return target_ts;
	}
	
	public int [] getPortMasks(LwirReaderParameters lrp) {
		boolean [] selected_channels = lrp.getSelectedChannels();
		int [] port_masks= new int [selected_channels.length / 4]; // {15,15,15,15,15}; // get
		for (int i = 0; i < port_masks.length; i++) {
			port_masks[i]=0;
			for (int j = 0; j < 4; j++) {
				if (selected_channels[i * 4 + j]) {
					port_masks[i] += 1 << j;
				}
			}
		}
		return port_masks;
	}
	/**
	 * 
	 * @param lrp
	 * @param delay
	 * @param duration -1 - start only, 0 - stop only, >0 - start+stop (LWIR cameras)
	 * @param duration_eo -1 - start only, 0 - stop only, >0 - start+stop (EO cameras). If stop only - will wait for both masters
	 * @return <0 - error, >0 last frame of the master camera/port in sequence
	 */
	public int captureRange(LwirReaderParameters lrp, double delay, int duration, int duration_eo) {
		boolean sync2eo = true;
		int lwir_master_index = 0;
		int eo_master_index = 4;
		String [] all_ips = lrp.getAllIPs();
		int num_lwir =      lrp.getLwirIPs().length;
		int [] port_masks= getPortMasks(lrp);
		String master_ip = all_ips[sync2eo? eo_master_index : lwir_master_index];
		double target_ts = getFutureTimestamnp (master_ip, delay);
		if (Double.isNaN(target_ts)) {
			return -1;
		}
		String [] urls = new String [all_ips.length];
		for (int i = 0; i < urls.length; i++) {
			int dur = (i>= num_lwir)?duration_eo:duration;
			urls[i]=String.format("http://%s/%s?sensor_port=0&ts=%f&m=%d&d=%d", all_ips[i], SCRIPT_CAPTURE, target_ts, port_masks[i], dur);
			if ((dur == 0) && ((i==lwir_master_index) || (i==eo_master_index))){
				urls[i] += "&wait";
			}
		}
		boolean OK = true;
   		Document [] docs = collectXmlResponses(urls);
   		if (lrp.getDebugLevel() > 0) {
   			for (int i = 0; i < urls.length; i++) {
   				System.out.println(i+": "+urls[i]);
   			}
   			for (Document doc:docs) {
   				printXmlDocument(doc);
   			}
   		}
   		for (int i = 0; i < urls.length; i++) {
   			if (docs[i] == null) {
   				OK=false;
   			} else {
   				if (!docs[i].getDocumentElement().getNodeName().equals("capture_range")) {
   					LOGGER.error("captureRange() in " + urls[i]+
   							": Root element: expected 'capture_range', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
   					OK = false;
   					// is there error?
   					NodeList allNodes=docs[i].getDocumentElement().getElementsByTagName("error");
   					if (allNodes.getLength()>0) {
   	   					String err = ((Node) 
   	   						(((Node) docs[i].getDocumentElement().getElementsByTagName("error").item(0)).getChildNodes().item(0))).getNodeValue();
   	   					LOGGER.error("captureRange() in " + urls[i]+
   	   							": Got error: \"" + err+"\"");
   	   					OK = false;
   					}
   				}
   			}
   		}
   		if (!OK) {
   			return -1;
   		}
		int master_frame = Integer.parseInt(((Node) 
				(((Node) docs[lwir_master_index].getDocumentElement().getElementsByTagName("frame").item(0)).getChildNodes().item(0))).getNodeValue());
		if ((lrp.getDebugLevel() > 0) && ((duration==0) || (duration_eo==0)) ) {
			System.out.println(String.format("master_frame=%d(0x%x)", master_frame,master_frame));
   			for (int i = 0; i < urls.length; i++) {
   				System.out.println(i+": "+urls[i]);
   			}
		}
//		System.out.println("mater_frame = "+master_frame+", target_ts="+target_ts);
		return master_frame;
	}
	
	public boolean startStopCompressor(LwirReaderParameters lrp, boolean start, double delay) {
//		if (lrp.getDebugLevel() > 0) System.out.println("startStopCompressor() start="+start);
		if (lrp.getDebugLevel() > -1) System.out.println("startStopCompressor() start="+start);
		int duration = start? -1: 0;
		int end_frame = captureRange(lrp, delay, duration, duration);
		if (end_frame < 0) return false;
		// end frame -> master_ts -> eo_frame -> increment -> wait eo_frame
		return true;
		
	}
	
	public boolean runFFC (LwirReaderParameters lrp) {
		int [] port_masks= getPortMasks(lrp);
		int ffc_groups = lrp.getFFCGroups();
		String [] lwir_IPs = lrp.getLwirIPs();
		boolean OK = true;
		int [] group_masks = new int [ffc_groups];
		if (ffc_groups == 1) {
			group_masks = new int [] {15};
		} else if (ffc_groups == 2) {
			group_masks = new int [] {5, 10};
		} else if (ffc_groups == 4) {
			group_masks = new int [] {1, 2, 4, 8};
		} else {
			throw new IllegalArgumentException ("Invalid number of FFC groups:"+ffc_groups+". Only 1,2 and 4 are supportded");
		}
		ArrayList<String> url_list = new ArrayList<String>();
		for (int ig = 0; ig < ffc_groups; ig++) {
			url_list.clear();
			int ip0 = -1;
			for (int i = 0; i < lwir_IPs.length; i++) { // start urls for 4 of lwir port0 and eo - ports 0-3 (total 8 URLs)
				int mask = port_masks[i] & group_masks[ig];
				if (mask != 0) {
					int p=0;
					for (;((1 << p) & mask) == 0; p++); // find lowest non-zero bit
					url_list.add(String.format("http://%s/parsedit.php?immediate&sensor_port=%d&%s=1&*%s=%d",lwir_IPs[i],p,REG_FFC_RUN,REG_FFC_RUN,mask));
					if (ip0 < 0) {
						ip0 = i; // first used ip index
					}
				}
			}
			if (!url_list.isEmpty()) {
				String[] urls = url_list.toArray(new String[url_list.size()]);
				// Print all urls
				if (lrp.getDebugLevel() > 0) {
					System.out.println("Performing FFC cycle with:");
					for (int i = 0; i < urls.length; i++) {
						System.out.println(i+": "+urls[i]);
					}
				}
				// run all programming not including TRIG=4
				Document [] docs = collectXmlResponses(urls);
				for (int i = 0; i < urls.length; i++) {
					if (docs[i] == null) {
						OK=false;
					} else {
						if (!docs[i].getDocumentElement().getNodeName().equals("parameters")) {
							LOGGER.error("runFFC() in " + urls[i]+
									": Root element: expected 'parameters', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
							OK = false;
						}
					}
				}
				int skip_frames = 2 * lrp.getFFCFrames() + lrp.getFFCWaitFrames();  
				skipMasterFrames(new String[] {lwir_IPs[ip0]},  skip_frames); // make sure all previous parameters are applied
			}
		}
		return OK;
	}
	
	/**
	 * Reset circbuf frames in each camera channel
	 * @param lrp
	 * @return
	 */
	public boolean resetFrameBuffers(LwirReaderParameters lrp) {
		int [] port_masks= getPortMasks(lrp);
		ArrayList<String> url_list = new ArrayList<String>();
		String [] all_IPs =  lrp.getAllIPs();
		for (int i = 0; i < all_IPs.length; i++) {
			for (int port = 0; port < IMGSRV_PORTS.length; port++) {
				if ((port_masks[i] & (1 << port)) != 0) {
					url_list.add(String.format("http://%s:%d%s", all_IPs[i],IMGSRV_PORTS[port],IMAGE_URL_RESET));
				}
			}
		}
		String[] urls = url_list.toArray(new String[url_list.size()]);
		boolean OK=true;
		if (lrp.getDebugLevel() > 0) {
			System.out.println("Resetting circbuf buffers with:");
			for (int i = 0; i < urls.length; i++) {
				System.out.println(i+": "+urls[i]);
			}
		}
		// run all programming not including TRIG=4
		Document [] docs = collectXmlResponses(urls);
		for (int i = 0; i < urls.length; i++) {
			if (docs[i] == null) {
				OK=false;
			} else {
				if (!docs[i].getDocumentElement().getNodeName().equals("frames")) {
					LOGGER.error("resetFrameBuffers() in " + urls[i]+
							": Root element: expected 'frames', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
					OK = false;
				}
			}
		}
		return OK;
	}
	
	public boolean isCompressorStopped(LwirReaderParameters lrp) {
		int [] port_masks= getPortMasks(lrp);
		String [] all_IPs = lrp.getAllIPs();
		int ip_index = 0;
		for (; (ip_index < port_masks.length) && (port_masks[ip_index] == 0);  ip_index++);
		if (ip_index >= port_masks.length) {
			throw new IllegalArgumentException ("No channels selected in isCompressing()");
		}
		String url = "";
		for (int port = 0; port < IMGSRV_PORTS.length; port++) {
			if ((port_masks[ip_index] & (1 << port)) != 0) {
				url = (String.format("http://%s:%d%s", all_IPs[ip_index],IMGSRV_PORTS[port],IMAGE_URL_STATUS));
				break;
			}
		}
		if (lrp.getDebugLevel() > 0) {
			System.out.println("Performing isCompressing");
			System.out.println(url);
		}
		Document doc = collectXmlResponse(url);
		String compressor_state = ((Node) 
				(((Node) doc.getDocumentElement().getElementsByTagName("compressor_state").item(0)).getChildNodes().item(0))).getNodeValue();
//		compressor_state = compressor_state.replace ("/(^\")|(\"$)/g", ""); // replace leading/trailing "
		compressor_state = compressor_state.replace ("\"", ""); // replace leading/trailing "
		
		return compressor_state.equals("COMPRESSOR_RUN_STOP");
	}
	
	public ImagePlus[] acquire (LwirReaderParameters lrp, String dirpath) {
		String std = "STD_";
		boolean scale = false;
		condProgramLWIRCamera();
		if (lrp.runFFCBeforeMeasurments()) {
			runFFC(lrp);
		}
		if (!isCompressorStopped(lrp)) {
			startStopCompressor(lrp, false, lrp.getAcquireDelay());
			resetFrameBuffers(lrp); //only needed after stopping
		}
		resetFrameBuffers(lrp); // actually only needed after stopping, here just in case
		// start capturing frames into each camera buffer
//		int end_frame = 
		captureRange(lrp, lrp.getAcquireDelay(), lrp.getAvgNumberLwir(),lrp.getAvgNumberEO());

//	public static final String IMAGE_URL_NEXT=  "/torp/next/wait/img/save"; // same image name, same image number/time
//	public static final String IMAGE_URL_NEXT=  "/torp/wait/img/next/save"; // same image name, same image number/time
//		int [] port_masks= getPortMasks(lrp);
//		String [] all_IPs = lrp.getAllIPs();
		final int num_lwir = lrp.getLwirIPs().length;
		final int avg_lwir = lrp.getAvgNumberLwir();
		final int avg_eo =   lrp.getAvgNumberEO();
		int avg_num = (avg_lwir > avg_eo)? avg_lwir : avg_eo; 
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
		final ImagePlus [] imgs = readImageSet(
				lrp,
				(avg_lwir> 0),      // boolean get_lwir,
				(avg_eo >  0),      // boolean get_eo,
				false, // lrp.isShowImages(), // boolean show,
				scale, // boolean scale,
				std, // String  std,
				""+0); // String dbgtxt)
   		final int num_channels = imgs.length;
		for (int nmeas = 1; nmeas < avg_num; nmeas++) {
			boolean use_lwir = nmeas < avg_lwir;
			boolean use_eo =   nmeas < avg_eo;
			final ImagePlus [] imgs1 = readImageSet(
					lrp,
					use_lwir,      // boolean get_lwir,
					use_eo,      // boolean get_eo,
					false, // lrp.isShowImages(), // boolean show,
					scale, // boolean scale,
					std, // String  std,
					""+nmeas); // String dbgtxt)
// add img1 to img in a multi-threaded way			
			indxAtomic.set(0);
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				@Override
					public void run() {
	   					for (int chn = indxAtomic.getAndIncrement(); chn < num_channels; chn = indxAtomic.getAndIncrement()) {
	   						if (imgs1[chn] != null) {
//	   							int width =  imgs[chn].getWidth();
//	   							int height = imgs[chn].getHeight();
	   							//	   						String title = sets[0][chn].getTitle()+"_average"+num_frames;
	   							float [] pixels_avg = (float []) imgs [chn].getProcessor().getPixels();
	   							float [] pixels =     (float []) imgs1[chn].getProcessor().getPixels();
	   							//	   						double scale = 1.0/num_frames;
	   							for (int i = 0; i < pixels_avg.length; i++) {
	   								pixels_avg[i] += pixels[i];
	   							}
//	   							imgs [chn].getProcessor().setPixels(pixels_avg); // or it is not needes as it is already in-place?
	   							/*
	   							ImageProcessor ip=new FloatProcessor(width,height);
	   							imps[chn].setProperty("average", ""+num_frames);
	   							if (motorsPosition!=null) for (int m=0;m<motorsPosition.length;m++ ) {
	   								imps_avg[chn].setProperty("MOTOR"+(m+1), ""+motorsPosition[m]);
	   							}
	   							ImagejJp4Tiff.encodeProperiesToInfo(imps_avg[chn]);
	   							*/
	   						}
	   					}
	   				}
	   			};
	   		}
	   		startAndJoin(threads);
		}
		// now divide average by number of measurements and finalize
		indxAtomic.set(0);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int chn = indxAtomic.getAndIncrement(); chn < num_channels; chn = indxAtomic.getAndIncrement()) {
   						if (imgs[chn] != null) {
   							boolean is_lwir = chn < (num_lwir*IMGSRV_PORTS.length);
//   							int num_frames = (chn < (num_lwir*IMGSRV_PORTS.length))? avg_lwir : avg_eo;
   							int num_frames = is_lwir? avg_lwir : avg_eo;
   							if (num_frames > 0) {
   								float [] pixels_avg = (float []) imgs[chn].getProcessor().getPixels();
   								double scale = 1.0/num_frames;
   								for (int i = 0; i < pixels_avg.length; i++) {
   									pixels_avg[i] *=scale;
   								}
   								imgs[chn].getProcessor().resetMinAndMax();
//   							imgs [chn].getProcessor().setPixels(pixels_avg); // or it is not needes as it is already in-place?
//   							ImageProcessor ip=new FloatProcessor(width,height);
   								imgs[chn].setProperty("average", ""+num_frames);
   								if (motorsPosition!=null) for (int m=0;m<motorsPosition.length;m++ ) {
   									imgs[chn].setProperty("MOTOR"+(m+1), ""+motorsPosition[m]);
   								}
								imgs[chn].setProperty(
										LwirReaderParameters.SENSOR_TYPE,
										LwirReaderParameters.SENSOR_TYPES[is_lwir?1:0]);
   								ImagejJp4Tiff.encodeProperiesToInfo(imgs[chn]);
   							}
   						}
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
   		
		if (dirpath != null) {
			File  f_dir = new File(dirpath);
			// get series path from the first channel file title
			String first_name = imgs[0].getTitle();
			String set_name = first_name.substring(0, first_name.lastIndexOf('_'));
			String set_path = dirpath+Prefs.getFileSeparator()+set_name;
			File  set_dir = new File(set_path);
			set_dir.mkdirs(); // including parent
			LOGGER.warn("Saving image set to: "+set_dir.getAbsolutePath());
			for (ImagePlus imp:imgs) {
				if (imp==null) { // channel was not selected?
					continue; 
				}
				String fname = imp.getTitle();
//				fname = fname.substring(0, fname.lastIndexOf('_')) + ".tiff"; // remove _average
//				fname = fname.substring(0, fname.lastIndexOf('.')) + ".tiff"; // remove _average
				fname = fname + ".tiff"; // remove _average
   				FileSaver fs=new FileSaver(imp);
   				String path=set_path+Prefs.getFileSeparator()+fname;
   				IJ.showStatus("Saving "+path);
   				LOGGER.info("LWIR_ACQUIRE: 'Saving "+path );
   				fs.saveAsTiff(path);
			}
		}
   		
		if (lrp.isShowImages()) {
			for (ImagePlus imp: imgs) {
				imp.show();
			}
		}
		return null;
	}
	
	public ImagePlus[] readImageSet(
			LwirReaderParameters lrp,
			boolean get_lwir,
			boolean get_eo,
			boolean show,
			boolean scale,
			String  std,
			String dbgtxt)
	{
		int [] port_masks= getPortMasks(lrp);
		ArrayList<String> url_list = new ArrayList<String>();
		String [] all_IPs =   lrp.getAllIPs();
		String [] lwir_IPs =  lrp.getLwirIPs();
		for (int i = 0; i < all_IPs.length; i++) {
			if ((i < lwir_IPs.length)? get_lwir:get_eo) {
				for (int port = 0; port < IMGSRV_PORTS.length; port++) if ((port_masks[i] & (1 << port)) != 0) {
					url_list.add(String.format("http://%s:%d/bchn%d%s",all_IPs[i],IMGSRV_PORTS[port],i*IMGSRV_PORTS.length,IMAGE_URL_NEXT));
				}
			}
		}
		String[] urls = url_list.toArray(new String[url_list.size()]);
		ImagePlus [] imps = new ImagePlus[urls.length];
		if (imagejJp4TiffMulti == null) {
			imagejJp4TiffMulti = new ImagejJp4TiffMulti();
		}
		LOGGER.info("---- Acquiring frame set "+dbgtxt);
		try {
			imagejJp4TiffMulti.getMultiImages(urls, imps,  scale,  std);
		} catch (IOException e) {
			LOGGER.error("readAllMultiple0: IOException, " + dbgtxt);
		} catch (FormatException e) {
			LOGGER.error("readAllMultiple0:FormatException, " + dbgtxt);
		}
		double [] img_seconds = new double [urls.length];
		int []    img_numbers = new int [urls.length];
		String []img_names = new String [urls.length];
		for (int i = 0; i < imps.length; i++) {
			String dt = null;
			try {
				dt = (String) imps[i].getProperty("DATE_TIME");
			} catch (Exception e) {
				LOGGER.error("Something is wrong!, no DATE_TIME in the image "+urls[i]);
			}
			if (dt==null) {
				System.out.println("i="+i+", dt is NULL!");
			}
			try {
				img_seconds[i] = Double.parseDouble(dt.substring(dt.lastIndexOf(":")+1));
			} catch (Exception e) {
				System.out.println("i="+i+", dt="+dt);
				System.out.println("i="+i+", dt="+dt);
			}
			
			try {
				img_numbers[i] = Integer.parseInt((String) imps[i].getProperty("STD_Image_Number"));
			} catch (Exception e) {
				img_numbers[i] = -1;
			}
			img_names[i] = (String) imps[i].getProperty("CONTENT_FILENAME");
			if (lrp.getDebugLevel() > 0) {
				LOGGER.warn("Seconds for "+i+" - "+img_seconds[i]+", number"+img_numbers[i]+", name "+img_names[i]);
			}
		}
		// Make a sparse array for all lwir+eo channels with nulls for channels that are not measured. 
		ImagePlus [] imps_all = new ImagePlus [all_IPs.length * IMGSRV_PORTS.length];
		int indx = 0;
		for (int i = 0; i < all_IPs.length; i++) {
			if ((i < lwir_IPs.length)? get_lwir:get_eo) {
				for (int port = 0; port < IMGSRV_PORTS.length; port++) if ((port_masks[i] & (1 << port)) != 0) {
					imps_all[i* IMGSRV_PORTS.length + port] = imps[indx++];
				}
			}
		}
		if (show) {
			for (ImagePlus imp: imps) {
				imp.show();
			}
		}
		return imps_all; // imps;
	}

	
	
	public boolean programLWIRCamera(LwirReaderParameters lrp) {
		return  programLWIRCamera(lrp, 0);
	}
	
	/**
	 * Program all 5 cameras for calibration mode
	 * @param lrp LWIR Reader parameters
	 * @param compressor_run 0 - exit with compressors stopped, 2 - keep compressors running
	 * @return true if OK, false otherwise
	 */
	
	public boolean programLWIRCamera(LwirReaderParameters lrp, int compressor_run) {
		String [] all_IPs =  lrp.getAllIPs();
		String [] lwir_IPs = lrp.getLwirIPs();
		String eo_ip =       lrp.getEOIP();
		String [] two_IPs=   {lwir_IPs[0], eo_ip};
		String [] urls;
//		int eo_master_port = 0;
		skipMasterFrames(two_IPs, 16); // make sure LWIR and EO (slow) flushed frame sequencer (overkill)
		boolean OK = resetAllFrames(lrp);
		skipMasterFrames(two_IPs,  1); 
		/// prepare trigger parameters (common for all)
		urls = new String[all_IPs.length+3]; // for EO - per port individual
		for (int i = 0; i < urls.length; i++) { // start urls for 4 of lwir port0 and eo - ports 0-3 (total 8 URLs)
			int ip = (i < all_IPs.length)?i:(all_IPs.length-1);
			int port = i-ip;
			urls[i] = "http://"+all_IPs[ip]+"/parsedit.php?immediate&sensor_port="+port;

		}
		for (int i = 0; i < all_IPs.length; i++) { // apply to all IPs port0 (LWIR and EO)
			urls[i] += "&TRIG_OUT=0x66555"+
				  	   "&TRIG_CONDITION=0x95555"+
//					   "&TRIG_DELAY="+lrp.lwir_trig_dly+
					   "&TRIG_BITLENGTH=31"+
					   "&EXTERN_TIMESTAMP=1"+
					   "&XMIT_TIMESTAMP=1";

		}
		for (int i = 0; i < lwir_IPs.length; i++) { // Apply to only 4 lwir cameras, port 0 (some to all ports) 
			urls[i] += "&TRIG_DELAY="+lrp.lwir_trig_dly+"&*TRIG_DELAY=15"+ // apply to all ports
					   "&BITS=16&*BITS=15"+
					   "&COLOR="+COLOR_RAW +"&*COLOR=15"+
					   "&WOI_HEIGHT="+(LWIR_HEIGHT + (lrp.lwir_telemetry?LWIR_TELEMETRY_LINES:0))+"&*WOI_HEIGHT=15"+
					   "&"+REG_FFC_FRAMES+"="+lrp.getFFCFrames() +"&*"+REG_FFC_FRAMES+"=15"; // apply to all channels
			urls[i] += "&COMPRESSOR_RUN=2&*COMPRESSOR_RUN=15"; // turn compressor and keep fort at least one frame after TRIG=4
		}
		
		for (int chn:lrp.eo_channels) { // add color settings to all EO channels			
	   		int minExposure=10; // usec
	   		int maxExposure=1000000; //usec
	   		int minGain=(int) (0x10000*1.0);
	   		int maxGain=(int) (0x10000*15.75);
	   		int minScale=0;
	   		int maxScale=(int) (0x10000*4.0);
	   		int exposure= (int) (Math.round(1000*lrp.eo_exposure_ms * lrp.eo_exp_corr[chn]));
	   		int autoExposureMax= (int) (Math.round(1000*lrp.eo_max_autoexp_ms));
	   		int gain=     (int) (Math.round(0x10000*lrp.eo_gain_g));
	   		int rScale=   (int) (Math.round(0x10000*lrp.eo_gain_rg*lrp.eo_gcorr_rbgb[3*chn+0]));
	   		int bScale=   (int) (Math.round(0x10000*lrp.eo_gain_bg*lrp.eo_gcorr_rbgb[3*chn+1]));
	   		int gScale=   (int) (Math.round(0x10000*                 lrp.eo_gcorr_rbgb[3*chn+2]));
	   		int autoExp=  lrp.eo_autoexp?1:0;
	   		int autoWB=   lrp.eo_whitebal?1:0;
	   		if (exposure<minExposure) exposure=minExposure; else if (exposure>maxExposure) exposure=maxExposure;
	   		if (autoExposureMax<minExposure) autoExposureMax=minExposure; else if (autoExposureMax>maxExposure) autoExposureMax=maxExposure;
	   		if (gain<minGain) gain= minGain ; else if (gain> maxGain) gain= maxGain;
	   		if (rScale<minScale) rScale= minScale ; else if (rScale> maxScale) rScale= maxScale;
	   		if (bScale<minScale) bScale= minScale ; else if (bScale> maxScale) bScale= maxScale;
	   		if (gScale<minScale) gScale= minScale ; else if (gScale> maxScale) gScale= maxScale;

	   		urls[all_IPs.length-1+chn] +=
	   				"&COLOR="+COLOR_JP4+ // "*1"+ // JP4 always
	   				"&QUALITY="+lrp.eo_quality+ // "*0"+
	   				"&EXPOS="+exposure+ // "*0"+
		   			"&AUTOEXP_EXP_MAX="+autoExposureMax+//"*0"+
		   			"&AUTOEXP_ON="+autoExp+//"*0"+
		   			"&GAING="+gain+//"*0"+
		   			"&RSCALE="+rScale+//"*0"+
		   			"&BSCALE="+bScale+//"*0"+
		   			"&GSCALE="+gScale+//"*0"+ // GB/G ratio
		   			"&WB_EN="+autoWB+//"*0"+
	   				"&DAEMON_EN_TEMPERATURE=1";//"*0";
	   		if (lrp.eo_full_window) {
	   			urls[all_IPs.length-1+chn] +=
	   					"&WOI_LEFT=0"+
	   					"&WOI_TOP=0"+
	   					"&WOI_WIDTH=2592"+
	   					"&WOI_HEIGHT=1936";
	   			if (chn == 0) {
	   				urls[all_IPs.length-1] += "&COMPRESSOR_RUN=2&*COMPRESSOR_RUN=15"; // turn compressor and keep fort at least one frame after TRIG=4
	   			}
	   		}
		}
		// run all programming not including TRIG=4
   		Document [] docs = collectXmlResponses(urls);
		// Print all urls
		if (lrp.getDebugLevel() > 0) {
			System.out.println("Configuring cameras with:");
			for (int i = 0; i < urls.length; i++) {
				System.out.println(i+": "+urls[i]);
			}
		}
   		for (int i = 0; i < urls.length; i++) {
   			if (docs[i] == null) {
   				OK=false;
   			} else {
   				if (!docs[i].getDocumentElement().getNodeName().equals("parameters")) {
   					LOGGER.error("programLWIRCamera() in " + urls[i]+
   							": Root element: expected 'parameters', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
   					OK = false;
   				}
   			}
   		}
		skipMasterFrames(two_IPs,  3); // make sure all previous parameters are applied
		
        // set external trigger mode for all LWIR and EO cameras
		urls = new String[all_IPs.length];
		for (int i = 0; i < urls.length; i++) {
			urls[i] = "http://"+all_IPs[i]+"/parsedit.php?immediate&sensor_port=0&TRIG=4&*TRIG=15"+
                    "&COMPRESSOR_RUN=" + compressor_run +"*5&*COMPRESSOR_RUN=15"; // delay turning off COMPRESSOR_RUN
		}
		
		docs = collectXmlResponses(urls);
		// Print all urls
		if (lrp.getDebugLevel() > 0) {
			System.out.println("\nSetting trigger mode = 4 on all cameras/ports");
			for (int i = 0; i < urls.length; i++) {
				System.out.println(i+": "+urls[i]);
			}
		}
   		for (int i = 0; i < urls.length; i++) {
   			if (docs[i] == null) {
   				OK=false;
   			} else {
   				if (!docs[i].getDocumentElement().getNodeName().equals("parameters")) {
   					LOGGER.error("resetAllFrames() in " + urls[i]+
   							": Root element: expected 'parameters', got \"" + docs[i].getDocumentElement().getNodeName()+"\"");
   					OK = false;
   				}
   			}
   		}
		skipMasterFrames(two_IPs,  4); // make sure all previous parameters are applied // waits for both LWIR and EO
		// second reset after cameras running synchronously
		OK &= resetAllFrames(lrp);
		skipMasterFrames(two_IPs,  1); 

		if (OK) {
			last_programmed = lrp.clone();
		}
		return OK;
	}
	
	public boolean programLWIRCamera_lepton(LwirReaderParameters lrp) {
		int lwir_master_port = 0;
		int eo_master_port = 0;
		int num_lwir = lrp.lwir_channels.length;
		int num_eo = lrp.eo_channels.length;
		final String [] urls = new String [num_lwir + num_eo];
		for (int chn:lrp.lwir_channels) {
			urls[chn] = "http://"+lrp.getLwirIP()+"/parsedit.php?immediate&sensor_port="+chn+
					"&BITS=16"+
					"&COLOR="+COLOR_RAW+ // +"*0"; // raw mode - delay 0 - breaks compressor
					"&WOI_HEIGHT="+(LWIR_HEIGHT + (lrp.lwir_telemetry?LWIR_TELEMETRY_LINES:0));
			if (chn == lwir_master_port) {
				urls[chn] +="&TRIG=0*0"+
						"&TRIG_DELAY="+lrp.lwir_trig_dly+"*0"+
						"&TRIG_OUT=419157*0"+
						"&TRIG_BITLENGTH=31*0"+
						"&EXTERN_TIMESTAMP=1*0"+
						"&XMIT_TIMESTAMP=1*0";
			}
		}
		for (int chn:lrp.eo_channels) {
	   		int minExposure=10; // usec
	   		int maxExposure=1000000; //usec
	   		int minGain=(int) (0x10000*1.0);
	   		int maxGain=(int) (0x10000*15.75);
	   		int minScale=0;
	   		int maxScale=(int) (0x10000*4.0);
	   		int exposure= (int) (Math.round(1000*lrp.eo_exposure_ms * lrp.eo_exp_corr[chn]));
	   		int autoExposureMax= (int) (Math.round(1000*lrp.eo_max_autoexp_ms));
	   		int gain=     (int) (Math.round(0x10000*lrp.eo_gain_g));
	   		int rScale=   (int) (Math.round(0x10000*lrp.eo_gain_rg*lrp.eo_gcorr_rbgb[3*chn+0]));
	   		int bScale=   (int) (Math.round(0x10000*lrp.eo_gain_bg*lrp.eo_gcorr_rbgb[3*chn+1]));
	   		int gScale=   (int) (Math.round(0x10000*                 lrp.eo_gcorr_rbgb[3*chn+2]));
	   		int autoExp=  lrp.eo_autoexp?1:0;
	   		int autoWB=   lrp.eo_whitebal?1:0;
	   		if (exposure<minExposure) exposure=minExposure; else if (exposure>maxExposure) exposure=maxExposure;
	   		if (autoExposureMax<minExposure) autoExposureMax=minExposure; else if (autoExposureMax>maxExposure) autoExposureMax=maxExposure;
	   		if (gain<minGain) gain= minGain ; else if (gain> maxGain) gain= maxGain;
	   		if (rScale<minScale) rScale= minScale ; else if (rScale> maxScale) rScale= maxScale;
	   		if (bScale<minScale) bScale= minScale ; else if (bScale> maxScale) bScale= maxScale;
	   		if (gScale<minScale) gScale= minScale ; else if (gScale> maxScale) gScale= maxScale;

	   		urls[num_lwir+chn] = "http://"+lrp.eo_ip+"/parsedit.php?immediate&sensor_port="+chn+
	   				"&COLOR="+COLOR_JP4+ // "*1"+ // JP4 always
	   				"&QUALITY="+lrp.eo_quality+ // "*0"+
	   				"&EXPOS="+exposure+ // "*0"+
		   			"&AUTOEXP_EXP_MAX="+autoExposureMax+//"*0"+
		   			"&AUTOEXP_ON="+autoExp+//"*0"+
		   			"&GAING="+gain+//"*0"+
		   			"&RSCALE="+rScale+//"*0"+
		   			"&BSCALE="+bScale+//"*0"+
		   			"&GSCALE="+gScale+//"*0"+ // GB/G ratio
		   			"&WB_EN="+autoWB+//"*0"+
	   				"&DAEMON_EN_TEMPERATURE=1";//"*0";
	   		if (lrp.eo_full_window) {
	   			urls[num_lwir+chn] +=
	   					"&WOI_LEFT=0"+
	   					"&WOI_TOP=0"+
	   					"&WOI_WIDTH=2592"+
	   					"&WOI_HEIGHT=1936";
	   		}


	   		if (chn == eo_master_port) {
	   			urls[num_lwir+chn] += "&TRIG=4*0"+
	   					"&TRIG_CONDITION=611669*0"+ // external input
	   					"&TRIG_BITLENGTH=31*0"+
	   					"&EXTERN_TIMESTAMP=1*0";
	   		}
		}
		for (int i = 0; i < urls.length; i++) {
				LOGGER.debug("programLWIRCamera(): reading url " + urls[i]);
		}
// multithreaded camera access:
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
   		final boolean [] get_success = new boolean [urls.length];
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					for (int indx = indxAtomic.getAndIncrement(); indx < urls.length; indx = indxAtomic.getAndIncrement()) {
   						Document dom=null;
   						try {
   							DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
   							DocumentBuilder db = dbf.newDocumentBuilder();
   							dom = db.parse(urls[indx]);
   							if (!dom.getDocumentElement().getNodeName().equals("parameters")) {
   								LOGGER.error("programLWIRCamera() in " + urls[indx]+
   										": Root element: expected 'parameters', got \"" + dom.getDocumentElement().getNodeName()+"\"");
   	   							get_success[indx] = false;
   	   							continue;
   							}
   				   		} catch(MalformedURLException e){
								LOGGER.error("programLWIRCamera() in " + urls[indx]+ ": " + e.toString());
	   							get_success[indx] = false;
	   							continue;
   				   		} catch(IOException  e1){
							LOGGER.error("programLWIRCamera() in " + urls[indx]+ " - camera did not respond: " + e1.toString());
   							get_success[indx] = false;
   							continue;
   				   		}catch(ParserConfigurationException pce) {
							LOGGER.error("programLWIRCamera() in " + urls[indx]+ " - PCE error: " + pce.toString());
   							get_success[indx] = false;
   							continue;
   				   		}catch(SAXException se) {
							LOGGER.error("programLWIRCamera() in " + urls[indx]+ " - SAX error: " + se.toString());
   							get_success[indx] = false;
   							continue;
   				   		}
   						get_success[indx] = true;
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
// See if there are any errors
   		boolean allOK = true;
   		for (boolean OK:get_success) {
   			allOK &= OK;
   		}
   		if (allOK){
   			lrp.reset_updated();

   			for (int i = 0; i < FRAMES_SKIP; i++) {
   				if (!skipFrame(lrp)) {
					LOGGER.error("programLWIRCamera():Failed to skip frame");
   				}
   			}
   		}
		return allOK;
	}

	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private Thread[] newThreadArray(int maxCPUs) {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try
		{
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}


}
