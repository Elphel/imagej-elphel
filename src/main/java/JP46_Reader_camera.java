/**
 ** -----------------------------------------------------------------------------**
 ** JP46_Reader_camera.java
 **
 ** Reads Elphel Camera JP46 files into ImageJ, un-applying gamma and gains
 **
 ** Copyright (C) 2010 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  JP46_Reader.java is free software: you can redistribute it and/or modify
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

import java.awt.Button;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.Panel;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.lang.reflect.Array;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

import javax.swing.JFileChooser;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GUI;
import ij.gui.GenericDialog;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.plugin.frame.PlugInFrame;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import loci.common.RandomAccessInputStream;
import loci.formats.tiff.IFD;
import loci.formats.tiff.IFDList;
import loci.formats.tiff.TiffParser;



/* This plugin opens images in Elphel JP4/JP46 format (opens as JPEG, reads MakerNote and converts). */
public class JP46_Reader_camera extends PlugInFrame implements ActionListener {

	/**
	 *
	 */
	private static final long serialVersionUID = 390855361964415147L;
	Panel panel1;
	Panel confpanel;
	Frame instance;

	String arg;

	static File dir;

	public String camera_url = "http://192.168.0.236:8081/";
	public String camera_img = "bimg";
	public String camera_img_new = "towp/wait/bimg"; // will always wait for the next image (repetitive acquisitions get new images)
	public String camera_jp46settings = "";
	public boolean IS_SILENT=true;
	public boolean ABSOLUTELY_SILENT=false;
    public boolean demux=true;
	public String imageTitle="cameraImage";
	private int ExifOffset=0x0c;

	public JP46_Reader_camera() {
		super("JP46 Reader Camera");
		if (IJ.versionLessThan("1.39t")) return;
		if (instance!=null) {
			instance.toFront();
			return;
		}
		instance = this;
		addKeyListener(IJ.getInstance());

		panel1 = new Panel();

		panel1.setLayout(new GridLayout(6, 1, 50, 5));

		addButton("Open JP4/JP46...",panel1);
		addButton("Open JP4/JP46 from camera",panel1);
		addButton("Configure...",panel1);
		addButton("Show image properties",panel1);
		addButton("Decode image info to properties",panel1);
		addButton("Split Bayer",panel1);


		add(panel1);

		pack();
		GUI.center(this);
		setVisible(true);
	}
	public JP46_Reader_camera(boolean showGUI) {
		super("JP46 Reader Camera");
		if (IJ.versionLessThan("1.39t")) return;
		if (instance!=null) {
			instance.toFront();
			return;
		}
		instance = this;
		addKeyListener(IJ.getInstance());

		panel1 = new Panel();

		panel1.setLayout(new GridLayout(6, 1, 50, 5));

		addButton("Open JP4/JP46...",panel1);
		addButton("Open JP4/JP46 from camera",panel1);
		addButton("Configure...",panel1);
		addButton("Show image properties",panel1);
		addButton("Decode image info to properties",panel1);
		addButton("Split Bayer",panel1);
		add(panel1);
		pack();
		GUI.center(this);
		setVisible(showGUI);
		if (!showGUI) {
			this.ABSOLUTELY_SILENT = true;
		}
	}

	void addButton(String label, Panel panel) {
		Button b = new Button(label);
		b.addActionListener(this);
		b.addKeyListener(IJ.getInstance());
		panel.add(b);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		String label = e.getActionCommand();

		/* nothing */
		if (label==null) return;

		/* button */
		if (label.equals("Open JP4/JP46...")) {
			read_jp46(arg,true);
		}else if (label.equals("Open JP4/JP46 (no scale)...")) {
			read_jp46(arg,false);
		}else if (label.equals("Configure...")) {
			showConfigDialog(); // open configure dialog
		}else if (label.equals("Open JP4/JP46 from camera")) {
			openURL(camera_url + camera_img_new + camera_jp46settings, arg, true);
		}else if (label.equals("Open JP4/JP46 from camera (no scale)")) {
			openURL(camera_url + camera_img_new + camera_jp46settings, arg, false);
		}else if (label.equals("Show image properties")) {
			ImagePlus imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","No images selected");
				return;
			}
			listImageProperties (imp_sel);
		}
		else if (label.equals("Decode image info to properties")) {
			ImagePlus imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","No images selected");
				return;
			}
			decodeProperiesFromInfo(imp_sel);
			listImageProperties (imp_sel);
		}
		else if (label.equals("Split Bayer")) {
			ImagePlus imp_sel = WindowManager.getCurrentImage();
			if (imp_sel==null){
				IJ.showMessage("Error","No images selected");
				return;
			}
			splitShowBayer(imp_sel);
		}
		IJ.showStatus("DONE");
	}
    public void splitShowBayer(ImagePlus imp){
    	float [] pixels= (float[]) imp.getProcessor().getPixels();
    	int height=imp.getHeight();
    	int width= imp.getWidth();
    	int halfHeight=height/2;
    	int halfWidth=width/2;
        float [][] bayerPixels=new float[4][halfHeight * halfWidth];
        for (int iy=0;iy<halfHeight;iy++) for (int ix=0;ix<halfWidth;ix++){
        	int oIndex=iy*halfWidth+ix;
        	int iIndex=iy*2*width+ix*2;
        	bayerPixels[0][oIndex]=pixels[iIndex];
        	bayerPixels[1][oIndex]=pixels[iIndex+1];
        	bayerPixels[2][oIndex]=pixels[iIndex+width];
        	bayerPixels[3][oIndex]=pixels[iIndex+width+1];
        }
        ImageStack array_stack=new ImageStack(halfWidth,halfHeight);
        for (int i=0;i<4;i++) array_stack.addSlice("chn-"+i, bayerPixels[i]);
        ImagePlus imp_stack = new ImagePlus(imp.getTitle()+"-BAYER", array_stack);
        imp_stack.getProcessor().resetMinAndMax();
        imp_stack.show();
    }

	public void read_jp46(String arg, boolean scale) {
		JFileChooser fc=null;
		//try {fc = new JFileChooser();}

		fc = new JFileChooser();

		//catch (Throwable e) {IJ.error("This plugin requires Java 2 or Swing."); return;}
		fc.setMultiSelectionEnabled(true);

		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);

		int returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return;
		File[] files = fc.getSelectedFiles();
		if (files.length==0) { // getSelectedFiles does not work on some JVMs
			files = new File[1];
			files[0] = fc.getSelectedFile();
		}
		String path = fc.getCurrentDirectory().getPath()+Prefs.getFileSeparator();
		dir = fc.getCurrentDirectory();
		for (int i=0; i<files.length; i++) {
			open(path, files[i].getName(), arg, scale);
		}
	}

	public boolean showConfigDialog() {
		GenericDialog gd = new GenericDialog("Configure");
		gd.addStringField ("Image title:    ", getTitle(), 20);
		gd.addStringField("Camera Address: ", getURL(), 20);
		gd.addStringField("Image: ",       camera_img, 20);
		gd.addStringField("Image (new): ", camera_img_new, 20);
		gd.addStringField("JP46 Parameters: ", camera_jp46settings, 50);
		gd.addCheckbox("Demux composite frame? ", demux);
		gd.addCheckbox("Silent? ", IS_SILENT);
		//      gd.addCheckbox("JP4 (not JP46)? ", IS_JP4);


		confpanel = new Panel();
		gd.addPanel(confpanel);
		addButton("Open JP4/JP46 (no scale)...", confpanel);
		addButton("Open JP4/JP46 from camera (no scale)", confpanel);

		//Vector textfields = gd.getStringFields();
		//((TextField)fields.elementAt(0)).SetWidth = 20;

		gd.showDialog();
		if (gd.wasCanceled()) return false;
		setTitle(gd.getNextString());
		setURL  (gd.getNextString());
		camera_img = gd.getNextString();
		camera_img_new = gd.getNextString();
		camera_jp46settings = gd.getNextString();
		demux=gd.getNextBoolean();
		IS_SILENT=gd.getNextBoolean();
		return true;
	}
	public ImagePlus open(String directory, String fileName, String arg, boolean scale) {
		return open(directory, fileName, arg, scale, null,true);
	}
	public ImagePlus open(String directory, String fileName, String arg, boolean scale, ImagePlus imp_src) {
		return open(directory, fileName, arg, scale, imp_src,true);
	}


	public ImagePlus open(
			String directory,
			String fileName,
			String arg,
			boolean scale,
			ImagePlus imp_src,
			boolean showImage) {
		long[] ElphelMakerNote=null;
		ImagePlus imp = null;
		boolean reuse_imp=false;
		boolean showDemux=showImage && demux;
		if (demux) showImage=false;
		double [] xtraExif=new double[1]; // ExposureTime
		double [] lla = null;
		try {
			imp = openJpegOrGif(directory, fileName);
			if (imp == null) {
				IJ.showMessage("JP46 Reader Error", "Could not open "+directory+"" + fileName + " as JPEG/JP46");
			} else {
				if ((imp_src==null)&& showImage) imp.show(); /* Shows before re-ordering*/
				ElphelMakerNote = readElphelMakerNote(directory, fileName, 16,xtraExif); /* after or 8.2.2 */
				if (ElphelMakerNote==null) ElphelMakerNote = readElphelMakerNote(directory, fileName, 14,xtraExif); /* after or 8.0.8.32 */
				if (ElphelMakerNote==null) ElphelMakerNote = readElphelMakerNote(directory, fileName, 12,xtraExif); /* after or 8.0.7.3 */
				if (ElphelMakerNote==null) ElphelMakerNote = readElphelMakerNote(directory, fileName, 8 ,xtraExif); /* before 8.0.7.3 */
				lla = readGpsLLA(directory, fileName);
				if (lla != null) {
					imp.setProperty("LATITUDE",   String.format("%f",lla[0]));
					imp.setProperty("LONGITUDE",  String.format("%f",lla[1]));
					imp.setProperty("ALTITUDE",   String.format("%f",lla[2]));
				}
			}
		} catch (IOException e) {
			IJ.showStatus("");
			String error = e.getMessage();
			if (error==null || error.equals("")) error = ""+e;
			IJ.showMessage("JP46 Reader", ""+error);
			return null;
		}
		if (imp!=null) {
			reuse_imp=jp46Reorder(imp, ElphelMakerNote, scale, imp_src);
			if (reuse_imp) {
				imp=imp_src;
			} else if ((imp_src!=null)&& showImage) { /* tried to reuse, but wrong size */
				imp.show(); /* never did that before */
			}
			if ((xtraExif!=null) && !Double.isNaN(xtraExif[0])){
				imp.setProperty("EXPOSURE",  String.format("%f",xtraExif[0]));

			}
			if (showImage) imp.updateAndDraw(); /* Redisplays final image*/

			if (showDemux) {
				if (!this.ABSOLUTELY_SILENT) System.out.println("demuxing...");
				ImagePlus imp_0 = demuxImage(imp, 0); if (imp_0!=null) imp_0.show();
				ImagePlus imp_1 = demuxImage(imp, 1); if (imp_1!=null) imp_1.show();
				ImagePlus imp_2 = demuxImage(imp, 2); if (imp_2!=null) imp_2.show();
				if ((imp_0==null) && (imp_0==null) && (imp_0==null)) imp.show(); // Show original image if demux failed (single original)
			}

			return imp;
		}
		return null;
	}

	public ImagePlus openURL(ImagePlus imp_src) {
		if (imp_src==null) return openURL(camera_url + camera_img_new + camera_jp46settings, arg, true);
		return openURL(camera_url + camera_img_new + camera_jp46settings, arg, true, imp_src,true);
	}

	public ImagePlus openURL() {
		return openURL(camera_url + camera_img_new + camera_jp46settings, arg, true);
	}

	public ImagePlus openURL(String url, String arg, boolean scale) {
		return openURL(url, arg, scale, null,true);
	}

	public ImagePlus openURL(
			String url,
			String arg,
			boolean scale,
			ImagePlus imp_src,
			boolean showImage) {
		long[] ElphelMakerNote=null;
		ImagePlus imp = null;
		boolean reuse_imp=false;
		boolean showDemux=showImage && demux;
		if (demux) showImage=false;
		double [] xtraExif=new double[1]; // ExposureTime


//		System.out.println("imp_src is "+((imp_src!=null)?"not ":"")+"null");
		try {
			imp = openJpegOrGifUsingURL(url);
			if (imp == null) {
				IJ.showMessage("JP46 Reader Error", "Could not open the URL: " + url + " as JPEG/JP46");
			} else {
				if ((imp_src==null) && showImage) {
//					System.out.println("show() 1");
					imp.show(); /* Shows before re-ordering*/
				}
				/// get rid of the "/towp/wait" if any - there is a chance to re-read the same image
				ElphelMakerNote = readElphelMakerNoteURL(url.replaceFirst("/towp/wait",""), 16,xtraExif); /* after or 8.2.2 */
				if (ElphelMakerNote==null) ElphelMakerNote = readElphelMakerNoteURL(url.replaceFirst("/towp/wait",""), 14,xtraExif); /* after or 8.0.8.32 */
				if (ElphelMakerNote==null) ElphelMakerNote = readElphelMakerNoteURL(url.replaceFirst("/towp/wait",""), 12,xtraExif); /* after or 8.0.7.3 */
				if (ElphelMakerNote==null) ElphelMakerNote = readElphelMakerNoteURL(url.replaceFirst("/towp/wait",""), 8 ,xtraExif ); /* before 8.0.7.3 */

			}
		} catch (IOException e) {
			IJ.showStatus("");
			String error = e.getMessage();
			if (error==null || error.equals(""))
				error = ""+e;
			IJ.showMessage("JP46 Reader", ""+error);
			return null;
		}
		if (imp!=null) {
			reuse_imp=jp46Reorder(imp, ElphelMakerNote, scale, imp_src);
			if (reuse_imp) {
				imp=imp_src;
			} else if ((imp_src!=null) && showImage) { /* tried to reuse, but wrong size */
//				System.out.println("show() 2");
				imp.show(); /* never did that before */
			}
			if ((xtraExif!=null) && !Double.isNaN(xtraExif[0])){
				imp.setProperty("EXPOSURE",  String.format("%f",xtraExif[0]));
			}
			if (showImage) imp.updateAndDraw(); /* Redisplays final image*/
			if (showDemux) {
				if (!this.ABSOLUTELY_SILENT) System.out.println("demuxing...");
				ImagePlus imp_0 = demuxImage(imp, 0); if (imp_0!=null) imp_0.show();
				ImagePlus imp_1 = demuxImage(imp, 1); if (imp_1!=null) imp_1.show();
				ImagePlus imp_2 = demuxImage(imp, 2); if (imp_2!=null) imp_2.show();
				if ((imp_0==null) && (imp_0==null) && (imp_0==null)) imp.show(); // Show original image if demux failed (single original)
			}

			return imp;
		}
		return null;
	}

	boolean jp46Reorder(ImagePlus imp, long[] MakerNote, boolean scale) {
		return jp46Reorder(imp, MakerNote, scale, null);
	}

	void swapArrayElements (double[]arr,int i, int j) {
		double tmp=arr[i];
		arr[i]=arr[j];
		arr[j]=tmp;
	}
	void swapArrayElements (long[]arr,int i, int j) {
		long tmp=arr[i];
		arr[i]=arr[j];
		arr[j]=tmp;
	}
	boolean  jp46Reorder(ImagePlus imp, long[] MakerNote, boolean scale, ImagePlus imp_src) {
		//    int MARGIN=2; // 2 pixels in JP4/JP46 mode around WOI
		double[] gains= new double[4];
		double[] blacks= new double[4];
		double[] blacks256= new double[4];
		double[] gammas= new double[4];
		long []   gamma_scales= new long[4]; /* now not used, was scale _after_ gamma is applied, 0x400(default) corresponds to 1.0 */
		int i;
		double[][] rgammas=new double[4][];
		double min_gain;
		long WOI_LEFT,WOI_WIDTH,WOI_TOP,WOI_HEIGHT,BAYER_MODE,DCM_HOR,DCM_VERT,BIN_HOR,BIN_VERT;
		long COLOR_MODE=0;
		long FLIPH=0;
		long FLIPV=0;
		long  HEIGHT1=0;
		long  HEIGHT2=0;
		long  HEIGHT3=0;
		long  BLANK1=0;
		long  BLANK2=0;
		boolean FLIPH1=false;
		boolean FLIPV1=false;
		boolean FLIPH2=false;
		boolean FLIPV2=false;
		boolean FLIPH3=false;
		boolean FLIPV3=false;
		boolean COMPOSITE=false;
		boolean PORTRAIT=false;
		boolean YTABLEFORC=false;
		long     QUALITY=0;
		long     CQUALITY=0;
		long     CORING_INDEX_Y=0;
		long     CORING_INDEX_C=0;
		double []   satValue={255.0, 255.0, 255.0, 255.0};
		if  (MakerNote !=null) {
			for (i=0;i<4;i++) { /* r,g,gb,b */
				gains[i]= MakerNote[i]/65536.0;
				blacks[i]=(MakerNote[i+4]>>24)/256.0;
				gammas[i]=((MakerNote[i+4]>>16)&0xff)/100.0;
				gamma_scales[i]=MakerNote[i+4] & 0xffff;
				imp.setProperty("gains_"+i,String.format("%f",gains[i]));
				imp.setProperty("blacks_"+i,String.format("%f",blacks[i]));
				imp.setProperty("gammas_"+i,String.format("%f",gammas[i]));
				imp.setProperty("gamma_scales_"+i,String.format("%d",gamma_scales[i]));
			}
			IJ.showStatus("R=" +(0.001*((int)(1000*gains[0])))+
					" G=" +IJ.d2s(gains[1],3)+
					" Gb="+IJ.d2s(gains[2],3)+
					" B=" +IJ.d2s(gains[3],3)+
					" Gamma[0]="+IJ.d2s(gammas[0],3)+
					" Black[0]="+((int) (256*blacks[0])));
			String info=new String();
			info+="Gain\t"+ IJ.d2s(gains[0],3) + "\t"+ IJ.d2s(gains[1],3) + "\t"+ IJ.d2s(gains[2],3) + "\t"+ IJ.d2s(gains[3],3) + "\n"+
			"Gamma\t"+IJ.d2s(gammas[0],3) + "\t"+IJ.d2s(gammas[1],3) + "\t"+IJ.d2s(gammas[2],3) + "\t"+IJ.d2s(gammas[3],3) + "\n"+
			"Black\t"+IJ.d2s(blacks[0],3) + "\t"+IJ.d2s(blacks[1],3) + "\t"+IJ.d2s(blacks[2],3) + "\t"+IJ.d2s(blacks[3],3) + "\n";
			if (MakerNote.length>=14) {
				COMPOSITE= ((MakerNote[10] & 0xc0000000)!=0);
				if (COMPOSITE) {
					HEIGHT1= MakerNote[11] & 0xffff;
					BLANK1= (MakerNote[11]>>16) & 0xffff;
					HEIGHT2= MakerNote[12] & 0xffff;
					BLANK2= (MakerNote[12]>>16) & 0xffff;
					//           HEIGHT3=(( MakerNote[9]>>16)+2*MARGIN) - HEIGHT1-BLANK1-HEIGHT2-BLANK2;
					HEIGHT3=( MakerNote[9]>>16) - HEIGHT1-BLANK1-HEIGHT2-BLANK2;

					FLIPH1= (((MakerNote[10] >> 24) & 1)!=0); // Same value as FLIP_H
					FLIPV1= (((MakerNote[10] >> 25) & 1)!=0); // Same value as FLIP_V
					FLIPH2= (((MakerNote[10] >> 26) & 1)!=0);
					FLIPV2= (((MakerNote[10] >> 27) & 1)!=0);
					FLIPH3= (((MakerNote[10] >> 28) & 1)!=0);
					FLIPV3= (((MakerNote[10] >> 29) & 1)!=0);
				}
				PORTRAIT=  (((MakerNote[13] >>  7) & 1)!=0);
				YTABLEFORC=(((MakerNote[13] >> 15) & 1)!=0);
				QUALITY=   (MakerNote[13] & 0x7f);
				CQUALITY=  ((MakerNote[13]>>8) & 0x7f);
				if (CQUALITY==0) CQUALITY=QUALITY;
				CORING_INDEX_Y= ((MakerNote[13]>>16) & 0x7f);
				CORING_INDEX_C= ((MakerNote[13]>>24) & 0x7f);
				if (CORING_INDEX_C==0) CORING_INDEX_C=CORING_INDEX_Y;

			}
			if (MakerNote.length>=12) {
				WOI_LEFT=   MakerNote[8]&0xffff;
				WOI_WIDTH=  MakerNote[8]>>16;
				WOI_TOP=    MakerNote[9]&0xffff;
				WOI_HEIGHT= MakerNote[9]>>16;
				FLIPH=      MakerNote[10]      & 1;
				FLIPV=     (MakerNote[10]>> 1) & 1;
				BAYER_MODE=(MakerNote[10]>> 2) & 3;
				COLOR_MODE=(MakerNote[10]>> 4) & 0x0f;
				DCM_HOR=   (MakerNote[10]>> 8) & 0x0f;
				DCM_VERT=  (MakerNote[10]>>12) & 0x0f;
				BIN_HOR=   (MakerNote[10]>>16) & 0x0f;
				BIN_VERT=  (MakerNote[10]>>20) & 0x0f;
				info+="WOI_LEFT\t" +  WOI_LEFT+"\t\t\t\n"+
				"WOI_WIDTH\t"+  WOI_WIDTH+"\t\t\t\n"+
				"WOI_TOP\t"  +  WOI_TOP+"\t\t\t\n"+
				"WOI_HEIGHT\t"+ WOI_HEIGHT+"\t\t\t\n"+
				"FLIP_HOR\t"+   (FLIPH!=0)+"\t\t\t\n"+
				"FLIP_VERT\t"+  (FLIPV!=0)+"\t\t\t\n"+
				"BAYER_MODE\t"+ BAYER_MODE+"\t\t\t\n"+
				"COLOR_MODE\t"+ COLOR_MODE+"\t"+ ((COLOR_MODE==2)?"JP46":((COLOR_MODE==5)?"JP4":((COLOR_MODE==0)?"MONO":"OTHER")))  +"\t\t\n"+
				"DECIM_HOR\t"+  DCM_HOR+"\t\t\t\n"+
				"DECIM_VERT\t"+ DCM_VERT+"\t\t\t\n"+
				"BIN_HOR\t"+    BIN_HOR+"\t\t\t\n"+
				"BIN_VERT\t"+   BIN_VERT+"\t\t\t\n";
				imp.setProperty("WOI_LEFT",  String.format("%d",WOI_LEFT));
				imp.setProperty("WOI_WIDTH", String.format("%d",WOI_WIDTH));
				imp.setProperty("WOI_TOP",   String.format("%d",WOI_TOP));
				imp.setProperty("WOI_HEIGHT",String.format("%d",WOI_HEIGHT));
				imp.setProperty("FLIPH",     String.format("%d",FLIPH));
				imp.setProperty("FLIPV",     String.format("%d",FLIPV));
				imp.setProperty("BAYER_MODE",String.format("%d",BAYER_MODE));
				imp.setProperty("COLOR_MODE",((COLOR_MODE==2)?"JP46":((COLOR_MODE==5)?"JP4":((COLOR_MODE==0)?"MONO":"OTHER"))));
				imp.setProperty("DCM_HOR",   String.format("%d",DCM_HOR));
				imp.setProperty("DCM_VERT",  String.format("%d",DCM_VERT));
				imp.setProperty("BIN_HOR",   String.format("%d",BIN_HOR));
				imp.setProperty("BIN_VERT",  String.format("%d",BIN_VERT));

			}
			if (MakerNote.length>=14) {
				info+="COMPOSITE\t" +  COMPOSITE+"\t\t\t\n";
				info+="ORIENTATION\t" +  (PORTRAIT?"PORTRAIT":"LANDSCAPE" )+"\t\t\t\n";
				info+="JPEG quality\t" +  QUALITY+"\t"+((CQUALITY!=QUALITY)?("("+CQUALITY+")"):"")+"\t"+(YTABLEFORC? "Color use Y table":"")+"\t\n";
				info+="Coring index\t" +  CORING_INDEX_Y+"\t"+((CORING_INDEX_C!=CORING_INDEX_Y)?("("+CORING_INDEX_C+")"):"")+"\t\t\n";
				imp.setProperty("COMPOSITE",String.format("%d",COMPOSITE?1:0));
				imp.setProperty("ORIENTATION",(PORTRAIT?"PORTRAIT":"LANDSCAPE" ));
				imp.setProperty("QUALITY",String.format("%d",QUALITY)); //not full
				imp.setProperty("CORING_INDEX_Y",String.format("%d",CORING_INDEX_Y));
				imp.setProperty("CORING_INDEX_C",String.format("%d",CORING_INDEX_C));
			}
			if (MakerNote.length>=16) {
				long [] iTemps={
						(MakerNote[14]>> 0) & 0xffff,
						(MakerNote[14]>>16) & 0xffff,
						(MakerNote[15]>> 0) & 0xffff,
						(MakerNote[15]>>16) & 0xffff};
				for (i=0;i<iTemps.length;i++) if (iTemps[i]!=0xffff){
					double temperature=(iTemps[i]&0xfff)/16.0;
					if (i==0) info+="SYSTEM TEMPERATURE\t" +  temperature+"\t\t\t\n";
					else  info+="SFE "+i+" TEMPERATURE\t" +  temperature+"\t\t\t\n";
					imp.setProperty("TEMPERATURE_"+i,""+temperature);

				}
			}
			if (COMPOSITE) {
				info+="SUB_FRAMES\t- 1 -\t- 2 -\t- 3 -\t\n";
				info+="HEIGHTS\t"+HEIGHT1+"\t"+HEIGHT2+"\t"+HEIGHT3+"\t\n";
				info+="BLANK_ROWS\t" +BLANK1+"\t"+BLANK2+"\t\t\n";
				info+="FLIPH\t"+FLIPH1+"\t"+FLIPH2+"\t"+FLIPH3+"\t\n";
				info+="FLIPV\t"+FLIPV1+"\t"+FLIPV2+"\t"+FLIPV3+"\t\n";

				imp.setProperty("HEIGHT1",String.format("%d",HEIGHT1));
				imp.setProperty("HEIGHT2",String.format("%d",HEIGHT2));
				imp.setProperty("HEIGHT3",String.format("%d",HEIGHT3));
				imp.setProperty("BLANK_ROWS1",String.format("%d",BLANK1));
				imp.setProperty("BLANK_ROWS2",String.format("%d",BLANK2));
				imp.setProperty("FLIPH1",FLIPH1?"1":"0");
				imp.setProperty("FLIPH2",FLIPH2?"1":"0");
				imp.setProperty("FLIPH3",FLIPH3?"1":"0");
				imp.setProperty("FLIPV1",FLIPV1?"1":"0");
				imp.setProperty("FLIPV2",FLIPV2?"1":"0");
				imp.setProperty("FLIPV3",FLIPV3?"1":"0");
			}
			if (!IS_SILENT) new TextWindow(imp.getTitle()+" info", "Parameter\tRed\tGreen(R)\tGreen(B)\tBlue",info, 400, COMPOSITE?(600):((MakerNote.length>=12)?600:160));

			//tw.setLocation(0,0);
// If there are FLIPH, FLIPV - swap gains, gammas, blacks accordingly. later the images will be also flipped
			if (FLIPV!=0) {
				swapArrayElements (gains,       1, 3);
				swapArrayElements (gains,       0, 2);
				swapArrayElements (blacks,      1, 3);
				swapArrayElements (blacks,      0, 2);
				swapArrayElements (gammas,      1, 3);
				swapArrayElements (gammas,      0, 2);
				swapArrayElements (gamma_scales,1, 3);
				swapArrayElements (gamma_scales,0, 2);
			}
			if (FLIPH!=0) {
				swapArrayElements (gains,       1, 0);
				swapArrayElements (gains,       3, 2);
				swapArrayElements (blacks,      1, 0);
				swapArrayElements (blacks,      3, 2);
				swapArrayElements (gammas,      1, 0);
				swapArrayElements (gammas,      3, 2);
				swapArrayElements (gamma_scales,1, 0);
				swapArrayElements (gamma_scales,3, 2);
			}
			for (i=0;i<4;i++) rgammas[i]=elphel_gamma_calc (gammas[i], blacks[i], gamma_scales[i]);
		} else {
			IJ.showMessage("WARNING", "MakerNote not found");
		}
		/**adjusting gains to have the result picture in the range 0..256 */
		min_gain=2.0*gains[0];
		for (i=0;i<4;i++) {
			if (min_gain > gains[i]*(1.0-blacks[i])) min_gain = gains[i]*(1.0-blacks[i]);
			if (!this.ABSOLUTELY_SILENT) System.out.println("gains["+i+"]="+gains[i]+" min_gain="+min_gain);
		}
		imp.setProperty("GAIN",String.format("%f",min_gain)); // common gain

		for (i=0;i<4;i++) gains[i]/=min_gain;
		for (i=0;i<4;i++) blacks256[i]=256.0*blacks[i];
		for (i=0;i<4;i++) {
			if (!this.ABSOLUTELY_SILENT) System.out.println("scaled gains["+i+"]="+gains[i]);
		}


		for (i=0;i<4;i++) {
			if  (MakerNote !=null) {
				if (scale) satValue[i]=((rgammas[i][255])-blacks256[i])/gains[i];
				else       satValue[i]=((rgammas[i][255])-blacks256[i]);
			}   else       satValue[i]=255.0;
			imp.setProperty("saturation_"+i,String.format("%f",satValue[i]));
			if (!this.ABSOLUTELY_SILENT) System.out.println("scaled gains["+i+"]="+gains[i]+" satValue["+i+"]="+satValue[i]);

		}
// swap satValue to match FLIPH,FLIPV again
		if (FLIPV!=0) {
			swapArrayElements (satValue,       1, 3);
			swapArrayElements (satValue,       0, 2);
		}
		if (FLIPH!=0) {
			swapArrayElements (satValue,       1, 0);
			swapArrayElements (satValue,       3, 2);
		}


		for (i=0;i<4;i++) {
			imp.setProperty("saturation_"+i,String.format("%f",satValue[i]));
//System.out.println("saturation_"+i+"\t"+String.format("%f",satValue[i]));

		}


		ImageProcessor ip = imp.getProcessor();
		//		if (FLIPH!=0) ip.flipHorizontal(); /* To correct Bayer */
		//		if (FLIPV!=0) ip.flipVertical(); /* To correct Bayer */

		int width = ip.getWidth();
		int height = ip.getHeight();
		int yb,y,xb,x,offset,nb,xbyr,ybyr;
		float [] pixels = (float[])ip.getPixels();
		float [][] macroblock=new float[16][16];
		float [] pixels1= null;
		boolean IS_JP4=(COLOR_MODE==5);
		boolean IS_JP46=(COLOR_MODE==2);
		if (IS_JP4) pixels1= pixels.clone(); ///JP4 mode

		boolean use_imp_src= (imp_src!=null) && (imp_src.getWidth()==imp.getWidth()) && (imp_src.getHeight()==imp.getHeight());
		ImageProcessor ip_src= use_imp_src? imp_src.getProcessor():ip;
		if (ip_src==null) {
			ip_src=ip;
			use_imp_src=false;
		}
		for (yb=0;yb<(height>>4); yb++) for (xb=0;xb<(width>>4); xb++) { /* iterating macroblocks */
			if (IS_JP4) {
				for (nb=0; nb<4;nb++) {
					xbyr=nb & 1;
					ybyr=(nb>>1) & 1;
					for (y=0;y<8;y++) {
						offset=((yb<<4)+y)*width+ (nb<<3) +((xb>=(width>>5))?(((xb<<5)-width)+(width<<3)):(xb<<5));
						for (x=0;x<8;x++) {
							macroblock[(y<<1) | ybyr][(x<<1) | xbyr]=pixels1[offset+x];
						}
					}
				}
			} else if (IS_JP46) {
				for (y=0;y<16;y++) {
					offset=((yb<<4)+y)*width+(xb<<4);
					for (x=0;x<16;x++) {
						macroblock[((y<<1)&0xe) | ((y>>3) & 0x1)][((x<<1)&0xe) | ((x>>3) & 0x1)]=pixels[offset+x];
					}
				}
			} else { /// mono and other non-processed
				for (y=0;y<16;y++) {
					offset=((yb<<4)+y)*width+(xb<<4);
					for (x=0;x<16;x++) {
						macroblock[y][x]=pixels[offset+x];
					}
				}
			}
			/* apply gammas here */
			if  (MakerNote !=null) {

				if (scale) {
					for (y=0;y<16;y+=2) for (x=0;x<16;x+=2) {
						i=(int) macroblock[y  ][x  ];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y  ][x  ]= (float) (((rgammas[1][i])-blacks256[1])/gains[1]);

						i=(int) macroblock[y  ][x+1];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y  ][x+1]= (float) (((rgammas[0][i])-blacks256[0])/gains[0]);

						i=(int) macroblock[y+1][x  ];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y+1][x  ]= (float) (((rgammas[3][i])-blacks256[3])/gains[3]);

						i=(int) macroblock[y+1][x+1];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y+1][x+1]= (float) (((rgammas[2][i])-blacks256[2])/gains[2]);
					}
				} else {
					for (y=0;y<16;y+=2) for (x=0;x<16;x+=2) {
						i=(int) macroblock[y  ][x  ];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y  ][x  ]= (float) ((rgammas[1][i])-blacks256[1]);

						i=(int) macroblock[y  ][x+1];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y  ][x+1]= (float) ((rgammas[0][i])-blacks256[0]);

						i=(int) macroblock[y+1][x  ];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y+1][x  ]= (float) ((rgammas[3][i])-blacks256[3]);

						i=(int) macroblock[y+1][x+1];  if (i<0) i=0 ; else if (i>255) i=255;
						macroblock[y+1][x+1]= (float) ((rgammas[2][i])-blacks256[2]);
					}
				}
			}
			if (ip_src==null)  System.out.println("ip_src is null");
			//   else if (ip_src.setf==null)  System.out.println("ip_src.setf is null");
			for (y=0;y<16;y++) {
				offset=(yb<<4)+y;
				for (x=0;x<16;x++) {
					ip_src.setf((xb<<4)+ x, offset, macroblock[y][x]); // here null pointer if image was closed
				}
			}
		}
		if (FLIPH!=0) ip_src.flipHorizontal(); /* To correct Bayer */
		if (FLIPV!=0) ip_src.flipVertical(); /* To correct Bayer */

		/* Is it needed here ? */
		/*    imp.draw();
    imp.show(); **/
		if (use_imp_src) copyProperties (imp, imp_src);
		return use_imp_src;
	}

	/* reverses gamma calculations in the camera
      returns double[] table , in the range 0.0..255.996
	 */
	double [] elphel_gamma_calc (double gamma, double black, long gamma_scale) {
		int i;
		double x, black256 ,k;
		int[] gtable = new int[257];
		double[] rgtable =new double[256];
		int ig;
		black256=black*256.0;
		k=1.0/(256.0-black256);
		if (gamma < 0.13) gamma=0.13;
		if (gamma >10.0)  gamma=10.0;
		for (i=0; i<257; i++) {
			x=k*(i-black256);
			if (x < 0.0 ) x=0.0;
			ig= (int) (0.5+65535.0*Math.pow(x,gamma));
			ig=(ig* (int) gamma_scale)/0x400;
			if (ig > 0xffff) ig=0xffff;
			gtable[i]=ig;
		}
		/* now gtable[] is the same as was used in the camera */
		/* FPGA was using linear interpolation between elements of the gamma table, so now we'll reverse that process */
		// double[] rgtable =new  double[256];
		int indx=0;
		double outValue;
		for (i=0; i<256; i++ ) {
			outValue=128+(i<<8);
			while ((gtable[indx+1]<outValue) && (indx<256)) indx++;
			if (indx>=256) rgtable[i]=65535.0/256;
			else if (gtable[indx+1]==gtable[indx]) rgtable[i]=i;
			else           rgtable[i]=indx+(1.0*(outValue-gtable[indx]))/(gtable[indx+1] - gtable[indx]);
		}
		return rgtable;
	}


	long[] readElphelMakerNote(String directory, String fileName, int len, double [] xtraExif) throws IOException  {
		byte [] sig=  {(byte) 0x92 ,0x7c, /* MakerNote*/
				0x00 ,0x04, /* type (long)*/
				0x00 ,0x00 ,0x00 ,0x08 }; /* number*/
		/* should always read all MakerNote - verify that format did not change (edit here when it does). */
		sig[7]=(byte) (len & 0xff);
		sig[6]=(byte) ((len>>8) & 0xff);
		sig[5]=(byte) ((len>>16) & 0xff);
		sig[4]=(byte) ((len>>24) & 0xff);

		RandomAccessFile in = new RandomAccessFile(directory + fileName, "r");
		byte[] head = new byte[4096]; /* just read the beginning of the file */
		in.readFully(head);
		in.close(); // was no close()! -? "too many open files"
		if ((head[this.ExifOffset]!=0x4d) || (head[this.ExifOffset+1]!=0x4d)) {
			IJ.showMessage("JP46 Reader", "Exif Header not found in " + directory + fileName);
			return null;
		}
		/* search for MakerNote */
		long [] note=getExifData (sig, head, len);
		if (xtraExif!=null){
			if (xtraExif.length>0){ // get exposure time
				byte [] exposureTime={
						(byte) 0x82,(byte) 0x9a,0x00,0x05,
						0x00,0x00,0x00,0x01};
				long [] nomDenom=getExifData (exposureTime, head, 2);
				if (nomDenom==null) xtraExif[0]=Double.NaN;
				else {
					xtraExif[0]=(1.0*nomDenom[0])/nomDenom[1];
				}
			}
		}
		return note;
	}

	double [] readGpsLLA(String directory, String fileName)  throws IOException
	{
		RandomAccessFile in = new RandomAccessFile(directory + fileName, "r");
		byte[] head = new byte[4096]; /* just read the beginning of the file */
		in.readFully(head);
		in.close(); // was no close()! -? "too many open files"
		if ((head[this.ExifOffset]!=0x4d) || (head[this.ExifOffset+1]!=0x4d)) {
			IJ.showMessage("JP46 Reader", "Exif Header not found in " + directory + fileName);
			return null;
		}
		// skip to gps header
		byte [] sig = {0x2, 0x2, 0x0, 0x0};
		int i = this.ExifOffset + 2;
		boolean match=false;
		for (i = this.ExifOffset + 2; i < (head.length - sig.length); i++ ) {
			match=true;
			for (int j=0;j<sig.length;j++)if (head[i+j]!=sig[j]){
				match=false;
				break;
			}
			if (match) break;
		}
		i += sig.length;
		if (i >= (head.length-4)) {
			/* IJ.showMessage("JP46 Reader", "MakerNote tag not found in "+directory + fileName+ ", finished at offset="+i); // re-enable with DEBUG_LEVEL*/
			return null;
		}
//		System.out.println(String.format("offset=0x%04x", i));
		double [] lla = new double[3];
		boolean [] lla_sign = new boolean[3];
		for (int ntag = 0; ntag < 10; ntag++) {
			// read 2 big endian bytesof tag, and 2 - type
			int tag =  ((head[i+0] << 8) & 0xff00) +     ((head[i+1] << 0) & 0x00ff);
			int typ =  ((head[i+2] << 8) & 0xff00) +     ((head[i+3] << 0) & 0x00ff);
			int len =  ((head[i+4] <<24) & 0xff000000) + ((head[i+5] <<16) & 0x00ff0000) + ((head[i+ 6] << 8) & 0x0000ff00) + ((head[i+ 7] << 0) & 0x000000ff);
			int val =  ((head[i+8] <<24) & 0xff000000) + ((head[i+9] <<16) & 0x00ff0000) + ((head[i+10] << 8) & 0x0000ff00) + ((head[i+11] << 0) & 0x000000ff);
//			System.out.println(String.format("ntag=%2d tag = 0x%04x type = 0x%04x len = 0x%08x value = 0x%08x", ntag, tag, typ, len, val));
			i+=12;
			switch (tag) {
			case 1: lla_sign[0] = ((val >> 24) & 0xff) == 0x53; break; // S GPSLatitudeRef
			case 2: int [] lat_arr = readExifBE(head, this.ExifOffset + val , 2 * len); //GPSLatitude
			lla[0] = 1.0*lat_arr[0]/lat_arr[1] + (1.0/60)* lat_arr[2]/lat_arr[3] + (1.0/3600)* lat_arr[4]/lat_arr[5];
			break;
			case 3: lla_sign[1] = ((val >> 24) & 0xff) == 0x57; break; // W GPSLongitudeRef
			case 4: int [] long_arr = readExifBE(head, this.ExifOffset + val , 2 * len); // GPSLongitude
			lla[1] = 1.0*long_arr[0]/long_arr[1] + (1.0/60)* long_arr[2]/long_arr[3] + (1.0/3600)* long_arr[4]/long_arr[5];
			break;
			case 5: lla_sign[2] = ((val >> 24) & 0xff) != 0x0; break; // GPSAltitudeRef
			case 6: int [] alt_arr = readExifBE(head, this.ExifOffset + val , 2 * len); // GPSLongitude
			lla[2] = 1.0*alt_arr[0]/alt_arr[1];
			break;
			}
		}
		for (int j = 0; j < lla.length; j++) {
			if (lla_sign[j]) lla[j] *= -1;
		}
//		System.out.println("Lat = "+lla[0]+" Long = "+lla[1]+ " alt = "+lla[2]+"m");
		return lla;

	}
	int [] readExifBE(byte [] head, int offset, int len) {
		int [] data = new int [len];
		for (int i = 0; i < len; i++) {
			int offs = offset + 4*i;
			data[i] = ((head[offs + 0] << 24) & 0xff000000) +
					  ((head[offs + 1] << 16) & 0x00ff0000) +
					  ((head[offs + 2] <<  8) & 0x0000ff00) +
					  ((head[offs + 3] <<  0) & 0x000000ff);
		}
		return data;
	}

	long[] readElphelMakerNoteURL(String url, int len, double [] xtraExif) throws IOException  {
		URL camURL = null;
		URLConnection urlConn = null;
		byte[] data = new byte[4096];
		//      System.out.println("loading exif from: " + url);

		try {
			camURL  = new URL(url);
			urlConn = camURL.openConnection();
			int contentLength = 4096; /* just read the beginning of the file */ //urlConn.getContentLength();
			//inStream = new InputStreamReader(urlConn.getInputStream());

			int bytesRead = 0;
			int offset = 0;
			InputStream raw = urlConn.getInputStream();
			InputStream in = new BufferedInputStream(raw);
			while (offset < contentLength) {
				bytesRead = in.read(data, offset, data.length - offset);
				if (bytesRead == -1)
					break;
				offset += bytesRead;
			}
			in.close();

		} catch(MalformedURLException e){
			System.out.println("Please check the URL:" + e.toString() );
		} catch(IOException  e1){
			System.out.println("Can't read  from the Internet: "+ e1.toString() );
		}
		byte [] sig=  {(byte) 0x92 ,0x7c, /* MakerNote*/
				0x00 ,0x04, /* type (long)*/
				0x00 ,0x00 ,0x00 ,0x08 }; /* number*/
		/* should always read all MakerNote - verify that format did not change (edit here when it does). */
		sig[7]=(byte) (len & 0xff);
		sig[6]=(byte) ((len>>8) & 0xff);
		sig[5]=(byte) ((len>>16) & 0xff);
		sig[4]=(byte) ((len>>24) & 0xff);

		//in = new RandomAccess(RandomAccessFactory.createRO(camURL), "r");

		byte[] head = new byte[4096]; /* just read the beginning of the file */
		head = data;
		//in.readFully(head);

		if ((head[this.ExifOffset] != 0x4d) || (head[this.ExifOffset+1] != 0x4d)) {
			IJ.showMessage("JP46 Reader", "Exif Header not found in " + url);
			return null;
		}
		/* search for MakerNote */
//		return getExifData (sig, head, len);
		long [] note=getExifData (sig, head, len);
		if (xtraExif!=null){
			if (xtraExif.length>0){ // get exposure time
				byte [] exposureTime={
						(byte) 0x82,(byte) 0x9a,0x00,0x05,
						0x00,0x00,0x00,0x01};
				long [] nomDenom=getExifData (exposureTime, head, 2);
				if (nomDenom==null) xtraExif[0]=Double.NaN;
				else {
					xtraExif[0]=(1.0*nomDenom[0])/nomDenom[1];
				}
			}
		}
		return note;

	}



	long [] getExifData (byte [] sig, byte [] head, int len){
		/* search for sig array */
		int i = this.ExifOffset + 2;
		boolean match=false;
		for (i = this.ExifOffset + 2; i < (head.length - sig.length); i++ ) {
			match=true;
			for (int j=0;j<sig.length;j++)if (head[i+j]!=sig[j]){
				match=false;
				break;
			}
			if (match) break;
		}
		i += sig.length;
		if (i >= (head.length-4)) {
			/* IJ.showMessage("JP46 Reader", "MakerNote tag not found in "+directory + fileName+ ", finished at offset="+i); // re-enable with DEBUG_LEVEL*/
			return null;
		}

		int offs=this.ExifOffset+(((head[i]<<24) & 0xff000000) |((head[i+1]<<16) & 0xff0000)| ((head[i+2]<<8) & 0xff00) | (head[i+3] & 0xff));

		// IJ.showMessage("JP46 Reader Debug", "MakerNote starts at offset "+offs);
		if (offs > (head.length-len) ) {
			IJ.showMessage("JP46 Reader", "Error: data (i.e.MakerNote)  starts too far - at offset "+offs+", while we read only "+ head.length+ "bytes");
			return null;
		}

		long[] note=new long[len];
		for (i=0; i<len; i++) note[i]=((head[offs+(i<<2)]&0xff) << 24) | ((head[offs+(i<<2)+1]&0xff) << 16)  | ((head[offs+(i<<2)+2]&0xff) << 8)  | (head[offs+(i<<2)+3]&0xff);
		return note;
	}

	/* Modified from Opener.java */
	ImagePlus openJpegOrGif(String dir, String name) {
		ImagePlus imp = null;
		boolean isTiff = false;
		Image img = Toolkit.getDefaultToolkit().createImage(dir+name);
		if (img!=null) {
			try {
				imp = new ImagePlus(name, img);
			} catch (IllegalStateException e) {
				//java.lang.IllegalStateException: Error loading image
				// Try TIFF
				Opener opener=new Opener(); // Reads Tiff images
				imp=opener.openImage("", dir+name);
				if (imp == null) return null; // error loading image
				isTiff = true;
			}

			if (imp.getType()==ImagePlus.COLOR_RGB) {
				checkGrayJpegTo32Bits(imp);
			}

			IJ.showStatus("Converting to 32-bits");
			new ImageConverter(imp).convertToGray32();
			FileInfo fi = imp.getFileInfo();
			fi.fileName = name;
			fi.directory = dir;
			fi.fileFormat = isTiff?FileInfo.TIFF: FileInfo.GIF_OR_JPG; // even if set originally, it is lost after convertToGray32
			imp.setFileInfo(fi);

			FileInfo ofi = imp.getOriginalFileInfo();
			fi = imp.getFileInfo();

// testing

			if ((ofi!=null) && (ofi.directory!=null) &&  (ofi.fileFormat ==FileInfo.TIFF)) {
				String path = ofi.directory + ofi.fileName;
				EyesisTiff ET = new EyesisTiff();
				ImagePlus imptiff = ET.readTiff(path);
				if (imptiff!=null) {
					imptiff.show();
				}

//				IJ.error("TIFF Dumper", "File path not available or not TIFF file");
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
		}


		return imp;
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
		for (int i=0;i<bstring.length;i++) {
			try {
				bstring[i]= (byte) Integer.parseInt(Array.get(value, i).toString());
			} catch (NumberFormatException e) {
				bstring[i] = 0;
			}
		}
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




	@Override
	public void setTitle (String title) {
		imageTitle=title;
	}
	@Override
	public String getTitle () {
		return imageTitle;
	}
	public void setURL (String url) {
		camera_url=url;
	}
	public String getURL () {
		return camera_url;
	}

	ImagePlus openJpegOrGifUsingURL (String cameraurl) {
		URL url = null;
		ImagePlus imp = null;
		Image img = null;

		/* Validate URL */
		try {
			url = new URL(cameraurl);
		} catch (MalformedURLException e) {
			System.out.println("Bad URL: " + cameraurl);
			return null;
		}

		img = Toolkit.getDefaultToolkit().createImage(url);
		if (!this.ABSOLUTELY_SILENT) System.out.println("loading image from: " + url);
		//      imp = new ImagePlus("test", img);
		imp = new ImagePlus(imageTitle, img);


		if (imp.getType() == ImagePlus.COLOR_RGB) {
			checkGrayJpegTo32Bits(imp);
		}

		IJ.showStatus("Converting to 32-bits");
		new ImageConverter(imp).convertToGray32();

		FileInfo fi = new FileInfo();
		fi.fileFormat = FileInfo.GIF_OR_JPG;
		fi.fileName = "aquired from camera";
		fi.directory = cameraurl;
		imp.setFileInfo(fi);

		return imp;
	}

	public static void checkGrayJpegTo32Bits(ImagePlus imp) {
		ImageProcessor ip = imp.getProcessor();
		int width = ip.getWidth();
		int height = ip.getHeight();
		int[] pixels = (int[])ip.getPixels();
		int c,r,g,b,offset;
		for (int y=0; y<(height-8); y++) {
			offset = y*width;
			for (int x=0; x<(width-8); x++) {
				c = pixels[offset+x];
				r = (c&0xff0000)>>16;
			g = (c&0xff00)>>8;
		b = c&0xff;
		if (!((r==g)&&(g==b))) {
			IJ.error("Warning: color image");
			return;
		}
			}
		}
		IJ.showStatus("Converting to 32-bits");
		new ImageConverter(imp).convertToGray32();
	}
	/* =====Other methods =================================================================== */
	/* ======================================================================== */
	public void listImageProperties (ImagePlus imp) {
		listImageProperties (imp,false);
	}
	public void listImageProperties (ImagePlus imp, boolean toConsole) {
		StringBuffer sb = new StringBuffer();
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
				sb.append(str+"\t"+jp4_prop.getProperty(str)+"\n");
				//				System.out.println(str + "=\t" +jp4_prop.getProperty(str));
			}
		}
		if (toConsole){
			System.out.println(imp.getTitle()+" properties\n"+sb.toString());
		} else {
		  new TextWindow(imp.getTitle()+" properties", "name\tvalue", sb.toString(),400,800);
		}
	}

	/* ======================================================================== */
	public double fracOverExposed(double [] map,   // map of overexposed pixels 0.0 - 0K, >0 (==1.0) - overexposed
			                      int  mapWidth,   // width of the map
			                      int        x0,   // X of the top left corner of the selection
			                      int        y0,   // Y of the top left corner of the selection
			                      int     width,   // selection width
			                      int     height){ // selection height
		int index,i,j;
		int y1=y0+height;
		int over=0;
		for (i=y0;i<y1;i++) {
			index=i*mapWidth+x0;
			for (j=0;j<width;j++) if (map[index++]>0.0) over++;

		}
		return (1.0*over)/width/height;
	}
	// returns 1.0 if there is overexposed pixel, 0.0 - if OK
	/* ======================================================================== */
	public double [] overexposedMap (ImagePlus imp) {
		return overexposedMap (imp, 0.999);
	}
	public double [] overexposedMap (ImagePlus imp, double relativeThreshold) {
		double [] satValues=new double[4];
		boolean noProperties=false;
		int i,j,index;
		for (i=0;i<4;i++) {
	//protect from Double.valueOf(null), move to a function
			if (imp.getProperty("saturation_"+i)!=null) satValues[i]= Double.valueOf((String) imp.getProperty("saturation_"+i)).doubleValue();
			else {
				noProperties=true;
				break;
			}
		}
		if (noProperties) return null;
	//0 - red, 1,2 - green (use Math.min()), 3 - blue
		for (i=0;i<4;i++) satValues[i]*=relativeThreshold;

		ImageProcessor ip=imp.getProcessor();
		int width=imp.getWidth();
		float []pixels=(float[]) ip.getPixels();
		double [] overexposed= new double [pixels.length];
	    for (index=0;index<overexposed.length;index++){
	    	i=(index / width) % 2;
	    	j=1-((index % width) % 2);
	    	overexposed[index]=(pixels[index]>=satValues[i*2+j])?1.0:0.0;
	    }
		return overexposed;
	}

	public ImagePlus demuxImageOrClone(ImagePlus imp, int numImg) {
		ImagePlus imp_new=demuxImage(imp, numImg);
		if (imp_new!=null) return imp_new;
		return demuxClone(imp);
	}

	public ImagePlus demuxClone(ImagePlus imp) {
		ImageProcessor ip=imp.getProcessor().duplicate();
		ImagePlus imp_new=new ImagePlus(imp.getTitle()+"-dup",ip);
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
				imp_new.setProperty(str,jp4_prop.getProperty(str));
		  }
		}
		return imp_new;
	}


	public ImagePlus demuxImage(ImagePlus imp, int numImg) {
		int width= imp.getWidth();
//		int height=imp.getHeight();
		int FLIPGV,FLIPGH;
		int [] FLIPV= new int[3];
		int [] FLIPH= new int[3];
		int [] HEIGHTS=new int[3];
		int [] BLANKS= new int[2];
		Object timestamp=null;
		if (imp.getProperty("FLIPV")!=null) FLIPGV= Integer.valueOf((String) imp.getProperty("FLIPV")).intValue(); else return null;
		if (imp.getProperty("FLIPH")!=null) FLIPGH= Integer.valueOf((String) imp.getProperty("FLIPH")).intValue(); else return null;
		int i;
		for (i=1;i<=3;i++) {
			if (imp.getProperty("FLIPV"+i)!=null)  FLIPV[i-1]=   Integer.valueOf((String) imp.getProperty("FLIPV"+i)).intValue(); else return null;
			if (imp.getProperty("FLIPH"+i)!=null)  FLIPH[i-1]=   Integer.valueOf((String) imp.getProperty("FLIPH"+i)).intValue(); else return null;
			if (imp.getProperty("HEIGHT"+i)!=null) HEIGHTS[i-1]= Integer.valueOf((String) imp.getProperty("HEIGHT"+i)).intValue(); else return null;
		}
		for (i=1;i<=2;i++) {
			if (imp.getProperty("BLANK_ROWS"+i)!=null) BLANKS[i-1]= Integer.valueOf((String) imp.getProperty("BLANK_ROWS"+i)).intValue(); else return null;
		}
		timestamp=imp.getProperty("timestamp");
		if (timestamp!=null);

/*
		System.out.println("FLIPV="+FLIPGV+" FLIPH="+FLIPGH);
		for (i=0;i<3;i++) System.out.println("FLIPV["+i+"]=  "+FLIPV[i]+" FLIPH["+i+"]=  "+FLIPH[i]);
		for (i=0;i<3;i++) System.out.println("HEIGHTS["+i+"]="+HEIGHTS[i]);
		for (i=0;i<2;i++) System.out.println("BLANKS["+i+"]= "+BLANKS[i]);
*/
		Rectangle [] r = new Rectangle[3];
		r[0]=new Rectangle(0, 0,                                        width,HEIGHTS[0]);
		r[1]=new Rectangle(0, HEIGHTS[0]+BLANKS[0],                     width,HEIGHTS[1]);
		r[2]=new Rectangle(0, HEIGHTS[0]+BLANKS[0]+HEIGHTS[1]+BLANKS[1],width,HEIGHTS[2]);
// assuming that 		(HEIGHTS[1]==0) && (HEIGHTS[2]!=0) == false
		if (FLIPGV!=0){
			if (HEIGHTS[1]>0) {
				if (HEIGHTS[2]>0) {
					Rectangle r_tmp=r[0];
					r[0]=r[2];
					r[2]=r_tmp;
				} else {
					Rectangle r_tmp=r[0];
					r[0]=r[1];
					r[1]=r_tmp;
				}
			}
		}
		if (FLIPGV>0) for (i=0;i<3;i++) FLIPV[i]=1-FLIPV[i];
		if (FLIPGH>0) for (i=0;i<3;i++) FLIPH[i]=1-FLIPH[i];

//		for (i=0;i<3;i++) System.out.println("Final: FLIPV["+i+"]=  "+FLIPV[i]+" FLIPH["+i+"]=  "+FLIPH[i]);

// if needed, we'll cut one pixel line. later can modify to add one extra, but then we need to duplicate the pre-last one (same Bayer),
// not just add zeros - later before sliding FHT the two border lines are repeated for 16 times to reduce border effects.
		for (i=0;i<3;i++) {
//			System.out.println("before r["+i+"].x=  "+r[i].x+" r["+i+"].width=  "+r[i].width);
//			System.out.println("before r["+i+"].y=  "+r[i].y+" r["+i+"].height=  "+r[i].height);
			if (((r[i].height & 1)==0 ) & (((r[i].y+FLIPV[i])&1)!=0)) r[i].height-=2;
			r[i].height &=~1;
			if (((r[i].y+FLIPV[i])&1)!=0) r[i].y+=1;

			if (((r[i].width & 1)==0 ) & (((r[i].x+FLIPH[i])&1)!=0)) r[i].width-=2;
			r[i].width &=~1;
			if (((r[i].x+FLIPH[i])&1)!=0) r[i].x+=1;



//			System.out.println("after r["+i+"].x=  "+r[i].x+" r["+i+"].width=  "+r[i].width);
//			System.out.println("after r["+i+"].y=  "+r[i].y+" r["+i+"].height=  "+r[i].height);
		}
		if (r[numImg].height<=0) return null;
//		ImageProcessor ip=imp.getProcessor();
//		ip.setRoi(r[numImg]);
//		ImageProcessor ip_individual=ip.crop().duplicate(); //java.lang.NegativeArraySizeException
/*
 * When using in multithreaded with (probably) the same composite image
Exception in thread "Thread-3564" java.lang.ArrayIndexOutOfBoundsException: 8970912
        at ij.process.FloatProcessor.crop(FloatProcessor.java:706)
        at JP46_Reader_camera.demuxImage(JP46_Reader_camera.java:1104)
        at CalibrationHardwareInterface$CamerasInterface$4.run(CalibrationHardwareInterface.java:1101)

 */
//		ImageProcessor ip=imp.getProcessor();
//		ip.setRoi(r[numImg]);
		ImageProcessor ip_individual=imp.getProcessor().duplicate(); //java.lang.NegativeArraySizeException
		ip_individual.setRoi(r[numImg]);
		ip_individual=ip_individual.crop();

		if (FLIPH[numImg]!=0) ip_individual.flipHorizontal();
		if (FLIPV[numImg]!=0) ip_individual.flipVertical();
		ImagePlus imp_result=new ImagePlus(imp.getTitle()+"-"+numImg,ip_individual);

		//copy all defined properties of the composite image
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
				imp_result.setProperty(str,jp4_prop.getProperty(str));
		  }
		}
// Replaced by copying all properties from the composite image
/*
		for (i=0;i<4;i++) {
			//protect from Double.valueOf(null), move to a function
			if (imp.getProperty("saturation_"+i)!=null) {
				imp_result.setProperty("saturation_"+i, imp.getProperty("saturation_"+i));
			}
		}
		if (timestamp!=null)imp_result.setProperty("timestamp", timestamp);
*/
		// fill in meta data
		return imp_result;
	}
	public void copyProperties (ImagePlus imp_src,ImagePlus imp_dst){
		// copy all the properties to the new image
   		Set<Object> set;
		Properties prop;
		Iterator<Object> itr;
		String str;
		prop=imp_src.getProperties();
		if (prop!=null) {
			set=prop.keySet();
			itr=set.iterator();
			while(itr.hasNext()) {
				str = (String) itr.next();
				imp_dst.setProperty(str,prop.getProperty(str));
			}
		}
	}


	public ImagePlus encodeProperiesToInfo(ImagePlus imp){
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
	public boolean decodeProperiesFromInfo(ImagePlus imp){
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
            String value=allNodes.item(i).getFirstChild().getNodeValue();
    		imp.setProperty(name, value);

	    }

		return true;
	}
	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 * Grabbed from https://github.com/imagej/minimal-ij1-plugin
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Aberration_Calibration.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		// start ImageJ
		new ImageJ();
		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}

}




