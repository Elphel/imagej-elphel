/**
 **
 ** CalibrationHardwareInterface.jave - hardware-related part (cameras, focusing motors,
 ** goniometer motors, lasers 
 **
 ** Copyright (C) 2010-2011 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  CalibrationHardwareInterface.java is free software: you can redistribute it and/or modify
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
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.StringWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.xml.sax.SAXException;

import Jama.LUDecomposition;
import Jama.Matrix;

public class CalibrationHardwareInterface {
	public static class CamerasInterface{
//		JP46_Reader_camera JP4_INSTANCE= new JP46_Reader_camera(false);
		public LaserPointers laserPointers=null;
		private int masterSubCamera=0; // "master" camera index of IP in the list 
		private int masterPort=     0; // "master" camera port (0..3) to apply trigger to 
		private JP46_Reader_camera [] jp4_Instances=null;
		private String [] resetURLs=null;
		private String [] imageURLs=null;
		private String [] metaURLs=null;
		private String    triggerURL="";
		private boolean [] flipImages=null;
		private ImagePlus [] images=   null; // to reuse same instances
		private ImagePlus [] imagesIP= null;
// TODO: when saving/restoring save cameraSubnet, iBaseIP, cameraIPs, so any IPs are OK through config, generate - sequential
		private String cameraSubnet="192.168.0.";
		private int iBaseIP=236;
        private String [] cameraIPs = null;
        private int imgsrvPort=8081;
        private String resetURLcmd="towp/save/pointers"; // advance buffer, next time will wait for the next frame acquired
// will return XML, just "trig" - 1x1 GIF
        private String triggerURLcmd="trig/pointers"; // TRIG=4 should be set in advance, this command will set in single-shot trigger mode
        private String imageURLcmd="torp/wait/bimg"; // will wait if needed. If repeated (as when reading Exif)- won't wait
        private String metaURLcmd="torp/wait/meta";  // will get XML, including timestamp
        private String lastTimestamp="";
		public int debugLevel=2;
		private double lastTemperature=Double.NaN;
		private int     colorMode=5; // JP4
		private boolean noWait=     true; // when false, IRQ_SMART=3 and the frame is available only 1 frame later, when true IRQ_SMART=6, frame is available after compression end
		private boolean nc393 =     false;
		private int     debugSensorNumber=-1; // increase debug level for this particular sensor
		private int     JPEGquality=99;  // JPEG quality
		private boolean cameraAutoExposure =    false;
		private boolean cameraAutoWhiteBalance =false;
		private String  cameraExtraURLCommon =  ""; // (should start with "&")
		private boolean setupTriggerMode=    false;
		private boolean externalTriggerCabling=false;
		private boolean noCabling=           false; // for single camera
		private double  cameraExposure=         5.0;
		private double  scaleExposureForLasers= 0.4;
		private double  scaleExposureForHeadLasers= 0.05;
		private double  cameraAutoExposureMax=  30.0; //Maximal autoexposure value
		private double  cameraGain =  2.0;
		private double  cameraRScale =1.25;
		private double  cameraBScale =1.5;
		private double  cameraGScale =1.0;
		private double [] cameraExposureCorr=null; // per-camera exposure correction1
		private double [] cameraGainCorr   = null;     // per-camera gain correction;
		private double [] cameraRScaleCorr = null;   // per-camera R/G scale correction;
		private double [] cameraBScaleCorr = null;   // per-camera B/G scale correction;
		private double [] cameraGScaleCorr = null;   // per-camera GB/G scale correction;
		private String [] cameraExtraURL =   null;   // per-camera extra URL (should start with "&")
		// these are initialized after being null, when the cameras are probed
		private int    []  cameraFrameNumber=null;
		private boolean [] triggeredMode=    null;   // true - triggered, false - free running
		private boolean [][] sensorPresent=  null;   // probe which sensors (of 3) are detected per system board
		// TODO - try if skipping setting TRIG_PERIOD=0, IRQ_SMART=6 (when they are already set) will fix hanging
		private int [] triggerPeriod=        null;
		private int [] irqSmart=             null;
		private int [] motorsPosition=      null; // motors steps when the images were acquired (for null)
		private long startTime=System.nanoTime();
		private long lastTime=startTime;
		private long thisTime=startTime;
		public boolean reportTiming=false;
		public int cameraBootTimeSeconds=100;
		public int connectionTimeoutMilliseconds=3000;
	
    	private void printTiming(String title){
    		if (this.reportTiming) {
    			this.thisTime=System.nanoTime();
    			System.out.println(title+ " done at "+IJ.d2s(0.000000001*(this.thisTime-this.startTime),3)+
    					" (+"+IJ.d2s(0.000000001*(this.thisTime-this.lastTime),3)+") sec");
    			this.lastTime=this.thisTime;
    		}
    	}
    	private void printTimingInit(){
    		this.startTime=System.nanoTime();
    		this.lastTime=this.startTime;
    		this.thisTime=this.startTime;
    	}

		private int [][] channelMap_23_393={ // ip index, channel number, port
				{0,0,3},{0,0,2},{0,0,0},{0,0,1},{0,0,3},{1,0,2},{1,0,0},{1,0,1},
				{0,1,3},{0,1,2},{0,1,0},{0,1,1},{0,1,3},{1,1,2},{1,1,0},{1,1,1},
				{0,2,3},{0,2,2},{0,2,0},{0,2,1},{0,2,3},{1,2,2},{1,2,0},{1,2,1},
				{2,0,2},{2,0,3}};
    	
    	
		private int [][] channelMap21={ // ip index, channel number
				{0,1,0},{0,0,0},{0,2,0},
				{1,1,0},{1,0,0},{1,2,0},
				{2,1,0},{2,0,0},{2,2,0},
				{3,1,0},{3,0,0},{3,2,0},
				{4,1,0},{4,0,0},{4,2,0},
				{5,1,0},{5,0,0},{5,2,0},
				{6,1,0},{6,0,0},{6,2,0}};
		private int [][] channelMap1={ // ip index, channel number
//				{0,-1}}; // negative channel - single camera
		{0,0,0}}; // Try with 0
		private int [][] channelMap2={ // ip index, channel number
		{0,0,0},{0,1,0}};
		private int [][] channelMap3={ // ip index, channel number
//				{0,-1}}; // negative channel - single camera
		{0,0,0},{1,0,0},{2,0,0}};
		private int [][] channelMap=null;
		public int maxNumberOfThreads=100;
		/**
		 * Initialize JP46_Reader_camera instances, one per sub-camera
		 */
		
        public CamerasInterface(int size, LaserPointers laserPointers){
        	this.laserPointers=laserPointers;
        	initDefaultMap(size);
        	initIPs();
        	initJP4();
        	initCamParsDefaultArrays(this.cameraIPs.length);
        }
        public CamerasInterface(int size){
        	initDefaultMap(size);
        	initIPs();
        	initJP4();
        	initCamParsDefaultArrays(this.cameraIPs.length);
        }
        public void setNumberOfThreads(int n){
        	this.maxNumberOfThreads=n;
        }
        public int getSubCamera (int channelNumber){
        	return ((channelNumber>=0)&& (channelNumber<this.channelMap.length))?this.channelMap[channelNumber][0]:-1;
        }
        public int getSubChannel (int channelNumber){
        	return ((channelNumber>=0)&& (channelNumber<this.channelMap.length))?this.channelMap[channelNumber][1]:-1;
        }

        public int getSensorPort (int channelNumber){
        	return ((channelNumber>=0)&& (channelNumber<this.channelMap.length))?this.channelMap[channelNumber][2]:-1;
        }
        
        // not used anywhere
        public int getChannel (int subCam, int subChn){
        	return getChannel (subCam, subChn, 0); // for compatibility with 353
        }

        public int getChannel (int subCam, int subChn, int port){
        	for (int channelNumber=0;channelNumber<this.channelMap.length;channelNumber++)
        		if ((this.channelMap[channelNumber][0]==subCam) &&
        				(this.channelMap[channelNumber][1]==subChn) &&
        				(this.channelMap[channelNumber][2]==port)) return  channelNumber;
        	return -1;
        }
        
        
		private void initJP4(){
			this.jp4_Instances=new JP46_Reader_camera[this.cameraIPs.length];
			this.resetURLs=new String [this.cameraIPs.length];
			this.imageURLs=new String [this.cameraIPs.length];
			this.metaURLs= new String [this.cameraIPs.length];
			this.triggerURL="http://"+this.cameraIPs[this.masterSubCamera]+":"+(this.imgsrvPort+this.masterPort)+"/"+triggerURLcmd;
			this.images=   new ImagePlus[this.channelMap.length];
			this.imagesIP= new ImagePlus[this.cameraIPs.length];
			
			for (int i=0; i<this.cameraIPs.length;i++){
				this.jp4_Instances[i]=new JP46_Reader_camera(false);// invisible
				this.jp4_Instances[i].camera_url="http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/";
				this.jp4_Instances[i].camera_img=    this.imageURLcmd; // not currently used
				this.jp4_Instances[i].camera_img_new=this.imageURLcmd; //"torp/wait/" will survive, only "towp/wait/" is removed for Exif re-read
				this.jp4_Instances[i].ABSOLUTELY_SILENT=true; 
				this.resetURLs[i]="http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/"+resetURLcmd;
				this.imageURLs[i]="http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/"+this.imageURLcmd;
				this.metaURLs[i]= "http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/"+metaURLcmd;
				this.imagesIP[i]= null;

			}
			for (int i=0; i<this.images.length;i++) this.images[i]= null;
		}
		
		private void initCamParsDefaultArrays(int num){
			this.cameraExposureCorr=new double[num];
			this.cameraGainCorr=    new double[num];
			this.cameraRScaleCorr=  new double[num];
			this.cameraBScaleCorr=  new double[num];
			this.cameraGScaleCorr=  new double[num];
			this.cameraExtraURL=    new String[num];
			for (int i=0; i<num;i++){
				this.cameraExposureCorr[i]=1.0;
				this.cameraGainCorr[i]=    1.0;
				this.cameraRScaleCorr[i]=  1.0;
				this.cameraBScaleCorr[i]=  1.0;
				this.cameraGScaleCorr[i]=  1.0;
				this.cameraExtraURL[i]=    "";
			}
		}

			
		/**
		 * Initialize cameraIPs from subNet and baseIP (sequentially)
		 */
		private void initIPs(){
			if (this.debugLevel>2) System.out.println("initIPs(): this.iBaseIP=" + this.iBaseIP );
			int size=0;
			for (int i=0;i<this.channelMap.length;i++) if (this.channelMap[i][0]>size) size=this.channelMap[i][0];
			size++;
			this.cameraIPs=new String [size];
			for (int i=0;i<size;i++) this.cameraIPs[i]=this.cameraSubnet+(this.iBaseIP+i);
//			this.masterSubCamera=0;
			
		}
		/**
		 * Initialize default subcamera map
		 * @param size number of subcameras
		 */
		private void initDefaultMap(int size){
			this.channelMap=new int [size][];
			this.flipImages=new boolean[size];
			int []port_seq={1,0,2,3};
			this.masterSubCamera=0;
			this.masterPort=0;
			if (this.nc393){
				this.imgsrvPort=2323;
				this.resetURLcmd="towp/save/pointers"; // advance buffer, next time will wait for the next frame acquired
				this.imageURLcmd="torp/wait/timestamp_name/bimg"; // will wait if needed. If repeated (as when reading Exif)- won't wait
				this.metaURLcmd="torp/wait/meta";  // will get XML, including timestamp
				if (size == 26) {
					for (int i=0;i<size;i++){
						this.channelMap[i]=channelMap_23_393[i].clone();
						this.flipImages[i]=false;
					}
					this.masterSubCamera = 2;
					this.masterPort=2;
				} else for (int i=0;i<size;i++){
					this.flipImages[i]=false;
					this.channelMap[i]=new int[3];
					this.channelMap[i][0]= i >> 2;
					this.channelMap[i][1]= 0;
					if (size <4) {
						this.channelMap[i][2]= i;
					} else {
						this.channelMap[i][2]= port_seq[i & 3];
					}
				}
			} else {
				if (size==1) { // single camera - old lens focusing
					this.channelMap[0]=channelMap1[0].clone();
					this.flipImages[0]=true;
				} else if (size==2){ // New lens focusing machine
					this.channelMap[0]=channelMap2[0].clone();
					this.flipImages[0]=true;  // main sensor under test
					this.channelMap[1]=channelMap2[1].clone();
					this.flipImages[1]=false; // extra sensor for location

				} else if (size==3){
					for (int i=0;i<size;i++){
						this.flipImages[i]=false;
						int i0=((i>=this.channelMap3.length)?(this.channelMap3.length-1):i);
						//					 this.channelMap[i]=this.channelMap21[i0].clone();
						this.channelMap[i]=this.channelMap3[i0].clone();
					}
				} else for (int i=0;i<size;i++){
					this.flipImages[i]=false;
					int i0=((i>=this.channelMap21.length)?(this.channelMap21.length-1):i);
					this.channelMap[i]=this.channelMap21[i0].clone();
				}
			}
		}

		public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"cameraSubnet",this.cameraSubnet);
    		properties.setProperty(prefix+"iBaseIP",this.iBaseIP+"");
    		properties.setProperty(prefix+"masterSubCamera",this.masterSubCamera+"");
    		properties.setProperty(prefix+"masterPort",this.masterPort+"");
    		properties.setProperty(prefix+"cameraBootTimeSeconds",this.cameraBootTimeSeconds+"");
    		properties.setProperty(prefix+"connectionTimeoutMilliseconds",this.connectionTimeoutMilliseconds+"");
    		properties.setProperty(prefix+"imgsrvPort",this.imgsrvPort+"");
    		properties.setProperty(prefix+"resetURLcmd",this.resetURLcmd);
    		properties.setProperty(prefix+"triggerURLcmd",this.triggerURLcmd);
    		properties.setProperty(prefix+"imageURLcmd",this.imageURLcmd);
    		properties.setProperty(prefix+"metaURLcmd",this.metaURLcmd);
    		properties.setProperty(prefix+"channelMap.length",this.channelMap.length+"");
    		for (int i=0;i<this.channelMap.length;i++) {
            		properties.setProperty(prefix+"channelMap_"+i+"_IPindex",   this.channelMap[i][0]+"");
            		properties.setProperty(prefix+"channelMap_"+i+"_subchannel",this.channelMap[i][1]+"");
            		properties.setProperty(prefix+"channelMap_"+i+"_port",      this.channelMap[i][2]+"");
            		properties.setProperty(prefix+"flipImages_"+i              ,this.flipImages[i]?"1":"0");
    		}
    		
    		properties.setProperty(prefix+"cameraIPs.length",this.cameraIPs.length+"");
    		properties.setProperty(prefix+"colorMode",this.colorMode+"");
    		properties.setProperty(prefix+"noWait",this.noWait+"");
    		properties.setProperty(prefix+"nc393",this.nc393+"");
    		properties.setProperty(prefix+"debugSensorNumber",this.debugSensorNumber+"");
    		properties.setProperty(prefix+"JPEGquality",this.JPEGquality+"");
    		properties.setProperty(prefix+"cameraAutoExposure",this.cameraAutoExposure+"");
    		properties.setProperty(prefix+"cameraAutoWhiteBalance",this.cameraAutoWhiteBalance+"");
    		properties.setProperty(prefix+"cameraExtraURLCommon","<![CDATA["+this.cameraExtraURLCommon+"]]>");
    		properties.setProperty(prefix+"setupTriggerMode",this.setupTriggerMode+"");
    		properties.setProperty(prefix+"externalTriggerCabling",this.externalTriggerCabling+"");
    		properties.setProperty(prefix+"noCabling",this.noCabling+"");
    		properties.setProperty(prefix+"cameraExposure",this.cameraExposure+"");
    		properties.setProperty(prefix+"scaleExposureForLasers",this.scaleExposureForLasers+"");
    		properties.setProperty(prefix+"scaleExposureForHeadLasers",this.scaleExposureForHeadLasers+"");
    		properties.setProperty(prefix+"cameraAutoExposureMax",cameraAutoExposureMax+"");
    		properties.setProperty(prefix+"cameraGain",this.cameraGain+"");
    		properties.setProperty(prefix+"cameraRScale",this.cameraRScale+"");
    		properties.setProperty(prefix+"cameraBScale",this.cameraBScale+"");
    		properties.setProperty(prefix+"cameraGScale",this.cameraGScale+"");
			for (int i=0;i<this.cameraIPs.length;i++){
        		properties.setProperty(prefix+"cameraIPs_"+i,this.cameraIPs[i]+"");
        		properties.setProperty(prefix+"cameraExposureCorr_"+i,this.cameraExposureCorr[i]+"");
        		properties.setProperty(prefix+"cameraGainCorr_"+i,   this.cameraGainCorr[i]+"");
        		properties.setProperty(prefix+"cameraRScaleCorr_"+i,   this.cameraRScaleCorr[i]+"");
        		properties.setProperty(prefix+"cameraBScaleCorr_"+i,   this.cameraBScaleCorr[i]+"");
        		properties.setProperty(prefix+"cameraGScaleCorr_"+i,   this.cameraGScaleCorr[i]+"");
        		properties.setProperty(prefix+"cameraExtraURL_"+i,   "<![CDATA["+this.cameraExtraURL[i]+"]]>");
        		
			}
    	}
     	public void getProperties(String prefix,Properties properties){
    		if (properties.getProperty(prefix+"cameraSubnet")!=null)
    			this.cameraSubnet=properties.getProperty(prefix+"cameraSubnet");
    		if (properties.getProperty(prefix+"iBaseIP")!=null)
    			this.iBaseIP=Integer.parseInt(properties.getProperty(prefix+"iBaseIP"));
    		if (properties.getProperty(prefix+"masterSubCamera")!=null)
    			this.masterSubCamera=Integer.parseInt(properties.getProperty(prefix+"masterSubCamera"));
    		if (properties.getProperty(prefix+"masterPort")!=null)
    			this.masterPort=Integer.parseInt(properties.getProperty(prefix+"masterPort"));
    		if (properties.getProperty(prefix+"cameraBootTimeSeconds")!=null)
    			this.cameraBootTimeSeconds=Integer.parseInt(properties.getProperty(prefix+"cameraBootTimeSeconds"));
    		if (properties.getProperty(prefix+"connectionTimeoutMilliseconds")!=null)
    			this.connectionTimeoutMilliseconds=Integer.parseInt(properties.getProperty(prefix+"connectionTimeoutMilliseconds"));
    		if (properties.getProperty(prefix+"imgsrvPort")!=null)
    			this.imgsrvPort=Integer.parseInt(properties.getProperty(prefix+"imgsrvPort"));

    		if (properties.getProperty(prefix+"resetURLcmd")!=null)
    			this.resetURLcmd=properties.getProperty(prefix+"resetURLcmd");
    		if (properties.getProperty(prefix+"triggerURLcmd")!=null)
    			this.triggerURLcmd=properties.getProperty(prefix+"triggerURLcmd");
    		if (properties.getProperty(prefix+"imageURLcmd")!=null)
    			this.imageURLcmd=properties.getProperty(prefix+"imageURLcmd");
    		if (properties.getProperty(prefix+"metaURLcmd")!=null)
    			this.metaURLcmd=properties.getProperty(prefix+"metaURLcmd");
    		if (properties.getProperty(prefix+"channelMap.length")!=null) {
    			initDefaultMap (Integer.parseInt(properties.getProperty(prefix+"channelMap.length")));
    			this.flipImages=new boolean[this.channelMap.length];
        		for (int i=0;i<this.channelMap.length;i++) {
            		if (properties.getProperty(prefix+"channelMap_"+i+"_IPindex")!=null) 
            			this.channelMap[i][0]=Integer.parseInt(properties.getProperty(prefix+"channelMap_"+i+"_IPindex"));
            		if (properties.getProperty(prefix+"channelMap_"+i+"_subchannel")!=null) 
            			this.channelMap[i][1]=Integer.parseInt(properties.getProperty(prefix+"channelMap_"+i+"_subchannel"));
            		if (properties.getProperty(prefix+"channelMap_"+i+"_port")!=null) 
            			this.channelMap[i][2]=Integer.parseInt(properties.getProperty(prefix+"channelMap_"+i+"_port"));
            		if (properties.getProperty(prefix+"flipImages_"+i)!=null) 
            			this.flipImages[i]=(Integer.parseInt(properties.getProperty(prefix+"flipImages_"+i))>0);
        			
        		}
   			
    		}
    		int numCams=0;
    		if (properties.getProperty(prefix+"cameraIPs.length")!=null) {
    			numCams=Integer.parseInt(properties.getProperty(prefix+"cameraIPs.length"));
    			this.cameraIPs=new String[numCams];
        		for (int i=0;i<numCams;i++) {
            		if (properties.getProperty(prefix+"cameraIPs_"+i)!=null) 
            			this.cameraIPs[i]=properties.getProperty(prefix+"cameraIPs_"+i);
        			
        		}
    		}
    		initCamParsDefaultArrays(numCams);
    		
			if (properties.getProperty(prefix+"colorMode")!=null)
				this.colorMode=Integer.parseInt(properties.getProperty(prefix+"colorMode"));
			if (properties.getProperty(prefix+"debugSensorNumber")!=null)
				this.debugSensorNumber=Integer.parseInt(properties.getProperty(prefix+"debugSensorNumber"));
			if (properties.getProperty(prefix+"JPEGquality")!=null)
				this.JPEGquality=Integer.parseInt(properties.getProperty(prefix+"JPEGquality"));
			if (properties.getProperty(prefix+"noWait")!=null)
				this.noWait=Boolean.parseBoolean(properties.getProperty(prefix+"noWait"));
			if (properties.getProperty(prefix+"nc393")!=null)
				this.nc393=Boolean.parseBoolean(properties.getProperty(prefix+"nc393"));
			if (properties.getProperty(prefix+"cameraAutoExposure")!=null)
				this.cameraAutoExposure=Boolean.parseBoolean(properties.getProperty(prefix+"cameraAutoExposure"));
			if (properties.getProperty(prefix+"cameraAutoWhiteBalance")!=null)
				this.cameraAutoWhiteBalance=Boolean.parseBoolean(properties.getProperty(prefix+"cameraAutoWhiteBalance"));
    		if (properties.getProperty(prefix+"cameraExtraURLCommon")!=null) {
    			this.cameraExtraURLCommon=properties.getProperty(prefix+"cameraExtraURLCommon");
				if ((this.cameraExtraURLCommon.length()>10) && this.cameraExtraURLCommon.substring(0,9).equals("<![CDATA["))
					this.cameraExtraURLCommon=this.cameraExtraURLCommon.substring(9,this.cameraExtraURLCommon.length()-3); 
    		}
			if (properties.getProperty(prefix+"setupTriggerMode")!=null)
				this.setupTriggerMode=Boolean.parseBoolean(properties.getProperty(prefix+"setupTriggerMode"));
			if (properties.getProperty(prefix+"externalTriggerCabling")!=null)
				this.externalTriggerCabling=Boolean.parseBoolean(properties.getProperty(prefix+"externalTriggerCabling"));
			if (properties.getProperty(prefix+"noCabling")!=null)
				this.noCabling=Boolean.parseBoolean(properties.getProperty(prefix+"noCabling"));
			if (properties.getProperty(prefix+"cameraExposure")!=null)
				this.cameraExposure=Double.parseDouble(properties.getProperty(prefix+"cameraExposure"));
			if (properties.getProperty(prefix+"scaleExposureForLasers")!=null)
				this.scaleExposureForLasers=Double.parseDouble(properties.getProperty(prefix+"scaleExposureForLasers"));
			if (properties.getProperty(prefix+"scaleExposureForHeadLasers")!=null)
				this.scaleExposureForHeadLasers=Double.parseDouble(properties.getProperty(prefix+"scaleExposureForHeadLasers"));
			
			if (properties.getProperty(prefix+"cameraAutoExposureMax")!=null)
				this.cameraAutoExposureMax=Double.parseDouble(properties.getProperty(prefix+"cameraAutoExposureMax"));
			if (properties.getProperty(prefix+"cameraGain")!=null)
				this.cameraGain=Double.parseDouble(properties.getProperty(prefix+"cameraGain"));
			if (properties.getProperty(prefix+"cameraRScale")!=null)
				this.cameraRScale=Double.parseDouble(properties.getProperty(prefix+"cameraRScale"));
			if (properties.getProperty(prefix+"cameraBScale")!=null)
				this.cameraBScale=Double.parseDouble(properties.getProperty(prefix+"cameraBScale"));
			if (properties.getProperty(prefix+"cameraGScale")!=null)
				this.cameraGScale=Double.parseDouble(properties.getProperty(prefix+"cameraGScale"));
    		
    		for (int i=0;i<numCams;i++) {
        		if (properties.getProperty(prefix+"cameraExposureCorr_"+i)!=null) 
        			this.cameraExposureCorr[i]=Double.parseDouble(properties.getProperty(prefix+"cameraExposureCorr_"+i));
        		if (properties.getProperty(prefix+"cameraGainCorr_"+i)!=null) 
        			this.cameraGainCorr[i]=Double.parseDouble(properties.getProperty(prefix+"cameraGainCorr_"+i));
        		if (properties.getProperty(prefix+"cameraRScaleCorr_"+i)!=null) 
        			this.cameraRScaleCorr[i]=Double.parseDouble(properties.getProperty(prefix+"cameraRScaleCorr_"+i));
        		if (properties.getProperty(prefix+"cameraBScaleCorr_"+i)!=null) 
        			this.cameraBScaleCorr[i]=Double.parseDouble(properties.getProperty(prefix+"cameraBScaleCorr_"+i));
        		if (properties.getProperty(prefix+"cameraGScaleCorr_"+i)!=null) 
        			this.cameraGScaleCorr[i]=Double.parseDouble(properties.getProperty(prefix+"cameraGScaleCorr_"+i));
        		if (properties.getProperty(prefix+"cameraExtraURL_"+i)!=null) {
        			this.cameraExtraURL[i]=properties.getProperty(prefix+"cameraExtraURL_"+i);
    				if ((this.cameraExtraURL[i].length()>10) && this.cameraExtraURL[i].substring(0,9).equals("<![CDATA["))
    					this.cameraExtraURL[i]=this.cameraExtraURL[i].substring(9,this.cameraExtraURL[i].length()-3); 
        		}
        		
    		}
    		initJP4();
     	}   	
    	
	   	public boolean showDialog(String title, int numCams, boolean askRegenerate) { // numCams<=0 -> do not initialize number
	   		if (numCams>0){
//				 initDefaultMap(1);
				 initDefaultMap(numCams);
				 initIPs();
		         initCamParsDefaultArrays(this.cameraIPs.length);
		    	 initJP4();
	   		}
//http://192.168.0.223/parsedit.php?title=PHASES&SENSOR_PHASE&MULTI_PHASE1&MULTI_PHASE2&MULTI_PHASE3&DEBUG&CABLE_TIM&FPGA_TIM0&FPGA_TIM1&DLY359_P1&DLY359_P2&DLY359_P3&DLY359_C1&DLY359_C2&DLY359_C3&SENSOR&FRAME_SIZE=@&refresh	   		
	   		//http://192.168.0.228/parsedit.php?immediate&EXPOS=25000*0&AUTOEXP_ON=0*0&WB_EN=0*0
	   		//http://192.168.0.221/parsedit.php?immediate&EXPOS=25000*0&AUTOEXP_ON=0*0&WB_EN=0*0&COLOR=1*0&GAINR=0x20000*0&GAING=0x20000*0&GAINB=0x36666*0&GAINGB=0x20000*0
	   		//http://192.168.0.221/parsedit.php?title=PHASES&SENSOR_PHASE&MULTI_PHASE1&MULTI_PHASE2&MULTI_PHASE3&FRAME_SIZE=@
//http://192.168.0.222/parsedit.php?title=day&immediate&EXPOS=20000*0&AUTOEXP_ON=0*0&WB_EN=0*0&COLOR=1*0&GAINR=0x270a4*0&GAING=0x20000*0&GAINB=0x2999a*0&GAINGB=0x20000*0
//http://192.168.0.236/parsedit.php?title=night&immediate&EXPOS=13000*0&AUTOEXP_ON=0*0&WB_EN=0*0&COLOR=1*0&GAINR=0x1fb35*0&GAING=0x20000*0&GAINB=0x3961c*0&GAINGB=0x20000*0
//http://192.168.0.221/parsedit.php?title=bright-day&immediate&EXPOS=10000*0&AUTOEXP_ON=0*0&WB_EN=0*0&COLOR=1*0&GAINR=0x270a4*0&GAING=0x20000*0&GAINB=0x2999a*0&GAINGB=0x20000*0	   		
//http://192.168.0.221/parsedit.php?immediate&TRIG&TRIG_PERIOD	   		
    		GenericDialog gd = new GenericDialog(title);
    		gd.addCheckbox    ("NC393 (unchecked - nc353)", this.nc393);
    		gd.addStringField ("Subnet of the cameras (3 first of the four IPv4 address ending with '.')",this.cameraSubnet,12);
    		gd.addNumericField("Last byte of the first sub-camera IP address)",this.iBaseIP,0);
    		gd.addNumericField("Camera boot time",this.cameraBootTimeSeconds,0,3,"sec");
    		gd.addNumericField("Network connection timeout",this.connectionTimeoutMilliseconds,0,3,"ms");
    		
    		gd.addNumericField("Index (in IP table) of the master camera (used for triggering)",this.masterSubCamera,0);
    		gd.addNumericField("Master (used for triggering), normally lowest connected port (0..3)",this.masterPort,0);
    		gd.addNumericField("Image server port number)",this.imgsrvPort,0);
    		gd.addStringField ("Image server command to reset image buffer",this.resetURLcmd,25);
    		gd.addStringField ("Image server command to trigger acquisition",this.triggerURLcmd,25);
    		gd.addStringField ("Image server command to acquire image (waits for the new one after reset)",this.imageURLcmd,25);
    		gd.addStringField ("Image server command to receive XML metadata (with timestamp)",this.metaURLcmd,25);
    		gd.addMessage("Configure each sub-camera - which IP index and channel does it use.");
    		if (askRegenerate) gd.addMessage("You may change number of subcameras and press REGENERATE below");
    		for (int i=0; i< this.channelMap.length;i++){
    			gd.addMessage("---------------------------------------------------------");
        		gd.addNumericField("Subcamera "+(i+1)+" IP index (starting from 0)", this.channelMap[i][0],0);
        		gd.addNumericField("Subcamera "+(i+1)+" port (0..3)",                this.channelMap[i][2],0);
        		gd.addNumericField("Subcamera "+(i+1)+" channel (0,1 or 2)",         this.channelMap[i][1],0);
        		gd.addCheckbox("Subcamera "+(i+1)+" - used mirror",                  this.flipImages[i]);
    		}
			gd.addMessage("---------------------------------------------------------");
    		if (askRegenerate) {
        		gd.addCheckbox ("Individually overwrite IP addresses",false);
    			gd.addNumericField("Number of subcameras (need to press REGENERATE button to change)", this.channelMap.length, 0);
    			gd.enableYesNoCancel("OK", "REGENERATE");
    		}
    	    WindowTools.addScrollBars(gd);
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
    	    this.nc393=                     gd.getNextBoolean();
    	    this.cameraSubnet=              gd.getNextString();
    	    if (!this.cameraSubnet.endsWith(".")) this.cameraSubnet += ".";
    	    int bip=                  (int) gd.getNextNumber(); 
    	    if ((bip>0) && (bip<255) && (bip!=this.iBaseIP)){
    	    	this.iBaseIP=bip;
    	    	initIPs();
    	    }
    		this.cameraBootTimeSeconds=(int) gd.getNextNumber();
    		this.connectionTimeoutMilliseconds=(int) gd.getNextNumber();
    	    this.masterSubCamera=     (int) gd.getNextNumber(); // should be after initIPs()!
    	    this.masterPort=          (int) gd.getNextNumber(); // should be after initIPs()!
    	    
    	    this.imgsrvPort=          (int) gd.getNextNumber();
    	    this.resetURLcmd=            gd.getNextString();
    	    this.triggerURLcmd=          gd.getNextString();
    	    this.imageURLcmd=            gd.getNextString();
    	    this.metaURLcmd=             gd.getNextString();
    		for (int i=0; i< this.channelMap.length;i++){
        		this.channelMap[i][0]=(int) gd.getNextNumber();
        		this.channelMap[i][2]=(int) gd.getNextNumber();
        		this.channelMap[i][1]=(int) gd.getNextNumber();
        		this.flipImages[i]=         gd.getNextBoolean();
    		}
    		if (askRegenerate) {
    			boolean overwriteIPs=gd.getNextBoolean();
    			int newSubCams=(int) gd.getNextNumber();
//    			if (!gd.wasOKed() &&(newSubCams!=this.channelMap.length)) {
       			if (!gd.wasOKed()) { // regenerate always if pressed, even if number is the same
//    				 initDefaultMap(newSubCams);
//    				 initIPs();
    				 int [][] backupMap=this.channelMap;
    				 if (!showDialog(title, newSubCams, false)) {
    					 this.channelMap=backupMap;
    					 return false; // channelMap is restored, but some other fields may be broken
    				 }
    			}
    			if (overwriteIPs && !editSubCamerasIPs()) return false;
    		}
	   		return true;
	   	}
	   	/**
	   	 * Wait for camera to boot
	   	 * @param chn - camera channel number (0)
	   	 * @return  -1 if camera is not yet set to trigger mode, else - last frame number acquired
	   	 * throws on timeout
	   	 */
	   	public int getCurrentFrameNumberWithTimeout(int chn, boolean showStatus, AtomicInteger stopRequested){
			long startTime=System.nanoTime();
			double dTimeOut=1E9* this.cameraBootTimeSeconds;
			long endTime=startTime+(long) dTimeOut;
			while (true){
				if (probeCameraState(chn, showStatus, true,this.connectionTimeoutMilliseconds)){
					if (showStatus) IJ.showProgress(0.0);
					if (!this.triggeredMode[chn]) return -1;
					if (this.triggerPeriod[chn]>1) return -1; // trigger period is not set to 1 - not yet programmed
					return this.cameraFrameNumber[chn];
				}
				long time=System.nanoTime();
				if (time>endTime){
					if (showStatus) IJ.showProgress(0.0);
					throw new IllegalArgumentException ("Timeout while waiting for the camera #"+chn+" to respond");
				}
				if (stopRequested.get()>0) {
					System.out.println("User requested stop");
					if (showStatus) IJ.showProgress(0.0);
					throw new IllegalArgumentException ("Waiting for camera #"+chn+" aborted by user request");
				}
				if (showStatus) IJ.showProgress((time-startTime)/dTimeOut);
			}
	   	}
	   	
	   	
	   	public String getSerialNumber(int chn, int EEPROM_chn){
	   		String url="http://"+this.cameraIPs[chn]+"/i2c.php?cmd=fromEEPROM0&EEPROM_chn="+EEPROM_chn;
	   			
	   	    	Document dom=null;
	   	    	String serial=null;
	   	    	try {
	   	    		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   	    		DocumentBuilder db = dbf.newDocumentBuilder();
	   	    		dom = db.parse(url);
	   	    		if (!dom.getDocumentElement().getNodeName().equals("board")) {
	   	    			String msg="Root element: expected 'board', got \"" + dom.getDocumentElement().getNodeName()+"\"";
	   	    			IJ.showMessage("Error",msg); 
    					throw new IllegalArgumentException (msg);
	   	    		}
	   	    		
	   	    		serial=((Node) (((Node) dom.getDocumentElement().getElementsByTagName("serial").item(0)).getChildNodes().item(0))).getNodeValue();
    				// remove opening and closing "
    				if (serial==null){
    					String msg="Could not read tag <serial>";
    					IJ.showMessage("Error",msg); 
    					System.out.println(msg);
    					return null;
    				}
    				if (serial.startsWith("\"")){
    					serial=serial.substring(1, serial.length()-1);
    				}
    			} catch(MalformedURLException e){
    				String msg="Please check the URL:" + e.toString();
					IJ.showMessage("Error",msg); 
					throw new IllegalArgumentException (msg);
    			} catch(IOException  e1){
    				String msg = e1.getMessage();
    				if (msg==null || msg.equals(""))  msg = ""+e1;
					IJ.showMessage("Error",msg); 
					throw new IllegalArgumentException (msg);
    			}catch(ParserConfigurationException pce) {
    				pce.printStackTrace();
    				throw new IllegalArgumentException ("PCE error");
    			}catch(SAXException se) {
    				se.printStackTrace(); 
    				throw new IllegalArgumentException ("SAX error");
    			}
    		return serial;
	   	}

	   	public double getSensorTemperature(int chn, int EEPROM_chn){
	   		// Need to configure camera to be able to read temperature
	   		if (sensorPresent==null) { // cameras were not probed/configured
				probeCameraState(); // testing detection
		   		printTiming("=== probeCameraState()");
				setupCameraAcquisition();
		   		printTiming("=== setupCameraAcquisition()");
	   		}
	   		if (!this.sensorPresent[chn][0] && !this.sensorPresent[chn][1] && !this.sensorPresent[chn][2]) EEPROM_chn=0; // no 10359 - null pointer while "lens center" if first
	   		String url="http://"+this.cameraIPs[chn]+"/i2c.php?cmd=fromEEPROM0&EEPROM_chn="+EEPROM_chn;
	   			this.lastTemperature=Double.NaN;
	   	    	Document dom=null;
	   	    	try {
	   	    		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   	    		DocumentBuilder db = dbf.newDocumentBuilder();
	   	    		dom = db.parse(url);
	   	    		if (!dom.getDocumentElement().getNodeName().equals("board")) {
	   	    			String msg="Root element: expected 'board', got \"" + dom.getDocumentElement().getNodeName()+"\"";
	   	    			IJ.showMessage("Error",msg); 
    					throw new IllegalArgumentException (msg);
	   	    		}
	   	    		
	   	    		String sTemperature=((Node) (((Node) dom.getDocumentElement().getElementsByTagName("sensorTemperature").item(0)).getChildNodes().item(0))).getNodeValue();
    				// remove opening and closing "
    				if (sTemperature==null){
    					String msg="Could not read sensor temperature";
//    					IJ.showMessage("Error",msg); 
    					System.out.println("Warning: "+msg);
    					return Double.parseDouble(sTemperature);
    				}
    				this.lastTemperature= Double.parseDouble(sTemperature);
    			} catch(MalformedURLException e){
    				String msg="Please check the URL:" + e.toString();
					IJ.showMessage("Error",msg); 
					throw new IllegalArgumentException (msg);
    			} catch(IOException  e1){
    				String msg = e1.getMessage();
    				if (msg==null || msg.equals(""))  msg = ""+e1;
					IJ.showMessage("Error",msg); 
					throw new IllegalArgumentException (msg);
    			}catch(ParserConfigurationException pce) {
    				pce.printStackTrace();
    				throw new IllegalArgumentException ("PCE error");
    			}catch(SAXException se) {
    				se.printStackTrace(); 
    				throw new IllegalArgumentException ("SAX error");
    			}
				return this.lastTemperature;
	   	}
	   	
	   	
//		private int [] motors=               null; // motors steps when the images were acquired (for null)
	   	public int [] getMotorsPosition() {return this.motorsPosition;} // may be null
	   	public void setMotorsPosition (int [] position) {this.motorsPosition=position;}
	   	public void resetInitialization(){
	   		System.out.println("resetInitialization()");
	   		this.sensorPresent=null;
	   	}
	   	public int probeCameraState(){
	   		int numOnline=0;
	   		int numSensors=0;
	   		int numChn=this.cameraIPs.length;
	   		initCameraArrays(numChn);
	   		/*
	   		this.cameraFrameNumber=new int[numChn];
	   		this.triggeredMode=new boolean[numChn];
	   		this.sensorPresent=new boolean[numChn][];
			// TODO - try if skipping setting TRIG_PERIOD=0, IRQ_SMART=6 (when they are already set) will fix hanging
	   		this.triggerPeriod=new int[numChn];
	   		this.irqSmart=     new int[numChn];
	   		*/
	   		
	   		
	   		for (int chn=0;chn<numChn;chn++){
		   		this.cameraFrameNumber[chn]=-1;
		   		this.triggeredMode[chn]=false;
		   		this.sensorPresent[chn]=null; // camera did not respond, {false,false,false} - single camera (no 10359)
		   		if (this.debugLevel>1) System.out.println("Probing camera "+chn+" ("+this.cameraIPs[chn]+")...");
		   		IJ.showStatus("Probing camera "+chn+" ("+this.cameraIPs[chn]+")...");
	   			if (probeCameraState(chn)) {
	   				numOnline++;
		   			printTiming("===== probing channel "+chn);
			   		if (this.debugLevel>1) System.out.println("Frame number: "+this.cameraFrameNumber[chn]+
			   				", Trigger mode:"+(this.triggeredMode[chn]?"ex":"in")+"ternal, "+
			   				" sensors attached:"+(this.sensorPresent[chn][0]?"1 ":"")+(this.sensorPresent[chn][1]?"2 ":"")+(this.sensorPresent[chn][2]?"3 ":"")+
			   				((!this.sensorPresent[chn][0] && !this.sensorPresent[chn][1] && !this.sensorPresent[chn][2])?"single-sensor, no multiplexer":"")
			   				);
			   		for (int i=0;i<this.sensorPresent[chn].length;i++) if (this.sensorPresent[chn][i]) numSensors++;
	   			} else {
	   				if (this.debugLevel>1) System.out.println("Camera did not respond");
	   			}
	   		}
	   		IJ.showStatus("Found "+numOnline+" camera IP, "+numSensors+" sensors");
			if (this.debugLevel>1) System.out.println("Found "+numOnline+" camera IP, "+numSensors+" sensors");
	   		return numOnline;
	   	}
	   	public void initCameraArrays(int numChn){
	   		if (this.debugLevel>1) System.out.println("initCameraArrays("+numChn+")");
	   		this.cameraFrameNumber=new int[numChn];
	   		this.triggeredMode=new boolean[numChn];
	   		this.sensorPresent=new boolean[numChn][];
			// TODO - try if skipping setting TRIG_PERIOD=0, IRQ_SMART=6 (when they are already set) will fix hanging
	   		this.triggerPeriod=new int[numChn];
	   		this.irqSmart=     new int[numChn];
	   	}
	   	public boolean probeCameraState(
	               int chn){
	                 return probeCameraState(chn, false, false,0);
	               }

		   	public boolean probeCameraState(
	               int chn,
	               boolean showStatus,
	               boolean doNotThrow,
	               int timeout // ms
	               ){
	   		//http://192.168.0.221/parsedit.php?immediate&TRIG&TRIG_PERIOD&FRAME
	   		String url;
	   		if (this.nc393){
	   			url="http://"+this.cameraIPs[chn]+"/parsedit.php?immediate&TRIG&TRIG_PERIOD&SENS_AVAIL&FRAME";
	   		} else {
	   			url="http://"+this.cameraIPs[chn]+"/parsedit.php?immediate&TRIG&TRIG_PERIOD&IRQ_SMART&SENS_AVAIL&FRAME";
	   		}
	   		if (this.debugLevel>1) System.out.println("url="+url);
	   		Document dom=null;
	   		if ((this.cameraFrameNumber==null) || ((this.cameraFrameNumber.length<(chn+1)))) initCameraArrays(chn+1);
	   		// did not find yet - why cameraFrameNumber is not null, but sensorPresent is? Added next lines until find
	   		if ((this.triggeredMode==null) || ((this.triggeredMode.length<(chn+1)))) initCameraArrays(chn+1);
	   		if ((this.sensorPresent==null) || ((this.sensorPresent.length<(chn+1)))) initCameraArrays(chn+1);
	   		if ((this.triggerPeriod==null) || ((this.triggerPeriod.length<(chn+1)))) initCameraArrays(chn+1);
	   		if (!this.nc393) {
	   			if ((this.irqSmart==null) || ((this.irqSmart.length<(chn+1)))) initCameraArrays(chn+1);
	   		}
	   		try {
	   			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   			DocumentBuilder db = dbf.newDocumentBuilder();
	   			if (timeout<=0){
	   			  dom = db.parse(url);
	   			} else {
	   	            URL uUrl = new URL(url);
	   	            URLConnection con = uUrl.openConnection();
	   	            con.setConnectTimeout(timeout);//The timeout in mills
	   	            dom = db.parse(con.getInputStream());
	   			}
	   			if (!dom.getDocumentElement().getNodeName().equals("parameters")) {
	   				String msg="Root element: expected 'parameters', got \"" + dom.getDocumentElement().getNodeName()+"\"";
	   				if (showStatus) IJ.showStatus(msg); 
	   				if (doNotThrow) return false;
	   				IJ.showMessage("Error",msg); 
	   				throw new IllegalArgumentException (msg);
	   			}
	   			//	   	    		String sNode;
	   			if (this.triggeredMode==null) System.out.println("this.triggeredMode==null"); 
	   			if (this.triggeredMode.length<(chn+1))  System.out.println("this.triggeredMode.length="+this.triggeredMode.length+" chn="+chn); 

	   			this.triggeredMode[chn]= (Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("TRIG").item(0)).getChildNodes().item(0))).getNodeValue())==4);
	   			int sensAvail= Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("SENS_AVAIL").item(0)).getChildNodes().item(0))).getNodeValue());

	   			if (this.sensorPresent==null) System.out.println("sensorPresent==null"); 
	   			if (this.sensorPresent.length<(chn+1))  System.out.println("sensorPresent="+this.sensorPresent.length+" chn="+chn);
	   			
	   			this.sensorPresent[chn]=new boolean[3];
	   			for (int i=0;i<this.sensorPresent[chn].length;i++) this.sensorPresent[chn][i]=(sensAvail & (1<<i))!=0;
	   			this.triggerPeriod[chn]=     Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("TRIG_PERIOD").item(0)).getChildNodes().item(0))).getNodeValue());
	   			if (!this.nc393) {
	   				this.irqSmart[chn]=          Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("IRQ_SMART").item(0)).getChildNodes().item(0))).getNodeValue());
	   			}
	   			this.cameraFrameNumber[chn]= Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("FRAME").item(0)).getChildNodes().item(0))).getNodeValue());
	   		} catch(MalformedURLException e){
	   			String msg="Please check the URL:" + e.toString();
	   			IJ.showMessage("Error",msg); 
   				if (showStatus) IJ.showStatus(msg); 
   				if (doNotThrow) return false;
	   			throw new IllegalArgumentException (msg);
	   		} catch(IOException  e1){
	   			String msg = e1.getMessage();
	   			if (msg==null || msg.equals(""))  msg = ""+e1;
	   			msg="Camera "+chn+" ("+this.cameraIPs[chn]+") did not respond: "+msg; // "Connection refused"
   				if (showStatus) IJ.showStatus(msg); 
   				if (doNotThrow) return false;
	   			IJ.showMessage("Error",msg); 
	   			System.out.println(msg);
	   			return false;
	   		}catch(ParserConfigurationException pce) {
	   			pce.printStackTrace();
	   			String msg="PCE error "+pce.getMessage();
   				if (showStatus) IJ.showStatus(msg); 
   				if (doNotThrow) return false;
	   			throw new IllegalArgumentException (msg);
	   		}catch(SAXException se) {
	   			se.printStackTrace();
	   			String msg="SAX error "+se.getMessage();
   				if (showStatus) IJ.showStatus(msg); 
   				if (doNotThrow) return false;
	   			throw new IllegalArgumentException (msg);
	   		}
	   		return true;
	   	}
	   	public boolean setupCameraAcquisition(){
	   		return setupCameraAcquisition(Double.NaN);
	   	}
	   	public boolean setupCameraAcquisition(final double exposureScale){
	   		final int ipLength=this.resetURLs.length;
	   		final boolean [] results=new boolean[ipLength];
	   		for (int chn=0;chn<ipLength;chn++) results[chn]=(this.sensorPresent[chn]!=null);
	   		final Thread[] threads = newThreadArray(this.maxNumberOfThreads);
	   		final AtomicInteger chnAtomic = new AtomicInteger(0);
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				public void run() {
	   					for (int chn=chnAtomic.getAndIncrement(); chn<ipLength;chn=chnAtomic.getAndIncrement()) if (results[chn]){
	   						results[chn]=setupCameraAcquisition(chn,exposureScale);
	   					}
	   				}

	   			};
	   		}
	   		startAndJoin(threads);
	   		if (Double.isNaN(exposureScale)){ // full init
	   			int nRepeat=this.setupTriggerMode?4:1;
	   			if (this.nc393) nRepeat++; // is it needed?
	   			for (int i=0;i<nRepeat;i++){
	   				if (this.debugLevel>0) System.out.println((i+1)+" of "+nRepeat+": Triggering cameras to give parameters a chance to propagate");
	   				trigger();
	   				try
	   				{
	   					Thread.sleep( 1000 ); // ms
	   				}
	   				catch ( InterruptedException e )
	   				{
	   					System.out.println( "awakened prematurely" );
	   				}
	   			}
	   			
	   		}
	   		boolean result=true;
	   		for (int chn=0;chn<ipLength;chn++) if (this.sensorPresent[chn]!=null) result &=results[chn];
	   		return result;
	   	}	   	
	   	public boolean setupCameraAcquisitionNoTHreads(double exposureScale){
	   		boolean result=true;
	   		for (int chn=0;chn<this.cameraIPs.length;chn++) if (this.sensorPresent[chn]!=null){
	   			result&=setupCameraAcquisition(chn,exposureScale);
	   		}
	   		// trigger ?
	   		if (Double.isNaN(exposureScale)){ // just a hack - if exposureScale is non-NaN only update exposure
			   trigger();
			//sleep (1000);
			  for ( long  endTime=System.nanoTime()+((long) 1000000000);endTime<System.nanoTime(););
	   		}
	   		return result;
	   	}
	   	// TODO: issue several TRIG pulses after setting parameters?
	   	public boolean setupCameraAcquisition(int chn, double exposureScale){
	   		boolean exposureOnly=!Double.isNaN(exposureScale);
	   		int minExposure=10; // usec
	   		int maxExposure=1000000; //usec
	   		int minGain=(int) (0x10000*1.0);
	   		int maxGain=(int) (0x10000*15.75);
	   		int minScale=0;
	   		int maxScale=(int) (0x10000*4.0);
	   		int exposure= (int) (Math.round(1000*this.cameraExposure*this.cameraExposureCorr[chn])*(exposureOnly?exposureScale:1.0));
	   		int autoExposureMax= (int) (Math.round(1000*this.cameraAutoExposureMax*this.cameraExposureCorr[chn]));
	   		int gain=     (int) (Math.round(0x10000*this.cameraGain*this.cameraGainCorr[chn]));
	   		int rScale=   (int) (Math.round(0x10000*this.cameraRScale*this.cameraRScaleCorr[chn]));
	   		int bScale=   (int) (Math.round(0x10000*this.cameraBScale*this.cameraBScaleCorr[chn]));
	   		int gScale=   (int) (Math.round(0x10000*this.cameraGScale*this.cameraGScaleCorr[chn]));
	   		int autoExp=  this.cameraAutoExposure?1:0;
	   		int autoWB=   this.cameraAutoWhiteBalance?1:0;
	   		String extraURL="";
	   		if (this.cameraExtraURLCommon.length()>0){
	   			extraURL+=(this.cameraExtraURLCommon.substring(0,1).equals("&"))?"":"&"+this.cameraExtraURLCommon;
	   		}
	   		if (this.cameraExtraURL[chn].length()>0){
	   			extraURL+=(this.cameraExtraURL[chn].substring(0,1).equals("&"))?"":"&"+this.cameraExtraURL[chn];
	   		}
	   		if (exposure<minExposure) exposure=minExposure; else if (exposure>maxExposure) exposure=maxExposure;
	   		if (autoExposureMax<minExposure) autoExposureMax=minExposure; else if (autoExposureMax>maxExposure) autoExposureMax=maxExposure;
	   		if (gain<minGain) gain= minGain ; else if (gain> maxGain) gain= maxGain;
	   		if (rScale<minScale) rScale= minScale ; else if (rScale> maxScale) rScale= maxScale;
	   		if (bScale<minScale) bScale= minScale ; else if (bScale> maxScale) bScale= maxScale;
	   		if (gScale<minScale) gScale= minScale ; else if (gScale> maxScale) gScale= maxScale;
	   		
	   		String triggerMode=this.setupTriggerMode?(
	   				"&TRIG_CONDITION="+(this.noCabling?(0):(this.externalTriggerCabling?"0x200000":"0x20000"))+"*0"+
	   				"&TRIG_OUT="+(this.noCabling?(0):(this.externalTriggerCabling?"0x800000":"0x80000"))+"*0"+
	   				"&TRIG=4*3"):"";
	   		if (this.nc393) {
	   			if (this.triggerPeriod[chn]>1)triggerMode+="&TRIG_PERIOD=0*1"; // just stop it if it wasn't already, imgsrv /trig does not set it, only FPGA register
	   		}else {
		   		if (this.triggerPeriod[chn]>1)triggerMode+="&TRIG_PERIOD=1*0"; // just imgsrv /trig does not set it, only FPGA register
	   		}
	   		String url="http://"+this.cameraIPs[chn]+"/parsedit.php?immediate";
	   		url+="&EXPOS="+exposure+"*0"; // always
	   		if (!exposureOnly){
	   			if (!this.nc393) {
	   				if (this.irqSmart[chn]!=(this.noWait?6:3)){
	   					url+="&IRQ_SMART="+(this.noWait?6:3)+"*0";
	   				}
	   			}
	   			url+="&COLOR="+this.colorMode+"*"+(this.nc393?"1":"0")+
	   			"&QUALITY="+this.JPEGquality+"*0"+
	   			"&EXPOS="+exposure+"*0"+
	   			"&AUTOEXP_EXP_MAX="+autoExposureMax+"*0"+
	   			"&AUTOEXP_ON="+autoExp+"*0"+
	   			"&GAING="+gain+"*0"+
	   			"&RSCALE="+rScale+"*0"+
	   			"&BSCALE="+bScale+"*0"+
	   			"&GSCALE="+gScale+"*0"+ // GB/G ratio
	   			"&WB_EN="+autoWB+"*0"+
	   			"&DAEMON_EN_TEMPERATURE=1"+"*0"+
	   			 triggerMode+ // needed here?
	   			"&FRAME"+
	   			extraURL;
	   		}
//	   		if (this.debugLevel>2) {
	   		if (this.debugLevel>1) { // temporary - debugging
	   			System.out.println("Setting camera "+chn+", URL="+url);
	   		}
	   		Document dom=null;
	   		try {
	   			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   			DocumentBuilder db = dbf.newDocumentBuilder();
	   			dom = db.parse(url);
	   			if (!dom.getDocumentElement().getNodeName().equals("parameters")) {
	   				String msg="Root element: expected 'parameters', got \"" + dom.getDocumentElement().getNodeName()+"\"";
	   				IJ.showMessage("Error",msg); 
	   				throw new IllegalArgumentException (msg);
	   			}
	   			//	   	    		String sNode;
	   			if (!exposureOnly) this.cameraFrameNumber[chn]= Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("FRAME").item(0)).getChildNodes().item(0))).getNodeValue());
	   		} catch(MalformedURLException e){
	   			String msg="Please check the URL:" + e.toString();
	   			IJ.showMessage("Error",msg); 
	   			throw new IllegalArgumentException (msg);
	   		} catch(IOException  e1){
	   			String msg = e1.getMessage();
	   			if (msg==null || msg.equals(""))  msg = ""+e1;
	   			msg="Camera "+chn+" ("+this.cameraIPs[chn]+") did not respond\n"+msg;
	   			IJ.showMessage("Error",msg); 
	   			System.out.println(msg);
	   			return false;
	   		}catch(ParserConfigurationException pce) {
	   			pce.printStackTrace();
	   			throw new IllegalArgumentException ("PCE error");
	   		}catch(SAXException se) {
	   			se.printStackTrace(); 
	   			throw new IllegalArgumentException ("SAX error");
	   		}
	   		return true;
	   	}
	   	
	   	
	   	
	   	private boolean editSubCamerasIPs() {
			GenericDialog gd = new GenericDialog("Edit IPs/hosts of the sub cameras");
			for (int i=0;i<this.cameraIPs.length;i++){
	    		gd.addStringField(i+": IP address/host of the subcamera",this.cameraIPs[i],20);
			}
    	    WindowTools.addScrollBars(gd);
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
			for (int i=0;i<this.cameraIPs.length;i++){
	    		this.cameraIPs[i]=gd.getNextString();
			}
			initJP4(); // initialize JP46_Reader_camera class instances
///        	initCamParsDefaultArrays(this.cameraIPs.length); - probably not needed as number of cameras is not changed
			return true;

	   	}
	   	
	   	public boolean editCameraSettings(String title){
			GenericDialog gd = new GenericDialog("title");
			gd.addCheckbox    ("NC393 (unchecked - nc353)", this.nc393);
    		gd.addNumericField("Camera exposure",this.cameraExposure,2,8,"ms");
    		gd.addNumericField("Scale camera exposure for target laser detection (4 lasers)",100*this.scaleExposureForLasers,1,5,"%");
    		gd.addNumericField("Scale camera exposure for optical head laser detection (2 lasers)",100*this.scaleExposureForHeadLasers,1,5,"%");
			gd.addCheckbox    ("Enable autoexposure", this.cameraAutoExposure);
    		gd.addNumericField("Maximal automatic exposure",this.cameraAutoExposureMax,2,8,"ms");
    		gd.addNumericField("Analog gain (2.0 recommended)",this.cameraGain,3,6,"x");
    		gd.addNumericField("Red-to-green color balance",this.cameraRScale,3,6,"x");
    		gd.addNumericField("Blue-to-green color balance",this.cameraBScale,3,6,"x");
    		gd.addNumericField("Green(blue row) to green (red row) color balance",this.cameraGScale,3,6,"x");
// TODO: Make one-time WB by averaging scales from all channels    		
			gd.addCheckbox    ("Enable auto white balance", this.cameraAutoWhiteBalance);
			
			
    		gd.addNumericField("Color mode (JP4=5)",this.colorMode,0,6,"");
    		gd.addNumericField("Compression quality",this.JPEGquality,0,6,"%");
			gd.addCheckbox    ("No wait for the next frame sync (IRQ_SMART=6, unchecked - IRQ_SMART=3)", this.noWait);
			
			gd.addCheckbox    ("Setup trigger mode", this.setupTriggerMode);
			gd.addCheckbox    ("Use external trigger wiring", this.externalTriggerCabling);
			gd.addCheckbox    ("Use internal FPGA trigger (no cables, for single camera)", this.noCabling);
			
//			this.setupTriggerMode
//			this.externalTriggerCabling

			
			
    		gd.addStringField ("Extra URL for all sub-cameras",this.cameraExtraURLCommon,40);
			for (int i=0;i<this.cameraIPs.length;i++){
				gd.addMessage("\n"+i+": camera "+this.cameraIPs[i]+" :");
	    		gd.addNumericField("     Exposure correction ",this.cameraExposureCorr[i],2,8,"x");
	    		gd.addNumericField("     Gain correction     ",this.cameraGainCorr[i],2,8,"x");
	    		gd.addNumericField("     R/G correction      ",this.cameraRScaleCorr[i],2,8,"x");
	    		gd.addNumericField("     B/G correction      ",this.cameraBScaleCorr[i],2,8,"x");
	    		gd.addNumericField("     G(B)/G(R) correction",this.cameraGScaleCorr[i],2,8,"x");
	    		gd.addStringField ("     Extra URL string    ",this.cameraExtraURL[i],40);
			}
    		gd.addNumericField("Debug sensor number (<0 - none)",this.debugSensorNumber,0,6,"");
			
    	    WindowTools.addScrollBars(gd);
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
    	    this.nc393=                          gd.getNextBoolean();
    	    this.cameraExposure=                 gd.getNextNumber();
    	    this.scaleExposureForLasers=    0.01*gd.getNextNumber();
    	    this.scaleExposureForHeadLasers=0.01*gd.getNextNumber();
    	    this.cameraAutoExposure=             gd.getNextBoolean();
    	    this.cameraAutoExposureMax=          gd.getNextNumber();
    	    this.cameraGain=                     gd.getNextNumber();
    	    this.cameraRScale=                   gd.getNextNumber();
    	    this.cameraBScale=                   gd.getNextNumber();
    	    this.cameraGScale=                   gd.getNextNumber();
    	    this.cameraAutoWhiteBalance=         gd.getNextBoolean();
    		this.colorMode=                (int) gd.getNextNumber();
    		this.JPEGquality=              (int) gd.getNextNumber();
    		this.noWait=                         gd.getNextBoolean();
			this.setupTriggerMode=               gd.getNextBoolean();
			this.externalTriggerCabling=         gd.getNextBoolean();
			this.noCabling=                      gd.getNextBoolean();

    	    this.cameraExtraURLCommon=           gd.getNextString();
			for (int i=0;i<this.cameraIPs.length;i++){
	    	    this.cameraExposureCorr[i]=      gd.getNextNumber();
	    	    this.cameraGainCorr[i]=          gd.getNextNumber();
	    	    this.cameraRScaleCorr[i]=        gd.getNextNumber();
	    	    this.cameraBScaleCorr[i]=        gd.getNextNumber();
	    	    this.cameraGScaleCorr[i]=        gd.getNextNumber();
	    	    this.cameraExtraURL[i]=          gd.getNextString();
			}
			this.debugSensorNumber=        (int) gd.getNextNumber();
	   		return true;
	   	}
	   	public void resetCameras(){
	   		resetCameras(selectAllSubcameras());
	   	}
	   	public double[] timestampCameras(){
	   		return timestampCameras(selectAllSubcameras());
	   	}
	   	public void resetIPs(){
	   		resetIPs(selectAllIPs());
	   	}
	   	public double[] timestampIPs(){
	   		return timestampIPs(selectAllIPs());
	   	}
	   	public boolean [] selectAllSubcameras(){
	   		boolean [] all= new boolean[this.channelMap.length];
	   		for (int i=0;i<all.length;i++) all[i]=true;
	   		return all;
	   	}
	   	public boolean [] selectNoneSubcameras(){
	   		boolean [] all= new boolean[this.channelMap.length];
	   		for (int i=0;i<all.length;i++) all[i]=false;
	   		return all;
	   	}
	   	public boolean [] selectAllIPs(){
	   		boolean [] all= new boolean[this.cameraIPs.length];
	   		for (int i=0;i<all.length;i++) all[i]=true;
	   		return all;
	   	}
	   	
	   	public boolean [] selectIPs(boolean [] selectCameras){
	   		boolean [] IPs= new boolean[this.cameraIPs.length];
	   		for (int i=0;i<IPs.length;i++) IPs[i]=false;
	   		for (int i=0;i<selectCameras.length;i++) if (i<this.channelMap.length) {
	   			IPs[this.channelMap[i][0]]=true;
	   		}
	   		return IPs;
	   	}

	   	
	   	public void trigger(){
			try {
				DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				DocumentBuilder db = dbf.newDocumentBuilder();
				String url=this.triggerURL;
				if (this.debugLevel>2) System.out.println(">>> trigger:" + url );
				db.parse(url); // should be some XML (currently discarded)
			} catch(MalformedURLException e){
				System.out.println("Please check the URL:" + e.toString() );
				return;
			} catch(IOException  e1){
				IJ.showStatus("");
				String error = e1.getMessage();
				if (error==null || error.equals(""))  error = ""+e1;
				IJ.showMessage("trigger() ERROR", ""+error);
				return;
			}catch(ParserConfigurationException pce) {
				pce.printStackTrace();
				return;
			}catch(SAXException se) {
				se.printStackTrace(); 
				return;
			}
		}
	   		
	   	public double[] timestampCameras(boolean [] selection){
	   		double [] ts_ip=timestampIPs(selectIPs(selection));
	   		double [] ts=new double [selection.length];
	   		for (int i=0;i<ts.length;i++){
	   			ts[i]=0.0;
	   			if (i<this.channelMap.length) ts[i]=ts_ip[this.channelMap[i][0]];
	   		}
	   		return ts;
	   	}
	   	
	   	public double[] timestampIPs(boolean [] ipSelection){

	   		final Thread[] threads = newThreadArray(this.maxNumberOfThreads);
	   		final AtomicInteger ipIndexAtomic = new AtomicInteger(0);
	   		final int ipLength=this.resetURLs.length;
	   		final int debugLevel=this.debugLevel;
			final String []	metaURLs=this.metaURLs;
	   		final double [] ts= new double[this.cameraIPs.length];
	   		for (int i=0;i<ts.length;i++) ts[i]=0.0;

	   		
	   	//TODO: Multithread the next cycle (per-sensor)
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				public void run() {
	   					for (int ipIndex=ipIndexAtomic.getAndIncrement(); ipIndex<ipLength;ipIndex=ipIndexAtomic.getAndIncrement()){

	   						//			for (int ipIndex=0;ipIndex<this.resetURLs.length;ipIndex++) if ((ipIndex<ipSelection.length) && ipSelection[ipIndex]){
	   						try {
	   							Document dom=null;
	   							DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   							DocumentBuilder db = dbf.newDocumentBuilder();
	   							String url=metaURLs[ipIndex];
	   							if (debugLevel>2) System.out.println("timestampIPs:" + url );
	   							dom = db.parse(url);
	   							if (!dom.getDocumentElement().getNodeName().equals("meta")) {
	   								System.out.println("Root element: expected 'meta', got'" + dom.getDocumentElement().getNodeName()+"'");
	   								IJ.showMessage("Error","Root element: expected 'meta', got'" + dom.getDocumentElement().getNodeName()+"'"); 
	   								continue;
	   							}
	   							if (dom.getDocumentElement().getElementsByTagName("frame").getLength()==0) {
	   								String sError="reading timestamp failed. Do you have camera firmware version >=8.1.1.1?";
	   								System.out.println("ERROR: "+sError );
	   								IJ.showMessage("Error",sError); 
	   								continue;
	   							}
	   							ts[ipIndex]=Double.parseDouble(
	   									((Node) (((Node) dom.getDocumentElement().getElementsByTagName("timestamp").item(0)).getChildNodes().item(0))).getNodeValue());

	   						} catch(MalformedURLException e){
	   							System.out.println("Please check the URL:" + e.toString() );
	   							continue;
	   						} catch(IOException  e1){
	   							IJ.showStatus("");
	   							String error = e1.getMessage();
	   							if (error==null || error.equals(""))  error = ""+e1;
	   							IJ.showMessage("timestampIPs() ERROR", ""+error);
	   							continue;
	   						}catch(ParserConfigurationException pce) {
	   							pce.printStackTrace();
	   							continue;
	   						}catch(SAXException se) {
	   							se.printStackTrace(); 
	   							continue;
	   						}
	   					}
	   				}
	   			};
	   		}
	   		startAndJoin(threads);
	   		return ts;
	   	}

	   	public void resetCameras(boolean [] selection){
	   		if (this.debugLevel>2) {
	   			System.out.println("resetCameras(...)");
	   		    for (int ii=0;ii<this.cameraIPs.length;ii++)System.out.println(ii+":  "+this.cameraIPs[ii]);
	   		}
	   		
	   		resetIPs(selectIPs(selection));
	   		
	   	}

	   	public void resetIPs(boolean [] ipSelection){
	   		final Thread[] threads = newThreadArray(this.maxNumberOfThreads);
	   		final AtomicInteger ipIndexAtomic = new AtomicInteger(0);
	   		final int ipLength=this.resetURLs.length;
	   		final int debugLevel=this.debugLevel;
			final String []	resetURLs=this.resetURLs;
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				public void run() {
	   					for (int ipIndex=ipIndexAtomic.getAndIncrement(); ipIndex<ipLength;ipIndex=ipIndexAtomic.getAndIncrement()){
	   						String url="";
	   						try {
	   							DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   							DocumentBuilder db = dbf.newDocumentBuilder();
	   							url=resetURLs[ipIndex];
	   							if (debugLevel>2) System.out.println("--- resetURLs:" + url );
//	   							if (debugLevel>1) System.out.println("--- resetURLs:" + url );
	   							db.parse(url); // should be some XML (currently discarded)
	   						} catch(MalformedURLException e){
	   							System.out.println("Please check the URL:" + e.toString() );
	   							continue;
	   						} catch(IOException  e1){
	   							IJ.showStatus("");
	   							String error = e1.getMessage();
	   							if (error==null || error.equals(""))  error = ""+e1;
	   							IJ.showMessage("resetIPs() ERROR", url+"\n"+error);
	   							continue;
	   						}catch(ParserConfigurationException pce) {
	   							pce.printStackTrace();
	   							continue;
	   						}catch(SAXException se) {
	   							se.printStackTrace(); 
	   							continue;
	   						}
	   					}
	   				}

	   			};
	   		}
	   		startAndJoin(threads);
	   	}
	   	
	   	public ImagePlus [] getImages(final boolean [] acquire, boolean resetAndTrigger, final boolean show){
		   	final boolean [] acquireIPs=selectIPs(acquire);
			if (this.debugLevel>2) {
				System.out.println("getImages(...) 2");
			    for (int ii=0;ii<this.cameraIPs.length;ii++)System.out.println(ii+":  "+this.cameraIPs[ii]);
			}

		   	if (resetAndTrigger) {
			  resetCameras();
			  trigger();
		   	}  
		   	final double [] timestamps= timestampIPs(acquireIPs);
			if (this.debugLevel>2) System.out.println("getImages(): this.imagesIP.length=" + this.imagesIP.length);
	   		final Thread[] threads = newThreadArray(this.maxNumberOfThreads);
	   		final AtomicInteger ipIndexAtomic = new AtomicInteger(0);
	   		final int ipLength=this.resetURLs.length;
	   		final int debugLevel=this.debugLevel;
			final String []	imageURLs=this.imageURLs;
			final ImagePlus []	imagesIP=this.imagesIP;
			final JP46_Reader_camera [] jp4_Instances=this.jp4_Instances;
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						for (int ipIndex=ipIndexAtomic.getAndIncrement(); ipIndex<ipLength;ipIndex=ipIndexAtomic.getAndIncrement()){
							//							for (int i=0;i<this.imagesIP.length;i++){
							if (debugLevel>2) System.out.println("getImages()3: ipIndex="+ipIndex+" acquireIPs.length=" +acquireIPs.length+
									((ipIndex<acquireIPs.length)? (" acquireIPs[ipIndex]="+acquireIPs[ipIndex]):""));
							if ((ipIndex<acquireIPs.length) && acquireIPs[ipIndex]) {
								if (debugLevel>2) System.out.println("getImages:" + imageURLs[ipIndex] );
								if ((debugLevel>3) &&(imagesIP[ipIndex]!=null) && show){
									System.out.println("=============== old image ==========");
									jp4_Instances[ipIndex].listImageProperties(imagesIP[ipIndex]);

								}
								imagesIP[ipIndex]=jp4_Instances[ipIndex].openURL( // gains inside were OK
										imageURLs[ipIndex],
										"",
										true, //scale
										imagesIP[ipIndex],
										false); //show); // show image
								imagesIP[ipIndex].setProperty("timestamp", IJ.d2s(timestamps[ipIndex],6));
								imagesIP[ipIndex].setProperty("MIRRORED","NO");
								if (debugLevel>2) {
									jp4_Instances[ipIndex].listImageProperties(imagesIP[ipIndex],true); // to console - properties old - fixed
								}


								if (show){
									imagesIP[ipIndex].updateAndDraw(); /// Redisplays final image
									if (debugLevel>2){
										System.out.println("=============== new image ==========");
										jp4_Instances[ipIndex].listImageProperties(imagesIP[ipIndex]);
									}
								}
							}
						}
					}
				};
			}
			startAndJoin(threads);
			//TODO: Multithread the next cycle (per-sensor)
			final ImagePlus []	images=this.images;
			final int [][] channelMap = this.channelMap;
	   		final AtomicInteger imageIndexAtomic = new AtomicInteger(0);
	   		final int [] motorsPosition=this.motorsPosition;
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				public void run() {
	   					// Next calls do not require image acquisition, so default settings for the new instance are OK
						JP46_Reader_camera jp4_Instance= new JP46_Reader_camera(false);
	   					for (int imageIndex=imageIndexAtomic.getAndIncrement(); imageIndex<images.length;imageIndex=imageIndexAtomic.getAndIncrement())
	   						if ((imageIndex<acquire.length) && acquire[imageIndex]){
	   							//		   	for (int i=0;i<this.images.length;i++) if ((i<acquire.length) && acquire[i]) {
	   							int iIP=channelMap[imageIndex][0];
	   							if (sensorPresent[iIP]==null) { // system board for this channel did not respond null pointer - check cameras were detected
	   								images[imageIndex]=null;
	   								continue;
	   							}
	   							boolean singleSensor=true;
	   							for (int s=0;s<sensorPresent[iIP].length;s++) singleSensor&= !sensorPresent[iIP][s];
	   							if (singleSensor){ // no 10359 multiplexor
	   								images[imageIndex]=jp4_Instance.demuxClone(imagesIP[iIP]); 
	   							} else {
	   								int subCam=channelMap[imageIndex][1];
	   								if (!sensorPresent[iIP][subCam]){ // requested sensor does not exist
	   									images[imageIndex]=null;
	   									continue;
	   								}
	   								// skip missing channels, then demux
	   								for (int s=0;s<channelMap[imageIndex][1];s++) if (!sensorPresent[iIP][s]) subCam--;
	   								images[imageIndex]=jp4_Instance.demuxImage(imagesIP[iIP],subCam);
	   							}
	   							//		   		this.images[i]=this.jp4_Instances[j].demuxImageOrClone(this.imagesIP[j],this.channelMap[i][1]);
	   							if (images[imageIndex]==null) continue;
	   							if (debugLevel>2) jp4_Instance.listImageProperties(images[imageIndex]);
	   							if (debugLevel>2) jp4_Instance.encodeProperiesToInfo(images[imageIndex]);
	   							if (debugLevel>2) {
	   								jp4_Instance.listImageProperties(images[imageIndex]);
	   								jp4_Instance.decodeProperiesFromInfo(images[imageIndex]);
	   								jp4_Instance.listImageProperties(images[imageIndex]);
	   							}
	   							String title=IJ.d2s(timestamps[iIP],6).replace('.','_')+String.format("-%02d.tiff", imageIndex); // sensor number
	   							images[imageIndex].setTitle(title);

	   							if (flipImages[imageIndex]) { // this is only for physical mirror
	   								flipKeepBayer(images[imageIndex]);
	   								if (debugLevel>2){
	   									System.out.println("=============== flipKeepBayer ==========");
	   									jp4_Instance.listImageProperties(images[imageIndex]);
	   								}
	   								if (show) images[imageIndex].updateAndDraw();
	   							}
	   							images[imageIndex].setProperty("channel", String.format("%02d", imageIndex));
	   							images[imageIndex].setProperty("subcamera",""+getSubCamera(imageIndex)); 
	   							images[imageIndex].setProperty("sensor_port",""+getSensorPort(imageIndex)); 
	   							images[imageIndex].setProperty("subchannel", ""+getSubChannel(imageIndex));
	   							//		private int [] motorsPosition=      null; // motors steps when the images were acquired (for null)
	   							if (motorsPosition!=null) for (int m=0;m<motorsPosition.length;m++ ) {
	   								images[imageIndex].setProperty("MOTOR"+(m+1), ""+motorsPosition[m]);
	   							}
	   							jp4_Instance.encodeProperiesToInfo(images[imageIndex]);
	   						}
	   				}
	   			};
	   		}
	   		startAndJoin(threads);
	   		return this.imagesIP;
	   	}

	   	
	// should be already triggered!   	
   	
	   	public ImagePlus [] getImages(
	   			UVLEDandLasers uvLEDandLasers, // or null - not null only for adjustment machine lasers
	   			boolean [] pre_acquire,
	   			boolean [] pre_lasers,
	   			boolean resetAndTrigger,
	   			final boolean show){
//	   		this.reportTiming=this.debugLevel>1; // Need to be set by caller
	   		printTimingInit();
//			long startTime = System.nanoTime();
	   		final int debugLevel=this.debugLevel;
	   		final boolean opticalHeadLasersMode=(uvLEDandLasers!=null);
//	   		boolean [] lasers2={true,true};
	   		final boolean [] lasers=opticalHeadLasersMode?pre_acquire:pre_lasers;
	   		if (this.sensorPresent==null) { // cameras were not probed/configured
				probeCameraState(); // testing detection
		   		printTiming("=== probeCameraState()");
				setupCameraAcquisition();
		   		printTiming("=== setupCameraAcquisition()");
	   		}
	   		boolean [] acquire= new boolean[this.channelMap.length];
	   		for (int i=0;i<acquire.length;i++) {
	   			acquire[i]= ((i<pre_acquire.length) && pre_acquire[i]) || ((lasers!=null) && (i<lasers.length) && lasers[i]); 
	   		}
	   		boolean [] opticalHeadLasersInitial=null;
	   		int[] opticalHeadSequence={0,1,2};
	   		final boolean [][] headLaserWasOn={{false,true,false},{false,false,true}};
	   		if (opticalHeadLasersMode){
	   			opticalHeadLasersInitial=uvLEDandLasers.getLasers();
	   		} else if ((lasers!=null) && (this.laserPointers!=null)) {
	   			if (this.debugLevel>2) System.out.println("************ turning lasers off *************");
	   			this.laserPointers.setLasers(0);// turn off all lasers
	   			//lasersSequence
	   			if (this.debugLevel>2) {
	   				System.out.println("getImages(...) 1");
	   				for (int ii=0;ii<this.cameraIPs.length;ii++)System.out.println(ii+":  "+this.cameraIPs[ii]);
	   			}
	   		}

	   		if 	((lasers==null) && !opticalHeadLasersMode) {
		   		getImages(acquire, resetAndTrigger, show); // get images w/o lasers, flip if needed
	   			printTiming("=== Image acquisition");
	   			return this.imagesIP; 
	   		}
//scaleExposureForLasers
	   		// reduce exposure for lasers (both 4 and 2)
	   		final double scaleExposureForLasers=opticalHeadLasersMode?
	   				(((this.scaleExposureForHeadLasers>0.0) && (this.scaleExposureForHeadLasers<1.0))?this.scaleExposureForHeadLasers:0.0):
	   				(((this.scaleExposureForLasers>0.0) && (this.scaleExposureForLasers<1.0))?this.scaleExposureForLasers:0.0);
	   		if (scaleExposureForLasers>0) setupCameraAcquisition(scaleExposureForLasers);
	   		final boolean [] lasersIPs= selectIPs(opticalHeadLasersMode?acquire:lasers);
	   		final int [] sequence=opticalHeadLasersMode?opticalHeadSequence:this.laserPointers.getSequence();
	   		if (debugLevel>2)	for (int i=0;i<sequence.length;i++) System.out.println(String.format("Laser sequence[%d]=0x%x", i,sequence[i]));
	   		final ImagePlus [][] laserImagesIP=new ImagePlus[sequence.length-1][this.imagesIP.length];
	   		if (this.debugLevel>2) System.out.println("this.imagesIP.length="+this.imagesIP.length+" sequence.length="+sequence.length);
	   		final String [] imageURLs=this.imageURLs;
			final JP46_Reader_camera [] jp4_Instances=this.jp4_Instances;
	   		final Thread[] threads = newThreadArray(this.maxNumberOfThreads);
	   		final AtomicInteger ipIndexAtomic = new AtomicInteger(0);
	   		for (int nSeqNum=1; nSeqNum<sequence.length;nSeqNum++){
	   			if (opticalHeadLasersMode) {
	   				boolean [] bLasers={(sequence[nSeqNum]&1)!=0,(sequence[nSeqNum]&2)!=0};
	   				uvLEDandLasers.setLasersAndUV(bLasers,null,null);
			   		if (debugLevel>2)	System.out.println("uvLEDandLasers.setLasersAndUV("+bLasers[0]+","+bLasers[1]+")");
	   			} else {
	   				this.laserPointers.setLasers (sequence[nSeqNum]);// turn on selected laser
			   		if (debugLevel>2)	System.out.println(String.format("this.laserPointers.setLasers (0x%x)", sequence[nSeqNum]));
	   			}
	   			resetIPs(selectIPs(lasersIPs)); // flush buffer
	   			trigger(); // trigger cameras
	   			ipIndexAtomic.set(0);
	   			final int fnSeqNum=nSeqNum;
	   			final boolean reportTiming=this.reportTiming;
	   			for (int ithread = 0; ithread < threads.length; ithread++) {
	   				threads[ithread] = new Thread() {
	   					public void run() {
	   						for (int ipIndex=ipIndexAtomic.getAndIncrement(); ipIndex<lasersIPs.length;ipIndex=ipIndexAtomic.getAndIncrement()){
	   							long st=System.nanoTime();
	   							if (debugLevel>1) {
	   								System.out.println("image url["+ipIndex+"]="+imageURLs[ipIndex]);
	   							}
	   							laserImagesIP[fnSeqNum-1][ipIndex]=jp4_Instances[ipIndex].openURL(
	   									imageURLs[ipIndex],
	   									"",
	   									true, //scale
	   									null, // new image, nothing to reuse
	   									false); // show image no timestamps here
	   							if (reportTiming){
	   								System.out.println("==== "+imageURLs[ipIndex]+ " opened in "+IJ.d2s(0.000000001*(System.nanoTime()-st),3)+" sec");
	   							}
	   							laserImagesIP[fnSeqNum-1][ipIndex].setProperty("MIRRORED","NO");
	   						}
	   					}
	   				};
	   			}
	   			startAndJoin(threads);
	   			printTiming("=== Image # "+nSeqNum+" acquisition");

	   		}
// acquire no-laser image after all others are acquired, so the camera is completely still at that time
   			if (opticalHeadLasersMode) {
   				boolean [] bLasers={(sequence[0]&1)!=0,(sequence[0]&2)!=0};
   				uvLEDandLasers.setLasersAndUV(bLasers,null,null);
		   		if (debugLevel>2)	System.out.println("uvLEDandLasers.setLasersAndUV("+bLasers[0]+","+bLasers[1]+")");
   			} else {
   				this.laserPointers.setLasers (sequence[0]);// turn on selected laser
		   		if (debugLevel>2)	System.out.println(String.format("this.laserPointers.setLasers (0x%x)", sequence[0]));
   			}
//	   		this.laserPointers.setLasers (0);// turn off all lasers
	   		if (scaleExposureForLasers>0.0) setupCameraAcquisition(1.0);
	   		
	   		getImages(acquire, resetAndTrigger, show); // get images w/o lasers, flip if needed
   			printTiming("=== Final (no-laser) image # 0 acquired");

   			if (opticalHeadLasersMode) {
   				uvLEDandLasers.setLasersAndUV(opticalHeadLasersInitial,null,null);
		   		if (debugLevel>2)	System.out.println("Restoring lasers: uvLEDandLasers.setLasersAndUV("+opticalHeadLasersInitial[0]+","+opticalHeadLasersInitial[1]+")");
   			}

	   		
	   		// now demux images (if composite) and process laser spots
	   		//		   	for (int sensorNum=0;sensorNum<this.images.length;sensorNum++) if (this.images[sensorNum]!=null){
//	   		final MatchSimulatedPattern matchSimulatedPattern = new MatchSimulatedPattern();
//	   		final ImagePlus imp_pointed=null;
	   		final int bayerG1=0;
	   		final int bayerG2=3;
	   		final int bayerR=1;
	   		final int bayerB=2;
	   		final boolean useOther=laserPointers.laserPointer.useOther;
	   		final boolean otherGreen=laserPointers.laserPointer.otherGreen;

	   		final int [][] channelMap=this.channelMap;
	   		final ImagePlus [] images=this.images;
			final boolean [] flipImages=this.flipImages;
			final LaserPointers laserPointers=this.laserPointers;
			final boolean [][] sensorPresent=this.sensorPresent;
			final int debugSensorNumber=this.debugSensorNumber;
   			printTiming("=== Acquisition done, starting multi-threaded laser pointer location processing");

//TODO: Multithread the next cycle (per-sensor)
			final AtomicInteger sensorNumAtomic = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
	   					// Next calls do not require image acquisition, so default settings for the new instance are OK
						JP46_Reader_camera jp4_Instance= new JP46_Reader_camera(false);
						MatchSimulatedPattern matchSimulatedPattern = new MatchSimulatedPattern();
						for (int sensorNum=sensorNumAtomic.getAndIncrement(); sensorNum<lasers.length;sensorNum=sensorNumAtomic.getAndIncrement()) // null pointer
							if (lasers[sensorNum] && (images[sensorNum]!=null)){
								//	   		for (int sensorNum=0;sensorNum<lasers.length;sensorNum++)  if (lasers[sensorNum] && (this.images[sensorNum]!=null)){ // lasers - here sensors to use lasers for
								int iIP=channelMap[sensorNum][0];
								double saturationRed=255.0;
								if (images[sensorNum].getProperty("saturation_0")!=null) saturationRed=Double.parseDouble((String)images[sensorNum].getProperty("saturation_0"));
								if (scaleExposureForLasers>0) saturationRed*=scaleExposureForLasers; // scaled to reduced exposure time  
								double[][] backgroundBayer =matchSimulatedPattern.splitBayer (images[sensorNum],  null, true); // full window
								double[][] pointedBayer;
								double [][] pontersXY=null;
								//		   			int iIP=this.channelMap[sensorNum][0];
								int len=backgroundBayer[bayerG1].length;
								double [][] greens=useOther?(new double [sequence.length][]):null;
								double [][] reds=  new double [sequence.length][];
								if (useOther){
									greens[0]=otherGreen?(new double[len]):backgroundBayer[bayerB];	
									if (otherGreen) for (int j=0; j<len;j++) greens[0][j]=0.5*(backgroundBayer[bayerG1][j]+backgroundBayer[bayerG2][j]);
								}
								reds[0]=backgroundBayer[bayerR]; // no need to clone?
								if (debugLevel>2)  System.out.println("getImages (): sensorNum="+sensorNum+", iIP="+iIP+" sequence.length="+sequence.length);
//								if ((debugLevel>1) && (sensorNum==13))  System.out.println("getImages (): sensorNum="+sensorNum+", iIP="+iIP+" sequence.length="+sequence.length);
								boolean singleSensor=true;
								for (int s=0;s<sensorPresent[iIP].length;s++) singleSensor&= !sensorPresent[iIP][s];
								ImagePlus imp_pointed=null;
								for (int nSeqNum=1; nSeqNum<sequence.length; nSeqNum++){
									if (debugLevel>2) System.out.println("getImages(): sensorNum="+sensorNum+", iIP="+iIP+" nSeqNum="+nSeqNum+" this.channelMap["+sensorNum+"][1]="+channelMap[sensorNum][1]);
									// here we are processing only same sensors that produced no-laser images, no need to check if they exist again
									if (singleSensor){ // no 10359 multiplexor
										//imp_pointed=jp4_Instances[iIP].demuxClone(laserImagesIP[nSeqNum-1][iIP]); 
										imp_pointed=jp4_Instance.demuxClone(laserImagesIP[nSeqNum-1][iIP]); 
									} else {
										int subCam=channelMap[sensorNum][1];
										// skip missing channels, then demux
										for (int s=0;s<channelMap[sensorNum][1];s++) if (!sensorPresent[iIP][s]) subCam--;
										//imp_pointed=jp4_Instances[iIP].demuxImage(laserImagesIP[nSeqNum-1][iIP],subCam);
										if (laserImagesIP[nSeqNum-1][iIP]==null){
											System.out.println("getImages(): laserImagesIP["+(nSeqNum-1)+"]["+iIP+"]==null");
										}
										imp_pointed=jp4_Instance.demuxImage(laserImagesIP[nSeqNum-1][iIP],subCam);
									}
									imp_pointed.setTitle(imp_pointed.getTitle()+"_"+sensorNum);
									if (flipImages[sensorNum]) {
										flipKeepBayer(imp_pointed);
										if (show) imp_pointed.updateAndDraw();
									}
									pointedBayer=matchSimulatedPattern.splitBayer (imp_pointed, null, true);
									if (useOther){
										greens[nSeqNum]=otherGreen?(new double[len]):pointedBayer[bayerB];	
										if (otherGreen) for (int j=0; j<len;j++) greens[nSeqNum][j]=0.5*(pointedBayer[bayerG1][j]+pointedBayer[bayerG2][j]);
									}
									reds[nSeqNum]=pointedBayer[bayerR].clone();
								}
								
//minimalIntensity		scaleExposureForLasers>0	saturationRed						
								pontersXY=laserPointers.laserPointer.getPointerXY( // returns x,y pair or null if pointer not detected
										opticalHeadLasersMode,
										saturationRed,
										scaleExposureForLasers, //>0.0,
										greens,        // combined Bayer greens for each image, starting with no-laser
										reds,          // red Bayer component for each image, starting with no-laser
										opticalHeadLasersMode?headLaserWasOn:laserPointers.laserWasOn(), // array specifying which image should have pointer on, for each laser
										imp_pointed.getWidth(),// image width in pixels
										imp_pointed.getTitle(),             // String title,
										debugLevel+  ((sensorNum==debugSensorNumber)?1:0)           // debug level (normal == 1)
//										debugLevel
								);
								int pointersDetected=0;
								if (!opticalHeadLasersMode){
									for (int nPointer=0; nPointer<laserPointers.getNumberOfLasers();nPointer++){
										if (pontersXY[nPointer]!=null){
											if (debugLevel>1) System.out.println("image:"+sensorNum+" pointer #"+(nPointer+1)+
													", X="+pontersXY[nPointer][0]+", Y="+pontersXY[nPointer][1]);
											images[sensorNum].setProperty("POINTER_X_"+nPointer, IJ.d2s(pontersXY[nPointer][0],1));
											images[sensorNum].setProperty("POINTER_Y_"+nPointer, IJ.d2s(pontersXY[nPointer][1],1));
											pointersDetected++;
										}

									}
									if ((pointersDetected>0) && (debugLevel>0)) {
										System.out.println("image:"+sensorNum+" - "+pointersDetected+" pointer"+((pointersDetected>1)?"s":"")+" detected");
									}									
								} else {
									for (int nPointer=0; nPointer<2;nPointer++){
										if (pontersXY[nPointer]!=null){
											if (debugLevel>1) System.out.println("image:"+sensorNum+" pointer #"+(nPointer+1)+
													", X="+pontersXY[nPointer][0]+", Y="+pontersXY[nPointer][1]);
											images[sensorNum].setProperty("HEAD_POINTER_X_"+nPointer, IJ.d2s(pontersXY[nPointer][0],1));
											images[sensorNum].setProperty("HEAD_POINTER_Y_"+nPointer, IJ.d2s(pontersXY[nPointer][1],1));
											pointersDetected++;
										}
									}
								}
								images[sensorNum].setProperty("HEAD_POINTERS", pointersDetected+"");
								//jp4_Instances[iIP].encodeProperiesToInfo(images[sensorNum]);		   		
								jp4_Instance.encodeProperiesToInfo(images[sensorNum]);		   		
							}
					}
				};
			}
			startAndJoin(threads);

//	   		this.laserPointers.setLasers (0);// turn off all lasers
   			printTiming("=== Image acquisition/laser pointers location");
	   		return this.imagesIP;
	   	}
	   	
	   	public double [][] getHeadPointers(ImagePlus imp){
	   		if (imp.getProperty("HEAD_POINTERS")==null) return null;
	   		int numHeadPointers=Integer.parseInt((String) imp.getProperty("HEAD_POINTERS")); 
	   		double [][] headPointers=new double[numHeadPointers][];
	   		for (int n=0;n<headPointers.length;n++){
	   			if ((imp.getProperty("HEAD_POINTER_X_"+n)!=null) && (imp.getProperty("HEAD_POINTER_Y_"+n)!=null)){
	   				headPointers[n]=new double[2];
	   				headPointers[n][0]=Double.parseDouble((String)imp.getProperty("HEAD_POINTER_X_"+n));
	   				headPointers[n][1]=Double.parseDouble((String)imp.getProperty("HEAD_POINTER_Y_"+n));
	   			} else headPointers[n]=null;
	   		}
	   		return headPointers;
	   	}
	   	

	   	
	   	
	   	//	   	return array of last acquired images
	   	public ImagePlus [] getImages(){
	   		return this.images;
	   	}
	   	public ImagePlus [] getImages(int numberOfPointers){ // return images with at least this number of pointers detected
	   		ImagePlus [] filteredImages=this.images.clone();
	   		for (int i=0;i<filteredImages.length;i++) if (getNumberOfPointers(i)<numberOfPointers) filteredImages[i]=null;
	   		return filteredImages;
	   	}

	   	public int getNumberOfPointers (int sensorNum){
	   		return getNumberOfPointers (this.images[sensorNum]);
	   	}

	   	public int getNumberOfPointers (ImagePlus imp){
	   		if (imp==null) return 0;
	   		if (imp.getProperty("POINTERS")!=null) return Integer.parseInt((String) imp.getProperty("POINTERS"));
	   		else return 0;
	   	}

	   	public int saveImages(boolean [] selection, String directory, boolean updateStatus){
	   		int numImg=0;
	   		boolean notYetPrinted=true;
	   		for (int i=0;(i<this.channelMap.length) && ((selection==null) || (i<selection.length));i++)
	   			if ((this.images[i]!=null) && ((selection==null) || selection[i])){
	   				FileSaver fs=new FileSaver(this.images[i]);
	   				String path=directory+Prefs.getFileSeparator()+this.images[i].getTitle();
	   				if (updateStatus) IJ.showStatus("Saving "+path);
	   				if (this.debugLevel>0){
	   					if (this.debugLevel>1) System.out.println("Saving "+path); // was >0
	   					else if (notYetPrinted) {
	   						System.out.println("Saving "+path+ " (and other with the same timestamp)" ); // was >0
	   						notYetPrinted=false;
	   					}
	   				}
	   				fs.saveAsTiff(path);
	   				numImg++;
	   			}
	   		return numImg;
	   	}
/*
								imp_psf = new ImagePlus(filePaths[dirNum][fileNum][1], stack);
								//							if (DEBUG_LEVEL>1) imp_psf.show();
								if (DEBUG_LEVEL>1) System.out.println("Saving result to"+filePaths[dirNum][fileNum][1]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
								FileSaver fs=new FileSaver(imp_psf);
								fs.saveAsTiffStack(filePaths[dirNum][fileNum][1]);
	   	
 */
	   	
	   	/**
	   	 * FLips image vertically and crops by one pixel at the top and bottom (to preserve Bayer pattern)
	   	 * @param imp
	   	 * @return
	   	 */
	   	public ImagePlus flipAndCrop(ImagePlus imp){
	   		String title=imp.getTitle();
	   		imp.setTitle(title+"-original");
	   		ImageProcessor ip=imp.getProcessor();
	   		ip.setRoi(0,1,ip.getWidth(),ip.getHeight()-2);
	   		ImageProcessor ip2=ip.crop();
	   		ip2.flipVertical();
	   		ImagePlus imp_new=new ImagePlus(title,ip2);
// copy all properties and add "MIRRORED=Y"	
			Set<Object> imp_set;
			Properties imp_prop;
			Iterator<Object> itr;
			String str;
			imp_prop=imp.getProperties();
			if (imp_prop!=null) {
				imp_set=imp_prop.keySet();
				itr=imp_set.iterator();
				while(itr.hasNext()) {
					str = (String) itr.next();
					imp_new.setProperty(str,imp_prop.getProperty(str));
				}
			}
			imp_new.setProperty("MIRRORED","YES");
	   		return imp_new;
	   	}
	   	public ImagePlus flipKeepBayer(ImagePlus imp){
	   		ImageProcessor ip=imp.getProcessor();
	   		float [] pixels=(float[]) ip.getPixels();
	   		float [] new_pixels=new float[pixels.length];
	   		int width=ip.getWidth();
	   		int height=ip.getHeight();
	   		for (int i=1;i<height-1;i++) for (int j=0;j<width;j++){
	   			new_pixels[i*width+j]=pixels[(height-i)*width+j];
	   		}
	   	    for (int j=0;j<width;j++) new_pixels[j]=new_pixels[2*width+j]; // copy second line to 0
	   	    ip.setPixels(new_pixels);
			imp.setProperty("MIRRORED","YES");
	   		return imp;
	   	}
	   	public ImagePlus flipKeepBayerNew(ImagePlus imp){
	   		float [] pixels=(float []) imp.getProcessor().getPixels();
	   		float [] new_pixels=new float[pixels.length];
	   		int width=imp.getWidth();
	   		int height=imp.getHeight();
	   		for (int i=1;i<height-1;i++) for (int j=0;j<width;j++){
	   			new_pixels[i*width+j]=pixels[(height-i)*width+j];
	   		}
	   		for (int j=0;j<width;j++) new_pixels[j]=new_pixels[2*width+j]; // copy second line to 0
	   		ImageProcessor ip=new FloatProcessor(imp.getWidth(), imp.getHeight());
	   		ip.setPixels(new_pixels);
	   		ip.resetMinAndMax();
	   		ImagePlus imp_flipped=  new ImagePlus(imp.getTitle()+"-flipped", ip);
	   		imp_flipped.setProperty("MIRRORED","YES");
	   		return imp_flipped;
	   	}

// in detection mode - reduce exp	   	
	   	
	   	public void test1(boolean useLasers) {
	   		getImages(
	   				null, // UVLEDLasers
	   				selectAllSubcameras(),
	   				(useLasers?selectAllSubcameras():null),
	   				true,
	   				this.debugLevel>1); // reset and trigger
			if (this.debugLevel>1) {
				if (this.debugLevel>2) System.out.println("++++++++++++++++++ Image Properies ++++++++++++++++++++++++++++");
				for (int i=0;i<this.images.length;i++) if (this.images[i]!=null) {
					int j=this.channelMap[i][0];
					if (this.debugLevel>2) {
						System.out.println("Image #"+i);
						this.jp4_Instances[j].listImageProperties(this.images[i]);
					}
					if ((this.debugLevel>2) || (getNumberOfPointers (i)>0)) {
						this.images[i].show();
						this.images[i].updateAndDraw();
					}
				}
			}
	   	}
	   	
	   	public void acquire(String directory, boolean useLasers, boolean updateStatus) {
			getImages(
	   				null, // UVLEDLasers
					selectAllSubcameras(),
					(useLasers?selectAllSubcameras():null),
					true,
					this.debugLevel>1); // reset and trigger
		   	int numImg=saveImages(selectAllSubcameras(), directory, updateStatus);
		   	if (this.debugLevel>1) System.out.println("Saved "+numImg+" images.");
	   	}
	    public ImagePlus acquireSingleImage (boolean useLasers, boolean updateStatus){
			getImages(
	   				null, // UVLEDLasers
					selectAllSubcameras(),
					(useLasers?selectAllSubcameras():null),
					true, 
					this.debugLevel>1); // reset and trigger
	    	this.lastTimestamp=(String) this.images[0].getProperty("timestamp");
	    	return this.images[0];
	    }
	    public ImagePlus [] acquireSeveralImages (boolean useLasers, boolean updateStatus){
			getImages(
	   				null, // UVLEDLasers
					selectAllSubcameras(),
					(useLasers?selectAllSubcameras():null),
					true, 
					this.debugLevel>1); // reset and trigger
	    	this.lastTimestamp=(String) this.images[0].getProperty("timestamp");
	    	return this.images;
	    }

	    public ImagePlus acquireSingleImage (UVLEDandLasers uvLEDLasers, boolean updateStatus){
			getImages(
					uvLEDLasers, // UVLEDLasers
					selectAllSubcameras(),
					null,
					true, 
					this.debugLevel>1); // reset and trigger
	    	this.lastTimestamp=(String) this.images[0].getProperty("timestamp");
	    	return this.images[0];
	    }

	    
	    public String getLastTimestampUnderscored(){
	    	return this.lastTimestamp.replace('.','_');
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
	
	public static class PowerControl{
		public boolean [] states={false,false,false,false,false};
		public boolean [] useControl={true,true,false,true,true};
//		public int [] lightsChannels={2,3,4}; // which lights to control on on/off
		public int [] lightsChannels={3,4}; // which lights to control on on/off
		public String [] groups={"heater","fan","light","light1","light2"};
    	public int debugLevel=1;
    	private String powerIP="192.168.0.80";
    	private double lightsDelay=1.0;
    	private final String urlFormat="http://%s/insteon/index.php?cmd=%s&group=%s&timestamp=%d";
    	private final String rootElement="Document";
    	public boolean powerConrtolEnabled=false;
    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"powerIP",this.powerIP+"");
    		properties.setProperty(prefix+"powerConrtolEnabled",this.powerConrtolEnabled+"");
    		properties.setProperty(prefix+"lightsDelay",this.lightsDelay+"");
    	}
    	//Integer.decode(string)
    	public void getProperties(String prefix,Properties properties){
    		if (properties.getProperty(prefix+"powerIP")!=null)
    			this.powerIP=properties.getProperty(prefix+"powerIP");
    		if (properties.getProperty(prefix+"powerConrtolEnabled")!=null)
    			this.powerConrtolEnabled=Boolean.parseBoolean(properties.getProperty(prefix+"powerConrtolEnabled"));
    		if (properties.getProperty(prefix+"lightsDelay")!=null)
    			this.lightsDelay=Double.parseDouble(properties.getProperty(prefix+"lightsDelay"));
    	}
    	
    	public boolean setPower (String group, String state){
    		if (!powerConrtolEnabled) {
    			System.out.println("=== Power control is disabled ===");
    			return false;
    		}
    		long thisTime=System.nanoTime();
    		if (debugLevel>0) System.out.println("=== Power control: "+group+":"+state+" ===");
    		String url=String.format(urlFormat,this.powerIP,state,group,thisTime);    	
    		if (this.debugLevel>1) System.out.println("setPower: "+url); 
    		Document dom=null;
    		try {
    			DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
    			DocumentBuilder db = dbf.newDocumentBuilder();
    			dom = db.parse(url);
    			if (!dom.getDocumentElement().getNodeName().equals(rootElement)) {
    				System.out.println("Root element: expected \""+rootElement+"\", got \"" + dom.getDocumentElement().getNodeName()+"\"");
    				IJ.showMessage("Error","Root element: expected \""+rootElement+"\", got \"" + dom.getDocumentElement().getNodeName()+"\""); 
    				return false;
    			}
    			//    				boolean responceError= (dom.getDocumentElement().getElementsByTagName("error").getLength()!=0);
    			//    				if (responceError) {
    			//    					System.out.println("ERROR: register write ("+url+") FAILED" );
    			//    					IJ.showMessage("Error","register write ("+url+") FAILED"); 
    			//    					return false;  
    			//    				}
//				System.out.println(dom);
//				System.out.println(dom.getDocumentElement().toString());
//				System.out.println(dom.getDocumentElement().getNodeName());
//				System.out.println(dom.getDocumentElement().getElementsByTagName("message"));
//				System.out.println(dom.getDocumentElement().getElementsByTagName("state"));
				if (dom.getDocumentElement().getElementsByTagName("state").getLength()>0){
					if (debugLevel>0) System.out.println("state="+((Node) (((Node) dom.getDocumentElement().getElementsByTagName("state").item(0)).getChildNodes().item(0))).getNodeValue());
				} else {
					System.out.println("*** Empty Document from Insteon");
					return false;
				}
    		} catch(MalformedURLException e){
    			System.out.println("Please check the URL:" + e.toString() );
    			return false;
    		} catch(IOException  e1){
    			IJ.showStatus("");
    			String error = e1.getMessage();
    			if (error==null || error.equals(""))  error = ""+e1;
    			IJ.showMessage("setPower ERROR", ""+error);
    			return false;
    		}catch(ParserConfigurationException pce) {
    			pce.printStackTrace();
    			return false;
    		}catch(SAXException se) {
    			se.printStackTrace(); 
    			return false;
    		}
    		if (debugLevel>1) {
    			// debugging
    			TransformerFactory tf = TransformerFactory.newInstance();
    			Transformer transformer=null;
    			try {
    				transformer = tf.newTransformer();
    			} catch (TransformerConfigurationException e) {
    				// TODO Auto-generated catch block
    				e.printStackTrace();
    			}
    			transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
    			StringWriter writer = new StringWriter();
    			try {
    				transformer.transform(new DOMSource(dom), new StreamResult(writer));
    			} catch (TransformerException e) {
    				// TODO Auto-generated catch block
    				e.printStackTrace();
    			}
    			String output = writer.getBuffer().toString().replaceAll("\n|\r", "");
    			System.out.println(output);
    		}
    		
    		if (this.debugLevel>1) System.out.println("=== Power control: OK ===");
    		for (int i=0;i<this.groups.length;i++) if (this.groups[i].equals(group)){
    			this.states[i]=state.equals("on");
    		}
    		return true;
    	}
    	
    	public boolean showDialog(String title, boolean control) {
    		GenericDialog gd = new GenericDialog(title);
    		boolean heaterOn=false, fanOn=false, lightOn=false, light1On=false, light2On=false;
			gd.addCheckbox("Enable power control (heater, fan, lights) ", this.powerConrtolEnabled);
    		gd.addStringField("IP address of the power control",this.powerIP,15);
			gd.addNumericField("Delay after lights on", this.lightsDelay,  1,4,"sec");
  		    		
    		if (control){
    			if(useControl[0]) gd.addCheckbox("Heater On", heaterOn);
    			if(useControl[1]) gd.addCheckbox("Fan On", fanOn);
    			if(useControl[2]) gd.addCheckbox("Lights On", lightOn);
    			if(useControl[3]) gd.addCheckbox("Lights Top On", light1On);
    			if(useControl[4]) gd.addCheckbox("Lights Bottom On", light2On);
    		}
    	    WindowTools.addScrollBars(gd);
			if (control) gd.enableYesNoCancel("OK", "Control Power");
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
    	    this.powerConrtolEnabled=gd.getNextBoolean();
    	    this.powerIP=gd.getNextString();
			this.lightsDelay=gd.getNextNumber();
    	    if (control){
    	    	if(useControl[0]) heaterOn=gd.getNextBoolean();
    	    	if(useControl[1]) fanOn=gd.getNextBoolean();
    	    	if(useControl[2]) lightOn=gd.getNextBoolean();
    	    	if(useControl[3]) light1On=gd.getNextBoolean();
    	    	if(useControl[4]) light2On=gd.getNextBoolean();
    	    	if (!gd.wasOKed()) {
    	    		if(useControl[0])setPower("heater",heaterOn?"on":"off"); // setPower("heater","state");
    	    		if(useControl[1])setPower("fan",fanOn?"on":"off");       // setPower("fan","state");
    	    		if(useControl[2])setPower("light",lightOn?"on":"off");   // setPower("light","state");
    	    		if(useControl[3])setPower("light1",light1On?"on":"off"); // setPower("light1","state");
    	    		if(useControl[4])setPower("light2",light2On?"on":"off"); //  setPower("light2","state");
    	    	}
    	    }
            return true;
    	}
    	public void lightsOnWithDelay(){
    		// turn on only new
    		boolean allOn=true;
    		for (int chn:this.lightsChannels) allOn&=this.states[chn];
    		if (allOn  || !this.powerConrtolEnabled) return; // already on
    		for (int chn:this.lightsChannels) if (!this.states[chn]) setPower(this.groups[chn],"on");
			System.out.print("Sleeping "+this.lightsDelay+" seconds to let lights stabilize on...");
    		try {
				TimeUnit.MILLISECONDS.sleep((long) (1000*this.lightsDelay));
			} catch (InterruptedException e) {
				System.out.println("Sleep was interrupted");
				// TODO Auto-generated catch block
			}
			System.out.println(" Done");

    	}
    	
    	public void lightsOff(){
    		if (!this.powerConrtolEnabled) return; // already on
    		for (int chn:this.lightsChannels) if (this.states[chn]) setPower(this.groups[chn],"off");
    	}
    	
    	public boolean isPowerControlEnabled(){
    		return this.powerConrtolEnabled;
    	}
		public void setDebugLevel(int debugLevel){
			this.debugLevel=debugLevel;
		}
	}
	
    public static class LaserPointers{
    	public int debugLevel=1;
    	public MatchSimulatedPattern.LaserPointer laserPointer=null; // class that defines geometry and extracts laser coordinates
    	private int state=-1; // bit-wise on/off, <0 - undefined
    	private String lasersIP="192.168.0.13";
    	private String lasersSequence="3c5a96"; // laser combinations to try, (non including 0)
    	private int [][][] registerData=
    	{
    			{{0x2600,0xff},{0x2700,0xff}}, // - - - - 0
    			{{0x2600,0xdf},{0x2700,0xff}}, // - - - + 1
    			{{0x2600,0xff},{0x2700,0xef}}, // - - + - 2
    			{{0x2600,0xdf},{0x2700,0xef}}, // - - + + 3
    			{{0x2600,0xef},{0x2700,0xff}}, // - + - - 4
    			{{0x2600,0xcf},{0x2700,0xff}}, // - + - + 5
    			{{0x2600,0xef},{0x2700,0xef}}, // - + + - 6
    			{{0x2600,0xcf},{0x2700,0xef}}, // - + + + 7
    			{{0x2600,0xff},{0x2700,0xdf}}, // + - - - 8
    			{{0x2600,0xdf},{0x2700,0xdf}}, // + - - + 9
    			{{0x2600,0xff},{0x2700,0xcf}}, // + - + - 10
    			{{0x2600,0xdf},{0x2700,0xcf}}, // + - + + 11
    			{{0x2600,0xef},{0x2700,0xdf}}, // + + - - 12
    			{{0x2600,0xcf},{0x2700,0xdf}}, // + + - + 13
    			{{0x2600,0xef},{0x2700,0xcf}}, // + + + - 14
    			{{0x2600,0xcf},{0x2700,0xcf}}  // + + + + 15
    	};
    	public LaserPointers(MatchSimulatedPattern.LaserPointer laserPointer){
    		this.laserPointer=laserPointer;
    	}
    	public int getNumberOfLasers(){
    		int numPointers=0;
			for (numPointers=0; (1<<numPointers)<this.registerData.length;numPointers++);
			if (this.debugLevel>3) System.out.println("getNumberOfLasers() first -> "+numPointers);

    		if (this.laserPointer==null) return numPointers;

			if (numPointers>this.laserPointer.getNumberOfPointers()) numPointers=this.laserPointer.getNumberOfPointers();
			if (this.debugLevel>3) System.out.println("getNumberOfLasers() second -> "+numPointers);
    		return numPointers;
    	}
    	//lasersSequence
    	
    	public int [] getSequence(){
    		int [] sequence = new int [this.lasersSequence.length()+1];
    		sequence[0]=0;
    		for (int i=0;i<this.lasersSequence.length();i++){
    			sequence[i+1]=Integer.parseInt(this.lasersSequence.substring(i, i+1),16);
    		}
    		return sequence;
    	}
    	/**
    	 * 
    	 * @param numLaser number of laser (0..1)
    	 * @return // array, starting with lasers-off image, showing if this laser was on on the corresponding image
    	 */
    	public boolean [] laserWasOn(int numLaser){
    		int [] sequence=getSequence();
    		boolean [] wasOn=new boolean [sequence.length];
    		for (int i=0;i<sequence.length;i++) wasOn[i]= ((sequence[i]>>numLaser) & 1)!=0; 
    		return wasOn;
    	}
    	
//getNumberOfLasers()    	
    	public boolean [][] laserWasOn(){
    		boolean [][] wasOn=new boolean [getNumberOfLasers()][];
    		for (int i=0;i<wasOn.length;i++) wasOn[i]= laserWasOn(i); 
    		return wasOn;
    	}
    	
    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"lasersIP",this.lasersIP+"");
    		properties.setProperty(prefix+"lasersSequence",this.lasersSequence+"");
    		properties.setProperty(prefix+"registerData.length",this.registerData.length+"");
    		for (int i=0;i<registerData.length;i++) {
        		properties.setProperty(prefix+"registerData_"+i+".length",this.registerData[i].length+"");
    			for (int j=0;j<registerData[i].length;j++) {
            		properties.setProperty(prefix+"registerData_"+i+"_"+j+".address","0x"+Integer.toHexString(this.registerData[i][j][0]));
            		properties.setProperty(prefix+"registerData_"+i+"_"+j+".data",   "0x"+Integer.toHexString(this.registerData[i][j][1]));
    			}
    		}
    	}
    	//Integer.decode(string)
    	public void getProperties(String prefix,Properties properties){
    		if (properties.getProperty(prefix+"lasersIP")!=null)
    			this.lasersIP=properties.getProperty(prefix+"lasersIP");
    		if (properties.getProperty(prefix+"lasersSequence")!=null)
    			this.lasersSequence=properties.getProperty(prefix+"lasersSequence");
    		if (properties.getProperty(prefix+"registerData.length")!=null) {
    			this.registerData=new int [Integer.parseInt(properties.getProperty(prefix+"registerData.length"))][][];
        		for (int i=0;i<registerData.length;i++) {
        			this.registerData[i]=new int [Integer.parseInt(properties.getProperty(prefix+"registerData_"+i+".length"))][2];
            		if (properties.getProperty(prefix+"registerData_"+i+".length")!=null) {
            			this.registerData[i]=new int [Integer.parseInt(properties.getProperty(prefix+"registerData_"+i+".length"))][2];
            			for (int j=0;j<registerData[i].length;j++) {
            				// or make it throw exceptions if some data is missing?
                    		if (properties.getProperty(prefix+"registerData_"+i+"_"+j+".address")!=null)
                    			this.registerData[i][j][0]=Integer.decode(properties.getProperty(prefix+"registerData_"+i+"_"+j+".address"));
                    		if (properties.getProperty(prefix+"registerData_"+i+"_"+j+".data")!=null)
                    			this.registerData[i][j][1]=Integer.decode(properties.getProperty(prefix+"registerData_"+i+"_"+j+".data"));
            			}
            		}
        		}    			
    		}
    	}
    	/**
    	 * 
    	 * @param title - Dialog title
    	 * @param askRegenerate (will not change number of laser pointers and registers to change if false)
    	 * @return true - OK, false - canceled
    	 */
    	public boolean showDialog(String title, boolean askRegenerate) {
    		GenericDialog gd = new GenericDialog(title);
    		gd.addStringField("IP address of the laser pointers",this.lasersIP,15);
    		gd.addStringField("Combinations of laser to try (excluding 0)",this.lasersSequence,15); // "3c5a96"
    		int numPointers;
    		int numRegisters=this.registerData[0].length;
			for (numPointers=0; (1<<numPointers)<this.registerData.length;numPointers++);
    		for (int i=0;i<this.registerData.length;i++) {
    			String onOff="";
    			for (int j=0; j<numPointers;j++) onOff+=(((i>>j) & 1) ==1)?" * ":" - ";
    			gd.addMessage(i+" Lasers:"+onOff);
    			for (int j=0;j<this.registerData[i].length;j++) {
    	    		gd.addStringField(j+": register address=  0x",Integer.toHexString(this.registerData[i][j][0]),4);
    	    		gd.addStringField(j+": register data=     0x",Integer.toHexString(this.registerData[i][j][1]),2);
    			}
    		}
    		if (askRegenerate) {
    			gd.addNumericField("Number of pointers (need to press REGENERATE button to change)", numPointers, 0);
    			gd.addNumericField("Number of registers (need to press REGENERATE button to change)", numRegisters, 0);
    			gd.enableYesNoCancel("OK", "REGENERATE");
    		}
    	    WindowTools.addScrollBars(gd);
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
    	    this.lasersIP=gd.getNextString();
    	    this.lasersSequence=gd.getNextString();
    	    
    		for (int i=0;i<this.registerData.length;i++) {
    			for (int j=0;j<this.registerData[i].length;j++) {
    				this.registerData[i][j][0]=Integer.parseInt(gd.getNextString(),16);
    				this.registerData[i][j][1]=Integer.parseInt(gd.getNextString(),16);
    			}
    		}
    		if (askRegenerate) {
    			int newNumPointers = (int) gd.getNextNumber();
    			int newNumRegisters = (int) gd.getNextNumber();
    			if (!gd.wasOKed() && ((newNumPointers!=numPointers) || (newNumRegisters!=numRegisters))) {
    				//modify
    				int [][][] oldRegisterData=this.registerData;
    				int numLines=oldRegisterData.length;
    				int newNumLines= 1 << newNumPointers;
    				this.registerData = new int [newNumLines][newNumRegisters][2];
    				for (int i=0;i<newNumLines;i++) {
    					int i0=(i<numLines)?i:(numLines-1);
    					for (int j=0;j<newNumRegisters;j++) {
    						int j0=(j<oldRegisterData[i0].length)?j:(oldRegisterData[i0].length-1);
    						this.registerData[i][j][0]=oldRegisterData[i0][j0][0];
    						this.registerData[i][j][1]=oldRegisterData[i0][j0][1];
    					}
    				}
    				return showDialog(title, false);
    			}
    		}
            return true;
    	}
    	public boolean setLasers (boolean [] lasers){
    		int intLasers=0;
    		for (int i=0;(i<lasers.length) && ((1<<i)<this.registerData.length);i++) intLasers |= lasers[i]?(1<<i):0;
    		return setLasers (intLasers);
    	}
    	public boolean setLasers (int intLasers){
    		if ((intLasers<0) || (intLasers>=this.registerData.length)) return false;
    		for (int i=0;i<this.registerData[intLasers].length; i++) {
    			String url="http://"+this.lasersIP+"/i2c.php?bus=1&raw=0x"+
    			Integer.toHexString(this.registerData[intLasers][i][0])+
    			"&data=0x"+Integer.toHexString(this.registerData[intLasers][i][1]);
    			if (this.debugLevel>2) System.out.println("setLasers: "+url); 
    			Document dom=null;
    			try {
    				DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
    				DocumentBuilder db = dbf.newDocumentBuilder();
    				dom = db.parse(url);
    				if (!dom.getDocumentElement().getNodeName().equals("i2c")) {
    					System.out.println("Root element: expected 'i2c', got'" + dom.getDocumentElement().getNodeName()+"'");
    					IJ.showMessage("Error","Root element: expected 'i2c', got'" + dom.getDocumentElement().getNodeName()+"'"); 
    					return false;
    				}

    				//Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor1").item(0)).getChildNodes().item(0))).getNodeValue());
    				boolean responceError= (dom.getDocumentElement().getElementsByTagName("error").getLength()!=0);
    				if (responceError) {
    					System.out.println("ERROR: register write ("+url+") FAILED" );
    					IJ.showMessage("Error","register write ("+url+") FAILED"); 
    					return false;  
    				}
    			} catch(MalformedURLException e){
    				System.out.println("Please check the URL:" + e.toString() );
    				return false;
    			} catch(IOException  e1){
    				IJ.showStatus("");
    				String error = e1.getMessage();
    				if (error==null || error.equals(""))  error = ""+e1;
    				IJ.showMessage("setLasers ERROR", ""+error);
    				return false;
    			}catch(ParserConfigurationException pce) {
    				pce.printStackTrace();
    				return false;
    			}catch(SAXException se) {
    				se.printStackTrace(); 
    				return false;
    			}
    		}
    		//http://192.168.0.13/i2c.php?bus=1&raw=0x2700&data=0xdf    		
    		this.state=intLasers;
    		return true;
    	}
    	public boolean manualSetLasers(){
    		int numPointers;
			for (numPointers=0; (1<<numPointers)<this.registerData.length;numPointers++);
            int oldState=(this.state>=0)?this.state:0; // now showing -1 (undefined) as 0 - off
    		boolean [] lasers = new boolean[numPointers];
    		for (int i=0;i<lasers.length;i++) lasers[i]= ((oldState>>i) & 1)!=0; 
    		GenericDialog gd = new GenericDialog("Turn laser pointers on/off");
    		for (int i=0; i<lasers.length;i++) {
        		gd.addCheckbox("Laser pointer "+(i+1),lasers[i]);
    		}
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
    		for (int i=0; i<lasers.length;i++) lasers[i]=gd.getNextBoolean();
    		return setLasers (lasers);
    	}
    	
    }
 
    public static class UVLEDandLasers{
    	public int debugLevel=1;
    	private int state=-1; // bit-wise on/off, <0 - undefined
    	private String uvLasersIP="192.168.0.236";
    	public int uvLasersBus=0;
    	public double [] uvLasersCurrents={0,0,0,0}; // will be overwritten
    	public boolean [] laserState=null;
    	public boolean [] uvState=null;
    	public double maxCurrent=100; //mA
    	public long [] lastExposureStart={-1,-1,-1,-1}; // timestamp the LED was last turned to this current
    	public double [] runningCurrent={0.0,0.0,0.0,0.0}; // currently running current, mA
    	public double [] ampsSeconds={0.0,0.0,0.0,0.0}; // cumulative Amps*seconds
    	public UVLEDandLasers(
    			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    	}
    	public void setParameters(
    			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		this.uvLasersIP=focusMeasurementParameters.uvLasersIP;
    		this.uvLasersBus=focusMeasurementParameters.uvLasersBus;
    		this.uvLasersCurrents=focusMeasurementParameters.uvLasersCurrents;
    		this.ampsSeconds=focusMeasurementParameters.ampsSeconds;
			for (int i=0;i<this.uvLasersCurrents.length;i++) {
				if (this.uvLasersCurrents[i]>this.maxCurrent) this.uvLasersCurrents[i]=this.maxCurrent;
				if (this.uvLasersCurrents[i]<0.0) this.uvLasersCurrents[i]=0.0;
			}
    	}
    	public void getParameters(
    			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		focusMeasurementParameters.uvLasersIP=this.uvLasersIP;
    		focusMeasurementParameters.uvLasersBus=this.uvLasersBus;
    		focusMeasurementParameters.uvLasersCurrents=this.uvLasersCurrents;
    		
    	}

    	public boolean uvLaserSettings(
    			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		boolean result=uvLaserSettings();
    		getParameters(focusMeasurementParameters);
    		return result;
    	}
    	
    	public boolean uvControl(
    			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		boolean result=uvControl();
    		getParameters(focusMeasurementParameters);
    		return result;
    	}
    	
    	public void updateCurrents(){
    		long thisTime=System.nanoTime();
    		for (int i=0;i<uvState.length;i++){
    			if (this.debugLevel>3) System.out.println("LED"+i+" runningCurrent="+runningCurrent[i]+" time passed="+
    					(thisTime-lastExposureStart[i])+" amps-sec="+ampsSeconds[i]);
    			if (runningCurrent[i]>0) ampsSeconds[i]+=(0.001*runningCurrent[i])*(0.000000001*(thisTime-lastExposureStart[i]));
    			lastExposureStart[i]=thisTime;
    			runningCurrent[i]=(uvState[i])?uvLasersCurrents[i]:0.0;
    			if (this.debugLevel>3) System.out.println("LED"+i+" new runningCurrent="+runningCurrent[i]);
    		}
    	}
    	public void resetCumulativeCurrents(){
    		long thisTime=System.nanoTime();
    		for (int i=0;i<uvState.length;i++){
    			ampsSeconds[i]=0.0;
    			lastExposureStart[i]=thisTime;
    			runningCurrent[i]=(uvState[i])?uvLasersCurrents[i]:0.0;
    		}    		
    	}
///55 504 690 774
//3 220 482 180 430 159
/**
 *     	returns true if any of the UV LEDs is (or was if dialog canceled) on
 */
    	public boolean uvControl(){
    		updateStatus();
    		GenericDialog gd = new GenericDialog("UV LED Control");
    		gd.addMessage("IP address="+this.uvLasersIP+" i2c bus="+this.uvLasersBus);
    		gd.addMessage("UV LEDs are referenced when viewd from the target to the camera");
    		gd.addCheckbox("UV LED 0 (left/near), "+this.uvLasersCurrents[0]+" mA (total "+IJ.d2s(ampsSeconds[0],3)+" Amp-sec)",this.uvState[0]);
    		gd.addCheckbox("UV LED 1 (right/near),"+this.uvLasersCurrents[1]+" mA (total "+IJ.d2s(ampsSeconds[1],3)+" Amp-sec)",this.uvState[1]);
    		gd.addCheckbox("UV LED 2 (right/far), "+this.uvLasersCurrents[2]+" mA (total "+IJ.d2s(ampsSeconds[2],3)+" Amp-sec)",this.uvState[2]);
    		gd.addCheckbox("UV LED 3 (left/far),  "+this.uvLasersCurrents[3]+" mA (total "+IJ.d2s(ampsSeconds[3],3)+" Amp-sec)",this.uvState[3]);
    		gd.addCheckbox("Turn them all",false);
    		gd.addCheckbox("Reset cumulative Amps-seconds",false);
    		gd.enableYesNoCancel("Apply", "Change Currents, then apply");
    		boolean anyOn=false;
    		for (int i=0;i<uvState.length;i++ ) anyOn|=uvState[i];
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return anyOn; 
    		if (!gd.wasOKed()) {
    			if (!uvLaserSettings()) return anyOn;
    		}
    		this.uvState[0]=gd.getNextBoolean();
    		this.uvState[1]=gd.getNextBoolean();
    		this.uvState[2]=gd.getNextBoolean();
    		this.uvState[3]=gd.getNextBoolean();
    		if (gd.getNextBoolean()){
        		this.uvState[0]=true;
        		this.uvState[1]=true;
        		this.uvState[2]=true;
        		this.uvState[3]=true;
    			
    		}
    		if (gd.getNextBoolean())resetCumulativeCurrents();
        	setLasersAndUV(
        			null, // lasers - do not touch
        			this.uvState,     // may be null
        			this.uvLasersCurrents// may be null
        			);
    		anyOn=false;
    		for (int i=0;i<uvState.length;i++ ) anyOn|=uvState[i];
        	return anyOn;
    	}
    	
    	
/*
 *     	public double [] uvLasersCurrents={0,0,0,0}; // will be overwritten
    	public boolean [] laserState=null;
    	public boolean [] uvState=null;
    	
 */
    	public boolean uvLaserSettings(){
    		GenericDialog gd = new GenericDialog("UV LED and laser configuration (103641 board)");
			gd.addStringField  ("IP address of the camera with 103641 board (UV LEDs and lasers) are attached",      this.uvLasersIP,40);
    		gd.addNumericField("I2C bus where LED/laser board is attached (0 - through 10359, 1 - through 10369)",   this.uvLasersBus,        0);
//    		gd.addMessage("New value for LED current will not be applied to the turned on device in this dialog");
    		gd.addNumericField("UV LED1 \"on\" current (left/near  when looking from the target)",                   this.uvLasersCurrents[0],  3,5,"mA");
    		gd.addNumericField("UV LED2 \"on\" current (right/near when looking from the target)",                   this.uvLasersCurrents[0],  3,5,"mA");
    		gd.addNumericField("UV LED3 \"on\" current (right/far  when looking from the target)",                   this.uvLasersCurrents[0],  3,5,"mA");
    		gd.addNumericField("UV LED4 \"on\" current (left/far   when looking from the target)",                   this.uvLasersCurrents[0],  3,5,"mA");
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
			this.uvLasersIP=                 gd.getNextString();
			this.uvLasersBus=          (int) gd.getNextNumber();
			this.uvLasersCurrents[0]=        gd.getNextNumber();
			this.uvLasersCurrents[1]=        gd.getNextNumber();
			this.uvLasersCurrents[2]=        gd.getNextNumber();
			this.uvLasersCurrents[3]=        gd.getNextNumber();
			for (int i=0;i<this.uvLasersCurrents.length;i++) {
				if (this.uvLasersCurrents[i]>this.maxCurrent) this.uvLasersCurrents[i]=this.maxCurrent;
				if (this.uvLasersCurrents[i]<0.0) this.uvLasersCurrents[i]=0.0;
			}
    		return true;
    	}
    	
    	
    	
    	public boolean allOff(LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		return allOff();
    	}
    	public boolean uvOff(LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		return uvOff();
    	}
    	public boolean lasersToggle (LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
        	return lasersToggle();
    	}
    	
    	public boolean lasersOff(LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		return lasersOff();
    	}
    	
    	public boolean lasersOn(LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		return lasersOn();
    	}
    	
    	public boolean setLasersAndUV(
    			boolean [] lasersOn, // may be null
    			boolean [] uvOn,     // may be null
    			double [] uvCurrents,// may be null
    			LensAdjustment.FocusMeasurementParameters focusMeasurementParameters){
    		setParameters(focusMeasurementParameters);
    		return setLasersAndUV(lasersOn,uvOn,uvCurrents);
    	}
    	
    	public boolean allOff(){
    		boolean [] lasersOn={false,false};
    		boolean [] uvOn={false,false,false,false};
    		return setLasersAndUV(lasersOn,uvOn,null);
    	}
    	public boolean uvOff(){
    		boolean [] uvOn={false,false,false,false};
    		return setLasersAndUV(null,uvOn,null);
    	}

    	public boolean lasersToggle(){
    		updateStatus();
    		for (int i=0;i<this.laserState.length;i++) this.laserState[i]=!this.laserState[i];
    		return setLasersAndUV(this.laserState,null,null);
    	}

    	public boolean lasersOff(){
    		boolean [] lasersOn={false,false};
    		return setLasersAndUV(lasersOn,null,null);
    	}
    	public boolean lasersOn(){
    		boolean [] lasersOn={true,true};
    		return setLasersAndUV(lasersOn,null,null);
    	}
    	public boolean setLasersAndUV(
    			boolean [] lasersOn, // may be null
    			boolean [] uvOn,     // may be null
    			double [] uvCurrents// may be null
    			){
			if (this.debugLevel>2) System.out.println("lasersOn="+((lasersOn==null)?"null":("{"+lasersOn[0]+","+lasersOn[1]+"}"))+
					" uvOn="+((uvOn==null)?"null":("{"+uvOn[0]+","+uvOn[1]+","+uvOn[2]+","+uvOn[3]+"} "))+
					" uvCurrents="+((uvCurrents==null)?"null":("{"+uvCurrents[0]+","+uvCurrents[1]+","+uvCurrents[2]+","+uvCurrents[3]+"} "))); 

    		String command="";
    		if (lasersOn!=null)	for (int i=0;i<lasersOn.length;i++) {
    			command+="&laser"+i+"="+(lasersOn[i]?1:0);
    		}
    		if (uvOn!=null)	for (int i=0;i<uvOn.length;i++) {
    			int data=(uvOn[i] && (uvCurrents!=null) && (uvCurrents.length>i))?((int) Math.round(uvCurrents[i])):0;
    			command+="&uv"+i+"="+data;
    		}
    		
    		boolean result=commandToDevice(command);
    		updateCurrents();
    		return result;

    	}
    	public boolean [] getLasers(){
    		updateStatus();
    		return this.laserState;
    	}
    	/**
    	 * Will send an empty command (just i2c bus number) and update current state of LEDs and lasers (if they were modified by other master)
    	 * @return true if OK, false - failure (not yet used, now all errors throw)
    	 */
    	public boolean updateStatus(){
    		return commandToDevice("");
    	}
    	public boolean commandToDevice(String command){
    			String url="http://"+this.uvLasersIP+"/103641.php?bus="+this.uvLasersBus+command; // should start with "&"
    			if (this.debugLevel>1) System.out.println("setDevice: "+url); 
    			Document dom=null;
    			try {
    				DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
    				DocumentBuilder db = dbf.newDocumentBuilder();
    				dom = db.parse(url);
    				if (!dom.getDocumentElement().getNodeName().equals("board_103641")) {
    					String msg="Root element: expected \"board_103641\", got\"" + dom.getDocumentElement().getNodeName()+"\"";
    					IJ.showMessage("Error",msg); 
    					throw new IllegalArgumentException (msg);
    				}
    				String quoted=(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("data9500").item(0)).getChildNodes().item(0))).getNodeValue());
    				// remove opening and closing "
    				if (quoted.startsWith("\"")){
    					quoted=quoted.substring(1, quoted.length()-1);
    				}
    				this.state=Integer.decode(quoted);
    				if (this.state==0){ // should normally not happen
    					String msg="Data read from the board is 0, that is not currently possible if everything works.";
    					IJ.showMessage("Error",msg); 
    					throw new IllegalArgumentException (msg);
    				}
    			} catch(MalformedURLException e){
    				String msg="Please check the URL:" + e.toString();
					IJ.showMessage("Error",msg); 
					throw new IllegalArgumentException (msg);
    			} catch(IOException  e1){
    				String msg = e1.getMessage();
    				if (msg==null || msg.equals(""))  msg = ""+e1;
					IJ.showMessage("Error",msg); 
					throw new IllegalArgumentException (msg);
    			}catch(ParserConfigurationException pce) {
    				pce.printStackTrace();
    				return false;
    			}catch(SAXException se) {
    				se.printStackTrace(); 
    				return false;
    			}
    			boolean [] laserState={
    				(this.state & 0x10)==0,
    				(this.state & 0x20)==0};
    			boolean [] uvState={
        				(this.state & 0x01)==0,
        				(this.state & 0x02)==0,
        				(this.state & 0x04)==0,
        				(this.state & 0x08)==0};
    		//http://192.168.0.13/i2c.php?bus=1&raw=0x2700&data=0xdf    	
    	    	this.laserState=laserState;
    	    	this.uvState=uvState;
    		return true;
    	}
    	
    }
    public static class GoniometerMotors{
        public double [][] motorsRange={null,{-100000,100000},{-500000,500000}};
        public double stepsPerSecond=     732.7; // motor speed limit in steps per second
        public double stepsPerDegreeTilt=-5682.48889; // minus that positive steps make negative elevation
        
//        gearmotor -64steps 1:19 - 1216 steps/turn. Worm 1:80 -> 97280/turn. roller/tube ~1.44 ->388.3/degree
//        public double stepsPerDegreeAxial=-388.3; // -36.0; // minus that positive steps make rotate CCW when looking from Eyesis top
        public double stepsPerDegreeAxial=-369.64; // measured 133070 steps per360
        
    	public String ipAddress="192.168.0.109";
        public int [] curpos=new int[3];
        public int [] targetPosition={0,0,0};
        public int    debugLevel=2;
        public int motorTolerance=5; // steps - disregard error less than that
        public int motorStuckTolerance=5; // steps - disregard error less than that
        
        public double coefficientETA=1.5;  // allow moving 1.5 longer than at maximal speed 
        
        private long nanoETA; 
        private long nanoReferenceTime; // last time the position was checked
        private int []referencePosition=null;
        private double motorsStuckTestTime=5.0; // seconds
        private boolean motorsInitialized=false;
        public int tiltMotor=2;  // 0-1-2
        public int axialMotor=1; // 0-1-2
        

        /**
         * Tries to initialize motors if needed. Assumes that non-initialized motors were not tried to be moved 
         * @param force unconditionally try to initialize
         * @return
         */
        public boolean tryInit(
        		boolean force,
        		boolean updateStatus){
        	int delta=2*this.motorTolerance;
        	if (!force) {
        		if (motorsInitialized) return true;
        		updateMotorsPosition();
        		if ((Math.abs(curpos[tiltMotor])>2) || (Math.abs(curpos[axialMotor])>2)){
        			motorsInitialized=true;
        			return true;
        		}
        		System.out.println("Trying to move axial motor to see if it is initialized ...");
        		motorsInitialized=true; // to avoid loop
        		int newpos=curpos[axialMotor]+delta;
        		if (moveMotorSetETA(axialMotor, newpos)) {
            		if (waitMotor(axialMotor, null, true, updateStatus)) {// no interaction
            			moveMotorSetETA(axialMotor, newpos-delta); // restore original
                		if (waitMotor(axialMotor, null, true, updateStatus)) {// no interaction
                			return true; // should be always true here 
                		}
            		}
        		}
        		motorsInitialized=false;
        	}
        	initMotors();
        	motorsInitialized=true;
        	return true;
        }
// TODO: Enable simultaneous motors        
        
        public boolean moveMotorSetETA(int motorNumber, int position){
        	if ((motorNumber<0) ||(motorNumber>=this.motorsRange.length) || (this.motorsRange[motorNumber]==null)){
        		String msg="Motor "+motorNumber+" is undefined";
        		IJ.showMessage("Error",msg);
        		throw new RuntimeException(msg);
        	}
        	if ((position<this.motorsRange[motorNumber][0]) || (position>this.motorsRange[motorNumber][1])){
        		String msg="Motor "+motorNumber+" requested position "+position+" is out of range ("+
        		this.motorsRange[motorNumber][0]+"..."+this.motorsRange[motorNumber][1]+")";
        		System.out.println("Error: "+msg);
        		IJ.showMessage("Error",msg);
        		//throw new RuntimeException(msg);
        		return false;
        		
        	}
        	if (!motorsInitialized) { // may be called from tryInit(), but with motorsInitialized set to true;
        		if (!tryInit(false, true)){ // do not force, update status
            		String msg="Motors failed to initialize";
                    		System.out.println("Error: "+msg);
                    		IJ.showMessage("Error",msg);
                    		//throw new RuntimeException(msg);
                    		return false;
        		}
        	}
        	// TODO: check init is needed
        	
        	updateMotorsPosition();
        	this.targetPosition[motorNumber]=position;
        	this.nanoReferenceTime=System.nanoTime();
        	this.referencePosition=this.curpos.clone();
//        	commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?m"+(motorNumber+1)+"="+this.targetPosition[motorNumber]+"&enable");
        	// set, wait, enable (to avoid continuing "old" movement that was aborted by disable
        	commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?m"+(motorNumber+1)+"="+this.targetPosition[motorNumber]+"sleep=1");
        	enableMotors(true); 
			nanoETA=System.nanoTime()+((long)(1E9*(Math.abs(this.targetPosition[motorNumber]-this.curpos[motorNumber])*(this.coefficientETA/this.stepsPerSecond))));
			return true;
        }
/*        
        public boolean [] checkGotTarget(){
    		updateMotorsPosition(0); // no wait here
    		boolean [] result=new boolean [this.curpos.length];
    		for (int motorNumber=0;motorNumber<result.length;motorNumber++){
    			result[motorNumber]=Math.abs(this.targetPosition[motorNumber]-this.curpos[motorNumber])<this.motorTolerance;
    		}
    		return result;
        }
*/
        public boolean checkGotTarget(int motorNumber, int position){
    		updateMotorsPosition(0); // no wait here
    		return Math.abs(position-this.curpos[motorNumber])<this.motorTolerance;
        }

        public boolean waitMotor(
        		int motorNumber,
        		AtomicInteger stopRequested, // or null
        		boolean quiet,
        		boolean updateStatus
        		){
        	if ((motorNumber<0) ||(motorNumber>=this.motorsRange.length) || (this.motorsRange[motorNumber]==null)){
        		String msg="Motor "+motorNumber+" is undefined";
        		IJ.showMessage("Error",msg);
        		throw new RuntimeException(msg);
        	}
        	while (true) {
//    			enableMotors(true); // just in case? - not here, test first
        		updateMotorsPosition(1); // wait one second before testing to decrease re-test frequency
        		if (updateStatus) {
        			double axialAngle=curpos[axialMotor]/this.stepsPerDegreeAxial;
        			double tiltAngle=curpos[tiltMotor]/this.stepsPerDegreeTilt;
        			double axialTargetAngle=targetPosition[axialMotor]/this.stepsPerDegreeAxial;
        			double tiltTargetAngle=targetPosition[tiltMotor]/this.stepsPerDegreeTilt;
        			IJ.showStatus("Goniometer: tilt="+IJ.d2s(tiltAngle,1)+"\u00b0 ("+IJ.d2s(tiltTargetAngle,1)+"\u00b0) "+
        					", axial="+IJ.d2s(axialAngle,1)+"\u00b0 ("+IJ.d2s(axialTargetAngle,1)+"\u00b0) "+
        					" Tilt steps="+curpos[tiltMotor]+" ("+this.targetPosition[tiltMotor]+") "+
        					", axial steps="+curpos[axialMotor]+" ("+this.targetPosition[axialMotor]+")");
        		}
        		int positionError= Math.abs(this.targetPosition[motorNumber]-this.curpos[motorNumber]);
        		if (positionError<this.motorTolerance){
        			updateMotorsPosition(1); // re-test
        			positionError= Math.abs(this.targetPosition[motorNumber]-this.curpos[motorNumber]);
        			enableMotors(false);
        			return true;
        		}
    			enableMotors(true); // does not need to be enabled before
        		long nanoNow=System.nanoTime();
        		if ((stopRequested!=null) && (stopRequested.get()>1)){
        			enableMotors(false);
					if (this.debugLevel>0) System.out.println("User interrupt");
					stopRequested.set(0);
					updateMotorsPosition(1);
        			double axialAngle=curpos[axialMotor]/this.stepsPerDegreeAxial;
        			double tiltAngle=curpos[tiltMotor]/this.stepsPerDegreeTilt;
        			double axialTargetAngle=targetPosition[axialMotor]/this.stepsPerDegreeAxial;
        			double tiltTargetAngle=targetPosition[tiltMotor]/this.stepsPerDegreeTilt;
        			String msg="Goniometer: tilt="+IJ.d2s(tiltAngle,1)+"\u00b0 ("+IJ.d2s(tiltTargetAngle,1)+"\u00b0) "+
        					", axial="+IJ.d2s(axialAngle,1)+"\u00b0 ("+IJ.d2s(axialTargetAngle,1)+"\u00b0) "+
        					" Tilt steps="+curpos[tiltMotor]+" ("+this.targetPosition[tiltMotor]+") "+
        					", axial steps="+curpos[axialMotor]+" ("+this.targetPosition[axialMotor]+")\n"+
        					"OK to continue, Cancel to abort movement";
        			if (!IJ.showMessageWithCancel("Goniometer interrupted", msg)) {
        				return false;
        			}
        			enableMotors(true);
        			this.nanoETA += System.nanoTime()-nanoNow;
        		}
        		if (nanoNow>this.nanoETA) {
        			enableMotors(false);
        			if (quiet) return false;
        			String msg="Motor "+motorNumber+" failed to reach destination "+this.targetPosition[motorNumber]+
        			", current position is "+this.curpos[motorNumber]+
        			"\nYou may try to manually fix the problem before hitting OK";
        			System.out.println ("Error:"+msg);
        			IJ.showMessage("Error",msg);
//        			if (IJ.showMessageWithCancel("Goniometer did not move", msg+"\n OK will try to initialize ");
// Give chance to manually fix the problem        			
        			updateMotorsPosition();
        			positionError= Math.abs(this.targetPosition[motorNumber]-this.curpos[motorNumber]);
        			if (positionError<this.motorTolerance){
            			enableMotors(false);
            			return true;
        			}
        			return false;
        			
        		}
        		if ((nanoNow-this.nanoReferenceTime)> ((long) (this.motorsStuckTestTime*1E9))){
        			if(Math.abs(this.referencePosition[motorNumber]-this.curpos[motorNumber]) < motorStuckTolerance){
            			enableMotors(false);
            			if (quiet) return false;
            			String msg="Motor "+motorNumber+" is stuck at "+this.curpos[motorNumber]+". "+
            			this.motorsStuckTestTime+" seconds ago it was at "+this.referencePosition[motorNumber]+
            			", target position is "+this.targetPosition[motorNumber]+
            			"\nYou may try to manually fix the problem before hitting OK";
            			System.out.println ("Error:"+msg);
            			IJ.showMessage("Error",msg);
            			// Give chance to manually fix the problem        			
            			updateMotorsPosition();
            			positionError= Math.abs(this.targetPosition[motorNumber]-this.curpos[motorNumber]);
            			if (positionError<this.motorTolerance){
                			enableMotors(false);
                			return true;
            			}
            			return false;
        			} else {
        	        	this.nanoReferenceTime=System.nanoTime();
        	        	this.referencePosition=this.curpos.clone();
        			}
        		}
        	}
        }

// first check tolerance, then - if motor is stuck        
        public int [] getCurrentPositions(){
        	updateMotorsPosition();
        	return this.curpos;
        }

        public int [] getTargetPositions(){
        	return this.targetPosition;
        }
        public int[] enableMotors(boolean enable)  {
        	return commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?"+(enable?"enable":"disable"));
        }
        public int[] initMotors()  {
        	return commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?init");
        }
        public int[] setHome()  {
        	return commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?m1=0&m2=0&m3=0&reset");
        }
        
        public int[] updateMotorsPosition()  {
        	return commandElphel10364Motors("http://"+this.ipAddress+"/10364.php");
        }
        public int[] updateMotorsPosition(int sleep)  {
        	return commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?sleep="+sleep);
        }

        private int[] commandElphel10364Motors(String url)  {
        	Document dom=null;
        	if (this.debugLevel>2)	System.out.println("commandElphel10364Motors("+url+")");
        	try {
        		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        		DocumentBuilder db = dbf.newDocumentBuilder();
        		dom = db.parse(url);
        		if (!dom.getDocumentElement().getNodeName().equals("motors")) {
        			IJ.showMessage("Root element: expected 'motors', got'" + dom.getDocumentElement().getNodeName()+"'"); 
        			return null;
        		}
        		this.curpos[0]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor1").item(0)).getChildNodes().item(0))).getNodeValue());
        		this.curpos[1]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor2").item(0)).getChildNodes().item(0))).getNodeValue());
        		this.curpos[2]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor3").item(0)).getChildNodes().item(0))).getNodeValue());
        	} catch(MalformedURLException e){
        		System.out.println("Please check the URL:" + e.toString() );
        		return null;
        	} catch(IOException  e1){
        		IJ.showStatus("");
        		String error = e1.getMessage();
        		if (error==null || error.equals(""))  error = ""+e1;
        		IJ.showMessage("commandElphel10364Motors ERRROR", ""+error);
        		return null;
        	}catch(ParserConfigurationException pce) {
        		pce.printStackTrace();
        		return null;
        	}catch(SAXException se) {
        		se.printStackTrace(); 
        		return null;
        	}
        	// prevent accumulating position errors
/*        	
        	if ((this.lastMovedTo!=null) && (Math.abs(this.lastMovedTo[0]-this.curpos[0])<=this.positionTolerance) &&
        			(Math.abs(this.lastMovedTo[1]-this.curpos[1])<=this.positionTolerance) &&
        			(Math.abs(this.lastMovedTo[2]-this.curpos[2])<=this.positionTolerance)) this.curpos=this.lastMovedTo.clone();  
*/
        	return this.curpos;
        }
        
    	
    }
    
    
    public static class FocusingMotors{ // should it be static?
    	
    public int motor3rd_signed_hysteresis=100;
//    public double [][] motorsRange={{-32000,32000},{-32000,32000},{-32000,32000}}; // hardware limited to 16 bits signed
    public double [][] motorsRange={{-32000,25000},{-32000,25000},{-32000,25000}}; // need longer spring, bad hysteresis at the center, starting at +12000
	public String ipAddress="192.168.0.236";
	public int stepsPerRevolution=3584; //112*32
	public double threadPitch=0.4; // mm/revolution
	public double linearReductionRatio=4.0/38.0; // sensor movement per 3 screws together movement
	public double waitPerStep=0.005; //~ 2 sec/revolution
	public double waitExtra=1.0;
	public double calmMotors=1.0; // extra time to wait to calm
	public int    positionTolerance=2; // actual position should not differ than this from the required
	public int [] lastMovedTo=null;
    public int[] curpos=new int[3];
    public boolean [] lastDirs= {false,false,false}; // true if motor was going in positive direction
    public FocusingHistory focusingHistory=new FocusingHistory(1000);
    public FocusingSharpness focusingSharpness=new FocusingSharpness(1000);
    public int debugLevel=1;
    public String motorsStatePath="lensAdjustmentMotorPosition.xml"; // will be saved in ~/.imagej/
    public String prefsDirectory=null; // prefs.getPrefsDir()+Prefs.getFileSeparator(); // "~/.imagej/" - too early !!! 
    private Properties motorProperties=new Properties();
    public boolean interactiveRestore=true; // ask for confirmation to restore motor position
    public boolean stateValid=false;
    public boolean motorHardwareInitialized=false; // check if motors need to be initialized, enabled and reset (cold start of the board)
    public boolean motorHardwareInitializationInProgress=false;
    public boolean motorsEnabled=false; // disable motors when UV was turned on
    
    public double [][] getRange(){return this.motorsRange;}
    public int getStepsPerRevolution(){return this.stepsPerRevolution;}
    public double getThreadPitch(){return this.threadPitch;}
    public double getLinearReductionRatio(){return this.linearReductionRatio;}
    public double getMicronsPerStep(){return  1000.0*this.threadPitch*this.linearReductionRatio/stepsPerRevolution;}
    public void   setLinearReductionRatio(double lrr){this.linearReductionRatio=lrr;}
    
    public String getPrefsDir(){
    	if (prefsDirectory==null) prefsDirectory=Prefs.getPrefsDir();
    	if (prefsDirectory!=null) prefsDirectory+=Prefs.getFileSeparator(); // "~/.imagej/"
    	return prefsDirectory;
    }

    public boolean isInitialized(){
    	return this.stateValid;
    }
    public boolean isEnabled(){
    	return this.motorsEnabled;
    }
    public void setEnable(boolean enable){
    	this.motorsEnabled=enable;
    }
    /**
     *  Run after the camera was turned off while program is still running
     */
    public void resetInitialization(){
    	this.stateValid=false;
        this.motorHardwareInitialized=false; // check if motors need to be initialized, enabled and reset (cold start of the board)
        this.motorHardwareInitializationInProgress=false;
        this.motorsEnabled=false; // disable motors when UV was turned on
    }

    
    public void setCalmMotors (double calmMotors){
    	this.calmMotors=calmMotors;
    }
    
    
    public void setHysteresis (int hysteresis){
    	this.motor3rd_signed_hysteresis=hysteresis;
    }
    public void setDebug (int debugLevel){
    	this.debugLevel=debugLevel;
    	System.out.println("Motor debug level is set to "+this.debugLevel);
    }

    public int[] moveElphel10364Motors(boolean wait, int [] positions, double sleep, boolean showStatus, boolean hysteresis)  {
    	 return moveElphel10364Motors(wait, positions, sleep, showStatus, "", hysteresis);
    }

    public void saveMotorState() throws IOException{
    	if (!this.stateValid) return; // do not save until allowed 
    	for (int i=0;i<curpos.length;i++) this.motorProperties.setProperty("motor"+(i+1),this.curpos[i]+"");
    	String path=this.getPrefsDir()+this.motorsStatePath;
    	OutputStream os;
    	try {
    		os = new FileOutputStream(path);
    	} catch (FileNotFoundException e1) {
    		String msg="Failed to open configuration file: "+path;
    		IJ.showMessage("Error",msg);
    		throw new FileNotFoundException (msg);
    	}
/*    	if (os==null) {
    		String msg="Failed to open configuration file for writing: "+path;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
*/    	
    	try {
    		this.motorProperties.storeToXML(os,
    				"last updated " + new java.util.Date(), "UTF8");
    	} catch (IOException e) {
    		String msg="Failed to write XML motor state file: "+path;
    		IJ.showMessage("Error",msg);
    		throw new IOException (msg);
    	}
    	try {
    		os.close();
    	} catch (IOException e) {
    		// TODO Auto-generated catch block
    		e.printStackTrace();
    	}
    }
    
    public int [] loadMotorState() throws IOException{
    	String path=this.getPrefsDir()+this.motorsStatePath;
    	InputStream is;
    	try {
    		is = new FileInputStream(path);
    	} catch (FileNotFoundException e) {
    		String msg="Failed to open configuration file: "+path;
    		System.out.println("Warning: "+msg);
    		return null; // file does not exist, use current motor position (if it is initialized)
    	}
    	try {
    		motorProperties.loadFromXML(is);

    	} catch (IOException e) {
    		String msg="Failed to read XML configuration file: "+path;
    		IJ.showMessage("Error",msg);
    		throw new IOException (msg);
    	}
    	try {
    		is.close();
    	} catch (IOException e) {
    		// TODO Auto-generated catch block
    		e.printStackTrace();
    	}
		int [] savedPosition=new int [this.curpos.length];
    	for (int i=0;i<curpos.length;i++){
    		String sValue=	this.motorProperties.getProperty("motor"+(i+1));
    		if (sValue!=null) {
    			savedPosition[i]=Integer.parseInt(sValue);
//    			if (this.debugLevel>2) System.out.println("Read savedPosition["+i+"]="+savedPosition[i]);
    			System.out.println("Read savedPosition["+i+"]="+savedPosition[i]);
    		} else {
    			String msg="motor"+(i+1)+" is undefined in "+path+". If the file is corrupted you may delete it.";
    			IJ.showMessage("Error",msg);
    			throw new IOException (msg);
    		}
    	}
    	return savedPosition;
    }
    public void checkEnabled(){
    	if (isEnabled()) return;
    	GenericDialog gd = new GenericDialog("Enable motors");
		gd.addMessage("Motors are currently disabled (i.e. SFE is being glued). Do you want to enable them?");
	    gd.showDialog();
	    if (gd.wasCanceled()) {
	    	throw new RuntimeException("Command canceled");
	    }
	    setEnable(true);
    }
    /**
     *  Check if the motor position is synchronized with disk, open dialog if needed (or throw) (use boolean flag - dialog or throw)
     */
    public void checkSynchronized(){
     	if (this.stateValid) return; // nothing to do
    	// first check if the hardware is initialized, if not - do that.
    	if (motorHardwareInitializationInProgress) return; 
    	motorHardwareInitializationInProgress=true;
		boolean wasEnabled= true; // innocent until proven guilty
    	if (!motorHardwareInitialized){
    		readElphel10364Motors(); // will update this.curpos
    		int [] testPosition=this.curpos.clone();
    		testPosition[0]+=2*this.positionTolerance;
    		double savedCalmMotors=this.calmMotors;
    		this.calmMotors=2;
        	moveElphel10364Motors(
        			true, // will not wait forever - will time out
        			testPosition, //int [] positions,
        			0.0, // double sleep,
        			true, // boolean showStatus,
        			"Looking if motors are initialized - trying to move motor1", // String message,
        			false); // , boolean hysteresis)
    		wasEnabled= maxPosDiff(testPosition)<=this.positionTolerance;
        	if (!wasEnabled){ // seems not to be enabled, and supposedly not initialized
        		String msg="Initializing motor driver";
    			if (this.debugLevel>0) System.out.println(msg);
    			IJ.showStatus("Initializing motor driver");
            	commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?init&enable");
            	moveElphel10364Motors(
            			true, //will timeout if noit initialized
            			testPosition, //int [] positions,
            			0.0, // double sleep,
            			true, // boolean showStatus,
            			"Trying to move motor1", // String message,
            			false); // , boolean hysteresis)
            	if (maxPosDiff(testPosition)>this.positionTolerance){//  failed to initialize
            		testPosition[0]-=2*this.positionTolerance; // restore back this.curpos may be wrong (probably not)
                	this.calmMotors=savedCalmMotors;
                	moveElphel10364Motors(
                			true, //boolean wait (Will be forever, if not initialized)
                			testPosition, //int [] positions,
                			0.0, // double sleep,
                			true, // boolean showStatus,
                			"Trying to move motor1", // String message,
                			false); // , boolean hysteresis)
                	motorHardwareInitializationInProgress=false;
            		String msg1="Failed to initialize motor driver board";
        			IJ.showMessage("Error",msg);
	    	    	throw new RuntimeException(msg1);
            	}
        	}
    		testPosition[0]-=2*this.positionTolerance; // restore back this.curpos may be wrong (probably not)
        	this.calmMotors=savedCalmMotors;
        	moveElphel10364Motors(
        			true, //boolean wait (Will be forever, if not initialized)
        			testPosition, //int [] positions,
        			0.0, // double sleep,
        			true, // boolean showStatus,
        			"Trying to move motor1", // String message,
        			false); // , boolean hysteresis)
        	motorHardwareInitializationInProgress=false;
    		motorHardwareInitialized=true;
    	}
    	
    	int [] savedPosition=null;
   // 	if (this.stateValid) return; //
    	if (!this.stateValid) try {
    		savedPosition=loadMotorState();
    	} catch (IOException e) {
    		// IJ.showMessage is already run
    		e.printStackTrace();
    	}
    	// is it the same position
    	GenericDialog gd;
    	if (this.debugLevel>2) System.out.println("savedPosition="+((savedPosition==null)?"null":("{"+savedPosition[0]+","+savedPosition[1]+","+savedPosition[2]+"}"))+
    			" motorHardwareInitialized="+motorHardwareInitialized+ " motorHardwareInitializationInProgress="+motorHardwareInitializationInProgress);
    	if (savedPosition!=null){
    		if (maxPosDiff(savedPosition)>this.positionTolerance){
    			if (this.interactiveRestore){
    	    		gd = new GenericDialog("Select motor position restoration");
    	    		gd.addMessage("Motor position saved at the previous run of this program {"+savedPosition[0]+":"+savedPosition[1]+":"+savedPosition[2]+"]");
    	    		gd.addMessage("differs from the curent position reported by the hardware: {"+this.curpos[0]+":"+this.curpos[1]+":"+this.curpos[2]+"]");
    	    		gd.addMessage("");
    	    		gd.addMessage("You may restore the saved position, keep new coordinates (and overwrite saved data) or abort current command if you select \"Cancel\"");
    	    		gd.addMessage("");
    	    		gd.addMessage("If the camera was powered down between the program runs, as it "+(wasEnabled?"does not seem to be":"seems to be")+",");
    	    		gd.addMessage("\"Restore\" should be used, if it was not - powered down and the motors were just moved by another application,");
    	    		gd.addMessage("\"Use current\" is appropriate. Of course, it could be powered down and moved - in that case it would be difficult to restore");
    	    		gd.enableYesNoCancel("Restore", "Use current");
    	    	    gd.showDialog();
    	    	    if (gd.wasCanceled()) {
    	    	    	throw new RuntimeException("Command canceled");
    	    	    }
    	    	    if (gd.wasOKed()) {
        	        	motorHardwareInitializationInProgress=true;
    	    	    	restorePosition(savedPosition); // otherwise just keep current
    	    	    }
    			} else { // !interactiveRestore
    	        	motorHardwareInitializationInProgress=true;
    				restorePosition(savedPosition);
    			}
    		} else {
    			if (this.debugLevel>0) System.out.println("Motors did not move since last time this program ran");
    		}
    	} else { // no saved position
        	if (this.debugLevel>2) System.out.println("--savedPosition="+((savedPosition==null)?"null":("{"+savedPosition[0]+","+savedPosition[1]+","+savedPosition[2]+"}")));
			if (this.interactiveRestore){
	    		gd = new GenericDialog("Missing motor position file");
	    		gd.addMessage("No motor position file ("+this.getPrefsDir()+this.motorsStatePath+") found.");
	    		gd.addMessage("If you run the program for the first time it is normal, you may click \"OK\" and the new file will be created");
	    		gd.addMessage("If you select \"Cancel\" the current command will be aborted, you may check the file and restart the program.");
	    	    gd.showDialog();
	    	    if (gd.wasCanceled()) {
	    	    	throw new RuntimeException("Command canceled");
	    	    }
			}    		
    	}
    	motorHardwareInitializationInProgress=false;
    	this.stateValid=true;
    }
    private void restorePosition(int [] position){
//    	int [] inversePosition=new int [position.length];
//    	for (int i=0;i<inversePosition.length;i++) inversePosition[i]=-position[i];
    	if (this.debugLevel>2) System.out.println("restorePosition {"+position[0]+","+position[1]+","+position[2]+"}, motorHardwareInitializationInProgress="+motorHardwareInitializationInProgress);
    	
    	commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?m1="+position[0]+"&m2="+position[1]+"&m3="+position[2]+"&reset");
    	
/*    	
 *      Code for old FPGA firmware
    	int [] wasCurPos=this.curpos.clone();
    	moveElphel10364Motors( // checks again and opens dialog
    			true, //boolean wait,
    			inversePosition, //int [] positions,
    			0.0, // double sleep,
    			true, // boolean showStatus,
    			"Restoring inversed saved position", // String message,
    			false); // , boolean hysteresis)
    	if (this.debugLevel>0) System.out.println("*** RESETING MOTORS at position {"+this.curpos[0]+","+this.curpos[1]+","+this.curpos[2]+"}");
    	commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?reset"); // that would zero hardware motor registers
    	// Move back to the original location
    	moveElphel10364Motors(
    			true, //boolean wait,
    			wasCurPos, //int [] positions,
    			0.0, // double sleep,
    			true, // boolean showStatus,
    			"Restoring saved position", // String message,
    			true); // , boolean hysteresis)
    			*/
	}
    
    //if (maxPosDiff(positions)>this.positionTolerance){
    
    public int[] moveElphel10364Motors(boolean wait, int [] positions, double sleep, boolean showStatus, String message, boolean hysteresis)  {
    	// global this.motor3rd_signed_hysteresis
    	int i;
    	if (this.debugLevel>2) System.out.println("moveElphel10364Motors(), before checks: motorHardwareInitialized="+motorHardwareInitialized+
    			" motorHardwareInitializationInProgress="+motorHardwareInitializationInProgress);

    	checkEnabled(); // will either enable or throw;
    	if (this.debugLevel>2) System.out.println("moveElphel10364Motors(), after checkEnabled(): motorHardwareInitialized="+motorHardwareInitialized+
    			" motorHardwareInitializationInProgress="+motorHardwareInitializationInProgress);
    	checkSynchronized();
    	if (this.debugLevel>2) System.out.println("moveElphel10364Motors(), after checkSynchronized(): motorHardwareInitialized="+motorHardwareInitialized+
    			" motorHardwareInitializationInProgress="+motorHardwareInitializationInProgress);
    	readElphel10364Motors(); // sets this.curpos
    	if (hysteresis) {
    		int [] hystPositions=positions.clone();
    		boolean hystRequired=false;
    		for (i=0;i<3;i++) if ((positions[i]-this.curpos[i])*this.motor3rd_signed_hysteresis<0){
    			hystPositions[i]-=this.motor3rd_signed_hysteresis;
    			hystRequired=true;
    		}
    		if (hystRequired) {
    			moveElphel10364Motors(true, hystPositions, 0, showStatus, message+((message.length()>0)?": c":"C")+"ompensating motors hysteresis -");
    		}
    	}
    	if (this.debugLevel>2) System.out.println("moveElphel10364Motors(), before actual move: motorHardwareInitialized="+motorHardwareInitialized+
    			" motorHardwareInitializationInProgress="+motorHardwareInitializationInProgress);
    	return moveElphel10364Motors(wait, positions, sleep, showStatus, message);
    }

    private int[] moveElphel10364Motors(boolean wait, int [] positions, double sleep, boolean showStatus, String message)  {
    	IJ.showStatus(message+"Moving m1="+positions[0]+" m2="+positions[1]+" m3="+positions[2]);
    	if (this.debugLevel>0) System.out.println(message+((message.length()>0)?"":"Moving motors")+" m1="+positions[0]+", m2="+positions[1]+", m3="+positions[2]+" wait="+wait);
    	return moveElphel10364Motors(wait, positions, sleep);
    }

    private int[] moveElphel10364Motors(boolean wait, int [] positions, double sleep)  {
    	moveElphel10364Motors(wait, positions);
    	return readElphel10364Motors(sleep);
    }

    private int[] moveElphel10364Motors(boolean wait, int [] positions)  {
    	int extraCalm=4; // sec
    	if (OK2moveElphel10364Motors(positions) || motorHardwareInitializationInProgress) { // OK to override ~OK
    		this.lastMovedTo=positions.clone();
    		String uri="http://"+this.ipAddress+"/10364.php?m1="+positions[0]+"&m2="+positions[1]+"&m3="+positions[2];
    		if (wait) {
    			double travel=0.0;

    			for (int i=0;i<positions.length;i++){
    				if (positions[i]!=this.curpos[i]) this.lastDirs[i]=(positions[i]>this.curpos[i]);
    				double t=Math.abs(positions[i]-this.curpos[i]);
    				if (travel<t) travel=t;
    			}
    			int waitSeconds=(int) Math.ceil(travel*this.waitPerStep+this.waitExtra);
    			uri+="&wait="+waitSeconds;
    			commandElphel10364Motors(uri+"&calm="+this.calmMotors);
    			if (maxPosDiff(positions)>this.positionTolerance){
    				System.out.println("Motors ["+this.curpos[0]+","+this.curpos[1]+","+this.curpos[2]+"] did not reach destination "+
    						"["+positions[0]+","+positions[1]+","+positions[2]+"], travel was="+travel+", waitSeconds="+waitSeconds+
    						". Waiting "+extraCalm+" seconds more...");
    				commandElphel10364Motors(uri+"&calm="+extraCalm);
    				if (maxPosDiff(positions)>this.positionTolerance){
    					System.out.println("Motors ["+this.curpos[0]+","+this.curpos[1]+","+this.curpos[2]+"] did not reach destination "+
    							"["+positions[0]+","+positions[1]+","+positions[2]+"] after extra waiting");
    				}
    			}
    		} else {
    			commandElphel10364Motors(uri);
    			if (this.debugLevel>0) System.out.println("Warning: send move command with no wait ["+positions[0]+","+positions[1]+","+positions[2]+"]");

    		}
// save position after each movement    		
    		try {
				saveMotorState();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    		return this.curpos;
    	} else {
			if (this.debugLevel>0) System.out.println("Warning: OK2moveElphel10364Motors("+positions[0]+","+positions[1]+","+positions[2]+") returned FALSE!");
    	}
    	return readElphel10364Motors(); // just read
    }
    private int maxPosDiff(int [] positions){
    	int max=0;
    	for (int i=0;i< positions.length;i++){
    		if ((positions[i]-this.curpos[i])>max) max= positions[i]-this.curpos[i];
    		if ((this.curpos[i]-positions[i])>max) max= this.curpos[i]-positions[i];
    	}
    	return max;
    }

    private boolean OK2moveElphel10364Motors( int [] positions)  {
    	// public double [][]    this.motorsRange{{-5000,5000},{-5000,5000},{-5000,5000}};
    	return (positions[0]>=this.motorsRange[0][0]) && (positions[0]<= this.motorsRange[0][1]) &&
    	(positions[1]>=this.motorsRange[1][0]) && (positions[1]<= this.motorsRange[1][1]) &&
    	(positions[2]>=this.motorsRange[2][0]) && (positions[2]<= this.motorsRange[2][1]);
    }


    public int[] readElphel10364Motors()  {
    	return readElphel10364Motors(0.0);
    }

    public int[] readElphel10364Motors(double sleep)  {
//    	if (isEnabled)
//   	if (isInitialized)
//    	checkSynchronized(); // Do not check synchronization here - just read current position (OK to do even if motorst are disabled) 
    	return commandElphel10364Motors("http://"+this.ipAddress+"/10364.php?sleep="+sleep);
    }

    private int[] commandElphel10364Motors(String url)  {
    	Document dom=null;
    	if (this.debugLevel>2)	System.out.println("commandElphel10364Motors("+url+")");
    	try {
    		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
    		DocumentBuilder db = dbf.newDocumentBuilder();
    		dom = db.parse(url);
    		if (!dom.getDocumentElement().getNodeName().equals("motors")) {
    			IJ.showMessage("Root element: expected 'motors', got'" + dom.getDocumentElement().getNodeName()+"'"); 
    			return null;
    		}
    		this.curpos[0]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor1").item(0)).getChildNodes().item(0))).getNodeValue());
    		this.curpos[1]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor2").item(0)).getChildNodes().item(0))).getNodeValue());
    		this.curpos[2]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor3").item(0)).getChildNodes().item(0))).getNodeValue());
    	} catch(MalformedURLException e){
    		System.out.println("Please check the URL:" + e.toString() );
    		return null;
    	} catch(IOException  e1){
    		IJ.showStatus("");
    		String error = e1.getMessage();
    		if (error==null || error.equals(""))  error = ""+e1;
    		IJ.showMessage("commandElphel10364Motors ERRROR", ""+error);
    		return null;
    	}catch(ParserConfigurationException pce) {
    		pce.printStackTrace();
    		return null;
    	}catch(SAXException se) {
    		se.printStackTrace(); 
    		return null;
    	}
    	// prevent accumulating position errors 
    	if ((this.lastMovedTo!=null) && (Math.abs(this.lastMovedTo[0]-this.curpos[0])<=this.positionTolerance) &&
    			(Math.abs(this.lastMovedTo[1]-this.curpos[1])<=this.positionTolerance) &&
    			(Math.abs(this.lastMovedTo[2]-this.curpos[2])<=this.positionTolerance)) this.curpos=this.lastMovedTo.clone();  
    	//    	    	  this.lastMovedTo=positions.clone();
    	//	public int    positionTolerance=2; // actual position should not differ than this from the required

    	return this.curpos;
    }

    public void clearHistory(){
    	this.focusingHistory.clear();
    }

    public void addToHistory( String sTimestamp,double temperature, double[][] psfMetrics, double [][][][] fullResults){ // null OK
    	this.focusingHistory.add(
    			sTimestamp,
    			temperature,
    			this.curpos,
    			this.lastDirs,
    			psfMetrics,
    			fullResults
    	);
    }
    
	public FocusingField.FocusingFieldMeasurement getThisFFMeasurement(
			FocusingField focusingField,
			String sTimestamp,
			double temperature,
			double[][] psfMetrics,
			double [][][][] fullResults){
		
		return focusingField.getFocusingFieldMeasurement(
				sTimestamp,   //focusingState.getTimestamp(),
				temperature, //focusingState.getTemperature(),
				this.curpos,      //focusingState.motorsPos,
				fullResults); //focusingState.getSamples());
	}    
    
    
    
/*
 * FocusingHistory.FocusingState focusingState
    public void addToHistory( String sTimestamp,double temperature, double[][] psfMetrics){ // null OK
    	this.focusingHistory.add(
    			sTimestamp,
    			temperature,
    			this.curpos,
    			this.lastDirs,
    			psfMetrics,
    			null
    	);
    }
*/
    public void addToHistory(String sTimestamp,double temperature){
    	this.focusingHistory.add(
    			sTimestamp,
    			temperature,
    			this.curpos,
    			this.lastDirs,
    			null,
    			null);
    }

    public void setLastMetrics(double[][] psfMetrics){
    	this.focusingHistory.setLastMetrics(psfMetrics);
    }

    public void listHistory(
    		boolean useLMA,
    		String path,
    		String lensSerial,
    		String comment,
    		boolean showIndividualComponents,
    		boolean showIndividualSamples,
    		double weightRatioRedToGreen,
    		double weightRatioBlueToGreen,
			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
			double weightY // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
			){
    	listHistory(
    			useLMA,
    			path,
    			"",
    			lensSerial,
    			comment,
    			showIndividualComponents,
    			showIndividualSamples,
    			weightRatioRedToGreen,
    			weightRatioBlueToGreen,
    			weightK,
    			weightY,
    			false,
    			Double.NaN,
    			Double.NaN,
    			Double.NaN,
    			Double.NaN);
    }
    
    public void listHistory(
    		boolean useLMA,
    		String path,
    		String serialNumber,
    		String lensSerial,
    		String comment,
    		boolean showIndividualComponents,
    		boolean showIndividualSamples,
    		double weightRatioRedToGreen,
    		double weightRatioBlueToGreen,
			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
			boolean justSummary,
			double pX0,
			double pY0,
			double result_lastKT,
			double result_allHistoryKT
			){

    	this.focusingHistory.list(
        		useLMA,
    			path,
        		serialNumber,
    			lensSerial,
    			comment,
    			showIndividualComponents,
    			showIndividualSamples,
    			weightRatioRedToGreen,
    			weightRatioBlueToGreen,
    			weightK,
    			weightY,
    			justSummary,
    			pX0,
    			pY0,
    			result_lastKT,
    			result_allHistoryKT);
    }
	public void saveHistoryXML(
			String path,
    		String serialNumber,
			String lensSerial, // if null - do not add average
    		String comment,
			double pX0,
			double pY0,
			double [][][] sampleCoord){ // x,y,r
		this.focusingHistory.saveXML(
    			path,
        		serialNumber,
    			lensSerial, // if null - do not add average
        		comment,
    			pX0,
    			pY0,
    			sampleCoord);
	}

	public FocusingField.FocusingFieldMeasurement getThisFFMeasurement(FocusingField focusingField){
		return getThisFFMeasurement(focusingField, -1);
	}    
	public FocusingField.FocusingFieldMeasurement getThisFFMeasurement(FocusingField focusingField, int index){
		return focusingField. getFocusingFieldMeasurement(
				historyGetTimestamp(index),   //focusingState.getTimestamp(),
				historyGetTemperature(index), //focusingState.getTemperature(),
				historyGetMotors(index),      //focusingState.motorsPos,
				historyGetSamples(index)); //focusingState.getSamples());
	}
	private FocusingHistory.FocusingState getFocusingState(int index){
		if (index<0) index+=historySize();
		if ((index>=0) && (index<historySize())) return this.focusingHistory.history.get(index);
		return null;
	}
	
	public void addCurrentHistoryToFocusingField(
			FocusingField focusingField,
			int start,
			int end){
		for (int i=start;i<end;i++){
			addCurrentHistoryToFocusingField(focusingField,i);
		}        	 
	} 

	public void addCurrentHistoryToFocusingField(FocusingField focusingField){
		addCurrentHistoryToFocusingField(
				focusingField,
				0,
				historySize()
				);
	} 

	public void addCurrentHistoryToFocusingField(
			FocusingField focusingField,
			int index){ // -1 - last (negative - from length)
		focusingField.addSample(
				historyGetTimestamp(index),   //focusingState.getTimestamp(),
				historyGetTemperature(index), //focusingState.getTemperature(),
				historyGetMotors(index),      //focusingState.motorsPos,
				historyGetSamples(index)); //focusingState.getSamples());
	}
	
	
	public String historyGetTimestamp(int index){
		return getFocusingState(index).getTimestamp();
	}

	public double historyGetTemperature(int index){
		return getFocusingState(index).getTemperature();
	}
	public int [] historyGetMotors(int index){
		return getFocusingState(index).getMotors();
	}
	
	public double [][][][] historyGetSamples(int index){
		return getFocusingState(index).getSamples();
	}

	public int historySize(){
    	return this.focusingHistory.history.size();
    }
    public void setLastProbed(){
    	this.focusingHistory.setLastProbed();
    }
    public double distFromProbed(){ // may return NaN
    	return this.focusingHistory.distFromProbed();
    }
    public double distFromProbed(int [] position){ // may return NaN
    	return this.focusingHistory.distFromProbed(position);
    }
    /**
     * Asks for the motor positions, moves motors, waits for the movement complete
     * @param newPosition array of 3 integers or null (will prompt current position
     * @param useHisteresis 1-element array (will be modified) - apply or not anti-hysteresis moves
     * @param showAuto if true, will show button "Auto", if it will be pressed - no movement
     * @return 0 - cancel, 1 - OK (moved), 2 - pressed "Auto"
     */
    public int dialogMove(
    		int [] newPosition, // null OK
    		boolean [] useHysteresis,
    		boolean showAuto
    ){
    	readElphel10364Motors();
    	GenericDialog gd = new GenericDialog("Moving Focus Motors");
    	int [] promptPosition=(newPosition==null)?this.curpos.clone():newPosition.clone();
    	for (int i=0;i<promptPosition.length;i++){
    		gd.addNumericField("Motor "+(i+1)+" "+(this.lastDirs[i]?"(+)":"(-)")+
    				((newPosition==null)?"":(" ("+this.curpos[i]+")"))+":", promptPosition[i], 0,5,"steps");
    	}
    	if (useHysteresis!=null) gd.addCheckbox("Use anti-hysteresis moves", useHysteresis[0]);
    	if (showAuto) gd.enableYesNoCancel("OK", "Auto");
    	gd.showDialog();
    	if (gd.wasCanceled()) return 0;
    	for (int i=0;i<promptPosition.length;i++){
    		promptPosition[i]= (int) gd.getNextNumber();
    	}
    	if (useHysteresis!=null)  useHysteresis[0]=gd.getNextBoolean();
    	else {
    		useHysteresis=new boolean[1];
    		useHysteresis[0]=false;
    	}
    	if (!gd.wasOKed()) return 2;
    	moveElphel10364Motors(
    			true, //boolean wait,
    			promptPosition,
    			0.0, //double sleep,
    			true, //boolean showStatus,
    			"",   //String message,
    			useHysteresis[0]); //boolean hysteresis)
    	// Ask something else
    	return 1;
    }

    public int dialogFocusing(
    		boolean isAdjusted,
    		int [] newPosition, // null OK
    		int dialogMode,  // 1 - manual, 2 auto first, 3 - auto normal, 4 fine focus (first), 5 fine focus (auto)
    		LensAdjustment.FocusMeasurementParameters focusMeasurementParameters
    ){
    	readElphel10364Motors();
    	GenericDialog gd = new GenericDialog("Moving Focus Motors  - mode="+dialogMode);
    	gd.addMessage("History size - "+historySize());
    	if (historySize()>0){
    		FocusingHistory.FocusingState focusingState=this.focusingHistory.history.get(historySize()-1);
			double fDistance=this.focusingHistory.getLensDistance(
					focusingState.getCenterResolutions(), // {R-sharpness,G-sharpness,B-sharpness}
	    			true, // boolean absolute, // return absolutely calibrated data
	    			focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
	    			focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
	    			this.debugLevel); // motors debug level
        	if (!Double.isNaN(fDistance)){
        		gd.addMessage("Relative focal distance - "+IJ.d2s(fDistance,3)+" um");
        	}
    		if (historySize()>1){
    			FocusingHistory.FocusingState prevFocusingState=this.focusingHistory.history.get(historySize()-2);
    			double [] xyz={
    					focusingState.motorsPos[0]-prevFocusingState.motorsPos[0],
    					focusingState.motorsPos[1]-prevFocusingState.motorsPos[1],
    					focusingState.motorsPos[2]-prevFocusingState.motorsPos[2]};
    			gd.addMessage("Last travel "+IJ.d2s(Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]),0)+" steps");
    		}
    		double [][] metrics=focusingState.getMetrics(focusMeasurementParameters.weightRatioRedToGreen,focusMeasurementParameters.weightRatioBlueToGreen);
    		for (int i=0;i<metrics.length;i++) {
    			gd.addMessage(focusingState.getName(i)+
    					" Far/Near="+IJ.d2s(metrics[i][0],3)+
    					" TiltX="+IJ.d2s(metrics[i][1],3)+
    					" TiltY="+IJ.d2s(metrics[i][2],3)+
    					" R50%="+IJ.d2s(metrics[i][3],3)+
    					" A50%="+IJ.d2s(metrics[i][4],3)+
    					" B50%="+IJ.d2s(metrics[i][5],3)+
    					" R50%Center="+IJ.d2s(metrics[i][6],3)+
    					((i==3)?(" Weight="+IJ.d2s(100*metrics[i][7],3)+"%"):""));
    		}
    	}
//    	int [] promptPosition=(newPosition==null)?this.curpos.clone():newPosition.clone();
    	
    	for (int i=0;i<newPosition.length;i++){
    		gd.addNumericField("Motor "+(i+1)+" "+(this.lastDirs[i]?"(+)":"(-)")+
    				((newPosition==null)?"":(" ("+this.curpos[i]+")"))+":", newPosition[i], 0,5,"steps");
    	}
    	gd.addNumericField("Weight ratio  red/green (use when combining defocusing data)",                    focusMeasurementParameters.weightRatioRedToGreen,    2,6,"");
    	gd.addNumericField("Weight ratio blue/green (use when combining defocusing data)",                    focusMeasurementParameters.weightRatioBlueToGreen,    2,6,"");
    	gd.addNumericField("Far/Near Focusing target (logariphm of average tangential-to-radial resolution)", focusMeasurementParameters.targetFarNear,    2,6,"");
		gd.addCheckbox    ("Use targetFarNear (radial/tangential resolution)  as a proxy for the distance",   focusMeasurementParameters.useRadialTangential); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
		gd.addNumericField("Needed lens center distance (away from the \"best focus\"",                       focusMeasurementParameters.targetMicrons,    2,6,"um");
		gd.addNumericField("Tolerance of the center focal distance",                                          focusMeasurementParameters.toleranceMicrons,    2,6,"um");
		gd.addNumericField("Tolerance of the tilt ",                                                          focusMeasurementParameters.toleranceTilt,    2,6,"");
		gd.addNumericField("When under this scaled tolerance, correction step is scaled down twice",          focusMeasurementParameters.toleranceThreshold,    2,6,"");
		gd.addNumericField("Adjust 3 motors by parallel movement if focal distance error in the center exceeds this value", focusMeasurementParameters.parallelAdjustThreshold,    2,6,"um");
		
    	gd.addNumericField("When fitting far/near and tilts, samples 'decay' with this sigma",                focusMeasurementParameters.motorsSigma,      1,7,"motors steps");
		gd.addNumericField("Same for the 3 motors moving together (center focusing)",                         focusMeasurementParameters.motorsSigma3,      1,7,"motors steps");
    	gd.addNumericField("Low limit for motor sigma when fitting is in the final stage",                    focusMeasurementParameters.motorsMinSigma,      1,7,"motors steps");
    	gd.addNumericField("When walk is getting smaller, sigma be reduced proportionally",                   focusMeasurementParameters.motorsVarSigmaToTravel,      1,7,"motors steps");
    	gd.addNumericField("After each step new sigma will have this part of the calculated from the travel", focusMeasurementParameters.motorsFadeSigma,      1,7,"motors steps");
    	gd.addNumericField("For quadratic maximum the correction will be increased by (1+motorsOverShootToBalance) if there are less samples on the other side", focusMeasurementParameters.motorsOverShootToBalance,      1,7,"motors steps");
    	gd.addNumericField("Maximal allowed single-step focusing adjustment",                                   focusMeasurementParameters.maxStep,      1,7,"motors steps");
    	gd.addNumericField("How far to go to probe around the current point to measure derivatives",            focusMeasurementParameters.probeStep,      1,7,"motors steps");
    	gd.addCheckbox    ("Probe 6 measurements, if unchecked - only 4 (tetrahedron)",                      focusMeasurementParameters.probeSymmetrical);
    	gd.addNumericField("Re-run probing around in orthogonal directions, if the current position moved farther from the last probing one", focusMeasurementParameters.reProbeDistance, 1,7,"motors steps");
    	gd.addNumericField("Bias towards the last measurement: 0% use best fit for all, 100% make planes throgh the last point",100*focusMeasurementParameters.believeLast,      1,5,"%");
    	gd.addCheckbox    ("Move motors in the same direction to compensate for the hysteresis",             focusMeasurementParameters.compensateHysteresis);
		gd.addNumericField("Motor anti-hysteresis travel",                                                   focusMeasurementParameters.motorHysteresis, 0,7,"motors steps");

		gd.addNumericField("Sensor travel to motors travel (all 3 together), by design it is 4/38~=0.105",   focusMeasurementParameters.linearReductionRatio, 5,7,"");

		gd.addNumericField("Motors debug (1 - show moves, 2 show moves+hysteresis mioves)",                  focusMeasurementParameters.motorDebug,        0);

		gd.addMessage("Parameters to calculate lens distance measurements from 3-color PSF measurements");
		gd.addNumericField("Number of points to tabulate center focus parameters vs. focal distance",        focusMeasurementParameters.lensDistanceNumPoints,        0);
		gd.addNumericField("polynomial degree to approximate center focus parameters vs. focal distance ",   focusMeasurementParameters.lensDistancePolynomialDegree,        0);
		gd.addNumericField("Normalize overall sharpness (that depends on the lens quality and/or PSF parameters) to differential ones",            focusMeasurementParameters.lensDistanceWeightY,        3,5,"");
		gd.addNumericField("Use derivartives-dependent weights for components (1.0), equal weight - 0.0",    focusMeasurementParameters.lensDistanceWeightK,        3,5,"");
		gd.addMessage("");
    	
    	gd.addNumericField("Minimal correction movement to initiate final series of corrections",            focusMeasurementParameters.minCorr,      1,5,"motors steps");
    	gd.addNumericField("Finish if this number of last corrections where below  minimum (previous input)",focusMeasurementParameters.numFinalCorr,      0);
    	gd.addNumericField("Unconditionally exit focusing adjustment after these number of itereations",     focusMeasurementParameters.maxAutoIterations,        0);
    	gd.addNumericField("Maximal allowed total motors travel for automatic correction",                   focusMeasurementParameters.maxAutoDistance, 1,7,"motors steps");
    	gd.addNumericField("If there are insufficient measurements to fit parabola - make this step",        focusMeasurementParameters.maxLinearStep, 1,7,"motors steps");
    	gd.addCheckbox    ("Run fine-focus adjustment for the center samples",                               (dialogMode==4)||(dialogMode==4));

    	// TODO: move to other dialog 
    	gd.addCheckbox    ("Scan motors for the focus in the center samples",                                false);
    	gd.addNumericField("Motor single movement (all 3 motors) in scan focus mode (signed value)",         focusMeasurementParameters.scanStep, 0,7,"motors steps");
    	gd.addNumericField("Number of scan steps during (center) focus scanning",                            focusMeasurementParameters.scanNumber,        0);
		gd.addCheckbox    ("Scan focus in 2 directions, after the calibration estimate hysteresis (play)",   focusMeasurementParameters.scanHysteresis);
		gd.addNumericField("Number of scan steps during hysteresis (play) measurement",                      focusMeasurementParameters.scanHysteresisNumber, 0);

    	if      (dialogMode==1) gd.enableYesNoCancel("OK", "Auto Suggest");
    	else if (dialogMode==2) gd.enableYesNoCancel("Apply and Continue", "Retry");
    	else if (dialogMode==3) {
    		if (isAdjusted) gd.enableYesNoCancel("OK/Done", "Single Step"); // here OK - same as cancel?
    		else            gd.enableYesNoCancel("Continue", "Single Step");
    	}
    	else if (dialogMode==4) gd.enableYesNoCancel("Apply and Continue", "Retry");
    	else if (dialogMode==5) gd.enableYesNoCancel("Continue", "Single Step");
    	else {
    		String msg="Mode should be 1,2,3,4 or 5. Got "+dialogMode;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
		WindowTools.addScrollBars(gd);
    	gd.showDialog();
    	if (gd.wasCanceled()) return 0;
    	int [] wasPosition=newPosition.clone();
    	boolean wasModified=false;
    	for (int i=0;i<newPosition.length;i++){
    		newPosition[i]= (int) gd.getNextNumber();
    		wasModified|=(newPosition[i]!=wasPosition[i]);
    	}
    	focusMeasurementParameters.weightRatioRedToGreen=      gd.getNextNumber();
    	focusMeasurementParameters.weightRatioBlueToGreen=     gd.getNextNumber();
    	focusMeasurementParameters.targetFarNear=              gd.getNextNumber();
    	focusMeasurementParameters.useRadialTangential=        gd.getNextBoolean();  
    	focusMeasurementParameters.targetMicrons=              gd.getNextNumber();
    	focusMeasurementParameters.toleranceMicrons=           gd.getNextNumber();
    	focusMeasurementParameters.toleranceTilt=              gd.getNextNumber();
    	focusMeasurementParameters.toleranceThreshold=         gd.getNextNumber();
    	focusMeasurementParameters.parallelAdjustThreshold=    gd.getNextNumber();

    	focusMeasurementParameters.motorsSigma=                gd.getNextNumber();
    	focusMeasurementParameters.motorsSigma3=               gd.getNextNumber();
    	focusMeasurementParameters.motorsMinSigma=             gd.getNextNumber(); 
    	focusMeasurementParameters.motorsVarSigmaToTravel=     gd.getNextNumber();
    	focusMeasurementParameters.motorsFadeSigma=            gd.getNextNumber(); 
    	focusMeasurementParameters.motorsOverShootToBalance=   gd.getNextNumber();

    	focusMeasurementParameters.maxStep=                    gd.getNextNumber();
    	focusMeasurementParameters.probeStep=                  gd.getNextNumber();
    	focusMeasurementParameters.probeSymmetrical=           gd.getNextBoolean();
    	focusMeasurementParameters.reProbeDistance=            gd.getNextNumber();
    	focusMeasurementParameters.believeLast=           0.01*gd.getNextNumber(); if (focusMeasurementParameters.believeLast>1.0) focusMeasurementParameters.believeLast=1.0; else if (focusMeasurementParameters.believeLast<0) focusMeasurementParameters.believeLast=0.0;
    	focusMeasurementParameters.compensateHysteresis=       gd.getNextBoolean();
    	focusMeasurementParameters.motorHysteresis=      (int) gd.getNextNumber();
    	focusMeasurementParameters.linearReductionRatio=       gd.getNextNumber();

    	focusMeasurementParameters.motorDebug=           (int) gd.getNextNumber();
    	focusMeasurementParameters.lensDistanceNumPoints=(int) gd.getNextNumber();
    	focusMeasurementParameters.lensDistancePolynomialDegree=(int) gd.getNextNumber();
    	focusMeasurementParameters.lensDistanceWeightY=        gd.getNextNumber();
    	focusMeasurementParameters.lensDistanceWeightK=        gd.getNextNumber();

    	focusMeasurementParameters.minCorr=                    gd.getNextNumber();
    	focusMeasurementParameters.numFinalCorr=         (int) gd.getNextNumber();
    	focusMeasurementParameters.maxAutoIterations=    (int) gd.getNextNumber();
    	focusMeasurementParameters.maxAutoDistance=            gd.getNextNumber();
    	focusMeasurementParameters.maxLinearStep=              gd.getNextNumber();
    	boolean fineFocus=                                     gd.getNextBoolean();

    	boolean scanMode=                                      gd.getNextBoolean();
    	focusMeasurementParameters.scanStep=             (int) gd.getNextNumber();
    	focusMeasurementParameters.scanNumber=           (int) gd.getNextNumber();
    	focusMeasurementParameters.scanHysteresis=            gd.getNextBoolean();
    	focusMeasurementParameters.scanHysteresisNumber= (int) gd.getNextNumber();

    	setHysteresis(focusMeasurementParameters.motorHysteresis);
    	setDebug(focusMeasurementParameters.motorDebug);

    	if (fineFocus){
    		if (gd.wasOKed()) return 7; // run fine focus
    		else return 6;              // run fine focus once, suggest result   
    	}
    	if (scanMode) return 5; // regardless of the non-cancel buttons
    	if (gd.wasOKed()) {
    		if      (dialogMode==1) return 1;
    		else if (dialogMode==2) return 3;
    		else if (dialogMode==3) {
    			if (!isAdjusted) return 4;
    			else if (wasModified) return 1; // manual mode
    			else return 0; // all done
    		}
    		else if (dialogMode==4) return 4;
    		else if (dialogMode==5) return 4;
    	} else {
    		if      (dialogMode==1) return 2;
    		else if (dialogMode==2) return 2;
    		else if (dialogMode==3) return 2;
    	}
    	return 0; // should never get here    		  
    }
    public int dialogFocusingSharpness
    (
    		boolean isAdjusted,
    		int [] newPosition, // null is not OK (may be modified)
    		int dialogMode, // 1 - manual, 2 auto first, 3 - auto normal
    		LensAdjustment.FocusMeasurementParameters focusMeasurementParameters
    ){
    	readElphel10364Motors();
    	GenericDialog gd = new GenericDialog("Moving Focus Motors - mode="+dialogMode);
    	if (isAdjusted ){ // should be true: (focusingSharpness.size()>0)
    		int [] currentPos=focusingSharpness.getPosition(); // last
    		double offCenter= 0.001*getMicronsPerStep()*(currentPos[0]+currentPos[1]+currentPos[2])/3.0;
    		String msg1="Lens focus is approximately "+IJ.d2s(Math.abs(offCenter),3)+" mm "+
			((offCenter>0)?"closer to":"farther from")+" the lens than the motor zero point.";
    		String msg2="You may screw the lens "+((offCenter>0)?"OUT":"IN")+ " by "+IJ.d2s(Math.abs(2.0*offCenter),3)+" turns ("+
			IJ.d2s(Math.abs(360*2.0*offCenter),0)+" degrees)";
    		if (this.debugLevel>0)System.out.println(msg1+"\n"+msg2);
    		gd.addMessage(msg1);
    		gd.addMessage(msg2);
    		
    	}
    	gd.addMessage("History size - "+focusingSharpness.size());
    	if (focusingSharpness.size()>0){
    		double sharpness=focusingSharpness.getSharpness(); // last
    		int [] currentPos=focusingSharpness.getPosition(); // last
    		if (focusingSharpness.size()>1){
    			int [] prevPos=focusingSharpness.getPosition(focusingSharpness.size()-2); // pre-last
    			double [] xyz={
    					currentPos[0]-prevPos[0],
    					currentPos[1]-prevPos[1],
    					currentPos[2]-prevPos[2]};
    			//    				  gd.addMessage("Last travel "+IJ.d2s(Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]),0)+" steps");
    			gd.addMessage("Last center travel "+IJ.d2s((xyz[0]+xyz[1]+xyz[2])/3.0,1)+" steps");
    		}
    		gd.addMessage("Focus sharpness="+sharpness);
    	}
    	for (int i=0;i<newPosition.length;i++){
    		gd.addNumericField("Motor "+(i+1)+" "+(this.lastDirs[i]?"(+)":"(-)")+
    				" ("+this.curpos[i]+")"+":", newPosition[i], 0,5,"steps");
    	}
    	gd.addCheckbox    ("Move motors in the same direction to compensate for the hysteresis",             focusMeasurementParameters.compensateHysteresis);
		gd.addNumericField("Motor anti-hysteresis travel",                                                   focusMeasurementParameters.motorHysteresis, 0,7,"motors steps");
		gd.addNumericField("Motors debug (1 - show moves, 2 show moves+hysteresis mioves)",                  focusMeasurementParameters.motorDebug,        0);
    	gd.addNumericField("Minimal correction movement to initiate final series of corrections",            focusMeasurementParameters.minCorr,      1,5,"motors steps");
    	gd.addNumericField("Finish if this number of last corrections where below  minimum (previous input)",focusMeasurementParameters.numFinalCorr,      0);
    	gd.addNumericField("Unconditionally exit focusing adjustment after these number of itereations",     focusMeasurementParameters.maxAutoIterations,        0);
    	gd.addNumericField("Maximal allowed total motors travel for automatic correction",                   focusMeasurementParameters.maxAutoDistance, 1,7,"motors steps");
    	gd.addNumericField("when fitting parabola for focusing sharpness in the center, far measurements decay with this sigma", focusMeasurementParameters.motorsPreSigma, 1,7,"motors steps");
    	gd.addNumericField("If there are insufficient measurements to fit parabola - make this step",        focusMeasurementParameters.maxLinearStep, 1,7,"motors steps");

    	if      (dialogMode==1) gd.enableYesNoCancel("OK", "Auto Suggest");
    	else if (dialogMode==2) gd.enableYesNoCancel("Apply and Continue", "Retry");
    	else if (dialogMode==3){
    		if (isAdjusted) gd.enableYesNoCancel("OK/Done", "Single Step"); // here OK - same as cancel?
    		else            gd.enableYesNoCancel("Continue", "Single Step");
    	}
    	else {
    		String msg="Mode shou8ld be 1,2 or 3. Got "+dialogMode;
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
    	}
    	gd.showDialog();
    	if (gd.wasCanceled()) return -1;
    	int [] wasPosition=newPosition.clone();
    	boolean wasModified=false;
    	for (int i=0;i<newPosition.length;i++){
    		newPosition[i]= (int) gd.getNextNumber();
    		wasModified|=(newPosition[i]!=wasPosition[i]);
    	}
    	
    	focusMeasurementParameters.compensateHysteresis=       gd.getNextBoolean();
    	focusMeasurementParameters.motorHysteresis=      (int) gd.getNextNumber();
    	focusMeasurementParameters.motorDebug=           (int) gd.getNextNumber();
    	focusMeasurementParameters.minCorr=                    gd.getNextNumber();
    	focusMeasurementParameters.numFinalCorr=         (int) gd.getNextNumber();
    	focusMeasurementParameters.maxAutoIterations=    (int) gd.getNextNumber();
    	focusMeasurementParameters.maxAutoDistance=            gd.getNextNumber();
    	focusMeasurementParameters.motorsPreSigma=             gd.getNextNumber();
    	focusMeasurementParameters.maxLinearStep=              gd.getNextNumber();
    	setHysteresis(focusMeasurementParameters.motorHysteresis);
    	setDebug(focusMeasurementParameters.motorDebug);

    	if (gd.wasOKed()) {
    		if      (dialogMode==1) return 1;
    		else if (dialogMode==2) return 3;
    		else if (dialogMode==3) {
    			if (!isAdjusted) return 4;
    			else if (wasModified) return 1; // manual mode
    			else return 0; // all done
    		}
    	} else {
    		if      (dialogMode==1) return 2;
    		else if (dialogMode==2) return 2;
    		else if (dialogMode==3) return 2;
    	}
    	return 0; // should never get here    		  
    }

    public void clearPreFocus(){
    	this.focusingSharpness.clear();
    }
    public void listPreFocus(String path, String comment){ // null - on screen
    	this.focusingSharpness.list(path,comment);
    }
    //   				  ,

    public void addPreFocus(
    		double sharpness){
    	addPreFocus(this.curpos, sharpness);
    }

    public void addPreFocus(
    		int [] motorsPos,
    		double sharpness){
    	this.focusingSharpness.add(motorsPos,sharpness);
    }
    public int sizePreFocus(){
    	return this.focusingSharpness.size();
    }
    public int sizeFocus(){
    	return this.focusingHistory.size();
    }
    public int [] solvePreFocus(
    		double motorsSigma,
    		double maxStep,        // Maximal motor move for each step
    		double linearStep,     // when no quadratic maximum
			double motorsOverShootToBalance,
    		int debugLevel
    ){
    	return focusingSharpness.solve(
    			getRange(),
    			motorsSigma,
    			maxStep,        // Maximal motor move for each step
    			linearStep,     // when no quadratic maximum
    			motorsOverShootToBalance,
    			debugLevel
    	);

    }

    public int [] solveSinglePoly(
    		boolean first,
    		double wR,  // weight of red in focus sharpness
    		double wG,
    		double wB,
    		boolean sameTiltOnly, // us only history where the difference between motors was the same as in the last sample
    		boolean centerOnly, // only use center samples
    		int polyDegree,      // polynomial order (just 2 always)
    		double motorsSigma,
    		double maxStep,        // Maximal motor move for each step
    		double linearStep,     // when no quadratic maximum
    		double motorsMinSigma, 
    		double motorsVarSigmaToTravel,
    		double motorsFadeSigma, 
    		double motorsOverShootToBalance,
    		int debugLevel
    ){
    	return focusingHistory.solveSinglePoly(
    			first,
    			wR,  // weight of red in focus sharpness
    			wG,
    			wB,
    			sameTiltOnly, // us only history where the difference between motors was the same as in the last sample
    			centerOnly, // only use center samples
    			polyDegree,      // polynomial order (just 2 always)
    			getRange(),
    			motorsSigma,
    			maxStep,        // Maximal motor move for each step
    			linearStep,     // when no quadratic maximum
    			motorsMinSigma, 
    			motorsVarSigmaToTravel,
    			motorsFadeSigma, 
    			motorsOverShootToBalance,
    			debugLevel
    	);

    }


    public class FocusingSharpness{
    	public ArrayList<FocusingSharpnessItem> history=null;
    	public FocusingSharpness(int size){
    		history=new ArrayList <FocusingSharpnessItem>(size);
    	}
    	public void clear(){
    		history.clear();
    	}
    	public int [] getPosition(){
    		if (this.history.size()<1) return null;
    		return this.history.get(this.history.size()-1).motorsPos;
    	}
    	public int [] getPosition(int index){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).motorsPos;
    	}
    	public int size(){
    		return this.history.size();
    	}
    	public double getSharpness(){
    		if (this.history.size()<1) return Double.NaN;
    		return this.history.get(this.history.size()-1).sharpness;
    	}
    	public double getSharpness(int index){
    		if ((index<0) || (index>=this.history.size())) return Double.NaN;
    		return this.history.get(index).sharpness;
    	}
    	public void add(
    			int [] motorsPos,
    			double sharpness
    	){
    		this.history.add(new FocusingSharpnessItem(motorsPos,sharpness));
    	}

    	// moves motors all together, using average as a single parameter.
    	// tries a second-degree polynomial, if not enough data - linear (maxSteop)
    	// TODO: modify to use PolunomialApproximation class    		  
    	public int [] solve(
    			double [][] motorRange,
    			double motorsSigma,
    			double maxStep,        // Maximal motor move for each step
    			double linearStep,     // when no quadratic maximum
    			double motorsOverShootToBalance,
    			int debugLevel
    	){
    		double S0=0.0,SX=0.0,SX2=0.0,SX3=0.0,SX4=0.0,SF=0.0,SFX=0.0,SFX2=0.0;
    		double rangeCenter=(
    				motorRange[0][0]+motorRange[0][1]+
    				motorRange[1][0]+motorRange[1][1]+
    				motorRange[2][0]+motorRange[2][1])/6.0;
    		int [] origin=getPosition();
    		if (
    				(origin[0]<=motorRange[0][0]) ||
    				(origin[0]>=motorRange[0][1]) ||
    				(origin[1]<=motorRange[1][0]) ||
    				(origin[1]>=motorRange[1][1]) ||
    				(origin[2]<=motorRange[2][0]) ||
    				(origin[2]>=motorRange[2][1])) {
    			String message="Motors are off range: ["+origin[0]+", "+origin[1]+", "+origin[2]+"], range is ["+
    			motorRange[0][0]+".."+motorRange[0][1]+", "+
    			motorRange[1][0]+".."+motorRange[1][1]+", "+
    			motorRange[2][0]+".."+motorRange[2][1]+"]";
    			System.out.println ("*** ERROR *** "+message);
    			IJ.showMessage("ERROR",message);
    			return null;
    		}
    		double x0=(origin[0]+origin[1]+origin[2])/3.0;
    		double move=0;
    		if (size()==1){
    			move=(x0>rangeCenter)?-linearStep:linearStep;
    		} else {

    			double kExp=0.5/motorsSigma/motorsSigma;
    			for (int i=0;i<size();i++){
    				int [] position=getPosition(i);
    				double f=getSharpness(i);
    				double x=((position[0]+position[1]+position[2])/3.0-x0);
    				double w=Math.exp(-kExp*(x*x));
    				double x2=x*x;
    				double x3=x2*x;
    				S0+=  w;
    				SX+=  w*x;
    				SX2+= w*x2;
    				SX3+= w*x3;
    				SX4+= w*x3*x;
    				SF+=  w*f;
    				SFX+= w*f*x;
    				SFX2+=w*f*x2;
    			}
    			Matrix M,B;
    			if (size()<3){
    				double [][] aM={
    						{SX2,SX },
    						{SX, S0 }};
    				double [][] aB={
    						{SFX},
    						{SF}};
    				M=new Matrix(aM);
    				B=new Matrix(aB);
    			} else {
    				double [][] aM={
    						{SX4,SX3,SX2},
    						{SX3,SX2,SX },
    						{SX2,SX, S0 }};
    				double [][] aB={
    						{SFX2},
    						{SFX},
    						{SF}};
    				M=new Matrix(aM);
    				B=new Matrix(aB);

    				if (!(new LUDecomposition(M)).isNonsingular()){
    					if (debugLevel>0){
    						System.out.println("Singular matrix - using linear approximation");
    					}
    					double [][] aM2={
    							{SX2,SX },
    							{SX, S0 }};
    					double [][] aB2={
    							{SFX},
    							{SF}};
    					M=new Matrix(aM2);
    					B=new Matrix(aB2);
    				}

    			}
    			if (debugLevel>1){
    				System.out.println("M:");
    				M.print(10, 5);
    				System.out.println("B:");
    				B.print(10, 5);
    			}
    			double [] abc= M.solve(B).getColumnPackedCopy();
    			if (abc.length<3){
    				double [] tmp={0.0,abc[0],abc[1]};
    				abc=tmp;
    			}
    			if (debugLevel>1){
    				System.out.println("A:");
    				M.solve(B).print(10, 5);
    			}
    			move=0.0;
    			if (abc[0]>=0.0){ // no maximum at all
    				double pos=abc[0]*linearStep*linearStep+abc[1]*linearStep+abc[2];
    				double neg=abc[0]*linearStep*linearStep-abc[1]*linearStep+abc[2];
    				move=((pos>neg)?1:-1)*linearStep;
    			} else {
    				move=-abc[1]/(2*abc[0]);
//    				System.out.println("move="+move+" SX="+SX);
					if (move*SX<0){ // only balance in one direction and if center is "behind"
						double balancedMove=move-(SX-move)*S0;
						if (debugLevel>0) {
							System.out.println("Balancing (first estimation): move="+move+" balancedMove="+balancedMove);
						}
						if ((balancedMove/move)>motorsOverShootToBalance) balancedMove=move*(1.0+motorsOverShootToBalance);
						move= (int) Math.round(balancedMove);
						if (debugLevel>0) {
							System.out.println("Balancing limited move="+move);
						}
					}
    				if (Math.abs(move)>maxStep) move*=maxStep/Math.abs(move);
    				
    			}
    		}
    		// limit to the range   			
    		for (int m=0; m<motorRange.length;m++) {
    			if (move>0){
    				if ((origin[m]+move)>motorRange[m][1]) move=motorRange[m][1]-origin[m];
    			} else {
    				if ((origin[m]+move)<motorRange[m][0]) move=motorRange[m][0]-origin[m];
    			}
    		}
    		int [] result = {
    				origin[0]+(int) move,
    				origin[1]+(int) move,
    				origin[2]+(int) move};
    		// see if it is an old result - then double the step?
    		boolean newPoint=true;
    		for (int i=0;i<size();i++){
    			int [] position=getPosition(i);
    			if ((position[0]==result[0])&&(position[1]==result[1])&&(position[2]==result[2])){
    				newPoint=false;
    				break;
    			}
    		}			
    		if (!newPoint) {
    			if (debugLevel>0) {
    				System.out.println("Got to the same point: ["+result[0]+", "+result[1]+", "+result[2]+"] again, doubling the step");
    			}
    			for (int m=0; m<motorRange.length;m++) {
    				if (move>0){
    					if ((result[m]+move)>motorRange[m][1]) move=motorRange[m][1]-result[m];
    				} else {
    					if ((result[m]+move)<motorRange[m][0]) move=motorRange[m][0]-result[m];
    				}
    			}
    			for (int m=0; m<motorRange.length;m++) result[m]+=move;
    		}
    		return result;
    	}


    	public void list(String path, String comment){
    		String header="#\tcenter travel\tCenter\tMotor1\tMotor2\tMotor3\tSharpness";
    		StringBuffer sb = new StringBuffer();
    		int [] prevPos=null;
    		for (int i=0;i<size();i++){
    			int [] motors=getPosition(i);
    			double sharpness=getSharpness(i);
    			sb.append((i+1)+"\t");
    			if (prevPos!=null) {
    				double [] xyz={motors[0]-prevPos[0],motors[1]-prevPos[1],motors[2]-prevPos[2]};
    				//	    	    		sb.append(""+IJ.d2s(Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]),0));
    				sb.append(""+IJ.d2s((xyz[0]+xyz[1]+xyz[2])/3.0,0)); // signed
    			}
    			sb.append("\t"+IJ.d2s((motors[0]+motors[1]+motors[2])/3.0,1)+"\t"+motors[0]+"\t"+motors[1]+"\t"+motors[2]+"\t"+IJ.d2s(sharpness,3)+"\n");
    			prevPos=motors;
    		}
    		if (path!=null) {
    			String footer="Comment: "+comment;
    			CalibrationFileManagement.saveStringToFile (
						path,
						header+"\n"+sb.toString()+"\n"+footer);
    		}else{
    			new TextWindow("Sharpness History", header, sb.toString(), 500,900);
    		}
    	}



    	public class FocusingSharpnessItem{
    		public int [] motorsPos=new int[3];
    		public double sharpness;
    		FocusingSharpnessItem(
    				int [] motorsPos,
    				double sharpness
    		){
    			this.motorsPos=motorsPos.clone();
    			this.sharpness=sharpness;
    		}
    	}
    }

    public class FocusingHistory{
    	private ArrayList<FocusingState> history=new ArrayList <FocusingState>(100);
    	private FocusingState focusingState = new FocusingState(); // just to access constant arrays
    	private double currentSigma=Double.NaN; // saves floating motor sigma between calls
    	private double [][]lensDistanceCalibration=null; // array to convert 3 color sharpnesses in the lens center to the distance
    	private double [][]polyCoeff=null; // polynomial approximations of the center color fractions and intensity
    	private double corrXMin;  // range of average motor coordinates to consider (should have the same tilt as last in history)
    	private double corrXMax;
    	private double micronsPerStep=Double.NaN;
    	private int    calibrationNumPoints; // number of calibration points in the table approximation;

    	private int [] calibrationHistoryFromTo=null; // history indices, to-from (inclusive) used to calculate calibration 
    	private int [] hysteresisHistoryFromTo=null; // history indices, to-from (inclusive) used to calculate hysteresis
    	private double measuredHysteresis=Double.NaN;
	    public int indexR=0; // index for Y data
	    public int indexB=1;
	    public int indexY=2;
	    public int indexRW=3;
	    public int indexBW=4;
	    public int indexYW=5;
	    public int indexSteps=6;
	    public int indexMicrons=7;

    	public FocusingHistory(int len){
    		this.history=new ArrayList <FocusingState>(len);
    	}
    	public double getMeasuredHysteresis(){ return this.measuredHysteresis;}
    	public void setCalibrationHistoryFromTo(int from, int to){
    		this.calibrationHistoryFromTo=new int [2];
    		this.calibrationHistoryFromTo[0]=from;
    		this.calibrationHistoryFromTo[1]=to;
    	}
    	
    	public void setHysteresisHistoryFromTo(int from, int to){
    		this.hysteresisHistoryFromTo=new int [2];
    		this.hysteresisHistoryFromTo[0]=from;
    		this.hysteresisHistoryFromTo[1]=to;
    	}
    	
    	public int [] getPosition(){
    		if (this.history.size()<1) return null;
    		return this.history.get(this.history.size()-1).getMotors();
    	}
    	public int [] getPosition(int index){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).getMotors();
    	}
    	public String getTimestamp(){
    		if (this.history.size()<1) return null;
    		return this.history.get(this.history.size()-1).getTimestamp();
    	}
    	public String getTimestamp(int index){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).getTimestamp();
    	}
    	
    	
    	public double [][]  getMetrics(double wR, double wG, double wB){
    		if (this.history.size()<1) return null;
    		return getMetrics(this.history.size()-1,wR,wG,wB);
    	}
    	public double [][]  getMetrics(double wRtoG, double wBtoG){
    		if (this.history.size()<1) return null;
    		return getMetrics(this.history.size()-1,wRtoG,wBtoG);
    	}

    	public double [][]  getMetrics(int index, double wR, double wG, double wB){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).getMetrics(wR,wG,wB);
    	}
    	public double [][]  getMetrics(int index, double wRtoG, double wBtoG){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).getMetrics(wRtoG,wBtoG);
    	}

    	public double getPositionAverage(){
    		if (this.history.size()<1) return Double.NaN;
    		return this.history.get(this.history.size()-1).getMotorsAverage();
    	}
    	public double  getPositionAverage(int index){
    		if ((index<0) || (index>=this.history.size())) return Double.NaN;
    		return this.history.get(index).getMotorsAverage();
    	}
    	public double []  getCenterResolutions(){
    		if (this.history.size()<1) return null;
    		return getCenterResolutions(this.history.size()-1);
    	}
    	public double [] getCenterResolutions(int index){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).getCenterResolutions();
    	}
    	public double [] getzTxTy(int index){
    		if ((index<0) || (index>=this.history.size())) return null;
    		return this.history.get(index).getzTxTy();
    	}
    	public double [] getzTxTy(){
    		if (this.history.size()<1) return null;
    		return getzTxTy(this.history.size()-1);
    	}
    	public int size(){
    		return this.history.size();
    	}
    	public int [] solveSinglePoly(
    			boolean first,
    			double wR,  // weight of red in focus sharpness
    			double wG,
    			double wB,
    			boolean sameTiltOnly, // us only history where the difference between motors was the same as in the last sample
    			boolean centerOnly, // only use center samples
    			int polyDegree,      // polynomial order (just 2 always)
    			double [][] motorRange,
    			double motorsSigma,
    			double maxStep,        // Maximal motor move for each step
    			double linearStep,     // when no quadratic maximum
    			double motorsMinSigma, 
    			double motorsVarSigmaToTravel,
    			double motorsFadeSigma, 
    			double motorsOverShootToBalance,
    			int debugLevel
    	){
    		if (first || Double.isNaN(this.currentSigma)) this.currentSigma=motorsSigma;
    		if (polyDegree>2)polyDegree=2; // for now limit
    		int indexR50=centerOnly?6:3; // 6 is average R50 in the center samples,3 - in all samples
    		double rangeCenter=(
    				motorRange[0][0]+motorRange[0][1]+
    				motorRange[1][0]+motorRange[1][1]+
    				motorRange[2][0]+motorRange[2][1])/6.0;
    		int [] origin=getPosition();
    		if (
    				(origin[0]<=motorRange[0][0]) ||
    				(origin[0]>=motorRange[0][1]) ||
    				(origin[1]<=motorRange[1][0]) ||
    				(origin[1]>=motorRange[1][1]) ||
    				(origin[2]<=motorRange[2][0]) ||
    				(origin[2]>=motorRange[2][1])) {
    			String message="Motors are off range: ["+origin[0]+", "+origin[1]+", "+origin[2]+"], range is ["+
    			motorRange[0][0]+".."+motorRange[0][1]+", "+
    			motorRange[1][0]+".."+motorRange[1][1]+", "+
    			motorRange[2][0]+".."+motorRange[2][1]+"]";
    			System.out.println ("*** ERROR *** "+message);
    			IJ.showMessage("ERROR",message);
    			return null;
    		}
    		boolean []sameTiltmask=new boolean[size()];
    		int sameTiltSize=0;
    		for (int i=0;i<size();i++){
    			if (!sameTiltOnly) {
    				sameTiltmask[i]=true;
    				sameTiltSize++;
    			} else {
    				int [] position=getPosition(i);
    				if (((position[0]-position[1])==(origin[0]-origin[1])) && ((position[2]-position[1])==(origin[2]-origin[1]))){
    					sameTiltmask[i]=true;
    					sameTiltSize++;
    				} else sameTiltmask[i]=false;
    			}
    		}
    		if (debugLevel>1) System.out.println("solveSinglePoly(): sameTiltSize="+sameTiltSize+", size()="+size());
    		double x0=(origin[0]+origin[1]+origin[2])/3.0;
    		int move=0;
    		double w=wR+wG+wB;
    		wR/=w;
    		wG/=w;
    		wB/=w;
    		int iLinearStep=(int) linearStep;
    		if (sameTiltSize==1){
    			move= (x0>rangeCenter)?-iLinearStep:iLinearStep;
    		} else {
    			// create a data set for the polynomial approximation
    			double [][] data=new double [size()][3]; //x,f,w
    			double kExp=0.5/this.currentSigma/this.currentSigma;
    			//sameDiffOnly
    			double SX=0,S0=0;
    			for (int i=0;i<size();i++) if (sameTiltmask[i]){
    				data[i][0]=getPositionAverage(i)- x0;
    				double [][] metrics= this.history.get(i).getMetrics(wR,wG,wB); // here we need to average 1/value, so have to calculate
    				data[i][1]=wR/metrics[0][indexR50]+ wG/metrics[1][indexR50]+wB/metrics[2][indexR50];// "sharpness"
    				data[i][2]=Math.exp(-kExp*(data [i][0]*data [i][0]));
    				SX+=data[i][0]*data[i][2];
    				S0+=data[i][2];
    			} else {
    				data[i][2]=0.0;
    			}
    			if (S0>0) SX/=S0; // center of gravity of the weighted samples used for parabola fitting
    			// the following line did not cause PolynomialApproximation class to recompile after it was modified

    			//        				double [] polyCoeff=(new PolynomialApproximation(debugLevel)).polynomialApproximation1d(data, polyDegree);
    			PolynomialApproximation pa= new PolynomialApproximation(debugLevel);
    			double [] polyCoeff=pa.polynomialApproximation1d(data, polyDegree);

    			move=0;
    			// history can consist of single point repeated several times
    			if ((polyCoeff[2]==0.0) && (polyCoeff[1]==0.0)){
    				move=(x0>rangeCenter)?-iLinearStep:iLinearStep;
    			} else {
    				if (polyCoeff[2]>=0.0){ // no maximum at all

    					double pos=polyCoeff[2]*linearStep*linearStep+polyCoeff[1]*linearStep+polyCoeff[0];
    					double neg=polyCoeff[2]*linearStep*linearStep-polyCoeff[1]*linearStep+polyCoeff[0];
    					move=((pos>neg)?1:-1)*iLinearStep;
    				} else {
    					move=(int) (-polyCoeff[1]/(2*polyCoeff[2]));
    					// modify move to balance around maximum when looking for the quadratic maximum        						
    					if (move*SX<0){ // only balance in one direction and if center is "behind"
    						double balancedMove=move-(SX-move)*S0;
    						if (debugLevel>0) {
    							System.out.println("Balancing (first estimation): move="+move+" balancedMove="+balancedMove);
    						}
    						if ((balancedMove/move)>motorsOverShootToBalance) balancedMove=move*(1.0+motorsOverShootToBalance);
    						move= (int) Math.round(balancedMove);
    						if (debugLevel>0) {
    							System.out.println("Balancing limited move="+move);
    						}
    					}
    					if (Math.abs(move)>maxStep) move*=maxStep/Math.abs(move);
    				}
    			}
    		} 
    		// limit to the range   			
    		for (int m=0; m<motorRange.length;m++) {
    			if (move>0){
    				if ((origin[m]+move)>motorRange[m][1]) move=(int)motorRange[m][1]-origin[m];
    			} else {
    				if ((origin[m]+move)<motorRange[m][0]) move=(int)motorRange[m][0]-origin[m];
    			}
    		}
    		int [] result = {
    				origin[0]+ move,
    				origin[1]+ move,
    				origin[2]+ move};
    		// see if it is an old result - then double the step?
    		boolean newPoint=true;
    		for (int i=0;i<size();i++){
    			int [] position=getPosition(i);
    			if ((position[0]==result[0])&&(position[1]==result[1])&&(position[2]==result[2])){
    				newPoint=false;
    				break;
    			}
    		}			
    		if (!newPoint) {
    			int move0=move;
    			if (debugLevel>0) {
    				System.out.println("Got to the same point: ["+result[0]+", "+result[1]+", "+result[2]+"] again, doubling the step");
    			}
    			for (int m=0; m<motorRange.length;m++) {
    				if (move>0){
    					if ((result[m]+move)>motorRange[m][1]) move=(int) motorRange[m][1]-result[m];
    				} else {
    					if ((result[m]+move)<motorRange[m][0]) move=(int) motorRange[m][0]-result[m];
    				}
    			}
    			for (int m=0; m<motorRange.length;m++) result[m]+=move;
    			move+=move0; // full movement;
    		}
    		// modify sigma
    		this.currentSigma=(1.0-motorsFadeSigma)*this.currentSigma+motorsFadeSigma*(motorsVarSigmaToTravel*Math.abs(move));
    		if      (this.currentSigma>motorsSigma)    this.currentSigma=motorsSigma;
    		else if (this.currentSigma<motorsMinSigma) this.currentSigma=motorsMinSigma;
    		if (debugLevel>0) {
    			System.out.println("Current motor sigma is "+IJ.d2s(this.currentSigma,2));
    		}
    		return result;
    	}

    	public void setProperties(String prefix,Properties properties){
    		if (this.polyCoeff!=null){
    			properties.setProperty(prefix+"corrXMin",this.corrXMin+"");
    			properties.setProperty(prefix+"corrXMax",this.corrXMax+"");
    			properties.setProperty(prefix+"polyNumChannels",(this.polyCoeff.length)+"");
    			properties.setProperty(prefix+"polyDegree",(this.polyCoeff[0].length-1)+"");
    			for (int i=0;i<this.polyCoeff.length;i++) for (int j=0;j< (this.polyCoeff[0].length-1);j++)
        			properties.setProperty(prefix+"polyCoeff_"+i+"_"+j,this.polyCoeff[i][j]+"");
    		}
    	}

		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"corrXMin")!=null)
				this.corrXMin=Double.parseDouble(properties.getProperty(prefix+"corrXMin"));
			if (properties.getProperty(prefix+"corrXMax")!=null)
				this.corrXMax=Double.parseDouble(properties.getProperty(prefix+"corrXMax"));
			if ((properties.getProperty(prefix+"polyNumChannels")!=null) && (properties.getProperty(prefix+"polyDegree")!=null)){
				this.polyCoeff=new double[Integer.parseInt(properties.getProperty(prefix+"polyNumChannels"))]
				                          [Integer.parseInt(properties.getProperty(prefix+"polyDegree"))+1];
    			for (int i=0;i<this.polyCoeff.length;i++) for (int j=0;j< (this.polyCoeff[0].length-1);j++){
    				if (properties.getProperty(prefix+"polyCoeff_"+i+"_"+j)==null){
    					String msg="Missing parameter \"polyCoeff_"+i+"_"+j+" in configuration file";
    					System.out.println("Error: "+msg);
    					IJ.showMessage("Error", msg);
    				} else {
    					this.polyCoeff[i][j]=Double.parseDouble(properties.getProperty(prefix+"polyCoeff_"+i+"_"+j));
    				}
    			}
			}

			// note: needs absolute calibration to convert steps to microns: optimalMotorPosition(
			// LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
			// double micronsPerStep,
			//  int debugLevel);
	}  
		public int [] optimalMotorPosition(
				LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
				double micronsPerStep,
				int debugLevel
		) {	
			if (this.polyCoeff==null) return null; // not calibrated at all
    		this.calibrationNumPoints=focusMeasurementParameters.lensDistanceNumPoints;
           	makeCalibrationTable(
        			this.polyCoeff,
        			this.calibrationNumPoints,
    				this.corrXMin,  // range of average motor coordinates to consider (should have the same tilt as last in history)
    				this.corrXMax,
        			debugLevel
        			);
			return absoluteCalibrate(
					micronsPerStep,
					focusMeasurementParameters.targetMicrons,
					focusMeasurementParameters.lensDistanceWeightY, 
					focusMeasurementParameters.lensDistanceWeightK,
					debugLevel);
		}

    	public int [] calibrateLensDistance(
    			String path, // or null
        		LensAdjustment.FocusMeasurementParameters focusMeasurementParameters,
    			boolean centerOnly,
    			double micronsPerStep,
    			  int debugLevel
    			){
    		if (this.calibrationHistoryFromTo==null){
    			return null; // not calibrated
    		}
    		if (focusMeasurementParameters.lensDistanceInteractive){
    			GenericDialog gd = new GenericDialog("Lens distance calculation");
    			gd.addNumericField("Number of division points",     focusMeasurementParameters.lensDistanceNumPoints, 0); // 1000
    			gd.addNumericField("Polynomial degree",                focusMeasurementParameters.lensDistancePolynomialDegree, 0); // 8
    			gd.addCheckbox    ("Center samples only",              centerOnly); // true
    			gd.addNumericField("Weight Y",                           focusMeasurementParameters.lensDistanceWeightY, 3,5,""); // 0.5 
    			gd.addNumericField("Weight K",                           focusMeasurementParameters.lensDistanceWeightK, 3,5,""); // 0.5
	    		gd.addCheckbox    ("Show results window from focal distance calibration", focusMeasurementParameters.lensDistanceShowResults);
	    		gd.addNumericField("Needed lens center distance (away from the \"best focus\"", focusMeasurementParameters.targetMicrons,    2,6,"um");

    			gd.addNumericField("Debug level",                       debugLevel,0);
    			gd.showDialog();
    			if (gd.wasCanceled()) return null;
    			focusMeasurementParameters.lensDistanceNumPoints=         (int) gd.getNextNumber();
    			focusMeasurementParameters.lensDistancePolynomialDegree=  (int) gd.getNextNumber();
    			centerOnly=                                                     gd.getNextBoolean();
    			focusMeasurementParameters.lensDistanceWeightY=                 gd.getNextNumber();
    			focusMeasurementParameters.lensDistanceWeightK=                 gd.getNextNumber();
    			focusMeasurementParameters.lensDistanceShowResults=             gd.getNextBoolean();
    			focusMeasurementParameters.targetMicrons=                       gd.getNextNumber();
    			debugLevel=                                               (int) gd.getNextNumber();
    		}
 //   		int sameTiltSize=0;
    		if (this.calibrationHistoryFromTo[0]<0)this.calibrationHistoryFromTo[0]=0;
    		if (this.calibrationHistoryFromTo[1]>=size())this.calibrationHistoryFromTo[1]=size()-1;
    		int [] origin=getPosition(this.calibrationHistoryFromTo[1]);
    		boolean []sameTiltmask=new boolean[size()];
    		for (int i=0;i<sameTiltmask.length;i++) sameTiltmask[i]=false;
    		for (int i=this.calibrationHistoryFromTo[0];i<=this.calibrationHistoryFromTo[1];i++){
    			int [] position=getPosition(i);
    			if (((position[0]-position[1])==(origin[0]-origin[1])) && ((position[2]-position[1])==(origin[2]-origin[1]))){
    				sameTiltmask[i]=true;
 //   				sameTiltSize++;
    			}
    		}
			// need to make sure the limits do not extend beyond actual samples (polynomials may go "crazy" there
    		this.corrXMin=getPositionAverage(this.calibrationHistoryFromTo[0]);
    		this.corrXMax=this.corrXMin;
    		for (int i=this.calibrationHistoryFromTo[0];i<=this.calibrationHistoryFromTo[1];i++) if (sameTiltmask[i]){
				double x=getPositionAverage(i);
				if (x<corrXMin) this.corrXMin=x;
				if (x>corrXMax) this.corrXMax=x;
			}
    		this.corrXMin=Math.floor(this.corrXMin)-1.0;
    		this.corrXMax=Math.ceil(this.corrXMax)+1.0;
    		this.calibrationNumPoints=focusMeasurementParameters.lensDistanceNumPoints;
    		if (debugLevel>0){
    			System.out.println( " Range is "+this.corrXMin+".."+this.corrXMax);
    		}
//    		xMin=corrXMin;
//    		xMax=corrXMax;

    		this.polyCoeff= calibrateFractions(
    				this.calibrationHistoryFromTo[0], // first history data to consider 
        			this.calibrationHistoryFromTo[1],   // last history data to consider
    				focusMeasurementParameters.lensDistancePolynomialDegree,  // polynomial degree to fit
					centerOnly, // only use center samples
					this.corrXMin,  // range of average motor coordinates to consider (should have the same tilt as last in history)
					this.corrXMax,
					debugLevel);
        	makeCalibrationTable(
        			this.polyCoeff,
        			this.calibrationNumPoints,
    				this.corrXMin,  // range of average motor coordinates to consider (should have the same tilt as last in history)
    				this.corrXMax,
        			debugLevel
        			);

    		String header="#\ttimestamp\tcenter\tResol-R\tResol-G\tResol-B\tfrac-R\tfrac-B\tY\tindex\tT-frac-R\tT-frac-B\tT-Y\tT-X";
//    		int indexR50=centerOnly?6:3; // 6 is average R50 in the center samples,3 - in all samples
    		StringBuffer sb = new StringBuffer();
    		double sumErrors=0.0;
    		double numSummedErrors=0;
    		
    		
    		
			for (int i=this.calibrationHistoryFromTo[0];i<=this.calibrationHistoryFromTo[1];i++)	if (sameTiltmask[i]){
				double x=getPositionAverage(i);
				String timestamp=getTimestamp(i);
				if (timestamp==null) timestamp="n/a";
				sb.append(i+"\t"+timestamp+"\t"+IJ.d2s(x,0));
				double [] sharp=this.history.get(i).getCenterResolutions();

				double l2=sharp[0]*sharp[0] + sharp[1]*sharp[1] + sharp[2]*sharp[2];
				sb.append("\t"+IJ.d2s(sharp[0],3)+"\t"+IJ.d2s(sharp[1],3)+"\t"+IJ.d2s(sharp[2],3)); // r,g,b here
				sb.append("\t"+IJ.d2s(sharp[0]*sharp[0]/l2,3));
				sb.append("\t"+IJ.d2s(sharp[2]*sharp[2]/l2,3));
				sb.append("\t"+IJ.d2s(l2,3));
				int j0=0;
				for (int j=0;j<this.lensDistanceCalibration.length;j++){
					if (this.lensDistanceCalibration[j][6]>x){
						j0=j;
						break;
					}
				}
				 sb.append("\t"+j0);
				if (j0==0) sb.append("\t---\t---\t---\t---\t---\t---\t---");
				else {
					double k= (x-this.lensDistanceCalibration[j0-1][6])/(this.lensDistanceCalibration[j0][6]-this.lensDistanceCalibration[j0-1][6]);
					for (int j=0;j<3;j++){
						double d=k*this.lensDistanceCalibration[j0][j]+(1.0-k)*this.lensDistanceCalibration[j0-1][j];
						sb.append("\t"+IJ.d2s(d,3));
					}
					double d=getLensDistance(
			    			sharp, // {R-sharpness,G-sharpness,B-sharpness}
			    			false, // boolean absolute, // return absolutely calibrated data
			    			focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
			    			focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
			    			debugLevel
			    			);
					sb.append("\t"+IJ.d2s(d,0));
					numSummedErrors++;
					sumErrors+=(d-x)*(d-x);

				}
				sb.append("\n");
			}
			// Add hysteresis if defined
			if (this.hysteresisHistoryFromTo!=null) {
				this.measuredHysteresis=0.0;
				for (int i=this.hysteresisHistoryFromTo[0];i<=this.hysteresisHistoryFromTo[1];i++)	{
					double x=getPositionAverage(i);
					String timestamp=getTimestamp(i);
					if (timestamp==null) timestamp="n/a";
					sb.append(i+"\t"+timestamp+"\t"+IJ.d2s(x,0));
					double [] sharp=this.history.get(i).getCenterResolutions();

					double l2=sharp[0]*sharp[0] + sharp[1]*sharp[1] + sharp[2]*sharp[2];
					sb.append("\t"+IJ.d2s(sharp[0],3)+"\t"+IJ.d2s(sharp[1],3)+"\t"+IJ.d2s(sharp[2],3)); // r,g,b here
					sb.append("\t"+IJ.d2s(sharp[0]*sharp[0]/l2,3));
					sb.append("\t"+IJ.d2s(sharp[2]*sharp[2]/l2,3));
					sb.append("\t"+IJ.d2s(l2,3));
					int j0=0;
					for (int j=0;j<this.lensDistanceCalibration.length;j++){
						if (this.lensDistanceCalibration[j][6]>x){
							j0=j;
							break;
						}
					}
					sb.append("\t"+j0);
					if (j0==0) sb.append("\t---\t---\t---\t---\t---\t---\t---");
					else {
						double k= (x-this.lensDistanceCalibration[j0-1][6])/(this.lensDistanceCalibration[j0][6]-this.lensDistanceCalibration[j0-1][6]);
						for (int j=0;j<3;j++){
							double d=k*this.lensDistanceCalibration[j0][j]+(1.0-k)*this.lensDistanceCalibration[j0-1][j];
							sb.append("\t"+IJ.d2s(d,3));
						}
						double d=getLensDistance(
								sharp, // {R-sharpness,G-sharpness,B-sharpness}
								false, // boolean absolute, // return absolutely calibrated data
								focusMeasurementParameters.lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
								focusMeasurementParameters.lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
								debugLevel
						);
						sb.append("\t"+IJ.d2s(d,0));
						if (this.measuredHysteresis<Math.abs((d-x))) this.measuredHysteresis=Math.abs((d-x));
					}
					sb.append("\n");
				}
			}
			double rms=Double.NaN;
			if (numSummedErrors>0) rms=Math.sqrt(sumErrors/numSummedErrors);
			String title= "Calibrating Lens Distance: RMS="+IJ.d2s(rms,3)+((this.hysteresisHistoryFromTo!=null)?(" hysteresis="+this.measuredHysteresis):"");

    		if (focusMeasurementParameters.showResults && focusMeasurementParameters.lensDistanceShowResults) new TextWindow(title, header, sb.toString(),1000,500);
    		//micronsPerStep
    		int [] optimalMotorPosition=absoluteCalibrate(
        			micronsPerStep,
        			focusMeasurementParameters.targetMicrons,
        			focusMeasurementParameters.lensDistanceWeightY, 
        			focusMeasurementParameters.lensDistanceWeightK,
        			debugLevel);
			double [] limits=absoluteCalibrationLimits() ; // here we know that calibration was made
			String footer=null;
    		if (optimalMotorPosition!=null) {
    			footer="Lens is calibrated from "+IJ.d2s(limits[0],3)+" microns to "+IJ.d2s(limits[1],3)+" microns from the \"best\" position\nBest offsets per color:";
// Add positions of color best focus    		
    			double [] bestColors=colorFocusDistances();
    			String [] colorNames={"red","green","blue"};
    			for (int c=0;c< bestColors.length;c++) footer+=" "+colorNames[c]+"="+(Double.isNaN(bestColors[c])?"out of range":(IJ.d2s(bestColors[c],2)+"um"));
    		} else {
    			footer="Scanned interval does not include best focus position.";
    			if (limits[0]==0.0) footer+=" Focus seems to be closer (smaller motor step values).";
    			if (limits[1]==0.0) footer+=" Focus seems to be farther (higher motor step values).";
    		}
    		if (this.hysteresisHistoryFromTo!=null){
    			footer+="Measured (maximal) motor hysteresis is "+IJ.d2s(this.measuredHysteresis,0)+" motor steps ("+IJ.d2s(micronsPerStep*this.measuredHysteresis,1)+" microns)";
    		}
    		footer+="\nComment: "+focusMeasurementParameters.comment;
			if (debugLevel>0) System.out.println(footer);
			if (focusMeasurementParameters.lensDistanceShowResults) IJ.showMessage("Success",footer);
    		// colorFocusDistances()
			if (focusMeasurementParameters.saveResults && (path!=null)){
    			CalibrationFileManagement.saveStringToFile (
						path,
						header+"\n"+sb.toString()+"\n"+footer);

			}
    		return optimalMotorPosition;
    	}

    	public double getLensDistance(
    			boolean absolute, // return absolutely calibrated data
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			int debugLevel
    			){
    		return getLensDistance(
    				getCenterResolutions(), // {R-sharpness,G-sharpness,B-sharpness}
        			absolute, // return absolutely calibrated data
        			weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
        			weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
        			debugLevel
        			);
    	}    	

    	public double getLensDistance(
    			int index,
    			boolean absolute, // return absolutely calibrated data
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			int debugLevel
    			){
    		return getLensDistance(
    				getCenterResolutions(index), // {R-sharpness,G-sharpness,B-sharpness}
        			absolute, // return absolutely calibrated data
        			weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
        			weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
        			debugLevel
        			);
    	}    	
    	
    	public double getLensDistance(
    			double [] sharpVector, // {R-sharpness,G-sharpness,B-sharpness}
    			boolean absolute, // return absolutely calibrated data
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			int debugLevel
    			){
 //   		System.out.println("getLensDistance(): debugLevel="+debugLevel);
 //   		System.out.println("getLensDistance(): weightK="+weightK);
 //   		System.out.println("getLensDistance(): weightY="+weightY);

    		int debugThreshold=2;
    		boolean useLogY=(weightY<0);
    		if (useLogY) weightY=-weightY;
    		// Caclulate yMmax
    		if (this.lensDistanceCalibration==null) return Double.NaN;
    		if (sharpVector==null) return Double.NaN;
    		double yMax=0.0;
    		for (int i=0;i<this.lensDistanceCalibration.length;i++) if (this.lensDistanceCalibration[i][2]>yMax)yMax=this.lensDistanceCalibration[i][2];
    		if (useLogY) yMax=Math.log(yMax);
    		double wY=weightY/yMax;
    		double best=Double.NaN;
    		int iBest=0;
    		double [] weights=new double[3];
    		double [] sample=new double[3];
			double s=0.0;
			for (int j=0;j<3;j++){
				sample[j]=sharpVector[j]*sharpVector[j];
				s+= sample[j];
			}
			sample[0]/=s;
			sample[1]=sample[2]/s;
			sample[2]=useLogY?Math.log(s):s;
if (debugLevel>=debugThreshold) System.out.println("sharpVector= {"+sharpVector[0]+" "+sharpVector[1]+" "+sharpVector[2]+"}, weightK="+weightK+" weightY="+weightY);
if (debugLevel>=debugThreshold) System.out.println("Sample="+sample[0]+" "+sample[1]+" "+sample[2]+" yMax="+yMax);
    		double [] diff=new double[3];
    		for (int i=0;i<this.lensDistanceCalibration.length;i++){
    			weights[0]=this.lensDistanceCalibration[i][3];
    			weights[1]=this.lensDistanceCalibration[i][4];
    			weights[2]=this.lensDistanceCalibration[i][5]*wY;
    			double w=weights[0]+weights[1]+weights[2];
    			for (int j=0;j<3;j++) weights[j]=weightK*weights[j]/w+(1.0-weightK)/3.0;
    			diff[0]=sample[0]-this.lensDistanceCalibration[i][0];
    			diff[1]=sample[1]-this.lensDistanceCalibration[i][1]; // sample[2], not 1 here (blue)
    			diff[2]=wY*(sample[2]-(useLogY?Math.log(this.lensDistanceCalibration[i][2]):this.lensDistanceCalibration[i][2]));

    			double err=0.0;
    			for (int j=0;j<3;j++) err+=weights[j]*diff[j]*diff[j];
    			if ((i==0) || (err<best)){
    				best=err;
    				iBest=i;
    			}
if (debugLevel>=debugThreshold) System.out.println(i+" "+diff[0]+" "+diff[1]+" "+diff[2]+" "+weights[0]+" "+weights[1]+" "+weights[2]+" "+err);
    		}
    		return this.lensDistanceCalibration[iBest][absolute?7:6];
    	}
    	/**
    	 * Applies absolute (from the "best") calibration to the distance between the lens and the center of the sensor
    	 * after the calibration ion motor steps for the range that includes "best focus" was performed.
    	 * "Best" here is the maximal sum of reversed areas of the focal spots (50% intensity) for each of the 3 colors,
    	 * measured for the "center" samples. After the calibration is performed, array lensDistanceCalibration can be used
    	 * to find the distance from the lens to the sensor center knowing the focal spot sizes for each of 3 colors.
    	 * 
    	 * @param micronsPerStep lens movement per step (all 3 motors together)
    	 * @param debugLevel - if >1, will produce debug output
    	 * @return motor position corresponding to zero shift ("best focus") or null if it is not inside
    	 * the scanned interval. May use absoluteCalibrationLimits() to find out which of the limits is wrong.
    	 */
    	public int[] absoluteCalibrate(
    			double micronsPerStep,
    			double targetMicrons,
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			int debugLevel){
// Temporary!
    		debugLevel=debugLevel+0; 
    		this.micronsPerStep=micronsPerStep;
    		if (this.lensDistanceCalibration==null) {
    			String msg="Lens travel is not yet measured/calibrated, impossible to apply absolute calibration";
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		
    		double bestY=0.0;
    		int i0=0;
    		for (int i=0;i<this.lensDistanceCalibration.length;i++){
    			if (this.lensDistanceCalibration[i][this.indexY]>bestY){
    				bestY=this.lensDistanceCalibration[i][this.indexY];
    				i0=i;
    			}
    		}
    		// use local maximum if the curve has N- shape
    		if ((i0==0) || (i0>=(this.lensDistanceCalibration.length-1))){
    			int i1=0;
    			int i2=this.lensDistanceCalibration.length-1;
    			while ((i1<(i2-1)) && (this.lensDistanceCalibration[i1+1][this.indexY]<this.lensDistanceCalibration[i1][this.indexY])) i1++; 
    			while ((i2>(i1+1)) && (this.lensDistanceCalibration[i2-1][this.indexY]<this.lensDistanceCalibration[i2][this.indexY])) i2--;
    			if (i1<(i2-1)){
    				bestY=0.0;
    				for (int i=i1;i<i2;i++){
    					if (this.lensDistanceCalibration[i][this.indexY]>bestY){
    						bestY=this.lensDistanceCalibration[i][this.indexY];
    						i0=i;
    					}
    				}
    			}
  				System.out.println("absoluteCalibrate(): ***** Warning: Absolute maximum of approximated sharpness function is not reached, using local maximum at "+i0+"/"+this.lensDistanceCalibration.length);
    		}
    		double x0=i0;
    		if ((i0>0) && (i0<(this.lensDistanceCalibration.length-1))){
    			double [][] data= {
    					{x0-1,this.lensDistanceCalibration[i0-1][this.indexY]},
    					{x0,  this.lensDistanceCalibration[i0  ][this.indexY]},
    					{x0+1,this.lensDistanceCalibration[i0+1][this.indexY]}
    			};
    			double [] coeff=(new PolynomialApproximation()).polynomialApproximation1d(data, 2);
    			if ((coeff!=null) && (coeff[2]<0.0)) {
    				x0=-coeff[1]/(2.0*coeff[2]);
    			}
    		}
			if (debugLevel>1) System.out.println ("absoluteCalibrate(): i0="+i0+", interpolated(?)  x0="+IJ.d2s(x0,3)+" motor average="+this.lensDistanceCalibration[i0][this.indexSteps]);
			double micronsPerElement=
				(this.lensDistanceCalibration[this.lensDistanceCalibration.length-1][indexSteps]-this.lensDistanceCalibration[0][indexSteps])/
				(this.lensDistanceCalibration.length-1)*this.micronsPerStep;
    		for (int i=0;i<this.lensDistanceCalibration.length;i++) {
    			this.lensDistanceCalibration[i][this.indexMicrons]= (i-x0)*micronsPerElement;
    			if (debugLevel>1) System.out.println ("absoluteCalibrate(): this.lensDistanceCalibration["+i+"]["+this.indexMicrons+"]="+this.lensDistanceCalibration[i][this.indexMicrons]);
    		}
//    		int [] position= parallelPositionMove(targetMicrons, weightK, weightY, debugLevel); // "optimal" motor position (or null)
    		int [] position= parallelPositionMove(0.0, weightK, weightY, debugLevel); // "optimal" motor position (or null) // here - disregard targetMicrons
    		if (debugLevel>0) if (position!=null) System.out.println("Optimal focus in the center at position ["+position[0]+","+position[1]+","+position[2]+"]");
    		else System.out.println("Optimal focus is undefined");
    		return position;
    	}

    	/**
    	 * 
    	 * @param targetMicrons requested distance from the lens to the center, 0 is the "best focus" (between green and red), in microns
    	 * @param weightK scale between derivative-driven (1.0) and equal (0.0) weights in fitting measured parameters to find distance
    	 * @param weightY relative weight of overall (not color-differential) "sharpness" in fitting
    	 * @param debugLevel
    	 * @return new motors position
    	 */
    	
    	public int [] parallelPositionMove(
    			double targetMicrons,
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			int debugLevel
    			){
    		if (this.lensDistanceCalibration==null){
    			if (debugLevel>0) System.out.println("parallelPositionMove(): this.lensDistanceCalibration==null"); 
    			return null;// not yet calibrated
    		}
    		int last=this.lensDistanceCalibration.length-1;

    		double [] limits=absoluteCalibrationLimits();
    		if ((limits==null) || ((limits[0]*limits[1])>0) || ((limits[0]==0) && (limits[1]==0))){
    			if ((debugLevel>0) && (limits==null)) System.out.println("parallelPositionMove(): limits==null"); 
    			if ((debugLevel>0) && ((limits[0]*limits[1])>=0)){
    				System.out.println("parallelPositionMove(): ((limits[0]*limits[1])>=0), limits[0]="+limits[0]+", limits[1]="+limits[1]); 
    			}
    			return null;
    		}

    		if (
    				(targetMicrons<this.lensDistanceCalibration[0][indexMicrons]) ||
    				(targetMicrons>this.lensDistanceCalibration[last][indexMicrons])) {
				System.out.println("parallelPositionMove(): out of range, this.lensDistanceCalibration[0]["+indexMicrons+"]="+
						this.lensDistanceCalibration[0][indexMicrons]+
						", this.lensDistanceCalibration["+last+"]["+indexMicrons+"]="+this.lensDistanceCalibration[last][indexMicrons]); 
    			return null ; // out of range
    		}
    		double x=
    			(targetMicrons - this.lensDistanceCalibration[0][indexMicrons])/
    			(this.lensDistanceCalibration[last][indexMicrons]-this.lensDistanceCalibration[0][indexMicrons]);
    		if (debugLevel>1) System.out.println("parallelPositionMove(): x="+x); 
    		double newAveragePos=
    			this.lensDistanceCalibration[0][this.indexSteps]+
    			x*(this.lensDistanceCalibration[last][indexSteps]-this.lensDistanceCalibration[0][indexSteps]);
    		if (debugLevel>1) System.out.println("parallelPositionMove(): newAveragePos="+newAveragePos); 
			int [] position=getPosition(); // last position from the history ** here was the bug (made returning clone()
//			double [] resolutions=getCenterResolutions();
			
			double predictedAverage=getLensDistance(
					getCenterResolutions(), // {R-sharpness,G-sharpness,B-sharpness}
	    			false, //true, // return absolutely calibrated data
	    			weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
	    			weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
	    			debugLevel);
			if (Double.isNaN(predictedAverage)) return null;
			int move = (int) (newAveragePos-predictedAverage);
    		if (debugLevel>1) System.out.println("parallelPositionMove(): calculated move of the center is "+move+" steps of 3 motors in parallel."); 
			for (int i=0;i<3;i++) position[i]+=move;
    		return position; // "optimal" motor position
    		
    	}
    	
    	/**
    	 * Calculates offsets (in microns) of each color best focus (if it is inside range) 
    	 * @return
    	 */
    	public double [] colorFocusDistances(){
    		int [] bestIndices={0,0,0};
    		double [] bestResolution={0.0,0.0,0.0};
    		for (int i=0;i<this.lensDistanceCalibration.length;i++) {
    			double l2=    this.lensDistanceCalibration[i][this.indexY];
    			double rFrac2=this.lensDistanceCalibration[i][this.indexR];
    			double bFrac2=this.lensDistanceCalibration[i][this.indexB];
    			double [] resolutions={
    					Math.sqrt(rFrac2*l2),
    					Math.sqrt((1-rFrac2-bFrac2)*l2),
    					Math.sqrt(bFrac2*l2)};
    			for (int c=0;c< resolutions.length;c++) if (resolutions[c]>bestResolution[c]){
    				bestResolution[c]=resolutions[c];
    				bestIndices[c]=i;
    			}
    		}
    		double [] bestFocus ={Double.NaN,Double.NaN,Double.NaN}; // maximum outside of range
    		for (int c=0;c<bestResolution.length;c++) if ((bestIndices[c]>0) && (bestIndices[c]<(this.lensDistanceCalibration.length-1))) {
    			double [][] data= {
    					{bestIndices[c]-1,0.0},
    					{bestIndices[c],  0.0},
    					{bestIndices[c]+1,0.0}
    			};
    			for (int j=-1;j<2;j++){
        			double l2=    this.lensDistanceCalibration[bestIndices[c]+j][this.indexY];
        			double rFrac2=this.lensDistanceCalibration[bestIndices[c]+j][this.indexR];
        			double bFrac2=this.lensDistanceCalibration[bestIndices[c]+j][this.indexB];
        			double [] resolutions={
        					Math.sqrt(rFrac2*l2),
        					Math.sqrt((1-rFrac2-bFrac2)*l2),
        					Math.sqrt(bFrac2*l2)};
        			data[j+1][1]=resolutions[c];
    			}
    			double [] coeff=(new PolynomialApproximation()).polynomialApproximation1d(data, 2);
    			if ((coeff!=null) && (coeff[2]<0.0)) {
    				double x=(-coeff[1]/(2.0*coeff[2]))-bestIndices[c];
    				
    				bestFocus[c]= this.lensDistanceCalibration[bestIndices[c]][this.indexMicrons]+
    				x*(this.lensDistanceCalibration[bestIndices[c]+1][this.indexMicrons]-this.lensDistanceCalibration[bestIndices[c]][this.indexMicrons]);
    			}
    		}
    		return bestFocus;
    	}
    	
     	/**
    	 *  
    	 * @return array of lower and upper limits of the calibration range. None should be 0.0, that would meant that maximum
    	 * was not covered by the scanned interval Returns null if there is no calibration 
    	 */
    	public double [] absoluteCalibrationLimits(){
    		if (this.lensDistanceCalibration==null)  return null; // no calibration at all
    		if (Double.isNaN(this.lensDistanceCalibration[0][this.indexMicrons])) return null; // no absolute calibration
    		double [] limits={
    				this.lensDistanceCalibration[0][this.indexMicrons],
    				this.lensDistanceCalibration[this.lensDistanceCalibration.length-1][this.indexMicrons]
    		};
    		return limits;

    	}    	
/*    	
    	public void calibrateFractions( // positive away from the lens!
    			int fromHistory, // first history data to consider
    			int toHistory,   // last history data to consider
    			int numPoints,
    			int polyDegree,  // polynomial degree to fit
//    			boolean sameTiltOnly,
    			boolean centerOnly, // only use center samples
    			double xMin,  // range of average motor coordinates to consider (should have the same tilt as last in history)
    			double xMax,
    			int debugLevel
    			){
    		if (fromHistory<0)fromHistory=0;          // already done this by the caller
    		if (toHistory>=size())toHistory=size()-1;

    		int indexR50=centerOnly?6:3; // 6 is average R50 in the center samples,3 - in all samples
    		int sameTiltSize=0;
//    		int [] origin=getPosition();
    		int [] origin=getPosition(toHistory);
    		boolean []sameTiltmask=new boolean[size()];
    		for (int i=0;i<sameTiltmask.length;i++) sameTiltmask[i]=false;
    		for (int i=fromHistory;i<=toHistory;i++){
    			sameTiltmask[i]=true;
    			sameTiltSize++;
    			int [] position=getPosition(i);
    			if (((position[0]-position[1])==(origin[0]-origin[1])) && ((position[2]-position[1])==(origin[2]-origin[1]))){
    				sameTiltmask[i]=true;
    				sameTiltSize++;
    			}
    		}

//    		System.out.println("solveSinglePoly(): sameTiltSize="+sameTiltSize+", size()="+size()+" debugLevel="+debugLevel);
    		if (debugLevel>1) System.out.println("solveSinglePoly(): sameTiltSize="+sameTiltSize+", size()="+size()+
    				" fromHistory="+fromHistory+" toHistory="+toHistory+" numPoints="+numPoints);
//    		double x0=(origin[0]+origin[1]+origin[2])/3.0;
			double x0=(xMax+xMin)/2;// center
			double scale=2.0/(xMax-xMin);
    		double [][] polyCoeff =new double[3][];
			PolynomialApproximation pa= new PolynomialApproximation(debugLevel);
    		for (int color=0;color<3;color++){
    			double [][] data=new double [size()][3]; //x,f,w
    			for (int i=0;i<size();i++) {
    				data[i][2]=0.0; // weight
    				if (sameTiltmask[i]){
    					double x=getPositionAverage(i);
    					if ((x<xMin) || (x>xMax)) continue; // out of range point
    					data[i][0]= scale*(x - x0); // interval center
    					double [][] metrics= this.history.get(i).getMetrics(0.0,0.0,0.0); // average is not used, any scales
    					double l2=
    						1.0/metrics[0][indexR50]/metrics[0][indexR50]+
    						1.0/metrics[1][indexR50]/metrics[1][indexR50]+
    						1.0/metrics[2][indexR50]/metrics[2][indexR50];
    					if      (color==0) data[i][1]=1.0/metrics[0][indexR50]/metrics[0][indexR50]/l2; // red 
    					else if (color==1) data[i][1]=1.0/metrics[2][indexR50]/metrics[2][indexR50]/l2; // blue
    					else               data[i][1]=l2; // Y-combined
 //   					data[i][1]=1.0/metrics[color][indexR50];// "sharpness"
    					data[i][2]=1.0;
    					if (debugLevel>1)	System.out.println("color=\t"+color+"\ti=\t"+i+"\tdata[i][0]=\t"+data[i][0]+"\tdata[i][1]=\t"+data[i][1]+"\tdata[i][2]=\t"+data[i][2]);
    				}
    			}
    			polyCoeff[color]=pa.polynomialApproximation1d(data, polyDegree);
    			if (debugLevel>1) {
    				String dbgStr="polyCoeff["+color+"]=";
    				for (int i=0;i<polyCoeff[color].length;i++) {
    					dbgStr+="\t"+polyCoeff[color][i];
    				}
    				System.out.println(dbgStr);
    			}
    		}
    		this.lensDistanceCalibration=new double [numPoints][8]; // fracR, fracG, fracB, weightR, weightG, weightB, centerSteps, microns (not assigned here)
			for (int n=0;n<numPoints;n++) {
				this.lensDistanceCalibration[n][6]=xMin+n*(xMax-xMin)/(numPoints-1); 
				double x=(2.0*n)/(numPoints-1)-1.0;
				double [] xp=new double[polyDegree+1];
				xp[0]=1;
				for (int i=1;i<=polyDegree;i++) xp[i]=xp[i-1]*x;
//				double sumDerSq=0.0;
				for (int color=0;color<3;color++){
					this.lensDistanceCalibration[n][color]=0;
					this.lensDistanceCalibration[n][color+3]=0; // derivatives^2
					for (int i=0;i<=polyDegree;i++){
						this.lensDistanceCalibration[n][color]+=polyCoeff[color][i]*xp[i];
					}
					for (int i=1;i<=polyDegree;i++){
						this.lensDistanceCalibration[n][color+3]+=polyCoeff[color][i]*xp[i-1]*i; // derivative
					}
					this.lensDistanceCalibration[n][color+3]*=this.lensDistanceCalibration[n][color+3]; //squared
//					sumDerSq+=this.lensDistanceCalibration[n][color+3];
				}
//				for (int color=0;color<3;color++)this.lensDistanceCalibration[n][color+3] /= sumDerSq;// should be non-zero, I hope 
				this.lensDistanceCalibration[n][7]=Double.NaN; // not yet absolutely calibrated
				if (debugLevel>1) System.out.println(n+":\t"+
						IJ.d2s(this.lensDistanceCalibration[n][0],4)+"\t"+
						IJ.d2s(this.lensDistanceCalibration[n][1],4)+"\t"+
						IJ.d2s(this.lensDistanceCalibration[n][2],4)+"\t"+
						IJ.d2s(this.lensDistanceCalibration[n][3],4)+"\t"+
						IJ.d2s(this.lensDistanceCalibration[n][4],4)+"\t"+
						IJ.d2s(this.lensDistanceCalibration[n][5],4)+"\t"+
						IJ.d2s(this.lensDistanceCalibration[n][6],4));
			}
			
			
			//    	private double lensDistanceCalibrationYScale=0.0; // R-frac, B-frac, Y*scale

//    		return calib; // should be stored between calls
    	}
 */   	
 
       	
    	public double [][] calibrateFractions( // positive away from the lens!
    			int fromHistory, // first history data to consider
    			int toHistory,   // last history data to consider
//    			int numPoints,
    			int polyDegree,  // polynomial degree to fit
    			boolean centerOnly, // only use center samples
    			double xMin,  // range of average motor coordinates to consider (should have the same tilt as last in history)
    			double xMax,
    			int debugLevel
    			){
    		if (fromHistory<0)fromHistory=0;          // already done this by the caller
    		if (toHistory>=size())toHistory=size()-1;

    		int indexR50=centerOnly?6:3; // 6 is average R50 in the center samples,3 - in all samples
    		int sameTiltSize=0;
    		int [] origin=getPosition(toHistory);
    		boolean []sameTiltmask=new boolean[size()];
    		for (int i=0;i<sameTiltmask.length;i++) sameTiltmask[i]=false;
    		for (int i=fromHistory;i<=toHistory;i++){
    			sameTiltmask[i]=true;
    			sameTiltSize++;
    			int [] position=getPosition(i);
    			if (((position[0]-position[1])==(origin[0]-origin[1])) && ((position[2]-position[1])==(origin[2]-origin[1]))){
    				sameTiltmask[i]=true;
    				sameTiltSize++;
    			}
    		}
    		if (debugLevel>1) System.out.println("solveSinglePoly(): sameTiltSize="+sameTiltSize+", size()="+size()+
    				" fromHistory="+fromHistory+" toHistory="+toHistory);
			double x0=(xMax+xMin)/2;// center
			double scale=2.0/(xMax-xMin);
    		double [][] polyCoeff =new double[3][];
			PolynomialApproximation pa= new PolynomialApproximation(debugLevel);
    		for (int color=0;color<3;color++){
    			double [][] data=new double [size()][3]; //x,f,w
    			for (int i=0;i<size();i++) {
    				data[i][2]=0.0; // weight
    				if (sameTiltmask[i]){
    					double x=getPositionAverage(i);
    					if ((x<xMin) || (x>xMax)) continue; // out of range point
    					data[i][0]= scale*(x - x0); // interval center
    					double [][] metrics= this.history.get(i).getMetrics(0.0,0.0,0.0); // average is not used, any scales
    					if (metrics==null) continue;
    					double l2=
    						1.0/metrics[0][indexR50]/metrics[0][indexR50]+
    						1.0/metrics[1][indexR50]/metrics[1][indexR50]+
    						1.0/metrics[2][indexR50]/metrics[2][indexR50];
    					if      (color==0) data[i][1]=1.0/metrics[0][indexR50]/metrics[0][indexR50]/l2; // red 
    					else if (color==1) data[i][1]=1.0/metrics[2][indexR50]/metrics[2][indexR50]/l2; // blue
    					else               data[i][1]=l2; // Y-combined
 //   					data[i][1]=1.0/metrics[color][indexR50];// "sharpness"
    					data[i][2]=1.0;
    					if (debugLevel>1)	System.out.println("color=\t"+color+"\ti=\t"+i+"\tdata[i][0]=\t"+data[i][0]+"\tdata[i][1]=\t"+data[i][1]+"\tdata[i][2]=\t"+data[i][2]);
    				}
    			}
    			polyCoeff[color]=pa.polynomialApproximation1d(data, polyDegree);
    			if (debugLevel>1) {
    				String dbgStr="polyCoeff["+color+"]=";
    				for (int i=0;i<polyCoeff[color].length;i++) {
    					dbgStr+="\t"+polyCoeff[color][i];
    				}
    				System.out.println(dbgStr);
    			}
    		}
    		return polyCoeff;
			
			
			//    	private double lensDistanceCalibrationYScale=0.0; // R-frac, B-frac, Y*scale

//    		return calib; // should be stored between calls
    	}
    	public void makeCalibrationTable(
    			double [][] polyCoeff,
    			int numPoints,
//    			int polyDegree,  // polynomial degree to fit
    			double xMin,  // range of average motor coordinates to consider (should have the same tilt as last in history)
    			double xMax,
    			int debugLevel
    			){
//    		debugLevel=2;    			
   			int polyDegree=polyCoeff[0].length-1;  // polynomial degree
    		this.lensDistanceCalibration=new double [numPoints][8]; // fracR, fracG, fracB, weightR, weightG, weightB, centerSteps, microns (not assigned here)
    		if (debugLevel>1) {
    			for (int c=0;c<3;c++) {
    				System.out.println("makeCalibrationTable(), c="+c);
    				for (int n=0;n<polyDegree;n++) System.out.println(n+": "+polyCoeff[c][n]);
    			}
    		}
    		for (int n=0;n<numPoints;n++) {
    			this.lensDistanceCalibration[n][6]=xMin+n*(xMax-xMin)/(numPoints-1); 
    			double x=(2.0*n)/(numPoints-1)-1.0;
    			double [] xp=new double[polyDegree+1];
    			xp[0]=1;
    			for (int i=1;i<=polyDegree;i++) xp[i]=xp[i-1]*x;
    			//		double sumDerSq=0.0;
    			for (int color=0;color<3;color++){
    				this.lensDistanceCalibration[n][color]=0;
    				this.lensDistanceCalibration[n][color+3]=0; // derivatives^2
    				for (int i=0;i<=polyDegree;i++){
    					this.lensDistanceCalibration[n][color]+=polyCoeff[color][i]*xp[i];
    				}
    				for (int i=1;i<=polyDegree;i++){
    					this.lensDistanceCalibration[n][color+3]+=polyCoeff[color][i]*xp[i-1]*i; // derivative
    				}
    				this.lensDistanceCalibration[n][color+3]*=this.lensDistanceCalibration[n][color+3]; //squared
    				//			sumDerSq+=this.lensDistanceCalibration[n][color+3];
    			}
    			//		for (int color=0;color<3;color++)this.lensDistanceCalibration[n][color+3] /= sumDerSq;// should be non-zero, I hope 
    			this.lensDistanceCalibration[n][7]=Double.NaN; // not yet absolutely calibrated
    			if (debugLevel>1) System.out.println(n+":\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][0],4)+"\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][1],4)+"\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][2],4)+"\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][3],4)+"\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][4],4)+"\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][5],4)+"\t"+
    					IJ.d2s(this.lensDistanceCalibration[n][6],4));
    		}

    	}

    	
    	/**
    	 *  Calculates worst ratio of the residual error to threshold. First parameters is focus or radial/tangential ratio, second and third - tilts
    	 *  Tilts should be 0, while the first is compared against "target"
    	 * @param weightRatioRedToGreen
    	 * @param weightRatioBlueToGreen
    	 * @param lensDistanceWeightK
    	 * @param lensDistanceWeightY
    	 * @param useRadialTangential
    	 * @param targetFarNear
    	 * @param targetMicrons
    	 * @param toleranceMicrons
    	 * @param toleranceTilt
    	 * @param toleranceThreshold
    	 * @param debugLevel
    	 * @return
    	 */
     	public double worstOverTolerance (
    			double weightRatioRedToGreen,
    			double weightRatioBlueToGreen,
    			double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
				double targetMicrons,
				double toleranceMicrons,
				double toleranceTilt,
    			int debugLevel
    	){ return worstOverTolerance (
    			weightRatioRedToGreen,
    			weightRatioBlueToGreen,
    			lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
				false, //boolean useRadialTangential,
    			0.0,   //double targetFarNear,
				targetMicrons,
				toleranceMicrons,
				toleranceTilt,
    			debugLevel);
     	}
    	public double worstOverTolerance (
    			double weightRatioRedToGreen,
    			double weightRatioBlueToGreen,
    			double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
				boolean useRadialTangential,
    			double targetFarNear,
				double targetMicrons,
				double toleranceMicrons,
				double toleranceTilt,
    			int debugLevel
    	){
    		if (this.history.size()<1)  return Double.NaN;
    		if (!useRadialTangential){
    			// Make sure calibration is done and targetMicrons is within limits
    			double [] limits=absoluteCalibrationLimits();
    			if ((limits==null) || Double.isNaN(limits[0]) || Double.isNaN(limits[1])) return Double.NaN; // nonexistent/invalid calibration
    		}
    		double targetValue=useRadialTangential?targetFarNear:targetMicrons;
    		double [] lastAverageMetrics=this.history.get(this.history.size()-1).getMetrics(weightRatioRedToGreen,weightRatioBlueToGreen)[3];
    		double [] lastResolutions= this.history.get(this.history.size()-1).getCenterResolutions();
    		if (!useRadialTangential){
    			lastAverageMetrics[0]=getLensDistance(
    					lastResolutions, // {R-sharpness,G-sharpness,B-sharpness}
		    			true, // boolean absolute, // return absolutely calibrated data
		    			lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
		    			lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
		    			1 //debugLevel
		    			);
    		}
    		lastAverageMetrics[0]-=targetValue;
    		lastAverageMetrics[0]/=useRadialTangential?toleranceTilt:toleranceMicrons; //
    		lastAverageMetrics[1]/=toleranceTilt;
    		lastAverageMetrics[2]/=toleranceTilt;
    		return Math.max(Math.abs(lastAverageMetrics[0]),Math.max(Math.abs(lastAverageMetrics[1]),Math.abs(lastAverageMetrics[2])));
     	}
     	
     	public int [] focusReadjustStep(
     			double targetMicrons, // target focal distance
     			double micronsFade, // reduce influence of far points
    			double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			int debugLevel
     			){
			double [] limits=absoluteCalibrationLimits();
    		int last=this.history.size()-1;
			if (limits==null){
    			String msg="Failed: lens focus quality per colors vs. distance is not yet measured/calibrated";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
			}
			if (Double.isNaN(limits[0]) || Double.isNaN(limits[1])){
    			String msg="Failed: lens focus quality per colors vs. distance is measured, but absolute calibration was not performed.";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
			}
			if ((targetMicrons<limits[0]) || (targetMicrons>limits[1])){
    			String msg="Failed: requested distance ("+targetMicrons+"um) is not withing calibrated range "+limits[0]+".."+limits[1]+"um.";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
			}
			if (this.history.size()<1){
    			String msg="Failed: history is empty.";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
				
			}
    		if (
    				(targetMicrons<this.lensDistanceCalibration[0][indexMicrons]) ||
    				(targetMicrons>this.lensDistanceCalibration[this.lensDistanceCalibration.length-1][indexMicrons])) {
				System.out.println("parallelPositionMove(): out of range, this.lensDistanceCalibration[0]["+indexMicrons+"]="+
						this.lensDistanceCalibration[0][indexMicrons]+
						", this.lensDistanceCalibration["+(this.lensDistanceCalibration.length-1)+"]["+indexMicrons+"]="+this.lensDistanceCalibration[this.lensDistanceCalibration.length-1][indexMicrons]); 
    			return null ; // out of range
    		}
			if (debugLevel>1) {
				System.out.println("focusReadjustStep(), last="+last);
			}

    		int [] lastPosition=getPosition();
    		// use only same difference between motors
    		int [] diffs= {lastPosition[0]-lastPosition[2],lastPosition[1]-lastPosition[2]};
			double [][] data = new double[this.history.size()][3];
    		double k=(micronsFade>0.0)?0.5/micronsFade/micronsFade:1.0;
    		int numSameDiff=0;
    		double [] resolutions;
			for (int nSample=0;nSample<this.history.size();nSample++){
				int [] motors=       getPosition(nSample);
				if ((motors[0]-motors[2]==diffs[0]) && (motors[1]-motors[2]==diffs[1])){
					resolutions= this.history.get(nSample).getCenterResolutions();
					data [nSample][0]=(motors[0]+motors[1]+motors[2])/3.0;
					data [nSample][1]=getLensDistance(
							resolutions, // {R-sharpness,G-sharpness,B-sharpness}
							true, // boolean absolute, // return absolutely calibrated data
							lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
							lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
							debugLevel
					);
					if (micronsFade<=0) data [nSample][2]=1.0;
					else {
						double diffDistance=data [nSample][1]-targetMicrons;
						data [nSample][2]=Math.exp(-k*diffDistance*diffDistance);
					}
					if (nSample==(last)) data [nSample][2]=1.0; // last measurement full weight?
					numSameDiff++;
				} else {
					for (int j=0; j<data[nSample].length;j++) data [nSample][j]=0.0;
					
				}
			}
			double [] coeff=(new PolynomialApproximation( debugLevel)).polynomialApproximation1d(data, 1); // just linear approximation
			int move=0;
			if ((coeff.length==1)|| (coeff[1]==0.0)){ // single sample or several for the same motor position
				move=(int) Math.round((targetMicrons-data[last][1])/this.micronsPerStep);
			} else {
				move=(int) Math.round((targetMicrons-coeff[0])/coeff[1] - data [last][0]); 
				// linear approximate
			}
			if (debugLevel>1) {
				System.out.println("focusReadjustStep(), numSameDiff="+numSameDiff+" targetMicrons="+targetMicrons+" data [last][0]="+data [last][0]);
				for (int ii=0;ii<coeff.length;ii++){
					System.out.println("focusReadjustStep(), coeff["+ii+"]="+coeff[ii]);
				}
				System.out.println("focusReadjustStep(), move="+move);
			}
			int [] result={lastPosition[0]+move,lastPosition[1]+move,lastPosition[2]+move};
			return result;
     	}

     	public int [][] probeSequence(
     			int [] motors, // should not be null
     			boolean symmetrical,
     			double probe_M1M2M3, // sigma for decay for average of the 3 motors: (M1+M2+M3)/3
     			double probe_M3_M1M2, // sigma for decay for M3 opposite to M1 and M2: M3-(M1+M2)/2
     			double probe_M2_M1){ // sigma for decay for  M2 opposite to M1:    M2-M1
     		double [] orthoCenter=motorsToOrtho(motors);
    		double [][] seqShort={
    				{probe_M1M2M3,          0,                0        },
    				{     0,           probe_M3_M1M2,         0        },
    				{     0,                0,           probe_M2_M1   },
    				{-probe_M1M2M3/3, -probe_M3_M1M2/3, -probe_M2_M1/3 },
    				{     0,                0,                0}};
    		double [][] seqSymm={
    				{ probe_M1M2M3,      0,             0},
    				{-probe_M1M2M3,      0,             0},
    				{      0,       probe_M3_M1M2,      0},
    				{      0,      -probe_M3_M1M2,      0},
    				{      0,            0,        probe_M2_M1},
    				{      0,            0,       -probe_M2_M1},
    				{      0,            0,             0}};

     		double [][] seq=symmetrical?seqSymm:seqShort;
     		int [][] motorSeq=new int [seq.length][];
     		for (int i =0; i<motorSeq.length;i++) {
     			double [] ortho= {
     					orthoCenter[0]+seq[i][0],
     					orthoCenter[1]+seq[i][1],
     					orthoCenter[2]+seq[i][2]};
     			
     			motorSeq[i]= motorsFromOrtho(ortho);
     		}
     		return motorSeq;
     	}
/**
 * Approximate focal distance versus sensor temperature as a linera function 
 * @param numSamples number of last samples from history to use
 * @param lensDistanceWeightK coefficients used for focal distance calculation
 * @param lensDistanceWeightY
 * @return array of 2 elements - {length at 0C, microns/degree}
 */
     	public double [] temperatureLinearApproximation(
     			boolean useLMA,
     			int numSamples,             // number of last samples from history to use, 0 - use all
    			double lensDistanceWeightK, // used in distance calculation 
    			double lensDistanceWeightY // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
     			){

     		if (numSamples<=0) numSamples=this.history.size();
     		if (numSamples>this.history.size()) numSamples=this.history.size();
     		double [][] data =new double [numSamples][2]; // no weights
    		double [] resolutions;
    		int firstSample=this.history.size()-numSamples;

     		for (int nSample=0;nSample<numSamples;nSample++){
				data[nSample][0]=this.history.get(firstSample+nSample).getTemperature();
				if (useLMA){
					data[nSample][1]=this.history.get(firstSample+nSample).getzTxTy()[0];
				} else {
					resolutions= this.history.get(firstSample+nSample).getCenterResolutions();
					data[nSample][1]=getLensDistance(
							resolutions, // {R-sharpness,G-sharpness,B-sharpness}
							true, // boolean absolute, // return absolutely calibrated data
							lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
							lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
							1 //debugLevel
							);
				}
     		}
			PolynomialApproximation pa= new PolynomialApproximation(debugLevel);
			double [] polyCoeff=pa.polynomialApproximation1d(data, 1); // just linear
			return polyCoeff;
     	}

     	public double [] temperatureLinearApproximation(
     			double [][] ZTT,            // Z, tXx, tY, m1,m2,m3 found with LMA - may contain NaN
     			int numSamples             // number of last samples from history to use, 0 - use all
     			){
//     		if (numSamples<=0)  this.history.size();
//     		if (numSamples>this.history.size()) numSamples=this.history.size();
//    		int firstSample=this.history.size()-numSamples;
    		
     		if (numSamples<=0) numSamples=ZTT.length;
     		if (numSamples>ZTT.length) numSamples=ZTT.length;
    		int firstSample=ZTT.length-numSamples;
    		int numGoodSamples=0;
     		for (int nSample=0;nSample<numSamples;nSample++){
     			if (ZTT[firstSample+nSample]!=null) numGoodSamples++;
     		}
     		double [][] data =new double [numGoodSamples][2]; // no weights
     		int index=0;
     		for (int nSample=0;nSample<numSamples;nSample++) if (ZTT[firstSample+nSample]!=null) {
//				data[index][0]=this.history.get(firstSample+nSample).getTemperature();
				data[index][0]=ZTT[firstSample+nSample][3];
				data[index++][1]=ZTT[firstSample+nSample][0]; // Z
     		}
			PolynomialApproximation pa= new PolynomialApproximation(debugLevel);
			double [] polyCoeff=pa.polynomialApproximation1d(data, 1); // just linear
			return polyCoeff;
     	}
     	
     	
/**
 *  will return same coordinates if tolerances are met     	
 * @param useTheBest
 * @param targetMicrons
 * @param micronsFade
 * @param tiltFade
 * @param sigma_M1M2M3
 * @param sigma_M3_M1M2
 * @param sigma_M2_M1
 * @param lensDistanceWeightK
 * @param lensDistanceWeightY
 * @param weightRatioRedToGreen
 * @param weightRatioBlueToGreen
 * @param toleranceMicrons
 * @param toleranceTilt
 * @param toleranceThreshold
 * @param maxStep
 * @param debugLevel
 * @return
 */
     	public int [] focusTiltStep(
     			boolean useTheBest, // start from the best sample, (false - from the last)
     			double targetMicrons, // target focal distance
     			double micronsFade, // reduce influence of far points
     			double tiltFade,    //goodTiltSigma
     			double sigma_M1M2M3, // sigma for decay for average of the 3 motors: (M1+M2+M3)/3
     			double sigma_M3_M1M2, // sigma for decay for M3 opposite to M1 and M2: M3-(M1+M2)/2
     			double sigma_M2_M1, // sigma for decay for  M2 opposite to M1:    M2-M1
    			double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			double weightRatioRedToGreen,  // used for tilt averaging
    			double weightRatioBlueToGreen, // used for tilt averaging
    			
    			double toleranceMicrons,
    			double toleranceTilt,
    			double toleranceThreshold,
    			double maxStep,        // Maximal motor move for each step

    			int debugLevel
     			){
     		int debugThreshold=0;
			double [] limits=absoluteCalibrationLimits();
    		int last=this.history.size()-1;
			if (limits==null){
    			String msg="Failed: lens focus quality per colors vs. distance is not yet measured/calibrated";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
			}
			if (Double.isNaN(limits[0]) || Double.isNaN(limits[1])){
    			String msg="Failed: lens focus quality per colors vs. distance is measured, but absolute calibration was not performed.";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
			}
			if ((targetMicrons<limits[0]) || (targetMicrons>limits[1])){
    			String msg="Failed: requested distance ("+targetMicrons+"um) is not withing calibrated range "+limits[0]+".."+limits[1]+"um.";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
			}
			if (this.history.size()<1){
    			String msg="Failed: history is empty.";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
				
			}
    		if (
    				(targetMicrons<this.lensDistanceCalibration[0][indexMicrons]) ||
    				(targetMicrons>this.lensDistanceCalibration[this.lensDistanceCalibration.length-1][indexMicrons])) {
				System.out.println("parallelPositionMove(): out of range, this.lensDistanceCalibration[0]["+indexMicrons+"]="+
						this.lensDistanceCalibration[0][indexMicrons]+
						", this.lensDistanceCalibration["+(this.lensDistanceCalibration.length-1)+"]["+indexMicrons+"]="+this.lensDistanceCalibration[this.lensDistanceCalibration.length-1][indexMicrons]); 
    			return null ; // out of range
    		}
			if (debugLevel>1) {
				System.out.println("focusTiltStep(), last="+last);
			}
// see if tolerances are met - return with unchanged position. If not met, but close - scale down step
    		int [] lastPosition=getPosition();
    		double  worstOver= worstOverTolerance (
        			weightRatioRedToGreen,
        			weightRatioBlueToGreen,
        			lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
        			lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    				targetMicrons,
    				toleranceMicrons,
    				toleranceTilt,
        			debugLevel);
	        	if (debugLevel>debugThreshold) {
	        		System.out.println("focusTiltStep() worstOver="+IJ.d2s(worstOver,5)+", m1="+lastPosition[0]+"m2="+lastPosition[1]+"m3="+lastPosition[2]);
	        	}

    		if (worstOver<1.0) return lastPosition; // no correction if worst is less than tolerance
    		
    		
    		boolean scaleDownCorrection=(worstOver<toleranceThreshold); 
    		double scaleDownValue=scaleDownCorrection?(worstOver/toleranceThreshold):1.0;

			
			double [][][] data = new double[this.history.size()][3][]; //{u,v,t},{f1(u,v,t), f2(u,v,t), f3(u,v,t)),{w}
    		double kF=(micronsFade>0.0)?0.5/micronsFade/micronsFade:1.0;
    		double kT=(tiltFade>0.0)?0.5/tiltFade/tiltFade:1.0;
    		double [] kP={
    				((sigma_M1M2M3>0.0)?(0.5/sigma_M1M2M3/sigma_M1M2M3):1.0),
    				((sigma_M3_M1M2>0.0)?(0.5/sigma_M3_M1M2/sigma_M3_M1M2):1.0),
    				((sigma_M2_M1>0.0)?(0.5/sigma_M2_M1/sigma_M2_M1):1.0)};
    		double [] resolutions;
    		double bestW=0.0;
    		int bestSample=0;
			for (int nSample=0;nSample<this.history.size();nSample++){
				int [] motors=       this.history.get(nSample).motorsPos;
				data[nSample][0]=motorsToOrtho(motors );
				data[nSample][1]=new double[3];
				data[nSample][2]=new double[1];
				resolutions= this.history.get(nSample).getCenterResolutions();
				data[nSample][1][0]=getLensDistance(
						resolutions, // {R-sharpness,G-sharpness,B-sharpness}
						true, // boolean absolute, // return absolutely calibrated data
						lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
						lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
						1 //debugLevel
				)-targetMicrons;
				double [] averageMetrics= this.history.get(nSample).getMetrics(weightRatioRedToGreen,weightRatioBlueToGreen)[3]; // 4-th line - weighted average for all colors (TODO: do averaging here!!)
				data [nSample][1][1]=averageMetrics[1];
				data [nSample][1][2]=averageMetrics[2];
				data [nSample][2][0]=1.0;
				if (micronsFade>0) data[nSample][2][0]*=Math.exp(-kF*data[nSample][1][0]*data[nSample][1][0]);
				if (tiltFade>0)    data[nSample][2][0]*=Math.exp(-kT*(data[nSample][1][1]*data[nSample][1][1]+data[nSample][1][2]*data[nSample][1][2]));
				if (data[nSample][2][0]>bestW){
					bestW=data[nSample][2][0];
					bestSample=nSample;
				}
				if (debugLevel>2){
					double [] oMotors=motorsToOrtho(motors );
					System.out.println ("#"+nSample+": ["+motors[0]+","+motors[1]+","+motors[2]+"]-> "+
							IJ.d2s(oMotors[0],1)+","+IJ.d2s(oMotors[1],1)+","+IJ.d2s(oMotors[2],1)); //+" == "+
//							IJ.d2s(data[nSample][0][0],1)+","+IJ.d2s(data[nSample][0][1],1)+","+IJ.d2s(data[nSample][0][2],1));
				}

			}
			if (!useTheBest) bestSample=last;
    		int [] bestPosition=getPosition(bestSample); // either last or best position 

			// Additionally decay by motor positions 
			double [] orthoBest=data[bestSample][0];
			for (int nSample=0;nSample<this.history.size();nSample++){
				double [] ortho=data[nSample][0].clone();
				double w=data[nSample][2][0]; // just for debugging
				for (int i=0;i<3;i++){
					ortho[i]-=orthoBest[i];
					if (kP[i]>0) data[nSample][2][0]*=Math.exp(-kP[i]*ortho[i]*ortho[i]);
				}
				if (debugLevel>2){
					System.out.println ("##"+nSample+": ["+IJ.d2s(data[nSample][0][0],1)+","+IJ.d2s(data[nSample][0][1],1)+","+IJ.d2s(data[nSample][0][2],1)+"] "+
							" fDistDiff="+IJ.d2s(data[nSample][1][0],3)+" tiltX="+IJ.d2s(data[nSample][1][1],3)+" tiltY="+IJ.d2s(data[nSample][1][2],3)+
							" weight="+data[nSample][2][0]+" ("+w+")");
				}
			}			
			double [] orthoSolution=(new PolynomialApproximation( debugLevel)).linear3d(data);
//			double [] orthoSolution=(new PolynomialApproximation(4)).linear3d(data);
			int [] newPosition=motorsFromOrtho(orthoSolution); 
			if (debugLevel>1){
				System.out.println ("Ortho solution:   "+IJ.d2s(orthoSolution[0],1)+" : "+IJ.d2s(orthoSolution[1],1)+" : "+IJ.d2s(orthoSolution[2],1));
				System.out.println ("Motors solution:  "+newPosition[0]+" : "+newPosition[1]+" : "+newPosition[2]);
			}
			if (debugLevel>0){
				System.out.println ("Best known position:"+bestPosition[0]+" : "+bestPosition[1]+" : "+bestPosition[2]);
			}
			if (toleranceThreshold<=0) return newPosition; // just testing
// Scale and limit movement if needed.
			double [] moveFromBest={
					(newPosition[0]-bestPosition[0]),
					(newPosition[1]-bestPosition[1]),
					(newPosition[2]-bestPosition[2])};
    		if (scaleDownCorrection){
    			if (debugLevel>0) System.out.println("Scaling down correction as the residual error is under scaled tolerance, multiplying by "+IJ.d2s(scaleDownValue,3));
    			for (int i=0;i<moveFromBest.length;i++) moveFromBest[i]*=scaleDownValue;
    		}
    		double corrLength=Math.sqrt(moveFromBest[0]*moveFromBest[0]+moveFromBest[1]*moveFromBest[1]+moveFromBest[2]*moveFromBest[2]);
    		if (debugLevel>debugThreshold) System.out.println("corrLength="+corrLength+" maxStep="+maxStep);
    		if (corrLength>maxStep){
    			if (debugLevel>0) System.out.println("Scaling down correction to match maximal allowed movement ("+maxStep+" steps), multiplying by "+IJ.d2s(maxStep/corrLength,3));
    			for (int i=0;i<moveFromBest.length;i++) moveFromBest[i]*=maxStep/corrLength;
    		}

    		if (debugLevel>debugThreshold) System.out.println("motorCorrections="+IJ.d2s(moveFromBest[0],0)+", "+IJ.d2s(moveFromBest[1],0)+", "+IJ.d2s(moveFromBest[2],0));
    		int [] iNewPos=new int [moveFromBest.length];
    		for (int i=0;i<moveFromBest.length;i++)  iNewPos[i]=(int) Math.round(moveFromBest[i])+bestPosition[i];
    		if (debugLevel>debugThreshold) System.out.println("new position="+iNewPos[0]+", "+iNewPos[1]+", "+iNewPos[2]);
    		return iNewPos;
     	}

     	public double [] motorsToOrtho(int [] motors ){
     		double [] ortho={
     				((motors[0]+motors[1]+motors[2])/3.0),
     				(motors[2]-0.5*(motors[0]+motors[1])),
     				(motors[1]-motors[0])};
     		return ortho;

     	}

     	public int [] motorsFromOrtho(double [] ortho ){
     		int [] motors={
     				(int) Math.round(ortho[0] -ortho[1]/3.0 -ortho[2]/2.0),
     				(int) Math.round(ortho[0] -ortho[1]/3.0 +ortho[2]/2.0),
     				(int) Math.round(ortho[0] +2*ortho[1]/3.0)};
     		return motors;
     	}
     	
     	
     	public int [] solveFocusing(
    			double weightRatioRedToGreen,
    			double weightRatioBlueToGreen,
    			double lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
				boolean useRadialTangential,
    			double targetFarNear,
				double targetMicrons,
				double toleranceMicrons,
				double toleranceTilt,
				double toleranceThreshold,
    			double motorsSigma,
    			double believeLast,   // 0 - use 'honest' best fit, 1.0 - make each plane go through the last sample 
    			double maxStep,        // Maximal motor move for each step
    			double goodDistanceSigma, // to use only samples with small distance errors (<=0.0 - do not use)

    			int debugLevel
    	){
     		int debugThreshold=0;
    		if (this.history.size()<4) {
    			String msg="Failed: too few samples ("+this.history.size()+") to solve";
    			System.out.println(msg);
    			IJ.showMessage("Error",msg);
    			return null;
    		}
//    		double scaleDownValue=0.5; // make a parameter? multiply correction by this if residual error is smaller than scaled tolerance
    		
    		if (!useRadialTangential){
    			// Make sure calibration is done and targetMicrons is within limits
    			double [] limits=absoluteCalibrationLimits();
    			if (limits==null){
        			String msg="Failed: lens focus quality per colors vs. distance is not yet measured/calibrated";
        			System.out.println(msg);
        			IJ.showMessage("Error",msg);
        			return null;
    			}
    			if (Double.isNaN(limits[0]) || Double.isNaN(limits[1])){
        			String msg="Failed: lens focus quality per colors vs. distance is measured, but absolute calibration was not performed.";
        			System.out.println(msg);
        			IJ.showMessage("Error",msg);
        			return null;
    			}
    			if ((targetMicrons<limits[0]) || (targetMicrons>limits[1])){
        			String msg="Failed: requested distance ("+targetMicrons+"um) is not withing calibrated range "+limits[0]+".."+limits[1]+"um.";
        			System.out.println(msg);
        			IJ.showMessage("Error",msg);
        			return null;
    			}
    		}
    		double targetValue=useRadialTangential?targetFarNear:targetMicrons;

    		int [] lastPosition=this.history.get(this.history.size()-1).motorsPos;
    		double [] lastAverageMetrics=this.history.get(this.history.size()-1).getMetrics(weightRatioRedToGreen,weightRatioBlueToGreen)[3];
    		double [] lastResolutions= this.history.get(this.history.size()-1).getCenterResolutions();

    		if (!useRadialTangential){
    			lastAverageMetrics[0]=getLensDistance(
    					lastResolutions, // {R-sharpness,G-sharpness,B-sharpness}
		    			true, // boolean absolute, // return absolutely calibrated data
		    			lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
		    			lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
		    			debugLevel
		    			);
    		}
    		double  worstOver= worstOverTolerance (
        			weightRatioRedToGreen,
        			weightRatioBlueToGreen,
        			lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
        			lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    				useRadialTangential,
        			targetFarNear,
    				targetMicrons,
    				toleranceMicrons,
    				toleranceTilt,
        			debugLevel);
	        	if (debugLevel>debugThreshold) {
	        		System.out.println("solveFocusing() worstOver="+IJ.d2s(worstOver,5)+", m1="+lastPosition[0]+"m2="+lastPosition[1]+"m3="+lastPosition[2]);
	        	}

    		if (worstOver<1.0) return lastPosition; // no correction if worst is less than tolerance
    		boolean scaleDownCorrection=(worstOver<toleranceThreshold); 
    		double scaleDownValue=scaleDownCorrection?1.0:(worstOver/toleranceThreshold);
    		/*
    		 * Approximating each of the 3 measured parameters (Far/near, tilt x and tilt y) with linear approximation in the vicinity of the last position
    		 * For each parameter F(x,y,z)=A*x + B*y +C*z + D, using Gaussian weight function with sigma= motorsSigma     	    		
    		 */
    		double kExp=0.5/motorsSigma/motorsSigma;
    		double kGoodDistance=(goodDistanceSigma>0.0)?0.5/goodDistanceSigma/goodDistanceSigma:1.0;
    		double [][] aM3=new double [3][3]; // matrix with lines of A, B, C for each parameter, first line - far/nea, second - tiltX, 3-rd tiltY
    		double [][] aB3=new double [3][1];    // right side:  {targetFarNear-D(Far/Near), -D(tiltX), -D(tiltY)}
    		for (int parNum=0;parNum<aM3.length;parNum++){
    			double S0=0.0,SX=0.0,SY=0.0,SZ=0.0,SXY=0.0,SXZ=0.0,SYZ=0.0,SX2=0.0,SY2=0.0,SZ2=0.0,SF=0.0,SFX=0.0,SFY=0.0,SFZ=0.0;
    			for (int nSample=0;nSample<this.history.size();nSample++){
    				int [] motors=       this.history.get(nSample).motorsPos;
    				double [] averageMetrics= this.history.get(nSample).getMetrics(weightRatioRedToGreen,weightRatioBlueToGreen)[3]; // 4-th line - weighted average for all colors (TODO: do averaging here!!)
    	    		double [] resolutions= this.history.get(nSample).getCenterResolutions();
    				double [] xyz={
    						motors[0]-lastPosition[0],
    						motors[1]-lastPosition[1],
    						motors[2]-lastPosition[2]};
    				double w=Math.exp(-kExp*(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]));
    	    		if ((!useRadialTangential) && (parNum>0)){ // do not apply this to focal distance, only to tilts
    	    			averageMetrics[0]=getLensDistance(
    	    					resolutions, // {R-sharpness,G-sharpness,B-sharpness}
    			    			true, // boolean absolute, // return absolutely calibrated data
    			    			lensDistanceWeightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			    			lensDistanceWeightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			    			debugLevel
    			    			);
    	    			if (goodDistanceSigma>0.0) {
    	    				double diffDistance=(averageMetrics[0]-targetValue);
    	    				w*=Math.exp(-kGoodDistance*diffDistance*diffDistance); // reduce weight of tilt measurements if di8stance is wrong
    	    			}
    	    		}
    	    		//TODO replace with PolynomialApproximation class 
    				S0+=w;
    				SX+=w*xyz[0];
    				SY+=w*xyz[1];
    				SZ+=w*xyz[2];
    				SXY+=w*xyz[0]*xyz[1];
    				SXZ+=w*xyz[0]*xyz[2];
    				SYZ+=w*xyz[1]*xyz[2];
    				SX2+=w*xyz[0]*xyz[0];
    				SY2+=w*xyz[1]*xyz[1];
    				SZ2+=w*xyz[2]*xyz[2];
    				SF+=w*averageMetrics[parNum];
    				SFX+=w*averageMetrics[parNum]*xyz[0];
    				SFY+=w*averageMetrics[parNum]*xyz[1];
    				SFZ+=w*averageMetrics[parNum]*xyz[2];
    			}
    			double [][] aM={
    					{SX2,SXY,SXZ,SX},
    					{SXY,SY2,SYZ,SY},
    					{SXZ,SYZ,SZ2,SZ},
    					{SX, SY, SZ, S0}};
    			double [][] aB={
    					{SFX},
    					{SFY},
    					{SFZ},
    					{SF}};
    			Matrix M=new Matrix(aM);
    			Matrix B=new Matrix(aB);
    			// Check for singular (near-singular) here
    			double [] abcd= M.solve(B).getColumnPackedCopy();
    			if (debugLevel>debugThreshold+1){
    				System.out.println(parNum+"M:");
    				M.print(10, 5);
    				System.out.println(parNum+"B:");
    				B.print(10, 5);
    				System.out.println(parNum+"A:");
    				M.solve(B).print(10, 5);
    			}
    			//believeLast
    			aM3[parNum][0]= abcd[0];
    			aM3[parNum][1]= abcd[1];
    			aM3[parNum][2]= abcd[2];
    			aB3[parNum][0]=believeLast*lastAverageMetrics[parNum] -abcd[3]*(1.0-believeLast);

    		}
    		aB3[0][0]+=targetValue;
    		Matrix M3=new Matrix(aM3);
    		Matrix B3=new Matrix(aB3);
    		double [] motorCorrections=M3.solve(B3).getColumnPackedCopy();
    		if (debugLevel>debugThreshold) {
    			System.out.println("M3:");
    			M3.print(10, 5);
    			System.out.println("B3:");
    			B3.print(10, 5);
    			System.out.println("A3:");
    			M3.solve(B3).print(10, 5);
    		}
    		
    		//Limit maximal step (proportionally?)
    		// Show correction here 
    		if (scaleDownCorrection){
    			if (debugLevel>0) System.out.println("Scaling down correction as the residual error is under scaled tolerance");
    			for (int i=0;i<motorCorrections.length;i++) motorCorrections[i]*=scaleDownValue;
    		}
    		
    		double corrLength=Math.sqrt(motorCorrections[0]*motorCorrections[0]+motorCorrections[1]*motorCorrections[1]+motorCorrections[2]*motorCorrections[2]);
    		if (debugLevel>debugThreshold) System.out.println("corrLength="+corrLength+" maxStep="+maxStep);
    		if (corrLength>maxStep){
    			// Show that correction is reduced
    			for (int i=0;i<motorCorrections.length;i++) motorCorrections[i]*=maxStep/corrLength;
    		}
    		if (debugLevel>debugThreshold) System.out.println("motorCorrections="+IJ.d2s(motorCorrections[0],0)+", "+IJ.d2s(motorCorrections[1],0)+", "+IJ.d2s(motorCorrections[2],0));
    		int [] iNewPos=new int [motorCorrections.length];
    		for (int i=0;i<motorCorrections.length;i++)  iNewPos[i]=(int) Math.round(motorCorrections[i])+lastPosition[i];
    		if (debugLevel>debugThreshold) System.out.println("new position="+iNewPos[0]+", "+iNewPos[1]+", "+iNewPos[2]);
    		return iNewPos;
    	}
    	// todo - add "probe around" - 6/3 points, +/- for each direction (fraction of sigma?)    	    	

    	public void list(
    			boolean useLMA,
    			String path,
    			String lensSerial, // if null - do not add average
        		String comment,
    			boolean showIndividualComponents,
    			boolean showIndividualSamples,
    			double weightRatioRedToGreen,
    			double weightRatioBlueToGreen,
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY){ // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    		list (
    				useLMA,
    				path,
    				"", // serial; number
    				lensSerial, // if null - do not add average
    				comment,
    				showIndividualComponents,
    				showIndividualSamples,
    				weightRatioRedToGreen,
    				weightRatioBlueToGreen,
    				weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    				weightY,
    				false,
    				Double.NaN,
    				Double.NaN,
    				Double.NaN,
    				Double.NaN);
    	}
 /*   	
    	public void addCurrentHistory(
    			FocusingField focusingField)
    	{
    		for (int i=0;i<this.history.size();i++){
    			FocusingState focusingState=this.history.get(i);
    			focusingField.addSample(
    					focusingState.getTimestamp(),
    					focusingState.getTemperature(),
    					focusingState.motorsPos,
    					focusingState.getSamples());
    		}        	 

    	}
    	
*/
    	public void saveXML(
    			String path,
        		String serialNumber,
    			String lensSerial, // if null - do not add average
        		String comment,
    			double pX0,
    			double pY0,
    			double [][][] sampleCoord){ // x,y,r
    		String sep=" "; //",";
    		ArrayList<String> nodeList=new ArrayList<String>(); 
        	XMLConfiguration hConfig=new XMLConfiguration();
        	hConfig.setRootElementName("focusingHistory");
        	if (comment!=null)      hConfig.addProperty("comment",comment);
        	if (serialNumber!=null) hConfig.addProperty("serialNumber",serialNumber);
        	if (lensSerial!=null)   hConfig.addProperty("lensSerial",lensSerial);
        	 hConfig.addProperty("lens_center_x",pX0);
        	 hConfig.addProperty("lens_center_y",pY0);
        	 if ((sampleCoord!=null) && (sampleCoord.length>0) && (sampleCoord[0] != null) && (sampleCoord[0].length>0)){
        		 hConfig.addProperty("samples_x",sampleCoord[0].length);
        		 hConfig.addProperty("samples_y",sampleCoord.length);
        		 for (int i=0;i<sampleCoord.length;i++)
        			 for (int j=0;j<sampleCoord[i].length;j++){
//                		 double coord[] = {sampleCoord[i][j][0],sampleCoord[i][j][1]};
                		 hConfig.addProperty("sample_"+i+"_"+j,sampleCoord[i][j][0]+sep+sampleCoord[i][j][1]);
        			 }
        	 }
        	 hConfig.addProperty("measurements",this.history.size());
        	 for (int i=0;i<this.history.size();i++){
        		 FocusingState focusingState=this.history.get(i);
        		 String prefix="measurement_"+i+".";
        		 if (focusingState.getTimestamp()!=null)
        			 hConfig.addProperty(prefix+"timestamp",focusingState.getTimestamp());
        		 hConfig.addProperty(prefix+"temperature",focusingState.getTemperature());
        		 hConfig.addProperty(prefix+"motors",focusingState.motorsPos[0]+sep+focusingState.motorsPos[1]+sep+focusingState.motorsPos[2]); // array of 3
        		 double [][][][] samples=focusingState.getSamples();
        		 nodeList.clear();
        		 if ((samples!=null) && (samples.length>0) && (samples[0]!=null) && (samples[0].length>0)){
        			 for (int ii=0;ii<samples.length;ii++) for (int jj=0;jj<samples[ii].length;jj++){
        				 if ((samples[ii][jj]!=null) && (samples[ii][jj].length>0)) {
        					 String sdata=ii+sep+jj;
        					 for (int cc=0;cc<samples[ii][jj].length;cc++) {
/*        						 
        						 double rt=samples[ii][jj][cc][0]/Math.sqrt((1+samples[ii][jj][cc][1]*samples[ii][jj][cc][1])/2.0); // tangential;
        						 double rs=rt*samples[ii][jj][cc][1]; // saggital
        						 sdata += sep+rs+ // saggital
        								  sep+rt; // tangential
*/        						 
        						 sdata += sep+samples[ii][jj][cc][0]+ // x2
        								 sep+samples[ii][jj][cc][1]+  // y2
        								 sep+samples[ii][jj][cc][2];  // xy
//        						 double [] samples_st={
//        								 samples[ii][jj][cc][0], // saggital
//        								 samples[ii][jj][cc][0]/samples[ii][jj][cc][1] // tangential (was ratio)
//        						 };
        					 }
							 nodeList.add(sdata);
        				 }
        			 }
					 hConfig.addProperty(prefix+"sample",nodeList);
        		 }        	 
        	 }
         	File file=new File (path);
         	BufferedWriter writer;
 			try {
 				writer = new BufferedWriter(new FileWriter(file));
 	        	hConfig.save(writer);
 			} catch (IOException e) {
 				// TODO Auto-generated catch block
 				e.printStackTrace();
 			} catch (ConfigurationException e) {
 				// TODO Auto-generated catch block
 				e.printStackTrace();
 			}
        	 
    	}
     	public void list(
        		boolean useLMA,
    			String path,
        		String serialNumber,
    			String lensSerial, // if null - do not add average
        		String comment,
    			boolean showIndividualComponents,
    			boolean showIndividualSamples,
    			double weightRatioRedToGreen,
    			double weightRatioBlueToGreen,
    			double weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    			double weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    			boolean justSummary,
    			double pX0,
    			double pY0,
    			double result_lastKT,
    			double result_allHistoryKT){
    		String shortHeader="timestamp\ttemperature\ttravel\tcenter\tMotor1\tMotor2\tMotor3\tDist\tFar/Near\tTilt X\tTilt Y\tR50%\tA50%\tB50%\tC-R50%\tResol\tResol-R\tResol-G\tResol-B\tC-Resol\tC-Resol-R\tC-Resol-G\tC-Resol-B";
    		String fullHeader= "timestamp\ttemperature\ttravel\tcenter\tdir1\tMotor1\tdir2\tMotor2\tdir3\tMotor3\tDist\tComponent\tFar/Near\tTilt X\tTilt Y\tR50%\tA50%\tB50%\tC-R50%\tResol\tC-Resol\tWeight(%)";
    		String header=showIndividualComponents?fullHeader:shortHeader;
    		if (justSummary){
    			header="Serial\tLens\t"+header+"\tCenter X\tCenter Y\tLast KT\tAll KT\tComment";
    		} else {
    			header="#\t"+header;
    		}
    		int sampleRows=0;
    		int sampleColumns=0;
    		if (showIndividualSamples && !justSummary) {
    			for (int i=0;i<this.history.size();i++){
    				FocusingState focusingState=this.history.get(i);
    				double [][][][] samples=focusingState.getSamples();
    				if (samples!=null){
    					sampleRows=samples.length;
    					sampleColumns=samples[0].length;
    					break;
    				}
    			}
    		}
    		if ((sampleRows==0) || (sampleColumns==0)) showIndividualSamples=false;
    		if (showIndividualSamples){
    			if (showIndividualComponents) { // additional 3 per-color rows
    				for (int i=0;i<sampleRows;i++) for (int j=0;j<sampleColumns;j++){
//    					header+="\tY"+i+"X"+j+"_R"+"\tY"+i+"X"+j+"_R/T";
    					header+="\tY"+i+"X"+j+"_X2"+"\tY"+i+"X"+j+"_Y2"+"\tY"+i+"X"+j+"_XY";
    				}

    			} else {
    				for (int i=0;i<sampleRows;i++) for (int j=0;j<sampleColumns;j++){ // single row fro all colors
//    					header+="\tY"+i+"X"+j+"_Red_R"+"\tY"+i+"X"+j+"_Red_R/T"+
//    					"\tY"+i+"X"+j+"_Green_R"+"\tY"+i+"X"+j+"_Green_R/T"+
//    					"\tY"+i+"X"+j+"_Blue_R"+"\tY"+i+"X"+j+"_Blue_R/T";
    					header+="\tY"+i+"X"+j+"_Red_X2"+"\tY"+i+"X"+j+"_Red_Y2"+"\tY"+i+"X"+j+"_Red_XY"+
    					"\tY"+i+"X"+j+"_Green_X2"+"\tY"+i+"X"+j+"_Green_Y2"+"\tY"+i+"X"+j+"_Green_XY"+
    					"\tY"+i+"X"+j+"_Blue_X2"+"\tY"+i+"X"+j+"_Blue_Y2"+"\tY"+i+"X"+j+"_Blue_XY";
    				}
    			}
    		}
    		StringBuffer sb = new StringBuffer();
    		int [] prevPos=null;
    		if (showIndividualComponents) {
    			for (int i=0;i<this.history.size();i++){
    				FocusingState focusingState=this.history.get(i);
    				double [] zTxTy=useLMA?focusingState.getzTxTy():null;
    				double [][] metrics=focusingState.getMetrics(weightRatioRedToGreen,weightRatioBlueToGreen);
    				double [][] resolution=focusingState.getSharpness(weightRatioRedToGreen,weightRatioBlueToGreen);
    				double dist= (zTxTy==null)?getLensDistance(
    						focusingState.getCenterResolutions(), // {R-sharpness,G-sharpness,B-sharpness}
    		    			true, //boolean absolute, // return absolutely calibrated data
    		    			weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
    		    			weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
    		    			1): //int debugLevel
    		    				zTxTy[0];
    				double []   averageMetrics=metrics[3];
    				if (zTxTy!=null){
        				averageMetrics=metrics[3].clone(); // to modify w/o changing original
    					averageMetrics[1]=zTxTy[1]; // tiltX
    					averageMetrics[2]=zTxTy[2]; // tiltY
    				}
    				double []   averageResolution=resolution[3];
    				sb.append((i+1)+"\t");
    				String timestamp=focusingState.getTimestamp();
    				if (timestamp==null) timestamp="n/a";
    				sb.append(timestamp+"\t");
    				
    				double sensorTemperature=focusingState.getTemperature();
    				sb.append(IJ.d2s(sensorTemperature,2)+"\t"); 
    				if (prevPos!=null) {
    					double [] xyz={focusingState.motorsPos[0]-prevPos[0],focusingState.motorsPos[1]-prevPos[1],focusingState.motorsPos[2]-prevPos[2]};
    					sb.append(""+IJ.d2s(Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]),0));
    				}
    				sb.append(  "\t"+IJ.d2s((focusingState.motorsPos[0]+focusingState.motorsPos[1]+focusingState.motorsPos[2])/3.0,1));
    				sb.append(  "\t"+(focusingState.motorsDir[0]?"(+)":"(-)")+"\t"+focusingState.motorsPos[0]);
    				sb.append(  "\t"+(focusingState.motorsDir[1]?"(+)":"(-)")+"\t"+focusingState.motorsPos[1]);
    				sb.append(  "\t"+(focusingState.motorsDir[2]?"(+)":"(-)")+"\t"+focusingState.motorsPos[2]);
    				if (Double.isNaN(dist)) sb.append("\t---"); else sb.append(  "\t"+IJ.d2s(dist,3));
    				sb.append("\t"+this.focusingState.getName(3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[0],3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[1],3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[2],3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[3],3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[4],3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[5],3));
    				sb.append(  "\t"+IJ.d2s(averageMetrics[6],3));
    				sb.append(  "\t"+IJ.d2s(averageResolution[0],3));
    				sb.append(  "\t"+IJ.d2s(averageResolution[1],3));
    				if (showIndividualSamples){
    					for (int ii=0;ii<sampleRows;ii++) for (int jj=0;jj<sampleColumns;jj++){
    						sb.append(  "\t-\t-");
    					}
    				}
					sb.append(  "\n");
    				for (int c=0;c<3;c++){
    					sb.append(  "\t\t\t\t\t\t\t\t\t\t\t"+focusingState.getName(c));
    					sb.append(  "\t"+IJ.d2s(metrics[c][0],3));
    					sb.append(  "\t"+IJ.d2s(metrics[c][1],3));
    					sb.append(  "\t"+IJ.d2s(metrics[c][2],3));
    					sb.append(  "\t"+IJ.d2s(metrics[c][3],3));
    					sb.append(  "\t"+IJ.d2s(metrics[c][4],3));
    					sb.append(  "\t"+IJ.d2s(metrics[c][5],3));
    					sb.append(  "\t"+IJ.d2s(metrics[c][6],3));
    					sb.append(  "\t"+IJ.d2s(resolution[c][0],3));
    					sb.append(  "\t"+IJ.d2s(resolution[c][1],3));
    					sb.append(  "\t"+IJ.d2s(100*metrics[c][7],1)+"%");
        				if (showIndividualSamples){
        					double [][][][] samples=focusingState.getSamples();
        					for (int ii=0;ii<sampleRows;ii++) for (int jj=0;jj<sampleColumns;jj++){
        						if ((samples[ii][jj]!=null) && (samples[ii][jj][c]!=null)) {
        							sb.append(  "\t"+IJ.d2s(samples[ii][jj][c][0],3));
        							sb.append(  "\t"+IJ.d2s(samples[ii][jj][c][1],3));
        							sb.append(  "\t"+IJ.d2s(samples[ii][jj][c][2],3));
        						} else {
        							sb.append(  "\t---\t---\t---");
        						}
        					}
        				}
    					sb.append(  "\n");
    				}
    				prevPos=focusingState.motorsPos;
    			}
    		} else { // calculate and add average values (not yet for individual above)
    			int numPars=23;
    			int numDistMeas=0;
    			double [] sums= new double [numPars];
    			for (int i=0;i<sums.length;i++) sums[i]=0.0;
    			for (int i=0;i<this.history.size();i++){
//    				int parIndex=0;
    				FocusingState focusingState=this.history.get(i);
    				double [] zTxTy=useLMA?focusingState.getzTxTy():null;
    				double [][] metrics=focusingState.getMetrics(weightRatioRedToGreen,weightRatioBlueToGreen);
    				double [][] resolution={{Double.NaN,Double.NaN},{Double.NaN,Double.NaN},{Double.NaN,Double.NaN},{Double.NaN,Double.NaN}};
    				try {
    					resolution=focusingState.getSharpness(weightRatioRedToGreen,weightRatioBlueToGreen);
    				} catch (Exception e){
    					System.out.println("Failed to get resolution for history("+i+")");
    					continue; // skip the whole line
    				}
    				double []   averageMetrics=metrics[3];
    				if (zTxTy!=null){
        				averageMetrics=metrics[3].clone(); // to modify w/o changing original
    					averageMetrics[1]=zTxTy[1]; // tiltX
    					averageMetrics[2]=zTxTy[2]; // tiltY
    				}
    				double []   averageResolution=resolution[3];
    				if (!justSummary) sb.append((i+1)+"\t");
    				String timestamp=focusingState.getTimestamp();
    				if (timestamp==null) timestamp="n/a";
    				else sums[0]+=Double.parseDouble(timestamp); //0
    				if (!justSummary) sb.append(timestamp+"\t");
    				double sensorTemperature=focusingState.getTemperature(); 
    				if (!justSummary) sb.append(IJ.d2s(sensorTemperature,2)+"\t");
    				sums[1]+=sensorTemperature;
    				if (prevPos!=null) {
    					double [] xyz={focusingState.motorsPos[0]-prevPos[0],focusingState.motorsPos[1]-prevPos[1],focusingState.motorsPos[2]-prevPos[2]};
    					if (!justSummary) sb.append(""+IJ.d2s(Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]),0));
    					sums[2]+=Math.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);//2
    				}
    				if (!justSummary) sb.append(  "\t"+IJ.d2s((focusingState.motorsPos[0]+focusingState.motorsPos[1]+focusingState.motorsPos[2])/3.0,1));
    				sums[3]+=(focusingState.motorsPos[0]+focusingState.motorsPos[1]+focusingState.motorsPos[2])/3.0; // 3
    				if (!justSummary) {
    					sb.append(  "\t"+focusingState.motorsPos[0]);
    					sb.append(  "\t"+focusingState.motorsPos[1]);
    					sb.append(  "\t"+focusingState.motorsPos[2]);
    				}
					sums[4]+=focusingState.motorsPos[0];
					sums[5]+=focusingState.motorsPos[1];
					sums[6]+=focusingState.motorsPos[2];
					double dist= (zTxTy==null)?getLensDistance(
							focusingState.getCenterResolutions(), // {R-sharpness,G-sharpness,B-sharpness}
							true, //boolean absolute, // return absolutely calibrated data
							weightK, // 0.0 - all 3 component errors are combined with the same weights. 1.0 - proportional to squared first derivatives 
							weightY, // R-frac, Y-frac have the same scale regardless of the sharpness, but not Y. This is to balance Y contribution
							1): //int debugLevel
								zTxTy[0];
    				if (Double.isNaN(dist)){
    					if (!justSummary) sb.append("\t---");
    				} else {
    					if (!justSummary) sb.append(  "\t"+IJ.d2s(dist,3));
    					numDistMeas++;
    					sums[7]+=dist; //may be NaN
    				}
    				if (!justSummary)  {
    					sb.append(  "\t"+IJ.d2s(averageMetrics[0],3));
    					sb.append(  "\t"+IJ.d2s(averageMetrics[1],3));
    					sb.append(  "\t"+IJ.d2s(averageMetrics[2],3));
    					sb.append(  "\t"+IJ.d2s(averageMetrics[3],3));
    					sb.append(  "\t"+IJ.d2s(averageMetrics[4],3));
    					sb.append(  "\t"+IJ.d2s(averageMetrics[5],3));
    					sb.append(  "\t"+IJ.d2s(averageMetrics[6],3));
    					sb.append(  "\t"+IJ.d2s(averageResolution[0],3));
    					sb.append(  "\t"+IJ.d2s(resolution[0][0],3));
    					sb.append(  "\t"+IJ.d2s(resolution[1][0],3));
    					sb.append(  "\t"+IJ.d2s(resolution[2][0],3));
    					sb.append(  "\t"+IJ.d2s(averageResolution[1],3));
    					sb.append(  "\t"+IJ.d2s(resolution[0][1],3));
    					sb.append(  "\t"+IJ.d2s(resolution[1][1],3));
    					sb.append(  "\t"+IJ.d2s(resolution[2][1],3));
        				if (showIndividualSamples){
        					double [][][][] samples=focusingState.getSamples();
        	    			for (int ii=0;ii<sampleRows;ii++) for (int jj=0;jj<sampleColumns;jj++){
        	    				for (int cc=0;cc<3;cc++) {
        	    					if ((samples[ii][jj]!=null) && (samples[ii][jj][cc]!=null)) {
        	        					sb.append(  "\t"+IJ.d2s(samples[ii][jj][cc][0],3));
        	        					sb.append(  "\t"+IJ.d2s(samples[ii][jj][cc][1],3));
        	        					sb.append(  "\t"+IJ.d2s(samples[ii][jj][cc][2],3));
        	    					} else {
        	    						sb.append(  "\t---\t---\t---");
        	    					}
        	    				}
        	    			}
        				}
    					sb.append(  "\n");
    				}
					sums[ 8]+=averageMetrics[0];
					sums[ 9]+=averageMetrics[1];
					sums[10]+=averageMetrics[2];
					sums[11]+=averageMetrics[3];
					sums[12]+=averageMetrics[4];
					sums[13]+=averageMetrics[5];
					sums[14]+=averageMetrics[6];
					sums[15]+=averageResolution[0];
					sums[16]+=resolution[0][0];
					sums[17]+=resolution[1][0];
					sums[18]+=resolution[2][0];
					sums[19]+=averageResolution[1];
					sums[20]+=resolution[0][1];
					sums[21]+=resolution[1][1];
					sums[22]+=resolution[2][1];
    				prevPos=focusingState.motorsPos;
    			}
    			String averageData="";
    			// how averages
    			for (int i=0;i<sums.length;i++) sums[i]/=(i==7)?numDistMeas:this.history.size();
    			averageData+=String.format("LENS%S",lensSerial)+"\t";
    			averageData+=IJ.d2s(sums[ 0],3)+"\t"; //timestamp
    			averageData+=IJ.d2s(sums[ 1],2)+"\t"; // temperature
    			averageData+=IJ.d2s(sums[ 2],0)+"\t"; // travel
    			averageData+=IJ.d2s(sums[ 3],0)+"\t"; // center
    			averageData+=IJ.d2s(sums[ 4],0)+"\t"; // motorsPos[0];
    			averageData+=IJ.d2s(sums[ 5],0)+"\t"; // motorsPos[1];
    			averageData+=IJ.d2s(sums[ 6],0)+"\t"; // motorsPos[2];
    			averageData+=((numDistMeas>0)?IJ.d2s(sums[ 7],3):"---")+"\t"; // distance;
    			averageData+=IJ.d2s(sums[ 8],3)+"\t"; // averageMetrics[0];
    			averageData+=IJ.d2s(sums[ 9],3)+"\t"; // averageMetrics[1];
    			averageData+=IJ.d2s(sums[10],3)+"\t"; // averageMetrics[2];
    			averageData+=IJ.d2s(sums[11],3)+"\t"; // averageMetrics[3];
    			averageData+=IJ.d2s(sums[12],3)+"\t"; // averageMetrics[4];
    			averageData+=IJ.d2s(sums[13],3)+"\t"; // averageMetrics[5];
    			averageData+=IJ.d2s(sums[14],3)+"\t"; // averageMetrics[6];
    			averageData+=IJ.d2s(sums[15],3)+"\t"; // averageResolution[0];
    			averageData+=IJ.d2s(sums[16],3)+"\t"; // resolution[0][0];
    			averageData+=IJ.d2s(sums[17],3)+"\t"; // resolution[1][0];
    			averageData+=IJ.d2s(sums[18],3)+"\t"; // resolution[2][0];
    			averageData+=IJ.d2s(sums[19],3)+"\t"; // averageResolution[1];
    			averageData+=IJ.d2s(sums[20],3)+"\t"; // resolution[0][1];
    			averageData+=IJ.d2s(sums[21],3)+"\t"; // resolution[1][1];
    			averageData+=IJ.d2s(sums[22],3); // resolution[2][1];
    			if (justSummary){
    				averageData=serialNumber+"\t"+averageData;
        			averageData+="\t"+ IJ.d2s(pX0,1); // 
        			averageData+="\t"+ IJ.d2s(pY0,1); //
        			averageData+="\t"+ IJ.d2s(result_lastKT,3);
        			averageData+="\t"+ IJ.d2s(result_allHistoryKT,3); 
        			averageData+="\t"+ comment; // 
    			}
    			averageData+="\n";
    			// TODO: add data for individula samples?
    			sb.append(averageData);
    			if (path!=null){
    				System.out.println("Averaging measurements:");
    				System.out.println(header);
    				System.out.println(averageData);
    			}
    		}
    		if (path!=null) {
    			String footer= justSummary?"":("\nComment: "+comment+"\tdX=\t"+ IJ.d2s(pX0,1)+"\tdY=\t"+IJ.d2s(pY0,1));
    			CalibrationFileManagement.saveStringToFile (
						path,
						header+"\n"+sb.toString()+footer);
    		}else{
        		new TextWindow("History", header, sb.toString(), showIndividualComponents?1500:1700,500);
    		}
    	}

    	public void setLastProbed(){
    		if (this.history.size()==0){
    			String msg="List is empty";
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		this.history.get(this.history.size()-1).isProbed=true;

    	}

    	public double distFromProbed(){
    		if (this.history.size()>0) return distFromProbed (null);
    		else return Double.NaN; 
    	}
    	public double distFromProbed(int [] position){
    		double minDist=Double.NaN;
    		if (position==null) position= this.history.get(this.history.size()-1).motorsPos;
    		if (this.history.size()>0){
    			for (int i=0;i<this.history.size();i++) if (this.history.get(i).isProbed) {
    				int [] probedPosition=this.history.get(i).motorsPos;
    				double dist=Math.sqrt(
    						(position[0]-probedPosition[0])*(position[0]-probedPosition[0])+
    						(position[1]-probedPosition[1])*(position[1]-probedPosition[1])+
    						(position[2]-probedPosition[2])*(position[2]-probedPosition[2]));
    				if (!(minDist<dist)){ //minDist may be NaN
    					minDist=dist;
    				}
    			}
    		}
    		return minDist; 
    	}

    	
    	public void clear(){
    		this.history.clear();
    		this.currentSigma=Double.NaN;
    		this.calibrationHistoryFromTo=null; // history indices, to-from (inclusive) used to calculate calibration 
    		this.hysteresisHistoryFromTo=null; // history indices, to-from (inclusive) used to calculate hysteresis
    	}
    	public void add(
    			String sTimestamp,
    			double temperature,
    			int [] motors, // FocusingMotors.curpos
    			boolean [] lastDirs, // FocusingMotors.lastDirs
    			double[][] psfMetrics, // null OK
    			double [][][][] fullResults // null OK
    	){
    		this.history.add(new FocusingState(sTimestamp,temperature,motors,lastDirs,psfMetrics,fullResults));
    	}
    	/*
    	      	public void add(
    			String sTimestamp,
    			double temperature,
    			int [] motors, // FocusingMotors.curpos
    			boolean [] lastDirs, // FocusingMotors.lastDirs
    			double[][] psfMetrics // null OK
    	){
    		this.history.add(new FocusingState(sTimestamp,temperature,motors,lastDirs,psfMetrics));

    	}
        */
    	
    	
    	public void setLastMetrics(
    			double[][] psfMetrics
    	){
    		if (this.history.size()==0){
    			String msg="List is empty";
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		this.history.get(this.history.size()-1).setMetrics(psfMetrics);
    	}


    	public class FocusingState{
    		public double temperature;
    		public int [] motorsPos=new int[3];
    		public boolean [] motorsDir=new boolean[3]; // true if was increasing
    		private double [][] psfMetricses=null;
    		boolean isProbed=false;
    		private int [] indices={1,5,2,6};
    		private String [] names={"red","green","blue","AVERAGE"};
    		public int indexFarNear=0;
    		public int indexTiltX=  1;
    		public int indexTiltY=  2;
    		public int indexR50All= 3;
    		public int indexA50All= 4;
    		public int indexB50All= 5;
    		public int indexR50Center= 6;
//    		public int indexA50Center= 7;
//    		public int indexB50Center= 8;
    		public String sTimestamp=null;
    		public double [][][][] fullResults=null;

    		public FocusingState(){}
    		public FocusingState(
    				String sTimestamp,
    				double temperature,
    				int [] motorsPos,
    				boolean [] motorsDir,
    				double [][] psfMetrics,
    				double [][][][] fullResults){
    			this.temperature=temperature;
    			this.motorsPos=motorsPos.clone();
    			this.motorsDir=motorsDir.clone();
    			this.fullResults=null;
    			if (fullResults!=null){ // deep clone
    				this.fullResults=new double [fullResults.length][][][];
    				for (int i=0;i<fullResults.length;i++){
    					if (fullResults[i]!=null){
    						this.fullResults[i]=new double [fullResults[i].length][][];			
    						for (int j=0;j<fullResults[i].length;j++){
    	    					if (fullResults[i][j]!=null){
    	    						this.fullResults[i][j]=new double [fullResults[i][j].length][];
    	    						for (int c=0;c<fullResults[i][j].length;c++){
    	    							if (fullResults[i][j][c]!=null){
    	    								this.fullResults[i][j][c]=new double [fullResults[i][j][c].length]; // 2
    	    								for (int n=0;n<fullResults[i][j][c].length;n++){
    	    									this.fullResults[i][j][c][n]=fullResults[i][j][c][n];
    	    								}
    	    							} else this.fullResults[i][j][c]=null;
    	    						}
    	    					} else this.fullResults[i][j]=null;
    						}
    					} else this.fullResults[i]=null; 
    				}
    			}

    			if (psfMetrics!=null) setMetrics(psfMetrics);
    			if (sTimestamp!=null) this.sTimestamp=sTimestamp;

    		}
/*    		
    		public FocusingState(
    				String sTimestamp,
    				double temperature,
    				int [] motorsPos,
    				boolean [] motorsDir,
    				double [][] psfMetrics){
    			this.temperature=temperature;
    			this.motorsPos=motorsPos.clone();
    			this.motorsDir=motorsDir.clone();
    			this.fullResults=null;
    			if (psfMetrics!=null) setMetrics(psfMetrics);
    			if (sTimestamp!=null) this.sTimestamp=sTimestamp;

    		}
*/
    		public FocusingState(
    				String sTimestamp,
    				double temperature,
    				int [] motorsPos,
    				boolean [] motorsDir){
    			this.temperature=temperature;
    			this.motorsPos=motorsPos.clone();
    			this.motorsDir=motorsDir.clone();
    			if (sTimestamp!=null) this.sTimestamp=sTimestamp;
    		}

    		public FocusingState(
    				String sTimestamp,
    				double temperature,
    				int [] motorsPos,
    				FocusingState previousState){
    			this.temperature=temperature;
    			this.motorsPos=motorsPos.clone();
    			setDirsFromPrevious(previousState);
    			if (sTimestamp!=null) this.sTimestamp=sTimestamp;
    		}

    		public FocusingState(
    				String sTimestamp,
    				double temperature,
    				int [] motorsPos,
    				FocusingState previousState,
    				double [][] psfMetrics){
    			this.temperature=temperature;
    			this.motorsPos=motorsPos.clone();
    			setDirsFromPrevious(previousState);
    			if (psfMetrics!=null) setMetrics(psfMetrics);
    			if (sTimestamp!=null) this.sTimestamp=sTimestamp;
    		}

    		public FocusingState(
    				String sTimestamp,
    				double temperature,
    				int m1,
    				int m2,
    				int m3,
    				FocusingState previousState,
    				double [][] psfMetrics){
    			this.temperature=temperature;
    			int [] motorsPos={m1,m2,m3};
    			this.motorsPos=motorsPos.clone();
    			setDirsFromPrevious(previousState);
    			if (psfMetrics!=null) setMetrics(psfMetrics);
    			if (sTimestamp!=null) this.sTimestamp=sTimestamp;

    		}
    		public void setMetrics(double [][] psfMetrics){
    			if (psfMetrics!=null) {
    				this.psfMetricses=new double [psfMetrics.length][];
    				for (int i=0;i<psfMetrics.length;i++){
    					if (psfMetrics[i]!=null){
    						this.psfMetricses[i]=psfMetrics[i].clone();
    					} else this.psfMetricses[i]=null;
    				}
    			} else this.psfMetricses=null;

    		}
    		public void setDirsFromPrevious(
    				FocusingState previousState){
    			for (int i=0;i<motorsPos.length;i++){
    				if (previousState.motorsPos[i]==this.motorsPos[i]) this.motorsDir[i]=previousState.motorsDir[i];
    				else if (previousState.motorsPos[i]<this.motorsPos[i]) this.motorsDir[i]=true;
    				else this.motorsDir[i]=false;
    			}
    		}
    		public double getTemperature(){
    			return this.temperature;
    		}
    		public String getTimestamp(){
    			return this.sTimestamp;
    		}

    		public String getName(int i){
    			return this.names[i];
    		}

    		public int [] getMotors(){
    			return this.motorsPos.clone();
    		}
    		public double getMotorsAverage(){
    			return (this.motorsPos[0]+this.motorsPos[1]+this.motorsPos[2])/3.0;
    		}
    		public double [][] getSharpness(
    				double weightRatioRedToGreen,
    				double weightRatioBlueToGreen){
    			return getSharpness(
    					weightRatioRedToGreen,
    					1.0,
    					weightRatioBlueToGreen);
    		}
    		public double [][] getSharpness( // first index - color+average, second index 0 - all, 1 - center
    				double weightRed,
    				double weightGreen,
    				double weightBlue){
    			double [][] metrics=getMetrics(weightRed, weightGreen, weightBlue);
    			double [][]sharpness={
    					{1.0/metrics[0][this.indexR50All],1.0/metrics[0][this.indexR50Center]}, //java.lang.NullPointerException
    					{1.0/metrics[1][this.indexR50All],1.0/metrics[1][this.indexR50Center]},
    					{1.0/metrics[2][this.indexR50All],1.0/metrics[2][this.indexR50Center]},
    					{0.0,0.0}};
    			double [] weights={weightRed,weightGreen,weightBlue};
    			double w=0.0;
    			for (int i=0;i<weights.length;i++) w+=weights[i];
    			for (int i=0;i<weights.length;i++) weights[i]/=w;
    			for (int i=0;i<3;i++)for (int j=0;j<2;j++) sharpness[3][j]+=weights[i]*sharpness[i][j];
    			return sharpness;
    		}
    		public double [] getCenterResolutions(){
    			try {
    			double [] resolutions={
    					1.0/this.psfMetricses[this.indices[0]][this.indexR50Center],  // red, r50% center
    					1.0/this.psfMetricses[this.indices[1]][this.indexR50Center],  // green, r50% center
    					1.0/this.psfMetricses[this.indices[2]][this.indexR50Center]}; // blue, r50% center
    			return resolutions;
    			} catch (Exception e){
    				return null;
    			}
    		}
    		
    		public double [][] getMetrics(
    				double weightRatioRedToGreen,
    				double weightRatioBlueToGreen){
    			return getMetrics(
    					weightRatioRedToGreen,
    					1.0,
    					weightRatioBlueToGreen);
    		}
    		// alternative mode (from aberration model) ising metrics[6][]
    		public double [] getzTxTy(){
    			if (this.psfMetricses==null){
    				return null;
    			}
    			if ((this.psfMetricses.length<7) || (this.psfMetricses[6]==null)) {
    				System.out.println("BUG? psfMetrics does not have line 6 with lens berration model z, tx, ty");
    				return null;
    			}
    			double [] zTxTy={this.psfMetricses[6][0],this.psfMetricses[6][1],this.psfMetricses[6][2]};
    			return zTxTy;

    		}
    		public double [][] getMetrics(
    				double weightRed,
    				double weightGreen,
    				double weightBlue){
    			boolean [] squared={false,false,false,false,true,true,false,false};
    			double [][] metrics=new double[4][];
    			double [] weights={weightRed,weightGreen,weightBlue};
    			double w=0.0;
    			for (int i=0;i<weights.length;i++) w+=weights[i];
    			for (int i=0;i<weights.length;i++) weights[i]/=w;
    			try {
    				metrics[3]=new double[this.psfMetricses[this.indices[0]].length];
    				for (int c=0;c<3;c++) metrics[c]= this.psfMetricses[this.indices[c]];
    			} catch (Exception e){
    				return null;
    			}
    			for (int i=0; i<metrics[3].length;i++) {
    				metrics[3][i]=0.0;
    				if (squared[i]){
    					for (int c=0;c<3;c++) metrics[3][i]+=metrics[c][i]*metrics[c][i]*weights[c];
    					metrics[3][i]=Math.sqrt(metrics[3][i]);
    				} else for (int c=0;c<3;c++) metrics[3][i]+=metrics[c][i]*weights[c];
    			}
    			return metrics;
    		}
    		public double [][][][] getSamples(){
    			return this.fullResults;
    		}

    	}
    }


    }

}
