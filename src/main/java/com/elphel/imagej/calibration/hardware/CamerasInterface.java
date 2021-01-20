package com.elphel.imagej.calibration.hardware;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import com.elphel.imagej.calibration.MatchSimulatedPattern;
import com.elphel.imagej.calibration.UVLEDandLasers;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.jp4.JP46_Reader_camera;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class CamerasInterface{
//		JP46_Reader_camera JP4_INSTANCE= new JP46_Reader_camera(false);
		public LaserPointersHardware laserPointers=null;
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
        private String [] cameraIPs = null; // since nc393 port is a part of cameraIPs[]
        private int [] channelIPPort =    null; // index in camareIPs (each IP/port combination) for each individual sensor
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
		private boolean [][] sensorPresent=  null;   // probe which sensors (of 3) are detected per system board (NC393 - per board/port)
		// TODO - try if skipping setting TRIG_PERIOD=0, IRQ_SMART=6 (when they are already set) will fix hanging
		private int [] triggerPeriod=        null;
		private int [] cameraMasterPort=     null;
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
				{0,0,3},{0,0,2},{0,0,0},{0,0,1},{1,0,3},{1,0,2},{1,0,0},{1,0,1},
				{0,1,3},{0,1,2},{0,1,0},{0,1,1},{1,1,3},{1,1,2},{1,1,0},{1,1,1},
				{0,2,3},{0,2,2},{0,2,0},{0,2,1},{1,2,3},{1,2,2},{1,2,0},{1,2,1},
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

        public CamerasInterface(int size, LaserPointersHardware laserPointers){
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
			// this.triggerURL is already defined
//			this.triggerURL="http://"+this.cameraIPs[this.masterSubCamera]+":"+(this.imgsrvPort+this.masterPort)+"/"+triggerURLcmd;

			this.images=   new ImagePlus[this.channelMap.length];
			this.imagesIP= new ImagePlus[this.cameraIPs.length];

			for (int i=0; i<this.cameraIPs.length;i++){
				this.jp4_Instances[i]=new JP46_Reader_camera(false);// invisible
//				this.jp4_Instances[i].camera_url="http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/";
				this.jp4_Instances[i].camera_url="http://"+this.cameraIPs[i]+"/";
				this.jp4_Instances[i].camera_img=    this.imageURLcmd; // not currently used
				this.jp4_Instances[i].camera_img_new=this.imageURLcmd; //"torp/wait/" will survive, only "towp/wait/" is removed for Exif re-read
				this.jp4_Instances[i].ABSOLUTELY_SILENT=true;
//				this.resetURLs[i]="http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/"+resetURLcmd;
//				this.imageURLs[i]="http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/"+this.imageURLcmd;
//				this.metaURLs[i]= "http://"+this.cameraIPs[i]+":"+this.imgsrvPort+"/"+metaURLcmd;
				this.resetURLs[i]="http://"+this.cameraIPs[i]+"/"+resetURLcmd;
				this.imageURLs[i]="http://"+this.cameraIPs[i]+"/"+this.imageURLcmd;
				this.metaURLs[i]= "http://"+this.cameraIPs[i]+"/"+metaURLcmd;
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
			ArrayList<Integer> ip_ports_list= new ArrayList<Integer>();
			for (int i=0;i<this.channelMap.length;i++) {
				Integer ip_port=(this.channelMap[i][0]<<2) + this.channelMap[i][2];
				if (!ip_ports_list.contains(ip_port)) ip_ports_list.add(ip_port);
			}
			Collections.sort(ip_ports_list);
			this.cameraIPs =   new String [ip_ports_list.size()];
			this.channelIPPort = new int [this.channelMap.length];
			for (int i = 0; i<this.cameraIPs.length; i++){
				int ip_index= ip_ports_list.get(i)>>2;
				int sensor_port = ip_ports_list.get(i) & 3;
				this.cameraIPs[i] = this.cameraSubnet+(this.iBaseIP + ip_index) + ":"+ (this.imgsrvPort+sensor_port);
				for (int j = 0; j<this.channelMap.length; j++){
					if ((this.channelMap[j][0] == ip_index) && (this.channelMap[j][2] == sensor_port)) {
						this.channelIPPort[j] = i;
					}
				}
			}

			this.triggerURL="http://"+this.cameraSubnet+(this.iBaseIP+this.masterSubCamera)+":"+
				(this.imgsrvPort+ (this.masterPort & 3))+"/"+triggerURLcmd;
			if (this.debugLevel>2) System.out.println("DEBUG393: initIPs(): this.triggerURL ="+this.triggerURL);
		}
/*
//pre nc393
		private void initIPs(){
			if (this.debugLevel>2) System.out.println("initIPs(): this.iBaseIP=" + this.iBaseIP );
			int size=0;
			for (int i=0;i<this.channelMap.length;i++) if (this.channelMap[i][0]>size) size=this.channelMap[i][0];
			size++;
			this.cameraIPs=new String [size];
			for (int i=0;i<size;i++) this.cameraIPs[i]=this.cameraSubnet+(this.iBaseIP+i);
//			this.masterSubCamera=0;
		}
 */
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
    		if (properties.getProperty(prefix+"channelMap.length")!=null) {
    			// next initializes default values, so it should be before reading them from saved properties
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


// both defaults were here



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
    		initIPs(); // was missing here?
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
	   		int colon_index = this.cameraIPs[chn].indexOf(":");
	   		int sensor_port = Integer.parseInt(this.cameraIPs[chn].substring(colon_index+1)) - this.imgsrvPort;
	   		String ip=this.cameraIPs[chn].substring(0, colon_index);
//	   		String url="http://"+this.cameraIPs[chn]+"/i2c.php?cmd=fromEEPROM0&EEPROM_chn="+EEPROM_chn;
	   		String url="http://"+ip+"/i2c.php?cmd=fromEEPROM" + sensor_port+ "&EEPROM_chn="+EEPROM_chn;
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

	   	    		serial=(dom.getDocumentElement().getElementsByTagName("serial").item(0).getChildNodes().item(0)).getNodeValue();
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

	   		int colon_index = this.cameraIPs[chn].indexOf(":");
	   		int sensor_port = Integer.parseInt(this.cameraIPs[chn].substring(colon_index+1)) - this.imgsrvPort;
	   		String ip=this.cameraIPs[chn].substring(0, colon_index);
//	   		String url="http://"+this.cameraIPs[chn]+"/i2c.php?cmd=fromEEPROM0&EEPROM_chn="+EEPROM_chn;
	   		String url="http://"+ip+"/i2c.php?cmd=fromEEPROM"+ sensor_port +"&EEPROM_chn="+EEPROM_chn;
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

	   			String sTemperature=null;

	   			if ((dom.getDocumentElement().getElementsByTagName("sensorTemperature").item(0))!=null){
	   				sTemperature=(dom.getDocumentElement().getElementsByTagName("sensorTemperature").item(0).getChildNodes().item(0)).getNodeValue();
	   			}else{
	   				sTemperature= "0.0";
	   			}

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
		   			boolean single_no_mux= (!this.sensorPresent[chn][0] && !this.sensorPresent[chn][1] && !this.sensorPresent[chn][2]);
			   		if (this.debugLevel>1) System.out.println("Frame number: "+this.cameraFrameNumber[chn]+
			   				", Trigger mode:"+(this.triggeredMode[chn]?"ex":"in")+"ternal, "+
			   				" sensors attached:"+(this.sensorPresent[chn][0]?"1 ":"")+(this.sensorPresent[chn][1]?"2 ":"")+(this.sensorPresent[chn][2]?"3 ":"")+
			   				(single_no_mux?"single-sensor, no multiplexer":"") +
			   				(" Master port "+this.cameraMasterPort[chn])
			   				);
			   		for (int i=0;i<this.sensorPresent[chn].length;i++) if (this.sensorPresent[chn][i]) numSensors++;
			   		if (single_no_mux) numSensors++;
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
	   		this.triggeredMode=    new boolean[numChn];
	   		this.sensorPresent=    new boolean[numChn][];
			// TODO - try if skipping setting TRIG_PERIOD=0, IRQ_SMART=6 (when they are already set) will fix hanging
	   		this.triggerPeriod=    new int[numChn];
	   		this.irqSmart=         new int[numChn];
	   		this.cameraMasterPort= new int[numChn]; // not sure if will use it - master port for each camera:port
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
	   		int colon_index = this.cameraIPs[chn].indexOf(":");
	   		int sensor_port = Integer.parseInt(this.cameraIPs[chn].substring(colon_index+1)) - this.imgsrvPort;
	   		String ip=this.cameraIPs[chn].substring(0, colon_index);

	   		if (this.nc393){
//	   			url="http://"+this.cameraIPs[chn]+"/parsedit.php?immediate&TRIG&TRIG_PERIOD&SENS_AVAIL&FRAME";
	   			url="http://"+ip+"/parsedit.php?sensor_port="+sensor_port+"&immediate&TRIG&TRIG_PERIOD&SENS_AVAIL&FRAME&TRIG_MASTER";
	   		} else {
//	   			url="http://"+this.cameraIPs[chn]+"/parsedit.php?immediate&TRIG&TRIG_PERIOD&IRQ_SMART&SENS_AVAIL&FRAME";
	   			url="http://"+ip+"/parsedit.php?immediate&TRIG&TRIG_PERIOD&IRQ_SMART&SENS_AVAIL&FRAME";
	   		}
	   		if (this.debugLevel>1) System.out.println("url="+url);
	   		Document dom=null;
	   		if ((this.cameraFrameNumber==null) ||    ((this.cameraFrameNumber.length<(chn+1)))) initCameraArrays(chn+1);
	   		// did not find yet - why cameraFrameNumber is not null, but sensorPresent is? Added next lines until find
	   		if ((this.triggeredMode==null) ||        ((this.triggeredMode.length<(chn+1)))) initCameraArrays(chn+1);
	   		if ((this.sensorPresent==null) ||        ((this.sensorPresent.length<(chn+1)))) initCameraArrays(chn+1);
	   		if ((this.triggerPeriod==null) ||        ((this.triggerPeriod.length<(chn+1)))) initCameraArrays(chn+1);
	   		if (!this.nc393) {
	   			if ((this.irqSmart==null) ||         ((this.irqSmart.length<(chn+1))))         initCameraArrays(chn+1);
	   		} else {
	   			if ((this.cameraMasterPort==null) || ((this.cameraMasterPort.length<(chn+1)))) initCameraArrays(chn+1);
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

	   			this.triggeredMode[chn]= (Integer.parseInt((dom.getDocumentElement().getElementsByTagName("TRIG").item(0).getChildNodes().item(0)).getNodeValue())==4);
	   			int sensAvail= Integer.parseInt((dom.getDocumentElement().getElementsByTagName("SENS_AVAIL").item(0).getChildNodes().item(0)).getNodeValue());

	   			if (this.sensorPresent==null) System.out.println("sensorPresent==null");
	   			if (this.sensorPresent.length<(chn+1))  System.out.println("sensorPresent="+this.sensorPresent.length+" chn="+chn);

	   			this.sensorPresent[chn]=new boolean[3];
	   			for (int i=0;i<this.sensorPresent[chn].length;i++) this.sensorPresent[chn][i]=(sensAvail & (1<<i))!=0;
	   			this.triggerPeriod[chn]=     Integer.parseInt((dom.getDocumentElement().getElementsByTagName("TRIG_PERIOD").item(0).getChildNodes().item(0)).getNodeValue());
	   			if (!this.nc393) {
	   				this.irqSmart[chn]=          Integer.parseInt((dom.getDocumentElement().getElementsByTagName("IRQ_SMART").item(0).getChildNodes().item(0)).getNodeValue());
	   				this.cameraMasterPort[chn]=  sensor_port;
	   			} else {
	   				this.cameraMasterPort[chn]=  Integer.parseInt((dom.getDocumentElement().getElementsByTagName("TRIG_MASTER").item(0).getChildNodes().item(0)).getNodeValue());
	   			}
	   			this.cameraFrameNumber[chn]= Integer.parseInt((dom.getDocumentElement().getElementsByTagName("FRAME").item(0).getChildNodes().item(0)).getNodeValue());
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
	   		final boolean mp14 = true;// FIXME: Modified for 14Mpix
	   		final int ipLength=this.resetURLs.length;
	   		final boolean [] results=new boolean[ipLength];
	   		for (int chn=0;chn<ipLength;chn++) results[chn]=(this.sensorPresent[chn]!=null);
	   		final Thread[] threads = newThreadArray(this.maxNumberOfThreads);
	   		final AtomicInteger chnAtomic = new AtomicInteger(0);
	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				@Override
					public void run() {
	   					for (int chn=chnAtomic.getAndIncrement(); chn<ipLength;chn=chnAtomic.getAndIncrement()) if (results[chn]){
	   						results[chn]=setupCameraAcquisition(chn,exposureScale);
	   					}
	   				}

	   			};
	   		}
	   		startAndJoin(threads);
	   		if (Double.isNaN(exposureScale)){ // full init
//	   			int nRepeat=this.setupTriggerMode?4:1;
	   			int nRepeat=this.setupTriggerMode?4:2; // FIXME: Modified for 14Mpix
	   			if (this.nc393) nRepeat++; // is it needed?
	   			for (int i=0;i<nRepeat;i++){
	   				if (this.debugLevel>0) System.out.println((i+1)+" of "+nRepeat+": Triggering cameras to give parameters a chance to propagate");
	   				trigger();
	   				try
	   				{
	   					Thread.sleep( 2000 ); // ms // FIXME: Modified for 14Mpix
	   				}
	   				catch ( InterruptedException e )
	   				{
	   					System.out.println( "awakened prematurely" );
	   				}
	   			}

	   		} else if (mp14) {
   				try
   				{
   					Thread.sleep( 2000 ); // ms // FIXME: Modified for 14Mpix
   				}
   				catch ( InterruptedException e )
   				{
   					System.out.println( "awakened prematurely" );
   				}
   				trigger();
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
			  for ( long  endTime=System.nanoTime()+(1000000000);endTime<System.nanoTime(););
	   		}
	   		return result;
	   	}
	   	// TODO: issue several TRIG pulses after setting parameters?
	   	public boolean setupCameraAcquisition(int chn, double exposureScale){
	   		int colon_index = this.cameraIPs[chn].indexOf(":");
	   		int sensor_port = Integer.parseInt(this.cameraIPs[chn].substring(colon_index+1)) - this.imgsrvPort;
	   		String ip=this.cameraIPs[chn].substring(0, colon_index);
	   		boolean isMasterSubcamera = ip.equals(this.cameraSubnet+(this.iBaseIP+this.masterSubCamera));
	   		boolean isMasterPort =  (this.cameraMasterPort[chn] == sensor_port); // for 353 should always be master port
	   		if (isMasterSubcamera && (this.cameraMasterPort[chn] != this.masterPort)){
	   			System.out.println("Master port mismatch for camera "+ip+" (master) - parameters master port = "+
	   					this.masterPort+", camera responded with "+this.cameraMasterPort[chn] );
	   		}
	   		if (this.debugLevel>1) System.out.println("DEBUG393: chn="+chn+" sensor_port="+sensor_port+" ip= "+ip+" isMasterSubcamera="+isMasterSubcamera+ " isMasterPort="+isMasterPort);
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

	   		String triggerMode="";


	   		if (this.nc393) {
	   			if (isMasterPort) {
	   				/*
	   				 * P_TRIG_OUT (outputs):
	   				 * off:      0x00000
	   				 * external: 0x02000
	   				 * internal: 0x20000
	   				 * both:     0x22000
	   				 * P_TRIG_IN (inputs):
	   				 * off:      0x00000
	   				 * external: 0x80000
	   				 * internal: 0x08000
	   				 * both:     0x88000
	   				 */
	   				if (this.setupTriggerMode){
	   					triggerMode+="&TRIG_CONDITION="+(this.noCabling?(0):(this.externalTriggerCabling?"0x80000":"0x08000"))+"*0"+
	   							"&TRIG_OUT="+(this.noCabling?(0):(this.externalTriggerCabling?"0x02000":"0x20000"))+"*0"+
	   							"&TRIG=4*3";
	   				}

	   				if ((this.triggerPeriod[chn]>1) || this.setupTriggerMode) {
	   					triggerMode+="&TRIG_PERIOD=0*1"; // just stop it if it wasn't already, imgsrv /trig does not set it, only FPGA register
	   				}
	   			}
	   		}else {
	   			if (this.setupTriggerMode){
	   				triggerMode+="&TRIG_CONDITION="+(this.noCabling?(0):(this.externalTriggerCabling?"0x200000":"0x20000"))+"*0"+
	   						"&TRIG_OUT="+(this.noCabling?(0):(this.externalTriggerCabling?"0x800000":"0x80000"))+"*0"+
	   						"&TRIG=4*3";
	   			}
	   			if ((this.triggerPeriod[chn]>1) || this.setupTriggerMode){
	   				triggerMode+="&TRIG_PERIOD=1*0"; // just imgsrv /trig does not set it, only FPGA register
	   			}
	   		}


	   		String url="http://"+ip+"/parsedit.php?immediate";

	   		if (this.nc393) {
	   			url += "&sensor_port="+sensor_port;

	   		}



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
	   			if (!exposureOnly) this.cameraFrameNumber[chn]= Integer.parseInt((dom.getDocumentElement().getElementsByTagName("FRAME").item(0).getChildNodes().item(0)).getNodeValue());
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
			GenericDialog gd = new GenericDialog("Edit IPs/hosts and ports of the sub cameras");
			for (int i=0;i<this.cameraIPs.length;i++){
	    		gd.addStringField(i+": IP address/host:sensor port of the subcamera",this.cameraIPs[i],25);
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
//	   			IPs[this.channelMap[i][0]]=true;
	   			IPs[this.channelIPPort[i]] = true;  // since NC393
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
	   				@Override
					public void run() {
	   					for (int ipIndex=ipIndexAtomic.getAndIncrement(); ipIndex<ipLength;ipIndex=ipIndexAtomic.getAndIncrement()){

	   						//			for (int ipIndex=0;ipIndex<this.resetURLs.length;ipIndex++) if ((ipIndex<ipSelection.length) && ipSelection[ipIndex]){
	   						try {
	   							Document dom=null;
	   							DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
	   							DocumentBuilder db = dbf.newDocumentBuilder();
	   							String url=metaURLs[ipIndex];
//	   							if (debugLevel>2) System.out.println("timestampIPs:" + url );
	   							if (debugLevel>0) System.out.println("timestampIPs:" + url );
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
	   									(dom.getDocumentElement().getElementsByTagName("timestamp").item(0).getChildNodes().item(0)).getNodeValue());

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
	   				@Override
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

//	   	public boolean preTrigger()	{ // mt9f002 fails if more than 1-3 sec between triggers 
//	   		return true;
//	   	}
	   	
	   	
	   	public ImagePlus [] getImages(final boolean [] acquire, boolean resetAndTrigger, final boolean show){
		   	final boolean dual_trig = false;
		   	final boolean [] acquireIPs=selectIPs(acquire);
			if (this.debugLevel>2) {
				System.out.println("getImages(...) 2");
			    for (int ii=0;ii<this.cameraIPs.length;ii++)System.out.println(ii+":  "+this.cameraIPs[ii]);
			}

			if (resetAndTrigger) {
				try { // FIXME: Make conditional for 14MPix
					Thread.sleep(2000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				resetCameras();
				trigger();
				if (dual_trig) {
					timestampIPs(acquireIPs); // ignoring results, just wait
					resetCameras();
					//				  trigger();
					jp4_Instances[0].readDummyURL("http://192.168.0.236/parsedit.php?immediate&TRIG_PERIOD=1*0"); 
				}
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
					@Override
					public void run() {
						for (int ipIndex=ipIndexAtomic.getAndIncrement(); ipIndex<ipLength;ipIndex=ipIndexAtomic.getAndIncrement()){
							//							for (int i=0;i<this.imagesIP.length;i++){
							if (debugLevel>1/*2*/) System.out.println("getImages()3: ipIndex="+ipIndex+" acquireIPs.length=" +acquireIPs.length+
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
			final int [] channelIPPort = this.channelIPPort;
	   		final AtomicInteger imageIndexAtomic = new AtomicInteger(0);
	   		final int [] motorsPosition=this.motorsPosition;

	   		for (int ithread = 0; ithread < threads.length; ithread++) {
	   			threads[ithread] = new Thread() {
	   				@Override
					public void run() {
	   					// Next calls do not require image acquisition, so default settings for the new instance are OK
						JP46_Reader_camera jp4_Instance= new JP46_Reader_camera(false);
	   					for (int imageIndex=imageIndexAtomic.getAndIncrement(); imageIndex<images.length;imageIndex=imageIndexAtomic.getAndIncrement())
	   						if ((imageIndex<acquire.length) && acquire[imageIndex]){
	   							//		   	for (int i=0;i<this.images.length;i++) if ((i<acquire.length) && acquire[i]) {
//	   							int iIP=channelMap[imageIndex][0];
	   							int iIP=channelIPPort[imageIndex]; // index in composite images (per ip/port)
	   							if (sensorPresent[iIP]==null) { // system board for this channel did not respond null pointer - check cameras were detected
	   								images[imageIndex]=null;
	   								continue;
	   							}
	   							boolean singleSensor=true;
	   							for (int s=0;s<sensorPresent[iIP].length;s++) singleSensor&= !sensorPresent[iIP][s];
	   							if (singleSensor){ // no 10359 multiplexor
	   								if (debugLevel>1) System.out.println("DEBUG393: demuxing single sensor imageIndex="+imageIndex+" iIP="+iIP);
	   								images[imageIndex]=jp4_Instance.demuxClone(imagesIP[iIP]);
	   							} else {
	   								int subCam=channelMap[imageIndex][1];
	   								if (!sensorPresent[iIP][subCam]){ // requested sensor does not exist
	   									images[imageIndex]=null;
	   									continue;
	   								}
	   								// skip missing channels, then demux
	   								for (int s=0;s<channelMap[imageIndex][1];s++) if (!sensorPresent[iIP][s]) subCam--;
	   								if (debugLevel>1) System.out.println("DEBUG393: demuxing imageIndex="+imageIndex+" iIP="+iIP+" subCam="+subCam);
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
	   			try { // FIXME: Make conditional for 14MPix
					Thread.sleep(2000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	   			trigger(); // trigger cameras
	   			ipIndexAtomic.set(0);
	   			final int fnSeqNum=nSeqNum;
	   			final boolean reportTiming=this.reportTiming;
	   			for (int ithread = 0; ithread < threads.length; ithread++) {
	   				threads[ithread] = new Thread() {
	   					@Override
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
//	   							laserImagesIP[fnSeqNum-1][ipIndex].setTitle("raw"+fnSeqNum);  //FIXME:Debugging 
//	   							laserImagesIP[fnSeqNum-1][ipIndex].show(); //FIXME:Debugging
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
			final int []   channelIPPort = this.channelIPPort;
	   		final ImagePlus [] images=this.images;
			final boolean [] flipImages=this.flipImages;
			final LaserPointersHardware laserPointers=this.laserPointers;
			final boolean [][] sensorPresent=this.sensorPresent;
			final int debugSensorNumber=this.debugSensorNumber;
   			printTiming("=== Acquisition done, starting multi-threaded laser pointer location processing");

//TODO: Multithread the next cycle (per-sensor)
			final AtomicInteger sensorNumAtomic = new AtomicInteger(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
	   					// Next calls do not require image acquisition, so default settings for the new instance are OK
						JP46_Reader_camera jp4_Instance= new JP46_Reader_camera(false);
						MatchSimulatedPattern matchSimulatedPattern = new MatchSimulatedPattern();
						for (int sensorNum=sensorNumAtomic.getAndIncrement(); sensorNum<lasers.length;sensorNum=sensorNumAtomic.getAndIncrement()) // null pointer
							if (lasers[sensorNum] && (images[sensorNum]!=null)){
								//	   		for (int sensorNum=0;sensorNum<lasers.length;sensorNum++)  if (lasers[sensorNum] && (this.images[sensorNum]!=null)){ // lasers - here sensors to use lasers for
//								int iIP=channelMap[sensorNum][0];
	   							int iIP=channelIPPort[sensorNum]; // index in composite images (per ip/port)
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
											
											if ((laserPointers.laserPointer.laserUVMap!=null) && (laserPointers.laserPointer.laserUVMap[nPointer]!=null)) {
												images[sensorNum].setProperty("POINTER_U_"+nPointer, IJ.d2s(laserPointers.laserPointer.laserUVMap[nPointer][0],1));
												images[sensorNum].setProperty("POINTER_V_"+nPointer, IJ.d2s(laserPointers.laserPointer.laserUVMap[nPointer][1],1));
												
											}
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
//					int j=this.channelMap[i][0];
					int j=channelIPPort[i]; // sincer NC393 adder port
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