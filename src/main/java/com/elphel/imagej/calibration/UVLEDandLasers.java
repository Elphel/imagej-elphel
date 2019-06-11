package com.elphel.imagej.calibration;

import java.io.IOException;
import java.net.MalformedURLException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import ij.IJ;
import ij.gui.GenericDialog;

public class UVLEDandLasers{
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
    		if (result)	updateCurrents();
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
    				String quoted=((dom.getDocumentElement().getElementsByTagName("data9500").item(0).getChildNodes().item(0)).getNodeValue());
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
					return false;
//					throw new IllegalArgumentException (msg);
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