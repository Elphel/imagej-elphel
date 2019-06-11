package com.elphel.imagej.calibration.hardware;

import java.io.IOException;
import java.net.MalformedURLException;
import java.util.Properties;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import com.elphel.imagej.calibration.LaserPointer;
import com.elphel.imagej.common.WindowTools;

import ij.IJ;
import ij.gui.GenericDialog;

public class LaserPointersHardware{
    	public int debugLevel=1;
    	public LaserPointer laserPointer=null; // class that defines geometry and extracts laser coordinates
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
    	public LaserPointersHardware(LaserPointer laserPointer){
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