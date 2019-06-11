package com.elphel.imagej.calibration.hardware;

import java.io.IOException;
import java.io.StringWriter;
import java.net.MalformedURLException;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

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

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import com.elphel.imagej.common.WindowTools;

import ij.IJ;
import ij.gui.GenericDialog;

public class PowerControl{
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
					if (debugLevel>0) System.out.println("state="+(dom.getDocumentElement().getElementsByTagName("state").item(0).getChildNodes().item(0)).getNodeValue());
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