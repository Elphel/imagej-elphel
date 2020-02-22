package com.elphel.imagej.calibration.hardware;

import java.io.IOException;
import java.net.MalformedURLException;
import java.util.concurrent.atomic.AtomicInteger;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import ij.IJ;

public class GoniometerMotors{
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
			
			nanoETA += (long)(1E9*this.motorsStuckTestTime); // =5.0; // seconds
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
        		this.curpos[0]=Integer.parseInt((dom.getDocumentElement().getElementsByTagName("motor1").item(0).getChildNodes().item(0)).getNodeValue());
        		this.curpos[1]=Integer.parseInt((dom.getDocumentElement().getElementsByTagName("motor2").item(0).getChildNodes().item(0)).getNodeValue());
        		this.curpos[2]=Integer.parseInt((dom.getDocumentElement().getElementsByTagName("motor3").item(0).getChildNodes().item(0)).getNodeValue());
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