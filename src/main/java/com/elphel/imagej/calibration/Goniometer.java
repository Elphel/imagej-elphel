package com.elphel.imagej.calibration;
/**
 ** -----------------------------------------------------------------------------**
 ** Goniometer.java
 **
 ** Measurements in the "goniometer" machine
 **
 **
 ** Copyright (C) 2010 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Goniometer.java is free software: you can redistribute it and/or modify
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

import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.calibration.hardware.CamerasInterface;
import com.elphel.imagej.calibration.hardware.GoniometerMotors;
import com.elphel.imagej.cameras.EyesisCameraParameters;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.lwir.LwirReader;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;

public class Goniometer {
	/*
horizontal axis:
131 * 244 * 64 = 2045696
244 - worm gear
131 - motor
64 - pulses per revolution
5682.48889 per degree
	 */
	private ShowDoubleFloatArrays sdfaInstance = new ShowDoubleFloatArrays(); // just
																				// for
																				// debugging
	public CamerasInterface cameras = null;
	public LwirReader lwirReader = null;
	// public CalibrationHardwareInterface.LaserPointers lasers = null;
	// public static CalibrationHardwareInterface.FocusingMotors motorsS=null;
	// public DistortionProcessConfiguration
	// distortionProcessConfiguration=null;

//	public LensAdjustment.FocusMeasurementParameters focusMeasurementParameters = null;
	// public PatternParameters patternParameters=null;
	// public LensDistortionParameters
	// lensDistortionParameters=null;
//	public MatchSimulatedPattern.DistortionParameters distortion = null;
	public MatchSimulatedPattern.DistortionParameters distortionParametersDefault=null;
	public EyesisCameraParameters eyesisCameraParameters = null;

	public MatchSimulatedPattern[] matchSimulatedPatterns = null; // =new
																	// MatchSimulatedPattern();
	public LaserPointer laserPointers = null;
	MatchSimulatedPattern.PatternDetectParameters patternDetectParameters=null;
	public SimulationPattern.SimulParameters simulParametersDefault=null;
	public Goniometer.GoniometerParameters goniometerParameters = null;
    public DistortionProcessConfiguration distortionProcessConfiguration=null;
    public int lastScanStep=-1;
	public int debugLevel = 2;
	public double bottomRollerTilt=60.0; // decrease scan step if tilt is above this
	public double bottomRollersClearance=36.0; // angular clearance between the two bottom rollers


	public Goniometer(CamerasInterface cameras,
			MatchSimulatedPattern.DistortionParameters distortionParametersDefault,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			EyesisCameraParameters eyesisCameraParameters,
			LaserPointer laserPointers,
			SimulationPattern.SimulParameters simulParametersDefault,
			Goniometer.GoniometerParameters goniometerParameters,
			DistortionProcessConfiguration distortionProcessConfiguration
			) {
		this.cameras = cameras;
		this.distortionParametersDefault = distortionParametersDefault;
//		this.distortion = distortion;
		this.patternDetectParameters=patternDetectParameters;
		this.eyesisCameraParameters = eyesisCameraParameters;
		this.laserPointers = laserPointers;
		this.simulParametersDefault=simulParametersDefault;
		this.goniometerParameters=goniometerParameters;
		this.distortionProcessConfiguration=distortionProcessConfiguration;
	}

	public Goniometer(
			LwirReader lwirReader,
//			CalibrationHardwareInterface.CamerasInterface cameras,
			MatchSimulatedPattern.DistortionParameters distortionParametersDefault,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			EyesisCameraParameters eyesisCameraParameters,
			LaserPointer laserPointers,
			SimulationPattern.SimulParameters simulParametersDefault,
			Goniometer.GoniometerParameters goniometerParameters,
			DistortionProcessConfiguration distortionProcessConfiguration
			) {
		this.lwirReader= lwirReader;
//		this.cameras = cameras;
		this.distortionParametersDefault = distortionParametersDefault;
//		this.distortion = distortion;
		this.patternDetectParameters=patternDetectParameters;
		this.eyesisCameraParameters = eyesisCameraParameters;
		this.laserPointers = laserPointers;
		this.simulParametersDefault=simulParametersDefault;
		this.goniometerParameters=goniometerParameters;
		this.distortionProcessConfiguration=distortionProcessConfiguration;
	}




	//goniometerMotors
	private enum MOT_ACT {
		MOVE_SPEC,
		HOME,
		TILT_P85,
		TILT_M85,
		RCW200,
		RCCW200,
		SET_HOME};
	private final String [] options={
			"Move to specified position",
			"Move home",
			"Tilt up (from target) 85 degrees",
			"Tilt down (to target) 85 degrees",
			"Rotate to clockwise 200 degrees",
			"Rotate to counter-clockwise 200 degrees",
			"set current position as new home"};


	public boolean manualMove(
			AtomicInteger stopRequested, // 1 - stop now, 2 - when convenient
			boolean updateStatus){
		int tiltMotor=2;  // 0-1-2
		int axialMotor=1; // 0-1-2
		boolean needsInit=false;

		int [] currentMotors=this.goniometerParameters.goniometerMotors.getCurrentPositions();
		double currentTilt=currentMotors[tiltMotor]/this.goniometerParameters.goniometerMotors.stepsPerDegreeTilt;
		double currentAxial=currentMotors[axialMotor]/this.goniometerParameters.goniometerMotors.stepsPerDegreeAxial;
		currentTilt=0.1*Math.round(10.0*currentTilt);
		currentAxial=0.1*Math.round(10.0*currentAxial);
		this.goniometerParameters.goniometerMotors.debugLevel=this.debugLevel+1;
		System.out.println(
				"Current position:\n"+
				"Tilt:  "+IJ.d2s(currentTilt,1)+" degrees ("+currentMotors[tiltMotor]+" steps)\n"+
				"Axial: "+IJ.d2s(currentAxial,1)+" degrees ("+currentMotors[axialMotor]+" steps)\n");

		GenericDialog gd = new GenericDialog("User interrupt");

		gd.addRadioButtonGroup("Select action", options, 7, 1, options[0]);
		gd.addNumericField("Goniometer tilt angle",    currentTilt, 1,5,"\u00b0 (positive - look up)");
		gd.addNumericField("Goniometer axial angle",    currentAxial, 1,5,"\u00b0 (positive - CCW)");
		gd.addCheckbox("Initialize goniometer motor driver",needsInit);
		gd.addNumericField("Debug level",    debugLevel, 0,1,"");
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		String selectedAction=gd.getNextRadioButton();
		currentTilt=  gd.getNextNumber();
		currentAxial= gd.getNextNumber();
		needsInit=    gd.getNextBoolean();
		debugLevel=(int)gd.getNextNumber();
		this.goniometerParameters.goniometerMotors.debugLevel=this.debugLevel; // +1;
		if (needsInit){
			if (!this.goniometerParameters.goniometerMotors.tryInit(true,updateStatus)){
				String msg="Failed to initialize goniometer motor driver";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false; // failed initialization
			}
		}
		if (selectedAction.equals(options[MOT_ACT.SET_HOME.ordinal()])){
			if (this.goniometerParameters.goniometerMotors.setHome()==null){
				String msg="Failed to set new home position";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false; // failed set home
			}
			return true;
		}
		MOT_ACT act=MOT_ACT.MOVE_SPEC;
		for (MOT_ACT a: MOT_ACT.values()){
			if (selectedAction.equals(options[a.ordinal()])){
				act=a;
				break;
			}
		}
		double targetTilt=currentTilt;
		double targetAxial=currentAxial;
		switch (act) {
		case HOME:
			targetTilt=0;
			targetAxial=0;
			act=MOT_ACT.MOVE_SPEC;
			break;
		case TILT_P85:
			targetTilt=85.0;
			act=MOT_ACT.MOVE_SPEC;
			break;
		case TILT_M85:
			targetTilt=-85.0;
			act=MOT_ACT.MOVE_SPEC;
			break;
		case RCW200:
			targetAxial=-200.0;
			act=MOT_ACT.MOVE_SPEC;
			break;
		case RCCW200:
			targetAxial=200.0;
			act=MOT_ACT.MOVE_SPEC;
			break;
		default:
			break;
		}
		if (act!=MOT_ACT.MOVE_SPEC){
			String msg="Program BUG, action "+selectedAction+" can not be processed";
			System.out.println("Error: "+msg);
			IJ.showMessage("Error",msg);
			return false;
		}
		int tiltMotorPosition= (int) Math.round(targetTilt*this.goniometerParameters.goniometerMotors.stepsPerDegreeTilt);
		int axialMotorPosition= (int) Math.round(targetAxial*this.goniometerParameters.goniometerMotors.stepsPerDegreeAxial);
		boolean OK= motorsMove(
				tiltMotor,
				axialMotor,
				tiltMotorPosition,
				axialMotorPosition,
				stopRequested,
				updateStatus);
		if (!OK) System.out.println("motorsMove()->false");
		return true; // OK; // So will re-open dialog even after abort
	}


	public boolean scanAndAcquire(
			double targetAngleHorizontal,
			double targetAngleVertical,
		    AtomicInteger stopRequested, // 1 - stop now, 2 - when convenient
			boolean updateStatus
			){

		int tiltMotor=2;  // 0-1-2
		int axialMotor=1; // 0-1-2
		int [] motors=this.goniometerParameters.goniometerMotors.updateMotorsPosition();
		double thisTilt= motors[tiltMotor]/this.goniometerParameters.goniometerMotors.stepsPerDegreeTilt;
		double thisAxial=motors[axialMotor]/this.goniometerParameters.goniometerMotors.stepsPerDegreeAxial;
		double scanOverlapVertical=this.goniometerParameters.scanOverlapVertical;
		double scanOverlapHorizontal=this.goniometerParameters.scanOverlapHorizontal;

		boolean  reverseAxial=(this.goniometerParameters.scanLimitAxialStart>this.goniometerParameters.scanLimitAxialEnd);
        double scanLimitAxialLow= reverseAxial?this.goniometerParameters.scanLimitAxialEnd:this.goniometerParameters.scanLimitAxialStart;
        double scanLimitAxialHigh=   reverseAxial?this.goniometerParameters.scanLimitAxialStart:this.goniometerParameters.scanLimitAxialEnd;


		boolean zenithToNadir=this.goniometerParameters.scanLatitudeHigh<this.goniometerParameters.scanLatitudeLow;
		double scanLatitudeHigh=zenithToNadir? this.goniometerParameters.scanLatitudeLow: this.goniometerParameters.scanLatitudeHigh;
		double scanLatitudeLow= zenithToNadir? this.goniometerParameters.scanLatitudeHigh:this.goniometerParameters.scanLatitudeLow;

		double scanStepTilt= targetAngleVertical*  (1.0-scanOverlapVertical);
		double scanStepAxial=targetAngleHorizontal*(1.0-scanOverlapHorizontal); // valid at equator
		if (this.debugLevel>1) System.out.println("scanStepTilt="+IJ.d2s(scanStepTilt,2)+", scanStepAxial="+IJ.d2s(scanStepAxial,2));

		int numTiltSteps=(int) Math.ceil((scanLatitudeHigh-scanLatitudeLow)/scanStepTilt); // includes first and last

		if (numTiltSteps>0){ // increase vertical overlap to make it same for all images
		   scanStepTilt=(scanLatitudeHigh-scanLatitudeLow)/numTiltSteps;
		   scanOverlapVertical=1.0-(scanStepTilt/targetAngleVertical);
		}
		if (this.debugLevel>1) System.out.println("Updated scanStepTilt="+IJ.d2s(scanStepTilt,2)+", scanOverlapVertical="+IJ.d2s(scanOverlapVertical,2));


		double [] tilts=new double [numTiltSteps];
		double [][] rots= new double [numTiltSteps][];
		boolean dirAxial=reverseAxial;
		int numStops=0;
		for (int i=0;i<numTiltSteps;i++){
			int tiltIndex= zenithToNadir?(numTiltSteps-i -1):i;
			double tilt=-(scanLatitudeLow+tiltIndex*scanStepTilt+ 0.5*(1.0-scanOverlapVertical)*targetAngleVertical);
			tilts [i]=tilt;
			double tiltL=tilt-0.5*targetAngleVertical;
			double tiltH=tilt+0.5*targetAngleVertical;
			if (tiltL<-90.0) tiltL=-180.0-tiltL;
			if (tiltH> 90.0) tiltH= 180.0-tiltH;
			double minAbsTilt=Math.min(Math.abs(tiltL),Math.abs(tiltH));
			if ((tiltL*tiltH)<0) minAbsTilt=0.0;
			double cosMinAbsTilt=Math.cos(minAbsTilt*Math.PI/180.0);
			if (this.debugLevel>2) System.out.println("tilt="+IJ.d2s(tilt,2)+"tiltL="+IJ.d2s(tiltL,2)+", tiltH="+IJ.d2s(tiltH,2)+
					", minAbsTilt="+IJ.d2s(minAbsTilt,2)+", cosMinAbsTilt="+IJ.d2s(cosMinAbsTilt,2));
			double scanStep=scanStepAxial;
			double overlap=scanOverlapHorizontal;
			int numAxialSteps=(int) Math.ceil((scanLimitAxialHigh-scanLimitAxialLow)*cosMinAbsTilt/scanStepAxial);

// Correction for bottom rollers that block view - if tilt is above 60degrees (positive - looking higher). The clear angle is ~36%
			if (tilts[i]>this.bottomRollerTilt){
				int numForRollers=(int) Math.ceil((scanLimitAxialHigh-scanLimitAxialLow)/(this.bottomRollersClearance*(1.0-scanOverlapHorizontal)));
				if (numForRollers>numAxialSteps){
					if (this.debugLevel>0){
						System.out.println("Increasing number of steps to mitigate occlusion by the bottom rollers. Original number of steps: "+
								numAxialSteps+", increased: "+numForRollers);
					}
					numAxialSteps=numForRollers;
				}
			}


			// spread evenly
			if (numAxialSteps>0){ // increase vertical overlap to make it same for all images
				scanStep=(scanLimitAxialHigh-scanLimitAxialLow)/numAxialSteps;
				overlap=1.0-(scanStep*cosMinAbsTilt/targetAngleHorizontal);
			}
			rots[i]=new double[numAxialSteps];
			for (int j=0;j<numAxialSteps;j++){
				int axialIndex= dirAxial?(numAxialSteps-j-1):j;
				double axial=(scanLimitAxialLow+axialIndex*scanStep+ 0.5*(1.0-overlap)*targetAngleHorizontal/cosMinAbsTilt);
				rots[i][j]=axial;
			}
			if (this.goniometerParameters.scanBidirectional) dirAxial=!dirAxial;
			numStops+=numAxialSteps;
		}
		if (this.debugLevel>0){
			System.out.println("First tilt: "+IJ.d2s(tilts[0],1)+" last tilt: "+IJ.d2s(tilts[numTiltSteps-1],1)+
					", number of tilt steps: "+numTiltSteps+", total number of measurements - "+numStops);
		}
		if (this.debugLevel>1){
			System.out.println("targetAngleHorizontal="+IJ.d2s(targetAngleHorizontal,2)+", targetAngleVertical="+IJ.d2s(targetAngleVertical,2));
		  for (int i=0;i<numTiltSteps;i++) {
			  System.out.println("Tilt # "+i+": "+tilts[i]+"\n axial ("+rots[i].length+" stops):");
			  for (int j=0;j<rots[i].length;j++) System.out.print(" "+IJ.d2s(rots[i][j],1));
			  System.out.println();
		  }
		}

		GenericDialog gd = new GenericDialog("Start scanning");
		gd.addMessage("About to start scanning and recording images, parameters are set in the \"Configure Goniometer\" dialog");
		gd.addMessage("Please make sure goniometer motors are set that 0,0 corresponds to the camera in the initial position -");
		gd.addMessage("vertical and cables are not entangled, camera exposure is set correctly.");
		gd.addMessage("Camera will start from the tilt "+tilts[0]+" degrees (0 is vertical, positive - away from the target) and move to "+
				tilts[tilts.length-1]+" degrees,");
		gd.addMessage("making "+numStops+" stops for image acquisirtion from all channels.");
		gd.addMessage("");
		gd.addMessage("Axial rotations are set to "+(this.goniometerParameters.scanBidirectional?"ALTERNATE directions.":"be all in THE SAME direction."));
		gd.addMessage("Simultaneous operation of the motors is "+(this.goniometerParameters.motorsSimultaneous?"ENABLED.":"DISABLED."));
		int startStep=0;
		if (this.lastScanStep>=0){
			gd.addMessage("Last scan finished at stop "+this.lastScanStep);
			if (this.lastScanStep<(numStops-1)) startStep= this.lastScanStep+1;

		}
		gd.addMessage("");
		gd.addMessage("Current position is:");
		gd.addMessage("tilt="+ IJ.d2s(thisTilt,2)+ " degrees ("+motors[tiltMotor]+ " motor steps)");
		gd.addMessage("axial="+IJ.d2s(thisAxial,2)+" degrees ("+motors[axialMotor]+" motor steps)");
		gd.addNumericField("Start from position (0.."+(numStops-1)+")", startStep, 0);
		gd.addCheckbox("Debug timing", false);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		startStep= (int) gd.getNextNumber();
		boolean debugTiming=gd.getNextBoolean();
		String src_dir=this.distortionProcessConfiguration.selectSourceDirectory(true, this.distortionProcessConfiguration.sourceDirectory, true);
		if (src_dir==null) {
			String msg="Failed to select directory "+this.distortionProcessConfiguration.sourceDirectory+" to save images ";
			System.out.println("Error: "+msg);
			IJ.showMessage("Error",msg);
			return false;
		}
		// overwrite some distortionProcessConfiguration parameters (show, save, ...?) from GoniometerParameters
		this.distortionProcessConfiguration.sourceDirectory=src_dir;
// just for now - setting motor debug 1 higher than this
		this.goniometerParameters.goniometerMotors.debugLevel=this.debugLevel+1;
		long startTime = System.nanoTime();
		String status;
		boolean OK;
		int startTilt=0;
		int startAxial=0;
		this.lastScanStep=0;
		while (startStep>=(this.lastScanStep+rots[startTilt].length)){
			this.lastScanStep+=rots[startTilt].length;
			startTilt++;
		}
		startAxial=startStep-this.lastScanStep;
		this.lastScanStep=startStep-1;
		 Runtime runtime = Runtime.getRuntime();
		for (int nTilt=startTilt;nTilt<numTiltSteps;nTilt++){
			double tilt=tilts[nTilt];
			tilt=0.1*Math.round(10*tilt); // is that needed?
			int tiltMotorPosition= (int) Math.round(tilt*this.goniometerParameters.goniometerMotors.stepsPerDegreeTilt);

			status=IJ.d2s(1E-9*(System.nanoTime()-startTime),3)+": Tilt run "+(nTilt+1)+" (of "+numTiltSteps+"), tilt angle "+
					IJ.d2s(tilt,2)+" degrees, motor steps: "+tiltMotorPosition;
			if (this.debugLevel>0) System.out.println(status);
			if (updateStatus) IJ.showStatus(status);
			this.goniometerParameters.motorsSimultaneous=false; // not yet implemented
			OK= this.goniometerParameters.goniometerMotors.moveMotorSetETA(tiltMotor, tiltMotorPosition);
			if (!OK) {
				String msg="Could not set motor "+(tiltMotor+1)+" to move to "+tiltMotorPosition+" - may be out of limit";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false;
			}
			OK=this.goniometerParameters.goniometerMotors.waitMotor(tiltMotor, null, false,updateStatus); // no interrupts during movement
			if (!OK) {
				String msg="Motor "+(tiltMotor+1)+" failed to reach "+tiltMotorPosition+".";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false;
			}
			//			}
			for (int nAxial=((nTilt==startTilt)?startAxial:0);nAxial<rots[nTilt].length;nAxial++){
				double axial=rots[nTilt][nAxial];
				axial=0.1*Math.round(10*axial);
				int axialMotorPosition= (int) Math.round(axial*this.goniometerParameters.goniometerMotors.stepsPerDegreeAxial);
				status=(this.lastScanStep+1)+"( last="+(numStops-1)+") "+IJ.d2s(1E-9*(System.nanoTime()-startTime),3)+" sec.  tilt:"+(nTilt+1)+"/"+numTiltSteps+" , "+
						IJ.d2s(tilt,2)+" deg., axial:"+
						(nAxial+1)+"/"+rots[nTilt].length+" , "+IJ.d2s(axial,2)+" deg.";
				runtime.gc();
				String memoryStatus="Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")";
				if (this.debugLevel>0) System.out.println(status+", axial motor steps: "+axialMotorPosition+" "+memoryStatus);
				if (updateStatus) IJ.showStatus(status);

				OK= this.goniometerParameters.goniometerMotors.moveMotorSetETA(axialMotor, axialMotorPosition);
				if (!OK) {
					String msg="Could not set motor "+(axialMotor+1)+" to move to "+axialMotorPosition+" - may be out of limit";
					System.out.println("Error: "+msg);
					IJ.showMessage("Error",msg);
					return false;
				}
				OK=this.goniometerParameters.goniometerMotors.waitMotor(axialMotor, null, false,updateStatus);
				if (!OK) {
					String msg="Motor "+(axialMotor+1)+" failed to reach "+axialMotorPosition+".";
					System.out.println("Error: "+msg);
					IJ.showMessage("Error",msg);
					return false;
				}
				// update motor positions in the image properties, acquire and save images.
				// TODO: Make acquisition/decoding/laser identification multi-threaded

				if (cameras != null) {
					this.cameras.setMotorsPosition(this.goniometerParameters.goniometerMotors.getTargetPositions()); // Used target, not current to prevent minor variations
					this.cameras.reportTiming=debugTiming;
					this.cameras.acquire(this.distortionProcessConfiguration.sourceDirectory,true, updateStatus); // true - use lasers, updateStatus - make false?
				} else if (lwirReader != null) {
					this.lwirReader.setMotorsPosition(this.goniometerParameters.goniometerMotors.getTargetPositions()); // Used target, not current to prevent minor variations
					this.lwirReader.reportTiming=debugTiming;
					this.lwirReader.acquire(this.distortionProcessConfiguration.sourceDirectory); // true - use lasers, updateStatus - make false?
				} else {
					System.out.println("Neignter traditional camera/rig, no LWIR rig are initialized, dry run");
				}

				this.lastScanStep++;
				if (stopRequested.get()>1){
					if (this.debugLevel>0) System.out.println("User interrupt");
					stopRequested.set(0);
					gd = new GenericDialog("User interrupt");
					//	    		this.serialNumber, // camera serial number string
					gd.addMessage("User requested interrupt after step "+this.lastScanStep+" (last would be "+(numStops-1)+")");
					System.out.println("User requested interrupt after step "+this.lastScanStep+" (last would be "+(numStops-1)+")");
					gd.addMessage("Cancel will terminate the scanning leaving motros where they are.");
					gd.enableYesNoCancel("Continue", "Motors Home");
					gd.showDialog();
					if (gd.wasCanceled()) return false;
					if (!gd.wasOKed()){
						motorsMove(
								tiltMotor,
								axialMotor,
								0, //tiltMotorPosition,
								0, //axialMotorPosition,
								null, // no interrupts
								updateStatus);
						return false;
					}
				}

			}

		}
		OK=motorsMove(
				tiltMotor,
				axialMotor,
				0, //tiltMotorPosition,
				0, //axialMotorPosition,
				null, // no interrupts
				updateStatus);
		if (this.debugLevel>0) System.out.println("Scan finished in "+IJ.d2s(1E-9*(System.nanoTime()-startTime),3)+" seconds.");

		return OK;

	}
	public boolean motorsMove(
			int tiltMotor,
			int axialMotor,
			int tiltMotorPosition,
			int axialMotorPosition,
			AtomicInteger stopRequested, // or null
			boolean updateStatus){
		String status;
		if (!this.goniometerParameters.goniometerMotors.checkGotTarget(axialMotor,axialMotorPosition)) {
			status="Moving axial motor to "+axialMotorPosition+"...";
			if (updateStatus) IJ.showStatus(status);
			boolean OK= this.goniometerParameters.goniometerMotors.moveMotorSetETA(axialMotor, axialMotorPosition);
			if (!OK) {
				String msg="Could not set motor "+(axialMotor+1)+" to move to "+axialMotorPosition+" - may be out of limit";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false;
			}
			OK=this.goniometerParameters.goniometerMotors.waitMotor(axialMotor, stopRequested, false,updateStatus);
			if (!OK) {
				String msg="Motor "+(axialMotor+1)+" failed to reach "+axialMotorPosition+".";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false;
			}
		}
		if (!this.goniometerParameters.goniometerMotors.checkGotTarget(tiltMotor,tiltMotorPosition)) {
			status="Moving tilt motor to "+tiltMotorPosition+"...";
			if (updateStatus) IJ.showStatus(status);
			boolean OK= this.goniometerParameters.goniometerMotors.moveMotorSetETA(tiltMotor, tiltMotorPosition);
			if (!OK) {
				String msg="Could not set motor "+(tiltMotor+1)+" to move to "+tiltMotorPosition+" - may be out of limit";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false;
			}
			OK=this.goniometerParameters.goniometerMotors.waitMotor(tiltMotor,stopRequested, false, updateStatus);
			if (!OK) {
				String msg="Motor "+(tiltMotor+1)+" failed to reach "+tiltMotorPosition+".";
				System.out.println("Error: "+msg);
				IJ.showMessage("Error",msg);
				return false;
			}
		}
		return true;
	}


	public boolean testHintedTarget (
			ImagePlus[] images,
			Distortions lensDistortions, // should not be null
			DistortionCalibrationData distortionCalibrationData,
			PatternParameters patternParameters, // should not be  null
			boolean equalizeGreens,
			int threadsMax,
			boolean updateStatus,
			int debug_level) {// debug level used inside loops
		if (lensDistortions == null) {
			String msg = "lensDistortions is not initialized";
			IJ.showMessage("Error", msg);
			throw new IllegalArgumentException(msg);
		}
		GenericDialog gd = new GenericDialog("Select image to try hint pattern");
		int numNonNull=0;
		for (int nImg=0;nImg<images.length;nImg++) if (images[nImg]!=null) numNonNull++;
		String [] choices=new String[numNonNull];
		int [] pointersArray=new int[images.length];
		int [] indices=new int [numNonNull];
		int index=0;
//		if (this.debugLevel>1)System.out.println ("images.length="+images.length);
		for (int nImg=0;nImg<images.length;nImg++) if (images[nImg]!=null){
			int subCam=distortionCalibrationData.getImageChannel(images[nImg]);
//			String imgTitle=images[nImg].getTitle();
			pointersArray[nImg]=0;
			if (images[nImg].getProperty("POINTERS")!=null) pointersArray[nImg]=Integer.parseInt((String) (images[nImg].getProperty("POINTERS")));

			choices[index]=index+" sub-camera:"+subCam+" "+images[nImg].getTitle()+" - "+pointersArray[nImg]+" pointers";
//			if (this.debugLevel>1)System.out.println ("Adding image, nImg="+nImg);
//			if (this.debugLevel>1)System.out.println ("Adding image "+images[nImg].getTitle()+": "+choices[index]);
			indices[index++]=nImg;
		}
		gd.addChoice("Image to process", choices,  choices[0]);
		double hintGridTolerance=0.0;
		gd.addNumericField("grid match tolerance (0 - only direction)", hintGridTolerance, 3);
		int stationNumber=0;
		if (distortionCalibrationData.getNumStations()>1) {
			gd.addNumericField("Station number (0.."+(distortionCalibrationData.getNumStations()-1)+")", stationNumber, 0);
		}

		gd.showDialog();
		if (gd.wasCanceled()) return false;
		index= gd.getNextChoiceIndex();
		hintGridTolerance=gd.getNextNumber();
		if (distortionCalibrationData.getNumStations()>1) {
			stationNumber=(int) gd.getNextNumber();
			if (stationNumber<0) stationNumber=0;
			else if (stationNumber>=distortionCalibrationData.getNumStations()) stationNumber = distortionCalibrationData.getNumStations()-1;
		}

		int nImg=indices[index];

		int subCam=        distortionCalibrationData.getImageChannel(images[nImg]);
//		int stationNumber= distortionCalibrationData.getImageStation(numGridImage), // station number
		double timeStamp=  distortionCalibrationData.getImageTimestamp(images[nImg]);

		int numPointers=pointersArray[nImg];
		double [] goniometerTiltAxial=distortionCalibrationData.getImagesetTiltAxial(timeStamp);
		if (goniometerTiltAxial==null){
			String msg= "No orientation data for timestamp="+IJ.d2s(timeStamp,6);
			System.out.println(msg);
			IJ.showMessage("Warning",msg);
			return false;
		}
		if (Double.isNaN(goniometerTiltAxial[0])){
			String msg= "Goniometer tilt angle is undefined for timestamp="+IJ.d2s(timeStamp,6);
			System.out.println(msg);
			IJ.showMessage("Warning",msg);
		}
		if (Double.isNaN(goniometerTiltAxial[1])){
			String msg= "Goniometer axial angle is undefined for timestamp="+IJ.d2s(timeStamp,6);
			System.out.println(msg);
			IJ.showMessage("Warning",msg);
		}
		if ((Double.isNaN(goniometerTiltAxial[0])) || (Double.isNaN(goniometerTiltAxial[1]))) return false;

		double [][][] hintGrid=lensDistortions.estimateGridOnSensor(
				stationNumber, // station number
				subCam,
				goniometerTiltAxial[0], // Tilt, goniometerHorizontal
				goniometerTiltAxial[1],  // Axial,goniometerAxial
				goniometerTiltAxial[2],  // inter-axis angle
				-1, // use camera parameters, not imageSet
				true // filter border
				);
		if (hintGrid==null){
			String msg= "Target is not visible";
			System.out.println(msg);
			IJ.showMessage("Warning",msg);
			return false;
		}
		if (this.debugLevel>1) lensDistortions.showHintGrid(hintGrid);


		MatchSimulatedPattern matchSimulatedPattern = new MatchSimulatedPattern(this.distortionParametersDefault.FFTSize); // new instance, all reset
		// next 2 lines are not needed for the new instance, but can be
		// used alternatively if keeping it
		matchSimulatedPattern.invalidateFlatFieldForGrid(); // Reset Flat Filed calibration - different image.
		matchSimulatedPattern.invalidateFocusMask();
		matchSimulatedPattern.debugLevel = debug_level;
		ImagePlus imp_eq = matchSimulatedPattern.applyFlatField(images[nImg]); // current image with grid flat-field  correction

//TODO: it shows always 4 pointers (if>0) - the array is sparse
//	   	public int getNumberOfPointers (int sensorNum){
//	   	public int getNumberOfPointers (ImagePlus imp){

		if (debug_level > 0){
			System.out.println("\n   ======= Looking for grid, matching pointers in image " +images[nImg].getTitle()+
					", initial number of pointers was "+numPointers);
		}
//matchSimulatedPatterns[numSensor].getChannel(images[numSensor])+" ");
		MatchSimulatedPattern.DistortionParameters distortionParameters = modifyDistortionParameters();
		SimulationPattern.SimulParameters simulParameters = modifySimulParameters();

		boolean noMessageBoxes=true;
		int numAbsolutePoints = matchSimulatedPattern.calculateDistortions(
				        null, // LwirReaderParameters lwirReaderParameters, // null is OK
						// allow more of grid around pointers?
						distortionParameters, //
						this.patternDetectParameters,
//						this.patternDetectParameters.minGridPeriod/2,
//						this.patternDetectParameters.maxGridPeriod/2,
						simulParameters,
						equalizeGreens, imp_eq,
						this.laserPointers, // null, //LASER_POINTERS, //
										// LaserPointer laserPointer, //
										// LaserPointer object or null
						true, // don't care -removeOutOfGridPointers
						hintGrid, //   double [][][] hintGrid, // predicted grid array (or null)
						hintGridTolerance,    //   double  hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
						threadsMax,
						updateStatus,
						debug_level,
						distortionParameters.loop_debug_level, // debug level
						noMessageBoxes);
		if (numAbsolutePoints < 0) { // no pointers in this image
			String msg = "*** No laser pointers matched for " + images[nImg].getTitle() + " - they are needed for absolute grid positioning";
			if (debug_level > 0) System.out.println("Warning: " + msg);
			if (debug_level > 2) IJ.showMessage("Warning", msg);
		}
		if (numAbsolutePoints == 0) { // no pointers in this image
			String msg = "*** No laser pointers matched for " + images[nImg].getTitle() + ", but the grid fit hinted one with specified tolerance";
			if (debug_level > 0) System.out.println("Warning: " + msg);
			if (debug_level > 2) IJ.showMessage("Warning", msg);
		}

		return true;

	}

	/*
	 * 	private showDoubleFloatArrays SDFA_INSTANCE= new showDoubleFloatArrays(); // just for debugging?
   		this.SDFA_INSTANCE.showArrays(gridXYZCorr, getGridWidth(), getGridHeight(),  true, "Grid corrections", titles);

	 *
					gd.addChoice( // ArrayIndexOutOfBoundsException: 21
							this.distortionCalibrationData.getParameterName(parIndex)+
							" ("+sValue+" "+
							this.distortionCalibrationData.getParameterUnits(parIndex)+")"+
							(this.distortionCalibrationData.isSubcameraParameter(parIndex)?(" s"+subCam):"com "),
							this.definedModes, this.definedModes[this.parameterMode[numSeries][i]]);

	 *
	 * 					this.parameterMode[numSeries][i]=gd.getNextChoiceIndex();
			PatternParameters patternParameters, // should not be  null
			boolean equalizeGreens,
			int threadsMax,
			boolean updateStatus,
			int debug_level) {// debug level used inside loops

	 *
	 *
	 */
	public MatchSimulatedPattern.DistortionParameters modifyDistortionParameters(){
		MatchSimulatedPattern.DistortionParameters distortionParameters = this.distortionParametersDefault.clone();
		distortionParameters.refineInPlace = false;
		distortionParameters.correlationMaxOffset = this.goniometerParameters.maxCorr;
		distortionParameters.correlationSize = this.goniometerParameters.correlationSize;
		distortionParameters.correlationGaussWidth = this.goniometerParameters.correlationGaussWidth;
		distortionParameters.refineCorrelations = false;
		distortionParameters.fastCorrelationOnFirstPass = true;
		distortionParameters.fastCorrelationOnFinalPass = true;
		distortionParameters.correlationAverageOnRefine = false;
		distortionParameters.minUVSpan = this.goniometerParameters.minUVSpan;
		distortionParameters.flatFieldCorrection = this.goniometerParameters.flatFieldCorrection;
		distortionParameters.flatFieldExpand = this.goniometerParameters.flatFieldExpand;
		distortionParameters.numberExtrapolated = 1; // measuring distortions -
		distortionParameters.correlationMinInitialContrast=this.goniometerParameters.correlationMinInitialContrast;
		distortionParameters.minimalPatternCluster=this.goniometerParameters.minimalPatternCluster;
		distortionParameters.scaleMinimalInitialContrast=this.goniometerParameters.scaleMinimalInitialContrast;
		distortionParameters.searchOverlap=this.goniometerParameters.searchOverlap;
		return distortionParameters;
	}

	public SimulationPattern.SimulParameters modifySimulParameters(){
		SimulationPattern.SimulParameters simulParameters = this.simulParametersDefault.clone();
		simulParameters.smallestSubPix = this.goniometerParameters.smallestSubPix;
		simulParameters.bitmapNonuniforityThreshold = this.goniometerParameters.bitmapNonuniforityThreshold;
		simulParameters.subdiv = this.goniometerParameters.subdiv;
		return simulParameters;
	}

	public double[] estimateOrientation(
			ImagePlus[] images, // last acquire images with number of pointers
								// detected>0
			DistortionCalibrationData distortionCalibrationData,
			PatternParameters patternParameters, // should not be  null
			Distortions lensDistortions, // should not be null
			LaserPointer laserPointers, // as a back up measure
			boolean equalizeGreens,
			int threadsMax,
			boolean updateStatus,
			int debug_level) {// debug level used inside loops
		long startTime = System.nanoTime();
		if (lensDistortions == null) {
			String msg = "lensDistortions is not initialized";
			IJ.showMessage("Error", msg);
			throw new IllegalArgumentException(msg);
		}

		// remove unneeded, copied from updateFocusGrid()
		SimulationPattern.SimulParameters simulParameters = modifySimulParameters();
		MatchSimulatedPattern.DistortionParameters distortionParameters = modifyDistortionParameters();

		int numImages = 0;
		for (int i = 0; i < images.length; i++)
			if (images[i] != null)
				numImages++;
		if (numImages == 0) {
			String msg = "No images with laser pointers";
			System.out.println("Error: " + msg);
			IJ.showMessage("Error", msg);
			return null;
		}
		if (this.matchSimulatedPatterns == null) {
			this.matchSimulatedPatterns = new MatchSimulatedPattern[images.length];
			for (int i = 0; i < this.matchSimulatedPatterns.length; i++)
				this.matchSimulatedPatterns[i] = null;
		}
		if (this.matchSimulatedPatterns.length < images.length) {
			MatchSimulatedPattern[] matchSimulatedPatternsTmp = matchSimulatedPatterns.clone();
			this.matchSimulatedPatterns = new MatchSimulatedPattern[images.length];
			for (int i = 0; i < matchSimulatedPatternsTmp.length; i++)
				this.matchSimulatedPatterns[i] = matchSimulatedPatternsTmp[i];
		}
		ImagePlus[] imp_calibrated = new ImagePlus[images.length];
		int [] numPointers=new int [images.length];
		for (int numSensor = 0; numSensor < imp_calibrated.length; numSensor++) {
			imp_calibrated[numSensor] = null;
			numPointers[numSensor] = 0;
		}
		boolean noMessageBoxes=true;
		for (int numSensor = 0; numSensor < images.length; numSensor++)
			if (images[numSensor] != null) {
				// reset matchSimulatedPattern, so it will start from scratch
				this.matchSimulatedPatterns[numSensor] = new MatchSimulatedPattern(
						this.distortionParametersDefault.FFTSize); // new instance, all reset
				// next 2 lines are not needed for the new instance, but can be
				// used alternatively if keeping it
				this.matchSimulatedPatterns[numSensor].invalidateFlatFieldForGrid(); // Reset Flat Filed calibration - different image.
				this.matchSimulatedPatterns[numSensor].invalidateFocusMask();
				if (MatchSimulatedPattern.getPointersXYUV(images[numSensor],this.laserPointers) == null) { // no pointers in this image
					String msg = "No laser pointers detected for "+ images[numSensor].getTitle()+ " - they are needed for absolute grid positioning";
					if (debug_level > 0) System.out.println("Warning: " + msg);
					IJ.showMessage("Warning", msg);
					continue;
				} else if (debug_level > 0) {
//					System.out.println("Image "+numSensor+" has "+ matchSimulatedPatterns[numSensor].getPointersXY(images[numSensor],this.laserPointers.laserUVMap.length).length+ " pointers");
				}
				this.matchSimulatedPatterns[numSensor].debugLevel = debug_level;
				ImagePlus imp_eq = this.matchSimulatedPatterns[numSensor].applyFlatField(images[numSensor]); // current image with grid flat-field  correction

//TODO: it shows always 4 pointers (if>0) - the array is sparse
//			   	public int getNumberOfPointers (int sensorNum){
//			   	public int getNumberOfPointers (ImagePlus imp){
				numPointers[numSensor]=0;
		   		if (images[numSensor].getProperty("POINTERS")!=null) numPointers[numSensor]= Integer.parseInt((String) images[numSensor].getProperty("POINTERS"));


				if (debug_level > 0){
					System.out.println("\n   ======= Looking for grid, matching pointers in image " +images[numSensor].getTitle()+
							", initial number of pointers was "+numPointers[numSensor]);
				}
//	matchSimulatedPatterns[numSensor].getChannel(images[numSensor])+" ");

				int numAbsolutePoints = this.matchSimulatedPatterns[numSensor].calculateDistortions(
						null, // LwirReaderParameters lwirReaderParameters, // null is OK
						// allow more of grid around pointers?
						distortionParameters, //
						this.patternDetectParameters,
//						this.patternDetectParameters.minGridPeriod/2,
//						this.patternDetectParameters.maxGridPeriod/2,
						simulParameters,
						equalizeGreens, imp_eq,
						this.laserPointers, // null, //LASER_POINTERS, //
						// LaserPointer laserPointer, //
						// LaserPointer object or null
						true, // don't care -removeOutOfGridPointers
						null, //   double [][][] hintGrid, // predicted grid array (or null)
						0,    //   double  hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only

						threadsMax,
						updateStatus,
						debug_level,
						distortionParameters.loop_debug_level, // debug level
						noMessageBoxes);
				if (numAbsolutePoints <= 0) { // no pointers in this image
					String msg = "*** No laser pointers matched for " + images[numSensor].getTitle() + " - they are needed for absolute grid positioning";
					if (debug_level > 0) System.out.println("Warning: " + msg);
					if (debug_level > 2) IJ.showMessage("Warning", msg);
					continue;
				}
				numPointers[numSensor]=numAbsolutePoints;

				// if (debug_level>1)
				// System.out.println("calculateDistortions() for channel "+numSensor+" finished at "+
				// IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				if (debug_level > 2) {
					double[] test_uv = new double[this.matchSimulatedPatterns[numSensor].UV_INDEX.length];
					for (int i = 0; i < this.matchSimulatedPatterns[numSensor].UV_INDEX.length; i++)
						test_uv[i] = this.matchSimulatedPatterns[numSensor].UV_INDEX[i];
					sdfaInstance.showArrays(test_uv,
							this.matchSimulatedPatterns[numSensor].getImageWidth(),
							this.matchSimulatedPatterns[numSensor].getImageHeight(), "UV_INDEX");
				}
				if (debug_level > 0)
					System.out.println("Matched "
									+ numAbsolutePoints
									+ " laser pointers, grid generated at "
									+ IJ.d2s(0.000000001 * (System.nanoTime() - startTime), 3));
				// TODO - here - multiple images possible, not just one!
				// First - create sparse, then remove nulls

				imp_calibrated[numSensor] = this.matchSimulatedPatterns[numSensor].getCalibratedPatternAsImage(images[numSensor],"grid-",numAbsolutePoints);
				if (this.goniometerParameters.showAcquiredImages)
					imp_calibrated[numSensor].show(); // DISTORTION_PROCESS_CONFIGURATION.showGridImages
			} // for (int numSensor=0;numSensor<images.length;numSensor++) if
				// (images[numSensor]!=null) {
		int numCalibrated = 0;
		for (int i = 0; i < imp_calibrated.length; i++)
			if (imp_calibrated[i] != null)
				numCalibrated++;
		if (numCalibrated == 0) {
			String msg = "No absolutely calibrated (with laser pointers) images available";
			if (debug_level > 0)
				System.out.println("Error: " + msg);
			IJ.showMessage("Error", msg);
			return null;
		} else {
			String msg = "Found " + numCalibrated + " images with the grid referenced by laser pointers";
			if (debug_level > 1) System.out.println(msg);
			if (debug_level > 2) IJ.showMessage("Info", msg);
		}
		int maxNumPointers=0;
		int imgWithMaxPointers=0;
		for (int numSensor=0;numSensor<numPointers.length;numSensor++){
			if (numPointers[numSensor]>maxNumPointers){
				maxNumPointers=numPointers[numSensor];
				imgWithMaxPointers=numSensor;
			}
		}
		if (numCalibrated < imp_calibrated.length) {
			ImagePlus[] imp_tmp = imp_calibrated.clone();
			imp_calibrated = new ImagePlus[numCalibrated];
			numCalibrated = 0;
			for (int i = 0; i < imp_tmp.length; i++) if (imp_tmp[i] != null){
				if (i==imgWithMaxPointers) imgWithMaxPointers= numCalibrated; // may only be decreased or stay the same, so
				// imgWithMaxPointes is a new index of the image with maximal number of recognized laser pointers
				imp_calibrated[numCalibrated++] = imp_tmp[i];
			}
		}
		/*
		 * DistortionCalibrationData distortionCalibrationData= new
		 * DistortionCalibrationData( imp_calibrated, //ImagePlus []
		 * images, // images in the memory patternParameters,
		 * //PatternParameters patternParameters, eyesisCameraParameters
		 * //EyesisCameraParameters eyesisCameraParameters );
		 */
		distortionCalibrationData.setImages(imp_calibrated, // ImagePlus [] images, // images in the memory
				patternParameters, // PatternParameters patternParameters);
				laserPointers);
		distortionCalibrationData.initImageSet(eyesisCameraParameters);

		// Set initial azimuth and elevation
		// Set initial heading and elevation
		double [] initialHeadEl=distortionCalibrationData.getHeadEl(imgWithMaxPointers);
        // set goniometer horizontal axis angle and goniometer axial angles in all images
		distortionCalibrationData.setGHGA(-initialHeadEl[1], -initialHeadEl[0]);
		if (debug_level > 1) System.out.println("Initial Heading and Elevation are set to heading="+IJ.d2s(-initialHeadEl[0],2)+", elvation="+IJ.d2s(-initialHeadEl[1],2));

		lensDistortions.copySensorConstants(eyesisCameraParameters); // copy from the first channel
		// lensDistortions.fittingStrategy will be defined later, no need to
		// update it with a reference to distortionCalibrationData now
//		if (debug_level > 1) {
//			System.out.println("distortionCalibrationData.setImages()");
//		}

		// fitting strategy
//		distortionCalibrationData.pathName=this.goniometerParameters.initialCalibrationFile;
		lensDistortions.debugLevel = this.debugLevel;
		lensDistortions.fittingStrategy = new FittingStrategy(true,
				this.goniometerParameters.strategyFile,
				distortionCalibrationData); // will use list of grid files
		if (lensDistortions.fittingStrategy.pathName == null) { // failed to select/open the file
			lensDistortions.fittingStrategy = null;
			IJ.showMessage("Error", "Failed to open fitting strategy file: "+this.goniometerParameters.strategyFile);
			return null;
		}
		if (debug_level > 1) System.out.println("Using fitting strategy template file: "+ lensDistortions.fittingStrategy.pathName);

		// TODO: modify fitting strategy to include all grid images
		lensDistortions.fittingStrategy.adjustNumberOfImages(imp_calibrated.length);
		this.goniometerParameters.strategyFile = lensDistortions.fittingStrategy.pathName;
		// saved fitting strategy maybe for different number of subcameras.
		// it will be adjusted here - that works only for simple strategies that use all subcameras at each step.


		// TODO: fix repeating subcamera parameters for all subcameras in each strategy
		lensDistortions.fittingStrategy.updateNumberOfSubcameras();
// enable azimuth adjust for all but the first camera? not here, later

		// Calculate Sensor Masks
		distortionCalibrationData.debugLevel = debug_level;
		distortionCalibrationData.updateStatus = updateStatus;
		distortionCalibrationData.calculateSensorMasks();

		if (debug_level > 0) System.out.println("Starting LMA at " + IJ.d2s(0.000000001 * (System.nanoTime() - startTime), 3));

		lensDistortions.seriesNumber = 0; // start from 0;
		lensDistortions.stopEachStep = false;
		lensDistortions.stopEachSeries = false;

		lensDistortions.thresholdFinish = this.goniometerParameters.thresholdFinish;
		lensDistortions.numIterations = this.goniometerParameters.numIterations;

		// TODO: Set initial values for the goniometer angles from the sensor
		// (channel) number, average them if there are several in the list

		lensDistortions.LevenbergMarquardt(
				false, // skip dialog
				false); // new: dry_run use it here?
		if (debug_level > 0)
			System.out.println("Finished LMA at "
					+ IJ.d2s(0.000000001 * (System.nanoTime() - startTime), 3));

		// Read camera parameters
		if (debug_level > 0)
			System.out.println("estimateOrientation finished at "
					+ IJ.d2s(0.000000001 * (System.nanoTime() - startTime), 3));
		// TODO: see if needs to be changed
        int stationNumber=0;
		double[] result = { this.eyesisCameraParameters.goniometerHorizontal[stationNumber], // goniometer rotation around "horizontal" axis (tilting from the target - positive)
				this.eyesisCameraParameters.goniometerAxial[stationNumber] // goniometer rotation around Eyesis axis (clockwise in plan - positive
		};
		return result;
	}
    public static class GoniometerParameters {
    	public String gridGeometryFile="";
    	public String initialCalibrationFile="";
    	public String strategyFile="";
    	public String resultsSuperDirectory=""; // directory with subdirectories named as serial numbers to store results
    	public String comment="no comments"; // Comment to add to the results

    	public int EEPROM_channel=1; // EEPROM channel to read serial number from
    	public boolean saveResults=true; // save focusing results
    	public boolean showResults=true; // show focusing (includingh intermediate) results
    	public String serialNumber=""; // camera serial number string
    	public double sensorTemperature=Double.NaN; // last measured sensor temperature
//other summary results to be saved with parameters
    	public double maxCorr=5;     // maximal grid correction between images allowed (larger will trigger full grid rebuild)
    	public boolean showHistoryDetails=false;   // show color info
    	public boolean showHistorySingleLine=true; // all parameters in a single line (easier to copy to spreadsheet)
        public boolean showAcquiredImages=false;
        public boolean showFittedParameters=true;
        // when approximating PSF with a second degree polynomial:
        public double psf_cutoffEnergy=0.5; // disregard pixels outside of this fraction of the total energy
        public double psf_cutoffLevel= 0.2; // disregard pixels below this fraction of the maximal value
        public int    psf_minArea    = 10;  // continue increasing the selected area, even if beyound psf_cutoffEnergy and psf_cutoffLevel,
        public double psf_blurSigma  = 0.0; // optionally blur the calculated mask

    	// the following  overwrite SimulParameters members
//TODO: Make initial pattern search more robust - if it first gets false positive, and number of detected cells is too low - increase threshold and re-try
		public double correlationMinInitialContrast=3.0;   // minimal contrast for the pattern of the center (initial point)
		public int    minimalPatternCluster=150;          //    minimal pattern cluster size (0 - disable retries)
		public double scaleMinimalInitialContrast=2.0;   // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
		public double searchOverlap=0.25;         // when searching for grid, step this amount of the FFTSize
    	public double smallestSubPix=0.3; // subdivide pixels down to that fraction when simulating
    	public double bitmapNonuniforityThreshold=0.1	; // subdivide pixels until difference between the corners is below this value
    	public int    subdiv=4;
    	// overwrites  	public static class MultiFilePSF.overexposedMaxFraction
    	public double overexposedMaxFraction=0.1; // allowed fraction of the overexposed pixels in the PSF kernel measurement area
    	// overwrites	public static class PSFParameters.minDefinedArea
    	public double minDefinedArea=0.75;    // minimal (weighted) fraction of the defined patter pixels in the FFT area
        public int PSFKernelSize=32;          // size of the detected PSF kernel
		public boolean approximateGrid=true; // approximate grid with polynomial
		public boolean centerPSF=true;       // Center PSF by modifying phase

		public double mask1_sigma=    1.0;
		public double mask1_threshold=0.25;
		public double gaps_sigma=     1.0;
		public double mask_denoise=   0.25;
    //  OTFFilterParameters
		public double deconvInvert=0.03; // with good focus can go to 0.015 or smaller
    //  DistortionParameters
        public int correlationSize=32;
        public double correlationGaussWidth=0.75;
        public double minUVSpan;           // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
        public boolean flatFieldCorrection=true;
        public double flatFieldExpand=4.0;

        public double thresholdFinish=0.001; // Stop iterations if 2 last steps had less improvement (but not worsening )

        public int    numIterations=  100; // maximal number of iterations

    	/*
        horizontal axis:
        131 * 244 * 64 = 2045696
        244 - worm gear
        131 - motor
        64 - pulses per revolution
        5682.48889 per degree
        	 */
        // motors rotate positive - look down, positive - CCW
        GoniometerMotors goniometerMotors=null;
//        public double stepsPerDegreeTilt=-5682.48889; // minus that positive steps make negative elevation
//        public double stepsPerDegreeAxial=-36.0; // minus that positive steps make rotate CCW when looking from Eyesis top
//        public double scanStepTilt=20.0;   // degrees (equal steps not larger than
//        public double scanStepAxial=10.0;        // degrees (equal steps not larger than
//        public double scanLimitTiltStart= 30.0;  // scan around horizontal axis from  that angle
//        public double scanLimitTiltEnd= -80.0;  // scan around horizontal axis to  that angle
        public double targetDistance=  5817; //mm - foer overlap calculation

        public double scanLatitudeLow=-90;   // lowest camera latitude to calibrate (nadir=-90)
        public double scanLatitudeHigh=90;  // highest camera latitude to calibrate (zenith=90)

        public double scanOverlapHorizontal= 0.5;
        public double scanOverlapVertical=   0.5;


        public double scanLimitAxialStart=     -200.0;  // scan around camera axis from  that angle
        public double scanLimitAxialEnd=     200.0;  // scan around camera axis  to  that angle
        public boolean scanBidirectional=     true;  // false - always move axial in the same direction, true - optimize movements
        public boolean motorsSimultaneous=     true;  // true - move motors simultaneously, false - one at a time



    	public GoniometerParameters(GoniometerMotors goniometerMotors){
    		    		this.goniometerMotors=goniometerMotors;
    	}
    	public GoniometerParameters(
    			GoniometerMotors goniometerMotors,
    	    	String gridGeometryFile,
    	    	String initialCalibrationFile, // not needed
    	    	String strategyFile,
    	    	String resultsSuperDirectory, // directory with subdirectories named as serial numbers to stro results
    	    	int EEPROM_channel, // EEPROM channel to read serial number from
    	    	boolean saveResults, // save focusing results
    	    	boolean showResults, // show focusing (includingh intermediate) results
    	    	String serialNumber, // camera serial number string
    	    	double sensorTemperature, // last measured sensor temperature
    	    	String comment, // Comment to add to the results
    			double maxCorr,     // maximal grid correction between images allowed (larger will trigger full grid rebuild)
    	    	boolean showHistoryDetails,
    	    	boolean showHistorySingleLine, // all parameters in a single line (easier to copy to spreadsheet)
    	    	boolean showAcquiredImages,
    	    	boolean showFittedParameters,
                double psf_cutoffEnergy, // disregard pixels outside of this fraction of the total energy
                double psf_cutoffLevel,  // disregard pixels below this fraction of the maximal value
                int    psf_minArea,      // continue increasing the selected area, even if beyound psf_cutoffEnergy and psf_cutoffLevel,
                                         // if the selected area is smaller than this (so approximation wpuld work)
                double psf_blurSigma,    // optionally blur the calculated mask
                double correlationMinInitialContrast,   // minimal contrast for the pattern of the center (initial point)
        		int    minimalPatternCluster,       //    minimal pattern cluster size (0 - disable retries)
        		double scaleMinimalInitialContrast, // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
        		double searchOverlap,         // when searching for grid, step this amount of the FFTSize
    	    	double smallestSubPix, // subdivide pixels down to that fraction when simulating
    	    	double bitmapNonuniforityThreshold, // subdivide pixels until difference between the corners is below this value
    	    	int    subdiv,
    	    	double overexposedMaxFraction, // allowed fraction of the overexposed pixels in the PSF kernel measurement area
    	    	double minDefinedArea, // minimal (weighted) fraction of the defined patter pixels in the FFT area
    	    	int PSFKernelSize,
				boolean approximateGrid, // approximate grid with polynomial
				boolean centerPSF,       // Center PSF by modifying phase
				double mask1_sigma,
				double mask1_threshold,
				double gaps_sigma,
				double mask_denoise,
				double deconvInvert,
                int correlationSize,
                double correlationGaussWidth,
                double minUVSpan,           // Minimal u/v span in correlation window that triggers increase of the correlation FFT size
                boolean flatFieldCorrection,
                double flatFieldExpand,
                double thresholdFinish,// (copied from series) stop iterations if 2 last steps had less improvement (but not worsening )
                int    numIterations, // maximal number of iterations
//                double stepsPerDegreeTilt, // minus that positive steps make negative elevation
//    	        double stepsPerDegreeAxial,      // minus that positive steps make rotate CCW when looking from Eyesis top
                double targetDistance,
                double scanLatitudeLow,   // lowest camera latitude to calibrate
                double scanLatitudeHigh,  // highest camera latitude to calibrate
                double scanOverlapHorizontal,
                double scanOverlapVertical,
                double scanLimitAxialStart,        // scan around camera axis from  that angle
    	        double scanLimitAxialEnd,       // scan around camera axis  to  that angle
    	        boolean scanBidirectional,
    	        boolean motorsSimultaneous
    			){
    		this.gridGeometryFile=gridGeometryFile;
    		this.initialCalibrationFile=initialCalibrationFile;
    		this.strategyFile=strategyFile;
    		this.resultsSuperDirectory=resultsSuperDirectory; // directory with subdirectories named as serial numbers to stro results
    		this.EEPROM_channel=EEPROM_channel; // EEPROM channel to read serial number from
    		this.saveResults=saveResults; // save focusing results
    		this.showResults=showResults; // show focusing (includingh intermediate) results
    		this.serialNumber=serialNumber; // camera serial number string
    		this.sensorTemperature=sensorTemperature; // last measured sensor temperature
    		this.comment=comment; // Comment to add to the results
			this.maxCorr=maxCorr;
			this.showHistoryDetails=showHistoryDetails;
			this.showHistorySingleLine=showHistorySingleLine; // all parameters in a single line (easier to copy to spreadsheet)
			this.showAcquiredImages=showAcquiredImages;
			this.showFittedParameters=showFittedParameters;
			this.psf_cutoffEnergy=psf_cutoffEnergy;
			this.psf_cutoffLevel= psf_cutoffLevel;
			this.psf_minArea=     psf_minArea;
			this.psf_blurSigma=   psf_blurSigma;
			this.correlationMinInitialContrast=correlationMinInitialContrast;   // minimal contrast for the pattern of the center (initial point)
			this.minimalPatternCluster=minimalPatternCluster;           //    minimal pattern cluster size (0 - disable retries)
			this.scaleMinimalInitialContrast=scaleMinimalInitialContrast; // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
			this.searchOverlap=searchOverlap;         // when searching for grid, step this amount of the FFTSize
			this.smallestSubPix=smallestSubPix;
			this.bitmapNonuniforityThreshold=bitmapNonuniforityThreshold;
			this.subdiv=subdiv;
			this.overexposedMaxFraction=overexposedMaxFraction;
			this.minDefinedArea=minDefinedArea;
			this.PSFKernelSize=PSFKernelSize;
			this.approximateGrid = approximateGrid; // approximate grid with polynomial
			this.centerPSF = centerPSF; // approximate grid with polynomial
			this.mask1_sigma = mask1_sigma;
			this.mask1_threshold = mask1_threshold;
			this.gaps_sigma=gaps_sigma;
			this.mask_denoise=mask_denoise;
			this.deconvInvert = deconvInvert;
			this.correlationSize=correlationSize;
			this.correlationGaussWidth=correlationGaussWidth;
			this.minUVSpan=minUVSpan;
			this.flatFieldCorrection=flatFieldCorrection;
			this.flatFieldExpand=flatFieldExpand;
			this.thresholdFinish=thresholdFinish;// (copied from series) stop iterations if 2 last steps had less improvement (but not worsening )
			this.numIterations=numIterations; // maximal number of iterations
			this.goniometerMotors=goniometerMotors;
//			this.goniometerMotors.stepsPerDegreeTilt=stepsPerDegreeTilt; // minus that positive steps make negative elevation
//			this.goniometerMotors.stepsPerDegreeAxial=stepsPerDegreeAxial;      // minus that positive steps make rotate CCW when looking from Eyesis top
			this.targetDistance=targetDistance;
			this.scanLatitudeLow=scanLatitudeLow;   // lowest camera latitude to calibrate
			this.scanLatitudeHigh=scanLatitudeHigh;  // highest camera latitude to calibrate
			this.scanOverlapHorizontal=scanOverlapHorizontal;
			this.scanOverlapVertical=scanOverlapVertical;
			this.scanLimitAxialStart=scanLimitAxialStart;        // scan around camera axis from  that angle
			this.scanLimitAxialEnd=scanLimitAxialEnd;       // scan around camera axis  to  that angle
			this.scanBidirectional=scanBidirectional;
			this.motorsSimultaneous=motorsSimultaneous;
    	}
    	@Override
		public GoniometerParameters clone(){
    		return new GoniometerParameters(
    				this.goniometerMotors,
    	    		this.gridGeometryFile,
    	    		this.initialCalibrationFile,
    	    		this.strategyFile,
    	    		this.resultsSuperDirectory, // directory with subdirectories named as serial numbers to stro results
    	    		this.EEPROM_channel,// EEPROM channel to read serial number from
    	    		this.saveResults, // save focusing results
    	    		this.showResults, // show focusing (includingh intermediate) results
    	    		this.serialNumber, // camera serial number string
    	    		this.sensorTemperature, // last measured sensor temperature
    	    		this.comment,
        			this.maxCorr,
    				this.showHistoryDetails,
    				this.showHistorySingleLine, // all parameters in a single line (easier to copy to spreadsheet)
    				this.showAcquiredImages,
    				this.showFittedParameters,
    				this.psf_cutoffEnergy,
    				this.psf_cutoffLevel,
    				this.psf_minArea,
    				this.psf_blurSigma,
    				this.correlationMinInitialContrast,   // minimal contrast for the pattern of the center (initial point)
    				this.minimalPatternCluster,           //    minimal pattern cluster size (0 - disable retries)
    				this.scaleMinimalInitialContrast, // increase/decrease minimal contrast if initial cluster is >0 but less than minimalPatternCluster
    				this.searchOverlap,         // when searching for grid, step this amount of the FFTSize
    				this.smallestSubPix,
    				this.bitmapNonuniforityThreshold,
    				this.subdiv,
    				this.overexposedMaxFraction,
    				this.minDefinedArea,
    				this.PSFKernelSize,
    				this.approximateGrid,
        			this.centerPSF,
        			this.mask1_sigma,
        			this.mask1_threshold,
        			this.gaps_sigma,
        			this.mask_denoise,
        			this.deconvInvert,
    				this.correlationSize,
    				this.correlationGaussWidth,
    				this.minUVSpan,
    				this.flatFieldCorrection,
    				this.flatFieldExpand,
    				this.thresholdFinish,
    				this.numIterations,
    				this.targetDistance,
    				this.scanLatitudeLow,
    				this.scanLatitudeHigh,
    				this.scanOverlapHorizontal,
    				this.scanOverlapVertical,
    				this.scanLimitAxialStart,        // scan around camera axis from  that angle
    				this.scanLimitAxialEnd,       // scan around camera axis  to  that angle
    				this.scanBidirectional,
    				this.motorsSimultaneous
        			);
    	}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"gridGeometryFile",this.gridGeometryFile+"");
			properties.setProperty(prefix+"initialCalibrationFile",this.initialCalibrationFile+"");
			properties.setProperty(prefix+"strategyFile",this.strategyFile+"");
			properties.setProperty(prefix+"resultsSuperDirectory",this.resultsSuperDirectory+"");
			properties.setProperty(prefix+"serialNumber",this.serialNumber);
			if (!Double.isNaN(this.sensorTemperature))properties.setProperty(prefix+"sensorTemperature",this.sensorTemperature+"");
			properties.setProperty(prefix+"EEPROM_channel",this.EEPROM_channel+"");
			properties.setProperty(prefix+"saveResults",this.saveResults+"");
			properties.setProperty(prefix+"showResults",this.showResults+"");
			properties.setProperty(prefix+"comment","<![CDATA["+this.comment+ "]]>");
			properties.setProperty(prefix+"maxCorr",this.maxCorr+"");
			properties.setProperty(prefix+"showHistoryDetails",this.showHistoryDetails+"");
			properties.setProperty(prefix+"showHistorySingleLine",this.showHistorySingleLine+"");
			properties.setProperty(prefix+"showAcquiredImages",this.showAcquiredImages+"");
			properties.setProperty(prefix+"showFittedParameters",this.showFittedParameters+"");
			properties.setProperty(prefix+"psf_cutoffEnergy",this.psf_cutoffEnergy+"");
			properties.setProperty(prefix+"psf_cutoffLevel",this.psf_cutoffLevel+"");
			properties.setProperty(prefix+"psf_minArea",this.psf_minArea+"");
			properties.setProperty(prefix+"psf_blurSigma",this.psf_blurSigma+"");
			properties.setProperty(prefix+"correlationMinInitialContrast",this.correlationMinInitialContrast+"");
			properties.setProperty(prefix+"minimalPatternCluster",this.minimalPatternCluster+"");
			properties.setProperty(prefix+"scaleMinimalInitialContrast",this.scaleMinimalInitialContrast+"");
			properties.setProperty(prefix+"searchOverlap",this.searchOverlap+"");
			properties.setProperty(prefix+"smallestSubPix",this.smallestSubPix+"");
			properties.setProperty(prefix+"bitmapNonuniforityThreshold",this.bitmapNonuniforityThreshold+"");
			properties.setProperty(prefix+"subdiv",this.subdiv+"");
			properties.setProperty(prefix+"overexposedMaxFraction",this.overexposedMaxFraction+"");
			properties.setProperty(prefix+"minDefinedArea",this.minDefinedArea+"");
			properties.setProperty(prefix+"PSFKernelSize",this.PSFKernelSize+"");
			properties.setProperty(prefix+"approximateGrid",this.approximateGrid+"");
			properties.setProperty(prefix+"centerPSF",this.centerPSF+"");
			properties.setProperty(prefix+"mask1_sigma",this.mask1_sigma+"");
			properties.setProperty(prefix+"mask1_threshold",this.mask1_threshold+"");
			properties.setProperty(prefix+"gaps_sigma",this.gaps_sigma+"");
			properties.setProperty(prefix+"mask_denoise",this.mask_denoise+"");
			properties.setProperty(prefix+"deconvInvert",this.deconvInvert+"");
			properties.setProperty(prefix+"correlationSize",this.correlationSize+"");
			properties.setProperty(prefix+"correlationGaussWidth",this.correlationGaussWidth+"");
			properties.setProperty(prefix+"minUVSpan",this.minUVSpan+"");
			properties.setProperty(prefix+"flatFieldCorrection",this.flatFieldCorrection+"");
			properties.setProperty(prefix+"flatFieldExpand",this.flatFieldExpand+"");
			properties.setProperty(prefix+"thresholdFinish",this.thresholdFinish+"");
			properties.setProperty(prefix+"numIterations",this.numIterations+"");
			properties.setProperty(prefix+"goniometerMotors_ipAddress",this.goniometerMotors.ipAddress+"");
			properties.setProperty(prefix+"goniometerMotors_stepsPerSecond",     this.goniometerMotors.stepsPerSecond+"");
			properties.setProperty(prefix+"goniometerMotors_stepsPerDegreeTilt", this.goniometerMotors.stepsPerDegreeTilt+"");
			properties.setProperty(prefix+"goniometerMotors_stepsPerDegreeAxial",this.goniometerMotors.stepsPerDegreeAxial+"");

			properties.setProperty(prefix+"targetDistance",this.targetDistance+"");
			properties.setProperty(prefix+"scanLatitudeLow",this.scanLatitudeLow+"");
			properties.setProperty(prefix+"scanLatitudeHigh",this.scanLatitudeHigh+"");
			properties.setProperty(prefix+"scanOverlapHorizontal",this.scanOverlapHorizontal+"");
			properties.setProperty(prefix+"scanOverlapVertical",this.scanOverlapVertical+"");

			properties.setProperty(prefix+"scanLimitAxialStart",this.scanLimitAxialStart+"");
			properties.setProperty(prefix+"scanLimitAxialEnd",this.scanLimitAxialEnd+"");
			properties.setProperty(prefix+"scanBidirectional",this.scanBidirectional+"");
			properties.setProperty(prefix+"motorsSimultaneous",this.motorsSimultaneous+"");


		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"gridGeometryFile")!=null)
				this.gridGeometryFile=properties.getProperty(prefix+"gridGeometryFile");
			if (properties.getProperty(prefix+"initialCalibrationFile")!=null)
				this.initialCalibrationFile=properties.getProperty(prefix+"initialCalibrationFile");
			if (properties.getProperty(prefix+"strategyFile")!=null)
				this.strategyFile=properties.getProperty(prefix+"strategyFile");
			if (properties.getProperty(prefix+"resultsSuperDirectory")!=null)
				this.resultsSuperDirectory=properties.getProperty(prefix+"resultsSuperDirectory");
			if (properties.getProperty(prefix+"serialNumber")!=null)
				this.serialNumber=properties.getProperty(prefix+"serialNumber");
			//	this.serialNumber is only written, but never read from the configuration file (only from device)

			if (properties.getProperty(prefix+"sensorTemperature")!=null) this.sensorTemperature=Double.parseDouble(properties.getProperty(prefix+"sensorTemperature"));
			else this.sensorTemperature=Double.NaN;
			if (properties.getProperty(prefix+"EEPROM_channel")!=null)
				this.EEPROM_channel=Integer.parseInt(properties.getProperty(prefix+"EEPROM_channel"));
			if (properties.getProperty(prefix+"saveResults")!=null)
				this.saveResults=Boolean.parseBoolean(properties.getProperty(prefix+"saveResults"));
			if (properties.getProperty(prefix+"showResults")!=null)
				this.showResults=Boolean.parseBoolean(properties.getProperty(prefix+"showResults"));
			if (properties.getProperty(prefix+"comment")!=null)
				this.comment=properties.getProperty(prefix+"comment");
			if ((this.comment.length()>10) && this.comment.substring(0,9).equals("<![CDATA[")) this.comment=this.comment.substring(9,this.comment.length()-3);
			if (properties.getProperty(prefix+"maxCorr")!=null)
				this.maxCorr=Double.parseDouble(properties.getProperty(prefix+"maxCorr"));
			if (properties.getProperty(prefix+"showHistoryDetails")!=null)
				this.showHistoryDetails=Boolean.parseBoolean(properties.getProperty(prefix+"showHistoryDetails"));
			if (properties.getProperty(prefix+"showHistorySingleLine")!=null)
				this.showHistorySingleLine=Boolean.parseBoolean(properties.getProperty(prefix+"showHistorySingleLine"));
			if (properties.getProperty(prefix+"showAcquiredImages")!=null)
				this.showAcquiredImages=Boolean.parseBoolean(properties.getProperty(prefix+"showAcquiredImages"));
			if (properties.getProperty(prefix+"showFittedParameters")!=null)
				this.showFittedParameters=Boolean.parseBoolean(properties.getProperty(prefix+"showFittedParameters"));
			if (properties.getProperty(prefix+"psf_cutoffEnergy")!=null)
				this.psf_cutoffEnergy=Double.parseDouble(properties.getProperty(prefix+"psf_cutoffEnergy"));
			if (properties.getProperty(prefix+"psf_cutoffLevel")!=null)
				this.psf_cutoffLevel=Double.parseDouble(properties.getProperty(prefix+"psf_cutoffLevel"));
			if (properties.getProperty(prefix+"psf_minArea")!=null)
				this.psf_minArea=Integer.parseInt(properties.getProperty(prefix+"psf_minArea"));
			if (properties.getProperty(prefix+"psf_blurSigma")!=null)
				this.psf_blurSigma=Double.parseDouble(properties.getProperty(prefix+"psf_blurSigma"));
			if (properties.getProperty(prefix+"correlationMinInitialContrast")!=null)
				this.correlationMinInitialContrast=Double.parseDouble(properties.getProperty(prefix+"correlationMinInitialContrast"));
			if (properties.getProperty(prefix+"minimalPatternCluster")!=null)
			    this.minimalPatternCluster=Integer.parseInt(properties.getProperty(prefix+"minimalPatternCluster"));
			if (properties.getProperty(prefix+"scaleMinimalInitialContrast")!=null)
			    this.scaleMinimalInitialContrast=Double.parseDouble(properties.getProperty(prefix+"scaleMinimalInitialContrast"));
			if (properties.getProperty(prefix+"searchOverlap")!=null)
			    this.searchOverlap=Double.parseDouble(properties.getProperty(prefix+"searchOverlap"));
			if (properties.getProperty(prefix+"smallestSubPix")!=null)
				this.smallestSubPix=Double.parseDouble(properties.getProperty(prefix+"smallestSubPix"));
			if (properties.getProperty(prefix+"bitmapNonuniforityThreshold")!=null)
				this.bitmapNonuniforityThreshold=Double.parseDouble(properties.getProperty(prefix+"bitmapNonuniforityThreshold"));
			if (properties.getProperty(prefix+"subdiv")!=null)
				this.subdiv=Integer.parseInt(properties.getProperty(prefix+"subdiv"));
			if (properties.getProperty(prefix+"overexposedMaxFraction")!=null)
				this.overexposedMaxFraction=Double.parseDouble(properties.getProperty(prefix+"overexposedMaxFraction"));
			if (properties.getProperty(prefix+"minDefinedArea")!=null)
				this.minDefinedArea=Double.parseDouble(properties.getProperty(prefix+"minDefinedArea"));
			if (properties.getProperty(prefix+"PSFKernelSize")!=null)
				this.PSFKernelSize=Integer.parseInt(properties.getProperty(prefix+"PSFKernelSize"));
			if (properties.getProperty(prefix+"approximateGrid")!=null)
				this.approximateGrid=Boolean.parseBoolean(properties.getProperty(prefix+"approximateGrid"));
			if (properties.getProperty(prefix+"centerPSF")!=null)
				this.centerPSF=Boolean.parseBoolean(properties.getProperty(prefix+"centerPSF"));
			if (properties.getProperty(prefix+"mask1_sigma")!=null)
				this.mask1_sigma=Double.parseDouble(properties.getProperty(prefix+"mask1_sigma"));
			if (properties.getProperty(prefix+"mask1_threshold")!=null)
				this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"mask1_threshold"));
			if (properties.getProperty(prefix+"gaps_sigma")!=null)
				this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"gaps_sigma"));
			if (properties.getProperty(prefix+"mask_denoise")!=null)
				this.mask_denoise=Double.parseDouble(properties.getProperty(prefix+"mask_denoise"));
			if (properties.getProperty(prefix+"deconvInvert")!=null)
				this.deconvInvert=Double.parseDouble(properties.getProperty(prefix+"deconvInvert"));
			if (properties.getProperty(prefix+"correlationSize")!=null)
				this.correlationSize=Integer.parseInt(properties.getProperty(prefix+"correlationSize"));
			if (properties.getProperty(prefix+"correlationGaussWidth")!=null)
				this.correlationGaussWidth=Double.parseDouble(properties.getProperty(prefix+"correlationGaussWidth"));
			if (properties.getProperty(prefix+"minUVSpan")!=null)
				this.minUVSpan=Double.parseDouble(properties.getProperty(prefix+"minUVSpan"));
			if (properties.getProperty(prefix+"flatFieldCorrection")!=null)
				this.flatFieldCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldCorrection"));
			if (properties.getProperty(prefix+"flatFieldExpand")!=null)
				this.flatFieldExpand=Double.parseDouble(properties.getProperty(prefix+"flatFieldExpand"));
			if (properties.getProperty(prefix+"thresholdFinish")!=null)
				this.thresholdFinish=Double.parseDouble(properties.getProperty(prefix+"thresholdFinish"));
			if (properties.getProperty(prefix+"numIterations")!=null)
				this.numIterations=Integer.parseInt(properties.getProperty(prefix+"numIterations"));
			if (properties.getProperty(prefix+"goniometerMotors_ipAddress")!=null)
				this.goniometerMotors.ipAddress=properties.getProperty(prefix+"goniometerMotors_ipAddress");
			if (properties.getProperty(prefix+"goniometerMotors_stepsPerSecond")!=null)
				this.goniometerMotors.stepsPerSecond=Double.parseDouble(properties.getProperty(prefix+"goniometerMotors_stepsPerSecond"));
			if (properties.getProperty(prefix+"goniometerMotors_stepsPerDegreeTilt")!=null)
				this.goniometerMotors.stepsPerDegreeTilt=Double.parseDouble(properties.getProperty(prefix+"goniometerMotors_stepsPerDegreeTilt"));
			if (properties.getProperty(prefix+"goniometerMotors_stepsPerDegreeAxial")!=null)
				this.goniometerMotors.stepsPerDegreeAxial=Double.parseDouble(properties.getProperty(prefix+"goniometerMotors_stepsPerDegreeAxial"));

			if (properties.getProperty(prefix+"targetDistance")!=null)
				this.targetDistance=Double.parseDouble(properties.getProperty(prefix+"targetDistance"));
			if (properties.getProperty(prefix+"scanLatitudeLow")!=null)
				this.scanLatitudeLow=Double.parseDouble(properties.getProperty(prefix+"scanLatitudeLow"));
			if (properties.getProperty(prefix+"scanLatitudeHigh")!=null)
				this.scanLatitudeHigh=Double.parseDouble(properties.getProperty(prefix+"scanLatitudeHigh"));
			if (properties.getProperty(prefix+"scanOverlapHorizontal")!=null)
				this.scanOverlapHorizontal=Double.parseDouble(properties.getProperty(prefix+"scanOverlapHorizontal"));
			if (properties.getProperty(prefix+"scanOverlapVertical")!=null)
				this.scanOverlapVertical=Double.parseDouble(properties.getProperty(prefix+"scanOverlapVertical"));
			if (properties.getProperty(prefix+"scanLimitAxialStart")!=null)
				this.scanLimitAxialStart=Double.parseDouble(properties.getProperty(prefix+"scanLimitAxialStart"));
			if (properties.getProperty(prefix+"scanLimitAxialEnd")!=null)
				this.scanLimitAxialEnd=Double.parseDouble(properties.getProperty(prefix+"scanLimitAxialEnd"));
			if (properties.getProperty(prefix+"scanBidirectional")!=null)
				this.scanBidirectional=Boolean.parseBoolean(properties.getProperty(prefix+"scanBidirectional"));
			if (properties.getProperty(prefix+"motorsSimultaneous")!=null)
				this.motorsSimultaneous=Boolean.parseBoolean(properties.getProperty(prefix+"motorsSimultaneous"));

		}
    	public boolean showDialog(String title) {
    		GenericDialog gd = new GenericDialog(title);
    		//	    		this.serialNumber, // camera serial number string
    		gd.addMessage("Sensor board serial number is "+(((this.serialNumber==null)||(this.serialNumber==""))?"not specified":this.serialNumber));
			gd.addStringField  ("Grid geometry file",                                 this.gridGeometryFile,40);
			gd.addStringField  ("Initial camera intrinsic/extrinsic parametres file", this.initialCalibrationFile,40);
			gd.addStringField  ("Levenberg-Marquardt algorithm strategy file",        this.strategyFile,40);
			gd.addStringField  ("Focusing results superdirectory (individual will be named by serial numbers)", this.resultsSuperDirectory,40);
			gd.addNumericField("EEPROM channel to read sensor serial number from",    this.EEPROM_channel, 0,4,"");
			gd.addCheckbox    ("Save goniometer results (including intermediate) ", this.saveResults);
			gd.addCheckbox    ("Show SFE focusing results (including intermediate) ", this.showResults);
			gd.addStringField  ("Comment to add to the result files",                 this.comment,80);
    		gd.addNumericField("Maximal grid correction between images",  this.maxCorr, 3,5,"pix");
    		gd.addCheckbox    ("Show history details (per color info)",   this.showHistoryDetails);
    		gd.addCheckbox    ("Show history details in a single line (for spreadheets)",  this.showHistorySingleLine);
    		gd.addCheckbox    ("Show acquired images",                    this.showAcquiredImages); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addCheckbox    ("Show LMA fitted parameters",              this.showFittedParameters); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addMessage("When approximating measured PSF for different areas/colors:");
    		gd.addNumericField("Disregard pixels outside of this fraction of the total energy",           100*this.psf_cutoffEnergy, 2,6,"%");
    		gd.addNumericField("Disregard pixels below this fraction of the maximal value",               100*this.psf_cutoffLevel,  2,6,"%");
    		gd.addNumericField("Minimal selection size (will continue even if previous conditions matched)",  this.psf_minArea,      0,3,"sub-pix");
    		gd.addNumericField("Optionally blur the calculated selection mask",                               this.psf_blurSigma,    2,6,"sub-pix");
    		gd.addMessage     ("The following parameters overwrite some defined for aberration measurements in other dialogs");
			gd.addNumericField("Correlation minimal contrast for initial search:", this.correlationMinInitialContrast, 3);
			gd.addNumericField("Minimal initial pattern cluster size (0 - disable retries)", this.minimalPatternCluster, 0);
			gd.addNumericField("Scale minimal contrast if the initial cluster is nonzero but smaller", this.scaleMinimalInitialContrast, 3);
			gd.addNumericField("Overlap of FFT areas when searching for pattern", this.searchOverlap, 3);
    		gd.addNumericField("Smallest fraction to subdivide pixels at simulation", this.smallestSubPix, 3,5,"sensor pix");
    		gd.addNumericField("Maximal difference of the pattern value in the corners that triggers subdivision", this.bitmapNonuniforityThreshold, 3);
    		gd.addNumericField("Subdivide simulated pattern by:",         this.subdiv, 0);
    		gd.addNumericField("Allowed overexposed pixels (fraction of the area) ",this.overexposedMaxFraction,3); //  0.005; // allowed fraction of the overexposed pixels in the PSF kernel measurement area
    		gd.addNumericField("Min fraction of the FFT square (weighted) to have defined pattern",  this.minDefinedArea, 3);
    		gd.addNumericField ("PSF kernel size",                        this.PSFKernelSize, 0);
    		gd.addCheckbox    ("Approximate pattern grid with a polynomial",this.approximateGrid); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addCheckbox    ("Center PSF by modifying phase",           this.centerPSF); // true; // ignore lateral chromatic aberration (center OTF to 0,0)
    		gd.addNumericField("Bluring power spectrum to remove pattern grid (in pattern base freq)",  this.mask1_sigma, 3);
    		gd.addNumericField("Threshold to supress spectral points not present in the pattern ",  this.mask1_threshold, 3);
    		gd.addNumericField("Sigma for filling the OTF ",               this.gaps_sigma, 3);
    		gd.addNumericField("Denoise mask ",                            this.mask_denoise, 3);
    		gd.addNumericField("Invert deconvolution if less than",        this.deconvInvert, 3);
			gd.addNumericField("Correlation size:",                        this.correlationSize, 0); // 64
			gd.addNumericField("Correlation Gauss width (relative):",      this.correlationGaussWidth, 3);
			gd.addNumericField("Minimal UV span in correlation window to trigger FFT size increase",this.minUVSpan, 3);
			gd.addCheckbox    ("Compensate uneven pattern intensity",      this.flatFieldCorrection);
			gd.addNumericField("Expand during extrapolation (relative to the average grid period)", this.flatFieldExpand, 3);
			gd.addNumericField("Threshold RMS to exit LMA",                this.thresholdFinish, 7,9,"pix");
			gd.addNumericField("Maximal number of LMA iterations per series",this.numIterations, 0);
			gd.addMessage("Parameters for scanning/acquisition");

			gd.addStringField ("Goniometer motors IP address",             this.goniometerMotors.ipAddress,40);
			gd.addNumericField("Motors rotation speed ",                   this.goniometerMotors.stepsPerSecond,6,12,"steps/second");
			gd.addNumericField("Motor steps per tilt angular degree (currently negative) ",this.goniometerMotors.stepsPerDegreeTilt,6,12,"steps");
			gd.addNumericField("Motor steps per axial angular degree (currently negative) ",this.goniometerMotors.stepsPerDegreeAxial,6,12,"steps");
//			gd.addNumericField("Tilt scan step (not larger than)",         this.scanStepTilt,2,6,"degrees");
//			gd.addNumericField("Axial scan step (not larger than)",        this.scanStepAxial,2,6,"degrees");



			gd.addNumericField("Distance to target (for overlap calculation)", this.targetDistance,2,7,"mm");

			gd.addNumericField("Horizontal overlap",                   100*this.scanOverlapHorizontal,2,5,"%");
			gd.addNumericField("Vertical overlap",                     100*this.scanOverlapVertical,2,5,"%");
			gd.addNumericField("Lowest camera view latitude to calibrate (nadir=-90)", this.scanLatitudeLow,2,6,"degrees");
			gd.addNumericField("Highest camera view latitude to calibrate (zenith=90)", this.scanLatitudeHigh,2,6,"degrees");
//			gd.addNumericField("Tilt scan start angle",                    this.scanLimitTiltStart,2,6,"degrees");
//			gd.addNumericField("Tilt scan end angle",                      this.scanLimitTiltEnd,2,6,"degrees");
			gd.addNumericField("Axial scan start angle",                   this.scanLimitAxialStart,2,6,"degrees");
			gd.addNumericField("Axial scan end angle",                     this.scanLimitAxialEnd,2,6,"degrees");
	        gd.addCheckbox    ("Axial scan bidirectional",                 this.scanBidirectional);
	        gd.addCheckbox    ("Allow simultaneous operation of motors",   this.motorsSimultaneous);


    		if (!Double.isNaN(this.sensorTemperature)) gd.addMessage("Last measured sensor temperature is "+this.sensorTemperature+" C");
    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
			this.gridGeometryFile=           gd.getNextString();
			this.initialCalibrationFile=     gd.getNextString();
			this.strategyFile=               gd.getNextString();
    		this.resultsSuperDirectory=      gd.getNextString();
    		this.EEPROM_channel=       (int) gd.getNextNumber();
			this.saveResults=                gd.getNextBoolean();
			this.showResults=                gd.getNextBoolean();
//    		this.comment=                    gd.getNextString().replace(' ','_'); //TODO: - add escape
    		this.comment=                    gd.getNextString();
    		this.maxCorr=                    gd.getNextNumber();
    		this.showHistoryDetails=         gd.getNextBoolean();
    		this.showHistorySingleLine=      gd.getNextBoolean();
    		this.showAcquiredImages=         gd.getNextBoolean();
    		this.showFittedParameters=       gd.getNextBoolean();
    		this.psf_cutoffEnergy=      0.01*gd.getNextNumber();
    		this.psf_cutoffLevel=       0.01*gd.getNextNumber();
    		this.psf_minArea=          (int) gd.getNextNumber();
    		this.psf_blurSigma=              gd.getNextNumber();
    		this.correlationMinInitialContrast=gd.getNextNumber();
    		this.minimalPatternCluster=(int) gd.getNextNumber();
    		this.scaleMinimalInitialContrast=gd.getNextNumber();
			this.searchOverlap=              gd.getNextNumber();
    		this.smallestSubPix=             gd.getNextNumber();
    		this.bitmapNonuniforityThreshold=gd.getNextNumber();
    		this.subdiv=               (int) gd.getNextNumber();
    		this.overexposedMaxFraction=     gd.getNextNumber();
    		this.minDefinedArea=             gd.getNextNumber();
    		this.PSFKernelSize=        (int) gd.getNextNumber();
    		this.approximateGrid=            gd.getNextBoolean();
    		this.centerPSF=                  gd.getNextBoolean();
    		this.mask1_sigma=                gd.getNextNumber();
    		this.mask1_threshold=            gd.getNextNumber();
    		this.gaps_sigma=                 gd.getNextNumber();
    		this.mask_denoise=               gd.getNextNumber();
    		this.deconvInvert=               gd.getNextNumber();
			this.correlationSize=      (int) gd.getNextNumber();
			this.correlationGaussWidth=      gd.getNextNumber();
			this.minUVSpan=                  gd.getNextNumber();
			this.flatFieldCorrection=        gd.getNextBoolean();
			this.flatFieldExpand=            gd.getNextNumber();
			this.thresholdFinish=            gd.getNextNumber();
			this.numIterations=        (int) gd.getNextNumber();
			this.goniometerMotors.ipAddress=          gd.getNextString();
			this.goniometerMotors.stepsPerSecond=     gd.getNextNumber();
			this.goniometerMotors.stepsPerDegreeTilt= gd.getNextNumber();
			this.goniometerMotors.stepsPerDegreeAxial=gd.getNextNumber();
			this.targetDistance=             gd.getNextNumber();
			this.scanOverlapHorizontal= 0.01*gd.getNextNumber();
			this.scanOverlapVertical=   0.01*gd.getNextNumber();
			this.scanLatitudeLow=            gd.getNextNumber();
			this.scanLatitudeHigh=           gd.getNextNumber();
			this.scanLimitAxialStart=        gd.getNextNumber();
			this.scanLimitAxialEnd=          gd.getNextNumber();
	        this.scanBidirectional=          gd.getNextBoolean();
	        this.motorsSimultaneous=         gd.getNextBoolean();
    		return true;
    	}
/* ======================================================================== */

    }

}
