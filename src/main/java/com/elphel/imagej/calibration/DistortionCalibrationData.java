package com.elphel.imagej.calibration;
/*
 **
 ** DistortionCalibrationData.java
 **
 ** Copyright (C) 2011-2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  DistortionCalibrationData.java is free software: you can redistribute it and/or modify
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

import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.configuration.XMLConfiguration;

import com.elphel.imagej.cameras.EyesisCameraParameters;
import com.elphel.imagej.cameras.EyesisSubCameraParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.PolynomialApproximation;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.common.WindowTools;
import com.elphel.imagej.jp4.JP46_Reader_camera;

import Jama.Matrix;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

// stores per pattern image camera/subcamera parameters, filenames, ...
// saves / restores them from a disk file
    public class DistortionCalibrationData{
    	public static final int INDEX_PX =       0;
    	public static final int INDEX_PY =       1;
    	public static final int INDEX_U =        2;
    	public static final int INDEX_V =        3;
    	public static final int INDEX_CONTRAST = 4;
    	public static final int INDEX_R =        5;
    	public static final int INDEX_G =        6;
    	public static final int INDEX_B =        7;
    	public static final double SMALL_FRACTION = 0.8; // consider sensor to be a "small" if average grid period < this fraction of the large

    	Goniometer.GoniometerParameters goniometerParameters = null;
    	public String pathName=null;
    	public EyesisCameraParameters eyesisCameraParameters; // has "cartesian"
        public int       numSubCameras=1;
        public int       numPointers=4;    // maximal number of pointers to look for
        public int       numMotors  =3;    // maximal number of motors to look for
        public GridImageParameters [] gIP= null; // per-grid image parameters
        public GridImageSet []        gIS= null; // sets of images with the same timestamp
        public boolean [] small_sensors =    null; // set by filter grids
        public double     small_period_frac =   0; // set by filter grids - ratio of small sensor period to large sensor period
        // keep for now?
    	public double [][] pars=null; // for each defined image: set of (22) parameters
    	public double [][] sensorMasks= null; // per-channel (not image) mask

    	//pixelsXY, pixelsUV should match, second dimension is variable
    	public boolean updateStatus=true;
    	public int     debugLevel=2;
    	private ShowDoubleFloatArrays SDFA_INSTANCE=null; // just for debugging
    	public int getNumStations(){
    		return (eyesisCameraParameters==null)?0:eyesisCameraParameters.getNumStations();
    	}

		public double getPixelSize(int station, int channel) {
			return this.eyesisCameraParameters.eyesisSubCameras[station][channel].getPixelSize();
		}
		public double getDistortionRadius(int station, int channel) {
			return this.eyesisCameraParameters.eyesisSubCameras[station][channel].getDistortionRadius();
		}
		public double getPixelSize(int imgNumber) {
			return getPixelSize(this.gIP[imgNumber].stationNumber, this.gIP[imgNumber].channel);
		}
		public double getDistortionRadius(int imgNumber) {
			return getDistortionRadius(this.gIP[imgNumber].stationNumber, this.gIP[imgNumber].channel);
		}
        public boolean isTripod() {
        	return (eyesisCameraParameters !=null) && this.eyesisCameraParameters.isTripod();
        }

        public boolean isCartesian() {
        	return (eyesisCameraParameters !=null) && this.eyesisCameraParameters.isCartesian();
        }


     	public class GridImageParameters{
    		public int         imgNumber=-1; // index of this image (for pars[][])
    		private int        setNumber=-1; // long overdue  - will be some inconsistency
    		GridImageSet       gridImageSet=null;
    		public int         stationNumber=0; // changes when camera/goniometer is moved to new position
    		public String      path=null;
    		public String      source_path = null; // Full path of the source image this grid was calculated from
    		public double [][] laserPixelCoordinates=null; // first index - absolute number of pointer. Each element may be either null or {x,y} pair
    		// moving for new files to have laser UV contained in the file
    		public double [][] laserUVCoordinates=null;    // first index - absolute number of pointer. Each element may be either null or {u,v} pair - never used??????
    		public int         matchedPointers=0;
    		public int         hintedMatch=-1; // -1 - not tried, 0 - no real grid (i.e. double reflection), applied orientation, applied orientation and shift
    		public boolean     enabled=true; //false;  // to mask out some images from all strategy steps (i.e w/o reliable absolute calibration)
    		public boolean     flatFieldAvailable=false; // grid files have flat field data
    		public boolean     newEnabled=false;
    		public int []      motors=null;
    		public ImagePlus   gridImage=null;
    		public double      timestamp=-1;
    		public int         channel=  -1;
    		public double []   intensityRange={255.0,255.0,255.0}; // r,g,b - used to normalize vign*
    		public double []   saturation={255.0,255.0,255.0}; // r,g,b - saturation range read from images
//    		public double  []  pars=null; // set of (22) parameters
    		public double [][] pixelsXY=   null; // for each image, each grid node - a set of of {px,py,contrast,vignR,vignG,vignB} vign* is in the 0..1.0 range
    		public double []   pixelsMask= null; // for each image, each grid node - weight function derived from contrast and 3 parameters
    		public int    [][] pixelsUV=  null; // for each image, each grid node - a pair of {gridU, gridV}
    		public boolean  [] badNodes=  null; // if not null, marks node with excessive errors
    		public double [][] pixelsXY_extra=  null; // extra data, for nodes that are out of the physical grid (may be needed after re-calibration)
    		public int    [][] pixelsUV_extra=  null;
    		private double      gridPeriod=0.0;  // average grid period, in pixels (to filter out (double-) reflected images
    		public boolean     noUsefulPSFKernels=false; // used to mark images w/o good PSF data
    		public double      diameter=0.0;
    		public int []      UVShiftRot={0,0,0}; // shift and rotation of the grid
    		public Rectangle   woi;
    		final int contrastIndex=2;



    		public double getGridPeriod() {	return gridPeriod;}
    		public void setGridPeriod(double v) {gridPeriod = v;}
    		public int getSetNumber(){return this.setNumber;}
    		public int getImageNumber() {return this.imgNumber;}
        	public GridImageParameters(int index){
        		this.imgNumber=index;
        	}
            public int [] getUVShiftRot(){
            	return this.UVShiftRot;
            }
            public void setUVShiftRot(int [] UVShiftRot){
            	this.UVShiftRot=UVShiftRot;
            }
        	public int getStationNumber(){ // TODO: make only a single station number - in GridImageSet?
        		return this.stationNumber;
        	}
        	public void setStationNumber(int stationNumber){ // TODO: make only a single station number - in GridImageSet?
        		this.stationNumber=stationNumber;
        	}
        	public double [] getGridWeight(){
        		return this.pixelsMask;
        	}
        	public void resetMask(){
        			this.pixelsMask=null;
    	    }
        	public void resetBadNodes(){
        		this.badNodes=null;
        	}
        	public void setBadNode(int index){
        		if (this.badNodes==null){
        			this.badNodes=new boolean[this.pixelsXY.length]; // let it throw if null
        			for (int i=0;i<this.badNodes.length;i++)this.badNodes[i]=false;
        		}
        		this.badNodes[index]=true;
        	}
        	public boolean isNodeBad(int index){
        		if (this.badNodes==null) return false;
        		if (index>=this.badNodes.length) {
        			System.out.println("### isNodeBad("+index+") - OOB, as this.badNodes="+this.badNodes.length);
        			return true;
        		}
        		return this.badNodes[index]; //OOB
        	}

        	public int getNumContrastNodes(double minContrast){
        	    int num=0;
        		for (int i=0;i<this.pixelsXY.length;i++) if (this.pixelsXY[i][contrastIndex]>=minContrast) num++;
        		return num;
        	}

        	public int getChannel() {
        		return channel;
        	}


        	/**
        	 * Calculate "diameter" of the image to be used for image weight
        	 * @param xc image center pixel X
        	 * @param yc image center pixel Y
        	 * @param r0 reference diameter
        	 * @param minContrast minimal contrast to count the node
        	 */

        	public void setImageDiameter( // need to get image center px,py. Maybe r0 - use to normalize result diameter
        			double xc,
        			double yc,
        			double r0,
        			double minContrast,
        			int dbgImgNum //>=0 - debug print with image number
        			){
        		boolean debug=(dbgImgNum>=0);
        		// find the farthest point from the center
        		double maxR2=-1;
        		int firstIndex=0;
        		for (int i=0;i<this.pixelsXY.length;i++) if (this.pixelsXY[i][contrastIndex]>=minContrast){
        			double dx=this.pixelsXY[i][0]-xc;
        			double dy=this.pixelsXY[i][1]-yc;
        			double r2=dx*dx+dy*dy;
        			if (r2>maxR2) {
        				maxR2=r2;
        				firstIndex=i;
        			}
        		}
        		if (maxR2<=0) {
        			this.diameter=0.0;
        			return;
        		}
        		double maxDx=this.pixelsXY[firstIndex][0]-xc;
        		double maxDy=this.pixelsXY[firstIndex][1]-yc;

        		if (debug) System.out.print("setImageDiameter("+IJ.d2s(xc,2)+","+IJ.d2s(yc,2)+","+IJ.d2s(r0,2)+","+IJ.d2s(minContrast,4)+","+dbgImgNum+") ---- > ");
        		if (debug) System.out.print(" maxR2="+IJ.d2s(maxR2,2)+" maxDx="+IJ.d2s(maxDx,2)+" maxDy="+IJ.d2s(maxDy,2));
        		double maxAamb=0;
        		double dbgDx=0.0,dbgDy=0.0;
        		for (int i=0;i<this.pixelsXY.length;i++) if (this.pixelsXY[i][contrastIndex]>=minContrast){
//        			double dx=maxDx-this.pixelsXY[i][0];
//        			double dy=maxDy-this.pixelsXY[i][1];
        			double dx=this.pixelsXY[firstIndex][0]-this.pixelsXY[i][0];
        			double dy=this.pixelsXY[firstIndex][1]-this.pixelsXY[i][1];
        			double aAmb=dx*maxDx+dy*maxDy;
        			if (aAmb>maxAamb) {
        				maxAamb=aAmb;
        				dbgDx=this.pixelsXY[i][0]; // debug only !
        				dbgDy=this.pixelsXY[i][1]; // debug only !
        			}
        		}
        		this.diameter=maxAamb/Math.sqrt(maxR2)/r0;
        		if (debug) System.out.println(" maxAamb="+IJ.d2s(maxAamb,2)+" dbgDx="+IJ.d2s(dbgDx,2)+" dbgDy="+IJ.d2s(dbgDy,2)+" --> "+IJ.d2s(this.diameter,2));
        	}
        	/**
        	 * Uzses data calculated by  setImageDiameter();
        	 * @return detected grid diameter (along the radius) to be uses as image weight (in r0 units)
        	 */
        	public double getGridDiameter(){
        		return this.diameter;
        	}

        	public void calculateMask(
        			double minContrast,
        			double shrinkBlurSigma,
        			double shrinkBlurLevel){
        		if (this.pixelsMask!=null) return; // need to reset to re-calculate
        		if (this.pixelsUV==null) {this.pixelsMask=null; return; }
        		if (this.pixelsUV.length==0){ this.pixelsMask=new double[0]; return; }

        		this.pixelsMask=new double [this.pixelsUV.length];
        		if (shrinkBlurSigma<=0){
            		for (int i=0;i<this.pixelsUV.length;i++){
            			this.pixelsMask[i]=(this.pixelsXY[i][contrastIndex]>=minContrast)?1.0:0.0;
            		}
            		return;
        		}
        		int minU=this.pixelsUV[0][0],minV=this.pixelsUV[0][1];
        		int maxU=minU,maxV=minV;
        		int margin=(int) (2*shrinkBlurSigma);
        		for (int i=0;i<this.pixelsUV.length;i++){
        			if (this.pixelsUV[i][0]>maxU) maxU=this.pixelsUV[i][0];
        			if (this.pixelsUV[i][0]<minU) minU=this.pixelsUV[i][0];
        			if (this.pixelsUV[i][1]>maxV) maxV=this.pixelsUV[i][1];
        			if (this.pixelsUV[i][1]<minV) minV=this.pixelsUV[i][1];
        		}

        		int U0=minU-margin;
        		int V0=minV-margin;
        		int width= (maxU-minU+1+2*margin);
        		int height=(maxV-minV+1+2*margin);
        		double [] mask = new double [width*height];
        		for (int i=0;i<mask.length;i++) mask[i]=-1.0;
        		for (int i=0;i<this.pixelsUV.length;i++){
        			int index=(this.pixelsUV[i][0]-U0)+width*(this.pixelsUV[i][1]-V0);
        			mask[index]=(this.pixelsXY[i][contrastIndex]>=minContrast)?1.0:-1.0; // java.lang.ArrayIndexOutOfBoundsException: 2230
        		}
        		(new DoubleGaussianBlur()).blurDouble(
							mask,
							width,
							height,
							shrinkBlurSigma,
							shrinkBlurSigma,
							0.01);
        		double k=1.0/(1.0-shrinkBlurLevel);
        		double dbgMax=0.0;
        		for (int i=0;i<this.pixelsUV.length;i++){
        			int index=(this.pixelsUV[i][0]-U0)+width*(this.pixelsUV[i][1]-V0);
        			double d=k*(mask[index]-shrinkBlurLevel);
        			this.pixelsMask[i]=(d>0.0)?(d*d):0.0;
        			if (this.pixelsMask[i]>dbgMax) dbgMax=this.pixelsMask[i];
        		}
 //      		System.out.print(" "+IJ.d2s(dbgMax,2)+" ");
        	}
    	}
    	public class GridImageSet{
    		private int numPars=53; // 27;
    		private int thisParsStartIndex=6;

    		public int         stationNumber=0; // changes when camera/goniometer is moved to new position
    		public GridImageParameters [] imageSet=null;
//    		public GridImageParameters firstImage=null; // first non-null image in the sert (update to have current parameters?)
    		public double timeStamp;
    		public int [] motors=null;
    		public double goniometerAxial=Double.NaN;
    		public double goniometerTilt=Double.NaN;
    		public double interAxisDistance;    // 8 distance in mm between two goniometer axes
    		public double interAxisAngle;       // 9 angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
    		public double horAxisErrPhi;        //10 angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
    		public double horAxisErrPsi;        //11 angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
    		public double entrancePupilForward; //12 common to all lenses - distance from the sensor to the lens entrance pupil
    		public double centerAboveHorizontal;//13 camera center distance along camera axis above the closest point to horizontal rotation axis (adds to height of each
    		public double [] GXYZ=new double [3];  //14 (12) coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
//			this.GXYZ[stationNumber][1],              //15 (13)  y
//			this.GXYZ[stationNumber][2],              //16 (14)  z
    		public boolean orientationEstimated=true; // orientation is estimated from other stes, not adjusted by LMA
    		public double setWeight=0.0; // weight of this set when calculating errors
    		public void setEstimatedFromNonNaN(){
    			this.orientationEstimated= Double.isNaN(this.goniometerTilt) ||  Double.isNaN(this.goniometerAxial);
    		}
    		public int getMinIndex(){
    			return this.thisParsStartIndex;
    		}
    		public int getMaxIndexPlusOne(){
    			return this.thisParsStartIndex+getSetVector().length;
    		}

    		public double [] getSetVector(){
    			double [] sv={
    		    		this.goniometerTilt,
    		    		this.goniometerAxial,
    		    		this.interAxisDistance,
    		    		this.interAxisAngle,
    		    		this.horAxisErrPhi,
    		    		this.horAxisErrPsi,
    		    		this.entrancePupilForward,
    		    		this.centerAboveHorizontal,
    		    		this.GXYZ[0],
    		    		this.GXYZ[1],
    		    		this.GXYZ[2]
    			};
    			return sv;
    		}
    		public void setSetVector(double [] vector){
    			if (vector.length!=getSetVector().length){
    				String msg="Wrong parameter vector length - got:"+vector.length+", expected: "+getSetVector().length;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
	    		this.goniometerTilt=       vector[ 0];
	    		this.goniometerAxial=      vector[ 1];
	    		this.interAxisDistance=    vector[ 2];
	    		this.interAxisAngle=       vector[ 3];
	    		this.horAxisErrPhi=        vector[ 4];
	    		this.horAxisErrPsi=        vector[ 5];
	    		this.entrancePupilForward= vector[ 6];
	    		this.centerAboveHorizontal=vector[ 7];
	    		this.GXYZ[0]=              vector[ 8];
	    		this.GXYZ[1]=              vector[ 9];
	    		this.GXYZ[2]=              vector[10];
    		}

    		public double getParameterValue(int index){
    			int thisIndex=index-this.thisParsStartIndex;
    			double [] sv=getSetVector();
    			if ((thisIndex<0) || (index >sv.length)) return Double.NaN;
    			return sv[thisIndex];

    		}
    		public void setParameterValue(int index,
    				double value,
    				boolean updateEstimated){
    			int thisIndex=index-this.thisParsStartIndex;
    			switch (thisIndex){
    			case  0:
    				this.goniometerTilt=       value;
    				setEstimatedFromNonNaN();
    				break;
    			case  1:
    				this.goniometerAxial=      value;
    				setEstimatedFromNonNaN();
    				break;
    			case  2: this.interAxisDistance=    value; break;
    			case  3: this.interAxisAngle=       value; break;
    			case  4: this.horAxisErrPhi=        value; break;
    			case  5: this.horAxisErrPsi=        value; break;
    			case  6: this.entrancePupilForward= value; break;
    			case  7: this.centerAboveHorizontal=value; break;
    			case  8: this.GXYZ[0]=              value; break;
    			case  9: this.GXYZ[1]=              value; break;
    			case 10: this.GXYZ[2]=              value; break;
    			}
    		}

    		public double [] updateParameterVectorFromSet(double [] vector){
    			if (vector==null){
    				vector=new double [this.numPars];
    				for (int i=0;i<vector.length;i++) vector[i]=Double.NaN;
    			}
    			if (vector.length!=this.numPars){
    				String msg="Wrong parameter vector length - got:"+vector.length+", expected: "+this.numPars;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
    			double [] sv=getSetVector();
    			for (int i=0;i<sv.length;i++) if (!Double.isNaN(sv[i])) vector[i+this.thisParsStartIndex]=sv[i];
    			return vector;
    		}
    		public double [] updateParameterVectorFromSet(double [] vector, boolean [] mask){
    			if (vector==null){
    				vector=new double [this.numPars];
    				for (int i=0;i<vector.length;i++) vector[i]=Double.NaN;
    			}
    			if (vector.length!=this.numPars){
    				String msg="Wrong parameter vector length - got:"+vector.length+", expected: "+this.numPars;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
    			double [] sv=getSetVector();
    			for (int i=0;i<sv.length;i++) if (!Double.isNaN(sv[i]) && mask[this.thisParsStartIndex+ i]) vector[i+this.thisParsStartIndex]=sv[i];
    			return vector;
    		}

    		public void updateSetFromParameterVector(double [] vector){
    			if (vector.length!=this.numPars){
    				String msg="Wrong parameter vector length - got:"+vector.length+", expected: "+this.numPars;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
	    		this.goniometerTilt=       vector[this.thisParsStartIndex+ 0];
	    		this.goniometerAxial=      vector[this.thisParsStartIndex+ 1];
	    		this.interAxisDistance=    vector[this.thisParsStartIndex+ 2];
	    		this.interAxisAngle=       vector[this.thisParsStartIndex+ 3];
	    		this.horAxisErrPhi=        vector[this.thisParsStartIndex+ 4];
	    		this.horAxisErrPsi=        vector[this.thisParsStartIndex+ 5];
	    		this.entrancePupilForward= vector[this.thisParsStartIndex+ 6];
	    		this.centerAboveHorizontal=vector[this.thisParsStartIndex+ 7];
	    		this.GXYZ[0]=              vector[this.thisParsStartIndex+ 8];
	    		this.GXYZ[1]=              vector[this.thisParsStartIndex+ 9];
	    		this.GXYZ[2]=              vector[this.thisParsStartIndex+10];
    		}

    		public void updateSetFromParameterVector(double [] vector, boolean [] mask){
    			if (vector.length!=this.numPars){
    				String msg="Wrong parameter vector length - got:"+vector.length+", expected: "+this.numPars;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
	    		if (mask[this.thisParsStartIndex+ 0]) this.goniometerTilt=       vector[this.thisParsStartIndex+ 0];
	    		if (mask[this.thisParsStartIndex+ 1]) this.goniometerAxial=      vector[this.thisParsStartIndex+ 1];
	    		if (mask[this.thisParsStartIndex+ 2]) this.interAxisDistance=    vector[this.thisParsStartIndex+ 2];
	    		if (mask[this.thisParsStartIndex+ 3]) this.interAxisAngle=       vector[this.thisParsStartIndex+ 3];
	    		if (mask[this.thisParsStartIndex+ 4]) this.horAxisErrPhi=        vector[this.thisParsStartIndex+ 4];
	    		if (mask[this.thisParsStartIndex+ 5]) this.horAxisErrPsi=        vector[this.thisParsStartIndex+ 5];
	    		if (mask[this.thisParsStartIndex+ 6]) this.entrancePupilForward= vector[this.thisParsStartIndex+ 6];
	    		if (mask[this.thisParsStartIndex+ 7]) this.centerAboveHorizontal=vector[this.thisParsStartIndex+ 7];
	    		if (mask[this.thisParsStartIndex+ 8]) this.GXYZ[0]=              vector[this.thisParsStartIndex+ 8];
	    		if (mask[this.thisParsStartIndex+ 9]) this.GXYZ[1]=              vector[this.thisParsStartIndex+ 9];
	    		if (mask[this.thisParsStartIndex+10]) this.GXYZ[2]=              vector[this.thisParsStartIndex+10];
    		}

    		public double getSetWeight(){return this.setWeight;}
        	public int getStationNumber(){ // TODO: make only a single station number - in GridImageSet?
        		return this.stationNumber;
        	}
        	public void setStationNumber(int stationNumber){ // TODO: make only a single station number - in GridImageSet?
        		this.stationNumber=stationNumber;
        	}
    	}

    	public int index_right; // =     getParameterIndexByName("subcamRight");          // 0
    	public int index_forward; // =   getParameterIndexByName("subcamForward");        // 1
    	public int index_azimuth; // =   getParameterIndexByName("subcamAzimuth");        // 0
    	public int index_heading; // =   getParameterIndexByName("subcamHeading");        // 3
    	public int index_elevation; // = getParameterIndexByName("subcamElevation");      // 4
    	public int index_gh; //=         getParameterIndexByName("goniometerHorizontal"); // 6
    	public int index_ga; //=         getParameterIndexByName("goniometerAxial");      // 7

        public String [][] parameterDescriptionsCartesian ={ // may be shorter, have null rows, shorter rows - will use parameterDescriptions for missing data
        		{"subcamRight",          "Subcamera distance from the vertical rotation axis, positive - right looking to the target","mm","S","E"},                  // 0
        		{"subcamForward",        "Subcamera distance from the vertical rotation axis, positive - towards the target","mm","S","E"},                               // 1
        		null,                                                   // 2
        		{"subcamHeading",        "Optical axis heading (0 - to the target, positive - CW looking from top)","degrees","S","E"}};               // 3

        public String [][] parameterDescriptions ={
        		{"subcamAzimuth",        "Subcamera azimuth, clockwise looking from top","degrees","S","E"},                                    // 0
        		{"subcamDistance",       "Subcamera distance from the axis","mm","S","E"},                                                      // 1
        		{"subcamHeight",         "Subcamera height from the 'equator'","mm","S","E"},                                                   // 2
        		{"subcamHeading",        "Optical axis heading (relative to azimuth)","degrees","S","E"},                                       // 3
        		{"subcamElevation",      "Optical axis elevation (up from equator)","degrees","S","E"},                                         // 4
        		{"subcamRoll",           "Subcamera roll, positive CW looking to the target","degrees","S","E"},                                // 5
    			{"goniometerHorizontal", "Goniometer rotation around 'horizontal' axis (tilting from the target - positive)","degrees","R","E"},// 6
    			{"goniometerAxial",      "Rotation around Eyesis main axis (clockwise in plan - positive)","degrees","R","E"},                  // 7
    			{"interAxisDistance",    "Distance between goniometer axes","mm","C","E"},                                                      // 8
    			{"interAxisAngle",       "Angle error between goniometer axes (<0 if vertical axis rotated CW )","degrees","C","E"},            // 9
    			{"horAxisErrPhi",        "Horizontal axis azimuth error (CW in plan)","degrees","C","E"},                                       //10
    			{"horAxisErrPsi",        "Horizontal axis roll error (CW looking to target)","degrees","C","E"},                                //11
    			{"entrancePupilForward", "Distance from the sensor to the lens entrance pupil","mm","C","E"},                              //12
    			{"centerAboveHorizontal","CenterAboveHorizontal","mm","C","E"},                                                            //13
    			{"GXYZ0",                "Goniometer reference point position X (target coordinates, left)","mm","T","E"},                      //14 (12)
    			{"GXYZ1",                "Goniometer reference point position Y (target coordinates, up)","mm","T","E"},                        //15 (13)
    			{"GXYZ2",                "Goniometer reference point position Z (target coordinates, away)","mm","T","E"} ,                     //16 (14)
    			{"subcamFocalLength",    "Lens focal length","mm","S","I"},                                                                     //17 (15)
    			{"subcamPX0",            "Lens axis on the sensor (horizontal, from left edge)","pixels","S","I"},                              //18 (16)
    			{"subcamPY0",            "Lens axis on the sensor (vertical, from top edge)","pixels","S","I"},                                 //19 (17)
    			{"subcamDistortionA8",   "Distortion A8(r^5)","relative","S","I"},                                                              //20 (18)
    			{"subcamDistortionA7",   "Distortion A7(r^5)","relative","S","I"},                                                              //21 (19)
    			{"subcamDistortionA6",   "Distortion A6(r^5)","relative","S","I"},                                                              //22 (20)
    			{"subcamDistortionA5",   "Distortion A5(r^5)","relative","S","I"},                                                              //23 (21)
    			{"subcamDistortionA",    "Distortion A (r^4)","relative","S","I"},                                                              //24 (22)
    			{"subcamDistortionB",    "Distortion B (r^3)","relative","S","I"},                                                              //25 (23)
    			{"subcamDistortionC",    "Distortion C (r^2)","relative","S","I"},                                                               //26 (24)

        		{"subcamElong_C_o",      "Orthogonal elongation for r^2","relative","S","I"},     // 27 39 (37)
        		{"subcamElong_C_d",      "Diagonal   elongation for r^2","relative","S","I"},     // 28 40 (38)

        		{"subcamEccen_B_x",      "Distortion center shift X for r^3","relative","S","I"}, // 29 27 (25)
        		{"subcamEccen_B_y",      "Distortion center shift Y for r^3","relative","S","I"}, // 30 28 (26)
        		{"subcamElong_B_o",      "Orthogonal elongation for r^3","relative","S","I"},     // 31 41 (39)
        		{"subcamElong_B_d",      "Diagonal   elongation for r^3","relative","S","I"},     // 32 42 (40)

        		{"subcamEccen_A_x",      "Distortion center shift X for r^4","relative","S","I"}, // 33 29 (27)
        		{"subcamEccen_A_y",      "Distortion center shift Y for r^4","relative","S","I"}, // 34 30 (28)
        		{"subcamElong_A_o",      "Orthogonal elongation for r^4","relative","S","I"},     // 35 43 (41)
        		{"subcamElong_A_d",      "Diagonal   elongation for r^4","relative","S","I"},     // 36 44 (42)

        		{"subcamEccen_A5_x",     "Distortion center shift X for r^5","relative","S","I"}, // 37 31 (29)
        		{"subcamEccen_A5_y",     "Distortion center shift Y for r^5","relative","S","I"}, // 38 32 (30)
        		{"subcamElong_A5_o",     "Orthogonal elongation for r^5","relative","S","I"},     // 39 45 (43)
        		{"subcamElong_A5_d",     "Diagonal   elongation for r^5","relative","S","I"},     // 40 46 (44)

        		{"subcamEccen_A6_x",     "Distortion center shift X for r^6","relative","S","I"}, // 41 33 (31)
        		{"subcamEccen_A6_y",     "Distortion center shift Y for r^6","relative","S","I"}, // 42 34 (32)
        		{"subcamElong_A6_o",     "Orthogonal elongation for r^6","relative","S","I"},     // 43 47 (45)
        		{"subcamElong_A6_d",     "Diagonal   elongation for r^6","relative","S","I"},     // 44 48 (46)

        		{"subcamEccen_A7_x",     "Distortion center shift X for r^7","relative","S","I"}, // 45 35 (33)
        		{"subcamEccen_A7_y",     "Distortion center shift Y for r^7","relative","S","I"}, // 46 36 (34)
        		{"subcamElong_A7_o",     "Orthogonal elongation for r^7","relative","S","I"},     // 47 49 (47)
        		{"subcamElong_A7_d",     "Diagonal   elongation for r^7","relative","S","I"},     // 48 50 (48)

        		{"subcamEccen_A8_x",     "Distortion center shift X for r^8","relative","S","I"}, // 49 37 (35)
        		{"subcamEccen_A8_y",     "Distortion center shift Y for r^8","relative","S","I"}, // 50 38 (36)
        		{"subcamElong_A8_o",     "Orthogonal elongation for r^8","relative","S","I"},     // 51 51 (49)
        		{"subcamElong_A8_d",     "Diagonal   elongation for r^8","relative","S","I"}      // 52 52 (50)
        };

        public String [] channelSuffixes={ // natural order (same as array indices, may be modified to camera/subcamera
        		"00","01","02","03","04","05","06","07","08","09",
        		"10","11","12","13","14","15","16","17","18","19",
        		"20","21","22","23","24","25","26","27","28","29"};


        public int getNumSets() {
        	return gIS.length;
        }

        public void setupIndices(){ // should be always called during initialization !
        	this.index_right =     getParameterIndexByName("subcamRight");          // 0 may be -1 if !cartesian
        	this.index_forward =   getParameterIndexByName("subcamForward");        // 1 may be -1 if !cartesian
        	this.index_azimuth =   getParameterIndexByName("subcamAzimuth");        // 0
        	this.index_heading =   getParameterIndexByName("subcamAzimuth");        // 3
        	this.index_elevation = getParameterIndexByName("subcamElevation");      // 4
        	this.index_gh=         getParameterIndexByName("goniometerHorizontal"); // 6
        	this.index_ga=         getParameterIndexByName("goniometerAxial");      // 7
        }

        public String descrField(int i,int j){
        	if (
        			(eyesisCameraParameters !=null) &&
        			eyesisCameraParameters.isCartesian() &&
        			(i < this.parameterDescriptionsCartesian.length) &&
        			(this.parameterDescriptionsCartesian[i]!=null) &&
        			(j<this.parameterDescriptionsCartesian[i].length)){
        		return this.parameterDescriptionsCartesian[i][j];
        	}
    		return this.parameterDescriptions[i][j];
        }

        public boolean isNonRadial(int index){
        	return parameterDescriptions[index][0].startsWith("subcamEccen_") || parameterDescriptions[index][0].startsWith("subcamElong_");
        }

        public int getParameterIndexByName(String name){
        	if (isCartesian()){
            	for (int i=0;i<this.parameterDescriptionsCartesian.length;i++) if ((this.parameterDescriptionsCartesian[i]!=null) && this.parameterDescriptionsCartesian[i][0].equals(name)){
            		return i;
            	}
        	}
        	for (int i=0;i<this.parameterDescriptions.length;i++) if (this.parameterDescriptions[i][0].equals(name)){
        		return i;
        	}
        	return -1;
        }

        public int getNumDescriptions(){
        	return this.parameterDescriptions.length;
        }

/**
 * Initialize data from scratch using filenames "grid-<timestamp-seconds>_<timestamp-microseconds>-<channel-number>.tiff
 * @param filenames List of grid filenames (2-slice TIFFs)
 */

        public DistortionCalibrationData (
        		String [] filenames,
        		PatternParameters patternParameters,
        		EyesisCameraParameters eyesisCameraParameters,
        		LaserPointer laserPointers,
        		Goniometer.GoniometerParameters goniometerParameters
        		) {
    	    String [][] stationFilenames={filenames};
        	setupDistortionCalibrationData(
        			stationFilenames,
            		patternParameters,
            		eyesisCameraParameters, // debugLevel
            		laserPointers,
            		goniometerParameters
            		);
        }
        public DistortionCalibrationData (
        		String []              filenames,
        		PatternParameters      patternParameters,
        		EyesisCameraParameters eyesisCameraParameters,
        		LaserPointer           laserPointers,
        		Goniometer.GoniometerParameters goniometerParameters,
        		int debugLevel
        		) {
        	    this.debugLevel=debugLevel;
        	    String [][] stationFilenames={filenames};
        	setupDistortionCalibrationData(
        			stationFilenames,
            		patternParameters,
            		eyesisCameraParameters, // debugLevel
            		laserPointers,
            		goniometerParameters
            		);
        }

        public DistortionCalibrationData (
        		String [][] stationFilenames,
        		PatternParameters patternParameters,
        		EyesisCameraParameters eyesisCameraParameters,
        		LaserPointer laserPointers,
        		Goniometer.GoniometerParameters goniometerParameters,
        		int debugLevel
        		) {
        	    this.debugLevel=debugLevel;
        	setupDistortionCalibrationData(
        			stationFilenames,
            		patternParameters,
            		eyesisCameraParameters, // debugLevel
            		laserPointers,
            		goniometerParameters
            		);
        }

        public DistortionCalibrationData (
        		String [][] stationFilenames,
        		String []                                              source_dirs,      // directories of the source files per station
        		MultipleExtensionsFileFilter gridFilter,
        		MultipleExtensionsFileFilter sourceFilter,
        		PatternParameters                                      patternParameters,
        		EyesisCameraParameters                                 eyesisCameraParameters,
        		LaserPointer                                           laserPointers, // as a backup if data is not available in the file
        		Goniometer.GoniometerParameters                        goniometerParameters,
        		boolean                                                read_grids,
        		int                                                    debugLevel
        		) {
        	    this.debugLevel=debugLevel;
        	    setupDirDistortionCalibrationData(
        			stationFilenames,
            		source_dirs,      // directories of the source files per station
        			gridFilter,
        			sourceFilter,
            		patternParameters,
            		eyesisCameraParameters, // debugLevel
            		laserPointers, // as a backup if data is not available in the file
            		goniometerParameters,
            		read_grids
            		);
        }


        public void setupDistortionCalibrationData (
        		String [][]            stationFilenames,
        		PatternParameters      patternParameters,
        		EyesisCameraParameters eyesisCameraParameters,
        		LaserPointer           laserPointers, // as a backup if data is not available in the file
        		Goniometer.GoniometerParameters goniometerParameters
        		) {
        	this.goniometerParameters = goniometerParameters;
        	setupIndices();
        	this.eyesisCameraParameters=eyesisCameraParameters;
        	int numSubCameras=(eyesisCameraParameters==null)?1:eyesisCameraParameters.eyesisSubCameras[0].length;
        	this.numSubCameras=numSubCameras;
        	this.eyesisCameraParameters.numStations=stationFilenames.length;
        	int numFiles=0;
        	for (int i=0;i<stationFilenames.length;i++) numFiles+=stationFilenames[i].length;
        	this.gIP=new GridImageParameters[numFiles];


        	int numFile=0;
        	for (int numStation=0;numStation<stationFilenames.length;numStation++){
        		for (int index=0;index<stationFilenames[numStation].length;index++){

        			System.out.println(numFile+" ("+numStation+":"+index+"): "+stationFilenames[numStation][index]);
        			this.gIP[numFile]=new GridImageParameters(numFile);
        			this.gIP[numFile].path=stationFilenames[numStation][index]; //Exception in thread "Run$_AWT-EventQueue-0" java.lang.NullPointerException at Distortions$DistortionCalibrationData.<init>(Distortions.java:5987)
        			this.gIP[numFile].setStationNumber(numStation);
        			int i1=stationFilenames[numStation][index].indexOf('-',stationFilenames[numStation][index].lastIndexOf(Prefs.getFileSeparator()));
        			int i2=stationFilenames[numStation][index].indexOf('-',i1+1);
        			int i3=stationFilenames[numStation][index].indexOf('.',i2+1);
        			// Extract timestamp from the filename
        			if ((i1<0) || (i2<0)) {
        				String msg="invalid file format - '"+stationFilenames[numStation][index]+"', should be '<timestamp-seconds>_<timestamp-microseconds>-<channel-number>.tiff'";
        				IJ.showMessage("Error",msg);
        				throw new IllegalArgumentException (msg);
        			}
        			// Extract channel number from the filename

        			this.gIP[numFile].timestamp=Double.parseDouble(stationFilenames[numStation][index].substring(i1+1,i2).replace('_','.'));
        			String channelSuffix=stationFilenames[numStation][index].substring(i2+1,i3);
        			this.gIP[numFile].channel=-1;
        			for (int j=0;j<this.channelSuffixes.length;j++) if (channelSuffix.equals(this.channelSuffixes[j])) {
        				this.gIP[numFile].channel=j;
        				break;
        			}
        			if (this.gIP[numFile].channel<0) {
        				String msg="invalid file format (channel suffix not recognized) - '"+stationFilenames[numStation][index]+"', should be '<timestamp-seconds>_<timestamp-microseconds>-<channel-number>.tiff'";
        				msg+="\nThis channel suffix is "+channelSuffix+", available channel suffixes are:\n";
        				for (int j=0;j<this.channelSuffixes.length;j++) msg+=this.channelSuffixes[j]+", ";
        				IJ.showMessage("Error",msg);
        				throw new IllegalArgumentException (msg);
        			}
        			numFile++;
        		}
        	}
// Create parameters array
        	initPars (this.gIP.length,parameterDescriptions.length);
        	if (this.debugLevel>1) System.out.println("setupDistortionCalibrationData(): Resetting this.gIS");
        	this.gIS=null; // so it will be initialized in readAllGrids()
        	readAllGrids(
        			patternParameters,
        			laserPointers, // prepare grid parameters for LMA
            		true); // boolean keep_images make it configurable parameter?

        	// no orientation
        }

// from data organized as image sets
        public void setupDirDistortionCalibrationData (
        		String [][]                                            stationFilenames, // per-station List of image set directories
        		String []                                              source_dirs,      // directories of the source files per station
        		MultipleExtensionsFileFilter gridFilter,
        		MultipleExtensionsFileFilter sourceFilter,
        		PatternParameters                                      patternParameters,
        		EyesisCameraParameters                                 eyesisCameraParameters,
        		LaserPointer                                           laserPointers, // as a backup if data is not available in the file
        		Goniometer.GoniometerParameters goniometerParameters,
        		boolean                                                read_grids
        		) {
        	class DirTs{
        		int       station;
        		String [] paths;
        		String [] spaths; // can be null
//        		String    dir;
        		double    ts;
        		int getStation() {return station;}
//        		String getDir() {return dir;}
        		double getTs() {return ts;}
        		String [] getPaths() {return paths;}
        		String [] getSourcePaths() {return spaths;} // may not be null
        		DirTs(int station,
        			  String dir,  // grid image set directory that contains channel files (may be different timestamps)
        			  String sdir, // source super directory that contains image set directories with files
        			  int num_chn,
        			  MultipleExtensionsFileFilter gridFilter,
        			  MultipleExtensionsFileFilter sourceFilter)
        		{
        			this.station = station;
//        			this.dir = dir;
        			int dot_index = dir.lastIndexOf("_");
        			String digits = "0123456789";
        			int ts_start = dot_index -1;
        			while ((ts_start >=0) && (digits.indexOf(dir.charAt(ts_start)) >= 0)) ts_start--;
        			this.ts = Double.parseDouble(dir.substring(ts_start + 1).replace('_','.'));
        			String [] files = (new File(dir)).list(gridFilter); // are these full files?
        			paths = new String[num_chn];
        			for (String path:files) {
        				int last_dash = path.lastIndexOf('-');
        				int last =      path.lastIndexOf('_');
        				if (last_dash >last) last = last_dash;
        				int last_dot = path.lastIndexOf('.');
        				if (last_dot < 0) {
        					last_dot = path.length();
        				}
        				int chn = Integer.parseInt(path.substring(last+1, last_dot));
        				paths[chn] = (new File(dir,path)).getPath();
        				//grid-elphelimg_1559195695_507621_4.tiff
        			}
    				spaths = new String[num_chn];
        			if (sdir != null) {
        				// construct source image set directory name
        				String set_name = (new File(dir)).getName();
        				File set_dir = new File(sdir, set_name );
        				String [] sfiles = set_dir.list(sourceFilter);
        				if (sfiles == null) {
        					System.out.println("sfiles == null");
        				}
        				for (String spath:sfiles) {
        					int last_dash = spath.lastIndexOf('-');
        					int last =      spath.lastIndexOf('_');
        					if (last_dash >last) last = last_dash;
        					int last_dot = spath.lastIndexOf('.');
        					if (last_dot < 0) {
        						last_dot = spath.length();
        					}
        					int chn = Integer.parseInt(spath.substring(last+1, last_dot));
        					spaths[chn] = (new File(set_dir,spath)).getPath();
        				}
        			}
        		}
        	}
    		this.goniometerParameters =  goniometerParameters;
        	boolean ignore_LWIR_pointers = true; // skip LWIR absolute marks, use them later
        	int max_lwir_width = 1023;  // use LWIR class
        	setupIndices();
        	this.eyesisCameraParameters=eyesisCameraParameters;
        	int numSubCameras=(eyesisCameraParameters==null)?1:eyesisCameraParameters.eyesisSubCameras[0].length;
        	this.numSubCameras=numSubCameras;
        	this.eyesisCameraParameters.numStations=stationFilenames.length;
//        	int numFiles=0;

//        	DirTs [][] dirTs = new DirTs [stationFilenames.length][];
    		ArrayList<DirTs> dirTsList = new ArrayList<DirTs>();
        	for (int numStation=0;numStation<stationFilenames.length;numStation++){
        		for (int is = 0; is < stationFilenames[numStation].length; is++) {
        			dirTsList.add(new DirTs(
        					numStation,
        					stationFilenames[numStation][is], // 	String dir,
        					source_dirs[numStation],
        					numSubCameras, // int num_chn,
        					gridFilter,
        					sourceFilter));
        		}
        	}
    		// sort list
    		Collections.sort(dirTsList, new Comparator<DirTs>() {
    		    @Override
    		    public int compare(DirTs lhs, DirTs rhs) {
    		        return rhs.ts > lhs.ts ? -1 : (rhs.ts < lhs.ts) ? 1 : 0;
    		    }
    		});
    		int numFiles = 0;
    		for (DirTs dt:dirTsList) {
    			String [] paths = dt.getPaths();
    			for (String p:paths) if (p!=null) numFiles++;
    		}
    		this.gIS=new GridImageSet[dirTsList.size()];
        	this.gIP=new GridImageParameters[numFiles];
        	int numFile=0;
    		Opener opener=new Opener();
    		JP46_Reader_camera jp4_reader= new JP46_Reader_camera(false);
    		ImagePlus imp_grid=null;
    		ImageStack stack;
        	boolean disableNoFlatfield=false;  // true only for processing transitional images - mixture of ff/ no-ff
    		int numOfGridNodes=        0;
    		int numOfGridNodes_extra=  0;
        	for (int nis = 0; nis<this.gIS.length; nis++) {
        		DirTs dt = dirTsList.get(nis);
        		this.gIS[nis]=new GridImageSet();
        		this.gIS[nis].timeStamp= dt.getTs();
        		this.gIS[nis].imageSet=new GridImageParameters [numSubCameras];
    			this.gIS[nis].setStationNumber(dt.getStation());
    			float [][][] set_pixels = new float [numSubCameras][][];
    			int []       set_widths = new int [numSubCameras];
        		String [] paths =  dt.getPaths();
        		String [] spaths = dt.getSourcePaths();
        		int with_pointers = -1;
        		for (int nc = 0; nc < paths.length; nc++) {
        			String p = paths[nc];
        			boolean first_in_set = true;
        			if (p != null) {
        				this.gIP[numFile] =              new GridImageParameters(numFile);
        				this.gIP[numFile].path =         p;
        				this.gIP[numFile].source_path =  spaths[nc];
        				this.gIP[numFile].setStationNumber(dt.getStation());
        				this.gIP[numFile].timestamp =    dt.getTs();
        				this.gIP[numFile].channel =      nc;
        				this.gIP[numFile].setNumber =    nis;
        				this.gIP[numFile].gridImageSet = this.gIS[nis];
            			this.gIS[nis].imageSet[nc]=this.gIP[numFile];

        				//numFile
        				if (first_in_set || read_grids) {
        					if (read_grids) {
        						if (this.updateStatus) IJ.showStatus("Reading grid file "+(numFile+1)+" (of "+(numFiles)+"): "+this.gIP[numFile].path);
        						if (this.debugLevel>-1) System.out.print(numFile+" ("+this.gIP[numFile].getStationNumber()+
        								":"+this.gIP[numFile].setNumber+":"+this.gIP[numFile].channel+"): "+this.gIP[numFile].path);
        					}
                			imp_grid=opener.openImage("", this.gIP[numFile].path);  // or (path+filenames[nFile])
                			if (imp_grid==null) {
                				String msg="Failed to read grid file "+this.gIP[numFile].path;
                				IJ.showMessage("Error",msg);
                				throw new IllegalArgumentException (msg);
                			}
                	        // TODO: here - need to decode properties
                			jp4_reader.decodeProperiesFromInfo(imp_grid);
                			this.gIP[numFile].woi = new Rectangle(
                					getImagePlusProperty(imp_grid,"WOI_LEFT",0),
                					getImagePlusProperty(imp_grid,"WOI_TOP",0),
                					getImagePlusProperty(imp_grid,"WOI_WIDTH",  eyesisCameraParameters.getSensorWidth(nc)),
                					getImagePlusProperty(imp_grid,"WOI_HEIGHT", eyesisCameraParameters.getSensorHeight(nc)));
//                			boolean woi_compensated = getImagePlusProperty(imp_grid,"WOI_COMPENSATED",false);
                			boolean woi_compensated = (this.gIP[numFile].woi.x == 0) && (this.gIP[numFile].woi.y == 0);

                			this.gIP[numFile].gridImage = imp_grid; // Save all images?
                    		this.gIP[numFile].laserPixelCoordinates = MatchSimulatedPattern.getPointersXYUV(imp_grid, laserPointers);
                    		this.gIP[numFile].motors =                getMotorPositions(imp_grid, this.numMotors);
                    		this.gIS[nis].motors=                     this.gIP[numFile].motors.clone();
                    		this.gIP[numFile].matchedPointers =       getUsedPonters(imp_grid);
                    		if (this.gIP[numFile].matchedPointers > 0) {
                    			// Not using LWIR pointers here!
                    			if (!ignore_LWIR_pointers || (getImagePlusProperty(imp_grid,"WOI_TOP",0) > max_lwir_width)) {
                    				with_pointers = numFile;
                    			}
                    		}
                    		double [] saturations=new double [4];
                    		for (int i=0;i<saturations.length;i++) {
                    			saturations[i]=Double.NaN;
                    			if (imp_grid.getProperty("saturation_" + i) !=null) saturations[i]=Double.parseDouble((String) imp_grid.getProperty("saturation_" + i));
                    		}
                    		if (!Double.isNaN(saturations[1])) this.gIP[numFile].saturation[0]=saturations[1];
                    		if (!Double.isNaN(saturations[2])) this.gIP[numFile].saturation[2]=saturations[2];
                    		if (!Double.isNaN(saturations[0]) && !Double.isNaN(saturations[3])) {
                    			this.gIP[numFile].saturation[1]=0.5*(saturations[0]+saturations[3]);
                    		} else {
                        		if (!Double.isNaN(saturations[0])) this.gIP[numFile].saturation[1]=saturations[0];
                        		if (!Double.isNaN(saturations[3])) this.gIP[numFile].saturation[1]=saturations[3];
                    		}
                    		if (read_grids) {
                    			stack=imp_grid.getStack();
                    			if ((stack==null) || (stack.getSize()<4)) {
                    				String msg="Expected a 8-slice stack in "+this.gIP[numFile].path;
                    				IJ.showMessage("Error",msg);
                    				throw new IllegalArgumentException (msg);
                    			}
                    			float [][] pixels=new float[stack.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB
                    			for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here
                    			if (!woi_compensated) {
                    				System.out.print(" woi_compensation-a "+numFile+" ");
                    				for (int i = 0; i < pixels[0].length; i++) {
                    					pixels[0][i] += this.gIP[numFile].woi.x;
                    					pixels[1][i] += this.gIP[numFile].woi.y;
                    				}
                    				this.gIP[numFile].woi.width += this.gIP[numFile].woi.x;
                    				this.gIP[numFile].woi.x = 0;
                    				this.gIP[numFile].woi.height += this.gIP[numFile].woi.y;
                    				this.gIP[numFile].woi.y = 0;
                    				woi_compensated = true;
                    			}
                    			set_pixels[nc] = pixels;
                    			set_widths[nc] = imp_grid.getWidth();


                    			int numBadNodes = 0;
                    			if (this.eyesisCameraParameters.badNodeThreshold>0.0){
                    				boolean thisDebug =false;
                    				//                            		thisDebug|=        (fileNumber== 720); // chn 25
                    				numBadNodes=fixBadGridNodes(
                    						pixels,
                    						stack.getWidth(),
                    						this.eyesisCameraParameters.badNodeThreshold,
                    						this.eyesisCameraParameters.maxBadNeighb,
                    						this.debugLevel+(thisDebug?3:0),
                    						thisDebug?("fixBad-"+numFile):null
                    						);
                    			}
                    			this.gIP[numFile].flatFieldAvailable=pixels.length>=8;
                            	if (disableNoFlatfield && !this.gIP[numFile].flatFieldAvailable) this.gIP[numFile].enabled=false; // just to use old mixed data

                            	int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIP[numFile].getUVShiftRot());
                            	int [] sizeSizeExtra=setGridsWithRemap(
                            			numFile,
                                		shiftRotMatrix, // int [][] reMap,
                                		pixels,
                                		patternParameters);
                            	numOfGridNodes+=sizeSizeExtra[0];
                            	numOfGridNodes_extra+=sizeSizeExtra[1];

                				if (this.debugLevel>-1) {
                					if (this.gIP[numFile].pixelsUV != null) {
                						System.out.print(" ["+ this.gIP[numFile].pixelsUV.length+"+"+this.gIP[numFile].pixelsUV_extra.length+"]");
                					} else {
                						System.out.print(" [null]");
                					}
                					if (numBadNodes>0)
                						System.out.print("  -- replaced "+numBadNodes+" bad grid nodes");
                					int [] uvrot=this.gIP[numFile].getUVShiftRot();
                					System.out.println(" shift:rot="+uvrot[0]+"/"+uvrot[1]+":"+uvrot[2]+
                							" enabled="+this.gIP[numFile].enabled+" hintedMatch="+this.gIP[numFile].hintedMatch);

                				}
                            	calcGridPeriod(numFile, true); // may be not absolutely calibrated, use_extra out-of-pattern nodes will be used to filter out reflections
                //System.out.println ("pixelsXY["+fileNumber+"]length="+pixelsXY[fileNumber].length);
                    		} //if (read_grids)
                    		// not reading the grid itself
                    		first_in_set = false;
        				}
        				numFile++;

        			}
        		}
        		if (with_pointers < 0) { // no matching pointers, will try to match selected channel with the pattern
        			int main_channel = 4; // one of the EO channels to match with the pattern
//        			boolean [] sensor_mask = null; // later may be used to limit scope to EO-only
        			int extra_search = 2;
//        			int base_channel = this.gIP[with_pointers].channel;
        			if (this.gIS[nis].imageSet[main_channel] != null) {
        				int imgNum =  this.gIS[nis].imageSet[main_channel].imgNumber;
        				boolean invert_color = (main_channel & 4) == 0; // first 4 - LWIR

						if (this.updateStatus) IJ.showStatus("Matching with the pattern, grid file "+(imgNum+1)+" (of "+(numFiles)+"): "+this.gIP[imgNum].path);
						if (this.debugLevel>-1) System.out.print(imgNum+">("+this.gIP[imgNum].getStationNumber()+
								":"+this.gIP[imgNum].setNumber+":"+this.gIP[imgNum].channel+"): "+this.gIP[imgNum].path);

    					double [] sensor_wh = {
    							this.gIP[imgNum].woi.width +  this.gIP[imgNum].woi.x,
    							this.gIP[imgNum].woi.height + this.gIP[imgNum].woi.y};

        				int [] uv_shift_rot = correlateWithPattern(
        		        		patternParameters,
    							set_widths[main_channel], // 		int        test_width,
    							set_pixels[main_channel], //		float [][] test_pixels,
        		        		invert_color,
        		        		extra_search,
        		        		5.0, // double     sigma,
        		        		sensor_wh, // test set pixels width/height pair to reduce weight near the margins (or null)
        		        		false // true // boolean    bdebug
        		        		);
    					System.out.print(" {"+uv_shift_rot[0]+":"+uv_shift_rot[1]+"]");

    					this.gIS[nis].imageSet[main_channel].setUVShiftRot(uv_shift_rot);

                    	int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIS[nis].imageSet[main_channel].getUVShiftRot());
                    	setGridsWithRemap( // null immediately
                    			imgNum,
                    			shiftRotMatrix, // int [][] reMap,
                    			set_pixels[main_channel],
                    			patternParameters);
                    	calcGridPeriod(imgNum, false);  // centered, can skip _extra

    					if (this.gIP[imgNum].pixelsUV != null) {
    						System.out.println(" ["+ this.gIP[imgNum].pixelsUV.length+"+"+this.gIP[imgNum].pixelsUV_extra.length+"]");
    					} else {
    						System.out.println(" [null]");
    					}

    					with_pointers = imgNum; // no adjust all other channels by this one
        			}
        		}

        		if (with_pointers >= 0) { // set initial grids offset from the grid files in the same image set that do not have absolute calibration
        			boolean [] sensor_mask = null; // later may be used to limit scope to EO-only
        			int extra_search = 1;
        			int base_channel = this.gIP[with_pointers].channel;
        			for (int nc = 0; nc < this.gIS[nis].imageSet.length; nc++) if ((sensor_mask == null) || sensor_mask[nc]) {
        				boolean invert_color = ((base_channel ^ nc) & 4) != 0;
        				if ((this.gIS[nis].imageSet[nc].matchedPointers <= 0) && (nc != base_channel)) { // Later add non-laser conditions
        					int imgNum =  this.gIS[nis].imageSet[nc].imgNumber; // with_pointers - base_channel + nc;
    						if (this.updateStatus) IJ.showStatus("Re-reading grid file "+(imgNum+1)+" (of "+(numFiles)+"): "+this.gIP[imgNum].path);
    						if (this.debugLevel>-1) System.out.print(imgNum+"*("+this.gIP[imgNum].getStationNumber()+
    								":"+this.gIP[imgNum].setNumber+":"+this.gIP[imgNum].channel+"): "+this.gIP[imgNum].path);

        					double [] sensor_wh = {
        							this.gIP[imgNum].woi.width +  this.gIP[imgNum].woi.x,
        							this.gIP[imgNum].woi.height + this.gIP[imgNum].woi.y};

        					int [] uv_shift_rot = correlateGrids(
        							set_widths[base_channel], // int        base_width,
        							set_pixels[base_channel], //		float [][] base_pixels,
        							set_widths[nc], // 		int        test_width,
        							set_pixels[nc], //		float [][] test_pixels,
        							invert_color,
        							extra_search,
        							5.0, // 2.0, // sigma
        							sensor_wh,
        							false); // true);
        					System.out.print(" {"+uv_shift_rot[0]+":"+uv_shift_rot[1]+"->");
                        	int [] combinedUVShiftRot=MatchSimulatedPattern.combineUVShiftRot(
                        			this.gIS[nis].imageSet[base_channel].getUVShiftRot(),
                        			uv_shift_rot);

        					this.gIS[nis].imageSet[nc].setUVShiftRot(combinedUVShiftRot); // uv_shift_rot);
        					System.out.print(combinedUVShiftRot[0]+":"+combinedUVShiftRot[1]+"}");

                        	int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIS[nis].imageSet[nc].getUVShiftRot());
                        	setGridsWithRemap( // null immediately
                        			imgNum,
                        			shiftRotMatrix, // int [][] reMap,
                        			set_pixels[nc],
                        			patternParameters);
                        	calcGridPeriod(imgNum, false); // centered, can skip _extra
        					if (this.gIP[imgNum].pixelsUV != null) {
        						System.out.println(" ["+ this.gIP[imgNum].pixelsUV.length+"+"+this.gIP[imgNum].pixelsUV_extra.length+"]");
        					} else {
        						System.out.println(" [null]");
        					}
                		}
        			}
        		} else {

        		}
        	} // for (int nis = 0; nis<this.gIS.length; nis++)


        	initPars (this.gIP.length,parameterDescriptions.length);
        	// readAllGrids(patternParameters); // prepare grid parameters for LMA

        	// Create parameters array
        	///        	initPars (this.gIP.length,parameterDescriptions.length);
        	///        	if (this.debugLevel>1) System.out.println("setupDistortionCalibrationData(): Resetting this.gIS");
        	///        	this.gIS=null; // so it will be initialized in readAllGrids()
        	///        	readAllGrids(patternParameters); // prepare grid parameters for LMA
        	// no orientation
        	if (read_grids) {
        		if (this.debugLevel>3) {
        			System.out.println("setupDirDistortionCalibrationData(), numFiles="+numFiles);
        			for (int n=0;n<this.gIP.length;n++) {
        				System.out.println(n+": length="+this.gIP[n].pixelsXY.length);
        				System.out.println("pixelsUV[][][0]/pixelsUV[][][1] pixelsXY[][][0]/pixelsXY[][][1]");
        				for (int i=0;i<this.gIP[n].pixelsXY.length;i++){
        					System.out.println(n+":"+i+"  "+
        							this.gIP[n].pixelsUV[i][0]+"/"+
        							this.gIP[n].pixelsUV[1][1]+"  "+
        							IJ.d2s(this.gIP[n].pixelsXY[i][0], 2)+"/"+
        							IJ.d2s(this.gIP[n].pixelsXY[i][1], 2)
        							);
        				}
        			}
        		}
        		if (this.debugLevel>0) {
        			System.out.println("setupDirDistortionCalibrationData(), numFiles="+numFiles+", total number of grid nodes="+numOfGridNodes+", unused nodes "+numOfGridNodes_extra);
        		}
        		/// buildImageSets(this.gIS!=null); // already done
        		 // probably - do not need to verify that this.gIS is null - should do that anyway. UPDATE: no, now reading config file creates gIS
        	}

        }

        public boolean initialSetLwirFromEO( //
        		int               num_set,
        		boolean           invert_unmarked_grid,
        		int               extra_search,         // 2
        		double            sigma,                 //5.0
        		PatternParameters patternParameters,
        		boolean           bdebug
        		) {
        	// see if there is any LWIR in the system is sensor and throw if there is none
        	if (!hasSmallSensors()) {
        		String msg="This system does not have any LWIR or other dependent sub-cameras";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if ((this.gIS[num_set] == null) || (this.gIS[num_set].imageSet == null)) {
        		return false;
        	}
        	// See if any of the LWIR subcameras has a mark (absolute grid) in this set
        	int lwir_mark = -1;
        	for (int ns = 0; ns <  this.gIS[num_set].imageSet.length; ns++) {
        		if (this.gIS[num_set].imageSet[ns]!=null) {
        			int imgNum = this.gIS[num_set].imageSet[ns].imgNumber;
        			if (isSmallSensor(imgNum) && (this.gIP[imgNum].matchedPointers > 0)) {
        				lwir_mark = ns;
        				break;
        			}
        		}
        	}
        	int master_sub = getEo0(); //  lowest number EO channel, use as a reference
        	if (lwir_mark >=0) { // some LWIR grid already has mark
        		master_sub = lwir_mark;
        		if (bdebug) {
        			System.out.println ("Aligning to marked LWIR grid image rather than to EO one");
        		}
//        		return false;
        	}
        	if (this.gIS[num_set].imageSet[master_sub] == null) {
        		return false; // master EO is not available // may try to search other channels
        	}
        	int imgMaster = this.gIS[num_set].imageSet[master_sub].imgNumber;
        	// get EO image and pixels
        	ImagePlus imp_master_grid = null;
    		if (this.gIP[imgMaster].gridImage!=null){ // use in-memory grid images instead of the files
    			imp_master_grid=this.gIP[imgMaster].gridImage;
    		} else if (this.gIP[imgMaster].path != null) {
    			imp_master_grid=(new Opener()).openImage("", this.gIP[imgMaster].path);
    			if (imp_master_grid==null) {
    				String msg="Failed to read grid file "+this.gIP[imgMaster].path;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
    			(new JP46_Reader_camera()).decodeProperiesFromInfo(imp_master_grid);
    		} else {
    			System.out.println("EO grid is not in memory, file path is not specified");
    			return false;
    		}
			ImageStack stack_master=imp_master_grid.getStack();
			if ((stack_master==null) || (stack_master.getSize()<4)) {
				String msg="Expected a 8-slice stack";
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
			}
//			boolean woi_compensated_master= getImagePlusProperty(imp_master_grid,"WOI_COMPENSATED",false);
			boolean woi_compensated_master = (getImagePlusProperty(imp_master_grid,"WOI_TOP",0) == 0) &&
					(getImagePlusProperty(imp_master_grid,"WOI_LEFT",0) == 0);


			float [][] pixels_master =new float[stack_master.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB
        	for (int i=0;i<pixels_master.length;i++) pixels_master[i]= (float[]) stack_master.getPixels(i+1); // pixel X : negative - no grid here
			if (!woi_compensated_master) {
				for (int i = 0; i < pixels_master[0].length; i++) {
					pixels_master[0][i] += this.gIP[imgMaster].woi.x;
					pixels_master[1][i] += this.gIP[imgMaster].woi.y;
				}
				this.gIP[imgMaster].woi.width += this.gIP[imgMaster].woi.x;
				this.gIP[imgMaster].woi.x = 0;
				this.gIP[imgMaster].woi.height += this.gIP[imgMaster].woi.y;
				this.gIP[imgMaster].woi.y = 0;
				woi_compensated_master = true;
			}

        	int width_master = imp_master_grid.getWidth();
        	// If master is LWIR - reset it's shift to zero
        	if (lwir_mark >=0) {
    			// reset master sub UV shift
    			int [] zero_uvr = {0,0,0};
	            this.gIS[num_set].imageSet[master_sub].setUVShiftRot(zero_uvr); // uv_shift_rot);
	            int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIS[num_set].imageSet[master_sub].getUVShiftRot());
	            setGridsWithRemap( // null immediately
	            		imgMaster,
	                    shiftRotMatrix, // int [][] reMap,
	                    pixels_master,
	                    patternParameters);
        	}


        	for (int ns = 0; ns <  this.gIS[num_set].imageSet.length; ns++) {
        		if ((this.gIS[num_set].imageSet[ns]!=null) &&(ns != lwir_mark)) {
        			int imgNum = this.gIS[num_set].imageSet[ns].imgNumber;
        			if (isSmallSensor(imgNum)) { // repeat for all target grid images in the image set
        	        	// get target EO image and pixels
        	        	ImagePlus imp_grid = null;
        	    		if (this.gIP[imgNum].gridImage!=null){ // use in-memory grid images instead of the files
        	    			imp_grid=this.gIP[imgNum].gridImage;
        	    		} else if (this.gIP[imgNum].path != null) {
        	    			imp_grid=(new Opener()).openImage("", this.gIP[imgNum].path);
        	    			if (imp_grid==null) {
        	    				String msg="Failed to read grid file "+this.gIP[imgNum].path;
        	    				IJ.showMessage("Error",msg);
        	    				throw new IllegalArgumentException (msg);
        	    			}
        	    			(new JP46_Reader_camera()).decodeProperiesFromInfo(imp_grid);
        	    		} else {
        	    			System.out.println("EO grid is not in memory, file path is not specified");
        	    			return false;
        	    		}
        				ImageStack stack = imp_grid.getStack();
        				if ((stack==null) || (stack.getSize()<4)) {
        					String msg="Expected a 8-slice stack";
        					IJ.showMessage("Error",msg);
        					throw new IllegalArgumentException (msg);
        				}
//        				boolean woi_compensated = getImagePlusProperty(imp_grid,"WOI_COMPENSATED",false);
            			boolean woi_compensated = (this.gIP[imgNum].woi.x == 0) && (this.gIP[imgNum].woi.y == 0);

        				float [][] pixels =new float[stack.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB
        	        	for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here
            			if (!woi_compensated) {
            				System.out.print(" woi_compensation-b "+imgNum+" ");
            				for (int i = 0; i < pixels[0].length; i++) {
            					pixels[0][i] += this.gIP[imgNum].woi.x;
            					pixels[1][i] += this.gIP[imgNum].woi.y;
            				}
            				this.gIP[imgNum].woi.width += this.gIP[imgNum].woi.x;
            				this.gIP[imgNum].woi.x = 0;
            				this.gIP[imgNum].woi.height += this.gIP[imgNum].woi.y;
            				this.gIP[imgNum].woi.y = 0;
            				woi_compensated = true;
            			}

        	        	int width_lwir = imp_grid.getWidth();
    					double [] sensor_wh = {
    							this.gIP[imgNum].woi.width +  this.gIP[imgNum].woi.x,
    							this.gIP[imgNum].woi.height + this.gIP[imgNum].woi.y};
        	            int [] uv_shift_rot = correlateGrids(
        	                    width_master, // int        base_width,
        	                    pixels_master, //		float [][] base_pixels,
        	                    width_lwir, // 		int        test_width,
        	                    pixels, //		float [][] test_pixels,
        	                    invert_unmarked_grid,
        	                    extra_search,
        	                    sigma,
        	                    sensor_wh,
        	                    false); // true); //bdebug

        	            // combined rotation - first this, next what is applied to EO channel

        	            int [] combinedUVShiftRot=MatchSimulatedPattern.combineUVShiftRot(
        	            		uv_shift_rot,
        	                    this.gIS[num_set].imageSet[master_sub].getUVShiftRot()
        	                    );
        	            if (bdebug) {
        	            	System.out.print(imgNum+": calculated uv_shift_rot= ["+uv_shift_rot[0]+":"+uv_shift_rot[1]+"], ");
        	            	System.out.print(imgNum+": EO uv_shift_rot= ["+this.gIS[num_set].imageSet[master_sub].getUVShiftRot()[0]+
        	            			":"+this.gIS[num_set].imageSet[master_sub].getUVShiftRot()[1]+"], ");
        	            	System.out.println(imgNum+": combined uv_shift_rot= ["+combinedUVShiftRot[0]+":"+combinedUVShiftRot[1]+"]");
        	            }

        	            this.gIS[num_set].imageSet[ns].setUVShiftRot(combinedUVShiftRot); // uv_shift_rot);
        	            int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIS[num_set].imageSet[ns].getUVShiftRot());
        	            setGridsWithRemap( // null immediately
        	                    imgNum,
        	                    shiftRotMatrix, // int [][] reMap,
        	                    pixels,
        	                    patternParameters);
        			}
        		}
        	}
        	return true; // OK
        }


        // provide image set index for the same station that has at least one marked image
        // non_estimated - disregard images with estimated orientation
        public int getMarkedSet(int num_set, boolean non_estimated) {
        	int station = this.gIS[num_set].stationNumber;
        	for (int ns = 0; ns < this.gIS.length; ns++) {
        		if (	(this.gIS[ns].stationNumber == station) &&
        				(!non_estimated || !this.gIS[ns].orientationEstimated)) {
                	for (int n=0;n<this.gIS[ns].imageSet.length;n++){
                		if ((this.gIS[ns].imageSet[n]!=null) && (this.gIS[ns].imageSet[n].matchedPointers > 0)){
                			return ns;
                		}
                	}
        		}

        	}
        	return -1; // none found
        }

        // get XYZ for this set's station from marked grid
        public double [] getXYZFromMarked(int num_set, boolean non_estimated) {
        	int ns = getMarkedSet(num_set, non_estimated);
        	if (ns < 0) return null;
        	return this.gIS[ns].GXYZ;
        }

        public double [] getXYZ(int num_img) {
        	int ns = gIP[num_img].getSetNumber();
        	return this.gIS[ns].GXYZ;
        }


        // suggest set grid offset by comparing with known (by mark) set.
        // Wrong Grid UV should cause parallel shift - same Z, different XY
        public int [] suggestOffset (
        		int num_img,
        		boolean non_estimated,
        		boolean even,
        		PatternParameters patternParameters) {
        	int num_set = this.gIP[num_img].getSetNumber();
        	double [] ref_xyz = getXYZFromMarked(num_set, non_estimated);
        	if (ref_xyz == null) {
        		System.out.println("Error: Could not find reference goniometer XYZ for set "+num_set);
        		return null;
        	}
        	double [] diff_xyz = this.gIS[num_set].GXYZ.clone();
        	for (int i = 0; i < diff_xyz.length; i++) diff_xyz[i]-=ref_xyz[i];
        	return suggestOffset (
        			num_img,
        			diff_xyz, // z is not used, may ne just[2]
        			even,
        			patternParameters);
        }

        public int [] suggestOffset (
        		int num_img,
        		double [] diff_xyz, // This XYZ minus reference XYZ  z is not used, may be just[2]
        		boolean even,
        		PatternParameters patternParameters) {
        	int num_set = this.gIP[num_img].setNumber;
        	int station = this.gIS[num_set].stationNumber;
        	int [][] pixelsUV =  this.gIP[num_img].pixelsUV ; // null; // for each image, each grid node - a pair of {gridU, gridV}
        	if ((pixelsUV == null) || ((pixelsUV.length <3 ))) {
        		System.out.println("No/too few pixelsUV data for image "+num_img);
        		return null;
        	}
        	double [][][] data =new double [pixelsUV.length][3][];
        	for (int i=0; i < pixelsUV.length; i++){
        		data[i][0]=new double[2];
        		data[i][1]=new double[2];
        		data[i][2]=new double[1];
        		double [] xyzm = patternParameters.getXYZM(
        				pixelsUV[i][0],
        				pixelsUV[i][1],
        				false, // boolean verbose,
        				station); // int station)
        		data[i][0][0]=xyzm[0];// pixelsXY[i][0];
        		data[i][0][1]=xyzm[1];// pixelsXY[i][1];
        		data[i][1][0]=pixelsUV[i][0];
        		data[i][1][1]=pixelsUV[i][1];
        		data[i][2][0]=xyzm[3];// mask
        	}
        	double [][] coeff=new PolynomialApproximation(this.debugLevel).quadraticApproximation(data, true); // force linear
        	double [] dUV = {
        			-(coeff[0][0]* diff_xyz[0] + coeff[0][1]* diff_xyz[1]),
        			-(coeff[1][0]* diff_xyz[0] + coeff[1][1]* diff_xyz[1])};
        	int [] idUV = {(int) Math.round(dUV[0]), (int) Math.round(dUV[1]), 0}; // 0 - no rot
        	int parity = (idUV[0]+idUV[1] + (even?0:1)) & 1;
        	double [] UV_err = {dUV[0]-idUV[0], dUV[1]-idUV[1]};
        	if (parity !=0) {
        		if (UV_err[1] > UV_err[0]) {
        			if (UV_err[1] > -UV_err[0]) idUV[1]++;
        			else             			idUV[0]--;
        		} else {
        			if (UV_err[1] > -UV_err[0]) idUV[0]++;
        			else	                    idUV[1]--;
        		}
        		UV_err[0] = dUV[0] - idUV[0];
        		UV_err[1] = dUV[1] - idUV[1];
        	}
        	System.out.println(String.format("Errors U/V = %.3f:%.3f",UV_err[0],UV_err[1]));
        	return idUV;
        }


        public int [] offsetGrid(
        		int img_num,
        		int [] uv_shift_rot,
        		PatternParameters patternParameters) {
        	ImagePlus imp_grid = null;
        	boolean woi_compensated = true;
    		if (this.gIP[img_num].gridImage!=null){ // use in-memory grid images instead of the files
    			int numGridImg=img_num;
    			if (numGridImg>=this.gIP.length) numGridImg=this.gIP.length-1;
    			imp_grid=this.gIP[numGridImg].gridImage;
    		}else if (this.gIP[img_num].path != null) {
    			imp_grid=(new Opener()).openImage("", this.gIP[img_num].path);
    			if (imp_grid==null) {
    				String msg="Failed to read grid file "+this.gIP[img_num].path;
    				IJ.showMessage("Error",msg);
    				throw new IllegalArgumentException (msg);
    			}
    			(new JP46_Reader_camera()).decodeProperiesFromInfo(imp_grid);
    			woi_compensated = (getImagePlusProperty(imp_grid,"WOI_TOP",0) == 0) &&
    					(getImagePlusProperty(imp_grid,"WOI_LEFT",0) == 0);
//    			woi_compensated = getImagePlusProperty(imp_grid,"WOI_COMPENSATED",false);
    		} else {
    			System.out.println("Grid is not in memory, file path is not specified");
    			return null;
    		}
			ImageStack stack=imp_grid.getStack();
			if ((stack==null) || (stack.getSize()<4)) {
				String msg="Expected a 8-slice stack";
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
			}
			float [][] pixels=new float[stack.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB

        	for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here

			if (!woi_compensated) { // in memory - always compensated
				System.out.print(" woi_compensation-c "+img_num+" ");
				for (int i = 0; i < pixels[0].length; i++) {
					pixels[0][i] += this.gIP[img_num].woi.x;
					pixels[1][i] += this.gIP[img_num].woi.y;
				}
				this.gIP[img_num].woi.width += this.gIP[img_num].woi.x;
				this.gIP[img_num].woi.x = 0;
				this.gIP[img_num].woi.height += this.gIP[img_num].woi.y;
				this.gIP[img_num].woi.y = 0;
				woi_compensated = true;
			}

        	int [] combinedUVShiftRot=MatchSimulatedPattern.combineUVShiftRot(
                    this.gIP[img_num].getUVShiftRot(),
                    uv_shift_rot);
            this.gIP[img_num].setUVShiftRot(combinedUVShiftRot); // uv_shift_rot);
            int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIP[img_num].getUVShiftRot());
            setGridsWithRemap( // null immediately
            		img_num,
                    shiftRotMatrix, // int [][] reMap,
                    pixels,
                    patternParameters);
            calcGridPeriod(img_num, false); // centered, can skip _extra
            return this.gIP[img_num].getUVShiftRot();
        }


        public static int getImagePlusProperty(ImagePlus imp, String name, int dflt) {
           	try {
           		dflt = Integer.parseInt((String) (imp.getProperty(name)));
           	} catch (Exception e) {

           	}
           	return dflt;
        }

        public static double getImagePlusProperty(ImagePlus imp, String name, double dflt) {
           	try {
           		dflt = Double.parseDouble((String) (imp.getProperty(name)));
           	} catch (Exception e) {

           	}
           	return dflt;
        }

        public static boolean getImagePlusProperty(ImagePlus imp, String name, boolean dflt) {
           	try {
           		dflt = Boolean.parseBoolean((String) (imp.getProperty(name)));
           	} catch (Exception e) {

           	}
           	return dflt;
        }

        public static String getImagePlusProperty(ImagePlus imp, String name, String dflt) {
        	Object obj = imp.getProperty(name);
        	if (obj != null) {
        		dflt =  (String) obj;
        	}
           	return dflt;
        }


        public DistortionCalibrationData (
        		EyesisCameraParameters eyesisCameraParameters,
        		Goniometer.GoniometerParameters goniometerParameters
        		) {
        	setupIndices();
        	this.goniometerParameters = goniometerParameters;
        	int numSubCameras=(eyesisCameraParameters==null)?1:eyesisCameraParameters.eyesisSubCameras[0].length;

        	this.numSubCameras=numSubCameras;
        	this.eyesisCameraParameters=eyesisCameraParameters;
        }


        public DistortionCalibrationData (
        		ImagePlus []           images, // images in the memory
        		PatternParameters      patternParameters,
        		EyesisCameraParameters eyesisCameraParameters,
        		LaserPointer           laserPointers,
        		Goniometer.GoniometerParameters goniometerParameters
        		) {
        	setupIndices();
        	this.goniometerParameters = goniometerParameters;
        	int numSubCameras=(eyesisCameraParameters==null)?1:eyesisCameraParameters.eyesisSubCameras[0].length;
        	this.numSubCameras=numSubCameras;
        	this.eyesisCameraParameters=eyesisCameraParameters;
        	setImages(
        			images,
        			patternParameters,
        			laserPointers);
        }

        public int get_gIS_index(int numImg){
        	if (this.gIS==null) return -1;
        	for (int i=0;i<this.gIS.length;i++)
        		if (this.gIS[i].imageSet!=null)
        			for (int j=0;j<this.gIS[i].imageSet.length;j++)
        				if ((this.gIS[i].imageSet[j]!=null) &&(this.gIS[i].imageSet[j].imgNumber==numImg)) return i;
        	return -1;

        }

        public void listCameraParameters(boolean xcam){
        	int numSubCameras=getNumSubCameras();
        	if (this.gIP!=null) {
        		int maxChn=0;
        		for (int i=0;i<this.gIP.length;i++) if ((this.gIP[i]!=null) && (this.gIP[i].channel>maxChn)){
        			maxChn=this.gIP[i].channel;
        		}
        		numSubCameras=maxChn+1;
        	}

        	if (xcam && (numSubCameras == 4)) {
//           	if (xcam) {
        		listCameraParametersXcam();
        	} else {
        		listCameraParameters();
        	}
        }
        public void listCameraParameters(){
            int numSubCameras=getNumSubCameras();
            if (this.gIP!=null) {
            	int maxChn=0;
            	for (int i=0;i<this.gIP.length;i++) if ((this.gIP[i]!=null) && (this.gIP[i].channel>maxChn)){
            		maxChn=this.gIP[i].channel;
            	}
            	numSubCameras=maxChn+1;
            }
        	String header="Name\tUnits";
        	StringBuffer sb = new StringBuffer();
        	for (int i=0;i<numSubCameras;i++) header+="\t"+i;
        	for (int stationNumber=0;stationNumber<this.eyesisCameraParameters.numStations;stationNumber++){
        		if (this.eyesisCameraParameters.numStations>1){
        			sb.append("Station "+stationNumber+" W="+(100*this.eyesisCameraParameters.stationWeight[stationNumber])+"%");  for (int i=-1;i<numSubCameras;i++) sb.append("\t===");  sb.append("\n");
        		}

        		int [] lensDistortionModels=new int [numSubCameras];
        		for (int i=0;i<numSubCameras;i++) lensDistortionModels[i]=eyesisCameraParameters.getLensDistortionModel(stationNumber,i);
        		sb.append("Lens Distortion Model\t");
        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+lensDistortionModels[i]);
        		sb.append("\n");

        		sb.append("Sensor width\t");
        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+eyesisCameraParameters.getSensorWidth(i));
        		sb.append("\n");

        		sb.append("Sensor height\t");
        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+eyesisCameraParameters.getSensorHeight(i));
        		sb.append("\n");


        		double [][] cameraPars=new double [numSubCameras][];

            	for (int i=0;i<numSubCameras;i++) cameraPars[i]=eyesisCameraParameters.getParametersVector(stationNumber,i);
            	// parameters same order as in this
            	for (int n=0;n<cameraPars[0].length;n++) if (isSubcameraParameter(n) && isIntrinsicParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(cameraPars[i][n],3));
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");

            	for (int n=0;n<cameraPars[0].length;n++) if (isSubcameraParameter(n) && !isIntrinsicParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(cameraPars[i][n],3));
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	for (int n=0;n<cameraPars[0].length;n++) if (
            			!isSubcameraParameter(n)&&
            			!isLocationParameter(n)&&
            			!isOrientationParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		sb.append("\t"+IJ.d2s(cameraPars[0][n],3));
            		for (int i=1;i<numSubCameras;i++) sb.append("\t---");
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	for (int n=0;n<cameraPars[0].length;n++) if (isLocationParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		sb.append("\t"+IJ.d2s(cameraPars[0][n],3));
            		for (int i=1;i<numSubCameras;i++) sb.append("\t---");
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	for (int n=0;n<cameraPars[0].length;n++) if (isOrientationParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		sb.append("\t"+IJ.d2s(cameraPars[0][n],3));
            		for (int i=1;i<numSubCameras;i++) sb.append("\t---");
            		sb.append("\n");
            	}
            }
     	    new TextWindow("Camera parameters", header, sb.toString(), 85*(numSubCameras+3),600);
         }

        public void listCameraParametersXcam(){ // getNumSubCameras() should be 4!
        	double rollDegPerTurn =  -0.45/33.5*180/Math.PI; //  -0.769644799429464 deg/turn, CW screw increases roll, degrees per 1 screw turn
        	double headDegPerTurn = 0.45/34.5*180/Math.PI; //  0.7473362545184652  deg/turn, both screws CW decreases heading (degree/turn)
        	double elevDegPerTurn = 0.45/14*180/Math.PI; //  1.8416500557776463 deg/turn, top CW, bottom CCW decreases elevation (degree/turn)


            int numSubCameras=getNumSubCameras();
            if (this.gIP!=null) {
            	int maxChn=0;
            	for (int i=0;i<this.gIP.length;i++) if ((this.gIP[i]!=null) && (this.gIP[i].channel>maxChn)){
            		maxChn=this.gIP[i].channel;
            	}
            	numSubCameras=maxChn+1;
            }
        	String header="Name\tUnits";
        	StringBuffer sb = new StringBuffer();
        	for (int i=0;i<numSubCameras;i++) header+="\t"+i;
        	for (int stationNumber=0;stationNumber<this.eyesisCameraParameters.numStations;stationNumber++){
        		if (this.eyesisCameraParameters.numStations>1){
        			sb.append("Station "+stationNumber+" W="+(100*this.eyesisCameraParameters.stationWeight[stationNumber])+"%");  for (int i=-1;i<numSubCameras;i++) sb.append("\t===");  sb.append("\n");
        		}

        		int [] lensDistortionModels=new int [numSubCameras];
        		for (int i=0;i<numSubCameras;i++) lensDistortionModels[i]=eyesisCameraParameters.getLensDistortionModel(stationNumber,i);
//        		sb.append("Lens Distortion Model\t");
//        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+lensDistortionModels[i]);
//        		sb.append("\n");

        		sb.append("Sensor width\t");
        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+eyesisCameraParameters.getSensorWidth(i));
        		sb.append("\n");

        		sb.append("Sensor height\t");
        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+eyesisCameraParameters.getSensorHeight(i));
        		sb.append("\n");

        		double [][] cameraPars=new double [numSubCameras][];

            	for (int i=0;i<numSubCameras;i++) cameraPars[i]=eyesisCameraParameters.getParametersVector(stationNumber,i);

            	// calculate average height average right
            	double [] subcamRight =      new double[numSubCameras];
            	double [] subcamHeight =     new double[numSubCameras];
            	double [] subcamCorrRight =  new double[numSubCameras]; // in rotated C.S.
            	double [] subcamCorrHeight = new double[numSubCameras]; // in rotated C.S.
            	double [] subcamHeading =    new double[numSubCameras];
            	double [] subcamElevation =  new double[numSubCameras];
            	double subcamRightCenter =     0.0;
            	double subcamHeightCenter =    0.0;
            	double subcamHeadingCenter =   0.0;
            	double subcamElevationCenter = 0.0;
            	double [] subcamRelRot =        new double[numSubCameras];
            	double [] subcamRelHeading =    new double[numSubCameras];
            	double [] subcamRelElevation =  new double[numSubCameras];

            	double [] rollCorrTurns =       new double[numSubCameras];
            	double [] topCorrTurns =        new double[numSubCameras];
            	double [] botCorrTurns =        new double[numSubCameras];
            	for (int i=0;i<numSubCameras;i++) {
            		subcamRight[i] =     cameraPars[i][getParameterIndexByName("subcamRight")];
            		subcamHeight[i] =    cameraPars[i][getParameterIndexByName("subcamHeight")];
            		subcamHeading[i] =   cameraPars[i][getParameterIndexByName("subcamHeading")];
            		subcamElevation[i] = cameraPars[i][getParameterIndexByName("subcamElevation")];
            		subcamRightCenter +=     subcamRight[i];
            		subcamHeightCenter +=    subcamHeight[i];
            		subcamHeadingCenter +=   subcamHeading[i];
            		subcamElevationCenter += subcamElevation[i];
            	}
            	subcamRightCenter /= numSubCameras;
            	subcamHeightCenter /= numSubCameras;
            	subcamHeadingCenter /= numSubCameras;
            	subcamElevationCenter /= numSubCameras;
            	double [] subcamNominalDirs = {135.0,45.0, -135.0, -45.0};
            	double [] subcamDirsDeg = new double[numSubCameras];
            	double commonRot = 0.0;
            	for (int i=0;i<numSubCameras;i++) {
            		subcamDirsDeg[i]=180.0/Math.PI*Math.atan2(subcamHeight[i]-subcamHeightCenter, subcamRight[i]-subcamRightCenter);
            		commonRot += subcamNominalDirs[i]-subcamDirsDeg[i];
            	}
            	commonRot /= numSubCameras;
            	for (int i=0;i<numSubCameras;i++) {
            		subcamRelRot[i] =       cameraPars[i][getParameterIndexByName("subcamRoll")] - commonRot;
            		subcamRelHeading[i] =   subcamHeading[i] - subcamHeadingCenter;
            		subcamRelElevation[i] =   subcamElevation[i] - subcamElevationCenter;

            		double r = Math.sqrt((subcamHeight[i]-subcamHeightCenter)*(subcamHeight[i]-subcamHeightCenter)+
            				(subcamRight[i]-subcamRightCenter)*(subcamRight[i]-subcamRightCenter));
            		subcamCorrRight[i] = r*Math.cos(Math.PI/180.0*(subcamDirsDeg[i]+commonRot));
            		subcamCorrHeight[i] = r*Math.sin(Math.PI/180.0*(subcamDirsDeg[i]+commonRot));

            		rollCorrTurns[i] = subcamRelRot[i] / rollDegPerTurn;
            		topCorrTurns[i] =  subcamRelHeading[i] / headDegPerTurn + subcamRelElevation[i] / elevDegPerTurn;
            		botCorrTurns[i] =  subcamRelHeading[i] / headDegPerTurn - subcamRelElevation[i] / elevDegPerTurn;

            	}
/*
            	// parameters same order as in this
            	for (int n=0;n<cameraPars[0].length;n++) if (isSubcameraParameter(n) && isIntrinsicParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(cameraPars[i][n],3));
            		sb.append("\n");
            	}
*/
//            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	int flindex =getParameterIndexByName("subcamFocalLength");
        		sb.append(getParameterName(flindex)+"\t"+getParameterUnits(flindex));
        		for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(cameraPars[i][flindex],3));
        		sb.append("\n");

            	sb.append("Camera roll"+"\t"+"degrees"+"\t"+IJ.d2s(commonRot,3));
            	for (int i=1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");

            	sb.append("Camera heading"+"\t"+"degrees"+"\t"+IJ.d2s(subcamHeadingCenter,3));
            	for (int i=1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");

            	sb.append("Camera elevation"+"\t"+"degrees"+"\t"+IJ.d2s(subcamElevationCenter,3));
            	for (int i=1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");

            	sb.append("Rel roll"+"\t"+"degrees");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(subcamRelRot[i],3));
            	sb.append("\n");


            	sb.append("Rel heading"+"\t"+"degrees");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(subcamRelHeading[i],3));
            	sb.append("\n");

            	sb.append("Rel elevation"+"\t"+"degrees");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(subcamRelElevation[i],3));
            	sb.append("\n");

            	sb.append("Corr right"+"\t"+"mm");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(subcamCorrRight[i],3));
            	sb.append("\n");

            	sb.append("Corr height"+"\t"+"mm");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(subcamCorrHeight[i],3));
            	sb.append("\n");

            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	sb.append("Screw roll"+"\t"+"turns CW");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(rollCorrTurns[i],2));
            	sb.append("\n");

            	sb.append("Screw top"+"\t"+"turns CW");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(topCorrTurns[i],2));
            	sb.append("\n");

            	sb.append("Screw bottom"+"\t"+"turns CW");
            	for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(botCorrTurns[i],2));
            	sb.append("\n");

            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");

            	for (int n=0;n<cameraPars[0].length;n++) if (isSubcameraParameter(n) && !isIntrinsicParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		for (int i=0;i<numSubCameras;i++) sb.append("\t"+IJ.d2s(cameraPars[i][n],3));
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	/*
            	for (int n=0;n<cameraPars[0].length;n++) if (
            			!isSubcameraParameter(n)&&
            			!isLocationParameter(n)&&
            			!isOrientationParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		sb.append("\t"+IJ.d2s(cameraPars[0][n],3));
            		for (int i=1;i<numSubCameras;i++) sb.append("\t---");
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	for (int n=0;n<cameraPars[0].length;n++) if (isLocationParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		sb.append("\t"+IJ.d2s(cameraPars[0][n],3));
            		for (int i=1;i<numSubCameras;i++) sb.append("\t---");
            		sb.append("\n");
            	}
            	sb.append("---");  for (int i=-1;i<numSubCameras;i++) sb.append("\t");  sb.append("\n");
            	for (int n=0;n<cameraPars[0].length;n++) if (isOrientationParameter(n)){
            		sb.append(getParameterName(n)+"\t"+getParameterUnits(n));
            		sb.append("\t"+IJ.d2s(cameraPars[0][n],3));
            		for (int i=1;i<numSubCameras;i++) sb.append("\t---");
            		sb.append("\n");
            	}
*/
            }
     	    new TextWindow("Camera parameters", header, sb.toString(), 85*(numSubCameras+3),600);
         }


        public void setImages(
        		ImagePlus [] images,  // images in the memory
        		PatternParameters patternParameters,
        		LaserPointer laserPointers){
        	this.gIP=new GridImageParameters[images.length];
        	for (int i=0;i<images.length;i++){
        		this.gIP[i]=new GridImageParameters(i);
        		this.gIP[i].path=      images[i].getTitle(); // not real path?
        		this.gIP[i].timestamp= getImageTimestamp(images[i]);
        		System.out.println(i+": "+this.gIP[i].path+" - timestamp="+this.gIP[i].timestamp);
        		this.gIP[i].channel=   getImageChannel(images[i]);
            	this.gIP[i].gridImage=images[i]; // free later?

        	}
// Create parameters array
        	initPars (this.gIP.length,parameterDescriptions.length);
        	this.gIS=null; // so it will be created in readAllGrids()
        	readAllGrids(
        			patternParameters, // prepare grid parameters for LMA
        			laserPointers, // prepare grid parameters for LMA
            		true); // boolean keep_images make it configurable parameter?
        	// no orientation
        }

        public void listImageSet(){
        	listImageSet(0, null,null, null,null,null);
        }

        public void listImageSet(
        		int    mode,
        		int [] numPoints, // All arrays may be twice long, then 1 - EO, second - LWIR
        		double [] setRMS,
        		boolean [] hasNaNInSet,
    			int [][] numImgPoints,
    			double [][] rmsPerImg
        		){
        	if ((this.gIS==null) || (this.gIS.length==0)){
        		return;
        	}
        	boolean showXYZ = true;
    		boolean hasLwir = false;
    		if (numPoints!=null) {
        		hasLwir = numPoints.length > this.gIS.length; // twice longer
    		} else if (setRMS != null) {
        		hasLwir = setRMS.length > this.gIS.length; // twice longer
    		} else if (hasNaNInSet != null) {
        		hasLwir = hasNaNInSet.length > this.gIS.length; // twice longer
    		}


        	String header="#\ttimestamp";
        	if (this.eyesisCameraParameters.numStations>1) header+="\tStation";
//        	header+="\tAxial\tTilt\thorPhi\thorPsi\tX\tY\tZ\tMotor2\tMotor3";
        	header+="\tAxial\tTilt\tdTilt\tInter";
        	if (showXYZ) {
            	header+="\tGXY0\tGXY1\tGXY2";
        	}
        	header+="\tMotor2\tMotor3";
        	if (numPoints!=null) {
        		if (hasLwir) {
            		header+="\tNumPointsEO\tNumPointsLWIR";
        		} else {
        			header+="\tNumPoints";
        		}
        	}
        	header+="\tEnabled\tMatched";
        	if (setRMS!=null) {
        		if (hasLwir) {
            		header+="\tRMS-EO\tRMS-LWIR\tWeight";
        		} else {
            		header+="\tRMS\tWeight";
        		}
        	}
        	for (int n=0;n<this.gIS[0].imageSet.length;n++) header+="\t"+n;
    		StringBuffer sb = new StringBuffer();

    		for (int i=0;i<this.gIS.length;i++){
    			double axial_corr_sign=this.gIS[i].goniometerAxial; // correct sign of rotation beyond +/-180 according to motor steps
    			if (this.gIS[i].motors != null) {
    				if (this.gIS[i].motors[1] > 0){
    					if (axial_corr_sign < -90.0) {
    						axial_corr_sign += 360.0;
    					}
    				} else {
    					if (axial_corr_sign > 90.0) {
    						axial_corr_sign -= 360.0;
    					}

    				}
    			}
    			// calculate average tilt for this tilt motor and difference of the current tilt from average
    			double dTilt=Double.NaN;
    			if (!Double.isNaN(this.gIS[i].goniometerTilt) && (this.gIS[i].motors != null)){
    				int i_low,i_high;
    				for (i_low=i-1;i_low>=0;i_low--){
    					if ((this.gIS[i_low].motors != null) && (this.gIS[i_low].motors[2] != this.gIS[i].motors[2])) break;
    				}
    				i_low++;
    				for (i_high=i+1;i_high < this.gIS.length;i_high++){
    					if ((this.gIS[i_high].motors != null) && (this.gIS[i_high].motors[2] != this.gIS[i].motors[2])) break;
    				}
    				int num_avg=0;
    				double sum_avg=0.0;
    				for (int i_avg=i_low;i_avg < i_high; i_avg++){
    					if (!Double.isNaN(this.gIS[i_avg].goniometerTilt)){
    						num_avg++;
    						sum_avg += this.gIS[i_avg].goniometerTilt;
    					}
    				}
    				if (num_avg>0) dTilt = this.gIS[i].goniometerTilt - (sum_avg/num_avg);

    			}
    			double firstInterAxisAngle=Double.NaN;
    			firstInterAxisAngle = this.gIS[i].interAxisAngle;

    			sb.append(i+"\t"+IJ.d2s(this.gIS[i].timeStamp,6));
    			if (this.eyesisCameraParameters.numStations>1)	sb.append("\t"+ this.gIS[i].getStationNumber());
    			sb.append("\t"+(Double.isNaN(this.gIS[i].goniometerAxial)?"---":((this.gIS[i].orientationEstimated?"(":"")+IJ.d2s(axial_corr_sign,3)+(this.gIS[i].orientationEstimated?")":""))));
    			sb.append("\t"+(Double.isNaN(this.gIS[i].goniometerTilt)?"---":((this.gIS[i].orientationEstimated?"(":"")+IJ.d2s(this.gIS[i].goniometerTilt,3)+(this.gIS[i].orientationEstimated?")":""))));

    			sb.append("\t"+(Double.isNaN(dTilt)?"---":IJ.d2s(dTilt,3)));
    			sb.append("\t"+(Double.isNaN(firstInterAxisAngle)?"---":IJ.d2s(firstInterAxisAngle,3)));

            	if (showXYZ) {
                	sb.append(String.format("\t%.1f\t%.1f\t%.1f", this.gIS[i].GXYZ[0], this.gIS[i].GXYZ[1], this.gIS[i].GXYZ[2]));
            	}



    			if (this.gIS[i].motors==null) {
    				sb.append("\t"+"bug"+"\t"+"bug");
    			} else {
    				sb.append("\t"+this.gIS[i].motors[1]+"\t"+this.gIS[i].motors[2]); // null pointer here????
    			}
            	if (numPoints!=null) {
            		if (hasLwir) {
                		sb.append("\t"+numPoints[2 * i + 0]);
                		sb.append("\t"+numPoints[2 * i + 1]);
            		} else {
                		sb.append("\t"+numPoints[i]);
            		}
            	}
            	int numEnImages=0;
            	for (int n=0;n<this.gIS[i].imageSet.length;n++)if (this.gIS[i].imageSet[n]!=null){
            		if (this.gIS[i].imageSet[n].enabled) numEnImages++;
            	}
            	sb.append("\t"+numEnImages);
            	int matchedPointersInSet=0;
            	for (int n=0;n<this.gIS[i].imageSet.length;n++){
            		if (this.gIS[i].imageSet[n]!=null){
            			matchedPointersInSet+=this.gIS[i].imageSet[n].matchedPointers;
            		}
            	}
            	sb.append("\t"+matchedPointersInSet);
            	if (setRMS!=null) {
            		if (hasLwir) {
                		sb.append("\t"+(((hasNaNInSet!=null) && hasNaNInSet[2 * i + 0])?"*":"")+IJ.d2s(setRMS[2 * i + 0],3));
                		sb.append("\t"+(((hasNaNInSet!=null) && hasNaNInSet[2 * i + 1])?"*":"")+IJ.d2s(setRMS[2 * i + 1],3)); //393
                		sb.append("\t"+IJ.d2s(this.gIS[i].setWeight,3));
            		} else {
                		sb.append("\t"+(((hasNaNInSet!=null) && hasNaNInSet[i])?"*":"")+IJ.d2s(setRMS[i],3));
                		sb.append("\t"+IJ.d2s(this.gIS[i].setWeight,3));
            		}
            	}
            	switch (mode) {
            	case 0:
            		for (int n=0;n<this.gIS[i].imageSet.length;n++){
            			sb.append("\t");
            			if (this.gIS[i].imageSet[n]!=null){
            				int numPointers=0; // count number of laser pointers
            				if (this.gIS[i].imageSet[n].laserPixelCoordinates!=null){
            					for (int j=0;j<this.gIS[i].imageSet[n].laserPixelCoordinates.length;j++) {
            						if (this.gIS[i].imageSet[n].laserPixelCoordinates[j]!=null) numPointers++;
            					}
            				}
            				if (!this.gIS[i].imageSet[n].enabled) sb.append("(");
            				sb.append(numPointers+"("+this.gIS[i].imageSet[n].matchedPointers+"):"+this.gIS[i].imageSet[n].hintedMatch +
            						" "+IJ.d2s(this.gIS[i].imageSet[n].getGridPeriod(),1));
            				if (!this.gIS[i].imageSet[n].enabled) sb.append(")");

            			}
            		}
            		break;
            	case 1:
            		for (int n=0;n<this.gIS[i].imageSet.length;n++){
            			sb.append("\t");
            			if (this.gIS[i].imageSet[n]!=null){
            				int [] uvrot = this.gIS[i].imageSet[n].getUVShiftRot();
            				sb.append(uvrot[0]+":"+uvrot[1]+"("+uvrot[2]+")");
            			} else {
            				sb.append("\t---");
            			}
            		}
            		break;
            	case 2:
            		for (int n=0;n<this.gIS[i].imageSet.length;n++){
            			sb.append("\t");
            			if (this.gIS[i].imageSet[n]!=null){
            				sb.append(this.gIS[i].imageSet[n].pixelsXY.length+"+"+this.gIS[i].imageSet[n].pixelsXY_extra.length);
            			} else {
            				sb.append("\t---");
            			}
            		}
            		break;
            	case 3:
            		for (int n=0;n<this.gIS[i].imageSet.length;n++){
            			sb.append("\t");
            			if ((this.gIS[i].imageSet[n]!=null) && (numImgPoints!=null) && (rmsPerImg !=null)){
            				if (Double.isNaN(rmsPerImg[i][n])) {
                				sb.append("NaN ("+numImgPoints[i][n]+")");
            				} else {
                				sb.append(String.format("%.3f (%d)",rmsPerImg[i][n],numImgPoints[i][n]));
            				}
            			} else {
            				sb.append("---");
            			}
            		}
            		break;
            	}
            	sb.append("\n");
    		}
			new TextWindow("Image calibration state (pointers/hinted state)", header, sb.toString(), 1400, 900);
        }


        /**
         * create list of image indices per image set
         * @return array of image indices for each image set
         */
        public int [][] listImages(boolean enabledOnly){
        	int [][] imageSets = new int [this.gIS.length][];
    		for (int i=0;i<this.gIS.length;i++){
    			int setSize=0;
    			for (int n=0;n<this.gIS[i].imageSet.length;n++) if ((this.gIS[i].imageSet[n]!=null) && (this.gIS[i].imageSet[n].imgNumber>=0) && (!enabledOnly || this.gIS[i].imageSet[n].enabled)) setSize++;
    			imageSets[i]=new int [setSize];
    		}
    		for (int i=0;i<this.gIS.length;i++){
    			int index=0;
    			for (int n=0;n<this.gIS[i].imageSet.length;n++) if ((this.gIS[i].imageSet[n]!=null) && (this.gIS[i].imageSet[n].imgNumber>=0) && (!enabledOnly || this.gIS[i].imageSet[n].enabled)) imageSets[i][index++]=this.gIS[i].imageSet[n].imgNumber;
    		}
        	return imageSets;
        }

        /**
         * Filter images (grids) by calibration status with laser pointers and "hinted" from the camera orientation
         * buildImageSets may be needed to be re-ran (if it was ran with all=false)
         * @param resetHinted - if true - reset status of "hinted" calibration to undefined
         * @param minPointers minimal number of laser pointers considered to be enough (usually 2, as mirror/non-mirror is aPriori known
         * @parame minGridPeriod - minimal detected grid period as a fraction of the maximal (filtering reflected grids)
         * @return number of enabled images
         */
        public int [] filterImages(
        		boolean resetHinted,
        		int minPointers,
        		double minGridPeriodFraction,
        		boolean disableNoVignetting,
        		int minGridNodes){
        	int notEnoughNodes=0;
        	int numEnabled=0;
        	int newEnabled=0;
        	int maxPeriod=100;
        	int periodSubdivide=10;
        	int numBins=maxPeriod*periodSubdivide;
        	double [] periodHistogram=new double[numBins];
        	double [] medianGridPeriod=new double [this.eyesisCameraParameters.numStations];
        	double [] maxGridPeriod=new double [this.eyesisCameraParameters.numStations];
        	double [] minGridPeriod=new double [this.eyesisCameraParameters.numStations];
        	// With different sensors different channels wil have different periods
        	int numChannels = getNumSubCameras(); // this.eyesisCameraParameters.getNumChannels();
        	double [][] sw = new double [numChannels][2];
    		for (int i=0;i<this.gIP.length;i++) if (getNumNodes(i, false) >0) {
    			double period = this.gIP[i].getGridPeriod();
    			if (!Double.isNaN(period) && !Double.isInfinite(period) ) {
    				int chn = this.gIP[i].getChannel();
    				sw[chn][0] += getNumNodes(i, false); // weight
    				sw[chn][1] += getNumNodes(i, false) * period; // weight
    			}
    		}

    		double [] avg_periods =  new double [numChannels];
    		double [] ravg_periods = new double [numChannels];
    		this.small_sensors =     new boolean [numChannels];
    		double max_per = 0;
    		for (int i = 0; i < numChannels; i++) {
    			avg_periods[i] =sw[i][1] / sw[i][0];
    			if (max_per < avg_periods[i])  max_per = avg_periods[i];
    		}
    		double [][] sw1 = new double [2][2];
    		for (int i = 0; i < numChannels; i++) {
    			ravg_periods[i] = avg_periods[i] / max_per;
    			small_sensors[i] = ravg_periods[i] < SMALL_FRACTION;
    			sw1[small_sensors[i]?1:0][0] += sw[i][0];
    			sw1[small_sensors[i]?1:0][1] += sw[i][0] * ravg_periods[i];
    		}
    		this.small_period_frac = (sw1[1][0] == 0) ? 0.0 : (sw1[1][1] * sw1[0][0] / (sw1[1][0] * sw1[0][1]));
    		if (small_period_frac > 0.0) {
    			System.out.println(String.format("2 types of sensors are detected, lowres has %5.2f%% resolution",100*small_period_frac));
    			System.out.print("Sensor map: ");
    			for (int i = 0; i < numChannels; i++) {
    				System.out.print(i+":"+(small_sensors[i]?"low-res":"high-res")+", ");
    			}
    			System.out.println();
    		}
        	for (int stationNumber=0;stationNumber<this.eyesisCameraParameters.numStations;stationNumber++){
        		for (int i=0;i<numBins;i++) periodHistogram[i]=0.0;
        		int numSamples=0;
        		for (int i=0;i<this.gIP.length;i++) if (this.gIP[i].getStationNumber()==stationNumber){
        			double period = getEffectivePeriod(i);
//        			if (!Double.isNaN(this.gIP[i].getGridPeriod())) {
        			if (!Double.isNaN(period)) {
//        				int iPeriod=(int) Math.round(this.gIP[i].getGridPeriod()*periodSubdivide);
        				int iPeriod=(int) Math.round(period*periodSubdivide);
        				if (iPeriod>=numBins) iPeriod=numBins-1;
        				else if (iPeriod<0) iPeriod=0; // does not count NaN
        				if (iPeriod>0) {
        					periodHistogram[iPeriod]++;
        					numSamples++;
        				}
        			}
        		}
        		int sumLess=0;
        		medianGridPeriod[stationNumber]=0.0;
        		for (int i=0;i<numBins;i++){
        			sumLess+=periodHistogram[i];
        			if (sumLess>(numSamples/2)) {
        				medianGridPeriod[stationNumber]=(1.0*i)/periodSubdivide;
        				break;
        			}
        		}

        		maxGridPeriod[stationNumber]=0.0;
        		for (int i=0;i<this.gIP.length;i++) if (this.gIP[i].getStationNumber()==stationNumber){
//        			if (this.gIP[i].getGridPeriod()>maxGridPeriod[stationNumber]) {
//        				maxGridPeriod[stationNumber]=this.gIP[i].getGridPeriod();
//        			}
        			if (getEffectivePeriod(i) > maxGridPeriod[stationNumber]) {
        				maxGridPeriod[stationNumber]=getEffectivePeriod(i);
        			}
        		}
        		minGridPeriod[stationNumber]=medianGridPeriod[stationNumber]*minGridPeriodFraction;
            	System.out.print("Station "+stationNumber+ ": maximal grid period="+maxGridPeriod[stationNumber]+" minimal grid period="+minGridPeriod[stationNumber]+" median grid period="+medianGridPeriod[stationNumber]+" numSamples="+numSamples);
            	if (minGridPeriodFraction>0.0) maxGridPeriod[stationNumber]=medianGridPeriod[stationNumber]/minGridPeriodFraction;
            	System.out.println(" new maximal grid period="+maxGridPeriod[stationNumber]);
        	}
        	// set which image set each image belongs
        	int [] gIS_index=new int [this.gIP.length];
        	for (int i=0;i<gIS_index.length;i++)gIS_index[i]=-1;
        	if (this.gIS!=null){
            	for (int i=0;i<this.gIS.length;i++) if (this.gIS[i].imageSet!=null)for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null){
            		gIS_index[this.gIS[i].imageSet[j].imgNumber]=i;
            	}
        	}
        	int numNoVignetting=0;
        	int disabledNoLaser=0;
        	for (int i=0;i<this.gIP.length;i++){
        		int stationNumber=this.gIP[i].getStationNumber();
        		boolean enableNoLaser=this.eyesisCameraParameters.getEnableNoLaser(stationNumber,this.gIP[i].channel);
        		boolean wasEnabled=this.gIP[i].enabled;
        		if (resetHinted) this.gIP[i].hintedMatch=-1; // undefined
///        		if (Double.isNaN(this.gIP[i].getGridPeriod()) ||
///        				((minGridPeriodFraction>0) && ((this.gIP[i].getGridPeriod()<minGridPeriod[stationNumber]) || (this.gIP[i].getGridPeriod()>maxGridPeriod[stationNumber])))){
           		if (Double.isNaN(getEffectivePeriod(i)) ||
           				((minGridPeriodFraction>0) && ((getEffectivePeriod(i)<minGridPeriod[stationNumber]) ||
           						(getEffectivePeriod(i) > maxGridPeriod[stationNumber])))){
        			this.gIP[i].hintedMatch=0; // is it needed?
        			this.gIP[i].enabled=false; // failed against minimal grid period (too far) - probably double reflection in the windows
        		}
        		if (this.gIP[i].hintedMatch==0) {
        			this.gIP[i].enabled=false; // failed against predicted grid
        		} else {
        			if (
        					(this.gIP[i].matchedPointers>=minPointers) ||
        					((this.gIP[i].matchedPointers>0) && (this.gIP[i].hintedMatch>0)) || // orientation and one pointer
        					((this.gIP[i].hintedMatch>1) && enableNoLaser)) { // do not use bottom images w/o matched pointers
        				// before enabling - copy orientation from gIS
        				if (!this.gIP[i].enabled && (gIS_index[i]>=0)){ // FIXME - is it correct to use set index 0?
        					if (!Double.isNaN(this.gIS[gIS_index[i]].goniometerTilt))	setGH(i,this.gIS[gIS_index[i]].goniometerTilt );
        					if (!Double.isNaN(this.gIS[gIS_index[i]].goniometerAxial))	setGA(i,this.gIS[gIS_index[i]].goniometerAxial );
        				}
        				this.gIP[i].enabled=true;
        			} else {
        				this.gIP[i].enabled=false;
        			}
        			if ((this.gIP[i].hintedMatch>1) && !enableNoLaser && (this.gIP[i].matchedPointers==0)){
        				disabledNoLaser++;
        			}
        		}

        		if (disableNoVignetting) {
        			if (this.gIP[i].enabled &!this.gIP[i].flatFieldAvailable) numNoVignetting++;
        			this.gIP[i].enabled &= this.gIP[i].flatFieldAvailable;
        		}
        		if (this.gIP[i].motors==null) { // only disable if any other set has motors
        			boolean hasMotors=false;
                	for (int j=0;j<this.gIP.length;j++){
                		if (this.gIP[j].motors != null) {
                			hasMotors=true;
                			break;
                		}
                	}
                	if (hasMotors) {
                		this.gIP[i].enabled=false; // got some no-motor images made without scanning
                	}
        		}

        		/* Disable no-pointer, new, number of points less than required */
        		if (this.gIP[i].enabled && !wasEnabled && (this.gIP[i].matchedPointers==0) && (this.gIP[i].pixelsXY.length<minGridNodes)){
        			this.gIP[i].enabled=false;
        			notEnoughNodes++;
        		}

            	if (this.gIP[i].enabled) numEnabled++;
            	this.gIP[i].newEnabled=this.gIP[i].enabled&&!wasEnabled;
            	if (this.gIP[i].newEnabled) newEnabled++;
        	}
        	// may need buildImageSets
        	int [] result={numEnabled,newEnabled,numNoVignetting,notEnoughNodes,disabledNoLaser};
        	return result;
        }
// TODO:
        // 1 -  Filter by lasers/hint state
        // 2 - recalculate hinted
        // connect "enabled" to strategies (not done yet)
      //  applyHintedGrids90 - moved to the parent class


        /**
         * Create array of image sets ("panoramas"), sorted by timestamps
         * @return number of sets
         */
        public int buildImageSets(boolean preserveSet){
        	if (this.debugLevel>0) {
        		System.out.println("buildImageSets("+preserveSet+")");
        	}
        	if (!preserveSet){
        		List <Double> timeStampList=new ArrayList<Double>(this.gIP.length);
        		int numChannels=0;
        		for (int i=0;i<this.gIP.length;i++) {
        			if (this.gIP[i].channel>numChannels) numChannels=this.gIP[i].channel;
        			int j=0;
        			Double ts=this.gIP[i].timestamp;
        			if (!timeStampList.contains(ts)){
        				for (;(j<timeStampList.size()) && (ts>timeStampList.get(j));j++);
        				timeStampList.add(j,ts);
        			}
        		}
        		numChannels++;
        		this.gIS=new GridImageSet[timeStampList.size()];
        		for (int i=0;i<this.gIS.length;i++){
        			this.gIS[i]=new GridImageSet();
        			this.gIS[i].timeStamp=timeStampList.get(i);
        			this.gIS[i].imageSet=new GridImageParameters [numChannels];
        			for (int j=0;j<numChannels;j++) this.gIS[i].imageSet[j]=null;

        		}
        		for (int i=0;i<this.gIP.length;i++) {
        			Double ts=this.gIP[i].timestamp;
        			int iIS=timeStampList.indexOf(ts);
        			this.gIS[iIS].setStationNumber(this.gIP[i].getStationNumber());
        			this.gIS[iIS].imageSet[this.gIP[i].channel]=this.gIP[i];
//        			if (this.gIP[i].motors!=null) this.gIS[iIS].motors=this.gIP[i].motors;
        			this.gIP[i].setNumber=iIS;
        			this.gIP[i].gridImageSet=this.gIS[iIS];
        		}
        		// verify that station number is the same for the same timestamp
        		for (int i=0;i<this.gIP.length;i++) {
        			Double ts=this.gIP[i].timestamp;
        			int iIS=timeStampList.indexOf(ts);
        			if (this.gIS[iIS].getStationNumber()!=this.gIP[i].getStationNumber()){
        				String msg="Inconsistent station number for timestamp "+ts+": this.gIS[iIS].getStationNumber()="+this.gIS[iIS].getStationNumber()+
        				" this.gIP[i].getStationNumber()="+this.gIP[i].getStationNumber()+", using "+this.gIS[iIS].getStationNumber();
        				System.out.println(msg);
        				IJ.showMessage("Error:",msg);
        				this.gIP[i].setStationNumber(this.gIS[iIS].getStationNumber());
        			}
        		}
        	}
    		for (int i=0;i<this.gIP.length;i++) {
    			int iIS=this.gIP[i].setNumber;
    			if (this.gIP[i].motors!=null) this.gIS[iIS].motors=this.gIP[i].motors.clone();
    		}
        	return this.gIS.length;
        }

        /**
         * Create array of image sets ("panoramas"), sorted by timestamps
         * @param all // use all images (false - only enabled)
         * @return number of sets
         */

        public int buildImageSetsOld(boolean all){
        	List <Double> timeStampList=new ArrayList<Double>(this.gIP.length);
        	int numChannels=0;
        	for (int i=0;i<this.gIP.length;i++) if (all || this.gIP[i].enabled){
        		if (this.gIP[i].channel>numChannels) numChannels=this.gIP[i].channel;
        		int j=0;
        		Double ts=this.gIP[i].timestamp;
        		if (!timeStampList.contains(ts)){
        			for (;(j<timeStampList.size()) && (ts>timeStampList.get(j));j++);
        			timeStampList.add(j,ts);
        		}
        	}
        	numChannels++;
        	this.gIS=new GridImageSet[timeStampList.size()];
        	for (int i=0;i<this.gIS.length;i++){
        		this.gIS[i]=new GridImageSet();
        		this.gIS[i].timeStamp=timeStampList.get(i);
        		this.gIS[i].imageSet=new GridImageParameters [numChannels];
        		for (int j=0;j<numChannels;j++) this.gIS[i].imageSet[j]=null;

        	}
        	for (int i=0;i<this.gIP.length;i++) if (all || this.gIP[i].enabled){
        		Double ts=this.gIP[i].timestamp;
        		int iIS=timeStampList.indexOf(ts);
        		this.gIS[iIS].imageSet[this.gIP[i].channel]=this.gIP[i];
        		if (this.gIP[i].motors!=null) this.gIS[iIS].motors=this.gIP[i].motors;
        	}
        	return this.gIS.length;
        }

        /**
         * Set goniometer initial orientation from the image with maximal number of laser pointers (make averaging later?)
         * Needed before LMA to have some reasonable initial orientation
         * @param overwriteAll if true, overwrite orientation data even if it is already not NaN, false -skip those that have orientation set
         */
        public void setInitialOrientation(
        		PatternParameters  patternParameters,
        		boolean            overwriteAll) {
			if (this.debugLevel>0) {
				System.out.println("setInitialOrientation("+overwriteAll+"), debugLevel= "+this.debugLevel);
			}

        	for (int i=0; i<this.gIS.length;i++){
        		int stationNumber=this.gIS[i].getStationNumber();
        		int bestRating=-1;
        		int bestChannel=-1;
        		for (int j=0;j<this.gIS[i].imageSet.length;j++) if ((this.gIS[i].imageSet[j]!=null) && this.gIS[i].imageSet[j].enabled){
        			int thisRating=this.gIS[i].imageSet[j].matchedPointers+((this.gIS[i].imageSet[j].hintedMatch>0)?1:0); // rate hintedMatch 2 higher?
        			if (thisRating>bestRating) {
        				bestRating=thisRating;
        				bestChannel=j;
        			}
        		}
        		if (bestRating>0){
        			EyesisSubCameraParameters esp = this.eyesisCameraParameters.eyesisSubCameras[stationNumber][bestChannel];
        			double [] uv_center = getGridUVfromXY(
        	        		esp.px0, // final double px,
        	        		esp.py0, // final double py,
        	        		this.gIS[i].imageSet[bestChannel].getImageNumber(), //  final int fileNumber,
        	        		true); // boolean use_extra)
// find UV of the center of the image getImageNumber
        			if (uv_center == null) {
        				if (this.debugLevel>0) {
        					System.out.println("Center UV = NULL");
        				}
        			} else {
        				if (this.debugLevel>0) {
        					System.out.println("Center UV = "+uv_center[0]+","+uv_center[1]);
        				}
        				double [] patt_xyz = 		patternParameters.getXYZ(
        						uv_center, // double [] uv,
        						false, // boolean verbose,
        						this.gIS[i].getStationNumber()); // int station); // u=0,v=0 - center!

            			if (patt_xyz == null) {
            				if (this.debugLevel>0) {
            					System.out.println("Center UV = NULL");
            				}
            			} else {
            				if (this.debugLevel>0) {
            					System.out.println("Center XYZ = "+patt_xyz[0]+","+patt_xyz[1]+","+patt_xyz[2]);
            				}
// Calculate position relative to the view point on the target
            				double [] aview = {
            						patt_xyz[0]- this.gIS[i].GXYZ[0],
            						-patt_xyz[1]+ (this.gIS[i].GXYZ[1] - this.gIS[i].centerAboveHorizontal),
            						-patt_xyz[2]+ this.gIS[i].GXYZ[2]
            				};
            				Matrix mview = new Matrix(aview,3);
            				double phi = -Math.PI/180.0*this.gIS[i].horAxisErrPhi;
            				double cp = Math.cos(phi);
            				double sp = Math.sin(phi);
            				double [][] aphi = {
            						{ cp, 0.0, sp},
            						{0.0, 1.0, 0.0},
            						{-sp, 0.0, cp}};
            				Matrix mphi = new Matrix(aphi);
            				Matrix mview_gon = mphi.times(mview); // view point on the target from the goniometer
            				double tilt = -Math.atan2(mview_gon.get(1, 0), mview_gon.get(2, 0)); // y pointed up
            				double ct =  Math.cos(tilt);
            				double st = Math.sin(tilt);

            				double [][] atilt = {
            						{1.0, 0.0, 0.0},
            						{0.0,  ct,  st},
            						{0.0, -st,  ct}};
            				Matrix mtilt = new Matrix(atilt);
            				Matrix mview_tilt = mtilt.times(mview_gon); // view point on the target from the tilted goniometer
            				double az = Math.atan2(mview_tilt.get(0, 0), mview_tilt.get(2, 0)); // x pointed right

            				double tilt_deg = tilt/Math.PI*180;
            				double az_deg =   az/Math.PI*180;
            				if (this.debugLevel>0) {
            					System.out.println("Tilt = "+tilt_deg+", az = "+az_deg);
            					System.out.print("");
            				}
            			}
        			}

        			if (overwriteAll || Double.isNaN(this.gIS[i].goniometerAxial)){
 //       				System.out.println("setInitialOrientation("+overwriteAll+"),  Double.isNaN(this.gIS["+i+"].goniometerAxial)="+Double.isNaN(this.gIS[i].goniometerAxial));

        				double subcam_heading = (esp.heading + (esp.cartesian? 0: esp.azimuth));
        				this.gIS[i].goniometerAxial=-subcam_heading;
        				for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) setGA(this.gIS[i].imageSet[j].imgNumber,this.gIS[i].goniometerAxial);
            			this.gIS[i].orientationEstimated=true;
            			if (this.debugLevel>1) {
            				System.out.print(String.format("Setting goniometerAxial for the image set #%4d (%18.6f) to ", i, this.gIS[i].timeStamp));
            				System.out.println(""+this.gIS[i].goniometerAxial+" +++++ orientationEstimated==true +++++");
//            				System.out.println("Setting goniometerAxial for the image set #"+i+" ("+this.gIS[i].timeStamp+") to "+this.gIS[i].goniometerAxial+" +++++ orientationEstimated==true +++++");
            			}
        			}
        			if (overwriteAll || Double.isNaN(this.gIS[i].goniometerTilt )){
//        				System.out.println("setInitialOrientation("+overwriteAll+"),  Double.isNaN(this.gIS["+i+"].goniometerTilt)="+Double.isNaN(this.gIS[i].goniometerTilt));
        				this.gIS[i].goniometerTilt= -esp.theta;
        				for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) setGH(this.gIS[i].imageSet[j].imgNumber,this.gIS[i].goniometerTilt);
            			this.gIS[i].orientationEstimated=true;
            			if (this.debugLevel>1) {
            				System.out.print(String.format("Setting goniometerTilt  for the image set #%4d (%18.6f) to ", i, this.gIS[i].timeStamp));
            				System.out.println(""+this.gIS[i].goniometerTilt+" ===== orientationEstimated==true =====");
            			}
        			}
        		}
        	}
        }
        /**
         * update image set (panorama, set of simultaneous images) goniometer orientation from the image parameters, do after running LMA
         * @param selectedImages boolean array of selected images (in current strategy) or null (all selected)
         */
// TODO: potential problem here if only some images were enabled in the strategy -- FIXED
// TODO: Add other extrinsic parameters here to sets?
        /**
         * Updated version - only flag as orientationEstimated if no enabled images exist in the set or any of the angles is NaN
         * Temporarily duplicate  image parameters from those of the set (should not be needed)
         * selectedImages will not be used
         */
        public void updateSetOrientation(boolean [] selectedImages){ // if selectedImages[] is not null will set orientationEstimated for unselected images
        	if (this.gIS==null){
        		String msg="Image set is not initilaized";
        		System.out.println(msg);
        		IJ.showMessage(msg);
        	}

        	for (int i=0; i<this.gIS.length;i++){
        		this.gIS[i].orientationEstimated=true;
        		if (!Double.isNaN(this.gIS[i].goniometerAxial) && !Double.isNaN(this.gIS[i].goniometerTilt)) {
        			for (int j=0;j<this.gIS[i].imageSet.length;j++) if ((this.gIS[i].imageSet[j]!=null) && this.gIS[i].imageSet[j].enabled){
        				if ((selectedImages==null) || selectedImages[this.gIS[i].imageSet[j].imgNumber]) {
        					this.gIS[i].goniometerAxial-=360.0*Math.floor((this.gIS[i].goniometerAxial+180.0)/360.0);
        					this.gIS[i].orientationEstimated=false;
        					break; // set from the first non-null, enabled image
        				}
        			}
        		}
        		if (!this.gIS[i].orientationEstimated){
        			// now fill that data to all disabled images of the same set (just for listing RMS errors and debugging)
        			for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) { // fill even those that are enabled
        				setGA(this.gIS[i].imageSet[j].imgNumber,this.gIS[i].goniometerAxial );
        				setGH(this.gIS[i].imageSet[j].imgNumber,this.gIS[i].goniometerTilt );
        			}
        		} else {
        			this.gIS[i].goniometerAxial=Double.NaN;
        			this.gIS[i].goniometerTilt= Double.NaN;
        			System.out.println("updateSetOrientation(): imageSet "+i+" orientationEstimated == true");
        		}
        	}
        }

        public void updateSetOrientationOld(boolean [] selectedImages){
        	if (this.gIS==null){
        		String msg="Image set is not initilaized";
        		System.out.println(msg);
        		IJ.showMessage(msg);
        	}
        	for (int i=0; i<this.gIS.length;i++){
        		if (selectedImages==null){ // if all selected - remove orientation if there are no enabled images (i.e. after removeOutliers)
    				this.gIS[i].goniometerAxial=Double.NaN;
    				this.gIS[i].goniometerTilt= Double.NaN;
    				this.gIS[i].orientationEstimated=true;

        		}
        		for (int j=0;j<this.gIS[i].imageSet.length;j++) if ((this.gIS[i].imageSet[j]!=null) && this.gIS[i].imageSet[j].enabled){
        			if ((selectedImages==null) || selectedImages[this.gIS[i].imageSet[j].imgNumber]) {
        				this.gIS[i].goniometerAxial=getGA(this.gIS[i].imageSet[j].imgNumber);  //update - most likely will do nothing (if set has non-NaN)
        				this.gIS[i].goniometerTilt= getGH(this.gIS[i].imageSet[j].imgNumber);
        				this.gIS[i].goniometerAxial-=360.0*Math.floor((this.gIS[i].goniometerAxial+180.0)/360.0);
        				this.gIS[i].orientationEstimated=false;
        				break; // set from the first non-null, enabled image
        			}
        		}
        		// now fill that data to all disabled images of the same set (just for listing RMS errors and debugging)
        		if (!Double.isNaN(this.gIS[i].goniometerAxial) && !Double.isNaN(this.gIS[i].goniometerTilt)){
        			for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) { // fill even those that are enabled
        				setGA(this.gIS[i].imageSet[j].imgNumber,this.gIS[i].goniometerAxial );
        				setGH(this.gIS[i].imageSet[j].imgNumber,this.gIS[i].goniometerTilt );
        			}
        		}
        	}
        }

        public boolean isEstimated(int imgNum){
        	if (this.gIS==null) {
            	String msg="Image sets are not initialized";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (this.gIP[imgNum].gridImageSet!=null) return this.gIP[imgNum].gridImageSet.orientationEstimated;
        	// should not get here
        	System.out.println("FIXME: isEstimated("+imgNum+"): this.gIP["+imgNum+"].gridImageSet==null");
        	for (int i=0;i<this.gIS.length;i++)	if (this.gIS[i].imageSet!=null){
        		for (int j=0;j<this.gIS[i].imageSet.length;j++) if ((this.gIS[i].imageSet[j]!=null) && (this.gIS[i].imageSet[j].imgNumber==imgNum)){
        			return this.gIS[i].orientationEstimated;
        		}
        	}
        	String msg="Image with index "+imgNum+" is not in the image set";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
        }
        public boolean isEstimatedOld(int imgNum){
        	if (this.gIS==null) {
            	String msg="Image sets are not initialized";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	for (int i=0;i<this.gIS.length;i++)	if (this.gIS[i].imageSet!=null){
        		for (int j=0;j<this.gIS[i].imageSet.length;j++) if ((this.gIS[i].imageSet[j]!=null) && (this.gIS[i].imageSet[j].imgNumber==imgNum)){
        			return this.gIS[i].orientationEstimated;
        		}
        	}
        	String msg="Image with index "+imgNum+" is not in the image set";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
        }
        public int getNumberOfEstimated(boolean enabledOnly) {
        	int numEstimated=0;
        	if (this.gIS==null) return 0;
        	for (int i=0;i<this.gIS.length;i++)	if (this.gIS[i].imageSet!=null){
        		for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) {
        			if ((!enabledOnly || this.gIS[i].imageSet[j].enabled) && this.gIS[i].orientationEstimated) numEstimated++;

        		}
        	}
        	return numEstimated;
        }

        public int [] getNumberOfEstimatedPerStation(boolean enabledOnly) {
        	int [] numEstimated=new int [this.eyesisCameraParameters.numStations];
        	for (int i=0;i<numEstimated.length;i++) numEstimated[i]=0;
        	if (this.gIS!=null){
        		for (int i=0;i<this.gIS.length;i++)	if (this.gIS[i].imageSet!=null){
        			for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) {
        				if ((!enabledOnly || this.gIS[i].imageSet[j].enabled) && this.gIS[i].orientationEstimated) numEstimated[this.gIS[i].getStationNumber()]++;
        			}
        		}
        	}
        	return numEstimated;
        }


        public int getNumEnabled(){
        	int num=0;
        	for (int i=0;i<this.gIP.length;i++) if ((this.gIP[i]!=null) && this.gIP[i].enabled) num++;
        	return num;
        }

        public int getNumNewEnabled(){
        	int num=0;
        	for (int i=0;i<this.gIP.length;i++) if ((this.gIP[i]!=null) && this.gIP[i].enabled && this.gIP[i].newEnabled) num++;
        	return num;
        }

        public int [] getNumNewEnabledPerStation(){
        	int [] numEnabled=new int [this.eyesisCameraParameters.numStations];
        	for (int i=0;i<numEnabled.length;i++) numEnabled[i]=0;
        	for (int i=0;i<this.gIP.length;i++) if ((this.gIP[i]!=null) && this.gIP[i].enabled && this.gIP[i].newEnabled) numEnabled[this.gIP[i].getStationNumber()]++;//  OOB 837
        	return numEnabled;
        }

        public int [] getStations(){
        	int [] result = new int [this.gIP.length];
        	for (int i=0;i<result.length;i++) result[i]=(this.gIP[i]!=null)?this.gIP[i].stationNumber:-1;
        	return result;
        }
        public int [] getChannels(){
        	int [] result = new int [this.gIP.length];
        	for (int i=0;i<result.length;i++) result[i]=(this.gIP[i]!=null)?this.gIP[i].channel:-1;
        	return result;
        }
        public int [] getMatchedPointers(){
        	int [] result = new int [this.gIP.length];
        	for (int i=0;i<result.length;i++) result[i]=(this.gIP[i]!=null)?this.gIP[i].matchedPointers:0;
        	return result;
        }
        public int [] getHintedMatch(){
        	int [] result = new int [this.gIP.length];
        	for (int i=0;i<result.length;i++) result[i]=(this.gIP[i]!=null)?this.gIP[i].hintedMatch:-1;
        	return result;
        }

        public boolean [] selectNewEnabled () {
        	boolean [] newEnabled=new boolean [this.gIP.length] ;
        	for (int i=0;i<this.gIP.length;i++) newEnabled[i]= (this.gIP[i]!=null) && this.gIP[i].enabled && this.gIP[i].newEnabled;
        	return newEnabled;
        }

        public boolean [] selectEnabled () {
        	boolean [] enabled=new boolean [this.gIP.length] ;
        	for (int i=0;i<this.gIP.length;i++) enabled[i]= (this.gIP[i]!=null) && this.gIP[i].enabled;
        	return enabled;
        }

        public boolean [] selectEstimated (boolean enabledOnly) {
        	boolean [] estimated=new boolean [getNumImages()];
        	if (this.gIS==null) {
            	String msg="Image sets are not initialized";
        		IJ.showMessage("Error",msg);
//        		throw new IllegalArgumentException (msg);
        		Arrays.fill(estimated, true);
            	return estimated;
        	}

        	for (int i=0;i<estimated.length;i++) estimated[i]=false;
        	for (int i=0;i<this.gIS.length;i++)	if (this.gIS[i].imageSet!=null){
        		for (int j=0;j<this.gIS[i].imageSet.length;j++) if (this.gIS[i].imageSet[j]!=null) {
        			if ((!enabledOnly || this.gIS[i].imageSet[j].enabled) ) estimated[this.gIS[i].imageSet[j].imgNumber]= this.gIS[i].orientationEstimated;

        		}
        	}
        	return estimated;
        }
        public void enableSelected(boolean [] selected){
        	for (int i=0;i<this.gIP.length  ;i++) if (this.gIP[i]!=null){
        		int i1=i;
        		if (i1>=selected.length) i1=selected.length-1;
        		this.gIP[i].enabled = selected[i1];
        	}
        }
        /**
         * Calculate goniometer orientation for one of the "known" images/grids
         * @param imgNum grid image number
         * @return pair of {goniometerHorizontal, goniometerAxial} (in angular degrees)
         */
        public double [] getImagesetTiltAxial(int imgNum){
        	return getImagesetTiltAxial(this.gIP[imgNum].timestamp);
        }
        /**
         * Return pair of {goniometerHorizontal, goniometerAxial} for the specified timestamp
         * updateSetOrientation() should be called after LMA or other updates to camera parameters
         * @param timeStamp - double timestamp identifying imageset (image does not need to be a part of selected grid files)
         * @return null if no images set has the specified timestamp, may contain Double.NaN if the orientation was not set.
         * Now 3-rd term - interAxisAngle - with goniometerTilt it is used for correction of non-pure axial movement of the camera.
         */
        public double [] getImagesetTiltAxial(double timeStamp){
        	int mAxial=1;     // m2
        	int mHorizontal=2;// m3
        	// this is probably already set
        	for (int i=0;i<this.gIS.length;i++){
        		if ((this.gIS[i].imageSet!=null) && (this.gIS[i].imageSet.length>0) && (this.gIS[i].imageSet[0]!=null)) this.gIS[i].setStationNumber(this.gIS[i].imageSet[0].getStationNumber());
            }
        	for (int i=0;i<this.gIS.length;i++)
        		if (this.gIS[i].timeStamp==timeStamp) {
    				int iBest=i;
        			if (Double.isNaN(this.gIS[i].goniometerTilt) || Double.isNaN(this.gIS[i].goniometerAxial)  || Double.isNaN(this.gIS[i].interAxisAngle)) {
// find the closest one (by motors)
        				if (this.gIS[i].motors==null) {
                			if (this.debugLevel>0) System.out.println("getImagesetTiltAxial("+timeStamp+"): No motor data");
                			if (this.debugLevel>0) System.out.println("Looking for closest timestamps in the same station, image set = "+i);
                			int early_set=-1;
                			int late_set=-1;
                			for (int j=0; j<this.gIS.length;j++) {
                    			if (Double.isNaN(this.gIS[j].goniometerTilt) || Double.isNaN(this.gIS[j].goniometerAxial)  || Double.isNaN(this.gIS[j].interAxisAngle)) continue;
                    			if (this.gIS[j].timeStamp > timeStamp){
                    				if ((late_set<0) || (this.gIS[j].timeStamp < this.gIS[late_set].timeStamp)) late_set = j;
                    			} else {
                    				if ((early_set<0) || (this.gIS[j].timeStamp >this.gIS[early_set].timeStamp)) early_set = j;
                    			}
                			}

                			if ((late_set <0) && (early_set<0)) {
                    			if (this.debugLevel>0) System.out.println("Failed to find any known orientation");
                    			return null;
                			}
                			if       (late_set <0) {
                				iBest= early_set;
                			} else if  (early_set <0) {
                				iBest= late_set;
                			} else {
                				// interpolate
            					double axialEarly=this.gIS[early_set].goniometerAxial;
            					double axialLate= this.gIS[late_set].goniometerAxial;
            					axialEarly-=360.0*Math.floor((axialEarly+180.0)/360.0); // convert to range +/-180
            					axialLate-= 360.0*Math.floor((axialLate+ 180.0)/360.0);
            					double axialCenter= 0.5*(axialEarly+axialLate);
            					if (Math.abs(axialEarly-axialLate)>180) {
            						if (axialCenter>0) axialCenter-=180.0;
            						else axialCenter+=180.0;
            					}

            					double interEarly=this.gIS[early_set].interAxisAngle;
            					double interLate= this.gIS[late_set].interAxisAngle;
            					interEarly-=360.0*Math.floor((interEarly+180.0)/360.0);
            					interLate-= 360.0*Math.floor((interLate+ 180.0)/360.0);
            					double interCenter= 0.5*(interEarly+interLate);
            					if (Math.abs(interEarly-interLate)>180) {
            						if (interCenter>0) interCenter-=180.0;
            						else interCenter+=180.0;
            					}

            					double tiltEarly=this.gIS[early_set].goniometerTilt;
            					double tiltLate= this.gIS[late_set].goniometerTilt;
            					tiltEarly-=360.0*Math.floor((tiltEarly+180.0)/360.0);
            					tiltLate-= 360.0*Math.floor((tiltLate+ 180.0)/360.0);
            					double tiltCenter= 0.5*(tiltEarly+tiltLate);
            					if (Math.abs(tiltEarly-tiltLate)>180) {
            						if (tiltCenter>0) tiltCenter-=180.0;
            						else tiltCenter+=180.0;
            					}
            					this.gIS[i].goniometerTilt= tiltCenter;
            					this.gIS[i].goniometerAxial=axialCenter;
            					this.gIS[i].interAxisAngle=interCenter;
                				if (this.debugLevel>2) System.out.println("getImagesetTiltAxial("+timeStamp+"):"+
                						" axialEarly - "+ axialEarly+
                						" axialLate - "+  axialLate+
                						" axialCenter - "+axialCenter+
                						" tiltEarly - "+  tiltEarly+
                						" tiltLate - "+   tiltLate+
                						" tiltCenter - "+ tiltCenter+
                						" interEarly - "+ interEarly+
                						" interLate - "+  interLate+
                						" interCenter - "+interCenter);
                			}
                			if (iBest!=i) {
                			// use closest
                				this.gIS[i].goniometerTilt= this.gIS[iBest].goniometerTilt;
                				this.gIS[i].goniometerAxial=this.gIS[iBest].goniometerAxial;
                				this.gIS[i].interAxisAngle=this.gIS[iBest].interAxisAngle;
                			}
            				this.gIS[i].orientationEstimated=true;
                			double [] result = {
//                					this.gIS[iBest].goniometerTilt,
//                					this.gIS[iBest].goniometerAxial,
///                					this.gIS[iBest].interAxisAngle
                					this.gIS[i].goniometerTilt,
                					this.gIS[i].goniometerAxial,
                					this.gIS[i].interAxisAngle
                			};
        					return result;
        				}
// Maybe later use both motors, for now - just the axial. It seems to have <0.5 degree error (but accumulates gradually as there are friction rollers involved).
        				int thisMotorHorizontal=this.gIS[i].motors[mHorizontal];
        				int thisMotorAxial=     this.gIS[i].motors[mAxial];
        				int stationNumber=      this.gIS[i].getStationNumber();
            			ArrayList<Integer> setList=new ArrayList<Integer>(100);
            			for (int j=0;j<this.gIS.length;j++) {
            				if (this.gIS[j]==null){
            					System.out.println("BUG?: getImagesetTiltAxial("+timeStamp+"): this.gIS["+j+"]==null");
            					continue;
            				}
            				if (this.gIS[j].motors==null){
            					System.out.println("BUG?: getImagesetTiltAxial("+timeStamp+"): this.gIS["+j+"].motors==null");
            					continue;
            				}
            				if ( //   (j!=i)  && // not needed - this set does not have orientation
            						(this.gIS[j].getStationNumber()==stationNumber) &&
            						(this.gIS[j].motors[mHorizontal]==thisMotorHorizontal) &&
            						!this.gIS[j].orientationEstimated && // new
            						!Double.isNaN(this.gIS[j].goniometerTilt) &&
            						!Double.isNaN(this.gIS[j].goniometerAxial) &&
            						!Double.isNaN(this.gIS[j].interAxisAngle)){
            					System.out.println("Found image set "+j+" with timestamp "+ timeStamp + " with known attitude to estimate this attitude");
            					setList.add(new Integer(j));
            				}
            			}
            			if (setList.size()>=2){
//            				if (this.debugLevel>2) System.out.println("getImagesetTiltAxial("+timeStamp+"): estimating orientation for set # "+i+": this.debugLevel="+this.debugLevel);
            				if (this.debugLevel > -1) System.out.println("getImagesetTiltAxial("+timeStamp+"): estimating orientation for set # "+i+": this.debugLevel="+this.debugLevel);
            				// find the closest one
            				int indexClosest=setList.get(0);
            				double dClosest=Math.abs(this.gIS[indexClosest].motors[mAxial]-thisMotorAxial);
            				for (int j=1;j<setList.size();j++) if (Math.abs(this.gIS[setList.get(j)].motors[mAxial]-thisMotorAxial)<dClosest){
            					indexClosest=setList.get(j);
            					dClosest=Math.abs(this.gIS[indexClosest].motors[mAxial]-thisMotorAxial);
            				}
            				// try to get the second on the other side than the closest first
            				int indexSecond=-1;
            				for (int j=0;j<setList.size();j++) {
            					if (((this.gIS[indexClosest].motors[mAxial]-thisMotorAxial)*
            							(this.gIS[setList.get(j)].motors[mAxial]-thisMotorAxial)<0) && // different side
            							((indexSecond<0) || (Math.abs(this.gIS[setList.get(j)].motors[mAxial]-thisMotorAxial)<dClosest))){
            						indexSecond=setList.get(j);
            						dClosest=Math.abs(this.gIS[indexSecond].motors[mAxial]-thisMotorAxial);
            					}
            				}
            				if (this.debugLevel>2) System.out.println("indexSecond="+indexSecond);
            				if (indexSecond<0){ // no sets on the opposite side from the indexClosest, use second closest on the same side as indexClosest
                				for (int j=0;j<setList.size();j++) {
                					if ((setList.get(j)!=indexClosest) &&
                							((indexSecond<0) || (Math.abs(this.gIS[setList.get(j)].motors[mAxial]-thisMotorAxial)<dClosest))){
                						indexSecond=setList.get(j);
                						dClosest=Math.abs(this.gIS[indexSecond].motors[mAxial]-thisMotorAxial);
                					}
                				}

            				}
            				if (indexSecond<0){ // no second sets at all
            					System.out.println("getImagesetTiltAxial("+timeStamp+") - this is a BUG ");
            				} else {
            					// now linear interpolate axail between theses two sets: indexClosest and indexSecond. (resolve/ guess crossing 360
            					double axialClosest=this.gIS[indexClosest].goniometerAxial;
            					double axialSecond= this.gIS[indexSecond].goniometerAxial;
            					double interClosest=this.gIS[indexClosest].interAxisAngle;
            					double interSecond= this.gIS[indexSecond].interAxisAngle;
            					axialClosest-=360.0*Math.floor((axialClosest+180.0)/360.0);
            					axialSecond-= 360.0*Math.floor((axialSecond+ 180.0)/360.0);
                				if (this.debugLevel>2) System.out.println("getImagesetTiltAxial("+timeStamp+"):"+
                						" same tilt - "+setList.size()+
                						" axialClosest="+axialClosest+
                						" axialSecond="+axialSecond+
                						" interClosest="+interClosest+
                						" interSecond="+interSecond+
                						" motor closest="+this.gIS[indexClosest].motors[mAxial]+
                						" motor second="+this.gIS[indexSecond].motors[mAxial]);
            					// axial motor has the same sign/direction as the axial angle
            					if (this.gIS[indexSecond].motors[mAxial]>this.gIS[indexClosest].motors[mAxial]){
            						if (axialSecond<axialClosest) axialSecond+=360.0;
            					} else {
            						if (axialSecond>axialClosest) axialClosest+=360.0;
            					}
            					this.gIS[i].goniometerAxial=
            						axialClosest+
            						(axialSecond-axialClosest)*
            						(thisMotorAxial-this.gIS[indexClosest].motors[mAxial])/
            						(this.gIS[indexSecond].motors[mAxial]-this.gIS[indexClosest].motors[mAxial]);
            					this.gIS[i].interAxisAngle=
            							interClosest+
                						(interSecond-interClosest)*
                						(thisMotorAxial-this.gIS[indexClosest].motors[mAxial])/
                						(this.gIS[indexSecond].motors[mAxial]-this.gIS[indexClosest].motors[mAxial]);
            					this.gIS[i].goniometerTilt=
            						this.gIS[indexClosest].goniometerTilt+
            						(this.gIS[indexSecond].goniometerTilt-this.gIS[indexClosest].goniometerTilt)*
            						(thisMotorAxial-this.gIS[indexClosest].motors[mAxial])/
            						(this.gIS[indexSecond].motors[mAxial]-this.gIS[indexClosest].motors[mAxial]);
            					// 06/06/2015 Andrey: Was missing setting estimated orientation. Was it a bug?
                    			this.gIS[i].orientationEstimated=true;
                    			if (this.debugLevel>0) System.out.println("Orientation for set # "+i+" timestamp "+IJ.d2s(this.gIS[i].timeStamp,6)+
                    					") is not defined, using interpolated between sets # "+indexClosest+" (timestamp "+IJ.d2s(this.gIS[indexClosest].timeStamp,6)+") "+
                    					"and # "+indexSecond+" (timestamp "+IJ.d2s(this.gIS[indexSecond].timeStamp,6)+")");
            				}
//            			} else if (setList.size() >= 1){
//            				if (this.debugLevel > -1) System.out.println("getImagesetTiltAxial("+timeStamp+
//            						"): estimating orientation for set # "+i+" from a single set "+setList.get(0)+": this.debugLevel="+this.debugLevel);

            			} else { // old way
            				// first try for the same station number only:
            				iBest = -1;
            				double d2Min=-1;
            				int station_number = this.gIS[i].getStationNumber();
            				for (int j=0;j<this.gIS.length;j++) if ((j!=i) &&
            						(this.gIS[j].getStationNumber() == station_number) &&
            						!this.gIS[j].orientationEstimated &&
            						(this.gIS[j].motors!=null) &&
            						!Double.isNaN(this.gIS[j].goniometerTilt) &&
            						!Double.isNaN(this.gIS[j].goniometerAxial )  &&
            						!Double.isNaN(this.gIS[j].interAxisAngle)) {
            					double d2=0;
            					for (int k=0;k<this.gIS[j].motors.length;k++){
            						d2+=1.0*(this.gIS[j].motors[k]-this.gIS[i].motors[k])*
            						(this.gIS[j].motors[k]-this.gIS[i].motors[k]);
            					}
            					if ((d2Min<0) || (d2Min>d2)) {
            						d2Min=d2;
            						iBest=j;
            					}
            				}
            				if (iBest < 0) {
                				d2Min=-1;
                				for (int j=0;j<this.gIS.length;j++) if ((j!=i) &&
                						(this.gIS[j].motors!=null) &&
                						!this.gIS[j].orientationEstimated &&
                						!Double.isNaN(this.gIS[j].goniometerTilt) &&
                						!Double.isNaN(this.gIS[j].goniometerAxial )  &&
                						!Double.isNaN(this.gIS[j].interAxisAngle)) {
                					double d2=0;
                					for (int k=0;k<this.gIS[j].motors.length;k++){
                						d2+=1.0*(this.gIS[j].motors[k]-this.gIS[i].motors[k])*
                						(this.gIS[j].motors[k]-this.gIS[i].motors[k]);
                					}
                					if ((d2Min<0) || (d2Min>d2)) {
                						d2Min=d2;
                						iBest=j;
                					}
                				}
                				if (iBest < 0) {
                    				d2Min=-1;
                    				for (int j=0;j<this.gIS.length;j++) if ((j!=i) &&
                    						(this.gIS[j].motors!=null) &&
                    						!Double.isNaN(this.gIS[j].goniometerTilt) &&
                    						!Double.isNaN(this.gIS[j].goniometerAxial )  &&
                    						!Double.isNaN(this.gIS[j].interAxisAngle)) {
                    					double d2=0;
                    					for (int k=0;k<this.gIS[j].motors.length;k++){
                    						d2+=1.0*(this.gIS[j].motors[k]-this.gIS[i].motors[k])*
                    						(this.gIS[j].motors[k]-this.gIS[i].motors[k]);
                    					}
                    					if ((d2Min<0) || (d2Min>d2)) {
                    						d2Min=d2;
                    						iBest=j;
                    					}
                    				}
                    				System.out.println("Used any station numer with even estimated orientation, iBest = "+iBest);
                				} else {
                    				System.out.println("Used different station numer, iBest = "+iBest);
                				}
            				} else {
                				System.out.println("Used the same station number, iBest = "+iBest);

            				}


            			}
        			}
//        			double [] result = {
//        					this.gIS[iBest].goniometerTilt,
//        					this.gIS[iBest].goniometerAxial,
//        					this.gIS[iBest].interAxisAngle
//        			};
        			if (iBest!=i){
        				boolean usable_tilt  = (this.gIS[i].motors != null);
        				boolean usable_axial = usable_tilt && (this.gIS[i].getStationNumber() == this.gIS[iBest].getStationNumber());
        				double diff_axial = usable_axial? (this.gIS[i].motors[mAxial]-this.gIS[iBest].motors[mAxial])/
        						goniometerParameters.goniometerMotors.stepsPerDegreeAxial : 0.0;
        				double diff_horizontal = usable_tilt? (this.gIS[i].motors[mHorizontal]-this.gIS[iBest].motors[mHorizontal])/
        						goniometerParameters.goniometerMotors.stepsPerDegreeTilt: 0.0;
        				this.gIS[i].goniometerTilt =  this.gIS[iBest].goniometerTilt + diff_horizontal;
        				this.gIS[i].goniometerAxial = this.gIS[iBest].goniometerAxial + diff_axial;

        				this.gIS[i].goniometerTilt-=360.0*Math.floor((this.gIS[i].goniometerTilt+180.0)/360.0);
        				this.gIS[i].goniometerAxial-=360.0*Math.floor((this.gIS[i].goniometerAxial+180.0)/360.0);

            			if (this.debugLevel>0) System.out.println("Orientation for set # "+i+" timestamp "+IJ.d2s(this.gIS[i].timeStamp,6)+
            					") is not defined, estimating from  # "+iBest+" (timestamp "+IJ.d2s(this.gIS[iBest].timeStamp,6)+")" );
            			this.gIS[i].orientationEstimated=true;
//    					this.gIS[i].goniometerTilt= this.gIS[iBest].goniometerTilt;
//    					this.gIS[i].goniometerAxial=this.gIS[iBest].goniometerAxial;
    					this.gIS[i].interAxisAngle=this.gIS[iBest].interAxisAngle;
        			}
        			double [] result = {
        					this.gIS[i].goniometerTilt,
        					this.gIS[i].goniometerAxial,
        					this.gIS[i].interAxisAngle
        			};

       				return result; // may have Double.NaN
        	}
        	return null;
        }

        public double getImageTimestamp(ImagePlus image){
        	if ((image.getProperty("timestamp")==null) || (((String) image.getProperty("timestamp")).length()==0)) {
        		(new JP46_Reader_camera(false)).decodeProperiesFromInfo(image);
        	}
        	return Double.parseDouble((String) image.getProperty("timestamp"));
        }

        public int getImageChannel(ImagePlus image){
        	if ((image.getProperty("channel")==null) || (((String) image.getProperty("channel")).length()==0)) {
        		(new JP46_Reader_camera(false)).decodeProperiesFromInfo(image);
        	}

        	String channelSuffix=(String) image.getProperty("channel");
        	int channel=-1;
        	for (int j=0;j<this.channelSuffixes.length;j++){
//        		System.out.println("== j="+j);
//        		System.out.println("channelSuffix="+channelSuffix);
//        		System.out.println("this.channelSuffixes[j]="+this.channelSuffixes[j]);
        		if (channelSuffix.equals(this.channelSuffixes[j])) {
        			channel=j;
        			break;
        		}
        	}
        	if (channel<0) {
        		String msg="Channel not recognized) - this channel suffix is "+channelSuffix+", available channel suffixes are:\n";
        		for (int j=0;j<this.channelSuffixes.length;j++) msg+=this.channelSuffixes[j]+", ";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return channel;
        }

        /**
         * initialize image data with camera defaults
         * @param distortionCalibrationData grid distortionCalibrationData
         * @param eyesisCameraParameters deafault camera parameters
         * @return
         */
        // Used in Goniometer
        public void initImageSet(
        		EyesisCameraParameters eyesisCameraParameters) {
        	for (int i=0;i<this.getNumImages();i++){
        		int subCam=this.getImageSubcamera(i);
        		int stationNumber=this.getImageStation(i);
        		this.setParameters(eyesisCameraParameters.getParametersVector(stationNumber,subCam), i);
        	}
        }



// constructor from XML file

        public DistortionCalibrationData (
        		boolean smart,
        		String defaultPath,
        		PatternParameters patternParameters,
        		EyesisCameraParameters eyesisCameraParameters,
    			EyesisAberrations.AberrationParameters aberrationParameters,
    			LaserPointer laserPointers,
    			Goniometer.GoniometerParameters goniometerParameters,
				ImagePlus[] gridImages  ){ // null - use specified files
        	this.goniometerParameters = goniometerParameters;
        	setupIndices();
			String [] extensions={".dcal-xml","-distcal.xml"};
			MultipleExtensionsFileFilter parFilter = new MultipleExtensionsFileFilter("",extensions,"Distortion calibration *.dcal-xml files");
			String pathname=CalibrationFileManagement.selectFile(
					smart,
					false,
					"Restore Calibration Parameters",
					"Restore",
					parFilter,
					defaultPath); //String defaultPath
			if ((pathname==null) || (pathname=="")) return;
//			setGridImages(gridImages);
//TODO: these images will be overwritten by setFromXML !!!!!!!!!
			this.gIS=null; // So readAllGrids will create it
        	setFromXML(
        			pathname,
            		eyesisCameraParameters,
        			aberrationParameters);
			if (gridImages!=null) {
//				this.pathName="";  // modified, keep the path anyway
// overwrite saved paths with the provided images, number of images{ should match
				if (this.gIP.length == gridImages.length){
					for (int i=0;i<this.gIP.length;i++){
						this.gIP[i].gridImage=gridImages[i];
						this.gIP[i].path=null; // not needed, just in case
						this.gIP[i].enabled=true;// enable all (actually just one) acquired images
					}
				} else {
					String msg="Number of provided images ("+gridImages.length+") does not match parameters restored from the "+pathname+" ("+this.gIP.length+")";
		    		IJ.showMessage("Error",msg);
//		    		throw new IllegalArgumentException (msg);
					for (int i=0; i<this.gIP.length ; i++){
						this.gIP[i].path=null; // not needed, just in case
						this.gIP[i].enabled=true;// enable all (actually just one) acquired images
						if (i < gridImages.length) {
							this.gIP[i].gridImage=gridImages[i];
						} else {
							this.gIP[i].gridImage=null;
						}
					}
				}
//				setGridImages(gridImages);
			}
        	readAllGrids(
        			patternParameters, // prepare grid parameters for LMA
        			laserPointers, // prepare grid parameters for LMA
            		true); // boolean keep_images make it configurable parameter?
			updateSetOrientation(null); // update orientation of image sets (built in readAllGrids() UPDATE - not anymore)

        }

        public void setFromXML(String pathname,
        		EyesisCameraParameters eyesisCameraParameters, // should have cartesian set
    			EyesisAberrations.AberrationParameters aberrationParameters) {
        	this.eyesisCameraParameters=eyesisCameraParameters;

        	XMLConfiguration hConfig=null;
        	try {
				hConfig=new XMLConfiguration(pathname);
			} catch (ConfigurationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    		this.numSubCameras=Integer.parseInt(hConfig.getString("subcameras","1"));
        	System.out.println("Number of subcameras is "+this.numSubCameras);
			int num=hConfig.getMaxIndex("file");
			num++;
        	this.gIP=new GridImageParameters[num];
        	this.pars=new double[num][parameterDescriptions.length];
        	System.out.println("Number of pattern grid images in "+pathname+" is "+num);

        	int numSets=hConfig.getMaxIndex("set")+1; // see if it returns -1 for none
        	System.out.println("Number of image sets in "+pathname+" is "+numSets);
        	if (numSets>0){
            	this.gIS=new GridImageSet[numSets];
            	for (int i=0;i<numSets;i++) {
            		HierarchicalConfiguration sub = hConfig.configurationAt("set("+i+")");
            		int index=Integer.parseInt(sub.getString("index"));
            		this.gIS[index]=new GridImageSet();
            		this.gIS[index].timeStamp=     Double.parseDouble(sub.getString("timestamp"));
            		this.gIS[index].stationNumber= Integer.parseInt(sub.getString("stationNumber"));
                	int minIndex=       this.gIS[index].getMinIndex();
                	int maxIndexPlusOne=this.gIS[index].getMaxIndexPlusOne();
                	for (int j=minIndex;j<maxIndexPlusOne;j++){
//                		if (sub.getString(parameterDescriptions[j][0])!=null) {
                   		if (sub.getString(descrField(j,0)) != null) {
                			this.gIS[index].setParameterValue(j,Double.parseDouble(sub.getString(descrField(j,0))), false);
                		}
                	}

            		if (sub.getString("orientationEstimated")!=null) {
            			this.gIS[i].orientationEstimated=Boolean.parseBoolean(sub.getString("orientationEstimated"));
/*                    	System.out.println(i+": restored orientationEstimated="+this.gIS[i].orientationEstimated+
                    			" tilt="+this.gIS[i].goniometerTilt+
                    			" axial="+this.gIS[i].goniometerAxial);*/
            		} else {
                		this.gIS[i].setEstimatedFromNonNaN();
/*                    	System.out.println(i+": guessed. orientationEstimated="+this.gIS[i].orientationEstimated+
                    			" tilt="+this.gIS[i].goniometerTilt+
                    			" axial="+this.gIS[i].goniometerAxial);*/
            		}

            	}

        	} else {
        		this.gIS=null; // has to be build later
        	}

        	for (int i=0;i<num;i++) {
        		this.gIP[i]=new GridImageParameters(i);
        		HierarchicalConfiguration sub = hConfig.configurationAt("file("+i+")");
        		this.gIP[i].imgNumber=i;
        		this.gIP[i].path=sub.getString("name");
        		this.gIP[i].source_path=sub.getString("source_path");
        		this.gIP[i].timestamp=Double.parseDouble(sub.getString("timestamp"));
        		this.gIP[i].channel=Integer.parseInt(sub.getString("channel"));
        		if (sub.getString("stationNumber")!=null) this.gIP[i].setStationNumber(Integer.parseInt(sub.getString("stationNumber")));
        		else this.gIP[i].setStationNumber(0);
        		if (sub.getString("enabled")!=null) this.gIP[i].enabled=Boolean.parseBoolean(sub.getString("enabled"));
        		if (sub.getString("noUsefulPSFKernels")!=null) this.gIP[i].noUsefulPSFKernels=Boolean.parseBoolean(sub.getString("noUsefulPSFKernels"));
        		this.gIP[i].setNumber=sub.getInt("setNumber",-1);
// new
        		this.gIP[i].hintedMatch=sub.getInt("hintedMatch",-1);
        		this.gIP[i].enabled=sub.getBoolean("enabled",false);
//        		if (aberrationParameters.trustEnabled && this.gIP[i].enabled) this.gIP[i].hintedMatch=2; // trusted
        		if (aberrationParameters.trustEnabled) this.gIP[i].hintedMatch= this.gIP[i].enabled?2:-1; // trusted and only trusted to enabled

        		for (int j=0;j<this.parameterDescriptions.length;j++){
//            		if (sub.getString(parameterDescriptions[j][0])!=null)
               		if (sub.getString(descrField(j,0))!=null)
        				this.pars[i][j] = Double.parseDouble(sub.getString(descrField(j,0)));
        			else
        				if (isNonRadial(j)){
        					this.pars[i][j] = 0.0; // old calibration files without non-radial parameters
        				} else {
        					this.pars[i][j] = Double.NaN;
        				}
        		}
        		int [] shiftRot={
        				sub.getInt("gridShiftX", 0),
        				sub.getInt("gridShiftY", 0),
        				sub.getInt("gridRotate", 0)};
        		this.gIP[i].setUVShiftRot(shiftRot);
//        		getInt(String key, int defaultValue)
        	}
        	if (this.gIS!=null){
            	System.out.println("Using stored image set data");
        		for (int is=0;is<this.gIS.length;is++){
            		this.gIS[is].imageSet=new GridImageParameters [this.numSubCameras];
            		for (int j=0;j<this.numSubCameras;j++) this.gIS[is].imageSet[j]=null;
        		}
        		for (int ip=0;ip<this.gIP.length;ip++) if (this.gIP[ip].setNumber>=0) {
        			this.gIS[this.gIP[ip].setNumber].imageSet[this.gIP[ip].channel]=this.gIP[ip];
        			this.gIP[ip].gridImageSet=this.gIS[this.gIP[ip].setNumber];
        			//this.gIP[i].channel
        		}

        	} else {
            	System.out.println("Re-creating image set data from individual images (old format)");
            	System.out.println("WARNING: Some parameters may get from unused images and so have wrong values");
            	buildImageSets(false); // from scratch
            	// copying only parameters that have the same values for all images in a set
            	for (int is=0;is<this.gIS.length;is++){
            		int minIndex=       this.gIS[is].getMinIndex();
            		int maxIndexPlusOne=this.gIS[is].getMaxIndexPlusOne();
            		for (int pi=minIndex;pi<maxIndexPlusOne;pi++){
                		double parVal=Double.NaN;
                		boolean differs=false;
            			for (int j=0;j<this.gIS[is].imageSet.length;j++) if (this.gIS[is].imageSet[j]!=null) {
            				int imgNum=this.gIS[is].imageSet[j].imgNumber;
            				if (!Double.isNaN(this.pars[imgNum][pi])){
            					if (!Double.isNaN(parVal) && (parVal!=this.pars[imgNum][pi])){
            						differs=true;
            						break;
            					} else {
            						parVal=this.pars[imgNum][pi];
            					}
            				}
            				if (!differs && !Double.isNaN(parVal)){
            					this.gIS[is].setParameterValue(pi,parVal,false);
            				}
            				if (differs){
//            					System.out.println("ImageSet #"+is+": "+parameterDescriptions[j][0] +" has different values for individual images, skipping");
            					System.out.println("ImageSet #"+is+": "+descrField(j,0) +" has different values for individual images, skipping");
            				}
            			}
            		}
            		this.gIS[is].setEstimatedFromNonNaN();
            		//orientationEstimated
            		System.out.println(is+": tilt="+this.gIS[is].goniometerTilt+" axial="+this.gIS[is].goniometerAxial+" estimated="+this.gIS[is].orientationEstimated);

            	}
            	System.out.println("setFromXML("+pathname+",eyesisCameraParameters) 1 -> this.gIS.length="+this.gIS.length);
        	}
//        	System.out.println("setFromXML("+pathname+",eyesisCameraParameters) 2 -> this.gIS.length="+((this.gIS==null)?"null":this.gIS.length));
        	this.pathName=pathname; // where this instance was created from
        }
  //http://commons.apache.org/configuration/userguide/howto_xml.html
        public String getPath(){
        	return this.pathName;
        }
        public String selectAndSaveToXML(boolean smart, String defaultPath){
        	return selectAndSaveToXML(smart, defaultPath, null);
        }
        public String selectAndSaveToXML(boolean smart, String defaultPath, String comment){
			String [] extensions={".dcal-xml","-distcal.xml"};
			MultipleExtensionsFileFilter parFilter = new MultipleExtensionsFileFilter("",extensions,"Distortion calibration *.dcal-xml files");
			String pathname=CalibrationFileManagement.selectFile(
					smart,
					true,
					"Save Calibration Parameters",
					"Save",
					parFilter,
					(defaultPath==null)?this.pathName:defaultPath); //String defaultPath
			if (pathname!=null) saveToXML(pathname,comment);
			return pathname;
        }
        public boolean saveTimestampedToXML(String pathname, String comment) {
        	return saveToXML(pathname+"_"+IJ.d2s(0.000001*(System.nanoTime()/1000),6).replace('.', '_')+".dcal-xml",      // full path or null
        			null);
        }
        public boolean saveToXML(String pathname) {
        	return saveToXML(pathname,null);
        }
        public boolean saveToXML(String pathname, String comment) {
        	XMLConfiguration hConfig=new XMLConfiguration();
        	if (comment!=null) hConfig.addProperty("comment",comment);
        	hConfig.setRootElementName("distortionCalibrationParameters");
        	hConfig.addProperty("subcameras",this.numSubCameras);
        	for (int i=0;i<this.gIP.length;i++){
            	hConfig.addProperty("file","");
            	hConfig.addProperty("file.setNumber",this.gIP[i].setNumber);
            	hConfig.addProperty("file.name",this.gIP[i].path);
            	hConfig.addProperty("file.source_path",this.gIP[i].source_path);
            	hConfig.addProperty("file.enabled",this.gIP[i].enabled);
            	hConfig.addProperty("file.hintedMatch",this.gIP[i].hintedMatch); // new
            	hConfig.addProperty("file.timestamp",IJ.d2s(this.gIP[i].timestamp,6));
            	hConfig.addProperty("file.channel",this.gIP[i].channel);
            	hConfig.addProperty("file.stationNumber",this.gIP[i].getStationNumber());
            	hConfig.addProperty("file.noUsefulPSFKernels",this.gIP[i].noUsefulPSFKernels);
            	int [] UVShiftRot=this.gIP[i].getUVShiftRot();
            	hConfig.addProperty("file.gridShiftX",UVShiftRot[0]);
            	hConfig.addProperty("file.gridShiftY",UVShiftRot[1]);
            	hConfig.addProperty("file.gridRotate",UVShiftRot[2]);
            	for (int j=0;j<this.parameterDescriptions.length;j++){
                	hConfig.addProperty("file."+descrField(j,0),this.pars[i][j]);
            	}
        	}
// save image sets
        	for (int i=0;i<this.gIS.length;i++){
            	hConfig.addProperty("set","");
            	hConfig.addProperty("set.index",i);
            	hConfig.addProperty("set.stationNumber",this.gIS[i].stationNumber);
            	hConfig.addProperty("set.timestamp",    IJ.d2s(this.gIS[i].timeStamp,6));
            	hConfig.addProperty("set.orientationEstimated",this.gIS[i].orientationEstimated);
            	double [] vector = this.gIS[i].updateParameterVectorFromSet(null); // unused parameters will be NaN
            	for (int j=0;j<vector.length;j++) if (!Double.isNaN(vector[j])){
            		hConfig.addProperty("set."+descrField(j,0),vector[j]);
            	}
        	}

//        	hConfig.addProperty("grids","");
        	File file=new File (pathname);
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
			this.pathName=pathname;
        	return true;
        }

//    		public double      gridPeriod=0.0;  // average grid period, in pixels (to filter out (double-) reflected images
        public double calcGridPeriod(
        		int fileNumber,
        		boolean use_extra){ // use out-of grid nodes (can be w/o absolute matching
        	use_extra &= (this.gIP[fileNumber].pixelsXY_extra !=null);
        	int len = (this.gIP[fileNumber].pixelsXY==null)? 0 : this.gIP[fileNumber].pixelsXY.length;
        	int len0 = len;
        	if (use_extra ) {
        		len += this.gIP[fileNumber].pixelsXY_extra.length;
        	}
        	if (len<3) {
        		this.gIP[fileNumber].setGridPeriod(Double.NaN);
        	} else {
//        		double [][][] data =new double [this.gIP[fileNumber].pixelsXY.length][2][2];
        		double [][][] data =new double [len][2][2];
        		// U(x,y), v(x,y)
        		for (int i=0; i < this.gIP[fileNumber].pixelsXY.length; i++){
        			data[i][0][0]=this.gIP[fileNumber].pixelsXY[i][0];
        			data[i][0][1]=this.gIP[fileNumber].pixelsXY[i][1];
        			data[i][1][0]=this.gIP[fileNumber].pixelsUV[i][0];
        			data[i][1][1]=this.gIP[fileNumber].pixelsUV[i][1];
        		}
        		if (use_extra) {
            		for (int i=0; i < this.gIP[fileNumber].pixelsXY_extra.length; i++){
            			data[i + len0][0][0]=this.gIP[fileNumber].pixelsXY_extra[i][0];
            			data[i + len0][0][1]=this.gIP[fileNumber].pixelsXY_extra[i][1];
            			data[i + len0][1][0]=this.gIP[fileNumber].pixelsUV_extra[i][0];
            			data[i + len0][1][1]=this.gIP[fileNumber].pixelsUV_extra[i][1];
            		}
        		}
        		if (this.debugLevel>3) {
        			System.out.println("calcGridPeriod("+fileNumber+"), debugLevel="+this.debugLevel+":");
            		for (int i=0;i<data.length;i++)System.out.println(i+": {{"+data[i][0][0]+","+data[i][0][1]+"},{"+data[i][1][0]+","+data[i][1][1]+"}}");
        		}
         	   double [][] coeff=new PolynomialApproximation(this.debugLevel).quadraticApproximation(data, true); // force linear
         	   if (coeff!=null) {
         	     this.gIP[fileNumber].setGridPeriod(2.0/Math.sqrt(coeff[0][0]*coeff[0][0]+coeff[0][1]*coeff[0][1]+coeff[1][0]*coeff[1][0]+coeff[1][1]*coeff[1][1]));
         	     if (this.debugLevel>3) {
         	    	System.out.println("coeff[][]={{"+coeff[0][0]+","+coeff[0][1]+"},{"+coeff[1][0]+","+coeff[1][1]+"}}");
         	     }
         	   } else {
        		  this.gIP[fileNumber].setGridPeriod(Double.NaN);
         	   }
        	}
    		if (this.debugLevel>3) {
    			System.out.println("calcGridPeriod("+fileNumber+") => "+this.gIP[fileNumber].getGridPeriod());
    		}
        	return this.gIP[fileNumber].getGridPeriod();

        }

        public double [] getGridUVfromXY(
        		final double px,
        		final double py,
        		final int fileNumber,
        		boolean use_extra){ // use out-of grid nodes (can be w/o absolute matching
        	use_extra &= (this.gIP[fileNumber].pixelsXY_extra !=null);
        	int len = (this.gIP[fileNumber].pixelsXY==null)? 0 : this.gIP[fileNumber].pixelsXY.length;
        	final int len0 = len;
        	if (use_extra ) {
        		len += this.gIP[fileNumber].pixelsXY_extra.length;
        	}
        	if (len<3) {
        		return null;
        	}
        	final double [][] all_xy = new double [len][2];
        	final int    [][] all_uv = new int   [len][2];
        	if (len0 > 0) {
        		System.arraycopy(this.gIP[fileNumber].pixelsXY, 0, all_xy, 0, len0);
        		System.arraycopy(this.gIP[fileNumber].pixelsUV, 0, all_uv, 0, len0);
        	}
        	if (len > len0) {
        		System.arraycopy(this.gIP[fileNumber].pixelsXY_extra, 0, all_xy, len0, len-len0);
        		System.arraycopy(this.gIP[fileNumber].pixelsUV_extra, 0, all_uv, len0, len-len0);
        	}

    		ArrayList<Integer> neibs = new ArrayList<Integer>();
    		for (int i = 0; i < len; i++) neibs.add(i);
    		Collections.sort(neibs, new Comparator<Integer>() {
    		    @Override
    		    public int compare(Integer lhs, Integer rhs) {
//    		    	double [] xy_lhs = (lhs < len0)? gIP[fileNumber].pixelsXY[lhs] : gIP[fileNumber].pixelsXY_extra[lhs-len0];
//    		    	double [] xy_rhs = (rhs < len0)? gIP[fileNumber].pixelsXY[rhs] : gIP[fileNumber].pixelsXY_extra[rhs-len0];
    		    	double x_lhs = all_xy[lhs][0] - px;
    		    	double y_lhs = all_xy[lhs][1] - py;
    		    	double x_rhs = all_xy[rhs][0] - px;
    		    	double y_rhs = all_xy[rhs][1] - py;
    		    	double l2_lhs = x_lhs*x_lhs + y_lhs*y_lhs;
    		    	double l2_rhs = x_rhs*x_rhs + y_rhs*y_rhs;
    		        return l2_rhs > l2_lhs ? -1 : (l2_rhs < l2_lhs) ? 1 : 0;
    		    }
    		});
    		// now list neibs start with closest to px,py node. Get first with non-collinearU,V
    		int dlen = -1;
    		int [] i01 = {neibs.get(0), neibs.get(1)};
    		int [] duv0 = {all_uv[i01[1]][0]-all_uv[i01[0]][0], all_uv[i01[1]][1]-all_uv[i01[0]][1]};
    		for (int ii = 2;ii < len; ii++) {
    			int i = neibs.get(ii);
//        		int [] duv = {all_uv[i][0]-all_uv[0][0], all_uv[i][1]-all_uv[0][1]};
        		int idet = (duv0[0] *  (all_uv[i][1]-all_uv[i01[0]][1]))
        				-  (duv0[1] *  (all_uv[i][0]-all_uv[i01[0]][0]));
        		if (idet != 0) {
        			dlen = i + 1;
        			break;
        		}
    		}
    		if (dlen < 0) {
    			return null;
    		}
    		double [][][] data =new double [dlen][2][2];
    		// U(x,y), v(x,y)
    		for (int i=0; i < dlen; i++){
    			int indx = neibs.get(i);
    			data[i][0][0]=all_xy[indx][0];
    			data[i][0][1]=all_xy[indx][1];
    			data[i][1][0]=all_uv[indx][0];
    			data[i][1][1]=all_uv[indx][1];
    		}
      	   double [][] coeff=new PolynomialApproximation(this.debugLevel).quadraticApproximation(data, true); // force linear
     	   if (coeff == null) {
     		   return null;
     	   }
     	   double [] uv0 =
     		   {       (coeff[0][0] * px + coeff[0][1] * py + coeff[0][2]),
     				   (coeff[1][0] * px + coeff[1][1] * py + coeff[1][2])};

     	   int [][] reMap= MatchSimulatedPattern.getRemapMatrix(this.gIP[fileNumber].getUVShiftRot());
//     	   double [] uv = { reMap[0][0]*uv0[0] + reMap[0][1]* uv0[1] + reMap[0][2], // u
//     			   (        reMap[1][0]*uv0[0] + reMap[1][1]* uv0[1] + reMap[1][2])}; // v;
//     	   Sign?
     	   double [] uv = { reMap[0][0]*uv0[0] + reMap[0][1]* uv0[1] - reMap[0][2], // u
     			   (        reMap[1][0]*uv0[0] + reMap[1][1]* uv0[1] - reMap[1][2])}; // v;

        	return uv;
        }



        public int [] setGridsWithRemap(
        		int fileNumber,
        		int [][] reMap,
        		float [][] pixels,
        		PatternParameters patternParameters){
    		int sensorWidth=this.eyesisCameraParameters.getSensorWidth(this.gIP[fileNumber].channel);
    		int sensorHeight=this.eyesisCameraParameters.getSensorHeight(this.gIP[fileNumber].channel);
        	int station=this.gIP[fileNumber].getStationNumber();
        	int size=0;
        	int size_extra=0;
        	for (int i=0;i<pixels[0].length;i++) if ((pixels[0][i]>=0) && (pixels[1][i]>=0) && (pixels[0][i]<sensorWidth) && (pixels[1][i]<sensorHeight)){
        		int u=Math.round(pixels[2][i]);
        		int v=Math.round(pixels[3][i]);
    			int u1= reMap[0][0]*u + reMap[0][1]*v + reMap[0][2]; // u
    			int v1= reMap[1][0]*u + reMap[1][1]*v + reMap[1][2]; // v;
        		if (patternParameters.getXYZM(u1,v1,false,station)!=null) size++; // already assumes correct uv?
        		else size_extra++;
        	}


        	this.gIP[fileNumber].resetMask();
        	this.gIP[fileNumber].pixelsXY=new double [size][6];
        	this.gIP[fileNumber].pixelsUV=new int    [size][2];
        	this.gIP[fileNumber].pixelsXY_extra=new double [size_extra][6];
        	this.gIP[fileNumber].pixelsUV_extra=new int    [size_extra][2];

        	int index=0;
        	int index_extra=0;
        	for (int i=0;i<pixels[0].length;i++) if ((pixels[0][i]>=0) && (pixels[1][i]>=0) && (pixels[0][i]<sensorWidth) && (pixels[1][i]<sensorHeight)) {
        		int u=Math.round(pixels[2][i]);
        		int v=Math.round(pixels[3][i]);
    			int u1= reMap[0][0]*u + reMap[0][1]*v + reMap[0][2]; // u
    			int v1= reMap[1][0]*u + reMap[1][1]*v + reMap[1][2]; // v;

        		if (patternParameters.getXYZM(u1,v1,false,station)!=null) {
        			this.gIP[fileNumber].pixelsXY[index][0]=pixels[0][i];
        			this.gIP[fileNumber].pixelsXY[index][1]=pixels[1][i];
        			this.gIP[fileNumber].pixelsUV[index][0]= u1; // u
        			this.gIP[fileNumber].pixelsUV[index][1]= v1; // v;
        			if (this.gIP[fileNumber].flatFieldAvailable){
        				this.gIP[fileNumber].pixelsXY[index][2]=pixels[4][i];
        				for (int n=0;n<3;n++) this.gIP[fileNumber].pixelsXY[index][n+3]=pixels[n+5][i]/this.gIP[fileNumber].intensityRange[n];
        			} else {
        				for (int n=0;n<4;n++)this.gIP[fileNumber].pixelsXY[index][n+2]=1.0;
        			}
        			index++;
        		} else {
        			this.gIP[fileNumber].pixelsXY_extra[index_extra][0]=pixels[0][i];
        			this.gIP[fileNumber].pixelsXY_extra[index_extra][1]=pixels[1][i];
        			this.gIP[fileNumber].pixelsUV_extra[index_extra][0]= u1; // u
        			this.gIP[fileNumber].pixelsUV_extra[index_extra][1]= v1; // v;
        			if (this.gIP[fileNumber].flatFieldAvailable){
        				this.gIP[fileNumber].pixelsXY_extra[index_extra][2]=pixels[4][i];
        				for (int n=0;n<3;n++){
        					this.gIP[fileNumber].pixelsXY_extra[index_extra][n+3]=pixels[n+5][i]/this.gIP[fileNumber].intensityRange[n];
        				}
        			} else {
        				for (int n=0;n<4;n++)this.gIP[fileNumber].pixelsXY_extra[index_extra][n+2]=1.0;
        			}
        			index_extra++;
        		}
        	}
        	int [] result = {size,size_extra};
        	return result;
        }


        public boolean readAllGrids(
        		PatternParameters patternParameters,
        		LaserPointer      laserPointers, // as a backup if data is not available in the file
        		boolean keep_images
            ){
        	boolean disableNoFlatfield=false;  // true only for processing transitional images - mixture of ff/ no-ff
			System.out.println("readAllGrids(), this.debugLevel="+this.debugLevel+" this.gIS is "+((this.gIS==null)?"null":"not null"));
        	int numImages=getNumImages();
    		Opener opener=new Opener();
    		JP46_Reader_camera jp4_reader= new JP46_Reader_camera(false);
    		ImagePlus imp_grid=null;
    		ImageStack stack;
    		int numOfGridNodes=0;
    		int numOfGridNodes_extra=0;
        	for (int fileNumber=0;fileNumber<numImages;fileNumber++){
        		boolean woi_compensated = true;
        		if (this.gIP[fileNumber].gridImage!=null){ // use in-memory grid images instead of the files
        			int numGridImg=fileNumber;
        			if (numGridImg>=this.gIP.length) numGridImg=this.gIP.length-1;
        			if (this.updateStatus) IJ.showStatus("Using in-memory grid image "+(fileNumber+1)+" (of "+(numImages)+"): "+
        					this.gIP[numGridImg].gridImage.getTitle());
        			if (this.debugLevel>1) System.out.print((fileNumber+1)+": "+this.gIP[numGridImg].gridImage.getTitle());
        			imp_grid=this.gIP[numGridImg].gridImage;
        			// in memory grid pixels are always compensated
        		} else {
        			if (this.updateStatus) IJ.showStatus("Reading grid file "+(fileNumber+1)+" (of "+(numImages)+"): "+this.gIP[fileNumber].path);
        			if (this.debugLevel>-1) System.out.print(fileNumber+" ("+this.gIP[fileNumber].getStationNumber()+"): "+this.gIP[fileNumber].path);
        			imp_grid=opener.openImage("", this.gIP[fileNumber].path);  // or (path+filenames[nFile])
        			if (imp_grid==null) {
        				String msg="Failed to read grid file "+this.gIP[fileNumber].path;
        				IJ.showMessage("Error",msg);
        				throw new IllegalArgumentException (msg);
        			}
// TODO: here - need to decode properties
        			jp4_reader.decodeProperiesFromInfo(imp_grid);
//        			woi_compensated = getImagePlusProperty(imp_grid,"WOI_COMPENSATED",false);
        			woi_compensated = (getImagePlusProperty(imp_grid,"WOI_TOP",0) == 0) &&
        					(getImagePlusProperty(imp_grid,"WOI_LEFT",0) == 0);

        			if (keep_images) {
        				this.gIP[fileNumber].gridImage = imp_grid;
        			}
        		}
    			this.gIP[fileNumber].woi = new Rectangle(
    					getImagePlusProperty(imp_grid,"WOI_LEFT",0),
    					getImagePlusProperty(imp_grid,"WOI_TOP",0),
    					getImagePlusProperty(imp_grid,"WOI_WIDTH",  eyesisCameraParameters.getSensorWidth(this.gIP[fileNumber].getChannel())),
    					getImagePlusProperty(imp_grid,"WOI_HEIGHT", eyesisCameraParameters.getSensorHeight(this.gIP[fileNumber].getChannel())));
        		this.gIP[fileNumber].laserPixelCoordinates=MatchSimulatedPattern.getPointersXYUV(imp_grid, laserPointers);
        		this.gIP[fileNumber].motors=getMotorPositions(imp_grid, this.numMotors);
        		this.gIP[fileNumber].matchedPointers=getUsedPonters(imp_grid);
        		double [] saturations=new double [4];
        		for (int i=0;i<saturations.length;i++) {
        			saturations[i]=Double.NaN;
        			if (imp_grid.getProperty("saturation_" + i) !=null) saturations[i]=Double.parseDouble((String) imp_grid.getProperty("saturation_" + i));
        		}
        		if (!Double.isNaN(saturations[1])) this.gIP[fileNumber].saturation[0]=saturations[1];
        		if (!Double.isNaN(saturations[2])) this.gIP[fileNumber].saturation[2]=saturations[2];
        		if (!Double.isNaN(saturations[0]) && !Double.isNaN(saturations[3])) this.gIP[fileNumber].saturation[1]=0.5*(saturations[0]+saturations[3]);
        		else {
            		if (!Double.isNaN(saturations[0])) this.gIP[fileNumber].saturation[1]=saturations[0];
            		if (!Double.isNaN(saturations[3])) this.gIP[fileNumber].saturation[1]=saturations[3];
        		}

                stack=imp_grid.getStack();
            	if ((stack==null) || (stack.getSize()<4)) {
            		String msg="Expected a 8-slice stack in "+this.gIP[fileNumber].path;
            		IJ.showMessage("Error",msg);
            		throw new IllegalArgumentException (msg);
            	}
        		float [][] pixels=new float[stack.getSize()][]; // now - 8 (x,y,u,v,contrast, vignR,vignG,vignB
            	for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1); // pixel X : negative - no grid here
    			if (!woi_compensated) {
    				System.out.print(" woi_compensation-d "+fileNumber+" ");
    				for (int i = 0; i < pixels[0].length; i++) {
    					pixels[0][i] += this.gIP[fileNumber].woi.x;
    					pixels[1][i] += this.gIP[fileNumber].woi.y;
    				}
    				this.gIP[fileNumber].woi.width += this.gIP[fileNumber].woi.x;
    				this.gIP[fileNumber].woi.x = 0;
    				this.gIP[fileNumber].woi.height += this.gIP[fileNumber].woi.y;
    				this.gIP[fileNumber].woi.y = 0;
    				woi_compensated = true;

    			}

            	if (this.eyesisCameraParameters.badNodeThreshold>0.0){
            		boolean thisDebug =false;
//            		thisDebug|=        (fileNumber== 720); // chn 25
                 int numBadNodes=fixBadGridNodes(
                		pixels,
                		stack.getWidth(),
                		this.eyesisCameraParameters.badNodeThreshold,
                		this.eyesisCameraParameters.maxBadNeighb,
                		this.debugLevel+(thisDebug?3:0),
                		thisDebug?("fixBad-"+fileNumber):null
                		);
                 if (this.debugLevel>-1) {
                  if (numBadNodes>0)
                	  System.out.print("  -- replaced "+numBadNodes+" bad grid nodes");
                  int [] uvrot=this.gIP[fileNumber].getUVShiftRot();
                  System.out.println(" shift:rot="+uvrot[0]+"/"+uvrot[1]+":"+uvrot[2]+
                		  " enabled="+this.gIP[fileNumber].enabled+" hintedMatch="+this.gIP[fileNumber].hintedMatch);
                 }
            	}

    			this.gIP[fileNumber].flatFieldAvailable=pixels.length>=8;
            	if (disableNoFlatfield && !this.gIP[fileNumber].flatFieldAvailable) this.gIP[fileNumber].enabled=false; // just to use old mixed data

            	int [][] shiftRotMatrix= MatchSimulatedPattern.getRemapMatrix(this.gIP[fileNumber].getUVShiftRot());
            	int [] sizeSizeExtra=setGridsWithRemap(
                		fileNumber,
                		shiftRotMatrix, // int [][] reMap,
                		pixels,
                		patternParameters);
            	numOfGridNodes+=sizeSizeExtra[0];
            	numOfGridNodes_extra+=sizeSizeExtra[1];

            	calcGridPeriod(fileNumber,true); // use _extra (out-of-pattern nodes) will be used to filter out reflections
//System.out.println ("pixelsXY["+fileNumber+"]length="+pixelsXY[fileNumber].length);
        	}
    		if (this.debugLevel>3) {
    			System.out.println("readAllGrids(), numImages="+numImages);
    			for (int n=0;n<this.gIP.length;n++) {
					System.out.println(n+": length="+this.gIP[n].pixelsXY.length);
	    			System.out.println("pixelsUV[][][0]/pixelsUV[][][1] pixelsXY[][][0]/pixelsXY[][][1]");
					for (int i=0;i<this.gIP[n].pixelsXY.length;i++){
    					System.out.println(n+":"+i+"  "+
    							this.gIP[n].pixelsUV[i][0]+"/"+
    							this.gIP[n].pixelsUV[1][1]+"  "+
    							IJ.d2s(this.gIP[n].pixelsXY[i][0], 2)+"/"+
    							IJ.d2s(this.gIP[n].pixelsXY[i][1], 2)
    					);
    				}
    			}
    		}
    		if (this.debugLevel>0) {
    			System.out.println("readAllGrids(), numImages="+numImages+", total number of grid nodes="+numOfGridNodes+", unused nodes "+numOfGridNodes_extra);
    		}
    		 // probably - do not need to verify that this.gIS is null - should do that anyway. UPDATE: no, now reading config file creates gIS
    		buildImageSets(this.gIS!=null); // with non-null just copies motors from any of the images that has it

        	return true;
        }
        /**
         * Sometimes "Process grid files" generates outliers (by 0.1..5 pixels) TODO: find the bug
         * This program replaces the "bad" ones with predicted by 8 neighbors using 2-nd order interpolation
         * @param fPixels stack of pX,pY,target-U,target-V,contrast (some bad pixels have low contrast), red,green,blue
         * @param width grid width
         * @param tolerance maximal tolerated difference between the predicted by 8 neigbors and center pixels
         * @parame maxBadNeighb - maximal number of bad cells among 8 neighbors
         * @parame gebugLevel debug level
         * @return number of fixed nodes
         * Neighbors of bad pixels can be reported bad, so they have to be re-tried with the worst removed
         */
        public int fixBadGridNodes(
        		float [][] fpixels,
        		int width,
        		double tolerance,
        		int maxBadNeighb,
        		int debugLevel,
        		String dbgTitle){
        	int debugThreshold=3;
        	double tolerance2=tolerance*tolerance;
        	double tolerance2Final=10.0*tolerance2; // final pass - fix even if the surronding are not that good
        	int [][] dirs8=   {{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}};
        	int [] dirs8Index={1,width+1,width,width-1,-1,-width-1,-width,-width+1};
        	double [] diffs2=new double [fpixels[0].length];
        	int height=diffs2.length/width;
        	for (int i=0;i<diffs2.length;i++) diffs2[i]=-1.0; // no nodes
        	double [][][] data=new double [8][3][];
        	for (int i=0;i<data.length;i++){
        		data[i][0]=new double[2];
        		data[i][1]=new double[2];
        		data[i][2]=new double[1];
        	}
        	PolynomialApproximation polynomialApproximation=new PolynomialApproximation(0); // do not report linear
        	double maxDiff2=0.0;
        	for (int y=1; y<(height-1);y++) for (int x=1;x<(width-1);x++) {
        		int index=y*width+x;
        		if (fpixels[0][index]>=0.0){
        			int numNonZero=0;
        			for (int iDir=0;iDir<dirs8.length;iDir++){
        				int index1=index+dirs8[iDir][1]*width+dirs8[iDir][0];
        				data[iDir][0][0]=dirs8[iDir][0];
        				data[iDir][0][1]=dirs8[iDir][1];
        				data[iDir][1][0]=fpixels[0][index1];
        				data[iDir][1][1]=fpixels[1][index1];
        				if ((fpixels[0][index1]<0) || (fpixels[1][index1]<0)){
        					data[iDir][2][0]=0.0;
        				} else {
        					data[iDir][2][0]=1.0;
        					numNonZero++;
        				}
        			}
        			if (numNonZero<8) continue; // should all be defined
        			double [][] coeff=polynomialApproximation.quadraticApproximation(
        					data,
        					false); // boolean forceLinear  // use linear approximation
        			if (coeff!=null) {
        				if ((coeff[0].length<6) || (coeff[1].length<6)){
        					if (debugLevel>0){
            					System.out.println("fixBadGridNodes() linear interpolate for x="+x+", y="+y);
            					for (int j=0;j<data.length;j++){
            						System.out.println(j+" "+data[j][0][0]+"/"+data[j][0][1]+" - "+data[j][1][0]+"/"+data[j][1][1]+" : "+data[j][2][0]);
            					}
            				}
        				}
        				double dx=coeff[0][coeff[0].length-1] - fpixels[0][index];
        				double dy=coeff[1][coeff[1].length-1] - fpixels[1][index];
        				diffs2[index]=dx*dx+dy*dy;
        				if (diffs2[index]>maxDiff2) maxDiff2=diffs2[index];
        			} else {
        				if (debugLevel>0){
        					System.out.println("fixBadGridNodes() failed for x="+x+", y="+y);
        				}
        			}
        		}
        	}
        	if (maxDiff2<=tolerance2) return 0; // nothing to fix
        	// here - first debug show?
        	boolean [] localWorst=new boolean[diffs2.length];
        	int numBad=0;
        	for (int i=0;i<localWorst.length;i++){
        		if (diffs2[i]<tolerance2){
        			localWorst[i]=false;
        		} else {
        			localWorst[i]=true;
        			for (int iDir=0;iDir<dirs8Index.length;iDir++) if (diffs2[i+dirs8Index[iDir]] > diffs2[i]){
        				localWorst[i]=false;
        				break;
        			}
        			if (localWorst[i]) numBad++;
        		}
        	}
        	if (numBad==0) {
				System.out.println("fixBadGridNodes() BUG - should not get here.");
        		return 0; // should not get here -
        	}
        	double [][] dbgData=null;
			if (debugLevel>debugThreshold){
				dbgData=new double[9][];
				dbgData[0]=diffs2.clone();
				dbgData[2]=dbgData[0].clone();
				for (int i=0;i< dbgData[2].length;i++) if (!localWorst[i]) dbgData[2][i]=-1.0;
//				(new showDoubleFloatArrays()).showArrays(diffs2, width, height,  "diffs2");
			}
        	// Trying to eliminate all non local worst (may that is just extra as there anot too many bad nodes)
        	int numStillBad=0;
        	for (int i=0;i<localWorst.length;i++) if (localWorst[i]){
        		for (int iDir0=0;iDir0<dirs8Index.length;iDir0++) if (diffs2[i+dirs8Index[iDir0]] > tolerance2){ // don't bother with not-so-bad
        			int index=i+dirs8Index[iDir0]; // will never be on the border as diffs2 is <=0.0 there
        			int numNonZero=0;
        			for (int iDir=0;iDir<dirs8.length;iDir++){
        				int index1=index+dirs8[iDir][1]*width+dirs8[iDir][0];
        				data[iDir][0][0]=dirs8[iDir][0];
        				data[iDir][0][1]=dirs8[iDir][1];
        				data[iDir][1][0]=fpixels[0][index1];
        				data[iDir][1][1]=fpixels[1][index1];
        				if ((data[iDir][1][0]<0) || (data[iDir][1][1]<0) || localWorst[index1]){
            				data[iDir][2][0]=0.0;
        				} else {
        					data[iDir][2][0]=1.0;
        					numNonZero++;
        				}

        			}
    				if (debugLevel>3){
    					System.out.print("+++ fixBadGridNodes() trying to fix for x="+(index%width)+", y="+(index/width)+", iDir0="+iDir0+" numNonZero="+numNonZero+" maxBadNeighb="+maxBadNeighb);
    				}

        			if (numNonZero<(data.length-maxBadNeighb-1)) continue;
        			double [][] coeff=polynomialApproximation.quadraticApproximation(
        					data,
        					false); // boolean forceLinear  // use linear approximation
        			if (coeff!=null) {
        				double dx=coeff[0][coeff[0].length-1] - fpixels[0][index];
        				double dy=coeff[1][coeff[1].length-1] - fpixels[1][index];
        				if (debugLevel>3){
        					System.out.print("fixBadGridNodes() old diffs2["+index+"]="+diffs2[index]);
        				}
        				diffs2[index]=dx*dx+dy*dy; // updated value
        				if (debugLevel>3){
        					System.out.print(" new diffs2["+index+"]="+diffs2[index]);
        				}
        				if (diffs2[index]>tolerance2) {
        					numStillBad++;
            				if (debugLevel>3){
            					System.out.print(" --- BAD");
            				}
        				} else if (debugLevel>3){
        					System.out.print(" --- GOOD");
        				}
        				if ((coeff[0].length<6) || (coeff[1].length<6)){
        					if (debugLevel>3){
            					System.out.print("fixBadGridNodes() 2 linear interpolate for x="+(index%width)+", y="+(index/width));
            					for (int j=0;j<data.length;j++){
            						System.out.println(j+" "+data[j][0][0]+"/"+data[j][0][1]+" - "+data[j][1][0]+"/"+data[j][1][1]+" : "+data[j][2][0]);
            					}
            				}
        				}
        			} else {
        				if (debugLevel>3){
        					System.out.println("fixBadGridNodes() failed for x="+(index%width)+", y="+(index/width)+", iDir0="+iDir0);
        				}
        			}
        			if (debugLevel>3) System.out.println();
        		}
        	}
        	if (numStillBad>0){
        		if (debugLevel>3){
        			System.out.println("fixBadGridNodes(): numStillBad="+numStillBad+" > 0 - probably near the border, just make sure  OK.");
        		}
        	}
			if (debugLevel>debugThreshold){
				dbgData[1]=diffs2.clone();
				for (int i=0;i< dbgData[1].length;i++) if (localWorst[i]) dbgData[1][i]=0.0;
				dbgData[3]=new double[dbgData[0].length];
				for (int i=0;i< dbgData[3].length;i++)  dbgData[3][i]=0.0;
				dbgData[4]=dbgData[3].clone();
				dbgData[5]=dbgData[3].clone();
				dbgData[6]=dbgData[3].clone();
				dbgData[7]=dbgData[3].clone();
				dbgData[8]=dbgData[3].clone();
				for (int i=0;i< dbgData[3].length;i++)  {
					dbgData[3][i]=fpixels[0][i];
					dbgData[4][i]=fpixels[1][i];
			    }
			}

// TODO - try to fix some around pixels first?

// Actually patching locally worst nodes
        	for (int index=0;index<localWorst.length;index++) if (localWorst[index]){
        		int numNonZero=0;
    			for (int iDir=0;iDir<dirs8.length;iDir++){
    				int index1=index+dirs8[iDir][1]*width+dirs8[iDir][0];
    				data[iDir][0][0]=dirs8[iDir][0];
    				data[iDir][0][1]=dirs8[iDir][1];
    				data[iDir][1][0]=fpixels[0][index1];
    				data[iDir][1][1]=fpixels[1][index1];
    				if (diffs2[index1]>tolerance2Final){ // increased tolerance for the final correction
    					data[iDir][2][0]=0.0; // do not count neighbors who are bad themselves
    				} else {
    					data[iDir][2][0]=1.0;
    					numNonZero++;
    				}
    			}
    			if (numNonZero<(data.length-maxBadNeighb)){
    				if (debugLevel>3){
    					System.out.println("fixBadGridNodes() failed x="+(index%width)+", y="+(index/width)+", number of good neighbors="+numNonZero);
    				}
    				continue; // do not fix anything
    			}
    			double [][] coeff=polynomialApproximation.quadraticApproximation(
    					data,
    					false); // boolean forceLinear  // use linear approximation
    			if (coeff!=null) {
    				if ((coeff[0].length<6) || (coeff[1].length<6)){
    					if (debugLevel>3){
        					System.out.println("fixBadGridNodes() linear interpolate for x="+(index%width)+", y="+(index/width));
        					for (int j=0;j<data.length;j++){
        						System.out.println(j+" "+data[j][0][0]+"/"+data[j][0][1]+" - "+data[j][1][0]+"/"+data[j][1][1]+" : "+data[j][2][0]);
        					}
        					for (int n=0;n<coeff.length;n++){
        						for (int j=0;j<coeff[n].length;j++){
        							System.out.print(coeff[n][j]+" ");
        						}
        						System.out.println();
        					}
        				}
    				} else if (debugLevel>3){
    					System.out.println("fixBadGridNodes() qudratic interpolate for x="+(index%width)+", y="+(index/width));
    					for (int j=0;j<data.length;j++){
    						System.out.println(j+" "+data[j][0][0]+"/"+data[j][0][1]+" - "+data[j][1][0]+"/"+data[j][1][1]+" : "+data[j][2][0]);
    					}
    					for (int n=0;n<coeff.length;n++){
    						for (int j=0;j<coeff[n].length;j++){
    							System.out.print(coeff[n][j]+" ");
    						}
    						System.out.println();
    					}
    					if (((index%width)==19) && ((index/width)==57)){
    						coeff=(new PolynomialApproximation(4)).quadraticApproximation(
    		    					data,
    		    					false);
    					}
    				}
    				fpixels[0][index]=(float) coeff[0][coeff[0].length-1];
    				fpixels[1][index]=(float) coeff[1][coeff[1].length-1];
    			} else {
    				if (debugLevel>3){
    					System.out.println("fixBadGridNodes() failed for x="+(index%width)+", y="+(index/width)+", last pass");
    				}
    			}
        	}
			if (debugLevel>debugThreshold){
				for (int i=0;i< dbgData[3].length;i++)  {
					dbgData[5][i]=fpixels[0][i];
					dbgData[6][i]=fpixels[1][i];
					dbgData[7][i]=dbgData[3][i]-fpixels[0][i];
					dbgData[8][i]=dbgData[4][i]-fpixels[1][i];
			    }

				String [] dbgTitles={"diff20","diff2Mod","localWorst", "old-X", "old-Y", "new-X", "new-Y","old-new-X","old-new-Y"};
				if (dbgTitle!=null) (new ShowDoubleFloatArrays()).showArrays(dbgData, width, height, true,  dbgTitle, dbgTitles);
			}
        	return numBad;
        }

// TODO: Move all custom image properties (including encode/decode from JP4_reader_camera) to a separate class.
// below is a duplicatie from MatchSimulatedPattern
        @Deprecated
        public double[][] getPointersXY(ImagePlus imp, int numPointers){
			   // read image info to properties (if it was not done yet - should it?
			   if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
				   JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
				   jp4_instance.decodeProperiesFromInfo(imp);
			   }
			   double [][] pointersXY=new double[numPointers][];
			   int numPointerDetected=0;
			   for (int i=0;i<pointersXY.length;i++) {
				   pointersXY[i]=null;
				   if ((imp.getProperty("POINTER_X_"+i)!=null) && (imp.getProperty("POINTER_Y_"+i)!=null)) {
					   pointersXY[i]=new double[2];
					   pointersXY[i][0]=Double.parseDouble((String) imp.getProperty("POINTER_X_"+i));
					   pointersXY[i][1]=Double.parseDouble((String) imp.getProperty("POINTER_Y_"+i));
					   numPointerDetected++;
				   }
			   }
			   if (numPointerDetected>0) return pointersXY;
			   else return null;
		   }

        public int [] getMotorPositions(ImagePlus imp, int numMotors){
        	// read image info to properties (if it was not done yet - should it?
        	if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
        		JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
        		jp4_instance.decodeProperiesFromInfo(imp);
        	}
        	int [] motorPos=new int [numMotors];
        	int numMotorsDetected=0;
        	for (int i=0;i<motorPos.length;i++) {
        		motorPos[i]=0;
        		if (imp.getProperty("MOTOR"+(i+1))!=null) {
        			motorPos[i]=Integer.parseInt((String) imp.getProperty("MOTOR"+(i+1)));
        			numMotorsDetected++;
        		}
        	}
        	if (numMotorsDetected>0) return motorPos;
        	else return null;
        }

        public int  getUsedPonters(ImagePlus imp){
        	// read image info to properties (if it was not done yet - should it?
        	if ((imp.getProperty("timestamp")==null) || (((String) imp.getProperty("timestamp")).length()==0)) {
        		JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
        		jp4_instance.decodeProperiesFromInfo(imp);
        	}
        	if (imp.getProperty("USED_POINTERS")!=null) {
        		return Integer.parseInt((String) imp.getProperty("USED_POINTERS"));
        	}
        	return 0;
        }


        // get "effective" grid period scaled for low-res (as LWIR) sensors
    	public double getEffectivePeriod(int numImg) {
    		double period = this.gIP[numImg].getGridPeriod();
    		int chn =       this.gIP[numImg].getChannel();
    		if ((this.small_sensors != null) && this.small_sensors[chn]) period /= small_period_frac;
    		return period;
    	}

    	public boolean hasSmallSensors() {
    		return small_period_frac > 0.0;
    	}
    	public boolean [] getSmallSensors() {
    		return small_sensors;
    	}
    	public boolean isSmallSensor(int numImg) {
    		if ((this.gIP != null) && (numImg >= 0) &&   (numImg < this.gIP.length) && (small_sensors != null)){
    			return small_sensors[this.gIP[numImg].getChannel()];
    		}
    		return false;
    	}


    	public double getSmallPeriodFrac() {
    		return small_period_frac;
    	}
//        public boolean [] small_sensors =    null; // set by filter grids
//        public double     small_period_frac =   0; // set by filter grids - ratio of small sensor period to large sensor period

    	// depending on camera type, return group, groups, group name
    	// camera type: eyesis26, lwir/eo (2 resolutions) , single, other
    	//getNumSubCameras()
    	public int getNumLwir() {
    		if (hasSmallSensors()) {
    			int n = 0;
    			for (int i = 0; i < small_sensors.length; i++) if (small_sensors[i]) n++;
    			return n;
    		} else {
    			return 0;
    		}
    	}
    	public int getNumEo() {
    		return getNumSubCameras() - getNumLwir();
    	}
    	public int getEo0() {
    		if (hasSmallSensors()) {
    			for (int i = 0; i < small_sensors.length; i++) if (!small_sensors[i]) return i;
    			return -1; // should not happen
    		}
    		return 0;
    	}

    	public int getLwir0() {
    		if (hasSmallSensors()) {
    			for (int i = 0; i < small_sensors.length; i++) if (small_sensors[i]) return i;
    			return -1; // should not happen
    		}
    		return -1;
    	}

    	// Get number of different subcameras for adjustments (to share adjustment types)
    	public int getSubGroups() {
    		int num_sub = getNumSubCameras();
    		if (num_sub == 1)  return 1; // single
    		if (num_sub == 26) return 3; // eyesis4pi-26
    		int n = 2;
    		if (hasSmallSensors()) {
    			if (getNumLwir() > 1) n++;
    			if (getNumEo() > 1) n++;
    		}
    		return n;
    	}

    	// Get subcamera adjustment group
    	// May be modified to use other subcamera as zero (eyesis in the middle row)
    	public int getSubGroup(int chn) {
    		int groups = getSubGroups();
    		if (groups == 1) return 0; // single camera
    		int num_sub = getNumSubCameras();
    		if (num_sub == 26) { // eyesis4pi-26
    			if (chn == 0) return 0;
    			if (chn < 24) return 1;
    			return 2;
    		}
			int n = 0;
    		if (hasSmallSensors()) {
    			if (small_sensors[chn]) { // current is LWIR
    				n = 1;
    				if (getNumEo() > 1) n = 2;
    				if (chn != getLwir0()) n++;
    			} else {  // current is EO
    				n = 0;
    				if (chn != getEo0()) n++;
    			}
    		}
			return n;
    	}
    	// check if the channel is the first in group that represents settings
    	public boolean firstInGroup(int chn) {
    		if (chn >= 24) return (chn == 24);
    		if ((chn == getEo0()) || (chn == getLwir0())) return true;
    		int [] num = {0,0};
    		for (int i = 0; i < chn; i++) {
    			if ((small_sensors != null) && small_sensors[i]) num[1]++;
    			else num[0]++;
    		}
    		if ((small_sensors != null) && small_sensors[chn]) return (num[1] == 1);
    		else                                               return (num[0] == 1);
    	}

    	// can be a clone of the previous channel (same value)
    	public int sourceToCopy(int chn) {
    		if ((chn == getEo0()) || (chn == getLwir0())) return -1; // no prototype to copy
    		return chn -1; // copy from previous
    	}



    	public String getSubName(int chn, boolean full) {
    		int groups = getSubGroups();
    		if (groups == 1) return full?"Single subcamera":"sub"; // single camera
    		int num_sub = getNumSubCameras();
    		if (num_sub == 26) { // eyesis4pi-26
    			if (chn == 0) return full?"Camera head subcamera 0":"sub-head-0";
    			if (chn < 24) return full?"Camera head other subcameras":"sub-head-other";
    			if (chn == 24) return full?"Camera bottom 0":"sub-bottom-0";
    			return                full?"Camera bottom other":"sub-bottom-other";
    		}
    		if (hasSmallSensors()) {
    			if (small_sensors[chn]) { // current is LWIR
    				if (chn != getLwir0()) return full?"Subcamera LWIR other":"sub-lwir-other";
    				if (getNumLwir() > 1) return  full?"Subcamera LWIR 0":"sub-lwir0";
    				return                        full?"Subcamera LWIR":"sub-lwir";
    			} else {  // current is EO
    				if (chn != getEo0()) return full?"Subcamera EO other":"sub-eo-other";
    				if (getNumEo() > 1) return  full?"Subcamera EO 0":"sub-eo0";
    				return                        full?"Subcamera EO":"sub-eo";
    			}
    		}
			if (chn != getEo0()) return full?"Subcamera other":"sub-other";
			return                        full?"Subcamera 0":"sub0";
    	}

    	// get subcamera index this one should copy parameters from (-1 if it is first in group
    	public int masterSub(int chn) {
    		if (firstInGroup(chn)) return -1;
    		for (int mchn = chn -1; mchn >= 0; mchn--) {
    			if (firstInGroup(mchn)) return mchn;
    		}
    		return -1; // should never get here
    	}



        public int getImageNumPoints(int numImg){
        	return this.gIP[numImg].pixelsUV.length;
        }

        public void initPars(int numImages, int numPars) {
        	this.pars=new double [numImages][numPars];
        	for (int i=0;i<numImages;i++) for (int j=0;j<numPars;j++) this.pars[i][j]=Double.NaN;
        }
        public double getParameterValue(int numImg, int numPar){
        	if ((numImg<0) || (numImg>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+numImg;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if ((numPar<0) || (numPar>=this.pars[numImg].length)) {
        		String msg="There are only "+this.pars[numImg].length+" parameters defined, requested #"+numPar;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	double par=(this.gIP[numImg].gridImageSet!=null)?this.gIP[numImg].gridImageSet.getParameterValue(numPar):Double.NaN;
        	if (Double.isNaN(par)) par=this.pars[numImg][numPar];
        	return par;
        }
        public void setParameterValue(int numImg, int numPar, double value, boolean updateEstimated){
        	if ((numImg<0) || (numImg>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+numImg;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if ((numPar<0) || (numPar>=this.pars[numImg].length)) {
        		String msg="There are only "+this.pars[numImg].length+" parameters defined, requested #"+numPar;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	this.pars[numImg][numPar]=value;
        	if (this.gIP[numImg].gridImageSet!=null) this.gIP[numImg].gridImageSet.setParameterValue(numPar,value,updateEstimated);
        }

        public void setParameters(double [] parameters, int numImg){
        	if ((numImg<0) || (numImg>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+numImg;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	this.pars[numImg]=parameters.clone();
        	if (this.gIP[numImg].gridImageSet!=null) this.gIP[numImg].gridImageSet.updateSetFromParameterVector(parameters);
        }

        public int getParametersLength(int numImg){
        	return this.pars[numImg].length;
        }
        public double [] getParameters(int numImg){
        	if ((numImg<0) || (numImg>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+numImg;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	double [] parameters=this.pars[numImg].clone();
        	if (this.gIP[numImg].gridImageSet!=null) this.gIP[numImg].gridImageSet.updateParameterVectorFromSet(parameters);
        	return parameters;
        }
        public int [] getUVShiftRot(int numImg){
        	return this.gIP[numImg].getUVShiftRot();
        }
        public GridImageParameters getGridImageParameters(int numImg){
        	return this.gIP[numImg];
        }

        // next is just for goniometer - use elevation and heading for cartesian mode?
        public double [] getHeadEl(int imgNum){ // get sensor heading +(azimuth) and elevation
        	if ((imgNum<0) || (imgNum>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+imgNum;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	double [] headel={this.pars[imgNum][index_heading],this.pars[imgNum][index_elevation]};
        	if (!isCartesian()) {
        		headel[0] += this.pars[imgNum][index_azimuth];
        	}
        	System.out.println("getHeadEl("+imgNum+") "+isCartesian()+" -> "+headel[0]+"/"+ headel[1]+", "+this.pars[imgNum][index_azimuth]+","+
        			this.pars[imgNum][index_heading]+", "+this.pars[imgNum][index_elevation]);
        	return headel;
        }
        // set goniometer horizontal axis angle and goniometer axial angles in all images
        public void setGHGA(double gh, double ga){
        	for (int imgNum=0;imgNum<this.pars.length;imgNum++) setGHGA( imgNum, gh,ga);
        }
        public void setGHGA(int imgNum, double gh, double ga){
        	setGH(imgNum, gh);
        	setGA(imgNum, ga);
        }
        public void setGH(int numImg, double gh){
        	this.pars[numImg][index_gh]=gh;
        	if (this.gIP[numImg].gridImageSet!=null) this.gIP[numImg].gridImageSet.goniometerTilt= gh;
        }
        public void setGA(int numImg,  double ga){
        	this.pars[numImg][index_ga]=ga;
        	if (this.gIP[numImg].gridImageSet!=null) this.gIP[numImg].gridImageSet.goniometerAxial= ga;
        }
        public double getGH(int numImg){
        	if (this.gIP[numImg].gridImageSet!=null) return this.gIP[numImg].gridImageSet.goniometerTilt;
        	return this.pars[numImg][index_gh];
        }

        public double getGA(int numImg){
        	if (this.gIP[numImg].gridImageSet!=null) return this.gIP[numImg].gridImageSet.goniometerAxial;
        	return this.pars[numImg][index_ga];
        }

        public void setParameters(double [] parameters, int numImg, boolean[] mask){
        	if ((numImg<0) || (numImg>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+numImg;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if ((this.pars[numImg].length!=parameters.length) || (this.pars[numImg].length!=mask.length)) {
        		String msg="Vector lengths for image #"+numImg+
        		" mismatch: this.pars["+numImg+"].length="+this.pars[numImg].length+
        		" parameters.length="+parameters.length+
        		" mask.length="+mask.length;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	for (int i=0;i<mask.length;i++) if (mask[i])this.pars[numImg][i]=parameters[i];
        	if (this.gIP[numImg].gridImageSet!=null) this.gIP[numImg].gridImageSet.updateSetFromParameterVector(parameters,mask);

        }

        public void setIntrinsicParameters(double [] parameters, int num){
        	if ((num<0) || (num>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (this.pars[num].length!=parameters.length) {
        		String msg="Vector lengths for image #"+num+
        		" mismatch: this.pars["+num+"].length="+this.pars[num].length+
        		" parameters.length="+parameters.length;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	for (int i=0;i<parameters.length;i++) if (isIntrinsicParameter(i))this.pars[num][i]=parameters[i];
        	// no need to update image sets
        }
        public void setSubcameraParameters(double [] parameters, int num){
        	if ((num<0) || (num>=this.pars.length)) {
        		String msg="There are only "+this.pars.length+" images defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (this.pars[num].length!=parameters.length) {
        		String msg="Vector lengths for image #"+num+
        		" mismatch: this.pars["+num+"].length="+this.pars[num].length+
        		" parameters.length="+parameters.length;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	for (int i=0;i<parameters.length;i++) if (isSubcameraParameter(i))this.pars[num][i]=parameters[i];
        	// no need to update image sets
        }


        public String getParameterName(int num){
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return descrField(num,0);

        }
        public String getParameterDescription(int num){
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return descrField(num,1);

        }
        public String getParameterUnits(int num){
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return descrField(num,2);
        }
        public boolean isSubcameraParameter(int num){
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return (descrField(num,3).equals("S"));

        }
        public boolean isLocationParameter(int num){ //X,Y or Z location of the camera
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return (descrField(num,3).equals("T"));
        }

        public boolean isOrientationParameter(int num){ //one of the 2 goniometer orientation angles
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return (descrField(num,3).equals("R"));
        }

        public boolean isIntrinsicParameter(int num){ // updated from image calibration file
        	if ((num<0) || (num>=this.parameterDescriptions.length)) {
        		String msg="There are only "+this.parameterDescriptions.length+" parameters defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return (descrField(num,4).equals("I"));

        }
        public String getImagePath(int num) {
        	if ((num<0) || (num>=this.gIP.length)) {
        		String msg="There are only "+this.gIP.length+" images defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return this.gIP[num].path;
        }
        public int getImageSubcamera(int num) {
        	if ((num<0) || (num>=this.gIP.length)) {
        		String msg="There are only "+this.gIP.length+" images defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return this.gIP[num].channel;
        }
        public int getImageStation(int num) {
        	if ((num<0) || (num>=this.gIP.length)) {
        		String msg="There are only "+this.gIP.length+" images defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return this.gIP[num].getStationNumber();
        }
        public double getImageTimestamp(int num) {
        	if ((num<0) || (num>=this.gIP.length)) {
        		String msg="There are only "+this.gIP.length+" images defined, requested #"+num;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return this.gIP[num].timestamp;
        }
        public int getNumImages() {
        	if (this.gIP == null) return 0;
        	return this.gIP.length;
        }
        public int getNumParameters() {
        	if (this.parameterDescriptions == null) return 0;
        	return this.parameterDescriptions.length;
        }
        public int getNumSubCameras() {
        	return this.numSubCameras;
        }
        public int getNumNodes(int num, boolean use_extra) {
        	if ((this.gIP == null) || (num >= this.gIP.length) || (this.gIP[num].pixelsXY == null)) return 0;
        	int len = this.gIP[num].pixelsXY.length;
        	if (use_extra && (this.gIP[num].pixelsXY_extra != null)) len += this.gIP[num].pixelsXY_extra.length;
        	return len;
        }

        public int getImageNumber(int set_number, int sub_number) {
        	if ((this.gIS== null) || (this.gIS.length <= set_number) || (this.gIS[set_number] == null)) return -1;
        	if (sub_number >= getNumChannels()) return -1;
        	if ((this.gIS[set_number].imageSet == null) || (this.gIS[set_number].imageSet[sub_number] == null)) return -1;
        	return this.gIS[set_number].imageSet[sub_number].imgNumber;
        }

        public int getImageSet(int img_number) {
        	if ((this.gIP== null) || (this.gIP.length <= img_number) || (img_number <  0) ||(this.gIP[img_number] == null)) return -1;
        	return this.gIP[img_number].setNumber;
        }


        /**
         *
         * @param imgNumber number of grid image to edit parameters (location, distortion) for
         * @return <2 - canceled, -1 - done, els - number of the next image to edit
         */
        public int editImageParameters(int imgNumber){
        	if ((this.gIP==null) || (imgNumber<0) ||(imgNumber>=this.gIP.length)) return -3;
			int sub=getImageSubcamera(imgNumber);
       	    String sTS=IJ.d2s(getImageTimestamp(imgNumber),6);
    		GenericDialog gd = new GenericDialog("Manually editing per-image parameters, timestamp="+sTS+
    				", subchannel-"+sub+" "+getImagePath(imgNumber));
    	    for (int i=0;i<getNumParameters();i++){
    	    	gd.addNumericField(
    	    			i+": "+getParameterDescription(i)+"["+ getParameterName(i)+"] "+
    	    			(isSubcameraParameter(i)?"S ":"  "),
//    	    			this.pars[imgNumber][i],5,10, getParameterUnits(i));
    	    			this.getParameterValue(imgNumber,i),5,10, getParameterUnits(i));
    	    }
    	    gd.addNumericField("Next image to edit (0.."+this.pars.length+", -1 - none) ", imgNumber+1,0);
   	        gd.enableYesNoCancel("OK", "Done");
    	    WindowTools.addScrollBars(gd);
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return -2;
    	    for (int i=0;i<getNumParameters();i++){
//    	    	this.pars[imgNumber][i]= gd.getNextNumber();
    	    	this.setParameterValue(imgNumber,i, gd.getNextNumber(),true);
    	    }
    	    imgNumber= (int) gd.getNextNumber();
    	    if ((imgNumber<0) || (imgNumber>=getNumImages())) return -1;
    		if (!gd.wasOKed()) return -1; // pressed Done (no need to ask for the next number)
            return imgNumber;
        }
        @Deprecated
        public void setMaskFromImageStack(String path){ // can not work with different size senors
    		Opener opener=new Opener();
			if (this.debugLevel>1) System.out.println("Opening "+path+" as a stack of sensor masks");
			ImagePlus imp=opener.openImage("", path);
        	if (imp==null) {
        		String msg="Failed to read sensors mask file "+path;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
        	if (imp.getProperty("shrinkGridForMask")!=null)
        		eyesisCameraParameters.shrinkGridForMask=Integer.parseInt((String) imp.getProperty("shrinkGridForMask"));
            	if (imp.getProperty("maskBlurSigma")!=null)
            		eyesisCameraParameters.maskBlurSigma=Double.parseDouble((String) imp.getProperty("maskBlurSigma"));
            	if (imp.getProperty("decimateMasks")!=null)
            		eyesisCameraParameters.setDecimateMasks(Integer.parseInt((String) imp.getProperty("decimateMasks")));
            	if (imp.getProperty("sensorWidth")!=null)
            		eyesisCameraParameters.setSensorWidth(Integer.parseInt((String) imp.getProperty("sensorWidth")));
            	if (imp.getProperty("sensorHeight")!=null)
            		eyesisCameraParameters.setSensorHeight(Integer.parseInt((String) imp.getProperty("sensorHeight")));
        	setMaskFromImageStack(imp);
        }
        /**
         * Find number of channels in this camera
         * @return maximal number of channel used plus one
         */
        public int getNumChannels(){
        	int nChn=-1;
        	for (int i=0;i<this.gIP.length;i++) if (this.gIP[i].channel>nChn) nChn=this.gIP[i].channel;
        	return nChn+1;
        }


        public double getMask(int chnNum, double px, double py){
        	int width= eyesisCameraParameters.getSensorWidth(chnNum)/eyesisCameraParameters.getDecimateMasks(chnNum);
        	int height=eyesisCameraParameters.getSensorHeight(chnNum)/eyesisCameraParameters.getDecimateMasks(chnNum);
        	int iPX= ((int) Math.round(px))/eyesisCameraParameters.getDecimateMasks(chnNum);
        	int iPY= ((int) Math.round(py))/eyesisCameraParameters.getDecimateMasks(chnNum);
        	if ((iPX<0) || (iPY<0) || (iPX>=width) || (iPY>=height)) return 0.0;
        	if ((this.sensorMasks==null) || (this.sensorMasks[chnNum]==null)) return 1.0;
        	return this.sensorMasks[chnNum][iPY*width+iPX];
        }
        @Deprecated
        public double getMask(double[] mask, double px, double py){ // problems with different size sensors
        	if (mask==null) return 0;
        	int width= eyesisCameraParameters.getSensorWidth()/eyesisCameraParameters.getDecimateMasks();
        	int height=eyesisCameraParameters.getSensorHeight()/eyesisCameraParameters.getDecimateMasks();
        	int iPX= ((int) Math.round(px))/eyesisCameraParameters.getDecimateMasks();
        	int iPY= ((int) Math.round(py))/eyesisCameraParameters.getDecimateMasks();
        	if ((iPX<0) || (iPY<0) || (iPX>=width) || (iPY>=height)) return 0.0;
        	return mask[iPY*width+iPX];  // null ponter
        }


        @Deprecated
        public void setMaskFromImageStack(ImagePlus imp){
        	if (imp == null){
        		String msg="sensors mask image is null";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (imp.getProperty("decimateMasks")!=null) eyesisCameraParameters.setDecimateMasks(Integer.parseInt((String) imp.getProperty("decimateMasks")));
        	eyesisCameraParameters.setSensorWidth(imp.getWidth()*eyesisCameraParameters.getDecimateMasks());
        	eyesisCameraParameters.setSensorHeight(imp.getHeight()*eyesisCameraParameters.getDecimateMasks());
        	if (imp.getProperty("sensorWidth")!=null)	eyesisCameraParameters.setSensorWidth(Integer.parseInt((String) imp.getProperty("sensorWidth")));
        	if (imp.getProperty("sensorHeight")!=null)  eyesisCameraParameters.setSensorHeight(Integer.parseInt((String) imp.getProperty("sensorHeight")));

    		if (this.sensorMasks==null) {
    			this.sensorMasks=new double[getNumChannels()][];
    			for (int i=0;i<this.sensorMasks.length;i++) this.sensorMasks[i]=null;
    		}
    		int numChannels=imp.getStackSize();
    		float [][] pixels =new float[numChannels][];
    		if (numChannels==1){
    			pixels[0]=(float[]) imp.getProcessor().getPixels();
    		} else {
        		ImageStack stack = imp.getStack();
            	if (stack==null) {
            		String msg="Expected a image stack with masks";
            		IJ.showMessage("Error",msg);
            		throw new IllegalArgumentException (msg);
            	}
            	for (int i=0;i<numChannels;i++) pixels[i]= (float[]) stack.getPixels(i+1);
    		}
    		for (int numChn=0;(numChn<numChannels) && (numChn<this.sensorMasks.length);numChn++){
    			//Make sure masks contain non-zero (>0.0) pixels, otherwise skip those
        		boolean defined=false;
        		for (int i=0;i<pixels[numChn].length;i++) if (pixels[numChn][i]>0.0){
        			defined=true;
        			break;
        		}
    			if (defined) {
    				this.sensorMasks[numChn]=new double [pixels[numChn].length];
    				for (int i=0;i<this.sensorMasks[numChn].length;i++) this.sensorMasks[numChn][i]=pixels[numChn][i];
    			}
    		}
        }
        public ImagePlus saveMaskAsImageStack(String title, String path){
        	ImagePlus imp=getMaskAsImageStack(title);
        	if (imp==null) return null;
	   				FileSaver fs=new FileSaver(imp);
	   				if (updateStatus) IJ.showStatus("Saving masks "+path);
	   				if (this.debugLevel>0) System.out.println("Saving masks "+path);
	   				if (imp.getStackSize()>1)
	   					fs.saveAsTiffStack(path);
	   				else
	   					fs.saveAsTiff(path);
        	return imp;
        }

        @Deprecated
        public ImagePlus getMaskAsImageStack(String title){
        	if (this.sensorMasks==null){
        		String msg="Sensor mask array does not exist, nothing to convert";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	int width= eyesisCameraParameters.getSensorWidth()/eyesisCameraParameters.getDecimateMasks();
        	int height=eyesisCameraParameters.getSensorHeight()/eyesisCameraParameters.getDecimateMasks();
        	float [][]pixels=new float [getNumChannels()][width*height];
        	ImagePlus imp=null;
        	for (int numChn=0;numChn<getNumChannels();numChn++){
        		if (this.sensorMasks[numChn]==null) for (int i=0;i<pixels[numChn].length;i++)pixels[numChn][i]=0.0F;
        		else for (int i=0;i<pixels[numChn].length;i++)pixels[numChn][i]=(float) this.sensorMasks[numChn][i];
        	}
        	if (this.sensorMasks.length>0){
        		ImageStack stack=new ImageStack(width,height);
        		for (int numChn=0;numChn<pixels.length;numChn++)  stack.addSlice("chn-"+numChn,    pixels[numChn]);
        		imp = new ImagePlus(title, stack);
        	} else {
        		ImageProcessor  ip =new FloatProcessor(width,height);
        		ip.setPixels(pixels[0]);
        		imp=new ImagePlus(title, ip);
        	}
// TODO: add more properties here (MAC+channel)? preserve other properties?
        	imp.setProperty("sensorWidth", ""+eyesisCameraParameters.getSensorWidth());
        	imp.setProperty("sensorHeight", ""+eyesisCameraParameters.getSensorHeight());
        	imp.setProperty("shrinkGridForMask", ""+eyesisCameraParameters.shrinkGridForMask);
        	imp.setProperty("maskBlurSigma", ""+eyesisCameraParameters.maskBlurSigma);
        	imp.setProperty("decimateMasks", ""+eyesisCameraParameters.getDecimateMasks());

        	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
        	imp.getProcessor().resetMinAndMax();
        	return imp;
        }
        /**
         * Generate low-vignetting sensor mask for flat-field calculation
         * @param sensorMask sensor mask, decimated array
         * @param width  sensor width, pixels
         * @param height sensor height, pixels
         * @param shrink shrink sensor mask by this amount (sensor, non-decimated pixels)
         * @param radius radial mask - zero if farther than radius, 0.5*(cos(pi*r/radius)+1.0) if less
         * @param minimalAlpha - zero mask below this threshold
         * @return returns arrray with the same size as sensorMask that corresponds to low-vignetting areas of the sensor/lens
         */

        public double [] nonVignettedMask(
        		double [] sensorMask,
        		int width,
        		int height,
        		double x0,     // lens center X (sensor, non-decimated pix)
        		double y0,     // lens center Y (sensor, non-decimated pix)
        		double shrink,
        		double radius,
        		double minimalAlpha){

        	int decimate= (int) Math.round(Math.sqrt(width*height/sensorMask.length));
        	int dcmWidth= width/decimate;
        	int dcmHeight=height/decimate;
        	double [] mask= sensorMask.clone();
        	if (shrink>0){
        		(new DoubleGaussianBlur() ).blurDouble(mask, dcmWidth, dcmHeight, shrink/decimate, shrink/decimate, 0.01);
        		for (int i=0;i<mask.length;i++){
        			double d=2*(mask[i]-0.5);
        			mask[i]=(d>0)?(d*d):(0.0);
        		}
        	}
        	if (radius>0.0){
        		int index=0;
        		for (int iy=0; iy<dcmHeight;iy++) for (int ix=0; ix<dcmWidth;ix++){
        			double r=Math.sqrt((iy*decimate-y0)*(iy*decimate-y0)+(ix*decimate-x0)*(ix*decimate-x0))/radius;
        			double k=(r>1.0)?0.0:(0.5*(Math.cos(Math.PI*r)+1.0));
        			mask[index++]*=k;
        		}
        	}
        	if (minimalAlpha>0.0) for (int i=0;i<mask.length;i++) if (mask[i]<minimalAlpha) mask[i]=0.0;
        	return mask;
        }
        @Deprecated
        public double [][] calculateSensorMasksOld() {
        	return calculateSensorMasks(
        			eyesisCameraParameters.getDecimateMasks(),
        			eyesisCameraParameters.getSensorWidth(),
        			eyesisCameraParameters.getSensorHeight(),
        			eyesisCameraParameters.shrinkGridForMask,
        			eyesisCameraParameters.maskBlurSigma);
        }

        public double [][] calculateSensorMasks() {
        	return calculateSensorMasks(
        			eyesisCameraParameters.shrinkGridForMask,
        			eyesisCameraParameters.maskBlurSigma);
        }

        /**
         *
         * @param width image width, in pixels (pixel X coordinates are between 0 and width-1, inclusive)
         * @param height image height, in pixels (pixel Y coordinates are between 0 and height-1, inclusive)
         * @param shrinkGridForMask shrink detected grids by this number of nodes in each direction before bluring
         * @param sigmaUV Gaussian sigma for blurring of the sensor mask (if negative - in grid inter-node distances)
         * @return array of pixel arrays (or nulls) for each camera subchannel (also keeps it in the class instance)
         */
        @Deprecated
        public double [][] calculateSensorMasks( int decimate, int width, int height, int shrinkGridForMask, double sigmaUV) {
        	int dWidth=  (width -1)/decimate+1;
        	int dHeight= (height-1)/decimate+1;
        	int numChannels=getNumChannels();
        	this.sensorMasks=new double [numChannels][];
        	DoubleGaussianBlur gb=new DoubleGaussianBlur();
        	if ((this.debugLevel>1) && (SDFA_INSTANCE==null)) SDFA_INSTANCE=new ShowDoubleFloatArrays();
			if (this.debugLevel>2)System.out.println("calculateSensorMasks("+width+","+height+","+shrinkGridForMask+","+sigmaUV+")");
        	for (int chNum=0;chNum<numChannels; chNum++){
        		this.sensorMasks[chNum]=new double[dWidth*dHeight];
        		for (int i=0;i<this.sensorMasks[chNum].length;i++) this.sensorMasks[chNum][i]=0.0;
        		double rAverage=0.0;
        		double rAverageNum=0.0;
        		for (int imgNum=0;imgNum<this.gIP.length;imgNum++) if (this.gIP[imgNum].channel==chNum){ // image is for this this channel
        	        double [][] preMask=preCalculateSingleImageMask(imgNum, decimate, width, height, shrinkGridForMask);
        	        if (preMask==null) continue; //nothing in this channel
            		rAverage+=preMask[0][0];
            		rAverageNum+=preMask[0][1];
       			    for (int i=0;i<this.sensorMasks[chNum].length;i++) if (preMask[1][i]>0.0) this.sensorMasks[chNum][i]=1.0;
        		}
        		if (rAverageNum==0.0) continue; // nothing to blur/process for this channel
        		rAverage/=rAverageNum; // average distance to the fartherst node from the current
        		double      sigma=sigmaUV;
        		if(sigma<0) sigma*=-rAverage;
        		gb.blurDouble(this.sensorMasks[chNum], dWidth, dHeight, sigma/decimate, sigma/decimate, 0.01);

        		//    this.sensorMasks[chNum] now contains 0.0/1.0 mask. Blur it
        		//	    		gb.blurDouble(pointedBayer[bayerR], halfWidth, halfHeight, this.lowpassSigma, this.lowpassSigma, 0.01);
        		//	    		if (debugLevel>2) sdfra_instance.showArrays(pointedBayer[bayerR].clone(), halfWidth, halfHeight, title+"-smooth");

        	}
        	return this.sensorMasks;
        }

        public double [][] calculateSensorMasks(int shrinkGridForMask, double sigmaUV) {
        	int numChannels=getNumChannels();
        	this.sensorMasks=new double [numChannels][];
        	DoubleGaussianBlur gb=new DoubleGaussianBlur();
        	if ((this.debugLevel>1) && (SDFA_INSTANCE==null)) SDFA_INSTANCE=new ShowDoubleFloatArrays();
			if (this.debugLevel>2)System.out.println("calculateSensorMasks("+shrinkGridForMask+","+sigmaUV+")");
        	for (int chNum=0;chNum<numChannels; chNum++){
        		int decimate = eyesisCameraParameters.getDecimateMasks(chNum);
        		int width = eyesisCameraParameters.getSensorWidth(chNum);
        		int height = eyesisCameraParameters.getSensorHeight(chNum);
            	int dWidth=  (width -1)/decimate+1;
            	int dHeight= (height-1)/decimate+1;

        		this.sensorMasks[chNum]=new double[dWidth*dHeight];
        		for (int i=0;i<this.sensorMasks[chNum].length;i++) this.sensorMasks[chNum][i]=0.0;
        		double rAverage=0.0;
        		double rAverageNum=0.0;
        		for (int imgNum=0;imgNum<this.gIP.length;imgNum++) if (this.gIP[imgNum].channel==chNum){ // image is for this this channel
        	        double [][] preMask=preCalculateSingleImageMask(imgNum, decimate, width, height, shrinkGridForMask);
        	        if (preMask==null) continue; //nothing in this channel
            		rAverage+=preMask[0][0];
            		rAverageNum+=preMask[0][1];
       			    for (int i=0;i<this.sensorMasks[chNum].length;i++) if (preMask[1][i]>0.0) this.sensorMasks[chNum][i]=1.0;
        		}
        		if (rAverageNum==0.0) continue; // nothing to blur/process for this channel
        		rAverage/=rAverageNum; // average distance to the fartherst node from the current
        		double      sigma=sigmaUV;
        		if(sigma<0) sigma*=-rAverage;
        		gb.blurDouble(this.sensorMasks[chNum], dWidth, dHeight, sigma/decimate, sigma/decimate, 0.01);
        	}
        	return this.sensorMasks;
        }





        /**
         * Create round mask inside the actual one, with the provided center. Blur result with the same sigma as original
         * @param chn sensor number
         * @param xCenter X of the center (before decimation)
         * @param yCenter Y of the center (before decimation)
         * @return this channel mask, also sets the round mask instead of the original
         */
        public double [] roundOffMask(int chn, double xCenter, double yCenter){
        	int dWidth=  (eyesisCameraParameters.getSensorWidth(chn) -1)/eyesisCameraParameters.getDecimateMasks(chn)+1;
        	int dHeight= (eyesisCameraParameters.getSensorHeight(chn)-1)/eyesisCameraParameters.getDecimateMasks(chn)+1;
        	DoubleGaussianBlur gb=new DoubleGaussianBlur();

        	int iXC=(int) Math.round(xCenter/eyesisCameraParameters.getDecimateMasks(chn));
        	int iYC=(int) Math.round(yCenter/eyesisCameraParameters.getDecimateMasks(chn));
        	double r0=iXC;
        	r0 = Math.min(r0, iYC);
        	r0 = Math.min(r0, dWidth - iXC);
        	r0 = Math.min(r0, dHeight - iYC);
        	int ir02=(int) Math.round(r0*r0);
        	System.out.println("iXC="+iXC);
        	System.out.println("iYC="+iYC);
        	System.out.println("initial ir02="+ir02+ "("+Math.sqrt(ir02)+")");
        	for (int i = 0; i < this.sensorMasks[chn].length; i++) if (this.sensorMasks[chn][i]<0.5) {
        		int ix = (i % dWidth) - iXC;
        		int iy = (i / dWidth) - iYC;
        		int ir2=ix*ix + iy*iy;
        		if (ir2 < ir02) ir02 = ir2;
        	}
        	System.out.println("second ir02="+ir02+ "("+Math.sqrt(ir02)+")");
        	double [] mask= new double[this.sensorMasks[chn].length];
        	for (int i = 0; i < mask.length; i++) {
        		int ix = (i % dWidth) - iXC;
        		int iy = (i / dWidth) - iYC;
        		mask[i]=((ix*ix + iy*iy) > ir02)?0.0:1.0;
        	}
        	// blur result
        	double [][] preMask=preCalculateSingleImageMask(chn,
        			eyesisCameraParameters.getDecimateMasks(chn),
        			eyesisCameraParameters.getSensorWidth(chn),
        			eyesisCameraParameters.getSensorHeight(chn),
        			eyesisCameraParameters.shrinkGridForMask);

        	if (preMask==null) return null; //nothing in this channel
        	double rAverage=preMask[0][0];
        	double rAverageNum=preMask[0][1];
        	if (rAverageNum==0.0) return null; // nothing to blur/process for this channel
        	rAverage/=rAverageNum; // average distance to the fartherst node from the current
        	double      sigma = eyesisCameraParameters.maskBlurSigma;
        	if(sigma<0) sigma*=-rAverage;
        	gb.blurDouble(mask, dWidth, dHeight, sigma/eyesisCameraParameters.getDecimateMasks(chn), sigma/eyesisCameraParameters.getDecimateMasks(chn), 0.01);
        	for (int i=0;i < mask.length;i++){
        		this.sensorMasks[chn][i] = Math.min(this.sensorMasks[chn][i],mask[i]);
        	}
        	return this.sensorMasks[chn];
        }


        public double [] calculateImageGridMask(int imgNum) {
        	int chn = this.gIP[imgNum].channel; // getChannel()
        	return calculateImageGridMask(
        			imgNum,
        			eyesisCameraParameters.getDecimateMasks(chn),
        			eyesisCameraParameters.getSensorWidth(chn),
        			eyesisCameraParameters.getSensorHeight(chn),
        			eyesisCameraParameters.shrinkGridForMask,
        			eyesisCameraParameters.maskBlurSigma);
        }

        public double [] calculateImageGridMask(int imgNum, int decimate, int width, int height, int shrinkGridForMask, double sigmaUV) {
        	int dWidth=  (width -1)/decimate+1;
        	int dHeight= (height-1)/decimate+1;
        	DoubleGaussianBlur gb=new DoubleGaussianBlur();
        	if ((this.debugLevel>1) && (SDFA_INSTANCE==null)) SDFA_INSTANCE=new ShowDoubleFloatArrays();
        	if (this.debugLevel>2)System.out.println("calculateSensorMasks("+width+","+height+","+shrinkGridForMask+","+sigmaUV+")");
        	double [][] preMask=preCalculateSingleImageMask(imgNum, decimate, width, height, shrinkGridForMask);
        	if (preMask==null) return null; //nothing in this channel
        	double rAverage=preMask[0][0];
        	double rAverageNum=preMask[0][1];
        	if (rAverageNum==0.0) return null; // nothing to blur/process for this channel
        	rAverage/=rAverageNum; // average distance to the fartherst node from the current
        	double      sigma=sigmaUV;
        	if(sigma<0) sigma*=-rAverage;
// old version, trying new - will influence all sensor masks!!
//        	gb.blurDouble(preMask[1], dWidth, dHeight, sigma/decimate, sigma/decimate, 0.01);
        	double [] mask0=preMask[1].clone();
        	gb.blurDouble(preMask[1], dWidth, dHeight, sigma/decimate, sigma/decimate, 0.01);
        	for (int i=0;i<preMask[1].length;i++){
				double d=2.0*(preMask[1][i]-0.5);
				preMask[1][i]=((mask0[i]>0) && (d>0))?(d*d):0.0;
        	}
        	return preMask[1];
        }



        /**
         *
         * @param imgNum number of image to process
         * @param decimate - reduce image resolution for the mask
         * @param width - image width (actual will be divided by decimate
         * @param height- image height (actual will be divided by decimate
         * @param shrinkGridForMask shrink defined grid before bluring
         * @return array of 2 rows - [0] has just rAverage and rAverageNum for the average radius of the grid [1] - mask (1.0/0.0)
         *         or null if there are no grid nodes at all;
         */
        public double [][] preCalculateSingleImageMask(
        		int imgNum,
        		int decimate,
        		int width,
        		int height,
        		int shrinkGridForMask){
        	if (!this.gIP[imgNum].enabled) return null; // this image is disabled, ignore it
        	int [][] dirs={{-1,0},{1,0},{0,-1},{0,1}};
        	double rAverage=0.0;
        	double rAverageNum=0.0;
        	int i;
        	int dWidth=  (width -1)/decimate+1;
        	int dHeight= (height-1)/decimate+1;
        	double [] mask = new double [dWidth*dHeight];
        	boolean hasGrid=false;
        	for (i=0;i<this.gIP[imgNum].pixelsXY.length;i++) if ((this.gIP[imgNum].pixelsXY[i]!=null) &&(this.gIP[imgNum].pixelsXY[i][0]>=0)) {
        		hasGrid=true;
        		break;
        	}
        	if (!hasGrid) return null; // image has no grid nodes
        	int minU=this.gIP[imgNum].pixelsUV[i][0];
        	int minV=this.gIP[imgNum].pixelsUV[i][1];
        	int maxU=minU;
        	int maxV=minV;
        	for (;i<this.gIP[imgNum].pixelsXY.length;i++) if ((this.gIP[imgNum].pixelsXY[i]!=null) &&(this.gIP[imgNum].pixelsXY[i][0]>=0)) {
        		if (this.gIP[imgNum].pixelsUV[i][0]<minU) minU=this.gIP[imgNum].pixelsUV[i][0];
        		if (this.gIP[imgNum].pixelsUV[i][1]<minV) minV=this.gIP[imgNum].pixelsUV[i][1];
        		if (this.gIP[imgNum].pixelsUV[i][0]>maxU) maxU=this.gIP[imgNum].pixelsUV[i][0];
        		if (this.gIP[imgNum].pixelsUV[i][1]>maxV) maxV=this.gIP[imgNum].pixelsUV[i][1];
        	}
        	if (this.debugLevel>2)System.out.println("calculateSensorMasks, imgNum="+imgNum+", minU="+minU+", maxU="+maxU+", minV="+minV+", maxV="+maxV);
        	// restore the grid rectangle for u,v ->pixel-x,pixel-y
        	double [][][] pXY=new double[maxV-minV+1][maxU-minU+1][2];
        	int [][]iMask=new int [pXY.length][pXY[0].length];
        	for (int v=0;v<pXY.length;v++) for (int u=0;u<pXY[0].length;u++) {
        		pXY[v][u][0]=-1;
        		iMask[v][u]=0;
        	}
        	if (this.debugLevel>2)System.out.println("calculateSensorMasks, pXY.length="+pXY.length+", pXY[0].length="+pXY[0].length);
        	for (i=0;i<this.gIP[imgNum].pixelsXY.length;i++) if ((this.gIP[imgNum].pixelsXY!=null) &&(this.gIP[imgNum].pixelsXY[i][0]>=0)) {
        		int v=this.gIP[imgNum].pixelsUV[i][1]-minV;
        		int u=this.gIP[imgNum].pixelsUV[i][0]-minU;
        		pXY[v][u][0]=this.gIP[imgNum].pixelsXY[i][0]; // out of bounds 22
        		pXY[v][u][1]=this.gIP[imgNum].pixelsXY[i][1];
        		//    			if (this.debugLevel>2)System.out.println("calculateSensorMasks, i="+i+", pXY["+v+"]["+u+"]={"+pXY[v][u][0]+","+pXY[v][u][1]+"}");
        		iMask[v][u]=1;
        	}
        	if (this.debugLevel>3){
        		double [][] testArray=new double[3][pXY.length*pXY[0].length];
        		int index=0;
        		for (int v=0;v<iMask.length;v++) for (int u=0;u<iMask[0].length;u++){
        			testArray[0][index]=pXY[v][u][0];
        			testArray[1][index]=pXY[v][u][1];
        			testArray[2][index++]=iMask[v][u];

        		}
        		String [] dbgTitles={"X","Y","iMask"};
        		this.SDFA_INSTANCE.showArrays(testArray, pXY[0].length, pXY.length,  true, "original", dbgTitles);

        	}
        	// shrink the grid
        	int vMax=iMask.length-1;
        	int uMax=iMask[0].length-1;
        	if (this.debugLevel>2)System.out.println("calculateSensorMasks, uMax="+uMax+", vMax="+vMax);
        	for (int n=0;n<shrinkGridForMask;n++) {
        		for (int v=0;v<iMask.length;v++) for (int u=0;u<iMask[0].length;u++) if (iMask[v][u]>0){

        			if ((v==0) || (v==vMax) || (u==0) || (u==uMax) ||
        					(iMask[v-1][u]==-n) || (iMask[v+1][u]==-n) ||(iMask[v][u-1]==-n) || (iMask[v][u+1]==-n)) {
        				iMask[v][u]=-n-1;
        			}
        		}
        	}
        	if (this.debugLevel>3){
        		double [][] testArray1=new double[3][pXY.length*pXY[0].length];
        		int index=0;
        		for (int v=0;v<iMask.length;v++) for (int u=0;u<iMask[0].length;u++){
        			testArray1[0][index]=pXY[v][u][0];
        			testArray1[1][index]=pXY[v][u][1];
        			testArray1[2][index++]=iMask[v][u];

        		}
        		String [] dbgTitles={"X","Y","iMask"};
        		this.SDFA_INSTANCE.showArrays(testArray1, pXY[0].length, pXY.length,  true, "shrank", dbgTitles);
        	}

        	// now in remaining grid nodes iMask[v][u]>0 (0 and negative - no grid)
        	// accumulate pixels around the grid points
        	for (int v=0;v<iMask.length;v++) for (int u=0;u<iMask[0].length;u++)if (iMask[v][u]>0){
        		// find the radius - distance to the fartherst of the 4 (existent) neighbors (if none exist - disregard the node)
        		double r2Max=0;
        		for (int d=0;d<dirs.length;d++){
        			int u1=u+dirs[d][0];
        			int v1=v+dirs[d][1];
        			double r2;
        			if ((v1>=0) && (v1<=vMax) && (u1>=0) && (u1<=uMax)){
        				r2=(pXY[v1][u][0]-pXY[v][u][0])*(pXY[v1][u][0]-pXY[v][u][0])+
        				(pXY[v][u1][0]-pXY[v][u][0])*(pXY[v][u1][0]-pXY[v][u][0]);
        				if (r2Max<r2) r2Max=r2;
        			}
        		}
        		if (r2Max==0.0) continue; // nothing around - skip this node
        		// calculate average radius (for bluring)
        		double r=Math.sqrt(r2Max);
        		rAverage+=r;
        		rAverageNum++;

        		int iR= (int) Math.round(r);
        		int iX0= (int) Math.round (pXY[v][u][0]);
        		int iY0= (int) Math.round (pXY[v][u][1]);
        		int xLowLim=iX0-iR;
        		int xHighLim=iX0+iR;
        		int yLowLim=iY0-iR;
        		int yHighLim=iY0+iR;
        		if (xLowLim<0)       xLowLim=0;
// decimation apply below
        		if (xHighLim>=width) xHighLim=width-1;
        		if (yLowLim<0)       yLowLim=0;
        		if (yHighLim>=height)yHighLim=height-1;
        		for (int iY=yLowLim;iY<=yHighLim;iY+=decimate)for (int iX=xLowLim;iX<=xHighLim;iX+=decimate){
        			double r2=(iX-pXY[v][u][0])*(iX-pXY[v][u][0])+(iY-pXY[v][u][1])*(iY-pXY[v][u][1]);
        			if (r2<=r2Max) {
        				if (decimate==1) mask[iY*width+iX]=1.0;
        				else mask[(iY/decimate)*dWidth+(iX/decimate)]=1.0;
        			}
        		}
        	}
        	double [][] result= {{rAverage,rAverageNum},mask};
        	return  result;

        }
// Methods related to Talon (instructor/student) systems
        // calculates UV offset (no rotations) for test grid by best contrast match to  known grid
//        int [] correlateGrids(int base_grid_num, int test_grid_num, boolean invert_color) {
// returns
        int [] correlateGrids(
        		int        base_width,
        		float [][] base_pixels,
        		int        test_width,
        		float [][] test_pixels,
        		boolean    invert_color,
        		int        extra_search) {
        	int base_height = base_pixels[0].length/base_width;
        	int test_height = test_pixels[0].length/test_width;
        	int search_rad = Math.max(
        			(Math.max(base_width,  test_width)- Math.min(base_width,  test_width) + 1) /2,
        			(Math.max(base_height, test_height)-Math.min(base_height, test_height) + 1) /2) + extra_search;
        	int offs_x = base_width/2 -  test_width/2;  // subtract from test.x
        	int offs_y = base_height/2 - test_height/2; // subtract from test.y
        	double [] corr = new double [(2*search_rad + 1)*(2*search_rad + 1)];
        	for (int dy = -search_rad; dy <= search_rad; dy++) {
        		for (int dx = -search_rad; dx <= search_rad; dx++) {
        			double sum = 0;
        			for (int y0 = 0; y0 < base_height; y0++) {
        				int y1 = y0 - offs_y - dy;
        				if ((y1 >= 0) && (y1 < test_height)) {

        					for (int x0 = 0; x0 < base_width; x0++) {
        						int x1 = x0 - offs_x - dx;
        						if ((x1 >= 0) && (x1 < test_width)) {
        							sum+= base_pixels[INDEX_CONTRAST][y0*base_width + x0] * test_pixels[INDEX_CONTRAST][y1*test_width + x1];
        						}
        					}
        				}
        			}
        			corr[(2*search_rad + 1)*(search_rad + dy) +(search_rad + dx)] = sum;
        		}
        	}
        	int [] indx_max_even_odd = {0,0};
        	for (int i = 1; i < corr.length; i++) {
        		int parity = ((i /(2*search_rad + 1)) +  (i %(2*search_rad + 1))) & 1;
        		if (corr[indx_max_even_odd[parity]] <corr[i]) {
        			indx_max_even_odd[parity] = i;
        		}
        	}
        	// find first non-zero matching cell
        	int indx0=-1,indx1=-1;
        	int [] rslt = new int[3];
        	for (int parity = 0; parity < 2; parity++) {
        		int dy = indx_max_even_odd[parity] / (2*search_rad + 1) - search_rad;
        		int dx = indx_max_even_odd[parity] % (2*search_rad + 1) - search_rad;

        		first_nonzero:

        			for (int y0 = 0; y0 < base_height; y0++) {
        				int y1 = y0 - offs_y - dy;
        				if ((y1 >= 0) && (y1 < test_height)) {
        					for (int x0 = 0; x0 < base_width; x0++) {
        						int x1 = x0 - offs_x - dx;
        						if ((x1 >= 0) && (x1 < test_width)) {
        							indx0 = y0*base_width + x0;
        							indx1 = y1*test_width + x1;
        							if ((base_pixels[INDEX_CONTRAST][indx0] > 0 )&&
        									(test_pixels[INDEX_CONTRAST][indx1] > 0)) {
        								break first_nonzero;
        							}
        						}
        					}
        				}
        			}
        		// test grid with index indx1 matches base grid with indx0
        		rslt[0] = Math.round(base_pixels[INDEX_U][indx0] - test_pixels[INDEX_U][indx1]);
        		rslt[1] = Math.round(base_pixels[INDEX_V][indx0] - test_pixels[INDEX_V][indx1]);
        		rslt[2] = 0; // rotation
        		if (((rslt[0] + rslt[1] + (invert_color?1:0)) & 1) == 0) break;
        	}
        	return rslt;
        }

        int [] correlateWithPattern(
        		PatternParameters patternParameters,
        		int        test_width,
        		float [][] test_pixels,
        		boolean    invert_color,
        		int        extra_search,
        		double     sigma,
        		double []  sensor_wh, // test set pixels width/height pair to reduce weight near the margins (or null)
        		boolean    bdebug
        		) {
        	int base_height = patternParameters.gridGeometry.length;
        	int base_width =  patternParameters.gridGeometry[0].length;
        	int index_mask =  3;
        	float [][] base_pixels = new float [INDEX_CONTRAST+1][base_height*base_width];
    		int indx = 0;
        	for (int iv = 0; iv < base_height; iv++) {
            	for (int iu = 0; iu < base_width; iu++) {
            		base_pixels[INDEX_U][indx] =        iu - patternParameters.U0;
            		base_pixels[INDEX_V][indx] =        iv - patternParameters.V0;
            		base_pixels[INDEX_CONTRAST][indx] = (float) patternParameters.gridGeometry[iv][iu][index_mask];
            		indx++;
            	}
        	}
        	return correlateGrids(
            		base_width,
            		base_pixels,
            		test_width,
            		test_pixels,
            		invert_color,
            		extra_search,
            		sigma,
            		sensor_wh, // test set pixels width/height pair to reduce weight near the margins (or null)
            		bdebug
            		) ;
        }

        int [] correlateGrids(
        		int        base_width,
        		float [][] base_pixels,
        		int        test_width,
        		float [][] test_pixels,
        		boolean    invert_color,
        		int        extra_search,
        		double     sigma,
        		double []  sensor_wh, // test set pixels width/height pair to reduce weight near the margins (or null)
        		boolean    bdebug
        		) { // Gaussian blur sigma to subtract
        	int base_height = base_pixels[0].length/base_width;
        	int test_height = test_pixels[0].length/test_width;
        	int search_rad = Math.max(
        			(Math.max(base_width,  test_width)- Math.min(base_width,  test_width) + 1) /2,
        			(Math.max(base_height, test_height)-Math.min(base_height, test_height) + 1) /2) + extra_search;
        	int offs_x = base_width/2 -  test_width/2;  // subtract from test.x
        	int offs_y = base_height/2 - test_height/2; // subtract from test.y
        	double [] corr = new double [(2*search_rad + 1)*(2*search_rad + 1)];
        	double [] base_contrast = new double [base_pixels[INDEX_U].length];
        	double [] test_contrast = new double [test_pixels[INDEX_U].length];
        	for (int i = 0; i < base_contrast.length; i++) {
        		base_contrast[i] = base_pixels[INDEX_CONTRAST][i];
        	}
        	for (int i = 0; i < test_contrast.length; i++) {
        		test_contrast[i] = test_pixels[INDEX_CONTRAST][i];
        	}
        	double [] tmp_b = base_contrast.clone();
        	DoubleGaussianBlur gb = new DoubleGaussianBlur();
        	gb.blurDouble(tmp_b, base_width, base_height, sigma, sigma, 0.01);
        	// subtract DC
        	double s = 0.0;
        	for (int i = 0; i < base_contrast.length; i++) {
        		base_contrast[i] -= tmp_b[i];
        		s+=base_contrast[i];
        	}
        	s/=base_contrast.length;
        	for (int i = 0; i < base_contrast.length; i++) {
        		base_contrast[i] -= s;
        	}

        	tmp_b = test_contrast.clone();
        	gb.blurDouble(tmp_b, test_width, test_height, sigma, sigma, 0.01);
        	// subtract DC
        	s = 0.0;
        	for (int i = 0; i < test_contrast.length; i++) {
        		test_contrast[i] -= tmp_b[i];
        		s+=test_contrast[i];
        	}
        	s/=test_contrast.length;
        	for (int i = 0; i < test_contrast.length; i++) {
        		test_contrast[i] -= s;
        	}

        	if (sensor_wh != null) {
        		double x0 = sensor_wh[0]/2;
        		double y0 = sensor_wh[1]/2;
        		for (int i = 0; i < test_contrast.length; i++) {
        			double x = (test_pixels[INDEX_PX][i] - x0)/x0;
        			double y = (test_pixels[INDEX_PY][i] - y0)/y0;
        			test_contrast[i] *= (1.0 - x*x)* (1.0 - y*y);
        		}
        	}
        	if (bdebug) {
        		(new ShowDoubleFloatArrays()).showArrays(base_contrast, base_width, base_height, "base_sigma-"+sigma);
        		(new ShowDoubleFloatArrays()).showArrays(test_contrast, test_width, test_height, "test_sigma-"+sigma);
        	}

        	for (int dy = -search_rad; dy <= search_rad; dy++) {
            	for (int dx = -search_rad; dx <= search_rad; dx++) {
            		double sum = 0;
            		for (int y0 = 0; y0 < base_height; y0++) {
            			int y1 = y0 - offs_y - dy;
            			if ((y1 >= 0) && (y1 < test_height)) {

                    		for (int x0 = 0; x0 < base_width; x0++) {
                    			int x1 = x0 - offs_x - dx;
                    			if ((x1 >= 0) && (x1 < test_width)) {
                    				sum+= base_contrast[y0*base_width + x0] * test_contrast[y1*test_width + x1];
                    			}
                    		}
            			}
            		}
            		corr[(2*search_rad + 1)*(search_rad + dy) +(search_rad + dx)] = sum;
            	}
        	}
           	if (bdebug) {
        		(new ShowDoubleFloatArrays()).showArrays(corr, "corr_sigma-"+sigma);
           	}
        	int [] indx_max_even_odd = {-1,-1};
        	for (int i = 1; i < corr.length; i++) {
        		int parity = ((i /(2*search_rad + 1)) +  (i %(2*search_rad + 1))) & 1;
        		if ((indx_max_even_odd[parity] < 0) || (corr[indx_max_even_odd[parity]] < corr[i])) {
        			indx_max_even_odd[parity] = i;
        		}
        	}
    		// find first non-zero matching cell
    		int indx0=-1,indx1=-1;
    		int [] rslt = new int[3];
        	for (int parity = 0; parity < 2; parity++) {
        		int dy = indx_max_even_odd[parity] / (2*search_rad + 1) - search_rad;
        		int dx = indx_max_even_odd[parity] % (2*search_rad + 1) - search_rad;

        		first_nonzero:

        			for (int y0 = 0; y0 < base_height; y0++) {
        				int y1 = y0 - offs_y - dy;
        				if ((y1 >= 0) && (y1 < test_height)) {
        					for (int x0 = 0; x0 < base_width; x0++) {
        						int x1 = x0 - offs_x - dx;
        						if ((x1 >= 0) && (x1 < test_width)) {
        							indx0 = y0*base_width + x0;
        							indx1 = y1*test_width + x1;
        							if ((base_pixels[INDEX_CONTRAST][indx0] > 0 )&&
        									(test_pixels[INDEX_CONTRAST][indx1] > 0)) {
        								break first_nonzero;
        							}
        						}
        					}
        				}
        			}
            	// test grid with index indx1 matches base grid with indx0
            	rslt[0] = Math.round(base_pixels[INDEX_U][indx0] - test_pixels[INDEX_U][indx1]);
            	rslt[1] = Math.round(base_pixels[INDEX_V][indx0] - test_pixels[INDEX_V][indx1]);
            	rslt[2] = 0; // rotation
            	if (((rslt[0] + rslt[1] + (invert_color?1:0)) & 1) == 0) break;
        	}
        	return rslt;
        }


        // Set initial shifts for the grids that do not have absolute match from the one in the same set that does

        // end of class DistortionCalibrationData


    // Set initial shifts for the grids that do not have absolute match from the one in the same set that does

    // end of class DistortionCalibrationData
}

