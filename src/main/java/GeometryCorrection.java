import ij.IJ;

/**
 **
 ** GeometryCorrection - geometry correction for multiple sensors sharing the same
 ** lens radial distortion model
 **
 ** Copyright (C) 2016 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  GeometryCorrection.java is free software: you can redistribute it and/or modify
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

public class GeometryCorrection {
//	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
//	public double radius;  // mm, distance from the rotation axis
//	public double height;       // mm, up - from the origin point
//	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise heading
//	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up elevation
//	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target roll
	public double focalLength=4.5;
	public double pixelSize=  2.2; //um
	public double distortionRadius=  2.8512; // mm - half width of the sensor
	public double distortionA8=0.0; //r^8 (normalized to focal length or to sensor half width?)
	public double distortionA7=0.0; //r^7 (normalized to focal length or to sensor half width?)
	public double distortionA6=0.0; //r^6 (normalized to focal length or to sensor half width?)
	public double distortionA5=0.0; //r^5 (normalized to focal length or to sensor half width?)
	public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	public double distortionB=0.0; // r^3
	public double distortionC=0.0; // r^2
	public double px0=1296.0;          // center of the lens on the sensor, pixels
	public double py0=968.0;           // center of the lens on the sensor, pixels

	private double [] rByRDist=null;
    private double    stepR=0.001;
    private double    maxR=2.0; // calculate up to this*distortionRadius
	

	
	
public void setDistortion()	{
//	imp.setProperty("distortion_formula",  "(normalized by distortionRadius in mm) Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
//	imp.setProperty("distortionRadius", ""+subCam.distortionRadius);
	
}

/*
	    			if (this.rByRDist==null){
	    				calcReverseDistortionTable();
	    			}
	    			double rND2R=getRByRDist(rD/this.distortionRadius,debugThis);
	    			x*= rND2R; // positive - right
	    			y*=-rND2R; // positive - up

 */


// Copied from PixelMapping
/**
 * Calculate reverse distortion table - from pixel radius to non-distorted radius	
 * Rdist/R=A5*R^4+A*R^3+B*R^2+C*R+(1-A5-A-B-C)    
 * @return false if distortion is too high
 */
public boolean calcReverseDistortionTable(){
	boolean debugThis=false; //true;
	double delta=1E-8;
	double minDerivative=0.1;
	int numIterations=1000;
	double drDistDr=1.0;
//	public double distortionA5=0.0; //r^5 (normalized to focal length or to sensor half width?)
//	public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
//	public double distortionB=0.0; // r^3
//	public double distortionC=0.0; // r^2
	boolean use8=(this.distortionA8!=0.0) || (this.distortionA7!=0.0) || (this.distortionA6!=0.0);
	double d=1.0-this.distortionA8-this.distortionA7-this.distortionA6-this.distortionA5-this.distortionA-this.distortionB-this.distortionC;
	double rPrev=0.0;
	this.rByRDist=new double [(int) Math.ceil(this.maxR/this.stepR)+1];
	for (int j=1;j<this.rByRDist.length;j++) this.rByRDist[j]=Double.NaN;
	this.rByRDist[0]=1.0/d;
	boolean bailOut=false;
	if (debugThis)	System.out.println("calcReverseDistortionTable()");

	for (int i=1;i<this.rByRDist.length;i++) {
		double rDist=this.stepR*i;
		double r=rPrev+this.stepR/drDistDr;
//		if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" rDist="+rDist+" r="+r+" rPrev="+rPrev);

		for (int iteration=0;iteration<numIterations;iteration++){
			double k;
			if (use8){
				k=(((((((this.distortionA8)*r+this.distortionA7)*r+this.distortionA6)*r+this.distortionA5)*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
    			drDistDr=(((((((8*this.distortionA8)*r + 7*this.distortionA7)*r + 6*this.distortionA6)*r + 5*this.distortionA5)*r + 4*this.distortionA)*r+3*this.distortionB)*r+2*this.distortionC)*r+d;
			} else {
				k=(((this.distortionA5*r + this.distortionA)*r+this.distortionB)*r+this.distortionC)*r+d;
    			drDistDr=(((5*this.distortionA5*r + 4*this.distortionA)*r+3*this.distortionB)*r+2*this.distortionC)*r+d;
			}
			double rD=r*k;
			if (drDistDr<minDerivative) {
				bailOut=true;
				break; // too high distortion
			}
			if (Math.abs(rD-rDist)<delta) break; // success
			r+=(rDist-rD)/drDistDr;
		}
		if (bailOut) {
			if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" Bailing out, drDistDr="+drDistDr);
			return false;
		} 
		rPrev=r;
		this.rByRDist[i]=r/rDist;
		if (debugThis)	System.out.println("calcReverseDistortionTable() i="+i+" rDist="+rDist+" r="+r+" rPrev="+rPrev+" this.rByRDist[i]="+this.rByRDist[i]);
	}
	return true;
}

public double getRByRDist(double rDist, boolean debug){
	// add exceptions;
	if (this.rByRDist==null) {
		if (debug)System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): this.rByRDist==null");
		return Double.NaN;
	}
	if (rDist<0) {
		if (debug)System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): rDist<0");
		return Double.NaN;
	}
	int index=(int) Math.floor(rDist/this.stepR);
	if (index>=(this.rByRDist.length-1)) {
		if (debug) System.out.println("getRByRDist("+IJ.d2s(rDist,3)+"): index="+index+">="+(this.rByRDist.length-1));
		return Double.NaN;
	}
	double result=this.rByRDist[index]+(this.rByRDist[index+1]-this.rByRDist[index])*(rDist/this.stepR-index);
	if (Double.isNaN(result)){
		if (debug) System.out.println("this.rByRDist["+index+"]="+this.rByRDist[index]);
		if (debug) System.out.println("this.rByRDist["+(index+1)+"]="+this.rByRDist[index+1]);
		if (debug) System.out.println("rDist="+rDist);
		if (debug) System.out.println("(rDist/this.stepR="+(rDist/this.stepR));
		
	}
	return result;
	
}



}
