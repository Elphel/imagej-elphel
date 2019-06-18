package com.elphel.imagej.calibration;
import java.util.Arrays;
import java.util.Properties;

import com.elphel.imagej.cameras.EyesisSubCameraParameters;
import com.elphel.imagej.common.WindowTools;

import Jama.Matrix;
import ij.IJ;
import ij.gui.GenericDialog;

//import EyesisCameraParameters;
//import Distortions.EyesisSubCameraParameters;
//import LensDistortionParameters;
//import PatternParameters;
//import DistortionCalibrationData.GridImageParameters;
//import DistortionCalibrationData.GridImageSet;

	public class LensDistortionParameters {
		/*
		 * Hugin Rsrc=a*Rdest^4+b*Rdest^3+c*Rdest^2+d*Rdest; d=1-(a+b+c)
		 */
		// lens parameters (add more later)

		final int numInputs= 53; //27; // with A8...// 24;   // parameters in subcamera+...
	    final int numOutputs=42; //16; // with A8...//13;  // parameters in a single camera
//		final
//		boolean cummulativeCorrection=true; // r_xy, r_od for higher terms are relative to lower ones
	    public int defaultLensDistortionModel=200;
	    public int lensDistortionModel=defaultLensDistortionModel;
	    public int lensDistortionModelType=0; // set from lensDistortionModel
		boolean cummulativeCorrection=false; //true; // r_xy, r_od for higher terms are relative to lower ones
		public double focalLength=4.5;
		public double pixelSize=  2.2; //um
		public double distortionRadius=  2.8512; // mm - half width of the sensor
		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
		public double distortionB=0.0; // r^3
		public double distortionC=0.0; // r^2
		// orientation/position parameters
		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
		public double x0=0;      // lens axis from pattern center, mm (to the right)
		public double y0=0;      // lens axis from pattern center, mm (down)
		public double z0=0;      // lens axis from pattern center, mm (away from the camera)
		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
		public double px0=1296.0;           // center of the lens on the sensor, pixels
		public double py0=968.0;           // center of the lens on the sensor, pixels
		public boolean flipVertical; // acquired image is mirrored vertically (mirror used)

		// new non-radial parameters
		final private double [][] r_xy_dflt={{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; // only 6, as for the first term delta x, delta y ==0
		final private double [][] r_od_dflt=   {{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; // ortho

		public double [][] r_xy=null; // only 6, as for the first term delta x, delta y ==0
		public double [][] r_od=null; // ortho

		public double [][] r_xyod=null; //{x0,y0,ortho, diagonal}

		// total number of new parameters = 6*2+7+7=26

		public int debugLevel=1; // was 2
/*
 Modifying to accommodate for eccentricity of different terms (2 parameters per term) and elliptical shape (another 2 terms). When all are
 zeroes, the old parameters are in effect:
		Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
		Rdist/R=A8*R6^7+A7*R5^6+A6*R4^5+A5*R3^4+A*R2^3+B*R1^2+C*R0+(1-A6-A7-A6-A5-A-B-C)");
		R[i] depends on px,py and r_xy[i][], r_o[i] (positive - "landscape", negative - "portrait"), r_d (positive - along y=x, negative - along y=-x)

		R[i] = r0[i]*(1+  (r_od[i][0]*(y[i]**2-x[i]**2)+ 2*r_od[i][1]*x[i]*y[i])/r0[i]**2;
		r0[i]=sqrt(x[i]**2+y[2]**2)
		x[i]=pixel_x-px0-((i>0)?r_xy[i-1][0]:0)
		y[i]=pixel_y-py0-((i>0)?r_xy[i-1][1]:0)
*/

// intermediate values
		public double phi, theta,psi,cPH,sPH,cTH,sTH,cPS,sPS;
		public double [][] rotMatrix=new double[3][3]; // includes mirroring for Y (target coordinates y- down, camera - y  up)
		public void setDistortionModelParameters(){
			if (lensDistortionModel<0){
				System.out.println("BUG:setDistortionModelParameters() - lensDistortionModel<0");
				lensDistortionModel=defaultLensDistortionModel;
			}
			this.cummulativeCorrection= (lensDistortionModel==101);
			if (lensDistortionModel<100){
				lensDistortionModelType=0;
			} else if (lensDistortionModel<200){
				lensDistortionModelType=1;
			} else if (lensDistortionModel<300){
				lensDistortionModelType=2;
			} else {
				lensDistortionModelType=0;
			}
		}
	    public LensDistortionParameters(
	    		boolean isTripod,
	    		boolean cartesian,
	    		double  pixelSize,
	    		double  distortionRadius,
	            double [][] interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
	    		double [] parVect,
	    		boolean [] mask, // calculate only selected derivatives (all parVect values are still
	    		int debugLevel
	    		){
	    	this.debugLevel=debugLevel;
	    	lensCalcInterParamers( // changed name to move calcInterParamers method from enclosing class
	    			this,
		    		isTripod,
		    		cartesian,
		    		pixelSize,
		    		distortionRadius,
		            interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
		    		parVect,
		    		mask // calculate only selected derivatives (all parVect values are still
		    		);
	    }

		public LensDistortionParameters(
				double focalLength,
				double pixelSize,  //um
				double distortionRadius, // mm
				double distortionA8, // r^8
				double distortionA7, // r^7
				double distortionA6, // r^6
				double distortionA5, // r^5
				double distortionA, // r^4
				double distortionB, // r^3
				double distortionC, // r^2
				// orientation/position parameters
				double yaw,   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
				double pitch, // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
				double roll,  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
				double x0,      // lens axis from pattern center, mm (to the right)
				double y0,      // lens axis from pattern center, mm (down)
				double z0,      // lens axis from pattern center, mm (away from the camera)
				double distance, // distance from the lens input pupil to the pattern plane along the camera axis, mm
				double px0,           // center of the lens on the sensor, pixels
				double py0,           // center of the lens on the sensor, pixels
				boolean flipVertical, // acquired image is mirrored vertically (mirror used)
				int lensDistortionModel,
				double [][] r_xy,
				double [][] r_od
		){
			setLensDistortionParameters(
			focalLength,
			pixelSize,
			distortionRadius,
			distortionA8,
			distortionA7,
			distortionA6,
			distortionA5,
			distortionA,
			distortionB,
			distortionC,
			yaw,
			pitch,
			roll,
			x0,
			y0,
			z0,
			distance,
			px0,
			py0,
			flipVertical,
			lensDistortionModel,
			r_xy,
			r_od);
		}

		public LensDistortionParameters(){
			setLensDistortionParameters(
			4.5,    // focalLength,
			2.2,    // pixelSize,
			2.8512, // distortionRadius,
			0.0,    // distortionA8,
			0.0,    // distortionA7,
			0.0,    // distortionA6,
			0.0,    // distortionA5,
			0.0,    // distortionA,
			0.0,    // distortionB,
			0.0,    // distortionC,
			0.0,    // yaw,
			0.0,    // pitch,
			0.0,    // roll,
			0.0,    // x0,
			0.0,    // y0,
			0.0,    // z0,
			2360.0, // distance,
			1296,   // px0,
			698,    // py0,
			true,   // flipVertical,
			-1,     // lensDistortionModel
			null,   // r_xy,
			null   // r_od,
			);
		}
		@Override
		public LensDistortionParameters clone() {
			return new LensDistortionParameters(
					this.focalLength,
					this.pixelSize,
					this.distortionRadius,
					this.distortionA8,
					this.distortionA7,
					this.distortionA6,
					this.distortionA5,
					this.distortionA,
					this.distortionB,
					this.distortionC,
					this.yaw,
					this.pitch,
					this.roll,
					this.x0,
					this.y0,
					this.z0,
					this.distance,
					this.px0,
					this.py0,
					this.flipVertical,
					this.lensDistortionModel,
					this.r_xy,
					this.r_od
			);
		}
		public void setLensDistortionParameters(
				double focalLength,
				double pixelSize,  //um
				double distortionRadius, // mm
				double distortionA8, // r^7
				double distortionA7, // r^6
				double distortionA6, // r^5
				double distortionA5, // r^4
				double distortionA, // r^4
				double distortionB, // r^3
				double distortionC, // r^2
				// orientation/position parameters
				double yaw,   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
				double pitch, // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
				double roll,  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
				double x0,      // lens axis from pattern center, mm (to the right)
				double y0,      // lens axis from pattern center, mm (down)
				double z0,      // lens axis from pattern center, mm (away)
				double distance, // distance from the lens input pupil to the pattern plane along the camera axis, mm
				double px0,           // center of the lens on the sensor, pixels
				double py0,           // center of the lens on the sensor, pixels
				boolean flipVertical, // acquired image is mirrored vertically (mirror used)
				int lensDistortionModel,
				double [][] r_xy,   // per polynomial term center x,y correction only 6, as for the first term delta x, delta y ==0
				double [][] r_od   // per polynomial term orthogonal+diagonal elongation
		){
			this.focalLength=focalLength;
			this.pixelSize=pixelSize;
			this.distortionRadius=distortionRadius;
			this.distortionA8=distortionA8;
			this.distortionA7=distortionA7;
			this.distortionA6=distortionA6;
			this.distortionA5=distortionA5;
			this.distortionA=distortionA;
			this.distortionB=distortionB;
			this.distortionC=distortionC;
			this.yaw=yaw;
			this.pitch=pitch;
			this.roll=roll;
			this.x0=x0;
			this.y0=y0;
			this.z0=z0;
			this.distance=distance;
			this.px0=px0;
			this.py0=py0;
			this.flipVertical=flipVertical;
			this.lensDistortionModel=(lensDistortionModel>=0)?lensDistortionModel:defaultLensDistortionModel;
			if (r_xy==null) r_xy=r_xy_dflt;
			if (r_od==null) r_od=r_od_dflt;
			this.r_xy=new double [r_xy.length][2];
			for (int i=0;i<r_xy.length;i++)this.r_xy[i]=r_xy[i].clone();
			this.r_od=new double [r_od.length][2];
			for (int i=0;i<r_od.length;i++)this.r_od[i]=r_od[i].clone();
			recalcCommons();
		}

		public void setIntrincicFromSubcamera(EyesisSubCameraParameters pars){
			setLensDistortionParameters(
					pars.focalLength,
					pars.pixelSize,  //um
					pars.distortionRadius, // mm
					pars.distortionA8, // r^7
					pars.distortionA7, // r^6
					pars.distortionA6, // r^5
					pars.distortionA5, // r^4
					pars.distortionA,  // r^4
					pars.distortionB,  // r^3
					pars.distortionC,  // r^2
					// orientation/position parameters
					this.yaw,          // (keep) angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
					this.pitch,        // (keep) angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
					this.roll,         // (keep) angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
					this.x0,           // (keep) lens axis from pattern center, mm (to the right)
					this.y0,           // (keep) lens axis from pattern center, mm (down)
					this.z0,           // (keep) lens axis from pattern center, mm (away)
					this.distance,     // (keep) distance from the lens input pupil to the pattern plane along the camera axis, mm
					pars.px0,          //        center of the lens on the sensor, pixels
					pars.py0,          // center of the lens on the sensor, pixels
					this.flipVertical, // (keep)  acquired image is mirrored vertically (mirror used)
					pars.lensDistortionModel,
					pars.r_xy,         // do not exist yet!
					pars.r_od          // do not exist yet!
			);
		}

		public void setLensDistortionParameters(LensDistortionParameters ldp
		){
			setLensDistortionParameters(
					ldp.focalLength,
					ldp.pixelSize,
					ldp.distortionRadius,
					ldp.distortionA8,
					ldp.distortionA7,
					ldp.distortionA6,
					ldp.distortionA5,
					ldp.distortionA,
					ldp.distortionB,
					ldp.distortionC,
					ldp.yaw,
					ldp.pitch,
					ldp.roll,
					ldp.x0,
					ldp.y0,
					ldp.z0,
					ldp.distance,
					ldp.px0,
					ldp.py0,
					ldp.flipVertical,
					ldp.lensDistortionModel,
					ldp.r_xy,
					ldp.r_od);
		}




// TODO: Fix for non-radial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// just for debugging
		public void setLensDistortionParameters(LensDistortionParameters ldp,
				int index, // parameter to add delta, 1..13->14->17
				double delta
		){
/*
			this.focalLength=ldp.focalLength+((index==7)?delta:0);
			this.pixelSize=ldp.pixelSize;
			this.distortionRadius=ldp.distortionRadius;
			this.distortionA8=ldp.distortionA8+((index==9)?delta:0);
			this.distortionA7=ldp.distortionA7+((index==10)?delta:0);
			this.distortionA6=ldp.distortionA6+((index==11)?delta:0);
			this.distortionA5=ldp.distortionA5+((index==12)?delta:0);
			this.distortionA=ldp.distortionA+((index==13)?delta:0);
			this.distortionB=ldp.distortionB+((index==14)?delta:0);
			this.distortionC=ldp.distortionC+((index==15)?delta:0);
			this.yaw=ldp.yaw+((index==1)?delta:0);
			this.pitch=ldp.pitch+((index==2)?delta:0);
			this.roll=ldp.roll+((index==3)?delta:0);
			this.x0=ldp.x0+((index==4)?delta:0);
			this.y0=ldp.y0+((index==5)?delta:0);
			this.z0=ldp.z0+((index==6)?delta:0);
			this.distance=ldp.distance+((index==8)?delta:0);
			this.px0=ldp.px0+((index==16)?delta:0);
			this.py0=ldp.py0+((index==17)?delta:0);
			this.flipVertical=ldp.flipVertical;
*/
			setLensDistortionParameters(ldp);
			final int index_r_xy=18;
			final int index_r_od=30;
			final int index_end=44;

			switch (index){
			case 1: this.yaw+=delta; break;
			case 2: this.pitch+=delta; break; //=ldp.pitch+((index==2)?delta:0);
			case 3: this.roll+=delta; break; // =ldp.roll+((index==3)?delta:0);
			case 4: this.x0+=delta; break; // =ldp.x0+((index==4)?delta:0);
			case 5: this.y0+=delta; break; // =ldp.y0+((index==5)?delta:0);
			case 6: this.z0+=delta; break; // =ldp.z0+((index==6)?delta:0);
			case 7: this.focalLength+=delta; break; // =ldp.focalLength+((index==7)?delta:0);
			case 8: this.distance+=delta; break; // =ldp.distance+((index==8)?delta:0);
			case 9: this.distortionA8+=delta; break; // =ldp.distortionA8+((index==9)?delta:0);
			case 10: this.distortionA7+=delta; break; // =ldp.distortionA7+((index==10)?delta:0);
			case 11: this.distortionA6+=delta; break; // =ldp.distortionA6+((index==11)?delta:0);
			case 12: this.distortionA5+=delta; break; // =ldp.distortionA5+((index==12)?delta:0);
			case 13: this.distortionA+=delta; break; // =ldp.distortionA+((index==13)?delta:0);
			case 14: this.distortionB+=delta; break; // =ldp.distortionB+((index==14)?delta:0);
			case 15: this.distortionC+=delta; break; // =ldp.distortionC+((index==15)?delta:0);
			case 16: this.px0+=delta; break; // =ldp.px0+((index==16)?delta:0);
			case 17: this.py0+=delta; break; // =ldp.py0+((index==17)?delta:0);
			default:
				if ((index>=index_r_xy) && (index<index_r_od)){
					this.r_xy[(index-index_r_xy)/2][(index-index_r_xy)%2]+=delta;
				} else if ((index>=index_r_od) && (index<index_end)){
					this.r_od[(index-index_r_od)/2][(index-index_r_od)%2]+=delta;
				}
			}
			recalcCommons();
		}

// recalculate common (point-invariant) intermediate values (cos, sin, rotation matrix)

		public void recalcCommons(){
//        	this.cummulativeCorrection=false; // just debugging
//		public double phi, theta,psi,cPH,sPH,cTH,sTH,cPS,sPS;
//    		public double [][] rotMatrix=new double[3][3];
        	this.phi=  this.yaw*Math.PI/180;
        	this.theta=this.pitch*Math.PI/180;
        	this.psi=  this.roll*Math.PI/180;
        	this.sPH=Math.sin(this.phi);
        	this.cPH=Math.cos(this.phi);
        	this.sTH=Math.sin(this.theta);
        	this.cTH=Math.cos(this.theta);
        	this.sPS=Math.sin(this.psi);
        	this.cPS=Math.cos(this.psi);
/*
| Xe |   |   0  |   | cPS*cPH+sPS*sTH*sPH   -sPS*cTH  -cPS*sPH+sPS*sTH*cPH |   | Xp |
| Ye | = |   0  | + | sPS*cPH-cPS*sTH*sPH    cPS*cTH  -sPS*sPH-cPS*sTH*cPH | * |-Yp |
| Ze |   | dist |   | cTH*sPH                    sTH       cTH*cPH         |   | Zp |

| PX | =(1000*f)/(Ze*Psz) * |  Xe | + | PX0 |
| PY | =                    | -Ye |   | PY0 |


 Xe =   (cPS*cPH+sPS*sTH*sPH)*Xp   +sPS*cTH*Yp  +(-cPS*sPH+sPS*sTH*cPH)*Zp
 Ye =   (sPS*cPH-cPS*sTH*sPH)*Xp   -cPS*cTH*Yp  +(-sPS*sPH-cPS*sTH*cPH)*Zp
 Ze =   (cTH*sPH)*Xp               -sTH*Yp      +( cTH*cPH)* Zp + dist


theta==0, psi==0:
 Xe =   (cPH)*Xp
 Ye =   Yp
 Ze =   cPH* Zp + dist

(4) PXmmc =f/(cPH* Zp + dist)*  (cPH)*Xp  // mm, left from the lens axis intersection with the sensor

dPXmmc/dphi=

 */
        	this.rotMatrix[0][0]= cPS*cPH+sPS*sTH*sPH;
        	this.rotMatrix[0][1]= sPS*cTH;
           	this.rotMatrix[0][2]=-cPS*sPH+sPS*sTH*cPH;
           	this.rotMatrix[1][0]= sPS*cPH-cPS*sTH*sPH;
           	this.rotMatrix[1][1]=-cPS*cTH;
           	this.rotMatrix[1][2]=-sPS*sPH-cPS*sTH*cPH;
           	this.rotMatrix[2][0]= cTH*sPH;
           	this.rotMatrix[2][1]=-sTH;
           	this.rotMatrix[2][2]= cTH*cPH;
           	if (this.debugLevel>2){
           		System.out.println("recalcCommons():this.rotMatrix:");
        		(new Matrix(this.rotMatrix)).print(10, 5);
           	}
           	//
    		this.r_xyod=new double [this.r_od.length][4]; //{x0,y0,ortho, diagonal}
    		this.r_xyod[0][0]=0.0; // this.px0; //
    		this.r_xyod[0][1]=0.0; // this.py0;
    		this.r_xyod[0][2]=this.r_od[0][0];
    		this.r_xyod[0][3]=this.r_od[0][1];
    		setDistortionModelParameters();
    		if (cummulativeCorrection){
    			for (int i=1;i<this.r_xyod.length;i++){
    	    		this.r_xyod[i][0]=this.r_xyod[i-1][0]+this.r_xy[i-1][0];
    	    		this.r_xyod[i][1]=this.r_xyod[i-1][1]+this.r_xy[i-1][1];
    	    		this.r_xyod[i][2]=this.r_xyod[i-1][2]+this.r_od[i][0];
    	    		this.r_xyod[i][3]=this.r_xyod[i-1][3]+this.r_od[i][1];
    			}
    		} else {
    			for (int i=1;i<this.r_xyod.length;i++){
    	    		this.r_xyod[i][0]=this.r_xy[i-1][0]; // this.px0+this.r_xy[i-1][0];
    	    		this.r_xyod[i][1]=this.r_xy[i-1][1]; // this.py0+this.r_xy[i-1][1];
    	    		this.r_xyod[i][2]=this.r_od[i][0];
    	    		this.r_xyod[i][3]=this.r_od[i][1];
    			}
    		}
        }

        private String [][] descriptions=
        {
        		{"distance",    "Distance from the intersection of the lens axis with z=0 target plane to the camera lens entrance pupil", "mm", "e"},
        		{"x0",          "Lens axis from pattern center, (to the right)", "mm","e"},
        		{"y0",          "Lens axis from pattern center, (down)", "mm","e"},
        		{"yaw",         "Angle from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top", "degrees","e"},
        		{"pitch",       "Angle from perpendicular to the pattern, 0 - towards wall, positive - up", "degrees","e"},
        		{"roll",        "Angle around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise", "degrees","e"},
        		{"focalLength", "Lens focal length", "mm","i"},
        		{"px0",         "horizontal (left to right) pixel number of the lens axis intersection with the sensor", "pix","i"},
        		{"py0",         "vertical (up to down) pixel number of the lens axis intersection with the sensor", "pix","i"},
        		{"distortionA8","Lens distortion coefficient for r^8 (r0=half sensor width)", "r0","i"},
        		{"distortionA7","Lens distortion coefficient for r^7 (r0=half sensor width)", "r0","i"},
        		{"distortionA6","Lens distortion coefficient for r^6 (r0=half sensor width)", "r0","i"},
        		{"distortionA5","Lens distortion coefficient for r^5 (r0=half sensor width)", "r0","i"},
        		{"distortionA", "Lens distortion coefficient for r^4 (r0=half sensor width)", "r0","i"},
        		{"distortionB", "Lens distortion coefficient for r^3 (r0=half sensor width)", "r0","i"},
        		{"distortionC", "Lens distortion coefficient for r^2 (r0=half sensor width)", "r0","i"},
        		{"elong_C_o",   "Orthogonal elongation for r^2", "rel","i"},
        		{"elong_C_d",   "Diagonal   elongation for r^2", "rel","i"},

        		{"eccen_B_x",   "Distortion center shift X for r^3", "rel","i"},
        		{"eccen_B_y",   "Distortion center shift Y for r^3", "rel","i"},
        		{"elong_B_o",   "Orthogonal elongation for r^3", "rel","i"},
        		{"elong_B_d",   "Diagonal   elongation for r^3", "rel","i"},

        		{"eccen_A_x",   "Distortion center shift X for r^4", "rel","i"},
        		{"eccen_A_y",   "Distortion center shift Y for r^4", "rel","i"},
        		{"elong_A_o",   "Orthogonal elongation for r^4", "rel","i"},
        		{"elong_A_d",   "Diagonal   elongation for r^4", "rel","i"},

        		{"eccen_A5_x",  "Distortion center shift X for r^5", "rel","i"},
        		{"eccen_A5_y",  "Distortion center shift Y for r^5", "rel","i"},
        		{"elong_A5_o",  "Orthogonal elongation for r^5", "rel","i"},
        		{"elong_A5_d",  "Diagonal   elongation for r^5", "rel","i"},

        		{"eccen_A6_x",  "Distortion center shift X for r^6", "rel","i"},
        		{"eccen_A6_y",  "Distortion center shift Y for r^6", "rel","i"},
        		{"elong_A6_o",  "Orthogonal elongation for r^6", "rel","i"},
        		{"elong_A6_d",  "Diagonal   elongation for r^6", "rel","i"},

        		{"eccen_A7_x",  "Distortion center shift X for r^7", "rel","i"},
        		{"eccen_A7_y",  "Distortion center shift Y for r^7", "rel","i"},
        		{"elong_A7_o",  "Orthogonal elongation for r^7", "rel","i"},
        		{"elong_A7_d",  "Diagonal   elongation for r^7", "rel","i"},

        		{"eccen_A8_x",  "Distortion center shift X for r^8", "rel","i"},
        		{"eccen_A8_y",  "Distortion center shift Y for r^8", "rel","i"},
        		{"elong_A8_o",  "Orthogonal elongation for r^8", "rel","i"},
        		{"elong_A8_d",  "Diagonal   elongation for r^8", "rel","i"}
        };
        private int numberIntExtrisic(String type){
        	int num=0;
        	for (int i=0;i<this.descriptions.length;i++) if (type.indexOf(descriptions[i][3].charAt(0))>=0) num++;
        	return num;
        }
        /**
         * Verifies that the camera is looking towards the target
         * @return true if looking to the target, false - if away
         */
        public boolean isTargetVisible(boolean verbose){
        	if (verbose) System.out.println("isTargetVisible(): this.distance="+this.distance+", this.yaw="+this.yaw+", this.pitch="+this.pitch);
        	if (this.distance <=0.0) return false;
        	if (Math.cos(this.yaw*Math.PI/180)<0.0) return false;
        	if (Math.cos(this.pitch*Math.PI/180)<0.0) return false;
        	return true;
        }

        public double [] getExtrinsicVector(){
        	double [] extVector = {
        			this.distance,
        			this.x0,
        			this.y0,
        			this.yaw,
        			this.pitch,
        			this.roll
        			};
        	return extVector;
        }
        public double [] getIntrinsicVector(){
        	double [] intVector = {
        			this.focalLength,
        			this.px0,
        			this.py0,
        			this.distortionA8,
        			this.distortionA7,
        			this.distortionA6,
        			this.distortionA5,
        			this.distortionA,
        			this.distortionB,
        			this.distortionC,
        			this.r_od[0][0],
        			this.r_od[0][1],

        			this.r_xy[0][0],
        			this.r_xy[0][1],
        			this.r_od[1][0],
        			this.r_od[1][1],

        			this.r_xy[1][0],
        			this.r_xy[1][1],
        			this.r_od[2][0],
        			this.r_od[2][1],

        			this.r_xy[2][0],
        			this.r_xy[2][1],
        			this.r_od[3][0],
        			this.r_od[3][1],

        			this.r_xy[3][0],
        			this.r_xy[3][1],
        			this.r_od[4][0],
        			this.r_od[4][1],

        			this.r_xy[4][0],
        			this.r_xy[4][1],
        			this.r_od[5][0],
        			this.r_od[5][1],

        			this.r_xy[5][0],
        			this.r_xy[5][1],
        			this.r_od[6][0],
        			this.r_od[6][1]
        	};
        	return intVector;
        }
        public double [] getAllVector(){
        	double [] extVector=getExtrinsicVector();
        	double [] intVector=getIntrinsicVector();
        	double [] allVector = new double[extVector.length+intVector.length];
        	int index=0;
        	for (int i=0;i<extVector.length;i++) allVector[index++]= extVector[i];
        	for (int i=0;i<intVector.length;i++) allVector[index++]= intVector[i];
        	return allVector;
        }
        public void setAllVector(double [] vector){
        	if (vector.length!=(getExtrinsicVector().length+getIntrinsicVector().length)){
        		String msg="Parameter vector should have exactly"+(getExtrinsicVector().length+getIntrinsicVector().length)+" elements";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	this.distance=    vector[ 0];
        	this.x0=          vector[ 1];
        	this.y0=          vector[ 2];
        	this.yaw=         vector[ 3];
        	this.pitch=       vector[ 4];
        	this.roll=        vector[ 5];
        	this.focalLength= vector[ 6];
        	this.px0=         vector[ 7];
        	this.py0=         vector[ 8];
        	this.distortionA8=vector[ 9];
        	this.distortionA7=vector[10];
        	this.distortionA6=vector[11];
        	this.distortionA5=vector[12];
        	this.distortionA= vector[13];
        	this.distortionB= vector[14];
        	this.distortionC= vector[15];
        	int index=16;
        	for (int i=0;i<this.r_od.length;i++){
        		if (i>0) {
            		this.r_xy[i-1][0]=vector[index++];
            		this.r_xy[i-1][1]=vector[index++];
        		}
        		this.r_od[i][0]=vector[index++];
        		this.r_od[i][1]=vector[index++];
        	}
        	/*
        	for (int i=0;i<this.r_xy.length;i++){
        		this.r_xy[i][0]=vector[index++];
        		this.r_xy[i][1]=vector[index++];
        	}
        	for (int i=0;i<this.r_od.length;i++){
        		this.r_od[i][0]=vector[index++];
        		this.r_od[i][1]=vector[index++];
        	}
        	*/
        	recalcCommons();
        }

        public String [] getExtrinsicNames()       {return getDescriptionStrings("e", 0);}
        public String [] getExtrinsicDescriptions(){return getDescriptionStrings("e", 1);}
        public String [] getExtrinsicUnits()       {return getDescriptionStrings("e", 2);}
        public String [] getIntrinsicNames()       {return getDescriptionStrings("i", 0);}
        public String [] getIntrinsicDescriptions(){return getDescriptionStrings("i", 1);}
        public String [] getIntrinsicUnits()       {return getDescriptionStrings("i", 2);}
        public String [] getAllNames()             {return getDescriptionStrings("ei", 0);}
        public String [] getAllDescriptions()      {return getDescriptionStrings("ei", 1);}
        public String [] getAllUnits()             {return getDescriptionStrings("ei", 2);}
        public String [] getAllFlags()             {return getDescriptionStrings("ei", 3);}

        private String [] getDescriptionStrings(String type, int var){
        	String [] s=new String [numberIntExtrisic(type)];
        	int num=0;
        	for (int i=0;i<this.descriptions.length;i++) if (type.indexOf(descriptions[i][3].charAt(0))>=0) s[num++]=this.descriptions[i][var];
        	return s;
        }


/*
 * Calculate pixel value of projection from the pattern point [xp,yp,zp] using current distortion/position parameters
(1) Xe =   (cPS*cPH+sPS*sTH*sPH)*(Xp-X0)   +sPS*cTH*(Yp-Y0)  +(-cPS*sPH+sPS*sTH*cPH)*(Zp-Z0)
(2) Ye =   (sPS*cPH-cPS*sTH*sPH)*(Xp-X0)   -cPS*cTH*(Yp-Y0)  +(-sPS*sPH-cPS*sTH*cPH)*(Zp-Z0)
(3) Ze =   (cTH*sPH)*(Xp-X0)               -sTH*(Yp-Y0)      +( cTH*cPH)* (Zp-Z0) + dist

(4) PXmmc =f/Ze*  Xe  // mm, left from the lens axis intersection with the sensor
(5) PYmmc =f/Ze*  Ye  // mm, up from the lens axis intersection with the sensor

(6) r=sqrt(PXmmc^2+PYmmc^2) // distance from the image point to the lens axis intersection with the sensor (pinhole model)
(7) kD=(Da*(r/r0)^3+Db*(r/r0)^2+Dc*(r/r0)^1+(1-Da-Db-Dc)) correction to the actual distance from the image point to the lens axis due to distortion
(8) xDist = kD *  PXmmc // horisontal distance (mm) from the lens axis on the sensor to the image point, mm (positive - right)
(9) yDist = kD *  PYmmc // vertical distance (mm) from the lens axis on the sensor to the image point, mm (positive - up)


(10) PX = 1000/Psz*( xDist) + PX0 // horizontal pixel of the image (positive - right)
(11) PY = 1000/Psz*(-yDist) + PY0 // vertical pixel of the image (positive - down)

 */
        public double [] patternToPixels(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp // target point horizontal, positive - away from camera,  mm
        		){
        	return calcPartialDerivatives(xp,yp,zp,false)[0];
        }
/*
0		public double x0=0;      // lens axis from pattern center, mm (to the right)
1		public double y0=0;      // lens axis from pattern center, mm (down)
2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
6		public double focalLength=4.5;
7		public double px0=1296.0;           // center of the lens on the sensor, pixels
8		public double py0=968.0;           // center of the lens on the sensor, pixels
 		public double distortionRadius=  2.8512; // mm - half width of the sensor
9		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
10		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
11		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
14		public double distortionB=0.0; // r^3
15		public double distortionC=0.0; // r^2

 */
        /**
         *  extract needed ones, and reorder partial derivatives to match calcInterParameters
         *  @param srcDerivatives - values and 15 derivatives for px, py
         */
        public double [][] reorderPartialDerivatives (double [][] srcDerivatives){
        	int [] order={
        			 4, // 0		public double x0=0;      // lens axis from pattern center, mm (to the right)
        			 5, // 1		public double y0=0;      // lens axis from pattern center, mm (down)
        			 8, // 2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
        			 1, // 3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
        			 2, // 4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
        			 3, // 5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
        			 7, // 6		public double focalLength=4.5;
        			16, // 7		public double px0=1296.0;           // center of the lens on the sensor, pixels
        			17, // 8		public double py0=968.0;           // center of the lens on the sensor, pixels
        			    // 	public double distortionRadius=  2.8512; // mm - half width of the sensor
       			     9, // 9		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
    			    10, // 10		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
    			    11, // 11		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			12, // 12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			13, // 13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
        			14, // 14		public double distortionB=0.0; // r^3
        			15, // 15		public double distortionC=0.0; // r^2
        			18, // 16		Orthogonal elongation for r^2"
        			19, // 17		Diagonal   elongation for r^2"
        			20, // 18		Distortion center shift X for r^3"
        			21, // 19		Distortion center shift Y for r^3"
        			22, // 20		Orthogonal elongation for r^3"
        			23, // 21		Diagonal   elongation for r^3"
        			24, // 22		Distortion center shift X for r^4"
        			25, // 23		Distortion center shift Y for r^4"
        			26, // 24		Orthogonal elongation for r^4"
        			27, // 25		Diagonal   elongation for r^4"
        			28, // 26		Distortion center shift X for r^5"
        			29, // 27		Distortion center shift Y for r^5"
        			30, // 28		Orthogonal elongation for r^5"
        			31, // 29		Diagonal   elongation for r^5"
        			32, // 30		Distortion center shift X for r^6"
        			33, // 31		Distortion center shift Y for r^6"
        			34, // 32		Orthogonal elongation for r^6"
        			35, // 33		Diagonal   elongation for r^6"
        			36, // 34		Distortion center shift X for r^7"
        			37, // 35		Distortion center shift Y for r^7"
        			38, // 36		Orthogonal elongation for r^6"
        			39, // 37		Diagonal   elongation for r^6"
        			40, // 38		Distortion center shift X for r^8"
        			41, // 29		Distortion center shift Y for r^8"
        			42, // 40		Orthogonal elongation for r^8"
        			43  // 41		Diagonal   elongation for r^8"
        	};
        	double [][] result = new double [order.length][2];
        	for (int i=0; i<order.length;i++){
        		result[i][0]=srcDerivatives[order[i]][0];
        		result[i][1]=srcDerivatives[order[i]][1];
        	}
        	return result;
        }
/*
 * result = {{srcDerivatives[4][0],srcDerivatives[4][1]}, - X0
 *           {srcDerivatives[5][0],srcDerivatives[5][1]}  - Y0
 *           {srcDerivatives[8][0],srcDerivatives[8][1]}  - dist
 *           ...
 */


/**
 * Reorder to match the sequence of names - seems to be different :-(
 * @param srcDerivatives  values and 15 derivatives for px, py
 * @return
 */
        public double [][] reorderPartialDerivativesAsNames (double [][] srcDerivatives){
        	int [] order={
        			 8, // 2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
        			 4, // 0		public double x0=0;      // lens axis from pattern center, mm (to the right)
        			 5, // 1		public double y0=0;      // lens axis from pattern center, mm (down)
        			 1, // 3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
        			 2, // 4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
        			 3, // 5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
        			 7, // 6		public double focalLength=4.5;
        			16, // 7		public double px0=1296.0;           // center of the lens on the sensor, pixels
        			17, // 8		public double py0=968.0;           // center of the lens on the sensor, pixels
        			 9, // 9		public double distortionA8=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			10, // 10		public double distortionA7=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			11, // 11		public double distortionA6=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			12, // 12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
        			13, // 13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
        			14, // 14		public double distortionB=0.0; // r^3
        			15, // 15		public double distortionC=0.0; // r^2
        			18, // 16		Orthogonal elongation for r^2"
        			19, // 17		Diagonal   elongation for r^2"

        			20, // 18		Distortion center shift X for r^3"
        			21, // 19		Distortion center shift Y for r^3"
        			22, // 20		Orthogonal elongation for r^3"
        			23, // 21		Diagonal   elongation for r^3"

        			24, // 22		Distortion center shift X for r^4"
        			25, // 23		Distortion center shift Y for r^4"
        			26, // 24		Orthogonal elongation for r^4"
        			27, // 25		Diagonal   elongation for r^4"

        			28, // 26		Distortion center shift X for r^5"
        			29, // 27		Distortion center shift Y for r^5"
        			30, // 28		Orthogonal elongation for r^5"
        			31, // 29		Diagonal   elongation for r^5"

        			32, // 30		Distortion center shift X for r^6"
        			33, // 31		Distortion center shift Y for r^6"
        			34, // 32		Orthogonal elongation for r^6"
        			35, // 33		Diagonal   elongation for r^6"

        			36, // 34		Distortion center shift X for r^7"
        			37, // 35		Distortion center shift Y for r^7"
        			38, // 36		Orthogonal elongation for r^6"
        			39, // 37		Diagonal   elongation for r^6"

        			40, // 38		Distortion center shift X for r^8"
        			41, // 29		Distortion center shift Y for r^8"
        			42, // 40		Orthogonal elongation for r^8"
        			43  // 41		Diagonal   elongation for r^8"

        	};
        	double [][] result = new double [order.length][2];
        	for (int i=0; i<order.length;i++){
        		result[i][0]=srcDerivatives[order[i]][0];
        		result[i][1]=srcDerivatives[order[i]][1];
        	}
        	return result;
        }
/**
 * Calculate lens center location in target coordinates
 * @return lens center coordinates
 */
        public double [] getLensCenterCoordinates(){
        	double [] p={
        			 this.x0-this.distance*this.rotMatrix[2][0],
        			-this.y0-this.distance*this.rotMatrix[2][1], // this.y0 - up?
        			 this.z0-this.distance*this.rotMatrix[2][2]};
/*
        	Matrix MR=new Matrix(this.rotMatrix);
        	Matrix MRT=MR.transpose();
        	Matrix E=MR.times(MRT);
        	System.out.println("===MR==");  	MR.print(8, 5);System.out.println("");
        	System.out.println("===MRT==");  	MRT.print(8, 5);System.out.println("");
        	System.out.println("===E===");  	E.print(8, 5);System.out.println("");

        	System.out.println("x0="+IJ.d2s(this.x0,1)+" y0="+IJ.d2s(this.y0,1)+" z0="+IJ.d2s(this.z0,1)+" this.distance="+IJ.d2s(this.distance,1));
        	System.out.println("phi="+IJ.d2s(this.phi,4)+" theta="+IJ.d2s(this.theta,4)+" psi="+IJ.d2s(this.psi,4));

        	System.out.println("|"+IJ.d2s(this.rotMatrix[0][0],4)+" "+IJ.d2s(this.rotMatrix[0][1],4)+" "+IJ.d2s(this.rotMatrix[0][2],4)+"|");
        	System.out.println("|"+IJ.d2s(this.rotMatrix[1][0],4)+" "+IJ.d2s(this.rotMatrix[1][1],4)+" "+IJ.d2s(this.rotMatrix[1][2],4)+"|");
        	System.out.println("|"+IJ.d2s(this.rotMatrix[2][0],4)+" "+IJ.d2s(this.rotMatrix[2][1],4)+" "+IJ.d2s(this.rotMatrix[2][2],4)+"|");
*/
        	return p;
        }

/*
 * Calculate pixel value and partial derivatives for different parameters
 * of projection from the pattern point [xp,yp,zp] using current distortion/position parameters
 * [0][0] - pixel x - NaN if looking away
 * [1][ 0] - pixel y- NaN if looking away
 * [*][ 1] - pixel x[or y] partial derivative for phi    (yaw)
 * [*][ 2] - pixel x[or y] partial derivative for theta  (pitch)
 * [*][ 3] - pixel x[or y] partial derivative for psi    (roll)
 * [*][ 4] - pixel x[or y] partial derivative for X0 (intersection of the lens axis with zp==Z0 plane of the target)
 * [*][ 5] - pixel x[or y] partial derivative for Y0 (intersection of the lens axis with zp==Z0 plane of the target)
 * [*][ 6] - pixel x[or y] partial derivative for Z0 (intersection of the lens axis with zp==Z0 plane of the target) - not used
 * [*][ 7] - pixel x[or y] partial derivative for f  (focal length)
 * [*][ 8] - pixel x[or y] partial derivative for dist (distance from [X0,Y0,Z0] to the lens entrance pupil
 * [*][ 9] - pixel x[or y] partial derivative for Da8 (distortion coefficient for r^8)
 * [*][10] - pixel x[or y] partial derivative for Da7 (distortion coefficient for r^7)
 * [*][11] - pixel x[or y] partial derivative for Da6 (distortion coefficient for r^6)
 * [*][12] - pixel x[or y] partial derivative for Da5 (distortion coefficient for r^5)
 * [*][13] - pixel x[or y] partial derivative for Da (distortion coefficient for r^4)
 * [*][14] - pixel x[or y] partial derivative for Db (distortion coefficient for r^3)
 * [*][15] - pixel x[or y] partial derivative for Dc (distortion coefficient for r^2)
 * [*][16] - pixel x[or y] partial derivative for PX0 (lens axis on the sensor, horizontal, right,in pixels)
 * [*][17] - pixel x[or y] partial derivative for PY0 (lens axis on the sensor, vertical, down, in pixels)
 * [*][18] - Orthogonal elongation for r^2"
 * [*][19] - Diagonal   elongation for r^2"
 * [*][20] - Distortion center shift X for r^3"
 * [*][21] - Distortion center shift Y for r^3"
 * [*][22] - Orthogonal elongation for r^3"
 * [*][23] - Diagonal   elongation for r^3"
 * [*][24] - Distortion center shift X for r^4"
 * [*][25] - Distortion center shift Y for r^4"
 * [*][26] - Orthogonal elongation for r^4"
 * [*][27] - Diagonal   elongation for r^4"
 * [*][28] - Distortion center shift X for r^5"
 * [*][29] - Distortion center shift Y for r^5"
 * [*][30] - Orthogonal elongation for r^5"
 * [*][31] - Diagonal   elongation for r^5"
 * [*][32] - Distortion center shift X for r^6"
 * [*][33] - Distortion center shift Y for r^6"
 * [*][34] - Orthogonal elongation for r^6"
 * [*][35] - Diagonal   elongation for r^6"
 * [*][36] - Distortion center shift X for r^7"
 * [*][37] - Distortion center shift Y for r^7"
 * [*][38] - Orthogonal elongation for r^7"
 * [*][39] - Diagonal   elongation for r^7"
 * [*][40] - Distortion center shift X for r^8"
 * [*][41] - Distortion center shift Y for r^8"
 * [*][42] - Orthogonal elongation for r^8"
 * [*][43] - Diagonal   elongation for r^8"
 */
        /*
         * TODO: minimaize calculations for individual {xp,yp,zp}
         */
        public double ipow(double a, int b){
        	switch (b) {
        	case 1: return a;
        	case 2: return a*a;
        	case 3: return a*a*a;
        	case 0: return 1.0;
        	default:
        		if (b<0) {
        			return 1.0/ipow(a,-b);
        		} else {
        			int b1=b>>1;
        		    return ipow(a,b1)*ipow(a,b-b1);
        		}
        	}
        }
        public double [][] calcPartialDerivatives(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp, // target point horizontal, positive - away from camera,  mm
        		boolean calculateAll){ // calculate derivatives, false - values only
        	double maxRelativeRadius=2.0;
        	return calcPartialDerivatives(xp,yp,zp,maxRelativeRadius,calculateAll);
        }
        public double [][] calcPartialDerivatives(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp, // target point horizontal, positive - away from camera,  mm
            	double maxRelativeRadius, // make configurable?
        		boolean calculateAll){ // calculate derivatives, false - values only
        	switch (this.lensDistortionModelType){
        	case 0:
        	case 1:
        		return calcPartialDerivatives_type1(xp,yp,zp,maxRelativeRadius,calculateAll);
        	case 2:
        		return calcPartialDerivatives_type2(xp,yp,zp,maxRelativeRadius,calculateAll);
        	default:
        		return calcPartialDerivatives_type1(xp,yp,zp,maxRelativeRadius,calculateAll);
        	}
        }


        public double [][] calcPartialDerivatives_type1(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp, // target point horizontal, positive - away from camera,  mm
        		double maxRelativeRadius,
        		boolean calculateAll){ // calculate derivatives, false - values only
//        	this.cummulativeCorrection=false; // just debugging

        	// TODO - add reduced calculations for less terms?
//        	final int numDerivatives=44; // 18+6*2+7*2; // 18  for radial and 26 more for non-radial
        	final int numRadialDerivatives=18;
        	final int numNonRadialDerivatives=26;
        	final int distortionCIndex=15; // a8->9, a7->10, a6->11, a5->12, a->13, b->14, c->15
        	double partDeriv[][] = new double [calculateAll?(numRadialDerivatives+numNonRadialDerivatives):1][2];
//        	double [] XYZ= {xp-this.x0, yp-this.y0, zp-this.z0};
        	double [] XYZ= {xp-this.x0, yp+this.y0, zp-this.z0}; // Y - "down" (as in images), not up
        	double [] XeYeZe={
        			this.rotMatrix[0][0]*XYZ[0] + this.rotMatrix[0][1]*XYZ[1] + this.rotMatrix[0][2]*XYZ[2],
        			this.rotMatrix[1][0]*XYZ[0] + this.rotMatrix[1][1]*XYZ[1] + this.rotMatrix[1][2]*XYZ[2],
        			this.rotMatrix[2][0]*XYZ[0] + this.rotMatrix[2][1]*XYZ[1] + this.rotMatrix[2][2]*XYZ[2]+this.distance
        	};
        	// PXYmmc - pinhole (non-distorted) projection coordinates. metric (in mm)
        	double [] PXYmmc={this.focalLength/XeYeZe[2]*XeYeZe[0],this.focalLength/XeYeZe[2]*XeYeZe[1]};
        	// now each term has individual radius
//        	double [] rr=new double [r_xyod.length];
        	// Geometric - get to pinhole coordinates on the sensor
            double [][] dXeYeZe=null; //[14];
            double [][] dPXYmmc=null;

    		/*
   		 Modifying to accommodate for eccentricity of different terms (2 parameters per term) and elliptical shape (another 2 terms). When all are
   		 zeroes, the old parameters are in effect:
   				Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
   				Rdist/R=A8*R6^7+A7*R5^6+A6*R4^5+A5*R3^4+A*R2^3+B*R1^2+C*R0+(1-A6-A7-A6-A5-A-B-C)");
   				R[i] depends on px,py and r_xy[i][], r_o[i] (positive - "landscape", negative - "portrait"), r_d (positive - along y=x, negative - along y=-x)
   				R[i] = r0[i]*(1+  (r_od[i][0]*(y[i]**2-x[i]**2)+ 2*r_od[i][1]*x[i]*y[i])/r0[i]**2;
   				r0[i]=sqrt(x[i]**2+y[2]**2)
   				x[i]=pixel_x-px0-((i>0)?r_xy[i-1][0]:0)
   				y[i]=pixel_y-py0-((i>0)?r_xy[i-1][1]:0)
   		*/
        	double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
        	double deltaX=0.0,deltaY=0.0; // difference between distorted and non-distorted as a fraction of this.distortionRadius
        	double xmmc=PXYmmc[0]/this.distortionRadius;
        	double ymmc=PXYmmc[1]/this.distortionRadius;
        	double dDeltaX_dxmmc=0.0, dDeltaX_dymmc=0.0, dDeltaY_dxmmc=0.0, dDeltaY_dymmc=0;
        	if (calculateAll){
            	for (double [] r:partDeriv) Arrays.fill(r, 0.0); // will fill the [0] (value) too, should be done before calculation of these 2 values
            	// Geometric - get to pinhole coordinates on the sensor
                dXeYeZe=new double[9][3]; //[14];
             // Geometric - get to pinhole coordinates on the sensor
                dXeYeZe[0][0]=XeYeZe[0];
                dXeYeZe[0][1]=XeYeZe[1];
                dXeYeZe[0][2]=XeYeZe[2];
    // /dphi
                dXeYeZe[1][0]=(-cPS*sPH + sPS*sTH*cPH)*XYZ[0]                    +(-cPS*cPH-sPS*sTH*sPH)*XYZ[2];
                dXeYeZe[1][1]=(-sPS*sPH-cPS*sTH*cPH)  *XYZ[0]                    +(-sPS*cPH+cPS*sTH*sPH)*XYZ[2];
                dXeYeZe[1][2]=(cTH*cPH)*XYZ[0]                                     -cTH*sPH* XYZ[2];
    // /dtheta
                dXeYeZe[2][0]=(sPS*cTH*sPH)*XYZ[0]              -sPS*sTH*XYZ[1]  +(sPS*cTH*cPH)*XYZ[2];
                dXeYeZe[2][1]=(-cPS*cTH*sPH)*XYZ[0]             +cPS*sTH*XYZ[1]  +(-cPS*cTH*cPH)*XYZ[2];
                dXeYeZe[2][2]=(-sTH*sPH)*XYZ[0]                 -cTH*XYZ[1]      -sTH*cPH* XYZ[2];
    // /dpsi
                dXeYeZe[3][0]=(-sPS*cPH+cPS*sTH*sPH)*XYZ[0]     +cPS*cTH*XYZ[1]  +(sPS*sPH+cPS*sTH*cPH)*XYZ[2];
                dXeYeZe[3][1]=(cPS*cPH+ sPS*sTH*sPH)*XYZ[0]     +sPS*cTH*XYZ[1]  +(-cPS*sPH+sPS*sTH*cPH)*XYZ[2];
                dXeYeZe[3][2]=0.0;
    // /dX0
//              dXeYeZe[4][0]=-cPS*cPH+sPS*sTH*sPH; // bad?
                dXeYeZe[4][0]=-cPS*cPH-sPS*sTH*sPH; // bad?
                dXeYeZe[4][1]=-sPS*cPH+cPS*sTH*sPH;
                dXeYeZe[4][2]=-cTH*sPH;
    // /dY0
//                dXeYeZe[5][0]=-sPS*cTH;
//                dXeYeZe[5][1]= cPS*cTH;
//                dXeYeZe[5][2]= sTH;
                dXeYeZe[5][0]=+sPS*cTH;
                dXeYeZe[5][1]=-cPS*cTH;
                dXeYeZe[5][2]=-sTH;
    // /dZ0
//              dXeYeZe[6][0]= cPS*sPH+sPS*sTH*cPH; //bad?
                dXeYeZe[6][0]= cPS*sPH-sPS*sTH*cPH; //bad?
                dXeYeZe[6][1]= sPS*sPH+cPS*sTH*cPH;
                dXeYeZe[6][2]=-cTH*cPH;
    // /df
                dXeYeZe[7][0]=0.0;
                dXeYeZe[7][1]=0.0;
                dXeYeZe[7][2]=0.0;
    // /ddist
                dXeYeZe[8][0]=0.0;
                dXeYeZe[8][1]=0.0;
                dXeYeZe[8][2]=1.0;


                dPXYmmc=new double[9][2]; //[14];
             // TODO: move up
//              double [][] dPXYmmc=new double[9][2]; //[14];
              dPXYmmc[0][0]=PXYmmc[0]; // is it needed? probably not used
              dPXYmmc[0][1]=PXYmmc[1];
  //(4) PXmmc =f/Ze*  Xe  // mm, left from the lens axis intersection with the sensor

  //dPXmmc/dphi   = f/Ze * dXe/dphi - f*Xe/Ze^2 * dZe/dphi
              dPXYmmc[1][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[1][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[1][2]);
  //dPXmmc/dtheta = f/Ze * dXe/dtheta - f*Xe/Ze^2 * dZe/dtheta
              dPXYmmc[2][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[2][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[2][2]);
  //dPXmmc/dpsi   = f/Ze * dXe/dpsi - f*Xe/Ze^2 * dZe/dpsi
              dPXYmmc[3][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[3][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[3][2]);
  //dPXmmc/dX0    = f/Ze * dXe/dX0 - f*Xe/Ze^2 * dZe/dX0
              dPXYmmc[4][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[4][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[4][2]);
  //dPXmmc/dY0    = f/Ze * dXe/dY0 - f*Xe/Ze^2 * dZe/dY0
              dPXYmmc[5][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[5][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[5][2]);
  //dPXmmc/dZ0    = f/Ze * dXe/dZ0 - f*Xe/Ze^2 * dZe/dZ0
              dPXYmmc[6][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[6][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[6][2]); //bad?
  //dPXmmc/df     = Xe/Ze
              dPXYmmc[7][0]=dXeYeZe[0][0]/dXeYeZe[0][2];
  //dPXmmc/ddist =  - f*Xe/Ze^2
              dPXYmmc[8][0]=-this.focalLength*dXeYeZe[0][0]/(dXeYeZe[0][2]*dXeYeZe[0][2]);

  //(5) PYmmc =f/Ze*  Ye  // mm, up from the lens axis intersection with the sensor
  //dPYmmc/dphi   = f/Ze * dYe/dphi - f*Ye/Ze^2 * dZe/dphi
              dPXYmmc[1][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[1][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[1][2]);
  //dPYmmc/dtheta = f/Ze * dYe/dtheta - f*Ye/Ze^2 * dZe/dtheta
              dPXYmmc[2][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[2][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[2][2]);
  //dPYmmc/dpsi   = f/Ze * dYe/dpsi - f*Ye/Ze^2 * dZe/dpsi
              dPXYmmc[3][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[3][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[3][2]);
  //dPYmmc/dX0    = f/Ze * dYe/dX0 - f*Ye/Ze^2 * dZe/dX0
              dPXYmmc[4][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[4][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[4][2]);
  //dPYmmc/dY0    = f/Ze * dYe/dY0 - f*Ye/Ze^2 * dZe/dY0
              dPXYmmc[5][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[5][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[5][2]);
  //dPYmmc/dZ0    = f/Ze * dYe/dZ0 - f*Ye/Ze^2 * dZe/dZ0
              dPXYmmc[6][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[6][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[6][2]); // good?
  //dPYmmc/df     = Ye/Ze
              dPXYmmc[7][1]=dXeYeZe[0][1]/dXeYeZe[0][2];
  //dPYmmc/ddist =  - f*Ye/Ze^2
              dPXYmmc[8][1]=-this.focalLength*dXeYeZe[0][1]/(dXeYeZe[0][2]*dXeYeZe[0][2]);

        	}
        	// conversion coefficient from relative (to distortionRadius) to pixels
        	// negate for y!
        	double rel_to_pix=this.distortionRadius*1000.0/this.pixelSize;
        	//TODO:  seems that rr[i] can be just a single running variable, not an array
        	for (int i=0;i<r_xyod.length;i++){

        		double x=xmmc-r_xyod[i][0]; // relative X-shift of this term center
        		double y=ymmc-r_xyod[i][1]; // relative X-shift of this term center
        		double x2=x*x;   // relative squared X-shift from this term center
        		double y2=y*y;   // relative squared Y-shift from this term center
        		double r2=x2+y2; // relative squared distance from this term center
        		// effective distance from this term center corrected for elongation
        		double rr=Math.sqrt(r2+ r_xyod[i][2]*(y2-x2)+ 2.0*r_xyod[i][3]*x*y );
                if (rr<0.00000001*this.distortionRadius) rr=0.00000001*this.distortionRadius; // avoid 1/0.0 -// Added 3 zeros (old comment)

        		double rr_pow_i=ipow(rr,i);
//        		double ki=a[i]*(Math.pow(rr,i+1)-1.0);
//        		double ki=a[i]*(ipow(rr,i+1)-1.0);
        		double ki=a[i]*(rr_pow_i*rr-1.0); // scaling of distorted distance from this term's center (does not include pinhole center shift itself)
        		deltaX+=x*ki; // relative distorted distance from the center
        		deltaY+=y*ki;
//        		if ((debugLevel>2) && ((r_xyod[i][0]!=0.0) || (r_xyod[i][0]!=0.0))){
//        			System.out.println ("i="+i+" r_xyod[i][0]="+r_xyod[i][0]+" r_xyod[i][1]="+r_xyod[i][1]+" a[0]="+a[0]+" a[1]="+a[1]+" a[2]="+a[2]+" a[3]="+a[3]);
//        		}
        		if (calculateAll) {
        			double csi=rel_to_pix*a[i]*(i+1)*rr_pow_i;

//        		double dx2_dxmmc=2*x, dy2_dymmc=2*y, dr2_dxmmc=2*x, dr2_dymmc=2*y;
// double drr_dxmmc=0.5/rr*(dr2_dxmmc+ r_xyod[i][2]*(-dx2_dxmmc)+2*r_xyod[i][3]*y)=0.5/rr*(2*x - 2*x*r_xyod[i][2]+2*r_xyod[i][3]*y)=
//                  =  (x*(1.0-r_xyod[i][2])+y*r_xyod[i][3])/rr
// double drr_dymmc=  (y*(1.0+r_xyod[i][2])+x*r_xyod[i][3])/rr
// double ki=a[i]*(rr^(i+1)-1.0);
// double dki_dxmmc=a[i]*(i+1)*rr^i*(x*(1.0-r_xyod[i][2])+y*r_xyod[i][3])/rr
// double dki_dymmc=a[i]*(i+1)*rr^i*(y*(1.0+r_xyod[i][2])+x*r_xyod[i][3])/rr

// double dki_dxmmc=a[i]*(i+1)*rr_pow_i*(x*(1.0-r_xyod[i][2])+y*r_xyod[i][3])/rr
// double dki_dymmc=a[i]*(i+1)*rr_pow_i*(y*(1.0+r_xyod[i][2])+x*r_xyod[i][3])/rr

        			double ai_iplus1_rr_pow_i= a[i]*(i+1)*rr_pow_i;
        			double dki_dxmmc=          ai_iplus1_rr_pow_i*(x*(1.0-r_xyod[i][2])+y*r_xyod[i][3])/rr;
        			double dki_dymmc=          ai_iplus1_rr_pow_i*(y*(1.0+r_xyod[i][2])+x*r_xyod[i][3])/rr;

                    // the following 4 shifts are "extra", on top of non-distorted (pinhole) - pinhole should be added to the per-term sum
        			double dDeltaXi_dxmmc = ki+x*dki_dxmmc; // here dDelta*_d*mmc are relative (both I/O are fractions of distortionRadius)
        			double dDeltaXi_dymmc =    x*dki_dymmc;
        			double dDeltaYi_dxmmc =    y*dki_dxmmc;
        			double dDeltaYi_dymmc = ki+y*dki_dymmc;
        			dDeltaX_dxmmc += dDeltaXi_dxmmc; // here dDelta*_d*mmc are relative (both I/O are fractions of distortionRadius)
        			dDeltaX_dymmc += dDeltaXi_dymmc;
        			dDeltaY_dxmmc += dDeltaYi_dxmmc;
        			dDeltaY_dymmc += dDeltaYi_dymmc;
// dependence on eccentricity
// dDeltaX_d_r_xyod0 =  dDeltaXi_dxmmc * dxmmc_d_r_xyod0 = -dDeltaXi_dxmmc
// dDeltaY_d_r_xyod0 =  dDeltaYi_dxmmc * dxmmc_d_r_xyod0 = -dDeltaYi_dxmmc
// dDeltaX_d_r_xyod1 =  dDeltaXi_dymmc * dymmc_d_r_xyod1 = -dDeltaXi_dymmc
// dDeltaY_d_r_xyod1 =  dDeltaYi_dymmc * dymmc_d_r_xyod1 = -dDeltaYi_dymmc
        			int index=numRadialDerivatives-2+4*i;
        			if (i>0){ // eccentricity is not applicabe to the first (C) term
        				partDeriv[index  ][0]= -rel_to_pix*(dDeltaXi_dxmmc-0.0); // dPx_dr_xyod0
        				partDeriv[index  ][1]=  rel_to_pix* dDeltaYi_dxmmc; // dPy_dr_xyod0
        				partDeriv[index+1][0]= -rel_to_pix* dDeltaXi_dymmc; // dPx_dr_xyod1
        				partDeriv[index+1][1]=  rel_to_pix*(dDeltaYi_dymmc-0.0); // dPy_dr_xyod1
        			}

        			// d/dai
        			// ki=a[i]*(rr_pow_i*rr-1.0);
        			// double dki_dai=(rr_pow_i*rr-1.0);
        			// double dpx_dai= rel_to_pix*x*(dki_dai)=rel_to_pix*x*(rr_pow_i*rr-1.0)
        			// double dpy_dai=-rel_to_pix*y*(dki_dai)=rel_to_pix*y*(rr_pow_i*rr-1.0)
        			//        	final int distortionCIndex=15; // a8->9, a7->10, a6->11, a5->12, a->13, b->14, c->15
        			index= distortionCIndex-i; // reverse order, // a8->9, a7->10, a6->11, a5->12, a->13, b->14, c->15
    				partDeriv[index][0]= rel_to_pix*x*(rr_pow_i*rr-1.0); // OK
    				partDeriv[index][1]=-rel_to_pix*y*(rr_pow_i*rr-1.0); // OK

        			// d/dr_xyod[0] (x shift of the center
        			// rr=Math.sqrt(r2+ r_xyod[i][2]*(y2-x2)+ 2.0*r_xyod[i][3]*x*y );
        			// ki=a[i]*(rr_pow_i*rr-1.0);
        			// dx_dr_xyod0=-1.0
        			// dy_dr_xyod1=-1.0
        			// dx2_dr_xyod0=-2*x
        			// dy2_dr_xyod1=-2*y
        			// dr2_dr_xyod0=-2*x
        			// dr2_dr_xyod1=-2*y
        			// dxy_dr_xyod0=-y
        			// dxy_dr_xyod1=-x
        			// dri_dr_xyod[0]=0.5*(dr2_dr_xyod0 -r_xyod[i][2]*dx2_dr_xyod0 + 2.0*r_xyod[i][3]*dxy_dr_xyod0)/rr=
        			//                0.5*(-2*x -r_xyod[i][2]*(-2*x) + 2.0*r_xyod[i][3]*(-y))/rr=
        			//                (-x +r_xyod[i][2]*x -r_xyod[i][3]*y)/rr=
        			//                (x*(r_xyod[i][2]-1.0) -r_xyod[i][3]*y)/rr
//        			double dri_dr_xyod0 = (x*(r_xyod[i][2]-1.0) -r_xyod[i][3]*y)/rr;
        			// dri_dr_xyod[1]=0.5*(dr2_dr_xyod1 -r_xyod[i][2]*dx2_dr_xyod1 + 2.0*r_xyod[i][3]*dxy_dr_xyod1)/rr=
        			//                0.5*(-2*y +r_xyod[i][2]*(-2*y) + 2.0*r_xyod[i][3]*(-x))/rr=
        			//                (-y -r_xyod[i][2]*y -r_xyod[i][3]*x)/rr=
        			//                (-y*(r_xyod[i][2]+1.0) -r_xyod[i][3]*x)/rr
//        			double dri_dr_xyod1 = (-y*(r_xyod[i][2]+1.0) -r_xyod[i][3]*x)/rr;
        			// dri_dr_xyod[2]=0.5*(y2-x2)/rr
        			double dri_dr_xyod2 = 0.5*(y2-x2)/rr;
        			// dri_dr_xyod[3]=(y*x)/rr
        			double dri_dr_xyod3 = (y*x)/rr;
        			//--------------------
        			// double ki=a[i]*(rr^(i+1)-1.0);
        			// double ki=a[i]*(rr_pow_i*rr-1.0);
        			// deltaX+=x*ki;
        			// deltaY+=y*ki;
        			// dki_dri=a[i]*(i+1)*rr^i = a[i]*(i+1)*rr_pow_i
        			// dDeltaX_dri=x*a[i]*(i+1)*rr_pow_i
        			// dDeltaY_dri=y*a[i]*(i+1)*rr_pow_i
        			// dPx_dri= rel_to_pix*x*a[i]*(i+1)*rr_pow_i =  csi*x
        			// dPy_dri=-rel_to_pix*y*a[i]*(i+1)*rr_pow_i = -csi*y
        			// dPx_dr_xyod0= dPx_dri*dri_dr_xyod0= csi*x*dri_dr_xyod0;
        			// dPx_dr_xyod1= dPx_dri*dri_dr_xyod1= csi*x*dri_dr_xyod1;
        			// dPx_dr_xyod2= dPx_dri*dri_dr_xyod2= csi*x*dri_dr_xyod2;
        			// dPx_dr_xyod3= dPx_dri*dri_dr_xyod3= csi*x*dri_dr_xyod3;
        			// dPy_dr_xyod0= dPy_dri*dri_dr_xyod0=-csi*y*dri_dr_xyod0;
        			// dPy_dr_xyod1= dPy_dri*dri_dr_xyod1=-csi*y*dri_dr_xyod1;
        			// dPy_dr_xyod2= dPy_dri*dri_dr_xyod2=-csi*y*dri_dr_xyod2;
        			// dPy_dr_xyod3= dPy_dri*dri_dr_xyod3=-csi*y*dri_dr_xyod3;
        			index=numRadialDerivatives-2+4*i;
        			/*
        			if (i>0){ // eccentricity is not applicabe to the first (C) term
        				partDeriv[index  ][0]= csi*x*dri_dr_xyod0; // dPx_dr_xyod0
        				partDeriv[index  ][1]=-csi*y*dri_dr_xyod0; // dPy_dr_xyod0
        				partDeriv[index+1][0]= csi*x*dri_dr_xyod1; // dPx_dr_xyod1
        				partDeriv[index+1][1]=-csi*y*dri_dr_xyod1; // dPy_dr_xyod1
        			}
        			*/
    				partDeriv[index+2][0]= csi*x*dri_dr_xyod2; // dPx_dr_xyod2
    				partDeriv[index+2][1]=-csi*y*dri_dr_xyod2; // dPy_dr_xyod2
    				partDeriv[index+3][0]= csi*x*dri_dr_xyod3; // dPx_dr_xyod3
    				partDeriv[index+3][1]=-csi*y*dri_dr_xyod3; // dPy_dr_xyod3

        		}
        	}

        	double [] xyDist={PXYmmc[0]+this.distortionRadius*deltaX,PXYmmc[1]+this.distortionRadius*deltaY};
        	// convert to sensor pixels coordinates
        	partDeriv[0][0]=  1000.0/this.pixelSize*xyDist[0] + this.px0;
        	partDeriv[0][1]= -1000.0/this.pixelSize*xyDist[1] + this.py0;

        	if (!calculateAll) { // TODO: how to deal with it when calculating Jacobian???
        		// TODO: Looking away from the target, trying only with no dervatives. Do the same for derivatives too?
            	if ((XeYeZe[2]<0.0) || ((xmmc*xmmc+ymmc*ymmc)>maxRelativeRadius*maxRelativeRadius)){ // non-distorted too far from the axis
            		partDeriv[0][0]=Double.NaN;
            		partDeriv[0][1]=Double.NaN;
            	}
        		return partDeriv;
        	}
        	// correct parameter derivatives to cumulative version
        	if (this.cummulativeCorrection){
            	for (int i=r_xyod.length-2;i>=0;i--){
        			int index=numRadialDerivatives-2+4*i;
        			if (i>0){ // eccentricity is not applicabe to the first (C) term
        				partDeriv[index  ][0]+=partDeriv[index+4][0]; // oob=36
        				partDeriv[index  ][1]+=partDeriv[index+4][1];
        				partDeriv[index+1][0]+=partDeriv[index+5][0];
        				partDeriv[index+1][1]+=partDeriv[index+5][1];
        			}
    				partDeriv[index+2][0]+=partDeriv[index+6][0];
    				partDeriv[index+2][1]+=partDeriv[index+6][1];
    				partDeriv[index+3][0]+=partDeriv[index+7][0];
    				partDeriv[index+3][1]+=partDeriv[index+7][1];
            	}
        	}

// convert    dDelta*_d*mmc from relative/relative to pix/mm (invert pixel Y direction)
// added 1.0 to account for non-distorted (pinhole) shift
        	double dPx_dPinholeX= (1.0+dDeltaX_dxmmc)*1000.0/this.pixelSize;
        	double dPx_dPinholeY=       dDeltaX_dymmc*1000.0/this.pixelSize;
        	double dPy_dPinholeX=      -dDeltaY_dxmmc*1000.0/this.pixelSize;
        	double dPy_dPinholeY=-(1.0+dDeltaY_dymmc)*1000.0/this.pixelSize;
        	if (this.debugLevel>2){
        		System.out.println(" deltaX="+deltaX+" deltaY="+deltaY+" xyDist[0]="+xyDist[0]+" xyDist[1]="+xyDist[1]);
        		System.out.println(" PXYmmc[0]="+PXYmmc[0]+" PXYmmc[1]="+PXYmmc[1]+" xmmc="+xmmc+" ymmc"+ymmc);
        		System.out.println(" dDeltaX_dxmmc="+dDeltaX_dxmmc+" dDeltaX_dymmc="+dDeltaX_dymmc+" dDeltaY_dxmmc="+dDeltaY_dxmmc+" dDeltaY_dymmc"+dDeltaY_dymmc);
        		System.out.println(" dPx_dPinholeX="+dPx_dPinholeX+" dPx_dPinholeY="+dPx_dPinholeY+" dPy_dPinholeX="+dPy_dPinholeX+" dPy_dPinholeY"+dPy_dPinholeY);
        	}

            double K=Math.PI/180; // multiply all derivatives my angles
// dPX/dphi   =  1000/Psz* dxDist/dphi
        	partDeriv[ 1][0]=  K*(dPx_dPinholeX*dPXYmmc[1][0]+dPx_dPinholeY*dPXYmmc[1][1]);
        	partDeriv[ 1][1]=  K*(dPy_dPinholeX*dPXYmmc[1][0]+dPy_dPinholeY*dPXYmmc[1][1]);

// dPX/dtheta =  1000/Psz* dxDist/dtheta
        	partDeriv[ 2][0]=  K*(dPx_dPinholeX*dPXYmmc[2][0]+dPx_dPinholeY*dPXYmmc[2][1]);
        	partDeriv[ 2][1]=  K*(dPy_dPinholeX*dPXYmmc[2][0]+dPy_dPinholeY*dPXYmmc[2][1]);
// dPX/dpsi   =  1000/Psz* dxDist/dpsi
        	partDeriv[ 3][0]=  K*(dPx_dPinholeX*dPXYmmc[3][0]+dPx_dPinholeY*dPXYmmc[3][1]);
        	partDeriv[ 3][1]=  K*(dPy_dPinholeX*dPXYmmc[3][0]+dPy_dPinholeY*dPXYmmc[3][1]);
// dPX/dX0    =  1000/Psz* dxDist/dX0
        	partDeriv[ 4][0]=  dPx_dPinholeX*dPXYmmc[4][0]+dPx_dPinholeY*dPXYmmc[4][1];
        	partDeriv[ 4][1]=  dPy_dPinholeX*dPXYmmc[4][0]+dPy_dPinholeY*dPXYmmc[4][1];
// dPX/dY0    =  1000/Psz* dxDist/dY0
        	partDeriv[ 5][0]=  dPx_dPinholeX*dPXYmmc[5][0]+dPx_dPinholeY*dPXYmmc[5][1];
        	partDeriv[ 5][1]=  dPy_dPinholeX*dPXYmmc[5][0]+dPy_dPinholeY*dPXYmmc[5][1];
// dPX/dZ0    =  1000/Psz* dxDist/dZ0
        	partDeriv[ 6][0]=  dPx_dPinholeX*dPXYmmc[6][0]+dPx_dPinholeY*dPXYmmc[6][1];
        	partDeriv[ 6][1]=  dPy_dPinholeX*dPXYmmc[6][0]+dPy_dPinholeY*dPXYmmc[6][1];
// dPX/df     =  1000/Psz* dxDist/df
        	partDeriv[ 7][0]=  dPx_dPinholeX*dPXYmmc[7][0]+dPx_dPinholeY*dPXYmmc[7][1];
        	partDeriv[ 7][1]=  dPy_dPinholeX*dPXYmmc[7][0]+dPy_dPinholeY*dPXYmmc[7][1];
// dPX/ddist  =  1000/Psz* dxDist/ddist
        	partDeriv[ 8][0]=  dPx_dPinholeX*dPXYmmc[8][0]+dPx_dPinholeY*dPXYmmc[8][1];
        	partDeriv[ 8][1]=  dPy_dPinholeX*dPXYmmc[8][0]+dPy_dPinholeY*dPXYmmc[8][1];

// dPX/dPX0   =  1
// dPY/dPX0   =  0
        	partDeriv[16][0]=  1.0;
        	partDeriv[16][1]=  0.0;
// dPX/dPY0   =  0
// dPY/dPY0   =  1
        	partDeriv[17][0]=  0.0;
        	partDeriv[17][1]=  1.0;
            return partDeriv;
        }

        public double [][] calcPartialDerivatives_type2(
        		double xp, // target point horizontal, positive - right,  mm
        		double yp, // target point vertical,   positive - down,  mm
        		double zp, // target point horizontal, positive - away from camera,  mm
        		double maxRelativeRadius,
        		boolean calculateAll){ // calculate derivatives, false - values only
//        	this.cummulativeCorrection=false; // just debugging

        	// TODO - add reduced calculations for less terms?
//        	final int numDerivatives=44; // 18+6*2+7*2; // 18  for radial and 26 more for non-radial
        	final int numRadialDerivatives=18;
        	final int numNonRadialDerivatives=26;
        	final int distortionCIndex=15; // a8->9, a7->10, a6->11, a5->12, a->13, b->14, c->15
        	double partDeriv[][] = new double [calculateAll?(numRadialDerivatives+numNonRadialDerivatives):1][2];
//        	double [] XYZ= {xp-this.x0, yp-this.y0, zp-this.z0};
        	double [] XYZ= {xp-this.x0, yp+this.y0, zp-this.z0}; // Y - "down" (as in images), not up
        	double [] XeYeZe={
        			this.rotMatrix[0][0]*XYZ[0] + this.rotMatrix[0][1]*XYZ[1] + this.rotMatrix[0][2]*XYZ[2],
        			this.rotMatrix[1][0]*XYZ[0] + this.rotMatrix[1][1]*XYZ[1] + this.rotMatrix[1][2]*XYZ[2],
        			this.rotMatrix[2][0]*XYZ[0] + this.rotMatrix[2][1]*XYZ[1] + this.rotMatrix[2][2]*XYZ[2]+this.distance
        	};
        	// PXYmmc - pinhole (non-distorted) projection coordinates. metric (in mm)
        	double [] PXYmmc={this.focalLength/XeYeZe[2]*XeYeZe[0],this.focalLength/XeYeZe[2]*XeYeZe[1]};
        	// now each term has individual radius
//        	double [] rr=new double [r_xyod.length];
        	// Geometric - get to pinhole coordinates on the sensor
            double [][] dXeYeZe=null; //[14];
            double [][] dPXYmmc=null;

/*
   		 Modifying to accommodate for eccentricity of different terms (2 parameters per term) and elliptical shape (another 2 terms). When all are
   		 zeroes, the old parameters are in effect:
   				Rdist/R=A8*R^7+A7*R^6+A6*R^5+A5*R^4+A*R^3+B*R^2+C*R+(1-A6-A7-A6-A5-A-B-C)");
   				Rdist/R=A8*R6^7+A7*R5^6+A6*R4^5+A5*R3^4+A*R2^3+B*R1^2+C*R0+(1-A6-A7-A6-A5-A-B-C)");
   				R[i] depends on px,py and r_xy[i][], r_o[i] (positive - "landscape", negative - "portrait"), r_d (positive - along y=x, negative - along y=-x)
   				R[i] = r0[i]*(1+  (r_od[i][0]*(y[i]**2-x[i]**2)+ 2*r_od[i][1]*x[i]*y[i])/r0[i]**2;
   				r0[i]=sqrt(x[i]**2+y[2]**2)
   				x[i]=pixel_x-px0-((i>0)?r_xy[i-1][0]:0)
   				y[i]=pixel_y-py0-((i>0)?r_xy[i-1][1]:0)
   		*/
        	double [] a={this.distortionC,this.distortionB,this.distortionA,this.distortionA5,this.distortionA6,this.distortionA7,this.distortionA8};
        	double deltaX=0.0,deltaY=0.0; // difference between distorted and non-distorted as a fraction of this.distortionRadius
        	double x=PXYmmc[0]/this.distortionRadius; // was xmmc
        	double y=PXYmmc[1]/this.distortionRadius; // was ymmc
        	double dDeltaX_dx=0.0, dDeltaX_dy=0.0, dDeltaY_dx=0.0, dDeltaY_dy=0;
        	if (calculateAll){
            	for (double [] r:partDeriv) Arrays.fill(r, 0.0); // will fill the [0] (value) too, should be done before calculation of these 2 values
            	// Geometric - get to pinhole coordinates on the sensor
                dXeYeZe=new double[9][3]; //[14];
             // Geometric - get to pinhole coordinates on the sensor
                dXeYeZe[0][0]=XeYeZe[0];
                dXeYeZe[0][1]=XeYeZe[1];
                dXeYeZe[0][2]=XeYeZe[2];
    // /dphi
                dXeYeZe[1][0]=(-cPS*sPH + sPS*sTH*cPH)*XYZ[0]                    +(-cPS*cPH-sPS*sTH*sPH)*XYZ[2];
                dXeYeZe[1][1]=(-sPS*sPH-cPS*sTH*cPH)  *XYZ[0]                    +(-sPS*cPH+cPS*sTH*sPH)*XYZ[2];
                dXeYeZe[1][2]=(cTH*cPH)*XYZ[0]                                     -cTH*sPH* XYZ[2];
    // /dtheta
                dXeYeZe[2][0]=(sPS*cTH*sPH)*XYZ[0]              -sPS*sTH*XYZ[1]  +(sPS*cTH*cPH)*XYZ[2];
                dXeYeZe[2][1]=(-cPS*cTH*sPH)*XYZ[0]             +cPS*sTH*XYZ[1]  +(-cPS*cTH*cPH)*XYZ[2];
                dXeYeZe[2][2]=(-sTH*sPH)*XYZ[0]                 -cTH*XYZ[1]      -sTH*cPH* XYZ[2];
    // /dpsi
                dXeYeZe[3][0]=(-sPS*cPH+cPS*sTH*sPH)*XYZ[0]     +cPS*cTH*XYZ[1]  +(sPS*sPH+cPS*sTH*cPH)*XYZ[2];
                dXeYeZe[3][1]=(cPS*cPH+ sPS*sTH*sPH)*XYZ[0]     +sPS*cTH*XYZ[1]  +(-cPS*sPH+sPS*sTH*cPH)*XYZ[2];
                dXeYeZe[3][2]=0.0;
    // /dX0
//              dXeYeZe[4][0]=-cPS*cPH+sPS*sTH*sPH; // bad?
                dXeYeZe[4][0]=-cPS*cPH-sPS*sTH*sPH; // bad?
                dXeYeZe[4][1]=-sPS*cPH+cPS*sTH*sPH;
                dXeYeZe[4][2]=-cTH*sPH;
    // /dY0
//                dXeYeZe[5][0]=-sPS*cTH;
//                dXeYeZe[5][1]= cPS*cTH;
//                dXeYeZe[5][2]= sTH;
                dXeYeZe[5][0]=+sPS*cTH;
                dXeYeZe[5][1]=-cPS*cTH;
                dXeYeZe[5][2]=-sTH;
    // /dZ0
//              dXeYeZe[6][0]= cPS*sPH+sPS*sTH*cPH; //bad?
                dXeYeZe[6][0]= cPS*sPH-sPS*sTH*cPH; //bad?
                dXeYeZe[6][1]= sPS*sPH+cPS*sTH*cPH;
                dXeYeZe[6][2]=-cTH*cPH;
    // /df
                dXeYeZe[7][0]=0.0;
                dXeYeZe[7][1]=0.0;
                dXeYeZe[7][2]=0.0;
    // /ddist
                dXeYeZe[8][0]=0.0;
                dXeYeZe[8][1]=0.0;
                dXeYeZe[8][2]=1.0;


                dPXYmmc=new double[9][2]; //[14];
             // TODO: move up
//              double [][] dPXYmmc=new double[9][2]; //[14];
              dPXYmmc[0][0]=PXYmmc[0]; // is it needed? probably not used
              dPXYmmc[0][1]=PXYmmc[1];
  //(4) PXmmc =f/Ze*  Xe  // mm, left from the lens axis intersection with the sensor

  //dPXmmc/dphi   = f/Ze * dXe/dphi - f*Xe/Ze^2 * dZe/dphi
              dPXYmmc[1][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[1][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[1][2]);
  //dPXmmc/dtheta = f/Ze * dXe/dtheta - f*Xe/Ze^2 * dZe/dtheta
              dPXYmmc[2][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[2][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[2][2]);
  //dPXmmc/dpsi   = f/Ze * dXe/dpsi - f*Xe/Ze^2 * dZe/dpsi
              dPXYmmc[3][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[3][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[3][2]);
  //dPXmmc/dX0    = f/Ze * dXe/dX0 - f*Xe/Ze^2 * dZe/dX0
              dPXYmmc[4][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[4][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[4][2]);
  //dPXmmc/dY0    = f/Ze * dXe/dY0 - f*Xe/Ze^2 * dZe/dY0
              dPXYmmc[5][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[5][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[5][2]);
  //dPXmmc/dZ0    = f/Ze * dXe/dZ0 - f*Xe/Ze^2 * dZe/dZ0
              dPXYmmc[6][0]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[6][0]-dXeYeZe[0][0]/dXeYeZe[0][2]*dXeYeZe[6][2]); //bad?
  //dPXmmc/df     = Xe/Ze
              dPXYmmc[7][0]=dXeYeZe[0][0]/dXeYeZe[0][2];
  //dPXmmc/ddist =  - f*Xe/Ze^2
              dPXYmmc[8][0]=-this.focalLength*dXeYeZe[0][0]/(dXeYeZe[0][2]*dXeYeZe[0][2]);

  //(5) PYmmc =f/Ze*  Ye  // mm, up from the lens axis intersection with the sensor
  //dPYmmc/dphi   = f/Ze * dYe/dphi - f*Ye/Ze^2 * dZe/dphi
              dPXYmmc[1][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[1][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[1][2]);
  //dPYmmc/dtheta = f/Ze * dYe/dtheta - f*Ye/Ze^2 * dZe/dtheta
              dPXYmmc[2][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[2][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[2][2]);
  //dPYmmc/dpsi   = f/Ze * dYe/dpsi - f*Ye/Ze^2 * dZe/dpsi
              dPXYmmc[3][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[3][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[3][2]);
  //dPYmmc/dX0    = f/Ze * dYe/dX0 - f*Ye/Ze^2 * dZe/dX0
              dPXYmmc[4][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[4][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[4][2]);
  //dPYmmc/dY0    = f/Ze * dYe/dY0 - f*Ye/Ze^2 * dZe/dY0
              dPXYmmc[5][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[5][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[5][2]);
  //dPYmmc/dZ0    = f/Ze * dYe/dZ0 - f*Ye/Ze^2 * dZe/dZ0
              dPXYmmc[6][1]=this.focalLength/dXeYeZe[0][2]*(dXeYeZe[6][1]-dXeYeZe[0][1]/dXeYeZe[0][2]*dXeYeZe[6][2]); // good?
  //dPYmmc/df     = Ye/Ze
              dPXYmmc[7][1]=dXeYeZe[0][1]/dXeYeZe[0][2];
  //dPYmmc/ddist =  - f*Ye/Ze^2
              dPXYmmc[8][1]=-this.focalLength*dXeYeZe[0][1]/(dXeYeZe[0][2]*dXeYeZe[0][2]);

        	}
        	// conversion coefficient from relative (to distortionRadius) to pixels
        	// negate for y!
        	double rel_to_pix=this.distortionRadius*1000.0/this.pixelSize;
    		double x2=x*x;   // relative squared X-shift from this term center
    		double y2=y*y;   // relative squared Y-shift from this term center
    		double r2=x2+y2; // relative squared distance from this term center
    		double xy=x*y;
    		double rr=Math.sqrt(r2);
        	//TODO:  seems that rr[i] can be just a single running variable, not an array
        	for (int i=0;i<r_xyod.length;i++){
        		// effective distance from this term center corrected for elongation
//        		double rr=Math.sqrt(r2+ r_xyod[i][2]*(y2-x2)+ 2.0*r_xyod[i][3]*x*y );
                if (rr<0.00000001*this.distortionRadius) rr=0.00000001*this.distortionRadius; // avoid 1/0.0 -// Added 3 zeros (old comment)
        		double rr_pow_i_minus_1,rr_pow_i,rr_pow_i_plus_1;
        		if (i==0){
        			rr_pow_i_minus_1=1.0/rr;
        			rr_pow_i=1.0;
        		}else{
        			rr_pow_i_minus_1=ipow(rr,i-1);
        			rr_pow_i=rr_pow_i_minus_1*rr;
        		}
        		rr_pow_i_plus_1=rr_pow_i*rr;
        		double ki1= r_xyod[i][0]*x+r_xyod[i][1]*y;
        		double ki2= r_xyod[i][2]*(y2-x2)+2.0*r_xyod[i][3]*xy;
        		double ki=a[i]*(rr_pow_i_plus_1-1.0)+ki1*rr_pow_i+ki2*rr_pow_i_minus_1; // scaling of distorted distance from this term's center (does not include pinhole center shift itself)
//        		double ki=a[i]*(rr_pow_i*rr-1.0); // scaling of distorted distance from this term's center (does not include pinhole center shift itself)
        		deltaX+=x*ki; // relative distorted distance from the center
        		deltaY+=y*ki;
//        		if ((debugLevel>2) && ((r_xyod[i][0]!=0.0) || (r_xyod[i][0]!=0.0))){
//        			System.out.println ("i="+i+" r_xyod[i][0]="+r_xyod[i][0]+" r_xyod[i][1]="+r_xyod[i][1]+" a[0]="+a[0]+" a[1]="+a[1]+" a[2]="+a[2]+" a[3]="+a[3]);
//        		}
        		if (calculateAll) {
//        			double csi=rel_to_pix*a[i]*(i+1)*rr_pow_i;
// drr_dx=x/rr, drr_dy=y/r
//                  ki=a[i]*(rr^(i+1-1.0) +ki1*rr^i+ki2*rr^(i-1)
//        			dki_dx=a[i]*(i+1)*(rr^i)*x/rr+ (dki1_dx*rr^i+ ki1*i* rr^(i-1)*x/r) + (dki2_dx+           ki2*(i-1)*rr^(i-2)*x/r)= // (WAS)
//        			dki_dx=a[i]*(i+1)*(rr^i)+      (dki1_dx*rr^i+ ki1*i* rr^(i-1)*x/r) + (dki2_dx*rr^(i-1) + ki2*(i-1)*rr^(i-2)*x/r)= // (CORRECTED)
//        			(a[i]*(i+1)*(rr^i)/rr)*x +                       // OK
//        			(axi*r^(i-1) + (axi*x+ayi*y)/r*i)) * x  +        // (WAS)
//        			(axi*r^i     + (axi*x+ayi*y)*i*rr^(i-2) * x) +   // (CORRECTED)
//        			2*(adi*y - aoi*x)*rr^(i-1)+                      // OK
//        			(aoi*(y^2-x^2)+2*adi*x*y)*(i-1)/r^3 * x          // (WAS)
//        			(aoi*(y^2-x^2)+2*adi*x*y)*(i-1)*r^(i-3) * x      // (CORRECTED)
//        			dki_dy=
//        			(a[i]*(i+1)*(rr^i)/rr)*y +                       // OK
//        			(ayi*r^(i-1)*y + (axi*x+ayi*y)/r*i)) * y +       // (WAS)
//        			(ayi*r^i       + (axi*x+ayi*y)*i*rr^(i-2) * x) + // (CORRECTED)
//        			2*(adi*x + aoi*y)*rr^(i-1)+                      // OK
//        			(aoi*(y^2-x^2)+2*adi*x*y)*(i-1)/r^3 * y          // (WAS)
//        			(aoi*(y^2-x^2)+2*adi*x*y)*(i-1)*r^(i-3) * x      // (CORRECTED)
// Below may speed up by branching for i==0, i==1 - starting with commonn code
        			double ai_iplus1_rr_pow_i_minus_1= a[i]*(i+1)*rr_pow_i_minus_1;
        			double rr_pow_i_minus_2=rr_pow_i_minus_1/rr;
//        			double ki1_div_r_mul_i=i*ki1/rr;
        			double ki1_r_pow_im2_mul_i=i*ki1*rr_pow_i_minus_2;       // (axi*x+ayi*y)*i*rr^(i-2)
//        			double ki2_div_r3_mul_im1=(i-1)*ki2/(rr*rr*rr);
        			double ki2_r_pow_im3_mul_im1=(i-1)*ki2*rr_pow_i_minus_2/rr; // (aoi*(y^2-x^2)+2*adi*x*y)* (i-1)*r^(i-3)
//        			double dki_dx= x*(ai_iplus1_rr_pow_i_minus_1 +
//        					r_xyod[i][0]*rr_pow_i_minus_1 + ki1_div_r_mul_i+ki2_div_r3_mul_im1)+
//        					2*(r_xyod[i][3]*y - r_xyod[i][2]*x)/rr_pow_i_minus_1;
        			double dki_dx= x*(ai_iplus1_rr_pow_i_minus_1 +
        					ki1_r_pow_im2_mul_i + ki2_r_pow_im3_mul_im1)+ // r_xyod[i][0]*rr_pow_i_minus_1 +
        					r_xyod[i][0]*rr_pow_i + 2*(r_xyod[i][3]*y - r_xyod[i][2]*x)*rr_pow_i_minus_1;
//        			double dki_dy= y*(ai_iplus1_rr_pow_i_minus_1 +
//        					r_xyod[i][1]*rr_pow_i_minus_1 + ki1_div_r_mul_i + ki2_div_r3_mul_im1)+
//        					2*(r_xyod[i][3]*x + r_xyod[i][2]*y)/rr_pow_i_minus_1;
        			double dki_dy= y*(ai_iplus1_rr_pow_i_minus_1 +
        					ki1_r_pow_im2_mul_i + ki2_r_pow_im3_mul_im1)+ // r_xyod[i][1]*rr_pow_i_minus_1 +
        					r_xyod[i][1]*rr_pow_i + 2*(r_xyod[i][3]*x + r_xyod[i][2]*y)*rr_pow_i_minus_1;
//        			double ai_iplus1_rr_pow_i= a[i]*(i+1)*rr_pow_i;
//        			double dki_dxmmc=          ai_iplus1_rr_pow_i*(x*(1.0-r_xyod[i][2])+y*r_xyod[i][3])/rr;
//        			double dki_dymmc=          ai_iplus1_rr_pow_i*(y*(1.0+r_xyod[i][2])+x*r_xyod[i][3])/rr;
                    // the following 4 shifts are "extra", on top of non-distorted (pinhole) - pinhole should be added to the per-term sum
        			double dDeltaXi_dx = ki+x*dki_dx; // mmc; // here dDelta*_d*mmc are relative (both I/O are fractions of distortionRadius)
        			double dDeltaXi_dy =    x*dki_dy; // mmc;
        			double dDeltaYi_dx =    y*dki_dx; // mmc;
        			double dDeltaYi_dy = ki+y*dki_dy; // mmc;
        			dDeltaX_dx += dDeltaXi_dx; // here dDelta*_d*mmc are relative (both I/O are fractions of distortionRadius)
        			dDeltaX_dy += dDeltaXi_dy;
        			dDeltaY_dx += dDeltaYi_dx;
        			dDeltaY_dy += dDeltaYi_dy;
//dpx_dai, dpi_dai - same as before
        			int index= distortionCIndex-i; // reverse order, // a8->9, a7->10, a6->11, a5->12, a->13, b->14, c->15
    				partDeriv[index][0]= rel_to_pix*x*(rr_pow_i_plus_1-1.0); // OK
    				partDeriv[index][1]=-rel_to_pix*y*(rr_pow_i_plus_1-1.0); // OK
// dpx_daxi,dpy_daxi,dpx_dayi,dpy_dayi
// dpx_daxi= x^2*r^i;
// dpy_daxi=-x*y*r^i;
// dpx_dayi= x*y*r^i;
// dpy_dayi=-y^2*r^i;
        			index=numRadialDerivatives-2+4*i;
        			double rel_to_pix_mul_rr_pow_i=rel_to_pix*rr_pow_i;
        			if (i>0){ // eccentricity is not applicabe to the first (C) term
        				partDeriv[index  ][0]=  rel_to_pix_mul_rr_pow_i*x2;
        				partDeriv[index  ][1]= -rel_to_pix_mul_rr_pow_i*xy;
        				partDeriv[index+1][0]=  rel_to_pix_mul_rr_pow_i*xy;
        				partDeriv[index+1][1]= -rel_to_pix_mul_rr_pow_i*y2;
        			}
// dpx_daoi,dpy_daoi,dpx_dadi,dpy_dadi
// dpx_daoi= x*(y^2-x^2)*r^(i-1);
// dpY_daoi=-y*(y^2-x^2)*r^(i-1);
// dpx_daDi= 2*x*xy*r^(i-1);
// dpy_daDi=-2*y*xy*r^(i-1);
    				double rel_to_pix_mul_y2_minus_x2_mul_rr_pow_i_minus_1=rel_to_pix* (y2-x2)*rr_pow_i_minus_1;
    				double rel_to_pix_mul_2xy_mul_rr_pow_i_minus_1=rel_to_pix*2.0*xy*rr_pow_i_minus_1;

        			index=numRadialDerivatives-2+4*i;
    				partDeriv[index+2][0]= x*rel_to_pix_mul_y2_minus_x2_mul_rr_pow_i_minus_1;
    				partDeriv[index+2][1]=-y*rel_to_pix_mul_y2_minus_x2_mul_rr_pow_i_minus_1;
    				partDeriv[index+3][0]= x*rel_to_pix_mul_2xy_mul_rr_pow_i_minus_1;
    				partDeriv[index+3][1]=-y*rel_to_pix_mul_2xy_mul_rr_pow_i_minus_1;
        		}
        	}

        	double [] xyDist={PXYmmc[0]+this.distortionRadius*deltaX,PXYmmc[1]+this.distortionRadius*deltaY};
        	// convert to sensor pixels coordinates
        	partDeriv[0][0]=  1000.0/this.pixelSize*xyDist[0] + this.px0;
        	partDeriv[0][1]= -1000.0/this.pixelSize*xyDist[1] + this.py0;

        	if (!calculateAll) { // TODO: how to deal with it when calculating Jacobian???
        		// TODO: Looking away from the target, trying only with no dervatives. Do the same for derivatives too?
            	if ((XeYeZe[2]<0.0) || (r2>maxRelativeRadius*maxRelativeRadius)){ // non-distorted too far from the axis
            		partDeriv[0][0]=Double.NaN;
            		partDeriv[0][1]=Double.NaN;
            	}
        		return partDeriv;
        	}
        	// correct parameter derivatives to cumulative version
        	if (this.cummulativeCorrection){
            	for (int i=r_xyod.length-2;i>=0;i--){
        			int index=numRadialDerivatives-2+4*i;
        			if (i>0){ // eccentricity is not applicabe to the first (C) term
        				partDeriv[index  ][0]+=partDeriv[index+4][0]; // oob=36
        				partDeriv[index  ][1]+=partDeriv[index+4][1];
        				partDeriv[index+1][0]+=partDeriv[index+5][0];
        				partDeriv[index+1][1]+=partDeriv[index+5][1];
        			}
    				partDeriv[index+2][0]+=partDeriv[index+6][0];
    				partDeriv[index+2][1]+=partDeriv[index+6][1];
    				partDeriv[index+3][0]+=partDeriv[index+7][0];
    				partDeriv[index+3][1]+=partDeriv[index+7][1];
            	}
        	}

// convert    dDelta*_d* from relative/relative to pix/mm (invert pixel Y direction)
// added 1.0 to account for non-distorted (pinhole) shift
        	double dPx_dPinholeX= (1.0+dDeltaX_dx)*1000.0/this.pixelSize;
        	double dPx_dPinholeY=       dDeltaX_dy*1000.0/this.pixelSize;
        	double dPy_dPinholeX=      -dDeltaY_dx*1000.0/this.pixelSize;
        	double dPy_dPinholeY=-(1.0+dDeltaY_dy)*1000.0/this.pixelSize;
        	if (this.debugLevel>2){
        		System.out.println(" deltaX="+deltaX+" deltaY="+deltaY+" xyDist[0]="+xyDist[0]+" xyDist[1]="+xyDist[1]);
        		System.out.println(" PXYmmc[0]="+PXYmmc[0]+" PXYmmc[1]="+PXYmmc[1]+" x="+x+" y"+y);
        		System.out.println(" dDeltaX_dxmmc="+dDeltaX_dx+" dDeltaX_dymmc="+dDeltaX_dy+" dDeltaY_dxmmc="+dDeltaY_dx+" dDeltaY_dymmc"+dDeltaY_dy);
        		System.out.println(" dPx_dPinholeX="+dPx_dPinholeX+" dPx_dPinholeY="+dPx_dPinholeY+" dPy_dPinholeX="+dPy_dPinholeX+" dPy_dPinholeY"+dPy_dPinholeY);
        	}

            double K=Math.PI/180; // multiply all derivatives my angles
// dPX/dphi   =  1000/Psz* dxDist/dphi
        	partDeriv[ 1][0]=  K*(dPx_dPinholeX*dPXYmmc[1][0]+dPx_dPinholeY*dPXYmmc[1][1]);
        	partDeriv[ 1][1]=  K*(dPy_dPinholeX*dPXYmmc[1][0]+dPy_dPinholeY*dPXYmmc[1][1]);

// dPX/dtheta =  1000/Psz* dxDist/dtheta
        	partDeriv[ 2][0]=  K*(dPx_dPinholeX*dPXYmmc[2][0]+dPx_dPinholeY*dPXYmmc[2][1]);
        	partDeriv[ 2][1]=  K*(dPy_dPinholeX*dPXYmmc[2][0]+dPy_dPinholeY*dPXYmmc[2][1]);
// dPX/dpsi   =  1000/Psz* dxDist/dpsi
        	partDeriv[ 3][0]=  K*(dPx_dPinholeX*dPXYmmc[3][0]+dPx_dPinholeY*dPXYmmc[3][1]);
        	partDeriv[ 3][1]=  K*(dPy_dPinholeX*dPXYmmc[3][0]+dPy_dPinholeY*dPXYmmc[3][1]);
// dPX/dX0    =  1000/Psz* dxDist/dX0
        	partDeriv[ 4][0]=  dPx_dPinholeX*dPXYmmc[4][0]+dPx_dPinholeY*dPXYmmc[4][1];
        	partDeriv[ 4][1]=  dPy_dPinholeX*dPXYmmc[4][0]+dPy_dPinholeY*dPXYmmc[4][1];
// dPX/dY0    =  1000/Psz* dxDist/dY0
        	partDeriv[ 5][0]=  dPx_dPinholeX*dPXYmmc[5][0]+dPx_dPinholeY*dPXYmmc[5][1];
        	partDeriv[ 5][1]=  dPy_dPinholeX*dPXYmmc[5][0]+dPy_dPinholeY*dPXYmmc[5][1];
// dPX/dZ0    =  1000/Psz* dxDist/dZ0
        	partDeriv[ 6][0]=  dPx_dPinholeX*dPXYmmc[6][0]+dPx_dPinholeY*dPXYmmc[6][1];
        	partDeriv[ 6][1]=  dPy_dPinholeX*dPXYmmc[6][0]+dPy_dPinholeY*dPXYmmc[6][1];
// dPX/df     =  1000/Psz* dxDist/df
        	partDeriv[ 7][0]=  dPx_dPinholeX*dPXYmmc[7][0]+dPx_dPinholeY*dPXYmmc[7][1];
        	partDeriv[ 7][1]=  dPy_dPinholeX*dPXYmmc[7][0]+dPy_dPinholeY*dPXYmmc[7][1];
// dPX/ddist  =  1000/Psz* dxDist/ddist
        	partDeriv[ 8][0]=  dPx_dPinholeX*dPXYmmc[8][0]+dPx_dPinholeY*dPXYmmc[8][1];
        	partDeriv[ 8][1]=  dPy_dPinholeX*dPXYmmc[8][0]+dPy_dPinholeY*dPXYmmc[8][1];

// dPX/dPX0   =  1
// dPY/dPX0   =  0
        	partDeriv[16][0]=  1.0;
        	partDeriv[16][1]=  0.0;
// dPX/dPY0   =  0
// dPY/dPY0   =  1
        	partDeriv[17][0]=  0.0;
        	partDeriv[17][1]=  1.0;
            return partDeriv;
        }



        public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"focalLength",this.focalLength+"");
			properties.setProperty(prefix+"pixelSize",this.pixelSize+"");
			properties.setProperty(prefix+"distortionRadius",this.distortionRadius+"");
			properties.setProperty(prefix+"distortionA8",this.distortionA8+"");
			properties.setProperty(prefix+"distortionA7",this.distortionA7+"");
			properties.setProperty(prefix+"distortionA6",this.distortionA6+"");
			properties.setProperty(prefix+"distortionA5",this.distortionA5+"");
			properties.setProperty(prefix+"distortionA",this.distortionA+"");
			properties.setProperty(prefix+"distortionB",this.distortionB+"");
			properties.setProperty(prefix+"distortionC",this.distortionC+"");
			properties.setProperty(prefix+"yaw",this.yaw+"");
			properties.setProperty(prefix+"pitch",this.pitch+"");
			properties.setProperty(prefix+"roll",this.roll+"");
			properties.setProperty(prefix+"x0",this.x0+"");
			properties.setProperty(prefix+"y0",this.y0+"");
			properties.setProperty(prefix+"z0",this.z0+"");
			properties.setProperty(prefix+"distance",this.distance+"");
			properties.setProperty(prefix+"px0",this.px0+"");
			properties.setProperty(prefix+"py0",this.py0+"");
			properties.setProperty(prefix+"flipVertical",this.flipVertical+"");
			for (int i=0;i<this.r_xy.length;i++){
				properties.setProperty(prefix+"r_xy_"+i+"_x",this.r_xy[i][0]+"");
				properties.setProperty(prefix+"r_xy_"+i+"_y",this.r_xy[i][1]+"");
			}
			for (int i=0;i<this.r_od.length;i++){
				properties.setProperty(prefix+"r_od_"+i+"_o",this.r_od[i][0]+"");
				properties.setProperty(prefix+"r_od_"+i+"_d",this.r_od[i][1]+"");
			}
		}
    	public void setDefaultNonRadial(){
			r_od=new double [r_od_dflt.length][2];
			for (int i=0;i<r_od.length;i++) r_od[i]=r_od_dflt[i].clone();
			r_xy=new double [r_xy_dflt.length][2];
			for (int i=0;i<r_xy.length;i++) r_xy[i]=r_xy_dflt[i].clone();
    	}

    	public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"focalLength")!=null)
				this.focalLength=Double.parseDouble(properties.getProperty(prefix+"focalLength"));
			if (properties.getProperty(prefix+"pixelSize")!=null)
				this.pixelSize=Double.parseDouble(properties.getProperty(prefix+"pixelSize"));
			if (properties.getProperty(prefix+"distortionRadius")!=null)
				this.distortionRadius=Double.parseDouble(properties.getProperty(prefix+"distortionRadius"));
			if (properties.getProperty(prefix+"distortionA8")!=null)
				this.distortionA8=Double.parseDouble(properties.getProperty(prefix+"distortionA8"));
			if (properties.getProperty(prefix+"distortionA7")!=null)
				this.distortionA7=Double.parseDouble(properties.getProperty(prefix+"distortionA7"));
			if (properties.getProperty(prefix+"distortionA6")!=null)
				this.distortionA6=Double.parseDouble(properties.getProperty(prefix+"distortionA6"));
			if (properties.getProperty(prefix+"distortionA5")!=null)
				this.distortionA5=Double.parseDouble(properties.getProperty(prefix+"distortionA5"));
			if (properties.getProperty(prefix+"distortionA")!=null)
				this.distortionA=Double.parseDouble(properties.getProperty(prefix+"distortionA"));
			if (properties.getProperty(prefix+"distortionB")!=null)
				this.distortionB=Double.parseDouble(properties.getProperty(prefix+"distortionB"));
			if (properties.getProperty(prefix+"distortionC")!=null)
				this.distortionC=Double.parseDouble(properties.getProperty(prefix+"distortionC"));
			if (properties.getProperty(prefix+"yaw")!=null)
				this.yaw=Double.parseDouble(properties.getProperty(prefix+"yaw"));
			if (properties.getProperty(prefix+"pitch")!=null)
				this.pitch=Double.parseDouble(properties.getProperty(prefix+"pitch"));
			if (properties.getProperty(prefix+"roll")!=null)
				this.roll=Double.parseDouble(properties.getProperty(prefix+"roll"));
			if (properties.getProperty(prefix+"x0")!=null)
				this.x0=Double.parseDouble(properties.getProperty(prefix+"x0"));
			if (properties.getProperty(prefix+"y0")!=null)
				this.y0=Double.parseDouble(properties.getProperty(prefix+"y0"));
			if (properties.getProperty(prefix+"z0")!=null)
				this.z0=Double.parseDouble(properties.getProperty(prefix+"z0"));
			if (properties.getProperty(prefix+"distance")!=null)
				this.distance=Double.parseDouble(properties.getProperty(prefix+"distance"));
			if (properties.getProperty(prefix+"px0")!=null)
				this.px0=Double.parseDouble(properties.getProperty(prefix+"px0"));
			if (properties.getProperty(prefix+"py0")!=null)
				this.py0=Double.parseDouble(properties.getProperty(prefix+"py0"));
			if (properties.getProperty(prefix+"flipVertical")!=null)
				this.flipVertical=Boolean.parseBoolean(properties.getProperty(prefix+"flipVertical"));

			setDefaultNonRadial();
			for (int i=0;i<this.r_xy.length;i++){
				if (properties.getProperty(prefix+"r_xy_"+i+"_x")!=null) this.r_xy[i][0]=Double.parseDouble(properties.getProperty(prefix+"r_xy_"+i+"_x"));
				if (properties.getProperty(prefix+"r_xy_"+i+"_y")!=null) this.r_xy[i][1]=Double.parseDouble(properties.getProperty(prefix+"r_xy_"+i+"_y"));
			}
			for (int i=0;i<this.r_od.length;i++){
				if (properties.getProperty(prefix+"r_od_"+i+"_o")!=null) this.r_od[i][0]=Double.parseDouble(properties.getProperty(prefix+"r_od_"+i+"_o"));
				if (properties.getProperty(prefix+"r_od_"+i+"_d")!=null) this.r_od[i][1]=Double.parseDouble(properties.getProperty(prefix+"r_od_"+i+"_d"));
			}
			recalcCommons();
		}
		public boolean showDialog() {
			GenericDialog gd = new GenericDialog("Lens distortion, location and orientation");
			gd.addNumericField("Lens focal length",               this.focalLength, 3,6,"mm");
			gd.addNumericField("Sensor pixel period",             this.pixelSize, 3,6,"um");
			gd.addNumericField("Distortion radius (halw width)",  this.distortionRadius, 5,8,"mm");
			gd.addNumericField("Distortion A8(r^5)",              this.distortionA8, 6,8,"");
			gd.addNumericField("Distortion A7(r^5)",              this.distortionA7, 6,8,"");
			gd.addNumericField("Distortion A6(r^5)",              this.distortionA6, 6,8,"");
			gd.addNumericField("Distortion A5(r^5)",              this.distortionA5, 6,8,"");
			gd.addNumericField("Distortion A (r^4)",              this.distortionA, 6,8,"");
			gd.addNumericField("Distortion B (r^3)",              this.distortionB, 6,8,"");
			gd.addNumericField("Distortion C (r^2)",              this.distortionC, 6,8,"");
			gd.addNumericField("Lens axis from perpendicular to the pattern, positive - clockwise (from top)", this.yaw, 2,6,"degrees");
			gd.addNumericField("Lens axis from perpendicular to the pattern, positive - up",                   this.pitch, 2,6,"degrees");
			gd.addNumericField("Rotation around lens axis, positive - clockwise (looking to pattern)",         this.roll, 2,6,"degrees");
			gd.addNumericField("Lens axis from the pattern center, (to the right)",                            this.x0, 1,6,"mm");
			gd.addNumericField("Lens axis from the pattern center, (down)",                                    this.y0, 1,6,"mm");
			gd.addNumericField("Lens axis from the pattern center, (away from camera, normally 0.0)",          this.z0, 1,6,"mm");
			gd.addNumericField("Distance from the lens input pupil to the pattern plane along the camera axis",this.distance, 1,6,"mm");
			gd.addNumericField("Lens axis on the sensor (horizontal, from left edge)",                         this.px0, 1,6,"pixels");
			gd.addNumericField("Lens axis on the sensor (vertical, from top edge)",                            this.py0, 1,6,"pixels");
			gd.addCheckbox    ("Camera looks through the mirror",                                              this.flipVertical);
			gd.addMessage("=== non-radial model parameters ===");
			gd.addMessage("For r^2 (Distortion C):");
			gd.addNumericField("Orthogonal elongation for r^2",   100*this.r_od[0][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^2",     100*this.r_od[0][1], 3,7,"%");
			gd.addMessage("For r^3 (Distortion B):");
			gd.addNumericField("Distortion center shift X for r^3", 100*this.r_xy[0][0], 1,6,"%");
			gd.addNumericField("Distortion center shift Y for r^3", 100*this.r_xy[0][1], 1,6,"%");
			gd.addNumericField("Orthogonal elongation for r^3",  100*this.r_od[1][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^3",    100*this.r_od[1][1], 3,7,"%");
			gd.addMessage("For r^4 (Distortion A):");
			gd.addNumericField("Distortion center shift X for r^4", 100*this.r_xy[1][0], 1,6,"%");
			gd.addNumericField("Distortion center shift Y for r^4", 100*this.r_xy[1][1], 1,6,"%");
			gd.addNumericField("Orthogonal elongation for r^4",  100*this.r_od[2][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^4",    100*this.r_od[2][1], 3,7,"%");
			gd.addMessage("For r^5 (Distortion A5):");
			gd.addNumericField("Distortion center shift X for r^5", 100*this.r_xy[2][0], 1,6,"%");
			gd.addNumericField("Distortion center shift Y for r^5", 100*this.r_xy[2][1], 1,6,"%");
			gd.addNumericField("Orthogonal elongation for r^5",  100*this.r_od[3][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^5",    100*this.r_od[3][1], 3,7,"%");
			gd.addMessage("For r^6 (Distortion A6:");
			gd.addNumericField("Distortion center shift X for r^6", 100*this.r_xy[3][0], 1,6,"%");
			gd.addNumericField("Distortion center shift Y for r^6", 100*this.r_xy[3][1], 1,6,"%");
			gd.addNumericField("Orthogonal elongation for r^6",  100*this.r_od[4][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^6",    100*this.r_od[4][1], 3,7,"%");
			gd.addMessage("For r^7 (Distortion A7):");
			gd.addNumericField("Distortion center shift X for r^7", 100*this.r_xy[4][0], 1,6,"%");
			gd.addNumericField("Distortion center shift Y for r^7", 100*this.r_xy[4][1], 1,6,"%");
			gd.addNumericField("Orthogonal elongation for r^7",  100*this.r_od[5][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^7",    100*this.r_od[5][1], 3,7,"%");
			gd.addMessage("For r^8 (Distortion A8):");
			gd.addNumericField("Distortion center shift X for r^8", 100*this.r_xy[5][0], 1,6,"%");
			gd.addNumericField("Distortion center shift Y for r^8", 100*this.r_xy[5][1], 1,6,"%");
			gd.addNumericField("Orthogonal elongation for r^8",  100*this.r_od[6][0], 3,7,"%");
			gd.addNumericField("Diagonal elongation for r^8",    100*this.r_od[6][1], 3,7,"%");
    		WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.focalLength=     gd.getNextNumber();
			this.pixelSize=       gd.getNextNumber();
			this.distortionRadius=gd.getNextNumber();
			this.distortionA8=    gd.getNextNumber();
			this.distortionA7=    gd.getNextNumber();
			this.distortionA6=    gd.getNextNumber();
			this.distortionA5=    gd.getNextNumber();
			this.distortionA=     gd.getNextNumber();
			this.distortionB=     gd.getNextNumber();
			this.distortionC=     gd.getNextNumber();
			this.yaw=             gd.getNextNumber();
			this.pitch=           gd.getNextNumber();
			this.roll=            gd.getNextNumber();
			this.x0=              gd.getNextNumber();
			this.y0=              gd.getNextNumber();
			this.z0=              gd.getNextNumber();
			this.distance=        gd.getNextNumber();
			this.px0=             gd.getNextNumber();
			this.py0=             gd.getNextNumber();
			this.flipVertical=    gd.getNextBoolean();

			this.r_od[0][0]= 0.01*gd.getNextNumber();
			this.r_od[0][1]= 0.01*gd.getNextNumber();
			this.r_xy[0][0]= 0.01*gd.getNextNumber();
			this.r_xy[0][1]= 0.01*gd.getNextNumber();
			this.r_od[1][0]= 0.01*gd.getNextNumber();
			this.r_od[1][1]= 0.01*gd.getNextNumber();
			this.r_xy[1][0]= 0.01*gd.getNextNumber();
			this.r_xy[1][1]= 0.01*gd.getNextNumber();
			this.r_od[2][0]= 0.01*gd.getNextNumber();
			this.r_od[2][1]= 0.01*gd.getNextNumber();
			this.r_xy[2][0]= 0.01*gd.getNextNumber();
			this.r_xy[2][1]= 0.01*gd.getNextNumber();
			this.r_od[3][0]= 0.01*gd.getNextNumber();
			this.r_od[3][1]= 0.01*gd.getNextNumber();
			this.r_xy[3][0]= 0.01*gd.getNextNumber();
			this.r_xy[3][1]= 0.01*gd.getNextNumber();
			this.r_od[4][0]= 0.01*gd.getNextNumber();
			this.r_od[4][1]= 0.01*gd.getNextNumber();
			this.r_xy[4][0]= 0.01*gd.getNextNumber();
			this.r_xy[4][1]= 0.01*gd.getNextNumber();
			this.r_od[5][0]= 0.01*gd.getNextNumber();
			this.r_od[5][1]= 0.01*gd.getNextNumber();
			this.r_xy[5][0]= 0.01*gd.getNextNumber();
			this.r_xy[5][1]= 0.01*gd.getNextNumber();
			this.r_od[6][0]= 0.01*gd.getNextNumber();
			this.r_od[6][1]= 0.01*gd.getNextNumber();
			return true;
		}
	    /**
	     * Calculate/set  this.lensDistortionParameters and this.interParameterDerivatives
     * UPDATE - Modifies lensDistortionParameters, not "this" for multi-threaded
	     * @param parVect 21-element vector for eyesis sub-camera, including common and individual parameters
	     * @param mask -mask - which partial derivatives are needed to be calculated (others will be null)
	     * @param calculateDerivatives calculate array of partial derivatives, if false - just the values
	     */
	    public void lensCalcInterParamers(
	    		LensDistortionParameters lensDistortionParameters,
	    		boolean isTripod,
	    		boolean cartesian,
	    		double  pixelSize,
	    		double  distortionRadius,
	            double [][] interParameterDerivatives, //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21) if null - just values, no derivatives
	    		double [] parVect,
	    		boolean [] mask // calculate only selected derivatives (all parVect values are still
	    		){
	    	boolean calculateDerivatives=(interParameterDerivatives!=null);  // calculate this.interParameterDerivatives -derivatives array (false - just this.values)
	    	 // change meaning of goniometerHorizontal (tripod vertical) and goniometerAxial (tripod horizontal)

	    	// Alternative variables for cartesian/cylindrical modes
	    	double azimuth_cyl=  cartesian? Double.NaN:  parVect[0];
	    	double right_cart=   cartesian? parVect[0] : Double.NaN;
	    	double radius_cyl=   cartesian? Double.NaN:  parVect[1];
	    	double forward_cart= cartesian? parVect[1]:  Double.NaN;
	    	double height=       parVect[2];
	    	double phi_cyl=      cartesian? Double.NaN:  parVect[3];
	    	double heading_cart= cartesian? parVect[3]:  Double.NaN;
	    	double theta=  parVect[4];
	    	double psi=    parVect[5];
	    	double goniometerHorizontal=parVect[6];
	    	double goniometerAxial=parVect[7];
	    	double interAxisDistance=parVect[8];
	    	double interAxisAngle=parVect[9];
	    	double horAxisErrPhi=parVect[10];
	    	double horAxisErrPsi=parVect[11];
	    	double entrancePupilForward=parVect[12];
	    	double centerAboveHorizontal=parVect[13];
	    	double GXYZ0=parVect[14];
	    	double GXYZ1=parVect[15];
	    	double GXYZ2=parVect[16];

	    	double cPS=   Math.cos(psi*Math.PI/180); // subCam.psi
	    	double sPS=   Math.sin(psi*Math.PI/180); // subCam.psi
	    	double cTH=   Math.cos(theta*Math.PI/180); // subCam.theta
	    	double sTH=   Math.sin(theta*Math.PI/180); // subCam.theta

	    	double cAZP=  Math.cos((cartesian?heading_cart:(azimuth_cyl+phi_cyl))*Math.PI/180); //subCam.azimuth+subCam.phi
	    	double sAZP=  Math.sin((cartesian?heading_cart:(azimuth_cyl+phi_cyl))*Math.PI/180); //subCam.azimuth+subCam.phi

	    	// renaming the following 2 to be replaced for cartesian coordinates
	    	double cAZ_cyl=   Math.cos(azimuth_cyl*Math.PI/180); //subCam.azimuth
	    	double sAZ_cyl=   Math.sin(azimuth_cyl*Math.PI/180); //subCam.azimuth

	    	double cGA=   Math.cos(goniometerAxial*Math.PI/180); //eyesisCameraParameters.goniometerAxial
	    	double sGA=   Math.sin(goniometerAxial*Math.PI/180); //eyesisCameraParameters.goniometerAxial
	    	double cGH=   Math.cos(goniometerHorizontal*Math.PI/180); //eyesisCameraParameters.goniometerHorizontal
	    	double sGH=   Math.sin(goniometerHorizontal*Math.PI/180); //eyesisCameraParameters.goniometerHorizontal
	    	double cGIAA= Math.cos(interAxisAngle*Math.PI/180); // eyesisCameraParameters.interAxisAngle
	    	double sGIAA= Math.sin(interAxisAngle*Math.PI/180); //eyesisCameraParameters.interAxisAngle
	    	double cHAEPH=Math.cos(horAxisErrPhi*Math.PI/180); //eyesisCameraParameters.horAxisErrPhi
	    	double sHAEPH=Math.sin(horAxisErrPhi*Math.PI/180); //eyesisCameraParameters.horAxisErrPhi
	    	double cHAEPS=Math.cos(horAxisErrPsi*Math.PI/180); //eyesisCameraParameters.horAxisErrPsi
	    	double sHAEPS=Math.sin(horAxisErrPsi*Math.PI/180); //eyesisCameraParameters.horAxisErrPsi\

	/*
	 0) Translate by distance to entrance pupil (lens center)
	    	| Xc0 |   | 0                     |   |Xc|
	    	| Yc0 | = | 0                     | + |Yc|
	    	| Zc0 |   | entrancePupilForward  |   |Zc|
	    	*/
	    	double [][] aT0={{0.0},{0.0},{entrancePupilForward}};
	    	Matrix T0=new Matrix(aT0);

	/*
	Converting from the sub-camera coordinates to the target coordinates
	1) rotate by -psi around CZ: Vc1= R1*Vc
	| Xc1 |   | cos(psi)  sin(psi)    0         |   |Xc0|
	| Yc1 | = |-sin(psi)  cos(psi)    0         | * |Yc0|
	| Zc1 |   |    0         0        1         |   |Zc0|
	*/
	    	double [][] aR1={{cPS,sPS,0.0},{-sPS,cPS,0.0},{0.0,0.0,1.0}};
	    	Matrix R1=new Matrix(aR1);
	/*
	2) rotate by - theta around C1X:Vc2= R2*Vc1
	| Xc2 |   |    1         0         0        |   |Xc1|
	| Yc2 | = |    0    cos(theta)   sin(theta) | * |Yc1|
	| Zc2 |   |    0   -sin(theta)   cos(theta) |   |Zc1|
	*/
	    	double [][] aR2={{1.0,0.0,0.0},{0.0,cTH,sTH},{0.0,-sTH,cTH}};
	    	Matrix R2=new Matrix(aR2);
	/*
	3) rotate by -(azimuth+phi) around C2Y:Vc3= R3*Vc2
	| Xc3 |   | cos(azimuth+phi)    0   sin(azimuth+phi)   |   |Xc2|
	| Yc3 | = |     0               1         0            | * |Yc2|
	| Zc3 |   | -sin(azimuth+phi)   0   cos(azimuth+phi)   |   |Zc2|
	*/
	    	double [][] aR3={{cAZP,0.0,sAZP},{0.0,1.0,0.0},{-sAZP,0.0,cAZP}};
	    	Matrix R3=new Matrix(aR3);
	/*
	4) Now axes are aligned, just translate to get to eyesis coordinates: Vey= T1+Vc3
	| Xey |   |      r * sin (azimuth)       |   |Xc3|
	| Yey | = | height+centerAboveHorizontal | + |Yc3|
	| Zey |   |      r * cos (azimuth)       |   |Zc3|
	*/
	    	double [][] aT1_cyl= {{radius_cyl*sAZ_cyl},{(height+centerAboveHorizontal)},{radius_cyl*cAZ_cyl}}; // {{subCam.radius*sAZ},{subCam.height},{subCam.radius*cAZ}};
	    	double [][] aT1_cart={{right_cart},        {(height+centerAboveHorizontal)},{forward_cart}}; // {{subCam.radius*sAZ},{subCam.height},{subCam.radius*cAZ}};

	    	double [][] aT1= cartesian ? aT1_cart : aT1_cyl;
	    	Matrix T1=new Matrix(aT1);

	/**
	5) rotate around moving goniometer axis, by (-goniometerAxial) - same as around EY: Vgm1=R4*Vey
	| Xgm1 |   |    1           0                     0              |   |Xey|
	| Ygm1 | = |    0   cos(goniometerAxial)  sin(goniometerAxial)   | * |Yey|
	| Zgm1 |   |    0   -sin(goniometerAxial) cos(goniometerAxial)   |   |Zey|
	*/
	    	double [][] aR4_tripod=    {{1.0,0.0,0.0},{0.0,cGA,sGA},{0.0,-sGA,cGA}};
	/*
	5) rotate around moving goniometer axis, by (-goniometerAxial) - same as around EY: Vgm1=R4*Vey
	| Xgm1 |   | cos(goniometerAxial)    0   sin(goniometerAxial)   |   |Xey|
	| Ygm1 | = |        0                1             0            | * |Yey|
	| Zgm1 |   |-sin(goniometerAxial)    0   cos(goniometerAxial)   |   |Zey|
	*/
			double [][] aR4_goniometer={{cGA,0.0,sGA},{0.0,1.0,0.0},{-sGA,0.0,cGA}};
			Matrix R4=new Matrix(isTripod?aR4_tripod:aR4_goniometer);
	/*
	6) move to the goniometer horizontal axis:Vgm2=T2+Vgm1
	| Xgm2 |   |                 0  |   |Xgm1|
	| Ygm2 | = |                 0  | + |Ygm1|
	| Zgm2 |   | interAxisDistance  |   |Zgm1|
	*/
	    	double [][] aT2={{0.0},{0.0},{interAxisDistance}}; //eyesisCameraParameters.interAxisDistance
	    	Matrix T2=new Matrix(aT2);
	/*
	7) rotate around Zgm2 by -interAxisAngle, so Xgm3 axis is the same as horizontal goniometer axis: Vgm3=R5*Vgm2
	| Xgm3 |   | cos(interAxisAngle)  sin(interAxisAngle)    0         |   |Xgm2|
	| Ygm3 | = |-sin(interAxisAngle)  cos(interAxisAngle)    0         | * |Ygm2|
	| Zgm3 |   |    0                         0              1         |   |Zgm2|
	*/
	    	double [][] aR5={{cGIAA,sGIAA,0.0},{-sGIAA,cGIAA,0.0},{0.0,0.0,1.0}};
	    	Matrix R5=new Matrix(aR5);
	/**
	8)  rotate around goniometer horizontal axis (Xgm3) by -goniometerHorizontal: Vgm4=R6*Vgm3
	| Xgm4 |   |    cos(goniometerHorizontal)  0   sin(goniometerHorizontal) |   |Xgm3|
	| Ygm4 | = |               0               1               0             | * |Ygm3|
	| Zgm4 |   |   -sin(goniometerHorizontal)  0   cos(goniometerHorizontal) |   |Zgm3|
	*/
	    	double [][] aR6_tripod=    {{cGH,0.0,sGH},{0.0,1.0,0.0},{-sGH,0.0,cGH}};
	/*
	8)  rotate around goniometer horizontal axis (Xgm3) by -goniometerHorizontal: Vgm4=R6*Vgm3
	| Xgm4 |   |    1                 0                           0            |   |Xgm3|
	| Ygm4 | = |    0    cos(goniometerHorizontal)   sin(goniometerHorizontal) | * |Ygm3|
	| Zgm4 |   |    0   -sin(goniometerHorizontal)   cos(goniometerHorizontal) |   |Zgm3|
	*/
	    	double [][] aR6_goniometer={{1.0,0.0,0.0},{0.0,cGH,sGH},{0.0,-sGH,cGH}};
			Matrix R6=new Matrix(isTripod?aR6_tripod:aR6_goniometer);
	/*
	9) Correct roll error of the goniometer horizontal axis - rotate by -horAxisErrPsi around Zgm4: Vgm5=R7*Vgm4
	| Xgm5 |   | cos(horAxisErrPsi)  sin(horAxisErrPsi)    0         |   |Xgm4|
	| Ygm5 | = |-sin(horAxisErrPsi)  cos(horAxisErrPsi)    0         | * |Ygm4|
	| Zgm5 |   |         0                   0             1         |   |Zgm4|
	*/
	    	double [][] aR7={{cHAEPS,sHAEPS,0.0},{-sHAEPS,cHAEPS,0.0},{0.0,0.0,1.0}};
	    	Matrix R7=new Matrix(aR7);
	/*
	10) Correct azimuth error of the goniometer hoirizontal axis - rotate by -horAxisErrPhi around Ygm5: Vgm6=R8*Vgm5

	| Xgm6 |   | cos(horAxisErrPhi)      0   sin(horAxisErrPhi)   |   |Xgm5|
	| Ygm6 | = |        0                1             0          | * |Ygm5|
	| Zgm6 |   |-sin(horAxisErrPhi)      0   cos(horAxisErrPhi)   |   |Zgm5|
	For Tripod - rotate around X-axis (like theta)
	| Xgm6 |   |   1           0                       0      |   |Xgm5|
	| Ygm6 | = |   0  cos(horAxisErrPhi)  sin(horAxisErrPhi)  | * |Ygm5|
	| Zgm6 |   |   0 -sin(horAxisErrPhi)  cos(horAxisErrPhi)  |   |Zgm5|

	*/
	    	double [][] aR8_tripod=    {{1.0,   0.0,   0.0},{0.0,cHAEPH,sHAEPH},{0.0,  -sHAEPH,cHAEPH}};
	    	double [][] aR8_goniometer={{cHAEPH,0.0,sHAEPH},{0.0,   1.0,   0.0},{-sHAEPH,  0.0,cHAEPH}};
	    	Matrix R8=new Matrix(isTripod?aR8_tripod:aR8_goniometer);
	/*
	11) translate to the target zero point: Vt=  T3+Vgm6
	| Xt |   | GXYZ[0]  |   |Xgm6|
	| Yt | = |-GXYZ[1]  | + |Ygm6| // Y - up positive
	| Zt |   |-GXYZ[2]  |   |Zgm6| // Z - away positive
	    	 */
//	    	double [][] aT3={{parVect[12]},{-parVect[13]},{-parVect[14]}};//{{eyesisCameraParameters.GXYZ[0]},{eyesisCameraParameters.GXYZ[1]},{eyesisCameraParameters.GXYZ[2]}};

	    	double [][] aT3={{GXYZ0},{-GXYZ1},{-GXYZ2}};//{{eyesisCameraParameters.GXYZ[0]},{eyesisCameraParameters.GXYZ[1]},{eyesisCameraParameters.GXYZ[2]}};
//	    	double [][] aT3={{parVect[12]},{ parVect[13]},{-parVect[14]}};//{{eyesisCameraParameters.GXYZ[0]},{eyesisCameraParameters.GXYZ[1]},{eyesisCameraParameters.GXYZ[2]}};
	    	Matrix T3=new Matrix(aT3);

	// MA=R8*R7*R6*R5*R4*R3*R2*R1;
	// MB=T3+(R8*R7*R6*R5*(T2+R4*T1)); - old
	// MB=T3+(R8*R7*R6*R5*(T2+R4*(T1+R3*R2*R1*T0)));
	        Matrix MA=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
//	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1)))))));
	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1.plus(R3.times(R2.times(R1.times(T0)))) )))))));
			if (this.debugLevel>2) {
				System.out.println("MA:");
				MA.print(10, 5);
				System.out.println("MB:");
				MB.print(10, 5);
				System.out.println("T3:");
				T3.print(10, 5);
//				System.out.println("interParameterDerivatives[0]="+sprintfArray(interParameterDerivatives[0]));
			}
	        double [] extrinsicParams=parametersFromMAMB(MA,MB); // all after 6 are 0;

	        lensDistortionParameters.pixelSize =        pixelSize;
	        lensDistortionParameters.distortionRadius = distortionRadius;

	        lensDistortionParameters.distance=extrinsicParams[2];
	        lensDistortionParameters.x0=extrinsicParams[0];
	        lensDistortionParameters.y0=extrinsicParams[1];
	        lensDistortionParameters.z0=0.0; // used here
	        lensDistortionParameters.pitch=extrinsicParams[4];
	        lensDistortionParameters.yaw=  extrinsicParams[3];
	        lensDistortionParameters.roll= extrinsicParams[5];

			lensDistortionParameters.focalLength=parVect[17]; //subCam.focalLength;
			lensDistortionParameters.px0=parVect[18]; //subCam.px0;
			lensDistortionParameters.py0=parVect[19]; //subCam.py0;
			lensDistortionParameters.distortionA8=parVect[20]; //subCam.distortion5;
			lensDistortionParameters.distortionA7=parVect[21]; //subCam.distortion5;
			lensDistortionParameters.distortionA6=parVect[22]; //subCam.distortion5;
			lensDistortionParameters.distortionA5=parVect[23]; //subCam.distortion5;
			lensDistortionParameters.distortionA=parVect[24]; //subCam.distortionA;
			lensDistortionParameters.distortionB=parVect[25]; //subCam.distortionB;
			lensDistortionParameters.distortionC=parVect[26]; //subCam.distortionC;

			lensDistortionParameters.r_xy=new double [6][2];
			lensDistortionParameters.r_od=new double [7][2];

			int index=27;
			for (int i=0;i<lensDistortionParameters.r_od.length;i++){
				if (i>0){
					lensDistortionParameters.r_xy[i-1][0]=parVect[index++];
					lensDistortionParameters.r_xy[i-1][1]=parVect[index++];
				}
				lensDistortionParameters.r_od[i][0]=parVect[index++];
				lensDistortionParameters.r_od[i][1]=parVect[index++];
			}

			lensDistortionParameters.recalcCommons(); // uncluding lensDistortionParameters.r_xyod


	        if (this.debugLevel>2){
	        	System.out.println("lensDistortionParameters.recalcCommons()");
	        	System.out.println("lensDistortionParameters.distance="+IJ.d2s(lensDistortionParameters.distance, 3));
	        	System.out.println("lensDistortionParameters.x0="+      IJ.d2s(lensDistortionParameters.x0, 3));
	        	System.out.println("lensDistortionParameters.y0="+      IJ.d2s(lensDistortionParameters.y0, 3));
	        	System.out.println("lensDistortionParameters.z0="+      IJ.d2s(lensDistortionParameters.z0, 3));
	        	System.out.println("lensDistortionParameters.pitch="+   IJ.d2s(lensDistortionParameters.pitch, 3));
	        	System.out.println("lensDistortionParameters.yaw="+IJ.d2s(lensDistortionParameters.yaw, 3));
	        	System.out.println("lensDistortionParameters.roll="+IJ.d2s(lensDistortionParameters.roll, 3));
	        	System.out.println("lensDistortionParameters.focalLength="+IJ.d2s(lensDistortionParameters.focalLength, 3));
	        	System.out.println("lensDistortionParameters.px0="+IJ.d2s(lensDistortionParameters.px0, 3));
	        	System.out.println("lensDistortionParameters.py0="+IJ.d2s(lensDistortionParameters.py0, 3));
	        	System.out.println("lensDistortionParameters.distortionA8="+IJ.d2s(lensDistortionParameters.distortionA8, 5));
	        	System.out.println("lensDistortionParameters.distortionA7="+IJ.d2s(lensDistortionParameters.distortionA7, 5));
	        	System.out.println("lensDistortionParameters.distortionA6="+IJ.d2s(lensDistortionParameters.distortionA6, 5));
	        	System.out.println("lensDistortionParameters.distortionA5="+IJ.d2s(lensDistortionParameters.distortionA5, 5));
	        	System.out.println("lensDistortionParameters.distortionA="+IJ.d2s(lensDistortionParameters.distortionA, 5));
	        	System.out.println("lensDistortionParameters.distortionB="+IJ.d2s(lensDistortionParameters.distortionB, 5));
	        	System.out.println("lensDistortionParameters.distortionC="+IJ.d2s(lensDistortionParameters.distortionC, 5));
	        	for (int i=0;i<lensDistortionParameters.r_xy.length;i++){
	        		System.out.println("lensDistortionParameters.r_xy["+i+"][0]"+IJ.d2s(lensDistortionParameters.r_xy[i][0], 5));
	        		System.out.println("lensDistortionParameters.r_xy["+i+"][1]"+IJ.d2s(lensDistortionParameters.r_xy[i][1], 5));

	        	}
	        	for (int i=0;i<lensDistortionParameters.r_od.length;i++){
	        		System.out.println("lensDistortionParameters.r_od["+i+"][0]"+IJ.d2s(lensDistortionParameters.r_od[i][0], 5));
	        		System.out.println("lensDistortionParameters.r_od["+i+"][1]"+IJ.d2s(lensDistortionParameters.r_od[i][1], 5));
	        	}
	        }
	        if (!calculateDerivatives) return;
	/* Calculate all derivativs as a matrix.
	 *  Input parameters (columns):
	0    	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
	1    	public double radius;  // mm, distance from the rotation axis
	2   	public double height;       // mm, up (was downwards?) - from the origin point
	3   	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise
	4   	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up
	5   	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
	6   	public double goniometerHorizontal; // goniometer rotation around "horizontal" axis (tilting from the target - positive)
	7   	public double goniometerAxial; // goniometer rotation around Eyesis axis (clockwise in plan - positive
	8   	public double interAxisDistance; // distance in mm between two goniometer axes
	9    	public double interAxisAngle;    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
	10   	public double horAxisErrPhi;   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)
	11   	public double horAxisErrPsi;   // angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
	Two new parameters
	12    	public double entrancePupilForward; // common to all lenses - distance from the sensor to the lens entrance pupil
	13    	public double centerAboveHorizontal; // camera center distance along camera axis above the closest point to horizontal rotation axis (adds to
	14(12)  x	public double [] GXYZ=null; // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
	15(13)  y
	16(14)  z
	17(15)	public double focalLength=4.5;
	18(16)	public double px0=1296.0;          // center of the lens on the sensor, pixels
	19(17)		public double py0=968.0;           // center of the lens on the sensor, pixels
	20(18)		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
	21(19)		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
	22(20)		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
	23(21)		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
	24(22)		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	25(23)		public double distortionB=0.0; // r^3
	26(24)		public double distortionC=0.0; // r^2
    27          Distortion center shift X for r^3"
    28          Distortion center shift Y for r^3"
    29          Distortion center shift X for r^4"
    30          Distortion center shift Y for r^4"
    31          Distortion center shift X for r^5"
    32          Distortion center shift Y for r^5"
    33          Distortion center shift X for r^6"
    34          Distortion center shift Y for r^6"
    35          Distortion center shift X for r^7"
    36          Distortion center shift Y for r^7"
    37          Distortion center shift X for r^8"
    38          Distortion center shift Y for r^8"
    39          Orthogonal elongation for r^2"
    40          Diagonal   elongation for r^2"
    41          Orthogonal elongation for r^3"
    42          Diagonal   elongation for r^3"
    43          Orthogonal elongation for r^4"
    44          Diagonal   elongation for r^4"
    45          Orthogonal elongation for r^5"
    46          Diagonal   elongation for r^5"
    47          Orthogonal elongation for r^6"
    48          Diagonal   elongation for r^6"
    49          Orthogonal elongation for r^7"
    50          Diagonal   elongation for r^7"
    51          Orthogonal elongation for r^8"
    52          Diagonal   elongation for r^8"

	 * Output parameters (rows):
	0		public double x0=0;      // lens axis from pattern center, mm (to the right)
	1		public double y0=0;      // lens axis from pattern center, mm (down)
	2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
	3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
	4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
	5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise
	6		public double focalLength=4.5;
	7		public double px0=1296.0;           // center of the lens on the sensor, pixels
	8		public double py0=968.0;           // center of the lens on the sensor, pixels
	 		public double distortionRadius=  2.8512; // mm - half width of the sensor
	9		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
	10		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
	11		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
	12		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
	13		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	14		public double distortionB=0.0; // r^3
	15		public double distortionC=0.0; // r^2
    16          Distortion center shift X for r^3"
    17          Distortion center shift Y for r^3"
    18          Distortion center shift X for r^4"
    19          Distortion center shift Y for r^4"
    20          Distortion center shift X for r^5"
    21          Distortion center shift Y for r^5"
    22          Distortion center shift X for r^6"
    23          Distortion center shift Y for r^6"
    24          Distortion center shift X for r^7"
    25          Distortion center shift Y for r^7"
    26          Distortion center shift X for r^8"
    27          Distortion center shift Y for r^8"
    28          Orthogonal elongation for r^2"
    29          Diagonal   elongation for r^2"
    30          Orthogonal elongation for r^3"
    31          Diagonal   elongation for r^3"
    32          Orthogonal elongation for r^4"
    33          Diagonal   elongation for r^4"
    34          Orthogonal elongation for r^5"
    35          Diagonal   elongation for r^5"
    36          Orthogonal elongation for r^6"
    37          Diagonal   elongation for r^6"
    38          Orthogonal elongation for r^7"
    39          Diagonal   elongation for r^7"
    40          Orthogonal elongation for r^8"
    41          Diagonal   elongation for r^8"
 */
//	        interParameterDerivatives=new double[getNumInputs()][]; //partial derivative matrix from subcamera-camera-goniometer to single camera (12x21)

	//calculate dMA_azimuth
	//calculate dMB_azimuth
	/*
	// MA=R8*R7*R6*R5*R4*R3*R2*R1;
	// MB=T3+(R8*R7*R6*R5*(T2+R4*T1)); - old
	// MB=T3+(R8*R7*R6*R5*(T2+R4*(T1+R3*R2*R1*T0)));
	        Matrix MA=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1)))))));  old
	        Matrix MB=T3.plus(R8.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1.plus(R3.times(R2.times(R1.times(T0)))) )))))));
	3) rotate by -(azimuth+phi) around C2Y:Vc3= R3*Vc2
	| Xc3 |   | cos(azimuth+phi)    0   sin(azimuth+phi)   |   |Xc2|
	| Yc3 | = |     0               1         0            | * |Yc2|
	| Zc3 |   | -sin(azimuth+phi)   0   cos(azimuth+phi)   |   |Zc2|
	    	double [][] aR3={{cAZP,0.0,sAZP},{0.0,1.0,0.0},{-sAZP,0.0,cAZP}};
	    	Matrix R3=new Matrix(aR3);
	4) Now axes are aligned, just translate to get to eyesis coordinates: Vey= T1+Vc3
	| Xey |   |  r * sin (azimuth)  |   |Xc3|
	| Yey | = |             height  | + |Yc3|
	| Zey |   |  r * cos (azimuth)  |   |Zc3|
	    	double [][] aT1={{subCam.radius*sAZ},{subCam.height},{subCam.radius*cAZ}};
	    	Matrix T1=new Matrix(aT1);

	        Matrix
	// Make a function MA, MB - >parameters (column) - reuse it above and for each    interParameterDerivatives row
	Matrix dMA_azimuth=R8*R7*R6*R5*R4*dR3_azimuth*R2*R1;
	Matrix dMB_azimuth=T3+(R8*R7*R6*R5*(R4*dT1_azimuth));

	Use extrinsicParams=parametersFromMAMB(dMA_azimuth,dMB_azimuth);
		       		if (this.debugLevel>3) {
		    			System.out.println("calculateFxAndJacobian->calcPartialDerivatives("+IJ.d2s(targetXYZ[fullIndex][0],2)+","+
		    					IJ.d2s(targetXYZ[fullIndex][1],2)+","+
		    					IJ.d2s(targetXYZ[fullIndex][2],2)+" ("+calcJacobian+") -> "+
		    					IJ.d2s(derivatives14[0][0],2)+"/"+IJ.d2s(derivatives14[0][1],2));
		    		}

	Which parameters affect which matrices
	                                               R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 || T0 | T1 | T2 | T3 |
	0    	public double azimuth_cyl;          //    |    | +  |    |    |    |    |    ||    | +  |    |    |
	1    	public double radius_cyl;           //    |    |    |    |    |    |    |    ||    | +  |    |    |
	0    	public double right_cart;           //    |    |    |    |    |    |    |    ||    | +  |    |    |
	1    	public double forward_cart;         //    |    |    |    |    |    |    |    ||    | +  |    |    |
	2   	public double height;               //    |    |    |    |    |    |    |    ||    | +  |    |    |
	3   	public double phi_cyl;              //    |    | +  |    |    |    |    |    ||    |    |    |    |
	3   	public double head_cart;            //    |    | +  |    |    |    |    |    ||    |    |    |    |
	4   	public double theta;                //    | +  |    |    |    |    |    |    ||    |    |    |    |
	5   	public double psi;                  // +  |    |    |    |    |    |    |    ||    |    |    |    |
	6   	public double goniometerHorizontal; //    |    |    |    |    | +  |    |    ||    |    |    |    |
	7   	public double goniometerAxial;      //    |    |    | +  |    |    |    |    ||    |    |    |    |
	8   	public double interAxisDistance;    //    |    |    |    |    |    |    |    ||    |    | +  |    |
	9    	public double interAxisAngle;       //    |    |    |    | +  |    |    |    ||    |    |    |    |
	10   	public double horAxisErrPhi;        //    |    |    |    |    |    |    | +  ||    |    |    |    |
	11   	public double horAxisErrPsi;        //    |    |    |    |    |    | +  |    ||    |    |    |    |
	12      public double entrancePupilForward  //    |    |    |    |    |    |    |    || +  |    |    |    |
	13      public double centerAboveHorizontal //    |    |    |    |    |    |    |    ||    | +  |    |    |
	14  x	public double [] GXYZ=null;         //    |    |    |    |    |    |    |    ||    |    |    | +  |
	15  y                                       //    |    |    |    |    |    |    |    ||    |    |    | +  |
	16  z                                       //    |    |    |    |    |    |    |    ||    |    |    | +  |

	 */
	// Below can be optimized with common intermediate results
	   		if (this.debugLevel>2) {
	   			for (int i=0;i<parVect.length;i++){
	   				System.out.println("calcInterParamers(): parVect["+i+"]="+parVect[i]);
	   			}
	   		}
	//0    	public double right_cart; // right displacement of the lens entrance pupil center, mm
	        if (mask[0]) {
	        	if (cartesian) {
	        		double [][] adT1_right_cart={{1.0},{0.0},{0.0}};
	        		Matrix dT1_right_cart=new Matrix(adT1_right_cart);
		        	Matrix dMA_right_cart=new Matrix(3,3,0.0); // zero
	        		Matrix dMB_right_cart=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_right_cart)))));
		        	interParameterDerivatives[0]=d_parametersFromMAMB(dMA_right_cart,dMB_right_cart,MA,MB,false);
		    		if (this.debugLevel>2) {
		    			System.out.println("dMA_right_cart:");
		    			dMA_right_cart.print(10, 5);
		    			System.out.println("dMB_right_cart:");
		    			dMB_right_cart.print(10, 5);
		    			System.out.println("interParameterDerivatives[0]="+sprintfArray(interParameterDerivatives[0]));
		    		}
	        	} else {
    //0    	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
	        		double [][] adR3_azimuth_cyl={{-sAZP,0.0,cAZP},{0.0,0.0,0.0},{-cAZP,0.0,-sAZP}};
	        		Matrix dR3_azimuth_cyl=new Matrix(adR3_azimuth_cyl);
	        		double [][] adT1_azimuth_cyl={{radius_cyl*cAZ_cyl},{0.0},{-radius_cyl*sAZ_cyl}}; //{{subCam.radius*cAZ},{subCam.height},{-subCam.radius*sAZ}}
	        		Matrix dT1_azimuth_cyl=new Matrix(adT1_azimuth_cyl);

	        		Matrix dMA_azimuth_cyl=R8.times(R7.times(R6.times(R5.times(R4.times(dR3_azimuth_cyl.times(R2.times(R1)))))));
	        		Matrix dMB0_azimuth_cyl=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_azimuth_cyl)))));
	        		Matrix dMB_azimuth_cyl=dMB0_azimuth_cyl.plus(dMA_azimuth_cyl.times(T0)); // new term
	        		interParameterDerivatives[0]=d_parametersFromMAMB(dMA_azimuth_cyl,dMB_azimuth_cyl,MA,MB,true); // all after 6 are 0;
	        		if (this.debugLevel>2) {
	        			System.out.println("dMA_azimuth_cyl:");
	        			dMA_azimuth_cyl.print(10, 5);
	        			System.out.println("dMB_azimuth_cyl:");
	        			dMB_azimuth_cyl.print(10, 5);
	        			System.out.println("interParameterDerivatives[0]="+sprintfArray(interParameterDerivatives[0]));
	        		}
	        	}
	        } else interParameterDerivatives[0]=null;
	//1    	public double radius;  // mm, distance from the rotation axis
	        if (mask[1]) {
	        	if (cartesian) {
	        		double [][] adT1_forward_cart={{0.0},{0.0},{1.0}}; //{{subCam.radius*sAZ},{0.0},{subCam.radius*cAZ}}
	        		Matrix dT1_forward_cart=new Matrix(adT1_forward_cart);
	        		Matrix dMA_forward_cart=new Matrix(3,3,0.0);
	        		Matrix dMB_forward_cart=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_forward_cart)))));
	        		interParameterDerivatives[1]=d_parametersFromMAMB(dMA_forward_cart,dMB_forward_cart,MA,MB,false); // all after 6 are 0;
	        		if (this.debugLevel>2) {
	        			System.out.println("dMA_forward_cart:");
	        			dMA_forward_cart.print(10, 5);
	        			System.out.println("dMB_forward_cart:");
	        			dMB_forward_cart.print(10, 5);
	        			System.out.println("interParameterDerivatives[1]="+sprintfArray(interParameterDerivatives[1]));
	        		}
	        	} else {
	        		double [][] adT1_radius_cyl={{sAZ_cyl},{0.0},{cAZ_cyl}}; //{{subCam.radius*sAZ},{0.0},{subCam.radius*cAZ}}
	        		Matrix dT1_radius_cyl=new Matrix(adT1_radius_cyl);
	        		Matrix dMA_radius_cyl=new Matrix(3,3,0.0);
	        		Matrix dMB_radius_cyl=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_radius_cyl)))));
	        		interParameterDerivatives[1]=d_parametersFromMAMB(dMA_radius_cyl,dMB_radius_cyl,MA,MB,false); // all after 6 are 0;
	        		if (this.debugLevel>2) {
	        			System.out.println("dMA_radius_cyl:");
	        			dMA_radius_cyl.print(10, 5);
	        			System.out.println("dMB_radius_cyl:");
	        			dMB_radius_cyl.print(10, 5);
	        			System.out.println("interParameterDerivatives[1]="+sprintfArray(interParameterDerivatives[1]));
	        		}
	        	}
	        } else interParameterDerivatives[1]=null;
	//2   	public double height;       // mm, downwards - from the origin point
	        if (mask[2]) {
	        	double [][] adT1_height={{0.0},{1.0},{0.0}};
	        	Matrix dT1_height=new Matrix(adT1_height);
	        	Matrix dMA_height=new Matrix(3,3,0.0);
	        	Matrix dMB_height=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_height)))));
	        	interParameterDerivatives[2]=d_parametersFromMAMB(dMA_height,dMB_height,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_height:");
	    			dMA_height.print(10, 5);
	    			System.out.println("dMB_height:");
	    			dMB_height.print(10, 5);
	    			System.out.println("interParameterDerivatives[2]="+sprintfArray(interParameterDerivatives[2]));
	    		}
	        } else interParameterDerivatives[2]=null;
	//3   	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise
	        if (mask[3]) { // here the same for cartesian
	        	double [][] adR3_phi={{-sAZP,0.0,cAZP},{0.0,0.0,0.0},{-cAZP,0.0,-sAZP}}; // same as adR3_azimuth
	        	Matrix dR3_phi=new Matrix(adR3_phi); // same as dR3_azimuth
	        	Matrix dMA_phi=R8.times(R7.times(R6.times(R5.times(R4.times(dR3_phi.times(R2.times(R1))))))); //same as dMA_azimuth
//	        	Matrix dMB_phi=new Matrix(3,1,0.0);
	        	Matrix dMB_phi=dMA_phi.times(T0); // new term
	        	interParameterDerivatives[3]=d_parametersFromMAMB(dMA_phi,dMB_phi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_phi:");
	    			dMA_phi.print(10, 5);
	    			System.out.println("dMB_phi:");
	    			dMB_phi.print(10, 5);
	    			System.out.println("interParameterDerivatives[3]="+sprintfArray(interParameterDerivatives[3]));
	    		}
	        } else interParameterDerivatives[3]=null;
	//4   	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up
	        if (mask[4]) {
	        	double [][] adR2_theta={{0.0,0.0,0.0},{0.0,-sTH,cTH},{0.0,-cTH,-sTH}};
	        	Matrix dR2_theta=new Matrix(adR2_theta);
	        	Matrix dMA_theta=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(dR2_theta.times(R1)))))));
//	        	Matrix dMB_theta=new Matrix(3,1,0.0);
	        	Matrix dMB_theta=dMA_theta.times(T0); // new term
	        	interParameterDerivatives[4]=d_parametersFromMAMB(dMA_theta,dMB_theta,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_theta:");
	    			dMA_theta.print(10, 5);
	    			System.out.println("dMB_theta:");
	    			dMB_theta.print(10, 5);
	    			System.out.println("interParameterDerivatives[4]="+sprintfArray(interParameterDerivatives[4]));
	    		}
	        } else interParameterDerivatives[4]=null;
	//5   	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
	        if (mask[5]) {
	        	double [][] adR1_psi={{-sPS,cPS,0.0},{-cPS,-sPS,0.0},{0.0,0.0,0.0}};
	        	Matrix dR1_psi=new Matrix(adR1_psi);
	        	Matrix dMA_psi=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(dR1_psi)))))));
//	        	Matrix dMB_psi=new Matrix(3,1,0.0);
	        	Matrix dMB_psi=dMA_psi.times(T0); // new term
	        	interParameterDerivatives[5]=d_parametersFromMAMB(dMA_psi,dMB_psi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	/*
	    			System.out.print("R1:");
	    			R1.print(10, 5);
	    			System.out.print("dR1_psi:");
	    			dR1_psi.print(10, 5);
	    			Matrix R82_psi=R8.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2))))));
	    			System.out.print("R82_psi:");
	    			R82_psi.print(10, 5);
	    			*/
	    			System.out.print("dMA_psi:");
	    			dMA_psi.print(10, 5);
	    			System.out.print("dMB_psi:");
	    			dMB_psi.print(10, 5);
	    			System.out.println("interParameterDerivatives[5]="+sprintfArray(interParameterDerivatives[5]));
	    		}
	        } else interParameterDerivatives[5]=null;
	//6   	public double goniometerHorizontal; // goniometer rotation around "horizontal" axis (tilting from the target - positive)
	        if (mask[6]) {
	/* define for isTripod */
	        	double [][] adR6_goniometerHorizontal_tripod=    {{-sGH,0.0,cGH},{0.0, 0.0,0.0},{-cGH, 0.0,-sGH}};
	        	double [][] adR6_goniometerHorizontal_goniometer={{ 0.0,0.0,0.0},{0.0,-sGH,cGH},{ 0.0,-cGH,-sGH}};
	        	Matrix dR6_goniometerHorizontal=new Matrix(isTripod?adR6_goniometerHorizontal_tripod:adR6_goniometerHorizontal_goniometer);
	        	Matrix dMA_goniometerHorizontal=R8.times(R7.times(dR6_goniometerHorizontal.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_goniometerHorizontal=R8.times(R7.times(dR6_goniometerHorizontal.times(R5.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[6]=d_parametersFromMAMB(dMA_goniometerHorizontal,dMB_goniometerHorizontal,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_goniometerHorizontal:");
	    			dMA_goniometerHorizontal.print(10, 5);
	    			System.out.println("dMB_goniometerHorizontal:");
	    			dMB_goniometerHorizontal.print(10, 5);
	    			System.out.println("interParameterDerivatives[6]="+sprintfArray(interParameterDerivatives[6]));
	    		}
	        } else interParameterDerivatives[6]=null;
	//7   	public double goniometerAxial; // goniometer rotation around Eyesis axis (clockwise in plan - positive
	        if (mask[7]) {
	// define for isTripod
	        	double [][] adR4_goniometerAxial_tripod=    {{ 0.0,0.0,0.0},{0.0,-sGA,cGA},{ 0.0,-cGA,-sGA}};
	        	double [][] adR4_goniometerAxial_goniometer={{-sGA,0.0,cGA},{0.0, 0.0,0.0},{-cGA, 0.0,-sGA}};
	        	Matrix dR4_goniometerAxial=new Matrix(isTripod?adR4_goniometerAxial_tripod:adR4_goniometerAxial_goniometer);
	        	Matrix dMA_goniometerAxial=R8.times(R7.times(R6.times(R5.times(dR4_goniometerAxial.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_goniometerAxial=R8.times(R7.times(R6.times(R5.times(dR4_goniometerAxial.times(T1)))));
	        	interParameterDerivatives[7]=d_parametersFromMAMB(dMA_goniometerAxial,dMB_goniometerAxial,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_goniometerAxial:");
	    			dMA_goniometerAxial.print(10, 5);
	    			System.out.println("dMB_goniometerAxial:");
	    			dMB_goniometerAxial.print(10, 5);
	    			System.out.println("interParameterDerivatives[7]="+sprintfArray(interParameterDerivatives[7]));
	    		}
	        } else interParameterDerivatives[7]=null;
	//8   	public double interAxisDistance; // distance in mm between two goniometer axes
	        if (mask[8]) {
	        	double [][] adT2_interAxisDistance={{0.0},{0.0},{1.0}};
	        	Matrix dT2_interAxisDistance=new Matrix(adT2_interAxisDistance);
	        	Matrix dMA_interAxisDistance=new Matrix(3,3,0.0);
	        	Matrix dMB_interAxisDistance=R8.times(R7.times(R6.times(R5.times(dT2_interAxisDistance))));
	        	interParameterDerivatives[8]=d_parametersFromMAMB(dMA_interAxisDistance,dMB_interAxisDistance,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_interAxisDistance:");
	    			dMA_interAxisDistance.print(10, 5);
	    			System.out.println("dMB_interAxisDistance:");
	    			dMB_interAxisDistance.print(10, 5);
	    			System.out.println("interParameterDerivatives[8]="+sprintfArray(interParameterDerivatives[8]));
	    		}
	        } else interParameterDerivatives[8]=null;
	//9    	public double interAxisAngle;    // angle in degrees between two goniometer axes minus 90. negative if "vertical" axis is rotated
	        if (mask[9]) {
	        	double [][] adR5_interAxisAngle={{-sGIAA,cGIAA,0.0},{-cGIAA,-sGIAA,0.0},{0.0,0.0,0.0}};
	        	Matrix dR5_interAxisAngle=new Matrix(adR5_interAxisAngle);
	        	Matrix dMA_interAxisAngle=R8.times(R7.times(R6.times(dR5_interAxisAngle.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_interAxisAngle=R8.times(R7.times(R6.times(dR5_interAxisAngle.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[9]=d_parametersFromMAMB(dMA_interAxisAngle,dMB_interAxisAngle,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_interAxisAngle:");
	    			dMA_interAxisAngle.print(10, 5);
	    			System.out.println("dMB_interAxisAngle:");
	    			dMB_interAxisAngle.print(10, 5);
	    			System.out.println("interParameterDerivatives[9]="+sprintfArray(interParameterDerivatives[9]));
	    		}
	        } else interParameterDerivatives[9]=null;
	//10   	public double horAxisErrPhi;   // angle in degrees "horizontal" goniometer axis is rotated around target Y axis from target X axis (CW)

	// change sHAEPH to rotate like theta /1 0 0 /0 c s/ 0 -s c/ (derivative -/0 0 0 /0 -s c/ 0 -c -s/
	        if (mask[10]) {
	        	double [][] adR8_horAxisErrPhi_tripod=    {{0.0,   0.0,    0.0},{0.0,-sHAEPH,cHAEPH},{    0.0,-cHAEPH,-sHAEPH}};
	        	double [][] adR8_horAxisErrPhi_goniometer={{-sHAEPH,0.0,cHAEPH},{0.0,    0.0,   0.0},{-cHAEPH,    0.0,-sHAEPH}};
	        	Matrix dR8_horAxisErrPhi=new Matrix(isTripod?adR8_horAxisErrPhi_tripod:adR8_horAxisErrPhi_goniometer);
	        	Matrix dMA_horAxisErrPhi=dR8_horAxisErrPhi.times(R7.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_horAxisErrPhi=dR8_horAxisErrPhi.times(R7.times(R6.times(R5.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[10]=d_parametersFromMAMB(dMA_horAxisErrPhi,dMB_horAxisErrPhi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_horAxisErrPhi:");
	    			dMA_horAxisErrPhi.print(10, 5);
	    			System.out.println("dMB_horAxisErrPhi:");
	    			dMB_horAxisErrPhi.print(10, 5);
	    			System.out.println("interParameterDerivatives[10]="+sprintfArray(interParameterDerivatives[10]));
	    		}
	        } else interParameterDerivatives[10]=null;
	//11   	public double horAxisErrPsi;   // angle in degrees "horizontal" goniometer axis is rotated around moving X axis (up)
	        if (mask[11]) {
	        	double [][] adR7_horAxisErrPsi={{-sHAEPS,cHAEPS,0.0},{-cHAEPS,-sHAEPS,0.0},{0.0,0.0,0.0}};
	        	Matrix dR7_horAxisErrPsi=new Matrix(adR7_horAxisErrPsi);
	        	Matrix dMA_horAxisErrPsi=R8.times(dR7_horAxisErrPsi.times(R6.times(R5.times(R4.times(R3.times(R2.times(R1)))))));
	        	Matrix dMB_horAxisErrPsi=R8.times(dR7_horAxisErrPsi.times(R6.times(R5.times(T2.plus(R4.times(T1))))));
	        	interParameterDerivatives[11]=d_parametersFromMAMB(dMA_horAxisErrPsi,dMB_horAxisErrPsi,MA,MB,true); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_horAxisErrPsi:");
	    			dMA_horAxisErrPsi.print(10, 5);
	    			System.out.println("dMB_horAxisErrPsi:");
	    			dMB_horAxisErrPsi.print(10, 5);
	    			System.out.println("interParameterDerivatives[11]="+sprintfArray(interParameterDerivatives[11]));
	    		}
	        } else interParameterDerivatives[11]=null;
//	                                          // R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 || T0 | T1 | T2 | T3 |
	//12   	public double entrancePupilForward  //    |    |    |    |    |    |    |    || +  |    |    |    |
	//13   	public double centerAboveHorizontal //    |    |    |    |    |    |    |    ||    | +  |    |    |

	        if (mask[12]) {
	        	double [][] adT0_entrancePupilForward={{0.0},{0.0},{1.0}};
	        	Matrix dT0_entrancePupilForward=new Matrix(adT0_entrancePupilForward);
	        	Matrix dMA_entrancePupilForward=new Matrix(3,3,0.0);
	        	Matrix dMB_entrancePupilForward=MA.times(dT0_entrancePupilForward);
	        	interParameterDerivatives[12]=d_parametersFromMAMB(dMA_entrancePupilForward,dMB_entrancePupilForward,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_centerAboveHorizontal:");
	    			dMA_entrancePupilForward.print(10, 5);
	    			System.out.println("dMB_centerAboveHorizontal:");
	    			dMB_entrancePupilForward.print(10, 5);
	    			System.out.println("interParameterDerivatives[12]="+sprintfArray(interParameterDerivatives[12]));
	    		}
	        } else interParameterDerivatives[12]=null;

	        if (mask[13]) {
	        	double [][] adT1_centerAboveHorizontal={{0.0},{1.0},{0.0}};
	        	Matrix dT1_centerAboveHorizontal=new Matrix(adT1_centerAboveHorizontal);
	        	Matrix dMA_centerAboveHorizontal=new Matrix(3,3,0.0);
	        	Matrix dMB_centerAboveHorizontal=R8.times(R7.times(R6.times(R5.times(R4.times(dT1_centerAboveHorizontal)))));
	        	interParameterDerivatives[13]=d_parametersFromMAMB(dMA_centerAboveHorizontal,dMB_centerAboveHorizontal,MA,MB,false); // all after 6 are 0;
	    		if (this.debugLevel>2) {
	    			System.out.println("dMA_centerAboveHorizontal:");
	    			dMA_centerAboveHorizontal.print(10, 5);
	    			System.out.println("dMB_centerAboveHorizontal:");
	    			dMB_centerAboveHorizontal.print(10, 5);
	    			System.out.println("interParameterDerivatives[13]="+sprintfArray(interParameterDerivatives[13]));
	    		}
	        } else interParameterDerivatives[13]=null;


	//14  x	public double [] GXYZ=null; // coordinates (in mm) of the goniometer horizontal axis closest to the moving one in target system
	        if (mask[14]) {
	        	double [][] adT3_GXYZ0={{1.0},{0.0},{0.0}};
	        	Matrix dMA_GXYZ0=new Matrix(3,3,0.0);
	        	Matrix dMB_GXYZ0=new Matrix(adT3_GXYZ0);
	        	interParameterDerivatives[14]=d_parametersFromMAMB(dMA_GXYZ0,dMB_GXYZ0,MA,MB,false); // all after 6 are 0;
	        	if (this.debugLevel>2) {
	        		System.out.println("dMA_GXYZ0:");
	        		dMA_GXYZ0.print(10, 5);
	        		System.out.println("dMA_GXYZ0:");
	        		dMB_GXYZ0.print(10, 5);
	        		System.out.println("interParameterDerivatives[14]="+sprintfArray(interParameterDerivatives[14]));
	        	}
	        } else interParameterDerivatives[14]=null;
	//13  y
	        if (mask[15]) {
//	    	double [][] adT3_GXYZ1={{0.0},{1.0},{0.0}};
	        	double [][] adT3_GXYZ1={{0.0},{-1.0},{0.0}}; // up - positive
//	        	double [][] adT3_GXYZ1={{0.0},{ 1.0},{0.0}}; // up - positive
	        	Matrix dMA_GXYZ1=new Matrix(3,3,0.0);
	        	Matrix dMB_GXYZ1=new Matrix(adT3_GXYZ1);
	        	interParameterDerivatives[15]=d_parametersFromMAMB(dMA_GXYZ1,dMB_GXYZ1,MA,MB,false); // all after 6 are 0;
	        	if (this.debugLevel>2) {
	        		System.out.println("dMA_GXYZ1:");
	        		dMA_GXYZ1.print(10, 5);
	        		System.out.println("dMB_GXYZ1:");
	        		dMB_GXYZ1.print(10, 5);
	        		System.out.println("interParameterDerivatives[15]="+sprintfArray(interParameterDerivatives[15]));
	        	}
	        } else interParameterDerivatives[15]=null;
	//14  z
	        if (mask[16]) {
//	    	double [][] adT3_GXYZ2={{0.0},{0.0},{1.0}};
	        	double [][] adT3_GXYZ2={{0.0},{0.0},{-1.0}}; // away - positive
	        	Matrix dMA_GXYZ2=new Matrix(3,3,0.0);
	        	Matrix dMB_GXYZ2=new Matrix(adT3_GXYZ2);
	        	interParameterDerivatives[16]=d_parametersFromMAMB(dMA_GXYZ2,dMB_GXYZ2,MA,MB,false); // all after 6 are 0;
	        	if (this.debugLevel>2) {
	        		System.out.println("dMA_GXYZ2:");
	        		dMA_GXYZ2.print(10, 5);
	        		System.out.println("dMB_GXYZ2:");
	        		dMB_GXYZ2.print(10, 5);
	        		System.out.println("interParameterDerivatives[16]="+sprintfArray(interParameterDerivatives[16]));
	        	}
	        } else interParameterDerivatives[16]=null;
	// now fill the rest, unchanged parameters
	/*
	    int numInputs=22;   //was 21  parameters in subcamera+...
	    int numOutputs=13;  //was 12  parameters in a single camera

	 */
	//17 (15)		public double focalLength=4.5;
	//18 (16)		public double px0=1296.0;          // center of the lens on the sensor, pixels
	//19 (17)		public double py0=968.0;           // center of the lens on the sensor, pixels
	//20 (18)		public double distortionA8=0.0; // r^8 (normalized to focal length or to sensor half width?)
	//21 (19)		public double distortionA7=0.0; // r^7 (normalized to focal length or to sensor half width?)
	//22 (20)		public double distortionA6=0.0; // r^6 (normalized to focal length or to sensor half width?)
	//23 (21)		public double distortionA5=0.0; // r^5 (normalized to focal length or to sensor half width?)
	//24 (22)		public double distortionA=0.0; // r^4 (normalized to focal length or to sensor half width?)
	//25 (23)		public double distortionB=0.0; // r^3
	//26 (24)		public double distortionC=0.0; // r^2
	//27..52 - non-radial parameters

//	        for (int inpPar=15; inpPar<getNumInputs(); inpPar++){
	        for (int inpPar=17; inpPar<getNumInputs(); inpPar++){ // OK with A8. Should be OK with non-radial parameters ((27..52) also
	        	if (mask[inpPar]){
	        		interParameterDerivatives[inpPar]=new double[getNumOutputs()];
	        		for (int outPar=0; outPar<getNumOutputs(); outPar++){
	        			interParameterDerivatives[inpPar][outPar]=((getNumOutputs()-outPar)==(getNumInputs()-inpPar))?1.0:0.0;
	            		if (this.debugLevel>2) {
	            			System.out.println("interParameterDerivatives["+inpPar+"]["+outPar+"]="+interParameterDerivatives[inpPar][outPar]);
	            		}
	        		}
	        	} else interParameterDerivatives[inpPar]=null;
	        }
	    }


	     public double [] parametersFromMAMB(Matrix MA, Matrix MB){
	    	double [] result= new double [getNumOutputs()]; // only first 6 are used
	    	for (int i=6; i<getNumOutputs();i++) result[i]=0.0;
	/*
	0		public double x0=0;      // lens axis from pattern center, mm (to the right)
	1		public double y0=0;      // lens axis from pattern center, mm (down)
	2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
	3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
	4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
	5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise

	        Vt=MB+MA*Vc
	        calculate X0,Y0 (Z0==0), dist, phi, theta, psi from MB and MA
	        | MA00 MA01 MA02 |
	        | MA10 MA11 MA12 |
	        | MA20 MA21 MA22 |
	        ===================
	        | MB0 |
	        | MB1 |
	        | MB2 |
	        Take point Pc on a sub-camera axis and on the target Zt=0 plane: [0,0,dist]
	        MA*Pc+MB= dist*[MA02,MA12,MA22]+MB =[Tx0,Ty0,0]
	        MA02*dist+MB0=Tx0
	        MA12*dist+MB1=Ty0
	        MA22*dist+MB2=0
	        dist=-MB2/MA22;
	        Tx0=MA02*(-MB2/MA22)+MB0
	        Ty0=MA12*(-MB2/MA22)+MB1
	        Tx0=MB0-MB2*MA02/MA22
	        Ty0=MB1-MB2*MA12/MA22
	        */
	    	result[2]=-MB.get(2,0)/MA.get(2,2); // distance
	    	result[0]= MB.get(0,0)-MB.get(2,0)*MA.get(0,2)/MA.get(2,2); // x0
	    	result[1]= MB.get(1,0)-MB.get(2,0)*MA.get(1,2)/MA.get(2,2); // y0
	        /*
	        // now find phi, theta, psi
	        MA - rotation from camera to target, transp(MA) - rotation from target to camera - same as rot(phi, theta,psi)
	        MT=transp(MA),
	        */
	        /*
	        MA[1,2]= sin(theta)
	        MA[0,2]= cos(theta)*sin(phi)
	        MA[2,2]= cos(theta)*cos(phi)
	        MA[1,0]=-cos(theta)*sin(psi)
	        MA[1,1]= cos(theta)*cos(psi)
	        theta=arcsin(MA[1,2])
	        phi=  atan2(MA[0,2],MA[2,2])
	        psi= -atan2(MA[1,0],MA[1,1])
	         */
	        result[4]=180.0/Math.PI* Math.asin(MA.get(1, 2));  //pitch
	        result[3]=  180.0/Math.PI* Math.atan2(MA.get(0,2),MA.get(2, 2));//yaw
	        result[5]=-180.0/Math.PI* Math.atan2(MA.get(1, 0),MA.get(1, 1));//roll
	    	return result;
	    }
	    /**
	     *
	     * @param d_MA differential of the rotational matrix MA by some parameter
	     * @param d_MB differential of the translational matrix MB by some parameter
	     * @param MA - rotation matrix
	     * @param MB - translation matrix
	     * @param isAngle - when true, the partial derivative is for angles, d_MA, d_MB should be divided by 180/pi
	     * @return differentials of the {x0,y0,dist,phi,theta,psi} by that parameter
	     */
	    public double [] d_parametersFromMAMB(Matrix d_MA, Matrix d_MB, Matrix MA, Matrix MB, boolean isAngle){
	    	double [] result= new double [getNumOutputs()]; // only first 6 are used, rest are 0
	    	Arrays.fill(result,6, result.length, 0.0);
//	    	for (int i=6; i<getNumOutputs();i++) result[i]=0.0;
	/*
	0		public double x0=0;      // lens axis from pattern center, mm (to the right)
	1		public double y0=0;      // lens axis from pattern center, mm (down)
	2		public double distance=2360; // distance from the lens input pupil to the pattern plane along the camera axis, mm
	3		public double yaw=0.0;   // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - clockwise from top
	4		public double pitch=0.0; // angle in degrees from perpendicular to the pattern, 0 - towards wall, positive - up
	5		public double roll=0.0;  // angle in degrees rotation around camera optical axis (perpendicular to pattern if yaw==0, pitch==0), positive - clockwise

	        Vt=MB+MA*Vc
	        calculate X0,Y0 (Z0==0), dist, phi, theta, psi from MB and MA
	        | MA00 MA01 MA02 |
	        | MA10 MA11 MA12 |
	        | MA20 MA21 MA22 |
	        ===================
	        | MB0 |
	        | MB1 |
	        | MB2 |
	        Take point Pc on a sub-camera axis and on the target Zt=0 plane: [0,0,dist]
	        MA*Pc+MB= dist*[MA02,MA12,MA22]+MB =[Tx0,Ty0,0]
	        MA02*dist+MB0=Tx0
	        MA12*dist+MB1=Ty0
	        MA22*dist+MB2=0
	        dist=-MB2/MA22;
	        Tx0=MB0-MB2*MA02/MA22
	        Ty0=MB1-MB2*MA12/MA22
	        */
//	    	result[2]=-MB.get(2,0)/MA.get(2,2); // distance
	/*
	        d_dist=-(d_MB2/MA22 - d_MA22*MB2/(MA22^2) ;
	 */
	    	double K=isAngle?(Math.PI/180):1.0;
	    	result[2]=K*(-d_MB.get(2,0)/MA.get(2,2)+
	    	d_MA.get(2,2)*MB.get(2,0)/(MA.get(2,2)*MA.get(2,2))); // d_distance
	    	/*
	        Tx0=MB0-MB2*MA02/MA22
	        d_Tx0=d_MB0 -  d_MB2*MA02/MA22 -d_MA02*MB2/MA22 +d_MA22*MB2*MA02/(MA22^2)
	    	 */
	    	result[0]= K*(d_MB.get(0,0) -
	    	d_MB.get(2,0)*MA.get(0,2)/MA.get(2,2) -
	    	d_MA.get(0,2)*MB.get(2,0)/MA.get(2,2) +
	    	d_MA.get(2,2)*MB.get(2,0)*MA.get(0,2)/(MA.get(2,2)*MA.get(2,2))); // x0

	    	/*
	        Ty0=MB0-MB2*MA02/MA22
	        d_Ty0=d_MB1 -  d_MB2*MA12/MA22 -d_MA12*MB2/MA22 +d_MA22*MB2*MA12/(MA22^2)
	    	 */
	    	result[1]= K*(d_MB.get(1,0) -
	    	d_MB.get(2,0)*MA.get(1,2)/MA.get(2,2) -
	    	d_MA.get(1,2)*MB.get(2,0)/MA.get(2,2) +
	    	d_MA.get(2,2)*MB.get(2,0)*MA.get(1,2)/(MA.get(2,2)*MA.get(2,2))); // y0
	        /*
	        // now find phi, theta, psi
	        MA - rotation from camera to target, transp(MA) - rotation from target to camera - same as rot(phi, theta,psi)
	        MT=transp(MA),
	        */
	        /*
	        MA[1,2]= sin(theta)
	        MA[0,2]= cos(theta)*sin(phi)
	        MA[2,2]= cos(theta)*cos(phi)
	        MA[1,0]=-cos(theta)*sin(psi)
	        MA[1,1]= cos(theta)*cos(psi)
	        theta=arcsin(MA[1,2])
	        phi=  atan2(MA[0,2],MA[2,2])
	        psi= -atan2(MA[1,0],MA[1,1])
	        arcsin(x)'= 1/sqrt(1-x^2)
	        arccos(x)'=-1/sqrt(1-x^2)
	        arctan(x)'= 1/sqrt(1+x^2)

	         */
	/*
	        theta=arcsin(MA12)
	        d_theta=d_MA12/sqrt(1-MA12^2)
	 */
	        result[4]=K*d_MA.get(1, 2)*180.0/Math.PI/Math.sqrt(1.0-MA.get(1, 2)*MA.get(1, 2));  //pitch
	/*
	        phi=  atan2(MA02,MA22)
	        d_phi=(d_MA02*MA22-d_MA22*MA02) / (MA22^2+MA02^2)
	 */

	        result[3]=  K*180.0/Math.PI*(d_MA.get(0,2)*MA.get(2,2) - d_MA.get(2,2)*MA.get(0,2)) /
	          ((MA.get(2,2)*MA.get(2,2)) + (MA.get(0,2)*MA.get(0,2)));//yaw
	/*
	        psi=  -atan2(MA10,MA11)
	        d_psi=-(d_MA10*MA11-d_MA11*MA10) / (MA10^2+MA11^2)
	*/
	        result[5]=  -K* 180.0/Math.PI*(d_MA.get(1,0)*MA.get(1,1) - d_MA.get(1,1)*MA.get(1,0)) /
	          ((MA.get(1,1)*MA.get(1,1)) + (MA.get(1,0)*MA.get(1,0)));//roll
	    	return result;
	    }

	    public String sprintfArray(double []arr){
	    	String result="";
	    	for (int i=0;i<arr.length;i++) result+= ((i>0)?", ":"")+arr[i];
	    	return "["+result+"]";
	    }

		public int getNumInputs(){return numInputs;}
		public int getNumOutputs(){return numOutputs;}
	}
