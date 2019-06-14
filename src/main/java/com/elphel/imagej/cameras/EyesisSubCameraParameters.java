package com.elphel.imagej.cameras;
/*
 **
 ** EyesisSubCameraParameters.java
 **
 ** Copyright (C) 2011-2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  EyesisSubCameraParameters.java is free software: you can redistribute it and/or modify
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

    public class EyesisSubCameraParameters{
    	public int sensorWidth=      2592; // moving here from the camera as a whole
    	public int sensorHeight=     1936;
    	public int decimateMasks=       1;
    	// origin is on the rotation axis of the tube body closest to the goniometer horizontal axis
    	public boolean cartesian = false; // cartesian coordinates mode (false - cylindrical)
    	public int    lensDistortionModel=0;
    	public boolean enableNoLaser=true; // enable images for this channel w/o matched laser pointer
    	public double right;   // distance to the right (radius*sin(azimuth))
    	public double forward; // distance forward (radius*cos(azimuth))
    	public double heading;   // absolute heading in degrees (used in cartesian mode)

    	public double azimuth; // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    	public double radius;  // mm, distance from the rotation axis
    	public double height;       // mm, up - from the origin point
    	public double phi;     // degrees, optical axis from azimuth/r vector, clockwise (absolute if cylindrical = false)
    	public double theta;   // degrees, optical axis from the eyesis horizon, positive - up
    	public double psi;     // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
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
		public double channelWeightDefault=1.0;
		public double channelWeightCurrent=1.0;
		public int [][] defectsXY=null; // pixel defects coordinates list (starting with worst)
		public double [] defectsDiff=null; // pixel defects value (diff from average of neighbors), matching defectsXY

		final private double [][] r_xy_dflt={{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; // only 6, as for the first term delta x, delta y ==0
		final private double [][] r_od_dflt=   {{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; // ortho

		public double [][] r_xy=null; // only 6, as for the first term delta x, delta y ==0
		public double [][] r_od=null; // ortho
		// was not used before LWIR and import subsystem
		public int subcamera =   0;
		public int sensor_port = 0;
		public int subchannel =  0;


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

    	public EyesisSubCameraParameters(
    	    	int sensorWidth,
    	    	int sensorHeight,
    	    	int decimateMasks,
    			boolean cartesian,
    			int lensDistortionModel,
    			boolean enableNoLaser,
    	    	double right,   // distance to the right (radius*sin(azimuth))
    	    	double forward, // distance forward (radius*cos(azimuth))
    	    	double heading, // used in cartesian mode
    			double azimuth, // azimuth of the lens entrance pupil center, degrees, clockwise looking from top
    			double radius,  // mm, distance from the rotation axis
    			double height,  // mm, up from the origin point
    			double phi,     // degrees, optical axis from azimuth/r vector, clockwise
    			double theta,   // degrees, optical axis from the eyesis horizon, positive - up
    			double psi,      // degrees, rotation (of the sensor) around the optical axis. Positive if camera is rotated clockwise looking to the target
    			double focalLength,
    			double pixelSize,//um
    			double distortionRadius, //mm - half width of the sensor
    			double distortionA8, // r^8
    			double distortionA7, // r^7
    			double distortionA6, // r^5
    			double distortionA5, // r^5
    			double distortionA, // r^4 (normalized to focal length or to sensor half width?)
    			double distortionB, // r^3
    			double distortionC, // r^2
    			double px0,           // center of the lens on the sensor, pixels
    			double py0,           // center of the lens on the sensor, pixels
				double [][] r_xy,     // eccentricity for b,a,a5,a6,a7,a8
				double [][] r_od,     // elongation for c,b,a,a5,a6,a7,a8
    			double channelWeightDefault,
    			int subcamera,
    			int sensor_port,
    			int subchannel
    	){
    		this.decimateMasks = decimateMasks;
        	this.sensorWidth= sensorWidth;
        	this.sensorHeight=sensorHeight;
    		this.cartesian = cartesian;
    		this.lensDistortionModel=lensDistortionModel;
    		this.enableNoLaser=enableNoLaser;
    		this.right =   right;
    		this.forward = forward;
    		this.heading = heading;

    		this.azimuth=azimuth;
    		this.radius=radius;
    		this.height=height;
    		this.phi=phi;
    		this.theta=theta;
    		this.psi=psi;
    		this.focalLength=focalLength;
    		this.pixelSize=  pixelSize; //um
    		this.distortionRadius=  distortionRadius; // mm - half width of the sensor
    		this.distortionA8=distortionA8;
    		this.distortionA7=distortionA7;
    		this.distortionA6=distortionA6;
    		this.distortionA5=distortionA5;
    		this.distortionA=distortionA; // r^4 (normalized to focal length or to sensor half width?)
    		this.distortionB=distortionB; // r^3
    		this.distortionC=distortionC; // r^2
    		this.px0=px0;
    		this.py0=py0;
			if (r_xy==null) r_xy=r_xy_dflt;
			if (r_od==null) r_od=r_od_dflt;
			this.r_xy=new double [r_xy.length][2];
			for (int i=0;i<r_xy.length;i++)this.r_xy[i]=r_xy[i].clone();
			this.r_od=new double [r_od.length][2];
			for (int i=0;i<r_od.length;i++)this.r_od[i]=r_od[i].clone();
    		this.channelWeightDefault=channelWeightDefault;
    		this.channelWeightCurrent=this.channelWeightDefault;
    		this.defectsXY=null; // pixel defects coordinates list (starting with worst)
    		this.defectsDiff=null; // pixel defects value (diff from average of neighbors), matching defectsXY
			this.subcamera = subcamera;
			this.sensor_port = sensor_port;
			this.subchannel = subchannel;

    		updateCartesian(); // set alternative
    	}
    	// defects are not cloned!
    	@Override
		public EyesisSubCameraParameters clone() {
    		return new EyesisSubCameraParameters(
    				this.sensorWidth,
    				this.sensorHeight,
    				this.decimateMasks,
    	    		this.cartesian,
    				this.lensDistortionModel,
    				this.enableNoLaser,
    	    		this.right,
    	    		this.forward,
    	    		this.heading,
    				this.azimuth,
    				this.radius,
    				this.height,
    				this.phi,
    				this.theta,
    				this.psi,
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
    	    		this.px0,
    	    		this.py0,
    				this.r_xy,
    				this.r_od,
    	    		this.channelWeightDefault,
    				this.subcamera,
    				this.sensor_port,
    				this.subchannel
    				);
    	}
    	public void setDefaultNonRadial(){
			r_od=new double [r_od_dflt.length][2];
			for (int i=0;i<r_od.length;i++) r_od[i]=r_od_dflt[i].clone();
			r_xy=new double [r_xy_dflt.length][2];
			for (int i=0;i<r_xy.length;i++) r_xy[i]=r_xy_dflt[i].clone();
    	}
// TODO: add/restore new properties
    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"sensorWidth",         this.sensorWidth+"");
    		properties.setProperty(prefix+"sensorHeight",        this.sensorHeight+"");
      		properties.setProperty(prefix+"decimateMasks",       this.decimateMasks+"");
    		properties.setProperty(prefix+"cartesian",           this.cartesian+"");
    		properties.setProperty(prefix+"lensDistortionModel", this.lensDistortionModel+"");
    		properties.setProperty(prefix+"enableNoLaser",       this.enableNoLaser+"");
    		properties.setProperty(prefix+"right",               this.right+"");
    		properties.setProperty(prefix+"heading",             this.heading+"");
    		properties.setProperty(prefix+"forward",             this.forward+"");
    		properties.setProperty(prefix+"azimuth",             this.azimuth+"");
    		properties.setProperty(prefix+"radius",              this.radius+"");
    		properties.setProperty(prefix+"height",              this.height+"");
    		properties.setProperty(prefix+"phi",                 this.phi+"");
    		properties.setProperty(prefix+"theta",               this.theta+"");
    		properties.setProperty(prefix+"psi",                 this.psi+"");
			properties.setProperty(prefix+"focalLength",         this.focalLength+"");
			properties.setProperty(prefix+"pixelSize",           this.pixelSize+"");
			properties.setProperty(prefix+"distortionRadius",    this.distortionRadius+"");
			properties.setProperty(prefix+"distortionA8",        this.distortionA8+"");
			properties.setProperty(prefix+"distortionA7",        this.distortionA7+"");
			properties.setProperty(prefix+"distortionA6",        this.distortionA6+"");
			properties.setProperty(prefix+"distortionA5",        this.distortionA5+"");
			properties.setProperty(prefix+"distortionA",         this.distortionA+"");
			properties.setProperty(prefix+"distortionB",         this.distortionB+"");
			properties.setProperty(prefix+"distortionC",         this.distortionC+"");
			properties.setProperty(prefix+"px0",                 this.px0+"");
			properties.setProperty(prefix+"py0",                 this.py0+"");
			for (int i=0;i<this.r_xy.length;i++){
				properties.setProperty(prefix+"r_xy_"+i+"_x",    this.r_xy[i][0]+"");
				properties.setProperty(prefix+"r_xy_"+i+"_y",    this.r_xy[i][1]+"");
			}
			for (int i=0;i<this.r_od.length;i++){
				properties.setProperty(prefix+"r_od_"+i+"_o",    this.r_od[i][0]+"");
				properties.setProperty(prefix+"r_od_"+i+"_d",    this.r_od[i][1]+"");
			}
			properties.setProperty(prefix+"channelWeightDefault",this.channelWeightDefault+"");
			properties.setProperty(prefix+"subcamera",           this.subcamera+"");
			properties.setProperty(prefix+"sensor_port",         this.sensor_port+"");
			properties.setProperty(prefix+"subchannel",          this.subchannel+"");


    	}
    	public void getProperties(String prefix,Properties properties){
    		getProperties(prefix,properties, -1);
    	}
    	public void getProperties(String prefix,Properties properties, int channel){
    		if (properties.getProperty(prefix+"sensorWidth")!=null)
    			this.sensorWidth=Integer.parseInt(properties.getProperty(prefix+"sensorWidth"));
    		if (properties.getProperty(prefix+"sensorHeight")!=null)
    			this.sensorHeight=Integer.parseInt(properties.getProperty(prefix+"sensorHeight"));
       		if (properties.getProperty(prefix+"decimateMasks")!=null)
       			this.decimateMasks=Integer.parseInt(properties.getProperty(prefix+"decimateMasks"));
    		if (properties.getProperty(prefix+"cartesian")!=null)
    			this.cartesian=Boolean.parseBoolean(properties.getProperty(prefix+"cartesian"));
    		if (properties.getProperty(prefix+"lensDistortionModel")!=null)
    			this.lensDistortionModel=Integer.parseInt(properties.getProperty(prefix+"lensDistortionModel"));
    		if (properties.getProperty(prefix+"right")!=null)
    			this.right=Double.parseDouble(properties.getProperty(prefix+"right"));
    		if (properties.getProperty(prefix+"forward")!=null)
    			this.forward=Double.parseDouble(properties.getProperty(prefix+"forward"));
    		if (properties.getProperty(prefix+"heading")!=null)
    			this.heading=Double.parseDouble(properties.getProperty(prefix+"heading"));

    		if (properties.getProperty(prefix+"azimuth")!=null)
    			this.azimuth=Double.parseDouble(properties.getProperty(prefix+"azimuth"));
    		if (properties.getProperty(prefix+"radius")!=null)
    			this.radius=Double.parseDouble(properties.getProperty(prefix+"radius"));
    		if (properties.getProperty(prefix+"height")!=null)
    			this.height=Double.parseDouble(properties.getProperty(prefix+"height"));
    		if (properties.getProperty(prefix+"phi")!=null)
    			this.phi=Double.parseDouble(properties.getProperty(prefix+"phi"));
    		if (properties.getProperty(prefix+"theta")!=null)
    			this.theta=Double.parseDouble(properties.getProperty(prefix+"theta"));
    		if (properties.getProperty(prefix+"psi")!=null)
    			this.psi=Double.parseDouble(properties.getProperty(prefix+"psi"));
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
			if (properties.getProperty(prefix+"px0")!=null)
				this.px0=Double.parseDouble(properties.getProperty(prefix+"px0"));
			if (properties.getProperty(prefix+"py0")!=null)
				this.py0=Double.parseDouble(properties.getProperty(prefix+"py0"));
			setDefaultNonRadial();
			for (int i=0;i<this.r_xy.length;i++){
				if (properties.getProperty(prefix+"r_xy_"+i+"_x")!=null) this.r_xy[i][0]=Double.parseDouble(properties.getProperty(prefix+"r_xy_"+i+"_x"));
				if (properties.getProperty(prefix+"r_xy_"+i+"_y")!=null) this.r_xy[i][1]=Double.parseDouble(properties.getProperty(prefix+"r_xy_"+i+"_y"));
			}
			for (int i=0;i<this.r_od.length;i++){
				if (properties.getProperty(prefix+"r_od_"+i+"_o")!=null) this.r_od[i][0]=Double.parseDouble(properties.getProperty(prefix+"r_od_"+i+"_o"));
				if (properties.getProperty(prefix+"r_od_"+i+"_d")!=null) this.r_od[i][1]=Double.parseDouble(properties.getProperty(prefix+"r_od_"+i+"_d"));
			}

			if (properties.getProperty(prefix+"channelWeightDefault")!=null) {
				this.channelWeightDefault=Double.parseDouble(properties.getProperty(prefix+"channelWeightDefault"));
				this.channelWeightCurrent=this.channelWeightDefault;
			}
//			/enableNoLaser
			if (properties.getProperty(prefix+"enableNoLaser")!=null) {
				this.enableNoLaser=Boolean.parseBoolean(properties.getProperty(prefix+"enableNoLaser"));
			} else {
				this.enableNoLaser=(channel<24);
			}
			if (properties.getProperty(prefix+"subcamera")!=null)
				this.subcamera=Integer.parseInt(properties.getProperty(prefix+"subcamera"));
			if (properties.getProperty(prefix+"sensor_port")!=null)
				this.sensor_port=Integer.parseInt(properties.getProperty(prefix+"sensor_port"));
			if (properties.getProperty(prefix+"subchannel")!=null)
				this.subchannel=Integer.parseInt(properties.getProperty(prefix+"subchannel"));

    	}
    	public void setChannelWeightCurrent(
    			double weight){
    		this.channelWeightCurrent=weight;
    	}
    	public double getChannelWeightCurrent(){
    		return this.channelWeightCurrent;
    	}
    	public double getChannelWeightDefault(){
    		return this.channelWeightDefault;
    	}

    	public void updateCartesian(){ // set alternative parameters
    		if (cartesian) {
    			this.azimuth = Math.atan2(this.right, this.forward)*180.0/Math.PI;
    			this.radius = Math.sqrt(this.forward*this.forward + this.right*this.right);
    			this.phi = this.heading - this.azimuth;
    		} else {
    			this.forward = this.radius * Math.cos(Math.PI*this.azimuth/180.0);
    			this.right =   this.radius * Math.sin(Math.PI*this.azimuth/180.0);
    			this.heading = this.phi +  this.azimuth;
    		}
    	}
    	public void setCartesian(boolean cartesian){
    		if (this.cartesian == cartesian) return; // already in this mode
    		updateCartesian();
    		this.cartesian = cartesian;
    	}
    }
