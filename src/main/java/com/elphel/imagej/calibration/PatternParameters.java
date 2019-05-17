package com.elphel.imagej.calibration;
/*
 **
 ** PatternParameters.java
 **
 ** Copyright (C) 2011-2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  PatternParameters.java is free software: you can redistribute it and/or modify
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
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.Opener;

import java.awt.Rectangle;
import java.util.Properties;

import com.elphel.imagej.dp.CalibrationFileManagement;
import com.elphel.imagej.dp.CalibrationFileManagement.MultipleExtensionsFileFilter;
import com.elphel.imagej.jp4.JP46_Reader_camera;
	/* gridGeometry:
	 * [v][u][0] - x(mm) of the node(u,v), right (looking to the wall) - positive
	 * [v][u][1] - y(mm) of the node(u,v), down - positive	
	 * [v][u][2] - z(mm) of the node(u,v), away (into the wall) - positive	
	 * [v][u][3] - mask 0.0 - outside of the pattern rectangle, 1.0 - on the pattern
	 * RGB moved to Photometric
	 * [v][u][4] - R-intensity (normalized to 1.0 full scale)
	 * [v][u][5] - G-intensity
	 * [v][u][6] - B-intensity
	 * ... repeat 3,4,5,6 for different views of the pattern (from different points)
	 */
	public  class PatternParameters {
		public double patternWidth;  // pattern full width in mm
		public double patternHeight; // pattern full height in mm
		public double patternHalfPeriod; // distance between opposite sign nodes
		public double patternTilt;       // pattern tilt (degrees) - U clockwise from X-right (V clockwise from Y-down)
		public int [] viewMap=null; // should have
		public int    numStations=1;
		public int    numberOfViews=1;
        public double [][][] gridGeometry=null; // [v][u]{x,y,z,{"alpha", R-intensity, G--intensity, B--intensity} alpha=0 - no ghrid, 1 - grid}
        public double [][][] stationZCorr=null; // [v][u]{station0... stationN} - per-station addition to the gridGeometry[][][2] 
        public double [] averageRGB={1.0,1.0,1.0}; 
    	public int numGeometricChannels=4;   // x,y,z,mask
    	public int numPhotometricChannels=4; // r,g,b,a
    	public int defaultNumberOfChannels=26;

//        public int numLayers=7; // number of layers in gridGeometry
        public int U0=0;  //39
        public int V0=0;  //36
        public int debugLevel=2;
        public boolean updateStatus=true;
        public String pathName=null; // path from which the grid was read/ where it was last saved to
        public double [][][][] photometricByView=null; // [numStation][numView]{a,r,g,b}[pixel]
//        public double [][][] photometricBySensor=null;
        public double [] patternErrors=null; // weighted RMS error for each pattern node, calculated
        public double [] patternErrorMask=null; // weighted RMS error for each pattern node, calculated
        public double [] patternErrorMaskSaved=null; // weighted RMS error for each pattern node, calculated
        
        public void initStationZCorr(){
        	this.stationZCorr=new double [this.gridGeometry.length][this.gridGeometry[0].length][this.numStations];
        	for (int v=0;v<this.stationZCorr.length;v++) for (int u=0;u<this.stationZCorr[0].length;u++) for (int s=0;s<this.stationZCorr[0][0].length;s++){
        		this.stationZCorr[v][u][s]=0.0;
        	}
        }
        public void resetStationZCorr(){
        	this.stationZCorr=null;
        }
        public Rectangle getUVDimensions(){
        	return new Rectangle(U0,V0,gridGeometry[0].length,gridGeometry.length);
        }
        public void setPatternErrors(double [] patternErrors){
        	this.patternErrors=patternErrors;
        }
        public double [] getPatternErrors(){
        	return this.patternErrors;
        }
        public double [] getPatternErrorMask(){
        	return this.patternErrorMask;
        }
        public double [] calculatePatternErrorMask(double maxRMS, double minRMS){
        	double [] mask=new double [this.patternErrors.length]; // null pointer
        	double a=1.0/(maxRMS-minRMS);
        	boolean binary=maxRMS==minRMS;
        	
        	for (int i=0;i<this.patternErrors.length;i++){
        		if (Double.isNaN(this.patternErrors[i])) mask[i]=0.0;
        		else {
        			if (binary) mask[i]=(this.patternErrors[i]>=maxRMS)?0.0:1.0;
        			else {
        				double x= a*(this.patternErrors[i]-minRMS);
        				if (x<=0) mask[i]=1.0;
        				else if (x>=1.0) mask[i]=0.0;
        				else mask[i]=1.0+2*x*x*x-3*x*x;
        			}
        		}
        	}
        	this.patternErrorMask=mask;
        	return mask;
        }
        public double [] expandPatternErrorMask(){
        	Rectangle dimensions=getUVDimensions();
        	int [] dirs={1,dimensions.width+1,dimensions.width,dimensions.width-1,-1,-dimensions.width-1,-dimensions.width,dimensions.width+1};
        	double [] mask=this.patternErrorMask.clone();
        	for (int v=1;v<(dimensions.height-1);v++) for (int u=1;u<(dimensions.width-1);u++){
        		int index=v*dimensions.width+u;
        		double min=this.patternErrorMask[index];
        		for (int iDir=0;iDir<dirs.length;iDir++){
        			if (this.patternErrorMask[index+dirs[iDir]]<min) min=this.patternErrorMask[index+dirs[iDir]];
        		}
        		mask[index]=min;
        	}
        	this.patternErrorMask=mask;
        	return this.patternErrorMask;
        }
        //patternErrorMaskSaved
        public void resetPatternErrorMask(){
        	this.patternErrorMask=null;
        }
        public void savePatternErrorMask(){
        	if (this.patternErrorMask!=null) this.patternErrorMaskSaved=this.patternErrorMask.clone();
        	else this.patternErrorMaskSaved=null;
        }
        public void restorePatternErrorMask(){
        	if (this.patternErrorMaskSaved!=null) this.patternErrorMask=this.patternErrorMaskSaved.clone();
        	else this.patternErrorMask=null;
        }
        public double [] getSavedPatternErrorMask(){
        	return this.patternErrorMaskSaved;
        }

		public PatternParameters(
				int [] viewMap,
				int numStations,
				double patternWidth,  // pattern full width in mm
				double patternHeight, // pattern full height in mm
				double patternHalfPeriod, // distance between opposite sign nodes
				double patternTilt        // pattern tilt (degrees) - U clockwise from X-right (V clockwise from Y-down)
		){
			this.numStations=numStations;
			this.patternWidth=patternWidth;
			this.patternHeight=patternHeight;
			this.patternHalfPeriod=patternHalfPeriod;
			this.patternTilt=patternTilt;
			setPhotometric(viewMap);
			calculateGridGeometryAndPhotometric(true);
		}
		public void updateNumStations(int numStations){
			this.numStations=numStations;
			if (numStations==this.photometricByView.length) {
				return;
			}
			System.out.println ("Updating pattern number of stations from "+this.photometricByView.length+" to "+numStations);
			double [][][][] photometricByViewCopy=this.photometricByView;
			this.photometricByView=new double [this.numStations][][][];
			for (int i=0;i<this.photometricByView.length;i++){
				if (i<photometricByViewCopy.length){
					this.photometricByView[i]= photometricByViewCopy[i];
				} else { // deep clone last station
					int iSource=photometricByViewCopy.length-1;
					this.photometricByView[i]=new double [photometricByViewCopy[iSource].length][][];
					for (int j=0;j<photometricByViewCopy[iSource].length;j++){
						if (photometricByViewCopy[iSource][j]==null) {
							this.photometricByView[i][j]=null;
						} else {
							this.photometricByView[i][j]=new double [photometricByViewCopy[iSource][j].length][];
							for (int k=0;k<photometricByViewCopy[iSource][j].length;k++){
								if (photometricByViewCopy[iSource][j][k]==null) {
									this.photometricByView[i][j][k]=null;
								} else {
									this.photometricByView[i][j][k]=photometricByViewCopy[iSource][j][k].clone();
								}								
							}
						}
					}
				}
			}
		}
		public void setPhotometric(int [] viewMap){
			if (viewMap==null){
				this.viewMap=null;
			} else{
				this.viewMap=viewMap.clone();
			}
			setPhotometric();
		}
		public void setPhotometric(){
			if (this.viewMap==null){
				this.viewMap=new int[1];
				this.viewMap[0]=0;
			}
			int maxView=0;
			for (int i=0;i<this.viewMap.length;i++) if (this.viewMap[i]>maxView) maxView=this.viewMap[i];
			this.numberOfViews=maxView+1;
			this.photometricByView=new double [this.numStations][this.numberOfViews][4][];
		}
		public void initDefaultChannels(int num){
			this.viewMap=new int[num];
			for (int i=0;i<num;i++) this.viewMap[i]=(i<24)?0:1;
			this.numberOfViews=2;
		}
		
		public double [][] getPhotometricByView(int stationNumber, int nView){
			if (stationNumber>=this.photometricByView.length) stationNumber=this.photometricByView.length-1;
			if (nView>=this.photometricByView[stationNumber].length) nView=this.photometricByView[stationNumber].length-1;
			return this.photometricByView[stationNumber][nView];
		}
		public double [][] getPhotometricBySensor(int stationNumber,int nSensor){
        	if (getNumStations()<=stationNumber) updateNumStations(stationNumber+1);
//			if (stationNumber>=this.photometricByView.length) stationNumber=this.photometricByView.length-1;
			int isens=nSensor;
			if (nSensor>=this.viewMap.length){
				System.out.println("nSensor="+nSensor+" this.viewMap.length="+this.viewMap.length+" this.photometricByView.length="+this.photometricByView.length);
				nSensor=this.viewMap.length-1;
			}
			int nView=this.viewMap[nSensor];
			if (nView>=this.photometricByView[stationNumber].length){
				System.out.println("nSensor was "+isens+", nView="+nView+" this.photometricByView["+stationNumber+"].length="+this.photometricByView[stationNumber].length);
				nView=this.photometricByView.length-1;
			}
			return this.photometricByView[stationNumber][nView];
		}
//		public int getNumViews(){return this.photometricByView[0].length;}
		public int getNumViews(){return this.numberOfViews;}
		public int getNumStations(){return this.numStations;}
		public void setNumStations(int numStations){this.numStations=numStations;}
		public int [] getViewMap(){return this.viewMap;}
		
		public int getNumGeometricChannels(){return this.numGeometricChannels;}
		public int getNumPhotometricChannels(){return this.numPhotometricChannels;}

		public PatternParameters(
				int [] viewMap,
				double patternWidth,  // pattern full width in mm
				double patternHeight, // pattern full height in mm
				double patternHalfPeriod, // distance between opposite sign nodes
				double patternTilt,        // pattern tilt (degrees) - U clockwise from X-right (V clockwise from Y-down)
				double [] averageRGB
		){
			this.patternWidth=patternWidth;
			this.patternHeight=patternHeight;
			this.patternHalfPeriod=patternHalfPeriod;
			this.patternTilt=patternTilt;
			this.averageRGB=averageRGB.clone();
			setPhotometric(viewMap);
			calculateGridGeometryAndPhotometric(true);
		}
		public PatternParameters clone() {
			PatternParameters patternParameters= new PatternParameters(
					this.viewMap,
					this.patternWidth,  // pattern full width in mm
					this.patternHeight, // pattern full height in mm
					this.patternHalfPeriod, // distance between opposite sign nodes
					this.patternTilt,        // pattern tilt (degrees) - U clockwise from X-right (V clockwise from Y-down)
					this.averageRGB);
			patternParameters.debugLevel=this.debugLevel;
			patternParameters.updateStatus=this.updateStatus;
			patternParameters.pathName=this.pathName;
			return patternParameters;
		}
        public ImagePlus saveGridAsImageStack(String title, String path){
        	ImagePlus imp=getGridAsImageStack(title);
        	if (imp==null) return null;
	   				FileSaver fs=new FileSaver(imp);
	   				if (updateStatus) IJ.showStatus("Saving grid "+path);
	   				if (imp.getStackSize()>1)
	   					fs.saveAsTiffStack(path);
	   				else
	   					fs.saveAsTiff(path);
	   				
	   				this.pathName=path;
	   				if (this.debugLevel>0) System.out.println("Pattern saved as "+this.pathName);
        	return imp;
        }
        
        public String selectAndSave(boolean smart, String defaultPath){
			String [] extensions={".grid-tiff","-grid.tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"Pattern grid *.grid-tiff files");
			String pathname=CalibrationFileManagement.selectFile(
					smart,
					true,
					"Save Pattern Grid Geometry",
					"Save",
					parFilter,
					(defaultPath==null)?this.pathName:defaultPath); //String defaultPath
			if ((pathname!=null)&& (pathname.length()!=0)) saveGridAsImageStack("Pattern Grid", pathname);
			return pathname;
        }

        
        public String selectAndRestore(boolean smart, String defaultPath, int numStations){
			String [] extensions={".grid-tiff","-grid.tiff"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"Pattern grid *.grid-tiff files");
			String pathname=CalibrationFileManagement.selectFile(
					smart,
					false,
					"Restore Pattern Grid Geometry",
					"Restore",
					parFilter,
					(defaultPath==null)?this.pathName:defaultPath); //String defaultPath
			if ((pathname==null) || (pathname=="")) return null;
			setNumStations(numStations);
			setGridFromImageStack(pathname);
			return pathname;
        }
        public void setGridFromImageStack(String path){
    		Opener opener=new Opener();
			ImagePlus imp=opener.openImage("", path);
        	if (imp==null) {
        		String msg="Failed to read grid geometry file "+path;
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	(new JP46_Reader_camera(false)).decodeProperiesFromInfo(imp);
        	setGridFromImageStack(imp);
        	this.pathName=path;
			if (this.debugLevel>0) System.out.println("Opened "+path+" as a stack of the pattern grid geometry");
        }
        /**
         * Loads grid geomnetry (data in mm, X - left, Y - down, Z - into the wall) frome mult-slice
         * ImagePlus. Properties should be set or decoded from the info in the tiff file
         * @param imp - ImagePlus stack, containg x,y,z,alpha slices
         */
        public void setGridFromImageStack(ImagePlus imp){
        	int indexAlpha=3;
        	if (imp == null){
        		String msg="Grid image is null";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (imp.getStackSize()<4){
        		String msg="Expected >=4 slice image with grid geometry";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (	(imp.getProperty("U0") == null) ||
        			(imp.getProperty("V0") == null)){
        		String msg="Properties \"U0\" and/or \"V0\" do not exist";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	this.U0=Integer.parseInt((String) imp.getProperty("U0")); 
        	this.V0=Integer.parseInt((String) imp.getProperty("V0"));
        	// Other properties that are needed only during pattern generation:
        	if (imp.getProperty("patternWidth")!=null) this.patternWidth=Double.parseDouble((String) imp.getProperty("patternWidth"));
        	if (imp.getProperty("patternHeight")!=null) this.patternHeight=Double.parseDouble((String) imp.getProperty("patternHeight"));
        	if (imp.getProperty("patternHalfPeriod")!=null) this.patternHalfPeriod=Double.parseDouble((String) imp.getProperty("patternHalfPeriod"));
        	if (imp.getProperty("patternTilt")!=null) this.patternTilt=Double.parseDouble((String) imp.getProperty("patternTilt"));
        	
        	if (imp.getProperty("AverageRed")!=null)   this.averageRGB[0]=Double.parseDouble((String) imp.getProperty("AverageRed"));
        	if (imp.getProperty("AverageGreen")!=null) this.averageRGB[1]=Double.parseDouble((String) imp.getProperty("AverageGreen"));
        	if (imp.getProperty("AverageBlue")!=null)  this.averageRGB[2]=Double.parseDouble((String) imp.getProperty("AverageBlue"));
        	int numZCorr=0;
        	if (imp.getProperty("numZCorr")!=null)  numZCorr=Integer.parseInt((String) imp.getProperty("numZCorr"));

 //       	if (imp.getProperty("numStations")!=null)  this.numStations=Integer.parseInt((String) imp.getProperty("numStations"));
        	
        	int width=imp.getWidth();
        	int height=imp.getHeight();
    		ImageStack stack = imp.getStack();
        	if (stack==null) {
        		String msg="Expected a image stack with grid geometry";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
    		float [][] pixels =new float[stack.getSize()][];
        	for (int i=0;i<pixels.length;i++) pixels[i]= (float[]) stack.getPixels(i+1);
        	if (this.debugLevel>3){
        		System.out.println("setGridFromImageStack() width="+width+", height="+height+", pixels[0].length="+pixels[0].length);
        	}
        	this.gridGeometry=new double[height][width][getNumGeometricChannels()]; // x,y,z, alpha
        	
        	
        	boolean geometryMaskOnly=(pixels.length==getNumGeometricChannels());
        	boolean singleNoAlpha=(pixels.length==7);
        	int totalNumViews=(geometryMaskOnly || singleNoAlpha)?1:((pixels.length-getNumGeometricChannels()-numZCorr)/getNumPhotometricChannels());
        	int fileNumStations=totalNumViews/getNumViews(); // keep current number of stations
        	int length=height*width;
        	for (int v=0;v<height;v++) for (int u=0;u<width;u++) for (int n=0;n<getNumGeometricChannels();n++){ // x,y,z, alpha
        		this.gridGeometry[v][u][n]=pixels[n][v*width+u];
        	}
        	
        	if (fileNumStations!=getNumStations()){
            	if (this.debugLevel>0){
            		System.out.println("File has "+totalNumViews+" photometric slices, expected "+(getNumStations()*getNumViews())+
            				" ("+getNumStations()+" stations, "+getNumViews()+" views ), skipping loading photometric data");
            	}
            	/// TODO: Fix me!
            	return; //
        	}
        	if (numZCorr>0) {
        		if (numZCorr==getNumStations()) {
        			if (this.debugLevel>0){
        				System.out.println("Loading zCorr data: "+getNumStations()+" slices");
        			}
        			this.stationZCorr=new double [height][width][numZCorr];
        			for (int v=0;v<height;v++) for (int u=0;u<width;u++) for (int n=0;n<numZCorr;n++){
        				this.stationZCorr[v][u][n]=pixels[n+getNumGeometricChannels()][v*width+u];
        			}
        		} else {
            		System.out.println("File has "+numZCorr+" ZCorr slices, current number of stations is "+getNumStations()+
            				", skipping loading zCorr data (per-station pattern Z-correction from the average Z)");
        			
        		}

        	}
        	
        	if (this.debugLevel>0){
        		System.out.println("Loading photometric data: "+(getNumStations()*getNumViews())+" slices "+
        				" ("+getNumStations()+" stations, "+getNumViews()+" vies )");
        	}
        	
        	this.photometricByView=new double [this.numStations][this.numberOfViews][getNumPhotometricChannels()][length]; // r,g,b,a
        	for (int numStation=0;numStation<this.numStations;numStation++) {
        		int useNumStation=numStation;
        		if ((useNumStation>0) && (useNumStation>=fileNumStations)) useNumStation=fileNumStations-1;
        		for (int nView=0;nView<getNumViews();nView++) {
        			for (int chn=0;chn<getNumPhotometricChannels();chn++){
        				int pixIndex=numGeometricChannels+numZCorr+(nView + getNumViews()*useNumStation)*getNumPhotometricChannels()+chn;
        				if ((pixIndex>=pixels.length) && (chn>2)) pixIndex=indexAlpha; // use mask for non-existent alpha
        				if (pixIndex>=pixels.length) {
        					for (int i=0;i<length;i++) {
        						this.photometricByView[useNumStation][nView][chn][i]=(pixels[indexAlpha][i]>0.5)?this.averageRGB[chn]:0.0; // 0<=chn<=2 here OOB =3
        					}
        				} else {
        					for (int i=0;i<length;i++) {
        						this.photometricByView[useNumStation][nView][chn][i]=pixels[pixIndex][i]; //OOB12
        					}
        				}
        			}
        		}
        	}	
        }
 //getNumPhotometricChannels()       
        public ImagePlus getGridAsImageStack(String title){
        	if (this.gridGeometry==null){
        		String msg="Grid geometry does not exist, nothing to convert";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	int width=this.gridGeometry[0].length;
        	int height=this.gridGeometry.length;
        	int length=height*width;
        	int numZCorr=0;
        	if (this.stationZCorr!=null) {
        		for (int v=0;(v<height) && (numZCorr==0);v++ ) for (int u=0;u<width;u++) if (this.stationZCorr[v][u]!=null){
        			numZCorr=this.stationZCorr[v][u].length;
        			break;
        		}
        	}
//        	float [][]pixels=new float [this.gridGeometry[0][0].length][width*height];
        	float [][] pixels=new float [getNumGeometricChannels()+numZCorr+getNumViews()*getNumStations()*getNumPhotometricChannels()] [length];
        	String [] titles=new String[pixels.length];
        	String [] geometricTitles=  {"X","Y","Z","Mask"};
        	String [] zCorrTitles=new String [numZCorr];
        	for (int i=0;i<numZCorr;i++) zCorrTitles[i]="dZ"+i;
        	String [] photometricTitles={"red","green","blue","alpha"};
        	
        	int index=0;
        	for (int v=0;v<height;v++) for (int u=0;u<width;u++){
        		for (int n=0;n<getNumGeometricChannels();n++) { // should be 4==numGeometricChannels
        			pixels[n][index]=  (float) this.gridGeometry[v][u][n];
        		}
        		for (int n=0;n<numZCorr;n++) { // should be 4==numGeometricChannels
        			pixels[n+getNumGeometricChannels()][index]=  (this.stationZCorr[v][u]!=null)?((float) this.stationZCorr[v][u][n]):0.0f;
        		}
        		index++;
        	}
        	for (int n=0;n<getNumGeometricChannels();n++) titles[n]=geometricTitles[n];
        	for (int numStation=0;numStation<this.numStations;numStation++) {
        		for (int nView=0;nView<getNumViews();nView++) {
        			for (int chn=0;chn<getNumPhotometricChannels();chn++){
        				//        		int pixIndex=getNumGeometricChannels()+chn*this.getNumPhotometricChannels()();
        				int pixIndex=getNumGeometricChannels()+numZCorr+(nView + numStation*getNumViews()) *getNumPhotometricChannels()+chn;
        				for (int i=0;i<length;i++) pixels[pixIndex][i]= (float) this.photometricByView[numStation][nView][chn][i]; //OOB 8 //oob 12
        				titles[pixIndex]=photometricTitles[chn]+nView;
        			}
        		}
        	}
        	for (int i=0;i<numZCorr;i++){
        		titles[getNumGeometricChannels()+i]=zCorrTitles[i];
        	}
        	ImagePlus imp=null;
      		ImageStack stack=new ImageStack(width,height);
       		for (int n=0;n<pixels.length;n++)  stack.addSlice(titles[n],    pixels[n]);
       		imp = new ImagePlus(title, stack);
        	imp.setProperty("patternWidth", ""+this.patternWidth);
        	imp.setProperty("patternHeight", ""+this.patternHeight);
        	imp.setProperty("patternHalfPeriod", ""+this.patternHalfPeriod);
        	imp.setProperty("patternTilt", ""+this.patternTilt);
        	imp.setProperty("U0", ""+this.U0);
        	imp.setProperty("V0", ""+this.V0);
        	imp.setProperty("AverageRed",   ""+this.averageRGB[0]);
        	imp.setProperty("AverageGreen", ""+this.averageRGB[1]);
        	imp.setProperty("AverageBlue",  ""+this.averageRGB[2]);
        	if (numZCorr>0) imp.setProperty("numZCorr",  ""+numZCorr);
        	(new JP46_Reader_camera(false)).encodeProperiesToInfo(imp);
        	imp.getProcessor().resetMinAndMax();
        	return imp;
        }
        
        public void applyGridCorrection(double [][] gridCorr){
        	applyGridCorrection(gridCorr, 1.0);
        }
        /**
         * Apply X,Y,Z correction to the current grid geometry
         * @param gridCorr array [4][width*height] of the corrections to grid geometry (x,y,z,weight).
         * Weight is not used, just for information (1.0 - one image file used with full weight)
         * @param scale scale correction
         */
        public void applyGridCorrection(double [][] gridCorr, double scale){
        	if (this.gridGeometry==null){
        		String msg="Grid geometry does not exist, nothing to apply correction to";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	int height=this.gridGeometry.length;
        	int width=this.gridGeometry[0].length;
        	if ((gridCorr==null) || (gridCorr.length!=4) || (gridCorr[0].length!=(width*height))){
        		String msg="Correction is null or does not match pattern grid geometry";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	for (int v=0;v<height;v++) for (int u=0;u<width;u++) if (this.gridGeometry[v][u][3]>0.0){ // only apply to defined grid
        		int vu=v*width+u;
        		for (int n=0;n<3;n++) this.gridGeometry[v][u][n]+= scale*gridCorr[n][vu];
        	}

        }
        
        public void applyZGridCorrection(
    			double [][] gridZCorr3d,
    			double scale){
        	int height=this.gridGeometry.length;
        	int width=this.gridGeometry[0].length;
        	if ((gridZCorr3d==null) || (gridZCorr3d.length!=this.numStations) || (gridZCorr3d[0].length!=(width*height))){
        		String msg="Correction is null or does not match pattern grid geometry";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	this.stationZCorr=new double [height][width][this.numStations];
        	for (int v=0;v<height;v++) for (int u=0;u<width;u++) {
        		int vu=v*width+u;
        		for (int s=0;s<this.numStations;s++){
        	        		this.stationZCorr[v][u][s]=scale*gridZCorr3d[s][vu];
        		}
        	}
    	}

        
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"patternWidth",this.patternWidth+"");
			properties.setProperty(prefix+"patternHeight",this.patternHeight+"");
			properties.setProperty(prefix+"patternHalfPeriod",this.patternHalfPeriod+"");
			properties.setProperty(prefix+"patternTilt",this.patternTilt+"");
			properties.setProperty(prefix+"averageRGB_0",this.averageRGB[0]+"");
			properties.setProperty(prefix+"averageRGB_1",this.averageRGB[1]+"");
			properties.setProperty(prefix+"averageRGB_2",this.averageRGB[2]+"");
			if (this.viewMap!=null){
				properties.setProperty(prefix+"viewMap_length",this.viewMap.length+"");
				for (int i=0;i<this.viewMap.length;i++) {
					properties.setProperty(prefix+"viewMap_"+i,this.viewMap[i]+"");
				}
			}
//			this.viewMap;
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"patternWidth")!=null)
				this.patternWidth=Double.parseDouble(properties.getProperty(prefix+"patternWidth"));
			if (properties.getProperty(prefix+"patternHeight")!=null)
				this.patternHeight=Double.parseDouble(properties.getProperty(prefix+"patternHeight"));
			if (properties.getProperty(prefix+"patternHalfPeriod")!=null)
				this.patternHalfPeriod=Double.parseDouble(properties.getProperty(prefix+"patternHalfPeriod"));
			if (properties.getProperty(prefix+"patternTilt")!=null)
				this.patternTilt=Double.parseDouble(properties.getProperty(prefix+"patternTilt"));
			if (properties.getProperty(prefix+"averageRGB_0")!=null)
				this.averageRGB[0]=Double.parseDouble(properties.getProperty(prefix+"averageRGB_0"));
			if (properties.getProperty(prefix+"averageRGB_1")!=null)
				this.averageRGB[1]=Double.parseDouble(properties.getProperty(prefix+"averageRGB_1"));
			if (properties.getProperty(prefix+"averageRGB_2")!=null)
				this.averageRGB[2]=Double.parseDouble(properties.getProperty(prefix+"averageRGB_2"));
			if (properties.getProperty(prefix+"viewMap_length")!=null) {
				this.viewMap=new int [Integer.parseInt(properties.getProperty(prefix+"viewMap_length"))];
				for (int i=0;i<this.viewMap.length;i++) {
					this.viewMap[i]=0;
					if (properties.getProperty(prefix+"viewMap_"+i)!=null)
						this.viewMap[i]=Integer.parseInt(properties.getProperty(prefix+"viewMap_"+i));
				}
			}
			setPhotometric();
			calculateGridGeometryAndPhotometric(true); // may need to setup this.photometricBySensor
		}
		public boolean showDialog() {
			GenericDialog gd = new GenericDialog("Initial Wall pattern parameters");
			gd.addNumericField("Pattern full width",                  this.patternWidth, 1,6,"mm"); // pattern full width in mm
			gd.addNumericField("Pattern full height",                 this.patternHeight, 1,6,"mm"); // pattern full width in mm
			gd.addNumericField("Distance between opposite sign nodes",this.patternHalfPeriod, 4,8,"mm"); // istance between opposite sign nodes in mm
			gd.addNumericField("Pattern tilt (clockwise)",            this.patternTilt, 1,5,"degrees"); //  pattern tilt (degrees) - U clockwise from X-right (V clockwise from Y-down)
			gd.addNumericField("Average grid RED   (1.0 for white)",  this.averageRGB[0], 3,5,"x"); //
			gd.addNumericField("Average grid GREEN (1.0 for white)",  this.averageRGB[1], 3,5,"x"); //
			gd.addNumericField("Average grid BLUE  (1.0 for white)",  this.averageRGB[2], 3,5,"x"); //
			
			gd.addNumericField("Number of sensors (>24 - two groups, 0 - do not change)",this.defaultNumberOfChannels,0); //
			
			gd.addMessage("Pressing OK will recalculate grid and clear current grid calibration");
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.patternWidth=         gd.getNextNumber();
			this.patternHeight=        gd.getNextNumber();
			this.patternHalfPeriod=    gd.getNextNumber();
			this.patternTilt=          gd.getNextNumber();
			this.averageRGB[0]=        gd.getNextNumber();
			this.averageRGB[1]=        gd.getNextNumber();
			this.averageRGB[2]=        gd.getNextNumber();
			int numberOfChannels= (int)gd.getNextNumber();
			
			if (numberOfChannels>0){
				initDefaultChannels(numberOfChannels);
//				setPhotometric();
			}
			setPhotometric();
			calculateGridGeometryAndPhotometric(true);
			return true;
		}
/**
* Calculate pattern x,y,z==0 and alpha (1.0 - inside, 0.0 - outside) for the grid
* @param resetAll - if true - reset all grid info, if false - only recalculate mask and reset flat field info, keep distortions
*/
		public void calculateGridGeometryAndPhotometric(boolean resetAll){
			// this.photometricByView should be initialized 
			int indexAlpha=3;
			int indexMask=3;
			double cosA=Math.cos(this.patternTilt/180*Math.PI);
			double sinA=Math.sin(this.patternTilt/180*Math.PI);
			double halfWidth=0.5*this.patternWidth;
			double halfHeight=0.5*this.patternHeight;
			double [][] uv={
					{( halfWidth*cosA+halfHeight*sinA)/this.patternHalfPeriod,
						(-halfWidth*sinA+halfHeight*cosA)/this.patternHalfPeriod},
						{( halfWidth*cosA-halfHeight*sinA)/this.patternHalfPeriod,
							(-halfWidth*sinA-halfHeight*cosA)/this.patternHalfPeriod}};
			double [] maxUV={Math.max(Math.abs(uv[0][0]), Math.abs(uv[1][0])),Math.max(Math.abs(uv[0][1]), Math.abs(uv[1][1]))};
			this.U0=(int)Math.ceil(maxUV[0]);
			this.V0=(int)Math.ceil(maxUV[1]);
			resetAll |= (this.gridGeometry==null);
			double x,y;
			// in any case
			int len=(2*this.U0+1)*(2*this.V0+1);
			for (int station=0;station<getNumStations();station++){
				for (int i=0;i<getNumViews();i++){
					for (int chn=0;chn<this.getNumPhotometricChannels();chn++){
						this.photometricByView[station][i][chn]=new double [len]; // r,g,b,alpha
					}
				}
			}

			if (resetAll) {
				for (int station=0;station<getNumStations();station++){
					for (int i=0;i<getNumViews();i++){
						for (int chn=0;chn<this.getNumPhotometricChannels();chn++){
							this.photometricByView[station][i][chn]=new double [len]; // r,g,b,alpha
							for (int j=0;j<len;j++){
								this.photometricByView[station][i][chn][j]=0.0;
							}
						}
					}
				}
				this.gridGeometry=new double[2*this.V0+1][2*this.U0+1][getNumGeometricChannels()]; // without resetAll all properties should be the same
				for (int v=-this.V0; v<=this.V0;v++) for (int u=-this.U0; u<=this.U0;u++){
					x=(u*cosA-v*sinA)*this.patternHalfPeriod;
					y=(u*sinA+v*cosA)*this.patternHalfPeriod;
					boolean inGrid=((x>=-halfWidth) && (x<=halfWidth) && (y>=-halfHeight) && (y<=halfHeight));
					this.gridGeometry[v+this.V0][u+this.U0][0]=x;
					this.gridGeometry[v+this.V0][u+this.U0][1]=y;
					this.gridGeometry[v+this.V0][u+this.U0][2]=0.0;
					this.gridGeometry[v+this.V0][u+this.U0][3]=inGrid?1.0:0.0; // common for all views - to be combined with per-view?
				}

			}
			// always - copy mask to alpha for each view
			int index=0;
			for (int v=-this.V0; v<=this.V0;v++) for (int u=-this.U0; u<=this.U0;u++){
				for (int station=0;station<getNumStations();station++){
					for (int nView=0;nView<getNumViews();nView++) {
						this.photometricByView[station][nView][indexAlpha][index]=this.gridGeometry[v+this.V0][u+this.U0][indexMask];
						for (int c=0;c<3;c++) {
							this.photometricByView[station][nView][c][index]=averageRGB[c];
						}
					}
				}
				index++;
			}			
		}

		
		public int [] uvIndicesToUV (int u1, int v1){
        	if ((v1<0) || (u1<0) ||
        			(v1 >= this.gridGeometry.length) ||
        			(u1 >= this.gridGeometry[0].length) ||
        			(this.gridGeometry[v1][u1][3]==0)) return null;
			int [] iUV={u1-this.U0, v1-this.V0};
			return iUV;
		}
		public double[] getXYZM(int u, int v, boolean verbose, int station){ // u=0,v=0 - center!
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length) ||
        			(this.gridGeometry[v1][u1][3]==0)) {
        		if ((this.debugLevel>1) && verbose && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))){
            		String msg="Requested (acquired) grid (point u="+u+",v="+v+") is outside of the physical grid ["+
            		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
            		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		  //IJ.showMessage("Error",msg);
        		  System.out.println(msg);
        		}
        		return null;
        	}
        	if (this.stationZCorr==null) return this.gridGeometry[v1][u1];
        	double [] result=this.gridGeometry[v1][u1].clone();
// use lower station if grid file does not have current        	
        	int useStation=(this.stationZCorr[v1][u1].length>station)?station:(this.stationZCorr[v1][u1].length-1);
        	result[2]+=this.stationZCorr[v1][u1][useStation];
			return result;
//			return this.gridGeometry[v1][u1];
		}
		public int getGridIndex(int u, int v){ // u=0,v=0 - center!
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length) ||
        			(this.gridGeometry[v1][u1][3]==0)) {
        		return -1;
        	}
			return u1+this.gridGeometry[0].length*v1;
		}
		
		public double[] getXYZM(int u, int v, int station){ // u=0,v=0 - center!
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((this.debugLevel >1) && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))) {
        		String msg="Requested (acquired) grid (point u="+u+",v="+v+") is outside of the physical grid ["+
        		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
        		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	if (this.stationZCorr==null) return this.gridGeometry[v1][u1];
        	if (station>this.stationZCorr[v1][u1].length) station = this.stationZCorr[v1][u1].length-1;
        	double [] result=this.gridGeometry[v1][u1].clone();
        	result[2]+=this.stationZCorr[v1][u1][station]; // if pattern did not have multi-station?
			return result;
		}
		public double[] getXYZMAverage(int u, int v){ // u=0,v=0 - center!
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((this.debugLevel >1) && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))) {
        		String msg="Requested (acquired) grid (point u="+u+",v="+v+") is outside of the physical grid ["+
        		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
        		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return this.gridGeometry[v1][u1];
		}
		public double getZCorr(int u, int v, int station){ // u=0,v=0 - center!
        	if (this.stationZCorr==null) return 0.0;
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((this.debugLevel >1) && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))) {
        		String msg="Requested (acquired) grid (point u="+u+",v="+v+") is outside of the physical grid ["+
        		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
        		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	 return this.stationZCorr[v1][u1][station];
		}
		public double[] getXYZMAverage(int vu){ // u=0,v=0 - center!
			int width=this.gridGeometry[0].length;
			int u1=vu%width;
			int v1=vu/width;
        	if ((this.debugLevel >1) && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))) {
        		String msg="Requested (acquired) grid (point vu="+vu+") is outside of the physical grid ["+
        		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
        		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	return this.gridGeometry[v1][u1];
		}
		public double getZCorr(int vu, int station){ // u=0,v=0 - center!
        	if (this.stationZCorr==null) return 0.0;
			int width=this.gridGeometry[0].length;
			int u1=vu%width;
			int v1=vu/width;
        	if ((this.debugLevel >1) && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))) {
        		String msg="Requested (acquired) grid (point vu="+vu+") is outside of the physical grid ["+
        		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
        		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
        	 return this.stationZCorr[v1][u1][station];
		}

		
		
		/**
		 * Return grid geometry and photometics: X,Y,Z,mask,R,G,B,Alpha
		 * @param u       signed grid U (0 in the center)
		 * @param v       signed grid V (0 in the center)
		 * @param station station number
		 * @param channel channel (sensor) number
		 * @param verbose report out of grid
		 * @return null if out of grid, otherwise X,Y,Z,mask (binary),R (~0.5..1.2),G,B,alpha (0.0..1.0)
		 */
		public double[] getXYZMP(
				int u,
				int v,
				int station,
				int channel,
				boolean verbose){ // u=0,v=0 - center!
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length) ||
        			(this.gridGeometry[v1][u1][3]==0)) {
        		if ((this.debugLevel>1) && verbose && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))){
            		String msg="Requested (acquired) grid (point u="+u+",v="+v+") is outside of the physical grid ["+
            		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
            		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		  //IJ.showMessage("Error",msg);
        		  System.out.println(msg);
        		}
        		return null;
        	}
        	int index=u1+v1*this.gridGeometry[0].length;
        	if (getNumStations()<=station) updateNumStations(station+1);
        	int nView=this.viewMap[channel];
			if (nView>=this.photometricByView[station].length){  // OOB 1// NUll pointer - need to run F-field first?
				nView=this.photometricByView.length-1;
			}

        	double [] result= { // null
        			this.gridGeometry[v1][u1][0],
        			this.gridGeometry[v1][u1][1],
        			this.gridGeometry[v1][u1][2]+((this.stationZCorr!=null)?this.stationZCorr[v1][u1][station]:0.0), // per-station correction
        			this.gridGeometry[v1][u1][3],
        			this.photometricByView[station][nView][0][index],
        			this.photometricByView[station][nView][1][index],
        			this.photometricByView[station][nView][2][index],
        			this.photometricByView[station][nView][3][index]
        	};
			return result;
		}
		public double[] getXYZMPE(
				int u,
				int v,
				int station,
				int channel,
				boolean verbose){ // u=0,v=0 - center!
			int u1=u+this.U0;
			int v1=v+this.V0;
        	if ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length) ||
        			(this.gridGeometry[v1][u1][3]==0)) {
        		if ((this.debugLevel>1) && verbose && ((v1<0) || (u1<0) || (v1 >= this.gridGeometry.length) || (u1 >= this.gridGeometry[0].length))){
            		String msg="Requested (acquired) grid (point u="+u+",v="+v+") is outside of the physical grid ["+
            		(-this.U0)+"..."+(this.gridGeometry[0].length-this.U0-1)+"]x["+
            		(-this.V0)+"..."+(this.gridGeometry.length-this.V0-1)+"]";
        		  //IJ.showMessage("Error",msg);
        		  System.out.println(msg);
        		}
        		return null;
        	}
        	int index=u1+v1*this.gridGeometry[0].length;
        	if (getNumStations()<=station) updateNumStations(station+1);
        	int nView=this.viewMap[channel];
        	int useStation=(this.stationZCorr!=null)?((this.stationZCorr[v1][u1].length>station)?station:(this.stationZCorr[v1][u1].length-1)):0;
        	if (nView>=this.photometricByView[station].length){  // OOB 1// NUll pointer - need to run F-field first? (oob1 when grid had less than now
				nView=this.photometricByView.length-1;
			}
        	if (index > this.photometricByView[station][nView][0].length){
        		System.out.println("getXYZMPE ("+u+","+v+","+station+","+channel+")");
        		System.out.println("this.photometricByView[station][nView][0].length="+this.photometricByView[station][nView][0].length);
        		System.out.println("index="+index+" vView="+nView);
        		System.out.println();
        	}

        	double [] result= { // null 
        			this.gridGeometry[v1][u1][0],
        			this.gridGeometry[v1][u1][1],
//        			this.gridGeometry[v1][u1][2]+((this.stationZCorr!=null)?this.stationZCorr[v1][u1][station]:0.0), // per-station correction
        			this.gridGeometry[v1][u1][2]+((this.stationZCorr!=null)?this.stationZCorr[v1][u1][useStation]:0.0), // per-station correction
        			this.gridGeometry[v1][u1][3],
        			this.photometricByView[station][nView][0][index], // java.lang.ArrayIndexOutOfBoundsException: 14215
        			this.photometricByView[station][nView][1][index],
        			this.photometricByView[station][nView][2][index],
        			this.photometricByView[station][nView][3][index],
        			(this.patternErrorMask==null)?1.0:this.patternErrorMask[index]
        	};
			return result;
		}
		public double [][][] getGeometry(){return this.gridGeometry;} 
	}
	
	
	
