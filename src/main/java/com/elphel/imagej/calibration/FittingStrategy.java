package com.elphel.imagej.calibration;
/*
 **
 ** FittingStrategy.java
 **
 ** Copyright (C) 2011-2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  FittingStrategy.java is free software: you can redistribute it and/or modify
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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.configuration.ConfigurationException;
import org.apache.commons.configuration.XMLConfiguration;

import com.elphel.imagej.calibration.CalibrationFileManagement.MultipleExtensionsFileFilter;
import com.elphel.imagej.common.WindowTools;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.text.TextWindow;



    /**
     *
     * Specifies images to process and parameters to adjust
     * Each parameter cab be:
     * 0 - "fixed" - use individual, per-image parameters, do not modify them
     * 1 - "common" - parameters are common for all selected images.
     *                When saving - save to all selected images.
     *                When loading (if different) use "master image" (or closest to it in time)
     * 2 - "super common" - same as common, but save to all images, not just selected
     * 3 - "individual"
     * 4 - "per group"
     * 5 - "per-station"
     * 6 - "per-station" save to all (super)
     * 7 - "weak common" - like common, but enable small individual variations (for a price)
     * 8 - "weak station" - like per-station, but enable individual (for a price)
     *
+====================+===========+===========+
|                    |  Same TS  |  Diff TS  |
|                    +-----+-----+-----+-----+
|                    |  C  |  I  |  C  |  I  |
+============+=======+=====+=====+=====+=====+
| Same       |Eyesis |  X  |  X  |  C  |  I  |
| Subcamera  +-------+-----+-----+-----+-----+
|            |Subcam |  X  |  X  |  C  |  I* |
+============+=======+=====+=====+=====+=====+
| Different  |Eyesis |  C  |  C  |  C  |  I  |
| Subcameras +-------+-----+-----+-----+-----+
|            |Subcam |  I  |  I  |  I  |  I  |
+============+=======+=====+=====+=====+=====+
I* - special case when the subcamera is being adjusted/replaced. How to deal with it?

     *
     */
    public  class FittingStrategy{
    	public String pathName=null; // path to XML file this instance was created from
    	public DistortionCalibrationData distortionCalibrationData=null;// per-image parameters
    	private boolean [][] selectedImages=null; // images selected for each step (will be masked with enabled images)
    	private boolean [][] selectedValidImages=null; // images selected for each step same as selected, but only if number of weight>0.0 nodes > threshold
    	public boolean [] stopAfterThis=null;
    	final public int modeFixed=0,modeCommon=1,modeSupercommon=2,modeIndividual=3,modeGroup=4;
    	final public int modeStation=5,modeSuperStation=6,modeWeakCommon=7,modeWeakStation=8,modeTiltEqualize=9;
    	public String[] definedModes={
    			"fixed",                         // modeFixed=0
    			"common",                        // modeCommon=1
    			"common, save to all",           // modeSupercommon=2
    			"individual",                    // modeIndividual=3
    			"per-group",                     // modeGroup=4
    			"per-station",                   // modeStation=5
    			"per-station, save to  all",     // modeSuperStation=6
    			"weak common, save to  all",     // modeWeakCommon=7
    			"weak per-station, save to  all" // modeWeakStation=8
    			};
    	public String[] definedModesNoWeak={     // show for parameters that can not use weak
    			"fixed",                         // modeFixed=0
    			"common",                        // modeCommon=1
    			"common, save to all",           // modeSupercommon=2
    			"individual",                    // modeIndividual=3
    			"per-group",                     // modeGroup=4
    			"per-station",                   // modeStation=5
    			"per-station, save to  all"      // modeSuperStation=6
//    			"weak common",              // modeWeakCommon=7
//    			"weak per-station",         // modeWeakStation=8
    			};
    	public String[] definedModesTiltEq={      // show for parameters that can not use weak
    			"fixed",                         // modeFixed=0
    			"common",                        // modeCommon=1
    			"common, save to all",           // modeSupercommon=2
    			"individual",                    // modeIndividual=3
    			"per-group",                     // modeGroup=4
    			"per-station",                   // modeStation=5
    			"per-station, save to  all",     // modeSuperStation=6
    			"weak common, save to  all",     // modeWeakCommon=7
    			"weak per-station, save to  all",// modeWeakStation=8
    			"tilt equalize"                 // modeTiltEqualize
    			};
    	public String[] definedModesAll= definedModesTiltEq;
    	public int [][] parameterMode=null; // per series, per-parameter
    	public int [][][] parameterGroups=null; // per series, per-parameter - null or array of group numbers (1 element per image)
    	public int [][] zGroups=null;
    	public boolean saveUnusedGroups=false; // purge groups for parameters when saving to XML, preserve if true
    	public double [] lambdas=null;   // LMA initial lambda for each step
    	public double defaultLambda=0.001;
    	public int [][] parameterList=null; // list of all parameters in the system, each has subcamera number and parameters index
    	public boolean [] parameterEnable=null; // select which parameters (of some 320 to display in selector)
    	public int [] masterImages=null;        // number of image no take "common" parameters from if they are different
    	public double defaultStepDone=1.0E-6;
    	public double [] stepDone=null;// delta_error to error ratio to consider series finished
//    	public double [] parameterVector;
    	public int [][] parameterMap=null;// ** valid for one strategy sequence *** for each parameter vector element hold image number/parameter number
    	public int [][] reverseParameterMap=null; // reversed map - for each [imageNumber][parameterNumber] -> vector index (-1 fixed)
    	public int currentSeriesNumber=-1;     // currently selected strategy step (for which parameterMap and reverseParameterMap)
    	public int debugLevel=2;
    	public int [] varianceModes=null;      // per-series: 0 - disabled, 1 - constant imageSets weight, 2 - variable imageSets weight
    	final public int varianceModeDisabled=0,varianceModeSameWeight=1,varianceModeVariableWeight=2;
    	public String [] definedVarianceModes={
    			"disabled",
    			"same weight",
    			"variable weight",
    	};

    	// next arrays will be initialized at buildVariancesMaps only if at least some parameters use variances, otherwise they will be null
    	public double [] variationsAverages=null; // holds per extrinsic parameter or per parameter/per station average values
    	public int    [] averageCellIndex=  null; // for each element in the parameters vector holds index in  variationsAverages array
    	public double [] weightOnAverage=   null; // how variation of this parameter influences average (for all or for currenrt station
    	public double [] weightVariance=    null; // weight for LMA over variance of parameters from average
    	public double [] varianceErrorsSquared=  null; // weighted squared errors

    	public String [] strategyComment = null;

    	public FittingStrategy(
    	    	DistortionCalibrationData distortionCalibrationData// per-image parameters
    			) {
    		this.distortionCalibrationData=distortionCalibrationData;
    		setDflt(0);
    	}
		public FittingStrategy(
    	    	DistortionCalibrationData distortionCalibrationData,// per-image parameters
    	    	int numSeries  // number of iteration series
    			){
    		this.distortionCalibrationData=distortionCalibrationData;
    		this.selectedImages=new boolean[numSeries][this.distortionCalibrationData.getNumImages()];
    		this.selectedValidImages=new boolean[numSeries][this.distortionCalibrationData.getNumImages()];

    		setDflt(numSeries);
    	}

	   	public FittingStrategy(
        		boolean smart,
        		String defaultPath,
    	    	DistortionCalibrationData distortionCalibrationData // per-image parameters
    			) {
			String [] extensions={".stg-xml","-strategy.xml"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"*.stg-xml files");
			String pathname=CalibrationFileManagement.selectFile(
					smart,
					false,
					"Restore Fitting Strategy",
					"Restore",
					parFilter,
					defaultPath); //String defaultPath
			if ((pathname==null) || (pathname=="")) return;
	   		setFromXML(
	    	    	distortionCalibrationData, // per-image parameters
	    	    	pathname);
			System.out.println("Restored fitting strategy from "+pathname);
	   	}


    	/**
    	 * Reads FittingStrategy from an XML file
    	 * @param distortionCalibrationData should be defined before this class!
    	 * @param pathname pathe to the saved data
    	 */
	   	public FittingStrategy(
    	    	DistortionCalibrationData distortionCalibrationData, // per-image parameters
    	    	String pathname
    			) {
	   		setFromXML(
	    	    	distortionCalibrationData, // per-image parameters
	    	    	pathname);
	   	}
    	public void setFromXML(
    	    	DistortionCalibrationData distortionCalibrationData, // per-image parameters
    	    	String pathname
    			) {
    		this.distortionCalibrationData=distortionCalibrationData;
        	XMLConfiguration hConfig=null;
        	try {
				hConfig=new XMLConfiguration(pathname);
				hConfig.setListDelimiter((char) 0);
				hConfig.setDelimiterParsingDisabled(true);
				hConfig.setAttributeSplittingDisabled(true);

			} catch (ConfigurationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			hConfig.setThrowExceptionOnMissing(false); // default value, will return null on missing
// read parameterList
			int len=Integer.parseInt(hConfig.getString("parameterMap.length","0"));
			if (len<=0) {
				String msg="No parameterMap specified in "+ pathname;
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
			}
	    	this.parameterList=new int [len][2];
	    	this.parameterEnable=new boolean[len];
    		for (int i=0; i<this.parameterList.length; i++){
    			this.parameterList[i][0]=Integer.parseInt(hConfig.getString("parameterMap.par_"+i+".subcamera"));
    			this.parameterList[i][1]=Integer.parseInt(hConfig.getString("parameterMap.par_"+i+".index"));
    			this.parameterEnable[i]= (Integer.parseInt(hConfig.getString("parameterMap.par_"+i+".visible"))>0);
    		}
			int nSer=Integer.parseInt(hConfig.getString("series.number","0"));
			if (nSer==0) return; // arrays will just be null
			this.selectedImages=new boolean[nSer][];
			this.selectedValidImages=new boolean[nSer][];
			this.masterImages=new int[nSer];
			this.stopAfterThis=new boolean[nSer];
			for (int i=0;i<nSer;i++) this.stopAfterThis[i]=true; // older configuration files did not have it
			this.parameterMode=new int[nSer][len];
			this.varianceModes=new int[nSer];
			this.zGroups=new int [nSer][];
			this.parameterGroups=new int[nSer][len][];
			this.lambdas=new double  [nSer];
			this.stepDone=new double [nSer];
			this.strategyComment=new String[nSer];
        	for (int i=0;i<nSer;i++){ // iterate through series
//this.stopAfterThis
// read selected images (no check here that it matches to the    distortionCalibrationData!
        		String sSeries="series.series_"+i;
        		String fs=hConfig.getString(sSeries+".selectedImages");
        		this.selectedImages[i]=new boolean[fs.length()];
        		this.selectedValidImages[i]=null;
        		for (int j=0;j<this.selectedImages[i].length;j++) this.selectedImages[i][j]=(fs.charAt(j)=='+');
    			this.masterImages[i]=Integer.parseInt(hConfig.getString(sSeries+".masterImage"));
    			if (hConfig.getString(sSeries+".varianceModes")!=null) {
    				this.varianceModes[i]=Integer.parseInt(hConfig.getString(sSeries+".varianceModes"));
    			} else {
    				this.varianceModes[i]=varianceModeDisabled;
    			}

        		for (int j=0;j<this.parameterList.length;j++){
        			int nPar=this.parameterList[j][1];
        			int nSub=this.parameterList[j][0];
        			String hconfigName=sSeries+".parameterMode."+
    				this.distortionCalibrationData.parameterDescriptions[nPar][0]+ // here intentionally left parameterDescriptions[][], not descrField
    				// so switching Cartesian on/off will not make make strategy invalid
    				(this.distortionCalibrationData.isSubcameraParameter(nPar)?("_sub"+nSub):"");
//        			System.out.println("Setting this.parameterMode["+i+"]["+j+"] from " +hconfigName);
        			if (hConfig.getString(hconfigName)!=null)  this.parameterMode[i][j]=Integer.parseInt(hConfig.getString(hconfigName));
        			else System.out.println("Failed to set this.parameterMode["+i+"]["+j+"] from " +hconfigName+" - maybe it is a new parameter not present in the file "+pathname);
        			this.parameterGroups[i][j]=null;
//Try to read series
        			for (int ni=0;ni<this.selectedImages.length;ni++){
        				String sGroup=hConfig.getString(sSeries+".parameterMode."+
                				this.distortionCalibrationData.parameterDescriptions[nPar][0]+
                				(this.distortionCalibrationData.isSubcameraParameter(nPar)?("_sub"+nSub):"")+
                				"_group"+ni);
        				if (sGroup!=null){
        					if (this.parameterGroups[i][j]==null){
        						this.parameterGroups[i][j]=new int[this.selectedImages.length];
        						for (int ni1=0;ni1<this.selectedImages.length;ni1++)this.parameterGroups[i][j][ni1]=0;
        					}
        					this.parameterGroups[i][j][ni]=Integer.parseInt(sGroup);
        				}
        			}
        		}
        		this.lambdas[i]=Double.parseDouble(hConfig.getString(sSeries+".lambdas"));
        		this.stepDone[i]=Double.parseDouble(hConfig.getString(sSeries+".stepDone"));
        		if (hConfig.getString(sSeries+".stopAfterThis")!=null){
        			this.stopAfterThis[i]=Boolean.parseBoolean(hConfig.getString(sSeries+".stopAfterThis"));
        		}
    			String zGroupsName=sSeries+".zGroups";
    			if (hConfig.getString(zGroupsName)!=null) {
    				int numZgroups=Integer.parseInt(hConfig.getString(zGroupsName));
    				this.zGroups[i]=new int [numZgroups];
    				for (int nZGroup=0;nZGroup<numZgroups;nZGroup++){
    					this.zGroups[i][nZGroup]=Integer.parseInt(hConfig.getString(sSeries+".zGroup_"+nZGroup));
    				}
    			} else {
    				this.zGroups[i]=null;
    			}

    			if (hConfig.getString(sSeries+".strategyComment")!=null) {
    				this.strategyComment[i]=hConfig.getString(sSeries+".strategyComment");
    				if ((this.strategyComment[i].length()>10) && this.strategyComment[i].substring(0,9).equals("<![CDATA["))
    					this.strategyComment[i]=this.strategyComment[i].substring(9,this.strategyComment[i].length()-3);
    			}
        	}
//        	if (hConfig.getString("variances")!=null){
//            if (!hConfig.configurationAt("variances").isEmpty()){
       		if (hConfig.configurationsAt("variances").size()>0){
//        		System.out.println("hConfig.configurationAt(\"variances\").isEmpty()=false");
        		this.distortionCalibrationData.eyesisCameraParameters.getCostsPropertiesXML("variances.",hConfig);
        	}
        	this.pathName=pathname;
//        	distortionCalibrationData.readAllGrids();
    	}

/**
 * Adjust  selectedImages and parameterGroups after number of images may change, old number of images should be >0
 * If the new number is larger, the first images attributes will be copied, the extra ones - use those of the last of the old ones
 * @param newNumberOfImages new number of images used by this series
 */
    	public void adjustNumberOfImages(int newNumberOfImages){
    		if ((this.selectedImages == null) ||(this.selectedImages[0].length==0)) {
    			String msg="selectedImages array is "+((this.selectedImages == null)?"null":"empty");
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		if (this.selectedImages[0].length==newNumberOfImages) return;
    		for (int nSer=0;nSer<this.selectedImages.length;nSer++){
    			boolean [] tmp=selectedImages[nSer].clone();
    			this.selectedImages[nSer]=new boolean[newNumberOfImages];
    			this.selectedValidImages[nSer]=null; // just invalidate
    			for (int i=0;i<newNumberOfImages;i++) this.selectedImages[nSer][i]=(i<tmp.length)?tmp[i]:tmp[tmp.length-1];
    			for (int nPar=0;nPar<this.parameterList.length;nPar++){
    				if (this.parameterGroups[nSer][nPar]!=null){
    					int [] tmpGroup=this.parameterGroups[nSer][nPar].clone();
    					this.parameterGroups[nSer][nPar]=new int[newNumberOfImages];
    	    			for (int i=0;i<newNumberOfImages;i++) this.parameterGroups[nSer][nPar][i]=(i<tmpGroup.length)?tmpGroup[i]:tmpGroup[tmpGroup.length-1];
    				}
    			}

    		}
    	}

    	public String selectAndSaveToXML(
    			boolean smart,
    			String defaultPath){
			String [] extensions={".stg-xml","-strategy.xml"};
			CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,"*.stg-xml files");
			if ((defaultPath==null) || (defaultPath.length()==0)){
				defaultPath=this.pathName;
			}
			String pathname=CalibrationFileManagement.selectFile(
					smart,
					true,
					"Save Fitting Strategy",
					"Save",
					parFilter,
					defaultPath); //String defaultPath
			if (pathname!=null) saveToXML(pathname);
			return pathname;

    	}
    	  //http://commons.apache.org/configuration/userguide/howto_xml.html
        public boolean saveToXML(String pathname) {
        	if (pathname==null) return false;
        	XMLConfiguration hConfig=new XMLConfiguration();

        	hConfig.setListDelimiter((char) 0);
			hConfig.setDelimiterParsingDisabled(true);
			hConfig.setAttributeSplittingDisabled(true);

        	hConfig.setRootElementName("FittingStrategy");
//        	write which of the overall parameter correspond to each global/subcamera one, and which parameters are hidden/shown
    		hConfig.addProperty("parameterMap","");
    		hConfig.addProperty("parameterMap.length",this.parameterList.length);
    		for (int i=0; i<this.parameterList.length; i++){
        		hConfig.addProperty("parameterMap.par_"+i,"");
        		hConfig.addProperty("parameterMap.par_"+i+".subcamera",this.parameterList[i][0]);
        		hConfig.addProperty("parameterMap.par_"+i+".index",this.parameterList[i][1]);
        		hConfig.addProperty("parameterMap.par_"+i+".visible",this.parameterEnable[i]?"1":"0");
    		}

    		hConfig.addProperty("series","");
    		hConfig.addProperty("series.number",this.selectedImages.length);

        	for (int i=0;i<this.selectedImages.length;i++){ // iterate through series
        		String sSeries="series.series_"+i;
        		hConfig.addProperty(sSeries,"");
        		String fs="";
        		for (int j=0;j<this.selectedImages[i].length;j++)fs+=this.selectedImages[i][j]?"+":"-";
        		hConfig.addProperty(sSeries+".selectedImages",fs);
        		hConfig.addProperty(sSeries+".masterImage",this.masterImages[i]);
        		if (this.varianceModes!=null) hConfig.addProperty(sSeries+".varianceModes",this.varianceModes[i]);
        		hConfig.addProperty(sSeries+".parameterMode","");
        		for (int j=0;j<this.parameterList.length;j++){
        			int nPar=this.parameterList[j][1];
        			int nSub=this.parameterList[j][0];
            		hConfig.addProperty(sSeries+".parameterMode."+
            				this.distortionCalibrationData.parameterDescriptions[nPar][0]+
            				(this.distortionCalibrationData.isSubcameraParameter(nPar)?("_sub"+nSub):""),
            				this.parameterMode[i][j]);
// cleaning up output - removing unused groups (may disable
            		if (((this.parameterMode[i][j]==this.modeGroup) || this.saveUnusedGroups ) &&  (this.parameterGroups[i][j]!=null)){
            			for (int ni=0;ni<this.parameterGroups[i][j].length;ni++){
                    		hConfig.addProperty(sSeries+".parameterMode."+
                    				this.distortionCalibrationData.parameterDescriptions[nPar][0]+
                    				(this.distortionCalibrationData.isSubcameraParameter(nPar)?("_sub"+nSub):"") +"_group"+ni,
                    				this.parameterGroups[i][j][ni]);
            			}
            		}
        		}
        		if (this.zGroups[i]!=null) {
        			hConfig.addProperty(sSeries+".zGroups", this.zGroups[i].length);
    				for (int nZGroup=0;nZGroup<this.zGroups[i].length;nZGroup++){
    					hConfig.addProperty(sSeries+".zGroup_"+nZGroup,this.zGroups[i][nZGroup]);
    				}
        		}
        		hConfig.addProperty(sSeries+".lambdas",lambdas[i]);
        		hConfig.addProperty(sSeries+".stepDone",stepDone[i]);
        		hConfig.addProperty(sSeries+".stopAfterThis",this.stopAfterThis[i]);

            	if (this.strategyComment[i]!=null){
            		String comment_esc=this.strategyComment[i].replace(",","\\,");
//            		comment_esc=comment_esc.replace("/","\\/");
            		//          hConfig.addProperty("comment","<![CDATA["+comment_esc+ "]]>");
            		hConfig.addProperty(sSeries+".strategyComment",comment_esc);
            	}
//        		hConfig.addProperty(sSeries+".strategyComment","<![CDATA["+this.strategyComment[i]+"]]>");
//        		hConfig.addPropertyDirect(sSeries+".strategyComment","<![CDATA["+this.strategyComment[i]+"]]>");
//        		hConfig.setProperty(sSeries+".strategyComment","<![CDATA["+this.strategyComment[i]+"]]>");
        	}
        	hConfig.addProperty("variances","");
       		this.distortionCalibrationData.eyesisCameraParameters.setCostsPropertiesXML("variances.",hConfig);
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
        public int findLastValidSeries(int numSeries){
        	for (;numSeries>=0;numSeries--) if (isSeriesValid(numSeries)) return numSeries;
        	return -1; // none valid
        }
        public boolean copySeries(int sourceSeries, int numSeries){
     		if (isSeriesValid(sourceSeries)){
    				invalidateSelectedImages(numSeries); // just in case -= will be recalculated
    				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++){
    					this.selectedImages[numSeries][i]=this.selectedImages[sourceSeries][i]; // copy for all, not only enabled
    				}
    				for (int i =0; i<this.parameterEnable.length;i++) if (this.parameterEnable[i]){
    					this.parameterMode[numSeries][i]=this.parameterMode[sourceSeries][i];
    					this.parameterGroups[numSeries][i]=(this.parameterGroups[sourceSeries][i]==null)?null:this.parameterGroups[sourceSeries][i].clone();
    				}
    				this.masterImages[numSeries]=this.masterImages[sourceSeries];
    				this.lambdas[numSeries]=this.lambdas[sourceSeries];
    				this.stepDone[numSeries]=this.stepDone[sourceSeries];
    				this.stopAfterThis[numSeries]=this.stopAfterThis[sourceSeries];
    				this.varianceModes[numSeries]=this.varianceModes[sourceSeries];
    				this.zGroups[numSeries]=(this.zGroups[sourceSeries]!=null)?this.zGroups[sourceSeries].clone():null;
    				this.strategyComment[numSeries] = this.strategyComment[sourceSeries]+" (copy)";
    				return true;
    		} else {
    			return false;
    		}
        }



        /**
         * Verifies that series exists and is not empty (valid for LMA)
         * @param num - series to verify
         * @return true if series is valid
         */
    	public boolean isSeriesValidVerbose(int num){
    		if ((num<0) || (num>=this.selectedImages.length)) return false;
    		boolean someImagesSelected=false;
    		for (int i=0;i<this.selectedImages[num].length;i++) if (this.selectedImages[num][i]){
    			someImagesSelected=true;
    			break;
    		}
    		if (!someImagesSelected) {
    			if (this.debugLevel>0) System.out.println("isSeriesValid("+num+"): no images selected");
    			return false;
    		}
    		boolean someParametersSelected=false;
    		for (int i=0;i<this.parameterMode[num].length;i++) if (this.parameterMode[num][i]>0){
    			someParametersSelected=true;
    			break;
    		}
    		if (!someParametersSelected) {
    			if (this.debugLevel>0) System.out.println("isSeriesValid("+num+"): no parameters selected");
    			return false;
    		}
    		return true;
    	}

    	public boolean isSeriesValid(int num){
    		if ((num<0) || (num>=this.selectedImages.length)) return false;
    		boolean someImagesSelected=false;
    		for (int i=0;i<this.selectedImages[num].length;i++) if (this.selectedImages[num][i]){
    			someImagesSelected=true;
    			break;
    		}
    		if (!someImagesSelected) {
    			return false;
    		}
    		boolean someParametersSelected=false;
    		for (int i=0;i<this.parameterMode[num].length;i++) if (this.parameterMode[num][i]>0){
    			someParametersSelected=true;
    			break;
    		}
    		if (!someParametersSelected) {
    			return false;
    		}
    		return true;
    	}
    	/**
    	 * Determines if LMA should stop after this series (either flag is set or next is invalid)
    	 * @param num number of series to test
    	 * @return true if LMA should stop
    	 */
    	public boolean isLastSeries(int num){
    		if (this.stopAfterThis[num]) return true;
    		if (!isSeriesValid(num+1)) return true;
    		return false;
    	}

        /**
         *
         * Specifies images to process and parameters to adjust
         * Each parameter cab be:
         * 0 - "fixed" - use individual, per-image parameters, do not modify them
         * 1 - "common" - parameters are common for all selected images.
         *                When saving - save to all selected images.
         *                When loading (if different) use "master image" (or closest to it in time)
         * 2 - "super common" - same as common, but save to all images, not just selected
         * 3 - "individual"
    +====================+===========+===========+
    |                    |  Same TS  |  Diff TS  |
    |                    +-----+-----+-----+-----+
    |                    |  C  |  I  |  C  |  I  |
    +============+=======+=====+=====+=====+=====+
    | Same       |Eyesis |  X  |  X  |  C  |  I  |
    | Subcamera  +-------+-----+-----+-----+-----+
    |            |Subcam |  X  |  X  |  C  |  I* |
    +============+=======+=====+=====+=====+=====+
    | Different  |Eyesis |  C  |  C  |  C  |  I  |
    | Subcameras +-------+-----+-----+-----+-----+
    |            |Subcam |  X  |  X  |  X  |  X  | subcameras for different subcameras are individual parameters
    +============+=======+=====+=====+=====+=====+
    I* - special case when the subcamera is being adjusted/replaced. How to deal with it?

         *
         */

 //    	public int [][] reverseParameterMap=null; // reversed map - for each [imageNumber][parameterNumber] -> vector index (-1 fixed)\
    	// selects all enabled images in the specified series (modifies strategy series!)
    	/**
    	 * selects all enabled images in the specified series (modifies strategy series!)
    	 * @param ser number of fitting strategy series to set
    	 * @return old value for series images selection
    	 */
    	public boolean [] selectAllImages(int ser){
 //   		this.selectedValidImages[ser]=null; // just invalidate - do nothing, so restore would work
    		boolean [] oldSelection=this.selectedImages[ser].clone();
			for (int i=0;i<this.selectedImages[ser].length;i++) this.selectedImages[ser][i]=true; // select all\
			return oldSelection;
    	}
    	/**
    	 * Sets image selection. Can be used to restore saved selection after selectAllImages(int ser)
    	 * @param ser number of fitting strategy series to set
    	 * @param selection array specifying selected images
    	 */
    	public void setImageSelection(int ser, boolean [] selection){
    		this.selectedImages[ser]=selection.clone();
    	}
    	public boolean [] selectedAllImages() {
//    		if (this.currentSeriesNumber<0) return null;
			return selectedAllImages(this.currentSeriesNumber);
		}
    	public boolean [] selectedAllImages(int ser) {
    		if ((ser<0) || (ser>=this.selectedImages.length)){
    			boolean [] allImages=new boolean[this.selectedImages[0].length];
    			for (int i=0;i<allImages.length;i++) allImages[i]=true; // select all
    			return allImages;
    		} else {
        		return this.selectedImages[ser];
    		}
		}
    	public boolean [] selectedImages() {
//    		if (this.currentSeriesNumber<0) return null;
			return selectedImages(this.currentSeriesNumber);
		}
    	public boolean [] selectedImagesNoBadKernels(int ser) {
    		boolean [] selected=selectedImages(ser);
    		for (int i=0;i<selected.length;i++) selected[i]&= !this.distortionCalibrationData.gIP[i].noUsefulPSFKernels;
    		return selected;
    	}
    	public void setNoUsefulPSFKernels(int i, boolean noUsefulPSFKernels){
    		this.distortionCalibrationData.gIP[i].noUsefulPSFKernels=noUsefulPSFKernels;
    	}
    	public void invalidateSelectedImages(int ser){
    		this.selectedValidImages[ser]=null;
    	}
    	public void initSelectedValidImages(int ser){
    		this.selectedValidImages[ser]=this.selectedImages[ser].clone();
    	}
    	public void invalidateSelectedImage(int ser,int nImg){
    		this.selectedValidImages[ser][nImg]=false;
    	}

    	//this.selectedValidImages
    	public boolean [] selectedImages(int ser) {
    		return selectedImages(ser,false);
    	}
    	public boolean [] selectedImages(int ser, boolean userSelection) {
    		boolean [] selectedMasked;
    		if ((ser<0) || (ser>=this.selectedImages.length)){
    			selectedMasked=new boolean[this.selectedImages[0].length];
    			for (int i=0;i<selectedMasked.length;i++) selectedMasked[i]=true; // select all
    		} else {
    			selectedMasked=this.selectedImages[ser].clone();
    		}
    		for (int i=0;i<selectedMasked.length;i++) selectedMasked[i] &= this.distortionCalibrationData.gIP[i].enabled;
    		// unselect empty
    		// TODO: add minimal number of nodes?
			if ((ser>=0) && (this.selectedValidImages[ser]!=null) && !userSelection){
	    		for (int i=0;i<selectedMasked.length;i++) selectedMasked[i] &= this.selectedValidImages[ser][i];
			} else {
				for (int i=0;i<selectedMasked.length;i++) selectedMasked[i] &= (this.distortionCalibrationData.gIP[i].pixelsXY.length>0);
			}
    		return selectedMasked; //this.selectedImages[ser];
    	}
    	public int getNumSeries(){
    		return this.selectedImages.length;
    	}

       	/**
    	 * Creates map from the parameter vector index to the {grid image number, parameter number}
    	 * When the parameter is shared by several images, the map points to the one which value will be used
    	 * (they might be different). Timestamp of the masterImages[] is used to determine which image to use.
    	 * Simultaneously creates this.reverseParameterMap that maps each of the image/parameter to the parameter vector
    	 * Needs to be run for each new strategy series
    	 * @param numSeries number of fitting strategy series
    	 * @return this.parameterMap
    	 */

    	public int  buildParameterMap (int numSeries){
    		int numPars=this.distortionCalibrationData.getNumParameters();
    		int numImg=this.distortionCalibrationData.getNumImages();
    		int numTPars=this.parameterMode[numSeries].length;
    		this.reverseParameterMap=new int [numImg][numPars];
//			boolean [] selectedEnabledImagesAll=selectedImages(numSeries,true); // strictly as in series, including no valid points ones
			boolean [] selectedEnabledImages=selectedImages(numSeries);
// set defaults - -1 - "fixed", use individual parameter from this image
    		for (int i=0;i<numImg;i++) for (int j=0;j<numPars;j++) this.reverseParameterMap[i][j]=-1;
    		int vectorIndex=0;
//    		int [][] tmpMap=new int[numTPars][3]; // temporary array parameterMap[][] (will be truncated)
    		int [][] tmpMap=new int[numPars*numImg][3]; // temporary array parameterMap[][] (will be truncated)
			double masterTS=this.distortionCalibrationData.getImageTimestamp(this.masterImages[numSeries]); // timestamp of the master image
			int [] masterVectorIndices = new int [numTPars];
// iterate through all global/subcamera parameters
    		for (int numTPar=0;numTPar<numTPars;numTPar++) if (this.parameterMode[numSeries][numTPar]!=this.modeFixed){ // skip "fixed"
    			boolean isCommon=
    				(this.parameterMode[numSeries][numTPar]==this.modeCommon) ||
    				(this.parameterMode[numSeries][numTPar]==this.modeSupercommon);
    			boolean isStation=
    				(this.parameterMode[numSeries][numTPar]==this.modeStation) ||
    				(this.parameterMode[numSeries][numTPar]==this.modeSuperStation);
    			boolean isGroup=(this.parameterMode[numSeries][numTPar]==this.modeGroup);
    			int numSub=this.parameterList[numTPar][0]; // number of sub-camera for this total parameter index
    			int numPar=this.parameterList[numTPar][1]; // number of per-image parameter for this total parameter index
    			boolean isSubCamera=this.distortionCalibrationData.isSubcameraParameter(numPar);
    			// is it "same as" ?
    			int sameAs = this.parameterMode[numSeries][numTPar] - this.definedModesAll.length; // number of subcamera to use
    			int numTParMaster = -1;
    			int vectorIndexMaster = -1;
    			if (sameAs < 0  ) sameAs = -1;
    			else {
    				for (numTParMaster = 0; numTParMaster < numTPar; numTParMaster++){
    					if ((this.parameterList[numTParMaster][0] == sameAs) && (this.parameterList[numTParMaster][1] == numPar)) break;
    				}
    				// find master index for this parameter
    				if (numTParMaster == numTPar) { // bug, reposrt and terat as fixed
    					System.out.println("Bug in buildParameterMap("+numSeries+"): could not find master parameter for numTPar="+numTPar+", sameAs="+sameAs);
    					System.out.println("Treaing as fixed");
    					continue;
    				}
    				vectorIndexMaster = masterVectorIndices[numTParMaster]; // use as vectorIndex for slaves
    			}

    			if (this.debugLevel>2) System.out.println("numTPar="+numTPar+" numSub="+numSub+" numPar="+numPar+" numTParMaster="+numTParMaster);
// iterate through available images
    			for (int numThisImg=0;numThisImg<numImg;numThisImg++) {
    				if ((selectedEnabledImages[numThisImg]) &&
    						(!isSubCamera || (numSub==this.distortionCalibrationData.getImageSubcamera(numThisImg))) &&
    						(this.reverseParameterMap[numThisImg][numPar]<0)){ // image used, this cell is not (yet) defined
    					if (this.debugLevel>2){
    						System.out.println("buildParameterMap("+numSeries+"): numThisImg="+numThisImg+", numPar="+numPar+", vectorIndex="+vectorIndex);
    					}
// is it "same as" ?
    					if (vectorIndexMaster >= 0){
    						this.reverseParameterMap[numThisImg][numPar] = vectorIndexMaster;
    						// does not look for same parameter in later images - leave it for the master?
    					} else {
// assign it a new parameter
    						this.reverseParameterMap[numThisImg][numPar]=vectorIndex;
    						masterVectorIndices[numTPar] = vectorIndex;
    						double thisTS=this.distortionCalibrationData.getImageTimestamp(numThisImg);
    						int thisStation=this.distortionCalibrationData.gIP[numThisImg].getStationNumber();
// set pointer to this first image
    						tmpMap[vectorIndex][0]=numThisImg; // vectorindex==22 > tmpMap.length?
    						tmpMap[vectorIndex][1]=numPar;
    						tmpMap[vectorIndex][2]=this.parameterMode[numSeries][numTPar];
    						double minDist=Math.abs(this.distortionCalibrationData.getImageTimestamp(numThisImg)-masterTS);
    						if (this.debugLevel>2) System.out.println("vectorIndex="+vectorIndex+" numThisImg="+numThisImg);
// see if same parameter in some other (later) image(s) is shared
    						for (int numOtherImg=numThisImg + 1; numOtherImg < numImg; numOtherImg++)
    							if ((selectedEnabledImages[numOtherImg]) && // OOB 1
    									(!isSubCamera || (numSub==this.distortionCalibrationData.getImageSubcamera(numOtherImg))) &&
    									(this.reverseParameterMap[numOtherImg][numPar]<0)){ // image used, this cell is not (yet) defined
    								if ((this.distortionCalibrationData.getImageTimestamp(numOtherImg)==thisTS) ||  // same parameter same timestamp - same group even if is set differently
    										(isStation && (this.distortionCalibrationData.gIP[numOtherImg].getStationNumber()==thisStation)) || // new
    										isCommon ||
    										(isGroup && (this.parameterGroups[numSeries][numTPar][numThisImg]==
    										this.parameterGroups[numSeries][numTPar][numOtherImg]))){
// assign it a the same parameter
    									this.reverseParameterMap[numOtherImg][numPar]=vectorIndex;
    									double thisDist=Math.abs(this.distortionCalibrationData.getImageTimestamp(numThisImg)-masterTS);
    									if (thisDist<minDist) {
    										minDist=thisDist;
    										tmpMap[vectorIndex][0]=numOtherImg;
    									}
    								}
    							}
    						vectorIndex++;
    					}
    				}
    			}

    		}
    		// reverseParameterMap built,   vectorIndex equals to the total number of parameters needed for fitting
    		//truncate tmpMap into this.parameterMap[][]
    		this.parameterMap=new int[vectorIndex][];
    		for (int i=0;i<vectorIndex;i++) this.parameterMap[i] =tmpMap[i];
    		this.currentSeriesNumber=numSeries;
		   	if (this.debugLevel>2) System.out.println("this.parameterMap.length="+this.parameterMap.length);
    		return this.parameterMap.length;
    	}

    	/**
    	 * Prepare data for calculating additional LMA terms for parameter variances
    	 * @param numSeries fitting series to use
    	 * @return number of parameters for which variance is considered
    	 */
    	public int  buildVariancesMaps (
    			int numSeries){
    		if (this.debugLevel>0){
    			System.out.println("buildVariancesMaps("+numSeries+")");
    		}
    		if ((this.varianceModes==null) || (this.varianceModes[numSeries]==varianceModeDisabled)){
    				this.averageCellIndex=null;
    				this.variationsAverages=null;
    				this.weightVariance=null;
    				this.varianceErrorsSquared=null;
    				return 0;

    		}
    		int debugThreshold=2;
    		boolean useSetWeights=(this.varianceModes[numSeries]==varianceModeVariableWeight);
    		int numPars=this.distortionCalibrationData.getNumParameters();
    		int numStations=this.distortionCalibrationData.eyesisCameraParameters.getNumStations();
//    		int numTPars=this.parameterMode[numSeries].length;
    		int numSeriesPars=this.parameterMap.length;
			// count tilt positions per station
			List<Integer>  tiltList= new ArrayList<Integer>(1000);
			Integer iStationTilt;
			int tiltMotorIndex=2;
			for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
				if (this.parameterMap[numSPar]==null){
					System.out.println("buildVariancesMaps() BUG - ,this.parameterMap["+numSPar+"]==null, numSeriesPars="+numSeriesPars);
					continue;
				}
				int imgNumber=this.parameterMap[numSPar][0]; // null pointer for triclops after adding variances for tilt
				int setNumber=this.distortionCalibrationData.gIP[imgNumber].getSetNumber();
				int station=  this.distortionCalibrationData.gIP[imgNumber].getStationNumber();
				if (this.distortionCalibrationData.gIS[setNumber]==null){
					System.out.println("buildVariancesMaps() BUG - ,this.distortionCalibrationData.gIS["+setNumber+"]==null, numSeriesPars="+numSeriesPars+
							" numSPar="+numSPar+" imgNumber="+imgNumber+" station="+station);
					continue;
				}
				if (this.distortionCalibrationData.gIS[setNumber].motors==null){
					System.out.println("buildVariancesMaps() BUG - ,this.distortionCalibrationData.gIS["+setNumber+"].motors==null, numSeriesPars="+numSeriesPars+
							" numSPar="+numSPar+" imgNumber="+imgNumber+" station="+station);
					continue;
				}


				int tiltMotor=this.distortionCalibrationData.gIS[setNumber].motors[tiltMotorIndex];// null pointer for triclops after adding variances for tilt


				iStationTilt=station+numStations*tiltMotor;
				if (!tiltList.contains(iStationTilt)) tiltList.add(iStationTilt);
			}
			int numStationTilts=tiltList.size();

			int [] numGroups=new int [numPars];
			for (int i=0;i<numGroups.length;i++)numGroups[i]=0;
			for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
				int parIndex=this.parameterMap[numSPar][1];
//				if ((numGroups[parIndex]==0) && (this.distortionCalibrationData.eyesisCameraParameters.isVarianceCostSet(parIndex))){
				if (numGroups[parIndex]==0){
					boolean isSet=this.distortionCalibrationData.eyesisCameraParameters.isVarianceCostSet(parIndex);
					switch (this.parameterMap[numSPar][2]){
					case modeWeakCommon:
						if (isSet) numGroups[parIndex]=1;
						if (this.debugLevel>debugThreshold) System.out.println(">>>1 "+numSPar+":"+parIndex +" - "+numGroups[parIndex]+" isSet="+isSet);
						break;
					case modeWeakStation:
						if (isSet) numGroups[parIndex]=numStations;
						if (this.debugLevel>debugThreshold) System.out.println(">>>2 "+numSPar+":"+parIndex +" - "+numGroups[parIndex]+" isSet="+isSet);
						break;
					case modeTiltEqualize:
						if (isSet) numGroups[parIndex]=numStationTilts;
						if (this.debugLevel>debugThreshold) System.out.println(">>>3 "+numSPar+":"+parIndex +" - "+numGroups[parIndex]+" isSet="+isSet);
						break;
					}
				}
			}
    		if (this.debugLevel>debugThreshold){
    			System.out.println("buildVariancesMaps() numGroups:");
    			for (int i=0;i<numGroups.length;i++) if (numGroups[i]>0) System.out.println("--- "+i+":"+numGroups[i]);
    		}

			int startIndex=0;
			for (int i=0;i<numGroups.length;i++){
				if (numGroups[i]>0){
					int num=numGroups[i];
					numGroups[i]=startIndex;
					startIndex+=num;
				} else {
					numGroups[i]=-1;
				}
			}
    		if (this.debugLevel>debugThreshold){
    			System.out.println("buildVariancesMaps() start indices:");
    			for (int i=0;i<numGroups.length;i++) if (numGroups[i]>=0) System.out.println("--- "+i+":"+numGroups[i]);
    		}

    		if (this.debugLevel>debugThreshold){
    			System.out.println("buildVariancesMaps() startIndex="+startIndex);
    		}

			if (startIndex==0){ // no variance parameters defined
				this.averageCellIndex=null;
				this.variationsAverages=null;
				this.weightVariance=null;
				this.varianceErrorsSquared=null;
				return 0;
			}

			// Assign average cell indices for each vector element, for each paramter it can be just 1 for global, numStations or number of tilt positions
			this.averageCellIndex=new int [numSeriesPars];
			for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
				this.averageCellIndex[numSPar]=-1; // not applicable
				int parIndex=this.parameterMap[numSPar][1];
//			for (int numTPar=0;numTPar<numTPars;numTPar++) {
//				int parIndex=this.parameterList[numTPar][1];
				if ((numGroups[parIndex]>=0) && (this.distortionCalibrationData.eyesisCameraParameters.isVarianceCostSet(parIndex))){
					int imgNumber=this.parameterMap[numSPar][0];
					int setNumber=this.distortionCalibrationData.gIP[imgNumber].getSetNumber();
					int station=  this.distortionCalibrationData.gIP[imgNumber].getStationNumber();

					switch (this.parameterMap[numSPar][2]){
					case modeWeakCommon:
						this.averageCellIndex[numSPar]=numGroups[parIndex];
						break;
					case modeWeakStation:
						this.averageCellIndex[numSPar]=numGroups[parIndex]+station;
						break;
					case modeTiltEqualize:
						int tiltMotor=this.distortionCalibrationData.gIS[setNumber].motors[tiltMotorIndex];
						iStationTilt=station+numStations*tiltMotor;
						this.averageCellIndex[numSPar]=numGroups[parIndex]+tiltList.indexOf(iStationTilt);
						break;
					}
				}
			}
    		if (this.debugLevel>1){
    			System.out.println("buildVariancesMaps() numStationTilts="+numStationTilts+" numSeriesPars="+numSeriesPars+ " useSetWeights="+useSetWeights+
    					" number of averages="+startIndex);
    		}

			this.variationsAverages=new double [startIndex]; // total number of different averages
			int [] numContributors=new int [startIndex];
			for (int i=0;i<this.variationsAverages.length;i++) {
				this.variationsAverages[i]=0.0;
				numContributors[i]=0;
			}
			// Calculate total weights for averages and number of contributors (to remove single-contributor items)
			for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
				int avIndex=this.averageCellIndex[numSPar];
				if (avIndex>=0){
					double w;
					if (useSetWeights){
						int imgNumber=this.parameterMap[numSPar][0];
						int setNumber=this.distortionCalibrationData.gIP[imgNumber].getSetNumber();
						w=this.distortionCalibrationData.gIS[setNumber].getSetWeight();
					} else {
						w=1.0;
					}
					if (w>0.0) { // do not count zero-contributors (are they possible?)
						this.variationsAverages[avIndex]+=w;
						numContributors[avIndex]++;
					}
				}
			}
    		if (this.debugLevel>debugThreshold){
    			System.out.println("buildVariancesMaps() numContributors:");
    			for (int i=0;i<numContributors.length;i++) if (numContributors[i]>0) System.out.println("--- "+i+":"+numContributors[i]);
    		}

			// calculate each parameter share of the average TODO: disable if share =100% (single contributor)
			int numVariancePars=0;
			this.weightOnAverage=new double [numSeriesPars];
			this.weightVariance= new double [numSeriesPars];
			this.varianceErrorsSquared= new double [numSeriesPars];

			for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
				this.weightOnAverage[numSPar]=0.0;
				this.weightVariance[numSPar]=0.0;
				this.varianceErrorsSquared[numSPar]=Double.NaN;
				int avIndex=this.averageCellIndex[numSPar];
				if (avIndex>=0){
					if (numContributors[avIndex]<2){
						this.averageCellIndex[numSPar]=-1; // single element - impossible to use
						if (this.debugLevel>0) System.out.println("buildVariancesMaps(): removed parameter #"+numSPar+" - single contributor to average");
					} else {
						numVariancePars++;
						double w;
						if (useSetWeights){
							int imgNumber=this.parameterMap[numSPar][0];
							int setNumber=this.distortionCalibrationData.gIP[imgNumber].getSetNumber();
							w=this.distortionCalibrationData.gIS[setNumber].getSetWeight();
						} else w=1.0;
						this.weightOnAverage[numSPar]=w/this.variationsAverages[avIndex];
						this.weightVariance[numSPar]= w;
					}
				}
			}
    		if (this.debugLevel>debugThreshold){
    			System.out.println("average indices/shares:");
    			for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
    				int avIndex=this.averageCellIndex[numSPar];
    				if (avIndex>=0){
    					int parIndex=this.parameterMap[numSPar][1];
    					System.out.println("==="+numSPar+" parIndex="+parIndex+" avIndex="+avIndex+
    							" weightOnAverage="+this.weightOnAverage[numSPar]+
    							" weightVariance="+this.weightVariance[numSPar]);
    				}
    			}
    		}
			if (numVariancePars==0){
				this.averageCellIndex=null;
				this.variationsAverages=null;
				this.weightVariance=null;
				this.varianceErrorsSquared=null;
				return 0;
			}
			// TODO this.variationsAverages needs to be reset to all zeros each calculations
			return numVariancePars;
    	}
    	/**
    	 * Ammend LMA arrays to pull individual values of the extrinsic parameters to their averages
    	 * Averages may be global, per-station or per same tilt motor position (same station)
    	 * The function to be minimized is
    	 * f(x) Scale*(1/p*Xeff+(1-1/p)*pow(Xeff,p)
    	 * Xeff=x/X0*(1-shareInAverage)
    	 * and x is difference between the poarameter and it's average over the specified set (all/station/tiltStation)
    	 * @param numSeries fitting series number
    	 * @param vector - parameters vector
    	 * @param jTByJ Jacobian transposed multiplied by Jacobian - some diagonal elements may be modified (weighted)  or null (will add to old values)
    	 * @param jTByDiff Jacobian transposed multiplied by difference vector (weighted)  or null (will add to old values)
    	 * @return true if parameter variances are applicable to the fitting series
    	 */
    	public boolean addVarianceToLMA(
    			int numSeries,
    			final double [] vector,
    	        double [][] jTByJ, // jacobian multiplied by Jacobian transposed
    	        double []   jTByDiff){ // jacobian multiplied difference vector
    		if (this.averageCellIndex==null) return false; // parameter varainces are not used
//    		int numTPars=this.parameterMode[numSeries].length;
    		int numSeriesPars=this.parameterMap.length;
    		int debugThreshold=2;

//    		double [] fX=   new double [numTPars];
//    		double [] diagJ=new double [numTPars]; // diagonal of Jacobian (here parameters only influence themselves)
    		double [] weights=new double [this.variationsAverages.length];
    		for (int i=0;i<this.variationsAverages.length;i++) {
    			this.variationsAverages[i]=0.0;
    			weights[i]=0.0;
    		}
    		// calculate appropriate averages
    		for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
    			this.varianceErrorsSquared[numSPar]=0.0;
				int avIndex=this.averageCellIndex[numSPar];
				if (avIndex>=0){
					this.variationsAverages[avIndex]+=this.weightOnAverage[numSPar]*vector[numSPar];
    			}
    		}
    		// now fX and diagJ
    		int numModded=0;
    		for (int numSPar=0;numSPar<numSeriesPars;numSPar++) {
				int avIndex=this.averageCellIndex[numSPar];
				int parIndex=this.parameterMap[numSPar][1];
//    			if ((avIndex>=0) &&  (this.distortionCalibrationData.eyesisCameraParameters.isVarianceCostSet(parIndex))) {
    			if (avIndex>=0) {
    				double scale=this.distortionCalibrationData.eyesisCameraParameters.varianceCostScale(parIndex);
    				double x0=   this.distortionCalibrationData.eyesisCameraParameters.varianceCostVariationAbs(parIndex);
    				double exp=  this.distortionCalibrationData.eyesisCameraParameters.varianceCostVariationExponent(parIndex);
    				double xEff=(vector[numSPar]-this.variationsAverages[avIndex])*(1.0-this.weightOnAverage[numSPar])/x0; //OOB-1
    				double axEff=Math.abs(xEff);
    				double sxEff=Math.signum(xEff);
    				double k=1.0/exp;
    				double fX=sxEff*scale*(k*axEff+(1.0-k)*Math.pow(axEff,exp));
    				double diagJ=scale*((1.0-this.weightOnAverage[numSPar])/x0)*(k+(exp-1)*Math.pow(axEff,exp-1.0));
    				this.varianceErrorsSquared[numSPar]=this.weightVariance[numSPar]*fX*fX;
    				if (jTByJ!=null) jTByJ[numSPar][numSPar]+=this.weightVariance[numSPar]*diagJ*diagJ;
    				if (jTByDiff!=null) jTByDiff[numSPar]-= this.weightVariance[numSPar]*fX*diagJ;
    				if (this.debugLevel>debugThreshold){
    					System.out.print  (numSPar+" "+parIndex+" "+avIndex+": "+" scale="+scale);
    					System.out.print  (" vector["+numSPar+"]="+vector[numSPar]+" ("+this.variationsAverages[avIndex]+") xEff="+xEff+" fX="+fX+" k="+k+" exp="+exp );
    					System.out.print  (" diagJ="+diagJ+" "+xEff+" jTByJ+"+ (this.weightVariance[numSPar]*diagJ*diagJ) );
    					System.out.print  (" jTByDiff+"+(-this.weightVariance[numSPar]*fX*diagJ)+" ves="+(this.weightVariance[numSPar]*fX*fX));
    					System.out.println(" weight="+this.weightVariance[numSPar]);
    				}
    				if (this.weightVariance[numSPar]!=0.0) numModded++;
    			}
    		}
			if (this.debugLevel>1)System.out.println ("addVarianceToLMA() - modified "+numModded+" elements");
    		return true;
    	}
    	public double [] getVarianceError2(){ return (this.averageCellIndex==null)?null:this.varianceErrorsSquared;}
    	public double [] getWeights(){ return (this.averageCellIndex==null)?null:this.weightVariance;}

    	/**
    	 *
    	 * @return number of the current fitting strategy series
    	 */
    	public int getCurrentSeries(){
    		return this.currentSeriesNumber;
    	}

    	/**
    	 * Calculate vector of the parameters used in LMA algorithm, extracted from the
    	 * individual data, using parameter map (calculated once after changing series)
    	 * @return vector of parameters used for the LMA
    	 */
    	public double [] getSeriesVector(){
        	if (this.parameterMap==null) {
        		String msg="Series parameters map is not calculated";
        		IJ.showMessage("Error",msg);
        		throw new IllegalArgumentException (msg);
        	}
    		double [] vector=new double[this.parameterMap.length];
    		for (int i=0;i<vector.length;i++) {
//    			vector[i]=this.distortionCalibrationData.pars[this.parameterMap[i][0]][this.parameterMap[i][1]];
    			vector[i]=this.distortionCalibrationData.getParameterValue(this.parameterMap[i][0],this.parameterMap[i][1]); // (numImg, numPar)
    		}
    		return vector;
    	}

    	/**
    	 * Saves data from the parameter vector to that of the images
    	 * @param vector vector of parameters (after LMA fitting)
    	 */
    	// TODO: Update the temporarily disabled images also, when possible (modify buildParameterMap also?
    	public void saveSeriesVector(double [] vector){
    		if ((this.parameterMap==null) ||(vector==null) || (vector.length!=this.parameterMap.length) ) {
    			String msg="Vector length does not match parameters length";
    			IJ.showMessage("Error",msg);
    			throw new IllegalArgumentException (msg);
    		}
    		// save for "individual" and "common" - and "station" (week - same as individual here)
    		boolean [] selectedEnabledImages=selectedImages(this.currentSeriesNumber);
    		for (int numImg=0;numImg<this.reverseParameterMap.length;numImg++)	if (selectedEnabledImages[numImg]){
    			for (int nPar=0;nPar<this.reverseParameterMap[numImg].length;nPar++) if (this.reverseParameterMap[numImg][nPar]>=0){
//    				this.distortionCalibrationData.pars[numImg][nPar]=vector[this.reverseParameterMap[numImg][nPar]];
    				this.distortionCalibrationData.setParameterValue(numImg,nPar,vector[this.reverseParameterMap[numImg][nPar]],true);

    				if (this.debugLevel>2){
    					System.out.println(" Updated image "+numImg+" "+distortionCalibrationData.getParameterName(nPar)+" = "+
    				//			this.distortionCalibrationData.pars[numImg][nPar]) ; //vector[this.reverseParameterMap[numImg][nPar]]);
    							this.distortionCalibrationData.getParameterValue(numImg,nPar));
    				}
    			}
    		}
    		// propagate to other (unselected) images for "super common" parameters
    		for (int vPar=0; vPar<this.parameterMap.length;vPar++) {
    			if (
    					(this.parameterMap[vPar][2]==modeSupercommon) ||
    					(this.parameterMap[vPar][2]==modeSuperStation)	){ // "super common" No "weak" here!

    				//        		int nSub=this.parameterMap[vPar][0]; // that's an image number, not a subcamera number
    				int nSub=this.distortionCalibrationData.gIP[this.parameterMap[vPar][0]].channel;
    				int nStation=this.distortionCalibrationData.gIP[this.parameterMap[vPar][0]].getStationNumber();
    				boolean superCommon=(this.parameterMap[vPar][2]==modeSupercommon);
    				int nPar=this.parameterMap[vPar][1];
    				//this.distortionCalibrationData.channels[i]
    				boolean isSubCamera=this.distortionCalibrationData.isSubcameraParameter(nPar);
    				for (int numImg=0;numImg<this.reverseParameterMap.length;numImg++) {
    					if (superCommon || (this.distortionCalibrationData.gIP[numImg].getStationNumber()==nStation)) {

    						if (this.debugLevel>2){
    							System.out.println("saveSeriesVector(): numImg="+numImg+" nPar="+nPar+" isSubCamera="+isSubCamera+
    									" nSub="+nSub+" distortionCalibrationData.getImageSubcamera("+numImg+")="+distortionCalibrationData.getImageSubcamera(numImg)+
    									" vector["+vPar+"]="+vector[vPar]);
    						}
    						if ((!isSubCamera) || (nSub==this.distortionCalibrationData.getImageSubcamera(numImg))) {
//    							this.distortionCalibrationData.pars[numImg][nPar]=vector[vPar];
    							this.distortionCalibrationData.setParameterValue(numImg,nPar,vector[vPar],true);
    						}
    					}
    				}
    			}
    		}
    	}

    	/**
    	 * Calculates current values of all parameters for the particular sensor - some ("fixed")
    	 * are taken from the data stored for this individual image, others - from the parameter
    	 * vector (used in fitting) UPDATE: works with null parameterVector
    	 * @param numImg number of image
    	 * @param vector parameters vector
    	 * @return vector used for the current image (parameters influencing the acquired grid
    	 * on the sensor (common parameters and those of the sensor's subchannel)
    	 */
    	public double [] getImageParametersVector(int numImg, double [] parameterVector){
//    		double [] vector = this.distortionCalibrationData.pars[numImg].clone();
    		double [] vector = this.distortionCalibrationData.getParameters(numImg); // returns a copy, no clone() is needed
    		if ((parameterVector!=null) && (this.reverseParameterMap!=null) && (this.reverseParameterMap[numImg]!=null)){
    			for (int i=0;i<vector.length;i++){
    				if (this.reverseParameterMap[numImg][i]>=0) vector[i]=parameterVector[this.reverseParameterMap[numImg][i]];
    			}
    		}
    		return vector;
    	}
    	/**
    	 * Calculates which of all parameters for the particular sensor are to be adjusted (to reduce calcualtions)
    	 * @param numImg number of image
    	 * @param vector parameters vector
    	 * @return mask vector to be used with the results of getImageParametersVector()
    	 */
    	public boolean [] getImageParametersVectorMask (int numImg){
    		if ((this.reverseParameterMap==null) || (this.reverseParameterMap[numImg]==null)) return null;
//    		boolean [] mask =new boolean [this.distortionCalibrationData.pars[numImg].length];
    		boolean [] mask =new boolean [this.distortionCalibrationData.getParametersLength(numImg)];
    		for (int i=0;i<mask.length;i++){
    			mask[i]=(this.reverseParameterMap[numImg][i]>=0);
    		}
    		return mask;
    	}
    	/**
    	 * Calculates index in the parameter vector corresponding to each image parameter vector element
    	 * @param numImg number of image
    	 * @param vector parameters vector
    	 * @return mask vector to be used with the results of getImageParametersVector()
    	 */
    	public int [] getImageParametersVectorReverseMap (int numImg){
    		if ((this.reverseParameterMap==null) || (this.reverseParameterMap[numImg]==null)) return null;
//    		int  [] map =new int [this.distortionCalibrationData.pars[numImg].length];
    		int  [] map =new int [this.distortionCalibrationData.getParametersLength(numImg)];
    		for (int i=0;i<map.length;i++){
    			map[i]=this.reverseParameterMap[numImg][i];
    		}
    		return map;
    	}
    	/**
    	 * Opens a text window with a table that shows map from each element of the parameters vector
    	 * to the image number used as a source of this parameter
    	 * @param title Window title
    	 */
    	public void showCurrentParameterMap(String title){
    	    String header="#\tMode\tName\tSubcamera\tImage Number\tMode";
    	    StringBuffer sb = new StringBuffer();
    	    for (int i=0; i<this.parameterMap.length;i++){
    	        sb.append(
    	        		i+"\t"+
    	        		this.definedModesAll[this.parameterMap[i][2]]+"\t"+
    	        		this.distortionCalibrationData.getParameterName(this.parameterMap[i][1])+"\t"+
    	        		this.distortionCalibrationData.getImageSubcamera(this.parameterMap[i][0])+"\t"+
    	        		this.parameterMap[i][0]+"\t"+
    	        		this.definedModesAll[this.parameterMap[i][2]]+"\n");
    	    }
    	    new TextWindow(title, header, sb.toString(), 800,600);
    	}
    	/**
    	 * Opens window with a table mapping each image (column) parameter (row) to parameter vector element
    	 * "-" (-1 in teh table) - fixed parameter, use one from the individual image
    	 * @param title Window title
    	 */
    	public void showCurrentReverseParameterMap(String title){
    	    String header="Image";
			boolean [] selectedEnabledImages=selectedImages(this.currentSeriesNumber);
    	    for (int i=0; i<this.reverseParameterMap.length;i++) if (selectedEnabledImages[i]){
    	    	header+="\t"+i;
    	    }
    	    StringBuffer sb = new StringBuffer();
    	    sb.append("timestmps");
    	    for (int i=0; i<this.reverseParameterMap.length;i++) if (selectedEnabledImages[i]){
    	    	sb.append("\t"+IJ.d2s(this.distortionCalibrationData.getImageTimestamp(i),6));
    	    }
    	    sb.append("\n");
    	    sb.append("subcamera");
    	    for (int i=0; i<this.reverseParameterMap.length;i++) if (selectedEnabledImages[i]){
    	    	sb.append("\t"+this.distortionCalibrationData.getImageSubcamera(i));
    	    }
    	    sb.append("\n");

    	    for (int k=0;k<this.reverseParameterMap[0].length;k++){
//    			boolean isSubCamera=this.distortionCalibrationData.isSubcameraParameter(k);
//        	    sb.append(this.distortionCalibrationData.getParameterName(this.parameterMap[k][1])); // name of parameter
        	    sb.append(this.distortionCalibrationData.getParameterName(k)); // name of parameter
        	    for (int i=0; i<this.reverseParameterMap.length;i++) if (selectedEnabledImages[i]){
        	    	int mode =this.reverseParameterMap[i][k];
        	    		sb.append("\t"+((mode>=0)?mode:"-"));

        	    }
        	    sb.append("\n");

    	    }
    	    new TextWindow(title, header, sb.toString(), 800,600);
    	}

    	/**
    	 * Add/reduce number of series in this fitting strategy
    	 * @param numSeries new number of series
    	 */
    	public void setLength(int numSeries){
    		if (numSeries<=0) {
    			setDflt(0);
    			return;
    		}
    		if ((this.selectedImages!=null) && (numSeries==this.selectedImages.length)) return;
    		if ((this.selectedImages==null) && (numSeries==0)) return;
        	boolean [][] oldSelectedImages=null; // images selected for each step
        	int [][] oldParameterMode=null;
        	int [] oldVarianceModes=null;
        	int [][] oldZGroups=null;
        	int [][][] oldParameterGroups=null;
        	double [] oldLambdas=null;   // LMA initial lambda for each step
    		int [] oldMasterImages=null;
        	double [] oldStepDone=null;
        	boolean [] oldStopAfterThis=null;
        	String [] oldStrategyComment=null;
            int oldNumSeries=(this.selectedImages==null)?0:this.selectedImages.length;
            int oldLength=0;
    		if (this.selectedImages!=null){
    			// deep clone this.* to old*;
    			oldSelectedImages=new boolean[oldNumSeries][];
    			oldParameterMode=new int[oldNumSeries][];
    			oldParameterGroups=new int[oldNumSeries][][];
    			oldZGroups=new int[oldNumSeries][];
        		oldStrategyComment = new String [oldNumSeries];
    			for (int i=0;i<oldNumSeries;i++){
    				oldSelectedImages[i]=this.selectedImages[i].clone();
        			oldParameterMode[i]= this.parameterMode[i].clone(); // out of bound - 2
        			oldParameterGroups[i]=new int[this.parameterGroups[i].length][];
        			for (int j=0;j<oldParameterGroups[i].length;j++)
        				oldParameterGroups[i][j]=(parameterGroups[i][j]==null)?null:parameterGroups[i][j].clone();
        			if (this.zGroups[i]!=null) oldZGroups[i]=this.zGroups[i].clone();
        			else this.zGroups[i]=oldZGroups[i];
        			oldStrategyComment[i]=this.strategyComment[i];
    			}
    			oldVarianceModes=(this.varianceModes==null)?null:this.varianceModes.clone();
    			oldLambdas=this.lambdas.clone();
    			oldLength=this.selectedImages.length;
        		oldMasterImages=this.masterImages.clone();
        		oldStepDone=this.stepDone.clone();
        		oldStopAfterThis=this.stopAfterThis.clone();
    		}
    		setDflt(numSeries);
    		if (this.selectedImages!=null) for (int i=0;(i<this.selectedImages.length) && (i<oldLength);i++){
    			this.selectedImages[i]=oldSelectedImages[i];
    			this.parameterMode[i]= oldParameterMode[i];
    			this.parameterGroups[i]=oldParameterGroups[i];
    			this.lambdas[i]=oldLambdas[i];
    			this.masterImages[i]=oldMasterImages[i];
        		this.stepDone[i]=oldStepDone[i];
        		this.stopAfterThis[i]=oldStopAfterThis[i];
    			this.varianceModes[i]=oldVarianceModes[i];
    			this.zGroups[i]=(oldZGroups[i]!=null)?oldZGroups[i].clone():null;
    			this.strategyComment[i] = oldStrategyComment[i];
    		}
    	}

    	public void updateNumberOfSubcameras(){  // can break "same as"?
    		int numPars=   this.distortionCalibrationData.getNumParameters();
    		int numSubCams=this.distortionCalibrationData.getNumSubCameras();
    		initParameterList();
    		int numSubPars=0;
    		int [] subParIndex=new int[numPars];
    		for (int i=0;i<numPars;i++){
    			if (this.distortionCalibrationData.isSubcameraParameter(i)) subParIndex[numSubPars++]=i;
    		}
    		int totalNumPars=numPars+(numSubCams-1)*numSubPars; // maximal number of parametes per timestamp?
    		int oldTotalPars=this.parameterEnable.length;
    		if (oldTotalPars==totalNumPars) return;
    		if (oldTotalPars>numPars) { // more than one subcamera, so last numSubPars are subcamera parameters
        		for (int i=0;i<numSubPars;i++)subParIndex[i]=(oldTotalPars-numSubPars)+i;
    		}
//       public boolean isSubcameraParameter(int num){

    		int numSeries=this.parameterMode.length;
    		if (oldTotalPars<totalNumPars){ // grow, repeat last subcamera
    			if (this.debugLevel>1) System.out.println("Increasing total number of parameters in fitting strategy: "+oldTotalPars+" -> "+totalNumPars);
    			for (int nSer=0;nSer<numSeries;nSer++){
    				int [] parameterMode_tmp=this.parameterMode[nSer]; //.clone();
    				int [][] parameterGroups_tmp=this.parameterGroups[nSer]; //.clone();
            		this.parameterMode[nSer]=new int[totalNumPars];
            		this.parameterGroups[nSer]=new int[totalNumPars][];
            		for (int nTPar=0;nTPar<oldTotalPars;nTPar++){
            			this.parameterMode[nSer][nTPar]=parameterMode_tmp[nTPar];
            			this.parameterGroups[nSer][nTPar]=parameterGroups_tmp[nTPar];
            		}
            		for (int nTPar=oldTotalPars;nTPar<totalNumPars;nTPar++){
            			int index= subParIndex[(nTPar-oldTotalPars) % numSubPars];
            			this.parameterMode[nSer][nTPar]=parameterMode_tmp[index];
            			this.parameterGroups[nSer][nTPar]=(parameterGroups_tmp[index]==null)?null:parameterGroups_tmp[index].clone();
            		}
    			}
				boolean []  parameterEnable_tmp=this.parameterEnable; //.clone();
        		this.parameterEnable= new boolean [totalNumPars];
        		for (int nTPar=0;nTPar<oldTotalPars;nTPar++){
        			this.parameterEnable[nTPar]=parameterEnable_tmp[nTPar];
        		}
        		for (int nTPar=oldTotalPars;nTPar<totalNumPars;nTPar++){
        			int index= subParIndex[(nTPar-oldTotalPars) % numSubPars];
        			this.parameterEnable[nTPar]=parameterEnable_tmp[index];
        		}
    		} else { // shrink
    			if (this.debugLevel>1) System.out.println("Reducing total number of parameters in fitting strategy: "+oldTotalPars+" -> "+totalNumPars);
    			for (int nSer=0;nSer<numSeries;nSer++){
    				int [] parameterMode_tmp=this.parameterMode[nSer]; //.clone();
    				int [][] parameterGroups_tmp=this.parameterGroups[nSer]; //.clone();
            		this.parameterMode[nSer]=new int[totalNumPars];
            		this.parameterGroups[nSer]=new int[totalNumPars][];
            		for (int nTPar=0;nTPar<totalNumPars;nTPar++){
            			this.parameterMode[nSer][nTPar]=parameterMode_tmp[nTPar];
            			this.parameterGroups[nSer][nTPar]=parameterGroups_tmp[nTPar];
            		}
    			}
				boolean []  parameterEnable_tmp=this.parameterEnable; //.clone();
        		this.parameterEnable= new boolean [totalNumPars];
        		for (int nTPar=0;nTPar<totalNumPars;nTPar++){
        			this.parameterEnable[nTPar]=parameterEnable_tmp[nTPar];
        		}
    		}
    	}

    	private void setDflt(int numSeries){
    		if (numSeries==0) {
        		this.selectedImages=null;
        		this.lambdas=null;
        		this.stepDone=null;
        		this.parameterMode=null;
        		this.parameterGroups=null;
        		this.varianceModes=null;
        		this.zGroups=null;
        		this.strategyComment=null;
                return;
    		}
    		this.stopAfterThis=new boolean[numSeries];
    		for (int i=0;i<numSeries;i++) this.stopAfterThis[i]=true;

// calculate total (potential) number of parameters in the system
//TODO: modify for groups
    		int numPars=   this.distortionCalibrationData.getNumParameters();
    		int numSubCams=this.distortionCalibrationData.getNumSubCameras();
    		int numSubPars=0;
    		for (int i=0;i<numPars;i++){
    			if (this.distortionCalibrationData.isSubcameraParameter(i)) numSubPars++;
    		}
    		int totalNumPars=numPars+(numSubCams-1)*numSubPars; // maximal number of parametes per timestamp?
    		// number of parameters to adjust may be greater, if they were changing between the images
    		this.selectedImages=new boolean[numSeries][this.distortionCalibrationData.getNumImages()];//
    		this.selectedValidImages=new boolean[numSeries][this.distortionCalibrationData.getNumImages()];//
    		this.masterImages=new int[numSeries];
    		this.lambdas=new double[numSeries];
    		this.stepDone=new double[numSeries];
    		this.parameterMode=new int[numSeries][totalNumPars];
    		this.varianceModes=new int[numSeries];
    		this.parameterGroups=new int[numSeries][totalNumPars][];
    		this.parameterEnable= new boolean [totalNumPars];
    		this.zGroups=new int[numSeries][];
    		this.strategyComment = new String [numSeries];
    		for (int i=0;i<this.selectedImages.length;i++){
    			invalidateSelectedImages(i);
    			for (int j=0;j<this.selectedImages[i].length;j++){
    				this.selectedImages[i][j]=false;
    			}
    			for (int j=0;j<this.parameterMode[i].length;j++){
    				this.parameterMode[i][j]=this.modeFixed; // fixed
    				this.parameterGroups[i][j]=null; // no groups yet defined
    			}
    			this.lambdas[i]=defaultLambda;
    			this.stepDone[i]=defaultStepDone;
    			this.masterImages[i]=0; // first image, supposingly the earliest timestamp
    			this.varianceModes[i]=varianceModeDisabled;
        		this.zGroups[i]=null;
        		this.strategyComment[i] = "";
    		}

    		initParameterList();

    		for (int i=0;i<numPars;i++){
    			this.parameterEnable[i]=true; // initially enable all parameters for the first camera
    		}
    		for (int i=1;i<numSubCams;i++){
    			int i1=numPars+numSubPars*(i-1);
    			int j1=0;
    			for (int j=0;j<numPars;j++) if (this.distortionCalibrationData.isSubcameraParameter(j)){
    				this.parameterEnable[i1+j1]=false; // initially disable all other subcamera parameters
        			j1++;
    			}
    		}
    	}

    	private void initParameterList(){
    		int numPars=   this.distortionCalibrationData.getNumParameters();
    		int numSubCams=this.distortionCalibrationData.getNumSubCameras();
    		int numSubPars=0;
    		for (int i=0;i<numPars;i++){
    			if (this.distortionCalibrationData.isSubcameraParameter(i)) numSubPars++;
    		}
    		int totalNumPars=numPars+(numSubCams-1)*numSubPars; // maximal number of parametes per timestamp?
    		this.parameterList=new int[totalNumPars][2];

    		for (int i=0;i<numPars;i++){
    			this.parameterList[i][0]=0;
    			this.parameterList[i][1]=i;
//    			this.parameterEnable[i]=true; // initially enable all parameters for the first camera
    		}

    		for (int i=1;i<numSubCams;i++){
    			int i1=numPars+numSubPars*(i-1);
    			int j1=0;
    			for (int j=0;j<numPars;j++) if (this.distortionCalibrationData.isSubcameraParameter(j)){
    				this.parameterList[i1+j1][0]=i;
    				this.parameterList[i1+j1][1]=j;
 //   				this.parameterEnable[i1+j1]=false; // initially disable all other subcamera parameters
        			j1++;
    			}
    		}

    	}
    	/**
    	 * Find parameter number from subcamera and index
    	 * @param nSub number of subcamera (channel)
    	 * @param index index of parameter
    	 * @return number of parameter that has specified sub-camera and index. -1 if such does not exist
    	 */
    	private int getParameterNumber(int nSub, int index){
    		for (int i=0;i<this.parameterList.length;i++) if ((this.parameterList[i][0]==nSub) && (this.parameterList[i][1]==index)) return i;
    		return -1;
    	}
    	private int [] arrayToSortedInt(int [] srcArray){
    		Set<Integer> set = new HashSet<Integer>();
    		for (Integer data:srcArray) set.add(data);
    		int [] array = new int [set.size()];
    		int i=0;
    		for (Integer val:set) array[i++]= val;
    		Arrays.sort(array);
    		return array;
    	}

    	public boolean selectIndividualImages(
    			boolean [] selection,
    			boolean allImages,
    			int startIndex,
    			int perPage){
    		boolean [] enabled=  this.distortionCalibrationData.selectEnabled();
    		int [] hintedMatch=this.distortionCalibrationData.getHintedMatch();
    		if (selection.length!=enabled.length){
    			System.out.println("BUG: selectIndividualImages(): selection.length!=enabled.length!");
    			return false;
    		}
    		int endIndex=startIndex;
    		int numImg=0;
    		for (endIndex=startIndex; (endIndex<enabled.length) && (numImg<perPage);endIndex++) if (allImages || enabled[endIndex]) numImg++;
    		for (;(endIndex<enabled.length) && !allImages && !enabled[endIndex];endIndex++); // advance over disabled images
    		GenericDialog gd=new GenericDialog("Select images "+startIndex+"..."+(endIndex-1));
    		for (int i=startIndex;i<endIndex;i++) if (allImages || enabled[i]){
				gd.addCheckbox    (i+" - "+(this.distortionCalibrationData.gIP[i].enabled?"":"(disabled) ")+
						IJ.d2s(this.distortionCalibrationData.gIP[i].timestamp,6)+
						": "+this.distortionCalibrationData.gIP[i].channel+
						" matched "+this.distortionCalibrationData.gIP[i].matchedPointers+" pointers"+
						", hinted state: "+((hintedMatch[i]<0)?"undefined":((hintedMatch[i]==0)?"failed":((hintedMatch[i]==1)?"orientation":"orientation and translation"))),
						selection[i]);
    		}
    		if (endIndex<enabled.length){
         		gd.enableYesNoCancel("Done", "Next");
    		}
			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
    		for (int i=startIndex;i<endIndex;i++) if (allImages || enabled[i]){
    			selection[i]=gd.getNextBoolean();
    		}
			if (gd.wasOKed()) return true;
			return selectIndividualImages(
	    			selection,
	    			allImages,
	    			endIndex,
	    			perPage);
    	}

    	public boolean selectImageSets(
    			boolean [] selection,
    			boolean allImages,
    			int startIndex,
    			int perPage){
    		boolean [] enabled=  this.distortionCalibrationData.selectEnabled();
    		int [] imageStations=this.distortionCalibrationData.getStations();
    		int [] imageChannels=this.distortionCalibrationData.getChannels();

    		if (selection.length!=enabled.length){
    			System.out.println("BUG: selectIndividualImages(): selection.length!=enabled.length!");
    			return false;
    		}
    	    int [][] imageSets=this.distortionCalibrationData.listImages(!allImages); // true - only enabled images
    	    boolean [] enabledSets=new boolean [imageSets.length];
    	    boolean [] selectedSets=new boolean [imageSets.length]; // at least one image selected in the series
    	    for (int i=0;i<imageSets.length;i++){
    	    	enabledSets[i]=false;
    	    	selectedSets[i]=false;
    	    	if (imageSets[i]!=null) for (int j=0;j<imageSets[i].length;j++){
    	    		enabledSets[i] |= enabled[imageSets[i][j]];
    	    		selectedSets[i] |= selection[imageSets[i][j]];
    	    	}
    	    }

    		int endIndex=startIndex;
    		int numSet=0;
    		for (endIndex=startIndex; (endIndex<enabledSets.length) && (numSet<perPage);endIndex++) if (enabledSets[endIndex]) numSet++;
    		for (;(endIndex<enabledSets.length) &&  !enabledSets[endIndex];endIndex++); // advance over disabled sets
    		GenericDialog gd=new GenericDialog("Select image Sets "+startIndex+"..."+(endIndex-1));
    		for (int i=startIndex;i<endIndex;i++) if (enabledSets[i]){
    			int station=-1;
    			String sImgList="";
    			int l=0;
    			for (int j=0;j<imageSets[i].length;j++) if (enabled[imageSets[i][j]]) {
    				station=imageStations[imageSets[i][j]];
    				int channel=imageChannels[imageSets[i][j]];
    				if (l > 0) sImgList+=", ";
    				sImgList+=imageSets[i][j]+" ["+channel+"] ";
    				l++;
    	    	}
				gd.addCheckbox    (i+" ("+sImgList+") === Station_"+(station+1), selectedSets[i]);
    		}
    		if (endIndex<enabledSets.length){
         		gd.enableYesNoCancel("Done", "Next");
    		}
			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
//			Arrays.fill(selection, false);
    		for (int i=startIndex;i<endIndex;i++) if (enabledSets[i]){
    			selectedSets[i]=gd.getNextBoolean();
    			for (int j=0;j<imageSets[i].length;j++) {
    				selection[imageSets[i][j]]=selectedSets[i];
    			}
    		}
			if (gd.wasOKed()) return true;
			return selectImageSets(
	    			selection,
	    			allImages,
	    			endIndex,
	    			perPage);
    	}


    	/**
    	 * Manage image selection for the current series
    	 * @param numSeries series number to manage
    	 * @return -1 - cancel, 0 - repeat again, 1 -  Done,
    	 */
    	public int manageSelection(
    			int numSeries){
//    		String [] firstOperand= {"Current","Inverted current","None","All"};
//    		String [] SecondOperand={"Selection","Inverted selection","None","All"};
//    		String [] operation={"And","Or"};
    		String [] actions={
    				"Replace current with new selection",
    				"Add selection to current",
    				"Add inverted selection to current",
    				"And selection with current",
    				"Remove selection from current",
    				"Invert current selection, disregard other settings"};
    		String [] selectionType={
    				"Select from all enabled images",                 // 0
    				"Select from new enabled images",                 // 1
    				"Select from images with estimated orientation",  // 2
    				"Individual images, start from empty selection",  // 3
    				"Individual images, start from current selection",// 4
    				"Image sets, start from empty selection",         // 5
    				"Image sets, start from current selection"        // 6
    				};
    		int numStations=this.distortionCalibrationData.eyesisCameraParameters.getNumStations();
    		int numChannels=this.distortionCalibrationData.eyesisCameraParameters.getNumChannels(0); // for station 0
    		boolean [] selected=          this.selectedImages[numSeries];
    		boolean [] enabled=           this.distortionCalibrationData.selectEnabled();
            boolean [] newEnabled=        this.distortionCalibrationData.selectNewEnabled ();
            boolean [] estimated=         this.distortionCalibrationData.selectEstimated(true); //boolean enabledOnly);
            boolean [] estimatedAll=      this.distortionCalibrationData.selectEstimated(false); //boolean enabledOnly);

    		if (selected.length!=enabled.length){
    			System.out.println("WARNING: manageSelection(): lengths (strategy and images) mismatch selected.length="+selected.length+" available images: "+enabled.length);
   				boolean [] newSelection=new boolean[enabled.length];
				for (int i=0;i<newSelection.length;i++) newSelection[i]=(i<enabled.length)?enabled[i]:false;
				selected=newSelection;
				this.selectedImages[numSeries]=selected;
    		}
    		int numSelected=0;
    		int [] numSelectedPerStation=new int [numStations];
    		Arrays.fill(numSelectedPerStation, 0);
    		int numEnabled  = 0; // this.distortionCalibrationData.getNumEnabled();
    		int [] numEnabledPerStation=new int [numStations];
    		Arrays.fill(numEnabledPerStation, 0);
    		int [] totalPerStation=new int [numStations];
    		Arrays.fill(totalPerStation, 0);
    		int [] imageStations=this.distortionCalibrationData.getStations();
    		int [] imageChannels=this.distortionCalibrationData.getChannels();

    		int numNewEnabled=0;
    		int [] numNewEnabledPerStation=new int [numStations];
    		Arrays.fill(numNewEnabledPerStation, 0);
    		int numNewEnabledSelected=0;
    		int [] numNewEnabledSelectedPerStation=new int [numStations];
    		Arrays.fill(numNewEnabledSelectedPerStation, 0);

    		int numEstimatedSelected=0;
    		int [] numEstimatedSelectedPerStation=new int [numStations];
    		Arrays.fill(numEstimatedSelectedPerStation, 0);

    		int numEstimated=0;
    		int [] numEstimatedPerStation=new int [numStations];
    		Arrays.fill(numEstimatedPerStation, 0);

    		int numEstimatedAll=0;
    		int [] numEstimatedAllPerStation=new int [numStations];
    		Arrays.fill(numEstimatedAllPerStation, 0);

    		int [] matchedPointers=this.distortionCalibrationData.getMatchedPointers();
    		int [] matchedPointersIndex=arrayToSortedInt(matchedPointers);
    		int [] hintedMatch=this.distortionCalibrationData.getHintedMatch();
    		int [] hintedMatchIndex=arrayToSortedInt(hintedMatch);
    		Map <Integer,Integer> mapMP=new HashMap<Integer,Integer>();
    		Map <Integer,Integer> mapHM=new HashMap<Integer,Integer>();
    		for (Integer index=0;index<matchedPointersIndex.length;index++) mapMP.put(new Integer(matchedPointersIndex[index]),index);
    		for (Integer index=0;index<hintedMatchIndex.length;index++) mapHM.put(new Integer(hintedMatchIndex[index]),index);

    		int []   numMatchedPointers=        new int [matchedPointersIndex.length];
    		int []   numMatchedPointersSelected=new int [matchedPointersIndex.length];
    		int []   numMatchedPointersEnabled= new int [matchedPointersIndex.length];
    		int [][] numMatchedPointersPerStation=        new int [matchedPointersIndex.length][];
    		int [][] numMatchedPointersSelectedPerStation=new int [matchedPointersIndex.length][];
    		int [][] numMatchedPointersEnabledPerStation= new int [matchedPointersIndex.length][];
    		for (int n=0;n<numMatchedPointers.length;n++){
    			numMatchedPointers[n]=0;
    			numMatchedPointersSelected[n]=0;
    			numMatchedPointersEnabled[n]=0;
    			numMatchedPointersPerStation[n]=new int [numStations];
    			numMatchedPointersSelectedPerStation[n]=new int [numStations];
    			numMatchedPointersEnabledPerStation[n]=new int [numStations];
        		Arrays.fill(numMatchedPointersPerStation[n], 0);
        		Arrays.fill(numMatchedPointersSelectedPerStation[n], 0);
        		Arrays.fill(numMatchedPointersEnabledPerStation[n], 0);
    		}

    		int []   numHintedMatch=        new int [matchedPointersIndex.length];
    		int []   numHintedMatchSelected=new int [matchedPointersIndex.length];
    		int []   numHintedMatchEnabled= new int [matchedPointersIndex.length];
    		int [][] numHintedMatchPerStation=        new int [matchedPointersIndex.length][];
    		int [][] numHintedMatchSelectedPerStation=new int [matchedPointersIndex.length][];
    		int [][] numHintedMatchEnabledPerStation= new int [matchedPointersIndex.length][];
    		for (int n=0;n<numHintedMatch.length;n++){
    			numHintedMatch[n]=0;
    			numHintedMatchSelected[n]=0;
    			numHintedMatchEnabled[n]=0;
    			numHintedMatchPerStation[n]=new int [numStations];
    			numHintedMatchSelectedPerStation[n]=new int [numStations];
    			numHintedMatchEnabledPerStation[n]=new int [numStations];
        		Arrays.fill(numHintedMatchPerStation[n], 0);
        		Arrays.fill(numHintedMatchSelectedPerStation[n], 0);
        		Arrays.fill(numHintedMatchEnabledPerStation[n], 0);
    		}


    		for (int i=0;i<selected.length;i++) if (imageStations[i]>=0){
    			int mpi=mapMP.get(new Integer(matchedPointers[i]));
    			int hmi=mapHM.get(new Integer(hintedMatch[i]));
    			if (enabled[i]) {
    				numEnabledPerStation[imageStations[i]]++;
    				numEnabled++;
        			if (selected[i]) {
        				numSelectedPerStation[imageStations[i]]++;
        				numSelected++;
        				numMatchedPointersSelectedPerStation[mpi][imageStations[i]]++;
        				numMatchedPointersSelected[mpi]++;

        				numHintedMatchSelectedPerStation[hmi][imageStations[i]]++;
        				numHintedMatchSelected[hmi]++;

        			}
        			if (newEnabled[i]){
        				numNewEnabledPerStation[imageStations[i]]++;
        				numNewEnabled++;
            			if (selected[i]) {
            				numNewEnabledSelectedPerStation[imageStations[i]]++;
            				numNewEnabledSelected++;
            			}
        			}
        			if (estimated[i]){
        				numEstimatedPerStation[imageStations[i]]++;
        				numEstimated++;
            			if (selected[i]) {
            				numEstimatedSelectedPerStation[imageStations[i]]++;
            				numEstimatedSelected++;
            			}
        			}
        			numMatchedPointersEnabledPerStation[mpi][imageStations[i]]++;
    				numMatchedPointersEnabled[mpi]++;

        			numHintedMatchEnabledPerStation[hmi][imageStations[i]]++;
    				numHintedMatchEnabled[hmi]++;


    			}
    			if (estimatedAll[i]){
    				numEstimatedAllPerStation[imageStations[i]]++;
    				numEstimatedAll++;
    			}
    			numMatchedPointersPerStation[mpi][imageStations[i]]++;
				numMatchedPointers[mpi]++;

				numHintedMatchPerStation[hmi][imageStations[i]]++;
				numHintedMatch[hmi]++;

				totalPerStation[imageStations[i]]++;

    		}
    		String sAvailable="["+numSelected+"/"+numEnabled+"/"+selected.length+"]   ";
    		String sNewEnabled="["+numNewEnabledSelected+"/"+numNewEnabled+"]   ";
    		String sEstimated="["+numEstimatedSelected+"/"+numEstimated+"/"+numEstimatedAll+"]   ";
    		String [] sMatchedPointers=new String [matchedPointersIndex.length];
    		String [] sHintedMatch=new String [hintedMatchIndex.length];
    		for (int n=0;n<sMatchedPointers.length;n++){
    			sMatchedPointers[n]="["+numMatchedPointersSelected[n]+"/"+numMatchedPointersEnabled[n]+"/"+numMatchedPointers[n]+"]   ";
    		}
    		for (int n=0;n<sHintedMatch.length;n++){
    			sHintedMatch[n]="["+numHintedMatchSelected[n]+" / "+numHintedMatchEnabled[n]+" / "+numHintedMatch[n]+"]   ";
    		}
    		if (numStations>1) for (int i=0;i<numStations;i++) {
    			sAvailable+= " station_"+(i+1)+": ["+numSelectedPerStation[i]+" / "+numEnabledPerStation[i]+" / "+totalPerStation[i]+"]";
    			sNewEnabled+=" station_"+(i+1)+": ["+numNewEnabledSelectedPerStation[i]+" / "+numNewEnabledPerStation[i]+"]";
    			sEstimated+=" station_"+(i+1)+": ["+numEstimatedSelectedPerStation[i]+" / "+numEstimatedPerStation[i]+" / "+numEstimatedAllPerStation[i]+"]";
    			for (int n=0;n<sMatchedPointers.length;n++){
    				sMatchedPointers[n]+=" station_"+(i+1)+": ["+numMatchedPointersSelectedPerStation[n][i]+" / "+
    						numMatchedPointersEnabledPerStation[n][i]+" / "+numMatchedPointersPerStation[n][i]+"]";
    			}
    			for (int n=0;n<sHintedMatch.length;n++){
    				sHintedMatch[n]+=" station_"+(i+1)+": ["+numHintedMatchSelectedPerStation[n][i]+" / "+
    						numHintedMatchEnabledPerStation[n][i]+" / "+numHintedMatchPerStation[n][i]+"]";
    			}

    		}

    		int operIndex=0,selectionTypeIndex=0;
//			boolean selectEstimated=false;
//			boolean selectNewEnabled=false;
    		boolean [] requiredMatchedPointers=new boolean [matchedPointersIndex.length];
    		boolean [] requiredHintedMatch=new boolean [hintedMatchIndex.length];
    		boolean [] requiredStations=new boolean [numStations];
    		boolean [] requiredChannels=new boolean [numChannels];
    		Arrays.fill(requiredMatchedPointers,true);
    		Arrays.fill(requiredHintedMatch,true);
    		Arrays.fill(requiredStations,true);
    		Arrays.fill(requiredChannels,true);
    		selectionType[0]+=" ("+numEnabled+")";
    		selectionType[1]+=" ("+numNewEnabled+")";
    		selectionType[2]+=" ("+numEstimated+")";
    		if (this.debugLevel>0){
    			System.out.println("Image statistics for series "+numSeries+": [currently selected/enabled/total]");
    			System.out.println("Grid images:"+sAvailable);
    			System.out.println("New enabled images: "+sNewEnabled);
    			System.out.println("Estimated orientation: "+sEstimated);
        		for (int n=0;n<sMatchedPointers.length;n++) System.out.println("Images with "+matchedPointersIndex[n]+" pointers: "+sMatchedPointers[n]);
        		for (int n=0;n<sHintedMatch.length;n++) System.out.println("Images with hinted match state=\""+hintedMatchIndex[n]+"\": "+sHintedMatch[n]);
        		System.out.println();
    		}

    		GenericDialog gd=new GenericDialog("Manage image selection for series "+numSeries);
    		gd.addMessage("Image statistics: [currently selected/enabled/total]");
    		gd.addMessage("Grid images:"+sAvailable);
    		gd.addMessage("New enabled images: "+sNewEnabled);
    		gd.addMessage("Estimated orientation: "+sEstimated);
    		for (int n=0;n<sMatchedPointers.length;n++) gd.addMessage("Images with "+matchedPointersIndex[n]+" pointers: "+sMatchedPointers[n]);
    		for (int n=0;n<sHintedMatch.length;n++) gd.addMessage("Images with hinted match state=\""+hintedMatchIndex[n]+"\": "+sHintedMatch[n]);
			gd.addChoice("Operation on selection", actions, actions[operIndex]);
			gd.addChoice("Selection type", selectionType, selectionType[selectionTypeIndex]);
    		//selectionType
//			gd.addCheckbox("Select images with estimated orientation", selectEstimated);
//			gd.addCheckbox("Select new enabled images", selectNewEnabled);
    		gd.addMessage("=== Filter selection by the number of matched laser pointers ===");
    		for (int i=0;i<matchedPointersIndex.length;i++){
    			gd.addCheckbox("Select images with "+matchedPointersIndex[i]+" pointers", requiredMatchedPointers[i]);
    		}
    		gd.addMessage("=== Filter selection by the hinted match state (-1 - none, 1 - orientation, 2 - position and orientation) ===");
    		for (int i=0;i<hintedMatchIndex.length;i++){
    			gd.addCheckbox("Select images hintedMatch="+hintedMatchIndex[i], requiredHintedMatch[i]);
    		}
    		gd.addMessage("=== Limit selection by the station ===");
    		for (int i=0;i<requiredStations.length;i++){
    			gd.addCheckbox("Select station "+(i+1), requiredStations[i]);
    		}
    		gd.addMessage("=== Limit selection by the channel ===");
    		for (int i=0;i<requiredChannels.length;i++){
    			gd.addCheckbox("Select channel "+i, requiredChannels[i]);
    		}
    		gd.enableYesNoCancel("OK", "More");
			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return -1;
			boolean more=!gd.wasOKed();
			operIndex=gd.getNextChoiceIndex();
			selectionTypeIndex=gd.getNextChoiceIndex(); // TODO:Implement!
//			selectEstimated=gd.getNextBoolean();
//			selectNewEnabled=gd.getNextBoolean();
    		for (int i=0;i<matchedPointersIndex.length;i++) requiredMatchedPointers[i]=gd.getNextBoolean();
    		for (int i=0;i<hintedMatchIndex.length;i++)     requiredHintedMatch[i]=    gd.getNextBoolean();
    		for (int i=0;i<requiredStations.length;i++)     requiredStations[i]=gd.getNextBoolean();
    		for (int i=0;i<requiredChannels.length;i++)     requiredChannels[i]=gd.getNextBoolean();
    		boolean [] selection=new boolean [enabled.length];
    		Arrays.fill(selection,false);
    		switch (selectionTypeIndex){
    		case 0:
    			selection=enabled.clone();
    			break;
    		case 1:
    			selection=newEnabled.clone();
    			break;
    		case 2:
    			selection=estimated.clone();
    			break;
    		case 4: // start from current selection
    			selection=selected.clone();
    		case 3: // start from empty selection
    			if (!selectIndividualImages(
	    			selection,
	    			false, // allImages,
	    			0, // star iIndex
	    			500)) return -1; //perPage))
    		case 6: // start from current selection
    			selection=selected.clone();
    			break;
    		case 5: // start from empty selection
    			if (!selectImageSets(
    	    			selection,
    	    			false, // allImages,
    	    			0, // star iIndex
    	    			500)) return -1; //perPage))
    			break;
    		}
			for (int i=0;i<selection.length;i++) selection[i] &= requiredMatchedPointers[mapMP.get(new Integer(matchedPointers[i]))];
			for (int i=0;i<selection.length;i++) selection[i] &= requiredHintedMatch[mapHM.get(new Integer(hintedMatch[i]))];

			for (int i=0;i<selection.length;i++) if (imageStations[i]>=0) selection[i] &= requiredStations[imageStations[i]];
			else selection[i] = false;

			for (int i=0;i<selection.length;i++) if (imageChannels[i]>=0) selection[i] &= requiredChannels[imageChannels[i]];
			else selection[i] = false;

			// now combine new/old selections
			switch (operIndex){
			case 0: // keep new selection
				break;
			case 1:
				for (int i=0;i<selection.length;i++) selection[i] |= selected[i]; // OR
				break;
			case 2:
				for (int i=0;i<selection.length;i++) selection[i] = selected[i] | !selection[i]; // OR-NOT
				break;
			case 3:
				for (int i=0;i<selection.length;i++) selection[i] = selected[i] && selection[i]; // AND
				break;
			case 4:
				for (int i=0;i<selection.length;i++) selection[i] = selected[i] && !selection[i]; // AND-NOT
				break;
			case 5:
				for (int i=0;i<selection.length;i++) selection[i]= !selected[i];
				break;
			}
			// Remove disabled
			for (int i=0;i<selection.length;i++) selection[i] &= enabled[i];
			// replace current selection
			this.selectedImages[numSeries]=selection;
			return more?0:1;
    	}
       	/**
   	 *
   	 * @param numSeries Number of series to edit
   	 * @param useImages Select images for this series
   	 * @param fromToImages - limit number of checkboxes, otherwise window does not show the bottom ones
   	 * @param useParameters Select parameters for this series
   	 * @param askNextSeries Ask for next series number
   	 * @param zeroAndOther use 2 channels 0 and "other", propagate settings for channel 1 to all the rest
   	 * @return -2 - cancel, -1, done, otherwise - number of step to edit
   	public int selectStrategyStep(
   			int numSeries,
   			boolean useImages,
				int [] fromToImages,
   			boolean allImages,
   			boolean useParameters,
   			boolean askLambdas,
   			boolean askNextSeries,
   			boolean zeroAndOther
   			){
			GenericDialog gd = new GenericDialog("Fitting Strategy Step Configuration, step "+numSeries+
					                             " number of enabled images="+this.distortionCalibrationData.getNumEnabled());
			gd.addStringField  ("Comment",   this.strategyComment[numSeries],80);

   	 */

    	public void listStrategies() {
    		if (this.selectedImages == null) {
    			String msg="No strategies defined";
    			IJ.showMessage("Error",msg);
    			return;
    		}
    		int fromSer = 0;
    		int toSer = this.selectedImages.length;
			boolean zeroAndOther=    true;
    		GenericDialog gd = new GenericDialog("List Fitting Strategies");
    		gd.addNumericField("First number of series to list", fromSer, 0);
    		gd.addNumericField("Last number of series to list", toSer - 1, 0);
    		gd.addCheckbox    ("Use only channel 0 and \"all other channels\"",zeroAndOther);
    		gd.showDialog();
    		if (gd.wasCanceled()) return;
    		fromSer=          (int) gd.getNextNumber();
    		toSer=            (int) gd.getNextNumber();
    		zeroAndOther=           gd.getNextBoolean();
    		listStrategies(fromSer, toSer, zeroAndOther);
    	}

    	public void listStrategies(
        		int     fromSer,
        		int     toSer,
       			boolean zeroAndOther
    			) {
    		if (this.selectedImages == null) {
    			String msg="No strategies defined";
    			IJ.showMessage("Error",msg);
    			return;
    		}
    		int [] choice_offsets = new int [this.parameterEnable.length];
    		//
    	    String header="Strategy #\tType";
    	    StringBuffer sb = new StringBuffer();
    		for (int i=fromSer;i<toSer;i++)  if (isSeriesValid(i)) header+="\t"+i;
			sb.append("Comment\t");
			for (int i=fromSer;i<toSer;i++) if (isSeriesValid(i)) sb.append("\t"+this.strategyComment[i]);sb.append("\n");
			sb.append("Number of selected images\t");
			for (int i=fromSer;i<toSer;i++) if (isSeriesValid(i)) {
				int ns = 0;
				for (int j=0; j<this.selectedImages[i].length;j++) if (this.selectedImages[i][j]) ns++;
				sb.append("\t"+ns);
			}
			sb.append("\n");
			// parameters
			for (int ipar =0; ipar<this.parameterEnable.length;ipar++) if (this.parameterEnable[ipar] &&
					(!zeroAndOther || (this.parameterList[ipar][0] <= 1) || (this.parameterList[ipar][0] ==24))){ // in "zeroAndOther" mode do not show other subcameras
				int parIndex=this.parameterList[ipar][1];
				int subCam=this.parameterList[ipar][0];
				boolean isSub=this.distortionCalibrationData.isSubcameraParameter(parIndex);
				String sChn=(zeroAndOther && (subCam>=1)&& (subCam<24))?"-head-other":
					((zeroAndOther && (subCam>=24))?"-bottom":("-"+subCam));
				boolean noWeak=!this.distortionCalibrationData.eyesisCameraParameters.isExtrinsic(parIndex);
				boolean isTilt=this.distortionCalibrationData.eyesisCameraParameters.isTilt(parIndex);
				// for non-subcamera or subcam 0 - no extra choices. For "zeroAndOther" - only "0". For other subCam - all less than this
				int subMaxNum = (this.distortionCalibrationData.isSubcameraParameter(parIndex) && (subCam > 0)) ? (zeroAndOther? 0 : (subCam-1)):-1;
				String [] commonChoices = (isTilt?this.definedModesTiltEq:(noWeak?this.definedModesNoWeak:this.definedModes));
				String [] theseChoices;
				String thisChoice = "";
				choice_offsets[ipar] = 0;

				sb.append(this.distortionCalibrationData.getParameterName(parIndex)+"\t"+(isSub?(" sub"+sChn):"com "));
				for (int numSeries=fromSer;numSeries<toSer;numSeries++) if (isSeriesValid(numSeries)) {
					if (subMaxNum < 0){
//						theseChoices = commonChoices;
						thisChoice = this.definedModesAll[this.parameterMode[numSeries][ipar]];
					} else { // only can happen if (!isTilt) && (noWeak)
						choice_offsets[ipar] = commonChoices.length;
						theseChoices = new String [commonChoices.length + subMaxNum + 1];
						for (int ch = 0; ch < commonChoices.length; ch++) theseChoices[ch] = commonChoices[ch];
						for (int ch = 0; ch <= subMaxNum; ch++)  theseChoices[ch + commonChoices.length] = "same as "+ch;
						// choice index for "same as ..." starts with  this.definedModesAll.length, but in the listbox index is lower
						int indx = this.parameterMode[numSeries][ipar];
						if (indx >= this.definedModesAll.length) {
							indx -= (this.definedModesAll.length - choice_offsets[ipar]);
						}
						thisChoice = theseChoices[indx];
//						System.out.println("selectStrategyStep(): this.parameterMode["+numSeries+"]["+ipar+"]=" + this.parameterMode[numSeries][ipar]+" indx = "+indx);
					}
					if (thisChoice.equals("fixed")) thisChoice = "-";
					sb.append("\t"+thisChoice);
				}
				sb.append("\n");

			}


		    new TextWindow("Strategies_List", header, sb.toString(), 800,1000);
			System.out.println("** Do not save as .csv - it will replace tabs with commas! ***");
    	}


    	/**
    	 *
    	 * @param numSeries Number of series to edit
    	 * @param useImages Select images for this series
    	 * @param fromToImages - limit number of checkboxes, otherwise window does not show the bottom ones
    	 * @param useParameters Select parameters for this series
    	 * @param askNextSeries Ask for next series number
    	 * @param zeroAndOther use 2 channels 0 and "other", propagate settings for channel 1 to all the rest
    	 * @return -2 - cancel, -1, done, otherwise - number of step to edit
    	 */
    	public int selectStrategyStep(
    			int numSeries,
    			boolean useImages,
				int [] fromToImages,
    			boolean allImages,
    			boolean useParameters,
    			boolean askLambdas,
    			boolean askNextSeries,
    			boolean zeroAndOther
    			){
    		boolean showDirectMap=false;
    		boolean showReverseMap=false;
    		boolean showAdvancedImageSelection=false;
    		int [] choice_offsets = new int [this.parameterEnable.length];

    		// if current series is not valid (probably just started a new one) - look for the last valid (if any)
    		// and copy it;
    		int numEstimated=this.distortionCalibrationData.getNumberOfEstimated(true); //(boolean enabledOnly
    		int [] numEstimatedPerStation=this.distortionCalibrationData.getNumberOfEstimatedPerStation(true); //(boolean enabledOnly
    		String sNumEstimatedPerStation="";
    		for (int i=0;i<numEstimatedPerStation.length;i++){
    			if (i>0) sNumEstimatedPerStation+=", ";
    			sNumEstimatedPerStation+=numEstimatedPerStation[i];
    		}
    		int numNewEnabled=this.distortionCalibrationData.getNumNewEnabled();
    		int [] numNewEnabledPerStation=this.distortionCalibrationData.getNumNewEnabledPerStation();
    		String sNumNewEnabledPerStation="";
    		for (int i=0;i<numNewEnabledPerStation.length;i++){
    			if (i>0) sNumNewEnabledPerStation+=", ";
    			sNumNewEnabledPerStation+=numNewEnabledPerStation[i];
    		}
     		if (!isSeriesValid(numSeries)){
    			int sourceSeries= findLastValidSeries(numSeries);
    			if (sourceSeries>=0)  copySeries(sourceSeries, numSeries);
    		}

			GenericDialog gd = new GenericDialog("Fitting Strategy Step Configuration, step "+numSeries+
					                             " number of enabled images="+this.distortionCalibrationData.getNumEnabled());
			gd.addStringField  ("Comment",   this.strategyComment[numSeries],80);
			gd.addCheckbox("Advanced image selection (disregard other fields)", showAdvancedImageSelection);
			gd.addCheckbox("Copy all from the series below, ignore all other fields", false);
			gd.addNumericField("Source series to copy from", (numSeries>0)?(numSeries-1):(numSeries+1), 0, 3, "");
			gd.addCheckbox("Remove all (but first) images, reopen dialog", false); // remove all will be invalid, copied from the previous
			gd.addCheckbox("Select all images, reopen dialog", false);
			if (numEstimated>0){
				gd.addMessage("There are "+numEstimated+" ("+sNumEstimatedPerStation+") enabled images that have estimated orientation");
				gd.addCheckbox("Select them and only them (and re-open dialog)", false);
			} else {
				gd.addMessage("There are no enabled images with estimated (from neighbors) orientation");
			}
			if (numNewEnabled>0){
				gd.addMessage("There are "+numNewEnabled+" ("+sNumNewEnabledPerStation+" new enabled images");
				gd.addCheckbox("Select them and only them (and re-open dialog)", false);
			} else {
				gd.addMessage("There are no new enabled images");
			}
			int numStations=this.distortionCalibrationData.eyesisCameraParameters.getNumStations();
			boolean [] constrainByStation=new boolean[numStations];
			for (int i=0;i<constrainByStation.length;i++) constrainByStation[i]=true;
			if (this.distortionCalibrationData.eyesisCameraParameters.numStations>1){
				gd.addMessage("Constrain by stations");
//				gd.addCheckbox("Remove images of unselected stations below", true);
				for (int i=0;i<this.distortionCalibrationData.eyesisCameraParameters.numStations;i++) 	gd.addCheckbox("Station "+i, constrainByStation[i]);
			}

			if (useImages) {
	    		gd.addNumericField("Image selection range, from", fromToImages[0], 0);
	    		gd.addNumericField("Image selection range, up to (including)", fromToImages[1], 0);
				gd.addMessage("Select files to include");
				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++)
					if ((allImages || this.distortionCalibrationData.gIP[i].enabled) && (i>=fromToImages[0]) && (i<=fromToImages[1])){
					int hm=this.distortionCalibrationData.gIP[i].hintedMatch;
					gd.addCheckbox    (i+" - "+(this.distortionCalibrationData.gIP[i].enabled?"":"(disabled) ")+
							IJ.d2s(this.distortionCalibrationData.gIP[i].timestamp,6)+
							": "+this.distortionCalibrationData.gIP[i].channel+
							" matched "+this.distortionCalibrationData.gIP[i].matchedPointers+" pointers"+
							", hinted state: "+((hm<0)?"undefined":((hm==0)?"failed":((hm==1)?"orientation":"orientation and translation"))),
							this.selectedImages[numSeries][i]);
				}
				if (allImages) gd.addCheckbox("Enable selected, disable deselected images", false);

				gd.addNumericField("The 'master' (used for common parameters)", this.masterImages[numSeries], 0);

			}
			if (useParameters) {
//				choice_offsets = new int [this.parameterEnable.length];
				gd.addMessage("Select parameters to fit");
				for (int i =0; i<this.parameterEnable.length;i++) if (this.parameterEnable[i] &&
						(!zeroAndOther || (this.parameterList[i][0] <= 1) || (this.parameterList[i][0] ==24))){ // in "zeroAndOther" mode do not show other subcameras
					int parIndex=this.parameterList[i][1];
					int subCam=this.parameterList[i][0];
					boolean isSub=this.distortionCalibrationData.isSubcameraParameter(parIndex);
					boolean defined=false;
					double min=0.0,max=0.0;
					for (int imgNumber=0;imgNumber<this.distortionCalibrationData.getNumImages(); imgNumber++)
						if (this.selectedImages[numSeries][imgNumber] && this.distortionCalibrationData.gIP[imgNumber].enabled){
							int sub=this.distortionCalibrationData.gIP[imgNumber].channel;
							if (!isSub || (sub==subCam) ||  // global or same subcamera
									(zeroAndOther && (subCam>=1) && (subCam<24) && (sub>=1) && (sub<24)) || // both head "other"
									(zeroAndOther && (subCam>=24) && (sub>=24) )
									) { // both subcameras are "other" subcameras
								double parValue=this.distortionCalibrationData.getParameterValue(imgNumber,parIndex);
								if (!defined) {
//									min=this.distortionCalibrationData.pars[imgNumber][parIndex];
									min=parValue;
									max=min;
									defined=true;
								}
//								if (this.distortionCalibrationData.pars[imgNumber][parIndex]<min) min=this.distortionCalibrationData.pars[imgNumber][parIndex];
//								if (this.distortionCalibrationData.pars[imgNumber][parIndex]>max) max=this.distortionCalibrationData.pars[imgNumber][parIndex];
								if (parValue<min) min=parValue;
								if (parValue>max) max=parValue;
							}
						}
//					System.out.println(i+": "+parIndex+":"+subCam+"defined="+defined+" min="+min+" max="+max);
					// undefined, min, max
					String sValue=(defined)?((min==max)?(min+""):(min+"..."+max)):"undefined";
					String sChn=(zeroAndOther && (subCam>=1)&& (subCam<24))?"-head-other":
						((zeroAndOther && (subCam>=24))?"-bottom":("-"+subCam));
					boolean noWeak=!this.distortionCalibrationData.eyesisCameraParameters.isExtrinsic(parIndex);
					boolean isTilt=this.distortionCalibrationData.eyesisCameraParameters.isTilt(parIndex);
					// for non-subcamera or subcam 0 - no extra choices. For "zeroAndOther" - only "0". For other subCam - all less than this
					int subMaxNum = (this.distortionCalibrationData.isSubcameraParameter(parIndex) && (subCam > 0)) ? (zeroAndOther? 0 : (subCam-1)):-1;
					String [] commonChoices = (isTilt?this.definedModesTiltEq:(noWeak?this.definedModesNoWeak:this.definedModes));
					String [] theseChoices;
					String thisChoice = "";
					choice_offsets[i] = 0;
					if (subMaxNum < 0){
						theseChoices = commonChoices;
						thisChoice = this.definedModesAll[this.parameterMode[numSeries][i]];
					} else { // only can happen if (!isTilt) && (noWeak)
						choice_offsets[i] = commonChoices.length;
						theseChoices = new String [commonChoices.length + subMaxNum + 1];
						for (int ch = 0; ch < commonChoices.length; ch++) theseChoices[ch] = commonChoices[ch];
						for (int ch = 0; ch <= subMaxNum; ch++)  theseChoices[ch + commonChoices.length] = "same as "+ch;
						// choice index for "same as ..." starts with  this.definedModesAll.length, but in the listbox index is lower
						int indx = this.parameterMode[numSeries][i];
						if (indx >= this.definedModesAll.length) {
							indx -= (this.definedModesAll.length - choice_offsets[i]);
						}
						thisChoice = theseChoices[indx];
						System.out.println("selectStrategyStep(): this.parameterMode["+numSeries+"]["+i+"]=" + this.parameterMode[numSeries][i]+" indx = "+indx);
					}
					gd.addChoice( // ArrayIndexOutOfBoundsException: 9
							this.distortionCalibrationData.getParameterName(parIndex)+
							" ("+sValue+" "+
							this.distortionCalibrationData.getParameterUnits(parIndex)+")"+
							(this.distortionCalibrationData.isSubcameraParameter(parIndex)?(" sub"+sChn):"com "),
							theseChoices,
							thisChoice);
				}
			}
			if (askLambdas) {
				gd.addNumericField("Initial lambda for the L-M algorithm", this.lambdas[numSeries], 6,8,"");
				gd.addStringField("Relative decrese in error to error ratio to consider series finished", ""+this.stepDone[numSeries], 8);
				gd.addCheckbox("Stop after this series", this.stopAfterThis[numSeries]);
			}
			if (this.varianceModes!=null) gd.addChoice(
					"Processing of selected parameters variances",
					this.definedVarianceModes,
					this.definedVarianceModes[this.varianceModes[numSeries]]);
			gd.addCheckbox("Edit variances costs", false);
			if (numStations>1){
				int oldZGroupsLength=(this.zGroups[numSeries]!=null)?this.zGroups[numSeries].length:0;
				int [] oldZGroups={};
				if (oldZGroupsLength>0) oldZGroups=this.zGroups[numSeries].clone();
				this.zGroups[numSeries]=new int [numStations];
				int nextZGroup=-1;
				for (int i=0;i<numStations;i++){
					if (i<oldZGroupsLength) {
						this.zGroups[numSeries][i]=oldZGroups[i];
						if (this.zGroups[numSeries][i]>nextZGroup) nextZGroup=this.zGroups[numSeries][i];
					} else {
						this.zGroups[numSeries][i]=++nextZGroup;
					}
				}
				String [] zGroupsChoices=new String [numStations+1];
				zGroupsChoices[0]="Not used";
				for (int i=0;i<numStations;i++) zGroupsChoices[i+1]="Group "+i;
				gd.addMessage ("Select groups of same pattern for different stations (pattern did not move between measurements). Used for pattern correction command");
				for (int i=0;i<numStations;i++){
					int zG=(this.zGroups[numSeries][i]<0)?0:(this.zGroups[numSeries][i]+1);
					gd.addChoice( // ArrayIndexOutOfBoundsException: 9
							"Same target group for station " +i,
							zGroupsChoices,
							zGroupsChoices[zG]); // definedModesAll - includes all others

				}
			} else {
				this.zGroups[numSeries]=new int [1];
				this.zGroups[numSeries][0]=0;
			}

			if (askNextSeries) {
				gd.addCheckbox("Rebuild/Show parameter Map", showDirectMap);
				gd.addCheckbox("Rebuild/Show reverse parameter Map", showReverseMap);
				gd.addNumericField("Next series to edit (<0 - done)", numSeries+1, 0);
			}

			if (askNextSeries) gd.enableYesNoCancel("OK", "Done");
			WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return -2;
			this.strategyComment[numSeries]= gd.getNextString();
			showAdvancedImageSelection=gd.getNextBoolean();
			if (showAdvancedImageSelection){
				int rslt=0;
				while (rslt==0) rslt=manageSelection(numSeries);
//				return (rslt<0)?-2:numSeries;
				return numSeries; // cancel from manageSelection will just exit that mode with no changes
			}
			boolean copyFromPrevious=gd.getNextBoolean();
			int sourceStrategy= (int) gd.getNextNumber();
			boolean removeAllImages=gd.getNextBoolean();
			boolean selectAllImages=gd.getNextBoolean();
			boolean selectEstimated=false;
			if (numEstimated>0) selectEstimated=gd.getNextBoolean();
			boolean selectNewEnabled=false;
			if (numNewEnabled>0) selectNewEnabled=gd.getNextBoolean();
			if (this.distortionCalibrationData.eyesisCameraParameters.numStations>1){
//				boolean removeUnselectedStations=gd.getNextBoolean();
				for (int i=0;i<constrainByStation.length; i++) constrainByStation[i]=gd.getNextBoolean();
			}
			if (selectNewEnabled) {
				this.selectedImages[numSeries]=this.distortionCalibrationData.selectNewEnabled();
				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++){
					this.selectedImages[numSeries][i]&=constrainByStation[this.distortionCalibrationData.gIP[i].getStationNumber()];
				}
				return numSeries; // caller will repeat with the same series
			}

			if (selectEstimated) {
				this.selectedImages[numSeries]=this.distortionCalibrationData.selectEstimated(true); //(boolean enabledOnly
				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++){
					this.selectedImages[numSeries][i]&=constrainByStation[this.distortionCalibrationData.gIP[i].getStationNumber()];
				}
				return numSeries; // caller will repeat with the same series
			}

			if (copyFromPrevious){
    			int sourceSeries= findLastValidSeries(sourceStrategy);
    			if (sourceSeries>=0) {
    				copySeries(sourceSeries, numSeries);
    			} else {
    				System.out.println("Could not copy from invalid series "+sourceSeries);
    			}
				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++){
					this.selectedImages[numSeries][i]&=constrainByStation[this.distortionCalibrationData.gIP[i].getStationNumber()];
				}
				return numSeries; // caller will repeat with the same series
			}
			if (removeAllImages || selectAllImages) {
//
				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++){
//					this.selectedImages[numSeries][i]=false; // invalidate - all, regardless of .enabled
					this.selectedImages[numSeries][i]=selectAllImages || ((i==0) && removeAllImages); // invalidate - all, regardless of .enabled
					this.selectedImages[numSeries][i]&=constrainByStation[this.distortionCalibrationData.gIP[i].getStationNumber()];
				}
				return numSeries; // caller will repeat with the same series
			}
			boolean enableDisableSelected=false;
			if (useImages) {
	    		fromToImages[0]=    (int) gd.getNextNumber();
	    		fromToImages[1]=    (int) gd.getNextNumber();
				for (int i =0; i<this.distortionCalibrationData.getNumImages();i++)
					if ((allImages || this.distortionCalibrationData.gIP[i].enabled) && (i>=fromToImages[0]) && (i<=fromToImages[1])){
					this.selectedImages[numSeries][i]=gd.getNextBoolean();
				}
				if (allImages) enableDisableSelected=gd.getNextBoolean();
				this.masterImages[numSeries]=(int) gd.getNextNumber();
				if (this.masterImages[numSeries]<0)this.masterImages[numSeries]=0;
				if (this.masterImages[numSeries]>=this.selectedImages[numSeries].length)this.masterImages[numSeries]=this.selectedImages[numSeries].length;
			}
			if (useParameters) {
				int [] lastGroups=null;
				for (int i =0; i<this.parameterEnable.length;i++) if (this.parameterEnable[i] &&
						(!zeroAndOther || (this.parameterList[i][0] <=1) || (this.parameterList[i][0]==24))){ // in "zeroAndOther" mode do not show other subcameras
					//				for (int i =0; i<this.parameterEnable.length;i++) if (this.parameterEnable[i]){
					this.parameterMode[numSeries][i]=gd.getNextChoiceIndex();

// make adjustment for "same as (other lower numbered subcamera)"
					if ((choice_offsets[i] > 0) && (this.parameterMode[numSeries][i] >= choice_offsets[i])){
						System.out.print("selectStrategyStep(): choice_offsets["+i+"]="+choice_offsets[i]+ " this.parameterMode["+numSeries+"]["+i+"]=" + this.parameterMode[numSeries][i]);
						this.parameterMode[numSeries][i] += (this.definedModesAll.length - choice_offsets[i]);
						System.out.println(", corrected=" + this.parameterMode[numSeries][i]);
					}

					if (this.parameterMode[numSeries][i]==this.modeGroup) {
						if (this.parameterGroups[numSeries][i]!=null) lastGroups=this.parameterGroups[numSeries][i]; // default groups
						else if (lastGroups!=null) this.parameterGroups[numSeries][i]=lastGroups.clone(); // may be null
						selectGroups(numSeries,i);
					}
				}
				if (zeroAndOther){
					for (int i =0; i<this.parameterEnable.length;i++) {
						if ((this.parameterList[i][0]>1) && (this.parameterList[i][0]!=24)){ // "other" subchannels - copy from subchannel1
							int refChannel=(this.parameterList[i][0]<24)?1:24;
							int iSub1=getParameterNumber(refChannel, this.parameterList[i][1]);
							if (this.parameterEnable[iSub1]){
								//					System.out.println(	"parameter number="+i+" this.parameterList[i][0]="+this.parameterList[i][0]+" this.parameterList[i][1]="+this.parameterList[i][1]+" iSub1="+iSub1);
								this.parameterMode[numSeries][i]=this.parameterMode[numSeries][iSub1];
								if (this.parameterMode[numSeries][i]==this.modeGroup) { // copy groups from channel 1
									if (this.parameterGroups[numSeries][i]!=null) this.parameterGroups[numSeries][i]=this.parameterGroups[numSeries][iSub1].clone(); // may be null
									else this.parameterGroups[numSeries][i]=null;
								}
							}
						}
					}
				}
			}
			if (askLambdas) {
				this.lambdas[numSeries]=gd.getNextNumber();
				this.stepDone[numSeries]=Double.parseDouble(gd.getNextString());
				this.stopAfterThis[numSeries]=gd.getNextBoolean();
			}
    		if (enableDisableSelected) {
    			this.distortionCalibrationData.enableSelected(this.selectedImages[numSeries]);
    		}

    		if (this.varianceModes!=null) this.varianceModes[numSeries]=gd.getNextChoiceIndex();

    		boolean editVariancesCosts=gd.getNextBoolean();
    		if (editVariancesCosts){
				for (int i =0; i<this.parameterList.length;i++) {
					int parIndex=this.parameterList[i][1];
					if ((this.parameterMode[numSeries][i]==modeWeakCommon) ||
							(this.parameterMode[numSeries][i]==modeWeakStation) ||
							(this.parameterMode[numSeries][i]==modeTiltEqualize)
					){
						if (!this.distortionCalibrationData.eyesisCameraParameters.isExtrinsic(parIndex)){
							System.out.println("BUG: this.parameterMode["+numSeries+"]["+i+"]="+this.parameterMode[numSeries][i]+
									", but this parameter ("+this.distortionCalibrationData.getParameterName(parIndex)+" is not valid for variances");
							continue;
						}
						this.distortionCalibrationData.eyesisCameraParameters.editCostProperties(
								parIndex,
								this.distortionCalibrationData.getParameterName(parIndex),
								this.distortionCalibrationData.getParameterDescription(parIndex),
								this.distortionCalibrationData.getParameterUnits(parIndex));
					}
				}
    		}

			if (numStations>1){
				for (int i=0;i<numStations;i++){
					this.zGroups[numSeries][i]=gd.getNextChoiceIndex()-1;
				}
			}

    		if (!gd.wasOKed()) return -1; // pressed Done (no need to ask for the next number)

			if (askNextSeries) {
				showDirectMap=gd.getNextBoolean();
				showReverseMap=gd.getNextBoolean();
				if (showDirectMap || showReverseMap){
					buildParameterMap(numSeries);
					if (showDirectMap)  showCurrentParameterMap       ("Parameter map");
					if (showReverseMap) showCurrentReverseParameterMap("Reverse parameter map");
				}

				int nextSeries= (int) gd.getNextNumber();
				if (nextSeries<-1) nextSeries=-1;
				return nextSeries;
			}
			return -1;
    	}

    	private boolean selectGroups(
				int numSeries,
				int numPar){
 //   		if (this.debugLevel>1){
 //   			System.out.println("selectGroups("+numSeries+", "+numPar+")");
 //   		}


			int parIndex=this.parameterList[numPar][1];
			int subCam=this.parameterList[numPar][0];
			String name=this.distortionCalibrationData.getParameterName(parIndex)+
			(this.distortionCalibrationData.isSubcameraParameter(parIndex)?(" s"+subCam):"com ");
    		GenericDialog gd = new GenericDialog("Select image groups for "+name);
    		gd.addMessage("Select which images share the same value of "+name);

    		if (this.parameterGroups[numSeries][numPar]==null) {
    			this.parameterGroups[numSeries][numPar]=new int [this.distortionCalibrationData.getNumImages()];
    			for (int i=0;i<this.parameterGroups[numSeries][numPar].length;i++)this.parameterGroups[numSeries][numPar][i]=0;
    		}
    		String [] choices=organizeGroups(this.parameterGroups[numSeries][numPar],false); // preserve group numbers
    		for (int i =0; i<this.distortionCalibrationData.getNumImages();i++) if (this.selectedImages[numSeries][i]){
    			gd.addChoice (i+" - "+IJ.d2s(this.distortionCalibrationData.gIP[i].timestamp,6)+
    					": "+this.distortionCalibrationData.gIP[i].channel,
    					choices,
    					choices[this.parameterGroups[numSeries][numPar][i]]
    					);
    		}
			WindowTools.addScrollBars(gd);

			gd.showDialog();
			if (gd.wasCanceled()) return false;
    		for (int i =0; i<this.distortionCalibrationData.getNumImages();i++) if (this.selectedImages[numSeries][i]){
    			this.parameterGroups[numSeries][numPar][i]=gd.getNextChoiceIndex();
    		}
			return true;

		}
		/**
		 * Organizes list of groups and creates a list of selection choices. If all members fit in the range
		 * of 0 (length-1) and (force==false), group numbers are preserved, otherwise they are renumbered
		 * @param groups array of integers - group numbers
		 * @param force force renumbering groups
		 * @return list of selection choices
		 */
		private String [] organizeGroups(int []groups, boolean force) {
			if ((groups==null) || (groups.length==0)){
				String msg="Cannot organize mull or empty group";
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
			}
			// See if 0.. length-1 groups is not enough to represent all numbers used
			for (int i=0;!force && (i<groups.length);i++) if ((groups[i]<0) || (groups[i]>=(groups.length-1))) force = true;
			if (force){
				int [] tmp=groups.clone();
				for (int i=0;i<groups.length;i++) groups[i]=-1;
				int groupNumber=0;
				int max=tmp[0]; for (int i=0;i<tmp.length;i++) if (max<tmp[i]) max=tmp[i];
				for (boolean organized=false; !organized;){
					int min=max+1;
					for (int i=0;i<tmp.length;i++) if ((groups[i]<0) && (min>tmp[i])) min=tmp[i];
					if (min>max) organized=true;
					else {
						for (int i=0;i<tmp.length;i++) if ((groups[i]<0) && (min == tmp[i])) {
							groups[i]=groupNumber;
						}
						groupNumber++;
					}
				}
			}
			String [] rslt= new String [groups.length];
			for (int i=0;i<rslt.length;i++) rslt[i]="Group "+(i+1);
			return rslt;

		}
    	public boolean selectStrategy(int startSerNumber){
    		int defaultLength=30;
    		boolean selectImages=    false;
    		boolean allImages=       false;
    		boolean selectParameters=true;
    		boolean askNextSeries=   true;
    		boolean askLambdas=      true;
    		boolean askParameterMask=false;
			boolean zeroAndOther=    true;
			int [] fromToImages={0,500};

    		int numSeries=startSerNumber;
    		int oldLength=(this.selectedImages==null)?0:this.selectedImages.length;

    		GenericDialog gd = new GenericDialog("Fitting Strategy Step Configuration");
    		if (oldLength<=0) gd.addNumericField("Number of series in this strategy", defaultLength, 0);
    		gd.addNumericField("Number of series to edit (<0 - none)", numSeries, 0);
//    		gd.addNumericField("Number of series in this strategy", (oldLength>0)?oldLength:1, 0); //
    		gd.addCheckbox    ("Select images",selectImages);
    		gd.addNumericField("Show image checkboxes from", fromToImages[0], 0);
    		gd.addNumericField("Show image checkboxes up to (including)", fromToImages[1], 0);
    		gd.addCheckbox    ("Select from all images (false - only enabled)",allImages);
    		gd.addCheckbox    ("Select parameters",selectParameters);
    		gd.addCheckbox    ("Ask for initial lambda",askLambdas);
    		gd.addCheckbox    ("Ask for parameter mask",askParameterMask);
    		gd.addCheckbox    ("Use only channel 0 and \"all other channels\"",zeroAndOther);
    		gd.addCheckbox    ("Ask for next series",askNextSeries);
    		if (oldLength>0)  gd.addNumericField("Increase number of series in this strategy", oldLength, 0);
   	        gd.enableYesNoCancel("OK", "Done");

    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		int numberOfSeries=0;
    		if (oldLength<=0) numberOfSeries= (int) gd.getNextNumber();
    		numSeries=          (int) gd.getNextNumber();
    		selectImages=             gd.getNextBoolean();
    		fromToImages[0]=    (int) gd.getNextNumber();
    		fromToImages[1]=    (int) gd.getNextNumber();
    		allImages=                gd.getNextBoolean();
    		selectParameters=         gd.getNextBoolean();
    		askLambdas=               gd.getNextBoolean();
    		askParameterMask=         gd.getNextBoolean();
    		zeroAndOther=             gd.getNextBoolean();
    		askNextSeries=            gd.getNextBoolean();
    		if (oldLength>0) numberOfSeries= (int) gd.getNextNumber();
    		if (numberOfSeries!=oldLength)  setLength(numberOfSeries);
    		if (!gd.wasOKed()) return true;
    		if (askParameterMask) setParameterSelectionMask(zeroAndOther);
    		while ((numSeries>=0) && (numSeries<this.selectedImages.length)){
    			numSeries=selectStrategyStep(
    					numSeries,
    					selectImages,
    					fromToImages,
    					allImages,
    					selectParameters,
    					askLambdas,
    					askNextSeries,
    					zeroAndOther);
    		}
    		return true;
    	}

    	public boolean setParameterSelectionMask(boolean zeroAndOther){
    		GenericDialog gd = new GenericDialog("Set Parameter Selection Mask");
    		gd.addMessage("Common parameters and sub-camera 0 parameters");
    		int subCam=0;
    		boolean showParameter=true;
    		for (int i=0;i<this.parameterList.length;i++){
    			if (this.parameterList[i][0]!=subCam){
    				subCam=this.parameterList[i][0];
    				if (zeroAndOther) {
    					showParameter=true;
    					if (subCam==1) gd.addMessage("Other sub-cameras parameters");
    					else if (subCam==24) gd.addMessage("Bottom sub-cameras parameters");
//    					else break;
    					else {
    						showParameter=false;
//    						continue;
    					}
    				} else {
    					gd.addMessage("Sub-camera "+subCam+" parameters");
    				}
    			}
        		if (showParameter) gd.addCheckbox (this.distortionCalibrationData.getParameterName(parameterList[i][1]),this.parameterEnable[i]);
    		}
      	     WindowTools.addScrollBars(gd);
//   	     gd.setBackground(Color.white);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			subCam=0;
			showParameter=true;
			//May add disablong all other subcameras, but try keeping it for later running with zeroAndOther==false
    		for (int i=0;i<this.parameterList.length;i++){
    			if (this.parameterList[i][0]!=subCam){
    				subCam=this.parameterList[i][0];
					showParameter=true;
//					if (subCam>1) break;
					if ((subCam>1) && (subCam!=24))	showParameter=false;
    			}
    			if (showParameter) this.parameterEnable[i]=gd.getNextBoolean();
    		}
    		return true;
    	}
//this.currentSeriesNumber
    	public double getLambda(){
    		return getLambda(this.currentSeriesNumber);
    	}
    	public double getLambda(int numSeries){
    		return this.lambdas[numSeries];

    	}
    	public double getStepDone(){
    		return getStepDone(this.currentSeriesNumber);
    	}
    	public double getStepDone(int numSeries){
    		return this.stepDone[numSeries];
    	}
    }
