/**
** -----------------------------------------------------------------------------**
** EyesisCorrectionParameters.java
**
** Parameter classes for aberration correction for Eyesis4pi
** 
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  EyesisCorrectionParameters.java is free software: you can redistribute it and/or modify
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
import ij.Prefs;
import ij.gui.GenericDialog;

import java.io.File;
import java.util.Properties;


public class EyesisCorrectionParameters {
    public static class CorrectionParameters{
    	public boolean swapSubchannels01=      true; // false; // (false: 0-1-2, true - 1-0-2)
  		public boolean split=                  true;
  		public boolean vignetting=             true;
  		public boolean pixelDefects=           true;
  		public double  pixelDefectsThreshold=  8.0; // normally none with less than 5.0 are stored?
  		public boolean debayer=                true;
  		public boolean showDebayerEnergy =     false;
  		public boolean saveDebayerEnergy =     true;
  		public boolean deconvolve =            true;
  		public boolean combine =               true;
  		public boolean showDenoiseMask =       false;
  		public boolean saveDenoiseMask =       true;
  		public boolean showChromaDenoiseMask = false;
  		public boolean saveChromaDenoiseMask = true;
  		public boolean showNoiseGains =        false;
  		public boolean saveNoiseGains =        false;
  		public boolean colorProc =             true;
  		public boolean blueProc =              true;
  		public boolean toRGB =                 true;
  		public boolean rotate =                true;
  		public boolean crop =                  true;  // crop to the sensor size 
  		public int     equirectangularFormat=     0;  // 0 - 8 bit RGBA, 1 - 16 bit RGBA, 2 (32 int or 16 float!) ?, 3 - 32-bit FP RGBA. only 0, 1 and 3 currently supported
  		public double  outputRangeInt=          0.25;  // 1.0 intensity will be mapped to 65535*0.25
  		public double  outputRangeFP=          255.0; // 1.0 intensity will be saved as 255.0 (in float 32-bit mode)
  		public boolean imageJTags=             false; // encode ImageJ info data to the TIFF output header 
  		
  		public boolean jpeg =                  true;  // convert to RGB and save JPEG (if save is true)
  		public boolean png =                   true;  // use PNG instead of TIFF for 32-bit ARGB
  		public boolean save =                  true;
  		public boolean save16 =                false; // save 16-bit tiff also if the end result is 8 bit 
  		public boolean save32 =                false; // save 32-bit tiff also if the end result is 8 or 16 bit
  		public boolean show =                  false ;
  		public int     JPEG_quality =          95;
  		public double  JPEG_scale =            0.5;
  		public boolean equirectangular=        true;
  		public boolean zcorrect=               true;
  		public boolean saveSettings =          true;

    	public String [] sourcePaths={};
    	public String sourceDirectory="";
    	public String sourcePrefix="";
    	public String sourceSuffix=".tiff"; //".jp4"
    	public int    firstSubCamera=1; // 0 or 1
    	public String sensorDirectory="";
    	public String sensorPrefix="sensor-";
    	public String sensorSuffix=".calib-tiff"; // fixed in PixelMapping
    	
    	public String sharpKernelDirectory="";
    	public String sharpKernelPrefix="sharpKernel-";
    	public String sharpKernelSuffix=".kernel-tiff";
    	public String smoothKernelDirectory="";
    	public String smoothKernelPrefix="smoothKernel-";
    	public String smoothKernelSuffix=".kernel-tiff";
    	public String dctKernelDirectory="";
    	public String dctKernelPrefix="dct-";
    	public String dctSymSuffix=".sym-tiff";
    	public String dctAsymSuffix=".asym-tiff";
    	public String equirectangularDirectory="";
    	public String equirectangularPrefix="";
    	public String equirectangularSuffix=".eqr-tiff";
    	public boolean equirectangularCut=true;
    	public String planeMapPrefix="";
    	public String planeMapSuffix=".plane-proj-tiff";
    	public boolean usePlaneProjection=false;  // 
    	public boolean planeAsJPEG=       true;   // save de-warped image as JPEG (only if equirectangularFormat==0)
//    	public String equirectangularSuffixA="A.eqr-tiff"; // or the roll-over part
    	public String resultsDirectory="";
    	public boolean removeUnusedSensorData=false;
    	public int exposureCorrectionMode=2; // - 0 - none, 1 - absolute, 2 - relative
    	public double referenceExposure=0.0003; // 3/10000 sec, used in absolute mode only
    	public double relativeExposure=0.5; // 0.0 - use shortest (darken), 1.0 - use longest (brighten)
    	
    	public String cltKernelDirectory="";
    	public String cltKernelPrefix="clt-";
    	public String cltSuffix=".clt-tiff";

    	public String x3dDirectory="";   	
    	
    	public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"split",this.split+"");
  			properties.setProperty(prefix+"vignetting",this.vignetting+"");
  			properties.setProperty(prefix+"pixelDefects",this.pixelDefects+"");
  			properties.setProperty(prefix+"pixelDefectsThreshold",this.pixelDefectsThreshold+"");
  			properties.setProperty(prefix+"debayer",this.debayer+"");
  			properties.setProperty(prefix+"showDebayerEnergy",this.showDebayerEnergy+"");
  			properties.setProperty(prefix+"saveDebayerEnergy",this.saveDebayerEnergy+"");
  			properties.setProperty(prefix+"deconvolve",this.deconvolve+"");
  			properties.setProperty(prefix+"combine",this.combine+"");
  			properties.setProperty(prefix+"showDenoiseMask",this.showDenoiseMask+"");
  			properties.setProperty(prefix+"saveDenoiseMask",this.saveDenoiseMask+"");
  			properties.setProperty(prefix+"showChromaDenoiseMask",this.showChromaDenoiseMask+"");
  			properties.setProperty(prefix+"saveChromaDenoiseMask",this.saveChromaDenoiseMask+"");
  			properties.setProperty(prefix+"showNoiseGains",this.showNoiseGains+"");
  			properties.setProperty(prefix+"saveNoiseGains",this.saveNoiseGains+"");
  			properties.setProperty(prefix+"colorProc",this.colorProc+"");
  			properties.setProperty(prefix+"blueProc",this.blueProc+"");
  			properties.setProperty(prefix+"toRGB",this.toRGB+"");
  			properties.setProperty(prefix+"rotate",this.rotate+"");
  			properties.setProperty(prefix+"crop",this.crop+"");
  			properties.setProperty(prefix+"equirectangularFormat",this.equirectangularFormat+"");
  			properties.setProperty(prefix+"outputRangeInt",this.outputRangeInt+"");
  			properties.setProperty(prefix+"outputRangeFP",this.outputRangeFP+"");
  			properties.setProperty(prefix+"imageJTags",this.imageJTags+"");
  			properties.setProperty(prefix+"jpeg",this.jpeg+"");
  			properties.setProperty(prefix+"png",this.png+"");
  			properties.setProperty(prefix+"save",this.save+"");
  			properties.setProperty(prefix+"save16",this.save16+"");
  			properties.setProperty(prefix+"save32",this.save32+"");
  			properties.setProperty(prefix+"show",this.show+"");
  			properties.setProperty(prefix+"JPEG_quality",this.JPEG_quality+"");
  			properties.setProperty(prefix+"JPEG_scale",this.JPEG_scale+"");
  			properties.setProperty(prefix+"equirectangular",this.equirectangular+"");
  			properties.setProperty(prefix+"zcorrect",this.zcorrect+"");
  			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");

    		properties.setProperty(prefix+"sourceDirectory",this.sourceDirectory);
    		properties.setProperty(prefix+"sourcePrefix",this.sourcePrefix);
    		properties.setProperty(prefix+"sourceSuffix",this.sourceSuffix);
    		properties.setProperty(prefix+"firstSubCamera",this.firstSubCamera+"");
    		
    		properties.setProperty(prefix+"sensorDirectory",this.sensorDirectory);
    		properties.setProperty(prefix+"sensorPrefix",this.sensorPrefix);
    		properties.setProperty(prefix+"sensorSuffix",this.sensorSuffix);
    		
    		
    		properties.setProperty(prefix+"sharpKernelDirectory",this.sharpKernelDirectory);
    		properties.setProperty(prefix+"sharpKernelPrefix",this.sharpKernelPrefix);
    		properties.setProperty(prefix+"sharpKernelSuffix",this.sharpKernelSuffix);
    		properties.setProperty(prefix+"smoothKernelDirectory",this.smoothKernelDirectory);
    		properties.setProperty(prefix+"smoothKernelPrefix",this.smoothKernelPrefix);
    		properties.setProperty(prefix+"smoothKernelSuffix",this.smoothKernelSuffix);
    		
    		properties.setProperty(prefix+"dctKernelDirectory",this.dctKernelDirectory);
    		properties.setProperty(prefix+"dctKernelPrefix",this.dctKernelPrefix);
    		properties.setProperty(prefix+"dctSymSuffix",this.dctSymSuffix);
    		properties.setProperty(prefix+"dctAsymSuffix",this.dctAsymSuffix);

    		properties.setProperty(prefix+"equirectangularDirectory",this.equirectangularDirectory);
    		properties.setProperty(prefix+"equirectangularPrefix",this.equirectangularPrefix);
    		properties.setProperty(prefix+"equirectangularSuffix",this.equirectangularSuffix);
    		properties.setProperty(prefix+"equirectangularCut",this.equirectangularCut+"");
    		
    		properties.setProperty(prefix+"planeMapPrefix",this.planeMapPrefix+"");
    		properties.setProperty(prefix+"planeMapSuffix",this.planeMapSuffix+"");
    		properties.setProperty(prefix+"usePlaneProjection",this.usePlaneProjection+"");
    		properties.setProperty(prefix+"planeAsJPEG",this.planeAsJPEG+"");
    		
    		properties.setProperty(prefix+"resultsDirectory",this.resultsDirectory);
    		properties.setProperty(prefix+"removeUnusedSensorData",this.removeUnusedSensorData+"");
    		if (this.sourcePaths!=null) {
        		properties.setProperty(prefix+"sourcePaths",this.sourcePaths.length+"");
        		for (int i=0;i<this.sourcePaths.length;i++){
        			properties.setProperty(prefix+"sourcePath"+i,this.sourcePaths[i]);
        		}
    		}
    		properties.setProperty(prefix+"exposureCorrectionMode",this.exposureCorrectionMode+"");
    		properties.setProperty(prefix+"referenceExposure",     this.referenceExposure+"");
    		properties.setProperty(prefix+"relativeExposure",      this.relativeExposure+"");
    		properties.setProperty(prefix+"swapSubchannels01",     this.swapSubchannels01+"");
    		
    		properties.setProperty(prefix+"cltKernelDirectory",    this.cltKernelDirectory);
    		properties.setProperty(prefix+"cltKernelPrefix",       this.cltKernelPrefix);
    		properties.setProperty(prefix+"cltSuffix",             this.cltSuffix);

    		properties.setProperty(prefix+"x3dDirectory",          this.x3dDirectory);

    	}

    	public void getProperties(String prefix,Properties properties){
  		    if (properties.getProperty(prefix+"split")!=null) this.split=Boolean.parseBoolean(properties.getProperty(prefix+"split"));
  		    if (properties.getProperty(prefix+"vignetting")!=null) this.vignetting=Boolean.parseBoolean(properties.getProperty(prefix+"vignetting"));
  		    if (properties.getProperty(prefix+"pixelDefects")!=null) this.pixelDefects=Boolean.parseBoolean(properties.getProperty(prefix+"pixelDefects"));
  		    if (properties.getProperty(prefix+"pixelDefectsThreshold")!=null) this.pixelDefectsThreshold=Double.parseDouble(properties.getProperty(prefix+"pixelDefectsThreshold"));
  		    if (properties.getProperty(prefix+"debayer")!=null) this.debayer=Boolean.parseBoolean(properties.getProperty(prefix+"debayer"));
  		    if (properties.getProperty(prefix+"showDebayerEnergy")!=null) this.showDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"showDebayerEnergy"));
  		    if (properties.getProperty(prefix+"saveDebayerEnergy")!=null) this.saveDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"saveDebayerEnergy"));
  		    if (properties.getProperty(prefix+"deconvolve")!=null) this.deconvolve=Boolean.parseBoolean(properties.getProperty(prefix+"deconvolve"));
  		    if (properties.getProperty(prefix+"combine")!=null) this.combine=Boolean.parseBoolean(properties.getProperty(prefix+"combine"));
  		    if (properties.getProperty(prefix+"showDenoiseMask")!=null) this.showDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveDenoiseMask")!=null) this.saveDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveDenoiseMask"));
  		    if (properties.getProperty(prefix+"showChromaDenoiseMask")!=null) this.showChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveChromaDenoiseMask")!=null) this.saveChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"showNoiseGains")!=null) this.showNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"showNoiseGains"));
  		    if (properties.getProperty(prefix+"saveNoiseGains")!=null) this.saveNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"saveNoiseGains"));
  		    if (properties.getProperty(prefix+"colorProc")!=null) this.colorProc=Boolean.parseBoolean(properties.getProperty(prefix+"colorProc"));
  		    if (properties.getProperty(prefix+"blueProc")!=null) this.blueProc=Boolean.parseBoolean(properties.getProperty(prefix+"blueProc"));
  		    if (properties.getProperty(prefix+"toRGB")!=null) this.toRGB=Boolean.parseBoolean(properties.getProperty(prefix+"toRGB"));
  		    if (properties.getProperty(prefix+"rotate")!=null) this.rotate=Boolean.parseBoolean(properties.getProperty(prefix+"rotate"));
  		    if (properties.getProperty(prefix+"crop")!=null) this.crop=Boolean.parseBoolean(properties.getProperty(prefix+"crop"));   // crop to the sensor size
  		    if (properties.getProperty(prefix+"equirectangularFormat")!=null) this.equirectangularFormat=Integer.parseInt(properties.getProperty(prefix+"equirectangularFormat"));
  		    if (properties.getProperty(prefix+"outputRangeInt")!=null) this.outputRangeInt=Double.parseDouble(properties.getProperty(prefix+"outputRangeInt"));
  		    if (properties.getProperty(prefix+"outputRangeFP")!=null) this.outputRangeFP=Double.parseDouble(properties.getProperty(prefix+"outputRangeFP"));
  		    if (properties.getProperty(prefix+"imageJTags")!=null) this.imageJTags=Boolean.parseBoolean(properties.getProperty(prefix+"imageJTags"));
  		    if (properties.getProperty(prefix+"jpeg")!=null) this.jpeg=Boolean.parseBoolean(properties.getProperty(prefix+"jpeg"));   // convert to RGB and save jpeg (if save is true)
  		    if (properties.getProperty(prefix+"png")!=null) this.png=Boolean.parseBoolean(properties.getProperty(prefix+"png"));   // convert to RGB and save jpeg (if save is true)
  		    if (properties.getProperty(prefix+"save")!=null) this.save=Boolean.parseBoolean(properties.getProperty(prefix+"save"));
  		    if (properties.getProperty(prefix+"save16")!=null) this.save16=Boolean.parseBoolean(properties.getProperty(prefix+"save16")); // save 16-bit tiff also if the end result is 8 bit 
  		    if (properties.getProperty(prefix+"save32")!=null) this.save32=Boolean.parseBoolean(properties.getProperty(prefix+"save32")); // save 32-bit tiff also if the end result is 8 or 16 bit
  		    if (properties.getProperty(prefix+"show")!=null) this.show=Boolean.parseBoolean(properties.getProperty(prefix+"show"));
  		    if (properties.getProperty(prefix+"JPEG_quality")!=null) this.JPEG_quality=Integer.parseInt(properties.getProperty(prefix+"JPEG_quality"));
  		    if (properties.getProperty(prefix+"JPEG_scale")!=null) this.JPEG_scale=Double.parseDouble(properties.getProperty(prefix+"JPEG_scale"));
  		    if (properties.getProperty(prefix+"equirectangular")!=null) this.equirectangular=Boolean.parseBoolean(properties.getProperty(prefix+"equirectangular"));
  		    if (properties.getProperty(prefix+"zcorrect")!=null) this.zcorrect=Boolean.parseBoolean(properties.getProperty(prefix+"zcorrect"));
  		    if (properties.getProperty(prefix+"saveSettings")!=null) this.saveSettings=Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
			if (properties.getProperty(prefix+"sourceDirectory")!=      null) this.sourceDirectory=properties.getProperty(prefix+"sourceDirectory");
			if (properties.getProperty(prefix+"sourcePrefix")!=         null) this.sourcePrefix=properties.getProperty(prefix+"sourcePrefix");
			if (properties.getProperty(prefix+"sourceSuffix")!=         null) this.sourceSuffix=properties.getProperty(prefix+"sourceSuffix");
  		    if (properties.getProperty(prefix+"firstSubCamera")!=null) this.firstSubCamera=Integer.parseInt(properties.getProperty(prefix+"firstSubCamera"));
			if (properties.getProperty(prefix+"sensorDirectory")!=      null) this.sensorDirectory=properties.getProperty(prefix+"sensorDirectory");
			if (properties.getProperty(prefix+"sensorPrefix")!=         null) this.sensorPrefix=properties.getProperty(prefix+"sensorPrefix");
			if (properties.getProperty(prefix+"sensorSuffix")!=         null) this.sensorSuffix=properties.getProperty(prefix+"sensorSuffix");
			
			if (properties.getProperty(prefix+"sharpKernelDirectory")!= null) this.sharpKernelDirectory=properties.getProperty(prefix+"sharpKernelDirectory");
			if (properties.getProperty(prefix+"sharpKernelPrefix")!=    null) this.sharpKernelPrefix=properties.getProperty(prefix+"sharpKernelPrefix");
			if (properties.getProperty(prefix+"sharpKernelSuffix")!=    null) this.sharpKernelSuffix=properties.getProperty(prefix+"sharpKernelSuffix");
			if (properties.getProperty(prefix+"smoothKernelDirectory")!=null) this.smoothKernelDirectory=properties.getProperty(prefix+"smoothKernelDirectory");
			if (properties.getProperty(prefix+"smoothKernelPrefix")!=   null) this.smoothKernelPrefix=properties.getProperty(prefix+"smoothKernelPrefix");
			if (properties.getProperty(prefix+"smoothKernelSuffix")!=   null) this.smoothKernelSuffix=properties.getProperty(prefix+"smoothKernelSuffix");

			if (properties.getProperty(prefix+"dctKernelDirectory")!=   null) this.dctKernelDirectory=properties.getProperty(prefix+"dctKernelDirectory");
			if (properties.getProperty(prefix+"dctKernelPrefix")!=      null) this.dctKernelPrefix=properties.getProperty(prefix+"dctKernelPrefix");
			if (properties.getProperty(prefix+"dctSymSuffix")!=         null) this.dctSymSuffix=properties.getProperty(prefix+"dctSymSuffix");
			if (properties.getProperty(prefix+"dctAsymSuffix")!=        null) this.dctAsymSuffix=properties.getProperty(prefix+"dctAsymSuffix");
			
			if (properties.getProperty(prefix+"equirectangularDirectory")!=null) this.equirectangularDirectory=properties.getProperty(prefix+"equirectangularDirectory");
			if (properties.getProperty(prefix+"equirectangularPrefix")!=null) this.equirectangularPrefix=properties.getProperty(prefix+"equirectangularPrefix");
			if (properties.getProperty(prefix+"equirectangularSuffix")!=null) this.equirectangularSuffix=properties.getProperty(prefix+"equirectangularSuffix");

			if (properties.getProperty(prefix+"equirectangularCut")!=null)
				this.equirectangularCut=Boolean.parseBoolean((String)properties.getProperty(prefix+"equirectangularCut"));
//			if (properties.getProperty(prefix+"equirectangularSuffixA")!=null) this.equirectangularSuffixA=properties.getProperty(prefix+"equirectangularSuffixA");

			if (properties.getProperty(prefix+"planeMapPrefix")!=null) this.planeMapPrefix=properties.getProperty(prefix+"planeMapPrefix");
			if (properties.getProperty(prefix+"planeMapSuffix")!=null) this.planeMapSuffix=properties.getProperty(prefix+"planeMapSuffix");
			if (properties.getProperty(prefix+"usePlaneProjection")!=null)
				this.usePlaneProjection=Boolean.parseBoolean((String)properties.getProperty(prefix+"usePlaneProjection"));
			if (properties.getProperty(prefix+"planeAsJPEG")!=null)
				this.planeAsJPEG=Boolean.parseBoolean((String)properties.getProperty(prefix+"planeAsJPEG"));
			if (properties.getProperty(prefix+"resultsDirectory")!=     null) this.resultsDirectory=properties.getProperty(prefix+"resultsDirectory");
			if (properties.getProperty(prefix+"removeUnusedSensorData")!= null)
				this.removeUnusedSensorData=Boolean.parseBoolean((String) properties.getProperty(prefix+"removeUnusedSensorData"));
			if (properties.getProperty(prefix+"sourcePaths")!=   null){
				int numFiles=Integer.parseInt(properties.getProperty(prefix+"sourcePaths"));
				this.sourcePaths=new String[numFiles];
				for (int i=0;i<this.sourcePaths.length;i++){
					this.sourcePaths[i]=properties.getProperty(prefix+"sourcePath"+i);
        		}
			}
  		    if (properties.getProperty(prefix+"exposureCorrectionMode")!=null) this.exposureCorrectionMode=Integer.parseInt(properties.getProperty(prefix+"exposureCorrectionMode"));
  		    if (properties.getProperty(prefix+"referenceExposure")     !=null) this.referenceExposure=   Double.parseDouble(properties.getProperty(prefix+"referenceExposure"));
  		    if (properties.getProperty(prefix+"relativeExposure")      !=null) this.relativeExposure=    Double.parseDouble(properties.getProperty(prefix+"relativeExposure"));
  		    if (properties.getProperty(prefix+"swapSubchannels01")!=null) this.swapSubchannels01=Boolean.parseBoolean(properties.getProperty(prefix+"swapSubchannels01"));
  		    
			if (properties.getProperty(prefix+"cltKernelDirectory")!=   null) this.cltKernelDirectory=properties.getProperty(prefix+"cltKernelDirectory");
			if (properties.getProperty(prefix+"cltKernelPrefix")!=      null) this.cltKernelPrefix=properties.getProperty(prefix+"cltKernelPrefix");
			if (properties.getProperty(prefix+"cltSuffix")!=            null) this.cltSuffix=properties.getProperty(prefix+"cltSuffix");
  		    
			if (properties.getProperty(prefix+"x3dDirectory")!=         null) this.x3dDirectory=properties.getProperty(prefix+"x3dDirectory");
  		    
    	}

    	public boolean showDialog(String title) { 
    		GenericDialog gd = new GenericDialog(title);
    		gd.addCheckbox ("Splt into Bayer stack (if false will exit)",       this.split);
    		gd.addCheckbox ("Apply vignetting/color correction to source files",this.vignetting);
    		gd.addCheckbox ("Replace hot/warm/cold pixels with average of neighbors",this.pixelDefects);
    		gd.addNumericField("Pixel difference thershold to consider it \"bad\" on 255.0 scale (0 - use all)", this.pixelDefectsThreshold, 2,6,"8.0");
			String [] choices={"none","absolute","relative"};
			if (this.exposureCorrectionMode<0) this.exposureCorrectionMode=0;
			else if (this.exposureCorrectionMode>=choices.length) this.exposureCorrectionMode=choices.length-1;
			gd.addChoice      ("Exposure correction",choices, choices[this.exposureCorrectionMode]);
    		gd.addNumericField("Reference exposure (effective only in \"absolute\" mode)", 1000.0*this.referenceExposure, 2,6,"ms");
    		gd.addNumericField("Exposure scale (effective only in \"relative\" mode) 0 - darken, 1 - lighten", this.relativeExposure, 3,5,"");
    		gd.addCheckbox ("De-mosaic (if false will exit)",                   this.debayer);
    		gd.addCheckbox ("Show de-mosaic middle-frequency 'energy",          this.showDebayerEnergy);
    		gd.addCheckbox ("Save de-mosaic middle-frequency 'energy",          this.saveDebayerEnergy);
    		gd.addCheckbox ("Sharpen (convolve with calibration kernels)",      this.deconvolve);
    		gd.addCheckbox ("Denoise (convolve with Gaussian in smooth areas)", this.combine);
    		gd.addCheckbox ("Show denoise mask (white - use hi-res, black - low-res)", this.showDenoiseMask);
    		gd.addCheckbox ("Save denoise mask (white - use hi-res, black - low-res)", this.saveDenoiseMask);
    		gd.addCheckbox ("Show kernel noise gains",                          this.showNoiseGains);
    		gd.addCheckbox ("Save kernel noise gains",                          this.saveNoiseGains);
    		gd.addCheckbox ("Convert colors",                                   this.colorProc);
    		gd.addCheckbox ("Fix blue leak",                                    this.blueProc);
    		gd.addCheckbox ("Show chroma denoise mask (white - use hi-res, black - low-res)", this.showChromaDenoiseMask);
    		gd.addCheckbox ("Save chroma denoise mask (white - use hi-res, black - low-res)", this.saveChromaDenoiseMask);
    		gd.addCheckbox ("Rotate result image",                              this.rotate);
    		gd.addCheckbox ("Crop result image to the original size",           this.crop);
			String [] equirectangularFormatChoices={"RGBA 8-bit","RGBA 16-bit","RGBA 32-bit integer","RGBA 32-bit float","ImageJ stack"};
			int [] equirectangularFormats={0,1,2,3,4};
			int equirectangularFormatIndex=0;
			for ( int i=0;i<equirectangularFormats.length;i++) if (equirectangularFormats[i]==this.equirectangularFormat){
				equirectangularFormatIndex=i;
				break;
			}
			gd.addChoice   ("Equirectangular output format",equirectangularFormatChoices, equirectangularFormatChoices[equirectangularFormatIndex]);
    		gd.addNumericField("Map 1.0 intensity to this fraction of the full range 8/16/32-bit integer mode output", 100*this.outputRangeInt, 2,6,"%");
    		gd.addNumericField("Map 1.0 intensity to this value in 32-bit floating point output mode", this.outputRangeFP, 2,6,"");
    		gd.addCheckbox ("Encode ImageJ specific Info metadata to the output file TIFF header", this.imageJTags);
   		
			gd.addCheckbox ("Convert to RGB48",                                 this.toRGB);
    		gd.addCheckbox ("Convert to 8 bit RGB (and save JPEG if save is enabled)", this.jpeg);
    		gd.addCheckbox ("Use PNG instead of TIFF for 32 bit (8 per color) RGBA", this.png);
    		gd.addCheckbox ("Save the result to file system",                   this.save);
    		gd.addCheckbox ("Save 16-bit tiff if the result is 8 bit",          this.save16);
    		gd.addCheckbox    ("Save 32-bit tiff if the result is 8 or 16 bit",    this.save32);
    		gd.addCheckbox    ("Show the result image",                            this.show);
    		gd.addNumericField("JPEG quality (%)",                                 this.JPEG_quality,0);
    		gd.addNumericField("JPEG scale   (%)",                            100* this.JPEG_scale,0);
    		gd.addCheckbox    ("Warp results to equirectangular",                  this.equirectangular);
    		gd.addCheckbox    ("Calculate distances in overlapping areas",         this.zcorrect);
    		gd.addCheckbox    ("Save current settings with results",               this.saveSettings);
    		gd.addStringField ("Source files directory",                           this.sourceDirectory, 60);
    		gd.addCheckbox    ("Select source directory",                          false);
    		gd.addStringField ("Sensor calibration directory",                     this.sensorDirectory, 60);
    		gd.addCheckbox    ("Select sensor calibration directory",              false);

    		gd.addStringField ("Aberration kernels (sharp) directory",             this.sharpKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels (sharp) directory",      false);
    		gd.addStringField ("Aberration kernels (smooth) directory",            this.smoothKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels (smooth) directory",     false);
    		
    		gd.addStringField ("Aberration kernels for DCT directory",             this.dctKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels for DCT directory",      false);

    		gd.addStringField ("Aberration kernels for CLT directory",             this.cltKernelDirectory, 60);
    		gd.addCheckbox    ("Select aberration kernels for CLT directory",      false);
    		
    		gd.addStringField ("x3d output directory",                             this.x3dDirectory, 60);
    		gd.addCheckbox    ("Select x3d output directory",                      false);

    		gd.addStringField("Equirectangular maps directory (may be empty)",     this.equirectangularDirectory, 60);
    		gd.addCheckbox("Select equirectangular maps directory",                false);
    		gd.addStringField("Results directory",                                 this.resultsDirectory, 40);
    		gd.addCheckbox("Select results directory",                             false);
    		gd.addStringField("Source files prefix",                               this.sourcePrefix, 40);
    		gd.addStringField("Source files suffix",                               this.sourceSuffix, 40);
    		gd.addNumericField("First subcamera (in the source filename)",         this.firstSubCamera, 0);
    		
    		gd.addStringField("Sensor files prefix",                               this.sensorPrefix, 40);
    		gd.addStringField("Sensor files suffix",                               this.sensorSuffix, 40);
    		gd.addStringField("Kernel files (sharp) prefix",                       this.sharpKernelPrefix, 40);
    		gd.addStringField("Kernel files (sharp) suffix",                       this.sharpKernelSuffix, 40);
    		gd.addStringField("Kernel files (smooth) prefix",                      this.smoothKernelPrefix, 40);
    		gd.addStringField("Kernel files (smooth) suffix",                      this.smoothKernelSuffix, 40);

    		gd.addStringField("DCT kernel files  prefix",                          this.dctKernelPrefix, 40);
    		gd.addStringField("DCT symmetical kernel files",                       this.dctSymSuffix, 40);
    		gd.addStringField("DCT asymmetrical kernel files suffix",              this.dctAsymSuffix, 40);
    		gd.addStringField("CLT kernel files  prefix",                          this.cltKernelPrefix, 40);
    		gd.addStringField("CLT symmetical kernel files",                       this.cltSuffix, 40);
    		
    		gd.addStringField("Equirectangular maps prefix",     this.equirectangularPrefix, 40);
    		gd.addStringField("Equirectangular maps suffix",     this.equirectangularSuffix, 40);
    		gd.addCheckbox("Cut rolling-over equirectangular images in two", this.equirectangularCut);
    		
    		gd.addStringField("Plane projection map prefix",     this.planeMapPrefix, 40);
    		gd.addStringField("Plane projection map suffix",     this.planeMapSuffix, 40);
    		gd.addCheckbox("Use projection to a common plane instead of the  equirectangular", this.usePlaneProjection);
    		gd.addCheckbox("Save de-warped images as JPEG instead of TIFF",  this.planeAsJPEG);

    		
//    		gd.addStringField("Suffix for the second part of rolled-over equirectangular images",  this.equirectangularSuffixA, 40);
			
    		gd.addCheckbox   ("Remove unused sensor data",       this.removeUnusedSensorData);
    		gd.addCheckbox   ("Swap top and equator images",     this.swapSubchannels01);
    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.split=                  gd.getNextBoolean();
    		this.vignetting=             gd.getNextBoolean();
    		this.pixelDefects=           gd.getNextBoolean();
    		this.pixelDefectsThreshold=  gd.getNextNumber();
			this.exposureCorrectionMode= gd.getNextChoiceIndex();
    		this.referenceExposure=0.001*gd.getNextNumber();
    		this.relativeExposure=       gd.getNextNumber();
    		this.debayer=           gd.getNextBoolean();
    		this.showDebayerEnergy= gd.getNextBoolean();
    		this.saveDebayerEnergy= gd.getNextBoolean();
    		this.deconvolve=        gd.getNextBoolean();
    		this.combine=           gd.getNextBoolean();
    		this.showDenoiseMask=   gd.getNextBoolean();
    		this.saveDenoiseMask=   gd.getNextBoolean();
    		this.showNoiseGains=    gd.getNextBoolean();
    		this.saveNoiseGains=    gd.getNextBoolean();
    		this.colorProc=         gd.getNextBoolean();
    		this.blueProc=         gd.getNextBoolean();
    		this.showChromaDenoiseMask=   gd.getNextBoolean();
    		this.saveChromaDenoiseMask=   gd.getNextBoolean();
    		this.rotate=            gd.getNextBoolean();
    		this.crop=              gd.getNextBoolean();
    		this.equirectangularFormat= equirectangularFormats[gd.getNextChoiceIndex()];
    		this.outputRangeInt=0.01*gd.getNextNumber();
    		this.outputRangeFP=     gd.getNextNumber();
    		this.imageJTags=        gd.getNextBoolean();
    		this.toRGB=             gd.getNextBoolean();
    		this.jpeg=              gd.getNextBoolean();
    		this.png=               gd.getNextBoolean();
    		this.save=              gd.getNextBoolean();
    		this.save16=            gd.getNextBoolean();
    		this.save32=            gd.getNextBoolean();
    		this.show=              gd.getNextBoolean();
    		this.JPEG_quality=(int) gd.getNextNumber();
    		this.JPEG_scale=   0.01*gd.getNextNumber();
    		this.equirectangular=   gd.getNextBoolean();
    		this.zcorrect=          gd.getNextBoolean();
    		this.saveSettings=      gd.getNextBoolean();

    		this.sourceDirectory=        gd.getNextString(); if (gd.getNextBoolean()) selectSourceDirectory(false, false); 
    		this.sensorDirectory=        gd.getNextString(); if (gd.getNextBoolean()) selectSensorDirectory(false, false); 
    		this.sharpKernelDirectory=   gd.getNextString(); if (gd.getNextBoolean()) selectSharpKernelDirectory(false, false); 
    		this.smoothKernelDirectory=  gd.getNextString(); if (gd.getNextBoolean()) selectSmoothKernelDirectory(false, true);
    		this.dctKernelDirectory=     gd.getNextString(); if (gd.getNextBoolean()) selectDCTKernelDirectory(false, true);
    		this.cltKernelDirectory=     gd.getNextString(); if (gd.getNextBoolean()) selectCLTKernelDirectory(false, true);
    		this.x3dDirectory=           gd.getNextString(); if (gd.getNextBoolean()) selectX3dDirectory(false, true);
    		this.equirectangularDirectory=  gd.getNextString(); if (gd.getNextBoolean()) selectEquirectangularDirectory(false, false); 
    		this.resultsDirectory=       gd.getNextString(); if (gd.getNextBoolean()) selectResultsDirectory(false, true); 
    		this.sourcePrefix=           gd.getNextString();
    		this.sourceSuffix=           gd.getNextString();
    		this.firstSubCamera=   (int) gd.getNextNumber();
    		this.sensorPrefix=           gd.getNextString();
    		this.sensorSuffix=           gd.getNextString();
    		this.sharpKernelPrefix=      gd.getNextString();
    		this.sharpKernelSuffix=      gd.getNextString();
    		this.smoothKernelPrefix=     gd.getNextString();
    		this.smoothKernelSuffix=     gd.getNextString();
    		this.dctKernelPrefix=        gd.getNextString();
    		this.dctSymSuffix=           gd.getNextString();
    		this.dctAsymSuffix=          gd.getNextString();
    		this.cltKernelPrefix=        gd.getNextString();
    		this.cltSuffix=              gd.getNextString();
    		this.equirectangularPrefix=  gd.getNextString();
    		this.equirectangularSuffix=  gd.getNextString();
    		this.equirectangularCut=     gd.getNextBoolean();
    		this.planeMapPrefix=         gd.getNextString();
    		this.planeMapSuffix=         gd.getNextString();
    		this.usePlaneProjection=     gd.getNextBoolean();
    		this.planeAsJPEG=            gd.getNextBoolean();
//    		this.equirectangularSuffixA= gd.getNextString();
    		this.removeUnusedSensorData= gd.getNextBoolean();
    		this.swapSubchannels01= gd.getNextBoolean();
    		return true;
    	}

    	
// TODO: extract timestamnp from JP4 or, at least combine movie timestamp+frame into a single filename string
    	public String [] getSourcePaths(){
    		String [] empty={};
    		return (this.sourcePaths!=null)?this.sourcePaths:empty;
    	}

    	public int getChannelFromTiff(String path, String suffix){
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) && 
    				(path.charAt(indexLastDash)!='_') && 
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		return Integer.parseInt(path.substring(indexLastDash+1,indexSuffix));

    	}
    	public String getNameFromTiff(String path, String suffix){
    		int indexSuffix=path.length()-suffix.length();
    		int indexLastDash=indexSuffix-1; // in jp4 it will be underscore, not dash? Or should we change that?
    		while ((indexLastDash>0) &&
    				(indexLastDash>(indexSuffix-3)) && 
    				(path.charAt(indexLastDash)!='_') && 
    				(path.charAt(indexLastDash)!='-')) indexLastDash--;
    		int nameStart=path.lastIndexOf(Prefs.getFileSeparator())+1;
    		return path.substring(nameStart,indexLastDash);
    	}
    	public boolean isJP4(){
			return this.sourceSuffix.equals(".jp4") || this.sourceSuffix.equals(".jp46");
    	}
    	public int getChannelFromSourceTiff(String path){ return getChannelFromTiff(path, this.sourceSuffix);	}
    	public String getNameFromSourceTiff(String path){ return getNameFromTiff(path, this.sourceSuffix);	}
    	
    	public int getChannelFromKernelTiff(String path, int type){return getChannelFromTiff(path, (type==0)?this.sharpKernelSuffix:this.smoothKernelSuffix);}
    	public String getNameFromKernelTiff(String path, int type){return getNameFromTiff(path, (type==0)?this.sharpKernelSuffix:this.smoothKernelSuffix);}
    	
    	public int getChannelFromDCTTiff(String path, int type){return getChannelFromTiff(path, (type==0)?this.dctSymSuffix:this.dctAsymSuffix);}
    	public String getNameFromDCTTiff(String path, int type){return getNameFromTiff(path, (type==0)?this.dctSymSuffix:this.dctAsymSuffix);}
    	
    	public int getChannelFromCLTTiff(String path){return getChannelFromTiff(path, this.cltSuffix);}
    	public String getNameFromCLTTiff(String path){return getNameFromTiff(path, this.cltSuffix);}
    	
    	public boolean selectSourceFiles(boolean allFiles) {
    		return selectSourceFiles(allFiles, 1); // debug level 1 - modify here
    	}
    	public boolean selectSourceFiles(boolean allFiles, int debugLevel) {
    		String [] defaultPaths=this.sourcePaths;
    		if ((defaultPaths==null) || (defaultPaths.length==0)){
    			defaultPaths = new String[1];
    			if ((this.sourceDirectory==null) || (this.sourceDirectory.length()==0)){
    				defaultPaths[0]="";
    			} else {
    				defaultPaths[0]=this.sourceDirectory+Prefs.getFileSeparator();
    			}
    		}
    		String [] extensions={this.sourceSuffix};
			CalibrationFileManagement.MultipleExtensionsFileFilter sourceFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(this.sourcePrefix,extensions,"Source files");
			String [] sourceFiles=null;
    		if (allFiles){
				File dir= new File (this.sourceDirectory);
				if (debugLevel>1) System.out.println("selectSourceFiles, dir="+this.sourceDirectory);
				if (!dir.exists()) {
					String error="Source directory "+this.sourceDirectory+" does not exist.";
					if (debugLevel>1) System.out.println("selectSourceFiles() ERROR:"+error);
					if (debugLevel>1) IJ.showMessage("No files selected");
					return false;
				}
				File [] fileList=dir.listFiles(sourceFilter);
				if (debugLevel>1) System.out.println("Source directory "+this.sourceDirectory+" has "+fileList.length+" files.");
				sourceFiles = new String[fileList.length];
				for (int i=0;i<sourceFiles.length;i++) sourceFiles[i]=fileList[i].getPath();
    		} else {
    			sourceFiles=CalibrationFileManagement.selectFiles(false,
    					"Select Source files, saved as "+ extensions[0],
    					"Select",
    					sourceFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    	       	if ((sourceFiles==null) || (sourceFiles.length==0)) {
					if (debugLevel>1) System.out.println("selectSourceFiles() ERROR: No files selected");
					if (debugLevel>1) IJ.showMessage("No files selected");
            		return false;
            	}
    		}
    		this.sourcePaths=sourceFiles;
    		if ((this.sourcePaths!=null) && (this.sourcePaths.length>0)){
    			this.sourceDirectory=this.sourcePaths[0].substring(0, this.sourcePaths[0].lastIndexOf(Prefs.getFileSeparator()));
    			//sourceNames
    		}
			return true;
    	}

    	public String [] selectSensorFiles(int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		if ((this.sensorDirectory==null) || (this.sensorDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=this.sensorDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.sensorSuffix};
			CalibrationFileManagement.MultipleExtensionsFileFilter sensorFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(this.sensorPrefix,extensions,this.sensorPrefix+"*"+extensions[0]+" Sensor calibration files");
			if (debugLevel>1) System.out.println("selectSensorFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+this.sensorPrefix+"*"+this.sensorSuffix);

			String [] sensorFiles=null;
// try reading all matching files
			File dir= new File (this.sensorDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(sensorFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				sensorFiles=CalibrationFileManagement.selectFiles(false,
    					"Select sensor calibration files",
    					"Select",
    					sensorFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((sensorFiles!=null) && (sensorFiles.length>0)){
    				this.sensorDirectory=sensorFiles[0].substring(0, sensorFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (this.sensorDirectory);
//    				if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
    				fileList=dir.listFiles(sensorFilter);
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
    		
			if (debugLevel>1) System.out.println("Sensor directory "+this.sensorDirectory+" has "+fileList.length+" matching sensor files.");
			sensorFiles = new String[fileList.length];
			for (int i=0;i<sensorFiles.length;i++) sensorFiles[i]=fileList[i].getPath();

			String directory=sensorFiles[0].substring(0, sensorFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=sensorFiles[0].substring(directory.length()+1, sensorFiles[0].length()-extensions[0].length()-2); // all but NN
			this.sensorDirectory=directory;
			this.sensorPrefix=prefix;
		
			return sensorFiles;
    	}

    	public String selectEquirectangularMapFile(
    			int channel,
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String defaultPath="";
    		if ((this.equirectangularDirectory!=null) && (this.equirectangularDirectory.length()>1)){ // empty or "/"
    			defaultPath=this.equirectangularDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={String.format("%02d",channel)+this.equirectangularSuffix}; // looking just for a single map
			CalibrationFileManagement.MultipleExtensionsFileFilter equirectangularFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(this.equirectangularPrefix,extensions,
						this.equirectangularPrefix+"*"+extensions[0]+" Equirectangular map for channel "+channel);
			if (debugLevel>1) System.out.println("selectEquirectangularMapFile("+debugLevel+"): defaultPath="+defaultPath+
					" "+this.equirectangularPrefix+"*"+this.equirectangularSuffix);

			String equirectangularFile=null;
// try reading all matching files
			File dir= new File (this.equirectangularDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(equirectangularFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				equirectangularFile=CalibrationFileManagement.selectFile(false,
    					"Select Equirectangular map file for channel "+channel,
    					"Select",
    					equirectangularFilter,
    					defaultPath); // String [] defaultPaths); //this.sourceDirectory // null
    			if (equirectangularFile!=null) {
    				this.equirectangularDirectory=equirectangularFile.substring(0, equirectangularFile.lastIndexOf(Prefs.getFileSeparator()));
    				this.equirectangularPrefix=equirectangularFile.substring(this.equirectangularDirectory.length()+1, equirectangularFile.length()-extensions[0].length()-2);
    				return equirectangularFile;
    			} else return null;
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			equirectangularFile=fileList[0].getPath();
			this.equirectangularDirectory=equirectangularFile.substring(0, equirectangularFile.lastIndexOf(Prefs.getFileSeparator()));
			this.equirectangularPrefix=equirectangularFile.substring(this.equirectangularDirectory.length()+1, equirectangularFile.length()-extensions[0].length()); // extensions[0] already includes channel
			if (fileList.length>1) {
				String msg = "Multiple files matched, prefix updated to match just the first one - "+equirectangularFile;
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
			}
			if (debugLevel>1) System.out.println("selectEquirectangularMapFile() -> "+ equirectangularFile);
			return equirectangularFile;
    	}

    	public String selectPlaneMapFile(
//    			int channel,
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String defaultPath="";
    		if ((this.equirectangularDirectory!=null) && (this.equirectangularDirectory.length()>1)){ // empty or "/"
    			defaultPath=this.equirectangularDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.planeMapSuffix}; // looking just for a single map
			CalibrationFileManagement.MultipleExtensionsFileFilter planeMapFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(this.planeMapPrefix,extensions,
						this.planeMapPrefix+"*"+extensions[0]+" Plane map (all channels)");
			if (debugLevel>1) System.out.println("selectPlaneMapFile("+debugLevel+"): defaultPath="+defaultPath+
					" "+this.planeMapPrefix+"*"+this.planeMapSuffix);
			String planeMapFile=null;
// try reading all matching files
			File dir= new File (this.equirectangularDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(planeMapFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				planeMapFile=CalibrationFileManagement.selectFile(false,
    					"SelectPlane map file for all channels",
    					"Select",
    					planeMapFilter,
    					defaultPath); // String [] defaultPaths); //this.sourceDirectory // null
    			if (planeMapFile!=null) {
    				this.equirectangularDirectory=planeMapFile.substring(0, planeMapFile.lastIndexOf(Prefs.getFileSeparator()));
    				this.planeMapPrefix=planeMapFile.substring(this.equirectangularDirectory.length()+1, planeMapFile.length()-extensions[0].length()-2);
    				return planeMapFile;
    			} else return null;
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			planeMapFile=fileList[0].getPath();
			this.equirectangularDirectory=planeMapFile.substring(0, planeMapFile.lastIndexOf(Prefs.getFileSeparator()));
			this.planeMapPrefix=planeMapFile.substring(this.equirectangularDirectory.length()+1, planeMapFile.length()-extensions[0].length()); // extensions[0] already includes channel
			if (fileList.length>1) {
				String msg = "Multiple files matched, prefix updated to match just the first one - "+planeMapFile;
				System.out.println("Warning: "+msg);
				IJ.showMessage("Warning",msg);
			}
			if (debugLevel>1) System.out.println("selectPlaneMapFile() -> "+ planeMapFile);
			return planeMapFile;
    	}
    	
    	
    	
    	public String [] selectKernelChannelFiles(
    			int type,  // 0 - sharp, 1 - smooth
    			int numChannels, // number of channels
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String [] kernelFiles= selectKernelFiles(
        			type,  // 0 - sharp, 1 - smooth
        			debugLevel);
    		if (kernelFiles==null) return null;
    		String [] channelPaths=new String[numChannels];
    		for (int i=0;i<channelPaths.length;i++)channelPaths[i]=null;
    		for (int fileNum=0;fileNum<kernelFiles.length;fileNum++){
    			int chn=getChannelFromKernelTiff(kernelFiles[fileNum], type);
    			if ((chn>=0) && (chn<numChannels)){
    				if (channelPaths[chn]==null){ // use first file for channel if there are multiple
    					channelPaths[chn]=kernelFiles[fileNum];
    				} else {
    					if (debugLevel>0) System.out.println("Multiple kernel files for channel "+
    							chn+": "+channelPaths[chn]+" and "+kernelFiles[fileNum]+". Using "+channelPaths[chn]);
    				}
    			}
    		}
    		return channelPaths;
    	}
    	
    	public String [] selectKernelFiles(
    			int type,  // 0 - sharp, 1 - smooth
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		String kernelDirectory=(type==0)?this.sharpKernelDirectory:this.smoothKernelDirectory;
    		if ((kernelDirectory==null) || (kernelDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=kernelDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={(type==0)?this.sharpKernelSuffix:this.smoothKernelSuffix};
    		String  kernelPrefix= (type==0)?this.sharpKernelPrefix:this.smoothKernelPrefix;
			CalibrationFileManagement.MultipleExtensionsFileFilter kernelFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(kernelPrefix,extensions,kernelPrefix+
						"*"+extensions[0]+" "+((type==0)?"Sharp":"Smooth")+" kernel files");
			if (debugLevel>1) System.out.println("selectKernelFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+kernelPrefix+"*"+extensions[0]);

			String [] kernelFiles=null;
// try reading all matching files
			File dir= new File (kernelDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(kernelFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				kernelFiles=CalibrationFileManagement.selectFiles(false,
    					"Select"+((type==0)?"sharp":"smooth")+" kernel files",
    					"Select",
    					kernelFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((kernelFiles!=null) && (kernelFiles.length>0)){
    				kernelDirectory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (kernelDirectory);
//    				if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
    				fileList=dir.listFiles(kernelFilter);
    				if (type==0) this.sharpKernelDirectory= kernelDirectory;
    				else         this.smoothKernelDirectory=kernelDirectory;
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			if (debugLevel>1) System.out.println(((type==0)?"Sharp":"Smooth")+" kernel directory "+kernelDirectory+" has "+fileList.length+" matching files.");
			kernelFiles = new String[fileList.length];
			for (int i=0;i<kernelFiles.length;i++) kernelFiles[i]=fileList[i].getPath();
			String directory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=kernelFiles[0].substring(directory.length()+1, kernelFiles[0].length()-extensions[0].length()-2); // all but NN
			if (type==0) this.sharpKernelDirectory=directory;
			else         this.smoothKernelDirectory=directory;
			if (type==0) this.sharpKernelPrefix=prefix;
			else         this.smoothKernelPrefix=prefix;
			return kernelFiles;
    	}
    	
    	public String [] selectDCTChannelFiles(
    			int numChannels, // number of channels
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String [] kernelFiles= selectDCTFiles(
        			debugLevel);
    		if (kernelFiles==null) return null;
    		String [] channelPaths=new String[numChannels];
    		for (int i=0;i<channelPaths.length;i++)channelPaths[i]=null;
    		for (int fileNum=0;fileNum<kernelFiles.length;fileNum++){
    			int chn=getChannelFromDCTTiff(kernelFiles[fileNum], 0); // 1 for asym files
    			if ((chn>=0) && (chn<numChannels)){
    				if (channelPaths[chn]==null){ // use first file for channel if there are multiple
    					channelPaths[chn]=kernelFiles[fileNum];
    				} else {
    					if (debugLevel>0) System.out.println("Multiple kernel files for channel "+
    							chn+": "+channelPaths[chn]+" and "+kernelFiles[fileNum]+". Using "+channelPaths[chn]);
    				}
    			}
    		}
    		return channelPaths;
    	}
    	
    	public String [] selectDCTFiles(
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		String kernelDirectory=this.dctKernelDirectory;
    		if ((kernelDirectory==null) || (kernelDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=kernelDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.dctSymSuffix};
    		String  kernelPrefix= this.dctKernelPrefix;
			CalibrationFileManagement.MultipleExtensionsFileFilter kernelFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(kernelPrefix,extensions,kernelPrefix+
						"*"+extensions[0]+" DCT symmetrical kernel files");
			if (debugLevel>1) System.out.println("selectKernelFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+kernelPrefix+"*"+extensions[0]);

			String [] kernelFiles=null;
// try reading all matching files
			File dir= new File (kernelDirectory);
//			if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(kernelFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				kernelFiles=CalibrationFileManagement.selectFiles(false,
    					"Select DCT symmetrical kernel files",
    					"Select",
    					kernelFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((kernelFiles!=null) && (kernelFiles.length>0)){
    				kernelDirectory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (kernelDirectory);
//    				if (debugLevel>1) System.out.println("selectSensorFiles, dir="+this.sensorDirectory);
    				fileList=dir.listFiles(kernelFilter);
    				this.dctKernelDirectory= kernelDirectory;
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			if (debugLevel>1) System.out.println("DCT kernel directory "+kernelDirectory+" has "+fileList.length+" matching files.");
			kernelFiles = new String[fileList.length];
			for (int i=0;i<kernelFiles.length;i++) kernelFiles[i]=fileList[i].getPath();
			String directory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=kernelFiles[0].substring(directory.length()+1, kernelFiles[0].length()-extensions[0].length()-2); // all but NN
			this.dctKernelDirectory=directory;
			this.dctKernelPrefix=prefix;
			return kernelFiles;
    	}
    	

    	
    	public String [] selectCLTChannelFiles(
    			int numChannels, // number of channels
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String [] kernelFiles= selectCLTFiles(
        			debugLevel);
    		if (kernelFiles==null) return null;
    		String [] channelPaths=new String[numChannels];
    		for (int i=0;i<channelPaths.length;i++)channelPaths[i]=null;
    		for (int fileNum=0;fileNum<kernelFiles.length;fileNum++){
    			int chn=getChannelFromCLTTiff(kernelFiles[fileNum]);
    			if ((chn>=0) && (chn<numChannels)){
    				if (channelPaths[chn]==null){ // use first file for channel if there are multiple
    					channelPaths[chn]=kernelFiles[fileNum];
    				} else {
    					if (debugLevel>0) System.out.println("Multiple kernel files for channel "+
    							chn+": "+channelPaths[chn]+" and "+kernelFiles[fileNum]+". Using "+channelPaths[chn]);
    				}
    			}
    		}
    		return channelPaths;
    	}
    	
    	public String [] selectCLTFiles(
    			int debugLevel) { // will only open dialog if directory or files are not found
    		String []defaultPaths = new String[1];
    		String kernelDirectory=this.cltKernelDirectory;
    		if ((kernelDirectory==null) || (kernelDirectory.length()<=1)){ // empty or "/"
    			defaultPaths[0]="";
    		} else {
    			defaultPaths[0]=kernelDirectory+Prefs.getFileSeparator();
    		}
    		String [] extensions={this.cltSuffix};
    		String  kernelPrefix= this.cltKernelPrefix;
			CalibrationFileManagement.MultipleExtensionsFileFilter kernelFilter =
				new CalibrationFileManagement.MultipleExtensionsFileFilter(kernelPrefix,extensions,kernelPrefix+
						"*"+extensions[0]+" CLT symmetrical kernel files");
			if (debugLevel>1) System.out.println("selectKernelFiles("+debugLevel+"): defaultPaths[0]="+defaultPaths[0]+" "+kernelPrefix+"*"+extensions[0]);

			String [] kernelFiles=null;
// try reading all matching files
			File dir= new File (kernelDirectory);
			File [] fileList=null;
			if (dir.exists()) {
				fileList=dir.listFiles(kernelFilter);
			}
			if ((fileList==null) || (fileList.length==0)){
				kernelFiles=CalibrationFileManagement.selectFiles(false,
    					"Select CLT kernel files",
    					"Select",
    					kernelFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    			if ((kernelFiles!=null) && (kernelFiles.length>0)){
    				kernelDirectory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
    				dir= new File (kernelDirectory);
    				fileList=dir.listFiles(kernelFilter);
    				this.cltKernelDirectory= kernelDirectory;
    			}
			}
			if ((fileList==null) || (fileList.length==0)) return null;
			if (debugLevel>1) System.out.println("CLT kernel directory "+kernelDirectory+" has "+fileList.length+" matching files.");
			kernelFiles = new String[fileList.length];
			for (int i=0;i<kernelFiles.length;i++) kernelFiles[i]=fileList[i].getPath();
			String directory=kernelFiles[0].substring(0, kernelFiles[0].lastIndexOf(Prefs.getFileSeparator()));
			String prefix=kernelFiles[0].substring(directory.length()+1, kernelFiles[0].length()-extensions[0].length()-2); // all but NN
			this.cltKernelDirectory=directory;
			this.cltKernelPrefix=prefix;
			return kernelFiles;
    	}
    	
    	public String selectSourceDirectory(boolean smart, boolean newAllowed) { // normally newAllowed=false
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Source (acquired from the camera) image directory", // title
    				"Select source directory", // button
    				null, // filter
    				this.sourceDirectory); // this.sourceDirectory);
    		if (dir!=null) this.sourceDirectory=dir;
    		return dir;
    	}
    	public String selectSensorDirectory(boolean smart,  boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Sensor calibration directory", // title
    				"Select sensor calibration directory", // button
    				null, // filter
    				this.sensorDirectory); //this.sourceDirectory);
    		if (dir!=null) this.sensorDirectory=dir;
    		return dir;
    	}
    	public String selectSharpKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Aberration kernels (sharp) directory", // title
    				"Select aberration kernels (sharp) directory", // button
    				null, // filter
    				this.sharpKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.sharpKernelDirectory=dir;
    		return dir;
    	}
    	
    	public String selectSmoothKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Aberration kernels (smooth) directory", // title
    				"Select aberration kernels (smooth) directory", // button
    				null, // filter
    				this.smoothKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.smoothKernelDirectory=dir;
    		return dir;
    	}
    	
    	public String selectDCTKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"DCT aberration kernels directory (sym and asym files)", // title
    				"Select DCT aberration kernel sdirectory", // button
    				null, // filter
    				this.dctKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.dctKernelDirectory=dir;
    		return dir;
    	}
    	
    	public String selectCLTKernelDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"CLT aberration kernels directory", // title
    				"Select CLT aberration kernels directory", // button
    				null, // filter
    				this.cltKernelDirectory); //this.sourceDirectory);
    		if (dir!=null) this.cltKernelDirectory=dir;
    		return dir;
    	}

    	public String selectX3dDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"x3d output directory", // title
    				"Select x3d output directory", // button
    				null, // filter
    				this.x3dDirectory); //this.sourceDirectory);
    		if (dir!=null) this.x3dDirectory=dir;
    		return dir;
    	}
    	
    	public String selectEquirectangularDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Equirectangular maps directory", // title
    				"Select equirectangular maps directory", // button
    				null, // filter
    				this.equirectangularDirectory);
    		if (dir!=null) this.equirectangularDirectory=dir;
    		return dir;
    	}
    	public String selectResultsDirectory(boolean smart, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Results directory", // title
    				"Select results directory", // button
    				null, // filter
    				this.resultsDirectory); //this.sourceDirectory);
    		if (dir!=null) this.resultsDirectory=dir;
    		return dir;
    	}
     }
    
    
    /* === Parameter classes === */
    public static class ProcessParameters {
  	    public int numEyesisChannels=3;
  	    public int numEyesisSubChannels=3;
  	    public boolean eyesisMode;
  		public boolean [][] frames=new boolean[3][3];
  		public boolean selectFile;
  		public boolean thisFileOnly;
  		public int     subChannelToProcess;
  		public boolean split;
  		public boolean debayer;
  		public boolean showDebayerEnergy;
  		public boolean saveDebayerEnergy;
  		public boolean deconvolve;
  		public boolean combine;
  		public boolean showDenoiseMask;
  		public boolean saveDenoiseMask;
  		public boolean showChromaDenoiseMask;
  		public boolean saveChromaDenoiseMask;
  		public boolean showNoiseGains;
  		public boolean saveNoiseGains;
  		public boolean colorProc;
  		public boolean blueProc;
  		public boolean toRGB;
  		public boolean rotate;
  		public boolean crop;   // crop to the sennor size 
  		public boolean jpeg;   // convert to RGB and save jpeg (if save is true)
  		public boolean save;
  		public boolean save16; // save 16-bit tiff also if the end result is 8 bit 
  		public boolean save32; // save 32-bit tiff also if the end result is 8 or 16 bit
  		public boolean show;
  		public int     JPEG_quality;
  		public double  JPEG_scale;
  		public boolean saveSettings;

  		public ProcessParameters(
  			boolean eyesisMode,	
  			boolean frames_11,
  			boolean frames_12,
  			boolean frames_13,
  			boolean frames_21,
  			boolean frames_22,
  			boolean frames_23,
  			boolean frames_31,
  			boolean frames_32,
  			boolean frames_33,
  			boolean selectFile, // ask for file(s) to process
  			boolean thisFileOnly,
  			int     subChannelToProcess,
  			boolean split,
  			boolean debayer,
  			boolean showDebayerEnergy,
  			boolean saveDebayerEnergy,
  			boolean deconvolve,
  			boolean combine,
  			boolean showDenoiseMask,
  			boolean saveDenoiseMask,
  			boolean showChromaDenoiseMask,
  			boolean saveChromaDenoiseMask,
  			boolean showNoiseGains,
  			boolean saveNoiseGains,
  			boolean colorProc,
  			boolean blueProc,
  			boolean toRGB,
  			boolean rotate,
  			boolean crop,   // crop to the sennor size 
  			boolean jpeg,   // convert to RGB and save jpeg (if save is true)
  			boolean save,
  			boolean save16, // save 16-bit tiff also if the end result is 8 bit 
  			boolean save32, // save 32-bit tiff also if the end result is 8 or 16 bit
  			boolean show,
  			int     JPEG_quality,
  			double  JPEG_scale,
  			boolean saveSettings
  		) {
  			this.eyesisMode=eyesisMode;
  			this.frames[0][0]=frames_11;
  			this.frames[0][1]=frames_12;
  			this.frames[0][2]=frames_13;
  			this.frames[1][0]=frames_21;
  			this.frames[1][1]=frames_22;
  			this.frames[1][2]=frames_23;
  			this.frames[2][0]=frames_31;
  			this.frames[2][1]=frames_32;
  			this.frames[2][2]=frames_33;
  			this.selectFile=selectFile;
  			this.thisFileOnly=thisFileOnly;
  			this.subChannelToProcess=subChannelToProcess;
  			this.split=split;
  			this.debayer=debayer;
  			this.showDebayerEnergy=showDebayerEnergy;
  			this.saveDebayerEnergy=saveDebayerEnergy;
  			this.deconvolve=deconvolve;
  			this.combine=combine;
  			this.showDenoiseMask=showDenoiseMask;
  			this.saveDenoiseMask=saveDenoiseMask;
  			this.showNoiseGains=showNoiseGains;
  			this.saveNoiseGains=saveNoiseGains;
  			this.showChromaDenoiseMask=showChromaDenoiseMask;
  			this.saveChromaDenoiseMask=saveChromaDenoiseMask;
  			this.colorProc=colorProc;
  			this.blueProc=blueProc;
  			this.toRGB=toRGB;
  			this.rotate=rotate;
  			this.crop=crop;
  			this.jpeg=jpeg;
  			this.save=save;
  			this.save16=save16;
  			this.save32=save32;
  			this.show=show;
  			this.JPEG_quality=JPEG_quality;
  			this.JPEG_scale=  JPEG_scale;
  			this.saveSettings=saveSettings;
  		}
  		public void setProperties(String prefix,Properties properties){
  			int i,j;
  			properties.setProperty(prefix+"numEyesisChannels",this.numEyesisChannels+"");
  			properties.setProperty(prefix+"numEyesisSubChannels",this.numEyesisSubChannels+"");
  			properties.setProperty(prefix+"eyesisMode",this.eyesisMode+"");
  		    for (i=0;i<this.frames.length;i++) for (j=0;j<this.frames[0].length;j++)
  				properties.setProperty(prefix+"frames_"+i+"_"+j,this.frames[i][j]+"");
  			properties.setProperty(prefix+"selectFile",this.selectFile+"");
  			properties.setProperty(prefix+"thisFileOnly",this.thisFileOnly+"");
  			properties.setProperty(prefix+"subChannelToProcess",this.subChannelToProcess+"");
  			properties.setProperty(prefix+"split",this.split+"");
  			properties.setProperty(prefix+"debayer",this.debayer+"");
  			properties.setProperty(prefix+"showDebayerEnergy",this.showDebayerEnergy+"");
  			properties.setProperty(prefix+"saveDebayerEnergy",this.saveDebayerEnergy+"");
  			properties.setProperty(prefix+"deconvolve",this.deconvolve+"");
  			properties.setProperty(prefix+"combine",this.combine+"");
  			properties.setProperty(prefix+"showDenoiseMask",this.showDenoiseMask+"");
  			properties.setProperty(prefix+"saveDenoiseMask",this.saveDenoiseMask+"");
  			properties.setProperty(prefix+"showChromaDenoiseMask",this.showChromaDenoiseMask+"");
  			properties.setProperty(prefix+"saveChromaDenoiseMask",this.saveChromaDenoiseMask+"");
  			properties.setProperty(prefix+"showNoiseGains",this.showNoiseGains+"");
  			properties.setProperty(prefix+"saveNoiseGains",this.saveNoiseGains+"");
  			properties.setProperty(prefix+"colorProc",this.colorProc+"");
			properties.setProperty(prefix+"blueProc",this.blueProc+"");
  			properties.setProperty(prefix+"toRGB",this.toRGB+"");
  			properties.setProperty(prefix+"rotate",this.rotate+"");
  			properties.setProperty(prefix+"crop",this.crop+"");
  			properties.setProperty(prefix+"jpeg",this.jpeg+"");
  			properties.setProperty(prefix+"save",this.save+"");
  			properties.setProperty(prefix+"save16",this.save16+"");
  			properties.setProperty(prefix+"save32",this.save32+"");
  			properties.setProperty(prefix+"show",this.show+"");
  			properties.setProperty(prefix+"JPEG_quality",this.JPEG_quality+"");
  			properties.setProperty(prefix+"JPEG_scale",this.JPEG_scale+"");
  			properties.setProperty(prefix+"saveSettings",this.saveSettings+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			int i,j;
  			if (properties.getProperty(prefix+"numEyesisChannels")!=null) this.numEyesisChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisChannels"));
  			if (properties.getProperty(prefix+"numEyesisSubChannels")!=null) this.numEyesisSubChannels=Integer.parseInt(properties.getProperty(prefix+"numEyesisSubChannels"));
  			if (properties.getProperty(prefix+"eyesisMode")!=null) this.eyesisMode=Boolean.parseBoolean(properties.getProperty(prefix+"eyesisMode"));
  		    for (i=0;i<this.frames.length;i++) for (j=0;j<this.frames[0].length;j++)
  		    	if (properties.getProperty(prefix+"frames_"+i+"_"+j)!=null) this.frames[i][j]=Boolean.parseBoolean(properties.getProperty(prefix+"frames_"+i+"_"+j));
  		    if (properties.getProperty(prefix+"selectFile")!=null) this.selectFile=Boolean.parseBoolean(properties.getProperty(prefix+"selectFile"));
  		    if (properties.getProperty(prefix+"thisFileOnly")!=null) this.thisFileOnly=Boolean.parseBoolean(properties.getProperty(prefix+"thisFileOnly"));
  		    if (properties.getProperty(prefix+"subChannelToProcess")!=null) this.subChannelToProcess=Integer.parseInt(properties.getProperty(prefix+"subChannelToProcess"));
  		    if (properties.getProperty(prefix+"split")!=null) this.split=Boolean.parseBoolean(properties.getProperty(prefix+"split"));
  		    if (properties.getProperty(prefix+"debayer")!=null) this.debayer=Boolean.parseBoolean(properties.getProperty(prefix+"debayer"));
  		    if (properties.getProperty(prefix+"showDebayerEnergy")!=null) this.showDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"showDebayerEnergy"));
  		    if (properties.getProperty(prefix+"saveDebayerEnergy")!=null) this.saveDebayerEnergy=Boolean.parseBoolean(properties.getProperty(prefix+"saveDebayerEnergy"));
  		    if (properties.getProperty(prefix+"deconvolve")!=null) this.deconvolve=Boolean.parseBoolean(properties.getProperty(prefix+"deconvolve"));
  		    if (properties.getProperty(prefix+"combine")!=null) this.combine=Boolean.parseBoolean(properties.getProperty(prefix+"combine"));
  		    if (properties.getProperty(prefix+"showDenoiseMask")!=null) this.showDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveDenoiseMask")!=null) this.saveDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveDenoiseMask"));
  		    if (properties.getProperty(prefix+"showChromaDenoiseMask")!=null) this.showChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"showChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"saveChromaDenoiseMask")!=null) this.saveChromaDenoiseMask=Boolean.parseBoolean(properties.getProperty(prefix+"saveChromaDenoiseMask"));
  		    if (properties.getProperty(prefix+"showNoiseGains")!=null) this.showNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"showNoiseGains"));
  		    if (properties.getProperty(prefix+"saveNoiseGains")!=null) this.saveNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"saveNoiseGains"));
  		    if (properties.getProperty(prefix+"colorProc")!=null) this.colorProc=Boolean.parseBoolean(properties.getProperty(prefix+"colorProc"));
  		    if (properties.getProperty(prefix+"blueProc")!=null) this.blueProc=Boolean.parseBoolean(properties.getProperty(prefix+"blueProc"));
  		    if (properties.getProperty(prefix+"toRGB")!=null) this.toRGB=Boolean.parseBoolean(properties.getProperty(prefix+"toRGB"));
  		    if (properties.getProperty(prefix+"rotate")!=null) this.rotate=Boolean.parseBoolean(properties.getProperty(prefix+"rotate"));
  		    if (properties.getProperty(prefix+"crop")!=null) this.crop=Boolean.parseBoolean(properties.getProperty(prefix+"crop"));   // crop to the sensor size 
  		    if (properties.getProperty(prefix+"jpeg")!=null) this.jpeg=Boolean.parseBoolean(properties.getProperty(prefix+"jpeg"));   // convert to RGB and save jpeg (if save is true)
  		    if (properties.getProperty(prefix+"save")!=null) this.save=Boolean.parseBoolean(properties.getProperty(prefix+"save"));
  		    if (properties.getProperty(prefix+"save16")!=null) this.save16=Boolean.parseBoolean(properties.getProperty(prefix+"save16")); // save 16-bit tiff also if the end result is 8 bit 
  		    if (properties.getProperty(prefix+"save32")!=null) this.save32=Boolean.parseBoolean(properties.getProperty(prefix+"save32")); // save 32-bit tiff also if the end result is 8 or 16 bit
  		    if (properties.getProperty(prefix+"show")!=null) this.show=Boolean.parseBoolean(properties.getProperty(prefix+"show"));
  		    if (properties.getProperty(prefix+"JPEG_quality")!=null) this.JPEG_quality=Integer.parseInt(properties.getProperty(prefix+"JPEG_quality"));
  		    if (properties.getProperty(prefix+"JPEG_scale")!=null) this.JPEG_scale=Double.parseDouble(properties.getProperty(prefix+"JPEG_scale"));
  		    if (properties.getProperty(prefix+"saveSettings")!=null) this.saveSettings=Boolean.parseBoolean(properties.getProperty(prefix+"saveSettings"));
  		}
  	}
    
  /* ======================================================================== */
    public static class FilesParameters {
  	  public String [][] rPSFNames=new String [3][3];
  	  public String [][] gaussianNames=new String [3][3];
  	  public String kernelDirectory;
  	  public String resultsDirectory;
  	  public String [] sourceFiles;
  	  public boolean useXML;
  	  public FilesParameters(
  			  String rPSFNames_11,
  			  String rPSFNames_12,
  			  String rPSFNames_13,
  			  String rPSFNames_21,
  			  String rPSFNames_22,
  			  String rPSFNames_23,
  			  String rPSFNames_31,
  			  String rPSFNames_32,
  			  String rPSFNames_33,
  			  String gaussianNames_11,
  			  String gaussianNames_12,
  			  String gaussianNames_13,
  			  String gaussianNames_21,
  			  String gaussianNames_22,
  			  String gaussianNames_23,
  			  String gaussianNames_31,
  			  String gaussianNames_32,
  			  String gaussianNames_33,
  			  String kernelDirectory,
    		      String resultsDirectory,
    		      boolean useXML
    		      ){
  		  this.rPSFNames[0][0]=rPSFNames_11;
  		  this.rPSFNames[0][1]=rPSFNames_12;
  		  this.rPSFNames[0][2]=rPSFNames_13;
  		  this.rPSFNames[1][0]=rPSFNames_21;
  		  this.rPSFNames[1][1]=rPSFNames_22;
  		  this.rPSFNames[1][2]=rPSFNames_23;
  		  this.rPSFNames[2][0]=rPSFNames_31;
  		  this.rPSFNames[2][1]=rPSFNames_32;
  		  this.rPSFNames[2][2]=rPSFNames_33;
  		  this.gaussianNames[0][0]=gaussianNames_11;
  		  this.gaussianNames[0][1]=gaussianNames_12;
  		  this.gaussianNames[0][2]=gaussianNames_13;
  		  this.gaussianNames[1][0]=gaussianNames_21;
  		  this.gaussianNames[1][1]=gaussianNames_22;
  		  this.gaussianNames[1][2]=gaussianNames_23;
  		  this.gaussianNames[2][0]=gaussianNames_31;
  		  this.gaussianNames[2][1]=gaussianNames_32;
  		  this.gaussianNames[2][2]=gaussianNames_33;
  		  this.kernelDirectory=    kernelDirectory;
  		  this.resultsDirectory=   resultsDirectory;
  		  this.useXML=useXML;
  	  }
  		public void setProperties(String prefix,Properties properties){
//  			properties.setProperty(prefix+"",this.+"");
  			int i,j;
  			for (i=0;i<this.rPSFNames.length;i++) for (j=0;j<this.rPSFNames[i].length;j++)
  				properties.setProperty(prefix+"rPSFNames_"+i+"_"+j,this.rPSFNames[i][j]);				
  			for (i=0;i<this.gaussianNames.length;i++) for (j=0;j<this.gaussianNames[i].length;j++)
  				properties.setProperty(prefix+"gaussianNames_"+i+"_"+j,this.gaussianNames[i][j]);
  			properties.setProperty(prefix+"kernelDirectory",this.kernelDirectory);
  			properties.setProperty(prefix+"resultsDirectory",this.resultsDirectory);	
  			properties.setProperty(prefix+"useXML",this.useXML+"");
  			j=(this.sourceFiles==null)?0:this.sourceFiles.length;
  			properties.setProperty(prefix+"sourceFiles_length",j+"");
  			for (i=0;i<j;i++)
  				properties.setProperty(prefix+"sourceFiles_"+i,this.sourceFiles[i]);
  		}
  		public void getProperties(String prefix,Properties properties){
  			int i,j;
  			for (i=0;i<this.rPSFNames.length;i++) for (j=0;j<this.rPSFNames[i].length;j++)
  				this.rPSFNames[i][j]=properties.getProperty(prefix+"rPSFNames_"+i+"_"+j);
  			for (i=0;i<this.gaussianNames.length;i++) for (j=0;j<this.gaussianNames[i].length;j++)
  				this.gaussianNames[i][j]=properties.getProperty(prefix+"gaussianNames_"+i+"_"+j);
  			this.kernelDirectory=properties.getProperty(prefix+"kernelDirectory");
  			this.resultsDirectory=properties.getProperty(prefix+"resultsDirectory");
  			this.useXML=Boolean.parseBoolean(properties.getProperty(prefix+"useXML"));
  			j=Integer.parseInt(properties.getProperty(prefix+"sourceFiles_length"));
  			this.sourceFiles=new String[j];
  			for (i=0;i<j;i++)
  				this.sourceFiles[i]=properties.getProperty(prefix+"sourceFiles_"+i);
  		}
    }

  /* ======================================================================== */
    public static class RGBParameters {
  		public double r_min = 0.075;
  		public double g_min = 0.075;
  		public double b_min = 0.075;
  		public double r_max = 1.0;
  		public double g_max = 1.0;
  		public double b_max = 1.0;
  		public double alpha_min = 0.0;
  		public double alpha_max = 1.0;

/*  		public RGBParameters(double r_min, double g_min, double b_min, double r_max, double g_max, double b_max) {
  			this.r_min = r_min;
  			this.g_min = g_min;
  			this.b_min = b_min;
  			this.r_max = r_max;
  			this.g_max = g_max;
  			this.b_max = b_max;
  		} */
  		public RGBParameters(double r_min, double g_min, double b_min, double r_max, double g_max, double b_max, double alpha_min, double alpha_max) {
  			this.r_min = r_min;
  			this.g_min = g_min;
  			this.b_min = b_min;
  			this.r_max = r_max;
  			this.g_max = g_max;
  			this.b_max = b_max;
  			this.alpha_min = alpha_min;
  			this.alpha_max = alpha_max;
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"r_min",this.r_min+"");
  			properties.setProperty(prefix+"g_min",this.g_min+"");
  			properties.setProperty(prefix+"b_min",this.b_min+"");
  			properties.setProperty(prefix+"r_max",this.r_max+"");
  			properties.setProperty(prefix+"g_max",this.g_max+"");
  			properties.setProperty(prefix+"b_max",this.b_max+"");
  			properties.setProperty(prefix+"alpha_min",this.alpha_min+"");
  			properties.setProperty(prefix+"alpha_max",this.alpha_max+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			this.r_min=Double.parseDouble(properties.getProperty(prefix+"r_min"));
  			this.g_min=Double.parseDouble(properties.getProperty(prefix+"g_min"));
  			this.b_min=Double.parseDouble(properties.getProperty(prefix+"b_min"));
  			this.r_max=Double.parseDouble(properties.getProperty(prefix+"r_max"));
  			this.g_max=Double.parseDouble(properties.getProperty(prefix+"g_max"));
  			this.b_max=Double.parseDouble(properties.getProperty(prefix+"b_max"));
  			if (properties.getProperty(prefix+"alpha_min")!=null)  this.alpha_min=Double.parseDouble(properties.getProperty(prefix+"alpha_min"));
  			if (properties.getProperty(prefix+"alpha_max")!=null)  this.alpha_max=Double.parseDouble(properties.getProperty(prefix+"alpha_max"));
  		}
  		
  	}
  /* ======================================================================== */
    
    public static class ColorProcParameters {
  		public double balanceRed;
  		public double balanceBlue;
  		public double gain;
  		public double weightScaleR;
  		public double weightScaleB;
//  		public double sigma;
  		public double gamma;
  		public double minLin;
  		public double kr;
  		public double kb;
  		public double saturationRed;
  		public double saturationBlue;
  		public boolean useFirstY;
  		public double maskSigma; // LPF for luma to calculate chroma mask
  		public double maskMin; // minimal level for the mask (absolute luma values)
  		public double maskMax; // maximal level for the mask (absolute luma values)
  		public boolean combineWithSharpnessMask; // combine chroma mask with sharpness mask to reduce color leak through borders
  		public double chromaBrightSigma; // LPF for chroma in the bright areas (determined by the mask)
  		public double chromaDarkSigma;   // LPF for chroma in the dark areas (determined by the mask)
	    public boolean corrBlueLeak; //Remove blue color leak in the darks near saturation
  		public int    satDetSquareSize;    // Size of sliding square do detect potential overexposure
  		public double satDetRelDiff;     // Most pixels in square should be within this difference from average
  		public double satDetPartInside;  // Fraction of all pixels in the square to fit inside
  		public double satDetMinFrac;     // minimal value for average compared to average over the whole picture
  		public double satDetFinRelDiff;  // maximal difference from average for the saturated tile to be considered saturated
  		public double satDetGrowRelDiff; // maximal difference from start tile average during growing of overexposed areas 
  		public double satDetNewWeight;   // weight of new pixel when expanding overexposed areas 
  		
  		public int    satDetExpSym; // number of overexposure expand steps, not allowing any brighter
  		public int    satDetExpOver;// number of overexposure expand steps, limited under, any over

  		public int    satDetExpCleanUp; // number of overexposure expand steps, not allowing any brighter (final to clean up oscillations)
  		public double satDetGrowRelDiffCleanUp; // maximal difference from start tile average during growing of overexposed areas (final to clean up oscillations) 
  		
  		public int    blueOverShrink; // shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)
  		public int    blueOverGrow; // grow blue overexposed area by this number of pixels
  		public double blueBandWidth; // average amount of blue leak in pixels
  		public double blueBandWidthDark; // average amount of blue leak in pixels (slope at dark)

  		public double blueNeutral; // Value of Yb/Yrg ratio for the small areas where safe color c an not be found
  		public double blueSolutionRadius;   // How far to trust blue color ratio from the found solution (in pixels)
  		public boolean blueLeakNoHint;      // use blueNeutral in the small areas that do not have reliable color sample
  		public boolean blueLeakNoBrighten; // Do not brighten corrected areas, only darken
  		
  		public boolean blueLeakFixWires; //Fix thin objects with saturated blue, but not R+G
  		public double  blueLeakWiresSize; //size (in pixels) of the small objects to fix blue
  		public double  blueLeakWiresThreshold; // relative to saturation level threshold of the small blue-flooded objects on red+green to process
  		
  		public boolean use8; // use 8 neighbors (false - only 4)
  		
  		public ColorProcParameters(
  				double balanceRed,
  				double balanceBlue,
  				double gain,
  				double weightScaleR,
  				double weightScaleB,
//  				double sigma,
  				double gamma,
  				double minLin,
  				double kr,
  				double kb,
  				double saturationRed,
  				double saturationBlue,
  				boolean useFirstY,
  				double maskSigma, // LPF for luma to calculate chroma mask
  				double maskMin, // minimal level for the mask (absolute luma values)
  				double maskMax, // maximal level for the mask (absolute luma values)
  				boolean combineWithSharpnessMask, // combine chroma mask with sharpness mask to reduce color leak through borders
  				double chromaBrightSigma, // LPF for chroma in the bright areas (determined by the mask)
  				double chromaDarkSigma,   // LPF for chroma in the dark areas (determined by the mask)
  				//----------------------Fixing blue near saturation ----------/
  				// Saturation detect parameters
  			    boolean corrBlueLeak, //Remove blue color leak in the darks near saturation
  				int satDetSquareSize,    // Size of sliding square do detect potential overexposure
  				double satDetRelDiff,    // Most pixels in square should be within this difference from average
  				double satDetPartInside, // Fraction of all pixels in the square to fit inside
  				double satDetMinFrac,    // minimal value for average compared to average over the whole picture
  				double satDetFinRelDiff, // maximal difference from average for the saturated tile to be considered saturated
  				double satDetGrowRelDiff, //maximal difference from start tile average during growing of overexposed areas
  		  		double satDetNewWeight,   // weight of new pixel when expanding overexposed areas 

  				int    satDetExpSym,     // number of overexposure expand steps, not allowing any brighter
  		  		int    satDetExpOver,     // number of overexposure expand steps, limited under, any over
  		  		int    satDetExpCleanUp, // number of overexposure expand steps, not allowing any brighter (final to clean up oscillations)
  		  		double satDetGrowRelDiffCleanUp, // maximal difference from start tile average during growing of overexposed areas (final to clean up oscillations) 
  		  		int    blueOverShrink, // shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)
  		  		int    blueOverGrow, // grow blue overexposed area by this number of pixels
  		  		double blueBandWidth,
  		  		double blueBandWidthDark, // average amount of blue leak in pixels (slope at dark)
  		  		double blueNeutral, // Value of Yb/Yrg ratio for the small areas where safe color c an not be found
  		  		double blueSolutionRadius, // How far to trust blue color ratio from the found solution (in pixels) 
  		  		boolean blueLeakNoHint,    // use blueNeutral in the small areas that do not have reliable color sample
  		  		boolean blueLeakNoBrighten, // Do not brighten corrected areas, only darken
  		  		boolean blueLeakFixWires, //Fix thin objects with saturated blue, but not R+G
  		  		double  blueLeakWiresSize, //size (in pixels) of the small objects to fix blue
  		  		double  blueLeakWiresThreshold, //size (in pixels) of the small objects to fix blue
  		  		boolean use8 // use 8 neighbors (false - only 4)
               ) {
  			this.balanceRed = balanceRed;
  			this.balanceBlue = balanceBlue;
  			this.gain = gain;
  			this.weightScaleR = weightScaleR;
  			this.weightScaleB = weightScaleB;
//  			this.sigma = sigma;
  			this.gamma = gamma;
  			this.minLin = minLin;
  			this.kr = kr;
  			this.kb = kb;
  			this.saturationRed = saturationRed;
  			this.saturationBlue = saturationBlue;
  			this.useFirstY=useFirstY;
  			this.maskSigma=maskSigma;
  			this.maskMin=maskMin;
  			this.maskMax=maskMax;
  			this.combineWithSharpnessMask=combineWithSharpnessMask;
  			this.chromaBrightSigma=chromaBrightSigma;
  			this.chromaDarkSigma=chromaDarkSigma;
  		    this.corrBlueLeak=corrBlueLeak; //Remove blue color leak in the darks near saturation
  			this.satDetSquareSize=satDetSquareSize;    // Size of sliding square do detect potential overexposure
  			this.satDetRelDiff=satDetRelDiff;    // Most pixels in square should be within this difference from average
  			this.satDetPartInside=satDetPartInside; // Fraction of all pixels in the square to fit inside
  			this.satDetMinFrac=satDetMinFrac;    // minimal value for average compared to average over the whole picture
  			this.satDetFinRelDiff=satDetFinRelDiff; // maximal difference from average for the saturated tile to be considered saturated
  			this.satDetGrowRelDiff=satDetGrowRelDiff; //maximal difference from start tile average during growing of overexposed areas
  			this.satDetNewWeight=satDetNewWeight;   // weight of new pixel when expanding overexposed areas 
  			this.satDetExpSym=satDetExpSym;     // number of overexposure expand steps, not allowing any brighter
  			this.satDetExpOver=satDetExpOver;    // number of overexposure expand steps, limited under, any over
  			this.satDetExpCleanUp=satDetExpCleanUp; // number of overexposure expand steps, not allowing any brighter (final to clean up oscillations)
  			this.satDetGrowRelDiffCleanUp=satDetGrowRelDiffCleanUp; // maximal difference from start tile average during growing of overexposed areas (final to clean up oscillations)
  			this.blueOverShrink=blueOverShrink; // shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)
 			this.blueOverGrow=blueOverGrow; // grow blue overexposed area by this number of pixels
  			this.blueBandWidth=blueBandWidth;
  			this.blueBandWidthDark=blueBandWidthDark; // average amount of blue leak in pixels (slope at dark)
  			this.blueNeutral=blueNeutral; // Value of Yb/Yrg ratio for the small areas where safe color c an not be found
  			this.blueSolutionRadius=blueSolutionRadius; // How far to trust blue color ratio from the found solution (in pixels)
  			this.blueLeakNoHint=blueLeakNoHint;    // use blueNeutral in the small areas that do not have reliable color sample
		  	this.blueLeakNoBrighten=blueLeakNoBrighten; // Do not brighten corrected areas, only darken
		  	
		  	this.blueLeakFixWires=blueLeakFixWires; //Fix thin objects with saturated blue, but not R+G
		  	this.blueLeakWiresSize=blueLeakWiresSize; //size (in pixels) of the small objects to fix blue
		  	this.blueLeakWiresThreshold=blueLeakWiresThreshold; //size (in pixels) of the small objects to fix blue
		  	
  			this.use8=use8; // use 8 neighbors (false - only 4)
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"balanceRed",this.balanceRed+"");
  			properties.setProperty(prefix+"balanceBlue",this.balanceBlue+"");
  			properties.setProperty(prefix+"gain",this.gain+"");
  			properties.setProperty(prefix+"weightScaleR",this.weightScaleR+"");
  			properties.setProperty(prefix+"weightScaleB",this.weightScaleB+"");
//  			properties.setProperty(prefix+"sigma",this.sigma+"");
  			properties.setProperty(prefix+"gamma",this.gamma+"");
  			properties.setProperty(prefix+"minLin",this.minLin+"");
  			properties.setProperty(prefix+"kr",this.kr+"");
  			properties.setProperty(prefix+"kb",this.kb+"");
  			properties.setProperty(prefix+"saturationRed",this.saturationRed+"");
  			properties.setProperty(prefix+"saturationBlue",this.saturationBlue+"");
  			properties.setProperty(prefix+"useFirstY",this.useFirstY+"");
  			properties.setProperty(prefix+"maskSigma",this.maskSigma+"");
  			properties.setProperty(prefix+"maskMin",this.maskMin+"");
  			properties.setProperty(prefix+"maskMax",this.maskMax+"");
  			properties.setProperty(prefix+"combineWithSharpnessMask",this.combineWithSharpnessMask+"");
  			properties.setProperty(prefix+"chromaBrightSigma",this.chromaBrightSigma+"");
  			properties.setProperty(prefix+"chromaDarkSigma",this.chromaDarkSigma+"");
  			properties.setProperty(prefix+"corrBlueLeak",this.corrBlueLeak+"");
  			properties.setProperty(prefix+"satDetSquareSize",this.satDetSquareSize+"");
  			properties.setProperty(prefix+"satDetRelDiff",this.satDetRelDiff+"");
  			properties.setProperty(prefix+"satDetPartInside",this.satDetPartInside+"");
  			properties.setProperty(prefix+"satDetMinFrac",this.satDetMinFrac+"");
  			properties.setProperty(prefix+"satDetFinRelDiff",this.satDetFinRelDiff+"");
  			properties.setProperty(prefix+"satDetGrowRelDiff",this.satDetGrowRelDiff+"");
  			
  			
  			properties.setProperty(prefix+"satDetNewWeight",this.satDetNewWeight+"");
  			
  			properties.setProperty(prefix+"satDetExpSym",this.satDetExpSym+"");
  			properties.setProperty(prefix+"satDetExpOver",this.satDetExpOver+"");

  			properties.setProperty(prefix+"satDetExpCleanUp",this.satDetExpCleanUp+"");
  			properties.setProperty(prefix+"satDetGrowRelDiffCleanUp",this.satDetGrowRelDiffCleanUp+"");
  			properties.setProperty(prefix+"blueOverShrink",this.blueOverShrink+"");
  			properties.setProperty(prefix+"blueOverGrow",this.blueOverGrow+"");
  			properties.setProperty(prefix+"blueBandWidth",this.blueBandWidth+"");
  			properties.setProperty(prefix+"blueBandWidthDark",this.blueBandWidthDark+"");
  			
  			properties.setProperty(prefix+"blueNeutral",this.blueNeutral+"");
  			properties.setProperty(prefix+"blueSolutionRadius",this.blueSolutionRadius+"");
 			
  			properties.setProperty(prefix+"blueLeakNoHint",this.blueLeakNoHint+"");
  			properties.setProperty(prefix+"blueLeakNoBrighten",this.blueLeakNoBrighten+"");
  			
  			properties.setProperty(prefix+"blueLeakFixWires",this.blueLeakFixWires+"");
  			properties.setProperty(prefix+"blueLeakWiresSize",this.blueLeakWiresSize+"");
  			properties.setProperty(prefix+"blueLeakWiresThreshold",this.blueLeakWiresThreshold+"");
  			
  			properties.setProperty(prefix+"use8",this.use8+"");
  			
  		}
  		public void getProperties(String prefix,Properties properties){
  			if (properties.getProperty(prefix+"balanceRed")!=null) this.balanceRed=Double.parseDouble(properties.getProperty(prefix+"balanceRed"));
  			if (properties.getProperty(prefix+"balanceBlue")!=null) this.balanceBlue=Double.parseDouble(properties.getProperty(prefix+"balanceBlue"));
  			if (properties.getProperty(prefix+"gain")!=null) this.gain=Double.parseDouble(properties.getProperty(prefix+"gain"));
  			if (properties.getProperty(prefix+"weightScaleR")!=null) this.weightScaleR=Double.parseDouble(properties.getProperty(prefix+"weightScaleR"));
  			if (properties.getProperty(prefix+"weightScaleB")!=null) this.weightScaleB=Double.parseDouble(properties.getProperty(prefix+"weightScaleB"));
//  			this.sigma=Double.parseDouble(properties.getProperty(prefix+"sigma"));
  			if (properties.getProperty(prefix+"gamma")!=null) this.gamma=Double.parseDouble(properties.getProperty(prefix+"gamma"));
  			if (properties.getProperty(prefix+"minLin")!=null) this.minLin=Double.parseDouble(properties.getProperty(prefix+"minLin"));
  			if (properties.getProperty(prefix+"kr")!=null) this.kr=Double.parseDouble(properties.getProperty(prefix+"kr"));
  			if (properties.getProperty(prefix+"kb")!=null) this.kb=Double.parseDouble(properties.getProperty(prefix+"kb"));
  			if (properties.getProperty(prefix+"saturationRed")!=null) this.saturationRed=Double.parseDouble(properties.getProperty(prefix+"saturationRed"));
  			if (properties.getProperty(prefix+"saturationBlue")!=null) this.saturationBlue=Double.parseDouble(properties.getProperty(prefix+"saturationBlue"));
  			if (properties.getProperty(prefix+"useFirstY")!=null) this.useFirstY=Boolean.parseBoolean(properties.getProperty(prefix+"useFirstY"));
  			if (properties.getProperty(prefix+"maskSigma")!=null) this.maskSigma=Double.parseDouble(properties.getProperty(prefix+"maskSigma"));
  			if (properties.getProperty(prefix+"maskMin")!=null) this.maskMin=Double.parseDouble(properties.getProperty(prefix+"maskMin"));
  			if (properties.getProperty(prefix+"maskMax")!=null) this.maskMax=Double.parseDouble(properties.getProperty(prefix+"maskMax"));
  			if (properties.getProperty(prefix+"combineWithSharpnessMask")!=null) this.combineWithSharpnessMask=Boolean.parseBoolean(properties.getProperty(prefix+"combineWithSharpnessMask"));
  			if (properties.getProperty(prefix+"chromaBrightSigma")!=null) this.chromaBrightSigma=Double.parseDouble(properties.getProperty(prefix+"chromaBrightSigma"));
  			if (properties.getProperty(prefix+"chromaDarkSigma")!=null) this.chromaDarkSigma=Double.parseDouble(properties.getProperty(prefix+"chromaDarkSigma"));

  			if (properties.getProperty(prefix+"corrBlueLeak")!=null) this.corrBlueLeak=Boolean.parseBoolean(properties.getProperty(prefix+"corrBlueLeak"));
  			if (properties.getProperty(prefix+"satDetSquareSize")!=null) this.satDetSquareSize=Integer.parseInt(properties.getProperty(prefix+"satDetSquareSize"));
  			if (properties.getProperty(prefix+"satDetRelDiff")!=null)    this.satDetRelDiff=Double.parseDouble(properties.getProperty(prefix+"satDetRelDiff"));
  			if (properties.getProperty(prefix+"satDetPartInside")!=null) this.satDetPartInside=Double.parseDouble(properties.getProperty(prefix+"satDetPartInside"));
  			if (properties.getProperty(prefix+"satDetMinFrac")!=null)    this.satDetMinFrac=Double.parseDouble(properties.getProperty(prefix+"satDetMinFrac"));
  			if (properties.getProperty(prefix+"satDetFinRelDiff")!=null) this.satDetFinRelDiff=Double.parseDouble(properties.getProperty(prefix+"satDetFinRelDiff"));
  			if (properties.getProperty(prefix+"satDetGrowRelDiff")!=null) this.satDetGrowRelDiff=Double.parseDouble(properties.getProperty(prefix+"satDetGrowRelDiff"));
  			if (properties.getProperty(prefix+"satDetNewWeight")!=null) this.satDetNewWeight=Double.parseDouble(properties.getProperty(prefix+"satDetNewWeight"));
  			if (properties.getProperty(prefix+"satDetExpSym")!=null) this.satDetExpSym=Integer.parseInt(properties.getProperty(prefix+"satDetExpSym"));
  			if (properties.getProperty(prefix+"satDetExpOver")!=null) this.satDetExpOver=Integer.parseInt(properties.getProperty(prefix+"satDetExpOver"));
  			if (properties.getProperty(prefix+"satDetExpCleanUp")!=null) this.satDetExpCleanUp=Integer.parseInt(properties.getProperty(prefix+"satDetExpCleanUp"));
  			if (properties.getProperty(prefix+"satDetGrowRelDiffCleanUp")!=null) this.satDetGrowRelDiffCleanUp=Double.parseDouble(properties.getProperty(prefix+"satDetGrowRelDiffCleanUp"));
  			if (properties.getProperty(prefix+"blueOverShrink")!=null) this.blueOverShrink=Integer.parseInt(properties.getProperty(prefix+"blueOverShrink"));
  			
//  			if (properties.getProperty(prefix+"blueOverGrow")!=null) this.blueOverGrow=Integer.parseInt(properties.getProperty(prefix+"blueOverGrow"));
  			if (properties.getProperty(prefix+"blueOverGrow")!=null) this.blueOverGrow=(int) Double.parseDouble(properties.getProperty(prefix+"blueOverGrow"));
  			
  			if (properties.getProperty(prefix+"blueBandWidth")!=null) this.blueBandWidth=Double.parseDouble(properties.getProperty(prefix+"blueBandWidth"));
  			if (properties.getProperty(prefix+"blueBandWidthDark")!=null) this.blueBandWidthDark=Double.parseDouble(properties.getProperty(prefix+"blueBandWidthDark"));
  			if (properties.getProperty(prefix+"blueNeutral")!=null) this.blueNeutral=Double.parseDouble(properties.getProperty(prefix+"blueNeutral"));
  			if (properties.getProperty(prefix+"blueSolutionRadius")!=null) this.blueSolutionRadius=Double.parseDouble(properties.getProperty(prefix+"blueSolutionRadius"));
  			if (properties.getProperty(prefix+"blueLeakNoHint")!=null) this.blueLeakNoHint=Boolean.parseBoolean(properties.getProperty(prefix+"blueLeakNoHint"));

  			if (properties.getProperty(prefix+"blueLeakNoBrighten")!=null) this.blueLeakNoBrighten=Boolean.parseBoolean(properties.getProperty(prefix+"blueLeakNoBrighten"));
  			
  			if (properties.getProperty(prefix+"blueLeakFixWires")!=null) this.blueLeakFixWires=Boolean.parseBoolean(properties.getProperty(prefix+"blueLeakFixWires"));
  			if (properties.getProperty(prefix+"blueLeakWiresSize")!=null) this.blueLeakWiresSize=Double.parseDouble(properties.getProperty(prefix+"blueLeakWiresSize"));
  			if (properties.getProperty(prefix+"blueLeakWiresThreshold")!=null) this.blueLeakWiresThreshold=Double.parseDouble(properties.getProperty(prefix+"blueLeakWiresThreshold"));
  			
  			if (properties.getProperty(prefix+"use8")!=null) this.use8=Boolean.parseBoolean(properties.getProperty(prefix+"use8"));
  		}
  	}
    /* ======================================================================== */
  // individual per-channel color balance and gain
    public static class ColorCalibParameters {
  		public double[][] gain=new double[3][3]; 
  		public double[][] balanceRed=new double[3][3];
  		public double[][] balanceBlue=new double[3][3];
  		
  	  	public ColorCalibParameters(
  	  			double gain_11,
  	  			double gain_12,
  	  			double gain_13,
  	  			double gain_21,
  	  			double gain_22,
  	  			double gain_23,
  	  			double gain_31,
  	  			double gain_32,
  	  			double gain_33,
  	  			double balanceRed_11,
  	  			double balanceRed_12,
  	  			double balanceRed_13,
  	  			double balanceRed_21,
  	  			double balanceRed_22,
  	  			double balanceRed_23,
  	  			double balanceRed_31,
  	  			double balanceRed_32,
  	  			double balanceRed_33,
  	  			double balanceBlue_11,
  	  			double balanceBlue_12,
  	  			double balanceBlue_13,
  	  			double balanceBlue_21,
  	  			double balanceBlue_22,
  	  			double balanceBlue_23,
  	  			double balanceBlue_31,
  	  			double balanceBlue_32,
  	  			double balanceBlue_33){
  	  		this.gain[0][0]=gain_11;
  	  		this.gain[0][1]=gain_12;
  	  		this.gain[0][2]=gain_13;
  	  		this.gain[1][0]=gain_21;
  	  		this.gain[1][1]=gain_22;
  	  		this.gain[1][2]=gain_23;
  	  		this.gain[2][0]=gain_31;
  	  		this.gain[2][1]=gain_32;
  	  		this.gain[2][2]=gain_33;
  	  		this.balanceRed[0][0]=balanceRed_11;
  	  		this.balanceRed[0][1]=balanceRed_12;
  	  		this.balanceRed[0][2]=balanceRed_13;
  	  		this.balanceRed[1][0]=balanceRed_21;
  	  		this.balanceRed[1][1]=balanceRed_22;
  	  		this.balanceRed[1][2]=balanceRed_23;
  	  		this.balanceRed[2][0]=balanceRed_31;
  	  		this.balanceRed[2][1]=balanceRed_32;
  	  		this.balanceRed[2][2]=balanceRed_33;
  	  		this.balanceBlue[0][0]=balanceBlue_11;
  	  		this.balanceBlue[0][1]=balanceBlue_12;
  	  		this.balanceBlue[0][2]=balanceBlue_13;
  	  		this.balanceBlue[1][0]=balanceBlue_21;
  	  		this.balanceBlue[1][1]=balanceBlue_22;
  	  		this.balanceBlue[1][2]=balanceBlue_23;
  	  		this.balanceBlue[2][0]=balanceBlue_31;
  	  		this.balanceBlue[2][1]=balanceBlue_32;
  	  		this.balanceBlue[2][2]=balanceBlue_33;
  	    }
  		public void setProperties(String prefix,Properties properties){
  			int i,j;
  			for (i=0;i<this.gain.length;i++) for (j=0;j<this.gain[i].length;j++)
  			  properties.setProperty(prefix+"gain_"+i+"_"+j,this.gain[i][j]+"");
  			for (i=0;i<this.balanceRed.length;i++) for (j=0;j<this.balanceRed[i].length;j++)
  				  properties.setProperty(prefix+"balanceRed_"+i+"_"+j,this.balanceRed[i][j]+"");
  			for (i=0;i<this.balanceBlue.length;i++) for (j=0;j<this.balanceBlue[i].length;j++)
  				  properties.setProperty(prefix+"balanceBlue_"+i+"_"+j,this.balanceBlue[i][j]+"");
  		}
  		public void getProperties(String prefix,Properties properties){
//  			this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
  			int i,j;
  			String s;
  			for (i=0;i<this.gain.length;i++) for (j=0;j<this.gain[i].length;j++) {
  				s=properties.getProperty(prefix+"gain_"+i+"_"+j);
  				if (s!=null) this.gain[i][j]=Double.parseDouble(s);
  			}
  			for (i=0;i<this.balanceRed.length;i++) for (j=0;j<this.balanceRed[i].length;j++) {
  				s=properties.getProperty(prefix+"balanceRed_"+i+"_"+j);
  				if (s!=null) this.balanceRed[i][j]=Double.parseDouble(s);
  			}
  			for (i=0;i<this.balanceBlue.length;i++) for (j=0;j<this.balanceBlue[i].length;j++) {
  				s=properties.getProperty(prefix+"balanceBlue_"+i+"_"+j);
  				if (s!=null) this.balanceBlue[i][j]=Double.parseDouble(s);
  			}
  		}

    }
    /* ======================================================================== */
    public static class NonlinParameters {
    	public boolean useRejectBlocksFilter;
    	public boolean combineBothModes; 
   	public int    maskFFTSize; // 256
   	public int    blockPeriod; // 32
    	public double rejectFreqSigma; // 1.0, frequency domain
    	public double lowPassSigma;    // 5.0, spatial domain
    	public double filtMin;
    	public double filtMax;
  	public double [][] thresholdCorrection=new double[3][3]; // apply to filtMin and filtMax
  	public double [] thresholdCorr={
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
  			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
  	};
    	public double threshold;
  	public boolean useDiffNoiseGains;
  	public double [] noiseGainWeights=new double[3]; 
  	double blurSigma;     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation 
  	public double  noiseGainPower;
    	public boolean showMask;
  // ring filter
    	public boolean useRingFilter;    // filter out spots on denoise mask
      public double minMaxValue;       // minimal value (relative to filtMax) of the local maximum to be processed
      public double overRingThreshold; // ratio of local max. and maximal value in the surrounding ring to trigger filter 
      public double overRingLimit;     // limit values in the center circle to scaled maximum in a ring
      public double ringIR;            // ring inner radius (center circle radius)
      public double ringOR;            // ring outer radius
      
    	public NonlinParameters(
    			boolean useRejectBlocksFilter,
    			boolean combineBothModes,
    			int maskFFTSize,
    			int blockPeriod,
    			double rejectFreqSigma,
    			double lowPassSigma,
    			double filtMin,
    			double filtMax,
    			double thresholdCorrection_11,
    			double thresholdCorrection_12,
    			double thresholdCorrection_13,
    			double thresholdCorrection_21,
    			double thresholdCorrection_22,
    			double thresholdCorrection_23,
    			double thresholdCorrection_31,
    			double thresholdCorrection_32,
    			double thresholdCorrection_33,
    			double threshold,
  			boolean useDiffNoiseGains,
  			double noiseGainWeights_0, // r
  			double noiseGainWeights_1, // b
  			double noiseGainWeights_2, // g
  			double blurSigma,     // blur sigma for mask calculation (blur convolution kernels for noise gain calculation 
  			double noiseGainPower,
  			boolean useRingFilter,    // filter out spots on denoise mask
    		    double minMaxValue,       // minimal value (relative to filtMax) of the local maximum to be processed
    		    double overRingThreshold, // ratio of local max. and maximal value in the surrounding ring to trigger filter 
    		    double overRingLimit,     // limit values in the center circle to scaled maximum in a ring
    		    double ringIR,            // ring inner radius (center circle radius)
    		    double ringOR             // ring outer radius
  			
  ) {
    		this.useRejectBlocksFilter=useRejectBlocksFilter;
    		this.combineBothModes=combineBothModes;
    		this.maskFFTSize = maskFFTSize;
    		this.blockPeriod = blockPeriod;
    		this.rejectFreqSigma = rejectFreqSigma;
    		this.lowPassSigma = lowPassSigma;
    		this.filtMin = filtMin;
    		this.filtMax = filtMax;
    		this.thresholdCorrection[0][0]=thresholdCorrection_11;
    		this.thresholdCorrection[0][1]=thresholdCorrection_12;
    		this.thresholdCorrection[0][2]=thresholdCorrection_13;
    		this.thresholdCorrection[1][0]=thresholdCorrection_21;
    		this.thresholdCorrection[1][1]=thresholdCorrection_22;
    		this.thresholdCorrection[1][2]=thresholdCorrection_23;
    		this.thresholdCorrection[2][0]=thresholdCorrection_31;
    		this.thresholdCorrection[2][1]=thresholdCorrection_32;
    		this.thresholdCorrection[2][2]=thresholdCorrection_33;
    		this.threshold = threshold;
  		this.useDiffNoiseGains=useDiffNoiseGains;
  		this.noiseGainWeights[0]=noiseGainWeights_0;
  		this.noiseGainWeights[1]=noiseGainWeights_1;
  		this.noiseGainWeights[2]=noiseGainWeights_2;
  		this.blurSigma=blurSigma;
  		this.noiseGainPower=noiseGainPower;
  		this.useRingFilter=useRingFilter;
  		this.minMaxValue=minMaxValue;
  		this.overRingThreshold=overRingThreshold; 
  		this.overRingLimit=overRingLimit;
  		this.ringIR=ringIR;
  		this.ringOR=ringOR;
  		
    	}
    	public void modifyNumChannels(int numChannels){
    		if ((numChannels>0) && (numChannels!=this.thresholdCorr.length)){ 
    			double [] thresholdCorr1=this.thresholdCorr;
    			this.thresholdCorr=  new double[numChannels];
    			for (int i=0;i<numChannels;i++) {
    				int j=i;
    				if (j>=thresholdCorr1.length) j=thresholdCorr1.length-1;
    				this.thresholdCorr[i]=thresholdCorr1[j];
    			}
    		}
    	}

  	public void setProperties(String prefix,Properties properties){
//  		properties.setProperty(prefix+"oversample",this.oversample+"");
  		properties.setProperty(prefix+"useRejectBlocksFilter",this.useRejectBlocksFilter+"");
  		properties.setProperty(prefix+"combineBothModes",this.combineBothModes+"");
  		properties.setProperty(prefix+"maskFFTSize",this.maskFFTSize+"");
  		properties.setProperty(prefix+"blockPeriod",this.blockPeriod+"");
  		properties.setProperty(prefix+"rejectFreqSigma",this.rejectFreqSigma+"");
  		properties.setProperty(prefix+"lowPassSigma",this.lowPassSigma+"");
  		properties.setProperty(prefix+"filtMin",this.filtMin+"");
  		properties.setProperty(prefix+"filtMax",this.filtMax+"");
  		for (int i=0;i<this.thresholdCorrection.length;i++) for (int j=0;j<this.thresholdCorrection[i].length;j++)
  		  properties.setProperty(prefix+"thresholdCorrection_"+i+"_"+j,this.thresholdCorrection[i][j]+"");
  		
  		
  		properties.setProperty(prefix+"threshold",this.threshold+"");
  		properties.setProperty(prefix+"useDiffNoiseGains",this.useDiffNoiseGains+"");
  		for (int i=0;i<this.noiseGainWeights.length;i++)
  		   properties.setProperty(prefix+"noiseGainWeights_"+i,this.noiseGainWeights[i]+"");
  		properties.setProperty(prefix+"blurSigma",this.blurSigma+"");
  		properties.setProperty(prefix+"noiseGainPower",this.noiseGainPower+"");
  		properties.setProperty(prefix+"useRingFilter",this.useRingFilter+"");
  		properties.setProperty(prefix+"minMaxValue",this.minMaxValue+"");
  		properties.setProperty(prefix+"overRingThreshold",this.overRingThreshold+""); 
  		properties.setProperty(prefix+"overRingLimit",this.overRingLimit+"");
  		properties.setProperty(prefix+"ringIR",this.ringIR+"");
  		properties.setProperty(prefix+"ringOR",this.ringOR+"");
  		properties.setProperty(prefix+"thresholdCorr",this.thresholdCorr.length+"");
  		for (int i =0;i<this.thresholdCorr.length;i++) properties.setProperty(prefix+"thresholdCorr_"+i,this.thresholdCorr[i]+"");
  	}
  	public void getProperties(String prefix,Properties properties){
  		//  		this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
  		String s;
  		this.useRejectBlocksFilter=Boolean.parseBoolean(properties.getProperty(prefix+"useRejectBlocksFilter"));
  		this.combineBothModes=Boolean.parseBoolean(properties.getProperty(prefix+"combineBothModes")); 
  		this.maskFFTSize=Integer.parseInt(properties.getProperty(prefix+"maskFFTSize"));
  		this.blockPeriod=Integer.parseInt(properties.getProperty(prefix+"blockPeriod"));
  		this.rejectFreqSigma=Double.parseDouble(properties.getProperty(prefix+"rejectFreqSigma"));
  		this.lowPassSigma=Double.parseDouble(properties.getProperty(prefix+"lowPassSigma"));
  		this.filtMin=Double.parseDouble(properties.getProperty(prefix+"filtMin"));
  		this.filtMax=Double.parseDouble(properties.getProperty(prefix+"filtMax"));
  		for (int i=0;i<this.thresholdCorrection.length;i++) for (int j=0;j<this.thresholdCorrection[i].length;j++)
  			this.thresholdCorrection[i][j]=Double.parseDouble(properties.getProperty(prefix+"thresholdCorrection_"+i+"_"+j));
  		this.threshold=Double.parseDouble(properties.getProperty(prefix+"threshold"));
  		this.useDiffNoiseGains=Boolean.parseBoolean(properties.getProperty(prefix+"useDiffNoiseGains"));
  		for (int i=0;i<this.noiseGainWeights.length;i++)
  			this.noiseGainWeights[i]=Double.parseDouble(properties.getProperty(prefix+"noiseGainWeights_"+i)); 
  		this.blurSigma=Double.parseDouble(properties.getProperty(prefix+"blurSigma"));
  		this.noiseGainPower=Double.parseDouble(properties.getProperty(prefix+"noiseGainPower"));
  		s=properties.getProperty(prefix+"useRingFilter");
  		if ((s==null) || (s=="")) return; // earlier revision
  		this.useRingFilter=Boolean.parseBoolean(properties.getProperty(prefix+"useRingFilter"));
  		this.minMaxValue=Double.parseDouble(properties.getProperty(prefix+"minMaxValue"));
  		this.overRingThreshold=Double.parseDouble(properties.getProperty(prefix+"overRingThreshold")); 
  		this.overRingLimit=Double.parseDouble(properties.getProperty(prefix+"overRingLimit"));
  		this.ringIR=Double.parseDouble(properties.getProperty(prefix+"ringIR"));
  		this.ringOR=Double.parseDouble(properties.getProperty(prefix+"ringOR"));
  		if (properties.getProperty(prefix+"thresholdCorr")!=null){
  			this.thresholdCorr=new double [Integer.parseInt((String) properties.getProperty(prefix+"thresholdCorr"))];
  			for (int i=0;i<this.thresholdCorr.length;i++) this.thresholdCorr[i]=Double.parseDouble((String)properties.getProperty(prefix+"thresholdCorr_"+i));
  		}
  	}

    }
  /* ======================================================================== */
    public static class SplitParameters {
  		public int oversample;
  		public int addLeft;
  		public int addTop;
  		public int addRight;
  		public int addBottom;

  		public SplitParameters(int oversample, int addLeft, int addTop,
  				int addRight, int addBottom) {
  			this.oversample = oversample;
  			this.addLeft = addLeft;
  			this.addTop = addTop;
  			this.addRight = addRight;
  			this.addBottom = addBottom;
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"oversample",this.oversample+"");
  			properties.setProperty(prefix+"addLeft",   this.addLeft+"");
  			properties.setProperty(prefix+"addTop",    this.addTop+"");
  			properties.setProperty(prefix+"addRight",  this.addRight+"");
  			properties.setProperty(prefix+"addBottom", this.addBottom+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			this.oversample=Integer.parseInt(properties.getProperty(prefix+"oversample"));
  			this.addLeft=Integer.parseInt(properties.getProperty(prefix+"addLeft"));
  			this.addTop=Integer.parseInt(properties.getProperty(prefix+"addTop"));
  			this.addRight=Integer.parseInt(properties.getProperty(prefix+"addRight"));
  			this.addBottom=Integer.parseInt(properties.getProperty(prefix+"addBottom"));
  		}
  	}
    public static class CLTParameters {
  		public int        transform_size =      8; //
  		public int        clt_window =          1; // currently only 3 types of windows - 0 (none), 1 and 2
  		public double     shift_x =           0.0;
  		public double     shift_y =           0.0;
  		public int        iclt_mask =          15;  // which transforms to combine
  		public int        tileX =             258;  // number of kernel tile (0..163) 
  		public int        tileY =             133;  // number of kernel tile (0..122)
  		public int        dbg_mode =            0;  // 0 - normal, +1 - no DCT/IDCT
  		public int        ishift_x =            0;  // debug feature - shift source image by this pixels left
  		public int        ishift_y =            0;  // debug feature - shift source image by this pixels down
  		public double     fat_zero =          0.0;  // modify phase correlation to prevent division by very small numbers
  		public double     corr_sigma =        0.8;  // LPF correlation sigma
  		public boolean    norm_kern =         true; // normalize kernels
  		public boolean    gain_equalize =     false;// equalize green channel gain
  		public boolean    colors_equalize =   true; // equalize R/G, B/G of the individual channels
  		public double     novignetting_r    = 0.2644; // reg gain in the center of sensor calibration R (instead of vignetting)
  		public double     novignetting_g    = 0.3733; // green gain in the center of sensor calibration G
  		public double     novignetting_b    = 0.2034; // blue gain in the center of sensor calibration B
  		public double     scale_r =           1.0; // extra gain correction after vignetting or nonvignetting, before other processing
  		public double     scale_g =           1.0;
  		public double     scale_b =           1.0;
  		public double     vignetting_max    = 0.4; // value in vignetting data to correspond to 1x in the kernel
  		public double     vignetting_range  = 5.0; // do not try to correct vignetting less than vignetting_max/vignetting_range
  		public int        kernel_step =       16;  // source kernels step in pixels (have 1 kernel margin on each side)  
  		public double     disparity  =        0.0; // nominal disparity between side of square cameras (pix)
  		public double     z_correction  =     0.0; // Inverse distance to infinity (misalignment cortrection)
  		public boolean    correlate =         true; // calculate correlation
  		public int        corr_mask =         15;  // bitmask of pairs to combine in the composite
  		public boolean    corr_sym =          false; // combine correlation with mirrored around disparity direction
  		public boolean    corr_keep =         true;  // keep all partial correlations (otherwise - only combined one)
  		public boolean    corr_show =         false; // Show combined correlations
  		public boolean    corr_mismatch=      false; // calculate per-pair X/Y variations of measured correlations 
// TODO: what to do if some occlusion is present (only some channels correlate)  		
  		public double     corr_offset =       0.1; //0.1;  // add to pair correlation before multiplying by other pairs (between sum and product)
  		                                            // negative - add, not mpy
  		public double     corr_red =          0.5;  // Red to green correlation weight 
  		public double     corr_blue =         0.2;  // Blue to green correlation weight
  		public boolean    corr_normalize =    false; // normalize each correlation tile by rms
  		public double     min_corr =          0.001; // minimal correlation value to consider valid 
  		public double     min_corr_normalized =  2.0; // minimal correlation value to consider valid when normalizing correlation results 
  		public double     max_corr_sigma =    1.5;  // weights of points around global max to find fractional
  		                                            // pixel location by quadratic approximation
  		public double     max_corr_radius =   3.5;  // maximal distance from int max to consider
  		
  		public int        enhortho_width =    2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
  		public double     enhortho_scale =    0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)
  		
  		public boolean    max_corr_double =   false; // NOT USED double pass when masking center of mass to reduce preference for integer values
  		public int        corr_mode =         2;    // which correlation mode to use: 0 - integer max, 1 - center of mass, 2 - polynomial
  		
          // pixel location by quadratic approximation
  		public double     corr_border_contrast = 0.01; // contrast of dotted border on correlation results
  		
  		public int        tile_task_op =      0xff;   // bitmask of operation modes applied to tiles (0 - nothing), bits TBD later
  		                                           // +(0..f) - images, +(00.f0) - process pairs + 256 - force disparity when combining images
  		// window to process tiles (later arbitrary masks will be generated to follow particular stages);
  		public int        tile_task_wl =      0;   // 
  		public int        tile_task_wt =      0;   // 
  		public int        tile_task_ww =      324; // 
  		public int        tile_task_wh =      242; //
  		public double     min_shot =          10.0;  // Do not adjust for shot noise if lower than
  		public double     scale_shot =        3.0;   // scale when dividing by sqrt
  		
  		public double     diff_sigma =        5.0;   // RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
  		public double     diff_threshold =    1.5;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
  		public boolean    diff_gauss =        true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
  		public double     min_agree =         3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
  		public boolean    dust_remove =       true;  // Do not reduce average weight when only one image differes much from the average
  		
  		public boolean    black_back =        true;  // use Black for backdrop outside of the FOV
  		public boolean    keep_weights =      true;  // add port weights to RGBA stack (debug feature)
  		public boolean    sharp_alpha =       false; // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
  		public double     alpha0 = 	          0.6; // > .525 Alpha channel 0.0 thereshold (lower - transparent) (watch for far objects)
  		public double     alpha1 = 	          0.8; // Alpha channel 1.0 threshold (higher - opaque) (watch for window dust)
  		
  		public boolean    gen_chn_stacks =    false; // generate shifted channel rgb stacks
  		public boolean    gen_chn_img =       true;  // generate shifted channel images
  		public boolean    gen_4_img =         true;  // Generate shifted channel images and save with the model
  		public boolean    show_nonoverlap =   true;  // show result RGBA before overlap combined (first channels, then RGBA combined?)
  		public boolean    show_overlap =      true;  // show result RGBA (first channels, then RGBA combined?)
  		public boolean    show_rgba_color =   true;  // show combined color image
  		public boolean    show_map =          true;  // show disparity maps
  		
  		public double     disp_scan_start =   0.0;   // disparity scan start value
  		public double     disp_scan_step =    1.0;   // disparity scan step
  		public int        disp_scan_count =   10;    // disparity scan number of measurements
  		
  		public boolean    fine_dbg =          false; // Debug infinity/lazy eye correction
  		public double     fine_corr_x_0 =     0.0;   // additionally shift image in port 0 in x direction
  		public double     fine_corr_y_0 =     0.0;   // additionally shift image in port 0 in y direction
  		public double     fine_corr_x_1 =     0.0;   // additionally shift image in port 1 in x direction
  		public double     fine_corr_y_1 =     0.0;   // additionally shift image in port 1 in y direction
  		public double     fine_corr_x_2 =     0.0;   // additionally shift image in port 2 in x direction
  		public double     fine_corr_y_2 =     0.0;   // additionally shift image in port 2 in y direction
  		public double     fine_corr_x_3 =     0.0;   // additionally shift image in port 3 in x direction
  		public double     fine_corr_y_3 =     0.0;   // additionally shift image in port 3 in y direction
  		
  		public double     fcorr_radius =       0.75 ; // Do not try to correct outside this fraction of width/hight
  		public double     fcorr_min_strength = 0.15 ; // 0.005 minimal correlation strength to apply fine correction
  		public double     fcorr_disp_diff =   1.5;   // consider only tiles with absolute residual disparity lower than
  		public boolean    fcorr_quadratic =   true;  // Use quadratic polynomial for fine correction (false - only linear)
  		public boolean    fcorr_ignore =      false; // Ignore currently calculated fine correction
  		public double     fcorr_inf_strength = 0.20 ; // Minimal correlation strength to use for infinity correction
  		public double     fcorr_inf_diff =    0.2;   // Disparity half-range for infinity
  		public boolean    fcorr_inf_quad =    true;  // Use quadratic polynomial for infinity correction (false - only linear)
  		public boolean    fcorr_inf_vert =    false; // Correct infinity in vertical direction (false - only horizontal)

//--
  		public boolean    inf_disp_apply =   true;   // Apply disparity correction to zero at infinity 
  		public int        inf_repeat =       5;      // Re run disparity correction at infinity multiple times 
//  		public boolean    inf_mism_apply =   true;   // Apply lazy eye correction at infinity
  		
  		public int        inf_iters =        20;     // Infinity extraction - maximum iterations
  		public double     inf_final_diff =   0.0001; // Coefficients maximal increment to exit iterations
  		public double     inf_far_pull =     0.0;    // include farther tiles than tolerance, but scale their weights
  		
  		// infinity filter
  		public double     inf_str_pow =      1.0;    // Strength power for infinity filtering
  		public int        inf_smpl_side =    3;      // Sample size (side of a square) for infinity filtering
  		public int        inf_smpl_num =     5;      // Number after removing worst (should be >1) for infinity filtering
  		public double     inf_smpl_rms =     0.1;    // Maximal RMS of the remaining tiles in a sample for infinity filtering

  		//Histogram infinity filter
  		public int        ih_smpl_step =     8;      // Square sample step (50% overlap) 
  		public double     ih_disp_min =     -1.0;    // Minimal disparity
  		public double     ih_disp_step =     0.05;   // Disparity step
  		public int        ih_num_bins =     40;      // Number of bins
  		public double     ih_sigma =         0.1;    // Gaussian sigma (in disparity pixels)
  		public double     ih_max_diff =      0.1;    // Keep samples within this difference from farthest maximum
  		public int        ih_min_samples =  10;      // Minimal number of remaining samples
  		public boolean    ih_norm_center =  true;    // Replace samples with a single average with equal weight

  		// Lazy eye parameters
  		public int        ly_smpl_side =    3;       // Sample size (side of a square)
  		public int        ly_smpl_num =     5;       // Number after removing worst (should be >1)
 		public double     ly_meas_disp =    1.5;     // Maximal measured relative disparity
 		public double     ly_smpl_rms =     0.2;     // 1;     // Maximal RMS of the remaining tiles in a sample
		public double     ly_disp_var =     0.5;     // 2;     // Maximal full disparity difference to 8 neighbors
 		public double     ly_inf_frac =     0.5;     // Relative weight of infinity calibration data
  		public boolean    ly_on_scan =      true;    // Calculate and apply lazy eye correction after disparity scan (poly or extrinsic) 
  		public boolean    ly_inf_en =       true;    // Simultaneously correct disparity at infinity (both poly and extrinsic) 
  		public boolean    ly_inf_force=     false;   // Force convergence correction during extrinsic, even with no infinity data 
  		public boolean    ly_poly =         false;   // Use polynomial correction, false - correct tilt/azimuth/roll of each sensor 
  		
  		

  		// old fcorr parameters, reuse?
// 		public int        fcorr_sample_size = 32;    // Use square this size side to detect outliers
// 		public int        fcorr_mintiles =    8;     // Keep tiles only if there are more in each square 
// 		public double     fcorr_reloutliers = 0.5;   // Remove this fraction of tiles from each sample
// 		public double     fcorr_sigma =       20.0;  // Gaussian blur channel mismatch data
  		
  		public double     corr_magic_scale =  0.85;  // reported correlation offset vs. actual one (not yet understood)
  		
  		// 3d reconstruction
  		public boolean    show_textures    = true;  // show generated textures
  		public boolean    debug_filters    = false;// show intermediate results of filtering
  		// not used anywhere so far
  		public double     min_smth         = 0.25;  // 0.25 minimal noise-normalized pixel difference in a channel to suspect something    
  		public double     sure_smth        = 2.0;   // reliable noise-normalized pixel difference in a channel to have something
  		public double     bgnd_range       = 0.3;   // disparity range to be considered background
  		public double     other_range      = 2.0;   // disparity difference from center (provided) disparity to trust
  		
  		public double     ex_strength      = 0.18;  // minimal 4-corr strength to trust tile
  		public double     ex_nstrength     = 0.4;   // minimal 4-corr strength divided by channel diff for new (border) tiles
  		
  		public boolean    ex_over_bgnd     = false; // Allow expansion over previously identified background (infinity)
  		public double     ex_min_over      = 1.0;   // When expanding over background, disregard lower disparity 
  		
  		public double     pt_super_trust   = 1.6;   // If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds
  		public boolean    pt_keep_raw_fg   = true;  // Do not replace raw tiles by the plates, if raw is closer (like poles)
  		public double     pt_scale_pre     = 1.5;   // Scale plates strength before comparing to raw strength
  		public double     pt_scale_post    = 2.5;   // Scale plates strength when replacing raw (plates d/s data is more reliable if it exists)
  		
  		public double     bgnd_sure        = 0.18;  // minimal strength to be considered definitely background
  		public double     bgnd_maybe       = 0.1; // maximal strength to ignore as non-background
//  		public double     bgnd_2diff       = 0.005; // maximal strength to ignore as non-background
  		public int        min_clstr_seed   = 4; //2;     // number of tiles in a cluster to seed (just background?)
  		public int        min_clstr_lone   = 4;     // number of tiles in a cluster not close to other clusters (more than 2 tiles apart)
  		public double     min_clstr_weight = 0.0;   // Minimal total strength of the cluster 
  		public double     min_clstr_max    = 0.25;  // Minimal maximal strength of the cluster 
  		
  		public int        fill_gaps        = 4;     // same as in grow - 1:  4 directions by 1 step, 2: 8 directions by 1 step. +2*n - alternating hor/vert
  		public int        fill_final       = 50;    // same as fill_gaps, on the final pass
  		public int        min_clstr_block  = 3;     // number of tiles in a cluster to block (just non-background?)
  		public int        bgnd_grow        = 2;     // number of tiles to grow (1 - hor/vert, 2 - hor/vert/diagonal)
  		
//  		public double     ortho_min        = 0.09;  // minimal strength of hor/vert correlation to be used instead of full 4-pair correlation
  		public boolean    ortho_old        = false; // use old ortho features processing (orth0_* parameters, false - use or_*)
  		public double     ortho_min_hor    = 0.07;  // minimal strength of hor correlation to be used instead of full 4-pair correlation - 
  		public double     ortho_min_vert   = 0.15;  // minimal strength of vert correlation to be used instead of full 4-pair correlation
  		public double     ortho_asym       = 1.2;   // vert/hor (or hor/vert) strength to be used instead of the full correlation
  		public double     ortho_over4      = 0.8;   // vert/hor (or hor/vert) strength exceeding scaled 4-pair strength
  		public double     ortho_sustain    = 0.05;  // minimal strength of hor/vert to bridge over 
  		public int        ortho_run        = 3;     // minimal run of hor/vert tiles to be considered (at least from one side)
  		public double     ortho_minmax     = 0.09;  // minimal maximal strength in an ortho run 
  		public int        ortho_bridge     = 10;    // number of tiles to bridge over hor/vert gaps
  		public double     ortho_rms        = 0.3;   // maximal disparity RMS in a run to replace by average
  		public int        ortho_half_length = 4;    // convolve hor/vert strength by 3*(2*l+1) kernels to detect multi-tile features 
  		public double     ortho_mix        = 0.5;   // Fraction ovf convolved ortho in a mix with raw
  		
// Alternative mixing of ortho disparity/strength
  		public boolean    or_hor           = true;  // Apply ortho correction to horizontal correlation (vertical features)
  		public boolean    or_vert          = true;  // Apply ortho correction to vertical correlation (horizontal features)
  		public double     or_sigma         = 2.0;   // Blur sigma: verically for horizontal correlation, horizontally - for vertically
  		public double     or_sharp         = 0.0; // 0.5;   // 3-point sharpening (-k, +2k+1, -k)
  		public double     or_scale         = 2.5;   // Scale ortho correletion strength relative to 4-directional one
  		public double     or_offset        = 0.1;   // Subtract from scaled correlation strength, limit by 0
  		public double     or_asym          = 1.5;   // Minimal ratio of orthogonal strengths required for dis[parity replacement
  		public double     or_threshold     = 0.3; // 1.5;   // Minimal scaled offset ortho strength to normal strength needed for replacement
  		public double     or_absHor        = 0.15;  // Minimal horizontal absolute scaled offset ortho strength needed for replacement
  		public double     or_absVert       = 0.19;  // Minimal vertical absolute scaled offset ortho strength needed for replacement
  		
  		public boolean    poles_fix        = true;  // Continue vertical structures to the ground
  		public int        poles_len        = 25;    // Number of tiles to extend over the poles bottoms
  		public double     poles_ratio      = 1.0;   // Maximal ratio of invisible to visible pole length
  		public double     poles_min_strength = 0.1; // Set new pole segment strength to max of horizontal correlation and this value
  		public boolean    poles_force_disp = true;  // Set disparity to that of the bottom of existing segment (false - use hor. disparity)
  		
  		public int        max_clusters     = 500;   // Maximal number of clusters to generate for one run
  		public boolean    remove_scans     = true;  // Remove all unneeded scans when generating x3d output to save memory
  		public boolean    output_x3d       = true;  // Generate x3d output
  		public boolean    output_obj       = true;  // Generate Wavefront obj output
  		
  		
  		public boolean    correct_distortions = false; // Correct lens geometric distortions in a model (will need backdrop to be corrected too)
  		public boolean    show_triangles =    true;  // Show generated triangles
  		public boolean    avg_cluster_disp =  false;  // Weight-average disparity for the whole cluster 
  		public double     maxDispTriangle   = 0.2;    // Maximal relative disparity difference in a triangle face
  		public double     infinityDistance  = 10000;  // Distance to generate backdrop (0 - use regular backdrop)
  		public boolean    shUseFlaps        = true;  // Split into shells with flaps
  		public boolean    shAggrFade        = true;  // Aggressive fade alpha (whole boundary)
  		public int        shMinArea         = 1;     // Minimal shell area (not counting flaps
  		public double     shMinStrength    = 0.2;   // Minimal value of the shell maximum strength

  		// Thin ice parameters
  		public double     tiRigidVertical   = 3.0;   // 2.0 relative disparity rigidity in vertical direction
  		public double     tiRigidHorizontal = 1.0;   // 2.0 relative disparity rigidity in horizontal direction
  		public double     tiRigidDiagonal   = 0.5;   // 0.5 relative disparity rigidity in diagonal direction
  		public double     tiStrengthOffset  = 0.1;   // 0.1 strength "floor" - subtract (limit 0) before applying
  		public double     tiDispScale       = 0.5;   // 0.5 divide actual (disparity*eff_Strength) by this  before Math.pow and applying to T.I.
  		public double     tiDispPow         = 0.0;   // 0.0 apply pow to disparity (restore sign) for disparity difference pressure on ice
  		public double     tiDispPull        =  .1;   // 10.0 tiDispPull: multiply strength*disparity difference to pull force
  		public double     tiDispPullPreFinal=  .1;   // 5.0 Scale tiDispPull for pre-final pass
  		public double     tiDispPullFinal   =  .01;  // 2.0 Scale tiDispPull for final pass
  		
  		public double     tiBreakNorm       =   .5;  // Normalize stresses to average disparity if it is above threshold
  		public double     tiBreak3          = 0.6;   // 0.5 TI break value of abs(d0-3d1+3d2-d3)
  		public double     tiBreak31         = 0.1;   // 0.04 TI break value of (d0-3d1+3d2-d3) * (d1 - d2)
  		public double     tiBreak21         = 0.1;   // 0.1 TI break value of (-d0+d1+d2-d3) * abs(d1 - d2)
  		public double     tiBreakFar        = 0.3;   // 0.3  TI disparity threshold to remove as too far tiles
  		public double     tiBreakNear       = 0.3;   // 0.3 TI disparity threshold to remove as too near tiles
  		public int        tiBreakMode       = 0;     // 1 TI break mode: +1: abs(3-rd derivative), +2: -(3-rd * 1-st), +4: -(2-nd * abs(1-st)) , +8 - remove far, +16 - remove near 
  		public double     tiBreakSame       = 0.5;   // 0.75 Amplify colinear breaks in neighbor tiles
  		public double     tiBreakTurn       = 0.125; // 0.125 Amplify 90-degree turnintg breaks in neighbor tiles
  		
  		public double     tiHealPreLast     = 0.1;    // 0.1 Heal disparity gap before pre-last smooth
  		public double     tiHealLast        = 0.05;   // 0.05 Heal disparity gap before last smooth
  		public int        tiHealSame        = 5;      // 10 Maximal length of an internal break in the cluster to heal 
  		
  		public int        tiIterations      = 1000;  // 300 maximal number of TI iterations for each step
  		public int        tiPrecision       = 6;     // 6 iteration maximal error (1/power of 10) 
  		public int        tiNumCycles       = 5;     // 5 Number of cycles break-smooth (after the first smooth)
  		
  		// FG/BG separation
  		public boolean    stUseRefine =       false; // Apply super-tiles during refine passes 
  		public boolean    stUsePass2 =        true;  // Apply super-tiles during pass2 
  		public boolean    stUseRender =       true;  // Apply super-tiles during render
  		
  		public boolean    stShow =            false; // Calculate and show supertiles histograms 
  		public int        stSize            = 8;     // Super tile size (square, in tiles)
  		public double     stStepFar         = 0.1;   // Disparity histogram step for far objects 
  		public double     stStepNear        = 0.5;   // Disparity histogram step for near objects
  		public double     stStepThreshold   = 1.0;   // Disparity threshold to switch from linear to logarithmic steps
  		public double     stMinDisparity    = 0.0;   // Minimal disparity (center of a bin)
//  		public double     stMaxDisparity    = 15.0;  // Maximal disparity (center of a bin)
  		public double     stFloor           = 0.15;  // Subtract from strength, discard negative  
  		public double     stPow             = 1.0;   // raise strength to this power 
  		public double     stSigma           = 1.5;   // Blur disparity histogram (sigma in bins) 
  		public double     stMinBgDisparity  = 0.0;   // Minimal backgroubnd disparity to extract as a maximum from the supertiles 
  		public double     stMinBgFract      = 0.1;   // Minimal fraction of the disparity histogram to use as background 
  		public double     stUseDisp         = 0.15;  // Use background disparity from supertiles if tile strength is less 
  		public double     stStrengthScale   = 50.0;  // Multiply st strength if used instead of regular strength

  		public boolean    stSmplMode        = true;   // Use sample mode (false - regular tile mode)
  		public int        stSmplSide        = 2;      // Sample size (side of a square)
  		public int        stSmplNum         = 3;      // Number after removing worst
  		public double     stSmplRms         = 0.1;    // Maximal RMS of the remaining tiles in a sample
		public boolean    stSmplWnd         = false;  // Use window function for the samples (TODO: change default to true after testing)
  		
  		public int        stGrowSel         = 2;     // Grow initial selection before processing supertiles, odd - ortho. <0 - use all tiles
  		public int        stMeasSel         = 1;     // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
  		
  		public double     stSmallDiff       = 0.4;   // Consider merging initial planes if disparity difference below
  		public double     stHighMix         = 0.4;   // Consider merging initial planes if jumps between ratio above
  		
  		
  		public double     outlayerStrength  = 0.3;   // Outlayer tiles weaker than this may be replaced from neighbors
  		public double     outlayerDiff      = 0.4;   // Replace weak outlayer tiles that do not have neighbors within this disparity difference
  		public double     outlayerDiffPos   = 1.0;   // Replace weak outlayer tiles that have higher disparity than weighted average
  		public double     outlayerDiffNeg   = 0.4;   // Replace weak outlayer tiles that have lower disparity than weighted average
  		
  		// TODO: Make refine skip if already good?
  		public boolean    combine_refine    = true; // combine with all previous after refine pass
  		public double     combine_min_strength = 0.12; // Disregard weaker tiles when combining scans
  		public double     combine_min_hor =      0.12; // Disregard weaker tiles when combining scans for horizontal correlation
  		public double     combine_min_vert =     0.12; // Disregard weaker tiles when combining scans for vertical correlation
  		
// TODO: Move together with similar parameters  		
//  		public double     unique_tolerance = 0.1; // Do not re-measure correlation if target disparity differs from some previous by this
  		
  		// Multi-pass growing disparity
  		public int        grow_sweep         = 8; // Try these number of tiles around known ones 
  		public double     grow_disp_max =   160.0;// Maximal disparity to try
  		public double     grow_disp_trust =  4.0; // Trust measured disparity within +/- this value 
  		public double     grow_disp_step =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?  
  		public double     grow_min_diff =    0.5; // Grow more only if at least one channel has higher variance from others for the tile  
  		
  		public boolean    grow_retry_far =  false; // Retry tiles around known foreground that have low max_tried_disparity
  		public boolean    grow_pedantic =   false; // Scan full range between max_tried_disparity of the background and known foreground
		public boolean    grow_retry_inf =  false;  // Retry border tiles that were identified as infinity earlier
  		
		// New for initial growing
		public boolean    gr_new_expand        =   true;
		public int        gr_max_expand        = 500; // 150; // 30;
		public double     gr_ovrbg_cmb         =   0.3; // 0.3; 
		public double     gr_ovrbg_cmb_hor     =   0.3; // 0.3;
		public double     gr_ovrbg_cmb_vert    =   0.3; // 0.3;
		public double     gr_ovrbg_filtered    =   0.3; // 0.3; 
		  
		public double     fds_str_floor        =   0.09; // Should be normally less than combine_min_strength 
		public double     fds_str_pow          =   1.0;
		public int        fds_smpl_side        =   5; // 3;      // Sample size (side of a square)
		public int        fds_smpl_num         =  13; // 13; // 5;      // Number after removing worst (should be >1)
		public double     fds_smpl_rms         =   0.15; // Maximal RMS of the remaining tiles in a sample
		public double     fds_smpl_rel_rms     =   0.01; // 05;  // Maximal RMS/disparity in addition to smplRms
		public boolean    fds_smpl_wnd         =   true; //
		public double     fds_abs_tilt         =   2.0; // pix per tile
		public double     fds_rel_tilt         =   0.2; // (pix / disparity) per tile
		
// Macro disparity scanning parameters		
		public double     mc_disp8_step        =   2.0;   // Macro disparity scan step (actual disparity step is 8x)
		public double     mc_disp8_trust       =   2.0;   //Trust measured disparity within +/- this value 
		public double     mc_strength          =   0.2;   // Minimal composite correlation to process (0.2..0.3)
 		public double     mc_unique_tol        =   0.05;  // Do not re-measure macro correlation if target disparity differs from some previous by this
 		public double     mc_trust_fin         =   0.3;   // When consolidating macro results, exclude high residual disparity
 		public double     mc_trust_sigma       =   0.2;   // Gaussian sigma to reduce weight of large residual disparity
 		public double     mc_ortho_weight      =   0.5;   // Weight from ortho neighbor supertiles
 		public double     mc_diag_weight       =   0.25;  // Weight from diagonal neighbor supertiles
 		public double     mc_gap               =   0.4;   // Do not remove measurements farther from the kept ones
		
//		  0x1e, // 0x1f, // final int         variants_mask,
		public int        gr_min_new           =  20;    // Discard variant if it requests too few tiles
		public boolean    gr_var_new_sngl      =   false;// Expand only unambiguous tiles over previously undefined  
		public boolean    gr_var_new_fg        =   true; // Expand unambiguous and foreground tiles over previously undefined
		public boolean    gr_var_all_fg        =   true;  
		public boolean    gr_var_new_bg        =   false;  
		public boolean    gr_var_all_bg        =   false;  
		public boolean    gr_var_next          =   false; // try next disparity range  TODO: add related statements
		public int        gr_num_steps         =   8;    // How far to extend over previously undefined disparity tiles
		public int        gr_steps_over        =   4;    // How far to extend over previously determined disparity tiles
		public int        gr_smpl_size         =   5;    // Extend sample square side
		public int        gr_min_pnts          =   3;    // Extend at least this number of the seed tiles
		public boolean    gr_use_wnd           =   true; // Use window function for square sample 
		public double     gr_tilt_damp         =   0.001;// Tilt cost for damping insufficient plane data 
		public double     gr_split_rng         =   5.0;  // When growing, range of disparities to be extended without far/near division
		public double     gr_same_rng          =   3.0;  // consider far/near tiles within that range from the farthest/closest
		public double     gr_diff_cont         =   2.0;  // Maximal difference from the old value when smoothing
		public double     gr_abs_tilt          =   2.0;  // Maximal filter disparity absolute tilt (pix per tile)
		public double     gr_rel_tilt          =   0.2;  // Maximal filter disparity tilt (pix / disparity) per tile
		public int        gr_smooth            =   50;   // Maximal number of smoothing steps (reduce if long?)
		public double     gr_fin_diff          =   0.01; // Maximal change to finish smoothing iterations 
 		public double     gr_unique_tol        =   0.15; // Do not re-measure correlation if target disparity differs from some previous by this
	 	public double     gr_unique_pretol     =   0.5;  // Larger tolerance for expanding (not refining)
		
  		public boolean    plPreferDisparity    =   false;// Always start with disparity-most axis (false - lowest eigenvalue)
  		public double     plDispNorm           =   5.0;  // Normalize disparities to the average if above (now only for eigenvalue comparison)
  		
  		public double     plBlurBinVert        =   1.2;  // Blur disparity histograms for constant disparity clusters by this sigma (in bins)
  		public double     plBlurBinHor         =   0.8;  // Blur disparity histograms for horizontal clusters by this sigma (in bins)
  		public double     plMaxDiffVert        =   0.4;  // Maximal normalized disparity difference when initially assigning to vertical plane
  		public double     plMaxDiffHor         =   0.2;  // Maximal normalized disparity difference when initially assigning to horizontal plane
  		public int        plInitPasses         =     3;  // Number of initial passes to assign tiles to vert (const disparity) and hor planes
  		
  		public int        plMinPoints          =     5;  // Minimal number of points for plane detection
  		public double     plTargetEigen        =   0.02; // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
  		public double     plFractOutliers      =   0.3;  // Maximal fraction of outliers to remove
  		public int        plMaxOutliers        =    20;  // Maximal number of outliers to remove
  		public double     plMinStrength        =   0.01; // Minimal total strength of a plane 
  		public double     plMaxEigen           =   0.06; // Maximal eigenvalue of a plane 
  		public double     plEigenFloor         =   0.005;// Add to eigenvalues of each participating plane and result to validate connections 
  		public double     plEigenStick         =   25.0; // Consider plane to be a "stick" if second eigenvalue is below 
  		public double     plBadPlate           =   0.2;  // Not a plate if sin^2 between normals from disparity and world exceeds this 
  		public boolean    plDbgMerge           =   true; // Combine 'other' plane with current
  		public double     plWorstWorsening     =   2.0;  // Worst case worsening after merge
  		public double     plWorstWorsening2    =   5.0;  // Worst case worsening for thin planes
  		public double     plWorstEq            =   1.0;  // Worst case worsening after merge with equal weights
  		public double     plWorstEq2           =   2.0;  // Worst case worsening for thin planes with equal weights
  		public double     plOKMergeEigen       =   0.03; // If result of the merged planes is below, OK to use thin planes (higher) threshold 
  		public double     plMaxWorldSin2       =   0.1;  // Maximal sine squared of the world angle between planes to merge. Set to >= 1.0 to disable
  		public double     pl2dForSin           =   2.5;  // Do not compare sin() between planes, if at least one has too small axis ratio 
  		public double     plWeakWorsening      =   1.0;  // Relax merge requirements for weaker planes
  		public double     plMaxOverlap         =   0.1;  // Maximal overlap between the same supertile planes to merge
  		// Merge same supetile planes if at least one is weak and they do not differ much
  		public double     plWeakWeight         =   0.2 ; // Maximal weight of the weak plane to merge
  		public double     plWeakEigen          =   0.1;  // Maximal eigenvalue of the result of non-weighted merge
  		public double     plWeakWeight2        =  10.0 ; // Maximal weight of the weak plane to merge (second variant)
  		public double     plWeakEigen2         =   0.05; // Maximal eigenvalue of the result of non-weighted merge  (second variant)
  		public double     plSumThick           =   1.6;  // Do not merge if any sqrt of merged eigenvalue exceeds scaled sum of components 
  		public double     plNeNeibCost         =   5.0;  // When calculating non-exclusive planes, do not use neighbors with high cost
  		public double     plNeOwn              =   5.0;  // When calculating non-exclusive planes, use center plane relative weight

  		public double     plExNeibCost         =   5.0;  // When calculating exclusive planes links, do not use neighbors with high cost
//  	    public double     plExNeibCostSngl     =  10.0;  // When calculating exclusive planes links, do not use no-link neighbors with high cost
  	    public double     plExNeibSmooth       =   0.5;  // Scale down maximal costs for smoothed planes (tighter requirements) 
  	    public double     plMergeCostStar      =   5.0;  // Cost threshold for merging same tile planes if the plane has connected neighbors
  	    public double     plMergeCost          =  10.0;  // Cost threshold for merging same tile planes if not connected

  	    public boolean    plConflMerge         =  true;  // Try to merge conflicting planes
  	    public double     plConflRelax         =   1.5;  // Scale parameters to relax planes fit for merging conflicting planes
  	    public boolean    plConflSngl          =  true;  // Only merge conflicting planes if this is the only conflicting pair in the supertile
  	    public boolean    plConflSnglPair      =  true;  // Only merge conflicting planes only if there are just two planes in the supertile
  	    
  	    public double     plWeakFgStrength     =   0.15; // Consider merging plane if it is foreground and maximal strength below this  
  	    public int        plWeakFgOutliers     =   1;    // Remove these strongest from foreground when determining the maximal strength
  	    public double     plWeakFgRelax        =   2.0;  // Relax cost requirements when merging with weak foreground   
  	    
  	    
  	    public double     plThickWorld         =   0.2;  // Maximal real-world thickness of merged overlapping planes (meters) 
  	    public double     plThickWorldConfl    =   0.4;  // Maximal real-world merged thickness for conflicting planes 
  	    public double     plRelaxComplete      =   1.5;  // Relax cost requirements when adding exclusive links to complete squares and triangles 
  	    public double     plRelaxComplete2     =   2.0;  // Relax cost requirements during the second pass 
  	    
  	    
  		public double     plMaxZRatio          =   2.0;  // Maximal ratio of Z to allow plane merging
  		public double     plMaxDisp            =   0.6;  // Maximal disparity of one of the planes to apply  maximal ratio
  		public double     plCutTail            =   1.4;  // When merging with neighbors cut the tail that is worse than scaled best
  		public double     plMinTail            =   0.015;// Set cutoff value level not less than
  		
  		
  		// parameters to recreate planes from tiles disparity/strengths using determined plane connections to neighbors
  		public boolean    plDiscrEn            =   true; // Enable planes tiles selection regeneration hinted by supertile neighbors
  		public double     plDiscrTolerance     =   0.4;  // Maximal disparity difference from the plane to consider tile 
  		public double     plDiscrDispRange     =   1.0;  // Parallel move known planes around original know value for the best overall fit
  		public int        plDiscrSteps         =   10;   // Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)
  		public int        plDiscrVariants      =   100;  // total number of variants to try (protect from too many planes) 
  		public int        plDiscrMode          =   3;    // 0 - weighted, 1 - equalized, 2 - best, 3 - combined
  		
  		public double     plDiscrVarFloor      =   0.03;  // Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)
  		public double     plDiscrSigma         =   0.05;  // Gaussian sigma to compare how measured data is attracted to planes
  		public double     plDiscrBlur          =   0.05;  // Sigma to blur histograms while re-discriminating
  		public double     plDiscrExclusivity   =   1.5;   // Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others
  		public double     plDiscrExclus2       =   0.8;   // For second pass if exclusivity > 1.0 - will assign only around strong neighbors 
  		public boolean    plDiscrStrict        =   false; // When growing selection do not allow any offenders around (false - more these than others)
  		public double     plDiscrCorrMax       =   0.7;   // Attraction to different planes correlation that is too high for re-discrimination.
  		public double     plDiscrCorrMerge     =   0.85;  // Attraction to different planes correlation that is high enough to merge planes				
  		public int        plDiscrSteal         =   4;     // If offender has this number of tiles (including center) the cell can not be used
  		public int        plDiscrGrown         =   4;     // Only use tiles within this range from original selection
  		
  		public double     plDiscrXMedian       =   1.5;   // Remove outliers from the final selection that have distance more than scaled median

  		// comparing merge quality for plane pairs
  		public double     plCostDist           =   4.0;  // Disparity (pix) - closer cost will use more of the real world, farther - disparity
  		public double     plCostKrq            =   0.8;  // Cost of merge quality sqrt(weighted*equal) in disparity space
  		public double     plCostKrqEq          =   0.2;  // Cost of merge quality averaje of weighted and equal weight in disparity space
  		public double     plCostWrq            =   0.8;  // Cost of merge quality sqrt(weighted*equal) in world space
  		public double     plCostWrqEq          =   0.2;  // Cost of merge quality average of weighted and equal weight in world space
  		public double     plCostSin2           =  10.0;  // Cost of sin squared between normals
  		public double     plCostRdist2         =1000.0;  // Cost of squared relative distances       
  		
  		
  		public boolean    plConflDualTri       =   false; // Resolve dual triangles conflict (odoodo)
  		public boolean    plConflMulti         =   false; // Resolve multiple odo triangles conflicts
  		public boolean    plConflDiag          =   false; // Resolve diagonal (ood) conflicts
  		public boolean    plConflStar          =   true;  // Resolve all conflicts around a supertile 

  		public int        plStarSteps          =   2;    // How far to look around when calculating connection cost
  		public double     plStarOrtho          =   0.5;  // When calculating cost for the connections scale 4 ortho neighbors
  		public double     plStarDiag           =   0.25; // When calculating cost for the connections scale 4 diagonal neighbors
  		public double     plStarPwr            =   0.5;  // Divide cost by number of connections to this power
  		public double     plStarWeightPwr      =   0.5;  // use this power of tile weight when calculating connection cost
  		public double     plWeightToDens       =   0.3;  // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
  		public double     plStarValPwr         =   1.0;  // Raise value of each tile before averaging
  		public double     plDblTriLoss         =   0.0001; // When resolving double triangles allow minor degradation (0.0 - strict)
  		public boolean    plNewConfl           =   false; // Allow more conflicts if overall cost is reduced
  		public int        plMaxChanges         =   0;     // Maximal number of simultaneous connection changes around one tile (0 - any)
  		
  		public boolean    plMutualOnly         =   true; // keep only mutual links, remove weakest if conflict
  		public boolean    plFillSquares        =   true; // Add diagonals to full squares
  		public boolean    plCutCorners         =   true; // Add ortho to 45-degree corners
  		public boolean    plHypotenuse         =   true; // Add hypotenuse connection if both legs exist

  		public double     plPull               =  5.0; // .3;   // Relative weight of original (measured) plane compared to average neighbor pull
  		                                                // when combing with neighbors
  		public int        plIterations         =  10;   // Maximal number of smoothing iterations for each step
  		public boolean    plStopBad            =  true; // Do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
  		public int        plPrecision          =  6;    // Maximal step difference (1/power of 10)
  		
  		public double     plSplitPull          =  .5;   // Relative weight of center plane when splitting into pairs
  		public double     plNormPow            =  .5;   // 0.0: 8 neighbors pull 8 times as 1, 1.0 - same as 1
  		public int        plSplitMinNeib       =  2;    // Minimal number of neighbors to split plane in pairs
  		public double     plSplitMinWeight     =  2.0;  // Minimal weight of split plains to show
  		public double     plSplitMinQuality    =  1.1;  // Maximal normalized disparity difference from the plane to consider
  		
  		public boolean    plSplitApply         = true;  // Apply plane split to pairs
  		public boolean    plNonExclusive       = true;  // Allow tiles to belong to both planes of the pair
  		public boolean    plUseOtherPlanes     = false; // Allow other tiles from the same supertile
  		public boolean    plAllowParallel      = true;  // Allow parallel shift of the specified planes before adding
  		public double     plMaxDiff            = 0.3;   // Maximal normalized tile disparity difference from the plane to consider
  		public double     plOtherDiff          = 1.4;   // Maximal difference of the added tile ratio to the average  disparity difference
  		
  		
  		public boolean    plSplitXY            = true;  // Separate tiles for split planes by X, Y 
  		public double     plSplitXYTolerance   = 0.2;   // Disparity tolerance when separating by X, Y 
  		
  		
  		public boolean    plFuse               =  true; // Fuse planes together (off for debug only)
  		public boolean    plKeepOrphans        =  true; // Keep unconnected supertiles
  		public double     plMinOrphan          =  2.0;  // Minimal strength unconnected supertiles to keep
  		
  		public double     plSnapDispAny        =  .2;   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength
  		public double     plSnapStrengthAny    =  .2;   // Maximal strength to fit any distance (if does not fit otherwise - treat as zero str4ength
  		public double     plSnapNegAny         =  .5;   // Maximal negative disparity difference from the best match
  		public double     plSnapDispMax        =  .5;   // Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength
  		public double     plSnapDispWeight     =  .5;   // Maximal disparity diff. by weight product to snap to plane
  		public int        plSnapZeroMode       =  1;    // Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest
  		
  		public boolean    msUseSel             =   true; // Use planes selection masks (generated when splitting to intersecting pairs  
  		public boolean    msDivideByArea       =   true; // Divide plane strengths by ellipsoid area
  		public double     msScaleProj          =   1.5;  // Scale projection of the plane ellipsoid 
  		public double     msFractUni           =   0.3;  // Spread this fraction of the ellipsoid weight among extended (double) supertile
  		
		public boolean    tsNoEdge             = true;   // Do not assign tiles to the surface edges (not having all 8 neighbors)
		public boolean    tsUseCenter          = true;   // Only assign outside of 8x8 center if no suitable alternative
  		public double     tsMaxDiff            = 0.3;    // Maximal disparity difference when assigning tiles
		public double     tsMinDiffOther       = 0.35;   // Minimal disparity difference to be considered as a competitor surface
		public double     tsMinStrength        = 0.05;   // Minimal tile correlation strength to be assigned
		public double     tsMaxStrength        = 10.0;   // Maximal tile correlation strength to be assigned
		public double     tsMinSurface         = 0.001;  // Minimal surface strength at the tile location
		public int        tsMoveDirs           = 3;      // Allowed tile disparity correction: 1 increase, 2 - decrease, 3 - both directions
		public double     tsSurfStrPow         = 0.0;    // Raise surface strengths ratio to this power when comparing candidates
		public double     tsAddStrength        = 0.01;   // Add to strengths when calculating pull of assigned tiles
		public double     tsSigma              = 2.0;    // Radius of influence (in tiles) of the previously assigned tiles
		public double     tsNSigma             = 2.0;    // Maximal relative to radius distance to calculate influence
		public double     tsMinPull            = 0.001;  // Additional pull of each surface 
		public double     tsMinAdvantage       = 3.0;    // Minimal ratio of the best surface candidate to the next one to make selection
		
		public int        tsClustSize          = 1;      // Minimal size of a cluster to keep
		public double     tsClustWeight        = 0.2;    // Minimal total weight of a cluster to keep
		
		public int        tsMinNeib            = 6;      // Minimal number of neighbors of unassigned tile to join (the farthest)
		public double     tsMaxSurStrength     = 0.05;   // Maximal strength of the surrounded unassigned tile to join
		public boolean    tsCountDis           = true;   // Include disabled tiles/borders when counting assigned neighbors
		
		public boolean    tsEnPlaneSeed        = true;   // Assign tiles that were used to generate planes
		public boolean    tsEnOnly             = true;   // Allow assignment only surface
		public boolean    tsEnGrow             = true;   // Grow the only surface assignments
		public double     tsGrowStrength       = 0.01;   // Maximal strength when growing the only surfaces
		public boolean    tsGrowStrong         = true;   // Grow over strong if disparity matches 
		public double     tsContStrength       = 0.1;    // Minimal strength to continue grow with disparity match 
		public double     tsContDiff           = 0.1;    // Maximal normalized disparity error to grow over strong tiles 
		
		
		public boolean    tsEnSingle           = true;   // Allow assignment to the nearest surface with no competitors
		public boolean    tsEnMulti            = true;   // Allow assignment when several surfaces fit

		public boolean    tsRemoveWeak1        = false;  // Remove weak clusters before growing
		public boolean    tsGrowSurround       = true;   // Assign tiles that have neighbors to the lowest disparity
		public boolean    tsRemoveWeak2        = true;   // Remove weak clusters after growing
		
		public boolean    tsLoopMulti          = true;   // Repeat multi-choice assignment while succeeding 
		public boolean    tsReset              = false;  // Reset tiles to surfaces assignment
		public boolean    tsShow               = false;  // Show results of tiles to surfaces assignment
		public int        tsNumClust           = 500;     // Number of clusters to keep

		public int        tsConsensMode        = 7;      // Which assignments to match +1 - combo, +2 grown single, +4 plane seeds 
		public int        tsConsensAgree       = 1;      // Minimal number of assignments to agree
		
		// Tile assignment parameters
		public double     taMinFgBg            = 0.1;    // Minimal foreground/ background separation to look for weak FG edge
		public double     taMinFgEdge          = 0.2;    // Minimal foreground edge strength (stronger edges will have proportionally smaller costs)
		public double     taMinColSep          = 0.05;   // Minimal surface separation that requires color change 
		public double     taMinColDiff         = 0.01;   // Minimal color variation (larger proportionally reduces cost) 
		public double     taOutlier            = 1.0;    // Disparity difference limit
		public double     taDiffPwr            = 0.25;   // Strength power when calculating disparity error
		public double     taBestPwr            = 0.0;    // Strength power when calculating disparity error over best
		public double     taDiff9Pwr           = 0.5;    // Strength power when calculating disparity error for group of 9
		public double     taColSigma           = 1.5;    // Gaussian sigma to blur color difference between tiles along each direction
		public double     taColFraction        = 0.3;    // Relative amount of the blurred color difference in the mixture  
		
		
  		public double     taCostEmpty          = 1.0;    // Cost of a tile that is not assigned
  		public double     taCostNoLink         = 1.0;    // Cost of a tile not having any neighbor in particular direction
  		public double     taCostSwitch         = 1.0;    // Cost of a tile switching to a neighbor that does not have a link
  		public double     taCostColor          = 1.0;    // Cost of a tile switching to a disconnected neighbor divided by a color mismatch
  		public double     taCostDiff           = 1.0;    // Cost of a weighted normalized tile disparity error
  		public double     taCostDiffBest       = 1.0;    // Cost of a weighted normalized tile disparity error above best surface
  		public double     taCostDiff9          = 1.0;    // Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)
  		public double     taCostWeakFgnd       = 1.0;    // Cost of a weak foreground edge
  		public double     taCostFlaps          = 1.0;    // Cost of using supertile "flaps" (not in the center 8x8 tiles area)
  		public double     taCostMismatch       = 1.0;    // Cost of a measurement layer not having same layer in the same location or near
  		
  		public boolean    taEnEmpty            = true;   // Enable cost of a tile that is not assigned
  		public boolean    taEnNoLink           = true;   // Enable cost of a tile not having any neighbor in particular direction
  		public boolean    taEnSwitch           = true;   // Enable cost of a tile switching to a neighbor that does not have a link
  		public boolean    taEnColor            = true;   // Enable cost of a tile switching to a disconnected neighbor divided by a color mismatch
  		public boolean    taEnDiff             = true;   // Enable cost of a weighted normalized tile disparity error
  		public boolean    taEnDiffBest         = true;   // Enable cost of a weighted normalized tile disparity error above best surface
  		public boolean    taEnDiff9            = true;   // Enable cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)
  		public boolean    taEnWeakFgnd         = true;   // Enable cost of a weak foreground edge
  		public boolean    taEnFlaps            = true;   // Enable cost of using supertile "flaps" (not in the center 8x8 tiles area)
  		public boolean    taEnMismatch         = false;  // Enable cost of a measurement layer not having same layer in the same location or near
  		
		
		
  		public boolean    replaceWeakOutlayers =   true; // false; 
  		
  		public boolean    dbg_migrate =            true; 
  		
  		// other debug images
  		public boolean    show_ortho_combine =     false; // Show 'ortho_combine' 
  		public boolean    show_refine_supertiles = false; // show 'refine_disparity_supertiles' 
  		public boolean    show_bgnd_nonbgnd =      false; // show 'bgnd_nonbgnd' 
  		public boolean    show_filter_scan =       false; // show 'FilterScan'
  		public boolean    show_combined =          false; // show 'combo_scan' (combined multiple scans)
  		public boolean    show_unique =            false; // show 'unique_scan' (removed already measured tiles with the same disparity)
  		public boolean    show_histograms =        false; // show supertile disparity histograms 
  		public boolean    show_init_refine =       false; // show debug images during initial refinement 
  		public boolean    show_expand =            false; // show debug images during disparity expansion
  		public boolean    show_variant =           false; // show prepareExpandVariant when elevating variant number
  		public boolean    show_retry_far =         false; // show debug images related to retrying far tiles near foreground
  		public boolean    show_macro =             false; // show debug images related to macro correlation
  		
  		public boolean    show_shells =            false; // show 'shells' 
  		public boolean    show_neighbors =         false; // show 'neighbors' 
  		public boolean    show_flaps_dirs =        false; // show 'flaps-dirs' 
  		public boolean    show_first_clusters =    false; // show 'first_N_clusters' 
  		public boolean    show_planes =            false; // show planes
  		public double []  vertical_xyz =           {0.0,1.0,0.0}; // real world up unit vector in camera CS (x - right, y - up, z - to camera};
  		
  		public CLTParameters(){}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"transform_size",   this.transform_size+"");
  			properties.setProperty(prefix+"clt_window",       this.clt_window+"");
  			properties.setProperty(prefix+"shift_x",          this.shift_x+"");
  			properties.setProperty(prefix+"shift_y",          this.shift_y+"");
  			properties.setProperty(prefix+"iclt_mask",        this.iclt_mask+"");
  			properties.setProperty(prefix+"tileX",            this.tileX+"");
  			properties.setProperty(prefix+"tileY",            this.tileY+"");
  			properties.setProperty(prefix+"dbg_mode",         this.dbg_mode+"");
  			properties.setProperty(prefix+"ishift_x",         this.ishift_x+"");
  			properties.setProperty(prefix+"ishift_y",         this.ishift_y+"");
  			properties.setProperty(prefix+"fat_zero",         this.fat_zero+"");
  			properties.setProperty(prefix+"corr_sigma",       this.corr_sigma+"");
			properties.setProperty(prefix+"norm_kern",        this.norm_kern+"");
			properties.setProperty(prefix+"gain_equalize",    this.gain_equalize+"");
			properties.setProperty(prefix+"colors_equalize",  this.colors_equalize+"");
  			properties.setProperty(prefix+"novignetting_r",   this.novignetting_r+"");
  			properties.setProperty(prefix+"novignetting_g",   this.novignetting_g+"");
  			properties.setProperty(prefix+"novignetting_b",   this.novignetting_b+"");
  			properties.setProperty(prefix+"scale_r",          this.scale_r+"");
  			properties.setProperty(prefix+"scale_g",          this.scale_g+"");
  			properties.setProperty(prefix+"scale_b",          this.scale_b+"");
  			properties.setProperty(prefix+"vignetting_max",   this.vignetting_max+"");
  			properties.setProperty(prefix+"vignetting_range", this.vignetting_range+"");
  			properties.setProperty(prefix+"kernel_step",      this.kernel_step+"");
  			properties.setProperty(prefix+"disparity",        this.disparity +"");
  			properties.setProperty(prefix+"z_correction",     this.z_correction +"");
			properties.setProperty(prefix+"correlate",        this.correlate+"");
  			properties.setProperty(prefix+"corr_mask",        this.corr_mask+"");
			properties.setProperty(prefix+"corr_sym",         this.corr_sym+"");
			properties.setProperty(prefix+"corr_keep",        this.corr_keep+"");
			properties.setProperty(prefix+"corr_show",        this.corr_show+"");
			properties.setProperty(prefix+"corr_mismatch",    this.corr_mismatch+"");
  			properties.setProperty(prefix+"corr_offset",      this.corr_offset +"");
  			properties.setProperty(prefix+"corr_red",         this.corr_red +"");
  			properties.setProperty(prefix+"corr_blue",        this.corr_blue +"");
			properties.setProperty(prefix+"corr_normalize",   this.corr_normalize+"");
  			properties.setProperty(prefix+"min_corr",         this.min_corr +"");
  			properties.setProperty(prefix+"min_corr_normalized",this.min_corr_normalized +"");
  			properties.setProperty(prefix+"max_corr_sigma",   this.max_corr_sigma +"");
  			properties.setProperty(prefix+"max_corr_radius",  this.max_corr_radius +"");
  			
  			properties.setProperty(prefix+"enhortho_width",   this.enhortho_width +"");
  			properties.setProperty(prefix+"enhortho_scale",   this.enhortho_scale +"");

  			
  			properties.setProperty(prefix+"max_corr_double",  this.max_corr_double+"");
  			properties.setProperty(prefix+"corr_mode",        this.corr_mode+"");
  			properties.setProperty(prefix+"corr_border_contrast", this.corr_border_contrast +"");
  			properties.setProperty(prefix+"tile_task_op",     this.tile_task_op+"");
  			properties.setProperty(prefix+"tile_task_wl",     this.tile_task_wl+"");
  			properties.setProperty(prefix+"tile_task_wt",     this.tile_task_wt+"");
  			properties.setProperty(prefix+"tile_task_ww",     this.tile_task_ww+"");
  			properties.setProperty(prefix+"tile_task_wh",     this.tile_task_wh+"");
  			properties.setProperty(prefix+"min_shot",       this.min_shot +"");
  			properties.setProperty(prefix+"scale_shot",       this.scale_shot +"");
  			properties.setProperty(prefix+"diff_sigma",       this.diff_sigma +"");
  			properties.setProperty(prefix+"diff_threshold",   this.diff_threshold +"");
			properties.setProperty(prefix+"diff_gauss",       this.diff_gauss+"");
  			properties.setProperty(prefix+"min_agree",        this.min_agree +"");
			properties.setProperty(prefix+"dust_remove",      this.dust_remove+"");
			properties.setProperty(prefix+"black_back",       this.black_back+"");
			properties.setProperty(prefix+"keep_weights",     this.keep_weights+"");
			properties.setProperty(prefix+"sharp_alpha",      this.sharp_alpha+"");

			properties.setProperty(prefix+"alpha0",           this.alpha0 +"");
  			properties.setProperty(prefix+"alpha1",           this.alpha1 +"");

  			properties.setProperty(prefix+"gen_chn_stacks",   this.gen_chn_stacks+"");
			properties.setProperty(prefix+"gen_chn_img",      this.gen_chn_img+"");
			properties.setProperty(prefix+"gen_4_img",        this.gen_4_img+"");
			properties.setProperty(prefix+"show_nonoverlap",  this.show_nonoverlap+"");
			properties.setProperty(prefix+"show_overlap",     this.show_overlap+"");
			properties.setProperty(prefix+"show_rgba_color",  this.show_rgba_color+"");
			properties.setProperty(prefix+"show_map",         this.show_map+"");
  			properties.setProperty(prefix+"disp_scan_start",  this.disp_scan_start +"");
  			properties.setProperty(prefix+"disp_scan_step",   this.disp_scan_step +"");
  			properties.setProperty(prefix+"disp_scan_count",  this.disp_scan_count+"");

			properties.setProperty(prefix+"fine_dbg",         this.fine_dbg+"");
  			properties.setProperty(prefix+"fine_corr_x_0",    this.fine_corr_x_0 +"");
  			properties.setProperty(prefix+"fine_corr_y_0",    this.fine_corr_y_0 +"");
  			properties.setProperty(prefix+"fine_corr_x_1",    this.fine_corr_x_1 +"");
  			properties.setProperty(prefix+"fine_corr_y_1",    this.fine_corr_y_1 +"");
  			properties.setProperty(prefix+"fine_corr_x_2",    this.fine_corr_x_2 +"");
  			properties.setProperty(prefix+"fine_corr_y_2",    this.fine_corr_y_2 +"");
  			properties.setProperty(prefix+"fine_corr_x_3",    this.fine_corr_x_3 +"");
  			properties.setProperty(prefix+"fine_corr_y_3",    this.fine_corr_y_3 +"");

  			properties.setProperty(prefix+"fcorr_radius",     this.fcorr_radius +"");
  			properties.setProperty(prefix+"fcorr_min_strength",this.fcorr_min_strength +"");
  			properties.setProperty(prefix+"fcorr_disp_diff",  this.fcorr_disp_diff +"");
			properties.setProperty(prefix+"fcorr_quadratic",  this.fcorr_quadratic+"");
			properties.setProperty(prefix+"fcorr_ignore",     this.fcorr_ignore+"");

  			properties.setProperty(prefix+"fcorr_inf_strength",this.fcorr_inf_strength +"");
  			properties.setProperty(prefix+"fcorr_inf_diff",   this.fcorr_inf_diff +"");
			properties.setProperty(prefix+"fcorr_inf_quad",   this.fcorr_inf_quad+"");
			properties.setProperty(prefix+"fcorr_inf_vert",   this.fcorr_inf_vert+"");
			
			properties.setProperty(prefix+"inf_disp_apply",   this.inf_disp_apply+"");
  			properties.setProperty(prefix+"inf_repeat",       this.inf_repeat+"");
  			
//			properties.setProperty(prefix+"inf_mism_apply",   this.inf_mism_apply+"");
  			properties.setProperty(prefix+"inf_iters",        this.inf_iters+"");
  			properties.setProperty(prefix+"inf_final_diff",   this.inf_final_diff +"");
  			properties.setProperty(prefix+"inf_far_pull",     this.inf_far_pull +"");
  			properties.setProperty(prefix+"inf_str_pow",      this.inf_str_pow +"");
  			properties.setProperty(prefix+"inf_smpl_side",    this.inf_smpl_side+"");
  			properties.setProperty(prefix+"inf_smpl_num",     this.inf_smpl_num+"");
  			properties.setProperty(prefix+"inf_smpl_rms",     this.inf_smpl_rms +"");

  			properties.setProperty(prefix+"ih_smpl_step",     this.ih_smpl_step+"");
  			properties.setProperty(prefix+"ih_disp_min",      this.ih_disp_min +"");
  			properties.setProperty(prefix+"ih_disp_step",     this.ih_disp_step +"");
  			properties.setProperty(prefix+"ih_num_bins",      this.ih_num_bins+"");
  			properties.setProperty(prefix+"ih_sigma",         this.ih_sigma +"");
  			properties.setProperty(prefix+"ih_max_diff",      this.ih_max_diff +"");
  			properties.setProperty(prefix+"ih_min_samples",   this.ih_min_samples+"");
			properties.setProperty(prefix+"ih_norm_center",   this.ih_norm_center+"");

			properties.setProperty(prefix+"ly_smpl_side",     this.ly_smpl_side+"");
			properties.setProperty(prefix+"ly_smpl_num",      this.ly_smpl_num+"");
			properties.setProperty(prefix+"ly_meas_disp",     this.ly_meas_disp +"");
			properties.setProperty(prefix+"ly_smpl_rms",      this.ly_smpl_rms +"");
			properties.setProperty(prefix+"ly_disp_var",      this.ly_disp_var +"");
			properties.setProperty(prefix+"ly_inf_frac",      this.ly_inf_frac +"");
			properties.setProperty(prefix+"ly_on_scan",       this.ly_on_scan+"");
			properties.setProperty(prefix+"ly_inf_en",        this.ly_inf_en+"");
			properties.setProperty(prefix+"ly_inf_force",     this.ly_inf_force+"");
			properties.setProperty(prefix+"ly_poly",          this.ly_poly+"");

			properties.setProperty(prefix+"corr_magic_scale", this.corr_magic_scale +"");
  			
			properties.setProperty(prefix+"show_textures",    this.show_textures+"");
			properties.setProperty(prefix+"debug_filters",    this.debug_filters+"");

			properties.setProperty(prefix+"min_smth",         this.min_smth +"");
			properties.setProperty(prefix+"sure_smth",        this.sure_smth +"");
			properties.setProperty(prefix+"bgnd_range",       this.bgnd_range +"");
			properties.setProperty(prefix+"other_range",      this.other_range +"");

			properties.setProperty(prefix+"ex_strength",      this.ex_strength +"");
			properties.setProperty(prefix+"ex_nstrength",     this.ex_nstrength +"");
			
			properties.setProperty(prefix+"ex_over_bgnd",     this.ex_over_bgnd+"");
			properties.setProperty(prefix+"ex_min_over",      this.ex_min_over +"");

			properties.setProperty(prefix+"pt_super_trust",   this.pt_super_trust +"");
			properties.setProperty(prefix+"pt_keep_raw_fg",   this.pt_keep_raw_fg+"");
			properties.setProperty(prefix+"pt_scale_pre",     this.pt_scale_pre +"");
			properties.setProperty(prefix+"pt_scale_post",    this.pt_scale_post +"");

			properties.setProperty(prefix+"bgnd_sure",        this.bgnd_sure +"");
			properties.setProperty(prefix+"bgnd_maybe",       this.bgnd_maybe +"");
  			properties.setProperty(prefix+"min_clstr_seed",   this.min_clstr_seed+"");
  			properties.setProperty(prefix+"min_clstr_lone",   this.min_clstr_lone+"");
			properties.setProperty(prefix+"min_clstr_weight", this.min_clstr_weight +"");
			properties.setProperty(prefix+"min_clstr_max",    this.min_clstr_max +"");
  			
  			properties.setProperty(prefix+"fill_gaps",        this.fill_gaps+"");
  			properties.setProperty(prefix+"fill_final",       this.fill_final+"");
  			properties.setProperty(prefix+"min_clstr_block",  this.min_clstr_block+"");
  			properties.setProperty(prefix+"bgnd_grow",        this.bgnd_grow+"");

			properties.setProperty(prefix+"ortho_old",        this.ortho_old+"");
  			properties.setProperty(prefix+"ortho_min_hor",    this.ortho_min_hor +"");
  			properties.setProperty(prefix+"ortho_min_vert",   this.ortho_min_vert +"");
			properties.setProperty(prefix+"ortho_asym",       this.ortho_asym +"");
			properties.setProperty(prefix+"ortho_over4",      this.ortho_over4 +"");
			properties.setProperty(prefix+"ortho_sustain",    this.ortho_sustain +"");
  			properties.setProperty(prefix+"ortho_run",        this.ortho_run+"");
			properties.setProperty(prefix+"ortho_minmax",     this.ortho_minmax +"");
  			properties.setProperty(prefix+"ortho_bridge",     this.ortho_bridge+"");
			properties.setProperty(prefix+"ortho_rms",        this.ortho_rms +"");
  			properties.setProperty(prefix+"ortho_half_length",this.ortho_half_length+"");
			properties.setProperty(prefix+"ortho_mix",        this.ortho_mix +"");
			

			properties.setProperty(prefix+"or_hor",           this.or_hor+"");
			properties.setProperty(prefix+"or_vert",          this.or_vert+"");
			properties.setProperty(prefix+"or_sigma",         this.or_sigma +"");
			properties.setProperty(prefix+"or_sharp",         this.or_sharp +"");
			properties.setProperty(prefix+"or_scale",         this.or_scale +"");
			properties.setProperty(prefix+"or_offset",        this.or_offset +"");
			properties.setProperty(prefix+"or_asym",          this.or_asym +"");
			properties.setProperty(prefix+"or_threshold",     this.or_threshold +"");
			properties.setProperty(prefix+"or_absHor",        this.or_absHor +"");
			properties.setProperty(prefix+"or_absVert",       this.or_absVert +"");
			
			properties.setProperty(prefix+"poles_fix",        this.poles_fix+"");
  			properties.setProperty(prefix+"poles_len",        this.poles_len+"");
			properties.setProperty(prefix+"poles_ratio",      this.poles_ratio +"");
			properties.setProperty(prefix+"poles_min_strength",this.poles_min_strength +"");
			properties.setProperty(prefix+"poles_force_disp", this.poles_force_disp+"");
			
			
			
			properties.setProperty(prefix+"max_clusters",     this.max_clusters+"");
			properties.setProperty(prefix+"remove_scans",     this.remove_scans+"");
			properties.setProperty(prefix+"output_x3d",       this.output_x3d+"");
			properties.setProperty(prefix+"output_obj",       this.output_obj+"");
			properties.setProperty(prefix+"correct_distortions",this.correct_distortions+"");
			properties.setProperty(prefix+"show_triangles",   this.show_triangles+"");
			properties.setProperty(prefix+"avg_cluster_disp", this.avg_cluster_disp+"");
			properties.setProperty(prefix+"maxDispTriangle",  this.maxDispTriangle +"");
			properties.setProperty(prefix+"infinityDistance", this.infinityDistance +"");
			properties.setProperty(prefix+"shUseFlaps",       this.shUseFlaps+"");
			properties.setProperty(prefix+"shAggrFade",       this.shAggrFade+"");
  			properties.setProperty(prefix+"shMinArea",        this.shMinArea+"");
			properties.setProperty(prefix+"shMinStrength",   this.shMinStrength +"");
			properties.setProperty(prefix+"tiRigidVertical",  this.tiRigidVertical +"");
			properties.setProperty(prefix+"tiRigidHorizontal",this.tiRigidHorizontal +"");
			properties.setProperty(prefix+"tiRigidDiagonal",  this.tiRigidDiagonal +"");
			properties.setProperty(prefix+"tiStrengthOffset", this.tiStrengthOffset +"");
			properties.setProperty(prefix+"tiDispScale",      this.tiDispScale +"");
			properties.setProperty(prefix+"tiDispPow",        this.tiDispPow +"");
			properties.setProperty(prefix+"tiDispPull",       this.tiDispPull +"");
			properties.setProperty(prefix+"tiDispPullPreFinal",this.tiDispPullPreFinal +"");
			properties.setProperty(prefix+"tiDispPullFinal",  this.tiDispPullFinal +"");
			properties.setProperty(prefix+"tiBreakNorm",      this.tiBreakNorm +"");
			properties.setProperty(prefix+"tiBreak3",         this.tiBreak3 +"");
			properties.setProperty(prefix+"tiBreak31",        this.tiBreak31 +"");
			properties.setProperty(prefix+"tiBreak21",        this.tiBreak21 +"");
			properties.setProperty(prefix+"tiBreakFar",       this.tiBreakFar +"");
			properties.setProperty(prefix+"tiBreakNear",      this.tiBreakNear +"");
			properties.setProperty(prefix+"tiBreakMode",      this.tiBreakMode +"");
			properties.setProperty(prefix+"tiBreakSame",      this.tiBreakSame +"");
			properties.setProperty(prefix+"tiBreakTurn",      this.tiBreakTurn +"");

			properties.setProperty(prefix+"tiHealPreLast",    this.tiHealPreLast +"");
			properties.setProperty(prefix+"tiHealLast",       this.tiHealLast +"");
			properties.setProperty(prefix+"tiHealSame",       this.tiHealSame+"");

			properties.setProperty(prefix+"tiIterations",     this.tiIterations+"");
  			properties.setProperty(prefix+"tiPrecision",      this.tiPrecision+"");
  			properties.setProperty(prefix+"tiNumCycles",      this.tiNumCycles+"");
			
			properties.setProperty(prefix+"stUseRefine",      this.stUseRefine+"");
			properties.setProperty(prefix+"stUsePass2",       this.stUsePass2+"");
			properties.setProperty(prefix+"stUseRender",      this.stUseRender+"");

			properties.setProperty(prefix+"stShow",           this.stShow+"");
  			properties.setProperty(prefix+"stSize",           this.stSize+"");
			properties.setProperty(prefix+"stStepFar",        this.stStepFar +"");
			properties.setProperty(prefix+"stStepNear",       this.stStepNear +"");
			properties.setProperty(prefix+"stStepThreshold",  this.stStepThreshold +"");
			properties.setProperty(prefix+"stMinDisparity",   this.stMinDisparity +"");
//			properties.setProperty(prefix+"stMaxDisparity",   this.stMaxDisparity +"");
			properties.setProperty(prefix+"stFloor",          this.stFloor +"");
			properties.setProperty(prefix+"stPow",            this.stPow +"");
			properties.setProperty(prefix+"stSigma",          this.stSigma +"");
			properties.setProperty(prefix+"stMinBgDisparity", this.stMinBgDisparity +"");
			properties.setProperty(prefix+"stMinBgFract",     this.stMinBgFract +"");
			properties.setProperty(prefix+"stUseDisp",        this.stUseDisp +"");
			properties.setProperty(prefix+"stStrengthScale",  this.stStrengthScale +"");

			properties.setProperty(prefix+"stSmplMode",       this.stSmplMode+"");
  			properties.setProperty(prefix+"stSmplSide",       this.stSmplSide+"");
  			properties.setProperty(prefix+"stSmplNum",        this.stSmplNum+"");
			properties.setProperty(prefix+"stSmplRms",        this.stSmplRms +"");
			properties.setProperty(prefix+"stSmplWnd",        this.stSmplWnd+"");
 
			properties.setProperty(prefix+"stGrowSel",        this.stGrowSel+"");
			properties.setProperty(prefix+"stMeasSel",        this.stMeasSel+"");
			properties.setProperty(prefix+"stSmallDiff",      this.stSmallDiff +"");
			properties.setProperty(prefix+"stHighMix",        this.stHighMix +"");
  			
			
			properties.setProperty(prefix+"outlayerStrength", this.outlayerStrength +"");
			properties.setProperty(prefix+"outlayerDiff",     this.outlayerDiff +"");
			properties.setProperty(prefix+"outlayerDiffPos",  this.outlayerDiffPos +"");
			properties.setProperty(prefix+"outlayerDiffNeg",  this.outlayerDiffNeg +"");

			properties.setProperty(prefix+"combine_refine",   this.combine_refine+"");

			properties.setProperty(prefix+"combine_min_strength", this.combine_min_strength +"");
			properties.setProperty(prefix+"combine_min_hor",  this.combine_min_hor +"");
			properties.setProperty(prefix+"combine_min_vert", this.combine_min_vert +"");
//			properties.setProperty(prefix+"unique_tolerance", this.unique_tolerance +"");
  			properties.setProperty(prefix+"grow_sweep",       this.grow_sweep+"");
			properties.setProperty(prefix+"grow_disp_max",    this.grow_disp_max +"");
			properties.setProperty(prefix+"grow_disp_trust",  this.grow_disp_trust +"");
			properties.setProperty(prefix+"grow_disp_step",   this.grow_disp_step +"");
			properties.setProperty(prefix+"grow_min_diff",    this.grow_min_diff +"");
			properties.setProperty(prefix+"grow_retry_far",   this.grow_retry_far+"");
			properties.setProperty(prefix+"grow_pedantic",    this.grow_pedantic+"");
			properties.setProperty(prefix+"grow_retry_inf",   this.grow_retry_inf+"");
			
			properties.setProperty(prefix+"gr_new_expand",    this.gr_new_expand+"");
  			properties.setProperty(prefix+"gr_max_expand",    this.gr_max_expand+"");
			properties.setProperty(prefix+"fds_str_floor",    this.fds_str_floor +"");
			properties.setProperty(prefix+"gr_ovrbg_cmb",     this.gr_ovrbg_cmb +"");
			properties.setProperty(prefix+"gr_ovrbg_cmb_hor", this.gr_ovrbg_cmb_hor +"");
			properties.setProperty(prefix+"gr_ovrbg_cmb_vert",this.gr_ovrbg_cmb_vert +"");
			properties.setProperty(prefix+"gr_ovrbg_filtered",this.gr_ovrbg_filtered +"");
			properties.setProperty(prefix+"fds_str_pow",      this.fds_str_pow +"");
  			properties.setProperty(prefix+"fds_smpl_side",    this.fds_smpl_side+"");
  			properties.setProperty(prefix+"fds_smpl_num",     this.fds_smpl_num+"");
			properties.setProperty(prefix+"fds_smpl_rms",     this.fds_smpl_rms +"");
			properties.setProperty(prefix+"fds_smpl_rel_rms", this.fds_smpl_rel_rms +"");
			properties.setProperty(prefix+"fds_smpl_wnd",     this.fds_smpl_wnd+"");
			properties.setProperty(prefix+"fds_abs_tilt",     this.fds_abs_tilt +"");
			properties.setProperty(prefix+"fds_rel_tilt",     this.fds_rel_tilt +"");
			
			properties.setProperty(prefix+"mc_disp8_step",    this.mc_disp8_step +"");
			properties.setProperty(prefix+"mc_disp8_trust",   this.mc_disp8_trust +"");
			properties.setProperty(prefix+"mc_strength",      this.mc_strength +"");
			properties.setProperty(prefix+"mc_unique_tol",    this.mc_unique_tol +"");
			
			properties.setProperty(prefix+"mc_trust_fin",     this.mc_trust_fin +"");
			properties.setProperty(prefix+"mc_trust_sigma",   this.mc_trust_sigma +"");
			properties.setProperty(prefix+"mc_ortho_weight",  this.mc_ortho_weight +"");
			properties.setProperty(prefix+"mc_diag_weight",   this.mc_diag_weight +"");
			properties.setProperty(prefix+"mc_gap",           this.mc_gap +"");

			properties.setProperty(prefix+"gr_min_new",       this.gr_min_new+"");
			properties.setProperty(prefix+"gr_var_new_sngl",  this.gr_var_new_sngl+"");
			properties.setProperty(prefix+"gr_var_new_fg",    this.gr_var_new_fg+"");
			properties.setProperty(prefix+"gr_var_all_fg",    this.gr_var_all_fg+"");
			properties.setProperty(prefix+"gr_var_new_bg",    this.gr_var_new_bg+"");
			properties.setProperty(prefix+"gr_var_all_bg",    this.gr_var_all_bg+"");
			properties.setProperty(prefix+"gr_var_next",      this.gr_var_next+"");
  			properties.setProperty(prefix+"gr_num_steps",     this.gr_num_steps+"");
  			properties.setProperty(prefix+"gr_steps_over",    this.gr_steps_over+"");
  			properties.setProperty(prefix+"gr_smpl_size",     this.gr_smpl_size+"");
  			properties.setProperty(prefix+"gr_min_pnts",      this.gr_min_pnts+"");
			properties.setProperty(prefix+"gr_use_wnd",       this.gr_use_wnd+"");
			properties.setProperty(prefix+"gr_tilt_damp",     this.gr_tilt_damp +"");
			properties.setProperty(prefix+"gr_split_rng",     this.gr_split_rng +"");
			properties.setProperty(prefix+"gr_same_rng",      this.gr_same_rng +"");
			properties.setProperty(prefix+"gr_diff_cont",     this.gr_diff_cont +"");
			properties.setProperty(prefix+"gr_abs_tilt",      this.gr_abs_tilt +"");
			properties.setProperty(prefix+"gr_rel_tilt",      this.gr_rel_tilt +"");
  			properties.setProperty(prefix+"gr_smooth",        this.gr_smooth+"");
			properties.setProperty(prefix+"gr_fin_diff",      this.gr_fin_diff +"");
			properties.setProperty(prefix+"gr_unique_tol",    this.gr_unique_tol +"");
			properties.setProperty(prefix+"gr_unique_pretol", this.gr_unique_pretol +"");

			properties.setProperty(prefix+"plPreferDisparity",this.plPreferDisparity+"");
			properties.setProperty(prefix+"plDispNorm",       this.plDispNorm +"");

			properties.setProperty(prefix+"plBlurBinVert",    this.plBlurBinVert +"");
			properties.setProperty(prefix+"plBlurBinHor",     this.plBlurBinHor +"");
			properties.setProperty(prefix+"plMaxDiffVert",    this.plMaxDiffVert +"");
			properties.setProperty(prefix+"plMaxDiffHor",     this.plMaxDiffHor +"");
  			properties.setProperty(prefix+"plInitPasses",     this.plInitPasses+"");
			
			properties.setProperty(prefix+"plMinPoints",      this.plMinPoints+"");
			properties.setProperty(prefix+"plTargetEigen",    this.plTargetEigen +"");
			properties.setProperty(prefix+"plFractOutliers",  this.plFractOutliers +"");
  			properties.setProperty(prefix+"plMaxOutliers",    this.plMaxOutliers+"");
			properties.setProperty(prefix+"plMinStrength",    this.plMinStrength +"");
			properties.setProperty(prefix+"plMaxEigen",       this.plMaxEigen +"");
			properties.setProperty(prefix+"plEigenFloor",     this.plEigenFloor +"");
			properties.setProperty(prefix+"plEigenStick",     this.plEigenStick +"");
			properties.setProperty(prefix+"plBadPlate",       this.plBadPlate +"");
			properties.setProperty(prefix+"plDbgMerge",       this.plDbgMerge+"");
			properties.setProperty(prefix+"plWorstWorsening", this.plWorstWorsening +"");
			properties.setProperty(prefix+"plWorstWorsening2",this.plWorstWorsening2 +"");
			properties.setProperty(prefix+"plWorstEq",        this.plWorstEq +"");
			properties.setProperty(prefix+"plWorstEq2",       this.plWorstEq2 +"");
			properties.setProperty(prefix+"plOKMergeEigen",   this.plOKMergeEigen +"");
			properties.setProperty(prefix+"plMaxWorldSin2",   this.plMaxWorldSin2 +"");
			properties.setProperty(prefix+"pl2dForSin",       this.pl2dForSin +"");
			properties.setProperty(prefix+"plWeakWorsening",  this.plWeakWorsening +"");
			properties.setProperty(prefix+"plMaxOverlap",     this.plMaxOverlap +"");

			properties.setProperty(prefix+"plWeakWeight",     this.plWeakWeight +"");
			properties.setProperty(prefix+"plWeakEigen",      this.plWeakEigen +"");
			properties.setProperty(prefix+"plWeakWeight2",    this.plWeakWeight2 +"");
			properties.setProperty(prefix+"plWeakEigen2",     this.plWeakEigen2 +"");
			properties.setProperty(prefix+"plSumThick",       this.plSumThick +"");
			properties.setProperty(prefix+"plNeNeibCost",     this.plNeNeibCost +"");
			properties.setProperty(prefix+"plNeOwn",          this.plNeOwn +"");

			properties.setProperty(prefix+"plExNeibCost",     this.plExNeibCost +"");
			properties.setProperty(prefix+"plExNeibSmooth",   this.plExNeibSmooth +"");
			properties.setProperty(prefix+"plMergeCostStar",  this.plMergeCostStar +"");
			properties.setProperty(prefix+"plMergeCost",      this.plMergeCost +"");
			
			properties.setProperty(prefix+"plConflMerge",     this.plConflMerge+"");
			properties.setProperty(prefix+"plConflRelax",     this.plConflRelax +"");
			properties.setProperty(prefix+"plConflSngl",      this.plConflSngl+"");
			properties.setProperty(prefix+"plConflSnglPair",  this.plConflSnglPair+"");

			properties.setProperty(prefix+"plWeakFgStrength", this.plWeakFgStrength +"");
  			properties.setProperty(prefix+"plWeakFgOutliers", this.plWeakFgOutliers+"");
			properties.setProperty(prefix+"plWeakFgRelax",    this.plWeakFgRelax +"");

			properties.setProperty(prefix+"plThickWorld",     this.plThickWorld +"");
			properties.setProperty(prefix+"plThickWorldConfl",this.plThickWorldConfl +"");
			properties.setProperty(prefix+"plRelaxComplete",  this.plRelaxComplete +"");
			properties.setProperty(prefix+"plRelaxComplete2", this.plRelaxComplete2 +"");

			properties.setProperty(prefix+"plMaxZRatio",      this.plMaxZRatio +"");
			properties.setProperty(prefix+"plMaxDisp",        this.plMaxDisp +"");
			properties.setProperty(prefix+"plCutTail",        this.plCutTail +"");
			properties.setProperty(prefix+"plMinTail",        this.plMinTail +"");

			properties.setProperty(prefix+"plDiscrEn",        this.plDiscrEn+"");
			properties.setProperty(prefix+"plDiscrTolerance", this.plDiscrTolerance +"");
			properties.setProperty(prefix+"plDiscrDispRange", this.plDiscrDispRange +"");
  			properties.setProperty(prefix+"plDiscrSteps",     this.plDiscrSteps+"");
  			properties.setProperty(prefix+"plDiscrVariants",  this.plDiscrVariants+"");
  			properties.setProperty(prefix+"plDiscrMode",      this.plDiscrMode+"");
			properties.setProperty(prefix+"plDiscrVarFloor",  this.plDiscrVarFloor +"");
			properties.setProperty(prefix+"plDiscrSigma",     this.plDiscrSigma +"");
			properties.setProperty(prefix+"plDiscrBlur",      this.plDiscrBlur +"");
			properties.setProperty(prefix+"plDiscrExclusivity",this.plDiscrExclusivity +"");
			properties.setProperty(prefix+"plDiscrExclus2",   this.plDiscrExclus2 +"");
			properties.setProperty(prefix+"plDiscrStrict",    this.plDiscrStrict+"");
			properties.setProperty(prefix+"plDiscrCorrMax",   this.plDiscrCorrMax +"");
			properties.setProperty(prefix+"plDiscrCorrMerge", this.plDiscrCorrMerge +"");
  			properties.setProperty(prefix+"plDiscrSteal",     this.plDiscrSteal+"");
  			properties.setProperty(prefix+"plDiscrGrown",     this.plDiscrGrown+"");
			properties.setProperty(prefix+"plDiscrXMedian",   this.plDiscrXMedian +"");

			properties.setProperty(prefix+"plCostDist",       this.plCostDist +"");
			properties.setProperty(prefix+"plCostKrq",        this.plCostKrq +"");
			properties.setProperty(prefix+"plCostKrqEq",      this.plCostKrqEq +"");
			properties.setProperty(prefix+"plCostWrq",        this.plCostWrq +"");
			properties.setProperty(prefix+"plCostWrqEq",      this.plCostWrqEq +"");
			properties.setProperty(prefix+"plCostSin2",       this.plCostSin2 +"");
			properties.setProperty(prefix+"plCostRdist2",     this.plCostRdist2 +"");

			properties.setProperty(prefix+"plConflDualTri",   this.plConflDualTri+"");
			properties.setProperty(prefix+"plConflMulti",     this.plConflMulti+"");
			properties.setProperty(prefix+"plConflDiag",      this.plConflDiag+"");
			properties.setProperty(prefix+"plConflStar",      this.plConflStar+"");
  			properties.setProperty(prefix+"plStarSteps",      this.plStarSteps+"");
			properties.setProperty(prefix+"plStarOrtho",      this.plStarOrtho +"");
			properties.setProperty(prefix+"plStarDiag",       this.plStarDiag +"");
			properties.setProperty(prefix+"plStarPwr",        this.plStarPwr +"");
			properties.setProperty(prefix+"plStarWeightPwr",  this.plStarWeightPwr +"");
			properties.setProperty(prefix+"plWeightToDens",   this.plWeightToDens +"");
			properties.setProperty(prefix+"plStarValPwr",     this.plStarValPwr +"");
			properties.setProperty(prefix+"plDblTriLoss",     this.plDblTriLoss +"");
			properties.setProperty(prefix+"plNewConfl",       this.plNewConfl+"");
  			properties.setProperty(prefix+"plMaxChanges",     this.plMaxChanges+"");

			properties.setProperty(prefix+"plMutualOnly",     this.plMutualOnly+"");
			properties.setProperty(prefix+"plFillSquares",    this.plFillSquares+"");
			properties.setProperty(prefix+"plCutCorners",     this.plCutCorners+"");
			properties.setProperty(prefix+"plHypotenuse",     this.plHypotenuse+"");
			properties.setProperty(prefix+"plPull",           this.plPull +"");
			properties.setProperty(prefix+"plNormPow",        this.plNormPow +"");
  			properties.setProperty(prefix+"plIterations",     this.plIterations+"");
			properties.setProperty(prefix+"plStopBad",        this.plStopBad+"");
  			properties.setProperty(prefix+"plPrecision",      this.plPrecision+"");

			properties.setProperty(prefix+"plSplitPull",      this.plSplitPull +"");
  			properties.setProperty(prefix+"plSplitMinNeib",   this.plSplitMinNeib+"");
			properties.setProperty(prefix+"plSplitMinWeight", this.plSplitMinWeight +"");
			properties.setProperty(prefix+"plSplitMinQuality",this.plSplitMinQuality +"");
			properties.setProperty(prefix+"plSplitApply",     this.plSplitApply+"");
			properties.setProperty(prefix+"plNonExclusive",   this.plNonExclusive+"");
			properties.setProperty(prefix+"plUseOtherPlanes", this.plUseOtherPlanes+"");
			properties.setProperty(prefix+"plAllowParallel",  this.plAllowParallel+"");
			properties.setProperty(prefix+"plMaxDiff",        this.plMaxDiff +"");
			properties.setProperty(prefix+"plOtherDiff",      this.plOtherDiff +"");
			properties.setProperty(prefix+"plSplitXY",        this.plSplitXY+"");
			properties.setProperty(prefix+"plSplitXYTolerance",this.plSplitXYTolerance +"");

  			properties.setProperty(prefix+"plFuse",           this.plFuse+"");
			properties.setProperty(prefix+"plKeepOrphans",    this.plKeepOrphans+"");
			properties.setProperty(prefix+"plMinOrphan",      this.plMinOrphan +"");
			
			properties.setProperty(prefix+"plSnapDispAny",    this.plSnapDispAny +"");
			properties.setProperty(prefix+"plSnapStrengthAny",this.plSnapStrengthAny +"");
			properties.setProperty(prefix+"plSnapNegAny",     this.plSnapNegAny +"");
			properties.setProperty(prefix+"plSnapDispMax",    this.plSnapDispMax +"");
			properties.setProperty(prefix+"plSnapDispWeight", this.plSnapDispWeight +"");
  			properties.setProperty(prefix+"plSnapZeroMode",   this.plSnapZeroMode+"");

			properties.setProperty(prefix+"msUseSel",         this.msUseSel+"");
			properties.setProperty(prefix+"msDivideByArea",   this.msDivideByArea+"");
			properties.setProperty(prefix+"msScaleProj",      this.msScaleProj +"");
			properties.setProperty(prefix+"msFractUni",       this.msFractUni +"");

			properties.setProperty(prefix+"tsNoEdge",         this.tsNoEdge+"");
			properties.setProperty(prefix+"tsUseCenter",      this.tsUseCenter+"");
			properties.setProperty(prefix+"tsMaxDiff",        this.tsMaxDiff +"");
			properties.setProperty(prefix+"tsMinDiffOther",   this.tsMinDiffOther +"");
			properties.setProperty(prefix+"tsMinStrength",    this.tsMinStrength +"");
			properties.setProperty(prefix+"tsMaxStrength",    this.tsMaxStrength +"");
			properties.setProperty(prefix+"tsMinSurface",     this.tsMinSurface +"");
  			properties.setProperty(prefix+"tsMoveDirs",       this.tsMoveDirs+"");
			properties.setProperty(prefix+"tsSurfStrPow",     this.tsSurfStrPow +"");
			properties.setProperty(prefix+"tsAddStrength",    this.tsAddStrength +"");
			properties.setProperty(prefix+"tsSigma",          this.tsSigma +"");
			properties.setProperty(prefix+"tsNSigma",         this.tsNSigma +"");
			properties.setProperty(prefix+"tsMinPull",        this.tsMinPull +"");
			properties.setProperty(prefix+"tsMinAdvantage",   this.tsMinAdvantage +"");
			
			properties.setProperty(prefix+"tsClustSize",      this.tsClustSize +"");
			properties.setProperty(prefix+"tsClustWeight",    this.tsClustWeight +"");
			properties.setProperty(prefix+"tsMinNeib",        this.tsMinNeib +"");
			properties.setProperty(prefix+"tsMaxSurStrength", this.tsMaxSurStrength +"");
			properties.setProperty(prefix+"tsCountDis",       this.tsCountDis +"");
			
			properties.setProperty(prefix+"tsEnPlaneSeed",    this.tsEnPlaneSeed+"");
			properties.setProperty(prefix+"tsEnOnly",         this.tsEnOnly+"");
			properties.setProperty(prefix+"tsEnGrow",         this.tsEnGrow+"");
			properties.setProperty(prefix+"tsGrowStrength",   this.tsGrowStrength +"");
			properties.setProperty(prefix+"tsGrowStrong",     this.tsGrowStrong+"");
			properties.setProperty(prefix+"tsContStrength",   this.tsContStrength +"");
			properties.setProperty(prefix+"tsContDiff",       this.tsContDiff +"");

			properties.setProperty(prefix+"tsEnSingle",       this.tsEnSingle+"");
			properties.setProperty(prefix+"tsEnMulti",        this.tsEnMulti+"");
			properties.setProperty(prefix+"tsRemoveWeak1",    this.tsRemoveWeak1+"");
			properties.setProperty(prefix+"tsGrowSurround",   this.tsGrowSurround+"");
			properties.setProperty(prefix+"tsRemoveWeak2",    this.tsRemoveWeak2+"");
			properties.setProperty(prefix+"tsLoopMulti",      this.tsLoopMulti+"");
			properties.setProperty(prefix+"tsShow",           this.tsShow+"");
			properties.setProperty(prefix+"tsNumClust",       this.tsNumClust +"");
			
			properties.setProperty(prefix+"tsConsensMode",    this.tsConsensMode +"");
			properties.setProperty(prefix+"tsConsensAgree",   this.tsConsensAgree +"");

			properties.setProperty(prefix+"taMinFgBg",        this.taMinFgBg +"");
			properties.setProperty(prefix+"taMinFgEdge",      this.taMinFgEdge +"");
			properties.setProperty(prefix+"taMinColSep",      this.taMinColSep +"");
			properties.setProperty(prefix+"taMinColDiff",     this.taMinColDiff +"");
			properties.setProperty(prefix+"taOutlier",        this.taOutlier +"");
			properties.setProperty(prefix+"taDiffPwr",        this.taDiffPwr +"");
			properties.setProperty(prefix+"taBestPwr",        this.taBestPwr +"");
			properties.setProperty(prefix+"taDiff9Pwr",       this.taDiff9Pwr +"");
			properties.setProperty(prefix+"taColSigma",       this.taColSigma +"");
			properties.setProperty(prefix+"taColFraction",    this.taColFraction +"");
			
			properties.setProperty(prefix+"taCostEmpty",      this.taCostEmpty +"");
			properties.setProperty(prefix+"taCostNoLink",     this.taCostNoLink +"");
			properties.setProperty(prefix+"taCostSwitch",     this.taCostSwitch +"");
			properties.setProperty(prefix+"taCostColor",      this.taCostColor +"");
			properties.setProperty(prefix+"taCostDiff",       this.taCostDiff +"");
			properties.setProperty(prefix+"taCostDiffBest",   this.taCostDiffBest +"");
			properties.setProperty(prefix+"taCostDiff9",      this.taCostDiff9 +"");
			properties.setProperty(prefix+"taCostWeakFgnd",   this.taCostWeakFgnd +"");
			properties.setProperty(prefix+"taCostFlaps",      this.taCostFlaps +"");
			properties.setProperty(prefix+"taCostMismatch",   this.taCostMismatch +"");

			properties.setProperty(prefix+"taEnEmpty",        this.taEnEmpty +"");
			properties.setProperty(prefix+"taEnNoLink",       this.taEnNoLink +"");
			properties.setProperty(prefix+"taEnSwitch",       this.taEnSwitch +"");
			properties.setProperty(prefix+"taEnColor",        this.taEnColor +"");
			properties.setProperty(prefix+"taEnDiff",         this.taEnDiff +"");
			properties.setProperty(prefix+"taEnDiffBest",     this.taEnDiffBest +"");
			properties.setProperty(prefix+"taEnDiff9",        this.taEnDiff9 +"");
			properties.setProperty(prefix+"taEnWeakFgnd",     this.taEnWeakFgnd +"");
			properties.setProperty(prefix+"taEnFlaps",        this.taEnFlaps +"");
			properties.setProperty(prefix+"taEnMismatch",     this.taEnMismatch +"");
			
			
			properties.setProperty(prefix+"dbg_migrate",            this.dbg_migrate+"");
  			
			properties.setProperty(prefix+"show_ortho_combine",     this.show_ortho_combine+"");
			properties.setProperty(prefix+"show_refine_supertiles", this.show_refine_supertiles+"");
			properties.setProperty(prefix+"show_bgnd_nonbgnd",      this.show_bgnd_nonbgnd+"");
			properties.setProperty(prefix+"show_filter_scan",       this.show_filter_scan+"");
			properties.setProperty(prefix+"show_combined",          this.show_combined+"");
			properties.setProperty(prefix+"show_unique",            this.show_unique+"");
			properties.setProperty(prefix+"show_histograms",        this.show_histograms+"");
			properties.setProperty(prefix+"show_init_refine",       this.show_init_refine+"");
			properties.setProperty(prefix+"show_expand",            this.show_expand+"");
			properties.setProperty(prefix+"show_variant",           this.show_variant+"");
			properties.setProperty(prefix+"show_retry_far",         this.show_retry_far+"");
			properties.setProperty(prefix+"show_macro",             this.show_macro+"");
			properties.setProperty(prefix+"show_shells",            this.show_shells+"");
			properties.setProperty(prefix+"show_neighbors",         this.show_neighbors+"");
			properties.setProperty(prefix+"show_flaps_dirs",        this.show_flaps_dirs+"");
			properties.setProperty(prefix+"show_first_clusters",    this.show_first_clusters+"");
			properties.setProperty(prefix+"show_planes",            this.show_planes+"");

			properties.setProperty(prefix+"vertical_xyz.x",         this.vertical_xyz[0]+"");
			properties.setProperty(prefix+"vertical_xyz.y",         this.vertical_xyz[1]+"");
			properties.setProperty(prefix+"vertical_xyz.z",         this.vertical_xyz[2]+"");
}
  		public void getProperties(String prefix,Properties properties){
  			if (properties.getProperty(prefix+"transform_size")!=null) this.transform_size=Integer.parseInt(properties.getProperty(prefix+"transform_size"));
  			if (properties.getProperty(prefix+"clt_window")!=null)     this.clt_window=Integer.parseInt(properties.getProperty(prefix+"clt_window"));
  			if (properties.getProperty(prefix+"shift_x")!=null)        this.shift_x=Double.parseDouble(properties.getProperty(prefix+"shift_x"));
  			if (properties.getProperty(prefix+"shift_y")!=null)        this.shift_y=Double.parseDouble(properties.getProperty(prefix+"shift_y"));
  			if (properties.getProperty(prefix+"iclt_mask")!=null)      this.iclt_mask=Integer.parseInt(properties.getProperty(prefix+"iclt_mask"));
  			if (properties.getProperty(prefix+"tileX")!=null)          this.tileX=Integer.parseInt(properties.getProperty(prefix+"tileX"));
  			if (properties.getProperty(prefix+"tileY")!=null)          this.tileY=Integer.parseInt(properties.getProperty(prefix+"tileY"));
  			if (properties.getProperty(prefix+"dbg_mode")!=null)       this.dbg_mode=Integer.parseInt(properties.getProperty(prefix+"dbg_mode"));
  			if (properties.getProperty(prefix+"ishift_x")!=null)       this.ishift_x=Integer.parseInt(properties.getProperty(prefix+"ishift_x"));
  			if (properties.getProperty(prefix+"ishift_y")!=null)       this.ishift_y=Integer.parseInt(properties.getProperty(prefix+"ishift_y"));
  			if (properties.getProperty(prefix+"fat_zero")!=null)       this.fat_zero=Double.parseDouble(properties.getProperty(prefix+"fat_zero"));
  			if (properties.getProperty(prefix+"corr_sigma")!=null)     this.corr_sigma=Double.parseDouble(properties.getProperty(prefix+"corr_sigma"));
  			if (properties.getProperty(prefix+"norm_kern")!=null)      this.norm_kern=Boolean.parseBoolean(properties.getProperty(prefix+"norm_kern"));
  			if (properties.getProperty(prefix+"gain_equalize")!=null)  this.gain_equalize=Boolean.parseBoolean(properties.getProperty(prefix+"gain_equalize"));
  			if (properties.getProperty(prefix+"colors_equalize")!=null)this.colors_equalize=Boolean.parseBoolean(properties.getProperty(prefix+"colors_equalize"));
  			if (properties.getProperty(prefix+"novignetting_r")!=null) this.novignetting_r=Double.parseDouble(properties.getProperty(prefix+"novignetting_r"));
  			if (properties.getProperty(prefix+"novignetting_g")!=null) this.novignetting_g=Double.parseDouble(properties.getProperty(prefix+"novignetting_g"));
  			if (properties.getProperty(prefix+"novignetting_b")!=null) this.novignetting_b=Double.parseDouble(properties.getProperty(prefix+"novignetting_b"));
  			if (properties.getProperty(prefix+"scale_r")!=null)        this.scale_r=Double.parseDouble(properties.getProperty(prefix+"scale_r"));
  			if (properties.getProperty(prefix+"scale_g")!=null)        this.scale_g=Double.parseDouble(properties.getProperty(prefix+"scale_g"));
  			if (properties.getProperty(prefix+"scale_b")!=null)        this.scale_b=Double.parseDouble(properties.getProperty(prefix+"scale_b"));
  			if (properties.getProperty(prefix+"vignetting_max")!=null) this.vignetting_max=Double.parseDouble(properties.getProperty(prefix+"vignetting_max"));
  			if (properties.getProperty(prefix+"vignetting_range")!=null) this.vignetting_range=Double.parseDouble(properties.getProperty(prefix+"vignetting_range"));
  			if (properties.getProperty(prefix+"kernel_step")!=null)    this.kernel_step=Integer.parseInt(properties.getProperty(prefix+"kernel_step"));
  			if (properties.getProperty(prefix+"disparity")!=null)      this.disparity=Double.parseDouble(properties.getProperty(prefix+"disparity"));
  			if (properties.getProperty(prefix+"z_correction")!=null)   this.z_correction=Double.parseDouble(properties.getProperty(prefix+"z_correction"));
  			if (properties.getProperty(prefix+"correlate")!=null)      this.correlate=Boolean.parseBoolean(properties.getProperty(prefix+"correlate"));
  			if (properties.getProperty(prefix+"corr_mask")!=null)      this.corr_mask=Integer.parseInt(properties.getProperty(prefix+"corr_mask"));
  			if (properties.getProperty(prefix+"corr_sym")!=null)       this.corr_sym=Boolean.parseBoolean(properties.getProperty(prefix+"corr_sym"));
  			if (properties.getProperty(prefix+"corr_keep")!=null)      this.corr_keep=Boolean.parseBoolean(properties.getProperty(prefix+"corr_keep"));
  			if (properties.getProperty(prefix+"corr_show")!=null)      this.corr_show=Boolean.parseBoolean(properties.getProperty(prefix+"corr_show"));
  			if (properties.getProperty(prefix+"corr_mismatch")!=null)  this.corr_mismatch=Boolean.parseBoolean(properties.getProperty(prefix+"corr_mismatch"));
  			if (properties.getProperty(prefix+"corr_offset")!=null)    this.corr_offset=Double.parseDouble(properties.getProperty(prefix+"corr_offset"));
  			if (properties.getProperty(prefix+"corr_red")!=null)       this.corr_red=Double.parseDouble(properties.getProperty(prefix+"corr_red"));
  			if (properties.getProperty(prefix+"corr_blue")!=null)      this.corr_blue=Double.parseDouble(properties.getProperty(prefix+"corr_blue"));
  			if (properties.getProperty(prefix+"corr_normalize")!=null) this.corr_normalize=Boolean.parseBoolean(properties.getProperty(prefix+"corr_normalize"));
  			if (properties.getProperty(prefix+"min_corr")!=null)       this.min_corr=Double.parseDouble(properties.getProperty(prefix+"min_corr"));
  			if (properties.getProperty(prefix+"min_corr_normalized")!=null)this.min_corr_normalized=Double.parseDouble(properties.getProperty(prefix+"min_corr_normalized"));
  			if (properties.getProperty(prefix+"max_corr_sigma")!=null) this.max_corr_sigma=Double.parseDouble(properties.getProperty(prefix+"max_corr_sigma"));
  			if (properties.getProperty(prefix+"max_corr_radius")!=null)this.max_corr_radius=Double.parseDouble(properties.getProperty(prefix+"max_corr_radius"));

  			if (properties.getProperty(prefix+"enhortho_width")!=null) this.enhortho_width=Integer.parseInt(properties.getProperty(prefix+"enhortho_width"));
  			if (properties.getProperty(prefix+"enhortho_scale")!=null) this.enhortho_scale=Double.parseDouble(properties.getProperty(prefix+"enhortho_scale"));
  			
  			if (properties.getProperty(prefix+"max_corr_double")!=null)this.max_corr_double=Boolean.parseBoolean(properties.getProperty(prefix+"max_corr_double"));
  			if (properties.getProperty(prefix+"corr_mode")!=null)      this.corr_mode=Integer.parseInt(properties.getProperty(prefix+"corr_mode"));
  			if (properties.getProperty(prefix+"corr_border_contrast")!=null) this.corr_border_contrast=Double.parseDouble(properties.getProperty(prefix+"corr_border_contrast"));
  			if (properties.getProperty(prefix+"tile_task_op")!=null)   this.tile_task_op=Integer.parseInt(properties.getProperty(prefix+"tile_task_op"));
  			if (properties.getProperty(prefix+"tile_task_wl")!=null)   this.tile_task_wl=Integer.parseInt(properties.getProperty(prefix+"tile_task_wl"));
  			if (properties.getProperty(prefix+"tile_task_wt")!=null)   this.tile_task_wt=Integer.parseInt(properties.getProperty(prefix+"tile_task_wt"));
  			if (properties.getProperty(prefix+"tile_task_ww")!=null)   this.tile_task_ww=Integer.parseInt(properties.getProperty(prefix+"tile_task_ww"));
  			if (properties.getProperty(prefix+"tile_task_wh")!=null)   this.tile_task_wh=Integer.parseInt(properties.getProperty(prefix+"tile_task_wh"));
  			if (properties.getProperty(prefix+"min_shot")!=null)       this.min_shot=Double.parseDouble(properties.getProperty(prefix+"min_shot"));
  			if (properties.getProperty(prefix+"scale_shot")!=null)     this.scale_shot=Double.parseDouble(properties.getProperty(prefix+"scale_shot"));
  			if (properties.getProperty(prefix+"diff_sigma")!=null)     this.diff_sigma=Double.parseDouble(properties.getProperty(prefix+"diff_sigma"));
  			if (properties.getProperty(prefix+"diff_threshold")!=null) this.diff_threshold=Double.parseDouble(properties.getProperty(prefix+"diff_threshold"));
  			if (properties.getProperty(prefix+"diff_gauss")!=null)     this.diff_gauss=Boolean.parseBoolean(properties.getProperty(prefix+"diff_gauss"));
  			if (properties.getProperty(prefix+"min_agree")!=null)      this.min_agree=Double.parseDouble(properties.getProperty(prefix+"min_agree"));
  			if (properties.getProperty(prefix+"dust_remove")!=null)    this.dust_remove=Boolean.parseBoolean(properties.getProperty(prefix+"dust_remove"));
  			if (properties.getProperty(prefix+"black_back")!=null)     this.black_back=Boolean.parseBoolean(properties.getProperty(prefix+"black_back"));
  			if (properties.getProperty(prefix+"keep_weights")!=null)   this.keep_weights=Boolean.parseBoolean(properties.getProperty(prefix+"keep_weights"));
  			if (properties.getProperty(prefix+"sharp_alpha")!=null)    this.sharp_alpha=Boolean.parseBoolean(properties.getProperty(prefix+"sharp_alpha"));
  			if (properties.getProperty(prefix+"alpha0")!=null)         this.alpha0=Double.parseDouble(properties.getProperty(prefix+"alpha0"));
  			if (properties.getProperty(prefix+"alpha1")!=null)         this.alpha1=Double.parseDouble(properties.getProperty(prefix+"alpha1"));
  			if (properties.getProperty(prefix+"gen_chn_stacks")!=null) this.gen_chn_stacks=Boolean.parseBoolean(properties.getProperty(prefix+"gen_chn_stacks"));
  			if (properties.getProperty(prefix+"gen_chn_img")!=null)    this.gen_chn_img=Boolean.parseBoolean(properties.getProperty(prefix+"gen_chn_img"));
  			if (properties.getProperty(prefix+"gen_4_img")!=null)      this.gen_4_img=Boolean.parseBoolean(properties.getProperty(prefix+"gen_4_img"));
  			if (properties.getProperty(prefix+"show_nonoverlap")!=null)this.show_nonoverlap=Boolean.parseBoolean(properties.getProperty(prefix+"show_nonoverlap"));
  			if (properties.getProperty(prefix+"show_overlap")!=null)   this.show_overlap=Boolean.parseBoolean(properties.getProperty(prefix+"show_overlap"));
  			if (properties.getProperty(prefix+"show_rgba_color")!=null)this.show_rgba_color=Boolean.parseBoolean(properties.getProperty(prefix+"show_rgba_color"));
  			if (properties.getProperty(prefix+"show_map")!=null)       this.show_map=Boolean.parseBoolean(properties.getProperty(prefix+"show_map"));
  			if (properties.getProperty(prefix+"disp_scan_start")!=null)this.disp_scan_start=Double.parseDouble(properties.getProperty(prefix+"disp_scan_start"));
  			if (properties.getProperty(prefix+"disp_scan_step")!=null) this.disp_scan_step=Double.parseDouble(properties.getProperty(prefix+"disp_scan_step"));
  			if (properties.getProperty(prefix+"disp_scan_count")!=null)this.disp_scan_count=Integer.parseInt(properties.getProperty(prefix+"disp_scan_count"));
  			
  			if (properties.getProperty(prefix+"fine_dbg")!=null)       this.fine_dbg=Boolean.parseBoolean(properties.getProperty(prefix+"fine_dbg"));
  			if (properties.getProperty(prefix+"fine_corr_x_0")!=null) this.fine_corr_x_0=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_0"));
  			if (properties.getProperty(prefix+"fine_corr_y_0")!=null) this.fine_corr_y_0=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_0"));
  			if (properties.getProperty(prefix+"fine_corr_x_1")!=null) this.fine_corr_x_1=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_1"));
  			if (properties.getProperty(prefix+"fine_corr_y_1")!=null) this.fine_corr_y_1=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_1"));
  			if (properties.getProperty(prefix+"fine_corr_x_2")!=null) this.fine_corr_x_2=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_2"));
  			if (properties.getProperty(prefix+"fine_corr_y_2")!=null) this.fine_corr_y_2=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_2"));
  			if (properties.getProperty(prefix+"fine_corr_x_3")!=null) this.fine_corr_x_3=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_3"));
  			if (properties.getProperty(prefix+"fine_corr_y_3")!=null) this.fine_corr_y_3=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_3"));
  			
  			if (properties.getProperty(prefix+"fcorr_radius")!=null)      this.fcorr_radius=Double.parseDouble(properties.getProperty(prefix+"fcorr_radius"));
  			if (properties.getProperty(prefix+"fcorr_min_strength")!=null) this.fcorr_min_strength=Double.parseDouble(properties.getProperty(prefix+"fcorr_min_strength"));
  			if (properties.getProperty(prefix+"fcorr_disp_diff")!=null)   this.fcorr_disp_diff=Double.parseDouble(properties.getProperty(prefix+"fcorr_disp_diff"));
  			if (properties.getProperty(prefix+"fcorr_quadratic")!=null)   this.fcorr_quadratic=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_quadratic"));
  			if (properties.getProperty(prefix+"fcorr_ignore")!=null)      this.fcorr_ignore=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_ignore"));

  			if (properties.getProperty(prefix+"fcorr_inf_strength")!=null) this.fcorr_inf_strength=Double.parseDouble(properties.getProperty(prefix+"fcorr_inf_strength"));
  			if (properties.getProperty(prefix+"fcorr_inf_diff")!=null)    this.fcorr_inf_diff=Double.parseDouble(properties.getProperty(prefix+"fcorr_inf_diff"));
  			if (properties.getProperty(prefix+"fcorr_inf_quad")!=null)    this.fcorr_inf_quad=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_inf_quad"));
  			if (properties.getProperty(prefix+"fcorr_inf_vert")!=null)    this.fcorr_inf_vert=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_inf_vert"));

  			
  			
  			if (properties.getProperty(prefix+"inf_disp_apply")!=null)    this.inf_disp_apply=Boolean.parseBoolean(properties.getProperty(prefix+"inf_disp_apply"));
  			if (properties.getProperty(prefix+"inf_repeat")!=null)        this.inf_repeat=Integer.parseInt(properties.getProperty(prefix+"inf_repeat"));
//  			if (properties.getProperty(prefix+"inf_mism_apply")!=null)    this.inf_mism_apply=Boolean.parseBoolean(properties.getProperty(prefix+"inf_mism_apply"));

  			if (properties.getProperty(prefix+"inf_iters")!=null)         this.inf_iters=Integer.parseInt(properties.getProperty(prefix+"inf_iters"));
  			if (properties.getProperty(prefix+"inf_final_diff")!=null)    this.inf_final_diff=Double.parseDouble(properties.getProperty(prefix+"inf_final_diff"));
  			if (properties.getProperty(prefix+"inf_far_pull")!=null)      this.inf_far_pull=Double.parseDouble(properties.getProperty(prefix+"inf_far_pull"));
  			
  			if (properties.getProperty(prefix+"inf_str_pow")!=null)       this.inf_str_pow=Double.parseDouble(properties.getProperty(prefix+"inf_str_pow"));
  			if (properties.getProperty(prefix+"inf_smpl_side")!=null)     this.inf_smpl_side=Integer.parseInt(properties.getProperty(prefix+"inf_smpl_side"));
  			if (properties.getProperty(prefix+"inf_smpl_num")!=null)      this.inf_smpl_num=Integer.parseInt(properties.getProperty(prefix+"inf_smpl_num"));
  			if (properties.getProperty(prefix+"inf_smpl_rms")!=null)      this.inf_smpl_rms=Double.parseDouble(properties.getProperty(prefix+"inf_smpl_rms"));
  			
  			if (properties.getProperty(prefix+"ih_smpl_step")!=null)      this.ih_smpl_step=Integer.parseInt(properties.getProperty(prefix+"ih_smpl_step"));
  			if (properties.getProperty(prefix+"ih_disp_min")!=null)       this.ih_disp_min=Double.parseDouble(properties.getProperty(prefix+"ih_disp_min"));
  			if (properties.getProperty(prefix+"ih_disp_step")!=null)      this.ih_disp_step=Double.parseDouble(properties.getProperty(prefix+"ih_disp_step"));
  			if (properties.getProperty(prefix+"ih_num_bins")!=null)       this.ih_num_bins=Integer.parseInt(properties.getProperty(prefix+"ih_num_bins"));
  			if (properties.getProperty(prefix+"ih_sigma")!=null)          this.ih_sigma=Double.parseDouble(properties.getProperty(prefix+"ih_sigma"));
  			if (properties.getProperty(prefix+"ih_max_diff")!=null)       this.ih_max_diff=Double.parseDouble(properties.getProperty(prefix+"ih_max_diff"));
  			if (properties.getProperty(prefix+"ih_min_samples")!=null)    this.ih_min_samples=Integer.parseInt(properties.getProperty(prefix+"ih_min_samples"));
  			if (properties.getProperty(prefix+"ih_norm_center")!=null)    this.ih_norm_center=Boolean.parseBoolean(properties.getProperty(prefix+"ih_norm_center"));
  			
			if (properties.getProperty(prefix+"ly_smpl_side")!=null)      this.ly_smpl_side=Integer.parseInt(properties.getProperty(prefix+"ly_smpl_side"));
			if (properties.getProperty(prefix+"ly_smpl_num")!=null)       this.ly_smpl_num=Integer.parseInt(properties.getProperty(prefix+"ly_smpl_num"));
			if (properties.getProperty(prefix+"ly_meas_disp")!=null)      this.ly_meas_disp=Double.parseDouble(properties.getProperty(prefix+"ly_meas_disp"));
			if (properties.getProperty(prefix+"ly_smpl_rms")!=null)       this.ly_smpl_rms=Double.parseDouble(properties.getProperty(prefix+"ly_smpl_rms"));
			if (properties.getProperty(prefix+"ly_disp_var")!=null)       this.ly_disp_var=Double.parseDouble(properties.getProperty(prefix+"ly_disp_var"));
			if (properties.getProperty(prefix+"ly_inf_frac")!=null)       this.ly_inf_frac=Double.parseDouble(properties.getProperty(prefix+"ly_inf_frac"));
  			if (properties.getProperty(prefix+"ly_on_scan")!=null)        this.ly_on_scan=Boolean.parseBoolean(properties.getProperty(prefix+"ly_on_scan"));
  			if (properties.getProperty(prefix+"ly_inf_en")!=null)         this.ly_inf_en=Boolean.parseBoolean(properties.getProperty(prefix+"ly_inf_en"));
  			if (properties.getProperty(prefix+"ly_inf_force")!=null)      this.ly_inf_force=Boolean.parseBoolean(properties.getProperty(prefix+"ly_inf_force"));
  			if (properties.getProperty(prefix+"ly_poly")!=null)           this.ly_poly=Boolean.parseBoolean(properties.getProperty(prefix+" "));
 			
  			if (properties.getProperty(prefix+"corr_magic_scale")!=null)  this.corr_magic_scale=Double.parseDouble(properties.getProperty(prefix+"corr_magic_scale"));

  			if (properties.getProperty(prefix+"show_textures")!=null)      this.show_textures=Boolean.parseBoolean(properties.getProperty(prefix+"show_textures"));
  			if (properties.getProperty(prefix+"debug_filters")!=null)      this.debug_filters=Boolean.parseBoolean(properties.getProperty(prefix+"debug_filters"));

  			if (properties.getProperty(prefix+"min_smth")!=null)          this.min_smth=Double.parseDouble(properties.getProperty(prefix+"min_smth"));
  			if (properties.getProperty(prefix+"sure_smth")!=null)         this.sure_smth=Double.parseDouble(properties.getProperty(prefix+"sure_smth"));
  			if (properties.getProperty(prefix+"bgnd_range")!=null)        this.bgnd_range=Double.parseDouble(properties.getProperty(prefix+"bgnd_range"));
  			if (properties.getProperty(prefix+"other_range")!=null)       this.other_range=Double.parseDouble(properties.getProperty(prefix+"other_range"));

  			if (properties.getProperty(prefix+"ex_strength")!=null)       this.ex_strength=Double.parseDouble(properties.getProperty(prefix+"ex_strength"));
  			if (properties.getProperty(prefix+"ex_nstrength")!=null)      this.ex_nstrength=Double.parseDouble(properties.getProperty(prefix+"ex_nstrength"));

  			if (properties.getProperty(prefix+"ex_over_bgnd")!=null)      this.ex_over_bgnd=Boolean.parseBoolean(properties.getProperty(prefix+"ex_over_bgnd"));
  			if (properties.getProperty(prefix+"ex_min_over")!=null)       this.ex_min_over=Double.parseDouble(properties.getProperty(prefix+"ex_min_over"));

  			if (properties.getProperty(prefix+"pt_super_trust")!=null)    this.pt_super_trust=Double.parseDouble(properties.getProperty(prefix+"pt_super_trust"));
  			if (properties.getProperty(prefix+"pt_keep_raw_fg")!=null)    this.pt_keep_raw_fg=Boolean.parseBoolean(properties.getProperty(prefix+"pt_keep_raw_fg"));
  			if (properties.getProperty(prefix+"pt_scale_pre")!=null)      this.pt_scale_pre=Double.parseDouble(properties.getProperty(prefix+"pt_scale_pre"));
  			if (properties.getProperty(prefix+"pt_scale_post")!=null)     this.pt_scale_post=Double.parseDouble(properties.getProperty(prefix+"pt_scale_post"));

  			if (properties.getProperty(prefix+"bgnd_sure")!=null)         this.bgnd_sure=Double.parseDouble(properties.getProperty(prefix+"bgnd_sure"));
  			if (properties.getProperty(prefix+"bgnd_maybe")!=null)        this.bgnd_maybe=Double.parseDouble(properties.getProperty(prefix+"bgnd_maybe"));
  			if (properties.getProperty(prefix+"min_clstr_seed")!=null)    this.min_clstr_seed=Integer.parseInt(properties.getProperty(prefix+"min_clstr_seed"));
  			if (properties.getProperty(prefix+"min_clstr_lone")!=null)    this.min_clstr_lone=Integer.parseInt(properties.getProperty(prefix+"min_clstr_lone"));
  			if (properties.getProperty(prefix+"min_clstr_weight")!=null)  this.min_clstr_weight=Double.parseDouble(properties.getProperty(prefix+"min_clstr_weight"));
  			if (properties.getProperty(prefix+"min_clstr_max")!=null)     this.min_clstr_max=Double.parseDouble(properties.getProperty(prefix+"min_clstr_max"));

  			
  			if (properties.getProperty(prefix+"fill_gaps")!=null)         this.fill_gaps=Integer.parseInt(properties.getProperty(prefix+"fill_gaps"));
  			if (properties.getProperty(prefix+"fill_final")!=null)        this.fill_final=Integer.parseInt(properties.getProperty(prefix+"fill_final"));
  			if (properties.getProperty(prefix+"min_clstr_block")!=null)   this.min_clstr_block=Integer.parseInt(properties.getProperty(prefix+"min_clstr_block"));
  			if (properties.getProperty(prefix+"bgnd_grow")!=null)         this.bgnd_grow=Integer.parseInt(properties.getProperty(prefix+"bgnd_grow"));

  			if (properties.getProperty(prefix+"ortho_old")!=null)         this.ortho_old=Boolean.parseBoolean(properties.getProperty(prefix+"ortho_old"));
  			if (properties.getProperty(prefix+"ortho_min_hor")!=null)     this.ortho_min_hor=Double.parseDouble(properties.getProperty(prefix+"ortho_min_hor"));
  			if (properties.getProperty(prefix+"ortho_min_vert")!=null)    this.ortho_min_vert=Double.parseDouble(properties.getProperty(prefix+"ortho_min_vert"));
  			if (properties.getProperty(prefix+"ortho_asym")!=null)        this.ortho_asym=Double.parseDouble(properties.getProperty(prefix+"ortho_asym"));
  			if (properties.getProperty(prefix+"ortho_over4")!=null)       this.ortho_over4=Double.parseDouble(properties.getProperty(prefix+"ortho_over4"));
  			if (properties.getProperty(prefix+"ortho_sustain")!=null)     this.ortho_sustain=Double.parseDouble(properties.getProperty(prefix+"ortho_sustain"));
  			if (properties.getProperty(prefix+"ortho_run")!=null)         this.ortho_run=Integer.parseInt(properties.getProperty(prefix+"ortho_run"));
  			if (properties.getProperty(prefix+"ortho_minmax")!=null)      this.ortho_minmax=Double.parseDouble(properties.getProperty(prefix+"ortho_minmax"));
  			if (properties.getProperty(prefix+"ortho_bridge")!=null)      this.ortho_bridge=Integer.parseInt(properties.getProperty(prefix+"ortho_bridge"));
  			if (properties.getProperty(prefix+"ortho_rms")!=null)         this.ortho_rms=Double.parseDouble(properties.getProperty(prefix+"ortho_rms"));
  			if (properties.getProperty(prefix+"ortho_half_length")!=null) this.ortho_half_length=Integer.parseInt(properties.getProperty(prefix+"ortho_half_length"));
  			if (properties.getProperty(prefix+"ortho_mix")!=null)         this.ortho_mix=Double.parseDouble(properties.getProperty(prefix+"ortho_mix"));
  			
  			if (properties.getProperty(prefix+"or_hor")!=null)            this.or_hor=Boolean.parseBoolean(properties.getProperty(prefix+"or_hor"));
  			if (properties.getProperty(prefix+"or_vert")!=null)           this.or_vert=Boolean.parseBoolean(properties.getProperty(prefix+"or_vert"));
  			if (properties.getProperty(prefix+"or_sigma")!=null)          this.or_sigma=Double.parseDouble(properties.getProperty(prefix+"or_sigma"));
  			if (properties.getProperty(prefix+"or_sharp")!=null)          this.or_sharp=Double.parseDouble(properties.getProperty(prefix+"or_sharp"));
  			if (properties.getProperty(prefix+"or_scale")!=null)          this.or_scale=Double.parseDouble(properties.getProperty(prefix+"or_scale"));
  			if (properties.getProperty(prefix+"or_offset")!=null)         this.or_offset=Double.parseDouble(properties.getProperty(prefix+"or_offset"));
  			if (properties.getProperty(prefix+"or_asym")!=null)           this.or_asym=Double.parseDouble(properties.getProperty(prefix+"or_asym"));
  			if (properties.getProperty(prefix+"or_threshold")!=null)      this.or_threshold=Double.parseDouble(properties.getProperty(prefix+"or_threshold"));
  			if (properties.getProperty(prefix+"or_absHor")!=null)         this.or_absHor=Double.parseDouble(properties.getProperty(prefix+"or_absHor"));
  			if (properties.getProperty(prefix+"or_absVert")!=null)        this.or_absVert=Double.parseDouble(properties.getProperty(prefix+"or_absVert"));

  			if (properties.getProperty(prefix+"poles_fix")!=null)         this.poles_fix=Boolean.parseBoolean(properties.getProperty(prefix+"poles_fix"));
  			if (properties.getProperty(prefix+"poles_len")!=null)         this.poles_len=Integer.parseInt(properties.getProperty(prefix+"poles_len"));
  			if (properties.getProperty(prefix+"poles_ratio")!=null)       this.poles_ratio=Double.parseDouble(properties.getProperty(prefix+"poles_ratio"));
  			if (properties.getProperty(prefix+"poles_min_strength")!=null)this.poles_min_strength=Double.parseDouble(properties.getProperty(prefix+"poles_min_strength"));
  			if (properties.getProperty(prefix+"poles_force_disp")!=null)  this.poles_force_disp=Boolean.parseBoolean(properties.getProperty(prefix+"poles_force_disp"));
  			
  			if (properties.getProperty(prefix+"max_clusters")!=null)      this.max_clusters=Integer.parseInt(properties.getProperty(prefix+"max_clusters"));
  			if (properties.getProperty(prefix+"remove_scans")!=null)      this.remove_scans=Boolean.parseBoolean(properties.getProperty(prefix+"remove_scans"));
  			if (properties.getProperty(prefix+"output_x3d")!=null)        this.output_x3d=Boolean.parseBoolean(properties.getProperty(prefix+"output_x3d"));
  			if (properties.getProperty(prefix+"output_obj")!=null)        this.output_obj=Boolean.parseBoolean(properties.getProperty(prefix+"output_obj"));
  			if (properties.getProperty(prefix+"correct_distortions")!=null) this.correct_distortions=Boolean.parseBoolean(properties.getProperty(prefix+"correct_distortions"));
  			if (properties.getProperty(prefix+"show_triangles")!=null)    this.show_triangles=Boolean.parseBoolean(properties.getProperty(prefix+"show_triangles"));
  			if (properties.getProperty(prefix+"avg_cluster_disp")!=null)  this.avg_cluster_disp=Boolean.parseBoolean(properties.getProperty(prefix+"avg_cluster_disp"));
  			if (properties.getProperty(prefix+"maxDispTriangle")!=null)   this.maxDispTriangle=Double.parseDouble(properties.getProperty(prefix+"maxDispTriangle"));
  			if (properties.getProperty(prefix+"infinityDistance")!=null)  this.infinityDistance=Double.parseDouble(properties.getProperty(prefix+"infinityDistance"));
  			if (properties.getProperty(prefix+"shUseFlaps")!=null)        this.shUseFlaps=Boolean.parseBoolean(properties.getProperty(prefix+"shUseFlaps"));
  			if (properties.getProperty(prefix+"shAggrFade")!=null)        this.shAggrFade=Boolean.parseBoolean(properties.getProperty(prefix+"shAggrFade"));
  			if (properties.getProperty(prefix+"shMinArea")!=null)         this.shMinArea=Integer.parseInt(properties.getProperty(prefix+"shMinArea"));
  			if (properties.getProperty(prefix+"shMinStrength")!=null)    this.shMinStrength=Double.parseDouble(properties.getProperty(prefix+"shMinStrength"));
  			if (properties.getProperty(prefix+"tiRigidVertical")!=null)   this.tiRigidVertical=Double.parseDouble(properties.getProperty(prefix+"tiRigidVertical"));
  			if (properties.getProperty(prefix+"tiRigidHorizontal")!=null) this.tiRigidHorizontal=Double.parseDouble(properties.getProperty(prefix+"tiRigidHorizontal"));
  			if (properties.getProperty(prefix+"tiRigidDiagonal")!=null)   this.tiRigidDiagonal=Double.parseDouble(properties.getProperty(prefix+"tiRigidDiagonal"));
  			if (properties.getProperty(prefix+"tiStrengthOffset")!=null)  this.tiStrengthOffset=Double.parseDouble(properties.getProperty(prefix+"tiStrengthOffset"));
  			if (properties.getProperty(prefix+"tiDispScale")!=null)       this.tiDispScale=Double.parseDouble(properties.getProperty(prefix+"tiDispScale"));
  			if (properties.getProperty(prefix+"tiDispPow")!=null)         this.tiDispPow=Double.parseDouble(properties.getProperty(prefix+"tiDispPow"));
  			if (properties.getProperty(prefix+"tiDispPull")!=null)        this.tiDispPull=Double.parseDouble(properties.getProperty(prefix+"tiDispPull"));
  			if (properties.getProperty(prefix+"tiDispPullPreFinal")!=null)this.tiDispPullPreFinal=Double.parseDouble(properties.getProperty(prefix+"tiDispPullPreFinal"));
  			if (properties.getProperty(prefix+"tiDispPullFinal")!=null)   this.tiDispPullFinal=Double.parseDouble(properties.getProperty(prefix+"tiDispPullFinal"));
  			if (properties.getProperty(prefix+"tiBreakNorm")!=null)       this.tiBreakNorm=Double.parseDouble(properties.getProperty(prefix+"tiBreakNorm"));
  			if (properties.getProperty(prefix+"tiBreak3")!=null)          this.tiBreak3=Double.parseDouble(properties.getProperty(prefix+"tiBreak3"));
  			if (properties.getProperty(prefix+"tiBreak31")!=null)         this.tiBreak31=Double.parseDouble(properties.getProperty(prefix+"tiBreak31"));
  			if (properties.getProperty(prefix+"tiBreak21")!=null)         this.tiBreak21=Double.parseDouble(properties.getProperty(prefix+"tiBreak21"));
  			if (properties.getProperty(prefix+"tiBreakFar")!=null)        this.tiBreakFar=Double.parseDouble(properties.getProperty(prefix+"tiBreakFar"));
  			if (properties.getProperty(prefix+"tiBreakNear")!=null)       this.tiBreakNear=Double.parseDouble(properties.getProperty(prefix+"tiBreakNear"));
  			if (properties.getProperty(prefix+"tiBreakMode")!=null)       this.tiBreakMode=Integer.parseInt(properties.getProperty(prefix+"tiBreakMode"));
  			if (properties.getProperty(prefix+"tiBreakSame")!=null)       this.tiBreakSame=Double.parseDouble(properties.getProperty(prefix+"tiBreakSame"));
  			if (properties.getProperty(prefix+"tiBreakTurn")!=null)       this.tiBreakTurn=Double.parseDouble(properties.getProperty(prefix+"tiBreakTurn"));

  			if (properties.getProperty(prefix+"tiHealPreLast")!=null)     this.tiHealPreLast=Double.parseDouble(properties.getProperty(prefix+"tiHealPreLast"));
  			if (properties.getProperty(prefix+"tiHealLast")!=null)        this.tiHealLast=Double.parseDouble(properties.getProperty(prefix+"tiHealLast"));
  			if (properties.getProperty(prefix+"tiHealSame")!=null)        this.tiHealSame=Integer.parseInt(properties.getProperty(prefix+"tiHealSame"));

  			if (properties.getProperty(prefix+"tiIterations")!=null)      this.tiIterations=Integer.parseInt(properties.getProperty(prefix+"tiIterations"));
  			if (properties.getProperty(prefix+"tiPrecision")!=null)       this.tiPrecision=Integer.parseInt(properties.getProperty(prefix+"tiPrecision"));
  			if (properties.getProperty(prefix+"tiNumCycles")!=null)       this.tiNumCycles=Integer.parseInt(properties.getProperty(prefix+"tiNumCycles"));

  			if (properties.getProperty(prefix+"stUseRefine")!=null)       this.stUseRefine=Boolean.parseBoolean(properties.getProperty(prefix+"stUseRefine"));
  			if (properties.getProperty(prefix+"stUsePass2")!=null)        this.stUsePass2=Boolean.parseBoolean(properties.getProperty(prefix+"stUsePass2"));
  			if (properties.getProperty(prefix+"stUseRender")!=null)       this.stUseRender=Boolean.parseBoolean(properties.getProperty(prefix+"stUseRender"));
  			
  			if (properties.getProperty(prefix+"stShow")!=null)            this.stShow=Boolean.parseBoolean(properties.getProperty(prefix+"stShow"));
  			if (properties.getProperty(prefix+"stSize")!=null)            this.stSize=Integer.parseInt(properties.getProperty(prefix+"stSize"));
  			if (properties.getProperty(prefix+"stStepFar")!=null)         this.stStepFar=Double.parseDouble(properties.getProperty(prefix+"stStepFar"));
  			if (properties.getProperty(prefix+"stStepNear")!=null)        this.stStepNear=Double.parseDouble(properties.getProperty(prefix+"stStepNear"));
  			if (properties.getProperty(prefix+"stStepThreshold")!=null)   this.stStepThreshold=Double.parseDouble(properties.getProperty(prefix+"stStepThreshold"));
  			if (properties.getProperty(prefix+"stMinDisparity")!=null)    this.stMinDisparity=Double.parseDouble(properties.getProperty(prefix+"stMinDisparity"));
//  			if (properties.getProperty(prefix+"stMaxDisparity")!=null)    this.stMaxDisparity=Double.parseDouble(properties.getProperty(prefix+"stMaxDisparity"));
  			if (properties.getProperty(prefix+"stFloor")!=null)           this.stFloor=Double.parseDouble(properties.getProperty(prefix+"stFloor"));
  			if (properties.getProperty(prefix+"stPow")!=null)             this.stPow=Double.parseDouble(properties.getProperty(prefix+"stPow"));
  			if (properties.getProperty(prefix+"stSigma")!=null)           this.stSigma=Double.parseDouble(properties.getProperty(prefix+"stSigma"));
  			if (properties.getProperty(prefix+"stMinBgDisparity")!=null)  this.stMinBgDisparity=Double.parseDouble(properties.getProperty(prefix+"stMinBgDisparity"));
  			if (properties.getProperty(prefix+"stMinBgFract")!=null)      this.stMinBgFract=Double.parseDouble(properties.getProperty(prefix+"stMinBgFract"));
  			if (properties.getProperty(prefix+"stUseDisp")!=null)         this.stUseDisp=Double.parseDouble(properties.getProperty(prefix+"stUseDisp"));
  			if (properties.getProperty(prefix+"stStrengthScale")!=null)   this.stStrengthScale=Double.parseDouble(properties.getProperty(prefix+"stStrengthScale"));

  			if (properties.getProperty(prefix+"stSmplMode")!=null)        this.stSmplMode=Boolean.parseBoolean(properties.getProperty(prefix+"stSmplMode"));
  			if (properties.getProperty(prefix+"stSmplSide")!=null)        this.stSmplSide=Integer.parseInt(properties.getProperty(prefix+"stSmplSide"));
  			if (properties.getProperty(prefix+"stSmplNum")!=null)         this.stSmplNum=Integer.parseInt(properties.getProperty(prefix+"stSmplNum"));
  			if (properties.getProperty(prefix+"stSmplRms")!=null)         this.stSmplRms=Double.parseDouble(properties.getProperty(prefix+"stSmplRms"));
  			if (properties.getProperty(prefix+"stSmplWnd")!=null)         this.stSmplWnd=Boolean.parseBoolean(properties.getProperty(prefix+"stSmplWnd"));

  			if (properties.getProperty(prefix+"stGrowSel")!=null)         this.stGrowSel=Integer.parseInt(properties.getProperty(prefix+"stGrowSel"));
  			if (properties.getProperty(prefix+"stMeasSel")!=null)         this.stMeasSel=Integer.parseInt(properties.getProperty(prefix+"stMeasSel"));
  			if (properties.getProperty(prefix+"stSmallDiff")!=null)       this.stSmallDiff=Double.parseDouble(properties.getProperty(prefix+"stSmallDiff"));
  			if (properties.getProperty(prefix+"stHighMix")!=null)         this.stHighMix=Double.parseDouble(properties.getProperty(prefix+"stHighMix"));

  			if (properties.getProperty(prefix+"outlayerStrength")!=null)  this.outlayerStrength=Double.parseDouble(properties.getProperty(prefix+"outlayerStrength"));
  			if (properties.getProperty(prefix+"outlayerDiff")!=null)      this.outlayerDiff=Double.parseDouble(properties.getProperty(prefix+"outlayerDiff"));
  			if (properties.getProperty(prefix+"outlayerDiffPos")!=null)   this.outlayerDiffPos=Double.parseDouble(properties.getProperty(prefix+"outlayerDiffPos"));
  			if (properties.getProperty(prefix+"outlayerDiffNeg")!=null)   this.outlayerDiffNeg=Double.parseDouble(properties.getProperty(prefix+"outlayerDiffNeg"));

  			if (properties.getProperty(prefix+"combine_refine")!=null)    this.combine_refine=Boolean.parseBoolean(properties.getProperty(prefix+"combine_refine"));

  			if (properties.getProperty(prefix+"combine_min_strength")!=null) this.combine_min_strength=Double.parseDouble(properties.getProperty(prefix+"combine_min_strength"));
  			if (properties.getProperty(prefix+"combine_min_hor")!=null)   this.combine_min_hor=Double.parseDouble(properties.getProperty(prefix+"combine_min_hor"));
  			if (properties.getProperty(prefix+"combine_min_vert")!=null)  this.combine_min_vert=Double.parseDouble(properties.getProperty(prefix+"combine_min_vert"));
//  			if (properties.getProperty(prefix+"unique_tolerance")!=null)  this.unique_tolerance=Double.parseDouble(properties.getProperty(prefix+"unique_tolerance"));
  			if (properties.getProperty(prefix+"grow_sweep")!=null)        this.grow_sweep=Integer.parseInt(properties.getProperty(prefix+"grow_sweep"));
  			if (properties.getProperty(prefix+"grow_disp_max")!=null)     this.grow_disp_max=Double.parseDouble(properties.getProperty(prefix+"grow_disp_max"));
  			if (properties.getProperty(prefix+"grow_disp_trust")!=null)   this.grow_disp_trust=Double.parseDouble(properties.getProperty(prefix+"grow_disp_trust"));
  			if (properties.getProperty(prefix+"grow_disp_step")!=null)    this.grow_disp_step=Double.parseDouble(properties.getProperty(prefix+"grow_disp_step"));
  			if (properties.getProperty(prefix+"grow_min_diff")!=null)     this.grow_min_diff=Double.parseDouble(properties.getProperty(prefix+"grow_min_diff"));
  			if (properties.getProperty(prefix+"grow_retry_far")!=null)    this.grow_retry_far=Boolean.parseBoolean(properties.getProperty(prefix+"grow_retry_far"));
  			if (properties.getProperty(prefix+"grow_pedantic")!=null)     this.grow_pedantic=Boolean.parseBoolean(properties.getProperty(prefix+"grow_pedantic"));
  			if (properties.getProperty(prefix+"grow_retry_inf")!=null)    this.grow_retry_inf=Boolean.parseBoolean(properties.getProperty(prefix+"grow_retry_inf"));
  			
  			if (properties.getProperty(prefix+"gr_new_expand")!=null)     this.gr_new_expand=Boolean.parseBoolean(properties.getProperty(prefix+"gr_new_expand"));
  			if (properties.getProperty(prefix+"gr_max_expand")!=null)     this.gr_max_expand=Integer.parseInt(properties.getProperty(prefix+"gr_max_expand"));
  			if (properties.getProperty(prefix+"fds_str_floor")!=null)     this.fds_str_floor=Double.parseDouble(properties.getProperty(prefix+"fds_str_floor"));
  			if (properties.getProperty(prefix+"gr_ovrbg_cmb")!=null)      this.gr_ovrbg_cmb=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_cmb"));
  			if (properties.getProperty(prefix+"gr_ovrbg_cmb_hor")!=null)  this.gr_ovrbg_cmb_hor=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_cmb_hor"));
  			if (properties.getProperty(prefix+"gr_ovrbg_cmb_vert")!=null) this.gr_ovrbg_cmb_vert=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_cmb_vert"));
  			if (properties.getProperty(prefix+"gr_ovrbg_filtered")!=null) this.gr_ovrbg_filtered=Double.parseDouble(properties.getProperty(prefix+"gr_ovrbg_filtered"));
  			if (properties.getProperty(prefix+"fds_str_pow")!=null)       this.fds_str_pow=Double.parseDouble(properties.getProperty(prefix+"fds_str_pow"));
  			if (properties.getProperty(prefix+"fds_smpl_side")!=null)     this.fds_smpl_side=Integer.parseInt(properties.getProperty(prefix+"fds_smpl_side"));
  			if (properties.getProperty(prefix+"fds_smpl_num")!=null)      this.fds_smpl_num=Integer.parseInt(properties.getProperty(prefix+"fds_smpl_num"));
  			if (properties.getProperty(prefix+"fds_smpl_rms")!=null)      this.fds_smpl_rms=Double.parseDouble(properties.getProperty(prefix+"fds_smpl_rms"));
  			if (properties.getProperty(prefix+"fds_smpl_rel_rms")!=null)  this.fds_smpl_rel_rms=Double.parseDouble(properties.getProperty(prefix+"fds_smpl_rel_rms"));
  			if (properties.getProperty(prefix+"fds_smpl_wnd")!=null)      this.fds_smpl_wnd=Boolean.parseBoolean(properties.getProperty(prefix+"fds_smpl_wnd"));
  			if (properties.getProperty(prefix+"fds_abs_tilt")!=null)      this.fds_abs_tilt=Double.parseDouble(properties.getProperty(prefix+"fds_abs_tilt"));
  			if (properties.getProperty(prefix+"fds_rel_tilt")!=null)      this.fds_rel_tilt=Double.parseDouble(properties.getProperty(prefix+"fds_rel_tilt"));

  			if (properties.getProperty(prefix+"mc_disp8_step")!=null)     this.mc_disp8_step=Double.parseDouble(properties.getProperty(prefix+"mc_disp8_step"));
  			if (properties.getProperty(prefix+"mc_disp8_trust")!=null)    this.mc_disp8_trust=Double.parseDouble(properties.getProperty(prefix+"mc_disp8_trust"));
  			if (properties.getProperty(prefix+"mc_strength")!=null)       this.mc_strength=Double.parseDouble(properties.getProperty(prefix+"mc_strength"));
  			if (properties.getProperty(prefix+"mc_unique_tol")!=null)     this.mc_unique_tol=Double.parseDouble(properties.getProperty(prefix+"mc_unique_tol"));
  			if (properties.getProperty(prefix+"mc_trust_fin")!=null)      this.mc_trust_fin=Double.parseDouble(properties.getProperty(prefix+"mc_trust_fin"));
  			if (properties.getProperty(prefix+"mc_trust_sigma")!=null)    this.mc_trust_sigma=Double.parseDouble(properties.getProperty(prefix+"mc_trust_sigma"));
  			if (properties.getProperty(prefix+"mc_ortho_weight")!=null)   this.mc_ortho_weight=Double.parseDouble(properties.getProperty(prefix+"mc_ortho_weight"));
  			if (properties.getProperty(prefix+"mc_diag_weight")!=null)    this.mc_diag_weight=Double.parseDouble(properties.getProperty(prefix+"mc_diag_weight"));
  			if (properties.getProperty(prefix+"mc_gap")!=null)            this.mc_gap=Double.parseDouble(properties.getProperty(prefix+"mc_gap"));

  			if (properties.getProperty(prefix+"gr_min_new")!=null)        this.gr_min_new=Integer.parseInt(properties.getProperty(prefix+"gr_min_new"));
  			if (properties.getProperty(prefix+"gr_var_new_sngl")!=null)   this.gr_var_new_sngl=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_new_sngl"));
  			if (properties.getProperty(prefix+"gr_var_new_fg")!=null)     this.gr_var_new_fg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_new_fg"));
  			if (properties.getProperty(prefix+"gr_var_all_fg")!=null)     this.gr_var_all_fg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_all_fg"));
  			if (properties.getProperty(prefix+"gr_var_new_bg")!=null)     this.gr_var_new_bg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_new_bg"));
  			if (properties.getProperty(prefix+"gr_var_all_bg")!=null)     this.gr_var_all_bg=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_all_bg"));
  			if (properties.getProperty(prefix+"gr_var_next")!=null)       this.gr_var_next=Boolean.parseBoolean(properties.getProperty(prefix+"gr_var_next"));
  			if (properties.getProperty(prefix+"gr_num_steps")!=null)      this.gr_num_steps=Integer.parseInt(properties.getProperty(prefix+"gr_num_steps"));
  			if (properties.getProperty(prefix+"gr_steps_over")!=null)     this.gr_steps_over=Integer.parseInt(properties.getProperty(prefix+"gr_steps_over"));
  			if (properties.getProperty(prefix+"gr_smpl_size")!=null)      this.gr_smpl_size=Integer.parseInt(properties.getProperty(prefix+"gr_smpl_size"));
  			if (properties.getProperty(prefix+"gr_min_pnts")!=null)       this.gr_min_pnts=Integer.parseInt(properties.getProperty(prefix+"gr_min_pnts"));
  			if (properties.getProperty(prefix+"gr_use_wnd")!=null)        this.gr_use_wnd=Boolean.parseBoolean(properties.getProperty(prefix+"gr_use_wnd"));
  			if (properties.getProperty(prefix+"gr_tilt_damp")!=null)      this.gr_tilt_damp=Double.parseDouble(properties.getProperty(prefix+"gr_tilt_damp"));
  			if (properties.getProperty(prefix+"gr_split_rng")!=null)      this.gr_split_rng=Double.parseDouble(properties.getProperty(prefix+"gr_split_rng"));
  			if (properties.getProperty(prefix+"gr_same_rng")!=null)       this.gr_same_rng=Double.parseDouble(properties.getProperty(prefix+"gr_same_rng"));
  			if (properties.getProperty(prefix+"gr_diff_cont")!=null)      this.gr_diff_cont=Double.parseDouble(properties.getProperty(prefix+"gr_diff_cont"));
  			if (properties.getProperty(prefix+"gr_abs_tilt")!=null)       this.gr_abs_tilt=Double.parseDouble(properties.getProperty(prefix+"gr_abs_tilt"));
  			if (properties.getProperty(prefix+"gr_rel_tilt")!=null)       this.gr_rel_tilt=Double.parseDouble(properties.getProperty(prefix+"gr_rel_tilt"));
  			if (properties.getProperty(prefix+"gr_smooth")!=null)         this.gr_smooth=Integer.parseInt(properties.getProperty(prefix+"gr_smooth"));
  			if (properties.getProperty(prefix+"gr_fin_diff")!=null)       this.gr_fin_diff=Double.parseDouble(properties.getProperty(prefix+"gr_fin_diff"));
  			if (properties.getProperty(prefix+"gr_unique_tol")!=null)     this.gr_unique_tol=Double.parseDouble(properties.getProperty(prefix+"gr_unique_tol"));
  			if (properties.getProperty(prefix+"gr_unique_pretol")!=null)  this.gr_unique_pretol=Double.parseDouble(properties.getProperty(prefix+"gr_unique_pretol"));
  			
  			if (properties.getProperty(prefix+"plPreferDisparity")!=null) this.plPreferDisparity=Boolean.parseBoolean(properties.getProperty(prefix+"plPreferDisparity"));
  			if (properties.getProperty(prefix+"plDispNorm")!=null)        this.plDispNorm=Double.parseDouble(properties.getProperty(prefix+"plDispNorm"));

  			if (properties.getProperty(prefix+"plBlurBinVert")!=null)     this.plBlurBinVert=Double.parseDouble(properties.getProperty(prefix+"plBlurBinVert"));
  			if (properties.getProperty(prefix+"plBlurBinHor")!=null)      this.plBlurBinHor=Double.parseDouble(properties.getProperty(prefix+"plBlurBinHor"));
  			if (properties.getProperty(prefix+"plMaxDiffVert")!=null)     this.plMaxDiffVert=Double.parseDouble(properties.getProperty(prefix+"plMaxDiffVert"));
  			if (properties.getProperty(prefix+"plMaxDiffHor")!=null)      this.plMaxDiffHor=Double.parseDouble(properties.getProperty(prefix+"plMaxDiffHor"));
  			if (properties.getProperty(prefix+"plInitPasses")!=null)      this.plInitPasses=Integer.parseInt(properties.getProperty(prefix+"plInitPasses"));
  			
  			if (properties.getProperty(prefix+"plMinPoints")!=null)       this.plMinPoints=Integer.parseInt(properties.getProperty(prefix+"plMinPoints"));
  			if (properties.getProperty(prefix+"plTargetEigen")!=null)     this.plTargetEigen=Double.parseDouble(properties.getProperty(prefix+"plTargetEigen"));
  			if (properties.getProperty(prefix+"plFractOutliers")!=null)   this.plFractOutliers=Double.parseDouble(properties.getProperty(prefix+"plFractOutliers"));
  			if (properties.getProperty(prefix+"plMaxOutliers")!=null)     this.plMaxOutliers=Integer.parseInt(properties.getProperty(prefix+"plMaxOutliers"));
  			if (properties.getProperty(prefix+"plMinStrength")!=null)     this.plMinStrength=Double.parseDouble(properties.getProperty(prefix+"plMinStrength"));
  			if (properties.getProperty(prefix+"plMaxEigen")!=null)        this.plMaxEigen=Double.parseDouble(properties.getProperty(prefix+"plMaxEigen"));
  			if (properties.getProperty(prefix+"plEigenFloor")!=null)      this.plEigenFloor=Double.parseDouble(properties.getProperty(prefix+"plEigenFloor"));
  			if (properties.getProperty(prefix+"plEigenStick")!=null)      this.plEigenStick=Double.parseDouble(properties.getProperty(prefix+"plEigenStick"));
  			if (properties.getProperty(prefix+"plBadPlate")!=null)        this.plBadPlate=Double.parseDouble(properties.getProperty(prefix+"plBadPlate"));
  			if (properties.getProperty(prefix+"plDbgMerge")!=null)        this.plDbgMerge=Boolean.parseBoolean(properties.getProperty(prefix+"plDbgMerge"));
  			if (properties.getProperty(prefix+"plWorstWorsening")!=null)  this.plWorstWorsening=Double.parseDouble(properties.getProperty(prefix+"plWorstWorsening"));
  			if (properties.getProperty(prefix+"plWorstWorsening2")!=null) this.plWorstWorsening2=Double.parseDouble(properties.getProperty(prefix+"plWorstWorsening2"));
  			if (properties.getProperty(prefix+"plWorstEq")!=null)         this.plWorstEq=Double.parseDouble(properties.getProperty(prefix+"plWorstEq"));
  			if (properties.getProperty(prefix+"plWorstEq2")!=null)        this.plWorstEq2=Double.parseDouble(properties.getProperty(prefix+"plWorstEq2"));
  			if (properties.getProperty(prefix+"plOKMergeEigen")!=null)    this.plOKMergeEigen=Double.parseDouble(properties.getProperty(prefix+"plOKMergeEigen"));
  			if (properties.getProperty(prefix+"plMaxWorldSin2")!=null)    this.plMaxWorldSin2=Double.parseDouble(properties.getProperty(prefix+"plMaxWorldSin2"));
  			if (properties.getProperty(prefix+"pl2dForSin")!=null)        this.pl2dForSin=Double.parseDouble(properties.getProperty(prefix+"pl2dForSin"));
  			if (properties.getProperty(prefix+"plWeakWorsening")!=null)   this.plWeakWorsening=Double.parseDouble(properties.getProperty(prefix+"plWeakWorsening"));
  			if (properties.getProperty(prefix+"plMaxOverlap")!=null)      this.plMaxOverlap=Double.parseDouble(properties.getProperty(prefix+"plMaxOverlap"));

  			if (properties.getProperty(prefix+"plWeakWeight")!=null)      this.plWeakWeight=Double.parseDouble(properties.getProperty(prefix+"plWeakWeight"));
  			if (properties.getProperty(prefix+"plWeakEigen")!=null)       this.plWeakEigen=Double.parseDouble(properties.getProperty(prefix+"plWeakEigen"));
  			if (properties.getProperty(prefix+"plWeakWeight2")!=null)     this.plWeakWeight2=Double.parseDouble(properties.getProperty(prefix+"plWeakWeight2"));
  			if (properties.getProperty(prefix+"plWeakEigen2")!=null)      this.plWeakEigen2=Double.parseDouble(properties.getProperty(prefix+"plWeakEigen2"));
  			if (properties.getProperty(prefix+"plSumThick")!=null)        this.plSumThick=Double.parseDouble(properties.getProperty(prefix+"plSumThick"));
  			if (properties.getProperty(prefix+"plNeNeibCost")!=null)      this.plNeNeibCost=Double.parseDouble(properties.getProperty(prefix+"plNeNeibCost"));
  			if (properties.getProperty(prefix+"plNeOwn")!=null)           this.plNeOwn=Double.parseDouble(properties.getProperty(prefix+"plNeOwn"));

  			if (properties.getProperty(prefix+"plExNeibCost")!=null)      this.plExNeibCost=Double.parseDouble(properties.getProperty(prefix+"plExNeibCost"));
  			if (properties.getProperty(prefix+"plExNeibSmooth")!=null)    this.plExNeibSmooth=Double.parseDouble(properties.getProperty(prefix+"plExNeibSmooth"));
  			if (properties.getProperty(prefix+"plMergeCostStar")!=null)   this.plMergeCostStar=Double.parseDouble(properties.getProperty(prefix+"plMergeCostStar"));
  			if (properties.getProperty(prefix+"plMergeCost")!=null)       this.plMergeCost=Double.parseDouble(properties.getProperty(prefix+"plMergeCost"));

  			if (properties.getProperty(prefix+"plConflMerge")!=null)      this.plConflMerge=Boolean.parseBoolean(properties.getProperty(prefix+"plConflMerge"));
  			if (properties.getProperty(prefix+"plConflRelax")!=null)      this.plConflRelax=Double.parseDouble(properties.getProperty(prefix+"plConflRelax"));
  			if (properties.getProperty(prefix+"plConflSngl")!=null)       this.plConflSngl=Boolean.parseBoolean(properties.getProperty(prefix+"plConflSngl"));
  			if (properties.getProperty(prefix+"plConflSnglPair")!=null)   this.plConflSnglPair=Boolean.parseBoolean(properties.getProperty(prefix+"plConflSnglPair"));

  			if (properties.getProperty(prefix+"plWeakFgStrength")!=null)  this.plWeakFgStrength=Double.parseDouble(properties.getProperty(prefix+"plWeakFgStrength"));
  			if (properties.getProperty(prefix+"plWeakFgOutliers")!=null)  this.plWeakFgOutliers=Integer.parseInt(properties.getProperty(prefix+"plWeakFgOutliers"));
  			if (properties.getProperty(prefix+"plWeakFgRelax")!=null)     this.plWeakFgRelax=Double.parseDouble(properties.getProperty(prefix+"plWeakFgRelax"));
  			
  			if (properties.getProperty(prefix+"plThickWorld")!=null)      this.plThickWorld=Double.parseDouble(properties.getProperty(prefix+"plThickWorld"));
  			if (properties.getProperty(prefix+"plThickWorldConfl")!=null) this.plThickWorldConfl=Double.parseDouble(properties.getProperty(prefix+"plThickWorldConfl"));
  			if (properties.getProperty(prefix+"plRelaxComplete")!=null)   this.plRelaxComplete=Double.parseDouble(properties.getProperty(prefix+"plRelaxComplete"));
  			if (properties.getProperty(prefix+"plRelaxComplete2")!=null)  this.plRelaxComplete2=Double.parseDouble(properties.getProperty(prefix+"plRelaxComplete2"));
  			
  			if (properties.getProperty(prefix+"plMaxZRatio")!=null)       this.plMaxZRatio=Double.parseDouble(properties.getProperty(prefix+"plMaxZRatio"));
  			if (properties.getProperty(prefix+"plMaxDisp")!=null)         this.plMaxDisp=Double.parseDouble(properties.getProperty(prefix+"plMaxDisp"));
  			if (properties.getProperty(prefix+"plCutTail")!=null)         this.plCutTail=Double.parseDouble(properties.getProperty(prefix+"plCutTail"));
  			if (properties.getProperty(prefix+"plMinTail")!=null)         this.plMinTail=Double.parseDouble(properties.getProperty(prefix+"plMinTail"));
  			
  			if (properties.getProperty(prefix+"plDiscrEn")!=null)         this.plDiscrEn=Boolean.parseBoolean(properties.getProperty(prefix+"plDiscrEn"));
  			if (properties.getProperty(prefix+"plDiscrTolerance")!=null)  this.plDiscrTolerance=Double.parseDouble(properties.getProperty(prefix+"plDiscrTolerance"));
  			if (properties.getProperty(prefix+"plDiscrDispRange")!=null)  this.plDiscrDispRange=Double.parseDouble(properties.getProperty(prefix+"plDiscrDispRange"));
  			if (properties.getProperty(prefix+"plDiscrSteps")!=null)      this.plDiscrSteps=Integer.parseInt(properties.getProperty(prefix+"plDiscrSteps"));
  			if (properties.getProperty(prefix+"plDiscrVariants")!=null)   this.plDiscrVariants=Integer.parseInt(properties.getProperty(prefix+"plDiscrVariants"));
  			if (properties.getProperty(prefix+"plDiscrMode")!=null)       this.plDiscrMode=Integer.parseInt(properties.getProperty(prefix+"plDiscrMode"));
  			if (properties.getProperty(prefix+"plDiscrVarFloor")!=null)   this.plDiscrVarFloor=Double.parseDouble(properties.getProperty(prefix+"plDiscrVarFloor"));
  			if (properties.getProperty(prefix+"plDiscrSigma")!=null)      this.plDiscrSigma=Double.parseDouble(properties.getProperty(prefix+"plDiscrSigma"));
  			if (properties.getProperty(prefix+"plDiscrBlur")!=null)       this.plDiscrBlur=Double.parseDouble(properties.getProperty(prefix+"plDiscrBlur"));
  			if (properties.getProperty(prefix+"plDiscrExclusivity")!=null)this.plDiscrExclusivity=Double.parseDouble(properties.getProperty(prefix+"plDiscrExclusivity"));
  			if (properties.getProperty(prefix+"plDiscrExclus2")!=null)    this.plDiscrExclus2=Double.parseDouble(properties.getProperty(prefix+"plDiscrExclus2"));
  			if (properties.getProperty(prefix+"plDiscrStrict")!=null)     this.plDiscrStrict=Boolean.parseBoolean(properties.getProperty(prefix+"plDiscrStrict"));
  			if (properties.getProperty(prefix+"plDiscrCorrMax")!=null)    this.plDiscrCorrMax=Double.parseDouble(properties.getProperty(prefix+"plDiscrCorrMax"));
  			if (properties.getProperty(prefix+"plDiscrCorrMerge")!=null)  this.plDiscrCorrMerge=Double.parseDouble(properties.getProperty(prefix+"plDiscrCorrMerge"));
  			if (properties.getProperty(prefix+"plDiscrSteal")!=null)      this.plDiscrSteal=Integer.parseInt(properties.getProperty(prefix+"plDiscrSteal"));
  			if (properties.getProperty(prefix+"plDiscrGrown")!=null)      this.plDiscrGrown=Integer.parseInt(properties.getProperty(prefix+"plDiscrGrown"));
  			if (properties.getProperty(prefix+"plDiscrXMedian")!=null)    this.plDiscrXMedian=Double.parseDouble(properties.getProperty(prefix+"plDiscrXMedian"));

  			if (properties.getProperty(prefix+"plCostDist")!=null)        this.plCostDist=Double.parseDouble(properties.getProperty(prefix+"plCostDist"));
  			if (properties.getProperty(prefix+"plCostKrq")!=null)         this.plCostKrq=Double.parseDouble(properties.getProperty(prefix+"plCostKrq"));
  			if (properties.getProperty(prefix+"plCostKrqEq")!=null)       this.plCostKrqEq=Double.parseDouble(properties.getProperty(prefix+"plCostKrqEq"));
  			if (properties.getProperty(prefix+"plCostWrq")!=null)         this.plCostWrq=Double.parseDouble(properties.getProperty(prefix+"plCostWrq"));
  			if (properties.getProperty(prefix+"plCostWrqEq")!=null)       this.plCostWrqEq=Double.parseDouble(properties.getProperty(prefix+"plCostWrqEq"));
  			if (properties.getProperty(prefix+"plCostSin2")!=null)        this.plCostSin2=Double.parseDouble(properties.getProperty(prefix+"plCostSin2"));
  			if (properties.getProperty(prefix+"plCostRdist2")!=null)      this.plCostRdist2=Double.parseDouble(properties.getProperty(prefix+"plCostRdist2"));

  			if (properties.getProperty(prefix+"plConflDualTri")!=null)    this.plConflDualTri=Boolean.parseBoolean(properties.getProperty(prefix+"plConflDualTri"));
  			if (properties.getProperty(prefix+"plConflMulti")!=null)      this.plConflMulti=Boolean.parseBoolean(properties.getProperty(prefix+"plConflMulti"));
  			if (properties.getProperty(prefix+"plConflDiag")!=null)       this.plConflDiag=Boolean.parseBoolean(properties.getProperty(prefix+"plConflDiag"));
  			if (properties.getProperty(prefix+"plConflStar")!=null)       this.plConflStar=Boolean.parseBoolean(properties.getProperty(prefix+"plConflStar"));
  			if (properties.getProperty(prefix+"plStarSteps")!=null)       this.plStarSteps=Integer.parseInt(properties.getProperty(prefix+"plStarSteps"));
  			if (properties.getProperty(prefix+"plStarOrtho")!=null)       this.plStarOrtho=Double.parseDouble(properties.getProperty(prefix+"plStarOrtho"));
  			if (properties.getProperty(prefix+"plStarDiag")!=null)        this.plStarDiag=Double.parseDouble(properties.getProperty(prefix+"plStarDiag"));
  			if (properties.getProperty(prefix+"plStarPwr")!=null)         this.plStarPwr=Double.parseDouble(properties.getProperty(prefix+"plStarPwr"));
  			if (properties.getProperty(prefix+"plStarWeightPwr")!=null)   this.plStarWeightPwr=Double.parseDouble(properties.getProperty(prefix+"plStarWeightPwr"));
  			if (properties.getProperty(prefix+"plWeightToDens")!=null)    this.plWeightToDens=Double.parseDouble(properties.getProperty(prefix+"plWeightToDens"));
  			if (properties.getProperty(prefix+"plStarValPwr")!=null)      this.plStarValPwr=Double.parseDouble(properties.getProperty(prefix+"plStarValPwr"));
  			if (properties.getProperty(prefix+"plDblTriLoss")!=null)      this.plDblTriLoss=Double.parseDouble(properties.getProperty(prefix+"plDblTriLoss"));
  			if (properties.getProperty(prefix+"plNewConfl")!=null)        this.plNewConfl=Boolean.parseBoolean(properties.getProperty(prefix+"plNewConfl"));
  			if (properties.getProperty(prefix+"plMaxChanges")!=null)      this.plMaxChanges=Integer.parseInt(properties.getProperty(prefix+"plMaxChanges"));
  			
  			if (properties.getProperty(prefix+"plMutualOnly")!=null)      this.plMutualOnly=Boolean.parseBoolean(properties.getProperty(prefix+"plMutualOnly"));
  			if (properties.getProperty(prefix+"plFillSquares")!=null)     this.plFillSquares=Boolean.parseBoolean(properties.getProperty(prefix+"plFillSquares"));
  			if (properties.getProperty(prefix+"plCutCorners")!=null)      this.plCutCorners=Boolean.parseBoolean(properties.getProperty(prefix+"plCutCorners"));
  			if (properties.getProperty(prefix+"plHypotenuse")!=null)      this.plHypotenuse=Boolean.parseBoolean(properties.getProperty(prefix+"plHypotenuse"));

  			if (properties.getProperty(prefix+"plPull")!=null)            this.plPull=Double.parseDouble(properties.getProperty(prefix+"plPull"));
  			if (properties.getProperty(prefix+"plNormPow")!=null)         this.plNormPow=Double.parseDouble(properties.getProperty(prefix+"plNormPow"));
  			if (properties.getProperty(prefix+"plIterations")!=null)      this.plIterations=Integer.parseInt(properties.getProperty(prefix+"plIterations"));
  			if (properties.getProperty(prefix+"plStopBad")!=null)         this.plStopBad=Boolean.parseBoolean(properties.getProperty(prefix+"plStopBad"));
  			if (properties.getProperty(prefix+"plPrecision")!=null)       this.plPrecision=Integer.parseInt(properties.getProperty(prefix+"plPrecision"));

  			if (properties.getProperty(prefix+"plSplitPull")!=null)       this.plSplitPull=Double.parseDouble(properties.getProperty(prefix+"plSplitPull"));
  			if (properties.getProperty(prefix+"plSplitMinNeib")!=null)    this.plSplitMinNeib=Integer.parseInt(properties.getProperty(prefix+"plSplitMinNeib"));
  			if (properties.getProperty(prefix+"plSplitMinWeight")!=null)  this.plSplitMinWeight=Double.parseDouble(properties.getProperty(prefix+"plSplitMinWeight"));
  			if (properties.getProperty(prefix+"plSplitMinQuality")!=null) this.plSplitMinQuality=Double.parseDouble(properties.getProperty(prefix+"plSplitMinQuality"));
  			if (properties.getProperty(prefix+"plSplitApply")!=null)      this.plSplitApply=Boolean.parseBoolean(properties.getProperty(prefix+"plSplitApply"));
  			if (properties.getProperty(prefix+"plNonExclusive")!=null)    this.plNonExclusive=Boolean.parseBoolean(properties.getProperty(prefix+"plNonExclusive"));
  			if (properties.getProperty(prefix+"plUseOtherPlanes")!=null)  this.plUseOtherPlanes=Boolean.parseBoolean(properties.getProperty(prefix+"plUseOtherPlanes"));
  			if (properties.getProperty(prefix+"plAllowParallel")!=null)   this.plAllowParallel=Boolean.parseBoolean(properties.getProperty(prefix+"plAllowParallel"));
  			if (properties.getProperty(prefix+"plMaxDiff")!=null)         this.plMaxDiff=Double.parseDouble(properties.getProperty(prefix+"plMaxDiff"));
  			if (properties.getProperty(prefix+"plOtherDiff")!=null)       this.plOtherDiff=Double.parseDouble(properties.getProperty(prefix+"plOtherDiff"));
  			if (properties.getProperty(prefix+"plSplitXY")!=null)         this.plSplitXY=Boolean.parseBoolean(properties.getProperty(prefix+"plSplitXY"));
  			if (properties.getProperty(prefix+"plSplitXYTolerance")!=null)this.plSplitXYTolerance=Double.parseDouble(properties.getProperty(prefix+"plSplitXYTolerance"));

  			if (properties.getProperty(prefix+"plFuse")!=null)            this.plFuse=Boolean.parseBoolean(properties.getProperty(prefix+"plFuse"));
  			if (properties.getProperty(prefix+"plKeepOrphans")!=null)     this.plKeepOrphans=Boolean.parseBoolean(properties.getProperty(prefix+"plKeepOrphans"));
  			if (properties.getProperty(prefix+"plMinOrphan")!=null)       this.plMinOrphan=Double.parseDouble(properties.getProperty(prefix+"plMinOrphan"));

  			if (properties.getProperty(prefix+"plSnapDispAny")!=null)     this.plSnapDispAny=Double.parseDouble(properties.getProperty(prefix+"plSnapDispAny"));
  			if (properties.getProperty(prefix+"plSnapStrengthAny")!=null) this.plSnapStrengthAny=Double.parseDouble(properties.getProperty(prefix+"plSnapStrengthAny"));
  			if (properties.getProperty(prefix+"plSnapNegAny")!=null)      this.plSnapNegAny=Double.parseDouble(properties.getProperty(prefix+"plSnapNegAny"));
  			if (properties.getProperty(prefix+"plSnapDispMax")!=null)     this.plSnapDispMax=Double.parseDouble(properties.getProperty(prefix+"plSnapDispMax"));
  			if (properties.getProperty(prefix+"plSnapDispWeight")!=null)  this.plSnapDispWeight=Double.parseDouble(properties.getProperty(prefix+"plSnapDispWeight"));
  			if (properties.getProperty(prefix+"plSnapZeroMode")!=null)    this.plPrecision=Integer.parseInt(properties.getProperty(prefix+"plSnapZeroMode"));

  			if (properties.getProperty(prefix+"msUseSel")!=null)          this.msUseSel=Boolean.parseBoolean(properties.getProperty(prefix+"msUseSel"));
  			if (properties.getProperty(prefix+"msDivideByArea")!=null)    this.msDivideByArea=Boolean.parseBoolean(properties.getProperty(prefix+"msDivideByArea"));
  			if (properties.getProperty(prefix+"msScaleProj")!=null)       this.msScaleProj=Double.parseDouble(properties.getProperty(prefix+"msScaleProj"));
  			if (properties.getProperty(prefix+"msFractUni")!=null)        this.msFractUni=Double.parseDouble(properties.getProperty(prefix+"msFractUni"));


  			if (properties.getProperty(prefix+"tsNoEdge")!=null)          this.tsNoEdge=Boolean.parseBoolean(properties.getProperty(prefix+"tsNoEdge"));
  			if (properties.getProperty(prefix+"tsUseCenter")!=null)       this.tsUseCenter=Boolean.parseBoolean(properties.getProperty(prefix+"tsUseCenter"));
  			if (properties.getProperty(prefix+"tsMaxDiff")!=null)         this.tsMaxDiff=Double.parseDouble(properties.getProperty(prefix+"tsMaxDiff"));
  			if (properties.getProperty(prefix+"tsMinDiffOther")!=null)    this.tsMinDiffOther=Double.parseDouble(properties.getProperty(prefix+"tsMinDiffOther"));
  			if (properties.getProperty(prefix+"tsMinStrength")!=null)     this.tsMinStrength=Double.parseDouble(properties.getProperty(prefix+"tsMinStrength"));
  			if (properties.getProperty(prefix+"tsMaxStrength")!=null)     this.tsMaxStrength=Double.parseDouble(properties.getProperty(prefix+"tsMaxStrength"));
  			if (properties.getProperty(prefix+"tsMinSurface")!=null)      this.tsMinSurface=Double.parseDouble(properties.getProperty(prefix+"tsMinSurface"));
  			if (properties.getProperty(prefix+"tsMoveDirs")!=null)        this.tsMoveDirs=Integer.parseInt(properties.getProperty(prefix+"tsMoveDirs"));
  			if (properties.getProperty(prefix+"tsSurfStrPow")!=null)      this.tsSurfStrPow=Double.parseDouble(properties.getProperty(prefix+"tsSurfStrPow"));
  			if (properties.getProperty(prefix+"tsAddStrength")!=null)     this.tsAddStrength=Double.parseDouble(properties.getProperty(prefix+"tsAddStrength"));
  			if (properties.getProperty(prefix+"tsSigma")!=null)           this.tsSigma=Double.parseDouble(properties.getProperty(prefix+"tsSigma"));
  			if (properties.getProperty(prefix+"tsNSigma")!=null)          this.tsNSigma=Double.parseDouble(properties.getProperty(prefix+"tsNSigma"));
  			if (properties.getProperty(prefix+"tsMinPull")!=null)         this.tsMinPull=Double.parseDouble(properties.getProperty(prefix+"tsMinPull"));
  			if (properties.getProperty(prefix+"tsMinAdvantage")!=null)    this.tsMinAdvantage=Double.parseDouble(properties.getProperty(prefix+"tsMinAdvantage"));
  			
  			if (properties.getProperty(prefix+"tsClustSize")!=null)       this.tsClustSize=Integer.parseInt(properties.getProperty(prefix+"tsClustSize"));
  			if (properties.getProperty(prefix+"tsClustWeight")!=null)     this.tsClustWeight=Double.parseDouble(properties.getProperty(prefix+"tsClustWeight"));
  			if (properties.getProperty(prefix+"tsMinNeib")!=null)         this.tsMinNeib=Integer.parseInt(properties.getProperty(prefix+"tsMinNeib"));
  			if (properties.getProperty(prefix+"tsMaxSurStrength")!=null)  this.tsMaxSurStrength=Double.parseDouble(properties.getProperty(prefix+"tsMaxSurStrength"));
  			if (properties.getProperty(prefix+"tsCountDis")!=null)        this.tsCountDis=Boolean.parseBoolean(properties.getProperty(prefix+"tsCountDis"));
  			
  			if (properties.getProperty(prefix+"tsEnPlaneSeed")!=null)     this.tsEnPlaneSeed=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnPlaneSeed"));
  			if (properties.getProperty(prefix+"tsEnOnly")!=null)          this.tsEnOnly=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnOnly"));
  			if (properties.getProperty(prefix+"tsEnGrow")!=null)          this.tsEnGrow=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnGrow"));
  			if (properties.getProperty(prefix+"tsGrowStrength")!=null)    this.tsGrowStrength=Double.parseDouble(properties.getProperty(prefix+"tsGrowStrength"));
  			if (properties.getProperty(prefix+"tsGrowStrong")!=null)      this.tsGrowStrong=Boolean.parseBoolean(properties.getProperty(prefix+"tsGrowStrong"));
  			if (properties.getProperty(prefix+"tsContStrength")!=null)    this.tsGrowStrength=Double.parseDouble(properties.getProperty(prefix+"tsContStrength"));
  			if (properties.getProperty(prefix+"tsContDiff")!=null)        this.tsContDiff=Double.parseDouble(properties.getProperty(prefix+"tsContDiff"));

  			if (properties.getProperty(prefix+"tsEnSingle")!=null)        this.tsEnSingle=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnSingle"));
  			if (properties.getProperty(prefix+"tsEnMulti")!=null)         this.tsEnMulti=Boolean.parseBoolean(properties.getProperty(prefix+"tsEnMulti"));
  			if (properties.getProperty(prefix+"tsRemoveWeak1")!=null)     this.tsRemoveWeak1=Boolean.parseBoolean(properties.getProperty(prefix+"tsRemoveWeak1"));
  			if (properties.getProperty(prefix+"tsGrowSurround")!=null)    this.tsGrowSurround=Boolean.parseBoolean(properties.getProperty(prefix+"tsGrowSurround"));
  			if (properties.getProperty(prefix+"tsRemoveWeak2")!=null)     this.tsRemoveWeak2=Boolean.parseBoolean(properties.getProperty(prefix+"tsRemoveWeak2"));

  			if (properties.getProperty(prefix+"tsLoopMulti")!=null)       this.tsLoopMulti=Boolean.parseBoolean(properties.getProperty(prefix+"tsLoopMulti"));
  			if (properties.getProperty(prefix+"tsShow")!=null)            this.tsShow=Boolean.parseBoolean(properties.getProperty(prefix+"tsShow"));
  			if (properties.getProperty(prefix+"tsNumClust")!=null)        this.tsNumClust=Integer.parseInt(properties.getProperty(prefix+"tsNumClust"));

  			if (properties.getProperty(prefix+"tsConsensMode")!=null)     this.tsConsensMode=Integer.parseInt(properties.getProperty(prefix+"tsConsensMode"));
  			if (properties.getProperty(prefix+"tsConsensAgree")!=null)    this.tsConsensAgree=Integer.parseInt(properties.getProperty(prefix+"tsConsensAgree"));
  			
  			if (properties.getProperty(prefix+"taMinFgBg")!=null)         this.taMinFgBg=Double.parseDouble(properties.getProperty(prefix+"taMinFgBg"));
  			if (properties.getProperty(prefix+"taMinFgEdge")!=null)       this.taMinFgEdge=Double.parseDouble(properties.getProperty(prefix+"taMinFgEdge"));
  			if (properties.getProperty(prefix+"taMinColSep")!=null)       this.taMinColSep=Double.parseDouble(properties.getProperty(prefix+"taMinColSep"));
  			if (properties.getProperty(prefix+"taMinColDiff")!=null)      this.taMinColDiff=Double.parseDouble(properties.getProperty(prefix+"taMinColDiff"));
  			if (properties.getProperty(prefix+"taOutlier")!=null)         this.taOutlier=Double.parseDouble(properties.getProperty(prefix+"taOutlier"));
  			if (properties.getProperty(prefix+"taDiffPwr")!=null)         this.taDiffPwr=Double.parseDouble(properties.getProperty(prefix+"taDiffPwr"));
  			if (properties.getProperty(prefix+"taBestPwr")!=null)         this.taBestPwr=Double.parseDouble(properties.getProperty(prefix+"taBestPwr"));
  			if (properties.getProperty(prefix+"taDiff9Pwr")!=null)        this.taDiff9Pwr=Double.parseDouble(properties.getProperty(prefix+"taDiff9Pwr"));
  			if (properties.getProperty(prefix+"taColSigma")!=null)        this.taDiff9Pwr=Double.parseDouble(properties.getProperty(prefix+"taColSigma"));
  			if (properties.getProperty(prefix+"taColFraction")!=null)     this.taColFraction=Double.parseDouble(properties.getProperty(prefix+"taColFraction"));

  			if (properties.getProperty(prefix+"taCostEmpty")!=null)       this.taCostEmpty=Double.parseDouble(properties.getProperty(prefix+"taCostEmpty"));
  			if (properties.getProperty(prefix+"taCostNoLink")!=null)      this.taCostNoLink=Double.parseDouble(properties.getProperty(prefix+"taCostNoLink"));
  			if (properties.getProperty(prefix+"taCostSwitch")!=null)      this.taCostSwitch=Double.parseDouble(properties.getProperty(prefix+"taCostSwitch"));
  			if (properties.getProperty(prefix+"taCostColor")!=null)       this.taCostColor=Double.parseDouble(properties.getProperty(prefix+"taCostColor"));
  			if (properties.getProperty(prefix+"taCostDiff")!=null)        this.taCostDiff=Double.parseDouble(properties.getProperty(prefix+"taCostDiff"));
  			if (properties.getProperty(prefix+"taCostDiffBest")!=null)    this.taCostDiffBest=Double.parseDouble(properties.getProperty(prefix+"taCostDiffBest"));
  			if (properties.getProperty(prefix+"taCostDiff9")!=null)       this.taCostDiff9=Double.parseDouble(properties.getProperty(prefix+"taCostDiff9"));
  			if (properties.getProperty(prefix+"taCostWeakFgnd")!=null)    this.taCostWeakFgnd=Double.parseDouble(properties.getProperty(prefix+"taCostWeakFgnd"));
  			if (properties.getProperty(prefix+"taCostFlaps")!=null)       this.taCostFlaps=Double.parseDouble(properties.getProperty(prefix+"taCostFlaps"));
  			if (properties.getProperty(prefix+"taCostMismatch")!=null)    this.taCostMismatch=Double.parseDouble(properties.getProperty(prefix+"taCostMismatch"));

  			if (properties.getProperty(prefix+"taEnEmpty")!=null)         this.taEnEmpty=Boolean.parseBoolean(properties.getProperty(prefix+"taEnEmpty"));
  			if (properties.getProperty(prefix+"taEnNoLink")!=null)        this.taEnNoLink=Boolean.parseBoolean(properties.getProperty(prefix+"taEnNoLink"));
  			if (properties.getProperty(prefix+"taEnSwitch")!=null)        this.taEnSwitch=Boolean.parseBoolean(properties.getProperty(prefix+"taEnSwitch"));
  			if (properties.getProperty(prefix+"taEnColor")!=null)         this.taEnColor=Boolean.parseBoolean(properties.getProperty(prefix+"taEnColor"));
  			if (properties.getProperty(prefix+"taEnDiff")!=null)          this.taEnDiff=Boolean.parseBoolean(properties.getProperty(prefix+"taEnDiff"));
  			if (properties.getProperty(prefix+"taEnDiffBest")!=null)      this.taEnDiffBest=Boolean.parseBoolean(properties.getProperty(prefix+"taEnDiffBest"));
  			if (properties.getProperty(prefix+"taEnDiff9")!=null)         this.taEnDiff9=Boolean.parseBoolean(properties.getProperty(prefix+"taEnDiff9"));
  			if (properties.getProperty(prefix+"taEnWeakFgnd")!=null)      this.taEnWeakFgnd=Boolean.parseBoolean(properties.getProperty(prefix+"taEnWeakFgnd"));
  			if (properties.getProperty(prefix+"taEnFlaps")!=null)         this.taEnFlaps=Boolean.parseBoolean(properties.getProperty(prefix+"taEnFlaps"));
  			if (properties.getProperty(prefix+"taEnMismatch")!=null)      this.taEnMismatch=Boolean.parseBoolean(properties.getProperty(prefix+"taEnMismatch"));
  			
  			
  			if (properties.getProperty(prefix+"dbg_migrate")!=null)       this.dbg_migrate=Boolean.parseBoolean(properties.getProperty(prefix+"dbg_migrate"));

  			if (properties.getProperty(prefix+"show_ortho_combine")!=null)     this.show_ortho_combine=Boolean.parseBoolean(properties.getProperty(prefix+"show_ortho_combine"));
  			if (properties.getProperty(prefix+"show_refine_supertiles")!=null) this.show_refine_supertiles=Boolean.parseBoolean(properties.getProperty(prefix+"show_refine_supertiles"));
  			if (properties.getProperty(prefix+"show_bgnd_nonbgnd")!=null)      this.show_bgnd_nonbgnd=Boolean.parseBoolean(properties.getProperty(prefix+"show_bgnd_nonbgnd"));
  			if (properties.getProperty(prefix+"show_filter_scan")!=null)       this.show_filter_scan=Boolean.parseBoolean(properties.getProperty(prefix+"show_filter_scan"));
  			if (properties.getProperty(prefix+"show_combined")!=null)          this.show_combined=Boolean.parseBoolean(properties.getProperty(prefix+"show_combined"));
  			if (properties.getProperty(prefix+"show_unique")!=null)            this.show_unique=Boolean.parseBoolean(properties.getProperty(prefix+"show_unique"));
  			if (properties.getProperty(prefix+"show_histograms")!=null)        this.show_histograms=Boolean.parseBoolean(properties.getProperty(prefix+"show_histograms"));
  			if (properties.getProperty(prefix+"show_init_refine")!=null)       this.show_init_refine=Boolean.parseBoolean(properties.getProperty(prefix+"show_init_refine"));
  			if (properties.getProperty(prefix+"show_expand")!=null)            this.show_expand=Boolean.parseBoolean(properties.getProperty(prefix+"show_expand"));
  			if (properties.getProperty(prefix+"show_variant")!=null)           this.show_variant=Boolean.parseBoolean(properties.getProperty(prefix+"show_variant"));
  			if (properties.getProperty(prefix+"show_retry_far")!=null)         this.show_retry_far=Boolean.parseBoolean(properties.getProperty(prefix+"show_retry_far"));
  			if (properties.getProperty(prefix+"show_macro")!=null)             this.show_macro=Boolean.parseBoolean(properties.getProperty(prefix+"show_macro"));
  			if (properties.getProperty(prefix+"show_shells")!=null)            this.show_shells=Boolean.parseBoolean(properties.getProperty(prefix+"show_shells"));
  			if (properties.getProperty(prefix+"show_neighbors")!=null)         this.show_neighbors=Boolean.parseBoolean(properties.getProperty(prefix+"show_neighbors"));
  			if (properties.getProperty(prefix+"show_flaps_dirs")!=null)        this.show_flaps_dirs=Boolean.parseBoolean(properties.getProperty(prefix+"show_flaps_dirs"));
  			if (properties.getProperty(prefix+"show_first_clusters")!=null)    this.show_first_clusters=Boolean.parseBoolean(properties.getProperty(prefix+"show_first_clusters"));
  			if (properties.getProperty(prefix+"show_planes")!=null)            this.show_planes=Boolean.parseBoolean(properties.getProperty(prefix+"show_planes"));
  			
  			
  			if (properties.getProperty(prefix+"vertical_xyz.x")!=null)     this.vertical_xyz[0]=Double.parseDouble(properties.getProperty(prefix+"vertical_xyz.x"));
  			if (properties.getProperty(prefix+"vertical_xyz.y")!=null)     this.vertical_xyz[1]=Double.parseDouble(properties.getProperty(prefix+"vertical_xyz.y"));
  			if (properties.getProperty(prefix+"vertical_xyz.z")!=null)     this.vertical_xyz[2]=Double.parseDouble(properties.getProperty(prefix+"vertical_xyz.z"));
  			
  		}
  		
  		public boolean showDialog() {
  			GenericDialog gd = new GenericDialog("Set CLT parameters");
  			gd.addNumericField("Transform size (default 8)",                                              this.transform_size,            0);
  			gd.addNumericField("Lapped transform window type (0- rectangular, 1 - sinus)",                this.clt_window,                0);
   			gd.addNumericField("shift_x",                                                                 this.shift_x,                   4);
   			gd.addNumericField("shift_y",                                                                 this.shift_y,                   4);
  			gd.addNumericField("Bit mask - which of 4 transforms to combine after iclt",                  this.iclt_mask,                 0);
  			gd.addNumericField("Tile X to extract (0..163)",                                              this.tileX,                     0);
  			gd.addNumericField("Tile Y to extract (0..122)",                                              this.tileY,                     0);
  			gd.addNumericField("dbg_mode: 0 - normal, +1 - no DCT/IDCT, just fold",                       this.dbg_mode,                  0);
  			gd.addNumericField("ishift_x: shift source image by this pixels left",                        this.ishift_x,                  0);
  			gd.addNumericField("ishift_y: shift source image by this pixels down",                        this.ishift_y,                  0);
   			gd.addNumericField("Modify phase correlation to prevent division by very small numbers",      this.fat_zero,                  4);
   			gd.addNumericField("LPF correlarion sigma ",                                                  this.corr_sigma,                3);
  			gd.addCheckbox    ("Normalize kernels ",                                                      this.norm_kern);
  			gd.addCheckbox    ("Equalize green channel gain of the individual cnannels",                  this.gain_equalize);
  			gd.addCheckbox    ("Equalize R/G, B/G balance of the individual channels",                    this.colors_equalize);
  			gd.addNumericField("Reg gain in the center of sensor calibration R (instead of vignetting)",  this.novignetting_r,   4);
  			gd.addNumericField("Green gain in the center of sensor calibration G (instead of vignetting)",this.novignetting_g, 4);
  			gd.addNumericField("Blue gain in the center of sensor calibration B (instead of vignetting)", this.novignetting_b,  4);
  			gd.addNumericField("Extra red correction to compensate for light temperature",                this.scale_r,  4);
  			gd.addNumericField("Extra green correction to compensate for light temperature",              this.scale_g,  4);
  			gd.addNumericField("Extra blue correction to compensate for light temperature",               this.scale_b,  4);
  			gd.addNumericField("Value (max) in vignetting data to correspond to 1x in the kernel",        this.vignetting_max,      3);
  			gd.addNumericField("Do not try to correct vignetting smaller than this fraction of max",      this.vignetting_range,  3);
  			gd.addNumericField("Kernel step in pixels (has 1 kernel margin on each side)",                this.kernel_step,            0);
  			gd.addNumericField("Nominal (rectilinear) disparity between side of square cameras (pix)",    this.disparity,  3);
  			gd.addNumericField("Inverse distance to infinity (misalignment cortrection)",                 this.z_correction,  6);
  			gd.addCheckbox    ("Perform correlation",                                                     this.correlate);
  			gd.addNumericField("Bitmask of pairs to combine in the composite (top, bottom, left,righth)", this.corr_mask,            0);
  			gd.addCheckbox    ("Combine correlation with mirrored around disparity direction",            this.corr_sym);
  			gd.addCheckbox    ("Keep all partial correlations (otherwise - only combined one)",           this.corr_keep);
  			gd.addCheckbox    ("Show combined correlations",                                              this.corr_show);
  			gd.addCheckbox    ("Calculate per-pair X/Y variations of measured correlations ",             this.corr_mismatch);
  			gd.addNumericField("Add to pair correlation before multiplying by other pairs (between sum and product)",    this.corr_offset,  6);
  			gd.addNumericField("Red to green correlation weight",                                         this.corr_red,  4);
  			gd.addNumericField("Blue to green correlation weight",                                        this.corr_blue,  4);
  			gd.addCheckbox    ("Normalize each correlation tile by rms",                                  this.corr_normalize);
  			gd.addNumericField("Minimal correlation value to consider valid",                             this.min_corr,  6);
  			gd.addNumericField("Minimal correlation value to consider valid when normalizing results",    this.min_corr_normalized,  6);
  			gd.addNumericField("Sigma for weights of points around global max to find fractional",        this.max_corr_sigma,  3);
  			gd.addNumericField("Maximal distance from int max to consider",                               this.max_corr_radius,  3);
  			
  			
  			gd.addMessage("--- Enhance detection of horizontal/vertical features (when enh_ortho is enabled for tile ---");
  			gd.addNumericField("Reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)",  this.enhortho_width,            0);
  			gd.addNumericField("Multiply center correlation pixels (inside enhortho_width) (1.0 - disables enh_orttho)",  this.enhortho_scale,  3);
  			
  			
  			gd.addCheckbox    ("Double pass when masking center of mass to reduce preference for integer values", this.max_corr_double);
  			gd.addNumericField("Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial",   this.corr_mode,            0);
  			gd.addNumericField("Contrast of dotted border on correlation results",                        this.corr_border_contrast,  6);
  			gd.addMessage("--- tiles tasks (current tile_task_op = "+this.tile_task_op+") ---");
  			gd.addCheckbox    ("Enhace ortho lines detection (enh_ortho)",                                               ImageDtt.getOrthoLines(this.tile_task_op));
  			gd.addCheckbox    ("Force disparity for image rendering (false - use found from tile correlation)",          ImageDtt.getForcedDisparity(this.tile_task_op));
  			gd.addNumericField("Bitmask of used images (1 - top left, 2 - top right, 4 - bottom left, 8 bottom right)",  ImageDtt.getImgMask(this.tile_task_op),            0);
  			gd.addNumericField("Bitmask of used pairs  (1 - top, 2 - bottom, 4 - left, 8 - right)",                      ImageDtt.getPairMask(this.tile_task_op),            0);
//  			gd.addNumericField("Tile operations bits: +(0..f) - images, +(00.f0) - process pairs +256, ... ",  this.tile_task_op,            0);
  			gd.addNumericField("Tile operations window left (in 8x8 tiles)",                              this.tile_task_wl,            0);
  			gd.addNumericField("Tile operations window top",                                              this.tile_task_wt,            0);
  			gd.addNumericField("Tile operations window width",                                            this.tile_task_ww,            0);
  			gd.addNumericField("Tile operations window height",                                           this.tile_task_wh,            0);
  			
  			gd.addNumericField("Do not adjust for shot noise (~sqrt) if lower than this",                 this.min_shot,  4);
  			gd.addNumericField("Scale when dividing by sqrt for shot noise compensation of pixel differences (<0 - disable)", this.scale_shot,  4);
  			
  			gd.addNumericField("RMS difference from average to reduce weights (255 full scale image)",    this.diff_sigma,  4);
  			gd.addNumericField("RMS difference from average in sigmas to discard channel",                this.diff_threshold,  4);
  			gd.addCheckbox    ("Gaussian as weight when averaging images (false - sharp all/nothing)",    this.diff_gauss);
  			gd.addNumericField("Minimal number of channels to agree on a point (real number to work with fuzzy averages)",   this.min_agree,  2);
  			gd.addCheckbox    ("Do not reduce average weight when only one image differes much from the average", this.dust_remove);
  			gd.addCheckbox    ("Use black for backdrop outside of the FOV",                               this.black_back);
  			gd.addCheckbox    ("Add port weights to RGBA stack (debug feature)",                          this.keep_weights);
  			gd.addCheckbox    ("Alpha channel: use center 8x8 (unchecked - treat same as RGB)",           this.sharp_alpha);
  			gd.addNumericField("Alpha channel 0.0 thereshold (lower - transparent)",                      this.alpha0,   3);
  			gd.addNumericField("Alpha channel 1.0 threshold (higher - opaque)",                           this.alpha1,   3);
  			gd.addCheckbox    ("Generate shifted channel linear RGB stacks",                              this.gen_chn_stacks);
  			gd.addCheckbox    ("Generate shifted channel color image stack",                              this.gen_chn_img);
  			gd.addCheckbox    ("Generate shifted channel images and save with the model 'CLT process corr'",this.gen_4_img);
  			gd.addCheckbox    ("Show result RGBA before overlap combined",                                this.show_nonoverlap);
  			gd.addCheckbox    ("Show result RGBA",                                                        this.show_overlap);
  			gd.addCheckbox    ("Show result color",                                                       this.show_rgba_color);
  			gd.addCheckbox    ("Show disparity maps",                                                     this.show_map);
  			gd.addNumericField("Disparity scan start value",                                              this.disp_scan_start,  2);
  			gd.addNumericField("Disparity scan step",                                                     this.disp_scan_step,  2);
  			gd.addNumericField("Disparity scan number of disparity values to scan",                       this.disp_scan_count,            0);

  			gd.addMessage("--- camera fine correction: X/Y for images 0..3  ---");
  			gd.addCheckbox    ("Debug infinity/lazy eye correction",                                      this.fine_dbg);
  			gd.addNumericField("X 0",                                                                     this.fine_corr_x_0,  3);
  			gd.addNumericField("Y 0",                                                                     this.fine_corr_y_0,  3);
  			gd.addNumericField("X 1",                                                                     this.fine_corr_x_1,  3);
  			gd.addNumericField("Y 1",                                                                     this.fine_corr_y_1,  3);
  			gd.addNumericField("X 2",                                                                     this.fine_corr_x_2,  3);
  			gd.addNumericField("Y 2",                                                                     this.fine_corr_y_2,  3);
  			gd.addNumericField("X 3",                                                                     this.fine_corr_x_3,  3);
  			gd.addNumericField("Y 4",                                                                     this.fine_corr_y_3,  3);

  			gd.addNumericField("Y 4",                                                                     this.fcorr_radius,  3);
  			gd.addNumericField("Do not try to correct outside this fraction of width/hight",              this.fcorr_min_strength,3);
  			gd.addNumericField("Consider only tiles with absolute residual disparity lower than",         this.fcorr_disp_diff,  3);
  			gd.addCheckbox    ("Use quadratic polynomial for fine correction (false - only linear)",      this.fcorr_quadratic);
  			gd.addCheckbox    ("Ignore current calculated fine correction (use manual only)",             this.fcorr_ignore);

  			gd.addNumericField("Minimal correlation strength to use for infinity correction",             this.fcorr_inf_strength,3);
  			gd.addNumericField("Disparity half-range for infinity",                                       this.fcorr_inf_diff,  3);
  			gd.addCheckbox    ("Use quadratic polynomial for infinity correction (false - only linear)",  this.fcorr_inf_quad);
  			gd.addCheckbox    ("Correct infinity in vertical direction (false - only horizontal)",        this.fcorr_inf_vert);
  			
  			
  			gd.addCheckbox    ("Apply disparity correction to zero at infinity",                          this.inf_disp_apply);
  			gd.addNumericField("Re run disparity correction at infinity multiple times",                  this.inf_repeat,      0);

  			
  			//  			gd.addCheckbox    ("Apply lazy eye correction at infinity",                                   this.inf_mism_apply);
  			gd.addNumericField("Infinity extraction - maximum iterations",                                this.inf_iters,       0);
  			gd.addNumericField("Coefficients maximal increment to exit iterations",                       this.inf_final_diff,  6);
  			gd.addNumericField("Include farther tiles than tolerance, but scale their weights",           this.inf_far_pull,    3);
  			gd.addMessage     ("--- Infinity filter ---");
  			gd.addNumericField("Strength power",                                                          this.inf_str_pow,     3);
  			gd.addNumericField("Sample size (side of a square)",                                          this.inf_smpl_side,   0);
  			gd.addNumericField("Number after removing worst (should be >1)",                              this.inf_smpl_num,    0);
  			gd.addNumericField("Maximal RMS of the remaining tiles in a sample",                          this.inf_smpl_rms,    3);
  			gd.addMessage     ("--- Infinity histogram filter ---");
  			gd.addNumericField("Square sample step (50% overlap)",                                        this.ih_smpl_step,    0);
  			gd.addNumericField("Histogram minimal disparity",                                             this.ih_disp_min,     3);
  			gd.addNumericField("Histogram disparity step",                                                this.ih_disp_step,    3);
  			gd.addNumericField("Histogram number of bins",                                                this.ih_num_bins,     0);
  			gd.addNumericField("Histogram Gaussian sigma (in disparity pixels)",                          this.ih_sigma,        3);
  			gd.addNumericField("Keep samples within this difference from farthest maximum",               this.ih_max_diff,     3);
  			gd.addNumericField("Minimal number of remaining samples",                                     this.ih_min_samples,  0);
  			gd.addCheckbox    ("Replace samples with a single average with equal weight",                 this.ih_norm_center);

  			gd.addMessage     ("--- Lazy eye parameters (disparity @ infinity should be adjusted first ---");
			gd.addNumericField("Sample size (side of a square)",                                          this.ly_smpl_side,  0);
			gd.addNumericField("Number after removing worst (should be >1)",                              this.ly_smpl_num,  0);
			gd.addNumericField("Maximal measured relative disparity",                                     this.ly_meas_disp,  3);
			gd.addNumericField("Maximal RMS of the remaining tiles in a sample",                          this.ly_smpl_rms,  3);
			gd.addNumericField("Maximal full disparity difference to 8 neighbors",                        this.ly_disp_var,  3);
			gd.addNumericField("Relative weight of infinity calibration data",                            this.ly_inf_frac,  3);
  			gd.addCheckbox    ("Calculate and apply lazy eye correction after disparity scan (poly or extrinsic), may repeat",this.ly_on_scan);
  			gd.addCheckbox    ("Use infinity disparity (disable if there is not enough of infinity data), both poly and extrinsic", this.ly_inf_en);
  			gd.addCheckbox    ("Force convergence correction during extrinsic, even with no infinity data", this.ly_inf_force);
  			gd.addCheckbox    ("*Use polynomial correction, false - correct tilt/azimuth/roll of each sensor)", this.ly_poly);
  			gd.addMessage     ("---");
//  			gd.addNumericField("Use square this size side to detect outliers",                            this.fcorr_sample_size,  0);
//  			gd.addNumericField("Keep tiles only if there are more in each square",                        this.fcorr_mintiles,     0);
//  			gd.addNumericField("Remove this fraction of tiles from each sample",                          this.fcorr_reloutliers,  3);
//  			gd.addNumericField("Gaussian blur channel mismatch data",                                     this.fcorr_sigma,        3);

  			gd.addNumericField("Calculated from correlation offset vs. actual one (not yet understood)",  this.corr_magic_scale,  3);
  			
  			gd.addMessage     ("--- 3D reconstruction ---");
  			gd.addCheckbox    ("Show generated textures",                                                      this.show_textures);
  			gd.addCheckbox    ("show intermediate results of filtering",                                       this.debug_filters);

  			gd.addNumericField("Minimal noise-normalized pixel difference in a channel to suspect something",  this.min_smth,  3);
  			gd.addNumericField("Reliable noise-normalized pixel difference in a channel to have something ",   this.sure_smth,  3);
  			gd.addNumericField("Disparity range to be considered background",                                  this.bgnd_range,  3);
  			gd.addNumericField("Disparity difference from the center (provided) disparity to trust",           this.other_range,  3);

  			gd.addNumericField("Minimal 4-corr strength to trust tile",                                        this.ex_strength,  3);
  			gd.addNumericField("Minimal 4-corr strength divided by channel diff for new (border) tiles",       this.ex_nstrength,  3);
  			
  			gd.addCheckbox    ("Allow expansion over previously identified background (infinity)",             this.ex_over_bgnd);
  			gd.addNumericField("When expanding over background, disregard lower disparity ",                   this.ex_min_over,  3);

  			gd.addMessage     ("********* Plates filetering when building initial z-map *********");
  			gd.addNumericField("If strength exceeds ex_strength * super_trust, do not apply ex_nstrength and plate_ds", this.pt_super_trust,  3);
  			gd.addCheckbox    ("Do not replace raw tiles by the plates, if raw is closer (like poles)",        this.pt_keep_raw_fg);
  			gd.addNumericField("Scale plates strength before comparing to raw strength",                       this.pt_scale_pre,  3);
  			gd.addNumericField("Scale plates strength when replacing raw (plates d/s data is more reliable if it exists)", this.pt_scale_post,  3);

  			gd.addNumericField("Minimal strength to be considered definitely background",                      this.bgnd_sure,  3);
  			gd.addNumericField("Maximal strength to ignore as non-background",                                 this.bgnd_maybe,  3);

  			gd.addNumericField("Number of tiles in a cluster to seed (just background?)",                      this.min_clstr_seed,   0);
  			gd.addNumericField("Number of tiles in a cluster not close to other clusters (more than 2 tiles apart)", this.min_clstr_lone,   0);
  			gd.addNumericField("Minimal total strength of the cluster",                                        this.min_clstr_weight,  3);
  			gd.addNumericField("Minimal maximal strength of the cluster",                                      this.min_clstr_max,  3);
  			
  			gd.addNumericField("Fill gaps betsween clusters, see comments for 'grow'",                         this.fill_gaps,   0);
  			gd.addNumericField("Same as fill_gaps above, on the final pass",                                   this.fill_final,   0);
  			gd.addNumericField("Number of tiles in a cluster to block (just non-background?)",                 this.min_clstr_block,   0);
  			gd.addNumericField("Number of tiles to grow tile selection (1 - hor/vert, 2 - hor/vert/diagonal)", this.bgnd_grow,   0);

  			gd.addCheckbox    ("Use old ortho features processing (orth0_* parameters, false - use or_*)",     this.ortho_old);

  			gd.addNumericField("Minimal strength of hor correlation to be used instead of full 4-pair correlation", this.ortho_min_hor,  3);
  			gd.addNumericField("Minimal strength of vert correlation to be used instead of full 4-pair correlation", this.ortho_min_vert,  3);
  			gd.addNumericField("Vert/hor (or hor/vert) strength to be used instead of the full correlation",   this.ortho_asym,  3);
  			gd.addNumericField("Vert/hor (or hor/vert) strength exceeding scaled 4-pair strength",             this.ortho_over4,  3);

  			gd.addNumericField("Minimal strength of hor/vert to bridge over",                                  this.ortho_sustain,  3);
  			gd.addNumericField("minimal run of hor/vert tiles to be considered (at least from one side)",      this.ortho_run,      0);
  			gd.addNumericField("Minimal maximal strength in an ortho run (max. in run >=...)",                 this.ortho_minmax,  3);
  			gd.addNumericField("Number of tiles to bridge over hor/vert gaps",                                 this.ortho_bridge,   0);
  			gd.addNumericField("Maximal disparity RMS in a run to replace by average)",                        this.ortho_rms,  3);
  			gd.addNumericField("Convolve hor/vert strength by 3*(2*l+1) kernels to detect multi-tile features",this.ortho_half_length,   0);
  			gd.addNumericField("Fraction of convolved ortho in a mix with raw",                                this.ortho_mix,  3);
  			
  			gd.addMessage     ("--- Combination of ortho and 4-pair correlations ---");
  			gd.addCheckbox    ("Apply ortho correction to horizontal correlation (vertical features)",            this.or_hor);
  			gd.addCheckbox    ("Apply ortho correction to vertical correlation (horizontal features)",            this.or_vert);
  			gd.addNumericField("Blur sigma: verically for horizontal correlation, horizontally - for vertically", this.or_sigma,  3);
  			gd.addNumericField("3-point sharpening (-k, +2k+1, -k)",                                              this.or_sharp,  3);
  			gd.addNumericField("Scale ortho correletion strength relative to 4-directional one",                  this.or_scale,  3);
  			gd.addNumericField("Subtract from scaled correlation strength, limit by 0",                           this.or_offset,  3);
  			gd.addNumericField("Minimal ratio of orthogonal strengths required for dis[parity replacement",       this.or_asym,  3);
  			gd.addNumericField("Minimal scaled offset ortho strength to normal strength needed for replacement",  this.or_threshold,  3);
  			gd.addNumericField("Minimal horizontal absolute scaled offset ortho strength needed for replacement", this.or_absHor,  3);
  			gd.addNumericField("Minimal vertical absolute scaled offset ortho strength needed for replacement",   this.or_absVert,  3);

  			gd.addMessage     ("--- Fix vertical structures, such as street poles ---");
  			gd.addCheckbox    ("Continue vertical structures to the ground",                                      this.poles_fix);
  			gd.addNumericField("Number of tiles to extend over the poles bottoms",                                this.poles_len,   0);
  			gd.addNumericField("Maximal ratio of invisible to visible pole length",                               this.poles_ratio,  3);
  			gd.addNumericField("Set new pole segment strength to max of horizontal correlation and this value",   this.poles_min_strength,  3);
  			gd.addCheckbox    ("Set disparity to that of the bottom of existing segment (false - use hor. disparity)",this.poles_force_disp);
  			
  			gd.addNumericField("Maximal number of clusters to generate for one run",                           this.max_clusters,   0);
  			gd.addCheckbox    ("Remove all unneeded scans when generating x3d output to save memory",          this.remove_scans);
  			gd.addCheckbox    ("Generate x3d output",                                                          this.output_x3d);
  			gd.addCheckbox    ("Generate Wavefront obj output",                                                this.output_obj);
  			gd.addCheckbox    ("Correct lens geometric distortions in a model (will need backdrop to be corrected too)", this.correct_distortions);
  			gd.addCheckbox    ("Show generated triangles",                                                     this.show_triangles);
  			gd.addCheckbox    ("Weight-average disparity for the whole cluster ",                              this.avg_cluster_disp);
  			gd.addNumericField("Maximal disparity difference in a triangle face to show",                      this.maxDispTriangle,  6);
  			gd.addNumericField("Distance to generate backdrop (0 - use regular backdrop)",                     this.infinityDistance,  8);

  			gd.addCheckbox    ("Split into shells with flaps",                                                 this.shUseFlaps);
  			gd.addCheckbox    ("Aggressive fade alpha (whole boundary)",                                       this.shAggrFade);
  			gd.addNumericField("Minimal shell area (not counting flaps",                                       this.shMinArea,   0);
  			gd.addNumericField("Minimal value of the shell maximum strength",                                  this.shMinStrength,  6);

  			gd.addMessage     ("--- Thin ice parameters ---");
  			gd.addNumericField("Relative disparity rigidity in vertical direction",                            this.tiRigidVertical,  6);
  			gd.addNumericField("Relative disparity rigidity in horizontal direction",                          this.tiRigidHorizontal,  6);
  			gd.addNumericField("Relative disparity rigidity in diagonal   direction",                          this.tiRigidDiagonal,  6);
  			gd.addNumericField("Strength floor - subtract (limit with 0) before applying",                     this.tiStrengthOffset,  6);
  			gd.addNumericField("Divide actual disparity by this  before Math.pow and applying to TI",          this.tiDispScale,  6);
  			gd.addNumericField("Apply pow to disparity (restore sign) for disparity difference pressure on TI",this.tiDispPow,  6);
  			gd.addNumericField("tiDispPull: multiply strength*disparity difference to pull force",             this.tiDispPull,  6);
  			gd.addNumericField("Scale tiDispPull for pre-final pass",                                          this.tiDispPullPreFinal,  6);
  			gd.addNumericField("Scale tiDispPull for final pass",                                              this.tiDispPullFinal,  6);
  			gd.addNumericField("Normalize stresses to average disparity if it is above threshold",             this.tiBreakNorm,  6);
  			gd.addNumericField("TI break value of abs(d0-3d1+3d2-d3)",                                         this.tiBreak3,  6);
  			gd.addNumericField("TI break value of (d0-3d1+3d2-d3) * (d1 - d2)",                                this.tiBreak31,  6);
  			gd.addNumericField("TI break value of (-d0+d1+d2-d3) * abs(d1 - d2)",                              this.tiBreak21,  6);
  			gd.addNumericField("TI disparity threshold to remove as too far tiles",                            this.tiBreakFar,  6);
  			gd.addNumericField("TI disparity threshold to remove as too near tiles",                           this.tiBreakNear,  6);
  			gd.addNumericField("TI break mode: +1: abs(3-rd derivative), +2: -(3-rd * 1-st), +4: -(2-nd * abs(1-st))", this.tiBreakMode,  0);
  			gd.addNumericField("Amplify colinear breaks in neighbor tiles",                                    this.tiBreakSame,  6);
  			gd.addNumericField("Amplify 90-degree turnintg breaks in neighbor tiles",                          this.tiBreakTurn,  6);

  			gd.addNumericField("Heal disparity gap before pre-last smooth",                                    this.tiHealPreLast,  6);
  			gd.addNumericField("Heal disparity gap before last smooth",                                        this.tiHealLast,  6);
  			gd.addNumericField("Maximal length of an internal break in the cluster to heal",                   this.tiHealSame, 0);

  			gd.addNumericField("Maximal number of TI iterations for each step",                                this.tiIterations, 0);
  			gd.addNumericField("Iteration maximal error (1/power of 10)",                                      this.tiPrecision,  0);
  			gd.addNumericField("Number of cycles break-smooth (after the first smooth)",                       this.tiNumCycles,  0);
  			gd.addMessage     ("--- Fg/Bg separation ---");
  			gd.addCheckbox    ("Apply super-tiles during refine passes",                                       this.stUseRefine);
  			gd.addCheckbox    ("Apply super-tiles during pass2 ",                                              this.stUsePass2);
  			gd.addCheckbox    ("Apply super-tiles during render",                                              this.stUseRender);
  			
  			gd.addCheckbox    ("Show supertiles histograms",                                                   this.stShow);
  			gd.addNumericField("Super tile size (square, in tiles)",                                           this.stSize,  0);
  			gd.addNumericField("Disparity histogram step for far objects",                                     this.stStepFar,  6);
  			gd.addNumericField("Disparity histogram step for near objects",                                    this.stStepNear,  6);
  			gd.addNumericField("Disparity threshold to switch from linear to logarithmic steps",               this.stStepThreshold,  6);
  			gd.addNumericField("Minimal disparity (center of a bin)",                                          this.stMinDisparity,  6);
  			
//  			gd.addNumericField("Maximal disparity (center of a bin)",                                          this.stMaxDisparity,  6);
  			gd.addMessage     ("aximal disparity (center of a bin) - using grow_disp_max="+this.grow_disp_max);
  			gd.addNumericField("Subtract from strength, discard negative",                                     this.stFloor,  6);
  			gd.addNumericField("Raise strength to this power ",                                                this.stPow,  6);
  			gd.addNumericField("Blur disparity histogram (sigma in bins)",                                     this.stSigma,  6);
  			gd.addNumericField("Minimal backgroubnd disparity to extract as a maximum from the supertiles",    this.stMinBgDisparity,  6);
  			gd.addNumericField("Minimal fraction of the disparity histogram to use as background",             this.stMinBgFract,  6);
  			gd.addNumericField("Use background disparity from supertiles if tile strength is less",            this.stUseDisp,  6);
  			gd.addNumericField("Multiply st strength if used instead of regular strength ",                    this.stStrengthScale,  6);

  			gd.addCheckbox    ("Use sample mode (false - regular tile mode)",                                  this.stSmplMode);
  			gd.addNumericField("Sample size (side of a square)",                                               this.stSmplSide,  0);
  			gd.addNumericField("Number after removing worst",                                                  this.stSmplNum,  0);
  			gd.addNumericField("Maximal RMS of the remaining tiles in a sample",                               this.stSmplRms,  6);
  			gd.addCheckbox    ("Use window function for the samples",                                          this.stSmplWnd);
  			
  			gd.addNumericField("Grow initial selection before processing supertiles, odd - ortho. <0 - use all tiles",this.stGrowSel,  0);
  			gd.addNumericField("Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert",this.stMeasSel,  0);
  			gd.addNumericField("Consider merging initial planes if disparity difference below",                this.stSmallDiff,  6);
  			gd.addNumericField("Consider merging initial planes if jumps between ratio above",                 this.stHighMix,  6);

  			gd.addNumericField("Outlier tiles weaker than this may be replaced from neighbors",               this.outlayerStrength,  6);
  			gd.addNumericField("Replace weak outlier tiles that do not have neighbors within this disparity difference", this.outlayerDiff,  6);
  			gd.addNumericField("Replace weak outlier tiles that have higher disparity than weighted average", this.outlayerDiffPos,  6);
  			gd.addNumericField("Replace weak outlier tiles that have lower disparity than weighted average",  this.outlayerDiffNeg,  6);

  			gd.addCheckbox    ("Combine with all previous after refine pass",                                  this.combine_refine);
  			gd.addNumericField("Disregard weaker tiles when combining scans",                                  this.combine_min_strength,  6);
  			gd.addNumericField("Disregard weaker tiles when combining scans  for horizontal correlation",      this.combine_min_hor,  6);
  			gd.addNumericField("Disregard weaker tiles when combining scans  for vertical correlation",        this.combine_min_vert,  6);
//  			gd.addNumericField("Do not re-measure correlation if target disparity differs from some previous by this",this.unique_tolerance,  6);
  			
  			gd.addMessage     ("========= Growing disparity range to scan ========");
  			gd.addNumericField("Try these number of tiles around known ones",                                         this.grow_sweep,  0);
  			gd.addNumericField("Maximal disparity to try",                                                            this.grow_disp_max,  6);
  			gd.addNumericField("Trust measured disparity within +/- this value",                                      this.grow_disp_trust,  6);
  			gd.addNumericField("Increase disparity (from maximal tried) if nothing found in that tile",               this.grow_disp_step,  6);
  			gd.addNumericField("Grow more only if at least one channel has higher variance from others for the tile", this.grow_min_diff,  6);
  			gd.addCheckbox    ("Retry tiles around known foreground that have low max_tried_disparity",               this.grow_retry_far);
  			gd.addCheckbox    ("Scan full range between max_tried_disparity of the background and known foreground",  this.grow_pedantic);
  			gd.addCheckbox    ("Retry border tiles that were identified as infinity earlier",                         this.grow_retry_inf);

  			gd.addMessage     ("--- more growing parameters ---");
  			
  			gd.addCheckbox    ("New expansion mode",                                                                  this.gr_new_expand);
  			gd.addNumericField("Expansion steps limit",                                                               this.gr_max_expand,  0);
  			gd.addNumericField("Strength floor for multi-tile (now 5x5) samples (normally < combine_min_strength) ",  this.fds_str_floor,  6);
  			gd.addNumericField("Over background extra reliable strength",                                             this.gr_ovrbg_cmb,  6);
  			gd.addNumericField("Over background extra reliable strength horizontal",                                  this.gr_ovrbg_cmb_hor,  6);
  			gd.addNumericField("Over background extra reliable strength vertical",                                    this.gr_ovrbg_cmb_vert,  6);
  			gd.addNumericField("Over background filtered extra reliable strength",                                    this.gr_ovrbg_filtered,  6);
  			
  			gd.addMessage     ("--- \"plate\" filtering when growing parameters ---");
  			gd.addNumericField("Strength power exponent for tilted plates growing",                                   this.fds_str_pow,  6);
  			gd.addNumericField("Sample size (side of a square) for tilted plates growing",                            this.fds_smpl_side,  0);
  			gd.addNumericField("Number of tiles in a square tilted plate (should be >1)",                             this.fds_smpl_num,  0);
  			gd.addNumericField("Maximal RMS for the tiles to the tilted plate",                                       this.fds_smpl_rms,  6);
  			gd.addNumericField("Maximal relative RMS for the tiles to the tilted plate - multiply by disparity and add", this.fds_smpl_rel_rms,  6);
  			gd.addCheckbox    ("Use window function for the square sample plates",                                    this.fds_smpl_wnd);
  			gd.addNumericField("Maximal growing plate tilt in disparity pix per tile",                                this.fds_abs_tilt,  6);
  			gd.addNumericField("Maximal relative growing plate tilt in disparity pix per tile per disaprity pixel",   this.fds_rel_tilt,  6);

  			gd.addMessage     ("--- Macro correlation parameters ---");

  			gd.addNumericField("Macro disparity scan step (actual disparity step is 8x)",                             this.mc_disp8_step,  6);
  			gd.addNumericField("Trust measured macro(8x)  disparity within +/- this value ",                          this.mc_disp8_trust,  6);
  			gd.addNumericField("Minimal macro correlation strength to process",                                       this.mc_strength,  6);
  			gd.addNumericField("Do not re-measure macro correlation if target disparity differs less",                this.mc_unique_tol,  6);

  			gd.addNumericField("When consolidating macro results, exclude high residual disparity",                   this.mc_trust_fin,  6);
  			gd.addNumericField("Gaussian sigma to reduce weight of large residual disparity",                         this.mc_trust_sigma,  6);
  			gd.addNumericField("Weight from ortho neighbor supertiles",                                               this.mc_ortho_weight,  6);
  			gd.addNumericField("Weight from diagonal neighbor supertiles",                                            this.mc_diag_weight,  6);
  			gd.addNumericField("Do not remove measurements farther from the kept ones",                               this.mc_gap,  6);
  			gd.addMessage     ("--- more growing parameters ---");
  			
  			gd.addNumericField("Discard variant if it requests too few tiles",                                        this.gr_min_new,  0);
  			gd.addCheckbox    ("Expand only unambiguous tiles over previously undefined",                             this.gr_var_new_sngl);
  			gd.addCheckbox    ("Expand unambiguous and FOREGROUND tiles over previously UNDEFINED",                   this.gr_var_new_fg);
  			gd.addCheckbox    ("Expand unambiguous and FOREGROUND tiles over already DEFINED",                        this.gr_var_all_fg);
  			gd.addCheckbox    ("Expand unambiguous and BACKGROUND tiles over previously UNDEFINED",                   this.gr_var_new_bg);
  			gd.addCheckbox    ("Expand unambiguous and BACKGROUND tiles over already DEFINED",                        this.gr_var_all_bg);
  			gd.addCheckbox    ("Try next disparity range that was not tried before",                                  this.gr_var_next);
  			gd.addNumericField("How far to extend over previously undefined disparity tiles",                         this.gr_num_steps,  0);
  			gd.addNumericField("How far to extend over previously determined disparity tiles",                        this.gr_steps_over,  0);
  			gd.addNumericField("Extend sample plate square side",                                                     this.gr_smpl_size,  0);
  			gd.addNumericField("Extend at least this number of the seed tiles",                                       this.gr_min_pnts,  0);
  			gd.addCheckbox    ("Use window function for square sample ",                                              this.gr_use_wnd);
  			gd.addNumericField("Tilt cost for damping insufficient plane data",                                       this.gr_tilt_damp,  6);
  			gd.addNumericField("When growing, range of disparities to be extended without far/near subdivision",      this.gr_split_rng,  6);
  			gd.addNumericField("Consider far/near tiles within that range from the farthest/closest",                 this.gr_same_rng,  6);
  			gd.addNumericField("Maximal difference from the old value when smoothing",                                this.gr_diff_cont,  6);
  			gd.addNumericField("Maximal filter disparity absolute tilt (pix per tile)",                               this.gr_abs_tilt,  6);
  			gd.addNumericField("Maximal filter disparity tilt (pix / disparity) per tile",                            this.gr_rel_tilt,  6);
  			gd.addNumericField("Maximal number of smoothing steps (reduce if long?)",                                 this.gr_smooth,  0);
  			gd.addNumericField("Maximal change to finish smoothing iterations",                                       this.gr_fin_diff,  6);
  			gd.addNumericField("Do not re-measure correlation if target disparity differs from some previous less",   this.gr_unique_tol,  6);
  			gd.addNumericField("Larger tolerance for expanding (not refining)",                                       this.gr_unique_pretol,  6);
  			

  			gd.addMessage     ("--- Planes detection ---");
  			gd.addCheckbox    ("Always start with disparity-most axis (false - lowest eigenvalue)",            this.plPreferDisparity);
  			gd.addNumericField("Normalize disparities to the average if above",                                this.plDispNorm,  6);

  			gd.addNumericField("Blur disparity histograms for constant disparity clusters by this sigma (in bins)",   this.plBlurBinVert,  6);
  			gd.addNumericField("Blur disparity histograms for horizontal clusters by this sigma (in bins)",           this.plBlurBinHor,  6);
  			gd.addNumericField("Maximal normalized disparity difference when initially assigning to vertical plane",  this.plMaxDiffVert,  6);
  			gd.addNumericField("Maximal normalized disparity difference when initially assigning to horizontal plane",this.plMaxDiffHor,  6);
  			gd.addNumericField("Number of initial passes to assign tiles to vert (const disparity) and hor planes",   this.plInitPasses,  0);
  			
  			gd.addNumericField("Minimal number of points for plane detection",                                 this.plMinPoints,  0);
  			gd.addNumericField("Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below", this.plTargetEigen,  6);
  			gd.addNumericField("Maximal fraction of outliers to remove",                                       this.plFractOutliers,  6);
  			gd.addNumericField("Maximal number of outliers to remove",                                         this.plMaxOutliers,  0);
  			gd.addNumericField("Minimal total strength of a plane",                                            this.plMinStrength,  6);
  			gd.addNumericField("Maximal eigenvalue of a plane",                                                this.plMaxEigen,  6);
  			gd.addNumericField("Add to eigenvalues of each participating plane and result to validate connections",this.plEigenFloor,  6);
  			gd.addNumericField("Consider plane to be a \"stick\" if second eigenvalue is below",               this.plEigenStick,  6);
  			gd.addNumericField("Not a plate if sin^2 between normals from disparity and world exceeds this",   this.plBadPlate,  6);
  			gd.addCheckbox    ("Combine 'other' plane with the current (unused)",                              this.plDbgMerge);
  			gd.addNumericField("Worst case worsening after merge",                                             this.plWorstWorsening,  6);
  			gd.addNumericField("Worst case worsening for thin planes",                                         this.plWorstWorsening2,  6);
  			gd.addNumericField("Worst case worsening after merge with equal weights",                          this.plWorstEq,  6);
  			gd.addNumericField("Worst case worsening for thin planes with equal weights",                      this.plWorstEq2,  6);
  			gd.addNumericField("If result of the merged planes is below, OK to use thin planes (higher) threshold ",this.plOKMergeEigen,  6);
  			gd.addNumericField("Maximal sine squared of the world angle between planes to merge. Set to >= 1.0 to disable", this.plMaxWorldSin2,  6);
  			gd.addNumericField("Do not compare sin() between planes, if at least one has too small axis ratio",this.pl2dForSin,  6);
  			gd.addNumericField("Relax merge requirements for weaker planes",                                   this.plWeakWorsening,  6);
  			gd.addNumericField("Maximal overlap between the same supertile planes to merge",                   this.plMaxOverlap,  6);

  			gd.addMessage     ("--- Merge same supetile planes if at least one is weak and they do not differ much ---");
  			gd.addNumericField("Maximal weight of the weak plane to merge (first variant)",                    this.plWeakWeight,  6);
  			gd.addNumericField("Maximal eigenvalue of the result of non-weighted merge (first variant)",       this.plWeakEigen,  6);
  			gd.addNumericField("Maximal weight of the weak plane to merge (second variant)",                   this.plWeakWeight2,  6);
  			gd.addNumericField("Maximal eigenvalue of the result of non-weighted merge (second variant)",      this.plWeakEigen2,  6);
  			gd.addNumericField("Do not merge if any sqrt of merged eigenvalue exceeds scaled sum of components", this.plSumThick,  6);
  			gd.addNumericField("When calculating non-exclusive planes, do not use neighbors with high cost",   this.plNeNeibCost,  6);
  			gd.addNumericField("When calculating non-exclusive planes, use cenrter plane relative weight",     this.plNeOwn,  6);
  			gd.addNumericField("When calculating exclusive planes links, do not use neighbors with high cost", this.plExNeibCost,  6);
  			gd.addNumericField("Scale down maximal costs for smoothed planes links (tighter requirements)",    this.plExNeibSmooth,  6);
  			gd.addNumericField("Cost threshold for merging same tile planes if the plane has connected neighbors", this.plMergeCostStar,  6);
  			gd.addNumericField("Cost threshold for merging same tile planes if not connected",                 this.plMergeCost,  6);

  			gd.addMessage     ("--- Merging planes with topological conflicts ---");
  			gd.addCheckbox    ("Try to merge conflicting planes",                                                     this.plConflMerge);
  			gd.addNumericField("Scale parameters to relax planes fit for merging conflicting planes",                 this.plConflRelax,  6);
  			gd.addCheckbox    ("Only merge conflicting planes if this is the only conflicting pair in the supertile", this.plConflSngl);
  			gd.addCheckbox    ("Only merge conflicting planes only if there are just two planes in the supertile",    this.plConflSnglPair);

  			gd.addNumericField("Consider merging plane if it is foreground and maximal strength below this",          this.plWeakFgStrength,  6);
  			gd.addNumericField("Remove these strongest from foreground when determining the maximal strength",        this.plWeakFgOutliers,  0);
  			gd.addNumericField("Relax cost requirements when merging with weak foreground",                           this.plWeakFgRelax,  6);
  			
  			gd.addNumericField("Maximal real-world thickness of merged overlapping planes (meters)",                  this.plThickWorld,  6);
  			gd.addNumericField("Maximal real-world merged thickness for conflicting planes",                          this.plThickWorldConfl,  6);
  			gd.addNumericField("Relax cost requirements when adding exclusive links to complete squares and triangles",this.plRelaxComplete,  6);
  			gd.addNumericField("Relax cost requirements more during the second pass",                                 this.plRelaxComplete2,  6);
  			
  			gd.addMessage     ("---  ---");
  			gd.addNumericField("Maximal ratio of Z to allow plane merging",                                    this.plMaxZRatio,  6);
  			gd.addNumericField("Maximal disparity of one of the planes to apply  maximal ratio",               this.plMaxDisp,  6);
  			gd.addNumericField("When merging with neighbors cut the tail that is worse than scaled best",      this.plCutTail,  6);
  			gd.addNumericField("Set cutoff value level not less than this",                                    this.plMinTail,  6);
  			
  			gd.addMessage     ("--- Parameters to regenerate planes using preliminary tile connections  ---");
  			gd.addCheckbox    ("Enable planes tiles selection regeneration hinted by supertile neighbors",     this.plDiscrEn);
  			gd.addNumericField("Maximal disparity difference from the plane to consider tile",                 this.plDiscrTolerance,  6);
  			gd.addNumericField("Parallel move known planes around original know value for the best overall fit",this.plDiscrDispRange,  6);
  			gd.addNumericField("Number of steps (each direction) for each plane to search for the best fit (0 - single, 1 - 1 each side)",this.plDiscrSteps,  0);
  			gd.addNumericField("Total number of variants to try (protect from too many planes)",               this.plDiscrVariants,  0);
  			gd.addNumericField("What plane to use as a hint: 0 - weighted, 1 - equalized, 2 - best, 3 - combined", this.plDiscrMode,  0);
  			gd.addNumericField("Squared add to variance to calculate reverse flatness (used mostly for single-cell clusters)",this.plDiscrVarFloor,  6);
  			gd.addNumericField("Gaussian sigma to compare how measured data is attracted to planes",           this.plDiscrSigma,  6);
  			gd.addNumericField("Sigma to blur histograms while re-discriminating",                             this.plDiscrBlur,  6);
  			gd.addNumericField("Tile exclusivity: 1.0 - tile belongs to one plane only, 0.0 - regardless of others",this.plDiscrExclusivity,  6);
  			gd.addNumericField("For second pass if exclusivity > 1.0 - will assign only around strong neighbors",this.plDiscrExclus2,  6);
  			gd.addCheckbox    ("When growing selection do not allow any offenders around (false - more these than others)", this.plDiscrStrict);
  			gd.addNumericField("Attraction to different planes correlation that is too high for re-discrimination",this.plDiscrCorrMax,  6);
  			gd.addNumericField("Attraction to different planes correlation that is high enough to merge",      this.plDiscrCorrMerge,  6);
  			gd.addNumericField("If offender has this number of tiles (including center) the cell can not be used", this.plDiscrSteal,  0);
  			gd.addNumericField("Only use tiles within this range from original selection",                     this.plDiscrGrown,  0);
  			gd.addNumericField("Remove outliers from the final selection that have distance more than scaled median",this.plDiscrXMedian,  6);

  			gd.addMessage     ("--- Planes merge costs ---");
  			gd.addNumericField(" Disparity (pix) - closer cost will use more of the real world, farther - disparity",this.plCostDist,  6);
  			gd.addNumericField("Cost of merge quality sqrt(weighted*equal) in disparity space",                this.plCostKrq,     6);
  			gd.addNumericField("Cost of merge quality averaje of weighted and equal weight in disparity space",this.plCostKrqEq,   6);
  			gd.addNumericField("Cost of merge quality sqrt(weighted*equal) in world space",                    this.plCostWrq,     6);
  			gd.addNumericField("Cost of merge quality average of weighted and equal weight in world space",    this.plCostWrqEq,   6);
  			gd.addNumericField("Cost of sin squared between normals",                                          this.plCostSin2,    6);
  			gd.addNumericField("Cost of squared relative plane-to-other-center distances",                     this.plCostRdist2,  6);

  			gd.addCheckbox    ("Resolve dual triangles conflict (odoodo)",                                     this.plConflDualTri);
  			gd.addCheckbox    ("Resolve multiple odo triangles conflicts",                                     this.plConflMulti);
  			gd.addCheckbox    ("Resolve diagonal (ood) conflicts",                                             this.plConflDiag);
  			gd.addCheckbox    ("Resolve all conflicts around a supertile",                                     this.plConflStar);
  			gd.addNumericField("How far to look around when calculationg connection cost",                     this.plStarSteps,  0);
  			gd.addNumericField("When calculating cost for the connections scale 4 ortho neighbors",            this.plStarOrtho,  6);
  			gd.addNumericField("When calculating cost for the connections scale 4 diagonal neighbors",         this.plStarDiag,  6);
  			gd.addNumericField("Divide cost by number of connections to this power",                           this.plStarPwr,  6);
  			gd.addNumericField("Use this power of tile weight when calculating connection cost",               this.plStarWeightPwr,  6);
  			gd.addNumericField("Balance weighted density against density. 0.0 - density, 1.0 - weighted density", this.plWeightToDens,  6);
  			gd.addNumericField("Raise value of each tile before averaging",                                    this.plStarValPwr,  6);
  			gd.addNumericField("When resolving double triangles allow minor degradation (0.0 - strict)",       this.plDblTriLoss,  6);
  			gd.addCheckbox    ("Allow more conflicts if overall cost is reduced",                              this.plNewConfl);
  			gd.addNumericField("aximal number of simultaneous connection changes around one tile (0 - any)",   this.plMaxChanges,  0);

  			gd.addCheckbox    ("Keep only mutual links, remove weakest if conflict",                           this.plMutualOnly);

  			gd.addCheckbox    ("Add diagonals to full squares",                                                this.plFillSquares);
  			gd.addCheckbox    ("Add ortho to 45-degree corners",                                               this.plCutCorners);
  			gd.addCheckbox    ("Add hypotenuse connection if both legs exist",                                 this.plHypotenuse);

  			gd.addNumericField("Relative (to average neighbor) weight of the measured plane when combing with neighbors",     this.plPull,  6);
  			gd.addNumericField("0.0: 8 neighbors pull 8 times as 1, 1.0 - same as 1",                          this.plNormPow,  6);
  			gd.addNumericField("Maximal number of smoothing iterations for each step",                         this.plIterations,  0);
  			gd.addCheckbox    ("Do not update supertile if any of connected is not good (false: just skip that neighbor)", this.plStopBad);
  			gd.addNumericField("Maximal step difference (1/power of 10)",                                      this.plPrecision,  0);

  			gd.addNumericField("Relative weight of center plane when splitting into pairs",                    this.plSplitPull,  6);
  			gd.addNumericField("Minimal number of neighbors to split plane in pairs",                          this.plSplitMinNeib,  0);
  			gd.addNumericField("Minimal weight of split plains to show",                                       this.plSplitMinWeight,  6);
  			gd.addNumericField("Minimal split quality to show",                                                this.plSplitMinQuality,  6);
  			gd.addCheckbox    ("Apply plane split to pairs",                                                   this.plSplitApply);
  			gd.addCheckbox    ("Allow tiles to belong to both planes of the pair",                             this.plNonExclusive);
  			gd.addCheckbox    ("Allow other tiles from the same supertile",                                    this.plUseOtherPlanes);
  			gd.addCheckbox    ("Allow parallel shift of the specified planes before adding",                   this.plAllowParallel);
  			gd.addNumericField("Maximal normalized tile disparity difference from the plane to consider",      this.plMaxDiff,  6);
  			gd.addNumericField("Maximal difference of the added tile ratio to the average  disparity difference",this.plOtherDiff,  6);
  			gd.addCheckbox    ("Separate tiles for split planes by X, Y",                                      this.plSplitXY);
  			gd.addNumericField("Disparity tolerance when separating by X, Y",                                  this.plSplitXYTolerance,  6);
  			
  			gd.addCheckbox    ("Fuse planes together (off for debug only)",                                    this.plFuse);
  			gd.addCheckbox    ("Keep unconnected supertiles",                                                  this.plKeepOrphans);
  			gd.addNumericField("Minimal strength unconnected supertiles to keep",                              this.plMinOrphan,  6);

  			gd.addMessage     ("--- Snap to planes ---");
  			gd.addNumericField("Maximal (scaled by plDispNorm) disparity difference to snap to plane at any strength",     this.plSnapDispAny,  6);
  			gd.addNumericField("Maximal strength to fit any distance (if does not fit otherwise - treat as zero strength", this.plSnapStrengthAny,  6);
  			gd.addNumericField("Maximal negative disparity difference from the best match",                    this.plSnapNegAny,  6);
  			gd.addNumericField("Maximal (scaled by plDispNorm) disparity difference to snap to plane at low strength",     this.plSnapDispMax,  6);
  			gd.addNumericField("Maximal disparity diff. by weight product to snap to plane",                   this.plSnapDispWeight,  6);
  			gd.addNumericField("Zero strength snap mode: 0: no special treatment, 1 - strongest, 2 - farthest",this.plSnapZeroMode,  0);

  			gd.addCheckbox    ("Use planes selection masks (generated when splitting to intersecting pairs",   this.msUseSel);
  			gd.addCheckbox    ("Divide plane strengths by ellipsoid area",                                     this.msDivideByArea);
  			gd.addNumericField("Scale projection of the plane ellipsoid",                                      this.msScaleProj,  6);
  			gd.addNumericField("Spread this fraction of the ellipsoid weight among extended (double) supertile",this.msFractUni,  6);

  			gd.addMessage     ("--- Tiles assignment ---");
  			 
  			gd.addCheckbox    ("Do not assign tiles to the surface edges (not having all 8 neighbors)",           this.tsNoEdge);
  			gd.addCheckbox    ("Only assign outside of 8x8 center if no suitable alternative",                    this.tsUseCenter);
  			gd.addNumericField("Maximal disparity difference when assigning tiles",                               this.tsMaxDiff,  6);
  			gd.addNumericField("Minimal disparity difference to be considered as a competitor surface",           this.tsMinDiffOther,  6);
  			gd.addNumericField("Minimal tile correlation strength to be assigned",                                this.tsMinStrength,  6);
  			gd.addNumericField("Maximal tile correlation strength to be assigned",                                this.tsMaxStrength,  6);
  			gd.addNumericField("Minimal surface strength at the tile location",                                   this.tsMinSurface,  6);
  			gd.addNumericField("Allowed tile disparity correction: 1 increase, 2 - decrease, 3 - both directions",this.tsMoveDirs,  0);
  			gd.addNumericField("Raise surface strengths ratio to this power when comparing candidates",           this.tsSurfStrPow,  6);
  			gd.addNumericField("Add to strengths when calculating pull of assigned tiles",                        this.tsAddStrength,  6);
  			gd.addNumericField("Radius of influence (in tiles) of the previously assigned tiles",                 this.tsSigma,  6);
  			gd.addNumericField("Maximal relative to radius distance to calculate influence",                      this.tsNSigma,  6);
  			gd.addNumericField(" Additional pull of each surface ",                                               this.tsMinPull,  6);
  			gd.addNumericField("Minimal ratio of the best surface candidate to the next one to make selection",   this.tsMinAdvantage,  6);
  			
  			gd.addNumericField("Minimal size of a cluster to keep",                                               this.tsClustSize,  0);
  			gd.addNumericField("Minimal total weight of a cluster to keep",                                       this.tsClustWeight,  6);
  			gd.addNumericField("Minimal number of neighbors of unassigned tile to join (the farthest)",           this.tsMinNeib,  0);
  			gd.addNumericField("Maximal strength of the surrounded unassigned tile to join",                      this.tsMaxSurStrength,  6);
  			gd.addCheckbox    ("Include disabled tiles/borders when counting assigned neighbors",                 this.tsCountDis);
  			
  			gd.addCheckbox    ("Assign tiles that were used to generate planes",                                  this.tsEnPlaneSeed);
  			gd.addCheckbox    ("Allow assignment only surface",                                                   this.tsEnOnly);
  			gd.addCheckbox    ("Grow the only surface assignments",                                               this.tsEnGrow);
  			gd.addNumericField("Maximal strength when growing the only surfaces",                                 this.tsGrowStrength,  6);
  			gd.addCheckbox    ("Grow over strong if disparity matches",                                           this.tsGrowStrong);
  			gd.addNumericField("Minimal strength to continue grow with disparity match",                          this.tsContStrength,  6);
  			gd.addNumericField("Maximal normalized disparity error to grow over strong tiles",                    this.tsContDiff,  6);
  			
  			gd.addCheckbox    ("Allow assignment to the nearest surface with no competitors",                     this.tsEnSingle);
  			gd.addCheckbox    ("Allow assignment when several surfaces fit",                                      this.tsEnMulti);
  			gd.addCheckbox    ("Remove weak clusters before growing",                                             this.tsRemoveWeak1);
  			gd.addCheckbox    ("Assign tiles that have neighbors to the lowest disparity",                        this.tsGrowSurround);
  			gd.addCheckbox    ("Remove weak clusters after growing",                                              this.tsRemoveWeak2);
  			
  			gd.addCheckbox    ("Repeat multi-choice assignment while succeeding",                                 this.tsLoopMulti);
  			gd.addCheckbox    ("Show results of tiles to surfaces assignment",                                    this.tsShow);
  			gd.addNumericField("Number of clusters to keep",                                                      this.tsNumClust,  0);
  			
  			gd.addNumericField("Which assignments to match +1 - combo, +2 grown single, +4 plane seeds",          this.tsConsensMode,  0);
  			gd.addNumericField("Minimal number of assignments to agree",                                          this.tsConsensAgree,  0);

  			gd.addMessage     ("--- Tile assignment parameters ---");
  			gd.addNumericField("Minimal foreground/ background separation to look for weak FG edge",              this.taMinFgBg,    6);
  			gd.addNumericField("Minimal foreground edge strength (stronger edges will have proportionally smaller costs)", this.taMinFgEdge,  6);
  			gd.addNumericField("Minimal surface separation that requires color change",                           this.taMinColSep,  6);
  			gd.addNumericField("Minimal color variation (larger proportionally reduces cost)",                    this.taMinColDiff, 6);
  			gd.addNumericField("Disparity difference limit (to handle outliers)",                                 this.taOutlier, 6);
  			gd.addNumericField("Strength power when calculating disparity error",                                 this.taDiffPwr, 6);
  			gd.addNumericField("Strength power when calculating disparity error over best",                       this.taBestPwr, 6);
  			gd.addNumericField("Strength power when calculating disparity error for group of 9",                  this.taDiff9Pwr, 6);
  			gd.addNumericField("Gaussian sigma to blur color difference between tiles along each direction",      this.taColSigma, 6);
  			gd.addNumericField("Relative amount of the blurred color difference in the mixture",                  this.taColFraction, 6);

  			gd.addNumericField("Cost of a tile that is not assigned",                                             this.taCostEmpty,  6);
  			gd.addNumericField("Cost of a tile not having any neighbor in particular direction",                  this.taCostNoLink,  6);
  			gd.addNumericField("Cost of a tile switching to a neighbor that does not have a link",                this.taCostSwitch,  6);
  			gd.addNumericField("Cost of a tile switching to a disconnected neighbor divided by a color",          this.taCostColor,  6);
  			gd.addNumericField("Cost of a weighted normalized tile disparity error",                              this.taCostDiff,  6);
  			gd.addNumericField("Cost of a weighted normalized tile disparity error above best surface",           this.taCostDiffBest,  6);
  			gd.addNumericField("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",this.taCostDiff9,  6);
  			gd.addNumericField("Cost of a weak foreground edge",                                                  this.taCostWeakFgnd,  6);
  			gd.addNumericField("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",            this.taCostFlaps,  6);
  			gd.addNumericField("Cost of a measurement layer not having same layer in the same location or near",  this.taCostMismatch,  6);

  			gd.addCheckbox    ("Cost of a tile that is not assigned",                                             this.taEnEmpty);
  			gd.addCheckbox    ("Cost of a tile not having any neighbor in particular direction",                  this.taEnNoLink);
  			gd.addCheckbox    ("Cost of a tile switching to a neighbor that does not have a link",                this.taEnSwitch);
  			gd.addCheckbox    ("Cost of a tile switching to a disconnected neighbor divided by a color",          this.taEnColor);
  			gd.addCheckbox    ("Cost of a weighted normalized tile disparity error",                              this.taEnDiff);
  			gd.addCheckbox    ("Cost of a weighted normalized tile disparity error above best surface",           this.taEnDiffBest);
  			gd.addCheckbox    ("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",this.taEnDiff9);
  			gd.addCheckbox    ("Cost of a weak foreground edge",                                                  this.taEnWeakFgnd);
  			gd.addCheckbox    ("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",            this.taEnFlaps);
  			gd.addCheckbox    ("Cost of a measurement layer not having same layer in the same location or near",  this.taEnMismatch);
  			
  			gd.addCheckbox    ("Test new mode after migration",                                                this.dbg_migrate);

  			gd.addMessage     ("--- Other debug images ---");
  			gd.addCheckbox    ("Show 'ortho_combine'",                                                         this.show_ortho_combine);
  			gd.addCheckbox    ("Show 'refine_disparity_supertiles'",                                           this.show_refine_supertiles);
  			gd.addCheckbox    ("Show 'bgnd_nonbgnd'",                                                          this.show_bgnd_nonbgnd);
  			gd.addCheckbox    ("Show 'FilterScan'",                                                            this.show_filter_scan);
  			gd.addCheckbox    ("Show 'combo_scan' (combined multiple scans)",                                  this.show_combined);
  			gd.addCheckbox    ("Show 'unique_scan' (removed already measured tiles with the same disparity)",  this.show_unique);
  			gd.addCheckbox    ("Show supertile disparity histograms ",                                         this.show_histograms);
  			gd.addCheckbox    ("Show debug images during initial refinement",                                  this.show_init_refine);
  			gd.addCheckbox    ("Show debug images during disparity expansion",                                 this.show_expand);
  			gd.addCheckbox    ("Show prepareExpandVariant when elevating variant",                             this.show_variant);
  			gd.addCheckbox    ("Show debug images related to retrying far tiles near foreground",              this.show_retry_far);
  			gd.addCheckbox    ("Show debug images related to macro correlation",                               this.show_macro);
  			gd.addCheckbox    ("Show 'shells'",                                                                this.show_shells);
  			gd.addCheckbox    ("show 'neighbors'",                                                             this.show_neighbors);
  			gd.addCheckbox    ("Show 'flaps-dirs'",                                                            this.show_flaps_dirs);
  			gd.addCheckbox    ("Show 'first_N_clusters'",                                                      this.show_first_clusters);
  			gd.addCheckbox    ("Show planes",                                                                  this.show_planes);
  			gd.addMessage     ("Unity up vector in camera coordinate system (x - right, y - up, z - to camera): {"+
  			this.vertical_xyz[0]+","+this.vertical_xyz[1]+","+this.vertical_xyz[2]+"}");
  			WindowTools.addScrollBars(gd);
  			gd.showDialog();
  			
  			if (gd.wasCanceled()) return false;
  			this.transform_size=  (int) gd.getNextNumber();
  			this.clt_window=      (int) gd.getNextNumber();
  			this.shift_x =              gd.getNextNumber();
  			this.shift_y =              gd.getNextNumber();
  			this.iclt_mask=       (int) gd.getNextNumber();
  			this.tileX=           (int) gd.getNextNumber();
  			this.tileY=           (int) gd.getNextNumber();
  			this.dbg_mode=        (int) gd.getNextNumber();
  			this.ishift_x=        (int) gd.getNextNumber();
  			this.ishift_y=        (int) gd.getNextNumber();
  			this.fat_zero =             gd.getNextNumber();
  			this.corr_sigma =           gd.getNextNumber();
  			this.norm_kern=             gd.getNextBoolean();
  			this.gain_equalize=         gd.getNextBoolean();
  			this.colors_equalize=       gd.getNextBoolean();
  			this.novignetting_r=        gd.getNextNumber();
  			this.novignetting_g=        gd.getNextNumber();
  			this.novignetting_b=        gd.getNextNumber();
  			this.scale_r=               gd.getNextNumber();
  			this.scale_g=               gd.getNextNumber();
  			this.scale_b=               gd.getNextNumber();
  			this.vignetting_max=        gd.getNextNumber();
  			this.vignetting_range=      gd.getNextNumber();
  			this.kernel_step=     (int) gd.getNextNumber();
  			this.disparity=             gd.getNextNumber();
  			this.z_correction=          gd.getNextNumber();
  			this.correlate=             gd.getNextBoolean();
  			this.corr_mask=       (int) gd.getNextNumber();
  			this.corr_sym=              gd.getNextBoolean();
  			this.corr_keep=             gd.getNextBoolean();
  			this.corr_show=             gd.getNextBoolean();
  			this.corr_mismatch=         gd.getNextBoolean();
  			this.corr_offset=           gd.getNextNumber();
  			this.corr_red=              gd.getNextNumber();
  			this.corr_blue=             gd.getNextNumber();
  			this.corr_normalize=        gd.getNextBoolean();
  			this.min_corr=              gd.getNextNumber();
  			this.min_corr_normalized=   gd.getNextNumber();
  			this.max_corr_sigma=        gd.getNextNumber();
  			this.max_corr_radius=       gd.getNextNumber();
  			
  			this.enhortho_width= (int) gd.getNextNumber();
  			this.enhortho_scale=        gd.getNextNumber();

  			this.max_corr_double=       gd.getNextBoolean();
  			this.corr_mode=       (int) gd.getNextNumber();
  			this.corr_border_contrast=  gd.getNextNumber();
  			
//  			this.tile_task_op=    (int) gd.getNextNumber();
  			this.tile_task_op = ImageDtt.setOrthoLines      (this.tile_task_op, gd.getNextBoolean());
  			this.tile_task_op = ImageDtt.setForcedDisparity (this.tile_task_op, gd.getNextBoolean());
  			this.tile_task_op = ImageDtt.setImgMask         (this.tile_task_op, (int) gd.getNextNumber());
  			this.tile_task_op = ImageDtt.setPairMask        (this.tile_task_op, (int) gd.getNextNumber());
  			
  			this.tile_task_wl=    (int) gd.getNextNumber();
  			this.tile_task_wt=    (int) gd.getNextNumber();
  			this.tile_task_ww=    (int) gd.getNextNumber();
  			this.tile_task_wh=    (int) gd.getNextNumber();
  			this.min_shot=              gd.getNextNumber();
  			this.scale_shot=            gd.getNextNumber();
  			this.diff_sigma=            gd.getNextNumber();
  			this.diff_threshold=        gd.getNextNumber();
  			this.diff_gauss=            gd.getNextBoolean();
  			this.min_agree=             gd.getNextNumber();
  			this.dust_remove=           gd.getNextBoolean();
  			this.black_back=            gd.getNextBoolean();
  			this.keep_weights=          gd.getNextBoolean();
  			this.sharp_alpha=           gd.getNextBoolean();
  			this.alpha0=                gd.getNextNumber();
  			this.alpha1=                gd.getNextNumber();
  			this.gen_chn_stacks=        gd.getNextBoolean();
  			this.gen_chn_img=           gd.getNextBoolean();
  			this.gen_4_img=             gd.getNextBoolean();
  			this.show_nonoverlap=       gd.getNextBoolean();
  			this.show_overlap=          gd.getNextBoolean();
  			this.show_rgba_color=       gd.getNextBoolean();
  			this.show_map=              gd.getNextBoolean();
  			this.disp_scan_start=       gd.getNextNumber();
  			this.disp_scan_step=        gd.getNextNumber();
  			this.disp_scan_count= (int) gd.getNextNumber();
  			
  			this.fine_dbg=              gd.getNextBoolean();
  			this.fine_corr_x_0=         gd.getNextNumber();
  			this.fine_corr_y_0=         gd.getNextNumber();
  			this.fine_corr_x_1=         gd.getNextNumber();
  			this.fine_corr_y_1=         gd.getNextNumber();
  			this.fine_corr_x_2=         gd.getNextNumber();
  			this.fine_corr_y_2=         gd.getNextNumber();
  			this.fine_corr_x_3=         gd.getNextNumber();
  			this.fine_corr_y_3=         gd.getNextNumber();
  			
  			this.fcorr_radius=         gd.getNextNumber();
  			this.fcorr_min_strength=    gd.getNextNumber();
  			this.fcorr_disp_diff=       gd.getNextNumber();
  			this.fcorr_quadratic=       gd.getNextBoolean();
  			this.fcorr_ignore=          gd.getNextBoolean();

  			this.fcorr_inf_strength=    gd.getNextNumber();
  			this.fcorr_inf_diff=        gd.getNextNumber();
  			this.fcorr_inf_quad=        gd.getNextBoolean();
  			this.fcorr_inf_vert=        gd.getNextBoolean();
  			
  			this.inf_disp_apply=        gd.getNextBoolean();
  			this.inf_repeat=      (int) gd.getNextNumber();
//  			this.inf_mism_apply=        gd.getNextBoolean();
  			this.inf_iters=       (int) gd.getNextNumber();
  			this.inf_final_diff=        gd.getNextNumber();
  			this.inf_far_pull=          gd.getNextNumber();
  			
  			this.inf_str_pow=           gd.getNextNumber();
  			this.inf_smpl_side=   (int) gd.getNextNumber();
  			this.inf_smpl_num=    (int) gd.getNextNumber();
  			this.inf_smpl_rms=          gd.getNextNumber();

  			this.ih_smpl_step=    (int) gd.getNextNumber();
  			this.ih_disp_min=           gd.getNextNumber();
  			this.ih_disp_step=          gd.getNextNumber();
  			this.ih_num_bins=     (int) gd.getNextNumber();
  			this.ih_sigma=              gd.getNextNumber();
  			this.ih_max_diff=           gd.getNextNumber();
  			this.ih_min_samples=  (int) gd.getNextNumber();
  			this.ih_norm_center=        gd.getNextBoolean();

			this.ly_smpl_side=    (int) gd.getNextNumber();
			this.ly_smpl_num=     (int) gd.getNextNumber();
			this.ly_meas_disp=          gd.getNextNumber();
			this.ly_smpl_rms=           gd.getNextNumber();
			this.ly_disp_var=           gd.getNextNumber();
			this.ly_inf_frac=           gd.getNextNumber();
  			this.ly_on_scan=            gd.getNextBoolean();
  			this.ly_inf_en=             gd.getNextBoolean();
  			this.ly_inf_force=          gd.getNextBoolean();
  			this.ly_poly=               gd.getNextBoolean();
  			
//  			this.fcorr_sample_size= (int)gd.getNextNumber();
//  			this.fcorr_mintiles= (int)  gd.getNextNumber();
//  			this.fcorr_reloutliers=     gd.getNextNumber();
//  			this.fcorr_sigma=           gd.getNextNumber();

  			this.corr_magic_scale=      gd.getNextNumber();

  			this.show_textures=         gd.getNextBoolean();
  			this.debug_filters=         gd.getNextBoolean();
  			this.min_smth=              gd.getNextNumber();
  			this.sure_smth=             gd.getNextNumber();
  			this.bgnd_range=            gd.getNextNumber();
  			this.other_range=           gd.getNextNumber();

  			this.ex_strength=           gd.getNextNumber();
  			this.ex_nstrength=          gd.getNextNumber();

  			this.ex_over_bgnd=          gd.getNextBoolean();
  			this.ex_min_over=           gd.getNextNumber();
  			
  			this.pt_super_trust=        gd.getNextNumber();
  			this.pt_keep_raw_fg=        gd.getNextBoolean();
  			this.pt_scale_pre=          gd.getNextNumber();
  			this.pt_scale_post=         gd.getNextNumber();

  			this.bgnd_sure=             gd.getNextNumber();
  			this.bgnd_maybe=            gd.getNextNumber();
  			this.min_clstr_seed=  (int) gd.getNextNumber();
  			this.min_clstr_lone=  (int) gd.getNextNumber();
  			this.min_clstr_weight=      gd.getNextNumber();
  			this.min_clstr_max=         gd.getNextNumber();
  			
  			this.fill_gaps=       (int) gd.getNextNumber();
  			this.fill_final=      (int) gd.getNextNumber();
  			this.min_clstr_block= (int) gd.getNextNumber();
  			this.bgnd_grow=       (int) gd.getNextNumber();

  			this.ortho_old=             gd.getNextBoolean();
  			this.ortho_min_hor=         gd.getNextNumber();
  			this.ortho_min_vert=        gd.getNextNumber();
  			this.ortho_asym=            gd.getNextNumber();
  			this.ortho_over4=           gd.getNextNumber();
  			
  			this.ortho_sustain=         gd.getNextNumber();
  			this.ortho_run=       (int) gd.getNextNumber();
  			this.ortho_minmax=          gd.getNextNumber();
  			this.ortho_bridge=    (int) gd.getNextNumber();
  			this.ortho_rms=             gd.getNextNumber();
  			this.ortho_half_length=(int)gd.getNextNumber();
  			this.ortho_mix=             gd.getNextNumber();

  			this.or_hor=                gd.getNextBoolean();
  			this.or_vert=               gd.getNextBoolean();
  			this.or_sigma=              gd.getNextNumber();
  			this.or_sharp=              gd.getNextNumber();
  			this.or_scale=              gd.getNextNumber();
  			this.or_offset=             gd.getNextNumber();
  			this.or_asym=               gd.getNextNumber();
  			this.or_threshold=          gd.getNextNumber();
  			this.or_absHor=             gd.getNextNumber();
  			this.or_absVert=            gd.getNextNumber();

  			this.poles_fix=             gd.getNextBoolean();
  			this.poles_len=        (int)gd.getNextNumber();
  			this.poles_ratio=           gd.getNextNumber();
  			this.poles_min_strength=    gd.getNextNumber();
  			this.poles_force_disp=      gd.getNextBoolean();

  			this.max_clusters=    (int) gd.getNextNumber();
  			this.remove_scans=          gd.getNextBoolean();
  			this.output_x3d=            gd.getNextBoolean();
  			this.output_obj=            gd.getNextBoolean();
  			this.correct_distortions=   gd.getNextBoolean();
  			this.show_triangles=        gd.getNextBoolean();
  			this.avg_cluster_disp=      gd.getNextBoolean();
  			this.maxDispTriangle=       gd.getNextNumber();
  			this.infinityDistance=      gd.getNextNumber();
  			this.shUseFlaps=            gd.getNextBoolean();
  			this.shAggrFade=            gd.getNextBoolean();
  			this.shMinArea=       (int) gd.getNextNumber();
  			this.shMinStrength=         gd.getNextNumber();
  			this.tiRigidVertical=       gd.getNextNumber();
  			this.tiRigidHorizontal=     gd.getNextNumber();
  			this.tiRigidDiagonal=       gd.getNextNumber();
  			this.tiStrengthOffset=      gd.getNextNumber();
  			this.tiDispScale=           gd.getNextNumber();
  			this.tiDispPow=             gd.getNextNumber();
  			this.tiDispPull=            gd.getNextNumber();
  			this.tiDispPullPreFinal=    gd.getNextNumber();
  			this.tiDispPullFinal=       gd.getNextNumber();
  			this.tiBreakNorm=           gd.getNextNumber();
  			this.tiBreak3=              gd.getNextNumber();
  			this.tiBreak31=             gd.getNextNumber();
  			this.tiBreak21=             gd.getNextNumber();
  			this.tiBreakFar=            gd.getNextNumber();
  			this.tiBreakNear=           gd.getNextNumber();
  			this.tiBreakMode=     (int) gd.getNextNumber();
  			this.tiBreakSame=           gd.getNextNumber();
  			this.tiBreakTurn=           gd.getNextNumber();

  			this.tiHealPreLast=         gd.getNextNumber();
  			this.tiHealLast=            gd.getNextNumber();
  			this.tiHealSame=      (int) gd.getNextNumber();
  			
  			this.tiIterations=    (int) gd.getNextNumber();
  			this.tiPrecision=     (int) gd.getNextNumber();
  			this.tiNumCycles=     (int) gd.getNextNumber();
  			
  			this.stUseRefine=           gd.getNextBoolean();
  			this.stUsePass2=            gd.getNextBoolean();
  			this.stUseRender=           gd.getNextBoolean();

  			this.stShow=                gd.getNextBoolean();
  			this.stSize=          (int) gd.getNextNumber();
  			this.stStepFar=             gd.getNextNumber();
  			this.stStepNear=            gd.getNextNumber();
  			this.stStepThreshold=       gd.getNextNumber();
  			this.stMinDisparity=        gd.getNextNumber();
//  			this.stMaxDisparity=        gd.getNextNumber();
  			this.stFloor=               gd.getNextNumber();
  			this.stPow=                 gd.getNextNumber();
  			this.stSigma=               gd.getNextNumber();
  			this.stMinBgDisparity=      gd.getNextNumber();
  			this.stMinBgFract=          gd.getNextNumber();
  			this.stUseDisp=             gd.getNextNumber();
  			this.stStrengthScale=       gd.getNextNumber();

  			this.stSmplMode=            gd.getNextBoolean();
  			this.stSmplSide=      (int) gd.getNextNumber();
  			this.stSmplNum=       (int) gd.getNextNumber();
  			this.stSmplRms=             gd.getNextNumber();
  			this.stSmplWnd=             gd.getNextBoolean();
  			
  			this.stGrowSel=       (int) gd.getNextNumber();
  			this.stMeasSel=       (int) gd.getNextNumber();
  			this.stSmallDiff=           gd.getNextNumber();
  			this.stHighMix=             gd.getNextNumber();

  			this.outlayerStrength=      gd.getNextNumber();
  			this.outlayerDiff=          gd.getNextNumber();
  			this.outlayerDiffPos=       gd.getNextNumber();
  			this.outlayerDiffNeg=       gd.getNextNumber();

  			this.combine_refine=        gd.getNextBoolean();

  			this.combine_min_strength=  gd.getNextNumber();
  			this.combine_min_hor=       gd.getNextNumber();
  			this.combine_min_vert=      gd.getNextNumber();
//  			this.unique_tolerance=      gd.getNextNumber();
  			this.grow_sweep=      (int) gd.getNextNumber();
  			this.grow_disp_max=         gd.getNextNumber();
  			this.grow_disp_trust=       gd.getNextNumber();
  			this.grow_disp_step=        gd.getNextNumber();
  			this.grow_min_diff=         gd.getNextNumber();
  			this.grow_retry_far=        gd.getNextBoolean();
  			this.grow_pedantic=         gd.getNextBoolean();
  			this.grow_retry_inf=        gd.getNextBoolean();
  			
  			this.gr_new_expand=         gd.getNextBoolean();
  			this.gr_max_expand=   (int) gd.getNextNumber();
  			this.fds_str_floor=         gd.getNextNumber();
  			this.gr_ovrbg_cmb=          gd.getNextNumber();
  			this.gr_ovrbg_cmb_hor=      gd.getNextNumber();
  			this.gr_ovrbg_cmb_vert=     gd.getNextNumber();
  			this.gr_ovrbg_filtered=     gd.getNextNumber();
  			this.fds_str_pow=           gd.getNextNumber();
  			this.fds_smpl_side=   (int) gd.getNextNumber();
  			this.fds_smpl_num=    (int) gd.getNextNumber();
  			this.fds_smpl_rms=          gd.getNextNumber();
  			this.fds_smpl_rel_rms=      gd.getNextNumber();
  			this.fds_smpl_wnd=          gd.getNextBoolean();
  			this.fds_abs_tilt=          gd.getNextNumber();
  			this.fds_rel_tilt=          gd.getNextNumber();

  			this.mc_disp8_step=         gd.getNextNumber();
  			this.mc_disp8_trust=        gd.getNextNumber();
  			this.mc_strength=           gd.getNextNumber();
  			this.mc_unique_tol=         gd.getNextNumber();
  			
  			this.mc_unique_tol=         gd.getNextNumber();
  			this.mc_unique_tol=         gd.getNextNumber();
  			this.mc_unique_tol=         gd.getNextNumber();
  			this.mc_unique_tol=         gd.getNextNumber();
  			this.mc_unique_tol=         gd.getNextNumber();

  			this.gr_min_new=      (int) gd.getNextNumber();
  			this.gr_var_new_sngl=       gd.getNextBoolean();
  			this.gr_var_new_fg=         gd.getNextBoolean();
  			this.gr_var_all_fg=         gd.getNextBoolean();
  			this.gr_var_new_bg=         gd.getNextBoolean();
  			this.gr_var_all_bg=         gd.getNextBoolean();
  			this.gr_var_next=           gd.getNextBoolean();
  			this.gr_num_steps=    (int) gd.getNextNumber();
  			this.gr_steps_over=   (int) gd.getNextNumber();
  			this.gr_smpl_size=    (int) gd.getNextNumber();
  			this.gr_min_pnts=     (int) gd.getNextNumber();
  			this.gr_use_wnd=            gd.getNextBoolean();
  			this.gr_tilt_damp=          gd.getNextNumber();
  			this.gr_split_rng=          gd.getNextNumber();
  			this.gr_same_rng=           gd.getNextNumber();
  			this.gr_diff_cont=          gd.getNextNumber();
  			this.gr_abs_tilt=           gd.getNextNumber();
  			this.gr_rel_tilt=           gd.getNextNumber();
  			this.gr_smooth=       (int) gd.getNextNumber();
  			this.gr_fin_diff=           gd.getNextNumber();
  			this.gr_unique_tol=         gd.getNextNumber();
  			this.gr_unique_pretol=      gd.getNextNumber();

  			this.plPreferDisparity=     gd.getNextBoolean();
  			this.plDispNorm=            gd.getNextNumber();

  			this.plBlurBinVert=         gd.getNextNumber();
  			this.plBlurBinHor=          gd.getNextNumber();
  			this.plMaxDiffVert=         gd.getNextNumber();
  			this.plMaxDiffHor=          gd.getNextNumber();
  			this.plInitPasses=    (int) gd.getNextNumber();
  			
  			this.plMinPoints=     (int) gd.getNextNumber();
  			this.plTargetEigen=         gd.getNextNumber();
  			this.plFractOutliers=       gd.getNextNumber();
  			this.plMaxOutliers=   (int) gd.getNextNumber();
  			this.plMinStrength=         gd.getNextNumber();
  			this.plMaxEigen=            gd.getNextNumber();
  			this.plEigenFloor=          gd.getNextNumber();
  			this.plEigenStick=          gd.getNextNumber();
  			this.plBadPlate=            gd.getNextNumber();
  			this.plDbgMerge=            gd.getNextBoolean();
  			this.plWorstWorsening=      gd.getNextNumber();
  			this.plWorstWorsening2=     gd.getNextNumber();
  			this.plWorstEq=             gd.getNextNumber();
  			this.plWorstEq2=            gd.getNextNumber();
  			this.plOKMergeEigen=        gd.getNextNumber();
  			this.plMaxWorldSin2=        gd.getNextNumber();
  			this.pl2dForSin=            gd.getNextNumber();
  			this.plWeakWorsening=       gd.getNextNumber();
  			this.plMaxOverlap=          gd.getNextNumber();

  			this.plWeakWeight=          gd.getNextNumber();
  			this.plWeakEigen=           gd.getNextNumber();
  			this.plWeakWeight2=         gd.getNextNumber();
  			this.plWeakEigen2=          gd.getNextNumber();
  			this.plSumThick=            gd.getNextNumber();
  			this.plNeNeibCost=          gd.getNextNumber();
  			this.plNeOwn=               gd.getNextNumber();

  			this.plExNeibCost=          gd.getNextNumber();
  			this.plExNeibSmooth=        gd.getNextNumber();
  			this.plMergeCostStar=       gd.getNextNumber();
  			this.plMergeCost=           gd.getNextNumber();

  			this.plConflMerge=          gd.getNextBoolean();
  			this.plConflRelax=          gd.getNextNumber();
  			this.plConflSngl=           gd.getNextBoolean();
  			this.plConflSnglPair=       gd.getNextBoolean();

  			this.plWeakFgStrength=      gd.getNextNumber();
  			this.plWeakFgOutliers=(int) gd.getNextNumber();
  			this.plWeakFgRelax=         gd.getNextNumber();
  			
  			this.plThickWorld=          gd.getNextNumber();
  			this.plThickWorldConfl=     gd.getNextNumber();
  			this.plRelaxComplete=       gd.getNextNumber();
  			this.plRelaxComplete2=      gd.getNextNumber();

  			this.plMaxZRatio=           gd.getNextNumber();
  			this.plMaxDisp=             gd.getNextNumber();
  			this.plCutTail=             gd.getNextNumber();
  			this.plMinTail=             gd.getNextNumber();

  			this.plDiscrEn=             gd.getNextBoolean();
  			this.plDiscrTolerance=      gd.getNextNumber();
  			this.plDiscrDispRange=      gd.getNextNumber();
  			this.plDiscrSteps=    (int) gd.getNextNumber();
  			this.plDiscrVariants= (int) gd.getNextNumber();
  			this.plDiscrMode=     (int) gd.getNextNumber();
  			this.plDiscrVarFloor=       gd.getNextNumber();
  			this.plDiscrSigma=          gd.getNextNumber();
  			this.plDiscrBlur=           gd.getNextNumber();
  			this.plDiscrExclusivity=    gd.getNextNumber();
  			this.plDiscrExclus2=        gd.getNextNumber();
  			this.plDiscrStrict=         gd.getNextBoolean();
  			this.plDiscrCorrMax=        gd.getNextNumber();
  			this.plDiscrCorrMerge=      gd.getNextNumber();
  			this.plDiscrSteal=    (int) gd.getNextNumber();
  			this.plDiscrGrown=    (int) gd.getNextNumber();
  			this.plDiscrXMedian=        gd.getNextNumber();

  			this.plCostDist=            gd.getNextNumber();
  			this.plCostKrq=             gd.getNextNumber();
  			this.plCostKrqEq=           gd.getNextNumber();
  			this.plCostWrq=             gd.getNextNumber();
  			this.plCostWrqEq=           gd.getNextNumber();
  			this.plCostSin2=            gd.getNextNumber();
  			this.plCostRdist2=          gd.getNextNumber();

  			this.plConflDualTri=        gd.getNextBoolean();
  			this.plConflMulti=          gd.getNextBoolean();
  			this.plConflDiag=           gd.getNextBoolean();
  			this.plConflStar=           gd.getNextBoolean();
  			this.plStarSteps=     (int) gd.getNextNumber();
  			this.plStarOrtho=           gd.getNextNumber();
  			this.plStarDiag=            gd.getNextNumber();
  			this.plStarPwr=             gd.getNextNumber();
  			this.plStarWeightPwr=       gd.getNextNumber();
  			this.plWeightToDens=        gd.getNextNumber();
  			this.plStarValPwr=          gd.getNextNumber();
  			this.plDblTriLoss=          gd.getNextNumber();
  			this.plNewConfl=            gd.getNextBoolean();
  			this.plMaxChanges=    (int) gd.getNextNumber();

  			this.plMutualOnly=          gd.getNextBoolean();

  			this.plFillSquares=         gd.getNextBoolean();
  			this.plCutCorners=          gd.getNextBoolean();
  			this.plHypotenuse=          gd.getNextBoolean();

  			this.plPull=                gd.getNextNumber();
  			this.plNormPow=             gd.getNextNumber();
  			this.plIterations=    (int) gd.getNextNumber();
  			this.plStopBad=             gd.getNextBoolean();
  			this.plPrecision=     (int) gd.getNextNumber();

  			this.plSplitPull=           gd.getNextNumber();
  			this.plSplitMinNeib=  (int) gd.getNextNumber();
  			this.plSplitMinWeight=      gd.getNextNumber();
  			this.plSplitMinQuality=     gd.getNextNumber();
  			this.plSplitApply=          gd.getNextBoolean();
  			this.plNonExclusive=        gd.getNextBoolean();
  			this.plUseOtherPlanes=      gd.getNextBoolean();
  			this.plAllowParallel=       gd.getNextBoolean();
  			this.plMaxDiff=             gd.getNextNumber();
  			this.plOtherDiff=           gd.getNextNumber();
  			this.plSplitXY=             gd.getNextBoolean();
  			this.plSplitXYTolerance=    gd.getNextNumber();

  			this.plFuse=                gd.getNextBoolean();
  			this.plKeepOrphans=         gd.getNextBoolean();
  			this.plMinOrphan=           gd.getNextNumber();

  			this.plSnapDispAny=         gd.getNextNumber();
  			this.plSnapStrengthAny=     gd.getNextNumber();
  			this.plSnapNegAny=          gd.getNextNumber();
  			this.plSnapDispMax=         gd.getNextNumber();
  			this.plSnapDispWeight=      gd.getNextNumber();
  			this.plSnapZeroMode=  (int) gd.getNextNumber();

  			this.msUseSel=              gd.getNextBoolean();
  			this.msDivideByArea=        gd.getNextBoolean();
  			this.msScaleProj=           gd.getNextNumber();
  			this.msFractUni=            gd.getNextNumber();

  			this.tsNoEdge=              gd.getNextBoolean();
  			this.tsUseCenter=           gd.getNextBoolean();
  			this.tsMaxDiff=             gd.getNextNumber();
  			this.tsMinDiffOther=        gd.getNextNumber();
  			this.tsMinStrength=         gd.getNextNumber();
  			this.tsMaxStrength=         gd.getNextNumber();
  			this.tsMinSurface=          gd.getNextNumber();
  			this.tsMoveDirs=      (int) gd.getNextNumber();
  			this.tsSurfStrPow=          gd.getNextNumber();
  			this.tsAddStrength=         gd.getNextNumber();
  			this.tsSigma=               gd.getNextNumber();
  			this.tsNSigma=              gd.getNextNumber();
  			this.tsMinPull=             gd.getNextNumber();
  			this.tsMinAdvantage=        gd.getNextNumber();
  			
  			this.tsClustSize=     (int) gd.getNextNumber();
  			this.tsClustWeight=         gd.getNextNumber();
  			this.tsMinNeib=       (int) gd.getNextNumber();
  			this.tsMaxSurStrength=      gd.getNextNumber();
  			this.tsCountDis=            gd.getNextBoolean();
  			this.tsReset              = false;  // Reset tiles to surfaces assignment

  			this.tsEnPlaneSeed=         gd.getNextBoolean();
  			this.tsEnOnly=              gd.getNextBoolean();
  			this.tsEnGrow=              gd.getNextBoolean();
  			this.tsGrowStrength=        gd.getNextNumber();
  			this.tsGrowStrong=          gd.getNextBoolean();
  			this.tsContStrength=        gd.getNextNumber();
  			this.tsContDiff=            gd.getNextNumber();

  			this.tsEnSingle=            gd.getNextBoolean();
  			this.tsEnMulti=             gd.getNextBoolean();
  			this.tsRemoveWeak1=         gd.getNextBoolean();
  			this.tsGrowSurround=        gd.getNextBoolean();
  			this.tsRemoveWeak2=         gd.getNextBoolean();
  			
  			this.tsLoopMulti=           gd.getNextBoolean();
  			this.tsShow               = gd.getNextBoolean();
  			this.tsNumClust =     (int) gd.getNextNumber();

  			this.tsConsensMode =  (int) gd.getNextNumber();
  			this.tsConsensAgree = (int) gd.getNextNumber();

  			this.taMinFgBg=             gd.getNextNumber();
  			this.taMinFgEdge=           gd.getNextNumber();
  			this.taMinColSep=           gd.getNextNumber();
  			this.taMinColDiff=          gd.getNextNumber();
  			this.taOutlier=             gd.getNextNumber();
  			this.taDiffPwr=             gd.getNextNumber();
  			this.taBestPwr=             gd.getNextNumber();
  			this.taDiff9Pwr=            gd.getNextNumber();
  			this.taColSigma=            gd.getNextNumber();
  			this.taColFraction=         gd.getNextNumber();
  			

  			this.taCostEmpty=           gd.getNextNumber();
  			this.taCostNoLink=          gd.getNextNumber();
  			this.taCostSwitch=          gd.getNextNumber();
  			this.taCostColor=           gd.getNextNumber();
  			this.taCostDiff=            gd.getNextNumber();
  			this.taCostDiffBest=        gd.getNextNumber();
  			this.taCostDiff9=           gd.getNextNumber();
  			this.taCostWeakFgnd=        gd.getNextNumber();
  			this.taCostFlaps=           gd.getNextNumber();
  			this.taCostMismatch=        gd.getNextNumber();

  			this.taEnEmpty=             gd.getNextBoolean();
  			this.taEnNoLink=            gd.getNextBoolean();
  			this.taEnSwitch=            gd.getNextBoolean();
  			this.taEnColor=             gd.getNextBoolean();
  			this.taEnDiff=              gd.getNextBoolean();
  			this.taEnDiffBest=          gd.getNextBoolean();
  			this.taEnDiff9=             gd.getNextBoolean();
  			this.taEnWeakFgnd=          gd.getNextBoolean();
  			this.taEnFlaps=             gd.getNextBoolean();
  			this.taEnMismatch=          gd.getNextBoolean();
  			
  			this.dbg_migrate=           gd.getNextBoolean();

  			this.show_ortho_combine=    gd.getNextBoolean();
  			this.show_refine_supertiles=gd.getNextBoolean();
  			this.show_bgnd_nonbgnd=     gd.getNextBoolean(); // first on second pass
  			this.show_filter_scan=      gd.getNextBoolean(); // first on refine
  			this.show_combined=         gd.getNextBoolean();
  			this.show_unique=           gd.getNextBoolean();
  			this.show_histograms=       gd.getNextBoolean();
  			this.show_init_refine=      gd.getNextBoolean();
  			this.show_expand=           gd.getNextBoolean();
  			this.show_variant=          gd.getNextBoolean();
  			this.show_retry_far=        gd.getNextBoolean();
  			this.show_macro=            gd.getNextBoolean();
  			this.show_shells=           gd.getNextBoolean();
  			this.show_neighbors=        gd.getNextBoolean();
  			this.show_flaps_dirs=       gd.getNextBoolean();
  			this.show_first_clusters=   gd.getNextBoolean();
  			this.show_planes=           gd.getNextBoolean();
  			
  			return true;
  		}
  		
  		public boolean showTsDialog() {
  			GenericDialog gd = new GenericDialog("Set CLT tiles to surfaces assignment parameters");
  			gd.addCheckbox    ("Do not assign tiles to the surface edges (not having all 8 neighbors)",           this.tsNoEdge);
  			gd.addCheckbox    ("Only assign outside of 8x8 center if no suitable alternative",                    this.tsUseCenter);
  			gd.addNumericField("Maximal disparity difference when assigning tiles",                               this.tsMaxDiff,  6);
  			gd.addNumericField("Minimal disparity difference to be considered as a competitor surface",           this.tsMinDiffOther,  6);
  			gd.addNumericField("Minimal tile correlation strength to be assigned",                                this.tsMinStrength,  6);
  			gd.addNumericField("Maximal tile correlation strength to be assigned",                                this.tsMaxStrength,  6);
  			gd.addNumericField("Minimal surface strength at the tile location",                                   this.tsMinSurface,  6);
  			gd.addNumericField("Allowed tile disparity correction: 1 increase, 2 - decrease, 3 - both directions",this.tsMoveDirs,  0);
  			gd.addNumericField("Raise surface strengths ratio to this power when comparing candidates",           this.tsSurfStrPow,  6);
  			gd.addNumericField("Add to strengths when calculating pull of assigned tiles",                        this.tsAddStrength,  6);
  			gd.addNumericField("Radius of influence (in tiles) of the previously assigned tiles",                 this.tsSigma,  6);
  			gd.addNumericField("Maximal relative to radius distance to calculate influence",                      this.tsNSigma,  6);
  			gd.addNumericField("Additional pull of each surface ",                                                this.tsMinPull,  6);
  			gd.addNumericField("Minimal ratio of the best surface candidate to the next one to make selection",   this.tsMinAdvantage,  6);
  			gd.addNumericField("Minimal size of a cluster to keep",                                               this.tsClustSize,  0);
  			gd.addNumericField("Minimal total weight of a cluster to keep",                                       this.tsClustWeight,  6);
  			gd.addNumericField("Minimal number of neighbors of unassigned tile to join (the farthest)",           this.tsMinNeib,  0);
  			gd.addNumericField("Maximal strength of the surrounded unassigned tile to join",                      this.tsMaxSurStrength,  6);
  			gd.addCheckbox    ("Include disabled tiles/borders when counting assigned neighbors",                 this.tsCountDis);
  			gd.addCheckbox    ("Reset tiles to surfaces assignment",                                              false);
  			
  			gd.addCheckbox    ("Assign tiles that were used to generate planes",                                  this.tsEnPlaneSeed);
  			gd.addCheckbox    ("Allow assignment only surface",                                                   this.tsEnOnly);
  			gd.addCheckbox    ("Grow the only surface assignments",                                               this.tsEnGrow);
  			gd.addNumericField("Maximal strength when growing the only surfaces",                                 this.tsGrowStrength,  6);
  			gd.addCheckbox    ("Grow over strong if disparity matches",                                           this.tsGrowStrong);
  			gd.addNumericField("Minimal strength to continue grow with disparity match",                          this.tsContStrength,  6);
  			gd.addNumericField("Maximal normalized disparity error to grow over strong tiles",                    this.tsContDiff,  6);

  			gd.addCheckbox    ("Allow assignment to the nearest surface with no competitors",                     this.tsEnSingle);
  			gd.addCheckbox    ("Allow assignment when several surfaces fit",                                      this.tsEnMulti);
  			gd.addCheckbox    ("Remove weak clusters before growing",                                             this.tsRemoveWeak1);
  			gd.addCheckbox    ("Assign tiles that have neighbors to the lowest disparity",                        this.tsGrowSurround);
  			gd.addCheckbox    ("Remove weak clusters after growing",                                              this.tsRemoveWeak2);
  			
  			gd.addCheckbox    ("Repeat multi-choice assignment while succeeding",                                 this.tsLoopMulti);
  			gd.addCheckbox    ("Show results of tiles to surfaces assignment",                                    this.tsShow);
  			gd.addNumericField("Number of clusters to keep",                                                      this.tsNumClust,  0);

  			gd.addNumericField("Which assignments to match +1 - combo, +2 grown single, +4 plane seeds",          this.tsConsensMode,  0);
  			gd.addNumericField("Minimal number of assignments to agree",                                          this.tsConsensAgree,  0);
  			
  			gd.addMessage     ("--- Tile assignment parameters ---");
  			gd.addNumericField("Minimal foreground/ background separation to look for weak FG edge",              this.taMinFgBg,    6);
  			gd.addNumericField("Minimal foreground edge strength (stronger edges will have proportionally smaller costs)", this.taMinFgEdge,  6);
  			gd.addNumericField("Minimal surface separation that requires color change",                           this.taMinColSep,  6);
  			gd.addNumericField("Minimal color variation (larger proportionally reduces cost)",                    this.taMinColDiff, 6);
  	  		gd.addNumericField("Disparity difference limit (to handle outliers)",                                 this.taOutlier, 6);
  	  		gd.addNumericField("Strength power when calculating disparity error",                                 this.taDiffPwr, 6);
  	  		gd.addNumericField("Strength power when calculating disparity error over best",                       this.taBestPwr, 6);
  	  		gd.addNumericField("Strength power when calculating disparity error for group of 9",                  this.taDiff9Pwr, 6);
  			gd.addNumericField("Gaussian sigma to blur color difference between tiles along each direction",      this.taColSigma, 6);
  			gd.addNumericField("Relative amount of the blurred color difference in the mixture",                  this.taColFraction, 6);

  			gd.addNumericField("Cost of a tile that is not assigned",                                             this.taCostEmpty,  6);
  			gd.addNumericField("Cost of a tile not having any neighbor in particular direction",                  this.taCostNoLink,  6);
  			gd.addNumericField("Cost of a tile switching to a neighbor that does not have a link",                this.taCostSwitch,  6);
  			gd.addNumericField("Cost of a tile switching to a disconnected neighbor divided by a color",          this.taCostColor,  6);
  			gd.addNumericField("Cost of a weighted normalized tile disparity error",                              this.taCostDiff,  6);
  			gd.addNumericField("Cost of a weighted normalized tile disparity error above best surface",           this.taCostDiffBest,  6);
  			gd.addNumericField("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",this.taCostDiff9,  6);
  			gd.addNumericField("Cost of a weak foreground edge",                                                  this.taCostWeakFgnd,  6);
  			gd.addNumericField("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",            this.taCostFlaps,  6);
  			gd.addNumericField("Cost of a measurement layer not having same layer in the same location or near",  this.taCostMismatch,  6);

  			gd.addCheckbox    ("Cost of a tile that is not assigned",                                             this.taEnEmpty);
  			gd.addCheckbox    ("Cost of a tile not having any neighbor in particular direction",                  this.taEnNoLink);
  			gd.addCheckbox    ("Cost of a tile switching to a neighbor that does not have a link",                this.taEnSwitch);
  			gd.addCheckbox    ("Cost of a tile switching to a disconnected neighbor divided by a color",          this.taEnColor);
  			gd.addCheckbox    ("Cost of a weighted normalized tile disparity error",                              this.taEnDiff);
  			gd.addCheckbox    ("Cost of a weighted normalized tile disparity error above best surface",           this.taEnDiffBest);
  			gd.addCheckbox    ("Cost of a weighted normalized tile disparity error for tile and 8 neighbors (DC)",this.taEnDiff9);
  			gd.addCheckbox    ("Cost of a weak foreground edge",                                                  this.taEnWeakFgnd);
  			gd.addCheckbox    ("Cost of using supertile \"flaps\" (not in the center 8x8 tiles area)",            this.taEnFlaps);
  			gd.addCheckbox    ("Cost of a measurement layer not having same layer in the same location or near",  this.taEnMismatch);
  			
  			WindowTools.addScrollBars(gd);
  			gd.showDialog();
  			if (gd.wasCanceled()) return false;
  			this.tsNoEdge=              gd.getNextBoolean();
  			this.tsUseCenter=           gd.getNextBoolean();
  			this.tsMaxDiff=             gd.getNextNumber();
  			this.tsMinDiffOther=        gd.getNextNumber();
  			this.tsMinStrength=         gd.getNextNumber();
  			this.tsMaxStrength=         gd.getNextNumber();
  			this.tsMinSurface=          gd.getNextNumber();
  			this.tsMoveDirs=      (int) gd.getNextNumber();
  			this.tsSurfStrPow=          gd.getNextNumber();
  			this.tsAddStrength=         gd.getNextNumber();
  			this.tsSigma=               gd.getNextNumber();
  			this.tsNSigma=              gd.getNextNumber();
  			this.tsMinPull=             gd.getNextNumber();
  			this.tsMinAdvantage=        gd.getNextNumber();

  			this.tsClustSize=     (int) gd.getNextNumber();
  			this.tsClustWeight=         gd.getNextNumber();
  			this.tsMinNeib=       (int) gd.getNextNumber();
  			this.tsMaxSurStrength=      gd.getNextNumber();
  			this.tsCountDis=            gd.getNextBoolean();
  			this.tsReset=               gd.getNextBoolean();

  			this.tsEnPlaneSeed=         gd.getNextBoolean();
  			this.tsEnOnly=              gd.getNextBoolean();
  			this.tsEnGrow=              gd.getNextBoolean();
  			this.tsGrowStrength=        gd.getNextNumber();
  			this.tsGrowStrong=          gd.getNextBoolean();
  			this.tsContStrength=        gd.getNextNumber();
  			this.tsContDiff=            gd.getNextNumber();
 			
  			this.tsEnSingle=            gd.getNextBoolean();
  			this.tsEnMulti=             gd.getNextBoolean();
  			this.tsRemoveWeak1=         gd.getNextBoolean();
  			this.tsGrowSurround=        gd.getNextBoolean();
  			this.tsRemoveWeak2=         gd.getNextBoolean();
  			
  			this.tsLoopMulti=           gd.getNextBoolean();
  			this.tsShow =               gd.getNextBoolean();
  			this.tsNumClust=      (int) gd.getNextNumber();
  			this.tsConsensMode =  (int) gd.getNextNumber();
  			this.tsConsensAgree = (int) gd.getNextNumber();

  			this.taMinFgBg=             gd.getNextNumber();
  			this.taMinFgEdge=           gd.getNextNumber();
  			this.taMinColSep=           gd.getNextNumber();
  			this.taMinColDiff=          gd.getNextNumber();
  			this.taOutlier=             gd.getNextNumber();
  			this.taDiffPwr=             gd.getNextNumber();
  			this.taBestPwr=             gd.getNextNumber();
  			this.taDiff9Pwr=            gd.getNextNumber();
  			this.taColSigma=            gd.getNextNumber();
  			this.taColFraction=         gd.getNextNumber();

  			this.taCostEmpty=           gd.getNextNumber();
  			this.taCostNoLink=          gd.getNextNumber();
  			this.taCostSwitch=          gd.getNextNumber();
  			this.taCostColor=           gd.getNextNumber();
  			this.taCostDiff=            gd.getNextNumber();
  			this.taCostDiffBest=        gd.getNextNumber();
  			this.taCostDiff9=           gd.getNextNumber();
  			this.taCostWeakFgnd=        gd.getNextNumber();
  			this.taCostFlaps=           gd.getNextNumber();
  			this.taCostMismatch=        gd.getNextNumber();

  			this.taEnEmpty=             gd.getNextBoolean();
  			this.taEnNoLink=            gd.getNextBoolean();
  			this.taEnSwitch=            gd.getNextBoolean();
  			this.taEnColor=             gd.getNextBoolean();
  			this.taEnDiff=              gd.getNextBoolean();
  			this.taEnDiffBest=          gd.getNextBoolean();
  			this.taEnDiff9=             gd.getNextBoolean();
  			this.taEnWeakFgnd=          gd.getNextBoolean();
  			this.taEnFlaps=             gd.getNextBoolean();
  			this.taEnMismatch=          gd.getNextBoolean();

  			return true;
  		}
    }

    public static class DCTParameters {
  		public int dct_size =            8; //
  		public int asym_size =          15; //
  		public int asym_pixels =         4; // maximal number of non-zero pixels in direct convolution kernel
  		public int asym_distance =       1; // how far to try a new asym kernel pixel from existing ones 
  		public int dct_window =          1; // currently only 3 types of windows - 0 (none), 1 and 2
  		public int LMA_steps =         100;
  		public double fact_precision=    0.003; // stop iterations if error rms less than this part of target kernel rms
  		public double compactness =      0.0;
  		public double sym_compactness =  0.5;
  		public double dc_weight =        1.0;  // importance of dc realtive to rms_pure
  		public int    asym_tax_free  =   0;    // "compactness" does not apply to pixels with |x|<=asym_tax_free  and |y| <= asym_tax_free
  		public int    seed_size =        1;    // number of initial cells in asym_kernel - should be 4*b + 1 (X around center cell) or 4*n + 0  (X around between cells)
  		public double asym_random  =    -1; // initialize asym_kernel with random numbers
  		public double dbg_x =            2.7;
  		public double dbg_y =            0.0;
  		public double dbg_x1 =          -1.3;
  		public double dbg_y1 =           2.0;
  		public double dbg_sigma =        0.8;
  		public double dbg_src_size =     8.0; // trying to slightly scale in dct space. == dct = 1:1, dct+1.0 - shrink dct(dct+1.0)
  		public double dbg_scale =        1.0; // Should ==DCT_PARAMETERS.dct_size / DCT_PARAMETERS.dbg_src_size
  		public double dbg_fold_scale =   1.0; // Modifies window during MDCT->DCT-IV folding
  		public String dbg_mask = ".........:::::::::.........:::::::::......*..:::::*:::.........:::::::::.........";
  		public int dbg_mode =            1; // 0 - old LMA, 1 - new LMA - *** not used anymore ***
  		public int dbg_window_mode =     1; // 0 - none, 1 - square, 2 - sin 3 - sin^2 Now _should_ be square !!!
  		public boolean centerWindowToTarget = true;
  		// parameters to extract a kernel from the kernel image file
  		public int    color_channel =    2; // green (<0 - use simulated kernel, also will use simulated if kernels are not set)
  		public int    decimation =       2; // decimate original kernel this much in each direction
  		public double decimateSigma =   -1.0; // special mode for 2:1 deciamtion 
  		public int    tileX =            82;  // number of kernel tile (0..163) 
  		public int    tileY =            62;  // number of kernel tile (0..122) 
  		public int    kernel_chn =      -1; // camera channel calibration to use for aberration correction ( < 0 - no correction)
  		public boolean normalize =       true; // normalize both sym and asym kernels (asym to have sum==1, sym to have sum = dct_size
  		public boolean normalize_sym =   true; // normalize sym kernels separately
  		public boolean antiwindow =      false; // divide symmetrical kernel by a window function
  		public boolean skip_sym =        false; // do not apply symmetrical correction
  		public boolean convolve_direct = false; // do not apply symmetrical correction
  		
  		// colors should be balanced before DCT color conversion!
  		public double novignetting_r    = 0.2644; // reg gain in the center of sensor calibration R (instead of vignetting)
  		public double novignetting_g    = 0.3733; // green gain in the center of sensor calibration G
  		public double novignetting_b    = 0.2034; // blue gain in the center of sensor calibration B
  		
  		public double scale_r =           1.0; // extra gain correction after vignetting or nonvignetting, before other processing
  		public double scale_g =           1.0;
  		public double scale_b =           1.0;
  		
  		public double vignetting_max    = 0.4; // value in vignetting data to correspond to 1x in the kernel
  		public double vignetting_range  = 5.0; // do not try to correct vignetting less than vignetting_max/vignetting_range
  		
  		public boolean post_debayer     = false; // perform de-bayer after aberrations in pixel domain
  		public boolean color_DCT        = true; // false - use old color processing mode
  		public double  sigma_rb =         0.9; // additional (to G) blur for R and B
  		public double  sigma_y =          0.7; // blur for G contribution to Y
  		public double  sigma_color =      0.7; // blur for Pb and Pr in addition to that of Y
  		public double  line_thershold =   1.0; // line detection amplitude to apply line enhancement - not used?
  		public boolean nonlin =           true; // enable nonlinear processing (including denoise
  		public double  nonlin_max_y =     1.0; // maximal amount of nonlinear line/edge emphasis for Y component
  		public double  nonlin_max_c =     1.0; // maximal amount of nonlinear line/edge emphasis for C component
  		public double  nonlin_y =         0.1; // amount of nonlinear line/edge emphasis for Y component
  		public double  nonlin_c =         0.01; // amount of nonlinear line/edge emphasis for C component
  		public double  nonlin_corn =      0.5;  // relative weight for nonlinear corner elements
  		public boolean denoise =          true; // suppress noise during nonlinear processing
  		public double  denoise_y =        1.0;  // maximal total smoothing of the Y post-kernel (will compete with edge emphasis)
  		public double  denoise_c =        1.0;  // maximal total smoothing of the color differences post-kernel (will compete with edge emphasis)
  		public double  denoise_y_corn =   0.3;  // weight of the 4 corner pixels during denoise y (straight - 1-denoise_y_corn)
  		public double  denoise_c_corn =   0.3;  // weight of the 4 corner pixels during denoise y (straight - 1-denoise_c_corn)
  		

  		public DCTParameters(){}
  		
  		public DCTParameters(
  				int dct_size,
  				int asym_size,
  				int asym_pixels,
  				int asym_distance,
  				int dct_window,
  				double compactness,
  				int asym_tax_free,
  				int seed_size) {
  			this.dct_size =       dct_size;
  			this.asym_size =      asym_size;
  			this.asym_pixels =    asym_pixels;
  			this.asym_distance =  asym_distance;
  			this.dct_window =     dct_window;
  			this.compactness =    compactness;
  			this.asym_tax_free =  asym_tax_free;
  			this.seed_size =      seed_size;
  		}
  		public void setProperties(String prefix,Properties properties){
  			properties.setProperty(prefix+"dct_size",this.dct_size+"");
  			properties.setProperty(prefix+"asym_size",this.asym_size+"");
  			properties.setProperty(prefix+"asym_pixels",this.asym_pixels+"");
  			properties.setProperty(prefix+"asym_distance",this.asym_distance+"");
  			properties.setProperty(prefix+"dct_window",   this.dct_window+"");
  			properties.setProperty(prefix+"compactness",  this.compactness+"");
  			properties.setProperty(prefix+"sym_compactness",  this.sym_compactness+"");
  			properties.setProperty(prefix+"dc_weight",  this.dc_weight+"");
  			properties.setProperty(prefix+"fact_precision",  this.fact_precision+"");
  			properties.setProperty(prefix+"asym_tax_free",  this.asym_tax_free+"");
  			properties.setProperty(prefix+"seed_size",  this.seed_size+"");
  			properties.setProperty(prefix+"asym_random",  this.asym_random+"");
  			properties.setProperty(prefix+"LMA_steps",  this.LMA_steps+"");
  			properties.setProperty(prefix+"dbg_x",      this.dbg_x+"");
  			properties.setProperty(prefix+"dbg_y",      this.dbg_y+"");
  			properties.setProperty(prefix+"dbg_x1",     this.dbg_x1+"");
  			properties.setProperty(prefix+"dbg_y1",     this.dbg_y1+"");
  			properties.setProperty(prefix+"dbg_sigma",  this.dbg_sigma+"");
  			properties.setProperty(prefix+"dbg_src_size",this.dbg_src_size+"");
  			properties.setProperty(prefix+"dbg_scale",  this.dbg_scale+"");
  			properties.setProperty(prefix+"dbg_fold_scale",  this.dbg_fold_scale+"");
  			properties.setProperty(prefix+"dbg_mask",   this.dbg_mask+"");
  			properties.setProperty(prefix+"dbg_mode",   this.dbg_mode+"");
  			properties.setProperty(prefix+"dbg_window_mode",    this.dbg_window_mode+"");
  			properties.setProperty(prefix+"centerWindowToTarget",   this.centerWindowToTarget+"");
  			properties.setProperty(prefix+"color_channel",      this.color_channel+"");
  			properties.setProperty(prefix+"decimation",         this.decimation+"");
  			properties.setProperty(prefix+"decimateSigma",      this.decimateSigma+"");
  			properties.setProperty(prefix+"tileX",              this.tileX+"");
  			properties.setProperty(prefix+"tileY",              this.tileY+"");
  			properties.setProperty(prefix+"kernel_chn",         this.kernel_chn+"");
  			properties.setProperty(prefix+"normalize",          this.normalize+"");
  			properties.setProperty(prefix+"normalize_sym",      this.normalize_sym+"");
  			properties.setProperty(prefix+"antiwindow",         this.antiwindow+"");
  			properties.setProperty(prefix+"skip_sym",           this.skip_sym+"");
  			properties.setProperty(prefix+"convolve_direct",    this.convolve_direct+"");
  			properties.setProperty(prefix+"novignetting_r",     this.novignetting_r+"");
  			properties.setProperty(prefix+"novignetting_g",     this.novignetting_g+"");
  			properties.setProperty(prefix+"novignetting_b",     this.novignetting_b+"");
  			properties.setProperty(prefix+"scale_r",            this.scale_r+"");
  			properties.setProperty(prefix+"scale_g",            this.scale_g+"");
  			properties.setProperty(prefix+"scale_b",            this.scale_b+"");
  			properties.setProperty(prefix+"vignetting_max",     this.vignetting_max+"");
  			properties.setProperty(prefix+"vignetting_range",   this.vignetting_range+"");
  			properties.setProperty(prefix+"post_debayer",       this.post_debayer+"");
  			properties.setProperty(prefix+"color_DCT",          this.color_DCT+"");
  			properties.setProperty(prefix+"sigma_rb",           this.sigma_rb+"");
  			properties.setProperty(prefix+"sigma_y",            this.sigma_y+"");
  			properties.setProperty(prefix+"sigma_color",        this.sigma_color+"");
  			properties.setProperty(prefix+"line_thershold",     this.line_thershold+"");
  			properties.setProperty(prefix+"nonlin",             this.nonlin+"");
  			properties.setProperty(prefix+"nonlin_max_y",       this.nonlin_max_y+"");
  			properties.setProperty(prefix+"nonlin_max_c",       this.nonlin_max_c+"");
  			properties.setProperty(prefix+"nonlin_y",           this.nonlin_y+"");
  			properties.setProperty(prefix+"nonlin_c",           this.nonlin_c+"");
  			properties.setProperty(prefix+"nonlin_corn",        this.nonlin_corn+"");
  			properties.setProperty(prefix+"denoise",            this.denoise+"");
  			properties.setProperty(prefix+"denoise_y",          this.denoise_y+"");
  			properties.setProperty(prefix+"denoise_c",          this.denoise_c+"");
  			properties.setProperty(prefix+"denoise_y_corn",     this.denoise_y_corn+"");
  			properties.setProperty(prefix+"denoise_c_corn",     this.denoise_c_corn+"");
  			
  			
  			
  		}
  		public void getProperties(String prefix,Properties properties){
  			if (properties.getProperty(prefix+"dct_size")!=null) this.dct_size=Integer.parseInt(properties.getProperty(prefix+"dct_size"));
  			if (properties.getProperty(prefix+"asym_size")!=null) this.asym_size=Integer.parseInt(properties.getProperty(prefix+"asym_size"));
  			if (properties.getProperty(prefix+"asym_pixels")!=null) this.asym_pixels=Integer.parseInt(properties.getProperty(prefix+"asym_pixels"));
  			if (properties.getProperty(prefix+"asym_distance")!=null) this.asym_distance=Integer.parseInt(properties.getProperty(prefix+"asym_distance"));
  			if (properties.getProperty(prefix+"dct_window")!=null) this.dct_window=Integer.parseInt(properties.getProperty(prefix+"dct_window"));
  			if (properties.getProperty(prefix+"compactness")!=null) this.compactness=Double.parseDouble(properties.getProperty(prefix+"compactness"));
  			if (properties.getProperty(prefix+"sym_compactness")!=null) this.sym_compactness=Double.parseDouble(properties.getProperty(prefix+"sym_compactness"));
  			if (properties.getProperty(prefix+"dc_weight")!=null) this.dc_weight=Double.parseDouble(properties.getProperty(prefix+"dc_weight"));
  			if (properties.getProperty(prefix+"fact_precision")!=null) this.fact_precision=Double.parseDouble(properties.getProperty(prefix+"fact_precision"));
  			if (properties.getProperty(prefix+"asym_tax_free")!=null) this.asym_tax_free=Integer.parseInt(properties.getProperty(prefix+"asym_tax_free"));
  			if (properties.getProperty(prefix+"seed_size")!=null) this.seed_size=Integer.parseInt(properties.getProperty(prefix+"seed_size"));
  			if (properties.getProperty(prefix+"asym_random")!=null) this.asym_random=Double.parseDouble(properties.getProperty(prefix+"asym_random"));
  			if (properties.getProperty(prefix+"LMA_steps")!=null) this.LMA_steps=Integer.parseInt(properties.getProperty(prefix+"LMA_steps"));
  			if (properties.getProperty(prefix+"dbg_x")!=null) this.dbg_x=Double.parseDouble(properties.getProperty(prefix+"dbg_x"));
  			if (properties.getProperty(prefix+"dbg_y")!=null) this.dbg_y=Double.parseDouble(properties.getProperty(prefix+"dbg_y"));
  			if (properties.getProperty(prefix+"dbg_x1")!=null) this.dbg_x1=Double.parseDouble(properties.getProperty(prefix+"dbg_x1"));
  			if (properties.getProperty(prefix+"dbg_y1")!=null) this.dbg_y1=Double.parseDouble(properties.getProperty(prefix+"dbg_y1"));
  			if (properties.getProperty(prefix+"dbg_sigma")!=null) this.dbg_sigma=Double.parseDouble(properties.getProperty(prefix+"dbg_sigma"));
  			if (properties.getProperty(prefix+"dbg_src_size")!=null) this.dbg_src_size=Double.parseDouble(properties.getProperty(prefix+"dbg_src_size"));
  			if (properties.getProperty(prefix+"dbg_scale")!=null) this.dbg_scale=Double.parseDouble(properties.getProperty(prefix+"dbg_scale"));
  			if (properties.getProperty(prefix+"dbg_fold_scale")!=null) this.dbg_fold_scale=Double.parseDouble(properties.getProperty(prefix+"dbg_fold_scale"));
  			if (properties.getProperty(prefix+"dbg_mask")!=null) this.dbg_mask=properties.getProperty(prefix+"dbg_mask");
  			if (properties.getProperty(prefix+"dbg_mode")!=null) this.dbg_mode=Integer.parseInt(properties.getProperty(prefix+"dbg_mode"));
  			if (properties.getProperty(prefix+"centerWindowToTarget")!=null) this.centerWindowToTarget=Boolean.parseBoolean(properties.getProperty(prefix+"centerWindowToTarget"));
  			if (properties.getProperty(prefix+"color_channel")!=null) this.color_channel=Integer.parseInt(properties.getProperty(prefix+"color_channel"));
  			if (properties.getProperty(prefix+"decimation")!=null) this.decimation=Integer.parseInt(properties.getProperty(prefix+"decimation"));
  			if (properties.getProperty(prefix+"decimateSigma")!=null) this.decimateSigma=Double.parseDouble(properties.getProperty(prefix+"decimateSigma"));
  			if (properties.getProperty(prefix+"tileX")!=null) this.tileX=Integer.parseInt(properties.getProperty(prefix+"tileX"));
  			if (properties.getProperty(prefix+"tileY")!=null) this.tileY=Integer.parseInt(properties.getProperty(prefix+"tileY"));
  			if (properties.getProperty(prefix+"dbg_window_mode")!=null) this.dbg_window_mode=Integer.parseInt(properties.getProperty(prefix+"dbg_window_mode"));
  			if (properties.getProperty(prefix+"kernel_chn")!=null) this.kernel_chn=Integer.parseInt(properties.getProperty(prefix+"kernel_chn"));
  			if (properties.getProperty(prefix+"normalize")!=null) this.normalize=Boolean.parseBoolean(properties.getProperty(prefix+"normalize"));
  			if (properties.getProperty(prefix+"normalize_sym")!=null) this.normalize_sym=Boolean.parseBoolean(properties.getProperty(prefix+"normalize_sym"));
  			if (properties.getProperty(prefix+"antiwindow")!=null) this.antiwindow=Boolean.parseBoolean(properties.getProperty(prefix+"antiwindow"));
  			if (properties.getProperty(prefix+"skip_sym")!=null) this.skip_sym=Boolean.parseBoolean(properties.getProperty(prefix+"skip_sym"));
  			if (properties.getProperty(prefix+"convolve_direct")!=null) this.convolve_direct=Boolean.parseBoolean(properties.getProperty(prefix+"convolve_direct"));
  			if (properties.getProperty(prefix+"novignetting_r")!=null) this.novignetting_r=Double.parseDouble(properties.getProperty(prefix+"novignetting_r"));
  			if (properties.getProperty(prefix+"novignetting_g")!=null) this.novignetting_g=Double.parseDouble(properties.getProperty(prefix+"novignetting_g"));
  			if (properties.getProperty(prefix+"novignetting_b")!=null) this.novignetting_b=Double.parseDouble(properties.getProperty(prefix+"novignetting_b"));
  			if (properties.getProperty(prefix+"scale_r")!=null)        this.scale_r=Double.parseDouble(properties.getProperty(prefix+"scale_r"));
  			if (properties.getProperty(prefix+"scale_g")!=null)        this.scale_g=Double.parseDouble(properties.getProperty(prefix+"scale_g"));
  			if (properties.getProperty(prefix+"scale_b")!=null)        this.scale_b=Double.parseDouble(properties.getProperty(prefix+"scale_b"));
  			if (properties.getProperty(prefix+"vignetting_max")!=null) this.vignetting_max=Double.parseDouble(properties.getProperty(prefix+"vignetting_max"));
  			if (properties.getProperty(prefix+"vignetting_range")!=null) this.vignetting_range=Double.parseDouble(properties.getProperty(prefix+"vignetting_range"));
  			if (properties.getProperty(prefix+"post_debayer")!=null)   this.post_debayer=Boolean.parseBoolean(properties.getProperty(prefix+"post_debayer"));
  			if (properties.getProperty(prefix+"color_DCT")!=null)      this.color_DCT=Boolean.parseBoolean(properties.getProperty(prefix+"color_DCT"));
  			if (properties.getProperty(prefix+"sigma_rb")!=null)       this.sigma_rb=Double.parseDouble(properties.getProperty(prefix+"sigma_rb"));
  			if (properties.getProperty(prefix+"sigma_y")!=null)        this.sigma_y=Double.parseDouble(properties.getProperty(prefix+"sigma_y"));
  			if (properties.getProperty(prefix+"sigma_color")!=null)    this.sigma_color=Double.parseDouble(properties.getProperty(prefix+"sigma_color"));
  			if (properties.getProperty(prefix+"line_thershold")!=null) this.line_thershold=Double.parseDouble(properties.getProperty(prefix+"line_thershold"));
  			if (properties.getProperty(prefix+"nonlin")!=null)         this.nonlin=Boolean.parseBoolean(properties.getProperty(prefix+"nonlin"));
  			if (properties.getProperty(prefix+"nonlin_max_y")!=null)   this.nonlin_max_y=Double.parseDouble(properties.getProperty(prefix+"nonlin_max_y"));
  			if (properties.getProperty(prefix+"nonlin_max_c")!=null)   this.nonlin_max_c=Double.parseDouble(properties.getProperty(prefix+"nonlin_max_c"));
  			if (properties.getProperty(prefix+"nonlin_y")!=null)       this.nonlin_y=Double.parseDouble(properties.getProperty(prefix+"nonlin_y"));
  			if (properties.getProperty(prefix+"nonlin_c")!=null)       this.nonlin_c=Double.parseDouble(properties.getProperty(prefix+"nonlin_c"));
  			if (properties.getProperty(prefix+"nonlin_corn")!=null)    this.nonlin_corn=Double.parseDouble(properties.getProperty(prefix+"nonlin_corn"));
  			if (properties.getProperty(prefix+"denoise")!=null)        this.denoise=Boolean.parseBoolean(properties.getProperty(prefix+"denoise"));
  			if (properties.getProperty(prefix+"denoise_y")!=null)      this.denoise_y=Double.parseDouble(properties.getProperty(prefix+"denoise_y"));
  			if (properties.getProperty(prefix+"denoise_c")!=null)      this.denoise_c=Double.parseDouble(properties.getProperty(prefix+"denoise_c"));
  			if (properties.getProperty(prefix+"denoise_y_corn")!=null) this.denoise_y_corn=Double.parseDouble(properties.getProperty(prefix+"denoise_y_corn"));
  			if (properties.getProperty(prefix+"denoise_c_corn")!=null) this.denoise_c_corn=Double.parseDouble(properties.getProperty(prefix+"denoise_c_corn"));
  			
  		}
  		public boolean showDialog() {
  			GenericDialog gd = new GenericDialog("Set DCT parameters");
  			gd.addNumericField("DCT size",                                                       this.dct_size,            0);
  			gd.addNumericField("Size of asymmetrical (non-DCT) kernel",                          this.asym_size,           0);
  			gd.addNumericField("Maximal number of non-zero pixels in direct convolution kernel", this.asym_pixels,         0);
  			gd.addNumericField("How far to try a new asym kernel pixel from existing ones",      this.asym_distance,       0);
  			gd.addNumericField("MDCT window type (0,1,2)",                                       this.dct_window,          0);
  			gd.addNumericField("LMA_steps",                                                      this.LMA_steps,           0);
  			gd.addNumericField("Compactness (punish off-center asym_kernel pixels (proportional to r^2)", this.compactness,2);
  			gd.addNumericField("Symmetrical kernel compactness (proportional to r^2)",           this.sym_compactness,     2);
  			gd.addNumericField("Relative importance of DC error to RMS",                         this.dc_weight,           2);
  			gd.addNumericField("Factorization target precision (stop if achieved)",              this.fact_precision,      4);
  			gd.addNumericField("Do not punish pixels in the square around center",               this.asym_tax_free,       0);
  			gd.addNumericField("Start asym_kernel with this number of pixels (0 - single, 4n+0 (X between cells), 4*n+1 - x around center cell",               this.seed_size,     0); //0..2
  			gd.addNumericField("Initialize asym_kernel with random numbers (amplitude)",         this.asym_random,         2);
  			gd.addNumericField("dbg_x",                                                          this.dbg_x,               2);
  			gd.addNumericField("dbg_y",                                                          this.dbg_y,               2);
  			gd.addNumericField("dbg_x1",                                                         this.dbg_x1,              2);
  			gd.addNumericField("dbg_y1",                                                         this.dbg_y1,              2);
  			gd.addNumericField("dbg_sigma",                                                      this.dbg_sigma,           3);
  			gd.addNumericField("== dct_size = 1:1, dct+1.0 - shrink dct(dct+1.0)",               this.dbg_src_size,        3);
  			gd.addNumericField("Should ==DCT_PARAMETERS.dct_size / DCT_PARAMETERS.dbg_src_size", this.dbg_scale,           3);
  			gd.addNumericField("Modifies window during MDCT->DCT-IV folding",                    this.dbg_fold_scale,      3);
  			gd.addStringField ("Debug mask (anything but * is false)",                           this.dbg_mask,          100);
  			gd.addNumericField("LMA implementation: 0 - old, 1 - new",                           this.dbg_mode,            0);
  			gd.addNumericField("Convolution window: 0 - none, [1 - square], 2 - sin, 3 - sin^2", this.dbg_window_mode,     0);
  			gd.addCheckbox    ("Center convolution window around target kernel center",          this.centerWindowToTarget);
  			gd.addNumericField("Color channel to extract kernel (<0 - use synthetic)",           this.color_channel,       0);
  			gd.addNumericField("Convolution kernel decimation (original is normally 2x)",        this.decimation,          0);
  			gd.addNumericField("Smooth convolution kernel before decimation",                    this.decimateSigma,       3);
  			gd.addNumericField("Tile X to extract (0..163)",                                     this.tileX,               0);
  			gd.addNumericField("Tile Y to extract (0..122)",                                     this.tileY,               0);
  			gd.addNumericField("Calibration channel to use for aberration ( <0 - no correction)",this.kernel_chn,          0);
  			gd.addCheckbox    ("Normalize both sym and asym kernels ",                           this.normalize);
  			gd.addCheckbox    ("Normalize sym kernels separately",                               this.normalize_sym);
  			gd.addCheckbox    ("Divide symmetrical kernel by a window function",                 this.antiwindow);
  			gd.addCheckbox    ("Do not apply symmetrical (DCT) correction ",                     this.skip_sym);
  			gd.addCheckbox    ("Convolve directly with symmetrical kernel (debug feature) ",     this.convolve_direct);
  			gd.addNumericField("Reg gain in the center of sensor calibration R (instead of vignetting)",this.novignetting_r,   4);
  			gd.addNumericField("Green gain in the center of sensor calibration G (instead of vignetting)",this.novignetting_g, 4);
  			gd.addNumericField("Blue gain in the center of sensor calibration B (instead of vignetting)",this.novignetting_b,  4);
  			gd.addNumericField("Extra red correction to compensate for light temperature",               this.scale_r,  4);
  			gd.addNumericField("Extra green correction to compensate for light temperature",             this.scale_g,  4);
  			gd.addNumericField("Extra blue correction to compensate for light temperature",              this.scale_b,  4);
  			gd.addNumericField("Value (max) in vignetting data to correspond to 1x in the kernel",    this.vignetting_max,      3);
  			gd.addNumericField("Do not try to correct vignetting smaller than this fraction of max",  this.vignetting_range,  3);
  			gd.addCheckbox    ("Perform de-bayer after aberrations in pixel domain",              this.post_debayer             );
  			gd.addCheckbox    ("Use DCT-based color conversion (false - just LPF RGB with dbg_sigma)",this.color_DCT             );
  			gd.addNumericField("Gaussian sigma to apply to R and B (in addition to G), pix",      this.sigma_rb,            3);
  			gd.addNumericField("Gaussian sigma to apply to Y in the MDCT domain, pix",            this.sigma_y,             3);
  			gd.addNumericField("Gaussian sigma to apply to Pr and Pb in the MDCT domain, pix",    this.sigma_color,         3);
  			gd.addNumericField("Threshold for line detection (not yet used)",                     this.line_thershold,      3);

  			gd.addCheckbox    ("Use non-linear line emphasis and denoise",                        this.nonlin             );
  			gd.addNumericField("Maximal amount of non-linear emphasis for linear edges for Y component",  this.nonlin_max_y,3);
  			gd.addNumericField("Maximal amount of non-linear emphasis for linear edges for color diffs.", this.nonlin_max_c,3);
  			gd.addNumericField("Sensitivity of non-linear emphasis for linear edges for Y component",  this.nonlin_y,       3);
  			gd.addNumericField("Sensitivity of non-linear emphasis for linear edges for color diffs.", this.nonlin_c,       3);
  			gd.addNumericField("Corretion for diagonal/corner emphasis elements",                 this.nonlin_corn,         3);
  			gd.addCheckbox    ("Suppress noise during nonlinear processing",                      this.denoise             );
  			gd.addNumericField("Maximal total smoothing of the Y post-kernel (will compete with edge emphasis)",  this.denoise_y,3);
  			gd.addNumericField("Maximal total smoothing of the color differences post-kernel (will compete with edge emphasis)",  this.denoise_c,3);
  			gd.addNumericField("Weight of the 4 corner pixels during denoising y (straight - 1.0-denoise_y_corn)",  this.denoise_y_corn,3);
  			gd.addNumericField("Weight of the 4 corner pixels during denoising color ((straight - 1.0-denoise_c_corn))",  this.denoise_c_corn,3);
  			
  			WindowTools.addScrollBars(gd);
  			gd.showDialog();
  			
  			if (gd.wasCanceled()) return false;
  			this.dct_size=        (int) gd.getNextNumber();
  			this.asym_size=       (int) gd.getNextNumber();
  			this.asym_pixels=     (int) gd.getNextNumber();
  			this.asym_distance=   (int) gd.getNextNumber();
  			this.dct_window=      (int) gd.getNextNumber();
  			this.LMA_steps =      (int) gd.getNextNumber();
  			this.compactness =          gd.getNextNumber();
  			this.sym_compactness =      gd.getNextNumber();
  			this.dc_weight =            gd.getNextNumber();
  			this.fact_precision =       gd.getNextNumber();
  			this.asym_tax_free =  (int) gd.getNextNumber();
  			this.seed_size =      (int) gd.getNextNumber();
  			this.asym_random=           gd.getNextNumber();
  			this.dbg_x=                 gd.getNextNumber();
  			this.dbg_y=                 gd.getNextNumber();
  			this.dbg_x1=                gd.getNextNumber();
  			this.dbg_y1=                gd.getNextNumber();
  			this.dbg_sigma=             gd.getNextNumber();
  			this.dbg_src_size=          gd.getNextNumber();
  			this.dbg_scale=             gd.getNextNumber();
  			this.dbg_fold_scale=        gd.getNextNumber();
  			this.dbg_mask=              gd.getNextString();
  			this.dbg_mode=        (int) gd.getNextNumber();
  			this.dbg_window_mode= (int) gd.getNextNumber();
  			this.centerWindowToTarget=  gd.getNextBoolean();
  			this.color_channel=   (int) gd.getNextNumber();
  			this.decimation=      (int) gd.getNextNumber();
  			this.decimateSigma=         gd.getNextNumber();
  			this.tileX=           (int) gd.getNextNumber();
  			this.tileY=           (int) gd.getNextNumber();
  			this.kernel_chn=      (int) gd.getNextNumber();
  			this.normalize=             gd.getNextBoolean();
  			this.normalize_sym=         gd.getNextBoolean();
  			this.antiwindow=            gd.getNextBoolean();
  			this.skip_sym=              gd.getNextBoolean();
  			this.convolve_direct=       gd.getNextBoolean();

  			this.novignetting_r=        gd.getNextNumber();
  			this.novignetting_g=        gd.getNextNumber();
  			this.novignetting_b=        gd.getNextNumber();
  			this.scale_r=               gd.getNextNumber();
  			this.scale_g=               gd.getNextNumber();
  			this.scale_b=               gd.getNextNumber();
  			this.vignetting_max=        gd.getNextNumber();
  			this.vignetting_range=      gd.getNextNumber();

  			this.post_debayer=          gd.getNextBoolean();
  			this.color_DCT=             gd.getNextBoolean();
  			this.sigma_rb=              gd.getNextNumber();
  			this.sigma_y=               gd.getNextNumber();
  			this.sigma_color=           gd.getNextNumber();
  			this.line_thershold=        gd.getNextNumber();
  			
  			this.nonlin=             gd.getNextBoolean();
  			this.nonlin_max_y=          gd.getNextNumber();
  			this.nonlin_max_c=          gd.getNextNumber();
  			this.nonlin_y=              gd.getNextNumber();
  			this.nonlin_c=              gd.getNextNumber();
  			this.nonlin_corn=           gd.getNextNumber();
  			
  			this.denoise=             gd.getNextBoolean();
  			this.denoise_y=          gd.getNextNumber();
  			this.denoise_c=          gd.getNextNumber();
  			this.denoise_y_corn=          gd.getNextNumber();
  			this.denoise_c_corn=          gd.getNextNumber();
  			
  			//  	    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
  			return true;
  		}  
  	}

    
    /* ======================================================================== */
    public static class DebayerParameters {
  		public int size;
  		public double polarStep;
  		public double debayerThreshold;
  		public double debayerRelativeWidthGreen;
  		public double debayerRelativeWidthRedblue;
  		public double debayerRelativeWidthRedblueMain;
  		public double debayerRelativeWidthRedblueClones;
  		public double debayerGamma;
  		public double debayerBonus;
  		public double mainToAlias;
  		public double debayerMaskBlur;
  		public boolean debayerUseScissors;
  		public boolean debug;
  		public int xDebug;
  		public int yDebug;
  		public boolean debayerStacks;
  		public DebayerParameters(int size, double polarStep,
  				double debayerThreshold, double debayerRelativeWidthGreen,
  				double debayerRelativeWidthRedblue,
  				double debayerRelativeWidthRedblueMain,
  				double debayerRelativeWidthRedblueClones, double debayerGamma,
  				double debayerBonus, double mainToAlias, double debayerMaskBlur,
  				boolean debayerUseScissors, 
  				boolean debug, int xDebug, int yDebug,
  				boolean debayerStacks) {
  			this.size = size;
  			this.polarStep = polarStep;
  			this.debayerThreshold = debayerThreshold;
  			this.debayerRelativeWidthGreen = debayerRelativeWidthGreen;
  			this.debayerRelativeWidthRedblue = debayerRelativeWidthRedblue;
  			this.debayerRelativeWidthRedblueMain = debayerRelativeWidthRedblueMain;
  			this.debayerRelativeWidthRedblueClones = debayerRelativeWidthRedblueClones;
  			this.debayerGamma = debayerGamma;
  			this.debayerBonus = debayerBonus;
  			this.mainToAlias = mainToAlias;
  			this.debayerMaskBlur = debayerMaskBlur;
  			this.debayerUseScissors = debayerUseScissors;
  			this.debug = debug;
  			this.xDebug = xDebug;
  			this.yDebug = yDebug;
  			this.debayerStacks = debayerStacks;
  		}
  		public void setProperties(String prefix,Properties properties){
//  			properties.setProperty(prefix+"oversample",this.oversample+"");
  			properties.setProperty(prefix+"size",this.size+"");
  			properties.setProperty(prefix+"polarStep",this.polarStep+"");
  			properties.setProperty(prefix+"debayerThreshold",this.debayerThreshold+"");
  			properties.setProperty(prefix+"debayerRelativeWidthGreen",this.debayerRelativeWidthGreen+"");
  			properties.setProperty(prefix+"debayerRelativeWidthRedblue",this.debayerRelativeWidthRedblue+"");
  			properties.setProperty(prefix+"debayerRelativeWidthRedblueMain",this.debayerRelativeWidthRedblueMain+"");
  			properties.setProperty(prefix+"debayerRelativeWidthRedblueClones",this.debayerRelativeWidthRedblueClones+"");
  			properties.setProperty(prefix+"debayerGamma",this.debayerGamma+"");
  			properties.setProperty(prefix+"debayerBonus",this.debayerBonus+"");
  			properties.setProperty(prefix+"mainToAlias",this.mainToAlias+"");
  			properties.setProperty(prefix+"debayerMaskBlur",this.debayerMaskBlur+"");
  			properties.setProperty(prefix+"debayerUseScissors",this.debayerUseScissors+"");
  			properties.setProperty(prefix+"debug",this.debug+"");
  			properties.setProperty(prefix+"xDebug",this.xDebug+"");
  			properties.setProperty(prefix+"yDebug",this.yDebug+"");
  			properties.setProperty(prefix+"debayerStacks",this.debayerStacks+"");
  		}
  		public void getProperties(String prefix,Properties properties){
  			this.size=                             Integer.parseInt(properties.getProperty(prefix+"size"));
  			this.polarStep=                        Double.parseDouble(properties.getProperty(prefix+"polarStep"));
  			this.debayerThreshold=                 Double.parseDouble(properties.getProperty(prefix+"debayerThreshold"));
  			this.debayerRelativeWidthGreen=        Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthGreen"));
  			this.debayerRelativeWidthRedblue=      Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblue"));
  			this.debayerRelativeWidthRedblueMain=  Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblueMain"));
  			this.debayerRelativeWidthRedblueClones=Double.parseDouble(properties.getProperty(prefix+"debayerRelativeWidthRedblueClones"));
  			this.debayerGamma=                     Double.parseDouble(properties.getProperty(prefix+"debayerGamma"));
  			this.debayerBonus=                     Double.parseDouble(properties.getProperty(prefix+"debayerBonus"));
  			this.mainToAlias=                      Double.parseDouble(properties.getProperty(prefix+"mainToAlias"));
  			this.debayerMaskBlur=                  Double.parseDouble(properties.getProperty(prefix+"debayerMaskBlur"));
  			this.debayerUseScissors=               Boolean.parseBoolean(properties.getProperty(prefix+"debayerUseScissors"));
  			this.debug=                            Boolean.parseBoolean(properties.getProperty(prefix+"debug"));
  			this.xDebug=                           Integer.parseInt(properties.getProperty(prefix+"xDebug"));
  			this.yDebug=                           Integer.parseInt(properties.getProperty(prefix+"yDebug"));
  			this.debayerStacks=                    Boolean.parseBoolean(properties.getProperty(prefix+"debayerStacks"));
  		}
  	}

    public static class EquirectangularParameters {
    	
		public double longitudeLeft=    -180.0;
		public double longitudeRight=    180.0;
		public double latitudeTop=        90.0;
		public double latitudeBottom=    -90.0;
		public int pixelsHorizontal=   14268;
		public int imageWidth=          2592;
		public int imageHeight=         1936;
		public double resolutionScale=     1.0;
		public double x0=                  0.0;
		public double y0=                  0.0;
		public int longitudeWidth=      3000; //pix
		public boolean clearFullMap=      true;
		public boolean clearAllMaps=      true;
		
		public boolean needRebuild=      false;
// common plane parameters (dual camera, triclope camera)
		public boolean generateCommonPlane =  false;
		public double projectionElevation=     0.0;
		public double projectionYaw=           0.0;
		public double projectionRoll=          0.0;
		
		public boolean matchPixelSize=    true; // disregard next value, calculate projectionPixelSize from the equirectangular map
		public double projectionPixelSize=0.00044036902;
		public int    projectionWidth=   2920;
		public int     projectionHeight= 2220;
		public double projectionCenterX= 0.5*this.projectionWidth;
		public double projectionCenterY= 0.5*this.projectionHeight;
		public double nominalHorizontalDisparity=60.0; // nominal distance between horizontal cameras, mm
		public boolean [] channelSelection=null;
    	

    	public EquirectangularParameters(){
    	}
    	public boolean isNeedRebuildSet(){
    		boolean result=this.needRebuild;
    		this.needRebuild=false;
    		return result;
    	}
    	
    	public boolean [] getChannelsToProcess(){ return this.channelSelection;}

    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"longitudeLeft",this.longitudeLeft+"");
    		properties.setProperty(prefix+"longitudeRight",this.longitudeRight+"");
    		properties.setProperty(prefix+"latitudeTop",this.latitudeTop+"");
    		properties.setProperty(prefix+"latitudeBottom",this.latitudeBottom+"");
    		
    		properties.setProperty(prefix+"pixelsHorizontal",this.pixelsHorizontal+"");
    		properties.setProperty(prefix+"imageWidth",this.imageWidth+"");
    		properties.setProperty(prefix+"imageHeight",this.imageHeight+"");
    		
    		properties.setProperty(prefix+"resolutionScale",this.resolutionScale+"");
    		properties.setProperty(prefix+"x0",this.x0+"");
    		properties.setProperty(prefix+"y0",this.y0+"");
    		
    		properties.setProperty(prefix+"longitudeWidth",this.longitudeWidth+"");
    		
    		properties.setProperty(prefix+"clearFullMap",this.clearFullMap+"");
    		properties.setProperty(prefix+"clearAllMaps",this.clearAllMaps+"");
    		properties.setProperty(prefix+"generateCommonPlane",this.generateCommonPlane+"");
    		properties.setProperty(prefix+"projectionElevation",this.projectionElevation+"");
    		properties.setProperty(prefix+"projectionYaw",this.projectionYaw+"");
    		properties.setProperty(prefix+"projectionRoll",this.projectionRoll+"");
    		properties.setProperty(prefix+"matchPixelSize",this.matchPixelSize+"");
    		properties.setProperty(prefix+"projectionPixelSize",this.projectionPixelSize+"");
    		properties.setProperty(prefix+"projectionWidth",this.projectionWidth+"");
    		properties.setProperty(prefix+"projectionHeight",this.projectionHeight+"");
    		properties.setProperty(prefix+"projectionCenterX",this.projectionCenterX+"");
    		properties.setProperty(prefix+"projectionCenterY",this.projectionCenterY+"");
    		properties.setProperty(prefix+"nominalHorizontalDisparity",this.nominalHorizontalDisparity+"");
    		if (this.channelSelection!=null){
    			properties.setProperty(prefix+"numberProjectedChannels",this.channelSelection.length+"");
    			for (int i=0;i<this.channelSelection.length;i++){
    				properties.setProperty(prefix+"projectedChannel"+i,this.channelSelection[i]+"");
    			}
    		}
    		
    	}
    	public void getProperties(String prefix,Properties properties){

    		if (properties.getProperty(prefix+"longitudeLeft")!=null) this.longitudeLeft=      Double.parseDouble(properties.getProperty(prefix+"longitudeLeft"));
    		if (properties.getProperty(prefix+"longitudeRight")!=null) this.longitudeRight=      Double.parseDouble(properties.getProperty(prefix+"longitudeRight"));
    		if (properties.getProperty(prefix+"latitudeTop")!=null)this.latitudeTop=       Double.parseDouble(properties.getProperty(prefix+"latitudeTop"));
    		if (properties.getProperty(prefix+"latitudeBottom")!=null)this.latitudeBottom=       Double.parseDouble(properties.getProperty(prefix+"latitudeBottom"));

    		if (properties.getProperty(prefix+"pixelsHorizontal")!=null)this.pixelsHorizontal=          Integer.parseInt(properties.getProperty(prefix+"pixelsHorizontal"));
    		if (properties.getProperty(prefix+"imageWidth")!=null)this.imageWidth=          Integer.parseInt(properties.getProperty(prefix+"imageWidth"));
    		if (properties.getProperty(prefix+"imageHeight")!=null)this.imageHeight=          Integer.parseInt(properties.getProperty(prefix+"imageHeight"));

    		if (properties.getProperty(prefix+"resolutionScale")!=null)this.resolutionScale=       Double.parseDouble(properties.getProperty(prefix+"resolutionScale"));
    		if (properties.getProperty(prefix+"x0")!=null)this.x0=       Double.parseDouble(properties.getProperty(prefix+"x0"));
    		if (properties.getProperty(prefix+"y0")!=null)this.y0=       Double.parseDouble(properties.getProperty(prefix+"y0"));
  			
    		if (properties.getProperty(prefix+"longitudeWidth")!=null)this.longitudeWidth=          Integer.parseInt(properties.getProperty(prefix+"longitudeWidth"));
    		
    		if (properties.getProperty(prefix+"clearFullMap")!=null)this.clearFullMap=   Boolean.parseBoolean(properties.getProperty(prefix+"clearFullMap"));
    		if (properties.getProperty(prefix+"clearAllMaps")!=null)this.clearAllMaps=   Boolean.parseBoolean(properties.getProperty(prefix+"clearAllMaps"));
      		
    		if (properties.getProperty(prefix+"generateCommonPlane")!=null)this.generateCommonPlane=       Boolean.parseBoolean(properties.getProperty(prefix+"generateCommonPlane"));
    		if (properties.getProperty(prefix+"projectionElevation")!=null)this.projectionElevation=       Double.parseDouble(properties.getProperty(prefix+"projectionElevation"));
    		if (properties.getProperty(prefix+"projectionYaw")!=null)this.projectionYaw=       Double.parseDouble(properties.getProperty(prefix+"projectionYaw"));
    		if (properties.getProperty(prefix+"projectionRoll")!=null)this.projectionRoll=       Double.parseDouble(properties.getProperty(prefix+"projectionRoll"));
    		if (properties.getProperty(prefix+"matchPixelSize")!=null)this.matchPixelSize=   Boolean.parseBoolean(properties.getProperty(prefix+"matchPixelSize"));
    		if (properties.getProperty(prefix+"projectionPixelSize")!=null)this.projectionPixelSize=       Double.parseDouble(properties.getProperty(prefix+"projectionPixelSize"));

    		if (properties.getProperty(prefix+"projectionWidth")!=null)this.projectionWidth=          Integer.parseInt(properties.getProperty(prefix+"projectionWidth"));
    		if (properties.getProperty(prefix+"projectionHeight")!=null)this.projectionHeight=          Integer.parseInt(properties.getProperty(prefix+"projectionHeight"));
    		
    		if (properties.getProperty(prefix+"projectionCenterX")!=null)this.projectionCenterX=       Double.parseDouble(properties.getProperty(prefix+"projectionCenterX"));
    		if (properties.getProperty(prefix+"projectionCenterY")!=null)this.projectionCenterY=       Double.parseDouble(properties.getProperty(prefix+"projectionCenterY"));
    		if (properties.getProperty(prefix+"nominalHorizontalDisparity")!=null)this.nominalHorizontalDisparity=       Double.parseDouble(properties.getProperty(prefix+"nominalHorizontalDisparity"));
    		if (properties.getProperty(prefix+"numberProjectedChannels")!=null){
    			int numberProjectedChannels=Integer.parseInt(properties.getProperty(prefix+"numberProjectedChannels"));
    			this.channelSelection=new boolean[numberProjectedChannels];
    			for (int i=0;i<this.channelSelection.length;i++){
    	    		if (properties.getProperty(prefix+"projectedChannel"+i)!=null)
    	    			this.channelSelection[i]=   Boolean.parseBoolean(properties.getProperty(prefix+"projectedChannel"+i));
    			}
    		}

    	}
    	
    	
    	public boolean showDialog() {
    		needRebuild=false;
			GenericDialog gd=new GenericDialog("Select parameters for equirectangular->sensor pixel mapping");
			gd.addMessage("Equirectangular area");
			gd.addNumericField("Longitude left", this.longitudeLeft, 1,6,"degrees" );
			gd.addNumericField("Longitude right", this.longitudeRight, 1,6,"degrees" );
			gd.addNumericField("Latitude top", this.latitudeTop, 1,6,"degrees" );
			gd.addNumericField("Latitude bottom", this.latitudeBottom, 1,6,"degrees" );
			gd.addNumericField("Pixels horizontal ", this.pixelsHorizontal,0,5,"image pix");
			gd.addMessage("Source image parameters");
			gd.addNumericField("Input image width", this.imageWidth,0,4,"image pix");
			gd.addNumericField("Input image height", this.imageHeight,0,4,"image pix");
			gd.addNumericField("Input image resolution scale (2.0 - twice resolution, 0.5 - half)", this.resolutionScale, 4,6,"x" );
			gd.addNumericField("Input image left margin", this.x0, 1,6,"sensor pix" );
			gd.addNumericField("Input image top margin", this.y0, 1,6,"sensor pix" );
			gd.addMessage("Reduction of memory usage");
			gd.addNumericField("Crop files horizontally to ", this.longitudeWidth,0,4,"longitude pix");
			gd.addCheckbox    ("Clear full map",    this.clearFullMap);
			gd.addCheckbox    ("Clear all data",    this.clearFullMap);
			gd.addMessage("Parameters for the common projection plane (binocular/trinocular cameras)");
			gd.addCheckbox    ("Generate common projection plane",    this.generateCommonPlane);
			gd.addNumericField("View axis elevation (orthogonal to projection plane)",this.projectionElevation,2,6,"degrees");
			gd.addNumericField("View axis heading   (orthogonal to projection plane)",this.projectionYaw,2,6,"degrees");
			gd.addNumericField("View plane rotation (roll) around the view axis",     this.projectionRoll,2,6,"degrees");
			gd.addCheckbox    ("Match projection pixel size to that of the equirectangular map",    this.matchPixelSize);
			gd.addNumericField("Projection pixel size (relative)     ",1000*this.projectionPixelSize,4,8,"x1/1000");
			gd.addNumericField("Projection plane width", this.projectionWidth,0,5,"pix");
			gd.addNumericField("Projection plane height",this.projectionHeight,0,5,"pix");
			gd.addNumericField("Projection plane Center X (point orthogonal to the view axis), right",   this.projectionCenterX,2,8,"pix");
			gd.addNumericField("Projection plane Center Y (point orthogonal to the view axis), down",    this.projectionCenterY,2,8,"pix");
			gd.addNumericField("Nominal distance between the 2 horizontal cameras",    this.nominalHorizontalDisparity,2,8,"mm");
			gd.enableYesNoCancel("OK", "Rebuild map files");
    		WindowTools.addScrollBars(gd);
			gd.showDialog();
			if (gd.wasCanceled()) return false;
			this.longitudeLeft=          gd.getNextNumber();
			this.longitudeRight=         gd.getNextNumber();
			this.latitudeTop=            gd.getNextNumber();
			this.latitudeBottom=         gd.getNextNumber();
			this.pixelsHorizontal= (int) gd.getNextNumber();
			this.imageWidth=      (int)  gd.getNextNumber();
			this.imageHeight=     (int)  gd.getNextNumber();
			this.resolutionScale=        gd.getNextNumber();
			this.x0=                     gd.getNextNumber();
			this.y0=                     gd.getNextNumber();
			this.longitudeWidth=   (int) gd.getNextNumber();
			this.clearFullMap=           gd.getNextBoolean();
			this.clearAllMaps=           gd.getNextBoolean();
    		this.generateCommonPlane =   gd.getNextBoolean();
			this.projectionElevation=    gd.getNextNumber();
			this.projectionYaw=          gd.getNextNumber();
			this.projectionRoll=         gd.getNextNumber();
			this.matchPixelSize =        gd.getNextBoolean();
			this.projectionPixelSize=    0.001*gd.getNextNumber();
			this.projectionWidth=  (int) gd.getNextNumber();
			this.projectionHeight= (int) gd.getNextNumber();
			this.projectionCenterX=      gd.getNextNumber();
			this.projectionCenterY=      gd.getNextNumber();
			this.nominalHorizontalDisparity=      gd.getNextNumber();
			
			if (!gd.wasOKed()) needRebuild=true;
    		return true;
    	}
    	public boolean selectChannelsToProcess(String title, int numChannels) {
    		if (numChannels<=0){
    			this.channelSelection=null;
    			return true;
    		}
    		boolean [] newSelecttion=new boolean [numChannels];
    		boolean lastChoice=true; // selected
    		for (int i=0;i<numChannels;i++){
    			if ((this.channelSelection!=null) && (i<this.channelSelection.length)){
    				newSelecttion[i]=this.channelSelection[i];
    				lastChoice=newSelecttion[i];
    			} else newSelecttion[i]=lastChoice;
    		}
    		while (true) {
    			GenericDialog gd = new GenericDialog(title);
    			for (int i=0;i<numChannels;i++) gd.addCheckbox("channel "+i, newSelecttion[i]);
    			gd.enableYesNoCancel("OK", "All like channel 0");
    			WindowTools.addScrollBars(gd);
    			gd.showDialog();
    			if (gd.wasCanceled()) return false; // but do not modify this.channelSelection
    			for (int i=0;i<numChannels;i++) newSelecttion[i]=gd.getNextBoolean();
    			if (gd.wasOKed()){
    				this.channelSelection=newSelecttion;
    				return true;
    			} else {
    				for (int i=1;i<numChannels;i++) newSelecttion[i]=newSelecttion[0];
    			}
    		}
    	}


    }

    
    
  /* ======================================================================== */

}
