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
  		public boolean crop =                  true;  // crop to the sennor size 
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
  		public boolean    show_nonoverlap =   true;  // show result RGBA before overlap combined (first channels, then RGBA combined?)
  		public boolean    show_overlap =      true;  // show result RGBA (first channels, then RGBA combined?)
  		public boolean    show_rgba_color =   true;  // show combined color image
  		public boolean    show_map =          true;  // show disparity maps
  		
  		public double     disp_scan_start =   0.0;   // disparity scan start value
  		public double     disp_scan_step =    1.0;   // disparity scan step
  		public int        disp_scan_count =   10;    // disparity scan number of measurements
  		
  		public double     fine_corr_x_0 =     0.0;   // additionally shift image in port 0 in x direction
  		public double     fine_corr_y_0 =     0.0;   // additionally shift image in port 0 in y direction
  		public double     fine_corr_x_1 =     0.0;   // additionally shift image in port 1 in x direction
  		public double     fine_corr_y_1 =     0.0;   // additionally shift image in port 1 in y direction
  		public double     fine_corr_x_2 =     0.0;   // additionally shift image in port 2 in x direction
  		public double     fine_corr_y_2 =     0.0;   // additionally shift image in port 2 in y direction
  		public double     fine_corr_x_3 =     0.0;   // additionally shift image in port 3 in x direction
  		public double     fine_corr_y_3 =     0.0;   // additionally shift image in port 3 in y direction
  		
  		public double     fcorr_min_stength = 0.005; // minimal correlation strength to apply fine correction
  		public double     fcorr_disp_diff =   3.0;   // consider only tiles with absolute residual disparity lower than
  		public boolean    fcorr_quadratic =   true;  // Use quadratic polynomial for fine correction (false - only linear)
  		public boolean    fcorr_ignore =      false; // Ignore currently calculated fine correction
  		
  		public double     corr_magic_scale =  0.85;  // reported correlation offset vs. actual one (not yet understood)
  		
  		// 3d reconstruction
  		public boolean    show_textures    = true;  // show generated textures
  		public boolean    debug_filters    = false;// show intermediate results of filtering
  		public double     min_smth         = 0.25;  // 0.25 minimal noise-normalized pixel difference in a channel to suspect something    
  		public double     sure_smth        = 2.0;   // reliable noise-normalized pixel difference in a channel to have something    
  		public double     bgnd_range       = 0.3;   // disparity range to be considered background
  		public double     other_range      = 2.0;   // disparity difference from center (provided) disparity to trust
  		
  		public double     bgnd_sure        = 0.18;  // minimal strength to be considered definitely background
  		public double     bgnd_maybe       = 0.1; // maximal strength to ignore as non-background
//  		public double     bgnd_2diff       = 0.005; // maximal strength to ignore as non-background
  		public int        min_clstr_seed   = 2;     // number of tiles in a cluster to seed (just background?)
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
  		
  		public int        max_clusters     = 300;    // Maximal number of clusters to generate for one run
  		public boolean    correct_distortions = false; // Correct lens geometric distortions in a model (will need backdrop to be corrected too)
  		public boolean    show_triangles =    true;  // Show generated triangles
  		public boolean    avg_cluster_disp =  false;  // Weight-average disparity for the whole cluster 
  		public double     maxDispTriangle   = 0.2;    // Maximal disparity difference in a triangle face to show
  		
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
  		public double     stStepFar         = 0.1;   // Disaprity histogram step for far objects 
  		public double     stStepNear        = 0.5;   // Disaprity histogram step for near objects
  		public double     stStepThreshold   = 1.0;   // Disaprity threshold to switch cfrom linear to logarithmic steps
  		public double     stMinDisparity    = 0.0;   // Minimal disparity (center of a bin)
  		public double     stMaxDisparity    = 15.0;  // Maximal disparity (center of a bin)
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
  		public double     unique_tolerance = 0.1; // Do not re-measure correlation if target disparity differs from some previous by this
  		
  		// Multi-pass growing disparity
  		public int        grow_sweep         = 8; // Try these number of tiles around known ones 
  		public double     grow_disp_max =   50.0; // Maximal disparity to try
  		public double     grow_disp_trust =  4.0; // Trust measured disparity within +/- this value 
  		public double     grow_disp_step =   6.0; // Increase disparity (from maximal tried) if nothing found in that tile // TODO: handle enclosed dips?  
  		public double     grow_min_diff =    0.5; // Grow more only if at least one channel has higher variance from others for the tile  
  		
  		public boolean    plPreferDisparity    =   false;// Always start with disparity-most axis (false - lowest eigenvalue)
  		public double     plDispNorm           =   3.0;  // Normalize disparities to the average if above (now only for eigenvalue comparison)
  		public int        plMinPoints          =     5;  // Minimal number of points for plane detection
  		public double     plTargetEigen        =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
  		public double     plFractOutliers      =   0.3;  // Maximal fraction of outliers to remove
  		public int        plMaxOutliers        =    20;  // Maximal number of outliers to remove
  		public double     plMinStrength        =   0.1;  // Minimal total strength of a plane 
  		public double     plMaxEigen           =   0.3;  // Maximal eigenvalue of a plane 
  		public boolean    plDbgMerge           =   true; // Combine 'other' plane with current
  		public double     plWorstWorsening     =   3.0;  // Worst case worsening after merge
  		public double     plWeakWorsening      =   1.0;  // Relax merge requirements for weaker planes
  		public boolean    plMutualOnly         =   true; // keep only mutual links, remove weakest if conflict
  		public boolean    plFillSquares        =   true; // Add diagonals to full squares
  		public boolean    plCutCorners         =   true; // Add ortho to 45-degree corners

  		public double     plPull               =  .3;   // Relative weight of original (measured) plane compared to average neighbor pull
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
  		
  		public boolean    replaceWeakOutlayers =   true; // false; 
  		
  		public boolean    dbg_migrate =            true; 
  		
  		// other debug images
  		public boolean    show_ortho_combine =     false; // Show 'ortho_combine' 
  		public boolean    show_refine_supertiles = false; // show 'refine_disparity_supertiles' 
  		public boolean    show_bgnd_nonbgnd =      false; // show 'bgnd_nonbgnd' 
  		public boolean    show_filter_scan =       false; // show 'FilterScan'
  		public boolean    show_combined =          false; // show 'combo_scan' (combined multiple scans)
  		public boolean    show_unique =            false; // show 'unique_scan' (removed already measured tiles with the same disparity)
  		public boolean    show_init_refine =       false; // show debug images during initial refinement 
  		public boolean    show_expand =            false; // show debug images during disparity expansion
  		
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
			properties.setProperty(prefix+"show_nonoverlap",  this.show_nonoverlap+"");
			properties.setProperty(prefix+"show_overlap",     this.show_overlap+"");
			properties.setProperty(prefix+"show_rgba_color",  this.show_rgba_color+"");
			properties.setProperty(prefix+"show_map",         this.show_map+"");
  			properties.setProperty(prefix+"disp_scan_start",  this.disp_scan_start +"");
  			properties.setProperty(prefix+"disp_scan_step",   this.disp_scan_step +"");
  			properties.setProperty(prefix+"disp_scan_count",  this.disp_scan_count+"");

  			properties.setProperty(prefix+"fine_corr_x_0",    this.fine_corr_x_0 +"");
  			properties.setProperty(prefix+"fine_corr_y_0",    this.fine_corr_y_0 +"");
  			properties.setProperty(prefix+"fine_corr_x_1",    this.fine_corr_x_1 +"");
  			properties.setProperty(prefix+"fine_corr_y_1",    this.fine_corr_y_1 +"");
  			properties.setProperty(prefix+"fine_corr_x_2",    this.fine_corr_x_2 +"");
  			properties.setProperty(prefix+"fine_corr_y_2",    this.fine_corr_y_2 +"");
  			properties.setProperty(prefix+"fine_corr_x_3",    this.fine_corr_x_3 +"");
  			properties.setProperty(prefix+"fine_corr_y_3",    this.fine_corr_y_3 +"");

  			properties.setProperty(prefix+"fcorr_min_stength",this.fcorr_min_stength +"");
  			properties.setProperty(prefix+"fcorr_disp_diff",  this.fcorr_disp_diff +"");
			properties.setProperty(prefix+"fcorr_quadratic",  this.fcorr_quadratic+"");
			properties.setProperty(prefix+"fcorr_ignore",     this.fcorr_ignore+"");

			properties.setProperty(prefix+"corr_magic_scale", this.corr_magic_scale +"");
  			
			properties.setProperty(prefix+"show_textures",    this.show_textures+"");
			properties.setProperty(prefix+"debug_filters",    this.debug_filters+"");

			properties.setProperty(prefix+"min_smth",         this.min_smth +"");
			properties.setProperty(prefix+"sure_smth",        this.sure_smth +"");
			properties.setProperty(prefix+"bgnd_range",       this.bgnd_range +"");
			properties.setProperty(prefix+"other_range",      this.other_range +"");
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
			properties.setProperty(prefix+"correct_distortions",this.correct_distortions+"");
			properties.setProperty(prefix+"show_triangles",   this.show_triangles+"");
			properties.setProperty(prefix+"avg_cluster_disp", this.avg_cluster_disp+"");
			properties.setProperty(prefix+"maxDispTriangle",  this.maxDispTriangle +"");
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
			properties.setProperty(prefix+"stMaxDisparity",   this.stMaxDisparity +"");
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
			properties.setProperty(prefix+"unique_tolerance", this.unique_tolerance +"");
  			properties.setProperty(prefix+"grow_sweep",       this.grow_sweep+"");
			properties.setProperty(prefix+"grow_disp_max",    this.grow_disp_max +"");
			properties.setProperty(prefix+"grow_disp_trust",  this.grow_disp_trust +"");
			properties.setProperty(prefix+"grow_disp_step",   this.grow_disp_step +"");
			properties.setProperty(prefix+"grow_min_diff",    this.grow_min_diff +"");

			properties.setProperty(prefix+"plPreferDisparity",this.plPreferDisparity+"");
			properties.setProperty(prefix+"plDispNorm",       this.plDispNorm +"");
  			properties.setProperty(prefix+"plMinPoints",      this.plMinPoints+"");
			properties.setProperty(prefix+"plTargetEigen",    this.plTargetEigen +"");
			properties.setProperty(prefix+"plFractOutliers",  this.plFractOutliers +"");
  			properties.setProperty(prefix+"plMaxOutliers",    this.plMaxOutliers+"");
			properties.setProperty(prefix+"plMinStrength",    this.plMinStrength +"");
			properties.setProperty(prefix+"plMaxEigen",       this.plMaxEigen +"");
			properties.setProperty(prefix+"plDbgMerge",       this.plDbgMerge+"");
			properties.setProperty(prefix+"plWorstWorsening", this.plWorstWorsening +"");
			properties.setProperty(prefix+"plWeakWorsening",  this.plWeakWorsening +"");
			properties.setProperty(prefix+"plMutualOnly",     this.plMutualOnly+"");

			properties.setProperty(prefix+"plFillSquares",    this.plFillSquares+"");
			properties.setProperty(prefix+"plCutCorners",     this.plCutCorners+"");

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

			properties.setProperty(prefix+"dbg_migrate",            this.dbg_migrate+"");
  			
			properties.setProperty(prefix+"show_ortho_combine",     this.show_ortho_combine+"");
			properties.setProperty(prefix+"show_refine_supertiles", this.show_refine_supertiles+"");
			properties.setProperty(prefix+"show_bgnd_nonbgnd",      this.show_bgnd_nonbgnd+"");
			properties.setProperty(prefix+"show_filter_scan",       this.show_filter_scan+"");
			properties.setProperty(prefix+"show_combined",          this.show_combined+"");
			properties.setProperty(prefix+"show_unique",            this.show_unique+"");
			properties.setProperty(prefix+"show_init_refine",       this.show_init_refine+"");
			properties.setProperty(prefix+"show_expand",            this.show_expand+"");
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
  			if (properties.getProperty(prefix+"show_nonoverlap")!=null)this.show_nonoverlap=Boolean.parseBoolean(properties.getProperty(prefix+"show_nonoverlap"));
  			if (properties.getProperty(prefix+"show_overlap")!=null)   this.show_overlap=Boolean.parseBoolean(properties.getProperty(prefix+"show_overlap"));
  			if (properties.getProperty(prefix+"show_rgba_color")!=null)this.show_rgba_color=Boolean.parseBoolean(properties.getProperty(prefix+"show_rgba_color"));
  			if (properties.getProperty(prefix+"show_map")!=null)       this.show_map=Boolean.parseBoolean(properties.getProperty(prefix+"show_map"));
  			if (properties.getProperty(prefix+"disp_scan_start")!=null)this.disp_scan_start=Double.parseDouble(properties.getProperty(prefix+"disp_scan_start"));
  			if (properties.getProperty(prefix+"disp_scan_step")!=null) this.disp_scan_step=Double.parseDouble(properties.getProperty(prefix+"disp_scan_step"));
  			if (properties.getProperty(prefix+"disp_scan_count")!=null)this.disp_scan_count=Integer.parseInt(properties.getProperty(prefix+"disp_scan_count"));
  			
  			if (properties.getProperty(prefix+"fine_corr_x_0")!=null) this.fine_corr_x_0=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_0"));
  			if (properties.getProperty(prefix+"fine_corr_y_0")!=null) this.fine_corr_y_0=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_0"));
  			if (properties.getProperty(prefix+"fine_corr_x_1")!=null) this.fine_corr_x_1=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_1"));
  			if (properties.getProperty(prefix+"fine_corr_y_1")!=null) this.fine_corr_y_1=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_1"));
  			if (properties.getProperty(prefix+"fine_corr_x_2")!=null) this.fine_corr_x_2=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_2"));
  			if (properties.getProperty(prefix+"fine_corr_y_2")!=null) this.fine_corr_y_2=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_2"));
  			if (properties.getProperty(prefix+"fine_corr_x_3")!=null) this.fine_corr_x_3=Double.parseDouble(properties.getProperty(prefix+"fine_corr_x_3"));
  			if (properties.getProperty(prefix+"fine_corr_y_3")!=null) this.fine_corr_y_3=Double.parseDouble(properties.getProperty(prefix+"fine_corr_y_3"));
  			
  			if (properties.getProperty(prefix+"fcorr_min_stength")!=null) this.fcorr_min_stength=Double.parseDouble(properties.getProperty(prefix+"fcorr_min_stength"));
  			if (properties.getProperty(prefix+"fcorr_disp_diff")!=null)   this.fcorr_disp_diff=Double.parseDouble(properties.getProperty(prefix+"fcorr_disp_diff"));
  			if (properties.getProperty(prefix+"fcorr_quadratic")!=null)   this.fcorr_quadratic=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_quadratic"));
  			if (properties.getProperty(prefix+"fcorr_ignore")!=null)      this.fcorr_ignore=Boolean.parseBoolean(properties.getProperty(prefix+"fcorr_ignore"));
  			if (properties.getProperty(prefix+"corr_magic_scale")!=null)  this.corr_magic_scale=Double.parseDouble(properties.getProperty(prefix+"corr_magic_scale"));

  			if (properties.getProperty(prefix+"show_textures")!=null)      this.show_textures=Boolean.parseBoolean(properties.getProperty(prefix+"show_textures"));
  			if (properties.getProperty(prefix+"debug_filters")!=null)      this.debug_filters=Boolean.parseBoolean(properties.getProperty(prefix+"debug_filters"));

  			if (properties.getProperty(prefix+"min_smth")!=null)          this.min_smth=Double.parseDouble(properties.getProperty(prefix+"min_smth"));
  			if (properties.getProperty(prefix+"sure_smth")!=null)         this.sure_smth=Double.parseDouble(properties.getProperty(prefix+"sure_smth"));
  			if (properties.getProperty(prefix+"bgnd_range")!=null)        this.bgnd_range=Double.parseDouble(properties.getProperty(prefix+"bgnd_range"));
  			if (properties.getProperty(prefix+"other_range")!=null)       this.other_range=Double.parseDouble(properties.getProperty(prefix+"other_range"));
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
  			if (properties.getProperty(prefix+"correct_distortions")!=null) this.correct_distortions=Boolean.parseBoolean(properties.getProperty(prefix+"correct_distortions"));
  			if (properties.getProperty(prefix+"show_triangles")!=null)    this.show_triangles=Boolean.parseBoolean(properties.getProperty(prefix+"show_triangles"));
  			if (properties.getProperty(prefix+"avg_cluster_disp")!=null)  this.avg_cluster_disp=Boolean.parseBoolean(properties.getProperty(prefix+"avg_cluster_disp"));
  			if (properties.getProperty(prefix+"maxDispTriangle")!=null)   this.maxDispTriangle=Double.parseDouble(properties.getProperty(prefix+"maxDispTriangle"));
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
  			if (properties.getProperty(prefix+"stMaxDisparity")!=null)    this.stMaxDisparity=Double.parseDouble(properties.getProperty(prefix+"stMaxDisparity"));
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
  			if (properties.getProperty(prefix+"unique_tolerance")!=null)  this.unique_tolerance=Double.parseDouble(properties.getProperty(prefix+"unique_tolerance"));
  			if (properties.getProperty(prefix+"grow_sweep")!=null)        this.grow_sweep=Integer.parseInt(properties.getProperty(prefix+"grow_sweep"));
  			if (properties.getProperty(prefix+"grow_disp_max")!=null)     this.grow_disp_max=Double.parseDouble(properties.getProperty(prefix+"grow_disp_max"));
  			if (properties.getProperty(prefix+"grow_disp_trust")!=null)   this.grow_disp_trust=Double.parseDouble(properties.getProperty(prefix+"grow_disp_trust"));
  			if (properties.getProperty(prefix+"grow_disp_step")!=null)    this.grow_disp_step=Double.parseDouble(properties.getProperty(prefix+"grow_disp_step"));
  			if (properties.getProperty(prefix+"grow_min_diff")!=null)     this.grow_min_diff=Double.parseDouble(properties.getProperty(prefix+"grow_min_diff"));

  			if (properties.getProperty(prefix+"plPreferDisparity")!=null) this.plPreferDisparity=Boolean.parseBoolean(properties.getProperty(prefix+"plPreferDisparity"));
  			if (properties.getProperty(prefix+"plDispNorm")!=null)        this.plDispNorm=Double.parseDouble(properties.getProperty(prefix+"plDispNorm"));
  			if (properties.getProperty(prefix+"plMinPoints")!=null)       this.plMinPoints=Integer.parseInt(properties.getProperty(prefix+"plMinPoints"));
  			if (properties.getProperty(prefix+"plTargetEigen")!=null)     this.plTargetEigen=Double.parseDouble(properties.getProperty(prefix+"plTargetEigen"));
  			if (properties.getProperty(prefix+"plFractOutliers")!=null)   this.plFractOutliers=Double.parseDouble(properties.getProperty(prefix+"plFractOutliers"));
  			if (properties.getProperty(prefix+"plMaxOutliers")!=null)     this.plMaxOutliers=Integer.parseInt(properties.getProperty(prefix+"plMaxOutliers"));
  			if (properties.getProperty(prefix+"plMinStrength")!=null)     this.plMinStrength=Double.parseDouble(properties.getProperty(prefix+"plMinStrength"));
  			if (properties.getProperty(prefix+"plMaxEigen")!=null)        this.plMaxEigen=Double.parseDouble(properties.getProperty(prefix+"plMaxEigen"));
  			if (properties.getProperty(prefix+"plDbgMerge")!=null)        this.plDbgMerge=Boolean.parseBoolean(properties.getProperty(prefix+"plDbgMerge"));
  			if (properties.getProperty(prefix+"plWorstWorsening")!=null)  this.plWorstWorsening=Double.parseDouble(properties.getProperty(prefix+"plWorstWorsening"));
  			if (properties.getProperty(prefix+"plWeakWorsening")!=null)  this.plWeakWorsening=Double.parseDouble(properties.getProperty(prefix+"plWeakWorsening"));
  			if (properties.getProperty(prefix+"plMutualOnly")!=null)      this.plMutualOnly=Boolean.parseBoolean(properties.getProperty(prefix+"plMutualOnly"));

  			if (properties.getProperty(prefix+"plFillSquares")!=null)     this.plFillSquares=Boolean.parseBoolean(properties.getProperty(prefix+"plFillSquares"));
  			if (properties.getProperty(prefix+"plCutCorners")!=null)      this.plCutCorners=Boolean.parseBoolean(properties.getProperty(prefix+"plCutCorners"));

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
 
  			if (properties.getProperty(prefix+"dbg_migrate")!=null)       this.dbg_migrate=Boolean.parseBoolean(properties.getProperty(prefix+"dbg_migrate"));

  			if (properties.getProperty(prefix+"show_ortho_combine")!=null)     this.show_ortho_combine=Boolean.parseBoolean(properties.getProperty(prefix+"show_ortho_combine"));
  			if (properties.getProperty(prefix+"show_refine_supertiles")!=null) this.show_refine_supertiles=Boolean.parseBoolean(properties.getProperty(prefix+"show_refine_supertiles"));
  			if (properties.getProperty(prefix+"show_bgnd_nonbgnd")!=null)      this.show_bgnd_nonbgnd=Boolean.parseBoolean(properties.getProperty(prefix+"show_bgnd_nonbgnd"));
  			if (properties.getProperty(prefix+"show_filter_scan")!=null)       this.show_filter_scan=Boolean.parseBoolean(properties.getProperty(prefix+"show_filter_scan"));
  			if (properties.getProperty(prefix+"show_combined")!=null)          this.show_combined=Boolean.parseBoolean(properties.getProperty(prefix+"show_combined"));
  			if (properties.getProperty(prefix+"show_unique")!=null)            this.show_unique=Boolean.parseBoolean(properties.getProperty(prefix+"show_unique"));
  			if (properties.getProperty(prefix+"show_init_refine")!=null)       this.show_init_refine=Boolean.parseBoolean(properties.getProperty(prefix+"show_init_refine"));
  			if (properties.getProperty(prefix+"show_expand")!=null)            this.show_expand=Boolean.parseBoolean(properties.getProperty(prefix+"show_expand"));
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
  			gd.addCheckbox    ("Perfcorm coorrelation",                                                   this.correlate);
  			gd.addNumericField("itmask of pairs to combine in the composite (top, bottom, left,righth)",  this.corr_mask,            0);
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
  			gd.addCheckbox    ("Show result RGBA before overlap combined",                                this.show_nonoverlap);
  			gd.addCheckbox    ("Show result RGBA",                                                        this.show_overlap);
  			gd.addCheckbox    ("Show result color",                                                       this.show_rgba_color);
  			gd.addCheckbox    ("Show disparity maps",                                                     this.show_map);
  			gd.addNumericField("Disparity scan start value",                                              this.disp_scan_start,  2);
  			gd.addNumericField("Disparity scan step",                                                     this.disp_scan_step,  2);
  			gd.addNumericField("Disparity scan number of disparity values to scan",                       this.disp_scan_count,            0);

  			gd.addMessage("--- camera fine correction: X/Y for images 0..3  ---");
  			gd.addNumericField("X 0",                                                                     this.fine_corr_x_0,  3);
  			gd.addNumericField("Y 0",                                                                     this.fine_corr_y_0,  3);
  			gd.addNumericField("X 1",                                                                     this.fine_corr_x_1,  3);
  			gd.addNumericField("Y 1",                                                                     this.fine_corr_y_1,  3);
  			gd.addNumericField("X 2",                                                                     this.fine_corr_x_2,  3);
  			gd.addNumericField("Y 2",                                                                     this.fine_corr_y_2,  3);
  			gd.addNumericField("X 3",                                                                     this.fine_corr_x_3,  3);
  			gd.addNumericField("Y 4",                                                                     this.fine_corr_y_3,  3);

  			gd.addNumericField("Minimal correlation strength to apply fine correction",                   this.fcorr_min_stength,  3);
  			gd.addNumericField("Consider only tiles with absolute residual disparity lower than",         this.fcorr_disp_diff,  3);
  			gd.addCheckbox    ("Use quadratic polynomial for fine correction (false - only linear)",      this.fcorr_quadratic);
  			gd.addCheckbox    ("Ignore current calculated fine correction (use manual only)",             this.fcorr_ignore);
  			gd.addNumericField("Calculated from correlation offset vs. actual one (not yet understood)",  this.corr_magic_scale,  3);
  			
  			gd.addMessage     ("--- 3D reconstruction ---");
  			gd.addCheckbox    ("Show generated textures",                                                      this.show_textures);
  			gd.addCheckbox    ("show intermediate results of filtering",                                       this.debug_filters);

  			gd.addNumericField("Minimal noise-normalized pixel difference in a channel to suspect something",  this.min_smth,  3);
  			gd.addNumericField("Reliable noise-normalized pixel difference in a channel to have something ",   this.sure_smth,  3);
  			gd.addNumericField("Disparity range to be considered background",                                  this.bgnd_range,  3);
  			gd.addNumericField("Disparity difference from the center (provided) disparity to trust",           this.other_range,  3);
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
  			gd.addCheckbox    ("Correct lens geometric distortions in a model (will need backdrop to be corrected too)", this.correct_distortions);
  			gd.addCheckbox    ("Show generated triangles",                                                     this.show_triangles);
  			gd.addCheckbox    ("Weight-average disparity for the whole cluster ",                              this.avg_cluster_disp);
  			gd.addNumericField("Maximal disparity difference in a triangle face to show",                      this.maxDispTriangle,  6);

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
  			gd.addNumericField("Disaprity histogram step for far objects",                                     this.stStepFar,  6);
  			gd.addNumericField("Disaprity histogram step for near objects",                                    this.stStepNear,  6);
  			gd.addNumericField("Disaprity threshold to switch cfrom linear to logarithmic steps",              this.stStepThreshold,  6);
  			gd.addNumericField("Minimal disparity (center of a bin)",                                          this.stMinDisparity,  6);
  			gd.addNumericField("Maximal disparity (center of a bin)",                                          this.stMaxDisparity,  6);
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
  			
  			gd.addNumericField("Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert",this.stMeasSel,  0);
  			gd.addNumericField("Consider merging initial planes if disparity difference below",                this.stSmallDiff,  6);
  			gd.addNumericField("Consider merging initial planes if jumps between ratio above",                 this.stHighMix,  6);

  			gd.addNumericField("Outlayer tiles weaker than this may be replaced from neighbors",               this.outlayerStrength,  6);
  			gd.addNumericField("Replace weak outlayer tiles that do not have neighbors within this disparity difference", this.outlayerDiff,  6);
  			gd.addNumericField("Replace weak outlayer tiles that have higher disparity than weighted average", this.outlayerDiffPos,  6);
  			gd.addNumericField("Replace weak outlayer tiles that have lower disparity than weighted average",  this.outlayerDiffNeg,  6);

  			gd.addCheckbox    ("Combine with all previous after refine pass",                                  this.combine_refine);
  			gd.addNumericField("Disregard weaker tiles when combining scans",                                  this.combine_min_strength,  6);
  			gd.addNumericField("Disregard weaker tiles when combining scans  for horizontal correlation",      this.combine_min_hor,  6);
  			gd.addNumericField("Disregard weaker tiles when combining scans  for vertical correlation",        this.combine_min_vert,  6);
  			gd.addNumericField("Do not re-measure correlation if target disparity differs from some previous by this",this.unique_tolerance,  6);
  			
  			gd.addMessage     ("--- Growing disparity range to scan ---");
  			gd.addNumericField("Try these number of tiles around known ones",                                  this.grow_sweep,  0);
  			gd.addNumericField("Maximal disparity to try",                                                     this.grow_disp_max,  6);
  			gd.addNumericField("Trust measured disparity within +/- this value",                               this.grow_disp_trust,  6);
  			gd.addNumericField("Increase disparity (from maximal tried) if nothing found in that tile",        this.grow_disp_step,  6);
  			gd.addNumericField("Grow more only if at least one channel has higher variance from others for the tile", this.grow_min_diff,  6);
  			gd.addMessage     ("--- Planes detection ---");
  			gd.addCheckbox    ("Always start with disparity-most axis (false - lowest eigenvalue)",            this.plPreferDisparity);
  			gd.addNumericField("Normalize disparities to the average if above",                                this.plDispNorm,  6);
  			gd.addNumericField("Minimal number of points for plane detection",                                  this.plMinPoints,  0);
  			gd.addNumericField("Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below", this.plTargetEigen,  6);
  			gd.addNumericField("Maximal fraction of outliers to remove",                                       this.plFractOutliers,  6);
  			gd.addNumericField("Maximal number of outliers to remove",                                         this.plMaxOutliers,  0);
  			gd.addNumericField("Minimal total strength of a plane",                                            this.plMinStrength,  6);
  			gd.addNumericField("Maximal eigenvalue of a plane",                                                this.plMaxEigen,  6);
  			gd.addCheckbox    ("Combine 'other' plane with the current (unused)",                              this.plDbgMerge);
  			gd.addNumericField("Worst case worsening after merge",                                             this.plWorstWorsening,  6);
  			gd.addNumericField("Relax merge requirements for weaker planes",                                   this.plWeakWorsening,  6);
  			gd.addCheckbox    ("Keep only mutual links, remove weakest if conflict",                           this.plMutualOnly);

  			gd.addCheckbox    ("Add diagonals to full squares",                                                this.plFillSquares);
  			gd.addCheckbox    ("Add ortho to 45-degree corners",                                               this.plCutCorners);

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
  			
  			gd.addCheckbox    ("Test new mode after migration",                                                this.dbg_migrate);

  			gd.addMessage     ("--- Other debug images ---");
  			gd.addCheckbox    ("Show 'ortho_combine'",                                                         this.show_ortho_combine);
  			gd.addCheckbox    ("Show 'refine_disparity_supertiles'",                                           this.show_refine_supertiles);
  			gd.addCheckbox    ("Show 'bgnd_nonbgnd'",                                                          this.show_bgnd_nonbgnd);
  			gd.addCheckbox    ("Show 'FilterScan'",                                                            this.show_filter_scan);
  			gd.addCheckbox    ("Show 'combo_scan' (combined multiple scans)",                                  this.show_combined);
  			gd.addCheckbox    ("Show 'unique_scan' (removed already measured tiles with the same disparity)",  this.show_unique);
  			gd.addCheckbox    ("Show debug images during initial refinement",                                  this.show_init_refine);
  			gd.addCheckbox    ("Show debug images during disparity expansion",                                 this.show_expand);
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
  			this.black_back=          gd.getNextBoolean();
  			this.keep_weights=          gd.getNextBoolean();
  			this.sharp_alpha=           gd.getNextBoolean();
  			this.alpha0=                gd.getNextNumber();
  			this.alpha1=                gd.getNextNumber();
  			this.gen_chn_stacks=        gd.getNextBoolean();
  			this.gen_chn_img=           gd.getNextBoolean();
  			this.show_nonoverlap=       gd.getNextBoolean();
  			this.show_overlap=          gd.getNextBoolean();
  			this.show_rgba_color=       gd.getNextBoolean();
  			this.show_map=              gd.getNextBoolean();
  			this.disp_scan_start=       gd.getNextNumber();
  			this.disp_scan_step=        gd.getNextNumber();
  			this.disp_scan_count= (int) gd.getNextNumber();
  			
  			this.fine_corr_x_0=         gd.getNextNumber();
  			this.fine_corr_y_0=         gd.getNextNumber();
  			this.fine_corr_x_1=         gd.getNextNumber();
  			this.fine_corr_y_1=         gd.getNextNumber();
  			this.fine_corr_x_2=         gd.getNextNumber();
  			this.fine_corr_y_2=         gd.getNextNumber();
  			this.fine_corr_x_3=         gd.getNextNumber();
  			this.fine_corr_y_3=         gd.getNextNumber();
  			
  			this.fcorr_min_stength=     gd.getNextNumber();
  			this.fcorr_disp_diff=       gd.getNextNumber();
  			this.fcorr_quadratic=       gd.getNextBoolean();
  			this.fcorr_ignore=          gd.getNextBoolean();
  			this.corr_magic_scale=      gd.getNextNumber();

  			this.show_textures=         gd.getNextBoolean();
  			this.debug_filters=         gd.getNextBoolean();
  			this.min_smth=              gd.getNextNumber();
  			this.sure_smth=             gd.getNextNumber();
  			this.bgnd_range=            gd.getNextNumber();
  			this.other_range=           gd.getNextNumber();
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
  			this.correct_distortions=   gd.getNextBoolean();
  			this.show_triangles=        gd.getNextBoolean();
  			this.avg_cluster_disp=      gd.getNextBoolean();
  			this.maxDispTriangle=       gd.getNextNumber();
  			this.shUseFlaps=            gd.getNextBoolean();
  			this.shAggrFade=            gd.getNextBoolean();
  			this.shMinArea=       (int) gd.getNextNumber();
  			this.shMinStrength=        gd.getNextNumber();
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
  			this.stMaxDisparity=        gd.getNextNumber();
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
  			this.unique_tolerance=      gd.getNextNumber();

  			this.grow_sweep=      (int) gd.getNextNumber();
  			this.grow_disp_max=         gd.getNextNumber();
  			this.grow_disp_trust=       gd.getNextNumber();
  			this.grow_disp_step=        gd.getNextNumber();
  			this.grow_min_diff=         gd.getNextNumber();

  			this.plPreferDisparity=     gd.getNextBoolean();
  			this.plDispNorm=            gd.getNextNumber();
  			this.plMinPoints=     (int) gd.getNextNumber();
  			this.plTargetEigen=         gd.getNextNumber();
  			this.plFractOutliers=       gd.getNextNumber();
  			this.plMaxOutliers=   (int) gd.getNextNumber();
  			this.plMinStrength=         gd.getNextNumber();
  			this.plMaxEigen=            gd.getNextNumber();
  			this.plDbgMerge=            gd.getNextBoolean();
  			this.plWorstWorsening=      gd.getNextNumber();
  			this.plWeakWorsening=       gd.getNextNumber();
  			this.plMutualOnly=          gd.getNextBoolean();

  			this.plFillSquares=         gd.getNextBoolean();
  			this.plCutCorners=          gd.getNextBoolean();

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

  			this.dbg_migrate=           gd.getNextBoolean();

  			this.show_ortho_combine=    gd.getNextBoolean();
  			this.show_refine_supertiles=gd.getNextBoolean();
  			this.show_bgnd_nonbgnd=     gd.getNextBoolean(); // first on second pass
  			this.show_filter_scan=      gd.getNextBoolean(); // first on refine
  			this.show_combined=         gd.getNextBoolean();
  			this.show_unique=           gd.getNextBoolean();
  			this.show_init_refine=      gd.getNextBoolean();
  			this.show_expand=           gd.getNextBoolean();
  			this.show_shells=           gd.getNextBoolean();
  			this.show_neighbors=        gd.getNextBoolean();
  			this.show_flaps_dirs=       gd.getNextBoolean();
  			this.show_first_clusters=   gd.getNextBoolean();
  			this.show_planes=           gd.getNextBoolean();
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
		
		public boolean matchPixelSize=    true; // disregard next value, calculate projectionPixelSize from teh equirectangular map
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
