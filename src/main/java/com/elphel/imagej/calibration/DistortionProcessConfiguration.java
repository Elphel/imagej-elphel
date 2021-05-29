package com.elphel.imagej.calibration;
/*
 **
 ** DistortionProcessConfiguration.java
 **
 ** Copyright (C) 2011-2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  DistortionProcessConfiguration.java is free software: you can redistribute it and/or modify
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

import java.io.File;
import java.io.FilenameFilter;
import java.util.Arrays;
import java.util.Properties;

import com.elphel.imagej.common.WindowTools;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;


    public class DistortionProcessConfiguration{
    	public String  sourceDirectory="";
    	public String  gridDirectory=  "";
    	public boolean useLaserPonters=true;
    	public boolean useNoPonters=true; // use images that do not have any lasre pointers
    	public boolean showAcquiredImages=true;
    	public boolean saveAcquiredImages=true;
    	public boolean selectSourceFiles = true;
    	public boolean removeOutOfGridPointers=true;
    	public boolean showGridImages=false;
    	public boolean saveGridImages=true;
    	public boolean overwriteResultFiles=false;
    	public int     debugLevel=1;

    	public String selectSourceDirectory(boolean smart, String defaultPath, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Source (acquired from the camera) image directory", // title
    				"Select source directory", // button
    				null, // filter
    				defaultPath); // this.sourceDirectory);
    		if (dir!=null) this.sourceDirectory=dir;
    		return dir;
    	}
    	public String selectGridFileDirectory(boolean smart, String defaultPath, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save
    				"Grid files directory (grid patterns extracted from the images)", // title
    				"Select grid files directory", // button
    				null, // filter
    				defaultPath); //this.sourceDirectory);
    		if (dir!=null) this.gridDirectory=dir;
    		return dir;
    	}

    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"sourceDirectory",        this.sourceDirectory);
    		properties.setProperty(prefix+"gridDirectory",          this.gridDirectory);
    		properties.setProperty(prefix+"useLaserPonters",        this.useLaserPonters+"");
    		properties.setProperty(prefix+"useNoPonters",           this.useNoPonters+"");
    		properties.setProperty(prefix+"showAcquiredImages",     this.showAcquiredImages+"");
    		properties.setProperty(prefix+"saveAcquiredImages",     this.saveAcquiredImages+"");
    		properties.setProperty(prefix+"selectSourceFiles",      this.selectSourceFiles+"");
    		properties.setProperty(prefix+"removeOutOfGridPointers",this.removeOutOfGridPointers+"");
    		properties.setProperty(prefix+"showGridImages",         this.showGridImages+"");
    		properties.setProperty(prefix+"saveGridImages",         this.saveGridImages+"");
    		properties.setProperty(prefix+"overwriteResultFiles",   this.overwriteResultFiles+"");
    		properties.setProperty(prefix+"debugLevel",             this.debugLevel+"");
    	}
       	public void getProperties(String prefix,Properties properties){
    		if (properties.getProperty(prefix+"sourceDirectory")!=null)
    			this.sourceDirectory=properties.getProperty(prefix+"sourceDirectory");
    		if (properties.getProperty(prefix+"gridDirectory")!=null)
    			this.gridDirectory=properties.getProperty(prefix+"gridDirectory");
    		if (properties.getProperty(prefix+"useLaserPonters")!=null)
    			this.useLaserPonters=Boolean.parseBoolean(properties.getProperty(prefix+"useLaserPonters"));
    		if (properties.getProperty(prefix+"useNoPonters")!=null)
    			this.useNoPonters=Boolean.parseBoolean(properties.getProperty(prefix+"useNoPonters"));
    		if (properties.getProperty(prefix+"showAcquiredImages")!=null)
    			this.showAcquiredImages=Boolean.parseBoolean(properties.getProperty(prefix+"showAcquiredImages"));
    		if (properties.getProperty(prefix+"saveAcquiredImages")!=null)
    			this.saveAcquiredImages=Boolean.parseBoolean(properties.getProperty(prefix+"saveAcquiredImages"));
    		if (properties.getProperty(prefix+"selectSourceFiles")!=null)
    			this.selectSourceFiles=Boolean.parseBoolean(properties.getProperty(prefix+"selectSourceFiles"));
    		if (properties.getProperty(prefix+"removeOutOfGridPointers")!=null)
    			this.removeOutOfGridPointers=Boolean.parseBoolean(properties.getProperty(prefix+"removeOutOfGridPointers"));
    		if (properties.getProperty(prefix+"showGridImages")!=null)
    			this.showGridImages=Boolean.parseBoolean(properties.getProperty(prefix+"showGridImages"));
    		if (properties.getProperty(prefix+"saveGridImages")!=null)
    			this.saveGridImages=Boolean.parseBoolean(properties.getProperty(prefix+"saveGridImages"));
    		if (properties.getProperty(prefix+"overwriteResultFiles")!=null)
    			this.overwriteResultFiles=Boolean.parseBoolean(properties.getProperty(prefix+"overwriteResultFiles"));
    		if (properties.getProperty(prefix+"debugLevel")!=null)
    			this.debugLevel=Integer.parseInt(properties.getProperty(prefix+"debugLevel"));
       	}
    	public boolean showDialog(String title) {
    		GenericDialog gd = new GenericDialog(title);
    		gd.addStringField("Source (acquired from the camera) image directory, blank will open selection window",this.sourceDirectory,40);
    		gd.addStringField("Grid files directory (grid patterns extracted from the images, blank will open selection window)",this.gridDirectory,40);
    		gd.addCheckbox   ("Locate laser pointers for each image", this.useLaserPonters);
    		gd.addCheckbox   ("Use images that do not contain laser pointers", this.useLaserPonters);
    		gd.addCheckbox   ("Show images after acquisition",        this.showAcquiredImages);
    		gd.addCheckbox   ("Save acquired images",                 this.saveAcquiredImages);
    		gd.addCheckbox   ("Individually select source image (false use all directory)", this.selectSourceFiles);
    		gd.addCheckbox   ("Remove detected laser pointers if they are outside of the grid", this.removeOutOfGridPointers);
    		gd.addCheckbox   ("Show grid files as images",            this.showGridImages);
    		gd.addCheckbox   ("Save grid files",                      this.saveGridImages);
    		gd.addCheckbox   ("Overwrite existing result files",      this.overwriteResultFiles);



    		gd.addNumericField("Debug level",                         this.debugLevel,0);
    	    WindowTools.addScrollBars(gd);
    	    gd.showDialog();
    	    if (gd.wasCanceled()) return false;
    	    String newSourceDirectory=   gd.getNextString();
    	    String newGridDirectory=     gd.getNextString();
    	    this.useLaserPonters=        gd.getNextBoolean();
    		this.useLaserPonters=        gd.getNextBoolean();
    	    this.showAcquiredImages=     gd.getNextBoolean();
    	    this.saveAcquiredImages=     gd.getNextBoolean();
    	    this.selectSourceFiles=      gd.getNextBoolean();
    	    this.removeOutOfGridPointers=gd.getNextBoolean();
    	    this.showGridImages=         gd.getNextBoolean();
    	    this.saveGridImages=         gd.getNextBoolean();
    	    this.overwriteResultFiles=   gd.getNextBoolean();
    	    this.debugLevel=       (int) gd.getNextNumber();
    	    System.out.println("1.newSourceDirectory = "+newSourceDirectory);
    	    System.out.println("1.newGridDirectory = "+  newGridDirectory);
    	    if ((newSourceDirectory.length()==0) || (newSourceDirectory.indexOf('?')>=0))
    	    	newSourceDirectory= selectSourceDirectory(false, this.sourceDirectory, true);
    	    else
    	    	newSourceDirectory= selectSourceDirectory(true, newSourceDirectory, true); // if matches, no dialog
    	    if (newSourceDirectory!=null) this.sourceDirectory=newSourceDirectory;
    	    if ((newGridDirectory.length()==0) || (newGridDirectory.indexOf('?')>=0))
    	    	newGridDirectory= selectGridFileDirectory(false, this.gridDirectory, true);
    	    else
    	    	newGridDirectory= selectGridFileDirectory(true, newGridDirectory, true);
    	    if (newGridDirectory!=null) this.gridDirectory=newGridDirectory;
//    	    System.out.println("2.newSourceDirectory = "+newSourceDirectory);
//    	    System.out.println("2.newGridDirectory = "+  newGridDirectory);
//    	    System.out.println("this.sourceDirectory = "+this.sourceDirectory);
//    	    System.out.println("this.gridDirectory = "+  this.gridDirectory);
    	    return true;
    	}

    	public String [] selectSourceFiles(){
    		return selectSourceFiles(!this.selectSourceFiles);
    	}
    	public String [] selectSourceFiles(boolean allFiles){
			String [] extensions={".tif",".tiff"};
			if (this.sourceDirectory.length()==0){
    	    	String newSourceDirectory= selectSourceDirectory(true, this.sourceDirectory, true);
    	    	if (newSourceDirectory!=null) this.sourceDirectory=newSourceDirectory;
			}
			String [] defaultPaths={this.sourceDirectory+Prefs.getFileSeparator()};
			if (this.sourceDirectory.length()==0){
				defaultPaths[0]="";
			}

			MultipleExtensionsFileFilter sourceFilter =
				new MultipleExtensionsFileFilter("",extensions,"Source files");
			String [] sourceFiles=null;
    		if (allFiles){
				File dir= new File (this.sourceDirectory);
				if (this.debugLevel>1) System.out.println("selectSourceFiles, dir="+this.sourceDirectory);
				if (!dir.exists()) {
					String error="Source directory "+this.sourceDirectory+" does not exist.";
            		IJ.showMessage("No files selected");
					if (this.debugLevel>1) System.out.println("selectSourceFiles() ERROR:"+error);
					return null;
				}
				File [] fileList=dir.listFiles(sourceFilter);
				if (this.debugLevel>1) System.out.println("Source directory "+this.sourceDirectory+" has "+fileList.length+" files.");
				sourceFiles = new String[fileList.length];
				for (int i=0;i<sourceFiles.length;i++) sourceFiles[i]=fileList[i].getPath();
    		} else {
    				new MultipleExtensionsFileFilter("",extensions,"Source files");
    			sourceFiles=CalibrationFileManagement.selectFiles(false,
    					"Select Source files, saved as TIFF",
    					"Select",
    					sourceFilter,
    					defaultPaths); // String [] defaultPaths); //this.sourceDirectory // null
    	       	if ((sourceFiles==null) || (sourceFiles.length==0)) {
					if (this.debugLevel>1) System.out.println("selectSourceFiles() ERROR: No files selected");
            		IJ.showMessage("No files selected");
            		return null;
            	}
    		}
	       	return sourceFiles;
    	}

        public String[] selectSourceSets() {
			File dir= new File (this.sourceDirectory);
			if (this.debugLevel>1) System.out.println("selectSourceSets, dir="+this.sourceDirectory);
			if (!dir.exists()) {
				String error="Source directory "+this.sourceDirectory+" does not exist.";
        		IJ.showMessage("No files selected");
				if (this.debugLevel>1) System.out.println("selectSourceFiles() ERROR:"+error);
				return null;
			}
			File [] sourceFileSets = dir.listFiles(new FilenameFilter() {
			  @Override
			  public boolean accept(File current, String name) {
			    return new File(current, name).isDirectory();
			  }
			});
			// TODO:Sort set list
			String [] sourceSets = new String[sourceFileSets.length];
			for (int i=0;i<sourceSets.length;i++) sourceSets[i]=sourceFileSets[i].getPath();
			Arrays.sort(sourceSets);
			return sourceSets;
        }

    }

