package com.elphel.imagej.calibration;

import java.io.File;
import java.io.FilenameFilter;

import javax.swing.filechooser.FileFilter;

/* ======================================================================== */
  public class MultipleExtensionsFileFilter extends FileFilter implements FilenameFilter {
	  protected String [] patterns; // case insensitive
	  protected String    description="JP4 files";
	  protected String [] prefixes= {""}; // case sensitive


	  // prefix/prefixeas can not be null - ambiguity
	  public MultipleExtensionsFileFilter (String prefix, String [] patterns,String description) {
		  this.prefixes=     new String[1];
		  this.prefixes[0] = prefix;
		  this.description=description;
		  this.patterns=   patterns.clone();
	  }
	  public MultipleExtensionsFileFilter (String [] prefixes, String [] patterns,String description) {
		  this.prefixes=   prefixes;
		  this.description=description;
		  this.patterns=   patterns.clone();
	  }

	  public MultipleExtensionsFileFilter (String [] patterns,String description) {
		  this.description=description;
		  this.patterns=patterns.clone();
	  }
	  public MultipleExtensionsFileFilter (String [] patterns) {
		  this.patterns=patterns.clone();
	  }
	  @Override
	  public boolean accept (File file) {
		  int i;
		  String name=file.getName();
		  if (file.isDirectory()) return true;
		  if (!testPrefix(name)) return false;
		  for (i=0;i<patterns.length;i++) {
			  if (name.toLowerCase().endsWith(patterns[i].toLowerCase())) return true;
		  }
		  return false;
	  }
	  @Override
	public boolean accept (File dir, String name) { // directory - don't care here, only name
		  if (!testPrefix(name)) return false;
		  for (int i=0;i<patterns.length;i++) {
			  if (name.toLowerCase().endsWith(patterns[i].toLowerCase())) return true;
		  }
		  return false;
	  }
	  @Override
	public String getDescription() {
		  return description;
	  }
	  private boolean testPrefix(String name) {
		  if ((prefixes == null) || (prefixes[0] == null)) return true;
		  for (String s : prefixes) {
			  if ((s == null) || name.startsWith(s)) {
				  return true;
			  }
		  }
		  return false; // prefix mismatch;
	  }
  }