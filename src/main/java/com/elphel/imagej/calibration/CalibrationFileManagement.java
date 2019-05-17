package com.elphel.imagej.calibration;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

import ij.IJ;
import ij.io.OpenDialog;


public class CalibrationFileManagement {
	public static String DEFAULT_DIRECTORY=null;
/*
	public String getLatestName(FileFilter filter, String dirPath){
		File dFile=new File(dirPath);
		if (!dFile.isDirectory()) {
			String msg="Could not find directory with saved SFE configuration: "+dirPath;
//			IJ.showMessage(msg);
			System.out.println("Error: "+msg);
			return null;
		}
		File [] files=dFile.listFiles(filter);


	}
	*/
/* ======================================================================== */
	  public static String selectDirectory(boolean save, String title, String button, FileFilter filter, String defaultPath) {
		  return selectDirectoryOrFile(false, save,true, title, button, filter,defaultPath);  // always open dialog
	  }
	  /**
	   *
	   * @param smart        if true, and defaultPath matches criteria, return it and do not open dialog
	   * @param save         file/directory for writing, new OK
	   * @param title        Dialog title
	   * @param button       Name on the Button
	   * @param filter       Selection filter
	   * @param defaultPath  default path
	   * @return             selected directory name
	   */
	  public static String selectDirectory(boolean smart, boolean save, String title, String button, FileFilter filter, String defaultPath) {
		  return selectDirectoryOrFile(smart, save,true, title, button, filter,defaultPath);
	  }
	  public static String [] selectDirectories(boolean save, String title, String button, FileFilter filter, String [] defaultPaths) {
		  return selectDirectoriesOrFiles(save,true, title, button, filter, defaultPaths);
	  }
	  public static String selectFile(boolean save,  String title, String button, FileFilter filter, String defaultPath) {
		  return  selectDirectoryOrFile(false, save,false, title, button, filter, defaultPath ); // always open dialog
	  }
	  public static String selectFile(boolean smart, boolean save,  String title, String button, FileFilter filter, String defaultPath) {
		  return  selectDirectoryOrFile(smart, save,false, title, button, filter, defaultPath );
	  }
	  public static String [] selectFiles(boolean save,  String title, String button, FileFilter filter, String [] defaultPaths) {
		  return  selectDirectoriesOrFiles(save,false, title, button, filter, defaultPaths );
	  }
/* ======================================================================== */
	  public static String [] selectDirectoriesOrFiles(boolean save,
			  boolean directory,
			  String title,
			  String button,
			  FileFilter filter,
			  String [] defaultPaths) {
		  File dir=null;
		  String defaultPath=null;
		  File [] files=null;
		  int fileNum;
		  if ((defaultPaths!=null) && (defaultPaths.length>0)) {
			  File [] tfiles=new File [defaultPaths.length];
			  int nf=defaultPaths.length;
			  for (fileNum=0;fileNum<defaultPaths.length; fileNum++) {
				  tfiles[fileNum]=new File(defaultPaths[fileNum]);
				  if ((!tfiles[fileNum].exists()) ||(!tfiles[fileNum].isFile())) {
					  tfiles[fileNum]=null;
					  nf--;
				  }
			  }
			  files=new File[nf];
			  nf=0;
			  for (fileNum=0;fileNum<defaultPaths.length; fileNum++) if (tfiles[fileNum]!=null){
				  files[nf++]=tfiles[fileNum];
			  }
		  }
		  if ((defaultPaths!=null) && (defaultPaths.length>0) &&  (!defaultPaths[0].equals(""))) {
			  defaultPath=defaultPaths[0];
			  dir = new File(defaultPath);
		  }
		  if ((dir==null) || (!dir.exists())) {
			  if (DEFAULT_DIRECTORY!=null) {
				  defaultPath = DEFAULT_DIRECTORY;
				  dir = new File(defaultPath);
			  }
		  }
		  if ((dir==null) || (!dir.exists())) {
			  defaultPath = OpenDialog.getDefaultDirectory();
			  if (defaultPath!=null) dir = new File(defaultPath);
		  }
		  if ((dir!=null) && (!dir.exists())) dir=null;
		  if ((dir!=null) && (!dir.isDirectory())){
			  dir=dir.getParentFile();
		  }
	//getSelectedFiles

		  JFileChooser fc= new JFileChooser();
		  fc.setFileSelectionMode(directory?JFileChooser.DIRECTORIES_ONLY:JFileChooser.FILES_ONLY);
		  fc.setMultiSelectionEnabled(true);
		  if ((title!=null)  && (title.length()>0)) fc.setDialogTitle(title);
		  if ((button!=null) && (button.length()>0)) fc.setApproveButtonText(button);
		  if (filter!=null) fc.setFileFilter(filter) ;
		  if (dir!=null) 	fc.setCurrentDirectory(dir);
		  fc.setSelectedFiles(files);
		  int returnVal = save?(fc.showSaveDialog(IJ.getInstance())):(fc.showOpenDialog(IJ.getInstance()));
		  if (returnVal!=JFileChooser.APPROVE_OPTION)	return null;
		  DEFAULT_DIRECTORY=fc.getCurrentDirectory().getPath();
		  files=fc.getSelectedFiles();
		  if (files.length<1) return null;
		  String [] filenames=new String[files.length];
//		  for (int nFile=0;nFile<files.length;nFile++) filenames[nFile]= files[nFile].getName();
		  for (int nFile=0;nFile<files.length;nFile++) filenames[nFile]= files[nFile].getPath();
		  return filenames;
	  }

	  public static String selectDirectoryOrFile(
			  boolean smart, // do not open dialog if defaultPath matches filter
			  boolean save,
			  boolean directory,
			  String title,
			  String button,
			  FileFilter filter,
			  String defaultPath) {
		  File dir=null;
		  if ((defaultPath!=null) &&  (defaultPath.length()>0)) {
			  dir = new File(defaultPath);
		  }
		  // If directory is specified, smart=true, save is enabled, but it does not exist - try to create it
		  if (smart &&
				  directory &&
				  (dir!=null) &&
				  (defaultPath.length()>1) && // skip "/"
				  save &&
				  !dir.exists()) dir.mkdirs();


		  // see if defaultPath matches
		  if (smart &&
				  (dir!=null) &&
				  (defaultPath.length()>1) && // skip "/"
				  (dir.exists()) &&
				  (dir.isDirectory() ^ (!directory)) && // directory if requested directory, file if requested file
				  (dir.isDirectory() || filter.accept(dir))){ // don't care for directory, match filter if file
			  return defaultPath;
		  }

		  if ((dir==null) || (!dir.exists())) {
			  if (DEFAULT_DIRECTORY!=null) {
				  defaultPath = DEFAULT_DIRECTORY;
				  dir = new File(defaultPath);
			  }
		  }
		  if ((dir==null) || (!dir.exists())) {
			  defaultPath = OpenDialog.getDefaultDirectory();
			  if (defaultPath!=null) dir = new File(defaultPath);
		  }
		  if ((dir!=null) && (!dir.exists())) dir=null;
		  if ((dir!=null) && (!dir.isDirectory())){
			  dir=dir.getParentFile();
		  }


		  JFileChooser fc= new JFileChooser();
		  fc.setFileSelectionMode(directory?JFileChooser.DIRECTORIES_ONLY:JFileChooser.FILES_ONLY);
		  fc.setMultiSelectionEnabled(false);
		  if ((title!=null)  && (title.length()>0)) fc.setDialogTitle(title);
		  if ((button!=null) && (button.length()>0)) fc.setApproveButtonText(button);
		  if (filter!=null) fc.setFileFilter(filter) ;
		  if (dir!=null) 	fc.setCurrentDirectory(dir);
		  int returnVal = save?(fc.showSaveDialog(IJ.getInstance())):(fc.showOpenDialog(IJ.getInstance()));
		  if (returnVal!=JFileChooser.APPROVE_OPTION)	return null;
		  DEFAULT_DIRECTORY=fc.getCurrentDirectory().getPath();
		  return fc.getSelectedFile().getPath();
	  }

	  public static void saveStringToFile (String path,String data){
		  BufferedWriter writer = null;
		  try  {
			  writer = new BufferedWriter( new FileWriter( path));
			  writer.write( data);

		  }  catch ( IOException e)  {
				String msg = e.getMessage();
				if (msg==null || msg.equals(""))  msg = ""+e;
				IJ.showMessage("Error",msg);
				throw new IllegalArgumentException (msg);
		  }
		  finally  {
			  try {
				  if ( writer != null)
					  writer.close( );
			  } catch ( IOException e) {

			  }
		  }
	  }

/* ======================================================================== */
	  public static class MultipleExtensionsFileFilter extends FileFilter implements FilenameFilter {
		  protected String [] patterns; // case insensitive
		  protected String    description="JP4 files";
		  protected String    prefix=""; // case sensitive

		  public MultipleExtensionsFileFilter (String prefix, String [] patterns,String description) {
			  this.prefix=     prefix;
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
			  if (!name.startsWith(this.prefix)) return false; // empty prefix OK
			  for (i=0;i<patterns.length;i++) {
				  if (name.toLowerCase().endsWith(patterns[i].toLowerCase())) return true;
			  }
			  return false;
		  }
		  @Override
		public boolean accept (File dir, String name) { // directory - don't care here, only name
			  if (!name.startsWith(this.prefix)) return false; // empty prefix OK
			  for (int i=0;i<patterns.length;i++) {
				  if (name.toLowerCase().endsWith(patterns[i].toLowerCase())) return true;
			  }
			  return false;
		  }

		  @Override
		public String getDescription() {
			  return description;
		  }
	  }

}
