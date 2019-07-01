package com.elphel.imagej.calibration;

import java.io.File;

// select directory that contains matching MultipleExtensionsFileFilter in specified min/max quantities
//https://stackoverflow.com/questions/22302199/java-filefilter-to-select-certain-directories
//	  public static class DirectoryContentsFilter extends FileFilter implements FilenameFilter {
	  public class DirectoryChoser extends javax.swing.JFileChooser {
		  private static final long serialVersionUID = 390855361964415146L;
		  protected MultipleExtensionsFileFilter multipleExtensionsFileFilter;
		  protected int    min_files;
		  protected int    max_files;
		  protected String description;

		  public DirectoryChoser (
				  MultipleExtensionsFileFilter multipleExtensionsFileFilter,
				  int                          min_files,
				  int                          max_files,
				  String                       description) // may be null;
		  {
			  this.multipleExtensionsFileFilter=   multipleExtensionsFileFilter;
			  this.min_files = min_files;
			  this.max_files = max_files;
			  this.description = description;
		  }
		  @Override
		  public boolean isDirectorySelectionEnabled() {
//			  setOpenButtonState(this, false);
			  File file = getSelectedFile();
//			  File [] files = getSelectedFiles();
			  if(file == null){
//				  setOpenButtonState(this, true);
				  return  true; // false;
			  }
//			  setOpenButtonState(this, false);
			  if (!file.isDirectory()) return false;
/*
	// Can not make it work correctly with multiple selection, giving up for now
			  // get a list of all matching files

			  int num_match = file.list(multipleExtensionsFileFilter).length;
			  if (num_match < min_files) return false;
			  if ((max_files > 0) && (num_match > max_files)) {
				  return false;
			  }
*/
			  setOpenButtonState(this, true);
			  return true;
		  }

		    private void setOpenButtonState(java.awt.Container c, boolean flag) {
		        int len = c.getComponentCount();
		        for (int i = 0; i < len; i++) {
		            java.awt.Component comp = c.getComponent(i);

		            if (comp instanceof javax.swing.JButton) {
		                javax.swing.JButton b = (javax.swing.JButton)comp;

		                if ( b != null && b.getText() != null && b.getText().equals("Select") ) {
		                    b.setEnabled(flag);
		                }

		            } else if (comp instanceof java.awt.Container) {
		                setOpenButtonState((java.awt.Container) comp, flag);
		            }
		        }
		    }


/*
		  @Override
		  public boolean accept (File file) {
			  if (!file.isDirectory()) return false;
			  // get a list of all matching files
			  int num_match = file.list(multipleExtensionsFileFilter).length;
			  if (num_match < min_files) return false;
			  if ((max_files > 0) && (num_match > max_files)) return false;
			  return true;
		  }
		  @Override
		  public boolean accept (File dir, String name) {
			  // No name filter, need to resolve
			  //			  Path dir_path = dir.toPath();
			  File target = Paths.get(dir.toString(),name).toFile();
			  return accept(target);
		  }

		  @Override
		  public String getDescription() { // if description is null - generate it
			  if (description == null) {
				  description = "Directories containing";
				  if (min_files > 0)  description +=" not less than "+min_files;
				  if (max_files > 0)  description +=" not more than "+max_files;
				  description += " "+multipleExtensionsFileFilter.description;

			  }
			  return description;
		  }
*/
	  }