package com.elphel.imagej.correction;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;

import com.elphel.imagej.calibration.CalibrationIllustration;
import com.elphel.imagej.calibration.CalibrationIllustrationParameters;
import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.common.GenericJTabbedDialog;

import ij.IJ;
import ij.Prefs;

public class FootageOrganize {
//	static String [] copy_modes = {"copy","move", "soft link", "hard link"};
//	int num_folders = 5;
//	int copy_mode =   0;
//	String [] sourceFolders = new String[num_folders];
//	String destinationFolder = "";
//	double inter_gap = 1.0; // seconds between series
	
	
	
	public boolean OrganizeSeries(
			CLTParameters   clt_parameters,
			CalibrationIllustration illustration) {
		CalibrationIllustrationParameters illustrationParameters = illustration.getIllustrationParameters();
		GenericJTabbedDialog gd = new GenericJTabbedDialog("Select directories to process",1000,1000);
		for (int si = 0; si < illustrationParameters.num_source_folders; si++) {
			gd.addStringField ("Source directory "+(si+1),     illustrationParameters.sourceFolders[si], 100,
					"Directory containing timestamped scene folders. Keep empty if not needed (e.g. multiple cameras already merged");
		}
		gd.addStringField ("Destination directory",     illustrationParameters.destinationFolder, 100,
				"Will contain folders for each sequence of (usually) 100 scene folders");
		gd.addChoice("Copy/move/link mode", CalibrationIllustrationParameters.copy_modes, CalibrationIllustrationParameters.copy_modes[illustrationParameters.copy_mode],
				"Copy mode. Hard links are possible only for the same disk");

  		gd.addNumericField("Minimal gap between series",           illustrationParameters.inter_gap, 4, 6,"s",
				"Series will be split if the time gap between scenes exceeds this value");
		gd.addNumericField("Sequence length",                             illustrationParameters.seq_len,     0,3,"scenes",
				"Limit sequence length if >=0 ");
    	gd.addCheckbox    ("Two-level directories",                       illustrationParameters.two_level_dirs,
    			"Consolidate scene files in sequences named by the first timestamp. Requires at least one of inter_gap or seq_len to be > 0");

    	gd.addCheckbox    ("Ignore scenes with partial LWIR channels",    illustrationParameters.isCaptures_all_lwir(), "Skip scenes where not all or none LWIR channels are available");
    	gd.addCheckbox    ("Ignore scenes with partial EO channels",      illustrationParameters.isCaptures_all_eo(),   "Skip scenes where not all or none EO channels are available");
    	gd.addCheckbox    ("Ignore partial scenes",                       illustrationParameters.isCaptures_all(),      "Skip scenes where not all (20) images are present");
    	
    	gd.addCheckbox    ("*** Check me to proceed with files ! ***",    false,
    			"Just a safety measure. Select correct link/move/copy mode and available disk space. Some files will not be processed (e.g.incomplete scenes)");
		
    	
		
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		for (int si = 0; si < illustrationParameters.num_source_folders; si++) {
			illustrationParameters.sourceFolders[si]=  gd.getNextString();
		}
		illustrationParameters.destinationFolder=      gd.getNextString();
		illustrationParameters.copy_mode =             gd.getNextChoiceIndex();
		illustrationParameters.inter_gap =             gd.getNextNumber();
		illustrationParameters.seq_len=          (int) gd.getNextNumber();
		illustrationParameters.two_level_dirs =        gd.getNextBoolean();
		
		illustrationParameters.setCaptures_all_lwir(gd.getNextBoolean());
		illustrationParameters.setCaptures_all_eo(gd.getNextBoolean());
		illustrationParameters.setCaptures_all(gd.getNextBoolean());
		
		boolean safety_OK =                            gd.getNextBoolean();
		if (!safety_OK) {
			IJ.showMessage("Warning","Operation aborted as safety checkbox is not checked");
			return false;
		}
		long startTime = System.nanoTime();
		//create directory lists
		
		double inter_gap = illustrationParameters.two_level_dirs ? illustrationParameters.inter_gap : 0.0;
		int seq_len =      illustrationParameters.two_level_dirs ? illustrationParameters.seq_len : 0;
		final CalibrationIllustration.CapturedScene [] captured_scenes = illustration.listCapturedScenes(
				illustrationParameters.destinationFolder,                  // String   result_path,
				illustrationParameters.sourceFolders, // String   captured_paths,
				// if both (seq_len <= 0) && (seq_gap <= 0) will create a single-level directories
				seq_len,                              // int      seq_len, // sequence fixed length (or 0)
				inter_gap, // 0.0,                    // double   seq_gap, // start a new sequence after gap longer than 
				illustrationParameters.getMin_ts(),// double   min_ts,
				illustrationParameters.getMax_ts(),// double   max_ts,
				illustrationParameters.isCaptures_all_lwir(),
				illustrationParameters.isCaptures_all_eo(),
				illustrationParameters.isCaptures_all());
		System.out.println("OrganizeSeries(): Directory tree tasks prepared in "
				+ IJ.d2s(0.000000001 * (System.nanoTime() - startTime), 3) + " sec");
		
		/*
		final CalibrationIllustration.CapturedScene [][] captured_series = new CalibrationIllustration.CapturedScene[illustrationParameters.num_source_folders][];
		for (int nfold = 0; nfold < illustrationParameters.num_source_folders; nfold++) if ((illustrationParameters.sourceFolders[nfold] != null) && !illustrationParameters.sourceFolders[nfold].equals("")) {
			System.out.println ("processing source folder "+nfold+": "+illustrationParameters.sourceFolders[nfold]);
			captured_series[nfold] = illustration.listCapturedScenes(
					illustrationParameters.sourceFolders[nfold], // String   captured_path,
					illustrationParameters.getMin_ts(),// double   min_ts,
					illustrationParameters.getMax_ts(),// double   max_ts,
					false, // illustrationParameters.isCaptures_all_lwir(),
					false, // illustrationParameters.isCaptures_all_eo(),
					false ); //  illustrationParameters.isCaptures_all());
		}
		*/
		// merge, then filter by all/not all
		System.out.println("Preparing to "+CalibrationIllustrationParameters.copy_modes[illustrationParameters.copy_mode]+" "+captured_scenes.length+" scenes");
		System.out.println();
		boolean relative_symlinks = true;
		for (CalibrationIllustration.CapturedScene cs:captured_scenes) {
			String scene_path = cs.getName();
			String [] srs_paths = cs.getImages();
			File dFile=new File(scene_path);
			if (!dFile.isDirectory() &&  !dFile.mkdirs()) {
				String msg="Failed to create directory "+scene_path;
				IJ.showMessage(msg);
				throw new IllegalArgumentException (msg);
			}
			for (int chn = 0; chn <srs_paths.length; chn++) if (srs_paths[chn]!= null){
				int basename_start = srs_paths[chn].lastIndexOf(Prefs.getFileSeparator());
				if (basename_start >= 0) {
					basename_start++;
				} else {
					basename_start = 0;
				}
				String basename = srs_paths[chn].substring(basename_start);
				Path link_name = Paths.get(scene_path+Prefs.getFileSeparator()+basename);
				Path link_target = Paths.get(srs_paths[chn]);
				if (CalibrationIllustrationParameters.copy_modes[illustrationParameters.copy_mode].equals("soft link")) {
					try {
						if (relative_symlinks) {
							Path targetRelative = link_name.getParent().relativize(link_target);	
							Files.createSymbolicLink(link_name, targetRelative);
						} else {
							Files.createSymbolicLink(link_name, link_target);
						}
					} catch (IOException x) {
						System.out.println("Failed to soft link "+link_target+" (target) to "+link_name+" (link name)");
						System.err.println(x);
					} catch (UnsupportedOperationException x) {
						System.out.println("Failed to soft link "+link_target+" (target) to "+link_name+" (link name)");
						// Some file systems do not support symbolic links.
						System.err.println(x);
					}
				} else if (CalibrationIllustrationParameters.copy_modes[illustrationParameters.copy_mode].equals("hard link")) {
					try {
					    Files.createLink(link_name, link_target);
					} catch (IOException x) {
						System.out.println("Failed to hardlink "+link_target+" (target) to "+link_name+" (link name)");
					    System.err.println(x);
					} catch (UnsupportedOperationException x) {
						System.out.println("Failed to hardlink "+link_target+" (target) to "+link_name+" (link name)");

					    // Some file systems do not
					    // support adding an existing
					    // file to a directory.
					    System.err.println(x);
					}
				} else if (CalibrationIllustrationParameters.copy_modes[illustrationParameters.copy_mode].equals("copy")) {
					try {
						Files.copy(link_target, link_name, StandardCopyOption.REPLACE_EXISTING);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						System.out.println("Failed to copy "+link_target+" to "+link_name);
						e.printStackTrace();
					}
				} else if (CalibrationIllustrationParameters.copy_modes[illustrationParameters.copy_mode].equals("copy")) {
					try {
						Files.move(link_target, link_name, StandardCopyOption.REPLACE_EXISTING);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						System.out.println("Failed to move "+link_target+" to "+link_name);
						e.printStackTrace();
					}
				}
			}
		}
		
		System.out.println("OrganizeSeries(): All done in "
				+ IJ.d2s(0.000000001 * (System.nanoTime() - startTime), 3) + " sec");
		
		return true;
	}
	

}
