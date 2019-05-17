package com.elphel.imagej.tensorflow;
/**
 * Copyright (C) 2018 Elphel, Inc.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;

import org.apache.ant.compress.taskdefs.Unzip;
import org.tensorflow.SavedModelBundle;
import org.tensorflow.Tensor;

import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.tileprocessor.ImageDtt;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;

/**

Note 0: articles & examples:

- https://divis.io/2017/11/enterprise-tensorflow-1/
- https://divis.io/2018/01/enterprise-tensorflow-code-examples/
- https://divis.io/2018/01/enterprise-tensorflow-2-saving-a-trained-model/
- https://divis.io/2018/01/enterprise-tensorflow-3-loading-a-savedmodel-in-java/
- https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/
- https://www.programcreek.com/java-api-examples/?api=org.tensorflow.SavedModelBundle
- https://github.com/imagej/imagej-tensorflow
- simple: http://www.riptutorial.com/tensorflow/example/32154/load-and-use-the-model-in-java-

Note 1: How to feed:
 a. https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/
 b. https://github.com/tensorflow/tensorflow/blob/master/tensorflow/java/src/main/java/org/tensorflow/examples/LabelImage.java



Note 2: https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/

...
Two types of objects will need manual closing for proper resource handling:
Sessions and Tensors
...



*/

public class TensorflowInferModel
{
	
    final static String TRAINED_MODEL_URL = "https://community.elphel.com/files/quad-stereo/ml/trained_model_v1.0.zip";
	
    final static String TRAINED_MODEL = "trained_model"; // /home/oleg/GIT/python3-imagej-tiff/data_sets/tf_data_5x5_main_13_heur/exportdir";
    final static String SERVING = "serve";
    
    final int         tilesX, tilesY, num_tiles, num_layers;
    final int         corr_side, corr_side2;
//    final long []     shape_corr2d;
//    final long []     shape_target_disparity;
//    final long []     shape_tiles_stage1;
//    final long []     shape_tiles_stage2;
    final FloatBuffer fb_corr2d;
    final FloatBuffer fb_target_disparity;
    final IntBuffer   fb_tiles_stage1;
    final IntBuffer   fb_tiles_stage2;
    final FloatBuffer fb_predicted;

/*
long[] shape = new long[] {batch, imageSize};
 */
//    public
    final SavedModelBundle bundle;

    // utils: download url
    private static Path download(String sourceURL, String targetDirectory) throws IOException
    {
        URL url = new URL(sourceURL);
        String fileName = sourceURL.substring(sourceURL.lastIndexOf('/') + 1, sourceURL.length());
        Path targetPath = new File(targetDirectory + File.separator + fileName).toPath();
        Files.copy(url.openStream(), targetPath, StandardCopyOption.REPLACE_EXISTING);

        return targetPath;
    }
    
    // utils: unpack zip to dir
    private static boolean download_and_unpack(String sourceURL, String targetDirectory) {
    	
    	Path zipped_model = null;
    	
    	System.out.println("Downloading "+sourceURL+". Please, wait...");
    	try {
			zipped_model = download(sourceURL, targetDirectory);
    	} catch (IOException e){
    		e.printStackTrace();
    		return false;
    	}
    	
		Unzip unzipper = new Unzip();
		unzipper.setSrc(zipped_model.toFile());
		unzipper.setDest(new File(targetDirectory));
		unzipper.execute();
		
		System.out.println(unzipper.getLocation());
    	
    	return true;
    }
    
    public TensorflowInferModel(int tilesX, int tilesY, int corr_side, int num_layers)
    {
    	this.tilesX = tilesX;
    	this.tilesY = tilesY;
    	this.num_tiles = tilesX*tilesY;
    	this.num_layers = num_layers;
    	this.corr_side = corr_side;
    	this.corr_side2 = corr_side * corr_side;
    	// allocate buffers to be used for tensors
    	this.fb_corr2d =           FloatBuffer.allocate(num_tiles * corr_side2*num_layers);
    	this.fb_target_disparity = FloatBuffer.allocate(num_tiles );
    	this.fb_tiles_stage1 =     IntBuffer.allocate(num_tiles );
    	this.fb_tiles_stage2 =     IntBuffer.allocate(num_tiles );
    	this.fb_predicted =        FloatBuffer.allocate(num_tiles );

    	//String resourceDir = System.getProperty("user.dir")+"/src/main/resources";
    	// ./target/classes/
    	String resourceDir = getClass().getClassLoader().getResource("").getFile();
    	String abs_model_path = null;
    	
    	try {
    		abs_model_path = getClass().getClassLoader().getResource(TRAINED_MODEL).getFile();
    	} catch (java.lang.NullPointerException e) {
    		//e.printStackTrace();
    		download_and_unpack(TRAINED_MODEL_URL, resourceDir);
    		// re-read
    		abs_model_path = getClass().getClassLoader().getResource(TRAINED_MODEL).getFile();
    		System.out.println("New downloaded path: "+abs_model_path);
    	}
    	
        System.out.println("TensorflowInferModel model path: "+abs_model_path);
        
        // this will load graph/data and open a session that does not need to be closed until the program is closed
////        bundle = null;
        bundle = SavedModelBundle.load(abs_model_path, SERVING);

//    	Operation opr = bundle.graph().operation("rv_stage1_out");
    }


	public void run_stage1(
			 FloatBuffer fb_corr2d,
			 FloatBuffer fb_target_disparity,
			 IntBuffer   fb_tiles_stage1)
	{
	    int ntiles = fb_tiles_stage1.limit(); // actual number of entries
	    long []     shape_corr2d =           new long [] {ntiles, corr_side2*num_layers};
	    long []     shape_target_disparity = new long [] {ntiles, 1};
	    long []     shape_tiles_stage1 =     new long [] {ntiles};

		final Tensor<Float>   t_corr2d =         Tensor.create(shape_corr2d,           fb_corr2d);
		final Tensor<Float>   target_disparity = Tensor.create(shape_target_disparity, fb_target_disparity);
		final Tensor<Integer> t_tiles_stage1 =   Tensor.create(shape_tiles_stage1,fb_tiles_stage1);

		// this run does not output any data, but maybe it should still be captured and disposed of?
		final Tensor<?> t_result_stage1 = bundle.session().runner()
        .feed("ph_corr2d",           t_corr2d)
        .feed("ph_target_disparity", target_disparity)
        .feed("ph_ntile",            t_tiles_stage1)
        .fetch("Disparity_net/stage1done:0")
        .run()
        .get(0);

		t_result_stage1.  close();
		t_corr2d.         close();
    	target_disparity. close();
    	t_tiles_stage1.   close();
	}

	public void run_stage2(
			 IntBuffer   fb_tiles_stage2,
			 FloatBuffer fb_disparity) {
	    int ntiles = fb_tiles_stage2.limit(); // actual number of entries
	    long []     shape_tiles_stage2 =         new long [] {ntiles};
		final Tensor<Integer> t_tiles_stage2 =   Tensor.create(shape_tiles_stage2, fb_tiles_stage2);
    	final Tensor<?> t_result_disparity = bundle.session().runner()
                .feed("ph_ntile_out",t_tiles_stage2)
                .fetch("Disparity_net/stage2_out_sparse:0")
                .run()
                .get(0);
    	fb_disparity.rewind();
    	t_result_disparity.writeTo(fb_disparity);

    	t_result_disparity. close();
    	t_tiles_stage2.     close();
	}


    //**************************
    // helper class to prepare TileProcessor task
	class TfInTile{
		float [] corr2d;
		float target_disparity;
		float gt_disparity;
		float gt_strength;
		int   tile;
		public void setFromCorr2d(
				float [][] data,
				int  tileX,
				int  tileY,
				int  corr_side,
				int  tilesX
				) {
			int corr_side2 = corr_side*corr_side;
			int layers = data.length -1;
			this.corr2d = new float [layers * corr_side2];
			float [] other = new float [corr_side2];

			int width = tilesX * corr_side;
			int tl = tileY * corr_side * width + tileX * corr_side; // index to the
			for (int nl = 0; nl <= layers; nl++) {
				if (nl < layers) {
					for (int row = 0; row < corr_side; row++) {
						System.arraycopy(data[nl], tl+width*row, corr2d, nl*corr_side2 + row*corr_side, corr_side);
					}
				} else {
					for (int row = 0; row < corr_side; row++) {
						System.arraycopy(data[layers], tl+width*row, other, row*corr_side, corr_side);
					}
					this.target_disparity = other[ImageDtt.ML_OTHER_TARGET];
					this.gt_disparity =     other[ImageDtt.ML_OTHER_GTRUTH];
					this.gt_strength =      other[ImageDtt.ML_OTHER_GTRUTH_STRENGTH];
				}
			}
			if (Float.isNaN(this.target_disparity)) {
				corr2d = new float[corr2d.length]; // zero them all
			}
			this.tile = tileY*tilesX + tileX;

		}
//System.arraycopy(sym_conv, 0, tile_in, 0, n2*n2)

	}
	public int test_tensorflow(
			boolean keep_empty) {

		int dbgX = 162;
		int dbgY = 121;
//		int dbgT = tilesX ^ dbgY + dbgX;

		int corr_side = 9;
		String [] slices = {"hor-pairs","vert-pairs","diagm-pair","vert-pairs","other"};
    	ImagePlus imp_src = WindowManager.getCurrentImage();
    	if (imp_src==null){
    		IJ.showMessage("Error","2D Correlation image stack required");
    		return -1;
    	}
    	ImageStack corr_stack = imp_src.getStack();
    	String [] labels =  corr_stack.getSliceLabels();
//    	for (int ii = 0; ii < labels.length; ii++) {
//    		System.out.println(ii+": "+labels[ii]);
//    	}
    	int tilesX= corr_stack.getWidth()  /  corr_side;
    	int tilesY= corr_stack.getHeight() /  corr_side;

    	float [][] corr_data = new float [slices.length][];
    	for (int nslice = 0; nslice < slices.length; nslice++ ) {
    		int ns = -1;
    		for (int i = 0; i < labels.length; i++) {
    			if (slices[nslice].equals(labels[i])) {
    				ns = i;
    				break;
    			}
    		}
    		if (ns < 0) {
    			System.out.println("Slice "+slices[nslice]+" is not found in the image");
    			return -1;
    		} else {
    			corr_data[nslice] = (float[]) corr_stack.getPixels(ns + 1) ;
    		}
    	}
    	ArrayList<TfInTile> tf_tiles = new ArrayList<TfInTile>();
    	for (int tileY = 0; tileY < tilesY; tileY++) {
    		for (int tileX = 0; tileX < tilesX; tileX++) {
    			if ((tileY==dbgY) && (tileX==dbgX)) {
    				System.out.println("tileY = "+tileY+", tileX = "+tileX+" tile = "+(tileY*tilesX + tileX));
    			}
    			TfInTile tf_tile = new TfInTile();
    			tf_tile.setFromCorr2d(
    					corr_data, // float [][] data,
    					tileX,     // int  tileX,
    					tileY,     // int  tileY,
    					corr_side, // int  corr_side,
    					tilesX);   // int  tilesX
    			if (keep_empty || !Float.isNaN(tf_tile.target_disparity)) {
//    				if (Float.isNaN(tf_tile.target_disparity)) {
//    					tf_tile.target_disparity = 0.0f;
//    					tf_tile.gt_strength =      0.0f;
//    				}
    				tf_tiles.add(tf_tile);
    			}
    		}
    	}

    	// sets the limit to the capacity and the position to zero.
    	fb_corr2d.clear();
    	fb_target_disparity.clear();
    	fb_tiles_stage1.clear();
    	fb_tiles_stage2.clear();
    	for (int i = 0; i < tf_tiles.size(); i++) {
    		TfInTile tf_tile = tf_tiles.get(i);
    		fb_corr2d.put          (tf_tile.corr2d);
    		float td = tf_tile.target_disparity;
    		if (Float.isNaN(td)) td = 0.0f;
    		fb_target_disparity.put(i, td);
    		fb_tiles_stage1.put(tf_tile.tile);
    		fb_tiles_stage2.put(tf_tile.tile);
    	}
//    	fb_target_disparity.limit(tf_tiles.size()); // put float absolute does not movew the pointer
    	fb_target_disparity.position(tf_tiles.size()); // put float absolute does not movew the pointer
    	fb_tiles_stage1.    position(tf_tiles.size()); // put float absolute does not movew the pointer
    	fb_tiles_stage2.    position(tf_tiles.size()); // put float absolute does not movew the pointer
    	//sets the limit to the current position and then sets the position to zero.
    	fb_corr2d.          flip();
    	fb_target_disparity.flip();
    	fb_tiles_stage1.    flip();
    	fb_tiles_stage2.    flip();
    	String [] titles = {"predicted", "target", "gt_disparity", "gt_strength", "nn_out", "nn_error","abs_err","abs_heur","clean_nn"};
    	fb_predicted.rewind(); // not needed for absolute get();
    	float [][] result = new float [titles.length][tilesX*tilesY];
/*
    	for (int i = 0; i < tf_tiles.size(); i++) {
    		TfInTile tf_tile = tf_tiles.get(i);
//    		result[0][i] = tf_tile.target_disparity + fb_predicted.get(i);
    		result[1][i] = tf_tile.target_disparity;
    		result[2][i] = tf_tile.gt_disparity;
    		result[3][i] = tf_tile.gt_strength;
    	}

		(new showDoubleFloatArrays()).showArrays(
				result,
				tilesX,
				tilesY,
				true,
				"NN_disparity-pre",
				titles);
*/


//if (corr_side > 0) return 0; // always
    	run_stage1(
    			fb_corr2d,           // FloatBuffer fb_corr2d,
    			fb_target_disparity, // FloatBuffer fb_target_disparity,
    			fb_tiles_stage1);    // IntBuffer   fb_tiles_stage1);
    	run_stage2(
    			fb_tiles_stage2, // IntBuffer   fb_tiles_stage2,
    			fb_predicted); // FloatBuffer fb_disparity)


    	if (!keep_empty) {
    		for (int i =0; i < result[0].length; i++) {
    			result[0][i]=Float.NaN;
    			result[1][i]=Float.NaN;
    			result[2][i]=Float.NaN;
    		}
    	}

    	fb_predicted.rewind(); // not needed for absolute get();
    	for (int i = 0; i < tf_tiles.size(); i++) {
    		TfInTile tf_tile = tf_tiles.get(i);
    		result[0][i] = tf_tile.target_disparity + fb_predicted.get(i);
    		result[1][i] = tf_tile.target_disparity;
    		result[2][i] = tf_tile.gt_disparity;
    		result[3][i] = tf_tile.gt_strength;
    		result[4][i] = fb_predicted.get(i);
    		result[5][i] = (result[3][i] > 0)?(result[0][i]-result[2][i]):Float.NaN;
    		result[6][i] = (result[3][i] > 0)?Math.abs(result[5][i]):Float.NaN;
    		result[7][i] = (result[3][i] > 0)?Math.abs(result[1][i] - result[2][i]):Float.NaN;
    		result[8][i] = (result[3][i] > 0)?result[0][i]:Float.NaN;

    	}
		(new ShowDoubleFloatArrays()).showArrays(
				result,
				tilesX,
				tilesY,
				true,
				"NN_disparity",
				titles);

    	return 0;
	}


}


















