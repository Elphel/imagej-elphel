package com.elphel.imagej.tileprocessor;
/**
 **
 ** QuadCLT - Process images with CLT-based methods (code specific to ImageJ plugin)
 ** Using CPU+GPU
 **
 ** Copyright (C) 2017-2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  QuadCLT.java is free software: you can redistribute it and/or modify
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

import java.awt.Rectangle;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.util.Properties;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.gpu.GPUTileProcessor;

import ij.ImagePlus;
import ij.ImageStack;

public class QuadCLT extends QuadCLTCPU {
	GPUTileProcessor.GpuQuad gpuQuad =                     null;	
	public QuadCLT(
			String                                          prefix,
			Properties                                      properties,
			EyesisCorrections                               eyesisCorrections,
			EyesisCorrectionParameters.CorrectionParameters correctionsParameters
			){
		super (prefix,
				properties,
				eyesisCorrections,
				correctionsParameters
				);
	}
	public void setGPU(GPUTileProcessor.GpuQuad gpuQuad) {
		this.gpuQuad = gpuQuad;
	}
	public GPUTileProcessor.GpuQuad getGPU() {
		return this.gpuQuad;
	}
	
	  public CLTPass3d CLTBackgroundMeasGPU( // measure background
			  CLTParameters       clt_parameters,
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  CLTPass3d scan_rslt = new CLTPass3d(tp);
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setPairMask(d,0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][]     tile_op =         tp.setSameTileOp(clt_parameters,  d, debugLevel);
		  double [][]  disparity_array = tp.setSameDisparity(0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		  // undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		  double [][][][]     clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][tilesY][tilesX][]; // will only be used inside?

		  double min_corr_selected = clt_parameters.min_corr; // 0.02

		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		  double [][] shiftXY = new double [4][2];
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }

		  double [][][][] texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()));
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  image_dtt.clt_aberrations_quad_corr(
				  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                      // per-tile operation bit codes
				  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
				  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
				  // correlation results - final and partial
				  clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null, // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				  null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  //	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				  disparity_map,    // [12][tp.tilesY * tp.tilesX]
				  texture_tiles,        // [tp.tilesY][tp.tilesX]["RGBA".length()][];
				  tilesX * image_dtt.transform_size, // imp_quad[0].getWidth(),       // final int width,
				  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				  clt_parameters.corr_sym,
				  clt_parameters.corr_offset,
				  clt_parameters.corr_red,
				  clt_parameters.corr_blue,
				  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
				  clt_parameters.corr_normalize, // normalize correlation results by rms
				  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
				  clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				  clt_parameters.max_corr_radius,
				  clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				  clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				  clt_parameters.kernel_step,
///				  image_dtt.transform_size,
				  clt_parameters.clt_window,
				  shiftXY, //
				  disparity_corr, // final double              disparity_corr, // disparity at infinity
				  (clt_parameters.fcorr_ignore? null: this.fine_corr),
				  clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85
				  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				  clt_parameters.tileX,         // final int               debug_tileX,
				  clt_parameters.tileY,         // final int               debug_tileY,
				  (clt_parameters.dbg_mode & 64) != 0, // no fract shift
				  (clt_parameters.dbg_mode & 128) != 0, // no convolve
				  //				  (clt_parameters.dbg_mode & 256) != 0, // transpose convolve
				  threadsMax,
				  debugLevel);

		  scan_rslt.disparity = disparity_array;
		  scan_rslt.tile_op = tile_op;
		  scan_rslt.disparity_map = disparity_map;
		  scan_rslt.texture_tiles = texture_tiles;
		  scan_rslt.is_measured =   true;
		  scan_rslt.is_combo =      false;
		  scan_rslt.resetProcessed();
		  return scan_rslt;
	  }
	
	
	
	public void processCLTQuadCorrGPU(
			ImagePlus []                                    imp_quad,
			boolean [][]                                    saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
			CLTParameters                                   clt_parameters,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			CorrectionColorProc.ColorGainsParameters        channelGainParameters,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			double []	                                    scaleExposures, // probably not needed here - restores brightness of the final image
			final int                                       threadsMax,  // maximal number of threads to launch
			final boolean                                   updateStatus,
			final int                                       debugLevel){
		if (gpuQuad == null) {
			System.out.println("GPU instance is not initialized, using CPU mode");
			processCLTQuadCorrCPU(
					imp_quad,              // ImagePlus []                  imp_quad, // should have properties "name"(base for saving results), "channel","path"
					saturation_imp,        // boolean [][]                  saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
					clt_parameters,        // CLTParameters                 clt_parameters,
					debayerParameters,     // EyesisCorrectionParameters.DebayerParameters     debayerParameters,
					colorProcParameters,   // ColorProcParameters                              colorProcParameters,
					channelGainParameters, // CorrectionColorProc.ColorGainsParameters         channelGainParameters,
					rgbParameters, // EyesisCorrectionParameters.RGBParameters             rgbParameters,
					scaleExposures, // double []	       scaleExposures, // probably not needed here
					false, // final boolean    apply_corr, // calculate and apply additional fine geometry correction
					false, // final boolean    infinity_corr, // calculate and apply geometry correction at infinity
					threadsMax,  // final int        threadsMax,  // maximal number of threads to launch
					updateStatus, // final boolean    updateStatus,
					debugLevel); // final int        debugLevel)			
			
			return;
		}
		
// GPU-specific
		boolean is_mono = isMonochrome();
		boolean is_lwir = isLwir();
		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		
/*
		double    fat_zero = clt_parameters.getGpuFatZero(is_mono); //   30.0;
		double [] scales = (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.gpu_weight_r, // 0.25
				clt_parameters.gpu_weight_b, // 0.25
				1.0 - clt_parameters.gpu_weight_r - clt_parameters.gpu_weight_b}); // 0.5
		double cwgreen = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
		double [] col_weights= (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.corr_red *  cwgreen,
				clt_parameters.corr_blue * cwgreen,
				cwgreen});
*/
		ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  is_mono,
				  is_lwir,
				  clt_parameters.getScaleStrength(isAux())); // 1.0);
		
		float [][] lpf_rgb = new float[][] {
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_r),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_b),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_g),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				!batch_mode);

		float [] lpf_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrSigma(is_mono));

		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				!batch_mode);

		float [] lpf_rb_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrRBSigma(is_mono));
		gpuQuad.setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				!batch_mode);
		final boolean use_aux = false; // currently GPU is configured for a single quad camera
		
//		String image_name = image_name; // correctionsParameters.getModelName((String) imp_quad_main[0].getProperty("name"));
//		String image_path = image_path; //  (String) imp_quad_main[0].getProperty("path"); // Only for debug output
		// now set to only 4 !
		if (debugLevel>1) System.out.println("processing: "+image_path);

		
// end of GPU-specific		

		boolean advanced=  correctionsParameters.zcorrect || correctionsParameters.equirectangular;
		boolean toRGB=     advanced? true: correctionsParameters.toRGB;
		
		ImagePlus [] results = new ImagePlus[imp_quad.length];
		for (int i = 0; i < results.length; i++) {
			results[i] = imp_quad[i];
			results[i].setTitle(results[i].getTitle()+"RAW");
		}
		if (debugLevel>1) System.out.println("processing: "+gpuQuad.quadCLT);

		double z_correction =  clt_parameters.z_correction;
		if (clt_parameters.z_corr_map.containsKey(image_name)){
			z_correction +=clt_parameters.z_corr_map.get(image_name);// not used in lwir
		}
		final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		// the following will be replaced by setFullFrameImages() with target disparities
		GPUTileProcessor.TpTask [] tp_tasks  = gpuQuad.setFullFrameImages( // when disparities array is known - use different arguments of setFullFrameImages
				null,                                 // Rectangle                 woi,
				clt_parameters.gpu_woi_round,         // boolean                   round_woi,
	    		(float) clt_parameters.disparity,     // float                     target_disparity, // apply same disparity to all tiles
	    		(float) disparity_corr,
	    		0xf,                                  // int                       out_image, // from which tiles to generate image (currently 0/1)
	    		0x3f,                                 // int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
	    		debugLevel);                          // final int                 debugLevel) - not yet used

		gpuQuad.setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_aux); // boolean use_aux)

		gpuQuad.execSetTilesOffsets();
		gpuQuad.execConvertDirect();
		gpuQuad.execImcltRbgAll(is_mono); 

		// get data back from GPU
		float [][][] iclt_fimg = new float [gpuQuad.getNumCams()][][];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			iclt_fimg[ncam] = gpuQuad.getRBG(ncam);
		}

		int out_width =  gpuQuad.getImageWidth()  + gpuQuad.getDttSize();
		int out_height = gpuQuad.getImageHeight() + gpuQuad.getDttSize();
//		int tilesX =     gpuQuad.getImageWidth()  / gpuQuad.getDttSize();
//		int tilesY =     gpuQuad.getImageHeight() / gpuQuad.getDttSize();
		
		/* Prepare 4-channel images*/
		ImagePlus [] imps_RGB = new ImagePlus[iclt_fimg.length];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
//			String title=image_name+"-"+String.format("%02d", ncam);
            String title=String.format("%s%s-%02d",image_name, sAux(), ncam);
			
			imps_RGB[ncam] = linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					title, // String name,
					"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
					toRGB,
					!correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					!batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // boolean saveShowFinal,        // save/show result (color image?)
					iclt_fimg[ncam],
					out_width,
					out_height,
					1.0, // scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
					debugLevel );
		}
		
		if (clt_parameters.gen_chn_img) { // save and show 4-slice image
			// combine to a sliced color image
			// assuming total number of images to be multiple of 4
			//			  int [] slice_seq = {0,1,3,2}; //clockwise
			int [] slice_seq = new int[results.length];
			for (int i = 0; i < slice_seq.length; i++) {
				slice_seq[i] = i ^ ((i >> 1) & 1); // 0,1,3,2,4,5,7,6, ...
			}
			int width = imps_RGB[0].getWidth();
			int height = imps_RGB[0].getHeight();
			ImageStack array_stack=new ImageStack(width,height);
			for (int i = 0; i<slice_seq.length; i++){
				if (imps_RGB[slice_seq[i]] != null) {
					array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
				} else {
					array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
				}
			}
			ImagePlus imp_stack = new ImagePlus(image_name+sAux()+"-SHIFTED-D"+clt_parameters.disparity, array_stack);
			imp_stack.getProcessor().resetMinAndMax();
			if (!batch_mode) {
				imp_stack.updateAndDraw();
			}
			eyesisCorrections.saveAndShowEnable(
					imp_stack,  // ImagePlus             imp,
					correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
					true, // boolean               enableSave,
					!batch_mode) ;// boolean               enableShow);
		}

		if (clt_parameters.gen_4_img) { // save 4 JPEG images
			// Save as individual JPEG images in the model directory
			String x3d_path= correctionsParameters.selectX3dDirectory(
					image_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					correctionsParameters.x3dModelVersion,
					true,  // smart,
					true);  //newAllowed, // save
			for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
				eyesisCorrections.saveAndShow(
						imps_RGB[sub_img],
						x3d_path,
						correctionsParameters.png && !clt_parameters.black_back,
						!batch_mode && clt_parameters.show_textures,
						correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
			}
			String model_path= correctionsParameters.selectX3dDirectory(
					image_name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					null,
					true,  // smart,
					true);  //newAllowed, // save

			createThumbNailImage(
					imps_RGB[0],
					model_path,
					"thumb"+sAux(),
					debugLevel);

		}

/**
       if (colorProcParameters.isLwir() && colorProcParameters.lwir_autorange) {
            double rel_low =  colorProcParameters.lwir_low;
            double rel_high = colorProcParameters.lwir_high;
            if (!Double.isNaN(getLwirOffset())) {
                rel_low -=  getLwirOffset();
                rel_high -= getLwirOffset();
            }
            double [] cold_hot =  autorange(
                    iclt_data, // double [][][] iclt_data, //  [iQuad][ncol][i] - normally only [][2][] is non-null
                    rel_low, // double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
                    rel_high, // double hard_hot, // matches data, DC (this.lwir_offset)  subtracted
                    colorProcParameters.lwir_too_cold, // double too_cold, // pixels per image
                    colorProcParameters.lwir_too_hot, // double too_hot,  // pixels per image
                    1024); // int num_bins)
            if (cold_hot != null) {
                if (!Double.isNaN(getLwirOffset())) {
                    cold_hot[0] += getLwirOffset();
                    cold_hot[1] += getLwirOffset();
                }
            }
            setColdHot(cold_hot); // will be used for shifted images and for texture tiles
        }
		
 */
		  
		
	}
	
	
	
	
	
	public static ImagePlus [] processCLTQuadCorrPairGpu(
			QuadCLT                                         quadCLT_main, // use gpuQuad_main.quadCLT
			QuadCLT                                         quadCLT_aux, // use gpuQuad_aux.quadCLT
			ImagePlus []                                    imp_quad_main,
			ImagePlus []                                    imp_quad_aux,
			boolean [][]                                    saturation_main, // (near) saturated pixels or null
			boolean [][]                                    saturation_aux, // (near) saturated pixels or null
			CLTParameters                                   clt_parameters,
			EyesisCorrectionParameters.CorrectionParameters ecp,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			ColorProcParameters                             colorProcParameters_aux,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			double []	                                    scaleExposures_main, // probably not needed here - restores brightness of the final image
			double []	                                    scaleExposures_aux, // probably not needed here - restores brightness of the final image
			boolean                                         notch_mode, // not used here: use pole-detection mode for inter-camera correlation
			final int                                       lt_rad,     // not used here: low texture mode - inter-correlation is averaged between the neighbors before argmax-ing, using
			final int        threadsMax,  // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel){
		boolean resetGC = false;
		boolean resetEV = false;
		
		if (resetGC) {
			if (quadCLT_main.getGPU() != null) quadCLT_main.getGPU().resetGeometryCorrection();
			if (quadCLT_aux.getGPU() != null)  quadCLT_aux.getGPU().resetGeometryCorrection();
		}

		if (resetEV) {
			if (quadCLT_main.getGPU() != null) quadCLT_main.getGPU().resetGeometryCorrectionVector();
			if (quadCLT_aux.getGPU() != null)  quadCLT_aux.getGPU().resetGeometryCorrectionVector();
		}
		
// get fat_zero (absolute) and color scales
		boolean is_mono = quadCLT_main.isMonochrome();
		boolean is_lwir = quadCLT_main.isLwir();

		double    fat_zero = clt_parameters.getGpuFatZero(is_mono); //   30.0;
		double [] scales = (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.gpu_weight_r, // 0.25
				clt_parameters.gpu_weight_b, // 0.25
				1.0 - clt_parameters.gpu_weight_r - clt_parameters.gpu_weight_b}); // 0.5
		double cwgreen = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
		double [] col_weights= (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.corr_red *  cwgreen,
				clt_parameters.corr_blue * cwgreen,
				cwgreen});

		ImageDtt image_dtt = new ImageDtt(
				  clt_parameters.transform_size,
				  is_mono,
				  is_lwir,
				  1.0);
		float [][] lpf_rgb = new float[][] {
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_r),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_b),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_g),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_m)
		};
		quadCLT_main.getGPU().setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				debugLevel > -1);

		float [] lpf_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrSigma(is_mono));

		quadCLT_main.getGPU().setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				debugLevel > -1);

		float [] lpf_rb_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrRBSigma(is_mono));
		quadCLT_main.getGPU().setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				debugLevel > -1);

		final boolean use_aux = false; // currently GPU is configured for a single quad camera

		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		boolean toRGB=     quadCLT_main.correctionsParameters.toRGB;

		// may use this.StartTime to report intermediate steps execution times
		String name = quadCLT_main.image_name; // correctionsParameters.getModelName((String) imp_quad_main[0].getProperty("name"));
		String path = quadCLT_main.image_path; //  (String) imp_quad_main[0].getProperty("path"); // Only for debug output
		// now set to only 4 !
		ImagePlus [] results = new ImagePlus[imp_quad_main.length]; // + imp_quad_aux.length];
		for (int i = 0; i < results.length; i++) {
			if (i< imp_quad_main.length) {
				results[i] = imp_quad_main[i];
			} else {
				results[i] = imp_quad_aux[i-imp_quad_main.length];
			}
			results[i].setTitle(results[i].getTitle()+"RAW");
		}
		if (debugLevel>1) System.out.println("processing: "+path);

		// Set task clt_parameters.disparity
		Rectangle twoi = clt_parameters.gpu_woi? (new Rectangle ( // normally - full window (0,0,324,242)
				clt_parameters.gpu_woi_tx,
				clt_parameters.gpu_woi_ty,
				clt_parameters.gpu_woi_twidth,
				clt_parameters.gpu_woi_theight)): null;
		
		// the following will be replaced by setFullFrameImages() with target disparities
		
		double z_correction =  clt_parameters.z_correction;
		if (clt_parameters.z_corr_map.containsKey(quadCLT_main.image_name)){
			z_correction +=clt_parameters.z_corr_map.get(quadCLT_main.image_name);// not used in lwir
		}
		final double disparity_corr = (z_correction == 0) ? 0.0 : quadCLT_main.geometryCorrection.getDisparityFromZ(1.0/z_correction);
		
		
		GPUTileProcessor.TpTask [] tp_tasks  = quadCLT_main.getGPU().setFullFrameImages( // when disparities array is known - use different arguments of setFullFrameImages
				twoi,                                 // Rectangle                 woi,
				clt_parameters.gpu_woi_round,         // boolean                   round_woi,
	    		(float) clt_parameters.disparity,     // float                     target_disparity, // apply same disparity to all tiles
	    		(float) disparity_corr,               // float                     disparity_corr, // add to disparity (at infinity)
	    		0xf,                                  // int                       out_image, // from which tiles to generate image (currently 0/1)
	    		0x3f,                                 // int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
	    		debugLevel);                          // final int                 debugLevel) - not yet used

		quadCLT_main.getGPU().setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_aux); // boolean use_aux)
/*
		quadCLT_main.getGPU().setGeometryCorrection( // copy Geometry correction data to the GPU memory (once per camera group?), allocate GPU memory for the rotation/deriv. matrices
				quadCLT_main.getGeometryCorrection(),
				false); // boolean use_java_rByRDist) { // false - use newer GPU execCalcReverseDistortions); // once
		
		quadCLT_main.getGPU().setExtrinsicsVector(quadCLT_main.getGeometryCorrection().getCorrVector()); // for each new image - copy geometry correction vector to the GPU
		
		// Optionally save offsets here?
		if (clt_parameters.gpu_save_ports_xy) {
			savePortsXY(
					ecp,        // EyesisCorrectionParameters.CorrectionParameters ecp,
					tp_tasks);  // GPUTileProcessor.TpTask [] tp_tasks
		}
*/		
		// All set, run kernel (correct and convert)
		int NREPEAT = 1; // 00;
		System.out.println("\n------------ Running GPU "+NREPEAT+" times ----------------");
		long startGPU=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Calculate reverse radial distortion table from gpu_geometry_correction
			// Needed once during initialization of GeometryCorrection
/// Following will be executed when/if needed			
///			quadCLT_main.getGPU().execCalcReverseDistortions();       

		}
		long startRotDerivs=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Calculate rotation matrices and their derivatives
			// Needed after extrinsics changed
/// Following will be executed when/if needed			
///			quadCLT_main.getGPU().execRotDerivs();
		}

		long startTasksSetup=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Calculate tiles offsets (before each direct conversion run)
			quadCLT_main.getGPU().execSetTilesOffsets();
		}

		long startDirectConvert=System.nanoTime();

		for (int i = 0; i < NREPEAT; i++ ) {
			// Direct CLT conversion and aberration correction
			quadCLT_main.getGPU().execConvertDirect();
		}

		long startIMCLT=System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ ) {
			// Generate corrected image(s) from the CLT representation created by execConvertDirect()
			quadCLT_main.getGPU().execImcltRbgAll(quadCLT_main.isMonochrome());
		}
		long endImcltTime = System.nanoTime();

		long startCorr2d=System.nanoTime();   // System.nanoTime();
		for (int i = 0; i < NREPEAT; i++ )
			//Generate 2D phase correlations from the CLT representation
			quadCLT_main.getGPU().execCorr2D(
	    		scales,// double [] scales,
	    		fat_zero, // double fat_zero);
	    		clt_parameters.gpu_corr_rad); // int corr_radius

		long endCorr2d = System.nanoTime();
// run textures
		long startTextures = System.nanoTime();   // System.nanoTime();
		boolean   calc_textures = clt_parameters.gpu_show_jtextures; //  true;
		boolean   calc_extra =    clt_parameters.gpu_show_extra; //  true;
		for (int i = 0; i < NREPEAT; i++ )
			//Generate non-overlapping (16x16) texture tiles, prepare 
			quadCLT_main.getGPU().execTextures(
				col_weights,                   // double [] color_weights,
				quadCLT_main.isLwir(),         // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change - never used in GPU ?
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
				calc_textures,                 // boolean   calc_textures,
				calc_extra);                   // boolean   calc_extra)

		long endTextures = System.nanoTime();
// run texturesRBGA
		long startTexturesRBGA = System.nanoTime();   // System.nanoTime();

		for (int i = 0; i < NREPEAT; i++ )
// Generate combined (overlapping) texture
			quadCLT_main.getGPU().execRBGA(
				col_weights,                   // double [] color_weights,
				quadCLT_main.isLwir(),         // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove);   // boolean   dust_remove,
		long endTexturesRBGA = System.nanoTime();
		long endGPUTime = System.nanoTime();

		long calcReverseTime=      (startRotDerivs-     startGPU)           /NREPEAT;
		long rotDerivsTime=        (startTasksSetup-    startRotDerivs)     /NREPEAT;
		long tasksSetupTime=       (startDirectConvert- startTasksSetup)    /NREPEAT;
		long firstGPUTime=         (startIMCLT-         startDirectConvert) /NREPEAT;
		long runImcltTime =        (endImcltTime -      startIMCLT)         /NREPEAT;
		long runCorr2DTime =       (endCorr2d -         startCorr2d)        /NREPEAT;
		long runTexturesTime =     (endTextures -       startTextures)      /NREPEAT;
		long runTexturesRBGATime = (endTexturesRBGA -   startTexturesRBGA)  /NREPEAT;
		long runGPUTime =          (endGPUTime -        startGPU)           /NREPEAT;
		// run corr2d
//RotDerivs
		System.out.println("\n------------ End of running GPU "+NREPEAT+" times ----------------");
		System.out.println("GPU run time ="+        (runGPUTime * 1.0e-6)+"ms");
		System.out.println(" - calc reverse dist.: "+(calcReverseTime*1.0e-6)+"ms");
		System.out.println(" - rot/derivs:         "+(rotDerivsTime*1.0e-6)+"ms");
		System.out.println(" - tasks setup:        "+(tasksSetupTime*1.0e-6)+"ms");
		System.out.println(" - direct conversion:  "+(firstGPUTime*1.0e-6)+"ms");
		System.out.println(" - imclt:              "+(runImcltTime*1.0e-6)+"ms");
		System.out.println(" - corr2D:             "+(runCorr2DTime*1.0e-6)+"ms");
		System.out.println(" - textures:           "+(runTexturesTime*1.0e-6)+"ms");
		System.out.println(" - RGBA:               "+(runTexturesRBGATime*1.0e-6)+"ms");
		// get data back from GPU
		float [][][] iclt_fimg = new float [quadCLT_main.getGPU().getNumCams()][][];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			iclt_fimg[ncam] = quadCLT_main.getGPU().getRBG(ncam);
		}

		int out_width =  quadCLT_main.getGPU().getImageWidth()  + quadCLT_main.getGPU().getDttSize();
		int out_height = quadCLT_main.getGPU().getImageHeight() + quadCLT_main.getGPU().getDttSize();
		int tilesX =     quadCLT_main.getGPU().getImageWidth()  / quadCLT_main.getGPU().getDttSize();
		int tilesY =     quadCLT_main.getGPU().getImageHeight() / quadCLT_main.getGPU().getDttSize();
		
		
		// Read extra data for macro generation: 4 DIFFs, 4 of R,  4 of B, 4 of G
		// Available after gpu_diff_rgb_combo is generated in execTextures
		if (calc_extra) {
			String [] extra_group_titles = {"DIFF","Red","Blue","Green"};
			String [] extra_titles = new String [extra_group_titles.length*quadCLT_main.getGPU().getNumCams()];
			for (int g = 0; g < extra_group_titles.length;g++) {
				for (int ncam=0; ncam < quadCLT_main.getGPU().getNumCams();ncam++) {
					extra_titles[g * quadCLT_main.getGPU().getNumCams() + ncam]= extra_group_titles[g]+"-"+ncam;
				}
			}
			float [][] extra = quadCLT_main.getGPU().getExtra();
			(new ShowDoubleFloatArrays()).showArrays(
					extra,
					tilesX,
					tilesY,
					true,
					name+"-EXTRA-D"+clt_parameters.disparity,
					extra_titles);
		}
		/* Prepare 4-channel images*/
		ImagePlus [] imps_RGB = new ImagePlus[iclt_fimg.length];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			String title=name+"-"+String.format("%02d", ncam);
			imps_RGB[ncam] = quadCLT_main.linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					title, // String name,
					"-D"+clt_parameters.disparity, //String suffix, // such as disparity=...
					toRGB,
					!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					!batch_mode, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // boolean saveShowFinal,        // save/show result (color image?)
					iclt_fimg[ncam],
					out_width,
					out_height,
					1.0, // scaleExposures[iAux][iSubCam], // double scaleExposure, // is it needed?
					debugLevel );
		}
		
		
		//Show 2D correlations
		int [] wh = new int[2];
		if (clt_parameters.show_corr) {
			int [] corr_indices = quadCLT_main.getGPU().getCorrIndices();
			float [][] corr2D = quadCLT_main.getGPU().getCorr2D(
					clt_parameters.gpu_corr_rad); //  int corr_rad);
			// convert to 6-layer image		 using tasks
			double [][] dbg_corr = GPUTileProcessor.getCorr2DView(
					tilesX,
					tilesY,
					corr_indices,
					corr2D,
					wh);
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_corr,
					wh[0],
					wh[1],
					true,
					name+"-CORR2D-D"+clt_parameters.disparity,
					GPUTileProcessor.getCorrTitles());
		}
// convert to overlapping and show
		if (clt_parameters.gen_chn_img) { // save and show 4-slice image
			// combine to a sliced color image
			// assuming total number of images to be multiple of 4
			//			  int [] slice_seq = {0,1,3,2}; //clockwise
			int [] slice_seq = new int[results.length];
			for (int i = 0; i < slice_seq.length; i++) {
				slice_seq[i] = i ^ ((i >> 1) & 1); // 0,1,3,2,4,5,7,6, ...
			}
			int width = imps_RGB[0].getWidth();
			int height = imps_RGB[0].getHeight();
			ImageStack array_stack=new ImageStack(width,height);
			for (int i = 0; i<slice_seq.length; i++){
				if (imps_RGB[slice_seq[i]] != null) {
					array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
				} else {
					array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
				}
			}
			ImagePlus imp_stack = new ImagePlus(name+"-SHIFTED-D"+clt_parameters.disparity, array_stack);
			imp_stack.getProcessor().resetMinAndMax();
			if (!batch_mode) {
				imp_stack.updateAndDraw();
			}
			quadCLT_main.eyesisCorrections.saveAndShowEnable(
					imp_stack,  // ImagePlus             imp,
					quadCLT_main.correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
					true, // boolean               enableSave,
					!batch_mode) ;// boolean               enableShow);
		}

		if (clt_parameters.gen_4_img) { // save 4 JPEG images
			// Save as individual JPEG images in the model directory
			String x3d_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
					name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					quadCLT_main.correctionsParameters.x3dModelVersion,
					true,  // smart,
					true);  //newAllowed, // save
			for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
				quadCLT_main.eyesisCorrections.saveAndShow(
						imps_RGB[sub_img],
						x3d_path,
						quadCLT_main.correctionsParameters.png && !clt_parameters.black_back,
						!batch_mode && clt_parameters.show_textures,
						quadCLT_main.correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
						(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
			}
			String model_path= quadCLT_main.correctionsParameters.selectX3dDirectory(
					name, // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
					null,
					true,  // smart,
					true);  //newAllowed, // save

			quadCLT_main.createThumbNailImage(
					imps_RGB[0],
					model_path,
					"thumb",
					debugLevel);

		}
		// Use GPU prepared RBGA
		if (clt_parameters.show_rgba_color) {
			Rectangle woi = new Rectangle();
			float [][] rbga = quadCLT_main.getGPU().getRBGA(
					(is_mono?1:3), // int     num_colors,
					woi);
			(new ShowDoubleFloatArrays()).showArrays( // show slices RBGA (colors - 256, A - 1.0)
					rbga,
					woi.width,
					woi.height,
					true,
					name+"-RGBA-STACK-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R"),
					new String[] {"R","B","G","A"}
					);



			// for now - use just RGB. Later add option for RGBA
			float [][] rgb_main = {rbga[0],rbga[1],rbga[2]};
			float [][] rgba_main = {rbga[0],rbga[1],rbga[2],rbga[3]};
			ImagePlus imp_rgba_main = quadCLT_main.linearStackToColor(
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					name+"-texture", // String name,
					"-D"+clt_parameters.disparity+"-MAINGPU", //String suffix, // such as disparity=...
					toRGB,
					!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // true, // boolean saveShowFinal,        // save/show result (color image?)
					((clt_parameters.alpha1 > 0)? rgba_main: rgb_main),
					woi.width, // clt_parameters.gpu_woi_twidth * image_dtt.transform_size, //  tilesX *  image_dtt.transform_size,
					woi.height, // clt_parameters.gpu_woi_theight *image_dtt.transform_size, //  tilesY *  image_dtt.transform_size,
					1.0,         // double scaleExposure, // is it needed?
					debugLevel );

			int width = imp_rgba_main.getWidth();
			int height =imp_rgba_main.getHeight();
			ImageStack texture_stack=new ImageStack(width,height);
			texture_stack.addSlice("main",      imp_rgba_main.getProcessor().getPixels()); // single slice
			ImagePlus imp_texture_stack = new ImagePlus(
					name+"-RGBA-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R"),
					texture_stack);
			imp_texture_stack.getProcessor().resetMinAndMax();
			String results_path= quadCLT_main.correctionsParameters.selectResultsDirectory( // selectX3dDirectory(
					true,  // smart,
					true);  //newAllowed, // save
			quadCLT_main.eyesisCorrections.saveAndShow( // save and show color RGBA texture
					imp_texture_stack,
					results_path,
					true, // quadCLT_main.correctionsParameters.png && !clt_parameters.black_back,
					true, // !batch_mode && clt_parameters.show_textures,
					0, // quadCLT_main.correctionsParameters.JPEG_quality, // jpegQuality); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
					(debugLevel > 0) ? debugLevel : 1); // int debugLevel (print what it saves)
		}
		

		// convert textures to RGBA in Java
//		if (clt_parameters.show_rgba_color && (debugLevel > 100)) { // disabling
		if (clt_parameters.show_rgba_color && clt_parameters.gpu_show_jtextures) {  // seem to have wrong colors
			
			int numcol = quadCLT_main.isMonochrome()?1:3;
			int ports = imp_quad_main.length;
			int [] texture_indices = quadCLT_main.getGPU().getTextureIndices();
			int          num_src_slices = numcol + 1 + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
			float [] flat_textures =  quadCLT_main.getGPU().getFlatTextures( // fatal error has been detected by the Java Runtime Environment:
					texture_indices.length,
		    		(is_mono?1:3), // int     num_colors,
		    		clt_parameters.keep_weights); // boolean keep_weights);
	    	int texture_slice_size = (2 * quadCLT_main.getGPU().getDttSize())* (2 * quadCLT_main.getGPU().getDttSize());
	    	int texture_tile_size = texture_slice_size * num_src_slices ;

			if (debugLevel > -1) {
		    	for (int indx = 0; indx < texture_indices.length; indx++) if ((texture_indices[indx] & (1 << GPUTileProcessor.LIST_TEXTURE_BIT)) != 0){
		    		int tile = texture_indices[indx] >> GPUTileProcessor.CORR_NTILE_SHIFT;
		    		int tileX = tile % tilesX;
		    		int tileY = tile / tilesX;
		    		if ((tileY == clt_parameters.tileY) && (tileX == clt_parameters.tileX)) {

		    			System.out.println("=== tileX= "+tileX+" tileY= "+tileY+" tile="+tile+" ===");

		    			for (int slice =0; slice < num_src_slices; slice++) {
		    				System.out.println("=== Slice="+slice+" ===");
		    				for (int i = 0; i < 2 * quadCLT_main.getGPU().getDttSize(); i++) {
		    					for (int j = 0; j < 2 * quadCLT_main.getGPU().getDttSize(); j++) {
		    						System.out.print(String.format("%10.4f ",
		    								flat_textures[indx*texture_tile_size + slice* texture_slice_size + 2 * quadCLT_main.getGPU().getDttSize() * i + j]));
		    					}
		    					System.out.println();
		    				}
		    			}
		    		}
		    	}
			}
			double [][][][] texture_tiles =     quadCLT_main.getGPU().doubleTextures(
		    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
		    		texture_indices,                  // int []       indices,
		    		flat_textures,                    // float [][][] ftextures,
		    		tilesX,                           // int          full_width,
		    		4, // rbga only /int          num_slices
		    		num_src_slices // int          num_src_slices
		    		);

			if ((debugLevel > -1) && (clt_parameters.tileX >= 0) && (clt_parameters.tileY >= 0) && (clt_parameters.tileX < tilesX) && (clt_parameters.tileY < tilesY)) {
				String [] rgba_titles = {"red","blue","green","alpha"};
				String [] rgba_weights_titles = {"red","blue","green","alpha","port0","port1","port2","port3","r-rms","b-rms","g-rms","w-rms"};
				double [][] texture_tile = texture_tiles[clt_parameters.tileY][clt_parameters.tileX];
				int tile = +clt_parameters.tileY * tilesX  +clt_parameters.tileX;
    			System.out.println("=== tileX= "+clt_parameters.tileX+" tileY= "+clt_parameters.tileY+" tile="+tile+" ===");

    			for (int slice =0; slice < texture_tile.length; slice++) {
    				System.out.println("\n=== Slice="+slice+" ===");
    				for (int i = 0; i < 2 * quadCLT_main.getGPU().getDttSize(); i++) {
    					for (int j = 0; j < 2 * quadCLT_main.getGPU().getDttSize(); j++) {
    						System.out.print(String.format("%10.4f ",
    								texture_tile[slice][2 * quadCLT_main.getGPU().getDttSize() * i + j]));
    					}
    					System.out.println();
    				}
    			}
    			(new ShowDoubleFloatArrays()).showArrays(
						texture_tile,
						2 * image_dtt.transform_size,
						2 * image_dtt.transform_size,
						true,
						name + "-TXTNOL-GPU-D"+clt_parameters.disparity+"-X"+clt_parameters.tileX+"-Y"+clt_parameters.tileY,
						(clt_parameters.keep_weights?rgba_weights_titles:rgba_titles));

			}


			int alpha_index = 3;
			// in monochrome mode only MONO_CHN == GREEN_CHN is used, R and B are null
			double [][] texture_overlap_main = image_dtt.combineRBGATiles(
					texture_tiles,                 // array [tp.tilesY][tp.tilesX][4][4*transform_size] or [tp.tilesY][tp.tilesX]{null}
//					image_dtt.transform_size,
					true,                         // when false - output each tile as 16x16, true - overlap to make 8x8
					clt_parameters.sharp_alpha,    // combining mode for alpha channel: false - treat as RGB, true - apply center 8x8 only
					threadsMax,                    // maximal number of threads to launch
					debugLevel);
			if (clt_parameters.alpha1 > 0){ // negative or 0 - keep alpha as it was
				double scale = (clt_parameters.alpha1 > clt_parameters.alpha0) ? (1.0/(clt_parameters.alpha1 - clt_parameters.alpha0)) : 0.0;
				for (int i = 0; i < texture_overlap_main[alpha_index].length; i++){
					double d = texture_overlap_main[alpha_index][i];
					if      (d >=clt_parameters.alpha1) d = 1.0;
					else if (d <=clt_parameters.alpha0) d = 0.0;
					else d = scale * (d- clt_parameters.alpha0);
					texture_overlap_main[alpha_index][i] = d;
				}
			}

			// for now - use just RGB. Later add option for RGBA
			double [][] texture_rgb_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2]};
			double [][] texture_rgba_main = {texture_overlap_main[0],texture_overlap_main[1],texture_overlap_main[2],texture_overlap_main[3]};
			ImagePlus imp_texture_main = quadCLT_main.linearStackToColor(
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					name+"-texture", // String name,
					"-D"+clt_parameters.disparity+"-MAINGPU", //String suffix, // such as disparity=...
					toRGB,
					!quadCLT_main.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // true, // boolean saveShowFinal,        // save/show result (color image?)
					((clt_parameters.alpha1 > 0)? texture_rgba_main: texture_rgb_main),
					tilesX *  image_dtt.transform_size,
					tilesY *  image_dtt.transform_size,
					1.0,         // double scaleExposure, // is it needed?
					debugLevel );

			int width = imp_texture_main.getWidth();
			int height =imp_texture_main.getHeight();
			ImageStack texture_stack=new ImageStack(width,height);
			texture_stack.addSlice("main",      imp_texture_main.getProcessor().getPixels()); // single slice
			ImagePlus imp_texture_stack = new ImagePlus(name+"-TEXTURES-D"+clt_parameters.disparity, texture_stack);
			imp_texture_stack.getProcessor().resetMinAndMax();
			imp_texture_stack.show();
		}

		return results;
	}
	
	
	public static boolean savePortsXY(
			EyesisCorrectionParameters.CorrectionParameters ecp,
			GPUTileProcessor.TpTask [] tp_tasks
			) {
		if ((ecp.tile_processor_gpu != null) && !ecp.tile_processor_gpu.isEmpty()) {
			int quad = 4;
			String file_prefix = ecp.tile_processor_gpu +"clt/main";
			for (int chn = 0; chn < quad; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".portsxy";
				FileOutputStream fos;
				try {
					fos = new FileOutputStream(img_path);
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					System.out.println("Could not write to "+img_path+" (file not found) port offsets");
					break;
				}
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				ByteBuffer bb = ByteBuffer.allocate(tp_tasks.length * 2 * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int i = 0; i <  tp_tasks.length; i++) {
					bb.putFloat((tp_tasks[i].xy[chn][0])); // x-offset
					bb.putFloat((tp_tasks[i].xy[chn][1])); // y-offset
				}
				bb.flip();
				try {
					channel.write(bb);
				} catch (IOException e) {
					System.out.println("Could not write to "+img_path+" port offsets");
					break;
				}
				try {
					dos.close();
				} catch (IOException e) {
					System.out.println("Could not close DataOutputStream for "+img_path+" port offsets");
				}
				System.out.println("Wrote port offsets to "+img_path+".");
			}
			return true;
		} else {
			return false;
		}
		
	}
	
}
