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

//import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;

import java.awt.Rectangle;
import java.io.DataOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;
import java.util.Arrays;
import java.util.Properties;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAccumulator;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.cameras.ColorProcParameters;
import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.correction.CorrectionColorProc;
import com.elphel.imagej.correction.EyesisCorrections;
import com.elphel.imagej.gpu.GPUTileProcessor;
import com.elphel.imagej.gpu.GpuQuad;
import com.elphel.imagej.gpu.TpTask;
import com.elphel.imagej.tileprocessor.QuadCLTCPU.SetChannels;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.io.FileSaver;
public class QuadCLT extends QuadCLTCPU {
	int dbg_lev = 1;
	private GpuQuad gpuQuad =              null; 	// use updateQuadCLT() to update after switching to a different
	// QuadCLT instance, this.getGPU().getQuadCLT() should be equal to this
	/*
			if (getGPU().getQuadCLT() != this) {
				getGPU().updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
			}

	 */
	
	public QuadCLT( // creates AUX, does not know numSensors?
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
	
	public GpuQuad getGPUQuad() {
		return gpuQuad;
	}
	
	public boolean hasGPU() {
		return (gpuQuad != null);
	}
	
	
	QuadCLT saveQuadClt() {
		if (gpuQuad == null) {
			return null;
		}
		QuadCLT savedQuadClt =  gpuQuad.getQuadCLT();
		if (savedQuadClt != this) {
			gpuQuad.updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
		} else {
			savedQuadClt = null;
		}
		return savedQuadClt;
	}
	void restoreQuadClt(QuadCLT savedQuadClt) {
		if (savedQuadClt !=  null) {
			gpuQuad.updateQuadCLT(savedQuadClt);
		}
	}

	
	public QuadCLT(QuadCLTCPU pq, String name) {
		super (pq, name);
		if (pq instanceof QuadCLT) {
			this.gpuQuad = ((QuadCLT) pq).gpuQuad; //  careful when switching - reset Geometry, vectors, bayer images. Kernels should be the same
		}
	}
	
	public static double [] removeDisparityOutliers(
			final double [][] ds0,
			final double      max_strength,  // do not touch stronger
			final int         nth_fromextrem, // 0 - compare to max/min. 1 - second max/min, ... 
			final double      tolerance_absolute,
			final double      tolerance_relative,
			final int         width,
			final int         max_iter,
			final boolean     fit_completely, // do not add tolerance when replacing
			final int         threadsMax,
			final int         debug_level)
	{
		final int tiles = ds0[0].length;
		final double [] disparity =     ds0[0].clone();
		final double [] disparity_out = ds0[0].clone();
		final double [] strength =      ds0[1]; // will not be updated
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_updated = new AtomicInteger(0);
		final int dbg_tile = 8309;
		for (int iter = 0; iter < max_iter; iter++) {
			ai.set(0);
			anum_updated.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] neibs = new double[8];
						for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
							if ((debug_level >0) && (nTile == dbg_tile)) {
								System.out.println("removeDisparityOutliers() nTile="+nTile);
							}
							if ((strength[nTile] < max_strength) && !Double.isNaN(disparity[nTile])) {
								Arrays.fill(neibs, Double.NaN);
								for (int dir = 0; dir < 8; dir++) {
									int ineib = tn.getNeibIndex(nTile, dir);
									if (ineib >= 0) {
										neibs[dir] = disparity[ineib];
									}
								}
								Arrays.sort(neibs);
								int num_defined = neibs.length;
								for (int i = 0; i < neibs.length; i++) {
									if (Double.isNaN(neibs[i])) {
										num_defined = i;
										break;
									}
								}
								if (num_defined > 0) {
									int nth_min = nth_fromextrem;
									if (nth_min >= num_defined) {
										nth_min = num_defined - 1;
									}
									int nth_max = num_defined - 1 - nth_min;
									double d_min = neibs[nth_min] - tolerance_absolute - neibs[nth_min]*tolerance_relative;
									double d_max = neibs[nth_max] + tolerance_absolute + neibs[nth_max]*tolerance_relative;
									if  (disparity[nTile] < d_min) {
										disparity_out[nTile] = fit_completely? neibs[nth_min]: d_min;
									} else if (disparity[nTile] > d_max) {
										disparity_out[nTile] = fit_completely? neibs[nth_max]: d_max;
									}
									if (disparity_out[nTile] != disparity[nTile]) {
										anum_updated.getAndIncrement();
									}
								}
							}
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			if (anum_updated.get() ==0) {
				break;
			}
			System.arraycopy(disparity_out,0,disparity,0,tiles);
		}
		return disparity_out;
		
	}
	public static double [][] fillDisparityStrength(
			final double [][] ds0,
			final double      min_disparity,
			final double      max_sym_disparity, // lower disparity - num_bottom = 8; 
			final int         num_bottom, // average this number of lowest disparity neighbors (of 8)
			final int         num_passes,
			final double      max_change,
			final int         width,
			final int         threadsMax,
			final int         debug_level)
	{
        final double   diagonal_weight = 0.5 * Math.sqrt(2.0); // relative to ortho
		double wdiag = 0.25 *diagonal_weight / (diagonal_weight + 1.0);
		double wortho = 0.25 / (diagonal_weight + 1.0);
		final double [] neibw = {wortho, wdiag, wortho, wdiag, wortho, wdiag, wortho, wdiag}; 

		final double max_change2 = max_change * max_change;
		final int min_defined = 3; // do not try to fill if less than this (>=4 - will never fill expand outwards)	
		final double [][] ds = new double [][] {ds0[0].clone(),ds0[1].clone()};
		final int tiles = ds[0].length;
		final TileNeibs tn =  new TileNeibs(width, tiles/width);
		
		final int [] tile_indices = new int [tiles];
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger anum_gaps = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) {
						if (ds[0][nTile] < min_disparity) {
							ds[0][nTile] = min_disparity;
						}
						if (Double.isNaN(ds[0][nTile]) || (ds[1][nTile] <= 0)) {
							ds[0][nTile] = Double.NaN;
							ds[1][nTile] = 0.0;
							int indx = anum_gaps.getAndIncrement();
							tile_indices[indx] = nTile;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);
		final int num_gaps = anum_gaps.get(); 
		if (num_gaps == 0) {
			return ds; // no gaps already
		}
		final boolean [] fill_all = {false};
		final double [] disp_new = ds[0].clone();
		DoubleAccumulator amax_diff =  new DoubleAccumulator (Double::max, Double.NEGATIVE_INFINITY);

		for (int npass = 0; npass < num_passes; npass+= fill_all[0]? 1:0 ) { // do not limit initial passes
			anum_gaps.set(0);
			amax_diff.reset();
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						double [] neibs = new double[8];
						double [] neibs_sorted = new double[8];
						for (int indx = ai.getAndIncrement(); indx < num_gaps; indx = ai.getAndIncrement()) {
							int nTile = tile_indices[indx];
							if (!fill_all[0] && !Double.isNaN(ds[0][nTile])) {
								continue; // fill only new
							}
							Arrays.fill(neibs, Double.NaN);
							for (int dir = 0; dir < 8; dir++) {
								int nt_neib = tn.getNeibIndex(nTile, dir);
								if (nt_neib >= 0) {
									neibs[dir] = ds[0][nt_neib];
								}
							}
							System.arraycopy(neibs, 0, neibs_sorted, 0, 8);
							Arrays.sort(neibs_sorted); // NaNs in the end
							if (Double.isNaN(neibs_sorted[min_defined - 1])) {
								continue; // too few defined neighbors
							}
							if (!fill_all[0]) {
								anum_gaps.getAndIncrement();
							}
							double max_disp = neibs_sorted[num_bottom - 1];
							if (!(ds[0][nTile] > max_sym_disparity)){
								max_disp = Double.NaN;
							}
							double swd = 0.0, sw = 0.0;
							for (int dir = 0; dir < 8; dir++) {
								if (!Double.isNaN(neibs[dir]) && !(neibs[dir] > max_disp)) { // handles NaNs
									sw +=  neibw[dir];
									swd += neibw[dir] * neibs[dir]; 
								}
							}
							swd /= sw;  // should never be 0.0
							//max_change
							if (fill_all[0]) {
								double d2 = max_change2 * 2;
								if (!Double.isNaN(disp_new[nTile])){
									double d = swd -  disp_new[nTile];
									d2 = d * d;
								}
								amax_diff.accumulate(d2);
							}
							disp_new[nTile] = swd;		
						}
					}
				};
			}		      
			ImageDtt.startAndJoin(threads);
			ai.set(0);
			System.arraycopy(disp_new, 0, ds[0], 0, tiles);
			if ((debug_level > 0) && fill_all[0]) {
				System.out.println("fillDisparityStrength() npass="+npass+", change="+Math.sqrt(amax_diff.get())+" ("+max_change+")");
			}
			if (fill_all[0] && (amax_diff.get() < max_change2)) {
				break;
			}
			fill_all[0] = anum_gaps.get() == 0; // no new tiles filled
		} // for (int npass = 0; npass < num_passes; npass+= fill_all[0]? 1:0 )
		
		return ds;
	}

	/*
    public double [][]  getDSRBG (){
    	return dsrbg;
    }

    public void setDSRBG(
			CLTParameters  clt_parameters,
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus,
			int            debugLevel)
    {
    	setDSRBG(
    			this.dsi[is_aux?TwoQuadCLT.DSI_DISPARITY_AUX:TwoQuadCLT.DSI_DISPARITY_MAIN],
    			this.dsi[is_aux?TwoQuadCLT.DSI_STRENGTH_AUX:TwoQuadCLT.DSI_STRENGTH_MAIN],
    			this.dsi[is_aux?TwoQuadCLT.DSI_DISPARITY_AUX_LMA:TwoQuadCLT.DSI_DISPARITY_MAIN_LMA],
    			clt_parameters,
    			threadsMax,
    	        updateStatus,
    			debugLevel);
    }
    public void setDSRBG(
    		double [] disparity,
    		double [] strength,
    		double [] disparity_lma,
			CLTParameters  clt_parameters,
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus,
			int            debugLevel)
    {
    	double[][] rbg = getTileRBG(
    			clt_parameters,
    			disparity,
    			strength,
    			disparity_lma,
    			threadsMax,  // maximal number of threads to launch
    			updateStatus,
    			debugLevel);
    	double [][] dsrbg = {
    			disparity,
    			strength,
    			disparity_lma,
    			rbg[0],rbg[1],rbg[2]};
    	this.dsrbg = dsrbg;
    }
	
	public double[][] getTileRBG(
			CLTParameters  clt_parameters,
			double []      disparity,
			double []      strength,
    		double []      disparity_lma,
			int            threadsMax,  // maximal number of threads to launch
			boolean        updateStatus,
			int            debugLevel)
	{
		CLTPass3d scan = new CLTPass3d(tp);
		scan.setCalcDisparityStrength(
				disparity,
				strength);
		if (disparity_lma == null) {
			scan.resetLMA();
		} else {
			scan.setLMA(disparity_lma);
		}
		boolean [] selection = new boolean [disparity.length];
		for (int i = 0; i < disparity.length; i++) {
			selection[i] = (!Double.isNaN(disparity[i]) && ((strength == null) || (strength[i] > 0)));
		}
		scan.setTileOpDisparity(selection, disparity);
		// will work only with GPU
		// reset bayer source, geometry correction/vector
		//this.new_image_data =     true;
		QuadCLT savedQuadClt =  gpuQuad.getQuadCLT();
		if (savedQuadClt != this) {
			gpuQuad.updateQuadCLT(this);
		} else {
			savedQuadClt = null;
		}
		setPassAvgRBGA( // get image from a single pass, return relative path for x3d // USED in lwir
				clt_parameters, // CLTParameters           CLTParameters           clt_parameters,,
				scan,
				threadsMax,  // maximal number of threads to launch
				updateStatus,
				debugLevel);
		double [][] rgba = scan.getTilesRBGA();
		if (debugLevel > -1) { // -2) {
			String title = image_name+"-RBGA";
			String [] titles = {"R","B","G","A"};
			(new ShowDoubleFloatArrays()).showArrays(
					rgba,
					tp.getTilesX(),
					tp.getTilesY(),
					true,
					title,
					titles);
		}
        // Maybe resotore y caller?
		if (savedQuadClt !=  null) {
			gpuQuad.updateQuadCLT(savedQuadClt);
		}
		return rgba;
	}
	
*/	
	
	
	@Deprecated
	public double [][] getOpticalFlow(
			double disparity_min,
			double disparity_max,
			double [][] ds0,
			double [][] ds1,
			int debugLevel){
//		double sigma = .3;
		int rel_num_passes = 10;
		int  n = 2* tp.getTileSize();
		double [][] ds0f = new double[2][];
		double [][] ds1f = new double[2][];
		for (int i = 0; i < 2; i++) {
//			ds0f[i] = fillNaNGaps(ds0[i], n,sigma);
//			ds1f[i] = fillNaNGaps(ds1[i], n,sigma);

			ds0f[i] = fillNaNGaps(ds0[i], n, rel_num_passes, 100);
			ds1f[i] = fillNaNGaps(ds1[i], n, rel_num_passes, 100);
		}
		for (int tile = 0; tile < ds0f[0].length; tile++) {
			double d = ds0f[0][tile];
			if ((d < disparity_min) || (d > disparity_max)) {
				ds0f[0][tile] = Double.NaN;
				ds0f[1][tile] = Double.NaN;
			}
			d = ds1f[0][tile];
			if ((d < disparity_min) || (d > disparity_max)) {
				ds1f[0][tile] = Double.NaN;
				ds1f[1][tile] = Double.NaN;
			}
		}
		double [][] ds0ff = new double[2][];
		double [][] ds1ff = new double[2][];
		for (int i = 0; i < 2; i++) {
//			ds0ff[i] = fillNaNGaps(ds0f[i], n,sigma);
//			ds1ff[i] = fillNaNGaps(ds1f[i], n,sigma);
			ds0ff[i] = fillNaNGaps(ds0f[i], n, rel_num_passes, 100);
			ds1ff[i] = fillNaNGaps(ds1f[i], n, rel_num_passes, 100);
		}
		
		
		if (debugLevel > 0) {
			double [][] dbg_img = {
//					ds0[1],ds0f[1],ds0ff[1],ds1[1],ds1f[1],ds1ff[1],
//					ds0[0],ds0f[0],ds0ff[0],ds1[0],ds1f[0],ds1ff[0]};
					ds0[1],ds1[1],ds0f[1],ds1f[1],ds0ff[1],ds1ff[1],
					ds0[0],ds1[0],ds0f[0],ds1f[0],ds0ff[0],ds1ff[0]};
			String title = "filled_nan_gaps";
//			String [] titles = {
//					"ds0[1]","ds0f[1]","ds0ff[1]","ds1[1]","ds1f[1]","ds1ff[1]",
//					"ds0[0]","ds0f[0]","ds0ff[0]","ds1[0]","ds1f[0]","ds1ff[0]"};
			String [] titles = {
					"ds0[1]","ds1[1]","ds0f[1]","ds1f[1]","ds0ff[1]","ds1ff[1]",
					"ds0[0]","ds1[0]","ds0f[0]","ds1f[0]","ds0ff[0]","ds1ff[0]"};

			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tp.getTilesX(),
					tp.getTilesY(),
					true,
					title,
					titles);
		}

		
		
		return null;
	}

	public static double[] blurWeak(
			double [][] ds,
			double      min_strength_blur,
			double      min_strength_replace,
			int         n,
			int         width,
			double      sigma) {
		double [] disparity = new double[ds[0].length];
		Arrays.fill(disparity, Double.NaN);
		for (int i = 0; i < disparity.length; i++) {
			if (Double.isNaN(ds[0][i]) || (ds[1][i] < min_strength_blur)) {
				disparity[i] = ds[0][i];
			}
		}
		for (int i = 0; i < n; i++) { 
			disparity =   (new DoubleGaussianBlur()).blurWithNaN(
					disparity,              // double[] pixels,
					null,                   // double [] in_weight, // or null
					width,                  // int width,
					disparity.length/width, // int height,
					sigma,                  // double sigmaX,
					sigma,                  // double sigmaY,
					0.00001);               // double accuracy);
		}
		for (int i = 0; i < disparity.length; i++) {
			if (ds[1][i] > min_strength_replace) {
				disparity[i] = ds[0][i];
			}
		}
		return disparity;
	}
	
	
	
	@Deprecated
	public double[] fillNaNGapsOld(
			double [] data,
			int       n,
			double sigma) {
		double [] d = data.clone();
		for (int i = 0; i < n; i++) { 
			d =   (new DoubleGaussianBlur()).blurWithNaN(
					d, // double[] pixels,
					null,  // double [] in_weight, // or null
					tp.getTilesX(), // int width,
					tp.getTilesY(), // int height,
					sigma, // double sigmaX,
					sigma, // double sigmaY,
					0.00001); // double accuracy);
		}
		for (int i = 0; i < data.length; i++) {
			if (!Double.isNaN(data[i])) {
				d[i] = data[i];
			}
		}
		return d;
	}
	@Deprecated	
	public double[] fillNaNGaps(
			double [] data,
			int       n,
			int       rel_num_passes,
			int       threadsMax)
	{
		int num_passes = n * rel_num_passes;
		
		return tp.fillNaNs(
				data,                 // double [] data,
				tp.getTilesX(),       //int       width, 
				2 * n,                // int       grow,
				0.5 * Math.sqrt(2.0), // double    diagonal_weight, // relative to ortho
				num_passes,           // int       num_passes,
				threadsMax); // final int threadsMax)      // maximal number of threads to launch                         
	}
	
	@Deprecated
	public double[][][]  get_pair( // not used
			double k_prev,
			QuadCLT qprev,
			double corr_scale, //  = 0.75
			int debug_level)
	{
		final int iscale = 8;
		double ts =     getTimeStamp();
		double ts_prev =   ts;
		double [] zero3 =  {0.0,0.0,0.0};	
		double [] camera_xyz0 = zero3.clone();
		double [] camera_atr0 = zero3.clone();
		
		ErsCorrection ersCorrection = (ErsCorrection) geometryCorrection;
		
		System.out.println("\n"+image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", image_name,
				ersCorrection.ers_wxyz_center[0], ersCorrection.ers_wxyz_center[1],ersCorrection.ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f",	image_name,
				ersCorrection.ers_wxyz_center_dt[0], ersCorrection.ers_wxyz_center_dt[1],ersCorrection.ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", image_name,
				ersCorrection.ers_wxyz_center_d2t[0], ersCorrection.ers_wxyz_center_d2t[1],ersCorrection.ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", image_name,
				ersCorrection.ers_watr_center_dt[0], ersCorrection.ers_watr_center_dt[1],ersCorrection.ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", image_name,
				ersCorrection.ers_watr_center_d2t[0], ersCorrection.ers_watr_center_d2t[1],ersCorrection.ers_watr_center_d2t[2] ));
		
		double dt = 0.0;
		if (qprev == null) {
			qprev = this;
		}
		if (qprev != null) {
			ts_prev = qprev.getTimeStamp();
			dt = ts-ts_prev;
			if (dt < 0) {
				k_prev = (1.0-k_prev);
			}
			if (Math.abs(dt) > 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
				k_prev = 0.5;
				System.out.println("Non-consecutive frames, dt = "+dt);
			}
			ErsCorrection ersCorrectionPrev = (ErsCorrection) (qprev.geometryCorrection);
			double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
			double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt;
			double [] wxyz_delta = new double[3];
			double [] watr_delta = new double[3];
			for (int i = 0; i <3; i++) {
				wxyz_delta[i] = - corr_scale * dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_wxyz_center_dt[i]);
				watr_delta[i] = - corr_scale * dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_watr_center_dt[i]);
			}
			camera_xyz0 = wxyz_delta;
			camera_atr0 = watr_delta;
		}
		
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		
		double [][] dsrbg = this.transformCameraVew(
				qprev,       // QuadCLT   camera_QuadClt,
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr, // camera orientation relative to world frame
				iscale);
		double [][][] pair = {getDSRBG(),dsrbg};
		
		// combine this scene with warped previous one
		if (debug_level > 0) {
			String [] rtitles = new String[2* dsrbg_titles.length];
			double [][] dbg_rslt = new double [rtitles.length][];
			for (int i = 0; i < dsrbg_titles.length; i++) {
				rtitles[2*i] =    dsrbg_titles[i]+"0";
				rtitles[2*i+1] =  dsrbg_titles[i];
				dbg_rslt[2*i] =   pair[0][i];
				dbg_rslt[2*i+1] = pair[1][i];
			}
			String title = image_name+"-"+qprev.image_name+"-dt"+dt;
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_rslt,
					tilesX,
					tilesY,
					true,
					title,
					rtitles);
		}		
		return pair;
	}

	@Deprecated
	public double [][] transformCameraVew(
			QuadCLT   camera_QuadClt,
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr, // camera orientation relative to world frame
			int       iscale)
	{
		double    line_err = 0.1; // 10.0; // 0.1; // BUG
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		int rel_num_passes = 10;
		int num_passes =    transform_size; // * 2;

		double [] zero3 =               {0.0,0.0,0.0};	
		int stilesX = iscale*tilesX; 
		int stilesY = iscale*tilesY;
		int stiles = stilesX*stilesY;
		double sigma = 0.5 * iscale;
		double scale =  1.0 * iscale/transform_size;
		double [][] dsrbg_camera = camera_QuadClt.getDSRBG();
		double [][] ds =        new double [dsrbg_camera.length][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
			}
		}
		
		ErsCorrection ersCorrection = getErsCorrection();
		ersCorrection.setupERS(); // just in case - setUP using instance paRAMETERS
		double [] zbuffer = new double [tiles];
		for (int tileY = 0; tileY < tilesY; tileY++) {
			int stileY = iscale * tileY + iscale/2;
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int stileX = iscale * tileX + iscale/2;
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg_camera[DSRBG_DISPARITY][nTile];
				if (disparity < 0) {
					disparity = 0.0;
				}
				// found that there are tiles with strength == 0.0, while disparity is not NaN
				if (!Double.isNaN(disparity) && (dsrbg_camera[DSRBG_STRENGTH][nTile] > 0.0)) {
					double [] pXpYD = ersCorrection.getImageCoordinatesERS(
							camera_QuadClt, // QuadCLT cameraQuadCLT, // camera station that got image to be to be matched 
							centerX,        // double px,                // pixel coordinate X in this camera view
							centerY,        //double py,                // pixel coordinate Y in this camera view
							disparity,      // double disparity,         // this view disparity 
							true,           // boolean distortedView,    // This camera view is distorted (diff.rect), false - rectilinear
							zero3,          // double [] reference_xyz,  // this view position in world coordinates (typically zero3)
							zero3,          // double [] reference_atr,  // this view orientation relative to world frame  (typically zero3)
							true,           // boolean distortedCamera,  // camera view is distorted (false - rectilinear)
							camera_xyz,     // double [] camera_xyz,     // camera center in world coordinates
							camera_atr,     // double [] camera_atr,     // camera orientation relative to world frame
							line_err);       // double    line_err)       // threshold error in scan lines (1.0)
					if (pXpYD != null) {
						int px = (int) Math.round(pXpYD[0]/transform_size);
						int py = (int) Math.round(pXpYD[1]/transform_size);
						int spx = (int) Math.round(pXpYD[0]*scale);
						int spy = (int) Math.round(pXpYD[1]*scale);
						if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
							//Z-buffer
							if (!(pXpYD[2] < zbuffer[px + py* tilesX])) {
								zbuffer[px + py* tilesX] = pXpYD[2];
								if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
									int sTile = spx + spy* stilesX;
									ds[DSRBG_DISPARITY][sTile] = pXpYD[2]; //reduce*
									for (int i = DSRBG_STRENGTH; i < dsrbg_camera.length; i++) {
										ds[i][sTile] = dsrbg_camera[i][nTile]; // reduce * 
									}
								}								
							}
						}
					}
				}
			}
		}
		
		//dsrbg_out[DSRBG_DISPARITY]
		for (int i = 0; i < ds.length; i++) {
			ds[i] = (new DoubleGaussianBlur()).blurWithNaN(
					ds[i], // double[] pixels,
					null,  // double [] in_weight, // or null
					stilesX, // int width,
					stilesY, // int height,
					sigma, // double sigmaX,
					sigma, // double sigmaY,
					0.01); // double accuracy);
		}
		double [][] dsrbg_out = new double [dsrbg_camera.length][tiles];
		int [][] num_non_nan = new int [dsrbg_out.length] [tiles];
		
		for (int stileY = 0; stileY < stilesY; stileY++) {
			int tileY = stileY / iscale; 
			for (int stileX = 0; stileX < stilesX; stileX++) {
				int tileX = stileX / iscale;
				int stile = stileX + stileY * stilesX;
				int tile =  tileX +  tileY *  tilesX;
				for (int i = 0; i < dsrbg_out.length; i++) {
					double d = ds[i][stile];
					if (!Double.isNaN(d)) {
						num_non_nan[i][tile] ++;
						dsrbg_out[i][tile] += d;
					}	
				}
			}
		}
		for (int i = 0; i < dsrbg_out.length; i++) {
			for (int j = 0; j < tiles; j++) {
				if (num_non_nan[i][j] == 0) {
					dsrbg_out[i][j] = Double.NaN;
				} else {
					dsrbg_out[i][j]/=num_non_nan[i][j];
				}
			}
		}

		if (num_passes > 0) {
			for (int i = 0; i < dsrbg_out.length; i++) {
				dsrbg_out[i] = fillNaNGaps(dsrbg_out[i], num_passes, rel_num_passes, 100); // threadsMax);
			}
		}
		return dsrbg_out;
	}

	@Deprecated
	public double[][][]  get_pair_old(
			double k_prev,
			QuadCLT qprev,
			double corr_scale, //  = 0.75
			int debug_level)
	{
		final int iscale = 8;
		double ts =     getTimeStamp();
		double ts_prev =   ts;
		double [] zero3 =  {0.0,0.0,0.0};	
		double [] camera_xyz0 = zero3.clone();
		double [] camera_atr0 = zero3.clone();
		
		ErsCorrection ersCorrection = (ErsCorrection) geometryCorrection;
		
		System.out.println("\n"+image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", image_name,
				ersCorrection.ers_wxyz_center[0], ersCorrection.ers_wxyz_center[1],ersCorrection.ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f",	image_name,
				ersCorrection.ers_wxyz_center_dt[0], ersCorrection.ers_wxyz_center_dt[1],ersCorrection.ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", image_name,
				ersCorrection.ers_wxyz_center_d2t[0], ersCorrection.ers_wxyz_center_d2t[1],ersCorrection.ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", image_name,
				ersCorrection.ers_watr_center_dt[0], ersCorrection.ers_watr_center_dt[1],ersCorrection.ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", image_name,
				ersCorrection.ers_watr_center_d2t[0], ersCorrection.ers_watr_center_d2t[1],ersCorrection.ers_watr_center_d2t[2] ));
		
		double dt = 0.0;
		if (qprev == null) {
			qprev = this;
		}
		if (qprev != null) {
			ts_prev = qprev.getTimeStamp();
			dt = ts-ts_prev;
			if (dt < 0) {
				k_prev = (1.0-k_prev);
			}
			if (Math.abs(dt) > 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
				k_prev = 0.5;
				System.out.println("Non-consecutive frames, dt = "+dt);
			}
			ErsCorrection ersCorrectionPrev = (ErsCorrection) (qprev.geometryCorrection);
			double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
			double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt;
			double [] wxyz_delta = new double[3];
			double [] watr_delta = new double[3];
			for (int i = 0; i <3; i++) {
				wxyz_delta[i] = - corr_scale * dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_wxyz_center_dt[i]);
				watr_delta[i] = - corr_scale * dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * ersCorrection.ers_watr_center_dt[i]);
			}
			camera_xyz0 = wxyz_delta;
			camera_atr0 = watr_delta;
		}
		
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] dsrbg_titles = {"d", "s", "r", "b", "g"};
		
		double [][] dsrbg0 = qprev.transformCameraVew( // previous camera view, ers only 
				zero3,       // double [] camera_xyz, // camera center in world coordinates
				zero3,       // double [] camera_atr,  // camera orientation relative to world frame
				iscale);

		if (debug_level > 1) {
			String title0 = String.format("%s_DSRBG_ers-only",qprev.image_name); // previous frame, ERS-corrected only, no shift/rot
			(new ShowDoubleFloatArrays()).showArrays(
					dsrbg0,
					tilesX,
					tilesY,
					true,
					title0,
					dsrbg_titles);
		}
		double [][] dsrbg = transformCameraVew(
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr,  // camera orientation relative to world frame
				iscale);      // 
		if (debug_level > 1) {
			String title = String.format("%s_%f:%f:%f_%f:%f:%f",image_name,
					camera_xyz0[0],camera_xyz0[1],camera_xyz0[2],camera_atr0[0],camera_atr0[1],camera_atr0[2]);
			(new ShowDoubleFloatArrays()).showArrays(
					dsrbg,
					tilesX,
					tilesY,
					true,
					title,
					dsrbg_titles);
		}
		
		double [][] dsrbg1 = this.transformCameraVew(
				qprev,       // QuadCLT   camera_QuadClt,
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr, // camera orientation relative to world frame
				iscale);
///		double [][][] pair = {dsrbg0,dsrbg};
///		double [][][] pair = {qprev.getDSRBG(),dsrbg};
///		double [][][] pair = {qprev.getDSRBG(),dsrbg1};
		double [][][] pair = {getDSRBG(),dsrbg1};
		
		// combine previous frame with this one
		if (debug_level > 0) {
			String [] rtitles = new String[2* dsrbg_titles.length];
			double [][] dbg_rslt = new double [rtitles.length][];
			for (int i = 0; i < dsrbg_titles.length; i++) {
				rtitles[2*i] =    dsrbg_titles[i]+"0";
				rtitles[2*i+1] =  dsrbg_titles[i];
				dbg_rslt[2*i] =   pair[0][i];
				dbg_rslt[2*i+1] = pair[1][i];
			}
			String title = image_name+"-"+qprev.image_name+"-dt"+dt;
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_rslt,
					tilesX,
					tilesY,
					true,
					title,
					rtitles);
		}		
		return pair;
	}
	
	
	@Deprecated
	public double [][] transformCameraVew(
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr, // camera orientation relative to world frame
			int       iscale)
	{
		double [][] dsrbg =   getDSRBG(); 
		return transformCameraVew(
				dsrbg,      // double [][] dsi, // 
				camera_xyz,  // double [] camera_xyz, // camera center in world coordinates
				camera_atr,  // double [] camera_atr,  // camera orientation relative to world frame
				iscale,
	// normally all 0; 			
				null, // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				null, // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				null, // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				null, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				null); // double [] watr_center_d2t) // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	}
	
	@Deprecated
	public double [][] transformCameraVew(
			double [][] dsrbg, // 
			double []   camera_xyz, // camera center in world coordinates
			double []   camera_atr,  // camera orientation relative to world frame
			final int   iscale,     // 8 
// normally all 0; 			
			double [] wxyz_center,     // world camera XYZ (meters) for the frame center
			double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
			double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
			double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
			double [] watr_center_d2t) // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	{
		double    line_err = 10.0; // 0.1; // BUG
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		int rel_num_passes = 10;
		int num_passes =    transform_size; // * 2;

		double [] zero3 =               {0.0,0.0,0.0};	
		int stilesX = iscale*tilesX; 
		int stilesY = iscale*tilesY;
		int stiles = stilesX*stilesY;
		double sigma = 0.5 * iscale;
		double scale =  1.0 * iscale/transform_size;
//		double reduce = 1.0 / iscale/iscale;
		double [][] ds =        new double [dsrbg.length][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
			}
		}

		ErsCorrection ersCorrection = (ErsCorrection) geometryCorrection;
		// save original
		double [] saved_ers_wxyz_center =     ersCorrection.ers_wxyz_center;
		double [] saved_ers_wxyz_center_dt =  ersCorrection.ers_wxyz_center_dt;
		double [] saved_ers_wxyz_center_d2t = ersCorrection.ers_wxyz_center_d2t;
		double [] saved_ers_watr_center_dt =  ersCorrection.ers_watr_center_dt;
		double [] saved_ers_watr_center_d2t = ersCorrection.ers_watr_center_d2t;
		
		if (wxyz_center ==     null) wxyz_center =     ersCorrection.ers_wxyz_center;
		if (wxyz_center_dt ==  null) wxyz_center_dt =  ersCorrection.ers_wxyz_center_dt;
		if (wxyz_center_d2t == null) wxyz_center_d2t = ersCorrection.ers_wxyz_center_d2t;
		if (watr_center_dt ==  null) watr_center_dt =  ersCorrection.ers_watr_center_dt;
		if (watr_center_d2t == null) watr_center_d2t = ersCorrection.ers_watr_center_d2t;
		
		ersCorrection.setupERS(
				wxyz_center,      // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				wxyz_center_dt,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				wxyz_center_d2t,  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				watr_center_dt,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		double [][] wxyz = new double [tiles][];
		for (int tileY = 0; tileY < tilesY; tileY++) {
//			int stileY = iscale * tileY + iscale/2;
			for (int tileX = 0; tileX < tilesX; tileX++) {
//				int stileX = iscale * tileX + iscale/2;
//				int stile = stileX + stilesX * stileY;
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				double disparity = dsrbg[DSRBG_DISPARITY][nTile];
				if (disparity < 0) {
					disparity = 0.0;
				}
				// found that there are tiles with strength == 0.0, while disparity is not NaN
				if (!Double.isNaN(disparity) && (dsrbg[DSRBG_STRENGTH][nTile] > 0.0)) {
					wxyz[nTile] = ersCorrection.getWorldCoordinatesERS(
							centerX,                                        // double px,
							centerY,                                        // double py,
							disparity, // double disparity,
							true,                                           // boolean correctDistortions,// correct distortion (will need corrected background too !)
							camera_xyz,  // double [] camera_xyz, // camera center in world coordinates
							camera_atr); // double [] camera_atr)  // camera orientation relative to world frame
				}
			}
		}
		
		ersCorrection.setupERS(
		zero3,   // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
		zero3,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
		zero3,   // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
		zero3,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		zero3);  // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		double [] zbuffer = new double [tiles];
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
//				if (!Double.isNaN(dbg_img[indx_xyz][nTile])) {
				if (wxyz[nTile] != null) {
					double [] pXpYD = ersCorrection.getImageCoordinatesERS(
							wxyz[nTile], // double [] wxyz, 
							true,                                           // boolean correctDistortions,// correct distortion (will need corrected background too !)
							zero3,       // double [] camera_xyz, // camera center in world coordinates
							zero3,       // double [] camera_atr)  // camera orientation relative to world frame
							line_err);   // double    line_err) // threshold error in scan lines (1.0)
					if (pXpYD != null) {
						int px = (int) Math.round(pXpYD[0]/transform_size);
						int py = (int) Math.round(pXpYD[1]/transform_size);
						int spx = (int) Math.round(pXpYD[0]*scale);
						int spy = (int) Math.round(pXpYD[1]*scale);
						if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
							//Z-buffer
							if (!(pXpYD[2] < zbuffer[px + py* tilesX])) {
								zbuffer[px + py* tilesX] = pXpYD[2];
								if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
									int sTile = spx + spy* stilesX;
									ds[DSRBG_DISPARITY][sTile] = pXpYD[2]; //reduce*
									for (int i = DSRBG_STRENGTH; i < dsrbg.length; i++) {
										ds[i][sTile] = dsrbg[i][nTile]; // reduce * 
									}
								}								
							}
						}
					}
				}
			}
		}
		
		//dsrbg_out[DSRBG_DISPARITY]
		for (int i = 0; i < ds.length; i++) {
//			int rel_num_passes = 10;
//			int num_passes =  2 * transform_size;
//			ds[i] = fillNaNGaps(ds[i], num_passes, rel_num_passes, 100); // threadsMax);
			/* */
			ds[i] = (new DoubleGaussianBlur()).blurWithNaN(
					ds[i], // double[] pixels,
					null,  // double [] in_weight, // or null
					stilesX, // int width,
					stilesY, // int height,
					sigma, // double sigmaX,
					sigma, // double sigmaY,
					0.01); // double accuracy);
				/*	*/
		}
		double [][] dsrbg_out = new double [dsrbg.length][tiles];
		int [][] num_non_nan = new int [dsrbg_out.length] [tiles];
		
		for (int stileY = 0; stileY < stilesY; stileY++) {
			int tileY = stileY / iscale; 
			for (int stileX = 0; stileX < stilesX; stileX++) {
				int tileX = stileX / iscale;
				int stile = stileX + stileY * stilesX;
				int tile =  tileX +  tileY *  tilesX;
				for (int i = 0; i < dsrbg_out.length; i++) {
					double d = ds[i][stile];
					if (!Double.isNaN(d)) {
						num_non_nan[i][tile] ++;
						dsrbg_out[i][tile] += d;
					}	
				}
			}
		}
		for (int i = 0; i < dsrbg_out.length; i++) {
			for (int j = 0; j < tiles; j++) {
				if (num_non_nan[i][j] == 0) {
					dsrbg_out[i][j] = Double.NaN;
				} else {
					dsrbg_out[i][j]/=num_non_nan[i][j];
				}
			}
		}

		if (num_passes > 0) {
			for (int i = 0; i < dsrbg_out.length; i++) {
				dsrbg_out[i] = fillNaNGaps(dsrbg_out[i], num_passes, rel_num_passes, 100); // threadsMax);
			}
		}
		
		// restore original
		ersCorrection.setupERS(
				saved_ers_wxyz_center,      // double [] wxyz_center,      // world camera XYZ (meters) for the frame center
				saved_ers_wxyz_center_dt,   // double [] wxyz_center_dt,   // world camera Vx, Vy, Vz (m/s)
				saved_ers_wxyz_center_d2t,  // double [] wxyz_center_d2t,  // world camera Vx, Vy, Vz (m/s^2)
				saved_ers_watr_center_dt,   // double [] watr_center_dt,   // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				saved_ers_watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		return dsrbg_out;
	}
	
	@Deprecated
	public double[][]  test_back_forth_debug(
			double k_prev,
			QuadCLT qprev,
			double corr_scale, //  = 0.75
			int debug_level)
	{
		double ts =     getTimeStamp();
		double ts_prev = ts;
		boolean force_manual = false;
		double [] zero3 =               {0.0,0.0,0.0};	
		double [] camera_xyz0 = {0.0,0.0,0.0};
		double [] camera_atr0 = {0.0,0.0,0.0};

		double [] watr_center_dt =  {0.0, 0.0, 0.0};
		double [] watr_center_d2t = {0.0, 0.0, 0.0};
		double [] wxyz_center =     {0.0, 0.0, 0.0};	
		double [] wxyz_center_dt =  {0.0, 0.0, 0.0};	
		double [] wxyz_center_d2t = {0.0, 0.0, 0.0};	
		
		ErsCorrection ersCorrection = (ErsCorrection) geometryCorrection;
		// save original
		double [] saved_ers_wxyz_center =     ersCorrection.ers_wxyz_center;
		double [] saved_ers_wxyz_center_dt =  ersCorrection.ers_wxyz_center_dt;
		double [] saved_ers_wxyz_center_d2t = ersCorrection.ers_wxyz_center_d2t;
		double [] saved_ers_watr_center_dt =  ersCorrection.ers_watr_center_dt;
		double [] saved_ers_watr_center_d2t = ersCorrection.ers_watr_center_d2t;
		
		System.out.println("\n"+image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", image_name,saved_ers_wxyz_center[0], saved_ers_wxyz_center[1],saved_ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f", image_name,saved_ers_wxyz_center_dt[0], saved_ers_wxyz_center_dt[1],saved_ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", image_name,saved_ers_wxyz_center_d2t[0], saved_ers_wxyz_center_d2t[1],saved_ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", image_name,saved_ers_watr_center_dt[0], saved_ers_watr_center_dt[1],saved_ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", image_name,saved_ers_watr_center_d2t[0], saved_ers_watr_center_d2t[1],saved_ers_watr_center_d2t[2] ));

		// Comment out to use manual values
		
		if (force_manual) {
			System.out.println("\nUsing manually set values:");
			System.out.println(String.format("%s: wxyz_center=     %f, %f, %f",     image_name,wxyz_center[0], wxyz_center[1],wxyz_center[2] ));
			System.out.println(String.format("%s: wxyz_center_dt=  %f, %f, %f",  image_name,wxyz_center_dt[0], wxyz_center_dt[1],wxyz_center_dt[2] ));
			System.out.println(String.format("%s: wxyz_center_d2t= %f, %f, %f", image_name,wxyz_center_d2t[0], wxyz_center_d2t[1],wxyz_center_d2t[2] ));
			System.out.println(String.format("%s: watr_center_dt=  %f, %f, %f",  image_name,watr_center_dt[0], watr_center_dt[1],watr_center_dt[2] ));
			System.out.println(String.format("%s: watr_center_d2t= %f, %f, %f", image_name,watr_center_d2t[0], watr_center_d2t[1],watr_center_d2t[2] ));
		}else {
			wxyz_center =     saved_ers_wxyz_center;
			wxyz_center_dt =  saved_ers_wxyz_center_dt;
			wxyz_center_d2t = saved_ers_wxyz_center_d2t;  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
			watr_center_dt =  saved_ers_watr_center_dt;
			watr_center_d2t = saved_ers_watr_center_d2t;
		}
		
		
		double dt = 0.0;
		if (qprev == null) {
			qprev = this;
		}
		if (qprev != null) {
			ts_prev = qprev.getTimeStamp();
			dt = ts-ts_prev;
			if (dt < 0) {
				k_prev = (1.0-k_prev);
			}
			if (Math.abs(dt) > 0.15) { // at least two frames TODO: use number of lines* line_time * ...? 
				k_prev = 0.5;
				System.out.println("Non-consecutive frames, dt = "+dt);
			}
			ErsCorrection ersCorrectionPrev = (ErsCorrection) (qprev.geometryCorrection);
			double [] wxyz_center_dt_prev =   ersCorrectionPrev.ers_wxyz_center_dt;
			double [] watr_center_dt_prev =   ersCorrectionPrev.ers_watr_center_dt;
			double [] wxyz_delta = new double[3];
			double [] watr_delta = new double[3];
			for (int i = 0; i <3; i++) {
				wxyz_delta[i] = - corr_scale * dt * (k_prev * wxyz_center_dt_prev[i] + (1.0-k_prev) * saved_ers_wxyz_center_dt[i]);
				watr_delta[i] = - corr_scale * dt * (k_prev * watr_center_dt_prev[i] + (1.0-k_prev) * saved_ers_watr_center_dt[i]);
			}
			camera_xyz0 = wxyz_delta;
			camera_atr0 = watr_delta;
		}
		
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		String [] titles = {"s0", "s1", "d0", "d1", "pX0","pX", "pY0","pY", "x","y","z","w"};
		double [][] dsi = {this.dsi[TwoQuadCLT.DSI_DISPARITY_MAIN], this.dsi[TwoQuadCLT.DSI_STRENGTH_MAIN]};
		double [][] dbg_img0 = new double [titles.length][]; // null; // previous camera view, ers only 
		double [][] dbg_img1 = new double [titles.length][]; // null; // this camera view, ers only
		
		if (!force_manual) { 
			dbg_img0 = qprev.transformCameraVewDebug( // previous camera view, ers only 
					zero3,       // double [] camera_xyz, // camera center in world coordinates
					zero3);       // double [] camera_atr,  // camera orientation relative to world frame
			dbg_img1 = transformCameraVewDebug(     // this camera view, ers only
					zero3,       // double [] camera_xyz, // camera center in world coordinates
					zero3);       // double [] camera_atr,  // camera orientation relative to world frame

			if (debug_level > 1) {
				String title0 = String.format("%s_ers_only",qprev.image_name); // previous frame, ERS-corrected only, no shift/rot
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img0,
						tilesX,
						tilesY,
						true,
						title0,
						titles);
			}
		}
		double [][] dbg_img = transformCameraVewDebug(
				dsi,         // 
				camera_xyz0, // double [] camera_xyz, // camera center in world coordinates
				camera_atr0, //double [] camera_atr,  // camera orientation relative to world frame
				wxyz_center,      // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				wxyz_center_dt,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				wxyz_center_d2t,  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				watr_center_dt,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		if (debug_level > 1) {
			String title = String.format("%s_%f:%f:%f_%f:%f:%f",image_name,
					camera_xyz0[0],camera_xyz0[1],camera_xyz0[2],camera_atr0[0],camera_atr0[1],camera_atr0[2]);
			if (!force_manual) {
				title += String.format ("_ers_%f:%f:%f_%f:%f:%f",wxyz_center_dt[0],wxyz_center_dt[1],wxyz_center_dt[2],
						watr_center_dt[0],watr_center_dt[1],watr_center_dt[2]);
			}
			(new ShowDoubleFloatArrays()).showArrays(
					dbg_img,
					tilesX,
					tilesY,
					true,
					title,
					titles);
		}
		int indx_s0 =  0;
		int indx_s1 =  1;
		int indx_d0 =  2;
		int indx_d1 =  3;
		// combine previous frame with this one
		String [] rtitles = {"s","s_ers", "s_prev","s_this",
				            "d","d_ers", "d_prev","d_this"};
		double [][] rslt = {
				dbg_img1[indx_s0], // this strength as is
				dbg_img1[indx_s1], // this strength - ers only
				dbg_img0[indx_s1], // previous strength - ers only to match to
				dbg_img [indx_s1], // this strength - ers and shift/rot
				dbg_img1[indx_d0], // this disparity as is
				dbg_img1[indx_d1], // this disparity - ers only
				dbg_img0[indx_d1], // previous disparity - ers only to match to
				dbg_img [indx_d1]};// this disparity - ers and shift/rot
		if (debug_level > 0) {
			String title = image_name+"-"+qprev.image_name+"-dt"+dt;
			(new ShowDoubleFloatArrays()).showArrays(
					rslt,
					tilesX,
					tilesY,
					true,
					title,
					rtitles);
		}		
		return rslt;
	}
	
	@Deprecated
	public double [][] transformCameraVewDebug(
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr)  // camera orientation relative to world frame
	{
		
		double [][] dsi = {this.dsi[TwoQuadCLT.DSI_DISPARITY_MAIN], this.dsi[TwoQuadCLT.DSI_STRENGTH_MAIN]};

		return transformCameraVewDebug(
				dsi, // 
				camera_xyz, // camera center in world coordinates
				camera_atr,  // camera orientation relative to world frame
	// normally all 0; 			
				null, // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				null, // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				null, // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				null, // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				null); // double [] watr_center_d2t) // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
	}

	@Deprecated
	public double [][] transformCameraVewDebug(
			double [][] dsi, // 
			double [] camera_xyz, // camera center in world coordinates
			double [] camera_atr,  // camera orientation relative to world frame
			
// normally all 0; 			
			double [] wxyz_center,     // world camera XYZ (meters) for the frame center
			double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
			double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
			double [] watr_center_dt,  // camera rotations (az, tilt, roll in radians/s, corresponding to the frame center)
			double [] watr_center_d2t) // camera rotations (az, tilt, roll in radians/s, corresponding to the frame center)
	{
		double    line_err = 10.0; // 0.1; // bug
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
//		double [] ers_watr_center_dt =  {0.168792, 0.037145, -0.060279};
		double [] zero3 =               {0.0,0.0,0.0};	
		int indx_s0 =  0;
		int indx_s1 =  1;
		int indx_d0 =  2;
		int indx_d1 =  3;
///		int indx_pX0 = 4;
		int indx_pX =  5;
///		int indx_pY0 = 6;
		int indx_pY =  7;
		int indx_xyz = 8;
		double [][] dbg_img = new double [12][tiles];
		for (int i = 0; i <dbg_img.length; i++) {
			for (int j = 0; j <dbg_img[i].length; j++) {
				dbg_img[i][j] = Double.NaN;
			}
		}
		final int iscale = 8;
		int stilesX = iscale*tilesX; 
		int stilesY = iscale*tilesY;
		int stiles = stilesX*stilesY;
		double sigma = 0.5 * iscale;
		double scale = 1.0*iscale/transform_size;
		double reduce = 1.0/iscale/iscale;
		double [][] ds = new double [2][stiles];
		double [][] ds0 = new double [2][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
				ds0[i][j] = Double.NaN;
			}
		}
		
		dbg_img [indx_d0] = dsi[0];
		dbg_img [indx_s0] = dsi[1];

		ErsCorrection ersCorrection = (ErsCorrection) geometryCorrection;
		// save original
		double [] saved_ers_wxyz_center =     ersCorrection.ers_wxyz_center;
		double [] saved_ers_wxyz_center_dt =  ersCorrection.ers_wxyz_center_dt;
		double [] saved_ers_wxyz_center_d2t = ersCorrection.ers_wxyz_center_d2t;
		double [] saved_ers_watr_center_dt =  ersCorrection.ers_watr_center_dt;
		double [] saved_ers_watr_center_d2t = ersCorrection.ers_watr_center_d2t;
		
		if (wxyz_center ==     null) wxyz_center =     ersCorrection.ers_wxyz_center;
		if (wxyz_center_dt ==  null) wxyz_center_dt =  ersCorrection.ers_wxyz_center_dt;
		if (wxyz_center_d2t == null) wxyz_center_d2t = ersCorrection.ers_wxyz_center_d2t;
		if (watr_center_dt ==  null) watr_center_dt =  ersCorrection.ers_watr_center_dt;
		if (watr_center_d2t == null) watr_center_d2t = ersCorrection.ers_watr_center_d2t;
		
		
		ersCorrection.setupERS(
				wxyz_center,      // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				wxyz_center_dt,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				wxyz_center_d2t,  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				watr_center_dt,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		
		for (int tileY = 0; tileY < tilesY; tileY++) {
			int stileY = iscale * tileY + iscale/2;
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int stileX = iscale * tileX + iscale/2;
				int stile = stileX + stilesX * stileY;
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				ds0[0][stile] = reduce* dbg_img [indx_d0][nTile]; //  this.dsi[TwoQuadCLT.DSI_DISPARITY_MAIN][nTile];
				ds0[1][stile] = reduce* dbg_img [indx_s0][nTile]; // this.dsi[TwoQuadCLT.DSI_STRENGTH_MAIN][nTile];
				
				double [] wxyz = ersCorrection.getWorldCoordinatesERS(
						centerX,                                        // double px,
						centerY,                                        // double py,
						dsi[0][nTile], // double disparity,
						true,                                           // boolean correctDistortions,// correct distortion (will need corrected background too !)
						camera_xyz,  // double [] camera_xyz, // camera center in world coordinates
						camera_atr); // double [] camera_atr)  // camera orientation relative to world frame
				if ((wxyz != null) && (wxyz[2] < 0.0)) {
					for (int i = 0; i < 4; i ++) {
						dbg_img[indx_xyz+i][nTile] = wxyz[i]; 
					}
				}
			}
		}
		
		ersCorrection.setupERS(
		zero3,   // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
		zero3,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
		zero3,   // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
		zero3,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		zero3);  // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				if (!Double.isNaN(dbg_img[indx_xyz][nTile])) {
					double [] wxyz = {
							dbg_img[indx_xyz][nTile],
							dbg_img[indx_xyz+1][nTile],
							dbg_img[indx_xyz+2][nTile],
							dbg_img[indx_xyz+3][nTile]};
					double [] pXpYD = ersCorrection.getImageCoordinatesERS(
							wxyz, // double [] wxyz, 
							true,                                           // boolean correctDistortions,// correct distortion (will need corrected background too !)
							zero3,  // double [] camera_xyz, // camera center in world coordinates
							zero3,  // double [] camera_atr)  // camera orientation relative to world frame
							line_err);   // double    line_err) // threshold error in scan lines (1.0)
					if (pXpYD != null) {
						dbg_img[indx_pX][nTile] = pXpYD[0];
						dbg_img[indx_pY][nTile] = pXpYD[1];
						
						int px = (int) Math.round(pXpYD[0]/transform_size);
						int py = (int) Math.round(pXpYD[1]/transform_size);
						int spx = (int) Math.round(pXpYD[0]*scale);
						int spy = (int) Math.round(pXpYD[1]*scale);
						if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
							if (!(pXpYD[2] < dbg_img[indx_d1][px + py* tilesX])) {
								dbg_img[indx_s1][px + py* tilesX] = dbg_img[indx_s0][nTile];
								dbg_img[indx_d1][px + py* tilesX] = pXpYD[2];
								if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
									ds[1][spx + spy* stilesX] = reduce*dbg_img[indx_s0][nTile];
									ds[0][spx + spy* stilesX] = reduce*pXpYD[2];
								}								
							}
						}
					}
				}
			}
		}
		ds[0] = (new DoubleGaussianBlur()).blurWithNaN(
				ds[0], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		ds[1] = (new DoubleGaussianBlur()).blurWithNaN(
				ds[1], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		ds0[0] = (new DoubleGaussianBlur()).blurWithNaN(
				ds0[0], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		ds0[1] = (new DoubleGaussianBlur()).blurWithNaN(
				ds0[1], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		dbg_img[indx_s1] = new double[tiles];
		dbg_img[indx_d1] = new double[tiles];
		dbg_img[indx_s0] = new double[tiles];
		dbg_img[indx_d0] = new double[tiles];
		for (int stileY = 0; stileY < stilesY; stileY++) {
			int tileY = stileY / iscale; 
			for (int stileX = 0; stileX < stilesX; stileX++) {
				int tileX = stileX / iscale;
				int stile = stileX + stileY * stilesX;
				int tile =  tileX +  tileY *  tilesX;
				dbg_img[indx_s1][tile] += ds[1][stile];
				dbg_img[indx_d1][tile] += ds[0][stile];
				dbg_img[indx_s0][tile] += ds0[1][stile];
				dbg_img[indx_d0][tile] += ds0[0][stile];
			}
		}
		
		// restore original
		ersCorrection.setupERS(
				saved_ers_wxyz_center,      // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				saved_ers_wxyz_center_dt,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				saved_ers_wxyz_center_d2t,  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				saved_ers_watr_center_dt,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				saved_ers_watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		return dbg_img;
		
		
	}
	
	@Deprecated	
	public double[][]  test_back_forth0(
			double [] camera_xyz,
			double [] camera_atr)
	{
		camera_atr[0] = 1.0;
		
		double [] camera_xyz0 = {0.0,0.0,0.0};
		double [] camera_atr0 = {0.0,0.0,0.0};
		double [] camera_xyz1 = {0.0,0.0,0.0};
		double [] camera_atr1 = {0.0,0.0,0.0};

		double    line_err = 10.0; // 0.1; // bug
		int tilesX = tp.getTilesX();
		int tilesY = tp.getTilesY();
		int tiles = tilesX*tilesY;
		int transform_size = tp.getTileSize();
		String [] titles = {"s0", "s1", "d0", "d1", "pX0","pX", "pY0","pY", "x","y","z","w"};
//		double [] ers_watr_center_dt =  {0.168792, 0.037145, -0.060279};
		double [] zero3 =               {0.0,0.0,0.0};	
//		double [] watr_center_dt =  {0.168792, 0.037145, -0.060279};
		double [] watr_center_dt =  {0.0, 0.0, 0.0};
		double [] watr_center_d2t = {0.0, 0.0, 0.0};
		double [] wxyz_center =     {0.0, 0.0, 0.0};	
		double [] wxyz_center_dt =  {0.0, 0.0, 0.0};	
		double [] wxyz_center_d2t = {0.0, 0.0, 0.0};	
		boolean force_manual = false; // true;
		
		int indx_s0 =  0;
		int indx_s1 =  1;
		int indx_d0 =  2;
		int indx_d1 =  3;
		int indx_pX0 = 4;
		int indx_pX =  5;
		int indx_pY0 = 6;
		int indx_pY =  7;
		int indx_xyz = 8;
		double [][] dbg_img = new double [titles.length][tiles];
		for (int i = 0; i <dbg_img.length; i++) {
			for (int j = 0; j <dbg_img[i].length; j++) {
				dbg_img[i][j] = Double.NaN;
			}
		}
		final int iscale = 8;
		int stilesX = iscale*tilesX; 
		int stilesY = iscale*tilesY;
		int stiles = stilesX*stilesY;
		double sigma = 0.5 * iscale;
		double scale = 1.0*iscale/transform_size;
		double reduce = 1.0/iscale/iscale;
		double [][] ds = new double [2][stiles];
		double [][] ds0 = new double [2][stiles];
		for (int i = 0; i <ds.length; i++) {
			for (int j = 0; j <ds[i].length; j++) {
				ds[i][j] = Double.NaN;
				ds0[i][j] = Double.NaN;
			}
		}

		dbg_img [indx_d0] = this.dsi[TwoQuadCLT.DSI_DISPARITY_MAIN];
		dbg_img [indx_s0] = this.dsi[TwoQuadCLT.DSI_STRENGTH_MAIN];
		ErsCorrection ersCorrection = (ErsCorrection) geometryCorrection;
		
		// save original
		double [] saved_ers_wxyz_center =     ersCorrection.ers_wxyz_center;
		double [] saved_ers_wxyz_center_dt =  ersCorrection.ers_wxyz_center_dt;
		double [] saved_ers_wxyz_center_d2t = ersCorrection.ers_wxyz_center_d2t;
		double [] saved_ers_watr_center_dt =  ersCorrection.ers_watr_center_dt;
		double [] saved_ers_watr_center_d2t = ersCorrection.ers_watr_center_d2t;
		System.out.println("\n"+image_name+":\n"+ersCorrection.extrinsic_corr.toString());
		System.out.println(String.format("%s: ers_wxyz_center=     %f, %f, %f", image_name,saved_ers_wxyz_center[0], saved_ers_wxyz_center[1],saved_ers_wxyz_center[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_dt=  %f, %f, %f", image_name,saved_ers_wxyz_center_dt[0], saved_ers_wxyz_center_dt[1],saved_ers_wxyz_center_dt[2] ));
		System.out.println(String.format("%s: ers_wxyz_center_d2t= %f, %f, %f", image_name,saved_ers_wxyz_center_d2t[0], saved_ers_wxyz_center_d2t[1],saved_ers_wxyz_center_d2t[2] ));
		System.out.println(String.format("%s: ers_watr_center_dt=  %f, %f, %f", image_name,saved_ers_watr_center_dt[0], saved_ers_watr_center_dt[1],saved_ers_watr_center_dt[2] ));
		System.out.println(String.format("%s: ers_watr_center_d2t= %f, %f, %f", image_name,saved_ers_watr_center_d2t[0], saved_ers_watr_center_d2t[1],saved_ers_watr_center_d2t[2] ));

		// Comment out to use manual values
		
		if (force_manual) {
			System.out.println("\nUsing manually set values:");
			System.out.println(String.format("%s: wxyz_center=     %f, %f, %f",     image_name,wxyz_center[0], wxyz_center[1],wxyz_center[2] ));
			System.out.println(String.format("%s: wxyz_center_dt=  %f, %f, %f",  image_name,wxyz_center_dt[0], wxyz_center_dt[1],wxyz_center_dt[2] ));
			System.out.println(String.format("%s: wxyz_center_d2t= %f, %f, %f", image_name,wxyz_center_d2t[0], wxyz_center_d2t[1],wxyz_center_d2t[2] ));
			System.out.println(String.format("%s: watr_center_dt=  %f, %f, %f",  image_name,watr_center_dt[0], watr_center_dt[1],watr_center_dt[2] ));
			System.out.println(String.format("%s: watr_center_d2t= %f, %f, %f", image_name,watr_center_d2t[0], watr_center_d2t[1],watr_center_d2t[2] ));
		}else {
			wxyz_center =     saved_ers_wxyz_center;
			wxyz_center_dt =  saved_ers_wxyz_center_dt;
			wxyz_center_d2t = saved_ers_wxyz_center_d2t;  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
			watr_center_dt =  saved_ers_watr_center_dt;
			watr_center_d2t = saved_ers_watr_center_d2t;
		}
		
		ersCorrection.setupERS(
				zero3,            // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				wxyz_center_dt,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				wxyz_center_d2t,  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				watr_center_dt,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		
		for (int tileY = 0; tileY < tilesY; tileY++) {
			int stileY = iscale * tileY + iscale/2;
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int stileX = iscale * tileX + iscale/2;
				int stile = stileX + stilesX * stileY;
				int nTile = tileX + tileY * tilesX;
				double centerX = tileX * transform_size + transform_size/2; //  - shiftX;
				double centerY = tileY * transform_size + transform_size/2; //  - shiftY;
				ds0[0][stile] = reduce* dbg_img [indx_d0][nTile]; //  this.dsi[TwoQuadCLT.DSI_DISPARITY_MAIN][nTile];
				ds0[1][stile] = reduce* dbg_img [indx_s0][nTile]; // this.dsi[TwoQuadCLT.DSI_STRENGTH_MAIN][nTile];
				
				dbg_img[indx_pX0][nTile] = centerX;
				dbg_img[indx_pY0][nTile] = centerY;
				double [] wxyz = ersCorrection.getWorldCoordinatesERS(
						centerX,                                        // double px,
						centerY,                                        // double py,
						this.dsi[TwoQuadCLT.DSI_DISPARITY_MAIN][nTile], // double disparity,
						true,                                           // boolean correctDistortions,// correct distortion (will need corrected background too !)
						camera_xyz0,  // double [] camera_xyz, // camera center in world coordinates
						camera_atr0); // double [] camera_atr)  // camera orientation relative to world frame
				if ((wxyz != null) && (wxyz[2] < 0.0)) {
					for (int i = 0; i < 4; i ++) {
						dbg_img[indx_xyz+i][nTile] = wxyz[i]; 
					}
				}
			}
		}
		
		ersCorrection.setupERS(
		zero3,   // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
		zero3,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
		zero3,   // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
		zero3,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		zero3);  // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		
		for (int tileY = 0; tileY < tilesY; tileY++) {
			for (int tileX = 0; tileX < tilesX; tileX++) {
				int nTile = tileX + tileY * tilesX;
				if (!Double.isNaN(dbg_img[indx_xyz][nTile])) {
					double [] wxyz = {
							dbg_img[indx_xyz][nTile],
							dbg_img[indx_xyz+1][nTile],
							dbg_img[indx_xyz+2][nTile],
							dbg_img[indx_xyz+3][nTile]};
					double [] pXpYD = ersCorrection.getImageCoordinatesERS(
							wxyz, // double [] wxyz, 
							true,                                           // boolean correctDistortions,// correct distortion (will need corrected background too !)
							camera_xyz1,  // double [] camera_xyz, // camera center in world coordinates
							camera_atr1,  // double [] camera_atr)  // camera orientation relative to world frame
							line_err);   // double    line_err) // threshold error in scan lines (1.0)
					if (pXpYD != null) {
						dbg_img[indx_pX][nTile] = pXpYD[0];
						dbg_img[indx_pY][nTile] = pXpYD[1];
						
						int px = (int) Math.round(pXpYD[0]/transform_size);
						int py = (int) Math.round(pXpYD[1]/transform_size);
						int spx = (int) Math.round(pXpYD[0]*scale);
						int spy = (int) Math.round(pXpYD[1]*scale);
						if ((px >= 0) && (py >= 0) && (px < tilesX) & (py < tilesY)) {
							if (!(pXpYD[2] < dbg_img[indx_d1][px + py* tilesX])) {
								dbg_img[indx_s1][px + py* tilesX] = dbg_img[indx_s0][nTile];
								dbg_img[indx_d1][px + py* tilesX] = pXpYD[2];
								if ((spx >= 0) && (spy >= 0) && (spx < stilesX) & (spy < stilesY)) {
									ds[1][spx + spy* stilesX] = reduce*dbg_img[indx_s0][nTile];
									ds[0][spx + spy* stilesX] = reduce*pXpYD[2];
								}								
							}
						}
					}
				}
			}
		}
		ds[0] = (new DoubleGaussianBlur()).blurWithNaN(
				ds[0], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		ds[1] = (new DoubleGaussianBlur()).blurWithNaN(
				ds[1], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		ds0[0] = (new DoubleGaussianBlur()).blurWithNaN(
				ds0[0], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		ds0[1] = (new DoubleGaussianBlur()).blurWithNaN(
				ds0[1], // double[] pixels,
				null,  // double [] in_weight, // or null
				stilesX, // int width,
				stilesY, // int height,
				sigma, // double sigmaX,
				sigma, // double sigmaY,
				0.01); // double accuracy);
		dbg_img[indx_s1] = new double[tiles];
		dbg_img[indx_d1] = new double[tiles];
		dbg_img[indx_s0] = new double[tiles];
		dbg_img[indx_d0] = new double[tiles];
		for (int stileY = 0; stileY < stilesY; stileY++) {
			int tileY = stileY / iscale; 
			for (int stileX = 0; stileX < stilesX; stileX++) {
				int tileX = stileX / iscale;
				int stile = stileX + stileY * stilesX;
				int tile =  tileX +  tileY *  tilesX;
				dbg_img[indx_s1][tile] += ds[1][stile];
				dbg_img[indx_d1][tile] += ds[0][stile];
				dbg_img[indx_s0][tile] += ds0[1][stile];
				dbg_img[indx_d0][tile] += ds0[0][stile];
			}
		}
		
		String title = String.format("%s_%f:%f:%f_%f:%f:%f",image_name,
				camera_xyz0[0],camera_xyz0[1],camera_xyz0[2],camera_atr0[0],camera_atr0[1],camera_atr0[2]);
		
		
		title += String.format ("_ers_%f:%f:%f_%f:%f:%f",wxyz_center_dt[0],wxyz_center_dt[1],wxyz_center_dt[2],
				watr_center_dt[0],watr_center_dt[1],watr_center_dt[2]);
		(new ShowDoubleFloatArrays()).showArrays(
				dbg_img,
				 tilesX,
				 tilesY,
				 true,
				 title,
				 titles);
		// restore original
		ersCorrection.setupERS(
				saved_ers_wxyz_center,      // double [] wxyz_center,     // world camera XYZ (meters) for the frame center
				saved_ers_wxyz_center_dt,   // double [] wxyz_center_dt,  // world camera Vx, Vy, Vz (m/s)
				saved_ers_wxyz_center_d2t,  // double [] wxyz_center_d2t, // world camera Vx, Vy, Vz (m/s^2)
				saved_ers_watr_center_dt,   // double [] watr_center_dt,  // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
				saved_ers_watr_center_d2t); // double [] watr_center_d2t); // camera rotaions (az, tilt, roll in radians/s, corresponding to the frame center)
		
		return dbg_img;
		
		
	}
	
	
	public void setGPU(GpuQuad gpuQuad) {
		this.gpuQuad = gpuQuad;
	}
	public GpuQuad getGPU() {
		return this.gpuQuad;
	}
	
	
	public void genSave4sliceImage(
			CLTParameters                                   clt_parameters,
			String                                          suffix,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			CorrectionColorProc.ColorGainsParameters        channelGainParameters,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			final int                                       threadsMax,  // maximal number of threads to launch
			final int                                       debugLevel){
		String x3d_path = getX3dDirectory();
		String file_name = image_name + suffix;
		String file_path = x3d_path + Prefs.getFileSeparator() + file_name + ".tiff";
		if ((getGPU() != null) && (getGPU().getQuadCLT() != this)) {
			getGPU().updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
		}

		ImagePlus img_noise = processCLTQuadCorrGPU(
				null,                  // ImagePlus []                        imp_quad, //null will be OK
				null,                  // boolean [][]                                    saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
				clt_parameters,        // CLTParameters                                   clt_parameters,
				debayerParameters,     // EyesisCorrectionParameters.DebayerParameters    debayerParameters,
				colorProcParameters,   // ColorProcParameters                             colorProcParameters,
				channelGainParameters, // CorrectionColorProc.ColorGainsParameters        channelGainParameters,
				rgbParameters,         // EyesisCorrectionParameters.RGBParameters        rgbParameters,
				null,                  // double []	                                      scaleExposures, // probably not needed here - restores brightness of the final image
				true,                  // boolean                                         only4slice,
				threadsMax,            // final int                                       threadsMax,  // maximal number of threads to launch
				false,                 // final boolean                                   updateStatus,
				debugLevel);           // final int                                       debugLevel);
//		FileSaver fs=new FileSaver(img_noise); // is null, will be saved inside to /home/elphel/lwir16-proc/proc1/results_cuda/1626032208_613623-AUX-SHIFTED-D0.0
//		fs.saveAsTiff(file_path);
	}	
	
	
	// added 12/21
	public void processCLTQuadCorrs( // not used in lwir
			CLTParameters           clt_parameters,
			EyesisCorrectionParameters.DebayerParameters     debayerParameters,
			ColorProcParameters colorProcParameters,
			CorrectionColorProc.ColorGainsParameters     channelGainParameters,
			EyesisCorrectionParameters.RGBParameters             rgbParameters,
			//			  int          convolveFFTSize, // 128 - fft size, kernel size should be size/2
			final boolean    apply_corr,    // calculate and apply additional fine geometry correction
			final boolean    infinity_corr, // calculate and apply geometry correction at infinity
			final int        threadsMax,    // maximal number of threads to launch
			final boolean    updateStatus,
			final int        debugLevel)
	{
		if (infinity_corr && (clt_parameters.z_correction != 0.0)){
			System.out.println(
					"****************************************\n"+
							"* Resetting manual infinity correction *\n"+
					"****************************************\n");
			clt_parameters.z_correction = 0.0;
		}

		this.startTime=System.nanoTime();
		String [] sourceFiles=correctionsParameters.getSourcePaths();
		SetChannels [] set_channels=setChannels(debugLevel);
		if ((set_channels == null) || (set_channels.length==0)) {
			System.out.println("No files to process (of "+sourceFiles.length+")");
			return;
		}
		// multiply each image by this and divide by individual (if not NaN)
		double [] referenceExposures = null;
		//		  if (!colorProcParameters.lwir_islwir) {
		if (!isLwir()) {
			referenceExposures=eyesisCorrections.calcReferenceExposures(debugLevel);
		}
		for (int nSet = 0; nSet < set_channels.length; nSet++){
			int [] channelFiles = set_channels[nSet].fileNumber();
			boolean [][] saturation_imp = (clt_parameters.sat_level > 0.0)? new boolean[channelFiles.length][] : null;
			double [] scaleExposures = new double[channelFiles.length];

			//			  ImagePlus [] imp_srcs = 
			conditionImageSet(
					clt_parameters,             // EyesisCorrectionParameters.CLTParameters  clt_parameters,
					colorProcParameters,
					sourceFiles,                // String []                                 sourceFiles,
					set_channels[nSet].name(),  // String                                    set_name,
					referenceExposures,         // double []                                 referenceExposures,
					channelFiles,               // int []                                    channelFiles,
					scaleExposures,   //output  // double [] scaleExposures
					saturation_imp,   //output  // boolean [][]                              saturation_imp,
					threadsMax,                 // int                                       threadsMax,
					debugLevel); // int                                       debugLevel);


			// once per quad here
			processCLTQuadCorrGPU(  // returns ImagePlus, but it already should be saved/shown
					null,           // imp_srcs, // [srcChannel], // should have properties "name"(base for saving results), "channel","path"
					saturation_imp, // boolean [][] saturation_imp, // (near) saturated pixels or null
					clt_parameters,
					debayerParameters,
					colorProcParameters,
					channelGainParameters,
					rgbParameters,
					scaleExposures,
					false, // true,            // boolean                                         only4slice,
					threadsMax,                // final int                                       threadsMax,  // maximal number of threads to launch
					false,                     // final boolean                                   updateStatus,
					debugLevel);               // final int                                       debugLevel);
			if (debugLevel >-1) System.out.println("Processing set "+(nSet+1)+" (of "+set_channels.length+") finished at "+
					IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
			if (eyesisCorrections.stopRequested.get()>0) {
				System.out.println("User requested stop");
				System.out.println("Processing "+(nSet + 1)+" file sets (of "+set_channels.length+") finished at "+
						IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
				return;
			}
		}
		System.out.println("processCLTQuadCorrs(): processing "+getTotalFiles(set_channels)+" files ("+set_channels.length+" file sets) finished at "+
				IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec, --- Free memory="+Runtime.getRuntime().freeMemory()+" (of "+Runtime.getRuntime().totalMemory()+")");
	}
	
	
	
	public ImagePlus processCLTQuadCorrGPU(
			ImagePlus []                                    imp_quad, //null will be OK
			boolean [][]                                    saturation_imp, // (near) saturated pixels or null // Not needed use this.saturation_imp
			CLTParameters                                   clt_parameters,
			EyesisCorrectionParameters.DebayerParameters    debayerParameters,
			ColorProcParameters                             colorProcParameters,
			CorrectionColorProc.ColorGainsParameters        channelGainParameters,
			EyesisCorrectionParameters.RGBParameters        rgbParameters,
			double []	                                    scaleExposures, // probably not needed here - restores brightness of the final image
			boolean                                         only4slice,
			final int                                       threadsMax,  // maximal number of threads to launch
			final boolean                                   updateStatus,
			final int                                       debugLevel){
		if (gpuQuad == null) {
			System.out.println("GPU instance is not initialized, using CPU mode");
			processCLTQuadCorrCPU(
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
			
			return null;
		}
		

// GPU-specific
		boolean showCoord = debugLevel > 1;
		boolean is_mono = isMonochrome();
		boolean is_lwir = isLwir();
		final boolean      batch_mode = clt_parameters.batch_run; //disable any debug images
		boolean test_execCorr2D =    batch_mode?false: false; // true;
		boolean try_lores =          batch_mode?false: true;
		boolean show_textures_rgba = batch_mode?false: clt_parameters.show_rgba_color;
		boolean try_textures =       batch_mode?false: true;
		
		ImageDtt image_dtt = new ImageDtt(
				getNumSensors(),
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				isAux(),
				is_mono,
				is_lwir,
				clt_parameters.getScaleStrength(isAux())); // 1.0);
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		float [][] lpf_rgb = new float[][] {
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_r),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_b),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_g),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_m)
		};
		gpuQuad.setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				true); // !batch_mode);

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
		
		ImagePlus [] results = null;
		if (imp_quad !=  null) {
			results = new ImagePlus[imp_quad.length];
			for (int i = 0; i < results.length; i++) {
				results[i] = imp_quad[i];
				results[i].setTitle(results[i].getTitle()+"RAW");
			}
		}
		if (debugLevel>1) System.out.println("processing: "+gpuQuad.quadCLT);

		double z_correction =  clt_parameters.z_correction;
		if (clt_parameters.z_corr_map.containsKey(image_name)){
			z_correction +=clt_parameters.z_corr_map.get(image_name);// not used in lwir
		}
		final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		// the following will be replaced by setFullFrameImages() with target disparities
		TpTask [] tp_tasks  = gpuQuad.setFullFrameImages( // when disparities array is known - use different arguments of setFullFrameImages
				null,                                 // Rectangle                 woi,
				clt_parameters.gpu_woi_round,         // boolean                   round_woi,
	    		(float) clt_parameters.disparity,     // float                     target_disparity, // apply same disparity to all tiles
	    		(float) disparity_corr,
	    		0xf,                                  // int                       out_image, // from which tiles to generate image (currently 0/1)
	    		0x3f,                                 // int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
	    		debugLevel);                          // final int                 debugLevel) - not yet used
		if (tp_tasks.length == 0) {
			System.out.println("Empty tasks - nothing to do");
			return null;
		}

		gpuQuad.setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_aux, // boolean use_aux)
				clt_parameters.img_dtt.gpu_verify); // boolean verify

		gpuQuad.execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers

		if (showCoord) {
			double max_r_by_rdist_err = gpuQuad.maxRbyRDistErr();
			System.out.println("Maximal R-by-RDist error = "+max_r_by_rdist_err);
			showCentersPortsXY("Tile_Virtual_Centers");
//			showCentersXY("Tile_Virtual_Centers");
//			showPortsXY("Tile_Physical_Sensor_Centers");
		}
		
		
		gpuQuad.execConvertDirect();
		int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt, getNumSensors());

		if (test_execCorr2D) {
			int [] i_mcorr_sel = Correlation2d.intCorrPairs(
					mcorr_sel,
					getNumSensors(),
					4); // int num_out); should be 4 int
			gpuQuad.setCorrMask(i_mcorr_sel);
			
			double fat_zero = 30.0; // 2000.0; // 30.00;
			int gpu_corr_rad =  7;
			int corr_size = 2* gpu_corr_rad + 1;

			final int numcol = isMonochrome()?1:3;
			final double [] col_weights= new double [numcol]; // colors are RBG
			if (isMonochrome()) {
				col_weights[0] = 1.00; // Was 0 ! 0;
			} else {
				col_weights[2] = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
				col_weights[0] = clt_parameters.corr_red *  col_weights[2];
				col_weights[1] = clt_parameters.corr_blue * col_weights[2];
			}

			//Generate 2D phase correlations from the CLT representation
			// Try to zero out memory before calculating?	
			gpuQuad.eraseGpuCorrs(); // no-diff
			gpuQuad.execCorr2D(
					gpuQuad.getCorrMask(),    // boolean [] pair_select,
					col_weights,              // double [] scales,
					fat_zero,                 // double fat_zero);
					gpu_corr_rad);            // int corr_radius
			// Should be done before execCorr2D_TD as corr_indices are shared to save memory
			int [] corr_indices = gpuQuad.getCorrIndices();
			// the following is not yet shared
			float [][] corr2D =  gpuQuad.getCorr2D(
					gpu_corr_rad); //  int corr_rad);
			final int tilesX=gpuQuad.getTilesX(); // width/transform_size;
			final int tilesY=gpuQuad.getTilesY(); // final int tilesY=height/transform_size;
//			int num_tiles = tilesX * tilesY;
			int sq = 16;
			int num_pairs = gpuQuad.getNumUsedPairs();
			float [][] corr_img = new float [num_pairs][tilesY * sq * tilesX * sq];
			for (int pair = 0; pair < num_pairs; pair++) {
				Arrays.fill(corr_img[pair], Float.NaN);
			}
			for (int ict = 0; ict < corr_indices.length; ict++){
				//    	int ct = cpu_corr_indices[ict];
				int ctt = ( corr_indices[ict] >>  GPUTileProcessor.CORR_NTILE_SHIFT);
				int cpair = corr_indices[ict] & ((1 << GPUTileProcessor.CORR_NTILE_SHIFT) - 1);
				int ty = ctt / tilesX;
				int tx = ctt % tilesX;
				int dst_offs0 = (ty * sq * tilesX * sq) + (tx * sq);
				for (int iy = 0; iy < corr_size; iy++){
					int dst_offs = dst_offs0 + iy * (tilesX * sq);
					System.arraycopy(corr2D[ict], iy * corr_size, corr_img[cpair], dst_offs, corr_size);
				}
			}
			(new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					corr_img,
					tilesX * sq,
					tilesY * sq,
					true,
					"test-corr1");
		}
		
		
		gpuQuad.execImcltRbgAll(is_mono); 

		// get data back from GPU
//		float [][][] iclt_fimg = new float [gpuQuad.getNumCams()][][];
		float [][][] iclt_fimg = new float [getNumSensors()][][];

		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			iclt_fimg[ncam] = gpuQuad.getRBG(ncam);
		}

		int out_width =  gpuQuad.getImageWidth()  + gpuQuad.getDttSize();
		int out_height = gpuQuad.getImageHeight() + gpuQuad.getDttSize();
		if (isLwir() && colorProcParameters.lwir_autorange) {
			double rel_low =  colorProcParameters.lwir_low;
			double rel_high = colorProcParameters.lwir_high;
			if (!Double.isNaN(getLwirOffset())) {
				rel_low -=  getLwirOffset();
				rel_high -= getLwirOffset();
			}
			double [] cold_hot =  autorange(
					iclt_fimg,                         // iclt_data, // double [][][] iclt_data, //  [iQuad][ncol][i] - normally only [][2][] is non-null
					rel_low,                           // double hard_cold,// matches data, DC (this.lwir_offset)  subtracted
					rel_high,                          // double hard_hot,   // matches data, DC (this.lwir_offset)  subtracted
					colorProcParameters.lwir_too_cold, // double too_cold, // pixels per image
					colorProcParameters.lwir_too_hot,  // double too_hot,  // pixels per image
					1024); // int num_bins)
			if (cold_hot != null) {
				if (!Double.isNaN(getLwirOffset())) {
					cold_hot[0] += getLwirOffset();
					cold_hot[1] += getLwirOffset();
				}
			}
			setColdHot(cold_hot); // will be used for shifted images and for texture tiles
		}
		
		/* Prepare 4-channel images*/
		ImagePlus [] imps_RGB = new ImagePlus[iclt_fimg.length];
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
            String title=String.format("%s%s-%02d",image_name, sAux(), ncam);
			imps_RGB[ncam] = linearStackToColor( // probably no need to separate and process the second half with quadCLT_aux (!)
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
		
		if (clt_parameters.gen_chn_img || only4slice) { // save and show 4-slice image
			// combine to a sliced color image
			int [] slice_seq = {0,1,3,2}; //clockwise
			if (imps_RGB.length > 4) {
				slice_seq = new int [imps_RGB.length];
				for (int i = 0; i < slice_seq.length; i++) {
					slice_seq[i] = i;
				}
			}
			int width = imps_RGB[0].getWidth();
			int height = imps_RGB[0].getHeight();
			ImageStack array_stack=new ImageStack(width,height);
			for (int i = 0; i<slice_seq.length; i++){
				///				if (imps_RGB[slice_seq[i]] != null) {
				array_stack.addSlice("port_"+slice_seq[i], imps_RGB[slice_seq[i]].getProcessor().getPixels());
				///				} else {
				///					array_stack.addSlice("port_"+slice_seq[i], results[slice_seq[i]].getProcessor().getPixels());
				///				}
			}
			ImagePlus imp_stack = new ImagePlus(image_name+sAux()+"-SHIFTED-D"+clt_parameters.disparity, array_stack);
			imp_stack.getProcessor().resetMinAndMax();
			if (only4slice) {
				return imp_stack;
			} else {
				if (!batch_mode) {
					imp_stack.updateAndDraw();
				}
				eyesisCorrections.saveAndShowEnable(
						imp_stack,  // ImagePlus             imp,
						correctionsParameters, // EyesisCorrectionParameters.CorrectionParameters  correctionsParameters,
						true, // boolean               enableSave,
						!batch_mode) ;// boolean               enableShow);
			}
		}

		if (clt_parameters.gen_4_img) { // save 4 JPEG images
			// Save as individual JPEG images in the model directory
			String x3d_path = getX3dDirectory();
			for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
				EyesisCorrections.saveAndShow(
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
// try textures here
		if (show_textures_rgba) {
			float [][] texture_img = new float [isMonochrome()?2:4][];
//			double [] col_weights = new double[3];
			double [] col_weights = new double[isMonochrome() ? 1 :3];
			if (isMonochrome()) {
				col_weights[0] = 1.0;
//			if (isMonochrome()) {
//				col_weights[0] = 1.0;
//				col_weights[1] = 0.0;
//				col_weights[2] = 0.0;// green color/mono
			} else {
				col_weights[2] = 1.0/(1.0 +  clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
				col_weights[0] = clt_parameters.corr_red *  col_weights[2];
				col_weights[1] = clt_parameters.corr_blue * col_weights[2];
			}
			Rectangle woi = new Rectangle(); // will be filled out to match actual available image
			gpuQuad.execRBGA(
					col_weights,                   // double [] color_weights,
					isLwir(),                      // boolean   is_lwir,
					clt_parameters.min_shot,       // double    min_shot,           // 10.0
					clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
					clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
					clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					clt_parameters.dust_remove);   // boolean   dust_remove,
			float [][] rbga = gpuQuad.getRBGA(
					(isMonochrome() ? 1 : 3), // int     num_colors,
					woi);
			for (int ncol = 0; ncol < texture_img.length; ncol++) if (ncol < rbga.length) {
				texture_img[ncol] = rbga[ncol];
			}
			// first try just multi-layer, then with palette
			
			(new ShowDoubleFloatArrays()).showArrays( // show slices RBGA (colors - 256, A - 1.0)
					rbga,
					woi.width,
					woi.height,
					true,
					getImageName()+"-RGBA-STACK-D"+clt_parameters.disparity+
					":"+clt_parameters.gpu_woi_tx+":"+clt_parameters.gpu_woi_ty+
					":"+clt_parameters.gpu_woi_twidth+":"+clt_parameters.gpu_woi_theight+
					":"+(clt_parameters.gpu_woi_round?"C":"R")
					//,new String[] {"R","B","G","A"}
					);
			float [][] rgb_text; //  = new float [(isMonochrome() ? 1 : 3)][]; // {rbga[0],rbga[1],rbga[2]};
			float [][] rgba_text; //  = new float [(isMonochrome() ? 2 : 4)][]; // {rbga[0],rbga[1],rbga[2],rbga[3]};
			if (isMonochrome()) {
				rgb_text = new float [][]{rbga[0]};
				rgba_text = new float [][]{rbga[0]};
			} else {
				rgb_text = new float [][]{rbga[0],rbga[2],rbga[1]};
				rgba_text = new float [][]{rbga[0],rbga[2],rbga[1],rbga[3]};
			}
			//isAux()
			ImagePlus imp_rgba = linearStackToColor(
					clt_parameters,
					colorProcParameters,
					rgbParameters,
					getImageName()+"-texture", // String name,
					"-D"+clt_parameters.disparity+"-"+(isAux()?"AUX":"MAIN")+ "GPU", //String suffix, // such as disparity=...
					toRGB,
					!correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
					false, // true, // boolean saveShowIntermediate, // save/show if set globally
					false, // true, // boolean saveShowFinal,        // save/show result (color image?)
					((clt_parameters.alpha1 > 0)? rgba_text: rgb_text),
					woi.width, // clt_parameters.gpu_woi_twidth * image_dtt.transform_size, //  tilesX *  image_dtt.transform_size,
					woi.height, // clt_parameters.gpu_woi_theight *image_dtt.transform_size, //  tilesY *  image_dtt.transform_size,
					1.0,         // double scaleExposure, // is it needed?
					debugLevel );
			imp_rgba.getProcessor().resetMinAndMax();
			imp_rgba.show();
		}
		if (try_lores) {
			//Generate non-overlapping (16x16) texture tiles, prepare
//			double [] col_weights = new double[3];
			double [] col_weights = new double[isMonochrome() ? 1 :3];
			if (isMonochrome()) {
				col_weights[0] = 1.0;
			} else {
				col_weights[2] = 1.0/(1.0 +  clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
				col_weights[0] = clt_parameters.corr_red *  col_weights[2];
				col_weights[1] = clt_parameters.corr_blue * col_weights[2];
			}
			gpuQuad.execTextures(
				col_weights,                   // double [] color_weights,
				isLwir(),                      // boolean   is_lwir,
				clt_parameters.min_shot,       // double    min_shot,           // 10.0
				clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
				clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
				clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
				clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // boolean   dust_remove,
				false,                         // boolean   calc_textures,
				true,                          // boolean   calc_extra
				false);                        // boolean   linescan_order) // TODO: use true to avoid reordering of the low-res output 
			float [][] extra = gpuQuad.getExtra(); // now 4*numSensors
			int num_cams = getNumSensors();
			
			(new ShowDoubleFloatArrays()).showArrays( // show slices RBGA (colors - 256, A - 1.0)
					extra,
					gpuQuad.img_width / GPUTileProcessor.DTT_SIZE,
					gpuQuad.img_height / GPUTileProcessor.DTT_SIZE,
					true,
					getImageName()+"-LOW-RES"
					//,new String[] {"R","B","G","A"}
					);
			for (int nc = 00; nc < (extra.length - num_cams); nc++) {
				int sindx = nc + num_cams;
				
			}
			if (try_textures) {
				//Generate non-overlapping (16x16) texture tiles, prepare 
				gpuQuad.execTextures(
					col_weights,                   // double [] color_weights,
					isLwir(),         // boolean   is_lwir,
					clt_parameters.min_shot,       // double    min_shot,           // 10.0
					clt_parameters.scale_shot,     // double    scale_shot,         // 3.0
					clt_parameters.diff_sigma,     // double    diff_sigma,         // pixel value/pixel change Used much larger sigma = 10.0 instead of 1.5
					clt_parameters.diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					clt_parameters.min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					clt_parameters.dust_remove,    // boolean   dust_remove,
					true,                          // boolean   calc_textures,
					false,                         // boolean   calc_extra
					false);                        // boolean   linescan_order) 

				int [] texture_indices = gpuQuad.getTextureIndices();
				int numcol = isMonochrome()? 1 : 3;
				int    num_src_slices = numcol + 1; //  + (clt_parameters.keep_weights?(ports + numcol + 1):0); // 12 ; // calculate
				float [] flat_textures =  gpuQuad.getFlatTextures( // fatal error has been detected by the Java Runtime Environment:
						texture_indices.length,
						numcol, // int     num_colors,
			    		false); // clt_parameters.keep_weights); // boolean keep_weights);
				int tilesX = gpuQuad.img_width / GPUTileProcessor.DTT_SIZE;
				int tilesY = gpuQuad.img_height / GPUTileProcessor.DTT_SIZE;
				double [][][][] texture_tiles = new double [tilesY][tilesX][][];
				gpuQuad.doubleTextures(
			    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
			    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
			    		texture_indices,                     // int []       indices,
			    		flat_textures,                       // float [][][] ftextures,
			    		tilesX,                              // int          full_width,
			    		isMonochrome()? 2: 4,                                   // rbga only /int          num_slices
			    		num_src_slices                       // int          num_src_slices
			    		);
				int num_out_slices = 0;
				for (int nt = 0; nt < tilesY * tilesX; nt++) {
					if (texture_tiles[nt / tilesX][nt % tilesX] != null) {
						num_out_slices = texture_tiles[nt / tilesX][nt % tilesX].length;
						break;
					}
				}
				if (num_out_slices > 0) {
					int ssize = 2*GPUTileProcessor.DTT_SIZE;
					int width =  tilesX * ssize;
					int height = tilesY * ssize;
					double [][] dbg_nonoverlap = new double[num_out_slices][width * height];
					for (int slice = 0; slice < num_out_slices; slice++) {
						Arrays.fill(dbg_nonoverlap[slice], Double.NaN);
					}
					for (int ty = 0; ty < tilesY; ty++) {
						for (int tx = 0; tx < tilesX; tx++) {
							if (texture_tiles[ty][tx] != null) {
								for (int slice = 0; slice < num_out_slices; slice++) {
									for (int row = 0; row < ssize; row++) {
										System.arraycopy(
												texture_tiles[ty][tx][slice],
												row * ssize,
												dbg_nonoverlap[slice],
												(ty * ssize + row) * width + (tx * ssize),
												ssize);
									}
								}
							}
						}
						
					}
					(new ShowDoubleFloatArrays()).showArrays( // show slices RBGA (colors - 256, A - 1.0)
							dbg_nonoverlap,
							width,
							height,
							true,
							getImageName()+"-textures"
							);
				}
				System.out.println("try_textures DONE");
			}
		}
// try low-res and non-overlap textures
		return null;	  
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
			if (quadCLT_main.getGPU() != null) quadCLT_main.gpuResetCorrVector(); // getGPU().resetGeometryCorrectionVector();
			if (quadCLT_aux.getGPU() != null)  quadCLT_aux.gpuResetCorrVector();   // .getGPU().resetGeometryCorrectionVector();
		}
		
// get fat_zero (absolute) and color scales
		boolean is_mono = quadCLT_main.isMonochrome();
		boolean is_lwir = quadCLT_main.isLwir();
		boolean is_aux =  quadCLT_main.isAux();


		double    fat_zero = clt_parameters.getGpuFatZero(is_mono); //   30.0;
/*		
		double [] scales = (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.gpu_weight_r, // 0.25
				clt_parameters.gpu_weight_b, // 0.25
				1.0 - clt_parameters.gpu_weight_r - clt_parameters.gpu_weight_b}); // 0.5
*/
		double corr_green = 1.0/(1.0 + clt_parameters.gpu_weight_r + clt_parameters.gpu_weight_b);    // green color
		// the same, as col_weights, for separate debugging
		double [] scales= (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.gpu_weight_r *  corr_green,
				clt_parameters.gpu_weight_b * corr_green,
				corr_green});
		
		
		double cwgreen = 1.0/(1.0 + clt_parameters.corr_red + clt_parameters.corr_blue);    // green color
		double [] col_weights= (is_mono) ? (new double [] {1.0}) :(new double [] {
				clt_parameters.corr_red *  cwgreen,
				clt_parameters.corr_blue * cwgreen,
				cwgreen});

		ImageDtt image_dtt = new ImageDtt(
				quadCLT_main.getNumSensors(),
				clt_parameters.transform_size,
				clt_parameters.img_dtt,
				is_aux,
				is_mono,
				is_lwir,
				1.0);
		image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		float [][] lpf_rgb = new float[][] {
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_r),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_b),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_g),
			image_dtt.floatGetCltLpfFd(clt_parameters.gpu_sigma_m)
		};
		quadCLT_main.getGPU().setLpfRbg( // constants memory - same for all cameras
				lpf_rgb,
				debugLevel > -1); // -3

		float [] lpf_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrSigma(is_mono));

		quadCLT_main.getGPU().setLpfCorr(// constants memory - same for all cameras
				"lpf_corr", // String const_name, // "lpf_corr"
				lpf_flat,
				debugLevel > -1); // -3

		float [] lpf_rb_flat = image_dtt.floatGetCltLpfFd(clt_parameters.getGpuCorrRBSigma(is_mono));
		quadCLT_main.getGPU().setLpfCorr(// constants memory - same for all cameras
				"lpf_rb_corr", // String const_name, // "lpf_corr"
				lpf_rb_flat,
				debugLevel > -1); // -3

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
		
		
		TpTask [] tp_tasks  = quadCLT_main.getGPU().setFullFrameImages( // when disparities array is known - use different arguments of setFullFrameImages
				twoi,                                 // Rectangle                 woi,
				clt_parameters.gpu_woi_round,         // boolean                   round_woi,
	    		(float) clt_parameters.disparity,     // float                     target_disparity, // apply same disparity to all tiles
	    		(float) disparity_corr,               // float                     disparity_corr, // add to disparity (at infinity)
	    		0xf,                                  // int                       out_image, // from which tiles to generate image (currently 0/1)
	    		0x3f,                                 // int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
	    		debugLevel);                          // final int                 debugLevel) - not yet used
		if (tp_tasks.length == 0) {
			System.out.println("--- Empty tasks - nothing to do ---");
			return null;
		}

		quadCLT_main.getGPU().setTasks( // copy tp_tasks to the GPU memory
				tp_tasks, // TpTask [] tile_tasks,
				use_aux, // boolean use_aux)
				clt_parameters.img_dtt.gpu_verify); // boolean verify
				
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
			quadCLT_main.getGPU().execSetTilesOffsets(true); // prepare tiles offsets in GPU memory // calculate tile centers
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

		// FIXME: Temporary setting all pairs to correlate
		int num_pairs = quadCLT_main.getGPU().getNumPairs();
		boolean [] sel_pairs = new boolean [num_pairs];
		Arrays.fill(sel_pairs, true);
		
		long startCorr2d=System.nanoTime();   // System.nanoTime();
		
		for (int i = 0; i < NREPEAT; i++ ) {
			//Generate 2D phase correlations from the CLT representation
			quadCLT_main.getGPU().execCorr2D(
					sel_pairs,                    // boolean [] pair_select,
					scales,                       // double [] scales,
					fat_zero,                     // double fat_zero);
					clt_parameters.gpu_corr_rad); // int corr_radius
		}
		long endCorr2d = System.nanoTime();
		
		// Should be done before execCorr2D_TD as corr_indices are shared to save memory
		int [] corr_indices = quadCLT_main.getGPU().getCorrIndices();
		// the following is not yet shared
		float [][] corr2D = quadCLT_main.getGPU().getCorr2D(
				clt_parameters.gpu_corr_rad); //  int corr_rad);
		
// calculate correlations, keep TD
		int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt,quadCLT_main.getNumSensors());
		int [] i_mcorr_sel = Correlation2d.intCorrPairs(
				mcorr_sel,
				quadCLT_main.getNumSensors(),
				4); // int num_out); should be 4 int

		quadCLT_main.getGPU().execCorr2D_TD(
	    		scales,i_mcorr_sel);// double [] scales,
		
		quadCLT_main.getGPU().execCorr2D_combine( // calculate cross pairs
		        true, // boolean init_corr,    // initialize output tiles (false - add to current)
		        6,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
		        0x0f); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
		
		quadCLT_main.getGPU().execCorr2D_normalize(
        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
	    		fat_zero, // double fat_zero);
				null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
	    		clt_parameters.gpu_corr_rad); // int corr_radius
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
				calc_extra,                    // boolean   calc_extra)
				false);                        // boolean   linescan_order) // TODO: use true to avoid reordering of the low-res output 
				

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
//		float [][][] iclt_fimg = new float [quadCLT_main.getGPU().getNumCams()][][];
		float [][][] iclt_fimg = new float [quadCLT_main.getNumSensors()][][];
		
		for (int ncam = 0; ncam < iclt_fimg.length; ncam++) {
			iclt_fimg[ncam] = quadCLT_main.getGPU().getRBG(ncam);
		}

		int out_width =  quadCLT_main.getGPU().getImageWidth()  + quadCLT_main.getGPU().getDttSize();
		int out_height = quadCLT_main.getGPU().getImageHeight() + quadCLT_main.getGPU().getDttSize();
		int tilesX =     quadCLT_main.getGPU().getImageWidth()  / quadCLT_main.getGPU().getDttSize();
		int tilesY =     quadCLT_main.getGPU().getImageHeight() / quadCLT_main.getGPU().getDttSize();
		
		if (clt_parameters.gpu_show_geometry) {
			//			GPUTileProcessor.TpTask []
			tp_tasks = quadCLT_main.getGPU().getTasks (false); // boolean use_aux)
			double [][] geom_dbg = new double [ImageDtt.GEOM_TITLES_DBG.length][tilesX*tilesY];
			int num_cams = quadCLT_main.getNumSensors(); // GPUTileProcessor.numSensors; // NUM_CAMS;
			for (int nt = 0; nt < tp_tasks.length; nt++) {
				for (int i = 0; i < num_cams; i++) {
					TpTask task = tp_tasks[nt];
					int nTile = task.ty * tilesX + task.tx;
					geom_dbg[2 * i + 0][nTile] = task.xy[i][0]; // x
					geom_dbg[2 * i + 1][nTile] = task.xy[i][1]; // y
					for (int j = 0; j < 4; j++) {
						geom_dbg[2 * num_cams + 4 * i + j][nTile] = task.disp_dist[i][j];
					}
				}
			}
		    (new ShowDoubleFloatArrays()).showArrays(
		            geom_dbg,
		            tilesX,
		            tilesY,
		            true,
		            name+"-GEOM-DBG-D"+clt_parameters.disparity,
		            ImageDtt.GEOM_TITLES_DBG);
			
		}

		/*
        public TpTask [] getTasks (boolean use_aux)
        {
        	float [] ftasks = new float [TPTASK_SIZE * num_task_tiles];
        	cuMemcpyDtoH(Pointer.to(ftasks), gpu_tasks, TPTASK_SIZE * num_task_tiles * Sizeof.FLOAT);
        	TpTask [] tile_tasks = new TpTask[num_task_tiles];
        	for (int i = 0; i < num_task_tiles; i++) {
        		tile_tasks[i] = new TpTask(ftasks, i* TPTASK_SIZE, use_aux);
        	}
        	return tile_tasks;
        }

		 */
		
		
		// Read extra data for macro generation: 4 DIFFs, 4 of R,  4 of B, 4 of G
		// Available after gpu_diff_rgb_combo is generated in execTextures
		if (calc_extra) {
			String [] extra_group_titles = {"DIFF","Red","Blue","Green"};
			int numSensors = quadCLT_main.getNumSensors();
			String [] extra_titles = new String [extra_group_titles.length* numSensors];
			for (int g = 0; g < extra_group_titles.length;g++) {
				for (int ncam=0; ncam < numSensors;ncam++) {
					extra_titles[g * numSensors + ncam]= extra_group_titles[g]+"-"+ncam;
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
			int [] corr_quad_indices = quadCLT_main.getGPU().getCorrComboIndices(); // get quad
			float [][] corr2D_quad = quadCLT_main.getGPU().getCorr2DCombo(clt_parameters.gpu_corr_rad);
// calculate and get cross here!			
			quadCLT_main.getGPU().execCorr2D_combine( // calculate cross pairs
			        true, // boolean init_corr,    // initialize output tiles (false - add to current)
			        6,    // int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			        0x30); // int     pairs_mask    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			
			quadCLT_main.getGPU().execCorr2D_normalize(
	        		true, // boolean combo, // normalize combo correlations (false - per-pair ones) 
		    		fat_zero, // double fat_zero);
					null, // float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
		    		clt_parameters.gpu_corr_rad); // int corr_radius
			
			int [] corr_cross_indices = quadCLT_main.getGPU().getCorrComboIndices(); // get quad
			float [][] corr2D_cross = quadCLT_main.getGPU().getCorr2DCombo(clt_parameters.gpu_corr_rad);
			
			double [][] dbg_corr_pairs = GPUTileProcessor.getCorr2DView(
					quadCLT_main.getNumSensors(),
					tilesX,
					tilesY,
					corr_indices,
					corr2D,
					wh);
			double [][] dbg_corr_quad = GPUTileProcessor.getCorr2DView(
					quadCLT_main.getNumSensors(),
					tilesX,
					tilesY,
					corr_quad_indices,
					corr2D_quad,
					wh);

			double [][] dbg_corr_cross = GPUTileProcessor.getCorr2DView(
					quadCLT_main.getNumSensors(),
					tilesX,
					tilesY,
					corr_cross_indices,
					corr2D_cross,
					wh);
			
			double [][] dbg_corr = {
					dbg_corr_pairs[0],
					dbg_corr_pairs[1],
					dbg_corr_pairs[2],
					dbg_corr_pairs[3],
					dbg_corr_pairs[4],
					dbg_corr_pairs[5],
					dbg_corr_quad[15],
					dbg_corr_cross[48]
			};

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
			String x3d_path= quadCLT_main.getX3dDirectory();
			for (int sub_img = 0; sub_img < imps_RGB.length; sub_img++){
				EyesisCorrections.saveAndShow(
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
			double [][][][] texture_tiles = new double [tilesY][tilesX][][];
			quadCLT_main.getGPU().doubleTextures(
		    		new Rectangle(0, 0, tilesX, tilesY), // Rectangle    woi,
		    		texture_tiles,                       // double [][][][] texture_tiles, // null or [tilesY][tilesX]
		    		texture_indices,                     // int []       indices,
		    		flat_textures,                       // float [][][] ftextures,
		    		tilesX,                              // int          full_width,
		    		4,                                   // rbga only /int          num_slices
		    		num_src_slices                       // int          num_src_slices
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
			TpTask [] tp_tasks
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
	
// overwrites super method	// Restored super - it now includes GPU
	  @Deprecated
	  public CLTPass3d CLTBackgroundMeas_Old( // measure background // USED in lwir
			  CLTParameters       clt_parameters,
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return super.CLTBackgroundMeas( // measure background // USED in lwir
					  clt_parameters,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
		  }
			  
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

		  // yes, needed  (for macro)
//		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)
		  double [][] disparity_map = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)
		  double [][][][] texture_tiles =     new double [tilesY][tilesX][][]; // ["RGBA".length()][];
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  
		  int disparity_modes = 
				  ImageDtt.BITS_ALL_DISPARITIES |
				  ImageDtt.BITS_ALL_DIFFS |
				  ImageDtt.BITS_OVEREXPOSED |
				  ImageDtt.BITS_TONE_RGB;
		  double [][][][][] clt_corr_partial = null;
		  if (clt_parameters.img_dtt.gpu_mode_debug) {
			  clt_corr_combo = null;
			  clt_corr_partial = new double [tilesY][tilesX][][][];
			  for (int i = 0; i < tilesY; i++){
				  for (int j = 0; j < tilesX; j++){
					  clt_corr_partial[i][j] = null;
				  }
			  }
		  }
		  
//		  float [][] texture_img = new float [isMonochrome()?2:4][];
//		  Rectangle  texture_woi = new Rectangle();
		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]	  fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  clt_corr_combo,                // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  clt_corr_partial,              // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  disparity_map,                 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  disparity_modes,               // disparity_modes, // bit mask of disparity_map slices to calculate/return
					
				  texture_tiles,                 // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					// new parameters, will replace some other?
//				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);

		  scan_rslt.disparity =     disparity_array;
		  scan_rslt.tile_op =       tile_op;
		  scan_rslt.disparity_map = disparity_map;
		  scan_rslt.texture_tiles = null; // texture_tiles;
//		  scan_rslt.texture_img =   texture_img;
//		  scan_rslt.texture_woi =   texture_woi;
		  scan_rslt.is_measured =   true;
		  scan_rslt.is_combo =      false;
		  scan_rslt.resetProcessed();
		  if (clt_corr_partial!=null){ // only to debug matching gpu/cpu
			  if (debugLevel > -1){ // -1
				  double [][] corr_rslt_partial = ImageDtt.corr_partial_dbg(
						  clt_corr_partial,
						  2*image_dtt.transform_size - 1,	//final int corr_size,
						  4,	// final int pairs,
						  4,    // final int colors,
						  clt_parameters.corr_border_contrast,
						  threadsMax,
						  debugLevel);
				  // titles.length = 15, corr_rslt_partial.length=16!
				  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						  corr_rslt_partial,
						  tilesX*(2*image_dtt.transform_size),
						  tilesY*(2*image_dtt.transform_size),
						  true,
						  image_name+sAux()+"-PART_CORR-GPU-D"+clt_parameters.disparity);
			  }
		  }
		  
		  return scan_rslt;
	  }

	  public ImagePlus getBackgroundImage(
			  boolean []                                bgnd_tiles, 
			  CLTParameters                             clt_parameters,
			  ColorProcParameters                       colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters  rgbParameters,
			  String     name,
//			  int        disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel
			  ) {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return super.getBackgroundImage( // measure background // USED in lwir
					  bgnd_tiles, 
					  clt_parameters,
					  colorProcParameters,
					  rgbParameters,
					  name,
//					  disparity_index, // index of disparity value in disparity_map == 2 (0,2 or 4)
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
		  }
		  
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  CLTPass3d bgnd_data = tp.clt_3d_passes.get(0);
		  boolean [] bgnd_tiles_grown2 = bgnd_data.getSelected().clone(); // only these have non 0 alpha
		  tp.growTiles(
				  2,      // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
				  bgnd_tiles_grown2,
				  null); // prohibit
		  //				  bgnd_data.getTextureSelection()
		  bgnd_data.setTextureSelection(bgnd_tiles_grown2);
		  Rectangle  texture_woi_pix = new Rectangle(); // in pixels
		  float [][] texture_img = GetTextureGPU( // returns texture CUDA error at com.elphel.imagej.gpu.GpuQuad.execRBGA_noDP(GpuQuad.java:2180)
				  clt_parameters,      // CLTParameters       clt_parameters,
				  texture_woi_pix,     // Rectangle  texture_woi, // = new Rectangle();
				  bgnd_data.disparity, // double [][]         disparity_array, // [tilesY][tilesX]
				  bgnd_tiles_grown2,   // bgnd_tiles_grown2,   // boolean []          selection,
				  threadsMax,          // final int           threadsMax,  // maximal number of threads to launch
				  updateStatus,        // final boolean       updateStatus,
				  debugLevel);         // final int           debugLevel);
		  int out_width =  tilesX *  clt_parameters.transform_size;
		  int out_height = tilesY *  clt_parameters.transform_size;
		  ImagePlus imp_texture_full = null;
		  if (texture_img == null) {
			  System.out.println("No background visible");
			  texture_woi_pix = new Rectangle(0,0,1,1);
			  texture_img = new float[][] {{0.0f},{0.0f}};
		  }

		  // for now - use just RGB. Later add option for RGBA
		  float [][] texture_rgb = isMonochrome() ? (new float [][] {texture_img[0]}): (new float [][] {texture_img[0],texture_img[1],texture_img[2]});
		  float [][] texture_rgba = isMonochrome() ? (new float [][] {texture_img[0],texture_img[1]}) : (new float [][] {texture_img[0],texture_img[1],texture_img[2],texture_img[3]});

			
		  ImagePlus imp_texture_bgnd = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture-bgnd", // String name,
				  "", //String suffix, // such as disparity=...
				  true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb),
				  texture_woi_pix.width, // tilesX *  clt_parameters.transform_size,
				  texture_woi_pix.height, // tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);
//		  imp_texture_bgnd.show();
		  // resize for backdrop here!
//	public double getFOVPix(){ // get ratio of 1 pixel X/Y to Z (distance to object)
		  imp_texture_full = resizeToFull(
				  out_width, // int       out_width,
				  out_height, // int       out_height,
				  texture_woi_pix.x, // int       x0, // image offset-x 
				  texture_woi_pix.y, // int       y0, // image offset-y
				  imp_texture_bgnd, // ImagePlus imp_texture_bgnd,
				  false, // boolean   fillBlack,
				  false, // boolean noalpha, // only with fillBlack, otherwise ignored
				  debugLevel);
//		  imp_texture_full.show();
		  return imp_texture_full;

	  }
/*
	  public CLTPass3d  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters     clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel) {
		  return  CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity // not used in lwir
				  clt_parameters,
				  scanIndex,
				  save_textures,
				  0, // clust_radius,
				  threadsMax,  // maximal number of threads to launch
				  updateStatus,
				  debugLevel);
	  }
*/	  
	  
	  public CLTPass3d  CLTMeasureCorr_replaced( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters     clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         clust_radius,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  // TODO: implement clust_radius for GPU !
			  return super.CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity
					  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					  scanIndex,      // final int         scanIndex,
					  save_textures,  // final boolean     save_textures,
					  clust_radius,   // final int         clust_radius,
					  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
					  updateStatus,   // final boolean     updateStatus,
					  debugLevel);    // final int         debugLevel);
		  }
		  // GPU will not use textures. Maybe set texture_selection instead?
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;
		  boolean [] tex_sel = null;
		  if (save_textures) {
			  tex_sel = new boolean [tilesY*tilesX];
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  if (tile_op[ty][tx] != 0){
						  tex_sel[ty * tilesX + tx] = true;
					  }
				  }
			  }
		  }
		  double [][] disparity_map = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][];
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  // TODO: find where BITS_ALL_DIFFS is not needed and remove - it takes much longer than correlation in the GPU
		  int disparity_modes = 
				  ImageDtt.BITS_ALL_DISPARITIES |
				  ImageDtt.BITS_ALL_DIFFS | // needs max_diff?
				  ImageDtt.BITS_OVEREXPOSED; //  |
		  if (clust_radius > 0) { 
			  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					  tile_op,                       // per-tile operation bit codes
					  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
					  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  //// Uses quadCLT from gpuQuad			
					  // correlation results - final and partial
					  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [type][tilesY][tilesX] should be set by caller
					  // types: 0 - selected correlation (product+offset), 1 - sum
					  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,                 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,               // disparity_modes, // bit mask of disparity_map slices to calculate/return

					  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
					  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
					  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
					  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					  // new parameters, will replace some other?
					  //				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
					  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
					  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
					  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
					  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.corr_red,       // +used
					  clt_parameters.corr_blue,      // +used
					  clt_parameters.max_corr_radius,// 3.9;
					  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
					  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_parameters.clt_window,
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
			  
		  } else {
			  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					  tile_op,                       // per-tile operation bit codes
					  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
					  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  //// Uses quadCLT from gpuQuad			
					  // correlation results - final and partial
					  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [type][tilesY][tilesX] should be set by caller
					  // types: 0 - selected correlation (product+offset), 1 - sum
					  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,                 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,               // disparity_modes, // bit mask of disparity_map slices to calculate/return

					  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
					  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
					  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
					  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					  // new parameters, will replace some other?
					  //				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
					  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
					  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
					  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
					  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.corr_red,       // +used
					  clt_parameters.corr_blue,      // +used
					  clt_parameters.max_corr_radius,// 3.9;
					  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
					  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_parameters.clt_window,
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
		  }	  
 		  scan.disparity_map = disparity_map;
		  scan.texture_tiles = null;
		  scan.setTextureSelection(tex_sel);
		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  scan.resetProcessed();
		  return scan;
	  }

	  /**
	   * Get GPU-calculated tiles ceneters coordinates for each sesnor
	   * @return per-tile null 2*num_sesnors array of alternating X and Y 
	   */
	  public double [][] getPortsXY() // run after GPU 
	  {
		  boolean use_aux = false;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  double [][] ports_xy = new double [tilesY*tilesX][];
		  TpTask [] tp_tasks = gpuQuad.getTasks(use_aux); // use_aux);
		  for (TpTask tp_task:tp_tasks) {
			  float [][] txy =  tp_task.getXY(); // use_aux);
			  double [] tile_ports_xy = new double [txy.length * txy[0].length];
			  int indx = 0;
			  for (int i = 0; i < txy.length; i++) {
				  for (int j = 0; j < txy[i].length; j++) {
					  tile_ports_xy[indx++] = txy[i][j];
				  }
			  }
			  ports_xy[tp_task.getTileY()*tilesX + tp_task.getTileX()] = tile_ports_xy;
		  }
		  return ports_xy;
	  }

	  /**
	   * Get GPU-calculated tiles centerXY (tile centers for the virtual camera)
	   * @return per-tile null or pair of X,Y
	   */
	  public double [][] getCenterXY() // run after GPU 
	  {
		  boolean use_aux = false;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

		  double [][] centers_xy = new double [tilesY*tilesX][];
		  TpTask [] tp_tasks = gpuQuad.getTasks(use_aux); // use_aux);
		  for (TpTask tp_task:tp_tasks) {
			  centers_xy[tp_task.getTileY()*tilesX + tp_task.getTileX()] = tp_task.getDoubleCenterXY(); // use_aux);
		  }
		  return centers_xy;
	  }

	  public void showPortsXY(String title) {
		  int tilesX = tp.getTilesX();
		  int tilesY = tp.getTilesY();
		  int num_sensors = getNumSensors();
		  String [] titles = new String[2*num_sensors];
		  for (int n = 0; n < num_sensors; n++) {
			  titles[n]= "x"+n;
			  titles[n+num_sensors]= "y"+n;
		  }
		  double [][] ports_xy =  getPortsXY();
		  double [][] data = new double [2*num_sensors][tilesX*tilesY];
		  for (int l = 0; l < data.length; l++)  Arrays.fill(data[l], Double.NaN);
		  for (int i = 0; i < ports_xy.length; i++) if (ports_xy[i]!= null) {
			  for (int j = 0; j < data.length; j++) {
				  data[j][i] = ports_xy[i][2 * (j % num_sensors)  + (j / num_sensors)];
			  }
		  }
		  (new ShowDoubleFloatArrays()).showArrays( 
				  data,
				  tilesX,
				  tilesY,
				  true,
				  title,
				  titles);
	  }

	  public void showCentersXY(String title) {
		  int tilesX = tp.getTilesX();
		  int tilesY = tp.getTilesY();
		  String [] titles = {"centerX","centerY"};
		  double [][] centers_xy =  getCenterXY();
		  double [][] data = new double [2][tilesX*tilesY];
		  for (int l = 0; l < data.length; l++)  Arrays.fill(data[l], Double.NaN);
		  for (int i = 0; i < centers_xy.length; i++) if (centers_xy[i]!= null) {
			  data[0][i] = centers_xy[i][0];
			  data[1][i] = centers_xy[i][1];
		  }
		  (new ShowDoubleFloatArrays()).showArrays( 
				  data,
				  tilesX,
				  tilesY,
				  true,
				  title,
				  titles);
	  }
	  
	  public void showCentersPortsXY(String title) {
		  int tilesX = tp.getTilesX();
		  int tilesY = tp.getTilesY();
		  int num_sensors = getNumSensors();
		  String [] titles = new String[2*num_sensors + 2];
		  titles[0]= "centerX" ;
		  titles[num_sensors + 1]= "centerY" ;
		  for (int n = 0; n < num_sensors; n++) {
			  titles[n + 1]= "x"+n ;
			  titles[n + 1 + (num_sensors + 1)]= "y"+n;
		  }
		  double [][] ports_xy =  getPortsXY();
		  double [][] centers_xy =  getCenterXY();

		  double [][] data = new double [2*(num_sensors + 1)][tilesX*tilesY];
		  for (int l = 0; l < data.length; l++)  Arrays.fill(data[l], Double.NaN);
		  for (int i = 0; i < ports_xy.length; i++) if (ports_xy[i]!= null) {
			  for (int j = 0; j < data.length; j++) {
				  int dir = j / (num_sensors + 1); // 0 - X, 1 - Y
				  int n = j - dir * (num_sensors + 1);
				  if (n == 0) {
					  data[j][i] = centers_xy[i][dir];
				  } else {
					  n --;
					  data[j][i] = ports_xy[i][2 * n  + dir];
				  }
			  }
		  }
		  (new ShowDoubleFloatArrays()).showArrays( 
				  data,
				  tilesX,
				  tilesY,
				  true,
				  title,
				  titles);
	  }

	  
	  
	  
	  public double [][] CLTCorrDisparityMap(
			  CLTParameters     clt_parameters,
			  double []         disparity, // do not measure NaN
			  double            disparity_corr,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return null; // no non-GPU 
		  }
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  int op = ImageDtt.setImgMask(0, 0xf);
		  op =     ImageDtt.setPairMask(op,0xf);
		  op =     ImageDtt.setForcedDisparity(op,true);
		  int [][] tile_op = new int [tilesY][tilesX];
		  double [][]  disparity_array = new double[tilesY][tilesX];
		  for (int ty = 0; ty < tilesY; ty++) {
			  for (int tx = 0; tx < tilesX; tx++) {
				  int indx = ty*tilesX+tx;
				  if (!Double.isNaN(disparity[indx])) {
					  tile_op        [ty][tx] = op; // 511; // d;
					  disparity_array[ty][tx] = disparity[indx];
					  
				  }
			  }
		  }
//		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][];
		  double [][] disparity_map = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][];

		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
//		  double z_correction =  clt_parameters.z_correction;
//		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
//			  z_correction +=clt_parameters.z_corr_map.get(image_name);
//		  }
		  // final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  // TODO: find where BITS_ALL_DIFFS is not needed and remove - it takes much longer than correlation in the GPU
		  int disparity_modes = 
				  ImageDtt.BITS_ALL_DISPARITIES |
				  ImageDtt.BITS_ALL_DIFFS | // needs max_diff?
				  ImageDtt.BITS_OVEREXPOSED; //  |
//				  ImageDtt.BITS_TONE_RGB;

		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  disparity_map,                 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  disparity_modes,               // disparity_modes, // bit mask of disparity_map slices to calculate/return
				  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  return disparity_map;
	  }
	  
	  @Deprecated
	  public double [][] CLTCorrDisparityMapCPUTasks(
			  String            suffix,
			  CLTParameters     clt_parameters,
			  final double []   debug_disparity_range,
			  double []         disparity, // do not measure NaN
			  double            disparity_corr,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  final boolean show_per_scene = false; // true; // false; //true; // false;
		  final boolean hv_lma =  (debugLevel > -100) && true; // use LMA with hor/vert pairs instead of the regular hor/vert method
		  final double fat_zero_cpu = (debugLevel > -100) ? clt_parameters.getGpuFatZero(isMonochrome()) : 0.0; // > 0 - use CPU normalization, DCT2 and unfolding. ==0 - GPU
		  //			final double fat_zero_pre = (debug_level>-100)?-1.0:0.0; //100000.0; // 10000.0; // 1000.0; //
		  //		  final double output_amplitude = (debugLevel < -100) ? 30000 : 0.0; // 10000; // 1000; // 200; // 50.0;
		  final double output_amplitude = 1.0; // 10000; // 1000; // 200; // 50.0;

		  final boolean opposite_mode = true; // false - parallel shift
		  final double offset_x = 0.0; // 1.0; // 0.0; // 0.5; // 0.0; // 0.0; // 0.1; // 0.0; // 0.1; // 0.1; // 0.0; // 0.1; // 0.0; // 0.1;
		  final double offset_y = 0.0; // 1.0; // 1.0; // 0.0; // 0.5; // 0.0; // 0.1; // -0.1; // 0.1; // 0.0; // 0.1; // 0.0;
//		  final Rectangle tile_woi = new Rectangle(250,114,32,24); // for visualizations
		  final Rectangle tile_woi = new Rectangle(250,114,50,24); // for visualizations
		  final int vis_gap = 2;

		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return null; // no non-GPU 
		  }
		  if (gpuQuad.getQuadCLT() != this) {
			  gpuQuad.updateQuadCLT(this); // to re-load new set of Bayer images to the GPU
		  }
		  final int margin = 8;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final int tileSize = tp.getTileSize();
		  final int tiles = tilesY*tilesX;
		  final int transform_size = tp.getTileSize();
		  final double [][]  pXpYD=new double[tilesY*tilesX][];            // per-tile array of pX,pY,disparity triplets
		  for (int nTile = 0; nTile < tiles; nTile++) if (!Double.isNaN(disparity[nTile])) {
			  int tileX = nTile % tilesX;
			  int tileY = nTile / tilesX;
			  double	centerX = tileX * tileSize + tileSize/2; //  - shiftX;
			  double	centerY = tileY * tileSize + tileSize/2; //  - shiftY;
			  pXpYD[nTile] = new double[] {centerX, centerY, disparity[nTile]};
		  }
		  float  [][][][]     fcorr_td =       new float[tilesY][tilesX][][];
//		  float  [][][][]     fcorr_combo_td = hv_lma ? null : new float[4][tilesY][tilesX][];
		  float  [][][][]     fcorr_combo_td = new float[4][tilesY][tilesX][];
		  final double gpu_sigma_corr =     clt_parameters.getGpuCorrSigma(isMonochrome());
		  final double gpu_sigma_rb_corr =  isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr;
		  final double gpu_sigma_log_corr = clt_parameters.getGpuCorrLoGSigma(isMonochrome());
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  TpTask[] tp_tasks =  gpuQuad.setInterTasks(
				  false, // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				  pXpYD,              // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				  geometryCorrection, // final GeometryCorrection  geometryCorrection,
				  disparity_corr,     // final double              disparity_corr,
				  margin,             // final int                 margin,      // do not use tiles if their centers are closer to the edges
				  null,               // final boolean []          valid_tiles,            
				  threadsMax);        // final int                 threadsMax)  // maximal number of threads to launch
		  // TODO: test-modify coordinates here
		  if (opposite_mode) {
			  for (TpTask tp_task:tp_tasks) {
				  tp_task.xy[0][0] -= offset_x;
				  tp_task.xy[1][0] += offset_x;
				  tp_task.xy[2][0] -= offset_x;
				  tp_task.xy[3][0] += offset_x;
				  tp_task.xy[0][1] -= offset_y;
				  tp_task.xy[1][1] -= offset_y;
				  tp_task.xy[2][1] += offset_y;
				  tp_task.xy[3][1] += offset_y;
			  }
		  } else {
			  for (TpTask tp_task:tp_tasks) {
				  tp_task.xy[0][0] += offset_x;
				  tp_task.xy[1][0] += offset_x;
				  tp_task.xy[2][0] += offset_x;
				  tp_task.xy[3][0] += offset_x;
				  tp_task.xy[0][1] += offset_y;
				  tp_task.xy[1][1] += offset_y;
				  tp_task.xy[2][1] += offset_y;
				  tp_task.xy[3][1] += offset_y;
			  }
		  }
		  image_dtt.quadCorrTD(
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  tp_tasks,                      // final GPUTileProcessor.TpTask[] tp_tasks,
				  fcorr_td,                      // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  fcorr_combo_td,                // final float  [][][][]     fcorr_combo_td,  // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
				  getErsCorrection(),            // final GeometryCorrection  geometryCorrection,
				  //					disparity_corr,                // final double              disparity_corr,   // disparity offset at infinity
				  margin,                        // final int                 margin,           // do not use tiles if their centers are closer to the edges
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  gpu_sigma_rb_corr,             // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  gpu_sigma_corr,                //  =    0.9;gpu_sigma_corr_m
				  gpu_sigma_log_corr,            // final double              gpu_sigma_log_corr,   // hpf to reduce dynamic range for correlations
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  threadsMax,                    // final int             threadsMax,       // maximal number of threads to launch
				  debugLevel);                   // final int                 globalDebugLevel)

		  double [][] hv_offsets = getPortsXY();
		  if (show_per_scene) {
			  String [] dbg_titles = {"disp_fract", "combo", "combo-fract", "disparity"};
			  double [][] dbg_img = new double [dbg_titles.length][tiles];
			  for (int i = 0; i < dbg_img.length; i++) {
				  Arrays.fill(dbg_img[i], Double.NaN);
			  }
			  for (int nt = 0; nt < tiles; nt++) if (!Double.isNaN(disparity[nt]) && (hv_offsets[nt] != null)){
				  dbg_img[0][nt] = disparity[nt]-Math.round(disparity[nt]);
				  dbg_img[3][nt] = disparity[nt];
				  double [] hv = hv_offsets[nt];
				  // vertical expansion minus horizontal expansion
				  dbg_img[1][nt] = (-hv[1] - hv[3] + hv[5] + hv[7]) -( -hv[0] + hv[2] -hv[4] + hv[6]);
				  dbg_img[2][nt] = dbg_img[1][nt] - Math.round(dbg_img[1][nt]);

			  }
			  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					  dbg_img,
					  tilesX,
					  tilesY,
					  true,
					  "vert-minus-hor-expansion"+suffix+"-CPU",
					  dbg_titles);
		  }
		  
		  
		  
		  int [][] offset_bins = new int [2][tiles];
		  Arrays.fill(offset_bins[0], -1);
		  Arrays.fill(offset_bins[1], -1);
		  double debug_hist_error_step =   0.025;
		  double debug_hist_max_error2 =   0.6; //
		  int err_steps2 =  2*((int) Math.ceil(debug_hist_max_error2/debug_hist_error_step)) - 1;
		  double err_offs2 = -err_steps2 * debug_hist_error_step / 2.0;
		  for (int nt = 0; nt < tiles; nt++) {
			  double td = disparity[nt];
			  if (!Double.isNaN(td) && (hv_offsets[nt]!= null) && (td >= debug_disparity_range[0]) && (td < debug_disparity_range[1])) {
				  double [] offsets = new double[2];
				  offsets[0] = -0.25 * (hv_offsets[nt][2] + hv_offsets[nt][6] - hv_offsets[nt][0] - hv_offsets[nt][4]); // null
				  offsets[1] = -0.25 * (hv_offsets[nt][5] + hv_offsets[nt][7] - hv_offsets[nt][1] - hv_offsets[nt][3]);

				  for (int mode = 0; mode < 2; mode++) { // hor, vert only
					  double offs = offsets[mode];
					  //	                    if (mode == 0) {
					  //	                        offs += 0.5; // experimental - visible shift by 0.5
					  //	                    }
					  offs -= Math.round(offs);
					  if (!Double.isNaN(offs)) {
						  int x = (int) Math.floor((offs - err_offs2)/debug_hist_error_step);
						  if      (x < 0) x = 0;
						  else if (x >= err_steps2) x = err_steps2 - 1;
						  offset_bins[mode][nt] = x;
					  }
				  }    						
			  }
		  }
		  int corr_len = 4 * transform_size * transform_size;
		  double [][][][] selected_td_tiles = new double [2][err_steps2][10][corr_len];
		  int [][]        num_acc = new int [2][err_steps2];
		  for (int nt = 0; nt < tiles; nt++) {
			  int tileX = nt % tilesX;
			  int tileY = nt / tilesX;
			  int ih = offset_bins[0][nt];
			  int iv = offset_bins[1][nt];
			  if ((fcorr_td[tileY][tileX] != null) && (ih >= 0) && (iv >= 0)) {
				  for (int s = 0; s < 6; s++) {
					  float [] ftile = fcorr_td[tileY][tileX][s];
					  for (int i = 0; i < corr_len; i++) {
						  selected_td_tiles[0][ih][s][i] += ftile[i];
						  selected_td_tiles[1][iv][s][i] += ftile[i];
					  }
				  }
				  if (fcorr_combo_td != null) {
					  for (int s = 0; s < 4; s++) if (fcorr_combo_td[s] != null) {
						  float [] ftile = fcorr_combo_td[s][tileY][tileX];
						  int s1 = s + 6;
						  for (int i = 0; i < corr_len; i++) {
							  selected_td_tiles[0][ih][s1][i] += ftile[i];
							  selected_td_tiles[1][iv][s1][i] += ftile[i];
						  }
					  }
				  }
				  num_acc[0][ih] ++;
				  num_acc[1][iv] ++;
			  }
		  }
		  for (int mode = 0; mode < 2; mode++) {
			  for (int ihv = 0; ihv < err_steps2; ihv++) {
				  if (num_acc[mode][ihv] > 0) {
					  for (int s = 0; s < 10; s ++) {
						  for (int i = 0; i < corr_len; i++) {
							  selected_td_tiles[mode][ihv][s][i] /= num_acc[mode][ihv];
						  }
					  }
				  }
			  }
		  }
		  // now inject accumulated into the image fcorr_td, fcorr_combo_td
		  for (int mode = 0; mode < 2; mode++) {
			  int ty = tile_woi.y + mode;
			  for (int ihv = 0; ihv < err_steps2; ihv++) {
				  int tx = tile_woi.x + ihv;
				  if (num_acc[mode][ihv] > 0) {
					  for (int s = 0; s < 6; s ++) {
						  for (int i = 0; i < corr_len; i++) {
							  fcorr_td[ty][tx][s][i] = (float) selected_td_tiles[mode][ihv][s][i];
						  }
					  }
					  for (int s = 0; s < 4; s ++) {
						  int s1 = s + 6;
						  if (fcorr_combo_td!=null) {
							  for (int i = 0; i < corr_len; i++) {
								  fcorr_combo_td[s][ty][tx][i] = (float) selected_td_tiles[mode][ihv][s1][i];
							  }
						  }
					  }
				  }
			  }
		  }

//		  double[][] disparity_map =      new double [ImageDtt.DISPARITY_TITLES.length][];
//		  double[][] disparity_map_fake = new double [ImageDtt.DISPARITY_TITLES.length][];
		  double [][] disparity_map = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][];
		  double [][] disparity_map_fake = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][];
		  
		  int disparity_modes = 
				  ImageDtt.BITS_ALL_DISPARITIES |
				  ImageDtt.BITS_ALL_DIFFS | // needs max_diff?
				  ImageDtt.BITS_OVEREXPOSED; //  |

		  // here accumulate TD tiles according to offsets (like in a histogram) using double array - convert to float in the very end and
		  // put them in the top left corner, then process as usual 


		  double [][][][][] corrs_pd = null;
		  double [][][][][] corrs_pd_combo = null;
		  double [][][][]   corrs = null;
		  double [][][][]   corrs_combo = null;
		  int [] dbg_wh = new int [2];

		  double [][] dbg_img = new double [4][];
		  float [] fdbg_img = ImageDtt.corr_td_wnd(
				  fcorr_td,             // final float [][][][] fcorr_td,       // float[tilesY][tilesX][num_slices][];
				  fcorr_combo_td,       // final float [][][][] fcorr_combo_td, // float[4][tilesY][tilesX][];
				  tile_woi,                 // final Rectangle      woi,
				  vis_gap,                  // final int            gap,
				  dbg_wh,                   // final int []         wh,
				  image_dtt.transform_size, // final int            transform_size,
				  threadsMax);              // final int            threadsMax)     // maximal number of threads to launch
		  dbg_img[0] = new double [fdbg_img.length];
		  for (int i = 0; i < fdbg_img.length; i++) {
			  dbg_img[0][i] = fdbg_img[i];
		  }
		  if (fat_zero_cpu > 0.0) {
			  // normalize fcorr_td, fcorr_combo_td in-place
			  int num_pairs = Correlation2d.getNumPairs(getNumSensors());
			  if (fcorr_td != null) {
				  ImageDtt.corr_td_normalize(
						  fcorr_td, // final float [][][][] fcorr_td, // will be updated
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  num_pairs, // final int            num_slices,
						  image_dtt.transform_size,   // final int            transform_size,
						  fat_zero_cpu, // final double         fat_zero_abs,
						  output_amplitude, // final double         output_amplitude,
						  threadsMax); // final int            threadsMax);     // maximal number of threads to launch
				  ImageDtt.corr_td_lpf(
						  fcorr_td, // final float [][][][] fcorr_td, // will be updated
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  num_pairs, // final int            num_slices,
						  image_dtt.transform_size,   // final int            transform_size,
						  gpu_sigma_corr,             // final double         corr_sigma,
						  threadsMax);                // final int            threadsMax);     // maximal number of threads to launch

				  fcorr_td =  ImageDtt.corr_transpose(
						  fcorr_td, // final float [][][][] fcorr_td, // will be updated
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][]; // last index - [4 * 64]
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  num_pairs, // final int            num_slices,
						  image_dtt.transform_size,   // final int            transform_size,
						  threadsMax);                // final int            threadsMax);     // maximal number of threads to launch

				  corrs_pd = ImageDtt.corr_td_dct2(
						  fcorr_td, // final float [][][][] fcorr_td, // normalized TD representation
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][]; // last index - [4 * 64]
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  num_pairs, // final int            num_slices,
						  image_dtt.transform_size,   // final int            transform_size,
						  threadsMax); // final int            threadsMax);     // maximal number of threads to launch

				  corrs = ImageDtt.corr_pd_unfold(
						  corrs_pd,                   // final double [][][][][] pd_data, // pixel-domain correlation representation
						  // if 0 - pd_data = new double[4][tilesY][tilesX][4][64]; // last index - [4 * 64]
						  // if > 0 - pd_data =       new double[tilesY][tilesX][num_slices][4][64];
						  num_pairs, // final int            num_slices,
						  image_dtt.transform_size,   // final int            transform_size,
						  threadsMax); // final int            threadsMax);     // maximal number of threads to launch

			  }
			  if (fcorr_combo_td != null) {
				  ImageDtt.corr_td_normalize(
						  fcorr_combo_td,              // final float [][][][] fcorr_td, // will be updated
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  0, // final int            num_slices,
						  image_dtt.transform_size,    // final int            transform_size,
						  fat_zero_cpu,                // final double         fat_zero_abs,
						  output_amplitude,            //final double         output_amplitude,
						  threadsMax);                 // final int            threadsMax);     // maximal number of threads to launch

				  ImageDtt.corr_td_lpf(
						  fcorr_combo_td, // final float [][][][] fcorr_td, // will be updated
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][];
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  0,                           // final int            num_slices,
						  image_dtt.transform_size,    // final int            transform_size,
						  gpu_sigma_corr,              // final double         corr_sigma,
						  threadsMax);                 // final int            threadsMax);     // maximal number of threads to launch

				  fcorr_combo_td =  ImageDtt.corr_transpose(
						  fcorr_combo_td, // final float [][][][] fcorr_td, // will be updated
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][]; // last index - [4 * 64]
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  0,                           // final int            num_slices,
						  image_dtt.transform_size,    // final int            transform_size,
						  threadsMax);                 // final int            threadsMax);     // maximal number of threads to launch

				  corrs_pd_combo = ImageDtt.corr_td_dct2(
						  fcorr_combo_td,              // final float [][][][] fcorr_td, // normalized TD representation
						  // if 0 - fcorr_combo_td = new float[4][tilesY][tilesX][]; // last index - [4 * 64]
						  // if > 0 - fcorr_td =       new float[tilesY][tilesX][num_slices][];
						  0,                           // final int            num_slices,
						  image_dtt.transform_size,    // final int            transform_size,
						  threadsMax);                 // final int            threadsMax);     // maximal number of threads to launch
				  corrs_combo = ImageDtt.corr_pd_unfold(
						  corrs_pd_combo,              // final double [][][][][] pd_data, // pixel-domain correlation representation
						  // if 0 - pd_data = new double[4][tilesY][tilesX][4][64]; // last index - [4 * 64]
						  // if > 0 - pd_data =       new double[tilesY][tilesX][num_slices][4][64];
						  0,                           // final int            num_slices,
						  image_dtt.transform_size,    // final int            transform_size,
						  threadsMax);                 // final int            threadsMax);     // maximal number of threads to launch
			  }
			  // after normalization
			  fdbg_img = ImageDtt.corr_td_wnd(
					  fcorr_td,             // final float [][][][] fcorr_td,       // float[tilesY][tilesX][num_slices][];
					  fcorr_combo_td,       // final float [][][][] fcorr_combo_td, // float[4][tilesY][tilesX][];
					  tile_woi,                 // final Rectangle      woi,
					  vis_gap,                  // final int            gap,
					  dbg_wh,                   // final int []         wh,
					  image_dtt.transform_size, // final int            transform_size,
					  threadsMax);              // final int            threadsMax)     // maximal number of threads to launch
			  dbg_img[1] = new double [fdbg_img.length];
			  for (int i = 0; i < fdbg_img.length; i++) {
				  dbg_img[1][i] = fdbg_img[i];
			  }
			  dbg_img[2] = ImageDtt.corr_td_wnd(
					  corrs_pd,                 // final float [][][][] fcorr_td,       // float[tilesY][tilesX][num_slices][];
					  corrs_pd_combo,           // final float [][][][] fcorr_combo_td, // float[4][tilesY][tilesX][];
					  tile_woi,                 // final Rectangle      woi,
					  vis_gap,                  // final int            gap,
					  dbg_wh,                   // final int []         wh,
					  image_dtt.transform_size, // final int            transform_size,
					  threadsMax);              // final int            threadsMax)     // maximal number of threads to launch


			  final double  [][][]       dclt_corr = new double [tilesX * tilesY][][];

			  ImageDtt.clt_process_pd_correlations(    // process correlations already prepared in dcorr and/or dcorr_combo
					  this,                            // QuadCLT                   quadCLT,
					  clt_parameters.img_dtt,			 // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  // both arrays should have same non-null tiles
					  corrs,                           // final double  [][][][]    dcorr,           // [tilesY][tilesX][pair][225] pixel domain representation of 6 corr pairs SHOULD not == null
					  corrs_combo,                     // final double  [][][][]    dcorr_combo,     // [4][tilesY][tilesX][225] PD of combo corrs: qud, cross, hor,vert MAY BE null
					  // each of the top elements may be null to skip particular combo type
					  null, 							 //	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					  null, 							 // final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // 											  [tilesY][tilesX] should be set by caller
					  dclt_corr,                       // final double  [][][]       dcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null, 							 // final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,					 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,                 // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					  63,                              // final int                 corr_pairs_mask, // which pairs to use in LMA
					  clt_parameters.max_corr_radius,	 // final double              max_corr_radius, // 3.9;
					  clt_parameters.tileX,      		 // final int               debug_tileX,
					  clt_parameters.tileY,          	 // final int               debug_tileY,
					  threadsMax,
					  debugLevel -1 );
			  if (show_per_scene) {
				  dbg_img[3] = ImageDtt.corr_partial_wnd( // not used in lwir
						  dclt_corr,                      // clt_corr_partial,               // final double [][][][][] corr_data,
						  tilesX,                         // final int               tilesX,
						  2*image_dtt.transform_size - 1, // final int               corr_size,
						  tile_woi,                       // final Rectangle      woi,
						  vis_gap,                        // final int            gap,
						  dbg_wh,                         // final int []         wh,
						  threadsMax);                    // final int               threadsMax)

				  String [] dbg_titles = {"td","td_norm", "pd", "corr"}; 
				  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						  dbg_img,
						  dbg_wh[0],
						  dbg_wh[1],
						  true,
						  "td-pd-dbg"+suffix+"-CPU",
						  dbg_titles);
			  }
			  if (!hv_lma) {
				  return disparity_map;
			  }
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_CM] = disparity_map[ImageDtt.DISPARITY_INDEX_CM].clone();
			  disparity_map_fake[ImageDtt.DISPARITY_STRENGTH_INDEX] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
			  // using LMA with horizontal pairs, pretend it is horizontal correlation
			  ImageDtt.clt_process_pd_correlations(    // process correlations already prepared in dcorr and/or dcorr_combo
					  this,                            // QuadCLT                   quadCLT,
					  clt_parameters.img_dtt,			 // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  // both arrays should have same non-null tiles
					  corrs,                           // final double  [][][][]    dcorr,           // [tilesY][tilesX][pair][225] pixel domain representation of 6 corr pairs SHOULD not == null
					  corrs_combo,                     // final double  [][][][]    dcorr_combo,     // [4][tilesY][tilesX][225] PD of combo corrs: qud, cross, hor,vert MAY BE null
					  // each of the top elements may be null to skip particular combo type
					  null, 							 //	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					  null, 							 // final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // 											  [tilesY][tilesX] should be set by caller
					  null,                            // final double  [][][]       dcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null, 							 // final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,					 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,                 // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					  3,                               // final int                 corr_pairs_mask, // which pairs to use in LMA
					  clt_parameters.max_corr_radius,	 // final double              max_corr_radius, // 3.9;
					  clt_parameters.tileX,      		 // final int               debug_tileX,
					  clt_parameters.tileY,          	 // final int               debug_tileY,
					  threadsMax,
					  debugLevel -1 );
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_HOR] =          disparity_map[ImageDtt.DISPARITY_INDEX_CM].clone();
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
			  // using LMA with vertical pairs, pretend it is vertical correlation
			  ImageDtt.clt_process_pd_correlations(    // process correlations already prepared in dcorr and/or dcorr_combo
					  this,                            // QuadCLT                   quadCLT,
					  clt_parameters.img_dtt,			 // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  // both arrays should have same non-null tiles
					  corrs,                           // final double  [][][][]    dcorr,           // [tilesY][tilesX][pair][225] pixel domain representation of 6 corr pairs SHOULD not == null
					  corrs_combo,                     // final double  [][][][]    dcorr_combo,     // [4][tilesY][tilesX][225] PD of combo corrs: qud, cross, hor,vert MAY BE null
					  // each of the top elements may be null to skip particular combo type
					  null, 							 //	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					  null, 							 // final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // 											  [tilesY][tilesX] should be set by caller
					  null,                            // final double  [][][]       dcorr_tiles,     // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null, 							 // final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,					 // final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,                 // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					  12,                               // final int                 corr_pairs_mask, // which pairs to use in LMA
					  clt_parameters.max_corr_radius,	 // final double              max_corr_radius, // 3.9;
					  clt_parameters.tileX,      		 // final int               debug_tileX,
					  clt_parameters.tileY,          	 // final int               debug_tileY,
					  threadsMax,
					  debugLevel -1 );
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_VERT] =          disparity_map[ImageDtt.DISPARITY_INDEX_CM].clone();
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();

		  } else {
			  // convert correlation arrays to disparity map
			  //				final float  [][][]       fclt_corr = null; // new float [tilesX * tilesY][][]; // only used for debugging
			  final float  [][][]       fclt_corr = new float [tilesX * tilesY][][]; // only used for debugging
			  image_dtt.clt_process_tl_correlations_GPU(	// convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					  clt_parameters.img_dtt,				// final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  // both arrays should have same non-null tiles
					  fcorr_td,				 	 		// final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  fcorr_combo_td, 			    	// final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  null, 								//	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					  null, 								// final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  fclt_corr, // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null, 								// final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,						// final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,                    // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					  63, // final int                 corr_pairs_mask, // which pairs to use in LMA
					  clt_parameters.gpu_corr_scale,      //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  image_dtt.transform_size - 1,			// clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.max_corr_radius,		// final double              max_corr_radius, // 3.9;
					  clt_parameters.clt_window,     		// final int                 window_type,     // GPU: will not be used
					  clt_parameters.tileX,      		    // final int               debug_tileX,
					  clt_parameters.tileY,          		// final int               debug_tileY,
					  threadsMax,
					  debugLevel -1 );

			  fdbg_img = ImageDtt.corr_partial_wnd( // not used in lwir
					  fclt_corr,                      // clt_corr_partial,               // final double [][][][][] corr_data,
					  tilesX,                         // final int               tilesX,
					  2*image_dtt.transform_size - 1, // final int               corr_size,
					  tile_woi,                       // final Rectangle      woi,
					  vis_gap,                        // final int            gap,
					  dbg_wh,                         // final int []         wh,
					  threadsMax);                    // final int               threadsMax)

			  dbg_img[3] = new double [fdbg_img.length];
			  for (int i = 0; i < fdbg_img.length; i++) {
				  dbg_img[3][i] = fdbg_img[i];
			  }
			  String [] dbg_titles = {"td","td_norm", "pd", "corr"}; 
			  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					  dbg_img,
					  dbg_wh[0],
					  dbg_wh[1],
					  true,
					  "td-pd-dbg"+suffix+"-GPU",
					  dbg_titles);

			  if (!hv_lma) {
				  return disparity_map;
			  }
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_CM] = disparity_map[ImageDtt.DISPARITY_INDEX_CM].clone();
			  disparity_map_fake[ImageDtt.DISPARITY_STRENGTH_INDEX] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
			  // using LMA with horizontal pairs, pretend it is horizontal correlation
			  image_dtt.clt_process_tl_correlations_GPU(	// convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					  clt_parameters.img_dtt,				// final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  // both arrays should have same non-null tiles
					  fcorr_td,				 	 		// final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  fcorr_combo_td, 			    	// final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  null, 								//	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					  null, 								// final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  fclt_corr, // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null, 								// final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,						// final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,                    // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					  3, // final int                 corr_pairs_mask, // which pairs to use in LMA
					  clt_parameters.gpu_corr_scale,      //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  image_dtt.transform_size - 1,			// clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.max_corr_radius,		// final double              max_corr_radius, // 3.9;
					  clt_parameters.clt_window,     		// final int                 window_type,     // GPU: will not be used
					  clt_parameters.tileX,      		    // final int               debug_tileX,
					  clt_parameters.tileY,          		// final int               debug_tileY,
					  threadsMax,
					  debugLevel -1 );
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_HOR] =          disparity_map[ImageDtt.DISPARITY_INDEX_CM].clone();
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_HOR_STRENGTH] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();
			  // using LMA with vertical pairs, pretend it is vertical correlation
			  image_dtt.clt_process_tl_correlations_GPU(	// convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					  clt_parameters.img_dtt,				// final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  // both arrays should have same non-null tiles
					  fcorr_td,				 	 		// final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  fcorr_combo_td, 			    	// final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  null, 								//	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					  null, 								// final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  fclt_corr, // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null, 								// final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  disparity_map,						// final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  disparity_modes,                    // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					  12, // final int                 corr_pairs_mask, // which pairs to use in LMA
					  clt_parameters.gpu_corr_scale,      //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  image_dtt.transform_size - 1,			// clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.max_corr_radius,		// final double              max_corr_radius, // 3.9;
					  clt_parameters.clt_window,     		// final int                 window_type,     // GPU: will not be used
					  clt_parameters.tileX,      		    // final int               debug_tileX,
					  clt_parameters.tileY,          		// final int               debug_tileY,
					  threadsMax,
					  debugLevel -1 );
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_VERT] =          disparity_map[ImageDtt.DISPARITY_INDEX_CM].clone();
			  disparity_map_fake[ImageDtt.DISPARITY_INDEX_VERT_STRENGTH] = disparity_map[ImageDtt.DISPARITY_STRENGTH_INDEX].clone();



		  }
		  return disparity_map_fake;
	  }	  
	  
	  
	  
	  
	  public CLTPass3d  CLTMeasureCorrTesting( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  CLTParameters     clt_parameters,
			  final int         scanIndex,
			  final boolean     save_textures,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return super.CLTMeasureCorr( // perform single pass according to prepared tiles operations and disparity
					  clt_parameters, // EyesisCorrectionParameters.CLTParameters           clt_parameters,
					  scanIndex,      // final int         scanIndex,
					  save_textures,  // final boolean     save_textures,
					  0,                 // final int         clust_radius,
					  threadsMax,     // final int         threadsMax,  // maximal number of threads to launch
					  updateStatus,   // final boolean     updateStatus,
					  debugLevel);    // final int         debugLevel);
		  }
		  // GPU will not use textures. Maybe set texture_selection instead?
//		  final int transform_size =clt_parameters.transform_size;
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;
		  boolean [] tex_sel = null;
		  if (save_textures) {
			  tex_sel = new boolean [tilesY*tilesX];
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  if (tile_op[ty][tx] != 0){
						  tex_sel[ty * tilesX + tx] = true;
					  }
				  }
			  }
		  }
//		  double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][];
		  double [][] disparity_map = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][];		  
		  /*
		  double [][] shiftXY = new double [4][2];
		  // not used
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }
          */
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  // TODO: find where BITS_ALL_DIFFS is not needed and remove - it takes much longer than correlation in the GPU
		  int disparity_modes = 
				  ImageDtt.BITS_ALL_DISPARITIES |
				  ImageDtt.BITS_ALL_DIFFS | // needs max_diff?
				  ImageDtt.BITS_OVEREXPOSED; //  |
//				  ImageDtt.BITS_TONE_RGB;
// testing - remove later
		  float  [][][][]     corr_td = new float[tilesY][tilesX][][];
		  float  [][][][]     corr_combo_td = new float [4][tilesY][tilesX][];
		  double [][][][][] clt_corr_partial0 = null;
		  clt_corr_partial0 = new double [tilesY][tilesX][][][];
		  for (int i = 0; i < tilesY; i++){
			  for (int j = 0; j < tilesX; j++){
				  clt_corr_partial0[i][j] = null;
			  }
		  }
		  
		  
		  
		  image_dtt.clt_aberrations_quad_corr_GPU(
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  corr_td,				 		 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  corr_combo_td, 				 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  clt_corr_partial0,             // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  disparity_map,                 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  disparity_modes,               // disparity_modes, // bit mask of disparity_map slices to calculate/return
					
				  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					// new parameters, will replace some other?
//				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  
 		  scan.disparity_map = disparity_map;
		  scan.texture_tiles = null;
		  scan.setTextureSelection(tex_sel);
		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  scan.resetProcessed();
// testing
//		  float  [][][][]     corr_td = new float[tilesY][tilesX][][];
//		  float  [][][][]     corr_combo_td = new float [4][tilesY][tilesX][];
		 (new ShowDoubleFloatArrays()).showArrays(
				 disparity_map,
				 tilesX,
				 tilesY,
				 true,
				 "original_disparity_map"
//				 ,dbg_titles
				 );
		  if (debugLevel > -10){ // -1
			  double [][] corr_rslt_partial = ImageDtt.corr_partial_dbg(
					  clt_corr_partial0,
					  2*image_dtt.transform_size - 1,	//final int corr_size,
					  4,	// final int pairs,
					  4,    // final int colors,
					  clt_parameters.corr_border_contrast,
					  threadsMax,
					  debugLevel);
			  // titles.length = 15, corr_rslt_partial.length=16!
			  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					  corr_rslt_partial,
					  tilesX*(2*image_dtt.transform_size),
					  tilesY*(2*image_dtt.transform_size),
					  true,
					  image_name+sAux()+"-ORIG-PART_CORR-D"+clt_parameters.disparity);
//					  titles);
		  }

		  double [][][][][] clt_corr_partial1 = null;
		  clt_corr_partial1 = new double [tilesY][tilesX][][][];
		  for (int i = 0; i < tilesY; i++){
			  for (int j = 0; j < tilesX; j++){
				  clt_corr_partial1[i][j] = null;
			  }
		  }
//		  double [][] disparity_map1 = new double [ImageDtt.DISPARITY_TITLES.length][];
		  double [][] disparity_map1 = new double [ImageDtt.getDisparityTitles(getNumSensors(), isMonochrome()).length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)
		  
		  
		  float [][][][] corr_td_blur = image_dtt.blur_corr_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
//					final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  corr_td, // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
//					final int                 debug_tileX,
//					final int                 debug_tileY,
				  threadsMax, // final int                 threadsMax,      // maximal number of threads to launch
				  debugLevel); // final int                 globalDebugLevel)
		  
			float [][][][]  corr_combo_blur = image_dtt.blur_corr_combo_GPU( // convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
//					final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					corr_combo_td, // final float  [][][][]     fcorr_combo_td,  // [n][tilesY][tilesX][4*64] transform domain representation of combined corr pairs
//					final int                 debug_tileX,
//					final int                 debug_tileY,
					  threadsMax, // final int                 threadsMax,      // maximal number of threads to launch
					  debugLevel); // final int                 globalDebugLevel)
		  
		  

			image_dtt.clt_process_tl_correlations_GPU(	// convert to pixel domain and process correlations already prepared in fcorr_td and/or fcorr_combo_td
					clt_parameters.img_dtt,				// final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					// both arrays should have same non-null tiles
					corr_td_blur,				 	 		// final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					corr_combo_blur, 						// final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
															// each of the top elements may be null to skip particular combo type
					null, // 	final double [][][][]     corr_tiles,      // [tilesY][tilesX][pair][] ([(2*gpu_corr_rad+1)*(2*gpu_corr_rad+1)]) or null
					clt_corr_partial1,           			// final double [][][][][]   clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// [tilesY][tilesX] should be set by caller
					null, // [tile][index][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					// When clt_mismatch is non-zero, no far objects extraction will be attempted
					null, 								// final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					// values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					disparity_map1,						// final double [][]         disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					// last 2 - contrast, avg/ "geometric average)
					disparity_modes,                      // final int                 disparity_modes, // bit mask of disparity_map slices to calculate/return
					clt_parameters.gpu_corr_scale,        //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					image_dtt.transform_size - 1,			// clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					clt_parameters.max_corr_radius,		// final double              max_corr_radius, // 3.9;
					clt_parameters.clt_window,     		// final int                 window_type,     // GPU: will not be used
					clt_parameters.tileX,      		    // final int               debug_tileX,
					clt_parameters.tileY,          		// final int               debug_tileY,
					threadsMax,
					debugLevel);
		  
			 (new ShowDoubleFloatArrays()).showArrays(
					 disparity_map1,
					 tilesX,
					 tilesY,
					 true,
					 "rebuilt_disparity_map"
//					 ,dbg_titles
					 );
			 
			 
			  if (debugLevel > -10){ // -1
				  double [][] corr_rslt_partial = ImageDtt.corr_partial_dbg(
						  clt_corr_partial1,
						  2*image_dtt.transform_size - 1,	//final int corr_size,
						  4,	// final int pairs,
						  4,    // final int colors,
						  clt_parameters.corr_border_contrast,
						  threadsMax,
						  debugLevel);
				  // titles.length = 15, corr_rslt_partial.length=16!
				  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
						  corr_rslt_partial,
						  tilesX*(2*image_dtt.transform_size),
						  tilesY*(2*image_dtt.transform_size),
						  true,
						  image_name+sAux()+"-TD-PART_CORR-D"+clt_parameters.disparity);
//						  titles);
			  }
			 
		  return scan;
	  }
	  
	  
	  
	  public String getPassImage( // get image from a single pass, return relative path for x3d // USED in lwir
			  CLTParameters           clt_parameters,
			  ColorProcParameters colorProcParameters,
			  EyesisCorrectionParameters.RGBParameters             rgbParameters,
			  String     name,
			  int        scanIndex,
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  return super.getPassImage( // get image from a single pass, return relative path for x3d // USED in lwir
					  clt_parameters,
					  colorProcParameters,
					  rgbParameters,
					  name,
					  scanIndex,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
		  }

		  //		  final int tilesX = tp.getTilesX();
		  //		  final int tilesY = tp.getTilesY();
		  final int transform_size =clt_parameters.transform_size;
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  Rectangle  texture_woi_pix = new Rectangle(); // in pixels
		  float [][] texture_img = GetTextureGPU( // returns texture
				  clt_parameters,             // CLTParameters       clt_parameters,
				  texture_woi_pix,                // Rectangle  texture_woi, // = new Rectangle();
				  scan.disparity,             // double [][]         disparity_array, // [tilesY][tilesX]
				  scan.getTextureSelection(), // bgnd_tiles_grown2,   // boolean []          selection,
				  threadsMax,                 // final int           threadsMax,  // maximal number of threads to launch
				  updateStatus,               // final boolean       updateStatus,
				  debugLevel);                // final int           debugLevel);

		  // for now - use just RGB. Later add option for RGBA
		  float [][] texture_rgb = {texture_img[0],texture_img[1],texture_img[2]};
		  float [][] texture_rgba = {texture_img[0],texture_img[1],texture_img[2],texture_img[3]};
		  float [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);
		  scan.texture_bounds = new Rectangle(
				  texture_woi_pix.x/transform_size,
				  texture_woi_pix.y/transform_size,
				  texture_woi_pix.width/transform_size,
				  texture_woi_pix.height/transform_size);
		  scan.setSelected(scan.getTextureSelection().clone()); // null
		  ImagePlus imp_texture_cluster = linearStackToColor(
				  clt_parameters,
				  colorProcParameters,
				  rgbParameters,
				  name+"-texture", // String name,
				  "", //String suffix, // such as disparity=...
				  true, // toRGB,
				  !this.correctionsParameters.jpeg, // boolean bpp16, // 16-bit per channel color mode for result
				  true, // boolean saveShowIntermediate, // save/show if set globally
				  false, //true, // boolean saveShowFinal,        // save/show result (color image?)
				  texture_rgbx,
				  scan.texture_bounds.width * transform_size, // tilesX *  clt_parameters.transform_size,
				  scan.texture_bounds.height * transform_size, // tilesY *  clt_parameters.transform_size,
				  1.0,         // double scaleExposure, // is it needed?
				  debugLevel);

		  String path= correctionsParameters.selectX3dDirectory(
				  //TODO: Which one to use - name or this.image_name ?
				  correctionsParameters.getModelName(this.image_name), // quad timestamp. Will be ignored if correctionsParameters.use_x3d_subdirs is false
				  correctionsParameters.x3dModelVersion,
				  true,  // smart,
				  true);  //newAllowed, // save
		  // only show/save original size if debug or debug_filters)
		  eyesisCorrections.saveAndShow(
				  imp_texture_cluster,
				  path,
				  correctionsParameters.png,
				  clt_parameters.show_textures,
				  -1); // jpegQuality){//  <0 - keep current, 0 - force Tiff, >0 use for JPEG
		  
		  scan.is_measured =   true; // should it be here?
		  scan.is_combo =      false;
		  return imp_texture_cluster.getTitle()+".png"; // imp_texture_cluster;
	  }	

	  public void setPassAvgRBGA( // get image from a single pass, return relative path for x3d // USED in lwir
			  CLTParameters           clt_parameters,
			  CLTPass3d  scan,			  
			  int        threadsMax,  // maximal number of threads to launch
			  boolean    updateStatus,
			  int        debugLevel)
	  {
		  if ((gpuQuad != null) && (isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  final int tilesX = scan.getTileProcessor().getTilesX();
			  final int tilesY = scan.getTileProcessor().getTilesY();
			  double [] disparity = scan.getDisparity();
			  double [] strength =  scan.getStrength();
			  boolean [] selection = null; // scan.getSelected();
			  if (selection == null) {
				  selection = new boolean[tilesX*tilesY];
				  for (int nTile = 0; nTile < selection.length; nTile++) {
					  selection[nTile] = !Double.isNaN(disparity[nTile]) && (strength[nTile] > 0.0);
				  }
				  scan.setSelected(selection);
			  }
			  if ((scan.disparity == null) || (scan.tile_op == null)) {
				  int d = ImageDtt.setImgMask(0, 0xf); // no correlations
				  d =     ImageDtt.setForcedDisparity(d,true);
				  scan.setTileOpDisparity(
						  d,
						  scan.getSelected(), // boolean [] selection,
						  scan.getDisparity()); // double []  disparity)
			  }
			  
			  scan.texture_tiles = getTextureTilesGPU( // returns texture tiles, as with CPU
					  clt_parameters,     // CLTParameters       clt_parameters,
					  scan.tile_op,       // int [][]            tile_op,
					  scan.disparity,             // double [][]         disparity_array, // [tilesY][tilesX]
					  threadsMax,                 // final int           threadsMax,  // maximal number of threads to launch
					  updateStatus,               // final boolean       updateStatus,
					  debugLevel);                // final int           debugLevel);
		  } else {
			  CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
					  clt_parameters,   // final CLTParameters           clt_parameters,
					  scan,             // final CLTPass3d   scan,
					  threadsMax,
					  updateStatus,
					  debugLevel);
			  
		  }
		  super.setPassAvgRBGA( // get image from a single pass, return relative path for x3d // USED in lwir
				  clt_parameters,
				  scan,         // scanIndex,
				  threadsMax,   // maximal number of threads to launch
				  updateStatus,
				  debugLevel);
	  }
	  
			  /*
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final int transform_size =clt_parameters.transform_size;
		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  Rectangle  texture_woi = new Rectangle();
// need	to set tasks first!	  
		  float [][] texture_img = GetTextureGPU( // returns texture 
				  clt_parameters,             // CLTParameters       clt_parameters,
				  texture_woi,                // Rectangle  texture_woi, // = new Rectangle();
				  scan.disparity,             // double [][]         disparity_array, // [tilesY][tilesX]
				  scan.getTextureSelection(), // bgnd_tiles_grown2,   // boolean []          selection,
				  threadsMax,                 // final int           threadsMax,  // maximal number of threads to launch
				  updateStatus,               // final boolean       updateStatus,
				  debugLevel);                // final int           debugLevel);

		  // for now - use just RGB. Later add option for RGBA
		  float [][] texture_rgb = {texture_img[0],texture_img[1],texture_img[2]};
		  float [][] texture_rgba = {texture_img[0],texture_img[1],texture_img[2],texture_img[3]};
		  float [][] texture_rgbx = ((clt_parameters.alpha1 > 0)? texture_rgba: texture_rgb);
		  scan.texture_bounds = new Rectangle(
				  texture_woi.x * transform_size,
				  texture_woi.y * transform_size,
				  texture_woi.width * transform_size,
				  texture_woi.height * transform_size);
		  scan.selected = scan.getTextureSelection().clone();
		  int numTiles = tilesX * tilesY; 
		  int num_layers = texture_rgbx.length; // here only 3/4
		  double [][] tileTones = new double [num_layers][numTiles];
		  double [] scales = new double [num_layers];
		  for (int n = 0; n < num_layers; n++){
			  if       (n < 3)  scales[n] = 1.0/255.0; // R,B,G
			  else if  (n == 3) scales[n] = 1.0; //alpha
			  else if  (n < 8)  scales[n] = 1.0; // ports 0..3
			  else              scales[n] = 1.0/255.0; // RBG rms, in 1/255 units, but small
		  }
		  int tx0 = scan.texture_bounds.x;
		  int ty0 = scan.texture_bounds.y;
		  int tw =  scan.texture_bounds.width;
		  int th =  scan.texture_bounds.height;
		  int text_step = transform_size * tw;
		  for (int ty = ty0; ty < (ty0 + th); ty++) {
			  for (int tx = tx0; tx < (tx0 + tw); tx++) {
				  int text_start = ((ty - ty0) * transform_size * tw  + (tx-tx0)) * transform_size;
				  double s = 0.0;
				  for (int n = 0; n < num_layers; n++) {
					  for (int i = 0; i < transform_size; i++) {
						  for (int j = 0; j < transform_size; j++) {
							 s += texture_rgbx[n][text_start + i * text_step + j]; 
						  }
					  }
					  tileTones[n][ty * tilesX + tx] = s * scales[n] /(transform_size * transform_size);
				  }
			  }
		  }
		  scan.setTilesRBGA(tileTones); 
		  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
				  tileTones,
				  tilesX,
				  tilesY,
				  true,
				  image_name+sAux()+"-TileTones");
	  }			  
*/


	  
	  public float [][] GetTextureGPU( // returns texture
			  CLTParameters       clt_parameters,
			  Rectangle           texture_woi_pix, // = new Rectangle(); in pixels!
			  double [][]         disparity_array, // [tilesY][tilesX]
			  boolean []          selection,
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
			  
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][] tile_op = new int [tilesY][tilesX];
		  int indx = 0;
		  if (selection != null) {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  if (selection[indx++]) {
						  tile_op[ty][tx] = d;
					  }
				  }
			  }
		  } else {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  tile_op[ty][tx] = d;
				  }
			  }
		  }

		  /*
		  double [][] shiftXY = new double [4][2];
		  // not used
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }
          */
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  float [][] texture_img = new float [isMonochrome()?2:4][];
		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]    fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  null,                          // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  0,                             // disparity_modes, // bit mask of disparity_map slices to calculate/return
				  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  texture_img,                   // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  texture_woi_pix,               // texture_woi_pix,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					// new parameters, will replace some other?
//				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  if ((texture_img[0] == null) || (texture_img[1]== null)) {
			  return null;
		  }
		  return texture_img;
	  }
	  
	  public double [][][][] getTextureTilesGPU( // returns texture tiles, as with CPU
			  CLTParameters       clt_parameters,
			  int [][]            tile_op,
			  double [][]         disparity_array, // [tilesY][tilesX]
//			  boolean []          selection,       // null - all
			  final int           threadsMax,      // maximal number of threads to launch
			  final boolean       updateStatus,
			  final int           debugLevel)
	  {
			  
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  /*
		  int d = ImageDtt.setImgMask(0, 0xf);
		  d =     ImageDtt.setForcedDisparity(d,true);
		  int [][] tile_op = new int [tilesY][tilesX];
		  int indx = 0;
		  if (selection != null) {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  if (selection[indx++]) {
						  tile_op[ty][tx] = d;
					  }
				  }
			  }
		  } else {
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  tile_op[ty][tx] = d;
				  }
			  }
		  }
         */
		  /*
		  double [][] shiftXY = new double [4][2];
		  // not used
		  if (!clt_parameters.fine_corr_ignore) {
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  shiftXY = shiftXY0;
		  }
          */
		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);
		  double [][][][] texture_tiles = new double [tilesY][tilesX][][];
		  image_dtt.clt_aberrations_quad_corr_GPU( // USED in LWIR
				  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
				  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				  tile_op,                       // per-tile operation bit codes
				  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
				  null,							 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
				  null, 						 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					                             // each of the top elements may be null to skip particular combo type
				  null,                          // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
				  null,                          //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
				  //// Uses quadCLT from gpuQuad			
				                                 // correlation results - final and partial
				  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [type][tilesY][tilesX] should be set by caller
				                                 // types: 0 - selected correlation (product+offset), 1 - sum
				  null,                          // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				                                 // [tilesY][tilesX] should be set by caller
				                                 // When clt_mismatch is non-zero, no far objects extraction will be attempted
				  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
				                                 // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
				  null,                          // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
				                                 // last 2 - contrast, avg/ "geometric average)
				  0,                             // disparity_modes, // bit mask of disparity_map slices to calculate/return
				  texture_tiles,                 // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
				  null,                          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
				  null,                          // texture_woi,     // null or generated texture location/size
				  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					// new parameters, will replace some other?
//				  clt_parameters.getFatZero(isMonochrome()),      // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
				  clt_parameters.getGpuFatZero(isMonochrome()), // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
				  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
				  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
				  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
				  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
				  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
				  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
				  clt_parameters.corr_red,       // +used
				  clt_parameters.corr_blue,      // +used
				  clt_parameters.max_corr_radius,// 3.9;
				  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
				  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
				  clt_parameters.clt_window,
				  disparity_corr,                // final double              disparity_corr, // disparity at infinity
				  clt_parameters.tileX,          // final int               debug_tileX,
				  clt_parameters.tileY,          // final int               debug_tileY,
				  threadsMax,
				  debugLevel);
		  
		  return texture_tiles;
	  }

	  
//	  public CLTPass3d  CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
	  public void CLTMeasureTextures( // perform single pass according to prepared tiles operations and disparity // not used in lwir
			  final CLTParameters           clt_parameters,
			  final CLTPass3d          scan,
			  final int         threadsMax,  // maximal number of threads to launch
			  final boolean     updateStatus,
			  final int         debugLevel)
	  {
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux : clt_parameters.gpu_use_main)) {
			  super.CLTMeasureTextures( // measure background // USED in lwir
					  clt_parameters,
					  scan,
					  threadsMax,  // maximal number of threads to launch
					  updateStatus,
					  debugLevel);
			  return;
		  }
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();

//		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  // maybe just rely on tile_op ?
		  boolean [] text_sel = new boolean [tilesX * tilesY];
		  for (int ty = 0; ty < tilesY; ty++) {
			  for (int tx = 0; tx < tilesX; tx++) {
				  text_sel[ty*tilesX + tx] =  scan.tile_op[ty][tx] != 0;
			  }
		  }
//		  scan.setTextureSelection(scan.selected); // .clone(); GPU will measure when the texture will be requested
		  scan.setTextureSelection(text_sel); // .clone(); GPU will measure when the texture will be requested
//		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  if (debugLevel > 10) {
			  tp.showScan(
					  scan,   // CLTPass3d   scan,
					  "CLTMeasureTextures->");
		  }
//		  return scan; // what about processed and other attributes?
		  
	  }
	
	  public void CLTMeasureLY( // perform single pass according to prepared tiles operations and disparity // USED in lwir
			  final CLTParameters clt_parameters,
			  final CLTPass3d     scan,
			  final int           bgIndex, // combine, if >=0
			  final int           threadsMax,  // maximal number of threads to launch
			  final boolean       updateStatus,
			  int           debugLevel)
	  {
		  final int dbg_x = -295-debugLevel;
		  final int dbg_y = -160-debugLevel;
		  final int tilesX = tp.getTilesX();
		  final int tilesY = tp.getTilesY();
		  final int cluster_size =clt_parameters.tileStep;
		  final int clustersX= (tilesX + cluster_size - 1) / cluster_size;
		  final int clustersY= (tilesY + cluster_size - 1) / cluster_size;

//		  CLTPass3d scan = tp.clt_3d_passes.get(scanIndex);
		  scan.setLazyEyeClusterSize(cluster_size);
		  boolean [] force_disparity= new boolean[clustersX * clustersY];
		  //		  scan.setLazyEyeForceDisparity(force_disparity);
		  if (bgIndex >= 0) {
			  CLTPass3d    bg_scan = tp.clt_3d_passes.get(bgIndex);
			  // if at least one tile in a cluster is BG, use BG for the whole cluster and set lazy_eye_force_disparity
			  for (int cY = 0; cY < clustersY; cY ++) {
				  for (int cX = 0; cX < clustersX; cX ++) {
					  boolean has_bg = false;
					  for (int cty = 0; (cty < cluster_size) && !has_bg; cty++) {
						  int ty = cY * cluster_size + cty;
						  if (ty < tilesY) for (int ctx = 0; ctx < cluster_size; ctx++) {
							  int tx = cX * cluster_size + ctx;
							  if ((tx < tilesX ) && (bg_scan.tile_op[ty][tx] > 0)) {
								  has_bg = true;
								  break;
							  }
						  }
					  }
					  if (has_bg) {
						  for (int cty = 0; cty < cluster_size; cty++) {
							  int ty = cY * cluster_size + cty;
							  if (ty < tilesY) for (int ctx = 0; ctx < cluster_size; ctx++) {
								  int tx = cX * cluster_size + ctx;
								  if (tx < tilesX ) {
									  scan.tile_op[ty][tx] =   bg_scan.tile_op[ty][tx];
									  scan.disparity[ty][tx] = bg_scan.disparity[ty][tx];
								  }
							  }
						  }
						  force_disparity[cY * clustersX + cX] = true;
					  }
				  }
			  }
			  scan.setLazyEyeForceDisparity(force_disparity);
		  }
		  int [][]     tile_op =         scan.tile_op;
		  double [][]  disparity_array = scan.disparity;

		  // Should not happen !
		  if (scan.disparity == null) { // not used in lwir
			  System.out.println ("** BUG: should not happen - scan.disparity == null ! **");
			  System.out.println ("Trying to recover");
			  double [] backup_disparity = scan.getDisparity(0);
			  if (backup_disparity == null) {
				  System.out.println ("** BUG: no disparity at all !");
				  backup_disparity = new double[tilesX*tilesY];
			  }
			  scan.disparity = new double[tilesY][tilesX];
			  for (int ty = 0; ty < tilesY; ty++) {
				  for (int tx = 0; tx < tilesX; tx++) {
					  scan.disparity[ty][tx] = backup_disparity[ty*tilesX + tx];
					  if (Double.isNaN(scan.disparity[ty][tx])) {
						  scan.disparity[ty][tx] = 0;
						  tile_op[ty][tx] = 0;
					  }
				  }
			  }
			  disparity_array = scan.disparity;
		  }

		  if (debugLevel > -1){
			  int numTiles = 0;
			  for (int ty = 0; ty < tile_op.length; ty ++) for (int tx = 0; tx < tile_op[ty].length; tx ++){
				  if (tile_op[ty][tx] != 0) numTiles ++;
			  }
//			  System.out.println("CLTMeasure("+scanIndex+"): numTiles = "+numTiles);
			  System.out.println("CLTMeasure(): numTiles = "+numTiles);
			  if ((dbg_y >= 0) && (dbg_x >= 0) && (tile_op[dbg_y][dbg_x] != 0)){
//				  System.out.println("CLTMeasure("+scanIndex+"): tile_op["+dbg_y+"]["+dbg_x+"] = "+tile_op[dbg_y][dbg_x]);
				  System.out.println("CLTMeasure(): tile_op["+dbg_y+"]["+dbg_x+"] = "+tile_op[dbg_y][dbg_x]);
			  }
		  }
		  double min_corr_selected = clt_parameters.min_corr;
		  double [][] shiftXY = new double [getNumSensors()][2];
		  if (!clt_parameters.fine_corr_ignore) {// invalid for AUX!
			  double [][] shiftXY0 = {
					  {clt_parameters.fine_corr_x_0,clt_parameters.fine_corr_y_0},
					  {clt_parameters.fine_corr_x_1,clt_parameters.fine_corr_y_1},
					  {clt_parameters.fine_corr_x_2,clt_parameters.fine_corr_y_2},
					  {clt_parameters.fine_corr_x_3,clt_parameters.fine_corr_y_3}};
			  // FIXME - only first 4 sensors have correction. And is it the same for aux and main? 
			  for (int i = 0; i < shiftXY0.length;i++) {
				  shiftXY[i] = shiftXY0[i];
			  }
		  }

		  ImageDtt image_dtt = new ImageDtt(
				  getNumSensors(),
				  clt_parameters.transform_size,
				  clt_parameters.img_dtt,
				  isAux(),
				  isMonochrome(),
				  isLwir(),
				  clt_parameters.getScaleStrength(isAux()),
				  gpuQuad);
		  image_dtt.getCorrelation2d(); // initiate image_dtt.correlation2d, needed if disparity_map != null  
		  
		  double z_correction =  clt_parameters.z_correction;
		  if (clt_parameters.z_corr_map.containsKey(image_name)){ // not used in lwir
			  z_correction +=clt_parameters.z_corr_map.get(image_name);
		  }
		  final double disparity_corr = (z_correction == 0) ? 0.0 : geometryCorrection.getDisparityFromZ(1.0/z_correction);

		  if (debugLevel > -3){ // now -3
			  tp.showScan(
					  scan,   // CLTPass3d   scan,
					  "LY-combo_scan-"+scan+"_post"); //String title)
		  }
		  // use new, LMA-based mismatch calculation
		  double [][] lazy_eye_data = null;
		  if ((gpuQuad == null) || !(isAux()?clt_parameters.gpu_use_aux_adjust : clt_parameters.gpu_use_main_adjust)) { // CPU
			  /*
			  lazy_eye_data = image_dtt.cltMeasureLazyEye ( // returns d,s lazy eye parameters
					  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  tile_op,                      // per-tile operation bit codes
					  disparity_array,              // clt_parameters.disparity,     // final double            disparity,
					  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
					  saturation_imp,               // boolean [][] saturation_imp, // (near) saturated pixels or null
					  null, // final double [][]         clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  null, // disparity_map,    // [12][tp.tilesY * tp.tilesX]
					  tilesX * image_dtt.transform_size, // imp_quad[0].getWidth(),       // final int width,
					  clt_parameters.getFatZero(isMonochrome()),      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					  clt_parameters.corr_red,
					  clt_parameters.corr_blue,
					  clt_parameters.getCorrSigma(image_dtt.isMonochrome()),
					  min_corr_selected, // 0.0001; // minimal correlation value to consider valid
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,
					  //					  image_dtt.transform_size,
					  clt_parameters.clt_window,
					  shiftXY, //
					  disparity_corr, // final double              disparity_corr, // disparity at infinity
					  clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileStep,         // 	final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  clt_parameters.tileX,         // final int               debug_tileX,
					  clt_parameters.tileY,         // final int               debug_tileY,
					  threadsMax,
					  debugLevel - 2);
			*/
			  lazy_eye_data = image_dtt.cltMeasureLazyEye ( // returns d,s lazy eye parameters 
					  clt_parameters.img_dtt,       // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  tile_op,                      // final int [][]            tile_op,         // [tilesY][tilesX] - what to do - 0 - nothing for this tile
					  disparity_array,              // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					  image_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
					  saturation_imp,               // final boolean [][]        saturation_imp, // (near) saturated pixels or null
					  tilesX * image_dtt.transform_size, // 	final int                 width,
					  clt_parameters.getFatZero(isMonochrome()),      // final double              corr_fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
					  clt_parameters.corr_red,      // final double              corr_red,
					  clt_parameters.corr_blue,     // final double              corr_blue,
					  clt_parameters.getCorrSigma(image_dtt.isMonochrome()), // final double              corr_sigma,
					  min_corr_selected, // 0.0001; //final double              min_corr,        // 0.02; // minimal correlation value to consider valid
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_kernels,                   // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
					  clt_parameters.kernel_step,    // final int                 kernel_step,
					  clt_parameters.clt_window,     // final int                 window_type,
					  shiftXY,                       // final double [][]         shiftXY, // [port]{shiftX,shiftY}
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.shift_x,        // final double              shiftX, // shift image horizontally (positive - right) - just for testing
					  clt_parameters.shift_y,        // final double              shiftY, // shift image vertically (positive - down)
					  clt_parameters.tileStep,       // final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  clt_parameters.img_dtt.getMcorrSelLY(getNumSensors()), //    final int                 mcorr_sel, // +1 - all, +2 - dia, +4 - sq, +8 - neibs, +16 - hor + 32 - vert

					  clt_parameters.img_dtt.mcorr_comb_width,					// final int                 mcorr_comb_width,  // combined correlation tile width
					  clt_parameters.img_dtt.mcorr_comb_height,					// final int                 mcorr_comb_height, // combined correlation tile full height
					  clt_parameters.img_dtt.mcorr_comb_offset,					// final int                 mcorr_comb_offset, // combined correlation tile height offset: 0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)
					  clt_parameters.img_dtt.mcorr_comb_disp,					// final double              mcorr_comb_disp,   // Combined tile per-pixel disparity for baseline == side of a square
					  
					  clt_parameters.tileX,        // final int                 debug_tileX,
					  clt_parameters.tileY,         // final int                 debug_tileY,
					  threadsMax, // final int                 threadsMax,  // maximal number of threads to launch
					  debugLevel - 2); // final int                 globalDebugLevel)
			  
		  } else {// GPU
			  // First - measure and get 	fcorr_td, fdisp_dist, fpxpy
			  float  [][][][]     fcorr_td =   new float[tilesY][tilesX][][];
			  float  [][][][]     fdisp_dist = new float[tilesY][tilesX][][];
			  float  [][][][]     fpxpy =      new float[tilesY][tilesX][][];
			  final double        gpu_fat_zero = clt_parameters.getGpuFatZero(isMonochrome()); 
			  
			  image_dtt.clt_aberrations_quad_corr_GPU( // measures 2D correlations
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  1,                             // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
					  tile_op,                       // per-tile operation bit codes
					  disparity_array,               // clt_parameters.disparity,     // final double            disparity,
					  fcorr_td,				 		 // final float  [][][][]     corr_td,         // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr pairs
					  null,			  				 //	final float  [][][][]     corr_combo_td,   // [4][tilesY][tilesX][pair][4*64] TD of combo corrs: qud, cross, hor,vert
					  // each of the top elements may be null to skip particular combo type
					  fdisp_dist,                    // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  fpxpy,	                     //	final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  //// Uses quadCLT from gpuQuad			
					  // correlation results - final and partial
					  null,                          // [type][tilesY][tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [type][tilesY][tilesX] should be set by caller
					  // types: 0 - selected correlation (product+offset), 1 - sum
					  null,     			         // clt_corr_partial,// [tilesY][tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
					  // [tilesY][tilesX] should be set by caller
					  // When clt_mismatch is non-zero, no far objects extraction will be attempted
					  null,                          // clt_mismatch,    // [12][tilesY * tilesX] // ***** transpose unapplied ***** ?. null - do not calculate
					  // values in the "main" directions have disparity (*_CM) subtracted, in the perpendicular - as is
					  null,              			 // disparity_map,   // [8][tilesY][tilesX], only [6][] is needed on input or null - do not calculate
					  // last 2 - contrast, avg/ "geometric average)
					  0, // disparity_modes,         // disparity_modes, // bit mask of disparity_map slices to calculate/return
					  null,                          // 	final double [][][][]     texture_tiles,   // compatible with the CPU ones      
					  null, // texture_img,          // texture_img,     // null or [3][] (RGB) or [4][] RGBA
					  null, // texture_woi,          // texture_woi,     // null or generated texture location/size
					  null,                          //  iclt_fimg,       // will return quad images or null to skip, use quadCLT.linearStackToColor 
					  clt_parameters.gpu_corr_scale, //.gpu_corr_scale, //  gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  gpu_fat_zero, // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  clt_parameters.gpu_sigma_r,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_b,    // 0.9, 1.1
					  clt_parameters.gpu_sigma_g,    // 0.6, 0.7
					  clt_parameters.gpu_sigma_m,    //  =       0.4; // 0.7;
					  (isMonochrome()? 1.0 : clt_parameters.gpu_sigma_rb_corr), // final double              gpu_sigma_rb_corr, //  = 0.5; // apply LPF after accumulating R and B correlation before G, monochrome ? 1.0 : gpu_sigma_rb_corr;
					  clt_parameters.gpu_sigma_corr, //  =    0.9;gpu_sigma_corr_m
					  image_dtt.transform_size - 1, // clt_parameters.gpu_corr_rad,   // = transform_size - 1 ?
					  clt_parameters.corr_red,       // +used
					  clt_parameters.corr_blue,      // +used
					  clt_parameters.max_corr_radius,// 3.9;
					  clt_parameters.corr_mode,      // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
					  clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
					  clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
					  clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
					  clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
					  clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					  clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  clt_parameters.clt_window,
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
if (debugLevel < -100) {
			  double [][] lazy_eye_data_dbg = image_dtt.cltMeasureLazyEyeGPU ( // processes 2D correlations, returns d,s lazy eye parameters
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  disparity_array,               // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					  fcorr_td,                      // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr
					  fdisp_dist,                    // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  fpxpy,                         // 						  final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  clt_parameters.gpu_corr_scale, // gpu_corr_scale,  //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  gpu_fat_zero,                  // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // geometryCorrection_main, // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  1, // clt_parameters.tileStep,       // final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  0, // clt_parameters.super_radius,   // final int                 super_radius, // 0 - none, 1 - 3x3, 2 - 5x5,

					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
			  ExtrinsicAdjustment ea_dbg = new ExtrinsicAdjustment (
					  geometryCorrection, // GeometryCorrection gc,
					  1, // int         clusterSize,
					  tilesX, // int         clustersX,
					  tilesY); // int         clustersY)
			  
//			  	double [][] dbg_img = new double[ExtrinsicAdjustment.DATA_TITLES.length][lazy_eye_data_dbg.length];
			  	double [][] dbg_img = new double[ea_dbg.data_titles.length][lazy_eye_data_dbg.length];
			  	for (int i = 0; i < dbg_img.length; i++) {
			  		Arrays.fill(dbg_img[i], Double.NaN);
			  	}
			  	for (int nt = 0; nt < lazy_eye_data_dbg.length; nt++) if (lazy_eye_data_dbg[nt] != null){
				  	for (int i = 0; i < dbg_img.length; i++) {
				  		dbg_img[i][nt] = lazy_eye_data_dbg[nt][i];
				  	}
			  	}
			  
				(new ShowDoubleFloatArrays()).showArrays(
						dbg_img,
						tp.getTilesX(),
						tp.getTilesY(),
						true,
						"LY-everytile",
						ea_dbg.data_titles); //	ExtrinsicAdjustment.DATA_TITLES);
			  
}			  
			  
			  lazy_eye_data = image_dtt.cltMeasureLazyEyeGPU ( // processes 2D correlations, returns d,s lazy eye parameters
					  clt_parameters.img_dtt,        // final ImageDttParameters  imgdtt_params,   // Now just extra correlation parameters, later will include, most others
					  disparity_array,               // final double [][]         disparity_array, // [tilesY][tilesX] - individual per-tile expected disparity
					  fcorr_td,                      // final float  [][][][]     fcorr_td,        // [tilesY][tilesX][pair][4*64] transform domain representation of 6 corr
					  fdisp_dist,                    // final float  [][][][]     fdisp_dist,      // [tilesY][tilesX][cams][4], // disparity derivatives vectors or null
					  fpxpy,                         // 						  final float  [][][][]     fpxpy,           // [tilesY][tilesX][cams][2], tile {pX,pY}
					  clt_parameters.gpu_corr_scale, // gpu_corr_scale,  //  final double              gpu_corr_scale,  //  0.75; // reduce GPU-generated correlation values
					  gpu_fat_zero,                  // final double              gpu_fat_zero,    // clt_parameters.getGpuFatZero(is_mono);absolute == 30.0\
					  geometryCorrection,            // final GeometryCorrection  geometryCorrection,
					  null,                          // geometryCorrection_main, // final GeometryCorrection  geometryCorrection_main, // if not null correct this camera (aux) to the coordinates of the main
					  disparity_corr,                // final double              disparity_corr, // disparity at infinity
					  clt_parameters.tileStep,       // final int                 tileStep, // process tileStep x tileStep cluster of tiles when adjusting lazy eye parameters
					  clt_parameters.super_radius,   // final int                 super_radius, // 0 - none, 1 - 3x3, 2 - 5x5,

					  clt_parameters.tileX,          // final int               debug_tileX,
					  clt_parameters.tileY,          // final int               debug_tileY,
					  threadsMax,
					  debugLevel);
			  
			  
		  }

		  scan.setLazyEyeData(lazy_eye_data);
		  scan.is_measured =   true; // but no disparity map/textures
		  scan.is_combo =      false;
		  scan.resetProcessed();
//		  return scan;
	  }
	  
	  public void gpuResetCorrVector() {
		  if (getGPU() != null) {
			  getGPU().resetGeometryCorrectionVector();
			  getGPU().resetGeometryCorrection(); // testing, just in case, should be removed
		  }
	  }
	  public void gpuResetGeometryCorrection() {
		  if (getGPU() != null) {
			  getGPU().resetGeometryCorrection();
		  }
	  }
	  
	  public void resetBayer() {
		  if (getGPU() != null) {
			  getGPU().resetBayer();
		  }
	  }
	  
		public QuadCLT spawnQuadCLT(
				String              set_name,
				CLTParameters       clt_parameters,
				ColorProcParameters colorProcParameters, //
				int                 threadsMax,
				int                 debugLevel)
		{
			QuadCLT quadCLT = new QuadCLT(this, set_name); //null
			
			quadCLT.restoreFromModel(
					clt_parameters,
					colorProcParameters,
					null,                 // double []    noise_sigma_level,
					-1,                   // noise_variant, // <0 - no-variants, compatible with old code
					null,                 // final QuadCLTCPU     ref_scene, // may be null if scale_fpn <= 0
					threadsMax,
					debugLevel);
			
			return quadCLT;
		}
	  
	  
}
