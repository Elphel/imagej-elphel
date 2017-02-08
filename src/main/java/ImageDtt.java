import java.util.concurrent.atomic.AtomicInteger;

import ij.ImageStack;

/**
 **
 ** ImageDtt - Process images with DTT-based methods
 **
 ** Copyright (C) 2016 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  ImageDtt.java is free software: you can redistribute it and/or modify
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

public class ImageDtt {
	public ImageDtt(){

	}

	public double [][][][] mdctStack(
			final ImageStack                                 imageStack,
			final int                                        subcamera, // 
			final EyesisCorrectionParameters.DCTParameters   dctParameters, //
			final EyesisDCT                                  eyesisDCT,
			final int                                        threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			final int                                        debugLevel,
			final boolean                                    updateStatus) // update status info

	{
	  	  if (imageStack==null) return null;
		  final int imgWidth=imageStack.getWidth();
		  final int nChn=imageStack.getSize();
		  double [][][][] dct_data = new double [nChn][][][];
		  float [] fpixels;
		  int i,chn; //tileX,tileY;
		  /* find number of the green channel - should be called "green", if none - use last */
		  // Extract float pixels from inage stack, convert each to double

		  EyesisDCT.DCTKernels dct_kernels = null;
		  dct_kernels = ((eyesisDCT==null) || (eyesisDCT.kernels==null))?null:eyesisDCT.kernels[subcamera];
		  if (dct_kernels == null){
			  System.out.println("No DCT kernels available for subcamera # "+subcamera);
		  } else if (debugLevel>0){
			  System.out.println("Using DCT kernels for subcamera # "+subcamera);
		  }
//		  if (dctParameters.kernel_chn >=0 ){
//			  dct_kernels = eyesisDCT.kernels[dctParameters.kernel_chn];
//		  }
		  
		  for (chn=0;chn<nChn;chn++) {
			  fpixels= (float[]) imageStack.getPixels(chn+1);
			  double[] dpixels = new double[fpixels.length];
			  for (i = 0; i <fpixels.length;i++) dpixels[i] = fpixels[i];
			  // convert each to DCT tiles
			  dct_data[chn] =lapped_dct(
						dpixels,
						imgWidth,
						dctParameters.dct_size,
						0, //     dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
						dctParameters.dct_window, // final int       window_type,
						chn,
						dct_kernels,
						dctParameters.skip_sym,
						dctParameters.convolve_direct,
						dctParameters.tileX,
						dctParameters.tileY,
						dctParameters.dbg_mode,
						threadsMax,  // maximal number of threads to launch                         
						debugLevel);
		  }
		return dct_data;
	}
	
	public double [][][] lapped_dct(
			final double [] dpixels,
			final int       width,
			final int       dct_size,
			final int       dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
			final int       window_type,
			final int       color,
			final EyesisDCT.DCTKernels dct_kernels,
			final boolean   skip_sym,
			final boolean   convolve_direct, // test feature - convolve directly with the symmetrical kernel
			final int       debug_tileX,
			final int       debug_tileY,
			final int       debug_mode,
			final int       threadsMax,  // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
		final int kernel_margin = 1; //move to parameters?
		final int height=dpixels.length/width;
		final int tilesX=width/dct_size-1;
		final int tilesY=height/dct_size-1;
		final int nTiles=tilesX*tilesY; 
		final double [][][] dct_data = new double[tilesY][tilesX][dct_size*dct_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int i=0; i<dct_data[tileY][tileX].length;i++) dct_data[tileY][tileX][i]= 0.0; // actually not needed, Java initializes arrays
			}
		}
		double [] dc = new double [dct_size*dct_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
		DttRad2 dtt0 = new DttRad2(dct_size);
		dtt0.set_window(window_type);
		final double [] dciii = dtt0.dttt_iii  (dc, dct_size);
		final double [] dciiie = dtt0.dttt_iiie  (dc, 0, dct_size);
		if ((globalDebugLevel > 0) && (color ==2)) {
			double [][]dcx = {dc,dciii,dciiie, dtt0.dttt_ii(dc, dct_size),dtt0.dttt_iie(dc, 0, dct_size)}; 
			showDoubleFloatArrays sdfa_instance0 = new showDoubleFloatArrays(); // just for debugging?
			sdfa_instance0.showArrays(dcx,  dct_size, dct_size, true, "dcx");
		}

		
		if (globalDebugLevel > 0) {
			System.out.println("lapped_dctdc(): width="+width+" height="+height);
		}

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*dct_size * dct_size];
					double [] sym_conv= null;
					if ((dct_kernels != null) && convolve_direct){ // debug feature - directly convolve with symmetrical kernel
						sym_conv = new double[4*dct_size * dct_size];
					}
					double [] tile_folded;
					double [] tile_out; // = new double[dct_size * dct_size];
					int tileY,tileX;
					int n2 = dct_size * 2;
					double [] tile_out_copy = null;
					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						int kernelTileY=0;
						int kernelTileX=0;
						//readDCTKernels() debugLevel = 1 kernels[0].size = 8 kernels[0].img_step = 16 kernels[0].asym_nonzero = 4 nColors = 3 numVert = 123 numHor =  164
						if (dct_kernels != null){ // convolve directly with asym_kernel
							int asym_center = dct_kernels.asym_size/2; // 7 for 15
							kernelTileY = kernel_margin + (tileY * dct_size) / dct_kernels.img_step;
							kernelTileX = kernel_margin + (tileX * dct_size) / dct_kernels.img_step;
							if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
								System.out.println("kernelTileY="+kernelTileY+" kernelTileX="+kernelTileX+" width="+width);
							}
							for (int i = 0; i < n2; i++){
								for (int j = 0; j < n2; j++){
									tile_in[i*n2 + j] = 0.0;
									// convolve list
									int [] asym_indx =   dct_kernels.asym_indx[color][kernelTileY][kernelTileX];
									double [] asym_val = dct_kernels.asym_val[color][kernelTileY][kernelTileX];
									for (int indx = 0; indx < asym_indx.length; indx++){
										int xy = asym_indx[indx];
										if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
											System.out.println("i="+i+" j="+j+" indx="+indx+" xy="+xy);
										}
										if (xy >= 0) {
											int dy = (xy / dct_kernels.asym_size) - asym_center;
											int dx = (xy % dct_kernels.asym_size) - asym_center;
											int y = tileY*dct_size - dy + i;
											int x = tileX*dct_size - dx + j;
											if (y < 0) y &= 1; 
											if (x < 0) x &= 1;
											if (y >= height) y = (height - 2) + (y & 1);
											if (x >= width)  x = (width - 2) +  (x & 1);
											tile_in[i*n2 + j] += asym_val[indx] * dpixels[ y * width + x]; 
											if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
												System.out.println("dy= "+dy+" dx="+dx+" x = "+x+" y="+y+" y*width + x="+(y*width + x));
												System.out.println("asym_val["+indx+"]="+asym_val[indx]+
														"  dpixels["+(y * width + x)+"]="+ dpixels[ y * width + x]+
														"tile_in["+(i*n2 + j)+"]="+tile_in[i*n2 + j]);
											}
										}
									}
								}
							}
							// directly convolve with symmetrical kernel (debug feature
							if ((dct_kernels != null) && convolve_direct){
								double [] dir_sym = dct_kernels.st_direct[color][kernelTileY][kernelTileX];
								double s0 = 0;
								for (int i = 0; i < n2; i++){
									for (int j = 0; j < n2; j++){
										int indx = i*n2+j;
										sym_conv[indx] = 0.0; // dir_sym[0]* tile_in[indx];
										for (int dy = -dct_size +1; dy < dct_size; dy++){
											int ady = (dy>=0)?dy:-dy;
											int sgny = 1;
											int y = i - dy;
											if (y < 0){
												y = -1 -y;
												sgny = -sgny;
											}
											if (y >= n2){
												y = 2*n2 - y -1;
												sgny = -sgny;
											}
											for (int dx = -dct_size +1; dx < dct_size; dx++){
												int adx = (dx >= 0)? dx:-dx;
												int sgn = sgny;
												int x = j - dx;
												if (x < 0){
													x = -1 -x;
													sgn = -sgn;
												}
												if (x >= n2){
													x = 2*n2 - x -1;
													sgn = -sgn;
												}
												sym_conv[indx] += sgn*dir_sym[ady * dct_size + adx] * tile_in[y * n2 + x];
												s0+=dir_sym[ady * dct_size + adx];
												if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2) &&
														(i == dct_size) && (j== dct_size)) {
													System.out.println("i="+i+" j="+j+" dy="+dy+" dx="+dx+" ady="+ady+" adx="+adx+
															" y="+y+" x="+x+" sgny="+sgny+" sgn="+sgn+
															"sym_conv["+indx+"] += "+sgn+"* dir_sym["+(ady * dct_size + adx)+"] * tile_in["+(y * n2 + x)+"] +="+
															sgn+"* "+ dir_sym[ady * dct_size + adx]+" * "+tile_in[y * n2 + x]+" +="+
															(sgn*dir_sym[ady * dct_size + adx] * tile_in[y * n2 + x])+" ="+sym_conv[indx]);
												}

											}
										}
									}
								}


								if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
									//								if ((tileY == debug_tileY) && (tileX == debug_tileX)) {
									double [][] pair = {tile_in, sym_conv};
									sdfa_instance.showArrays(pair,  n2, n2, true, "dconv-X"+tileX+"Y"+tileY+"C"+color);
									sdfa_instance.showArrays(dir_sym,  dct_size, dct_size, "dk-X"+tileX+"Y"+tileY+"C"+color);
									double s1=0,s2=0;
									for (int i = 0; i<tile_in.length; i++){
										s1 +=tile_in[i];
										s2 +=sym_conv[i];
									}
									double s3 = 0.0;
									for (int i=0; i<dct_size;i++){
										for (int j=0; j<dct_size;j++){
											double d = dir_sym[i*dct_size+j];
											if (i > 0) d*=2;
											if (j > 0) d*=2;
											s3+=d;
										}
									}
									System.out.println("s1="+s1+" s2="+s2+" s1/s2="+(s1/s2)+" s0="+s0+" s3="+s3);
								}								
//								tile_in = sym_conv.clone(); 
								System.arraycopy(sym_conv, 0, tile_in, 0, n2*n2);
							}
						} else { // no aberration correction, just copy data
							for (int i = 0; i < n2;i++){
								System.arraycopy(dpixels, (tileY*width+tileX)*dct_size + i*width, tile_in, i*n2, n2);
							}
						}
						tile_folded=dtt.fold_tile(tile_in, dct_size, 0); // DCCT
						tile_out=dtt.dttt_iv  (tile_folded, dct_mode, dct_size);
						if ((dct_kernels != null) && !skip_sym){ // convolve in frequency domain with sym_kernel
							double s0 =0;

							if (debug_mode == 2){
								for (int i=0;i<dct_kernels.st_kernels[color][kernelTileY][kernelTileX].length; i++){
									s0+=dct_kernels.st_kernels[color][kernelTileY][kernelTileX][i];
								}
								s0 = dct_size*dct_size/s0;
							} else if (debug_mode == 3){
								for (int i=0;i<dct_size;i++){
									double scale0 = (i>0)?2.0:1.0; 
									for (int j=0;j<dct_size;j++){
										double scale = scale0*((j>0)?2.0:1.0);
										int indx = i*dct_size+j;
										s0+=scale*dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
									}
								}
								s0 = (2*dct_size-1)*(2*dct_size-1)/s0;
							}else if (debug_mode == 4){
								//dciii								
								for (int i=0;i<dct_kernels.st_kernels[color][kernelTileY][kernelTileX].length; i++){
									s0+=dciii[i]* dct_kernels.st_kernels[color][kernelTileY][kernelTileX][i];
								}
								s0 = dct_size*dct_size/s0;
							} else s0 = 1.0;
							
							for (int i = 0; i < tile_out.length; i++){
								tile_out[i] *= s0;
							}
						}
						
						if ((tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
							tile_out_copy = tile_out.clone();
						}
						
						
						if ((dct_kernels != null) && !skip_sym){ // convolve in frequency domain with sym_kernel
							for (int i = 0; i < tile_out.length; i++){
								tile_out[i] *=dct_kernels.st_kernels[color][kernelTileY][kernelTileX][i];
							}
						}						

						
						if ((dct_kernels!=null) && (tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
							double [][] dbg_tile = {
									dct_kernels.st_direct[color][kernelTileY][kernelTileX],
									dct_kernels.st_kernels[color][kernelTileY][kernelTileX],
									tile_out_copy,
									tile_out};
							if (globalDebugLevel > 0){
								sdfa_instance.showArrays(tile_in,  n2, n2, "tile_in-X"+tileX+"Y"+tileY+"C"+color);
								sdfa_instance.showArrays(dbg_tile,  dct_size, dct_size, true, "dbg-X"+tileX+"Y"+tileY+"C"+color);
								System.out.println("tileY="+tileY+" tileX="+tileX+" kernelTileY="+kernelTileY+" kernelTileX="+kernelTileX);
								double s0=0.0, s1=0.0, s2=0.0, s3=0.0;
								for (int i=0;i<dct_size;i++){
									double scale0 = (i>0)?2.0:1.0; 
									for (int j=0;j<dct_size;j++){
										double scale = scale0*((j>0)?2.0:1.0);
										int indx = i*dct_size+j;
										s0+=scale*dct_kernels.st_direct[color][kernelTileY][kernelTileX][indx];
										s1+=scale*dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
										s2+=      dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
										s3+=dciii[indx]*dct_kernels.st_kernels[color][kernelTileY][kernelTileX][indx];
									}
								}
								System.out.println("s0="+s0+" s1="+s1+" s2="+s2+" s3="+s3);
							}
						}
						System.arraycopy(tile_out, 0, dct_data[tileY][tileX], 0, tile_out.length);
					}
				}
			};
		}		      
		startAndJoin(threads);
		return dct_data;
	}
	
	// extract DCT transformed parameters in linescan order (for visualization)
	public double [] lapped_dct_dbg(
			final double [][][] dct_data,
			final int           threadsMax,     // maximal number of threads to launch                         
			final int           globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] dct_data_out = new double[tilesY*tilesX*dct_len];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int i=0; i<dct_data_out.length;i++) dct_data_out[i]= 0;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < dct_size;i++){
							System.arraycopy(dct_data[tileY][tileX], dct_size* i, dct_data_out, ((tileY*dct_size + i) *tilesX + tileX)*dct_size , dct_size);
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return dct_data_out;
	}
	
	public void dct_lpf(
			final double sigma,
			final double [][][] dct_data,
			final int       threadsMax,     // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] filter_direct= new double[dct_len];
		if (sigma == 0) {
			filter_direct[0] = 1.0; 
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0; 
		} else {
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					filter_direct[i*dct_size+j] = Math.exp(-(i*i+j*j)/(2*sigma));
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	filter_direct[i*dct_size+j];
				d*=Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}
		
		if (globalDebugLevel > 0) {
			for (int i=0; i<filter_direct.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter_direct[i]); 
			}
		}
		DttRad2 dtt = new DttRad2(dct_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		final double [] dbg_filter= dtt.dttt_ii(filter);
		
//		for (int i=0; i < filter.length;i++) filter[i] *= dct_size;  
		for (int i=0; i < filter.length;i++) filter[i] *= 2*dct_size;  
		
		if (globalDebugLevel > 0) {
			for (int i=0; i<filter.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter[i]); 
			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filter_direct,filter,dbg_filter};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filter_lpf");
		}
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < filter.length; i++){
							dct_data[tileY][tileX][i] *= filter[i];
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
	}
	
	public double [][][][] dct_color_convert(
			final double [][][][] dct_data,
			final double kr,
			final double kb,
			final double sigma_rb,        // blur of channels 0,1 (r,b) in addition to 2 (g)
			final double sigma_y,         // blur of Y from G
			final double sigma_color,     // blur of Pr, Pb in addition to Y
			final int       threadsMax,     // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data[0].length;
		final int tilesX=dct_data[0][0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [][][][] yPrPb = new double [3][tilesY][tilesX][dct_len];
		final double [][][] filters = new double [3][3][dct_len];
		final double kg = 1.0 - kr - kb;
		final double [][] filters_proto_direct = new double[3][dct_len];
		final double [][] filters_proto = new double[3][];
		System.out.println("dct_color_convert(): kr="+kr+" kg="+kg+" kb="+kb);
		final double [] sigmas = {sigma_rb,sigma_y,sigma_color};
		double [] norm_sym_weights = new double [dct_size*dct_size];
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				norm_sym_weights[i*dct_size+j] = d;
			}
		}
		
		for (int n = 0; n<3; n++) {
			
			double s = 0.0;
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					double d;
					if (sigmas[n] == 0.0)   d = ((i == 0) && (j==0))? 1.0:0.0;
					else                    d = Math.exp(-(i*i+j*j)/(2*sigmas[n]));
					filters_proto_direct[n][i*dct_size+j] = d;
				}
				
			}
			for (int i = 0; i< dct_len; i++){
				s += norm_sym_weights[i]*filters_proto_direct[n][i];
			}
			
			if (globalDebugLevel>0) System.out.println("dct_color_convert(): sigmas["+n+"]="+sigmas[n]+", sum="+s);
			for (int i = 0; i < dct_len; i++){
				filters_proto_direct[n][i] /=s;
			}
		}
		
		DttRad2 dtt = new DttRad2(dct_size);
		for (int i = 0; i < filters_proto.length; i++){
			filters_proto[i] = dtt.dttt_iiie(filters_proto_direct[i]);
			if (globalDebugLevel > 0)  System.out.println("filters_proto.length="+filters_proto.length+" filters_proto["+i+"].length="+filters_proto[i].length+" dct_len="+dct_len+" dct_size="+dct_size);
			for (int j=0; j < dct_len; j++) filters_proto[i][j] *= 2*dct_size;  

		}
		if (globalDebugLevel > 0) {
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filters_proto_direct[0],filters_proto_direct[1],filters_proto_direct[2],filters_proto[0],filters_proto[1],filters_proto[2]};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filters_proto");
		}

		double [][] coeff_arr ={
				{ kr,              kb,               kg             },  // Y = R*Kr+G*Hg+B*Kb
				{ 0.5,            -kb/(2.0*(1-kr)), -kg/(2.0*(1-kr))},  // Pr =  R* 0.5  - G* Kg/(2.0*(1-Kr)) - B *Kb/(2.0*(1-Kr))
				{-kr/(2.0*(1-kb)), 0.5,             -kg/(2.0*(1-kb))}}; // Pb =  B* 0.5  - G* Kg/(2.0*(1-Kb)) - R *Kr/(2.0*(1-Kb))
		for (int k = 0; k < dct_len; k++){
			for (int i = 0; i < coeff_arr.length; i++){
				for (int j = 0; j < coeff_arr.length; j++){
					filters[i][j][k] = coeff_arr[i][j]* filters_proto[1][k];      // minimal blur - for all sigma_y
					if (i > 0){
						filters[i][j][k] *= filters_proto[2][k]; // for Pr, Pb sigma_color
					}
					if (j <2){ // all but green 
						filters[i][j][k] *= filters_proto[0][k]; // for R,B sigma_rb
					}
					
				}
			}
		}
		if (globalDebugLevel > 0) {
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			double [][] ff = {
					filters[0][0], filters[0][1], filters[0][2],
					filters[1][0], filters[1][1], filters[1][2],
					filters[2][0], filters[2][1], filters[2][2]};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filters");
		}
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < filters.length; i++){
							for (int k = 0; k <dct_len; k++){
								yPrPb[i][tileY][tileX][k]=0.0;
								for (int j = 0; j < filters[i].length; j++){
									yPrPb[i][tileY][tileX][k] += filters[i][j][k] * dct_data[j][tileY][tileX][k];
								}
							}
							
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return yPrPb;
	}
	

	
	
	
	
	public double [] lapped_idct(
//			final double [][][] dctdc_data,  // array [tilesY][tilesX][dct_size*dct_size+1] - last element is DC value  
			final double [][][] dct_data,  // array [tilesY][tilesX][dct_size*dct_size]  
			final int       dct_size,
			final int       window_type,
			final int       threadsMax,  // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
//		final int tilesX=dct_width/dct_size;
//		final int tilesY=dct_data.length/(dct_width*dct_size);
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;

		final int width=  (tilesX+1)*dct_size;
		final int height= (tilesY+1)*dct_size;
		if (globalDebugLevel > 0) {
			System.out.println("lapped_idct():tilesX=   "+tilesX);
			System.out.println("lapped_idct():tilesY=   "+tilesY);
			System.out.println("lapped_idct():width=    "+width);
			System.out.println("lapped_idct():height=   "+height);
		}
		final double [] dpixels = new double[width*height];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger nser = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int i=0; i<dpixels.length;i++) dpixels[i]= 0;
		for (int n=0; n<4; n++){
			nser.set(n);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						DttRad2 dtt = new DttRad2(dct_size);
						dtt.set_window(window_type);
						double [] tile_in = new double[dct_size * dct_size];
						double [] tile_dct; // = new double[dct_size * dct_size];
						double [] tile_out; //  = new double[4*dct_size * dct_size];
						int tileY,tileX;
						int n2 = dct_size * 2;
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];
							System.arraycopy(dct_data[tileY][tileX], 0, tile_in, 0, tile_in.length);
							tile_dct=dtt.dttt_iv  (tile_in, 0, dct_size);
							tile_out=dtt.unfold_tile(tile_dct, dct_size, 0); // mpode=0 - DCCT
							for (int i = 0; i < n2;i++){
								int start_line = ((tileY*dct_size + i) *(tilesX+1) + tileX)*dct_size; 
								for (int j = 0; j<n2;j++) {
									dpixels[start_line + j] += tile_out[n2 * i + j]; //  +1.0; 
								}
							}
						}
					}
				};
			}		      
			startAndJoin(threads);
		}
		return dpixels;
	}
	
// perform 2d clt, result is [tileY][tileX][cc_sc_cs_ss][index_in_tile]
	public double [][][][] clt_2d(
			final double [] dpixels,
			final int       width,
			final int       dct_size,
			final int       window_type,
			final int       shiftX, // shift image horizontally (positive - right)
			final int       shiftY, // shift image vertically (positive - down)
			final int       debug_tileX,
			final int       debug_tileY,
			final int       debug_mode,
			final int       threadsMax,  // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
		final int height=dpixels.length/width;
		final int tilesX=width/dct_size-1;
		final int tilesY=height/dct_size-1;
		final int nTiles=tilesX*tilesY; 
		final double [][][][] dct_data = new double[tilesY][tilesX][4][dct_size*dct_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int dct_mode = 0; dct_mode < dct_data[tileY][tileX].length; dct_mode++){
					for (int i=0; i<dct_data[tileY][tileX][dct_mode].length;i++) {
						dct_data[tileY][tileX][dct_mode][i]= 0.0; // actually not needed, Java initializes arrays
					}
				}
			}
		}
		double [] dc = new double [dct_size*dct_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
		DttRad2 dtt0 = new DttRad2(dct_size);
		dtt0.set_window(window_type);
		if (globalDebugLevel > 0) {
			System.out.println("clt_2d(): width="+width+" height="+height+" dct_size="+dct_size+
					" debug_tileX="+debug_tileX+" debug_tileY="+debug_tileY+" globalDebugLevel="+globalDebugLevel);
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*dct_size * dct_size];
					double [][] tile_folded = new double[4][];
					double [][] tile_out =    new double[4][]; // = new double[dct_size * dct_size];
					int tileY,tileX;
					int n2 = dct_size * 2;
//					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						if ((shiftX == 0) && (shiftY == 0)){
							for (int i = 0; i < n2;i++){
								System.arraycopy(dpixels, (tileY*width+tileX)*dct_size + i*width, tile_in, i*n2, n2);
							}
						} else {
							int x0 = tileX * dct_size - shiftX;
							if      (x0 < 0)             x0 = 0; // first/last will be incorrect
							else if (x0 >= (width - n2)) x0 = width - n2;  
							for (int i = 0; i < n2;i++){
								int y0 = tileY * dct_size + i - shiftY;
								if      (y0 < 0)       y0 = 0;
								else if (y0 >= height) y0 = height -1;
								System.arraycopy(dpixels, y0 * width+ x0, tile_in, i*n2, n2);
							}
						}
						for (int dct_mode = 0; dct_mode <4; dct_mode++) {
							tile_folded[dct_mode] = dtt.fold_tile(tile_in, dct_size, dct_mode); // DCCT, DSCT, DCST, DSST
							if ((debug_mode & 1) != 0) {
								tile_out[dct_mode] = tile_folded[dct_mode];
							} else {
								tile_out[dct_mode] =    dtt.dttt_iv  (tile_folded[dct_mode], dct_mode, dct_size);
							}
							System.arraycopy(tile_out[dct_mode], 0, dct_data[tileY][tileX][dct_mode], 0, tile_out[dct_mode].length);
						}
						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
							sdfa_instance.showArrays(tile_in,  n2, n2, "tile_in_x"+tileX+"_y"+tileY);
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(tile_folded,  dct_size, dct_size, true, "folded_x"+tileX+"_y"+tileY, titles);
							if (globalDebugLevel > 0) {
								sdfa_instance.showArrays(tile_out,     dct_size, dct_size, true, "clt_x"+tileX+"_y"+tileY, titles);
							}
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return dct_data;
	}
	
	public double [] iclt_2d(
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]  
			final int             dct_size,
			final int             window_type,
			final int             debug_mask, // which transforms to combine
			final int             debug_mode, // skip idct - just unfold
			final int             threadsMax,  // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;

		final int width=  (tilesX+1)*dct_size;
		final int height= (tilesY+1)*dct_size;
		final double debug_scale = 1.0 /((debug_mask & 1) + ((debug_mask >> 1) & 1) + ((debug_mask >> 2) & 1) + ((debug_mask >> 3) & 1));
		if (globalDebugLevel > 0) {
			System.out.println("iclt_2d():tilesX=        "+tilesX);
			System.out.println("iclt_2d():tilesY=        "+tilesY);
			System.out.println("iclt_2d():width=         "+width);
			System.out.println("iclt_2d():height=        "+height);
			System.out.println("iclt_2d():debug_mask=    "+debug_mask);
			System.out.println("iclt_2d():debug_scale=   "+debug_scale);
		}
		final double [] dpixels = new double[width*height];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger nser = new AtomicInteger(0);
		final int [][][] tiles_list = new int[4][][];
		for (int n=0; n<4; n++){
			int nx = (tilesX + 1 - (n &1)) / 2;
			int ny = (tilesY + 1 - ((n>>1) & 1)) / 2;
			tiles_list[n] = new int [nx*ny][2];
			int indx = 0;
			for (int i = 0;i < ny; i++) for (int j = 0; j < nx; j++){
				tiles_list[n][indx][0]=2*j+(n &1);
				tiles_list[n][indx++][1]=2*i+((n>>1) & 1);
			}
		}
		for (int i=0; i<dpixels.length;i++) dpixels[i]= 0;
		for (int n=0; n<4; n++){
			nser.set(n);
			ai.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					public void run() {
						DttRad2 dtt = new DttRad2(dct_size);
						dtt.set_window(window_type);
						double [] tile_in =   new double [dct_size * dct_size];
						double [] tile_dct;
						double [] tile_mdct;
						int tileY,tileX;
						int n2 = dct_size * 2;
						for (int nTile = ai.getAndIncrement(); nTile < tiles_list[nser.get()].length; nTile = ai.getAndIncrement()) {
							tileX = tiles_list[nser.get()][nTile][0];
							tileY = tiles_list[nser.get()][nTile][1];
							for (int dct_mode = 0; dct_mode < 4; dct_mode++) if (((1 << dct_mode) & debug_mask) != 0) {
								System.arraycopy(dct_data[tileY][tileX][dct_mode], 0, tile_in, 0, tile_in.length);
								if ((debug_mode & 1) != 0) {
									tile_dct = tile_in;
								} else {
									// IDCT-IV should be in reversed order: CC->CC, SC->CS, CS->SC, SS->SS 
									int idct_mode = ((dct_mode << 1) & 2) | ((dct_mode >> 1) & 1); 
									tile_dct = dtt.dttt_iv  (tile_in, idct_mode, dct_size);
								}
								tile_mdct = dtt.unfold_tile(tile_dct, dct_size, dct_mode); // mode=0 - DCCT
								for (int i = 0; i < n2;i++){
									int start_line = ((tileY*dct_size + i) *(tilesX+1) + tileX)*dct_size; 
									for (int j = 0; j<n2;j++) {
										dpixels[start_line + j] += debug_scale * tile_mdct[n2 * i + j]; // add (cc+sc+cs+ss)/4 
									}
								}
							}
						}
					}
				};
			}		      
			startAndJoin(threads);
		}
		return dpixels;
	}

	public double [][][][] clt_shiftXY(
			final double [][][][] dct_data,  // array [tilesY][tilesX][4][dct_size*dct_size]  
			final int             dct_size,
			final double          shiftX,
			final double          shiftY,
			final int             dbg_swap_mode,
			final int             threadsMax,  // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles = tilesY* tilesX; 
		if (globalDebugLevel > 0) {
			System.out.println("clt_shift():tilesX=        "+tilesX);
			System.out.println("clt_shift():tilesY=        "+tilesY);
			System.out.println("clt_shift():shiftX=        "+shiftX);
			System.out.println("clt_shift():shiftY=        "+shiftY);
		}
		final double [] cos_hor =  new double [dct_size*dct_size];
		final double [] sin_hor =  new double [dct_size*dct_size];
		final double [] cos_vert = new double [dct_size*dct_size];
		final double [] sin_vert = new double [dct_size*dct_size];
		for (int i = 0; i < dct_size; i++){
			double ch = Math.cos((i+0.5)*Math.PI*shiftX/dct_size);
			double sh = Math.sin((i+0.5)*Math.PI*shiftX/dct_size);
			double cv = Math.cos((i+0.5)*Math.PI*shiftY/dct_size);
			double sv = Math.sin((i+0.5)*Math.PI*shiftY/dct_size);
			for (int j = 0; j < dct_size; j++){
				int iv = dct_size * j + i; // 2d DTT results are stored transposed! 
				int ih = dct_size * i + j; 
				cos_hor[ih] = ch; 
				sin_hor[ih] = sh; 
				cos_vert[iv] = cv; 
				sin_vert[iv] = sv; 
			}
			
		}

		if (globalDebugLevel > 0){
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			String [] titles = {"cos_hor","sin_hor","cos_vert","sin_vert"};
			double [][] cs_dbg = {cos_hor, sin_hor, cos_vert, sin_vert};
			sdfa_instance.showArrays(cs_dbg,  dct_size, dct_size, true, "shift_cos_sin", titles);
		}
		
		final double [][][][] rslt = new double[dct_data.length][dct_data[0].length][dct_data[0][0].length][dct_data[0][0][0].length];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						// Horizontal shift CLT tiled data is stored in transposed way (horizontal - Y, vertical X) 
						for (int i = 0; i < cos_hor.length; i++) {
							rslt[tileY][tileX][0][i] = dct_data[tileY][tileX][0][i] * cos_hor[i] - dct_data[tileY][tileX][1][i] * sin_hor[i];
							rslt[tileY][tileX][1][i] = dct_data[tileY][tileX][1][i] * cos_hor[i] + dct_data[tileY][tileX][0][i] * sin_hor[i] ;

							rslt[tileY][tileX][2][i] = dct_data[tileY][tileX][2][i] * cos_hor[i]  - dct_data[tileY][tileX][3][i] * sin_hor[i];
							rslt[tileY][tileX][3][i] = dct_data[tileY][tileX][3][i] * cos_hor[i]  + dct_data[tileY][tileX][2][i] * sin_hor[i] ;
						}
						// Vertical shift (in-place)
						for (int i = 0; i < cos_hor.length; i++) {
							double tmp =               rslt[tileY][tileX][0][i] * cos_vert[i] - rslt[tileY][tileX][2][i] * sin_vert[i];
							rslt[tileY][tileX][2][i] = rslt[tileY][tileX][2][i] * cos_vert[i] + rslt[tileY][tileX][0][i] * sin_vert[i];
							rslt[tileY][tileX][0][i] = tmp;

							tmp =                      rslt[tileY][tileX][1][i] * cos_vert[i] - rslt[tileY][tileX][3][i] * sin_vert[i];
							rslt[tileY][tileX][3][i] = rslt[tileY][tileX][3][i] * cos_vert[i] + rslt[tileY][tileX][1][i] * sin_vert[i];
							rslt[tileY][tileX][1][i] = tmp;
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return rslt;
	}

	public double [][][][] clt_correlate(
			final double [][][][] data1,  // array [tilesY][tilesX][4][dct_size*dct_size]  
			final double [][][][] data2,  // array [tilesY][tilesX][4][dct_size*dct_size]  
			final int             dct_size,
			final double          fat_zero,    // add to denominator to modify phase correlation (same units as data1, data2)
			final int             debug_tileX,
			final int             debug_tileY,
			final int             threadsMax,  // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=(data1.length > data2.length)?data2.length:data1.length;
		final int tilesX=(data1[0].length > data2[0].length)?data2[0].length:data1[0].length;
		final int nTiles = tilesY* tilesX;
		if (globalDebugLevel > 0) {
			System.out.println("clt_shift():tilesX= "+tilesX);
			System.out.println("clt_shift():tilesY= "+tilesY);
		}
		/* Direct matrix Z1: X2 ~= Z1 * Shift   
		 * {{+cc  -sc  -cs  +ss},
		 *  {+sc  +cc  -ss  -cs},
		 *  {+cs  -ss  +cc  -sc},
		 *  {+ss  +cs  +sc  +cc}}
		 *  
		 * T= transp({cc, sc, cs, ss}) 
		 */
		/*
		final int [][] zi = 
			{{ 0, -1, -2,  3},
			 { 1,  0, -3, -2},
			 { 2, -3,  0, -1},
			 { 3,  2,  1,  0}};
		*/
		final int [][] zi = 
			{{ 0,  1,  2,  3},
			 {-1,  0, -3,  2},
			 {-2, -3,  0,  1},
			 { 3, -2, -1,  0}};
		
		final int dct_len = dct_size * dct_size;
		final double [][][][] rslt = new double[tilesY][tilesX][4][dct_len];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < dct_len; i++) {
							double s1 = 0.0, s2=0.0;
							for (int n = 0; n< 4; n++){
								s1+=data1[tileY][tileX][n][i] * data1[tileY][tileX][n][i];
								s2+=data2[tileY][tileX][n][i] * data2[tileY][tileX][n][i];
							}
							double scale = 1.0 / (Math.sqrt(s1*s2) + fat_zero*fat_zero); // squared to match units
							for (int n = 0; n<4; n++){
								/*
								if (
										(tileY >= rslt.length) ||
										(tileX >= rslt[tileY].length) ||
										(n >= rslt[tileY][tileX].length) ||
										(i >= rslt[tileY][tileX][n].length)) {
									
									System.out.println("===== tileY="+tileY+" ("+tilesY+") tileX="+tileX+" ("+tilesX+") n="+n+" i="+i);
								
									System.out.println(
											" rslt.length="+rslt.length+
											" rslt.length[tileY]="+rslt[tileY].length+
											" rslt.length[tileY][tileX]="+rslt[tileY][tileX].length+
											" rslt.length[tileY][tileX][n]="+rslt[tileY][tileX][n].length);
								System.out.println("===== tileY="+tileY+" ("+tilesY+") tileX="+tileX+" ("+tilesX+") n="+n+" i="+i);
								}
								*/
								rslt[tileY][tileX][n][i] = 0;
								for (int k=0; k<4; k++){
									if (zi[n][k] < 0)
										rslt[tileY][tileX][n][i] -= 
											data1[tileY][tileX][-zi[n][k]][i] * data2[tileY][tileX][k][i];
									else
										rslt[tileY][tileX][n][i] += 
										data1[tileY][tileX][zi[n][k]][i] * data2[tileY][tileX][k][i];
								}
								rslt[tileY][tileX][n][i] *= scale;
							}
						}
						if ((globalDebugLevel > 0) && (debug_tileX == tileX) && (debug_tileY == tileY)) {
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
							String [] titles = {"CC","SC","CS","SS"};
							sdfa_instance.showArrays(data1[tileY][tileX], dct_size, dct_size, true, "data1_x"+tileX+"_y"+tileY, titles);
							sdfa_instance.showArrays(data2[tileY][tileX], dct_size, dct_size, true, "data2_x"+tileX+"_y"+tileY, titles);
							sdfa_instance.showArrays(rslt[tileY][tileX],  dct_size, dct_size, true, "rslt_x"+ tileX+"_y"+tileY, titles);
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return rslt;
	}

	public void clt_lpf(
			final double          sigma,
			final double [][][][] clt_data,
			final int             threadsMax,     // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=clt_data.length;
		final int tilesX=clt_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(clt_data[0][0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] filter_direct= new double[dct_len];
		if (sigma == 0) {
			filter_direct[0] = 1.0; 
			for (int i= 1; i<filter_direct.length;i++) filter_direct[i] =0; 
		} else {
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					filter_direct[i*dct_size+j] = Math.exp(-(i*i+j*j)/(2*sigma));
				}
			}
		}
		// normalize
		double sum = 0;
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	filter_direct[i*dct_size+j];
				d*=Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				sum +=d;
			}
		}
		for (int i = 0; i<filter_direct.length; i++){
			filter_direct[i] /= sum;
		}
		
		if (globalDebugLevel > 0) {
			for (int i=0; i<filter_direct.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter_direct[i]); 
			}
		}
		DttRad2 dtt = new DttRad2(dct_size);
		final double [] filter= dtt.dttt_iiie(filter_direct);
		final double [] dbg_filter= dtt.dttt_ii(filter);
		for (int i=0; i < filter.length;i++) filter[i] *= 2*dct_size;  
		
		if (globalDebugLevel > 0) {
			for (int i=0; i<filter.length;i++){
				System.out.println("dct_lpf_psf() "+i+": "+filter[i]); 
			}
			showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
			double [][] ff = {filter_direct,filter,dbg_filter};
			sdfa_instance.showArrays(ff,  dct_size,dct_size, true, "filter_lpf");
		}
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int n = 0; n < 4; n++){
							for (int i = 0; i < filter.length; i++){
								clt_data[tileY][tileX][n][i] *= filter[i];
							}
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
	}
	
	public void clt_dtt2( // transform dcct2, dsct2, dcst2, dsst2
			final double [][][][] data,
			final boolean         transpose, // when doing inverse transform, the data comes in transposed form, so CS <->SC
			final int             threadsMax,     // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=data.length;
		final int tilesX=data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(data[0][0][0].length));
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int quadrant = 0; quadrant < 4; quadrant++){
							int mode = transpose ? (((quadrant << 1) & 2) | ((quadrant >> 1) & 1)) : quadrant;
							data[tileY][tileX][quadrant] = dtt.dttt_iie(data[tileY][tileX][quadrant], mode, dct_size);
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
	}
	
	public double [][][] clt_corr_quad( // combine 4 correlation quadrants after DTT2
			final double [][][][] data,
			final int             threadsMax,     // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=data.length;
		final int tilesX=data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(data[0][0][0].length));
		final int rslt_size=dct_size*2-1;
		
		final double [][][] rslt = new double[tilesY][tilesX][rslt_size*rslt_size];
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					double scale = 0.25;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						rslt[tileY][tileX][rslt_size*dct_size - dct_size] = scale * data[tileY][tileX][0][0]; // center
						for (int j = 1; j < dct_size; j++) { //  for i == 0
							rslt[tileY][tileX][rslt_size*dct_size - dct_size + j] = scale * (data[tileY][tileX][0][j] + data[tileY][tileX][1][j-1]); 
							rslt[tileY][tileX][rslt_size*dct_size - dct_size - j] = scale * (data[tileY][tileX][0][j] - data[tileY][tileX][1][j-1]); 
						}						
						for (int i = 1; i < dct_size; i++) {
							rslt[tileY][tileX][rslt_size*(dct_size + i) - dct_size] =
									scale * (data[tileY][tileX][0][i*dct_size] + data[tileY][tileX][2][(i-1)*dct_size]); 
							rslt[tileY][tileX][rslt_size*(dct_size - i) - dct_size] =
									scale * (data[tileY][tileX][0][i*dct_size] - data[tileY][tileX][2][(i-1)*dct_size]); 
							for (int j = 1; j < dct_size; j++) {
								rslt[tileY][tileX][rslt_size*(dct_size + i) - dct_size + j] =
										scale * (data[tileY][tileX][0][i*    dct_size + j] + 
												 data[tileY][tileX][1][i*    dct_size + j - 1] +
												 data[tileY][tileX][2][(i-1)*dct_size + j] +
												 data[tileY][tileX][3][(i-1)*dct_size + j - 1]); 
								
								rslt[tileY][tileX][rslt_size*(dct_size + i) - dct_size - j] =
										scale * ( data[tileY][tileX][0][i*    dct_size + j] + 
												 -data[tileY][tileX][1][i*    dct_size + j - 1] +
												  data[tileY][tileX][2][(i-1)*dct_size + j] +
												 -data[tileY][tileX][3][(i-1)*dct_size + j - 1]); 
								rslt[tileY][tileX][rslt_size*(dct_size - i) - dct_size + j] =
										scale * (data[tileY][tileX][0][i*    dct_size + j] + 
												 data[tileY][tileX][1][i*    dct_size + j - 1] +
												 -data[tileY][tileX][2][(i-1)*dct_size + j] +
												 -data[tileY][tileX][3][(i-1)*dct_size + j - 1]);
								rslt[tileY][tileX][rslt_size*(dct_size - i) - dct_size - j] =
										scale * (data[tileY][tileX][0][i*    dct_size + j] + 
												 -data[tileY][tileX][1][i*    dct_size + j - 1] +
												 -data[tileY][tileX][2][(i-1)*dct_size + j] +
												 data[tileY][tileX][3][(i-1)*dct_size + j - 1]); 
							}
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return rslt;
	}
	
	// extract correlation result  in linescan order (for visualization)
	public double [] corr_dbg(
			final double [][][] corr_data,
			final int           threadsMax,     // maximal number of threads to launch                         
			final int           globalDebugLevel)
	{
		final int tilesY=corr_data.length;
		final int tilesX=corr_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int corr_size = (int) Math.round(Math.sqrt(corr_data[0][0].length));
		
		final int tile_size = corr_size+1;
		final int corr_len = corr_size*corr_size;
		
		final double [] corr_data_out = new double[tilesY*tilesX*tile_size*tile_size];
		
		System.out.println("corr_dbg(): tilesY="+tilesY+", tilesX="+tilesX+", corr_size="+corr_size+", corr_len="+corr_len);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int i=0; i<corr_data_out.length;i++) corr_data_out[i]= 0;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int i = 0; i < corr_size;i++){
							System.arraycopy(corr_data[tileY][tileX], corr_size* i, corr_data_out, ((tileY*tile_size + i) *tilesX + tileX)*tile_size , corr_size);
							corr_data_out[((tileY*tile_size + i) *tilesX + tileX)*tile_size+corr_size] = (i & 1) - 0.5;
						}
						for (int i = 0; i < tile_size; i++){
							corr_data_out[((tileY*tile_size + corr_size) *tilesX + tileX)*tile_size+i] = (i & 1) - 0.5;
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return corr_data_out;
	}
	
	
	
	public double [][][][][] cltStack(
			final ImageStack                                 imageStack,
			final int                                        subcamera, // 
			final EyesisCorrectionParameters.CLTParameters   cltParameters, //
			final int                                        shiftX, // shift image horizontally (positive - right)
			final int                                        shiftY, // shift image vertically (positive - down)
			final int                                        threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			final int                                        debugLevel,
			final boolean                                    updateStatus) // update status info

	{
	  	  if (imageStack==null) return null;
		  final int imgWidth=imageStack.getWidth();
		  final int nChn=imageStack.getSize();
		  double [][][][][] dct_data = new double [nChn][][][][];
		  float [] fpixels;
		  int i,chn; //tileX,tileY;
		  /* find number of the green channel - should be called "green", if none - use last */
		  // Extract float pixels from inage stack, convert each to double

//		  EyesisDCT.DCTKernels dct_kernels = null;
//		  dct_kernels = eyesisDCT.kernels[subcamera];
//		  if (dct_kernels == null){
//			  System.out.println("No DCT kernels available for subcamera # "+subcamera);
//		  } else if (debugLevel>0){
//			  System.out.println("Using DCT kernels for subcamera # "+subcamera);
//		  }
//		  if (dctParameters.kernel_chn >=0 ){
//			  dct_kernels = eyesisDCT.kernels[dctParameters.kernel_chn];
//		  }
		  
		  for (chn=0;chn<nChn;chn++) {
			  fpixels= (float[]) imageStack.getPixels(chn+1);
			  double[] dpixels = new double[fpixels.length];
			  for (i = 0; i <fpixels.length;i++) dpixels[i] = fpixels[i];
			  // convert each to DCT tiles
			  dct_data[chn] = clt_2d(
						dpixels,
						imgWidth,
						cltParameters.transform_size,
						cltParameters.clt_window,
						shiftX,
						shiftY,
						cltParameters.tileX,    //       debug_tileX,
						cltParameters.tileY,    //       debug_tileY,
						cltParameters.dbg_mode, //       debug_mode,
						threadsMax,  // maximal number of threads to launch                         
						debugLevel);
		  }
		return dct_data;
	}
	
	
	
	
	// extract DCT transformed parameters in linescan order (for visualization)
	public double [][] clt_dbg(
			final double [][][][] dct_data,
			final int             threadsMax,     // maximal number of threads to launch                         
			final int             globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [][] dct_data_out = new double[4][tilesY*tilesX*dct_len];
		
		System.out.println("clt_dbg(): tilesY="+tilesY+", tilesX="+tilesX+", dct_size="+dct_size+", dct_len="+dct_len+", dct_data_out[0].length="+dct_data_out[0].length);
		
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int n=0; n<dct_data_out.length;n++) for (int i=0; i<dct_data_out[n].length;i++) dct_data_out[n][i]= 0;

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						for (int n=0; n<dct_data_out.length;n++) {
							for (int i = 0; i < dct_size;i++){
								System.arraycopy(dct_data[tileY][tileX][n], dct_size* i, dct_data_out[n], ((tileY*dct_size + i) *tilesX + tileX)*dct_size , dct_size);
							}
						}
					}
				}
			};
		}		      
		startAndJoin(threads);
		return dct_data_out;
	}

	void clt_convert_double_kernel( // converts double resolution kernel
			double []   src_kernel, //
			double []   dst_kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) +2 size - kernel and dx, dy to the nearest 1/2 pixels
			int src_size, // 64
			int dtt_size) // 8
	{
		
		int [] indices = {0,-src_size,-1,1,src_size,-src_size-1,-src_size+1,src_size-1,src_size+1};
		double [] weights = {0.25,0.125,0.125,0.125,0.125,0.0625,0.0625,0.0625,0.0625};
		int src_center = src_size / 2; // 32
		// Find center
		double sx=0.0, sy = 0.0, s = 0.0;
		int indx = 0;
		for (int i= -src_center; i < src_center; i++){
			for (int j = -src_center; j < src_center; j++){
				double d = src_kernel[indx++];
				sx+= j*d;
				sy+= i*d;
				s += d;
			}
		}
		int src_x = (int) Math.round(sx / s) + src_center;
		int src_y = (int) Math.round(sy / s) + src_center;
		// make sure selected area (2*dst_size-1) * (2*dst_size-1) fits into src_kernel, move center if not
		if      (src_x < 2 * dtt_size)             src_x = 2 * dtt_size - 1; // 15
		else if (src_x > (src_size - 2* dtt_size)) src_x = src_size - 2* dtt_size;
		
		if      (src_y < 2 * dtt_size)             src_y = 2 * dtt_size - 1; // 15
		else if (src_y > (src_size - 2* dtt_size)) src_y = src_size - 2* dtt_size;
		indx = 0;
		// downscale, copy
		for (int i = -dtt_size + 1; i < dtt_size; i++){
			int src_i = (src_y + 2 * i) * src_size  + src_x; 
			for (int j = -dtt_size + 1; j < dtt_size; j++){
				double d = 0.0;
				for (int k = 0; k < indices.length; k++){
					d += weights[k]*src_kernel[src_i + 2 * j + indices[k]];
				}
				dst_kernel[indx++] = d;
			}			
		}
		dst_kernel[indx++] = 0.5*(src_x - src_center);
		dst_kernel[indx++] = 0.5*(src_y - src_center);
	}

	void clt_normalize_kernel( // 
			double []   kernel, // should be (2*dtt_size-1) * (2*dtt_size-1) +2 size (last 2 are not modified)
			double []   window, // normalizes result kernel * window to have sum of elements == 1.0 
			final int dtt_size) // 8
	{
		double s = 0.0;
		int indx = 0;
		for (int i = -dtt_size + 1; i < dtt_size; i++){
			int ai = (i < 0)? -i: i;
			for (int j = -dtt_size + 1; j < dtt_size; j++){
				int aj = (j < 0)? -j: j;
				s += kernel[indx++] * window[ai*dtt_size+aj];
			}
		}
		s = 1.0/s;
		indx = 0;
		for (int i = 0; i < (dtt_size * dtt_size); i++) {
			kernel[indx++] *= s;
		}
 	}

	void clt_symmetrize_kernel( // 
			double []     kernel,      // should be (2*dtt_size-1) * (2*dtt_size-1) +2 size (last 2 are not modified)
			double [][]   sym_kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift  
			final int     dtt_size) // 8
	{
		int in_size = 2*dtt_size-1;
		int dtt_size_m1 = dtt_size - 1;
		int center = dtt_size_m1 * in_size + dtt_size_m1;
		
		for (int i = 0; i < dtt_size; i++){
			for (int j = 0; j < dtt_size; j++){
				int indx0 = center - i * in_size - j;  
				int indx1 = center - i * in_size + j;  
				int indx2 = center + i * in_size - j;  
				int indx3 = center + i * in_size + j;  
				sym_kernels[0][i*dtt_size+j] =                             0.25*( kernel[indx0] + kernel[indx1] + kernel[indx2] + kernel[indx3]);
				if (j > 0)              sym_kernels[1][i*dtt_size+j-1] =   0.25*(-kernel[indx0] + kernel[indx1] - kernel[indx2] + kernel[indx3]);
				if (i > 0)              sym_kernels[2][(i-1)*dtt_size+j] = 0.25*(-kernel[indx0] - kernel[indx1] + kernel[indx2] + kernel[indx3]);
				if ((i > 0) && (j > 0)) sym_kernels[3][i*dtt_size+j] =     0.25*(-kernel[indx0] + kernel[indx1] - kernel[indx2] + kernel[indx3]);
			}
			sym_kernels[1][i*dtt_size + dtt_size_m1] = 0.0;   
			sym_kernels[2][dtt_size_m1*dtt_size + i] = 0.0;   
			sym_kernels[3][i*dtt_size + dtt_size_m1] = 0.0;   
			sym_kernels[3][dtt_size_m1*dtt_size + i] = 0.0;   
		}
 	}

	void clt_dtt3_kernel( // 
			double [][]   kernels, // set of 4 SS, AS, SA, AA kdernels, each dtt_size * dtt_size (may have 5-th with center shift  
			final int     dtt_size, // 8
			DttRad2       dtt)
	{
		if (dtt == null) dtt = new DttRad2(dtt_size);
		for (int quad = 0; quad < 4; quad ++){
			kernels[quad] = dtt.dttt_iiie(kernels[quad], quad, dtt_size);
		}
 	}
	
	
	
	
	public double [][][][] mdctScale(
			final ImageStack                                 imageStack,
			final int                                        subcamera, // 
			final EyesisCorrectionParameters.DCTParameters   dctParameters, //
			final int                                        threadsMax, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			final int                                        debugLevel,
			final boolean                                    updateStatus) // update status info

	{
	  	  if (imageStack==null) return null;
		  final int imgWidth=imageStack.getWidth();
		  final int nChn=imageStack.getSize();
		  double [][][][] dct_data = new double [nChn][][][];
		  float [] fpixels;
		  int i,chn; //tileX,tileY;
		  /* find number of the green channel - should be called "green", if none - use last */
		  // Extract float pixels from inage stack, convert each to double

		  
		  for (chn=0;chn<nChn;chn++) {
			  fpixels= (float[]) imageStack.getPixels(chn+1);
			  double[] dpixels = new double[fpixels.length];
			  for (i = 0; i <fpixels.length;i++) dpixels[i] = fpixels[i];
			  // convert each to DCT tiles
			  dct_data[chn] =lapped_dct_scale(
						dpixels,
						imgWidth,
						dctParameters.dct_size,
						(int) Math.round(dctParameters.dbg_src_size),
						dctParameters.dbg_fold_scale,
						dctParameters.dbg_fold_scale,
						0, //     dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
						dctParameters.dct_window, // final int       window_type,
						chn,
						dctParameters.tileX,
						dctParameters.tileY,
						dctParameters.dbg_mode,
						threadsMax,  // maximal number of threads to launch                         
						debugLevel);
		  }
		return dct_data;
	}
	
	
	
	public double [][][] lapped_dct_scale( // scale image to 8/9 size in each direction
			final double [] dpixels,
			final int       width,
			final int       dct_size,
			final int       src_size,    // source step (== dct_size - no scale, == 9 - shrink, ==7 - expand
			final double    scale_hor,
			final double    scale_vert,
			final int       dct_mode,    // 0: dct/dct, 1: dct/dst, 2: dst/dct, 3: dst/dst
			final int       window_type,
			final int       color,
			final int       debug_tileX,
			final int       debug_tileY,
			final int       debug_mode,
			final int       threadsMax,  // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
		final int height=dpixels.length/width;
		final int n2 = dct_size * 2;

		final int tilesX = (width - n2) / src_size + 1;
		final int tilesY = (height - n2) / src_size + 1;

		final int nTiles=tilesX*tilesY; 
		final double [][][] dct_data = new double[tilesY][tilesX][dct_size*dct_size];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int tileY = 0; tileY < tilesY; tileY++){
			for (int tileX = 0; tileX < tilesX; tileX++){
				for (int i=0; i<dct_data[tileY][tileX].length;i++) dct_data[tileY][tileX][i]= 0.0; // actually not needed, Java initializes arrays
			}
		}
		double [] dc = new double [dct_size*dct_size];
		for (int i = 0; i<dc.length; i++) dc[i] = 1.0;
//		DttRad2 dtt0 = new DttRad2(dct_size);
//		dtt0.set_window(window_type);
//		final double [] dciii = dtt0.dttt_iii  (dc, dct_size);
//		final double [] dciiie = dtt0.dttt_iiie  (dc, 0, dct_size);

		
		if (globalDebugLevel > 0) {
			System.out.println("lapped_dctdc(): width="+width+" height="+height);
		}

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					DttRad2 dtt = new DttRad2(dct_size);
					dtt.set_window(window_type);
					double [] tile_in = new double[4*dct_size * dct_size];
					double [] tile_folded;
					double [] tile_out; // = new double[dct_size * dct_size];
					int tileY,tileX;
					double [][][] fold_k =  dtt.get_fold_2d(
//							int n,
							scale_hor,
							scale_vert
							);
//					double [] tile_out_copy = null;
//					showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						//readDCTKernels() debugLevel = 1 kernels[0].size = 8 kernels[0].img_step = 16 kernels[0].asym_nonzero = 4 nColors = 3 numVert = 123 numHor =  164
						// no aberration correction, just copy data
						for (int i = 0; i < n2;i++){
							System.arraycopy(dpixels, (tileY*width+tileX)*src_size + i*width, tile_in, i*n2, n2);
						}
						tile_folded=dtt.fold_tile(tile_in, dct_size, 0, fold_k); // DCCT
						tile_out=dtt.dttt_iv  (tile_folded, dct_mode, dct_size);

//						if ((tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
//							tile_out_copy = tile_out.clone();
//						}
						
						System.arraycopy(tile_out, 0, dct_data[tileY][tileX], 0, tile_out.length);
					}
				}
			};
		}		      
		startAndJoin(threads);
		return dct_data;
	}
	
	public void dct_scale(
			final double  scale_hor,  // < 1.0 - enlarge in dct domain (shrink in time/space)
			final double  scale_vert, // < 1.0 - enlarge in dct domain (shrink in time/space)
			final boolean normalize, // preserve weighted dct values
			final double [][][] dct_data,
			final int       debug_tileX,
			final int       debug_tileY,
			final int       threadsMax,     // maximal number of threads to launch                         
			final int       globalDebugLevel)
	{
		final int tilesY=dct_data.length;
		final int tilesX=dct_data[0].length;
		final int nTiles=tilesX*tilesY;
		final int dct_size = (int) Math.round(Math.sqrt(dct_data[0][0].length));
		final int dct_len = dct_size*dct_size;
		final double [] norm_sym_weights = new double [dct_size*dct_size];
		for (int i = 0; i < dct_size; i++){
			for (int j = 0; j < dct_size; j++){
				double d = 	Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
				if (i > 0) d*= 2.0;
				if (j > 0) d*= 2.0;
				norm_sym_weights[i*dct_size+j] = d;
			}
		}
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int tileY,tileX;
					double [] dct1 = new double [dct_size*dct_size];
					double [] dct;
					double [][] bidata = new double [2][2];
					int dct_m1 = dct_size - 1;
					int dct_m2 = dct_size - 2;
					for (int nTile = ai.getAndIncrement(); nTile < nTiles; nTile = ai.getAndIncrement()) {
						tileY = nTile/tilesX;
						tileX = nTile - tileY * tilesX;
						dct = dct_data[tileY][tileX];
						double sum_orig=0;
						if (normalize) {
							for (int i = 0; i < dct_len; i++){
								sum_orig += dct[i] *norm_sym_weights[i];
							}
						}
						for (int i = 0; i < dct_size; i++){
							double fi = i * scale_vert;
							int i0 = (int) fi;
							fi -= i0;
							for (int j = 0; j < dct_size; j++){
								double fj = j * scale_hor;
								int j0 = (int) fj;
								fj -= j0;
								int indx = i0*dct_size+j0;
								if ((i0 > dct_m1) || (j0 > dct_m1)){
									bidata[0][0] = 0.0;
									bidata[0][1] = 0.0;
									bidata[1][0] = 0.0;
									bidata[1][1] = 0.0;
								} else {
									bidata[0][0] = dct[indx];
									if (i0 > dct_m2) {
										bidata[1][0] = 0.0;
										bidata[1][1] = 0.0;
										if (j0 > dct_m2) {
											bidata[0][1] = 0.0;
										} else {
											bidata[0][1] = dct[indx + 1];
										}
									} else {
										bidata[1][0] = dct[indx+dct_size];
										if (j0 > dct_m2) {
											bidata[0][1] = 0.0;
											bidata[1][1] = 0.0;
										} else {
											bidata[0][1] = dct[indx + 1];
											bidata[1][1] = dct[indx + dct_size + 1];
										}
									}
									
								}
								// bilinear interpolation 
								dct1[i*dct_size+j] =
										bidata[0][0] * (1.0-fi) * (1.0-fj) +
										bidata[0][1] * (1.0-fi) *      fj  +
										bidata[1][0] *      fi  * (1.0-fj) +
										bidata[1][1] *      fi *       fj;
								if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX)) {
									System.out.println(i+":"+j+" {"+bidata[0][0]+","+bidata[0][1]+","+bidata[1][0]+","+bidata[1][1]+"}, ["+fi+","+fj+"] "+bidata[1][1]);
								}
							}
						}
						if (normalize) {
							double sum=0;
							for (int i = 0; i < dct_len; i++){
								sum += dct1[i] *norm_sym_weights[i];
							}
							if (sum >0.0) {
								double k = sum_orig/sum;
								for (int i = 0; i < dct_len; i++){
									dct1[i] *= k;
								}
							}
						}
//						if ((tileY == debug_tileY) && (tileX == debug_tileX) && (color == 2)) {
						if ((globalDebugLevel > 0) && (tileY == debug_tileY) && (tileX == debug_tileX)) {
							double [][] scaled_tiles = {dct, dct1};
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
							String [] titles = {"orig","scaled"};
							sdfa_instance.showArrays(scaled_tiles,  dct_size, dct_size, true, "scaled_tile", titles);
						}
						
						System.arraycopy(dct1, 0, dct, 0, dct_len); // replace original data
					}
				}
			};
		}		      
		startAndJoin(threads);
	}
	
	
	
	
	
	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	static Thread[] newThreadArray(int maxCPUs) {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	public static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try
		{   
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}
	
}
