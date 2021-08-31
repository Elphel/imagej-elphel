package com.elphel.imagej.tileprocessor;

import java.util.Arrays;
import java.util.Random;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

/**
 **
 ** Clt1d.java - Testing 1D CLT
 **
 ** Copyright (C) 2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Clt1d.java is free software: you can redistribute it and/or modify
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

public class Clt1d {

	public int    transform_size = 8;
	int           transform_len = 64;
	int           nSens =          4;
	DttRad2       dtt;
	ImageDtt      image_dtt;
	Correlation2d corr2d;
	
	// it is just for testing, so use existing 2D for 1D
	public Clt1d(int transform_size) {
		this.transform_size = transform_size;
		this.transform_len = transform_size * transform_size;
		dtt = new DttRad2(transform_size);
		dtt.set_window(1);
		image_dtt = new ImageDtt(
				nSens,
				transform_size,
				null, // FIXME: needs ImageDttParameters (clt_parameters.img_dtt), 
				false,
				false,
				1.0);
		corr2d = new Correlation2d(
				nSens,
				transform_size,             // int transform_size,
				false,                      // boolean monochrome,
				false);                     //   boolean debug)
		
	}
	public void test1d(CLTParameters clt_parameters) {
		  double fat_zero =   clt_parameters.getGpuFatZero(false); // isMonochrome());
		  double sigma_corr = clt_parameters.getGpuCorrSigma(false); // isMonochrome());
		  double corr_scale = 8; // 32.0;
		  double scale =      1.0; // 0.8; // 1.25; // 1; // 1-st channel frequency relative to the second
//		  double freq =       1.3; // periods per 8 pixels
		  double shift_vert =  0.0;
		  double shift_hor1 =  0.0;
		  double shift_hor2 =  0.0;
		  long    seed =       1; // 0;
		  double  sigma =      0.5; // 1.5;
		  int     oversample = 16;
		  double  offset0 =    -0.5;
		  double  offset1 =    0.0;
// fill in signal 1
		  int n = transform_size;
		  int n2 = 2 * n;
		  /*
		  double [] signal1 = new double [n2];
		  double [] signal2 = new double [n2];
		  for (int i = 0; i < n2; i++) {
			  signal1[i] = Math.sin((i - n) * scale * freq * 2 * Math.PI / n);
			  signal2[i] = Math.sin((i - n) *         freq * 2 * Math.PI / n);
		  }

		  double [][] signals_2d = {setTile2d(signal1), setTile2d(signal2)};
		  */
		  double [][] signals_2d = getRandomStereo(
				  seed,       // long    seed,
				  sigma,      // double  sigma,
				  oversample, // int     oversample,
				  offset0,    // double  offset0, // -0.5...+0.5
				  offset1,    // double  offset1, // -0.5...+0.5
				  scale,      // double  scale0,   // horizontal "reciprocal zoom" for view 0
				  true);      // boolean debug);		
		  
		  double [][][] signals_td = new double [2][][]; 
		  signals_td[0] = foldConvertShiftTile(
				  signals_2d[0],  // double [] tile_in,
				  shift_hor1,  // double shift_hor,
				  shift_vert); // double shift_vert
		  signals_td[1] = foldConvertShiftTile(
				  signals_2d[1],  // double [] tile_in,
				  shift_hor2,  // double shift_hor,
				  shift_vert); // double shift_vert
		  double [] corr= correlate(
				  signals_td[0],  // double [][] tile1,
				  signals_td[1],  // double [][] tile2,
				  2.0 , // sigma_corr,  // 	double      corr_sigma,
				  5); //fat_zero);   // double      fat_zero_abs)
		  double [][] dbg_img = new double[3][];
		  dbg_img[0] = signals_2d[0].clone();
		  dbg_img[1] = signals_2d[1].clone();
		  dbg_img[2] = new double [n2*n2];
		  Arrays.fill(dbg_img[2], Double.NaN);
		  for (int y = 0; y < (n2 - 1); y++) {
			  for (int x = 0; x < (n2 - 1); x++) {
				  dbg_img[2][y * n2 + x] = corr_scale * corr [ y * (n2-1) + x];
			  }
		  }
		  double [] corr1d= new double[n2-1];
		  System.arraycopy(corr, (n-1)*(n2-1), corr1d, 0, n2-1);
		  double ym = corr_scale * corr1d[n-2];
		  double y0 = corr_scale * corr1d[n-1];
		  double yp = corr_scale * corr1d[n-0];
		  double w = Math.cos(Math.PI/n);
		  ym /= w;
		  yp /= w;
		  //		  double dx = (yp-ym)/(2*(ym+yp+2*y0));
		  double dx = (yp-ym)/(2*(2*y0 - ym-yp));
		  System.out.println(String.format("ym=%6f y0=%6f yp=%6f dx=%6f",ym, y0,yp, dx));
		  //=(K41-K43)/(2*(K41+K43+2*K42))
		  //  double dx = (yp-ym)/(2*(ym+yp+2*y0));

		  String [] dbg_titles = {"first", "second", "corr"};
		  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
				  dbg_img,
				  n2,
				  n2,
				  true,
				  "corr_2d",
				  dbg_titles);
	}
	
	public double [][] getRandomStereo(
			long    seed,
			double  sigma,
			int     oversample,
			double  offset0, // -0.5...+0.5
			double  offset1, // -0.5...+0.5
			double  scale0,   // horizontal "reciprocal zoom" for view 0
			boolean debug)
	{
		int n = transform_size;
		int n2 = n*2;
		int n2n2 = n2*n2;
		double hrange1 = n+1;
		if (scale0 > 1) hrange1 *= scale0;
		int ihrange1 = (int) Math.ceil(hrange1);
		int swidth = 2* oversample * ihrange1;
		int sheight = n2 * oversample;
		int slength = swidth * sheight;
		int margin = (int) Math.ceil(sigma * oversample);
		int swidth2 = swidth + 2*margin;
		int sheight2 = sheight + 2*margin;
		int slength2 = swidth2 * sheight2;
		double [] sdata2 =new double[slength2];
		Random random = new Random(seed);
		for (int i = 0; i < slength; i++) {
			sdata2[i] = random.nextGaussian();
		}
		(new DoubleGaussianBlur()).blurDouble(sdata2,  swidth2, sheight2, sigma * oversample, sigma * oversample, 0.01);
		
		double [] sdata =new double[slength];
		for (int row = 0; row < sheight; row++) {
			System.arraycopy(sdata2, (row+margin)*swidth2 + margin, sdata, row * swidth, swidth);
		}
		double s2 = 0.0;
		for (int row = 0; row < sheight; row++) {
			for (int col = 0; col < swidth; col++) {
				double d = sdata[row * swidth + col];
				s2 += d*d;
			}
		}
		double std = Math.sqrt(s2/(swidth * sheight));
		for (int row = 0; row < sheight; row++) {
			for (int col = 0; col < swidth; col++) {
				double d = sdata[row * swidth + col];
				sdata[row * swidth + col] = 0.5*(d/(2*std) + 1.0);
			}
		}
		

		
		
		
		if (debug) {
			  (new ShowDoubleFloatArrays()).showArrays( // out of boundary 15
					  new double [][] {sdata},
					  swidth,
					  sheight,
					  true,
					  "random-source-data");
		}
		double [][] pair = new double [2][n2n2];
		for (int chn = 0; chn < 2; chn++) {
			double offset = (chn > 0) ? offset1 : offset0;
			double scale =  (chn > 0) ? 1.0 : scale0;
			for (int row = 0; row < n2; row++) {
				int srow = row * oversample;
				for (int col = 0; col < n2; col++) {
					double dcol = ((col - n) * scale - offset) * oversample + (swidth/2); // offset sign?
					// linear interpolate
					int    icol = (int) Math.floor(dcol);
					double fcol = dcol-icol;
					int sindx = srow * swidth + icol;
					pair[chn][row * n2 + col] = (1.0 - fcol) * sdata[sindx] + fcol * sdata[sindx + 1]; 
				}
			}
		}
		return pair;
	}
	
	
	double [] setTile2d (double [] x){ // x should be 16-long
		int l = x.length;
		double [] tile = new double [l * l];
		for (int i = 0; i < l; i++) {
			if (i == l/2) {
				System.arraycopy(x, 0, tile, i*l, l);
			}
		}
		return tile;
	}

	
	
	
	double [][] foldConvertShiftTile(
			double [] tile_in,
			double shift_hor,
			double shift_vert
			) {
		double [][][] fold_coeff = dtt.get_shifted_fold_2d ( // get_shifted_fold_2d(
				transform_size,
				shift_hor,
				shift_vert,
				0); // debug level
		double [][] tile = new double [4][];
		for (int dct_mode = 0; dct_mode < 4; dct_mode++) {
			tile[dct_mode] = dtt.fold_tile (tile_in,        transform_size, dct_mode, fold_coeff);
			tile[dct_mode] = dtt.dttt_iv   (tile[dct_mode], dct_mode,       transform_size);
		}
		// Apply shift
        image_dtt.fract_shift( // fractional shift in transform domain. Currently uses sin/cos - change to tables with 2? rotations
        		tile,          // double  [][]  clt_tile,
        		shift_hor,     // double        shiftX,
        		shift_vert,    // double        shiftY,
                false);        // debug);
		return tile;
	}
	//		final double [] filter =     image_dtt.doubleGetCltLpfFd(corr_sigma);
	double [] correlate(
			double [][] tile1,
			double [][] tile2,
			double      corr_sigma,
			double      fat_zero_abs) {
		double [] filter =     image_dtt.doubleGetCltLpfFd(corr_sigma);
		/*
		double[][] corr_td = corr2d.correlateSingleColorFD( // USED in lwir
				tile1,      // double [][] clt_data1,
				tile2,      // double [][] clt_data2,
	    		null,       // double [][] tcorr, // null or initialized to [4][transform_len]
	    		fat_zero);  // double      fat_zero);
		*/
		double[][] corr_td = corr2d.correlateSingleColorFD( // USED in lwir
				tile1,      // double [][] clt_data1,
				tile2,      // double [][] clt_data2,
	    		null);       // double [][] tcorr, // null or initialized to [4][transform_len]
		// normalize as in GPU
		final double fat_zero_abs2 = fat_zero_abs * fat_zero_abs; 
        for (int i = 0; i < transform_len; i++) {
            double s2 = fat_zero_abs2;
            for (int q = 0; q < 4; q++) {
                double d = corr_td[q][i]; 
                s2 += d*d;
            }
            double k = 1.0/Math.sqrt(s2);
            for (int q = 0; q < 4; q++) {
            	corr_td[q][i] *= k;
            }
        }
        // apply lpf:
    	if (filter != null) {
    		for (int n = 0; n < 4; n++) {
    			for (int i = 0; i < transform_len; i++) {
    				corr_td[n][i] *= filter[i];
    			}
    		}
    	}
    	for (int quadrant = 0; quadrant < 4; quadrant++){
    		int mode = ((quadrant << 1) & 2) | ((quadrant >> 1) & 1); // transpose
    		corr_td[quadrant] = dtt.dttt_iie(corr_td[quadrant], mode, transform_size);
    	}
		// convert from 4 quadrants to 15x15 centered tiles (only composite)
    	double [] corr_pd =  dtt.corr_unfold_tile(corr_td,	transform_size);
    	return corr_pd;
	}
	
}
