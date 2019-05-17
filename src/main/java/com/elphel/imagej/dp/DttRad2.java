package com.elphel.imagej.dp;
/**
 **
 ** DttRad2 - Calculate DCT types II and IV for n=2^t and n*n 2-d
 ** also DST-IV and some other related transforms
 **
 ** Uses algorithm described in
 ** Plonka, Gerlind, and Manfred Tasche. "Fast and numerically stable algorithms for discrete cosine transforms."
 ** Linear algebra and its applications 394 (2005): 309-345.
 **
 ** Copyright (C) 2016 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  DttRad2.java is free software: you can redistribute it and/or modify
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

public class DttRad2 {
	int N = 0;
	double [][][] CII=    null;
	double [][][] CIIe=   null; // alternative matrix with all coefficients the same (non-orthogonal, but matching DFT)
	double [][][] CIIIe=  null; // alternative matrix with k0=1/2, k(n-1) = 1/2 (non-orthogonal, but matching DFT)
	double [][][] SIIe=   null; // alternative matrix with all coefficients the same (non-orthogonal, but matching DFT)
	double [][][] SIIIe=  null; // alternative matrix with k0=1/2, k(n-1) = 1/2 (non-orthogonal, but matching DFT)

	double [][][] CIV=  null;
	double [][][] SIV=  null;
	double [][] CN1=null;
	double [][] SN1=null;
	double      COSPI_1_8_SQRT2 = Math.cos(Math.PI/8)*Math.sqrt(2.0);
	double      COSPI_3_8_SQRT2 = Math.cos(3*Math.PI/8)*Math.sqrt(2.0);
	double sqrt2 = Math.sqrt(2.0);
	double sqrt1_2 = 1/sqrt2;
	double [] hwindow = null; // half window
	int    [][] fold_index = null; // index of the source item in 2nx2n array input to mdct_2d.
	                               // First index (0..n^2-1) index in the folded array (dct-iV input)
	                               // Second index(0..3) - item to add (2 vertical, 2 - horizontal)
	double [][][] fold_k = null; // First index - mode: 0 - CC 1: SC, 2: CS, 3: SS. Other indices matching fold_index items. Each is a product of 2 window coefficients and sign
	int    [] unfold_index = null;  // index  for each element of idct(2nx2n)
	double [][] unfold_k = null; // First index - mode: 0 - CC 1: SC, 2: CS, 3: SS. Other indices matching unfold_index items. Each is a product of 2 window coefficients and sign

	public int [][] getFoldIndex(){
		return fold_index;
	}
	public int [] getUnfoldIndex(){
		return unfold_index;
	}

	public double [][][] getFoldK(){
		return fold_k;
	}

	public DttRad2 (int maxN){ // n - maximal
		setup_arrays(maxN); // always setup arrays for fast calculations
	}

	public double [] dct_ii(double[] x){
		if (x.length > N){
			N = x.length;
		}
		double [] y=  _dctii_recurs(x);
		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) y[i] *= scale;
		return y;
	}

	public double [] dct_iv(double[] x){
		double [] y=  _dctiv_recurs(x);
		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) y[i] *= scale;
		return y;
	}

	public double [] dst_iv(double[] x){
		double [] xr= new double[x.length];
		int j= x.length-1;
		for (int i=0; i < x.length;i++) xr[i] = x[j--];
		double [] y=  _dctiv_recurs(xr);
		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) {
			y[i] *= scale;
			scale = -scale;
		}
		return y;
	}


	// For index in dct-iv input (0..n-1) get 2 variants of index in mdct input array (0..2*n-1)
	// second index : 0 - index in X array 2*n long
	//                1 - window index (0..n-1), [0] - minimal, [n-1] - max
	//                2 - sign of the C term  (-~c - d, +a -~b)
	//                3 - sign of the S term  (+~c - d, +a +~b)
	// Added for shift (+/- 0.5), using window derivative for correction:
	//                4 - window index (0..n-1), [0] - minimal, [n-1] - max (cos: same table, just different order
	//                5 - sign of derivative (to multiply by shift)
	private int [][] get_fold_indices(int x, int n){
		int n1 = n>>1;
		int [][] ind = new int[2][6];
		if (x <n1) {
			ind[0][0] = n + n1 - x - 1; // C: -cR, S: +cR
			ind[0][1] = n1     + x;
			ind[0][2] = -1;
			ind[0][3] =  1;
			ind[0][4] = n1     - x -1;
			ind[0][5] =  -1; // c - window derivative over shift is negative

			ind[1][0] = n + n1 + x;     // C: -d,  S: -d
			ind[1][1] = n1     - x - 1;
			ind[1][2] = -1;
			ind[1][3] = -1;

			ind[1][4] = n1     + x;
			ind[1][5] =  -1; // d - window derivative over shift is negative
		} else {
			x-=n1;
			ind[0][0] =          x;     // C: +a, S: +a
			ind[0][1] =          x;
			ind[0][2] =  1;
			ind[0][3] =  1;
			ind[0][4] =  n     - x - 1;
			ind[0][5] =  1;   // a - window derivative over shift is positive

			ind[1][0] = n      - x - 1; // C: -bR, S: +bR
			ind[1][1] = n      - x - 1;
			ind[1][2] = -1;
			ind[1][3] =  1;
			ind[1][4] =          x;
			ind[1][5] =  1;   // b - window derivative over shift is positive
		}
		return ind;
	}
	// is called when window is set
	private void set_fold_2d(int n){ // n - DCT and window size
		if ((fold_index != null) && (fold_index.length == n*n)) return;
		fold_index = new int[n*n][4];
		fold_k =     new double[4][n*n][4];
		int []    vert_ind = new int[2];
		double [][] vert_k = new double[2][2];
		int    [] hor_ind = new int[2];
		double [][] hor_k = new double[2][2];
		int [][] fi;
		int n2 = 2*n;
		for (int i = 0; i < n; i++ ){
			fi = get_fold_indices(i,n);
			vert_ind[0] = fi[0][0];
			vert_ind[1] = fi[1][0];
			vert_k[0][0] =   fi[0][2] * hwindow[fi[0][1]]; // use cosine sign
			vert_k[0][1] =   fi[1][2] * hwindow[fi[1][1]]; // use cosine sign
			vert_k[1][0] =   fi[0][3] * hwindow[fi[0][1]]; // use sine sign
			vert_k[1][1] =   fi[1][3] * hwindow[fi[1][1]]; // use sine sign
			for (int j = 0; j < n; j++ ){
				fi = get_fold_indices(j,n);
				hor_ind[0] = fi[0][0];
				hor_ind[1] = fi[1][0];
				hor_k[0][0] =   fi[0][2] * hwindow[fi[0][1]]; // use cosine sign
				hor_k[0][1] =   fi[1][2] * hwindow[fi[1][1]]; // use cosine sign
				hor_k[1][0] =   fi[0][3] * hwindow[fi[0][1]]; // use sine sign
				hor_k[1][1] =   fi[1][3] * hwindow[fi[1][1]]; // use sine sign
				int indx = n*i + j;
				for (int k = 0; k<4;k++) {
					fold_index[indx][k] = n2 * vert_ind[(k>>1) & 1] + hor_ind[k & 1];
				}
				for (int mode = 0; mode<4; mode++){
					for (int k = 0; k<4;k++) {
						fold_k[mode][indx][k] =     vert_k[(mode>>1) &1][(k>>1) & 1] * hor_k[mode &1][k & 1];
					}
				}
			}
		}
		if (n < 8) {
			for (int i = 0; i < n; i++ ){
				fi = get_fold_indices(i,n);
				System.out.println(i+"->"+String.format("?[%2d %2d %2d %2d] [%2d %2d %2d %2d] %f %f",
						fi[0][0],fi[0][1],fi[0][2],fi[0][3],
						fi[1][0],fi[1][1],fi[1][2],fi[1][3], hwindow[fi[0][1]], hwindow[fi[1][1]]));
			}
			for (int i = 0; i < n*n; i++){
				System.out.println(String.format("%3x:   %6x   %6x   %6x   %6x",i,fold_index[i][0],fold_index[i][1],fold_index[i][2],fold_index[i][3]));
//				System.out.println(String.format("   : %8.5f %8.5f %8.5f %8.5f", fold_k[0][i][0],  fold_k[0][i][1], fold_k[0][i][2], fold_k[0][i][3]));
//				System.out.println(String.format("   : %8.5f %8.5f %8.5f %8.5f", fold_k[1][i][0],  fold_k[1][i][1], fold_k[1][i][2], fold_k[1][i][3]));
//				System.out.println(String.format("   : %8.5f %8.5f %8.5f %8.5f", fold_k[2][i][0],  fold_k[2][i][1], fold_k[2][i][2], fold_k[2][i][3]));
//				System.out.println(String.format("   : %8.5f %8.5f %8.5f %8.5f", fold_k[3][i][0],  fold_k[3][i][1], fold_k[3][i][2], fold_k[3][i][3]));
				System.out.println(String.format("   :   %2d   %2d   %2d   %2d", (fold_k[0][i][0]<0)?-1:1,  (fold_k[0][i][1]<0)?-1:1, (fold_k[0][i][2]<0)?-1:1, (fold_k[0][i][3]<0)?-1:1));
				System.out.println(String.format("   :   %2d   %2d   %2d   %2d", (fold_k[1][i][0]<0)?-1:1,  (fold_k[1][i][1]<0)?-1:1, (fold_k[1][i][2]<0)?-1:1, (fold_k[1][i][3]<0)?-1:1));
				System.out.println(String.format("   :   %2d   %2d   %2d   %2d", (fold_k[2][i][0]<0)?-1:1,  (fold_k[2][i][1]<0)?-1:1, (fold_k[2][i][2]<0)?-1:1, (fold_k[2][i][3]<0)?-1:1));
				System.out.println(String.format("   :   %2d   %2d   %2d   %2d", (fold_k[3][i][0]<0)?-1:1,  (fold_k[3][i][1]<0)?-1:1, (fold_k[3][i][2]<0)?-1:1, (fold_k[3][i][3]<0)?-1:1));
			}

		}
	}

	public double [][][] get_fold_2d(
			double scale_hor,
			double scale_vert)
	{
		return get_fold_2d( this.N, scale_hor, scale_vert);
	}
	public double [][][] get_fold_2d(
			int n,
			double scale_hor,
			double scale_vert)
	{ // n - DCT and window size

			double [] hwindow_h = new double[n];
			double [] hwindow_v = new double[n];
			double f = Math.PI/(2.0*n);
			for (int i = 0; i < n; i++ ) {
				double ah = f*scale_hor * (n-i-0.5);
				double av = f*scale_vert * (n-i-0.5);
				hwindow_h[i] = (ah > (Math.PI/2))? 0.0: Math.cos(ah);
				hwindow_v[i] = (av > (Math.PI/2))? 0.0: Math.cos(av);
			}
		double [][][] fold_sk = new double[4][n*n][4]; // was [2][n*n][4]
		int []    vert_ind =    new int[2];
		double [][] vert_k =    new double[2][2];
		int    [] hor_ind =     new int[2];
		double [][] hor_k =     new double[2][2];
		int [][] fi;
		for (int i = 0; i < n; i++ ){
			fi = get_fold_indices(i,n);
			vert_ind[0] = fi[0][0];
			vert_ind[1] = fi[1][0];
			vert_k[0][0] =   fi[0][2] * hwindow_v[fi[0][1]]; // use cosine sign
			vert_k[0][1] =   fi[1][2] * hwindow_v[fi[1][1]]; // use cosine sign
			vert_k[1][0] =   fi[0][3] * hwindow_v[fi[0][1]]; // use sine sign
			vert_k[1][1] =   fi[1][3] * hwindow_v[fi[1][1]]; // use sine sign

			for (int j = 0; j < n; j++ ){
				fi = get_fold_indices(j,n);
				hor_ind[0] = fi[0][0];
				hor_ind[1] = fi[1][0];
				hor_k[0][0] =   fi[0][2] * hwindow_h[fi[0][1]]; // use cosine sign
				hor_k[0][1] =   fi[1][2] * hwindow_h[fi[1][1]]; // use cosine sign
				hor_k[1][0] =   fi[0][3] * hwindow_h[fi[0][1]]; // use sine sign
				hor_k[1][1] =   fi[1][3] * hwindow_h[fi[1][1]]; // use sine sign

				int indx = n*i + j;

				for (int mode = 0; mode<4; mode++){
					for (int k = 0; k<4;k++) {
						fold_sk[mode][indx][k] =     vert_k[(mode>>1) &1][(k>>1) & 1] * hor_k[mode &1][k & 1];
					}
				}
			}
		}
		return fold_sk;
	}

	// Generate (slightly - up to +/- 0.5) shifted window for standard sin window (derivative uses same table)
	public double [][][] get_shifted_fold_2d(
			int n,
			double shift_hor,
			double shift_vert,
			int debugLevel)
	{ // n - DCT and window size
		int n2 = 2* n;
		double [][][] fold_sk = new double[4][n*n][4];
		int []    vert_ind =    new int[2];
		double [][] vert_k =    new double[2][2];
		int    [] hor_ind =     new int[2];
		double [][] hor_k =     new double[2][2];
		double ahc = Math.cos(Math.PI/n2*shift_hor);
		double ahs = Math.sin(Math.PI/n2*shift_hor);
		double avc = Math.cos(Math.PI/n2*shift_vert);
		double avs = Math.sin(Math.PI/n2*shift_vert);

		int [][] fi;
		for (int i = 0; i < n; i++ ){
			fi = get_fold_indices(i,n);
			vert_ind[0] = fi[0][0];
			vert_ind[1] = fi[1][0];
			double vw0 = avc*hwindow[fi[0][1]] + avs*hwindow[fi[0][4]]*fi[0][5];
			double vw1 = avc*hwindow[fi[1][1]] + avs*hwindow[fi[1][4]]*fi[1][5];
			vert_k[0][0] =   fi[0][2] * vw0; // use cosine sign
			vert_k[0][1] =   fi[1][2] * vw1; // use cosine sign
			vert_k[1][0] =   fi[0][3] * vw0; // use sine sign
			vert_k[1][1] =   fi[1][3] * vw1; // use sine sign

			for (int j = 0; j < n; j++ ){
				fi = get_fold_indices(j,n);
				hor_ind[0] = fi[0][0];
				hor_ind[1] = fi[1][0];
				double hw0 = ahc*hwindow[fi[0][1]] + ahs*hwindow[fi[0][4]]*fi[0][5];
				double hw1 = ahc*hwindow[fi[1][1]] + ahs*hwindow[fi[1][4]]*fi[1][5];

				hor_k[0][0] =   fi[0][2] * hw0; // use cosine sign
				hor_k[0][1] =   fi[1][2] * hw1; // use cosine sign
				hor_k[1][0] =   fi[0][3] * hw0; // use sine sign
				hor_k[1][1] =   fi[1][3] * hw1; // use sine sign

				int indx = n*i + j;

				for (int mode = 0; mode<4; mode++){
					for (int k = 0; k<4;k++) {
						fold_sk[mode][indx][k] =     vert_k[(mode>>1) &1][(k>>1) & 1] * hor_k[mode &1][k & 1];
					}
				}
			}
		}
		return fold_sk;
	}

	// Generate (slightly - up to +/- 0.5) shifted window for standard sin window (derivative uses same table)
	// only for mode=1 (sin window)
	public double [][][] get_shifted_fold_2d_direct(
			int n,
			double shift_hor,
			double shift_vert,
			int debugLevel) // fpga scale (1 << 17 -1)
	{ // n - DCT and window size

		int n2 = 2 * n;
		double [][][] fold_sk = new double[4][n*n][4];
		int []    vert_ind =    new int[2];
		double [][] vert_k =    new double[2][2];
		int    [] hor_ind =     new int[2];
		double [][] hor_k =     new double[2][2];
//		double ahc = Math.cos(Math.PI/n2*shift_hor);
//		double ahs = Math.sin(Math.PI/n2*shift_hor);
//		double avc = Math.cos(Math.PI/n2*shift_vert);
//		double avs = Math.sin(Math.PI/n2*shift_vert);
		double [] wnd_hor = new double[n2];
		double [] wnd_vert = new double[n2];
		double f = Math.PI/n2;
		for (int i = 0; i < n2; i++ ){
			wnd_hor[i] =  Math.sin(f * (i+ 0.5 +shift_hor));
			wnd_vert[i] = Math.sin(f * (i+ 0.5 +shift_vert));
		}
		if (debugLevel > 0){
			System.out.println("Window (index,hor,vert) ");
			for (int i = 0; i <16; i++) System.out.print(String.format("%5x ", i)); System.out.println();
			for (int i = 0; i <16; i++) System.out.print(String.format("%5x ", (int) Math.round(debugLevel * wnd_hor[i]))); System.out.println();
			for (int i = 0; i <16; i++) System.out.print(String.format("%5x ", (int) Math.round(debugLevel * wnd_vert[i]))); System.out.println();
			System.out.println();
		}
		int [][] fi;
		for (int i = 0; i < n; i++ ){
			fi = get_fold_indices(i,n);
			vert_ind[0] = fi[0][0];
			vert_ind[1] = fi[1][0];
//			double vw0 = avc*hwindow[fi[0][1]] + avs*hwindow[fi[0][4]]*fi[0][5];
//			double vw1 = avc*hwindow[fi[1][1]] + avs*hwindow[fi[1][4]]*fi[1][5];
			double vw0 = wnd_vert[vert_ind[0]];
			double vw1 = wnd_vert[vert_ind[1]];
			vert_k[0][0] =   fi[0][2] * vw0; // use cosine sign
			vert_k[0][1] =   fi[1][2] * vw1; // use cosine sign
			vert_k[1][0] =   fi[0][3] * vw0; // use sine sign
			vert_k[1][1] =   fi[1][3] * vw1; // use sine sign

			for (int j = 0; j < n; j++ ){
				fi = get_fold_indices(j,n);
				hor_ind[0] = fi[0][0];
				hor_ind[1] = fi[1][0];
//				double hw0 = ahc*hwindow[fi[0][1]] + ahs*hwindow[fi[0][4]]*fi[0][5];
//				double hw1 = ahc*hwindow[fi[1][1]] + ahs*hwindow[fi[1][4]]*fi[1][5];
				double hw0 = wnd_hor[hor_ind[0]];
				double hw1 = wnd_hor[hor_ind[1]];

				hor_k[0][0] =   fi[0][2] * hw0; // use cosine sign
				hor_k[0][1] =   fi[1][2] * hw1; // use cosine sign
				hor_k[1][0] =   fi[0][3] * hw0; // use sine sign
				hor_k[1][1] =   fi[1][3] * hw1; // use sine sign

				int indx = n*i + j;

				for (int mode = 0; mode<4; mode++){
					for (int k = 0; k<4;k++) {
						fold_sk[mode][indx][k] =     vert_k[(mode>>1) &1][(k>>1) & 1] * hor_k[mode &1][k & 1];
					}
				}
			}
		}
		return fold_sk;
	}


	// return index and two signs (c,s) for 1-d imdct. x is index (0..2*n-1) of the imdct array, value is sign * (idct_index+1),
	// where idct_index (0..n-1) is index in the dct-iv array
	private int [] get_unfold_index_signs(int x, int n){
		int n1 = n>>1;
		int segm = x / n1;
		x = x % n1;
		int [] is2  = new int[3];
		switch (segm){
		case 0: is2[0] =  x + n1;     is2[1]=  1; is2[2] =  1; return is2; // return 1+ (x + n1);
		case 1: is2[0] =  n -  x - 1; is2[1]= -1; is2[2] =  1; return is2; // return -(n - x);
		case 2: is2[0] =  n1 - x - 1; is2[1]= -1; is2[2] =  1; return is2; // return -(n1 - x);
		case 3: is2[0] =  x;          is2[1]= -1; is2[2] = -1; return is2; // return -(1 + x);
		}
		return null; //should never happen
	}
	private void set_unfold_2d(int n){ // n - DCT size
		if ((unfold_index != null) && (unfold_index.length == 4*n*n)) return;
		unfold_index = new int[4*n*n];
		unfold_k = new double[4][4*n*n];
		int n2 = 2*n;
		for (int i = 0; i < 2*n; i++ ){
			int [] is2_vert = get_unfold_index_signs(i,n);
			double [] k_vert = {is2_vert[1]*hwindow[(i < n)?i:n2 -i -1], is2_vert[2]*hwindow[(i < n)?i:n2 -i -1]};
			int index_vert = is2_vert[0] * n;

			for (int j = 0; j < 2*n; j++ ){
				int [] is2_hor = get_unfold_index_signs(j,n);
				int index_hor = is2_hor[0];
				double [] k_hor = {is2_hor[1] * hwindow[(j < n)?j:n2 -j -1], is2_hor[2] * hwindow[(j < n)?j:n2 -j -1]};
				unfold_index[n2*i+j]=(index_vert+index_hor);
				for (int mode = 0; mode < 4; mode++){
					unfold_k[mode][n2*i+j] = k_vert[(mode>>1) &1]*k_hor[mode &1];
				}

				if (n < 8) System.out.print(String.format("%4d", unfold_index[n2*i+j]));

			}
			if (n < 8) System.out.println();
		}
		if (n < 8) {
			for (int i = 0; i < 2*n; i++ ){
				System.out.println(i+"=>"+get_unfold_index_signs(i,n)[0]+", "+get_unfold_index_signs(i,n)[1]+", "+get_unfold_index_signs(i,n)[2]);
			}
		}
	}

	/**
	 * Unfolds 2d correlation tile in pixel domain
	 * @param qdata 4-quadrant result of 2cd DCT2
	 * @param transform_size DCT transformation size
	 * @return packed array [(2*transform_size-1) * (2*transform_size-1)]
	 */

	public  double [] corr_unfold_tile(
			double [][]  qdata, // [4][transform_size*transform_size] data after DCT2 (pixel domain)
			int          transform_size
			)
	{
		int corr_pixsize = transform_size * 2 - 1;
		double corr_pixscale = 0.25;
		double [] rslt = new double [corr_pixsize*corr_pixsize];
		rslt[corr_pixsize*transform_size - transform_size] = corr_pixscale * qdata[0][0]; // center
		for (int j = 1; j < transform_size; j++) { //  for i == 0
			rslt[corr_pixsize*transform_size - transform_size + j] = corr_pixscale * (qdata[0][j] + qdata[1][j-1]);
			rslt[corr_pixsize*transform_size - transform_size - j] = corr_pixscale * (qdata[0][j] - qdata[1][j-1]);
		}
		for (int i = 1; i < transform_size; i++) {
			rslt[corr_pixsize*(transform_size + i) - transform_size] =
					corr_pixscale * (qdata[0][i*transform_size] + qdata[2][(i-1)*transform_size]);
			rslt[corr_pixsize*(transform_size - i) - transform_size] =
					corr_pixscale * (qdata[0][i*transform_size] - qdata[2][(i-1)*transform_size]);
			for (int j = 1; j < transform_size; j++) {
				rslt[corr_pixsize*(transform_size + i) - transform_size + j] =
						corr_pixscale * (qdata[0][i*    transform_size + j] +
								 qdata[1][i*    transform_size + j - 1] +
								 qdata[2][(i-1)*transform_size + j] +
								 qdata[3][(i-1)*transform_size + j - 1]);

				rslt[corr_pixsize*(transform_size + i) - transform_size - j] =
						corr_pixscale * ( qdata[0][i*    transform_size + j] +
								 -qdata[1][i*    transform_size + j - 1] +
								  qdata[2][(i-1)*transform_size + j] +
								 -qdata[3][(i-1)*transform_size + j - 1]);
				rslt[corr_pixsize*(transform_size - i) - transform_size + j] =
						corr_pixscale * (qdata[0][i*    transform_size + j] +
								 qdata[1][i*    transform_size + j - 1] +
								 -qdata[2][(i-1)*transform_size + j] +
								 -qdata[3][(i-1)*transform_size + j - 1]);
				rslt[corr_pixsize*(transform_size - i) - transform_size - j] =
						corr_pixscale * (qdata[0][i*    transform_size + j] +
								 -qdata[1][i*    transform_size + j - 1] +
								 -qdata[2][(i-1)*transform_size + j] +
								 qdata[3][(i-1)*transform_size + j - 1]);
			}
		}


		return rslt;

	}






	public double [] dttt_iv(double [] x){
		return dttt_iv(x, 0, 1 << (ilog2(x.length)/2));
	}

	public double [] dttt_iv(double [] x, int mode){
		return dttt_iv(x, mode, 1 << (ilog2(x.length)/2));
	}
	public double [] dttt_iv(double [] x, int mode, int n){ // mode 0 - dct,dct 1:dst,dct, 2: dct, dst, 3: dst,dst
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = ((mode & 1)!=0)? dst_iv(line):dct_iv(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = ((mode & 2)!=0)? dst_iv(line):dct_iv(line);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}

	public double [] dttt_iv(double [] x, int mode, int n, double scale, int mask){ // mode 0 - dct,dct 1:dst,dct, 2: dct, dst, 3: dst,dst
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		System.out.println("dttt_iv, mode="+mode);
		System.out.println("horizontal pass "+(((mode & 1)!=0)? "dst_iv":"dct_iv"));
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = ((mode & 1)!=0)? dst_iv(line):dct_iv(line);
			System.out.print(String.format("%02x: ", i));
			for (int j=0; j < n;j++) {
				int di = ((int) Math.round(x[n*i+j]*scale)) & mask;
				System.out.print(String.format("%07x ", di));
			}
			System.out.print("  ");
			for (int j=0; j < n;j++) {
				int di = ((int) Math.round(line[j] * scale)) & mask;
				System.out.print(String.format("%07x ", di));
			}
			System.out.println();
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		System.out.println("vertical pass "+(((mode & 2)!=0)? "dst_iv":"dct_iv")+" (after transpose)");
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = ((mode & 2)!=0)? dst_iv(line):dct_iv(line);
			System.out.print(String.format("%02x: ", i));
			for (int j=0; j < n;j++) {
				int di = ((int) Math.round(y[n*i+j]*scale*0.5)) & mask;
				System.out.print(String.format("%07x ", di));
			}
			System.out.print("  ");
			for (int j=0; j < n;j++) {
				int di = ((int) Math.round(line[j] * scale)) & mask;
				System.out.print(String.format("%07x ", di));
			}
			System.out.println();
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}



	public double [] dttt_ii(double [] x){
		return dttt_ii(x, 1 << (ilog2(x.length)/2));
	}

	public double [] dttt_ii(double [] x, int n){
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = dctii_direct(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = dctii_direct(line);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}

	public double [] dttt_iie(double [] x){
		return dttt_iie(x,0, 1 << (ilog2(x.length)/2));
	}
	public double [] dttt_iie(double [] x, int mode){
		return dttt_iie(x, mode, 1 << (ilog2(x.length)/2));
	}

	public double [] dttt_iie(double [] x, int mode, int n){

		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
//			line = dctiie_direct(line);
			line = ((mode & 1)!=0)? dstiie_direct(line):dctiie_direct(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
//			line = dctiie_direct(line);
			line = ((mode & 2)!=0)? dstiie_direct(line):dctiie_direct(line);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}




	public double [] dttt_iii(double [] x){
		return dttt_iii(x, 1 << (ilog2(x.length)/2));
	}

	public double [] dttt_iii(double [] x, int n){
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = dctiii_direct(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = dctiii_direct(line);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}

	public double [] dttt_iiie(double [] x){
		return dttt_iiie(x,0, 1 << (ilog2(x.length)/2));
	}
	public double [] dttt_iiie(double [] x, int mode){
		return dttt_iiie(x, mode, 1 << (ilog2(x.length)/2));
	}

	public double [] dttt_iiie(double [] x, int mode, int n){
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
//			line = dctiiie_direct(line);
			line = ((mode & 1)!=0)? dstiiie_direct(line):dctiiie_direct(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
//			line = dctiiie_direct(line);
			line = ((mode & 2)!=0)? dstiiie_direct(line):dctiiie_direct(line);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}

	public void set_window(){
		set_window(0);
	}
	public void set_window(int mode){
		set_window(mode, N);
	}
	public void set_window(int mode, int len){ // using mode==1
		// for N=8: sin (pi/32), sin (3*pi/32), sin (5*pi/32), sin (7*pi/32), sin (9*pi/32), sin (11*pi/32), sin (13*pi/32), sin (15*pi/32)
		hwindow = new double[len];
		double f = Math.PI/(2.0*len);
		double sqrt1_2=Math.sqrt(0.5);
		if      (mode < 0) mode =0;
		else if (mode > 2) mode = 2;
		if (mode ==0){
			for (int i = 0; i < len; i++ ) hwindow[i] = sqrt1_2;
		} else if (mode ==1){
			for (int i = 0; i < len; i++ ) hwindow[i] = Math.sin(f*(i+0.5));
		} else if (mode ==2){
			double s;
			for (int i = 0; i < len; i++ ) {
				s = Math.sin(f*(i+0.5));
				hwindow[i] = Math.sin(Math.PI*s*s/2);
			}
		}
		set_fold_2d(len);
		set_unfold_2d(len);
	}

	// get current LT window as a 2d tile (2*size * 2*size)
	public double [] getWin2d(){
		int size = this.hwindow.length;
		int size2 = 2*size;
		double [] rslt = new double [4*size*size];
		for (int i=0; i< size; i++){
			int ia0 = i * size2;
			int ia1 = (size2 - i -1) * size2;
			for (int j=0; j< size; j++){
				double d = this.hwindow[i] * this.hwindow[j];
				rslt [ia0 + j] = d;
				rslt [ia1 + j] = d;
				rslt [ia0 + size2 - j -1] = d;
				rslt [ia1 + size2 - j -1] = d;
			}
		}
		return rslt;
	}


	// Convert 2nx2n overlapping tile to n*n for dct-iv
	public double [] fold_tile(double [] x, int mode) { // x should be 2n*2n
		return fold_tile(x, 1 << (ilog2(x.length/4)/2));
	}
	public double [] fold_tile(double [] x, int n, int mode) { // x should be 2n*2n
		return fold_tile(x,n, mode,this.fold_k);
	}

	public double [] fold_tile(
			double [] x,
			int n,
			int mode, //////
			double [][][] fold_k
			) { // x should be 2n*2n
		double [] y = new double [n*n];
		for (int i = 0; i<y.length;i++) {
			y[i] = 0;
			for (int k = 0; k < 4; k++){
				y[i] += x[fold_index[i][k]] * fold_k[mode][i][k];
			}
		}
		return y;
	}

	public double [] fold_tile_debug(
			double [] x,
			int n,
			int mode, //////
			double [][][] fold_k
			) { // x should be 2n*2n
		System.out.println("fold_tile_debug, mode = "+mode);
		double [] y = new double [n*n];
		for (int i = 0; i<y.length;i++) {
			y[i] = 0;
			System.out.print(String.format("%2d: ",i));
			for (int k = 0; k < 4; k++){
				y[i] += x[fold_index[i][k]] * fold_k[mode][i][k];
				System.out.print(String.format("(%f * %f) ", x[fold_index[i][k]], fold_k[mode][i][k]));
				if (k < 3) System.out.print("+ ");
			}
			System.out.println(String.format("= %f", y[i]));
		}
		return y;
	}



	public double [] unfold_tile(
			double [] x,  // x should be n*n
			int mode)
	{
		return unfold_tile(x, 1 << (ilog2(x.length)/2), mode);

	}
	public double [] unfold_tile(
			double [] x,  // x should be 2n*2n
			int n,
			int mode)
	{
		double [] y = new double [4*n*n];
		for (int i = 0; i<y.length;i++) {
			y[i] = unfold_k[mode][i]* x[unfold_index[i]];
		}
		return y;
	}


	public double [] dctii_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (CII==null){
			setup_CII(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= CII[t][i][j]*x[j];
			}
		}
		return y;
	}

	public double [] dctiie_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (CIIe==null){
			setup_CIIe(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= CIIe[t][i][j]*x[j];
			}
		}
		return y;
	}


	public double [] dctiii_direct(double[] x){
		// CIII=transp(CII)
		int n = x.length;
		int t = ilog2(n)-1;
		if (CII==null){
			setup_CII(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= CII[t][j][i]*x[j];
			}
		}
		return y;
	}

	public double [] dctiiie_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (CIIIe==null){
			setup_CIIIe(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= CIIIe[t][i][j]*x[j];
			}
		}
		return y;
	}

	public double [] dctiv_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (CIV==null){
			setup_CIV(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= CIV[t][i][j]*x[j];
			}
		}
		return y;
	}

	public double [] dstiie_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (SIIe==null){
			setup_SIIe(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= SIIe[t][i][j]*x[j];
			}
		}
		return y;
	}

	public double [] dstiiie_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (SIIIe==null){
			setup_SIIIe(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= SIIIe[t][i][j]*x[j];
			}
		}
		return y;
	}

	public double [] dstiv_direct(double[] x){
		int n = x.length;
		int t = ilog2(n)-1;
		if (SIV==null){
			setup_SIV(N); // just full size
		}
		double [] y = new double[n];
		for (int i = 0; i<n; i++) {
			y[i] = 0.0;
			for (int j = 0; j< n; j++){
				y[i]+= SIV[t][i][j]*x[j];
			}
		}
		return y;
	}
	private void setup_arrays(int maxN){
		if (N >= maxN) return;
		N = maxN;
		int l = ilog2(N)-1;
		CN1 = new double[l][];
		SN1 = new double[l][];

		for (int t = 0; t<CN1.length; t++) {
			int n1 = 2 << t; // for N==3: 2, 4, 8
			double pi_4n=Math.PI/(8*n1); // n1 = n/2
			CN1[t] = new double[n1];
			SN1[t] = new double[n1];
			for (int k=0; k<n1; k++){
				CN1[t][k] = Math.cos((2*k+1)*pi_4n);
				SN1[t][k] = Math.sin((2*k+1)*pi_4n);
			}
		}
	}

	private void setup_CII(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(CII==null) && (CII.length >= l)) return;
		CII = new double[l][][]; // only needed for direct? Assign only when needed?
		for (int t = 0; t<CII.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			CII[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double ej;
			double pi_2n=Math.PI/(2*n);
			for (int j=0;j<n; j++){
				if (j==0) ej= Math.sqrt(0.5);
				else ej = 1.0;
				for (int k = 0; k<n; k++){
					CII[t][j][k] = scale * ej * Math.cos(j*(2*k+1)*pi_2n);
				}
			}
		}
	}

	private void setup_CIIe(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(CIIe==null) && (CIIe.length >= l)) return;
		CIIe = new double[l][][]; // only needed for direct? Assign only when needed?
		for (int t = 0; t<CIIe.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			CIIe[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double pi_2n=Math.PI/(2*n);
			for (int j=0;j<n; j++){
				for (int k = 0; k<n; k++){
					CIIe[t][j][k] = scale * Math.cos(j*(2*k+1)*pi_2n);
				}
			}
		}
	}

	private void setup_CIIIe(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(CIIIe==null) && (CIIIe.length >= l)) return;
		CIIIe = new double[l][][]; // only needed for direct? Assign only when needed?
		for (int t = 0; t<CIIIe.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			CIIIe[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double ej;
			double pi_2n=Math.PI/(2*n);
			for (int j=0;j < n; j++){
				if ((j==0) || (j == (n-1))) ej= 0.5; // Math.sqrt(0.5);
//				if (j==0) ej= 0.5; // Math.sqrt(0.5); Should it be this? https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-III
				else ej = 1.0;
				for (int k = 0; k<n; k++){
					CIIIe[t][k][j] = scale * ej * Math.cos(j*(2*k+1)*pi_2n);
				}
			}
		}
	}


	private void setup_CIV(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(CIV==null) && (CIV.length >= l)) return;
		CIV = new double[l][][];
		for (int t = 0; t<CIV.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			CIV[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double pi_4n=Math.PI/(4*n);
			for (int j=0;j<n; j++){
				for (int k = 0; k < j; k++){
					CIV[t][j][k] = CIV[t][k][j];
				}
				for (int k = j; k<n; k++){
					CIV[t][j][k] = scale *  Math.cos((2*j+1)*(2*k+1)*pi_4n);
				}
			}
		}
	}

	private void setup_SIIe(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(SIIe==null) && (SIIe.length >= l)) return;
		SIIe = new double[l][][]; // only needed for direct? Assign only when needed?
		for (int t = 0; t<SIIe.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			SIIe[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double pi_2n=Math.PI/(2*n);
			for (int j=0;j<n; j++){
				for (int k = 0; k<n; k++){
					SIIe[t][j][k] = scale * Math.sin((j+1)*(2*k+1)*pi_2n);
				}
			}
		}
	}

	private void setup_SIIIe(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(SIIIe==null) && (SIIIe.length >= l)) return;
		SIIIe = new double[l][][]; // only needed for direct? Assign only when needed?
		for (int t = 0; t<SIIIe.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			SIIIe[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double ej;
			double pi_2n=Math.PI/(2*n);
			for (int j=0;j < n; j++){
//				if ((j==0) || (j == (n-1))) ej= 0.5; // Math.sqrt(0.5);
				if (j == (n-1)) ej= 0.5; // Math.sqrt(0.5);
				else ej = 1.0;
				for (int k = 0; k<n; k++){
					SIIIe[t][k][j] = scale * ej * Math.sin((j+1)*(2*k+1)*pi_2n);
				}
			}
		}
	}

	private void setup_SIV(int maxN){
		if (maxN > N) setup_arrays(maxN);
		int l = ilog2(N);
		if (!(SIV==null) && (SIV.length >= l)) return;
		SIV = new double[l][][];
		for (int t = 0; t<SIV.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
			SIV[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double pi_4n=Math.PI/(4*n);
			for (int j=0;j<n; j++){
				for (int k = 0; k < j; k++){
					SIV[t][j][k] = SIV[t][k][j];
				}
				for (int k = j; k<n; k++){
					SIV[t][j][k] = scale *  Math.sin((2*j+1)*(2*k+1)*pi_4n);
				}
			}
		}
	}


	private int ilog2(int n){
		int i;
		for (i=0; n>1; n= n >> 1) i++;
		return i;
	}


	private double [] _dctii_recurs(double[] x){
		int n = x.length;
		if (n ==2) {
			double [] y= {x[0]+x[1],x[0]-x[1]};
			return y;
		}
		int n1 = n >> 1;
		double [] u0 = new double [n1];
		double [] u1 = new double [n1];
		// u = sqrt(2)*Tn(0) * x
		for (int j = 0; j< n1; j++){
			u0[j]=            (x[j] +    x[n-j-1]);
			u1[j]=            (x[j] -    x[n-j-1]);
		}
		double [] v0 = _dctii_recurs(u0);
		double [] v1 = _dctiv_recurs(u1);
		double [] y = new double[n];
		for (int j = 0; j< n1; j++){
			y[2*j] =     v0[j];
			y[2*j+1] =   v1[j];
		}
		return y;
	}

	private double [] _dctiv_recurs(double[] x){
		int n = x.length;
		if (n ==2) {
			double [] y= {COSPI_1_8_SQRT2*x[0] + COSPI_3_8_SQRT2*x[1],
					      COSPI_3_8_SQRT2*x[0] - COSPI_1_8_SQRT2*x[1]};
			return y;
		}

		int n1 = n >> 1;
		int t = ilog2(n1)-1;

		double [] u0 = new double [n1];
		double [] u1 = new double [n1];
		// u = sqrt(2)*Tn(1) * x
		for (int j = 0; j< n1; j++){
			u0[j]=                             ( CN1[t][j] *      x[j] +      SN1[t][j] *         x[n -  j - 1]);
			u1[j]=            ( 1 - 2*(j & 1))*(-SN1[t][n1-j-1] * x[n1-j-1] + CN1[t][n1 - j -1] * x[n1 + j    ]);
		}
		double [] v0 = _dctii_recurs(u0);
		double [] v1 = _dctii_recurs(u1); //both cos-II


		double [] w0 = new double [n1];
		double [] w1 = new double [n1];

		w0[0] =     sqrt2 * v0[0];
		w1[n1-1] =     sqrt2 * v1[0];
		for (int j = 0; j< n1; j++){
			int sgn = (1 - 2* (j & 1));
			if (j > 0)	    w0[j] = v0[j]   - sgn * v1[n1 - j];
			if (j < (n1-1))	w1[j] = v0[j+1] - sgn * v1[n1 - j -1];
		}

		double [] y = new double[n];
		for (int j = 0; j< n1; j++){
			y[2*j] =     w0[j];
			y[2*j+1] =   w1[j];
		}
		return y;
	}

// ======================= preparing for GPU code, replacing *_recurs ===========================
// tested
	public double [] dttt_ivn(double [] x, int mode, int n, boolean gpu){ // mode 0 - dct,dct 1:dst,dct, 2: dct, dst, 3: dst,dst
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = ((mode & 1)!=0)? dst_ivn(line, gpu):dct_ivn(line, gpu);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = ((mode & 2)!=0)? dst_ivn(line, gpu):dct_ivn(line, gpu);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}


	public double [] dct_iin(double[] x, boolean gpu){
		if (x.length > N){
			N = x.length;
		}
		double [] y;
		if (gpu) {
			y=  new double[8];
			_dctii_nrecurs8(x,y);
		} else {
			y=  _dctii_nrecurs8(x);
		}

		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) y[i] *= scale;
		return y;
	}

	public double [] dct_ivn(double[] x, boolean gpu){
		double [] y;
		if (gpu) {
			y=  new double[8];
			_dctiv_nrecurs8(x,y);
		} else {
			y=  _dctiv_nrecurs8(x);
		}

		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) y[i] *= scale;
		return y;
	}

	public double [] dst_ivn(double[] x, boolean gpu){
		double [] xr= new double[x.length];
		int j= x.length-1;
		for (int i=0; i < x.length;i++) xr[i] = x[j--];

		double [] y;
		if (gpu) {
			y=  new double[8];
			_dctiv_nrecurs8(xr,y);
		} else {
			y=  _dctiv_nrecurs8(xr);
		}

		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) {
			y[i] *= scale;
			scale = -scale;
		}
		return y;
	}




	private double [] _dctii_nrecurs2(double[] x){
		double [] y= {x[0]+x[1],x[0]-x[1]};
		System.out.println(String.format("_dctii_nrecurs2: x={%8f,%8f}",x[0],x[1]));
		System.out.println(String.format("_dctii_nrecurs2: y={%8f,%8f}",y[0],y[1]));

		return y;
	}

	private double [] _dctii_nrecurs4(double[] x){
		double [] u0 = new double [2];
		double [] u1 = new double [2];
		for (int j = 0; j < 2; j++){
			u0[j]=            (x[j] +    x[3-j]);
			u1[j]=            (x[j] -    x[3-j]);
		}
		double [] v0 = _dctii_nrecurs2(u0);
		double [] v1 = _dctiv_nrecurs2(u1);
		double [] y = new double[4];
		for (int j = 0; j< 2; j++){
			y[2*j] =     v0[j];
			y[2*j+1] =   v1[j];
		}

		System.out.println(String.format("_dctii_nrecurs4: x={%8f,%8f,%8f,%8f}",x[0],x[1],x[2],x[3]));
		System.out.println(String.format("_dctii_nrecurs4: u0={%8f,%8f}, u1={%8f,%8f}",u0[0],u0[1],u1[0],u1[1]));
		System.out.println(String.format("_dctii_nrecurs4: v0={%8f,%8f}, v1={%8f,%8f}",v0[0],v0[1],v1[0],v1[1]));
		System.out.println(String.format("_dctii_nrecurs4: y={%8f,%8f,%8f,%8f}",y[0],y[1],y[2],y[3]));

		return y;
	}

	private double [] _dctii_nrecurs8(double[] x){
//		int n = x.length; // = 8
//		int n1 = n >> 1; // 4
		double [] u0 = new double [4];
		double [] u1 = new double [4];
		for (int j = 0; j< 4; j++){
			u0[j]=            (x[j] +    x[7-j]);
			u1[j]=            (x[j] -    x[7-j]);
		}
		double [] v0 = _dctii_nrecurs4(u0);
		double [] v1 = _dctiv_nrecurs4(u1);
		double [] y = new double[8];
		for (int j = 0; j< 4; j++){
			y[2*j] =     v0[j];
			y[2*j+1] =   v1[j];
		}
		System.out.println(String.format("_dctii_nrecurs8: x={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]));
		System.out.println(String.format("_dctii_nrecurs8: u0={%8f,%8f,%8f,%8f}, u1={%8f,%8f,%8f,%8f}",u0[0],u0[1],u0[2],u0[3],u1[0],u1[1],u1[2],u1[3]));
		System.out.println(String.format("_dctii_nrecurs8: v0={%8f,%8f,%8f,%8f}, v1={%8f,%8f,%8f,%8f}",v0[0],v0[1],v0[2],v0[3],v1[0],v1[1],v1[2],v1[3]));
		System.out.println(String.format("_dctii_nrecurs8: y={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]));

		return y;
	}


	private double [] _dctiv_nrecurs2(double[] x){
//		int n = x.length; == 2
		double [] y= {COSPI_1_8_SQRT2*x[0] + COSPI_3_8_SQRT2*x[1],
				      COSPI_3_8_SQRT2*x[0] - COSPI_1_8_SQRT2*x[1]};
		System.out.println(String.format("_dctiv_nrecurs2: x={%8f,%8f}",x[0],x[1]));
		System.out.println(String.format("_dctiv_nrecurs2: y={%8f,%8f}",y[0],y[1]));

		return y;
	}

	private double [] _dctiv_nrecurs4(double[] x){
		double [] u0 = new double [2];
		double [] u1 = new double [2];
		for (int j = 0; j< 2; j++){
			u0[j]=                             ( CN1[0][j] *   x[j] +   SN1[0][j] *   x[3 - j]);
			u1[j]=            ( 1 - 2*(j & 1))*(-SN1[0][1-j] * x[1-j] + CN1[0][1-j] * x[2 + j]);
		}
		double [] v0 = _dctii_nrecurs2(u0);
		double [] v1 = _dctii_nrecurs2(u1); //both cos-II


		double [] w0 = new double [2];
		double [] w1 = new double [2];

		w0[0] =     sqrt2 * v0[0];
		w1[1] =     sqrt2 * v1[0];
		for (int j = 0; j< 2; j++){
			int sgn = (1 - 2* (j & 1));
			if (j > 0)	    w0[j] = v0[j]   - sgn * v1[2 - j];
			if (j < (1))	w1[j] = v0[j+1] - sgn * v1[2 - j -1];
		}

		double [] y = new double[4];
		for (int j = 0; j< 2; j++){
			y[2*j] =     w0[j];
			y[2*j+1] =   w1[j];
		}

		System.out.println(String.format("_dctiv_nrecurs4: x={%8f,%8f,%8f,%8f}",x[0],x[1],x[2],x[3]));
		System.out.println(String.format("_dctiv_nrecurs4: u0={%8f,%8f}, u1={%8f,%8f}",u0[0],u0[1],u1[0],u1[1]));
		System.out.println(String.format("_dctiv_nrecurs4: v0={%8f,%8f}, v1={%8f,%8f}",v0[0],v0[1],v1[0],v1[1]));
		System.out.println(String.format("_dctiv_nrecurs4: w0={%8f,%8f}, w1={%8f,%8f}",w0[0],w0[1],w1[0],w1[1]));
		System.out.println(String.format("_dctiv_nrecurs4: y={%8f,%8f,%8f,%8f}",y[0],y[1],y[2],y[3]));
		return y;
	}

	private double [] _dctiv_nrecurs8(double[] x){

		double [] u0 = new double [4];
		double [] u1 = new double [4];
		for (int j = 0; j< 4; j++){
			u0[j]=                             ( CN1[1][j] *   x[j] +     SN1[1][j] *     x[7 - j]);
			u1[j]=            ( 1 - 2*(j & 1))*(-SN1[1][3-j] * x[3 - j] + CN1[1][3 - j] * x[4 + j]);
		}

		double [] v0 = _dctii_nrecurs4(u0);
		double [] v1 = _dctii_nrecurs4(u1); //both cos-II

		double [] w0 = new double [4];
		double [] w1 = new double [4];

		w0[0] =     sqrt2 * v0[0];
		w1[3] =     sqrt2 * v1[0];
		for (int j = 0; j< 4; j++){
			int sgn = (1 - 2* (j & 1)); //TODO:  optimize it out
			if (j > 0)	w0[j] = v0[j]   - sgn * v1[4 - j];
			if (j < 3)	w1[j] = v0[j+1] - sgn * v1[3 - j];
		}

		double [] y = new double[8];
		for (int j = 0; j< 4; j++){
			y[2*j] =     w0[j];
			y[2*j+1] =   w1[j];
		}
		System.out.println(String.format("_dctiv_nrecurs8: x={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]));
		System.out.println(String.format("_dctiv_nrecurs8: u0={%8f,%8f,%8f,%8f}, u1={%8f,%8f,%8f,%8f}",u0[0],u0[1],u0[2],u0[3],u1[0],u1[1],u1[2],u1[3]));
		System.out.println(String.format("_dctiv_nrecurs8: v0={%8f,%8f,%8f,%8f}, v1={%8f,%8f,%8f,%8f}",v0[0],v0[1],v0[2],v0[3],v1[0],v1[1],v1[2],v1[3]));
		System.out.println(String.format("_dctiv_nrecurs8: w0={%8f,%8f,%8f,%8f}, w1={%8f,%8f,%8f,%8f}",w0[0],w0[1],w0[2],w0[3],w1[0],w1[1],w1[2],w1[3]));
		System.out.println(String.format("_dctiv_nrecurs8: y={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]));
		return y;
	}


// ++++++++++++++++++++++++++++ More comparison with GPU ++++++++++++++++++++++++++++++++++
//
//	double COSPI_1_8_SQRT2 = 1.306563;
//	double COSPI_3_8_SQRT2 = 0.541196;
	double SQRT_2 = 1.414214;
	double SQRT1_2 = 0.707107;
	double SQRT1_8 = 0.353553;
	double [] COSN1 = {0.980785,0.831470};
	double [] COSN2 = {0.995185,0.956940,0.881921,0.773010};
	double [] SINN1 = {0.195090,0.555570};
	double [] SINN2 = {0.098017,0.290285,0.471397,0.634393};


	void _dctii_nrecurs8( double [] x, double [] y) // x,y point to 8-element arrays each
	{
		double u00= (x[0] + x[7]);
		double u10= (x[0] - x[7]);

		double u01= (x[1] + x[6]);
		double u11= (x[1] - x[6]);

		double u02= (x[2] + x[5]);
		double u12= (x[2] - x[5]);

		double u03= (x[3] + x[4]);
		double u13= (x[3] - x[4]);

//		double v00, v01, v02, v03, v10, v11, v12, v13;
//		_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03);

		double w00= u00 + u03;
		double w10= u00 - u03;

		double w01= (u01 + u02);
		double w11= (u01 - u02);

		double v00= w00 + w01;
		double v02= w00 - w01;
		double v01= COSPI_1_8_SQRT2 * w10 + COSPI_3_8_SQRT2 * w11;
		double v03= COSPI_3_8_SQRT2 * w10 - COSPI_1_8_SQRT2 * w11;

//		_dctiv_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

/*
		double [] u0 = new double [2];
		double [] u1 = new double [2];
		for (int j = 0; j< 2; j++){
			u0[j]=                             ( CN1[0][j] *   x[j] +   SN1[0][j] *   x[3 - j]);
			u1[j]=            ( 1 - 2*(j & 1))*(-SN1[0][1-j] * x[1-j] + CN1[0][1-j] * x[2 + j]);
		}
		double [] v0 = _dctii_nrecurs2(u0);
		double [] v1 = _dctii_nrecurs2(u1); //both cos-II


		double [] w0 = new double [2];
		double [] w1 = new double [2];

		w0[0] =     sqrt2 * v0[0];
		w1[1] =     sqrt2 * v1[0];
		for (int j = 0; j< 2; j++){
			int sgn = (1 - 2* (j & 1));
			if (j > 0)	    w0[j] = v0[j]   - sgn * v1[2 - j];
			if (j < (1))	w1[j] = v0[j+1] - sgn * v1[2 - j -1];
		}

		double [] y = new double[4];
		for (int j = 0; j< 2; j++){
			y[2*j] =     w0[j];
			y[2*j+1] =   w1[j];
		}
 */


		double w20=            ( COSN1[0] * u10 + SINN1[0] * u13);
		double w30=            (-SINN1[1] * u11 + COSN1[1] * u12);

		double w21=            ( COSN1[1] * u11 + SINN1[1] * u12);
		double w31=           -(-SINN1[0] * u10 + COSN1[0] * u13);
	//u->t, v ->z, t0->w2, t1->w3
//		double z00, z01, z10,z11;
// TODO: change w20<->w02
//		_dctii_nrecurs2(u00, u01, &v00, &v01);
		double z00= w20 + w21;
		double z01= w20 - w21;

//		_dctii_nrecurs2(u10, u11, &v10, &v11);
		double z10= w30 + w31;
		double z11= w30 - w31;

		double v10 = SQRT_2 * z00;
		double v11 = z01 - z11;

		double v12 = z01 + z11;
		double v13 = SQRT_2 * z10;

		y[0] =   v00;
		y[1] =   v10;

		y[2] =   v01;
		y[3] =   v11;

		y[4] =   v02;
		y[5] =   v12;

		y[6] =   v03;
		y[7] =   v13;
		System.out.println(String.format("_dctii_nrecurs8-gpu: x={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]));
		System.out.println(String.format("_dctii_nrecurs8-gpu: u0={%8f,%8f,%8f,%8f}, u1={%8f,%8f,%8f,%8f}",u00,u01,u02,u03,u10,u11,u12,u13));
		System.out.println(String.format("_dctii_nrecurs8-gpu: v0={%8f,%8f,%8f,%8f}, v1={%8f,%8f,%8f,%8f}",v00,v01,v02,v03,v10,v11,v12,v13));
		System.out.println(String.format("_dctii_nrecurs8-gpu: w0={%8f,%8f,%8f,%8f}, w1={%8f,%8f,%8f,%8f}",w00,w10,w20,w30,w01,w11,w21,w31));
		System.out.println(String.format("_dctii_nrecurs8-gpu: y={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]));

	}

	void _dctiv_nrecurs8( double [] x, double [] y) // x,y point to 8-element arrays each
	{
		double u00=            ( COSN2[0] * x[0] + SINN2[0] * x[7]);
		double u10=            (-SINN2[3] * x[3] + COSN2[3] * x[4]);

		double u01=            ( COSN2[1] * x[1] + SINN2[1] * x[6]);
		double u11=           -(-SINN2[2] * x[2] + COSN2[2] * x[5]);

		double u02=            ( COSN2[2] * x[2] + SINN2[2] * x[5]);
		double u12=            (-SINN2[1] * x[1] + COSN2[1] * x[6]);

		double u03=            ( COSN2[3] * x[3] + SINN2[3] * x[4]);
		double u13=           -(-SINN2[0] * x[0] + COSN2[0] * x[7]);

//		_dctii_nrecurs4(u00, u01, u02, u03, &v00, &v01, &v02, &v03); // ua - u inside first _dctii_nrecurs4, ub - inside second
		double ua00= u00 + u03;
		double ua10= u00 - u03;

		double ua01= u01 + u02;
		double ua11= u01 - u02;

		double v00= ua00 + ua01;
		double v02= ua00 - ua01;

		double v01= COSPI_1_8_SQRT2 * ua10 + COSPI_3_8_SQRT2 * ua11;
		double v03= COSPI_3_8_SQRT2 * ua10 - COSPI_1_8_SQRT2 * ua11;
//		_dctii_nrecurs4(u10, u11, u12, u13, &v10, &v11, &v12, &v13);

		System.out.println(String.format("_dctiv_nrecurs8-gpu (_dctii_nrecurs4(x): x={%8f,%8f,%8f,%8f}",u10, u11, u12, u13));

		double ub00= u10 + u13;
		double ub10= u10 - u13;

		double ub01= u11 + u12;
		double ub11= u11 - u12;
		System.out.println(String.format("_dctiv_nrecurs8-gpu: ub0={%8f,%8f}, ub1={%8f,%8f}",ub00,ub01,ub10,ub11));

		double vb00= ub00 + ub01;
		double vb01= ub00 - ub01;

		double vb10= COSPI_1_8_SQRT2*ub10 + COSPI_3_8_SQRT2*ub11;
		double vb11= COSPI_3_8_SQRT2*ub10 - COSPI_1_8_SQRT2*ub11;

		System.out.println(String.format("_dctiv_nrecurs8-gpu: vb0={%8f,%8f}, vb1={%8f,%8f}",vb00,vb01,vb10,vb11));

		double v10 = vb00;
		double v11 = vb10;
		double v12 = vb01;
		double v13 = vb11;
		System.out.println(String.format("_dctiv_nrecurs8-gpu: yb={%8f,%8f,%8f,%8f}",vb00, vb10, vb01, vb11));

		// j == 0
		y[0] =  SQRT_2 * v00;    // w0[0];
		y[1] =  v01 -  v13;    // w1[0];
		// j == 1
		y[2] =  v01 +  v13;    // w0[1];
		y[3] =  v02 +  v12;    // w1[1];
		// j == 2
		y[4] =  v02 -  v12;    // w0[2];
		y[5] =  v03 -  v11;    // w1[2]; - same as y[3]
		// j == 3
		y[6] =  v03 +  v11;    // w0[3];
		y[7] =  SQRT_2 * v10;    // w1[3];
		System.out.println(String.format("_dctiv_nrecurs8-gpu: x={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]));
		System.out.println(String.format("_dctiv_nrecurs8-gpu: u0={%8f,%8f,%8f,%8f}, u1={%8f,%8f,%8f,%8f}",u00,u01,u02,u03,u10,u11,u12,u13));
		System.out.println(String.format("_dctiv_nrecurs8-gpu: v0={%8f,%8f,%8f,%8f}, v1={%8f,%8f,%8f,%8f}",v00,v01,v02,v03,v10,v11,v12,v13));
//		System.out.println(String.format("_dctiv_nrecurs8-gpu: w0={%8f,%8f,%8f,%8f}, w1={%8f,%8f,%8f,%8f}",w00,w10,w20,w30,w01,w11,w21,w31));
		System.out.println(String.format("_dctiv_nrecurs8-gpu: y={%8f,%8f,%8f,%8f,%8f,%8f,%8f,%8f}",y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]));
	}



}
