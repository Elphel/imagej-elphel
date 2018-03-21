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


}
