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
	                               // Second index(0..3) - item to add (2 vertiacl, 2 - horizontal)  
	double [][] fold_k = null; // Matching fold_index items. Each is a product of 2 window coefficients and sign  
	int    [] unfold_index = null;  // index  for each element of idct(2nx2n)
	double [] unfold_k = null; // Matching unfold_index items. Each is a product of 2 window coefficients and sign  
	
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
	//                2 - sign of the term
	private int [][] get_fold_indices(int x, int n){
		int n1 = n>>1;
		int [][] ind = new int[2][3];
		if (x <n1) {
			ind[0][0] = n + n1 - x - 1; // -cR
			ind[0][1] = n1     + x;
			ind[0][2] = -1;
			ind[1][0] = n + n1 + x;     // -d
			ind[1][1] = n1     - x - 1;
			ind[1][2] = -1;
		} else {
			x-=n1;
			ind[0][0] =          x;     // +a
			ind[0][1] =          x;
			ind[0][2] = 1;
			ind[1][0] = n      - x - 1; // -bR
			ind[1][1] = n      - x - 1;
			ind[1][2] = -1;
		}
		
		return ind;
	}
	// is called when window is set
	private void set_fold_2d(int n){ // n - DCT and window size
		if ((fold_index != null) && (fold_index.length == n*n)) return;
		fold_index = new int[n*n][4];
		fold_k =     new double[n*n][4];
		int []    vert_ind = new int[2];
		double [] vert_k = new double[2];
		int    [] hor_ind = new int[2];
		double [] hor_k = new double[2];
		int [][] fi;
		int n2 = 2*n;
		for (int i = 0; i < n; i++ ){
			fi = get_fold_indices(i,n);
			vert_ind[0] = fi[0][0];
			vert_ind[1] = fi[1][0];
			vert_k[0] =   fi[0][2] * hwindow[fi[0][1]];
			vert_k[1] =   fi[1][2] * hwindow[fi[1][1]];
			for (int j = 0; j < n; j++ ){
				fi = get_fold_indices(j,n);
				hor_ind[0] = fi[0][0];
				hor_ind[1] = fi[1][0];
				hor_k[0] =   fi[0][2] * hwindow[fi[0][1]];
				hor_k[1] =   fi[1][2] * hwindow[fi[1][1]];
				int indx = n*i + j;
				for (int k = 0; k<4;k++) {
					fold_index[indx][k] = n2 * vert_ind[(k>>1) & 1] + hor_ind[k & 1];
					fold_k[indx][k] =     vert_k[(k>>1) & 1] * hor_k[k & 1]; 
				}
			}
		}
		if (n < 8) {
			for (int i = 0; i < n; i++ ){
				fi = get_fold_indices(i,n);
				System.out.println(i+"->"+String.format("[%2d % 2d % 2d] [%2d %2d %2d] %f %f",
						fi[0][0],fi[0][1],fi[0][2],
						fi[1][0],fi[1][1],fi[1][2], hwindow[fi[0][1]], hwindow[fi[1][1]]));
			}
		}
		
	}
	// return index+1 and sign for 1-d imdct. x is index (0..2*n-1) of the imdct array, value is sign * (idct_index+1),
	// where idct_index (0..n-1) is index in the dct-iv array
	private int get_unfold_index(int x, int n){
		int n1 = n>>1;
		int segm = x / n1;
		x = x % n1;
		switch (segm){
		case 0: return 1+ (x + n1);
		case 1: return -(n - x);
		case 2: return -(n1 - x);
		case 3: return -(1 + x);
		}
		return 0; //should never happen
	}
	private void set_unfold_2d(int n){ // n - DCT size
		if ((unfold_index != null) && (unfold_index.length == 4*n*n)) return;
		unfold_index = new int[4*n*n];
		unfold_k = new double[4*n*n];
		int n2 = 2*n;
		for (int i = 0; i < 2*n; i++ ){
			int index_vert = get_unfold_index(i,n);
			double k_vert = hwindow[(i < n)?i:n2 -i -1];
			if (index_vert <0 ){
				k_vert = -k_vert;
				index_vert = -index_vert;
			}
			index_vert --;
			index_vert *= n;
			
			for (int j = 0; j < 2*n; j++ ){
				int index_hor = get_unfold_index(j,n);
				double k_hor = hwindow[(j < n)?j:n2 -j -1];
				if (index_hor <0 ){
					k_hor = -k_hor;
					index_hor = -index_hor;
				}
				index_hor --; // pass 1 to next
//				unfold_index1[n2*i+j]=sgn_vert*sgn_hor*(index_vert+index_hor); // should never be 0
				unfold_index[n2*i+j]=(index_vert+index_hor);
				unfold_k[n2*i+j]=k_vert*k_hor;
				
				if (n < 8) System.out.print(String.format("%4d", unfold_index[n2*i+j]));
				
			}
			if (n < 8) System.out.println();
		}
		if (n < 8) {
			for (int i = 0; i < 2*n; i++ ){
				System.out.println(i+"->"+get_unfold_index(i,n));
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
		return dttt_iie(x, 1 << (ilog2(x.length)/2)); 
	}

	public double [] dttt_iie(double [] x, int n){
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = dctiie_direct(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose 
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = dctiie_direct(line);
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
		return dttt_iiie(x, 1 << (ilog2(x.length)/2)); 
	}

	public double [] dttt_iiie(double [] x, int n){
		double [] y = new double [n*n];
		double [] line = new double[n];
		// first (horizontal) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = dctiiie_direct(line);
			for (int j=0; j < n;j++) y[j*n+i] =line[j]; // transpose 
		}
		// second (vertical) pass
		for (int i = 0; i<n; i++){
			System.arraycopy(y, n*i, line, 0, n);
			line = dctiiie_direct(line);
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
	public void set_window(int mode, int len){
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
	
	// Convert 2nx2n overlapping tile to n*n for dct-iv
	public double [] fold_tile(double [] x) { // x should be 2n*2n
		return fold_tile(x, 1 << (ilog2(x.length/4)/2));
	}
	public double [] fold_tile(double [] x, int n) { // x should be 2n*2n
		double [] y = new double [n*n];
		for (int i = 0; i<y.length;i++) {
			y[i] = 0;
			for (int k = 0; k < 4; k++){
				y[i] += x[fold_index[i][k]] * fold_k[i][k];
			}
		}
		return y;
	}

	public double [] unfold_tile(double [] x) { // x should be n*n
		return unfold_tile(x, 1 << (ilog2(x.length)/2));
		
	}
	public double [] unfold_tile(double [] x, int n) { // x should be 2n*2n
		double [] y = new double [4*n*n];
		for (int i = 0; i<y.length;i++) {
			y[i] = unfold_k[i]* x[unfold_index[i]];
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
//			double ej;
			double pi_2n=Math.PI/(2*n);
			for (int j=0;j<n; j++){
//				if (j==0) ej= Math.sqrt(0.5);
//				else ej = 1.0;
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
//				if (j==0) ej= 0.5; // Math.sqrt(0.5);
				else ej = 1.0;
				for (int k = 0; k<n; k++){
//					CIIIe[t][j][k] = scale * ej * Math.cos(j*(2*k+1)*pi_2n);  
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
