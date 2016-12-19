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
	double [][][] CII=null;
	double [][][] CIV=null;
	double [][][] SIV=null;
	double [][] CN1=null; 
	double [][] SN1=null; 
	double      COSPI_1_8_SQRT2 = Math.cos(Math.PI/8)*Math.sqrt(2.0);
	double      COSPI_3_8_SQRT2 = Math.cos(3*Math.PI/8)*Math.sqrt(2.0);
	double sqrt2 = Math.sqrt(2.0);
	double sqrt1_2 = 1/sqrt2;
	double [] hwindow = null; // half window
	
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
		double [] y=  _dctiv_recurs(x);
		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) {
			y[i] *= scale;
			scale = -scale; 
		}
		return y;
	}
	
	public double [] mdct_fold (double [] x ) { // x is 2n long
		int n = x.length/2;
		int n1 = n/2;
		int n3 = n+n1;
		double [] x1 = new double [n];
		for (int i=0; i<n1; i++){
			x1[i]=    -hwindow[n1 + i]* x[n3 -1 -i] - hwindow[n1 - i -1] * x[n3     + i]; 
			x1[i+n1]=  hwindow[     i]* x[       i] - hwindow[n  - i -1] * x[n  - 1 - i]; 
		}
		return x1;
	}

	public double [] mdct_unfold (double [] x) { // x is 2n long
		int n = x.length;
		int n1 = n/2;
		int n3 = n+n1;
		double [] x1 = new double [2*n];
		for (int i=0; i<n1; i++){
			x1[i]=           x[n1+i];
			x1[n - i - 1] =  x1[i]; 
			x1[n3 + i]=     -x[i];
			x1[n1 - i - 1] =-x1[i]; 
		}
		return x1;
	}
	
	
	
	public double [] mdct_2d(double [] x){
		return mdct_2d(x, 0, 1 << (ilog2(x.length/4)/2)); 
	}

	public double [] mdct_2d(double [] x, int mode){
		return mdct_2d(x, mode, 1 << (ilog2(x.length/4)/2)); 
	}

	public double [] mdct_2d (double [] x, int mode, int n) { // x is 2n*2n long
//		int n = 1 << (ilog2(x.length/4)/2);
		int n2 = 2*n;
		double [] transp = new double [2*n*n];
		double [] y =      new double [n*n];
		double [] line2 = new double[n*2];
		double [] line = new double[n*2];
		// first (horizontal) pass
		for (int i = 0; i < n2; i++){
			System.arraycopy(x, n2*i, line2, 0, n2);
			line = mdct_fold(line2);
			line = ((mode & 1)!=0)? dst_iv(line):dct_iv(line);
			for (int j=0; j < n;j++) transp[j*n2+i] =line[j]; // transpose 
		}
		// second (vertical) pass
		for (int i = 0; i < n2; i++){
			System.arraycopy(transp, n2*i, line2, 0, n2);
			line = mdct_fold(line2);
			line = ((mode & 2)!=0)? dst_iv(line):dct_iv(line);
			System.arraycopy(line, 0, y, n*i, n);
		}
		return y;
	}

	public double [] imdct_2d(double [] x){
		return imdct_2d(x, 1 << (ilog2(x.length)/2)); 
	}

	public double [] imdct_2d (double [] x, int n) { // x is n*n long
		int n2 = 2*n;
		double [] transp = new double [2*n*n];
		double [] y =      new double [4*n*n];
		double [] line2 = new double[n*2];
		double [] line = new double[n*2];
		// first (horizontal) pass
		for (int i = 0; i < n2; i++){
			System.arraycopy(x, n*i, line, 0, n);
			line = dct_iv(line);
			line2 = mdct_unfold(line);
			for (int j=0; j < n2;j++) transp[j*n+i] =line2[j]; // transpose 
		}
		// second (vertical) pass
		for (int i = 0; i < n2; i++){
			System.arraycopy(transp, n*i, line, 0, n);
			line = dct_iv(line);
			line2 = mdct_unfold(line);
			System.arraycopy(line2, 0, y, n2*i, n2);
		}
		return y;
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
	
	public void set_window(){
		set_window(0);
	}
	public void set_window(int mode){
		set_window(mode, N);
	}
	public void set_window(int mode, int len){
		hwindow = new double[len];
		double f = Math.PI/(2.0*len);

		if (mode ==0){
			for (int i = 0; i < len; i++ ) hwindow[i] = Math.sin(f*(i+0.5));
		} else { // add more types?
			double s;
			for (int i = 0; i < len; i++ ) {
				s = Math.sin(f*(i+0.5));
				hwindow[i] = Math.sin(Math.PI*s*s);
			}
		}
	}
	// Convert 2nx2n overlapping tile to n*n for dct-iv
	public double [] fold_tile(double [] x) { // x should be 2n*2n
		return fold_tile(x, 1 << (ilog2(x.length/4)/2));
	}
	public double [] fold_tile(double [] x, int n) { // x should be 2n*2n
		double [] y = new double [n*n];
		for (int i = 0; i<y.length;i++) y[i] = 0;
		int n1 = n/2;
		int n2 = 2*n;
		for (int tile_y_v=0; tile_y_v<2; tile_y_v++){     // 2 rows of y tiles
			for (int tile_y_h=0; tile_y_h<2; tile_y_h++){ // 2 columns of y tiles
				int start_y_addr = n*n1*tile_y_v + n1*tile_y_h; //atart address in the aoutput array
				int start_x_tl_addr = (2*n*n) * (1 - tile_y_v) + n * (1 - tile_y_h); // address of the top left corner of a group of 4 tiles 
				for (int tile_x_v=0; tile_x_v<2; tile_x_v++){     // 2 rows of x tiles (contributing to the same y tile)
					for (int tile_x_h=0; tile_x_h<2; tile_x_h++){ // 2 columns of x tiles (contributing to the same y tile)
						int dir_x = ((tile_y_h ^ tile_x_h) !=0)? 1 : -1;
						int dir_y = ((tile_y_v ^ tile_x_v) !=0)? 1 : -1;
						int start_x_addr =  start_x_tl_addr +
								            tile_x_v * n *n + // 2n * n/2
								            tile_x_h * n1 +
								            ((dir_y < 0)?(n1-1)*2*n : 0)+
								            ((dir_x < 0)?(n1-1) : 0);
						int dir_window_vert =  (tile_x_v > 0)? -1 : +1; // same for any tile_y_*
						int dir_window_hor =   (tile_x_h > 0)? -1 : +1;
						int start_window_vert =  (tile_y_v > 0)? ((tile_x_v > 0)? n-1: 0 ):((tile_x_v > 0)? n1-1: n1);
						int start_window_hor =   (tile_y_h > 0)? ((tile_x_h > 0)? n-1: 0 ):((tile_x_h > 0)? n1-1: n1);
						
						for (int i = 0; i < n1; i++){             // n1 rows in each y tile
							for (int j = 0; j < n1; j++){         // n1 columns in each y tile
								y[start_y_addr+ n*i+j] += hwindow[start_window_vert + dir_window_vert * i] *
														  hwindow[start_window_hor +  dir_window_hor * j] *
														  x[start_x_addr + n2 * dir_y * i + dir_x * j];
							}
						}
					}
				}
			}

		}
		return y;
	}

	public double [] unfold_tile(double [] x) { // x should be n*n
		return fold_tile(x, 1 << (ilog2(x.length)/2));
		
	}
	public double [] unfold_tile(double [] x, int n) { // x should be 2n*2n
		double [] y = new double [4*n*n];
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
					SIV[t][j][k] = scale *  Math.cos((2*j+1)*(2*k+1)*pi_4n);
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
		System.out.println("_dctii_recurs: n="+n);
		if (n ==2) {
			double [] y= {x[0]+x[1],x[0]-x[1]};
			System.out.println("_dctii_recurs(2): "); 
			for (int j = 0; j< n; j++){
				System.out.print("--x["+j+"]="+x[j]+" ");
			}		
			System.out.println();
			for (int j = 0; j< n; j++){
				System.out.print("--y["+j+"]="+y[j]+" ");
			}		
			System.out.println();
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
		
		for (int j = 0; j< n1; j++){
			System.out.println("u0["+j+"]="+u0[j]+", u1["+j+"]="+u1[j]);
		}		
		
		
		double [] v0 = _dctii_recurs(u0);
		double [] v1 = _dctiv_recurs(u1);
		for (int j = 0; j< n1; j++){
			System.out.println("_dctii_recurs(): v0["+j+"]="+v0[j]+", v1["+j+"]="+v1[j]);
		}		
		
		double [] y = new double[n];
		for (int j = 0; j< n1; j++){
			y[2*j] =     v0[j];
			y[2*j+1] =   v1[j];
		}
		
		
		for (int j = 0; j< n; j++){
			System.out.println("_dctii_recurs(): y["+j+"]="+y[j]);
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
