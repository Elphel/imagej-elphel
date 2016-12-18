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
	int N = 8;
	double [][][] CII;
	double [][][] CIV;
	double [][] CN1; 
	double [][] SN1; 
//	double [][] C_N1M1; 
//	double [][] S_N1M1; 
	double      COSPI_1_8_SQRT2 = Math.cos(Math.PI/8)*Math.sqrt(2.0);
	double      COSPI_3_8_SQRT2 = Math.cos(3*Math.PI/8)*Math.sqrt(2.0);
	double sqrt2 = Math.sqrt(2.0);
	double sqrt1_2 = 1/sqrt2;
	public DttRad2 (int maxN){ // n - maximal
		N=maxN;
		CII = null; // only needed for direct transforms. Assign when first used
		CIV = null; // same

		CN1 = new double[ilog2(N)-1][];
		SN1 = new double[ilog2(N)-1][];
//		C_N1M1 = new double[ilog2(N)-1][];
//		S_N1M1 = new double[ilog2(N)-1][];
		
		for (int t = 0; t<CN1.length; t++) {
			int n1 = 2 << t; // for N==3: 2, 4, 8
			double pi_4n=Math.PI/(8*n1); // n1 = n/2
//			double pi_2n=Math.PI/(4*n1); // n1 = n/2
			CN1[t] = new double[n1];
			SN1[t] = new double[n1];
//			C_N1M1[t] = new double[n1-1];
//			S_N1M1[t] = new double[n1-1];
			for (int k=0; k<n1; k++){
				CN1[t][k] = Math.cos((2*k+1)*pi_4n);
				SN1[t][k] = Math.sin((2*k+1)*pi_4n);
			}
//			for (int k=1; k<n1; k++){
//				C_N1M1[t][k-1] = Math.cos(k*pi_2n);
//				S_N1M1[t][k-1] = Math.sin(k*pi_2n);
//			}
		} 

	}
	
	public void setup_CII(int maxN){
		CII = new double[ilog2(N)][][]; // only needed for direct? Assign only when needed?
		for (int t = 0; t<CII.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
//			System.out.println("t="+t+", n="+n);
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
	public void setup_CIV(int maxN){
		CIV = new double[ilog2(N)][][];
		for (int t = 0; t<CIV.length; t++) {
			int n = 2 << t; // for N==3: 2, 4, 8
//			System.out.println("t="+t+", n="+n);
			CIV[t] = new double[n][n];
			double scale = Math.sqrt(2.0/n);
			double pi_4n=Math.PI/(4*n);
			for (int j=0;j<n; j++){
				for (int k = 0; k < j; k++){
					CIV[t][j][k] = CIV[t][k][j];  
				}
				for (int k = j; k<n; k++){
					CIV[t][j][k] = scale *  Math.cos((2*j+1)*(2*k+1)*pi_4n);
//                    if (t<2) System.out.println("CIV["+t+"]["+j+"]["+k+"]="+CIV[t][j][k]);
				}
			}
		}
	}	
	
	public int ilog2(int n){
		int i;
		for (i=0; n>1; n= n >> 1) i++;
		return i;
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
	public double [] dctii_recurs(double[] x){
		double [] y=  _dctii_recurs(x);
		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) y[i] *= scale;
		return y;
	}

	public double [] dctiv_recurs(double[] x){
		double [] y=  _dctiv_recurs(x);
		double scale = 1.0/Math.sqrt(x.length);
		for (int i = 0; i < y.length ; i++) y[i] *= scale;
		return y;
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
			for (int j = 0; j< n; j++){
//				System.out.println("_dctiv_recurs(2): y["+j+"]="+y[j]);
			}		
			
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

		for (int j = 0; j< n1; j++){
//			System.out.println("v0["+j+"]="+v0[j]+", v1["+j+"]="+v1[j]);
		}		
		
		w0[0] =     sqrt2 * v0[0];
		w1[n1-1] =     sqrt2 * v1[0];
		for (int j = 0; j< n1; j++){
			int sgn = (1 - 2* (j & 1));
			if (j > 0)	    w0[j] = v0[j]   - sgn * v1[n1 - j]; 
			if (j < (n1-1))	w1[j] = v0[j+1] - sgn * v1[n1 - j -1]; 
		}		
		for (int j = 0; j< n1; j++){
//			System.out.println("w0["+j+"]="+w0[j]+", w1["+j+"]="+w1[j]);
		}		
		double [] y = new double[n];
		for (int j = 0; j< n1; j++){
			y[2*j] =     w0[j];
			y[2*j+1] =   w1[j];
		}
		for (int j = 0; j< n; j++){
//			System.out.println("y["+j+"]="+y[j]);
		}		
		return y;
	}
}
