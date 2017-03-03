import Jama.LUDecomposition;
import Jama.Matrix;

public class PolynomialApproximation {
	public int debugLevel=1;
// TODO Move other methods here
	public PolynomialApproximation(){}
	public PolynomialApproximation(int debugLevel){
		this.debugLevel=debugLevel;
	}
	public double [] polynomialApproximation1d(double [][] data, int N){
		double [] S=new double [2*N+1];
		double [] SF=new double [N+1];
		for (int i=0;i<=2*N;i++) S[i]=0.0;
		for (int i=0;i<=N;i++)  SF[i]=0.0;
		for (int i=0;i<data.length;i++){
			double wxn=(data[i].length>2)?data[i][2]:1.0;
			if (wxn>0.0){ // save time on 0.0 that can be used to mask out some samples
				double f=data[i][1];
				double x=data[i][0];
				for (int j=0;j<=N;j++){
					S[j]+=wxn;
					SF[j]+=wxn*f;
					wxn*=x;
				}
				for (int j=N+1;j<2*N;j++){
					S[j]+=wxn;
					wxn*=x;
				}
				S[2*N]+=wxn;
				if (this.debugLevel>1){
					if (data[i].length > 2)
						System.out.println("polynomialApproximation1d() |"+i+"|: x=|"+data[i][0]+"| f(x)=|"+data[i][1]+"| (w=\t|"+data[i][2]+"|\t)");
					else
						System.out.println("polynomialApproximation1d() |"+i+"|: x=|"+data[i][0]+"| f(x)=|"+data[i][1]+"|\t)");
				}
			}
		}
		double [][] aM=new double [N+1][N+1];
		double [][] aB=new double [N+1][1];
		for (int i=0;i<=N; i++) {
			aB[i][0]=SF[i];
			for (int j=0;j<=N;j++) aM[i][j]=S[i+j];
		}
		Matrix M=new Matrix(aM,N+1,N+1);
		Matrix B=new Matrix(aB,N+1,1);
		int N1=N;
		// TODO: use try/catch with solve
		if (this.debugLevel>1){
			System.out.println("polynomialApproximation1d(data,"+N+") M:");
			M.print(10, 5);
			System.out.println("polynomialApproximation1d() B:");
			B.print(10, 5);
		}
//		while (!(new LUDecomposition(M)).isNonsingular() && (N1>0)){
		while (!(new LUDecomposition(M)).isNonsingular() && (N1>=0)){ // make N=0 legal ?
			aM=new double [N1][N1];
			aB=new double [N1][1];
			N1--;
			for (int i=0;i<=N1; i++) {
				aB[i][0]=B.getArray()[i][0];
				for (int j=0;j<=N1;j++) aM[i][j]=M.getArray()[i][j];
			}
			M=new Matrix(aM,N1+1,N1+1);
			B=new Matrix(aB,N1+1,1);
			if (this.debugLevel>1){
				System.out.println("polynomialApproximation1d() Reduced degree: M:");
				M.print(10, 5);
				System.out.println("polynomialApproximation1d() Reduced degree: B:");
				B.print(10, 5);
			}
		}
		double [][] aR=M.solve(B).getArray();
		if (this.debugLevel>1){
			System.out.println("polynomialApproximation1d() solution=");
			M.solve(B).print(10, 12);
		}
		double []result=new double [N+1];
		
		for (int i=0;i<=N;i++) result[i]=(i<=N1)?aR[i][0]:0.0;
		return result;
	}
/**
 * Linear approximates each of 3 functions of 3 variables and finds where they are all zero	
 * @param data: for each sample (1-st index):
 *               0 - {x,y,z}
 *               1 - {f1, f2, f3},
 *               2 - {weight}
 * @return {x0, y0, z0} where A1*x0+B1*y0+C1*z0+D1=0, A2*x0+B2*y0+C2*z0+D2=0, A3*x0+B3*y0+C3*z0+D3=0
 */
	public double [] linear3d(double [][][] data){
		/*
		 * Approximating each of the 3 measured parameters (Far/near, tilt x and tilt y) with linear approximation in the vicinity of the last position
		 * For each parameter F(x,y,z)=A*x + B*y +C*z + D, using Gaussian weight function with sigma= motorsSigma     	    		
		 */
		double [][] aM3=new double [3][3];
		double [][] aB3=new double [3][1];
		for (int parNum=0;parNum<aM3.length;parNum++){
			double S0=0.0,SX=0.0,SY=0.0,SZ=0.0,SXY=0.0,SXZ=0.0,SYZ=0.0,SX2=0.0,SY2=0.0,SZ2=0.0,SF=0.0,SFX=0.0,SFY=0.0,SFZ=0.0;
			for (int nSample=0; nSample< data.length;nSample++){
				if (this.debugLevel>3){
					System.out.println(
							parNum+","+data[nSample][0][0]+","+data[nSample][0][1]+","+data[nSample][0][2]+","+
							data[nSample][1][0]+","+data[nSample][1][1]+","+data[nSample][1][2]+","+data[nSample][2][0]);
				}

	    		//TODO replace with PolynomialApproximation class
				double w=(data[nSample].length>2)?data[nSample][2][0]:1.0;
				double [] xyz=data[nSample][0];
				S0+=w;
				SX+=w*xyz[0];
				SY+=w*xyz[1];
				SZ+=w*xyz[2];
				SXY+=w*xyz[0]*xyz[1];
				SXZ+=w*xyz[0]*xyz[2];
				SYZ+=w*xyz[1]*xyz[2];
				SX2+=w*xyz[0]*xyz[0];
				SY2+=w*xyz[1]*xyz[1];
				SZ2+=w*xyz[2]*xyz[2];
				SF+=w*data[nSample][1][parNum];
				SFX+=w*data[nSample][1][parNum]*xyz[0];
				SFY+=w*data[nSample][1][parNum]*xyz[1];
				SFZ+=w*data[nSample][1][parNum]*xyz[2];
			}
			double [][] aM={
					{SX2,SXY,SXZ,SX},
					{SXY,SY2,SYZ,SY},
					{SXZ,SYZ,SZ2,SZ},
					{SX, SY, SZ, S0}};
			double [][] aB={
					{SFX},
					{SFY},
					{SFZ},
					{SF}};
			Matrix M=new Matrix(aM);
			Matrix B=new Matrix(aB);
			// Check for singular (near-singular) here
			double [] abcd= M.solve(B).getColumnPackedCopy();
			if (this.debugLevel>2){
				System.out.println(parNum+"M:");
				M.print(10, 5);
				System.out.println(parNum+"B:");
				B.print(10, 5);
				System.out.println(parNum+"A:");
				M.solve(B).print(10, 7);
			}
			//believeLast
			aM3[parNum][0]= abcd[0];
			aM3[parNum][1]= abcd[1];
			aM3[parNum][2]= abcd[2];
			aB3[parNum][0]=-abcd[3];
			if (this.debugLevel>1) System.out.println("** "+parNum+": A="+abcd[0]+" B="+abcd[1]+" C="+abcd[2]+" D="+abcd[3]);
		}
		Matrix M3=new Matrix(aM3);
		Matrix B3=new Matrix(aB3);
		double [] result=M3.solve(B3).getColumnPackedCopy();
		if (this.debugLevel>2) {
			System.out.println("M3:");
			M3.print(10, 7);
			System.out.println("B3:");
			B3.print(10, 7);
			System.out.println("A3:");
			M3.solve(B3).print(10, 5);
		}
		return result;

		
	}

	  public double[] quadraticMax2d (double [][][] data){
		  return quadraticMax2d (data,1.0E-15);
	  }

	  public double[] quadraticMax2d (double [][][] data, double thresholdQuad){
		  return quadraticMax2d (data,thresholdQuad, debugLevel);
	  }

	  public double[] quadraticMax2d (double [][][] data,double thresholdQuad, int debugLevel){
		  double [][] coeff=quadraticApproximation(data, false, 1.0E-20,  thresholdQuad,debugLevel);
		  if (coeff==null) return null;
		  if (coeff[0].length<6) return null;
		  double [][] aM={
				  {2*coeff[0][0],  coeff[0][2]},  // | 2A,  C |
				  {  coeff[0][2],2*coeff[0][1]}}; // |  C, 2B |
		   Matrix M=(new Matrix(aM));
 	   	   double nmQ=normMatix(aM);
	   	   if (debugLevel>3) System.out.println("M.det()="+M.det()+" normMatix(aM)="+nmQ+" data.length="+data.length);
		   if ((nmQ==0.0) || (Math.abs(M.det())/normMatix(aM)<thresholdQuad)) {
			   if (debugLevel>3) System.out.println("quadraticMax2d() failed: M.det()="+M.det()+" normMatix(aM)="+normMatix(aM));
			  return  null;
		   }
		  double [][] aB={
				  {-coeff[0][3]},  // | - D |
				  {-coeff[0][4]}}; // | - E |
		  double [] xy=M.solve(new Matrix(aB)).getColumnPackedCopy();
		  return xy;
	  }
	
	  
	/* ======================================================================== */
	/**
	 * Approximate function z(x,y) as a second degree polynomial (or just linear)
	 * f(x,y)=A*x^2+B*y^2+C*x*y+D*x+E*y+F or f(x,y)=D*x+E*y+F 
	 * data array consists of lines of either 2 or 3 vectors:
	 *  2-element vector x,y
	 *  variable length vector z (should be the same for all samples)
	 *  optional 1- element vector w (weight of the sample)
	 * 
	 * returns array of vectors or null
	 * each vector (one per each z component) is either 6-element-  (A,B,C,D,E,F) if quadratic is possible and enabled
	 * or 3-element - (D,E,F) if linear is possible and quadratic is not possible or disabled
	 * returns null if not enough data even for the linear approximation
	 
	 */
	   public double [][] quadraticApproximation(
			   double [][][] data,
			   boolean forceLinear)  // use linear approximation
			   {
		   return  quadraticApproximation(
				   data,
				   forceLinear,  // use linear approximation
				   1.0E-10,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail) 11.0E-10 failed where it shouldn't?
				   1.0E-15,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				   this.debugLevel);
			   }
	   public double [][] quadraticApproximation(
			   double [][][] data,
			   boolean forceLinear,  // use linear approximation
			   int debugLevel
			   ){
		   return  quadraticApproximation(
				   data,
				   forceLinear,  // use linear approximation
				   1.0E-10,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail) 11.0E-10 failed where it shouldn't?
				   1.0E-15,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				   debugLevel);
/*				   
				   1.0E-12,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail) 11.0E-10 failed where it shouldn't?
				   1.0E-20);  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
*/				   
	   }
	   public double [][] quadraticApproximation(
			   double [][][] data,
			   boolean forceLinear,  // use linear approximation
			   double thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
			   double thresholdQuad)  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
			   {
		   return  quadraticApproximation(
				   data,
				   forceLinear,  // use linear approximation
				   1.0E-10,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail) 11.0E-10 failed where it shouldn't?
				   1.0E-15,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				   this.debugLevel);
	   }
	   
	   public double [][] quadraticApproximation(
				   double [][][] data,
				   boolean forceLinear,  // use linear approximation
				   double thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
				   double thresholdQuad,  // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				   int debugLevel
				   ){
			   if (debugLevel>3) System.out.println("quadraticApproximation(...), debugLevel="+debugLevel+":");
	/* ix, iy - the location of the point with maximal value. We'll approximate the vicinity of that maximum using a
	 * second degree polynomial:
	   Z(x,y)~=A*x^2+B*y^2+C*x*y+D*x+E*y+F
	   by minimizing sum of squared differenceS00between the actual (Z(x,uy)) and approximated values. 
	   and then find the maximum on the approximated surface. Here iS00the math:
	  	   		
	Z(x,y)~=A*x^2+B*y^2+C*x*y+D*x+E*y+F
	minimizing squared error, using W(x,y) aS00weight function

	error=Sum(W(x,y)*((A*x^2+B*y^2+C*x*y+D*x+E*y+F)-Z(x,y))^2)

	error=Sum(W(x,y)*(A^2*x^4 + 2*A*x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)) +(...) )
	0=derror/dA=Sum(W(x,y)*(2*A*x^4 + 2*x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)))
	0=Sum(W(x,y)*(A*x^4 + x^2*(B*y^2+C*x*y+D*x+E*y+F-Z(x,y)))

	S40=Sum(W(x,y)*x^4), etc

	(1) 0=A*S40 + B*S22 + C*S31 +D*S30 +E*S21 +F*S20 - SZ20

	derror/dB:

	error=Sum(W(x,y)*(B^2*y^4 + 2*B*y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)) +(...) )
	0=derror/dB=Sum(W(x,y)*(2*B*y^4 + 2*y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)))
	0=Sum(W(x,y)*(B*y^4 + y^2*(A*x^2+C*x*y+D*x+E*y+F-Z(x,y)))

	(2) 0=B*S04 + A*S22 + C*S13 +D*S12 +E*S03 +F*SY2 - SZ02
	(2) 0=A*S22 + B*S04 + C*S13 +D*S12 +E*S03 +F*SY2 - SZ02

	derror/dC:

	error=Sum(W(x,y)*(C^2*x^2*y^2 + 2*C*x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) +(...) )
	0=derror/dC=Sum(W(x,y)*(2*C*x^2*y^2 + 2*x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) )
	0=Sum(W(x,y)*(C*x^2*y^2 + x*y*(A*x^2+B*y^2+D*x+E*y+F-Z(x,y)) )

	(3) 0= A*S31 +  B*S13 +  C*S22 + D*S21 + E*S12 + F*S11 - SZ11

	derror/dD:

	error=Sum(W(x,y)*(D^2*x^2 + 2*D*x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) +(...) )
	0=derror/dD=Sum(W(x,y)*(2*D*x^2 + 2*x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) )
	0=Sum(W(x,y)*(D*x^2 + x*(A*x^2+B*y^2+C*x*y+E*y+F-Z(x,y)) )

	(4) 0= A*S30 +   B*S12 +  C*S21 + D*S20  + E*S11 +  F*S10  - SZ10

	derror/dE:

	error=Sum(W(x,y)*(E^2*y^2 + 2*E*y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) +(...) )
	0=derror/dE=Sum(W(x,y)*(2*E*y^2 + 2*y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) )
	0=Sum(W(x,y)*(E*y^2 + y*(A*x^2+B*y^2+C*x*y+D*x+F-Z(x,y)) )
	(5) 0= A*S21 +  B*S03 +   C*S12 + D*S11 +  E*SY2  + F*SY  - SZ01

	derror/dF:

	error=Sum(W(x,y)*(F^2 +  2*F*(A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) +(...) )
	0=derror/dF=Sum(W(x,y)*(2*F +  2*(A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) )
	0=Sum(W(x,y)*(F +  (A*x^2+B*y^2+C*x*y+D*x+E*y-Z(x,y)) )
	(6) 0= A*S20 +   B*SY2 +   C*S11 +  D*S10 +   E*SY   + F*S00  - SZ00


	(1) 0= A*S40 + B*S22 + C*S31 + D*S30 + E*S21 + F*S20 - SZ20
	(2) 0= A*S22 + B*S04 + C*S13 + D*S12 + E*S03 + F*S02 - SZ02
	(3) 0= A*S31 + B*S13 + C*S22 + D*S21 + E*S12 + F*S11 - SZ11
	(4) 0= A*S30 + B*S12 + C*S21 + D*S20 + E*S11 + F*S10 - SZ10
	(5) 0= A*S21 + B*S03 + C*S12 + D*S11 + E*S02 + F*S01 - SZ01
	(6) 0= A*S20 + B*S02 + C*S11 + D*S10 + E*S01 + F*S00 - SZ00
	*/
			   int zDim=data[0][1].length;   

			   double w,z,x,x2,x3,x4,y,y2,y3,y4,wz;
			   int i,j,n=0;
			   double S00=0.0,
			   S10=0.0,S01=0.0,
			   S20=0.0,S11=0.0,S02=0.0,
			   S30=0.0,S21=0.0,S12=0.0,S03=0.0,
			   S40=0.0,S31=0.0,S22=0.0,S13=0.0,S04=0.0;
			   double [] SZ00=new double [zDim];
			   double [] SZ01=new double [zDim];
			   double [] SZ10=new double [zDim];
			   double [] SZ11=new double [zDim];
			   double [] SZ02=new double [zDim];
			   double [] SZ20=new double [zDim];
			   for (i=0;i<zDim;i++) {
				   SZ00[i]=0.0;
				   SZ01[i]=0.0;
				   SZ10[i]=0.0;
				   SZ11[i]=0.0;
				   SZ02[i]=0.0;
				   SZ20[i]=0.0;
			   }
			   for (i=0;i<data.length;i++)  {
				   w=(data[i].length>2)? data[i][2][0]:1.0;
				   if (w>0) {
					   n++;
					   x=data[i][0][0];
					   y=data[i][0][1];
					   x2=x*x;
					   y2=y*y;
					   S00+=w;
					   S10+=w*x;
					   S01+=w*y;
					   S11+=w*x*y;
					   S20+=w*x2;
					   S02+=w*y2;
					   if (!forceLinear) {
						   x3=x2*x;
						   x4=x3*x;
						   y3=y2*y;
						   y4=y3*y;
						   S30+=w*x3;
						   S21+=w*x2*y;
						   S12+=w*x*y2;
						   S03+=w*y3;
						   S40+=w*x4;
						   S31+=w*x3*y;
						   S22+=w*x2*y2;
						   S13+=w*x*y3;
						   S04+=w*y4;
					   }
					   for (j=0;j<zDim;j++) {
						   z=data[i][1][j];
						   wz=w*z;
						   SZ00[j]+=wz;
						   SZ10[j]+=wz*x;
						   SZ01[j]+=wz*y;
						   if (!forceLinear) {
							   SZ20[j]+=wz*x2;
							   SZ11[j]+=wz*x*y;
							   SZ02[j]+=wz*y2;
						   }
					   }
					   
				   }
			   }
			   //need to decide if there is enough data for linear and quadratic
			   double [][] mAarrayL= {
					   {S20,S11,S10},
					   {S11,S02,S01},
					   {S10,S01,S00}};
			   Matrix M=new Matrix (mAarrayL);
			   Matrix Z;
	 	   	   if (debugLevel>3) System.out.println(">>> n="+n+" det_lin="+M.det()+" norm_lin="+normMatix(mAarrayL));
	 	   	   double nmL=normMatix(mAarrayL);
			   if ((nmL==0.0) || (Math.abs(M.det())/nmL<thresholdLin)){
// return average value for each channel
				   if (S00==0.0) return null; // not even average
				   double [][] ABCDEF=new double[zDim][3];
				   for (i=0;i<zDim;i++) {
					   ABCDEF[i][0]=0.0;
					   ABCDEF[i][1]=0.0;
					   ABCDEF[i][2]=SZ00[i]/S00;
				   }
				   return ABCDEF;
			   }
			   double []zAarrayL=new double [3];
			   double [][] ABCDEF=new double[zDim][];
//			   double [] zAarrayL={SZ10,SZ01,SZ00};
			   for (i=0;i<zDim;i++) {
				   zAarrayL[0]=SZ10[i];
				   zAarrayL[1]=SZ01[i];
				   zAarrayL[2]=SZ00[i];
			       Z=new Matrix (zAarrayL,3);
			       ABCDEF[i]= M.solve(Z).getRowPackedCopy();
			   }
			   if (forceLinear) return ABCDEF;
			   // quote try quadratic approximation            
			   double [][] mAarrayQ= {
					   {S40,S22,S31,S30,S21,S20},
					   {S22,S04,S13,S12,S03,S02},
					   {S31,S13,S22,S21,S12,S11},
					   {S30,S12,S21,S20,S11,S10},
					   {S21,S03,S12,S11,S02,S01},
					   {S20,S02,S11,S10,S01,S00}};
			   M=new Matrix (mAarrayQ);
	 	   	   if (debugLevel>3) {
	 	   		   System.out.println("    n="+n+" det_quad="+M.det()+" norm_quad="+normMatix(mAarrayQ)+" data.length="+data.length);
	 	   		   M.print(10,5);
	 	   	   }
	 	   	   double nmQ=normMatix(mAarrayQ);
			   if ((nmQ==0.0) || (Math.abs(M.det())/normMatix(mAarrayQ)<thresholdQuad)) {
				   if (debugLevel>0) System.out.println("Using linear approximation, M.det()="+M.det()+
						   " normMatix(mAarrayQ)="+normMatix(mAarrayQ)+
						   ", thresholdQuad="+thresholdQuad+
						   ", nmQ="+nmQ+
						   ", Math.abs(M.det())/normMatix(mAarrayQ)="+(Math.abs(M.det())/normMatix(mAarrayQ))); //did not happen
				   return ABCDEF; // not enough data for the quadratic approximation, return linear
			   }
//			   double [] zAarrayQ={SZ20,SZ02,SZ11,SZ10,SZ01,SZ00};
			   double [] zAarrayQ=new double [6];
			   for (i=0;i<zDim;i++) {
				   zAarrayQ[0]=SZ20[i];
				   zAarrayQ[1]=SZ02[i];
				   zAarrayQ[2]=SZ11[i];
				   zAarrayQ[3]=SZ10[i];
				   zAarrayQ[4]=SZ01[i];
				   zAarrayQ[5]=SZ00[i];
				   Z=new Matrix (zAarrayQ,6);
				   ABCDEF[i]= M.solve(Z).getRowPackedCopy();
			   }
			   return ABCDEF;
		   }
//			calcualte "volume" made of the matrix row-vectors, placed orthogonally
		// to be compared to determinant	   
			public double normMatix(double [][] a) {
		        double d,norm=1.0;
		        for (int i=0;i<a.length;i++) {
		        	d=0;
		        	for (int j=0;j<a[i].length;j++) d+=a[i][j]*a[i][j];
		        	norm*=Math.sqrt(d);
		        }
				return norm;
			}


//RuntimeException
}
