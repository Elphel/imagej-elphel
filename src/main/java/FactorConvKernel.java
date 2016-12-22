/**
 **
 ** FactorConvKernel Split convolution kernel into small asymmetrical 
 ** (to be applied directly to Bayer data) and large symmetrical to use
 ** with MDCT
 **
 ** Copyright (C) 2016 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  FactorConvKernel.java is free software: you can redistribute it and/or modify
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
import Jama.LUDecomposition;
import Jama.Matrix;
import ij.IJ;

public class FactorConvKernel {
	public int asym_size =           6;
	public int sym_radius =          8; // 2*2^n  - for DCT
	public double[]  target_kernel = null; // should be expanded to 2*(sym_radius)+asym_size- 1 in each direction 
//	public double[]  asym_kernel =   null;
//	public double[]  sym_kernel =    null;
	public int       debugLevel =    3;
	public double    init_lambda =   0.001;
	
	public double    compactness_weight = 1.0; // realtive "importance of asymmetrical kernel being compact"
	public double    lambdaStepUp=   8.0; // multiply lambda by this if result is worse
	public double    lambdaStepDown= 0.5; // multiply lambda by this if result is better
	public double    thresholdFinish=0.001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
	public int       numIterations=  100; // maximal number of iterations
	public double    maxLambda=      100.0;  // max lambda to fail
	public boolean   stopOnFailure = true;
	
	
	public double    sym_kernel_scale =1.0;//  to scale sym kernel to match original one
	public int       center_i0;
	public int       center_j0;
	
	public double    lambda =     init_lambda;
	public double    currentRMS ;
	public double    firstRMS;
	public double    nextRMS ;

	public double    currentRMSPure=-1.0; // calculated RMS for the currentVector->currentfX
    public double    nextRMSPure=   -1.0; // calculated RMS for the nextVector->nextfX
    public double    firstRMSPure=  -1.0; // RMS before current series of LMA started
	
	
	public double [] currentfX = null; // conv
	public double [] nextfX =  null;
	public double [] currentVector=null; //kern_vector;
	public double [] nextVector;
    public double [][] jacobian=null; // partial derivatives of fX (above) by parameters to be adjusted (rows)
	public int       iterationStepNumber=0;
	public long      startTime=0;
	public LMAArrays lMAArrays=null;
	public LMAArrays  savedLMAArrays=null;
    public double [] lastImprovements= {-1.0,-1.0}; // {last improvement, previous improvement}. If both >0 and < thresholdFinish - done
	
	
	
	public class LMAArrays {
		public double [][] jTByJ=  null; // jacobian multiplied by Jacobian transposed
		public double []   jTByDiff=null; // jacobian multiplied difference vector
		public LMAArrays clone() {
			LMAArrays lma=new LMAArrays();
			lma.jTByJ = this.jTByJ.clone();
			for (int i=0;i<this.jTByJ.length;i++) lma.jTByJ[i]=this.jTByJ[i].clone();
			lma.jTByDiff=this.jTByDiff.clone();
			return lma;
		}
	}

	public FactorConvKernel(){
	}

	public FactorConvKernel(int asym_size, int sym_radius){
		this.asym_size = asym_size;
		this.sym_radius =  sym_radius;
	}

	public boolean calcKernels(
			double []target_kernel,
			int asym_size,
			int sym_radius){
		this.asym_size = asym_size;
		this.sym_radius =  sym_radius;
		this.target_kernel = target_kernel;
		return levenbergMarquardt();
	}
	
	public double [] getSymKernel(){
		return  getSymKernel(currentVector,sym_kernel_scale);
	}

	public double [] getAsymKernel(){
		return  getAsymKernel(currentVector,1.0/sym_kernel_scale);
	}

	public double [] getConvolved(){ // check that it matches original
		return this.currentfX;
	}
	public void setDebugLevel(int debugLevel){
		this.debugLevel =debugLevel;
	}
	public void setCompactnessWeight(double compactness_weight){
		this.compactness_weight =compactness_weight;
	}
	
	
	private double [] getSymKernel(
			double [] kvect,
			double    scale){
		int sym_len= sym_radius * sym_radius;
		double [] sym_kernel = new double [sym_len];
		sym_kernel[0] = scale;
		for (int i =1; i<  sym_kernel.length; i++) sym_kernel[i] = scale * kvect[i-1];
		if (debugLevel>1){
			System.out.println("getSymKernel(): sym_kernel.length="+sym_kernel.length);
		}
		return sym_kernel;
	}
	private double [] getAsymKernel(
			double [] kvect,
			double    scale){ // should be 1/scale of the sym_vector
		int sym_len= sym_radius * sym_radius;
		double [] asym_kernel = new double [asym_size*asym_size];
		int indx = sym_len -1;
		for (int i =0; i<  asym_kernel.length; i++) asym_kernel[i] = scale * kvect[indx++];
		if (debugLevel>1){
			System.out.println("getAsymKernel(): asym_kernel.length="+asym_kernel.length);
		}
		return asym_kernel;
	}
	
	// initial estimation
	private double [] setInitialVector(
			double [] target_kernel) // should be (asym_size + 2*sym_radius-1)**2
	{
		int conv_size = asym_size + 2*sym_radius-1;
//		int sym_size = 2*sym_radius - 1;
		int sym_rad_m1 = sym_radius - 1; // 7
		// find center of the target kernel squared value
		double s0=0.0,sx=0.0,sy=0.0;
		for (int i = 0; i < conv_size; i++){
			for (int j = 0; j < conv_size; j++){
				double d = target_kernel[conv_size*i+j];
				d  *= d;
				s0 += d;
				sx += d*j;
				sy += d*i;
			}
		}
//		int async_center = asym_size / 2; // will be 2 for 5x5 and 3 for 6x6 - more room in negative
		
		
		int j0= (int) Math.round(sx/s0 - sym_rad_m1); // should be ~ async_center 
		int i0= (int) Math.round(sy/s0 - sym_rad_m1); // should be ~ async_center
		if (debugLevel>1){
			System.out.println("setInitialVector(): fi0 = "+(sy/s0 - sym_radius)+" fj0 = "+(sx/s0 - sym_radius));
			System.out.println("setInitialVector(): i0 = "+i0+" j0 = "+j0 + " asym_size="+asym_size);
		}
		// fit i0,j0 to asym_kernel (it should be larger)
//		i0 += asym_size/2; 
//		j0 += asym_size/2;
		if ((i0<0) || (i0>=asym_size) || (i0<0) || (i0>=asym_size)){
			System.out.println("Warning: kernel center too far : i="+(i0 - asym_size/2)+", j= "+(j0 - asym_size/2));
			if (i0 <  0) i0=0;
			if (j0 <  0) j0=0;
			if (i0 >= asym_size) i0=asym_size-1;
			if (j0 >= asym_size) j0=asym_size-1;
		}
		double [] asym_kernel = new double [asym_size * asym_size];
		for (int i = 0; i < asym_kernel.length; i++) asym_kernel[i] = 0.0;
		asym_kernel[asym_size*i0+j0] = 1.0;
		
		center_i0 = i0; // center to calcualte compatess odf asymmetrical kernel
		center_j0 = j0; // center to calcualte compatess odf asymmetrical kernel

		if (debugLevel>2){
			for (int i = 0; i < asym_kernel.length; i++) {
				System.out.println("asym_kernel["+i+"] = "+asym_kernel[i]);
			}
		}
		
		
		//		i0 -= asym_size/2;  
//		j0 -= asym_size/2;
		double [] sym_kernel = new double [sym_radius * sym_radius];
		int []    sym_kernel_count = new int [sym_radius * sym_radius];
		for (int i = 0; i < sym_kernel.length; i++){
			sym_kernel[i] = 0.0;
			sym_kernel_count[i] = 0;
		}
		for (int i=0; i< (2*sym_radius -1); i++ ){
			for (int j=0; j< (2*sym_radius -1); j++ ){
				int indx =((i >= sym_rad_m1)? (i-sym_rad_m1): (sym_rad_m1 -i)) * sym_radius +
						  ((j >= sym_rad_m1)? (j-sym_rad_m1): (sym_rad_m1 -j));
				sym_kernel[indx] += target_kernel[conv_size*(i+i0) + (j+j0)];
				sym_kernel_count[indx]++;
			}
			
		}
		for (int i = 0; i < sym_kernel.length; i++){
			if (sym_kernel_count[i] >0) sym_kernel[i] /= sym_kernel_count[i];
			else sym_kernel[i] = 0.0;
		}		
		return setVectorFromKernels(sym_kernel, asym_kernel);
	}
	
	private double [] setVectorFromKernels(
			double [] sym_kernel,
			double [] asym_kernel){
		if (this.debugLevel>2){
			System.out.println("setVectorFromKernels(): sym_kernel.length=    "+sym_kernel.length);
			System.out.println("setVectorFromKernels(): asym_kernel.length=   "+asym_kernel.length);
		}
		sym_kernel_scale = sym_kernel[0];
		double [] kvect = new double [sym_kernel.length+asym_kernel.length -1];
		int indx = 0;
		for (int i =1; i<  sym_kernel.length; i++) kvect[indx++] = sym_kernel[i]/sym_kernel[0];
		for (int i =0; i< asym_kernel.length; i++) kvect[indx++] = asym_kernel[i]*sym_kernel[0];
		if (debugLevel>2){
			for (int i = 0; i < kvect.length; i++) {
				System.out.println("kvect["+i+"] = "+kvect[i]);
			}
		}

		return kvect;
	}

	private double [] convKernels(
		double [] kvect) // first - all elements of sym kernel but [0] (replaced by 1.0), then - asym ones
		{
		int conv_size = asym_size + 2*sym_radius-1;
		double [] conv_data = new double [conv_size*conv_size];
//		int asym_start= (2 * sym_radius - 1) * (2 * sym_radius - 1)-1;
		int asym_start=  sym_radius * sym_radius - 1;
		if (this.debugLevel>2){
			System.out.println("convKernels(): vector_length=    "+kvect.length);
			System.out.println("convKernels(): sym_radius=       "+sym_radius);
			System.out.println("convKernels(): asym_size=        "+asym_size);
			System.out.println("convKernels(): conv_size=        "+conv_size);
			System.out.println("convKernels(): conv_data.length= "+conv_data.length);
			System.out.println("convKernels(): asym_start=       "+asym_start);
		}
		for (int i = 0; i < conv_data.length; i++) conv_data[i] = 0.0;
		for (int i = -sym_radius+1; i < sym_radius; i++) {
			for (int j = -sym_radius+1; j < sym_radius; j++) {
				int indx = ((i < 0)? -i : i) * sym_radius + ((j < 0)? -j : j);
				double sd = (indx >0)? kvect[indx -1] : 1.0;
				int base_indx = conv_size * (i + sym_radius -1) + (j + sym_radius -1);  
				for (int ia=0; ia < asym_size; ia++) {
					for (int ja=0; ja < asym_size; ja++) {
						conv_data[base_indx + conv_size * ia + ja] += sd* kvect[asym_size * ia + ja + asym_start];
					}
				}
			}
		}
		return conv_data;
	}

	private double getRms(
			double [] conv_data,
			double [] target_kernel){
		double rms = 0;
		for (int i=0; i < conv_data.length; i++){
			double d = target_kernel[i] - conv_data[i]; 
			rms += d * d;
		}
		return Math.sqrt(rms/conv_data.length);
	}

	private double getCompactRms(double [] kvect){ // uses global compactness_weight, sym_kernel_scale, center_i0, center_j0 and asym_kernel
		double cw = getCompactWeight();
		double rms = 0;
		int indx = sym_radius * sym_radius -1;
		for (int ia=0; ia <asym_size;ia++){
			for (int ja=0;ja< asym_size;ja++){
				int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
				double d = kvect[indx++]*ir2; 
			rms += d * d;
			}
		}
		return cw*Math.sqrt(rms)/(asym_size*asym_size); // square for area, another - to normalize by r^2
	}
	
/*
		int sym_len= sym_radius * sym_radius;
		double [] asym_kernel = new double [asym_size*asym_size];
		int indx = sym_radius * sym_radius -1;
		for (int i =0; i<  asym_kernel.length; i++) asym_kernel[i] = scale * kvect[indx++];
		if (debugLevel>1){
			System.out.println("getAsymKernel(): asym_kernel.length="+asym_kernel.length);
		}
	
 */
	private double getCompactWeight(){
//		return compactness_weight*sym_kernel_scale/(asym_size*asym_size); // use
		return compactness_weight*sym_kernel_scale/(asym_size*asym_size*asym_size*asym_size); // use
	}

	private double [][] getJacobian(
			double [] kvect)
	{
		int conv_size =       asym_size + 2*sym_radius-1;
//		int asym_start =     (2 * sym_radius - 1) * (2 * sym_radius - 1)-1;
		int asym_start=  sym_radius * sym_radius - 1;
		double cw = getCompactWeight();
		double [][] jacob = new double [kvect.length][conv_size*conv_size + asym_size*asym_size];
		int asym_terms_start=conv_size*conv_size;
		for (int i = -sym_radius+1; i < sym_radius; i++) {
			for (int j = -sym_radius+1; j < sym_radius; j++) {
				int indx = ((i < 0)? -i : i) * sym_radius + ((j < 0)? -j : j);
				double sd = (indx >0)? kvect[indx -1] : 1.0;
				int base_indx = conv_size * (i + sym_radius -1) + (j + sym_radius -1);  
				for (int ia=0; ia < asym_size; ia++) {
					for (int ja=0; ja < asym_size; ja++) {
						int async_index = asym_size*ia + ja+ asym_start;
						int conv_index = base_indx + conv_size * ia + ja;
						if (indx > 0) jacob[indx-1][conv_index] += kvect[async_index];
						jacob[async_index][conv_index] += sd;
//						conv_data[base_indx + conv_size * ia + ja] += sd* kvect[asym_size * ia + ja + asym_start];

					}
				}
			}
		}
		for (int i = 0; i < jacob.length; i++) for (int j=asym_terms_start; j < jacob[i].length;j++) jacob[i][j] = 0.0;
		for (int ia=0; ia <asym_size;ia++){
			for (int ja=0;ja< asym_size;ja++){
				int async_index = asym_size*ia + ja;
				int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
				jacob[async_index + asym_start][async_index+asym_terms_start] = ir2 * cw;
			}
		}
		
		return jacob;
	}
	
	private double [][] getJTByJ(double [][] jacob){
		double [][] jTByJ = new double [jacob.length][jacob.length];
		for (int i = 0; i < jacob.length; i++ ){
			for (int j = 0; j < jacob.length; j++ ){
				if (j<i){
					jTByJ[i][j] = jTByJ[j][i];
				} else {
					jTByJ[i][j] = 0;
					for (int k=0; k< jacob[i].length; k++){
						jTByJ[i][j] += jacob[i][k] * jacob[j][k];
					}
				}
			}
		}
		return jTByJ;
	}

	private double [] getJTByDiff(
			double [][] jacob,       // jacobian
			double [] target_kernel, // target kernel
			double [] conv_data,     // current convolution result of async_kernel (*) sync_kernel
			double [] kvect          // parameter vector - used asym values
			){
		double [] jTByDiff = new double [jacob.length];
		for (int i=0; i < jTByDiff.length; i++){
			jTByDiff[i] = 0;
//			for (int k=0; k< jacob[i].length; k++){
			for (int k=0; k< target_kernel.length; k++){
				jTByDiff[i] += jacob[i][k]*(target_kernel[k]-conv_data[k]);
			}
		}
		int asym_start=       sym_radius * sym_radius - 1;
		int conv_size =       asym_size + 2*sym_radius-1;
		int asym_terms_start= conv_size*conv_size;
		double cw = getCompactWeight();
		for (int ia=0; ia <asym_size;ia++){
			for (int ja=0;ja< asym_size;ja++){
				int async_index = asym_size*ia + ja;
				int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
//				jacob[async_index + asym_start][async_index+asym_terms_start] = ir2 * cw;
				jTByDiff[async_index + asym_start] += jacob[async_index + asym_start][async_index+asym_terms_start]*
						ir2 * cw* kvect[asym_start + async_index];
			}
		}
		
		return jTByDiff;
	}

	private double [] solveLMA(
		    LMAArrays lMAArrays,
		    double lambda,
			int debugLevel){
		this.debugLevel=debugLevel;
		double [][] JtByJmod= lMAArrays.jTByJ.clone();
		int numPars=JtByJmod.length;
		for (int i=0;i<numPars;i++){
			JtByJmod[i]=lMAArrays.jTByJ[i].clone();
			JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
		}
		//     M*Ma=Mb
		Matrix M=new Matrix(JtByJmod);
		if (debugLevel>2) {
			System.out.println("Jt*J -lambda* diag(Jt*J), lambda="+lambda+":");
			M.print(10, 5);
		}
		Matrix Mb=new Matrix(lMAArrays.jTByDiff,numPars); // single column
		if (!(new LUDecomposition(M)).isNonsingular()){
			double [][] arr=M.getArray();
			System.out.println("Singular Matrix "+arr.length+"x"+arr[0].length);
			// any rowsx off all 0.0?
			for (int n=0;n<arr.length;n++){
				boolean zeroRow=true;
				for (int i=0;i<arr[n].length;i++) if (arr[n][i]!=0.0){
					zeroRow=false;
					break;
				}
				if (zeroRow){
					System.out.println("Row of all zeros: "+n);
				}
			}
			//            M.print(10, 5);
			return null;
		}
		Matrix Ma=M.solve(Mb); // singular
		return Ma.getColumnPackedCopy();
	}

	private boolean levenbergMarquardt(){
    	this.startTime=System.nanoTime();
    	this.firstRMS=-1; //undefined
    	this.currentVector = setInitialVector(target_kernel); // should be (asym_size + 2*sym_radius-1)**2
    	if (this.numIterations < 0){
    		this.currentfX=convKernels  (this.currentVector); // to be able to read back
            return true;    		
    	}
    	
    	
    	while (true) { // loop for the same series
    		boolean [] state=stepLevenbergMarquardtFirst();
    		if (this.debugLevel>1) System.out.println(this.iterationStepNumber+": stepLevenbergMarquardtFirst()==>"+state[1]+":"+state[0]);
    		//    		boolean cont=true;
    		// Make it success if this.currentRMS<this.firstRMS even if LMA failed to converge
    		if (state[1] && !state[0] && (this.firstRMS>this.currentRMS)){
    			if (this.debugLevel>1) System.out.println("LMA failed to converge, but RMS improved from the initial value ("+this.currentRMS+" < "+this.firstRMS+")");
    			state[0]=true;
    		}
    		if ((this.stopOnFailure && state[1] && !state[0])){
    			if (this.debugLevel>0){
    				System.out.println("LevenbergMarquardt(): step ="+this.iterationStepNumber+
    						", RMS="+IJ.d2s(this.currentRMS,8)+
    						" ("+IJ.d2s(this.firstRMS,8)+") "+
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    			}
    			long startDialogTime=System.nanoTime();
    			this.startTime+=(System.nanoTime()-startDialogTime); // do not count time used by the User.
    			//    			if (this.showThisImages) showDiff (this.currentfX, "fit-"+this.iterationStepNumber);
    			//    			if (this.showNextImages) showDiff (this.nextfX,    "fit-"+(this.iterationStepNumber+1));
    		} else if (this.debugLevel>1){
    			System.out.println("==> LevenbergMarquardt(): before action step ="+this.iterationStepNumber+
    					", RMS="+IJ.d2s(this.currentRMS,8)+
    					" ("+IJ.d2s(this.firstRMS,8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    		}
    		stepLevenbergMarquardtAction(); // apply step - in any case?
    		if ((this.debugLevel>0) && ((this.debugLevel>1) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
    			System.out.println("--> LevenbergMarquardt(): series:step = "+this.iterationStepNumber+
    					", RMS="+IJ.d2s(this.currentRMS,8)+
    					" ("+IJ.d2s(this.firstRMS,8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    		}
    		//stepLevenbergMarquardtAction();    			
    		if (state[1]) {
    			if (!state[0]) return false; // sequence failed
    			break; // while (true), proceed to the next series
    		}
    	}
    	//    		if (this.fittingStrategy.isLastSeries(this.seriesNumber)) break;
    	if (this.debugLevel>0) System.out.println("LevenbergMarquardt(): RMS="+this.currentRMS+
    			" ("+this.firstRMS+") "+
    			") at "+ IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3));
    	return true; // all series done
    }
	
	private boolean [] stepLevenbergMarquardtFirst(){
		double [] deltas=null;
		if (this.currentVector==null) {
			this.currentRMS=-1;
			this.currentRMSPure=-1;
			this.currentfX=null; // invalidate
			this.jacobian=null;  // invalidate
			this.lMAArrays=null;
			lastImprovements[0]=-1.0;
			lastImprovements[1]=-1.0;
		}
		// calculate  this.currentfX, this.jacobian if needed
		if (this.debugLevel>2) {
			System.out.println("this.currentVector");
			for (int i=0;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		//    	if ((this.currentfX==null)|| ((this.jacobian==null) && !this.threadedLMA )) {
		if ((this.currentfX==null)|| (this.lMAArrays==null)) {
			this.currentfX=convKernels  (this.currentVector); // first - all elements of sym kernel but [0] (replaced by 1.0), then - asym ones
			this.jacobian = getJacobian (this.currentVector);
			this.lMAArrays= new LMAArrays();
			lMAArrays.jTByJ=getJTByJ(
					this.jacobian);
			lMAArrays.jTByDiff=getJTByDiff(
					this.jacobian,
					this.target_kernel,   // target kernel to factor
					this.currentfX,      // current convolution result of async_kernel (*) sync_kernel
					this.currentVector); // used to control compactness of asym_kernel


			this.currentRMSPure= getRms(this.currentfX,this.target_kernel);
			this.currentRMS=     getCompactRms(this.currentVector) + this.currentRMSPure;
			
			if (this.debugLevel>1) {
				System.out.println("initial RMS="+IJ.d2s(this.currentRMS,8)+
						" ("+IJ.d2s(this.currentRMSPure,8)+")"+
						". Calculating next Jacobian. Points:"+this.target_kernel.length+" Parameters:"+this.currentVector.length);
			}
		} else {
			this.currentRMSPure= getRms(this.currentfX,this.target_kernel); 
			this.currentRMS=     getCompactRms(this.currentVector) + this.currentRMSPure;
		}
		//		this.currentRMS= calcError(calcYminusFx(this.currentfX));
		if (this.firstRMS<0) {
			this.firstRMS=this.currentRMS;
    		this.firstRMSPure=this.currentRMSPure;
		}
		// calculate deltas

		deltas=solveLMA(this.lMAArrays,	this.lambda, this.debugLevel);

		boolean matrixNonSingular=true;
		if (deltas==null) {
			deltas=new double[this.currentVector.length];
			for (int i=0;i<deltas.length;i++) deltas[i]=0.0;
			matrixNonSingular=false;
		}
		if (this.debugLevel>2) {
			System.out.println("deltas");
			for (int i=0;i<deltas.length;i++){
				System.out.println(i+": "+ deltas[i]);
			}
		}
		// apply deltas    	
		this.nextVector=this.currentVector.clone();
		for (int i=0;i<this.nextVector.length;i++) this.nextVector[i]+=deltas[i];
		// another option - do not calculate J now, just fX. and late - calculate both if it was improvement    	
		//    	save current Jacobian
		if (this.debugLevel>2) {
			System.out.println("this.nextVector");
			for (int i=0;i<this.nextVector.length;i++){
				System.out.println(i+": "+ this.nextVector[i]);
			}
		}

		//        this.savedJacobian=this.jacobian;
		this.savedLMAArrays=lMAArrays.clone();
		this.jacobian=null; // not needed, just to catch bugs
		// calculate next vector and Jacobian  (this.jacobian)  	
		//    	this.nextfX=calculateFxAndJacobian(this.nextVector,true); //=========== OLD

		//	    	this.nextfX=calculateFxAndJacobian(this.nextVector,true);
		//			this.lMAArrays=calculateJacobianArrays(this.nextfX);
		this.nextfX=convKernels(this.nextVector);
		this.jacobian = getJacobian(this.currentVector);

		this.lMAArrays= new LMAArrays();
		lMAArrays.jTByJ=getJTByJ(
				this.jacobian);
		lMAArrays.jTByDiff=getJTByDiff(
				this.jacobian,
				this.target_kernel,    // target kernel to factor
				this.nextfX,           // next convolution result of async_kernel (*) sync_kernel
				this.nextVector);      // used to control compactness of asym_kernel
				

		this.nextRMSPure= getRms(this.nextfX,this.target_kernel); 
		this.nextRMS=     getCompactRms(this.nextVector) + this.nextRMSPure;
		
		
		
		

		this.lastImprovements[1]=this.lastImprovements[0];
		this.lastImprovements[0]=this.currentRMS-this.nextRMS;
		if (this.debugLevel>2) {
			System.out.println("stepLMA this.currentRMS="+this.currentRMS+
					", this.nextRMS="+this.nextRMS+
					", delta="+(this.currentRMS-this.nextRMS));
		}
		boolean [] status={matrixNonSingular && (this.nextRMS<=this.currentRMS),!matrixNonSingular};
		// additional test if "worse" but the difference is too small, it was be caused by computation error, like here:
		//stepLevenbergMarquardtAction() step=27, this.currentRMS=0.17068403807026408,   this.nextRMS=0.1706840380702647

		if (!status[0] && matrixNonSingular) {
			if (this.nextRMS<(this.currentRMS+this.currentRMS*this.thresholdFinish*0.01)) {
				this.nextRMS=this.currentRMS;
				status[0]=true;
				status[1]=true;
				this.lastImprovements[0]=0.0;
				if (this.debugLevel>1) {
					System.out.println("New RMS error is larger than the old one, but the difference is too small to be trusted ");
					System.out.println(
							"stepLMA this.currentRMS="+this.currentRMS+
							", this.currentRMSPure="+this.currentRMSPure+
							", this.nextRMS="+this.nextRMS+
							", this.nextRMSPure="+this.nextRMSPure+
							", delta="+(this.currentRMS-this.nextRMS)+
							", deltaPure="+(this.currentRMSPure-this.nextRMSPure));
				}
			}
		}
		if (status[0] && matrixNonSingular) { //improved
			status[1]=(this.iterationStepNumber>this.numIterations) || ( // done
					(this.lastImprovements[0]>=0.0) &&
					(this.lastImprovements[0]<this.thresholdFinish*this.currentRMS) &&
					(this.lastImprovements[1]>=0.0) &&
					(this.lastImprovements[1]<this.thresholdFinish*this.currentRMS));
		} else if (matrixNonSingular){
			//    		this.jacobian=this.savedJacobian;// restore saved Jacobian
			this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff

			status[1]=(this.iterationStepNumber>this.numIterations) || // failed
					((this.lambda*this.lambdaStepUp)>this.maxLambda);
		}
		///this.currentRMS    	
		//TODO: add other failures leading to result failure?    	
		if (this.debugLevel>2) {
			System.out.println("stepLevenbergMarquardtFirst()=>"+status[0]+","+status[1]);
		}
		return status;
	}
	
    /**
     * Apply fitting step 
     */
	private void stepLevenbergMarquardtAction(){// 
    	this.iterationStepNumber++;
// apply/revert,modify lambda    	
		if (this.debugLevel>1) {
			System.out.println(
					"stepLevenbergMarquardtAction() step="+this.iterationStepNumber+
					", this.currentRMS="+this.currentRMS+
					", this.currentRMSPure="+this.currentRMSPure+
					", this.nextRMS="+this.nextRMS+
					", this.nextRMSPure="+this.nextRMSPure+
					" lambda="+this.lambda+" at "+IJ.d2s(0.000000001*(System.nanoTime()-this.startTime),3)+" sec");
		}
    	if (this.nextRMS<this.currentRMS) { //improved
    		this.lambda*=this.lambdaStepDown;
    		this.currentRMS=this.nextRMS;
    		this.currentfX=this.nextfX;
    		this.currentVector=this.nextVector;
    	} else {
    		this.lambda*=this.lambdaStepUp;
//    		this.jacobian=this.savedJacobian;// restore saved Jacobian
    		this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff
    	}
    }

}
