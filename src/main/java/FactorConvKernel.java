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
import java.util.ArrayList;
import java.util.Random;

import Jama.LUDecomposition;
import Jama.Matrix;
import ij.IJ;

/*
 * TODO: Add options for rotation-symmetrical (2-fold) kernels (with the same center as the symmetrical one for DCT) and a
 * sub-pixel shift (4 (2x2) pixel kernel. As the target kernel is rather smooth, half-pixel shift widening can be absorbed
 * by corresponding "sharpening" of the symmetrical kernel.
 * 
 * Try to minimize (or limit) the number of non-zero elements of asymmetrical kernels - they are convolved directly with the
 * Bayer mosaic data.
 * 1. Do as now, then select N of the highest absolute value asymmetrical elements, mask out (and zero) all others
 * 2. Optionally try to improve: Remove some from step 1, then add one-by-one: select from neighbors (or neighbors of neighbors?)
 * add the one that gets the best improvement.
 * 3. Try "jumping" (one remove, one add)?
 */


public class FactorConvKernel {
	public int       asym_size =     6;
	public int       sym_radius =    8; // 2*2^n  - for DCT
	public double[]  target_kernel = null; // should be expanded to 2*(sym_radius)+asym_size- 1 in each direction
	public double    target_rms;   // Target kernel rma (to compare with residual error)
	public int       debugLevel =    3;
	public double    init_lambda =   0.001;
	public int       asym_pixels =  10;      // maximal number of non-zero pixels in asymmmetrical kernel
	public int       asym_distance = 2;      // how far to seed a new pixel
	
	public double    compactness_weight = 1.0; // realtive "importance of asymmetrical kernel being compact"
	public int       asym_tax_free = 1; // do not apply compactness_weight for pixels close to the center
	
	
	public double    lambdaStepUp=   8.0; // multiply lambda by this if result is worse
	public double    lambdaStepDown= 0.5; // multiply lambda by this if result is better
	public double    thresholdFinish=0.001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
	public int       numIterations=  100; // maximal number of iterations
	public double    maxLambda=      1000.0;  // max lambda to fail
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
	
	
//    public boolean [] vector_mask = null;
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
    public double    goal_rms_pure;	
	
	
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
			int sym_radius,
			double fact_precision,
			boolean [] mask){
		this.asym_size = asym_size;
		this.sym_radius =  sym_radius;
		this.target_kernel = target_kernel;
		initLevenbergMarquardt(fact_precision);
		if (mask != null){
			if (this.debugLevel>2) {
				System.out.println("calcKernels(): this.currentVector 0");
				for (int i=63;i<this.currentVector.length;i++){
					System.out.println(i+": "+ this.currentVector[i]);
				}
			}
			
			int asym_start=  sym_radius * sym_radius - 1;
			for (int i = 0; i<mask.length; i++) if (!mask[i]) this.currentVector[asym_start+i] = Double.NaN;
			if (debugLevel>0){
				System.out.println("mask.length="+mask.length+" asym_start="+asym_start+" this.currentVector.length="+this.currentVector.length);
				System.out.println("asym mask: ");
				for (int ii=0;ii < asym_size;ii++){
					System.out.print(ii+": ");
					for (int jj=0;jj < asym_size;jj++){
						System.out.print((mask[ii * asym_size + jj]?" X":" .")+" ");
					}
					System.out.println();	
				}
			}
		}
		if (this.debugLevel>2) {
			System.out.println("calcKernels(): this.currentVector 1");
			for (int i=63;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		return levenbergMarquardt();
	}

	public int calcKernels(
			double []target_kernel,
			int asym_size,
			int sym_radius,
			double fact_precision,
			int asym_pixels,          // maximal number of non-zero pixels in asymmmetrical kernel
			int asym_distance){       // how far to seed a new pixel
		this.asym_size =     asym_size;
		this.sym_radius =    sym_radius;
		this.target_kernel = target_kernel;
		this.asym_pixels =   asym_pixels;
		this.asym_distance = asym_distance;
		int asym_start=  sym_radius * sym_radius - 1;
		int num_lma = 0;
		double s = 0.0;
		for (int i = 0; i<target_kernel.length; i++){
			s+= target_kernel[i]*target_kernel[i];
		}
		this.goal_rms_pure = Math.sqrt(s/target_kernel.length)*fact_precision;
		
		double [] bestRms =  new double [asym_pixels];
		int    [] enPixels = new int    [asym_pixels];
//		int       numPixels = 0;
    	this.startTime=System.nanoTime();
    	
    	double[] initialVector = setInitialVector(target_kernel, null); // should be (asym_size + 2*sym_radius-1)**2
    	
		enPixels[0] = center_i0*asym_size+center_j0;
    	boolean [] mask = new boolean [asym_size*asym_size];
    	for (int i = 0; i<mask.length; i++) mask[i] = false;
    	mask[enPixels[0]] = true;
    	this.currentVector = initialVector.clone();
    	System.out.println("mask.length="+mask.length+" asym_start="+asym_start+" this.currentVector.length="+this.currentVector.length);
    	for (int i = 0; i<mask.length; i++) if (!mask[i]) this.currentVector[asym_start+i] = Double.NaN;
    	boolean OK =       levenbergMarquardt();
    	num_lma++;
    	bestRms[0] =       this.currentRMSPure; 
/*    	
    	for (int i = 0; i<numPixels; i++) mask[enPixels[i]] = true;
    	this.currentVector = initialVector.clone();
    	for (int i = 0; i<mask.length; i++) if (!mask[i]) this.currentVector[asym_start+1] = Double.NaN;
    	boolean OK = levenbergMarquardt();
    	if (numPixels == 1){
        	bestRms[numPixels-1] = this.currentRMSPure; 
    	}
*/		
		int numPixels = 0;
    	for (numPixels = 1; numPixels < asym_pixels; numPixels++) {
    		if (debugLevel >0) {
    			System.out.println("calcKernels() numPixels="+numPixels);
    		}
    		
        	for (int i = 0; i < mask.length; i++) mask[i] = false;
        	for (int i = 0; i < numPixels; i++) mask[enPixels[i]] = true;
        	boolean[] mask0 = mask.clone();
        	for (int n = 0; n < asym_distance; n++){
            	boolean[] mask1 = mask0.clone();
            	if (asym_size < 0){
            		asym_size = -asym_size;
            		for (int i = 0; i <asym_size; i++){
            			for (int j = 1; j<asym_size; j++){
            				mask1[asym_size*i + j] |=                mask0[asym_size*i + j-1];
            				mask1[asym_size*i + (asym_size-j-1)] |=  mask0[asym_size*i + (asym_size-j)];
            				mask1[asym_size*j + i] |=                mask0[asym_size*(j-1) + i];
            				mask1[asym_size*(asym_size-j-1) + i] |=  mask0[asym_size*(asym_size-j) + i];
            			}
            		}
            	} else {
            		// hor
            		for (int i = 0; i <asym_size; i++){
            			for (int j = 1; j<asym_size; j++){
            				mask1[asym_size*i + j] |=                mask0[asym_size*i + j-1];
            				mask1[asym_size*i + (asym_size-j-1)] |=  mask0[asym_size*i + (asym_size-j)];
            			}
            		}
            		// vert
                	mask0 = mask1.clone();
            		for (int i = 0; i <asym_size; i++){
            			for (int j = 1; j<asym_size; j++){
            				mask1[asym_size*j + i] |=                mask0[asym_size*(j-1) + i];
            				mask1[asym_size*(asym_size-j-1) + i] |=  mask0[asym_size*(asym_size-j) + i];
            			}
            		}
            	}
            	mask0 = mask1.clone();
        		if (debugLevel > 1) {
        			System.out.println("mask0: (n="+n+"), asym_size="+asym_size+" mask0.length="+mask0.length);
        			for (int i=0;i<asym_size;i++){
            			System.out.print(i+": ");
            			for (int j=0;j<asym_size;j++){
                			System.out.print((mask0[i*asym_size+j]?" X":" .")+" ");
            			}
            			System.out.println();	
        				
        			}
        			System.out.println("calcKernels() numPixels="+numPixels);
        		}
        	}
    		if (debugLevel >0) {
    			System.out.println("mask/mask0: , asym_size="+asym_size+" mask0.length="+mask0.length);
    			for (int i=0;i<asym_size;i++){
        			System.out.print(i+": ");
        			for (int j=0;j<asym_size;j++){
            			System.out.print((mask[i*asym_size+j]?" O":(mask0[i*asym_size+j]?" X":" ."))+" ");
        			}
        			System.out.println();	
    				
    			}
    			System.out.println("calcKernels() numPixels="+numPixels);
    		}
        	
        	
        	
        	ArrayList<Integer> asym_candidates = new ArrayList<Integer>();
        	for (int i = 0; i<mask.length; i++) if (mask0[i] && !mask[i]) asym_candidates.add(new Integer(i));
        	double [] results= new double[asym_candidates.size()];
        	
			System.out.println("asym_candidates.size()="+asym_candidates.size()+" asym_start="+asym_start);
			
        	for (int ncand = 0; ncand<asym_candidates.size();ncand++){
        		mask0 = mask.clone();
        		mask0[asym_candidates.get(ncand)] = true;
            	this.currentVector = initialVector.clone();
            	for (int i = 0; i<mask0.length; i++) if (!mask0[i]) this.currentVector[asym_start+i] = Double.NaN;

        		if (debugLevel >0) {
        			System.out.println("+++ mask0: asym_size="+asym_size+" mask0.length="+mask0.length);
        			for (int i=0;i<asym_size;i++){
            			System.out.print(i+": ");
            			for (int j=0;j<asym_size;j++){
                			System.out.print((mask0[i*asym_size+j]?" X":" .")+" ");
            			}
            			System.out.println();	
        				
        			}
        			System.out.println("calcKernels() numPixels="+numPixels);
        		}
            	
            	
        		if (debugLevel >1) {
        			System.out.println("Before calling levenbergMarquardt  numPixels = "+numPixels+" ncand="+ncand);
        		}
            	
            	
            	OK = levenbergMarquardt();
            	num_lma++;
        		if (debugLevel >1) {
        			System.out.println("After calling levenbergMarquardt, OK="+OK+" numPixels = "+numPixels+" ncand="+ncand);
        			
        		}
        		if (debugLevel >2) {
        			for (int i=asym_start;i< this.currentVector.length;i++){
            			System.out.println("currentVector["+i+"]="+currentVector[i]);
        			}
        		}
            	results[ncand] = this.currentRMSPure;
        	}
        	int best_cand=0;
        	for (int ncand = 1; ncand<asym_candidates.size();ncand++){
        		if (results[ncand] < results[best_cand]){
        			best_cand = ncand;
        		}
        	}
        	bestRms[numPixels] = results[best_cand];
        	enPixels[numPixels] =asym_candidates.get(best_cand);
        	if (results[best_cand] < this.goal_rms_pure) {
        		if (debugLevel >0) {
        			System.out.println("Reached goal at numPixels="+(numPixels+1)+ " results["+best_cand+"]= "+results[best_cand]+" < "+this.goal_rms_pure);
        		}
        		numPixels++;
        		break;
        	}
        	
    	}
    	
// Re-run with the best settings
    	for (int i = 0; i < mask.length; i++) mask[i] = false;
    	for (int i = 0; i < numPixels; i++) mask[enPixels[i]] = true;
    	this.currentVector = initialVector.clone();
    	for (int i = 0; i<mask.length; i++) if (!mask[i]) this.currentVector[asym_start+i] = Double.NaN;

		if (debugLevel >0) {
			System.out.println("Final mask");
			for (int i=0;i<asym_size;i++){
    			System.out.print(i+": ");
    			for (int j=0;j<asym_size;j++){
        			System.out.print((mask[i*asym_size+j]?" X":" .")+" ");
    			}
    			System.out.println();	
				
			}
			System.out.println("calcKernels() numPixels="+numPixels);
		}
    	
    	OK = levenbergMarquardt();
   	
		if (debugLevel >0) {
			for (int i = 0; i< numPixels; i++){
				System.out.println(i+": "+bestRms[i]+" i="+(enPixels[i]/asym_size)+" j="+(enPixels[i]%asym_size));
			}
			System.out.println("Target pure rms is="+this.goal_rms_pure);
		}
		System.out.println("Number of LMA runs = "+num_lma+", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
    	return numPixels;
		
	}
	
	
	
/*
  			properties.setProperty(prefix+"asym_pixels",this.asym_pixels+"");
  			properties.setProperty(prefix+"asym_distance",this.asym_distance+"");
	
 */
	
	
	
	
	
	public double [] getSymKernel(){
		return  getSymKernel(currentVector,sym_kernel_scale);
	}

	public double [] getAsymKernel(){
		return  getAsymKernel(currentVector,1.0/sym_kernel_scale);
	}

	public double [] getConvolved(){ // check that it matches original
		double [] convolved = new double [target_kernel.length];
		System.arraycopy(currentfX, 0, convolved, 0, convolved.length);
		return convolved;
	}
	
	public void setDebugLevel(int debugLevel){
		this.debugLevel =debugLevel;
	}
	
	public void setAsymCompactness(
			double compactness_weight,
			int asym_tax_free){
		this.compactness_weight = compactness_weight;
		this.asym_tax_free =      asym_tax_free;
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
			double [] target_kernel, // should be (asym_size + 2*sym_radius-1)**2
			boolean [] asym_mask)    // which of the asymmetrical kernel to use
	{
		int conv_size = asym_size + 2*sym_radius-2;
		int sym_rad_m1 = sym_radius - 1; // 7
		// find center of the target kernel squared value
		double s0=0.0,sx=0.0,sy=0.0;
		double scx=0.0,scy=0.0;
		if (asym_mask == null) {
			asym_mask = new boolean[asym_size*asym_size];
			for (int i = 0; i<asym_mask.length; i ++ ) asym_mask[i] = true;
		}
		for (int i = 0; i < conv_size; i++){
			for (int j = 0; j < conv_size; j++){
				double d = target_kernel[conv_size*i+j];
				d  *= d;
				s0 += d;
				sx += d*j;
				sy += d*i;
				scx += d*(j-sym_rad_m1-asym_size/2);
				scy += d*(i-sym_rad_m1-asym_size/2);
			}
		}
		
		int j0= (int) Math.round(sx/s0 - sym_rad_m1); // should be ~ async_center 
		int i0= (int) Math.round(sy/s0 - sym_rad_m1); // should be ~ async_center
		if (debugLevel>1){
			System.out.println("setInitialVector(): scx="+(scx/s0) + " scy="+(scy/s0));
			System.out.println("setInitialVector(): x="+(sx/s0) + " y="+(sy/s0));

			System.out.println("setInitialVector(): conv_size = "+conv_size + " asym_size="+asym_size);
			System.out.println("setInitialVector(): fj0 = "+(sx/s0 - sym_rad_m1)+" j0 = "+j0 );
			System.out.println("setInitialVector(): fi0 = "+(sy/s0 - sym_rad_m1)+" i0 = "+i0 );
		}
		// fit i0,j0 to asym_kernel (it should be larger)
		if ((i0<0) || (i0>=asym_size) || (i0<0) || (i0>=asym_size)){
			System.out.println("Warning: kernel center too far : i="+(i0 - asym_size/2)+", j= "+(j0 - asym_size/2));
			if (i0 <  0) i0=0;
			if (j0 <  0) j0=0;
			if (i0 >= asym_size) i0=asym_size-1;
			if (j0 >= asym_size) j0=asym_size-1;
		}
		Random rand = new Random();
		double [] asym_kernel = new double [asym_size * asym_size];
		for (int i = 0; i < asym_kernel.length; i++) asym_kernel[i] = 0; // 0.001*(rand.nextDouble()-0.5)/(asym_size*asym_size);
		asym_kernel[asym_size*i0+j0] = 1.0;
		
		center_i0 = i0; // center to calcualte compatess odf asymmetrical kernel
		center_j0 = j0; // center to calcualte compatess odf asymmetrical kernel

		if (debugLevel>2){
			for (int i = 0; i < asym_kernel.length; i++) {
				System.out.println("asym_kernel["+i+"] = "+asym_kernel[i]);
			}
		}
		for (int i = 0; i < asym_kernel.length; i++) {
			if (!asym_mask[i]) asym_kernel[i] = Double.NaN;
		}
		
		
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
	/*
	// Remove masked out elements of the parameters vector
	private double [] compressVector(
			double [] kvect){
		int comp_len = kvect.length;
		for (int i = 0; i<kvect.length; i++) if (Double.isNaN(kvect[i])) comp_len--;
		double [] cvect = new double[comp_len];
		this.vector_mask = new boolean[kvect.length];
		int indx = 0;
		for (int i = 0; i<kvect.length; i++){
			this.vector_mask[i] = !Double.isNaN(kvect[i]);
			if (this.vector_mask[i]) cvect[indx++] = kvect[i];
		}
		return cvect;
	}
	private double [] expandVector(
			double [] cvect,
			double [] mvect){ // 'old' parameter vectors where Double.NaN means disabled 
		double [] kvect = new double [mvect.length];
		int indx = 0;
		for (int i = 0; i < mvect.length; i++){
			if (!Double.isNaN(mvect[i])) kvect[i] = cvect[indx++];
			else  kvect[i] = Double.NaN;
		}
		return kvect;
	}
	*/
	private void maskAsymPoint(
			double [] kvect,
			int indx) {
		kvect[indx + sym_radius * sym_radius -1] = Double.NaN;
	}
	
	private void unMaskAsymPoint(
			double [] kvect,
			int indx) {
		kvect[indx + sym_radius * sym_radius -1] = 0.0;
	}
	private int [] getVectorLengths(double [] kvect){ // full lengths, some elements may be Double.NaN (disabled)
		int [] lengths = {0,0,0};// [0] - full length; [1] - asym start, [2] - asym length;
		for (int i = 0; i<sym_radius*sym_radius - 1; i++){
			if (!Double.isNaN(kvect[i])) lengths[1] ++;
		}
		for (int i = sym_radius*sym_radius - 1; i < kvect.length; i++){
			if (!Double.isNaN(kvect[i])) lengths[2] ++;
		}
		lengths[0] = lengths[1] + lengths[2];
		return lengths;
	}
	

			

	private double [] getDerivDelta( // calcualte approximate partial derivative as delta
			double [] kvect,      // parameters vector
			int cindx,             // index of parameter to calculate approximate derivative
			double delta
			){
		int indx=0;
		int j = 0;
		for (indx =0;indx < kvect.length;indx++) if (!Double.isNaN(kvect[indx])){
			if (j== cindx) break;
			j++;
		}
		double [] kvect_inc = kvect.clone();
		double [] fx = getFX( kvect);
		kvect_inc[indx] += delta;
		double [] fx1 = getFX( kvect_inc);
//		if (indx == 63) {
//			System.out.println("----- getDerivDelta(): indx="+indx+" delta="+delta+" kvect["+indx+"]="+kvect[indx]+" kvect_inc["+indx+"]="+kvect_inc[indx]);
//		}
		for (int i = 0; i < fx.length; i++){
			fx[i] = (fx1[i]-fx[i])/delta;
		}
		return fx;
	}
	
	public double [] compareDerivative(
			int indx,
			double delta,          // value to increment parameter by for derivative calculation
			boolean verbose){
		double [] rslt = {0.0,0.0};
		double [][] jacob = getJacobian(this.currentVector);
		double [] deriv = getDerivDelta(this.currentVector, indx, delta);
		double [] fx = getFX(this.currentVector);
		for (int i = 0; i< fx.length; i++) {
			if (verbose) {
				System.out.println(i+": "+(jacob[indx][i]-deriv[i])+ " jacob["+indx+"]["+i+"] = "+jacob[indx][i]+
						" deriv["+i+"]="+deriv[i]+" f["+i+"]="+fx[i]);
			}
			rslt[0]+=(jacob[indx][i]-deriv[i])*(jacob[indx][i]-deriv[i]);
			rslt[1]+=jacob[indx][i]*jacob[indx][i];
		}
		rslt[0] = Math.sqrt(rslt[0]/fx.length);
		rslt[1] = Math.sqrt(rslt[1]/fx.length);
		if (debugLevel>3) {
			System.out.println("rms(jacob["+indx+"][]) = "+rslt[1]+", rms(diff) = "+rslt[0]);
		}
		return rslt;
	}
	
	private double [] getFX(
			double [] kvect) // first - all elements of sym kernel but [0] (replaced by 1.0), then - asym ones. May have Double NaN
	{
		int conv_size =       asym_size + 2*sym_radius-2;
		int asym_start=       sym_radius * sym_radius - 1;
		int sym_radius_m1 =   sym_radius -1;
		int asym_terms_start= conv_size*conv_size;
		double cw =           getCompactWeight();
		int num_pars = 0;
		int [] cind = new int [kvect.length];
		for (int i = 0; i < asym_start ; i++){
			if (Double.isNaN(kvect[i])) cind[i] = -1;
			else cind[i] = num_pars++;
		}
		int num_asym  = -num_pars;
		for (int i = asym_start; i<kvect.length ; i++){
			if (Double.isNaN(kvect[i])) cind[i] = -1;
			else cind[i] = num_pars++;
		}
		num_asym +=num_pars;
//		double [][] jacob = new double [num_pars][asym_terms_start + num_asym];

//		double [] fx = new double [conv_size*conv_size+asym_size*asym_size];
		double [] fx = new double [asym_terms_start + num_asym];

		if (this.debugLevel>3){
			System.out.println("fx(): vector_length=    "+kvect.length);
			System.out.println("fx(): sym_radius=       "+sym_radius);
			System.out.println("fx(): asym_size=        "+asym_size);
			System.out.println("fx(): conv_size=        "+conv_size);
			System.out.println("fx(): fx.length= "+fx.length);
			System.out.println("fx(): asym_start=       "+asym_start);
		}
		for (int i = 0; i < fx.length; i++) fx[i] = 0.0;

		for (int i = -sym_radius_m1; i <= sym_radius_m1; i++) {
			for (int j = -sym_radius_m1; j <=sym_radius_m1; j++) {
				double sd = 1;
				int indx = (((i < 0)? -i : i) * sym_radius + ((j < 0)? -j : j)) - 1;
				if ((indx < 0) || (cind[indx] >= 0)) {// <0 - only one, center of symmetrical == 1
					if (indx >=0 ) { // so cind[indx]>=0
						sd =   kvect[indx];
						indx = cind[indx]; // it can be used only as jacobian first index
						if (indx <0) {
							System.out.println("getJacobian() BUG ");
							continue; // should never happen?
						}
						// now sd - value of the parameter vector, indx - index of the compacted (for only enabled parameters) jacobian
					}
					int base_indx = conv_size * (i + sym_radius_m1) + (j + sym_radius_m1);  
					for (int ia=0; ia < asym_size; ia++) {
						for (int ja=0; ja < asym_size; ja++) {
							int asym_index = asym_size*ia + ja + asym_start; // index of the parameter space, full, not compacted
							if (cind[asym_index] >= 0) {
								int conv_index = base_indx + conv_size * ia + ja;
								fx[conv_index] += sd* kvect[asym_index];

//								if ((cind[asym_index]==63) && (conv_index==74)) {
//									System.out.println("cind["+asym_index+"]="+cind[asym_index]+
//											" conv_index ="+conv_index+" kvect["+asym_index+"]=" +kvect[asym_index]+
//											" fx["+conv_index+"]="+fx[conv_index]); 
//									
//								}
								
							}
						}
					}
				}
			}
		}
		
		int asym_term = asym_terms_start;
		for (int ia=0; ia <asym_size;ia++){
			for (int ja=0;ja< asym_size;ja++){
				int asym_index = asym_size*ia + ja;
				int casym_index = cind[asym_index+asym_start];
				if (casym_index >= 0) {
					int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
					if ((ia - center_i0 <= asym_tax_free) &&
							(center_i0 - ia <= asym_tax_free) &&
							(ja - center_j0 <= asym_tax_free) &&
							(center_j0 - ja <= asym_tax_free)) ir2 = 0;
					fx[asym_term++] = ir2 * cw* kvect[asym_start + asym_index];
				}
			}
		}
		return fx;
	}
	
	
	private double [] getFXOld(
			double [] kvect) // first - all elements of sym kernel but [0] (replaced by 1.0), then - asym ones
	{
		int conv_size =       asym_size + 2*sym_radius-2;
		int asym_start=       sym_radius * sym_radius - 1;
		int sym_radius_m1 =   sym_radius -1;
		int asym_terms_start= conv_size*conv_size;
		double cw =           getCompactWeight();

		double [] fx = new double [conv_size*conv_size+asym_size*asym_size];

		if (this.debugLevel>2){
			System.out.println("fx(): vector_length=    "+kvect.length);
			System.out.println("fx(): sym_radius=       "+sym_radius);
			System.out.println("fx(): asym_size=        "+asym_size);
			System.out.println("fx(): conv_size=        "+conv_size);
			System.out.println("fx(): fx.length= "+fx.length);
			System.out.println("fx(): asym_start=       "+asym_start);
		}
		for (int i = 0; i < fx.length; i++) fx[i] = 0.0;
		for (int i = -sym_radius_m1; i <= sym_radius_m1; i++) {
			for (int j = -sym_radius_m1; j <= sym_radius_m1; j++) {
				int indx = ((i < 0)? -i : i) * sym_radius + ((j < 0)? -j : j);
				double sd = (indx >0)? kvect[indx -1] : 1.0;
				int base_indx = conv_size * (i + sym_radius_m1) + (j + sym_radius_m1);  
				for (int ia=0; ia < asym_size; ia++) {
					for (int ja=0; ja < asym_size; ja++) {
						int asym_index = asym_size*ia + ja + asym_start;
						fx[base_indx + conv_size * ia + ja] += sd* kvect[asym_index];
					}
				}
			}
		}
		for (int ia=0; ia <asym_size;ia++){
			for (int ja=0;ja< asym_size;ja++){
				int asym_index = asym_size*ia + ja;
				int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
				if ((ia - center_i0 <= asym_tax_free) &&
						(center_i0 - ia <= asym_tax_free) &&
						(ja - center_j0 <= asym_tax_free) &&
						(center_j0 - ja <= asym_tax_free)) ir2 = 0;
				fx[asym_index+asym_terms_start] = ir2 * cw* kvect[asym_start + asym_index];
			}
		}
		return fx;
	}
	
	
	private double getCompactWeight(){
		return compactness_weight*sym_kernel_scale; // (asym_size*asym_size*asym_size*asym_size); // use
	}

	private double [][] getJacobian(
			double [] kvect) // some entries may be Double.NaN - skip them  as well as asym_kernel entries in the end
	{
		int conv_size = asym_size + 2*sym_radius-2;
		int asym_start=  sym_radius * sym_radius - 1;
		int sym_radius_m1 = sym_radius -1;
		double cw = getCompactWeight();
		int num_pars = 0;
		int [] cind = new int [kvect.length];
		for (int i = 0; i < asym_start ; i++){
			if (Double.isNaN(kvect[i])) cind[i] = -1;
			else cind[i] = num_pars++;
		}
		int num_asym  = -num_pars;
		for (int i = asym_start; i<kvect.length ; i++){
			if (Double.isNaN(kvect[i])) cind[i] = -1;
			else cind[i] = num_pars++;
		}
		num_asym +=num_pars;
		int asym_terms_start=conv_size*conv_size;
		double [][] jacob = new double [num_pars][asym_terms_start + num_asym];
		if (debugLevel > 3) {
			System.out.println("getJacobian(): asym_terms_start="+asym_terms_start+" num_asym="+num_asym);
			for (int i = asym_start; i<kvect.length ; i++){
				System.out.println("getJacobian(): cind["+i+"]="+cind[i]+", kvect["+i+"]="+kvect[i]);
			}
		}
		for (int i = -sym_radius_m1; i <= sym_radius_m1; i++) {
			for (int j = -sym_radius_m1; j <=sym_radius_m1; j++) {
				double sd = 1;
				int indx = (((i < 0)? -i : i) * sym_radius + ((j < 0)? -j : j)) - 1;
				if ((indx < 0) || (cind[indx] >= 0)) {// <0 - only one, center of symmetrical == 1
					if (indx >=0 ) { // so cind[indx]>=0
						sd =   kvect[indx];
						indx = cind[indx]; // it can be used only as jacobian first index
						if (indx <0) {
							System.out.println("getJacobian() BUG ");
							continue; // should never happen?
						}
						// now sd - value of the parameter vector, indx - index of the compacted (for only enabled parameters) jacobian
					}
					int base_indx = conv_size * (i + sym_radius_m1) + (j + sym_radius_m1); // base index in the convolved space 
					for (int ia=0; ia < asym_size; ia++) {
						for (int ja=0; ja < asym_size; ja++) {
							int asym_index = asym_size*ia + ja + asym_start; // index of the parameter space, full, not compacted
							if (cind[asym_index] >= 0) {
								int conv_index = base_indx + conv_size * ia + ja;
								if (indx >= 0) jacob[indx][conv_index] += kvect[asym_index];
								jacob[cind[asym_index]][conv_index] += sd;
								if (debugLevel > 3) {
									System.out.println("getJacobian: indx="+indx+" asym_index="+asym_index+" cind[asym_index]="+cind[asym_index]+
											" conv_index="+conv_index+" kvect["+asym_index+"], sd="+sd);
								}
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < jacob.length; i++) for (int j=asym_terms_start; j < jacob[i].length;j++) jacob[i][j] = 0.0;
		int asym_term = asym_terms_start;
		for (int ia=0; ia <asym_size;ia++){
			for (int ja=0;ja< asym_size;ja++){
				int asym_index = asym_size*ia + ja;
				int casym_index = cind[asym_index+asym_start];
				if (casym_index >= 0) {
					int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
					if ((ia - center_i0 <= asym_tax_free) &&
							(center_i0 - ia <= asym_tax_free) &&
							(ja - center_j0 <= asym_tax_free) &&
							(center_j0 - ja <= asym_tax_free)) ir2 = 0;
					jacob[casym_index][asym_term++] = ir2 * cw;
				}
			}
		}
		
		return jacob;
	}

	private double [][] getJacobianOld(
			double [] kvect){ // some entries may be Double.NaN - skip them 

	int conv_size = asym_size + 2*sym_radius-2;
	int asym_start=  sym_radius * sym_radius - 1;
	int sym_radius_m1 = sym_radius -1;
	double cw = getCompactWeight();
	double [][] jacob = new double [kvect.length][conv_size*conv_size + asym_size*asym_size];
	int asym_terms_start=conv_size*conv_size;
	for (int i = -sym_radius_m1; i <= sym_radius_m1; i++) {
		for (int j = -sym_radius_m1; j <=sym_radius_m1; j++) {
			int indx = ((i < 0)? -i : i) * sym_radius + ((j < 0)? -j : j);
			double sd = (indx >0)? kvect[indx -1] : 1.0;
			int base_indx = conv_size * (i + sym_radius -1) + (j + sym_radius -1);  
			for (int ia=0; ia < asym_size; ia++) {
				for (int ja=0; ja < asym_size; ja++) {
					int asym_index = asym_size*ia + ja+ asym_start;
					int conv_index = base_indx + conv_size * ia + ja;
					if (indx > 0) jacob[indx-1][conv_index] += kvect[asym_index];
					jacob[asym_index][conv_index] += sd;
				}
			}
		}
	}
	for (int i = 0; i < jacob.length; i++) for (int j=asym_terms_start; j < jacob[i].length;j++) jacob[i][j] = 0.0;
	for (int ia=0; ia <asym_size;ia++){
		for (int ja=0;ja< asym_size;ja++){
			int asym_index = asym_size*ia + ja;
			int ir2 = (ia-center_i0)*(ia-center_i0)+(ja-center_j0)*(ja-center_j0);
			if ((ia - center_i0 <= asym_tax_free) &&
					(center_i0 - ia <= asym_tax_free) &&
					(ja - center_j0 <= asym_tax_free) &&
					(center_j0 - ja <= asym_tax_free)) ir2 = 0;
			jacob[asym_index + asym_start][asym_index+asym_terms_start] = ir2 * cw;
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
			double [] fx)            // current convolution result of async_kernel (*) sync_kernel, extended by asym_kernel components
			{
		double [] jTByDiff = new double [jacob.length];
		for (int i=0; i < jTByDiff.length; i++){
			jTByDiff[i] = 0;
			for (int k = 0; k< target_kernel.length; k++){
				jTByDiff[i] += jacob[i][k]*(target_kernel[k]-fx[k]);
			}
			for (int k = target_kernel.length; k< fx.length; k++){
				jTByDiff[i] += jacob[i][k]*(-fx[k]);
			}

		}

		return jTByDiff;
	}
	private double[] getDiffByDiff(
			double [] target_kernel, // target kernel
			double [] fx)            // current convolution result of async_kernel (*) sync_kernel, extended async kernel components
			{
		double [] diffByDiff = {0.0,0.0};
		for (int k=0; k< target_kernel.length; k++){
			double d = target_kernel[k]-fx[k];
			diffByDiff[0] += d*d;
		}
		diffByDiff[1] = diffByDiff[0]; // actual squared error, without compactness 
		for (int k=target_kernel.length; k< fx.length; k++){
			double d = fx[k];
			diffByDiff[0] += d*d;
		}
		if (this.debugLevel > 2){
			System.out.println("getDiffByDiff: diffByDiff[0]="+diffByDiff[0]+" diffByDiff[1]="+diffByDiff[1]);
		}
		return diffByDiff;
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

	private void initLevenbergMarquardt(double fact_precision){
		double s = 0.0;
		for (int i = 0; i<target_kernel.length; i++){
			s+= target_kernel[i]*target_kernel[i];
		}
		this.goal_rms_pure = Math.sqrt(s/target_kernel.length)*fact_precision;
    	this.currentVector = setInitialVector(target_kernel, null); // should be (asym_size + 2*sym_radius-1)**2
	}
	
	private boolean levenbergMarquardt(){
    	long startTime=System.nanoTime();
    	this.firstRMS=-1; //undefined
    	this.iterationStepNumber = 0;
    	this.lambda = this.init_lambda;
// New
    	this.currentfX=null;
		this.lMAArrays=null;
		this.currentRMS=-1;
		this.currentRMSPure=-1;
		this.jacobian=null;  // invalidate
//System.out.println("Setting both lastImprovements(); to -1");			
		lastImprovements[0]=-1.0;
		lastImprovements[1]=-1.0;
    	
    	
    	if (this.numIterations < 0){
			this.currentfX=getFX  (this.currentVector);			
            return true;    		
    	}
    	
		if (this.debugLevel>2) {
			System.out.println("this.currentVector 0");
			for (int i=63;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}

    	while (true) { // loop for the same series
    		boolean [] state=stepLevenbergMarquardtFirst(goal_rms_pure);
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
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()- startTime),3));
    			}
    			long startDialogTime=System.nanoTime();
    			startTime+=(System.nanoTime()-startDialogTime); // do not count time used by the User.
    			//    			if (this.showThisImages) showDiff (this.currentfX, "fit-"+this.iterationStepNumber);
    			//    			if (this.showNextImages) showDiff (this.nextfX,    "fit-"+(this.iterationStepNumber+1));
    		} else if (this.debugLevel>1){
    			System.out.println("==> LevenbergMarquardt(): before action step ="+this.iterationStepNumber+
    					", RMS="+IJ.d2s(this.currentRMS,8)+
    					" ("+IJ.d2s(this.firstRMS,8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
    		}
    		stepLevenbergMarquardtAction(startTime); // apply step - in any case?
//    		if ((this.debugLevel>0) && ((this.debugLevel>1) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
       		if ((this.debugLevel>0) && ((this.debugLevel>0) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
    			System.out.println("--> LevenbergMarquardt(): step = "+this.iterationStepNumber+
    					", RMS="+IJ.d2s(this.currentRMS,8)+
    					" ("+IJ.d2s(this.firstRMS,8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
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
    			") at "+ IJ.d2s(0.000000001*(System.nanoTime()- startTime),3));
    	if (this.debugLevel>2) {
    		double worstRatio = 0;
    		int worstIndex = -1;
    		int param_index=0;
    		for (int i = 0; i<currentVector.length;i++) {
    			if (!Double.isNaN(currentVector[i])){
    				double [] r=compareDerivative(
    						param_index++,
    						0.0000001, // delta,          // value to increment parameter by for derivative calculation
    						false); // param_index>61); //false); // verbose)
    				if (r[1] > 0){
    					if (r[0]/r[1] > worstRatio){
    						worstRatio = r[0]/r[1];
    						worstIndex = i;
    					}

    				}
    			}
    		}
    		System.out.println(" rms(relative diff["+worstIndex+"]) = >>>>> "+worstRatio+" <<<<<");
    	}
    	return true; // all series done
	}
	
	private boolean [] stepLevenbergMarquardtFirst(double goal_rms_pure){
		int deltas_indx;
		double [] deltas=null;
		double [] rmses; // [0]: full rms, [1]:pure rms
		// moved to caller
/*		
		if (this.currentVector==null) {
			this.currentRMS=-1;
			this.currentRMSPure=-1;
			this.currentfX=null; // invalidate
			this.jacobian=null;  // invalidate
			this.lMAArrays=null;
System.out.println("Setting both lastImprovements(); to -1");			
			lastImprovements[0]=-1.0;
			lastImprovements[1]=-1.0;
		}
*/		
		// calculate  this.currentfX, this.jacobian if needed
		if (this.debugLevel>3) {
			System.out.println("this.currentVector 1");
			for (int i=0;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		if (this.debugLevel>2) {
			System.out.println("this.currentVector 2");
			for (int i=63;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		if ((this.currentfX==null)|| (this.lMAArrays==null)) {
			
			this.currentfX=getFX  (this.currentVector);			
			this.jacobian = getJacobian (this.currentVector);
			
//System.out.println("stepLevenbergMarquardtFirst(): this.jacobian.length="+this.jacobian.length);
//for (int i = 0;i<this.currentVector.length;i++){
//	System.out.println(i+": "+this.currentVector[i]);
//}	
			if (debugLevel > 2){
				int conv_size =       asym_size + 2*sym_radius-2;
				int asym_terms_start= conv_size*conv_size;
				for (int i=asym_terms_start; i<currentfX.length; i++){
					System.out.println("this.currentfX["+i+"]="+this.currentfX[i]);
				}
				for (int n=63; n<this.jacobian.length;n++){
//					for (int i=asym_terms_start; i<currentfX.length; i++){
					for (int i=0; i<currentfX.length; i++){
						System.out.println("this.jacobian["+n+"]["+i+"]="+this.jacobian[n][i]+" this.currentfX["+i+"]="+this.currentfX[i]);
					}
				}
			}			
			this.lMAArrays= new LMAArrays();
			lMAArrays.jTByJ=getJTByJ(
					this.jacobian);
			
			lMAArrays.jTByDiff=getJTByDiff(
					this.jacobian,
					this.target_kernel,   // target kernel to factor
					this.currentfX); // used to control compactness of asym_kernel

			if (debugLevel >2) {
				for (int n=63; n < this.lMAArrays.jTByJ.length;n++){
					for (int i=0; i < this.lMAArrays.jTByJ.length; i++){
						System.out.println("this.lMAArrays.jTByJ["+n+"]["+i+"]="+this.lMAArrays.jTByJ[n][i]);
					}
					System.out.println("this.lMAArrays.jTByDiff["+n+"]="+this.lMAArrays.jTByDiff[n]);
				}
			}
			rmses = getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX);
			this.currentRMSPure=  Math.sqrt(rmses[1] / target_kernel.length);
			this.currentRMS =     Math.sqrt(rmses[0] / (asym_size*asym_size+target_kernel.length));
			if (debugLevel > 1){
				System.out.println("currentRMSPure= "+ currentRMSPure + " getDiffByDiff[1] = "+Math.sqrt(getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX)[1] / target_kernel.length));
			}
			
			if (this.debugLevel>1) {
				System.out.println("initial RMS="+IJ.d2s(this.currentRMS,8)+
						" ("+IJ.d2s(this.currentRMSPure,8)+")"+
						". Calculating next Jacobian. Points:"+this.target_kernel.length+" Parameters:"+this.currentVector.length);
			}
		} else {
			rmses = getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX);
			this.currentRMSPure=  Math.sqrt(rmses[1] / target_kernel.length);
			this.currentRMS =     Math.sqrt(rmses[0] / (asym_size*asym_size+target_kernel.length));
			
			if (debugLevel > 2){
				System.out.println("this.currentRMS=" + this.currentRMS+ " getDiffByDiff[1] = "+Math.sqrt(getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX)[1] / target_kernel.length));
			}
			
			
		}
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
			System.out.println("--- deltas ---");
			for (int i=63;i<deltas.length;i++){
				System.out.println("deltas["+i+"]="+ deltas[i]+" this.currentRMS="+this.currentRMS);
			}
		}

		if (this.debugLevel>2) {
			System.out.println("deltas");
			for (int i=0;i<deltas.length;i++){
				System.out.println(i+": "+ deltas[i]);
			}
		}

		// apply deltas    	
		this.nextVector=this.currentVector.clone();
//		System.out.println("this.nextVector.length="+this.nextVector.length+"  deltas.length="+deltas.length);
//		for (int i=0;i<this.nextVector.length;i++) {
//			System.out.println(i+": "+this.nextVector[i]);
//		}
		
		deltas_indx = 0;
		for (int i=0;i<this.nextVector.length;i++) {
			if (!Double.isNaN(this.nextVector[i])){
				if (this.debugLevel>2) {
					if (i>=63) System.out.println(i+": "+this.nextVector[i]+" deltas["+deltas_indx+"]="+deltas[deltas_indx]+
							" this.currentRMS="+this.currentRMS+" this.currentRMSPure="+this.currentRMSPure);
				}
//				System.out.println(deltas[deltas_indx]);
				this.nextVector[i]+=deltas[deltas_indx++];
			}
		}
		
		// another option - do not calculate J now, just fX. and late - calculate both if it was improvement    	
		//    	save current Jacobian
		if (this.debugLevel>2) {
			System.out.println("this.nextVector");
			for (int i=0;i<this.nextVector.length;i++){
				System.out.println(i+": "+ this.nextVector[i]);
			}
		}

		this.savedLMAArrays=lMAArrays.clone();
		this.jacobian=null; // not needed, just to catch bugs
		this.nextfX=getFX  (this.nextVector);			
		
		this.jacobian = getJacobian(this.currentVector);

		this.lMAArrays= new LMAArrays();
		lMAArrays.jTByJ=getJTByJ(
				this.jacobian);
		lMAArrays.jTByDiff=getJTByDiff(
				this.jacobian,
				this.target_kernel,    // target kernel to factor
				this.nextfX);           // next convolution result of async_kernel (*) sync_kernel

		rmses = getDiffByDiff(
				this.target_kernel, // target kernel
				this.nextfX);     // current convolution result of async_kernel (*) sync_kernel
		this.nextRMSPure=  Math.sqrt(rmses[1] / target_kernel.length);
		this.nextRMS =     Math.sqrt(rmses[0] / (asym_size*asym_size+target_kernel.length));

		if (debugLevel > 2){
			System.out.println("nextRMSPure= "+ nextRMSPure + " target_kernel.length = "+target_kernel.length+" getDiffByDiff[1] = "+Math.sqrt(getDiffByDiff(
				this.target_kernel, // target kernel
				this.nextfX)[1] / target_kernel.length));     // current convolution result of async_kernel (*) sync_kernel
		}
		this.lastImprovements[1]=this.lastImprovements[0];
		this.lastImprovements[0]=this.currentRMS-this.nextRMS;
//		System.out.println("Setting this.lastImprovements[1]="+this.lastImprovements[1]+" new this.lastImprovements[0]="+this.lastImprovements[0]);
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
//			System.out.println("this.lastImprovements[0]="+this.lastImprovements[0]+" this.lastImprovements[1]="+this.lastImprovements[1]);
//			System.out.println("this.thresholdFinish="+this.thresholdFinish+
//					" this.thresholdFinish*this.currentRMS="+(this.thresholdFinish*this.currentRMS));
			status[1]=(this.iterationStepNumber>this.numIterations) || ( // done
					(this.lastImprovements[0]>=0.0) &&
					(this.lastImprovements[0]<this.thresholdFinish*this.currentRMS) &&
					(this.lastImprovements[1]>=0.0) &&
					(this.lastImprovements[1]<this.thresholdFinish*this.currentRMS));
			if (!status[1] && (this.currentRMSPure < goal_rms_pure)) {
				status[1] = true;
				if (this.debugLevel>1) {
					System.out.println("Improvent is possible, but the factorization precision reached its goal");
					System.out.println(
							"stepLMA this.currentRMS="+this.currentRMS+
							", this.currentRMSPure="+this.currentRMSPure+
							", this.nextRMS="+this.nextRMS+
							", this.nextRMSPure="+this.nextRMSPure+
							", delta="+(this.currentRMS-this.nextRMS)+
							", deltaPure="+(this.currentRMSPure-this.nextRMSPure));
				}				
			}
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
	private void stepLevenbergMarquardtAction(long startTime){// 
    	this.iterationStepNumber++;
// apply/revert,modify lambda    	
		if (this.debugLevel>1) {
			System.out.println(
					"stepLevenbergMarquardtAction() step="+this.iterationStepNumber+
					", this.currentRMS="+this.currentRMS+
					", this.currentRMSPure="+this.currentRMSPure+
					", this.nextRMS="+this.nextRMS+
					", this.nextRMSPure="+this.nextRMSPure+
					" lambda="+this.lambda+" at "+IJ.d2s(0.000000001*(System.nanoTime()- startTime),3)+" sec");
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
