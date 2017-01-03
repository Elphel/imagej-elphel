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
	public boolean   new_mode =      false; // trying new version of LMA
	public int       asym_size =     6;
	public int       sym_radius =    8; // 2*2^n  - for DCT
	public double[]  target_kernel = null; // should be expanded to 2*(sym_radius)+asym_size- 1 in each direction
	public double    target_rms;   // Target kernel rma (to compare with residual error)
	public int       debugLevel =    3;
	public double    init_lambda =   0.001;
	public int       asym_pixels =  10;      // maximal number of non-zero pixels in asymmmetrical kernel
	public int       asym_distance = 2;      // how far to seed a new pixel
	
	public double    compactness_weight = 1.0; // relative "importance of asymmetrical kernel being compact"
	public int       asym_tax_free = 1; // do not apply compactness_weight for pixels close to the center
	
	
	public double    lambdaStepUp=   8.0; // multiply lambda by this if result is worse
	public double    lambdaStepDown= 0.5; // multiply lambda by this if result is better
//	public double    thresholdFinish=0.001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
	public double    thresholdFinish=0.00001; // (copied from series) stop iterations if 2 last steps had less improvement (but not worsening ) 
	public int       numIterations=  100; // maximal number of iterations
	public double    maxLambda=      1000.0;  // max lambda to fail
	public boolean   stopOnFailure = true;
	
	
	public double    sym_kernel_scale =1.0;//  to scale sym kernel to match original one
	public int       center_i0;
	public int       center_j0;
	public double    center_y;
	public double    center_x;
	
	public double    lambda =     init_lambda;
	public double    currentRMS ;
	public double    firstRMS;
	public double    nextRMS ;
	public double [] RMSes = {-1.0,-1.0};
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
	public LMAData   lMAData = null;
	public int       numLMARuns = 0;
	public int       target_window_mode = 2; // 0 - none, 1 - square, 2 - sin, 3 - sin^2
	public boolean   centerWindowToTarget = true; // center convolution weights around target kernel center
	
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
	public class LMAData{
		public  int          sym_radius=      0;
		public  int          asym_size=       0;
		private double  [][] kernels =        {null,null}; // 0 - sym kernel (including 0 [(sym_radius*2-1)*(sym_radius*2-1)]), 1 - asym kernel
		private double  [][] savedKernels =   null;
		private boolean [][] kernel_masks =   {null,null};
		private int     [][] map_from_pars =  null; // first index - number in (compressed) parameter vector,  value - a pair of {kernel_type, index}  
		private int     [][] map_to_pars =    null; // first index - kernel type, second - index in kernel, value index in parameter vector
		public  double  []   fX =             null; // make always fixed length (convolution plus asym), get rid of map_*_fX
		public  double  []   weight =         null; // same length as fX - combineds weights (used while calculating JTByJ, JTByDiff, DiffByDiff
		                                            // when calculating weight2 - use kernel_masks[1] (zero if false), weights[1] contains full
													// normalized so that sum == 1.0
		public  double  [][] weights =        {null, null}; // same length as fX [0] - weighs for the convolution array, [1] - for asym_kernel
/*?*/   public  double       weight_pure =    0.0;  // 0.. 1.0 - fraction of weights for the convolution part
		
		public  double  []   saved_fX =       null;
		public  double  [][] jacobian =       null;
		private double  []   target_kernel =  null;
		private boolean []   fx_mask =        null;
//		private double  []   asym_weights =   null; // multiply by asym_kernel elements to enforce compactness (proportional to r2 from the center)
//		private int     [][] map_from_fx =    null; // first index - element in fX vector, second [0] - where to look  (0 - target, 1 - asym),
		                                         // second - index in target kernel or asym_kernel (kernels[1])
//		private int     [][] map_to_fx =      null;
		private double  []   par_vector =     null;
		private LMAArrays    lMAArrays =      null;
		public  int          debugLevel =        1;
		public double   []   rmses =         {-1.0, -1.0}; // [0] - composite, [1] - "pure" (without extra "punishment")
		public double   []   saved_rmses =   null;
		public double   []   first_rmses =   null;
		public int           target_window_mode = 2; // 0 - none, 1 - square, 2 - sin
		public boolean       centerWindowToTarget = true; // center convolution weights around target kernel center
		
		
		public LMAData(){
		}
		public LMAData( int debugLevel){
			this.debugLevel = debugLevel;
		}
		public void setTargetWindowMode(int mode, boolean   centerWindowToTarget){
			this.target_window_mode =    mode;
			this.centerWindowToTarget = 	centerWindowToTarget;
		}

		public void initLMAData(){
		}
		public void initLMAData( int debugLevel){
			this.debugLevel = debugLevel;
		}
		public void invalidate_maps_pars(){
			map_from_pars = null; // will need to be re-generated
			map_to_pars = null;
//			nvalidate_weight(); // not needed - should be done for asym only 
		}

		public void invalidate_weight(){
			weight = null;
		}

		public void setSymKernel (double [] sym_kernel){ // does not set mask! Use ( , null) to set default one
			kernels[0] = sym_kernel.clone();
			sym_radius = (int) Math.round(Math.sqrt(sym_kernel.length));
			boolean hasNaN = false;
			for (int i = 0; i < sym_kernel.length; i++) {
				hasNaN =  Double.isNaN(kernels[0][i]);
				if (hasNaN) break;
			}
            if ((kernel_masks[0]==null) || (kernel_masks[0].length != kernels[0].length) || hasNaN){
            	kernel_masks[0] = new boolean [kernels[0].length];
            	kernel_masks[0][0] = false; // do not adjust center element
				for (int i=1;i< kernel_masks[0].length;i++) kernel_masks[0][i] =  hasNaN? (!Double.isNaN(kernels[0][i])): true;
				invalidate_maps_pars();
				invalidate_weight();
			}
            if (hasNaN){
				for (int i=0; i< kernels[0].length; i++) if (Double.isNaN(kernels[0][i])) kernels[0][i] = 0.0;
            }
            
		}
		public void updateSymKernel (double [] sym_kernel){ // just change values (and reset fX)
			kernels[0] = sym_kernel.clone();
			fX=null;
		}

		public void updateAsymKernel (double [] asym_kernel){ // just change values (and reset fX)
			kernels[1] = asym_kernel.clone();
			fX=null;
		}

		public void setSymKernel (double [] sym_kernel, // if null - set mask only
				                  boolean [] mask)
		{
			if (mask == null)      kernel_masks[0] = null; // setSymKernel() will generate default  
			if (sym_kernel !=null) setSymKernel(sym_kernel);
			if (mask != null)	   {
				kernel_masks[0] = mask.clone();
				invalidate_maps_pars();
				invalidate_weight();
			}
		}

		public double[] getSymKernel(){
			return kernels[0];
		}
		
		public double[] getSymKernel(double scale){
			return getScaledKernel(0,scale);
		}

		public double[] getScaledKernel(int kernType, double scale){
			double [] kernel = new double [kernels[kernType].length];
			for (int i=0; i< kernel.length; i++) kernel[i] = Double.isNaN(kernels[kernType][i])? 0.0: scale*kernels[kernType][i];
			return kernel;
		}
		
		public boolean[] getSymKernelMask(){
			return kernel_masks[0]; 
		}

		public void setAsymKernel (double [] asym_kernel){ // does not set mask! Use ( , null) to set default one
			kernels[1] = asym_kernel.clone();
			asym_size = (int) Math.round(Math.sqrt(asym_kernel.length));
            if (debugLevel > 3){
    			System.out.println("setAsymKernel(): kernel_masks[1] is "+((kernel_masks[1]==null)?"":"not ")+"null");
    			if (kernel_masks[1]!=null) System.out.println("kernel_masks[1].length= "+kernel_masks[1].length);
            }
			boolean hasNaN = false;
			for (int i = 0; i < asym_kernel.length; i++) {
				hasNaN =  Double.isNaN(kernels[1][i]);
				if (hasNaN) break;
			}
            if ((kernel_masks[1]==null) || (kernel_masks[1].length != kernels[1].length) || hasNaN){
            	kernel_masks[1] = new boolean [kernels[1].length];
				for (int i=0;i< kernel_masks[1].length;i++) kernel_masks[1][i] = hasNaN? (!Double.isNaN(kernels[1][i])): true;
				invalidate_maps_pars();
				invalidate_weight();
			}
            if (hasNaN){
				for (int i=0; i< kernels[1].length; i++) if (Double.isNaN(kernels[1][i])) kernels[1][i] = 0.0;
            }
            
            if ((weights[1] == null) || (weights[1].length != kernels[1].length)){
            	weights[1] = new double[kernels[1].length];
            	for (int i=0; i<weights[1].length; i++){
            		weights[1][i] = 0.0; // disable
            	}
            }
            
		}

		public void setAsymKernel (double [] asym_kernel, // if null - set mask only
				                   boolean [] mask)
		{
			if (mask == null)       kernel_masks[1] = null; // setSymKernel() will generate default  
			if (asym_kernel !=null) setAsymKernel(asym_kernel);
			if (mask != null){
				kernel_masks[1] = mask.clone();
				invalidate_weight();
			}
			invalidate_maps_pars();
		}		

		public void updateAsymKernelMask (boolean [] mask) // new mask, should not be null
		{
			for (int i = 0; i<mask.length; i++){
				if (!mask[i] || !kernel_masks[1][i]) kernels[1][i] = 0.0;
			}
		    kernel_masks[1] = mask.clone();
			invalidate_maps_pars();
			invalidate_weight();
			if (debugLevel > 2) System.out.println("updateAsymKernelMask(): invalidated map_from_pars and map_to_pars");
		}		
		
		public double[] getAsymKernel(double scale){
			return getScaledKernel(1,scale);
		}
		
		public double[] getAsymKernel(){
			return kernels[1]; 
		}

		public boolean[] getAsymKernelMask(){
			return kernel_masks[1]; 
		}
		
		public void setTarget(double [] target_kernel){
			this.target_kernel = target_kernel.clone();
			fx_mask = new boolean[this.target_kernel.length];
			for (int i= 0;i< fx_mask.length;i++) fx_mask[i] = true;
			invalidate_weight();
		}

		public double [] getTarget(){
			return this.target_kernel;
		}

		public void setTargetMask(boolean [] mask){
			fx_mask=mask.clone();
			invalidate_weight();
		}
		public boolean[] getTargetMask(){
			return fx_mask; 
		}
		public void setAsymWeights(double [] weights){
			this.weights[1] = weights.clone(); // should have the same dimensions
		}

		public void setTargetWeights(double [] weights){
			this.weights[0] = weights.clone(); // should have the same dimensions
		}
		public double[] getTargetWeights(){
			if (this.weights[0] == null) getWeight(false); // generate
			return this.weights[0];
		}
		
		public void rebuildMapsPars(boolean force)
		{
			if (force || (map_from_pars == null)){
				
				par_vector = null; // invalidate
				int numPars=0;
				for (int n = 0; n < kernel_masks.length; n++){
					for (int i = 0; i<kernel_masks[n].length; i++){ // will throw if masks are not initialized
						if (kernel_masks[n][i]) numPars++;
					}
				}
				map_from_pars = new int [numPars][2];
				int indx = 0;
				for (int n = 0; n < kernel_masks.length; n++){
					for (int i = 0; i<kernel_masks[n].length; i++){
						if (kernel_masks[n][i]){
							map_from_pars[indx  ][0] = n;
							map_from_pars[indx++][1] = i;
							
						}
					}
					
				}
				if (debugLevel > 4){
					System.out.println("rebuildMapsPars("+force+"), map_from_pars");
					for (int i = 0; i< map_from_pars.length; i++){
						System.out.println(i+": ("+map_from_pars[i][0]+","+map_from_pars[i][1]+")");
						
					}
				}
			}
			
			if (force || (map_to_pars == null)){
				map_to_pars = new int [kernel_masks.length][];
				int numPar = 0;
				for (int n = 0; n < map_to_pars.length; n++){
					map_to_pars[n] = new int [kernel_masks[n].length];
					for (int i = 0; i < map_to_pars[n].length; i++ ){
						if (kernel_masks[n][i]){
							map_to_pars[n][i] = numPar++;
						} else {
							map_to_pars[n][i] = -1;
						}
						
					}
				}
				if (debugLevel > 4){
					System.out.println("rebuildMapsPars("+force+"), map_to_pars");
					for (int n = 0; n < map_to_pars.length; n++){
						for (int i = 0; i < map_to_pars[n].length; i++ ){
							System.out.println(n+","+i+": "+map_to_pars[n][i]);
						}
					}
				}
				
			}
		}
		
		
		
		public double [] getWeight(boolean force)
		{
			if (force || (weight == null)){
				if (debugLevel > 4){
					System.out.println("getWeight(), delete jacobian");
				}
				if ((weights[0] == null) || (weights[0].length != target_kernel.length)){
					int xc = 0;
					int yc = 0;
					int conv_size = asym_size + 2*sym_radius-2;
					int cc = conv_size/2;
					if (this.centerWindowToTarget) {
						double s0=0.0,sx=0.0,sy=0.0;
						for (int i = 0; i < conv_size; i++){
							for (int j = 0; j < conv_size; j++){
								double d = target_kernel[conv_size*i+j];
								d  *= d;
								s0 += d;
								sx += d* (j - cc);
								sy += d* (i - cc);
							}
						}
						xc = (int) Math.round(sx/s0);
						yc = (int) Math.round(sy/s0);
					}
					double [] sins = new double [2*sym_radius-1];
					for (int i = 1; i< 2*sym_radius; i++) sins[i-1] = Math.sin(Math.PI*i/(2.0 * sym_radius));
					weights[0] = new double [target_kernel.length];
					int left_margin = ((asym_size-1)/2) +xc; // inclusive
					int top_margin = ((asym_size-1)/2) + yc; // inclusive
					int right_margin = left_margin + (2 * sym_radius - 1); // exclusive  
					int bottom_margin = top_margin + (2 * sym_radius - 1); // exclusive  
					for (int i = 0; i < conv_size; i++) {
						for (int j = 0; j < conv_size; j++){
							int cindx = i * conv_size + j;
							if ((i >= top_margin) &&  (i < bottom_margin) && (j >= left_margin) &&  (j < right_margin)) {
								weights[0][cindx] = (target_window_mode>=2)?(sins[i-top_margin]*sins[j-left_margin]):1.0;
								if (target_window_mode == 3) {
									weights[0][cindx]*=weights[0][cindx];
								}
							} else {
								weights[0][cindx] = (target_window_mode == 0)? 1.0: 0.0;
							}
						}
					}
				}
//	public boolean   centerWindowToTarget = true; // center convolution weights around target kernel center
				
			//target_window_mode
				rebuildMapsPars(false); // only if it does not exist
				double sw = 0.0;
				for (int i=0; i< weights[0].length; i++) sw+= weights[0][i];
				weight_pure = sw;
				for (int i=0; i< weights[1].length; i++) sw+= weights[1][i];
				weight_pure /= sw;
				this.weight = new double[weights[0].length+weights[1].length];
				for (int i=0; i< weights[0].length; i++) weight[i] = weights[0][i]/sw;
				for (int i=0; i< weights[1].length; i++) weight[i + weights[0].length] = weights[1][i]/sw;
				fX =       null; // invalidate
				jacobian = null;
			}
			return this.weight;
		}

		// just for debugging
		public int [][] getMapToPars()     { return map_to_pars;	}
		public int [][] getMapFromPars()   { return map_from_pars;}
		public int      getNumPars()       { return map_from_pars.length;}
		public int      getNumPoints()     { return weights[0].length+weights[1].length;}
		public int      getNumPurePoints() { return weights[0].length;}
		

		public double [] getVector(){
			rebuildMapsPars(false); // create maps if not current, invalidates par_vector if rebuilds maps
			if (par_vector == null) {
				par_vector = new double [map_from_pars.length];
				for (int i = 0; i < par_vector.length; i++) {
					par_vector[i] = kernels[map_from_pars[i][0]][map_from_pars[i][1]];
				}
			}
			return par_vector;
		}

		public double [] getConvolved(boolean skip_disabled_asym) // consider all masked out asym_kernel elements to be 0
		{
			return getFX(skip_disabled_asym, true);
		}		

		
		public double [] getFX(boolean skip_disabled_asym) // consider all masked out asym_kernel elements to be 0
		{
			return getFX(skip_disabled_asym, false);
		}		

		public double [] getFX(boolean skip_disabled_asym, boolean justConvolved) // consider all masked out asym_kernel elements to be 0
		{
			getWeight(false); // will invalidate (make null) fX if data is not current, otherwise just return last calculated value.
			if ((fX == null) || justConvolved) {
				int conv_size = asym_size + 2*sym_radius-2;
				int conv_len = conv_size * conv_size;
				int sym_rad_m1 = sym_radius - 1; // 7
				int sym_rad2 = 2*sym_radius; // 16
				int sym_rad4 = 4*sym_radius; // 32
				double [] fX = new double [justConvolved? conv_len:  this.weight.length];
				// calculate convolution, for kernels - regardless of kernels enabled/disabled
				// calculate convolution part
				for (int ci =0; ci < conv_size; ci++) for (int cj =0; cj < conv_size; cj++){
					int cindx = ci*conv_size + cj;
//					if (!justConvolved && (this.weight[cindx] == 0.0)){
					if (this.weight[cindx] == 0.0){
						fX[cindx] = 0.0;
					} else { // calculate convolution for ci, cj
						fX[cindx] = 0;
						for (int ai = 0; ai < asym_size; ai ++){
							int si = (ci - ai) - sym_rad_m1;
							if (si < 0) si = -si;
							int sgni = 1;
							if (si > sym_rad2) si = sym_rad4 - si;
							if (si > sym_radius) {
								sgni = -1;
								si = sym_rad2 - si;
							}
							if (si < sym_radius) { // skip si == sym_radius, coefficient is 0 (WA)
								for (int aj = 0; aj < asym_size; aj ++){
									int aindx = ai*asym_size + aj;
									if (!skip_disabled_asym || kernel_masks[1][aindx]){
										int sj = (cj - aj) - sym_rad_m1;
										if (sj < 0) sj = -sj;
										int sgn = sgni;
										if (sj > sym_rad2) sj = sym_rad4 - sj;
										if (sj > sym_radius) {
											sgn = -sgn;
											sj = sym_rad2 - sj;
										}
										if (sj < sym_radius) { // skip sj == sym_radius, coefficient is 0 (WA)
											int sindx = si * sym_radius + sj;
											fX[cindx] += sgn*kernels[0][sindx] * kernels[1][aindx];
										}
									}
								}	
							}
						}
					}
				}
				if (justConvolved) {
					return fX; // do not copy to this.fX
				}
				// calculate asym kernel elements "handicaps"
				if (!justConvolved) {
					for (int ai = 0; ai < asym_size; ai ++){
						for (int aj = 0; aj < asym_size; aj ++){
							int aindx = ai*asym_size + aj;
							int fx_indx = conv_len + aindx; // from index in the asym_kernel to fX vector
							if ((this.weight[fx_indx] != 0.0)) {
								fX[fx_indx] += kernels[1][aindx]; 
							}
						}
					}
				}
				this.fX = fX;
				double [] diffs = getDiffByDiffW();
				this.rmses = new double[2]; 
				this.rmses[0] = Math.sqrt(diffs[0]);
				this.rmses[1] = Math.sqrt(diffs[1]/this.weight_pure);
				if (first_rmses == null) first_rmses = rmses.clone();
			}
			return fX;
		}
		
		public double [][] getJacobian(boolean recalculate, boolean skip_disabled_asym)
		{
			getWeight(false); // will invalidate (make null) fX and jacobian if data is not current
			if (recalculate || (jacobian == null)){
				jacobian = new double [map_from_pars.length][this.weight.length];
				// zero elements?
				// calculate convolution parts, for kernels - regardless of kernels enabled/disabled
				int conv_size = asym_size + 2*sym_radius-2;
				int conv_len = conv_size * conv_size;
				int sym_rad_m1 = sym_radius - 1; // 7
				int sym_rad2 = 2*sym_radius; // 16
				int sym_rad4 = 4*sym_radius; // 32
				// calculate convolution part
				for (int ci =0; ci < conv_size; ci++) for (int cj =0; cj < conv_size; cj++){
					int cindx = ci*conv_size + cj;
					
					if (this.weight[cindx] != 0.0){ // calculate convolution for ci, cj (skip masked out for speed)
						for (int ai = 0; ai < asym_size; ai ++){
							int si = (ci - ai) - sym_rad_m1;
							if (si < 0) si = -si;
							int sgni = 1;
							if (si > sym_rad2) si = sym_rad4 - si;
							if (si > sym_radius) {
								sgni = -1;
								si = sym_rad2 - si;
							}
							if (si < sym_radius) {  // skip si == sym_radius, coefficient is 0 (WA)
								for (int aj = 0; aj < asym_size; aj ++){
									int aindx = ai*asym_size + aj;
									int apar_indx = map_to_pars[1][aindx];
									int sj = (cj - aj) - sym_rad_m1;
									if (sj < 0) sj = -sj;
									int sgn = sgni;
									if (sj > sym_rad2) sj = sym_rad4 - sj;
									if (sj > sym_radius) {
										sgn = -sgn;
										sj = sym_rad2 - sj;
									}
									if (sj < sym_radius) { // skip sj == sym_radius, coefficient is 0 (WA)
										int sindx = si * sym_radius + sj;
//										if (sindx <0){
//											System.out.println(
//													"ci="+ci+" cj="+cj+" si="+si+" sj="+sj+" ai="+ai+" aj="+aj+
//													" sgni="+sgni+" sgn="+sgn+" sym_rad2="+sym_rad2);
//										}
										if (apar_indx >= 0){
											jacobian[apar_indx][cindx] += sgn*kernels[0][sindx];
										}
										int spar_indx = map_to_pars[0][sindx];
										if ((spar_indx>=0) && (!skip_disabled_asym || kernel_masks[1][aindx])){
											jacobian[spar_indx][cindx] += sgn*kernels[1][aindx];
										}
									}	
								}
							}
						}
					}
				}
				// calculate asym kernel elements "handicaps"
				for (int ai = 0; ai < asym_size; ai ++){
					for (int aj = 0; aj < asym_size; aj ++){
						int aindx = ai*asym_size + aj;
						int par_indx = map_to_pars[1][aindx];
						if (par_indx >=0) { 
							jacobian[par_indx][conv_len + aindx] = 1.0; // asym_weights[aindx];
						}
					}
				}
			}

			return jacobian;
		}
		
		public void invalidateLMAArrays(){
			lMAArrays = null;
			savedLMAArrays = null;
		}
		
		public boolean isValidLMAArrays(){
			return lMAArrays != null;
		}
	
		
		public void save(){
			savedKernels = new double[kernels.length][];
			for (int i=0; i<kernels.length;i++) savedKernels[i] = kernels[i].clone();
			saved_fX = fX.clone();
			saved_rmses = rmses.clone();
		}
		
		public void restore(){
			kernels = new double[savedKernels.length][];
			for (int i=0; i<savedKernels.length;i++) kernels[i] = savedKernels[i].clone();
			fX = saved_fX; // no need to clone()
			rmses = saved_rmses; // no need to clone()
		}
		
		public double [] getRMSes(){
			return rmses;
		}

		public double [] getSavedRMSes(){
			if (saved_rmses == null) {
				double [] m1 = {-1.0,-1.0};
				return m1;
			}
			return saved_rmses;
		}
		public double [] getFirstRMSes(){
			if (first_rmses == null) {
				double [] m1 = {-1.0,-1.0};
				return m1;
			}
			return first_rmses;
		}
		public void resetRMSes(){
			first_rmses = null;
			saved_rmses = null;
		}
		
		
		public double [][] getJTByJW(){
			return getJTByJW(true);
		}
		public double [][] getJTByJW(boolean recalculate){
			if (recalculate) {
				if (lMAArrays==null) lMAArrays = new LMAArrays();
				lMAArrays.jTByJ = new double [jacobian.length][jacobian.length];
				for (int i = 0; i < jacobian.length; i++ ){
					for (int j = 0; j < jacobian.length; j++ ){
						if (j<i){
							lMAArrays.jTByJ[i][j] = lMAArrays.jTByJ[j][i];
						} else {
							lMAArrays.jTByJ[i][j] = 0;
							for (int k=0; k< jacobian[i].length; k++){
								lMAArrays.jTByJ[i][j] += jacobian[i][k] * jacobian[j][k]*weight[k];
							}
						}
					}
				}
			}
			return lMAArrays.jTByJ;
		}
		
		//getFX should be ran
		public double [] getJTByDiffW()
		{
			return getJTByDiffW(true);
		}
		
		public double [] getJTByDiffW(boolean recalculate) // current convolution result of async_kernel (*) sync_kernel, extended by asym_kernel components
		{
			if (recalculate) {
				int conv_size = asym_size + 2*sym_radius-2;
				int conv_len = conv_size * conv_size;
				
				if (lMAArrays==null) lMAArrays = new LMAArrays();
				lMAArrays.jTByDiff = new double [jacobian.length];
				for (int i=0; i < lMAArrays.jTByDiff.length; i++){
					lMAArrays.jTByDiff[i] = 0;
					for (int k = 0; k< jacobian[i].length; k++){
						if (k<conv_len) {
							lMAArrays.jTByDiff[i] += jacobian[i][k]*(target_kernel[k]-fX[k])*weight[k];
						} else {
							lMAArrays.jTByDiff[i] += jacobian[i][k]*(-fX[k])*weight[k];
						}
					}
				}
			}
			return lMAArrays.jTByDiff;
		}

		private double[] getDiffByDiffW()
		{
			double [] diffByDiff = {0.0,0.0};
			for (int i = 0; i< this.weights[0].length; i++){
				double d = target_kernel[i]-fX[i];
				diffByDiff[0] += d*d*weight[i];
			}
			diffByDiff[1] = diffByDiff[0];
			for (int i = 0; i < this.weights[1].length; i++){
				double d =                 -fX[i];
				diffByDiff[0] += d*d*weight[i];
			}
			return diffByDiff;
		}
		
		public double getTargetRMSW(){
			double [] w = getTargetWeights(); // will generate if null
			double trms = 0.0;
			double sw = 0.0;
			for (int i = 0; i< w.length; i++){
				double d = target_kernel[i];
				trms += d*d*w[i];
				sw+=   w[i];
			}
			return Math.sqrt(trms/sw);
		}

		
		public double [] solveLMA(
				double lambda,
				int debugLevel){
			double [][] JtByJmod= lMAArrays.jTByJ.clone();

			int numPars=JtByJmod.length;
			for (int i=0;i<numPars;i++){
				JtByJmod[i]=lMAArrays.jTByJ[i].clone();
				JtByJmod[i][i]+=lambda*JtByJmod[i][i]; //Marquardt mod
			}
			//     M*Ma=Mb
			Matrix M=new Matrix(JtByJmod);
			if (debugLevel > 3) {
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
			return Ma.getColumnPackedCopy(); // deltas
		}
		
		public void applyDeltas(double [] deltas)
		{
			for (int i = 0; i< deltas.length; i++){
				kernels[map_from_pars[i][0]][map_from_pars[i][1]] += deltas[i];
			}
			fX = null; //needs to be recalculated
		}
	} // class LMAData
	
	public FactorConvKernel(){
		this.new_mode =true;
	}

	public FactorConvKernel(boolean new_mode){
		this.new_mode =new_mode;
	}
	
	public FactorConvKernel(int asym_size, int sym_radius){
		this.asym_size = asym_size;
		this.sym_radius =  sym_radius;
	}

	public void setTargetWindowMode(int mode, boolean centerWindowToTarget){
		target_window_mode = mode;
		this.centerWindowToTarget = centerWindowToTarget;
	}

	public double[] getRMSes(){
		return this.RMSes; // lMAData.getRMSes();
	}
	
	public double getTargetRMS(){
		return lMAData.getTargetRMSW();
	}
	
	public double [] generateAsymWeights(
			int asym_size,
			double scale,
			double xc,
			double yc,
			int taxfree){
		double [] weight = new double [asym_size*asym_size];
		int ixc = (int) Math.round(xc);
		int iyc = (int) Math.round(yc);
		for (int i=0; i<asym_size; i++){
			for (int j=0; j<asym_size; j++) {
				double r2 = (i-yc)*(i-yc)+(j-xc)*(j-xc);
				if (	(i - iyc <= taxfree) &&
						(j - ixc <= taxfree) &&
						(iyc - i <= taxfree) &&
						(ixc - j <= taxfree)) r2 = 0.0;
				weight [i*asym_size + j] = r2 * scale;
			}
		}
		return weight;
	}
	
	
	public boolean calcKernels(
			double []target_kernel,
			int asym_size,
			int sym_radius,
			double fact_precision,
			boolean [] mask,
			int seed_size,
			double asym_random){
		this.asym_size = asym_size;
		this.sym_radius =  sym_radius;
		this.target_kernel = target_kernel;
    	this.startTime=System.nanoTime();
    	if (new_mode) {
    		initLevenbergMarquardt(fact_precision, seed_size, asym_random);
    		if (mask != null) {
//    			lMAData.setAsymKernel (
//   					null, // double [] asym_kernel, // if null - set mask only
//  					mask);
    			lMAData.updateAsymKernelMask(mask); // this will zero out asym_kernel elements that are masked out
    		}
    		
    		if (debugLevel > 1){
    			if (mask == null){
        			System.out.println("mask is null, retrieving from the kernels");
    				mask = lMAData.getAsymKernelMask();
    			}
    			System.out.println("mask.length="+mask.length);
    			System.out.println("asym mask: ");
    			for (int ii=0;ii < asym_size;ii++){
    				System.out.print(ii+": ");
    				for (int jj=0;jj < asym_size;jj++){
    					System.out.print((mask[ii * asym_size + jj]?" X":" .")+" ");
    				}
    				System.out.println();	
    			}
    		}
    		double [] asym_weights = generateAsymWeights(
    				asym_size,
    				this.compactness_weight * this.sym_kernel_scale, // double scale,
    				this.center_j0, //double xc,
    				this.center_i0, // double yc,
    				this.asym_tax_free); // int taxfree);
    		lMAData.setAsymWeights (asym_weights);
    		if (this.debugLevel > 3) {
    			double [][] kernels = {lMAData.getSymKernel(),lMAData.getAsymKernel()};
    			System.out.println("calcKernels(): kernels data:");
    			for (int n=0;n<kernels.length;n++) for (int i=0;i<kernels[n].length;i++){
    				System.out.println(n+"/"+i+": "+ kernels[n][i]);
    			}
    		}
    		// first - freeze sym_kernel, then unfreeze
    		boolean [] sym_mask = lMAData.getSymKernelMask();
    		boolean [] sym_mask_frozen =new boolean[sym_mask.length];
    		for (int i = 0; i<sym_mask_frozen.length; i++) sym_mask_frozen[i] = false;
    		lMAData.setSymKernel(null, sym_mask_frozen);
    		levenbergMarquardt();
    		lMAData.setSymKernel(null, sym_mask);
    		boolean OK = levenbergMarquardt();
    		this.RMSes = lMAData.getRMSes().clone();
    		return OK;
		} else{
			initLevenbergMarquardt_old(fact_precision, seed_size, asym_random);
			if (mask != null){
				if (this.debugLevel > 3) {
					System.out.println("calcKernels(): this.currentVector 0");
					for (int i=63;i<this.currentVector.length;i++){
						System.out.println(i+": "+ this.currentVector[i]);
					}
				}

				int asym_start=  sym_radius * sym_radius - 1;
				for (int i = 0; i<mask.length; i++) {
					if (Double.isNaN(this.currentVector[asym_start+i])) this.currentVector[asym_start+i] = 0.0; // replace all NaN with 0.0;
					if (!mask[i]) this.currentVector[asym_start+i] = Double.NaN;
				}
				if (debugLevel > 1){
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
			if (this.debugLevel > 3) {
				System.out.println("calcKernels(): this.currentVector 1");
				for (int i=63;i<this.currentVector.length;i++){
					System.out.println(i+": "+ this.currentVector[i]);
				}
			}
			return levenbergMarquardt_old();
		}
	}

	public int calcKernels(
			double []target_kernel,
			int asym_size,
			int sym_radius,
			double fact_precision,
			int asym_pixels,         // maximal number of non-zero pixels in asymmmetrical kernel
			int asym_distance,       // how far to seed a new pixel
			int seed_size){
		if (!new_mode) return calcKernels_old(
				target_kernel,
				asym_size,
				sym_radius,
				fact_precision,
				asym_pixels,          // maximal number of non-zero pixels in asymmmetrical kernel
				asym_distance,       // how far to seed a new pixel
				seed_size);

		this.asym_size =     asym_size;
		this.sym_radius =    sym_radius;
		this.target_kernel = target_kernel;
		this.asym_pixels =   asym_pixels;
		this.asym_distance = asym_distance;
		
		double [] RMSes = null;
    	this.startTime=System.nanoTime(); // need local?
		int numWeakest = 0;
		int numAny=0;

		initLevenbergMarquardt(fact_precision, seed_size,0.0);
		this.goal_rms_pure = lMAData.getTargetRMSW()*fact_precision;
		this.target_rms = lMAData.getTargetRMSW(); //Math.sqrt(s/target_kernel.length);
		
		if (debugLevel > 1){
			boolean []	mask = lMAData.getAsymKernelMask();
			System.out.println("mask.length="+mask.length);
			System.out.println("asym mask: ");
			for (int ii=0;ii < asym_size;ii++){
				System.out.print(ii+": ");
				for (int jj=0;jj < asym_size;jj++){
					System.out.print((mask[ii * asym_size + jj]?" X":" .")+" ");
				}
				System.out.println();	
			}
		}
		double [] asym_weights = generateAsymWeights(
				asym_size,
				this.compactness_weight * this.sym_kernel_scale, // double scale,
				this.center_j0, //double xc,
				this.center_i0, // double yc,
				this.asym_tax_free); // int taxfree);
		lMAData.setAsymWeights (asym_weights);
		if (!levenbergMarquardt()){
			boolean [] asym_mask = lMAData.getAsymKernelMask();
			int numAsym = 0;
			for (int i = 0; i < asym_mask.length; i++) if (asym_mask[i]) numAsym++;
			System.out.println("===== calcKernels(): failed to run first LMA ======");
			return numAsym; 
		}
		RMSes = lMAData.getRMSes();
		// Add points until reached asym_pixels
		int numAsym = 0;
		while (true){
			boolean [] asym_mask = lMAData.getAsymKernelMask();
			numAsym = 0;
			for (int i = 0; i < asym_mask.length; i++) if (asym_mask[i]) numAsym++;
			if (numAsym >= asym_pixels) break;
			if ((RMSes!=null) && (RMSes[1] < this.goal_rms_pure)){
				if (debugLevel > 1){
					System.out.println("calcKernels(): reached goal: "+RMSes[1]+" < "+this.goal_rms_pure);
				}
				break;
			}
			RMSes =addCell(
					asym_distance, // how far to seed a new pixel (hor/vert
					RMSes[0],       // if Double.NaN - will not compare against it, return the best one
					-1); // no skip points
			if (RMSes == null) {
				if (debugLevel > 1) {
					System.out.println("calcKernels(): failed to add cell");
				}
				return numAsym;
			}
		}
		if (debugLevel > 0){
			System.out.println(
					"Finished adding cells, number of LMA runs = "+getLMARuns()+
					", RMS = "+RMSes[0]+
					", RMSPure = "+RMSes[1]+
					", relRMSPure = "+(RMSes[1]/this.target_rms)+
					", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
		}
		if (debugLevel > 0){
			boolean [] am = lMAData.getAsymKernelMask();
			for (int ii=0;ii < asym_size;ii++){
				System.out.print(ii+": ");
				for (int jj=0;jj < asym_size;jj++){
					System.out.print((am[ii * asym_size + jj]?" X":" .")+" ");
				}
				System.out.println();	
			}
		}

		if ((RMSes!=null) && (RMSes[1] < this.goal_rms_pure)){
			if (debugLevel > 0){
				System.out.println("calcKernels(): reached goal: "+RMSes[1]+" < "+this.goal_rms_pure);
			}
		} else  {
			// Replace weakest
			while (true){
				double [] newRMSes = replaceWeakest(
						asym_distance,       // how far to seed a new pixel (hor/vert
						RMSes[0]); // double   was_rms
				if (newRMSes == null){
					break;
				}
				RMSes = newRMSes;
				numWeakest++;
				if (debugLevel>0){
					System.out.println(
							"Replaced weakiest cell ("+numWeakest+"), number of LMA runs = "+getLMARuns()+
							", RMS = "+RMSes[0]+
							", RMSPure = "+RMSes[1]+
							", relRMSPure = "+(RMSes[1]/this.target_rms)+
							", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
				}
				if (debugLevel>0){
					boolean [] am = lMAData.getAsymKernelMask();
					for (int ii=0;ii < asym_size;ii++){
						System.out.print(ii+": ");
						for (int jj=0;jj < asym_size;jj++){
							System.out.print((am[ii * asym_size + jj]?" X":" .")+" ");
						}
						System.out.println();	
					}
				}
				if ((RMSes!=null) && (RMSes[1] < this.goal_rms_pure)){
					if (debugLevel > 0){
						System.out.println("calcKernels(): reached goal: "+RMSes[1]+" < "+this.goal_rms_pure);
					}
					break;
				}
				

			}

			if (debugLevel > 0){
				System.out.println(
						"Finished replacing weakiest cells, number of LMA runs = "+getLMARuns()+ ", number of weakest replaced = "+numWeakest+
						", RMS = "+RMSes[0]+
						", RMSPure = "+RMSes[1]+
						", relRMSPure = "+(RMSes[1]/this.target_rms)+
						", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
			}
			if (debugLevel > 0){
				boolean [] am = lMAData.getAsymKernelMask();
				for (int ii=0;ii < asym_size;ii++){
					System.out.print(ii+": ");
					for (int jj=0;jj < asym_size;jj++){
						System.out.print((am[ii * asym_size + jj]?" X":" .")+" ");
					}
					System.out.println();	
				}
			}
			/*
			if ((RMSes!=null) && (RMSes[1] < this.goal_rms_pure)){
				if (debugLevel > 0){
					System.out.println("calcKernels(): reached goal: "+RMSes[1]+" < "+this.goal_rms_pure);
				}
			} else  {
				// replace any cell
				while (true){
					double [] newRMSes = replaceAny(
							asym_distance,       // how far to seed a new pixel (hor/vert
							RMSes[0]); // double   was_rms
					if (newRMSes == null){
						break;
					}
					RMSes = newRMSes;
					numAny++;
					if (debugLevel > 0){
						System.out.println(
								"Replaced a cell (not the weakest), total number "+numAny+", number of LMA runs = "+getLMARuns()+
								", RMS = "+RMSes[0]+
								", RMSPure = "+RMSes[1]+
								", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
					}
					if (debugLevel > 0){
						boolean [] am = lMAData.getAsymKernelMask();
						for (int ii=0;ii < asym_size;ii++){
							System.out.print(ii+": ");
							for (int jj=0;jj < asym_size;jj++){
								System.out.print((am[ii * asym_size + jj]?" X":" .")+" ");
							}
							System.out.println();	
						}
					}
					if ((RMSes!=null) && (RMSes[1] < this.goal_rms_pure)){
						if (debugLevel > 0){
							System.out.println("calcKernels(): reached goal: "+RMSes[1]+" < "+this.goal_rms_pure);
						}
						break;
					}
					
				}
			}
			*/
		}
		if (debugLevel > 0){
			System.out.println(
					"Finished replacing (any) cells, number of LMA runs = "+getLMARuns()+", number of cells replaced - "+numAny+
					", RMS = "+RMSes[0]+
					", RMSPure = "+RMSes[1]+
					", relRMSPure = "+(RMSes[1]/this.target_rms)+
					", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
		}
		if (debugLevel > 0){
			boolean [] am = lMAData.getAsymKernelMask();
			for (int ii=0;ii < asym_size;ii++){
				System.out.print(ii+": ");
				for (int jj=0;jj < asym_size;jj++){
					System.out.print((am[ii * asym_size + jj]?" X":" .")+" ");
				}
				System.out.println();	
			}
		}
		this.RMSes = RMSes.clone();
		return numAsym;		
	}

	
	public int calcKernels_old(
			double []target_kernel,
			int asym_size,
			int sym_radius,
			double fact_precision,
			int asym_pixels,          // maximal number of non-zero pixels in asymmmetrical kernel
			int asym_distance,       // how far to seed a new pixel
			int seed_size){
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
    	
///    	double[] initialVector = setInitialVector(target_kernel, null); // should be (asym_size + 2*sym_radius-1)**2
    	double [][]kernels = setInitialVector(target_kernel, seed_size, 0.0); // should be (asym_size + 2*sym_radius-1)**2
    	double[] initialVector = setVectorFromKernels(kernels[0], kernels[1]);
    	
		enPixels[0] = center_i0*asym_size+center_j0;
    	boolean [] mask = new boolean [asym_size*asym_size];
    	for (int i = 0; i<mask.length; i++) mask[i] = false;
    	mask[enPixels[0]] = true;
    	this.currentVector = initialVector.clone();
    	System.out.println("mask.length="+mask.length+" asym_start="+asym_start+" this.currentVector.length="+this.currentVector.length);
    	for (int i = 0; i<mask.length; i++) if (!mask[i]) this.currentVector[asym_start+i] = Double.NaN;
    	boolean OK =       levenbergMarquardt_old();
    	num_lma++;
    	bestRms[0] =       this.currentRMSPure; 
/*    	
    	for (int i = 0; i<numPixels; i++) mask[enPixels[i]] = true;
    	this.currentVector = initialVector.clone();
    	for (int i = 0; i<mask.length; i++) if (!mask[i]) this.currentVector[asym_start+1] = Double.NaN;
    	boolean OK = levenbergMarquardt_old();
    	if (numPixels == 1){
        	bestRms[numPixels-1] = this.currentRMSPure; 
    	}
*/		
		int numPixels = 0;
    	for (numPixels = 1; numPixels < asym_pixels; numPixels++) {
    		if (debugLevel > 0) {
    			System.out.println("calcKernels() numPixels="+numPixels);
    		}
    		
        	for (int i = 0; i < mask.length; i++) mask[i] = false;
        	for (int i = 0; i < numPixels; i++) mask[enPixels[i]] = true;
        	boolean diagonal = false;
        	if (asym_distance < 0){
        		asym_distance = -asym_distance;
        		diagonal = true;
        	}
        	boolean[] mask0 = mask.clone();
        	for (int n = 0; n < asym_distance; n++){
            	boolean[] mask1 = mask0.clone();
            	if (diagonal){
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
    		if (debugLevel > 0) {
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

        		if (debugLevel > 0) {
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
            	
            	
        		if (debugLevel > 1) {
        			System.out.println("Before calling levenbergMarquardt_old  numPixels = "+numPixels+" ncand="+ncand);
        		}
            	
            	
            	OK = levenbergMarquardt_old();
            	num_lma++;
        		if (debugLevel > 1) {
        			System.out.println("After calling levenbergMarquardt_old, OK="+OK+" numPixels = "+numPixels+" ncand="+ncand);
        			
        		}
        		if (debugLevel > 2) {
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
        		if (debugLevel > 0) {
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

		if (debugLevel > 0) {
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
    	
    	OK = levenbergMarquardt_old();
   	
		if (debugLevel > 0) {
			for (int i = 0; i< numPixels; i++){
				System.out.println(i+": "+bestRms[i]+" i="+(enPixels[i]/asym_size)+" j="+(enPixels[i]%asym_size));
			}
			System.out.println("Target pure rms is="+this.goal_rms_pure);
		}
		System.out.println("Number of LMA runs = "+num_lma+", spent "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)+" sec");
    	return numPixels;
		
	}
	
	// LMAData contains current kernels/masks (but not RMSes!)
	// if improved - will leave with the new kernels/masks and return RMSes, if failed - return null and restore
	public double [] addCell(
			int      asym_distance,       // how far to seed a new pixel (hor/vert
			double   was_rms, // if Double.NaN - will not compare against it, return the best one
			int      skip_index // to prevent just removed to appear again 
			){
		// save state
		double  [][] saved_kernels = {lMAData.getSymKernel().clone(),lMAData.getAsymKernel().clone()};
		boolean [][] saved_kernel_masks = {lMAData.getSymKernelMask().clone(),lMAData.getAsymKernelMask().clone()};
		boolean [] sym_mask_frozen =new boolean[saved_kernel_masks[0].length];
		for (int i = 0; i<sym_mask_frozen.length; i++) sym_mask_frozen[i] = false;
		
//		double  []   saved_RMSes = lMAData.getRMSes();
		// create list of neighbors
		int [] candidates = expandMask(
				saved_kernel_masks[1],
				asym_distance);
		double [][] results = new double[candidates.length][];
		double [][][] result_kernels = new double[candidates.length][][];
		for (int n = 0; n < candidates.length; n++) if (candidates[n] != skip_index){
			lMAData.setSymKernel  (saved_kernels[0], saved_kernel_masks[0]);
			lMAData.setAsymKernel (saved_kernels[1], saved_kernel_masks[1]);
			boolean [] asym_mask = saved_kernel_masks[1].clone();
			asym_mask[candidates[n]] = true;
			lMAData.updateAsymKernelMask(asym_mask); // will set new asym_kernel values to 0;
// First - just adjust asym
			if (debugLevel > 1) System.out.println("Running asym-only LMA");
    		lMAData.setSymKernel(null, sym_mask_frozen);
    		levenbergMarquardt();
    		if (debugLevel > 1) System.out.println("Running full LMA");
    		lMAData.setSymKernel(null, saved_kernel_masks[0]);
			if ( levenbergMarquardt()) {
				results[n] = lMAData.getRMSes().clone();
				result_kernels[n] = new double[2][];
				result_kernels[n][0] = lMAData.getSymKernel().clone();
				result_kernels[n][1] = lMAData.getAsymKernel().clone();
			} else {
				results[n]= null;
				result_kernels[n] = null;
			}
		}
		int best_index=-1;
		for (int n = 0; n<results.length; n++){
			if ((results[n] !=null) && (Double.isNaN(was_rms) || (results[n][0] < was_rms)) && ((best_index < 0) || (results[n][0] < results[best_index][0]))){
				best_index = n;
			}
		}
		if (best_index <0) { // failed
			lMAData.setSymKernel  (saved_kernels[0], saved_kernel_masks[0]);
			lMAData.setAsymKernel (saved_kernels[1], saved_kernel_masks[1]);
			return null;
		} else {
			boolean [] asym_mask = saved_kernel_masks[1].clone();
			asym_mask[candidates[best_index]] = true;
			lMAData.setSymKernel  (result_kernels[best_index][0], saved_kernel_masks[0]);
			lMAData.setAsymKernel (result_kernels[best_index][1], asym_mask);
			return results[best_index];
		}
	}
	
	// LMAData contains current kernels/masks (but not RMSes!)
	// if improved - will leave with the new kernels/masks and return RMSes, if failed - return null and restore
	// Try to replace the cell that has the lowest absolute value with the different one
	public double [] replaceWeakest(
			int      asym_distance,       // how far to seed a new pixel (hor/vert
			double   was_rms
			){
		// save state
		double  [][] saved_kernels = {lMAData.getSymKernel().clone(),lMAData.getAsymKernel().clone()};
		boolean [][] saved_kernel_masks = {lMAData.getSymKernelMask().clone(),lMAData.getAsymKernelMask().clone()};
		int weakestIndex =-1;
		for (int i = 0; i < saved_kernel_masks[1].length; i++) if (saved_kernel_masks[1][i]){
			if ((weakestIndex < 0) || (Math.abs(saved_kernels[1][i]) < Math.abs(saved_kernels[1][weakestIndex]))) {
				weakestIndex = i;
			}
		}
		if (weakestIndex < 0) return null ; //should not happen
		boolean [] asym_mask = saved_kernel_masks[1].clone();
		asym_mask[weakestIndex] = false;
		lMAData.setAsymKernel (null, asym_mask); // set mask only
		double [] RMSes = addCell( asym_distance, Double.NaN, weakestIndex); // if Double.NaN - will not compare against it, return the best one
//		double [] RMSes = addCell( asym_distance, Double.NaN, -1); // if Double.NaN - will not compare against it, return the best one
		if (RMSes == null){
			System.out.println("replaceWeakest(): Failed to find any replacements at all");
		} else if (RMSes[0] > was_rms){
			if (this.debugLevel > 1){
				System.out.println("replaceWeakest(): Failed to find a better replacemnet ("+ RMSes[0]+">"+was_rms+")");
			}
			RMSes = null;
		}
		if (RMSes == null) { // failed in any way - restore state
			lMAData.setSymKernel  (saved_kernels[0], saved_kernel_masks[0]);
			lMAData.setAsymKernel (saved_kernels[1], saved_kernel_masks[1]);
			return null;
		}
		return RMSes; // keep new kernels, return RMSes 
	}
	
	// LMAData contains current kernels/masks (but not RMSes!)
	// if improved - will leave with the new kernels/masks and return RMSes, if failed - return null and restore
	// Try to replace each of the already enabled cells the different ones
	public double [] replaceAny(
			int      asym_distance,       // how far to seed a new pixel (hor/vert
			double   was_rms
			){
		// save state
		double  [][] saved_kernels = {lMAData.getSymKernel().clone(),lMAData.getAsymKernel().clone()};
		boolean [][] saved_kernel_masks = {lMAData.getSymKernelMask().clone(),lMAData.getAsymKernelMask().clone()};
		int numCells = 0;
		for (int i = 0; i < saved_kernel_masks[1].length; i++) if (saved_kernel_masks[1][i]) numCells++;
		int [] cells = new int [numCells];
		int indx = 0;
		for (int i = 0; i < saved_kernel_masks[1].length; i++) if (saved_kernel_masks[1][i]) cells[indx++] = i;
		double [][] results = new double[cells.length][];
		double [][][] result_kernels = new double[cells.length][][];
		boolean [][] result_asym_mask = new boolean[cells.length][];
		
		for (int removedCell = 0; removedCell < cells.length; removedCell++){
			boolean [] asym_mask = saved_kernel_masks[1].clone();
			asym_mask[cells[removedCell]] = false;
			lMAData.setAsymKernel (null, asym_mask); // set mask only
			results[removedCell] = addCell( asym_distance, Double.NaN,cells[removedCell]); // if Double.NaN - will not compare against it, return the best one
			if (results[removedCell] == null){
				System.out.println("replaceAny(): Failed to find any replacements at all");
			} else if (results[removedCell][0] > was_rms){
				if (this.debugLevel > 2){
					System.out.println("replaceAny(): Failed to find a better replacemnet for cell "+removedCell+" ("+ results[removedCell][0]+">"+was_rms+")");
				}
				results[removedCell] = null;
			}
			if (results[removedCell] != null) {
				result_kernels[removedCell] = new double[2][];
				result_kernels[removedCell][0] = lMAData.getSymKernel().clone();
				result_kernels[removedCell][1] = lMAData.getAsymKernel().clone();
				result_asym_mask[removedCell] =  lMAData.getAsymKernelMask().clone();
			}
		}
		// See if any of the replacements improved result
		int bestToRemove = -1;
		for (int n = 0; n < results.length; n++){
			if ((results[n] != null) && (results[n][0] < was_rms) && ((bestToRemove < 0) || (results[n][0] < results[bestToRemove][0]))) {
				bestToRemove = n;
			}
		}
		if (bestToRemove < 0){
			lMAData.setSymKernel  (saved_kernels[0], saved_kernel_masks[0]);
			lMAData.setAsymKernel (saved_kernels[1], saved_kernel_masks[1]);
			return null;
		}
		
		lMAData.setSymKernel  (result_kernels[bestToRemove][0], saved_kernel_masks[0]);
		lMAData.setAsymKernel (result_kernels[bestToRemove][1], result_asym_mask[bestToRemove]);
		return results[bestToRemove];
	}
	
	
	private int [] expandMask(
			boolean [] mask,
			int asym_distance){
		boolean diagonal = false;
		if (asym_distance < 0){
			asym_distance = -asym_distance;
			diagonal = true;
		}
		boolean[] mask0 = mask.clone();
		for (int n = 0; n < asym_distance; n++){
			boolean[] mask1 = mask0.clone();
			if (diagonal){
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
				System.out.println("expandMask(): (n="+n+"), asym_size="+asym_size+" mask0.length="+mask0.length);
				for (int i=0;i<asym_size;i++){
					System.out.print(i+": ");
					for (int j=0;j<asym_size;j++){
						System.out.print((mask0[i*asym_size+j]?(mask[i*asym_size+j]?" +":" X"):" .")+" ");
					}
					System.out.println();	

				}
			}
		}
    	ArrayList<Integer> asym_candidates = new ArrayList<Integer>();
    	for (int i = 0; i<mask.length; i++) if (mask0[i] && !mask[i]) asym_candidates.add(new Integer(i));
    	int [] front = new int[asym_candidates.size()];
    	for (int i = 0; i<front.length;i++) front[i] = asym_candidates.get(i);
    	return front;
	}
	
	public double [] getSymKernel(){
		if (new_mode) return lMAData.getSymKernel(sym_kernel_scale);
		return  getSymKernel(currentVector, sym_kernel_scale);
	}

	public double [] getAsymKernel(){
		if (new_mode) return lMAData.getAsymKernel(1.0/sym_kernel_scale);
		return  getAsymKernel(currentVector,1.0/sym_kernel_scale);
	}

	public double [] getTargetWeights(){
		if (new_mode) return lMAData.getTargetWeights();
		return  null;
	}

	public double [] getConvolved(){ // check that it matches original
		if (new_mode) {
			return lMAData.getConvolved(true); //boolean skip_disabled_asym
		}
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
	private double [][] setInitialVector(
			double [] target_kernel, // should be (asym_size + 2*sym_radius-1)**2
			int seed_size,    // 4*n +1 - center and 45-degree, 4*n - 4 center and 45-degree cross
			double asym_random)
	{
		int conv_size = asym_size + 2*sym_radius-2;
		int sym_rad_m1 = sym_radius - 1; // 7
		// find center of the target kernel squared value
		double s0=0.0,sx=0.0,sy=0.0;
		double scx=0.0,scy=0.0;
		boolean [] asym_mask = null;
		if (asym_mask == null) {
			asym_mask = new boolean[asym_size*asym_size];
			for (int i = 0; i<asym_mask.length; i ++ ) asym_mask[i] =  seed_size == 0;
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
		double xc = sx/s0 - sym_rad_m1;
		double yc = sy/s0 - sym_rad_m1;
		int j0= (int) Math.round(xc); // should be ~ async_center 
		int i0= (int) Math.round(yc); // should be ~ async_center
		if (debugLevel>0){
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
		double [] asym_kernel = new double [asym_size * asym_size];
		Random rand = new Random();
//		for (int i = 0; i < asym_kernel.length; i++) asym_kernel[i] = 0; // 0.001*(rand.nextDouble()-0.5)/(asym_size*asym_size);
//		for (int i = 0; i < asym_kernel.length; i++) asym_kernel[i] = 0.001*(rand.nextDouble()-0.5)/(asym_size*asym_size);
		if (asym_random > 0) {
			for (int i = 0; i < asym_kernel.length; i++) asym_kernel[i] = asym_random*(rand.nextDouble()-0.5)/(asym_size*asym_size);
		}
		asym_kernel[asym_size*i0+j0] = 1.0;
		if (asym_random < 0) {
			DoubleGaussianBlur gb = new DoubleGaussianBlur();
	    	gb.blurDouble(asym_kernel, asym_size, asym_size, -asym_random, -asym_random, 0.01);
		}
		
		
		asym_mask  [asym_size*i0+j0] = true;

		if (seed_size > 0 ){ // 0 - just center, for compatibility with the old code
			// the following code assumes that the asym_kernel can accommodate all initial cells
			if ((seed_size & 1) == 1) { // around single center pixel
				for (int n = 1; n <= (seed_size-1)/4; n++){
					asym_mask [asym_size * (i0 + n) + (j0 + n)] = true;
					asym_mask [asym_size * (i0 + n) + (j0 - n)] = true;
					asym_mask [asym_size * (i0 - n) + (j0 + n)] = true;
					asym_mask [asym_size * (i0 - n) + (j0 - n)] = true;
				}
			} else {
				int j00 = (xc < j0) ? (j0-1):j0;
				int i00 = (yc < i0) ? (i0-1):i0;
				for (int n = 0; n < seed_size/4; n++){
					asym_mask [asym_size * (i00 + n +1) + (j00 + n +1)] = true;
					asym_mask [asym_size * (i00 + n +1) + (j00 - n +0)] = true;
					asym_mask [asym_size * (i00 - n +0) + (j00 + n +1)] = true;
					asym_mask [asym_size * (i00 - n +0) + (j00 - n +0)] = true;
				}
			}
		}
		if (debugLevel>2){
			System.out.println("setInitialVector(target_kernel,"+seed_size+"): asym_mask.length="+asym_mask.length);
			System.out.println("asym mask: ");
			for (int ii=0;ii < asym_size;ii++){
				System.out.print(ii+": ");
				for (int jj=0;jj < asym_size;jj++){
					System.out.print((asym_mask[ii * asym_size + jj]?" X":" .")+" ");
				}
				System.out.println();	
			}
		}
		
		center_x = xc; // not limited by asym_size
		center_y = yc; // not limited by asym_size
		center_j0 = j0; // center to calcualte compatess odf asymmetrical kernel
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
				if ((debugLevel > 2) && (indx == 0)) {
					if (debugLevel>0)System.out.println("1.sym_kernel[0] = "+sym_kernel[0]+" i="+i+" j="+j+" i0="+i0+" j0="+j0+" conv_size="+conv_size+
				" target_kernel["+(conv_size*(i+i0) + (j+j0))+"] = " + target_kernel[(conv_size*(i+i0) + (j+j0))]);
				}
				sym_kernel_count[indx]++;
			}
			
		}
		if (debugLevel > 2)System.out.println("sym_kernel[0] = "+sym_kernel[0]+" sym_kernel_count[0] = " +sym_kernel_count[0]);
		for (int i = 0; i < sym_kernel.length; i++){
			if (sym_kernel_count[i] >0) sym_kernel[i] /= sym_kernel_count[i];
			else sym_kernel[i] = 0.0;
		}
		if (debugLevel > 2)System.out.println("sym_kernel[0] = "+sym_kernel[0]);
		double [][] kernels = {sym_kernel, asym_kernel};
		return kernels;
//		return setVectorFromKernels(sym_kernel, asym_kernel);
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
	
*/	

			

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
		double [] fX = getFX( kvect);
		kvect_inc[indx] += delta;
		double [] fx1 = getFX( kvect_inc);
//		if (indx == 63) {
//			System.out.println("----- getDerivDelta(): indx="+indx+" delta="+delta+" kvect["+indx+"]="+kvect[indx]+" kvect_inc["+indx+"]="+kvect_inc[indx]);
//		}
		for (int i = 0; i < fX.length; i++){
			fX[i] = (fx1[i]-fX[i])/delta;
		}
		return fX;
	}
	
	public double [] compareDerivative(
			int indx,
			double delta,          // value to increment parameter by for derivative calculation
			boolean verbose){
		double [] rslt = {0.0,0.0};
		double [][] jacob = getJacobian(this.currentVector);
		double [] deriv = getDerivDelta(this.currentVector, indx, delta);
		double [] fX = getFX(this.currentVector);
		for (int i = 0; i< fX.length; i++) {
			if (verbose) {
				System.out.println(i+": "+(jacob[indx][i]-deriv[i])+ " jacob["+indx+"]["+i+"] = "+jacob[indx][i]+
						" deriv["+i+"]="+deriv[i]+" f["+i+"]="+fX[i]);
			}
			rslt[0]+=(jacob[indx][i]-deriv[i])*(jacob[indx][i]-deriv[i]);
			rslt[1]+=jacob[indx][i]*jacob[indx][i];
		}
		rslt[0] = Math.sqrt(rslt[0]/fX.length);
		rslt[1] = Math.sqrt(rslt[1]/fX.length);
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

//		double [] fX = new double [conv_size*conv_size+asym_size*asym_size];
		double [] fX = new double [asym_terms_start + num_asym];

		if (this.debugLevel>3){
			System.out.println("fX(): vector_length=    "+kvect.length);
			System.out.println("fX(): sym_radius=       "+sym_radius);
			System.out.println("fX(): asym_size=        "+asym_size);
			System.out.println("fX(): conv_size=        "+conv_size);
			System.out.println("fX(): fX.length= "+fX.length);
			System.out.println("fX(): asym_start=       "+asym_start);
		}
		for (int i = 0; i < fX.length; i++) fX[i] = 0.0;

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
								fX[conv_index] += sd* kvect[asym_index];

//								if ((cind[asym_index]==63) && (conv_index==74)) {
//									System.out.println("cind["+asym_index+"]="+cind[asym_index]+
//											" conv_index ="+conv_index+" kvect["+asym_index+"]=" +kvect[asym_index]+
//											" fX["+conv_index+"]="+fX[conv_index]); 
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
					fX[asym_term++] = ir2 * cw* kvect[asym_start + asym_index];
				}
			}
		}
		return fX;
	}
	
	
	private double getCompactWeight(){
		return compactness_weight*sym_kernel_scale; // (asym_size*asym_size*asym_size*asym_size); // use
	}

	public int getNumPars(double [] kvect){
		int num_pars = 0;
		for (int i = 0; i<kvect.length ; i++){
			if (!Double.isNaN(kvect[i]))num_pars++;
		}
		return num_pars;
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
//								if ((debugLevel > 1) && (indx == 0)) {
//									System.out.println("                  getJacobian: indx="+indx+" asym_index="+asym_index+" cind[asym_index]="+cind[asym_index]+
//											" conv_index="+conv_index+" kvect["+asym_index+"], sd="+sd);
//									System.out.println("jacob["+indx+"]["+conv_index+"] += kvect["+asym_index+"] +="+kvect[asym_index]+" = "+ jacob[indx][conv_index]);
//								}
								
								jacob[cind[asym_index]][conv_index] += sd;
//								if ((debugLevel > 4) || (cind[asym_index] == 0)) {
//									System.out.println("+++++++++++++++++ getJacobian: indx="+indx+" asym_index="+asym_index+" cind[asym_index]="+cind[asym_index]+
//											" conv_index="+conv_index+" kvect["+asym_index+"], sd="+sd+ " jacob["+cind[asym_index]+"]["+conv_index+"]="+jacob[cind[asym_index]][conv_index]);
//								}
							}
						}
					}
				}
			}
		}
		
//		System.out.println("3:jacobian[0][167]="+jacob[0][167]);
//		System.out.println("3:jacobian[0][168]="+jacob[0][168]);
//		System.out.println("3:jacobian[0][169]="+jacob[0][169]);
//		System.out.println("3:jacobian[0][170]="+jacob[0][170]);
		
		
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
		
//		System.out.println("4:jacobian[0][167]="+jacob[0][167]);
//		System.out.println("4:jacobian[0][168]="+jacob[0][168]);
//		System.out.println("4:jacobian[0][169]="+jacob[0][169]);
//		System.out.println("4:jacobian[0][170]="+jacob[0][170]);
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
			double [] fX)            // current convolution result of async_kernel (*) sync_kernel, extended by asym_kernel components
			{
		double [] jTByDiff = new double [jacob.length];
		for (int i=0; i < jTByDiff.length; i++){
			jTByDiff[i] = 0;
			for (int k = 0; k< target_kernel.length; k++){
				jTByDiff[i] += jacob[i][k]*(target_kernel[k]-fX[k]);
			}
			for (int k = target_kernel.length; k< fX.length; k++){
				jTByDiff[i] += jacob[i][k]*(-fX[k]);
			}

		}

		return jTByDiff;
	}
	private double[] getDiffByDiff(
			double [] target_kernel, // target kernel
			double [] fX)            // current convolution result of async_kernel (*) sync_kernel, extended async kernel components
			{
		double [] diffByDiff = {0.0,0.0};
		for (int k=0; k< target_kernel.length; k++){
			double d = target_kernel[k]-fX[k];
			diffByDiff[0] += d*d;
		}
		diffByDiff[1] = diffByDiff[0]; // actual squared error, without compactness 
		for (int k=target_kernel.length; k< fX.length; k++){
			double d = fX[k];
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

	public void resetLMARuns(){
		numLMARuns = 0;
	}
	public int getLMARuns(){
		return numLMARuns;
	}
	

	private void initLevenbergMarquardt(double fact_precision, int seed_size, double asym_random){
		lMAData = new LMAData(debugLevel);
		lMAData.setTarget(target_kernel);
		lMAData.setTargetWindowMode(target_window_mode, centerWindowToTarget);		
    	double [][]kernels = setInitialVector(target_kernel, seed_size, asym_random); // should be (asym_size + 2*sym_radius-1)**2
    	sym_kernel_scale = kernels[0][0];
    	for (int i=0;i<kernels[0].length;i++) if (!Double.isNaN(kernels[0][i])) kernels[0][i] /=sym_kernel_scale;
    	for (int i=0;i<kernels[1].length;i++) if (!Double.isNaN(kernels[1][i])) kernels[1][i] *=sym_kernel_scale;
    	lMAData.setSymKernel (kernels[0]);
    	lMAData.setAsymKernel(kernels[1]);
		this.goal_rms_pure = lMAData.getTargetRMSW()*fact_precision;
		this.target_rms = lMAData.getTargetRMSW(); //Math.sqrt(s/target_kernel.length);
    	resetLMARuns();
	}
	
	private boolean levenbergMarquardt(){
    	long startTime=System.nanoTime();
    	this.iterationStepNumber = 0;
    	this.lambda = this.init_lambda;
    	lMAData.resetRMSes(); // both first and saved
		lastImprovements[0]=-1.0;
		lastImprovements[1]=-1.0;
    	lMAData.invalidateLMAArrays(); // should be in each run , not at init
    	
    	if (this.numIterations < 0){
			lMAData.getFX(true); // try false too			
            return true;    		
    	}
    	
		if (this.debugLevel > 3) {
			System.out.println("this.currentVector 0");
			double [] dbg_vector = lMAData.getVector();
			for (int i=63;i<dbg_vector.length;i++){
				System.out.println(i+": "+ dbg_vector[i]);
			}
		}
		if (this.debugLevel>2){
			System.out.println("LMA before loop, RMS = "+lMAData.getRMSes()[0]+", pure="+lMAData.getRMSes()[1]);
		}
    	while (true) { // loop for the same series
    		boolean [] state=stepLevenbergMarquardt(goal_rms_pure);
    		if ((this.iterationStepNumber==1) && (this.debugLevel > 2)){
    			System.out.println("LMA after first stepLevenbergMarquardt, RMS = "+lMAData.getRMSes()[0]+", pure="+lMAData.getRMSes()[1]);
    		}

    		
    		// state[0] - better, state[1] - finished
    		if (this.debugLevel > 2) System.out.println(this.iterationStepNumber+": stepLevenbergMarquardt()==>"+state[1]+":"+state[0]);
    		//    		boolean cont=true;
    		// Make it success if this.currentRMS<this.firstRMS even if LMA failed to converge
    		if (state[1] && !state[0] && (lMAData.getFirstRMSes()[0] > lMAData.getRMSes()[0])) {
    			if (this.debugLevel > 2) System.out.println("LMA failed to converge, but RMS improved from the initial value ("+lMAData.getRMSes()[0]+" < "+lMAData.getFirstRMSes()[0]+")");
    			state[0]=true;
    		}
    		if (this.debugLevel > 1){
    			if (state[1] && !state[0]){ // failure, show at debugLevel >0
    				System.out.println("LevenbergMarquardt(): failed step ="+this.iterationStepNumber+
    						", RMS="+IJ.d2s(lMAData.getRMSes()[0], 8)+
    						" ("+IJ.d2s(lMAData.getFirstRMSes()[0], 8)+") "+
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()- startTime), 3));
    			} else if (this.debugLevel > 2){ // success: show only if debugLevel > 1
    				System.out.println("==> LevenbergMarquardt(): before action step ="+this.iterationStepNumber+
    						", RMS="+IJ.d2s(lMAData.getRMSes()[0], 8)+
    						" ("+IJ.d2s(lMAData.getFirstRMSes()[0], 8)+") "+
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
    			}
    		}
//    		stepLevenbergMarquardtAction(startTime); // apply step - in any case?
    		this.iterationStepNumber++;
    		if (this.debugLevel > 2) {
    			System.out.println(
    					"stepLevenbergMarquardtAction() step="+this.iterationStepNumber+
    					", this.currentRMS="+lMAData.getSavedRMSes()[0]+
    					", this.currentRMSPure="+lMAData.getSavedRMSes()[1]+
    					", this.nextRMS="+lMAData.getRMSes()[0]+
    					", this.nextRMSPure="+lMAData.getRMSes()[1]+
    					" lambda="+this.lambda+" at "+IJ.d2s(0.000000001*(System.nanoTime()- startTime),3)+" sec");
    		}
//    		if (this.nextRMS   <  this.currentRMS) { //improved
    		if (state[0]) { // improved
    			this.lambda    *= this.lambdaStepDown;
    			this.currentRMS = this.nextRMS;
    		} else {
    			this.lambda*=     this.lambdaStepUp;
    		}
    		
    		if ((this.debugLevel > 1) && ((this.debugLevel > 3) || ((System.nanoTime()-this.startTime)>30000000000.0))){ // > 10 sec
//       		if ((this.debugLevel>0) && ((this.debugLevel>0) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
    			System.out.println("--> long wait: LevenbergMarquardt(): step = "+this.iterationStepNumber+
    					", RMS="+IJ.d2s(lMAData.getRMSes()[0],8)+
    					" ("+IJ.d2s(lMAData.getFirstRMSes()[0],8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
    		}
    		if (state[1]) { // finished
    			if (!state[0]) {
    		    	numLMARuns++;
    				return false; // sequence failed
    			}
    			break; // while (true), proceed to the next series
    		}
    	}
    	if (this.debugLevel > 1) System.out.println("LevenbergMarquardt() finished in "+this.iterationStepNumber+
    			" steps, RMS="+lMAData.getRMSes()[0]+
    			" ("+lMAData.getFirstRMSes()[0]+") "+
    			", RMSpure="+lMAData.getRMSes()[1]+
    			" ("+lMAData.getFirstRMSes()[1]+") "+
    			", Relative RMS="+ (lMAData.getRMSes()[1]/this.target_rms)+
    			" at "+ IJ.d2s(0.000000001*(System.nanoTime()- startTime),3));
    	numLMARuns++;
    	return true; // all series done
	}
	
	private boolean [] stepLevenbergMarquardt(double goal_rms_pure){
		double [] deltas=null;
//		System.out.println("lMAData.isValidLMAArrays()="+lMAData.isValidLMAArrays());
		if (!lMAData.isValidLMAArrays()) { // First time and after improvement
			lMAData.getFX(true); // Will calculate fX and rms (both composite and pure) if null (only when firs called) (try with false too?)
			lMAData.getJacobian(
					true, // boolean recalculate,
					true); //boolean skip_disabled_asym)
			lMAData.getJTByJW(true);    // force recalculate (it is null anyway)
			lMAData.getJTByDiffW(true); // force recalculate
			if (debugLevel >2) {
				double [][] jTByJ =   lMAData.getJTByJW(false);    // do not recalculate
				double [] jTByDiff =  lMAData.getJTByDiffW(false); // do not recalculate
				int    [][] map_from_pars = lMAData.getMapFromPars();
				
				for (int n=0; n < jTByJ.length;n++) if ((debugLevel > 3) || (map_from_pars[n][0]==1)){
					for (int i=0; i < jTByJ.length; i++){
						System.out.println("jTByJ["+n+"]["+i+"]="+jTByJ[n][i]);
					}
					System.out.println("jTByDiff["+n+"]="+jTByDiff[n]);
				}
			}
			if (this.debugLevel > 2) {
				System.out.println("initial RMS="+IJ.d2s(lMAData.getRMSes()[0],8)+
						" ("+IJ.d2s(lMAData.getRMSes()[1],8)+")"+
						". Calculating next Jacobian. Points:"+lMAData.getNumPoints()+"("+lMAData.getNumPurePoints()+")"+
						" Parameters:"+lMAData.getNumPars());
			}
			lMAData.save(); // save kernels (vector), fX and RMS values
			
		} else { // LMA arrays already calculated after previous failure
			if (debugLevel > 3){
				System.out.println("existing data: this.currentRMS="+IJ.d2s(lMAData.getRMSes()[0], 8)+
						" ("+IJ.d2s(lMAData.getRMSes()[1], 8)+")");
			}
		}

		// calculate deltas
		deltas=lMAData.solveLMA(this.lambda, this.debugLevel);
		
		if (deltas==null) {
			if (this.debugLevel > 0) {
				System.out.println("--- Singular matrix - failed to compute deltas ---");
			}
			deltas=new double[lMAData.getNumPars()];
			for (int i=0;i<deltas.length;i++) deltas[i]=0.0;
			boolean [] status={false, true}; // done / bad
			return status; // done, bad . No need to restore - nothing was changed
		}
		if (this.debugLevel > 3 ){
			System.out.println("--- deltas ---"+" this.currentRMS="+lMAData.getRMSes()[0]);
			int    [][] map_from_pars = lMAData.getMapFromPars();
			for (int i=0; i < deltas.length; i++)  if ((debugLevel > 3) || (map_from_pars[i][0]==1)) {
				System.out.println("deltas["+i+"]="+ deltas[i]);
			}
		}

		// apply deltas
		lMAData.applyDeltas(deltas); 
		if (this.debugLevel > 3) {
			System.out.println("modified Vector");
			double [] dbg_vector = lMAData.getVector();
			int    [][] map_from_pars = lMAData.getMapFromPars();
			for (int i= 0; i<dbg_vector.length;i++) if ((debugLevel > 3) || (map_from_pars[i][0]==1)){
//			for (int i= 0; i<dbg_vector.length;i++) if ((debugLevel > 1) || (map_from_pars[i][0]==1)){
				System.out.println(i+": "+ dbg_vector[i]+" ("+deltas[i]+")");
			}
		}
		
		// calculate for new (modified vector) 
		lMAData.getFX(true); // also calculates RMSes  (try false too)

		this.lastImprovements[1]=this.lastImprovements[0];
		this.lastImprovements[0]=lMAData.getSavedRMSes()[0] - lMAData.getRMSes()[0]; // this.currentRMS-this.nextRMS;
		if (this.debugLevel > 3) {
			System.out.println("stepLMA currentRMS="+lMAData.getSavedRMSes()[0]+
					", nextRMS="+lMAData.getRMSes()[0]+
					", delta="+(lMAData.getSavedRMSes()[0]-lMAData.getRMSes()[0]));
		}
		boolean [] status={lMAData.getSavedRMSes()[0] > lMAData.getRMSes()[0], false};
		// additional test if "worse" but the difference is too small, it was be caused by computation error, like here:
		//stepLevenbergMarquardtAction() step=27, this.currentRMS=0.17068403807026408,   this.nextRMS=0.1706840380702647

		if (!status[0]) { // worse
			if (lMAData.getRMSes()[0] < (lMAData.getSavedRMSes()[0] * (1.0+this.thresholdFinish*0.01))) {
//				this.nextRMS=this.currentRMS; // no need to modify ?
				status[0]=true;
				status[1]=true;
				this.lastImprovements[0]=0.0;
				if (this.debugLevel > 2) {
					System.out.println("New RMS error is larger than the old one, but the difference is too small to be trusted ");
					System.out.println(
							"stepLMA this.currentRMS="+lMAData.getSavedRMSes()[0]+
							", this.currentRMSPure="+lMAData.getSavedRMSes()[1]+
							", this.nextRMS="+lMAData.getRMSes()[0]+
							", this.nextRMSPure="+lMAData.getRMSes()[1]+
							", delta="+(lMAData.getSavedRMSes()[0] - lMAData.getRMSes()[0])+
							", deltaPure="+(lMAData.getSavedRMSes()[1] - lMAData.getRMSes()[1]));
				}
			}
		}
		if (status[0]) { // improved
			status[1]=(this.iterationStepNumber>this.numIterations) || ( // done
					(this.lastImprovements[0]>=0.0) &&
					(this.lastImprovements[0]<this.thresholdFinish*lMAData.getSavedRMSes()[0]) &&
					(this.lastImprovements[1]>=0.0) &&
					(this.lastImprovements[1]<this.thresholdFinish*lMAData.getSavedRMSes()[0]));
			if ( lMAData.getRMSes()[1] < goal_rms_pure) {
				status[1] = true;
				if (this.debugLevel > 2) {
					System.out.println("Improvent is possible, but the factorization precision reached its goal");
					System.out.println(
							"stepLMA this.currentRMS="+lMAData.getSavedRMSes()[0]+
							", this.currentRMSPure="+lMAData.getSavedRMSes()[1]+
							", this.nextRMS="+lMAData.getRMSes()[0]+
							", this.nextRMSPure="+lMAData.getRMSes()[1]+
							", delta="+(lMAData.getSavedRMSes()[0]-lMAData.getRMSes()[0])+
							", deltaPure="+(lMAData.getSavedRMSes()[1]-lMAData.getRMSes()[1]));
				}				
			}
			lMAData.invalidateLMAArrays();
		} else { // did not improve - roll back
			lMAData.restore();       // roll back: kernels, fX, RMSes
			status[1]=(this.iterationStepNumber>this.numIterations) || // failed
					((this.lambda*this.lambdaStepUp)>this.maxLambda);
		}
		///this.currentRMS    	
		//TODO: add other failures leading to result failure?    	
		if (this.debugLevel > 3) {
			System.out.println("stepLevenbergMarquardt()=>"+status[0]+","+status[1]);
		}
		return status;
	}
	

	
	

	
// Old version
	
	private void initLevenbergMarquardt_old(double fact_precision, int seed_size, double asym_random){
		double s = 0.0;
		for (int i = 0; i<target_kernel.length; i++){
			s+= target_kernel[i]*target_kernel[i];
		}
		this.goal_rms_pure = Math.sqrt(s/target_kernel.length)*fact_precision;
//    	this.currentVector = setInitialVector(target_kernel, 1); // should be (asym_size + 2*sym_radius-1)**2
    	double [][]kernels = setInitialVector(target_kernel, seed_size, asym_random); // should be (asym_size + 2*sym_radius-1)**2
    	this.currentVector = setVectorFromKernels(kernels[0], kernels[1]);
    	resetLMARuns();
    	
	}
	
	private boolean levenbergMarquardt_old(){
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
    	
		if (this.debugLevel > 3) {
			System.out.println("this.currentVector 0");
			for (int i=63;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}

    	while (true) { // loop for the same series
    		boolean [] state=stepLevenbergMarquardtFirst_old(goal_rms_pure);
    		if (this.debugLevel > 2) System.out.println(this.iterationStepNumber+": stepLevenbergMarquardtFirst_old()==>"+state[1]+":"+state[0]);
    		//    		boolean cont=true;
    		// Make it success if this.currentRMS<this.firstRMS even if LMA failed to converge
    		if (state[1] && !state[0] && (this.firstRMS>this.currentRMS)){
    			if (this.debugLevel>1) System.out.println("LMA failed to converge, but RMS improved from the initial value ("+this.currentRMS+" < "+this.firstRMS+")");
    			state[0]=true;
    		}
    		if ((this.stopOnFailure && state[1] && !state[0])){
    			if (this.debugLevel > 1){
    				System.out.println("LevenbergMarquardt(): step ="+this.iterationStepNumber+
    						", RMS="+IJ.d2s(this.currentRMS,8)+
    						" ("+IJ.d2s(this.firstRMS,8)+") "+
    						") at "+ IJ.d2s(0.000000001*(System.nanoTime()- startTime),3));
    			}
    			long startDialogTime=System.nanoTime();
    			startTime+=(System.nanoTime()-startDialogTime); // do not count time used by the User.
    			//    			if (this.showThisImages) showDiff (this.currentfX, "fit-"+this.iterationStepNumber);
    			//    			if (this.showNextImages) showDiff (this.nextfX,    "fit-"+(this.iterationStepNumber+1));
    		} else if (this.debugLevel > 2){
    			System.out.println("==> LevenbergMarquardt(): before action step ="+this.iterationStepNumber+
    					", RMS="+IJ.d2s(this.currentRMS,8)+
    					" ("+IJ.d2s(this.firstRMS,8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
    		}
    		stepLevenbergMarquardtAction_old(startTime); // apply step - in any case?
    		if ((this.debugLevel > 1) && ((this.debugLevel > 2) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
//       		if ((this.debugLevel>0) && ((this.debugLevel>0) || ((System.nanoTime()-this.startTime)>10000000000.0))){ // > 10 sec
    			System.out.println("--> LevenbergMarquardt(): step = "+this.iterationStepNumber+
    					", RMS="+IJ.d2s(this.currentRMS,8)+
    					" ("+IJ.d2s(this.firstRMS,8)+") "+
    					") at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
    		}
    		//stepLevenbergMarquardtAction_old();    			
    		if (state[1]) {
    			if (!state[0]) {
    		    	numLMARuns++;
    				return false; // sequence failed
    			}
    			break; // while (true), proceed to the next series
    		}
    	}
    	//    		if (this.fittingStrategy.isLastSeries(this.seriesNumber)) break;
    	if (this.debugLevel > 1) System.out.println(
    			"LevenbergMarquardt_old() finished in "+this.iterationStepNumber+
    			" steps, RMS="+this.currentRMS+
    			" ("+this.firstRMS+")"+
    			", RMSPure="+this.currentRMSPure+
    			" ("+this.firstRMSPure+") "+
    			"at "+ IJ.d2s(0.000000001*(System.nanoTime()- startTime),3));
    	if (this.debugLevel > 3) {
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
    	numLMARuns++;
    	return true; // all series done
	}
	
	private boolean [] stepLevenbergMarquardtFirst_old(double goal_rms_pure){
		int deltas_indx;
		double [] deltas=null;
		double [] rmses; // [0]: full rms, [1]:pure rms
		// moved to caller
		// calculate  this.currentfX, this.jacobian if needed
		if (this.debugLevel > 4) {
			System.out.println("this.currentVector 1");
			for (int i=0;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		if (this.debugLevel > 3) {
			System.out.println("this.currentVector 2");
			for (int i=63;i<this.currentVector.length;i++){
				System.out.println(i+": "+ this.currentVector[i]);
			}
		}
		if ((this.currentfX==null)|| (this.lMAArrays==null)) {
			
			this.currentfX=getFX  (this.currentVector);			
			this.jacobian = getJacobian (this.currentVector);
			
			if (debugLevel > 2){
				int conv_size =       asym_size + 2*sym_radius-2;
				int asym_terms_start= conv_size*conv_size;
				for (int i=asym_terms_start; i<currentfX.length; i++){
					System.out.println("this.currentfX["+i+"]="+this.currentfX[i]);
				}
				for (int n=(debugLevel > 3)?0:63; n<this.jacobian.length;n++){
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
			if (debugLevel > 3) {
				System.out.println("Calculated new lMAArrays, lMAArrays.jTByJ[0][0] = "+lMAArrays.jTByJ[0][0]);
			}
			if (debugLevel > 3) {
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
			this.currentRMS =     Math.sqrt(rmses[0] / currentfX.length);
			if (debugLevel > 2){
				System.out.println("currentRMSPure= "+ currentRMSPure + " getDiffByDiff[1] = "+Math.sqrt(getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX)[1] / target_kernel.length));
			}
			
			if (this.debugLevel > 3) {
				System.out.println("initial RMS="+IJ.d2s(this.currentRMS,8)+
						" ("+IJ.d2s(this.currentRMSPure,8)+")"+
//						". Calculating next Jacobian. Points:"+(asym_size*asym_size+target_kernel.length)+" Parameters:"+this.currentVector.length);
				". Calculating next Jacobian. Points:"+currentfX.length+" Parameters:"+getNumPars(this.currentVector));
				
			}
		} else {
			rmses = getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX);
			this.currentRMSPure=  Math.sqrt(rmses[1] / target_kernel.length);
			this.currentRMS =     Math.sqrt(rmses[0] / currentfX.length);
			
			if (debugLevel > 3){
				System.out.println("this.currentRMS=" + this.currentRMS+ " getDiffByDiff[1] = "+Math.sqrt(getDiffByDiff(
					this.target_kernel, // target kernel
					this.currentfX)[1] / target_kernel.length));
			}
			if (debugLevel > 2) {
				System.out.println("Reused existing lMAArrays (after previous failure), lMAArrays.jTByJ[0][0] = "+lMAArrays.jTByJ[0][0]);
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
		if (this.debugLevel > 3) {
			System.out.println("--- deltas.1 ---");
			for (int i=63;i<deltas.length;i++){
				System.out.println("deltas["+i+"]="+ deltas[i]+" this.currentRMS="+this.currentRMS);
			}
		}

		if (this.debugLevel > 3){ // 2!!! ) {
			System.out.println("--- deltas ---"+" this.currentRMS="+this.currentRMS);
//			for (int i=0; i < deltas.length; i++)  if ((debugLevel > 3) || (map_from_pars[i][0]==1)) {
			for (int i=0; i < deltas.length; i++)  if ((this.debugLevel > 4) || (i >= 63)) {
				System.out.println("deltas["+i+"]="+ deltas[i]);
			}
		}

		// apply deltas    	
		this.nextVector=this.currentVector.clone();
		
		deltas_indx = 0;
		for (int i=0;i<this.nextVector.length;i++) {
			if (!Double.isNaN(this.nextVector[i])){
				this.nextVector[i]+=deltas[deltas_indx++];
			}
		}
		
		// another option - do not calculate J now, just fX. and late - calculate both if it was improvement    	
		//    	save current Jacobian
		if (this.debugLevel > 3) { // change to 2 later
			int indx = 0;
			System.out.println("modified Vector");
			for (int i=0;i<this.nextVector.length;i++) if (!Double.isNaN(this.nextVector[i])){
				System.out.println(indx+": "+ this.nextVector[i]+" ("+deltas[indx]+")");
				indx++;
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
		this.nextRMS =     Math.sqrt(rmses[0] / nextfX.length);

		if (debugLevel > 3){
			System.out.println("nextRMSPure= "+ nextRMSPure + " target_kernel.length = "+target_kernel.length+" getDiffByDiff[1] = "+Math.sqrt(getDiffByDiff(
				this.target_kernel, // target kernel
				this.nextfX)[1] / target_kernel.length));     // current convolution result of async_kernel (*) sync_kernel
		}
		this.lastImprovements[1]=this.lastImprovements[0];
		this.lastImprovements[0]=this.currentRMS-this.nextRMS;
//		System.out.println("Setting this.lastImprovements[1]="+this.lastImprovements[1]+" new this.lastImprovements[0]="+this.lastImprovements[0]);
		if (this.debugLevel > 3) {
			System.out.println("stepLMA this.currentRMS="+this.currentRMS+
					", this.nextRMS="+this.nextRMS+
					", delta="+(this.currentRMS-this.nextRMS));
		}
		boolean [] status={matrixNonSingular && (this.nextRMS<=this.currentRMS),!matrixNonSingular};
		// additional test if "worse" but the difference is too small, it was be caused by computation error, like here:
		//stepLevenbergMarquardtAction_old() step=27, this.currentRMS=0.17068403807026408,   this.nextRMS=0.1706840380702647

		if (!status[0] && matrixNonSingular) {
			if (this.nextRMS<(this.currentRMS+this.currentRMS*this.thresholdFinish*0.01)) {
				this.nextRMS=this.currentRMS;
				status[0]=true;
				status[1]=true;
				this.lastImprovements[0]=0.0;
				if (this.debugLevel > 2) {
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
				if (this.debugLevel > 2) {
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
			this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff - it is done in stepLevenbergMarquardtAction_old

			status[1]=(this.iterationStepNumber>this.numIterations) || // failed
					((this.lambda*this.lambdaStepUp)>this.maxLambda);
		}
		///this.currentRMS    	
		//TODO: add other failures leading to result failure?    	
		if (this.debugLevel > 3) {
			System.out.println("stepLevenbergMarquardt()=>"+status[0]+","+status[1]);
		}
		return status;
	}
	
    /**
     * Apply fitting step 
     */
	private void stepLevenbergMarquardtAction_old(long startTime){// 
    	this.iterationStepNumber++;
// apply/revert,modify lambda    	
		if (this.debugLevel>1) {
			System.out.println(
					"stepLevenbergMarquardtAction_old() step="+this.iterationStepNumber+
					", this.currentRMS="+this.currentRMS+
					", this.currentRMSPure="+this.currentRMSPure+
					", this.nextRMS="+this.nextRMS+
					", this.nextRMSPure="+this.nextRMSPure+
					" lambda="+this.lambda+" at "+IJ.d2s(0.000000001*(System.nanoTime()- startTime),3)+" sec");
		}
    	if (this.nextRMS<this.currentRMS) { //improved
    		if (this.debugLevel > 2) {
    			System.out.println("Using new fX and vector");
    		}
    		this.lambda*=this.lambdaStepDown;
    		this.currentRMS=this.nextRMS;
    		this.currentfX=this.nextfX;
    		this.currentVector=this.nextVector;
    		this.lMAArrays = null; // need to calculate new ones
    	} else {
    		this.lambda*=this.lambdaStepUp;
    		this.lMAArrays=this.savedLMAArrays; // restore Jt*J and Jt*diff
    		if (this.debugLevel > 2) {
    			System.out.println("Re-using saved saved LMAArrays");
    		}
    	}
    }

}
