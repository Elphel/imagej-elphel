package com.elphel.imagej.dp;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Rectangle;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;

import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.showDoubleFloatArrays;
import com.elphel.imagej.jp4.JP46_Reader_camera;

public class EyesisAberrations {
	public double [][][][] pdfKernelMap=null;
	JP46_Reader_camera JP4_INSTANCE=       new JP46_Reader_camera(false);
	showDoubleFloatArrays SDFA_INSTANCE=   new showDoubleFloatArrays();
    public AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
	public Distortions distortions=null;
	public AberrationParameters aberrationParameters=null;
    
    public EyesisAberrations (AtomicInteger stopRequested,
    		AberrationParameters aberrationParameters){
    	this.stopRequested=stopRequested;
    	this.aberrationParameters=aberrationParameters;
    }
   	public void setDistortions(Distortions distortions){
		this.distortions=distortions;
	}
   	
   	int countExistentFiles(String directory,String [] paths, boolean remove){
   		int numFiles=0;
   		for (int i=0;i<paths.length;i++) if (paths[i]!=null){
   			String path= ( (directory!=null)?(directory+Prefs.getFileSeparator()):"")+paths[i];
   			if ((new File(path)).exists()) {
   				numFiles++;
   			} else if (remove) {
   				paths[i]=null;
   			}
   		}
   		return numFiles;
   	}

  	public boolean reverseKernels(
		    AtomicInteger stopRequested, // 1 - stop now, 2 - when convenient
  			EyesisAberrations.InverseParameters inverseParameters, // size (side of square) of direct PSF kernel
   			boolean                saveResult,
   			boolean                showResult,
   			boolean                updateStatus,          // UPDATE_STATUS
   			int                    threadsMax,
   			int                    globalDebugLevel
   	){
   		if ((this.aberrationParameters.aberrationsKernelDirectory==null) || (this.aberrationParameters.aberrationsKernelDirectory.length()==0)){
   			if (aberrationParameters.selectAberrationsKernelDirectory(true, this.aberrationParameters.aberrationsKernelDirectory, false)==null) {
   				String msg = "Nothing selected";
   				System.out.println("Warning"+msg);
   				IJ.showMessage("Warning",msg);
   				return false;
   			}
   		}
   		int numChannels=distortions.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
    	boolean [] selectedChannels=this.aberrationParameters.getChannelSelection(distortions);
		String [] srcPaths=    new String[numChannels];
		String [] resultPaths= new String[numChannels];
		int numToProcess=0;
    	for (int nChn=0;nChn<selectedChannels.length;nChn++){
    		if (!selectedChannels[nChn]){
    			srcPaths[nChn]=null;
    			resultPaths[nChn]=null;
    		} else {
    			srcPaths[nChn]=this.aberrationParameters.psfKernelDirectory+Prefs.getFileSeparator()+
       			this.aberrationParameters.interpolatedPSFPrefix+String.format("%02d", nChn)+
       			this.aberrationParameters.interpolatedPSFSuffix;
    			resultPaths[nChn]=this.aberrationParameters.aberrationsKernelDirectory+Prefs.getFileSeparator()+
       			this.aberrationParameters.aberrationsPrefix+String.format("%02d", nChn)+
       			this.aberrationParameters.aberrationsSuffix;
    			if (!this.aberrationParameters.overwriteResultFiles && (new File(resultPaths[nChn])).exists()) {
    				srcPaths[nChn]=null;
    				if (globalDebugLevel>0) System.out.println("File "+resultPaths[nChn]+" already exists and overwrite is disabled in configuration, channel "+nChn+" will be skipped");
    			}
    			numToProcess++;
    		}
    	}
    	if (numToProcess==0){
				String msg = "No kernels to process";
   				System.out.println("Warning"+msg);
   				IJ.showMessage("Warning",msg);
   				return false;
    	}
   		int numProcessed=0;
   		Opener opener=new Opener();
   		ImagePlus impSpsf;
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
   		for (int nChn=0;nChn<numChannels;nChn++) if (srcPaths[nChn]!=null){
   			if (!(new File(srcPaths[nChn])).exists()) {
   				String msg = "Interpolated PSF kernel stack for channel #"+nChn+": "+srcPaths[nChn]+" does not exist";
   				System.out.println("Warning"+msg);
   				continue;
   			}
   			impSpsf=opener.openImage("", srcPaths[nChn]);
   			if (impSpsf==null) {
   				System.out.println("Failed to open interpolated PSF kernel stack "+srcPaths[nChn]);
   				continue;
   			}
   			if (impSpsf.getStackSize()<3) {
   				System.out.println("Need a 3-layer stack with interpolated PSF kernels");
   				continue;
   			}
   			ImageStack stack= reversePSFKernelStack(
   					impSpsf.getStack(), // stack of 3 32-bit (float) images, made of square kernels
   					inverseParameters, // size (side of square) of direct PSF kernel
   					threadsMax, // size (side of square) of reverse PSF kernel
   					updateStatus,
   					globalDebugLevel); 

   			ImagePlus impInvertedPSF = new ImagePlus("interpolated kernel stack", stack);
   			JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
   			jp4_instance.decodeProperiesFromInfo(impSpsf);
   			// copy properties from the source image
   			jp4_instance.copyProperties (impSpsf,impInvertedPSF);
   			inverseParameters.setProperties("INVERSE.", impInvertedPSF);
   			jp4_instance.encodeProperiesToInfo(impInvertedPSF);
   			if (showResult) {
   				impInvertedPSF.getProcessor().resetMinAndMax(); // imp_psf will be reused
   				impInvertedPSF.show();
   			} 
   			if (saveResult){
   				if (globalDebugLevel>0) System.out.println((numProcessed+1)+" of "+numToProcess+": saving invered (of the file"+srcPaths[nChn]+") kernel to "+resultPaths[nChn]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
   				FileSaver fs=new FileSaver(impInvertedPSF);
   				fs.saveAsTiffStack(resultPaths[nChn]);
   			}
   			numProcessed++;
    		if 	(stopRequested.get()>0) {				
				if (globalDebugLevel>0) System.out.println("User requested stop");
				break;
    		}
   		}
   		if (numProcessed>0){
			if (globalDebugLevel>0) {
				System.out.println("Inverted "+numProcessed+" kernel stacks at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			}
   		} else {
				String msg = "No kernel files to invert";
   				System.out.println("Warning"+msg);
   				return false;
   		}
   		return true;
   	}
   	
	public ImageStack  reversePSFKernelStack(
			final ImageStack            PSFStack, // stack of 3 32-bit (float) images, made of square kernels
			final EyesisAberrations.InverseParameters inverseParameters, // size (side of square) of direct PSF kernel
			final int                 threadsMax, // size (side of square) of reverse PSF kernel
			final boolean           updateStatus,
			final int globalDebugLevel){  // update status info
		if (PSFStack==null) return null;
		final int tilesX=PSFStack.getWidth()/inverseParameters.dSize;
		final int tilesY=PSFStack.getHeight()/inverseParameters.dSize;
		final int nChn=PSFStack.getSize();
		final double [] sigmas={inverseParameters.blurIndividual,inverseParameters.blurIndividual,inverseParameters.blurChecker};
		final float [][] outPixels=new float[nChn][tilesX*inverseParameters.rSize*tilesY*inverseParameters.rSize];
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int numberOfKernels=     tilesY*tilesX*nChn;
		final int numberOfKernelsInChn=tilesY*tilesX;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					float [] pixels=null;
					double [] kernel= new double[inverseParameters.dSize*inverseParameters.dSize];
					double [] rKernel=new double[inverseParameters.rSize*inverseParameters.rSize];
					int  [][]selection;
					double [] ellipse_coeff;
					double [] variableSigmas;
					int chn,tileY,tileX;
					int chn0=-1;
					DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					for (int nTile = ai.getAndIncrement(); nTile < numberOfKernels; nTile = ai.getAndIncrement()) {
						chn=nTile/numberOfKernelsInChn;
						tileY =(nTile % numberOfKernelsInChn)/tilesX;
						tileX = nTile % tilesX;
						if (updateStatus) IJ.showStatus("Invertinging PSF, channel "+(chn+1)+" of "+nChn+", row "+(tileY+1)+" of "+tilesY);
						if (chn!=chn0) {
							pixels=(float[]) PSFStack.getPixels(chn+1);
							chn0=chn;
						}
						extractOneKernel( pixels, //  array of combined square kernels, each 
								kernel, // will be filled, should have correct size before call
								tilesX, // number of kernels in a row
								tileX, // horizontal number of kernel to extract
								tileY); // vertical number of kernel to extract
						/* Find direct kernel approximation ellipse, increase it, mirror center around 0,0 and use it as a mask for the reversed kernel */
						selection=    findClusterOnPSF(kernel,  inverseParameters.psfCutoffEnergy, "",globalDebugLevel);
						ellipse_coeff=findEllipseOnPSF(kernel,  selection, "",globalDebugLevel); // coefficients for direct PSF, for rPSF [0] and [1] need to be opposite size

						rKernel=resizeForFFT(kernel,inverseParameters.rSize);
/* Apply variable blur to direct kernel using it's center X,Y */
						if (inverseParameters.filterDirect) {
							variableSigmas= createSigmasFromCenter(inverseParameters.rSize, // side of square
									inverseParameters.sigmaToRadiusDirect, // variable blurring - sigma will be proportional distance from the center
									sigmas[chn]*inverseParameters.sigmaScaleDirect, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
									ellipse_coeff[0], // coordinates of the center (0:0 - size/2: size/2)
									ellipse_coeff[1]);
							rKernel=variableGaussBlurr(          rKernel, // input square pixel array, preferably having many exact zeros (they will be skipped)
									variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
									3.5, // drop calculatin if farther then nSigma
									0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
									0, // int WOICenterY, // 
									inverseParameters.rSize, //int WOIWidth, reduce later
									inverseParameters.rSize, //int WOIHeight)
									globalDebugLevel);
						}
						
/* reverse PSF kernel */
						rKernel= cleanupAndReversePSF (rKernel,  // input pixels
								inverseParameters,
								//    						  false,  // fold high frequency into low, when false - use Hamming to cut off high frequencies
								fht_instance,
						"",
						globalDebugLevel); // just for the plot names
/*  mask  the reversed kernel */
						rKernel= maskReversePSFKernel(rKernel, // reversed psf, square array
								ellipse_coeff, // ellipse coefficients from _direct_ kernel
								inverseParameters.psfEllipseScale,
								inverseParameters.rpsfMinMaskThreshold); // zero output element if elliptical Gauss mask is below this threshold

						normalizeKernel(rKernel); // in-place
/* Apply variable blur to inversed kernel, using (and reversing sign) the center X,Y from the direct kernel */
						if (inverseParameters.filter) {
							variableSigmas= createSigmasFromCenter(inverseParameters.rSize, // side of square
									inverseParameters.sigmaToRadius, // variable blurring - sigma will be proportional distance from the center
									sigmas[chn]*inverseParameters.sigmaScale, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
									-ellipse_coeff[0], // coordinates of the center (0:0 - size/2: size/2)
									-ellipse_coeff[1]);
							rKernel=variableGaussBlurr(          rKernel, // input square pixel array, preferrably having many exact zeros (they will be skipped)
									variableSigmas, // array of sigmas to be used for each pixel, matches pixels[]
									3.5, // drop calculation if farther then nSigma
									0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
									0, // int WOICenterY, // 
									inverseParameters.rSize, //int WOIWidth, reduce later
									inverseParameters.rSize,
									globalDebugLevel); //int WOIHeight)

						}
						//TODO: verify if it is OK that sum changed (was 10.5) after variableGaussBlurr(), for now - just re-calibrate
						normalizeKernel(rKernel); // in-place

						storeOneKernel( outPixels[chn], // float [] array of combined square kernels - will be filled
								rKernel, // square kernel to store
								tilesX, // number of kernels in a row
								tileX, // horizontal number of kernel to store
								tileY); // vertical number of kernel to store

					}
				}
			};
		}
		startAndJoin(threads);
		//	  System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
/* prepare result stack to return */
		final ImageStack outStack=new ImageStack(tilesX*inverseParameters.rSize,tilesY*inverseParameters.rSize);
		for (int chn=0;chn<nChn;chn++) {
			outStack.addSlice(PSFStack.getSliceLabel(chn+1), outPixels[chn]);
		}
		return outStack;
	}
  	
	/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
	private double [] maskReversePSFKernel( double []rpsf_pixels, // reversed psf, square array
			double [] ellipse_coeff, // ellipse coefficients from _direct_ kernel
			double ellipse_scale,
			double min_mask_threshold) // zero output element if elliptical Gauss mask is below this threshold
	{
		int rpsf_size=(int)Math.sqrt(rpsf_pixels.length);
		double [] masked_rpsf=new double[rpsf_size*rpsf_size];
		int ix,iy;
		double x,y,r2;
		int indx=0;
		double k2=1/ellipse_scale/ellipse_scale;
		double m;
		for (iy=0;iy<rpsf_size;iy++) {
			y=iy-rpsf_size/2+ellipse_coeff[1];  // move center opposite to that of direct kernel (psf)
			for (ix=0;ix<rpsf_size;ix++) {
				x=ix -rpsf_size/2 +ellipse_coeff[0]; //  move center opposite to that of direct kernel (psf)
				r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
				m=Math.exp(-k2*r2);
				masked_rpsf[indx]=(m>=min_mask_threshold)?(rpsf_pixels[indx]*Math.exp(-k2*r2)):0.0;
				indx++;
			}
		}
		return masked_rpsf;
	}

	/* ======================================================================== */
	
	
	private  double [] createSigmasFromCenter(
			int               size, // side of square
			double sigma_to_radius, // variable blurring - sigma will be proportional distance from the center
			double    center_sigma, //blurring in the center sigma(r)=sqrt((sigma_to_radius*r)^2)+center_sigma^2)
			double         centerX, // coordinates of the center (0:0 - size/2: size/2)
			double         centerY) {
		double [] sigmas = new double [size*size];
		int i,j;
		double x,y;
		double center_sigma2=center_sigma*center_sigma;
		double sigma_to_radius2=sigma_to_radius*sigma_to_radius;
		for (i=0;i<size;i++) for (j=0;j<size;j++) {
			y=i-size/2-centerY;
			x=j-size/2-centerX;
			sigmas[i*size+j]=Math.sqrt((x*x+y*y)*sigma_to_radius2+ center_sigma2);
		}
		return sigmas;
	}

	
	
	/* ======================================================================== */
	public double [] cleanupAndReversePSF (double []   psf_pixels,  // input pixels
			EyesisAberrations.InverseParameters inverseParameters, // size (side of square) of direct PSF kernel
			DoubleFHT fht_instance,  // provide DoubleFHT instance to save on initializations (or null)
			String           title,   // just for the plot names
			int debugLevel
	) {
		int size=(int) Math.sqrt(psf_pixels.length);
		double[][][] fft_complex;
		int i,j,ix,iy;
		double a,k,r,r2,k2;

		double [] cpixels=psf_pixels.clone();
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
/* Swapping quadrants, so the center will be 0,0 */
		fht_instance.swapQuadrants(cpixels);
/* get to frequency domain */
		fht_instance.transform(cpixels);
/* Convert from FHT to complex FFT - avoid that in the future, process FHT directly*/
		fft_complex= FHT2FFTHalf (cpixels,size);
		double [][]fft_energy=new double[(size/2)+1][size];
		for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
			fft_energy[i][j]=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
		}
		int  [][] clusterPS = findClusterOnPS(fft_energy, inverseParameters.otfCutoffEnergy,title,debugLevel);
		double [] ellipse_coeff = findEllipseOnPS(fft_energy, clusterPS, title,debugLevel);
/* create ellipse window using Hamming */
/* TODO: scale radius */
		double [][] ellipseMask=new double [size/2+1][size];
		k2=1/inverseParameters.otfEllipseScale/inverseParameters.otfEllipseScale;
		for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
			iy=(i==size/2)?-i:i;
			ix=(j>=(size/2))?(j-size):j;
			if (iy<0) ix=-ix;
			r2=ellipse_coeff[0]*ix*ix+ellipse_coeff[1]*iy*iy+ellipse_coeff[2]*ix*iy;
			if (inverseParameters.otfEllipseGauss){
				ellipseMask[i][j]=Math.exp(-k2*r2);
			} else {
				r=Math.sqrt(r2)/inverseParameters.otfEllipseScale;
				ellipseMask[i][j]=(r>1.0)?0.0:(0.54+0.46*Math.cos(r*Math.PI));
			}
		}
/* optionally display selection */
		if (debugLevel>2) {
			ImageProcessor ip_ellipse = new FloatProcessor(size,size);
			float [] ellipsePixels = new float [size*size];
			for (i=0;i<ellipsePixels.length;i++) {
				iy=i/size-size/2;
				ix=i%size-size/2;
				if (iy<0) {
					ix=-ix;
					iy=-iy;
				}
				ix= (ix+size) % size;
				ellipsePixels[i]= (float) ellipseMask[iy][ix];
			}
			ip_ellipse.setPixels(ellipsePixels);
			ip_ellipse.resetMinAndMax();
			ImagePlus imp_ellipse= new ImagePlus(title+"_EL-MASK_"+ inverseParameters.otfCutoffEnergy+"-"+inverseParameters.otfEllipseScale, ip_ellipse);
			imp_ellipse.show();
		}

/* inverse fft_complex */
		if (inverseParameters.invertRange>0.0) {
			/// Invert Z for large values, but make them Z - for small ones. So it will be a mixture of correlation and deconvolution
			//here the targets are round, but what will th\be the colrrect way fo assymmetrical ones?
			/// First - find maximal value
			double fft_max=0;
			for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
				r2=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
				if (r2>fft_max) fft_max=r2;
			}
			k=Math.sqrt(fft_max)*inverseParameters.invertRange;
			k2=k*k;
			for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
				r=Math.sqrt(fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1]);
				a=-Math.atan2(fft_complex[i][j][1],fft_complex[i][j][0]); /// was zero for circular targets)
				r=r/(r*r+k2);
				fft_complex[i][j][0]=r*Math.cos(a);
				fft_complex[i][j][1]=r*Math.sin(a);
			}
/* multiply by ellipse window */
			for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) {
				fft_complex[i][j][0]*=ellipseMask[i][j];
				fft_complex[i][j][1]*=ellipseMask[i][j];
			}
		} else { // Do just the division (low power frequencies will be masked out by ellipse window)
			for (i=0;i<fft_complex.length; i++) for (j=0;j<fft_complex[0].length;j++) if (ellipseMask[i][j]>=0.0){
				r2=fft_complex[i][j][0]*fft_complex[i][j][0]+fft_complex[i][j][1]*fft_complex[i][j][1];
				fft_complex[i][j][0]*= ellipseMask[i][j]/r2;
				fft_complex[i][j][1]*=-ellipseMask[i][j]/r2;
			} else {
				fft_complex[i][j][0]=0.0;
				fft_complex[i][j][1]=0.0;
			}
		}

		double [] pixels=null;
/* convert back original dimension array if there was no decimation or debug is set (in that case both sizes arrays will be converted) */
/* Convert fft array back to fht array and 
    set fht pixels with new values */
	    pixels=FFTHalf2FHT (fft_complex,size);
/* optionally show the result FHT*/
/* transform to space */
		fht_instance.inverseTransform(pixels);
		fht_instance.swapQuadrants(pixels);
/*   return inverted psf pixels */
		return pixels;
	}
	
	/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

	/* finds cluster (with the center at DC)  by flooding from DC, so total energy is cutoff_energy fraction
	returns integer array (same dimensions as input) with 1 - selected, 0 - not selected */
		private int [][] findClusterOnPS(
				double [][]       ps, // half power spectrum, starting from 0.0 (DC)
				double cutoff_energy, // fraction of energy in the pixels to be used
				String         title,
				int       debugLevel) {
			int i,j;
			List <Integer> pixelList=new ArrayList<Integer>(100);
			Integer Index;
			int size=ps[0].length;
			int [][]clusterMap=new int[size/2+1][size];
			double full_energy=0.0;
			int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
			for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
				full_energy+=((i%(size/2))==0)?ps[i][j]:(2*ps[i][j]); /* first and last line are counted once, others - twice */
				clusterMap[i][j]=0;
			}
			double threshold=full_energy*cutoff_energy;
			double cluster_energy=0.0;
			double maxValue;
			int ix,iy,ix1,iy1,maxX, maxY;
			int clusterSize=0;
			ix=0;
			iy=0;
			maxX=0;
			maxY=0;
			int listIndex;
			Index=iy*size + ix;
			pixelList.clear();
			pixelList.add (Index);
			clusterSize++;
			clusterMap[iy][ix]=1;
			cluster_energy+=ps[iy][ix];
			boolean noNew=true;
			while ((pixelList.size()>0) &&  (cluster_energy<threshold)) {
	/* Find maximal new neighbor */
				maxValue=0.0;
				listIndex=0;
				while (listIndex<pixelList.size()) {
					Index=pixelList.get(listIndex);
					iy=Index/size;
					ix=Index%size;
					noNew=true;
					for (j=0;j<8;j++) if (((iy > 0 ) || (dirs[j][1]>=0)) && ((iy < (size/2) ) || (dirs[j][1]<=0))){
						ix1=(ix+dirs[j][0]+size) % size;
						iy1= iy+dirs[j][1];
						if (clusterMap[iy1][ix1]==0) {
							noNew=false;
							if (ps[iy1][ix1]>maxValue) {
								maxValue= ps[iy1][ix1];
								maxX=ix1;
								maxY=iy1;
							}
						}
					}
					if (noNew) pixelList.remove(listIndex);  //  remove current list element
					else       listIndex++;     // increase list index
				}
				if (maxValue==0.0) { // Should
					System.out.println("findClusterOnPS: - should not get here - no points around >0, and threshold is not reached yet.");
					break;
				}
	/* Add this new point to the list */
				Index=maxY*size + maxX;
				pixelList.add (Index);
				clusterSize++;
				clusterMap[maxY][maxX]=1;
				cluster_energy+=((maxY%(size/2))==0)?ps[maxY][maxX]:(2*ps[maxY][maxX]);
			} // end of while ((pixelList.size()>0) &&  (cluster_energy<threshold))
			if (debugLevel>3)   System.out.println("findClusterOnPS: cluster size is "+clusterSize);
			if (debugLevel>6) {
				ImageProcessor ip2 = new FloatProcessor(size,size/2+1);
				float [] floatPixels = new float [size*(size/2+1)];
				for (i=0;i<floatPixels.length;i++) {
					floatPixels[i]=(float) ps[i/size][i%size];
				}
				ip2.setPixels(floatPixels);
				ip2.resetMinAndMax();
				ImagePlus imp2= new ImagePlus(title+"_PS1_"+cutoff_energy, ip2);
				imp2.show();
			}
			if (debugLevel>6) {
				ImageProcessor ip1 = new FloatProcessor(size,size);
				float [] floatPixels = new float [size*size];
				for (i=0;i<floatPixels.length;i++) {
					iy=i/size-size/2;
					ix=i%size-size/2;
					if (iy<0) {
						ix=-ix;
						iy=-iy;
					}
					ix= (ix+size) % size;
					floatPixels[i]=(float) ps[iy][ix];
				}
				ip1.setPixels(floatPixels);
				ip1.resetMinAndMax();
				ImagePlus imp1= new ImagePlus(title+"_PS_"+cutoff_energy, ip1);
				imp1.show();
			}

			if (debugLevel>5) {
				ImageProcessor ip = new FloatProcessor(size,size);
				float [] floatPixels = new float [size*size];
				for (i=0;i<floatPixels.length;i++) {
					iy=i/size-size/2;
					ix=i%size-size/2;
					if (iy<0) {
						ix=-ix;
						iy=-iy;
					}
					ix= (ix+size) % size;
					floatPixels[i]=(float) clusterMap[iy][ix];
				}
				ip.setPixels(floatPixels);
				ip.resetMinAndMax();
				ImagePlus imp= new ImagePlus(title+"_SEL_"+cutoff_energy, ip);
				imp.show();
			}
			return clusterMap;
		}

	/* calculates ellipse (with the center at DC) that interpolates area of the points defined by flooding from DC, so total energy is cutoff_energy fraction
	returns {a,b,c} , where a*x^2+b*y^2 + c*x*y=r^2 , so r^2 can be used for a window that removes high frequancy components that are too low to be useful*/

		private double [] findEllipseOnPS(
				double [][]        ps,   // half power spectrum, starting from 0.0 (DC)
				int    [][] selection, // 0/1 - selected/not selected
				String          title,
				int debugLevel) {
			int i,j;
			double x,y;
			int size=ps[0].length;
			double SX2=0.0;
			double SY2=0.0;
			double SXY=0.0;
			double S0=0.0;
			double k=2.0;
			double d;
			double area=0; // selection area
			for (i=0;i<(size/2+1);i++) {
				k=((i%(size/2))==0)?1.0:2.0;
				y=i;
				for (j=0;j<size;j++) if (selection[i][j]>0){
					x=(j>(size/2))?(j-size):j;
					d=k*ps[i][j];
					S0+=d;
					SX2+=x*x*d;
					SY2+=y*y*d;
					SXY+=x*y*d;
					area+=1.0;
				}
			}
			if (debugLevel>5) {
				System.out.println("findEllipseOnPS: title="+title+" area="+area+" S0="+S0+" SX2="+SX2+" SY2="+SY2+" SXY="+SXY);
			}
			//k=Math.PI*Math.PI/(2.0*S0*S0*area*area);
			//double [] result = {k*SY2,k*SX2,2*k*SXY};
			k=Math.PI*Math.PI/(2.0*S0*area*area);
			double [] result = {k*SY2,k*SX2,-2*k*SXY};
			if (debugLevel>3) {
				System.out.println("findEllipseOnPS: title="+title+" a="+result[0]+" b="+result[1]+" c="+result[2]);
			}
			return result;
		}

	
	
	
	/* ======================================================================== */
	/* TODO: REPLACE doubleFHT  */
	/* converts FHT results (frequency space) to complex numbers of [fftsize/2+1][fftsize] */

		private double[][][] FHT2FFTHalf (double [] fht_pixels, int fftsize) {
			double[][][] fftHalf=new double[(fftsize>>1)+1][fftsize][2];
			int row1,row2,col1,col2;

			for (row1=0;row1<=(fftsize>>1);row1++) {
				row2=(fftsize-row1) %fftsize;
				for (col1=0;col1<fftsize;col1++) {
					col2=(fftsize-col1) %fftsize;
					fftHalf[row1][col1][0]=   0.5*(fht_pixels[row1*fftsize+col1] + fht_pixels[row2*fftsize+col2]);
					fftHalf[row1][col1][1]=   0.5*(fht_pixels[row2*fftsize+col2] - fht_pixels[row1*fftsize+col1]);
				}
			}
			return fftHalf;
		}


		private double[] FFTHalf2FHT (double [][][] fft, int fftsize) {
			double[] fht_pixels=new double [fftsize*fftsize];
			int row1,row2,col1,col2;
			for (row1=0;row1<=(fftsize>>1);row1++) {
				row2=(fftsize-row1) %fftsize;
				for (col1=0;col1 < fftsize;col1++) {
					col2=(fftsize-col1) %fftsize;
					fht_pixels[row1*fftsize+col1]=(double) (fft[row1][col1][0]-fft[row1][col1][1]);
					fht_pixels[row2*fftsize+col2]=(double) (fft[row1][col1][0]+fft[row1][col1][1]);
				}
			}
			return fht_pixels;
		}


	
   	
   	/* interpolate kernels minimizing memory image - use directly the image stack (32-bit, float) with kernels.
   	  Add kernels around by either replication or extrapolation to compensate for "margins" in the original; kernels */
 //TODO: FIXME: Does not work if overwrite is disabled
   	public boolean interpolateKernels( 
		    AtomicInteger stopRequested, // 1 - stop now, 2 - when convenient
			EyesisAberrations.InterpolateParameters  interpolateParameters, // INTERPOLATE
			EyesisAberrations.MultiFilePSF           multiFilePSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
   			boolean                saveResult,
   			boolean                showResult,
   			boolean                updateStatus,          // UPDATE_STATUS
   			int                    globalDebugLevel
   	){
   		if ((this.aberrationParameters.psfKernelDirectory==null) || (this.aberrationParameters.psfKernelDirectory.length()==0)){
   			if (aberrationParameters.selectPSFKernelDirectory(true, this.aberrationParameters.partialKernelDirectory, false)==null) {
   				String msg = "Nothing selected";
   				System.out.println("Warning"+msg);
   				IJ.showMessage("Warning",msg);
   				return false;
   			}
   		}
   		int numChannels=distortions.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
    	boolean [] selectedChannels=this.aberrationParameters.getChannelSelection(distortions);
		String [] srcPaths=    new String[numChannels];
		String [] resultPaths= new String[numChannels];
		int numToProcess=0;
		for (int nChn=0;nChn<selectedChannels.length;nChn++){
			if (!selectedChannels[nChn]){
				srcPaths[nChn]=null;
				resultPaths[nChn]=null;
			} else {
				srcPaths[nChn]=this.aberrationParameters.psfKernelDirectory+Prefs.getFileSeparator()+
				this.aberrationParameters.psfPrefix+String.format("%02d", nChn)+
				this.aberrationParameters.psfSuffix;
				resultPaths[nChn]=this.aberrationParameters.psfKernelDirectory+Prefs.getFileSeparator()+
				this.aberrationParameters.interpolatedPSFPrefix+String.format("%02d", nChn)+
				this.aberrationParameters.interpolatedPSFSuffix;
				if (!this.aberrationParameters.overwriteResultFiles && (new File(resultPaths[nChn])).exists()) {
					srcPaths[nChn]=null;
					if (globalDebugLevel>0) System.out.println("File "+resultPaths[nChn]+" already exists and overwrite is disabled in configuration, channel "+nChn+" will be skipped");
					continue;
				}
				numToProcess++;
			}
		}
    	if (numToProcess==0){
				String msg = "No kernels to process";
   				System.out.println("Warning"+msg);
   				IJ.showMessage("Warning",msg);
   				return false;
    	}
   		int numProcessed=0;
   		Opener opener=new Opener();;
   		ImagePlus impSpsf;
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
   		for (int nChn=0;nChn<numChannels;nChn++) if (srcPaths[nChn]!=null){
   			if (!(new File(srcPaths[nChn])).exists()) {
   				String msg = "Combined PSF kernel stack for channel #"+nChn+": "+srcPaths[nChn]+" does not exist";
   				System.out.println("Warning"+msg);
   				continue;
   			}
   			impSpsf=opener.openImage("", srcPaths[nChn]);
   			if (impSpsf==null) {
   				System.out.println("Failed to open raw PSF kernel stack "+srcPaths[nChn]);
   				continue;
   			}
   			if (impSpsf.getStackSize()<3) {
   				System.out.println("Need a 3-layer stack with raw PSF kernels");
   				continue;
   			}
   			ImageStack stack= interpolateKernelStack(
   					impSpsf.getStack(), // Image stack, each slice consists of square kernels of one channel
   					interpolateParameters,
   					updateStatus,
   					globalDebugLevel); // update status info

   			ImagePlus impInterpolatedPSF = new ImagePlus("interpolated kernel stack", stack);
   			JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
   			jp4_instance.decodeProperiesFromInfo(impSpsf);
   			// copy properties from the source image
   			jp4_instance.copyProperties (impSpsf,impInterpolatedPSF);
   			interpolateParameters.setProperties("INTERPOLATE.", impInterpolatedPSF);
   			jp4_instance.encodeProperiesToInfo(impInterpolatedPSF);
   			if (showResult) {
   				impInterpolatedPSF.getProcessor().resetMinAndMax(); // imp_psf will be reused
   				impInterpolatedPSF.show();
   			} 
   			if (saveResult){
   				if (globalDebugLevel>0) System.out.println((numProcessed+1)+" of "+numToProcess+": saving interpolation result (of the file"+srcPaths[nChn]+") to "+
   						resultPaths[nChn]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
   				FileSaver fs=new FileSaver(impInterpolatedPSF);
   				fs.saveAsTiffStack(resultPaths[nChn]);
   			}
   			numProcessed++;
    		if 	(stopRequested.get()>0) {				
				if (globalDebugLevel>0) System.out.println("User requested stop");
				break;
    		}
   		}
   		if (numProcessed>0){
			if (globalDebugLevel>0) {
				System.out.println("Interpolated "+numProcessed+" kernel stacks at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
			}
   		} else {
				String msg = "No kernel files to interpolate";
   				System.out.println("Warning"+msg);
   				return false;
   		}
   		return true;
   	}




   		public ImageStack interpolateKernelStack(
   				ImageStack kernelStack, // Image stack, each slice consists of square kernels of one channel
   				EyesisAberrations.InterpolateParameters interpolateParameters,
   				boolean   updateStatus,
   				int globalDebugLevel) // update status info
   		{
   			DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
   			if (kernelStack==null) return null;
   			int inTilesX=kernelStack.getWidth()/interpolateParameters.size;
   			int inTilesY=kernelStack.getHeight()/interpolateParameters.size;
   			int outTilesX= (inTilesX-1)*interpolateParameters.step +1 + interpolateParameters.add_left + interpolateParameters.add_right;
   			int outTilesY= (inTilesY-1)*interpolateParameters.step +1 + interpolateParameters.add_top + interpolateParameters.add_bottom;
   			int nChn=kernelStack.getSize();
   			float [][] outPixels=new float[nChn][outTilesX*interpolateParameters.size*outTilesY*interpolateParameters.size];
   			float [] pixels;
   			int i,j,chn;
   			int xTile0=(interpolateParameters.add_left>0)?-1:0;
   			int xTile1=inTilesX+((interpolateParameters.add_right>0)?0:-1);
   			int yTile0=(interpolateParameters.add_top>0)?-1:0;
   			int yTile1=inTilesY+((interpolateParameters.add_bottom>0)?0:-1);
   			int tileY,tileX; //,subTileY,subTileX;

   			int tileWidth, tileHeight; // for inner cells (interpolateParameters.step+1)*(interpolateParameters.step+1), for outer includes exte row/column fro extrapolation
   			//  int maxTileWidth= Math.max(interpolateParameters.step,1+Math.max(interpolateParameters.add_right,interpolateParameters.add_left));
   			//  int maxTileHeight=Math.max(interpolateParameters.step,1+Math.max(interpolateParameters.add_bottom,interpolateParameters.add_top));
   			boolean lastColumn=false;  //last column - inverse convert and copy the last column of rectangleFHT to the result array
   			boolean lastRow=false;     //last row - interpolate, inverse convert and copy the last row of rectangleFHT to the result array

   			double [] pointsVert;
   			double [] pointsHor;
   			double [][] fhtLine;
   			double extraScale=interpolateParameters.extrapolate/interpolateParameters.step;
   			int [] outTopLeft=new int [2]; // top left kernel in the output array
   			int [] inTopLeft=new int [2]; // top left kernel in the input array
   			double [][] firstFHTColumn=null;
   			double [][] secondFHTColumn=null;
   			double [][][] cornerFHT=new double[2][2][interpolateParameters.size*interpolateParameters.size]; //[y][x][pixel] 
   			double [] swapArray=null;

   			for (chn=0;chn<nChn;chn++) {
   				pixels=(float[]) kernelStack.getPixels(chn+1);
   				for (tileY=yTile0;tileY<yTile1;tileY++) {
   					if (updateStatus) IJ.showStatus("Interpolating kernels, channel "+kernelStack.getSliceLabel(chn+1)+", row "+(tileY-yTile0+1)+" of "+(yTile1-yTile0));
   					lastRow=(tileY==(yTile1-1));
   					if (tileY<0) {
   						inTopLeft[1]=0;
   						tileHeight=interpolateParameters.add_top;
   						outTopLeft[1]=0;
   						pointsVert=new double[tileHeight];
   						for (i=0;i<tileHeight;i++)  pointsVert[i]=(i-tileHeight)*extraScale; // negative values
   					} else if (tileY>=(inTilesY-1)){
   						inTopLeft[1]=tileY-1;
   						tileHeight=interpolateParameters.add_bottom+1; // always last row, if got here at all (interpolateParameters.add_bottom>0)
   						outTopLeft[1]=interpolateParameters.add_top+interpolateParameters.step*tileY;
   						pointsVert=new double[tileHeight];
   						for (i=0;i<tileHeight;i++)  pointsVert[i]=1.0+i*extraScale;
   					} else {
   						inTopLeft[1]=tileY;
   						tileHeight=interpolateParameters.step+ (lastRow?1:0); // last tile row includes bottom outpout kernel row
   						outTopLeft[1]=interpolateParameters.add_top+interpolateParameters.step*tileY;
   						pointsVert=new double[tileHeight];
   						for (i=0;i<tileHeight;i++) pointsVert[i]=(1.0*i)/tileHeight;
   					}
   					firstFHTColumn=null; // invalidate
   					secondFHTColumn=null; // invalidate
   					for (tileX=xTile0;tileX<xTile1;tileX++) {
   						if (globalDebugLevel>2)  System.out.println(" interpolateKernelStack(): chn="+chn+" tileY="+tileY+" tileX="+tileX);

   						lastColumn=(tileX==(xTile1-1));
   						if (tileX<0) {
   							inTopLeft[0]=0;
   							tileWidth=interpolateParameters.add_left;
   							outTopLeft[0]=0;
   							pointsHor=new double[tileWidth];
   							for (i=0;i<tileWidth;i++)  pointsHor[i]=(i-tileWidth)*extraScale; // negative values
   						} else if (tileX>=(inTilesX-1)){
   							inTopLeft[0]=tileX-1;
   							tileWidth=interpolateParameters.add_right+1; // always last columnw, if got here at all (interpolateParameters.add_right>0)
   							outTopLeft[0]=interpolateParameters.add_left+interpolateParameters.step*tileX;
   							pointsHor=new double[tileWidth];
   							for (i=0;i<tileWidth;i++)  pointsHor[i]=1.0+ i*extraScale;
   							// else keep both firstFHTColumn and secondFHTColumn
   							if (globalDebugLevel>2)  System.out.println("last column: tileX="+tileX);
   						} else {
   							inTopLeft[0]=tileX;
   							tileWidth=interpolateParameters.step+ (lastColumn?1:0); // last tile column includes rightmost outpout kernel column
   							outTopLeft[0]=interpolateParameters.add_left+interpolateParameters.step*tileX;
   							pointsHor=new double[tileWidth];
   							for (i=1;i<tileWidth;i++)  pointsHor[i]=(1.0*i)/tileWidth;
   							//  if (DEBUG_LEVEL>2)  System.out.println("else: tileX="+tileX);
   							if (tileX!=0) {
   								firstFHTColumn=secondFHTColumn;
   								secondFHTColumn=null; // invalidate
   								//  if (DEBUG_LEVEL>2)  System.out.println(" secondFHTColumn==null");
   	/* swap columns, the new second one will be just reused */
   								swapArray=cornerFHT[0][0];
   								cornerFHT[0][0]=cornerFHT[0][1];
   								cornerFHT[0][1]=swapArray;
   								swapArray=cornerFHT[1][0];
   								cornerFHT[1][0]=cornerFHT[1][1];
   								cornerFHT[1][1]=swapArray;

   							} // else keep both firstFHTColumn and secondFHTColumn
   						}
   						if (globalDebugLevel>2)  System.out.println(" interpolateKernelStack(): tileHeight="+tileHeight+" tileWidth="+tileWidth+" inTopLeft[0]="+inTopLeft[0]+" inTopLeft[1]="+inTopLeft[1]+
   								" outTopLeft[0]="+outTopLeft[0]+" outTopLeft[1]="+outTopLeft[1]);

   						if (firstFHTColumn==null) { /* First colum needs to be input and calculated*/
   							extractOneKernel(          pixels, //  array of combined square kernels, each 
   									cornerFHT[0][0], // will be filled, should have correct size before call
   									inTilesX, // number of kernels in a row
   									inTopLeft[0], // horizontal number of kernel to extract
   									inTopLeft[1]); // vertical number of kernel to extract
   							extractOneKernel(          pixels, //  array of combined square kernels, each 
   									cornerFHT[1][0], // will be filled, should have correct size before call
   									inTilesX, // number of kernels in a row
   									inTopLeft[0], // horizontal number of kernel to extract
   									inTopLeft[1]+1); // vertical number of kernel to extract
   	/* convert to frequency domain */
   							fht_instance.swapQuadrants(cornerFHT[0][0]);
   							fht_instance.transform(    cornerFHT[0][0]);
   							fht_instance.swapQuadrants(cornerFHT[1][0]);
   							fht_instance.transform(    cornerFHT[1][0]);
   	/* inter/extrapolate the column */
   							firstFHTColumn=fht_instance.interpolateFHT (cornerFHT[0][0],    // first FHT array
   									cornerFHT[1][0],    // second FHT array
   									pointsVert,    // array of interpolation points - 0.0 - fht0, 1.0 - fht1
   									false);   // OK not to clone, so corners will be referenced?
   							if (globalDebugLevel>2)  System.out.println(" firstFHTColumn.length="+firstFHTColumn.length);
   						}
   						if (secondFHTColumn==null) { /* Last colum needs to be input and calculated*/
   							extractOneKernel(          pixels, //  array of combined square kernels, each 
   									cornerFHT[0][1], // will be filled, should have correct size before call
   									inTilesX, // number of kernels in a row
   									inTopLeft[0]+1, // horizontal number of kernel to extract
   									inTopLeft[1]); // vertical number of kernel to extract
   							extractOneKernel(          pixels, //  array of combined square kernels, each 
   									cornerFHT[1][1], // will be filled, should have correct size before call
   									inTilesX, // number of kernels in a row
   									inTopLeft[0]+1, // horizontal number of kernel to extract
   									inTopLeft[1]+1); // vertical number of kernel to extract
   	/* convert to frequency domain */
   							fht_instance.swapQuadrants(cornerFHT[0][1]);
   							fht_instance.transform(    cornerFHT[0][1]);
   							fht_instance.swapQuadrants(cornerFHT[1][1]);
   							fht_instance.transform(    cornerFHT[1][1]);
   	/* inter/extrapolate the column */
   							secondFHTColumn=fht_instance.interpolateFHT (cornerFHT[0][1],    // first FHT array
   									cornerFHT[1][1],    // second FHT array
   									pointsVert,    // array of interpolation points - 0.0 - fht0, 1.0 - fht1
   									false);   // OK not to clone, so corners will be referenced?

   							if (globalDebugLevel>2)  {
   								System.out.println(" secondFHTColumn.length="+secondFHTColumn.length);
   								for (i=0;i<pointsVert.length;i++) System.out.println(""+pointsVert[i]);
   								System.out.println("");
   							}
   						}
   	/* interpolate horizontally */
   	/* TODO: calculate top-left corner in output array */
   						/*
   	   if ((DEBUG_LEVEL>1) &&(tileY==0)) {
   	      SDFA_instance.showArrays(firstFHTColumn,size,size, "firstFHTColumn");
   	      SDFA_instance.showArrays(secondFHTColumn,size,size, "secondFHTColumn");
   	      DEBUG_LEVEL=4;
   	      return null;
   	   }
   						 */
   						for (i=0;i<tileHeight;i++) {
   							if (globalDebugLevel>2)  System.out.print("i="+i);

   							fhtLine=fht_instance.interpolateFHT ( firstFHTColumn[i],    // first FHT array
   									secondFHTColumn[i],    // second FHT array
   									pointsHor,    // array of interpolation points - 0.0 - fht0, 1.0 - fht1
   									true); //clone ends
   							if (globalDebugLevel>2)  System.out.print(": ");
   							for (j=0;j<tileWidth;j++) {
   								if (globalDebugLevel>2)  System.out.print(j);
   								fht_instance.inverseTransform(fhtLine[j]);
   								fht_instance.swapQuadrants   (fhtLine[j]);
   								storeOneKernel(           outPixels[chn], // float [] array of combined square kernels - will be filled
   										fhtLine[j], // square kernel to store
   										outTilesX, // number of kernels in a row
   										outTopLeft[0]+j, // horizontal number of kernel to store
   										outTopLeft[1]+i); // vertical number of kernel to store
   							}
   							if (globalDebugLevel>2)  System.out.println("");

   						}
   					}
   				}
   			}    
   	/* prepare result stack to return */
   			ImageStack outStack=new ImageStack(outTilesX*interpolateParameters.size,outTilesY*interpolateParameters.size);
   			for (chn=0;chn<nChn;chn++) {
   				outStack.addSlice(kernelStack.getSliceLabel(chn+1), outPixels[chn]);
   			}
   			return outStack;
   		}
   	/* ======================================================================== */
   	/* Used in interpolateKernelStack() */  
   		private void storeOneKernel(
   				float []  pixels, // float [] array of combined square kernels - will be filled
   				double [] kernel, // square kernel to store
   				int       numHor, // number of kernels in a row
   				int        xTile, // horizontal number of kernel to store
   				int        yTile) { // vertical number of kernel to store
   			int length=kernel.length;
   			int size=(int) Math.sqrt(length);
   			int i,j;
   			int pixelsWidth=numHor*size;
   			int base=(yTile*pixelsWidth+xTile)*size;
   			for (i=0;i<size;i++) for (j=0;j<size;j++) pixels[base+i*pixelsWidth+j]= (float) kernel[i*size+j];
   		}

   	/* ======================================================================== */
   	
   	
   	public String [][]  preparePartialKernelsFilesList(
   			int debugLevel){
   		DistortionCalibrationData distortionCalibrationData= distortions.fittingStrategy.distortionCalibrationData;
   		boolean [] selectedImages=distortions.fittingStrategy.selectedImagesNoBadKernels(this.aberrationParameters.allImages?-1:this.aberrationParameters.seriesNumber); // negative series number OK - will select all enabled
   		int num=0;
   		for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) num++;
   		if (debugLevel>1) {
   			System.out.println("Enabled "+num+" source files");
   		}
   		if (num==0){
				String msg="No enabled files selected. Command aborted";
   				System.out.println("Warning"+msg);
   				IJ.showMessage("Warning",msg);
   				return null;
   		}
   		String [] partialKernelPaths=new String [selectedImages.length];
   		for (int imgNum=0;imgNum<partialKernelPaths.length;imgNum++){
   			if (!selectedImages[imgNum]) {
   				partialKernelPaths[imgNum]=null;
   			} else {
   				partialKernelPaths[imgNum]=this.aberrationParameters.partialPrefix+IJ.d2s(distortionCalibrationData.gIP[imgNum].timestamp,6).replace('.','_')+
   				String.format("-%02d"+this.aberrationParameters.partialSuffix, distortionCalibrationData.gIP[imgNum].channel); // sensor number
   				//   			partialKernelPaths[imgNum]=this.aberrationParameters.sourceDirectory+Prefs.getFileSeparator()+filename;
   				
   		   		if (debugLevel>2) System.out.println("preparePartialKernelsFilesList() "+imgNum+": "+partialKernelPaths[imgNum]);

   			}
   		}
   		if (countExistentFiles(this.aberrationParameters.partialKernelDirectory,partialKernelPaths,false)==0){ // keep non-existent
   			if (aberrationParameters.selectPartialKernelDirectory(true, this.aberrationParameters.partialKernelDirectory, false)==null) {
   				String msg="No partial kernel directory selected. Command aborted";
   				System.out.println("Warning"+msg);
   				IJ.showMessage("Warning",msg);
   				return null;
   			}
   		}
   		if (countExistentFiles(this.aberrationParameters.partialKernelDirectory,partialKernelPaths,true)==0){ // will remove all non-existent files
   			String msg="No partial kernel files found. Command aborted";
   			System.out.println("Warning"+msg);
   			IJ.showMessage("Warning",msg);
   			return null;
   		}
   		
   		
   		int numChannels=distortions.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
   		String [][] fileList=new String[numChannels][];
   		for (int numChn=0;numChn<numChannels;numChn++){
   			int n=0;
   			for (int  imgNum=0;imgNum<partialKernelPaths.length;imgNum++)
   				if ((partialKernelPaths[imgNum]!=null) && (distortionCalibrationData.gIP[imgNum].channel==numChn))n++;
   			if (n==0) {
   				fileList[numChn]=null;
   			} else {
   				fileList[numChn]=new String[n];
   				n=0;
   				for (int  imgNum=0;imgNum<partialKernelPaths.length;imgNum++)
   					if ((partialKernelPaths[imgNum]!=null) && (distortionCalibrationData.gIP[imgNum].channel==numChn)){
   						fileList[numChn][n++]=this.aberrationParameters.partialKernelDirectory+Prefs.getFileSeparator()+partialKernelPaths[imgNum];
   					}
   			}
   		}
   		if (debugLevel>0) {
   			System.out.println("Partial kernel files available:");
   			for (int numChn=0;numChn<numChannels;numChn++){
   				if (fileList[numChn]!=null) {
   					System.out.println("   channel "+numChn+": "+fileList[numChn].length+" files");

   					if (debugLevel>1) {
   						for (int i=0;i<fileList[numChn].length;i++) {
   							System.out.println(numChn+":"+i+": "+fileList[numChn][i]);
   						}
   					}
   				}
   			}
   		}
   		return fileList;
   	}
/*
    	int numChannels=this.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
  
 */
	public boolean createPartialKernels(
		    AtomicInteger stopRequested, // 1 - stop now, 2 - when convenient
			int             mapFFTsize, // scanImageForPatterns:FFT size
			int            fft_overlap,
			int               fft_size,
			int           PSF_subpixel, 
			OTFFilterParameters otfFilterParameters,
			PSFParameters psfParameters,
			int          PSFKernelSize, // size of square used in the new map (should be multiple of map step)
			double       gaussWidth,  // ** NEW
			MultiFilePSF multiFilePSF,
			MatchSimulatedPattern.DistortionParameters distortionParameters, //
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			SimulationPattern.SimulParameters  simulParameters,
			ColorComponents colorComponents,
			boolean resetBadKernels, // ignore and reset noUsefulKernels mark for selected channel
			int threadsMax,
			boolean updateStatus,
			int loopDebugLevel, // debug level used inside loops
			int debugLevel
			){
    	DistortionCalibrationData distortionCalibrationData= distortions.fittingStrategy.distortionCalibrationData;
    	boolean partialToReprojected=this.aberrationParameters.partialToReprojected;
    	boolean applySensorCorrection=this.aberrationParameters.partialCorrectSensor;
    	// this.distortions is set to top level LENS_DISTORTIONS
		if (distortions==null){
    		String msg="Distortions instance does not exist, exiting";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);

		}
		if (distortions.fittingStrategy==null){
    		String msg="Fitting strategy does not exist, exiting";
    		IJ.showMessage("Error",msg);
    		throw new IllegalArgumentException (msg);
		}
		long startTime=System.nanoTime(); // restart timer after possible interactive dialogs
//		long tmpTime;
		//resetBadKernels
		int serNumber=this.aberrationParameters.allImages?-1:this.aberrationParameters.seriesNumber;
    	boolean [] selectedImages=resetBadKernels?
    			distortions.fittingStrategy.selectedImages(serNumber):
    			distortions.fittingStrategy.selectedImagesNoBadKernels(serNumber); // negative series number OK - will select all enabled
    	boolean [] selectedChannels=this.aberrationParameters.getChannelSelection(distortions);
    	int numSelected=0;
    	int numDeselected=0;
    	if (debugLevel>2){
    		for (int i=0;i<selectedChannels.length;i++){
    			System.out.println("Channel "+i+" is "+(selectedChannels[i]?"Enabled":"Disabled"));
    		}
    	}
    	for (int imgNum=0;imgNum<selectedImages.length;imgNum++) if (selectedImages[imgNum]) {
    		int numChannel=distortionCalibrationData.gIP[imgNum].channel;
        	if (debugLevel>2){
        		System.out.println("Image "+imgNum+" channel "+numChannel+" is "+(selectedChannels[numChannel]?"ENABLED":"DISABLED"));
        	}    		
    		if (!selectedChannels[numChannel]){
    			selectedImages[imgNum]=false;
    			numDeselected++;
    		}else{
    			distortions.fittingStrategy.setNoUsefulPSFKernels(imgNum,false); // reset noUsefulKernels mark (if it was not set - OK)
    			numSelected++;
    		}
    	} else if (debugLevel>2){
    		System.out.println("Skipping disabled image "+imgNum);
    	}
    	if (debugLevel>0)System.out.println("Enabled "+numSelected+" source files ("+numDeselected+") were removed by channel selection. partialToReprojected="+partialToReprojected);

    	String [] sourcePaths=new String [selectedImages.length];
    	// Set/verify source paths
    	int numFiles=0;
    	boolean skipMissing=false;
    	for (int imgNum=0;imgNum<sourcePaths.length;imgNum++){
    		if (!selectedImages[imgNum]) sourcePaths[imgNum]=null;
    		else {
    			String filename=this.aberrationParameters.sourcePrefix+IJ.d2s(distortionCalibrationData.gIP[imgNum].timestamp,6).replace('.','_')+
    			String.format("-%02d"+this.aberrationParameters.sourceSuffix, distortionCalibrationData.gIP[imgNum].channel); // sensor number
    			sourcePaths[imgNum]=this.aberrationParameters.sourceDirectory+Prefs.getFileSeparator()+filename;
    			File srcFile=new File(sourcePaths[imgNum]);	
    			if (!srcFile.exists()){
    				if (skipMissing) {
    					if (debugLevel>0) System.out.println("Skipping missing file: "+sourcePaths[imgNum]);
    	    	    	sourcePaths[imgNum]=null;
    					continue;
    				}
    				GenericDialog gd=new GenericDialog("Missing source file(s)");
    				gd.addMessage ("Source file "+sourcePaths[imgNum]+" does not exist");
    				gd.addMessage ("This may be a file from the different source directory (acquired at different station).");
    				gd.addMessage ("You may change source and destination directories and re-run this command for another station.");
    				gd.enableYesNoCancel("Find file", "Skip  this and other missing files");
    	    	    gd.showDialog();
    	    	    if (gd.wasCanceled()) return false;
    	    	    if (!gd.wasOKed()){
    	    	    	skipMissing=true;
    	    	    	sourcePaths[imgNum]=null;
    	    	    	continue;
    	    	    }
    				String [] extensions={filename}; // just this one file
    				CalibrationFileManagement.MultipleExtensionsFileFilter parFilter = new CalibrationFileManagement.MultipleExtensionsFileFilter("",extensions,filename);
    				String pathname=CalibrationFileManagement.selectFile(
    						false,
    						false,
    						"Find source file",
    						"Select",
    						parFilter,
    						this.aberrationParameters.sourceDirectory); //String defaultPath
    				if ((pathname==null) || (pathname=="")) return false;
    				sourcePaths[imgNum]=pathname;
    				this.aberrationParameters.sourceDirectory=pathname.substring(0, pathname.lastIndexOf(Prefs.getFileSeparator()));
    			}
				numFiles++;
    		}
    	}
    	if (numFiles==0 ){
    		String msg="createPartialKernels(): No files selected";
    		System.out.println("Warning: "+msg);
    		if (!aberrationParameters.noMessageBoxes)IJ.showMessage("Warning",msg);
    		return true;
    	}
    	if (debugLevel>0) {
    		System.out.println("Selected "+numFiles+" source files");
    	}

    	if (aberrationParameters.selectPartialKernelDirectory(true, this.aberrationParameters.partialKernelDirectory, true)==null){
    		String msg="createPartialKernels(): No partial kernel directory selected";
    		System.out.println("Warning: "+msg);
    		if (!aberrationParameters.noMessageBoxes)IJ.showMessage("Warning",msg);
    		return true;
    	}
    	String [] partialKernelsPaths=new String [selectedImages.length];
    	for (int imgNum=0;imgNum<sourcePaths.length;imgNum++){
    		if (sourcePaths[imgNum]==null){
    			partialKernelsPaths[imgNum]=null;
    		} else {
    			String filename=this.aberrationParameters.partialPrefix+IJ.d2s(distortionCalibrationData.gIP[imgNum].timestamp,6).replace('.','_')+
    			String.format("-%02d", distortionCalibrationData.gIP[imgNum].channel)+this.aberrationParameters.partialSuffix;
    			partialKernelsPaths[imgNum]=this.aberrationParameters.partialKernelDirectory+Prefs.getFileSeparator()+filename;
    		}
    	}
    	int numOld=0;
    	if (!this.aberrationParameters.overwriteResultFiles){
        	for (int imgNum=0;imgNum<sourcePaths.length;imgNum++) if (partialKernelsPaths[imgNum]!=null){
            	if (debugLevel>1){
            		System.out.println(imgNum+": "+partialKernelsPaths[imgNum]+((new File(partialKernelsPaths[imgNum]).exists())?" EXISTS":" DOES NOT EXIST"));
            	}
        		if (new File(partialKernelsPaths[imgNum]).exists()){
        			numOld++;
        			numFiles--;
        			partialKernelsPaths[imgNum]=null;
        			sourcePaths[imgNum]=null;
        		}
        	}
    	}
    	if (debugLevel>0){
    		System.out.println((numFiles+numOld)+" source files selected, "+((numOld>0)?(numOld+" existent files skipped, "):"")+numFiles+" to process");
    	}
    	if (numFiles<=0){
    		String msg="createPartialKernels(): No files to process";
    		System.out.println("Warning: "+msg);
    		if (!aberrationParameters.noMessageBoxes)IJ.showMessage("Warning",msg);
    		return true;
    	}
    	// reorder in the ascending channel number order
    	String [][] files=new String [numFiles][2]; // 0 - source, 1 - result
    	int [] fileIndices =new int [numFiles]; // needed to mark bad kernels (and also to reference grid parameters to replace extracted)
    	int numListedFiles=0;
    	int numChannel=0;
    	while (numListedFiles<numFiles) {
    		for (int imgNum=0;imgNum<sourcePaths.length;imgNum++) if ((sourcePaths[imgNum]!=null) && (distortionCalibrationData.gIP[imgNum].channel<=numChannel)){
    			if (debugLevel>1) System.out.println("numListedFiles="+numListedFiles+" numFiles"+numFiles+" sourcePaths["+imgNum+"]="+sourcePaths[imgNum]+" numChannel="+numChannel);
    			files[numListedFiles][0]=sourcePaths[imgNum];
    			files[numListedFiles][1]=partialKernelsPaths[imgNum];
    			fileIndices[numListedFiles++]=imgNum;
    			sourcePaths[imgNum]=null;
    		}
    		numChannel++;
    	}
    	startTime=System.nanoTime(); // restart timer after possible interactive dialogs
    	for (int imgNum=0;imgNum<files.length;imgNum++){ // add stopRequested
			if (debugLevel>0) System.out.println("Processing file #"+(imgNum+1)+ " ( of "+files.length+") :"+files[imgNum][0]);
        	ImagePlus imp=new ImagePlus(files[imgNum][0]); // read source file
        	JP4_INSTANCE.decodeProperiesFromInfo(imp);
// TODO: Add vignetting correction ?
        	MatchSimulatedPattern matchSimulatedPattern= new MatchSimulatedPattern(distortionParameters.FFTSize);
			boolean [] correlationSizesUsed=null;
			float [][] simArray=         	null;

        	int MaxRetries=4;
        	int iRetry=0;
        	for (iRetry=0;iRetry<MaxRetries;iRetry++){ // is this retry needed?
        		try {
        			
        			double [][][] projectedGrid=null;
        			double hintTolerance=0.0; 
        			if (partialToReprojected){ // replace px, py with projected values form the grid
        				// this.distortions is set to the global LENS_DISTORTIONS
        				int numGridImage=fileIndices[imgNum];
        				projectedGrid=distortions.estimateGridOnSensor( // return grid array [v][u][0- x,  1 - y, 2 - u, 3 - v] 
        						distortions.fittingStrategy.distortionCalibrationData.getImageStation(numGridImage), // station number,
        						distortions.fittingStrategy.distortionCalibrationData.gIP[numGridImage].channel, // subCamera,
        						Double.NaN, // goniometerHorizontal, - not used
        						Double.NaN, // goniometerAxial, - not used
        						Double.NaN, // inter-axis angle, - not used ?
        						distortions.fittingStrategy.distortionCalibrationData.gIP[numGridImage].getSetNumber(), //imageSet,
        						true); //filterBorder)
        				hintTolerance=5.0; // TODO:set from configurable parameter
        				if (applySensorCorrection){
        					boolean applied=distortions.correctGridOnSensor(
        							projectedGrid,
        							distortions.fittingStrategy.distortionCalibrationData.gIP[numGridImage].channel);
                			if (debugLevel>0) {
                				if (applied) System.out.println("Applied sensor correction to the projected grid");
                				else System.out.println("No sensor correction available to apply to the projected grid");
                			}
        				}
        			}
        			
        			int rslt=matchSimulatedPattern.calculateDistortions(
        					distortionParameters, //
        					patternDetectParameters,
        					simulParameters,
        					colorComponents.equalizeGreens,
        					imp,
        					null, // LaserPointer laserPointer, // LaserPointer object or null
        					true, // don't care -removeOutOfGridPointers
        					projectedGrid, // null, //   double [][][] hintGrid, // predicted grid array (or null)
        					hintTolerance, // 0,    //   double  hintGridTolerance, // allowed mismatch (fraction of period) or 0 - orientation only
        					threadsMax,
        					updateStatus,
        					debugLevel,
        					loopDebugLevel, // debug level
        					aberrationParameters.noMessageBoxes);
        			if (rslt<0){
            			if (debugLevel>0) System.out.println("calculateDistortions failed, returned error code "+rslt+" iRetry="+iRetry+" (of "+MaxRetries+")");
            			continue;
        			}
        			// now replace extracted grid X,Y with projected (need to add sensor correction)
        			if (projectedGrid!=null){
        				int numReplaced= matchSimulatedPattern.replaceGridXYWithProjected(projectedGrid,(debugLevel>1)?imp.getTitle():null);
            			if (debugLevel>0) System.out.println("Replaced extracted XY with projected ones for "+numReplaced+" nodes");
        			}
        			correlationSizesUsed=matchSimulatedPattern.getCorrelationSizesUsed();
        			simArray=	(new SimulationPattern(simulParameters)).simulateGridAll (
        					imp.getWidth(),
        					imp.getHeight(),
        					matchSimulatedPattern,
        					2, // gridFrac, // number of grid steps per pattern full period
        					simulParameters,
        					threadsMax,
        					updateStatus,
        					debugLevel,
        					loopDebugLevel); // debug level
        			
        			createPSFMap(
        					matchSimulatedPattern,
        					matchSimulatedPattern.applyFlatField (imp), // if grid is flat-field calibrated, apply it (may throw here)
        					null,     //  int [][][] sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  PSF_KERNEL_MAP structure)  
        					multiFilePSF.overexposedMaxFraction, //MULTIFILE_PSF.overexposedMaxFraction,
        					simulParameters, //SIMUL, //simulation parameters
        					mapFFTsize, // MAP_FFT_SIZE, // scanImageForPatterns:FFT size int             mapFFTsize, // scanImageForPatterns:FFT size
        					patternDetectParameters, //PATTERN_DETECT, //MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
        					fft_overlap, //FFT_OVERLAP, // int            fft_overlap,
        					fft_size, // FFT_SIZE, // int               fft_size,
        					colorComponents, //COMPONENTS,   // ColorComponents colorComponents,
        					PSF_subpixel, //PSF_SUBPIXEL, // int           PSF_subpixel, 
        					otfFilterParameters, // OTF_FILTER, // OTFFilterParameters otfFilterParameters,
        					psfParameters, //PSF_PARS, // final PSFParameters psfParameters
        					psfParameters.minDefinedArea , //PSF_PARS.minDefinedArea, // final double       minDefinedArea,
        					PSFKernelSize, //// int          PSFKernelSize, // size of square used in the new map (should be multiple of map step)
        					gaussWidth, //gaussWidth
        					simArray, // simArray
        					threadsMax,   // threadsMax,
        					updateStatus, // updateStatus,
        					debugLevel, //masterDebugLevel
        					debugLevel, //globalDebugLevel
        					loopDebugLevel);// debug level used inside loops
        			break; // success
        		} catch (Exception e) {
        			if (debugLevel>0) System.out.println("Attempt "+(iRetry+1)+" of "+MaxRetries+"Failed to find initial pattern in file #"+
        					(imgNum+1)+ " ( of "+files.length+") :"+files[imgNum][0]);
        			e.printStackTrace();
        			continue;
        		}
        	}
        	if (iRetry==MaxRetries) {
				System.out.println("File "+files[imgNum][1]+ " has problems - finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
	    		if 	(stopRequested.get()>0) {				
					if (debugLevel>0) System.out.println("User requested stop");
					break;
	    		}
        		continue;
        	}
        	
			
			ImageStack stack=mergeKernelsToStack(this.pdfKernelMap);
			
			// TODO: Add properties,
			// Save configuration (filename with timestamp?) before files from the top class, test directory is writable
			
			
			if (stack!=null) {
				if (debugLevel>0) System.out.println("Saving result to"+files[imgNum][1]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
				 savePartialKernelStack(
						 files[imgNum][1],
							stack,
							imp,
							psfParameters,
							correlationSizesUsed);
			} else {
				System.out.println("File "+files[imgNum][1]+ " has no useful PSF kernels - at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3)); 
				distortions.fittingStrategy.setNoUsefulPSFKernels(fileIndices[imgNum], true); // mark (need to save configuration) not to try them next time
// todo - write a placeholder file (different suffix/prefix) instead of using		setNoUsefulPSFKernels()?		
			}
    		if 	(stopRequested.get()>0) {				
				if (debugLevel>0) System.out.println("User requested stop");
				break;
    		}
    	}
		return true;
	}

	public boolean combinePSFKernels ( // save configuration to combined kernels directory before calling this method
		    AtomicInteger stopRequested, // 1 - stop now, 2 - when convenient
			EyesisAberrations.InterpolateParameters  interpolateParameters, // INTERPOLATE
			EyesisAberrations.MultiFilePSF           multiFilePSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
//			showDoubleFloatArrays  sdfa_instance,        // SDFA_INSTANCE
			boolean                saveResult,
			boolean                showResult,
			boolean                updateStatus,          // UPDATE_STATUS
			int                    thisDebugLevel,
			int                    globalDebugLevel
	){
		String [][] fileList=  preparePartialKernelsFilesList(
				globalDebugLevel);
    	boolean [] selectedChannels=this.aberrationParameters.getChannelSelection(distortions);
		String [] resultPaths= new String[fileList.length];
    	for (int nChn=0;nChn<selectedChannels.length;nChn++){
    		if (!selectedChannels[nChn]){
    			fileList[nChn]=null;
    		} else {
    			resultPaths[nChn]=this.aberrationParameters.psfKernelDirectory+Prefs.getFileSeparator()+
    			this.aberrationParameters.psfPrefix+String.format("%02d", nChn)+
    			this.aberrationParameters.psfSuffix;
    	    	if (!this.aberrationParameters.overwriteResultFiles){
    	    		if ((new File(resultPaths[nChn])).exists()) {
    	    			fileList[nChn]=null;
    					if (globalDebugLevel>0) System.out.println("File "+resultPaths[nChn]+" already exists and overwrite is disabled in configuration, channel "+nChn+" will be skipped");
    	    		}
    	    	}
    		}
    	}
	
		showDoubleFloatArrays sdfa_instance=new showDoubleFloatArrays();
	
		
		ImagePlus              impShow=new ImagePlus("CombinedKernels");              // just to show in the same window?
		long 	  startTime=System.nanoTime();
		for (int nChn=0; nChn<fileList.length;nChn++) if (fileList[nChn]!=null){
			// TODO: add parameters to kernel files			
			boolean OK=combinePSFKernels (
					interpolateParameters, // INTERPOLATE
					multiFilePSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
					fileList[nChn],
					resultPaths[nChn],
					sdfa_instance,        // SDFA_INSTANCE
					impShow, // just to show in the same window?
					saveResult,
					showResult,
					updateStatus,          // UPDATE_STATUS
					thisDebugLevel,
					globalDebugLevel);
			if (OK && (globalDebugLevel>0)) System.out.println("Saved combined kernel for channel "+nChn+" to"+resultPaths[nChn]+ " at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
    		if 	(stopRequested.get()>0) {				
				if (globalDebugLevel>0) System.out.println("User requested stop");
				break;
    		}

		}
		return true;
	}
	
	
	public boolean combinePSFKernels(
			EyesisAberrations.InterpolateParameters  interpolateParameters, // INTERPOLATE
			EyesisAberrations.MultiFilePSF           multiFilePSF ,         // MULTIFILE_PSF = new EyesisAberrations.MultiFilePSF(
			String []              filenames,
			String                 resultPath,
			showDoubleFloatArrays  sdfa_instance,        // SDFA_INSTANCE
			ImagePlus              imp_sel, // just to show in the same window?
			boolean                saveResult,
			boolean                showResult,
			boolean                updateStatus,          // UPDATE_STATUS
			int                    thisDebugLevel,
			int                    globalDebugLevel
	){	
		double [][][][] psfKernelMap; // will be lost - do we need it outside
		double [][][][][] kernelsElllipsePars = new double[filenames.length][][][][]; //x0,y0,a,b,c,area
		if (thisDebugLevel>0){
			System.out.println("combinePSFKernels(): filenames.length="+filenames.length);
		}
//		int i;
//		int nFile;
		int impProtoIndex=-1; // image index to copy all properties from (add combine? - i.e. 32/64 correlation)
		JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
		Opener opener=new Opener();
		for (int nFile=0;nFile<filenames.length;nFile++) {
			if (updateStatus) IJ.showStatus("Scanning file "+(nFile+1)+" (of "+(filenames.length)+"): "+filenames[nFile]);
			if (thisDebugLevel>1) System.out.println((nFile+1)+": "+filenames[nFile]);
			imp_sel=opener.openImage("", filenames[nFile]);  // or (path+filenames[nFile])
			// see if it has any usable properties and impProto is not set yet
			if (impProtoIndex<0){
				if ((imp_sel.getProperty("timestamp")==null) || (((String) imp_sel.getProperty("timestamp")).length()==0)) {
					jp4_instance.decodeProperiesFromInfo(imp_sel);
					if ((imp_sel.getProperty("timestamp")!=null) && (((String) imp_sel.getProperty("timestamp")).length()>0)) {
						impProtoIndex=nFile;
					}
				}
			}


			kernelsElllipsePars[nFile]= kernelStackToEllipseCoefficients( // null pointer
					imp_sel.getStack(), // Image stack, each slice consists of square kernels of one channel
					interpolateParameters.size, // size of each kernel (should be square)
					multiFilePSF.validateThreshold,
					globalDebugLevel);               //      threshold) // to find ellipse
		}

		// Visualize the array as stacks
		int nFiles=kernelsElllipsePars.length;
		int kHeight=kernelsElllipsePars[0].length;
		int kWidth=kernelsElllipsePars[0][0].length;
		int kLength=kHeight*kWidth;
		int nChn=imp_sel.getStack().getSize();
		int numResults=7;
		double [][][][] c= new double[numResults][nChn][nFiles+1][kLength];
		double [][][] numVals=new double[numResults][nChn][kLength];
//		int chn, tileY,tileX;
		boolean [] channels=new boolean[nChn];
		double a;
		if (thisDebugLevel>1) { 
			System.out.println("nFiles="+nFiles);
			System.out.println("kWidth="+kWidth);
			System.out.println("kHeight="+kHeight);
			System.out.println("nChn="+nChn);
		}
		Double D;
		int nOut;
		for (int chn=0;chn<nChn;chn++) {
			channels[chn]=false;
			for (int nFile=0;nFile<nFiles;nFile++) for (int tileY=0;tileY<kHeight;tileY++) for (int tileX=0;tileX<kWidth;tileX++) {
				//   			  System.out.println("nChn="+nChn+" nFile="+nFile+" tileY="+tileY+" tileX="+tileX);
				if (kernelsElllipsePars[nFile][tileY][tileX][chn]!=null) {
					channels[chn]=true;
					c[0][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][0]; // x0
					c[1][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][1]; // y0
					c[2][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][2]; // a
					c[3][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][3]; // b
					c[4][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][4]; // c
					a=1/Math.sqrt(kernelsElllipsePars[nFile][tileY][tileX][chn][2]*kernelsElllipsePars[nFile][tileY][tileX][chn][3]-
							kernelsElllipsePars[nFile][tileY][tileX][chn][4]*kernelsElllipsePars[nFile][tileY][tileX][chn][4]/4);
					c[5][chn][nFile+1][tileY*kWidth+tileX]= Math.sqrt(a); // radius
					c[6][chn][nFile+1][tileY*kWidth+tileX]=kernelsElllipsePars[nFile][tileY][tileX][chn][5]; // area

				} else {
					c[0][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[1][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[2][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[3][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[4][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[5][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
					c[6][chn][nFile+1][tileY*kWidth+tileX]=Double.NaN;
				}

			}
		}
		/* 
		 * Combine files - now just average all that are not NaN
		 */
		int [][] dirs={{-1,-1},{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1}};
//		int yn,xn,index;
		// remove any tiles that are not OK in all channels
		double [][] weights=new double[nFiles+1][kLength];
		for (int nFile=0;nFile<nFiles;nFile++) {
			for (int i=0;i<kLength;i++){
				weights[nFile+1][i]=1.0;
				for (int chn=0;chn<nChn;chn++) {
					D=c[0][chn][nFile+1][i];
					if (D.isNaN()) weights[nFile+1][i]=0.0;
				}
			}
			// Set weight to 0.5 if it has zero cells around        	
			for (int tileY=0;tileY<kHeight;tileY++) for (int tileX=0;tileX<kWidth;tileX++) {
				int index=tileY*kWidth+tileX;
				if ( weights[nFile+1][index]>0.0){
					for (int i=0;i<dirs.length;i++) {
						int yn=tileY+dirs[i][1];
						int xn=tileX+dirs[i][0];
						if ((yn>=0) && (yn<kHeight) && (xn>=0) && (xn<kWidth) && (weights[nFile+1][yn*kWidth+xn]==0.0)){
							weights[nFile+1][index]=0.5; // multiFilePSF.weightOnBorder; //0.5->0.01;
						}
					}
				}  
				weights[0][index]+=weights[nFile+1][index]; 
			}
		}
		if (thisDebugLevel>1) sdfa_instance.showArrays(weights, kWidth, kHeight,  true, "weights0");

		// remove any border ones if non-border is present in the same cell
		double [][] weightsMasked=new double[weights.length][];
		for (int i=0;i<weights.length;i++) weightsMasked[i]=weights[i].clone();
		double [][] weightsNotMasked=new double[weights.length][];
		for (int i=0;i<weights.length;i++) weightsNotMasked[i]=weights[i].clone();
		for (int i=0;i<kLength;i++){
			weightsMasked[0][i]=0.0;
			double maxWeight=0.0;
			for (int nFile=0;nFile<nFiles;nFile++)if (weightsMasked[nFile+1][i]>maxWeight) maxWeight=weightsMasked[nFile+1][i];
			for (int nFile=0;nFile<nFiles;nFile++) if ((weightsMasked[nFile+1][i]<1.0) && (maxWeight >= 1.0)) weightsMasked[nFile+1][i]=0.0; // do not count  half-weights if full one(s) are present
			for (int nFile=0;nFile<nFiles;nFile++)  weightsMasked[0][i]+=weightsMasked[nFile+1][i];
		}
		if (thisDebugLevel>1) sdfa_instance.showArrays(weightsMasked, kWidth, kHeight,  true, "weightsMasked");
		
		double [][][] psfRadius=c[5]; // later may remove all other calculations for c[i]?
		double [][][] pxfCenterX=c[0]; // some outlayer kernels have large x/y shift with normal radius - remove them too
		double [][][] pxfCenterY=c[1];
		if (thisDebugLevel>1) {
			for (int color=0;color<nChn;color++) sdfa_instance.showArrays(psfRadius[color], kWidth, kHeight,  true, "psfRadius-"+color);
			for (int color=0;color<nChn;color++) sdfa_instance.showArrays(pxfCenterX[color], kWidth, kHeight,  true, "pxfCenterX-"+color);
			for (int color=0;color<nChn;color++) sdfa_instance.showArrays(pxfCenterY[color], kWidth, kHeight,  true, "pxfCenterY-"+color);
		}
		double [][][] radiusRatio=new double[nChn][nFiles+1][kLength];
		for (int tileY=0;tileY<kHeight;tileY++) for (int tileX=0;tileX<kWidth;tileX++) {
			int index=tileY*kWidth+tileX;
			int totalNumSamples=0;
			for (int nFile=0;nFile<nFiles;nFile++) if ( weights[nFile+1][index]>0.0) totalNumSamples++;
			int samplesAfterWorse=totalNumSamples - ((int) Math.floor(multiFilePSF.maxFracDiscardWorse*totalNumSamples));
			int samplesAfterAll=  totalNumSamples - ((int) Math.floor(multiFilePSF.maxFracDiscardAll*totalNumSamples));

			int numSamples=totalNumSamples;
			
			while (numSamples>samplesAfterAll){ // calculate and remove worst sample until it is close enough to the average or too few samples are left
				boolean removeOnlyWorse=numSamples>samplesAfterWorse;
				for (int color=0;color<nChn;color++){
					radiusRatio[color][0][index]=0.0;
					pxfCenterX[color][0][index]=0.0; // same as c[0], zero before calculating average
					pxfCenterY[color][0][index]=0.0;
				}
				double sumWeights=0.0;
				for (int nFile=0;nFile<nFiles;nFile++) 	if ((weightsMasked[nFile+1][index]>0.0) && (weights[nFile+1][index]>0.0)){ // both, with outlayers removed
					for (int i=0;i<dirs.length;i++) {
						int yn=tileY+dirs[i][1];
						int xn=tileX+dirs[i][0];
						if ((yn>=0) && (yn<kHeight) && (xn>=0) && (xn<kWidth) && (weightsNotMasked[nFile+1][yn*kWidth+xn]>0.0)){ // including removed outlayers
							int indexNeib=xn+ kWidth*yn;
							double weight=weightsMasked[nFile+1][indexNeib];
							if (multiFilePSF.sharpBonusPower>0) {
								weight/=Math.pow(psfRadius[2][nFile+1][indexNeib],multiFilePSF.sharpBonusPower); // use green color
							}
							sumWeights+=weight;
							for (int color=0;color<nChn;color++) {
								radiusRatio[color][0][index]+=weight*psfRadius[color][nFile+1][indexNeib];
								pxfCenterX[color][0][index]+=weight*pxfCenterX[color][nFile+1][indexNeib];
								pxfCenterY[color][0][index]+=weight*pxfCenterY[color][nFile+1][indexNeib];
							}
						}
					}
				}
				/*
			System.out.println(tileY+":"+tileX+" - "+IJ.d2s(radiusRatio[0][0][index],3)+
					" "+IJ.d2s(radiusRatio[1][0][index],3)+
					" "+IJ.d2s(radiusRatio[2][0][index],3)+
					" sumWeights="+IJ.d2s(sumWeights,3));
				 */
				if (sumWeights>0.0) for (int color=0;color<nChn;color++) {
					radiusRatio[color][0][index]/=sumWeights; // average radius, without border-over-non-border cells 
					pxfCenterX[color][0][index]/=sumWeights; 
					pxfCenterY[color][0][index]/=sumWeights; 
				}
				double [] diffs=new double[nFiles];
				double [] diffsXY2=new double[nFiles];
//				for (int nFile=0;nFile<nFiles;nFile++) 	if ( weights[nFile+1][index]>0.0){ // here all , not just masked - why?
				for (int nFile=0;nFile<nFiles;nFile++) 	if ((weights[nFile+1][index]>0.0) && (weightsMasked[nFile+1][index]>0.0)){ //no outlayers, no masked - find worst
					diffs[nFile]=0;
					diffsXY2[nFile]=0;
					for (int color=0;color<nChn;color++) {
						radiusRatio[color][nFile+1][index]=psfRadius[color][nFile+1][index]/radiusRatio[color][0][index];
						double diff=(radiusRatio[color][nFile+1][index]>1.0)?(radiusRatio[color][nFile+1][index]-1.0):(1.0/radiusRatio[color][nFile+1][index]-1.0);
						if (diff>diffs[nFile]) diffs[nFile]=diff; // worst of 3 colors
						double diffX=pxfCenterX[color][nFile+1][index]-pxfCenterX[color][0][index];
						double diffY=pxfCenterY[color][nFile+1][index]-pxfCenterY[color][0][index];
						double diffXY2=diffX*diffX+diffY*diffY;
						if (diffXY2>diffsXY2[nFile]) diffsXY2[nFile]=diffXY2; // worst of 3 colors
						if (removeOnlyWorse && (psfRadius[color][nFile+1][index]<radiusRatio[color][0][index])) diffXY2=0.0; // only remove if radius is greater than average
						diffs[nFile]+=multiFilePSF.shiftToRadiusContrib*(Math.sqrt(diffXY2)/radiusRatio[color][0][index]); // now difference combines size and position
					}
				}
				// mask out outlayers
				//weightsMasked[0]

				// TODO: when averaging, divide by r^2 to some power to give bonus low-radius samples.
				//       also - combine dX^2+dY2+dR^2 when selecting outlayers
				if (thisDebugLevel>1) {
					System.out.print("\n === "+tileY+":"+tileX);
					for (int nFile=0;nFile<nFiles;nFile++) if ( weights[nFile+1][index]>0.0){
						System.out.print(" "+nFile);
					}
					System.out.println();
				}

				//			while (true){
//				int numSamples=0;
				numSamples=0;
				double worstDiff=0.0;
				int worstFile=0;
				for (int nFile=0;nFile<nFiles;nFile++) if ( weights[nFile+1][index]>0.0){
					numSamples++;
					int numNeib=0;
					for (int i=0;i<dirs.length;i++) {
						int yn=tileY+dirs[i][1];
						int xn=tileX+dirs[i][0];
						if ((yn>=0) && (yn<kHeight) && (xn>=0) && (xn<kWidth) && (weights[nFile+1][yn*kWidth+xn]>0.0))numNeib++;
					}
					double scale= 8.0/(8.0+  multiFilePSF.internalBonus*numNeib); // make internal cells look better
					double diff=diffs[nFile]*scale;
					if (diff>worstDiff){
						worstDiff=diff;
						worstFile=nFile;
					}
					if (thisDebugLevel>2) {
						System.out.println(tileY+":"+tileX+" - "+
								" nFile="+nFile+
								" scale=="+IJ.d2s(scale,3)+
								" diff="+diff+
								" numNeib="+numNeib);
						}

				}
				if (thisDebugLevel>1) {
				System.out.println(tileY+":"+tileX+" - "+
						" numSamples="+numSamples+
						" worstDiff="+IJ.d2s(worstDiff,3)+
						" removeOnlyWorse="+removeOnlyWorse+
						" worstFile="+worstFile+" ("+filenames[worstFile]+")");
				
				System.out.println(
						" this radius  ={"+psfRadius[0][worstFile+1][index]+","+psfRadius[1][worstFile+1][index]+","+psfRadius[2][worstFile+1][index]+"}\n"+
						" mean radius  ={"+radiusRatio[0][0][index]+","+radiusRatio[1][0][index]+","+radiusRatio[2][0][index]+"}\n"+
						
						" this center X={"+pxfCenterX[0][worstFile+1][index]+","+pxfCenterX[1][worstFile+1][index]+","+pxfCenterX[2][worstFile+1][index]+"}\n"+
						" mean center X={"+pxfCenterX[0][0][index]+","+pxfCenterX[1][0][index]+","+pxfCenterX[2][0][index]+"}\n"+
						
						" this center Y={"+pxfCenterY[0][worstFile+1][index]+","+pxfCenterY[1][worstFile+1][index]+","+pxfCenterY[2][worstFile+1][index]+"}\n"+
						" mean center Y={"+pxfCenterY[0][0][index]+","+pxfCenterY[1][0][index]+","+pxfCenterY[2][0][index]+"}");
				}
				if (numSamples<1) break; // nothing left 
				if (worstDiff>multiFilePSF.radiusDiffHigh){
					if ( (numSamples==1) && (globalDebugLevel>0)){
						System.out.println("PSF size for the cell "+tileX+":"+tileY+", file# "+worstFile+" varies too much from the neighbor cells, so it is removed, creating a gap");
					}
					weights[worstFile+1][index]=0.0;
					continue;
				} else if ((worstDiff>multiFilePSF.radiusDiffLow) && (numSamples>1)){
					weights[worstFile+1][index]=0.0;
					continue;
				}
				break;
			}
			// recalculate sum of weights;
			weights[0][index]=0.0;
			for (int nFile=0;nFile<nFiles;nFile++) if ( weights[nFile+1][index]>0.0){
				weights[0][index]+=weights[nFile+1][index];
			}
		}
		
		
		// for each channel, each cell - compare radius calculated for neighbors (use masked weights) and the cell
		
		
		// TODO: Filter out outlayers: Add bonus to cells surrounded by others?
		
		
		//    	double [][][][] c= new double[numResults][nChn][nFiles+1][kLength];
		//     	double [][][] numVals=new double[numResults][nChn][kLength];
		for (int chn=0;chn<nChn;chn++) {

			for (nOut=0;nOut<c.length;nOut++) {
				c[nOut][chn][0]=null;
				for (int i=0;i<kLength;i++) {
					numVals[nOut][chn][i]=0.0;
				}
			}
			if (channels[chn]) {
				for (nOut=0;nOut<c.length;nOut++) {
					c[nOut][chn][0]=new double [kLength];
					for (int nFile=0;nFile<nFiles;nFile++) {
						for (int i=0;i<kLength;i++){
							D=c[nOut][chn][nFile+1][i];
							if (!D.isNaN()){
								numVals[nOut][chn][i]+=1.0;
								c[nOut][chn][0][i]+=D*weights[nFile+1][i]/weights[0][i];
							}
						}

					}
					for (int i=0;i<kLength;i++){
						if (numVals[nOut][chn][i]==0.0 )c[nOut][chn][0][i]=Double.NaN;
						//    			  else c[nOut][chn][0][i]/=numVals[nOut][chn][i];
					}       	    	
				}
				
				if (multiFilePSF.validateShowEllipse) {
						sdfa_instance.showArrays(radiusRatio[chn], kWidth, kHeight,  true, "ratio-"+chn);
						sdfa_instance.showArrays(c[5][chn],kWidth, kHeight,  true, "radius-"+chn);
					if (thisDebugLevel>1) {
						sdfa_instance.showArrays(c[0][chn], kWidth, kHeight,  true, "x-shift-"+chn);
						sdfa_instance.showArrays(c[1][chn], kWidth, kHeight,  true, "y-shift-"+chn);
						sdfa_instance.showArrays(c[2][chn], kWidth, kHeight,  true, "x2-"+chn);
						sdfa_instance.showArrays(c[3][chn], kWidth, kHeight,  true, "y2-"+chn);
						sdfa_instance.showArrays(c[4][chn], kWidth, kHeight,  true, "xy-"+chn);
						sdfa_instance.showArrays(c[6][chn], kWidth, kHeight,  true, "area-"+chn);
					}  
				}
			}

		}
		if (multiFilePSF.showWeights) sdfa_instance.showArrays(weights, kWidth, kHeight,  true, "weights");
		//    	double [][] weights=new double[nFiles+1][kLength];
		for (int i=0;i<kLength;i++) weights[0][i]=0.0;
		psfKernelMap=new double [kHeight][kWidth][nChn][];
		for (int tileY=0;tileY<kHeight;tileY++) for (int tileX=0;tileX<kWidth;tileX++) for (int chn=0;chn<nChn;chn++){
			psfKernelMap[tileY][tileX][chn]=null;
		}
		String [] originalSliceLabels=null;
		for (int nFile=0;nFile<nFiles;nFile++) {
			if (updateStatus) IJ.showStatus("Accumulating file "+(nFile+1)+" (of "+nFiles+"): "+filenames[nFile]);
			if (thisDebugLevel>1) System.out.println("Accumulating file "+nFile+": "+filenames[nFile]);
			imp_sel=opener.openImage("", filenames[nFile]);  // or (path+filenames[nFile])
			if (originalSliceLabels==null) {
				originalSliceLabels=imp_sel.getStack().getSliceLabels();
			}
			accumulatePartialKernelStack(
					psfKernelMap,
					imp_sel.getStack(), // Image stack with partial array of kernels, each slice consists of square kernels of one channel
					interpolateParameters.size, // size of each kernel (should be square)
					weights[nFile+1], // weights of the kernel tiles in the current stack
					weights[0],
					globalDebugLevel);// weights of the kernel tiles already accumulated (will be updated)

		}
// optionally fill in blanks from nearest neighbors
		int filledMissing=0;
		//Finalize accumulated kernels - transform them from frequency to space domain
		inverseTransformKernels(psfKernelMap);
// should be done after inversion, because filled in kernels are just pointers to original ones		
		if (multiFilePSF.fillMissing) filledMissing=fillMissingKernels (psfKernelMap);
		int numMissing=0;
		ImageStack mergedStack= mergeKernelsToStack(psfKernelMap,originalSliceLabels);
		System.out.println("mergedStack.getSize()= "+mergedStack.getSize());
		System.out.println("mergedStack.getWidth()= "+mergedStack.getWidth()  );
		System.out.println("mergedStack.getHeight()= "+mergedStack.getHeight()  );
		System.out.println("psfKernelMaplength= "+psfKernelMap.length  );
		System.out.println("psfKernelMap[0].length= "+psfKernelMap[0].length  );
		System.out.println("mergedStack= "+((mergedStack==null)?"null":"not null"));

		if (mergedStack.getSize()==0) {
			System.out.println("*** Error - result is empty");
			return false;
		}

		for (int tileY=0;tileY<kHeight;tileY++) for (int tileX=0;tileX<kWidth;tileX++) if ((psfKernelMap[tileY][tileX]==null) || (psfKernelMap[tileY][tileX][0]==null)) numMissing++;
        ImagePlus imp_psf = new ImagePlus(resultPath, mergedStack);
        
        if (impProtoIndex>=0){
			imp_sel=opener.openImage("", filenames[impProtoIndex]); 
			jp4_instance.decodeProperiesFromInfo(imp_sel);
			// copy properties from the source image
			jp4_instance.copyProperties (imp_sel,imp_psf);
        }
        multiFilePSF.setProperties("MULTIFILE_PSF.", imp_psf);
  // other properties
		jp4_instance.encodeProperiesToInfo(imp_psf);
        if (showResult) {
        	imp_psf.getProcessor().resetMinAndMax();
        	imp_psf.show();
        }
		if (saveResult) {
			if (numMissing==0) {
			  if (thisDebugLevel>1) System.out.println("Saving result to "+resultPath);
			  FileSaver fs=new FileSaver(imp_psf);
			  fs.saveAsTiffStack(resultPath);
			  if (multiFilePSF.fillMissing && (filledMissing>0)) {
					System.out.println("*** Warning "+filledMissing+" kernel tiles are missing from the results (insufficient overlap) \n"+
					"You may disable filling missing kernels from neighbors in Conf. Multifile");
/*
  					IJ.showMessage("Warning",filledMissing+" kernel tiles were missing from the results\n"+
 							"(i.e.insufficient overlap) and filled from neighbors, it is OK only for the fisheye lens.\n"+
					        "You may disable filling missing kernels from neighbors in Conf. Multifile");*/
			  }
			  return true;
			} else {
				System.out.println("*** Error "+numMissing+" kernel tiles are missing from the results (insufficient overlap), result is not saved\n"+
				"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");

				IJ.showMessage("Error",numMissing+" kernel tiles are missing from the results\n (insufficient overlap), result file is not saved\n"+
				"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");

				if (!showResult) { // not yet shown
		        	imp_psf.getProcessor().resetMinAndMax();
		        	imp_psf.show();
				}
				return false;
			}
		}
		if (numMissing>0) {
			System.out.println("*** Error "+numMissing+" kernel tiles are missing from the results (insufficient overlap) \n"+
					"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");
			IJ.showMessage("Error",numMissing+" kernel tiles are missing from the results\n (insufficient overlap)\n"+
					"You may enable filling missing kernels from neighbors if it is a fisheye lens (in Conf. Multifile)");
			return false;

		} else if (multiFilePSF.fillMissing && (filledMissing>0)) {
			System.out.println("*** Warning "+filledMissing+" kernel tiles are missing from the results (insufficient overlap) \n"+
			"You may disable filling missing kernels from neighbors in Conf. Multifile");
			IJ.showMessage("Warning",filledMissing+" kernel tiles were missing from the results\n"+
					"(i.e.insufficient overlap) and filled from neighbors, it is OK only for the fisheye lens.\n"+
			        "You may disable filling missing kernels from neighbors in Conf. Multifile");
		}
		return true;
	}
	
	private int fillMissingKernels(double [][][][] kernels){
		int [][] dirs={{-1,0},{1,0},{0,-1},{0,1}};
		List <Integer> kernelList=new ArrayList<Integer>(100);
		Integer Index;
		kernelList.clear();
		int tileY,tileX,newTileY,newTileX,nDir,numMissing=0;
		int width= kernels[0].length;
		int height=kernels.length;
		for (tileY=0;tileY<height;tileY++) for (tileX=0;tileX<width;tileX++) {
			if ((kernels[tileY][tileX]==null) || (kernels[tileY][tileX][0]==null)) {
				Index=tileY*width+tileX;
				for (nDir=0;nDir<dirs.length;nDir++) {
					newTileX=tileX+dirs[nDir][0];
					newTileY=tileY+dirs[nDir][1];
					if ((newTileX>=0) && (newTileY>=0) && (newTileX<width) && (newTileY<height) &&
							(kernels[newTileY][newTileX]!=null) && (kernels[newTileY][newTileX][0]!=null) ) {
						kernelList.add(Index);
					}
				}				
				numMissing++;
			}
		}
		System.out.println("fillMissingKernels: numMissing="+numMissing);
		System.out.println("fillMissingKernels: kernelList.size()="+kernelList.size());

		while (kernelList.size()>0) {
			Index=kernelList.get(0);
			kernelList.remove(0);
			tileY=Index/width;
			tileX=Index%width;
			if ((kernels[tileY][tileX]==null) || (kernels[tileY][tileX][0]==null)) {// may be duplicates (added several times)
//TODO: - change order of directions?
				for (nDir=0;nDir<dirs.length;nDir++) {
					newTileX=tileX+dirs[nDir][0];
					newTileY=tileY+dirs[nDir][1];
					if ((newTileX>=0) && (newTileY>=0) && (newTileX<width) && (newTileY<height)) {
						if ((kernels[newTileY][newTileX]==null) || (kernels[newTileY][newTileX][0]==null)) {
							Index=newTileY*width+newTileX;
							kernelList.add(Index);
						} else if ((kernels[tileY][tileX]==null) || (kernels[tileY][tileX][0]==null)) { // may be already added
// need to copy - they will be subject to reverse fht	?						
							kernels[tileY][tileX]=kernels[newTileY][newTileX];
							System.out.println("fillMissingKernels: filled "+tileX+"/"+tileY);
						}
					}
				}
			}
		}
		return numMissing;
	}
	
	
	
	/* ======================================================================== */
	//Finalize accumulated kernels - transform them from frequency to space domain
	public void inverseTransformKernels(
			double [][][][] psfKernelMap){
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		int tilesX=psfKernelMap[0].length;
		int tilesY=psfKernelMap.length;
		int tileY,tileX, chn; //,subTileY,subTileX;
		for (tileY=0;tileY<tilesY;tileY++) for (tileX=0;tileX<tilesX;tileX++) for (chn=0; chn<psfKernelMap[tileY][tileX].length; chn++){
			if (psfKernelMap[tileY][tileX][chn]!=null) {
				fht_instance.inverseTransform(psfKernelMap[tileY][tileX][chn]);
				fht_instance.swapQuadrants   (psfKernelMap[tileY][tileX][chn]);
			}
		}
	}


	
	
	
	// Will build global PSF_KERNEL_MAP (each [][][]element should be set to null?
	// kernels are supposed to be normalized?
	public void accumulatePartialKernelStack(
			double [][][][] psfKernelMap,
			ImageStack   kernelStack, // Image stack with partial array of kernels, each slice consists of square kernels of one channel
			int                 size, // size of each kernel (should be square)
			double []   theseWeights, // weights of the kernel tiles in the current stack
			double []   accumWeights,// weights of the kernel tiles already accumulated (will be updated)
			int debugLevel){
		DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		int tilesX=kernelStack.getWidth()/size;
		int tilesY=kernelStack.getHeight()/size;
		int nChn=kernelStack.getSize();
		int tileY,tileX, chn; //,subTileY,subTileX;
		float [] pixels;
		int length=size*size;
		double [] kernel=new double[length];
		int index;
		int debugTileX=18;
		int debugTileY=22;
		int debugIndex=debugTileY*tilesX+debugTileX;
		boolean lastChn;
		for (chn=0;chn<nChn;chn++) {
			pixels=(float[]) kernelStack.getPixels(chn+1);
			lastChn= (chn==(nChn-1));
			for (tileY=0;tileY<tilesY;tileY++) for (tileX=0;tileX<tilesX;tileX++) {
				index=tileY*tilesX+tileX;
				boolean debugThis= (index==debugIndex) && (debugLevel>1);
				if (theseWeights[index]>0.0){
					extractOneKernel(
							pixels, //  array of combined square kernels, each 
							kernel, // will be filled, should have correct size before call
							tilesX, // number of kernels in a row
							tileX, // horizontal number of kernel to extract
							tileY); // vertical number of kernel to extract
					// convert to frequency domain (interpolation is for FHT)
					fht_instance.swapQuadrants(kernel);
					fht_instance.transform(    kernel);
					if (debugThis) System.out.println("tileY="+tileY+" tileX="+tileX+" chn= "+chn+
							" theseWeights["+index+"]="+theseWeights[index]+" accumWeights["+index+"]="+accumWeights[index]);
					if (!(accumWeights[index]>0.0)) { // nothing yet in this tile
						psfKernelMap[tileY][tileX][chn]=kernel.clone();
						if (lastChn) accumWeights[index]=theseWeights[index];
					} else { // "accumulate" - interpolate between existent and new kernel, using/updating weights
						if (debugLevel>5) {
							System.out.println("tileY="+tileY+" tileX="+tileX+" chn= "+chn);
							System.out.println("PSF_KERNEL_MAP[tileY][tileX][chn].length= "+psfKernelMap[tileY][tileX][chn].length);
							System.out.println("kernel.length= "+kernel.length);
						}

//						kernel=fht_instance.interpolateFHT (
						psfKernelMap[tileY][tileX][chn]=fht_instance.interpolateFHT (
								psfKernelMap[tileY][tileX][chn],    // first FHT array
								kernel,    // second FHT array
								theseWeights[index]/accumWeights[index]);    //interpolation ratio - 0.0 - fht0, 1.0 - fht1
						if (lastChn) accumWeights[index]+=theseWeights[index];
					}
				}
			}
		}
	}

	
	
	public double [][][][] kernelStackToEllipseCoefficients(
			ImageStack kernelStack, // Image stack, each slice consists of square kernels of one channel
			int               size, // size of each kernel (should be square)
			double       threshold,
			int         debugLevel) // to find ellipse
	// update status info
	{
		//	  DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
		if (kernelStack==null) return null;
		int tilesX=kernelStack.getWidth()/size;
		int tilesY=kernelStack.getHeight()/size;
		int nChn=kernelStack.getSize();
		float [] pixels;
		int i,j;
		int tileY,tileX, chn; //,subTileY,subTileX;
		double [][][][] ellipseCoeffs=new double [tilesY][tilesX][nChn][];
		int length=size*size;
		double [] kernel=new double[length]; 
		double max;
		int  [][]selection;
		double [] ec;
		int l;
		for (chn=0;chn<nChn;chn++) {
			pixels=(float[]) kernelStack.getPixels(chn+1);
			for (tileY=0;tileY<tilesY;tileY++) for (tileX=0;tileX<tilesX;tileX++) {
				extractOneKernel(
						pixels, //  array of combined square kernels, each 
						kernel, // will be filled, should have correct size before call
						tilesX, // number of kernels in a row
						tileX, // horizontal number of kernel to extract
						tileY); // vertical number of kernel to extract
				max=0.0;
				for (i=0;i<length;i++) if (max<kernel[i]) max=kernel[i];
				if (max<=0.0) ellipseCoeffs[tileY][tileX][chn]=null;
				else {
					selection= findClusterOnPSF(
							kernel, // PSF function, square array
							threshold, // fraction of energy in the pixels to be used
					"",
					debugLevel);
					//				  ellipseCoeffs[tileY][tileX][chn]=findEllipseOnPSF(kernel,  selection,   "");
					ec=findEllipseOnPSF(kernel,  selection,   "", debugLevel); // x0,y0,a,b,c (r2= a* x^2*+b*y^2+c*x*y)

					l=ec.length;
					ellipseCoeffs[tileY][tileX][chn]=new double[l+1];
					for (i=0;i<ec.length;i++) ellipseCoeffs[tileY][tileX][chn][i]=ec[i];
					ellipseCoeffs[tileY][tileX][chn][l]=0;
					for (i=0;i<selection.length;i++) for (j=0;j<selection[0].length;j++) ellipseCoeffs[tileY][tileX][chn][l]+=selection[i][j];
				}
			}
		}
		return ellipseCoeffs;

	} 
	
	/* ======================================================================== */
	private void extractOneKernel(
			float []  pixels, //  array of combined square kernels, each 
			double [] kernel, // will be filled, should have correct size before call
			int       numHor, // number of kernels in a row
			int        xTile, // horizontal number of kernel to extract
			int        yTile) { // vertical number of kernel to extract
		int length=kernel.length;
		int size=(int) Math.sqrt(length);
		int i,j;
		int pixelsWidth=numHor*size;
		int pixelsHeight=pixels.length/pixelsWidth;
		int numVert=pixelsHeight/size;
/* limit tile numbers - effectively add margins around the known kernels */
		if (xTile<0) xTile=0;
		else if (xTile>=numHor) xTile=numHor-1;
		if (yTile<0) yTile=0;
		else if (yTile>=numVert) yTile=numVert-1;
		int base=(yTile*pixelsWidth+xTile)*size;
		for (i=0;i<size;i++) for (j=0;j<size;j++) kernel [i*size+j]=pixels[base+i*pixelsWidth+j];
	}

	
	
	
	
	
//=======================================================	
	public void savePartialKernelStack(
			String path,
			ImageStack stack,
			ImagePlus impSrc, // properties - decoded
			PSFParameters psfParameters,
			boolean [] correlationSizesUsed
	){
		int [] corrSizes={};
		if (correlationSizesUsed!=null) {
			int numDifferentFFT=0;
			for (int i=0;i<correlationSizesUsed.length;i++) if (correlationSizesUsed[i]) {
				numDifferentFFT++;
			}
			corrSizes=new int [numDifferentFFT];
			int index=0;
			for (int i=0;i<correlationSizesUsed.length;i++) if (correlationSizesUsed[i]) {
				corrSizes[index++]=1<<i;
			}
		}
		ImagePlus impPsf = new ImagePlus(path, stack);

		JP46_Reader_camera jp4_instance= new JP46_Reader_camera(false);
		if ((impSrc.getProperty("timestamp")==null) || (((String) impSrc.getProperty("timestamp")).length()==0)) {
			jp4_instance.decodeProperiesFromInfo(impSrc);
		}
		// copy properies from the source image
		jp4_instance.copyProperties (impSrc,impPsf);
		// save psf parameters (same as in configuration file)
		psfParameters.setProperties("PSF_PARS.", impPsf);
		for (int i=0;i<corrSizes.length;i++){
			impPsf.setProperty("corr_size_"+corrSizes[i], true+"");
		}
//TODO:  Add more properties?

		jp4_instance.encodeProperiesToInfo(impPsf);
		FileSaver fs=new FileSaver(impPsf);
		fs.saveAsTiffStack(path);
	}
	
	
	
	private  ImageStack mergeKernelsToStack(double [][][][] kernels) {
		return mergeKernelsToStack(kernels,null);

	}
	private  ImageStack mergeKernelsToStack(double [][][][] kernels,String [] names) { // use oldStack.getSliceLabels() to get names[]
		if (kernels==null) return null;
		int tilesY=kernels.length;
		int tilesX=kernels[0].length;
		int i,j,k,nChn, chn,x,y,index;
		double [][]kernel=null;
		for (i=0;(i<tilesY) && (kernel==null);i++)  for (j=0;(j<tilesX) && (kernel==null);j++)  kernel=kernels[i][j];
		if (kernel==null) return null;
		int length=0;
		for (i=0;i<kernel.length;i++) if (kernel[i]!=null){
			length=kernel[i].length;
			break;
		}
		if (length==0){
			System.out.println("mergeKernelsToStack(): no non-null kernels");
			return null;
		}
		int [] channelsMask = new int [kernel.length];
		for (i=0;i<kernel.length;i++) channelsMask[i]=0;
		for (i=0;i<tilesY ;i++)  for (j=0;j<tilesX;j++) if (kernels[i][j]!=null) {
			for (k=0;(k<kernel.length)&& (k<channelsMask.length);k++) if (kernels[i][j][k]!=null) {
				channelsMask[k]=1;
				if (kernels[i][j][k].length>length) length=kernels[i][j][k].length;
			}
		}

		nChn=0;
		for (i=0;i<channelsMask.length;i++) if (channelsMask[i]!=0) nChn++;
		int [] channels = new int [nChn];
		nChn=0;
		for (i=0;i<channelsMask.length;i++) if (channelsMask[i]!=0) channels[nChn++]=i;

		//	    for (i=0;i<kernel.length;i++) if (kernel[i]!=null)  if (nChn<channels.length) channels[nChn++]=i;
		int size=(int) Math.sqrt(length);
		int outWidth= size*tilesX;
		int outHeight=size*tilesY;

		ImageStack stack=new ImageStack(outWidth,outHeight);
		int numKernels=0;
		float [] fpixels;
		for (chn=0;chn<nChn;chn++) {
			fpixels= new float [outWidth*outHeight];
			k=channels[chn];
			for (i=0; i<tilesY;i++)  for (j=0;j<tilesX;j++) {
				for (y=0;y<size;y++) for (x=0;x<size;x++) {
					index=((i*size+y)*outWidth)+(j*size+x);
					if ((kernels[i][j]==null || (kernels[i][j][k]==null))) fpixels[index]=0.0f;
					else {
						fpixels[index]= (float) kernels[i][j][k][y*size+x];
					}
				}
				if ((kernels[i][j]!=null && (kernels[i][j][k]!=null))) numKernels++;
			}
			if (names==null) stack.addSlice("channel"+k, fpixels);
			else             stack.addSlice(names[chn], fpixels);
		}
		if (numKernels==0){
			System.out.println("mergeKernelsToStack(): all kernels are empty");
			return null;
		}
		System.out.println("mergeKernelsToStack(): got "+numKernels +" non-null kernels");
		return stack;
	}


	public double [][][][] createPSFMap(
			final MatchSimulatedPattern commonMatchSimulatedPattern, // to be cloned in threads, using common data
			final ImagePlus          imp_sel, // linearized Bayer mosaic image from the camera, GR/BG
			final int [][][]        sampleList, // optional (or null) 2-d array: list of coordinate pairs (2d - to match existent  pdfKernelMap structure)  
			final double  overexposedAllowed, // fraction of pixels OK to be overexposed
			final SimulationPattern.SimulParameters simulParameters,
			final int             mapFFTsize, // scanImageForPatterns:FFT size
			final MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			final int            fft_overlap,
			final int               fft_size,
			final ColorComponents colorComponents,
			final int           PSF_subpixel, 
			final OTFFilterParameters otfFilterParameters,
			final PSFParameters psfParameters,
			final double       minDefinedArea,
			final int          PSFKernelSize, // size of square used in the new map (should be multiple of map step)
			final double       gaussWidth,  // ** NEW
			final float [][]   simArray,    // ** NEW
			final int          threadsMax,
			final boolean      updateStatus,          // UPDATE_STATUS
			final int          masterDebugLevel, // get rid of it? // ** NEW
			final int          globalDebugLevel,// ** NEW
			final int          debug_level){// debug level used inside loops
		final boolean debugLateralShifts=(globalDebugLevel>1);
		System.out.println("createPSFMap(): masterDebugLevel="+masterDebugLevel+" globalDebugLevel="+globalDebugLevel+" debug_level="+debug_level); // 2 2 0
		final long startTime = System.nanoTime();
		  Runtime runtime = Runtime.getRuntime();
		  runtime.gc();
		  if (globalDebugLevel>1) System.out.println("--- Free memory="+runtime.freeMemory()+" (of "+runtime.totalMemory()+")");

		// Generate hi-res pattern bitmap (one cell)	
		SimulationPattern simulationPattern= new SimulationPattern();
		simulationPattern.debugLevel=globalDebugLevel;
		final double [] bitmaskPattern= simulationPattern.patternGenerator(simulParameters);
		int nTileX,nTileY;
		int numPatternCells=0;
/* Filter results based on correlation with the actual pattern */
		boolean [][]   PSFBooleanMap; // map of 2*fft_size x 2*fft_size squares with 2*fft_overlap step, true means that that tile can be used for PSF
		if (sampleList==null){
				PSFBooleanMap= mapFromPatternMask ( // count number of defined cells
						commonMatchSimulatedPattern,
						imp_sel.getWidth(), // image (mask) width
						fft_size*2,
						fft_overlap*2,
						fft_size,   // backward compatibility margin==tileSize/2
						gaussWidth,
						//	psfParameters.minDefinedArea);
						minDefinedArea,
						globalDebugLevel);
		} else {
			PSFBooleanMap= new boolean[sampleList.length][sampleList[0].length];
			for (int i=0;i<sampleList.length;i++) for (int j=0;j<sampleList[0].length;j++) PSFBooleanMap[i][j]=(sampleList[i][j][0]>=0); // all with positive X
		}
		if (PSFBooleanMap==null) return null;
		numPatternCells=0;
		for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++) if (PSFBooleanMap[nTileY][nTileX]) numPatternCells++;
		if (globalDebugLevel>1) {
			System.out.println("Remapped for PSF measurment- converted to an array["+PSFBooleanMap.length+"]["+PSFBooleanMap[0].length+"], "+
					numPatternCells+" cells (of "+(PSFBooleanMap.length*PSFBooleanMap[0].length)+") with pattern detected");
		}
		pdfKernelMap=new double[PSFBooleanMap.length][PSFBooleanMap[0].length][][]; //pdfKernelMap - global (or final)
//		int saved_globalDebugLevel=globalDebugLevel;
//		globalDebugLevel=debug_level;
		simulationPattern.debugLevel=globalDebugLevel;
		int ncell=0;
/* Create array of coordinates of cells to process, fill result array with zeros (to be actually written by threads */       
		final int [][] tilesToProcessXY=new int [numPatternCells][4];

		for (nTileY=0;nTileY<PSFBooleanMap.length;nTileY++) for (nTileX=0;nTileX<PSFBooleanMap[0].length;nTileX++){
			if (PSFBooleanMap[nTileY][nTileX]) {
				tilesToProcessXY[ncell  ][0]=nTileX;
				tilesToProcessXY[ncell  ][1]=nTileY;
				tilesToProcessXY[ncell  ][2]=(sampleList==null)?(fft_overlap*2*nTileX):sampleList[nTileY][nTileX][0];
				tilesToProcessXY[ncell++][3]=(sampleList==null)?(fft_overlap*2*nTileY):sampleList[nTileY][nTileX][1];
				pdfKernelMap[nTileY][nTileX]=new double[colorComponents.colorsToCorrect.length][];
			} else pdfKernelMap[nTileY][nTileX]=null;
		}
		final double [][][][] debugLateral=new double [PSFBooleanMap.length][PSFBooleanMap[0].length][][];
		for (nTileY=0;nTileY<debugLateral.length;nTileY++) for (nTileX=0;nTileX<debugLateral[0].length;nTileX++){
			debugLateral[nTileY][nTileX]=null;
		}
		final Thread[] threads = newThreadArray(threadsMax);
		if (globalDebugLevel>1) System.out.println("Starting "+threads.length+" threads: "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		final AtomicInteger ai = new AtomicInteger(0);
		final int patternCells=numPatternCells;
		//	  final double []   overexposedMap, // map of overexposed pixels in the image (may be null)
		final double [] overexposed=(overexposedAllowed>0)?JP4_INSTANCE.overexposedMap (imp_sel):null;
		final int mapWidth=imp_sel.getWidth();
   		final AtomicInteger tilesFinishedAtomic = new AtomicInteger(1); // first finished will be 1
   		final int debugNumColors=6; 
		for (int ithread = 0; ithread < threads.length; ithread++) {
			// Concurrently run in as many threads as CPUs
			threads[ithread] = new Thread() {
				public void run() {

					// Each thread processes a few items in the total list
					// Each loop iteration within the run method has a unique 'i' number to work with
					// and to use as index in the results array:
					//	double [] sum_kern_el=new double[6]; // just testing					
					int x0,y0,nTX,nTY,nChn;
					double [][] kernels;
					// change to true (first 2 only?) to separate memory arrays for threads
				    MatchSimulatedPattern matchSimulatedPattern=commonMatchSimulatedPattern.cloneDeep(
				    		false, // boolean clonePATTERN_GRID,
				    		false, // boolean cloneTargetUV,
				    		false, // boolean clonePixelsUV,
				    		false, // boolean cloneFlatFieldForGrid,
				    		false  // boolean cloneFocusMask
				    		);
				    matchSimulatedPattern.debugLevel=globalDebugLevel;
					SimulationPattern simulationPattern= new SimulationPattern(bitmaskPattern);
					simulationPattern.debugLevel=globalDebugLevel;
					double [] windowFFTSize=    matchSimulatedPattern.initWindowFunction(fft_size,gaussWidth); //=initHamming( fft_size) calculate once
					double [] windowFullFFTSize=matchSimulatedPattern.initWindowFunction(fft_size*PSF_subpixel,gaussWidth); //=initHamming( fft_size*subpixel);
					DoubleFHT fht_instance =new DoubleFHT(); // provide DoubleFHT instance to save on initializations (or null)
					double over;
// individual per-thread - will be needed when converted to doubleFHT					
//				    MatchSimulatedPattern matchSimulatedPattern=new MatchSimulatedPattern(FFT_SIZE);
					for (int nTile = ai.getAndIncrement(); nTile < patternCells; nTile = ai.getAndIncrement()) {
						nTX=tilesToProcessXY[nTile][0];
						nTY=tilesToProcessXY[nTile][1];
						y0=tilesToProcessXY[nTile][3];
						x0=tilesToProcessXY[nTile][2];
						if (updateStatus) IJ.showStatus("Processing tile["+nTY+"]["+nTX+"] ("+(nTile+1)+" of "+patternCells+")");
						if (masterDebugLevel>1) System.out.println("#!# "+x0+":"+y0+" Processing tile["+nTY+"]["+nTX+"] ("+(nTile+1)+" of "+patternCells+") : "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
						if (overexposed!=null){
							over=JP4_INSTANCE.fracOverExposed(overexposed,   // map of overexposed pixels 0.0 - 0K, >0 (==1.0) - overexposed
									mapWidth,    // width of the map
									x0,          // X of the top left corner of the selection
									y0,          // Y of the top left corner of the selection
									2*fft_size,  // selection width
									2*fft_size); // selection height
							//						  if ((globalDebugLevel>0) && (over>0.0)) System.out.println("Overexposed fraction of "+over+" at x0="+x0+" y0="+y0+" width"+(2*fft_size));
						} else over=-1.0;
						if ( over > overexposedAllowed) {
							pdfKernelMap[nTY][nTX]=null;
							if (globalDebugLevel>0) System.out.println("Overexposed fraction of "+over+" at x0="+x0+" y0="+y0+" width"+(2*fft_size));
						} else {
							if (debugLateralShifts) {
								debugLateral[nTY][nTX]= new double [debugNumColors][];    // X/Y shift of the PSF array, in Bayer component pixel coordinates (same as PSF arrays)
							} 
							kernels=getPSFKernels(imp_sel,
									simArray, //simulation image, scaled PSF_subpixel/2
									2*fft_size,       // size in pixels (twice fft_size)
									x0,               // top left corner X (pixels)
									y0,               // top left corner Y (pixels)
									simulationPattern,
									matchSimulatedPattern,
									patternDetectParameters,
									windowFFTSize,    //=initHamming( fft_size) calculate once
									windowFullFFTSize,//=initHamming( fft_size*subpixel);
									PSF_subpixel,     // use finer grid than actual pixels 
									simulParameters,
									colorComponents,  // color channels to process, equalizeGreens
									otfFilterParameters,
									5,                // int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
									psfParameters,
									fht_instance,      // provide DoubleFHT instance to save on initializations (or null)
									debug_level,// ((x0<512)&& (y0<512))?3:debug_level DEBUG during "focusing"
									masterDebugLevel, // get rid of it? // ** NEW
									globalDebugLevel,// ** NEW
									debugLateralShifts?debugLateral[nTY][nTX]:null
							);
							if (kernels!=null) {
								if (kernelLength(kernels)>(PSFKernelSize*PSFKernelSize)) kernels=resizeForFFT(kernels,PSFKernelSize); // shrink before normalizing
								normalizeKernel(kernels); // in-place
								if (kernelLength(kernels)<(PSFKernelSize*PSFKernelSize)) kernels=resizeForFFT(kernels,PSFKernelSize); // expand after normalizing
								for (nChn=0;nChn<kernels.length;nChn++) if (kernels[nChn]!=null){
									pdfKernelMap[nTY][nTX][nChn]=kernels[nChn]; // not .clone()?
								}
//(new showDoubleFloatArrays()).showArrays(kernels, "***kernels-"+nTX+"-"+nTY);
							} else {
								if (masterDebugLevel>1) System.out.println("Empty kernel for tile["+nTY+"]["+nTX+"]");
							}
							//save results into common array
							//pdfKernelMap[nTY][nTX]
						}
   						final int numFinished=tilesFinishedAtomic.getAndIncrement();
   						SwingUtilities.invokeLater(new Runnable() {
   							public void run() {
   								IJ.showProgress(numFinished,patternCells);
   							}
   						});

					}
				}
			};
		}
		startAndJoin(threads);
		if (globalDebugLevel>1) System.out.println("Threads done at "+IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
//		globalDebugLevel=saved_globalDebugLevel;
		if (debugLateralShifts) {
		    boolean [] debugUsedColors=new boolean [debugNumColors];
		    Arrays.fill(debugUsedColors,false);
			for (double [][][] dbgRow:debugLateral) for (double [][] dbgTile:dbgRow) if (dbgTile!=null) for (int c=0;c<dbgTile.length;c++){
				if (dbgTile[c]!=null) debugUsedColors[c]=true;
			}
			int numUsedColors=0;
			for (boolean b:debugUsedColors) if (b) numUsedColors++;
			if (numUsedColors>0) {
				int [] cIndex=new int [numUsedColors];
				int index=0;
				for (int i=0;i<debugUsedColors.length;i++) if (debugUsedColors[i]) cIndex[index++]=i;
				String [] dbgComponets={"latChromX","latChromY","shftX","shftY","centX","centY"};
				String [] debugTitles=new String[dbgComponets.length*numUsedColors];
				index=0;
				for (int i=0;i<dbgComponets.length;i++) for (int j=0;j<cIndex.length;j++){
					debugTitles[index++]=dbgComponets[i]+"-"+cIndex[j];
				}
				double [][] dbgLat= new double [debugTitles.length][debugLateral.length*debugLateral[0].length];
				int layer=0;
				for (int i=0;i<dbgComponets.length;i++) for (int j=0;j<cIndex.length;j++){
					Arrays.fill(dbgLat[layer],Double.NaN);
					for (nTileY=0;nTileY<debugLateral.length;nTileY++) for (nTileX=0;nTileX<debugLateral[0].length;nTileX++) if (debugLateral[nTileY][nTileX]!=null){
						if (debugLateral[nTileY][nTileX][cIndex[j]]!=null) dbgLat[layer][nTileY*debugLateral[0].length+nTileX]=debugLateral[nTileY][nTileX][cIndex[j]][i];
					} layer++;
				}				
				SDFA_INSTANCE.showArrays(dbgLat, debugLateral[0].length, debugLateral.length, true, "lateral"+imp_sel.getTitle(), debugTitles);
				/*
					debugLateralTile[i][0]=lateralChromatic[i][0];
					debugLateralTile[i][1]=lateralChromatic[i][1];
					debugLateralTile[i][2]=PSF_shifts[i][0];
					debugLateralTile[i][3]=PSF_shifts[i][1];
					debugLateralTile[i][4]=PSF_centroids[i][0];
					debugLateralTile[i][5]=PSF_centroids[i][1];
						    			sdfra_instance.showArrays(pointedBayer.clone(), halfWidth, halfHeight, true, title+"-bayer", subtitles);

				 */
			}
		}
		return pdfKernelMap;
	}
	/* Combine both greens as a checkerboard pattern (after oversampleFFTInput()) */
	private  double [][] combineCheckerGreens (double[][] input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
			int ratio) { // same as used in oversampleFFTInput() - oversampling ratio
		int width=(int) Math.sqrt(input_pixels[0].length);
		return combineCheckerGreens (input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
				width,   // width of the image
				ratio);
	}

	private  double [][] combineCheckerGreens (double[][] input_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
			int width,   // width of the image
			int ratio) { // same as used in oversampleFFTInput() - oversampling ratio
		if ((ratio<2) ||
				(input_pixels==null) ||
				((input_pixels.length>5) && (input_pixels[5]!=null)) ||
				(input_pixels.length<4) ||
				(input_pixels[0]==null) ||
				(input_pixels[3]==null)) return input_pixels;
		int height=input_pixels[0].length/width;
		int i,j;
		double [][] pixels={null,null,null,null,null,null};
		for (i=0;i<input_pixels.length;i++) pixels[i]=input_pixels[i];
		pixels[5]= new double[input_pixels[0].length];
		int index=0;
		int index_diff=(width+1)*ratio/2;
		double d;
		for (i=0;i<height;i++) for (j=0;j<width;j++) {
			d=input_pixels[0][index];
			if ((i>=ratio) && (j>=ratio)) d=0.5*(d+input_pixels[3][index-index_diff]);
			pixels[5][index++]=d;
		}
		return pixels;
	}
/*	
	private  double [] combineDiagonalGreens (double [] green0, double []green3, int half_width, int half_height) {
		int y,x,base;
		int base_b=0;
		double [] result= new double [green0.length];
		for (y=0;y<half_height/2; y++){
			base=half_height*half_width/2+ y* (half_width+1);
			for (x=0; x<half_width/2; x++) {
				result[base_b++]=green0[base];
				base-=half_width;
				result[base_b++]=green3[base++];
			}
			base=half_height*half_width/2+ y* (half_width+1);
			for (x=0; x<half_width/2; x++) {
				//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
				result[base_b++]=green3[base++];
				result[base_b++]=green0[base];
				base-=half_width;
			}
		}
		return result;
	}
	*/
	/* ======================================================================== */
	private  double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction) {
		return normalizeAndWindow (pixels, windowFunction, true);
	}
	/*
	private  double[] normalizeAndWindow (double [] pixels, double [] windowFunction) {
		return normalizeAndWindow (pixels, windowFunction, true);
	}
	*/
	private  double[][] normalizeAndWindow (double [][] pixels, double [] windowFunction, boolean removeDC) {
		int i;
		for (i=0;i<pixels.length;i++)  if (pixels[i]!=null) pixels[i]=normalizeAndWindow (pixels[i],  windowFunction, removeDC);
		return pixels;
	}
	private  double[] normalizeAndWindow (double [] pixels, double [] windowFunction, boolean removeDC) {
		int j;
		double s=0.0;
		if (pixels==null) return null;
		if (removeDC) {
			for (j=0;j<pixels.length;j++) s+=pixels[j];
			s/=pixels.length;
		}
		for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
		return pixels;
	}

	/* inserts zeros between pixels */ 
	private  double [][] oversampleFFTInput (double[][] input_pixels,
			int ratio) {
		double [][] pixels=new double[input_pixels.length][];
		int i;
		for (i=0;i<pixels.length;i++) pixels[i]= oversampleFFTInput (input_pixels[i], ratio);
		return pixels;
	}


	private  double [] oversampleFFTInput (double[] input_pixels, int ratio) {
		if (input_pixels==null) return null;
		int width=(int) Math.sqrt(input_pixels.length);
		return oversampleFFTInput (input_pixels,
				width,   // width of the image
				ratio);
	}

	private  double [] oversampleFFTInput (double[] input_pixels,
			int width,   // width of the image
			int ratio) {
		if (input_pixels==null) return null;
		double [] pixels=new double[input_pixels.length*ratio*ratio];
		int i,j,x,y;
		int height=input_pixels.length/width;
		for (i=0;i<pixels.length;i++) pixels[i]=0.0;
		j=0;
		for (y=0;y<height;y++) {
			i=width*ratio*ratio*y;
			for (x=0;x<width;x++) {
				pixels[i]=input_pixels[j++];
				i+=ratio;
			}
		}
		return pixels;
	}


/* ======================================================================== */
	private  void normalizeKernel(double [][] kernel) {
		int i;
		for (i=0;i<kernel.length;i++) if (kernel[i]!=null) normalizeKernel(kernel[i]);
	}

	private  void normalizeKernel(double [] kernel) {
		//	    if (kernel==null) return null;
		int i;
		double s=0;
		for (i=0;i<kernel.length;i++) s+= kernel[i];
		s=1.0/s;
		for (i=0;i<kernel.length;i++) kernel[i]*=s;
	}

/* ======================================================================== */
/* extends/shrinks image to make it square for FFT */
	private double[][] resizeForFFT (double[][]kernels, int size) {
		if (kernels==null) return null;
		double [][]result=new double [kernels.length][];
		for (int i=0;i<kernels.length;i++) {
			if (kernels[i]!=null) result[i]=resizeForFFT(kernels[i],size);
			else result[i]=null;
		}
		return result;
	}

	private double[] resizeForFFT (double[]kernel, int size) {
		int ksize=(int) Math.sqrt(kernel.length);
		double [] kernelForFFT = new double[size*size];
		int i,j,index, full_index;
//		if (DEBUG_LEVEL>10) System.out.println("resizeForFFT: new size="+size+" old size="+ksize);
		index=0;
		if (size==ksize) {
			return kernel.clone();
		} else if (size>ksize) {
			for (full_index=0;full_index<kernelForFFT.length; full_index++) kernelForFFT [full_index]=0.0;
			for (i=0;i<ksize; i++) {
				full_index=size* (size/2- ksize/2 + i) +size/2-ksize/2;
				for (j=0;j<ksize; j++) kernelForFFT[full_index++]=kernel[index++];
			}
		} else {
			for (i=0; i<size; i++) {
				full_index= ksize* (ksize/2-(size/2) +i) + (ksize/2-(size/2));
				for (j=0; j<size; j++) kernelForFFT[index++]=kernel[full_index++];
			}
		}
		return kernelForFFT;
	}
/* ======================================================================== */

	private  int kernelLength(double[][]kernels) {
		if (kernels==null) return 0;
		for (int i=0; i<kernels.length;i++) if (kernels[i]!=null) return kernels[i].length;
		return 0;
	}


	public double [][] getPSFKernels ( ImagePlus imp,
			float [][] simArray, //simulation image, scaled PSF_subpixel/2 (or null), [0] - main pixels, [1] - shifted diagonally by 0.5 pixels (for checker greens)
			int size,   // size in pixels (twice FFT_SIZE)
			int x0,      // top left corner X (pixels)
			int y0,      // top left corner Y (pixels)
			SimulationPattern     simulationPattern,
		    MatchSimulatedPattern matchSimulatedPattern,
			MatchSimulatedPattern.PatternDetectParameters patternDetectParameters,
			double [] Hamming, //=initHamming( fft_size) calculate once
			double [] fullHamming, //=initHamming( fft_size*subpixel);
			int subpixel, // use finer grid than actual pixels 
			SimulationPattern.SimulParameters  simulParameters,
			EyesisAberrations.ColorComponents colorComponents,
			EyesisAberrations.OTFFilterParameters otfFilterParameters,
			int referenceComp, // number of color component to reference lateral chromatic aberration to (now 4 - checkerboard greens)
			EyesisAberrations.PSFParameters psfParameters,
			DoubleFHT fht_instance, // provide DoubleFHT instance to save on initializations (or null)
			int masterDebugLevel, // get rid of it? // ** NEW
			int globalDebugLevel,// ** NEW
			int debug,
			double [][] debugLateralTile
	){
		boolean debugThis=false; //(y0==384) && ((x0==448) || (x0==512));
		if (globalDebugLevel>1){
			System.out.println("getPSFKernels(), simArray is "+((simArray==null)?"":"not ")+"null"); 
		}

		if (imp==null) return null; // Maybe convert to double pixel array once to make it faster?
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
		double [][] kernels=         new double[6][];  // was global
		String title=imp.getTitle()+"X"+x0+"Y"+y0;
		Rectangle PSFCell=new Rectangle (x0,y0,size,size);
		int fft_size=size/2;
		double [][] input_bayer=matchSimulatedPattern.splitBayer(imp,PSFCell,colorComponents.equalizeGreens); // does it work trhe same?
		//int greensToProcess=4;
		int i,j,l;
		double [][] simul_pixels;
		double [][]wVectors=new double[2][2];
		int imgWidth=imp.getWidth();
		
		double [][] dbgSimPix=null;
		
		double [] localBarray;
		
		if ((simArray==null) || (psfParameters.approximateGrid)){ // just for testing
			/* Calculate pattern parameters, including distortion */
			if (matchSimulatedPattern.PATTERN_GRID==null) {
				double[][] distortedPattern= matchSimulatedPattern.findPatternDistorted(input_bayer, // pixel array to process (no windowing!)
						patternDetectParameters,
						true, //(greensToProcess==4), // boolean greens, // this is a pattern for combined greens (diagonal), adjust results accordingly
						title); // title prefix to use for debug  images

				if (distortedPattern==null) return null;
				if (globalDebugLevel>3){
					System.out.println(
							" W0x="+     IJ.d2s(distortedPattern[0][0],4)+
							" W0y="+     IJ.d2s(distortedPattern[0][1],4)+
							" W0_phase="+IJ.d2s(distortedPattern[0][2],2)+
							" W1x="+     IJ.d2s(distortedPattern[1][0],4)+
							" W1y="+     IJ.d2s(distortedPattern[1][1],4)+
							" W1_phase="+IJ.d2s(distortedPattern[1][2],2));

				}
//				simulationPattern.simulatePatternFullPattern( // Not thread safe!
						localBarray=simulationPattern.simulatePatternFullPatternSafe(
						distortedPattern[0][0],
						distortedPattern[0][1],
						distortedPattern[0][2],
						distortedPattern[1][0],
						distortedPattern[1][1],
						distortedPattern[1][2],
						distortedPattern[2], //
						simulParameters.subdiv,
						fft_size,
						simulParameters.center_for_g2);
				wVectors[0][0]=2.0*distortedPattern[0][0]/subpixel;
				wVectors[0][1]=2.0*distortedPattern[0][1]/subpixel;
				wVectors[1][0]=2.0*distortedPattern[1][0]/subpixel;
				wVectors[1][1]=2.0*distortedPattern[1][1]/subpixel;
			} else { // approximate pattern grid and simulate it
				double[][] distPatPars= matchSimulatedPattern.findPatternFromGrid(
						x0, // top-left pixel of the square WOI
						y0,
						size, // size of square (pixels)
						Hamming, // only half-window!
						false,  // use linear approximation (instead of quadratic)
						1.0E-10,  // thershold ratio of matrix determinant to norm for linear approximation (det too low - fail)
						1.0E-20);  // thershold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
				int [] iUV={(int) Math.floor(distPatPars[0][2]),(int) Math.floor(distPatPars[1][2])};
				boolean negative=((iUV[0]^iUV[1])&1)!=0;
				double [] simCorr={
						distPatPars[0][3]/4,
						distPatPars[0][4]/4,
						distPatPars[0][5]/4,
						distPatPars[1][3]/4,
						distPatPars[1][4]/4,
						distPatPars[1][5]/4,
						0.0,0.0,0.0,0.0};
				double [] phases={
						1.0*Math.PI*(distPatPars[0][2]-iUV[0]+(negative?(-0.5):0.5)), // measured from the center of white
						1.0*Math.PI*(distPatPars[1][2]-iUV[1]+0.5)};
//				double [][]wVectors={{distPatPars[0][0],distPatPars[0][1]},{distPatPars[1][0],distPatPars[1][1]}};
				wVectors[0][0]=distPatPars[0][0];
				wVectors[0][1]=distPatPars[0][1];
				wVectors[1][0]=distPatPars[1][0];
				wVectors[1][1]=distPatPars[1][1];
//		    	simulationPattern.simulatePatternFullPattern( // Not thread safe!
				localBarray=simulationPattern.simulatePatternFullPatternSafe(
						wVectors[0][0],
						wVectors[0][1],
						phases[0],
						wVectors[1][0],
						wVectors[1][1],
						phases[1],
						simCorr, //
						simulParameters.subdiv,
						fft_size,
						simulParameters.center_for_g2);
			}
//			simul_pixels= simulationPattern.extractSimulPatterns (
			simul_pixels= simulationPattern.extractSimulPatterns (
					localBarray,		// this version is thread safe
					simulParameters,
					subpixel, // subdivide pixels
					fft_size*subpixel, // number of Bayer cells in width of the square selection (half number of pixels)
					0.0,    // selection center, X (in pixels)
					0.0);   // selection center, y (in pixels)
			if (subpixel>1) {
				if (colorComponents.colorsToCorrect[5])  simul_pixels=combineCheckerGreens (simul_pixels,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
						subpixel); // same as used in oversampleFFTInput() - oversampling ratio
			}
			for (i=0;i<simul_pixels.length; i++) {
				if (!colorComponents.colorsToCorrect[i]) simul_pixels[i]=null; // removed unused
			}
			simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
			
		} else {
			Rectangle PSFCellSim=new Rectangle (x0*subpixel/2,y0*subpixel/2,size*subpixel/2,size*subpixel/2);
			
			simul_pixels=new double[6][];
// simulationPattern.debugLevel=globalDebugLevel;
			for (i=0;i<simul_pixels.length; i++) {
				if (colorComponents.colorsToCorrect[i]) simul_pixels[i]=simulationPattern.extractBayerSim (
						simArray, // [0] - regular pixels, [1] - shifted by 1/2 diagonally, for checker greens
						imgWidth*subpixel/2,
						PSFCellSim,
						subpixel, // 4
						i);
				else simul_pixels[i]=null;
			}
//System.out.println("PSFCell.y="+PSFCell.y+" PSFCell.height="+PSFCell.height+" imgWidth="+imgWidth+" PSFCell.x="+PSFCell.x+" PSFCell.width="+PSFCell.width+" matchSimulatedPattern.UV_INDEX.length="+matchSimulatedPattern.UV_INDEX.length);
			int index=matchSimulatedPattern.getUVIndex((PSFCell.y+PSFCell.height/2)*imgWidth+(PSFCell.x+PSFCell.width/2));
//			int index=matchSimulatedPattern.getUVIndex((PSFCell.y+PSFCell.height/2)*matchSimulatedPattern.getWOI().width+(PSFCell.x+PSFCell.width/2));
			
			if (index<0) {
				System.out.println ("Error, No UV pattern @ x="+(PSFCell.x+PSFCell.width/2)+", y="+(PSFCell.y+PSFCell.height/2));
				return null;
			}
			//			int [] iUV={index % matchSimulatedPattern.getDArrayHeight(), index / matchSimulatedPattern.getDArrayHeight()}; // TODO: make sure it is correct?
			int [] iUV={index % matchSimulatedPattern.getDArrayWidth(), index / matchSimulatedPattern.getDArrayWidth()}; // TODO: make sure it is correct?
			if (matchSimulatedPattern.getDArray(iUV[1],iUV[0])==null) {
				if (globalDebugLevel>0){
					System.out.println ( "Tried to extract wave vectors from non-existent node "+iUV[0]+"/"+iUV[1]);
					System.out.println ( "index="+index+"  matchSimulatedPattern.getDArrayHeight()"+ matchSimulatedPattern.getDArrayHeight());
					System.out.println("PSFCell.y="+PSFCell.y+" PSFCell.height="+PSFCell.height+" imgWidth="+imgWidth+" PSFCell.x="+PSFCell.x+" PSFCell.width="+PSFCell.width+
							" matchSimulatedPattern.UV_INDEX.length="+matchSimulatedPattern.UV_INDEX.length);
				}
				return null;
			}
			if (matchSimulatedPattern.getDArray(iUV[1],iUV[0],1)==null) {
				if (globalDebugLevel>0) System.out.println ( "Tried to extract non-existent wave vectors from "+iUV[0]+"/"+iUV[1]);
				return null;
			}
			//TODO:  Need to define wave vectors here - how?
			wVectors[0]=matchSimulatedPattern.getDArray(iUV[1],iUV[0],1); //null pointer
			wVectors[1]=matchSimulatedPattern.getDArray(iUV[1],iUV[0],2);
			// should it be averaged WV?			
			if (globalDebugLevel>2) System.out.println ( " x0="+x0+" y0="+y0);
			if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(input_bayer, true, title+"-in");
			if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(simul_pixels, true, title+"-S");
			//if (globalDebugLevel>2) System.out.println (simArray[0][-1]); // cause error
			if (masterDebugLevel>1){
				dbgSimPix=new double[simul_pixels.length][];
				for (int ii=0;ii<dbgSimPix.length;ii++)
					if (simul_pixels[ii]!=null) dbgSimPix[ii]=simul_pixels[ii].clone();
					else dbgSimPix[ii]=null;

			}
			simul_pixels= normalizeAndWindow (simul_pixels, fullHamming);
		}
		
		input_bayer= normalizeAndWindow (input_bayer, Hamming);
		if (subpixel>1) {
			input_bayer= oversampleFFTInput (input_bayer,subpixel);
			if (colorComponents.colorsToCorrect[5])  input_bayer=combineCheckerGreens (input_bayer,   // pixel arrays after oversampleFFTInput() or extractSimulPatterns())
					subpixel); // same as used in oversampleFFTInput() - oversampling ratio
		}
		for (i=0;i<4;i++) if (!colorComponents.colorsToCorrect[i]) input_bayer[i]=null; // leave composite greens even if disabled
		if (debugThis) {
			SDFA_INSTANCE.showArrays(input_bayer, fft_size*subpixel, fft_size*subpixel, title);
		}
		if (globalDebugLevel>2) System.out.println ( " input_bayer.length="+input_bayer.length+" simul_pixels.length="+simul_pixels.length+" fft_size*subpixel="+fft_size*subpixel);
		for (i=0;(i<input_bayer.length) && (i<simul_pixels.length);i++) if ((colorComponents.colorsToCorrect[i]) && (input_bayer[i]!=null)){
			if (globalDebugLevel>2) System.out.println ( "input_bayer["+i+"].length="+input_bayer[i].length+" simul_pixels["+i+"].length="+simul_pixels[i].length);
		}
		if (debugThis) SDFA_INSTANCE.showArrays(input_bayer, true, title+"-input");
		if (debugThis) SDFA_INSTANCE.showArrays(simul_pixels, true, title+"-SIM");

if (globalDebugLevel>2)globalDebugLevel=0; //************************************************************
		double [][] inverted=new double[colorComponents.colorsToCorrect.length][];
		double wvAverage=Math.sqrt(0.5*(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1]+
				wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1]));

		for (i=0;(i<input_bayer.length) && (i<simul_pixels.length);i++) if ((colorComponents.colorsToCorrect[i]) && (input_bayer[i]!=null)){
			if (globalDebugLevel>2) System.out.println ( "Color "+colorComponents.getColorName(i)+" is re-calculated into bayer pixels ");
			if (globalDebugLevel>2) System.out.println ( "input_bayer["+i+"].length="+input_bayer[i].length+" simul_pixels["+i+"].length="+simul_pixels[i].length);
			inverted[i]=limitedInverseOfFHT(input_bayer[i],
					simul_pixels[i],
					fft_size*subpixel,
					(i==5),     //    boolean checker // checkerboard pattern in the source file (use when filtering)
					true, //      forwardOTF,
					subpixel,
					otfFilterParameters,
					fht_instance,
					psfParameters.mask1_sigma*size*wvAverage,      // normalize to wave vectors!
					psfParameters.mask1_threshold,
					psfParameters.gaps_sigma*size*wvAverage,
					psfParameters.mask_denoise,
					debug,
					globalDebugLevel,
					title+"-"+i);
		}
		int debugThreshold=1;
		if (debugThis) SDFA_INSTANCE.showArrays(inverted, fft_size*subpixel, fft_size*subpixel, title+"_Combined-PSF");
/* correct composite greens */
/* Here we divide wave vectors by subpixel as the pixels are already added */
		double [][] wVrotMatrix= {{0.5,0.5},{-0.5,0.5}};
		double [][]wVectors4= new double [2][2];
		for (i=0;i<2;i++) for (j=0;j<2;j++) {
			wVectors4[i][j]=0.0;
			for (l=0;l<2;l++) wVectors4[i][j]+=wVectors[i][l]*wVrotMatrix[l][j];
		}
		double [][] PSF_shifts=         new double [input_bayer.length][];    // X/Y shift of the PSF array, in Bayer component pixel coordinates (same as PSF arrays)
		double [][] PSF_centroids=      new double [input_bayer.length][];    // X/Y coordinates of the centroids of PSF in Bayer component pioxel coordinates (same as PSF arrays) (after they were optionally shifted)
		double [][] lateralChromatic=   new double [input_bayer.length][]; // X/Y coordinates of the centroids of Bayer component PSF in sensor pixel coordinates
		double [][] kernelsForFFT=      new double [input_bayer.length][];
		double [][] psf_inverted=       new double [input_bayer.length][];
		double [][] psf_inverted_masked=new double [input_bayer.length][];
		double [] lateralChromaticAbs=new double [input_bayer.length];
		double [] zeroVector={0.0,0.0};
		for (i=input_bayer.length-1;i>=0;i--) {
			if (colorComponents.colorsToCorrect[i]) {
				PSF_shifts[i]=       zeroVector.clone();
				PSF_centroids[i]=    zeroVector.clone();
				lateralChromatic[i]= zeroVector.clone();
			} else {
				PSF_shifts[i]=       null;
				PSF_centroids[i]=    null;
				lateralChromatic[i]= null;
			}
			lateralChromaticAbs[i]=0.0;
			kernelsForFFT[i]=null;
			psf_inverted[i]=null;
			psf_inverted_masked[i]=null;
		}
		//int [][]  clusterMask;
/* Start with referenceComp */
		i= referenceComp;
		if (globalDebugLevel>debugThreshold) {
			System.out.println(x0+":"+y0+"1-PSF_shifts.length= "+PSF_shifts.length+" i="+i+" input_bayer.length="+input_bayer.length);
			System.out.println("Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
					" PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
		}  


		kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
				!psfParameters.absoluteCenter, //true, // master, force ignoreChromatic
				PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if PSF_PARS.ignoreChromatic is true
				PSF_centroids[i], // will return array of XY coordinates of the result centroid
				(i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
						psfParameters,
						fht_instance,
						title+"_"+i,    // reduce the PSF cell size to this part of the area connecting first negative clones
						(globalDebugLevel>4),
						globalDebugLevel
						);
		if (globalDebugLevel>debugThreshold)     System.out.println(x0+":"+y0+"After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
		if (globalDebugLevel>debugThreshold)     System.out.println(x0+":"+y0+"After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));

		if (!psfParameters.ignoreChromatic && !psfParameters.absoluteCenter) { /* Recalculate center to pixels from greens (diagonal)) and supply it to other colors (lateral chromatic aberration correction) */
			for (j=0;j<input_bayer.length;j++) if ((colorComponents.colorsToCorrect[j]) && (j!=referenceComp)) {
				PSF_shifts[j]=shiftSensorToBayer (shiftBayerToSensor(PSF_shifts[referenceComp],referenceComp,subpixel),j,subpixel);
				if (globalDebugLevel>debugThreshold)       System.out.println(x0+":"+y0+"After-2 (recalc): color Component "+j+" PSF_shifts["+j+"][0]="+IJ.d2s(PSF_shifts[j][0],3)+" PSF_shifts["+j+"][1]="+IJ.d2s(PSF_shifts[j][1],3));
			}
		}

		lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
				PSF_shifts[i][1]+PSF_centroids[i][1],
				i,
				subpixel);
		lateralChromaticAbs[i]=Math.sqrt(lateralChromatic[i][0]*lateralChromatic[i][0]+lateralChromatic[i][1]*lateralChromatic[i][1]);
/* Now process all the other components */
		for (i=0; i<input_bayer.length;i++) if ((i!=referenceComp) && (colorComponents.colorsToCorrect[i])) {
			if (globalDebugLevel>debugThreshold) {
				System.out.println(x0+":"+y0+"2-PSF_shifts.length= "+PSF_shifts.length+" i="+i+" input_bayer.length="+input_bayer.length);

				System.out.println(x0+":"+y0+"Before: color Component "+i+" PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
						" PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3));
			}  
			kernels[i]=combinePSF (inverted[i], // Square array of pixels with multiple repeated PSF (alternating sign)
					false, // !master, use ignoreChromatic
					PSF_shifts[i],  // centerXY[] - will be modified inside combinePSF() if psfParameters.ignoreChromatic is true
					PSF_centroids[i], // will return array of XY coordinates of the result centroid
					(i==4)?wVectors4:wVectors, // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
							psfParameters,
							fht_instance,
							title+"_"+i,    // reduce the PSF cell size to this part of the area connecting first negative clones
							(globalDebugLevel>4),
							globalDebugLevel);
			if (globalDebugLevel>debugThreshold)     System.out.println(x0+":"+y0+"After-1: color Component "+i+"    PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts   [i][0],3)+"    PSF_shifts["+i+"][1]="+IJ.d2s(   PSF_shifts[i][1],3));
			if (globalDebugLevel>debugThreshold)     System.out.println(x0+":"+y0+"After-1: color Component "+i+" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
			lateralChromatic[i]=shiftBayerToSensor ( PSF_shifts[i][0]+PSF_centroids[i][0],
					PSF_shifts[i][1]+PSF_centroids[i][1],
					i,
					subpixel);
			lateralChromaticAbs[i]=Math.sqrt((lateralChromatic[i][0]-lateralChromatic[referenceComp][0])*(lateralChromatic[i][0]-lateralChromatic[referenceComp][0])+
					(lateralChromatic[i][1]-lateralChromatic[referenceComp][1])*(lateralChromatic[i][1]-lateralChromatic[referenceComp][1]));
		}
		if (globalDebugLevel>1) { //1
			for (i=0;i<PSF_shifts.length;i++) if (colorComponents.colorsToCorrect[i]){
				if (globalDebugLevel>debugThreshold) { //2
					System.out.println(x0+":"+y0+" Color Component "+i+" subpixel="+subpixel+
							" psfParameters.ignoreChromatic="+psfParameters.ignoreChromatic+
							" psfParameters.absoluteCenter="+psfParameters.absoluteCenter+
							" psfParameters.symm180="+psfParameters.symm180);
					System.out.println(x0+":"+y0+                     " PSF_shifts["+i+"][0]="+IJ.d2s(PSF_shifts[i][0],3)+
							" PSF_shifts["+i+"][1]="+IJ.d2s(PSF_shifts[i][1],3)+
							" PSF_centroids["+i+"][0]="+IJ.d2s(PSF_centroids[i][0],3)+
							" PSF_centroids["+i+"][1]="+IJ.d2s(PSF_centroids[i][1],3));
					System.out.println(x0+":"+y0+"  lateralChromatic["+i+"][0]="+IJ.d2s(lateralChromatic[i][0],3)+
							"  lateralChromatic["+i+"][1]="+IJ.d2s(lateralChromatic[i][1],3));
				}
			}
			if (colorComponents.colorsToCorrect[referenceComp]) for (i=0;i<colorComponents.colorsToCorrect.length;i++) if ((colorComponents.colorsToCorrect[i])&& (i!=referenceComp)){
				System.out.println("#!# "+x0+":"+y0+" "+colorComponents.getColorName(i)+" lateral chromatic (from green) "+IJ.d2s(lateralChromaticAbs[i],3)+"pix(sensor):  ["+i+"][0]="+IJ.d2s(lateralChromatic[i][0]-lateralChromatic[referenceComp][0],3)+
						"  ["+i+"][1]="+IJ.d2s(lateralChromatic[i][1]-lateralChromatic[referenceComp][1],3));
			}
			System.out.println("#!# "+x0+":"+y0+" "+"Lateral shift green from simulation "+IJ.d2s(lateralChromaticAbs[referenceComp],3)+"pix(sensor):  ["+referenceComp+"][0]="+IJ.d2s(lateralChromatic[referenceComp][0],3)+
					"  ["+referenceComp+"][1]="+IJ.d2s(lateralChromatic[referenceComp][1],3));
		}
		if (debugLateralTile!=null)	for (i=0;i<PSF_shifts.length;i++) {
				if (colorComponents.colorsToCorrect[i]){
					debugLateralTile[i]=new double [6];
					debugLateralTile[i][0]=lateralChromatic[i][0];
					debugLateralTile[i][1]=lateralChromatic[i][1];
					debugLateralTile[i][2]=PSF_shifts[i][0];
					debugLateralTile[i][3]=PSF_shifts[i][1];
					debugLateralTile[i][4]=PSF_centroids[i][0];
					debugLateralTile[i][5]=PSF_centroids[i][1];
			} else {
				debugLateralTile[i]=null;
			}
			
			
		}
		if (debugThis && (kernels!=null)){
			int debugSize=0;
			for (int ii=0;ii<kernels.length;ii++) if (kernels[ii]!=null){
				debugSize=(int)Math.sqrt(kernels[ii].length);
				break;
			}
			if (debugSize>0) SDFA_INSTANCE.showArrays(kernels, debugSize, debugSize, title+"_KERNELS");
		}
		return kernels;
	}
	/* ======================================================================== */
	/* shift (like lateral chromatic aberration) in Bayer component to sensor pixels */

		private  double [] shiftBayerToSensor ( double [] dxy,
				int color,
				int subPixel) {
			return shiftBayerToSensor (dxy[0], dxy[1], color, subPixel);
		}

		private  double [] shiftBayerToSensor ( double dx,
				double dy,
				int color,
				int subPixel) {
			double [] dxy=new double[2];
			switch (color) {
			case 5:
			case 0:
			case 1:
			case 2:
			case 3:dxy[0]=2.0*dx/subPixel;  dxy[1]= 2.0*dy/subPixel;  break;
			case 4:dxy[0]=(dx+dy)/subPixel; dxy[1]= (dy-dx)/subPixel; break;
			}
//			if (DEBUG_LEVEL>2)  System.out.println("shiftBayerToSensor(), color="+color+" subPixel="+subPixel+" ["+IJ.d2s(dx,3)+"/"+IJ.d2s(dy,3)+"] ->["+IJ.d2s(dxy[0],3)+"/"+IJ.d2s(dxy[1],3)+"]");
			return dxy;
		}

		private  double [] shiftSensorToBayer ( double [] dxy,
				int color,
				int subPixel) {
			return shiftSensorToBayer (dxy[0], dxy[1], color, subPixel);
		}
		private  double [] shiftSensorToBayer ( double dx,
				double dy,
				int color,
				int subPixel) {
			double [] dxy=new double[2];
			switch (color) {
			case 5:
			case 0:
			case 1:
			case 2:
			case 3:dxy[0]=0.5*dx*subPixel;      dxy[1]=0.5*dy*subPixel; break;
			case 4:dxy[0]=0.5*(dx-dy)*subPixel; dxy[1]=0.5*(dx+dy)*subPixel; break;
			}
//			if (DEBUG_LEVEL>2)  System.out.println("shiftSensorToBayer(), color="+color+" subPixel="+subPixel+" ["+IJ.d2s(dx,3)+"/"+IJ.d2s(dy,3)+"] ->["+IJ.d2s(dxy[0],3)+"/"+IJ.d2s(dxy[1],3)+"]");

			return dxy;
		}


	/* ======================================================================== */

	private double[] limitedInverseOfFHT(double [] measuredPixels,  // measured pixel array
			double [] modelPixels,  // simulated (model) pixel array)
			int size,  // FFT size
			boolean checker,  // checkerboard pattern in the source file (use when filtering)
			boolean forward_OTF,  // divide measured by simulated when true, simulated by measured - when false
			int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
			EyesisAberrations.OTFFilterParameters filterOTFParameters,  //  fraction of the maximal value to be used to limit zeros
			DoubleFHT fht_instance,  // add rejection of zero frequency (~2-3pix)
			double mask1_sigma,
			double mask1_threshold,
			double gaps_sigma,
			double mask_denoise,
			int debug,
			int globalDebugLevel,
			String title){ // title base for optional plots names
		return limitedInverseOfFHT(measuredPixels,
				modelPixels,
				size,
				checker,
				forward_OTF,
				oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
				filterOTFParameters.deconvInvert,
				filterOTFParameters.zerofreqSize,  // add rejection of zero frequency (~2-3pix)
				filterOTFParameters.smoothPS,       // 0 - none, otherwise Gauss width
				filterOTFParameters.thresholdHigh,  // reject completely if energy is above this part of maximal
				filterOTFParameters.thresholdLow,  // leave intact if energy is below this part of maximal
				-1.0, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
				0.0, // low-pass result with low pass filter (should be later defined automatically)
				fht_instance,
				mask1_sigma,
				mask1_threshold,
				gaps_sigma,
				mask_denoise,
				debug,
				globalDebugLevel,
				title);
	}
// TODO: It now selects a single PSF, so combinePSF() and binPSF() can be simplified and eliminated
	private double[] limitedInverseOfFHT(double [] measuredPixels,  // measured pixel array
			double [] modelPixels,  // simulated (model) pixel array)
			int size,  // FFT size
			boolean checker,  // checkerboard pattern in the source file (use when filtering)
			boolean forward_OTF,  // divide measured by simulated when true, simulated by measured - when false
			int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
			double deconvInvert,  //  fraction of the maximal value to be used to limit zeros
			double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
			double smoothPS,       // 0 - none, otherwise Gauss width = FFT size/2/smoothPS
			double threshold_high,  // reject completely if energy is above this part of maximal
			double threshold_low,  // leave intact if energy is below this part of maximal
			double threshold, // if 0 use normalize amplitude, if 0..1 - make binary: 1.0 if > threshold, 0.0 - otherwise -1 - disable mask
			double radius, // low-pass result with low pass filter (should be later defined automatically)
			DoubleFHT fht_instance, // provide DoubleFHT instance to save on initializations (or null)
			double mask1_sigma,
			double mask1_threshold,
			double gaps_sigma,
			double mask_denoise,
			int debug,
			int globalDebugLevel,
			String title){
        
		double [] denominatorPixels= forward_OTF? modelPixels.clone():    measuredPixels.clone();
		double [] nominatorPixels=   forward_OTF? measuredPixels.clone(): modelPixels.clone();
		if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
		int i;
		fht_instance.swapQuadrants(denominatorPixels);
		fht_instance.transform(denominatorPixels);
		double [] mask= null;
		double [] mask1=null;
		DoubleGaussianBlur gb=new DoubleGaussianBlur();
		if ((oversample>1) && (threshold_low<1.0)) {
			double [] ps=fht_instance.calculateAmplitude2(denominatorPixels);
/* create mask */
			mask= maskAliases (denominatorPixels,   // complex spectrum, [size/2+1][size]
					checker, // checkerboard pattern in the source file (use when filtering)
					oversample,   // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
					zerofreq_size,   // add rejection of zero frequency (~2-3pix)
					smoothPS,
					deconvInvert,
					threshold_high,   // reject completely if energy is above this part of maximal
					threshold_low,  // leave intact if energy is below this part of maximal
					fht_instance,
					globalDebugLevel);
/* debug show the mask */
			if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(mask, title+"-MASK");
			}
			for (int ii=0;ii<ps.length;ii++) ps[ii]=Math.log(ps[ii]); // can be twice faster
			if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(ps, "LOG-"+title);
			}
			double [] ps_smooth=ps.clone();
			gb.blurDouble(ps_smooth, size, size, mask1_sigma, mask1_sigma, 0.01);
			if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(ps_smooth, "SM-"+title);
			}
			double threshold1=Math.log(2.0*mask1_threshold);
			mask1=new double [ps.length];
			for (int ii=0;ii<ps.length;ii++) mask1[ii]= ps[ii]-ps_smooth[ii]-threshold1;
			if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(mask1, "M1-"+title);
			}
			fht_instance.swapQuadrants(mask1); // zero in the corner
			for (int ii=0;ii<mask1.length;ii++){
				if (mask1[ii]<0) {
					//				mask[ii]=0.0;
					mask1[ii]=0.0;
				}
				mask1[ii]*=mask[ii];
			}
			if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(mask1, "M1A-"+title);
			}
		}
/* Mask already includes zeros on ps, so we can just use divisions of FHT*/		  
		//Swapping quadrants of the nominator, so the center will be 0,0
		fht_instance.swapQuadrants(nominatorPixels);
		//get to frequency domain
		fht_instance.transform(nominatorPixels);
		if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug evel later */ // was 3
			SDFA_INSTANCE.showArrays(nominatorPixels, title+"-NOM-FHT");
			SDFA_INSTANCE.showArrays(denominatorPixels, title+"-DENOM-FHT");
		}			    
		double [] pixels=fht_instance.divide(nominatorPixels,denominatorPixels);
		if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug evel later */ // was 3
			SDFA_INSTANCE.showArrays(pixels, title+"-DECONV");
		}			    
		for (i=0;i<pixels.length;i++) {
			if (mask[i]==0.0) pixels[i]=0.0; // preventing NaN*0.0
			else pixels[i]*=mask[i];
		}
		if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
			SDFA_INSTANCE.showArrays(pixels, title+"-MASKED");
			double [][] aphase=fht_instance.fht2AmpHase(pixels,true);
			SDFA_INSTANCE.showArrays(aphase, true,"AP="+title+"-MASKED");
			
		}
		if (gaps_sigma>0.0){
			double [][] fft_reIm_centered=fht_instance.fht2ReIm(pixels, true); //0 in the center, full square
			fht_instance.swapQuadrants(mask1); // zero in the center
			for (int ii=0;ii<2;ii++) for (int jj=0;jj<mask1.length;jj++) fft_reIm_centered[ii][jj]*=mask1[jj];
			gb.blurDouble(mask1, size, size, gaps_sigma, gaps_sigma, 0.01);
			gb.blurDouble(fft_reIm_centered[0], size, size, gaps_sigma, gaps_sigma, 0.01);
			gb.blurDouble(fft_reIm_centered[1], size, size, gaps_sigma, gaps_sigma, 0.01);
			for (int ii=0;ii<2;ii++) for (int jj=0;jj<mask1.length;jj++)
				if (mask1[jj]>mask_denoise) fft_reIm_centered[ii][jj]/=mask1[jj];
				else if (mask1[jj]>=0.0) fft_reIm_centered[ii][jj]/=mask_denoise;
				else  fft_reIm_centered[ii][jj]=0.0;
			if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
				SDFA_INSTANCE.showArrays(fft_reIm_centered, true,"ReIm-"+title);
			}
			fht_instance.swapQuadrants(fft_reIm_centered[0]); // zero in the corner
			fht_instance.swapQuadrants(fft_reIm_centered[1]); // zero in the corner
			pixels=fht_instance.FFTHalf2FHT(fft_reIm_centered, size);
		//mask_denoise	
		}
		/// transform to space
		fht_instance.inverseTransform(pixels);
		fht_instance.swapQuadrants(pixels);
		if ((debug>2) ||((globalDebugLevel>2) && (title!=""))) { /* Increase debug level later */ // was 3
			SDFA_INSTANCE.showArrays(pixels, "PSF-"+title);
		}
		return pixels;
	}

	/* ======================================================================== */
	/* Trying to remove aliasing artifacts when the decimated (pixel resolution) image is deconvolved with full resolution (sub-pixel resolution)
	model pattern. This effect is also easily visible if the decimated model is deconvolved with the same one art full resolution.
	Solution is to clone the power spectrum of the full resolution model with the shifts to match oversampling (15 clones for the 4x oversampling),
	And add them together (adding also zero frequerncy point - it might be absent on the model) but not include the original (true one) and
	use the result to create a rejectiobn mask - if the energy was high, (multiplicative) mask should be zero at those points. */

		private double [] maskAliases (double [] fht, // complex spectrum, [size/2+1][size]
				boolean checker, // checkerboard pattern in the source file (use when filtering)
				int oversample,  // measured array is sampled at 1/oversample frequency than model (will add more parameters later)
				double zerofreq_size,  // add rejection of zero frequency (~2-3pix)
				double sigma,
				double deconvInvert,
				double threshold_high,  // reject completely if energy is above this part of maximal
				double threshold_low,  // leave intact if energy is below this part of maximal
				DoubleFHT fht_instance,
				int globalDebugLevel){ // provide DoubleFHT instance to save on initializations (or null)
			double th=threshold_high*threshold_high;
			double tl=threshold_low*threshold_low;

			int length=fht.length;
			int size=(int) Math.sqrt(fht.length);
			//	double [][] ps=new double [size/2+1][size];
			int i,ix,iy, cloneNx, cloneNy, cloneX, cloneY;
			int cloneStep=size/oversample;
	/* generating power spectrum for the high-res complex spectrum, find maximum value and normalize */
			if (fht_instance==null) fht_instance=new DoubleFHT(); // move upstream to reduce number of initializations
			double [] ps=fht_instance.calculateAmplitude2(fht);
			double psMax=0.0;
			for (i=0;i<length; i++) if (psMax<ps[i]) psMax=ps[i];
			double k=1.0/psMax;
			for (i=0;i<length; i++) ps[i]*=k;
			if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(ps, "PS");
	/* Add maximum at (0,0) */
			double [] psWithZero=ps;
			if (zerofreq_size>0.0) {
				psWithZero=ps.clone();
				int zs=(int) (4*zerofreq_size);
				int base=size*(size+1)/2;
				k=0.5/(zerofreq_size*zerofreq_size);
				if (zs>=size/2) zs =size/2;
				for (iy=-zs;iy<=zs;iy++) for (ix=-zs; ix <= zs; ix++) {
					psWithZero[base+iy*size+ix]+=Math.exp(-k*(iy*iy+ix*ix));
				}
			}
	/* put zero in the center */
			double [] mask=new double [length];
			for (i=0;i<length; i++) mask[i]=0.0;
	/* clone spectrums */
			for (iy=0;iy<size;iy++) for (ix=0;ix<size;ix++){
				for (cloneNy=0;cloneNy<oversample;cloneNy++) for (cloneNx=0;cloneNx<oversample;cloneNx++)
					if (((cloneNy!=0) || (cloneNx!=0)) && // not a zero point
							(!checker ||                      // use all if it is not a checkerboard pattren
									(((cloneNx ^ cloneNy) & 1)==0) )) { // remove clones in a checker pattern
						cloneY=(iy+cloneNy*cloneStep)%size;
						cloneX=(ix+cloneNx*cloneStep)%size;
						mask[cloneY*size+cloneX]+=psWithZero[iy*size+ix];
					}			
			}
	/* debug show the mask */
			if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(mask, "PS-cloned");
			if (sigma>0) {
				DoubleGaussianBlur gb = new DoubleGaussianBlur();
				gb.blurDouble(mask,size,size,sigma,sigma, 0.01);
				if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(mask, "PS-smooth");
			}

	/* make mask of cloned power spectrums */
			double a;
			double k2=deconvInvert*deconvInvert;
			double min=0.01*k2; // less than 1/10 of that value - mask=0.0
			if (globalDebugLevel>2) System.out.println("maskAliases() threshold_high="+threshold_high+" threshold_low="+threshold_low+" th="+th+" tl="+tl+" k2="+k2+" min="+min);
			for (i=0;i<length;i++) {
				if      (mask[i]<tl)  mask[i]=1.0;
				else if (mask[i]>th) mask[i]=0.0;
				else { // make smooth transition
					a=(2.0 * mask[i] - th - tl)/(th - tl);
					mask[i]=0.5*(1.0-a*a*a);
				}
				// now mask out zeros on the ps		
				if (ps[i]<min) mask[i]=0.0;
				else {
					mask[i]*=ps[i]/(ps[i]+k2);
				}
			}
			if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(mask, "mask-all");
			/* zeros are now for FHT - in the top left corner */	
			fht_instance.swapQuadrants(mask);
			return mask;
		}


	private boolean [][] mapFromPatternMask (
			MatchSimulatedPattern matchSimulatedPattern, // to use windowFunction
			int width, // image (mask) width
			int tileSize,
			int tileStep,
			int margin,   // backward compatibility margin==tileSize/2
			double gaussWidth,
			double threshold,
			int debugLevel){
		int [] uvIndex=matchSimulatedPattern.getUVIndex(); // int array, >=0 - uv exist, <0 - empty
		if (uvIndex==null) return null;
		double[] windowFunction= matchSimulatedPattern.initWindowFunction(tileSize, gaussWidth);
		int height =uvIndex.length/width;
		int tileHeight=(height-2*margin)/tileStep+1;
		int tileWidth= (width- 2*margin)/tileStep+1;
		boolean [][] result = new boolean [tileHeight][tileWidth];
		int index;
		int len=tileSize*tileSize;
		double absThresh=0.0, sum;
		for (index=0;index<len;index++) absThresh+=windowFunction[index];
		absThresh*=threshold;
		if (debugLevel>1) System.out.println(" threshold="+threshold+" absThresh="+absThresh);
		
		int y,x,y0,x0;
		for (int tileY=0;tileY<tileHeight;tileY++) for (int tileX=0;tileX<tileWidth;tileX++) {
			y0=-tileSize/2+margin+tileStep*tileY;
			x0=-tileSize/2+margin+tileStep*tileX;
			sum=0;
			for (index=0;index<len;index++) {
				y=index/tileSize+y0;
				x=index%tileSize+x0;
				if ((y>=0) && (x>=0) && (y<height) && (x<width) && (uvIndex[y*width+x]>=0)) sum+= windowFunction[index];
//				if ((globalDebugLevel>0) && (tileY==22) && (tileX==32)) System.out.println(" x="+x+" y="+y);
//				if ((globalDebugLevel>0) && (tileY==22) && (tileX==32) && (y>=0) && (x>=0) && (y<height) && (x<width))System.out.println(" uvIndex["+(y*width+x)+"]="+uvIndex[y*width+x]);
			}
			result[tileY][tileX]=(sum>absThresh);
			if (debugLevel>1) System.out.println(" tileY="+tileY+" tileX="+tileX+" x0="+x0+" y0="+y0+" sum="+sum+" threshold="+threshold+" absThresh="+absThresh+" rrsult="+result[tileY][tileX]);
		}
		return result;
	}
	
	/* ======================================================================== 
	/**
	 * Mostly done, need to move where szis\
	 * TODO: currently the shift of the PSF during binning is done with the integer steps. If ignoreChromatic - to all colors
	 * independently, if it is false - all components are moved in sync, but again - with integer steps. That causes
	 * mis-match between the PSF calculated in nearly identical runs (i.e. use the data shifted by 2 pixels) caused by 1 pixel shift.
	 * That can be improved if PSF are shifted smoothly (not so easy though). It is probably already handled when averaging PSF - 
	 * amplitude and phase is handled separately so shift should be OK.
	 * 	
	 */
		
		
		double [] combinePSF (double []pixels,         // Square array of pixels with multiple repeated PSF (alternating sign)
				boolean   master,          // force ignoreChromatic
				double[] centerXY,         // coordinates (x,y) of the center point (will update if ignoreChromatic is true)
				double [] centroid_xy,    // RETURNS centroid of the result array (should be small) if ignoreChromatic is true
				double [][] wVectors,    // two wave vectors, lengths in cycles/pixel (pixels match pixel array)
				EyesisAberrations.PSFParameters psfParameters,    // minimal instance contrast to use in binning
				DoubleFHT fht_instance, // provide DoubleFHT instance to save on initializations (or null) // used for sub-pixel shift, null OK
				String title,     // reduce the PSF cell size to this part of the area connecting first negative clones
				boolean debug,
				int debugLevel)
		{
			if (pixels==null) return null;
			//    double [] contrastCache=new double[pixelSize*pixelSize];
			int i,j;

			if (debugLevel>2) {
				System.out.println("combinePSF title="+title+" wV[0][0]="+IJ.d2s(wVectors[0][0],4)+" wV[0][1]="+IJ.d2s(wVectors[0][1],4));
				System.out.println("combinePSF title="+title+" wV[1][0]="+IJ.d2s(wVectors[1][0],4)+" wV[1][1]="+IJ.d2s(wVectors[1][1],4));
			}

	/* vectors perpendicular to the checkerboard edges, lengths equal to the periods */
			double [][] f= {{wVectors[0][0]/(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1]),
				wVectors[0][1]/(wVectors[0][0]*wVectors[0][0]+wVectors[0][1]*wVectors[0][1])},
				{wVectors[1][0]/(wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1]),
					wVectors[1][1]/(wVectors[1][0]*wVectors[1][0]+wVectors[1][1]*wVectors[1][1])}};
			if (debugLevel>2) {
				System.out.println("combinePSF title="+title+" f[0][0]="+IJ.d2s(f[0][0],4)+" f[0][1]="+IJ.d2s(f[0][1],4));
				System.out.println("combinePSF title="+title+" f[1][0]="+IJ.d2s(f[1][0],4)+" f[1][1]="+IJ.d2s(f[1][1],4));
			}

	/* vectors parallel to checkerboard edges, lenghs equal to the period along those lines */
			double l2f1=   f[0][0]*f[0][0]+f[0][1]*f[0][1];
			double l2f2=   f[1][0]*f[1][0]+f[1][1]*f[1][1];
			double pf1f2  =f[0][1]*f[1][0]-f[1][1]*f[0][0];
			double [][]g0= {{f[0][1]*l2f2/pf1f2,  -f[0][0]*l2f2/pf1f2},
					{f[1][1]*l2f1/pf1f2,  -f[1][0]*l2f1/pf1f2}};
			if (debugLevel>2) {
				System.out.println("combinePSF title="+title+" g0[0][0]="+IJ.d2s(g0[0][0],4)+" g[0][1]="+IJ.d2s(g0[0][1],4));
				System.out.println("combinePSF title="+title+" g0[1][0]="+IJ.d2s(g0[1][0],4)+" g[1][1]="+IJ.d2s(g0[1][1],4));
			}
	/* calculate vectors connecting centers of the "positive" PSF copies */

			double [][] g= {{0.5*(g0[0][0]+g0[1][0]), 0.5*(g0[0][1]+g0[1][1])},
					{0.5*(g0[0][0]-g0[1][0]), 0.5*(g0[0][1]-g0[1][1])}};

			if (debugLevel>2) {
				System.out.println("combinePSF title="+title+" g[0][0]="+IJ.d2s(g[0][0],4)+" g[0][1]="+IJ.d2s(g[0][1],4));
				System.out.println("combinePSF title="+title+" g[1][0]="+IJ.d2s(g[1][0],4)+" g[1][1]="+IJ.d2s(g[1][1],4));
			}
			/// =================

	/* calculate outSize to be able to use FFT here */
			double sizeNegatives= Math.max(Math.max(Math.abs(g[0][0]+ g[1][0]),Math.abs(g[0][1]+ g[1][1])),
					Math.max(Math.abs(g[0][0]- g[1][0]),Math.abs(g[0][1]- g[1][1])));
			double scaleSize=2.5; /// Will include next positive centers and overlap
			int outSize;
			for (outSize=8;outSize<scaleSize*sizeNegatives; outSize<<=1);
			int halfOutSize=outSize/2;
			if (debugLevel>2) {
				System.out.println("sizeNegatives="+sizeNegatives+ " scaled="+ (scaleSize*sizeNegatives)+" outSize="+outSize+" halfOutSize="+halfOutSize);
			}

			double [] pixelsPSF= binPSF(pixels,
					g,
					outSize,
					psfParameters.minContrast,  // minimal contrast of PSF clones
					centerXY,  //  coordinates (x,y) of the center point
					null,  // coordinates of the center of symmetry - not applicable
					1, // pass 1
					title,
					debug,
					debugLevel);
			//                   true);
			
			if (!master && !psfParameters.ignoreChromatic && !psfParameters.absoluteCenter && psfParameters.centerPSF && (centerXY!=null)){
//				System.out.println("1:pixelsPSF.length="+pixelsPSF.length+" outSize+"+outSize);

				// TODO: Shift +/- 0.5 Pix here {centerXY[0]-Math.round(centerXY[0]),centerXY[1]-Math.round(centerXY[1])}	
				if (fht_instance==null) fht_instance=new DoubleFHT();
//				fht_instance.debug=(centerXY[0]-Math.round(centerXY[0]))<-0.4; // just reducing number
//				double dx=centerXY[0]-Math.round(centerXY[0]);
//				double dy=centerXY[1]-Math.round(centerXY[1]);
//				if (dx<-0.4) SDFA_INSTANCE.showArrays(pixelsPSF.clone(), "before:"+dx+":"+dy);
	        
				pixelsPSF=fht_instance.translateSubPixel (
						 pixelsPSF,
						 -(centerXY[0]-Math.round(centerXY[0])),
						 -(centerXY[1]-Math.round(centerXY[1])));
//				fht_instance.debug=false;
//				if (dx<-0.4) SDFA_INSTANCE.showArrays(pixelsPSF.clone(), "after:"+dx+":"+dy);

			}

			double distToNegativeClones=0.5*Math.sqrt(((g[0][0]+g[1][0])*(g[0][0]+g[1][0])+
					(g[0][1]+g[1][1])*(g[0][1]+g[1][1])+
					(g[0][0]-g[1][0])*(g[0][0]-g[1][0])+
					(g[0][1]-g[1][1])*(g[0][1]-g[1][1]))/2.0);
			if (debugLevel>2) {
				System.out.println("distToNegativeClones="+distToNegativeClones+ " gaussWidth="+ distToNegativeClones*psfParameters.smoothSeparate);
			}
			double smoothSigma=distToNegativeClones*psfParameters.smoothSeparate;

			//	double [] smoothPixelsPSF= lowPassGauss(pixelsPSF, smoothSigma, true);
			double [] smoothPixelsPSF= pixelsPSF.clone();
			DoubleGaussianBlur gb=new DoubleGaussianBlur();
			gb.blurDouble(smoothPixelsPSF, outSize, outSize, smoothSigma, smoothSigma, 0.01);

	/* find amplitude of smoothed pixel array */
			double smoothMin=0.0;
			double smoothMax=0.0;
			for (i=0;i<smoothPixelsPSF.length;i++) {
				if      (smoothPixelsPSF[i] > smoothMax) smoothMax=smoothPixelsPSF[i];
				else if (smoothPixelsPSF[i] < smoothMin) smoothMin=smoothPixelsPSF[i];
			}
			int [][]  clusterMask = findClusterOnPSF(smoothPixelsPSF, // PSF function, square array (use smooth array)
					-psfParameters.topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
					outSize/2,  // location of a start point, x-coordinate
					outSize/2,  // location of a start point, y-coordinate
					title,
					debugLevel);
			double [] centroidXY=       calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
					//--centroidXY is in function call arguments
					//centroidXY=            calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
					clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
					psfParameters.topCenter);// subtract level below topCenter*max
			double [] centroidXY_smooth=calcCentroidFromCenter(smoothPixelsPSF, // use smooth - not final, just for clones rejection
					clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
					psfParameters.topCenter);// subtract level below topCenter*max

			if (debugLevel>2) System.out.println("Centroid after first binPSF: x="+IJ.d2s(centroidXY[0],3)+" y="+IJ.d2s(centroidXY[1],3)+" center was at x="+IJ.d2s(centerXY[0],3)+" y="+IJ.d2s(centerXY[1],3));

	/* Re-bin results with the new center if ignoreChromatic is true, update centerXY[](shift of the result PSF array) and centroidXY[] (center of the optionally shifted PDF array) */
			if (!psfParameters.absoluteCenter && (master || psfParameters.ignoreChromatic)) {
				if (centerXY!=null) {
					centerXY[0]+=centroidXY[0];
					centerXY[1]+=centroidXY[1];
				}
				pixelsPSF= binPSF(   pixels,
						g,
						outSize,
						psfParameters.minContrast,  // minimal contrast of PSF clones
						centerXY,  // now includes centroid from the pass 1
						psfParameters.symm180?centroidXY:null,
								2, // pass2
								title,
								debug,
								debugLevel);
				if (psfParameters.centerPSF && (centerXY!=null)){
//					System.out.println("2:pixelsPSF.length="+pixelsPSF.length+" outSize+"+outSize);
					// TODO: Shift +/- 0.5 Pix here {centerXY[0]-Math.round(centerXY[0]),centerXY[1]-Math.round(centerXY[1])}	
					if (fht_instance==null) fht_instance=new DoubleFHT();
//					fht_instance.debug=(centerXY[0]-Math.round(centerXY[0]))<-0.4; // just reducing number
					pixelsPSF=fht_instance.translateSubPixel (
							 pixelsPSF,
							 -(centerXY[0]-Math.round(centerXY[0])),
							 -(centerXY[1]-Math.round(centerXY[1])));
//					fht_instance.debug=false;
				}
	/*  recalculate centroids  */
				smoothPixelsPSF= pixelsPSF.clone();
				gb.blurDouble(smoothPixelsPSF, outSize, outSize, smoothSigma, smoothSigma, 0.01);
				smoothMin=0.0;
				smoothMax=0.0;
				for (i=0;i<smoothPixelsPSF.length;i++) {
					if      (smoothPixelsPSF[i] > smoothMax) smoothMax=smoothPixelsPSF[i];
					else if (smoothPixelsPSF[i] < smoothMin) smoothMin=smoothPixelsPSF[i];
				}
				clusterMask = findClusterOnPSF(smoothPixelsPSF, // PSF function, square array (use smooth array)
						-psfParameters.topCenter, // fraction of energy in the pixels to be used (or minimal level if it is negative)
						outSize/2,  // location of a start point, x-coordinate
						outSize/2,  // location of a start point, y-coordinate
						title,
						debugLevel);
				centroidXY= calcCentroidFromCenter(pixelsPSF, // use original array (mask from the smoothed one)
						clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
						psfParameters.topCenter);// subtract level below topCenter*max
				// seems it is not used anymore
				centroidXY_smooth=calcCentroidFromCenter(smoothPixelsPSF, // use smooth - not final, just for clones rejection
						clusterMask, // integer mask -0 - don't use this pixel, 1 - use it
						psfParameters.topCenter);// subtract level below topCenter*max
				if (debugLevel>2) System.out.println("Centroid after second binPSF: x="+IJ.d2s(centroidXY[0],3)+" y="+IJ.d2s(centroidXY[1],3)+" center was at x="+IJ.d2s(centerXY[0],3)+" y="+IJ.d2s(centerXY[1],3));

			}
			
			
			
	/* compensate center point and/or add center-symmetrical points if enabled */
			double [] rejectedClonesPixels=null;
			double [][] modelPSFVectors={{0.5*(g[0][0]+g[1][0]),0.5*(g[0][1]+g[1][1])},
					{0.5*(g[0][0]-g[1][0]),0.5*(g[0][1]-g[1][1])}};
	/********* removed subtraction of clones *****************************************************************/		
			rejectedClonesPixels=pixelsPSF; // Maybe fo the opposite?
			maskClonesPSF(rejectedClonesPixels, // square pixel array where the model PSF is added
					psfParameters.windowFrac, // multiply window by this value
					centroidXY[0], // Center of the remaining single PSF
					centroidXY[1], // same for Y
					modelPSFVectors, // vectors that connect center of PSF with two oppositre sign clones
					psfParameters.useWindow);  // use Hamming window, if false - just cut sharp

			if (psfParameters.wingsEnergy>0.0) {
				rejectedClonesPixels=cutPSFWings (rejectedClonesPixels, // direct PSF function, square array, may be proportionally larger than reversed
						psfParameters.wingsEnergy, // fraction of energy in the pixels to be used
						psfParameters.wingsEllipseScale,
						0.003, // wings_min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
						title+"-w",
						debugLevel);
			}
			double [] sigmas=createSigmasRadius(rejectedClonesPixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
					psfParameters.sigmaToRadius, // sigma is proportional to the distance from the center
					centroidXY[0], // model PSF center X-coordinate (in pixels[] units, from the center of the array )
					centroidXY[1], // same for Y
					0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
					0, // int WOICenterY, // 
					outSize, //int WOIWidth, reduce later
					outSize); //int WOIHeight)

			double max1=0;
			for (i=0;i<smoothPixelsPSF.length;i++) if (smoothPixelsPSF[i]>max1) max1=smoothPixelsPSF[i];
			double minSigma=0.5;
			double varSigmaTop=1.0 ; //0.7;
			double kk;

			for (i=0;i<sigmas.length;i++) {
				kk=smoothPixelsPSF[i]/max1;
				if (kk>varSigmaTop) sigmas[i]=minSigma;  
				else                sigmas[i] = minSigma+ sigmas[i]*((varSigmaTop-kk)*(varSigmaTop-kk)/varSigmaTop/varSigmaTop);
			}
			double [] varFilteredPSF=variableGaussBlurr(rejectedClonesPixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
					sigmas, // array of sigmas to be used for each pixel, matches pixels[]
					3.5, // drop calculatin if farther then nSigma
					0, // int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
					0, // int WOICenterY, // 
					outSize, //int WOIWidth, reduce later
					outSize,
					debugLevel); //int WOIHeight)


			if (debugLevel>2) {
	/* Sigmas are 0 here ??? */
				if (psfParameters.sigmaToRadius>0.0) {
					float [] floatPixelsSigmas=new float[sigmas.length];
					for (j=0;j<sigmas.length;j++) floatPixelsSigmas[j]=(float) sigmas[j];
					ImageProcessor ip_Sigmas=new FloatProcessor(outSize,outSize);
					ip_Sigmas.setPixels(floatPixelsSigmas);
					ip_Sigmas.resetMinAndMax();
					ImagePlus imp_Sigmas=  new ImagePlus(title+"_Sigmas", ip_Sigmas);
					imp_Sigmas.show();
				}

				System.out.println("title="+title+" center X(pix)="+centroidXY_smooth[0]+"(smooth) center Y(pix)="+centroidXY_smooth[1]+"(smooth)");
				System.out.println("title="+title+" center X(pix)="+centroidXY[0]+"          center Y(pix)="+centroidXY[1]);
			}
			centroid_xy[0]=centroidXY[0];
			centroid_xy[1]=centroidXY[1];
			return  varFilteredPSF;
		}
		/* ======================================================================== */
		public double [][] matrix2x2_invert(double [][] m ){
			double det=m[0][0]*m[1][1]-m[0][1]*m[1][0];
			double [][] rslt= {{ m[1][1]/det,  -m[0][1]/det},
					{-m[1][0]/det,   m[0][0]/det}};
			return rslt;
		}
		public double [][] matrix2x2_mul(double [][] a, double [][] b ){
			double [][] rslt={{a[0][0]*b[0][0]+a[0][1]*b[1][0], a[0][0]*b[0][1]+a[0][1]*b[1][1]},
					{a[1][0]*b[0][0]+a[1][1]*b[1][0], a[1][0]*b[0][1]+a[1][1]*b[1][1]}};
			return rslt;
		}
		public double [] matrix2x2_mul(double [][] a, double [] b ){
			double [] rslt={a[0][0]*b[0]+a[0][1]*b[1],
					a[1][0]*b[0]+a[1][1]*b[1]};
			return rslt;
		}
		public double [][] matrix2x2_scale(double [][] a, double  b ){
			double [][] rslt={{a[0][0]*b, a[0][1]*b},
					{a[1][0]*b, a[1][1]*b}};
			return rslt;
		}

		public double [][] matrix2x2_add(double [][] a, double [][] b ){
			double [][] rslt={{a[0][0]+b[0][0], a[0][1]+b[0][1]},
			         		  {a[1][0]+b[1][0], a[1][1]+b[1][1]}};
			return rslt;
		}

		public double [] matrix2x2_add(double [] a, double [] b ){
			double [] rslt={a[0]+b[0], a[1]+b[1]};
			return rslt;
		}
		
		public double [][] matrix2x2_transp(double [][] m ){
			double [][] rslt= {{ m[0][0],  m[1][0]},
			            	   { m[0][1],  m[1][1]}};
			return rslt;
		}
		
		/* ======================================================================== */
		/* zeroes out area outside of the area bound by 4 negative clones (or a fraction of it), either sharp or with Hamming */
			private double [] maskClonesPSF(double [] pixels, // square pixel array where the model PSF is added
					double windowPart, // multiply window by this value
					double xc, // Center of the remaining single PSF
					double yc, // same for Y
					double[][] vectors, // vectors that connect center of PSF with two oppositre sign clones
					boolean  useHamming  // use Hamming window, if false - just cut sharp
			) {
				int ix,iy;
				int size = (int) Math.sqrt (pixels.length);
				double [] xy= new double[2];
				double [] uv;
		/* matrix that converts u,v (lengths along the) 2 input vectors connecting opposite sign PSFs into x,y coordinates */
				double [][] uv2xy= {{vectors[0][0]*windowPart,vectors[1][0]*windowPart},
						{vectors[0][1]*windowPart,vectors[1][1]*windowPart}};
				double [][] xy2uv=  matrix2x2_invert(uv2xy);
				for (iy=0;iy<size;iy++) {
					xy[1]=(iy-size/2)-yc;
					for (ix=0;ix<size;ix++) {
						xy[0]=(ix-size/2)-xc;
						uv=matrix2x2_mul(xy2uv, xy);
						if ((Math.abs(uv[0])>1.0) || (Math.abs(uv[1])>1.0)) pixels[iy*size+ix]=0.0;
						else if (useHamming) {
							pixels[iy*size+ix]*=(0.54+0.46*Math.cos(uv[0]*Math.PI))*(0.54+0.46*Math.cos(uv[1]*Math.PI));
						}
					}
				}
				return pixels;
			}
			/* ======================================================================== */
			private double [] variableGaussBlurr (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
					double []sigmas, // array of sigmas to be used for each pixel, matches pixels[]
					double nSigma, // drop calculatin if farther then nSigma
					int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
					int WOICenterY, // 
					int WOIWidth, //
					int WOIHeight,
					int globalDebugLevel){ //
				int size = (int) Math.sqrt(pixels.length);
				double [] result =new double [size*size];
				double [] gauss= new double [2*size];
				int x0= (size-WOIWidth)/2 +WOICenterX;
				int y0= (size-WOIHeight)/2+WOICenterY;
				int x1=x0+WOIWidth;
				int y1=x0+WOIHeight;
				int i,ix,iy,max_i;
				double sum,k,sigma,d,gy,scale,g;
				int xk0,xk1,yk0,yk1, ikx,iky, index;
				for (i=0;i<result.length;i++) result[i]=0.0;
				if (globalDebugLevel>2) {
					System.out.println(" variableGaussBlurr(), x0="+x0+" y0="+y0+" x1="+x1+" y1="+y1);
				}
				if (x0<0) x0=0; if (x1>size) x1=size; if (y0<0) y0=0; if (y1>size) y1=size;
				for (iy=0;iy<size;iy++) {
					for (ix=0;ix<size;ix++) {
						d=pixels[iy*size+ix];
						if (d!=0.0) {
							sigma=sigmas[iy*size+ix];
							if (sigma==0.0) {
								result[iy*size+ix]+=d; // just copy input data, no convolving
							} else {
		/* opposite to "normal" convolution we have diffrent kernel for each point, so we need to make sure that two points with the same values but
		  diffrent sigma values will not move "energy" from one to another. For this we can do accumulation both ways - from the source point to all
		   points "reachable" by the kernel (proportional to the pixel value) and also in opposite direction - from those other points to the current
		   pointer (where kernel is centered) with the value proportional to that othre point  */

								max_i= (int) (sigma*nSigma+1);
								k=1.0/(2.0*sigma*sigma);
								if (max_i>=gauss.length) max_i=gauss.length-1;
								sum=-0.5; // 0 is counted twice
								for (i=0; i<=max_i; i++) {
									gauss[i]=Math.exp(-k*i*i);
									sum+= gauss[i]; // could use - more errors for small values of gamma 1/Math.sqrt(2*Math.PI*sigma*sigma)
								}
								scale=0.5/sum;
								for (i=0; i<=max_i; i++) gauss[i]*=scale;
								yk0=-max_i; if (yk0<(y0-iy)) yk0=y0-iy;
								yk1= max_i; if (yk1>=(y1-iy)) yk1=y1-iy-1;
								xk0=-max_i; if (xk0<(x0-ix)) xk0=x0-ix;
								xk1= max_i; if (xk1>=(x1-ix)) xk1=x1-ix-1;

								for (iky=yk0;iky<=yk1;iky++) {
									gy=gauss[Math.abs(iky)]/2; // Extra /2 because we'll calculate the convolution twice from the [ix,iy] and to [ix,iy]
									for (ikx=xk0;ikx<=xk1;ikx++) {
										index=(iy+iky)*size+ix+ikx;
										g=gy*gauss[Math.abs(ikx)];
										result[index]+=d*g;
										result[iy*size+ix]+=pixels[index]*g;

									}
								}
							}
						}
					}
				}
				return result;
			}
			
			
			
		/* ======================================================================== */
		/* find ellipse approximating section of the PSF, scale ellipse and use it as a mask to remove PSF far wings */
			private double [] cutPSFWings (double [] psf_pixels, // direct PSF function, square array, may be proportionally larger than reversed
					double cutoff_energy, // fraction of energy in the pixels to be used
					double ellipse_scale,
					double min_mask_threshold, // zero output element if elliptical Gauss mask is below this threshold
					String title,
					int globalDebugLevel)
			{
				int psf_size=(int)Math.sqrt(psf_pixels.length);
				double [] masked_psf=new double[psf_size*psf_size];
				int  [][]selection=   findClusterOnPSF(psf_pixels, cutoff_energy, title, globalDebugLevel);
				double [] ellipse_coeff=findEllipseOnPSF(psf_pixels,  selection,    title, globalDebugLevel);
				int ix,iy;
				double x,y,r2;
				int indx=0;
				double k2=1/ellipse_scale/ellipse_scale;
				double m;

				for (iy=0;iy<psf_size;iy++) {
					y=(iy-psf_size/2)-ellipse_coeff[1];  // scale to the original psf (and ellipse_coeff)
					for (ix=0;ix<psf_size;ix++) {
						x=(ix-psf_size/2)-ellipse_coeff[0]; // scale to the original psf (and ellipse_coeff)
						r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
						m=Math.exp(-k2*r2);
						masked_psf[indx]=(m>=min_mask_threshold)?(psf_pixels[indx]*Math.exp(-k2*r2)):0.0;
						indx++;
					}
				}

				if (globalDebugLevel>2) {
					ImageProcessor ip_ellipse = new FloatProcessor(psf_size,psf_size);
					float [] ellipsePixels = new float [psf_size*psf_size];
					indx=0;
					for (iy=0;iy<psf_size;iy++) {
						y=(iy-psf_size/2)+ellipse_coeff[1];  // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
						for (ix=0;ix<psf_size;ix++) {
							x=(ix-psf_size/2)+ellipse_coeff[0]; // scale to the original psf (and ellipse_coeff), move center opposite to that of direct kernel (psf)
							r2=ellipse_coeff[2]*x*x+ellipse_coeff[3]*y*y+ellipse_coeff[4]*x*y;
							m=Math.exp(-k2*r2);
							ellipsePixels[indx++]=(float)((m>=min_mask_threshold)?(Math.exp(-k2*r2)):0.0);
						}
					}
					ip_ellipse.setPixels(ellipsePixels);
					ip_ellipse.resetMinAndMax();
					ImagePlus imp_ellipse= new ImagePlus(title+"_PSFWINGS-MASK_"+cutoff_energy+"-"+ellipse_scale, ip_ellipse);
					imp_ellipse.show();
				}
				return masked_psf;
			}


		/* ======================================================================== */
			private double PSFAtXY(double [] pixels, int size, double x, double y) {
				int ix=(int) Math.round(x);
				int iy=(int) Math.round(y);
				if      (ix <  -size/2) ix=-size/2;
				else if (ix >=  size/2) ix= size/2-1;
				if      (iy <  -size/2) iy=-size/2;
				else if (iy >=  size/2) iy= size/2-1;
				int index=size* (size/2 + iy)+ size/2 + ix;
				if ((index<0) || (index > pixels.length)) {
					System.out.println("PSFAtXY error, x="+IJ.d2s(x,0)+" y="+IJ.d2s(y,0)+ " index="+(size*(size/2 + (int) Math.round(y))+ size/2 + (int) Math.round(x))+ " pixels.length="+pixels.length);
				}
				return pixels[index];
			}
		/* ======================================================================== */

			private double contrastAtXY(int sign, double [] pixels, int size, double x, double y, double [][] g, double [] cache) {
				int ir= (int) Math.round(0.2*Math.min(Math.max(Math.abs(g[0][0]),Math.abs(g[1][0])),Math.max(Math.abs(g[0][1]),Math.abs(g[1][1])))); // sample at square 1 1/2x1/2 of the grid "square"

				int ix=(int) Math.round(x);
				int iy=(int) Math.round(y);
				if      (ix <  -size/2) ix=-size/2;
				else if (ix >=  size/2) ix= size/2-1;
				if      (iy <  -size/2) iy=-size/2;
				else if (iy >=  size/2) iy= size/2-1;
				int index= size* (size/2 + iy)+ size/2 + ix;
				//  if ((cache!=null) && (cache[index]>=0)) return sign*cache[index];
				if ((cache!=null) && (cache[index]>=0)) return cache[index];
				double rslt=0.0;
				int i,j;
				for (i=-ir;i<=ir;i++) for (j=-ir;j<=ir;j++) {
					rslt+=     PSFAtXY(pixels,size,j+ix,i+iy) -
					0.25* (PSFAtXY(pixels,size,j+ix+(g[0][0]+ g[1][0])/2  ,i+iy+(g[0][1]+ g[1][1])/2)+
							PSFAtXY(pixels,size,j+ix+(g[0][0]- g[1][0])/2  ,i+iy+(g[0][1]- g[1][1])/2)+
							PSFAtXY(pixels,size,j+ix-(g[0][0]+ g[1][0])/2  ,i+iy-(g[0][1]+ g[1][1])/2)+
							PSFAtXY(pixels,size,j+ix-(g[0][0]- g[1][0])/2  ,i+iy-(g[0][1]- g[1][1])/2));

				}
				rslt=rslt*sign;
				cache[index] = (rslt>0.0)?rslt:0.0;
				return rslt/ir/ir;
			}


		/* ======================================================================== */
		/* create aray (to be used with variableGaussBlurr() ) of per-pixel sigma values for gauss blur, proportional to distance from the specified center */
			private double [] createSigmasRadius (double []pixels, // input square pixel array, preferrably having many exact zeros (they will be skipped)
					double sigmaToRadius, // sigma is proportional to the distance from the center
					double xc, // model PSF center X-coordinate (in pixels[] units, from the center of the array )
					double yc, // same for Y
					int WOICenterX, // window of interest in pixels[] array - do not generate data outside it
					int WOICenterY, // 
					int WOIWidth, //
					int WOIHeight) {
				int size = (int) Math.sqrt(pixels.length);
				double [] sigmas =new double [size*size];
				int x0= (size-WOIWidth)/2 +WOICenterX;
				int y0= (size-WOIHeight)/2+WOICenterY;
				int x1=x0+WOIWidth;
				int y1=x0+WOIHeight;
				int i,ix,iy;
				double r,x,y;
				for (i=0;i<sigmas.length;i++) sigmas[i]=0.0;
				if (x0<0) x0=0; if (x1>size) x1=size; if (y0<0) y0=0; if (y1>size) y1=size;
				for (iy=0;iy<size;iy++) {
					y=(iy-size/2)-yc;
					for (ix=0;ix<size;ix++) {
						x=(ix-size/2)-xc;
						r=Math.sqrt(x*x+y*y);
						//        sigma=r*sigmaToRadius;
						//        sigma=r*r/radiusSigma;
						//        sigmas[iy*size+ix]=(r*sigmaToRadius)+1;
						sigmas[iy*size+ix]=(r*sigmaToRadius);
					}
				}


				return sigmas;
			}

			/* calculates ellipse (with the center at DC) that interpolates area of the points defined by flooding from the initial center,
			so total energy is cutoff_energy fraction
			returns {x0,y0,a,b,c} , where a*x^2+b*y^2 + c*x*y=r^2 , so r^2 can be used for a window that removes high far pixels
			distribute the whol mass at the ends of short and long ellipse axis
			u^2/Ru^2+V^2/Rv^2=1, u=cos(a)*x+sin(a)*y, v=-sin(a)*x+cos(a)*y
			c=cos(a), s=sin(a), S0=sum(f(x,y), SX2=sum(f(x,y)*(x-x0)*(x-x0)),SY2=sum(f(x,y)*(y-y0)*(y-y0)), SXY=sum(f(x,y)*(x-x0)*(y-y0))
			"effective" squared radius (to be used in Gaussian)
			r2= u^2/Ru^2+V^2/Rv^2
			r2= 1/Ru^2 * 1/Rv^2 * (x^2*(c^2*Rv^2+s^2*Ru^2)+y^2*(c^2*Ru^2+s^2*Rv^2)+2*x*y*c*s*(Rv^2-Ru^2)

			SX2/S0=1/2* ((c*Ru)^2 + (s*Rv)^2)         =1/2*(c^2*Ru^2 + s^2*Rv^2)
			SY2/S0=1/2* ((s*Ru)^2 + (c*Rv)^2)         =1/2*(c^2*Rv^2 + s^2*Ru^2)
			SXY/S0=1/2* ((c*Ru)*(s*Ru)-(c*Rv)*(s*rv)) =1/2*(c*s*(Ru^2 -Rv^2))

			r2= 1/Ru^2 * 1/Rv^2 * (x^2*(2*SY2/S0))+y^2*(2*SX2/S0)-2*2*x*y*(SXY/S0)

			SX2/S0+SY2/S0= 1/2*(Ru^2 + Rv^2)
			Ru^2+Rv^2= 2*(SX2+SY2)/S0
			Ru^2-Rv^2= 2* SXY /S0

			Ru^2=(SX2+SY2+SXY)/S0
			Rv^2=(SX2+SY2-SXY)/S0

			r2= a* x^2*+b*y^2+c*x*y
			a=  1/Ru^2 * 1/Rv^2 * (2*SY2/S0)
			b=  1/Ru^2 * 1/Rv^2 * (2*SX2/S0)
			c= -1/Ru^2 * 1/Rv^2 * (4*SXY/S0)
				 */
				private double [] findEllipseOnPSF(
						double []         psf,   // Point Spread Function (may be off-center)
						int    [][] selection, // 0/1 - selected/not selected
						String          title,
						int globalDebugLevel) {
					int i,j;
					double x,y;
					int size=(int) Math.sqrt(psf.length);
					double SX=0.0;
					double SY=0.0;
					double SX2=0.0;
					double SY2=0.0;
					double SXY=0.0;
					double S0=0.0;
					double d; //,k;
					//	double area=0; // selection area
			/* find centyer */

					for (i=0;i<size;i++) {
						y=i-size/2;
						for (j=0;j<size;j++) if (selection[i][j]>0){
							x=j-size/2;
							d=psf[i*size+j];
							S0+=d;
							SX+=x*d;
							SY+=y*d;
							//			area+=1.0;
						}
					}
					double centerX=SX/S0;
					double centerY=SY/S0;
					if (globalDebugLevel>5) {
						//		System.out.println("findEllipseOnPSF: title="+title+" area="+area+" S0="+S0+" SX="+SX+" SY="+SY+" centerX="+centerX+" centerY="+centerY);
						System.out.println("findEllipseOnPSF: title="+title+" S0="+S0+" SX="+SX+" SY="+SY+" centerX="+centerX+" centerY="+centerY);
					}

			/* second pass (could all be done in a single) */
					SX2=0.0;
					SY2=0.0;
					SXY=0.0;
					for (i=0;i<size;i++) {
						y=i-size/2-centerY;
						for (j=0;j<size;j++) if (selection[i][j]>0){
							x=j-size/2-centerX;
							d=psf[i*size+j];
							SX2+=x*x*d;
							SY2+=y*y*d;
							SXY+=x*y*d;
						}
					}
					if (globalDebugLevel>5) {
						System.out.println("findEllipseOnPXF: title="+title+" SX2="+SX2+" SY2="+SY2+" SXY="+SXY);
					}
					/*
			Ru^2=(SX2+SY2+SXY)/S0
			Rv^2=(SX2+SY2-SXY)/S0

			r2= a* x^2*+b*y^2+c*x*y
			a=  1/Ru^2 * 1/Rv^2 * (2*SY2/S0)
			b=  1/Ru^2 * 1/Rv^2 * (2*SX2/S0)
			c= -1/Ru^2 * 1/Rv^2 * (4*SXY/S0)
					 */
					double Ru2=(SX2+SY2+SXY)/S0;
					double Rv2=(SX2+SY2-SXY)/S0;
					double [] result = {centerX,
							centerY,
							1/Ru2 * 1/Rv2 * (2*SY2/S0),
							1/Ru2 * 1/Rv2 * (2*SX2/S0),
							-1/Ru2 * 1/Rv2 * (4*SXY/S0)};
					//	k=Math.PI*Math.PI/(2.0*S0*area*area);
					//	double [] result = {centerX,centerY,k*SY2,k*SX2,-2*k*SXY};
					if (globalDebugLevel>3) {
						System.out.println("findEllipseOnPS: title="+title+" x0="+result[0]+" y0="+result[1]+" a="+result[2]+" b="+result[3]+" c="+result[4]);
					}
					return result;
				}

		
		/* ======================================================================== */
		/* finds cluster on the PSF (with the center at specidfied point)  by flooding from the specified center, so total energy is cutoff_energy fraction
		returns integer array (same dimensions as input) with 1 - selected, 0 - not selected
		cutoff_energy: if positive - specifies fraction of total energy, if negative -cutoff_energy is the minimal value of the pixel to be included 
		UPDATE: follows gradient from the start point to a local maximum if "cutoff_energy" is negative" */
			private int [][] findClusterOnPSF(
					double []        psf, // PSF function, square array
					double cutoff_energy, // fraction of energy in the pixels to be used
					String         title,
					int            globalDebugLevel) {
				int size=(int) Math.sqrt(psf.length);
				return findClusterOnPSF(psf,          // PSF function, square array
						cutoff_energy, // fraction of energy in the pixels to be used
						size/2,        // X0
						size/2,        // Y0
						title,
						globalDebugLevel);
			}



			private int [][] findClusterOnPSF(
					double []        psf, // PSF function, square array
					double cutoff_energy, // fraction of energy in the pixels to be used (or minimal level if it is negative)
					int           startX,  // location of a start point, x-coordinate
					int           startY,  // location of a start point, y-coordinate
					String         title,
					int globalDebugLevel) {
				int i,j;
				int ix,iy,ix1,iy1,maxX, maxY;
				List <Integer> pixelList=new ArrayList<Integer>(100);
				Integer Index;
				int size=(int) Math.sqrt(psf.length);
				int [][]clusterMap=new int[size][size];
				double full_energy=0.0;
				int [][] dirs={{-1,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1}};
				ix=startX;
				iy=startY;
				Index=iy*size + ix;
				double maxValue=psf[Index];
		/* Make ix,iy to start from the maximal value on PSF */
				Index=0;
				for (i=0;i<size;i++) for (j=0;j<size;j++) {
					full_energy+=psf[Index];
					clusterMap[i][j]=0;
					if (psf[Index]>maxValue){
						maxValue=psf[Index];
						ix=j;
						iy=i;
					}
					Index++;
				}
				boolean noThreshold=(cutoff_energy<=0);
				double threshold=full_energy*((cutoff_energy>0)?cutoff_energy:1.0); // no limit for negative values of cutoff_energy
				double minValue=0.0; // no limit if total energy is controlled
				double cluster_energy=0.0;
				int clusterSize=0;
				boolean noNew=true;
				if (cutoff_energy<=0) { // find nearest local maximum following gradient
					ix=startX;
					iy=startY;
					maxValue=psf[iy*size + ix];
					for (noNew=false;noNew==false;){
						noNew=true;
						for (j=0;j<dirs.length;j++) if (((iy > 0 )        || (dirs[j][1]>=0)) &&
								((iy < (size-1) ) || (dirs[j][1]<=0)) &&
								((ix > 0 )        || (dirs[j][0]>=0)) &&
								((ix < (size-1) ) || (dirs[j][0]<=0))){
							ix1= ix+dirs[j][0];
							iy1= iy+dirs[j][1];
							if (psf[iy1*size+ix1]>maxValue) {
								noNew=false;
								maxValue= psf[iy1*size+ix1];
								ix=ix1;
								iy=iy1;
								break;
							}
						}
					}
					minValue=maxValue*(-cutoff_energy);
				}
		//
		if (globalDebugLevel>1)		System.out.println("findClusterOnPSF: full_energy="+full_energy+" minValue="+minValue+" maxValue="+maxValue);
		if (globalDebugLevel>1)		System.out.println("findClusterOnPSF: ix="+ix+" iy="+iy);
				maxX=0;
				maxY=0;
				int listIndex;
				Index=iy*size + ix;
				pixelList.clear();
				pixelList.add (Index);
				clusterSize++;
				clusterMap[iy][ix]=1;
				cluster_energy+=psf[Index];
				noNew=true;
				while ((pixelList.size()>0) &&  (noThreshold || (cluster_energy<threshold) )) { // will break from the loop if  (psf[Index] <minValue)
		/* Find maximal new neighbor */
					maxValue=0.0;
					listIndex=0;
					while (listIndex<pixelList.size()) {
						Index=pixelList.get(listIndex);
						iy=Index/size;
						ix=Index%size;
						noNew=true;
						for (j=0;j<8;j++) if (((iy > 0 ) || (dirs[j][1]>=0)) && ((iy < (size-1) ) || (dirs[j][1]<=0))){
							ix1=(ix+dirs[j][0]+size) % size;
							iy1= iy+dirs[j][1];
							if (clusterMap[iy1][ix1]==0) {
								noNew=false;
								if (psf[iy1*size+ix1]>maxValue) {
									maxValue= psf[iy1*size+ix1];
									maxX=ix1;
									maxY=iy1;
								}
							}
						}
						if (noNew) pixelList.remove(listIndex);  //  remove current list element
						else       listIndex++;     // increase list index
					}
					if (maxValue==0.0) { // Should
						if (!noThreshold) System.out.println("findClusterOnPSF: - should not get here - no points around >0, and threshold is not reached yet.");
						break;
					}
		/* Add this new point to the list */
					if (psf[Index]<minValue) break; // break if the condition was value, not total energy
					Index=maxY*size + maxX;
					pixelList.add (Index);
					clusterSize++;
					clusterMap[maxY][maxX]=1;
					cluster_energy+=psf[Index];

				} // end of while ((pixelList.size()>0) &&  (cluster_energy<threshold))
				if (globalDebugLevel>3)   System.out.println("findClusterOnPSF: cluster size is "+clusterSize);
				if (globalDebugLevel>6) {
					ImageProcessor ip2 = new FloatProcessor(size,size);
					float [] floatPixels = new float [size*size];
					for (i=0;i<floatPixels.length;i++) {
						floatPixels[i]=(float) psf[i];
					}
					ip2.setPixels(floatPixels);
					ip2.resetMinAndMax();
					ImagePlus imp2= new ImagePlus(title+"_PSF1_"+cutoff_energy, ip2);
					imp2.show();
				}
				if (globalDebugLevel>5) {
					ImageProcessor ip = new FloatProcessor(size,size);
					float [] floatPixels = new float [size*size];
					for (i=0;i<floatPixels.length;i++) {
						floatPixels[i]=(float) clusterMap[i/size][i%size];
					}
					ip.setPixels(floatPixels);
					ip.resetMinAndMax();
					ImagePlus imp= new ImagePlus(title+"_PSF-SEL_"+cutoff_energy, ip);
					imp.show();
				}
				return clusterMap;
			}

			/* ======================================================================== */
			/* ======================================================================== */
			/* calculates 2x2 matrix that converts two pairs of vectors: u2=M*u1, v2=M*v1*/

			/* ======================================================================== */

			/* ======================================================================== */

				private  int [] convert2d_1d(int [][] pixels){
					int i,j;
					int width=pixels[0].length;
					int [] rslt=new int[pixels.length*pixels[0].length];
					for (i=0;i<pixels.length;i++) for (j=0;j<width;j++) rslt[i*width+j]=pixels[i][j];
					return rslt;
				}

			/* pixels should be a square array, zero is in the center (/center+0.5 for even dimensions) */
//				private  double [] calcCentroidFromCenter(double [] pixels) {return calcCentroidFromCenter(pixels, (int[]) null, 0.0);}
				private  double [] calcCentroidFromCenter(double [] pixels, // square pixel array
						int[][] mask, // integer mask -0 - don't use this pixel, 1 - use it
						double refLevel) { // subtract this fraction of maximal level from all pixels
					return calcCentroidFromCenter(pixels, convert2d_1d(mask), refLevel);
				}
				private  double [] calcCentroidFromCenter(double [] pixels, // square pixel array
						int[] mask, // integer mask -0 - don't use this pixel, 1 - use it
						double refLevel) { // subtract this fraction of maximal leve from all pixels
					int size = (int) Math.sqrt ( pixels.length);
					int c= size/2;
					double S0=0.0;
					double SX=0.0;
					double SY=0.0;
					double x,y,p;
					int i,j,indx;
					double maxValue = 0.0;
					if (refLevel>0.0) for (i=0;i<pixels.length;i++) if (((mask==null) || (mask[i]>0)) && (pixels[i] > maxValue)) maxValue=pixels[i];

					double minValue=refLevel*maxValue;

					for (i=0;i<size;i++) {
						y=i-c;
						for (j=0;j<size;j++) {
							indx=i*size+j;
							if ((mask==null) || (mask[indx]>0)) {
								x=j-c;
								p=pixels[indx]-minValue;
								if (p>0.0) { // with mask mis-match there could be negative total mask
									S0+=p;
									SX+=p*x;
									SY+=p*y;
								}
							}
						}
					}
					double [] result={SX/S0,SY/S0};
					return result;
				}

		/* ======================================================================== */
		private double [] binPSF(double [] pixels,
				double [][] g,
				int outSize,
				//		int      decimate,     // sub-pixel decimation 
				double minContrast,
				double [] centerXY,    // coordinates (x,y) of the center point (will be alway subtracted)
				double[] symmXY,       // coordinates (x,y) of the center of symmetry (to combine with 180 if enabled by symm180)
				int pass,              // mostly for debug purposes
				String title,
				boolean debug,
				int globalDebugLevel) {
			int multiple=2;         // 0 - use each pixel once, 1 - add first negatives (4), 2 - second positives()4)
			int pixelSize=(int) Math.sqrt(pixels.length);
			int halfOutSize=outSize/2;
			int indx,i,j,outIndex,ix,iy;
			double x,y,xc,yc,uc,vc,u,v,p,q,d, du, dv, dp,dq, xr,yr, overThreshold;
			int np,nq;
			int PSF_sign=1;
			double [] contrastCache=new double[pixelSize*pixelSize];
			double [] debugPixels=null;
			if (debug)  debugPixels=new double[pixelSize*pixelSize];

			double det_g=g[0][0]*g[1][1]-g[0][1]*g[1][0];
			double [][] xy2uv= {{-2.0*g[0][1]/det_g,  2.0*g[0][0]/det_g},
					{-2.0*g[1][1]/det_g,  2.0*g[1][0]/det_g}};
			double [][] uv2xy= matrix2x2_scale(matrix2x2_invert(xy2uv),2); // real pixels are twice
			double [] pixelsPSF       =new double [outSize*outSize];  
			int    [] pixelsPSFCount  =new int    [outSize*outSize];
			double [] pixelsPSFWeight =new double [outSize*outSize];  
			double [] center=centerXY;
			for (i=0;i<contrastCache.length;i++) {
				contrastCache[i]=-1.0;
			}
			double threshold=minContrast*contrastAtXY(1, pixels, pixelSize, 0.0, 0.0,  g, contrastCache);
			if (debug)  {
				System.out.println("binPSF title="+title+" g[0][0]="+IJ.d2s(g[0][0],4)+" g[0][1]="+IJ.d2s(g[0][1],4));
				System.out.println("binPSF title="+title+" g[1][0]="+IJ.d2s(g[1][0],4)+" g[1][1]="+IJ.d2s(g[1][1],4));
				System.out.println("  center[0]="+center[0]+"  center[1]="+center[1]);
				//		System.out.println("  decimate="+decimate+"  threshold="+threshold);
				System.out.println("  threshold="+threshold);
			}

			if (center==null) {
				center = new double[2];
				center[0]=0.0;
				center[1]=0.0;
			}
			for (i=0;i<pixelsPSF.length;i++) {
				pixelsPSF[i]=0.0;
				pixelsPSFCount[i]=0;
				pixelsPSFWeight[i]=0.0;
			}

			for (indx=0;indx<pixels.length;indx++) {
				y= indx / pixelSize- pixelSize/2;
				x= indx % pixelSize- pixelSize/2;
				u= xy2uv[0][0]*x + xy2uv[0][1]*y;
				v= xy2uv[1][0]*x + xy2uv[1][1]*y;
				p=u+v;
				q=u-v;
				np=(int)Math.floor((1+p)/2);
				nq=(int)Math.floor((1+q)/2);
				//if (debug)  debugPixels[indx]=(int)Math.floor((1+q)/2);
	/* see if the point is in the cell of positive or negative OTF instance */
				PSF_sign= (((np + nq) & 1)==0)?1:-1;
	/* find x,y coordinates of the center of the cell */
				uc=0.5*(np+nq);
				vc=0.5*(np-nq);
				//xc=g[0][0]*uc + g[1][0]*vc;
				//yc=g[0][1]*uc + g[1][1]*vc;

				yc=-g[0][0]*uc - g[1][0]*vc;
				xc= g[0][1]*uc + g[1][1]*vc;


				//if (debug) debugPixels[indx]=p/2-Math.round(p/2);

	/* See if this cell has enough contrast */
				overThreshold=contrastAtXY(PSF_sign,pixels, pixelSize, xc,yc,  g, contrastCache);
				//if (debug) debugPixels[indx]=overThreshold;
				if (overThreshold<threshold) {
					if (debug) debugPixels[indx]=0.0;
					//if (debug) debugPixels[indx]=yc;
				} else {
					//if (debug) debugPixels[indx]=yc;

	/* Do binning itself here */
					d=PSF_sign*PSFAtXY(pixels, pixelSize, x,y);

	/* map to the segment around 0,0 */        
					dp=p/2-Math.round(p/2);
					dq=q/2-Math.round(q/2);
	/* dp, dq are between +/- 0.5 - use them for Hamming windowing -NOT HERE, moved later*/
					du=(dp+dq)/2;
					dv=(dp-dq)/2;

	/* bin this point to the center and some (positive) duplicates if enabled */
					for (i=-(multiple/2); i<=(multiple/2); i++) for (j=-(multiple/2); j<=(multiple/2); j++) {
						xr= uv2xy[0][0]*(j+du) + uv2xy[0][1]*(i+dv);
						yr= uv2xy[1][0]*(j+du) + uv2xy[1][1]*(i+dv);
						xr= Math.round(xr-center[0]);
						yr= Math.round(yr-center[1]);
	/* does it fit into output array ? */
						if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
							outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
							pixelsPSFCount[outIndex]++;
							pixelsPSF[outIndex]+=d*overThreshold;
							pixelsPSFWeight[outIndex]+=overThreshold;
						}
					}
	/* bin this to center-symmetrical point if enabled */
					if (symmXY!=null) {
						for (i=-(multiple/2); i<=(multiple/2); i++) for (j=-(multiple/2); j<=(multiple/2); j++) {
							xr= uv2xy[0][0]*(j+du) + uv2xy[0][1]*(i+dv);
							yr= uv2xy[1][0]*(j+du) + uv2xy[1][1]*(i+dv);
							xr= Math.round(symmXY[0]*2.0-xr-center[0]);
							yr= Math.round(symmXY[1]*2.0-yr-center[1]);
							//does it fit into output array ?
							if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
								outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
								pixelsPSFCount[outIndex]++;
								pixelsPSF[outIndex]+=d*overThreshold;
								pixelsPSFWeight[outIndex]+=overThreshold;
							}
						}
					}
	/* Now bin this point to the negative duplicates if enabled (debug feature). Normally it will be skipped */
					if (multiple>0) for (i=-((multiple+1)/2); i<((multiple+1)/2); i++) for (j=-((multiple+1)/2); j<((multiple+1)/2); j++) {
						xr= uv2xy[0][0]*(j+du+0.5) + uv2xy[0][1]*(i+dv+0.5);
						yr= uv2xy[1][0]*(j+du+0.5) + uv2xy[1][1]*(i+dv+0.5);
						xr= Math.round(xr-center[0]);
						yr= Math.round(yr-center[1]);
						//does it fit into output array ?
						if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
							outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
							pixelsPSFCount[outIndex]++;
							pixelsPSF[outIndex]-=d*overThreshold;
							pixelsPSFWeight[outIndex]+=overThreshold;
						}
					}
	/* bin this to center-symmetrical point if enabled */
	/* Now bin this point to the negative duplicates if enabled (debug feature). Normally it will be skipped */
					if (symmXY!=null) {
						if (multiple>0) for (i=-((multiple+1)/2); i<((multiple+1)/2); i++) for (j=-((multiple+1)/2); j<((multiple+1)/2); j++) {
							xr= uv2xy[0][0]*(j+du+0.5) + uv2xy[0][1]*(i+dv+0.5);
							yr= uv2xy[1][0]*(j+du+0.5) + uv2xy[1][1]*(i+dv+0.5);
							xr= Math.round(symmXY[0]*2.0-xr-center[0]);
							yr= Math.round(symmXY[1]*2.0-yr-center[1]);
							//does it fit into output array ?
							if ((yr>=-halfOutSize) && (yr<halfOutSize) && (xr>=-halfOutSize) && (xr<halfOutSize)) {
								outIndex=outSize*(outSize/2+ ((int) yr))+(outSize/2)+((int) xr);
								pixelsPSFCount[outIndex]++;
								pixelsPSF[outIndex]+=d*overThreshold;
								pixelsPSFWeight[outIndex]+=overThreshold;
							}
						}
					}
				}
			}


			for (i=0;i<pixelsPSF.length;i++) {
				if (pixelsPSFWeight[i]>0.0) pixelsPSF[i]/=pixelsPSFWeight[i];
			}
	/* Interpolate  missing points (pixelsPSFCount[i]==0) */

			for (i=0;i<pixelsPSF.length;i++) if (pixelsPSFWeight[i]==0.0){
				iy=i/outSize;
				ix=i%outSize;
				if ((ix>0)&&(ix<(outSize-1))&&(iy>0)&&(iy<(outSize-1))) {
					if ((pixelsPSFWeight[(iy-1)*outSize+ix  ]>0.0) &&
							(pixelsPSFWeight[(iy+1)*outSize+ix  ]>0.0) &&
							(pixelsPSFWeight[(iy  )*outSize+ix-1]>0.0) &&
							(pixelsPSFWeight[(iy  )*outSize+ix+1]>0.0)) {
						if (globalDebugLevel>5) System.out.println("Interpolating missing OTF point at x="+ix+" y="+iy);
						pixelsPSF[i]=
							0.25*(pixelsPSF[(iy-1)*outSize+ix  ]+
									pixelsPSF[(iy+1)*outSize+ix  ]+
									pixelsPSF[(iy  )*outSize+ix-1]+
									pixelsPSF[(iy  )*outSize+ix+1]);
					}
				}
			}
	/* optionally show original array with masked out low-contrast cells */
			if ((globalDebugLevel>2) && (pass==1))  SDFA_INSTANCE.showArrays(pixelsPSF, title+"_Used-PSF");
			if (debug) {
				SDFA_INSTANCE.showArrays(debugPixels, title+"_mask_PSF");
				double [] doublePixelsPSFCount=new double [pixelsPSF.length];
				for (j=0;j<doublePixelsPSFCount.length;j++) doublePixelsPSFCount[j]=(double)pixelsPSFCount[j];
				SDFA_INSTANCE.showArrays(doublePixelsPSFCount, title+"_PSF_bin_count");
				SDFA_INSTANCE.showArrays(pixelsPSFWeight,      title+"_PSF_bin_weight");
				double [] doubleContrastCache=new double [contrastCache.length];
				for (j=0;j<doubleContrastCache.length;j++) doubleContrastCache[j]=(double)((contrastCache[j]>=0.0)?contrastCache[j]:-0.00001);
				SDFA_INSTANCE.showArrays(doubleContrastCache,  title+"_ContrastCache");
			}
			return pixelsPSF;
		}

	
	
	
	
	/* ======================================================================== */
	/* Create a Thread[] array as large as the number of processors available.
		 * From Stephan Preibisch's Multithreading.java class. See:
		 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
		 */
		private Thread[] newThreadArray(int maxCPUs) {
			int n_cpus = Runtime.getRuntime().availableProcessors();
			if (n_cpus>maxCPUs)n_cpus=maxCPUs;
			return new Thread[n_cpus];
		}
	/* Start all given threads and wait on each of them until all are done.
		 * From Stephan Preibisch's Multithreading.java class. See:
		 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
		 */
		public static void startAndJoin(Thread[] threads)
		{
			for (int ithread = 0; ithread < threads.length; ++ithread)
			{
				threads[ithread].setPriority(Thread.NORM_PRIORITY);
				threads[ithread].start();
			}

			try
			{   
				for (int ithread = 0; ithread < threads.length; ++ithread)
					threads[ithread].join();
			} catch (InterruptedException ie)
			{
				throw new RuntimeException(ie);
			}
		}
		/* === Parameter classes === */
		public static class MultiFilePSF {
			public double  overexposedMaxFraction; // allowed fraction of the overexposed pixels in the PSF kernel measurement area
			public double  weightOnBorder=0.5;
			public double  radiusDiffLow= 0.1; // do not remove partial kernel cell if radius differs from average less than by this fraction
			public double  radiusDiffHigh=0.25;  // remove this cell even if it is the only one
			public double  shiftToRadiusContrib=1.0; // Center shift (in pixels) addition to the difference relative to radius difference (in pixels)
			public double  sharpBonusPower=2.0; // increase weight of the "sharp" kernels by dividing weight by radius to this power
			public double  maxFracDiscardWorse=0.1; // discard up to this fraction of samples that have larger radius (i.e. falling on the target seam that may only make PSF larger)
			public double  maxFracDiscardAll=0.5; // continue removing outlayers (combined radius and shift), removing not more that this fraction (including maxFracDiscardWorse)
			public double  internalBonus=1.0;    // cell having 8 around will "seem" twice better than having none (radiusDiff* twice higher)
			public double  validateThreshold;      // fraction of full PSF "energy"
			public boolean validateShowEllipse;    // show ellipse parameters of partial PSF arrays
			public boolean showWeights;            // show image indicating frame coverage
			public boolean fillMissing;            // replace missing kernels with neighbors
			public MultiFilePSF (
					double  overexposedMaxFraction,
					double  weightOnBorder,          
					double  radiusDiffLow, // do not remove partial kernel cell if radius differs from average less than by this fraction
					double  radiusDiffHigh,  // remove this cell even if it is the only one
					double  shiftToRadiusContrib, // Center shift (in pixels) addition to the difference relative to radius difference (in pixels)
					double  sharpBonusPower, // increase weight of the "sharp" kernels by dividing weight by radius to this power
					double  maxFracDiscardWorse, // discard up to this fraction of samples that have larger radius (i.e. falling on the target seam that may only make PSF larger)
					double  maxFracDiscardAll,  // continue removing outlayers (combined radius and shift), removing not more that this fraction (including maxFracDiscardWorse)
					double  internalBonus,
					double  validateThreshold,
					boolean validateShowEllipse,
					boolean showWeights,
					boolean fillMissing
			) {
				this.overexposedMaxFraction=overexposedMaxFraction;
				this.weightOnBorder=weightOnBorder;          
				this.radiusDiffLow=radiusDiffLow; // do not remove partial kernel cell if radius differs from average less than by this fraction
				this.radiusDiffHigh=radiusDiffHigh;  // remove this cell even if it is the only one
				this.shiftToRadiusContrib=shiftToRadiusContrib; // Center shift (in pixels) addition to the difference relative to radius difference (in pixels)
				this.sharpBonusPower=sharpBonusPower; // increase weight of the "sharp" kernels by dividing weight by radius to this power
				this.maxFracDiscardWorse=maxFracDiscardWorse; // discard up to this fraction of samples that have larger radius (i.e. falling on the target seam that may only make PSF larger)
				this.maxFracDiscardAll=maxFracDiscardAll; // continue removing outlayers (combined radius and shift), removing not more that this fraction (including maxFracDiscardWorse)
				this.internalBonus=internalBonus;
				this.validateThreshold=validateThreshold;
				this.validateShowEllipse=validateShowEllipse;
				this.showWeights=showWeights;
				this.fillMissing=fillMissing;
			}

			public void setProperties(String prefix,Properties properties){
				properties.setProperty(prefix+"overexposedMaxFraction",this.overexposedMaxFraction+"");
				properties.setProperty(prefix+"weightOnBorder",this.weightOnBorder+"");
				properties.setProperty(prefix+"radiusDiffLow",this.radiusDiffLow+"");
				properties.setProperty(prefix+"radiusDiffHigh",this.radiusDiffHigh+"");
				properties.setProperty(prefix+"shiftToRadiusContrib",this.shiftToRadiusContrib+"");
				properties.setProperty(prefix+"sharpBonusPower",this.sharpBonusPower+"");
				properties.setProperty(prefix+"maxFracDiscardWorse",this.maxFracDiscardWorse+"");
				properties.setProperty(prefix+"maxFracDiscardAll",this.maxFracDiscardAll+"");
				properties.setProperty(prefix+"internalBonus",this.internalBonus+"");
				properties.setProperty(prefix+"validateThreshold",this.validateThreshold+"");
				properties.setProperty(prefix+"validateShowEllipse",this.validateShowEllipse+"");
				properties.setProperty(prefix+"showWeights",this.showWeights+"");
				properties.setProperty(prefix+"fillMissing",this.fillMissing+"");
			}

			public void setProperties(String prefix,ImagePlus properties){
				properties.setProperty(prefix+"overexposedMaxFraction",this.overexposedMaxFraction+"");
				properties.setProperty(prefix+"weightOnBorder",this.weightOnBorder+"");
				properties.setProperty(prefix+"radiusDiffLow",this.radiusDiffLow+"");
				properties.setProperty(prefix+"radiusDiffHigh",this.radiusDiffHigh+"");
				properties.setProperty(prefix+"shiftToRadiusContrib",this.shiftToRadiusContrib+"");
				properties.setProperty(prefix+"sharpBonusPower",this.sharpBonusPower+"");
				properties.setProperty(prefix+"maxFracDiscardWorse",this.maxFracDiscardWorse+"");
				properties.setProperty(prefix+"maxFracDiscardAll",this.maxFracDiscardAll+"");
				properties.setProperty(prefix+"internalBonus",this.internalBonus+"");
				properties.setProperty(prefix+"validateThreshold",this.validateThreshold+"");
				properties.setProperty(prefix+"validateShowEllipse",this.validateShowEllipse+"");
				properties.setProperty(prefix+"showWeights",this.showWeights+"");
				properties.setProperty(prefix+"fillMissing",this.fillMissing+"");
			}
			
			public void getProperties(String prefix,Properties properties){
				if (properties.getProperty(prefix+"overexposedMaxFraction")!=null) this.overexposedMaxFraction=Double.parseDouble(properties.getProperty(prefix+"overexposedMaxFraction"));
				if (properties.getProperty(prefix+"weightOnBorder")!=null) this.weightOnBorder=Double.parseDouble(properties.getProperty(prefix+"weightOnBorder"));
				if (properties.getProperty(prefix+"radiusDiffLow")!=null) this.radiusDiffLow=Double.parseDouble(properties.getProperty(prefix+"radiusDiffLow"));
				if (properties.getProperty(prefix+"radiusDiffHigh")!=null) this.radiusDiffHigh=Double.parseDouble(properties.getProperty(prefix+"radiusDiffHigh"));
				if (properties.getProperty(prefix+"shiftToRadiusContrib")!=null)
					this.shiftToRadiusContrib=Double.parseDouble(properties.getProperty(prefix+"shiftToRadiusContrib"));
				if (properties.getProperty(prefix+"sharpBonusPower")!=null)
					this.sharpBonusPower=Double.parseDouble(properties.getProperty(prefix+"sharpBonusPower"));
				if (properties.getProperty(prefix+"maxFracDiscardWorse")!=null)
					this.maxFracDiscardWorse=Double.parseDouble(properties.getProperty(prefix+"maxFracDiscardWorse"));
				if (properties.getProperty(prefix+"maxFracDiscardAll")!=null)
					this.maxFracDiscardAll=Double.parseDouble(properties.getProperty(prefix+"maxFracDiscardAll"));
				if (properties.getProperty(prefix+"internalBonus")!=null) this.internalBonus=Double.parseDouble(properties.getProperty(prefix+"internalBonus"));
				if (properties.getProperty(prefix+"validateThreshold")!=null)this.validateThreshold=Double.parseDouble(properties.getProperty(prefix+"validateThreshold"));
				if (properties.getProperty(prefix+"validateShowEllipse")!=null)this.validateShowEllipse=Boolean.parseBoolean(properties.getProperty(prefix+"validateShowEllipse"));
				if (properties.getProperty(prefix+"showWeights")!=null)this.showWeights=Boolean.parseBoolean(properties.getProperty(prefix+"showWeights"));
				if (properties.getProperty(prefix+"fillMissing")!=null)this.fillMissing=Boolean.parseBoolean(properties.getProperty(prefix+"fillMissing"));
				
			}

		}

    public static class AberrationParameters{
    	public String sourceDirectory="";
    	public String partialKernelDirectory="";
    	public String psfKernelDirectory="";
    	public String aberrationsKernelDirectory="";
    	public boolean autoRestore;
    	public String calibrationPath="";
    	public String strategyPath="";
    	public String gridPath="";
    	public String sensorsPath="";
    	public boolean autoRestoreSensorOverwriteOrientation=false; // dangerous! true; // overwrite camera parameters from sensor calibration files
    	public boolean autoRestoreSensorOverwriteDistortion= false; // dangerous! true; // overwrite camera parameters from sensor calibration files
		public boolean autoReCalibrate=true; // Re-calibrate grids on autoload
		public boolean autoReCalibrateIgnoreLaser=false; // "Ignore laser pointers on recalibrate"
    	public boolean autoFilter=false; // true;
    	public boolean trustEnabled= true; // mark all enabled images as hintedMatch=2 on Read Calibration
    	public boolean noMessageBoxes=true;
    	public boolean overwriteResultFiles=false;
    	public boolean partialToReprojected=true; // Use reprojected grid for partial kernel calculation (false - use extracted)
    	public boolean partialCorrectSensor=true; // Apply sensor correction to the projected grid
    	public int     seriesNumber=0;
    	public boolean allImages;
    	public String sourcePrefix="";
    	public String sourceSuffix=".tiff";
    	public String partialPrefix="partial-";
    	public String partialSuffix=".ppsf-tiff";
    	public String psfPrefix="direct-psf-";
    	public String psfSuffix=".psf-tiff";
    	public String interpolatedPSFPrefix="interpolated-psf-";
    	public String interpolatedPSFSuffix=".ipsf-tiff";
    	public String aberrationsPrefix="kernel-";
    	public String aberrationsSuffix=".kernel-tiff";
    	public boolean [] selectedChannels=null;
    	
    	
    	
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"sourceDirectory",this.sourceDirectory);
			properties.setProperty(prefix+"partialKernelDirectory",this.partialKernelDirectory);
			properties.setProperty(prefix+"psfKernelDirectory",this.psfKernelDirectory);
			properties.setProperty(prefix+"aberrationsKernelDirectory",this.aberrationsKernelDirectory);
			
			properties.setProperty(prefix+"autoRestore",this.autoRestore+"");
			properties.setProperty(prefix+"calibrationPath",this.calibrationPath);
			properties.setProperty(prefix+"strategyPath",this.strategyPath);
			properties.setProperty(prefix+"gridPath",this.gridPath);
			properties.setProperty(prefix+"sensorsPath",this.sensorsPath);
			properties.setProperty(prefix+"autoRestoreSensorOverwriteOrientation",this.autoRestoreSensorOverwriteOrientation+"");
			properties.setProperty(prefix+"autoRestoreSensorOverwriteDistortion",this.autoRestoreSensorOverwriteDistortion+"");
			properties.setProperty(prefix+"autoReCalibrate",this.autoReCalibrate+"");
			properties.setProperty(prefix+"autoReCalibrateIgnoreLaser",this.autoReCalibrateIgnoreLaser+"");
			properties.setProperty(prefix+"autoFilter",this.autoFilter+"");
			properties.setProperty(prefix+"trustEnabled",this.trustEnabled+"");
			properties.setProperty(prefix+"noMessageBoxes",this.noMessageBoxes+"");
			properties.setProperty(prefix+"overwriteResultFiles",this.overwriteResultFiles+"");
			properties.setProperty(prefix+"partialToReprojected",this.partialToReprojected+"");
			properties.setProperty(prefix+"partialCorrectSensor",this.partialCorrectSensor+"");
			
			
			properties.setProperty(prefix+"seriesNumber",this.seriesNumber+"");
			properties.setProperty(prefix+"allImages",this.allImages+"");

			properties.setProperty(prefix+"sourcePrefix",this.sourcePrefix);
			properties.setProperty(prefix+"sourceSuffix",this.sourceSuffix);
			properties.setProperty(prefix+"partialPrefix",this.partialPrefix);
			properties.setProperty(prefix+"partialSuffix",this.partialSuffix);
			properties.setProperty(prefix+"psfPrefix",this.psfPrefix);
			properties.setProperty(prefix+"psfSuffix",this.psfSuffix);
			properties.setProperty(prefix+"interpolatedPSFPrefix",this.interpolatedPSFPrefix);
			properties.setProperty(prefix+"interpolatedPSFSuffix",this.interpolatedPSFSuffix);
			properties.setProperty(prefix+"aberrationsPrefix",this.aberrationsPrefix);
			properties.setProperty(prefix+"aberrationsSuffix",this.aberrationsSuffix);
			if (this.selectedChannels!=null){
				String sSelectedChannels="";
				for (int i=0;i<this.selectedChannels.length;i++) sSelectedChannels+= selectedChannels[i]?"+":"-";
				properties.setProperty(prefix+"selectedChannels",sSelectedChannels);
			}

		
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"sourceDirectory")!=null)	           this.sourceDirectory=properties.getProperty(prefix+"sourceDirectory");
			if (properties.getProperty(prefix+"partialKernelDirectory")!=null)     this.partialKernelDirectory=properties.getProperty(prefix+"partialKernelDirectory");
			if (properties.getProperty(prefix+"psfKernelDirectory")!=null)         this.psfKernelDirectory=properties.getProperty(prefix+"psfKernelDirectory");
			if (properties.getProperty(prefix+"aberrationsKernelDirectory")!=null) this.aberrationsKernelDirectory=properties.getProperty(prefix+"aberrationsKernelDirectory");
			
			if (properties.getProperty(prefix+"autoRestore")!=null)                this.autoRestore=Boolean.parseBoolean(properties.getProperty(prefix+"autoRestore"));
			if (properties.getProperty(prefix+"calibrationPath")!=null)            this.calibrationPath=properties.getProperty(prefix+"calibrationPath");
			if (properties.getProperty(prefix+"strategyPath")!=null)               this.strategyPath=properties.getProperty(prefix+"strategyPath");
			if (properties.getProperty(prefix+"gridPath")!=null)                   this.gridPath=properties.getProperty(prefix+"gridPath");
			if (properties.getProperty(prefix+"sensorsPath")!=null)                this.sensorsPath=properties.getProperty(prefix+"sensorsPath");
			if (properties.getProperty(prefix+"autoRestoreSensorOverwriteOrientation")!=null)
				this.autoRestoreSensorOverwriteOrientation=Boolean.parseBoolean(properties.getProperty(prefix+"autoRestoreSensorOverwriteOrientation"));
			if (properties.getProperty(prefix+"autoRestoreSensorOverwriteDistortion")!=null)
				this.autoRestoreSensorOverwriteDistortion=Boolean.parseBoolean(properties.getProperty(prefix+"autoRestoreSensorOverwriteDistortion"));
			
			if (properties.getProperty(prefix+"autoReCalibrate")!=null)            this.autoReCalibrate=Boolean.parseBoolean(properties.getProperty(prefix+"autoReCalibrate"));
			if (properties.getProperty(prefix+"autoReCalibrateIgnoreLaser")!=null) this.autoReCalibrateIgnoreLaser=Boolean.parseBoolean(properties.getProperty(prefix+"autoReCalibrateIgnoreLaser"));
			if (properties.getProperty(prefix+"autoFilter")!=null)                 this.autoFilter=Boolean.parseBoolean(properties.getProperty(prefix+"autoFilter"));
			if (properties.getProperty(prefix+"trustEnabled")!=null)               this.trustEnabled=Boolean.parseBoolean(properties.getProperty(prefix+"trustEnabled"));
			if (properties.getProperty(prefix+"noMessageBoxes")!=null)             this.noMessageBoxes=Boolean.parseBoolean(properties.getProperty(prefix+"noMessageBoxes"));
			if (properties.getProperty(prefix+"overwriteResultFiles")!=null)       this.overwriteResultFiles=Boolean.parseBoolean(properties.getProperty(prefix+"overwriteResultFiles"));
			if (properties.getProperty(prefix+"partialToReprojected")!=null)       this.partialToReprojected=Boolean.parseBoolean(properties.getProperty(prefix+"partialToReprojected"));
			if (properties.getProperty(prefix+"partialCorrectSensor")!=null)       this.partialCorrectSensor=Boolean.parseBoolean(properties.getProperty(prefix+"partialCorrectSensor"));
			
			
			if (properties.getProperty(prefix+"seriesNumber")!=null)               this.seriesNumber=Integer.parseInt(properties.getProperty(prefix+"seriesNumber"));
			if (properties.getProperty(prefix+"allImages")!=null)                  this.allImages=Boolean.parseBoolean(properties.getProperty(prefix+"allImages"));
			if (properties.getProperty(prefix+"sourcePrefix")!=null)	      this.sourcePrefix=properties.getProperty(prefix+"sourcePrefix");
			if (properties.getProperty(prefix+"sourceSuffix")!=null)	      this.sourceSuffix=properties.getProperty(prefix+"sourceSuffix");
			if (properties.getProperty(prefix+"partialPrefix")!=null)	      this.partialPrefix=properties.getProperty(prefix+"partialPrefix");
			if (properties.getProperty(prefix+"partialSuffix")!=null)         this.partialSuffix=properties.getProperty(prefix+"partialSuffix");
			if (properties.getProperty(prefix+"psfPrefix")!=null)	          this.psfPrefix=properties.getProperty(prefix+"psfPrefix");
			if (properties.getProperty(prefix+"psfSuffix")!=null)	          this.psfSuffix=properties.getProperty(prefix+"psfSuffix");
			if (properties.getProperty(prefix+"interpolatedPSFPrefix")!=null) this.interpolatedPSFPrefix=properties.getProperty(prefix+"interpolatedPSFPrefix");
			if (properties.getProperty(prefix+"interpolatedPSFSuffix")!=null) this.interpolatedPSFSuffix=properties.getProperty(prefix+"interpolatedPSFSuffix");
			if (properties.getProperty(prefix+"aberrationsPrefix")!=null)     this.aberrationsPrefix=properties.getProperty(prefix+"aberrationsPrefix");
			if (properties.getProperty(prefix+"aberrationsSuffix")!=null)     this.aberrationsSuffix=properties.getProperty(prefix+"aberrationsSuffix");
			
			
			if (properties.getProperty(prefix+"selectedChannels")!=null){
				String sSelectedChannels=properties.getProperty(prefix+"selectedChannels");
				this.selectedChannels=new boolean[sSelectedChannels.length()];
				for (int i=0;i<this.selectedChannels.length;i++) selectedChannels[i]= sSelectedChannels.charAt(i)=='+';
			}
		}
		/**
		 * calibration files paths
		 * @param distortions Distortion class instance
		 * @param combine when true - return the configured paths if current is not set, false - return null
		 *  for the paths that are not set or did not change from configured
		 * @return array of 4 paths
		 */
		
		public String [] currentConfigPaths(Distortions distortions, boolean combine){
    		String currentCalibrationPath=null;
    		String currentStrategyPath=null;
    		String currentGridPath=null;
    		String currentSensorsPath=null;
    		if (distortions!=null) {
    			if (distortions.fittingStrategy!=null){
    				currentStrategyPath=distortions.fittingStrategy.pathName;
//        			System.out.println("currentConfigPaths():currentStrategyPath="+((currentStrategyPath==null)?"null":currentStrategyPath));
        			if (distortions.fittingStrategy.distortionCalibrationData!=null) {
        				currentCalibrationPath=distortions.fittingStrategy.distortionCalibrationData.pathName;
        			} else {
//            			System.out.println("currentConfigPaths():distortions.fittingStrategy.distortionCalibrationData==null");
        			}
    			} else {
//        			System.out.println("currentConfigPaths():distortions.fittingStrategy==null");
    				
    			}
    			currentGridPath=distortions.patternParameters.pathName;
//    			System.out.println("currentConfigPaths():currentGridPath="+((currentGridPath==null)?"null":currentGridPath));
    			currentSensorsPath=distortions.getSensorPath(-1);
//    			System.out.println("currentConfigPaths():currentSensorsPath="+((currentSensorsPath==null)?"null":currentSensorsPath));
    		} else {
    			System.out.println("currentConfigPaths():distortions==null");
    		}
    		if ((currentCalibrationPath==null) || (currentCalibrationPath.equals(this.calibrationPath)) || (currentCalibrationPath.length()==0))
    			currentCalibrationPath=combine?this.calibrationPath:null;
    		if ((currentStrategyPath==   null) || (currentStrategyPath.equals   (this.strategyPath))    || (currentStrategyPath.length()==0))
    			currentStrategyPath=combine?this.strategyPath:null;
    		if ((currentGridPath==       null) || (currentGridPath.equals       (this.gridPath))        || (currentGridPath.length()==0))
    			currentGridPath=combine?this.gridPath:null;
    		if ((currentSensorsPath==    null) || (currentSensorsPath.equals    (this.sensorsPath))     || (currentSensorsPath.length()==0))
    			currentSensorsPath=combine?this.sensorsPath:null;
    		String [] result={
    				currentCalibrationPath,
    				currentStrategyPath,
    				currentGridPath,
    				currentSensorsPath
    		};
			return result;
		}
		public String [] autoLoadPaths(){
    		String [] result={
    				this.calibrationPath,
    				this.strategyPath,
    				this.gridPath,
    				this.sensorsPath
    		};
			return result;
		}
		public boolean [] getChannelSelection(Distortions distortions){
			if (distortions==null) return null;
		   	int numChannels=distortions.fittingStrategy.distortionCalibrationData.getNumChannels(); // number of used channels
    		if (this.selectedChannels==null) {
    			this.selectedChannels=new boolean[1];
    			this.selectedChannels[0]=true;
    		}
    		if (this.selectedChannels.length!=numChannels){
    			boolean [] tmp=this.selectedChannels;
    			this.selectedChannels=new boolean[numChannels];
    			for (int i=0;i<numChannels;i++){
    				this.selectedChannels[i]=(i<tmp.length)?tmp[i]:tmp[tmp.length-1];
    			}
    		}
    		return this.selectedChannels;
			
		}
		public boolean selectChannelsToProcess(String title, Distortions distortions) {
    		boolean [] newSelecttion=getChannelSelection(distortions).clone(); //java.lang.NullPointerException
			if (newSelecttion==null) return false;
			int numChannels=newSelecttion.length;
    		while (true) {
    			GenericDialog gd = new GenericDialog(title);
    			for (int i=0;i<numChannels;i++) gd.addCheckbox("channel "+i, newSelecttion[i]);
    			gd.enableYesNoCancel("OK", "All like channel 0");
    			WindowTools.addScrollBars(gd);
    			gd.showDialog();
    			if (gd.wasCanceled()) return false;
    			for (int i=0;i<numChannels;i++) newSelecttion[i]=gd.getNextBoolean();
    			if (gd.wasOKed()){
    				for (int i=0;i<numChannels;i++) this.selectedChannels[i]=newSelecttion[i];
    				return true;
    			} else {
    				for (int i=1;i<numChannels;i++) newSelecttion[i]=newSelecttion[0];
    			}
    		}
		}

    	public boolean showDialog(String title, Distortions distortions) { 
    		String [] currentConfigs;
    		String []nulls={null,null,null,null};
    		currentConfigs=(distortions!=null)?currentConfigPaths(distortions, false):nulls;
    		GenericDialog gd = new GenericDialog(title);
    		gd.addStringField("Source files directory", this.sourceDirectory, 60);
    		gd.addCheckbox("Select source directory", false);
    		gd.addStringField("Partial kernels directory", this.partialKernelDirectory, 60);
    		gd.addCheckbox("Select partial kernels directory", false);
    		gd.addStringField("Combined kernels directory", this.psfKernelDirectory, 60);
    		gd.addCheckbox("Select combined kernsls directory", false);
    		gd.addStringField("Aberrations kernels directory", this.aberrationsKernelDirectory, 60);
    		gd.addCheckbox("Select aberrations kernels directory", false);
    		gd.addCheckbox("Supress non-essential message boxes", this.noMessageBoxes);
    		gd.addCheckbox("Overwrite result files if they exist", this.overwriteResultFiles);
    		gd.addCheckbox("Use reprojected grids for partial kernel calculation (false - extracted grids)", this.partialToReprojected);
    		gd.addCheckbox("Apply sensor correction during for partial kernel calculation", this.partialCorrectSensor);
    		    		
    		gd.addNumericField("Fitting series number to use for image selection", this.seriesNumber,0);
    		gd.addCheckbox("Process all enabled image files (false - use selected fitting series)", this.allImages);
    		gd.addMessage("===== Autoload options (when restoring configuration) =====");
    		gd.addCheckbox("Autoload additional files on \"Restore\"", this.autoRestore);
    		gd.addCheckbox("Overwrite all (including position/orientation) SFE parameters from the sensor calibration files (at auto-load) DANGEROUS!", this.autoRestoreSensorOverwriteOrientation);
    		gd.addCheckbox("Overwrite SFE distortion parameters from the sensor calibration files (at auto-load) DANGEROUS!", this.autoRestoreSensorOverwriteDistortion);
    		
    		gd.addCheckbox("Re-calibrate grids on autoload", this.autoReCalibrate);
    		gd.addCheckbox("Ignore laser pointers on recalibrate", this.autoReCalibrateIgnoreLaser);
    		gd.addCheckbox("Filter grids after restore", this.autoFilter);
    		gd.addCheckbox("Trust enabled images on input (mark as hintedGrid=2)", this.trustEnabled);
   		
    		gd.addMessage("Calibration: "+(((this.calibrationPath==null) || (this.calibrationPath.length()==0))?"not configured ":(this.calibrationPath+" "))+
    				((currentConfigs[0]!=null)?("(current: "+currentConfigs[0]+")"):("") ));
    		gd.addMessage("Strategy: "+(((this.strategyPath==null) || (this.strategyPath.length()==0))?"not configured ":(this.strategyPath+" "))+
    				((currentConfigs[1]!=null)?("(current: "+currentConfigs[1]+")"):("") ));
    		gd.addMessage("Pattern grid: "+(((this.gridPath==null) || (this.gridPath.length()==0))?"not configured ":(this.gridPath+" "))+
    				((currentConfigs[2]!=null)?("(current: "+currentConfigs[2]+")"):("") ));
    		gd.addMessage("Sensors(one of): "+(((this.sensorsPath==null) || (this.sensorsPath.length()==0))?"not configured ":(this.sensorsPath+" "))+
    				((currentConfigs[3]!=null)?("(current: "+currentConfigs[3]+")"):("") ));
    		gd.addCheckbox("Update configured (for auto-load) paths from current ones", false);
    		gd.addMessage("Filename prefixes/suffixes:");
    		gd.addStringField("Source files prefix",             this.sourcePrefix, 40);
    		gd.addStringField("Source files suffix",             this.sourceSuffix, 40);
    		gd.addStringField("Partial kernels prefix",          this.partialPrefix, 40);
    		gd.addStringField("Partial kernels suffix",          this.partialSuffix, 40);
    		gd.addStringField("Combined kernels prefix",         this.psfPrefix, 40);
    		gd.addStringField("Combined kernels suffix",         this.psfSuffix, 40);
    		gd.addStringField("Interpolated kernels prefix",     this.interpolatedPSFPrefix, 40);
    		gd.addStringField("Interpolated kernels suffix",     this.interpolatedPSFSuffix, 40);
    		gd.addStringField("Inverted (final) kernels prefix", this.aberrationsPrefix, 40);
    		gd.addStringField("Inverted (final) kernels suffix", this.aberrationsSuffix, 40);
    		
    		gd.addCheckbox("Select channels to process", true);
    		
    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		this.sourceDirectory=       gd.getNextString();
    		if (gd.getNextBoolean()) selectSourceDirectory(false, this.sourceDirectory, false); 
    		this.partialKernelDirectory=gd.getNextString();
    		if (gd.getNextBoolean()) selectPartialKernelDirectory(false, this.partialKernelDirectory, false); 
    		this.psfKernelDirectory=gd.getNextString();
    		if (gd.getNextBoolean()) selectPSFKernelDirectory(false, this.psfKernelDirectory, false); 
    		this.aberrationsKernelDirectory=gd.getNextString();
    		if (gd.getNextBoolean()) selectAberrationsKernelDirectory(false, this.aberrationsKernelDirectory, false);
    		this.noMessageBoxes=        gd.getNextBoolean();
    		this.overwriteResultFiles=  gd.getNextBoolean();
    		this.partialToReprojected=  gd.getNextBoolean();
    		this.partialCorrectSensor=  gd.getNextBoolean();
    		this.seriesNumber=    (int) gd.getNextNumber();
    		this.allImages=             gd.getNextBoolean();
    		this.autoRestore=           gd.getNextBoolean();
    		this.autoRestoreSensorOverwriteOrientation= gd.getNextBoolean();
    		this.autoRestoreSensorOverwriteDistortion= gd.getNextBoolean();
    		this.autoReCalibrate=           gd.getNextBoolean();
    		this.autoReCalibrateIgnoreLaser=gd.getNextBoolean();
    		this.autoFilter=            gd.getNextBoolean();
    		this.trustEnabled=          gd.getNextBoolean();
    		
    		if (gd.getNextBoolean()) {
    			if (currentConfigs[0]!=null) this.calibrationPath=currentConfigs[0];
    			if (currentConfigs[1]!=null) this.strategyPath=   currentConfigs[1];
    			if (currentConfigs[2]!=null) this.gridPath=       currentConfigs[2];
    			if (currentConfigs[3]!=null) this.sensorsPath=    currentConfigs[3];
    		}
    		
    		this.sourcePrefix=          gd.getNextString();
    		this.sourceSuffix=          gd.getNextString();
    		this.partialPrefix=         gd.getNextString();
    		this.partialSuffix=         gd.getNextString();
    		this.psfPrefix=             gd.getNextString();
    		this.psfSuffix=             gd.getNextString();
    		this.interpolatedPSFPrefix= gd.getNextString();
    		this.interpolatedPSFSuffix= gd.getNextString();
    		this.aberrationsPrefix=     gd.getNextString();
    		this.aberrationsSuffix=     gd.getNextString();
    		if (gd.getNextBoolean()) selectChannelsToProcess("Select channels to process", distortions);
    		return true;
    	}
    	public String selectSourceDirectory(boolean smart, String defaultPath, boolean newAllowed) { // normally newAllowed=false
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Source (acquired from the camera) image directory", // title
    				"Select source directory", // button
    				null, // filter
    				defaultPath); // this.sourceDirectory);
    		if (dir!=null) this.sourceDirectory=dir;
    		return dir;
    	}
    	public String selectPartialKernelDirectory(boolean smart, String defaultPath, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Partial PSF directory", // title
    				"Select partial PSF files directory", // button
    				null, // filter
    				defaultPath); //this.sourceDirectory);
    		if (dir!=null) this.partialKernelDirectory=dir;
    		return dir;
    	}
    	public String selectPSFKernelDirectory(boolean smart, String defaultPath, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Combined direct PSF kernel directory", // title
    				"Select combined kernel directory", // button
    				null, // filter
    				defaultPath); //this.sourceDirectory);
    		if (dir!=null) this.psfKernelDirectory=dir;
    		return dir;
    	}
    	public String selectAberrationsKernelDirectory(boolean smart, String defaultPath, boolean newAllowed) {
    		String dir= CalibrationFileManagement.selectDirectory(
    				smart,
    				newAllowed, // save  
    				"Aberrations kernel directory", // title
    				"Select aberrations kernel directory", // button
    				null, // filter
    				defaultPath); //this.sourceDirectory);
    		if (dir!=null) this.aberrationsKernelDirectory=dir;
    		return dir;
    	}

    	
    }
	
	
	
	public static class ColorComponents {
		public boolean [] colorsToCorrect=    new boolean[6];
		public int        referenceComponent; // component to calculate lateral chromatic from (0 - G1, 1 - R, 2 - B, 3 - G2,4 - diagonal greens, 5 - checker greens)
		public boolean    equalizeGreens;   // equalize 2 greens in Bayer mosaic
		public static String [] componentColorNames={"green1","red","blue","green2", "greens (diagonal)", "greens (checker)"};
		public static String [] stackColorNames={"red","green","blue"};
		public String getColorName(int i) {return  componentColorNames[i];}
		public String getStackColorName(int i) {return  stackColorNames[i];}

		public ColorComponents (
				boolean green1,
				boolean red,
				boolean blue,
				boolean green2,
				boolean diagonal, // both greens combined in a 45-degree rotated array 
				boolean checker,   // both greens combined in a checkerboard pattern
				int        referenceComponent,
				boolean    equalizeGreens
		) {
			this.colorsToCorrect[0]=green1;
			this.colorsToCorrect[1]=red;
			this.colorsToCorrect[2]=blue;
			this.colorsToCorrect[3]=green2;
			this.colorsToCorrect[4]=diagonal;
			this.colorsToCorrect[5]=checker;
			this.referenceComponent=referenceComponent;
			this.equalizeGreens=equalizeGreens;
		}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"green1",this.colorsToCorrect[0]+"");
			properties.setProperty(prefix+"red",this.colorsToCorrect[1]+"");
			properties.setProperty(prefix+"blue",this.colorsToCorrect[2]+"");
			properties.setProperty(prefix+"green2",this.colorsToCorrect[3]+"");
			properties.setProperty(prefix+"diagonal",this.colorsToCorrect[4]+"");
			properties.setProperty(prefix+"checker",this.colorsToCorrect[5]+"");
			properties.setProperty(prefix+"referenceComponent",this.referenceComponent+"");
			properties.setProperty(prefix+"equalizeGreens",this.equalizeGreens+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.colorsToCorrect[0]=Boolean.parseBoolean(properties.getProperty(prefix+"green1"));
			this.colorsToCorrect[1]=Boolean.parseBoolean(properties.getProperty(prefix+"red"));
			this.colorsToCorrect[2]=Boolean.parseBoolean(properties.getProperty(prefix+"blue"));
			this.colorsToCorrect[3]=Boolean.parseBoolean(properties.getProperty(prefix+"green2"));
			this.colorsToCorrect[4]=Boolean.parseBoolean(properties.getProperty(prefix+"diagonal"));
			this.colorsToCorrect[5]=Boolean.parseBoolean(properties.getProperty(prefix+"checker"));
			this.referenceComponent=Integer.parseInt(properties.getProperty(prefix+"referenceComponent"));
			this.equalizeGreens=Boolean.parseBoolean(properties.getProperty(prefix+"equalizeGreens"));
		}
	}

	public static class OTFFilterParameters {
		public double deconvInvert;
		public double zerofreqSize;
		public double smoothPS;
		public double thresholdHigh;
		public double thresholdLow;

		public OTFFilterParameters(
				double deconvInvert,
				double zerofreqSize,
				double smoothPS,
				double thresholdHigh,
				double thresholdLow) {
			this.deconvInvert = deconvInvert;
			this.zerofreqSize = zerofreqSize;
			this.smoothPS = smoothPS;
			this.thresholdHigh = thresholdHigh;
			this.thresholdLow = thresholdLow;
		}
        public OTFFilterParameters clone() {
        	return new OTFFilterParameters(
        			this.deconvInvert,
        			this.zerofreqSize,
        			this.smoothPS,
        			this.thresholdHigh,
        			this.thresholdLow);
        }
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"deconvInvert",this.deconvInvert+"");
			properties.setProperty(prefix+"zerofreqSize",this.zerofreqSize+"");
			properties.setProperty(prefix+"smoothPS",this.smoothPS+"");
			properties.setProperty(prefix+"thresholdHigh",this.thresholdHigh+"");
			properties.setProperty(prefix+"thresholdLow",this.thresholdLow+"");
		}
		public void getProperties(String prefix,Properties properties){
			this.deconvInvert=Double.parseDouble(properties.getProperty(prefix+"deconvInvert"));
			this.zerofreqSize=Double.parseDouble(properties.getProperty(prefix+"zerofreqSize"));
			this.smoothPS=Double.parseDouble(properties.getProperty(prefix+"smoothPS"));
			this.thresholdHigh=Double.parseDouble(properties.getProperty(prefix+"thresholdHigh"));
			this.thresholdLow=Double.parseDouble(properties.getProperty(prefix+"thresholdLow"));
		}

	}


	public static class PSFParameters {
		public double minContrast;
		public double windowFrac;
		public boolean useWindow;
		public boolean symm180;
		public boolean ignoreChromatic;
		public boolean absoluteCenter;
		public double smoothSeparate;
		public double topCenter;
		public double sigmaToRadius;
		public double wingsEnergy;
		public double wingsEllipseScale;
		public double minDefinedArea;   // minimal (weighted) fraction of the defined patter pixels in the FFT area
		public boolean approximateGrid; // approximate grid with polynomial 
		public boolean centerPSF;       // Center PSF by modifying phase 
		public double mask1_sigma;
		public double mask1_threshold;
		public double gaps_sigma;
		public double mask_denoise;
		

		public PSFParameters(double minContrast,
				double windowFrac,
				boolean useWindow,
				boolean symm180,
				boolean ignoreChromatic,
				boolean absoluteCenter,
				double smoothSeparate,
				double topCenter,
				double sigmaToRadius,
				double wingsEnergy,
				double wingsEllipseScale,
				double minDefinedArea, // minimal (weighted) fraction of the defined patter pixels in the FFT area
				boolean approximateGrid, // approximate grid with polynomial
				boolean centerPSF,       // Center PSF by modifying phase
				double mask1_sigma,
				double mask1_threshold,
				double gaps_sigma,
				double mask_denoise

		) {
			this.minContrast = minContrast;
			this.windowFrac = windowFrac;
			this.useWindow = useWindow;
			this.symm180 = symm180;
			this.ignoreChromatic = ignoreChromatic;
			this.absoluteCenter=absoluteCenter;
			this.smoothSeparate = smoothSeparate;
			this.topCenter = topCenter;
			this.sigmaToRadius = sigmaToRadius;
			this.wingsEnergy = wingsEnergy;
			this.wingsEllipseScale = wingsEllipseScale;
			this.minDefinedArea = minDefinedArea; // minimal (weighted) fraction of the defined patter pixels in the FFT area
			this.approximateGrid = approximateGrid; // approximate grid with polynomial 
			this.centerPSF = centerPSF; // approximate grid with polynomial 
			this.mask1_sigma = mask1_sigma;
			this.mask1_threshold = mask1_threshold;
			this.gaps_sigma=gaps_sigma;
			this.mask_denoise=mask_denoise;


		}
        public PSFParameters clone(){
        	return new PSFParameters(
        			this.minContrast,
        			this.windowFrac,
        			this.useWindow,
        			this.symm180,
        			this.ignoreChromatic,
        			this.absoluteCenter,
        			this.smoothSeparate,
        			this.topCenter,
        			this.sigmaToRadius,
        			this.wingsEnergy,
        			this.wingsEllipseScale,
        			this.minDefinedArea, // minimal (weighted) fraction of the defined patter pixels in the FFT area
        			this.approximateGrid, // approximate grid with polynomial
        			this.centerPSF, // approximate grid with polynomial
        			this.mask1_sigma,
        			this.mask1_threshold,
        			this.gaps_sigma,
        			this.mask_denoise
                      );
        }
        public void setProperties(String prefix,Properties properties){
        	properties.setProperty(prefix+"minContrast",this.minContrast+"");
        	properties.setProperty(prefix+"windowFrac",this.windowFrac+"");
        	properties.setProperty(prefix+"useWindow",this.useWindow+"");
        	properties.setProperty(prefix+"symm180",this.symm180+"");
        	properties.setProperty(prefix+"ignoreChromatic",this.ignoreChromatic+"");
        	properties.setProperty(prefix+"absoluteCenter",this.absoluteCenter+"");
        	properties.setProperty(prefix+"smoothSeparate",this.smoothSeparate+"");
        	properties.setProperty(prefix+"topCenter",this.topCenter+"");
        	properties.setProperty(prefix+"sigmaToRadius",this.sigmaToRadius+"");
        	properties.setProperty(prefix+"wingsEnergy",this.wingsEnergy+"");
        	properties.setProperty(prefix+"wingsEllipseScale",this.wingsEllipseScale+"");
        	properties.setProperty(prefix+"minDefinedArea",this.minDefinedArea+"");
        	properties.setProperty(prefix+"approximateGrid",this.approximateGrid+"");
        	properties.setProperty(prefix+"centerPSF",this.centerPSF+"");
        	properties.setProperty(prefix+"mask1_sigma",this.mask1_sigma+"");
        	properties.setProperty(prefix+"mask1_threshold",this.mask1_threshold+"");
        	properties.setProperty(prefix+"gaps_sigma",this.gaps_sigma+"");
        	properties.setProperty(prefix+"mask_denoise",this.mask_denoise+"");
        }
        public void setProperties(String prefix, ImagePlus properties){
        	properties.setProperty(prefix+"minContrast",this.minContrast+"");
        	properties.setProperty(prefix+"windowFrac",this.windowFrac+"");
        	properties.setProperty(prefix+"useWindow",this.useWindow+"");
        	properties.setProperty(prefix+"symm180",this.symm180+"");
        	properties.setProperty(prefix+"ignoreChromatic",this.ignoreChromatic+"");
        	properties.setProperty(prefix+"absoluteCenter",this.absoluteCenter+"");
        	properties.setProperty(prefix+"smoothSeparate",this.smoothSeparate+"");
        	properties.setProperty(prefix+"topCenter",this.topCenter+"");
        	properties.setProperty(prefix+"sigmaToRadius",this.sigmaToRadius+"");
        	properties.setProperty(prefix+"wingsEnergy",this.wingsEnergy+"");
        	properties.setProperty(prefix+"wingsEllipseScale",this.wingsEllipseScale+"");
        	properties.setProperty(prefix+"minDefinedArea",this.minDefinedArea+"");
        	properties.setProperty(prefix+"approximateGrid",this.approximateGrid+"");
        	properties.setProperty(prefix+"centerPSF",this.centerPSF+"");
        	properties.setProperty(prefix+"mask1_sigma",this.mask1_sigma+"");
        	properties.setProperty(prefix+"mask1_threshold",this.mask1_threshold+"");
        	properties.setProperty(prefix+"gaps_sigma",this.gaps_sigma+"");
        	properties.setProperty(prefix+"mask_denoise",this.mask_denoise+"");
        }

		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"minContrast")!=null)       this.minContrast=Double.parseDouble(properties.getProperty(prefix+"minContrast"));
			if (properties.getProperty(prefix+"windowFrac")!=null)        this.windowFrac=Double.parseDouble(properties.getProperty(prefix+"windowFrac"));
			if (properties.getProperty(prefix+"useWindow")!=null)         this.useWindow=Boolean.parseBoolean(properties.getProperty(prefix+"useWindow"));
			if (properties.getProperty(prefix+"symm180")!=null)           this.symm180=Boolean.parseBoolean(properties.getProperty(prefix+"symm180"));
			if (properties.getProperty(prefix+"ignoreChromatic")!=null)   this.ignoreChromatic=Boolean.parseBoolean(properties.getProperty(prefix+"ignoreChromatic"));
			if (properties.getProperty(prefix+"absoluteCenter")!=null)   this.absoluteCenter=Boolean.parseBoolean(properties.getProperty(prefix+"absoluteCenter"));
			if (properties.getProperty(prefix+"smoothSeparate")!=null)    this.smoothSeparate=Double.parseDouble(properties.getProperty(prefix+"smoothSeparate"));
			if (properties.getProperty(prefix+"topCenter")!=null)         this.topCenter=Double.parseDouble(properties.getProperty(prefix+"topCenter"));
			if (properties.getProperty(prefix+"sigmaToRadius")!=null)     this.sigmaToRadius=Double.parseDouble(properties.getProperty(prefix+"sigmaToRadius"));
			if (properties.getProperty(prefix+"wingsEnergy")!=null)       this.wingsEnergy=Double.parseDouble(properties.getProperty(prefix+"wingsEnergy"));
			if (properties.getProperty(prefix+"wingsEllipseScale")!=null) this.wingsEllipseScale=Double.parseDouble(properties.getProperty(prefix+"wingsEllipseScale"));
			if (properties.getProperty(prefix+"minDefinedArea")!=null)    this.minDefinedArea=Double.parseDouble(properties.getProperty(prefix+"minDefinedArea"));
			if (properties.getProperty(prefix+"approximateGrid")!=null)   this.approximateGrid=Boolean.parseBoolean(properties.getProperty(prefix+"approximateGrid"));
			if (properties.getProperty(prefix+"centerPSF")!=null)         this.centerPSF=Boolean.parseBoolean(properties.getProperty(prefix+"centerPSF"));
			if (properties.getProperty(prefix+"mask1_sigma")!=null)       this.mask1_sigma=Double.parseDouble(properties.getProperty(prefix+"mask1_sigma"));
			if (properties.getProperty(prefix+"mask1_threshold")!=null)   this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"mask1_threshold"));
			if (properties.getProperty(prefix+"gaps_sigma")!=null)        this.mask1_threshold=Double.parseDouble(properties.getProperty(prefix+"gaps_sigma"));
			if (properties.getProperty(prefix+"mask_denoise")!=null)        this.mask_denoise=Double.parseDouble(properties.getProperty(prefix+"mask_denoise"));
		}
	}

	public static class InverseParameters {
		public int dSize;
		public int rSize;
		public double invertRange;
		public double otfCutoffEnergy;
		public double otfEllipseScale;
		public boolean otfEllipseGauss;
		public double psfCutoffEnergy;
		public double psfEllipseScale;
		public double rpsfMinMaskThreshold;
		public boolean filter;
		public double blurIndividual;
		public double blurDiagonal;
		public double blurChecker;
		public double gaussianSigmaIndividual;
		public double gaussianSigmaDiagonal;
		public double gaussianSigmaChecker;
		public double sigmaScale;
		public double sigmaToRadius;
		public boolean filterDirect;
		public double sigmaScaleDirect;
		public double sigmaToRadiusDirect;

		public InverseParameters(int dSize, int rSize, double invertRange,
				double otfCutoffEnergy, double otfEllipseScale,
				boolean otfEllipseGauss, double psfCutoffEnergy,
				double psfEllipseScale, double rpsfMinMaskThreshold,
				boolean filter, double blurIndividual, double blurDiagonal, double blurChecker,
				double gaussianSigmaIndividual, double gaussianSigmaDiagonal, double gaussianSigmaChecker,
				double sigmaScale, double sigmaToRadius,
				boolean filterDirect, double sigmaScaleDirect, double sigmaToRadiusDirect
				) {
			this.dSize = dSize;
			this.rSize = rSize;
			this.invertRange = invertRange;
			this.otfCutoffEnergy = otfCutoffEnergy;
			this.otfEllipseScale = otfEllipseScale;
			this.otfEllipseGauss = otfEllipseGauss;
			this.psfCutoffEnergy = psfCutoffEnergy;
			this.psfEllipseScale = psfEllipseScale;
			this.rpsfMinMaskThreshold = rpsfMinMaskThreshold;
			this.filter = filter;
			this.blurIndividual = blurIndividual;
			this.blurDiagonal = blurDiagonal;
			this.blurChecker = blurChecker;
			this.gaussianSigmaIndividual = gaussianSigmaIndividual;
			this.gaussianSigmaDiagonal = gaussianSigmaDiagonal;
			this.gaussianSigmaChecker = gaussianSigmaChecker;
			this.sigmaScale = sigmaScale;
			this.sigmaToRadius = sigmaToRadius;
			this.filterDirect=filterDirect;
			this.sigmaScaleDirect=sigmaScaleDirect;
			this.sigmaToRadiusDirect=sigmaToRadiusDirect;
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"dSize",this.dSize+"");
			properties.setProperty(prefix+"rSize",this.rSize+"");
			properties.setProperty(prefix+"invertRange",this.invertRange+"");
			properties.setProperty(prefix+"otfCutoffEnergy",this.otfCutoffEnergy+"");
			properties.setProperty(prefix+"otfEllipseScale",this.otfEllipseScale+"");
			properties.setProperty(prefix+"otfEllipseGauss",this.otfEllipseGauss+"");
			properties.setProperty(prefix+"psfCutoffEnergy",this.psfCutoffEnergy+"");
			properties.setProperty(prefix+"psfEllipseScale",this.psfEllipseScale+"");
			properties.setProperty(prefix+"rpsfMinMaskThreshold",this.rpsfMinMaskThreshold+"");
			properties.setProperty(prefix+"filter",this.filter+"");
			properties.setProperty(prefix+"blurIndividual",this.blurIndividual+"");
			properties.setProperty(prefix+"blurDiagonal",this.blurDiagonal+"");
			properties.setProperty(prefix+"blurChecker",this.blurChecker+"");
			properties.setProperty(prefix+"gaussianSigmaIndividual",this.gaussianSigmaIndividual+"");
			properties.setProperty(prefix+"gaussianSigmaDiagonal",this.gaussianSigmaDiagonal+"");
			properties.setProperty(prefix+"gaussianSigmaChecker",this.gaussianSigmaChecker+"");
			properties.setProperty(prefix+"sigmaScale",this.sigmaScale+"");
			properties.setProperty(prefix+"sigmaToRadius",this.sigmaToRadius+"");
			properties.setProperty(prefix+"filterDirect",this.filterDirect+"");
			properties.setProperty(prefix+"sigmaScaleDirect",this.sigmaScaleDirect+"");
			properties.setProperty(prefix+"sigmaToRadiusDirect",this.sigmaToRadiusDirect+"");
		}

		public void setProperties(String prefix,ImagePlus properties){
			properties.setProperty(prefix+"dSize",this.dSize+"");
			properties.setProperty(prefix+"rSize",this.rSize+"");
			properties.setProperty(prefix+"invertRange",this.invertRange+"");
			properties.setProperty(prefix+"otfCutoffEnergy",this.otfCutoffEnergy+"");
			properties.setProperty(prefix+"otfEllipseScale",this.otfEllipseScale+"");
			properties.setProperty(prefix+"otfEllipseGauss",this.otfEllipseGauss+"");
			properties.setProperty(prefix+"psfCutoffEnergy",this.psfCutoffEnergy+"");
			properties.setProperty(prefix+"psfEllipseScale",this.psfEllipseScale+"");
			properties.setProperty(prefix+"rpsfMinMaskThreshold",this.rpsfMinMaskThreshold+"");
			properties.setProperty(prefix+"filter",this.filter+"");
			properties.setProperty(prefix+"blurIndividual",this.blurIndividual+"");
			properties.setProperty(prefix+"blurDiagonal",this.blurDiagonal+"");
			properties.setProperty(prefix+"blurChecker",this.blurChecker+"");
			properties.setProperty(prefix+"gaussianSigmaIndividual",this.gaussianSigmaIndividual+"");
			properties.setProperty(prefix+"gaussianSigmaDiagonal",this.gaussianSigmaDiagonal+"");
			properties.setProperty(prefix+"gaussianSigmaChecker",this.gaussianSigmaChecker+"");
			properties.setProperty(prefix+"sigmaScale",this.sigmaScale+"");
			properties.setProperty(prefix+"sigmaToRadius",this.sigmaToRadius+"");
			properties.setProperty(prefix+"filterDirect",this.filterDirect+"");
			properties.setProperty(prefix+"sigmaScaleDirect",this.sigmaScaleDirect+"");
			properties.setProperty(prefix+"sigmaToRadiusDirect",this.sigmaToRadiusDirect+"");
		}

		public void getProperties(String prefix,Properties properties){
			this.dSize=Integer.parseInt(properties.getProperty(prefix+"dSize"));
			this.rSize=Integer.parseInt(properties.getProperty(prefix+"rSize"));
			this.invertRange=Double.parseDouble(properties.getProperty(prefix+"invertRange"));
			this.otfCutoffEnergy=Double.parseDouble(properties.getProperty(prefix+"otfCutoffEnergy"));
			this.otfEllipseScale=Double.parseDouble(properties.getProperty(prefix+"otfEllipseScale"));
			this.otfEllipseGauss=Boolean.parseBoolean(properties.getProperty(prefix+"otfEllipseGauss"));
			this.psfCutoffEnergy=Double.parseDouble(properties.getProperty(prefix+"psfCutoffEnergy"));
			this.psfEllipseScale=Double.parseDouble(properties.getProperty(prefix+"psfEllipseScale"));
			this.rpsfMinMaskThreshold=Double.parseDouble(properties.getProperty(prefix+"rpsfMinMaskThreshold"));
			this.filter=Boolean.parseBoolean(properties.getProperty(prefix+"filter"));
			this.blurIndividual=Double.parseDouble(properties.getProperty(prefix+"blurIndividual"));
			this.blurDiagonal=Double.parseDouble(properties.getProperty(prefix+"blurDiagonal"));
			this.blurChecker=Double.parseDouble(properties.getProperty(prefix+"blurChecker"));
			this.gaussianSigmaIndividual=Double.parseDouble(properties.getProperty(prefix+"gaussianSigmaIndividual"));
			this.gaussianSigmaDiagonal=Double.parseDouble(properties.getProperty(prefix+"gaussianSigmaDiagonal"));
			this.gaussianSigmaChecker=Double.parseDouble(properties.getProperty(prefix+"gaussianSigmaChecker"));
			this.sigmaScale=Double.parseDouble(properties.getProperty(prefix+"sigmaScale"));
			this.sigmaToRadius=Double.parseDouble(properties.getProperty(prefix+"sigmaToRadius"));
			this.filterDirect=Boolean.parseBoolean(properties.getProperty(prefix+"filterDirect"));
			this.sigmaScaleDirect=Double.parseDouble(properties.getProperty(prefix+"sigmaScaleDirect"));
			this.sigmaToRadiusDirect=Double.parseDouble(properties.getProperty(prefix+"sigmaToRadiusDirect"));
		}
	}
	
	
	public static class InterpolateParameters {
		public int    size;        // size of each kernel (should be square)
		public int    step;        // number of subdivisions from input to output
		public int    add_top;     // add this number of kernel rows to the output above the existent/interpolated
		public int    add_left;    // add this number of kernel columns to the output on the left of the existent/interpolated
		public int    add_right;   // add this number of kernel columns to the output on the right of the existent/interpolated
		public int    add_bottom;  // add this number of kernel rows to the output below the existent/interpolated
		public double extrapolate; // 0 - duplicate, 1.0 - extrapolate outside of the known kernels

		public InterpolateParameters(
				int    size,
				int    step,
				int    add_top,
				int    add_left,
				int    add_right,
				int    add_bottom,
				double extrapolate
		) {
			this.size=size;
			this.step=step;
			this.add_top=add_top;
			this.add_left=add_left;
			this.add_right=add_right;
			this.add_bottom=add_bottom;
			this.extrapolate=extrapolate;
		}
		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"size",this.size+"");
			properties.setProperty(prefix+"step",this.step+"");
			properties.setProperty(prefix+"add_top",this.add_top+"");
			properties.setProperty(prefix+"add_left",this.add_left+"");
			properties.setProperty(prefix+"add_right",this.add_right+"");
			properties.setProperty(prefix+"add_bottom",this.add_bottom+"");
			properties.setProperty(prefix+"extrapolate",this.extrapolate+"");
		}
		public void setProperties(String prefix,ImagePlus properties){
			properties.setProperty(prefix+"size",this.size+"");
			properties.setProperty(prefix+"step",this.step+"");
			properties.setProperty(prefix+"add_top",this.add_top+"");
			properties.setProperty(prefix+"add_left",this.add_left+"");
			properties.setProperty(prefix+"add_right",this.add_right+"");
			properties.setProperty(prefix+"add_bottom",this.add_bottom+"");
			properties.setProperty(prefix+"extrapolate",this.extrapolate+"");
		}
		
		public void getProperties(String prefix,Properties properties){
			this.size=Integer.parseInt(properties.getProperty(prefix+"size"));
			this.step=Integer.parseInt(properties.getProperty(prefix+"step"));
			this.add_top=Integer.parseInt(properties.getProperty(prefix+"add_top"));
			this.add_left=Integer.parseInt(properties.getProperty(prefix+"add_left"));
			this.add_right=Integer.parseInt(properties.getProperty(prefix+"add_right"));
			this.add_bottom=Integer.parseInt(properties.getProperty(prefix+"add_bottom"));
			this.extrapolate=Double.parseDouble(properties.getProperty(prefix+"extrapolate"));
		}
		
	}

	
	

}
