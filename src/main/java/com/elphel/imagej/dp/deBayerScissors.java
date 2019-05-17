package com.elphel.imagej.dp;
/**
 ** -----------------------------------------------------------------------------**
 ** deBayerScissors.java
 **
 ** Frequency-domain de-mosoaic filters generation
 ** 
 **
 ** Copyright (C) 2010 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  focus_tuning.java is free software: you can redistribute it and/or modify
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

import ij.process.*;
import ij.plugin.filter.GaussianBlur;
import java.util.HashSet;

import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
public class deBayerScissors {
	private PolarSpectrums        pol_instace=null;
	private double [][][]         lopass=null;
	private int                   size;
	private double                lastMidEnergy; // last midrange spectral energy
	private ShowDoubleFloatArrays SDFA_instance; // just for debugging?
	private DoubleFHT             fht_instance;
	private int [][] aliasMapRedBlue={
			{-2,-2},{-2,-1},{-2,0},{-2,1},
			{-1,-2},{-1,-1},{-1,0},{-1,1},
			{ 0,-2},{ 0,-1},       { 0,1},
			{ 1,-2},{ 1,-1},{ 1,0},{ 1,1}};
    private int [][][][] speedTable = null;

	public double getMidEnergy() {return lastMidEnergy; } // instead of the  DOUBLE_DEBUG_RESULT
	public deBayerScissors(
			int isize, // size of the square array, centar is at size/2, size/2, only top half+line will be used
			double polarStep, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
			double debayer_width_green, // result green mask mpy by scaled default (diamond)
			double debayer_width_redblue, // result red/blue mask mpy by scaled default (square)
			double debayer_width_redblue_main, // green mask when applied to red/blue, main (center)
			double debayer_width_redblue_clones){// green mask when applied to red/blue, clones 
		size=isize;
		fht_instance=       new DoubleFHT();
		SDFA_instance=      new ShowDoubleFloatArrays();
		pol_instace=new PolarSpectrums(size, // size of the square array, centar is at size/2, size/2, only top half+line will be used
				Math.PI,
				size/2-2, // width of the polar array - should be <= size/2-2
				polarStep, //0.5, //0.75, //2.0, //0.5, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
				4);// angular symmetry - 0- none,1 - pi corresponds to integer, 2 - pi/2 corresponds to integer, n - pi/n corresponds to integer angular step 

		lopass=  createAliasFilters (debayer_width_green, // result green mask mpy by scaled default (diamond)
				debayer_width_redblue, // result red/blue mask mpy by scaled default (square)
				debayer_width_redblue_main, // green mask when applied to red/blue, main (center)
				debayer_width_redblue_clones, // green mask when applied to red/blue, clones 
				size, // side of the square
				4); // should be 4 now

		int hsize=size/2;
		int subpixel=4; // hardwired - when changing it will need to change alias maps
		int aliasX=size/subpixel;
	    speedTable = new int [hsize+1][size][aliasMapRedBlue.length][3];
        int i,j,nAlias,x,y;
		for (i = 0;i <= hsize;i++) for (j = 0;j < size;j++) {
			speedTable[i][j][0][0] = i*size+j;
			speedTable[i][j][0][1] = ((size-i) % size) * size + ((size-j) % size);
			for(nAlias=0;nAlias<aliasMapRedBlue.length; nAlias++) {
				y = (i-aliasX*aliasMapRedBlue[nAlias][0]+size) % size;
				x = (j-aliasX*aliasMapRedBlue[nAlias][1]+size) % size;
				if (y > hsize) {
					y = size - y;
					x = (size - x) % size;
				}
				if (y>hsize) {
					y=size-y;
					x=(size-x)%size;
				}
				speedTable[i][j][nAlias][2]=y*size+x;
			}
		}
	}

	/* returns 2 masks (0:0 in the top left corner, match fht) [0] - for greens, [1] - for red/blue */
	/* Possible improvements: - 1 make the initial green mask (or actually "fan"-like image) to have sharper ends.
                             2. detect periodic (line of spots) on the spectrum amplitudes (corresponds to thin lines) and use this
                                info to confirm this area to belong to the main spectrum */

	public double [][] aliasScissors(double [] green_fht, // fht array for green, will be masked in-place
			double debayer_threshold, // no high frequencies - use default uniform filter
			double debayer_gamma, // power function applied to the amplitudes before generating spectral masks
			double debayer_bonus, // scale far pixels as (1.0+bonus*r/rmax)
			double mainToAlias,// relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
			double debayer_mask_blur, // for both masks  sigma for Gaussian blur of the binary masks (<0 -do not use "scissors")
			boolean debayer_use_scissors, // use "scissors", if false - just apply "diamond" ands "square" with DEBAYER_WIDTH_GREEN and DEBAYER_WIDTH_REDBLUE
			int this_debug){ // internal debug level
		int length=green_fht.length;
		int size=(int) Math.sqrt(length);
		double [] green_mask;
		double [] red_blue_mask;
		double [] green_amp=fht_instance.calculateAmplitude(green_fht);
		int i,j;
		/**normalize amplitudes, apply gamma */
		double dmax=0.0;
		for (i=0;i<green_amp.length;i++) if (green_amp[i]>dmax) dmax=green_amp[i];
		dmax=1.0/dmax;
		for (i=0;i<green_amp.length;i++) green_amp[i]= Math.pow(green_amp[i]*dmax,debayer_gamma);
		if (this_debug>2)   SDFA_instance.showArrays(green_amp,  "DT-gam"); // only top half+1 will be used
		double midRangeSpectral=pol_instace.maxAmpInRing (green_amp);
		boolean useFancyDebayer=(midRangeSpectral>=debayer_threshold);

		lastMidEnergy= midRangeSpectral; // for optional monitoring outside of this class

		if (useFancyDebayer && debayer_use_scissors) { /* calculate and apply "scissors" masks */
			green_mask= calcGreensAliasMaskRays (green_amp, // normalized amplitude spectrum, (0,0) in the center
					pol_instace, // initialized instance
					debayer_bonus, // hack - here it is "bonus"
					this_debug);//

			if (this_debug>3) SDFA_instance.showArrays(green_mask,  "G-raw");
			if (debayer_mask_blur>0) {
				blurDouble(green_mask,   size, debayer_mask_blur, debayer_mask_blur, 0.01);
				if (this_debug>3) SDFA_instance.showArrays(green_mask,  "G-blurred");
			}
			double [] green_mask_for_redblue_main=  green_mask.clone();
			double [] green_mask_for_redblue_clones=green_mask.clone();
			for (i=0;i<green_mask.length;i++) {
				green_mask_for_redblue_main[i]*=  lopass[2][0][i];
				green_mask_for_redblue_clones[i]*=lopass[2][1][i];
			}

			if (this_debug>2) {
				SDFA_instance.showArrays(green_mask_for_redblue_main,  "MAIN");
				SDFA_instance.showArrays(green_mask_for_redblue_main,  "CLONES");
			}

			/* Maybe here we need to unmasked (wide bandwidth) green_amp? */
			red_blue_mask= calcRedBlueAliasMaskRays (
					green_amp, // both halves are needed ??
					green_mask_for_redblue_main, // may be null if amp_pixels is already masked
					green_mask_for_redblue_clones,
					pol_instace, // initialized instance (if null - skip rays processing)
					mainToAlias,// relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
					debayer_bonus, // scale far pixels as (1.0+bonus*r/rmax)
					this_debug);// relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

			/* add    double mainToAlias){// relative main/alias amplitudes to enable pixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out) */

			if (this_debug>3) SDFA_instance.showArrays(red_blue_mask,  "RB-raw");
			if (debayer_mask_blur>0) {
				blurDouble(red_blue_mask, size,debayer_mask_blur, debayer_mask_blur, 0.01);
				if (this_debug>3) SDFA_instance.showArrays(red_blue_mask,  "RB-blurred");
			}
			for (i=0;i<red_blue_mask.length;i++) red_blue_mask[i]*=lopass[1][1][i]; //  scaled, red-blue - was red_blue_lopass[i];
		} else { // debayer_mask_blur<0 : use default masks
			green_mask=lopass[1][0].clone(); //green_lopass.clone(); variable (wide) filter here)
			red_blue_mask=lopass[1][1].clone(); //red_blue_lopass.clone();
			if (!useFancyDebayer) for (i=0;i<green_mask.length;i++) { // no high-frequency componets detected - reduce noise by extra (narrow) filtering
				green_mask[i]*=   lopass[0][0][i]; // *=   green_lopass[i];
				red_blue_mask[i]*=lopass[0][1][i]; // *=red_blue_lopass[i];
			}
		}
		/* Swap quadrants in the masks to match FHT arrays (0:0 in the top left corner) */
		fht_instance.swapQuadrants(green_mask);
		fht_instance.swapQuadrants(red_blue_mask);
		/* return both masks */
		double [][] result =new double [2][];
		result[0]= green_mask;
		result[1]= red_blue_mask;
		//    if (this_debug>3) SDFA_instance.showArrays(result,  "before_norm_masks");


		/* normalize masks to have exactly 1.0 at 0:0 - it can be reduced by blurring */
		for (i=0;i<result.length;i++) {
			dmax=1.0/result[i][0];
			for (j=0;j<result[i].length;j++) result[i][j]*=dmax;
		}
		//    if (this_debug>3) SDFA_instance.showArrays(result,  "masks");
		return result;
	}

	public double [] calcRedBlueAliasMaskRays (
			double [] green_amp, // both halves are needed ??
			double [] green_mask, // may be null if amp_pixels is already masked
			double [] green_mask_clones, // mask (more inclusive than just green_mask) to be used with clones
			PolarSpectrums pol_instace, // initialized instance (if null - skip rays processing)
			double mainToAlias,// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
			double bonus, // scale far pixels as (1.0+bonus*r/rmax)
			int this_debug){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)

		int length=green_amp.length;
		int size = (int) Math.sqrt(length);
		int hsize=size/2;
//		int subpixel=4; // hardwired - when changing it will need to change alias maps
//		int aliasX=size/subpixel;
		int i,j,index,index_back; //,x,y;
		double [] amp=       green_amp.clone();
		double [] amp_clones=green_amp.clone();
		if (green_mask!=null)        for (i=0;i<amp.length;i++)        amp[i]*=green_mask[i];
		if (green_mask_clones!=null) for (i=0;i<amp_clones.length;i++) amp_clones[i]*=green_mask_clones[i];
		double [] mask= new double [length];
		for (i=0;i<length;i++) mask[i]=0.0;
		/* Combine into mask by comparing pixels[] from the zero and 7 aliases */
		double d;
		int nAlias;
		/*
		int [][] aliasMapRedBlue={
				{-2,-2},{-2,-1},{-2,0},{-2,1},
				{-1,-2},{-1,-1},{-1,0},{-1,1},
				{ 0,-2},{ 0,-1},       { 0,1},
				{ 1,-2},{ 1,-1},{ 1,0},{ 1,1}};

		// First step - mask out all the pixels where at least one of the alias amplitude is above the main one
		if (this_debug>2) SDFA_instance.showArrays(amp.clone(),  "amp");
		if (this_debug>2) SDFA_instance.showArrays(amp_clones,  "amp_clones");

		for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
			index=i*size+j;
			index_back=((size-i) % size) * size + ((size-j) % size);
			d=amp[index]*mainToAlias;
			if (d>0.0) {
				mask[index]=1.0;
				mask[index_back]=1.0;
				//        isGreater=true;
				for(nAlias=0;nAlias<aliasMapRedBlue.length; nAlias++) {
					y=(i-aliasX*aliasMapRedBlue[nAlias][0]+size) % size;
					x=(j-aliasX*aliasMapRedBlue[nAlias][1]+size) % size;
					if (y>hsize) {
						y=size-y;
						x=(size-x)%size;
					}
					if (amp_clones[y*size+x]>d) {
						mask[index]=-1.0;
						mask[index_back]=-1.0;
						break;
					}
					
				}
			}
		}
		*/
		/*
	    speedTable = new int [hsize][size][aliasMapRedBlue.length][3];
        int i,j,nAlias,x,y;
		for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
			speedTable[i][j][0][0] = i*size+j;
			speedTable[i][j][0][1] = ((size-i) % size) * size + ((size-j) % size);
			for(nAlias=0;nAlias<aliasMapRedBlue.length; nAlias++) {
				y = (i-aliasX*aliasMapRedBlue[nAlias][0]+size) % size;
				x = (j-aliasX*aliasMapRedBlue[nAlias][1]+size) % size;
				if (y > hsize) {
					y = size - y;
					x = (size - x) % size;
				}
				if (y>hsize) {
					y=size-y;
					x=(size-x)%size;
				}
				speedTable[i][j][nAlias][2]=y*size+x;
			}
		}
		*/

		// First step - mask out all the pixels where at least one of the alias amplitude is above the main one
		if (this_debug>2) SDFA_instance.showArrays(amp.clone(),  "amp");
		if (this_debug>2) SDFA_instance.showArrays(amp_clones,  "amp_clones");

		for (i=0;i<=hsize;i++) for (j=0;j<size;j++) {
			index =      speedTable[i][j][0][0]; // i*size+j;
			index_back = speedTable[i][j][0][1]; //((size-i) % size) * size + ((size-j) % size);
			d=amp[index]*mainToAlias;
			if (d>0.0) {
				mask[index]=1.0;
				mask[index_back]=1.0;
				for(nAlias=0;nAlias<aliasMapRedBlue.length; nAlias++) {
					if (amp_clones[speedTable[i][j][nAlias][2]]>d) {
						mask[index]=-1.0;
						mask[index_back]=-1.0;
						break;
					}
				}
			}
		}
		
// End of replacement code
		
		if (this_debug>2)  SDFA_instance.showArrays(mask,  "mask");

		if (pol_instace==null) return mask;
		/* Now apply mask to amplitudes and use ray processing (same as with greens)*/
		for (i=0;i<amp.length;i++) amp[i]*=mask[i];
		if (this_debug>2) SDFA_instance.showArrays(amp,  "amp-mask");
		double [] polar_amp=pol_instace.cartesianToPolar(amp);
		int width = pol_instace.getWidth();
		int height = pol_instace.getHeight();
		if (this_debug>2)   SDFA_instance.showArrays(polar_amp.clone(),width, height,  "RB-polar-amp");
		double k= bonus/width;
		for (i=0;i<pol_instace.getHeight();i++) for (j = 0; j < width; j++) polar_amp[i * width + j]*=1.0+k*j;
		double [] polar_mask_pixels=pol_instace.genPolarRedBlueMask(polar_amp,0); // 0 - just 1.0/0.0, 1 - "analog"
		double [] cart_mask_pixels= pol_instace.polarToCartesian (polar_mask_pixels,size,0.0);
		if (this_debug>2) {
			SDFA_instance.showArrays(polar_amp,         width, height,     "RB-amp-bonus");
			SDFA_instance.showArrays(polar_mask_pixels, width, height, "pRBm");
			SDFA_instance.showArrays(cart_mask_pixels,  size,  size,   "cRBm");
		}
		if (this_debug>2) {
			double [] polar_mask_pixels1=pol_instace.genPolarRedBlueMask(polar_amp,1);
			double [] cart_mask_pixels1= pol_instace.polarToCartesian (polar_mask_pixels1,size,0.0);
			SDFA_instance.showArrays(polar_mask_pixels1, width, height, "pRBm1");
			SDFA_instance.showArrays(cart_mask_pixels1,  size,  size,   "cRBm1");
		}
		return cart_mask_pixels;

	}


	public double [] calcGreensAliasMaskRays (double [] amp_pixels, // normalized amplitude spectrum, (0,0) in the center
			PolarSpectrums pol_instace, // initialized instance
			double bonus, // scale far pixels as (1.0+bonus*r/rmax)
			int this_debug){// relative main/alias amplitudes to enable lixels (i.e. 0.5 means that if alias is >0.5*main, the pixel will be masked out)
		int length=amp_pixels.length;
		int size = (int) Math.sqrt(length);
		double [] polar_amp_pixels=pol_instace.cartesianToPolar(amp_pixels);
		if (this_debug>2)   SDFA_instance.showArrays(polar_amp_pixels.clone(),pol_instace.getWidth(),pol_instace.getHeight(),  "polar-amp");

		double k= bonus/pol_instace.getWidth();
		for (int i=0;i<pol_instace.getHeight();i++) for (int j=0;j<pol_instace.getWidth();j++) polar_amp_pixels[i*pol_instace.getWidth()+j]*=1.0+k*j;
		double [] polar_green_mask_pixels=pol_instace.genPolarGreenMask(polar_amp_pixels,0); // 0 - just 1.0/0.0, 1 - "analog"
		double [] cart_green_mask_pixels= pol_instace.polarToCartesian (polar_green_mask_pixels,size,0.0);
		if (this_debug>2) {
			SDFA_instance.showArrays(polar_amp_pixels,  pol_instace.getWidth(),pol_instace.getHeight(),     "amp-bonus");
			SDFA_instance.showArrays(polar_green_mask_pixels,pol_instace.getWidth(),pol_instace.getHeight(), "pgm");
			SDFA_instance.showArrays(cart_green_mask_pixels,size,size,   "cgm");
		}

		if (this_debug>2) {
			double [] polar_green_mask_pixels1=pol_instace.genPolarGreenMask(polar_amp_pixels,1);
			double [] cart_green_mask_pixels1= pol_instace.polarToCartesian (polar_green_mask_pixels1,size,0.0);
			SDFA_instance.showArrays(polar_green_mask_pixels1,pol_instace.getWidth(),pol_instace.getHeight(), "PGM1");
			SDFA_instance.showArrays(cart_green_mask_pixels1,size,size,   "CGM1");
		}
		return cart_green_mask_pixels;

	}


	double [][][] createAliasFilters (double debayer_width_green, // result green mask mpy by scaled default (diamond)
			double debayer_width_redblue, // result red/blue mask mpy by scaled default (square)
			double debayer_width_redblue_main, // green mask when applied to red/blue, main (center), square
			double debayer_width_redblue_clones, // green mask when applied to red/blue, clones , square
			int size, // side of the square
			int subpixel){ // should be 4 now
		int i;
		double [] cosMask= createCosMask (size,  subpixel);   // oversampling
		double [][] [] lopass =new double [3][2][];
		lopass[0][0]=new double [size*size];
		for (i=0;i<lopass[0][0].length;i++) lopass[0][0][i]=1.0;
		lopass[0][1]=lopass[0][0].clone();
		lopass[1][0]=lopass[0][0].clone();
		lopass[1][1]=lopass[0][0].clone();
		lopass[2][0]=lopass[0][0].clone();
		lopass[2][1]=lopass[0][0].clone();
		maskBayerAliases (lopass[0][0],   // FHT array to be filtered
				cosMask,   // cosine mask array
				true); // this fht array is for the checkerboard greens
		maskBayerAliases (lopass[0][1],   // FHT array to be filtered
				cosMask,   // cosine mask array
				false); // this fht array is for the checkerboard greens

		cosMask= createCosMask ((int) Math.round(size*debayer_width_green),  subpixel);   // oversampling
		maskBayerAliases (lopass[1][0],   // FHT array to be filtered
				cosMask,   // cosine mask array
				true); // this fht array is for the checkerboard greens
		cosMask= createCosMask ((int) Math.round(size*debayer_width_redblue),  subpixel);   // oversampling
		maskBayerAliases (lopass[1][1],   // FHT array to be filtered
				cosMask,   // cosine mask array
				false); // this fht array is for the checkerboard greens

		cosMask= createCosMask ((int) Math.round(size*debayer_width_redblue_main),  subpixel);   // oversampling
		maskBayerAliases (lopass[2][0],   // FHT array to be filtered
				cosMask,   // cosine mask array
				false); // this fht array is for the checkerboard greens
		cosMask= createCosMask ((int) Math.round(size*debayer_width_redblue_clones),  subpixel);   // oversampling
		maskBayerAliases (lopass[2][1],   // FHT array to be filtered
				cosMask,   // cosine mask array
				false); // this fht array is for the checkerboard greens


		fht_instance.swapQuadrants(lopass[0][0]);
		fht_instance.swapQuadrants(lopass[0][1]);
		fht_instance.swapQuadrants(lopass[1][0]);
		fht_instance.swapQuadrants(lopass[1][1]);
		fht_instance.swapQuadrants(lopass[2][0]);
		fht_instance.swapQuadrants(lopass[2][1]);
		return lopass;
	}

	void maskBayerAliases (double [] fht,   // FHT array to be filtered
			double [] cosMask,   // cosine mask array
			boolean isChecker) { // this fht array is for the checkerboard greens
		int size= (int) Math.sqrt(fht.length);
		int iy,ix, ix1,iy1;
		int tsize= (cosMask.length-1)/(isChecker?1:2);
		int index=0;
		int hsizeM1=(size/2)-1;
		if (isChecker) {

			for (iy=0;iy<size;iy++) {
				iy1=(iy+hsizeM1)%size -hsizeM1;
				for (ix=0;ix<size;ix++) {
					ix1=(ix+hsizeM1)%size -hsizeM1;
					if (((ix1+iy1)>-tsize) &&
							((ix1-iy1)>-tsize) &&
							((ix1+iy1)<=tsize) &&
							((ix1-iy1)<=tsize)) fht[index++]*=cosMask[Math.abs(ix1+iy1)]*cosMask[Math.abs(ix1-iy1)];
					else fht[index++]=0.0;
				}
			}

		} else {
			for (iy=0;iy<size;iy++) {
				iy1=(iy+hsizeM1)%size -hsizeM1;
				for (ix=0;ix<size;ix++) {
					ix1=(ix+hsizeM1)%size -hsizeM1;
					if ((iy1>-tsize) && (iy1<=tsize) && (ix1>-tsize) && (ix1<=tsize)) fht[index++]*=cosMask[2*Math.abs(iy1)]*cosMask[2*Math.abs(ix1)];
					else fht[index++]=0.0;
				}
			}
		}
	}


	double [] createCosMask (int fftsize,   // FHT array to be filtered - just length is used
			int subdiv   // oversampling
			) { // this fht array is for the checkerboard greens
		int size= 2*fftsize/subdiv;
		double [] cosMask=new double [size+1];
		for (int i=0;i<=size;i++) cosMask[i]=0.5*(1.0+Math.cos(i*Math.PI/size));
		return cosMask; 
	}

	// temporary using float implementation in ImageJ - re-write to directly use double [] arrays
	public void  blurDouble(double[] pixels,
			int width,
			double sigmaX,
			double sigmaY,
			double precision) {
		//  public void  blurFloat(red_blue_mask, DEBAYER_MASK_BLUR, DEBAYER_MASK_BLUR, 0.01);
		int i;
		int height = pixels.length/width;
		float [] fpixels=new float [pixels.length];
		for (i=0;i<pixels.length;i++) fpixels[i]= (float) pixels[i];
		FloatProcessor fp = new FloatProcessor(width, height, fpixels, null);
		GaussianBlur gb = new GaussianBlur();
		gb.blurFloat(fp, sigmaX, sigmaY, precision);
		for (i=0;i<pixels.length;i++) pixels[i]=fpixels[i];
	}
	/* ====================================================== */
	public class PolarSpectrums  {
		public int radius=0;
		public int iRadiusPlus1; // number of radius steps
		public int iAngle;
		public double aStep;
		public double rStep;
		public int size;
		public int length;
		// Make them private later, after debugging
		private int    [][] polar2CartesianIndices;   // for each polar angle/radius (angle*iRadiusPlus1+radius) - 4 interpolation corners (0:0, dx:0, 0:dy, dx:dy), the first (0:0) being the closest to the polar point
		private double [][] polar2CartesianFractions; //  a pair of dx, dy for interpolations (used with ) polar2CartesianIndices[][]]
		private int    [][]   cartesian2PolarIndices;   // each per-pixel array is a list of indices in polar array pointing to this cell (may be empty)
		private int    []     cartesian2PolarIndex;     // Cartesian->polar array index (cell closest to the center). Is it possible that cartesian2PolarIndices does not include this one?
		private int    [][]   polarGreenMap=null ;     // each element is a variable length integer array with a list of the alias indices
		private int    [][]   polarRedBlueMap=null ;    // each element is a variable length integer array with a list of the alias indices
		private int    [][]   sameCartesian=null ;     // each element is a variable length integer array with a list of indices of the other polar cells that belong (point to) the same cartesian cell
		private int    []     cartAmpList = null;      // list of indices of the elements of the cartesian array (symmetrical around the center) so the distance is between ampRMinMax[0] and ampRMinMax[1]
		private double []     ampRMinMax  ={0.0,0.0};
		public PolarSpectrums() { }  // so "Compile and Run" will be happy
		/* Convert cartesian to polar array, dimensions are set in the class constructor. Uses bi-linear interpolation */
		public double [] cartesianToPolar (double [] cartPixels ) {
			double [] polPixels=new double[iRadiusPlus1*(iAngle+1)];
			int i;
			for (i=0;i<polPixels.length;i++) {
				polPixels[i]=(1-polar2CartesianFractions[i][1])*( (1-polar2CartesianFractions[i][0])*cartPixels[polar2CartesianIndices[i][0]]  +
						polar2CartesianFractions[i][0]*cartPixels[polar2CartesianIndices[i][1]])+
						polar2CartesianFractions[i][1] *( (1-polar2CartesianFractions[i][0])*cartPixels[polar2CartesianIndices[i][2]]  +
								polar2CartesianFractions[i][0]*cartPixels[polar2CartesianIndices[i][3]]) ;
			}
			return polPixels;
		}
		public double [] polarToCartesian (double [] polPixels , int height, double undefined ) { return polarToCartesian (polPixels , height, undefined, height==size); }
		public double [] polarToCartesian (double [] polPixels , double undefined)              { return polarToCartesian (polPixels ,size, undefined,false); }
		public double [] polarToCartesian (double [] polPixels , int height )                   { return polarToCartesian (polPixels , height, Double.NaN,height==size); }
		public double [] polarToCartesian (double [] polPixels )                                { return polarToCartesian (polPixels , size, Double.NaN,false);  }
		public double [] polarToCartesian (double [] polPixels,
				int height, // for partial arrays
				double undefined, // use this value in the undefined areas
				boolean   symmHalf){ // add center-symmetrical top to the bottom(spectrums of real signals)
			int length=size*height;
			double [] cartPixels=new double[length];
			int i,j;
			int [] sameCartCell;
			double d;
			int l=symmHalf?((size+1)*size/2+1) :(size*height);
			int l2=(size+1)*size;
			for (i=0;i<l;i++) {
				sameCartCell=cartesian2PolarIndices[i];
				if (sameCartCell==null) {
					if (cartesian2PolarIndex[i]>=0) cartPixels[i]=polPixels[cartesian2PolarIndex[i]];
					else cartPixels[i]=undefined;
				} else {
					d=0;
					for (j=0;j<sameCartCell.length;j++) d+=polPixels[sameCartCell[j]];
					cartPixels[i]=d/sameCartCell.length;
				}
				if (symmHalf) {
					j=l2-i;
					if (j<length) cartPixels[j] = cartPixels[i];
				}
			}
			return cartPixels;
		}
		/* Caculates maximal value of a center-symmetrical array of the amplitudes in a ring. Uses cached table of indices, recalculates if it changed */
		public double maxAmpInRing ( double []amps ){ return  maxAmpInRing (amps,size*0.118,size*0.236);} // ~=1/3* (Math.sqrt(2)/4), 2/3* (Math.sqrt(2)/4) (center 1/3 ring between center and the closest alias for greens)
		public double maxAmpInRing ( double []amps,
				double rMin,
				double rMax
				){
			int i,j,x,y;
			if ((cartAmpList==null) || (rMin!=ampRMinMax[0]) || (rMax!=ampRMinMax[1])) {
				ampRMinMax[0]=rMin;
				ampRMinMax[1]=rMax;
				double rMin2=rMin*rMin;
				double rMax2=rMax*rMax;
				double r2;
				// pass 1 - count number of elements
				int numMembers=0;
				for (i=0;i<=size/2;i++) {
					y=i-(size/2);
					for (j=0;j<size;j++) {
						x=j-(size/2);
						r2=x*x+y*y;
						if ((r2>=rMin2) && (r2<=rMax2)) numMembers++;
					}
				}
				cartAmpList=new int [numMembers];
				// pass 2 - count number of elements fill in the array
				numMembers=0;
				for (i=0;i<=size/2;i++) {
					y=i-(size/2);
					for (j=0;j<size;j++) {
						x=j-(size/2);
						r2=x*x+y*y;
						if ((r2>=rMin2) && (r2<=rMax2)) cartAmpList[numMembers++]=i*size+j;
					}
				}
			}
			if (cartAmpList.length<1) return Double.NaN;
			double max=amps[cartAmpList[0]];
			for (i=1;i<cartAmpList.length;i++) if (max<amps[cartAmpList[i]]) max=amps[cartAmpList[i]];
			return max;
		}
		/* return polar array width (== radius+1) */
		public int getWidth() { return iRadiusPlus1; } 
		public int getHeight() { return iAngle+1; } 
		public double [] genPolarGreenMask(double [] polarAmps, // polar array of amplitude values, <0 - stop
				int mode){ // result mode - 0: output mask as 0/1, 1 -output proportional, positive - pass, negative - rejected
			return genPolarMask(polarAmps,0,mode);
		}

		public double [] genPolarRedBlueMask(
				double [] polarAmps, // polar array of amplitude values, <0 - stop
				int mode)
		{ // result mode - 0: output mask as 0/1, 1 -output proportional, positive - pass, negative - rejected
			return genPolarMask(polarAmps, 1, mode);
		}

// **** Seems to be most time-critical ****
		public double [] genPolarMask(double [] polarAmps, // polar array of amplitude values, <0 - stop
				int type, // 0 - green, 1 red/blue
				int mode){ // result mode - 0: output mask as 0/1, 1 -output proportional, positive - pass, negative - rejected
			int [][] polarMap=(type>0)?polarRedBlueMap: polarGreenMap;
			int length=iRadiusPlus1*(iAngle+1);
			int [] intMap= new int[length];
			int i,ia;
			for (i=0;i<intMap.length;i++) intMap[i]=0;
			int []    rayLength=new int[iAngle+1]; // index (radius)) of the current ray end for this angle
			boolean [] rayOpen= new boolean[iAngle+1]; // this ray can grow (not blocked)
			for (i=0;i<iAngle;i++) {
				rayLength[i]=0;
				rayOpen[i]=true;
			}
			int lastIndex;
			int base;
			int l;
			double max;
			int newVal;
			int step=0;
			int iMax=0; // index of the best ray
			int index=0;
			boolean good=true;
			while (iMax>=0) {
				step++;
				/* add polar point index */
				newVal=good?step:-step;
				intMap[index]=newVal;
				if (sameCartesian[index]!=null)
					for (i=0;i<sameCartesian[index].length;i++)
						intMap[sameCartesian[index][i]]=newVal;
				/* add aliases of point index (as negative values) */
				if ((good) &&(polarMap[index]!=null))
					for (i=0;i<polarMap[index].length;i++)
						intMap[polarMap[index][i]]=-step;
				/* update ray lengths and status */
				max=-1.0;
				iMax=-1;
				for  (ia=0;ia<=iAngle;ia++) if (rayOpen[ia]) {
					base=ia*iRadiusPlus1+1;  // second for this angle
					l=base+rayLength[ia]; // first after the pointer
					lastIndex=base+iRadiusPlus1; // first in the next row
					while ((l<lastIndex) && (intMap[l]>0)) l++;
					rayLength[ia]=l-base; // last "good" ( >0 and in the same row)
					if ((l==lastIndex) || (intMap[l]<0) || (polarAmps[l]<0.0) )
						rayOpen[ia]=false;
					else if (polarAmps[l]>max) {
						max=polarAmps[l];
						iMax=ia;
					}
					
				}
				if (iMax>=0) {
					rayLength[iMax]++;
					index=iMax*iRadiusPlus1+rayLength[iMax];
					/* See if any of the aliases of the new point  hit the positive value, then this point is prohibited (good=false). Otherwise add it with good=true */
					good=true;
					if (polarMap[index]!=null) for (i=0;i<polarMap[index].length;i++) {
						if (intMap[polarMap[index][i]]>0) {
							good=false;
							break;
						}
					}
				}
				/* index is set if (iMax>=0) */
			}
			double [] result=new double [intMap.length];
			if (mode==0) {
				for (i=0;i<length;i++) result[i]=(intMap[i]>0)?1.0:0.0;
			} else {
				for (i=0;i<length;i++) result[i]=(intMap[i]>0)?(step-intMap[i]):-(step+intMap[i]);
			}
			return result;
		}
		
		public PolarSpectrums(
				int isize, // size of the square array, center is at size/2, size/2, only top half+line will be used
				double fullAngle, // i.e. Math.PI, 2*Math.PI
				int    maxRadius, // width of the polar array - should be <= size/2-2
				double outerStep, // maximal step in pixels on the maxRadius for 1 angular step (i.e. 0.5)
				int         symm // angular symmetry - 0- none,1 - pi corresponds to integer, 2 - pi/2 corresponds to integer, n - pi/n corresponds to integer angular step     
				) {
			size= isize;
			length=size*size;
			if (maxRadius>(size/2-2)) maxRadius=(size/2-2);
			radius=maxRadius;
			if (symm==0) aStep=fullAngle/Math.ceil(fullAngle*radius/outerStep);
			else        aStep=Math.PI/symm/Math.ceil(Math.PI*radius/outerStep/symm);
			iRadiusPlus1=(int) Math.ceil(radius/outerStep)+1;
			rStep=radius/(iRadiusPlus1-1.0);
			iAngle=(int) Math.round(fullAngle/aStep);
			polar2CartesianIndices=  new int    [(iAngle+1)*iRadiusPlus1][4]; // [0] - closest one
			polar2CartesianFractions=new double [(iAngle+1)*iRadiusPlus1][2];
			int ia,ir,y,x,i,j; //,PolarIndex;
			double a,r,cos,sin,dy,dx;
			cartesian2PolarIndex=  new int[length];
			cartesian2PolarIndices=new int[length][];
			@SuppressWarnings("unchecked")
			HashSet  <Integer> [] polarList=  (HashSet  <Integer> []) new HashSet[length];
			for (i=0;i<length;i++) {
				polarList[i]=new HashSet <Integer>(); // 16, 0.75
			}
			Integer PolarIndex,CartesianIndex;
			for (ia=0;ia<=iAngle;ia++) {
				a=ia*aStep;
				cos=Math.cos(a);
				sin=Math.sin(a);
				for (ir=0;ir<iRadiusPlus1;ir++) {
					PolarIndex=ia*iRadiusPlus1+ir;
					r=ir*rStep;
					dy=r*sin;
					y=(int) Math.round(dy);
					dy-=y;
					i=size/2-y;
					dx=r*cos;
					x=(int) Math.round(dx);
					dx-=x;
					j=x+size/2;
					CartesianIndex=i*size+j;
					polar2CartesianIndices[PolarIndex][0]=CartesianIndex;
					polarList[CartesianIndex].add(PolarIndex);
					if (dx<0) {
						polar2CartesianIndices[PolarIndex][1]=polar2CartesianIndices[PolarIndex][0]-1;
						dx=-dx;
					} else {
						polar2CartesianIndices[PolarIndex][1]=polar2CartesianIndices[PolarIndex][0]+1;
					}
					if (dy<0) {
						polar2CartesianIndices[PolarIndex][2]=polar2CartesianIndices[PolarIndex][0]+size;
						polar2CartesianIndices[PolarIndex][3]=polar2CartesianIndices[PolarIndex][1]+size;
						dy=-dy;
					} else {
						polar2CartesianIndices[PolarIndex][2]=polar2CartesianIndices[PolarIndex][0]-size;
						polar2CartesianIndices[PolarIndex][3]=polar2CartesianIndices[PolarIndex][1]-size;
					}
					polar2CartesianFractions[PolarIndex][0]=dx;
					polar2CartesianFractions[PolarIndex][1]=dy;
				}
			}
			for (i=0;i<length;i++) {
				y=size/2-(i/size);
				x=(i % size)- size/2;
				r=Math.sqrt(x*x+y*y);
				a=Math.atan2(y,x);
				if (a<0) a+=2*Math.PI;
				ir =(int) Math.round(r/rStep);
				ia= (int) Math.round(a/aStep);
				if ((ir>=0) && (ir<iRadiusPlus1) && (ia>=0) && (ia<=iAngle)) {
					cartesian2PolarIndex[i]=ia*iRadiusPlus1+ir;  
					if (polarList[i].size()==0) cartesian2PolarIndices[i]=null;
					else {
						cartesian2PolarIndices[i]=new int[polarList[i].size()];
						j=0;
						for (Integer val : polarList[i]) cartesian2PolarIndices[i][j++]=val;
					}
				} else {
					cartesian2PolarIndex[i]=-1; // invalid  
					cartesian2PolarIndices[i]=null;
				}
			}
			initSameCartesian();
			polarGreenMap=  new int [iRadiusPlus1*(iAngle+1)][];
			initAliasMaps(0);
			polarRedBlueMap=new int [iRadiusPlus1*(iAngle+1)][];
			initAliasMaps(1);
		}

		public double [] testMapsLengths(int mode) { // 0 - return lengths of "sameCartesian[]", 1 - same for polarGreenMap
			int [][] map= (mode==0)?sameCartesian:((mode==1)?polarGreenMap:polarRedBlueMap);
			double [] result = new double [map.length];
			for (int i=0; i<map.length;i++) {
				result[i]=(map[i]==null)?0.0:map[i].length;
			}
			return result;
		}

		public double [] testGreenMap(int ia, int ir) {
			double [] result = new double [polarGreenMap.length];
			int i;
			for ( i=0; i<result.length;i++) result[i]=0.0;
			int index=ia*iRadiusPlus1+ir;
			if (polarGreenMap[index]!=null){
				for (i=0;i<polarGreenMap[index].length;i++) result [polarGreenMap[index][i]]+=1.0;
				System.out.println("testGreenMap("+ia+","+ir+"): polarGreenMap["+index+"].length="+polarGreenMap[index].length);
			} else {
				System.out.println("testGreenMap("+ia+","+ir+"): polarGreenMap["+index+"]=null");
			}
			result [index]=-1.0;
			return result;
		}

		public double [] testRedBlueMap(int ia, int ir) {
			double [] result = new double [polarRedBlueMap.length];
			int i;
			for ( i=0; i<result.length;i++) result[i]=0.0;
			int index=ia*iRadiusPlus1+ir;
			if (polarRedBlueMap[index]!=null){
				for (i=0;i<polarRedBlueMap[index].length;i++) result [polarRedBlueMap[index][i]]+=1.0;
				System.out.println("testRedBlueMap("+ia+","+ir+"): polarRedBlueMap["+index+"].length="+polarRedBlueMap[index].length);
			} else {
				System.out.println("testRedBlueMap("+ia+","+ir+"): polarRedBlueMap["+index+"]=null");
			}
			result [index]=-1.0;
			return result;
		}


		/* Create per-polar pixel list of aliases for green Bayer. For each polar point it shows the polar coordinates of the same (and rotated by pi) point of aliases */
		/* current implementation - us cartesian (original) pixels as all/nothing, maybe it makes sense to go directly polar-polar, but then it may leave gaps */
		public void initAliasMaps (int type) { // 0 - green, 1 - Red/Blue
			int [][] aliasMapGreen=  {{-2,-2},{-2,0},            // using rollover, so only unique aliases are needed
					{-1,-1},{-1,1},
					{ 0,-2},
					{ 1,-1},{ 1,1}};
			int [][] aliasMapRedBlue={{-2,-2},{-2,-1},{-2,0},{-2,1},
					{-1,-2},{-1,-1},{-1,0},{-1,1},
					{ 0,-2},{ 0,-1},       { 0,1},
					{ 1,-2},{ 1,-1},{ 1,0},{ 1,1}};
			int [][] aliasMap=(type>0)?aliasMapRedBlue:aliasMapGreen;
			int [][] polarMap=(type>0)?polarRedBlueMap: polarGreenMap;
			HashSet  <Integer> aliasList=new HashSet <Integer>();
			int j,ix,iy, nAlias, dirAlias,ixa,iya, index, polarIndex;
			for (polarIndex=0;polarIndex<polarMap.length;polarIndex++) {
				iy= size/2- (polar2CartesianIndices[polarIndex][0] / size);
				ix= (polar2CartesianIndices[polarIndex][0] % size)-size/2 ;
				aliasList.clear();
				for (nAlias=0;nAlias<aliasMap.length;nAlias++) for (dirAlias=-1;dirAlias<2;dirAlias+=2) {
					ixa=(size+ size/2+ aliasMap[nAlias][0]*size/4+ dirAlias*ix) % size;
					iya=(size+ size/2- aliasMap[nAlias][1]*size/4- dirAlias*iy) % size;
					index=iya*size + ixa;
					if (cartesian2PolarIndices[index]==null) {
						if (cartesian2PolarIndex[index]>=0) {
							aliasList.add (new Integer(cartesian2PolarIndex[index]));
						}
					} else {
						for (j=0;j<cartesian2PolarIndices[index].length;j++) {
							aliasList.add (new Integer(cartesian2PolarIndices[index][j]));
						}
					}
				}
				/*  convert set to int[] */
				if (aliasList.size()==0) polarMap[polarIndex]=null;
				else {
					polarMap[polarIndex]=new int[aliasList.size()];
					j=0;
					for (Integer val : aliasList) polarMap[polarIndex][j++]=val;
				}
			}
		}

		public void initSameCartesian () {
			int i,j, polarIndex, cartesianIndex;
			sameCartesian=new int [iRadiusPlus1*(iAngle+1)][];
			for (polarIndex=0;polarIndex<sameCartesian.length;polarIndex++) {
				cartesianIndex=polar2CartesianIndices[polarIndex][0];
				if ((cartesian2PolarIndices[cartesianIndex]==null) || (cartesian2PolarIndices[cartesianIndex].length<=1)) sameCartesian[polarIndex]=null;
				else {
					sameCartesian[polarIndex]=new int [cartesian2PolarIndices[cartesianIndex].length-1];
					j=0;
					/* copy all elements but this one - out of bounds may mean that it was not included - bug */
					for (i=0;i<cartesian2PolarIndices[cartesianIndex].length;i++) if (cartesian2PolarIndices[cartesianIndex][i]!=polarIndex) sameCartesian[polarIndex][j++]=cartesian2PolarIndices[cartesianIndex][i];
				}
			}
		}
	}

}

