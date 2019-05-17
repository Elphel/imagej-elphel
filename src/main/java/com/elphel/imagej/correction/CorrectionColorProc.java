package com.elphel.imagej.correction;
/**
** -----------------------------------------------------------------------------**
** CorrectionColorProc.java
**
** Color conversion methods used in aberration correction for Eyesis4pi
**
**
** Copyright (C) 2012 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  CorrectionColorProc.java is free software: you can redistribute it and/or modify
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

import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.GenericJTabbedDialog;
import com.elphel.imagej.common.ShowDoubleFloatArrays;
import com.elphel.imagej.common.WindowTools;

import ij.ImageStack;
import ij.gui.GenericDialog;


public class CorrectionColorProc {
	ShowDoubleFloatArrays SDFA_INSTANCE=   new ShowDoubleFloatArrays();
    public AtomicInteger stopRequested=null; // 1 - stop now, 2 - when convenient
    double [] denoiseMaskChroma;
    int       denoiseMaskChromaWidth;

    int debugLevel=1;
    public CorrectionColorProc(AtomicInteger stopRequested){
    	this.stopRequested=stopRequested;
    }
    double [] getDenoiseMaskChroma() {return this.denoiseMaskChroma;}
    int       getDenoiseMaskChromaWidth() {return this.denoiseMaskChromaWidth;}
    void setDebug(int debugLevel){this.debugLevel=debugLevel;}


    public void  processColorsWeights(ImageStack stack, // does not like NaN (in GaussianBlur
  		  double scale,     // initial maximal pixel value (16))
  		  EyesisCorrectionParameters.ColorProcParameters  colorProcParameters,
  		  CorrectionColorProc.ColorGainsParameters channelGainParameters,
  		  int channel,
  		  double [] denoiseMask,
  		  boolean blueProc,
  		  int debugLevel
    		) {
    	double thisGain=       colorProcParameters.gain;
    	double thisBalanceRed= colorProcParameters.balanceRed;
    	double thisBalanceBlue=colorProcParameters.balanceBlue;
    	if (channelGainParameters!=null) {
    		thisGain*=       channelGainParameters.gain[channel];
    		thisBalanceRed*= channelGainParameters.balanceRed[channel];
    		thisBalanceBlue*=channelGainParameters.balanceBlue[channel];
    	}
    	float [] fpixels_r= (float[]) stack.getPixels(1);
    	float [] fpixels_g= (float[]) stack.getPixels(2);
    	float [] fpixels_b= (float[]) stack.getPixels(3);
    	boolean useWeights=(stack.getSize()>=5);
    	if (!useWeights) {
    		stack.addSlice("dummy1",  fpixels_r);
    		stack.addSlice("dummy2",  fpixels_g);
    	}
    	float [] fpixels_wr=(float[]) stack.getPixels(4);
    	float [] fpixels_wb=(float[]) stack.getPixels(5);
    	int length=fpixels_r.length;
    	int width= stack.getWidth();
    	int height=stack.getHeight();
    	/* Scale colors, gamma-convert */
    	int i;
    	double gain_red= thisBalanceRed* thisGain/scale;
    	double gain_blue=thisBalanceBlue*thisGain/scale;
    	double gain_green=thisGain/scale;
    	double gamma_a=Math.pow(colorProcParameters.minLin,colorProcParameters.gamma)*(1.0-colorProcParameters.gamma);
    	gamma_a=gamma_a/(1.0-gamma_a);
    	double gamma_linK=(1.0+gamma_a)*colorProcParameters.gamma*Math.pow(colorProcParameters.minLin,colorProcParameters.gamma)/colorProcParameters.minLin;


    	double Kg=1.0-colorProcParameters.kr-colorProcParameters.kb;

    	// Correct color saturation for gamma
    	double Ar=colorProcParameters.kr*gain_red;
    	double Ag=Kg*gain_green;
    	double Ab=colorProcParameters.kb*gain_blue;
    	if (debugLevel>1) {
    		System.out.println ( " processColorsWeights() Ar="+Ar+" Ag="+Ag+" Ab="+Ab);
    		System.out.println ( " processColorsWeights() colorProcParameters.saturationBlue="+colorProcParameters.saturationBlue+
    				" colorProcParameters.saturationRed="+colorProcParameters.saturationRed);
    	}


    	//public void showArrays(double[][] pixels, int width, int height,  boolean asStack, String title) {




    	for (i=0;i<length;i++) {
    		double Y=Ar*fpixels_r[i]+Ag*fpixels_g[i]+Ab*fpixels_b[i];
    		Y= (Y !=0.0) ? (linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, Y)/Y): 0.0;
    		fpixels_r[i]*=Y*gain_red;
    		fpixels_g[i]*=Y*gain_green;
    		fpixels_b[i]*=Y*gain_blue;
    	}

    	if (colorProcParameters.corrBlueLeak && blueProc) {
    		double [][] blueLeakRgb =new double[3][length];
    		for (int px=0;px<length;px++){
    			blueLeakRgb[0][px]=fpixels_r[px];
    			blueLeakRgb[1][px]=fpixels_g[px];
    			blueLeakRgb[2][px]=fpixels_b[px];
    		}
    		BlueLeak blueLeak = new BlueLeak(
    				colorProcParameters,
    				blueLeakRgb,
    				width,
    				SDFA_INSTANCE,
    				null, // "blue_corr", //TODO: add normal generation/saving of the intermediate images
    				debugLevel); // debug level
    		double [][] blueRemovedRGB= blueLeak.process(); // will later return corrected RGB to use

    		for (int px=0;px<length;px++){
    			fpixels_r[px]=(float) blueRemovedRGB[0][px];
    			fpixels_g[px]=(float) blueRemovedRGB[1][px];
    			fpixels_b[px]=(float) blueRemovedRGB[2][px];
    		}
    	}

    	/* Convert to YPbPr */
    	double Y,Pb,Pr;
    	//        double Kg=1.0-colorProcParameters.kr-colorProcParameters.kb;
    	double Sb=0.5/(1.0-colorProcParameters.kb)*colorProcParameters.saturationBlue;
    	double Sr=0.5/(1.0-colorProcParameters.kr)*colorProcParameters.saturationRed;
    	double Yr,Yg,Yb,Wr,Wg,Wb,S;
    	/* coefficients to find Y from Pb, Pr and a color (R,G or B)
   Yr = R- Pr*KPrR
   Yb = B- Pb*KPbB
   Yg = G+ Pr*KPrG  + Pb*KPbG
    	 */
    	double KPrR= -(2.0*(1-colorProcParameters.kr))/colorProcParameters.saturationRed;
    	double KPbB= -(2.0*(1-colorProcParameters.kb))/colorProcParameters.saturationBlue;
    	double KPrG=  2.0*colorProcParameters.kr*(1-colorProcParameters.kr)/Kg/colorProcParameters.saturationRed;
    	double KPbG=  2.0*colorProcParameters.kb*(1-colorProcParameters.kb)/Kg/colorProcParameters.saturationBlue;
    	if (debugLevel> 1) {
    		System.out.println ( " processColorsWeights() gain_red="+gain_red+" gain_green="+gain_green+" gain_blue="+gain_blue);
    		System.out.println ( " processColorsWeights() gamma="+colorProcParameters.gamma+      " minLin="+colorProcParameters.minLin+" gamma_a="+gamma_a+" gamma_linK="+gamma_linK);
    		System.out.println ( " processColorsWeights() Kr="+colorProcParameters.kr+" Kg="+Kg+" Kb="+colorProcParameters.kb+" Sr="+Sr+" Sb="+Sb);
    		System.out.println ( " processColorsWeights() KPrR="+KPrR+" KPbB="+KPbB+" KPrG="+KPrG+" KPbG="+KPbG);
    	}

    	float [] fpixels_pb= new float [length];
    	float [] fpixels_pr= new float [length];
    	float [] fpixels_y0= new float [length];
    	float [] fpixels_y= fpixels_y0;

    	float [] fpixels_yR=null;
    	float [] fpixels_yG=null;
    	float [] fpixels_yB=null;

    	if (debugLevel>2) {
    		fpixels_yR= new float [length];
    		fpixels_yG= new float [length];
    		fpixels_yB= new float [length];
    	}
    	for (i=0;i<length;i++) {
    		Y=colorProcParameters.kr*fpixels_r[i]+Kg*fpixels_g[i]+colorProcParameters.kb*fpixels_b[i];
    		fpixels_pb[i] = (float) (Sb*(fpixels_b[i]-Y));
    		fpixels_pr[i] = (float) (Sr*(fpixels_r[i]-Y));
    		fpixels_y0[i]=(float) Y;
    	}
    	/* calculate Y from weighted colors, weights derived from how good each color component predicts signal in each subpixel of Bayer pattern */
    	if (useWeights) {
    		fpixels_y=  new float [length];
    		for (i=0;i<length;i++) {
    			Pb=fpixels_pb[i];
    			Pr=fpixels_pr[i];
    			Yr = fpixels_r[i]- Pr*KPrR;
    			Yb = fpixels_b[i]- Pb*KPbB;
    			Yg = fpixels_g[i]+ Pr*KPrG  + Pb*KPbG;
    			Wr=fpixels_wr[i];
    			Wb=fpixels_wb[i];
    			Wg=1.0-Wr-Wb;
    			S=1.0/(Wr*(colorProcParameters.weightScaleR-1.0)+Wb*(colorProcParameters.weightScaleB-1.0)+1.0);
    			Wr*=S*colorProcParameters.weightScaleR;
    			Wb*=S*colorProcParameters.weightScaleB;
    			Wg*=S;
    			Y=Yr*Wr+Yg*Wg+Yb*Wb;
    			fpixels_y[i]=(float) Y;
    			if (debugLevel>2) {
    				fpixels_yR[i]= (float) Yr;
    				fpixels_yG[i]= (float) Yg;
    				fpixels_yB[i]= (float) Yb;
    			}
    		}
    	}
    	/* Low-pass filter Pb and Pr */
    	DoubleGaussianBlur gb=new DoubleGaussianBlur();
    	double [] dpixels_pr=new double[fpixels_pr.length];
    	double [] dpixels_pb=new double[fpixels_pb.length];
    	for (i=0;i<dpixels_pr.length;i++) {
    		dpixels_pr[i]=fpixels_pr[i];
    		dpixels_pb[i]=fpixels_pb[i];
    	}
    	if (colorProcParameters.maskMax>0.0) {
    		double [] dmask=new double[fpixels_y0.length];
    		for (i=0;i<dpixels_pr.length;i++)  dmask[i]=fpixels_y0[i];
    		double [] dpixels_pr_dark=dpixels_pr.clone();
    		double [] dpixels_pb_dark=dpixels_pb.clone();
    		gb.blurDouble(dmask, width, height, colorProcParameters.maskSigma, colorProcParameters.maskSigma, 0.01);
    		gb.blurDouble(dpixels_pr, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    		gb.blurDouble(dpixels_pb, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    		gb.blurDouble(dpixels_pr_dark, width, height, colorProcParameters.chromaDarkSigma, colorProcParameters.chromaDarkSigma, 0.01);
    		gb.blurDouble(dpixels_pb_dark, width, height, colorProcParameters.chromaDarkSigma, colorProcParameters.chromaDarkSigma, 0.01);
    		if (debugLevel> 2) {
    			SDFA_INSTANCE.showArrays(dmask, width, height,"dmask");
    			SDFA_INSTANCE.showArrays(dpixels_pr, width, height,"dpixels_pr");
    			SDFA_INSTANCE.showArrays(dpixels_pb, width, height,"dpixels_pb");
    			SDFA_INSTANCE.showArrays(dpixels_pr_dark, width, height,"dpixels_pr_dark");
    			SDFA_INSTANCE.showArrays(dpixels_pb_dark, width, height,"dpixels_pb_dark");
    		}
    		double mp;
    		double k =1.0/(colorProcParameters.maskMax-colorProcParameters.maskMin);
    		for (i=0;i<dmask.length;i++) {
    			mp=dmask[i];
    			if (mp < colorProcParameters.maskMin) {
    				dmask[i]=0.0;
    			} else if (mp< colorProcParameters.maskMax) {
    				dmask[i]= k*(mp-colorProcParameters.maskMin);
    			} else dmask[i]=1.0;
    		}
    		//TODO: null DENOISE_MASK if it is not calculated
    		if (colorProcParameters.combineWithSharpnessMask) {
    			if (denoiseMask==null) {
    				System.out.println ( "Can not combine masks as denoiseMask is null (i.e. no denoise was performed)"); // here
    			} else if (denoiseMask.length!=dmask.length) {
    				System.out.println ( "Can not combine masks as denoiseMask length is different from that of dmask");
    			} else {
    				for (i=0;i<dmask.length;i++) {
    					dmask[i]+=denoiseMask[i];
    					if (dmask[i]>1.0) dmask[i]=1.0;
    				}
    			}

    		}
    		for (i=0;i<dmask.length;i++) {
    			mp=dmask[i];
    			dpixels_pb[i]= (1.0-mp)*dpixels_pb_dark[i]+ mp* dpixels_pb[i];
    			dpixels_pr[i]= (1.0-mp)*dpixels_pr_dark[i]+ mp* dpixels_pr[i];
    		}
    		this.denoiseMaskChroma=dmask; // (global, used to return denoise mask to save/show
    		this.denoiseMaskChromaWidth=width; // width of the this.denoiseMaskChroma image
    	} else {
    		gb.blurDouble(dpixels_pr, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    		gb.blurDouble(dpixels_pb, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
    		this.denoiseMaskChroma=null; // (global, used to return denoise mask to save/show
    	}
    	for (i=0;i<dpixels_pr.length;i++) {
    		fpixels_pr[i]=(float) dpixels_pr[i];
    		fpixels_pb[i]=(float) dpixels_pb[i];
    	}
    	stack.addSlice("Pr", fpixels_pr);
    	stack.addSlice("Pb", fpixels_pb);
    	stack.addSlice("Y",  fpixels_y);
    	stack.addSlice("Y0", fpixels_y0); // not filtered by low-pass, preliminary (for comparison only)
    	if (debugLevel>2) {
    		stack.addSlice("Yr",fpixels_yR);
    		stack.addSlice("Yg",fpixels_yG);
    		stack.addSlice("Yb",fpixels_yB);
    	}

    }

// old versiion
    public void  processColorsWeightsOld(ImageStack stack,
    		  double scale,     // initial maximal pixel value (16))
    		  EyesisCorrectionParameters.ColorProcParameters  colorProcParameters,
    		  CorrectionColorProc.ColorGainsParameters channelGainParameters,
    		  int channel,
    		  double [] denoiseMask,
    		  int debugLevel
      ) {
    	  double thisGain=       colorProcParameters.gain;
    	  double thisBalanceRed= colorProcParameters.balanceRed;
    	  double thisBalanceBlue=colorProcParameters.balanceBlue;
    	  if (channelGainParameters!=null) {
    		  thisGain*=       channelGainParameters.gain[channel];
    		  thisBalanceRed*= channelGainParameters.balanceRed[channel];
    		  thisBalanceBlue*=channelGainParameters.balanceBlue[channel];
    	  }
          float [] fpixels_r= (float[]) stack.getPixels(1);
          float [] fpixels_g= (float[]) stack.getPixels(2);
          float [] fpixels_b= (float[]) stack.getPixels(3);
          boolean useWeights=(stack.getSize()>=5);
          if (!useWeights) {
            stack.addSlice("dummy1",  fpixels_r);
            stack.addSlice("dummy2",  fpixels_g);
          }
          float [] fpixels_wr=(float[]) stack.getPixels(4);
          float [] fpixels_wb=(float[]) stack.getPixels(5);
          int length=fpixels_r.length;
          int width= stack.getWidth();
          int height=stack.getHeight();
    /* Scale colors, gamma-convert */
          int i;
          double gain_red= thisBalanceRed* thisGain/scale;
          double gain_blue=thisBalanceBlue*thisGain/scale;
          double gain_green=thisGain/scale;
          double gamma_a=Math.pow(colorProcParameters.minLin,colorProcParameters.gamma)*(1.0-colorProcParameters.gamma);
          gamma_a=gamma_a/(1.0-gamma_a);
          double gamma_linK=(1.0+gamma_a)*colorProcParameters.gamma*Math.pow(colorProcParameters.minLin,colorProcParameters.gamma)/colorProcParameters.minLin;

          for (i=0;i<length;i++) {
            fpixels_r[i]=(float) linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, fpixels_r[i]*gain_red);
            fpixels_g[i]=(float) linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, fpixels_g[i]*gain_green);
            fpixels_b[i]=(float) linGamma(colorProcParameters.gamma, gamma_a, gamma_linK, colorProcParameters.minLin, fpixels_b[i]*gain_blue);
          }

    /* Convert to YPbPr */
          double Y,Pb,Pr;
          double Kg=1.0-colorProcParameters.kr-colorProcParameters.kb;
          double Sb=0.5/(1.0-colorProcParameters.kb)*colorProcParameters.saturationBlue;
          double Sr=0.5/(1.0-colorProcParameters.kr)*colorProcParameters.saturationRed;
          double Yr,Yg,Yb,Wr,Wg,Wb,S;
    /* coefficients to find Y from Pb, Pr and a color (R,G or B)
     Yr = R- Pr*KPrR
     Yb = B- Pb*KPbB
     Yg = G+ Pr*KPrG  + Pb*KPbG
     */
          double KPrR= -(2.0*(1-colorProcParameters.kr))/colorProcParameters.saturationRed;
          double KPbB= -(2.0*(1-colorProcParameters.kb))/colorProcParameters.saturationBlue;
          double KPrG=  2.0*colorProcParameters.kr*(1-colorProcParameters.kr)/Kg/colorProcParameters.saturationRed;
          double KPbG=  2.0*colorProcParameters.kb*(1-colorProcParameters.kb)/Kg/colorProcParameters.saturationBlue;
          if (debugLevel > 1) {
              System.out.println ( " processColorsWeights() gain_red="+gain_red+" gain_green="+gain_green+" gain_blue="+gain_blue);
              System.out.println ( " processColorsWeights() gamma="+colorProcParameters.gamma+      " minLin="+colorProcParameters.minLin+" gamma_a="+gamma_a+" gamma_linK="+gamma_linK);
              System.out.println ( " processColorsWeights() Kr="+colorProcParameters.kr+" Kg="+Kg+" Kb="+colorProcParameters.kb+" Sr="+Sr+" Sb="+Sb);
              System.out.println ( " processColorsWeights() KPrR="+KPrR+" KPbB="+KPbB+" KPrG="+KPrG+" KPbG="+KPbG);
          }

          float [] fpixels_pb= new float [length];
          float [] fpixels_pr= new float [length];
          float [] fpixels_y0= new float [length];
          float [] fpixels_y= fpixels_y0;

          float [] fpixels_yR=null;
          float [] fpixels_yG=null;
          float [] fpixels_yB=null;

          if (debugLevel>2) {
            fpixels_yR= new float [length];
            fpixels_yG= new float [length];
            fpixels_yB= new float [length];
          }
          for (i=0;i<length;i++) {
            Y=colorProcParameters.kr*fpixels_r[i]+Kg*fpixels_g[i]+colorProcParameters.kb*fpixels_b[i];
            fpixels_pb[i] = (float) (Sb*(fpixels_b[i]-Y));
            fpixels_pr[i] = (float) (Sr*(fpixels_r[i]-Y));
            fpixels_y0[i]=(float) Y;
          }
    /* calculate Y from weighted colors, weights derived from how good each color component predicts signal in each subpixel of Bayer pattern */
          if (useWeights) {
            fpixels_y=  new float [length];
            for (i=0;i<length;i++) {
              Pb=fpixels_pb[i];
              Pr=fpixels_pr[i];
              Yr = fpixels_r[i]- Pr*KPrR;
              Yb = fpixels_b[i]- Pb*KPbB;
              Yg = fpixels_g[i]+ Pr*KPrG  + Pb*KPbG;
              Wr=fpixels_wr[i];
              Wb=fpixels_wb[i];
              Wg=1.0-Wr-Wb;
              S=1.0/(Wr*(colorProcParameters.weightScaleR-1.0)+Wb*(colorProcParameters.weightScaleB-1.0)+1.0);
              Wr*=S*colorProcParameters.weightScaleR;
              Wb*=S*colorProcParameters.weightScaleB;
              Wg*=S;
              Y=Yr*Wr+Yg*Wg+Yb*Wb;
              fpixels_y[i]=(float) Y;
              if (debugLevel>2) {
                fpixels_yR[i]= (float) Yr;
                fpixels_yG[i]= (float) Yg;
                fpixels_yB[i]= (float) Yb;
              }
            }
          }
    /* Low-pass filter Pb and Pr */
    	  DoubleGaussianBlur gb=new DoubleGaussianBlur();
          double [] dpixels_pr=new double[fpixels_pr.length];
          double [] dpixels_pb=new double[fpixels_pb.length];
          for (i=0;i<dpixels_pr.length;i++) {
        	  dpixels_pr[i]=fpixels_pr[i];
        	  dpixels_pb[i]=fpixels_pb[i];
          }
          if (colorProcParameters.maskMax>0.0) {
              double [] dmask=new double[fpixels_y0.length];
              for (i=0;i<dpixels_pr.length;i++)  dmask[i]=fpixels_y0[i];
              double [] dpixels_pr_dark=dpixels_pr.clone();
              double [] dpixels_pb_dark=dpixels_pb.clone();
        	  gb.blurDouble(dmask, width, height, colorProcParameters.maskSigma, colorProcParameters.maskSigma, 0.01);
        	  gb.blurDouble(dpixels_pr, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
        	  gb.blurDouble(dpixels_pb, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
        	  gb.blurDouble(dpixels_pr_dark, width, height, colorProcParameters.chromaDarkSigma, colorProcParameters.chromaDarkSigma, 0.01);
        	  gb.blurDouble(dpixels_pb_dark, width, height, colorProcParameters.chromaDarkSigma, colorProcParameters.chromaDarkSigma, 0.01);
              if (debugLevel>2) {
            	  SDFA_INSTANCE.showArrays(dmask, width, height,"dmask");
            	  SDFA_INSTANCE.showArrays(dpixels_pr, width, height,"dpixels_pr");
            	  SDFA_INSTANCE.showArrays(dpixels_pb, width, height,"dpixels_pb");
            	  SDFA_INSTANCE.showArrays(dpixels_pr_dark, width, height,"dpixels_pr_dark");
            	  SDFA_INSTANCE.showArrays(dpixels_pb_dark, width, height,"dpixels_pb_dark");
              }
              double mp;
              double k =1.0/(colorProcParameters.maskMax-colorProcParameters.maskMin);
              for (i=0;i<dmask.length;i++) {
            	  mp=dmask[i];
            	  if (mp < colorProcParameters.maskMin) {
            		  dmask[i]=0.0;
            	  } else if (mp< colorProcParameters.maskMax) {
            		  dmask[i]= k*(mp-colorProcParameters.maskMin);
            	  } else dmask[i]=1.0;
              }
    //TODO: null DENOISE_MASK if it is not calculated
              if (colorProcParameters.combineWithSharpnessMask) {
            	  if (denoiseMask==null) {
                      System.out.println ( "Can not combine masks as denoiseMask is null (i.e. no denoise was performed)");
            	  } else if (denoiseMask.length!=dmask.length) {
                      System.out.println ( "Can not combine masks as denoiseMask length is different from that of dmask");
            	  } else {
                      for (i=0;i<dmask.length;i++) {
                    	  dmask[i]+=denoiseMask[i];
                    	  if (dmask[i]>1.0) dmask[i]=1.0;
                      }
            	  }

              }
              for (i=0;i<dmask.length;i++) {
            	  mp=dmask[i];
          		  dpixels_pb[i]= (1.0-mp)*dpixels_pb_dark[i]+ mp* dpixels_pb[i];
          		  dpixels_pr[i]= (1.0-mp)*dpixels_pr_dark[i]+ mp* dpixels_pr[i];
              }
              this.denoiseMaskChroma=dmask; // (global, used to return denoise mask to save/show
              this.denoiseMaskChromaWidth=width; // width of the this.denoiseMaskChroma image
          } else {
        	  gb.blurDouble(dpixels_pr, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
        	  gb.blurDouble(dpixels_pb, width, height, colorProcParameters.chromaBrightSigma, colorProcParameters.chromaBrightSigma, 0.01);
       	      this.denoiseMaskChroma=null; // (global, used to return denoise mask to save/show
          }
          for (i=0;i<dpixels_pr.length;i++) {
        	  fpixels_pr[i]=(float) dpixels_pr[i];
        	  fpixels_pb[i]=(float) dpixels_pb[i];
          }
          stack.addSlice("Pr", fpixels_pr);
          stack.addSlice("Pb", fpixels_pb);
          stack.addSlice("Y",  fpixels_y);
          stack.addSlice("Y0", fpixels_y0); // not filtered by low-pass, preliminary (for comaprison only)
          if (debugLevel>2) {
            stack.addSlice("Yr",fpixels_yR);
            stack.addSlice("Yg",fpixels_yG);
            stack.addSlice("Yb",fpixels_yB);
          }

      }


    /* ======================================================================== */
    public double linGamma(double gamma, double a, double k, double x0, double x) {
    	if (x<0) return 0.0;

    	if (x<=x0) return k*x;
    	return (1.0+a)*Math.pow(x,gamma)-a;
    	//  return x;

    	// individual per-channel color balance and gain
    }
    public void  YPrPbToRGB(ImageStack stack,
    		double Kr,        // 0.299;
    		double Kb,        // 0.114;
    		int sliceY,
    		int slicePr,
    		int slicePb
    ) {
    	float [] fpixels_r= (float[]) stack.getPixels(1);
    	float [] fpixels_g= (float[]) stack.getPixels(2);
    	float [] fpixels_b= (float[]) stack.getPixels(3);
    	float [] fpixels_Y= (float[]) stack.getPixels(sliceY);
    	float [] fpixels_Pr=(float[]) stack.getPixels(slicePr);
    	float [] fpixels_Pb=(float[]) stack.getPixels(slicePb);
    	int length=fpixels_r.length;
    	double Kg=1.0-Kr-Kb;
    	int i;
    	/**
R= Y+ Pr*2.0*(1-Kr)
B= Y+ Pb*2.0*(1-Kb)
G= Y  +Pr*(- 2*Kr*(1-Kr))/Kg + Pb*(-2*Kb*(1-Kb))/Kg

    	 */
    	double KPrR=  2.0*(1-Kr);
    	double KPbB=  2.0*(1-Kb);
    	double KPrG= -2.0*Kr*(1-Kr)/Kg;
    	double KPbG= -2.0*Kb*(1-Kb)/Kg;
    	double Y,Pr,Pb;
    	for (i=0;i<length;i++) {
    		Pb=fpixels_Pb[i];
    		Pr=fpixels_Pr[i];
    		Y =fpixels_Y [i];
    		fpixels_r[i]=(float) (Y+ Pr*KPrR);
    		fpixels_b[i]=(float) (Y+ Pb*KPbB);
    		fpixels_g[i]=(float) (Y+ Pr*KPrG + Pb*KPbG);
    	}
    }




    public static class ColorGainsParameters {
    	public static final String AUX_PREFIX = "AUX-";
    	public double[] gain={
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    	public double[] balanceRed={
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
    	public double[] balanceBlue={
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    			1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};


    	public ColorGainsParameters(){
    	}
    	public void setProperties(String prefix,Properties properties){
    		properties.setProperty(prefix+"channels",this.gain.length+"");
    		for (int i=0; i<this.gain.length;i++) {
    			properties.setProperty(prefix+"gain_"+i,this.gain[i]+"");
    			properties.setProperty(prefix+"balanceRed_"+i,this.balanceRed[i]+"");
    			properties.setProperty(prefix+"balanceBlue_"+i,this.balanceBlue[i]+"");
    		}
    	}
    	public boolean getProperties(String prefix,Properties properties){
    		if (properties.getProperty(prefix+"channels")!=null) {
        		int numChannels=Integer.parseInt(properties.getProperty(prefix+"channels"));
        		this.gain=       new double[numChannels];
        		this.balanceRed= new double[numChannels];
        		this.balanceBlue=new double[numChannels];
        		for (int i=0;i<numChannels;i++) {
        			this.gain[i]=       Double.parseDouble(properties.getProperty(prefix+"gain_"+i));
        			this.balanceRed[i]= Double.parseDouble(properties.getProperty(prefix+"balanceRed_"+i));
        			this.balanceBlue[i]=Double.parseDouble(properties.getProperty(prefix+"balanceBlue_"+i));
        		}
        		return true;
    		}
    		return false;
    	}

    	public void modifyNumChannels(int numChannels){
    		if ((numChannels>0) && (numChannels!=this.gain.length)) {
    			double [] gain1=this.gain;
    			double [] balanceRed1=this.balanceRed;
    			double [] balanceBlue1=this.balanceBlue;
    			this.gain=       new double[numChannels];
    			this.balanceRed= new double[numChannels];
    			this.balanceBlue=new double[numChannels];
    			for (int i=0;i<numChannels;i++) {
    				int j=i;
    				if (j>=gain1.length) j=gain1.length-1;
    				this.gain[i]=gain1[j];
    				this.balanceRed[i]=balanceRed1[j];
    				this.balanceBlue[i]=balanceBlue1[j];
    			}
    		}
    	}

    	public boolean showDialog() {
    		GenericDialog gd = new GenericDialog("Individual channels colors/gains");
    		for (int i =0; i<this.gain.length;i++){
    			gd.addMessage(String.format("=== CHANNEL %02d ===",i));
    			gd.addNumericField(String.format("%02d: Gain (brightness)",i), this.gain[i], 3);
    			gd.addNumericField(String.format("%02d: Balance Red/Green",i), this.balanceRed[i], 3);
    			gd.addNumericField(String.format("%02d: Balance Blue/Green",i), this.balanceBlue[i], 3);
    		}
    		WindowTools.addScrollBars(gd);
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		for (int i =0; i<this.gain.length;i++){
    			this.gain[i]=       gd.getNextNumber();
    			this.balanceRed[i]= gd.getNextNumber();
    			this.balanceBlue[i]=gd.getNextNumber();
    		}
    		return true;
    	}

    	public boolean showDialog(ColorGainsParameters aux) {

    		GenericJTabbedDialog gd = new GenericJTabbedDialog("Individual channels colors/gains", 600,1000);
    		if (aux != null) gd.addTab("Main camera");
    		for (int i =0; i<this.gain.length;i++){
    			gd.addMessage(String.format("=== CHANNEL %02d ===",i));
    			gd.addNumericField(String.format("%02d: Gain (brightness)",i), this.gain[i], 3);
    			gd.addNumericField(String.format("%02d: Balance Red/Green",i), this.balanceRed[i], 3);
    			gd.addNumericField(String.format("%02d: Balance Blue/Green",i), this.balanceBlue[i], 3);
    		}
    		if (aux != null) {
    			gd.addTab("Auxiliary camera");
        		for (int i =0; i<this.gain.length;i++){
        			gd.addMessage(String.format("=== CHANNEL %02d ===",i));
        			gd.addNumericField(String.format("%02d: Gain (brightness)",i), aux.gain[i], 3);
        			gd.addNumericField(String.format("%02d: Balance Red/Green",i), aux.balanceRed[i], 3);
        			gd.addNumericField(String.format("%02d: Balance Blue/Green",i), aux.balanceBlue[i], 3);
        		}

    		}
    		gd.showDialog();
    		if (gd.wasCanceled()) return false;
    		for (int i =0; i<this.gain.length;i++){
    			this.gain[i]=       gd.getNextNumber();
    			this.balanceRed[i]= gd.getNextNumber();
    			this.balanceBlue[i]=gd.getNextNumber();
    		}
       		if (aux != null) {
        		for (int i =0; i<this.gain.length;i++){
        			aux.gain[i]=       gd.getNextNumber();
        			aux.balanceRed[i]= gd.getNextNumber();
        			aux.balanceBlue[i]=gd.getNextNumber();
        		}
       		}
    		return true;
    	}


    }




}
