package com.elphel.imagej.calibration;
import java.util.Arrays;
import java.util.Random;

import com.elphel.imagej.common.DoubleFHT;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/*
 **
 ** WavePatternGenerator.java
 **
 ** Copyright (C) 2014 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  WavePatternGenerator.java is free software: you can redistribute it and/or modify
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

public class WavePatternGenerator {
	String DEFAULT_DIRECTORY=null;
//	/    		GenericDialog gd=new GenericDialog("Select images "+startIndex+"..."+(endIndex-1));
	public boolean selectAndGenerate(){
   		GenericDialog gd=new GenericDialog("Select wave pattern parameters");
   		boolean generateCirc=false;
   		boolean generateFFT=false;
   		int width=4096;
   		int height=4096;
   		double period=10.0;
   		double x0=width/2;
   		double y0=height/2;
   		boolean binary=false;

		int fft_size= 4096;
		boolean useSelectedImage=false;
		double fxc=0.15; // center frequency X (0..0.5)
		double fyc=0.05; // center frequency Y (0..0.5)
		double sigma_long=0.001; // coeff for x^2, relative to full size
		double sigma_lat= 0.001; //


   		gd.addCheckbox("Halftone",!binary);
   		gd.addMessage("=== Circular waves generation ===");
   		gd.addCheckbox("Generate circular wave pattern",generateCirc);
		gd.addNumericField("Image width", width, 0, 4, "pixels");
		gd.addNumericField("Image height", height, 0, 4, "pixels");
		gd.addNumericField("Wave period", period, 2, 6, "pixels");
		gd.addNumericField("Wave center X", x0, 2, 10, "pixels");
		gd.addNumericField("Wave center Y", y0, 2, 10, "pixels");
   		gd.addMessage("=== FFT filter ===");
   		gd.addCheckbox("Generate FFT-filtered wave pattern",generateFFT);
   		gd.addCheckbox("Use selected image (otherwise use random data)",useSelectedImage);
		gd.addNumericField("FFT size", fft_size, 0, 4, "pixels");
		gd.addNumericField("Filter center X (frequency domain)", fxc, 5, 7, "fraction");
		gd.addNumericField("Filter center Y (frequency domain)", fyc, 5, 7, "fraction");
		gd.addNumericField("Filter half-width longitudinal", sigma_long, 5, 7, "fraction");
		gd.addNumericField("Filter half-width lateral", sigma_lat, 5, 7, "fraction");
//		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		binary=        !gd.getNextBoolean();
		generateCirc=   gd.getNextBoolean();
		width=    (int) gd.getNextNumber();
		height=   (int) gd.getNextNumber();
		period=         gd.getNextNumber();
		x0=             gd.getNextNumber();
		y0=             gd.getNextNumber();
		generateFFT=    gd.getNextBoolean();
		fft_size= (int) gd.getNextNumber();
		fxc=            gd.getNextNumber();
		fyc=            gd.getNextNumber();
		sigma_long=       gd.getNextNumber();
		sigma_lat=       gd.getNextNumber();

		if (generateCirc) {
			circularPatternGenerator(
					width,
					height,
					period,
					x0,
					y0,
					binary);
		}
		if (generateFFT) {
		FHTFilter (
				fft_size,
				useSelectedImage?WindowManager.getCurrentImage():null, // input image or null (will use random)
				fxc, // center frequency X (0..0.5)
				fyc, // center frequency Y (0..0.5)
				sigma_long, // coeff for x^2, relative to full size
				sigma_lat, //
				binary);
		}
		return true;
	}
	public ImagePlus FHTFilter (
			int size,
			ImagePlus imp_in, // input image or null (will use random)
			double fxc, // center frequency X (0..0.5)
			double fyc, // center frequency Y (0..0.5)
			double sigma_long, // coeff for x^2, relative to full size
			double sigma_lat, //
			boolean binary
			){
// Random random=new Random();
		// double rnd=random.nextDouble();
		double kSigma=5.0;
		double vLen=Math.sqrt(fxc*fxc+fyc*fyc);
		double [] vLong={fxc/vLen, fyc/vLen};
		double [] vLat= {-fyc/vLen, fxc/vLen};
		double [] data=new double[size*size];
		if (imp_in!=null){
			int iWidth=imp_in.getWidth();
			int iHeight=imp_in.getHeight();
			int dy=iHeight/2-size/2;
			int dx=iWidth/2-size/2;
			float [] fpixels= (float[]) imp_in.getProcessor().getPixels();
			for (int y=0;y<size;y++){
				int y_in=y+dy;
				if ((y_in<0) || (y_in>=iHeight)) y_in=-1;
				for (int x=0;x<size;x++){
					if (y_in<0) data[y*size+x]=0.0;
					else {
						int x_in=x+dx;
						if ((x_in<0) || (x_in>=iWidth)) data[y*size+x]=0.0;
						else {

							data[y*size+x]=fpixels[y_in*iWidth+x_in];
						}
					}
				}
			}
		} else {
			Random random=new Random();
			for (int i=0;i<data.length;i++) data[i]=random.nextDouble();
		}
		// generate frequency mask
		double [] mask=new double[size*size];
		Arrays.fill(mask,0.0);
		int ifxc=(int) Math.round(fxc*size);
		int ifyc=(int) Math.round(fyc*size);
		int iRangeX= (int) Math.round((sigma_long*Math.abs(vLong[0]) + sigma_lat*Math.abs(vLat[0]))* size*kSigma);
		int iRangeY= (int) Math.round((sigma_long*Math.abs(vLong[1]) + sigma_lat*Math.abs(vLat[1]))* size*kSigma);
		double kLong= -0.5/(sigma_long*sigma_long*size*size);
		double kLat=  -0.5/(sigma_lat* sigma_lat* size*size);
		for (int y=ifyc-iRangeY;y<ifyc+iRangeY;y++){
			int iy=(y+size)%size;
			int base=size*iy;
			for (int x=ifxc-iRangeX;x<ifxc+iRangeX;x++){
				int ix=(x+size)%size;
				double cLong=(x-ifxc)*vLong[0]+(y-ifyc)*vLong[1];
				double cLat= (x-ifxc)*vLat[0]+ (y-ifyc)*vLat[1];
				mask[base+ix]=Math.exp(kLong*cLong*cLong+kLat*cLat*cLat);
			}
		}
		ShowDoubleFloatArrays.showArrays(data, "input data");
		ShowDoubleFloatArrays.showArrays(mask, "mask");
		DoubleFHT fht =new DoubleFHT();
		fht.swapQuadrants(data);
		fht.transform(    data);
		ShowDoubleFloatArrays.showArrays(data, "FHT data");
		for (int i=0;i<data.length;i++) data[i]*=mask[i];
		ShowDoubleFloatArrays.showArrays(data, "masked FHT data");
		fht.inverseTransform(data);
		fht.swapQuadrants   (data);
//		ShowDoubleFloatArrays.showArrays(data, "restored data");
		float [] pixels = new float [data.length];
		if (binary) for (int i=0;i<data.length;i++) pixels[i] = (data[i]>0)?1.0f:0.0f;
		else        for (int i=0;i<data.length;i++) pixels[i] = (float) data[i];
		ImageProcessor ip_wave = new FloatProcessor(size,size);
		ip_wave.setPixels(pixels);
		ip_wave.resetMinAndMax();
		ImagePlus imp_circle= new ImagePlus("wave_fxc"+fxc+"_fyc"+fyc+"_sigma_long"+sigma_long+"_sigma_lat"+sigma_lat,ip_wave);
		imp_circle.show();
		return imp_circle;
	}

	public ImagePlus circularPatternGenerator(
			int width,
			int height,
			double period,
			double centerX,
			double centerY,
			boolean binary
			){
		ImageProcessor ip_cirWave = new FloatProcessor(width,height);
		float [] pixels = new float [width*height];
		double l=period/Math.PI;
		for (int y=0;y<height;y++){
			int base=y*width;
			for (int x=0;x<width;x++){

			double dx=x-centerX;
			double dy=y-centerY;

			double r=Math.sqrt(dx*dx+dy*dy);
			float d=(float) Math.sin(r/l);
			pixels[base+x]=binary?((d>0)?1.0f:0.0f):d;
			}
		}
		ip_cirWave.setPixels(pixels);
		ip_cirWave.resetMinAndMax();
		ImagePlus imp_circle= new ImagePlus("circularWaves_period"+period+"_xc"+centerX+"_yc"+centerY,ip_cirWave);
		imp_circle.show();
		return imp_circle;
	}

}
