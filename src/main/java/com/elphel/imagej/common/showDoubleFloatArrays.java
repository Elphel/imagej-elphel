package com.elphel.imagej.common;
import java.awt.Rectangle;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

  public class showDoubleFloatArrays {
	  // defaults for color conversion
	  public int sliceRminusG=1;
	  public int sliceBminusG=2;
	  public double colorSpan=-.25;
	  public double brightnessModulate=1.0;
	  public double Kr=0.299;   // 0.299;
	  public double Kb=0.114;   // 0.114;
	  public double brightness=0.5;



/* For square arrays */

  public void showArrays(double[][] pixels,  String title) { showArrays(pixels, false, title);}
  public void showArrays(double[][] pixels,  boolean asStack, String title) {
    int width=0;
    int i;
    if (pixels==null) return;
    for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
      width= (int) Math.sqrt(pixels[i].length);
      break;
    }
    showArrays(pixels, width, width, asStack, title);
  }

  public void showArrays(float[][] pixels,  String title) { showArrays(pixels, false, title);}
  public void showArrays(float[][] pixels,  boolean asStack, String title) {
    int width=0;
    int i;
    if (pixels==null) return;
    for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
      width= (int) Math.sqrt(pixels[i].length);
      break;
    }
    showArrays(pixels, width, width, asStack, title);
  }


  public void showArrays(double[] pixels,  String title) {
   int width;
    if (pixels!=null) {
      width=(int) Math.sqrt(pixels.length);
      showArrays(pixels,  width,width, title);
    }
  }

  public void showArrays(float[] pixels,  String title) {
   int width;
    if (pixels!=null) {
      width=(int) Math.sqrt(pixels.length);
      showArrays(pixels,  width,width, title);
    }
  }


  public void showArrays(double[][] pixels, int width, int height,  boolean asStack, String title) {
    int i,j;
    if (asStack) {
      float [] fpixels;
      ImageStack array_stack=new ImageStack(width,height);
      for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
          if (pixels[i].length!=(width*height)){
        	  System.out.println("showArrays(): pixels["+i+"].length="+pixels[i].length+" != width (+"+width+") * height("+height+")="+(width*height));
        	  return;
          }
        fpixels=new float[pixels[i].length];
        for (j=0;j<fpixels.length;j++) fpixels[j]=(float)pixels[i][j];
        array_stack.addSlice("chn-"+i,    fpixels);
        if (pixels[i].length!=(width*height)){

        }
      }
      ImagePlus imp_stack = new ImagePlus(title, array_stack);
      imp_stack.getProcessor().resetMinAndMax();
      imp_stack.show();
      return;
    } else showArrays(pixels, width, height, title);
  }
  public void showArrays(double[][] pixels, int width, int height,  boolean asStack, String title, String [] titles) {
	    int i,j;
	    if (pixels==null) {
	    	System.out.println("showDoubleFloatArrays.showArrays(): - pixel array is null");
	    }
	    if (asStack) {
	      float [] fpixels;
	      ImageStack array_stack=new ImageStack(width,height);
	      boolean not_empty = false;
	      for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
	    	  not_empty = true;
	    	  fpixels=new float[pixels[i].length];
	    	  for (j=0;j<fpixels.length;j++) fpixels[j]=(float)pixels[i][j];
	    	  array_stack.addSlice(titles[i], fpixels);
	      }
	      if (not_empty) {
	    	  ImagePlus imp_stack = new ImagePlus(title, array_stack);
	    	  imp_stack.getProcessor().resetMinAndMax();
	    	  imp_stack.show();
	      }
	      return;
	    } else showArrays(pixels, width, height, titles);
	  }

  public void showArrays(float[][] pixels, int width, int height,  boolean asStack, String title, String [] titles) {
	    int i,j;
	    if (asStack) {
	      float [] fpixels;
	      ImageStack array_stack=new ImageStack(width,height);
	      for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
	        fpixels=new float[pixels[i].length];
	        for (j=0;j<fpixels.length;j++) fpixels[j]=pixels[i][j];
	        array_stack.addSlice(titles[i], fpixels);
	      }
	      ImagePlus imp_stack = new ImagePlus(title, array_stack);
	      imp_stack.getProcessor().resetMinAndMax();
	      imp_stack.show();
	      return;
	    } else showArrays(pixels, width, height, titles);
	  }
  public void showArraysSparse(double[][] pixels, int width, int height,  boolean asStack, String title, String [] titles) {
	  int i,j;
	  float [] fpixels;
	  if (asStack) {
		  ImageStack array_stack=new ImageStack(width,height);
		  for (i=0;i<titles.length;i++) {
			  fpixels=new float[pixels.length];
			  for (j=0;j<fpixels.length;j++) fpixels[j]= (pixels[j]==null)?0.0f: ((float)pixels[j][i]);
			  array_stack.addSlice(titles[i], fpixels);
		  }
		  ImagePlus imp_stack = new ImagePlus(title, array_stack);
		  imp_stack.getProcessor().resetMinAndMax();
		  imp_stack.show();
	  } else {
		  ImagePlus[]      imp= new ImagePlus[titles.length];
		  ImageProcessor[] ip= new ImageProcessor[titles.length];
		  for (i=0;i<titles.length;i++) {
			  fpixels=new float[pixels.length];
			  for (j=0;j<fpixels.length;j++) fpixels[j]= (pixels[j]==null)?0.0f: ((float)pixels[j][i]);
			  ip[i]=new FloatProcessor(width,height);
			  ip[i].setPixels(fpixels);
			  ip[i].resetMinAndMax();
			  imp[i]=  new ImagePlus(title+"_"+i, ip[i]);
			  imp[i].show();
		  }
	  }
  }

  public void showArraysSparse(float[][] pixels, int width, int height,  boolean asStack, String title, String [] titles) {
	  int i,j;
	  float [] fpixels;
	  if (asStack) {
		  ImageStack array_stack=new ImageStack(width,height);
		  for (i=0;i<titles.length;i++) {
			  fpixels=new float[pixels.length];
			  for (j=0;j<fpixels.length;j++) fpixels[j]= (pixels[j]==null)?0.0f: ((float)pixels[j][i]);
			  array_stack.addSlice(titles[i], fpixels);
		  }
		  ImagePlus imp_stack = new ImagePlus(title, array_stack);
		  imp_stack.getProcessor().resetMinAndMax();
		  imp_stack.show();
	  } else {
		  ImagePlus[]      imp= new ImagePlus[titles.length];
		  ImageProcessor[] ip= new ImageProcessor[titles.length];
		  for (i=0;i<titles.length;i++) {
			  fpixels=new float[pixels.length];
			  for (j=0;j<fpixels.length;j++) fpixels[j]= (pixels[j]==null)?0.0f: ((float)pixels[j][i]);
			  ip[i]=new FloatProcessor(width,height);
			  ip[i].setPixels(fpixels);
			  ip[i].resetMinAndMax();
			  imp[i]=  new ImagePlus(title+"_"+i, ip[i]);
			  imp[i].show();
		  }
	  }
  }
  public ImageStack makeStack(double[][] pixels, int width, int height) {
	  return makeStack(pixels, width,height, null);
  }

  public ImageStack makeStack(double[][] pixels, int width, int height, String [] titles) {
	  return makeStack( pixels, width, height, titles, false);
  }
  public ImageStack makeStack(double[][] pixels, int width, int height, String [] titles, boolean noNaN) {
      float [] fpixels;
      ImageStack array_stack=new ImageStack(width,height);
      for (int i=0;i<pixels.length;i++) if (pixels[i]!=null) {
          if (pixels[i].length!=(width*height)){
        	  System.out.println("showArrays(): pixels["+i+"].length="+pixels[i].length+" != width (+"+width+") * height("+height+")="+(width*height));
        	  return null;
          }
        fpixels=new float[pixels[i].length];
        if (noNaN){
        	for (int j=0;j<fpixels.length;j++) fpixels[j]= Double.isNaN(pixels[i][j])? 0.0F: ((float)pixels[i][j]);
        } else {
        	for (int j=0;j<fpixels.length;j++) fpixels[j]=(float)pixels[i][j];
        }
        if (titles!=null){
            array_stack.addSlice(titles[i],    fpixels);
        } else {
        	array_stack.addSlice("chn-"+i,    fpixels);
        }
      }
      return array_stack;
  }

  public ImageStack makeStack(float[][] pixels, int width, int height, String [] titles, boolean noNaN) {
      float [] fpixels;
      ImageStack array_stack=new ImageStack(width,height);
      for (int i=0;i<pixels.length;i++) if (pixels[i]!=null) {
          if (pixels[i].length!=(width*height)){
        	  System.out.println("showArrays(): pixels["+i+"].length="+pixels[i].length+" != width (+"+width+") * height("+height+")="+(width*height));
        	  return null;
          }
        if (noNaN){
            fpixels=new float[pixels[i].length];
        	for (int j=0;j<fpixels.length;j++) fpixels[j]= Float.isNaN(pixels[i][j])? 0.0F: ((float)pixels[i][j]);
        } else {
        	fpixels=pixels[i];
        }
        if (titles!=null){
            array_stack.addSlice(titles[i],    fpixels);
        } else {
        	array_stack.addSlice("chn-"+i,    fpixels);
        }
      }
      return array_stack;
  }


  public ImagePlus [] makeArrays(double[][] pixels, int width, int height, String title) {
	  int i,j;
	  float [] fpixels;
	  ImageProcessor[] ip= new ImageProcessor[pixels.length];
	  ImagePlus[]      imp=new ImagePlus[pixels.length];
	  for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
		  fpixels=new float[pixels[i].length];
		  for (j=0;j<fpixels.length;j++) fpixels[j]=(float)pixels[i][j];
		  ip[i]=new FloatProcessor(width,height);
		  ip[i].setPixels(fpixels);
		  ip[i].resetMinAndMax();
		  imp[i]=  new ImagePlus(title+"_"+i, ip[i]);
	  } else imp[i]=null;
	  return imp;
  }
  public ImagePlus makeArrays(double[][] pixels, int width, int height,  String title, String [] titles) {
	  int i,j;
	  if (pixels==null) {
		  System.out.println("showDoubleFloatArrays.makeArrays(): - pixel array is null");
	  }
	  float [] fpixels;
	  ImageStack array_stack=new ImageStack(width,height);
	  boolean not_empty = false;
	  for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
		  not_empty = true;
		  fpixels=new float[pixels[i].length];
		  for (j=0;j<fpixels.length;j++) fpixels[j]=(float)pixels[i][j];
		  array_stack.addSlice(titles[i], fpixels);
	  }
	  if (not_empty) {
		  ImagePlus imp_stack = new ImagePlus(title, array_stack);
		  imp_stack.getProcessor().resetMinAndMax();
//		  imp_stack.show();
		  return imp_stack;
	  }
	  return null;
  }

  public ImagePlus [] makeArrays(double[][] pixels, int width, int height, String [] titles) {
	  int i,j;
	  float [] fpixels;
	  ImageProcessor[] ip= new ImageProcessor[pixels.length];
	  ImagePlus[]      imp=new ImagePlus[pixels.length];
	  for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
		  fpixels=new float[pixels[i].length];
		  for (j=0;j<fpixels.length;j++) fpixels[j]=(float)pixels[i][j];
		  ip[i]=new FloatProcessor(width,height);
		  ip[i].setPixels(fpixels);
		  ip[i].resetMinAndMax();
		  imp[i]=  new ImagePlus(titles[i], ip[i]);
	  } else imp[i]=null;
	  return imp;
  }
  public ImagePlus [] makeArrays(float[][] pixels, int width, int height, String [] titles) {
	  int i,j;
	  float [] fpixels;
	  ImageProcessor[] ip= new ImageProcessor[pixels.length];
	  ImagePlus[]      imp=new ImagePlus[pixels.length];
	  for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
		  fpixels=new float[pixels[i].length];
		  for (j=0;j<fpixels.length;j++) fpixels[j]=pixels[i][j];
		  ip[i]=new FloatProcessor(width,height);
		  ip[i].setPixels(fpixels);
		  ip[i].resetMinAndMax();
		  imp[i]=  new ImagePlus(titles[i], ip[i]);
	  } else imp[i]=null;
	  return imp;
  }

  public void showArrays(double[][] pixels, int width, int height, String title) {
	  int i;
	  ImagePlus[] imp=makeArrays(pixels, width, height, title);
	  if (imp==null) return;
	  for (i=0;i<imp.length;i++) if (imp[i]!=null) {
		  imp[i].show();
	  }
  }

  public void showArrays(double[][] pixels, int width, int height, String [] titles) {
	  int i;
	  ImagePlus[] imp=makeArrays(pixels, width, height, titles);
	  if (imp==null) return;
	  for (i=0;i<imp.length;i++) if (imp[i]!=null) {
		  imp[i].show();
	  }
  }

  public void showArrays(float[][] pixels, int width, int height, String [] titles) {
	  int i;
	  ImagePlus[] imp=makeArrays(pixels, width, height, titles);
	  if (imp==null) return;
	  for (i=0;i<imp.length;i++) if (imp[i]!=null) {
		  imp[i].show();
	  }
  }


  public ImagePlus makeArrays(double[] pixels, int width, int height, String title) {
	  int j;
	  float [] fpixels;
	  if (pixels!=null) {
		  fpixels=new float[pixels.length];
		  for (j=0;j<pixels.length;j++) fpixels[j]=(float)pixels[j];
		  ImageProcessor ip=new FloatProcessor(width,height);
		  ip.setPixels(fpixels);
		  ip.resetMinAndMax();
		  ImagePlus imp=  new ImagePlus(title, ip);
		  return imp;
	  }
	  return null;
  }

  public ImagePlus makeArrays(int[] pixels, int width, int height, String title) {
	  int j;
	  float [] fpixels;
	  if (pixels!=null) {
		  fpixels=new float[pixels.length];
		  for (j=0;j<pixels.length;j++) fpixels[j]=pixels[j];
		  ImageProcessor ip=new FloatProcessor(width,height);
		  ip.setPixels(fpixels);
		  ip.resetMinAndMax();
		  ImagePlus imp=  new ImagePlus(title, ip);
		  return imp;
	  }
	  return null;
  }

  public ImagePlus makeArrays(boolean[] pixels, int width, int height, String title) {
	  int j;
	  float [] fpixels;
	  if (pixels!=null) {
		  fpixels=new float[pixels.length];
		  for (j=0;j<pixels.length;j++) fpixels[j]=pixels[j]?1.0f:0.0f;
		  ImageProcessor ip=new FloatProcessor(width,height);
		  ip.setPixels(fpixels);
		  ip.resetMinAndMax();
		  ImagePlus imp=  new ImagePlus(title, ip);
		  return imp;
	  }
	  return null;
  }


  public ImagePlus makeArrays(float[] pixels, int width, int height, String title) {
	  int j;
	  float [] fpixels;
	  if (pixels!=null) {
		  fpixels=new float[pixels.length];
		  for (j=0;j<pixels.length;j++) fpixels[j]=pixels[j];
		  ImageProcessor ip=new FloatProcessor(width,height);
		  ip.setPixels(fpixels);
		  ip.resetMinAndMax();
		  ImagePlus imp=  new ImagePlus(title, ip);
		  return imp;
	  }
	  return null;
  }

  public void showArrays(double[] pixels, int width, int height, String title) {
	  ImagePlus imp= makeArrays(pixels, width, height, title);
	  if (imp!=null) imp.show();
   }

  public void showArrays(int[] pixels, int width, int height, String title) {
	  ImagePlus imp= makeArrays(pixels, width, height, title);
	  if (imp!=null) imp.show();
   }

  public void showArrays(boolean[] pixels, int width, int height, String title) {
	  ImagePlus imp= makeArrays(pixels, width, height, title);
	  if (imp!=null) imp.show();
   }

  public void showArrays(float[][] pixels, int width, int height, boolean asStack, String title) {
    int i;
    if (asStack) {
      ImageStack array_stack=new ImageStack(width,height);
      for (i=0;i<pixels.length;i++) if (pixels[i]!=null)  array_stack.addSlice("chn-"+i,    pixels[i]);
      ImagePlus imp_stack = new ImagePlus(title, array_stack);
      imp_stack.getProcessor().resetMinAndMax();
      imp_stack.show();
      return;
    } else showArrays(pixels, width, height, title);
  }

  public void showArrays(float[][] pixels, int width, int height, String title) {
    int i;
    ImageProcessor[] ip= new ImageProcessor[pixels.length];
    ImagePlus[]      imp=new ImagePlus[pixels.length];
    for (i=0;i<pixels.length;i++) if (pixels[i]!=null) {
      ip[i]=new FloatProcessor(width,height);
      ip[i].setPixels(pixels[i]);
      ip[i].resetMinAndMax();
      imp[i]=  new ImagePlus(title+"_"+i, ip[i]);
      imp[i].show();
    }
  }

  public void showArrays(float[] pixels, int width, int height, String title) {
    if (pixels!=null) {
      ImageProcessor ip=new FloatProcessor(width,height);
      ip.setPixels(pixels);
      ip.resetMinAndMax();
      ImagePlus imp=  new ImagePlus(title, ip);
      imp.show();
    }
  }
  public ImagePlus showImageStack(ImageStack stack, String title) {
	  if (stack==null) return null;
      ImagePlus imp_stack = new ImagePlus(title, stack);
      imp_stack.getProcessor().resetMinAndMax();
      imp_stack.show();
      return imp_stack;
  }
  public void showImageStackThree(ImageStack stack, String title) {
	  if (stack==null) return;
      ImageStack stack3=new ImageStack(stack.getWidth(),stack.getHeight());
      float [] fpixels_r= (float[]) stack.getPixels(1);
      float [] fpixels_g= (float[]) stack.getPixels(2);
      float [] fpixels_b= (float[]) stack.getPixels(3);
      stack3.addSlice("red",  fpixels_r);
      stack3.addSlice("green",fpixels_g);
      stack3.addSlice("blue", fpixels_b);

      ImagePlus imp_stack = new ImagePlus(title, stack3);
      imp_stack.getProcessor().resetMinAndMax();
      imp_stack.show();
  }

  // additional methods to show 2-d "flows" in color
  /*
  	  public int sliceRminusG=1;
	  public int sliceBminusG=2;
	  public double colorSpan=-1.0;
	  public double brightnessModulate=0.5;
	  public double Kr=0.299;   // 0.299;
	  public double Kb=0.114;   // 0.114;
  */
  public ImagePlus showFlowFromSlices(ImagePlus imp){
	  if ((imp==null) || (imp.getStackSize()<2))	{
		  String msg="Source image with at least two slices is required";
		  IJ.showMessage("Error",msg);
		  throw new IllegalArgumentException (msg);
	  }

	    GenericDialog gd = new GenericDialog("Select parameters to convert 2 slices (representing a vector field) to color");
	    gd.addMessage("Maximal values for normalization are calculated inside the selected area (or the full image if there is no selection");
		gd.addNumericField("Slice number (1.."+imp.getStackSize()+" to convert to R/G", this.sliceRminusG, 0);
		gd.addNumericField("Slice number (1.."+imp.getStackSize()+" to convert to B/G", this.sliceBminusG, 0);
		gd.addNumericField("Color span - if positive - absolute value of the source data to get saturated color, -1 - calculate maximum", this.colorSpan, 3);
		gd.addNumericField("Modulate brightness: 1.0 - use sqrt(x^2+y^2), 0.0 - constant brightness ", this.brightnessModulate, 3);
		gd.addNumericField("Color conversion coefficient Kr (default =0.299) ", this.Kr, 3);
		gd.addNumericField("Color conversion coefficient Kb (default =0.114) ", this.Kb, 3);
		gd.addNumericField("Brightness (0..1.0)", this.brightness, 3);

	    gd.showDialog();
	    if (gd.wasCanceled()) return null;
	    this.sliceRminusG=     (int) gd.getNextNumber();
	    this.sliceBminusG=     (int) gd.getNextNumber();
	    this.colorSpan=              gd.getNextNumber();
	    this.brightnessModulate=     gd.getNextNumber();
	    this.Kr=                     gd.getNextNumber();
	    this.Kb=                     gd.getNextNumber();
	    this.brightness=             gd.getNextNumber();
		  return showFlowFromSlices(
				  imp,
				  this.sliceRminusG,
				  this.sliceBminusG,
				  this.colorSpan,
				  this.brightnessModulate,
				  this.brightness,
				  this.Kr,
				  this.Kb);
  }
	/**
	 *
	 * @param imp image containing at least 2 slices, rectangular selection will be used to find min/max
	 * @param sliceRminusG slice number to be interpreted as R-G
	 * @param sliceBminusG slice number to be interpreted as B-G
	 * @param colorSpan    when positive - absolute value in slice to have full color saturation, if negative -1 - full saturation
	 * @param brightnessModulate when 0 - same brightness, 1.0 - equals to sqrt(slice1^2 + slice2^2), normalized
	 * @return RGB24 color image
	 */
  public ImagePlus showFlowFromSlices(
		  ImagePlus imp,
		  int sliceRminusG,
		  int sliceBminusG,
		  double colorSpan,
		  double brightnessModulate){
	  return showFlowFromSlices(
			  imp,
			  sliceRminusG,
			  sliceBminusG,
			  colorSpan,
			  brightnessModulate,0.5,0.299,0.114);

  }
  public ImagePlus showFlowFromSlices(
		  ImagePlus imp,
		  int sliceRminusG,
		  int sliceBminusG,
		  double colorSpan,
		  double brightnessModulate,
		  double brightness,
		  double Kr,        // 0.299;
          double Kb        // 0.114;

  ){

	  if (	(sliceRminusG<1) || (sliceRminusG>imp.getStackSize()) ||
			  (sliceBminusG<1) || (sliceBminusG>imp.getStackSize()))	{
		  String msg="Source image does not contain specified slices";
		  IJ.showMessage("Error",msg);
		  throw new IllegalArgumentException (msg);
	  }
	  Roi roi= imp.getRoi();
	  Rectangle selection;
	  if (roi==null){
		  selection=new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
	  } else {
		  selection=roi.getBounds();
	  }
	  float [][] pixels = new float[2][];
	  ImageStack stack = imp.getStack();
	  pixels[0]= (float[]) stack.getPixels(sliceRminusG);
	  pixels[1]= (float[]) stack.getPixels(sliceBminusG);
      double maxR=0.0; //Math.abs(pixels[0][0]);
      double maxB=0.0; //Math.abs(pixels[1][0]);
      double maxA=0.0;
      int width=imp.getWidth();
      for (int y=selection.y;y<selection.y+selection.height;y++)
      for (int x=selection.x;x<selection.x+selection.width;x++) {
    	  int index=y*width+x;
    	  if (maxR<Math.abs(pixels[0][index])) maxR=Math.abs(pixels[0][index]);
    	  if (maxB<Math.abs(pixels[1][index])) maxB=Math.abs(pixels[1][index]);
    	  double a=Math.sqrt(pixels[0][index]*pixels[0][index]+pixels[1][index]*pixels[1][index]);
    	  if (maxA<a) maxA=a;

      }
      double maxRB=Math.max(maxR,maxB);
      double Kg=1.0-Kr-Kb;
      double colorScale=1.0;
      if (colorSpan>0){
    	  colorScale=1.0/colorSpan;
      } else {
    	  colorScale=(-colorSpan)* 1.0/maxRB;
      }
/**
R= Y+ Pr*2.0*(1-Kr)
B= Y+ Pb*2.0*(1-Kb)
G= Y  +Pr*(- 2*Kr*(1-Kr))/Kg + Pb*(-2*Kb*(1-Kb))/Kg

*/
      double KPrR=  2.0*(1-Kr);
      double KPbB=  2.0*(1-Kb);
      double KPrG= -2.0*Kr*(1-Kr)/Kg;
      double KPbG= -2.0*Kb*(1-Kb)/Kg;
      System.out.println("maxR="+maxRB+" maxB="+maxB+" maxRB="+maxRB+" maxA="+maxA);
      System.out.println("KPrR="+KPrR+" KPbB="+KPbB+" KPrG="+KPrG+" KPbG="+KPbG);
      int [] pixelsRGB= new int [pixels[0].length];
      double [] rgb=new double[3];
      int debugIndex=selection.y*width+selection.x;
      for (int index=0;index<pixelsRGB.length;index++){
    	  double a=Math.sqrt(pixels[0][index]*pixels[0][index]+pixels[1][index]*pixels[1][index])/maxA;
    	  double Y= brightness*((1.0-brightnessModulate)+brightnessModulate*a);
    	  double Pr=colorScale*pixels[0][index];
    	  double Pb=colorScale*pixels[1][index];
    	  rgb[0]= (Y+ Pr*KPrR);
    	  rgb[1]= (Y+ Pb*KPbB);
    	  rgb[2]= (Y+ Pr*KPrG + Pb*KPbG);
    	  pixelsRGB[index]=0;
		  for (int c=0;c<3;c++) {
			  int d= (int) Math.round(255.0* rgb[c]);
			  if (d>255) d=255;
			  else if (d<0) d=0;
			  pixelsRGB[index]|= d<<(8*(2-c));
		  }
		  if (index==debugIndex){
		      System.out.println("a="+a+" Y="+Y+" Pr="+Pr+" Pb="+Pb);
		      System.out.println("rgb[0]="+rgb[0]+" rgb[1]="+rgb[1]+" rgb[2]="+rgb[2]);
		      System.out.println("pixelsRGB["+index+"]="+pixelsRGB[index]);
		  }
      }
      String title=imp.getTitle()+".png";
	  ColorProcessor cp=new ColorProcessor(imp.getWidth(),imp.getHeight());
	  cp.setPixels(pixelsRGB);
	  ImagePlus imp_color=new ImagePlus(title,cp);
	  return imp_color;

  }


}

