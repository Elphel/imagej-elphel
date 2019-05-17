package com.elphel.imagej.calibration;
/**
** -----------------------------------------------------------------------------**
** Crosstalk_Deconv.java
**
** Calculates 4 kernels for convolution with 4 color components to
** compensate for inter-pixel crosstalk
** Applies the calculated (or direct) kernels to selected image
**
** Copyright (C) 2010 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  Crosstalk_Deconv.java is free software: you can redistribute it and/or modify
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

import ij.*;
import ij.process.*;
import ij.gui.*;

import java.awt.*;
import java.awt.event.*;

import ij.plugin.frame.*;

import java.util.Random;
import java.util.Arrays;

import ij.text.*;

public class Crosstalk_Deconv extends PlugInFrame implements ActionListener {
  /**
	 * 
	 */
	private static final long serialVersionUID = 650738637136446153L;
Panel panel;
  int previousID;
  static Frame instance;
  String crossTalkName="Crosstalk";
  String simulName="Simulation";

 public static int DEBUG_LEVEL = 1;
 public static int MASTER_DEBUG_LEVEL = 1;
 public static boolean in_place = true; /* replace image with the processed */
 public static int BayerMask=1; /* Which components to keep during "Mask Bayer" command (+1 - g1, +2 - r, +4 - b, +8 - g2*/
 public static boolean kernels_valid=false;
 public static boolean monoRed=false; /* use crosstalk for red in all channels */
 public static boolean monoBlue=false; /* use crosstalk for blue in all channels */
 public static double [][] crossCoeff= { /**[red.green,blue][N,NE,E,...NW] */
                                        {0.061,  0.0,  0.091,  0.0,  0.061,  0.0,  0.061,  0.0},
                                        {0.023,  0.0,  0.038,  0.0,  0.023,  0.0,  0.038,  0.0},
                                        {0.0086, 0.0,  0.019,  0.0,  0.0086, 0.0,  0.019,  0.0}};
 public static int FFTSize=16;
 public static int imageFFTSize=128;
 public static int [][][][] invertSeq ; /* once calculated sequence of elements for 4x4 matrix inversion i,j,i,j,i,j,sign */
 public static int [][]     determinantSeq ;      /* once calculated sequence of elements for 4x4 matrix determinant: j,j,j,j,sign */
 public static int invertDimension=0;                           /* determins if 4x4 matrix inversion arrays are initialized */

 public static double    crosstalkThreshold=0.02; /* ignore crosstalk coefficients when abs(ip_RminusB)< is less than crosstalkThreshold */

 public static int       simulWidth=1200;
 public static int       simulHeight=800;
 public static double    simulSlant=0.033; /* pixel shift between rows */
 public static double    simulPeriod=64.1; /* stripes period */
 public static double    simulContrast=0.9; /* blackLevel=simulContrast*whiteLevel  */
 public static double    simulLSFWidth=1.1;
 public static double [] simulBayer={30.0,120.0,1.0,30.0}; /* Gr,R,B,Gb - as with 600x10nm filter - need correction. that was with crosstalk */
 public static double    simulElectrons=66.7; /* 120 - 8000e- */

 String                  patternName="Pattern";
// public boolean          patternHorizontal;
 public static int       patternType=2; /* 0 - vertical, 1 - horizontal, 2 - H+V, 3 - H*V */
 public static int       patternWidth=1920;
 public static int       patternHeight=1200;
 public static int       patternColor=1; /* 1 - red, 2 - green, 4 - blue and their combinartions */
 public static double    patternContrast=1.0; /* stripes contrast */
 public static double    patternPeriod=12.96; /* stripes period */
 public static double    patternNonLinear=0.0; /* non-linearity of the pattern - 0.0 - linear, 1.0 square*/
 public static double    patternSlant=0.033; /* pixel shift between rows */
 public static double    patternScale=0.303; /* sensor pixels per LCD screen pixels */
 public static double    patternV2H=1.0; /* contrast ratio (>1 - V more than H), to simulate astigmatism */
 public static double    pattern3rdHA=0.0; /* simulation: amount of third harmonic, 1.0 - 100% (horizontal)*/
 public static double    pattern3rdHP=0.0; /* simulation: 3-rd harmonic phase,  */
 public static double    pattern3rdVA=0.0; /* simulation: amount of third harmonic, 1.0 - 100% (horizontal)*/
 public static double    pattern3rdVP=0.0; /* simulation: 3-rd harmonic phase,  */
 public static int       patternSimSize=300; /* size of simulated pattern (square, poweer of 2) */
 public static double [] patternSimBayer={30.0,120.0,1.0,30.0}; /* Gr,R,B,Gb as acquiered by simulated sensor */
 public static double    patternSimElectrons=66.7; /* 120 - 8000e- */

 public static double    DC_RminusB,DC_RBminusGG;
 public static double [] DC_bayer;  
 
 private ImageProcessor ip_kr,ip_kg,ip_kb,ip_simul,ip_pattern,ip_RminusB, ip_GbminusGr, ip_patternSim;
 private ImagePlus imp_kr,imp_kg,imp_kb,imp_simul,imp_pattern,imp_RminusB, imp_GbminusGr, imp_patternSim, imp_src;
 
 private FHT fht_kr, fht_kg,fht_kb;
 private FHT fht_RminusB, fht_GbminusGr,fht_crosstalk;
 private static double [] accum_crosstalk; /* accumulating fht pixles here */ 
 private static double [] accum_crosstalk_weight; /* accumulating fht pixles weights (abs of denominators)  here */ 
 private static  double [][][] FKr,FKrX,FKrY,FKrXY,  FKg,FKgX,FKgY,FKgXY,  FKb,FKbX,FKbY,FKbXY;
 private static  double [][][] FRslt_g1, FRslt_r, FRslt_b, FRslt_g2; 
 private static  float[] direct_kr,direct_kg,direct_kb; 
 private static  float[] reverse_kg1,reverse_kr,reverse_kb,reverse_kg2;
 private static float [][][] kernels; // Four kernels for convolution 
// private double [][]    Rslt_g1,  Rslt_r,  Rslt_b,  Rslt_g2; 
// private ImagePlus     Rslt_g1,  Rslt_r,  Rslt_b,  Rslt_g2;
 private FHT     Rslt_g1,  Rslt_r,  Rslt_b,  Rslt_g2;
/*
 private double[][] MS= {{1 ,  1,  1,  1},
                         {1 , -1,  1, -1},
                         {1 ,  1, -1, -1},
                         {1 , -1, -1,  1}};
*/
 private double[][] MS= {{0.5 ,  0.5,  0.5,  0.5},
                         {0.5 , -0.5,  0.5, -0.5},
                         {0.5 ,  0.5, -0.5, -0.5},
                         {0.5 , -0.5, -0.5,  0.5}};

// FHT fht_kr, fht_kg, fht_kb;
// public void run(String arg) {
 public Crosstalk_Deconv() {
    super("Crosstalk Compensation");
    if (IJ.versionLessThan("1.39t")) return;
    if (instance!=null) {
      instance.toFront();
      return;
    }
    instance = this;
    addKeyListener(IJ.getInstance());

    setLayout(new FlowLayout());
    panel = new Panel();
    addButton("Configure");
    addButton("Remove Crosstalk");
    addButton("Add Crosstalk");
    addButton("Create kernels");
    addButton("Mask Bayer");
    addButton("Create Simulation");
    addButton("Create Color Pattern");
    addButton("Simulate Color Pattern");

    addButton("Measure Crosstalk");

    add(panel);
    pack();
    GUI.center(this);
    setVisible(true);
 }



  void addButton(String label) {
    Button b = new Button(label);
    b.addActionListener(this);
    b.addKeyListener(IJ.getInstance());
    panel.add(b);
  }
  public void actionPerformed(ActionEvent e) {
    String label = e.getActionCommand();
    if (label==null) return;
    if (label.equals("Configure")) {
      showDialog();
      return;
    } else if (label.equals("Create kernels")) {
       DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
       createReverseKernels();
       return;
    } else if (label.equals("Create Simulation")) {
       createSimulation();
       return;
    } else if (label.equals("Create Color Pattern")) {
       createColorPattern();
       return;
    } else if (label.equals("Simulate Color Pattern")) {
       simulatePatternCapture();
       return;
    }
/* other commands need current image (float) */
    imp_src = WindowManager.getCurrentImage();
    if (imp_src==null) {
       IJ.showStatus("No image");
       IJ.showMessage("Crosstalk_Deconv Error","Image required");
       return;
    }
    if (imp_src.getType()!=ImagePlus.GRAY32) {
       IJ.showStatus("Converting source image to gray 32 bits (float)");
       new ImageConverter(imp_src).convertToGray32();
    }
    ImageProcessor ip=imp_src.getProcessor();
    ImagePlus imp=imp_src;
    String newTitle= imp_src.getTitle();
    if (label.equals("Measure Crosstalk")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
/*
      Rectangle r=prepareCrosstalkRegion(imp_src);
      if (r.width==0) return;
      imp_src.setRoi(r);
      preprocessCrosstalk(ip,newTitle, false);
      measureCrosstalk(ip_RminusB, ip_GbminusGr,newTitle);
*/
      accumulateCrosstalk (imp_src, imageFFTSize, newTitle, crosstalkThreshold);
      return;
    }
/*
public Rectangle prepareCrosstalkRegion(ImageProcessor ip, String title) {
//  imp_src.setRoi(r);
public Rectangle prepareCrosstalkRegion(ImagePlus imp) {

  Rectangle r=roi_src.getBounds();
public boolean preprocessCrosstalk(ImageProcessor ip, String title) {

*/

/* calculate kernels if they are not current */
    if (!kernels_valid) {
       DEBUG_LEVEL=0;
       createReverseKernels();
    }
   if (label.equals("Remove Crosstalk")) {
//      IJ.showMessage("Crosstalk_Deconv","Remove Crosstalk");
       set2removeCrosstalk();
       newTitle+="_removed_"+crossTalkName;
       IJ.showStatus("Convolving source image with the inverse crosstalk kernel");
    } else if (label.equals("Add Crosstalk")) {
//    IJ.showMessage("Crosstalk_Deconv","Add Crosstalk");
       set2addCrosstalk();
//  String crossTalkName;

       newTitle+="_added_"+crossTalkName;
       IJ.showStatus("Convolving source image with the crosstalk kernel");
    } else if (label.equals("Mask Bayer")) {
//    IJ.showMessage("Crosstalk_Deconv","Add Crosstalk");
       newTitle+="_"+(((BayerMask&1)!=0)?"G1":"")+(((BayerMask&2)!=0)?"R":"")+(((BayerMask&2)!=0)?"B":"")+(((BayerMask&8)!=0)?"G2":"");
       IJ.showStatus("Masking out Bayer components");
    } else return; /* add more options later, if needed */
    if (!in_place) {
      ip=ip.duplicate(); /* create a copy of the image before convolving */
      imp= new ImagePlus(newTitle, ip);
      imp.show();
    } else imp.setTitle(newTitle);
    if (label.equals("Mask Bayer")) {
      maskBayer(ip);
    } else {
      convolveBayerKernel(ip);
    }

    imp.updateAndDraw(); /* Redisplays final image*/
    IJ.showStatus("Crosstalk_Deconv DONE");

  }

  public void processWindowEvent(WindowEvent e) {
    super.processWindowEvent(e);
    if (e.getID()==WindowEvent.WINDOW_CLOSING) {
      instance = null;	
    }
  }

public boolean showColorPatternDialog() {
   double a;
   String [] PatternTypeName={"Vertical","Horizontal","Sum V+H" /*, "Product V*H"*/};

   GenericDialog gd = new GenericDialog("Simulation parameters");
   gd.addStringField ("Pattern name: ",                                patternName, 32);
   gd.addChoice      ("Pattern type", PatternTypeName, PatternTypeName[patternType]);
   gd.addNumericField("Nonlenearity (0 - sine, 1 - square, <0 - odd harmonics):",           patternNonLinear, 3);
   gd.addNumericField("Image width:",                                  patternWidth,  0);
   gd.addNumericField("Image height:",                                 patternHeight, 0);
   gd.addNumericField("Pattern Color (1-red, 2-green, 4 - blue):",     patternColor, 0);
   gd.addNumericField("Contrast (%):",                                 patternContrast*100, 1);
   gd.addNumericField("Stripes period (pixels):",                      patternPeriod, 3);
//   gd.addCheckbox    ("Use square wave instead of sine wave? ", patternSquareWave);

   gd.addNumericField("Pixel shift per row:",                          patternSlant,  3);
   gd.addNumericField("Simulated pattern size (square):",  patternSimSize,0);
   gd.addNumericField("Simulated scale (sensor pixels per LCD pixel):",patternScale,3);
   gd.addNumericField("Relative contrast in vertical direction (to horizontal):",patternV2H,3);
   gd.addNumericField("Horizontal 3-rd harmonic amount:",             pattern3rdHA,3);
   gd.addNumericField("Horizontal 3-rd harmonic phase (degrees):",    pattern3rdHP,3);
   gd.addNumericField("Vertical 3-rd harmonic amount:",               pattern3rdVA,3);
   gd.addNumericField("Vertical 3-rd harmonic phase (degrees):",      pattern3rdVP,3);
   gd.addNumericField("Bayer R  scale :",        patternSimBayer[1], 3);
   gd.addNumericField("Bayer Gr scale :",        patternSimBayer[0], 3);
   gd.addNumericField("Bayer Gb scale :",        patternSimBayer[3], 3);
   gd.addNumericField("Bayer B  scale :",        patternSimBayer[2], 3);
   gd.addNumericField("Electrons per 1.0 in the output (0 - no noise) :", patternSimElectrons, 3);

   gd.showDialog();
   if (gd.wasCanceled()) return false;
   patternName =         gd.getNextString();
//   patternHorizontal=    gd.getNextBoolean();
   patternType=gd.getNextChoiceIndex();
   patternNonLinear =        gd.getNextNumber();
   if (patternNonLinear<-.99) patternNonLinear=-.99;
   else if (patternNonLinear>1.0) patternNonLinear=1.0;
   patternWidth =  (int) gd.getNextNumber();
   patternHeight = (int) gd.getNextNumber();
   patternColor =  (int) gd.getNextNumber();
   a=gd.getNextNumber(); if ((a>=0.0) && (a<=100.0)) patternContrast=0.01*a;
   patternPeriod =       gd.getNextNumber();
//   patternSquareWave=    gd.getNextBoolean();
   patternSlant  =       gd.getNextNumber();
   patternSimSize =  (int) gd.getNextNumber();
   patternScale =        gd.getNextNumber();
   patternV2H =        gd.getNextNumber();
   pattern3rdHA =        gd.getNextNumber();
   pattern3rdHP =        gd.getNextNumber();
   pattern3rdVA =        gd.getNextNumber();
   pattern3rdVP =        gd.getNextNumber();
   patternSimBayer[1] =       gd.getNextNumber();
   patternSimBayer[0] =       gd.getNextNumber();
   patternSimBayer[3] =       gd.getNextNumber();
   patternSimBayer[2] =       gd.getNextNumber();
   patternSimElectrons =       gd.getNextNumber();
   return true;
  }


public void createColorPattern () {
   int i,j,iy,ipix;
   double x,p,wh,wv;
   if (!showColorPatternDialog()) return;
   ip_pattern=new ColorProcessor(patternWidth, patternHeight);
   int[] pixels = (int[])ip_pattern.getPixels();
   ipix=0;
   p=0.0;
   for (i = 0; i < patternHeight; i++) {
     for (j = 0; j < patternWidth; j++) {
      x=j+patternSlant*i;
      x/=patternPeriod;
      x-=Math.floor(x);
      wv=Math.sin(2*Math.PI*x);
      x=(patternHeight-i-1)+patternSlant*j;
      x/=patternPeriod;
      x-=Math.floor(x);
      wh=Math.sin(2*Math.PI*x);
      switch (patternType) {
        case 0:  p=wv; break;
        case 1:  p=wh; break;
        case 2:  p=(wv+wh)/2.0; break;
        case 3:  p=wv*wh; break;
      }
      if (patternNonLinear > 0.0) {
        if (patternNonLinear >= 1.0) p= (p>0.0)?1.0:-1.0;
        else if (p>0.0) p=Math.pow(p,1.0-patternNonLinear);
        else if (p<0.0) p=-Math.pow(-p,1.0-patternNonLinear);
      } else if (patternNonLinear < 0.0) { /* odd harmonics also */
        p= 2.0*Math.pow(0.5*(p+1.0),1.0+patternNonLinear)-1.0;
      }

      iy= (int) (127.5*(1.0+patternContrast* p));
      iy&= 255; // just in case
      pixels[ipix++] =  (255 << 24) | (((patternColor & 1)!=0)? (iy<<16):0) | (((patternColor & 2)!=0)? (iy<< 8):0) | (((patternColor & 4)!=0)? iy:0) ;
//   if ((i==10) && (j==10))   IJ.showMessage("Debug","ipix="+ipix+"\npixels[ipix-1]="+(pixels[ipix-1]&0xffffff)+"\npixels.length="+pixels.length);

     }
   }
   ip_pattern.setPixels(pixels);
   imp_pattern=  new ImagePlus(patternName+"-"+patternType+"-"+patternPeriod+"-"+patternNonLinear, ip_pattern);
   imp_pattern.show();
   IJ.showStatus("Pattren image done.");

}

public void simulatePatternCapture() {
   int i,j,index;
   double x,p,wv,wh;
   double noise=(patternSimElectrons>0)?(1.0/patternSimElectrons):0.0;
   if (!showColorPatternDialog()) return;
   Random generator = new Random( 123456 );
   ip_patternSim= new FloatProcessor(patternSimSize,patternSimSize);
   float[] pixels = (float[])ip_patternSim.getPixels();
   double simPeriod=patternPeriod*patternScale;

   index=0;
   p=0.0;
   double v2h=1.0/(patternV2H + 1.0);
   for (i = 0; i < patternSimSize; i++) {
     for (j = 0; j < patternSimSize; j++) {
      x=j+patternSlant*i;
      x/=simPeriod;
      x-=Math.floor(x);
      wv=Math.sin(2*Math.PI*x);
      if (pattern3rdHA>0.0)  wv=wv*(1-pattern3rdHA)+pattern3rdHA*Math.sin(2*Math.PI*(3*x+pattern3rdHP/360.0)); /* may be less than 1.0 amplitude */
      x=(patternHeight-i-1)+patternSlant*j;
      x/=simPeriod;
      x-=Math.floor(x);
      wh=Math.sin(2*Math.PI*x);
      if (pattern3rdVA>0.0)  wh=wh*(1-pattern3rdVA)+pattern3rdVA*Math.sin(2*Math.PI*(3*x+pattern3rdVP/360.0)); /* may be less than 1.0 amplitude */
      switch (patternType) {
        case 0:  p=wv; break;
        case 1:  p=wh; break;
        case 2:
//               p=0.5*(wv+wh);
                p=v2h*(wh*patternV2H+wv);

                 break;
        case 3:  p=wv*wh; break;
      }
      if (patternNonLinear > 0.0) { /* only even harmonics */
        if (patternNonLinear >= 1.0) p= (p>0.0)?1.0:-1.0;
        else if (p>0.0) p=Math.pow(p,1.0-patternNonLinear);
        else if (p<0.0) p=-Math.pow(-p,1.0-patternNonLinear);
      } else if (patternNonLinear < 0.0) { /* odd harmonics also */
        p= 2.0*Math.pow(0.5*(p+1.0),1.0+patternNonLinear)-1.0;
      }
      p= 0.5*(1.0+patternContrast*p);
      p*=patternSimBayer[((i&1)<<1) +(j&1)];
      if (noise>0) {
        p+=Math.sqrt(noise*p)*generator.nextGaussian();
      }
      pixels[index++]=(float) p;
     }
   }
   ip_patternSim.setPixels(pixels);
   ip_patternSim.resetMinAndMax();

   imp_patternSim=  new ImagePlus("simulated_"+patternName+"-"+patternType+"-"+patternPeriod+"-"+patternNonLinear, ip_patternSim);

   imp_patternSim.show();
   IJ.showStatus("Simulated pattern image done.");
}

public boolean accumulateCrosstalk (ImagePlus imp, int size, String title, double threshold) {
  ImageProcessor ip=imp.getProcessor();
  Roi roi_src= imp.getRoi();
  if (roi_src==null){
    imp.setRoi(0, 0, imp.getWidth(), imp.getHeight());
    roi_src= imp.getRoi();    
  }
  Rectangle r=roi_src.getBounds();
  if ((r.width<size) || (r.height<size)) {
    IJ.showMessage("Error","Selection is too small");
    return false;
  }
/* align to Bayer */
  if ((r.x & 1) !=0) {
    r.width+=1;
    r.x--;
  }
  if ((r.y & 1) !=0) {
    r.height+=1;
    r.y--;
  }
  int dw,dh,nw,nh,iw,ih,i;
  int row, base, l, c;
  dw=r.width-size;
  dh=r.height-size;
  nw=(dw + (size >> 2)) / (size >>1);
  nh=(dh + (size >> 2)) / (size >>1);

/* initialize arrays here */
  accum_crosstalk = new double[size*size];
  accum_crosstalk_weight= new double[size*size];
  for (i=0;i<accum_crosstalk.length;i++) accum_crosstalk[i]=0.0 ; /* is there arrayFill() ? */
  for (i=0;i<accum_crosstalk_weight.length;i++) accum_crosstalk_weight[i]=0.0 ; /* is there arrayFill() ? */
/* silent if (nw>0) || (nh>0) or just iw>0 or ih>0 */
  Rectangle rs= new Rectangle();
  float [] pixels_RminusB;
//  float [] pixels_GbminusGr;
  float [] pixels_crosstalk;
  double p, dn;
  dn=1.0/(nh+1)/(nw+1);
  rs.width=size;
  rs.height=size;
  pixels_crosstalk=null;
  pixels_RminusB=null; 
  for (ih=0;ih<=nh;ih++) for (iw=0;iw<=nw;iw++) {
    rs.x=(r.x+ ((nw>0)?((dw*iw)/nw):0)) & (~1);
    rs.y=(r.y+ ((nh>0)?((dh*ih)/nh):0)) & (~1);
    imp_src.setRoi(rs);
    preprocessCrosstalk(ip, title+"-"+ih+":"+iw, ((ih!=0) || (iw!=0)));
//      measureCrosstalk(ip_RminusB, ip_GbminusGr,newTitle);
// private FHT fht_RminusB, fht_GbminusGr,fht_crosstalk;
    fht_RminusB=new FHT(ip_RminusB);
    fht_GbminusGr=new FHT(ip_GbminusGr);
    fht_RminusB.transform();
    fht_GbminusGr.transform();
    fht_crosstalk=fht_GbminusGr.divide(fht_RminusB); /* values are twice the crosstalk */
    pixels_RminusB=(float[])fht_RminusB.getPixels();
//    pixels_GbminusGr=(float[])fht_GbminusGr.getPixels();
    pixels_crosstalk=(float[]) fht_crosstalk.getPixels();
    for (row=0; row<size; row++) {
      base=row*size;
      for (c=0; c<size; c++) {
        l = ((size-row)%size) * size + (size-c)%size;
        p=Math.sqrt(pixels_RminusB[base+c]*pixels_RminusB[base+c] + pixels_RminusB[l]*pixels_RminusB[l]);
        accum_crosstalk_weight[base+c]+=p;
        accum_crosstalk[base+c]+=p*pixels_crosstalk[base+c];
      }
    }
  }
/* restore, mask*/
  double max =0.0;
  for (row=0; row<size; row++) {
    base=row*size;
    for (c=0; c<size; c++) {
      p=accum_crosstalk_weight[base+c];
      if (p>max) max=p;
    }
  }
/* divide arrays here */
  max*=threshold;
  for (row=0; row<size; row++) {
    base=row*size;
    for (c=0; c<size; c++) {
      l = ((size-row)%size) * size + (size-c)%size;
      p=accum_crosstalk_weight[base+c];
      if (p>=max) {
        pixels_crosstalk[base+c]=(float)(accum_crosstalk[base+c]/accum_crosstalk_weight[base+c]);
        pixels_crosstalk[l]=(float)(accum_crosstalk[l]/accum_crosstalk_weight[l]);
      } else {
        pixels_crosstalk[base+c]=0f;
        pixels_crosstalk[l]=0f;
      }
      pixels_RminusB[base+c]=(float) (accum_crosstalk_weight[base+c]*dn); 
    }
  }
  fht_crosstalk.setPixels(pixels_crosstalk);
  fht_RminusB.setPixels(pixels_RminusB);

/* find results */
  finalizeCrosstalk(fht_RminusB, fht_crosstalk, size, title);
  imp_src.setRoi(r); /* restore origional ROI */
  return true;
}


public void finalizeCrosstalk(FHT fht_RminusB, FHT fht_crosstalk, int size, String title) {
  int precision=3;
  int region=16;
  double [][] K=new double[4][];
  int [][] maxXY=new int[4][2];
  double [] Kswap;
  double [] c;
  double e1,e2;
  int [] maxXYswap;
  float [] pixels_RminusB=(float[]) fht_RminusB.getPixels();
  float [] pixels_crosstalk=(float[]) fht_crosstalk.getPixels();
  if (DEBUG_LEVEL>3) printComplexSubArray(pixels_RminusB, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "RminusB");
  if (DEBUG_LEVEL>2) printComplexSubArray(pixels_crosstalk, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "crosstalk");
/* find two maximums in R-B near vertical lines pattern. One is real, other - alias (gets high energy in Gb-Gr) */
  maxXY[0]=findTwoMax(pixels_RminusB, -(size>>2) , 0, region,  region, size);

  maxXY[1][0]=-(size>>1)-maxXY[0][0];
  maxXY[1][1]=-maxXY[0][1];
  K[0]=getFHTComplexPixel(pixels_crosstalk, maxXY[0][0], maxXY[0][1], size);
  K[1]=getFHTComplexPixel(pixels_crosstalk, maxXY[1][0], maxXY[1][1], size);
/* find two maximums in R-B near horizontal lines pattern. One is real, other - alias (gets high energy in Gb-Gr) */
/* Do it independently, so will work with just vertical or just horizontal lines pattern */
  maxXY[2]=findTwoMax(pixels_RminusB, 0, -(size>>2) , region,  region, size);
  maxXY[3][0]=-maxXY[2][0];
  maxXY[3][1]=-(size>>1)-maxXY[2][1];
  K[2]=getFHTComplexPixel(pixels_crosstalk, maxXY[2][0], maxXY[2][1], size);
  K[3]=getFHTComplexPixel(pixels_crosstalk, maxXY[3][0], maxXY[3][1], size);

/* Find which is real, which is fake */
  if ((K[0][0]*K[0][0]+K[0][1]*K[0][1])>(K[1][0]*K[1][0]+K[1][1]*K[1][1])) { /* swap */
    Kswap=K[0];
    K[0]=K[1];
    K[1]=Kswap;
    maxXYswap=maxXY[0];
    maxXY[0]=maxXY[1];
    maxXY[1]=maxXYswap;
  }
/* Find which is real, which is fake (make it independently) */
  if ((K[2][0]*K[2][0]+K[2][1]*K[2][1])>(K[3][0]*K[3][0]+K[3][1]*K[3][1])) { /* swap */
    Kswap=K[2];
    K[2]=K[3];
    K[3]=Kswap;
    maxXYswap=maxXY[2];
    maxXY[2]=maxXY[3];
    maxXY[3]=maxXYswap;
  }
  c=getFHTComplexPixel(pixels_RminusB, maxXY[0][0], maxXY[0][1], size);
  e1=(Math.sqrt(c[0]*c[0]+c[1]*c[1])/Math.abs(DC_RminusB))/size/size*8.0; /* 8.0 - 8 maximums*/
  c=getFHTComplexPixel(pixels_RminusB, maxXY[2][0], maxXY[2][1], size);
  e2=(Math.sqrt(c[0]*c[0]+c[1]*c[1])/Math.abs(DC_RminusB))/size/size*8.0;
  new TextWindow(title+"_crosstalk_raw_data", "parameter\tvalue",
         "File\t"+title+"\n"+
         "Color\t"+ ((DC_RBminusGG<0)?"Green/white":((DC_RminusB>0)?"Red":"Blue"))+"\n"+
       ((DEBUG_LEVEL>4)?
         ("DC_RBminusGG\t"+DC_RBminusGG+"\n"+
         "DC_RminusB\t"+DC_RminusB+"\n"+
         "DC_bayer[0]\t"+DC_bayer[0]+"\n"+
         "DC_bayer[1]\t"+DC_bayer[1]+"\n"+
         "DC_bayer[2]\t"+DC_bayer[2]+"\n"+
         "DC_bayer[3]\t"+DC_bayer[3]+"\n"):"")+
       ((DEBUG_LEVEL>2)?
         ("K[0]("+maxXY[0][0]+","+maxXY[0][1]+")\t"+IJ.d2s(K[0][0],precision)+ ((K[0][1]>=0)?"+":"")+IJ.d2s(K[0][1],precision)+"i\n"+
         "(K[1]("+maxXY[1][0]+","+maxXY[1][1]+")\t"+IJ.d2s(K[1][0],precision)+ ((K[1][1]>=0)?"+":"")+IJ.d2s(K[1][1],precision)+"i)\n"+
         "K[2]("+maxXY[2][0]+","+maxXY[2][1]+")\t"+IJ.d2s(K[2][0],precision)+ ((K[2][1]>=0)?"+":"")+IJ.d2s(K[2][1],precision)+"i\n"+
         "(K[3]("+maxXY[3][0]+","+maxXY[3][1]+")\t"+IJ.d2s(K[3][0],precision)+ ((K[3][1]>=0)?"+":"")+IJ.d2s(K[3][1],precision)+"i)\n"):"")+
       ((DEBUG_LEVEL>1)?
/* K[0][0] is inverted as y in the image is upside down */
         ("S+N\t"+IJ.d2s(-K[0][0],precision)+"\n"+
         "E-W\t"+IJ.d2s(K[0][1],precision)+"\n"+
         "E+W\t"+IJ.d2s(K[2][0],precision)+"\n"+
         "N-S\t"+IJ.d2s(K[2][1],precision)+"\n"):"")+
         "N\t"+IJ.d2s(0.5*(-K[0][0]+K[2][1]),precision)+"\n"+
         "E\t"+IJ.d2s(0.5*(K[0][1]+K[2][0]),precision)+"\n"+
         "S\t"+IJ.d2s(0.5*(-K[0][0]-K[2][1]),precision)+"\n"+
         "W\t"+IJ.d2s(0.5*(-K[0][1]+K[2][0]),precision)+"\n"+
         "Eh\t"+IJ.d2s(e1,precision)+"\n"+
         "Ev\t"+IJ.d2s(e2,precision)+"\n"+
         "",
 500,300);
  

  if (crosstalkThreshold>0.0) maskFHT (fht_crosstalk, fht_RminusB, crosstalkThreshold); 
//  ImagePlus imp= new ImagePlus("crosstalk_fht", fht_crosstalk);
//  imp.show();
  float [] pixels_crosstalk_masked=(float[]) fht_crosstalk.getPixels();
  if (DEBUG_LEVEL>2)  printComplexSubArray(pixels_crosstalk_masked, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "crosstalk");
  if (DEBUG_LEVEL>1) {
    ImageStack stack_crosstalk =   fht_crosstalk.getComplexTransform();
        ImagePlus imp2 = new ImagePlus(title+"_crosstalk_fht_stack", stack_crosstalk);
        imp2.getProcessor().resetMinAndMax();
        imp2.show();
  }
//fht_RminusB
  if (DEBUG_LEVEL>2) {
    ImagePlus imp3 = new ImagePlus(title+"_accumulated_abs(R-B)", fht_RminusB);
    imp3.getProcessor().resetMinAndMax();
    imp3.show();
  }

}






public Rectangle prepareCrosstalkRegion(ImagePlus imp) {
  Roi roi_src= imp.getRoi();
  if (roi_src==null){
    imp.setRoi(0, 0, imp.getWidth(), imp.getHeight());
    roi_src= imp.getRoi();    
  }
  Rectangle r=roi_src.getBounds();
  if ((r.width<8) || (r.height<8)) {
    IJ.showMessage("Error","Selection is too small");
    r.x=0;r.y=0;r.width=0;r.height=0;
    return r;
  }
/* find largest square power of 2 roi inside the current one */
  if ((r.x & 1) !=0) {
    r.width+=1;
    r.x--;
  }
  if ((r.y & 1) !=0) {
    r.height+=1;
    r.y--;
  }
  int size;
  for (size=1; (size<=r.height) && (size<=r.width); size<<=1);
  size >>=1;
  r.x+=((r.width-size)>>2) <<1;
  r.y+=((r.height-size)>>2) <<1;
  r.width=size;
  r.height=size;
//  imp.setRoi(r);
  return r;
}
  

public void preprocessCrosstalk(ImageProcessor ip, String title, boolean silent) {
  int n,row,col,index, brow,bcol;
  Rectangle r=ip.getRoi();
  int size=r.width;
  float [] pixels=(float[])ip.crop().getPixels();
/* remove DC, apply Hamming */
  DC_bayer=new double[4];
  index=0;
  for (row=0;row<size;row++) {
      brow=row & 1;
      for (col=0;col<size;col++) {
        bcol=col & 1;
        DC_bayer[(brow<<1)+bcol]+=pixels[index++];
      }
  }
  for (n=0;n<4;n++) DC_bayer[n]/=(size*size/4);
  DC_RminusB=DC_bayer[1]-DC_bayer[2];
  DC_RBminusGG=DC_bayer[1]+DC_bayer[2]-DC_bayer[0]-DC_bayer[3];
  double [] hamming=new double[size];
  double a=2.0*Math.PI/size;
  for (col=0;col< size; col++ ) hamming [col]= 0.54-0.46*Math.cos(a*col);
  
  float [] pixels_RminusB=   new float[pixels.length];
  float [] pixels_GbminusGr= new float[pixels.length];
  ip_RminusB=   new FloatProcessor(size,size);
  ip_GbminusGr= new FloatProcessor(size,size);
  index=0;
  for (row=0;row<size;row++) {
      brow=row & 1;
      for (col=0;col<size;col++) {
        bcol=col & 1;
        switch ((brow<<1)+bcol) {
          case 0: /* Gr */
                pixels_GbminusGr[index] = - (float) ((pixels[index]-DC_bayer[0])*hamming[row]*hamming[col]);
                pixels_RminusB[index] = 0.0f;
                break;
          case 1: /* R */
                pixels_GbminusGr[index] = 0.0f;
                pixels_RminusB[index] =  (float) ((pixels[index]-DC_bayer[1])*hamming[row]*hamming[col]);
                break;
          case 2: /* B */
                pixels_GbminusGr[index] = 0.0f;
                pixels_RminusB[index] = - (float) ((pixels[index]-DC_bayer[2])*hamming[row]*hamming[col]);
                break;
          case 3: /* Gb */
                pixels_GbminusGr[index] =   (float) ((pixels[index]-DC_bayer[3])*hamming[row]*hamming[col]);
                pixels_RminusB[index] = 0.0f;
                break;
        
        }
        index++;
      }
    }
    ip_RminusB.setPixels(pixels_RminusB);
    ip_GbminusGr.setPixels(pixels_GbminusGr);
    ip_RminusB.resetMinAndMax();
    ip_GbminusGr.resetMinAndMax();
    if (!silent && (DEBUG_LEVEL>1)) {
      imp_RminusB=   new ImagePlus(title+"_RminusB", ip_RminusB);
      imp_GbminusGr= new ImagePlus(title+"_GbminusGr", ip_GbminusGr);
      imp_RminusB.show();
      imp_GbminusGr.show();
    }
}

public void measureCrosstalk(ImageProcessor ip_RminusB,ImageProcessor ip_GbminusGr, String title) {
  int precision=3;
  fht_RminusB=new FHT(ip_RminusB);
  fht_GbminusGr=new FHT(ip_GbminusGr);
  fht_RminusB.transform();
  fht_GbminusGr.transform();
  int size=ip_RminusB.getWidth();
  int region=16;
  double [][] K=new double[4][];
  int [][] maxXY=new int[4][2];
  double [] Kswap;
  double [] c;
  double e1,e2;
  int [] maxXYswap;
  float [] pixels_RminusB=(float[])fht_RminusB.getPixels();
  float [] pixels_GbminusGr=(float[])fht_GbminusGr.getPixels();
  if (DEBUG_LEVEL>3) printComplexSubArray(pixels_RminusB, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "RminusB");
  if (DEBUG_LEVEL>3) printComplexSubArray(pixels_GbminusGr, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "GrminusGb");

//public FHT divide(FHT fht) {
  fht_crosstalk=fht_GbminusGr.divide(fht_RminusB); /* values are twice the crosstalk */
  float [] pixels_crosstalk=(float[]) fht_crosstalk.getPixels();

  if (DEBUG_LEVEL>2) printComplexSubArray(pixels_crosstalk, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "crosstalk");
/* find two maximums in R-B near vertical lines pattern. One is real, other - alias (gets high energy in Gb-Gr) */

//  maxXY[0]=findTwoMax(pixels_RminusB, -(size>>2) , - (region>>1), region,  region, size); // worked?
  maxXY[0]=findTwoMax(pixels_RminusB, -(size>>2) , 0, region,  region, size);

  maxXY[1][0]=-(size>>1)-maxXY[0][0];
  maxXY[1][1]=-maxXY[0][1];
  K[0]=getFHTComplexPixel(pixels_crosstalk, maxXY[0][0], maxXY[0][1], size);
  K[1]=getFHTComplexPixel(pixels_crosstalk, maxXY[1][0], maxXY[1][1], size);
/* find two maximums in R-B near horizontal lines pattern. One is real, other - alias (gets high energy in Gb-Gr) */
/* Do it independently, so will work with just vertical or just horizontal lines pattern */
  maxXY[2]=findTwoMax(pixels_RminusB, 0, -(size>>2) , region,  region, size);
  maxXY[3][0]=-maxXY[2][0];
  maxXY[3][1]=-(size>>1)-maxXY[2][1];
  K[2]=getFHTComplexPixel(pixels_crosstalk, maxXY[2][0], maxXY[2][1], size);
  K[3]=getFHTComplexPixel(pixels_crosstalk, maxXY[3][0], maxXY[3][1], size);

/* Find which is real, which is fake */
  if ((K[0][0]*K[0][0]+K[0][1]*K[0][1])>(K[1][0]*K[1][0]+K[1][1]*K[1][1])) { /* swap */
    Kswap=K[0];
    K[0]=K[1];
    K[1]=Kswap;
    maxXYswap=maxXY[0];
    maxXY[0]=maxXY[1];
    maxXY[1]=maxXYswap;
  }
/* Find which is real, which is fake (make it independently) */
  if ((K[2][0]*K[2][0]+K[2][1]*K[2][1])>(K[3][0]*K[3][0]+K[3][1]*K[3][1])) { /* swap */
    Kswap=K[2];
    K[2]=K[3];
    K[3]=Kswap;
    maxXYswap=maxXY[2];
    maxXY[2]=maxXY[3];
    maxXY[3]=maxXYswap;
  }
  c=getFHTComplexPixel(pixels_RminusB, maxXY[0][0], maxXY[0][1], size);
  e1=(Math.sqrt(c[0]*c[0]+c[1]*c[1])/Math.abs(DC_RminusB))/size/size*8.0; /* 8.0 - 8 maximums*/
  c=getFHTComplexPixel(pixels_RminusB, maxXY[2][0], maxXY[2][1], size);
  e2=(Math.sqrt(c[0]*c[0]+c[1]*c[1])/Math.abs(DC_RminusB))/size/size*8.0;
//DC_RminusB,DC_RBminusGG
  new TextWindow(title+"_crosstalk_raw_data", "parameter\tvalue",
         "File\t"+title+"\n"+
         "Color\t"+ ((DC_RBminusGG<0)?"Green/white":((DC_RminusB>0)?"Red":"Blue"))+"\n"+
       ((DEBUG_LEVEL>4)?
         ("DC_RBminusGG\t"+DC_RBminusGG+"\n"+
         "DC_RminusB\t"+DC_RminusB+"\n"+
         "DC_bayer[0]\t"+DC_bayer[0]+"\n"+
         "DC_bayer[1]\t"+DC_bayer[1]+"\n"+
         "DC_bayer[2]\t"+DC_bayer[2]+"\n"+
         "DC_bayer[3]\t"+DC_bayer[3]+"\n"):"")+
       ((DEBUG_LEVEL>2)?
         ("K[0]("+maxXY[0][0]+","+maxXY[0][1]+")\t"+IJ.d2s(K[0][0],precision)+ ((K[0][1]>=0)?"+":"")+IJ.d2s(K[0][1],precision)+"i\n"+
         "(K[1]("+maxXY[1][0]+","+maxXY[1][1]+")\t"+IJ.d2s(K[1][0],precision)+ ((K[1][1]>=0)?"+":"")+IJ.d2s(K[1][1],precision)+"i)\n"+
         "K[2]("+maxXY[2][0]+","+maxXY[2][1]+")\t"+IJ.d2s(K[2][0],precision)+ ((K[2][1]>=0)?"+":"")+IJ.d2s(K[2][1],precision)+"i\n"+
         "(K[3]("+maxXY[3][0]+","+maxXY[3][1]+")\t"+IJ.d2s(K[3][0],precision)+ ((K[3][1]>=0)?"+":"")+IJ.d2s(K[3][1],precision)+"i)\n"):"")+
       ((DEBUG_LEVEL>1)?
         ("S+N\t"+IJ.d2s(-K[0][0],precision)+"\n"+
         "E-W\t"+IJ.d2s(K[0][1],precision)+"\n"+
         "E+W\t"+IJ.d2s(K[2][0],precision)+"\n"+
         "N-S\t"+IJ.d2s(K[2][1],precision)+"\n"):"")+
         "N\t"+IJ.d2s(0.5*(-K[0][0]+K[2][1]),precision)+"\n"+
         "E\t"+IJ.d2s(0.5*(K[0][1]+K[2][0]),precision)+"\n"+
         "S\t"+IJ.d2s(0.5*(-K[0][0]-K[2][1]),precision)+"\n"+
         "W\t"+IJ.d2s(0.5*(-K[0][1]+K[2][0]),precision)+"\n"+
         "Eh\t"+IJ.d2s(e1,precision)+"\n"+
         "Ev\t"+IJ.d2s(e2,precision)+"\n"+
         "",
 500,300);
//      sb.append("\t"+IJ.d2s(c[0],precision)+ ((c[1]>0)?"+":"")+ IJ.d2s(c[1],precision)+"i");
  

  if (crosstalkThreshold>0.0) maskFHT (fht_crosstalk, fht_RminusB, crosstalkThreshold); 
//  ImagePlus imp= new ImagePlus("crosstalk_fht", fht_crosstalk);
//  imp.show();
  float [] pixels_crosstalk_masked=(float[]) fht_crosstalk.getPixels();
  if (DEBUG_LEVEL>2)  printComplexSubArray(pixels_crosstalk_masked, -(size>>2) - (region>>1) , - (region>>1), region,  region, size, "crosstalk");
  if (DEBUG_LEVEL>1) {
    ImageStack stack_crosstalk =   fht_crosstalk.getComplexTransform();
        ImagePlus imp2 = new ImagePlus(title+"_crosstalk_fht_stack", stack_crosstalk);
        imp2.getProcessor().resetMinAndMax();
        imp2.show();
 }
}
/* return one pixel as complex, -size/2<=x<=+size/2,-size/2<=y<=+size/2, y negative at row=0 */

public double[] getFHTComplexPixel(float [] fht_pixels, int x, int y, int size) {
   double [] rslt = {0.0,0.0};
   int row= y+(size>>1);
   int col= x+(size>>1);
   int row1= (size-row) %size;
   int col1= (size-col) %size;
//   IJ.showMessage("DEBUG520","x="+x+"\ny="+y+"\nsize="+size+"\nrow="+row+"\ncol="+col+"\nrow1="+row1+"\ncol1="+col1);

   rslt[0]= 0.5*(fht_pixels[row1*size+col1] + fht_pixels[row*size+col]);
   rslt[1]= 0.5*(fht_pixels[row1*size+col1] - fht_pixels[row*size+col]);
   return rslt;
}

public void printComplexSubArray(float [] fht_pixels, int x, int y, int w, int h, int size, String title) {
  int ix,iy;
  int precision=3;
//  int size=Math.sqrt(fht_pixels.length);
  String header=new String();
  StringBuffer sb = new StringBuffer();
  double [] c;
  header+="y\\x";
  for (ix=x; ix< (x+w); ix++ ) header +="\t"+ix;
  header +="\n";
  for (iy=y; iy< (y+h); iy++ ) {
    sb.append(IJ.d2s(iy,0));
    for (ix=x; ix< (x+w); ix++ ) {
      c=getFHTComplexPixel(fht_pixels, ix, iy, size);
      sb.append("\t"+IJ.d2s(c[0],precision)+ ((c[1]>0)?"+":"")+ IJ.d2s(c[1],precision)+"i");
    }
    sb.append("\n");
  }
  new TextWindow(title+"_(csv)", header, sb.toString(), 800,600);
}

public int[] findTwoMax(float [] fht_pixels, int x0, int y0, int w, int h, int size) {
  int [] rslt= new int[2];
  int ix,iy;
  int ix0= x0-(w>>1);
  int ix1= x0+(w>>1);
  int iy1= y0+(w>>1);
  double max=0.0;
  double [] c1,c2;
  double a;
  for (iy=y0; iy<= iy1; iy++ ) {
    for (ix=ix0; ix<= ix1; ix++ ) {
      c1=getFHTComplexPixel(fht_pixels, ix, iy, size);
      c2=getFHTComplexPixel(fht_pixels, (x0<<1)-ix, (y0<<1)-iy, size);
      a=c1[0]*c1[0]+c1[1]*c1[1]+c2[0]*c2[0]+c2[1]*c2[1];
      if (a>max) {
        max=a;
        rslt[0]=ix;
        rslt[1]=iy;
      }
    }
  }
  return rslt;
}






/*
    StringBuffer sb = new StringBuffer();
    for (i=0;i<xValues.length;i++)  {
       sb.append(IJ.d2s(xValues[i],2));
       for (n=0;n<yValues.length;n+=2)  sb.append("\t"+IJ.d2s(yValues[n][i],4)+((yValues[n+1][i]>=0)?"+":"")+IJ.d2s(yValues[n+1][i],4)+"i");
       sb.append("\n");
   }
   TextWindow tw = new TextWindow(title+"_(csv)", hs, sb.toString(), (Width>0)?Width:(84*(yValues.length+1)), (Height>0)?Height:600);

*/


public void maskFHT (FHT fht, FHT fht_mask, double threshold) {
  int row, base, l, c;
  int size=fht.getWidth();
  float max = 0f;
  float p;
  float[] pixels= (float[])fht.getPixels();
  float[] mask=   (float[])fht_mask.getPixels();
  for (row=0; row<size; row++) {
    base=row*size;
    for (c=0; c<size; c++) {
      l = ((size-row)%size) * size + (size-c)%size;
      p=(mask[base+c]*mask[base+c] + mask[l]*mask[l])/2f;
      if (p>max) max=p;
    }
  }
  max*=threshold*threshold;
  for (row=0; row<size; row++) {
    base=row*size;
    for (c=0; c<size; c++) {
      l = ((size-row)%size) * size + (size-c)%size;
      p=(mask[base+c]*mask[base+c] + mask[l]*mask[l])/2f;
      if (p<max) {
        pixels[base+c]=0f;
        pixels[l]=0f;
      }
    }
  }
  fht.setPixels(pixels);
}

 public boolean showSimulationDialog() {
   double a;
   GenericDialog gd = new GenericDialog("Simulation parameters");
   gd.addStringField("Simulation name: ",        simulName, 32);
   gd.addNumericField("Image width:",            simulWidth,  0);
   gd.addNumericField("Image height:",           simulHeight, 0);
   gd.addNumericField("SE pixel shift per row:", simulSlant,  3);
   gd.addNumericField("Stripes period (pixels):",simulPeriod, 3);
   gd.addNumericField("Contrast (%):",           simulContrast*100, 1);
   gd.addNumericField("LSF width (pixels):",     simulLSFWidth, 3);
   gd.addNumericField("Bayer R  scale :",        simulBayer[1], 3);
   gd.addNumericField("Bayer Gr scale :",        simulBayer[0], 3);
   gd.addNumericField("Bayer Gb scale :",        simulBayer[3], 3);
   gd.addNumericField("Bayer B  scale :",        simulBayer[2], 3);
   gd.addNumericField("Electrons per 1.0 in the output (0 - no noise) :", simulElectrons, 3);
   gd.showDialog();
   if (gd.wasCanceled()) return false;
   simulName =         gd.getNextString();
   simulWidth =  (int) gd.getNextNumber();
   simulHeight = (int) gd.getNextNumber();
   simulSlant        = gd.getNextNumber();
   simulPeriod        = gd.getNextNumber();
   a=gd.getNextNumber(); if ((a>=0.0) && (a<=100.0)) simulContrast=0.01*a;
   simulLSFWidth      = gd.getNextNumber();
   simulBayer[1]      = gd.getNextNumber();
   simulBayer[0]      = gd.getNextNumber();
   simulBayer[3]      = gd.getNextNumber();
   simulBayer[2]      = gd.getNextNumber();
   simulElectrons     = gd.getNextNumber();
   return true;
  }




 public void createSimulation() {
   int i,j,index,ix;
   double x,a,b;
   int pattDiv=16;
   double min,max;
   b=2*Math.sqrt(Math.log(2))/simulLSFWidth/pattDiv;
   if (!showSimulationDialog()) return;
//   IJ.showMessage("Debug","Creating simulation image");
   int pattLength= (int)(pattDiv*simulPeriod)+2;
   int pattLength4=pattLength/4;
   double [] pattern=new double [pattLength];
   for (i=0; i<pattLength; i++) pattern[i]=0;
   for (i=pattLength4; i< pattLength- pattLength4; i++ ) {
      for (j=0;j<pattLength4;j++) {
        a=b*j;
        a=Math.exp(-(a*a));
        pattern[i-j]+=a;
        if (j>0) pattern[i+j]+=a;
      }
   }
   min=pattern[0];
   max=pattern[1];
   for (i=1; i<pattLength; i++) {
     if (pattern[i]>max) max=pattern[i];
     if (pattern[i]<max) min=pattern[i];
   }
   a=1.0/(max-min);
   for (i=0; i<pattLength; i++) {
     pattern[i]=(pattern[i]-min)*a;
   }   
   ip_simul = new FloatProcessor(simulWidth,simulHeight);
   float [] pixels= new float [simulWidth *simulHeight];
/**create left side "horizontal" stripes */
   IJ.showStatus("Creating horizontal stripes...");
//   IJ.showMessage("Debug","simulPeriod*pattDiv="+simulPeriod*pattDiv+"\npattern.length="+pattern.length);

   for (i=0;i<(simulWidth>>1); i++ ) {
     index=simulWidth*(simulHeight-1)+i;
     for (j=0;j<simulHeight;j++) {
       x=j+simulSlant*i;
       x-=simulPeriod*Math.floor(x/simulPeriod);
       x*=pattDiv;
       ix= (int) x;
       x-=ix;
//if (ix>=(pattern.length-1)) IJ.showMessage("Debug","ix="+ix+"\nx="+x);
       pixels[index]=(float)(pattern[ix]*(1.0- x) +  x*pattern[ix+1]);
//       pixels[index]=(float)(ix+x);
       index-=simulWidth;
     }
   }
/**create right side "vertical" stripes */
   IJ.showStatus("Creating verical stripes...");

   for (i=0; i<simulHeight; i++ ) {
     index=simulWidth*i+ (simulWidth>>1);
     for (j=0;j<(simulWidth>>1);j++) {
       x=j+simulSlant*i;
       x-=simulPeriod*Math.floor(x/simulPeriod);
       x*=pattDiv;
       ix= (int) x;
       x-=ix;
       pixels[index+j]=(float)(pattern[ix]*(1.0- x) +  x*pattern[ix+1]);
     }
   }
   IJ.showStatus("Applying Bayer colors...");
   for (i=0;i<simulHeight; i++ ) {
     index=simulWidth*i;
     for (j=0;j<simulWidth;j++) {
       ix=((i&1)<<1) +(j&1);
       a=pixels[index+j];
       pixels[index+j]=(float) ((1.0-simulContrast)+simulContrast * a*simulBayer[ix]);
     }
   }

   Random generator = new Random( 123456 );
   if (simulElectrons>0.0) {
     IJ.showStatus("Aplying noise...");
     a=1.0/simulElectrons;
     for (i=0;i<simulHeight; i++ ) {
       index=simulWidth*i;
       for (j=0;j<simulWidth;j++) {
         ix=((i&1)<<1) +(j&1);
         pixels[index+j]+= Math.sqrt(a*pixels[index+j])*generator.nextGaussian();
       }
     }
   }
   ip_simul.setPixels(pixels);
   ip_simul.resetMinAndMax();
   imp_simul= new ImagePlus(simulName, ip_simul);
   imp_simul.show();
   IJ.showStatus("Simulation image done.");
 }




 public boolean showDialog() {
   int i,c;
   String []threeColors= {"Red","Green","Blue"};
   String []dirs={"N__","NE","__E","SE","S__","SW","__W","NW"};
   boolean rf,bf;
   GenericDialog gd = new GenericDialog("Pixel Crosstalk parameters");
   gd.addStringField("Image suffix: ", crossTalkName, 32);
   gd.addCheckbox("Convert in-place?", in_place);
// public static boolean monoRed=false; /* use crosstalk for red in all channels */
// public static boolean monoBlue=false; /* use crosstalk for blue in all channels */
   gd.addCheckbox("Use red filter?",  monoRed && !monoBlue);
   gd.addCheckbox("Use blue filter?", !monoRed && monoBlue);
   for (c=0;c<threeColors.length;c++) for (i=0;i<dirs.length;i++) {
     gd.addNumericField(threeColors[c]+"_"+dirs[i],   crossCoeff[c][i], 4);
   }
   gd.addNumericField("Crosstalk_kernel_FFT_Size:",           FFTSize, 0);
   gd.addCheckbox("Bayer_mask:_keep_G1?", (BayerMask & 1)!=0);
   gd.addCheckbox("Bayer_mask:_keep_R?",  (BayerMask & 2)!=0);
   gd.addCheckbox("Bayer_mask:_keep_B?",  (BayerMask & 4)!=0);
   gd.addCheckbox("Bayer_mask:_keep_G2?", (BayerMask & 8)!=0);
   gd.addNumericField("Crosstalk_measurement FFT_Size:",           imageFFTSize, 0);
   gd.addNumericField("Crosstalk processing threshold:", crosstalkThreshold, 3);

   gd.addNumericField("Debug Level:",      MASTER_DEBUG_LEVEL, 0);



   gd.showDialog();
   if (gd.wasCanceled()) return false;
   crossTalkName = gd.getNextString();
     in_place=gd.getNextBoolean();
     rf=gd.getNextBoolean();
     bf=gd.getNextBoolean();
     monoRed=rf;
     monoBlue=bf && !rf;
     for (c=0;c<threeColors.length;c++) for (i=0;i<dirs.length;i++) {
       crossCoeff[c][i]=gd.getNextNumber();
     }
     FFTSize=1;
     for (i=(int) gd.getNextNumber(); i >1; i>>=1) FFTSize <<=1; /* make FFTSize to be power of 2*/

     BayerMask  = (gd.getNextBoolean())?1:0;
     BayerMask |= (gd.getNextBoolean())?2:0;
     BayerMask |= (gd.getNextBoolean())?4:0;
     BayerMask |= (gd.getNextBoolean())?8:0;
     imageFFTSize=1;
     for (i=(int) gd.getNextNumber(); i >1; i>>=1) imageFFTSize <<=1; /* make FFTSize to be power of 2*/
     crosstalkThreshold= gd.getNextNumber();
     MASTER_DEBUG_LEVEL=  (int) gd.getNextNumber();
     kernels_valid=false;
     return true;
  }
  public void maskBayer(ImageProcessor ip) {
//BayerMask
    float [] pixels=(float[])ip.getPixels();
    int row,col,c;
    int width=ip.getWidth();
    int height=ip.getHeight();
    for (row=0;row<height;row++) {
      for (col=0;col<width; col++) {
        c= (((row & 1) <<1) | (col &1));
        if (((1<<c) & BayerMask) ==0) pixels[row*width+col]=0.0f;
      }
    }
   ip.setPixels(pixels);
  }

  public void convolveBayerKernel(ImageProcessor ip) {
    float [] src_pixels=(float[])ip.getPixels();
    float [] dst_pixels=new float[src_pixels.length];
//arraycopy(Object src, int srcPos, Object dest, int destPos, int length) 
//    System.arraycopy(src_pixels, 0, dst_pixels, 0, src_pixels.length); /* for the borders closer to 1/2 kernel size*/
//    java.util.Arrays.fill(mollyarray,0);
    Arrays.fill(dst_pixels,0.0f);
    int row,col,i,j,c,c0;
    int size=(FFTSize>>1)-1;
    int fullSize=FFTSize-1;
    float d;
    int width=ip.getWidth();
    int height=ip.getHeight();
    c0=(size & 1) | ((size & 1)<<1); //usually 3
    for (row=0;row<(height-fullSize);row++) {
      for (col=0;col<(width-fullSize); col++) {
        c= (((row & 1) <<1) | (col &1)) ^ c0;
        d=src_pixels[(row+size)*width+col+size];
        for (i=0;i<fullSize;i++) for (j=0;j<fullSize;j++) {
          dst_pixels[(row+i)*width+col+j]+=d*kernels[c][i][j];
        }
      }
    }

/**TODO:  Add margins */
   ip.setPixels(dst_pixels);
  }


 public void set2addCrosstalk() {
   int i,j,l;
   kernels=new float [4][FFTSize-1][FFTSize-1];
   for (i=1;i<FFTSize;i++) for (j=1;j<FFTSize;j++) {
      l=i*FFTSize+j;
      kernels[0][i-1][j-1]=direct_kg[l];
      kernels[1][i-1][j-1]=direct_kr[l];
      kernels[2][i-1][j-1]=direct_kb[l];
      kernels[3][i-1][j-1]=direct_kg[l];
   }
 }

 public void set2removeCrosstalk() {
   int i,j,l;
   kernels=new float [4][FFTSize-1][FFTSize-1];
   for (i=1;i<FFTSize;i++) for (j=1;j<FFTSize;j++) {
      l=i*FFTSize+j;
      kernels[0][i-1][j-1]=reverse_kg1[l];
      kernels[1][i-1][j-1]=reverse_kr[l];
      kernels[2][i-1][j-1]=reverse_kb[l];
      kernels[3][i-1][j-1]=reverse_kg2[l];
   }
 }


 public void createReverseKernels() {
//   showDialog();
   ip_kr = new FloatProcessor(FFTSize,FFTSize);
   ip_kg = new FloatProcessor(FFTSize,FFTSize);
   ip_kb = new FloatProcessor(FFTSize,FFTSize);
   IJ.showStatus("Creating crosstalk kernels...");
   initCrosstalkKernels();
   if (DEBUG_LEVEL>1) showCrosstalkKernels();
   IJ.showStatus("Calculating FFHT of crosstalk kernels...");
   FHTCrosstalkKernels();
 //  showCrosstalkKernelsFHT(); ///Not needed - spectra are displayed during transfrorm (depending on FFT options)
   setupFK();
   if (DEBUG_LEVEL>3) testShowAllFHT();
    IJ.showStatus("Reversing crosstalk kernels...");
   crossTalkResolve();
   if (DEBUG_LEVEL>1) showFHT (FFTHalf2FHT(FRslt_g1),  "FRslt_g1");
   if (DEBUG_LEVEL>1) showFHT (FFTHalf2FHT(FRslt_r),   "FRslt_r");
   if (DEBUG_LEVEL>1) showFHT (FFTHalf2FHT(FRslt_b),   "FRslt_b");
   if (DEBUG_LEVEL>1) showFHT (FFTHalf2FHT(FRslt_g2),  "FRslt_g2");
   Rslt_g1=iFHT(FFTHalf2FHT(FRslt_g1),  "Rslt_g1");
   Rslt_r= iFHT(FFTHalf2FHT(FRslt_r),   "Rslt_r");
   Rslt_b= iFHT(FFTHalf2FHT(FRslt_b),   "Rslt_b");
   Rslt_g2=iFHT(FFTHalf2FHT(FRslt_g2),  "Rslt_g2");
   kernels_valid=true;
// private static  float[] reverse_kg,reverse_kr,reverse_kb,reverse_kg1;
   reverse_kg1=(float[])Rslt_g1.getPixels();
   reverse_kr= (float[])Rslt_r.getPixels();
   reverse_kb= (float[])Rslt_b.getPixels();
   reverse_kg2=(float[])Rslt_g2.getPixels();
   IJ.showStatus("Calculating reverse kernels DONE.");

 }





public void crossTalkResolve () {
  int row, col;
  double[][][] MK=new double [4][4][2]; /* elements are initialized from FK** */
  double[][][] RK;
  FRslt_g1=new double[(FFTSize>>1)+1][FFTSize][2];
  FRslt_r= new double[(FFTSize>>1)+1][FFTSize][2];
  FRslt_b= new double[(FFTSize>>1)+1][FFTSize][2];
  FRslt_g2=new double[(FFTSize>>1)+1][FFTSize][2];
  for (row=0;row<=(FFTSize>>1); row++ ) {
    for (col = 0; col<FFTSize; col++) {
      MK[0][0][0]= FKg[row][col][0];   MK[0][1][0]=  FKr[row][col][0];   MK[0][2][0]=  FKb[row][col][0];   MK[0][3][0]=  FKg[row][col][0];
      MK[0][0][1]= FKg[row][col][1];   MK[0][1][1]=  FKr[row][col][1];   MK[0][2][1]=  FKb[row][col][1];   MK[0][3][1]=  FKg[row][col][1];
      MK[1][0][0]= FKgX[row][col][0];  MK[1][1][0]= -FKrX[row][col][0];  MK[1][2][0]=  FKbX[row][col][0];  MK[1][3][0]= -FKgX[row][col][0];
      MK[1][0][1]= FKgX[row][col][1];  MK[1][1][1]= -FKrX[row][col][1];  MK[1][2][1]=  FKbX[row][col][1];  MK[1][3][1]= -FKgX[row][col][1];
      MK[2][0][0]= FKgY[row][col][0];  MK[2][1][0]=  FKrY[row][col][0];  MK[2][2][0]= -FKbY[row][col][0];  MK[2][3][0]= -FKgY[row][col][0];
      MK[2][0][1]= FKgY[row][col][1];  MK[2][1][1]=  FKrY[row][col][1];  MK[2][2][1]= -FKbY[row][col][1];  MK[2][3][1]= -FKgY[row][col][1];
      MK[3][0][0]= FKgXY[row][col][0]; MK[3][1][0]= -FKrXY[row][col][0]; MK[3][2][0]= -FKbXY[row][col][0]; MK[3][3][0]=  FKgXY[row][col][0];
      MK[3][0][1]= FKgXY[row][col][1]; MK[3][1][1]= -FKrXY[row][col][1]; MK[3][2][1]= -FKbXY[row][col][1]; MK[3][3][1]=  FKgXY[row][col][1];
/* RK:= ~(MS*MK*MS)  */
      RK=complexMatrixInvert4(realComplexMatrixMultiply(MS,complexRealMatrixMultiply(MK,MS)));
/* Save to resuls */
/* lost 1/4 somewhere */
      FRslt_g1[row][col][0]=4.0*RK[0][0][0]; FRslt_r[row][col][0]=4.0*RK[0][1][0];  FRslt_b[row][col][0]=4.0*RK[0][2][0]; FRslt_g2[row][col][0]=4.0*RK[0][3][0]; 
      FRslt_g1[row][col][1]=4.0*RK[0][0][1]; FRslt_r[row][col][1]=4.0*RK[0][1][1];  FRslt_b[row][col][1]=4.0*RK[0][2][1]; FRslt_g2[row][col][1]=4.0*RK[0][3][1]; 
    }
  }
}
/*
  | Kg    Kr    Kb    Kg   |
  | KgX  -KrX   KbX  -KgX  | = MK
  | KgY   KrY  -KbY  -KgY  |
  | KgXY -KrXY -KbXY  KgXY |

 private double [][][] FKr,FKrX,FKrY,FKrXY,  FKg,FKgX,FKgY,FKgXY,  FKb,FKbX,FKbY,FKbXY;
 private double [][][] FRslt_g1, FRslt_r, FRslt_b, FRslt_g2; 
 private double [][]    Rslt_g1,  Rslt_r,  Rslt_b,  Rslt_g2; 

MS[4][4]
  private double[][][] complexRealMatrixMultiply   (double[][][] a, double[][] b  ) {
  private double[][][] realComplexMatrixMultiply   (double[][] a, double[][][] b  ) {
  private double[][][] complexMatrixInvert4    (double[][][] a ) {
 private FHT fht_rslt_g1, fht_rslt_r, fht_rslt_b, fht_rslt_g2;

	public FHT getCopy() {
		ImageProcessor ip = super.duplicate();
		FHT fht = new FHT(ip);
		fht.isFrequencyDomain = isFrequencyDomain;
		fht.quadrantSwapNeeded = quadrantSwapNeeded;
		fht.rgb = rgb;
		fht.originalWidth = originalWidth;
		fht.originalHeight = originalHeight;
		fht.originalBitDepth = originalBitDepth;		
		fht.originalColorModel = originalColorModel;		
		return fht;
	}


fht_kr
inverseTransform()

*/
/* make FHT from array of pixels, perform inverse FHT */
//public ImagePlus iFHT(float[] fht_pixels, String title) {
public FHT iFHT(float[] fht_pixels, String title) {
   FHT fht=fht_kr.getCopy(); /* could not find a way to create a frequency-domain FHT, had to copy one */
//   ImageProcessor ip_fht = new FloatProcessor(FFTSize,FFTSize);
/*
   ip_fht.setPixels(fht_pixels);
   ip_fht.resetMinAndMax();
   ImagePlus imp= new ImagePlus(title, ip_fht);
*/
   fht.setPixels(fht_pixels);
   fht.inverseTransform();
   fht.swapQuadrants(); /* put 0,0 in the center */

   float[] pixels=(float[])fht.getPixels();
   int i;
   double S=0.0;
   for (i=0;i<pixels.length;i++) S+=pixels[i];
   fht.resetMinAndMax();
//   ImagePlus imp= new ImagePlus(title, ip_fht);
   if (DEBUG_LEVEL>0) {
     ImagePlus imp= new ImagePlus(title+" "+S, fht);
     imp.show();
   }
//   return imp;
   return fht;
}

public void FHTCrosstalkKernels() {
    fht_kr =  new FHT(ip_kr);
    fht_kg =  new FHT(ip_kg);
    fht_kb =  new FHT(ip_kb);
/* Swapping quadrants, so the center will be 0,0 */

    fht_kr.swapQuadrants();
    fht_kg.swapQuadrants();
    fht_kb.swapQuadrants();

    fht_kr.transform();
    fht_kg.transform();
    fht_kb.transform();
}

private void showFHT (float [] fht_pixels, String title) {
   ImageProcessor ip_fht = new FloatProcessor(FFTSize,FFTSize);
   ip_fht.setPixels(fht_pixels);
   ip_fht.resetMinAndMax();
   ImagePlus imp= new ImagePlus(title, ip_fht);
   imp.show();
}




private void testShowAllFHT () {
   showFHT (FFTHalf2FHT(FKr),  "FKr");
   showFHT (FFTHalf2FHT(FKrX), "FKrX");
   showFHT (FFTHalf2FHT(FKrY), "FKrY");
   showFHT (FFTHalf2FHT(FKrXY),"FKrXY");

   showFHT (FFTHalf2FHT(FKg),  "FKg");
   showFHT (FFTHalf2FHT(FKgX), "FKgX");
   showFHT (FFTHalf2FHT(FKgY), "FKgY");
   showFHT (FFTHalf2FHT(FKgXY),"FKgXY");

   showFHT (FFTHalf2FHT(FKb),  "FKb");
   showFHT (FFTHalf2FHT(FKbX), "FKbX");
   showFHT (FFTHalf2FHT(FKbY), "FKbY");
   showFHT (FFTHalf2FHT(FKbXY),"FKbXY");
}


// private double [][][] FKr,FKrX,FKrY,FKrXY,  FKg,FKgX,FKgY,FKgXY,  FKb,FKbX,FKbY,FKbXY;

private void setupFK () {
  FKr  =FHT2FFTHalf (fht_kr);
  FKrX =shiftX(FKr);
  FKrY =shiftY(FKr);
  FKrXY=shiftY(FKrX);
  FKg  =FHT2FFTHalf (fht_kg);
  FKgX =shiftX(FKg);
  FKgY =shiftY(FKg);
  FKgXY=shiftY(FKgX);
  FKb  =FHT2FFTHalf (fht_kb);
  FKbX =shiftX(FKb);
  FKbY =shiftY(FKb);
  FKbXY=shiftY(FKbX);
}

/* shifts sides horizontally by 1/2 of the range */
private double[][][] shiftX (double[][][] k) {
  double[][][] result = new double[(FFTSize>>1)+1][FFTSize][2];
  int row, col;
   for (row=0;row<=(FFTSize>>1);row++) for (col=0;col<FFTSize;col++) {
     result[row][(col + (FFTSize >> 1)) % FFTSize][0]=k[row][col][0];
     result[row][(col + (FFTSize >> 1)) % FFTSize][1]=k[row][col][1];
   }
  return result;
}

/* shifts sides vertically by 1/2 of the range */
private double[][][] shiftY (double[][][] k) {
  double[][][] result = new double[(FFTSize>>1)+1][FFTSize][2];
  int row1, row2, col;
   for (row1=0;row1<=(FFTSize>>1);row1++) {
     row2=(FFTSize>>1)-row1;
     for (col=0;col< FFTSize;col++) {
       result[row1][col][0]= k[row2][col][0];
       result[row1][col][1]=-k[row2][col][1];
     }
   }
  return result;
}

/* converts FHT results (frequency space) to complex numbers of [FFTSize/2+1][FFTSize] */
private double[][][] FHT2FFTHalf (FHT fht) {
   float[] fht_pixels=(float[])fht.getPixels();
   double[][][] fftHalf=new double[(FFTSize>>1)+1][FFTSize][2];
   int row1,row2,col1,col2;

 //  double dbg1, dbg2,dbg3, dbg4;
   for (row1=0;row1<=(FFTSize>>1);row1++) {
     row2=(FFTSize-row1) %FFTSize;
     for (col1=0;col1<FFTSize;col1++) {
       col2=(FFTSize-col1) %FFTSize;
       fftHalf[row1][col1]=   complex( 0.5*(fht_pixels[row1*FFTSize+col1] + fht_pixels[row2*FFTSize+col2]),
                                       0.5*(fht_pixels[row2*FFTSize+col2] - fht_pixels[row1*FFTSize+col1]));
     }
   }
   return fftHalf;
}

/* converts FFT arrays of complex numbers of [FFTSize/2+1][FFTSize] to FHT arrays */
private float[] FFTHalf2FHT (double [][][] fft) {
   float[] fht_pixels=new float [FFTSize*FFTSize];
   int row1,row2,col1,col2;
   for (row1=0;row1<=(FFTSize>>1);row1++) {
     row2=(FFTSize-row1) %FFTSize;
     for (col1=0;col1 < FFTSize;col1++) {
       col2=(FFTSize-col1) %FFTSize;
/* out of bounds */
       fht_pixels[row1*FFTSize+col1]=(float)(fft[row1][col1][0]+fft[row1][col1][1]);
       fht_pixels[row2*FFTSize+col2]=(float)(fft[row1][col1][0]-fft[row1][col1][1]);
     }
   }
   return fht_pixels;
}


//public static boolean monoRed=false; /* use crosstalk for red in all channels */
// public static boolean monoBlue=false; /* use crosstalk for blue in all channels */


public void initCrosstalkKernels() {
   double [] cScales = new double[3];
   int center=FFTSize>>1;
   int [][]dirsVectors={{ 0,-1}, /* N*/
                        { 1,-1}, /* NE*/ 
                        { 1, 0}, /*  E*/ 
                        { 1, 1}, /* SE*/ 
                        { 0, 1}, /* S*/ 
                        {-1, 1}, /* SW*/ 
                        {-1, 0}, /*  W*/ 
                        {-1,-1}};/* NW*/ 
   int c,i;

   for (c=0;c<3;c++) {
     cScales[c]=1.0;
     for (i=0;i<8;i++) cScales[c]+=crossCoeff[monoRed?0:(monoBlue?2:c)][i];
     cScales[c]=1.0/cScales[c];
   }
   ip_kr.putPixelValue(center,  center  ,cScales[monoRed?0:(monoBlue?2:0)]);  
   ip_kg.putPixelValue(center,  center  ,cScales[monoRed?0:(monoBlue?2:1)]);  
   ip_kb.putPixelValue(center,  center  ,cScales[monoRed?0:2]);  
   for (i=0;i<8;i++) {
     ip_kr.putPixelValue(center+dirsVectors[i][0],  center+dirsVectors[i][1]  ,cScales[0]*crossCoeff[monoRed?0:(monoBlue?2:0)][i]);  
     ip_kg.putPixelValue(center+dirsVectors[i][0],  center+dirsVectors[i][1]  ,cScales[1]*crossCoeff[monoRed?0:(monoBlue?2:1)][i]);  
     ip_kb.putPixelValue(center+dirsVectors[i][0],  center+dirsVectors[i][1]  ,cScales[2]*crossCoeff[monoRed?0:2][i]);  
   }

// private static  double [][] direct_kr,direct_kg,direct_kb,; 
   direct_kr=(float[])ip_kr.getPixels();
   direct_kg=(float[])ip_kg.getPixels();
   direct_kb=(float[])ip_kb.getPixels();
   ip_kr.resetMinAndMax();
   ip_kg.resetMinAndMax();
   ip_kb.resetMinAndMax();

}


public void printComplexMatrix(double[][][] a, String title) {
    int i,j;
    String s=new String();
    for (i=0;i<a.length;i++) {
       for (j=0;j<a[0].length;j++) {
         s+= String.format(" %7.3f(%7.3f)",a[i][j][0],a[i][j][1]);
       }
       s+="\n";
    }
    IJ.showMessage(title+" ("+ a.length+"x"+a[0].length+")",s);
  }
  public double [][][] complexRandom(int h, int w) {
    Random generator = new Random( 123456 );
    int i,j;
    double[][][] result =new double[h][w][2];
    for (i=0;i<h;i++) for (j=0;j<w;j++) {
       result[i][j][0]=2.0*(generator.nextDouble()-0.5);
       result[i][j][1]=2.0*(generator.nextDouble()-0.5);
    }
   return result;
  }

  public double[][][] complexMatrixAdd        (double[][][] a, double[][][] b ) {
   int h=a.length;
   int w=a[0].length;
   int i,j;
   if ((b[0].length!=w) || (b.length!=h)) return null;
   double[][][] result=new double[h][w][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[i][j]=complexAdd(a[i][j], b[i][j]);
   }
   return result;
  }
  public double[][][] complexMatrixScale      (double[][][] a, double k ) {
   int h=a.length;
   int w=a[0].length;
   int i,j;
   double[][][] result=new double[h][w][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[i][j]=complexScale(a[i][j], k);
   }
   return result;
  }
  public double[][][] complexMatrixConjugate  (double[][][] a ) {
   int h=a.length;
   int w=a[0].length;
   int i,j;
   double[][][] result=new double[h][w][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[i][j]=complexConjugate(a[i][j]);
   }
   return result;
  }
  public double[][][] complexMatrixTranspose (double[][][] a ) {
   int h=a.length;
   int w=a[0].length;
   int i,j;
   double[][][] result=new double[w][h][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[j][i]=a[i][j];
   }
   return result;
  }
  
  public double[][][] complexMatrixMultiply   (double[][][] a, double[][][] b  ) {
   int h=a.length;
   int n=a[0].length;
   int w=b[0].length;
   int i,j,k;
   if (b.length!=n) return null;
   double[][][] result=new double[h][w][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[i][j][0]=0.0;
     result[i][j][1]=0.0;
     for (k=0;k<n;k++) result[i][j]=complexAdd(result[i][j],complexMultiply(a[i][k],b[k][j])); /* out of boud */
   }
   return result;
  }

  public double[][][] complexRealMatrixMultiply   (double[][][] a, double[][] b  ) {
   int h=a.length;
   int n=a[0].length;
   int w=b[0].length;
   int i,j,k;
   if (b.length!=n) return null;
   double[][][] result=new double[h][w][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[i][j][0]=0.0;
     result[i][j][1]=0.0;
     for (k=0;k<n;k++) result[i][j]=complexAdd(result[i][j],complexScale(a[i][k],b[k][j]));
   }
   return result;
  }

  public double[][][] realComplexMatrixMultiply   (double[][] a, double[][][] b  ) {
   int h=a.length;
   int n=a[0].length;
   int w=b[0].length;
   int i,j,k;
   if (b.length!=n) return null;
   double[][][] result=new double[h][w][2];
   for (i=0;i<h;i++) for (j=0;j<w;j++) {
     result[i][j][0]=0.0;
     result[i][j][1]=0.0;
     for (k=0;k<n;k++) result[i][j]=complexAdd(result[i][j],complexScale(b[k][j],a[i][k]));
   }
   return result;
  }


  public void matrix4InvertInit() { /* May be extended to different dimensions */
    int row,k,l,l1,s;
    int[] i ={0,0,0,0,0};
    boolean[] t= {true,true,true,true};
    int[] seq= {0,0,0,0};
    if ( invertDimension!=4) {
    invertSeq=     new int [4][4][6][7] ; /* once calculated sequence of elements for 4x4 matrix inversion i,j,i,j,i,j,sign */
    determinantSeq=new int [24][5] ;      /* once calculated sequence of elements for 4x4 matrix determinant: j,j,j,j,sign */
      k=0;
      for (i[1]=0;i[1]<4;i[1]++) {
        t[i[1]]=false;
        for (i[2]=0;i[2]<4;i[2]++) if (t[i[2]]) {
          t[i[2]]=false;
          for (i[3]=0;i[3]<4;i[3]++) if (t[i[3]]) {
            t[i[3]]=false;
            for (i[4]=0;i[4]<4;i[4]++) if (t[i[4]]) { /* now all i1,i2,i3,i4 are different */
              determinantSeq[k][0]=i[1];
              determinantSeq[k][1]=i[2];
              determinantSeq[k][2]=i[3];
              determinantSeq[k][3]=i[4];
              seq[0]=i[1];seq[1]=i[2];seq[2]=i[3];seq[3]=i[4];
              determinantSeq[k][4]=1;
              while ((seq[0]>seq[1]) || (seq[1]>seq[2])  || (seq[2]>seq[3])) {
                if      (seq[0]>seq[1]) {s=seq[0];seq[0]=seq[1];seq[1]=s;}
                else if (seq[1]>seq[2]) {s=seq[1];seq[1]=seq[2];seq[2]=s;}
                else if (seq[2]>seq[3]) {s=seq[2];seq[2]=seq[3];seq[3]=s;}
                determinantSeq[k][4]*=-1;
              }
              k++;
            }
            t[i[3]]=true;
          }
          t[i[2]]=true;
        }
        t[i[1]]=true;
      }
      for (row=0; row<4; row++) for (i[1]=0; i[1]<4; i[1]++) {
        t[i[1]]=false;
        k=0;
        for (i[2]=0;i[2]<4;i[2]++) if (t[i[2]]) {
          t[i[2]]=false;
          for (i[3]=0;i[3]<4;i[3]++) if (t[i[3]]) {
            t[i[3]]=false;
            for (i[4]=0;i[4]<4;i[4]++) if (t[i[4]]) { /* now all i1,i2,i3,i4 are different */
              seq[row]=i[1];
              l1=0;
              for (l=0;l<4;l++) if (l!=row) {
                seq[l]=(l>row)?i[l+1]:i[l+2];
                invertSeq[row][i[1]][k][l1++]=l;
                invertSeq[row][i[1]][k][l1++]=seq[l];
              }
/* calculate sign of the term */
              invertSeq[row][i[1]][k][l1]=1;
              while ((seq[0]>seq[1]) || (seq[1]>seq[2])  || (seq[2]>seq[3])) {
                if      (seq[0]>seq[1]) {s=seq[0];seq[0]=seq[1];seq[1]=s;}
                else if (seq[1]>seq[2]) {s=seq[1];seq[1]=seq[2];seq[2]=s;}
                else if (seq[2]>seq[3]) {s=seq[2];seq[2]=seq[3];seq[3]=s;}
                invertSeq[row][i[1]][k][l1]*=-1;
              }
              k++;/* next term */
            }
            t[i[3]]=true;
          }
          t[i[2]]=true;
        }
        t[i[1]]=true;
      }
      invertDimension=4;
    }
  }

/* to troubleshoot - try real matrices */
  private double[][][] complexMatrixInvert4    (double[][][] a ) {
   if ( invertDimension!=4) matrix4InvertInit();
   if ((a.length!=4) ||(a[0].length!=4) || (a[0][0].length!=2)) {
     IJ.showMessage("ERROR in complexMatrixInvert4()","complexMatrixInvert4() requires complex matrices [4,4]" );
     return null;
   }
   int i,j,k;
   double[][][] result = new double[4][4][2];
   double[] det=new double[2];
   double[] ddet=new double[2];
   double[] d=new double[2];
   double[] d1=new double[2];
//String debug_string= new String();
   det=complex (0.0,0.0);
   for (k=0;k<determinantSeq.length;k++) {
     d=complexMultiply(a[0][determinantSeq[k][0]],a[1][determinantSeq[k][1]]);
     d=complexMultiply(d,a[2][determinantSeq[k][2]]);
     d=complexMultiply(d,a[3][determinantSeq[k][3]]);
     d=complexScale   (d,     determinantSeq[k][4]); //sign
     det=complexAdd(det,d);
//     debug_string+=d[0]+"("+d[1]+")\n";
   }
//  IJ.showMessage("Debug:complexMatrixInvert4()","Determinant="+det[0]+"("+det[1]+")\n"+debug_string);
   if ((det[0]==0.0) && (det[1]==0.0)) det[0]=0.0000001; // to make sure determinant is non-zero
   ddet=complexDivide(complex(1.0,0.0),det);
//  IJ.showMessage("Debug:complexMatrixInvert4()","1/Determinant="+ddet[0]+"("+ddet[1]+")" );
/*   d=complexMultiply(det,ddet);
  IJ.showMessage("Debug:complexMatrixInvert4()","1/Determinant*Determinant="+d[0]+"("+d[1]+")" );*/
   for (i=0;i<4;i++) for (j=0;j<4;j++) {
     d1=complex(0.0,0.0);
     for (k=0;k<invertSeq[i][j].length;k++) {
       d=complexMultiply(a[invertSeq[i][j][k][0]][invertSeq[i][j][k][1]] , a[invertSeq[i][j][k][2]][invertSeq[i][j][k][3]] );
       d=complexMultiply(d , a[invertSeq[i][j][k][4]][invertSeq[i][j][k][5]]);
       if (invertSeq[i][j][k][6] <0)  {d[0]=-d[0];d[1]=-d[1];}
       d1[0]+=d[0];
       d1[1]+=d[1]; 
     }
     result[j][i]=complexMultiply(d1,ddet);
   }
   return result;
  }

  public double[] complex(double re, double im) {
   double[] result = {re,im};
   return result;
  }

  public double[] complexAdd(double[] a, double [] b) {
   double[] result = {a[0]+b[0],a[1]+b[1]};
   return result;
  }
  public double[] complexSubtract(double[] a, double [] b) {
   double[] result = {a[0]-b[0],a[1]-b[1]};
   return result;
  }
  public double[] complexScale(double[] a, double k) {
   double[] result = {a[0]*k, a[1]*k};
   return result;
  }
  public double[] complexConjugate(double[] a) {
   double[] result = {a[0],-a[1]};
   return result;
  }
  public double[] complexMultiply(double[] a, double [] b) {
   double[] result = {a[0]*b[0]-a[1]*b[1],a[1]*b[0]+a[0]*b[1]};
   return result;
  }
  public double[] complexDivide(double[] a, double [] b) {
   double l2=b[0]*b[0] + b[1]*b[1];
   double[] result = {(a[0]*b[0]+a[1]*b[1])/l2,(a[1]*b[0]-a[0]*b[1])/l2};
   return result;
  }

/* Still do not understand - how to open scaled image that is all painted */
public void showCrosstalkKernels() {
   int drawScale=16;

   imp_kr= new ImagePlus("Red Crosstalk Kernel", ip_kr);
   imp_kr.show();
   ImageWindow win_kr = imp_kr.getWindow();
   win_kr.getCanvas().setMagnification(drawScale);
   ImageCanvas ic_kr = imp_kr.getCanvas();
   ic_kr.setSourceRect(new Rectangle(0, 0, drawScale*FFTSize, drawScale*FFTSize));
   ic_kr.setDrawingSize(drawScale*FFTSize, drawScale*FFTSize);
   win_kr.pack();
   win_kr.repaint();



 //  imp_kr.updateAndDraw();
   imp_kg= new ImagePlus("Green Crosstalk Kernel", ip_kg);
   imp_kg.show();
   ImageWindow win_kg = imp_kg.getWindow();
   win_kg.getCanvas().setMagnification(drawScale);
   ImageCanvas ic_kg = imp_kg.getCanvas();
   ic_kg.setSourceRect(new Rectangle(0, 0, drawScale*FFTSize, drawScale*FFTSize));
   ic_kg.setDrawingSize(drawScale*FFTSize, drawScale*FFTSize);
//imp_kg.updateAndDraw();
   win_kg.pack();
   win_kg.repaint();


   imp_kb= new ImagePlus("Blue Crosstalk Kernel", ip_kb);
   imp_kb.show();

   ImageWindow win_kb = imp_kb.getWindow();
   win_kb.getCanvas().setMagnification(drawScale);
   ImageCanvas ic_kb = imp_kb.getCanvas();
   ic_kb.setSourceRect(new Rectangle(0, 0, drawScale*FFTSize, drawScale*FFTSize));
   ic_kb.setDrawingSize(drawScale*FFTSize, drawScale*FFTSize);
   win_kb.pack();
//   win_kb.repaint();
//   win_kb.updateImage(imp_kb); /**makes it magnification==1 *
   imp_kb.updateAndRepaintWindow();
   imp_kr.updateAndRepaintWindow();

  }
public void showCrosstalkKernelsFHT() {
   fht_kr.getPowerSpectrum ();
   fht_kg.getPowerSpectrum ();
   fht_kb.getPowerSpectrum ();
  }
/**
 * Main method for debugging.
 *
 * For debugging, it is convenient to have a method that starts ImageJ, loads an
 * image and calls the plugin, e.g. after setting breakpoints.
 * Grabbed from https://github.com/imagej/minimal-ij1-plugin
 * @param args unused
 */
public static void main(String[] args) {
	// set the plugins.dir property to make the plugin appear in the Plugins menu
	Class<?> clazz = Aberration_Calibration.class;
	String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
	String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
	System.setProperty("plugins.dir", pluginsDir);
	// start ImageJ
	new ImageJ();
	// run the plugin
	IJ.runPlugIn(clazz.getName(), "");
}

}
/**
.
Interpixel crosstalk causes some "leak" between pixels, leak that happens after the color filter array.
And photons of different wavelengths get to the different depth, so red pixels cause more crosstalk than
the green or blue ones. This program calculates the convolution kernels to correct that effect.
To use those kernels you need to split acquired pixel (Bayer) array:
G1 R  G1 R  G1 R..
B  G2 B  G2 B  G2..
G1 R  G1 R  G1 R..
B  G2 B  G2 B  G2..
...
into 4 sub-arrays:
G1 0  G1 0  G1 0..
0  0  0  0  0  0
G1 0  G1 0  G1 0..
0  0  0  0  0  0

0  R  0  R  0  R.
0  0  0  0  0  0
0  R  0  R  0  R.
0  0  0  0  0  0

0  0  0  0  0  0
B  0  B  0  0  G2..
0  0  0  0  0  0
B  0  B  0  0  G2..

and
0  0  0  0  0  0
0  G2 0  G2 0  G2..
0  0  0  0  0  0
0  G2 0  G2 0  G2..

Then convolve each of them with the corresponding reverse kernel calculated by this program and add together.
The result should be image with eth crosstalk compensated.

Below is the explanation to the reverse kernels calculation.

Bayer arrays (repeat through all the image)
Bg1=1 0
    0 0 

Br= 0 1
    0 0 

Bb= 0 0
    1 0 

Bg2=0 0
    0 1

Inp (x,y) - input light (before crosstalk)
Meas(x,y) - measured pixels (suffered by crosstalk)

Kr, Kg, Kb - forward crosstalk kernels 

@ - convolution

Meas=(Inp*Br) @ Kr + (Inp*(Bg1+Bg2)) @ Kg + (Inp*Bb) @ Kb

In the frequency domain (after applying forward 2d Fourier transform F[])
F[Meas]= (F[Inp]@F[Br])  * F[Kr] + 
         (F[Inp]@F[Bg1]) * F[Kg] +
         (F[Inp]@F[Bg2]) * F[Kg] +
         (F[Inp]@F[Bb])  * F[Kb]
So multiplication (masking) became convolution, convolution - multiplication. Luckily convolution with F[Br] (and other Bayer masks)
is rather simple. Convolutions result that the source is replicated 4 times, each replica is shifted by half of the full range
vertically and/or horizontally (i.e. for 16x16 FFT the shift is by 8 pixels), the replicas are combined with coeffients of +/- 1

FKr=F[Kr], FKg=F[Kg], FKb=F[Kb],...
FX  = F[Inp] shifted by 1/2 X range (horizontal quadrant swap: 12/34 -> 21/43)
FY  = F[Inp] shifted by 1/2 Y range (vertical quadrant swap: 12/34 -> 34/12)
FXY = F[Inp] shifted by 1/2 X and 1/2 Y range (diagonal quadrant swap: 12/34 -> 43/21)

Here are the signs of the replicas:

F[Inp]@F[Br] =F-FX+FY-FXY
F[Inp]@F[Bg1]=F+FX+FY+FXY
F[Inp]@F[Bg2]=F-FX-FY+FXY
F[Inp]@F[Bb] =F+FX-FY-FXY

This can be reprfesented with matrices:

  | F[Inp]@F[Bg1] |     | 1  1  1  1 |     | F  |
  | F[Inp]@F[Br]  |  =  | 1 -1  1 -1 | *   | FX |
  | F[Inp]@F[Bb]  |     | 1  1 -1 -1 |     | FY |
  | F[Inp]@F[Bg2] |     | 1 -1 -1  1 |     | FXY|

  | F[Inp]@F[Bg1] |
  | F[Inp]@F[Br]  |  =  VFIC
  | F[Inp]@F[Bb]  |
  | F[Inp]@F[Bg2] |

  | 1  1  1  1 |
  | 1 -1  1 -1 |     = MS
  | 1  1 -1 -1 |
  | 1 -1 -1  1 |

  | F  |
  | FX |  = VF
  | FY |
  | FXY|

  VFIC = MS *VF

F[Meas]= (F+FX+FY+FXY) * FKg +
         (F-FX+FY-FXY) * FKr + 
         (F+FX-FY-FXY) * FKb +
         (F-FX-FY+FXY) * FKg

  | Kg |
  | Kr | = VK
  | Kb |
  | Kg |

FM =|Kg Kr Kb Kg| * VFIC
where FM =F[Meas]

FM= transpose(VK) * VFIC = transpose(VK) * MS * VF

Applying same quadrant swap to FM (measured data) we can get 3 more equations (to be used to represent the
result as a convolution with 4 masked sub-arrays of the image)

FM =F[Meas], FMX=..., FMY=..., FMXY=... (quadrant swap)

  | FMBg1 |    | F[Meas*Bg1] |    | F[Meas]@F[Bg1] |    | FM+FMX+FMY+FMXY  |    |  1  1  1  1  |    | FM   |
  | FMBr  | =  | F[Meas*Br]  | =  | F[Meas]@F[Br]  |  = | FM-FMX+FMY-FMXY  | =  |  1 -1  1 -1  | *  | FMX  |
  | FMBb  |    | F[Meas*Bb]  |    | F[Meas]@F[Br]  |    | FM+FMX-FMY-FMXY  |    |  1  1 -1 -1  |    | FMY  |
  | FMBg2 |    | F[Meas*Bg2] |    | F[Meas]@F[Bg2] |    | FM-FMX-FMY+FMXY  |    |  1 -1 -1  1  |    | FMXY |

  | FMBg1 |
  | FMBr  | =  VFMB // vector of F[] of the pixels multiplied (masked) by Bayer (i.e. FMBg1 - Fourier of all but G1 pixels set to zero)
  | FMBb  |
  | FMBg2 |

  | FM   |
  | FMX  | =  VFM // Vector of Fourier of measured pixels and with quadrants swapped
  | FMY  |
  | FMXY |

VFMB= MS* VFM



================

FM= transpose(VK) * VFIC = transpose(VK) * MS * VF


VFX=MX*VF
  | 0  1  0  0 |
  | 1  0  0  0 | = MX
  | 0  0  0  1 |
  | 0  0  1  0 |

  | 0  0  1  0 |
  | 0  0  0  1 | = MY
  | 1  0  0  0 |
  | 0  1  0  0 |

  | 0  0  0  1 |
  | 0  0  1  0 | = MXY
  | 0  1  0  0 |
  | 1  0  0  0 |


FM=   transp(VK)   * MS *       VF
FMX=  transp(VKX)  * MS * MX *  VF
FMY=  transp(VKY)  * MS * MY *  VF
FMXY= transp(VKXY) * MS * MXY * VF

FMX=  transp(VKX)  * (MS * MX  * ~MS) * MS * VF
FMY=  transp(VKY)  * (MS * MY  * ~MS) * MS * VF
FMXY= transp(VKXY) * (MS * MXY * ~MS) * MS * VF

MSX=  (MS * MX  * ~MS)
MSY=  (MS * MY  * ~MS)
MSXY= (MS * MXY * ~MS)

  | 1  0  0  0 |
  | 0 -1  0  0 | = MSX=  (MS * MX  * ~MS)
  | 0  0  1  0 |
  | 0  0  0 -1 |

  | 1  0  0  0 |
  | 0  1  0  0 | = MSY=  (MS * MY  * ~MS)
  | 0  0 -1  0 |
  | 0  0  0 -1 |

  | 1  0  0  0 |
  | 0 -1  0  0 | = MSXY= (MS * MXY  * ~MS)
  | 0  0 -1  0 |
  | 0  0  0  1 |

FM=   transp(VK)   *        MS * VF
FMX=  transp(VKX)  * MSX  * MS * VF
FMY=  transp(VKY)  * MSY  * MS * VF
FMXY= transp(VKXY) * MSXY * MS * VF

  | Kg |
  | Kr | = VK
  | Kb |
  | Kg |

  | KgX  |
  | KrX  | = VKX
  | KbX  |
  | KgX  |

  | KgY  |
  | KrY  | = VKY
  | KbY  |
  | KgY  |

  | KgXY |
  | KrXY | = VKXY
  | KbXY |
  | KgXY |


  | Kg    Kr    Kb    Kg   |
  | KgX  -KrX   KbX  -KgX  | = MK
  | KgY   KrY  -KbY  -KgY  |
  | KgXY -KrXY -KbXY  KgXY |

  | FM   |
  | FMX  | =  VFM // Vector of Fourier of measured pixels and with quadrants swapped
  | FMY  |
  | FMXY |

  | FMBg1 |
  | FMBr  | =  VFMB // vector of F[] of the pixels multiplied by bayer (i.e. FMBg1 - Fourier of all but G1 pixels set to zero)
  | FMBb  |
  | FMBg2 |

VFM=MK*MS*VF  //! Fourie of the measured picture (differently swapped quadrants) from "ideal" picture and kernels
VFMB= MS* VFM //! Fourie of the measured picture masked by Bayer patterns from Fourie of the (not masked) measured picture, shifted by quadrants

VFMB= MS*MK*MS* VF
VF =~(MS*MK*MS) * VFMB

MI=~(MS*MK*MS) //!, then

VF =MI * VFMB //! We need just the first line of MI, make inverse FFT and it will produce the 4 kernel for convolution
              //! with measured picture, masked with Bayer patterns (result is the sum of them)



*/
