package com.elphel.imagej.dp;
/**
** -----------------------------------------------------------------------------**
** MTF_Bayer.java
**
** Calculates per-color, MTF for raw Bayer color images
** Based on original SE_MTF plugin
**
** Copyright (C) 2010 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  MTF_Bayer.java is free software: you can redistribute it and/or modify
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
import ij.text.*;


//public class MTF_Bayer extends PlugInFrame implements ActionListener, DialogListener {
public class MTF_Bayer extends PlugInFrame implements ActionListener {
  /**
	 * 
	 */
	private static final long serialVersionUID = 2575167919068149071L;
Panel panel1,panel2;
  int previousID;
  static Frame instance;
  Plot plotResult;

 private static int DEBUG_LEVEL = 1;
 private static int MASTER_DEBUG_LEVEL = 1;
 private static int sSize=32; /* change name to something meaningful */
 private static boolean IS_MONOCHROMATIC=true; /* in monochromatic mode chromatic aberrations are eliminated, edge position is common for all Bayer components */
 private static int EMODE=2; /* trying different modes to determine the edge position on individual lines. 2 seems to be the best so far */
 private static boolean TRUNCATE_ROWS=true; /* truncate number of processed rows to equlaize number of phase rotations */
 private static double LSFHamming=0.5; /**Hamming window width for LSF calculation, fraction of the full window. 0.0 - no window */
 private static boolean USE_COMPLEX=true;
 private static boolean SHOW_edgeinfo=false;
 private static boolean [] SHOW_ESF={false,true,true,false};
 private static boolean [] SHOW_LSF={false,true,true,false};
 private static boolean [] SHOW_MTF={true,true,false,false,false,false,false};
 private static boolean [] SHOW_GGRB={false,true,true,true,true,true,true};
 private static boolean SHOW_CROSSTALK=false;

 private static boolean ALLOW_FULL_IMAGE=false;
 private static double GB2GR_CORR=1.0; /* difference in Green sensitivity Gb/Gr*/
 private int maxJump=8;
 private ImagePlus imp_src;
 private boolean horizontalEdges=false;
 String title_src;
 int  selecWidth,selecHeight;
 private double [][]   ESFArray;
 private int [][]       posEdges;
 private double[][][]   color_coeff_abc;
 private double[][][] binDataBoth;
 private double[][][] binDataDiffBoth;
 private double[][]   MTFVectorsBoth =    new double [8][];   /* Amplitudes */
 private double[][]   MTFPVectorsBoth =   new double [8][];  /**Phases */
 private double[][]   MTFCVectorsBoth =   new double [16][]; /**Complex */
 private Complex[][]  MTFComplex=         new Complex [8][];
 private double[][]   MTFCRatios =        new double [16][]; /**Complex, each MTF divided  by maximal (separately for each direction) */
 private double[][]   MTFCRatiosA =       new double [8][];  /* Amplitude */
 private double[][]   MTFCRatiosRe=       new double [8][];  /* Real part */
 private double[][]   MTFCRatiosAMasked=  new double [8][];  /* Amplitude, plots for main color zeroed */
 private double[][]   MTFCRatiosReMasked= new double [8][];  /* Real part, Masked */
 private double[][]   MTFCRatiosP =       new double [8][];  /* Phase */
 private Complex[][]  MTFRatioComplex=    new Complex [8][];
 private Complex[][]  MTF_GG_RB=          new Complex [2][];
 private double [][]  MTF_GG_RB_REIM    = new double  [4][];

 private int          MTFMaximalChannel;                   /* Color number that has highest DC value - averaged for both directions */
 private double[]     MTF_DC          =   new double [8]; /* DC coefficinets of MTF - calculated separately from ESF to avoid window function skew */



 private int [][][] bayer_loc={{{0,0}, {0,1},{1,0},{1,1}},
                              {{0,1}, {1,1},{0,0},{1,0}}};
 Color [] bayer_colors=  { Color.green, Color.red, Color.blue, Color.cyan, Color.black };
 Color [] bayerColorsTwice= { Color.green, Color.green, Color.red, Color.red,
                              Color.blue,  Color.blue,  Color.cyan, Color.cyan, Color.black, Color.black };
 Color [] redGreen_colors=    { Color.red, Color.green, Color.black };
 Color [] redGreenColorsTwice= { Color.red, Color.green, Color.red, Color.green, Color.black  };
 Color [] redBlue_colors=      { Color.red, Color.blue, Color.black };
 Color [] redBlueColorsTwice=  { Color.red, Color.blue, Color.red, Color.blue, Color.black };

 Color [] two_colors=  { Color.black, Color.blue};
 String colorCSVHeaders=    "\tGreen (R)\tRed\tBlue\tGreen (B)";
 String colorCSVHeadersBoth="\tGreen (R) _/\tGreen (R) \\_\tRed _/\tRed \\_\tBlue _/\tBlue \\_\tGreen (B) _/\tGreen (B) \\_";
 String colorCSVHeadersBothComplex="\tRe(Gr _/)\tIm(Gr _/)\tRe(Gr \\_)\tIm(Gr \\_)\tRe(R _/)\tIm(R _/)\tRe(R \\_)\tIm(R \\_)\tRe(B _/)\tIm(B _/)\tRe(B \\_)\tIm(B \\_)\tRe(Gb _/)\tIm(Gb _/)\tRe(Gb \\_)\tIm(Gb \\_)";
// no yMax
 String CSVHeadersBoth="\t _/\t \\_";
 String CSVHeadersBothComplex="\tRe(_/)\tIm(_/)\tRe(\\_)\tIm(\\_)";


 public MTF_Bayer() {
    super("OTF for Bayer mosaic images");
    if (IJ.versionLessThan("1.39t")) return;
    if (instance!=null) {
      instance.toFront();
      return;
    }
    instance = this;
    addKeyListener(IJ.getInstance());

    setLayout(new FlowLayout());
    panel1 = new Panel();
    panel1.setLayout(new GridLayout(4, 1, 50, 5));
    addButton("Configure",panel1);
    addButton("Process",panel1);
    addButton("Calculate Crosstalk",panel1);
    addButton("OTF Ratio",panel1);
    add(panel1);
    panel2 = new Panel();
    panel2.setLayout(new GridLayout(4, 3, 5, 5));
    addButton("Show ESF",panel2);
    addButton("Show ESF Normalized",panel2);
    addButton("Show ESF Text",panel2);
    addButton("Show LSF",panel2);
    addButton("Show LSF Normalized",panel2);
    addButton("Show LSF Text",panel2);
    addButton("Show MTF",panel2);
    addButton("Show MTF Text",panel2);
    addButton("Show OTF Text",panel2);
    addButton("Show Approximation",panel2);
    addButton("Help",panel2);
    add(panel2);
    pack();
    GUI.center(this);
    setVisible(true);
 }


  void addButton(String label, Panel panel) {
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
    }
    DEBUG_LEVEL=MASTER_DEBUG_LEVEL; /* Move to particular command if needed */
    if (label.equals("Process")) {
      if (!openImage()) {
//        IJ.showMessage("Error","openImage() failed");
        return;
      }
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      double [][] esfV,esfH;
      int [][] posEdgesV=null;
      int [][] posEdgesH=null;

/* try both vertical and horizontal edges, use one that has more edges (if any) */
      esfV = generateSelection(imp_src,  false);
      if (esfV!=null) posEdgesV=findEdges(esfV, maxJump);

      esfH = generateSelection(imp_src,  true);
      if (esfH!=null) posEdgesH=findEdges(esfH, maxJump);
/* Both directions OK - use the best one */
      if ((posEdgesV!=null) && (posEdgesH!=null)) {
        if (posEdgesH[0].length>posEdgesV[0].length) {
          ESFArray=esfH;
          posEdges=posEdgesH;
          esfV=null;
          posEdgesV=null;
          horizontalEdges=true;
        } else {
          ESFArray=esfV;
          posEdges=posEdgesV;
          esfH=null;
          posEdgesH=null;
          horizontalEdges=false;
        }
      } else if (posEdgesH!=null) {
          ESFArray=esfH;
          posEdges=posEdgesH;
          horizontalEdges=true;
      } else if (posEdgesV!=null) {

          ESFArray=esfV;
          posEdges=posEdgesV;
          horizontalEdges=false;
      } else {
        IJ.showMessage("Error","No edges that run side-to-side found.");
        return;
      }
      if (DEBUG_LEVEL>4)  IJ.showMessage("Debug","Using "+(horizontalEdges?"horizontal":"vertical")+" edges, numer of edges="+posEdges[0].length);
//Restore selecWidth,selecHeight
    selecWidth=ESFArray[0].length;
    selecHeight=ESFArray.length;

      color_coeff_abc=approximateEdges(posEdges);
      if (SHOW_edgeinfo) showEdgesCoeff(color_coeff_abc, (ESFArray.length)>>1);
      binESF();
      if (SHOW_ESF[0]) showESF(SHOW_ESF[1],SHOW_ESF[2],SHOW_ESF[3], false,false) ;
      binDataDiffBoth=diffColorVector4Both(binDataBoth, LSFHamming);
      if (SHOW_LSF[0]) showLSF(SHOW_LSF[1],SHOW_LSF[2],SHOW_LSF[3], false,false) ;
      makeMTF();
      if (SHOW_MTF[0]) showMTF(SHOW_MTF[1],SHOW_MTF[2],SHOW_MTF[3],SHOW_MTF[4],SHOW_MTF[5],SHOW_MTF[6]);
      if (SHOW_CROSSTALK) crosstalk();
      return;
    } else if (label.equals("Show ESF"))             showESF(true,  false, false,  false, false) ;
    else   if (label.equals("Show ESF Normalized"))  showESF(false, true,  false,  false, false) ;
    else   if (label.equals("Show ESF Text"))        showESF(false, false, true,   false, false) ;
    else   if (label.equals("Show LSF"))             showLSF(true,  false, false,  false, false) ;
    else   if (label.equals("Show LSF Normalized"))  showLSF(false, true,  false,  false, false) ;
    else   if (label.equals("Show LSF Text"))        showLSF(false, false, true,   false, false) ;
    else   if (label.equals("Show MTF"))             showMTF(true,  SHOW_MTF[2], SHOW_MTF[3], SHOW_MTF[4],false, false) ;
    else   if (label.equals("Show MTF Text"))        showMTF(false, false, false, false, true,  false) ;
    else   if (label.equals("Show OTF Text"))        showMTF(false, false, false, false, false, true) ;
    else   if (label.equals("Show Approximation")) {
      if ((color_coeff_abc==null) || (ESFArray==null)) {
        IJ.showMessage("Error","No data available, please run \"Process\" first");
        return;
      }
      showEdgesCoeff(color_coeff_abc, (ESFArray.length)>>1);
    }    else   if (label.equals("Calculate Crosstalk")) {
        crosstalk();
    } else  if (label.equals("OTF Ratio")) {
       showMTFRatios(true,true,true,true);
    } else  if (label.equals("Help")) {
       showOTFHelp();
    }
  }

  public void crosstalk() {
    String [] Bayer={"green","red","blue","green"};
    int i,qf;
    qf= sSize/4;
    Complex [] quaterFreq =new Complex[8];
    double[] av0= new double [8];
    double[] av4= new double [8];
    double[] aav0= new double [4];
    Complex quaterFreqMC; /* average of two directions of main component */
    if (MTFVectorsBoth[0]==null) {
        IJ.showMessage("Error","No data available, please run \"Process\" first");
        return;
    }
   for (i=0;i<8;i++) {
      quaterFreq[i]=  new Complex(MTFCVectorsBoth[i<<1][qf],MTFCVectorsBoth[(i<<1)+1][qf]);
      av0[i]=MTFCVectorsBoth[i<<1][0]; //always real
      av4[i]=quaterFreq[i].abs();
   }
   quaterFreqMC=new Complex(quaterFreq[MTFMaximalChannel<<1]);
   quaterFreqMC=quaterFreqMC.plus(quaterFreq[(MTFMaximalChannel<<1)+1]);
   quaterFreqMC=quaterFreqMC.scale(1.0/(av0[MTFMaximalChannel<<1]+av0[(MTFMaximalChannel<<1)+1]));



   for (i=0;i<4;i++)  aav0[i]=0.5*(av0[i<<1]+av0[(i<<1)+1]);
   Complex GBmGR1 = quaterFreq[6].minus(quaterFreq[0]);
   Complex GBmGR2 = quaterFreq[7].minus(quaterFreq[1]);
   Complex RmB1 =   quaterFreq[2].minus(quaterFreq[4]);
   Complex RmB2 =   quaterFreq[3].minus(quaterFreq[5]);
   if ((RmB1.abs()==0.0) || (RmB2.abs()==0.0)) {
      IJ.showMessage("Error","Divide by 0");
      return;
   }
   Complex Kv1 =   GBmGR1.divideBy(RmB1);
   Complex Kv2 =   GBmGR2.divideBy(RmB2);
   Kv1=Kv1.scale(horizontalEdges?(-0.5):0.5);   
   Kv2=Kv2.scale(horizontalEdges?(-0.5):0.5);
   Complex Kv=Kv1.plus(Kv2);
   Kv=Kv.scale(0.5);
   double KHV, Khv, Khv1, Khv2; /**Kh-Kv, considering Kd==0, Kv-Kh= (GR(0)-GB(0))/(R(0)-B(0))/2/(K0-4*Kd)  */
/* next direction is always absolute, does not depend on horizontalEdges */
   Khv1=0.5*(av0[0]-av0[6])/(av0[2]-av0[4]); 
   Khv2=0.5*(av0[1]-av0[7])/(av0[3]-av0[5]);
   Khv=0.5*(Khv1+Khv2);
   KHV=0.5*((MTF_DC[0]+MTF_DC[1])-(MTF_DC[6]+MTF_DC[7]))/((MTF_DC[2]+MTF_DC[3])-(MTF_DC[4]+MTF_DC[5]));
//(GB(0)-GR(0))/(R(0)-B(0))=2*(Kv-Kh)/(K0-4*kd)
//so KV/K0 = (Gb(1/4)-Gr(1/4))/(R(1/4)-B(1/4))/2
// private double[][]   MTFCVectorsBoth =   new double [16][]; /**Complex */
  int precision=3;
  new TextWindow(title_src+"_"+"pixel crosstalk", "Parameter                  \tValue",
  ((DEBUG_LEVEL>0)?("File\t"+title_src+"\n"):"")+
  "Edge mode\t"+(IS_MONOCHROMATIC?"monochrome":"color")+"\n"+
  ((GB2GR_CORR!=1.0)?("Gb/Gr corr.\t"+IJ.d2s(GB2GR_CORR,precision+2)+"\n"):"")+
  ((DEBUG_LEVEL>2)?("Red[0]\t"+     IJ.d2s(aav0[1],precision)+"\n"):"")+
  ((DEBUG_LEVEL>2)?("GreenR[0]\t"+  IJ.d2s(aav0[0],precision)+"\n"):"")+
  ((DEBUG_LEVEL>2)?("GreenB[0]\t"+  IJ.d2s(aav0[3],precision)+"\n"):"")+
  ((DEBUG_LEVEL>2)?("Blue[0]\t"+    IJ.d2s(aav0[2],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_GR _/\t"+  IJ.d2s(MTF_DC[0],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_GR \\_\t"+  IJ.d2s(MTF_DC[1],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_R  _/\t"+  IJ.d2s(MTF_DC[2],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_R  \\_\t"+  IJ.d2s(MTF_DC[3],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_B  _/\t"+  IJ.d2s(MTF_DC[4],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_B  \\_\t"+  IJ.d2s(MTF_DC[5],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_GB _/\t"+  IJ.d2s(MTF_DC[6],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("MTF_DC_GB \\_\t"+  IJ.d2s(MTF_DC[7],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("R["+qf+"]  _/\t"+ stringComplexComplex(quaterFreq[2],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("R["+qf+"] \\_\t"+ stringComplexComplex(quaterFreq[3],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("GR["+qf+"]  _/\t"+stringComplexComplex(quaterFreq[0],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("GR["+qf+"] \\_\t"+stringComplexComplex(quaterFreq[1],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("GB["+qf+"]  _/\t"+stringComplexComplex(quaterFreq[6],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("GB["+qf+"] \\_\t"+stringComplexComplex(quaterFreq[7],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("B["+qf+"]  _/\t"+ stringComplexComplex(quaterFreq[4],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("B["+qf+"] \\_\t"+ stringComplexComplex(quaterFreq[5],precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("GB-GR  _/\t"+     stringComplexComplex(GBmGR1,precision)+"\n"):"")+ 
  ((DEBUG_LEVEL>1)?("GB-GR \\_\t"+     stringComplexComplex(GBmGR2,precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("R-B  _/\t"+       stringComplexComplex(RmB1,precision)+"\n"):"")+ 
  ((DEBUG_LEVEL>1)?("R-B  \\_\t"+       stringComplexComplex(RmB2,precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("Kv  _/\t"+        stringComplexComplex(Kv1,precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("Kv  \\_\t"+        stringComplexComplex(Kv2,precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("Kv average\t"+         stringComplexComplex(Kv,precision)+"\n"):"")+
  ((DEBUG_LEVEL>1)?("dG["+qf+"]/G[0] _/\t"+IJ.d2s(GBmGR1.abs()/aav0[0],precision)+"\n"):"")+ /* quality factor */
  ((DEBUG_LEVEL>1)?("dG["+qf+"]/G[0] \\_\t"+IJ.d2s(GBmGR2.abs()/aav0[0],precision)+"\n"):"")+ /* quality factor */
  ((DEBUG_LEVEL>2)?("Kh-Kv _/1\t"+    IJ.d2s(Khv1,precision+1)+"\n"):"")+
  ((DEBUG_LEVEL>2)?("Kh-Kv \\_\t"+    IJ.d2s(Khv1,precision+1)+"\n"):"")+
  "Filter\t"+     (((aav0[1] > aav0[0]) && (aav0[1] > aav0[2]) && (aav0[1]>aav0[3]))?"Red":(((aav0[2] > aav0[0]) && (aav0[2] > aav0[1]) && (aav0[2]>aav0[3]))? "Blue":"Green/None"))+"\n"+
  ((DEBUG_LEVEL>0)?(Bayer[MTFMaximalChannel]+"[1/4].abs\t"+IJ.d2s(quaterFreqMC.abs(),precision) +"\n"):"")+
  ((DEBUG_LEVEL>0)?(Bayer[MTFMaximalChannel]+"[1/4].phase\t"+IJ.d2s(quaterFreqMC.phase(),precision)+"\n"):"")+
//  (horizontalEdges?"Horizontal (E+W)/2":"Vertical (S+N)/2")+" crosstalk\t"+(USE_COMPLEX? stringComplexComplex(Kv,precision+1): IJ.d2s(KV,precision+1))+"\n"+
  (horizontalEdges?"Horizontal (E+W)/2":"Vertical (N+S)/2")+" crosstalk\t"+IJ.d2s(Kv.Re(),precision+1)+"\n"+
  ((DEBUG_LEVEL>2)?("Diff. H-V old\t"+IJ.d2s(Khv,precision+1)+"\n"):"")+
  "Diff. H-V\t"    +IJ.d2s(KHV,precision+1)+"\n"+
//  (horizontalEdges?"Vertical (S+N)/2":"Horizontal (E+W)/2")+" crosstalk\t" + (USE_COMPLEX? stringComplexComplex(horizontalEdges?Kv.minus(new Complex(KHV,0)):Kv.plus(new Complex(KHV,0)) ,precision+1):IJ.d2s(KV +(horizontalEdges?-1:1)*KHV,precision+1))+"\n",
  (horizontalEdges?"Vertical (N+S)/2":"Horizontal (E+W)/2")+" crosstalk\t" + IJ.d2s(Kv.Re() +(horizontalEdges?-1:1)*KHV,precision+1)+"\n"+
  (horizontalEdges?"Vertical (N-S)/2":"Horizontal (E-W)/2")+" asymmetry\t" + IJ.d2s(Kv.Im() ,precision+1)+"\n",
  500, (DEBUG_LEVEL>2)?925:((DEBUG_LEVEL>1)?825:325));
  if (SHOW_GGRB[0]) {
    showESF(false,  false, false,  SHOW_GGRB[1], SHOW_GGRB[2]) ;
    showLSF(false,  false, false,  SHOW_GGRB[3], SHOW_GGRB[4]) ;
    showMTFGGRB(SHOW_GGRB[5], SHOW_GGRB[6]);
  }
}
//
  public String stringComplexDouble(double re, double im, int precision) {
    return IJ.d2s(re,precision)+((im>=0)?"+":"")+IJ.d2s(im,precision)+"i";
  }
  public String stringComplexComplex(Complex c, int precision) {
    return stringComplexDouble(c.Re(), c.Im(), precision);
  }
  public void makeMTF() {
    int n,i, n1;
    int N=sSize*4;
    int M=(N>>3)+1; /* Only quarter (because of bins) of the first half is needed */
    Complex[] ArrayComplex = new Complex[N];
    Complex[] VectorFFTC = new Complex[N];
    boolean masked=false;
    for  (n=0;n<8;n++) {
      MTFVectorsBoth[n]=     new double[M];
      MTFPVectorsBoth[n]=    new double[M];
      MTFCVectorsBoth[2*n]=  new double[M];
      MTFCVectorsBoth[2*n+1]=new double[M];
      MTFComplex[n]=new Complex[M];
      for (i=0;i<N;i++) {
         ArrayComplex[(i+(N>>1))%N] = new Complex(binDataDiffBoth[n & 1][n>>1][i], 0); 
      }
      /* Calculate DC coefficients of MTF using ESF (not LSF multiplied by window function. Precise DC coefficients are needed to determine crosstalk
          assymmentry Horizontal-Vertical (may be caused by limited bandwidth of an ammplifier/ADC).
          Use difference between average ESF in the first and the last quaters */
      MTF_DC[n]=0.0;
      for (i=0;i<sSize; i++) {
        MTF_DC[n]+=binDataBoth[n & 1][n>>1][N-i-1];
        MTF_DC[n]-=binDataBoth[n & 1][n>>1][i];
      }
      MTF_DC[n]/=sSize;
      if ((n & 1) !=0) MTF_DC[n]=-MTF_DC[n];
      VectorFFTC = fft(ArrayComplex); 
      for (i = 0; i < M; i++) {
        MTFVectorsBoth[n][i]=VectorFFTC[i].abs();
        MTFPVectorsBoth[n][i]=VectorFFTC[i].phase();
        MTFCVectorsBoth[2*n][i]=  VectorFFTC[i].Re(); 
        MTFCVectorsBoth[2*n+1][i]=VectorFFTC[i].Im();
        MTFComplex[n][i]=new Complex(VectorFFTC[i]);
      }
    }
    for (n=0;n<2;n++) {
      MTF_GG_RB[n]=new Complex[M];
      MTF_GG_RB_REIM[(n<<1)]=  new double[M];
      MTF_GG_RB_REIM[(n<<1)+1]=new double[M];
      for (i = 0; i < M; i++) {
        MTF_GG_RB[n][i]=MTFComplex[6+n][i].minus(MTFComplex[0+n][i]).divideBy(MTFComplex[2+n][i].minus(MTFComplex[4+n][i]));
        MTF_GG_RB_REIM[n<<1][i]=(horizontalEdges?(-1.0):1.0)*MTF_GG_RB[n][i].Re();
        MTF_GG_RB_REIM[(n<<1)+1][i]=MTF_GG_RB[n][i].Im();
      }
    }
// private Complex[][]  MTF_GG_RB=          new Complex [2][];
// private double [][]  MTF_GG_RB_REIM    = new double  [4][];
//divideBy   

    /* calculate results of division of each channel MTF by the channel that has the maximal total energy.
    Based on assumption that:
    1 -  image was acquired with monochromatic light (no chromatic aberrations), and
    2 - one channel has much higher signal, so it can be considered as "actual" light input
    the other (lower) channels are considered to be a sum of the same signal as the main channel (just attenuated by mosaic filter) and
    the crosstalk signal from the main channel pixels.

    */
    double max=0.0;
    MTFMaximalChannel=0;
    for (n=0;n<4;n++) if ((MTF_DC[2*n]+MTF_DC[2*n+1])>max) {
      max=MTF_DC[2*n]+MTF_DC[2*n+1];
      MTFMaximalChannel=n;
    }
    for (n=0;n<8;n++) {
      switch (MTFMaximalChannel) {
        case 0: /* Green red line */
        case 3: /* Green blue line */
          masked= ((n >> 1) == MTFMaximalChannel) || ((n >> 1) == (MTFMaximalChannel^3));
          break;
        default: masked= ((n >> 1) == MTFMaximalChannel);
      }

      n1=(n & 1) + (MTFMaximalChannel << 1);
      MTFCRatios[n<<1]=      new double[M];
      MTFCRatios[(n<<1)+1]=  new double[M];
      MTFCRatiosA[n] =       new double [M];  /* Amplitude */
      MTFCRatiosAMasked[n] = new double [M];  /* Amplitude, masked */
      MTFCRatiosRe[n]=       new double [M];  /* Real part */
      MTFCRatiosReMasked[n]= new double [M];  /* Real part, Masked */


      MTFCRatiosP[n] = new double [M];  /* Phase */
      MTFRatioComplex[n]=new Complex[M];
      for (i=0;i<M;i++) {

//        d= (MTFCVectorsBoth[n1][i]*MTFCVectorsBoth[n1][i])+(MTFCVectorsBoth[n1+1][i]*MTFCVectorsBoth[n1+1][i]);
        if (MTFComplex[n1][i].abs()>0) MTFRatioComplex[n][i]=  MTFComplex[n][i].divideBy(MTFComplex[n1][i]);
        else  MTFRatioComplex[n][i]=  new Complex(0.0,0.0);
        MTFCRatios[n<<1][i]=      MTFRatioComplex[n][i].Re();
        MTFCRatios[(n<<1)+1][i]=  MTFRatioComplex[n][i].Im();
        MTFCRatiosA[n][i] =       MTFRatioComplex[n][i].abs();
        MTFCRatiosP[n][i] =       MTFRatioComplex[n][i].phase();
        MTFCRatiosAMasked[n][i] = masked?0.0:MTFCRatiosA[n][i];
        MTFCRatiosRe[n][i]=       MTFRatioComplex[n][i].Re();
        MTFCRatiosReMasked[n][i]= masked?0.0:MTFCRatiosRe[n][i];

      }
    }


  }

  public void showMTFRatios(boolean plot_amplitude, boolean plot_phase, boolean text, boolean complex) {
    if (MTFVectorsBoth[0]==null) {
      IJ.showMessage("Error","No data available, please run \"Process\" first");
      return;
    }
    if (plot_amplitude)  generatePlots(MTFCRatiosReMasked,"RE","OTF_Ratio_Re",        bayerColorsTwice,false);
    if (plot_phase)      generatePlots(MTFCRatiosP,"PHASE","OTF_Ratio_phase",bayerColorsTwice,false);
    if (text)    generatePlots(MTFCRatiosA, "MTF","OTF_Ratio",bayerColorsTwice,true,"Frequency"+colorCSVHeadersBoth,0,0);
    if (complex) {
        if (USE_COMPLEX) generatePlots(MTFCRatios,  "MTF","OTF_Ratio_complex",bayerColorsTwice,true,"Frequency"+colorCSVHeadersBoth,0,0);
        else             generatePlots(MTFCRatios,  "MTF","OTF_Ratio_complex",bayerColorsTwice,true,"Frequency"+colorCSVHeadersBothComplex,0,0);
    }
 }

  public void showMTFGGRB(boolean plot_reim, boolean text) {
    if (MTFVectorsBoth[0]==null) {
      IJ.showMessage("Error","No data available, please run \"Process\" first");
      return;
    }
    if (plot_reim)         generatePlots(MTF_GG_RB_REIM,"REIM", "Crosstalk_Re_Im",         redBlueColorsTwice,false);
    if (text ) {
        if (USE_COMPLEX) generatePlots(MTF_GG_RB_REIM,  "MTF","Crosstalk_complex",redBlueColorsTwice,true,"Frequency"+CSVHeadersBoth,0,0);
        else             generatePlots(MTF_GG_RB_REIM,  "MTF","Crosstalk_complex",redBlueColorsTwice,true,"Frequency"+CSVHeadersBothComplex,0,0);
    }
 }



//MTF_GG_RB_REIM
// private double[][]   MTFPVectorsBoth =   new double [8][];  /**Phases */
// private double[][]   MTFCVectorsBoth =   new double [16][]; /**Complex */

  public void showMTF(boolean plot_amplitude, boolean plot_phase, boolean plot_reim, boolean plot_reim_scaled, boolean text, boolean complex) {
   int n,i;
    if (MTFVectorsBoth[0]==null) {
        IJ.showMessage("Error","No data available, please run \"Process\" first");
        return;
    }

    if (plot_amplitude)  generatePlots(normalizePlots(MTFVectorsBoth,false),"MTF","MTF",        bayerColorsTwice,false);
    if (plot_phase)      generatePlots(MTFPVectorsBoth,"PHASE","OTF_phase",bayerColorsTwice,false);
    if (plot_reim || plot_reim_scaled){
        double[][] MTF_Re =   new double [8][MTFVectorsBoth[0].length];
        double[][] MTF_Im =   new double [8][MTFVectorsBoth[0].length];
        for (n=0;n<8;n++) for (i=0;i<MTFVectorsBoth[0].length;i++) {
          MTF_Re[n][i]=MTFCVectorsBoth[2*n][i];
          MTF_Im[n][i]=MTFCVectorsBoth[2*n+1][i];
        }
        if (plot_reim) {
          generatePlots(MTF_Re,"RE","re(OTF)",        bayerColorsTwice,false);
          generatePlots(MTF_Im,"IM","im(OTF)",        bayerColorsTwice,false);
        }
        if (plot_reim_scaled) {
          generatePlots(MTF_Re,"RES","re(OTF)_scaled",        bayerColorsTwice,false);
          generatePlots(MTF_Im,"IMS","im(OTF)_scaled",        bayerColorsTwice,false);
        }
    }
    if (text)    generatePlots(MTFVectorsBoth, "MTF","MTF",bayerColorsTwice,true,"Frequency"+colorCSVHeadersBoth,0,0);
    if (complex) {
        if (USE_COMPLEX) generatePlots(MTFCVectorsBoth,"MTF","OTF",bayerColorsTwice,true,"Frequency"+colorCSVHeadersBoth,0,0);
        else             generatePlots(MTFCVectorsBoth,"MTF","OTF",bayerColorsTwice,true,"Frequency"+colorCSVHeadersBothComplex,0,0);
    }
 }


  public void showESF(boolean actual, boolean normalized, boolean text, boolean actual_diff, boolean normalized_diff) {
    int n, bin;
    int sSize4=sSize*4;
    double[][] ESFBin      = new double[8][sSize*4];
    double[][] ESFBinRmBGmG   = new double[4][sSize*4];
    if (binDataBoth==null) {
        IJ.showMessage("Error","No data available, please run \"Process\" first");
        return;
    }
    for (n=0;n<4;n++) for (bin=0; bin<sSize4; bin++) {
      ESFBin[(n<<1)  ][bin]=binDataBoth[0][n][bin];
      ESFBin[(n<<1)+1][bin]=binDataBoth[1][n][bin];
    }
    for (bin=0; bin<sSize4; bin++) {
      ESFBinRmBGmG[0][bin]= binDataBoth[0][1][bin]-binDataBoth[0][2][bin];
      ESFBinRmBGmG[1][bin]= binDataBoth[0][3][bin]-binDataBoth[0][0][bin];
      ESFBinRmBGmG[2][bin]= binDataBoth[1][1][bin]-binDataBoth[1][2][bin];
      ESFBinRmBGmG[3][bin]= binDataBoth[1][3][bin]-binDataBoth[1][0][bin];
    }
    if (actual) generatePlots(ESFBin,"DV4","ESF",bayerColorsTwice,false);
    if (normalized) generatePlots(normalizePlots(ESFBin,true),"DV4","ESF_normalized",bayerColorsTwice,false);
    if (text) generatePlots(ESFBin,"DV4","ESF",bayerColorsTwice,true,"Pixel"+colorCSVHeadersBoth,0,0);
    if (actual_diff) generatePlots(ESFBinRmBGmG,"DV4","ESF_diff",redGreenColorsTwice,false);
    if (normalized_diff) generatePlots(normalizePlots(ESFBinRmBGmG,true),"DV4","ESF_normalized_diff",redGreenColorsTwice,false);
 }

  public void showLSF(boolean actual, boolean normalized, boolean text, boolean actual_diff, boolean normalized_diff) {
    int n, bin,edgeSign;
    int sSize4=sSize*4;
    double[][] LSFBoth   = new double[8][sSize4];
    double[][] LSFRmBGmG   = new double[4][sSize*4];
    if (binDataDiffBoth==null) {
        IJ.showMessage("Error","No data available, please run \"Process\" first");
        return;
    }
    for (edgeSign=0; edgeSign<2; edgeSign++) for (n=0;n<4;n++)  for (bin=0; bin<sSize*4; bin++) {
      LSFBoth[(n<<1)+edgeSign][bin]=binDataDiffBoth[edgeSign][n][bin];
    }
    for (bin=0; bin<sSize4; bin++) {
      LSFRmBGmG[0][bin]= binDataDiffBoth[0][1][bin]-binDataDiffBoth[0][2][bin];
//      LSFRmBGmG[1][bin]= binDataDiffBoth[1][1][bin]-binDataDiffBoth[1][2][bin];
//      LSFRmBGmG[2][bin]= binDataDiffBoth[0][3][bin]-binDataDiffBoth[0][0][bin];
      LSFRmBGmG[1][bin]= binDataDiffBoth[0][3][bin]-binDataDiffBoth[0][0][bin];
      LSFRmBGmG[2][bin]= binDataDiffBoth[1][1][bin]-binDataDiffBoth[1][2][bin];
      LSFRmBGmG[3][bin]= binDataDiffBoth[1][3][bin]-binDataDiffBoth[1][0][bin];
    }
    if (actual) generatePlots(LSFBoth,"DV4","LSF",bayerColorsTwice,false);
    if (normalized) generatePlots(normalizePlots(LSFBoth,false),"DV4","LSF_normalized",bayerColorsTwice,false);
    if (text) generatePlots(LSFBoth,"DV4","LSF",bayerColorsTwice,true,"Pixel"+colorCSVHeadersBoth,0,0);
    if (actual_diff) generatePlots(LSFRmBGmG,"DV4","LSF_diff",redGreenColorsTwice,false);
    if (normalized_diff) generatePlots(normalizePlots(LSFRmBGmG,false),"DV4","LSF_normalized_diff",redGreenColorsTwice,false);
  }

/*
 Color [] redGreen_colors=    { Color.red, Color.green, Color.black };
 Color [] redGreenColorsTwice= { Color.red, Color.green, Color.red, Color.green, Color.black  };

*/

  public void processWindowEvent(WindowEvent e) {
    super.processWindowEvent(e);
    if (e.getID()==WindowEvent.WINDOW_CLOSING) {
      instance = null;
    }
  }

 public boolean showDialog() {
   int i;
// private static boolean [] SHOW_GGRB={false,true,true,true,true,true,true};

   String[] ESF_labels={"Show_ESF","Show_ESF_actual?","Show_ESF_normalized?","Show_ESF_as_text?"};
   String[] LSF_labels={"Show_LSF","Show_LSF_actual?","Show_LSF_normalized?","Show_LSF_as_text?"};
   String[] MTF_labels={"Show_MTF","Show_MTF_plot?","Show_OTF_phase?","Show_OTF_re/im?","Show_OTF_re/im scaled?","Show_MTF_as_text?","Show_OTF_as_text?"};
   String[] GGRB_labels={"Show_crosstalk_plots?","GGRB_ESF?","GGRB_ESF_normalized?","GGRB_LSF?","GGRB_LSF_normalized?","F[(Gb-Gr)]/F[R-B]?","F[(Gb-Gr)]/F[R-B]_as_text?"};
   String [] Edge_Algorithms={"Centroid","Two ESF crossing points","Four ESF crossing points", "Centroid with double windowing"};
   GenericDialog gd = new GenericDialog("Pixel Crosstalk parameters");
   gd.addStringField("Image title (updated automatically): ", title_src, 80);
   gd.addNumericField("Conversion_strip_width:",         sSize, 0);
   gd.addCheckbox("Monochrome (color filtered) image? ", IS_MONOCHROMATIC);
//   gd.addNumericField("Edge_location_algorithm_number:", EMODE, 0);
   gd.addChoice("Edge_location_algorithm", Edge_Algorithms, Edge_Algorithms[EMODE]);
   gd.addCheckbox("Truncate rows to equalize phase rotations? ", TRUNCATE_ROWS);

   gd.addNumericField("LSF_calculation_window (fraction of full, 0.0 - no window):",    LSFHamming, 2);
   gd.addCheckbox("Output complex numbers?",  USE_COMPLEX);
   gd.addCheckbox("Show_edge_approximation_info?",  SHOW_edgeinfo);
   gd.addCheckboxGroup(1,4,ESF_labels, SHOW_ESF);
   gd.addCheckboxGroup(1,4,LSF_labels, SHOW_LSF);

   gd.addCheckboxGroup(1,7,MTF_labels, SHOW_MTF);
   gd.addCheckboxGroup(1,7,GGRB_labels, SHOW_GGRB);
   gd.addCheckbox("Show calculated crosstalk (monochrome mode, red or blue filter required)",SHOW_CROSSTALK);
   gd.addCheckbox("Allow_full_image?", ALLOW_FULL_IMAGE);
// private static double bg2gr_corr=1.0; /* difference in Green sensitivity Gb/Gr*/
   gd.addNumericField("Gb/Gr sensitivity correction:",GB2GR_CORR, 5);
 
   gd.addNumericField("Debug_Level:",      MASTER_DEBUG_LEVEL, 0);
   gd.showDialog();
   if (gd.wasCanceled()) return false;
   title_src = gd.getNextString();
   sSize=1;
   for (i=(int) gd.getNextNumber(); i >1; i>>=1) sSize <<=1; /* make sSize to be power of 2*/

   IS_MONOCHROMATIC=gd.getNextBoolean();
//   EMODE=(int) gd.getNextNumber();
   EMODE=gd.getNextChoiceIndex();
   TRUNCATE_ROWS=gd.getNextBoolean();
   LSFHamming=  gd.getNextNumber();
   USE_COMPLEX=gd.getNextBoolean();
   SHOW_edgeinfo=gd.getNextBoolean();
   SHOW_ESF[0]=gd.getNextBoolean();
   SHOW_ESF[1]=gd.getNextBoolean();
   SHOW_ESF[2]=gd.getNextBoolean();
   SHOW_ESF[3]=gd.getNextBoolean();

   SHOW_LSF[0]=gd.getNextBoolean();
   SHOW_LSF[1]=gd.getNextBoolean();
   SHOW_LSF[2]=gd.getNextBoolean();
   SHOW_LSF[3]=gd.getNextBoolean();

   SHOW_MTF[0]=gd.getNextBoolean();
   SHOW_MTF[1]=gd.getNextBoolean();
   SHOW_MTF[2]=gd.getNextBoolean();
   SHOW_MTF[3]=gd.getNextBoolean();
   SHOW_MTF[4]=gd.getNextBoolean();
   SHOW_MTF[5]=gd.getNextBoolean();
   SHOW_MTF[6]=gd.getNextBoolean();

   SHOW_GGRB[0]=gd.getNextBoolean();
   SHOW_GGRB[1]=gd.getNextBoolean();
   SHOW_GGRB[2]=gd.getNextBoolean();
   SHOW_GGRB[3]=gd.getNextBoolean();
   SHOW_GGRB[4]=gd.getNextBoolean();
   SHOW_GGRB[5]=gd.getNextBoolean();
   SHOW_GGRB[6]=gd.getNextBoolean();

   SHOW_CROSSTALK=gd.getNextBoolean();
   ALLOW_FULL_IMAGE=gd.getNextBoolean();
   GB2GR_CORR=gd.getNextNumber();
   MASTER_DEBUG_LEVEL=  (int) gd.getNextNumber();
   return true;
 }
/*
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
       IJ.showMessage("Debug","dialogItemChanged");
        return false;
    }
*/
  private void showEdgesCoeff(double[][][] color_coeff_abc, int halfHeight) {
    int nValid=color_coeff_abc.length;
    int i,n, edge;
    double rmsMax=0.0;
    StringBuffer sb = new StringBuffer();
    double [][][] data_total=new double[2][4][6]; /* calculate totals for each edge, each Bayer component */
    for (edge=0;edge<2;edge++) for (n=0;n<4;n++) for (i=0;i<6;i++) data_total[edge][n][i]=0.0;

    for (i=0;i<nValid;i++) for (n=0;n<4;n++) {
      edge=posEdges[halfHeight][i]; /* 0 - black-> white, 1 white->black*/
      data_total[edge][n][0]+=color_coeff_abc[i][n][0]*color_coeff_abc[i][n][5];
      data_total[edge][n][1]+=color_coeff_abc[i][n][1]*color_coeff_abc[i][n][5];
      data_total[edge][n][2]+=color_coeff_abc[i][n][2]*color_coeff_abc[i][n][5];
      data_total[edge][n][3]+=color_coeff_abc[i][n][3]*color_coeff_abc[i][n][5];
      data_total[edge][n][4]+=color_coeff_abc[i][n][4]*color_coeff_abc[i][n][5];
      data_total[edge][n][5]+=                         color_coeff_abc[i][n][5];
      if (color_coeff_abc[i][n][3]>rmsMax) rmsMax=color_coeff_abc[i][n][3];
    }




    for (edge=0;edge<2;edge++) for (n=0;n<4;n++) {
      data_total[edge][n][0]/=data_total[edge][n][5];
      data_total[edge][n][1]/=data_total[edge][n][5];
      data_total[edge][n][2]/=data_total[edge][n][5];
      data_total[edge][n][3]/=data_total[edge][n][5];
      data_total[edge][n][4]/=data_total[edge][n][5];
    }
    for (edge=0;edge<2;edge++) {
      for (n=0;n<4;n++) {
        sb.append("total"+
                "\t"+((edge>0)?" \\_":" _/")+
                "\t"+n+
                "\t"+IJ.d2s((data_total[edge][n][2]*halfHeight*halfHeight),3)+
                "\t"+IJ.d2s((data_total[edge][n][1]*halfHeight),3)+
                "\t"+IJ.d2s((data_total[edge][n][0]),3)+
                "\t"+IJ.d2s((data_total[edge][n][3]),3)+
               (IS_MONOCHROMATIC?( "\t"+IJ.d2s((data_total[edge][n][4]),3)):"")+
                "\t"+IJ.d2s((data_total[edge][n][5]),0)+
                "\n"); 
      }
      sb.append("-\t-\t-\t-\t-\t-\t-\t\n");
    }

    for (i=0;i<nValid;i++) {
      for (n=0;n<4;n++) {
        sb.append( i+
                "\t"+((posEdges[halfHeight][i]>0)?" \\_":" _/")+
                "\t"+n+
                "\t"+IJ.d2s((color_coeff_abc[i][n][2]*halfHeight*halfHeight),3)+
                "\t"+IJ.d2s((color_coeff_abc[i][n][1]*halfHeight),3)+
                "\t"+IJ.d2s((color_coeff_abc[i][n][0]),3)+
                "\t"+IJ.d2s((color_coeff_abc[i][n][3]),3)+
               (IS_MONOCHROMATIC?( "\t"+IJ.d2s((color_coeff_abc[i][n][4]),3)):"")+
                "\t"+IJ.d2s((color_coeff_abc[i][n][5]),0)+
                "\n"); 
      }
      sb.append("\t\t\t\t\t\t\t\n");
    }
    new TextWindow(title_src+" "+(horizontalEdges?"H-":"V-")+"Edges, RMS<"+IJ.d2s(rmsMax,3), "Edge#\tSlope\tBayer\tA\tB\tC\tRMS"+
                                   (IS_MONOCHROMATIC?"\tSkew":"")+"\tUsed rows", sb.toString(), 600, 600);
  }

//ALLOW_FULL_IMAGE
  public boolean openImage(){
    imp_src = WindowManager.getCurrentImage();
    if (imp_src==null){
      IJ.showMessage("Error","There are no images open\nProcess canceled");
      return false;
    }
    title_src=imp_src.getTitle();
    Roi roi_src= imp_src.getRoi();    
    if (imp_src.getType() !=ImagePlus.GRAY32 ) {
      if ((imp_src.getType() ==ImagePlus.GRAY8 ) ||
          (imp_src.getType() ==ImagePlus.GRAY16) ) {
         IJ.showStatus("Converting source image to gray 32 bits (float)");
         new ImageConverter(imp_src).convertToGray32();
      } else {
        IJ.showMessage("Error","Image should be Bayer array as a grayscale (8,16 or 32 bits)");
        return false;
      }
    }
    if (roi_src==null){
      if (ALLOW_FULL_IMAGE) {
        imp_src.setRoi(0, 0, imp_src.getWidth(), imp_src.getHeight());
        roi_src= imp_src.getRoi();    
      } else {
        IJ.showMessage("Error","No selection (line or rectangle) in the source image.\n"+
                               "You may allow processing of the full image in \"Configure\"");
        return false; /* Full image selected */
      }
    }
    Rectangle r=roi_src.getBounds();
/* align ROI to Bayer */
    r.width &= ~1;
    r.height &= ~1;
    r.x &= ~1;
    r.y &= ~1;
    selecWidth=r.width;
    selecHeight=r.height;
    imp_src.setRoi(r);
    return true;
  }
  

  double [][] generateSelection(ImagePlus imp,  boolean horizontal) {
    int i,i1,j,base;
    Roi roi_src= imp.getRoi();    
    Rectangle r=roi_src.getBounds();
    ImageProcessor ip=imp.getProcessor();
    float [] pixels =(float[]) ip.getPixels();
    selecWidth=horizontal?r.height:r.width;
    selecHeight=horizontal?r.width:r.height;
    int w=ip.getWidth();
    if (selecWidth<sSize) return null; //sample size is bigger than selection width
/* selection already aligned to color mosaic */
    double [][] esf=new double[selecHeight][selecWidth];
    if (DEBUG_LEVEL>4) IJ.showMessage("Debug","selecWidth="+selecWidth+"\tselecHeight="+selecHeight);

    if (horizontal) { /* (nearly) horizontal edges */
      for (i=0; i< selecWidth; i++) {
         base=((r.y+i) * w+r.x);
         i1=selecWidth-i-1;
         for (j=0;j< selecHeight; j++) {
            if (((i&1)==0) && ((j&1)==0)) esf[j][i1]= GB2GR_CORR* pixels[base+j]; /* Gr*=GB2GR_CORR */
            else                          esf[j][i1]= pixels[base+j];
         }
//GB2GR_CORR
      }
    } else {  /* (nearly) vertical edges */
      for (i=0; i< selecHeight; i++) {
         base=((r.y+i) * w+r.x);
         for (j=0;j< selecWidth; j++) {
            if (((i&1)==0) && ((j&1)==0)) esf[i][j]= GB2GR_CORR* pixels[base+j]; /* Gr*=GB2GR_CORR */
            else                          esf[i][j]= pixels[base+j];
         }
      }
    }
    return esf;
  }

/**
  Calculates 2-d integer array of vertical (in esf) edges that run through all height (esf[][] first index)
  maxJump - maximal difference between edge crossing subsequent lines.
  Returns array of edge positions, null on failure
*/

  int [][] findEdges(double[][] esf, int maxJump) {
    int halfHeight=selecHeight/2;
    int k,i,j,k2,k21, sign,pos;
    int [][] intEdges=new int [halfHeight][selecWidth]; /* used later !!!*/
    int [] initialSign=new int[halfHeight];
    int [] numEdges=new int[halfHeight];
    double [] smoothLine=new double [selecWidth]; // 0,1,selecWidth-1 - undefined
    double mx,mn, thresh;
    int halfSize=sSize/2;
//    int quaterSize=sSize/4;

    for (k=0;k<halfHeight;k++) for (i=0;i<selecWidth;i++) intEdges[k][i]=0;
    for (k=0;k<halfHeight;k++) {
      initialSign[k]=0;
      numEdges[k]=0;
      k2=2*k;
      k21=k*2+1;
      mx=0.0;
      mn=0.0;

      for (i=2; i<(selecWidth-1); i++ ) {
        smoothLine[i]=esf[k2][i+1]+esf[k21][i+1]+esf[k2][i]+esf[k21][i]-esf[k2][i-2]-esf[k21][i-2]-esf[k2][i-1]-esf[k21][i-1];
        if (smoothLine[i] > mx) mx=smoothLine[i];
        if (smoothLine[i] < mn) mn=smoothLine[i];
      }
/*  mx and mn are supposed to be of opposite sign here, set thereshold 1/4 of the minimal of the absolute values of mx and mn */
      if ((mn*mx)>=0.0) return null; /* no edges found */
      thresh=mx/4;
      if (thresh > -(mn/4)) thresh = -(mn/4);
      sign=0; // 0 - unknown yet, 1 - last was max, -1 - last was min
      mx=0.0;
      pos=0;
      for (i=2;i<(selecWidth-1);i++) {
        if (smoothLine[i] > thresh) {
          if (sign>0) {
            if (smoothLine[i] > mx) {
              mx=smoothLine[i];
              pos=i;
            }
          } else  {
            if (sign<0) {
              intEdges[k][pos]=-1;
              numEdges[k]++;
            } else initialSign[k]=1;
            sign=1;
            mx=smoothLine[i];
            pos=i;
          }
        } else if (smoothLine[i] < -thresh) {
          if (sign<0) {
            if (smoothLine[i] < mx) {
              mx=smoothLine[i];
              pos=i;
            }
          } else {
            if (sign>0) {
              intEdges[k][pos]=1;
              numEdges[k]++;
            } else initialSign[k]=-1;
            sign=-1;
            mx=smoothLine[i];
            pos=i;
          }
        }
      }
      if  (sign<0) {
        intEdges[k][pos]=-1;
        numEdges[k]++;
      } else if (sign>0) {
        intEdges[k][pos]=1;
        numEdges[k]++;
      }
    }
//    IJ.showMessage("calcEdges() debug 1  halfHeight="+halfHeight,"First line has "+numEdges[0]+" edges");
/* for now - considering only edges that run the full height */
    boolean []   valid=new boolean [numEdges[0]]; // if the edge runs all the way from top to bottom
//    int [][] posEdges=new int [halfHeight][numEdges[0]]; // horizontal positions of edges
    int [][] posEdges=new int [halfHeight+1][numEdges[0]]; // horizontal positions of edges, [halfHeight][i] - signs of edges (0/1)
    for (i=0; i<numEdges[0];i++)  valid[i]=true;
    j=0;
    for (i=0; i<(selecWidth-1);i++) if  (intEdges[0][i] !=0) posEdges[0][j++]=i; // set first row

    for (k=0; k<halfHeight;k++) {

/* Verify current line does not get out of halfSize from the ends (enough room on both sides of the edge) */
      for (i=0; i<numEdges[0]; i++) if (valid[i]) {
        if ((posEdges[k][i]<halfSize) || (posEdges[k][i]>= (selecWidth-halfSize))) valid[i] = false;
      }

      if (k< (halfHeight-1)) {
/* See if the edge is continued to the next line */
        for (i=0; i<numEdges[0]; i++) if (valid[i]) {
          posEdges[k+1][i]=-1;
          for (j=0; j<=maxJump; j=-j+((j<=0)?1:0)) {
            if (intEdges[k][posEdges[k][i]]==intEdges[k+1][posEdges[k][i]+j]) {
              posEdges[k+1][i]=posEdges[k][i]+j;
              break;
            }
          }
          if (posEdges[k+1][i]<0) valid[i]=false;
        }
      }

    }
/* compact posEdges by removing invalid columns, calculate nValid; */
    int nValid=0;
    for (i=0; i<numEdges[0]; i++) if (valid[i]) {
      if (i!=nValid) for (k=0; k< halfHeight; k++) posEdges[k][nValid]=posEdges[k][i];
      nValid++;
    }
    for (i=0;i<nValid; i++ )  posEdges[halfHeight][i]= (intEdges[0][posEdges[0][i]]>0)?0:1; /* 0 - black-> white, 1 white->black*/
    if (DEBUG_LEVEL>4) IJ.showMessage("calcEdges() debug 2","First line has "+numEdges[0]+" edges, of the "+nValid+" run all the way");
    if (nValid==0) return null;

/* Copy posEdges to a new array that has all columns valid, use this array as return value of this function.
    Last row contains signes of the edges (0/1)*/
    int [][] posEdges1=new int [halfHeight+1][nValid];
//  if (DEBUG_LEVEL>10)IJ.showMessage("Debug-22","posEdges1.length="+((posEdges1==null)?"null":posEdges1.length)+"\nhalfHeight="+halfHeight );
    for (k=0; k<=halfHeight; k++) for (i=0;i<nValid; i++ ) posEdges1[k][i]=posEdges[k][i];
    return posEdges1;
  }
  
/**
 Find second order polinomial approximation for each edge.
 returns polinomial coefficient (a,b,c), rms for each line, for each Bayer component
*/
  private double[][][] approximateEdges(int [][] posEdges) {
    int nValid=posEdges[0].length;
    int halfHeight=posEdges.length-1;
    int halfSize=sSize/2;
    int quarterSize=sSize/4;
    int HammingHalfWidth=(int) (0.5*LSFHamming*halfSize);
    double [] LSFColorVector= new double [halfSize+1];
    double [][][] posColorEdges=new double [halfHeight][4][nValid]; // horizontal positions of edges
    double ax,a, cx,aw,s,s1,s2,y0,y1,y2,y3,yy,pa,pb,pc,pd,x0,dx;
//    IJ.showMessage("Debug","approximateEdges: HammingHalfWidth="+HammingHalfWidth);
    int i, k, k2, n, n1, pos, byr_h,byr_v, byr_index; //bayer position hor/vert (0/1)
    for (k=0; k<halfHeight; k++) {
      for (byr_v=0;byr_v<2; byr_v++) {
        k2=2*k+byr_v;
        for (byr_h=0;byr_h<2; byr_h++) {
          byr_index=horizontalEdges?(((byr_h^1)<<1) + byr_v):((byr_v<<1)+byr_h);
          for (n=0; n<nValid; n++) {
            pos=2*(posEdges[k][n]/2)-halfSize+byr_h;
            switch (EMODE) {
              case 0: /* centroid -> multiply by window -> recalcualte centroid */
                for (i=0; i<halfSize; i++) LSFColorVector[i]=ESFArray[k2][pos+2*i+2]-ESFArray[k2][pos+2*i];
                a=0.0;
                
                ax=0.0;
                for(i=0;i<halfSize;i++){
                  a +=LSFColorVector[i];
                  ax+=LSFColorVector[i]*i;
                }
                cx=ax/a+0.5; /* center position without window */
                if (HammingHalfWidth> 0) {
                  a=0.0;
                  ax=0.0;
                  for(i=0;i<halfSize;i++) if (Math.abs(i-cx)<=HammingHalfWidth){
                    aw=LSFColorVector[i]*(0.54+0.46*Math.cos(Math.PI*(i-cx)/HammingHalfWidth));
                    a +=aw;
                   ax+=aw*i;
                 }
                }
                posColorEdges[k][byr_index][n]= 2*(ax/a)+pos+0.5; ///pixel position (in original pixels) of the maximum
                break;
              case 3: /* centroid -> multiply by window -> recalcualte centroid, reduce window twice - recalculate again */
                double HammingQuaterWidth=HammingHalfWidth/1;
                for (i=0; i<halfSize; i++) LSFColorVector[i]=ESFArray[k2][pos+2*i+2]-ESFArray[k2][pos+2*i];
                a=0.0;
                ax=0.0;
                for(i=0;i<halfSize;i++){
                  a +=LSFColorVector[i];
                  ax+=LSFColorVector[i]*i;
                }
                cx=ax/a+0.5; /* center position without window */
                if (HammingHalfWidth> 0) {
                  a=0.0;
                  ax=0.0;
                  for(i=0;i<halfSize;i++) if (Math.abs(i-cx)<=HammingQuaterWidth){
                    aw=LSFColorVector[i]*(0.54+0.46*Math.cos(Math.PI*(i-cx)/HammingQuaterWidth));
                    a +=aw;
                    ax+=aw*i;
                  }
                  cx=ax/a+0.5; /* center position with full Hamming window */
                  a=0.0;
                  ax=0.0;
                  for(i=0;i<halfSize;i++) if (Math.abs(i-cx)<=HammingHalfWidth){
                    aw=LSFColorVector[i]*(0.54+0.46*Math.cos(Math.PI*(i-cx)/HammingHalfWidth));
                    a +=aw;
                    ax+=aw*i;
                  }
                }
                posColorEdges[k][byr_index][n]= 2*(ax/a)+pos+0.5; ///pixel position (in original pixels) of the maximum
                break;



              case 1: /* average leftmost quarter, rightmost quarter, find middle between them, find crossing point with linear interpolation.*/
              case 2: /* average leftmost quarter, rightmost quarter, find middle between them, find crossing point with 3-rd polynome interpolation. */
                s1=0.0;
                s2=0.0;
                for (i=0;i<quarterSize;i++) {
                	
                  s1+=ESFArray[k2][pos+2*i];
                  s2+=ESFArray[k2][pos+2*i+halfSize];
                }
                s1/=quarterSize;
                s2/=quarterSize;
                s=0.5*(s1+s2);
                for (i=quarterSize; (i< (sSize-quarterSize) && ((s1-s)*(ESFArray[k2][pos+2*i]-s) >0)); i++) ;
                /* edge between i-1 and i*/
                y1=ESFArray[k2][pos+2*i-2]-s;
                y2=ESFArray[k2][pos+2*i]-s;
                if (y2==0.0) {
                  posColorEdges[k][byr_index][n]= pos+2*i-0.5; /* unlikely exact match */
                  break;
                }
                if ((y1)*(y2)>=0) { /* Should not be so, retry wider */
                  for (i=2; (i< (sSize-2) && ((s1-s)*(ESFArray[k2][pos+2*i]-s) >0)); i++) ;
                  /* edge between i-1 and i*/
                  y1=ESFArray[k2][pos+2*i-2]-s;
                  y2=ESFArray[k2][pos+2*i]-s;
                }
                posColorEdges[k][byr_index][n]= pos+2*((i-1)+(y1)/(y1-y2))-0.5; /* linear, uase a backup for cubic */
                if (EMODE==1) break;
                y0=ESFArray[k2][pos+2*i-4]-s;
                y3=ESFArray[k2][pos+2*i+2]-s;
/* y= pa*x^3+pb*x^2+pc*x+pd*/
                pd=          y0;
                pa= -(1.0/6)*y0 +0.5*y1 -0.5*y2 +(1.0/6)*y3;
                pb=          y0 -2.5*y1   +2*y2     -0.5*y3;
                pc=-(11.0/6)*y0   +3*y1 -1.5*y2 +(1.0/3)*y3;
/* solve cubic equation and select solution between 1 and 2, if none - use linear */
/* simplified solution - add the point from linear interpolation, linear interpolate before/after */
                x0=y1/(y1-y2) +1.0;
                yy=pa*x0*x0*x0+pb*x0*x0+pc*x0+pd;
                if (y1*yy == 0.0) {
                  posColorEdges[k][byr_index][n]= pos+2*(i-2+x0)-0.5; /* unlikely exact match */
                } else if (y1*yy > 0.0) { /* interpolate between yy and y2 */
                  posColorEdges[k][byr_index][n]= pos+2*((i-2+x0)+ (2-x0)*yy/(yy-y2))-0.5;
                } else { /* interpolate between y1 and yy */
                  posColorEdges[k][byr_index][n]= pos+2*(i-1+ (x0-1)*y1/(y1-yy))-0.5;
                }
                break;
            }
          }
        }
      }
    }
//    color_coeff_abc=new double [nValid][4][4];
    color_coeff_abc=new double [nValid][4][6]; /* [4] (monochromatic only) - average difference between approximated edge and this color component,
                                                   [5] - number of rows used */
    double SX4,SYX2,SX3,SX2,SYX,SY,SX,S0;
    double x,rms, skew, pmax,err;
    for (i=0; i<nValid; i++) for (n=0;n<4;n++) if ((!IS_MONOCHROMATIC) || (n==3)) {
      SX4=0.0;
      SX3=0.0;
      SX2=0.0;
      SX=0.0;
      SYX2=0.0;
      SYX=0.0;
      SY=0.0;
      S0=0.0;
/* TODO: Make weighted contributions of different components (or even the main one) */
      if (IS_MONOCHROMATIC) { /* gets here with (n==3) only, combine all Bayer components, taking care of shifts */
        for (n1=0;n1<4;n1++) {
          dx=0.5*bayer_loc[horizontalEdges?1:0][n1][0];
          for(k=0;k<halfHeight;k++){
            pmax=posColorEdges[k][n1][i];
            x=k+dx;
            SX4+=x*x*x*x;
            SYX2+=pmax*x*x;
            SX3+=x*x*x;
            SX2+=x*x;
            SYX+=pmax*x;
            SY+=pmax;
            SX+=x;
            S0+=1.0;
          }
        }
      } else {
        for(k=0;k<halfHeight;k++){
          pmax=posColorEdges[k][n][i];
          x=k;
          SX4+=x*x*x*x;
          SYX2+=pmax*x*x;
          SX3+=x*x*x;
          SX2+=x*x;
          SYX+=pmax*x;
          SY+=pmax;
          SX+=x;
          S0+=1.0;
        }
      }
/**
(1)         a*SX4 +b*SX3 + c*SX2 -SYX2 =0
(2)         a*SX3 +b*SX2 + c*SX  -SYX  =0
(3)         a*SX2 +b*SX  + c*S0  -SY   =0
(1a)        a*SX4*SX  +b*SX3*SX  + c*SX2*SX  -SYX2*SX =0
(2a)        a*SX3*SX2 +b*SX2*SX2 + c*SX*SX2  -SYX*SX2  =0
(1a-2a)     a*(SX4*SX-SX3*SX2) + b* (SX3*SX-SX2*SX2) - (SYX2*SX-SYX*SX2)
(1a2a)      b=  ((SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2)
(2b)        a*SX3*S0 +b*SX2*S0 + c*SX*S0  -SYX*S0  =0
(3b)        a*SX2*SX +b*SX*SX  + c*S0*SX  -SY*SX   =0
(2b-3b)     a*(SX3*S0-SX2*SX) +b*(SX2*S0-SX*SX)  -(SYX*S0-SY*SX)  =0
(2b3b1a2a)  a*(SX3*S0-SX2*SX) +((SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX)  -(SYX*S0-SY*SX)  =0
            a*(SX3*S0-SX2*SX) +(SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) -  a* (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX)  - (SYX*S0-SY*SX)  =0
            a*((SX3*S0-SX2*SX) - (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX)) + (SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) - (SYX*S0-SY*SX)  =0
            a=((SYX*S0-SY*SX)-(SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) )/((SX3*S0-SX2*SX) - (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX))
*/
      color_coeff_abc[i][n][2]=((SYX*S0-SY*SX)-(SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) )/((SX3*S0-SX2*SX) - (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX));
/**(1a2a)   b=  (SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) */
      color_coeff_abc[i][n][1]=((SYX2*SX-SYX*SX2) -  color_coeff_abc[i][n][2]* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2);
/**(3)      a*SX2 +b*SX  + c*S0  -SY   =0
            c =(SY -  a*SX2 - b*SX)/S0 */
      color_coeff_abc[i][n][0]=(SY -  color_coeff_abc[i][n][2]*SX2 - color_coeff_abc[i][n][1]*SX)/S0;

      if (IS_MONOCHROMATIC) { /* gets here with (n==3) only, combine all Bayer components, taking care of shifts */
        for (n1=0;n1<3;n1++) for (k=0;k<3;k++)  color_coeff_abc[i][n1][k]=color_coeff_abc[i][3][k]; /* copy coefficients to all other components */
        for (n1=0;n1<4;n1++) {
          dx=0.5*bayer_loc[horizontalEdges?1:0][n1][0];
          rms=0;
          skew=0.0;
          for(k=0;k<halfHeight;k++){
            pmax=posColorEdges[k][n1][i];
            x=k+dx;
            err=pmax-color_coeff_abc[i][n1][2]*x*x-color_coeff_abc[i][n1][1]*x-color_coeff_abc[i][n1][0];
            skew+=err;
            rms+=err*err;
          }
          color_coeff_abc[i][n1][3]=Math.sqrt(rms/selecHeight);
          color_coeff_abc[i][n1][4]=skew/selecHeight;
        }
      } else {
        rms=0;
        skew=0.0;
        for(k=0;k<halfHeight;k++){
          x=k;
          pmax=posColorEdges[k][n][i];
          err=pmax-color_coeff_abc[i][n][2]*x*x-color_coeff_abc[i][n][1]*x-color_coeff_abc[i][n][0];
          skew+=err;
          rms+=err*err;
        }
        color_coeff_abc[i][n][3]=Math.sqrt(rms/selecHeight);
        color_coeff_abc[i][n][4]=skew/selecHeight; /* should be 0 */
      } 
    }
    return color_coeff_abc;
  }

  private void binESF() {
    int nValid=color_coeff_abc.length;
    binDataBoth   = new double[2][4][sSize*4+1]; /* global double[][][]  */
    int [][][] binNumberBoth = new int [2][4][sSize*4+1];
/**
Scale data so all the ESF curves would have the same amplitude before adding to bins - that will decrease fluctuations caused by binNumber modulation
*/
    double[][] binMinBoth  = new double[2][4]; 
    double[][] binAmpBoth  = new double[2][4]; 
    int[][]    binScaleBothNumber  = new int[2][4]; 
    double FirstQuaterTotal,LastQuaterTotal,ThisESFMin,ThisESFAmp;
    int    FirstQuaterNumber, LastQuaterNumber;
    int i,n,k,edgeSign,edgeNumber, binShift, hor, hor_min, hor_max, bin;
    double dk;
    int sSize4=sSize*4;
    int halfHeight=selecHeight/2;
    int []rb={1,3,0,2};
    int nb;
    double ph0,ph1,phLim;
    int thisHeight;
    int slantDir;
    for (edgeSign=0;edgeSign<2;edgeSign++) for (n=0;n<4;n++){
       binMinBoth[edgeSign][n]=0.0;
       binAmpBoth[edgeSign][n]=0.0;
       binScaleBothNumber[edgeSign][n]=0;
       for (i=0;i<=sSize4;i++) {
         binDataBoth[edgeSign][n][i]=0.0;
         binNumberBoth[edgeSign][n][i]=0;
       }
    }
    for (edgeNumber=0; edgeNumber<nValid; edgeNumber++) {
      for (n=0;n<4;n++) {
        thisHeight=halfHeight; /* Will backup to this if there are too few edges */
        edgeSign=posEdges[halfHeight][edgeNumber]; /* 0 - black-> white, 1 white->black*/
/* Equalize number of phases, taking into account that we are binning each Bayer component separately,
    so number of pixel crossings should be (approximately) multiple of 2
*/
        dk=thisHeight-1;
        ph0=color_coeff_abc[edgeNumber][n][0];
        ph1=color_coeff_abc[edgeNumber][n][2]*dk*dk +  color_coeff_abc[edgeNumber][n][1]*dk +  color_coeff_abc[edgeNumber][n][0];
        if (TRUNCATE_ROWS && (Math.abs(ph1-ph0)>2.0)) {/* reduce  thisHeight */
          slantDir=(ph1>ph0)?1:-1;
          phLim=ph0+slantDir*2.0*Math.floor(0.5*Math.abs(ph1-ph0));
          while (slantDir*(ph1-phLim) >0) {
            thisHeight--;
            dk=thisHeight-1;
            ph1=color_coeff_abc[edgeNumber][n][2]*dk*dk +  color_coeff_abc[edgeNumber][n][1]*dk +  color_coeff_abc[edgeNumber][n][0];
          }
        }
        color_coeff_abc[edgeNumber][n][5]=thisHeight; /* for statistics */
        for(k=0;k<thisHeight;k++){
          nb=horizontalEdges?rb[n]:n; 
          dk=k+0.5*(nb>>1);
          binShift= (int) (4*(color_coeff_abc[edgeNumber][n][2]*dk*dk +
                              color_coeff_abc[edgeNumber][n][1]*dk +
                              color_coeff_abc[edgeNumber][n][0]+0.5)-2*sSize); /* where is the real center? */
//                              color_coeff_abc[edgeNumber][n][0]) +4.0 -2*sSize);
          hor_min=(binShift)/4;
          hor_max=binShift/4+sSize+4;
          if (hor_min < 0) hor_min=0;
          if (hor_max > selecWidth) hor_max=selecWidth;
          if (hor_min > hor_max) hor_min=hor_max; /// or just fail??
//          if (((hor_min ^ n) & 1) !=0) hor_min++;
          if (((hor_min ^ nb) & 1) !=0) hor_min++;
          FirstQuaterNumber= 0;
          LastQuaterNumber=  0;
          FirstQuaterTotal=0.0;
          LastQuaterTotal= 0.0;
          ThisESFMin=0.0;
          ThisESFAmp=0.0;
          for (hor=hor_min;hor<=hor_max;hor+=2) {
            bin=hor*4-binShift;
            if ((bin>=0) && (bin<=sSize)) {
               FirstQuaterTotal+=ESFArray[2*k+(nb>>1)][hor];
               FirstQuaterNumber++;
            }
            if ((bin>=(sSize4-sSize)) && (bin<=sSize4)) {
               LastQuaterTotal+=ESFArray[2*k+(nb>>1)][hor];
               LastQuaterNumber++;
            }
          }
          if (edgeSign>0) {
            ThisESFMin=LastQuaterTotal/LastQuaterNumber;
            ThisESFAmp=FirstQuaterTotal/FirstQuaterNumber-ThisESFMin;
          } else {
            ThisESFMin=FirstQuaterTotal/FirstQuaterNumber;
            ThisESFAmp=LastQuaterTotal/LastQuaterNumber-ThisESFMin;
          }
          binMinBoth[edgeSign][n]+=ThisESFMin; /* out of bounds here*/
          binAmpBoth[edgeSign][n]+=ThisESFAmp;
          binScaleBothNumber[edgeSign][n]++;
          for (hor=hor_min;hor<=hor_max;hor+=2) {
            bin=hor*4-binShift;
            if ((bin>=0) && (bin<=sSize4)) {
//              binDataBoth[edgeSign][n][bin]  +=ESFArray[2*k+(n>>1)][hor];
              binDataBoth[edgeSign][n][bin]  +=(ESFArray[2*k+(nb>>1)][hor]-ThisESFMin)/ThisESFAmp;
              binNumberBoth[edgeSign][n][bin]++;
            }
          }
        }
      }
    }
    for (edgeSign=0; edgeSign<2; edgeSign++)
      for (n=0;n<4;n++)
        for (bin=0; bin<=sSize4; bin++)
          if (binNumberBoth[edgeSign][n][bin]>0)
            binDataBoth[edgeSign][n][bin]/=binNumberBoth[edgeSign][n][bin];
/**
binDataBoth[edgeSign][n][bin] is now normalized, so average of the first quater is 0.0, last quater 1.0 (or opposite for white-to-black edges)
Next is restoration of the original shifts and scales applied before accumulating bins
*/
    for (edgeSign=0;edgeSign<2;edgeSign++) for (n=0;n<4;n++){
      binMinBoth[edgeSign][n]/=binScaleBothNumber[edgeSign][n];
      binAmpBoth[edgeSign][n]/=binScaleBothNumber[edgeSign][n];
      for (bin=0; bin<=sSize4; bin++)
        binDataBoth[edgeSign][n][bin] = (binDataBoth[edgeSign][n][bin] * binAmpBoth[edgeSign][n]) + binMinBoth[edgeSign][n];

    }
  }


  public double[][][] diffColorVector4Both(double[][][] dataArray, double window){ /**window = 0.0  - as is, 1.0 full width Hamming, 0.5 - middle half  */
    double[][][] diffData   = new double[2][4][sSize*4];
    int i,n;
    double a;
    for (n=0;n<4;n++) for (i=0;i<sSize*4;i++) diffData[0][n][i]=dataArray[0][n][i+1]-dataArray[0][n][i];
    for (n=0;n<4;n++) for (i=0;i<sSize*4;i++) diffData[1][n][i]=dataArray[1][n][i]-dataArray[1][n][i+1];
    if (window>0.0)  for (n=0;n<4;n++) {
      for (i=0;i<sSize*4;i++) {
        a=(i-2*sSize)/(2*sSize*window);
        if ((a< -1.0) || (a>1.0)) {
          diffData[0][n][i]=0.0;
          diffData[1][n][i]=0.0;
        } else {
          a=(0.54+0.46*Math.cos(Math.PI*a));
          diffData[0][n][i]*=a;
          diffData[1][n][i]*=a;
        }
      }
   }
   return diffData;
  }



/* ============================== */
  private double[][] normalizePlots (double[][] Vectors, boolean both ) {
     double [][] result= new double [Vectors.length][Vectors[0].length];
     int n,i;
     double min,max, scale;
     for (n=0;n<Vectors.length; n++) {
       min=Vectors[n][0];
       max=min;
       for (i=0;i<Vectors[n].length;i++) {
         if (Vectors[n][i]<min) min=Vectors[n][i];
         if (Vectors[n][i]>max) max=Vectors[n][i];
       }
       if (!both) min=0;
       if (max<=min) scale =1.0;
       else scale = 1.0/(max-min);
       for (i=0;i<Vectors[n].length;i++) result[n][i]=scale*(Vectors[n][i]-min);
     }
     return result;
  }
  public void generatePlots(double[][] Vectors, String plotType, String plotName, Color [] colors, String headers) {
    generatePlots(Vectors, plotType, plotName, colors, true, headers);
  }
  public void generatePlots(double[][] Vectors, String plotType, String plotName, Color [] colors, boolean asText) {
    generatePlots(Vectors, plotType, plotName, colors, asText, "");
 
  }
  public void generatePlots(double[][] Vectors, String plotType, String plotName, Color [] colors, boolean asText, String headers) {
    generatePlots(Vectors, plotType, plotName, colors, asText, headers, 0,0);
  }
  public void generatePlots(double[][] Vectors, String plotType, String plotName, Color [] colors, boolean asText, String headers, int textWidth, int TextHeight) {
    double[]xValues;
    String ejeX="pixel";
    String ejeY="";
    String allTitle="";
    boolean scaled=false;
    double[][] finalVectors;
///        String plotType=plot.substring(0,3);
//    ImageProcessor imgProc;
    double thisYmin, thisYmax,thisXmin,thisXmax;
    int i,n;
    if (colors == null) colors=bayer_colors;
    //If MTF plot, calculate the scale of cycles per pixel for x-axis values
    if (plotType.equals("MTF")){
      ejeY="Modulation Factor";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],plotType);
    } else if (plotType.equals("PHASE")){
      ejeY="Phase";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
    } else if (plotType.equals("AMPLITUDE")){
      ejeY="Amplitude";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
    } else if (plotType.equals("RE")){
      ejeY="Re";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
    } else if (plotType.equals("REIM")){
      ejeY="Re/Im";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
    } else if (plotType.equals("REIMS")){
      ejeY="Re/Im (scaled)";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
    } else if (plotType.equals("RES")){
      ejeY="Re (scaled)";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
      scaled=true;
    } else if (plotType.equals("IM")){
      ejeY="Im";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
    } else if (plotType.equals("IMS")){
      ejeY="Im (scaled)";
      ejeX="cycles / pixel";
      xValues=calculateXValues(Vectors[0],"MTF");
      scaled=true;
    } else if (plotType.equals("EDG")){
      xValues=calculateXValues(Vectors[0],"EDG");
    } else {
      xValues=calculateXValues(Vectors[0],"DV4"); // divide by 4 - pixels were subdivided into bins
    }

    //plot titles         
    if (plotType.equals("ESF")){
      ejeY="Grey Value";
    }
    if (plotType.equals("LSF")){ 
      ejeY="Grey Value / pixel";
    }
    if (plotType.equals("MTF")){
      ejeY="Modulation Factor";
      ejeX="cycles / pixel";
    }
    finalVectors=(scaled)?scaleReIm(Vectors):Vectors;

    allTitle=plotName + "_" + title_src;        
    //plot limits
    thisYmin=finalVectors[0][0];
    thisYmax=finalVectors[0][0];
    thisXmin=xValues[0];
    thisXmax=xValues[0];
    for (i=0; i<Vectors[0].length; i++) {
      if      (xValues[i] < thisXmin) thisXmin = xValues[i];
      else if (xValues[i] > thisXmax) thisXmax = xValues[i];
      for (n=0;n<Vectors.length; n++) {
        if      (finalVectors[n][i] < thisYmin) thisYmin = finalVectors[n][i];
        else if (finalVectors[n][i] > thisYmax) thisYmax = finalVectors[n][i];
      }
    }
    if (asText) {
      int headerLength=0;
      for (i=headers.indexOf("\t");i>=0;i=headers.indexOf("\t",i+1))  headerLength++;
//      IJ.showMessage(headers,"headerLength="+headerLength+" i="+i+" finalVectors.length="+finalVectors.length);
      if ((headerLength==(finalVectors.length>>1)) && USE_COMPLEX)       createComplexCSV(allTitle, xValues, finalVectors, headers, textWidth, TextHeight);
      else createCSV(allTitle, xValues, finalVectors, headers, textWidth, TextHeight);
      return;
    }
    plotResult = new Plot(allTitle, ejeX, ejeY, (double []) null, (double []) null);
    if (plotType.equals("ESF")){
      plotResult.setLimits(1,finalVectors[0].length,0,thisYmax);
    } else if (plotType.equals("LSF")){
      plotResult.setLimits(1,finalVectors[0].length,0,thisYmax);
    } else if (plotType.equals("MTF")){
      plotResult.setLimits(0,0.5,0,1);
    } else if (plotType.equals("PHASE")){
      plotResult.setLimits(0,0.5,-Math.PI,Math.PI);
    } else if (plotType.equals("AMPLITUDE")){
      plotResult.setLimits(0,0.5,0,thisYmax);
    } else if ((plotType.equals("RE")) || (plotType.equals("IM")) || (plotType.equals("REIM"))){
      plotResult.setLimits(0,0.5,thisYmin,thisYmax);
    } else {
      plotResult.setLimits(thisXmin,thisXmax,thisYmin,thisYmax);
    }
    for (n=0; n<finalVectors.length; n++) {
      plotResult.setColor(colors[n]);
      plotResult.addPoints(xValues, finalVectors[n], Plot.LINE);
    }
    plotResult.draw();
    plotResult.show();
  }

  double [][] scaleReIm (double [][] data) { /* Scale each plot, preserving 0 */
    double [][] result = new double [data.length][data[0].length];
    double mx;
    int i,n;
    for (n=0;n<data.length;n++) {
      mx=Math.abs(data[n][0]);
      for (i=1;i<data[0].length;i++) if (Math.abs(data[n][i]) > mx ) mx=Math.abs(data[n][i]);
      if (mx>0) mx=1.0/mx;
      for (i=0;i<data[0].length;i++) result[n][i]=mx*data[n][i];
    }
    return result;
  }
 
  private void createCSV(String title, double [] xValues, double[][] yValues, String headers, int Width, int Height) {
    if ((xValues==null) || (xValues.length==0)) return;
    if ((yValues==null) || (yValues.length==0)) return;

    String hs=headers;
    int i,n;
    if (hs.equals("")) {
       hs=new String();
       hs+="X";
       for (i=0;i<yValues.length; i++)  hs+="\t"+(i+1);
    }
    StringBuffer sb = new StringBuffer();
    for (i=0;i<xValues.length;i++)  {
       sb.append(IJ.d2s(xValues[i],2));
       for (n=0;n<yValues.length;n++)  sb.append("\t"+IJ.d2s(yValues[n][i],4));
       sb.append("\n");
   }
   new TextWindow(title+"_(csv)", hs, sb.toString(), (Width>0)?Width:(84*(yValues.length+1)), (Height>0)?Height:600);
  }

  private void createComplexCSV(String title, double [] xValues, double[][] yValues, String headers, int Width, int Height) {
    if ((xValues==null) || (xValues.length==0)) return;
    if ((yValues==null) || (yValues.length==0)) return;

    String hs=headers;
    int i,n;
    if (hs.equals("")) {
       hs=new String();
       hs+="X";
       for (i=0;i<(yValues.length>>1); i++)  hs+="\t"+(i+1);
    }
    StringBuffer sb = new StringBuffer();
    for (i=0;i<xValues.length;i++)  {
       sb.append(IJ.d2s(xValues[i],2));
       for (n=0;n<yValues.length;n+=2)  sb.append("\t"+IJ.d2s(yValues[n][i],4)+((yValues[n+1][i]>=0)?"+":"")+IJ.d2s(yValues[n+1][i],4)+"i");
       sb.append("\n");
   }
   new TextWindow(title+"_(csv)", hs, sb.toString(), (Width>0)?Width:(84*(yValues.length+1)), (Height>0)?Height:600);
  }



  public double[] calculateXValues(double[] Vector, String plot){
    int i;
    int N=Vector.length;
    double[]xValues = new double[N];
    if(plot.substring(0,3).equals("MTF")){
      xValues[0]=0;
      //Scale of values for x-axis
      for(i=1;i<N;i++) xValues[i]=xValues[i-1]+(0.5/(N-1));
    } else if(plot.substring(0,3).equals("DV4")){
      for(i=0;i<N;i++) xValues[i]=0.25*(i+1); /* why is it from 1, not from 0 ? */
    } else { 
      for(i=0;i<N;i++) xValues[i]=i+1; /* why is it from 1, not from 0 ? */
    }
    return xValues;
  }




/* ============================== */

  public class Complex {
    private final double re;   // the real part
    private final double im;   // the imaginary part
    // create a copy olf a complex
    public Complex(Complex c) {
      this.re = c.Re();
      this.im = c.Im();
    }
    // create a new object with the given real and imaginary parts
    public Complex(double real, double imag) {
      this.re = real;
      this.im = imag;
    }
    // return a string representation of the invoking object
    public String toString()  {return re + " + " + im + "i"; }

    // return a new object whose value is (this - b)
    public Complex plus(Complex b) { 
      Complex a = this;
      double real = a.re + b.re;
      double imag = a.im + b.im;
      Complex rslt = new Complex(real, imag);
      return rslt;
    }

    // return a new object whose value is (this - b)
    public Complex minus(Complex b) { 
      Complex a = this;   
      double real = a.re - b.re;
      double imag = a.im - b.im;
      Complex rslt = new Complex(real, imag);
      return rslt;
    }        
    // return a new object whose value is (this * b)
    public Complex times(Complex b) {
      Complex a = this;
      double real = a.re * b.re - a.im * b.im;
      double imag = a.re * b.im + a.im * b.re;
      Complex rslt = new Complex(real, imag);
      return rslt;
    }
    public Complex scale(double b) {
      Complex a = this;
      Complex rslt = new Complex(a.re *b, a.im *b );
      return rslt;
    }
    // return a new object whose value is (this / b)
    public Complex divideBy(Complex b) {
      Complex a = this;
      double div=b.re*b.re+b.im*b.im;
      double real = (a.im * b.im + a.re * b.re)/div;
      double imag = (a.im * b.re - a.re * b.im)/div;
      Complex rslt = new Complex(real, imag);
      return rslt;
    }
      // return |this|
      public double abs() { return Math.sqrt(re*re + im*im);  }
      public double phase() { return Math.atan2(im,re); }
      public double Re() {return this.re; }
      public double Im() {return this.im; }
    }    

    public Complex[] fft(Complex[] x) {
      int N = x.length;
      Complex[] y = new Complex[N];
    // base case
      if (N == 1) {
        y[0] = x[0];
        return y;
      }
     // radix 2 Cooley-Tukey FFT
      if (N % 2 != 0) throw new RuntimeException("N is not a power of 2");
      Complex[] even = new Complex[N/2];
      Complex[] odd  = new Complex[N/2];
      for (int k = 0; k < N/2; k++) even[k] = x[2*k];
      for (int k = 0; k < N/2; k++) odd[k]  = x[2*k + 1];
      Complex[] q = fft(even);
      Complex[] r = fft(odd);
      for (int k = 0; k < N/2; k++) {
        double kth = -2 * k * Math.PI / N;
        Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
        y[k]       = q[k].plus(wk.times(r[k]));
        y[k + N/2] = q[k].minus(wk.times(r[k]));
      }
      return y;
    }
  public void showOTFHelp() {

  new TextWindow(title_src+"_"+"pixel crosstalk", "Item                       \tDescription",
  "LICENSE"+
  "\tMTF_Bayer.java is free software: you can redistribute it and/or modify\n"+
  "\tit under the terms of the GNU General Public License as published by\n"+
  "\tthe Free Software Foundation, either version 3 of the License, or\n"+
  "\t(at your option) any later version.\n"+
  "\t\n"+
  "\tThis program is distributed in the hope that it will be useful,\n"+
  "\tbut WITHOUT ANY WARRANTY; without even the implied warranty of\n"+
  "\tMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"+
  "\tGNU General Public License for more details.\n"+
  "\t\n"+
  "\tYou should have received a copy of the GNU General Public License\n"+
  "\talong with this program.  If not, see <http://www.gnu.org/licenses/>.\n"+

  "\t\n"+
  "\t\n"+
  "Program\tThis program is designed to calculate optical transfer function (OTF)  using slanted edge (SE) method\n"+
  "\tbased on (but not exactly implemented) ISO standard 12233, it reuses parts of the code of SE_MTF plugin\n"+
  "\tfor ImageJ (http://rsbweb.nih.gov/ij/plugins/se-mtf/index.html)\n"+
  "\t\n"+
  "\tThis program is designed to work only with the images acquired from the color sensors with Bayer mosaic\n"+
  "\tfilters, the images should be presented as monochrome pixel array with repetitive 2x2 pixel pattern of GR/BG\n"+
  "\tThat means that all pixels in even rows, even columns (starting from 0/0) correspond to green sensor pixels\n"+
  "\t(referenced in the program as Gr - Green in \"red\" row), pixels in even rows/odd columns correspond to red (R) ones\n"+
  "\tpixels in odd rows/even columns - blue (B) and odd rows/odd columns - green (Gb - green in \"blue\" row)."+
  "\tImages that have different Bayer patterns can be converted to GR/BG by shifting by one pixel horizontally\n"+
  "\tand/or vertically by ImageJ.\n"+
  "\t\n"+
  "\tAdditionally images are assumed to be linear, and to have the same signal gain in each of the color channels,\n"+
  "\tyou may use JP46_reader plugin that can process jp46 format images from Elphel cameras and un-apply gamma and\n"+
  "\tchannel gain settings (JPEG quality is recommended to be set to exactly 100%, gamma to be set to 0.5, color\n"+
  "\tgains adjusted for approximate white balance (if possible), having both green channel gains set exactly\n"+
  "\tto the same value.\n"+
  "\t\n"+
  "\tThe program is run when you select the desired area in a valid image window with rectangular selection tool of ImageJ\n"+
  "\tand press the \"Process\" button. The program detects a number of near-vertical or near horizontal edges that should\n"+
  "\trun all the width (or height) of the selection, approximates them with second-degree polynomials and then accumulates\n"+
  "\tinto quarter-pixel bins (similarly to ISO 12233), separately for each edge polarity (black-to-white and white-to-black)\n"+
  "\tand separately for each Bayer component (two greens are treated individually)\n"+
  "\t\n"+
  "\tThis processing depends on the settings that can be customized in \"Configure\" menu, you may select which of the\n"+
  "\tplots and/or tables (exportable to spreadsheet programs) to be output each time the processing takes place.\n"+
  "\t\n"+
  "\tAdditionally you may use \"Show ...\" buttons, they display various plots/tables for the data being last processed,\n"+
  "\tthey will apply to the same image selection and parameters until you explicitly press \"Process\" again.\n"+
  "\t\n"+
  "Configure"+
  "\tSet parameters that modify image data processing\n"+
  "\t\n"+
  "Conversion strip"+
  "\tFull width of the pixel row across the edges that are being processed. Should be a power of 2.\n"+
  "width\t\n"+
  "\t\n"+
  "Monochrome"+
  "\tThis box should be checked for the images acquired with near-monochromatic light, when the chromatic aberrations\n"+
  "(color filtered)"+
  "\tare negligible and the edge is approximated for all color component pixels together. Without the color filtering \n"+
  "image"+
  "\tlens lateral chromatic aberrations may require different edge approximation coefficients for each color separately.\n"+
  "\t\n"+
  "Edge location"+
  "\t This selection box switches between different algorithms that are used to detect edge location in each processed\n"+
  "algorithm"+
  "\trow with sub-pixel accuracy. This selection is also modified by the \"LSF calculation window\" explained below.\n"+
  "\t\n"+
  "Centroid"+
  "\tCalculate centroid of each line LSF obtained as a digital derivative, similarly to ISO 12233 recommendation.\n"+
  "\tIf the \"LSF calculation window\" is set to a non-zero value, the centroid calculation takes place in two steps:\n"+
  "\tfirst all the data is used, then the LSF data is multiplied by the Hamming window function around the first iteration\n"+
  "\tposition, and centroid is re-calculated for the LSF modified by window function. This second step is designed to\n"+
  "\treduce the influence of the far from the edge pixels in the edge locating, that also helps to reduce the skew caused\n"+
  "\tby non-uniform illumination of the target\n"+
  "\t\n"+
  "Centroid with"+
  "\tSimilarly to \"Centroid\" but with additional iteration step with a window function twice narrower than the selected\n"+
  "double windowing"+
  "\tone. So first centroid is calculated for all the data, second - multiplied by the Hamming of specified width and\n"+
  "\tthe third one uses half-width window.\n"+
  "\t\n"+
  "Two ESF"+
  "\tThe edge threshold value is calculated as the average value of the first and last quarter pixels (away from the edge).\n"+
  "crossing points"+
  "\tThen the edge location is determined by linear approximation between the two points immediately above and below the\n"+
  "\tthreshold. In this method the far points have very little influence on the edge position\n"+
  "\t\n"+
  "Four ESF"+
  "\tSimilar to the \"Two ESF crossing points\", but with the cubic interpolation between the 4 points - two from each\n"+
  "crossing points"+
  "\tside of the threshold crossing. Cubic equation is solved approximately using two linear iterations.\n"+
  "\t\n"+
  "Truncate rows"+
  "\tWhen this box is checked, the program tries to follow ISO 12233 by truncating the bottom rows from the processing\n"+
  "to equalize"+
  "\tso the edges cross multiple of pixel columns (in this program for Bayer pixels it is multiple of two-pixel columns)\n"+
  "phase rotations"+
  "\tso during binning each phase gets fair representation. Actually in this program this has rather small effect because\n"+
  "\tof an additional deviation from the ISO algorithm - each row ESF data is normalized (with normalization coefficient\n"+
  "\tsaved) so the first and last quarter pixels average value map to 0.0 and 1.0, respectively. That reduces artifacts\n"+
  "\tcaused by different number of samples accumulated in each bin when row data is subject to non-uniform illumination.\n"+
  "\t\n"+
  "LSF calculation"+
  "\tWindow function (Hamming) width, applied to the LSF before calculating Discrete Fourier Transform to determine OTF.\n"+
  "window"+
  "\tWindow width is entered as a fraction of the \"Conversion strip width\", value 0.0 turns windowing off. This wimndow\n"+
  "\talso applies to centroid calculations if selected in \"Edge location algorithm\".\n"+
  "\t\n"+
  "Show ..."+
  "\tThe following check boxes determine which of the plots/tables are displayed automatically when you press \"Process\".\n"+
  "\tThose plots/tables can be shown later by pressing \"Show ...\" buttons. \"Show OTF phase\" modifies \"Show MTF\"\n"+
  "\tbox and button, so each time two plots are shown - MTF (OTF amplitude) and OTF phase.\n"+
  "\t\n"+
  "Show"+
  "\tDisplay edge approximation coefficients, RMS error for each edge, as well as total number of rows used for ESF/LSF/OTF\n"+
  "Approximation"+
  "\tmeasurements. Worst of each edge RMS error is included in the window title bar.\n"+
  "\t\n"+
  "Gb/Gr sensitivity"+
  "\tdebugging feature, allowing manual input for entering difference in the two green channels analog gain. Rather small\n"+
  "correction"+
  "\tunequity of these gains can influence crosstalk measurements\n"+
  "\t\n"+
  "Debug level"+
  "\tHigher integer numbers increase amount of text output. Currently the highest value is 3\n"+
  "",
  1000, 800);

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