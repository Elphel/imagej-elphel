package com.elphel.imagej.calibration;
/**
** -----------------------------------------------------------------------------**
** target_points.java
**
** Measures focus sharpnes at orthogonal directions, differenct colors,
** Displays results for manual focusing/image plane tilting.
** NOTE: Requires special targets !
**
** Copyright (C) 2010 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  target_points.java is free software: you can redistribute it and/or modify
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

import java.util.List;

import com.elphel.imagej.jp4.JP46_Reader_camera;

import java.util.ArrayList;

import ij.text.*;

import java.lang.Integer;


public class target_points extends PlugInFrame implements ActionListener {
private static final long serialVersionUID = -3057496866952930812L;
JP46_Reader_camera jp4_instance;
//  MTF_Bayer MTF_Bayer_instance;
  Panel panel;
  static Frame instance;

 public static int DEBUG_LEVEL = 1;
 public static int MASTER_DEBUG_LEVEL = 1;

 public static int FFTSize=64;
 public static int FFTScanStep=8;
 
 public static int test_x=FFTSize;
 public static int test_y=FFTSize;
 public static int     displayWidth=800;
 public static int     displayHeight=600;
 public static float [][] input_bayer=null;
 public static float [][] convolved_bayer=null;
 public static float [][] normalized_convolved_bayer=null;
 
 public static float   [] target_kernel=null;
 public static int   [][] clusterMaps=null; 
 public static double [][][] targetCoordinates;

 public static double targetOuterDMin =38; // minimal outer diameter of the target image , in pixels
 public static double targetOuterDMax =47; // maximal outer diameter of the target image , in pixels
 public static int    numTargetRings  = 2;  // number of target black rings (notg counting center black circle)
 public static int    pixelsSubdivide =10;  // subdivide pixels by this number (each direction) when generating targets
 public static double deconvInvert =   0.1; // when FFT component is lass than this fraction of the maximal value, replace 1/z with Z
 public static double filteredRadius=   7.0;  //pix - search maximums after convolution in (2*filteredRadius+1) squares
 public static double backgroundRadius= 15.0; //pix - consider ring area between backgroundRadius and filteredRadius as reference
 public static double clusterThreshold= 1.0; //1.2; // tested with 0.8 - many extras, but filtered out
 public static int    clusterSize=      20;  // cluster size (will be expanded/shrank before finding centroid

 public static int    discrAngularFreq=     2 ;  // pixels on FFT image of tragets converted polar (the smaller, the less angular variations)

/// these parameters are dependent on targets, use debug mode and manula fft for 64x64 polar coordinates target areas
 public static int    discrRadialMinFreq=  7 ; // pixels on FFT image of targets converted polar (radial component)
 public static int    discrRadialMaxFreq=  9 ; // pixels on FFT image of targets converted polar (radial component)
 public static double discrThreshold=      0.1; // FFT energy fraction in selecter area should be > this threshold to pass the test
 public static double maxChromaticDistance= 10.0; // Maximal distance between the same target on different color copmponents


 public static double[][][] targets; // For each target {averageX, averageY, num_non_zero_components},{X1,Y1,Qulaity1},...,{X4,Y4,Qulaity4}
/**
discrRadialMinFreq=(size*(2* numTargetRings +1)/targetOuterDMax)-1
discrRadialMaxFreq=(size*(2* numTargetRings +1)/targetOuterDMin)+1
*/

 private ImagePlus     imp_src;
 public ImageProcessor ip_display;
 public ImagePlus      imp_display;
 public ImagePlus      imp_camera=null;
 
  Plot plotResult;
 

  public target_points() {
    super("target_points");
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
    addButton("Split Bayer");
    addButton("Create Target");
    addButton("Split&Convolve");
    add(panel);
    pack();
    GUI.center(this);
    setVisible(true);
    initHamming(FFTSize);
//    initDisplay();
    jp4_instance=       new JP46_Reader_camera();
  }
  void addButton(String label) {
    Button b = new Button(label);
    b.addActionListener(this);
    b.addKeyListener(IJ.getInstance());
    panel.add(b);
  }
  public void actionPerformed(ActionEvent e) {
    int i,j,ir2, size;
    String label = e.getActionCommand();
    if (label==null) return;
    if (label.equals("Configure")) {
      if (showDialog()) {
        initHamming(FFTSize);
      }
      return;
    }
    if (label.equals("Create Target")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      target_kernel=createTargetDialog ();
      return;
    }
    imp_src = WindowManager.getCurrentImage();
    String newTitle= imp_src.getTitle();
    Rectangle r=new Rectangle(imp_src.getWidth(),imp_src.getHeight());

    if (label.equals("Split Bayer")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      input_bayer=splitBayer (imp_src);
      showBayers(input_bayer, r.width>>1, r.height>>1, newTitle);
      return;
    }
    if (label.equals("Split&Convolve")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      input_bayer=splitBayer (imp_src);
      if (DEBUG_LEVEL>5) showBayers(input_bayer, r.width>>1, r.height>>1, newTitle);
       
      if (target_kernel==null) {
        IJ.showMessage("Error","Target kernel does not exist, please generate one");
        return;
      }
      size=(int) Math.sqrt(target_kernel.length);
//target_kernel
/* Convolve Bayer slices with prepared target template  */
      convolved_bayer=new float[input_bayer.length][];
      for (i=0; i<input_bayer.length; i++) {
        IJ.showStatus("Convolving Bayer "+i);
/* Double in convolution works twice faster than float!*/
        convolved_bayer[i]=doubleConvolveWithTarget(input_bayer[i], target_kernel, r.width>>1, r.height>>1, size);

      }
      if (DEBUG_LEVEL>2) showBayers(convolved_bayer, r.width>>1, r.height>>1, newTitle+"_"+deconvInvert);
//
//    filteredRadius   =(int) gd.getNextNumber();
//    backgroundRadius =(int) gd.getNextNumber();
/* normalize convolved Bayer slices   */
/* prepare pixel mask for the normalization (ring) */
      int filtSize= (int) filteredRadius;
      boolean [][] mask=new boolean[2*filtSize+1][2*filtSize+1];
      for (i=0;i<(2*filtSize+1);i++) for (j=0;j<(2*filtSize+1);j++) {
        ir2=(i-filtSize)*(i-filtSize) + (j-filtSize)*(j-filtSize);
        mask[i][j]=(ir2>(filteredRadius*filteredRadius)) && (ir2 < (backgroundRadius*backgroundRadius));
      }
// public static float [][] normalized_convolved_bayer=null;
      normalized_convolved_bayer=new float[input_bayer.length][];
      for (i=0; i<input_bayer.length; i++) {
        IJ.showStatus("Normalizing Bayer "+i);
        normalized_convolved_bayer[i]=normalizeAtRing(convolved_bayer[i], r.width>>1, r.height>>1, mask);
      }
      if (DEBUG_LEVEL>1) showBayers(normalized_convolved_bayer, r.width>>1, r.height>>1, newTitle+"_"+deconvInvert+"_normalized");
// public static int   [][] clusterMaps=null; 
//clusterThreshold
      clusterMaps= new int[4][];
      targetCoordinates=new double[4][][];
      for (i=0; i<input_bayer.length; i++) {
        IJ.showStatus("Clusterizing Bayer "+i);
//        clusterMaps[i]=clusteriseTargets(normalized_convolved_bayer[i], r.width>>1, r.height>>1,clusterThreshold,clusterSize);
        targetCoordinates[i]=clusteriseTargets(normalized_convolved_bayer[i], r.width>>1, r.height>>1,clusterThreshold,clusterSize);
      }
      float[][]clusterPixels=new float[4][input_bayer[0].length];
      if (DEBUG_LEVEL>2) {
//         float[][]clusterPixels=new float[4][input_bayer[0].length];
         for (i=0; i<clusterPixels.length; i++) {
          for (j=0; j<clusterPixels[0].length; j++) clusterPixels[i][j]=0.0f;
          for (j=0;j<targetCoordinates[i].length;j++) {
             clusterPixels[i][((int)(Math.round(targetCoordinates[i][j][1])*(r.width>>1))) +((int)Math.round(targetCoordinates[i][j][0]))]=1.0f;
          }
        }
        showBayers(clusterPixels, r.width>>1, r.height>>1, newTitle+"_"+deconvInvert+"_clusters");
      }
      float [] rectTarget;
      ImageProcessor ip_dbg;
      ImagePlus      imp_dbg;
      double          likely;
      double [] likelyness;
      int numGoodTargets;
      double [][] goodTargets; /// x/y/above Threshold
      double r0=(targetOuterDMin+targetOuterDMax)/(2* numTargetRings +1)/8; /// copmpensate for the center circle twice wider than rings
/**
discrRadialMinFreq=(size*(2* numTargetRings +1)/targetOuterDMax)-1
discrRadialMaxFreq=(size*(2* numTargetRings +1)/targetOuterDMin)+1
*/
      for (i=0; i<clusterPixels.length; i++) {
        likelyness= new double[targetCoordinates[i].length];
        numGoodTargets=0;
        for (j=0;j<targetCoordinates[i].length;j++) {
          rectTarget=circle2DoubleRect (input_bayer[i], r.width>>1, r.height>>1, size,targetCoordinates[i][j][0],targetCoordinates[i][j][1], r0);
          likely=likelyTarget(rectTarget, size,discrAngularFreq,discrRadialMinFreq,discrRadialMaxFreq);
          likelyness[j]=likely;
          if (DEBUG_LEVEL>2) {
   System.out.println("Cluster="+j+" x="+targetCoordinates[i][j][0]+" y="+targetCoordinates[i][j][1]+" likely="+likely);

          }
          if (DEBUG_LEVEL>3) {
            if (i==0) {  // just to reduce debug clutter
              ip_dbg=new FloatProcessor(size,size);
              ip_dbg.setPixels(rectTarget);
              ip_dbg.resetMinAndMax();
              imp_dbg=  new ImagePlus(newTitle+"_rect_"+j, ip_dbg);
              imp_dbg.show();
            }
          }
          if (likely >=discrThreshold) numGoodTargets++;
        }
        goodTargets=new double[numGoodTargets][3];
        numGoodTargets=0;
        for (j=0;j<targetCoordinates[i].length;j++) {
          if (likelyness[j]>=discrThreshold) {
            goodTargets[numGoodTargets][2]=likelyness[j]/discrThreshold;
            goodTargets[numGoodTargets][0]=2*targetCoordinates[i][j][0]+(i&1);      // convert to full image, use Bayer shift
            goodTargets[numGoodTargets][1]=2*targetCoordinates[i][j][1]+((i>>1)&1); // convert to full image, use Bayer shift
            numGoodTargets++;
          } 
        }
        targetCoordinates[i]=goodTargets;
      }

      if (DEBUG_LEVEL>1) {
         for (i=0; i<clusterPixels.length; i++) {
          for (j=0; j<clusterPixels[0].length; j++) clusterPixels[i][j]=0.0f;
          for (j=0;j<targetCoordinates[i].length;j++) {
             clusterPixels[i][((int)(Math.round((targetCoordinates[i][j][1] - ((i>>1)&1))/2) *(r.width>>1))) +
                              ((int)(Math.round((targetCoordinates[i][j][0] - ( i    &1))/2)))]=1.0f;
   System.out.println("Bayer="+i+" Target="+(j+1)+" x="+targetCoordinates[i][j][0]+" y="+targetCoordinates[i][j][1]+" quality="+targetCoordinates[i][j][2]);
          }
        }
        showBayers(clusterPixels, r.width>>1, r.height>>1, newTitle+"_"+deconvInvert+"_clusters");
      }
/* TODO: Verify that all Bayer components have the same targets (build composite table) */
      targets= combineTargets(targetCoordinates, maxChromaticDistance);
      showTargetsLocationTable(targets,  newTitle, 2,  (DEBUG_LEVEL>1));
      return;
    }
  }

/* Combine target locations from 4 Bayer components */
  double [][][] combineTargets(double[][][] targetCoord, ///[bayer][number][x,y,q>1]
                               double maxDistance) {     /// maximal distance between the same target in different Bayer components

   int [][] targetNumbers=new int[targetCoord.length][];
   int i,i1,j,k,n,l,maxLen, bNum, numTargets;
   double d2=maxDistance*maxDistance;
   double [][][] targets;
   double avX,avY;
   maxLen=0;
   bNum=-1;
   for (i=0;i<targetCoord.length;i++) {
      l=targetCoord[i].length;
      targetNumbers[i]=new int [l];
      for (j=0;j<l;j++) targetNumbers[i][j]=0;
      if (maxLen<l){
        maxLen=l;
        bNum=i;
      }
   }
 /* assign target number according to the component that has most of the targets (does not mean others do not have that this one is missing */
   for (j=0;j<maxLen;j++) targetNumbers[bNum][j]=j+1;
   numTargets=maxLen; // may increase later
 /* compare all other color components with the coordinates in the seslected one (not too many to bother with good guess) */
   for (i=0;i<targetNumbers.length;i++) if (i!=bNum) {
     for (j=0;j<targetNumbers[i].length;j++) if (targetNumbers[i][j]==0){
        for (k=0;k<targetNumbers[bNum].length;k++) {
          if (((targetCoord[bNum][k][0]-targetCoord[i][j][0])*(targetCoord[bNum][k][0]-targetCoord[i][j][0])+
               (targetCoord[bNum][k][1]-targetCoord[i][j][1])*(targetCoord[bNum][k][1]-targetCoord[i][j][1]))<=d2) {
              targetNumbers[i][j]=targetNumbers[bNum][k];
            break;
          }
        }
     }
   }
/* See if any targets are missing, add them */
   for (i=0;i<targetNumbers.length;i++) if (i!=bNum) {
     for (j=0;j<targetNumbers[i].length;j++) if (targetNumbers[i][j]==0){
       numTargets++;
       targetNumbers[i][j]=numTargets;
       for (i1=i+1;i1<targetNumbers.length;i1++) if (i1!=bNum) {
         for (k=0;k<targetNumbers[i1].length;k++) if (targetNumbers[i1][k]==0) {
           targetNumbers[i1][k]=numTargets;
         }
       }
     }
   }
   if (DEBUG_LEVEL>2) {
     System.out.println("numTargets="+numTargets);
     for (i=0;i<targetCoord.length;i++) for (j=0;j<targetCoord[i].length;j++) {
       System.out.println("["+targetNumbers[i][j]+"] "+i+":"+j+" "+targetCoord[i][j][0]+","+targetCoord[i][j][1]+" :"+targetCoord[i][j][2]);
     }
   }




   targets = new double [numTargets][targetNumbers.length+1][3];
   for (i=0;i<numTargets;i++) for (j=0;j<targets[i].length;j++) for (k=0;k<3;k++) targets [i][j][k]=0.0;
   for (i=0;i<targetNumbers.length;i++) for (j=0;j<targetNumbers[i].length;j++) if (targetNumbers[i][j]!=0){ // should be always non-zero
     k=targetNumbers[i][j]-1;
     targets[k][i+1][0]=targetCoord[i][j][0]; // x
     targets[k][i+1][1]=targetCoord[i][j][1]; // y
     targets[k][i+1][2]=targetCoord[i][j][2]; // quality >1.0
   }
/* Calculate average values*/
   for (i=0;i<numTargets;i++) {
     avX=0.0;
     avY=0.0;
     n=0;
     for (j=1;j< targets[i].length; j++) if (targets[i][j][2]>0){
       avX+=targets[i][j][0];
       avY+=targets[i][j][1];
       n++;
     }
     if (n>0) {
       targets[i][0][0]=avX/n;
       targets[i][0][1]=avY/n;
       targets[i][0][2]=n;
     }
   }
   if (DEBUG_LEVEL>2) {
     System.out.println("targets");
     for (i=0;i<targets.length;i++) {
       System.out.println(i+" | "+targets[i][0][0]+","+targets[i][0][1]+" :"+targets[i][0][2]+
                            " | "+targets[i][1][0]+","+targets[i][1][1]+" :"+targets[i][1][2]+
                            " | "+targets[i][2][0]+","+targets[i][2][1]+" :"+targets[i][2][2]+
                            " | "+targets[i][3][0]+","+targets[i][3][1]+" :"+targets[i][3][2]+
                            " | "+targets[i][4][0]+","+targets[i][4][1]+" :"+targets[i][4][2]);
     }
   }

   return targets;
  }

  public void showTargetsLocationTable(double[][][] targets, String title, int precision, boolean showQuality) {
    int i,n;
    
    String header="#\tX\tY";
    for (i=1;i<targets[0].length;i++) header+="\tdX"+i+"\tdY"+i+(showQuality?("\tQ"+i):"");
    StringBuffer sb = new StringBuffer();
    for (n=0;n<targets.length;n++) {
      sb.append((n+1)+
                  "\t"+IJ.d2s(targets[n][0][0],precision)+  // Average X 
                  "\t"+IJ.d2s(targets[n][0][1],precision)); // Average Y 
      for (i=1;i<targets[0].length;i++) {
        if (targets[n][i][2]>0) {
          sb.append(  "\t"+(((targets[n][i][0]-targets[n][0][0])>0)?"+":"")+IJ.d2s(targets[n][i][0]-targets[n][0][0],precision)+  // X
                      "\t"+(((targets[n][i][1]-targets[n][0][1])>0)?"+":"")+IJ.d2s(targets[n][i][1]-targets[n][0][1],precision)); // Y
        } else {
          sb.append(  "\t---\t---");
        }
        if (showQuality) sb.append("\t"+((targets[n][0][2]>0)?IJ.d2s(targets[n][i][2],precision):"---"));
      }
      sb.append(  "\n");
    }
    new TextWindow(title+"_Target_Locations_Table", header, sb.toString(), showQuality?900:700,500);
  }



/* determines if it was likely a target of concentric circles, after convertion to polar coordinates expect nearly horizontal b/w bands */
double likelyTarget(float[] pixels, // pixel array
                          int size, // image size (should be square for FFT
                          int hor,  // horizontal selection area (half width)
                          int vertMin,
                          int vertMax) // vertical selection area (half height > half width for horizointal bands)
 {

      ImageProcessor ip;
      FHT            fht;
      double[][][]   fft;
      double         s1,s2,e;
      int i,j;
      ip=new FloatProcessor(size,size);
      ip.setPixels(pixels);
      fht =  new FHT(ip);
// Swapping quadrants, so the center will be 0,0
      fht.swapQuadrants();
// get to frequency
      fht.transform();
// Convert from FHT to complex FFT
      fft= FHT2FFTHalf (fht);
      s1=0;
      s2=0;
     
      for (i=0;i<(size/2+1);i++) for (j=0;j<size;j++) {
       if ((i>0) || (j>0)) { // skip DC
         e=fft[i][j][0]*fft[i][j][0]+fft[i][j][1]*fft[i][j][1];
         if ((i>=vertMin) && (i<=vertMax) && ((j<=hor) || (j>=size-hor))) s1+=e;
         else s2+=e;
       }
      }
      return s1/(s1+s2); // fraction inside selected area, use as likelyhood of the needed target
 }

float [] circle2DoubleRect (float [] pixels, int width, int height, int size, double x0, double y0, double r0) {
  float [] outPixels=new float[size*size];
  int x,y;
  double a,r;
  int px,py;
  for (y=0;y<size;y++) for (x=0;x<size;x++) {
    if ((y>(size>>1)) ||  ((y==(size>>1)) && (x>=(size>>1)))){
      r=y-(size>>1); // +0.5?
      a=size-x-1;
    } else {
      r=(size>>1)-y; // -0.5?
      a=size+x;
    }
    r+=r0; /// to match periodic pattern on both sides of zero (center circle is twice wider)
    a*=Math.PI/size;
    px=((int)Math.round(x0+r*Math.cos(a)) + width ) % width;
    py=((int)Math.round(y0+r*Math.sin(a)) + height) % height;
    outPixels[y*size+x]=pixels[py*width+px];    
  }
  return outPixels;
}

  double [][] clusteriseTargets(float [] pixels,int width, int height, double threshold, int clusterSize) {
    if ((width*height) != pixels.length) {
      IJ.showMessage("Error in  clasteriseTargets","pixels.length ("+pixels.length+") does not match width ("+width+") x height ("+height+") = "+(width*height));
      return null;
    }
    
    int x,y,i,j;
    Integer Index, NewIndex, NextIndex;
    int clusterNumber=1;
    int []clusterMap=new int[width*height];
    List <Integer> pixelList=new ArrayList<Integer>(100);
    int [] dirs={1,-width+1,-width,-width-1,-1,+width-1,width,width+1};
    int listIndex;
    float f;
    boolean first;
    double cx,cy,cm; // for centroid calculation;
    int ix,iy;
    Double [] cxy=null;
    for (i=0;i<clusterMap.length;i++) clusterMap[i]=0; /// 0 - unused, -1 - "do not use"
    List <Double[]> Centroids=new ArrayList<Double[]>(100);
    for (y=0;y<height;y++) for (x=0;x<width;x++) {
      if ((pixels[y*width+x]>=threshold) && (clusterMap[y*width+x]==0)) {
/// mark all connected pixels above the threshold
         Index=y*width+x;
         pixelList.clear();
         pixelList.add (Index);
         clusterMap[Index]=clusterNumber;
         listIndex=0;
         while (listIndex<pixelList.size() ) {
           Index=pixelList.get(listIndex++);
           for (i=0;i<dirs.length;i++) {
             NewIndex=Index+dirs[i];
             if ((NewIndex>=0) && (NewIndex<clusterMap.length) && (clusterMap[NewIndex]==0) && (pixels[NewIndex]>=threshold)) {
               pixelList.add (NewIndex);
               clusterMap[NewIndex]=clusterNumber;
             }
           }
         } //  while (!pixelList.isEmpty() )
    if (DEBUG_LEVEL>2) {
   System.out.println("Cluster="+clusterNumber+", n="+pixelList.size()+" x="+x+" y="+y);
    }
	 if (clusterSize>0) { // 0 - leave as is
           if (pixelList.size()>clusterSize) { // shrink
             while (pixelList.size()>clusterSize) {
               i=0;
               f=pixels[pixelList.get(i)];
               for (j=1;j<pixelList.size();j++) if (pixels[pixelList.get(j)]<f){
                 i=j;
                 f=pixels[pixelList.get(j)];
               }
               clusterMap[pixelList.get(i)]=-1; // Do not use looking for the next cluster
               pixelList.remove(i);
             }
           } else if (pixelList.size()<clusterSize) { // expand
             while (pixelList.size()<clusterSize) {
               first=true;
               f=0.0f;
               NextIndex=0;
               for (j=0;j<pixelList.size();j++) {
                 Index= pixelList.get(j);
                 for (i=0;i<dirs.length;i++){
                   NewIndex=Index+dirs[i];
                   if ((NewIndex>=0) && (NewIndex<clusterMap.length) && (clusterMap[NewIndex]==0) && (first || (pixels[NewIndex]>f))) {
                     f=pixels[NewIndex];
                     NextIndex=NewIndex;
                     first=false;
                   }
                 }
               }
               pixelList.add (NextIndex);
               clusterMap[NextIndex]=clusterNumber;
             }
           }
/* calculate centroid */
           cx=0.0; cy=0.0; cm=0.0;
           for (i=0;i<pixelList.size();i++) {
             j=pixelList.get(i);
             iy=j/width;
             ix=j-width*iy;
//   System.out.println("j="+j+" x="+ix+" y="+iy);

             f=pixels[j];
             cm+=f;
             cx+=f*ix;
             cy+=f*iy;
           }
           cx/=cm;
           cy/=cm;
           cxy=new Double[2];
           cxy[0]=cx;
           cxy[1]=cy;
           Centroids.add(cxy);
//   System.out.println("New cluster size="+pixelList.size()+" x="+cx+" y="+cy);

         }
         clusterNumber++;
      }
    }
    double[][] coordList=new double[Centroids.size()][2];
    for (i=0;i<coordList.length;i++) {
      coordList[i][0]=Centroids.get(i)[0];
      coordList[i][1]=Centroids.get(i)[1];
    }
//    Double[][] coordList= (Double[][]) Centroids.toArray();
  if (DEBUG_LEVEL>2) {
    for (i=0;i<coordList.length;i++) System.out.println(i+": x="+coordList[i][0]+" y="+coordList[i][1]);
  }
//    return clusterMap;
    return coordList;
  }
////    System.out.println("measureTargets(), N="+N);

/* Normalize pixels values as ratios of difference to average in the surrounding ring to variation in the ring*/
/* TODO: don't roll over, limit */
/// BUG: Seems something wrong - if convolution kernel had DC component - generated all "1.0"

  float [] normalizeAtRing(float [] pixels, int width, int height, boolean[][] mask ) {
    if ((width*height) != pixels.length) {
      IJ.showMessage("Error in  normalizeAtRing","pixels.length ("+pixels.length+") does not match width ("+width+") x height ("+height+") = "+(width*height));
      return null;
    }
    int i,j,x,y,x1,y1,x2,y2, pre_x,pre_y;
    int nFiltPix=0;
    float [] result=new float [width*height];
    double s,s2,d, mean,sigma, meang,sigmag;
    double min=0.0;
    double max=0.0;
    boolean first;
    int ir=(mask.length-1)>>1;
    for (i=0;i<mask.length; i++) for (j=0;j<mask[0].length;j++) if (mask[i][j]) nFiltPix++;

    s= 0.0;
    s2=0.0;
    for (y=0;y<height;y++) for (x=0;x<width;x++) {
      d=pixels[y*width+x];
      s+=d;
      s2+=d*d;
    }
    meang= s/(width+height);
    sigmag=Math.sqrt(s2/(width+height)-meang*meang);

    for (y=0;y<height;y++) for (x=0;x<width;x++) {
      s= 0.0;
      s2=0.0;
      pre_y=y+ir+height; // preparing for "%", making sure it will be positive 
      pre_x=x+ir+width; // preparing for "%", making sure it will be positive 
      first=true;
      for (y1=0;y1 < mask.length; y1++) {
        y2=(pre_y-y1)%height;
        for (x1=0;x1<mask[0].length;x1++) if (mask[y1][x1]) {
          x2=(pre_x-x1)%width;
          d=pixels[y2*width+x2];
          if (first) {
            min=d;
            max=d;
            first=false;
          }
          if (d>max) max=d;
          if (d<min) min=d;
          s+=d;
          s2+=d*d;
        }
      }
      mean= s/nFiltPix;
      sigma=Math.sqrt(s2/nFiltPix-mean*mean);
      sigma=Math.sqrt(sigma*sigmag); // average with image-global sigma
      
//      mean=max;
//      sigma=max-min;
      d=pixels[y*width+x];
      if (sigma>0) { // should always be so
        result[y*width+x]= (float) ((d-mean)/sigma);
      } else {
        result[y*width+x]=1.0f; // any number?
      }
    }
    return result;
  }

/* Convolve image (one Bayer slice) with the inverted target kernel
    Center should be at size/2, size/2 - will convolve only (size-1)*(size-1) */
/**Which is faster - double or float? Double i TWICE faster!*/
  float [] doubleConvolveWithTarget(float [] pixels, float [] kernel_full, int width, int height, int size) {
    int hsize=size/2;
    double [] kernel;
    
    if ((width*height) != pixels.length) {
      IJ.showMessage("Error","pixels.length ("+pixels.length+") does not match width ("+width+") x height ("+height+") = "+(width*height));
      return null;
    }
    if ((size*size) != kernel_full.length) {
      IJ.showMessage("Error","kernel.length ("+kernel_full.length+") does not match size ("+size+") ^2  = "+(size*size));
      return null;
    }
    int i,j;
    double d; // is float faster than double? or opposite (then it makes sesne to convert everything to double first
/* if kernel has even dimensions - ignore first (0) row and first (0) column */

    if ((size & 1)!=0) {
      kernel= new double[size*size];
      for (i=0;i<kernel.length;i++) kernel[i]=kernel_full[i];
    } else {
      size-=1;
      hsize-=1;
      d=0.0f;
      kernel= new double[size*size];
      for (i=0;i<size;i++) for (j=0;j<size;j++) {
        kernel[i*size+j]=kernel_full[(i+1)*(size+1)+(j+1)];
        d+=kernel[i*size+j];
      }
      d/=size*size;
//    System.out.println("Subtracting average value ("+d+") from the convolution kernel");
      for (i=0;i<kernel.length;i++) kernel[i]-=d;
    }
    double [] dPixels=new double[pixels.length];
    for (i=0;i<pixels.length;i++) dPixels[i]=pixels[i];
    
    if (DEBUG_LEVEL>10) IJ.showMessage("Debug doubleConvolveWithTarget","pixels.length="+pixels.length+"\nwidth="+width+"\nheight="+height+"\nkernel.length="+kernel.length+"\nsize="+size);
    float [] result=new float [width*height]; /* this is still float - one conversion on tghe output*/
    int x,y,x1,y1, x2, y2, pre_y,pre_x;
//    double d;
    boolean yMiddle=false;
    int index_kernel, index_source;
    for (y=0;y<height;y++) {


/**/
      yMiddle= (y>=hsize) && (y<(height-hsize));
      if (yMiddle) { // calculate faster when no need to care about borders
        for (x=hsize;x<width-hsize;x++) {
          d=0;
          index_kernel=0;
          index_source=(y+hsize)*width+x+hsize;
          for (y1=0;y1<size;y1++) {
            for (x1=0;x1<size;x1++) {
//   if (index_source<0)   System.out.println("index_source="+index_source+" index_kernel="+index_kernel+" x="+x+" y="+y+" x1="+x1+" y1="+y1);
              d+=dPixels[index_source--]*kernel[index_kernel++]; ///out of bounds: -834
            }
            index_source-=width-size; 
          }
          result[y*width+x]= (float) d;
        }
      }
/**/
/* now finish beginnings/ends of the middle lines and process complete first/last lines*/
      pre_y=y+hsize+height; // preparing for "%", making sure it will be positive 
      for (x=0;x<width;x++) if ((x<hsize) || (x>=(width-hsize)) || !yMiddle){
        d=0;
        pre_x=x+hsize+width; // preparing for "%", making sure it will be positive 
        for (y1=0;y1<size;y1++) {
          y2=(pre_y-y1)%height;
          for (x1=0;x1<size;x1++) {
            x2=(pre_x-x1)%width;
            d+=dPixels[y2*width+x2]*kernel[y1*size+x1];
          }
          result[y*width+x]= (float) d;
        }
      }
    }
    return result;
  } 
  





  public void processWindowEvent(WindowEvent e) {
    super.processWindowEvent(e);
    if (e.getID()==WindowEvent.WINDOW_CLOSING) {
      instance = null;	
    }
  }

  public boolean showDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Target Points parameters");
    gd.addStringField ("Filename prefix:                   ", jp4_instance.getTitle(), 20);
    gd.addNumericField("FFT_Size:                          ", FFTSize, 0);
//    gd.addNumericField("Target minimal outer diameter (pix)", targetOuterDMin, 2);
//    gd.addNumericField("Target maximal outer diameter (pix)", targetOuterDMax, 2);
//    gd.addNumericField("Number of target rings             ", numTargetRings, 0);
//    gd.addNumericField("Subdivide pixels for target generation ", pixelsSubdivide, 0);

    gd.addNumericField("Filtered radius (pix)             ", filteredRadius,   2); //3;  //pix - search maximums after convolution in (2*filteredRadius+1) squares
    gd.addNumericField("Background radius (pix)           ", backgroundRadius, 3); //25; //pix - consider ring area between backgroundRadius and filteredRadius as reference
    gd.addNumericField("Cluster threshold                 ", clusterThreshold, 3); //1.5
    gd.addNumericField("Cluster size (pix)                ", clusterSize, 0); //20

    gd.addNumericField("Target discriminator angular freq.   ", discrAngularFreq,  0); //2 ;  // pixels on FFT image of tragets converted polar (the smaller, the less angular variations)
    gd.addNumericField("Target discriminator radial min freq ", discrRadialMinFreq, 0); //10 ; // pixels on FFT image of tragets converted polar (radial component)
    gd.addNumericField("Target discriminator radial max freq ", discrRadialMaxFreq, 0); //10 ; // pixels on FFT image of tragets converted polar (radial component)
    gd.addNumericField("Target discriminator threshold    ", discrThreshold, 3); //0.3; // FFT energy fraction in selecter area should be > this threshold to pass the test
    gd.addNumericField("Max chromatic aberration (pix)    ", maxChromaticDistance, 1); //10.0; // Maximal distance between the same target on different color copmponents


    gd.addNumericField("Debug Level:                       ", MASTER_DEBUG_LEVEL, 0);


    gd.showDialog();
    if (gd.wasCanceled()) return false;
    jp4_instance.setTitle(gd.getNextString());
    FFTSize=1;
    for (i=             (int) gd.getNextNumber(); i >1; i>>=1) FFTSize <<=1; /* make FFTSize to be power of 2*/
//    targetOuterDMin =         gd.getNextNumber(); // minimal outer diameter of the target image , in pixels
//    targetOuterDMax =         gd.getNextNumber(); // maximal outer diameter of the target image , in pixels
//    numTargetRings  =   (int) gd.getNextNumber();  // number of target black rings (notg counting center black circle)
//    pixelsSubdivide  =   (int) gd.getNextNumber();  // Subdivide pixels for target generation

    filteredRadius   =      gd.getNextNumber();
    backgroundRadius =      gd.getNextNumber();
    clusterThreshold=       gd.getNextNumber();
    clusterSize=      (int) gd.getNextNumber();

    discrAngularFreq=    (int) gd.getNextNumber();
    discrRadialMinFreq=  (int) gd.getNextNumber();
    discrRadialMaxFreq=  (int) gd.getNextNumber();
    discrThreshold=            gd.getNextNumber();
    maxChromaticDistance=      gd.getNextNumber();
    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
    return true;
  }

  public float []createTargetDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Target template parameters");
    gd.addNumericField("FFT_Size:                          ", FFTSize, 0);
    gd.addNumericField("Target minimal outer diameter (pix)", targetOuterDMin, 2);
    gd.addNumericField("Target maximal outer diameter (pix)", targetOuterDMax, 2);
    gd.addNumericField("Number of target rings             ", numTargetRings, 0);
    gd.addNumericField("Subdivide pixels for target generation ", pixelsSubdivide, 0);
    gd.addNumericField("Invert deconvolution if less than",   deconvInvert, 3);

//    gd.addNumericField("Debug Level:                       ", MASTER_DEBUG_LEVEL, 0);


    gd.showDialog();
    if (gd.wasCanceled()) return null;
    FFTSize=1;
    for (i=             (int) gd.getNextNumber(); i >1; i>>=1) FFTSize <<=1; /* make FFTSize to be power of 2*/
    targetOuterDMin =         gd.getNextNumber(); // minimal outer diameter of the target image , in pixels
    targetOuterDMax =         gd.getNextNumber(); // maximal outer diameter of the target image , in pixels
    numTargetRings  =   (int) gd.getNextNumber();  // number of target black rings (notg counting center black circle)
    pixelsSubdivide  =  (int) gd.getNextNumber();  // Subdivide pixels for target generation
    deconvInvert=             gd.getNextNumber();  //0.05; // when FFT component is lass than this fraction of the maximal value, replace 1/z with Z
//    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
    return createTarget(FFTSize,pixelsSubdivide,targetOuterDMin,targetOuterDMax,numTargetRings,deconvInvert);
  }

  public float [] createTarget(int size, int subdiv, double DMin, double DMax, int nRings, double deconvInvert) {
   ImageProcessor ip_target;
   FHT fht_target;
   double[][][] fft_target;
   int hsizeP1= (size>>1)+1;
   double [][] dpixels=new double [hsizeP1][hsizeP1];
   double [] rMinIn=   new double[nRings+1];
   double [] rMaxIn=   new double[nRings+1];
   double [] rMinOut=  new double[nRings+1];
   double [] rMaxOut=  new double[nRings+1];
   int i,j,i1,j1,n;
   double x,y,r,ks,ke;
   double subFraction=1.0/(subdiv*subdiv);
   double DCLevel=0.0;
   double a,k,r2,k2;
   if (DMin>DMax) {
     x=DMin;
     DMin=DMax;
     DMax=x;
   }
   for (n=0;n<=nRings;n++) {
     rMinIn[n]= DMin*(n*2  )/(2*(2*nRings+1));
     rMinOut[n]=DMin*(n*2+1)/(2*(2*nRings+1));
     rMaxIn[n]= DMax*(n*2  )/(2*(2*nRings+1));
     rMaxOut[n]=DMax*(n*2+1)/(2*(2*nRings+1));
   }

   for (i=0;i<hsizeP1; i++) for (j=0;j<hsizeP1; j++) {
     dpixels[i][j]=0.0;
     for (i1=0;i1<subdiv; i1++) for (j1=0;j1<subdiv; j1++) {
       x=j+0.1*j1;
       y=i+0.1*i1;
       r=Math.sqrt(x*x+y*y);
       for (n=0;n<=nRings;n++) if ((rMinIn[n] <= r)&& (rMaxOut[n]>r)){
         if (rMaxIn[n]>r)   ke=(r-rMinIn[n])/(rMaxIn[n]-rMinIn[n]);
         else ke=1.0;
         if (rMinOut[n]<=r) ks=(r-rMinOut[n])/(rMaxOut[n]-rMinOut[n]);
         else ks=0.0;
///         dpixels[i][j]+=subFraction*(ke-ks);
         dpixels[i][j]-=subFraction*(ke-ks);
       }
     }
     r=dpixels[i][j];
// some piuxels will appear once, some - twice, most - four times
     if ((i>0) &&(i<(size>>1))) r*=2.0; 
     if ((j>0) &&(j<(size>>1))) r*=2.0;
     DCLevel+=r;
   }
   DCLevel/=(size*size);
   for (i=0;i<hsizeP1; i++) for (j=0;j<hsizeP1; j++)  dpixels[i][j]-=DCLevel;

   ip_target = new FloatProcessor(FFTSize,FFTSize);
   for (i=0;i<size; i++) for (j=0;j<size; j++) {
//     ip_target.putPixelValue(j,i, (float) dpixels[(i>=hsizeP1)?(size-i):i][(j>=hsizeP1)?(size-j):j]);
     ip_target.putPixelValue(j ^ (size>>1),i ^ (size>>1), (float) dpixels[(i>=hsizeP1)?(size-i):i][(j>=hsizeP1)?(size-j):j]);
   }
   ip_target.resetMinAndMax();
   if (DEBUG_LEVEL>5) {
     ImagePlus imp_target=  new ImagePlus("Target_direct_"+deconvInvert, ip_target);
     imp_target.show();
   }

   fht_target =  new FHT(ip_target);
// Swapping quadrants, so the center will be 0,0
   fht_target.swapQuadrants();
// get to frequency
   fht_target.transform();
   float [] fht_target_pixels=(float []) fht_target.getPixels();
   
   if (DEBUG_LEVEL>5) {
     ImageProcessor ip_fht_target = new FloatProcessor(size,size);
     ip_fht_target.setPixels(fht_target_pixels);
     ip_fht_target.resetMinAndMax();
     ImagePlus imp_fht_target= new ImagePlus("FHT_"+deconvInvert, ip_fht_target);
     imp_fht_target.show();
   }


// Convert from FHT to complex FFT
   fft_target= FHT2FFTHalf (fht_target);
/* */

/// deconvInvert
/// Now tricky thing. Invert Z for large values, but make them Z - for small ones. So it will be a mixture of correlation and deconvolution
// here the targets are round, but what will the the correct way fo assymmetrical ones?


/// First - find maximal value
   
//   double[][][] fft_target;
   double fft_max=0;
   for (i=0;i<fft_target.length; i++) for (j=0;j<fft_target[0].length;j++) {
     r2=fft_target[i][j][0]*fft_target[i][j][0]+fft_target[i][j][1]*fft_target[i][j][1];
     if (r2>fft_max) fft_max=r2;
   }
   k=Math.sqrt(fft_max)*deconvInvert;
   k2=k*k;

   for (i=0;i<fft_target.length; i++) for (j=0;j<fft_target[0].length;j++) {
     r=Math.sqrt(fft_target[i][j][0]*fft_target[i][j][0]+fft_target[i][j][1]*fft_target[i][j][1]);
     a=-Math.atan2(fft_target[i][j][1],fft_target[i][j][0]); /// will be zero for these targets)
     r=r/(r*r+k2);
     fft_target[i][j][0]=r*Math.cos(a);
     fft_target[i][j][1]=r*Math.sin(a);
   }   

// Convert fft array back to fht array
/**/

   fht_target_pixels= FFTHalf2FHT (fft_target);
// set fht_target pixels with new values 
   fht_target.setPixels (fht_target_pixels);
/// optionally show the result
   if (DEBUG_LEVEL>5) {
     ImageProcessor ip_fht_target1 = new FloatProcessor(size,size);
     ip_fht_target1.setPixels(fht_target_pixels);
     ip_fht_target1.resetMinAndMax();
     ImagePlus imp_fht_target1= new ImagePlus("Inverted_FHT_"+deconvInvert, ip_fht_target1);
     imp_fht_target1.show();
   }
/// transform 

   fht_target.inverseTransform();

   fht_target.swapQuadrants();


   fht_target.resetMinAndMax();
//   ImagePlus imp= new ImagePlus(title, ip_fht);
   if (DEBUG_LEVEL>1) {
     ImagePlus imp_target_inverted= new ImagePlus("Inverted_"+deconvInvert, fht_target);
     imp_target_inverted.show();
   }
//   return direct_target;
   return (float[])fht_target.getPixels();
  }




 /* ignore ROI, use whole image */
  public float[][] splitBayer (ImagePlus imp) {
    ImageProcessor ip=imp.getProcessor();
    Rectangle r=new Rectangle(imp.getWidth(),imp.getHeight());
    float [] pixels;
    pixels=(float[])ip.getPixels();    
    if (DEBUG_LEVEL>10) IJ.showMessage("splitBayer","r.width="+r.width+
                              "\nr.height="+r.height+
                              "\nlength="+pixels.length);
    float [][] bayer_pixels=new float[4][pixels.length>>2];
    int x,y,base,base_b,bv;
    int half_height=r.height>>1;
    int half_width=r.width>>1;
    for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
      base=r.width*((y<<1)+bv);
      base_b=half_width*y;
      if (bv==0) for (x=0; x<half_width; x++) {
        bayer_pixels[0][base_b]=  pixels[base++];
        bayer_pixels[1][base_b++]=pixels[base++];
      } else  for (x=0; x<half_width; x++) {
        bayer_pixels[2][base_b]=  pixels[base++];
        bayer_pixels[3][base_b++]=pixels[base++];
      }
    }
    return bayer_pixels;
  }

  public void showBayers(float[][] bayer_pixels, int width, int height, String title) {
    int i;
   if (DEBUG_LEVEL>10) IJ.showMessage("showBayers","width="+width+
                              "\nheight="+height+
                              "\nlength="+bayer_pixels[0].length);

    ImageProcessor[] ip= new ImageProcessor[4];
    ImagePlus[]      imp=new ImagePlus[4];
    for (i=0;i<4;i++) {
      ip[i]=new FloatProcessor(width,height);
      ip[i].setPixels(bayer_pixels[i]);
      ip[i].resetMinAndMax();
      imp[i]=  new ImagePlus(title+"_"+i, ip[i]);
      imp[i].show();
    }
  }



  public float[] initHamming(int size) {
    float [] hamming =new float [size*size];
    float [] hamming_line=new float [size];
    int i,j;
    for (i=0; i<size; i++) hamming_line[i]=(float) (0.54-0.46*Math.cos((i*2.0*Math.PI)/size));
    for (i=0; i<size; i++) for (j=0; j<size; j++){
       hamming[size*i+j]=hamming_line[i]*hamming_line[j];
    }
    return hamming;
  }


  public void initDisplay() {
    if ((imp_display==null) || (imp_display.getWidth()!=displayWidth) || (imp_display.getHeight()!=displayHeight)) {
      if (imp_display!=null)     imp_display.close();
      ip_display=   new ColorProcessor (displayWidth,displayHeight);
      imp_display=  new ImagePlus("Target Points", ip_display);
      imp_display.show();

    }
  }


/* converts FHT results (frequency space) to complex numbers of [FFTSize/2+1][FFTSize] */
private double[][][] FHT2FFTHalf (FHT fht) {
   float[] fht_pixels=(float[])fht.getPixels();
   double[][][] fftHalf=new double[(FFTSize>>1)+1][FFTSize][2];
   int row1,row2,col1,col2;

   for (row1=0;row1<=(FFTSize>>1);row1++) {
     row2=(FFTSize-row1) %FFTSize;
     for (col1=0;col1<FFTSize;col1++) {
       col2=(FFTSize-col1) %FFTSize;
       fftHalf[row1][col1][0]=   0.5*(fht_pixels[row1*FFTSize+col1] + fht_pixels[row2*FFTSize+col2]);
       fftHalf[row1][col1][1]=   0.5*(fht_pixels[row2*FFTSize+col2] - fht_pixels[row1*FFTSize+col1]);
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

