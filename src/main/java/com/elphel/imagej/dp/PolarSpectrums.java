package com.elphel.imagej.dp;
/**
** -----------------------------------------------------------------------------**
** PolarSpectrums.java
**
** Used in "scissors" frequency-domain demosaic
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
import java.util.HashSet;
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
//  private int    [][][] polarAliases;          // a,r,number triplet of the "master" polar cell and number of aliases (including this))- the one that returns by the decart-> polar  conversion of the polar->decart
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
//   System.out.println("polar2CartesianFractions["+i+"].length="+polar2CartesianFractions[i].length);
//    System.out.println("polar2CartesianIndices["+i+"].length="+    polar2CartesianIndices[i].length);

      polPixels[i]=(1-polar2CartesianFractions[i][1])*( (1-polar2CartesianFractions[i][0])*cartPixels[polar2CartesianIndices[i][0]]  + polar2CartesianFractions[i][0]*cartPixels[polar2CartesianIndices[i][1]])+
                      polar2CartesianFractions[i][1] *( (1-polar2CartesianFractions[i][0])*cartPixels[polar2CartesianIndices[i][2]]  + polar2CartesianFractions[i][0]*cartPixels[polar2CartesianIndices[i][3]]) ;
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

  public double [] genPolarRedBlueMask(double [] polarAmps, // polar array of amplitude values, <0 - stop
                                                 int mode){ // result mode - 0: output mask as 0/1, 1 -output proportional, positive - pass, negative - rejected
       return genPolarMask(polarAmps,1,mode);
  }


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
//      index=iMax*iRadiusPlus1+rayLength[iMax]; // rayLength[iMax] should point to a new cell (with intMap[]==0) may ommit - set in the end of the loop and before the loop?
      intMap[index]=newVal;
      if (sameCartesian[index]!=null) for (i=0;i<sameCartesian[index].length;i++) intMap[sameCartesian[index][i]]=newVal;
/* add aliases of point index (as negative values) */
      if ((good) &&(polarMap[index]!=null)) for (i=0;i<polarMap[index].length;i++) intMap[polarMap[index][i]]=-step;
/* update ray lengths and status */
      max=-1.0;
      iMax=-1;
      for  (ia=0;ia<=iAngle;ia++) if (rayOpen[ia]) {
        base=ia*iRadiusPlus1+1;  // second for this angle
        l=base+rayLength[ia]; // first after the pointer
        lastIndex=base+iRadiusPlus1; // first in the next row
        while ((l<lastIndex) && (intMap[l]>0)) l++;
        rayLength[ia]=l-base; // last "good" ( >0 and in the same row)
        if ((l==lastIndex) || (intMap[l]<0) || (polarAmps[l]<0.0) ) rayOpen[ia]=false;
        else {
          if (polarAmps[l]>max) {
            max=polarAmps[l];
            iMax=ia;
          }
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
      int isize, // size of the square array, centar is at size/2, size/2, only top half+line will be used
//      double  mask, //1d - array of pixels, maybe size*(size/2+1) or full
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
//    HashSet  [] polarList=  new HashSet [length];
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
//if (PolarIndex<5)  System.out.println(">>>>> x="+x+" y="+y+" dx="+dx+" dy="+dy+" i="+i+" j="+j+" CartesianIndex="+CartesianIndex+" polar2CartesianIndices["+PolarIndex+"][0]="+polar2CartesianIndices[PolarIndex][0]);
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
//   System.out.println("testMapsLengths("+mode+").length="+result.length);
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
/*    int [][] aliasMapGreen=  {{-2,-2},{-2,0},{-2,2},
                                  {-1,-1},{-1,1},
                              { 0,-2},       { 0,2},
                                  { 1,-1},{ 1,1},
                              { 2,-2},{ 2,0},{ 2,2}};*/
    int [][] aliasMapGreen=  {{-2,-2},{-2,0},            // using rollover, so only unique aliases are needed
                                  {-1,-1},{-1,1},
                              { 0,-2},
                                  { 1,-1},{ 1,1}};
/*    int [][] aliasMapRedBlue={{-1,-1},{-1,0},{-1,1},
                              { 0,-1},       { 0,1},
                              { 1,-1},{ 1,0},{ 1,1}};*/

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
//if (polarIndex<5)  System.out.println("ix="+ix+" iy="+iy+" polar2CartesianIndices["+polarIndex+"][0]="+polar2CartesianIndices[polarIndex][0]);
//if (polarIndex<5)  System.out.println("ix="+ix+" iy="+iy+" polar2CartesianIndices["+polarIndex+"][1]="+polar2CartesianIndices[polarIndex][1]);
//if (polarIndex<5)  System.out.println("ix="+ix+" iy="+iy+" polar2CartesianIndices["+polarIndex+"][2]="+polar2CartesianIndices[polarIndex][2]);
//if (polarIndex<5)  System.out.println("ix="+ix+" iy="+iy+" polar2CartesianIndices["+polarIndex+"][3]="+polar2CartesianIndices[polarIndex][3]);
      aliasList.clear();
      for (nAlias=0;nAlias<aliasMap.length;nAlias++) for (dirAlias=-1;dirAlias<2;dirAlias+=2) {
        ixa=(size+ size/2+ aliasMap[nAlias][0]*size/4+ dirAlias*ix) % size;
        iya=(size+ size/2- aliasMap[nAlias][1]*size/4- dirAlias*iy) % size;
        index=iya*size + ixa;
        if (cartesian2PolarIndices[index]==null) {
//if (polarIndex<5)  System.out.println("cartesian2PolarIndices["+index+"]=null");
          if (cartesian2PolarIndex[index]>=0) {
            aliasList.add (new Integer(cartesian2PolarIndex[index]));
//if (polarIndex<5)  System.out.println("cartesian2PolarIndex["+index+"]="+cartesian2PolarIndex[index]+ " (ir="+(cartesian2PolarIndex[index] % iRadiusPlus1)+" ia="+(cartesian2PolarIndex[index] % iRadiusPlus1)+")");
          }
        } else {
          for (j=0;j<cartesian2PolarIndices[index].length;j++) {
             aliasList.add (new Integer(cartesian2PolarIndices[index][j]));
//if (polarIndex<5)  System.out.println("cartesian2PolarIndices["+index+"]["+j+ "]="+cartesian2PolarIndices[index][j]+ " (ir="+(cartesian2PolarIndices[index][j] % iRadiusPlus1)+" ia="+(cartesian2PolarIndices[index][j] % iRadiusPlus1)+")");
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
//     System.out.println("polarIndex="+polarIndex+" cartesianIndex="+cartesianIndex+ " polar2CartesianIndices["+polarIndex+"][0]"+polar2CartesianIndices[polarIndex][0]);
//     System.out.println("polar2CartesianIndices["+polarIndex+"]="+polar2CartesianIndices[polarIndex].length);

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
