package com.elphel.imagej.calibration;
import java.awt.Rectangle;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;

import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

import ij.IJ;

/*
 **
 ** SimulationPattern.java - Generate simulated pattern
 **
 ** Copyright (C) 2010-2011 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  SimulationPattern.java is free software: you can redistribute it and/or modify
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

public class SimulationPattern {
	//	private double []  bPattern; // pattern bitmap (does not change)
	//	private double[][] barray;   // high resolution boolean pattern array (specific to distortions in each area)
	public double []  bPattern=null; // pattern bitmap (does not change)
	public int        bPatternSize=0;
	///	public double[][] barray;   // high resolution boolean pattern array (specific to distortions in each area)
	public double[] barray;   // high resolution boolean pattern array (specific to distortions in each area)
	public double bPatternSigma=0.0;
	public double barraySigma=0.0;
	public  int         debugLevel=2;
	private DoubleGaussianBlur gb = new DoubleGaussianBlur();
	private ShowDoubleFloatArrays SDFA_INSTANCE= new ShowDoubleFloatArrays(); // just for debugging?

	public SimulationPattern (){
		this.bPattern=null;
	}
	public SimulationPattern (double [] bPattern){ // reuse the same barray
		this.bPattern=bPattern;
	}
	public SimulationPattern(
			SimulParameters simulParameters ) {
		this.bPatternSigma=simulParameters.bPatternSigma;
		this.barraySigma=simulParameters.barraySigma;
		patternGenerator(simulParameters);
	}
	public SimulationPattern(
			int size,
			int patternNumber,
			double patternModifier) {
		patternGenerator(size,patternNumber,patternModifier);
	}
	/* ======================================================================== */
	public double [] patternGenerator(
			SimulParameters simulParameters ) {
		return patternGenerator(
				simulParameters.patternSize,
				simulParameters.pattern_type,
				simulParameters.pattern_modifier);
	}
	public double [] patternGenerator(int size,
			int patternNumber,
			double patternModifier) {
		this.bPattern=new double [size*size];
		this.bPatternSize=size;
		int i,j,index,k;
		double p;
		double a,r,r2,h;
		double qSize=size/4;
		switch (patternNumber) {
		case 1:
			a=patternModifier*(Math.sqrt(2)-1.0);
			r=(a*a+1)/(2*a)*qSize;
			r2=r*r;
			h=Math.sqrt(r2-qSize*qSize);
			if (a>1.0) h=-h;
			double [][] pattern1Centers={{qSize,         -h},
					{ size+h,       qSize},
					{ size-qSize,   size+h},
					{-h,            size-qSize}};
			index=0;
			for (i=0;i<size;i++) for (j=0;j<size;j++) {
				p=1.0;
				for (k=0;k<pattern1Centers.length;k++) if ((((i-pattern1Centers[k][1])*(i-pattern1Centers[k][1])+(j-pattern1Centers[k][0])*(j-pattern1Centers[k][0])))<r2) p=0.0;
				this.bPattern[index++]=p;
			}
			break;
		case 2:
			index=0;
			for (i=0;i<size;i++) for (j=0;j<size;j++) {
				p= ((i>=0.3*size) && (i<0.7*size) && (j>=0.3*size) && (j<0.7*size))?1.0:0.0;
				this.bPattern[index++]=p;
			}
			break;
		case 3:
			index=0;
			for (i=0;i<size;i++) for (j=0;j<size;j++) {
				p= ((i>=0.1*size) && (i<0.9*size) && (j>=0.1*size) && (j<0.9*size))?1.0:0.0;
				this.bPattern[index++]=p;
			}
			break;
		default: for (index=0;index<this.bPattern.length;index++) this.bPattern[index]=1.0;
		}
		// blur pattern
		if (this.bPatternSigma>0) {
			if (this.bPatternSigma>0.25) this.bPatternSigma=0.25;
			// 1 - add margins around the pattern
			int i1,j1;
			int margin= (int) Math.ceil(3*size*this.bPatternSigma);
			int sizeM=size+2*margin;
			boolean invertY,invertX;
			double [] bPatternM=new double [sizeM*sizeM];
			for (i=0;i<sizeM;i++) {
				i1= (i+size-margin)%size;
				invertY=(((i+size-margin)/size)&1)==0;
				for (j=0;j<sizeM;j++) {
					invertX=(((j+size-margin)/size)&1)==0;
					j1= (j+size-margin)%size;
					bPatternM[i*sizeM+j]= (invertX ^ invertY)?(1.0-this.bPattern[i1*size+j1]):this.bPattern[i1*size+j1];
				}
			}
			// apply blur
			if (this.debugLevel>3) SDFA_INSTANCE.showArrays(bPatternM,sizeM,sizeM, "bPatternM");
			this.gb.blurDouble(bPatternM,sizeM,sizeM,size*this.bPatternSigma,size*this.bPatternSigma, 0.01);
			if (this.debugLevel>3) SDFA_INSTANCE.showArrays(bPatternM,sizeM,sizeM, "bPatternM-blured");
			// remove margins
			for (i=0;i<size;i++) for (j=0;j<size;j++) {
				this.bPattern[i*size+j]= bPatternM[(i+margin)*sizeM+(j+margin)];
			}
		}
		return this.bPattern;
	}

	/* ======================================================================== */
	public void simulatePatternFullPattern(
			double freqX1,
			double freqY1,
			double phase1,
			double freqX2,
			double freqY2,
			double phase2,
			double [] corr,
			int subdiv,
			int size,
			boolean center_for_g2,
			boolean mono) {// takes precedence over center_for_g2, makes symmetrical around center of array (between pixels)
		this.barray=simulatePatternFullPatternSafe(
				freqX1,
				freqY1,
				phase1,
				freqX2,
				freqY2,
				phase2,
				corr,
				subdiv,
				size,
				center_for_g2,
				mono);

	}


	public double [] simulatePatternFullPatternSafe(
			double freqX1,
			double freqY1,
			double phase1,
			double freqX2,
			double freqY2,
			double phase2,
			double [] corr,
			int subdiv,
			int size,
			boolean center_for_g2,
			boolean mono) {// takes precedence over center_for_g2, makes symmetrical around center of array (between pixels)
		int patternSize= (this.bPattern!=null)?((int) Math.sqrt(this.bPattern.length)):0;
		double twicePatternSize=2*patternSize;
		int i,j;
		int fullSize=subdiv*(size+4)*2;
		double [] localBarray=new double [fullSize*fullSize];
		double xl,yl; //,x,y;//,p1,p2;


		double [][] xy2uv= {{freqX1,freqY1},
				{freqX2,freqY2}};

		if (this.debugLevel>2) {
			System.out.println("simulatePatternFullPattern:");
			System.out.println(" Ax="+IJ.d2s(corr[0],5)+" Bx="+IJ.d2s(corr[1],5)+" Cx="+IJ.d2s(corr[2],5)+" Dx="+IJ.d2s(corr[6],5)+" Ex="+IJ.d2s(corr[7],5));
			System.out.println(" Ay="+IJ.d2s(corr[3],5)+" By="+IJ.d2s(corr[4],5)+" Cy="+IJ.d2s(corr[5],5)+" Dy="+IJ.d2s(corr[8],5)+" Ey="+IJ.d2s(corr[9],5));
		}

		if (this.debugLevel>2) {
			System.out.println("simulatePatternFullPattern:  xy2uv[0][0]="+IJ.d2s(xy2uv[0][0],4)+" xy2uv[0][1]="+IJ.d2s(xy2uv[0][1],4));
			System.out.println("simulatePatternFullPattern:  xy2uv[1][0]="+IJ.d2s(xy2uv[1][0],4)+" xy2uv[1][1]="+IJ.d2s(xy2uv[1][1],4));
		}

		double []uv, xy;
		xy=new double [2];

		double [] phases={phase1/(Math.PI*2)+0.25,phase2/(Math.PI*2)+0.25}; // period=1.0;
		int iu,iv;
		boolean invert;

		double offset = mono? 0.0: ((center_for_g2?0.5:1.0));
		double center = 0.5*fullSize; //  +(mono?-0.5:  0.0); // see if it is needed!
		// in mono mode result barray is symmetrical around [fullSize/2,fullSize/2], but
		// after extractSimulMono result will be symmetrical around [size/2-0.5, size/2 - 0.5]
		for (i=0; i<fullSize; i++) {
			yl=(i - center) / subdiv - offset; // center in the middle of Bayer
			for (j=0; j<fullSize; j++) {
				xl=(j - center) / subdiv - offset; // center in the middle of Bayer
				/* apply second order polynomial correction to x,y
    x=xl+Ax*xl^2+Bx*yl^2+2*Cx*xl*yl;
    y=xl+Ay*xl^2+By*yl^2+2*Cy*xl*yl; */
				if (corr==null) {
					xy[0]=xl;
					xy[1]=yl;
				} else {
					xy[0]=xl + corr[0]*xl*xl + corr[1]*yl*yl + 2* corr[2]*xl*yl + corr[6]*xl + corr[7]*yl;
					xy[1]=yl + corr[3]*xl*xl + corr[4]*yl*yl + 2* corr[5]*xl*yl + corr[8]*xl + corr[9]*yl;
				}
				uv= matrix2x2_mul(xy2uv, xy);
				uv= vector_add(uv,phases);
				uv[0]-=Math.floor(uv[0]);
				uv[1]-=Math.floor(uv[1]);
				invert=false;
				if (uv[0]>=0.5){
					invert=!invert;
					uv[0]-=0.5;
				}
				if (uv[1]>=0.5){
					invert=!invert;
					uv[1]-=0.5;
				}
				if (this.bPattern==null) {
					localBarray[i*fullSize+j]=invert?0.0:1.0; //!invert;
				} else {
					iu= (int) Math.round(uv[0]*twicePatternSize);
					iv= (int) Math.round(uv[1]*twicePatternSize);
					if ((iu<0) || (iu>=patternSize)) {
						invert=!invert;
						iu=(iu+patternSize)% patternSize;
					}
					if ((iv<0) || (iv>=patternSize)) {
						invert=!invert;
						iv=(iv+patternSize)% patternSize;
					}
					localBarray[i*fullSize+j]=invert?(1.0-this.bPattern[iv*patternSize + iu]): this.bPattern[iv*patternSize + iu];
				}
			}
		}
		// Blur barray pattern if sigma >0
		if (this.barraySigma>0) {
			double sigma=this.barraySigma*subdiv; //*/ 2?
			if (this.debugLevel>3) SDFA_INSTANCE.showArrays(localBarray, "localBarray");
			this.gb.blurDouble(localBarray,fullSize,fullSize,sigma,sigma, 0.01);
			if (this.debugLevel>3) SDFA_INSTANCE.showArrays(localBarray, "localBarray-blured");
		}
		return localBarray;
	}
	/* ======================================================================== */
	public double [] recursiveFillPixels ( // invert pattern in the caller, return signed value (-1..1 - pattern is 0..1)
			SimulParameters  simulParameters,
			double [] xy,  // top-left corner
			double [] dxy, // increments to other corners
			double [][][] cornersXY, // xy pairs for the 4 corners of the square in UV (pattern) coordinates (u0v0,u1v0,u0v1,u1v1)
			double [] uv,  // UV value for the top-left corner (matching cornersXY[0][0])
			double [] duv,  // distances to the opposite corner in UV
			//			   final boolean   maskOnly, // just mark defined cells
			int debug
			){ //use this.bPattern, this.bPatternSize (side of the square)

		double [][][] cornersUV=new double [2][2][];
		double [] xy4=new double[2];
		double [] result ={0.0,0.0};
		int numInside=0;
		for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
			xy4[0]=xy[0]+j*dxy[0];
			xy4[1]=xy[1]+i*dxy[1];
			cornersUV[i][j]=bilinearXY2UV(cornersXY,xy4,debug);
			if ((cornersUV[i][j][0]>=0.0) && (cornersUV[i][j][0]<=1.0) && (cornersUV[i][j][1]>=0.0) && (cornersUV[i][j][1]<=1.0)) numInside++;
		}
		if (debug>21){
			String dbgStr="";
			//			 IJ.d2s(quarter_patterns[iq][0][0],4)
			dbgStr+="xy={"+IJ.d2s(xy[0],2)+","+IJ.d2s(xy[1],2)+"} ";
			dbgStr+=" dxy={"+IJ.d2s(dxy[0],2)+","+IJ.d2s(dxy[1],2)+"} ";
			dbgStr+=" uv={"+IJ.d2s(uv[0],2)+","+IJ.d2s(uv[1],2)+"} ";
			dbgStr+=" duv={"+IJ.d2s(duv[0],2)+","+IJ.d2s(duv[1],2)+"} ";
			dbgStr+=" cornersXY={{{"+IJ.d2s(cornersXY[0][0][0],5)+","+IJ.d2s(cornersXY[0][0][1],5)+"},";
			dbgStr+=             "{"+IJ.d2s(cornersXY[0][1][0],5)+","+IJ.d2s(cornersXY[0][1][1],5)+"}},";
			dbgStr+=            "{{"+IJ.d2s(cornersXY[1][0][0],5)+","+IJ.d2s(cornersXY[1][0][1],5)+"},";
			dbgStr+=             "{"+IJ.d2s(cornersXY[1][1][0],5)+","+IJ.d2s(cornersXY[1][1][1],5)+"}}}";
			dbgStr+=" cornersUV={{{"+IJ.d2s(cornersUV[0][0][0],3)+","+IJ.d2s(cornersUV[0][0][1],3)+"},";
			dbgStr+=             "{"+IJ.d2s(cornersUV[0][1][0],3)+","+IJ.d2s(cornersUV[0][1][1],3)+"}},";
			dbgStr+=            "{{"+IJ.d2s(cornersUV[1][0][0],3)+","+IJ.d2s(cornersUV[1][0][1],3)+"},";
			dbgStr+=             "{"+IJ.d2s(cornersUV[1][1][0],3)+","+IJ.d2s(cornersUV[1][1][1],3)+"}}}";
			dbgStr+=" numInside="+numInside;
			System.out.println(dbgStr);
		}

		if (numInside==0) return result; // all corners outside of the (sub)pattern cell
		//		if (maskOnly) {
		if (simulParameters==null) {
			result[1]=dxy[0]*dxy[1];
			result[0]=result[1];
			return result;
		}
		// recalculate to the full uv
		boolean cornersInvert;
		double  [][] cornerValue=new double [2][2];
		int  []   iPat=new int [2];
		double min=1.0,max=-1.0;
		for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
			//			cornersUV[i][j][0]=uv[0]+j*cornersUV[i][j][0]*duv[0];
			//			cornersUV[i][j][1]=uv[1]+i*cornersUV[i][j][1]*duv[1];
			cornersUV[i][j][0]=uv[0]+cornersUV[i][j][0]*duv[0];
			cornersUV[i][j][1]=uv[1]+cornersUV[i][j][1]*duv[1];
			cornersInvert=false;
			for (int k=0;k<2;k++) {
				//				iPat[k] = (int) Math.floor(cornersUV[i][j][k]*this.bPatternSize*2.0); // 0.5 ->bPatternSize
				iPat[k] = (int) Math.floor(cornersUV[i][j][k]*this.bPatternSize); // 1.0 ->bPatternSize
				if (iPat[k]<0){
					iPat[k]+=this.bPatternSize;
					cornersInvert=!cornersInvert;
				} else if (iPat[k]>=this.bPatternSize){
					iPat[k]-=this.bPatternSize;
					cornersInvert=!cornersInvert;
				}
				if ((iPat[k]<0) || (iPat[k]>=this.bPatternSize)) {
					if (debug>0) System.out.println("Too far, cornersUV["+i+"]["+j+"]["+k+"]="+cornersUV[i][j][k]);
					return result; // {0,0} here
				}
			}
			cornerValue[i][j]=(2.0*this.bPattern[this.bPatternSize*iPat[1]+iPat[0]]-1.0)*(cornersInvert?1.0:-1.0);
			if (max<cornerValue[i][j]) max=cornerValue[i][j];
			if (min>cornerValue[i][j]) min=cornerValue[i][j];
		}
		if (((max-min)>simulParameters.bitmapNonuniforityThreshold) &&
				(dxy[0]>simulParameters.smallestSubPix) &&
				(dxy[1]>simulParameters.smallestSubPix)) {

			// divide this square into 4 quadrants, return sum of the recursively called method on them
			double [][] quadrants={{0.0,0.0},{0.5,0.0},{0.0,0.5},{0.5,0.5}};
			double [] subResult;
			double [] subxy=new double [2];
			double [] subdxy={0.5*dxy[0],0.5*dxy[1]};
			if (debug>1){
				System.out.println("---> Subdividing into "+subdxy[0]+"x"+subdxy[0]+" (max="+IJ.d2s(max,3)+" min="+IJ.d2s(min,3)+
						" max-min="+IJ.d2s(max-min,3)+
						" cornerValue: [0][0]="+IJ.d2s(cornerValue[0][0],3)+
						" [0][1]="+IJ.d2s(cornerValue[0][1],3)+
						" [1][0]="+IJ.d2s(cornerValue[1][0],3)+
						" [1][1]="+IJ.d2s(cornerValue[1][1],3)+

						" cornersUV={{{"+IJ.d2s(cornersUV[0][0][0],3)+","+IJ.d2s(cornersUV[0][0][1],3)+"},"+
						"{"+IJ.d2s(cornersUV[0][1][0],3)+","+IJ.d2s(cornersUV[0][1][1],3)+"}},"+
						"{{"+IJ.d2s(cornersUV[1][0][0],3)+","+IJ.d2s(cornersUV[1][0][1],3)+"},"+
						"{"+IJ.d2s(cornersUV[1][1][0],3)+","+IJ.d2s(cornersUV[1][1][1],3)+"}}}");
			}
			for (int i=0;i<quadrants.length;i++) {
				subxy[0]=xy[0]+dxy[0]*quadrants[i][0];
				subxy[1]=xy[1]+dxy[1]*quadrants[i][1];
				subResult=	recursiveFillPixels ( // invert pattern in the caller, return signed value (-1..1 - pattern is 0..1)
						simulParameters,
						subxy,  // top-left corner
						subdxy, // increments to other corners
						cornersXY, // xy pairs for the 4 corners of the square in UV (pattern) coordinates (u0v0,u1v0,u0v1,u1v1)
						uv,  // UV value for the top-left corner (matching cornersXY[0][0])
						duv,  // distances to the opposite corner in UV
						//						maskOnly, // just mark defined cells //always false - will never get here
						debug
						);
				result[0]+=subResult[0];
				result[1]+=subResult[1];
			}
			if (debug>1){
				if (result[1]==0.0) System.out.println("<--- Combined results "+IJ.d2s(result[0],3)+" / "+IJ.d2s(result[1],3));
				else                System.out.println("<--- Combined results "+IJ.d2s(result[0],3)+" / "+IJ.d2s(result[1],3)+"="+IJ.d2s(result[0]/result[1],3));
			}

		} else { // no more subdivisions - calculate average value, taking into account partial pixels
			for (int i=0;i<2;i++) for (int j=0;j<2;j++)  result[0]+=0.25*cornerValue[i][j];
			result[1]=dxy[0]*dxy[1];
			result[0]*=result[1];
			if (numInside <4) {
				double f=((double) numInside)/4; // estimate fraction of the pixel - start with simple number of corners
				result[0]*=f;
				result[1]*=f;
			}
			if (debug>1){
				if (result[1]==0.0)System.out.println("< === Returnimg "+IJ.d2s(result[0],3)+" / "+IJ.d2s(result[1],3)+" ("+numInside+" corners inside)");
				else System.out.println("< === Returning "+IJ.d2s(result[0],5)+" / "+IJ.d2s(result[1],5)+"="+IJ.d2s(result[0]/result[1],3)+" ("+numInside+" corners inside)");
			}
		}
		return result;
	}
	/* ======================================================================== */
	/**
	 * @param cornersXY first index V, second index U, third index:0 - x, 1-y
	 * @param xy 0-x,1-y of the point, for which UV should be generated
	 * @return UV pair
	 */
	public double [] bilinearXY2UV(
			double [][][] cornersXY, // first index V, second index U, third index:0 - x, 1-y
			double []     xy,         // 0-x,1-y of the point, for which
			int debug
			) {
		return bilinearXY2UV(
				cornersXY, // first index V, second index U, third index:0 - x, 1-y
				xy,         // 0-x,1-y of the point, for which
				1E-6,
				debug);
	}
	/**
	 * @param cornersXY first index V, second index U, third index:0 - x, 1-y
	 * @param xy 0-x,1-y of the point, for which UV should be generated
	 * @param quadThreshold if abs(4*a*c)/b^2 is less than this, use linear, not quadratic equations
	 * @return UV pair
	 */
	public double [] bilinearXY2UV(
			double [][][] cornersXY, // first index V, second index U, third index:0 - x, 1-y
			double []     xy,         // 0-x,1-y of the point, for which
			double quadThreshold,      // if abs(4*a*c)/b^2 is less than this, use linear, not quadratic equations
			int debug
			) {
		/*
x,y -> u,v

(1) x= v*u*Ax + v*Bx + u*Cx + Dx
(2) y= v*u*Ay + v*By + u*Cy + Dy

Ax=x11-x10-x01+x00
Bx=x10-x00
Cx=x01-x00
Dx=x00

Ay=y11-y10-y01+y00
By=y10-y00
Cy=y01-y00
Dy=y00

u*u*(Cy*Ax-Cx*Ay)+u*((-Cx*By-Ay*Dx+Cy*Bx+Ax*Dy)+Ay*x-Ax*y)+(By*x-Bx*y)+(-By*Dx+Bx*Dy)=0

Au*u*u+Bu*u+Cu=0
Av*v*v+Bv*v+Cv=0

Au=(Cy*Ax-Cx*Ay)
Bu=((-Cx*By-Ay*Dx+Cy*Bx+Ax*Dy)+Ay*x-Ax*y)
Cu=(By*x-Bx*y)+(-By*Dx+Bx*Dy)

Av=(By*Ax-Bx*Ay)
Bv=((-Bx*Cy-Ay*Dx+By*Cx+Ax*Dy)+Ay*x-Ax*y)
Cv=(Cy*x-Cx*y)+(-Cy*Dx+Cx*Dy)


		 */
		double Ax=cornersXY[1][1][0]-cornersXY[1][0][0]-cornersXY[0][1][0]+cornersXY[0][0][0];
		double Bx=cornersXY[1][0][0]-cornersXY[0][0][0];
		double Cx=cornersXY[0][1][0]-cornersXY[0][0][0];
		double Dx=cornersXY[0][0][0];

		double Ay=cornersXY[1][1][1]-cornersXY[1][0][1]-cornersXY[0][1][1]+cornersXY[0][0][1];
		double By=cornersXY[1][0][1]-cornersXY[0][0][1];
		double Cy=cornersXY[0][1][1]-cornersXY[0][0][1];
		double Dy=cornersXY[0][0][1];

		double Au=(Cy*Ax-Cx*Ay);
		double Bu=((-Cx*By-Ay*Dx+Cy*Bx+Ax*Dy)+Ay*xy[0]-Ax*xy[1]);
		double Cu=(By*xy[0]-Bx*xy[1])+(-By*Dx+Bx*Dy);

		double Av=(By*Ax-Bx*Ay);
		double Bv=((-Bx*Cy-Ay*Dx+By*Cx+Ax*Dy)+Ay*xy[0]-Ax*xy[1]);
		double Cv=(Cy*xy[0]-Cx*xy[1])+(-Cy*Dx+Cx*Dy);
		//		double [] UV={-Cv/Bv,-Cu/Bu}; // linear solution - use for linear grid
		double [] UV={-Cu/Bu,-Cv/Bv}; // linear solution - use for linear grid
		double au=0.0,bu=0.0,av=0.0,bv=0.0;
		if (Math.abs(Au*Cu)/(Bu*Bu)>quadThreshold) { // use quadratic equation for U
			au=-Bu/(2*Au);
			bu=Math.sqrt(Bu*Bu-4*Au*Cu)/Math.abs(2*Au);
			// Use solution that is closer to linear one
			if (UV[0]>au) UV[0]=au+bu;
			else          UV[0]=au-bu;
		}
		if (Math.abs(Av*Cv)/(Bv*Bv)>quadThreshold) { // use quadratic equation for V
			av=-Bv/(2*Av);
			bv=Math.sqrt(Bv*Bv-4*Av*Cv)/Math.abs(2*Av);
			// Use solution that is closer to linear one
			if (UV[1]>av) UV[1]=av+bv;
			else          UV[1]=av-bv;
		}
		if (debug>2){
			String dbgStr="";
			//			 IJ.d2s(quarter_patterns[iq][0][0],4)
			dbgStr+=" Ax="+IJ.d2s(Ax,5)+", Bx="+IJ.d2s(Bx,5)+", Cx="+IJ.d2s(Cx,5)+", Dx="+IJ.d2s(Dx,5);
			dbgStr+=" Ay="+IJ.d2s(Ay,5)+", By="+IJ.d2s(By,5)+", Cy="+IJ.d2s(Cy,5)+", Dy="+IJ.d2s(Dy,5);
			dbgStr+=" Au="+IJ.d2s(Au,5)+", Bu="+IJ.d2s(Bu,5)+", Cu="+IJ.d2s(Cu,5);
			dbgStr+=" Av="+IJ.d2s(Av,5)+", Bv="+IJ.d2s(Bv,5)+", Cv="+IJ.d2s(Cv,5);
			dbgStr+=" LinU="+IJ.d2s(-Cu/Bu,3)+", LinV="+IJ.d2s(-Cv/Bv,5);
			dbgStr+=" au="+IJ.d2s(au,5)+", bu="+IJ.d2s(bu,5);
			dbgStr+=" av="+IJ.d2s(av,5)+", bv="+IJ.d2s(bv,5);
			System.out.println(dbgStr);
		}

		return UV;
	}
	/* ======================================================================== */
	private boolean isCellValid(
			double [][][][] grid,
			int [] uv){
		if ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length)) {
			double [][] cell = grid[uv[1]][uv[0]];
			return ((cell!=null) && (cell.length>1));
		}
		return false;
	}
	/*
	   private boolean isCellDefined(
    		   double [][][][] grid,
    		   int [] uv){
    	   return ((uv[1]>=0) && (uv[0]>=0) && (uv[1]<grid.length) && (uv[0]<grid[uv[1]].length) &&
    			   (grid[uv[1]][uv[0]]!=null) && (grid[uv[1]][uv[0]][0]!=null));
       }
	 */
	public float [] combineWithCanvas(
			double canvasFill,
			int width,
			int height,
			Rectangle woi,
			float [] selection ){
		float []canvas=new float[width*height];
		for (int i=0;i<canvas.length;i++)canvas[i]= (float) canvasFill;
		return combineWithCanvas(canvas, width, woi, selection );
	}
	public float [] combineWithCanvas(
			float [] canvas,
			int width,
			Rectangle woi,
			float [] selection ){
		// debug
		if (selection==null) System.out.println("combineWithCanvas(): selection==null");
		if (woi==null) System.out.println("combineWithCanvas(): woi==null");
		if (selection.length!=(woi.width*woi.height)) throw new IllegalArgumentException ("selection.length="+selection.length+", woi.width="+woi.width+", woi.height="+woi.height);
		int i0=0;
		int i1=width*woi.y+woi.x;
		for (int y=0;y<woi.height;y++){
			for (int x=0;x<woi.width;x++){
				if ((i1>canvas.length) ||(i0>=selection.length)){
					System.out.println("canvas.length="+canvas.length+" width="+width+" selection.length="+selection.length+" y="+y+" x="+x+" i0="+i0+" i1="+i1+
							" woi.x="+woi.x+" woi.y="+woi.y+" woi.width="+woi.width+" woi.height="+woi.height);
				}
				canvas[i1++]=selection[i0++]; // OOB 18720

			}
			i1+=(width-woi.width);
		}
		return canvas;
	}

	public double [] combineWithCanvas(
			double canvasFill,
			int width,
			int height,
			Rectangle woi,
			double [] selection ){
		double []canvas=new double[width*height];
		for (int i=0;i<canvas.length;i++)canvas[i]= canvasFill;
		return combineWithCanvas(canvas, width, woi, selection );
	}
	public double [] combineWithCanvas(
			double [] canvas,
			int width,
			Rectangle woi,
			double [] selection ){
		if (selection.length!=(woi.width*woi.height)) throw new IllegalArgumentException ("selection.length="+selection.length+", woi.width="+woi.width+", woi.height="+woi.height);
		int i0=0;
		int i1=width*woi.y+woi.x;
		for (int y=0;y<woi.height;y++){
			for (int x=0;x<woi.width;x++) canvas[i1++]=selection[i0++];
			i1+=(width-woi.width);
		}
		return canvas;
	}
	//===================== Moved from Aberration_Calibration
	public float[][] simulateGridAll (
			int width, // extend to full image, width, height - original (not scaled) image size
			int height,
			MatchSimulatedPattern matchSimulatedPattern,
			//				double [][][][] patternGrid, // should be aligned to gridFrac
			int gridFrac, // number of grid steps per pattern full period
			SimulParameters  simulParameters,
			int       threadsMax,
			boolean   updateStatus,
			int globalDebugLevel,
			int debug_level){// debug level used inside loops
		//			SimulationPattern simulationPattern=new SimulationPattern(simulParameters);
		float [][] simArray0=simulateGridAll (
				matchSimulatedPattern,
				//					patternGrid, // should be aligned to gridFrac
				gridFrac, // number of grid steps per pattern full period
				simulParameters,
				//					simulationPattern,
				threadsMax,
				updateStatus,
				globalDebugLevel,
				debug_level);
		Rectangle woi=matchSimulatedPattern.getWOI();
		if ((woi.x==0) && (woi.y==0) && (woi.width==width) && (woi.height==height)) return simArray0;
		int k=simulParameters.subdiv/2;
		Rectangle scaledWoi=new Rectangle(k*woi.x, k*woi.y, k*woi.width, k*woi.height);
		float [][] simArray=new float [2][];
		simArray[0]=(new SimulationPattern(simulParameters)).combineWithCanvas(0.0,  k*width, k*height, scaledWoi,simArray0[0]);
		simArray[1]=(new SimulationPattern(simulParameters)).combineWithCanvas(0.0,  k*width, k*height, scaledWoi,simArray0[1]);
		if (globalDebugLevel>1) SDFA_INSTANCE.showArrays(simArray,width*k,height*k,true, "full-simulation");
		return simArray;
	}

	public float[][] simulateGridAll (
			MatchSimulatedPattern matchSimulatedPattern,
			//				double [][][][] patternGrid, // should be aligned to gridFrac
			int gridFrac, // number of grid steps per pattern full period
			SimulationPattern.SimulParameters  simulParameters,
			//				SimulationPattern simulationPattern, // or null
			int       threadsMax,
			boolean   updateStatus,
			int globalDebugLevel,
			int debug_level){// debug level used inside loops
		long 	  startTime=System.nanoTime();
		double [][] xy0={{simulParameters.offsetX,simulParameters.offsetY},{simulParameters.offsetX-0.5,simulParameters.offsetY-0.5}} ;
		//			if (simulationPattern==null) simulationPattern=new SimulationPattern(simulParameters);
		float[][] simArray=new float[2][];
		simArray[0]=  simulateGrid (
				matchSimulatedPattern.getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				simulParameters,
				matchSimulatedPattern.getWOI(),
				simulParameters.subdiv/2,
				xy0[0],    // add to patterGrid xy
				threadsMax,
				updateStatus,
				debug_level); // debug level
		simArray[1]=  simulateGrid (
				matchSimulatedPattern.getDArray(),
				2, // gridFrac, // number of grid steps per pattern full period
				simulParameters,
				matchSimulatedPattern.getWOI(),
				simulParameters.subdiv/2,
				xy0[1],    // add to patterGrid xy
				threadsMax,
				updateStatus,
				debug_level); // debug level
		if (globalDebugLevel>2) SDFA_INSTANCE.showArrays(simArray,matchSimulatedPattern.getWOI().width*simulParameters.subdiv/2,matchSimulatedPattern.getWOI().height*simulParameters.subdiv/2,true, "a-simulation");
		if (globalDebugLevel>1) System.out.println("Grid simulation is finished at "+ IJ.d2s(0.000000001*(System.nanoTime()-startTime),3));
		return simArray;
	}

	//========================
	public float [] simulateGrid (
			final double [][][][] patternGrid, // should be aligned to gridFrac
			final int gridFrac, // number of grid steps per pattern full period: black+white
			final SimulParameters  simulParameters, // Try to use null here for maskOnly
			final Rectangle woi,
			final int       subdiv,    // subdivide output array from woi (normally 2)
			double[]  shift_xy,    // add to patterGrid xy, null OK
			//			   final boolean   maskOnly, // just mark defined cells
			final int       threadsMax,
			final boolean   updateStatus,
			final int debug_level){// debug level used inside loops
		double []xy_zero={0.0,0.0};
		if (patternGrid==null) return null;
		final double [] xy0=(shift_xy==null)?xy_zero:shift_xy;

		if ((simulParameters!=null) && (this.bPattern==null)){
			System.out.println("simulateGrid(), running patternGenerator(simulParameters )");
			patternGenerator(simulParameters ); // generate bPattern if it was not done yet
		}
		final Rectangle woiOut=new Rectangle(subdiv*woi.x,subdiv*woi.y,subdiv*woi.width,subdiv*woi.height);
		final float [] pixels=new float[woiOut.width*woiOut.height];
		final float [] pixelsDenom=new float[woiOut.width*woiOut.height];
		for (int i=0;i<pixels.length;i++ ){
			pixels[i]= 0.0F;
			pixelsDenom[i]= 0.0F;
		}
		final Thread[] threads = newThreadArray(threadsMax);
		final AtomicInteger cellNum = new AtomicInteger(0);
		final int [] series =  new int[1];
		final int uvhwidth=(patternGrid[0].length+1)/2;
		final int uvwidth=  patternGrid[0].length;
		final int uvhheight=(patternGrid.length+1)/2;
		final int numInSeries=(uvhwidth-1)*(uvhheight-1);
		final AtomicInteger debugCellNum = new AtomicInteger(0);
		final AtomicInteger finishedAtomic = new AtomicInteger(1);
		final int cellsToProcess=numInSeries*4;
		final int debugCellNum0=0;
		IJ.showStatus("Generating simulated pattern...");
		for (series[0]=0;series[0]<4;series[0]++) { // split processing in 4 series (odd/even row/column to avoid races between threads
			if (debug_level>2)System.out.println("**** series[0]="+series[0]);
			cellNum.set(0);
			for (int ithread = 0; ithread < threads.length; ithread++) {
				threads[ithread] = new Thread() {
					@Override
					public void run() {
						//	    					 String dbgStr="";
						int [][][] iUV=new int [2][2][2];
						double [][][] dUV=new double [2][2][2];
						double [][][] xy=new double [2][2][2];
						boolean invPattern;
						double [] pixDXY={1.0,1.0};
						for (int ncell=cellNum.getAndIncrement(); ncell<numInSeries;ncell=cellNum.getAndIncrement()){
							iUV[0][0][0]=2*(ncell%(uvhwidth-1))+ (series[0]      & 1);
							iUV[0][0][1]=2*(ncell/(uvhwidth-1))+ ((series[0]>>1) & 1);
							if ((updateStatus) && (debugLevel>1)) IJ.showStatus("Generating simulated pattern, series "+series[0]+" (of 4), row "+(iUV[0][0][1]/2+1)+"(of "+(uvhheight-1)+")");
							if (debugLevel>2) System.out.println("Generating pattern, series "+series[0]+" (of 4), row "+(iUV[0][0][1]/2+1)+"(of "+(uvhheight-1)+")");
							iUV[0][1][0]=iUV[0][0][0]+1;
							iUV[0][1][1]=iUV[0][0][1];
							iUV[1][0][0]=iUV[0][0][0];
							iUV[1][0][1]=iUV[0][0][1]+1;
							iUV[1][1][0]=iUV[0][0][0]+1;
							iUV[1][1][1]=iUV[0][0][1]+1;
							if ((isCellValid(patternGrid,iUV[0][0])) &&
									(isCellValid(patternGrid,iUV[0][1])) &&
									(isCellValid(patternGrid,iUV[1][0])) &&
									(isCellValid(patternGrid,iUV[1][1]))){



								// All 4 corners are valid
								invPattern=((iUV[0][0][0]%gridFrac)>=(gridFrac/2))^((iUV[0][0][1]%gridFrac)>=(gridFrac/2));
								if (debug_level>2)System.out.println("iUV[0][0][1]="+iUV[0][0][1]+" iUV[0][0][0]="+iUV[0][0][0]+"    invert="+invPattern);


								for (int i=0;i<2;i++) for (int j=0;j<2;j++) for (int k=0;k<2;k++) {
									xy[i][j][k]=subdiv*patternGrid[iUV[i][j][1]][iUV[i][j][0]][0][k]+xy0[k];
								}
								for (int k=0;k<2;k++) {
									dUV[0][0][k]=((double) (iUV[0][0][k]%(gridFrac/2)))/(gridFrac/2);
								}
								dUV[0][1][0]=dUV[0][0][0]+1.0/(gridFrac/2);
								dUV[0][1][1]=dUV[0][0][1];
								dUV[1][0][0]=dUV[0][0][0];
								dUV[1][0][1]=dUV[0][0][1]+1.0/(gridFrac/2);
								dUV[1][1][0]=dUV[0][1][0];
								dUV[1][1][1]=dUV[1][0][1];

								double [] minXY={xy[0][0][0],xy[0][0][1]};
								double [] maxXY={xy[0][0][0],xy[0][0][1]};
								for (int i=0;i<2;i++) for (int j=0;j<2;j++) for (int k=0;k<2;k++) {
									if (minXY[k]>xy[i][j][k]) minXY[k]=xy[i][j][k];
									if (maxXY[k]<xy[i][j][k]) maxXY[k]=xy[i][j][k];
								}
								Rectangle rcell=new Rectangle((int)minXY[0], // contains all pixels
										(int)minXY[1],
										((int) Math.ceil(maxXY[0]))-((int)minXY[0]) ,
										((int) Math.ceil(maxXY[1]))-((int)minXY[1]));
								double [] cornersDUV={1.0/(gridFrac/2),1.0/(gridFrac/2)};
								boolean debugNow=debugCellNum.getAndIncrement()==debugCellNum0;
								if (woiOut.intersects (rcell)) // do not bother if no
									for (int iy=rcell.y;iy<(rcell.y+rcell.height);iy++) for (int ix=rcell.x;ix<(rcell.x+rcell.width);ix++)
										if (woiOut.contains (ix,iy))
										{
											double [] pixXY={ix,iy};
											double [] pixData=recursiveFillPixels ( // invert pattern in the caller, return signed value (-1..1 - pattern is 0..1)
													simulParameters,
													pixXY,  // top-left corner
													pixDXY, // increments to other corners
													xy, // xy pairs for the 4 corners of the square in UV (pattern) coordinates (u0v0,u1v0,u0v1,u1v1)
													dUV[0][0],  // UV value for the top-left corner (matching cornersXY[0][0])
													cornersDUV,  // distances to the opposite corner in UV
													//	    												  maskOnly, // just mark defined cells
													debugNow?debug_level: debug_level-2
													);
											//	    										 int index=woiOut.width*iy+ix;
											int index=woiOut.width*(iy-woiOut.y)+(ix-woiOut.x);
											//	    										 pixels[index]+=((invPattern || maskOnly)?1.0:-1.0)*pixData[0];
											if (index>pixels.length){
												//											         final float [] pixels=new float[woiOut.width*woiOut.height];

												System.out.println("simulateGrid(), pixels.length="+pixels.length+
														" index="+index+
														" iy="+iy+" ix="+ix+
														" ncell="+ncell+
														" subdiv="+subdiv+
														" woiOut.x="+woiOut.x+
														" woiOut.y="+woiOut.y+
														" woiOut.width="+woiOut.width+
														" woiOut.height="+woiOut.height+
														" woi.x="+woi.x+
														" woi.y="+woi.y+
														" woi.width="+woi.width+
														" woi.height="+woi.height+
														" rcell.x="+rcell.x+
														" rcell.y="+rcell.y+
														" rcell.width="+rcell.width+
														" rcell.height="+rcell.height);

											}
											//	    										 if (maskOnly) pixels[index]=patternGrid.length*iUV[0][0][1]+iUV[0][0][0]; // out of bounds
											// OLD NASTY BUG!
											//	    										 if (maskOnly) pixels[index]=uvwidth*iUV[0][0][1]+iUV[0][0][0]; // Should be width, not height!
											if (simulParameters==null) pixels[index]=uvwidth*iUV[0][0][1]+iUV[0][0][0]; // Should be width, not height!
											else          pixels[index]+=(invPattern?1.0:-1.0)*pixData[0];

											pixelsDenom[index]+=pixData[1];
										}
							}
							final int numFinished=finishedAtomic.getAndIncrement();
							SwingUtilities.invokeLater(new Runnable() {
								@Override
								public void run() {
									IJ.showProgress(numFinished,cellsToProcess);
								}
							});

						}
					}
				};
			}
			startAndJoin(threads);
		}
		for (int i=0;i<pixels.length;i++ ) {
			//				/(simulParameters!=null)
			//				if (maskOnly) {
			if (simulParameters==null) {
				if (pixelsDenom[i]==0.0F) pixels[i]=-1;
			} else {
				if (pixelsDenom[i]!=0.0F){
					pixels[i]/= pixelsDenom[i];
					pixels[i]=(float) ((pixels[i]+1.0)/2); // convert from -1..+1 to 0..1.0
				}
			}
		}

		return pixels;
	}
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
	private static void startAndJoin(Thread[] threads)
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

	/* ======================================================================== */
	/* make it faster when outSubdiv =2*n (usually so) */
	/* TODO: cleanup shifts - they seem now to work correctly */
	public double [][] extractSimulPatterns (
			SimulParameters  simulParameters,
			int outSubdiv,  // subdivide output pixels
			int size,    // number of Bayer cells in width of the square selection (half number of pixels)
			double x0,    // selection center, X (in pixels)
			double y0) {
		return extractSimulPatterns (
				this.barray,
				simulParameters,
				outSubdiv,  // subdivide output pixels
				size,    // number of Bayer cells in width of the square selection (half number of pixels)
				x0,    // selection center, X (in pixels)
				y0);
	}

	public double [][] extractSimulPatterns (
			double [] localbArray,
			SimulParameters  simulParameters,
			int outSubdiv,  // subdivide output pixels
			int size,    // number of Bayer cells in width of the square selection (half number of pixels)
			double x0,    // selection center, X (in pixels)
			double y0) {
		int sampleWidth=(int) (Math.sqrt(simulParameters.fill)*simulParameters.subdiv);
		int sampleN=sampleWidth*sampleWidth;
		if      (sampleWidth<1)     sampleWidth=1;
		else if (sampleWidth>simulParameters.subdiv)sampleWidth=simulParameters.subdiv;
		double sampleAverage=0.5*sampleN;

		int n,i,j;
		int fullSize=(int) Math.sqrt(localbArray.length);
		double [][] simul_pixels=new double [5][size*size];
		int ix,iy, iy0,ix0,px,py;
		double bx,by;
		double s;
		double span=((double) size)/outSubdiv;
		int sampLow=-sampleWidth/2;
		int sampHigh=sampLow+sampleWidth;
		for (n=0;n<4;n++) {
			bx=(n&1)-0.5+0.5;        // last 0.5 to make same center as for dual greens
			by=((n>>1) & 1)-0.5-0.5;// last 0.5 to make same center as for dual greens
			for (iy=0;iy<size;iy++) {
				iy0=(fullSize/2) + (int) ((-span+y0+by  +1.5  +2.0*iy/outSubdiv)*simulParameters.subdiv);
				for (ix=0;ix<size;ix++) {
					ix0=(fullSize/2) + (int) ((-span+x0+bx+0.5   +2.0*ix/outSubdiv)*simulParameters.subdiv);
					s=0.0;
					for (py=iy0+sampLow;py<iy0+sampHigh;py++) for (px=ix0+sampLow;px<ix0+sampHigh;px++) {
						try {
							s+=localbArray[py*fullSize+px];
						} catch (Exception e){
							System.out.println("Bug in extractSimulPatterns(): px="+px+" py="+py+" fullSize="+fullSize+" size="+size+" x0="+x0+" y0="+y0);
							e.printStackTrace();
							return null;
						}
					}
					simul_pixels[n][iy*size+ix]= (s-sampleAverage)/sampleAverage;
				}
			}

		}
		if (outSubdiv>1) {
			if (this.debugLevel>2)System.out.println("Generating combined greens pattern greens from scratch");
			n=4;
			bx=0.0;
			by=0.0;
			for (iy=0;iy<size;iy++) {
				for (ix=0;ix<size;ix++) {
					iy0=(fullSize/2) + (int) ((-span+y0+by-1+1.5  +1.0*(size+iy-ix)/outSubdiv)*simulParameters.subdiv);
					ix0=(fullSize/2) + (int) ((-span+x0+bx  +0.5  +1.0*(iy+ix)/outSubdiv)*simulParameters.subdiv);
					s=0.0;
					for (py=iy0+sampLow;py<iy0+sampHigh;py++) for (px=ix0+sampLow;px<ix0+sampHigh;px++) {
						s+=localbArray[py*fullSize+px];
					}
					simul_pixels[n][iy*size+ix]= (s-sampleAverage)/sampleAverage;
				}

			}
		} else { // just reuse available greens
			if (this.debugLevel>2)System.out.println("Generating combined greens pattern from individual greens");
			/* now combine greens - same as in splitBayer() */

			int base, base_b;
			base_b=0;
			for (i=0;i<size/2; i++){
				base=size*size/2+ i* (size+1);
				for (j=0; j<size/2; j++) {
					simul_pixels[4][base_b++]=simul_pixels[0][base];
					base-=size;
					simul_pixels[4][base_b++]=simul_pixels[3][base++];
				}
				base=size*size/2+ i* (size+1);
				for (j=0; j<size/2; j++) {
					//System.out.println("2:y="+y+" x="+x+" base_b="+base_b+" base="+base);
					simul_pixels[4][base_b++]=simul_pixels[3][base++];
					simul_pixels[4][base_b++]=simul_pixels[0][base];
					base-=size;
				}
			}
		}
		if (this.debugLevel>2) {
			System.out.println("extractSimulPatterns, x0="+x0+" y0="+y0+" fullSize="+fullSize+" size="+size+" subdiv="+simulParameters.subdiv+" outSubdiv="+outSubdiv);
			System.out.println(" sampLow="+sampLow+" sampHigh="+sampHigh+" span="+span+" size="+size);
			for (n=0;n<simul_pixels.length;n++) {
				s=0.0;
				for (i=0;i<simul_pixels[n].length;i++) s+=simul_pixels[n][i];
				System.out.println(" component="+i+" sum of pixels="+s);
			}
		}

		if (this.debugLevel>2) SDFA_INSTANCE.showArrays(simul_pixels,size,size, "SIMUL");

		return simul_pixels;
	}

	public double [] extractSimulMono (
			SimulParameters  simulParameters,
			int outSubdiv,  // subdivide output pixels
			int size,    // number of Bayer cells in width of the square selection (half number of pixels)
			double x0,    // selection center, X (in pixels)
			double y0) {
		return extractSimulMono (
				this.barray,
				simulParameters,
				outSubdiv,  // subdivide output pixels
				size,    // number of Bayer cells in width of the square selection (half number of pixels)
				x0,    // selection center, X (in pixels)
				y0);
	}


	public double [] extractSimulMono ( // TODO: can use twice smaller barray
			double [] localbArray,
			SimulParameters  simulParameters,
			int outSubdiv,  // subdivide output pixels - now 4
			int size,    // number of Bayer cells in width of the square selection (half number of pixels)
			double x0,    // selection center, X (in pixels)
			double y0) {
		int sampleWidth=(int) (Math.sqrt(simulParameters.fill)*simulParameters.subdiv);
		int sampleN=sampleWidth*sampleWidth;
		if      (sampleWidth<1)     sampleWidth=1;
		else if (sampleWidth>simulParameters.subdiv)sampleWidth=simulParameters.subdiv;
		double sampleAverage=0.5*sampleN;
		int fullSize=(int) Math.sqrt(localbArray.length);
		double [] simul_pixels=new double [size*size];
		int ix,iy, iy0,ix0,px,py;
		double s;
		double span=((double) size)/outSubdiv;
		int sampLow=-sampleWidth/2;
		int sampHigh=sampLow+sampleWidth;
		for (iy=0;iy<size;iy++) {
			iy0=(fullSize/2) + (int) ((    -span/2 + y0 + 0.5 + 1.0 * iy / outSubdiv)*simulParameters.subdiv);
			for (ix=0;ix<size;ix++) {
				ix0=(fullSize/2) + (int) ((-span/2 + x0 + 0.5 + 1.0 * ix / outSubdiv)*simulParameters.subdiv);
				s=0.0;
				for (py=iy0+sampLow;py<iy0+sampHigh;py++) for (px=ix0+sampLow;px<ix0+sampHigh;px++) {
					try {
						s+=localbArray[py*fullSize+px];
					} catch (Exception e){
						System.out.println("Bug in extractSimulMono(): px="+px+" py="+py+" fullSize="+fullSize+" size="+size+" x0="+x0+" y0="+y0);
						e.printStackTrace();
						return null;
					}
				}
				simul_pixels[iy*size+ix]= (s-sampleAverage)/sampleAverage;
			}
		}
		if (this.debugLevel>2) {
			System.out.println("extractSimulMono, x0="+x0+" y0="+y0+" fullSize="+fullSize+" size="+size+" subdiv="+simulParameters.subdiv+" outSubdiv="+outSubdiv);
			System.out.println(" sampLow="+sampLow+" sampHigh="+sampHigh+" span="+span+" size="+size);
			s=0.0;
			for (int i=0;i<simul_pixels.length;i++) s+=simul_pixels[i];
			System.out.println("sum of pixels="+s);
		}
		if (this.debugLevel>2) SDFA_INSTANCE.showArrays(simul_pixels,size,size, "SIMUL");
		return simul_pixels;
	}

	private double [] matrix2x2_mul(double [][] a, double [] b ){
		double [] rslt={a[0][0]*b[0]+a[0][1]*b[1],
				a[1][0]*b[0]+a[1][1]*b[1]};
		return rslt;
	}
	private double []    vector_add(double [] a, double [] b ){
		double [] rslt= new double [a.length];
		int i;
		for (i=0;i<rslt.length;i++) rslt[i]=a[i]+b[i];
		return rslt;
	}

	public double[] extractBayerSim (
			float [][] spixels, // [0] - regular pixels, [1] - shifted by 1/2 diagonally, for checker greens
			int full_width,
			Rectangle woi,
			int bayerPeriod, // 4
			int colorComp) {
		Rectangle r=new Rectangle(woi); // clone
		int full_height=spixels[0].length/full_width; // full image height
		if (debugLevel>10) IJ.showMessage("splitBayer","r.width="+r.width+
				"\nr.height="+r.height+
				"\nr.x="+r.x+
				"\nr.y="+r.y+
				"\nlength="+spixels[0].length);
		if ((debugLevel>2) && ((r.x<0) || (r.y<0) || ((r.x+r.width)>=full_width) || ((r.y+r.height)>=full_height))) System.out.println("r.width="+r.width+
				" r.height="+r.height+
				" r.x="+r.x+
				" r.y="+r.y);
		if (colorComp==5) colorComp=0; // for compatibility, combined grees and green 0 generate the same result
		double []result=new double[r.width*r.height];
		int index;
		if (colorComp==4) { // checkerboard greens
			r.y+=r.width/2; // now it is the "top left" corner of the diagonal greens
			for (index=0;index<result.length;index++){
				//				int iy=r.y+(index / r.width);
				//				int ix=r.x+(index % r.width);
				int iyi=index / r.width;
				int ixi=index % r.width;
				int iy=r.y+(iyi+ixi)/2 -ixi;
				int ix=r.x+(iyi+ixi)/2;
				if (iy<0) iy=0;
				else if (iy>=full_height) iy=full_height-1;
				if (ix<0) ix=0;
				else if (ix>=full_width) iy=full_width-1;
				result[index]=spixels[(ixi+iyi) & 1][iy*full_width+ix];
			}
		} else { // components 0..3
			r.x+=(bayerPeriod/2)*(colorComp &1);
			r.y+=(bayerPeriod/2)*((colorComp>>1) &1);
			if (debugLevel>2) System.out.println(">>> r.width="+r.width+
					" r.height="+r.height+
					" r.x="+r.x+
					" r.y="+r.y+
					" colorComp="+colorComp);
			for (index=0;index<result.length;index++){
				int iy=r.y+(index / r.width);
				int ix=r.x+(index % r.width);
				if (iy<0) iy=0;
				else if (iy>=full_height) iy=full_height-1;
				if (ix<0) ix=0;
				else if (ix>=full_width) iy=full_width-1;
				result[index]=spixels[0][iy*full_width+ix];
			}
		}
		return result;
	}

	//=====================
	public static class SimulParameters {
		public int    patternSize;
		public int    pattern_type;
		public double pattern_modifier;
		public double freq_x1;
		public double freq_y1;
		public double phase1;
		public double freq_x2;
		public double freq_y2;
		public double phase2;
		public int    subdiv;
		public double fill;
		public boolean center_for_g2;
		public double bPatternSigma; // blur bPattern with this sigma
		public double barraySigma; // blur barray with this sigma, multiplied by subdiv
		public double smallestSubPix; // subdivide pixels down to that fraction when simulating
		public double bitmapNonuniforityThreshold; // subdivide pixels until difference between the corners is below this value
		public double offsetX; // debug - add to X during simulation, in pixels
		public double offsetY; // debug - add to Y during simulation, in pixels


		public SimulParameters(
				int    patternSize,
				int    pattern_type,
				double pattern_modifier,
				double freq_x1,
				double freq_y1,
				double phase1,
				double freq_x2,
				double freq_y2,
				double phase2,
				int    subdiv,
				double fill,
				boolean center_for_g2,
				double bPatternSigma, // blur bPattern with this sigma
				double barraySigma, // blur barray with this sigma, multiplied by subdiv
				double smallestSubPix, // subdivide pixels down to that fraction when simulating
				double bitmapNonuniforityThreshold, // subdivide pixels until difference between the corners is below this value
				double offsetX, // debug - add to X during simulation, in pixels
				double offsetY // debug - add to Y during simulation, in pixels
				) {

			this.patternSize=    patternSize;
			this.pattern_type=    pattern_type;
			this.pattern_modifier=pattern_modifier;
			this.freq_x1=         freq_x1;
			this.freq_y1=         freq_y1;
			this.phase1=          phase1;
			this.freq_x2=         freq_x2;
			this.freq_y2=         freq_y2;
			this.phase2=          phase2;
			this.subdiv=          subdiv;
			this.fill=            fill;
			this.center_for_g2=   center_for_g2;
			this.bPatternSigma=   bPatternSigma;
			this.barraySigma=barraySigma;
			this.smallestSubPix=  smallestSubPix;
			this.bitmapNonuniforityThreshold=bitmapNonuniforityThreshold;
			this.offsetX=         offsetX;
			this.offsetY=         offsetY;

		}
		@Override
		public SimulParameters clone() {
			return new SimulParameters(
					this.patternSize,
					this.pattern_type,
					this.pattern_modifier,
					this.freq_x1,
					this.freq_y1,
					this.phase1,
					this.freq_x2,
					this.freq_y2,
					this.phase2,
					this.subdiv,
					this.fill,
					this.center_for_g2,
					this.bPatternSigma,
					this.barraySigma,
					this.smallestSubPix,
					this.bitmapNonuniforityThreshold,
					this.offsetX,
					this.offsetY
					);
		}

		public void setProperties(String prefix,Properties properties){
			properties.setProperty(prefix+"patternSize",this.patternSize+"");
			properties.setProperty(prefix+"pattern_type",this.pattern_type+"");
			properties.setProperty(prefix+"pattern_modifier",this.pattern_modifier+"");
			properties.setProperty(prefix+"freq_x1",this.freq_x1+"");
			properties.setProperty(prefix+"freq_y1",this.freq_y1+"");
			properties.setProperty(prefix+"phase1",this.phase1+"");
			properties.setProperty(prefix+"freq_x2",this.freq_x2+"");
			properties.setProperty(prefix+"freq_y2",this.freq_y2+"");
			properties.setProperty(prefix+"phase2",this.phase2+"");
			properties.setProperty(prefix+"subdiv",this.subdiv+"");
			properties.setProperty(prefix+"fill",this.fill+"");
			properties.setProperty(prefix+"center_for_g2",this.center_for_g2+"");
			properties.setProperty(prefix+"bPatternSigma",this.bPatternSigma+"");
			properties.setProperty(prefix+"barraySigma",this.barraySigma+"");
			properties.setProperty(prefix+"smallestSubPix",this.smallestSubPix+"");
			properties.setProperty(prefix+"bitmapNonuniforityThreshold",this.bitmapNonuniforityThreshold+"");
			properties.setProperty(prefix+"offsetX",this.offsetX+"");
			properties.setProperty(prefix+"offsetY",this.offsetY+"");
		}
		public void getProperties(String prefix,Properties properties){
			if (properties.getProperty(prefix+"patternSize")!=null) this.patternSize=Integer.parseInt(properties.getProperty(prefix+"patternSize"));
			if (properties.getProperty(prefix+"pattern_type")!=null) this.pattern_type=Integer.parseInt(properties.getProperty(prefix+"pattern_type"));
			if (properties.getProperty(prefix+"pattern_modifier")!=null) this.pattern_modifier=Double.parseDouble(properties.getProperty(prefix+"pattern_modifier"));
			if (properties.getProperty(prefix+"freq_x1")!=null) this.freq_x1=Double.parseDouble(properties.getProperty(prefix+"freq_x1"));
			if (properties.getProperty(prefix+"freq_y1")!=null) this.freq_y1=Double.parseDouble(properties.getProperty(prefix+"freq_y1"));
			if (properties.getProperty(prefix+"phase1")!=null) this.phase1=Double.parseDouble(properties.getProperty(prefix+"phase1"));
			if (properties.getProperty(prefix+"freq_x2")!=null) this.freq_x2=Double.parseDouble(properties.getProperty(prefix+"freq_x2"));
			if (properties.getProperty(prefix+"freq_y2")!=null) this.freq_y2=Double.parseDouble(properties.getProperty(prefix+"freq_y2"));
			if (properties.getProperty(prefix+"phase2")!=null) this.phase2=Double.parseDouble(properties.getProperty(prefix+"phase2"));
			if (properties.getProperty(prefix+"subdiv")!=null) this.subdiv=Integer.parseInt(properties.getProperty(prefix+"subdiv"));
			if (properties.getProperty(prefix+"fill")!=null) this.fill=Double.parseDouble(properties.getProperty(prefix+"fill"));
			if (properties.getProperty(prefix+"center_for_g2")!=null) this.center_for_g2=Boolean.parseBoolean(properties.getProperty(prefix+"center_for_g2"));
			if (properties.getProperty(prefix+"bPatternSigma")!=null) this.bPatternSigma=Double.parseDouble(properties.getProperty(prefix+"bPatternSigma"));
			if (properties.getProperty(prefix+"barraySigma")!=null) this.barraySigma=Double.parseDouble(properties.getProperty(prefix+"barraySigma"));
			if (properties.getProperty(prefix+"smallestSubPix")!=null) this.smallestSubPix=Double.parseDouble(properties.getProperty(prefix+"smallestSubPix"));
			if (properties.getProperty(prefix+"bitmapNonuniforityThreshold")!=null) this.bitmapNonuniforityThreshold=Double.parseDouble(properties.getProperty(prefix+"bitmapNonuniforityThreshold"));
			if (properties.getProperty(prefix+"offsetX")!=null) this.offsetX=Double.parseDouble(properties.getProperty(prefix+"offsetX"));
			if (properties.getProperty(prefix+"offsetY")!=null) this.offsetY=Double.parseDouble(properties.getProperty(prefix+"offsetY"));
		}

	}
}
