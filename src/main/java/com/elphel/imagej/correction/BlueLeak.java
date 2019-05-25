package com.elphel.imagej.correction;
/**
** -----------------------------------------------------------------------------**
** BlueLeak.java
**
** Mitigating lens poor performance for blue light, especially noticeable in
** the dark areas near overexposed ones.
**
** Copyright (C) 2014 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  BlueLeak.java is free software: you can redistribute it and/or modify
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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.elphel.imagej.cameras.EyesisCorrectionParameters;
import com.elphel.imagej.common.DoubleGaussianBlur;
import com.elphel.imagej.common.ShowDoubleFloatArrays;

public class BlueLeak {
	private double [][] rgb_in;
	private double [] yrg;
	private double [] blue;
	private int width;
	private int height;
	private int length;
	EyesisCorrectionParameters.ColorProcParameters  colorProcParameters;
	private ShowDoubleFloatArrays SDFA_INSTANCE;
	private int [] dirs20;
	private double [] weights20;
	private int [] dirs8;
	private int [] dirs4;
	private int debugLevel;
	private String dbgImgTitle;
	private int margin=50; // make configurable?
	private int xMax,yMax;
	
	public BlueLeak(
			EyesisCorrectionParameters.ColorProcParameters  colorProcParameters,
			double [][] rgb,
			int width,
			ShowDoubleFloatArrays SDFA_INSTANCE,
			String dbgImgTitle,
			int debugLevel){
		length=rgb[0].length;
		rgb_in=rgb;
		this.width=width;
		height=length/width;
		this.SDFA_INSTANCE=SDFA_INSTANCE;
		this.colorProcParameters=colorProcParameters;
		if (debugLevel>2) System.out.println("width="+width+" height="+height+" length="+length);
		int [] dirs8 = {-width,-width+1, 1, width+1, width,width-1,-1,-width-1};
		this.dirs8=dirs8;
		int [] dirs4 = {-width, 1,  width,-1};
		this.dirs4=dirs4;
		int [] dirs20 = {
				-width,-width+1, 1, width+1, width,width-1,-1,-width-1,
				-2*width-1,-2*width,-2*width+1,
				-width+2, 2, width+2,
				2*width-1,2*width,2*width+1,
				-width-2, -2, width-2
				};
		this.dirs20=dirs20;
		double [] weights20={
				1.0, 0.75, 1.0, 0.75, 1.0, 0.75, 1.0, 0.75,
				0.5, 0.5, 0.5, 0.5,
				0.5, 0.5, 0.5, 0.5,
				0.5, 0.5, 0.5, 0.5,
				0.5, 0.5, 0.5, 0.5};
		this.weights20=weights20;
		this.debugLevel=debugLevel; //+1;  /*** CHANGE later ***/
		this.dbgImgTitle=dbgImgTitle;
		this.xMax=width-margin;
		this.yMax=height-margin;
		this.yrg=getYrg();

		
	}

	public enum BLUE_STATE{
		UNDEF,
		OVEREXP,
		EXPANDED,
		PROBLEM,
		SOLUTION,
		SOLVED,
		TMP
	}

	private BLUE_STATE[] markBlueOverexp(double [] blueOverExp){
		BLUE_STATE[] result= new BLUE_STATE[blueOverExp.length];
		for (int i=0;i<result.length;i++){
			result[i]=Double.isNaN(blueOverExp[i])?BLUE_STATE.UNDEF:BLUE_STATE.OVEREXP;
		}
		return result;
	}
	
	private List<Integer> createInitialList(
			BLUE_STATE[] state,
			BLUE_STATE thisState,
			BLUE_STATE neighborState,
			boolean use8
			){
		List <Integer> pixelList=new ArrayList<Integer>();
		int [] dirs=use8?dirs8:dirs4;
		for (int i=0;i<length;i++){
			if (state[i]!=thisState) continue;
			for (int dir=0; dir<dirs.length; dir++){
				int index=i+dirs[dir];
				int y=index/width;
				if ((y<margin) || (y>=yMax)) continue;
				int x=index%width;
				if ((x<margin) || (x>=xMax)) continue;
				if (state[index]==neighborState){
					pixelList.add(new Integer(i));
					break;
				}
			}
		}
		return pixelList;
	}

	private double[] findBlueSolutions(
			double [] blueOverExp,
			boolean blueLeakNoHint,
			boolean blueLeakNoBrighten,
			boolean use8
			){
		int [] dirs=use8?dirs8:dirs4;
		BLUE_STATE[] state= markBlueOverexp(blueOverExp);
		List<Integer> pixelList= createInitialList(
				state,
				BLUE_STATE.UNDEF,
				BLUE_STATE.OVEREXP,
				use8);
		
		// Shrink overexposed areas by blueOverShrink steps, mark them as UNDEF
		for (int nShrink=0;nShrink<colorProcParameters.blueOverShrink; nShrink++){
			int len=pixelList.size();
			for (int i=0;i<len;i++){
				int oldIndex=pixelList.remove(0);
				for (int dir=0; dir<dirs.length; dir++){
					int index=oldIndex+dirs[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if (state[index]==BLUE_STATE.OVEREXP){
						pixelList.add(new Integer(index));
						state[index]=BLUE_STATE.UNDEF;
					}
				}
			}
			if (debugLevel>2) System.out.println("Shrinking blue overexposed area, step "+nShrink+", number of pixels in the list: "+pixelList.size());
		}
		double [] yRef=new double [length];
		double [] leak=new double [length];
//		double [] yrg=getYrg();
		for (int i=0;i<length;i++){
			yRef[i] = ((state[i]==BLUE_STATE.OVEREXP) || (state[i]==BLUE_STATE.EXPANDED)) ?yrg[i]:Double.NaN;
			leak[i] = yRef[i];
		}
		List<Integer> solutionList=new ArrayList<Integer>();
		int problemPixels=0;
		double rate1=1.0/colorProcParameters.blueBandWidth;
		double rate0=1.0/colorProcParameters.blueBandWidthDark;
		
		if (debugLevel>1){
			String [] channel_blue_titles={"orig. blue","yRef","state"};
			double [][] dbgImg= new double [3][];
			dbgImg[0]=rgb_in[2];
			dbgImg[1]=yRef;
			dbgImg[2]=dbgStateToDouble(state);
			SDFA_INSTANCE.showArrays(dbgImg, width, height, true, "blue_overexp_expanded", channel_blue_titles);
		}
		
		int numStep=0;
		while (!pixelList.isEmpty()){
			boolean lookForSolution=numStep>=colorProcParameters.blueOverGrow; // first passes - just expand, do not look for solution
			//pixelList has pixels still UNDEFINED near OVEREXP (or problem/solution)
			if (debugLevel>2) System.out.println("Starting step "+numStep+", number of pixels in the list: "+pixelList.size()+", in solutions: "+solutionList.size());
			int [] solutionFront=new int[pixelList.size()]; // 0-problem, 1 solution, -1 - none (overexposed area shrank completely) 
			for (int i=0;i<solutionFront.length;i++){
				int oldIndex=pixelList.get(i);
				if (debugLevel>2) if (i==0) System.out.println ("For i==0 oldIndex="+oldIndex+", state["+oldIndex+"]="+state[oldIndex]);
				if (state[oldIndex]!=BLUE_STATE.UNDEF){
					System.out.println("BUG: State is not UNDEF for i="+i+", oldIndex="+oldIndex+" x="+(oldIndex%width)+", y="+oldIndex/width);
					return null; //rgb_in[2];
				}
				
				double w=0.0;
				double avgRef=0.0;
				double avgLeak=0.0;
				
				for (int dir=0; dir<dirs20.length; dir++){
					int index=oldIndex+dirs20[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if (state[index]!=BLUE_STATE.UNDEF){ // can be problem/solution or overexp?
						avgRef+= weights20[dir]*yRef[index];
						avgLeak+=weights20[dir]*leak[index];
						w+=weights20[dir];
					}
				}
				if (w>0.0){
					avgRef/=w;
					avgLeak/=w;
					yRef[oldIndex]=avgRef;
					if (lookForSolution){
//						double rate=rate0+ yrg[oldIndex]/avgRef*rate1;
						double rate=avgRef*rate0+ yrg[oldIndex]*rate1;
						leak[oldIndex]=avgLeak-rate;
						solutionFront[i]=(leak[oldIndex]<=yrg[oldIndex])?1:0;
					} else {
						leak[oldIndex]=	avgRef;
						solutionFront[i]=0;
					}
				} else { // now after shrinking it is OK
//					System.out.println("BUG: Could not find any neighbors.");
					solutionFront[i]=-1;
				}
			}
// Second pass - apply new states
			for (int i=0;i<solutionFront.length;i++){
				int index=pixelList.remove(0);
				if (solutionFront[i]>0){
					state[index]=BLUE_STATE.SOLUTION;
					solutionList.add(new Integer(index));
				} if (solutionFront[i]==0){
					state[index]=BLUE_STATE.PROBLEM;
					pixelList.add(new Integer(index));
					problemPixels++;
				} // Do nothing if <0 - keep undefined, do not add to list
			}
// Find fresh cells near pixelList members (destructive)
			int len=pixelList.size();
			for (int i=0;i<len;i++) {
				int oldIndex=pixelList.remove(0);
				for (int dir=0; dir<dirs.length; dir++){
					int index=oldIndex+dirs[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if (state[index]==BLUE_STATE.UNDEF){
						pixelList.add(new Integer(index));
						state[index]=BLUE_STATE.TMP; // will be removed before return
					}
				}
			}
// restore state to undefined			
			for (Iterator<Integer> iter= pixelList.iterator(); iter.hasNext();){
				int i=iter.next();
				assert (state[i]==BLUE_STATE.TMP):"BLUE_STATE.TMP expected, got "+state[i]+" for i="+i+", x="+(i%width)+" y="+(i/width);
				state[i]=BLUE_STATE.UNDEF;
			}
			numStep++;
		}
		if (debugLevel>2) System.out.println("Found "+problemPixels+" blue problem pixels, "+solutionList.size()+" pixels on solutions front");
		// just verify they are all OK
		for (Iterator<Integer> iter= solutionList.iterator(); iter.hasNext();){
			int index=iter.next();
			assert(state[index]==BLUE_STATE.SOLUTION):"Expected SOLUTION, got "+state[index];
		}
		
		if ((debugLevel>1)){
			String [] channel_blue_titles={"orig. blue","Yrg", "yRef","leak","state"};
			double [][] dbgImg= new double [5][];
			dbgImg[0]=rgb_in[2];
			dbgImg[1]=yrg;
			dbgImg[2]=yRef;
			dbgImg[3]=leak;
			dbgImg[4]=dbgStateToDouble(state);
			SDFA_INSTANCE.showArrays(dbgImg, width, height, true, "blue_overexp_expanded__"+numStep, channel_blue_titles);
		}
		blue=new double[length];
		for (int i=0;i<length;i++) blue[i]=rgb_in[2][i]; // just copy
//		if ((debugLevel>1)) return blue;

		// re-use yRef to hold blue/yrg ratio, deducted from solution. Actual blue/Yrg ration will gradually shift to
		// blue_neutral_ratio farther from solution to prevent colored spikes far from solution origin
		for (Iterator<Integer> iter= solutionList.iterator(); iter.hasNext();){
			int index=iter.next();
			yRef[index]=blue[index]/yrg[index];
		}
		// repeat filling problem areas
		int numSolStep=0;
		int numSolved=0;
		double fadeRate=(colorProcParameters.blueSolutionRadius>0)?(1.0/colorProcParameters.blueSolutionRadius):0.0;
		while (!solutionList.isEmpty()){
			int numSolvedPass=0;
			// Find problem pixels around, add to list, mark as TMP, add to list (while removing original) //EXPANDED
			int len=solutionList.size();
			for (int i=0;i<len;i++){
				int oldIndex=solutionList.remove(0);
				for (int dir=0; dir<dirs.length; dir++){
					int index=oldIndex+dirs[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if (state[index]==BLUE_STATE.PROBLEM){
						state[index]=BLUE_STATE.TMP;
						solutionList.add(new Integer(index));
					}
				}
			}
			if (debugLevel>2) System.out.println(numSolStep+": new front length="+solutionList.size()+", old length="+len);
			
			
			// first non-destructive pass through the list - calculate yRef[] (applying fadeRate), blue[]
			for (Iterator<Integer> iter= solutionList.iterator(); iter.hasNext();){
				int oldIndex=iter.next();
				assert(state[oldIndex]==BLUE_STATE.TMP):"expected TMP, got"+state[iter.next()];

				double w=0.0;
				double avg=0.0;
				for (int dir=0; dir<dirs20.length; dir++){
					int index=oldIndex+dirs20[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if ((state[index]==BLUE_STATE.SOLUTION) || (state[index]==BLUE_STATE.SOLVED)){
						avg+= weights20[dir]*yRef[index];
						w+=weights20[dir];
					}
				}
				if (w>0){
					avg/=w;
					avg=(1-fadeRate)*avg+fadeRate*colorProcParameters.blueNeutral;
					yRef[oldIndex]=avg;
					// blue[oldIndex]=yrg[oldIndex]*avg;
					double newBlue=yrg[oldIndex]*avg;
					if (!blueLeakNoBrighten || (blue[oldIndex] > newBlue)) blue[oldIndex]=newBlue;
					numSolvedPass++;
					numSolved++;
				} else {
					System.out.println("BUG: Could not find any neighbors when applying solution, state[oldIndex]="+state[oldIndex]); // Got here many times
				}
			}			
			// second non-destructive pass through the list - change state to "SOLVED" 
			for (Iterator<Integer> iter= solutionList.iterator(); iter.hasNext();){
				state[iter.next()]=BLUE_STATE.SOLVED;
			}
			if (debugLevel>2) System.out.println("Applying solution, step #"+(numSolStep++)+" solved on this pass "+numSolvedPass+" pixels.");
		}
		if (debugLevel>1) System.out.println("Appled solution for blue color leak, "+(numSolStep)+" steps,  solved "+numSolved+" pixels.");
		// fix small "windows" where there are no blue color reference pixels in the closed area
		//blueLeakNoBrighten
		
		if (blueLeakNoHint) {
			int numSmallFixed=0;
			for (int i=0;i<length;i++) if (state[i]== BLUE_STATE.PROBLEM){
				double newBlue=yrg[i]*colorProcParameters.blueNeutral;
				if (!blueLeakNoBrighten || (blue[i] > newBlue)) blue[i]=newBlue;
				numSmallFixed++;
			}
			if (debugLevel>1) System.out.println("Corrected"+numSmallFixed+" blue pixels near overexposed blue");
		}
		if (dbgImgTitle!=null){
			String [] channel_blue_titles={"orig. blue","corrected blue","state"};
			double [][] dbgImg= new double [3][];
			dbgImg[0]=rgb_in[2];
			dbgImg[1]=blue;
			dbgImg[2]=dbgStateToDouble(state);
			SDFA_INSTANCE.showArrays(dbgImg, width, height, true, dbgImgTitle, channel_blue_titles);
		}
		
		return blue; // corrected blue channel
	}
	
	private double [] dbgStateToDouble(BLUE_STATE[] state){
		double [] result = new double[length];
		for (int i=0;i<length;i++){
			if      (state[i] == BLUE_STATE.OVEREXP)  result[i]=0.2;
			else if (state[i] == BLUE_STATE.EXPANDED) result[i]=0.6;
			else if (state[i] == BLUE_STATE.PROBLEM)  result[i]=0.4;
			else if (state[i] == BLUE_STATE.SOLUTION) result[i]=0.8;
			else if (state[i] == BLUE_STATE.SOLVED)   result[i]=1.0;
			else if (state[i] == BLUE_STATE.TMP)      result[i]=1.2;
		}
		return result;
		
	}
	
	private double [] getYrg (){
		double kr=colorProcParameters.kr;
		double kb=colorProcParameters.kb;
		double kg=1.0-kr-kb;
		double krg=kr+kg;
		kr/=krg;
		kg/=krg;
		double [] yrg=new double [rgb_in[0].length];
		for (int i=0;i<yrg.length;i++) {
			yrg[i]=kr*rgb_in[0][i]+kg*rgb_in[1][i];
		}
		return yrg;
	}


	private double [] overexposed(
			EyesisCorrectionParameters.ColorProcParameters  colorProcParameters,
			double [] dpixels,
			int width
			){
		int size=colorProcParameters.satDetSquareSize;
		double [] mask=overexposed(
				colorProcParameters,
				dpixels,
				width,
				0,
				0);
		for (int shftY=0; shftY<size; shftY+=size/2)
			for (int shftX=0; shftX<size; shftX+=size/2)
				if ((shftY!=0) || (shftX!=0)){
					double [] mask1=overexposed(
							colorProcParameters,
							dpixels,
							width,
							shftY,
							shftY);
					for (int i=0;i<mask.length;i++){
						if (Double.isNaN(mask[i])){
							mask[i]=mask1[i]; 
						} else if (!Double.isNaN(mask1[i])){
							if (Math.abs(mask[i]-dpixels[i])>Math.abs(mask1[i]-dpixels[i])){
								mask[i]=mask1[i];
							}
						}
					}
				}
		return mask;
	}


	private double [] overexposed(
			EyesisCorrectionParameters.ColorProcParameters  colorProcParameters,
			double [] dpixels,
			int width,
			int shftX,
			int shftY
			){
		int height=dpixels.length/width;
		int size=colorProcParameters.satDetSquareSize;
		double [] mask=new double[dpixels.length];
		for (int i=0;i<dpixels.length;i++) mask[i]=Double.NaN;
		double avrg=0.0;
		for (int i=0;i<dpixels.length;i++) avrg+=dpixels[i];
		avrg/=dpixels.length;
		double [] tile=new double [size*size];
		//		  int dbgCount=0;
		for (int y0=shftY; y0<(height-size/2); y0+=size){
			int ymax=size;
			if ((y0+ymax) > height) ymax=height-y0;

			for (int x0=shftX;x0<(width-size/2);x0+=size){
				int len=0;
				int indx;
				int xmax=size;
				if ((x0+xmax) > width) xmax=width-x0;
				for (int y=0;y<ymax;y++){
					indx=(y+y0)*width+x0;
					for (int x=0;x<xmax;x++){
						tile[len++]=dpixels[indx++];
					}
				}
				double tileAvg=0.0;
				for (int i=0;i<len;i++) tileAvg+=tile[i];
				tileAvg/=len;
				if (tileAvg<(avrg*colorProcParameters.satDetMinFrac)) continue; // too small value to be an overexposed tile
				int numFit=0;
				double wnd_min=tileAvg - avrg*colorProcParameters.satDetRelDiff;
				double wnd_max=tileAvg + avrg*+colorProcParameters.satDetRelDiff;

				for (int i=0;i<len;i++) if ((tile[i]>=wnd_min) && (tile[i]<=wnd_max)) numFit++;
				if (numFit< (len*colorProcParameters.satDetPartInside)) continue; // too small number of pixels fit close to average;
				len=0;
				double wnd_min_fin=tileAvg-avrg*colorProcParameters.satDetFinRelDiff;
				double wnd_max_fin=tileAvg+avrg*colorProcParameters.satDetFinRelDiff;

				for (int y=0;y<ymax;y++){
					indx=(y+y0)*width+x0;
					for (int x=0;x<xmax;x++){
						double d=tile[len++];
						double dpx=dpixels[indx];
						mask[indx++] =  ((d>=wnd_min_fin) && (d<=wnd_max_fin))?dpx:Double.NaN;
					}
				}
			}
		}
		return mask;
	}
	
	private void overexposedGrow(
			double [] dpixels,
			double [] ovrexp,
			BLUE_STATE[] state,
			List<Integer> pixelList,
			double relDiff,
			boolean addBrighter,
			double newWeight,
			int numSteps,
			boolean use8) {
		int [] dirs=use8?dirs8:dirs4;
		List<Integer> laterCheckList=new ArrayList<Integer>(); // put pixels that do not match current filter, but may be used later
// TODO: move markBlueOverexp and createInitialList outside as this method can be called several times
		if (debugLevel>1) System.out.println("overexposedGrow, relDiff="+relDiff+", addBrighter="+addBrighter);
		for (int step=0;(step<numSteps) && !pixelList.isEmpty();step++){
			// first pass - consume pixelList, re-add satisfying ones, mark them as TMP
			int len=pixelList.size();
			for (int i=0;i<len;i++){
				int oldIndex=pixelList.remove(0);
// first find if this pixel falls withing overexposed classification, then use weighted average of known overexposed pixels
// and pixel value (if it is inside range)
				for (int dir=0; dir<dirs.length; dir++){
					int index=oldIndex+dirs[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if (state[index]!=BLUE_STATE.OVEREXP) continue;
					double d=ovrexp[index]-dpixels[oldIndex];
					double abs_d=Math.abs(d);
					if ((abs_d>(relDiff*ovrexp[index])) && (!addBrighter || (d>0.0))) continue; // not eligible
					state[oldIndex]=BLUE_STATE.TMP;
					break;
				}
				// added, now calculate the overexposed value
				if (state[oldIndex]==BLUE_STATE.TMP){
					pixelList.add(new Integer(oldIndex)); // re-add to the list
					double w=0.0;
					double avgRef=0.0;
					for (int dir=0; dir<dirs20.length; dir++){
						int index=oldIndex+dirs20[dir];
						int y=index/width;
						if ((y<margin) || (y>=yMax)) continue;
						int x=index%width;
						if ((x<margin) || (x>=xMax)) continue;
						if (state[index]!=BLUE_STATE.UNDEF){ // can be problem/solution or overexp?
							avgRef+= weights20[dir]*ovrexp[index];
							w+=weights20[dir];
						}
					}
					assert (w>0):"marked with TMP are supposed to have at least 1 eligible neighbor";
					double avgOther=avgRef/w;
					if (Math.abs(dpixels[oldIndex]-avgOther)<=(relDiff*avgOther)){ // overshoots should not contribute to the average saturation level
						avgRef+=newWeight*dpixels[oldIndex];
						w+=newWeight;
					}
					ovrexp[oldIndex]= avgRef/=w;
				} else {
					laterCheckList.add(new Integer(oldIndex)); // check if it will be valid later when filter will be loosened
				}
			}
			// mark list (non-destructively) with OVEREXP (was TMP)
			for (Iterator<Integer> iter= pixelList.iterator(); iter.hasNext();){
				state[iter.next()]=BLUE_STATE.OVEREXP;
			}
			if (debugLevel>2) System.out.print("Step "+step+" pixelList.size()="+pixelList.size());
			// Consume list, mark UNDEFINED neighbors as TMP, add to list
			len=pixelList.size();
			for (int i=0;i<len;i++){
				int oldIndex=pixelList.remove(0);
				for (int dir=0; dir<dirs.length; dir++){
					int index=oldIndex+dirs[dir];
					int y=index/width;
					if ((y<margin) || (y>=yMax)) continue;
					int x=index%width;
					if ((x<margin) || (x>=xMax)) continue;
					if (state[index]!=BLUE_STATE.UNDEF) continue;
					pixelList.add(new Integer(index));
					state[index]=BLUE_STATE.TMP;
				}
			}
			
			
			// change TMP back to UNDEF (nondestructively)
			for (Iterator<Integer> iter= pixelList.iterator(); iter.hasNext();){
				state[iter.next()]=BLUE_STATE.UNDEF;
			}
			if (debugLevel>2) System.out.println(" --> "+pixelList.size());
		}
		// add abandoned pixels - they may become valid again with different filter
		// first mark current pixels
		if (debugLevel>1) System.out.print("pixelList.size()="+pixelList.size()+",  laterCheckList.size()="+laterCheckList.size());
		for (Iterator<Integer> iter= pixelList.iterator(); iter.hasNext();){
			state[iter.next()]=BLUE_STATE.TMP;
		}
		for (Iterator<Integer> iter= laterCheckList.iterator(); iter.hasNext();){
			Integer index=iter.next();
			if (state[index]==BLUE_STATE.UNDEF){
				state[index]=BLUE_STATE.TMP;
				pixelList.add(index);
			}
		}
		// change TMP back to UNDEF (nondestructively)
		for (Iterator<Integer> iter= pixelList.iterator(); iter.hasNext();){
			state[iter.next()]=BLUE_STATE.UNDEF;
		}
		if (debugLevel>1) System.out.println(" ==> combined pixelList.size()="+pixelList.size());
	}	
	
	public double [] detailsMask(
			double [] pixels,
			double sigma,
			double relThreshold){ // relative value, 
		double [] fpix=new double [length];
		for (int i=0;i<length;i++) fpix[i]=pixels[i];
		DoubleGaussianBlur gb = new DoubleGaussianBlur();
		gb.blurDouble(fpix,width,height,sigma,sigma, 0.01);
		for (int i=0;i<length;i++) fpix[i]-=(pixels[i]+relThreshold*fpix[i]);
		return fpix;
	}
	
	
	
	public double [][] process(){
		boolean debugThis=debugLevel>1;
		boolean calcRGOver=false; // currently red and green overexposure is not used, only for debug images - sabve on calculations

		double [][] colorDiff=null;
		double [][] dbg_proc=null;
//		double [][] overExpPix = new double [3][length]; // may be needed all (not yet)
		double [][] overExpPix = new double [3][]; // may be needed all (not yet)
		overExpPix[2]=new double[length];
		if (calcRGOver || debugThis){
			overExpPix[0]=new double[length];
			overExpPix[1]=new double[length];
		}
		if (debugThis) {
			colorDiff=new double[2][length];
			for (int px=0;px<length;px++){
				colorDiff[0][px]=rgb_in[0][px]/rgb_in[1][px];
				colorDiff[1][px]=rgb_in[2][px]/rgb_in[1][px];
			}
			dbg_proc =new double[11][]; // Only some are currently needed, rest is just for debugging
			dbg_proc[0]=rgb_in[0];
			dbg_proc[1]=rgb_in[1];
			dbg_proc[2]=rgb_in[2];
			dbg_proc[3]=colorDiff[0];
			dbg_proc[4]=colorDiff[1];
		}
		

		for (int chn=(calcRGOver || debugThis)?0:2;chn<3;chn++) {
			Runtime runtime = Runtime.getRuntime();
			runtime.gc();
			overExpPix[chn] =    overexposed(
					colorProcParameters,
					rgb_in[chn],
					width);
			if (debugThis) {
				dbg_proc[5+chn]=new double[length];
					for (int i=0;i<length;i++) dbg_proc[5+chn][i]= overExpPix[chn][i];
			}
			
			BLUE_STATE[] state= markBlueOverexp(overExpPix[chn]);
			List<Integer>	pixelList= createInitialList(
						state,
						BLUE_STATE.UNDEF,
						BLUE_STATE.OVEREXP,
						colorProcParameters.use8);
			double avrg=0.0;
			for (int i=0;i<length;i++) avrg+=rgb_in[chn][i];
			avrg/=length;
			if (debugLevel>1) System.out.println("*** chn="+chn+", avrg="+avrg);
			overexposedGrow(
					rgb_in[chn], //double [] dpixels,
					overExpPix[chn], // double [] ovrexp,
					state, // BLUE_STATE[] state,
					pixelList,
					colorProcParameters.satDetGrowRelDiff, // double maxDiff,
					false,  // boolean addBrighter,
					colorProcParameters.satDetNewWeight,
					colorProcParameters.satDetExpSym, // int numSteps,
					colorProcParameters.use8); // boolean use8) {
	
			overexposedGrow(
					rgb_in[chn], //double [] dpixels,
					overExpPix[chn], // double [] ovrexp,
					state, // BLUE_STATE[] state,
					pixelList,
					colorProcParameters.satDetGrowRelDiff, // double maxDiff,
					true,  // boolean addBrighter,
					colorProcParameters.satDetNewWeight,
					colorProcParameters.satDetExpOver, // int numSteps,
					colorProcParameters.use8); // boolean use8) {
			
			overexposedGrow(
					rgb_in[chn], //double [] dpixels,
					overExpPix[chn], // double [] ovrexp,
					state, // BLUE_STATE[] state,
					pixelList,
					colorProcParameters.satDetGrowRelDiffCleanUp, // double maxDiff,
					true,  // boolean addBrighter,
					colorProcParameters.satDetNewWeight,
					colorProcParameters.satDetExpCleanUp, // int numSteps,
					colorProcParameters.use8); // boolean use8) {
		}
		
		if (colorProcParameters.blueLeakFixWires){
			double [] details=	detailsMask(
					yrg,
					colorProcParameters.blueLeakWiresSize, // double sigma,
					colorProcParameters.blueLeakWiresThreshold); //double relThreshold, // relative value,
			int numDet=0;
			for (int i=0;i<length;i++)
				if (!(overExpPix[2][i]==Double.NaN) && (details[i]>0.0)){
				overExpPix[2][i]=Double.NaN; // mark as if not over-exposed
				numDet++;
			}
			if (debugLevel>1) System.out.println("Marked "+numDet+" blue pixels as if not overexposed (trees, wires)");
			if (debugLevel>1) {
				SDFA_INSTANCE.showArrays(details, width, height,"Details");
			}
		}

		
		
		if (debugThis){
			String [] channel_titles={"red","green","blue", "r/g", "b/g","oe_R","oeG","oe_B","filt_oe_R","filt_oeG","filt_oe_B"};
			dbg_proc[8]=overExpPix[0];
			dbg_proc[9]=overExpPix[1];
			dbg_proc[10]=overExpPix[2];

			SDFA_INSTANCE.showArrays(dbg_proc, width, height, true, "test_chn",channel_titles);
		}
		/// public void showArraysSparse(double[][] pixels, int width, int height,  boolean asStack, String title, String [] titles) {
		double [] blue = findBlueSolutions(
				overExpPix[2],
				colorProcParameters.blueLeakNoHint,
				colorProcParameters.blueLeakNoBrighten,
				colorProcParameters.use8); // use8);
		double [][] rgb_out={rgb_in[0],rgb_in[1],blue};
		return rgb_out;
		//TODO: add processing of R/G/b overexposed areas using the least sensitive channel. 
		
	}
}
