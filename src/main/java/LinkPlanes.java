/**
 **
 ** LinkPlanes - manage links between supertile planes
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  LinkPlanes.java is free software: you can redistribute it and/or modify
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
import java.awt.Point;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;

public class LinkPlanes {
	SuperTiles st;
	public boolean    plPreferDisparity; //    =   false;// Always start with disparity-most axis (false - lowest eigenvalue)
	public double     plDispNorm; //           =   3.0;  // Normalize disparities to the average if above (now only for eigenvalue comparison)
	public double     plMinStrength; //        =   0.1;  // Minimal total strength of a plane 
	public double     plMaxEigen; //           =   0.05; // Maximal eigenvalue of a plane 
	public double     plEigenFloor; //         =   0.01; // Add to eigenvalues of each participating plane and result to validate connections 
	public double     plEigenStick; //         =   25.0; // Consider plane to be a "stick" if second eigenvalue is below 
	public double     plBadPlate; //           =   0.2;  // Not a plate if sin^2 between normals from disparity and world exceeds this 

	public double     plWorstWorsening; //     =   2.0;  // Worst case worsening after merge
	public double     plWorstWorsening2; //    =   5.0;  // Worst case worsening for thin planes
	public double     plWorstEq; //            =   1.0;  // Worst case worsening after merge with equal weights
	public double     plWorstEq2; //           =   2.0;  // Worst case worsening for thin planes with equal weights

	public double     plOKMergeEigen; //       =   0.03; // If result of the merged planes is below, OK to use thin planes (higher) threshold 
	public double     plMaxWorldSin2; //       =   0.1;  // Maximal sine squared of the world angle between planes to merge. Set to >= 1.0 to disable
    public double     pl2dForSin; //           =   2.5;  // Do not compare sin() between planes, if at least one has too small axis ratio 

	public double     plWeakWorsening; //      =   1.0;  // Relax merge requirements for weaker planes
	public double     plMaxOverlap; //         =   0.1;  // Maximal overlap between the same supertile planes to merge
		// Merge same supetile planes if at least one is weak and they do not differ much
	public double     plWeakWeight; //         =   0.2; // Maximal weight of the weak plane to merge
	public double     plWeakEigen; //          =   0.1; // Maximal eigenvalue of the result of non-weighted merge
    public double     plWeakWeight2; //        =  10.0 ; // Maximal weight of the weak plane to merge (second variant)
	public double     plWeakEigen2; //         =   0.05; // Maximal eigenvalue of the result of non-weighted merge  (second variant)
    public double     plSumThick; //           =   1.2;  // Do not merge if any sqrt of merged eigenvalue exceeds scaled sum of components
    public double     plNeNeibCost; //         =   5.0;  // When calculating non-exclusive planes, do not use neighbors with high cost
    public double     plNeOwn; //              =   5.0;  // When calculating non-exclusive planes, use center plane relative weight
    public double     plExNeibCost; //         =   5.0;  // When calculating non-exclusive planes, do not use neighbors with high cost
    public double     plExNeibSmooth; //       =   0.5;  // Scale down maximal costs for smoothed planes (tighter requirements) 
    public double     plMergeCostStar; //      =   5.0;  // Cost threshold for merging same tile planes if the plane has connected neighbors
    public double     plMergeCost; //          =  10.0;  // Cost threshold for merging same tile planes if not connected
    
    public boolean    plConflMerge; //         =  true;  // Try to merge conflicting planes
    public double     plConflRelax; //         =   1.5;  // Scale parameters to relax planes fit for merging conflicting planes
    public boolean    plConflSngl; //          =  true;  // Only merge conflicting planes if this is the only conflicting pair in the supertile
    public boolean    plConflSnglPair; //      =  true;  // Only merge conflicting planes only if there are only two planes in the supertile 

    public double     plWeakFgStrength; //     =   0.15; // Consider merging plane if it is foreground and maximal strength below this  
    public int        plWeakFgOutliers; //     =   1;    // Remove these strongest from foreground when determining the maximal strength
    public double     plWeakFgRelax; //        =   2.0;  // Relax cost requirements when merging with weak foreground   
    
    public double     plThickWorld; //         =   0.2;  // Maximal real-world thickness of merged overlapping planes (meters) 
    public double     plThickWorldConfl; //    =   0.4;  // Maximal real-world merged thickness for conflicting planes 
    public double     plRelaxComplete;  //     =   1.5;  // Relax cost requirements when adding exclusive links to complete squares and triangles 
    public double     plRelaxComplete2; //     =   2.0;  // Relax cost requirements during the second pass
    
    public boolean    plFillSquares; //        =   true; // Add diagonals to full squares
    public boolean    plCutCorners; //         =   true; // Add ortho to 45-degree corners
    public boolean    plHypotenuse; //         =   true; // Add hypotenuse connection if both legs exist
    

		// comparing merge quality for plane pairs
	public double     plCostDist; //           =   4.0;  // Disparity (pix) - closer cost will use more of the real world, farther - disparity
	public double     plCostKrq; //            =   0.8;  // cost of merge quality weighted in disparity space
	public double     plCostKrqEq; //          =   0.2;  // cost of merge quality equal weight in disparity space
	public double     plCostWrq; //            =   0.8;  // cost of merge quality weighted in world space
	public double     plCostWrqEq; //          =   0.2;  // cost of merge quality equal weight  in world space
	public double     plCostSin2; //           =  10.0;  // cost of sin squared between normals
	public double     plCostRdist2; //         =1000.0;  // cost of squared relative distances
	
	public double     plMaxZRatio; //          =   2.5;  // Maximal ratio of Z to allow plane merging
	public double     plMaxDisp; //            =   0.5;  // Maximal disparity of one of the planes to apply  maximal ratio
	public double     plCutTail; //            =   1.4;  // When merging with neighbors cut the tail that is worse than scaled best
	public double     plMinTail; //            =   0.015;// Set cutoff value level not less than

	public int        plMinPoints; //          =     5;  // Minimal number of points for plane detection
	public double     plTargetEigen; //        =   0.02; // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
	public double     plFractOutliers; //      =   0.3;  // Maximal fraction of outliers to remove
	public int        plMaxOutliers; //        =    20;  // Maximal number of outliers to remove

	public double     plPull               =  5.0; // .3;   // Relative weight of original (measured) plane compared to average neighbor pull
	public int        plIterations         =  10;   // Maximal number of smoothing iterations for each step
	public boolean    plStopBad            =  true; // Do not update supertile if any of connected neighbors is not good (false: just skip that neighbor)
	public int        plPrecision          =  6;    // Maximal step difference (1/power of 10)
	public double     plNormPow            =  .5;   // 0.0: 8 neighbors pull 8 times as 1, 1.0 - same as 1
	
	
	
	public int        dbg_tileX;
	public int        dbg_tileY;


	public LinkPlanes (
			EyesisCorrectionParameters.CLTParameters           clt_parameters,
			SuperTiles st)
	{
		plWorstWorsening =  clt_parameters.plWorstWorsening;
		plWorstWorsening2 = clt_parameters.plWorstWorsening2;
		plWorstEq =         clt_parameters.plWorstEq;
		plWorstEq2 =        clt_parameters.plWorstEq2;

		plWeakWorsening =   clt_parameters.plWeakWorsening;
		plOKMergeEigen =    clt_parameters.plOKMergeEigen;
		plMaxWorldSin2 =    clt_parameters.plMaxWorldSin2;
  		pl2dForSin =        clt_parameters.pl2dForSin;

		plDispNorm =        clt_parameters.plDispNorm;
		plMaxEigen =        clt_parameters.plMaxEigen;
		plEigenFloor =      clt_parameters.plEigenFloor; 
		plEigenStick =      clt_parameters.plEigenStick; 
		plBadPlate =        clt_parameters.plBadPlate;
		plMinStrength =     clt_parameters.plMinStrength;
		plMaxOverlap =      clt_parameters.plMaxOverlap;
		
		plWeakWeight =      clt_parameters.plWeakWeight;
		plWeakEigen =       clt_parameters.plWeakEigen;
		plWeakWeight2 =     clt_parameters.plWeakWeight2;
		plWeakEigen2 =      clt_parameters.plWeakEigen2;
		plSumThick =        clt_parameters.plSumThick;
		plNeNeibCost =      clt_parameters.plNeNeibCost;
		plNeOwn =           clt_parameters.plNeOwn;
		plExNeibCost =      clt_parameters.plExNeibCost;
		plExNeibSmooth =    clt_parameters.plExNeibSmooth;
		plMergeCostStar =   clt_parameters.plMergeCostStar;
		plMergeCost =       clt_parameters.plMergeCost;
		
		plConflMerge =      clt_parameters.plConflMerge;
		plConflRelax =      clt_parameters.plConflRelax;
		plConflSngl =       clt_parameters.plConflSngl;
		plConflSnglPair =   clt_parameters.plConflSnglPair;

		plWeakFgStrength =  clt_parameters.plWeakFgStrength;
		plWeakFgOutliers =  clt_parameters.plWeakFgOutliers;
		plWeakFgRelax =     clt_parameters.plWeakFgRelax;

		plThickWorld =      clt_parameters.plThickWorld;
		plThickWorldConfl = clt_parameters.plThickWorldConfl;
		plRelaxComplete =   clt_parameters.plRelaxComplete;
		plRelaxComplete2 =  clt_parameters.plRelaxComplete2;

		plFillSquares =     clt_parameters.plFillSquares;
		plCutCorners =      clt_parameters.plCutCorners;
		plHypotenuse =      clt_parameters.plHypotenuse;

		plCostDist =        clt_parameters.plCostDist;
		plCostKrq =         clt_parameters.plCostKrq;
		plCostKrqEq =       clt_parameters.plCostKrqEq;
		plCostWrq =         clt_parameters.plCostWrq;
		plCostWrqEq =       clt_parameters.plCostWrqEq;
		plCostSin2 =        clt_parameters.plCostSin2;
		plCostRdist2 =      clt_parameters.plCostRdist2;       

  		plMaxZRatio =       clt_parameters.plMaxZRatio;
  		plMaxDisp =         clt_parameters.plMaxDisp;
  		plCutTail =         clt_parameters.plCutTail;
  		plMinTail =         clt_parameters.plMinTail;
  		
  		plMinPoints =       clt_parameters.plMinPoints;
  		plTargetEigen =     clt_parameters.plTargetEigen;
  		plFractOutliers =   clt_parameters.plFractOutliers;
  		plMaxOutliers =     clt_parameters.plMaxOutliers;
  		
  		plPull =            clt_parameters.plPull;
  		plIterations =      clt_parameters.plIterations;
  		plStopBad =         clt_parameters.plStopBad;
  		plPrecision =       clt_parameters.plPrecision;
  		plNormPow =         clt_parameters.plNormPow;
		
		dbg_tileX =         clt_parameters.tileX;
		dbg_tileY =         clt_parameters.tileY;
		this.st =           st;
	}
	
	public double getExNeibCost(){
		return plExNeibCost;
	}
	public double getMergeCostStar(){
		return plMergeCostStar;
	}
	public double getMergeCostNoStar(){
		return plMergeCost;
	}
	public double getExNeibSmooth(){
		return plExNeibSmooth;
	}

	public double getConflRelax(){
		return plConflRelax;
	}
	
	public boolean getConflMerge(){
		return plConflMerge;
	}

	public double getRelaxComplete(){
		return plRelaxComplete;
	}
	public double getRelaxComplete2(){
		return plRelaxComplete2;
	}

	public boolean areWeakSimilar(
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			String prefix,
			int debugLevel)
	{
		return planesFit(
				plane1,        // should belong to the same supertile (or be converted for one)
				plane2,
				1.0,           // double               relax,
				true,          // use for same supertile merge
				true,          // boolean              check_is_weak_only, // only verify if the two planes are close and one is weak 
				Double.NaN,    // if NaN will calculate assuming the same supertile
				merged_ev_eq,  // if NaN will calculate assuming the same supertile
				Double.NaN,    // if NaN will calculate assuming the same supertile - for world
				Double.NaN,    // if NaN will calculate assuming the same supertile - for world
				prefix,
				debugLevel);
		
	}
	
	public boolean planesFit(
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               relax,
			boolean              merge_weak,    // use for same supertile merge
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
			double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
			String prefix,
			int debugLevel)
	{
		return planesFit(
				plane1, // should belong to the same supertile (or be converted for one)
				plane2,
				relax,         // double               relax,
				merge_weak,    // use for same supertile merge
				false, // boolean              check_is_weak_only, // only verify if the two planes are close and one is weak 
				merged_ev,    // if NaN will calculate assuming the same supertile
				merged_ev_eq, // if NaN will calculate assuming the same supertile
				merged_wev,    // if NaN will calculate assuming the same supertile - for world
				merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
				prefix,
				debugLevel);
	}
	
	public boolean planesFit(
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               relax,
			boolean              merge_weak,    // use for same supertile merge
			boolean              check_is_weak_only, // only verify if the two planes are close and one is weak 
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
			double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
			String prefix,
			int debugLevel)
	{
		double plSumThick =   this.plSumThick * relax;
		double plWeakEigen =  this.plWeakEigen * relax;
		double plWeakEigen2 = this.plWeakEigen2 * relax;
		merge_weak |= check_is_weak_only;
		if ((plane1 == null) || (plane2 == null)) return false;
		boolean dbg = debugLevel > 1;
		if (debugLevel > 1){
			System.out.println("planesFit() debug:");
		}
		TilePlanes.PlaneData merged_pd = null;
		TilePlanes.PlaneData merged_pd_eq = null;
		
		double disp1 = plane1.getZxy()[0];
		double disp2 = plane2.getZxy()[0];
		if ((disp1 <= 0) || (disp2 <= 0) || (((disp1 <= plMaxDisp) || (disp2 <= plMaxDisp)) && ((disp1/disp2 > plMaxZRatio)  || (disp2/disp1 > plMaxZRatio)))){
			if (debugLevel > 0)	System.out.println(prefix+" planes have too high disparity ratio ("+disp1+":"+disp2+" > plMaxZRatio="+plMaxZRatio+
					", at least one of them < plMaxDisp="+plMaxDisp);
			return false;
		} else {
			if (debugLevel > 1)	System.out.println(prefix+" disparity ratio ("+disp1+":"+disp2+" is OK, <= plMaxZRatio="+plMaxZRatio);
		}
		
		if (!merge_weak && (plane1.getWeight() < plMinStrength)) {
			if (debugLevel > 1)	System.out.println(prefix+" plane1 is too weak ("+plane1.getWeight()+" < plMinStrength="+plMinStrength+")");
			return false;
		}
		if (!merge_weak && (plane2.getWeight() < plMinStrength)) {
			if (debugLevel > 1)	System.out.println(prefix+" plane2 is too weak ("+plane2.getWeight()+" < plMinStrength="+plMinStrength+")");
			return false;
		}
		double corr_max_eigen = corrMaxEigen(
				plMaxEigen,
				plDispNorm,
				plane1);
		if ((plMaxEigen != 0.0) &&
				(plane1.getValue() > corr_max_eigen)){
			if (debugLevel > 1)	System.out.println(prefix+" plane1 is too thick ("+plane1.getValue()+" > corr_max_eigen="+corr_max_eigen+")");
			return false;
		}

		if ((plMaxEigen != 0.0) &&
				(plane2.getValue() > corr_max_eigen)){
			if (debugLevel > 1)	System.out.println(prefix+" plane2 is too thick ("+plane2.getValue()+" > corr_max_eigen="+corr_max_eigen+")");
			return false;
		}
		
		if (!check_is_weak_only && (Double.isNaN(merged_ev) || Double.isNaN(merged_wev))) {
			merged_pd = plane1.mergePlaneToThis(
					plane2, // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					false,       // boolean   ignore_weights,
					true, // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 2); // int       debugLevel)
			merged_ev = merged_pd.getValue();
			merged_wev = merged_pd.getWValue();
		}
		if (Double.isNaN(merged_ev_eq) || (!check_is_weak_only && Double.isNaN(merged_wev_eq))) {
			merged_pd_eq = plane1.mergePlaneToThis(
					plane2,      // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					true,        // boolean   ignore_weights,
					true,        // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 2); // int       debugLevel)
			merged_ev_eq = merged_pd_eq.getValue();
			merged_wev_eq = merged_pd_eq.getWValue();
			
		}
		double w1 = plane1.getWeight();
		double w2 = plane2.getWeight();
		double weakest = (w1 > w2) ? w2 : w1;
		double this_rq_eq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue() , // double L2,
				merged_ev_eq, // double L,
				1.0, // double w1,
				1.0, // double w2)
				plEigenFloor);//			double eigen_floor)
//		if (this_rq_eq == 0.01) System.out.println("planesFit(): this_rq_eq was negative");
		double this_rq_eq_norm = this_rq_eq;
		if ((w1 + w2) < plWeakWorsening) this_rq_eq_norm *= (w1 + w2) / plWeakWorsening; // forgive more for weak planes 
		boolean weak_and_close = false;
		
		if (merge_weak && (weakest <= plWeakWeight) && (merged_ev_eq <= plWeakEigen )){
			weak_and_close = true;
			if (dbg) System.out.println(prefix+": same supertile planes are weak and close: the weakest ("+weakest+") is below "+plWeakWeight+
					" and merged non-weighted eigenvalue ("+merged_ev_eq+") is below "+plWeakEigen);
		}
		if (merge_weak && (weakest <= plWeakWeight2) && (merged_ev_eq <= plWeakEigen2 )){
			weak_and_close = true;
			if (dbg) System.out.println(prefix+": same supertile planes are weak and close (variant 2): the weakest ("+weakest+") is below "+plWeakWeight2+
					" and merged non-weighted eigenvalue ("+merged_ev_eq+") is below "+plWeakEigen2);
		}
//		if (check_is_weak_only) return weak_and_close;
		
		
		double this_rq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue() , // double L2,
				merged_ev, // double L,
				w1, // double w1,
				w2, // double w2)
				plEigenFloor);//			double eigen_floor)
//		if (this_rq == 0.01) System.out.println("planesFit(): this_rq was negative");

		double this_rq_norm = this_rq;
		if ((w1 + w2) < plWeakWorsening) this_rq_norm *= (w1 + w2) / plWeakWorsening; // forgive more for weak planes

		
		double this_wrq = mergeRQuality(
				plane1.getWValue(), // double L1,
				plane2.getWValue(), // double L2,
				merged_wev, // double L,
				w1, // double w1,
				w2, // double w2)
				0.0);//			double eigen_floor)
//		if (this_wrq == 0.01) System.out.println("planesFit(): this_wrq was negative");

		double this_wrq_norm = this_wrq;
		if ((w1 + w2) < plWeakWorsening) this_wrq_norm *= (w1 + w2) / plWeakWorsening; // forgive more for weak planes

		double this_wrq_eq = mergeRQuality(
				plane1.getWValue(), // double L1,
				plane2.getWValue(), // double L2,
				merged_wev_eq, // double L,
				1.0, // double w1,
				1.0, // double w2)
				0.0);//			double eigen_floor)
///		this_wrq_eq /= (w1 + w2); // for comparison reduce this value for stronger planes
//		if (this_wrq_eq == 0.01) {
//			System.out.println("planesFit(): this_wrq_eq was negative");
//		}
		double this_wrq_eq_norm = this_wrq_eq;
		if ((w1 + w2) < plWeakWorsening) this_wrq_eq_norm *= (w1 + w2) / plWeakWorsening; // forgive more for weak planes 
		
		boolean OK_to_merge = false;
		boolean notOK_to_merge = false;
		// plSumThick should override weak_and_close
		double sum_thick = Math.sqrt(plane1.getValue()) + Math.sqrt(plane2.getValue());
		double sum_wthick = Math.sqrt(plane1.getWValue()) + Math.sqrt(plane2.getWValue());
		if (Math.sqrt(merged_ev) > plSumThick * sum_thick){
			notOK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes do not fit as weighted pixel thickness = "+Math.sqrt(merged_ev)+
					" > " +plSumThick+" * "+ sum_thick+" = "+(plSumThick * sum_thick) );
		}
		if (Math.sqrt(merged_ev_eq) > plSumThick * sum_thick){
			notOK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes do not fit as equalized pixel thickness = "+Math.sqrt(merged_ev_eq)+
					" > " +plSumThick+" * "+ sum_thick+" = "+(plSumThick * sum_thick) );
		}
		if (Math.sqrt(merged_wev) > plSumThick * sum_wthick){
			notOK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes do not fit as weighted world thickness = "+Math.sqrt(merged_wev)+
					" > " +plSumThick+" * "+ sum_wthick+" = "+(plSumThick * sum_wthick) );
		}
		if (Math.sqrt(merged_wev_eq) > plSumThick*sum_wthick){
			notOK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes do not fit as equalized world thickness = "+Math.sqrt(merged_wev_eq)+
					" > " +plSumThick+" * "+ sum_wthick+" = "+(plSumThick * sum_wthick) );
		}
		
		if (notOK_to_merge) {
			weak_and_close = false;
		}
		
		
		
		if (!notOK_to_merge && (this_rq_norm <= plWorstWorsening)){
			OK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes may fit  (this_rq_norm="+this_rq_norm+"  <= plWorstWorsening="+plWorstWorsening+")");
		}
		if (!notOK_to_merge && (!OK_to_merge || dbg) && (merged_ev <= plOKMergeEigen) && (this_rq_norm <= plWorstWorsening2)){
			OK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes may fit  (merged_ev="+merged_ev+
					" <= plOKMergeEigen="+plOKMergeEigen+") && (this_rq_norm="+this_rq_norm+
					" <= plWorstWorsening2="+plWorstWorsening2+")");
		}
		if (!notOK_to_merge && (!OK_to_merge || dbg) && (this_rq_eq_norm <= plWorstEq)){
			OK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes may fit (this_rq_eq_norm="+this_rq_eq_norm+" <= plWorstEq="+plWorstEq+")");
		}
		if (!notOK_to_merge && (!OK_to_merge || dbg) && (merged_ev_eq <= plOKMergeEigen) && (this_rq_eq_norm <= plWorstEq2)){
			OK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes may fit (merged_ev_eq="+merged_ev_eq+" <= plOKMergeEigen) && (this_rq_eq_norm="+
								this_rq_eq_norm+" <= plWorstEq2="+plWorstEq2+")");
		}
		
		if (!notOK_to_merge && (!OK_to_merge || dbg) && (this_wrq_norm <= plWorstWorsening)){
			OK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes may fit  (this_wrq_norm="+this_wrq_norm+"  <= plWorstWorsening="+plWorstWorsening+")");
		}
		if (!notOK_to_merge && (!OK_to_merge || dbg) && (this_wrq_eq_norm <= plWorstEq)){
			OK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes may fit (this_wrq_eq_norm="+this_wrq_eq_norm+" <= plWorstEq="+plWorstEq+")");
		}
		// do not apply sin2 to weak planes or if they are not flat enough
		// TODO: Compare what can be compared (for 2 linear features or plane + linear?
		
		double min2d = Math.min(plane1.get2dRatio(), plane2.get2dRatio());
		if (Double.isNaN(min2d)) min2d = pl2dForSin;
		if ((plMaxWorldSin2 < 1.0) &&
				(min2d >= pl2dForSin) &&
				(!merge_weak || (weakest > plWeakWeight)) &&
				(plane1.checkBadPlate(true) < plBadPlate) &&
				(plane2.checkBadPlate(true) < plBadPlate) &&
				(plane1.getWorldSin2(plane2) > plMaxWorldSin2)) {
			notOK_to_merge = true;
			if (dbg) System.out.println(prefix+": planes do not fit as sin2 > "+plMaxWorldSin2+" weakest="+weakest+
					" or minimal 2d = "+min2d+" >= pl2dForSin=" + pl2dForSin +" and none of the planes is \"bad\"");
		}
		if (weak_and_close){
			OK_to_merge = true;
			notOK_to_merge = false; // weak can have large angles
			if (dbg) System.out.println(prefix+": same supertile planes fit as the weakest ("+weakest+") is below "+plWeakWeight+
					" and merged non-weighted eigenvalue ("+merged_ev_eq+") is below "+plWeakEigen);
		}
		
		
		if (check_is_weak_only) {
			if (debugLevel>1){
				System.out.println(prefix+" weak_and_close = " + weak_and_close+
						" this_rq="+this_rq+
						" this_rq_eq="+this_rq_eq+
						" this_wrq=" + (this_wrq) +
						" this_wrq_eq=" + (this_wrq_eq) +
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+
						" L="+merged_ev+" L_eq="+merged_ev_eq+
						" L1W="+plane1.getWValue()+" L2W="+plane2.getWValue()+" LW="+merged_wev+
						" L_eqW="+merged_wev_eq);
				System.out.println(prefix+" (fit) world sin2 ="+
						plane1.getWorldSin2(plane2));
				System.out.println(prefix+ " (fit)" +
						" world rdist this="+ Math.sqrt(plane1.getWorldPlaneRDist2(plane2))+
						", world rdist other="+Math.sqrt(plane2.getWorldPlaneRDist2(plane1))+
						", world rdist sum="+Math.sqrt(plane1.getWorldPlaneRDist2(plane2)+
								plane2.getWorldPlaneRDist2(plane1)));
				
			}
			return weak_and_close;
		}

		
		
		if (OK_to_merge && !notOK_to_merge) {
			if (debugLevel > 0){
				System.out.println(prefix+": planes FIT");
				if (debugLevel > 1){
					System.out.println(prefix+" (fit) this_rq="+this_rq+
							" this_rq_eq="+this_rq_eq+
							" this_wrq=" + (this_wrq) +
							" this_wrq_eq=" + (this_wrq_eq) +
							" w1="+w1+" w2="+w2+
							" L1="+plane1.getValue()+" L2="+plane2.getValue()+
							" L="+merged_ev+" L_eq="+merged_ev_eq+
							" L1W="+plane1.getWValue()+" L2W="+plane2.getWValue()+" LW="+merged_wev+
							" L_eqW="+merged_wev_eq+
							" relax=" +relax);
					System.out.println(prefix+" (fit) world sin2 ="+
							plane1.getWorldSin2(plane2));
					System.out.println(prefix+ " (fit)" +
							" world rdist this="+ Math.sqrt(plane1.getWorldPlaneRDist2(plane2))+
							", world rdist other="+Math.sqrt(plane2.getWorldPlaneRDist2(plane1))+
							", world rdist sum="+Math.sqrt(plane1.getWorldPlaneRDist2(plane2)+
									plane2.getWorldPlaneRDist2(plane1)));
				}
			}
			return true;
		}
		if (debugLevel > 0) {
			System.out.println(prefix+": planes DO NOT FIT");
			if (debugLevel > 1){
				System.out.println(prefix+" (do not fit) this_rq="+this_rq+
						", this_rq_eq="+this_rq_eq+
						" this_wrq=" + (this_wrq) +
						" this_wrq_eq=" + (this_wrq_eq) +
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+
						" L="+merged_ev+" L_eq="+merged_ev_eq+
						" L1W="+plane1.getWValue()+" L2W="+plane2.getWValue()+" LW="+merged_wev+
						" L_eqW="+merged_wev_eq +
						" relax=" +relax);
				System.out.println(prefix+" (do not fit) world sin2 ="+
						plane1.getWorldSin2(plane2));
				System.out.println(prefix+" (do not fit)"+
						", world dist this="+ Math.sqrt(plane1.getWorldPlaneRDist2(plane2))+
						", world dist other="+Math.sqrt(plane2.getWorldPlaneRDist2(plane1))+
						", world dist sum="+Math.sqrt(plane1.getWorldPlaneRDist2(plane2)+
								plane2.getWorldPlaneRDist2(plane1)));
			}
		}
		return false;
	}
	
	// 0 - this_rq,
	// 1 - this_rq_eq
	// 2 - composite
	// 3..7 contribution of each of the factor to the overall cost
	public double getLinkCost(
			boolean              en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
			double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
			String prefix,
			int debugLevel){
		int cost_index = 2;
		
		double [] costs = getFitQualities(
				en_sticks, // boolean              en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
				plane1, // should belong to the same supertile (or be converted for one)
				plane2,
				merged_ev,    // if NaN will calculate assuming the same supertile
				merged_ev_eq, // if NaN will calculate assuming the same supertile
				merged_wev,    // if NaN will calculate assuming the same supertile - for world
				merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
				prefix,
				debugLevel); // [cost_index];
		if (costs == null) return Double.NaN;
		else return costs[cost_index];
	}
	public double [] getFitQualities(
			boolean              en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
			double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
			String prefix,
			int debugLevel)
	{
		if ((plane1 == null) || (plane2 == null)) return null;
		boolean plate1 = true;
		boolean plate2 = true;
		if (en_sticks) {
			plate1 &= (plane1.getValues()[1] >=     plEigenStick);
			plate2 &= (plane1.getValues()[1] >=     plEigenStick);
			plate1 &= (plane1.checkBadPlate(true) < plBadPlate);
			plate2 &= (plane2.checkBadPlate(true) < plBadPlate);
		}
		if ((debugLevel > 0) && (!plate1 || !plate2)) {
			System.out.println(prefix+ " plate1="+plate1+", plate1="+plate2+" ("+(plane1.getValues()[1])+", "+plane2.getValues()[1]+" < "+plEigenStick+")"+
		   " or bad_plate1="+plane1.checkBadPlate(true)+" bad_plate2="+plane2.checkBadPlate(true)+" >= plBadPlate="+plBadPlate);
		}
		TilePlanes.PlaneData merged_pd = null;
		TilePlanes.PlaneData merged_pd_eq = null;
		if (Double.isNaN(merged_ev) || Double.isNaN(merged_wev)) {
			merged_pd = plane1.mergePlaneToThis(
					plane2, // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					false,       // boolean   ignore_weights,
					true, // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 2); // int       debugLevel)
			merged_ev =  merged_pd.getValue();
			merged_wev = merged_pd.getWValue();
		}
		if (Double.isNaN(merged_ev_eq) || Double.isNaN(merged_wev_eq)) {
			merged_pd_eq = plane1.mergePlaneToThis(
					plane2,      // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					true,        // boolean   ignore_weights,
					true,        // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 2); // int       debugLevel)
			merged_ev_eq =  merged_pd_eq.getValue();
			merged_wev_eq = merged_pd_eq.getWValue();
		}
		double w1 = plane1.getWeight();
		double w2 = plane2.getWeight();
		double this_rq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue(), // double L2,
				merged_ev, // double L,
				w1, // double w1,
				w2, // double w2)
				plEigenFloor);//			double eigen_floor)
//		if (this_rq == 0.01) System.out.println("getFitQualities(): this_rq was negative");
		double this_rq_nofloor = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue(), // double L2,
				merged_ev, // double L,
				w1, // double w1,
				w2, // double w2)
				0); // eigenFloor);//			double eigen_floor)
//		if (this_rq_nofloor == 0.01) System.out.println("getFitQualities(): this_rq_nofloor was negative");
		double this_rq_eq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue(), // double L2,
				merged_ev_eq, // double L,
				1.0, // double w1,
				1.0, // double w2)
				plEigenFloor);//			double eigen_floor)
///		this_rq /= (w1 + w2); // for comparison reduce this value for stronger planes
//		if (this_rq_eq == 0.01) System.out.println("getFitQualities(): this_rq_eq was negative");
		double this_wrq = mergeRQuality(
				plane1.getWValue(), // double L1,
				plane2.getWValue(), // double L2,
				merged_wev, // double L,
				w1, // double w1,
				w2, // double w2)
				0.0);//			double eigen_floor)
//		if (this_wrq == 0.01) System.out.println("getFitQualities(): this_wrq was negative");
///		this_wrq /= (w1 + w2); // for comparison reduce this value for stronger planes
		double this_wrq_eq = mergeRQuality(
				plane1.getWValue(), // double L1,
				plane2.getWValue(), // double L2,
				merged_wev_eq, // double L,
				1.0, // double w1,
				1.0, // double w2)
				0.0);//			double eigen_floor)
///		this_wrq_eq /= (w1 + w2); // for comparison reduce this value for stronger planes
//		if (this_wrq_eq == 0.01) System.out.println("getFitQualities(): this_wrq_eq was negative");
		double sin2 = plane1.getWorldSin2(plane2);
		double rdist2 = plane1.getWorldPlaneRDist2(plane2) + plane2.getWorldPlaneRDist2(plane1);
		if (!plate2 && plate1){
			rdist2 = 2 * plane2.getWorldPlaneRDist2(plane1);
		} else if (plate2 && plate1){
			rdist2 = 2 * plane1.getWorldPlaneRDist2(plane2);
		}
		// reduce cost contribution for disparity space parameters for near objects, increase - for far 
		
		double average_disp = 0.5 * (plane1.getZxy()[0] + plane2.getZxy()[0]);
		double cost_near = (this.plCostDist > 0.0) ?
				(2.0 * average_disp * average_disp / (average_disp * average_disp + this.plCostDist * this.plCostDist)) :
					1.0; 
		double plCostKrq =   this.plCostKrq    * (2.0 - cost_near);
		double plCostKrqEq = this.plCostKrqEq  * (2.0 - cost_near);
		double plCostWrq =   this.plCostWrq    *  cost_near;
		double plCostWrqEq = this.plCostWrqEq  *  cost_near;
		
		
		double [] costs = {
				Math.sqrt(this_rq  * this_rq_eq) *  plCostKrq,
				0.5 *    (this_rq  + this_rq_eq) *  plCostKrqEq,
				Math.sqrt(this_wrq * this_wrq_eq) * plCostWrq,
				0.5 *    (this_wrq + this_wrq_eq) * plCostWrqEq,
				sin2 *        plCostSin2,
				rdist2 *      plCostRdist2};
		double cost = costs[0]+costs[1]+costs[2]+costs[3];
		if (plate1 && plate2){
			cost += costs[4]; // sin
		}
		if (plate1 || plate2){
			cost += costs[5]; // distance 
		}
		double [] qualities = {
				this_rq,
				this_rq_eq,
				cost,
				costs[0]/cost,
				costs[1]/cost,
				costs[2]/cost,
				costs[3]/cost,
				costs[4]/cost,
				costs[5]/cost};
		
		if (debugLevel > 0){
			System.out.println(prefix+" cost="+cost);
			if (debugLevel > 1){
				System.out.println(prefix+ " cost contributions: "+
						((int)(100*qualities[3]))+"%, "+
						((int)(100*qualities[4]))+"%, "+
						((int)(100*qualities[5]))+"%, "+
						((int)(100*qualities[6]))+"%, "+
						((int)(100*qualities[7]))+"%, "+
						((int)(100*qualities[8]))+"%");
				System.out.println(prefix+
						" this_rq=" + this_rq+
						" this_rq_eq="+(this_rq_eq) +
						" this_rq_nofloor="+(this_rq_nofloor) +
						" this_wrq=" + (this_wrq) +
						" this_wrq_eq=" + (this_wrq_eq) +
//						" this_wrq_raw=" + (this_wrq * (w1+w2)) +
//						" this_wrq_eq_raw=" + (this_wrq_eq * (w1+w2)) +
//						" this_rq_raw="+(this_rq * (w1+w2)) +
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+" L="+merged_ev+
						" L_eq="+merged_ev_eq+
						" L1W="+plane1.getWValue()+" L2W="+plane2.getWValue()+" LW="+merged_wev+
						" L_eqW="+merged_wev_eq);
				System.out.println(prefix + ", world sin2 =" + plane1.getWorldSin2(plane2));
				System.out.println(prefix+
						", world rdist this="+ Math.sqrt(plane1.getWorldPlaneRDist2(plane2))+
						", world rdist other="+Math.sqrt(plane2.getWorldPlaneRDist2(plane1))+
						", world rdist sum="+Math.sqrt(plane1.getWorldPlaneRDist2(plane2)+
								plane2.getWorldPlaneRDist2(plane1)));
			}
		}
		return qualities; // TODO: add modes to select what is output
	}
   /**
    * Calculate what be thickness eigenvalues for merging this plane with individual neighbors
    * @param planes array of per supertile, per layer plane instances
    * @param debugLevel
    * @param dbg_X
    * @param dbg_Y
    */
	
	public void matchPlanes(
			final TilePlanes.PlaneData [][] planes,			
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		//				final int debug_stile = 20 * stilesX + 27;
		final int debug_stile = dbg_Y * stilesX + dbg_X;

		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		// Select best symmetrical match, consider only N, NE, E, SE - later opposite ones will be copied
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					TilePlanes.PlaneData [] dbg_planes = null; 
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
//						int dl = ((debugLevel > 0) && (nsTile0 == debug_stile)) ? 1:0;
                        int dl = ((debugLevel > 1) && (nsTile0 == debug_stile)) ? 3: debugLevel;
						
						if ( planes[nsTile0] != null) {
							if (dl > 1){
								System.out.println("matchPlanes(): nsTile0 ="+nsTile0);
//								dbg_planes = planes[nsTile0];
							}
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								//										planes[nsTile0][np0].initNeibBest(); // 
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
									this_plane.initMergedValue();
									for (int dir = 0; dir < 4; dir++){ // just half directions - relations are symmetrical
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										//											if ((sty < stilesY) && (sty > 0) && (stx < 0)) {
										if ((stx < stilesX) && (sty < stilesY) && (sty >= 0)) {
//											if ((stx < stilesX) && (sty < stilesY) && (sty > 0)) {
											int nsTile = sty * stilesX + stx; // from where to get
											if (nsTile >= planes.length){
												System.out.println("BUG!!!!");
											} else {
												TilePlanes.PlaneData [] other_planes = planes[nsTile];
												if (other_planes != null) {

													this_plane.initMergedValue(dir,other_planes.length); // filled with NaN
													for (int np = 0; np < other_planes.length; np ++){
														if (other_planes[np] != null) {
															TilePlanes.PlaneData other_plane = this_plane.getPlaneToThis(
																	other_planes[np],
																	dl - 2); // debugLevel);
															if (other_plane !=null) { // now always, but may add later
																TilePlanes.PlaneData merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
																		false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		plPreferDisparity, 
																		dl-2); // int       debugLevel)

																if (merged_pd !=null) { // now always, but may add later
																	///															merged_pd.scaleWeight(0.5);
																	this_plane.setNeibMatch (dir, np, merged_pd.getValue()); // smallest eigenValue
																	this_plane.setNeibWMatch(dir, np, merged_pd.getWValue()); // smallest eigenValue
																}
																if (dl > 1){
																	System.out.println("matchPlanes(): nsTile0 ="+nsTile0+":"+np0+"-("+dir+")->"+nsTile+":"+np+" (ignore_weights=true)");
																}
																merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
																		true, // false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		plPreferDisparity, 
																		dl-1); // int       debugLevel)

																if (merged_pd !=null) { // now always, but may add later
																	///															merged_pd.scaleWeight(0.5);
																	this_plane.setNeibMatchEq (dir, np, merged_pd.getValue()); // smallest eigenValue
																	this_plane.setNeibWMatchEq(dir, np, merged_pd.getWValue()); // smallest eigenValue
																}
																if (dl > 1){
																	System.out.println("matchPlanes(): nsTile0 ="+nsTile0+":"+np0+"-("+dir+")->"+nsTile+":"+np+"...DONE, merged_pd.getWValue()="+merged_pd.getWValue());
																}

															}
														}
													}
												}
											}
										}
									}
								}
							}
							if (dl > 1){
								System.out.println("matchPlanes(): nsTile0 ="+nsTile0+ " Done.");
							}

						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		// copy symmetrical relations
		ai.set(0);				
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					TilePlanes.PlaneData [][] dbg_planes = planes;
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if (dl > 1) {
							System.out.println("matchPlanes() nsTile0="+nsTile0);
						}
						if ( planes[nsTile0] != null) {
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
									for (int dir = 4; dir < 8; dir++){ // other half - copy from opposite
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
//										if ((sty < stilesY) && (sty > 0) && (stx > 0)) {
										if ((sty < stilesY) && (sty >= 0) && (stx >= 0)) {
											int nsTile = sty * stilesX + stx; // from where to get
											TilePlanes.PlaneData [] other_planes = planes[nsTile];
											if (other_planes !=null) {
												this_plane.initMergedValue(dir,other_planes.length); // filled with NaN
												for (int np = 0; np < other_planes.length; np ++){
													if (other_planes[np] != null) { // && (other_planes[np].getMergedValue(dir-4) != null)) {
														double [] nm = other_planes[np].getMergedValue(dir-4);
														if (nm != null) {
															this_plane.setNeibMatch(dir,np, nm[np0]); //
														}
														nm = other_planes[np].getMergedValueEq(dir-4);
														if (nm != null) {
															this_plane.setNeibMatchEq(dir,np, nm[np0]); //
														}
														nm = other_planes[np].getMergedWValue(dir-4);
														if (nm != null) {
															this_plane.setNeibWMatch(dir,np, nm[np0]); //
														}
														nm = other_planes[np].getMergedWValueEq(dir-4);
														if (nm != null) {
															this_plane.setNeibWMatchEq(dir,np, nm[np0]); //
														}
														
													}

												}
											}
										}
									}
								}
							}
						}
						if (dl > 1) {
							System.out.println("matchPlanes() nsTile0="+nsTile0);
						}

					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	/**
	 * Create matrix of connection costs between layer pairs for each direction. Even as currently costs seem to be symmetrical
	 * this method calculates all 8 directions, from the PoV of the current supertile.
	 * Intended to be used after first round of smoothing to re-evaluate which connections are valid and run more smoothing with remaining
	 * connections. No-Link supertiles will be set to original measured data
	 * @param planes
	 * @param debugLevel
	 * @param dbg_X
	 * @param dbg_Y
	 */
	
	public void interPlaneCosts(
			final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
			final TilePlanes.PlaneData [][] planes,			
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		//				final int debug_stile = 20 * stilesX + 27;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);

		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
//						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
                        int dl = ((debugLevel > 1) && (nsTile0 == debug_stile)) ? 3: debugLevel;
						
						if ( planes[nsTile0] != null) {
							if (dl > 1){
								System.out.println("interPlaneCosts(): nsTile0 ="+nsTile0);
							}
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								//										planes[nsTile0][np0].initNeibBest(); // 
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
									this_plane.initLinkCosts();
									for (int dir = 0; dir < 8; dir++){ // All 8 directions if the costs are not symmetrical
										int nsTile = tnSurface.getNeibIndex(nsTile0, dir);
										if (nsTile >= 0){
											TilePlanes.PlaneData [] other_planes = planes[nsTile];
											if (other_planes != null) {
												this_plane.initLinkCosts(dir,other_planes.length); // filled with NaN
												for (int np = 0; np < other_planes.length; np ++){
													if (other_planes[np] != null) {
														TilePlanes.PlaneData other_plane = this_plane.getPlaneToThis(
																other_planes[np],
																dl-2); // debugLevel);
														if (other_plane !=null) { // now always, but may add later
															double link_cost = getLinkCost(
																	en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
																	this_plane,  // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																	other_plane, // TilePlanes.PlaneData plane2,
																	Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
																	Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
																	Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
																	Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
																	"interPlaneCosts(): "+nsTile0+":"+np0+"-("+dir+")-->"+nsTile+":"+np, // String prefix,
																	dl - 2);
															this_plane.setLinkCosts(dir, np, link_cost);
														}
													}
												}
											}
										}
									}
								}
							}
							if (dl > 1){
								System.out.println("interPlaneCosts(): nsTile0 ="+nsTile0+ " Done.");
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
/*
	public void setExclusiveLinks(
			final TilePlanes.PlaneData [][] planes,
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		setExclusiveLinks(
				planes,
				plExNeibCost,
				debugLevel,
				dbg_X,
				dbg_Y);
	}
*/
	/**
	 * Set links between planes in 8 directions, in such a way that no start/end planes can be shared
	 * @param planes per supertile, per plane - array of supertile instances
	 * @param max_cost maximal composite connection cost allowed
	 * @param debugLevel
	 * @param dbg_X
	 * @param dbg_Y
	 */
	public void setExclusiveLinks(
			final TilePlanes.PlaneData [][] planes,
			final double                    max_cost,
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		//				final int debug_stile = 20 * stilesX + 27;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);

		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		// Reset neighbors
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						if ( planes[nsTile0] != null) {
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
									this_plane.initNeibBest();
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);				
		final class LinkPair{
			int src;
			int dst;
			double cost;
			LinkPair (int src, int dst, double cost){
				this.src = src;
				this.dst = dst;
				this.cost = cost;
			}
			boolean hasConflict(LinkPair lp){
				return (src == lp.src) || (dst == lp.dst);
			}
			public String toString(){
				return String.format("%d-(%6.3f)->%d", src,cost,dst);
			}
			
		}
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
//						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
                        int dl = ((debugLevel > 1) && (nsTile0 == debug_stile)) ? 3: debugLevel;
						
						if ( planes[nsTile0] != null) {
							int num_sp = planes[nsTile0].length; 
							if (dl > 1){
								System.out.println("setExclusiveLinks(): nsTile0 ="+nsTile0);
							}
							for (int dir = 0; dir < 4; dir++){//
								int dir_back = (dir+4) % 8;
								int nsTile = tnSurface.getNeibIndex(nsTile0, dir);
								if ((nsTile >= 0) && (planes[nsTile] != null)){
									int num_dp = planes[nsTile].length;
									double [][] inter_costs = new double [num_sp][num_dp];
									for (int np0 = 0; np0 < num_sp; np0++) if (planes[nsTile0][np0]!=null){
										for (int np = 0; np < num_dp; np++)if (planes[nsTile][np]!=null){
											inter_costs[np0][np] = 0.5* (
													planes[nsTile0][np0].getLinkCosts(dir, np)+
													planes[nsTile][np].getLinkCosts(dir_back, np0));
										}
									}
									HashSet<LinkPair> pairs_set = new HashSet<LinkPair>();
									for (int np0 = 0; np0 < num_sp; np0++) if (planes[nsTile0][np0]!=null){
										for (int np = 0; np < num_dp; np++)if (planes[nsTile][np]!=null){
											if (!Double.isNaN(inter_costs[np0][np]) && (inter_costs[np0][np] <= max_cost)){
												pairs_set.add(new LinkPair(np0, np, inter_costs[np0][np]));
											}
										}										
									}
									HashSet<LinkPair> links_set = new HashSet<LinkPair>();
									while (!pairs_set.isEmpty()){
										LinkPair link =Collections.min(pairs_set, new Comparator<LinkPair>() {
											@Override
											public int compare(LinkPair lhs, LinkPair rhs) { // ascending
												return (lhs.cost < rhs.cost) ? -1 : (lhs.cost > rhs.cost ) ? 1 : 0;
											}
										});
										links_set.add(link);
										pairs_set.remove(link);
										HashSet<LinkPair> conflict_links = new HashSet<LinkPair>();
										for (LinkPair lp:pairs_set){
											if (link.hasConflict(lp)) {
												conflict_links.add(lp);
											}
										}
										pairs_set.removeAll(conflict_links);
									}
									for (LinkPair lp:links_set){
										planes[nsTile0][lp.src].setNeibBest(dir, lp.dst);
										planes[nsTile ][lp.dst].setNeibBest(dir_back, lp.src);
									}
								}
							}							
							if (dl > 1){
								System.out.println("setExclusiveLinks(): nsTile0 ="+nsTile0+ " Done.");
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
	
	/**
	 * Mark which links between neighbor planes are valid
	 * @param debugLevel debug level
	 * @param dbg_X debug supertile X coordinate
	 * @param dbg_Y debug supertile Y coordinate
	 */

	public void filterNeighborPlanes(
			final TilePlanes.PlaneData [][] planes,
			final boolean merge_low_eigen,
			final int     debugLevel,
			final int     dbg_X,
			final int     dbg_Y)
	{
		final double               relax = 1.0;

		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		//				final int debug_stile = 20 * stilesX + 27;
		//				final int debug_stile = 17 * stilesX + 27;
		//				final int debug_stile = 9 * stilesX + 26;
		final int debug_stile = dbg_Y * stilesX + dbg_X;

		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					TilePlanes.PlaneData [][] dbg_planes = planes;
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = (nsTile0 == debug_stile) ? debugLevel:0;
						if ( planes[nsTile0] != null) {
							if (dl > 0){
								System.out.println("filterNeighborPlanes() nsTile0="+nsTile0);
							}
							int np0_min = (planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
							for (int dir = 0; dir < 4; dir++){ //
								int stx = stx0 + dirsYX[dir][1];
								int sty = sty0 + dirsYX[dir][0];
								int nsTile = sty * stilesX + stx; // from where to get

								for (int np0 = np0_min; np0 < planes[nsTile0].length; np0++){
									if ((planes[nsTile0][np0] != null) && (planes[nsTile0][np0].getMergedValue(dir) != null)){
										double [] merge_ev =     planes[nsTile0][np0].getMergedValue(dir);
										double [] merge_ev_eq =  planes[nsTile0][np0].getMergedValueEq(dir);
										double [] merge_wev =    planes[nsTile0][np0].getMergedWValue(dir);
										double [] merge_wev_eq = planes[nsTile0][np0].getMergedWValueEq(dir);
										
										if (	(merge_ev != null) &&
												(merge_ev_eq != null)) {
											int np_min = SuperTiles.LOWEST_PLANE(merge_ev.length);
											for (int np = np_min; np < merge_ev.length; np++){
												if (	(planes[nsTile][np] != null) &&
														!Double.isNaN(merge_ev[np])) {
													if (dl > 0){
														System.out.println("filterNeighborPlanes() nsTile0="+nsTile0+" np0="+np0+"-("+dir+")-> nsTile="+nsTile+" np="+np);
													}
													String prefix = "filterNeighborPlanes() nsTile0="+nsTile0+" np0="+np0+" dir="+dir+" nsTile="+nsTile+" np="+np;
													if (planesFit(
															planes[nsTile0][np0], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
															planes[nsTile][np],   // 			TilePlanes.PlaneData plane2,
															relax, // double               relax = 1.0;
															true, // merge_low_eigen, // false,                // boolean              merge_weak,   // use for same supertile merge 
															merge_ev[np],         // double               merged_ev,    // if NaN will calculate assuming the same supertile
															merge_ev_eq[np],      //double               merged_ev_eq, // if NaN will calculate assuming the same supertile
															merge_wev[np],        // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
															merge_wev_eq[np],     // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
															prefix,               // String prefix,
															dl)                   // int debugLevel)
															){
														planes[nsTile0][np0].setMergedValid(dir, np, true, planes[nsTile].length);
														planes[nsTile][np].setMergedValid((dir + 4) %8, np0, true, planes[nsTile0].length);
														if (planesFit(
																planes[nsTile0][np0], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																planes[nsTile][np],   // 			TilePlanes.PlaneData plane2,
																relax, // double               relax = 1.0;
																false, // merge_low_eigen, // false,                // boolean              merge_weak,   // use for same supertile merge 
																merge_ev[np],         // double               merged_ev,    // if NaN will calculate assuming the same supertile
																merge_ev_eq[np],      //double               merged_ev_eq, // if NaN will calculate assuming the same supertile
																merge_wev[np],        // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
																merge_wev_eq[np],     // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
																prefix,               // String prefix,
																dl)                   // int debugLevel)
																){
															planes[nsTile0][np0].setMergedStrongValid(dir, np, true, planes[nsTile].length);
															planes[nsTile][np].setMergedStrongValid((dir + 4) %8, np0, true, planes[nsTile0].length);
														}
													}
												}															
											}
										}
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	/**
	 * Merge the supertile planes with agreeing neighbors, non-exclusively (no considering other planes
	 * of the same supertile. Start with the best fit, then goes to lower quality, until the individual
	 * merge quality falls below scaled quality of the best, pre-set minimum or the merged plane becomes
	 * too thick
	 * Separately calculates merged weighted plane and with equal weights of the neighbors
	 * TODO: Maybe just switch to the connection costs instead? 
	 * @param planes array of plane instances for the same supertile
	 * @param debugLevel
	 * @param dbg_X
	 * @param dbg_Y
	 */

	public void setNonExclusive(
			final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
			final TilePlanes.PlaneData [][] planes,
//			final double  center_weight,
			final int     debugLevel,
			final int     dbg_X,
			final int     dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					TilePlanes.PlaneData [][] dbg_planes = planes;
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int dl = (nsTile0 == debug_stile) ? debugLevel:0;
						if ( planes[nsTile0] != null) {
							if (dl > 0){
								System.out.println("setNonExclusive() nsTile0="+nsTile0);
							}
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++) if (planes[nsTile0][np0] != null) {
								// find best connection quality for each direction
								int [] neibs = {-1,-1,-1,-1,-1,-1,-1,-1};
								double [][] costs = new double[neibs.length][];
								double [] weights = new double[neibs.length];
								TilePlanes.PlaneData [] neib_pd = new TilePlanes.PlaneData[neibs.length];

								int cost_index = 2;
								double sum_rcosts = 0.0;
								int non_zero = 0;
								for (int dir = 0; dir < neibs.length; dir++) if (planes[nsTile0][np0].hasMergedValid(dir)){
									int nsTile1 = tnSurface.getNeibIndex(nsTile0, dir);
									if (nsTile1 >= 0){
										boolean [] merged_valid = planes[nsTile0][np0].getMergedValid(dir);
										double [][] costs_dir = new double[merged_valid.length][];
										int best_np = -1;
										for (int np = 0; np < merged_valid.length; np++) if (merged_valid[np] && (planes[nsTile1][np] != null)){
											TilePlanes.PlaneData other_plane = planes[nsTile0][np0].getPlaneToThis(
													planes[nsTile1][np], // neighbor, previous value
													debugLevel - 2); // debugLevel);
											if (other_plane != null) {
												costs_dir[np] = getFitQualities(
														en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
														planes[nsTile0][np0], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
														other_plane,          // TilePlanes.PlaneData plane2,
														Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
														Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
														Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
														Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
														nsTile0+":"+np0+"-"+dir+":"+np,  // String prefix,
														0); // int debugLevel)
												if (costs_dir[np] != null) {
													for (int i = 3; i < costs_dir[np].length; i++) costs_dir[np][i] *= costs_dir[np][2];
													if ((best_np < 0) || (costs_dir[np][cost_index] < costs_dir[best_np][cost_index])){
														best_np = np;
													}
												}
											}
										}
										if ((best_np >= 0) && (costs_dir[best_np][cost_index] < plNeNeibCost)) {
											neibs[dir] = best_np;
											costs[dir] = costs_dir[best_np];
											non_zero ++;
											weights[dir] = 1.0/costs[dir][cost_index];
											sum_rcosts += weights[dir];
											neib_pd[dir] = planes[nsTile1][best_np];
										}
									}
								}
								
								
								
								for (int dir = 0; dir < neibs.length; dir++) if (neibs[dir] >= 0) {
									weights[dir] *= non_zero/sum_rcosts; // average weight for active links will be 1.0 
								}								
								if (dl > 0) {
									for (int dir = 0; dir <8; dir++)if (costs[dir] != null){
										System.out.print(nsTile0+":"+np0+"-"+dir+"  "+String.format(" weight= %6.3f ",weights[dir]));
										for (int i = 0; i < costs[dir].length; i++){
											System.out.print(String.format("%8.3f", costs[dir][i]));
										}
										System.out.println();
									}
									System.out.println();
								}
								// maybe remove outliers here?
								if (non_zero > 0) {
									TilePlanes.PlaneData merged_pd = mergeWithNeibs(
											planes[nsTile0][np0], // TilePlanes.PlaneData   center_pd,
											neib_pd,              // TilePlanes.PlaneData[] neibs_pd,
											plMaxEigen,           // double  maxValue, // do not combine with too bad planes, do not merge if result is above
											plNeOwn,              // double                 center_weight,
											false,                // boolean                stopBadNeib,
											1.0,                  // final double           normPow,
											weights,              // double []              weights,
											false, // boolean                ignore_weights,
											dl); // int debugLevel)
									
									merged_pd.getWorldXYZ(0); // debugLevel); // just to recalculate world data for debugging
									merged_pd.setNeibBest(neibs);
									planes[nsTile0][np0].setNonexclusiveStar(merged_pd);
									TilePlanes.PlaneData merged_pd_eq = mergeWithNeibs(
											planes[nsTile0][np0], // TilePlanes.PlaneData   center_pd,
											neib_pd,              // TilePlanes.PlaneData[] neibs_pd,
											plMaxEigen,           // double  maxValue, // do not combine with too bad planes, do not merge if result is above
											plNeOwn,              // double                 center_weight,
											false,                // boolean                stopBadNeib,
											1.0,                  // final double           normPow,
											weights,              // double []              weights,
											true, // boolean                ignore_weights,
											dl); // int debugLevel)
									
									merged_pd_eq.getWorldXYZ(0); // debugLevel); // just to recalculate world data for debugging
									merged_pd_eq.setNeibBest(neibs);
									planes[nsTile0][np0].setNonexclusiveStarEq(merged_pd_eq);
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
	
	// TODO: reuse in planesSmoothStep
	TilePlanes.PlaneData mergeWithNeibs(
			TilePlanes.PlaneData   center_pd,
			TilePlanes.PlaneData[] neibs_pd,
			double                 maxValue, // do not combine with too bad planes, do not merge if result is above
			double                 center_weight,
			boolean                stopBadNeib,
			final double           normPow,
			double []              weights,
			boolean                ignore_weights,
			int debugLevel)
	{
		TilePlanes.PlaneData this_new_plane = center_pd.clone();

		this_new_plane.setWeight(0.0); //
		double num_merged = 0.0; // double to add fractional pull weight of the center
		double true_num_merged = 0.0;
		for (int dir = 0; dir < neibs_pd.length; dir++){
			if (neibs_pd[dir] != null) {
				TilePlanes.PlaneData other_plane = this_new_plane.getPlaneToThis(
						neibs_pd[dir], // neighbor, previous value
						debugLevel - 2); // debugLevel);
				//				if (dl > 0) dbg_img[ 2 + dir] = other_plane.getSinglePlaneDisparity(false);
				if ((other_plane != null) && ((other_plane.getValue() <= maxValue) || (maxValue == 0))) { // TODO:
					if (this_new_plane.getWeight() > 0.0){
						this_new_plane = this_new_plane.mergePlaneToThis(
								other_plane, // PlaneData otherPd,
								weights[dir], // 1.0,         // double    scale_other,
								// here it should be no power function for the weights
								1.0,               // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
								ignore_weights,    // boolean   ignore_weights,
								true,              // boolean   sum_weights,
								center_pd.getPreferDisparity(),
								debugLevel - 2);   // int       debugLevel)

					} else {
						this_new_plane = other_plane; // should increment num_merged
						this_new_plane.scaleWeight(weights[dir]);
					}
					if (this_new_plane != null){
						num_merged += 1.0;
						true_num_merged += 1.0;
						// just for debug / calculate
						this_new_plane.getWorldXYZ(0);
					}
					//					new_planes[nsTile0][np0] = this_new_plane;
					//					if (dl > 0) dbg_img[10 + dir] = this_new_plane.getSinglePlaneDisparity(false);
				} else if (stopBadNeib){ // just skip, not abandon? (set to false)
					this_new_plane = null;
					break;
				}
			}
		}

		if (this_new_plane != null) {
			// average weight over participating directions, so the relative pull
			// does not depend on number of neighbors
			if ((num_merged > 0.0) && (this_new_plane != null) && (normPow > 1.0)){
				double scale = Math.pow(num_merged, normPow);
				this_new_plane.scaleWeight(1.0/scale);
				num_merged /=scale;
				true_num_merged /= scale;
			}
			
			
			if (    (center_weight > 0.0) &&
					((center_pd.getValue() < maxValue) || (maxValue == 0) || 
							(this_new_plane == null) || (this_new_plane.getWeight() == 0.0)) // keep measured if impossible to merge
					){

				if (this_new_plane.getWeight() > 0.0){
					this_new_plane = this_new_plane.mergePlaneToThis(
							center_pd, // PlaneData otherPd,
							center_weight,                     // double    scale_other,
							1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost													
							false,                         // boolean   ignore_weights,
							true, // boolean   sum_weights,
							center_pd.getPreferDisparity(),
							debugLevel - 2);                           // int       debugLevel)
					if (this_new_plane != null){
						num_merged += center_weight;	// num_merged was 1.0 and weight is averaged over all neighbors
						true_num_merged += 1.0;
					}
				} else {
					this_new_plane = center_pd.clone();
					num_merged = 1.0;
					true_num_merged = 1.0;
				}
//				if (dl > 0) dbg_img[18] =  this_new_plane.getSinglePlaneDisparity(false);
			}
		}
		if ((num_merged > 0.0) && (this_new_plane != null)){
			// just for debug / calculate
			this_new_plane.getWorldXYZ(0);
			this_new_plane.scaleWeight(1.0/num_merged);
//			double true_num_merged = num_merged - center_weight + 1;
			this_new_plane.setNumPoints((int) Math.round(this_new_plane.getNumPoints()/true_num_merged));
		}
		
		return this_new_plane;
	}
	
	
	
	/**
	 * Find mutual links between multi-layer planes for supertiles. requires that for each plane there are calculated smalles eigenvalues
	 * for merging with each plane for each of 8 neighbors 
	 * @param debugLevel debug level
	 * @param dbg_X debug supertile X coordinate
	 * @param dbg_Y debug supertile Y coordinate
	 */
	public double [][]  selectNeighborPlanesMutual(
			final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
			final TilePlanes.PlaneData [][] planes,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final double [] nan_plane = new double [superTileSize*superTileSize];
		for (int i = 0; i < nan_plane.length; i++) nan_plane[i] = Double.NaN;
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final int numThreads = threads.length;
		final double [][][] all_quality_stats = new double  [numThreads][2][6]; // contributions of all [0] and winners [2], 6 parameters 
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						if ( planes[nsTile0] != null) {
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
									this_plane.initNeibBest();
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		ai.set(0);				
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
					double [][] quality_stats = all_quality_stats[numThread];
//					TilePlanes.PlaneData [][] dbg_planes = planes;
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = (nsTile0 == debug_stile) ? debugLevel:0;
						if ( planes[nsTile0] != null) {
							if (dl > 0){
								System.out.println(" ===== selectNeighborPlanesMutual() nsTile0="+nsTile0+" =====");
							}
							int np0_min = (planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
							for (int dir = 0; dir < 4; dir++){ //
								int stx = stx0 + dirsYX[dir][1];
								int sty = sty0 + dirsYX[dir][0];
								int nsTile = sty * stilesX + stx; // from where to get

								int num_other_planes = 0;
								for (int np0 = np0_min; np0 < planes[nsTile0].length; np0++){
									if ((planes[nsTile0][np0] != null) && (planes[nsTile0][np0].getMergedValue(dir) != null)){
										int l = planes[nsTile0][np0].getMergedValue(dir).length;
										if (l > num_other_planes)num_other_planes = l;
									}
								}
								if (num_other_planes > 0){ // will eliminate bad margins
									int np_min = SuperTiles.LOWEST_PLANE(num_other_planes);
									
									boolean [] this_matched = new boolean [planes[nsTile0].length];
									boolean [] other_matched = new boolean [num_other_planes];
									int num_pairs = this_matched.length - np0_min;
									if ((other_matched.length - np_min) < num_pairs) num_pairs = other_matched.length - np_min;
									double [][][] qualities = new double [this_matched.length][][];

									for (int pair = 0; pair < num_pairs; pair ++){
										int [] best_pair = {-1,-1};
										double best_rqual = Double.NaN;
										for (int np0 = np0_min; np0 < this_matched.length; np0++) if (planes[nsTile0][np0] != null){
											double [] merge_ev =    planes[nsTile0][np0].getMergedValue(dir);
											double [] merge_ev_eq = planes[nsTile0][np0].getMergedValueEq(dir);
											double [] merge_wev =    planes[nsTile0][np0].getMergedWValue(dir);
											double [] merge_wev_eq = planes[nsTile0][np0].getMergedWValueEq(dir);
											boolean [] merge_valid = planes[nsTile0][np0].getMergedValid(dir);
											qualities[np0] = new double[ merge_ev.length][];
											
											if (!this_matched[np0] &&(merge_valid != null)) {
												for (int np = np_min; np < merge_ev.length; np++){
													if (!other_matched[np] && merge_valid[np]) {
														String prefix = "selectNeighborPlanesMutual(): nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np;
														qualities[np0][np] = getFitQualities( // {this_rq, this_rq_eq}; 
																en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
																planes[nsTile0][np0], //TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																planes[nsTile][np], //TilePlanes.PlaneData plane2,
																merge_ev[np],       // double merged_ev,    // if NaN will calculate assuming the same supertile
																merge_ev_eq[np],    // double merged_ev_eq, // if NaN will calculate assuming the same supertile
																merge_wev[np],      // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
																merge_wev_eq[np],   // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
																prefix, // String prefix,
																dl); // int debugLevel)
														if (qualities != null) {
//															double this_rq = qualities[np0][np][0]; // just for old compatibility - CHANGE later
															double this_rq = qualities[np0][np][2]; // just for old compatibility - CHANGE later
															// statistics for all
															for (int i = 0; i < quality_stats[0].length; i++){
																quality_stats[0][i] += qualities[np0][np][3 + i];
															}
															
															if (Double.isNaN(best_rqual) || (this_rq < best_rqual)){ // OK if Double.isNaN(this_rq[np])
																if (dl > 0){
																	System.out.println(" ===== selectNeighborPlanesMutual) : nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np+", this_rq="+this_rq);
																}
																best_rqual = this_rq;
																best_pair[0]= np0;
																best_pair[1]= np;
															}
														}
													}															
												}														
											}
										}
										// calculate statistics
//										for (int pair = 0; pair < num_pairs; pair ++){
//											
//										}
										
										if (Double.isNaN(best_rqual)){
											if (dl >0) {
												System.out.println("selectNeighborPlanesMutual - nothing more found for "+nsTile0+" dir = " + dir);
											}
											break; // nothing (more) found
										}
										// statistics for the winner
										for (int i = 0; i < quality_stats[0].length; i++){
											quality_stats[1][i] += qualities[best_pair[0]][best_pair[1]][3 + i];
										}

										
										this_matched[best_pair[0]] = true;
										other_matched[best_pair[1]] = true;
										// neib_best should be initialized as  int [8];
										//												planes[nsTile0][0].initNeibBest();
										if (dl >0) {
											System.out.println("1. planes["+nsTile0+"]["+best_pair[0]+"].setNeibBest("+dir+","+best_pair[1]+")");
										}
										if (dl >0) {
											System.out.println("2. planes["+nsTile+"]["+best_pair[1]+"].setNeibBest("+(dir+4)+","+best_pair[0]+")");
										}

										planes[nsTile0][best_pair[0]].setNeibBest(dir,best_pair[1]);
										planes[nsTile][best_pair[1]].setNeibBest(dir + 4,best_pair[0]);
									}
									// disable remaining neighbors
									for (int np = 0; np < this_matched.length; np++) if (planes[nsTile0][np] != null){
										if (dl >0) {
											System.out.println("this_matched["+np+"]="+this_matched[np]);
										}
										if (!this_matched[np]){
											planes[nsTile0][np].setNeibBest(dir,-1);
										}
									}
//									for (int np = 0; np < other_matched.length; np++) if (planes[nsTile][np] != null){
									for (int np = 0; np < other_matched.length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
										if (!other_matched[np]){
											if (dl >0) {
												System.out.println("other_matched["+np+"]="+other_matched[np]);
											}
											planes[nsTile][np].setNeibBest(dir + 4,-1);
										}
									}
									if (dl >0) {
										for (int np = 0; np < this_matched.length; np++) if (planes[nsTile0][np] != null){
											int [] bn = planes[nsTile0][np].getNeibBest();
											System.out.println("nsTile0="+nsTile0+":"+np+" : ["+
											bn[0]+","+bn[1]+","+bn[2]+","+bn[3]+","+bn[4]+","+bn[5]+","+bn[6]+","+bn[7]+"]");
										}
//										for (int np = 0; np < other_matched.length; np++) if (planes[nsTile][np] != null){
										for (int np = 0; np < other_matched.length; np++) if ((planes[nsTile][np] != null) && (planes[nsTile][np].getWeight() > 0.0)){ // disregard 0-weight planes
											int [] bn = planes[nsTile][np].getNeibBest();
											System.out.println("nsTile="+nsTile+":"+np+" best neighbors : ["+
											bn[0]+","+bn[1]+","+bn[2]+","+bn[3]+","+bn[4]+","+bn[5]+","+bn[6]+","+bn[7]+"]");
										}
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		double [][] quality_stats = new double  [2][6]; // contributions of all [0] and winners [2], 4 parameters
		for (int n = 0; n < all_quality_stats.length; n++){
			for (int i = 0; i < all_quality_stats[n].length; i++){
				for (int j = 0; j < all_quality_stats[n][i].length; j++){
					quality_stats[i][j] += all_quality_stats[n][i][j];
				}
			}
		}
		// normalize
		for (int i = 0; i < quality_stats.length; i++){
			double s=0.0;
			for (int j = 0; j < quality_stats[i].length; j++){
				s += quality_stats[i][j];
			}
			if (s != 0) {
				for (int j = 0; j < quality_stats[i].length; j++){
					quality_stats[i][j] /= s;
				}
			}
		}
		if (debugLevel > -1){
			System.out.println("Contribution of various factors for all considered pairs and the winners:");
			System.out.println("    weighted quality (disp) : "+quality_stats[0][0]+" (all), "+quality_stats[1][0]+" (winners)");
			System.out.println("non-weighted quality (disp) : "+quality_stats[0][1]+" (all), "+quality_stats[1][1]+" (winners)");
			System.out.println("    weighted quality (world): "+quality_stats[0][2]+" (all), "+quality_stats[1][2]+" (winners)");
			System.out.println("non-weighted quality (world): "+quality_stats[0][3]+" (all), "+quality_stats[1][3]+" (winners)");
			System.out.println("                        sin2: "+quality_stats[0][4]+" (all), "+quality_stats[1][4]+" (winners)");
			System.out.println("                      rdist2: "+quality_stats[0][5]+" (all), "+quality_stats[1][5]+" (winners)");
		}
		return quality_stats;

		
	}
	
	/**
	 * Select initial pairs of plane merge candidates (same supertile plane). All pairs are ordered, first is lower than second
	 * Candidates for merging share some non-exclusive neighbors
	 * @param planes per supertile, per plane - array of supertile instances
	 * @param debugLevel debug level
	 * @param dbg_X tile x-index for detailed debug data
	 * @param dbg_Y tile y-index for detailed debug data
	 * @return array of per-supertile, per pair {start plane, end_plane} pairs
	 */
	public int [][][] getMergeSameTileCandidates(
			final TilePlanes.PlaneData [][] planes,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final int [][][] merge_candidates = new int [nStiles][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		class LayersLinks{
			int nl1, nl2, links1,  links2, shared;
			LayersLinks (int nl1, int nl2, int links1, int links2, int shared){
				this.nl1 =    nl1;
				this.nl2 =    nl2;
				this.links1 = links1;
				this.links2 = links2;
				this.shared = shared;
			}
			int [] toArray()
			{
				int [] data = {nl1, nl2, links1,  links2, shared};
				return data;
			}
		}
		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( planes[nsTile] != null) {
						ArrayList<LayersLinks> links_list = new ArrayList<LayersLinks>();
						int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? 1:0;
						if (dl > 0){
							System.out.println("getMergeCandidates(): nsTile="+nsTile);
						}
						for (int np1 = 0; np1 < planes[nsTile].length; np1++) if (planes[nsTile][np1] != null){ // nu
							boolean [][] merged_valid1 = planes[nsTile][np1].getMergedValid();
							if (merged_valid1 != null){
								for (int np2 = np1 + 1; np2 < planes[nsTile].length; np2++) if (planes[nsTile][np2] != null){ // nu
									boolean [][] merged_valid2 = planes[nsTile][np2].getMergedValid();
									if (merged_valid2 != null){
										int num_links1 = 0;
										int num_links2 = 0;
										int num_shared = 0;
										for (int dir = 0; dir < 8; dir++){
											if (merged_valid1[dir] != null){
												for (int nl = 0; nl < merged_valid1[dir].length; nl++){
													if (merged_valid1[dir][nl]) num_links1++;
												}
											}
											if (merged_valid2[dir] != null){
												for (int nl = 0; nl < merged_valid2[dir].length; nl++){
													if (merged_valid2[dir][nl]) num_links2++;
												}
											}
											if ((merged_valid1[dir] != null) && (merged_valid2[dir] != null)) { // should be the same length
												for (int nl = 0; nl < merged_valid2[dir].length; nl++){
													if (merged_valid1[dir][nl] && merged_valid2[dir][nl]) num_shared++;
												}
											}
										}
										if (num_shared > 0) links_list.add(new LayersLinks(np1, np2, num_links1, num_links2, num_shared));
									}
								}
							}
						}
						if (!links_list.isEmpty()){
							merge_candidates[nsTile] = new int [links_list.size()][];
							int indx = 0;
							for (LayersLinks ll : links_list){
								merge_candidates[nsTile][indx++] = ll.toArray();
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debugLevel > 1){
			System.out.println("Supertile planes that are candidates for merging:");
			for (int nsTile = 0; nsTile < nStiles; nsTile++){
				int stx = nsTile % stilesX;
				int sty = nsTile / stilesX;
				if (merge_candidates[nsTile] != null){
					for (int i = 0 ; i < merge_candidates[nsTile].length; i++){
						double sharedRatio = 2.0 * merge_candidates[nsTile][i][4] / (merge_candidates[nsTile][i][2] + merge_candidates[nsTile][i][3]);
						System.out.println(nsTile+" ["+stx+":"+sty+"] ("+merge_candidates[nsTile][i][0]+", "+merge_candidates[nsTile][i][1]+")"+
						" shared "+(((int) (sharedRatio * 1000)) / 10) + "%" + 
						" links1 = "+merge_candidates[nsTile][i][2]+
						" links2 = "+merge_candidates[nsTile][i][3]+
						" shared links = "+merge_candidates[nsTile][i][4]);
					}
				}
			}
		}
		return merge_candidates;
	}

	public  int [][][] filterMergeSameTileCandidates(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
//		final boolean [][] merge_pairs = new boolean [nStiles][];
		final int [][][]  filtered_merge_candidates = new int [nStiles][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( merge_candidates[nsTile] != null) {
//						merge_pairs[nsTile] = new boolean [merge_candidates[nsTile].length];
						ArrayList<Integer> filtered_pairs_list = new ArrayList<Integer>();
						int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? 2: ((debugLevel > 1) ? 1:0);
						if (dl > 1){
							System.out.println("filterMergeSameTileCandidates(): nsTile="+nsTile);
						}
						int n_planes = planes[nsTile].length;
						double [][] plane_strengths = new double [n_planes * 2][];
						for (int np = 0; np < planes[nsTile].length; np++) if (planes[nsTile][np] != null){
							double [][] ds = planes[nsTile][np].getDoublePlaneDisparityStrength(
									false, // boolean   useWorld,
									null, // double [] window,
									-1, // int       dir,
									false, // boolean   use_sel,
									false, // boolean   divide_by_area,
									1.0, // double    scale_projection,
									0.0, // double    fraction_uni,
									0); // int       debugLevel)
							double weight = planes[nsTile][np].getWeight();
							if (weight != 0){ // can it be? probably not
								for (int i = 0; i < ds[1].length; i++){
									ds[1][i]/=weight;
								}
							}
							plane_strengths[np] = ds[1];
							plane_strengths[np + n_planes] = new double [ds[1].length];
							for (int i = 0; i < ds[1].length; i++){
								plane_strengths[np + n_planes][i] = (ds[1][i] > 0.5)? 1.0 : 0.0;
							}
						}
						if ((dl > 1) && (nsTile == debug_stile)){
							System.out.println("filterMergeSameTileCandidates().2: nsTile="+nsTile);
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
							sdfa_instance.showArrays(plane_strengths, 2 * superTileSize, 2* superTileSize, true, "filterMergeSameTileCandidates_"+nsTile);
						}
						
						
						
						for (int pair = 0; pair < merge_candidates[nsTile].length; pair ++){
							int np1 =  merge_candidates[nsTile][pair][0];
							int np2 =  merge_candidates[nsTile][pair][1];
							if ((plane_strengths[np1] != null) && (plane_strengths[np2] != null)){
								int [] counts = {0, 0, 0 };
								for (int i = 0; i < plane_strengths[np1].length; i++){
									if (plane_strengths[np1][i] > 0.5) counts[0]++;
									if (plane_strengths[np2][i] > 0.5) counts[1]++;
									if ((plane_strengths[np1][i] > 0.5) && (plane_strengths[np2][i] > 0.5)) counts[2]++;
								}
								double [] overlaps = {1.0 * counts[2]/counts[0], 1.0 * counts[2]/counts[1]};
								if ((overlaps[0] < plMaxOverlap) && (overlaps[0] < plMaxOverlap)) {
									filtered_pairs_list.add(pair);
									if (debugLevel > 1){
										System.out.println("filterMergeSameTileCandidates(): VERIFIED pair nsTile="+nsTile+":"+np1+":"+np2+
												" as it has LOW enough overlap: "+
												" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
												" overlap2="+ ((int) (100 *overlaps[1]))+"% ");
									}
									
								} else {
									if (debugLevel > 0){
										System.out.println("filterMergeSameTileCandidates(): REMOVED pair nsTile="+nsTile+":"+np1+":"+np2+
												" as it has HIGH overlap: "+
												" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
												" overlap2="+ ((int) (100 *overlaps[1]))+"% ");
									}
								}
							}
						}
						if (!filtered_pairs_list.isEmpty()){
							filtered_merge_candidates[nsTile] = new int [filtered_pairs_list.size()][];
							int indx = 0;
							for (Integer pair:filtered_pairs_list){
								filtered_merge_candidates[nsTile][indx++] = merge_candidates[nsTile][pair];
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return filtered_merge_candidates;
	}

	/**
	 * Verify merge candidates by evaluating overlaps. Close planes that have significant overlap on the image are mofre likely
	 * to be actually different planes, if they do not overlap it is likely to be the same one, just multiple separate parts of
	 * it mistakenly discriminated during initial plane generation from the tiles data.
	 * 
	 * Merging is also permitted if the planes are "weak and similar" - low strength and fit for merging 
	 * 
	 * Added a 'hack' - allow merging even overlapping tiles if they are close enough in the real world space
	 * 
	 * @param planes per supertile, per plane - array of supertile instances
	 * @param merge_candidates array of merge pairs (per-supertile, per pair {start plane, end_plane} pairs) 
	 * @param min_distance minimal distance (in meters) to keep merging overlapping pair
	 * @param debugLevel debug level
	 * @param dbg_X tile x-index for detailed debug data
	 * @param dbg_Y tile y-index for detailed debug data
	 * @return per super tile , per first plane index, per second plane index boolean array of permitted to merge
	 */
	public  boolean [][][] overlapSameTileCandidates(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			// may be a hack - has problems with not merging some close portions of the (horizontal) pavement - 
			// maybe just some dirty window glitches
			// ignore overlap if the planes are close in the real world. Only comparing distances from good plates
			final double     min_distance, 
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final double min_distance2 = min_distance * min_distance;
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
//		final boolean [][] merge_pairs = new boolean [nStiles][];
		final boolean [][][]  overlap_merge_candidates = new boolean [nStiles][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( merge_candidates[nsTile] != null) {
						int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? 2: ((debugLevel > 1) ? 1:0);
						if (dl > 1){
							System.out.println("overlapSameTileCandidates(): nsTile="+nsTile);
						}
						int n_planes = planes[nsTile].length;
						overlap_merge_candidates[nsTile] = new boolean [n_planes][n_planes];
						double [][] plane_strengths = new double [n_planes * 2][];
						for (int np = 0; np < planes[nsTile].length; np++) if (planes[nsTile][np] != null){
							double [][] ds = planes[nsTile][np].getDoublePlaneDisparityStrength(
									false, // boolean   useWorld,
									null, // double [] window,
									-1, // int       dir,
									false, // boolean   use_sel,
									false, // boolean   divide_by_area,
									1.0, // double    scale_projection,
									0.0, // double    fraction_uni,
									0); // int       debugLevel)
							double weight = planes[nsTile][np].getWeight();
							if (weight != 0){ // can it be? probably not
								for (int i = 0; i < ds[1].length; i++){
									ds[1][i]/=weight;
								}
							}
							plane_strengths[np] = ds[1];
							plane_strengths[np + n_planes] = new double [ds[1].length];
							for (int i = 0; i < ds[1].length; i++){
								plane_strengths[np + n_planes][i] = (ds[1][i] > 0.5)? 1.0 : 0.0;
							}
						}
						if ((dl > 1) && (nsTile == debug_stile)){
							System.out.println("overlapSameTileCandidates().2: nsTile="+nsTile);
							showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
							sdfa_instance.showArrays(plane_strengths, 2 * superTileSize, 2* superTileSize, true, "overlapSameTileCandidates_"+nsTile);
						}
						for (int np1 = 0; np1 < planes[nsTile].length; np1++) if (planes[nsTile][np1] != null){
							for (int np2 = np1 + 1; np2 < planes[nsTile].length; np2++) if (planes[nsTile][np2] != null){
								int [] counts = {0, 0, 0 };
								for (int i = 0; i < plane_strengths[np1].length; i++){
									if (plane_strengths[np1][i] > 0.5) counts[0]++;
									if (plane_strengths[np2][i] > 0.5) counts[1]++;
									if ((plane_strengths[np1][i] > 0.5) && (plane_strengths[np2][i] > 0.5)) counts[2]++;
								}
								double [] overlaps = {1.0 * counts[2]/counts[0], 1.0 * counts[2]/counts[1]};
								if ((overlaps[0] < plMaxOverlap) && (overlaps[0] < plMaxOverlap)) {
									overlap_merge_candidates[nsTile][np1][np2] = true;
									overlap_merge_candidates[nsTile][np2][np1] = true;
									if (debugLevel > 1){
										System.out.println("overlapSameTileCandidates(): ACCEPTED pair nsTile="+nsTile+":"+np1+":"+np2+
												" as it has LOW enough overlap: "+
												" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
												" overlap2="+ ((int) (100 *overlaps[1]))+"% ");
									}
								} else {
									// maybe one of the planes is very weak and they are close by disparity?									
									// planes[nsTile][np1]
									if (areWeakSimilar(
											planes[nsTile][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
											planes[nsTile][np2], // TilePlanes.PlaneData plane2,
											Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
											"overlapSameTileCandidates() "+nsTile+":"+np1+":"+np2, //  String prefix,
											dl)// int debugLevel)
											){
										overlap_merge_candidates[nsTile][np1][np2] = true;
										overlap_merge_candidates[nsTile][np2][np1] = true;
										if (debugLevel > 0){
											System.out.println("overlapSameTileCandidates(): ACCEPTED pair nsTile="+nsTile+":"+np1+":"+np2+
													" even as it has HIGH overlap: "+
													" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
													" overlap2="+ ((int) (100 *overlaps[1]))+"% "+
													"because at least one plane is weak and they have small disparity difference");
										}
									// see if the planes are close enough to still be merged	
									} else if (min_distance > 0){
										boolean bad_plate1 = (planes[nsTile][np1].checkBadPlate(true) > plBadPlate);
										boolean bad_plate2 = (planes[nsTile][np2].checkBadPlate(true) > plBadPlate);
										if (bad_plate1 && bad_plate2) {
/*											
											overlap_merge_candidates[nsTile][np1][np2] = true;
											overlap_merge_candidates[nsTile][np2][np1] = true;
											if (debugLevel > 0){
												System.out.println("overlapSameTileCandidates(): ACCEPTED pair nsTile="+nsTile+":"+np1+":"+np2+
														" even as it has HIGH overlap: "+
														" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
														" overlap2="+ ((int) (100 *overlaps[1]))+"% "+
														"because both are \"bad plates\" (should it be?)");
											}
*/											
											if (debugLevel > 0){
												System.out.println("overlapSameTileCandidates(): REJECTED pair nsTile="+nsTile+":"+np1+":"+np2+
														" even as it has HIGH overlap: "+
														" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
														" overlap2="+ ((int) (100 *overlaps[1]))+"% "+
														"because both are \"bad plates\"");
											}

										} else {
											// from plane2 to 1 center 
											double wd2_1 = planes[nsTile][np1].getWorldPlaneDist2(planes[nsTile][np2]);
											double wd2_2 = planes[nsTile][np2].getWorldPlaneDist2(planes[nsTile][np1]);
											double wd2 = (!bad_plate1 && !bad_plate2) ? Math.min(wd2_1,wd2_2): (bad_plate1 ? wd2_1 : wd2_2);
//											wd2 = Math.sqrt(wd2); // just fro printing
											if (wd2 < min_distance2) {
												overlap_merge_candidates[nsTile][np1][np2] = true;
												overlap_merge_candidates[nsTile][np2][np1] = true;
												if (debugLevel > 0){
													System.out.println("overlapSameTileCandidates(): ACCEPTED pair nsTile="+nsTile+":"+np1+":"+np2+
															" even as it has HIGH overlap: "+
															" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
															" overlap2="+ ((int) (100 *overlaps[1]))+"% "+
															"because theay are close to each other: "+
															Math.sqrt(wd2) + " < " + min_distance + " meters.");
												}
											} else {
												if (debugLevel > 0){
													System.out.println("overlapSameTileCandidates(): REJECTED pair nsTile="+nsTile+":"+np1+":"+np2+
															" as it has HIGH overlap: "+
															" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
															" overlap2="+ ((int) (100 *overlaps[1]))+"% "+
															" and no other excuses - strong enough and the distance "+
															Math.sqrt(wd2) + " >= " + min_distance + " meters.");
												}
											}
										}
//										if ()
//										double d
										
									} else {
										if (debugLevel > 0){
											System.out.println("overlapSameTileCandidates(): REJECTED pair nsTile="+nsTile+":"+np1+":"+np2+
													" as it has HIGH overlap: "+
													" overlap1="+ ((int) (100 *overlaps[0]))+"% "+
													" overlap2="+ ((int) (100 *overlaps[1]))+"% ");
										}
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return overlap_merge_candidates;
	}
	
	
	// verify that after merging a pair composite plane  will have connections valid in each of the directions the pair had valid
	public  boolean [][][] keepSameTileConnections(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final boolean [][][]   valid_candidates, // will be updated
			final boolean    merge_low_eigen,
			final boolean useNonExcl, // consider only directions available for non-exclusive merges
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final double relax = 1.0;
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
//		final boolean [][] merge_pairs = new boolean [nStiles][];
//		final boolean [][][]  overlap_merge_candidates = new boolean [nStiles][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) if ( merge_candidates[nsTile0] != null) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > 0) && (nsTile0 == debug_stile)) ? 3: ((debugLevel > 1) ? 1:0);
						if (dl > 1){
							System.out.println("overlapSameTileCandidates(): nsTile="+nsTile0);
						}
//						int n_planes = planes[nsTile0].length;
						//						overlap_merge_candidates[nsTile] = new boolean [n_planes][n_planes];
						// get original directions
						for (int np1 = 0; np1 < planes[nsTile0].length; np1++) if (planes[nsTile0][np1] != null){
							for (int np2 = np1 + 1; np2 < planes[nsTile0].length; np2++) if (planes[nsTile0][np2] != null){
								if (valid_candidates[nsTile0][np1][np2]) { // only check pair considered valid
									boolean [] old_valid = new boolean[8];
									for (int dir = 0; dir < 8; dir++){
										if (planes[nsTile0][np1].hasMergedStrongValid(dir) || planes[nsTile0][np2].hasMergedStrongValid(dir)){
											if (!useNonExcl) {
												old_valid[dir] = true;
											} else {
												if ((planes[nsTile0][np1].getNonexclusiveStar() == null) ||
														(planes[nsTile0][np1].getNonexclusiveStar().getNeibBest(dir) >=0)){
													old_valid[dir] = true;
												}
												if ((planes[nsTile0][np1].getNonexclusiveStarEq() == null) ||
														(planes[nsTile0][np1].getNonexclusiveStarEq().getNeibBest(dir) >=0)){
													old_valid[dir] = true;
												}
												if ((planes[nsTile0][np2].getNonexclusiveStar() == null) ||
														(planes[nsTile0][np2].getNonexclusiveStar().getNeibBest(dir) >=0)){
													old_valid[dir] = true;
												}
												if ((planes[nsTile0][np2].getNonexclusiveStarEq() == null) ||
														(planes[nsTile0][np2].getNonexclusiveStarEq().getNeibBest(dir) >=0)){
													old_valid[dir] = true;
												}
											}
										}
//										old_valid[dir] = planes[nsTile0][np1].hasMergedStrongValid(dir) || planes[nsTile0][np2].hasMergedStrongValid(dir);
										//getNonexclusiveStar
									}
									// should be merged same way as later actually. Does it need to be recalculated from the original tiles?
									TilePlanes.PlaneData merged_pd =  planes[nsTile0][np1].mergePlaneToThis(
											planes[nsTile0][np2],      // PlaneData otherPd,
											1.0,         // double    scale_other,
											1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
											false,       // boolean   ignore_weights,
											true,        // boolean   sum_weights,
											plPreferDisparity, 
											debugLevel - 2); // int       debugLevel)
									// is the merge too bad already?
									double corr_max_eigen = corrMaxEigen(
											plMaxEigen,
											plDispNorm,
											merged_pd);
									if ((plMaxEigen != 0.0) &&
											(merged_pd.getValue() > corr_max_eigen)){
										valid_candidates[nsTile0][np1][np2] = false;
										valid_candidates[nsTile0][np2][np1] = false;
										if (debugLevel > 0){
											System.out.println("keepSameTileConnections(): REMOVING pair nsTile0="+nsTile0+":"+np1+":"+np2+
													" as the merge would have high eigenvalue = "+merged_pd.getValue()+" > " + corr_max_eigen);
										}
										continue; // to the next pair
									}

									// now verify that the merged plane can be connected in each of the original directions
									ArrayList<Integer> debug_dirs_list = new ArrayList<Integer>();
									for (int dir = 0; dir < 8; dir++) if (old_valid[dir]){ //
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										int nsTile = sty * stilesX + stx; // from where to get
										boolean fit = false;
										for (int np = 0; np < planes[nsTile].length; np++){
											if (planes[nsTile][np] != null){
												String prefix = "keepSameTileConnections() nsTile0="+nsTile0+":"+np1+":"+np2+" dir="+dir+" nsTile="+nsTile+" np="+np;
												TilePlanes.PlaneData other_plane = merged_pd.getPlaneToThis(
														planes[nsTile][np],
														dl-3); // debugLevel);

												if (planesFit(
														merged_pd,    // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
														other_plane,  // 			TilePlanes.PlaneData plane2,
														relax,        // double               relax = 1.0;
														merge_low_eigen, // false,                // boolean              merge_weak,   // use for same supertile merge 
														Double.NaN,         // double               merged_ev,    // if NaN will calculate assuming the same supertile
														Double.NaN,      //double               merged_ev_eq, // if NaN will calculate assuming the same supertile
														Double.NaN,        // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
														Double.NaN,     // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
														prefix,               // String prefix,
														dl-1)){                   // int debugLevel)
													fit = true;
													break;
												}
											}															
										}
										if (!fit){
											valid_candidates[nsTile0][np1][np2] = false;
											valid_candidates[nsTile0][np2][np1] = false;
											if (debugLevel > 0){
												debug_dirs_list.add(dir);
											}
											if (debugLevel < 2){
												break; // no need to check other directions, keep just for debug
											}
										}
									}
									
									if (debugLevel > 0){
										if (!debug_dirs_list.isEmpty()){
											System.out.println("keepSameTileConnections(): REMOVING pair nsTile0="+nsTile0+":"+np1+":"+np2+
													" as the merge would break previous connection in directions "+ debug_dirs_list);
										} else {
										System.out.println("keepSameTileConnections(): KEEPING pair nsTile0="+nsTile0+":"+np1+":"+np2+
												" as the merge would keep connection in each of the previously connected directions");
										}
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return valid_candidates;
	}
	
	// TODO: also evaluate background hiding behind hole in the with no edge?
	public  boolean[][] weakForeground(
			final TilePlanes.PlaneData [][] planes,
			final int [][][]                merge_candidates,
			final boolean [][][]            valid_candidates, // will be updated
//			final int                       stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
			final double                    min_fg_strength,
			final int                       outliers,  // remove any glitches?
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
//		final double [][][][][][]  merged_neib_ev = new double [nStiles][][][][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final boolean [][] weak_fgnd = new boolean [valid_candidates.length][];
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( merge_candidates[nsTile] != null) {
                        int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
						if (dl > 2){
							System.out.println("weakForeground(): nsTile="+nsTile);
						}
						//int n_planes = planes[nsTile].length;
						weak_fgnd[nsTile] = new boolean [merge_candidates[nsTile].length];
						for (int nPair = 0; nPair < merge_candidates[nsTile].length; nPair++) {
							int [] pair = merge_candidates[nsTile][nPair];
							if (valid_candidates[nsTile][pair[0]][pair[1]]){
								boolean weakfg = planes[nsTile][pair[0]].isWeakForeground(
										planes[nsTile][pair[1]], // final PlaneData fg_plane,
										//stMeasSel,               // final int stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
										outliers,                // final int outliers,  // remove any glitches?
										min_fg_strength,        // double min_fg_strength);
										dl);
								
								weakfg |= planes[nsTile][pair[1]].isWeakForeground(
										planes[nsTile][pair[0]], // final PlaneData fg_plane,
										//stMeasSel,               // final int stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
										outliers,                // final int outliers,  // remove any glitches?
										min_fg_strength,        // double min_fg_strength);
										dl);
								weak_fgnd[nsTile][nPair] = weakfg;
								if (debugLevel > 0) {
									if (weakfg && (debugLevel > 0)) {
										System.out.println(nsTile+":"+pair[0]+":"+pair[1]+" is a WEAK foreground pair");
									}else if (debugLevel > 1) {
										System.out.println(nsTile+":"+pair[0]+":"+pair[1]+" is NOT a weak foreground pair");
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return weak_fgnd;
	}
	
	/**
	 * Possible problem is that "normalizing" merge quality for low weights is not applicable for "star" plane that include neighbors
	 * Also calculates and outputs some of other costs - just for debugging
	 * TODO: Switch to a single "cost" function (costSameTileConnectionsAlt())
	 * @param planes per supertile, per plane - array of supertile instances
	 * @param merge_candidates array of merge pairs (per-supertile, per pair {start plane, end_plane} pairs) 
	 * @param valid_candidates per super tile , per first plane index, per second plane index boolean array of permitted to merge.
	 *                         valid_candidates array is updated as a result of this method                 
	 * @param relax - multiply thresholds by this value, so when relax > 1.0, the fitting requirements are loosened (used for
	 * conflicting planes)
	 * @param debugLevel debug level
	 * @param dbg_X tile x-index for detailed debug data
	 * @param dbg_Y tile y-index for detailed debug data
	 */
	public  void costSameTileConnections(
			final TilePlanes.PlaneData [][] planes,
			final int [][][]                merge_candidates,
			final boolean [][][]            valid_candidates, // will be updated
			final boolean [][]              weak_fg_pairs, 
			final double                    relax,
			final double                    relax_weak_fg,
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final double [][][][][][]  merged_neib_ev = new double [nStiles][][][][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) if ( merge_candidates[nsTile0] != null) {
//						int dl = ((debugLevel > 0) && (nsTile0 == debug_stile)) ? 3: ((debugLevel > 1) ? 2:0);
                        int dl = ((debugLevel > 1) && (nsTile0 == debug_stile)) ? 3: debugLevel;

						if (dl > 2){
							System.out.println("costSameTileConnections(): nsTile="+nsTile0);
						}
						int n_planes = planes[nsTile0].length;
						boolean [][] weak_fg = new boolean [n_planes][n_planes];
						if ((weak_fg_pairs != null) && (weak_fg_pairs[nsTile0] != null)) {
							for (int i = 0; i < merge_candidates[nsTile0].length; i++) {
								weak_fg[merge_candidates[nsTile0][i][0]][merge_candidates[nsTile0][i][1]] = true;
								weak_fg[merge_candidates[nsTile0][i][1]][merge_candidates[nsTile0][i][0]] = true;
							}
						}
						//						overlap_merge_candidates[nsTile] = new boolean [n_planes][n_planes];
						merged_neib_ev[nsTile0] = new double [n_planes][n_planes][4][][];
						// get original directions
						for (int np1 = 0; np1 < planes[nsTile0].length; np1++) if (planes[nsTile0][np1] != null){
							for (int np2 = np1 + 1; np2 < planes[nsTile0].length; np2++) if (planes[nsTile0][np2] != null){
								if (valid_candidates[nsTile0][np1][np2]) { // only check pair considered valid
									double var_relax = relax;
									if (weak_fg[np1][np2]){
										var_relax *= relax_weak_fg;
									}
									String prefix = "costSameTileConnections() fit weighted: nsTile0="+nsTile0+" np1="+np1+" np2="+np2;
									boolean fit1 = 	planesFit(
											planes[nsTile0][np1].getNonexclusiveStarFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
											planes[nsTile0][np2].getNonexclusiveStarFb(), // 			TilePlanes.PlaneData plane2,
											var_relax, //relax,      // double                    relax,
											true,       // boolean              merge_weak,   // use for same supertile merge 
											Double.NaN, // double merged_ev,    // if NaN will calculate assuming the same supertile
											Double.NaN, // double merged_ev_eq, // if NaN will calculate assuming the same supertile
											Double.NaN, // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
											Double.NaN, // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
											prefix, // String prefix,
											dl - 2); // int debugLevel)
									prefix = "costSameTileConnections() fit equal weight: nsTile0="+nsTile0+" np1="+np1+" np2="+np2;
									boolean fit2 = 	planesFit(
											planes[nsTile0][np1].getNonexclusiveStarEqFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
											planes[nsTile0][np2].getNonexclusiveStarEqFb(), // 			TilePlanes.PlaneData plane2,
											var_relax, // relax,      // double                    relax,
											true,       // boolean              merge_weak,   // use for same supertile merge 
											Double.NaN, // double merged_ev,    // if NaN will calculate assuming the same supertile
											Double.NaN, // double merged_ev_eq, // if NaN will calculate assuming the same supertile
											Double.NaN, // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
											Double.NaN, // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
											prefix, // String prefix,
											dl - 2); // int debugLevel)
									//										if (!fit1 || !fit2){
									if (!fit1 && !fit2){
										valid_candidates[nsTile0][np1][np2] = false;
										valid_candidates[nsTile0][np2][np1] = false;
										if (dl > 0){
											System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+
													(weak_fg[np1][np2]? " WEAK_FG_PAIR":"") +
													" REMOVING PAIR, fit1="+fit1+" fit2="+fit2);
										}

									} else {
										if (dl > 0){
											System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+
													(weak_fg[np1][np2]? " WEAK_FG_PAIR":"") +
													" KEEPING PAIR, fit1="+fit1+" fit2="+fit2);
										}

									}
									if (dl > 1){
										double [][] costs = new double[6][]; 
										costs[0] = getFitQualities(
												true, // en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs[1] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarEqFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarEqFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs[2] = getFitQualities(
												true, // en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2],          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs[3] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarEqFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2],          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs[4] = getFitQualities(
												true, // en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs[5] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarEqFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)

										System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs:");
										System.out.print("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted: ");
										for (int i = 0; i < costs[0].length;i++) {
											System.out.print(costs[0][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized: ");
										for (int i = 0; i < costs[1].length;i++) {
											System.out.print(costs[1][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted, 1-st star, 2-nd measured: ");
										for (int i = 0; i < costs[2].length;i++) {
											System.out.print(costs[2][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized, 1-st star, 2-nd measured: ");
										for (int i = 0; i < costs[3].length;i++) {
											System.out.print(costs[3][i]+" ");
										}
										System.out.println();
										System.out.println();
										System.out.print("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted, 1-st measured, 2-nd star: ");
										for (int i = 0; i < costs[4].length;i++) {
											System.out.print(costs[4][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized, 1-st measured, 2-nd star: ");
										for (int i = 0; i < costs[5].length;i++) {
											System.out.print(costs[5][i]+" ");
										}
										System.out.println();
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	public int [][][] filterPairsByConflicts(
			final TilePlanes.PlaneData [][] planes,
			final int [][][]                merge_candidates,
			final int [][][]                conflicts)
	{
		final int [][][]          filtered_pairs = new int [merge_candidates.length][][];
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if ((merge_candidates[nsTile] != null) && (conflicts[nsTile] != null)) {
							if (plConflSnglPair){
								int num_planes = 0;
								for (int i = 0; i < planes[nsTile].length; i++){
									if (planes[nsTile][i] != null) num_planes++;
								}
								if (num_planes > 2) {
									continue;
								}
							}
							if (plConflSngl && (conflicts[nsTile].length > 1)){
								continue;
							}

// TODO: maybe additionally check link costs correlation							
							HashSet<Point> pairs_set = new HashSet<Point>();
							for (int i = 0; i < merge_candidates[nsTile].length; i++ ){
								if (merge_candidates[nsTile][i] != null){
									Point p = new Point(merge_candidates[nsTile][i][0],merge_candidates[nsTile][i][1]);
									if (p.x > p.y){ // currently they are always ordered, just in case
										int tmp = p.x;
										p.x = p.y;
										p.y = tmp;
									}
									pairs_set.add(p);
								}
							}
							HashSet<Point> conflicts_set = new HashSet<Point>();
							for (int i = 0; i < conflicts[nsTile].length; i++ ){
								if (conflicts[nsTile][i] != null){
									Point p = new Point(conflicts[nsTile][i][0],conflicts[nsTile][i][1]);
									if (p.x > p.y){ // currently they are always ordered, just in case
										int tmp = p.x;
										p.x = p.y;
										p.y = tmp;
									}
									conflicts_set.add(p);
								}
							}
							pairs_set.retainAll(conflicts_set);
							if (!pairs_set.isEmpty()){
								filtered_pairs[nsTile] = new int [pairs_set.size()][2];
								int indx = 0;
								for (Point p: pairs_set){
									filtered_pairs[nsTile][indx][0] =   p.x;	
									filtered_pairs[nsTile][indx++][1] = p.y;	
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return filtered_pairs;
	}
	
	/**
	 * Calculate costs of merging planes of the same supertile and remove those that are exceed threshold
	 * Also calculates and outputs some of other costs - just for debugging
	 * @param threshold cost threshold for planes that do have valid non-exclusive neighbors that help to determine
	 * plane orientation
	 * @param threshold_nostar cost threshold for the planes that do not have valid neighbors
	 * @param planes per supertile, per plane - array of supertile instances
	 * @param merge_candidates array of merge pairs (per-supertile, per pair {start plane, end_plane} pairs) 
	 * @param valid_candidates per super tile , per first plane index, per second plane index boolean array of permitted to merge.
	 *                         valid_candidates array is updated as a result of this method                 
	 * @param debugLevel debug level
	 * @param dbg_X tile x-index for detailed debug data
	 * @param dbg_Y tile y-index for detailed debug data
	 */
	public  void costSameTileConnectionsAlt(
			final double                    threshold,
			final double                    threshold_nostar,
			final TilePlanes.PlaneData [][] planes,
			final int [][][]                merge_candidates,
			final boolean [][][]            valid_candidates, // will be updated
			final boolean [][]              weak_fg_pairs,
			final double                    relax_weak_fg,
			final int                       debugLevel,
			final int                       dbg_X,
			final int                       dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final double [][][][][][]  merged_neib_ev = new double [nStiles][][][][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) if ( merge_candidates[nsTile0] != null) {
//						int dl = ((debugLevel > 0) && (nsTile0 == debug_stile)) ? 3: ((debugLevel > 1) ? 2:0);
                        int dl = ((debugLevel > 1) && (nsTile0 == debug_stile)) ? 3: debugLevel;
						
						if (dl > 2){
							System.out.println("costSameTileConnectionsAlt(): nsTile="+nsTile0);
						}
						int n_planes = planes[nsTile0].length;
						boolean [][] weak_fg = new boolean [n_planes][n_planes];
						if ((weak_fg_pairs != null) && (weak_fg_pairs[nsTile0] != null)) {
							for (int i = 0; i < merge_candidates[nsTile0].length; i++) {
								weak_fg[merge_candidates[nsTile0][i][0]][merge_candidates[nsTile0][i][1]] = true;
								weak_fg[merge_candidates[nsTile0][i][1]][merge_candidates[nsTile0][i][0]] = true;
							}
						}
						merged_neib_ev[nsTile0] = new double [n_planes][n_planes][4][][];
						// get original directions
						for (int np1 = 0; np1 < planes[nsTile0].length; np1++) if (planes[nsTile0][np1] != null){
							for (int np2 = np1 + 1; np2 < planes[nsTile0].length; np2++) if (planes[nsTile0][np2] != null){
								if ((valid_candidates[nsTile0] == null) || (valid_candidates[nsTile0][np1] == null)){
									System.out.println("BUG in costSameTileConnectionsAlt()");
								}
								if (valid_candidates[nsTile0][np1][np2]) { // only check pair considered valid
									String prefix = "costSameTileConnectionsAlt() fit weighted: nsTile0="+nsTile0+" np1="+np1+" np2="+np2;
									double [] costs = new double[2]; 
									costs[0] = getLinkCost(
											true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
											planes[nsTile0][np1].getNonexclusiveStarFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
											planes[nsTile0][np2].getNonexclusiveStarFb(),          // TilePlanes.PlaneData plane2,
											Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
											Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
											Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
											Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
											prefix,  // String prefix,
											dl - 1); // int debugLevel)
									prefix = "costSameTileConnectionsAlt() fit equal weight: nsTile0="+nsTile0+" np1="+np1+" np2="+np2;
									costs[1] = getLinkCost(
											true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
											planes[nsTile0][np1].getNonexclusiveStarEqFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
											planes[nsTile0][np2].getNonexclusiveStarEqFb(),          // TilePlanes.PlaneData plane2,
											Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
											Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
											Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
											Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
											prefix,  // String prefix,
											dl - 1); // int debugLevel)
									boolean star1 = (planes[nsTile0][np1].getNonexclusiveStar() != null) && (planes[nsTile0][np2].getNonexclusiveStar() != null);
									boolean star2 = (planes[nsTile0][np1].getNonexclusiveStarEq() != null) && (planes[nsTile0][np2].getNonexclusiveStarEq() != null);
									double ts = threshold;
									double tns = threshold_nostar;
									if (weak_fg[np1][np2]){
										ts *= relax_weak_fg;
										tns *= relax_weak_fg;
									}
//									boolean fit1 = 	costs[0] < (star1 ? threshold : threshold_nostar);
//									boolean fit2 = 	costs[1] < (star2 ? threshold : threshold_nostar);
									
									boolean fit1 = 	costs[0] < (star1 ? ts : tns);
									boolean fit2 = 	costs[1] < (star2 ? ts : tns);
									
//									if (!fit1 || !fit2){
									if (!fit1 && !fit2){
										valid_candidates[nsTile0][np1][np2] = false;
										valid_candidates[nsTile0][np2][np1] = false;
										if (dl > 0){
											System.out.println("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+
													(weak_fg[np1][np2]? " WEAK_FG_PAIR":"") +
													" REMOVING PAIR, fit1="+fit1+" fit2="+fit2+ " (star1="+star1+", star2="+star2+")");
										}
									} else {
										if (dl > 0){
											System.out.println("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+
													(weak_fg[np1][np2]? " WEAK_FG_PAIR":"") +
													" KEEPING PAIR, fit1="+fit1+" fit2="+fit2+ " (star1="+star1+", star2="+star2+")");
										}
									}
									if (dl > 1){
										double [][] costs0 = new double[6][]; 
										costs0[0] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // int debugLevel)
										costs0[1] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarEqFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarEqFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // int debugLevel)
										
										costs0[2] = getFitQualities(
												true, // en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2],          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs0[3] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1].getNonexclusiveStarEqFb(), // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2],          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs0[4] = getFitQualities(
												true, // en_sticks, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)
										costs0[5] = getFitQualities(
												true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
												planes[nsTile0][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
												planes[nsTile0][np2].getNonexclusiveStarEqFb(),          // TilePlanes.PlaneData plane2,
												Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
												Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
												Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
												nsTile0+":"+np1+":"+np2,  // String prefix,
												dl); // 0); // int debugLevel)

										System.out.println("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs:");
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted: ");
										for (int i = 0; i < costs0[0].length;i++) {
											System.out.print(costs0[0][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized: ");
										for (int i = 0; i < costs0[1].length;i++) {
											System.out.print(costs0[1][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted, 1-st star, 2-nd measured: ");
										for (int i = 0; i < costs0[2].length;i++) {
											System.out.print(costs0[2][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized, 1-st star, 2-nd measured: ");
										for (int i = 0; i < costs0[3].length;i++) {
											System.out.print(costs0[3][i]+" ");
										}
										System.out.println();
										System.out.println();
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted, 1-st measured, 2-nd star: ");
										for (int i = 0; i < costs0[4].length;i++) {
											System.out.print(costs0[4][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized, 1-st measured, 2-nd star: ");
										for (int i = 0; i < costs0[5].length;i++) {
											System.out.print(costs0[5][i]+" ");
										}
										System.out.println();
										
										
										
										
										System.out.println("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs:");
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs weighted: ");
										for (int i = 0; i < costs0[0].length;i++) {
											System.out.print(costs0[0][i]+" ");
										}
										System.out.println();
										System.out.print("costSameTileConnectionsAlt(): nsTile="+nsTile0+":"+np1+":"+np2+" costs equalized: ");
										for (int i = 0; i < costs0[1].length;i++) {
											System.out.print(costs0[1][i]+" ");
										}
										System.out.println();
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public  double  [][][][][][] costSameTileConnectionsOld(
			final boolean ignore_weights,
			final double threshold_worst,
			final double threshold_world_worst,
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final boolean [][][]   valid_candidates, // will be updated
//			final boolean    merge_low_eigen,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
//		final boolean [][] merge_pairs = new boolean [nStiles][];
		final double [][][][][][]  merged_neib_ev = new double [nStiles][][][][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) if ( merge_candidates[nsTile0] != null) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > 0) && (nsTile0 == debug_stile)) ? 3: ((debugLevel > 1) ? 1:0);
						if (dl > 1){
							System.out.println("costSameTileConnections(): nsTile="+nsTile0);
						}
						int n_planes = planes[nsTile0].length;
						//						overlap_merge_candidates[nsTile] = new boolean [n_planes][n_planes];
						merged_neib_ev[nsTile0] = new double [n_planes][n_planes][4][][];
						// get original directions
						for (int np1 = 0; np1 < planes[nsTile0].length; np1++) if (planes[nsTile0][np1] != null){
							for (int np2 = np1 + 1; np2 < planes[nsTile0].length; np2++) if (planes[nsTile0][np2] != null){
								if (valid_candidates[nsTile0][np1][np2]) { // only check pair considered valid
									boolean [] old_valid = new boolean[8];
									for (int dir = 0; dir < 8; dir++){
										old_valid[dir] = planes[nsTile0][np1].hasMergedValid(dir) || planes[nsTile0][np2].hasMergedValid(dir);
									}
									// should be merged same way as later actually. Does it need to be recalculated from the original tiles?
									TilePlanes.PlaneData this_plane =  planes[nsTile0][np1].mergePlaneToThis(
											planes[nsTile0][np2],      // PlaneData otherPd,
											1.0,            // double    scale_other,
											1.0,            // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
											ignore_weights, // boolean   ignore_weights,
											true,           // boolean   sum_weights,
											plPreferDisparity, 
											dl - 3); // int       debugLevel)
									// is the merge too bad already - should be already tested in keepSameTileConnections()
									// now for each of the valid directions calculate similar to  matchPlanes(), but in all 8 directions
									double [][] merged_ev =     new double [8][];
									double [][] merged_ev_eq =  new double [8][];
									double [][] merged_wev =    new double [8][];
									double [][] merged_wev_eq = new double [8][];

									
									for (int dir = 0; dir < 8; dir++) if (old_valid[dir]){ // all 8 here
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										//											if ((sty < stilesY) && (sty > 0) && (stx < 0)) {
										if ((stx < stilesX) && (sty < stilesY) && (sty > 0)) {
											int nsTile = sty * stilesX + stx; // from where to get
											if (nsTile >= planes.length){
												System.out.println("BUG!!!!");
											} else {
												TilePlanes.PlaneData [] other_planes = planes[nsTile];
												if (other_planes != null) {
													merged_ev[dir] =     new double [other_planes.length];
													merged_ev_eq[dir] =  new double [other_planes.length];
													merged_wev[dir] =    new double [other_planes.length];
													merged_wev_eq[dir] = new double [other_planes.length];
//													this_plane.initMergedValue(dir,other_planes.length); // filled with NaN
													for (int np = 0; np < other_planes.length; np ++){
														if (other_planes[np] != null) {
															TilePlanes.PlaneData other_plane = this_plane.getPlaneToThis(
																	other_planes[np],
																	dl - 3); // debugLevel);
															if (other_plane !=null) { // now always, but may add later
																TilePlanes.PlaneData merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
																		false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		plPreferDisparity, 
																		dl - 3); // int       debugLevel)

																if (merged_pd !=null) { // now always, but may add later
																	merged_ev[dir][np] =  merged_pd.getValue(); // smallest eigenValue
																	merged_wev[dir][np] = merged_pd.getWValue(); // smallest eigenValue
																	if (Double.isNaN(merged_ev[dir][np]) || Double.isNaN(merged_wev[dir][np]) ){
																		System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" NaN");
																	}
																}
																
																merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
																		true, // false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		plPreferDisparity, 
																		dl - 3); // int       debugLevel)

																if (merged_pd !=null) { // now always, but may add later
																	merged_ev_eq[dir][np] =  merged_pd.getValue(); // smallest eigenValue
																	merged_wev_eq[dir][np] = merged_pd.getWValue(); // smallest eigenValue
																	if (Double.isNaN(merged_ev_eq[dir][np]) || Double.isNaN(merged_wev_eq[dir][np]) ){
																		System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" NaN2");
																	}
																}
															}
														}
													}
												}
											}
										}
									}
									merged_neib_ev[nsTile0][np1][np2][0] = merged_ev;
									merged_neib_ev[nsTile0][np1][np2][1] = merged_ev_eq;
									merged_neib_ev[nsTile0][np1][np2][2] = merged_wev;
									merged_neib_ev[nsTile0][np1][np2][3] = merged_wev_eq;
// calculate here, later move to a separate method
									double [] weighted_costs= new double[8]; // first plane with its connections * startWeight, second, then composite with same neibs
									double sw1 = 0.0, sw2 = 0.0;
									for (int dir = 0; dir < 8; dir++) { // all 8 here
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										if ((stx < stilesX) && (sty < stilesY) && (sty > 0)) {
											int nsTile = sty * stilesX + stx; // from where to get
											int nnp1 =  planes[nsTile0][np1].getNeibBest(dir);
											if (nnp1 >= 0){
												double sw = planes[nsTile][nnp1].getStarValueWeight()[1];
												weighted_costs[0] += sw * planes[nsTile0][np1].getMergedValueEq(dir, nnp1);
												weighted_costs[4] += sw * merged_ev_eq[dir][nnp1];
												weighted_costs[2] += sw * planes[nsTile0][np1].getMergedWValueEq(dir, nnp1);
												weighted_costs[6] += sw * merged_wev_eq[dir][nnp1];
												sw1 += sw;
												if (Double.isNaN(weighted_costs[0]) || Double.isNaN(weighted_costs[2]) ){
													System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" NaN3");
												}
											}

											int nnp2 =  planes[nsTile0][np2].getNeibBest(dir);
											if (nnp2 >= 0){
												double sw = planes[nsTile][nnp2].getStarValueWeight()[1];
												weighted_costs[1] += sw * planes[nsTile0][np2].getMergedValueEq(dir, nnp2);
												weighted_costs[5] += sw * merged_ev_eq[dir][nnp2];
												weighted_costs[3] += sw * planes[nsTile0][np2].getMergedValueEq(dir, nnp2);
												weighted_costs[7] += sw * merged_ev_eq[dir][nnp2];
												sw2 += sw;
												if (Double.isNaN(weighted_costs[1]) || Double.isNaN(weighted_costs[3]) ){
													System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+" NaN3");
												}
											}
										}
									}
									if ((sw1 > 0.0) && (sw2 > 0.0)) {
										weighted_costs[0] /= sw1;
										weighted_costs[2] /= sw1;
										weighted_costs[4] /= sw1;
										weighted_costs[6] /= sw1;
										weighted_costs[1] /= sw2;
										weighted_costs[3] /= sw2;
										weighted_costs[5] /= sw2;
										weighted_costs[7] /= sw2;
										double k1 = weighted_costs[4]/weighted_costs[0];
										double k2 = weighted_costs[5]/weighted_costs[1];
										double k1w = weighted_costs[6]/weighted_costs[2];
										double k2w = weighted_costs[7]/weighted_costs[3];
										double worst_k = (k1 > k2)? k1: k2;
										double worst_kw = (k1w > k2w)? k1w: k2w;
										if (dl > -1){
											System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+
													" worst_k="+worst_k+", sum="+(k1+k2)+", k1 = "+k1+", k2 = "+k2+
													" weighted costs = ["+weighted_costs[0]+", "+weighted_costs[1]+", "+weighted_costs[4]+", "+weighted_costs[5]+"]");
											System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+
													" worst_kw="+worst_kw+", sum="+(k1w+k2w)+", k1w = "+k1w+", k2w = "+k2w+
													" weighted costs = ["+weighted_costs[2]+", "+weighted_costs[3]+", "+weighted_costs[6]+", "+weighted_costs[7]+"]");
										}
										if ((worst_k > threshold_worst) || (worst_kw > threshold_world_worst)){
											valid_candidates[nsTile0][np1][np2] = false;
											valid_candidates[nsTile0][np2][np1] = false;
											if (dl > -1){
												System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+
														" REMOVING PAIR");
											}
											
										}
										
///		final double threshold_worst,
//		final double threshold_world_worst,
										
									} else {
										if (dl > -1){
											System.out.println("costSameTileConnections(): nsTile="+nsTile0+":"+np1+":"+np2+
													" one of the tiles was not connected");
										}
										
									}
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return merged_neib_ev;
	}
	
	
	public boolean [][] mergeSameTileEvaluate(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final boolean [][][] plane_nooverlaps,
			final double                    relax,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final boolean [][] merge_pairs = new boolean [nStiles][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( merge_candidates[nsTile] != null) {
						merge_pairs[nsTile] = new boolean [merge_candidates[nsTile].length];
						int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? 2: ((debugLevel > 1) ? 1:0);
						if (dl > 0){
							System.out.println("mergeSameTileEvaluate(): nsTile="+nsTile);
						}
						
						for (int pair = 0; pair < merge_candidates[nsTile].length; pair ++){
							int np1 =  merge_candidates[nsTile][pair][0];
							int np2 =  merge_candidates[nsTile][pair][1];
							if ((plane_nooverlaps == null) || (plane_nooverlaps[nsTile] == null) || plane_nooverlaps[nsTile][np1][np2]) {
								String prefix = "mergeSameTileEvaluate() pair="+pair+" nsTile="+nsTile+" np1="+np1+" np2="+np2;
								if (planesFit(
										planes[nsTile][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
										planes[nsTile][np2], // 			TilePlanes.PlaneData plane2,
										relax,      // double                    relax,
										true,       // boolean              merge_weak,   // use for same supertile merge 
										Double.NaN, // double merged_ev,    // if NaN will calculate assuming the same supertile
										Double.NaN, // double merged_ev_eq, // if NaN will calculate assuming the same supertile
										Double.NaN, // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
										Double.NaN, // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
										prefix, // String prefix,
										dl) // int debugLevel)
										){
									merge_pairs[nsTile][pair] = true;
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return merge_pairs;
	}
	/**
	 * Extract non-conflicting groups of the same supertile planes that are OK to merge. For multiple intersecting pairs
	 * this method starts for from the lowest cost merging and adds planes (from merge_candidates) until they do not cause
	 * conflicts with already selected ones 
	 * @param merge_candidates array of merge pairs (per-supertile, per pair {start plane, end_plane} pairs) 
	 * @param valid_candidates per super tile , per first plane index, per second plane index boolean array of permitted to merge.
	 *                         valid_candidates array is updated as a result of this method                 
	 * @param relax - multiply thresholds by this value, so when relax > 1.0, the fitting requirements are loosened (used for
	 * conflicting planes)
	 * @param debugLevel debug level
	 * @param dbg_X tile x-index for detailed debug data
	 * @param dbg_Y tile y-index for detailed debug data
	 * @return array of per supertile group sets for merging. Each group set is an array of groups. Each group is an array
	 * of plane indices. 
	 */
	// TODO: get rid of quality_index use getLinkCost()
	public int [][][] extractMergeSameTileGroups(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final boolean [][][] valid_candidates, // will be updated
//			final boolean [][][] plane_nooverlaps,
			final double                    relax,

			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY;
		final int [][][] merge_groups = new int [nStiles][][];
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int quality_index = 3; // 0 - using strengths, 1 - equal strengths, 2 - composite
		
// TODO Make nooverlaps be overriden if ne of the planes is very weak and they are close by disparity		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( merge_candidates[nsTile] != null) {
						if (planes[nsTile] == null){
							System.out.println("extractMergeSameTileGroups() planes["+nsTile+"] = null");
							continue;
						}
						boolean [][] merge_pairs = new boolean [planes[nsTile].length][planes[nsTile].length]; 
//						int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 2: ((debugLevel > 1) ? 1:0);
						int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
						if (dl > 1){
							System.out.println("extractMergeSameTileGroups(): nsTile="+nsTile);
						}
						HashSet<Integer> yet_to_merge = new HashSet<Integer>(); 
						for (int pair = 0; pair < merge_candidates[nsTile].length; pair ++){
							int np1 =  merge_candidates[nsTile][pair][0];
							int np2 =  merge_candidates[nsTile][pair][1];
							if ((valid_candidates == null) || (valid_candidates[nsTile] == null) || valid_candidates[nsTile][np1][np2]) {
								yet_to_merge.add(np1);
								yet_to_merge.add(np2);
								String prefix = "extractMergeSameTileGroups() pair="+pair+" nsTile="+nsTile+" np1="+np1+" np2="+np2;
								if (planesFit(
										planes[nsTile][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
										planes[nsTile][np2], // 			TilePlanes.PlaneData plane2,
										relax,      // double                    relax,
										true,       // boolean              merge_weak,   // use for same supertile merge 
										Double.NaN, // double merged_ev,    // if NaN will calculate assuming the same supertile
										Double.NaN, // double merged_ev_eq, // if NaN will calculate assuming the same supertile
										Double.NaN, // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
										Double.NaN, // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
										prefix, // String prefix,
										dl-1) // int debugLevel)
										){
									merge_pairs[np1][np2] = true; 
									merge_pairs[np2][np1] = true; 
								}
							}
						}
						if (!yet_to_merge.isEmpty()){
							ArrayList<HashSet<Integer>> groups_list = new ArrayList<HashSet<Integer>>();
							HashSet<Integer> rejected = new HashSet<Integer>();
							while (!yet_to_merge.isEmpty()){
								//- see if there are simple sets - all connected to each other - they do not need any recalculation of the merge fitness
								// extract group of somehow connected to each other, disregarding all conflicts
								HashSet<Integer> group = new HashSet<Integer>();
								for (Integer el: yet_to_merge){
									if (group.isEmpty()) {
										group.add(el);
									} else {
										// see if it is connected to any in the group
										for (Integer el1: group){
											if (merge_pairs[el][el1]){
												group.add(el);
												break;
											}
										} // nothing special to do if not connected
									}
								}
								// either add this group if it is all connected with no conflicts, are subdivide
								// in any case remove for the yet_to_merge
								yet_to_merge.removeAll(group);
								// verify that group has more than one element and they are all connected to each other
								if (group.size() == 1){
									rejected.addAll(group); // is it needed?
									if (debugLevel>1) {
										System.out.println("extractMergeSameTileGroups() nsTile = "+nsTile+" : single,  not connected to others: "+group); 
									}
									continue; // searching for the next group 
								} else {
									
									boolean all_connected = true;
									for (Integer el1: group){
										for (Integer el2: group){
											if ((el2 > el1) && !merge_pairs[el2][el1]){
												all_connected = false;
												break;
											}
										}
										if (!all_connected){
											break;
										}
									}
									if (all_connected){
										if (debugLevel>1) {
											System.out.println("extractMergeSameTileGroups() nsTile = "+nsTile+" : found interconnected group: "+group); 
										}
										groups_list.add(group);
										continue;
									} else {
										if (debugLevel>1) {
											System.out.println("extractMergeSameTileGroups() nsTile = "+nsTile+" : found incomplete group: "+group); 
										}
										while (group.size() > 1){ // maybe several sub-groups can be generated
											// here resolve that incomplete group:
											// First - find the
											Point best_pair = null;
											double best_quality = Double.NaN; // lower - better
											
											// create a seed - merge first pair (with strongest affiliation)
											for (Integer np1: group){
												for (Integer np2: group){
													if ((np2 > np1) && merge_pairs[np1][np2]){
														String prefix = "extractMergeSameTileGroups() nsTile="+nsTile+" np1="+np1+" np2="+np2;
														double [] qualities = getFitQualities( // {this_rq, this_rq_eq};
																true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
																planes[nsTile][np1], //TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																planes[nsTile][np2], //TilePlanes.PlaneData plane2,
																Double.NaN, // double merged_ev,    // if NaN will calculate assuming the same supertile
																Double.NaN, // double merged_ev_eq, // if NaN will calculate assuming the same supertile
																Double.NaN, // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
																Double.NaN, // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
																prefix, // String prefix,
																dl-1); // int debugLevel)
														if (qualities != null) {
															double this_rq = qualities[quality_index];
															if ((best_pair == null) || (best_quality > this_rq)){
																best_quality = this_rq;
																best_pair = new Point(np1, np2); // np2 > np1
															}
														}
													}
												}
											}
											// best_pair == null - can it be? What to do?
											if (best_pair != null){
												// merge np1, np2, add
												TilePlanes.PlaneData merged_pd = planes[nsTile][best_pair.x].mergePlaneToThis(
														planes[nsTile][best_pair.y], // PlaneData otherPd,
														1.0,         // double    scale_other,
														1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
														false,       // boolean   ignore_weights,
														true, // boolean   sum_weights,
														plPreferDisparity, 
														debugLevel - 2); // int       debugLevel)
												//													merged_ev = merged_pd.getValue();
												HashSet<Integer> sub_group =    new HashSet<Integer>();
												sub_group.add(best_pair.x);
												sub_group.add(best_pair.y);
												group.removeAll(sub_group);
												// now grow group while possible and no conflicts with existing members
												label_grow:
												{
													while (!group.isEmpty()) {
														Integer bestMatch = null;
														best_quality = Double.NaN;
														for (Integer np: group){
															// make sure it does not overlap with any of existing
															boolean nooverlap = true;
															for (Integer np1: sub_group){
																if ((valid_candidates != null) &&
																		(valid_candidates[nsTile] != null) &&
																		!valid_candidates[nsTile][np][np1]) {
																	nooverlap = false;
																	break;
																}
															}
															if (nooverlap){
																String prefix = "extractMergeSameTileGroups().2 nsTile="+nsTile+" np="+np;
																double [] qualities = getFitQualities( // {this_rq, this_rq_eq}; 
																		true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
																		merged_pd, //TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																		planes[nsTile][np], //TilePlanes.PlaneData plane2,
																		Double.NaN, // double merged_ev,    // if NaN will calculate assuming the same supertile
																		Double.NaN, // double merged_ev_eq, // if NaN will calculate assuming the same supertile
																		Double.NaN, // double merged_wev,    // if NaN will calculate assuming the same supertile - for world
																		Double.NaN, // double merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
																		prefix, // String prefix,
																		dl-1); // int debugLevel)
																if (qualities != null) {
																	double this_rq = qualities[quality_index];
																	if ((bestMatch == null) || (best_quality > this_rq)){
																		best_quality = this_rq;
																		bestMatch = np;
																	}
																}
															}
														}
														if (bestMatch != null){
															sub_group.add(bestMatch);
															group.remove(bestMatch);
														} else {
															break label_grow;
														}
													}
												}
												if (debugLevel>1) {
													System.out.println("extractMergeSameTileGroups() nsTile = "+nsTile+" : extracted sub-group: "+sub_group); 
												}
												groups_list.add(sub_group);
											} else {
												System.out.println("===== All remaining planes are incompatible to each other: extractMergeSameTileGroups() nsTile = "+nsTile+" : found incomplete group: "+group+""
														+ " best_pair = null");
												break; // all remaining planes are incompatible to each other 

											}
										}

									}
									// see if all are connected
								}
							}
							if (!groups_list.isEmpty()){
								//		final int [][][] merge_groups = new int [nStiles][][];
								merge_groups[nsTile] = new int [groups_list.size()][];
								int ng = 0;
								//								for (int ng = 0; ng < groups_list.size(); ng++) {
								for (HashSet<Integer> group: groups_list){
									merge_groups[nsTile][ng] = new int[group.size()];
									// sort to ascending plane indices 
									ArrayList<Integer> sorted_group = new ArrayList<Integer>(group);
									Collections.sort(sorted_group);
									int indx = 0;
									for (Integer np:sorted_group){
										merge_groups[nsTile][ng][indx++] = np;
									}
									ng++;
								}
								//							}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		if (debugLevel > 0){
			System.out.println("extractMergeSameTileGroups():");
			for (int nsTile = 0; nsTile < merge_groups.length; nsTile++) if (merge_groups[nsTile] != null){
				for (int ng = 0; ng < merge_groups[nsTile].length; ng++){
					System.out.print("nsTile="+nsTile+" ("+(nsTile % stilesX)+":"+(nsTile / stilesX)+"): [");
					for (int i = 0; i < merge_groups[nsTile][ng].length; i++){
						if (i > 0) System.out.print(", "); 
						System.out.print(merge_groups[nsTile][ng][i]);
					}
					System.out.println("]");
				}
			}
		}
		return merge_groups;
	}
	
	
	public double corrMaxEigen(
			double maxEigen,
			double dispNorm,
			TilePlanes.PlaneData pd)
	{
		return corrMaxEigen(
				maxEigen,
				dispNorm,
				pd.getZxy()[0]);
	}

	public double corrMaxEigen(
			double maxEigen,
			double dispNorm,
			double disparity)
	{
		double corrV = maxEigen;
		if ((dispNorm > 0.0) && (disparity > dispNorm)) {
			//					double dd = (dispNorm + z0)/ dispNorm; // > 1
			double dd = disparity/ dispNorm; // > 1
			corrV *= dd * dd; // > original
		}
		return corrV;

	}
	
	/**
	 * Measure how merging of two plains degrades individual flatness (smaller - better). For comparing divide by (w1+w2) to make strong
	 * planes score better  			
	 * @param L1 smallest eigenvalue of the first plane
	 * @param L2 smallest eigenvalue of the second plane
	 * @param L  smallest eigenvalue of the merged plane
	 * @param w1 weight of the first plane
	 * @param w2 weight of the second plane
	 * @param eigen_floor add to each L
	 * @return degrading by merging measure. 0 if both are co-planar, is supposed to be positive. very "bad" planes do produce negative results - 
	 * not yet clear why (related to non-linear coordinate transformation?) 
	 */
	
	public double mergeRQuality(
			double L1_in,
			double L2_in,
			double L_in,
			double w1,
			double w2,
			double eigen_floor)
	{
		double minVal = 0.00001; // Sometimes get small negative (precision) when merging 2 equal planes 
		double L1 = L1_in + eigen_floor;
		double L2 = L2_in + eigen_floor;
		double L =  L_in + eigen_floor;
		//				double Lav = Math.sqrt((L1*L1*w1 + L2*L2*w2)/(w1+w2));
		double Lav = (L1*w1 + L2*w2)/(w1+w2);
		///				double wors = (L - Lav)*(w1+w2)*(w1+w2) /(Lav*w1*w2);
		///				double rquality = (L - Lav)*(w1+w2) /(Lav*w1*w2); // ==wors/(w1+w2) to amplify stronger planes
		double rquality =  (L - Lav)*(w1+w2)*(w1+w2) /(Lav*w1*w2); 
		
		if (rquality < minVal){
//			System.out.println("Precision limitation for the same tile merging: mergeRQuality("+L1_in+", "+L2_in+", "+L_in+", "+w1+", "+w2+", "+eigen_floor+") -> "+rquality);
			rquality=minVal;
		}
		return rquality;
	}
	
	public void calcStarValueStrength(
			final boolean        set_start_planes,
			final double         orthoWeight,
			final double         diagonalWeight,
			final double         starPwr,    // Divide cost by number of connections to this power
			final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
			final double         weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			final double         starValPwr, //  Raise value of each tile before averaging
			final int            steps,
			final TilePlanes.PlaneData [][] planes,
			final boolean        preferDisparity,
			final int            debugLevel)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final int nStiles =       stilesX * stilesY; 
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					ConnectionCosts connectionCosts = new ConnectionCosts(
							orthoWeight,    // double         orthoWeight,
							diagonalWeight, // double         diagonalWeight,
							starPwr,        // double         starPwr,    // Divide cost by number of connections to this power
							starWeightPwr,  // double         starWeightPwr,    // Use this power of tile weight when calculating connection cost
							weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
							starValPwr,     //double          starValPwr, //  Raise value of each tile before averaging
							steps,          // int            steps,
							planes,         // TilePlanes.PlaneData [][] planes,
							tnSurface,      // TileNeibs tnSurface,
							preferDisparity); // boolean preferDisparity)
					int [] mod_supertiles = new int[1];
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) {
						if ( planes[nsTile] != null) {
							mod_supertiles[0] = nsTile;
							connectionCosts.initConnectionCosts(
									set_start_planes,									
									mod_supertiles,
									debugLevel);
							double [][][] val_weights = connectionCosts.getValWeights();
							for (int np = 0; np < planes[nsTile].length; np++){ // nu
								if (planes[nsTile][np] != null) {
									planes[nsTile][np].setStarValueWeight(val_weights[0][np]);
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}

	public void updateStarValueStrength(
			final int []         mod_supertiles,
			final double         orthoWeight,
			final double         diagonalWeight,
			final double         starPwr,    // Divide cost by number of connections to this power
			final double         starWeightPwr,    // Use this power of tile weight when calculating connection cost			
			final double         weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
			final double         starValPwr, //  Raise value of each tile before averaging
			final int            steps,
			final TilePlanes.PlaneData [][] planes,
			final boolean        preferDisparity,
			final int            debugLevel)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					ConnectionCosts connectionCosts = new ConnectionCosts(
							orthoWeight,    // double         orthoWeight,
							diagonalWeight, // double         diagonalWeight,
							starPwr,        // double         starPwr,    // Divide cost by number of connections to this power
							starWeightPwr,    // Use this power of tile weight when calculating connection cost							
							weightToDens,    // Balance weighted density against density. 0.0 - density, 1.0 - weighted density
							starValPwr,      //double     starValPwr, //  Raise value of each tile before averaging
							steps,          // int            steps,
							planes,         // TilePlanes.PlaneData [][] planes,
							tnSurface,      // TileNeibs tnSurface,
							preferDisparity); // boolean preferDisparity)
					int [] supertiles = new int[1];
					for (int isTile = ai.getAndIncrement(); isTile < mod_supertiles.length; isTile = ai.getAndIncrement()) {
						int nsTile = mod_supertiles[isTile];
						if ((nsTile >= 0) && ( planes[nsTile] != null)) {
							supertiles[0] = nsTile;
							connectionCosts.initConnectionCosts(supertiles, debugLevel - 2);
							double [][][] val_weights = connectionCosts.getValWeights();
							for (int np = 0; np < planes[nsTile].length; np++){ // nu
								if (planes[nsTile][np] != null) {
									planes[nsTile][np].setStarValueWeight(val_weights[0][np]);
								}
							}
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
	}
	
	
	// result is needed for debug only
	public TilePlanes.PlaneData [][][] conditionSuperTiles(
			final TilePlanes.PlaneData [][] planes,
			final int   max_num_merge_try, 
			final int         debugLevel)
	{
		// try to merge multiple times
		TilePlanes.PlaneData [][][] dbg_orig_planes = null;
		if (debugLevel > 0) {
			dbg_orig_planes = new TilePlanes.PlaneData [max_num_merge_try][][];
		}
		
		for (int num_merge_try = 0; num_merge_try < max_num_merge_try; num_merge_try++){
			// Calculate what be thickness eigenvalues for merging this plane with individual neighbors
			matchPlanes(
					planes, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			// Create matrix of connection costs between layer pairs for each direction. Even as currently costs seem to be symmetrical
			// this method calculates all 8 directions, from the PoV of the current supertile.
			interPlaneCosts( //
					true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
					planes, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			// Mark which links between neighbor planes are valid. Uses planesFit method, not costs
			// results used by resolveStarConflict(), getMergeSameTileCandidates(), setNonExclusive()
			// TODO: remove and switch to use connection costs
			// Sets both setMergedValid and setMergedStrongValid - should they be removed later
			filterNeighborPlanes(
					planes, // final TilePlanes.PlaneData [][] planes,
					true, // final boolean merge_low_eigen,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			// Set links between planes in 8 directions, in such a way that no start/end planes can be shared
			setExclusiveLinks(
					planes, // final TilePlanes.PlaneData [][] planes,
					getExNeibCost(), // final double                    max_cost,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			// Merge the supertile planes with agreeing neighbors, non-exclusively (no considering other planes
			// of the same supertile. Start with the best fit, then goes to lower quality, until the individual
			// merge quality falls below scaled quality of the best, pre-set minimum or the merged plane becomes
			// too thick.
			// Separately calculates merged weighted plane and with equal weights of the neighbors
			// TODO: Maybe just switch to the connection costs instead? 
			setNonExclusive(
					true, // final boolean                   en_sticks, // allow merging with bad plates
					planes, // final TilePlanes.PlaneData [][] planes,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			// Select initial pairs of plane merge candidates (same supertile plane). All pairs are ordered, first is lower than second
			// Candidates for merging share some non-exclusive neighbors
			int [][][] merge_candidates =  getMergeSameTileCandidates(
					planes, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			 // Verify merge candidates by evaluating overlaps. Close planes that have significant overlap on the image are mofre likely
			 // to be actually different planes, if they do not overlap it is likely to be the same one, just multiple separate parts of
			 // it mistakenly discriminated during initial plane generation from the tiles data.
			 // 
			 // Merging is also permitted if the planes are "weak and similar" - low strength and fit for merging 
			 // 
			 // Added a 'hack' - allow merging even overlapping tiles if they are close enough in the real world space
			boolean [][][] valid_candidates = overlapSameTileCandidates (
					planes, // final TilePlanes.PlaneData [][] planes,			
					merge_candidates,       // final int [][][] merge_candidates,
					plThickWorld, // 0.2, // final double     min_distance,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

		  boolean[][] weak_pairs = weakForeground(
					planes,           // final TilePlanes.PlaneData [][] planes,
					merge_candidates, // final int [][][]                merge_candidates,
					valid_candidates, // final boolean [][][]            valid_candidates,
					//					final int                       stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					plWeakFgStrength, // final double                    min_fg_strength,
					plWeakFgOutliers, // final int                       outliers,  // remove any glitches?
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			
//			 * Possible problem is that "normalizing" merge quality for low weights is not applicable for "star" plane that includes neighhbors
//			 * Switch to a single "cost" function (costSameTileConnectionsAlt())
// Still - how to merge stray tiles that do not have neighbors/star? Still merge them "old way"	(costSameTileConnections()) if at least 1 does not
//			have a "star"
			
			// Possible problem is that "normalizing" merge quality for low weights is not applicable for "star" plane that include neighbors
			// ALso calculates and outputs some of other costs - just for debugging
			// TODO: Switch to a single "cost" function (costSameTileConnectionsAlt())
			costSameTileConnections(
					planes, // final TilePlanes.PlaneData [][] planes,
					merge_candidates, // final int [][][] merge_candidates,
					valid_candidates, // final boolean [][][]   valid_candidates, // will be updated
					weak_pairs,       // final boolean [][]              weak_fg_pairs, 
					1.0,              // double                    relax,
					plWeakFgRelax,    // final double                    relax_weak_fg,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			// Calculate costs of merging planes of the same supertile and remove those that are exceed threshold
			// Also calculates and outputs some of other costs - just for debugging
			costSameTileConnectionsAlt(
					//5.0,  // final double         threshold,
					//10.0, // final double         threshold_nostar,
					getMergeCostStar(),   // relax_for_conflicts * 5.0,  // final double         threshold, // 
					getMergeCostNoStar(), //relax_for_conflicts * 10.0, // final double         threshold_nostar,

					planes, // final TilePlanes.PlaneData [][] planes,
					merge_candidates,       // final int [][][] merge_candidates,
					valid_candidates, // final boolean [][][]   valid_candidates, // will be updated
					weak_pairs, // final boolean [][]              weak_fg_pairs, 
					plWeakFgRelax,               // final double                    relax_weak_fg,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			// Extract non-conflicting groups of the same supertile planes that are OK to merge. For multiple intersecting pairs
			// this method starts for from the lowest cost merging and adds planes (from merge_candidates) until they do not cause
			// conflicts with already selected ones 
					
			int [][][] merge_groups = extractMergeSameTileGroups(
					planes,              // final TilePlanes.PlaneData [][] planes,
					merge_candidates,       // final int [][][] merge_candidates,
					valid_candidates, // boolean [][][] plane_overlaps,
					1.0,         // double                    relax,
					debugLevel + 1,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			// Save current planes for debugging, before merging

			if (dbg_orig_planes!=null) {
				dbg_orig_planes[num_merge_try] = planes.clone();
				for (int nsTile=0; nsTile < planes.length; nsTile++) if (planes[nsTile] != null){
					dbg_orig_planes[num_merge_try][nsTile] = planes[nsTile].clone();
					for (int np = 0; np < planes[nsTile].length; np++ ) if (planes[nsTile][np] != null){
						dbg_orig_planes[num_merge_try][nsTile][np] = planes[nsTile][np].clone(); 	
					}
				}
			}
			// Actually merge the planes by merging corresponding planes tile selections and re-building planes
			// Apply same supertile planes merge by combining tiles and re-generating ellipsoids by diagonalizing
			// covariance matrices. Some outliers may be removed after merge
			
			int num_removed_by_merging = st.applyMergePlanes(
					planes,    // final TilePlanes.PlaneData[][]   planes,
					merge_groups, // final int [][][]                 merge_groups,			
					// parameters to generate ellipsoids			
					0.0, // 3,                       // final double                     disp_far, // minimal disparity to select (or NaN)
					Double.NaN,                      // final double                     disp_near, // maximal disparity to select (or NaN)
					plDispNorm,       // final double                     dispNorm,   //  Normalize disparities to the average if above
					0.0,                             // final double                     min_weight,
					plMinPoints,      // final int                        min_tiles,
					// parameters to reduce outliers			
					plTargetEigen,    // final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					plFractOutliers,  // final double                     fractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
					plMaxOutliers,    // final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			if ( debugLevel > -1) {
				System.out.println("conditionSuperTiles(): try "+num_merge_try+ ": removed "+num_removed_by_merging+((num_removed_by_merging>0)?" planes by merging, recalculating connections":""));
			}
			if (num_removed_by_merging == 0){ // re-calculate all links
				break;

			}
		}
		return dbg_orig_planes;
	}
	


	public TilePlanes.PlaneData[][] copyPlanes(
			TilePlanes.PlaneData[][] src_planes)
	{
		TilePlanes.PlaneData[][] dst_planes = new TilePlanes.PlaneData[src_planes.length][];
		return copyPlanes(src_planes, dst_planes);
	}

	public TilePlanes.PlaneData[][] copyPlanes(
			final TilePlanes.PlaneData[][] src_planes,
			final TilePlanes.PlaneData[][] dst_planes)
	{
		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < src_planes.length; nsTile = ai.getAndIncrement()) {
						if (src_planes[nsTile] != null){
							dst_planes[nsTile] = new TilePlanes.PlaneData[src_planes[nsTile].length];
							for (int np = 0; np < src_planes[nsTile].length; np++){
								if (src_planes[nsTile][np] != null){
									dst_planes[nsTile][np] = src_planes[nsTile][np].clone(); 
								} else {
									dst_planes[nsTile][np] = null;
								}
							}
						} else {
							dst_planes[nsTile] = null;
						}
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return dst_planes;
	}

	
	
	
	public TilePlanes.PlaneData[][] planesSmooth(
			final TilePlanes.PlaneData[][] planes,
			TilePlanes.PlaneData[][] planes_mod, // should start with a copy of planes
			final int         debugLevel,
			final int         dbg_X,
			final int         dbg_Y)
	{
//		if (planes_mod == null){
//			planes_mod =copyPlanes(planes); // make always (for now) *********************
//		}
		final double maxDiff = Math.pow(10.0,  -plPrecision);
		for (int pass = 0; pass < plIterations; pass++){
			double diff = planesSmoothStep(
					planes,     // final TilePlanes.PlaneData[][] measured_planes,
					planes_mod, // final TilePlanes.PlaneData[][] mod_planes,
					true,            // final boolean calc_diff,
					(pass < 10)? debugLevel: 0,
							dbg_X,
							dbg_Y);
			if (diff < maxDiff){
				if (debugLevel > -1){
					System.out.println("planesSmooth(): pass:"+pass+" (of "+plIterations+"), rms = "+diff+" < "+maxDiff);
					break;
				}
			}
			if (debugLevel > -1){
				System.out.println("planesSmooth() -  pass:"+pass+" (of "+plIterations+"), rms = "+diff);
			}
		}
		return planes_mod;
	}			

	public double planesSmoothStep(
			final TilePlanes.PlaneData[][] measured_planes,
			final TilePlanes.PlaneData[][] mod_planes,
			final boolean calc_diff,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		
		final int [][] dirsYX = {{-1, 0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};
		final int tilesX =        st.tileProcessor.getTilesX();
//		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		final int stilesX = (tilesX + superTileSize -1)/superTileSize;  
//		final int stilesY = (tilesY + superTileSize -1)/superTileSize;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		final TilePlanes.PlaneData[][] new_planes = copyPlanes(mod_planes);
//		final Thread[] threads = ImageDtt.newThreadArray(tileProcessor.threadsMax);
		final Thread[] threads = ImageDtt.newThreadArray((debugLevel > 1)? 1 : st.tileProcessor.threadsMax);
		final int numThreads = threads.length;
		final double [] rslt_diffs = calc_diff ? new double [numThreads] : null; // all 0;
		final AtomicInteger ai_numThread = new AtomicInteger(0);
		final AtomicInteger ai = new AtomicInteger(0);
		final String [] titles= {
				"orig",     // 0
				"measured", // 1
				"n0","n1","n2","n3","n4","n5","n6","n7", // 2 .. 9
				"m0","m1","m2","m3","m4","m5","m6","m7","mm", // 10..18
		"diff"}; // 19

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					double [][] dbg_img=null;
					int numThread = ai_numThread.getAndIncrement(); // unique number of thread to write to rslt_diffs[numThread]
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < mod_planes.length; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
//						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
                        int dl = ((debugLevel > 1) && (nsTile0 == debug_stile)) ? 3: debugLevel;
						if ( new_planes[nsTile0] != null) {
							if (dl > 0){
								System.out.println("planesSmoothStep nsTile0="+nsTile0);
								dbg_img = new double [titles.length][];
							}
							int np0_min = (new_planes[nsTile0].length > 1) ? 1:0; // Modify if overall plane will be removed
							for (int np0 = np0_min; np0 < new_planes[nsTile0].length; np0 ++){
								TilePlanes.PlaneData this_new_plane = new_planes[nsTile0][np0];
								if (this_new_plane == null){
									System.out.println("Bug?  new_planes["+nsTile0+"]["+np0+"] == null");
									continue;
								}
								if (dl > 0) dbg_img[ 0] = this_new_plane.getSinglePlaneDisparity(false);
								if (dl > 0) dbg_img[ 1] = measured_planes[nsTile0][np0].getSinglePlaneDisparity(false);
								int [] neibs = this_new_plane.getNeibBest();
								double [][] costs = new double[neibs.length][];
								double [] weights = new double[neibs.length];
								int cost_index = 2; // overall cost
								double sum_rcosts = 0.0;
								int non_zero = 0;
								for (int dir = 0; dir < neibs.length; dir++) if (neibs[dir] >= 0) {
									int stx = stx0 + dirsYX[dir][1];
									int sty = sty0 + dirsYX[dir][0];
									int nsTile = sty * stilesX + stx; // from where to get
									TilePlanes.PlaneData other_plane = this_new_plane.getPlaneToThis(
											mod_planes[nsTile][neibs[dir]], // neighbor, previous value
											dl - 1); // debugLevel);
									costs[dir] = getFitQualities(
											// there probably be no sticks here after merging
											false, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
											this_new_plane, // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
											other_plane,   // TilePlanes.PlaneData plane2,
											Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
											Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
											Double.NaN, // double               merged_wev,    // if NaN will calculate assuming the same supertile - for world
											Double.NaN, // double               merged_wev_eq, // if NaN will calculate assuming the same supertile - for world
											nsTile0+":"+np0+"-"+dir,  // String prefix,
											0); // int debugLevel)
									for (int i = 3; i < costs[dir].length; i++) costs[dir][i] *= costs[dir][2];
									non_zero ++;
									weights[dir] = 1.0/costs[dir][cost_index];
									sum_rcosts += weights[dir]; 
								}
								for (int dir = 0; dir < neibs.length; dir++) if (neibs[dir] >= 0) {
									weights[dir] *= non_zero/sum_rcosts; // average weight for active links will be 1.0 
								}								
								if (dl > 0) {
									for (int dir = 0; dir <8; dir++)if (costs[dir] != null){
										System.out.print(nsTile0+":"+np0+"-"+dir+"  "+String.format(" weight= %6.3f ",weights[dir]));
										for (int i = 0; i < costs[dir].length; i++){
											System.out.print(String.format("%8.3f", costs[dir][i]));
//											if (i < 3)	System.out.print(String.format("%8.3f", costs[dir][i]));
//											else System.out.print(String.format("%8.3f", costs[dir][i]*costs[dir][2]));
										}
										System.out.println();
									}
									System.out.println();
								}
								
								this_new_plane =this_new_plane.clone(); // not to change weight!
								this_new_plane.setWeight(0.0); //
								double num_merged = 0.0; // double to add fractional pull weight of the center
								double true_num_merged = 0.0;
								for (int dir = 0; dir < neibs.length; dir++){
									if (neibs[dir] >= 0) {
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										int nsTile = sty * stilesX + stx; // from where to get
										TilePlanes.PlaneData other_plane = this_new_plane.getPlaneToThis(
												mod_planes[nsTile][neibs[dir]], // neighbor, previous value
												dl - 1); // debugLevel);
										if (dl > 0) dbg_img[ 2 + dir] = other_plane.getSinglePlaneDisparity(false);
										if ((other_plane != null) && ((other_plane.getValue() <= plMaxEigen) || (plMaxEigen == 0))) { // TODO:
											if (this_new_plane.getWeight() > 0.0){
												this_new_plane = this_new_plane.mergePlaneToThis(
														other_plane, // PlaneData otherPd,
														weights[dir], // 1.0,         // double    scale_other,
														// here it should be no power function for the weights
														1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost
														false,       // boolean   ignore_weights,
														true,        // boolean   sum_weights,
														plPreferDisparity,
														dl - 1); // int       debugLevel)
												
											} else {
												this_new_plane.copyNeib(this_new_plane, other_plane); // keep neighbors of the original center plane
												this_new_plane.copyStar(this_new_plane, other_plane);
												this_new_plane = other_plane; // should increment num_merged
												this_new_plane.scaleWeight(weights[dir]);
												this_new_plane.invalidateCalculated();
											}
											if (this_new_plane != null){
												num_merged += 1.0;
												true_num_merged += 1.0;

												// just for debug / calculate
												this_new_plane.getWorldXYZ(0);
											}
											new_planes[nsTile0][np0] = this_new_plane;
											if (dl > 0) dbg_img[10 + dir] = this_new_plane.getSinglePlaneDisparity(false);
										} else if (plStopBad){ // just skip, not abandon? (set to false)
											this_new_plane = null;
											break;
										}
									}
								}
								if (this_new_plane != null) {
									// average weight over participating directions, so the relative pull
									// does not depend on number of neighbors
									if ((num_merged > 0.0) && (this_new_plane != null) && (plNormPow > 1.0)){
										double scale = Math.pow(num_merged, plNormPow);
										this_new_plane.scaleWeight(1.0/scale);
										num_merged /=scale;
										true_num_merged /= scale;

									}
									
									if (    (plPull > 0.0) &&
											(measured_planes != null) &&
											(measured_planes[nsTile0] != null) &&
											(measured_planes[nsTile0][np0] != null) &&
//											((measured_planes[nsTile0][np0].getValue() < plMaxEigen) || (plMaxEigen == 0))
											((measured_planes[nsTile0][np0].getValue() < plMaxEigen) || (plMaxEigen == 0) || 
													(this_new_plane == null) || (this_new_plane.getWeight() == 0.0)) // keep measured if impossible to merge
											){ // merge with "measured"

										if (this_new_plane.getWeight() > 0.0){
											this_new_plane = this_new_plane.mergePlaneToThis(
													measured_planes[nsTile0][np0], // PlaneData otherPd,
													plPull,                     // double    scale_other,
													1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost													
													false,                         // boolean   ignore_weights,
													true, // boolean   sum_weights,
													plPreferDisparity,
													dl - 1);                           // int       debugLevel)
											if (this_new_plane != null){
												num_merged += plPull;	// num_merged was 1.0 and weight is averaged over all neighbors
												true_num_merged += 1.0;
											}
										} else {
											this_new_plane = measured_planes[nsTile0][np0].clone();
											num_merged = 1.0;
											true_num_merged = 1.0;
										}
										new_planes[nsTile0][np0] = this_new_plane;
										if (dl > 0) dbg_img[18] =  this_new_plane.getSinglePlaneDisparity(false);
									}
									if ((num_merged > 0.0) && (this_new_plane != null)){
										this_new_plane.scaleWeight(1.0/num_merged);
//										double true_num_merged = num_merged - plPull + 1;
										this_new_plane.setNumPoints((int) (this_new_plane.getNumPoints()/true_num_merged));
									}
									// Revert if the result value is higher than imposed maximum
									if ((this_new_plane.getValue() > plMaxEigen)  && (plMaxEigen != 0)){ // TODO: Set more relaxed here?
										if (dl > 0){
											System.out.println("planesSmoothStep nsTile0="+nsTile0+" smoothed plane is too thick, using previous one");
											dbg_img = new double [titles.length][];
										}
										this_new_plane = mod_planes[nsTile0][np0].clone();
										new_planes[nsTile0][np0] = this_new_plane;
									}
//									Use non-exclusive
									
									
									// just for debug / calculate 
									this_new_plane.getWorldXYZ(0);
									// calculate largest disparity difference between old and new plane
									if (rslt_diffs != null){ // filter out outliers here?
										// get plane for both old and new, calc rms of diff
										double [] oldPlane = mod_planes[nsTile0][np0].getSinglePlaneDisparity(
												false); // use_NaN)
										double [] newPlane = new_planes[nsTile0][np0].getSinglePlaneDisparity(
												false); // use_NaN)
										double s = 0.0;
										for (int i = 0; i < oldPlane.length; i++){
											double d = newPlane[i] - oldPlane[i];
											s+= d*d;
										}
										s= Math.sqrt(s/oldPlane.length);
										if (s > rslt_diffs[numThread]){
											rslt_diffs[numThread] = s;
										}
										if (dl > 0) {
											dbg_img[19] =  new double[oldPlane.length];
											for (int i = 0; i < oldPlane.length; i++){
												dbg_img[19][i] = newPlane[i] - oldPlane[i];
											}
											if (debugLevel > 0) {
												System.out.println("planesSmoothStep() nsTile0="+nsTile0+" rms = "+s);
											}
										}
										if (debugLevel > -1) {
											//										if ((s > 5.0) || (dl > 0)){
											if ((s > 5.0) || (dl > 0)){
												System.out.println("planesSmoothStep() nsTile0="+nsTile0+":"+np0+
														" num_merged="+num_merged + " rms = "+s+
														" new_weight = "+new_planes[nsTile0][np0].getWeight()+
														" old_weight = "+mod_planes[nsTile0][np0].getWeight()+
														" new_value = "+new_planes[nsTile0][np0].getValue()+
														" old_value = "+mod_planes[nsTile0][np0].getValue()+
														" measured_weight = "+measured_planes[nsTile0][np0].getWeight()+
														" measured_value = "+measured_planes[nsTile0][np0].getValue());
											}
										}

									}
									if ((dl > 0) && (debugLevel > 0)){
										showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
										sdfa_instance.showArrays(dbg_img, superTileSize, superTileSize, true, "smooth_step_x"+stx0+"_y"+sty0, titles);
									}
								} else { // if (this_new_plane != null)
									this_new_plane = mod_planes[nsTile0][np0].clone();
									new_planes[nsTile0][np0] = this_new_plane;
								}
							}
						}								
					}
				}
			};
		}		      
		ImageDtt.startAndJoin(threads);
		copyPlanes (new_planes, mod_planes); // copy back
		if (rslt_diffs == null){
			return Double.NaN;
		}
		double diff = 0.0;
		for (int i = 0; (i < numThreads) ; i++) if (diff < rslt_diffs[i]) diff = rslt_diffs[i];
		return diff; // return maximal difference
	}

	
	public TilePlanes.PlaneData[][] planesSmoothAndMerge(
			final TilePlanes.PlaneData[][] planes,
			final int                      max_num_tries,
			final int                      debugLevel)
	{
		TilePlanes.PlaneData[][] planes_mod = null;
		for (int num_merge_try = 0; num_merge_try < max_num_tries; num_merge_try ++ ) { // smooth and merge
			
			planes_mod =copyPlanes(planes); // Start with fresh copied planes extracted measure data, nort smoothed 
			
			planes_mod = planesSmooth(
					planes,      // final TilePlanes.PlaneData[][] planes,
					planes_mod,  // TilePlanes.PlaneData[][] planes_mod,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			// create costs for the modified planes
			interPlaneCosts(
					true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
					planes_mod, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			setExclusiveLinks(
					planes_mod, // final TilePlanes.PlaneData [][] planes,
					//						2.5, //final double                      max_cost
					getExNeibCost() * getExNeibSmooth(), // final double                    max_cost,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			// once more after updating exclusive links
			planes_mod = planesSmooth(
					planes,      // final TilePlanes.PlaneData[][] planes,
					planes_mod,  // TilePlanes.PlaneData[][] planes_mod,
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			interPlaneCosts(
					true,          // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
					planes_mod, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,    // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			// recalculate links? more smooth?
			setExclusiveLinks(
					planes_mod,                          // final TilePlanes.PlaneData [][] planes,
					getExNeibCost() * getExNeibSmooth(), // final double                    max_cost,
					debugLevel,                          // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			// just in case? Not yet needed			
			setNonExclusive(
					true,                                // final boolean                   en_sticks, // allow merging with bad plates
					planes_mod,                          // final TilePlanes.PlaneData [][] planes,
					debugLevel,                          // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			// see if some modified planes need to be merged (but merge originals)
			// TODO: Stricter requirements for merging here than for original planes?				

			int [][][] merge_candidates =  getMergeSameTileCandidates(
					planes_mod, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,                          // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			boolean [][][] valid_candidates = overlapSameTileCandidates (
					planes_mod, // final TilePlanes.PlaneData [][] planes,			
					merge_candidates,       // final int [][][] merge_candidates,
					plThickWorld, // 0.2, // final double     min_distance,
					debugLevel,                          // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			// Consider supertiles with conflicts, merge conflicting layers with relaxed requirements
			
			Conflicts iconflicts0 = new Conflicts(st); 
			int [][][] conflicts0 = iconflicts0.detectTriangularConflicts(
					debugLevel); // 1); // final int debugLevel)

			// Will be used later, now just use before merge candidates are modified
			int [][][] conflicting_candidates = filterPairsByConflicts(
					planes_mod,    // final TilePlanes.PlaneData [][] planes,			
					merge_candidates, // final int [][][]                merge_candidates,
					conflicts0);      // final int [][][]                conflicts)

			boolean[][] weak_pairs = weakForeground(
					planes,           // final TilePlanes.PlaneData [][] planes,
					merge_candidates, // final int [][][]                merge_candidates,
					valid_candidates, // final boolean [][][]            valid_candidates,
					//					final int                       stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
					plWeakFgStrength, // final double                    min_fg_strength,
					plWeakFgOutliers, // final int                       outliers,  // remove any glitches?
					debugLevel,                  // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			

			//				 * Possible problem is that "normalizing" merge quality for low weights is not applicable for "star" plane that include neighhbors
			//				 * Switch to a single "cost" function (costSameTileConnectionsAlt())
			// Still - how to merge stray tiles that do not have neighbors/star? Still merge them "old way"	(costSameTileConnections()) if at least 1 does not
			//				have a "star"
			costSameTileConnections(
					planes_mod, // final TilePlanes.PlaneData [][] planes,
					merge_candidates,       // final int [][][] merge_candidates,
					valid_candidates,       // final boolean [][][]   valid_candidates, // will be updated
					weak_pairs, // final boolean [][]              weak_fg_pairs, 
					1.0,                    // final double                    relax,
					plWeakFgRelax,          // final double                    relax_weak_fg,
					debugLevel,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			costSameTileConnectionsAlt(
					getMergeCostStar(),   // relax_for_conflicts * 5.0,  // final double         threshold, // 
					getMergeCostNoStar(), //relax_for_conflicts * 10.0, // final double         threshold_nostar,
					planes_mod,           // final TilePlanes.PlaneData [][] planes,
					merge_candidates,     // final int [][][] merge_candidates,
					valid_candidates,     // final boolean [][][]   valid_candidates, // will be updated
					weak_pairs,           // final boolean [][]              weak_fg_pairs, 
					plWeakFgRelax,        // final double                    relax_weak_fg,
					debugLevel,           // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			int [][][] merge_groups = extractMergeSameTileGroups(
					planes_mod,              // final TilePlanes.PlaneData [][] planes,
					merge_candidates,       // final int [][][] merge_candidates,
					valid_candidates, // boolean [][][] plane_overlaps,
					1.0,              // final double                    relax,
					debugLevel +1,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			// Apply merging to the original planes, not the smoothed ones
			int num_removed_by_merging = st.applyMergePlanes(
					planes,    // final TilePlanes.PlaneData[][]   planes,
					merge_groups, // final int [][][]                 merge_groups,			
					// parameters to generate ellipsoids			
					0.0, // 3,        // final double                     disp_far, // minimal disparity to select (or NaN)
					Double.NaN,       // final double                     disp_near, // maximal disparity to select (or NaN)
					plDispNorm,       // final double                     dispNorm,   //  Normalize disparities to the average if above
					0.0,              // final double                     min_weight,
					plMinPoints,      // final int                        min_tiles,
					// parameters to reduce outliers			
					plTargetEigen,    // final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
					plFractOutliers,  // final double                     fractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
					plMaxOutliers,    // final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
					debugLevel,       // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
// TODO: above is very similar to conditionSuperTiles(). Combine
			
			System.out.println("planesSmoothAndMerge(): try "+num_merge_try+ ": removed "+num_removed_by_merging+
					((num_removed_by_merging>0)?" planes by merging, recalculating connections":""));





			if (num_removed_by_merging == 0){ // re-calculate all links
				// Consider supertiles with conflicts, merge conflicting layers with relaxed requirements

				// just to list them
				if (debugLevel > 0){
					Conflicts conflicts0_stats =  new Conflicts(
							conflicts0,
							st,
							-1); // debugLevel);
				}
				System.out.println("Trying relaxed merging for conflicting plane pairs");

				valid_candidates = overlapSameTileCandidates (
						planes_mod, // final TilePlanes.PlaneData [][] planes,			
						conflicting_candidates,       // final int [][][] merge_candidates,\\
						plThickWorldConfl, // 0.4, // final double     min_distance,
						debugLevel,       // final int        debugLevel)
						dbg_tileX,
						dbg_tileY);

				// again same sequence
				boolean [][] confl_weak_pairs = weakForeground(
							planes,                 // final TilePlanes.PlaneData [][] planes,
							conflicting_candidates, // final int [][][]                merge_candidates,
							valid_candidates,       // final boolean [][][]            valid_candidates,
							//					final int                       stMeasSel, //            = 1;      // Select measurements for supertiles : +1 - combo, +2 - quad +4 - hor +8 - vert
							plWeakFgStrength, // final double                    min_fg_strength,
							plWeakFgOutliers, // final int                       outliers,  // remove any glitches?
							debugLevel,                  // final int        debugLevel)
							dbg_tileX,
							dbg_tileY);
				
				// try to merge original (measured) planes, not smoothed ones ??
				costSameTileConnections(
						planes,               // final TilePlanes.PlaneData [][] planes,
						conflicting_candidates,  // final int [][][] merge_candidates,
						valid_candidates,        // final boolean [][][]   valid_candidates, // will be updated
						confl_weak_pairs,        // final boolean [][]              weak_fg_pairs, 
						getConflRelax(),         // relax_for_conflicts,         // final double                    relax,
						plWeakFgRelax,           // final double                    relax_weak_fg,
						debugLevel, // 2,                       // -1, // debugLevel,                  // final int        debugLevel)
						dbg_tileX,
						dbg_tileY);

				costSameTileConnectionsAlt(
						getConflRelax() * getMergeCostStar(),   // relax_for_conflicts * 5.0,  // final double         threshold, // 
						getConflRelax() * getMergeCostNoStar(), //relax_for_conflicts * 10.0, // final double         threshold_nostar,

						planes, // final TilePlanes.PlaneData [][] planes,
						conflicting_candidates, // final int [][][] merge_candidates,
						valid_candidates,       // final boolean [][][]   valid_candidates, // will be updated
						confl_weak_pairs,       // final boolean [][]              weak_fg_pairs, 
						plWeakFgRelax,          // final double                    relax_weak_fg,
						debugLevel,             // final int        debugLevel)
						dbg_tileX,
						dbg_tileY);

				merge_groups = extractMergeSameTileGroups(
						planes,                 // final TilePlanes.PlaneData [][] planes,
						conflicting_candidates, // final int [][][] merge_candidates,
						valid_candidates,       // boolean [][][] plane_overlaps,
						getConflRelax(),        // relax_for_conflicts,         // final double                    relax,
						debugLevel + 1,          // final int        debugLevel)
						dbg_tileX,
						dbg_tileY);

				num_removed_by_merging = st.applyMergePlanes(
						planes,    // final TilePlanes.PlaneData[][]   planes,
						merge_groups, // final int [][][]                 merge_groups,			
						// parameters to generate ellipsoids			
						0.0, // 3,                       // final double                     disp_far, // minimal disparity to select (or NaN)
						Double.NaN,                      // final double                     disp_near, // maximal disparity to select (or NaN)
						plDispNorm,       // final double                     dispNorm,   //  Normalize disparities to the average if above
						0.0,                             // final double                     min_weight,
						plMinPoints,      // final int                        min_tiles,
						// parameters to reduce outliers			
						plTargetEigen,    // final double                     targetEigen,   //     =   0.1;  // Remove outliers until main axis eigenvalue (possibly scaled by plDispNorm) gets below
						plFractOutliers,  // final double                     fractOutliers, //     =   0.3;  // Maximal fraction of outliers to remove
						plMaxOutliers,    // final int                        maxOutliers,   //     =   20;  // Maximal number of outliers to remove
						debugLevel,             // final int        debugLevel)
						dbg_tileX,
						dbg_tileY);
				
				System.out.println("planesSmoothAndMerge(): try "+num_merge_try+ ": removed "+num_removed_by_merging+
						((num_removed_by_merging>0)?" conflicting planes by merging, recalculating connections":""));
				
				if ( num_merge_try >= max_num_tries) {
					System.out.println("Exceeded maximal number of iterations, beaking anyway...");
					break;
				}
				if (num_removed_by_merging == 0){ // re-calculate all links
					break;
				}
			}

			// the following is not actually needed, just to keep measured (not smoothed) plane data consistent ?
			
			// Do the same as in conditionSuperTiles before smoothing again
			// Some actions below may be duplicate?

			matchPlanes(
					planes, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
			
			interPlaneCosts( //
					true, // final boolean                   en_sticks, // treat planes with second eigenvalue below plEigenStick as "sticks"
					planes, // final TilePlanes.PlaneData [][] planes,			
					debugLevel,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			filterNeighborPlanes(
					planes, // final TilePlanes.PlaneData [][] planes,
					true, // final boolean merge_low_eigen,
					debugLevel,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			setExclusiveLinks( // stricter? why?
					planes, // final TilePlanes.PlaneData [][] planes,
					getExNeibCost()*getExNeibSmooth(), // final double                    max_cost,
					debugLevel,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);

			setNonExclusive(
					true, // final boolean                   en_sticks, // allow merging with bad plates
					planes, // final TilePlanes.PlaneData [][] planes,
					debugLevel,             // final int        debugLevel)
					dbg_tileX,
					dbg_tileY);
		}
		return planes_mod;
	}	

	public boolean goodLink(
			double threshold,
			int nsTile,
			int np,
			int dir,
			int dest_plane,
			TilePlanes.PlaneData [][] planes,
			TileNeibs tnSurface)
	{
		if (threshold == 0.0) return true;
		double cost = getMutualLinkCost(
				nsTile,
				np,
				dir,
				dest_plane,
				planes,
				tnSurface);
		return cost <= threshold; // cost may be Double.NaN, will return false
	}
	
	public double getMutualLinkCost(
			int nsTile,
			int np,
			int dir,
			int dest_plane,
			TilePlanes.PlaneData [][] planes,
			TileNeibs tnSurface)
	{
		double cost = planes[nsTile][np].getLinkCosts(dir, dest_plane);
		if (!Double.isNaN(cost)){
			int back_dir = (dir + 4) % 8;
			cost = 0.5 * (cost + planes[tnSurface.getNeibIndex(nsTile, dir)][dest_plane].getLinkCosts(back_dir, np));
		}
		return cost;
	}
	
	public int fillSquares(
			final TilePlanes.PlaneData[][] planes,
			final double                   threshold,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final int debug_stile = dbg_Y * stilesX + dbg_X;
//		final int nStiles =       stilesX * stilesY;
		int num_added = 0;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				int [] nsTiles = { nsTile, nsTile + 1, nsTile + stilesX, nsTile + stilesX + 1};
                int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
				if (dl > 2){
					System.out.println("fillSquares(): nsTile="+nsTile);
				}
				TilePlanes.PlaneData [][] planes4 = {planes[nsTiles[0]],planes[nsTiles[1]],planes[nsTiles[2]],planes[nsTiles[3]]};
				if ((planes4[0] != null) && (planes4[1] != null) && (planes4[2] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if ((neibs0[2] < 0) || (neibs0[4] < 0)) continue;
						if (planes4[1][neibs0[2]] == null) continue;
						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
						int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
						if ((neibs1[4] < 0) || (neibs2[2] < 0) || (neibs1[4] != neibs2[2])) continue;
						int [] neibs3 = planes4[3][neibs1[4]].getNeibBest();
						// got full square, does it already have diagonals
						if ((neibs0[3] < 0) && (neibs3[7] < 0)){
							// See if the link cost is below threshold
							// from this to right-down and back
							if (goodLink(threshold,nsTiles[0], np, 3, neibs1[4], planes, tnSurface )) {
								neibs0[3] = neibs1[4];
								neibs3[7] = np;
								num_added++;
								if (dl > 0){
									System.out.println("fillSquares().1: added link "+
											nsTiles[0]+":" + np+ "-(3)->"+neibs1[4]+", " + nsTiles[3]+":" + neibs1[4]+"-(7)->" + np+
											" cost = "+getMutualLinkCost(nsTiles[0], np, 3, neibs1[4], planes, tnSurface ));												
								}
							}
						}
						if ((neibs1[5] < 0) && (neibs2[1] < 0)){
							// See if the link cost is below threshold
							// from right to down and back
							if (goodLink(threshold, nsTiles[1], neibs0[2], 5, neibs0[4], planes, tnSurface )) {
								neibs1[5] = neibs0[4];
								neibs2[1] = neibs0[2];
								num_added++;
								if (dl > 0){
									System.out.println("fillSquares().2: added link "+
											nsTiles[1]+":" + neibs0[2]+"-(5)->"+neibs0[4]+", " + nsTiles[2]+":" + neibs0[4]+"-(1)->" + neibs0[2]+
											" cost = " + getMutualLinkCost(nsTiles[1], neibs0[2], 5, neibs0[4], planes, tnSurface ));												
								}
							}
						}
					}
				}
			}
		}
		return num_added;
	}
	
	public int fillHypotenuse(
			final TilePlanes.PlaneData[][] planes,
			final double                   threshold,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		final int debug_stile = dbg_Y * stilesX + dbg_X;
//		final int nStiles =       stilesX * stilesY;
		int num_added = 0;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				int [] nsTiles = { nsTile, nsTile + 1, nsTile + stilesX, nsTile + stilesX + 1};
                int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
				if (dl > 2){
					System.out.println("fillHypotenuse(): nsTile="+nsTile);
				}
				TilePlanes.PlaneData [][] planes4 = {planes[nsTiles[0]],planes[nsTiles[1]],planes[nsTiles[2]],planes[nsTiles[3]]};
				// when right and down links exist
				if ((planes4[0] != null) && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if ((neibs0[2] < 0) || (neibs0[4] < 0)) continue;
//						if (planes4[1][neibs0[2]] == null) continue; // should never happen 
//						if (planes4[2][neibs0[4]] == null) continue; // should never happen 
						int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
						int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
						if ((neibs1[5] < 0) && (neibs2[1] < 0)){
							// See if the link cost is below threshold
							// from right to down and back
							if (goodLink(threshold, nsTiles[1], neibs0[2], 5, neibs0[4], planes, tnSurface )) {
								neibs1[5] = neibs0[4];
								neibs2[1] = neibs0[2];
								num_added++;
								if (dl > 0){
									System.out.println("fillHypotenuse().1: added link "+
											nsTiles[1]+":" + neibs0[2]+"-(5)->"+neibs0[4]+", " + nsTiles[2]+":" + neibs0[4]+"-(1)->" + neibs0[2]+
											" cost = " + getMutualLinkCost(nsTiles[1], neibs0[2], 5, neibs0[4], planes, tnSurface ));												
								}
							}
						}
					}
				}
				// when right and down from right links exist
				if ((planes4[0] != null) && (planes4[1] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if (neibs0[2] < 0) continue;
//						if (planes4[1][neibs0[2]] == null) continue;
//						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
						if (neibs1[4] < 0) continue;
						int [] neibs3 = planes4[3][neibs1[4]].getNeibBest();
						// got right then down, connect if possible
						if ((neibs0[3] < 0) && (neibs3[7] < 0)){
							// See if the link cost is below threshold
							// from this to right-down and back
							if (goodLink(threshold,nsTiles[0], np, 3, neibs1[4], planes, tnSurface )) {
								neibs0[3] = neibs1[4];
								neibs3[7] = np;
								num_added++;
								if (dl > 0){
									System.out.println("fillHypotenuse().2: added link "+
											nsTiles[0]+":" + np+ "-(3)->"+neibs1[4]+", " + nsTiles[3]+":" + neibs1[4]+"-(7)->" + np+
											" cost = "+getMutualLinkCost(nsTiles[0], np, 3, neibs1[4], planes, tnSurface ));												
								}
							}
						}
					}
				}
				// when down and right from down links exist
				if ((planes4[0] != null) && (planes4[2] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if (neibs0[4] < 0) continue;
//						if (planes4[1][neibs0[2]] == null) continue;
//						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
						if (neibs2[2] < 0) continue;
						int [] neibs3 = planes4[3][neibs2[2]].getNeibBest();
						// got down then right, connect if possible
						if ((neibs0[3] < 0) && (neibs3[7] < 0)){
							// See if the link cost is below threshold
							// from this to right-down and back
							if (goodLink(threshold,nsTiles[0], np, 3, neibs2[2], planes, tnSurface )) {
								neibs0[3] = neibs2[2];
								neibs3[7] = np;
								num_added++;
								if (dl > 0){
									System.out.println("fillHypotenuse().3: added link "+
											nsTiles[0]+":" + np+ "-(3)->"+neibs2[2]+", " + nsTiles[3]+":" + neibs2[2]+"-(7)->" + np+
											" cost = "+getMutualLinkCost(nsTiles[0], np, 3, neibs2[2], planes, tnSurface ));												
								}
							}
						}
					}
				}
				// when both ortho links exist from the opposite ([3]) tile
				if ((planes4[1] != null) && (planes4[2] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[3].length; np++){
						if (planes4[3][np] == null) continue;
						int [] neibs3 = planes4[3][np].getNeibBest();
						if ((neibs3[6] < 0) || (neibs3[0] < 0)) continue;
//						if (planes4[1][neibs0[2]] == null) continue;
//						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs1 = planes4[1][neibs3[0]].getNeibBest();
						int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();
						if ((neibs1[5] < 0) && (neibs2[1] < 0)){
							// See if the link cost is below threshold
							// from right to down and back
							if (goodLink(threshold, nsTiles[1], neibs3[0], 5, neibs3[6], planes, tnSurface )) {
								neibs1[5] = neibs3[6];
								neibs2[1] = neibs3[0];
								num_added++;
								if (dl > 0){
									System.out.println("fillSquares().2: added link "+
											nsTiles[1]+":" + neibs3[0]+"-(5)->"+neibs3[6]+", " + nsTiles[2]+":" + neibs3[6]+"-(1)->" + neibs3[0]+
											" cost = " + getMutualLinkCost(nsTiles[1], neibs3[0], 5, neibs3[6], planes, tnSurface ));												
								}
							}
						}
					}
				}

				
				
				if ((planes4[0] != null) && (planes4[1] != null) && (planes4[2] != null) && (planes4[3] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] == null) continue;
						int [] neibs0 = planes4[0][np].getNeibBest();
						if ((neibs0[2] < 0) || (neibs0[4] < 0)) continue;
						if (planes4[1][neibs0[2]] == null) continue;
						if (planes4[2][neibs0[4]] == null) continue;
						int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
						int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
						if ((neibs1[4] < 0) || (neibs2[2] < 0) || (neibs1[4] != neibs2[2])) continue;
						int [] neibs3 = planes4[3][neibs1[4]].getNeibBest();
						// got full square, does it already have diagonals
						if ((neibs0[3] < 0) && (neibs3[7] < 0)){
							// See if the link cost is below threshold
							// from this to right-down and back
							if (goodLink(threshold,nsTiles[0], np, 3, neibs1[4], planes, tnSurface )) {
								neibs0[3] = neibs1[4];
								neibs3[7] = np;
								num_added++;
								if (dl > 0){
									System.out.println("fillSquares().1: added link "+
											nsTiles[0]+":" + np+ "-(3)->"+neibs1[4]+", " + nsTiles[3]+":" + neibs1[4]+"-(7)->" + np+
											" cost = "+getMutualLinkCost(nsTiles[0], np, 3, neibs1[4], planes, tnSurface ));												
								}
							}
						}
						if ((neibs1[5] < 0) && (neibs2[1] < 0)){
							// See if the link cost is below threshold
							// from right to down and back
							if (goodLink(threshold, nsTiles[1], neibs0[2], 5, neibs0[4], planes, tnSurface )) {
								neibs1[5] = neibs0[4];
								neibs2[1] = neibs0[2];
								num_added++;
								if (dl > 0){
									System.out.println("fillSquares().2: added link "+
											nsTiles[1]+":" + neibs0[2]+"-(5)->"+neibs0[4]+", " + nsTiles[2]+":" + neibs0[4]+"-(1)->" + neibs0[2]+
											" cost = " + getMutualLinkCost(nsTiles[1], neibs0[2], 5, neibs0[4], planes, tnSurface ));												
								}
							}
						}
					}
				}
			}
		}
		return num_added;
	}

	
	public int cutCorners(
			final TilePlanes.PlaneData[][] planes,
			final double                   threshold,
			final int debugLevel,
			final int dbg_X,
			final int dbg_Y)
	{
		final int tilesX =        st.tileProcessor.getTilesX();
		final int tilesY =        st.tileProcessor.getTilesY();
		final int superTileSize = st.tileProcessor.getSuperTileSize();
		//				final int tileSize =      tileProcessor.getTileSize();
		final int stilesX =       (tilesX + superTileSize -1)/superTileSize;  
		final int stilesY =       (tilesY + superTileSize -1)/superTileSize;
//		final int nStiles =       stilesX * stilesY;
		final TileNeibs tnSurface = new TileNeibs(stilesX, stilesY);
		int num_added = 0;
		final int debug_stile = dbg_Y * stilesX + dbg_X;
		for (int stY = 0; stY < (stilesY - 1); stY++ ) {
			for (int stX = 0; stX < (stilesX - 1); stX++ ) {
				int nsTile = stY * stilesX + stX;
				int dl = ((debugLevel > 1) && (nsTile == debug_stile)) ? 3: debugLevel;
				if (dl > 2){
					System.out.println("cutCorners(): nsTile="+nsTile);
				}
				int [] nsTiles = { nsTile, nsTile + 1, nsTile + stilesX, nsTile + stilesX + 1};
				TilePlanes.PlaneData [][] planes4 = {planes[nsTiles[0]],planes[nsTiles[1]],planes[nsTiles[2]],planes[nsTiles[3]]};
				if (planes4[0] != null) { // && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[0].length; np++){
						if (planes4[0][np] != null) {
							int [] neibs0 = planes4[0][np].getNeibBest();
							if ((neibs0[2] >= 0) && (neibs0[4] < 0)){
								if (planes4[1][neibs0[2]] != null) {
									int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
									if (neibs1[5] >= 0) {
										int [] neibs2 = planes4[2][neibs1[5]].getNeibBest();
										// from this to down and back
										if (neibs2[0] < 0) {
											if (goodLink(threshold, nsTiles[0], np, 4, neibs1[5], planes, tnSurface )) {
												neibs0[4] = neibs1[5];
												neibs2[0] = np;  //?
												num_added++;
												if (dl > 0){
													System.out.println("cutCorners().1: added link "+
															nsTiles[0]+":" + np +"-(4)->"+neibs1[5]+", "+nsTiles[2]+":" + neibs1[5]+"-(0)->"+np+
															" cost = " + getMutualLinkCost(nsTiles[0], np, 4, neibs1[5], planes, tnSurface ));												
												}
											}
										}
									}
								}
							}

							if ((neibs0[2] < 0) && (neibs0[4] >= 0)){
								if (planes4[2][neibs0[4]] != null) {
									int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
									if (neibs2[1] >= 0) {
										int [] neibs1 = planes4[1][neibs2[1]].getNeibBest();
										// from this to right and back
										if (neibs1[6] < 0) {
											if (goodLink(threshold, nsTiles[0], np, 2, neibs2[1], planes, tnSurface )) {
												neibs0[2] = neibs2[1];
												neibs1[6] = np;
												num_added++;
												if (dl > 0){
													System.out.println("cutCorners().2: added link "+
															nsTiles[0]+":" + np+"-(2)->"+neibs2[1]+", "+nsTiles[1]+":" + neibs2[1]+"-(6)->"+np+
															" cost = " + getMutualLinkCost(nsTiles[0], np, 2, neibs2[1], planes, tnSurface ));												
												}
											}
										}
									}
								}
							}

							if ((neibs0[2] >= 0) && (neibs0[3] >= 0)){
								int [] neibs1 = planes4[1][neibs0[2]].getNeibBest();
								if (neibs1[4] < 0) {
									int [] neibs3 = planes4[3][neibs0[3]].getNeibBest();
									// from right to right-down and back
									if (neibs3[0] < 0) {
										if (goodLink(threshold, nsTiles[1], neibs0[2], 4, neibs0[3], planes, tnSurface )) {
											neibs1[4] = neibs0[3];
											neibs3[0] = neibs0[2];
											num_added++;
											if (dl > 0){
												System.out.println("cutCorners().3: added link "+
														nsTiles[1]+":" + neibs0[2]+"-(4)->"+neibs0[3]+", " + nsTiles[3]+":" + neibs0[3]+"-(0)->" + neibs0[2]+
														" cost = " + getMutualLinkCost(nsTiles[1], neibs0[2], 4, neibs0[3], planes, tnSurface ));												
											}
										}
									}
								}
							}
							
							if ((neibs0[3] >= 0) && (neibs0[4] >= 0)){
								int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
								if (neibs2[2] < 0) {
									int [] neibs3 = planes4[3][neibs0[3]].getNeibBest();
									// from down to right-down and back
									if (neibs3[6] < 0) {
										if (goodLink(threshold, nsTiles[2], neibs0[4], 2, neibs0[3], planes, tnSurface )) {
											neibs2[2] = neibs0[3];
											neibs3[6] = neibs0[4];
											num_added++;
											if (dl > 0){
												System.out.println("cutCorners().4: added link "+
														nsTiles[2]+":" + neibs0[4]+"-(2)->"+neibs0[3]+", " + nsTiles[3]+":" + neibs0[3]+"-(6)->" + neibs0[4]+
														" cost = " + getMutualLinkCost(nsTiles[2], neibs0[4], 2, neibs0[3], planes, tnSurface ));												
											}
										}
									}
								}
							}
						}
					}
				}
				if (planes4[3] != null) { // && (planes4[1] != null) && (planes4[2] != null)){
					for (int np = 0; np < planes4[3].length; np++){
						if (planes4[3][np] != null){
							int [] neibs3 = planes4[3][np].getNeibBest();
							if ((neibs3[0] >= 0) && (neibs3[7] >= 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								if (neibs0[2] < 0) {
									int [] neibs1 = planes4[1][neibs3[0]].getNeibBest();
									// from this to right and back
									if ( neibs1[6] < 0) {
										if (goodLink(threshold, nsTiles[0], neibs3[7], 2, neibs3[0], planes, tnSurface )) {
											neibs0[2] = neibs3[0];
											neibs1[6] = neibs3[7];
											num_added++;
											if (dl > 0){
												System.out.println("cutCorners().5: added link "+
														nsTiles[0]+":" + neibs3[7]+"-(2)->"+neibs3[0]+", " + nsTiles[1]+":" + neibs3[0]+"-(6)->" + neibs3[7]+
														" cost = " + getMutualLinkCost(nsTiles[0], neibs3[7], 2, neibs3[0], planes, tnSurface ));												
											}
										}
									}
								}
							}

							if ((neibs3[6] >= 0) && (neibs3[7] >= 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								if (neibs0[4] < 0) {
									int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();
									// from this to down and back
									if (neibs2[0] < 0) {
										if (goodLink(threshold, nsTiles[0], neibs3[7], 4, neibs3[6], planes, tnSurface )) {
											neibs0[4] = neibs3[6];
											neibs2[0] = neibs3[7];
											num_added++;
											if (dl > 0){
												System.out.println("cutCorners().6: added link "+
														nsTiles[0]+":" + neibs3[7]+"-(4)->"+neibs3[6]+", " + nsTiles[2]+":" + neibs3[6]+"-(0)->" + neibs3[7]+
														" cost = " + getMutualLinkCost(nsTiles[0], neibs3[7], 4, neibs3[6], planes, tnSurface ));												
											}
										}
									}
								}
							}

							if ((neibs3[7] >= 0) && (neibs3[6] < 0)){
								int [] neibs0 = planes4[0][neibs3[7]].getNeibBest();
								if (neibs0[4] >= 0) {
									int [] neibs2 = planes4[2][neibs0[4]].getNeibBest();
									// from down to down-right and back
									if (neibs2[2] < 0) {
										if (goodLink(threshold, nsTiles[2], neibs0[4], 2, np, planes, tnSurface )) {
											neibs3[6] = neibs0[4];
											neibs2[2] = np;  //?
											num_added++;
											if (dl > 0){
												System.out.println("cutCorners().7: added link "+
														nsTiles[3]+":" + np+"-(6)->"+neibs0[4]+", " + nsTiles[2]+":" + neibs0[4]+"-(2)->" + np+
														" cost = " + getMutualLinkCost(nsTiles[2], neibs0[4], 2, np, planes, tnSurface ));												
											}
										}
									}
								}
							}
							
							if ((neibs3[6] >= 0) && (neibs3[0] < 0)){
								int [] neibs2 = planes4[2][neibs3[6]].getNeibBest();
								if (neibs2[1] >= 0) {
									int [] neibs1 = planes4[1][neibs2[1]].getNeibBest();
									// from right to down-right and back
									if (neibs1[4] < 0) {
										if (goodLink(threshold, nsTiles[1], neibs2[1], 4, np, planes, tnSurface )) {
											neibs3[0] = neibs2[1]; //?
											neibs1[4] = np;
											num_added++;
											if (dl > 0){
												System.out.println("cutCorners().8: added link "+
														nsTiles[3]+":" + np+"-(0)->"+neibs2[1]+", " + nsTiles[1]+":" + neibs2[1]+"-(4)->" + np+
														" cost = " + getMutualLinkCost(nsTiles[1], neibs2[1], 4, np, planes, tnSurface ));												
											}
										}
									}
								}
							}
						}
					}					
				}
			}
		}
		return num_added;

	}

	public int addFinalLinks (
			final TilePlanes.PlaneData[][] planes,
			final double                   threshold_sq,
			final double                   threshold_corn,
			final double                   threshold_hyp,
			final boolean                  en_sq,
			final boolean                  en_corn,
			final boolean                  en_hyp,
			final int                      debugLevel)
			
	{
		int total_added = 0;
		while (true) {
			while (true) {
				int num_added = 0;
				if (en_sq){
					num_added += fillSquares(
							planes,
							threshold_sq,
							debugLevel,
							dbg_tileX,
							dbg_tileY);

				}
				if (debugLevel > -1) {
					System.out.println("after fillSquares() added "+num_added);
				}
				if (en_corn){
					num_added += cutCorners(
							planes,
							threshold_corn,
							debugLevel,
							dbg_tileX,
							dbg_tileY);
				}
				if (debugLevel > -1) {
					System.out.println("after cutCorners() added (cumulative) "+num_added);
				}
				if (num_added == 0) break;
				total_added += num_added;
			}
			int num_added = 0;
			if (en_hyp){
				num_added += fillHypotenuse(
						planes,
						threshold_hyp,
						debugLevel,
						dbg_tileX,
						dbg_tileY);
			}
			if (debugLevel > -1) {
				System.out.println("after fillHypotenuse() added "+num_added);
			}
			if (num_added == 0) break;
			total_added += num_added;
		}
		return total_added;
	}

}
