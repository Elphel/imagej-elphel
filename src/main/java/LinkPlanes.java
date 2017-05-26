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
/*
				clt_parameters.plWorstWorsening, // final double worst_worsening,
				clt_parameters.plWorstWorsening2,// final double worst_worsening2  Worst case worsening for thin planes,
				clt_parameters.plWorstEq,        // final double worstEq,        // Worst case worsening after merge with equal weights
				clt_parameters.plWorstEq2,       // final double worstEq2,       // Worst case worsening for thin planes with equal weights
				clt_parameters.plWeakWorsening,  // final double worst_worsening,
				clt_parameters.plOKMergeEigen,   // final double okMergeEigen, f result of the merged planes is below, OK to use thin planes (higher) threshold
				clt_parameters.plMaxWorldSin2,   // final double maxWorldSin2,
				clt_parameters.plDispNorm,
				clt_parameters.plMaxEigen,
				clt_parameters.plEigenFloor,    // final double eigenFloor, // Add to eigenvalues of each participating plane and result to validate connections 
				clt_parameters.plMinStrength,
				0, // final int debugLevel)
				clt_parameters.tileX,
				clt_parameters.tileY);

 */

public class LinkPlanes {
	SuperTiles st;
	public boolean    plPreferDisparity; //    =   false;// Always start with disparity-most axis (false - lowest eigenvalue)
	public double     plDispNorm; //           =   3.0;  // Normalize disparities to the average if above (now only for eigenvalue comparison)
	public double     plMinStrength; //        =   0.1;  // Minimal total strength of a plane 
	public double     plMaxEigen; //           =   0.05; // Maximal eigenvalue of a plane 
	public double     plEigenFloor; //         =   0.01; // Add to eigenvalues of each participating plane and result to validate connections 

	public double     plWorstWorsening; //     =   2.0;  // Worst case worsening after merge
	public double     plWorstWorsening2; //    =   5.0;  // Worst case worsening for thin planes
	public double     plWorstEq; //            =   1.0;  // Worst case worsening after merge with equal weights
	public double     plWorstEq2; //           =   2.0;  // Worst case worsening for thin planes with equal weights

	public double     plOKMergeEigen; //       =   0.03; // If result of the merged planes is below, OK to use thin planes (higher) threshold 
	public double     plMaxWorldSin2; //       =   0.1;  // Maximal sine squared of the world angle between planes to merge. Set to >= 1.0 to disable

	public double     plWeakWorsening; //      =   1.0;  // Relax merge requirements for weaker planes
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

		plDispNorm =        clt_parameters.plDispNorm;
		plMaxEigen =        clt_parameters.plMaxEigen;
		plEigenFloor =      clt_parameters.plEigenFloor;    // final double eigenFloor, // Add to eigenvalues of each participating plane and result to validate connections 
		plMinStrength =     clt_parameters.plMinStrength;
		dbg_tileX =         clt_parameters.tileX;
		dbg_tileY =         clt_parameters.tileY;
	}
	public boolean planesFit(
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			String prefix,
			int debugLevel)
	{
		if ((plane1 == null) || (plane2 == null)) return false;
		TilePlanes.PlaneData merged_pd = null;
		TilePlanes.PlaneData merged_pd_eq = null;
		if (plane1.getWeight() < plMinStrength) {
			if (debugLevel > 1)	System.out.println(prefix+" plane1 is too weak ("+plane1.getWeight()+" < plMinStrength="+plMinStrength+")");
			return false;
		}
		if (plane2.getWeight() < plMinStrength) {
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
		
		if (Double.isNaN(merged_ev)) {
			merged_pd = plane1.mergePlaneToThis(
					plane2, // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					false,       // boolean   ignore_weights,
					true, // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 1); // int       debugLevel)
			merged_ev = merged_pd.getValue();
		}
		if (Double.isNaN(merged_ev_eq)) {
			merged_pd_eq = plane1.mergePlaneToThis(
					plane2,      // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					true,        // boolean   ignore_weights,
					true,        // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 1); // int       debugLevel)
			merged_ev_eq = merged_pd_eq.getValue();
		}
		double w1 = plane1.getWeight();
		double w2 = plane2.getWeight();
		double this_rq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue() , // double L2,
				merged_ev, // double L,
				w1, // double w1,
				w2, // double w2)
				plEigenFloor);//			double eigen_floor)
		double this_rq_norm = this_rq;
		if ((w1 + w2) < plWeakWorsening) this_rq_norm *= (w1 + w2) / plWeakWorsening; // forgive more for weak planes
		double this_rq_eq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue() , // double L2,
				merged_ev_eq, // double L,
				1.0, // double w1,
				1.0, // double w2)
				plEigenFloor);//			double eigen_floor)
		
		double this_rq_eq_norm = this_rq_eq;
		if ((w1 + w2) < plWeakWorsening) this_rq_eq_norm *= (w1 + w2) / plWeakWorsening; // forgive more for weak planes 
		
		if ((this_rq_norm <= plWorstWorsening) ||
				((merged_ev <= plOKMergeEigen) && (this_rq_norm <= plWorstWorsening2)) || // use higher threshold
				(this_rq_eq_norm <= plWorstEq) ||
				((merged_ev_eq <= plOKMergeEigen) && (this_rq_eq_norm <= plWorstEq2)) // use higher threshold
				) {
			if ((plMaxWorldSin2 >= 1.0) || (plane1.getWorldSin2(plane2) <= plMaxWorldSin2)) {
				if (debugLevel > 0){
					System.out.print(prefix+": planes FIT");
					if (this_rq_norm <= plWorstWorsening)
						System.out.print(" (this_rq_norm="+this_rq_norm+"  <= plWorstWorsening="+plWorstWorsening+")");
					if ((merged_ev <= plOKMergeEigen) && (this_rq_norm <= plWorstWorsening2))
						System.out.print(" merged_ev="+merged_ev+" <= plOKMergeEigen="+plOKMergeEigen+") && (this_rq_norm="+this_rq_norm+
								" <= plWorstWorsening2="+plWorstWorsening2+")");
					if (this_rq_eq_norm <= plWorstEq)
						System.out.print(" this_rq_eq_norm="+this_rq_eq_norm+" <= plWorstEq="+plWorstEq);
					if ((merged_ev_eq <= plOKMergeEigen) && (this_rq_eq_norm <= plWorstEq2))
						System.out.print(" ((merged_ev_eq="+merged_ev_eq+" <= plOKMergeEigen) && (this_rq_eq_norm="+
								this_rq_eq_norm+" <= plWorstEq2="+plWorstEq2+")");
					System.out.println();
					if (debugLevel > 1){
						System.out.println(prefix+" this_rq="+this_rq+
								", this_rq_eq="+this_rq_eq+
								" w1="+w1+" w2="+w2+
								" L1="+plane1.getValue()+" L2="+plane2.getValue()+
								" L="+merged_ev+" L_eq="+merged_ev_eq);
						System.out.println(prefix+", world sin2 ="+
								plane1.getWorldSin2(plane2));
						System.out.println(prefix+
								", world dist this="+ Math.sqrt(plane1.getWorldPlaneDist2(plane2))+
								", world dist other="+Math.sqrt(plane2.getWorldPlaneDist2(plane1))+
								", world dist sum="+Math.sqrt(plane1.getWorldPlaneDist2(plane2)+
										plane2.getWorldPlaneDist2(plane1)));
					}
				}
				return true;
			}
		}
		if (debugLevel > 0) {
			System.out.print(prefix+": planes DO NOT FIT");
			if (debugLevel > 1){
				System.out.println(prefix+" this_rq="+this_rq+
						", this_rq_eq="+this_rq_eq+
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+
						" L="+merged_ev+" L_eq="+merged_ev_eq);
				System.out.println(prefix+", world sin2 ="+
						plane1.getWorldSin2(plane2));
				System.out.println(prefix+
						", world dist this="+ Math.sqrt(plane1.getWorldPlaneDist2(plane2))+
						", world dist other="+Math.sqrt(plane2.getWorldPlaneDist2(plane1))+
						", world dist sum="+Math.sqrt(plane1.getWorldPlaneDist2(plane2)+
								plane2.getWorldPlaneDist2(plane1)));
			}
		}
		return false;
	}
	
	public double getFitQuality(
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			String prefix,
			int debugLevel)
	{
		if ((plane1 == null) || (plane2 == null)) return 0.0;
		TilePlanes.PlaneData merged_pd = null;
		TilePlanes.PlaneData merged_pd_eq = null;
		if (Double.isNaN(merged_ev)) {
			merged_pd = plane1.mergePlaneToThis(
					plane2, // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					false,       // boolean   ignore_weights,
					true, // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 1); // int       debugLevel)
			merged_ev = merged_pd.getValue();
		}
		if (Double.isNaN(merged_ev_eq)) {
			merged_pd_eq = plane1.mergePlaneToThis(
					plane2,      // PlaneData otherPd,
					1.0,         // double    scale_other,
					1.0,         // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
					true,        // boolean   ignore_weights,
					true,        // boolean   sum_weights,
					plPreferDisparity, 
					debugLevel - 1); // int       debugLevel)
			merged_ev_eq = merged_pd_eq.getValue();
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
		double this_rq_nofloor = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue(), // double L2,
				merged_ev, // double L,
				w1, // double w1,
				w2, // double w2)
				0); // eigenFloor);//			double eigen_floor)
		double this_rq_eq = mergeRQuality(
				plane1.getValue(), // double L1,
				plane2.getValue(), // double L2,
				merged_ev_eq, // double L,
				1.0, // double w1,
				1.0, // double w2)
				plEigenFloor);//			double eigen_floor)
		this_rq /= (w1 + w2); // for comparison reduce this value for stronger planes
		
		if (debugLevel > 0){
			System.out.println(prefix+", this_rq="+this_rq);
			if (debugLevel > 1){
				System.out.println(prefix+", this_rq="+this_rq+
						" this_rq_raw="+(this_rq * (w1+w2)) +
						" this_rq_eq="+(this_rq_eq) +
						" this_rq_nofloor="+(this_rq_nofloor) +
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+" L="+merged_ev+
						" L_eq="+merged_ev_eq);
				System.out.println(prefix + ", world sin2 =" + plane1.getWorldSin2(plane2));
				System.out.println(prefix+
						", world dist this="+ Math.sqrt(plane1.getWorldPlaneDist2(plane2))+
						", world dist other="+Math.sqrt(plane2.getWorldPlaneDist2(plane1))+
						", world dist sum="+Math.sqrt(plane1.getWorldPlaneDist2(plane2)+
								plane2.getWorldPlaneDist2(plane1)));
			}
		}
		return this_rq; // TODO: add modes to select what is output
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
	public double mergeRQuality(
			double L1_in,
			double L2_in,
			double L_in,
			double w1,
			double w2,
			double eigen_floor)
	{
		double L1 = L1_in + eigen_floor;
		double L2 = L2_in + eigen_floor;
		double L =  L_in + eigen_floor;
		//				double Lav = Math.sqrt((L1*L1*w1 + L2*L2*w2)/(w1+w2));
		double Lav = (L1*w1 + L2*w2)/(w1+w2);
		///				double wors = (L - Lav)*(w1+w2)*(w1+w2) /(Lav*w1*w2);
		///				double rquality = (L - Lav)*(w1+w2) /(Lav*w1*w2); // ==wors/(w1+w2) to amplify stronger planes
		double rquality =  (L - Lav)*(w1+w2)*(w1+w2) /(Lav*w1*w2); 
		return rquality;
	}
	
}
