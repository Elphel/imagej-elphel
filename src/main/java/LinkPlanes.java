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
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicInteger;


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
	public double     plMaxOverlap; //         =   0.1;  // Maximal overlap between the same supertile planes to merge

		// comparing merge quality for plane pairs
	public double     plCostKrq; //            =   0.8;  // Cost of merge quality weighted
	public double     plCostKrqEq; //          =   0.2;  // Cost of merge quality equal weight
	public double     plCostSin2; //           = 100.0;  // Cost of sin squared between normals
	public double     plCostRdist2; //         =1000.0;  // Cost of squared relative distances       
	
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
		plMaxOverlap =      clt_parameters.plMaxOverlap;
		
		plCostKrq =         clt_parameters.plCostKrq;
		plCostKrqEq =       clt_parameters.plCostKrqEq;
		plCostSin2 =        clt_parameters.plCostSin2;
		plCostRdist2 =      clt_parameters.plCostRdist2;       

		
		
		
		dbg_tileX =         clt_parameters.tileX;
		dbg_tileY =         clt_parameters.tileY;
		this.st =           st;
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
		if (debugLevel > 1){
			System.out.println("planesFit() debug:");
		}
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
					debugLevel - 2); // int       debugLevel)
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
					debugLevel - 2); // int       debugLevel)
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
						System.out.println(prefix+" (fit) this_rq="+this_rq+
								", this_rq_eq="+this_rq_eq+
								" w1="+w1+" w2="+w2+
								" L1="+plane1.getValue()+" L2="+plane2.getValue()+
								" L="+merged_ev+" L_eq="+merged_ev_eq);
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
		}
		if (debugLevel > 0) {
			System.out.println(prefix+": planes DO NOT FIT");
			if (debugLevel > 1){
				System.out.println(prefix+" (do not fit) this_rq="+this_rq+
						", this_rq_eq="+this_rq_eq+
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+
						" L="+merged_ev+" L_eq="+merged_ev_eq);
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
	public double [] getFitQualities(
			TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
			TilePlanes.PlaneData plane2,
			double               merged_ev,    // if NaN will calculate assuming the same supertile
			double               merged_ev_eq, // if NaN will calculate assuming the same supertile
			String prefix,
			int debugLevel)
	{
		if ((plane1 == null) || (plane2 == null)) return null;
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
		
		double sin2 = plane1.getWorldSin2(plane2);
		double rdist2 = plane1.getWorldPlaneRDist2(plane2) + plane2.getWorldPlaneRDist2(plane1);
		double [] costs = {
				this_rq *    plCostKrq,
				this_rq_eq * plCostKrqEq,
				sin2 *       plCostSin2,
				rdist2 *     plCostRdist2};
		double cost = costs[0]+costs[1]+costs[2]+costs[3];
		double [] qualities = {
				this_rq,
				this_rq_eq,
				cost,
				costs[0]/cost,
				costs[1]/cost,
				costs[2]/cost,
				costs[3]/cost};
		
		if (debugLevel > 0){
			System.out.println(prefix+" cost="+cost);
			if (debugLevel > 1){
				System.out.println(prefix+ "cost contributions: "+
						((int)(100*qualities[3]))+"%, "+
						((int)(100*qualities[4]))+"%, "+
						((int)(100*qualities[5]))+"%, "+
						((int)(100*qualities[6]))+"%");
				System.out.println(prefix+
						" this_rq="+this_rq+
						" this_rq=" +   (this_rq) +
						" this_rq_raw="+(this_rq * (w1+w2)) +
						" this_rq_eq="+(this_rq_eq) +
						" this_rq_nofloor="+(this_rq_nofloor) +
						" w1="+w1+" w2="+w2+
						" L1="+plane1.getValue()+" L2="+plane2.getValue()+" L="+merged_ev+
						" L_eq="+merged_ev_eq);
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

		final Thread[] threads = ImageDtt.newThreadArray(st.tileProcessor.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		// Select best symmetrical match, consider only N, NE, E, SE - later opposite ones will be copied
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
//					TilePlanes.PlaneData [] dbg_planes = null; 
					for (int nsTile0 = ai.getAndIncrement(); nsTile0 < nStiles; nsTile0 = ai.getAndIncrement()) {
						int sty0 = nsTile0 / stilesX;  
						int stx0 = nsTile0 % stilesX;
						int dl = ((debugLevel > -1) && (nsTile0 == debug_stile)) ? 1:0;
						if ( planes[nsTile0] != null) {
							if (dl > 0){
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
										if ((stx < stilesX) && (sty < stilesY) && (sty > 0)) {
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
																	dl-1); // debugLevel);
															if (other_plane !=null) { // now always, but may add later
																TilePlanes.PlaneData merged_pd = this_plane.mergePlaneToThis(
																		other_plane, // PlaneData otherPd,
																		1.0,         // double    scale_other,
																		1.0, // double     starWeightPwr,    // Use this power of tile weight when calculating connection cost																		
																		false,       // boolean   ignore_weights,
																		true, // boolean   sum_weights,
																		plPreferDisparity, 
																		dl-1); // int       debugLevel)

																if (merged_pd !=null) { // now always, but may add later
																	///															merged_pd.scaleWeight(0.5);
																	this_plane.setNeibMatch(dir, np, merged_pd.getValue()); // smallest eigenValue
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
																	this_plane.setNeibMatchEq(dir, np, merged_pd.getValue()); // smallest eigenValue
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
							if (dl > 0){
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
						if (dl>0) {
							System.out.println("matchPlanes() nsTile0="+nsTile0);
						}
						if ( planes[nsTile0] != null) {
							for (int np0 = 0; np0 < planes[nsTile0].length; np0++){ // nu
								TilePlanes.PlaneData this_plane = planes[nsTile0][np0];
								if (this_plane != null) {
									for (int dir = 4; dir < 8; dir++){ // other half - copy from opposite
										int stx = stx0 + dirsYX[dir][1];
										int sty = sty0 + dirsYX[dir][0];
										if ((sty < stilesY) && (sty > 0) && (stx > 0)) {
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
													}

												}
											}
										}
									}
								}
							}
						}
						if (dl>0) {
							System.out.println("matchPlanes() nsTile0="+nsTile0);
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
										double [] merge_ev = planes[nsTile0][np0].getMergedValue(dir);
										double [] merge_ev_eq = planes[nsTile0][np0].getMergedValueEq(dir);
										if (	(merge_ev != null) &&
												(merge_ev_eq != null)) {
											int np_min = SuperTiles.LOWEST_PLANE(merge_ev.length);
											for (int np = np_min; np < merge_ev.length; np++){
												if (	(planes[nsTile][np] != null) &&
														!Double.isNaN(merge_ev[np])) {
													String prefix = "filterNeighborPlanes() nsTile0="+nsTile0+" np0="+np0+" dir="+dir+" nsTile="+nsTile+" np="+np;
													if (planesFit(
															planes[nsTile0][np0], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
															planes[nsTile][np],   // 			TilePlanes.PlaneData plane2,
															merge_ev[np],         // double               merged_ev,    // if NaN will calculate assuming the same supertile
															merge_ev_eq[np],      //double               merged_ev_eq, // if NaN will calculate assuming the same supertile
															prefix,               // String prefix,
															dl)                   // int debugLevel)
															){
														planes[nsTile0][np0].setMergedValid(dir, np, true, planes[nsTile].length);
														planes[nsTile][np].setMergedValid((dir + 4) %8, np0, true, planes[nsTile0].length);
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
	 * Find mutual links between multi-layer planes for supertiles. requires that for each plane there are calculated smalles eigenvalues
	 * for merging with each plane for each of 8 neighbors 
	 * @param debugLevel debug level
	 * @param dbg_X debug supertile X coordinate
	 * @param dbg_Y debug supertile Y coordinate
	 */
	public double [][]  selectNeighborPlanesMutual(
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
		final double [][][] all_quality_stats = new double  [numThreads][2][4]; // contributions of all [0] and winners [2], 4 parameters 
		
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
											boolean [] merge_valid = planes[nsTile0][np0].getMergedValid(dir);
											qualities[np0] = new double[ merge_ev.length][];
											
											if (!this_matched[np0] &&(merge_valid != null)) {
												for (int np = np_min; np < merge_ev.length; np++){
													if (!other_matched[np] && merge_valid[np]) {
														String prefix = "selectNeighborPlanesMutual(): nsTile0="+nsTile0+":"+np0+", nsTile="+nsTile+":"+np;
														qualities[np0][np] = getFitQualities( // {this_rq, this_rq_eq}; 
																planes[nsTile0][np0], //TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																planes[nsTile][np], //TilePlanes.PlaneData plane2,
																merge_ev[np], // double               merged_ev,    // if NaN will calculate assuming the same supertile
																merge_ev_eq[np], // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
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
		double [][] quality_stats = new double  [2][4]; // contributions of all [0] and winners [2], 4 parameters
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
			System.out.println("    weighted quality: "+quality_stats[0][0]+" (all), "+quality_stats[1][0]+" (winners)");
			System.out.println("non-weighted quality: "+quality_stats[0][1]+" (all), "+quality_stats[1][1]+" (winners)");
			System.out.println("                sin2: "+quality_stats[0][2]+" (all), "+quality_stats[1][2]+" (winners)");
			System.out.println("              rdist2: "+quality_stats[0][3]+" (all), "+quality_stats[1][3]+" (winners)");
		}
		return quality_stats;

		
	}
	
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
							sdfa_instance.showArrays(plane_strengths, 2 * superTileSize, 2* superTileSize, true, "planes_"+nsTile);
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

	public  boolean [][][] overlapSameTileCandidates(
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
							sdfa_instance.showArrays(plane_strengths, 2 * superTileSize, 2* superTileSize, true, "planes_"+nsTile);
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
			};
		}		      
		ImageDtt.startAndJoin(threads);
		return overlap_merge_candidates;
	}
	
	public boolean [][] mergeSameTileEvaluate(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final boolean [][][] plane_nooverlaps,			
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
										Double.NaN, // calculate double               merged_ev,    // if NaN will calculate assuming the same supertile
										Double.NaN, // calculate double               merged_ev_eq, // if NaN will calculate assuming the same supertile
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
	
	public int [][][] extractMergeSameTileGroups(
			final TilePlanes.PlaneData [][] planes,
			final int [][][] merge_candidates,
			final boolean [][][] plane_nooverlaps,			
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
		
		
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				public void run() {
					for (int nsTile = ai.getAndIncrement(); nsTile < nStiles; nsTile = ai.getAndIncrement()) if ( merge_candidates[nsTile] != null) {
						boolean [][] merge_pairs = new boolean [planes[nsTile].length][planes[nsTile].length]; 
						int dl = ((debugLevel > 0) && (nsTile == debug_stile)) ? 2: ((debugLevel > 1) ? 1:0);
						if (dl > 1){
							System.out.println("extractMergeSameTileGroups(): nsTile="+nsTile);
						}
						HashSet<Integer> yet_to_merge = new HashSet<Integer>(); 
						for (int pair = 0; pair < merge_candidates[nsTile].length; pair ++){
							int np1 =  merge_candidates[nsTile][pair][0];
							int np2 =  merge_candidates[nsTile][pair][1];
							if ((plane_nooverlaps == null) || (plane_nooverlaps[nsTile] == null) || plane_nooverlaps[nsTile][np1][np2]) {
								yet_to_merge.add(np1);
								yet_to_merge.add(np2);
								String prefix = "mergeSameTileEvaluate() pair="+pair+" nsTile="+nsTile+" np1="+np1+" np2="+np2;
								if (planesFit(
										planes[nsTile][np1], // TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
										planes[nsTile][np2], // 			TilePlanes.PlaneData plane2,
										Double.NaN, // calculate double               merged_ev,    // if NaN will calculate assuming the same supertile
										Double.NaN, // calculate double               merged_ev_eq, // if NaN will calculate assuming the same supertile
										prefix, // String prefix,
										dl) // int debugLevel)
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
																planes[nsTile][np1], //TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																planes[nsTile][np2], //TilePlanes.PlaneData plane2,
																Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
																Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
																prefix, // String prefix,
																dl); // int debugLevel)
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
																if ((plane_nooverlaps != null) &&
																		(plane_nooverlaps[nsTile] != null) &&
																		!plane_nooverlaps[nsTile][np][np1]) {
																	nooverlap = false;
																	break;
																}
															}
															if (nooverlap){
																String prefix = "extractMergeSameTileGroups().2 nsTile="+nsTile+" np="+np;
																double [] qualities = getFitQualities( // {this_rq, this_rq_eq}; 
																		merged_pd, //TilePlanes.PlaneData plane1, // should belong to the same supertile (or be converted for one)
																		planes[nsTile][np], //TilePlanes.PlaneData plane2,
																		Double.NaN, // double               merged_ev,    // if NaN will calculate assuming the same supertile
																		Double.NaN, // double               merged_ev_eq, // if NaN will calculate assuming the same supertile
																		prefix, // String prefix,
																		dl); // int debugLevel)
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
												System.out.println("====== BUG: extractMergeSameTileGroups() nsTile = "+nsTile+" : found incomplete group: "+group+""
														+ "best_pair = null");
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
