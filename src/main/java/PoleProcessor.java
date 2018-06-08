/**
 ** PoleProcessor - handling poles like street lights and other vertical objects
 ** common in artificial environments
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  PoleProcessor.java is free software: you can redistribute it and/or modify
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
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;


public class PoleProcessor {
	BiCamDSI biCamDSI;
	ArrayList<PoleCluster> pole_clusters = new ArrayList<PoleCluster>();
	int tilesX;
	int tilesY;

	class PoleCluster{
		double mean_disp = 0.0;
		double mean_x =    0.0;
		double mean_y =    0.0;
		double strength =  0.0;
		ArrayList <Integer> tiles = new ArrayList<Integer>();
		Rectangle bBox = new Rectangle();
		Rectangle eBox; // not calculated, assigned by caller
		int    layer =     -1;
		double  getMeanDisp() {return mean_disp;}
		double  getStrength() {return strength;} // to get average strength - need to divide by getNumTiles()
		double  getMeanX()    {return mean_x;}
		double  getMeanY()    {return mean_y;}
		int     getNumTiles() {return tiles.size();}
		ArrayList<Integer>  getTiles() {return tiles;}
		public Rectangle getBBox() {return bBox;}
		public Rectangle getEBox() {return eBox;}
		public void trimEBoxBottom(int height) {eBox.height = height;}

		public void setEBox(
				int ext_left,
				int ext_right,
				int ext_up,
				int ext_down)
		{
			this.eBox = new Rectangle(
					bBox.x - ext_left,
					bBox.y - ext_up,
					bBox.width +  ext_right + ext_left,
					bBox.height + ext_up +   ext_down);
			if (eBox.x < 0) {
				eBox.width+=eBox.x;
				eBox.x = 0;
			}
			if (eBox.y < 0) {
				eBox.height +=eBox.y;
				eBox.y = 0;
			}
			if ((eBox.x + eBox.width) >= tilesX) {
				eBox.width = tilesX - eBox.x;
			}
			if ((eBox.y + eBox.height) >= tilesY) {
				eBox.height = tilesY - eBox.y;
			}
		}
		public boolean intersects(PoleCluster cluster) {
			return this.eBox.intersects(cluster.eBox);
		}

		public void setLayer(int layer) {this.layer = layer;}

		public int getLayer() {return layer;}

		void addCluster(
				double [][] norm_ds,
				PoleCluster cluster)
		{
			for (int otherTile:cluster.getTiles()) {
				addTile(norm_ds, otherTile);
			}
			// merge eBox-es if any
			if ((eBox != null) && (cluster.eBox !=null)) {
				eBox.add(cluster.eBox);
			}
		}


		void addTile(
				double [][] norm_ds,
				int nTile)
		{
			if (tiles.indexOf(nTile) >=0) {
				System.out.println("*** BUG in PoleCluster.addTile("+nTile+") - tile already exists here ***");
				return; //
			}
			if (norm_ds[1][nTile] <= 0) {
				System.out.println("*** BUG in PoleCluster.addTile("+nTile+") - strength= "+norm_ds[1][nTile]+" ***");
				return; //
			}
			TileNeibs  tnImage = biCamDSI.tnImage;
			int tileX = nTile % tnImage.sizeX;
			int tileY = nTile / tnImage.sizeX;
			mean_disp = (mean_disp * strength + norm_ds[0][nTile] * norm_ds[1][nTile])/(strength + norm_ds[1][nTile]);
			mean_x =    (mean_x *    strength + tileX *             norm_ds[1][nTile])/(strength + norm_ds[1][nTile]);
			mean_y =    (mean_y *    strength + tileY *             norm_ds[1][nTile])/(strength + norm_ds[1][nTile]);
			strength += norm_ds[1][nTile];
			if (tiles.isEmpty() || ( tileX < bBox.x))  bBox.x = tileX;
			if (tiles.isEmpty() || ((tileX - bBox.x) > (bBox.width -1)))  bBox.width = tileX - bBox.x + 1;
			if (tiles.isEmpty() || ( tileY < bBox.y))  bBox.y = tileY;
			if (tiles.isEmpty() || ((tileY - bBox.y) > (bBox.height -1)))  bBox.height = tileY - bBox.y + 1;
			tiles.add(nTile);
		}
		// remove would be trikier - will need to recalculate bBox
		void removeTile(
				double [][] norm_ds,
				int nTile)
		{
			int indx = tiles.indexOf(nTile);
			if (indx < 0) {
				System.out.println("*** BUG in PoleCluster.removeTile("+nTile+") - tile does not belong here ***");
				return; //
			}
			TileNeibs  tnImage = biCamDSI.tnImage;
			int tileX = nTile % tnImage.sizeX;
			int tileY = nTile / tnImage.sizeX;
			mean_disp = (mean_disp * strength - norm_ds[0][nTile] * norm_ds[1][nTile])/(strength - norm_ds[1][nTile]);
			mean_x =    (mean_x *    strength - tileX *             norm_ds[1][nTile])/(strength - norm_ds[1][nTile]);
			mean_y =    (mean_y *    strength - tileY *             norm_ds[1][nTile])/(strength - norm_ds[1][nTile]);
			strength -= norm_ds[1][nTile];
			tiles.remove(indx);
			// re-calculate bounding box
			for (int i = 0; i < tiles.size(); i++ ) {
				int nt = tiles.get(i);
				tileX = nt % tnImage.sizeX;
				tileY = nt / tnImage.sizeX;
				if ((i==0) || ( tileX < bBox.x))  bBox.x = tileX;
				if ((i==0) || ((tileX - bBox.x) > (bBox.width -1)))  bBox.width = tileX - bBox.x + 1;
				if ((i==0) || ( tileY < bBox.y))  bBox.y = tileX;
				if ((i==0) || ((tileY - bBox.y) > (bBox.height -1)))  bBox.height = tileY - bBox.y + 1;
			}
		}


	}



	public PoleProcessor (
			BiCamDSI   biCamDSI,
			int tilesX,
			int tilesY)
	{
		this.biCamDSI = biCamDSI;
		this.tilesX = tilesX;
		this.tilesY = tilesY;
	}

	public double [][] conditionDisparityStrength(
			final BiScan biScan,
			final double     trusted_strength, // trusted correlation strength
			final double     strength_rfloor,
			final double     strength_pow
			) {
		final double     strength_floor = trusted_strength * strength_rfloor;
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final double [][] cond_ds = new double [2] [num_tiles];
		final double [][] disparity_strength = biScan.getDisparityStrength(
				false, // only_strong,
				false, // only_trusted,
				true); // only_enabled);
//		final AtomicInteger ai_num_seeds = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) {
						if (!Double.isNaN(disparity_strength[0][nTile]) &&(disparity_strength[1][nTile] > strength_floor)){
							double w = disparity_strength[1][nTile] - strength_floor;
							if (strength_pow != 1.0) {
								w = Math.pow(w, strength_pow);
							}
							cond_ds[0][nTile] = disparity_strength[0][nTile];
							cond_ds[1][nTile] = w;

						} else {
							cond_ds[0][nTile] = Double.NaN;
							cond_ds[1][nTile] = 0.0;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return cond_ds;
	}




	public boolean [] findSeeds(
			final double [][] norm_ds,
			final double      min_strength,     // after subtracting floor and possible applying pow
			final int         seed_down,
			final double      seed_aover,
			final double      seed_rover,
			final double      max_disparity
			) {
		final TileNeibs  tnImage = biCamDSI.tnImage;
		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final boolean [] seeds = new boolean [num_tiles];
//		final AtomicInteger ai_num_seeds = new AtomicInteger(0);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < num_tiles; nTile = ai.getAndIncrement()) if (norm_ds[1][nTile] > min_strength){
						double d_head = norm_ds[0][nTile];
						if ((d_head <= max_disparity) && (d_head > seed_aover)) {
							double sw = 0.0;
							double swd = 0.0;
							int nTile1 = nTile; // tnImage.getNeibIndex(indx, TileNeibs.DIR_DOWN);
							for (int i = 0; i < seed_down; i++) {
								nTile1 = tnImage.getNeibIndex(nTile1, TileNeibs.DIR_DOWN);
								if ((nTile1 >= 0) && (norm_ds[1][nTile1] > 0.0) && !Double.isNaN(norm_ds[0][nTile1])) {
									double w = norm_ds[1][nTile1];
									sw += w;
									swd += w * norm_ds[0][nTile1];
								}
							}
							if (sw > 0.0){
								swd /= sw; // weighted average of the disparities right below
								double threshold = swd + seed_aover + seed_rover * d_head;
								if (d_head >= threshold) {
									seeds[nTile] = true;
//									ai_num_seeds.getAndIncrement();
								}
							}
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		return seeds;
	}

	public ArrayList<PoleCluster> initPoleClusters(
			final double max_dx,
			final double max_dy,
			final double max_dd,
			boolean []   bseeds,
			double [][]  norm_ds,
			final int    debugLevel)
	{
		final TileNeibs  tnImage = biCamDSI.tnImage;
		ArrayList<PoleCluster> clusters = new ArrayList<PoleCluster>();
		for (int nTile = 0; nTile < bseeds.length; nTile++) if (bseeds[nTile]){
			double this_disp = norm_ds[0][nTile];
			int tileX = nTile % tnImage.sizeX;
			int tileY = nTile / tnImage.sizeX;
			label_gotit: {
				for (PoleCluster cluster: clusters) {
					if (
							!(Math.abs(cluster.getMeanDisp() - this_disp) > max_dd) && // to allow Double.NaN
							!(Math.abs(cluster.getMeanX() - tileX) > max_dx) &&
							!(Math.abs(cluster.getMeanY() - tileY) > max_dy)) {
						cluster.addTile(norm_ds, nTile);
						break label_gotit;
					}
				}
				// did not find, need to create a new cluster
				PoleCluster cluster = new PoleCluster();
				cluster.addTile(norm_ds, nTile);
				clusters.add(cluster);
			}
		}
		sortClustersByDisparity(clusters);
		return clusters;
	}

	public void extendClusters(
			int ext_left,
			int ext_right,
			int ext_up,
			int ext_down,
			ArrayList<PoleCluster> clusters)
	{
		for (PoleCluster cluster: clusters) {
			// Assign extended bounding box to each cluster
			cluster.setEBox(
					ext_left,
					ext_right,
					ext_up,
					ext_down);
		}
	}


	public void cutClusterBottoms(
			final double [][] norm_ds,
			final double      pedestal_strength,     // bottom tiles should be at least this strong (normalized strength)
			final double      pedestal_disp_over,         // bottom tiles should have disparity at least by this above cluster
			final ArrayList<PoleCluster> clusters)
	{

		// can be multithreaded
		final TileNeibs  tnImage = biCamDSI.tnImage;
//		final int num_tiles = tnImage.sizeX * tnImage.sizeY;
		final Thread[] threads = ImageDtt.newThreadArray(biCamDSI.threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final int num_clust = clusters.size();
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nClust = ai.getAndIncrement(); nClust < num_clust; nClust = ai.getAndIncrement()) {
						PoleCluster cluster = clusters.get(nClust);
						Rectangle bBox = cluster.getBBox(); // never extends beyond ?
						Rectangle eBox = cluster.getEBox();
						int nTile = eBox.x + tilesX * eBox.y;
						double disp_marg = cluster.getMeanDisp() + pedestal_disp_over;
						for (int dy = bBox.y - eBox.y + bBox.height; dy < eBox.height; dy++) {
							// label
							label_still_some_far_weak:{
							for ( int dx = 0; dx < eBox.width;  dx++) {
								int nTile1 = tnImage.getNeibIndex(nTile, dx, dy); // should always be inside image
								if (    (norm_ds[1][nTile1] < pedestal_strength) ||
										(norm_ds[0][nTile1] < disp_marg)) {
									break label_still_some_far_weak;
								}
							}
							// found pedestal - trim eBox
							cluster.trimEBoxBottom(dy);
							break;
						}

						}


					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}



	public int assignClustersToLayers(
			ArrayList<PoleCluster> clusters)
	{

		ArrayList<ArrayList<PoleCluster>> layer_lists = new ArrayList<ArrayList<PoleCluster>>();

		for (PoleCluster cluster: clusters) {
			// Assign extended bounding box to each cluster
			// see if this cluster fits to any of the existing layers
			label_layers:{
				for (int nlayer = 0; nlayer < layer_lists.size(); nlayer++) {
					ArrayList<PoleCluster> layer = layer_lists.get(nlayer);
					// see if current cluster intersects with any on this layer
					label_clusters: {
						for (PoleCluster other_cluster: layer){
							if (cluster.intersects(other_cluster)) {
								break label_clusters; // this layer won't work
							}
						}
						// Layer OK
						cluster.setLayer(nlayer);
						layer.add(cluster);
						break label_layers;
					}
				}
				// need new layer
				cluster.setLayer(layer_lists.size());
				ArrayList<PoleCluster> new_layer = new ArrayList<PoleCluster>();
				new_layer.add(cluster);
				layer_lists.add(new_layer);
			}
		}
		return layer_lists.size();
	}

	public int mergeOverlappingClusters(
			boolean     extended,
			double      disp_tolerance,
			double [][] norm_ds,
			ArrayList<PoleCluster> clusters,
			int debugLevel)
	{
		int num_merge = 0;
		for (int nCluster = 0; nCluster < clusters.size(); nCluster++) {
			PoleCluster cluster = clusters.get(nCluster);
			double this_disp = cluster.getMeanDisp();
			Rectangle this_box = extended?cluster.getEBox():cluster.getBBox();
			for (int nOther = 0; nOther < nCluster; nOther++) {
				PoleCluster other_cluster = clusters.get(nOther);
				double other_disp = other_cluster.getMeanDisp();
				Rectangle other_box = extended?other_cluster.getEBox():other_cluster.getBBox();
				if (this_box.intersects(other_box) && (Math.abs(other_disp - this_disp) < disp_tolerance)) {
					if (debugLevel > -2) {
						System.out.println("Merging cluster  " +nCluster+ " to cluster " +nOther);
						System.out.println("Merging cluster   ("+this_box.x+","+this_box.y+","+this_box.width+","+this_box.height+") disp="+this_disp+
								" to cluster ("+other_box.x+","+other_box.y+","+other_box.width+","+other_box.height+") disp="+other_disp);
					}
					other_cluster.addCluster(norm_ds, cluster);
					clusters.remove(nCluster);
					nCluster--; // so when incremented by for() it will stay the same
					num_merge++;
					if (debugLevel > -2) {
						other_disp = other_cluster.getMeanDisp();
						other_box =  extended ? other_cluster.getEBox(): other_cluster.getBBox();
						System.out.println("Combined cluster: ("+other_box.x+","+other_box.y+","+other_box.width+","+other_box.height+") disp="+other_disp);
					}
					break;
				}
			}
		}
		return num_merge;
	}



	double [][] dbgClusterLayers( // layer and eBox should be set
			boolean show_bbox,
			boolean show_ebox,
			ArrayList<PoleCluster> clusters)
	{
		final TileNeibs  tnImage = biCamDSI.tnImage;
		int tilesX = tnImage.sizeX;
		int tilesY = tnImage.sizeY;
		int num_tiles =  tilesY*tilesX;
		int num_layers = 0;
		for (PoleCluster cluster: clusters){
			int cl = cluster.getLayer();
			if (cl > num_layers) num_layers = cl;
		}
		num_layers++;
		double [][] dbg_layers = new double [num_layers][num_tiles];
		for (int i = 0; i < num_layers; i++) {
			for (int j=0; j < num_tiles; j++) {
				dbg_layers[i][j] = Double.NaN;
			}
		}
		for (PoleCluster cluster: clusters){
			int layer = cluster.getLayer();
			double disp = cluster.getMeanDisp();
			Rectangle [] rectangles = {show_bbox?cluster.getBBox():null,show_ebox?cluster.getEBox():null};
			for (int nbox = 0; nbox< rectangles.length; nbox++) {
				Rectangle box = new Rectangle(rectangles[nbox]);
				if (box.x < 0) {
					box.width+=box.x;
					box.x = 0;
				}
				if (box.y < 0) {
					box.height +=box.y;
					box.y = 0;
				}
				if ((box.x + box.width) >= tilesX) {
					box.width = tilesX - box.x;
				}
				if ((box.y + box.height) >= tilesY) {
					box.height = tilesY - box.y;
				}
				int nTile = box.y * tilesX + box.x;
				for (int dx = 0; dx < box.width;dx++) {
					if (tnImage.getNeibIndex(nTile, dx, 0) < 0 ) {
						System.out.println("dbgClusterLayers() bug 1: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dx="+dx);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, 0)] = disp; // should never be out of bounds
					if (tnImage.getNeibIndex(nTile, dx, box.height -1) < 0 ) {
						System.out.println("dbgClusterLayers() bug 2: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dx="+dx);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, box.height -1)] = disp; // should never be out of bounds
				}
				for (int dy = 1; dy < box.height - 1; dy++) {
					if (tnImage.getNeibIndex(nTile, 0, dy) < 0 ) {
						System.out.println("dbgClusterLayers() bug 3: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dy="+dy);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, 0,            dy)] = disp; // should never be out of bounds
					if (tnImage.getNeibIndex(nTile, box.width -1, dy) < 0 ) {
						System.out.println("dbgClusterLayers() bug 4: box.x ="+box.x+", box.width="+box.width+"box.y ="+box.y+", box.height="+box.height+"dy="+dy);
						continue;
					}
					dbg_layers[layer][tnImage.getNeibIndex(nTile, box.width -1, dy)] = disp; // should never be out of bounds
				}
			}
		}
		return dbg_layers;
	}


	double [][] getClusterLayers( // layer and eBox should be set
			ArrayList<PoleCluster> clusters)
	{
		final TileNeibs  tnImage = biCamDSI.tnImage;
		int tilesX = tnImage.sizeX;
		int tilesY = tnImage.sizeY;
		int num_tiles =  tilesY*tilesX;
		int num_layers = 0;
		for (PoleCluster cluster: clusters){
			int cl = cluster.getLayer();
			if (cl > num_layers) num_layers = cl;
		}
		num_layers++;
		double [][] dbg_layers = new double [num_layers][num_tiles];
		for (int i = 0; i < num_layers; i++) {
			for (int j=0; j < num_tiles; j++) {
				dbg_layers[i][j] = Double.NaN;
			}
		}
		for (PoleCluster cluster: clusters){
			int layer = cluster.getLayer();
			double disp = cluster.getMeanDisp();
			Rectangle box = new Rectangle(cluster.getEBox());
			if (box.x < 0) {
				box.width+=box.x;
				box.x = 0;
			}
			if (box.y < 0) {
				box.height +=box.y;
				box.y = 0;
			}
			if ((box.x + box.width) >= tilesX) {
				box.width = tilesX - box.x;
			}
			if ((box.y + box.height) >= tilesY) {
				box.height = tilesY - box.y;
			}
			int nTile = box.y * tilesX + box.x;
			for (int dy = 0; dy < box.height; dy++) {
				for (int dx = 0; dx < box.width;dx++) {
					dbg_layers[layer][tnImage.getNeibIndex(nTile, dx, dy)] = disp; // should never be out of bounds
				}
			}
		}
		return dbg_layers;
	}


	public void sortClustersByDisparity(
			ArrayList<PoleCluster> clusters)
	{
		Collections.sort(clusters, new Comparator<PoleCluster>() {
			@Override
			public int compare(PoleCluster lhs, PoleCluster rhs) {
				// -1 - less than, 1 - greater than, 0 - equal
				return (rhs.mean_disp > lhs.mean_disp) ? -1 : (rhs.mean_disp < lhs.mean_disp ) ? 1 : 0;
			}
		});

	}

}