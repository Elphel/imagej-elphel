/**
 ** MeasuredLayer - per-tile measured disparity/strength pairs,
 ** multiple layers can be used for 4-disparity and 2 (hor/vert) pairs
 ** separately 
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  MeasuredLayer.java is free software: you can redistribute it and/or modify
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

public class MeasuredLayers {
	private DispStrength [][] layers;
	private int tilesX;
	private int tilesY;
	private int superTileSize;
	private double [][] lapWeight = null;
	public class DispStrength {
		public double disparity = Double.NaN;
		public double strength =  0.0;
		public DispStrength (
				double disparity,
				double strength)
		{
			this.disparity = disparity;
			this.strength = strength;
			if (Double.isNaN(disparity)){
				this.strength = 0.0;
			}
		}
		public double [] getDispStrength()
		{
			double [] ds = {disparity, strength};
			return ds;
		}
		public void setDisparity(double disparity)
		{
			this.disparity = disparity;
			if (Double.isNaN(disparity)){
				this.strength = 0.0;
			}
		}
		public void setStrength(double strength)
		{
			this.strength = strength;
		}
		public double getDisparity()
		{
			return disparity;
		}
		public double getStrength()
		{
			return strength;
		}
	}
	/**
	 * Create a set of measured layers. Several layers can be used to represent alternative measurements
	 * (with different preset disparity) or for different modes (4 pair correlation, vertical, horizontal)
	 * @param num_layers number of alternative layers used (3 for quad, hor, vert)
	 * @param tilesX number of tiles horizontally in the image
	 * @param tilesY number of tiles vertically in the image
	 * @param superTileSize size of the supertile square
	 */
	public MeasuredLayers(
			int num_layers,
			int tilesX,
			int tilesY,
			int superTileSize)
	{
		layers = new DispStrength [num_layers][];
		this.tilesX = tilesX;
		this.tilesY = tilesY;
		this.superTileSize = superTileSize;
		this.lapWeight =  getLapWeights();
	}
	
	/**
	 * Get number of tiles in the image horizontally
	 * @return tilesX
	 */
	public int getTilesX()
	{
			return tilesX;
	}

	/**
	 * Get number of tiles in the image vertically
	 * @return tilesY
	 */
	public int getTilesY()
	{
			return tilesY;
	}
	/**
	 * Get number of tiles in each supertile in eac direction
	 * @return superTileSize
	 */

	public int getSuperTileSize()
	{
			return superTileSize;
	}
	
	/**
	 * Set measured layer 
	 * @param num_layer number of the layer to set
	 * @param disparity array of per-tile disparity values (linescan order),
	 *  NaN deletes the tile 
	 * @param strength array of per-tile disparity values (linescan order),
	 *  0.0 strength is OK - does not delete the tile
	 * @param selection optional tile selection (or null). If unselected, tile is deleted 
	 */
	public void setLayer (
			int       num_layer,
			double [] disparity,
			double [] strength,
			boolean [] selection) // may be null
	{	
		if (layers[num_layer] == null) {
			int llen = disparity.length;
			if (layers[num_layer] == null) {
				layers[num_layer] = new DispStrength [llen];
			}
		}
		for (int i = 0; i < disparity.length; i++){
			if ((selection != null) && !selection[i]){
				layers[num_layer][i] = null;
			} else {
				if (Double.isNaN(disparity[i])){
					layers[num_layer][i] = null;
				} else {
					if (layers[num_layer][i] != null){
						layers[num_layer][i].setDisparity(disparity[i]);
						layers[num_layer][i].setStrength(strength[i]);
					} else {
						layers[num_layer][i] = new DispStrength(disparity[i], strength[i]);
					}
				}
			}
		}
	}
	
	/**
	 * Set disparity values for selected measurement layer 
	 * @param num_layer number of the layer to set
	 * @param disparity array of per-tile disparity values (linescan order),
	 *  NaN deletes the tile. New tile set strength = 0.0 
	 */
	public void setDisparity (
			int       num_layer,
			double [] disparity)
	{	
		if (layers[num_layer] == null) {
			int llen = disparity.length;
			if (layers[num_layer] == null) {
				layers[num_layer] = new DispStrength [llen];
			}
		}
		for (int i = 0; i < disparity.length; i++){
			if (Double.isNaN(disparity[i])){
				layers[num_layer][i] = null;
			} else {
				if (layers[num_layer][i] != null){
					layers[num_layer][i].setDisparity(disparity[i]);
				} else {
					layers[num_layer][i] = new DispStrength(disparity[i], 0.0);
				}
			}
		}
	}

	/**
	 * Set strength values for selected measurement layer. Does not allocate new tiles
	 * if that element was null  
	 * @param num_layer number of the layer to set
	 * @param strength array of per-tile strength values (linescan order)
	 */
	public void setStrength (
			int       num_layer,
			double [] strength)
	{	
		if (layers[num_layer] != null) {
			for (int i = 0; i < strength.length; i++){
				if (layers[num_layer][i] != null){
					layers[num_layer][i].setStrength(strength[i]);
				}
			}
		}
	}
	
	/**
	 * Change tile selection for the layer. Can only unselect existing tiles, not add new ones
	 * @param num_layer number of the layer to set
	 * @param selection array of per-tile boolean selection values (linescan order),
	 */

	public void setSelection (
			int       num_layer,
			boolean [] selection)
	{	
		for (int i = 0; i < selection.length; i++){
			if (!selection[i]){
				layers[num_layer][i] = null;
			}
		}
	}
	
	/**
	 * Get total number of the measurement layers in this object
	 * @return number of layers (some may be uninitialized
	 */
	
	public int getNumLayers()
	{
		if (layers == null){
			return 0;
		} else {
			return layers.length;
		}
	}

	/**
	 * Get array of disparity values for the selected layer
	 * @param num_layer number of the layer to read
	 * @return array of per-tile disparity values (in linescan order),
	 *  Double.NaN for missing tiles 
	 */
	public double [] getDisparity(
			int num_layer)
	{
		double [] disparity = new double [layers[num_layer].length];
		for (int i = 0; i < disparity.length; i++){
			if (layers[num_layer][i] == null){
				disparity[i] = Double.NaN;
			} else {
				disparity[i] = layers[num_layer][i].getDisparity();
			}
		}
		return disparity;
	}
	
	/**
	 * Get array of correlation strength values for the selected layer
	 * @param num_layer number of the layer to read
	 * @return array of per-tile strength values (in linescan order),
	 * 0.0 for missing tiles 
	 */

	public double [] getStrength(
			int num_layer)
	{
		double [] strength = new double [layers[num_layer].length];
		for (int i = 0; i < strength.length; i++){
			if (layers[num_layer][i] == null){
				strength[i] = 0.0;
			} else {
				strength[i] = layers[num_layer][i].getStrength();
			}
		}
		return strength;
	}
	
	/**
	 * Get array of existing tiles for the selected layer
	 * @param num_layer number of the layer to read
	 * @return boolean (per tile) array of existing tiles
	 */
	public boolean [] getSelection(
			int num_layer)
	{
		boolean [] selection = new boolean [layers[num_layer].length];
		for (int i = 0; i < selection.length; i++){
			selection[i] = layers[num_layer][i] != null; 
		}
		return selection;
	}
	
	/**
	 * Calculate weights for overlapping supertiles to multiply correlation strengths
	 * @return square 2-d array of weights 0.0 ... 1.0, each dimension twice the
	 * supertile size
	 */
	private double [][] getLapWeights(){
		final double [][] lapWeight = new double [2 * superTileSize][2 * superTileSize];
		final double [] lapWeight1d = new double [superTileSize];
		final int superTileSize2 = 2 * superTileSize;
		for (int i = 0; i < superTileSize; i++){
			lapWeight1d[i] = 0.5*(1.0 - Math.cos((i + 0.5)* Math.PI/superTileSize));
		}
		for (int i = 0; i < superTileSize; i++){
			for (int j = 0; j < superTileSize; j++){
				lapWeight[i]                     [                     j] = lapWeight1d[i]*lapWeight1d[j]; 
				lapWeight[superTileSize2 - 1 - i][                     j] = lapWeight[i][j];
				lapWeight[i]                     [superTileSize2 - 1 - j] = lapWeight[i][j]; 
				lapWeight[superTileSize2 - 1 - i][superTileSize2 - 1 - j] = lapWeight[i][j]; 
			}
		}
		double s = 0.0;
		for (int i = 0; i < superTileSize2; i++){
			for (int j = 0; j < superTileSize2; j++){
				s+=lapWeight[i][j];
			}
		}
		System.out.println("getLapWeights: sum = "+s);
		return lapWeight;
	}
	
	public double [] getLapWeights1d(){
		int indx = 0;
		final int superTileSize2 = 2 * superTileSize;
		double [] weights = new double [superTileSize2 * superTileSize2];
		for (int i = 0; i < superTileSize2; i++){
			for (int j = 0; j < superTileSize2; j++){
				weights[indx++] = lapWeight[i][j];
			}
		}
		return weights;
	}
	
	
	/**
	 * Get selection for the specific measurement layer and supertile X,Y coordinates
	 * in the image. Combined with input selection
	 * @param num_layer number of the measurement layer to process
	 * @param stX supertile horizontal position in the image
	 * @param stY supertile vertical position in the image
	 * @param sel_in optional selection for this supertile (linescan, 4 * supetile size)
	 * @param null_if_none return null if there are no selected tiles in teh result selection
	 * @return boolean array [4*superTileSize] of the selected tiles on the specified
	 * measurement layer
	 */
	public boolean [] getSupertileSelection(
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			double     strength_floor,
			boolean null_if_none)
	{
		if ((layers[num_layer] == null) && null_if_none){
			return null;
		}
		int st2 = 2 * superTileSize;
		int st_half = superTileSize/2;
		boolean [] selection = new boolean [st2 * st2];
		int num_selected = 0;
		if (layers[num_layer] != null) {
			for (int dy = 0; dy < st2; dy ++){
				int y = superTileSize * stY -st_half + dy;
				if ((y >= 0) && (y < tilesY)) {
					for (int dx = 0; dx < st2; dx ++){
						int x = superTileSize * stX -st_half + dx;
						if ((x >= 0) && (x < tilesX)) {
							int indx = y * tilesX + x;
							int indx_st = dy * st2 + dx;
							if (((sel_in == null) || sel_in[indx_st]) && 
									(layers[num_layer][indx] != null) &&
									(layers[num_layer][indx].getStrength() >= strength_floor)){
								selection[indx_st] = true;
								num_selected ++;
							}
						}
					}
				}
			}
		}
		if (null_if_none && (num_selected == 0)) return null;
		return selection;
	}
	

	/**
	 * Get selection for the specific measurement layer and supertile X,Y coordinates
	 * in the image and combine with disparity far/near limits. Alse combined with
	 * input selection
	 * @param num_layer number of the measurement layer to process
	 * @param stX supertile horizontal position in the image
	 * @param stY supertile vertical position in the image
	 * @param sel_in optional selection for this supertile (linescan, 4 * supetile size)
	 * @param disp_far lowest acceptable disparity value. Double.NaN - do not check.
	 * @param disp_near highest acceptable disparity value. Double.NaN - do not check.
	 * @param null_if_none return null if there are no selected tiles in the result selection
	 * @return boolean array [4*superTileSize] of the selected tiles on the specified
	 * measurement layer
	 */
	public boolean [] getSupertileSelection(
			int        num_layer,
			int        stX,
			int        stY,
			boolean [] sel_in,
			double     disp_far,
			double     disp_near,
			double     strength_floor,
			boolean    null_if_none)
	{
		if ((layers[num_layer] == null) && null_if_none){
			return null;
		}
		int st2 = 2 * superTileSize;
		int st_half = superTileSize/2;
		boolean [] selection = new boolean [st2 * st2];
		int num_selected = 0;
		if (layers[num_layer] != null) {
			for (int dy = 0; dy < st2; dy ++){
				int y = superTileSize * stY -st_half + dy;
				if ((y >= 0) && (y < tilesY)) {
					for (int dx = 0; dx < st2; dx ++){
						int x = superTileSize * stX -st_half + dx;
						if ((x >= 0) && (x < tilesX)) {
							int indx = y * tilesX + x;
							int indx_st = dy * st2 + dx;
							if (((sel_in == null) || sel_in[indx_st]) && (layers[num_layer][indx] != null)){
								if (    (Double.isNaN(disp_far)  || (layers[num_layer][indx].getDisparity() >= disp_far)) &&
										(Double.isNaN(disp_near) || (layers[num_layer][indx].getDisparity() <= disp_near)) &&
										(layers[num_layer][indx].getStrength() >= strength_floor)){
									selection[indx_st] = true;
									num_selected ++;
								}
							}
						}
					}
				}
			}
		}
		if (null_if_none && (num_selected == 0)) return null;
		return selection;
	}
	
	/**
	 * Get selection when disparity/strength is already calculated. Useful when
	 * disparity/strength are not just measured values, by sampled over an area
	 * @param disparityStrength array of {disparity, strength} where each is
	 *  a linescan data of the (2 * superTileSize) * (2 * superTileSize)   
	 * @param sel_in optional input selection array or null - use all
	 * @param null_if_none return null if selection has no tiles enabled
	 * @return boolean array of per-tile enabled/disabled values in linescan
	 *  order, (2 * superTileSize) * (2 * superTileSize)
	 */

	public boolean [] getSupertileSelection(
			double [][] disparityStrength,
			boolean [] sel_in,
			boolean null_if_none)
	{
		int st2 = 2 * superTileSize;
		boolean [] selection = new boolean [st2 * st2];
		int num_selected = 0;
		for (int i = 0; i < disparityStrength[1].length; i++){
			if (	(disparityStrength[1][i] > 0.0) &&
					((sel_in == null) || (sel_in[i]))){
				selection[i] = true;
				num_selected ++;
			}
		}
		
		if (null_if_none && (num_selected == 0)) return null;
		return selection;
	}
	
	/**
	 * Get selection when disparity/strength is already calculated. Useful when
	 * disparity/strength are not just measured values, by sampled over an area
	 * Includes low/high limits for disparity values
	 * @param disparityStrength array of {disparity, strength} where each is
	 *  a linescan data of the (2 * superTileSize) * (2 * superTileSize)   
	 * @param sel_in optional input selection array or null - use all
	 * @param disp_far low limit for disparity, Double.NaN - do not check
	 * @param disp_near high limit for disparity, Double.NaN - do not check
	 * @param null_if_none return null if selection has no tiles enabled
	 * @return boolean array of per-tile enabled/disabled values in linescan
	 *  order, (2 * superTileSize) * (2 * superTileSize)
	 * @return
	 */
	public boolean [] getSupertileSelection(
			double [][] disparityStrength,
			boolean [] sel_in,
			double     disp_far,
			double     disp_near,
			boolean null_if_none)
	{
		int st2 = 2 * superTileSize;
		boolean [] selection = new boolean [st2 * st2];
		int num_selected = 0;
		for (int i = 0; i < disparityStrength[1].length; i++){
			if (	(disparityStrength[1][i] > 0.0) &&
					((sel_in == null) || (sel_in[i]))){
				if (    (Double.isNaN(disp_far)  || (disparityStrength[0][i] >= disp_far)) &&
						(Double.isNaN(disp_near) || (disparityStrength[0][i] <= disp_near))) {
					
				}
				selection[i] = true;
				num_selected ++;
			}
		}
		
		if (null_if_none && (num_selected == 0)) return null;
		return selection;
	}
	
	
	
	
	/**
	 * Get number of "true" elements in a boolean array. Null is OK, it results in 0
	 * @param selected boolean array to count set elements in
	 * @return number of selected elements
	 */
	public static int getNumSelected(
			boolean [] selected)
	{
		if (selected == null) return 0;
		int num_selected = 0;
		for (int i = 0; i < selected.length; i++){
			if (selected[i]) num_selected ++;
		}
		return num_selected;
	}

	public static double getSumStrength(
			double [][] disp_strength)
	{
		if (disp_strength == null) return 0.0;
		double sw = 0.0;
		for (int i = 0; i < disp_strength[1].length; i++){
			sw += disp_strength[1][i];
		}
		return sw;
	}
	public static double getSumStrength(
			double [][] disp_strength,
			boolean [] selected)
	{
		if (disp_strength == null) return 0.0;
		double sw = 0.0;
		for (int i = 0; i < disp_strength[1].length; i++){
			if ((selected == null) || selected[i]) {
				sw += disp_strength[1][i];
			}
		}
		return sw;
	}
	
	/**
	 * Get disparity and correlation strength for the specific measurement layer and
	 * supertile X,Y coordinates in the image. Combined with input selection
	 * @param num_layer number of the measurement layer to process
	 * @param stX supertile horizontal position in the image
	 * @param stY supertile vertical position in the image
	 * @param sel_in optional selection for this supertile (linescan, 4 * supetile size)
	 * @param strength_floor subtract from the correlation strength, limit by 0
	 * @param strength_pow Non-linear treatment of correlation strength. Raise data to the strength_pow power
	 * @param null_if_none return null if there are no selected tiles in the result selection
	 * @return double [2][4*superTileSize] {disparity[4*superTileSize], strength [4*superTileSize]}
	 */

	public double[][] getDisparityStrength (
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			double strength_floor,
			double strength_pow,
			boolean null_if_none)
	{
		if ((layers[num_layer] == null) && null_if_none){
			return null;
		}
		int st2 = 2 * superTileSize;
		int st_half = superTileSize/2;
		double [][] ds = new double [2][st2*st2];
		int num_selected = 0;
		if (layers[num_layer] != null) {
			for (int dy = 0; dy < st2; dy ++){
				int y = superTileSize * stY -st_half + dy;
				if ((y >= 0) && (y < tilesY)) {
					for (int dx = 0; dx < st2; dx ++){
						int x = superTileSize * stX -st_half + dx;
						if ((x >= 0) && (x < tilesX)) {
							int indx = y * tilesX + x;
							int indx_st = dy * st2 + dx;
							if (((sel_in == null) || sel_in[indx_st]) && (layers[num_layer][indx] != null)){
								ds[0][indx_st] = layers[num_layer][indx].getDisparity();
								double w = layers[num_layer][indx].getStrength() - strength_floor;
								if (w > 0) {
									if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
									w *= lapWeight[dy][dx];
									ds[0][indx_st] = layers[num_layer][indx].getDisparity();
									ds[1][indx_st] = w;
									num_selected ++;
								}
							}
						}
					}
				}
			}
		}
		if (null_if_none && (num_selected == 0)) return null;
		return ds;
	}

	/**
	 * Get disparity and correlation strength for the specific measurement layer.
	 * Combined with input selection
	 * @param num_layer number of the measurement layer to process
	 * @param strength_floor subtract from the correlation strength, limit by 0
	 * @param strength_pow Non-linear treatment of correlation strength. Raise data to the strength_pow power
	 * @return double {disparity[tilesX * tilesY], strength[tilesX * tilesY]}
	 */

	public double[][] getDisparityStrength (
			int num_layer,
			double strength_floor,
			double strength_pow)
	{
		if (layers[num_layer] == null){
			return null;
		}
		double [][] ds = {getDisparity(num_layer),getStrength(num_layer)};
		for (int i = 0; i < ds[1].length; i++){
			double w = ds[1][i] - strength_floor;
			if (w > 0) {
				if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
				ds[1][i] = w;
			} else {
				ds[1][i] = 0.0;
			}
		}

		return ds;
	}
	
	
	
	
	/**
	 * Get double-size  (for overlapping) array of disparities and strengths for the supertile
	 * Using best (producing lowest disparity variance) subset of neighbor tiles
	 * @param num_layer number of measurement layer (currently 0 - composite, 1 - quad, 2 - horizontal
	 * and 3 - vertical pairs correlation
	 * @param stX supertile horizontal index
	 * @param stY supertile vertical index
	 * @param sel_in input selection of the output data samples (or null)
	 * @param strength_floor subtract from the correlation strength (limit by 0) before using as a
	 *  sample weight
	 * @param strength_pow raise correlation strength (after subtracting strength_floor) to this power
	 *  to use as a sample weight
	 * @param smplSide size of the square sample side 
	 * @param smplNum number of averaged samples (should be <= smplSide * smplSide and > 1) 
	 * @param smplRms maximal square root of variance (in disparity pixels) to accept the result
	 * @param null_if_none return null if there are no usable tiles in the result
	 * @return a pair of arrays (disparity and strengths) in line-scan order each
	 */
	public double[][] getDisparityStrength (
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			double     strength_floor,
			double     strength_pow,
			int        smplSide, //        = 2;      // Sample size (side of a square)
			int        smplNum, //         = 3;      // Number after removing worst (should be >1)
			double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
			boolean null_if_none)
	{
		if ((layers[num_layer] == null) && null_if_none){
			return null;
		}
		int st2 = 2 * superTileSize;
		int st_half = superTileSize/2;
		double [][] ds = new double [2][st2*st2];
		int num_selected = 0;
		int smpl_center = smplSide /2;
		int st2e = st2 + smplSide;
		int smplLen = smplSide*smplSide;
		double [] disp =     new double [st2e * st2e]; 
		double [] weight = new double [st2e * st2e]; 
		int st_halfe = st_half + smpl_center;
		double smlVar = smplRms * smplRms; // maximal variance (weighted average of the squared difference from the mean)
		
		if (layers[num_layer] != null) {
			for (int dy = 0; dy < st2e; dy ++){
				int y = superTileSize * stY -st_halfe + dy;
				if ((y >= 0) && (y < tilesY)) {
					for (int dx = 0; dx < st2e; dx ++){
						int x = superTileSize * stX -st_halfe + dx;
						if ((x >= 0) && (x < tilesX)) {
							int indx = y * tilesX + x;
							int indx_ste = dy * st2e + dx;
							if (layers[num_layer][indx] != null){ // apply sel_in later
								disp[indx_ste] = layers[num_layer][indx].getDisparity();
								double w = layers[num_layer][indx].getStrength() - strength_floor;
								if (w > 0) {
									if (strength_pow != 1.0) w = Math.pow(w, strength_pow);
//									w *= lapWeight[dy][dx];
									disp[indx_ste] = layers[num_layer][indx].getDisparity();
									weight[indx_ste] = w;
									num_selected ++;
								}
							}
						}
					}
				}
			}
		}
		if (null_if_none && (num_selected == 0)) return null;
		// now work with disp, strength [st2e*st2de] and filter results to ds[2][st2*st2], applying sel_in
		num_selected = 0;
		for (int dy = 0; dy < st2; dy ++){
			for (int dx = 0; dx < st2; dx ++){
				int indx = dy * st2 + dx;
				if (((sel_in == null) || sel_in[indx])){
					int num_in_sample = 0;
					boolean [] smpl_sel = new boolean [smplLen];
					double [] smpl_d =  new double [smplLen];
					double [] smpl_w =  new double [smplLen];
					for (int sy = 0; sy < smplSide; sy++){
						int y = dy + sy; //  - smpl_center;
						for (int sx = 0; sx < smplSide; sx++){
							int x = dx + sx; // - smpl_center;
							int indxe = y * st2e + x;
							if (weight[indxe] > 0.0){
								int indxs = sy * smplSide + sx;
								smpl_sel[indxs] = true;
								smpl_d[indxs] = disp[indxe]; 
								smpl_w[indxs] = weight[indxe]; 
								num_in_sample ++;
							}
						}
					}
					if (num_in_sample >= smplNum){ // try, remove worst
						// calculate 
						double sd=0.0, sd2 = 0.0, sw = 0.0;
						for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
							double dw = smpl_d[i] * smpl_w[i];
							sd += dw;
							sd2 += dw * smpl_d[i];
							sw +=       smpl_w[i];
						}
						// remove worst, update sd2, sd and sw
						while ((num_in_sample > smplNum) && (sw > 0)){ // try, remove worst
							double d_mean = sd/sw;
							int iworst = -1;
							double dworst2 = 0.0;
							for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
								double d2 = (smpl_d[i] - d_mean);
								d2 *=d2;
								if (d2 > dworst2) {
									iworst = i;
									dworst2 = d2;
								}
							}
							if (iworst < 0){
								System.out.println("**** this is a BUG in getDisparityStrength() ****");
								break;
							}
							// remove worst sample
							smpl_sel[iworst] = false;
							double dw = smpl_d[iworst] * smpl_w[iworst];
							sd -= dw; 
							sd2 -= dw * smpl_d[iworst];
							sw -=       smpl_w[iworst];
							num_in_sample --;
						}
						// calculate variance of the remaining set
						if (sw > 0.0) {
							sd /= sw;
							sd2 /= sw;
							double var = sd2 - sd * sd;
							if (var < smlVar) { // good, save in the result array
								ds[0][indx] = sd;
								ds[1][indx] = sw * lapWeight[dy][dx] /num_in_sample; // average weights, multiply by window
							}
						} else {
							num_in_sample = 0;
							System.out.println("**** this is a BUG in getDisparityStrength(), shoud not happen ? ****");
						}
					}
				}				
			}
		}
		return ds;
	}
	

	public void growSelection(
			int        grow,           // grow tile selection by 1 over non-background tiles 1: 4 directions, 2 - 8 directions, 3 - 8 by 1, 4 by 1 more
			boolean [] tiles,
			boolean [] prohibit)
	{
		boolean [] src_tiles = tiles.clone(); // just in case
		// grow
		boolean hor = true;
		for (; grow > 0; grow--){
			boolean single = (grow ==1) && hor;
			src_tiles = tiles.clone();
			if (hor){
				for (int tileY = 0; tileY < tilesY; tileY++){
					for (int tileX = 0; tileX < (tilesX - 1); tileX++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx + 1])) {
							tiles[tindx + 1] |= src_tiles[tindx];  
						}

					}
					for (int tileX = 1; tileX < tilesX; tileX++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx - 1])) {
							tiles[tindx - 1] |= src_tiles[tindx];
						}
					}
				}
			}
			if (!hor || single){ // do vertically, but from previous state
				for (int tileX = 0; tileX < tilesX; tileX++){
					for (int tileY = 0; tileY < (tilesY - 1); tileY++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx + tilesX])) {
							tiles[tindx + tilesX] |= src_tiles[tindx];
						}

					}
					for (int tileY = 1; tileY < tilesY; tileY++){
						int tindx = tileY * tilesX + tileX;
						if ((prohibit == null) || (!prohibit[tindx] && !prohibit[tindx - tilesX])) {
							tiles[tindx - tilesX] |= src_tiles[tindx];
						}
					}
				}
			}
			hor = !hor;
		}
	}
	
}
