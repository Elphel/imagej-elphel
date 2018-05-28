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
	 * Get window function for tile samples (currently just 3x3) in a line-scan order
	 * @param smplSide square sample side
	 * @return [smplSide * smplSide]array of weights, 1.0 in the center
	 */

	static double [] getSampleWindow(int smplSide, boolean all1){
		double [] weights = new double [smplSide * smplSide];
		for (int sy = 0; sy < smplSide; sy++){
			for (int sx = 0; sx < smplSide; sx++){
				weights[sy * smplSide + sx] = all1? 1.0 : Math.sin(Math.PI*(sy+0.5)/smplSide) * Math.sin(Math.PI*(sx+0.5)/smplSide);
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

	public double[][] getDisparityStrengthML (
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			MeasuredLayersFilterParameters mlfp,
			boolean null_if_none){
		return getDisparityStrengthML(
				num_layer,
				stX,
				stY,
				sel_in,
				mlfp.strength_floor,
				mlfp.strength_pow,
				null_if_none);
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

	public double[][] getDisparityStrengthML (
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

	public double[][] getDisparityStrengthML (
			int num_layer,
			MeasuredLayersFilterParameters mlfp){

		return getDisparityStrengthML(
				num_layer,
				mlfp.strength_floor,
				mlfp.strength_pow);
	}


	/**
	 * Get disparity and correlation strength for the specific measurement layer.
	 * Combined with input selection
	 * @param num_layer number of the measurement layer to process
	 * @param strength_floor subtract from the correlation strength, limit by 0
	 * @param strength_pow Non-linear treatment of correlation strength. Raise data to the strength_pow power
	 * @return double {disparity[tilesX * tilesY], strength[tilesX * tilesY]}
	 */

	public double[][] getDisparityStrengthML (
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
	 * @param smplWnd multiply samples weights by a window function
	 * @param debugLevel debug level
	 * @return a pair of arrays (disparity and strengths) in line-scan order each
	 */
	public double[][] getDisparityStrength_old (
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			MeasuredLayersFilterParameters mlfp,
//			double     strength_floor,
//			double     strength_pow,
//			int        smplSide, //        = 2;      // Sample size (side of a square)
//			int        smplNum, //         = 3;      // Number after removing worst (should be >1)
//			double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//			boolean    smplWnd, //
			boolean    null_if_none,
			int        debugLevel)
	{
		if ((layers[num_layer] == null) && null_if_none){
			return null;
		}
		int st2 = 2 * superTileSize;
		int st_half = superTileSize/2;
		double [][] ds = new double [2][st2*st2];
		final int dbg_tile = -1; // = ((stX == 22) && (stY == 19)) ? (5 + 7*16) : -1;// 50397;

		int num_selected = 0;
		int smpl_center = mlfp.smplSide /2;
		int st2e = st2 + mlfp.smplSide;
		int smplLen = mlfp.smplSide*mlfp.smplSide;
		final double [] smpl_weights = getSampleWindow(mlfp.smplSide, !mlfp.smplWnd);

		double [] disp =     new double [st2e * st2e];
		double [] weight = new double [st2e * st2e];
		int st_halfe = st_half + smpl_center;
		double smplVar = mlfp.smplRms * mlfp.smplRms; // maximal variance (weighted average of the squared difference from the mean)


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
								double w = layers[num_layer][indx].getStrength() - mlfp.strength_floor;
								if (w > 0) {
									if (mlfp.strength_pow != 1.0) w = Math.pow(w, mlfp.strength_pow);
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
				if (indx == dbg_tile){
					System.out.println("getDisparityStrengthML(): stX="+stX+" stY="+stY+" dx="+dx+" dy="+dy);
				}
				if (((sel_in == null) || sel_in[indx])){
					int num_in_sample = 0;
					double sum_wnd = 0.0;
					boolean [] smpl_sel = new boolean [smplLen];
					double [] smpl_d =  new double [smplLen];
					double [] smpl_w =  new double [smplLen];
					for (int sy = 0; sy < mlfp.smplSide; sy++){
						int y = dy + sy; //  - smpl_center;
						for (int sx = 0; sx < mlfp.smplSide; sx++){
							int x = dx + sx; // - smpl_center;
							int indxe = y * st2e + x;
							if (weight[indxe] > 0.0){
								int indxs = sy * mlfp.smplSide + sx;
								smpl_sel[indxs] = true;
								smpl_d[indxs] = disp[indxe];
								smpl_w[indxs] = weight[indxe] * smpl_weights[indxs];
								sum_wnd += smpl_weights[indxs];
								num_in_sample ++;
							}
						}
					}
					if (num_in_sample >= mlfp.smplNum){ // try, remove worst
						// calculate
						double sd=0.0, sd2 = 0.0, sw = 0.0;
						for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
							double dw = smpl_d[i] * smpl_w[i];
							sd += dw;
							sd2 += dw * smpl_d[i];
							sw +=       smpl_w[i];
						}
						// remove worst, update sd2, sd and sw
						while ((num_in_sample > mlfp.smplNum) && (sw > 0)){ // try, remove worst
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
								System.out.println("**** this is a BUG1 in getDisparityStrengthML() ****");
								break;
							}
							// remove worst sample
							smpl_sel[iworst] = false;
							double dw = smpl_d[iworst] * smpl_w[iworst];
							sd -= dw;
							sd2 -= dw * smpl_d[iworst];
							sw -=       smpl_w[iworst];
							sum_wnd -= smpl_weights[iworst];
							num_in_sample --;
						}
						// calculate variance of the remaining set
						if (sw > 0.0) {
							sd /= sw;
							sd2 /= sw;
							double var = sd2 - sd * sd;
							if (var < smplVar) { // good, save in the result array
								ds[0][indx] = sd;
//								ds[1][indx] = sw * lapWeight[dy][dx] /num_in_sample; // average weights, multiply by window //** TODO: change
								ds[1][indx] = sw * lapWeight[dy][dx] /sum_wnd; // average weights, multiply by window //** TODO: change
							}
						} else {
							num_in_sample = 0;
							System.out.println("**** this is a BUG in getDisparityStrengthML(), shoud not happen ? ****");
						}
					}
				}
			}
		}
		return ds;
	}

	/**
	 * Verify that selected points are not all on the same line
	 * @param sel 2-d sample selection in linescan order
	 * @param side square samples side
	 * @return true if there are enough samples for plane extraction, false otherwise
	 */
	public boolean notColinear (
			boolean [] sel,
			int side)
	{
		int indx0, indx1;
		for (indx0 = 0; indx0 < sel.length; indx0++){
			if (sel[indx0]) break;
		}
		for (indx1 = indx0+1; indx1 < sel.length; indx1++){
			if (sel[indx1]) break;
		}
		if (indx1 >= sel.length) return false; // too few points;
		int sx0 = indx0 % side;
		int sy0 = indx0 / side;
		int sx1 = indx1 % side;
		int sy1 = indx1 / side;
		for (int indx = indx1 +1; indx < sel.length; indx++){
			int sx = indx % side;
			int sy = indx / side;
			if ((sx - sx0) * (sy - sy1) != (sx - sx1) * (sy - sy0)){
				return true;
			}
		}
		return false;
	}

	/**
	 * Verify that selected points are not all on the same line, even if the specified one is removed
	 * @param indx index of the point to be removed
	 * @param sel 2-d sample selection in linescan order
	 * @param side square samples side
	 * @return true if there are enough samples for plane extraction, false otherwise
	 */

	public boolean notColinearWithout (
			int        indx,
			boolean [] sel,
			int side)
	{
		if (!sel[indx]){
			throw new IllegalArgumentException ("notCoplanarWithout(): specified is the non existing index");
		}
		sel[indx] = false;
		boolean rslt = notColinear ( sel, side);
		sel[indx] = true;
		return rslt;
	}

	// testing - redirecting all existing requests to this one with floating planes
	public double[][] getDisparityStrengthMLTilted (
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			MeasuredLayersFilterParameters mlfp,     // filter parameters
//			double     strength_floor,
//			double     strength_pow,
//			int        smplSide, //        = 2;      // Sample size (side of a square)
//			int        smplNum, //         = 3;      // Number after removing worst (should be >1)
//			double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//			boolean    smplWnd, //

//  			double     max_abs_tilt,  //  2.0;   // pix per tile
//			double     max_rel_tilt,  //  0.2;   // (pix / disparity) per tile
//			double     damp_tilt,     //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//			double     min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//			double     transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//			int        far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//			double     far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects

			boolean    null_if_none,
			int        debugLevel)
	{
		boolean use_new = true; // false;

		if (use_new) {
		return getDisparityStrengthMLTilted (
				num_layer,      // int num_layer,
				stX,            // int stX,
				stY,            // int stY,
				sel_in,         // boolean [] sel_in,
				null,           // double []  tiltXY, // null - free with limit on both absolute (2.0?) and relative (0.2) values
				mlfp,     // MeasuredLayersFilterParameters mlfp,     // filter parameters
//				strength_floor, // double     strength_floor,
//				strength_pow,   // double     strength_pow,
//				smplSide,       // int        smplSide, //        = 2;      // Sample size (side of a square)
//				smplNum,        // int        smplNum, //         = 3;      // Number after removing worst (should be >1)
//				smplRms,        // double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//				smplWnd,        // boolean    smplWnd, //
//				max_abs_tilt,   // double     max_abs_tilt, //  = 2.0; // pix per tile
//				max_rel_tilt,   // double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
//				damp_tilt,      //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//				min_tilt_disp,  //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//				transition,     //  1.0;   // Mode transition range (between tilted and maximal disparity)
//				far_mode,       //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//				far_power,      //  1.0;   // Raise disparity to this power before averaging for far objects
				null_if_none,   // boolean    null_if_none,
				debugLevel);    // int        debugLevel)
		} else {
			return getDisparityStrength_old (
					num_layer,      // int num_layer,
					stX,            // int stX,
					stY,            // int stY,
					sel_in,         // boolean [] sel_in,
//					null,           // double []  tiltXY, // null - free with limit on both absolute (2.0?) and relative (0.2) values
					mlfp,     // MeasuredLayersFilterParameters mlfp,     // filter parameters
//					strength_floor, // double     strength_floor,
//					strength_pow,   // double     strength_pow,
//					smplSide,       // int        smplSide, //        = 2;      // Sample size (side of a square)
//					smplNum,        // int        smplNum, //         = 3;      // Number after removing worst (should be >1)
//					smplRms,        // double     smplRms, //         = 0.1;    // Maximal RMS of the remaining tiles in a sample
//					smplWnd,        // boolean    smplWnd, //
//					2.0,            // double     max_abs_tilt, //  = 2.0; // pix per tile
//					0.2,            // double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
					null_if_none,   // boolean    null_if_none,
					debugLevel);    // int        debugLevel)

		}
	}

// Try two modes of operation: Far tiles (disparity <~4?) use power of disparity, no tilts. Near tiles - current.
// For far tiles power can use closest (lowest disparity) - either power or just census or smth.
// Detect strong background?


	public double[][] getDisparityStrengthMLTilted (
			int num_layer,
			int stX,
			int stY,
			boolean [] sel_in,
			double [][] atiltXY, // null - free with limit on both absolute (2.0?) and relative (0.2) values
			                     // if it has just one element - apply to all tiles, otherwise it is supposed
			                     // to be per-tile of a supertile array
			MeasuredLayersFilterParameters mlfp,
//			double     strength_floor,
//			double     strength_pow,
//			int        smplSide, //        = 3;      // Sample size (side of a square)
//			int        smplNum, //         = 5;      // Number after removing worst (should be >1)
//			double     smplRms, //         = 0.3;    // Maximal RMS of the remaining tiles in a sample
//			boolean    smplWnd, //
//			double     max_abs_tilt,     //  2.0;   // pix per tile
//			double     max_rel_tilt,     //  0.2;   // (pix / disparity) per tile
//			double     damp_tilt,        //  0.001; // Damp tilt to handle insufficient  (co-linear)data
//			double     min_tilt_disp, //  4.0;   // Disparity switch between filtering modes - near objects use tilts, far - use max disparity
//			double     transition,    //  1.0;   // Mode transition range (between tilted and maximal disparity)
//			int        far_mode,      //  1;     // Far objects filtering mode (0 - off, 1 - power of disparity)
//			double     far_power,     //  3.0;   // Raise disparity to this power before averaging for far objects
			boolean    null_if_none,
			int        debugLevel)
	{

		if ((layers[num_layer] == null) && null_if_none){
			return null;
		}
		final double [] damping = {mlfp.damp_tilt, mlfp.damp_tilt, 0.0}; // 0.0 will be applied to average value, tilt_cost - to both tilts

		int st2 = 2 * superTileSize;
		int st_half = superTileSize/2;
		double [][] ds = new double [2][st2*st2];
		final int dbg_tile =-1; //  ((stX == 22) && (stY == 19)) ? (5 + 7*16) : -1;// 50397;

		int num_selected = 0;
		int smpl_center = mlfp.smplSide /2;
		double smpl_dcenter = (mlfp.smplSide -1.0) /2;
		int st2e = st2 + mlfp.smplSide;
		int smplLen = mlfp.smplSide*mlfp.smplSide;
		final double [] smpl_weights = getSampleWindow(mlfp.smplSide, !mlfp.smplWnd);

		double [] disp =     new double [st2e * st2e];
		double [] weight = new double [st2e * st2e];
		int st_halfe = st_half + smpl_center;
		double smplVar = mlfp.smplRms * mlfp.smplRms; // maximal variance (weighted average of the squared difference from the mean)
		PolynomialApproximation pa = new PolynomialApproximation();
		double thresholdLin = 1.0E-20;  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
		double thresholdQuad = 1.0E-30; // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)

		double  disp_pwr_far =  mlfp.min_tilt_disp - mlfp.transition/2.0;
		double  disp_pwr_near = mlfp.min_tilt_disp + mlfp.transition/2.0;
		double  disp_pwr_near2 = disp_pwr_near + mlfp.smplRms;
		double [] disp_pwr = null;
		boolean far_mode_en = (mlfp.far_mode == 1) || (mlfp.far_mode == 2);
		boolean remove_far_only = (mlfp.far_mode == 2);
		if (far_mode_en) {
			disp_pwr = new double [st2e * st2e];
		}
		int num_far = 0;
		boolean [] smpl_far_sel = null;
		double []  smpl_pwr_d =   null;
		double []  smpl_far_d =   null;
		double []  smpl_far_w =   null;
		int num_in_sample_far = 0;

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
								double w = layers[num_layer][indx].getStrength() - mlfp.strength_floor;
								if (w > 0) {
									if (mlfp.strength_pow != 1.0) w = Math.pow(w, mlfp.strength_pow);
//									w *= lapWeight[dy][dx];
									disp[indx_ste] = layers[num_layer][indx].getDisparity(); // same?
									weight[indx_ste] = w;
									num_selected ++;
									if (far_mode_en && (disp[indx_ste] < disp_pwr_near2)  && (disp[indx_ste] > 0)) {
										disp_pwr[indx_ste] = Math.pow(disp[indx_ste], mlfp.far_power);
										num_far++;
									}
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
				if (indx == dbg_tile){
					System.out.println("getDisparityStrengthML(): stX="+stX+" stY="+stY+" dx="+dx+" dy="+dy);
				}
				if (((sel_in == null) || sel_in[indx])){
					// check if it need filtering
					int indxe_center = (dy + smpl_center)*st2e + (dx + smpl_center);
					if (weight[indxe_center] >= mlfp.strength_sure) { // tile does not need filtering
						ds[0][indx] = disp[indxe_center];
						ds[1][indx] = weight[indxe_center];
						if (Double.isNaN(ds[0][indx])){
							System.out.println("**** this is a BUG6 in getDisparityStrengthML() ****"); // all smpl_d are NaNs
//							break;
						}
					} else {
						int num_in_sample = 0;
						double sum_wnd = 0.0;
						boolean [] smpl_sel = new boolean [smplLen];
						double [] smpl_d =  new double [smplLen];
						double [] smpl_p =  new double [smplLen];
						double [] smpl_w =  new double [smplLen];
						if (far_mode_en ) {
							smpl_far_sel =  new boolean [smplLen];
							smpl_pwr_d =  new double [smplLen];
							smpl_far_w =  new double [smplLen];
							num_in_sample_far = 0;
						}
						for (int sy = 0; sy < mlfp.smplSide; sy++){
							int y = dy + sy; //  - smpl_center;
							for (int sx = 0; sx < mlfp.smplSide; sx++){
								int x = dx + sx; // - smpl_center;
								int indxe = y * st2e + x;
								if (weight[indxe] > 0.0){
									int indxs = sy * mlfp.smplSide + sx;
									smpl_sel[indxs] = true;
									smpl_d[indxs] = disp[indxe];
									if (Double.isNaN(smpl_d[indxs])){
										System.out.println("**** this is a BUG5 in getDisparityStrengthML() ****"); // all smpl_d are NaNs
//										break;
									}
									smpl_w[indxs] = weight[indxe] * smpl_weights[indxs];
									sum_wnd += smpl_weights[indxs];
									num_in_sample ++;
									if (far_mode_en && (disp_pwr[indxe]>0.0)) {
										smpl_far_sel[indxs] = true;
										smpl_pwr_d[indxs] =   disp_pwr[indxe];
										smpl_far_w[indxs] = smpl_w[indxs];
										num_in_sample_far++;
									}
								}
							}
						}
						if (far_mode_en ) {
							smpl_far_d = smpl_d.clone();
						}
						if (num_in_sample >= mlfp.smplNum){ // try, remove worst
							sample_loop:
							{
							boolean en_tilt = (atiltXY == null);
							if (en_tilt) { // make sure there are enough samples and not all of them are on the same line
								if ((mlfp.damp_tilt == 0.0) && !notColinear(smpl_sel, mlfp.smplSide)){
									en_tilt = false;
								}
							}
							if (en_tilt) { // enable floating tilt
								double sd2 = 0.0, d_center = 0.0,  sw = 0.0;

								// TODO: making simple - recalculate after removing. Can be done more efficient.
								while (num_in_sample >= mlfp.smplNum) { // try, remove worst
									sd2 = 0.0;
									d_center = 0.0;
									sw = 0.0;

									double [][][] mdata = new double [num_in_sample][3][];
									int mindx = 0;
									for (int sy = 0; sy < mlfp.smplSide; sy++){
										for (int sx = 0; sx < mlfp.smplSide; sx++){
											int indxs = sy * mlfp.smplSide + sx;
											if (smpl_sel[indxs]) {
												mdata[mindx][0] = new double [2];
												mdata[mindx][0][0] =  sx - smpl_dcenter;
												mdata[mindx][0][1] =  sy - smpl_dcenter;
												mdata[mindx][1] = new double [1];
												mdata[mindx][1][0] = smpl_d[indxs];
												mdata[mindx][2] = new double [1];
												mdata[mindx][2][0] =  smpl_w[indxs];
												mindx ++;
											}
										}
									}
									double[][] approx2d = pa.quadraticApproximation(
											mdata,
											true,          // boolean forceLinear,  // use linear approximation
											damping,       // double [] damping,
											thresholdLin,  // threshold ratio of matrix determinant to norm for linear approximation (det too low - fail)
											thresholdQuad, // threshold ratio of matrix determinant to norm for quadratic approximation (det too low - fail)
											debugLevel);
									if (approx2d == null){
										if (debugLevel > -1){
											System.out.println("getDisparityStrengthML(): can not find linear approximation");
										}
										break sample_loop;
									}
									// limit tilt to be within range
									//											double     max_abs_tilt, //  = 2.0; // pix per tile
									//											double     max_rel_tilt, //  = 0.2; // (pix / disparity) per tile
									double max_tilt = Math.min(mlfp.max_abs_tilt, mlfp.max_rel_tilt * approx2d[0][2]);
									boolean overlimit = (Math.abs(approx2d[0][0]) > max_tilt) || (Math.abs(approx2d[0][1]) > max_tilt);
									if (overlimit) {
										approx2d[0][0] = Math.min(approx2d[0][0],  max_tilt);
										approx2d[0][1] = Math.min(approx2d[0][1],  max_tilt);
										approx2d[0][0] = Math.max(approx2d[0][0], -max_tilt);
										approx2d[0][1] = Math.max(approx2d[0][1], -max_tilt);
									}
									// subtract tilt from disparity
									for (int sy = 0; sy < mlfp.smplSide; sy++){
										for (int sx = 0; sx < mlfp.smplSide; sx++){
											int indxs = sy * mlfp.smplSide + sx;
											if (smpl_sel[indxs]) {
												smpl_p[indxs] = approx2d[0][0] * (sx - smpl_dcenter) + approx2d[0][1] * (sy - smpl_dcenter) + approx2d[0][2];
											}
										}
									}

									if (overlimit){ // re-calculate disparity average (in the center)
										double sd=0.0;
										for (int indxs = 0; indxs < smplLen;indxs++) if (smpl_sel[indxs]) {
											double d = smpl_d[indxs] - smpl_p[indxs];
											double dw = d * smpl_w[indxs];
											sd += dw;
											sw += smpl_w[indxs];
										}
										sd /= sw;
										for (int indxs = 0; indxs < smplLen;indxs++) if (smpl_sel[indxs]) {
											smpl_p[indxs] += sd;
										}
										approx2d[0][2] += sd;
									}
									d_center = approx2d[0][2];
									sw = 0.0;
									for (int indxs = 0; indxs < smplLen;indxs++) if (smpl_sel[indxs]) {
										double d = smpl_d[indxs] - smpl_p[indxs];
										double dw = d * smpl_w[indxs];
										//									sd += dw;
										sd2 += dw * smpl_d[indxs];
										sw +=       smpl_w[indxs];
									}


									// remove worst - it should not make remaining set
									if (num_in_sample > mlfp.smplNum) { // remove worst if it is not the last run where only calculations are needed
										//									double d_mean = sd/sw;
										int iworst = -1;
										double dworst2 = 0.0;
										for (int indxs = 0; indxs < smplLen; indxs++) if (smpl_sel[indxs]) {
											//										double d2 = (smpl_d[i] - d_mean);
											double d2 = smpl_d[indxs] - smpl_p[indxs];
											d2 *=d2;
											if (d2 > dworst2) {
												if ((mlfp.damp_tilt !=0.0) || notColinearWithout (
														indxs, // int        indx,
														smpl_sel, // boolean [] sel,
														mlfp.smplSide)) { // int side))
													iworst = indxs;
													dworst2 = d2;
												}
											}
										}
										if (iworst < 0){
											if (debugLevel > 0) {
												System.out.println("**** this may be BUG in getDisparityStrengthML() can not find the worst sample  - all tiles fit perfectly ****");
											}
											// this can happen if some samples are the same and all the pixels fit exactly - use all of them
											break;
										}
										// remove worst sample
										smpl_sel[iworst] = false;
										//									double dw = smpl_d[iworst] * smpl_w[iworst];
										//									sd -= dw;
										//									sd2 -= dw * smpl_d[iworst];
										//									sw -=       smpl_w[iworst];
										sum_wnd -= smpl_weights[iworst];
										num_in_sample --;
									} else {
										break;
									}

								} // removing worst tiles, all done,
								// calculate variance of the remaining set
								if (sw > 0.0) {
									//								sd /= sw;
									//								sd2 /= sw;
									double var = sd2/sw; //   - sd * sd;
									if (var < smplVar) { // good, save in the result array
										ds[0][indx] = d_center;
										//								ds[1][indx] = sw * lapWeight[dy][dx] /num_in_sample; // average weights, multiply by window //** TODO: change
										ds[1][indx] = sw * lapWeight[dy][dx] /sum_wnd; // average weights, multiply by window //** TODO: change
									}
								}

							} else { // fixed tilt
								// tilt around center
								double [] tiltXY= {0.0,0.0};
								if (atiltXY != null){ // if null ( free but failed co-linear test - keep both tilts == 0.0
									if (atiltXY.length == 1){
										tiltXY = atiltXY[0]; // common for all tiles (such as constant disparity)
									} else if (atiltXY[indx] != null){
										tiltXY = atiltXY[indx];
										if (Double.isNaN(tiltXY[0]) || Double.isNaN(tiltXY[1])){
											System.out.println("**** this is a BUG4A in getDisparityStrengthML() ****");
											System.out.println("atilt["+indx+"]= {"+tiltXY[0]+","+tiltXY[1]+"}");
											tiltXY[0] = 0.0;
											tiltXY[1] = 0.0;
										}

									}
								}

								if ((tiltXY != null) && (tiltXY[0] != 0.0) && (tiltXY[1] != 0.0)){ // {NaN,NaN}
									for (int sy = 0; sy < mlfp.smplSide; sy++){
										for (int sx = 0; sx < mlfp.smplSide; sx++){
											int indxs = sy * mlfp.smplSide + sx;
											if (smpl_w[indxs] > 0.0) {
												smpl_d[indxs] -= tiltXY[0]* (sx - smpl_dcenter) + tiltXY[1]* (sy - smpl_dcenter);
												if (Double.isNaN(smpl_d[indxs])){
													System.out.println("**** this is a BUG4 in getDisparityStrengthML() ****"); // all smpl_d are NaNs
//													break;
												}

											}
										}
									}
								}


								// calculate
								double sd=0.0, sd2 = 0.0, sw = 0.0;
								for (int i = 0; i < smplLen; i++) if (smpl_sel[i]) {
									if (Double.isNaN(smpl_d[i])){
										System.out.println("**** this is a BUG3 in getDisparityStrengthML() ****"); // all smpl_d are NaNs
//										break;
									}

									double dw = smpl_d[i] * smpl_w[i];
									sd += dw;
									sd2 += dw * smpl_d[i];
									sw +=       smpl_w[i];
								}
								// remove worst, update sd2, sd and sw
								while ((num_in_sample > mlfp.smplNum) && (sw > 0)){ // try, remove worst
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
										//remove_far_only
									}
									if (iworst < 0){
										System.out.println("**** this is a BUG2 in getDisparityStrengthML() ****"); // all smpl_d are NaNs
										break;
									}
									// remove worst sample
									smpl_sel[iworst] = false;
									double dw = smpl_d[iworst] * smpl_w[iworst];
									sd -= dw;
									sd2 -= dw * smpl_d[iworst];
									sw -=       smpl_w[iworst];
									sum_wnd -= smpl_weights[iworst];
									num_in_sample --;
								}
								// calculate variance of the remaining set
								if (sw > 0.0) {
									sd /= sw;
									sd2 /= sw;
									double var = sd2 - sd * sd;
									if (var < smplVar) { // good, save in the result array
										ds[0][indx] = sd;
										//								ds[1][indx] = sw * lapWeight[dy][dx] /num_in_sample; // average weights, multiply by window //** TODO: change
										ds[1][indx] = sw * lapWeight[dy][dx] /sum_wnd; // average weights, multiply by window //** TODO: change
									}
								} else {
									num_in_sample = 0;
									System.out.println("**** this is a BUG in getDisparityStrengthML(), shoud not happen ? ****");
								}
							}
							}
						}
						// now apply far if needed
						if (num_in_sample_far > 0) { // calculate small far disparity, then merge with normal one
							// remove extra
							// smpl_pwr_d - disparity to power
							// smpl_far_d - just disparity
							if (num_in_sample_far >= mlfp.smplNum){ // try, remove worst
								double sd=0.0, sd2 = 0.0, sw = 0.0;
								double sdp =0.0; // , sdp2 = 0.0;
								for (int i = 0; i < smplLen; i++) if (smpl_far_sel[i]) {
									double dw = smpl_far_d[i] * smpl_far_w[i];
									sd += dw;
									sd2 += dw * smpl_far_d[i];
									sw +=       smpl_far_w[i];
									double dpw = smpl_pwr_d[i] * smpl_far_w[i];
									sdp += dpw;
									//								sdp2 += dpw * smpl_pwr_d[i];
								}


								// remove worst (before power), update sd2, sd and sw
								while ((num_in_sample_far > mlfp.smplNum) && (sw > 0)){ // try, remove worst
									double d_mean = sd/sw;
									int iworst = -1;
									double dworst2 = 0.0;
									for (int i = 0; i < smplLen; i++) if (smpl_far_sel[i]) {
										double d2 = (smpl_far_d[i] - d_mean);
										d2 *=d2;
										if ((d2 > dworst2) && (!remove_far_only ||(smpl_far_d[i] < d_mean))) {
											iworst = i;
											dworst2 = d2;
										}
									}
									if (iworst < 0){
										System.out.println("**** this is a BUG3 in getDisparityStrengthML() ****");
										break;
									}
									// remove worst sample
									smpl_far_sel[iworst] = false;
									double dw = smpl_far_d[iworst] * smpl_far_w[iworst];
									sd -= dw;
									sd2 -= dw * smpl_far_d[iworst];
									sw -=       smpl_far_w[iworst];
									sum_wnd -= smpl_weights[iworst];
									num_in_sample_far --;

									// Remove from the power average too

									sdp -= smpl_pwr_d[iworst] * smpl_far_w[iworst];
								}
								// calculate variance of the remaining set
								if (sw > 0.0) {
									sd /= sw;
									sd2 /= sw;
									double var = sd2 - sd * sd;
									if (var < smplVar) { // good, save in the result array (here use linear disparity
										// here need to combine two filtering modes - "tilted" and "far"
										sdp = Math.pow(sdp/sw, 1.0/mlfp.far_power);
										if (sdp < disp_pwr_near) { // otherwise keep tilted mode result
											double far_w = sw * lapWeight[dy][dx] /sum_wnd; // average weights, multiply by window //** TODO: change
											if (sdp < disp_pwr_far) {
												ds[0][indx] = sdp;
												ds[1][indx] = far_w;
											} else { // interpolate
												double k = (sdp - disp_pwr_far) / mlfp.transition;
												ds[0][indx] = (1.0 - k) * sdp +   k * ds[0][indx];
												ds[1][indx] = (1.0 - k) * far_w + k * ds[1][indx];
											}
										}
										ds[0][indx] = sd;
										ds[1][indx] = sw * lapWeight[dy][dx] /sum_wnd; // average weights, multiply by window //** TODO: change
									}
								} else {
									num_in_sample = 0;
									System.out.println("**** this is a BUG4 in getDisparityStrengthML(), shoud not happen ? ****");
								}
							}
						}
					}
				} // end of tile
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
