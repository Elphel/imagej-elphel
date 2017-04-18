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
	
	public int getnumLayers()
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
	 * @param null_if_none return null if there are no selected tiles in teh result selection
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
	
	/**
	 * Get disparity and correlation strength for the specific measurement layer and
	 * supertile X,Y coordinates in the image. Combined with input selection
	 * @param num_layer number of the measurement layer to process
	 * @param stX supertile horizontal position in the image
	 * @param stY supertile vertical position in the image
	 * @param sel_in optional selection for this supertile (linescan, 4 * supetile size)
	 * @param null_if_none return null if there are no selected tiles in the result selection
	 * @param strength_pow Non-linear treatment of correlation strength. Raise data to the strength_pow power
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
}
