/**
 **
 ** MacroCorrelation handle "correlations of correlations" - tiles instead of pixels
 ** to find the ranges for pixel correlations
 **
 ** Copyright (C) 2017 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **  
 **  MacroCorrelation.java is free software: you can redistribute it and/or modify
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

public class MacroCorrelation {
	TileProcessor tp; // pixel tile processor
	TileProcessor mtp; // macro tile processor
	public MacroCorrelation(
			TileProcessor tp,
			double trusted_correlation){
		this.tp = tp;
		  final int pTilesX = tp.getTilesX();
		  final int pTilesY = tp.getTilesY();
		  final int tileSize = tp.getTileSize(); // 
		  final int mTilesX = (pTilesX + tileSize - 1) / tileSize;   // clt_aberrations_quad_corr truncates
		  final int mTilesY = (pTilesY + tileSize - 1) / tileSize;
		  this.mtp = new TileProcessor(
					mTilesX,   // int tilesX,
					mTilesY,   // int tilesY,
					tileSize,  // int tileSize,   
					tp.superTileSize, // int superTileSize,
					tp.getMagicScale(), // double scale,
					trusted_correlation, // double trustedCorrelation,
					tp.getThreadsMax()); // int threadsMax)
	}
	public TileProcessor  CLTMacroScan( // perform single pass according to prepared tiles operations and disparity
			final CLTPass3d                          src_scan, // results of the normal correlations (now expecting infinity)
//			final double [][][]                      other_channels, // other channels to correlate, such as average RGB (first index - subcamera, 2-nd - channel, 3-rd - pixel)
			EyesisCorrectionParameters.CLTParameters clt_parameters,
			GeometryCorrection                       geometryCorrection,			  
			final double                             macro_disparity_low,
			final double                             macro_disparity_high,
			final double                             macro_disparity_step,
			final int                                debugLevel){

		double [][][] input_data = CLTMacroSetData( // perform single pass according to prepared tiles operations and disparity
				src_scan,           // final CLTPass3d      src_scan, // results of the normal correlations (now expecting infinity)
				null);               // final double [][][]  other_channels, // other channels to correlate, such as average RGB (first index - subcamera, 2-nd - channel, 3-rd - pixel)
				
		mtp.resetCLTPasses();
		for (double mdisp = macro_disparity_low; mdisp < (macro_disparity_high - macro_disparity_step); mdisp +=macro_disparity_step){
			final CLTPass3d macro_scan = new CLTPass3d(mtp);

			CLTMacroSetMeasure( // perform single pass according to prepared tiles operations and disparity
					macro_scan,         // final CLTPass3d                          macro_scan, // new, will be filled out
					clt_parameters,     // EyesisCorrectionParameters.CLTParameters clt_parameters,
					geometryCorrection, // GeometryCorrection   geometryCorrection,			  
					mdisp,              // final double         macro_disparity,
					false,              // final boolean                            show_corr_partial,
					false,              // final boolean                            show_corr_combo,
					debugLevel);        // final int                                debugLevel)
			
			CLTMacroMeasure( // perform single pass according to prepared tiles operations and disparity
					macro_scan, // final CLTPass3d                          macro_scan, // 
					input_data, // final double [][][]                      input_data,
					clt_parameters,     // EyesisCorrectionParameters.CLTParameters clt_parameters,
					geometryCorrection, // GeometryCorrection   geometryCorrection,	
					""+mdisp,           // final String                             suffix,
					false,              // final boolean                            show_corr_partial,
					false,              // final boolean                            show_corr_combo,
					debugLevel);        // final int                                debugLevel)
			
			mtp.clt_3d_passes.add(macro_scan);
		}
		return mtp;
		
	}
	
	public double [][][] CLTMacroSetData( // perform single pass according to prepared tiles operations and disparity
			final CLTPass3d                          src_scan, // results of the normal correlations (now expecting infinity)
			final double [][][]                      other_channels // other channels to correlate, such as average RGB (first index - subcamera, 2-nd - channel, 3-rd - pixel)
			)
	{
		final boolean tmp_input = true; // use 3 correlation channels made of the same data

		//		  final int dbg_x = 295;
		//		  final int dbg_y = 160;
		final int pTilesX = tp.getTilesX();
		final int pTilesY = tp.getTilesY();
		final int tileSize = tp.getTileSize(); // 
		final int mTilesX = (pTilesX + tileSize - 1) / tileSize;   // clt_aberrations_quad_corr truncates
		final int mTilesY = (pTilesY + tileSize - 1) / tileSize;
		final int mTiles = mTilesX * mTilesY; 
		final int num_chn = tmp_input? 3: (1 + ((other_channels == null) ? 0 : (other_channels[0].length)));
		final double [][][] input_data = new double [ImageDtt.QUAD][num_chn][mTiles*tileSize*tileSize];
		//			double [][] tiles_tone = src_scan.getTileRBGA( No, we need individual subcameras
		//					4); // int num_layers); 
		for (int sub_cam =0; sub_cam < input_data.length; sub_cam++){
			for (int pty = 0; pty < pTilesY; pty++){
				for (int ptx = 0; ptx < pTilesX; ptx++){
					int pTile = ptx + pty * pTilesX;
					int mTile = ptx + pty * (mTilesX * tileSize);
					input_data[sub_cam][0][mTile]= src_scan.disparity_map[ImageDtt.IMG_DIFF0_INDEX + sub_cam][pTile];
					if (tmp_input) {
						input_data[sub_cam][1][mTile]= input_data[sub_cam][0][mTile];
						input_data[sub_cam][2][mTile]= input_data[sub_cam][0][mTile];
					} else { 
						for (int other_chn = 1; other_chn < num_chn; other_chn++){
							input_data[sub_cam][other_chn][mTile] = other_channels[sub_cam][other_chn-1][pTile];
						}
					}
				}
			}

		}
		return input_data;
	}
	
	
	
	public CLTPass3d CLTMacroSetMeasure( // perform single pass according to prepared tiles operations and disparity
			final CLTPass3d                          macro_scan, // new, will be filled out
			EyesisCorrectionParameters.CLTParameters clt_parameters,
			GeometryCorrection                       geometryCorrection,			  
			final double                             macro_disparity,
			final boolean                            show_corr_partial,
			final boolean                            show_corr_combo,
			final int                                debugLevel)
	{
		final int pTilesX = tp.getTilesX();
		final int pTilesY = tp.getTilesY();
		final int tileSize = tp.getTileSize(); // 
		final int mTilesX = (pTilesX + tileSize - 1) / tileSize;   // clt_aberrations_quad_corr truncates
		final int mTilesY = (pTilesY + tileSize - 1) / tileSize;
		
		int op = ImageDtt.setImgMask(0, 0xf);
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);
		final int    [][]    tile_op =         new int [mTilesY][mTilesX];
		final double [][]    disparity_array =  new double [mTilesY][mTilesX];
		for (int mty = 0; mty < mTilesY; mty ++){
			for (int mtx = 0; mtx < mTilesX; mtx ++){
				tile_op[mty][mtx] = op;
				disparity_array[mty][mtx] = macro_disparity;
			}
		}
//		final CLTPass3d macro_scan = new CLTPass3d(mtp);
		macro_scan.tile_op =   tile_op;
		macro_scan.disparity = disparity_array;

		return macro_scan;
	}

	public CLTPass3d  CLTMacroMeasure( // perform single pass according to prepared tiles operations and disparity
			final CLTPass3d                          macro_scan, // 
			final double [][][]                      input_data,
			EyesisCorrectionParameters.CLTParameters clt_parameters,
			GeometryCorrection                       geometryCorrection,			  
			final String                             suffix,
			final boolean                            show_corr_partial,
			final boolean                            show_corr_combo,
			final int                                debugLevel)
	{
		int mTilesY = macro_scan.tile_op.length;
		int mTilesX = macro_scan.tile_op[0].length;
		
		// undecided, so 2 modes of combining alpha - same as rgb, or use center tile only
		if (debugLevel > -1){
			int numTiles = 0;
			for (int ty = 0; ty < macro_scan.tile_op.length; ty ++) for (int tx = 0; tx < macro_scan.tile_op[ty].length; tx ++){
				if (macro_scan.tile_op[ty][tx] != 0) numTiles ++;
			}
			System.out.println("CLTMacroMeasure(): numTiles = "+numTiles);
		}
		double min_corr_selected = clt_parameters.min_corr;

		double [][] disparity_map = new double [ImageDtt.DISPARITY_TITLES.length][]; //[0] -residual disparity, [1] - orthogonal (just for debugging)

		double [][] shiftXY = {{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}};

		double [][][][] clt_corr_combo =    new double [ImageDtt.TCORR_TITLES.length][mTilesY][mTilesX][]; // needed always
		double [][][][][]   clt_corr_partial = null; // [tp.tilesY][tp.tilesX][pair][color][(2*transform_size-1)*(2*transform_size-1)]

		if (show_corr_partial) {
			clt_corr_partial = new double [mTilesY][mTilesX][][][];
			for (int i = 0; i < mTilesY; i++){
				for (int j = 0; j < mTilesX; j++){
					clt_corr_partial[i][j] = null;
				}
			}
		}

		if (show_corr_combo) {
			for (int i = 0; i < mTilesY; i++){
				for (int j = 0; j < mTilesX; j++){
					for (int k = 0; k<clt_corr_combo.length; k++){
						clt_corr_combo[k][i][j] = null;
					}
				}
			}
		}


		//		  double [][][][] texture_tiles =   save_textures ? new double [tilesY][tilesX][][] : null; // ["RGBA".length()][];
		ImageDtt image_dtt = new ImageDtt();
		image_dtt.clt_aberrations_quad_corr(
				8,                            // final int  macro_scale, // to correlate tile data instead of the pixel data: 1 - pixels, 8 - tiles
				macro_scan.tile_op,                      // per-tile operation bit codes
				macro_scan.disparity,         // clt_parameters.disparity,     // final double            disparity,
				input_data,                   // final double [][][]      imade_data, // first index - number of image in a quad
				// correlation results - final and partial          
				clt_corr_combo,               // [tp.tilesY][tp.tilesX][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				clt_corr_partial,             // clt_corr_partial,    // [tp.tilesY][tp.tilesX][quad]color][(2*transform_size-1)*(2*transform_size-1)] // if null - will not calculate
				null,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				//	Use it with disparity_maps[scan_step]?		  clt_mismatch,    // [tp.tilesY][tp.tilesX][pair]{dx,dy,weight}[(2*transform_size-1)*(2*transform_size-1)] // transpose unapplied. null - do not calculate
				disparity_map,    // [12][tp.tilesY * tp.tilesX]
				null, // [tp.tilesY][tp.tilesX]["RGBA".length()][]; 			  
				mTilesX * clt_parameters.transform_size, // imp_quad[0].getWidth(),       // final int width,
				clt_parameters.fat_zero,      // add to denominator to modify phase correlation (same units as data1, data2). <0 - pure sum
				clt_parameters.corr_sym,
				clt_parameters.corr_offset,
				clt_parameters.corr_red,
				clt_parameters.corr_blue,
				clt_parameters.corr_sigma,
				clt_parameters.corr_normalize, // normalize correlation results by rms 
				min_corr_selected, // 0.0001; // minimal correlation value to consider valid 
				clt_parameters.max_corr_sigma,// 1.5;  // weights of points around global max to find fractional
				clt_parameters.max_corr_radius,
				clt_parameters.enhortho_width,  // 2;    // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
				clt_parameters.enhortho_scale,  // 0.2;  // multiply center correlation pixels (inside enhortho_width)
				clt_parameters.max_corr_double, // Double pass when masking center of mass to reduce preference for integer values
				clt_parameters.corr_mode,     // Correlation mode: 0 - integer max, 1 - center of mass, 2 - polynomial
				clt_parameters.min_shot,       // 10.0;  // Do not adjust for shot noise if lower than
				clt_parameters.scale_shot,     // 3.0;   // scale when dividing by sqrt ( <0 - disable correction)
				clt_parameters.diff_sigma,     // 5.0;//RMS difference from average to reduce weights (~ 1.0 - 1/255 full scale image)
				clt_parameters.diff_threshold, // 5.0;   // RMS difference from average to discard channel (~ 1.0 - 1/255 full scale image)
				clt_parameters.diff_gauss,     // true;  // when averaging images, use gaussian around average as weight (false - sharp all/nothing)
				clt_parameters.min_agree,      // 3.0;   // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				clt_parameters.dust_remove,    // Do not reduce average weight when only one image differes much from the average
				clt_parameters.keep_weights,   // Add port weights to RGBA stack (debug feature)
				geometryCorrection,           // final GeometryCorrection  geometryCorrection,
				null,     // clt_kernels,                  // final double [][][][][][] clt_kernels, // [channel_in_quad][color][tileY][tileX][band][pixel] , size should match image (have 1 tile around)
				clt_parameters.kernel_step,
				clt_parameters.transform_size,
				clt_parameters.clt_window,
				shiftXY, // 
				0.0, // disparity_corr, // final double              disparity_corr, // disparity at infinity
				null, // (clt_parameters.fcorr_ignore? null: this.fine_corr),
				clt_parameters.corr_magic_scale, // still not understood coefficient that reduces reported disparity value.  Seems to be around 0.85 
				clt_parameters.shift_x,       // final int               shiftX, // shift image horizontally (positive - right) - just for testing
				clt_parameters.shift_y,       // final int               shiftY, // shift image vertically (positive - down)
				31, // clt_parameters.tileX,         // final int               debug_tileX,
				10, // clt_parameters.tileY,         // final int               debug_tileY,
				(clt_parameters.dbg_mode & 64) != 0, // no fract shift
				true,                        // no convolve
				mtp.getThreadsMax(),
				debugLevel - 2);
		macro_scan.disparity_map = disparity_map;
		macro_scan.is_measured =   true;
		macro_scan.is_combo =      false;
		macro_scan.resetProcessed();
		if (debugLevel > 1) {
			double [][]dbg_img = new double [input_data.length][];
			for (int i = 0; i < input_data.length; i++) {
				dbg_img[i] = input_data[i][0];
			}
			(new showDoubleFloatArrays()).showArrays(
					dbg_img,
					tp.getTilesX(),
					tp.getTilesY(),
					true,
					"DBG_INPUT_MACRO-D"+suffix);
		}

		if (debugLevel > 0) mtp.showScan(
				macro_scan, // CLTPass3d   scan,
				"macro_scan-D"+suffix); //String title)

		if (show_corr_partial) {
			String [] allColorNames = {"red","blue","green","combo"};
			String [] titles = new String[clt_corr_partial.length];
			for (int i = 0; i < titles.length; i++){
				titles[i]=allColorNames[i % allColorNames.length]+"_"+(i / allColorNames.length);
			}
			double [][] corr_rslt_partial = image_dtt.corr_partial_dbg(
					clt_corr_partial,
					2*clt_parameters.transform_size - 1,	//final int corr_size,
					4,	// final int pairs,
					4,    // final int colors,
					clt_parameters.corr_border_contrast,
					mtp.getThreadsMax(),
					debugLevel);
			(new showDoubleFloatArrays()) .showArrays(
					corr_rslt_partial,
					mTilesX*(2*clt_parameters.transform_size),
					mTilesY*(2*clt_parameters.transform_size),
					true,
					"MACRO-PART_CORR-D"+suffix,
					titles);
		}
		if (show_corr_combo) {
			double [][] corr_rslt = new double [clt_corr_combo.length][];
			String [] titles = new String[clt_corr_combo.length]; // {"combo","sum"};
			for (int i = 0; i< titles.length; i++) titles[i] = ImageDtt.TCORR_TITLES[i];
			for (int i = 0; i<corr_rslt.length; i++) {
				corr_rslt[i] = image_dtt.corr_dbg(
						clt_corr_combo[i],
						2*clt_parameters.transform_size - 1,
						clt_parameters.corr_border_contrast,
						mtp.getThreadsMax(),
						debugLevel);
			}

			(new showDoubleFloatArrays()).showArrays(
					corr_rslt,
					mTilesX*(2*clt_parameters.transform_size),
					mTilesY*(2*clt_parameters.transform_size),
					true,
					"MACRO-CORR-D"+suffix,
					titles );

		}
		return macro_scan;
	}
	
	public CLTPass3d refineMacro(
			final double [][][]                      input_data,
			EyesisCorrectionParameters.CLTParameters clt_parameters,
			GeometryCorrection                       geometryCorrection,			  
			final double                             trustedCorrelation,
			final double                             disp_far,   // limit results to the disparity range
			final double                             disp_near,
			final double                             minStrength,
			final double                             unique_tolerance,
			final int                                debugLevel)
	{
		CLTPass3d refined_macro = mtp.RefineQuadMulti(
				mtp.clt_3d_passes,           // final ArrayList <CLTPass3d> passes,
				0, // final int                   firstPass,
				mtp.clt_3d_passes.size(), // final int                   lastPassPlus1,
				trustedCorrelation, // final double                trustedCorrelation,
				disp_far, // final double                disp_far,   // limit results to the disparity range
				disp_near, // final double                disp_near,
				minStrength, // final double                minStrength,
				unique_tolerance, // final double                unique_tolerance,

				// TODO: when useCombo - pay attention to borders (disregard)
				false, // final boolean               usePoly,  // use polynomial method to find max), valid if useCombo == false
				1); // final int                   debugLevel)
		if (refined_macro == null) return null;
		
		CLTMacroMeasure( // perform single pass according to prepared tiles operations and disparity
				refined_macro,        // final CLTPass3d                          macro_scan, // 
				input_data,           // final double [][][]                      input_data,
				clt_parameters,       // EyesisCorrectionParameters.CLTParameters clt_parameters,
				geometryCorrection,   // GeometryCorrection                       geometryCorrection,			  
				""+mtp.clt_3d_passes, // final String                             suffix,
				false, // final boolean                            show_corr_partial,
				false, // final boolean                            show_corr_combo,
				0); // final int                                debugLevel)
		return refined_macro;
		
	}
	
	
	
	
}
