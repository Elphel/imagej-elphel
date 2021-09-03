package com.elphel.imagej.tileprocessor;
/**
 **
 ** ImageDttParameters - parameters defining TP operations (at first extra)
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ImageDttParameters.java is free software: you can redistribute it and/or modify
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

import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;

public class ImageDttParameters {
	public boolean gpu_mode_debug =         true;
	public boolean gpu_verify =             false; // verify tasks/ input data
	public boolean corr_mode_debug =        true;
	public boolean mix_corr_poly =          true;
	public boolean corr_poly_only =         false; // if LMA fails, discard tile
	public double  min_poly_strength =      0.2; /// 0.1
	public double  max_poly_hwidth =        2.5; // Maximal polynomial approximation half-width (in both directions)
	public double  poly_corr_scale =        2.0; /// 0.0 // Shift value if correlation maximum is wide in X than in Y to detect near objects (negative - far ones)

	public double  poly_pwr =               1.0;
	public double  poly_vasw_pwr =          2.0; // raise value to this power and apply as weight (0 - disable)
	public double  corr_magic_scale_cm =    1.0; /// 1.0 //0.85;  // reported correlation offset vs. actual one (not yet understood)
	public double  corr_magic_scale_poly =  1.0; /// 1.0 // 0.95;  // reported correlation offset vs. actual one (not yet understood)

	public int     ortho_height =             7;   // height of non-zero weights for hor/vert correlation to compensate borders
	public double  ortho_eff_height =         2.58; // effective correlation stripe height to match strengths

	public int     ortho_nsamples =           5; // number of samples to fit parabola
	public double  ortho_vasw_pwr =         2.0; // use data as weights when fitting parabola (high value samples are more important (when false use 3 samples only)

	private int     enhortho_width =        2;   // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
	private int     enhortho_width_aux =    1;   // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
	private double  enhortho_scale =        0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)
	private double  enhortho_scale_aux =    0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)
	public boolean ly_poly =                false; // Use polynomial when measuring mismatch (false - use center of mass)
	public double  ly_crazy_poly =          1.0; ///2.0  // Maximal allowed mismatch difference calculated as polynomial maximum
	public boolean ly_poly_backup =         true; // Use CM offset measuremets if poly failed

	public boolean fo_correct =             true; // correct far objects by comparing orthogonal correlations
	public boolean fo_far =                 false;// Use smallest disparity (false - largest)
	public double  fo_min_strength =        0.1;  // minimal strength for all correlations (full, hor, vert) to try to correct
	public double  fo_min_eff =             0.25; // minimal effective strength (corrected for width/hight)
	public double  fo_min_eff_ratio =       0.65; // minimal effective strength ration (this to orthogonal)
	public double  fo_max_hwidth =          3.0;  // maximal half-width disparity  direction to try to correct far objects
	public double  fo_min_diff =            0.2;  // minimal disparity difference to attempt extraction of foreground objects
	public boolean fo_ortho =               true; // require min/max to be orthogonal
	public double  fo_overcorrection =      0.0;  // add scaled hor/vert difference to the largest of hor/vert disparity
	public double  fo_lim_overcorr =        2.0;  // limit full correction with respect to largest - fullcorr difference

	public double  mismatch_max_diff =      1.0;  // maximal disparity difference to calculate mismatch. Too low value will disqualify
	                                              // large mismatches, too high - treat multi-z disparity as mismatch

	public double  corr_offset =            0.1; //0.1;  // add to pair correlation before multiplying by other pairs
	public boolean twice_diagonal =         true;  // twice_diagonal diagonal pairs provide twice denser samples, when true it doubles the weight of diagonal pairs
	public double  min_corr =               0.02; // minimal correlation value to consider valid

	public int     dbg_pair_mask =          0x3f;  // which pairs to combine
	public int     corr_strip_hight =       9;     // number of rows to calculate

	// Extracting bi-convex (convex in both orthogonal directions) cells and allowing non-convex on the selection border only
	public int     cnvx_hwnd_size =         4;     // half window size (both horizontal and vertical to extract bi-convex cells
	public double  cnvx_weight =            0.5;   // relative weight of non-convex (border) cell
	public boolean cnvx_add3x3 =            true;  // always select 3x3 cells around integer maximum
	public int     cnvx_min_samples =       5;     // minimal number of used samples on a 2D correlation pair
	public boolean cnvx_non_coll =          true;  // require non-collinear samples in a valid correlation pair
	public int     cnvx_min_pairs =         5;     // minimal number 2D correlation pairs per tile

	public int     mcorr_multi =            4;     // apply *_multi if number of sensors >= mcorr_multi
	public boolean mcorr_all =              true;  // correlate all N*(N-1)/2 pairs
	public boolean mcorr_all_multi =        false; // correlate all N*(N-1)/2 pairs for multi
	public boolean mcorr_dia =              true;  // all diameters (2 pairs for quad, 8 - for lwir16)
	public boolean mcorr_dia_multi =        true;  // all diameters
	public boolean mcorr_sq =               true;  // all squares (4 pairs for quad, 16 - for lwir16)
	public boolean mcorr_sq_multi =         true;  // all squares
	public boolean mcorr_neib =             true;  // all neighbors (4 pairs for quad, 16 - for lwir16)
	public boolean mcorr_neib_multi =       true;  // all neighbors
	public boolean mcorr_hor =              true;  // all horizontal (2 pairs for quad, 7 - for lwir16)
	public boolean mcorr_hor_multi =        true;  // all horizontal
	public boolean mcorr_vert =             true;  // all vertical (2 pairs for quad, 8 - for lwir16)
	public boolean mcorr_vert_multi =       true;  // all vertical

	
	public boolean mcorr_cons_all =         true;  // consolidate all available pairs
	public boolean mcorr_cons_dia =         true;  // consolidate all having length of a diameter
	public boolean mcorr_cons_sq =          true;  // consolidate all having length of a square side
	public boolean mcorr_cons_neib =        true;  // consolidate all having shortest length
	public boolean mcorr_cons_hor =         true;  // consolidate all horizontal pairs
	public boolean mcorr_cons_vert =        true;  // consolidate all vertical pairs
	
	/// these are just for testing, actual will be generated specifically for different applications (LY will use closest pairs)  
	public int     mcorr_comb_width =       15;
	public int     mcorr_comb_height =      15;
	public int     mcorr_comb_offset =      0;
	public double  mcorr_comb_disp =        1.0;   // per-pixel disparity relative to diameter
	public boolean mcorr_comb_dbg =         false; // Generate/show correlation pairs scaled/rotated for combining
	
	// pole window will be inverted
	public int     corr_strip_notch =       11;     // number of rows to calculate for vertical poles
	public double  corr_notch_hwidth =      2.0; // 6.0;   // 50% window cutoff height
	public double  corr_notch_blur =        2.0; // 5.0;   // 100% to 0 % vertical transition range


	// Used for common interpolated stripes
	public int     corr_wndy_size =         2;   // 9;     // number of rows to calculate CM disparity
	public double  corr_wndy_hwidth =       1.0; // 6.0;   // 50% window cutoff height
	public double  corr_wndy_blur =         1.0; // 5.0;   // 100% to 0 % vertical transition range
	public int     corr_wndx_size =         9;     // half number of columns to calculate CM disparity (each row has only odd/even columns, so disparity range is smaller
	public double  corr_wndx_hwidth =       6.0;   // 50% window cutoff width
	public double  corr_wndx_blur =         5.0;   // 100% to 0 % vertical transition range
	
	public boolean pcorr_use =              false; // use group phase correlation for disparity calculation in existing code
	public boolean pcorr_use_hv =           true;  // use group phase correlation to combine hor and vert pairs
	public double  pcorr_sigma_mono =       0.15;  // correlation sigma for monochrome images  (after normalization)
	public double  pcorr_sigma =            0.8;  // correlation sigma for Bayer images  (after normalization)
	public double  pcorr_sigma_rb =         0.5;  // correlation sigma extra for R and B (before normalization)
	public double  pcorr_fat_zero =         0.05; // correlation relative fat zero for Bayer images
	public double  pcorr_fat_zero_mono =    0.03; // correlation relative fat zero for monochrome images
	public double  pcorr_dbg_offsx =        0.0;  // X-offset correlation in TD
	public double  pcorr_dbg_offsy =        0.0;  // Y-offset correlation in TD

// LMA parameters
	public double  lma_disp_range =          5.0;   // disparity range to combine in one cluster (to mitigate ERS
// LMA single parameters
	public boolean lmas_gaussian =           false; // model correlation maximum as a Gaussian (false - as a parabola)
	public boolean lmas_adjust_wm =          true;  // used in new for width
	public boolean lmas_adjust_wy =          true;  // adjust non-circular
	public boolean lmas_adjust_ag =          true;  // adjust gains gains
	// Pre-lma poly
	public double  lmas_poly_str_scale =     1.0;   // scale pre-lma poly strength
	public double  lmas_poly_str_min =       0.05;  // ignore tiles with poly strength (scaled) below

	public double  lmas_lambda_initial =     0.03;   //
	public double  lmas_rms_diff =           0.0003; //
	public int     lmas_num_iter =          20;     ///10
	// Filtering and strength calculation
	public double  lmas_max_rel_rms =        0.3;   // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
	public double  lmas_min_strength =       0.7;   // minimal composite strength (sqrt(average amp squared over absolute RMS)
	public double  lmas_min_ac =             0.02;  // minimal of a and C coefficients maximum (measures sharpest point/line)
	public double  lmas_min_min_ac =         0.007; // minimal of a and C coefficients minimum (measures sharpest point)
	public double  lmas_max_area =           0.0;  // maximal half-area (if > 0.0)

	public boolean lma_gaussian =           false; // model correlation maximum as a Gaussian (false - as a parabola)
	public boolean lma_second =             true;  // re-run LMA after removing weak/failed tiles
	public boolean lma_second_gaussian =    false;  // re-run after removing weal/failed in Gaussian mode
	public boolean lma_adjust_wm =          true;  // used in new for width
	public boolean lma_adjust_wy =          true;  // false; // used in new for ellipse
	public boolean lma_adjust_wxy =         true;  // used in new for lazy eye adjust parallel-to-disparity correction
	public boolean lma_adjust_ly1 =         true;  // adjust ortho-to-disparity correction
	public boolean lma_adjust_ag =          true;  // used in new for gains

// new LMA parameters
	public double  lma_wnd =                1.0; //1.5;   // raise cosine window to this power (1.0 - just 2D cosine)
	public double  lma_min_wnd =            0.4;   // divide values by the 2D correlation window if it is >= this value for finding maximums and convex areas
	public double  lma_wnd_pwr =            0.8;   // Raise window for finding a maximum and a convex region to this power
	public int     lma_hard_marg =          1;     // Zero out this width margins before blurring
	public int     lma_soft_marg =          2;     // Do not look for maximums inside this width margins
	public double  lma_sigma =              0.7;   // Blur correlation before finding maximum and convex region



	// maybe try using sqrt (corr_wnd) ? or variable power?


	public double  lma_half_width =         2.0;   //
	public double  lma_cost_wy =            0.0;  /// 0.00050 // cost of parallel-to-disparity correction
	public double  lma_cost_wxy =           0.0;  /// 0.00100 // cost of ortho-to-disparity correction

	public double  lma_lambda_initial =     0.03;   //
	public double  lma_lambda_scale_good =  0.5;   //
	public double  lma_lambda_scale_bad =   8.0;   //
	public double  lma_lambda_max =       100.0;   //
	public double  lma_rms_diff =           0.003; //
	public int     lma_num_iter =          10;     //
	// Filtering and strength calculation
	public double  lma_max_rel_rms =        0.25;  // maximal relative (to average max/min amplitude LMA RMS) // May be up to 0.3)
	public double  lma_min_strength =       1.0;   // minimal composite strength (sqrt(average amp squared over absolute RMS)
	public double  lma_min_ac =             0.05;  // minimal of a and C coefficients maximum (measures sharpest point/line)
	public double  lma_min_min_ac =         0.015; // minimal of a and C coefficients minimum (measures sharpest point)
	public double  lma_max_area =           30.0;  //45.0;  // maximal half-area (if > 0.0)

	public double  lma_str_scale =          0.2;   // convert lma-generated strength to match previous ones - scale
	public double  lma_str_offset =         0.05;  // convert lma-generated strength to match previous ones - add to result


	// Lazy eye results interpretation
	public boolean lma_diff_xy =            true;  // convert dd/nd to x,y
	public double  lma_diff_minw =          0.07;   // minimal weight to keep
	public double  lma_diff_sigma =         2.0;   // blur differential data (relative to the cluster linear size)

	public int     lma_debug_level =        0;     //
	public int     lma_debug_level1 =       0;     //
	public boolean lma_debug_graphic =      false; //
	@Deprecated
	public boolean corr_var_cam =           true;  // New correlation mode compatible with 8 subcameras 
	public double  cm_max_normalization =   0.55; // fraction of correlation maximum radius, being squared multiplied by maximum to have the same total mass
	public double  lma_dbg_scale =          0.0;  // scale lma_dbg_offset
	public double  [][] lma_dbg_offset =    {{ 1.0, 0.0},{ -1.0, 0.0},{-1.0, 0.0},{ 1.0, 0.0}}; //  new double [4][2];

	public double getCorrSigma(boolean monochrome) {
		return monochrome ? pcorr_sigma_mono : pcorr_sigma;
	}
	public double getFatZero(boolean monochrome) {
		return monochrome ? pcorr_fat_zero_mono : pcorr_fat_zero;
	}
	
	public int getEnhOrthoWidth(boolean aux) {
		return aux ? enhortho_width_aux : enhortho_width;
	}

	public double getEnhOrthoScale(boolean aux) {
		return aux ? enhortho_scale_aux : enhortho_scale;
	}

	// next 2 only used to read old config files
	public void setEnhOrthoWidth (int w) {
		enhortho_width = w;
	}
	public void setEnhOrthoScale (double s) {
		enhortho_scale = s;
	}
	
	
	public boolean getMcorrAll(int numSens) {
		return (numSens > mcorr_multi) ? mcorr_all_multi : mcorr_all;
	}

	public boolean getMcorrDia(int numSens) {
		return (numSens > mcorr_multi) ? mcorr_dia_multi : mcorr_dia;
	}

	public boolean getMcorrSq(int numSens) {
		return (numSens > mcorr_multi) ? mcorr_sq_multi :  mcorr_sq;
	}

	public boolean getMcorrNeib(int numSens) {
		return (numSens > mcorr_multi) ? mcorr_neib_multi : mcorr_neib;
	}

	public boolean getMcorrHor(int numSens) {
		return (numSens > mcorr_multi) ? mcorr_hor_multi : mcorr_hor;
	}

	public boolean getMcorrVert(int numSens) {
		return (numSens > mcorr_multi) ? mcorr_vert_multi : mcorr_vert;
	}

	public void dialogQuestions(GenericJTabbedDialog gd) {
		
		    gd.addCheckbox    ("Debug CPU->GPU matching",                                         this.gpu_mode_debug,
				"output clt_corr_partial");
		    gd.addCheckbox    ("Verify GPU input",                                                this.gpu_verify,
				"Check tp_tasks and fix NaN");
			gd.addCheckbox    ("Enable ImageDtt correlation debug layers",                        this.corr_mode_debug,
					"false - return (old) per-coord correlations, true - replace them with more pairs correlation (new)");
			gd.addCheckbox    ("Replace CM layer with mixed/new poly one",                        this.mix_corr_poly);
			gd.addCheckbox    ("If LMA fails, discard tile",                                      this.corr_poly_only);
			gd.addNumericField("Use poly mode if strength is greater than",                       this.min_poly_strength,  3,6,"", "AND condition");
			gd.addNumericField("Maximal polynomial approximation half-width",                     this.max_poly_hwidth,  3,6,"pix", "Maximal polynomial approximation half-width (in both directions), Most now are ~2.0");
			gd.addNumericField("Polynomial argmax correction (positive - near, negative - far)",  this.poly_corr_scale,  3,6,"×", "Shift value if correlation maximum is wide in X than in Y to detect near objects (negative - far ones)");

			gd.addNumericField("When calculating poly raise correlation value to this power",     this.poly_pwr,  3,6,"", "Trying to reduce sticking to integer values");
			gd.addNumericField("When calculating poly multiply weight by correlation value",      this.poly_vasw_pwr,3,6,"",
					"Raise value to this power and apply as weight");
			gd.addNumericField("Magic scale for CM correlation",                                  this.corr_magic_scale_cm,  3,6,"",
					"Reported center of mass correlation value is this part of actual");
			gd.addNumericField("Magic scale for Poly correlation",                                this.corr_magic_scale_poly,  3,6,"",
					"Reported polynomial correlation value is this part of actual");

			gd.addNumericField("Height of non-zero weights for hor/vert correlation",             this.ortho_height,  3);
			gd.addNumericField("Effective correlation stripe height to match strengths",          this.ortho_eff_height,  3,6,"", "For matching 2-d and 1-d correlation maximums");

			gd.addNumericField("Number of samples to fit parabola in ortho mode",                 this.ortho_nsamples,  3);
			gd.addNumericField("Use data as weights when fitting parabola for ortho mode",        this.ortho_vasw_pwr,3,6,"",
					"Raise value to this power and apply as weight.  Reduce width to 3 samples if false, 5 OK when true");

			gd.addNumericField("Reduce weight of center correlation pixels from center for the main camera",this.enhortho_width, 0, 3, "pix",
					"Reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)");
			gd.addNumericField("Reduce weight of center correlation pixels from center for the aux camera", this.enhortho_width_aux, 0, 3, "pix",
					"Reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)");
			gd.addNumericField("Multiply center correlation pixels for the main camera",                   this.enhortho_scale, 0,3, "pix",
					"Multiply center correlation pixels (inside enhortho_width) (1.0 - disables enh_ortho) - main camera");
			gd.addNumericField("Multiply center correlation pixels for the aux camera",                    this.enhortho_scale_aux, 0,3, "pix",
					"Multiply center correlation pixels (inside enhortho_width) (1.0 - disables enh_ortho) - aux camera");

			gd.addCheckbox    ("Use polynomial when measuring mismatch (false - use center of mass)",      this.ly_poly);
			gd.addNumericField("Maximal allowed mismatch difference calculated as polynomial maximum",     this.ly_crazy_poly,3,6,"px",
					"Use CM method as a back-up, if polynomial maximum is far off");
			gd.addCheckbox    ("Use CM offset measuremets if poly failed",                                  this.ly_poly_backup);

			gd.addMessage("Far objects correction");
			gd.addCheckbox    ("Try to correct far objects (make them closer) by hor/vert comparison",      this.fo_correct);
			gd.addCheckbox    ("Use smallest disparity/foreground objects (false - largest)",               this.fo_far,
					"When disparity is different for different directios, use smallest (longest distance)");
			gd.addNumericField("Minimal strength for all correlations (full, hor, vert) to try to correct", this.fo_min_strength,  3,6,"",
					"Do not correct if any of the full, hor or vert strength is lower than this");

			gd.addNumericField("Minimal effective strength (corrected for width/hight)",                    this.fo_min_eff,  3,6,"",
					"Minimal strength for the direction divided by width in disparity direction and multiplied by width in orthogonal direction");
			gd.addNumericField("Minimal ratio of the effective strength to that of the orthogonal one",     this.fo_min_eff_ratio,  3,6,"",
					"Ratio of the effective strength in the disparity max/min direction to that of the orthogonal pair(s)");

			gd.addNumericField("Maximal correlation half-width",                                            this.fo_max_hwidth,  3,6,"pix",
					"Do not correct if any of the full, hor or vert half-width is above this");
			gd.addNumericField("Minimal disparity difference to attempt extraction of foreground objects",  this.fo_min_diff,  3,6,"pix",
					"Keep all-pair disparity if directional difference is below this value");
			gd.addCheckbox    ("Require min/max to be orthogonal",                                          this.fo_ortho,
					"Correct fore/background only if the smallest/largest disparity is for orthogonal directions");
			gd.addNumericField("Overcorrection scale",                                            this.fo_overcorrection,  3,6,"",
					"Add scaled hor/vert difference to the largest of hor/vert disparity. Use 0 to just select largest of disparities");
			gd.addNumericField("Limit overcorrection",                                            this.fo_lim_overcorr,  3,6,"",
					"Limit full correction with respect to largest - fullcorr difference. 1.0 does not allow overcorrection, < 1.0 - the result will be closer to full correction");

			gd.addNumericField("Maximal disparity difference to calculate mismatch",             this.mismatch_max_diff,  3,6,"",
					"Too low value will disqualify large mismatches, too high - treat multi-z disparity as mismatch");

	  		gd.addNumericField("Add to pair correlation before multiplying by other pairs",       this.corr_offset,  6,8,"",
	  				"0.0 - pure product (false-positive tolerant), higher the value - closer to the sum (linear, better S/N");

			gd.addCheckbox    ("Double weight of diagonal pairs when combining",                  this.twice_diagonal,
					"Diagonal pairs provide twice denser samples, when true it doubles the weight of diagonal pairs");

		    gd.addNumericField("Extract disparity max/argmax if maximal value is above",          this.min_corr,  6,8,"",
		    		"skip unreliable correlations");


			gd.addNumericField("Debug: which pairs to combine",                                   this.dbg_pair_mask,  0, 3, "",
					"Bits: 0, 1 - horizontal pairs, 2,3 - verical pairs, 4,5 - diagonal pairs");
			gd.addNumericField("Number of correlation rows to combine (strip height)",            this.corr_strip_hight,  0, 3, "",
					"Number of rows to combine/interpolate correlation results. Rows are twice denser than pixels correponding to largest baseline disparity");


			gd.addNumericField("Half window size to extract bi-convex cells",                     this.cnvx_hwnd_size,  0, 3, "pix",
					"Create selection mask for quadratic approximation inside square around initial maximum position, specify distance from the center");
		    gd.addNumericField("Relative weight of non-convex (border) cell",                     this.cnvx_weight,  6,8,"",
		    		"Weight of the bi-convex points is 1.0, this value specifies reduced weight of the border (non-bi-convex) points");
			gd.addCheckbox    ("Always select 3x3 cells around integer maximum",                  this.cnvx_add3x3,
					"Add 3x3 cells selection around the original argmax, regardless of bi-convex property");

			gd.addNumericField("Min samples per pair",                                            this.cnvx_min_samples,  0, 3, "pix",
					"Minimal number of used samples on a 2D correlation pair");
			gd.addCheckbox    ("Always select 3x3 cells around integer maximum",                  this.cnvx_non_coll,
					"Must have non-collinear samples in each correlation pair");
			gd.addNumericField("Minimal correlation pairs per tile",                              this.cnvx_min_pairs,  0, 3, "",
					"Minimal number 2D correlation pairs per tile");

			gd.addMessage("Correlation pairs selection for multicamera confifurations");
			gd.addNumericField("Min number of sensors for \"multi\"",                             this.mcorr_multi,  0, 3, "sensors",
					"Cameras with this or larger number of sensor will use a \"_multi\" version of parameters (when available)");

			gd.addCheckbox    ("Calculate all correlation pairs for small cameras",               this.mcorr_all,
					"N*(N-1)/2: 6 for quad, 120 for lwir16 ");
			gd.addCheckbox    ("Calculate all correlation pairs for multi cameras",               this.mcorr_all_multi,
					"N*(N-1)/2: 6 for quad, 120 for lwir16");

			gd.addCheckbox    ("Calculate correlation pairs - diameters for small cameras",       this.mcorr_dia,
					"All diameter pairs. N/2: 2 for quad, 8 for lwir16");
			gd.addCheckbox    ("Calculate correlation pairs - diameters for multi cameras",       this.mcorr_dia_multi,
					"All diameter pairs. N/2: 2 for quad, 8 for lwir16");

			gd.addCheckbox    ("Calculate correlation pairs of squares for small cameras",        this.mcorr_sq,
					"All pairs with the lengths equal to side of a square. N: 4 for quad, 16 for lwir16");
			gd.addCheckbox    ("Calculate correlation pairs of squares for multi cameras",        this.mcorr_sq_multi,
					"All pairs with the lengths equal to side of a square. N: 4 for quad, 16 for lwir16");
			
			gd.addCheckbox    ("Calculate neighbor pairs for small cameras",                     this.mcorr_neib,
					"All neighbor pairs. N: 4 for quad, 16 for lwir16");
			gd.addCheckbox    ("Calculate neighbor pairs for multi cameras",                     this.mcorr_neib_multi,
					"All neighbor pairs. N: 4 for quad, 16 for lwir16");
			
			gd.addCheckbox    ("Calculate horizontal pairs for small cameras",                   this.mcorr_hor,
					"All horizontal pairs. N: 2 for quad, 7 for lwir16");
			gd.addCheckbox    ("Calculate horizontal pairs for multi cameras",                   this.mcorr_hor_multi,
					"All horizontal pairs. N: 2 for quad, 7 for lwir16");

			gd.addCheckbox    ("Calculate vertical pairs for small cameras",                     this.mcorr_vert,
					"All vertical pairs. N: 2 for quad, 8 for lwir16");
			gd.addCheckbox    ("Calculate vertical pairs for multi cameras",                     this.mcorr_vert_multi,
					"All vertical pairs. N: 2 for quad, 8 for lwir16");

			

			gd.addCheckbox    ("Combine all pairs",                                               this.mcorr_cons_all,
					"Combine all calculated correlation pairs");
			gd.addCheckbox    ("Combine diameters",                                               this.mcorr_cons_dia,
					"Combine all pairs with the length of a diameter");
			gd.addCheckbox    ("Combine squares",                                                 this.mcorr_cons_sq,
					"Combine all pairs with the length of a side of an inscribed square");
			gd.addCheckbox    ("Combine neighbors",                                               this.mcorr_cons_neib,
					"Combine all neighbor pairs ");
			gd.addCheckbox    ("Combine horizontal pairs",                                        this.mcorr_cons_hor,
					"Combine all calculated horizontal pairs regardless of length");
			gd.addCheckbox    ("Combine vertical pairs",                                          this.mcorr_cons_vert,
					"Combine all calculated vertical pairs regardless of length");

			gd.addMessage("Generating grid for combining visualization, actual will be provided programmatically");
			gd.addNumericField("Width of a combined correlation tile",                            this.mcorr_comb_width,  0, 3, "pix",
					"Width of a tile to combine correlations after rotation/scaling");
			gd.addNumericField("Height of a combined correlation tile",                           this.mcorr_comb_height,  0, 3, "pix",
					"Full height of a tile to combine correlations after rotation/scaling");
			gd.addNumericField("Height offset of a combined correlation tile",                    this.mcorr_comb_offset,  0, 3, "pix",
					"0 - centered (-height/2 to height/2), height/2 - only positive (0 to height)");
			gd.addNumericField("Relative per-pixel disparity",                                    this.mcorr_comb_disp,  3,6,"pix",
					"Combined tile per-pixel disparity for baseline == diameter");
			gd.addCheckbox    ("Debug correlation scaling/rotation for combining",                this.mcorr_comb_dbg,
					"Generate/show correlation pairs scaled/rotated for combining");
			
			gd.addMessage("Window for pole detection mode");
			gd.addNumericField("Strip height for pole detection",                                 this.corr_strip_notch,  0, 3, "half-pix",
					"Number of rows to combine/interpolate correlation results. Rows are twice denser than pixels correponding to largest baseline disparity");

			gd.addNumericField("50% correlation window cutoff height for poles (0 in the center)",this.corr_notch_hwidth,  3, 6, "half-pix",
					"Correlation window height argument for 50% value");
			gd.addNumericField("0% to 100 % transition range for poles",                          this.corr_notch_blur,  3,6,"half-pix",
					"Transition range, shifted sine is used");

			gd.addMessage("Window for normal correlations");
			gd.addNumericField("Number of rows to calculate CM disparity",                        this.corr_wndy_size,  0, 3, "",
					"Number of rows to calculate maximum. Normally should be equal to the previous parameter");

			gd.addNumericField("50% correlation window cutoff height",                            this.corr_wndy_hwidth,  3, 6, "",
					"Correlation window height argument for 50% value");
			gd.addNumericField("100% to 0 % correlation vertical window transition range",        this.corr_wndy_blur,  3,6,"",
					"Transition range, shifted sine is used");

			gd.addNumericField("Half-number of columns to calculate CM disparity",                this.corr_wndx_size,  0, 3, "",
					"Horizontal size of the window (symmetrical around integer argmax, only odd/even values are used in odd/even rows");

			gd.addNumericField("50% correlation window cutoff width",                             this.corr_wndx_hwidth,  3, 6, "",
					"Correlation window width argument for 50% value");
			gd.addNumericField("100% to 0 % correlation horizontal window transition range",      this.corr_wndx_blur,  3, 6, "",
					"Transition range, shifted sine is used");

			gd.addTab("Corr Intra","Parameters Group 2D Phase correlation");

			gd.addCheckbox    ("Use group phase correlation for disparity calculation",                this.pcorr_use,
					"Default false, for compatibility with existing code/ old configs");
			gd.addCheckbox    ("Use group phase correlation to combine hor and vert pairs separately", this.pcorr_use_hv,
					"Use even when hor are not combined with vert");
			gd.addNumericField("Correlation sigma for monochrome images",                         this.pcorr_sigma_mono,  3, 6, "pix",
					"after normalization");
			gd.addNumericField("Correlation sigma for Bayer images",                              this.pcorr_sigma,  3, 6, "pix",
					"after normalization");
			gd.addNumericField("Extra correlation sigma for red and blue",                        this.pcorr_sigma_rb,  3, 6, "pix",
					"before normalization");
			gd.addNumericField("Correlation relative fat zero for Bayer images",                  this.pcorr_fat_zero,  3, 6, "",
					"Normalized to average correlation absolute value");
			gd.addNumericField("Correlation relative fat zero for monochrome images",             this.pcorr_fat_zero_mono,  3, 6, "",
					"Normalized to average correlation absolute value");
			gd.addNumericField("Debug feature - shift correlation result X",                      this.pcorr_dbg_offsx,  3, 6, "",
					"Rotate in Transform Domain");
			gd.addNumericField("Debug feature - shift correlation result Y",                      this.pcorr_dbg_offsy,  3, 6, "",
					"Rotate in Transform Domain");
			
			gd.addTab("Corr LMA","Parameters for LMA fitting of the correlation maximum parameters");
			gd.addMessage("Single-tile (no lazy eye) only parameters (some are common");
			gd.addNumericField("Cluster disparity range",                                         this.lma_disp_range,  3, 6, "pix",
					"Disparity range to combine in one cluster (to mitigate ERS");
			gd.addCheckbox    ("Correlation maximum as gaussian",                                 this.lmas_gaussian,
					"Model correlation maximum as a Gaussian exp(-r^2)  (false - as a parabola - 1-r^2)");
			gd.addCheckbox    ("Fit correlation defined half-width",                              this.lmas_adjust_wm,
					"Allow fitting of the half-width common for all pairs, defined by the LPF filter of the phase correlation");
			gd.addCheckbox    ("Adjust ellipse parameters (was Fit extra vertical half-width)",   this.lmas_adjust_wy,
					"Adjust ellipse (non-circular) of the correlation maximum (was Fit extra perpendicular to disparity half-width (not used? and only possible with multi-baseline cameras))");
			gd.addCheckbox    ("Adjust per-pair scale (was Adjust per-group amplitudes)",         this.lmas_adjust_ag,
					"Each correlation pair gain (was Each correlation type's amplitude (now always needed))");

			gd.addMessage("pre-LMA (polynomial) filtering");
			gd.addNumericField("Scale pre-LMA poly strength",                                     this.lmas_poly_str_scale,  3, 6, "",
					"Calculated as maximal value over average radius");
			gd.addNumericField("Minimal pre-LMA poly strength (scaled)",                          this.lmas_poly_str_min,  3, 6, "",
					"Ignore tiles with pre-LMA poly strength (scaled with above) below this value");
		    gd.addMessage("LMA (single) LMA fitting parameters");
		    gd.addNumericField("Initial value of LMA lambda",                                     this.lmas_lambda_initial,  3, 6, "",
		            "The higher the lambda the more close it will be to the gradient descent (slower/safer)");
		    gd.addNumericField("Relative RMS improvement to exit LMA",                            this.lmas_rms_diff,  6, 8, "",
		            "LMA will report success when realtive RMS improvements fall below this value");
		    gd.addNumericField("LMA maximal iterations",                                          this.lmas_num_iter,  0, 3, "",
		            "Limit LMA cycles, so it will exit after certain number of small improvements");
		    gd.addMessage("LMA (single) results filtering");
		    gd.addNumericField("Maximal relative RMS ",                                           this.lmas_max_rel_rms,  6, 8, "",
		            "Discard tile if ratio of RMS to average of min and max amplitude exceeds this value");
		    gd.addNumericField("Minimal composite strength",                                      this.lmas_min_strength,  6, 8, "",
		            "Discard tile if composite strength (average amplitude over SQRT of RMS) is below");
		    gd.addNumericField("Minimal max (A,C)",                                               this.lmas_min_ac,  6, 8, "",
		            "Minimal value of max (A,C) coefficients to keep the tile (measures sharpest point/line correlation maximum)");
			gd.addNumericField("Minimal min (A,C)",                                               this.lmas_min_min_ac,  6, 8, "",
					"Minimal value of min (A,C) coefficients to keep the tile (measures sharpest point correlation maximum)");
		    gd.addNumericField("Maximal area",                                                    this.lmas_max_area,  6, 8, "sq.pix",
		            "Maximal product of maximum half-width by half-height, ignore check if <=0");

			gd.addMessage("Multi-tile (for lazy eye) LMA (some are used for with single-tile mode too)");
			gd.addCheckbox    ("Correlation maximum as gaussian",                                 this.lma_gaussian,
					"Model correlation maximum as a Gaussian exp(-r^2)  (false - as a parabola - 1-r^2)");
			gd.addCheckbox    ("Re-run LMA after removing weak/failed tiles",                     this.lma_second,
					"Re-run LMA with filtered tiles (see Correlation strength calculation section below)");
			gd.addCheckbox    ("Gaussian mode during LMA re-run",                                 this.lma_second_gaussian,
					"Parabola is more stable when using with un-filtered tiles, so it makes sense to use Gaussina only on filtered tiles");

			gd.addCheckbox    ("Fit correlation defined half-width",                              this.lma_adjust_wm,
					"Allow fitting of the half-width common for all pairs, defined by the LPF filter of the phase correlation");
			gd.addCheckbox    ("Adjust ellipse parameters (was Fit extra vertical half-width)",   this.lma_adjust_wy,
					"Adjust ellipse (non-circular) of the correlation maximum (was Fit extra perpendicular to disparity half-width (not used? and only possible with multi-baseline cameras))");
			gd.addCheckbox    ("Adjust \"lazy eye\" parameters parallel to disparity",            this.lma_adjust_wxy,
					"(was Fit extra half-width along disparity) Increased width in disparity direction caused by multi-distance objects in the tile");
			gd.addCheckbox    ("Adjust \"lazy eye\" parameters orthogonal CW to disparity",       this.lma_adjust_ly1,
					"Increased width in disparity direction caused by multi-distance objects in the tile");

			gd.addCheckbox    ("Adjust per-pair scale (was Adjust per-group amplitudes)",         this.lma_adjust_ag,
					"Each correlation pair gain (was Each correlation type's amplitude (now always needed))");

			gd.addNumericField("LMA window power",                                                this.lma_wnd,  3, 6, "",
					"Raise cosine window to this power (1.0 - plane 2D cosine");
			gd.addNumericField("Minimal window value for normalization during max/convex",        this.lma_min_wnd,  3, 6, "",
					"divide values by the 2D correlation window if it is >= this value for finding maximums and convex areas");
			gd.addNumericField("LMA window power for convex region",                              this.lma_wnd_pwr,  3, 6, "",
					"Raise window for finding a maximum and a convex region to this power");
			gd.addNumericField("LMA hard margin",                                                 this.lma_hard_marg,  0, 3, "",
					"Zero out this width margins before blurring");
			gd.addNumericField("LMA soft margins iterations",                                     this.lma_soft_marg,  0, 3, "",
					"Do not look for maximums inside this width margins");
			gd.addNumericField("LMA blur sigma",                                                  this.lma_sigma,  3, 6, "",
					"Blur correlation before finding maximum and convex region");


			gd.addNumericField("Initial/expected half-width of the correlation maximum in both directions", this.lma_half_width,  3, 6, "pix",
					"With LPF sigma = 0.9 it seems to be ~= 2.0. Used both as initial parameter and the fitted value difference from this may be penalized");
			gd.addNumericField("Lazy eye cost parallel to disparity (was Cost of the difference of the actual half-width...)",  this.lma_cost_wy,  5, 8, "",
					"The higher this cost, the more close the fitted half-width will be to the expected one");
			gd.addNumericField("Lazy eye cost ortho to disparity (was Cost of hor / vert widths difference)",  this.lma_cost_wxy,  5, 8, "",
					"Tries to enforce equal width and hight of the correlation maximum");

			gd.addNumericField("Initial value of LMA lambda",                                     this.lma_lambda_initial,  3, 6, "",
					"The higher the lambda the more close it will be to the gradient descent (slower/safer)");
			gd.addNumericField("Scale (decrease) LMA lambda on success (usually 0.5)",            this.lma_lambda_scale_good,  3, 6, "",
					"Make it smaller for unsafe but faster converging (if all goes well)");
			gd.addNumericField("Scale (increase) LMA lambda after failure (usually 8.0)",         this.lma_lambda_scale_bad,  3, 6, "",
					"Bad convergence usually means errors in Jacobian (derivatives)");
			gd.addNumericField("Lambda value to give up increasing it",                           this.lma_lambda_max,  3, 6, "",
					"Gives up LMA if increased lambda still does not improve RMS");
			gd.addNumericField("Relative RMS improvement to exit LMA",                            this.lma_rms_diff,  6, 8, "",
					"LMA will report success when realtive RMS improvements fall below this value");

			gd.addNumericField("LMA maximal iterations",                                          this.lma_num_iter,  0, 3, "",
					"Limit LMA cycles, so it will exit after certain number of small improvements");

			gd.addMessage("LMA results filtering");
			gd.addNumericField("Maximal relative RMS ",                                           this.lma_max_rel_rms,  6, 8, "",
					"Discard tile if ratio of RMS to average of min and max amplitude exceeds this value");
			gd.addNumericField("Minimal composite strength",                                      this.lma_min_strength,  6, 8, "",
					"Discard tile if composite strength (average amplitude over SQRT of RMS) is below");
			gd.addNumericField("Minimal max (A,C)",                                               this.lma_min_ac,  6, 8, "",
					"Minimal value of max (A,C) coefficients to keep the tile (measures sharpest point/line correlation maximum)");
			gd.addNumericField("Minimal min (A,C)",                                               this.lma_min_min_ac,  6, 8, "",
					"Minimal value of min (A,C) coefficients to keep the tile (measures sharpest point correlation maximum)");
			gd.addNumericField("Maximal area",                                                    this.lma_max_area,  6, 8, "sq.pix",
					"Maximal product of maximum half-width by half-height, ignore check if <=0");

			gd.addMessage("Correlation strength calculation (match legacy)");
			gd.addNumericField("Composite correlation strength scale",                            this.lma_str_scale,  6, 8, "",
					"Multiply LMA composite correlation strength to match older CM value");
			gd.addNumericField("Composite correlation strength offset",                           this.lma_str_offset,  6, 8, "",
					"Add to scaled composite correlation strength to match older CM value");

			gd.addMessage("Lazy eye results interpretation");
			gd.addCheckbox    ("Convert dd, nd -> dx, dx",                                        this.lma_diff_xy,
					"Convert per camera disparity correction (and orthogonal to disparity)  to X,Y correction");
			gd.addNumericField("Minimal per-tile weight for lazy eye correction",                 this.lma_diff_minw,  6, 8, "",
					"Minimal average per-tile weight to use this tile cluster for lazy eye correction");
			gd.addNumericField("Lazy eye data blur sigma",                                        this.lma_diff_sigma,  6, 8, "",
					"Lazy eye data blur sigma relative to tile cluster side");


			gd.addMessage("Debug LMA parameters");
			gd.addNumericField("LMA debug level",                                                 this.lma_debug_level,  0, 3, "",
					"Debug/verbosity level for the LMA correaltion maximum fitting");
			gd.addNumericField("LMA debug level1",                                                this.lma_debug_level1,  0, 3, "",
					"Debug/verbosity level for the new LMA correaltion maximum fitting");
			gd.addCheckbox    ("Enable LMA debug images",                                         this.lma_debug_graphic,
					"If false, no debug images generated regardless of debug levels");

			gd.addCheckbox    ("Use new correlation methods compatible with x8 camera",           this.corr_var_cam,
					"Debug feature to compare old/new methods");

			gd.addNumericField("Normalization for the CM correlation strength (old)",             this.cm_max_normalization,  6, 8, "",
					"Fraction of correlation maximum radius, being squared multiplied by maximum to have the same total mass. ~= 0.5, the lower the value, the higher strength reported by the CM");
			gd.addMessage("Cameras offsets in the disparity direction and orthogonal to disparity (debugging LMA)");
			gd.addNumericField("LMA debug offsets scale",                                         this.lma_dbg_scale,  6, 8, "",
					"Scale the following offsets by this value");
			gd.addNumericField("LMA debug offset: camera0, parallel",                   this.lma_dbg_offset[0][0],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera0, ortho",                      this.lma_dbg_offset[0][1],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera1, parallel",                   this.lma_dbg_offset[1][0],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera1, ortho",                      this.lma_dbg_offset[1][1],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera2, parallel",                   this.lma_dbg_offset[2][0],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera2, ortho",                      this.lma_dbg_offset[2][1],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera3, parallel",                   this.lma_dbg_offset[3][0],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");
			gd.addNumericField("LMA debug offset: camera3, ortho",                      this.lma_dbg_offset[3][1],  6, 8, "pix",
					"Add camera offset in the direction of disparity (to/from center)");


			//	public double  cm_max_normalization =   0.55; //

	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
			this.gpu_mode_debug =        gd.getNextBoolean();
			this.gpu_verify =            gd.getNextBoolean();
			this.corr_mode_debug=        gd.getNextBoolean();
			this.mix_corr_poly=          gd.getNextBoolean();
			this.corr_poly_only=         gd.getNextBoolean();
			this.min_poly_strength=      gd.getNextNumber();
			this.max_poly_hwidth=        gd.getNextNumber();
			this.poly_corr_scale=        gd.getNextNumber();

			this.poly_pwr=               gd.getNextNumber();
			this.poly_vasw_pwr=          gd.getNextNumber();
			this.corr_magic_scale_cm=    gd.getNextNumber();
			this.corr_magic_scale_poly=  gd.getNextNumber();

			this.ortho_height =    (int) gd.getNextNumber();
			this.ortho_eff_height=       gd.getNextNumber();
			this.ortho_nsamples =  (int) gd.getNextNumber();
			this.ortho_vasw_pwr =        gd.getNextNumber();

  			this.enhortho_width=   (int) gd.getNextNumber();
  			this.enhortho_width_aux=(int)gd.getNextNumber();
  			this.enhortho_scale=         gd.getNextNumber();
  			this.enhortho_scale_aux=     gd.getNextNumber();

  			this.ly_poly =               gd.getNextBoolean();
  			this.ly_crazy_poly=          gd.getNextNumber();
  			this.ly_poly_backup =        gd.getNextBoolean();

  			this.fo_correct =            gd.getNextBoolean();
  			this.fo_far =                gd.getNextBoolean();
  			this.fo_min_strength =       gd.getNextNumber();
  			this.fo_min_eff =            gd.getNextNumber();
  			this.fo_min_eff_ratio =      gd.getNextNumber();
  			this.fo_max_hwidth =         gd.getNextNumber();
  			this.fo_min_diff =           gd.getNextNumber();
  			this.fo_ortho =              gd.getNextBoolean();
  			this.fo_overcorrection =     gd.getNextNumber();
  			this.fo_lim_overcorr =       gd.getNextNumber();

  			this.mismatch_max_diff =     gd.getNextNumber();

  			this.corr_offset =           gd.getNextNumber();
  			this.twice_diagonal =        gd.getNextBoolean();
		    this.min_corr =              gd.getNextNumber();


  			this.dbg_pair_mask=    (int) gd.getNextNumber();
  			this.corr_strip_hight= (int) gd.getNextNumber();


  			this.cnvx_hwnd_size=   (int) gd.getNextNumber();
			this.cnvx_weight =           gd.getNextNumber();
  			this.cnvx_add3x3 =           gd.getNextBoolean();

  			this.cnvx_min_samples= (int) gd.getNextNumber();
  			this.cnvx_non_coll =         gd.getNextBoolean();
  			this.cnvx_min_pairs=   (int) gd.getNextNumber();

  			this.mcorr_multi=      (int) gd.getNextNumber();
  			this.mcorr_all =             gd.getNextBoolean();
  			this.mcorr_all_multi =       gd.getNextBoolean();
  			this.mcorr_dia =             gd.getNextBoolean();
  			this.mcorr_dia_multi =       gd.getNextBoolean();
  			this.mcorr_sq =              gd.getNextBoolean();
  			this.mcorr_sq_multi =        gd.getNextBoolean();
  			this.mcorr_neib =            gd.getNextBoolean();
  			this.mcorr_neib_multi =      gd.getNextBoolean();
  			
  			this.mcorr_hor =             gd.getNextBoolean();
  			this.mcorr_hor_multi =       gd.getNextBoolean();
  			this.mcorr_vert =            gd.getNextBoolean();
  			this.mcorr_vert_multi =      gd.getNextBoolean();
  			
  			this.mcorr_cons_all =        gd.getNextBoolean();
  			this.mcorr_cons_dia =        gd.getNextBoolean();
  			this.mcorr_cons_sq =         gd.getNextBoolean();
  			this.mcorr_cons_neib =       gd.getNextBoolean();
  			this.mcorr_cons_hor =        gd.getNextBoolean();
  			this.mcorr_cons_vert =       gd.getNextBoolean();
  			
  			this.mcorr_comb_width= (int) gd.getNextNumber();
  			this.mcorr_comb_height=(int) gd.getNextNumber();
  			this.mcorr_comb_offset=(int) gd.getNextNumber();
  			this.mcorr_comb_disp=        gd.getNextNumber();
  			this.mcorr_comb_dbg =        gd.getNextBoolean();
  			
  			this.corr_strip_notch= (int) gd.getNextNumber();
  			this.corr_notch_hwidth=      gd.getNextNumber();
  			this.corr_notch_blur=        gd.getNextNumber();

  			this.corr_wndy_size=   (int) gd.getNextNumber();

			this.corr_wndy_hwidth =      gd.getNextNumber();
			this.corr_wndy_blur =        gd.getNextNumber();

  			this.corr_wndx_size=   (int) gd.getNextNumber();

			this.corr_wndx_hwidth =      gd.getNextNumber();
			this.corr_wndx_blur =        gd.getNextNumber();

			this.pcorr_use=              gd.getNextBoolean();
			this.pcorr_use_hv=           gd.getNextBoolean();
			this.pcorr_sigma_mono =      gd.getNextNumber();
			this.pcorr_sigma =           gd.getNextNumber();
			this.pcorr_sigma_rb =        gd.getNextNumber();
			this.pcorr_fat_zero =        gd.getNextNumber();
			this.pcorr_fat_zero_mono =   gd.getNextNumber();
			this.pcorr_dbg_offsx =       gd.getNextNumber();
			this.pcorr_dbg_offsy =       gd.getNextNumber();
//LMA tab
			this.lma_disp_range =        gd.getNextNumber();
			this.lmas_gaussian=          gd.getNextBoolean();
			this.lmas_adjust_wm=         gd.getNextBoolean();
			this.lmas_adjust_wy=         gd.getNextBoolean();
			this.lmas_adjust_ag=         gd.getNextBoolean();

			this.lmas_poly_str_scale =   gd.getNextNumber();
			this.lmas_poly_str_min =     gd.getNextNumber();

			this.lmas_lambda_initial =   gd.getNextNumber();
			this.lmas_rms_diff =         gd.getNextNumber();
  			this.lmas_num_iter=    (int) gd.getNextNumber();
			this.lmas_max_rel_rms =      gd.getNextNumber();
			this.lmas_min_strength =     gd.getNextNumber();
			this.lmas_min_ac =           gd.getNextNumber();
			this.lmas_min_min_ac =           gd.getNextNumber();
			this.lmas_max_area =         gd.getNextNumber();

			this.lma_gaussian=           gd.getNextBoolean();
			this.lma_second=             gd.getNextBoolean();
			this.lma_second_gaussian=    gd.getNextBoolean();
			this.lma_adjust_wm=          gd.getNextBoolean();
			this.lma_adjust_wy=          gd.getNextBoolean();
			this.lma_adjust_wxy=         gd.getNextBoolean();
			this.lma_adjust_ly1=         gd.getNextBoolean();
			this.lma_adjust_ag=          gd.getNextBoolean();

			this.lma_wnd =               gd.getNextNumber();
			this.lma_min_wnd =           gd.getNextNumber();
			this.lma_wnd_pwr =           gd.getNextNumber();
  			this.lma_hard_marg=    (int) gd.getNextNumber();
  			this.lma_soft_marg=    (int) gd.getNextNumber();
			this.lma_sigma =             gd.getNextNumber();


			this.lma_half_width =        gd.getNextNumber();
			this.lma_cost_wy =           gd.getNextNumber();
			this.lma_cost_wxy =          gd.getNextNumber();

			this.lma_lambda_initial =    gd.getNextNumber();
			this.lma_lambda_scale_good = gd.getNextNumber();
			this.lma_lambda_scale_bad =  gd.getNextNumber();
			this.lma_lambda_max =        gd.getNextNumber();
			this.lma_rms_diff =          gd.getNextNumber();
  			this.lma_num_iter=     (int) gd.getNextNumber();

			this.lma_max_rel_rms =       gd.getNextNumber();
			this.lma_min_strength =      gd.getNextNumber();
			this.lma_min_ac =            gd.getNextNumber();
			this.lma_min_min_ac =        gd.getNextNumber();
			this.lma_max_area =          gd.getNextNumber();

			this.lma_str_scale =         gd.getNextNumber();
			this.lma_str_offset =        gd.getNextNumber();

  			this.lma_diff_xy =           gd.getNextBoolean();
  			this.lma_diff_minw=          gd.getNextNumber();
  			this.lma_diff_sigma=         gd.getNextNumber();

  			this.lma_debug_level=  (int) gd.getNextNumber();
  			this.lma_debug_level1= (int) gd.getNextNumber();
  			this.lma_debug_graphic =     gd.getNextBoolean();
  			this.corr_var_cam =          gd.getNextBoolean();
  			this.cm_max_normalization=   gd.getNextNumber();
  			this.lma_dbg_scale=          gd.getNextNumber();

  			for (int i = 0; i < 4; i++) for (int j=0; j < 2; j++) {
  				this.lma_dbg_offset[i][j]=   gd.getNextNumber();
  			}
	}


	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"gpu_mode_debug",       this.gpu_mode_debug+"");
		properties.setProperty(prefix+"gpu_verify",           this.gpu_verify+"");
		properties.setProperty(prefix+"corr_mode_debug",      this.corr_mode_debug+"");
		properties.setProperty(prefix+"mix_corr_poly",        this.mix_corr_poly+"");
		properties.setProperty(prefix+"corr_poly_only",       this.corr_poly_only+"");
		properties.setProperty(prefix+"min_poly_strength",    this.min_poly_strength+"");
		properties.setProperty(prefix+"max_poly_hwidth",      this.max_poly_hwidth+"");
		properties.setProperty(prefix+"poly_corr_scale",      this.poly_corr_scale+"");

		properties.setProperty(prefix+"poly_pwr",             this.poly_pwr+"");
		properties.setProperty(prefix+"poly_vasw_pwr",        this.poly_vasw_pwr+"");
		properties.setProperty(prefix+"corr_magic_scale_cm",  this.corr_magic_scale_cm+"");
		properties.setProperty(prefix+"corr_magic_scale_poly",this.corr_magic_scale_poly+"");

		properties.setProperty(prefix+"ortho_height",         this.ortho_height+"");
		properties.setProperty(prefix+"ortho_eff_height",     this.ortho_eff_height+"");
		properties.setProperty(prefix+"ortho_nsamples",       this.ortho_nsamples+"");
		properties.setProperty(prefix+"ortho_vasw_pwr",       this.ortho_vasw_pwr+"");

		properties.setProperty(prefix+"enhortho_width",       this.enhortho_width +"");
		properties.setProperty(prefix+"enhortho_width_aux",   this.enhortho_width_aux +"");
		properties.setProperty(prefix+"enhortho_scale",       this.enhortho_scale +"");
		properties.setProperty(prefix+"enhortho_scale_aux",   this.enhortho_scale_aux +"");

		properties.setProperty(prefix+"corr_offset",          this.corr_offset +"");
		properties.setProperty(prefix+"twice_diagonal",       this.twice_diagonal +"");
		properties.setProperty(prefix+"min_corr",             this.min_corr +"");

		properties.setProperty(prefix+"ly_poly",              this.ly_poly +"");
		properties.setProperty(prefix+"ly_crazy_poly",        this.ly_crazy_poly +"");
		properties.setProperty(prefix+"ly_poly_backup",       this.ly_poly_backup +"");

		properties.setProperty(prefix+"fo_correct",           this.fo_correct +"");
		properties.setProperty(prefix+"fo_far",               this.fo_far +"");
		properties.setProperty(prefix+"fo_min_strength",      this.fo_min_strength +"");
		properties.setProperty(prefix+"fo_min_eff",           this.fo_min_eff +"");
		properties.setProperty(prefix+"fo_min_eff_ratio",     this.fo_min_eff_ratio +"");
		properties.setProperty(prefix+"fo_max_hwidth",        this.fo_max_hwidth +"");
		properties.setProperty(prefix+"fo_min_diff",          this.fo_min_diff +"");
		properties.setProperty(prefix+"fo_ortho",             this.fo_ortho +"");
		properties.setProperty(prefix+"fo_overcorrection",    this.fo_overcorrection +"");
		properties.setProperty(prefix+"fo_lim_overcorr",      this.fo_lim_overcorr +"");

		properties.setProperty(prefix+"mismatch_max_diff",    this.mismatch_max_diff +"");

		properties.setProperty(prefix+"dbg_pair_mask",        this.dbg_pair_mask +"");
		properties.setProperty(prefix+"corr_strip_hight",     this.corr_strip_hight +"");

		properties.setProperty(prefix+"cnvx_hwnd_size",       this.cnvx_hwnd_size +"");
		properties.setProperty(prefix+"cnvx_weight",          this.cnvx_weight +"");
		properties.setProperty(prefix+"cnvx_add3x3",          this.cnvx_add3x3 +"");

		properties.setProperty(prefix+"cnvx_min_samples",     this.cnvx_min_samples +"");
		properties.setProperty(prefix+"cnvx_non_coll",        this.cnvx_non_coll +"");
		properties.setProperty(prefix+"cnvx_min_pairs",       this.cnvx_min_pairs +"");

		properties.setProperty(prefix+"mcorr_multi",          this.mcorr_multi +"");
		properties.setProperty(prefix+"mcorr_all",            this.mcorr_all +"");
		properties.setProperty(prefix+"mcorr_all_multi",      this.mcorr_all_multi +"");
		properties.setProperty(prefix+"mcorr_dia",            this.mcorr_dia +"");
		properties.setProperty(prefix+"mcorr_dia_multi",      this.mcorr_dia_multi +"");
		properties.setProperty(prefix+"mcorr_sq",             this.mcorr_sq +"");
		properties.setProperty(prefix+"mcorr_sq_multi",       this.mcorr_sq_multi +"");
		properties.setProperty(prefix+"mcorr_neib",           this.mcorr_neib +"");
		properties.setProperty(prefix+"mcorr_neib_multi",     this.mcorr_neib_multi +"");

		properties.setProperty(prefix+"mcorr_hor",            this.mcorr_hor +"");
		properties.setProperty(prefix+"mcorr_hor_multi",      this.mcorr_hor_multi +"");
		properties.setProperty(prefix+"mcorr_vert",           this.mcorr_vert +"");
		properties.setProperty(prefix+"mcorr_vert_multi",     this.mcorr_vert_multi +"");
		
		properties.setProperty(prefix+"mcorr_cons_all",       this.mcorr_cons_all +"");
		properties.setProperty(prefix+"mcorr_cons_dia",       this.mcorr_cons_dia +"");
		properties.setProperty(prefix+"mcorr_cons_sq",        this.mcorr_cons_sq +"");
		properties.setProperty(prefix+"mcorr_cons_neib",      this.mcorr_cons_neib +"");
		properties.setProperty(prefix+"mcorr_cons_hor",       this.mcorr_cons_hor +"");
		properties.setProperty(prefix+"mcorr_cons_vert",      this.mcorr_cons_vert +"");
		
		properties.setProperty(prefix+"mcorr_comb_width",     this.mcorr_comb_width +"");
		properties.setProperty(prefix+"mcorr_comb_height",    this.mcorr_comb_height +"");
		properties.setProperty(prefix+"mcorr_comb_offset",    this.mcorr_comb_offset +"");
		properties.setProperty(prefix+"mcorr_comb_disp",      this.mcorr_comb_disp +"");
		properties.setProperty(prefix+"mcorr_comb_dbg",       this.mcorr_comb_dbg +"");
		
		properties.setProperty(prefix+"corr_strip_notch",     this.corr_strip_notch +"");
		properties.setProperty(prefix+"corr_notch_hwidth",    this.corr_notch_hwidth +"");
		properties.setProperty(prefix+"corr_notch_blur",      this.corr_notch_blur +"");

		properties.setProperty(prefix+"corr_wndy_size",       this.corr_wndy_size +"");
		properties.setProperty(prefix+"corr_wndy_hwidth",     this.corr_wndy_hwidth +"");
		properties.setProperty(prefix+"corr_wndy_blur",       this.corr_wndy_blur +"");

		properties.setProperty(prefix+"corr_wndx_size",       this.corr_wndx_size +"");

		properties.setProperty(prefix+"corr_wndx_hwidth",     this.corr_wndx_hwidth +"");
		properties.setProperty(prefix+"corr_wndx_blur",       this.corr_wndx_blur +"");

		properties.setProperty(prefix+"pcorr_use",            this.pcorr_use +"");
		properties.setProperty(prefix+"pcorr_use_hv",         this.pcorr_use_hv +"");
		properties.setProperty(prefix+"pcorr_sigma_mono",     this.pcorr_sigma_mono +"");
		properties.setProperty(prefix+"pcorr_sigma",          this.pcorr_sigma +"");
		properties.setProperty(prefix+"pcorr_sigma_rb",       this.pcorr_sigma_rb +"");
		properties.setProperty(prefix+"pcorr_fat_zero",       this.pcorr_fat_zero +"");
		properties.setProperty(prefix+"pcorr_fat_zero_mono",  this.pcorr_fat_zero_mono +"");
		properties.setProperty(prefix+"pcorr_dbg_offsx",      this.pcorr_dbg_offsx +"");
		properties.setProperty(prefix+"pcorr_dbg_offsy",      this.pcorr_dbg_offsy +"");

		properties.setProperty(prefix+"lma_disp_range",       this.lma_disp_range +"");
		properties.setProperty(prefix+"lmas_gaussian",        this.lmas_gaussian +"");
		properties.setProperty(prefix+"lmas_adjust_wm",       this.lmas_adjust_wm +"");
		properties.setProperty(prefix+"lmas_adjust_wy",       this.lmas_adjust_wy +"");
		properties.setProperty(prefix+"lmas_adjust_ag",       this.lmas_adjust_ag +"");

		properties.setProperty(prefix+"lmas_poly_str_scale",  this.lmas_poly_str_scale +"");
		properties.setProperty(prefix+"lmas_poly_str_min",    this.lmas_poly_str_min +"");

		properties.setProperty(prefix+"lmas_lambda_initial",  this.lmas_lambda_initial +"");
		properties.setProperty(prefix+"lmas_rms_diff",        this.lmas_rms_diff +"");
		properties.setProperty(prefix+"lmas_num_iter",        this.lmas_num_iter +"");
		properties.setProperty(prefix+"lmas_max_rel_rms",     this.lmas_max_rel_rms +"");
		properties.setProperty(prefix+"lmas_min_strength",    this.lmas_min_strength +"");
		properties.setProperty(prefix+"lmas_min_ac",          this.lmas_min_ac +"");
		properties.setProperty(prefix+"lmas_min_min_ac",      this.lmas_min_min_ac +"");
		properties.setProperty(prefix+"lmas_max_area",        this.lmas_max_area +"");

		properties.setProperty(prefix+"lma_gaussian",         this.lma_gaussian +"");
		properties.setProperty(prefix+"lma_second",           this.lma_second +"");
		properties.setProperty(prefix+"lma_second_gaussian",  this.lma_second_gaussian +"");
		properties.setProperty(prefix+"lma_adjust_wm",        this.lma_adjust_wm +"");
		properties.setProperty(prefix+"lma_adjust_wy",        this.lma_adjust_wy +"");
		properties.setProperty(prefix+"lma_adjust_wxy",       this.lma_adjust_wxy +"");
		properties.setProperty(prefix+"lma_adjust_ly1",       this.lma_adjust_ly1 +"");

		properties.setProperty(prefix+"lma_adjust_ag",        this.lma_adjust_ag +"");

		properties.setProperty(prefix+"lma_wnd",              this.lma_wnd +"");
		properties.setProperty(prefix+"lma_min_wnd",          this.lma_min_wnd +"");
		properties.setProperty(prefix+"lma_wnd_pwr",          this.lma_wnd_pwr +"");
		properties.setProperty(prefix+"lma_hard_marg",        this.lma_hard_marg +"");
		properties.setProperty(prefix+"lma_soft_marg",        this.lma_soft_marg +"");
		properties.setProperty(prefix+"lma_sigma",            this.lma_sigma +"");


		properties.setProperty(prefix+"lma_half_width",       this.lma_half_width +"");
		properties.setProperty(prefix+"lma_cost_wy",          this.lma_cost_wy +"");
		properties.setProperty(prefix+"lma_cost_wxy",         this.lma_cost_wxy +"");

		properties.setProperty(prefix+"lma_lambda_initial",   this.lma_lambda_initial +"");
		properties.setProperty(prefix+"lma_lambda_scale_good",this.lma_lambda_scale_good +"");
		properties.setProperty(prefix+"lma_lambda_scale_bad", this.lma_lambda_scale_bad +"");
		properties.setProperty(prefix+"lma_lambda_max",       this.lma_lambda_max +"");
		properties.setProperty(prefix+"lma_rms_diff",         this.lma_rms_diff +"");
		properties.setProperty(prefix+"lma_num_iter",         this.lma_num_iter +"");

		properties.setProperty(prefix+"lma_max_rel_rms",      this.lma_max_rel_rms +"");
		properties.setProperty(prefix+"lma_min_strength",     this.lma_min_strength +"");
		properties.setProperty(prefix+"lma_min_ac",           this.lma_min_ac +"");
		properties.setProperty(prefix+"lma_min_min_ac",       this.lma_min_min_ac +"");
		properties.setProperty(prefix+"lma_max_area",         this.lma_max_area +"");
		properties.setProperty(prefix+"lma_str_scale",        this.lma_str_scale +"");
		properties.setProperty(prefix+"lma_str_offset",       this.lma_str_offset +"");

		properties.setProperty(prefix+"lma_diff_xy",          this.lma_diff_xy +"");
		properties.setProperty(prefix+"lma_diff_minw",        this.lma_diff_minw +"");
		properties.setProperty(prefix+"lma_diff_sigma",       this.lma_diff_sigma +"");

		properties.setProperty(prefix+"lma_debug_level",      this.lma_debug_level +"");
		properties.setProperty(prefix+"lma_debug_level1",     this.lma_debug_level1 +"");
		properties.setProperty(prefix+"lma_debug_graphic",    this.lma_debug_graphic +"");

		properties.setProperty(prefix+"corr_var_cam",         this.corr_var_cam +"");

		properties.setProperty(prefix+"cm_max_normalization", this.cm_max_normalization +"");

		properties.setProperty(prefix+"lma_dbg_scale",        this.lma_dbg_scale +"");
		for (int i = 0; i < 4; i++) for (int j=0; j < 2; j++) {
			properties.setProperty(prefix+"lma_dbg_offset_"+i+"_"+j, this.lma_dbg_offset[i][j] +"");
	    }


	}

	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"gpu_mode_debug")!=null)        this.gpu_mode_debug=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_mode_debug"));
		if (properties.getProperty(prefix+"gpu_verify")!=null)            this.gpu_verify=Boolean.parseBoolean(properties.getProperty(prefix+"gpu_verify"));
		if (properties.getProperty(prefix+"corr_mode_debug")!=null)       this.corr_mode_debug=Boolean.parseBoolean(properties.getProperty(prefix+"corr_mode_debug"));
		if (properties.getProperty(prefix+"mix_corr_poly")!=null)         this.mix_corr_poly=Boolean.parseBoolean(properties.getProperty(prefix+"mix_corr_poly"));
		if (properties.getProperty(prefix+"corr_poly_only")!=null)        this.corr_poly_only=Boolean.parseBoolean(properties.getProperty(prefix+"corr_poly_only"));
		if (properties.getProperty(prefix+"min_poly_strength")!=null)     this.min_poly_strength=Double.parseDouble(properties.getProperty(prefix+"min_poly_strength"));
		if (properties.getProperty(prefix+"max_poly_hwidth")!=null)       this.max_poly_hwidth=Double.parseDouble(properties.getProperty(prefix+"max_poly_hwidth"));
		if (properties.getProperty(prefix+"poly_corr_scale")!=null)       this.poly_corr_scale=Double.parseDouble(properties.getProperty(prefix+"poly_corr_scale"));

		if (properties.getProperty(prefix+"poly_pwr")!=null)              this.poly_pwr=Double.parseDouble(properties.getProperty(prefix+"poly_pwr"));
		if (properties.getProperty(prefix+"poly_vasw_pwr")!=null)         this.poly_vasw_pwr=Double.parseDouble(properties.getProperty(prefix+"poly_vasw_pwr"));
		if (properties.getProperty(prefix+"corr_magic_scale_cm")!=null)   this.corr_magic_scale_cm=Double.parseDouble(properties.getProperty(prefix+"corr_magic_scale_cm"));
		if (properties.getProperty(prefix+"corr_magic_scale_poly")!=null) this.corr_magic_scale_poly=Double.parseDouble(properties.getProperty(prefix+"corr_magic_scale_poly"));

		if (properties.getProperty(prefix+"ortho_height")!=null)          this.ortho_height=Integer.parseInt(properties.getProperty(prefix+"ortho_height"));
		if (properties.getProperty(prefix+"ortho_eff_height")!=null)      this.ortho_eff_height=Double.parseDouble(properties.getProperty(prefix+"ortho_eff_height"));
		if (properties.getProperty(prefix+"ortho_nsamples")!=null)        this.ortho_nsamples=Integer.parseInt(properties.getProperty(prefix+"ortho_nsamples"));
		if (properties.getProperty(prefix+"ortho_vasw_pwr")!=null)        this.ortho_vasw_pwr=Double.parseDouble(properties.getProperty(prefix+"ortho_vasw_pwr"));

		if (properties.getProperty(prefix+"enhortho_width")!=null)        this.enhortho_width=Integer.parseInt(properties.getProperty(prefix+"enhortho_width"));
		if (properties.getProperty(prefix+"enhortho_width_aux")!=null)    this.enhortho_width_aux=Integer.parseInt(properties.getProperty(prefix+"enhortho_width_aux"));
		if (properties.getProperty(prefix+"enhortho_scale")!=null)        this.enhortho_scale=Double.parseDouble(properties.getProperty(prefix+"enhortho_scale"));
		if (properties.getProperty(prefix+"enhortho_scale_aux")!=null)    this.enhortho_scale_aux=Double.parseDouble(properties.getProperty(prefix+"enhortho_scale_aux"));

		if (properties.getProperty(prefix+"fo_correct")!=null)            this.fo_correct=Boolean.parseBoolean(properties.getProperty(prefix+"fo_correct"));

		if (properties.getProperty(prefix+"ly_poly")!=null)               this.ly_poly=Boolean.parseBoolean(properties.getProperty(prefix+"ly_poly"));
		if (properties.getProperty(prefix+"ly_crazy_poly")!=null)         this.ly_crazy_poly=Double.parseDouble(properties.getProperty(prefix+"ly_crazy_poly"));
		if (properties.getProperty(prefix+"ly_poly_backup")!=null)        this.ly_poly_backup=Boolean.parseBoolean(properties.getProperty(prefix+"ly_poly_backup"));

		if (properties.getProperty(prefix+"fo_far")!=null)                this.fo_far=Boolean.parseBoolean(properties.getProperty(prefix+"fo_far"));
		if (properties.getProperty(prefix+"fo_min_strength")!=null)       this.fo_min_strength=Double.parseDouble(properties.getProperty(prefix+"fo_min_strength"));
		if (properties.getProperty(prefix+"fo_min_eff")!=null)            this.fo_min_eff=Double.parseDouble(properties.getProperty(prefix+"fo_min_eff"));
		if (properties.getProperty(prefix+"fo_min_eff_ratio")!=null)      this.fo_min_eff_ratio=Double.parseDouble(properties.getProperty(prefix+"fo_min_eff_ratio"));
		if (properties.getProperty(prefix+"fo_max_hwidth")!=null)         this.fo_max_hwidth=Double.parseDouble(properties.getProperty(prefix+"fo_max_hwidth"));
		if (properties.getProperty(prefix+"fo_min_diff")!=null)           this.fo_min_diff=Double.parseDouble(properties.getProperty(prefix+"fo_min_diff"));
		if (properties.getProperty(prefix+"fo_ortho")!=null)              this.fo_ortho=Boolean.parseBoolean(properties.getProperty(prefix+"fo_ortho"));
		if (properties.getProperty(prefix+"fo_overcorrection")!=null)     this.fo_overcorrection=Double.parseDouble(properties.getProperty(prefix+"fo_overcorrection"));
		if (properties.getProperty(prefix+"fo_lim_overcorr")!=null)       this.fo_lim_overcorr=Double.parseDouble(properties.getProperty(prefix+"fo_lim_overcorr"));

		if (properties.getProperty(prefix+"mismatch_max_diff")!=null)     this.mismatch_max_diff=Double.parseDouble(properties.getProperty(prefix+"mismatch_max_diff"));

		if (properties.getProperty(prefix+"corr_offset")!=null)           this.corr_offset=Double.parseDouble(properties.getProperty(prefix+"corr_offset"));
		if (properties.getProperty(prefix+"twice_diagonal")!=null)        this.twice_diagonal=Boolean.parseBoolean(properties.getProperty(prefix+"twice_diagonal"));
		if (properties.getProperty(prefix+"min_corr")!=null)              this.min_corr=Double.parseDouble(properties.getProperty(prefix+"min_corr"));


		if (properties.getProperty(prefix+"dbg_pair_mask")!=null)        this.dbg_pair_mask=Integer.parseInt(properties.getProperty(prefix+"dbg_pair_mask"));
		if (properties.getProperty(prefix+"corr_strip_hight")!=null)     this.corr_strip_hight=Integer.parseInt(properties.getProperty(prefix+"corr_strip_hight"));

		if (properties.getProperty(prefix+"cnvx_hwnd_size")!=null)       this.cnvx_hwnd_size=Integer.parseInt(properties.getProperty(prefix+"cnvx_hwnd_size"));
		if (properties.getProperty(prefix+"cnvx_weight")!=null)          this.cnvx_weight=Double.parseDouble(properties.getProperty(prefix+"cnvx_weight"));
		if (properties.getProperty(prefix+"cnvx_add3x3")!=null)          this.cnvx_add3x3=Boolean.parseBoolean(properties.getProperty(prefix+"cnvx_add3x3"));

		if (properties.getProperty(prefix+"cnvx_min_samples")!=null)     this.cnvx_min_samples=Integer.parseInt(properties.getProperty(prefix+"cnvx_min_samples"));
		if (properties.getProperty(prefix+"cnvx_non_coll")!=null)        this.cnvx_non_coll=Boolean.parseBoolean(properties.getProperty(prefix+"cnvx_non_coll"));
		if (properties.getProperty(prefix+"cnvx_min_pairs")!=null)       this.cnvx_min_pairs=Integer.parseInt(properties.getProperty(prefix+"cnvx_min_pairs"));
		
		if (properties.getProperty(prefix+"mcorr_multi")!=null)          this.mcorr_multi=Integer.parseInt(properties.getProperty(prefix+"mcorr_multi"));
		if (properties.getProperty(prefix+"mcorr_all")!=null)            this.mcorr_all=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_all"));
		if (properties.getProperty(prefix+"mcorr_all_multi")!=null)      this.mcorr_all_multi=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_all_multi"));
		if (properties.getProperty(prefix+"mcorr_dia")!=null)            this.mcorr_dia=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_dia"));
		if (properties.getProperty(prefix+"mcorr_dia_multi")!=null)      this.mcorr_dia_multi=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_dia_multi"));
		if (properties.getProperty(prefix+"mcorr_sq")!=null)             this.mcorr_sq=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_sq"));
		if (properties.getProperty(prefix+"mcorr_sq_multi")!=null)       this.mcorr_sq_multi=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_sq_multi"));
		if (properties.getProperty(prefix+"mcorr_neib")!=null)           this.mcorr_neib=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_neib"));
		if (properties.getProperty(prefix+"mcorr_neib_multi")!=null)     this.mcorr_neib_multi=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_neib_multi"));
		if (properties.getProperty(prefix+"mcorr_hor")!=null)            this.mcorr_hor=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_hor"));
		if (properties.getProperty(prefix+"mcorr_hor_multi")!=null)      this.mcorr_hor_multi=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_hor_multi"));
		if (properties.getProperty(prefix+"mcorr_vert")!=null)           this.mcorr_vert=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_vert"));
		if (properties.getProperty(prefix+"mcorr_vert_multi")!=null)     this.mcorr_vert_multi=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_vert_multi"));
		
		if (properties.getProperty(prefix+"mcorr_cons_all")!=null)       this.mcorr_cons_all=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_cons_all"));
		if (properties.getProperty(prefix+"mcorr_cons_dia")!=null)       this.mcorr_cons_dia=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_cons_dia"));
		if (properties.getProperty(prefix+"mcorr_cons_sq")!=null)        this.mcorr_cons_sq=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_cons_sq"));
		if (properties.getProperty(prefix+"mcorr_cons_neib")!=null)      this.mcorr_cons_neib=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_cons_neib"));
		if (properties.getProperty(prefix+"mcorr_cons_hor")!=null)       this.mcorr_cons_hor=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_cons_hor"));
		if (properties.getProperty(prefix+"mcorr_cons_vert")!=null)      this.mcorr_cons_vert=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_cons_vert"));
		
		if (properties.getProperty(prefix+"mcorr_comb_width")!=null)     this.mcorr_comb_width=Integer.parseInt(properties.getProperty(prefix+"mcorr_comb_width"));
		if (properties.getProperty(prefix+"mcorr_comb_height")!=null)    this.mcorr_comb_height=Integer.parseInt(properties.getProperty(prefix+"mcorr_comb_height"));
		if (properties.getProperty(prefix+"mcorr_comb_offset")!=null)    this.mcorr_comb_offset=Integer.parseInt(properties.getProperty(prefix+"mcorr_comb_offset"));
		if (properties.getProperty(prefix+"mcorr_comb_disp")!=null)      this.mcorr_comb_disp=Double.parseDouble(properties.getProperty(prefix+"mcorr_comb_disp"));
		if (properties.getProperty(prefix+"mcorr_comb_dbg")!=null)       this.mcorr_comb_dbg=Boolean.parseBoolean(properties.getProperty(prefix+"mcorr_comb_dbg"));
		
		if (properties.getProperty(prefix+"corr_strip_notch")!=null)     this.corr_strip_notch=Integer.parseInt(properties.getProperty(prefix+"corr_strip_notch"));
		if (properties.getProperty(prefix+"corr_notch_hwidth")!=null)    this.corr_notch_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_notch_hwidth"));
		if (properties.getProperty(prefix+"corr_notch_blur")!=null)      this.corr_notch_blur=Double.parseDouble(properties.getProperty(prefix+"corr_notch_blur"));

		if (properties.getProperty(prefix+"corr_wndy_size")!=null)       this.corr_wndy_size=Integer.parseInt(properties.getProperty(prefix+"corr_wndy_size"));

		if (properties.getProperty(prefix+"corr_wndy_hwidth")!=null)     this.corr_wndy_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_wndy_hwidth"));
		if (properties.getProperty(prefix+"corr_wndy_blur")!=null)       this.corr_wndy_blur=Double.parseDouble(properties.getProperty(prefix+"corr_wndy_blur"));

		if (properties.getProperty(prefix+"corr_wndx_size")!=null)       this.corr_wndx_size=Integer.parseInt(properties.getProperty(prefix+"corr_wndx_size"));

		if (properties.getProperty(prefix+"corr_wndx_hwidth")!=null)     this.corr_wndx_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_wndx_hwidth"));
		if (properties.getProperty(prefix+"corr_wndx_blur")!=null)       this.corr_wndx_blur=Double.parseDouble(properties.getProperty(prefix+"corr_wndx_blur"));

		if (properties.getProperty(prefix+"pcorr_use")!=null)            this.pcorr_use=Boolean.parseBoolean(properties.getProperty(prefix+"pcorr_use"));
		if (properties.getProperty(prefix+"pcorr_use_hv")!=null)         this.pcorr_use_hv=Boolean.parseBoolean(properties.getProperty(prefix+"pcorr_use_hv"));
		if (properties.getProperty(prefix+"pcorr_sigma_mono")!=null)     this.pcorr_sigma_mono=Double.parseDouble(properties.getProperty(prefix+"pcorr_sigma_mono"));
		if (properties.getProperty(prefix+"pcorr_sigma")!=null)          this.pcorr_sigma=Double.parseDouble(properties.getProperty(prefix+"pcorr_sigma"));
		if (properties.getProperty(prefix+"pcorr_sigma_rb")!=null)       this.pcorr_sigma_rb=Double.parseDouble(properties.getProperty(prefix+"pcorr_sigma_rb"));
		if (properties.getProperty(prefix+"pcorr_fat_zero")!=null)       this.pcorr_fat_zero=Double.parseDouble(properties.getProperty(prefix+"pcorr_fat_zero"));
		if (properties.getProperty(prefix+"pcorr_fat_zero_mono")!=null)  this.pcorr_fat_zero_mono=Double.parseDouble(properties.getProperty(prefix+"pcorr_fat_zero_mono"));
		if (properties.getProperty(prefix+"pcorr_dbg_offsx")!=null)      this.pcorr_dbg_offsx=Double.parseDouble(properties.getProperty(prefix+"pcorr_dbg_offsx"));
		if (properties.getProperty(prefix+"pcorr_dbg_offsy")!=null)      this.pcorr_dbg_offsy=Double.parseDouble(properties.getProperty(prefix+"pcorr_dbg_offsy"));

		if (properties.getProperty(prefix+"lma_disp_range")!=null)       this.lma_disp_range=Double.parseDouble(properties.getProperty(prefix+"lma_disp_range"));
		if (properties.getProperty(prefix+"lmas_gaussian")!=null)        this.lmas_gaussian=Boolean.parseBoolean(properties.getProperty(prefix+"lmas_gaussian"));
		if (properties.getProperty(prefix+"lmas_adjust_wm")!=null)       this.lmas_adjust_wm=Boolean.parseBoolean(properties.getProperty(prefix+"lmas_adjust_wm"));
		if (properties.getProperty(prefix+"lmas_adjust_wy")!=null)       this.lmas_adjust_wy=Boolean.parseBoolean(properties.getProperty(prefix+"lmas_adjust_wy"));
		if (properties.getProperty(prefix+"lmas_adjust_ag")!=null)       this.lmas_adjust_ag=Boolean.parseBoolean(properties.getProperty(prefix+"lmas_adjust_ag"));

		if (properties.getProperty(prefix+"lmas_poly_str_scale")!=null)  this.lmas_poly_str_scale=Double.parseDouble(properties.getProperty(prefix+"lmas_poly_str_scale"));
		if (properties.getProperty(prefix+"lmas_poly_str_min")!=null)    this.lmas_poly_str_min=Double.parseDouble(properties.getProperty(prefix+"lmas_poly_str_min"));

		if (properties.getProperty(prefix+"lmas_lambda_initial")!=null)  this.lmas_lambda_initial=Double.parseDouble(properties.getProperty(prefix+"lmas_lambda_initial"));
		if (properties.getProperty(prefix+"lmas_rms_diff")!=null)        this.lmas_rms_diff=Double.parseDouble(properties.getProperty(prefix+"lmas_rms_diff"));
		if (properties.getProperty(prefix+"lmas_num_iter")!=null)        this.lmas_num_iter=Integer.parseInt(properties.getProperty(prefix+"lmas_num_iter"));
		if (properties.getProperty(prefix+"lmas_max_rel_rms")!=null)     this.lmas_max_rel_rms=Double.parseDouble(properties.getProperty(prefix+"lmas_max_rel_rms"));
		if (properties.getProperty(prefix+"lmas_min_strength")!=null)    this.lmas_min_strength=Double.parseDouble(properties.getProperty(prefix+"lmas_min_strength"));
		if (properties.getProperty(prefix+"lmas_min_ac")!=null)          this.lmas_min_ac=Double.parseDouble(properties.getProperty(prefix+"lmas_min_ac"));
		if (properties.getProperty(prefix+"lmas_min_min_ac")!=null)      this.lmas_min_min_ac=Double.parseDouble(properties.getProperty(prefix+"lmas_min_min_ac"));
		if (properties.getProperty(prefix+"lmas_max_area")!=null)        this.lmas_max_area=Double.parseDouble(properties.getProperty(prefix+"lmas_max_area"));

		if (properties.getProperty(prefix+"lma_gaussian")!=null)         this.lma_gaussian=Boolean.parseBoolean(properties.getProperty(prefix+"lma_gaussian"));
		if (properties.getProperty(prefix+"lma_second")!=null)           this.lma_second=Boolean.parseBoolean(properties.getProperty(prefix+"lma_second"));
		if (properties.getProperty(prefix+"lma_second_gaussian")!=null)  this.lma_second_gaussian=Boolean.parseBoolean(properties.getProperty(prefix+"lma_second_gaussian"));
		if (properties.getProperty(prefix+"lma_adjust_wm")!=null)        this.lma_adjust_wm=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_wm"));
		if (properties.getProperty(prefix+"lma_adjust_wy")!=null)        this.lma_adjust_wy=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_wy"));
		if (properties.getProperty(prefix+"lma_adjust_wxy")!=null)       this.lma_adjust_wxy=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_wxy"));
		if (properties.getProperty(prefix+"lma_adjust_ly1")!=null)       this.lma_adjust_ly1=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_ly1"));
		if (properties.getProperty(prefix+"lma_adjust_ag")!=null)        this.lma_adjust_ag=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_ag"));

		if (properties.getProperty(prefix+"lma_wnd")!=null)              this.lma_wnd=Double.parseDouble(properties.getProperty(prefix+"lma_wnd"));
		if (properties.getProperty(prefix+"lma_min_wnd")!=null)          this.lma_min_wnd=Double.parseDouble(properties.getProperty(prefix+"lma_min_wnd"));
		if (properties.getProperty(prefix+"lma_wnd_pwr")!=null)          this.lma_wnd_pwr=Double.parseDouble(properties.getProperty(prefix+"lma_wnd_pwr"));
		if (properties.getProperty(prefix+"lma_hard_marg")!=null)        this.lma_hard_marg=Integer.parseInt(properties.getProperty(prefix+"lma_hard_marg"));
		if (properties.getProperty(prefix+"lma_soft_marg")!=null)        this.lma_soft_marg=Integer.parseInt(properties.getProperty(prefix+"lma_soft_marg"));
		if (properties.getProperty(prefix+"lma_sigma")!=null)            this.lma_sigma=Double.parseDouble(properties.getProperty(prefix+"lma_sigma"));

		if (properties.getProperty(prefix+"lma_half_width")!=null)       this.lma_half_width=Double.parseDouble(properties.getProperty(prefix+"lma_half_width"));
		if (properties.getProperty(prefix+"lma_cost_wy")!=null)          this.lma_cost_wy=Double.parseDouble(properties.getProperty(prefix+"lma_cost_wy"));
		if (properties.getProperty(prefix+"lma_cost_wxy")!=null)         this.lma_cost_wxy=Double.parseDouble(properties.getProperty(prefix+"lma_cost_wxy"));

		if (properties.getProperty(prefix+"lma_lambda_initial")!=null)   this.lma_lambda_initial=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_initial"));
		if (properties.getProperty(prefix+"lma_lambda_scale_good")!=null)this.lma_lambda_scale_good=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_scale_good"));
		if (properties.getProperty(prefix+"lma_lambda_scale_bad")!=null) this.lma_lambda_scale_bad=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_scale_bad"));
		if (properties.getProperty(prefix+"lma_lambda_max")!=null)       this.lma_lambda_max=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_max"));
		if (properties.getProperty(prefix+"lma_rms_diff")!=null)         this.lma_rms_diff=Double.parseDouble(properties.getProperty(prefix+"lma_rms_diff"));
		if (properties.getProperty(prefix+"lma_num_iter")!=null)         this.lma_num_iter=Integer.parseInt(properties.getProperty(prefix+"lma_num_iter"));

		if (properties.getProperty(prefix+"lma_max_rel_rms")!=null)      this.lma_max_rel_rms=Double.parseDouble(properties.getProperty(prefix+"lma_max_rel_rms"));
		if (properties.getProperty(prefix+"lma_min_strength")!=null)     this.lma_min_strength=Double.parseDouble(properties.getProperty(prefix+"lma_min_strength"));
		if (properties.getProperty(prefix+"lma_min_ac")!=null)           this.lma_min_ac=Double.parseDouble(properties.getProperty(prefix+"lma_min_ac"));
		if (properties.getProperty(prefix+"lma_min_min_ac")!=null)       this.lma_min_min_ac=Double.parseDouble(properties.getProperty(prefix+"lma_min_min_ac"));
		if (properties.getProperty(prefix+"lma_max_area")!=null)         this.lma_max_area=Double.parseDouble(properties.getProperty(prefix+"lma_max_area"));
		if (properties.getProperty(prefix+"lma_str_scale")!=null)        this.lma_str_scale=Double.parseDouble(properties.getProperty(prefix+"lma_str_scale"));
		if (properties.getProperty(prefix+"lma_str_offset")!=null)       this.lma_str_offset=Double.parseDouble(properties.getProperty(prefix+"lma_str_offset"));

		if (properties.getProperty(prefix+"lma_diff_xy")!=null)         this.lma_diff_xy=Boolean.parseBoolean(properties.getProperty(prefix+"lma_diff_xy"));
		if (properties.getProperty(prefix+"lma_diff_minw")!=null)       this.lma_diff_minw=Double.parseDouble(properties.getProperty(prefix+"lma_diff_minw"));
		if (properties.getProperty(prefix+"lma_diff_sigma")!=null)      this.lma_diff_sigma=Double.parseDouble(properties.getProperty(prefix+"lma_diff_sigma"));

		if (properties.getProperty(prefix+"lma_debug_level")!=null)      this.lma_debug_level=Integer.parseInt(properties.getProperty(prefix+"lma_debug_level"));
		if (properties.getProperty(prefix+"lma_debug_level1")!=null)     this.lma_debug_level1=Integer.parseInt(properties.getProperty(prefix+"lma_debug_level1"));
		if (properties.getProperty(prefix+"lma_debug_graphic")!=null)    this.lma_debug_graphic=Boolean.parseBoolean(properties.getProperty(prefix+"lma_debug_graphic"));

		if (properties.getProperty(prefix+"corr_var_cam")!=null)         this.corr_var_cam=Boolean.parseBoolean(properties.getProperty(prefix+"corr_var_cam"));

		if (properties.getProperty(prefix+"cm_max_normalization")!=null) this.cm_max_normalization=Double.parseDouble(properties.getProperty(prefix+"cm_max_normalization"));
		if (properties.getProperty(prefix+"lma_dbg_scale")!=null)        this.lma_dbg_scale=Double.parseDouble(properties.getProperty(prefix+"lma_dbg_scale"));
		for (int i = 0; i < 4; i++) for (int j=0; j < 2; j++) {
			if (properties.getProperty(prefix+"lma_dbg_offset_"+i+"_"+j)!=null) this.lma_dbg_offset[i][j]=Double.parseDouble(properties.getProperty(prefix+"lma_dbg_offset_"+i+"_"+j));
	    }
	}

	@Override
	public ImageDttParameters clone() throws CloneNotSupportedException {
        ImageDttParameters idp =     new ImageDttParameters();
		idp.gpu_mode_debug =         this.gpu_mode_debug;
		idp.gpu_verify =             this.gpu_verify;
		idp.corr_mode_debug =        this.corr_mode_debug;
		idp.mix_corr_poly =          this.mix_corr_poly;
		idp.corr_poly_only =         this.corr_poly_only;
		idp.min_poly_strength =      this.min_poly_strength;
		idp.max_poly_hwidth =        this.max_poly_hwidth;
		idp.poly_corr_scale =        this.poly_corr_scale;

		idp.poly_pwr =               this.poly_pwr;
		idp.poly_vasw_pwr =          this.poly_vasw_pwr;
		idp.corr_magic_scale_cm =    this.corr_magic_scale_cm;
		idp.corr_magic_scale_poly =  this.corr_magic_scale_poly;

		idp.ortho_height =           this.ortho_height;
		idp.ortho_eff_height =       this.ortho_eff_height;
		idp.ortho_nsamples =         this.ortho_nsamples;
		idp.ortho_vasw_pwr =         this.ortho_vasw_pwr;

		idp.enhortho_width =         this.enhortho_width;
		idp.enhortho_width_aux =     this.enhortho_width_aux;
		idp.enhortho_scale =         this.enhortho_scale;
		idp.enhortho_scale_aux =     this.enhortho_scale_aux;

		idp.ly_poly =                this.ly_poly;
		idp.ly_crazy_poly =          this.ly_crazy_poly;
		idp.ly_poly_backup =         this.ly_poly_backup;

		idp.fo_correct =             this.fo_correct;
		idp.fo_far =                 this.fo_far;
		idp.fo_min_strength =        this.fo_min_strength;
		idp.fo_min_eff =             this.fo_min_eff;
		idp.fo_min_eff_ratio =       this.fo_min_eff_ratio;
		idp.fo_max_hwidth =          this.fo_max_hwidth;
		idp.fo_min_diff =            this.fo_min_diff;
		idp.fo_ortho =               this.fo_ortho;
		idp.fo_overcorrection =      this.fo_overcorrection;
		idp.fo_lim_overcorr =        this.fo_lim_overcorr;

		idp.mismatch_max_diff =      this.mismatch_max_diff;

		idp.corr_offset=             this.corr_offset;
		idp.twice_diagonal=          this.twice_diagonal;
		idp.min_corr=                this.min_corr;

		idp.dbg_pair_mask=           this.dbg_pair_mask;
		idp.corr_strip_hight=        this.corr_strip_hight;

		idp.cnvx_hwnd_size=          this.cnvx_hwnd_size;
		idp.cnvx_weight=             this.cnvx_weight;
		idp.cnvx_add3x3=             this.cnvx_add3x3;

		idp.cnvx_min_samples=        this.cnvx_min_samples;
		idp.cnvx_non_coll=           this.cnvx_non_coll;
		idp.cnvx_min_pairs=          this.cnvx_min_pairs;

		idp.mcorr_multi=             this.mcorr_multi;
		idp.mcorr_all=               this.mcorr_all;
		idp.mcorr_all_multi=         this.mcorr_all_multi;
		idp.mcorr_dia=               this.mcorr_dia;
		idp.mcorr_dia_multi=         this.mcorr_dia_multi;
		idp.mcorr_sq=                this.mcorr_sq;
		idp.mcorr_sq_multi=          this.mcorr_sq_multi;
		idp.mcorr_neib=              this.mcorr_neib;
		idp.mcorr_neib_multi=        this.mcorr_neib_multi;
		
		idp.mcorr_hor=               this.mcorr_hor;
		idp.mcorr_hor_multi=         this.mcorr_hor_multi;
		idp.mcorr_vert=              this.mcorr_vert;
		idp.mcorr_vert_multi=        this.mcorr_vert_multi;
		
		idp.mcorr_cons_all=          this.mcorr_cons_all;
		idp.mcorr_cons_dia=          this.mcorr_cons_dia;
		idp.mcorr_cons_sq=           this.mcorr_cons_sq;
		idp.mcorr_cons_neib=         this.mcorr_cons_neib;
		idp.mcorr_cons_hor=          this.mcorr_cons_hor;
		idp.mcorr_cons_vert=         this.mcorr_cons_vert;
		
		idp.mcorr_comb_width=        this.mcorr_comb_width;
		idp.mcorr_comb_height=       this.mcorr_comb_height;
		idp.mcorr_comb_offset=       this.mcorr_comb_offset;
		idp.mcorr_comb_disp=         this.mcorr_comb_disp;
		idp.mcorr_comb_dbg=          this.mcorr_comb_dbg;
		
		idp.corr_strip_notch=        this.corr_strip_notch;
		idp.corr_notch_hwidth=       this.corr_notch_hwidth;
		idp.corr_notch_blur=         this.corr_notch_blur;

		idp.corr_wndy_size=          this.corr_wndy_size;

		idp.corr_wndy_hwidth =       this.corr_wndy_hwidth;
		idp.corr_wndy_blur =         this.corr_wndy_blur;

		idp.corr_wndx_size=          this.corr_wndx_size;

		idp.corr_wndx_hwidth =       this.corr_wndx_hwidth;
		idp.corr_wndx_blur =         this.corr_wndx_blur;

		idp.pcorr_use =              this.pcorr_use;
		idp.pcorr_use_hv =           this.pcorr_use_hv;
		idp.pcorr_sigma_mono =       this.pcorr_sigma_mono;
		idp.pcorr_sigma =            this.pcorr_sigma;
		idp.pcorr_sigma_rb =         this.pcorr_sigma_rb;
		idp.pcorr_fat_zero =         this.pcorr_fat_zero;
		idp.pcorr_fat_zero_mono =    this.pcorr_fat_zero_mono;
		idp.pcorr_dbg_offsx =        this.pcorr_dbg_offsx;
		idp.pcorr_dbg_offsy =        this.pcorr_dbg_offsy;

		idp.lma_disp_range=          this.lma_disp_range;
		idp.lmas_gaussian =          this.lmas_gaussian;
		idp.lmas_adjust_wm =         this.lmas_adjust_wm;
		idp.lmas_adjust_wy =         this.lmas_adjust_wy;
		idp.lmas_adjust_ag =         this.lmas_adjust_ag;

		idp.lmas_poly_str_scale =    this.lmas_poly_str_scale;
		idp.lmas_poly_str_min =      this.lmas_poly_str_min;

		idp.lmas_lambda_initial =    this.lmas_lambda_initial;
		idp.lmas_rms_diff =          this.lmas_rms_diff;
		idp.lmas_num_iter =          this.lmas_num_iter;
		idp.lmas_max_rel_rms=        this.lmas_max_rel_rms;
		idp.lmas_min_strength=       this.lmas_min_strength;
		idp.lmas_min_ac=             this.lmas_min_ac;
		idp.lmas_min_min_ac=         this.lmas_min_min_ac;
		idp.lmas_max_area=           this.lmas_max_area;

		idp.lma_gaussian =           this.lma_gaussian;
		idp.lma_second =             this.lma_second;
		idp.lma_second_gaussian =    this.lma_second_gaussian;
		idp.lma_adjust_wm =          this.lma_adjust_wm;
		idp.lma_adjust_wy =          this.lma_adjust_wy;
		idp.lma_adjust_wxy =         this.lma_adjust_wxy;
		idp.lma_adjust_ly1 =         this.lma_adjust_ly1;
		idp.lma_adjust_ag =          this.lma_adjust_ag;

		idp.lma_wnd =                this.lma_wnd;
		idp.lma_min_wnd =            this.lma_min_wnd;
		idp.lma_wnd_pwr =            this.lma_wnd_pwr;
		idp.lma_hard_marg =          this.lma_hard_marg;
		idp.lma_soft_marg =          this.lma_soft_marg;
		idp.lma_sigma =              this.lma_sigma;

		idp.lma_half_width =         this.lma_half_width;
		idp.lma_cost_wy =            this.lma_cost_wy;
		idp.lma_cost_wxy =           this.lma_cost_wxy;

		idp.lma_lambda_initial =     this.lma_lambda_initial;
		idp.lma_lambda_scale_good =  this.lma_lambda_scale_good;
		idp.lma_lambda_scale_bad =   this.lma_lambda_scale_bad;
		idp.lma_lambda_max =         this.lma_lambda_max;
		idp.lma_rms_diff =           this.lma_rms_diff;
		idp.lma_num_iter =           this.lma_num_iter;

		idp.lma_max_rel_rms=         this.lma_max_rel_rms;
		idp.lma_min_strength=        this.lma_min_strength;
		idp.lma_min_ac=              this.lma_min_ac;
		idp.lma_min_min_ac=          this.lma_min_min_ac;
		idp.lma_max_area=            this.lma_max_area;
		idp.lma_str_scale=           this.lma_str_scale;
		idp.lma_str_offset=          this.lma_str_offset;

		idp.lma_diff_xy =            this.lma_diff_xy;
		idp.lma_diff_minw =          this.lma_diff_minw;
		idp.lma_diff_sigma =         this.lma_diff_sigma;

		idp.lma_debug_level =        this.lma_debug_level;
		idp.lma_debug_level1 =       this.lma_debug_level1;
		idp.lma_debug_graphic =      this.lma_debug_graphic;

		idp.corr_var_cam =           this.corr_var_cam;
		idp.cm_max_normalization=    this.cm_max_normalization;

		idp.lma_dbg_scale=           this.lma_dbg_scale;
		idp.lma_dbg_offset=          new double [this.lma_dbg_offset.length][];
		for (int i = 0; i < idp.lma_dbg_offset.length; i++) {
			idp.lma_dbg_offset[i] =  this.lma_dbg_offset[i].clone();
		}

		return idp;
	}
}
