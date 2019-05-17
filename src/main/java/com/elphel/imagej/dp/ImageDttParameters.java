package com.elphel.imagej.dp;
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
	public boolean corr_mode_debug =        true;
	public boolean mix_corr_poly =          true;
	public double  min_poly_strength =      0.2;
	public double  max_poly_hwidth =        2.5; // Maximal polynomial approximation half-width (in both directions)
	public double  poly_corr_scale =        2.0; // Shift value if correlation maximum is wide in X than in Y to detect near objects (negative - far ones)

	public double  poly_pwr =               1.0;
	public double  poly_vasw_pwr =          2.0; // raise value to this power and apply as weight (0 - disable)
	public double  corr_magic_scale_cm =    1.0; //0.85;  // reported correlation offset vs. actual one (not yet understood)
	public double  corr_magic_scale_poly =  1.0; // 0.95;  // reported correlation offset vs. actual one (not yet understood)

	public int     ortho_height =             7;   // height of non-zero weights for hor/vert correlation to compensate borders
	public double  ortho_eff_height =         2.58; // effective correlation stripe height to match strengths

	public int     ortho_nsamples =           5; // number of samples to fit parabola
	public double  ortho_vasw_pwr =         2.0; // use data as weights when fitting parabola (high value samples are more important (when false use 3 samples only)

	public int     enhortho_width =         2;   // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
	public double  enhortho_scale =         0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)
	public boolean ly_poly =                false; // Use polynomial when measuring mismatch (false - use center of mass)
	public double  ly_crazy_poly =          1.0;  // Maximal allowed mismatch difference calculated as polynomial maximum
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

// LMA parameters
	public boolean lma_adjust_wm =          true;
	public boolean lma_adjust_wy =          false;
	public boolean lma_adjust_wxy =         true;
	public boolean lma_adjust_ag =          true;

	public double  lma_half_width =         2.0;   //
	public double  lma_cost_wy =            0.1;   //
	public double  lma_cost_wxy =           0.1;   //

	public double  lma_lambda_initial =     0.1;   //
	public double  lma_lambda_scale_good =  0.5;   //
	public double  lma_lambda_scale_bad =   8.0;   //
	public double  lma_lambda_max =       100.0;   //
	public double  lma_rms_diff =           0.001; //
	public int     lma_num_iter =          20;     //
	public int     lma_debug_level =        3;     //
	public boolean corr_var_cam =           true;  // New correlation mode compatible with 8 subcameras
	public double  cm_max_normalization =   0.55; // fraction of correlation maximum radius, being squared multiplied by maximum to have the same total mass


	public void dialogQuestions(GenericJTabbedDialog gd) {
			gd.addCheckbox    ("Enable ImageDtt correlation debug layers",                        this.corr_mode_debug,
					"false - return (old) per-coord correlations, true - replace them with more pairs correlation (new)");
			gd.addCheckbox    ("Replace CM layer with mixed/new poly one",                        this.mix_corr_poly);
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

			gd.addNumericField("Reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)",  this.enhortho_width,            0);
			gd.addNumericField("Multiply center correlation pixels (inside enhortho_width) (1.0 - disables enh_ortho)",  this.enhortho_scale,  3);

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

			gd.addMessage("Window for pole detection mode");
			gd.addNumericField("Strip height for pole detection",                                 this.corr_strip_notch,  0, 3, "half-pix",
					"Number of rows to combine/interpolate correlation results. Rows are twice denser than pixels correponding to largest baseline disparity");

			gd.addNumericField("50% correlation window cutoff height for poles (0 in the center)",this.corr_notch_hwidth,  3, 6, "half-pix",
					"Correlation window height argument for 50% value");
			gd.addNumericField("0% to 100 % transition range for poles",                          this.corr_notch_blur,  3,6,"half-pix",
					"Transition range, shifted sine is used");

			gd.addMessage("Window for niormal correlations");
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

			gd.addTab("Corr LMA","Parameters for LMA fitting of the correlation maximum parameters");
			gd.addCheckbox    ("Fit correlation defined half-width",                              this.lma_adjust_wm,
					"Allow fitting of the half-width common for all pairs, defined by the LPF filter of the phase correlation");
			gd.addCheckbox    ("Fit extra vertical half-width",                                   this.lma_adjust_wy,
					"Fit extra perpendicular to disparity half-width (not used? and only possible with multi-baseline cameras)");
			gd.addCheckbox    ("Fit extra half-width along disparity direction",                  this.lma_adjust_wxy,
					"Increased width in disparity direction caused by multi-distance objects in the tile");
			gd.addCheckbox    ("Adjust per-group amplitudes",                                     this.lma_adjust_ag,
					"Each correlation type's amplitude (now always needed)");

			gd.addNumericField("Initial/expected half-width of the correlation maximum in both directions", this.lma_half_width,  3, 6, "pix",
					"With LPF sigma = 0.9 it seems to be ~= 2.0. Used both as initial parameter and the fitted value difference from this may be penalized");
			gd.addNumericField("Cost of the difference of the actual half-width from the expected one",  this.lma_cost_wy,  5, 8, "",
					"The higher this cost, the more close the fitted half-width will be to the expected one");
			gd.addNumericField("Cost of the difference between horizontal and vertical widths",  this.lma_cost_wxy,  5, 8, "",
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
			gd.addNumericField("LMA debug level",                                                 this.lma_debug_level,  0, 3, "",
					"Debug/verbosity level for the LMA correaltion maximum fitting");
			gd.addCheckbox    ("Use new correlation methods compatible with x8 camera",           this.corr_var_cam,
					"Debug feature to compare old/new methods");

			gd.addNumericField("Normalization for the CM correlation strength",                  this.cm_max_normalization,  6, 8, "",
					"Fraction of correlation maximum radius, being squared multiplied by maximum to have the same total mass. ~= 0.5, the lower the value, the higher strength reported by the CM");
//	public double  cm_max_normalization =   0.55; //

	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
			this.corr_mode_debug=        gd.getNextBoolean();
			this.mix_corr_poly=          gd.getNextBoolean();
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
  			this.enhortho_scale=         gd.getNextNumber();

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

  			this.corr_strip_notch= (int) gd.getNextNumber();
  			this.corr_notch_hwidth=      gd.getNextNumber();
  			this.corr_notch_blur=        gd.getNextNumber();

  			this.corr_wndy_size=   (int) gd.getNextNumber();

			this.corr_wndy_hwidth =      gd.getNextNumber();
			this.corr_wndy_blur =        gd.getNextNumber();

  			this.corr_wndx_size=   (int) gd.getNextNumber();

			this.corr_wndx_hwidth =      gd.getNextNumber();
			this.corr_wndx_blur =        gd.getNextNumber();

//LMA tab
			this.lma_adjust_wm=          gd.getNextBoolean();
			this.lma_adjust_wy=          gd.getNextBoolean();
			this.lma_adjust_wxy=         gd.getNextBoolean();
			this.lma_adjust_ag=          gd.getNextBoolean();

			this.lma_half_width =        gd.getNextNumber();
			this.lma_cost_wy =           gd.getNextNumber();
			this.lma_cost_wxy =          gd.getNextNumber();

			this.lma_lambda_initial =    gd.getNextNumber();
			this.lma_lambda_scale_good = gd.getNextNumber();
			this.lma_lambda_scale_bad =  gd.getNextNumber();
			this.lma_lambda_max =        gd.getNextNumber();
			this.lma_rms_diff =          gd.getNextNumber();

  			this.lma_num_iter=     (int) gd.getNextNumber();
  			this.lma_debug_level=  (int) gd.getNextNumber();
  			this.corr_var_cam =          gd.getNextBoolean();
  			this.cm_max_normalization=   gd.getNextNumber();

	}


	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"corr_mode_debug",      this.corr_mode_debug+"");
		properties.setProperty(prefix+"mix_corr_poly",        this.mix_corr_poly+"");
		properties.setProperty(prefix+"min_poly_strength",    this.min_poly_strength+"");
		properties.setProperty(prefix+"max_poly_hwidth",      this.max_poly_hwidth+"");
		properties.setProperty(prefix+"poly_corr_scale",      this.poly_corr_scale+"");

		properties.setProperty(prefix+"poly_pwr",             this.poly_pwr+"");
		properties.setProperty(prefix+"poly_value_to_weight", this.poly_vasw_pwr+"");
		properties.setProperty(prefix+"corr_magic_scale_cm",  this.corr_magic_scale_cm+"");
		properties.setProperty(prefix+"corr_magic_scale_poly",this.corr_magic_scale_poly+"");

		properties.setProperty(prefix+"ortho_height",         this.ortho_height+"");
		properties.setProperty(prefix+"ortho_eff_height",     this.ortho_eff_height+"");
		properties.setProperty(prefix+"ortho_nsamples",       this.ortho_nsamples+"");
		properties.setProperty(prefix+"ortho_vasw",           this.ortho_vasw_pwr+"");

		properties.setProperty(prefix+"enhortho_width",       this.enhortho_width +"");
		properties.setProperty(prefix+"enhortho_scale",       this.enhortho_scale +"");

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


		properties.setProperty(prefix+"corr_strip_notch",     this.corr_strip_notch +"");
		properties.setProperty(prefix+"corr_notch_hwidth",    this.corr_notch_hwidth +"");
		properties.setProperty(prefix+"corr_notch_blur",      this.corr_notch_blur +"");

		properties.setProperty(prefix+"corr_wndy_size",       this.corr_wndy_size +"");
		properties.setProperty(prefix+"corr_wndy_hwidth",     this.corr_wndy_hwidth +"");
		properties.setProperty(prefix+"corr_wndy_blur",       this.corr_wndy_blur +"");

		properties.setProperty(prefix+"corr_wndx_size",       this.corr_wndx_size +"");

		properties.setProperty(prefix+"corr_wndx_hwidth",     this.corr_wndx_hwidth +"");
		properties.setProperty(prefix+"corr_wndx_blur",       this.corr_wndx_blur +"");



		properties.setProperty(prefix+"lma_adjust_wm",        this.lma_adjust_wm +"");
		properties.setProperty(prefix+"lma_adjust_wy",        this.lma_adjust_wy +"");
		properties.setProperty(prefix+"lma_adjust_wxy",       this.lma_adjust_wxy +"");
		properties.setProperty(prefix+"lma_adjust_ag",        this.lma_adjust_ag +"");

		properties.setProperty(prefix+"lma_half_width",       this.lma_half_width +"");
		properties.setProperty(prefix+"lma_cost_wy",          this.lma_cost_wy +"");
		properties.setProperty(prefix+"lma_cost_wxy",         this.lma_cost_wxy +"");

		properties.setProperty(prefix+"lma_lambda_initial",   this.lma_lambda_initial +"");
		properties.setProperty(prefix+"lma_lambda_scale_good",this.lma_lambda_scale_good +"");
		properties.setProperty(prefix+"lma_lambda_scale_bad", this.lma_lambda_scale_bad +"");
		properties.setProperty(prefix+"lma_lambda_max",       this.lma_lambda_max +"");
		properties.setProperty(prefix+"lma_rms_diff",         this.lma_rms_diff +"");

		properties.setProperty(prefix+"lma_num_iter",         this.lma_num_iter +"");
		properties.setProperty(prefix+"lma_debug_level",      this.lma_debug_level +"");
		properties.setProperty(prefix+"corr_var_cam",         this.corr_var_cam +"");

		properties.setProperty(prefix+"cm_max_normalization", this.cm_max_normalization +"");

	}

	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"corr_mode_debug")!=null)       this.corr_mode_debug=Boolean.parseBoolean(properties.getProperty(prefix+"corr_mode_debug"));
		if (properties.getProperty(prefix+"mix_corr_poly")!=null)         this.mix_corr_poly=Boolean.parseBoolean(properties.getProperty(prefix+"mix_corr_poly"));
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
		if (properties.getProperty(prefix+"enhortho_scale")!=null)        this.enhortho_scale=Double.parseDouble(properties.getProperty(prefix+"enhortho_scale"));

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



		if (properties.getProperty(prefix+"corr_strip_notch")!=null)     this.corr_strip_notch=Integer.parseInt(properties.getProperty(prefix+"corr_strip_notch"));
		if (properties.getProperty(prefix+"corr_notch_hwidth")!=null)    this.corr_notch_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_notch_hwidth"));
		if (properties.getProperty(prefix+"corr_notch_blur")!=null)      this.corr_notch_blur=Double.parseDouble(properties.getProperty(prefix+"corr_notch_blur"));

		if (properties.getProperty(prefix+"corr_wndy_size")!=null)       this.corr_wndy_size=Integer.parseInt(properties.getProperty(prefix+"corr_wndy_size"));

		if (properties.getProperty(prefix+"corr_wndy_hwidth")!=null)     this.corr_wndy_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_wndy_hwidth"));
		if (properties.getProperty(prefix+"corr_wndy_blur")!=null)       this.corr_wndy_blur=Double.parseDouble(properties.getProperty(prefix+"corr_wndy_blur"));

		if (properties.getProperty(prefix+"corr_wndx_size")!=null)       this.corr_wndx_size=Integer.parseInt(properties.getProperty(prefix+"corr_wndx_size"));

		if (properties.getProperty(prefix+"corr_wndx_hwidth")!=null)     this.corr_wndx_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_wndx_hwidth"));
		if (properties.getProperty(prefix+"corr_wndx_blur")!=null)       this.corr_wndx_blur=Double.parseDouble(properties.getProperty(prefix+"corr_wndx_blur"));



		if (properties.getProperty(prefix+"lma_adjust_wm")!=null)        this.lma_adjust_wm=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_wm"));
		if (properties.getProperty(prefix+"lma_adjust_wy")!=null)        this.lma_adjust_wy=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_wy"));
		if (properties.getProperty(prefix+"lma_adjust_wxy")!=null)       this.lma_adjust_wxy=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_wxy"));
		if (properties.getProperty(prefix+"lma_adjust_ag")!=null)        this.lma_adjust_ag=Boolean.parseBoolean(properties.getProperty(prefix+"lma_adjust_ag"));

		if (properties.getProperty(prefix+"lma_half_width")!=null)       this.lma_half_width=Double.parseDouble(properties.getProperty(prefix+"lma_half_width"));
		if (properties.getProperty(prefix+"lma_cost_wy")!=null)          this.lma_cost_wy=Double.parseDouble(properties.getProperty(prefix+"lma_cost_wy"));
		if (properties.getProperty(prefix+"lma_cost_wxy")!=null)         this.lma_cost_wxy=Double.parseDouble(properties.getProperty(prefix+"lma_cost_wxy"));

		if (properties.getProperty(prefix+"lma_lambda_initial")!=null)   this.lma_lambda_initial=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_initial"));
		if (properties.getProperty(prefix+"lma_lambda_scale_good")!=null)this.lma_lambda_scale_good=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_scale_good"));
		if (properties.getProperty(prefix+"lma_lambda_scale_bad")!=null) this.lma_lambda_scale_bad=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_scale_bad"));
		if (properties.getProperty(prefix+"lma_lambda_max")!=null)       this.lma_lambda_max=Double.parseDouble(properties.getProperty(prefix+"lma_lambda_max"));
		if (properties.getProperty(prefix+"lma_rms_diff")!=null)         this.lma_rms_diff=Double.parseDouble(properties.getProperty(prefix+"lma_rms_diff"));

		if (properties.getProperty(prefix+"lma_num_iter")!=null)         this.lma_num_iter=Integer.parseInt(properties.getProperty(prefix+"lma_num_iter"));
		if (properties.getProperty(prefix+"lma_debug_level")!=null)      this.lma_debug_level=Integer.parseInt(properties.getProperty(prefix+"lma_debug_level"));
		if (properties.getProperty(prefix+"corr_var_cam")!=null)         this.corr_var_cam=Boolean.parseBoolean(properties.getProperty(prefix+"corr_var_cam"));

		if (properties.getProperty(prefix+"cm_max_normalization")!=null) this.cm_max_normalization=Double.parseDouble(properties.getProperty(prefix+"cm_max_normalization"));
	}

	@Override
	public ImageDttParameters clone() throws CloneNotSupportedException {
        ImageDttParameters idp =     new ImageDttParameters();
		idp.corr_mode_debug =        this.corr_mode_debug;
		idp.mix_corr_poly =          this.mix_corr_poly;
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
		idp.enhortho_scale =         this.enhortho_scale;

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

		idp.corr_strip_notch=        this.corr_strip_notch;
		idp.corr_notch_hwidth=       this.corr_notch_hwidth;
		idp.corr_notch_blur=         this.corr_notch_blur;

		idp.corr_wndy_size=          this.corr_wndy_size;

		idp.corr_wndy_hwidth =       this.corr_wndy_hwidth;
		idp.corr_wndy_blur =         this.corr_wndy_blur;

		idp.corr_wndx_size=          this.corr_wndx_size;

		idp.corr_wndx_hwidth =       this.corr_wndx_hwidth;
		idp.corr_wndx_blur =         this.corr_wndx_blur;

		idp.lma_adjust_wm =          this.lma_adjust_wm;
		idp.lma_adjust_wy =          this.lma_adjust_wy;
		idp.lma_adjust_wxy =         this.lma_adjust_wxy;
		idp.lma_adjust_ag =          this.lma_adjust_ag;

		idp.lma_half_width =         this.lma_half_width;
		idp.lma_cost_wy =            this.lma_cost_wy;
		idp.lma_cost_wxy =           this.lma_cost_wxy;

		idp.lma_lambda_initial =     this.lma_lambda_initial;
		idp.lma_lambda_scale_good =  this.lma_lambda_scale_good;
		idp.lma_lambda_scale_bad =   this.lma_lambda_scale_bad;
		idp.lma_lambda_max =         this.lma_lambda_max;
		idp.lma_rms_diff =           this.lma_rms_diff;

		idp.lma_num_iter =           this.lma_num_iter;
		idp.lma_debug_level =        this.lma_debug_level;
		idp.corr_var_cam =           this.corr_var_cam;

		idp.cm_max_normalization=    this.cm_max_normalization;

		return idp;
	}
}
