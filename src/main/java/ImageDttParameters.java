/**
 **
 ** ImageDttParameters - parameters defining TP operations (at first extra)
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ImageDtt.java is free software: you can redistribute it and/or modify
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

public class ImageDttParameters {
	public boolean corr_mode_debug =        true;
	public boolean mix_corr_poly =          false;
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
	public double  ortho_vasw_pwr =         2.0; // use data as weights when fitting parabola (high vale samples are more important (when false use 3 samples only)

	public int     enhortho_width =         2;   // reduce weight of center correlation pixels from center (0 - none, 1 - center, 2 +/-1 from center)
	public double  enhortho_scale =         0.0; // 0.2;  // multiply center correlation pixels (inside enhortho_width)

	public boolean far_object_correct =     true; // correct far objects by comparing orthogonal correlations
	public double  fo_min_strength =        0.3;  // minimal strength for all correlations (full, hor, vert) to try to correct
	public double  fo_max_hwidth =          2.4;  // maximal half-width for all directions to try to correct far objects
	public double  fo_overcorrection =      0.0;  // add scaled hor/vert difference to the largest of hor/vert disparity
	public double  fo_lim_overcorr =        2.0;  // limit full correction with respect to largest - fullcorr difference

	public double  corr_offset =            0.1; //0.1;  // add to pair correlation before multiplying by other pairs
	public double  min_corr =               0.02; // minimal correlation value to consider valid

	public int     dbg_pair_mask =          0x3ff; // which pairs to combine
	public int     corr_strip_hight =       9;     // number of rows to calculate
	public int     corr_wndy_size =         9;     // number of rows to calculate CM disparity
	public double  corr_wndy_hwidth =       6.0;   // 50% window cutoff height
	public double  corr_wndy_blur =         5.0;   // 100% to 0 % vertical transition range
	public int     corr_wndx_size =         9;     // half number of columns to calculate CM disparity (each row has only odd/even columns, so disparity range is smaller
	public double  corr_wndx_hwidth =       6.0;   // 50% window cutoff width
	public double  corr_wndx_blur =         5.0;   // 100% to 0 % vertical transition range





	public void dialogQuestions(GenericJTabbedDialog gd) {
			gd.addCheckbox    ("Enable ImageDtt correlation debug layers",                        this.corr_mode_debug);
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
			gd.addNumericField("Multiply center correlation pixels (inside enhortho_width) (1.0 - disables enh_orttho)",  this.enhortho_scale,  3);

			gd.addMessage("Far objects correction");
			gd.addCheckbox    ("Try to correct far objects (make them closer) by hor/vert comparison",      this.far_object_correct);
			gd.addNumericField("Minimal strength for all correlations (full, hor, vert) to try to correct", this.fo_min_strength,  3,6,"",
					"Do not correct if any of the full, hor or vert strength is lower than this");
			gd.addNumericField("Maximal correlation half-width",                                            this.fo_max_hwidth,  3,6,"",
					"Do not correct if any of the full, hor or vert half-width is above this");
			gd.addNumericField("Overcorrection scale",                                            this.fo_overcorrection,  3,6,"",
					"Add scaled hor/vert difference to the largest of hor/vert disparity. Use 0 to just select largest of disparities");
			gd.addNumericField("Limit overcorrection",                                            this.fo_lim_overcorr,  3,6,"",
					"Limit full correction with respect to largest - fullcorr difference. 1.0 does not allow overcorrection, < 1.0 - the result will be closer to full correction");


	  		gd.addNumericField("Add to pair correlation before multiplying by other pairs",       this.corr_offset,  6,8,"",
	  				"0.0 - pure product (false-positive tolerant), higher the value - closer to the sum (linear, better S/N");

		    gd.addNumericField("Extract disparity max/argmax if maximal value is above",          this.min_corr,  6,8,"",
		    		"skip unreliable correlations");


			gd.addNumericField("Debug: which pairs to combine",                                   this.dbg_pair_mask,  0, 3, "",
					"Bits: 0, 1 - horizontal pairs, 2,3 - verical pairs, 4,5 - diagonal pairs");
			gd.addNumericField("Number of correlation rows to combine (strip height)",            this.corr_strip_hight,  0, 3, "",
					"Number of rows to combine/interpolate correlation results. Rows are twice denser than pixels correponding to largest baseline disparity");
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

  			this.far_object_correct =    gd.getNextBoolean();
  			this.fo_min_strength =       gd.getNextNumber();
  			this.fo_max_hwidth =         gd.getNextNumber();
  			this.fo_overcorrection =     gd.getNextNumber();
  			this.fo_lim_overcorr =       gd.getNextNumber();

	  		this.corr_offset =           gd.getNextNumber();
		    this.min_corr =              gd.getNextNumber();


  			this.dbg_pair_mask=    (int) gd.getNextNumber();
  			this.corr_strip_hight= (int) gd.getNextNumber();
  			this.corr_wndy_size=   (int) gd.getNextNumber();

			this.corr_wndy_hwidth =      gd.getNextNumber();
			this.corr_wndy_blur =        gd.getNextNumber();

  			this.corr_wndx_size=   (int) gd.getNextNumber();

			this.corr_wndx_hwidth =      gd.getNextNumber();
			this.corr_wndx_blur =        gd.getNextNumber();

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
		properties.setProperty(prefix+"min_corr",             this.min_corr +"");

		properties.setProperty(prefix+"far_object_correct",   this.far_object_correct +"");
		properties.setProperty(prefix+"fo_min_strength",      this.fo_min_strength +"");
		properties.setProperty(prefix+"fo_max_hwidth",        this.fo_max_hwidth +"");
		properties.setProperty(prefix+"fo_overcorrection",    this.fo_overcorrection +"");
		properties.setProperty(prefix+"fo_lim_overcorr",      this.fo_lim_overcorr +"");

		properties.setProperty(prefix+"dbg_pair_mask",        this.dbg_pair_mask +"");
		properties.setProperty(prefix+"corr_strip_hight",     this.corr_strip_hight +"");
		properties.setProperty(prefix+"corr_wndy_size",       this.corr_wndy_size +"");

		properties.setProperty(prefix+"corr_wndy_hwidth",     this.corr_wndy_hwidth +"");
		properties.setProperty(prefix+"corr_wndy_blur",       this.corr_wndy_blur +"");

		properties.setProperty(prefix+"corr_wndx_size",       this.corr_wndx_size +"");

		properties.setProperty(prefix+"corr_wndx_hwidth",     this.corr_wndx_hwidth +"");
		properties.setProperty(prefix+"corr_wndx_blur",       this.corr_wndx_blur +"");

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

		if (properties.getProperty(prefix+"far_object_correct")!=null)    this.far_object_correct=Boolean.parseBoolean(properties.getProperty(prefix+"far_object_correct"));
		if (properties.getProperty(prefix+"fo_min_strength")!=null)       this.fo_min_strength=Double.parseDouble(properties.getProperty(prefix+"fo_min_strength"));
		if (properties.getProperty(prefix+"fo_max_hwidth")!=null)         this.fo_max_hwidth=Double.parseDouble(properties.getProperty(prefix+"fo_max_hwidth"));
		if (properties.getProperty(prefix+"fo_overcorrection")!=null)     this.fo_overcorrection=Double.parseDouble(properties.getProperty(prefix+"fo_overcorrection"));
		if (properties.getProperty(prefix+"fo_lim_overcorr")!=null)       this.fo_lim_overcorr=Double.parseDouble(properties.getProperty(prefix+"fo_lim_overcorr"));


		if (properties.getProperty(prefix+"corr_offset")!=null)           this.corr_offset=Double.parseDouble(properties.getProperty(prefix+"corr_offset"));
		if (properties.getProperty(prefix+"min_corr")!=null)              this.min_corr=Double.parseDouble(properties.getProperty(prefix+"min_corr"));


		if (properties.getProperty(prefix+"dbg_pair_mask")!=null)        this.dbg_pair_mask=Integer.parseInt(properties.getProperty(prefix+"dbg_pair_mask"));
		if (properties.getProperty(prefix+"corr_strip_hight")!=null)     this.corr_strip_hight=Integer.parseInt(properties.getProperty(prefix+"corr_strip_hight"));
		if (properties.getProperty(prefix+"corr_wndy_size")!=null)       this.corr_wndy_size=Integer.parseInt(properties.getProperty(prefix+"corr_wndy_size"));

		if (properties.getProperty(prefix+"corr_wndy_hwidth")!=null)     this.corr_wndy_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_wndy_hwidth"));
		if (properties.getProperty(prefix+"corr_wndy_blur")!=null)       this.corr_wndy_blur=Double.parseDouble(properties.getProperty(prefix+"corr_wndy_blur"));

		if (properties.getProperty(prefix+"corr_wndx_size")!=null)       this.corr_wndx_size=Integer.parseInt(properties.getProperty(prefix+"corr_wndx_size"));

		if (properties.getProperty(prefix+"corr_wndx_hwidth")!=null)     this.corr_wndx_hwidth=Double.parseDouble(properties.getProperty(prefix+"corr_wndx_hwidth"));
		if (properties.getProperty(prefix+"corr_wndx_blur")!=null)       this.corr_wndx_blur=Double.parseDouble(properties.getProperty(prefix+"corr_wndx_blur"));

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

		idp.far_object_correct =     this.far_object_correct;
		idp.fo_min_strength =        this.fo_min_strength;
		idp.fo_max_hwidth =          this.fo_max_hwidth;
		idp.fo_overcorrection =      this.fo_overcorrection;
		idp.fo_lim_overcorr =        this.fo_lim_overcorr;

		idp.corr_offset=             this.corr_offset;
		idp.min_corr=                this.min_corr;

		idp.dbg_pair_mask=           this.dbg_pair_mask;
		idp.corr_strip_hight=        this.corr_strip_hight;
		idp.corr_wndy_size=          this.corr_wndy_size;

		idp.corr_wndy_hwidth =       this.corr_wndy_hwidth;
		idp.corr_wndy_blur =         this.corr_wndy_blur;

		idp.corr_wndx_size=          this.corr_wndx_size;

		idp.corr_wndx_hwidth =       this.corr_wndx_hwidth;
		idp.corr_wndx_blur =         this.corr_wndx_blur;

		return idp;
	}

	// TODO move 2 dialog methods here
}
