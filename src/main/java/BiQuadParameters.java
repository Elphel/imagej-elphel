import java.util.Properties;

/**
 **
 ** BiQuadParameters - parameters defining Operation of a two quad camera rig
 **
 ** Copyright (C) 2018 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  BiQuadParameters.java is free software: you can redistribute it and/or modify
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

public class BiQuadParameters {
	public boolean rig_mode_debug =            true;
	public boolean use_poly =                  true;
	public boolean use_xy_poly =               true;

	public double  min_poly_strength =         0.2;
	public double  min_xy_poly_strength =      1.0; // never
	public double  inf_min_strength_main =     0.12;
	public double  inf_min_strength_aux =      0.12;
	public double  inf_min_strength_rig =      0.25;
	public double  inf_max_disp_main =         0.15;
	public double  inf_max_disp_aux =          0.15;
	public double  inf_max_disp_rig =          0.5; // maybe even higher (2.0) to lock to initially high mismatch

	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addCheckbox    ("Debug rig/bi-camera functionality ",                              this.rig_mode_debug,"Enable debugging of the methods related to dual camera rig");
		gd.addCheckbox    ("Use poly for main/aux/rig correlations (false - CM)",             this.use_poly,"Use LMA/polynomial if correlation is strong enough");
		gd.addCheckbox    ("Use poly for rig X/Y mismatch measurements",                      this.use_xy_poly,"Use polynomial forX/Y offset measurements if correlation is strong enough");

		gd.addNumericField("Use poly mode for disparity if strength is greater than",                             this.min_poly_strength,  3,6,"",
				"Minimal correlation strength to use poly mode (below use CM)");
		gd.addNumericField("Use poly mode for X/Y rig correction  if strength is greater than",                   this.min_xy_poly_strength,  3,6,"",
				"Minimal correlation strength to use poly mode (below use CM) for measureking X/Y offset at infinity");
		gd.addNumericField("Minimal strength for main camera correlation to use for infinity rig adjustment",     this.inf_min_strength_main,  3,6,"",
				"Do not use tile for infinity adjustment if main correlation strength is less than");
		gd.addNumericField("Minimal strength for aux camera correlation to use for infinity rig adjustment",      this.inf_min_strength_aux,  3,6,"",
				"Do not use tile for infinity adjustment if aux correlation strength is less than");
		gd.addNumericField("Minimal strength for inter-camera correlation to use for infinity rig adjustment",    this.inf_min_strength_rig,  3,6,"",
				"Do not use tile for infinity adjustment if inter-camera correlation strength is less than");
		gd.addNumericField("Maximal absolute value of main camera disparity to use for infinity rig adjustment",  this.inf_max_disp_main,  3,6,"pix",
				"Do not use tile for infinity adjustment if absolute value of the main camera disparity is too high");
		gd.addNumericField("Maximal absolute value of aux camera disparity to use for infinity rig adjustment",   this.inf_max_disp_aux,  3,6,"pix",
				"Do not use tile for infinity adjustment if absolute value of the main camera disparity is too high");
		gd.addNumericField("Maximal absolute value of inter-camera disparity to use for infinity rig adjustment", this.inf_max_disp_rig,  3,6,"pix",
				"Do not use tile for infinity adjustment if absolute value of the inter-camera disparity is too high");
	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
		this.rig_mode_debug=         gd.getNextBoolean();
		this.use_poly=               gd.getNextBoolean();
		this.use_xy_poly=            gd.getNextBoolean();
		this.min_poly_strength=      gd.getNextNumber();
		this.min_xy_poly_strength=   gd.getNextNumber();
		this.inf_min_strength_main=  gd.getNextNumber();
		this.inf_min_strength_aux=   gd.getNextNumber();
		this.inf_min_strength_rig=   gd.getNextNumber();
		this.inf_max_disp_main=      gd.getNextNumber();
		this.inf_max_disp_aux=       gd.getNextNumber();
		this.inf_max_disp_rig=       gd.getNextNumber();
	}
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"rig_mode_debug",        this.rig_mode_debug+"");
		properties.setProperty(prefix+"use_poly",              this.use_poly+"");
		properties.setProperty(prefix+"use_xy_poly",           this.use_xy_poly+"");
		properties.setProperty(prefix+"min_poly_strength",     this.min_poly_strength+"");
		properties.setProperty(prefix+"min_xy_poly_strength",  this.min_xy_poly_strength+"");
		properties.setProperty(prefix+"inf_min_strength_main", this.inf_min_strength_main+"");
		properties.setProperty(prefix+"inf_min_strength_aux",  this.inf_min_strength_aux+"");
		properties.setProperty(prefix+"inf_min_strength_rig",  this.inf_min_strength_rig+"");
		properties.setProperty(prefix+"inf_max_disp_main",     this.inf_max_disp_main+"");
		properties.setProperty(prefix+"inf_max_disp_aux",      this.inf_max_disp_aux+"");
		properties.setProperty(prefix+"inf_max_disp_rig",      this.inf_max_disp_rig+"");
	}
	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"rig_mode_debug")!=null)        this.rig_mode_debug=Boolean.parseBoolean(properties.getProperty(prefix+"rig_mode_debug"));
		if (properties.getProperty(prefix+"use_poly")!=null)              this.use_poly=Boolean.parseBoolean(properties.getProperty(prefix+"use_poly"));
		if (properties.getProperty(prefix+"use_xy_poly")!=null)           this.use_xy_poly=Boolean.parseBoolean(properties.getProperty(prefix+"use_xy_poly"));
		if (properties.getProperty(prefix+"min_poly_strength")!=null)     this.min_poly_strength=Double.parseDouble(properties.getProperty(prefix+"min_poly_strength"));
		if (properties.getProperty(prefix+"min_xy_poly_strength")!=null)  this.min_xy_poly_strength=Double.parseDouble(properties.getProperty(prefix+"min_xy_poly_strength"));
		if (properties.getProperty(prefix+"inf_min_strength_main")!=null) this.inf_min_strength_main=Double.parseDouble(properties.getProperty(prefix+"inf_min_strength_main"));
		if (properties.getProperty(prefix+"inf_min_strength_aux")!=null)  this.inf_min_strength_aux=Double.parseDouble(properties.getProperty(prefix+"inf_min_strength_aux"));
		if (properties.getProperty(prefix+"inf_min_strength_rig")!=null)  this.inf_min_strength_rig=Double.parseDouble(properties.getProperty(prefix+"inf_min_strength_rig"));
		if (properties.getProperty(prefix+"inf_max_disp_main")!=null)     this.inf_max_disp_main=Double.parseDouble(properties.getProperty(prefix+"inf_max_disp_main"));
		if (properties.getProperty(prefix+"inf_max_disp_aux")!=null)      this.inf_max_disp_aux=Double.parseDouble(properties.getProperty(prefix+"inf_max_disp_aux"));
		if (properties.getProperty(prefix+"inf_max_disp_rig")!=null)      this.inf_max_disp_rig=Double.parseDouble(properties.getProperty(prefix+"inf_max_disp_rig"));
	}
	@Override
	public BiQuadParameters clone() throws CloneNotSupportedException {
		BiQuadParameters bqp =       new BiQuadParameters();
		bqp.rig_mode_debug=          this.rig_mode_debug;
		bqp.use_poly=                this.use_poly;
		bqp.use_xy_poly=             this.use_xy_poly;
		bqp.min_poly_strength =      this.min_poly_strength;
		bqp.min_xy_poly_strength =   this.min_xy_poly_strength;
		bqp.inf_min_strength_main =  this.inf_min_strength_main;
		bqp.inf_min_strength_aux =   this.inf_min_strength_aux;
		bqp.inf_min_strength_rig =   this.inf_min_strength_rig;
		bqp.inf_max_disp_main =      this.inf_max_disp_main;
		bqp.inf_max_disp_aux =       this.inf_max_disp_aux;
		bqp.inf_max_disp_rig =       this.inf_max_disp_rig;
		return bqp;


	}
}
