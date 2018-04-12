import java.util.Properties;

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

public class ImageDttParameters {
	public boolean corr_mode_debug =        true;
	public boolean mix_corr_poly =          false;
	public double  min_poly_strength =      0.2;
	public double  max_poly_diff =          0.6;
	public double  poly_pwr =               1.0;
	public boolean poly_value_to_weight =   true;

	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"corr_mode_debug",   this.corr_mode_debug+"");
		properties.setProperty(prefix+"mix_corr_poly",     this.mix_corr_poly+"");
		properties.setProperty(prefix+"min_poly_strength", this.min_poly_strength+"");
		properties.setProperty(prefix+"max_poly_diff",     this.max_poly_diff+"");
		properties.setProperty(prefix+"poly_pwr",          this.poly_pwr+"");
		properties.setProperty(prefix+"mix_corr_poly",     this.mix_corr_poly+"");
	}

	public void getProperties(String prefix,Properties properties){
		if (properties.getProperty(prefix+"corr_mode_debug")!=null)      this.corr_mode_debug=Boolean.parseBoolean(properties.getProperty(prefix+"corr_mode_debug"));
		if (properties.getProperty(prefix+"mix_corr_poly")!=null)        this.mix_corr_poly=Boolean.parseBoolean(properties.getProperty(prefix+"mix_corr_poly"));
		if (properties.getProperty(prefix+"min_poly_strength")!=null)    this.min_poly_strength=Double.parseDouble(properties.getProperty(prefix+"min_poly_strength"));
		if (properties.getProperty(prefix+"max_poly_diff")!=null)        this.max_poly_diff=Double.parseDouble(properties.getProperty(prefix+"max_poly_diff"));
		if (properties.getProperty(prefix+"poly_pwr")!=null)             this.poly_pwr=Double.parseDouble(properties.getProperty(prefix+"poly_pwr"));
		if (properties.getProperty(prefix+"poly_value_to_weight")!=null) this.poly_value_to_weight=Boolean.parseBoolean(properties.getProperty(prefix+"poly_value_to_weight"));
	}
	@Override
	public ImageDttParameters clone() {
		ImageDttParameters idp = new ImageDttParameters();
		idp.corr_mode_debug =        corr_mode_debug;
		idp.mix_corr_poly =          mix_corr_poly;
		idp.min_poly_strength =      min_poly_strength;
		idp.max_poly_diff =          max_poly_diff;
		idp.poly_pwr =               poly_pwr;
		idp.poly_value_to_weight =   poly_value_to_weight;
		return idp;
	}
}
