package com.elphel.imagej.tileprocessor;

import java.util.Properties;

import com.elphel.imagej.common.GenericJTabbedDialog;

public class IntersceneLmaParameters {
	public boolean [] ilma_lma_select =             new boolean [ErsCorrection.DP_NUM_PARS]; // first three will not be used
	public double  [] ilma_regularization_weights = new double  [ErsCorrection.DP_NUM_PARS]; // first three will not be used
	public double     ilma_lambda = 0.1;
	public double     ilma_lambda_scale_good = 0.5;
	public double     ilma_lambda_scale_bad =  8.0;
	public double     ilma_lambda_max =      100;
	public double     ilma_rms_diff =          0.001;
	public int        ilma_num_iter =          20;
	public int        ilma_debug_level =       1;

	public IntersceneLmaParameters() {
		ilma_lma_select[ErsCorrection.DP_DSAZ]= true;
		ilma_lma_select[ErsCorrection.DP_DSTL]= true;
		ilma_lma_select[ErsCorrection.DP_DSRL]= true;
		ilma_lma_select[ErsCorrection.DP_DSX]=  true;
		ilma_lma_select[ErsCorrection.DP_DSY]=  true;
		ilma_lma_select[ErsCorrection.DP_DSZ]=  true;
	}
	
	public void dialogQuestions(GenericJTabbedDialog gd) {
		gd.addMessage("Interframe LMA parameters selection");
		for (int i = ErsCorrection.DP_DVAZ; i < ErsCorrection.DP_NUM_PARS; i++) {
		    gd.addCheckbox    (ErsCorrection.DP_DERIV_NAMES[i],           this.ilma_lma_select[i],
		    		"Adjust parameter "+ErsCorrection.DP_DERIV_NAMES[i]+" with interscene LMA" );
		}
		gd.addMessage("Regularization parameters - pull strength to the initial values");
		for (int i = ErsCorrection.DP_DVAZ; i < ErsCorrection.DP_NUM_PARS; i++) {
			gd.addNumericField(ErsCorrection.DP_DERIV_NAMES[i],           this.ilma_regularization_weights[i], 6,8,"",
			"Weight of "+ErsCorrection.DP_DERIV_NAMES[i]+" pull, 1.0 means that the paramter offset from initial corresponding to 1 image pixel\n"+
			" will cause error equal to all reprojection ones");
		}
		gd.addMessage("LMA other parameters");
		gd.addNumericField("LMA lambda",                                  this.ilma_lambda, 6,8,"",
				"Initial value of the LMA lambda");
		gd.addNumericField("Scale lambda after successful LMA iteration", this.ilma_lambda_scale_good, 3,5,"",
				"Scale lambda (reduce) if the new RMSE is lower than the previous one.");
		gd.addNumericField("Scale lambda after failed LMA iteration",     this.ilma_lambda_scale_bad, 3,5,"",
				"Scale lambda (increase) if the new RMSE is higher than the previous one.");
		gd.addNumericField("Maximal value of lambda to try",              this.ilma_lambda_max, 2,7,"",
				"Fail LMA if the result is still worse than before parameters were updates.");
		gd.addNumericField("Minimal relative RMSE improvement",           this.ilma_rms_diff, 5,7,"",
				"Exit LMA iterations if relative RMSE improvement drops below this value.");
		gd.addNumericField("Maximal number of LMA iterations",            this.ilma_num_iter, 0,3,"",
				"A hard limit on LMA iterations.");
		gd.addNumericField("Debug level",                                 this.ilma_debug_level, 0,3,"",
				"Debug level of interscene LMA operation.");
	}
	public void dialogAnswers(GenericJTabbedDialog gd) {
		for (int i = ErsCorrection.DP_DVAZ; i < ErsCorrection.DP_NUM_PARS; i++) {
		    this.ilma_lma_select[i] =             gd.getNextBoolean();
		}
		for (int i = ErsCorrection.DP_DVAZ; i < ErsCorrection.DP_NUM_PARS; i++) {
			this.ilma_regularization_weights[i] = gd.getNextNumber();
		}
		this.ilma_lambda =                        gd.getNextNumber();
		this.ilma_lambda_scale_good =             gd.getNextNumber();
		this.ilma_lambda_scale_bad =              gd.getNextNumber();
		this.ilma_lambda_max =                    gd.getNextNumber();
		this.ilma_rms_diff =                      gd.getNextNumber();
		this.ilma_num_iter =                (int) gd.getNextNumber();
		this.ilma_debug_level =             (int) gd.getNextNumber();
	}
	public void setProperties(String prefix,Properties properties){
		for (int i = ErsCorrection.DP_DVAZ; i < ErsCorrection.DP_NUM_PARS; i++) {
			properties.setProperty(prefix+ErsCorrection.DP_DERIV_NAMES[i]+"_sel",           this.ilma_lma_select[i]+"");	
			properties.setProperty(prefix+ErsCorrection.DP_DERIV_NAMES[i]+"_regweight",     this.ilma_regularization_weights[i]+"");	
		}
		properties.setProperty(prefix+"ilma_lambda",            this.ilma_lambda+"");	
		properties.setProperty(prefix+"ilma_lambda_scale_good", this.ilma_lambda_scale_good+"");	
		properties.setProperty(prefix+"ilma_lambda_scale_bad",  this.ilma_lambda_scale_bad+"");	
		properties.setProperty(prefix+"ilma_lambda_max",        this.ilma_lambda_max+"");	
		properties.setProperty(prefix+"ilma_rms_diff",          this.ilma_rms_diff+"");	
		properties.setProperty(prefix+"ilma_num_iter",          this.ilma_num_iter+"");	
		properties.setProperty(prefix+"ilma_debug_level",       this.ilma_debug_level+"");	
	}
	public void getProperties(String prefix,Properties properties){
		for (int i = ErsCorrection.DP_DVAZ; i < ErsCorrection.DP_NUM_PARS; i++) {
			String pn_sel = prefix+ErsCorrection.DP_DERIV_NAMES[i]+"_sel";
			if (properties.getProperty(pn_sel)!=null) this.ilma_lma_select[i]=Boolean.parseBoolean(properties.getProperty(pn_sel));
			pn_sel = prefix+ErsCorrection.DP_DERIV_NAMES[i]+"_regweight";
			if (properties.getProperty(pn_sel)!=null) this.ilma_regularization_weights[i]=Double.parseDouble(properties.getProperty(pn_sel));
		}

		if (properties.getProperty(prefix+"ilma_lambda")!=null)            this.ilma_lambda=Double.parseDouble(properties.getProperty(prefix+"ilma_lambda"));
		if (properties.getProperty(prefix+"ilma_lambda_scale_good")!=null) this.ilma_lambda_scale_good=Double.parseDouble(properties.getProperty(prefix+"ilma_lambda_scale_good"));
		if (properties.getProperty(prefix+"ilma_lambda_scale_bad")!=null)  this.ilma_lambda_scale_bad=Double.parseDouble(properties.getProperty(prefix+"ilma_lambda_scale_bad"));
		if (properties.getProperty(prefix+"ilma_lambda_max")!=null)        this.ilma_lambda_max=Double.parseDouble(properties.getProperty(prefix+"ilma_lambda_max"));
		if (properties.getProperty(prefix+"ilma_rms_diff")!=null)          this.ilma_rms_diff=Double.parseDouble(properties.getProperty(prefix+"ilma_rms_diff"));
		if (properties.getProperty(prefix+"ilma_num_iter")!=null)          this.ilma_num_iter=Integer.parseInt(properties.getProperty(prefix+"ilma_num_iter"));
		if (properties.getProperty(prefix+"ilma_debug_level")!=null)       this.ilma_debug_level=Integer.parseInt(properties.getProperty(prefix+"ilma_debug_level"));
	}
	
	@Override
	public IntersceneLmaParameters clone() throws CloneNotSupportedException {
		IntersceneLmaParameters ilp =     new IntersceneLmaParameters();
		System.arraycopy(this.ilma_lma_select,             0, ilp.ilma_lma_select,             0, ilma_lma_select.length);
		System.arraycopy(this.ilma_regularization_weights, 0, ilp.ilma_regularization_weights, 0, ilma_regularization_weights.length);
		ilp.ilma_lambda =            this.ilma_lambda;
		ilp.ilma_lambda_scale_good = this.ilma_lambda_scale_good;
		ilp.ilma_lambda_scale_bad =  this.ilma_lambda_scale_bad;
		ilp.ilma_lambda_max =        this.ilma_lambda_max;
		ilp.ilma_rms_diff =          this.ilma_rms_diff;
		ilp.ilma_num_iter =          this.ilma_num_iter;
		ilp.ilma_debug_level =       this.ilma_debug_level;
		return ilp;
	}	
}
