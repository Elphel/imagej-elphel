package com.elphel.imagej.cameras;

import java.util.Properties;

import com.elphel.imagej.common.WindowTools;

import ij.gui.GenericDialog;

public class ColorProcParameters {
	public static final String AUX_PREFIX = "AUX-";
	public boolean lwir_islwir =        false;
	public double  lwir_low =           27000;
	public double  lwir_high    =       31000;
	public int     lwir_palette =           0; // 0 - white - hot, 1 - black - hot, 2+ - colored
	public boolean lwir_subtract_dc   = false;
	public boolean lwir_eq_chn   =      true;  // adjust average temperature between channels
	public boolean correct_vignetting = true;

	public double balanceRed;
	public double balanceBlue;
	public double gain;
	public double weightScaleR;
	public double weightScaleB;
	public double gamma;
	public double minLin;
	public double kr;
	public double kb;
	public double saturationRed;
	public double saturationBlue;
	public boolean useFirstY;
	public double maskSigma; // LPF for luma to calculate chroma mask
	public double maskMin; // minimal level for the mask (absolute luma values)
	public double maskMax; // maximal level for the mask (absolute luma values)
	public boolean combineWithSharpnessMask; // combine chroma mask with sharpness mask to reduce color leak through borders
	public double chromaBrightSigma; // LPF for chroma in the bright areas (determined by the mask)
	public double chromaDarkSigma;   // LPF for chroma in the dark areas (determined by the mask)
	public boolean corrBlueLeak; //Remove blue color leak in the darks near saturation
	public int    satDetSquareSize;    // Size of sliding square do detect potential overexposure
	public double satDetRelDiff;     // Most pixels in square should be within this difference from average
	public double satDetPartInside;  // Fraction of all pixels in the square to fit inside
	public double satDetMinFrac;     // minimal value for average compared to average over the whole picture
	public double satDetFinRelDiff;  // maximal difference from average for the saturated tile to be considered saturated
	public double satDetGrowRelDiff; // maximal difference from start tile average during growing of overexposed areas
	public double satDetNewWeight;   // weight of new pixel when expanding overexposed areas

	public int    satDetExpSym; // number of overexposure expand steps, not allowing any brighter
	public int    satDetExpOver;// number of overexposure expand steps, limited under, any over

	public int    satDetExpCleanUp; // number of overexposure expand steps, not allowing any brighter (final to clean up oscillations)
	public double satDetGrowRelDiffCleanUp; // maximal difference from start tile average during growing of overexposed areas (final to clean up oscillations)

	public int    blueOverShrink; // shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)
	public int    blueOverGrow; // grow blue overexposed area by this number of pixels
	public double blueBandWidth; // average amount of blue leak in pixels
	public double blueBandWidthDark; // average amount of blue leak in pixels (slope at dark)

	public double blueNeutral; // Value of Yb/Yrg ratio for the small areas where safe color c an not be found
	public double blueSolutionRadius;   // How far to trust blue color ratio from the found solution (in pixels)
	public boolean blueLeakNoHint;      // use blueNeutral in the small areas that do not have reliable color sample
	public boolean blueLeakNoBrighten; // Do not brighten corrected areas, only darken

	public boolean blueLeakFixWires; //Fix thin objects with saturated blue, but not R+G
	public double  blueLeakWiresSize; //size (in pixels) of the small objects to fix blue
	public double  blueLeakWiresThreshold; // relative to saturation level threshold of the small blue-flooded objects on red+green to process

	public boolean use8; // use 8 neighbors (false - only 4)

	public boolean isMonochrome() {
		return lwir_islwir; // for now it is the only reason to be monochrome
	}

	private ColorProcParameters() {}
	public ColorProcParameters(
			boolean lwir_islwir,        // false;
			double  lwir_low,           // 27000;
			double  lwir_high,          // 31000;
			int     lwir_palette,       //  0 - white - hot, 1 - black - hot, 2+ - colored
			boolean lwir_subtract_dc,   //   = false;
			boolean lwir_eq_chn,        // true
			boolean correct_vignetting, //  = true;

			double balanceRed,
			double balanceBlue,
			double gain,
			double weightScaleR,
			double weightScaleB,
			double gamma,
			double minLin,
			double kr,
			double kb,
			double saturationRed,
			double saturationBlue,
			boolean useFirstY,
			double maskSigma, // LPF for luma to calculate chroma mask
			double maskMin, // minimal level for the mask (absolute luma values)
			double maskMax, // maximal level for the mask (absolute luma values)
			boolean combineWithSharpnessMask, // combine chroma mask with sharpness mask to reduce color leak through borders
			double chromaBrightSigma, // LPF for chroma in the bright areas (determined by the mask)
			double chromaDarkSigma,   // LPF for chroma in the dark areas (determined by the mask)
			//----------------------Fixing blue near saturation ----------/
			// Saturation detect parameters
			boolean corrBlueLeak, //Remove blue color leak in the darks near saturation
			int satDetSquareSize,    // Size of sliding square do detect potential overexposure
			double satDetRelDiff,    // Most pixels in square should be within this difference from average
			double satDetPartInside, // Fraction of all pixels in the square to fit inside
			double satDetMinFrac,    // minimal value for average compared to average over the whole picture
			double satDetFinRelDiff, // maximal difference from average for the saturated tile to be considered saturated
			double satDetGrowRelDiff, //maximal difference from start tile average during growing of overexposed areas
			double satDetNewWeight,   // weight of new pixel when expanding overexposed areas

			int    satDetExpSym,     // number of overexposure expand steps, not allowing any brighter
			int    satDetExpOver,     // number of overexposure expand steps, limited under, any over
			int    satDetExpCleanUp, // number of overexposure expand steps, not allowing any brighter (final to clean up oscillations)
			double satDetGrowRelDiffCleanUp, // maximal difference from start tile average during growing of overexposed areas (final to clean up oscillations)
			int    blueOverShrink, // shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)
			int    blueOverGrow, // grow blue overexposed area by this number of pixels
			double blueBandWidth,
			double blueBandWidthDark, // average amount of blue leak in pixels (slope at dark)
			double blueNeutral, // Value of Yb/Yrg ratio for the small areas where safe color c an not be found
			double blueSolutionRadius, // How far to trust blue color ratio from the found solution (in pixels)
			boolean blueLeakNoHint,    // use blueNeutral in the small areas that do not have reliable color sample
			boolean blueLeakNoBrighten, // Do not brighten corrected areas, only darken
			boolean blueLeakFixWires, //Fix thin objects with saturated blue, but not R+G
			double  blueLeakWiresSize, //size (in pixels) of the small objects to fix blue
			double  blueLeakWiresThreshold, //size (in pixels) of the small objects to fix blue
			boolean use8 // use 8 neighbors (false - only 4)
			) {
		this.lwir_islwir = lwir_islwir;
		this.lwir_low = lwir_low;
		this.lwir_high = lwir_high;
		this.lwir_palette = lwir_palette;
		this.lwir_subtract_dc = lwir_subtract_dc;
		this.lwir_eq_chn =lwir_eq_chn;
		this.correct_vignetting = correct_vignetting;

		this.balanceRed = balanceRed;
		this.balanceBlue = balanceBlue;
		this.gain = gain;
		this.weightScaleR = weightScaleR;
		this.weightScaleB = weightScaleB;
		//  			this.sigma = sigma;
		this.gamma = gamma;
		this.minLin = minLin;
		this.kr = kr;
		this.kb = kb;
		this.saturationRed = saturationRed;
		this.saturationBlue = saturationBlue;
		this.useFirstY=useFirstY;
		this.maskSigma=maskSigma;
		this.maskMin=maskMin;
		this.maskMax=maskMax;
		this.combineWithSharpnessMask=combineWithSharpnessMask;
		this.chromaBrightSigma=chromaBrightSigma;
		this.chromaDarkSigma=chromaDarkSigma;
		this.corrBlueLeak=corrBlueLeak; //Remove blue color leak in the darks near saturation
		this.satDetSquareSize=satDetSquareSize;    // Size of sliding square do detect potential overexposure
		this.satDetRelDiff=satDetRelDiff;    // Most pixels in square should be within this difference from average
		this.satDetPartInside=satDetPartInside; // Fraction of all pixels in the square to fit inside
		this.satDetMinFrac=satDetMinFrac;    // minimal value for average compared to average over the whole picture
		this.satDetFinRelDiff=satDetFinRelDiff; // maximal difference from average for the saturated tile to be considered saturated
		this.satDetGrowRelDiff=satDetGrowRelDiff; //maximal difference from start tile average during growing of overexposed areas
		this.satDetNewWeight=satDetNewWeight;   // weight of new pixel when expanding overexposed areas
		this.satDetExpSym=satDetExpSym;     // number of overexposure expand steps, not allowing any brighter
		this.satDetExpOver=satDetExpOver;    // number of overexposure expand steps, limited under, any over
		this.satDetExpCleanUp=satDetExpCleanUp; // number of overexposure expand steps, not allowing any brighter (final to clean up oscillations)
		this.satDetGrowRelDiffCleanUp=satDetGrowRelDiffCleanUp; // maximal difference from start tile average during growing of overexposed areas (final to clean up oscillations)
		this.blueOverShrink=blueOverShrink; // shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)
		this.blueOverGrow=blueOverGrow; // grow blue overexposed area by this number of pixels
		this.blueBandWidth=blueBandWidth;
		this.blueBandWidthDark=blueBandWidthDark; // average amount of blue leak in pixels (slope at dark)
		this.blueNeutral=blueNeutral; // Value of Yb/Yrg ratio for the small areas where safe color c an not be found
		this.blueSolutionRadius=blueSolutionRadius; // How far to trust blue color ratio from the found solution (in pixels)
		this.blueLeakNoHint=blueLeakNoHint;    // use blueNeutral in the small areas that do not have reliable color sample
		this.blueLeakNoBrighten=blueLeakNoBrighten; // Do not brighten corrected areas, only darken

		this.blueLeakFixWires=blueLeakFixWires; //Fix thin objects with saturated blue, but not R+G
		this.blueLeakWiresSize=blueLeakWiresSize; //size (in pixels) of the small objects to fix blue
		this.blueLeakWiresThreshold=blueLeakWiresThreshold; //size (in pixels) of the small objects to fix blue

		this.use8=use8; // use 8 neighbors (false - only 4)
	}
	public void setProperties(String prefix,Properties properties){
		properties.setProperty(prefix+"lwir_islwir",        this.lwir_islwir+"");
		properties.setProperty(prefix+"lwir_low",           this.lwir_low+"");
		properties.setProperty(prefix+"lwir_high",          this.lwir_high+"");
		properties.setProperty(prefix+"lwir_palette",       this.lwir_palette+"");
		properties.setProperty(prefix+"lwir_subtract_dc",   this.lwir_subtract_dc+"");
		properties.setProperty(prefix+"lwir_eq_chn",        this.lwir_eq_chn+"");

		properties.setProperty(prefix+"correct_vignetting", this.correct_vignetting+"");

		properties.setProperty(prefix+"balanceRed",this.balanceRed+"");
		properties.setProperty(prefix+"balanceBlue",this.balanceBlue+"");
		properties.setProperty(prefix+"gain",this.gain+"");
		properties.setProperty(prefix+"weightScaleR",this.weightScaleR+"");
		properties.setProperty(prefix+"weightScaleB",this.weightScaleB+"");
		//  			properties.setProperty(prefix+"sigma",this.sigma+"");
		properties.setProperty(prefix+"gamma",this.gamma+"");
		properties.setProperty(prefix+"minLin",this.minLin+"");
		properties.setProperty(prefix+"kr",this.kr+"");
		properties.setProperty(prefix+"kb",this.kb+"");
		properties.setProperty(prefix+"saturationRed",this.saturationRed+"");
		properties.setProperty(prefix+"saturationBlue",this.saturationBlue+"");
		properties.setProperty(prefix+"useFirstY",this.useFirstY+"");
		properties.setProperty(prefix+"maskSigma",this.maskSigma+"");
		properties.setProperty(prefix+"maskMin",this.maskMin+"");
		properties.setProperty(prefix+"maskMax",this.maskMax+"");
		properties.setProperty(prefix+"combineWithSharpnessMask",this.combineWithSharpnessMask+"");
		properties.setProperty(prefix+"chromaBrightSigma",this.chromaBrightSigma+"");
		properties.setProperty(prefix+"chromaDarkSigma",this.chromaDarkSigma+"");
		properties.setProperty(prefix+"corrBlueLeak",this.corrBlueLeak+"");
		properties.setProperty(prefix+"satDetSquareSize",this.satDetSquareSize+"");
		properties.setProperty(prefix+"satDetRelDiff",this.satDetRelDiff+"");
		properties.setProperty(prefix+"satDetPartInside",this.satDetPartInside+"");
		properties.setProperty(prefix+"satDetMinFrac",this.satDetMinFrac+"");
		properties.setProperty(prefix+"satDetFinRelDiff",this.satDetFinRelDiff+"");
		properties.setProperty(prefix+"satDetGrowRelDiff",this.satDetGrowRelDiff+"");


		properties.setProperty(prefix+"satDetNewWeight",this.satDetNewWeight+"");

		properties.setProperty(prefix+"satDetExpSym",this.satDetExpSym+"");
		properties.setProperty(prefix+"satDetExpOver",this.satDetExpOver+"");

		properties.setProperty(prefix+"satDetExpCleanUp",this.satDetExpCleanUp+"");
		properties.setProperty(prefix+"satDetGrowRelDiffCleanUp",this.satDetGrowRelDiffCleanUp+"");
		properties.setProperty(prefix+"blueOverShrink",this.blueOverShrink+"");
		properties.setProperty(prefix+"blueOverGrow",this.blueOverGrow+"");
		properties.setProperty(prefix+"blueBandWidth",this.blueBandWidth+"");
		properties.setProperty(prefix+"blueBandWidthDark",this.blueBandWidthDark+"");

		properties.setProperty(prefix+"blueNeutral",this.blueNeutral+"");
		properties.setProperty(prefix+"blueSolutionRadius",this.blueSolutionRadius+"");

		properties.setProperty(prefix+"blueLeakNoHint",this.blueLeakNoHint+"");
		properties.setProperty(prefix+"blueLeakNoBrighten",this.blueLeakNoBrighten+"");

		properties.setProperty(prefix+"blueLeakFixWires",this.blueLeakFixWires+"");
		properties.setProperty(prefix+"blueLeakWiresSize",this.blueLeakWiresSize+"");
		properties.setProperty(prefix+"blueLeakWiresThreshold",this.blueLeakWiresThreshold+"");

		properties.setProperty(prefix+"use8",this.use8+"");

	}
	public boolean getProperties(String prefix,Properties properties){

		if (properties.getProperty(prefix+"lwir_islwir")!=null) this.lwir_islwir=Boolean.parseBoolean(properties.getProperty(prefix+"lwir_islwir"));
		if (properties.getProperty(prefix+"lwir_low")!=null) this.lwir_low=Double.parseDouble(properties.getProperty(prefix+"lwir_low"));
		if (properties.getProperty(prefix+"lwir_high")!=null) this.lwir_high=Double.parseDouble(properties.getProperty(prefix+"lwir_high"));
		if (properties.getProperty(prefix+"lwir_palette")!=null) this.lwir_palette=Integer.parseInt(properties.getProperty(prefix+"lwir_palette"));
		if (properties.getProperty(prefix+"lwir_subtract_dc")!=null) this.lwir_subtract_dc=Boolean.parseBoolean(properties.getProperty(prefix+"lwir_subtract_dc"));
		if (properties.getProperty(prefix+"lwir_eq_chn")!=null) this.lwir_eq_chn=Boolean.parseBoolean(properties.getProperty(prefix+"lwir_eq_chn"));
		if (properties.getProperty(prefix+"correct_vignetting")!=null) this.correct_vignetting=Boolean.parseBoolean(properties.getProperty(prefix+"correct_vignetting"));
		if (properties.getProperty(prefix+"balanceRed")!=null) this.balanceRed=Double.parseDouble(properties.getProperty(prefix+"balanceRed"));
		if (properties.getProperty(prefix+"balanceBlue")!=null) this.balanceBlue=Double.parseDouble(properties.getProperty(prefix+"balanceBlue"));
		if (properties.getProperty(prefix+"gain")!=null) this.gain=Double.parseDouble(properties.getProperty(prefix+"gain"));
		if (properties.getProperty(prefix+"weightScaleR")!=null) this.weightScaleR=Double.parseDouble(properties.getProperty(prefix+"weightScaleR"));
		if (properties.getProperty(prefix+"weightScaleB")!=null) this.weightScaleB=Double.parseDouble(properties.getProperty(prefix+"weightScaleB"));
		//  			this.sigma=Double.parseDouble(properties.getProperty(prefix+"sigma"));
		if (properties.getProperty(prefix+"gamma")!=null) this.gamma=Double.parseDouble(properties.getProperty(prefix+"gamma"));
		if (properties.getProperty(prefix+"minLin")!=null) this.minLin=Double.parseDouble(properties.getProperty(prefix+"minLin"));
		if (properties.getProperty(prefix+"kr")!=null) this.kr=Double.parseDouble(properties.getProperty(prefix+"kr"));
		if (properties.getProperty(prefix+"kb")!=null) this.kb=Double.parseDouble(properties.getProperty(prefix+"kb"));
		if (properties.getProperty(prefix+"saturationRed")!=null) this.saturationRed=Double.parseDouble(properties.getProperty(prefix+"saturationRed"));
		if (properties.getProperty(prefix+"saturationBlue")!=null) this.saturationBlue=Double.parseDouble(properties.getProperty(prefix+"saturationBlue"));
		if (properties.getProperty(prefix+"useFirstY")!=null) this.useFirstY=Boolean.parseBoolean(properties.getProperty(prefix+"useFirstY"));
		if (properties.getProperty(prefix+"maskSigma")!=null) this.maskSigma=Double.parseDouble(properties.getProperty(prefix+"maskSigma"));
		if (properties.getProperty(prefix+"maskMin")!=null) this.maskMin=Double.parseDouble(properties.getProperty(prefix+"maskMin"));
		if (properties.getProperty(prefix+"maskMax")!=null) this.maskMax=Double.parseDouble(properties.getProperty(prefix+"maskMax"));
		if (properties.getProperty(prefix+"combineWithSharpnessMask")!=null) this.combineWithSharpnessMask=Boolean.parseBoolean(properties.getProperty(prefix+"combineWithSharpnessMask"));
		if (properties.getProperty(prefix+"chromaBrightSigma")!=null) this.chromaBrightSigma=Double.parseDouble(properties.getProperty(prefix+"chromaBrightSigma"));
		if (properties.getProperty(prefix+"chromaDarkSigma")!=null) this.chromaDarkSigma=Double.parseDouble(properties.getProperty(prefix+"chromaDarkSigma"));

		if (properties.getProperty(prefix+"corrBlueLeak")!=null) this.corrBlueLeak=Boolean.parseBoolean(properties.getProperty(prefix+"corrBlueLeak"));
		if (properties.getProperty(prefix+"satDetSquareSize")!=null) this.satDetSquareSize=Integer.parseInt(properties.getProperty(prefix+"satDetSquareSize"));
		if (properties.getProperty(prefix+"satDetRelDiff")!=null)    this.satDetRelDiff=Double.parseDouble(properties.getProperty(prefix+"satDetRelDiff"));
		if (properties.getProperty(prefix+"satDetPartInside")!=null) this.satDetPartInside=Double.parseDouble(properties.getProperty(prefix+"satDetPartInside"));
		if (properties.getProperty(prefix+"satDetMinFrac")!=null)    this.satDetMinFrac=Double.parseDouble(properties.getProperty(prefix+"satDetMinFrac"));
		if (properties.getProperty(prefix+"satDetFinRelDiff")!=null) this.satDetFinRelDiff=Double.parseDouble(properties.getProperty(prefix+"satDetFinRelDiff"));
		if (properties.getProperty(prefix+"satDetGrowRelDiff")!=null) this.satDetGrowRelDiff=Double.parseDouble(properties.getProperty(prefix+"satDetGrowRelDiff"));
		if (properties.getProperty(prefix+"satDetNewWeight")!=null) this.satDetNewWeight=Double.parseDouble(properties.getProperty(prefix+"satDetNewWeight"));
		if (properties.getProperty(prefix+"satDetExpSym")!=null) this.satDetExpSym=Integer.parseInt(properties.getProperty(prefix+"satDetExpSym"));
		if (properties.getProperty(prefix+"satDetExpOver")!=null) this.satDetExpOver=Integer.parseInt(properties.getProperty(prefix+"satDetExpOver"));
		if (properties.getProperty(prefix+"satDetExpCleanUp")!=null) this.satDetExpCleanUp=Integer.parseInt(properties.getProperty(prefix+"satDetExpCleanUp"));
		if (properties.getProperty(prefix+"satDetGrowRelDiffCleanUp")!=null) this.satDetGrowRelDiffCleanUp=Double.parseDouble(properties.getProperty(prefix+"satDetGrowRelDiffCleanUp"));
		if (properties.getProperty(prefix+"blueOverShrink")!=null) this.blueOverShrink=Integer.parseInt(properties.getProperty(prefix+"blueOverShrink"));

		//  			if (properties.getProperty(prefix+"blueOverGrow")!=null) this.blueOverGrow=Integer.parseInt(properties.getProperty(prefix+"blueOverGrow"));
		if (properties.getProperty(prefix+"blueOverGrow")!=null) this.blueOverGrow=(int) Double.parseDouble(properties.getProperty(prefix+"blueOverGrow"));

		if (properties.getProperty(prefix+"blueBandWidth")!=null) this.blueBandWidth=Double.parseDouble(properties.getProperty(prefix+"blueBandWidth"));
		if (properties.getProperty(prefix+"blueBandWidthDark")!=null) this.blueBandWidthDark=Double.parseDouble(properties.getProperty(prefix+"blueBandWidthDark"));
		if (properties.getProperty(prefix+"blueNeutral")!=null) this.blueNeutral=Double.parseDouble(properties.getProperty(prefix+"blueNeutral"));
		if (properties.getProperty(prefix+"blueSolutionRadius")!=null) this.blueSolutionRadius=Double.parseDouble(properties.getProperty(prefix+"blueSolutionRadius"));
		if (properties.getProperty(prefix+"blueLeakNoHint")!=null) this.blueLeakNoHint=Boolean.parseBoolean(properties.getProperty(prefix+"blueLeakNoHint"));

		if (properties.getProperty(prefix+"blueLeakNoBrighten")!=null) this.blueLeakNoBrighten=Boolean.parseBoolean(properties.getProperty(prefix+"blueLeakNoBrighten"));

		if (properties.getProperty(prefix+"blueLeakFixWires")!=null) this.blueLeakFixWires=Boolean.parseBoolean(properties.getProperty(prefix+"blueLeakFixWires"));
		if (properties.getProperty(prefix+"blueLeakWiresSize")!=null) this.blueLeakWiresSize=Double.parseDouble(properties.getProperty(prefix+"blueLeakWiresSize"));
		if (properties.getProperty(prefix+"blueLeakWiresThreshold")!=null) this.blueLeakWiresThreshold=Double.parseDouble(properties.getProperty(prefix+"blueLeakWiresThreshold"));

		if (properties.getProperty(prefix+"use8")!=null) this.use8=Boolean.parseBoolean(properties.getProperty(prefix+"use8"));

		return (properties.getProperty(prefix+"gain")!=null);  // gain should always be present
	}

	public boolean showColorProcessDialog() {
		return showColorProcessDialog(null);
	}

	public boolean showColorProcessDialog(ColorProcParameters cp) {
		GenericDialog gd = new GenericDialog("Color processing parameters");
		if (cp != null) {
			gd.addCheckbox    ("Copy from master system, disregard other fields and reopen dialog", false);
		}
		gd.addMessage("--- Parameters related to thermal imaging ---");
		gd.addCheckbox    ("These sensors are thermal vision with absolute temperature",    this.lwir_islwir);
		gd.addNumericField("Lowest value (temperature) to display",                         this.lwir_low,     3); //0.53
		gd.addNumericField("Highest value (temperature) to display",                        this.lwir_high,     3); //0.53
		gd.addNumericField("LWIR pallet (0-white hot, 1-black hot, 2+ - pseudo colors ",    this.lwir_palette,    0);
		gd.addCheckbox    ("Subtract each image DC when conditioning",                      this.lwir_subtract_dc);
		gd.addCheckbox    ("Adjust average temperature between cameras",                    this.lwir_eq_chn);

		gd.addMessage("--- Common parameters and parameters, related to color imaging ---");
		gd.addCheckbox    ("Apply vignetting correctionUse ffrom sensor calibration files", this.correct_vignetting);
		gd.addNumericField("YCbCr gamma",                             this.gamma,     3); //0.53
		gd.addNumericField("Color saturation, red.green",             this.saturationRed,3); //2.0
		gd.addNumericField("Color saturation, blue/green",            this.saturationBlue,3); //2.0
		gd.addNumericField("Color balance, red-to-green",             this.balanceRed,     3); //1.8
		gd.addNumericField("Color balance, blue-to-green",            this.balanceBlue,    3); //1.8
		gd.addNumericField("Gain green",                              this.gain,      3); //1.8
		gd.addNumericField("Weight scale RED  (which color to use)",  this.weightScaleR,  3); //1.8
		gd.addNumericField("Weight scale BLUE (which color to use)",  this.weightScaleB,  3); //1.8
		gd.addNumericField("Minimal linear value to apply gammaY",    this.minLin,    3); //0.53
		gd.addNumericField("YCbCr Kb",                                this.kb,        3); //0.114
		gd.addNumericField("YCbCr Kr",                                this.kr,        3); //0.299
		gd.addCheckbox    ("Use first approximation for Y",           this.useFirstY);
		gd.addMessage("Color denoise parameters");
		gd.addNumericField("Low-pass sigma for luma to calculate mask",this.maskSigma,        3);
		gd.addNumericField("Filtered luma minimal level for mask transition",this.maskMin,        3);
		gd.addNumericField("Filtered luma maximal level for mask transition",this.maskMax,        3);
		gd.addCheckbox    ("Combine chroma mask with sharpness one (reduce color leak)",  this.combineWithSharpnessMask);
		gd.addNumericField("Low-pass sigma for chroma in bright areas",this.chromaBrightSigma,        3);
		gd.addNumericField("Low-pass sigma for chroma in dark areas",this.chromaDarkSigma,        3);
		gd.addMessage("Parameters to remove blue near overexposed areas");
		gd.addCheckbox ("Remove blue color leak in the darks near saturation",this.corrBlueLeak);

		gd.addNumericField("Size of sliding square do detect potential overexposure",this.satDetSquareSize,    0);
		gd.addNumericField("Most pixels in square should be within this difference from average",this.satDetRelDiff,       4);
		gd.addNumericField("Fraction of all pixels in the square to fit inside",this.satDetPartInside,    4);
		gd.addNumericField("Minimal value for average compared to average over the whole picture",this.satDetMinFrac,       4);
		gd.addNumericField("Maximal difference from average for the saturated tile to be considered saturated",this.satDetFinRelDiff,    4);
		gd.addNumericField("Maximal difference from the start tile average during growing",this.satDetGrowRelDiff,    4);


		gd.addNumericField("Weight of new pixel when expanding overexposed areas",this.satDetNewWeight,    4);


		gd.addNumericField("Number of overexposure expand steps, nothing special for brighter",this.satDetExpSym,    0);
		gd.addNumericField("Number of overexposure expand steps, limited under, any over",this.satDetExpOver,    0);

		gd.addNumericField("Maximal difference from the start tile average during growing (final to clean up oscillations)",this.satDetGrowRelDiffCleanUp,    4);
		gd.addNumericField("Number of overexposure expand steps, final to clean up oscillations)",this.satDetExpCleanUp,    0);
		gd.addNumericField("Shrink blue overexposed area by this number of pixels (to get to undisturbed R/G)",this.blueOverShrink,    0);
		gd.addNumericField("Grow blue overexposed area by this number of pixels",this.blueOverGrow,    0);
		gd.addNumericField("Average amount of blue leak in pixels",this.blueBandWidth,    4);
		gd.addNumericField("Average amount of blue leak in pixels (slope at dark)",this.blueBandWidthDark,    4);


		gd.addNumericField("Value of Yb/Yrg ratio for the small areas where safe color can not be found",this.blueNeutral,    4);
		gd.addNumericField("How far to trust blue color ratio from the found solution (in pixels)",this.blueSolutionRadius,    4);


		gd.addCheckbox    ("Use blueNeutral in the small areas that do not have reliable color sample",  this.blueLeakNoHint);
		gd.addCheckbox    ("Do not brighten blue in corrected areas, only darken",  this.blueLeakNoBrighten);

		gd.addCheckbox    ("Fix thin objects with saturated blue, but not R+G",  this.blueLeakFixWires);
		gd.addNumericField("Width (in pixels) of the small objects to fix blue flood",this.blueLeakWiresSize,    4);
		gd.addNumericField("Relative amplitude (fraction of saturation) of samll flooded by blue objects to process",this.blueLeakWiresThreshold,    4);

		gd.addCheckbox    ("Use 8 neighbors (false - only 4)",  this.use8);

		//  	    gd.addNumericField("Debug Level:",                          MASTER_DEBUG_LEVEL,      0);
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return false;
		if (cp != null) {
			boolean copy_from_master = gd.getNextBoolean();
			if (copy_from_master) {
				ColorProcParameters cp_backup = clone();
				copyFrom(cp);
				boolean copy_accepted = showColorProcessDialog(null);
				if (!copy_accepted) {
					copyFrom(cp_backup);
					return false;
				}
				return true;
			}
		}
		this.lwir_islwir=              gd.getNextBoolean();
		this.lwir_low=                 gd.getNextNumber();
		this.lwir_high=                gd.getNextNumber();
		this.lwir_palette=       (int) gd.getNextNumber();
		this.lwir_subtract_dc=         gd.getNextBoolean();
		this.lwir_eq_chn=              gd.getNextBoolean();
		this.correct_vignetting=       gd.getNextBoolean();

		this.gamma=                    gd.getNextNumber();
		this.saturationRed=            gd.getNextNumber();
		this.saturationBlue=           gd.getNextNumber();
		this.balanceRed=               gd.getNextNumber();
		this.balanceBlue=              gd.getNextNumber();
		this.gain=                     gd.getNextNumber();
		this.weightScaleR=             gd.getNextNumber();
		this.weightScaleB=             gd.getNextNumber();
		this.minLin=                   gd.getNextNumber();
		this.kb=                       gd.getNextNumber(); //---
		this.kr=                       gd.getNextNumber();
		this.useFirstY=                gd.getNextBoolean();
		this.maskSigma=                gd.getNextNumber();
		this.maskMin=                  gd.getNextNumber();
		this.maskMax=                  gd.getNextNumber();
		this.combineWithSharpnessMask= gd.getNextBoolean();
		this.chromaBrightSigma=        gd.getNextNumber();
		this.chromaDarkSigma=          gd.getNextNumber();
		this.corrBlueLeak=             gd.getNextBoolean();
		this.satDetSquareSize=   (int) gd.getNextNumber();
		this.satDetRelDiff=            gd.getNextNumber();
		this.satDetPartInside=         gd.getNextNumber();
		this.satDetMinFrac=            gd.getNextNumber();
		this.satDetFinRelDiff=         gd.getNextNumber();
		this.satDetGrowRelDiff=        gd.getNextNumber();
		this.satDetNewWeight=          gd.getNextNumber();

		this.satDetExpSym=       (int) gd.getNextNumber();
		this.satDetExpOver=      (int) gd.getNextNumber();

		this.satDetGrowRelDiffCleanUp= gd.getNextNumber();
		this.satDetExpCleanUp=   (int) gd.getNextNumber();
		this.blueOverShrink=     (int) gd.getNextNumber();
		this.blueOverGrow=       (int) gd.getNextNumber();
		this.blueBandWidth=            gd.getNextNumber();
		this.blueBandWidthDark=        gd.getNextNumber();
		this.blueNeutral=              gd.getNextNumber();
		this.blueSolutionRadius=       gd.getNextNumber();

		this.blueLeakNoHint=           gd.getNextBoolean();
		this.blueLeakNoBrighten=       gd.getNextBoolean();

		this.blueLeakFixWires=         gd.getNextBoolean();
		this.blueLeakWiresSize=        gd.getNextNumber();
		this.blueLeakWiresThreshold=   gd.getNextNumber();

		this.use8=                     gd.getNextBoolean();
		return true;
	}

	@Override
	public ColorProcParameters clone() {
		ColorProcParameters cp = new ColorProcParameters();
		cp.lwir_islwir =             this.lwir_islwir;
		cp.lwir_low =                this.lwir_low;
		cp.lwir_high =               this.lwir_high;
		cp.lwir_palette =            this.lwir_palette;
		cp.lwir_subtract_dc =        this.lwir_subtract_dc;
		cp.lwir_eq_chn =             this.lwir_eq_chn;
		cp.correct_vignetting =      this.correct_vignetting;

		cp.gamma=                    this.gamma;
		cp.saturationRed=            this.saturationRed;
		cp.saturationBlue=           this.saturationBlue;
		cp.balanceRed=               this.balanceRed;
		cp.balanceBlue=              this.balanceBlue;
		cp.gain=                     this.gain;
		cp.weightScaleR=             this.weightScaleR;
		cp.weightScaleB=             this.weightScaleB;
		cp.minLin=                   this.minLin;
		cp.kb=                       this.kb; //---
		cp.kr=                       this.kr;
		cp.useFirstY=                this.useFirstY;
		cp.maskSigma=                this.maskSigma;
		cp.maskMin=                  this.maskMin;
		cp.maskMax=                  this.maskMax;
		cp.combineWithSharpnessMask= this.combineWithSharpnessMask;
		cp.chromaBrightSigma=        this.chromaBrightSigma;
		cp.chromaDarkSigma=          this.chromaDarkSigma;
		cp.corrBlueLeak=             this.corrBlueLeak;
		cp.satDetSquareSize=         this.satDetSquareSize;
		cp.satDetRelDiff=            this.satDetRelDiff;
		cp.satDetPartInside=         this.satDetPartInside;
		cp.satDetMinFrac=            this.satDetMinFrac;
		cp.satDetFinRelDiff=         this.satDetFinRelDiff;
		cp.satDetGrowRelDiff=        this.satDetGrowRelDiff;
		cp.satDetNewWeight=          this.satDetNewWeight;

		cp.satDetExpSym=             this.satDetExpSym;
		cp.satDetExpOver=            this.satDetExpOver;

		cp.satDetGrowRelDiffCleanUp= this.satDetGrowRelDiffCleanUp;
		cp.satDetExpCleanUp=         this.satDetExpCleanUp;
		cp.blueOverShrink=           this.blueOverShrink;
		cp.blueOverGrow=             this.blueOverGrow;
		cp.blueBandWidth=            this.blueBandWidth;
		cp.blueBandWidthDark=        this.blueBandWidthDark;
		cp.blueNeutral=              this.blueNeutral;
		cp.blueSolutionRadius=       this.blueSolutionRadius;

		cp.blueLeakNoHint=           this.blueLeakNoHint;
		cp.blueLeakNoBrighten=       this.blueLeakNoBrighten;

		cp.blueLeakFixWires=         this.blueLeakFixWires;
		cp.blueLeakWiresSize=        this.blueLeakWiresSize;
		cp.blueLeakWiresThreshold=   this.blueLeakWiresThreshold;

		cp.use8=                     this.use8;
		return cp;
	}

	public void  copyFrom(ColorProcParameters cp) {
		this.lwir_islwir =             cp.lwir_islwir;
		this.lwir_low =                cp.lwir_low;
		this.lwir_high =               cp.lwir_high;
		this.lwir_palette =            cp.lwir_palette;
		this.lwir_subtract_dc =        cp.lwir_subtract_dc;
		this.lwir_eq_chn =             cp.lwir_eq_chn;
		this.correct_vignetting =      cp.correct_vignetting;

		this.gamma=                    cp.gamma;
		this.saturationRed=            cp.saturationRed;
		this.saturationBlue=           cp.saturationBlue;
		this.balanceRed=               cp.balanceRed;
		this.balanceBlue=              cp.balanceBlue;
		this.gain=                     cp.gain;
		this.weightScaleR=             cp.weightScaleR;
		this.weightScaleB=             cp.weightScaleB;
		this.minLin=                   cp.minLin;
		this.kb=                       cp.kb; //---
		this.kr=                       cp.kr;
		this.useFirstY=                cp.useFirstY;
		this.maskSigma=                cp.maskSigma;
		this.maskMin=                  cp.maskMin;
		this.maskMax=                  cp.maskMax;
		this.combineWithSharpnessMask= cp.combineWithSharpnessMask;
		this.chromaBrightSigma=        cp.chromaBrightSigma;
		this.chromaDarkSigma=          cp.chromaDarkSigma;
		this.corrBlueLeak=             cp.corrBlueLeak;
		this.satDetSquareSize=         cp.satDetSquareSize;
		this.satDetRelDiff=            cp.satDetRelDiff;
		this.satDetPartInside=         cp.satDetPartInside;
		this.satDetMinFrac=            cp.satDetMinFrac;
		this.satDetFinRelDiff=         cp.satDetFinRelDiff;
		this.satDetGrowRelDiff=        cp.satDetGrowRelDiff;
		this.satDetNewWeight=          cp.satDetNewWeight;
		this.satDetExpSym=             cp.satDetExpSym;
		this.satDetExpOver=            cp.satDetExpOver;
		this.satDetGrowRelDiffCleanUp= cp.satDetGrowRelDiffCleanUp;
		this.satDetExpCleanUp=         cp.satDetExpCleanUp;
		this.blueOverShrink=           cp.blueOverShrink;
		this.blueOverGrow=             cp.blueOverGrow;
		this.blueBandWidth=            cp.blueBandWidth;
		this.blueBandWidthDark=        cp.blueBandWidthDark;
		this.blueNeutral=              cp.blueNeutral;
		this.blueSolutionRadius=       cp.blueSolutionRadius;
		this.blueLeakNoHint=           cp.blueLeakNoHint;
		this.blueLeakNoBrighten=       cp.blueLeakNoBrighten;
		this.blueLeakFixWires=         cp.blueLeakFixWires;
		this.blueLeakWiresSize=        cp.blueLeakWiresSize;
		this.blueLeakWiresThreshold=   cp.blueLeakWiresThreshold;
		this.use8=                     cp.use8;
	}



}