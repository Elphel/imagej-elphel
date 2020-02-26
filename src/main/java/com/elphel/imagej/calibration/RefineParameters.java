/*
 **
 ** RefineParameters.java - Parameters for sensor residual correction
 **
 ** Copyright (C) 2011-2020 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Distortions.java is free software: you can redistribute it and/or modify
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

package com.elphel.imagej.calibration;

import java.util.Properties;

import com.elphel.imagej.common.WindowTools;

import ij.gui.GenericDialog;

public class RefineParameters{
	// New parameters 2020

	public double center_fract =       0.5;  // 0.5 of half-height
	public double transit_fract =      0.2;  // 0.2 of half-height - transition from center ortho to outer polar
	public double gaus_ang =           0.2;  // in radians
	public double gaus_rad =           0.05; // in fractions of the full radius
	public double max_diff_err_geom =  0.25; // before second pass linearly fade R/T and RGB where high-frequency error nears thios value
	public double max_diff_err_photo = 0.25;

	@Deprecated public boolean extrapolate=true;         // extrapolate sensor distortion correction
	@Deprecated public double  alphaThreshold =0.8;      // ignore sensor correction pixels with mask value below this
	@Deprecated public double  fatZero=0.01;             // when extrapolatging color transfer coefficients (flat field) use this for logariphm
	@Deprecated public double  extrapolationSigma=30.0;  // sigmna for Gaussian weight function when fittinga plane to known pixels
	// calculated for non-decimated pixels
	@Deprecated public double  extrapolationKSigma=2.0;  // consider pixels in 2*extrapolationSigma*extrapolationKSigma square when fitting
	@Deprecated public boolean smoothCorrection=true;    // apply Gaussian blur to calculated pixel correction field
	@Deprecated public double  smoothSigma=50.0; // sigma for Gaussian weight function when fittinga plane to known pixels
	public double  correctionScale=1.0; // scale correction when accumulating;
	public boolean showCumulativeCorrection=false; // show correction afther this one is applied
	public boolean showUnfilteredCorrection=true; // show this (additional) correction before extrapolation and/or smoothing
	public boolean showExtrapolationCorrection=false; // show Extrapolation
	public boolean showThisCorrection=false; // show this (additional) correction separately
	public boolean showPerImage=false;     // show residuals for each individual image
	public int     showIndividualNumber=0; // which image to show (-1 - all)
	public boolean applyCorrection=true;   // apply calculated corerction
	public boolean applyFlatField=true;   // apply calculated flat-field
	public boolean grid3DCorrection=true; // Correct patetrn grid node locations in 3d (false - in 2d only)
	public boolean rotateCorrection=true; // old value - did not yet understand why is it needed
	public double  grid3DMaximalZCorr=20.0; // Maximal Z-axis correction (if more will fall back to 2d correction algorithm)

	public boolean useVariations=  false; // allow different Z for different stations (for not a wall/stable pattern)
	public double  variationPenalty=0.001; // "stiffness" of individual (per-station) Z-values of the target pattern
	public boolean fixXY=          false; // adjust only Z of the target pattern, keep X and Y
	public boolean resetVariations=false;
	public boolean noFallBack=     true; // may have bugs - not tested yet


	@Deprecated public boolean usePatternAlpha= true;  // use pattern grid alpha data, false - old calculation

	// New individual parameters for modify pattern grid
	public boolean targetShowPerImage=false;
	public boolean targetShowThisCorrection=false;
	public boolean targetApplyCorrection=true;
	public double  targetCorrectionScale=1.0; // scale correction when accumulating;
	// New parameters for new sensor correction
	@Deprecated public boolean sensorExtrapolateDiff =      false; // true - extrapolate correction, false - composite
	@Deprecated public double sensorShrinkBlurComboSigma =  50.0;
	@Deprecated public double sensorShrinkBlurComboLevel =  0.25;
	@Deprecated public double sensorAlphaThreshold =        0.1;
	@Deprecated public double sensorStep =                  5;
	@Deprecated public double sensorInterpolationSigma=     100;
	@Deprecated public double sensorTangentialRadius=       0.5;
	@Deprecated public int    sensorScanDistance=           200;
	@Deprecated public int    sensorResultDistance=         500;
	@Deprecated public int    sensorInterpolationDegree=    2;
	//New parameters for Flat field correction
	public int    flatFieldSerNumber=          -1;
	public int    flatFieldReferenceStation=    0;
	public double flatFieldShrink=              100.0;
	public double flatFieldNonVignettedRadius = 1000.0;
	public double flatFieldMinimalAlpha =       0.01; // use %
	public double flatFieldMinimalContrast=     0.1;
	@Deprecated public double flatFieldMinimalAccumulate =  0.01; // use %
	public double flatFieldShrinkForMatching =  2.0;
	public double flatFieldMaxRelDiff =         0.1;  // use %
	public int    flatFieldShrinkMask=          2;
	public double flatFieldFadeBorder =         2.0;

	public boolean flatFieldResetMask=      true;
	public boolean flatFieldShowSensorMasks=false;
	public boolean flatFieldShowIndividual= false;
	public boolean flatFieldShowResult=     true;
	public boolean flatFieldApplyResult=    true;
	public boolean flatFieldUseInterpolate= true;
	public double  flatFieldMaskThresholdOcclusion=0.15; // use %
	public int     flatFieldShrinkOcclusion= 2;
	public double  flatFieldFadeOcclusion=   2.0;

	public boolean flatFieldIgnoreSensorFlatField= false;
	// Other
	public int     repeatFlatFieldSensor=10; // TODO: add stop !

	public double  specularHighPassSigma=            10.0;
	public double  specularLowPassSigma=              2.0;
	public double  specularDiffFromAverageThreshold= 0.01;
	public int     specularNumIter=                  5;
	public boolean specularApplyNewWeights=          true;
	public boolean specularPositiveDiffOnly=         true;
	public int     specularShowDebug=                1; // 0 - do not show, 1 - show on last iteration only, 2 - show always
	public RefineParameters refineParametersSmall; // same parameters for a small sensor
	public boolean is_small = false;

// will keep @Deprecated parameters, but remove them from dialogs

	public RefineParameters(){}
	public RefineParameters(
		    //new in 2020
			double center_fract,
			double transit_fract,
			double gaus_ang,
			double gaus_rad,
			double max_diff_err_geom,
			double max_diff_err_photo,

			boolean extrapolate,
			double alphaThreshold,
			double fatZero,
			double extrapolationSigma,
			double extrapolationKSigma,
			boolean smoothCorrection,
			double smoothSigma,
			double correctionScale,
			boolean showCumulativeCorrection,
			boolean showUnfilteredCorrection,
			boolean showExtrapolationCorrection,
			boolean showThisCorrection,
			boolean showPerImage,
			int     showIndividualNumber, // which image to show (-1 - all)
			boolean applyCorrection,
			boolean applyFlatField,   // apply calculated flat-field
			boolean grid3DCorrection, // Correct patetrn grid node locations in 3d (false - in 2d only)
			boolean rotateCorrection, // not clear
			double  grid3DMaximalZCorr, // Maximal Z-axis correc tion (if more will fall back to 2d correction algorithm)
			boolean useVariations,
			double  variationPenalty, // "stiffness" of individual (per-station) Z-values of the target pattern
			boolean  fixXY,
			boolean  resetVariations,
			boolean  noFallBack, // may have bugs - not tested yet
			boolean usePatternAlpha,
			boolean  targetShowPerImage,
			boolean  targetShowThisCorrection,
			boolean  targetApplyCorrection,
			double   targetCorrectionScale,
			boolean sensorExtrapolateDiff,
			double sensorShrinkBlurComboSigma,
			double sensorShrinkBlurComboLevel,
			double sensorAlphaThreshold,
			double sensorStep,
			double sensorInterpolationSigma,
			double sensorTangentialRadius,
			int    sensorScanDistance,
			int    sensorResultDistance,
			int    sensorInterpolationDegree,
			int flatFieldSerNumber,
			int flatFieldReferenceStation,
			double flatFieldShrink,
			double flatFieldNonVignettedRadius,
			double flatFieldMinimalAlpha,
			double flatFieldMinimalContrast,
			double flatFieldMinimalAccumulate,
			double flatFieldShrinkForMatching,
			double flatFieldMaxRelDiff,
			int    flatFieldShrinkMask,
			double flatFieldFadeBorder,
			boolean flatFieldResetMask,
			boolean flatFieldShowSensorMasks,
			boolean flatFieldShowIndividual,
			boolean flatFieldShowResult,
			boolean flatFieldApplyResult,
			boolean flatFieldUseInterpolate,
			double  flatFieldMaskThresholdOcclusion,
			int     flatFieldShrinkOcclusion,
			double  flatFieldFadeOcclusion,
			boolean flatFieldIgnoreSensorFlatField,
			int    repeatFlatFieldSensor,
			double  specularHighPassSigma,
			double  specularLowPassSigma,
			double  specularDiffFromAverageThreshold,
			int     specularNumIter,
			boolean specularApplyNewWeights,
			boolean specularPositiveDiffOnly,
			int     specularShowDebug,
			RefineParameters refineParametersSmall, // same parameters for a small sensor
			boolean is_small
			){
	    //new in 2020
		this.center_fract =  center_fract;
		this.transit_fract = transit_fract;
		this.gaus_ang =      gaus_ang;
		this.gaus_rad =      gaus_rad;
		this.max_diff_err_geom = max_diff_err_geom;
		this.max_diff_err_photo = max_diff_err_photo;

		this.extrapolate=extrapolate;
		this.alphaThreshold=alphaThreshold;
		this.fatZero=fatZero;        // when extrapolatging color transfer coefficients (flat field) use this for logariphm
		this.extrapolationSigma=extrapolationSigma;
		this.extrapolationKSigma=extrapolationKSigma;
		this.smoothCorrection=smoothCorrection;
		this.smoothSigma=smoothSigma;
		this.correctionScale=correctionScale;
		this.showCumulativeCorrection=showCumulativeCorrection;
		this.showUnfilteredCorrection=showUnfilteredCorrection;
		this.showExtrapolationCorrection=showExtrapolationCorrection;
		this.showThisCorrection=showThisCorrection;
		this.showPerImage=showPerImage;
		this.showIndividualNumber=showIndividualNumber; // which image to show (-1 - all)
		this.applyCorrection=applyCorrection;
		this.applyFlatField=applyFlatField;
		this.grid3DCorrection=grid3DCorrection;
		this.rotateCorrection=rotateCorrection; // not clear
		this.grid3DMaximalZCorr=grid3DMaximalZCorr; // Maximal Z-axis correc tion (if more will fall back to 2d correction algorithm)
		this.useVariations=useVariations;
		this.variationPenalty=variationPenalty; // "stiffness" of individual (per-station) Z-values of the target pattern
		this.fixXY=fixXY;
		this.resetVariations=resetVariations;
		this.noFallBack=     noFallBack; // may have bugs - not tested yet
		this.usePatternAlpha=usePatternAlpha;
		this.targetShowPerImage=       targetShowPerImage;
		this.targetShowThisCorrection=  targetShowThisCorrection;
		this.targetApplyCorrection=     targetApplyCorrection;
		this.targetCorrectionScale=     targetCorrectionScale;
		this.sensorExtrapolateDiff=     sensorExtrapolateDiff;
		this.sensorShrinkBlurComboSigma=sensorShrinkBlurComboSigma;
		this.sensorShrinkBlurComboLevel=sensorShrinkBlurComboLevel;
		this.sensorAlphaThreshold=sensorAlphaThreshold;
		this.sensorStep=sensorStep;
		this.sensorInterpolationSigma=sensorInterpolationSigma;
		this.sensorTangentialRadius=sensorTangentialRadius;
		this.sensorScanDistance=sensorScanDistance;
		this.sensorResultDistance=sensorResultDistance;
		this.sensorInterpolationDegree=sensorInterpolationDegree;
		this.flatFieldSerNumber=flatFieldSerNumber;
		this.flatFieldReferenceStation=flatFieldReferenceStation;
		this.flatFieldShrink=flatFieldShrink;
		this.flatFieldNonVignettedRadius=flatFieldNonVignettedRadius;
		this.flatFieldMinimalAlpha=flatFieldMinimalAlpha;
		this.flatFieldMinimalContrast=flatFieldMinimalContrast;
		this.flatFieldMinimalAccumulate=flatFieldMinimalAccumulate;
		this.flatFieldShrinkForMatching=flatFieldShrinkForMatching;
		this.flatFieldMaxRelDiff=flatFieldMaxRelDiff;
		this.flatFieldShrinkMask=flatFieldShrinkMask;
		this.flatFieldFadeBorder=flatFieldFadeBorder;
		this.flatFieldResetMask=flatFieldResetMask;
		this.flatFieldShowSensorMasks=flatFieldShowSensorMasks;
		this.flatFieldShowIndividual=flatFieldShowIndividual;
		this.flatFieldShowResult=flatFieldShowResult;
		this.flatFieldApplyResult=flatFieldApplyResult;
		this.flatFieldUseInterpolate=flatFieldUseInterpolate;
		this.flatFieldMaskThresholdOcclusion=flatFieldMaskThresholdOcclusion;
		this.flatFieldShrinkOcclusion=flatFieldShrinkOcclusion;
		this.flatFieldFadeOcclusion=flatFieldFadeOcclusion;
		this.flatFieldIgnoreSensorFlatField=flatFieldIgnoreSensorFlatField;
		//         	  this.flatFieldUseSelectedChannels=flatFieldUseSelectedChannels;
		this.repeatFlatFieldSensor=repeatFlatFieldSensor;

		this.specularHighPassSigma=specularHighPassSigma;
		this.specularLowPassSigma=specularLowPassSigma;
		this.specularDiffFromAverageThreshold=specularDiffFromAverageThreshold;
		this.specularNumIter=specularNumIter;
		this.specularApplyNewWeights=specularApplyNewWeights;
		this.specularPositiveDiffOnly=specularPositiveDiffOnly;
		this.specularShowDebug=specularShowDebug;
		this.refineParametersSmall=refineParametersSmall; // same parameters for a small sensor
		this.is_small = is_small;
	}
	@Override
	public RefineParameters clone(){
		return new RefineParameters(
			    //new in 2020
			    this.center_fract,
			    this.transit_fract,
			    this.gaus_ang,
			    this.gaus_rad,
			    this.max_diff_err_geom,
			    this.max_diff_err_photo,

				this.extrapolate,
				this.alphaThreshold,
				this.fatZero,
				this.extrapolationSigma,
				this.extrapolationKSigma,
				this.smoothCorrection,
				this.smoothSigma,
				this.correctionScale,
				this.showCumulativeCorrection,
				this.showUnfilteredCorrection,
				this.showExtrapolationCorrection,
				this.showThisCorrection,
				this.showPerImage,
				this.showIndividualNumber,
				this.applyCorrection,
				this.applyFlatField,
				this.grid3DCorrection,
				this.rotateCorrection, // not clear
				this.grid3DMaximalZCorr, // Maximal Z-axis correc tion (if more will fall back to 2d correction algorithm)
				this.useVariations,
				this.variationPenalty, // "stiffness" of individual (per-station) Z-values of the target pattern
				this.fixXY,
				this.resetVariations,
				this.noFallBack, // may have bugs - not tested yet
				this.usePatternAlpha,
				this.targetShowPerImage,
				this.targetShowThisCorrection,
				this.targetApplyCorrection,
				this.targetCorrectionScale,
				this.sensorExtrapolateDiff,
				this.sensorShrinkBlurComboSigma,
				this.sensorShrinkBlurComboLevel,
				this.sensorAlphaThreshold,
				this.sensorStep,
				this.sensorInterpolationSigma,
				this.sensorTangentialRadius,
				this.sensorScanDistance,
				this.sensorResultDistance,
				this.sensorInterpolationDegree,
				this.flatFieldSerNumber,
				this.flatFieldReferenceStation,
				this.flatFieldShrink,
				this.flatFieldNonVignettedRadius,
				this.flatFieldMinimalAlpha,
				this.flatFieldMinimalContrast,
				this.flatFieldMinimalAccumulate,
				this.flatFieldShrinkForMatching,
				this.flatFieldMaxRelDiff,
				this.flatFieldShrinkMask,
				this.flatFieldFadeBorder,
				this.flatFieldResetMask,
				this.flatFieldShowSensorMasks,
				this.flatFieldShowIndividual,
				this.flatFieldShowResult,
				this.flatFieldApplyResult,
				this.flatFieldUseInterpolate,
				this.flatFieldMaskThresholdOcclusion,
				this.flatFieldShrinkOcclusion,
				this.flatFieldFadeOcclusion,
				this.flatFieldIgnoreSensorFlatField,
				this.repeatFlatFieldSensor,
				this.specularHighPassSigma,
				this.specularLowPassSigma,
				this.specularDiffFromAverageThreshold,
				this.specularNumIter,
				this.specularApplyNewWeights,
				this.specularPositiveDiffOnly,
				this.specularShowDebug,
				this.refineParametersSmall,
				this.is_small);
	}
	public void setProperties(String prefix,Properties properties){
		//new in 2020
		properties.setProperty(prefix+"center_fract",      this.center_fract+"");
		properties.setProperty(prefix+"transit_fract",     this.transit_fract+"");
		properties.setProperty(prefix+"gaus_ang",          this.gaus_ang+"");
		properties.setProperty(prefix+"gaus_rad",          this.gaus_rad+"");
		properties.setProperty(prefix+"max_diff_err_geom", this.max_diff_err_geom+"");
		properties.setProperty(prefix+"max_diff_err_photo",this.max_diff_err_photo+"");

		properties.setProperty(prefix+"extrapolate",this.extrapolate+"");
		properties.setProperty(prefix+"alphaThreshold",this.alphaThreshold+"");
		properties.setProperty(prefix+"fatZero",this.fatZero+"");
		properties.setProperty(prefix+"extrapolationSigma",this.extrapolationSigma+"");
		properties.setProperty(prefix+"extrapolationKSigma",this.extrapolationKSigma+"");
		properties.setProperty(prefix+"smoothCorrection",this.smoothCorrection+"");
		properties.setProperty(prefix+"smoothSigma",this.smoothSigma+"");
		properties.setProperty(prefix+"correctionScale",this.correctionScale+"");
		properties.setProperty(prefix+"showCumulativeCorrection",this.showCumulativeCorrection+"");
		properties.setProperty(prefix+"showUnfilteredCorrection",this.showUnfilteredCorrection+"");
		properties.setProperty(prefix+"showExtrapolationCorrection",this.showExtrapolationCorrection+"");
		properties.setProperty(prefix+"showThisCorrection",this.showThisCorrection+"");
		properties.setProperty(prefix+"showPerImage",this.showPerImage+"");
		properties.setProperty(prefix+"showIndividualNumber",this.showIndividualNumber+"");
		properties.setProperty(prefix+"applyCorrection",this.applyCorrection+"");
		properties.setProperty(prefix+"applyFlatField",this.applyFlatField+"");
		properties.setProperty(prefix+"grid3DCorrection",this.grid3DCorrection+"");
		properties.setProperty(prefix+"rotateCorrection",this.rotateCorrection+"");
		properties.setProperty(prefix+"grid3DMaximalZCorr",this.grid3DMaximalZCorr+"");

		properties.setProperty(prefix+"useVariations",this.useVariations+"");
		properties.setProperty(prefix+"variationPenalty",this.variationPenalty+"");
		properties.setProperty(prefix+"fixXY",this.fixXY+"");

		properties.setProperty(prefix+"resetVariations",this.resetVariations+"");
		properties.setProperty(prefix+"noFallBack",this.noFallBack+"");

		properties.setProperty(prefix+"usePatternAlpha",this.usePatternAlpha+"");

		properties.setProperty(prefix+"targetShowPerImage",this.targetShowPerImage+"");
		properties.setProperty(prefix+"targetShowThisCorrection",this.targetShowThisCorrection+"");
		properties.setProperty(prefix+"targetApplyCorrection",this.targetApplyCorrection+"");
		properties.setProperty(prefix+"targetCorrectionScale",this.targetCorrectionScale+"");

		properties.setProperty(prefix+"sensorExtrapolateDiff",this.sensorExtrapolateDiff+"");
		properties.setProperty(prefix+"sensorShrinkBlurComboSigma",this.sensorShrinkBlurComboSigma+"");

		properties.setProperty(prefix+"sensorShrinkBlurComboLevel",this.sensorShrinkBlurComboLevel+"");
		properties.setProperty(prefix+"sensorAlphaThreshold",this.sensorAlphaThreshold+"");
		properties.setProperty(prefix+"sensorStep",this.sensorStep+"");
		properties.setProperty(prefix+"sensorInterpolationSigma",this.sensorInterpolationSigma+"");
		properties.setProperty(prefix+"sensorTangentialRadius",this.sensorTangentialRadius+"");

		properties.setProperty(prefix+"sensorScanDistance",this.sensorScanDistance+"");
		properties.setProperty(prefix+"sensorResultDistance",this.sensorResultDistance+"");
		properties.setProperty(prefix+"sensorInterpolationDegree",this.sensorInterpolationDegree+"");
		properties.setProperty(prefix+"flatFieldSerNumber",this.flatFieldSerNumber+"");
		properties.setProperty(prefix+"flatFieldReferenceStation",this.flatFieldReferenceStation+"");

		properties.setProperty(prefix+"flatFieldShrink",this.flatFieldShrink+"");
		properties.setProperty(prefix+"flatFieldNonVignettedRadius",this.flatFieldNonVignettedRadius+"");
		properties.setProperty(prefix+"flatFieldMinimalAlpha",this.flatFieldMinimalAlpha+"");

		properties.setProperty(prefix+"flatFieldMinimalContrast",this.flatFieldMinimalContrast+"");
		properties.setProperty(prefix+"flatFieldMinimalAccumulate",this.flatFieldMinimalAccumulate+"");
		properties.setProperty(prefix+"flatFieldShrinkForMatching",this.flatFieldShrinkForMatching+"");

		properties.setProperty(prefix+"flatFieldMaxRelDiff",this.flatFieldMaxRelDiff+"");
		properties.setProperty(prefix+"flatFieldShrinkMask",this.flatFieldShrinkMask+"");
		properties.setProperty(prefix+"flatFieldFadeBorder",this.flatFieldFadeBorder+"");
		properties.setProperty(prefix+"flatFieldResetMask",this.flatFieldResetMask+"");
		properties.setProperty(prefix+"flatFieldShowSensorMasks",this.flatFieldShowSensorMasks+"");

		properties.setProperty(prefix+"flatFieldShowIndividual",this.flatFieldShowIndividual+"");
		properties.setProperty(prefix+"flatFieldShowResult",this.flatFieldShowResult+"");
		properties.setProperty(prefix+"flatFieldApplyResult",this.flatFieldApplyResult+"");
		properties.setProperty(prefix+"flatFieldUseInterpolate",this.flatFieldUseInterpolate+"");
		properties.setProperty(prefix+"flatFieldMaskThresholdOcclusion",this.flatFieldMaskThresholdOcclusion+"");

		properties.setProperty(prefix+"flatFieldShrinkOcclusion",this.flatFieldShrinkOcclusion+"");
		properties.setProperty(prefix+"flatFieldFadeOcclusion",this.flatFieldFadeOcclusion+"");
		properties.setProperty(prefix+"flatFieldIgnoreSensorFlatField",this.flatFieldIgnoreSensorFlatField+"");
		properties.setProperty(prefix+"repeatFlatFieldSensor",this.repeatFlatFieldSensor+"");

		properties.setProperty(prefix+"specularHighPassSigma",this.specularHighPassSigma+"");
		properties.setProperty(prefix+"specularLowPassSigma", this.specularLowPassSigma+"");
		properties.setProperty(prefix+"specularDiffFromAverageThreshold",this.specularDiffFromAverageThreshold+"");
		properties.setProperty(prefix+"specularNumIter",this.specularNumIter+"");
		properties.setProperty(prefix+"specularApplyNewWeights",this.specularApplyNewWeights+"");
		properties.setProperty(prefix+"specularPositiveDiffOnly",this.specularPositiveDiffOnly+"");
		properties.setProperty(prefix+"specularShowDebug",this.specularShowDebug+"");
		properties.setProperty(prefix+"is_small",this.is_small+"");
		if (this.refineParametersSmall != null) {
			/*
			if (this.refineParametersSmall.refineParametersSmall != null) {
				System.out.println("Multiple-recursive RefineParameters - should only be main and small");
				System.out.println("*********** This is a bug! ***************");
				this.refineParametersSmall.refineParametersSmall = null;
				return;
			}
			*/
			this.refineParametersSmall.setProperties(prefix+"SMALL.", properties);
		}
	}

	public void getProperties(String prefix,Properties properties){
		//new in 2020
		if (properties.getProperty(prefix+"center_fract")!=null)       this.center_fract=Double.parseDouble(properties.getProperty(prefix+"center_fract"));
		if (properties.getProperty(prefix+"transit_fract")!=null)      this.transit_fract=Double.parseDouble(properties.getProperty(prefix+"transit_fract"));
		if (properties.getProperty(prefix+"gaus_ang")!=null)           this.gaus_ang=Double.parseDouble(properties.getProperty(prefix+"gaus_ang"));
		if (properties.getProperty(prefix+"gaus_rad")!=null)           this.gaus_rad=Double.parseDouble(properties.getProperty(prefix+"gaus_rad"));
		if (properties.getProperty(prefix+"max_diff_err_geom")!=null)  this.max_diff_err_geom=Double.parseDouble(properties.getProperty(prefix+"max_diff_err_geom"));
		if (properties.getProperty(prefix+"max_diff_err_photo")!=null) this.max_diff_err_photo=Double.parseDouble(properties.getProperty(prefix+"max_diff_err_photo"));

		if (properties.getProperty(prefix+"extrapolate")!=null)
			this.extrapolate=Boolean.parseBoolean(properties.getProperty(prefix+"extrapolate"));
		if (properties.getProperty(prefix+"alphaThreshold")!=null)
			this.alphaThreshold=Double.parseDouble(properties.getProperty(prefix+"alphaThreshold"));
		if (properties.getProperty(prefix+"fatZero")!=null)
			this.fatZero=Double.parseDouble(properties.getProperty(prefix+"fatZero"));
		if (properties.getProperty(prefix+"extrapolationSigma")!=null)
			this.extrapolationSigma=Double.parseDouble(properties.getProperty(prefix+"extrapolationSigma"));
		if (properties.getProperty(prefix+"extrapolationKSigma")!=null)
			this.extrapolationKSigma=Double.parseDouble(properties.getProperty(prefix+"extrapolationKSigma"));
		if (properties.getProperty(prefix+"smoothCorrection")!=null)
			this.smoothCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"smoothCorrection"));
		if (properties.getProperty(prefix+"smoothSigma")!=null)
			this.smoothSigma=Double.parseDouble(properties.getProperty(prefix+"smoothSigma"));
		if (properties.getProperty(prefix+"correctionScale")!=null)
			this.correctionScale=Double.parseDouble(properties.getProperty(prefix+"correctionScale"));
		if (properties.getProperty(prefix+"showCumulativeCorrection")!=null)
			this.showCumulativeCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"showCumulativeCorrection"));
		if (properties.getProperty(prefix+"showUnfilteredCorrection")!=null)
			this.showUnfilteredCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"showUnfilteredCorrection"));
		if (properties.getProperty(prefix+"showExtrapolationCorrection")!=null)
			this.showExtrapolationCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"showExtrapolationCorrection"));
		if (properties.getProperty(prefix+"showThisCorrection")!=null)
			this.showThisCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"showThisCorrection"));
		if (properties.getProperty(prefix+"showPerImage")!=null)
			this.showPerImage=Boolean.parseBoolean(properties.getProperty(prefix+"showPerImage"));
		if (properties.getProperty(prefix+"showIndividualNumber")!=null)
			this.showIndividualNumber=Integer.parseInt(properties.getProperty(prefix+"showIndividualNumber"));
		if (properties.getProperty(prefix+"applyCorrection")!=null)
			this.applyCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"applyCorrection"));
		if (properties.getProperty(prefix+"applyFlatField")!=null)
			this.applyFlatField=Boolean.parseBoolean(properties.getProperty(prefix+"applyFlatField"));
		if (properties.getProperty(prefix+"grid3DCorrection")!=null)
			this.grid3DCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"grid3DCorrection"));
		if (properties.getProperty(prefix+"rotateCorrection")!=null)
			this.rotateCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"rotateCorrection"));
		if (properties.getProperty(prefix+"grid3DMaximalZCorr")!=null)
			this.grid3DMaximalZCorr=Double.parseDouble(properties.getProperty(prefix+"grid3DMaximalZCorr"));

		if (properties.getProperty(prefix+"useVariations")!=null)
			this.useVariations=Boolean.parseBoolean(properties.getProperty(prefix+"useVariations"));
		if (properties.getProperty(prefix+"variationPenalty")!=null)
			this.variationPenalty=Double.parseDouble(properties.getProperty(prefix+"variationPenalty"));
		if (properties.getProperty(prefix+"fixXY")!=null)
			this.fixXY=Boolean.parseBoolean(properties.getProperty(prefix+"fixXY"));
		if (properties.getProperty(prefix+"resetVariations")!=null)
			this.resetVariations=Boolean.parseBoolean(properties.getProperty(prefix+"resetVariations"));
		if (properties.getProperty(prefix+"noFallBack")!=null)
			this.noFallBack=Boolean.parseBoolean(properties.getProperty(prefix+"noFallBack"));
		if (properties.getProperty(prefix+"usePatternAlpha")!=null)
			this.usePatternAlpha=Boolean.parseBoolean(properties.getProperty(prefix+"usePatternAlpha"));
		if (properties.getProperty(prefix+"targetShowPerImage")!=null)
			this.targetShowPerImage=Boolean.parseBoolean(properties.getProperty(prefix+"targetShowPerImage"));
		if (properties.getProperty(prefix+"targetShowThisCorrection")!=null)
			this.targetShowThisCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"targetShowThisCorrection"));
		if (properties.getProperty(prefix+"targetApplyCorrection")!=null)
			this.targetApplyCorrection=Boolean.parseBoolean(properties.getProperty(prefix+"targetApplyCorrection"));
		if (properties.getProperty(prefix+"targetCorrectionScale")!=null)
			this.targetCorrectionScale=Double.parseDouble(properties.getProperty(prefix+"targetCorrectionScale"));
		if (properties.getProperty(prefix+"sensorExtrapolateDiff")!=null)
			this.sensorExtrapolateDiff=Boolean.parseBoolean(properties.getProperty(prefix+"sensorExtrapolateDiff"));
		if (properties.getProperty(prefix+"sensorShrinkBlurComboSigma")!=null)
			this.sensorShrinkBlurComboSigma=Double.parseDouble(properties.getProperty(prefix+"sensorShrinkBlurComboSigma"));
		if (properties.getProperty(prefix+"sensorShrinkBlurComboLevel")!=null)
			this.sensorShrinkBlurComboLevel=Double.parseDouble(properties.getProperty(prefix+"sensorShrinkBlurComboLevel"));
		if (properties.getProperty(prefix+"sensorAlphaThreshold")!=null)
			this.sensorAlphaThreshold=Double.parseDouble(properties.getProperty(prefix+"sensorAlphaThreshold"));
		if (properties.getProperty(prefix+"sensorStep")!=null)
			this.sensorStep=Double.parseDouble(properties.getProperty(prefix+"sensorStep"));
		if (properties.getProperty(prefix+"sensorInterpolationSigma")!=null)
			this.sensorInterpolationSigma=Double.parseDouble(properties.getProperty(prefix+"sensorInterpolationSigma"));
		if (properties.getProperty(prefix+"sensorTangentialRadius")!=null)
			this.sensorTangentialRadius=Double.parseDouble(properties.getProperty(prefix+"sensorTangentialRadius"));
		if (properties.getProperty(prefix+"sensorScanDistance")!=null)
			this.sensorScanDistance=Integer.parseInt(properties.getProperty(prefix+"sensorScanDistance"));
		if (properties.getProperty(prefix+"sensorResultDistance")!=null)
			this.sensorResultDistance=Integer.parseInt(properties.getProperty(prefix+"sensorResultDistance"));
		if (properties.getProperty(prefix+"sensorInterpolationDegree")!=null)
			this.sensorInterpolationDegree=Integer.parseInt(properties.getProperty(prefix+"sensorInterpolationDegree"));
		if (properties.getProperty(prefix+"flatFieldSerNumber")!=null)
			this.flatFieldSerNumber=Integer.parseInt(properties.getProperty(prefix+"flatFieldSerNumber"));
		if (properties.getProperty(prefix+"flatFieldReferenceStation")!=null)
			this.flatFieldReferenceStation=Integer.parseInt(properties.getProperty(prefix+"flatFieldReferenceStation"));
		if (properties.getProperty(prefix+"flatFieldShrink")!=null)
			this.flatFieldShrink=Double.parseDouble(properties.getProperty(prefix+"flatFieldShrink"));
		if (properties.getProperty(prefix+"flatFieldNonVignettedRadius")!=null)
			this.flatFieldNonVignettedRadius=Double.parseDouble(properties.getProperty(prefix+"flatFieldNonVignettedRadius"));
		if (properties.getProperty(prefix+"flatFieldMinimalAlpha")!=null)
			this.flatFieldMinimalAlpha=Double.parseDouble(properties.getProperty(prefix+"flatFieldMinimalAlpha"));
		if (properties.getProperty(prefix+"flatFieldMinimalContrast")!=null)
			this.flatFieldMinimalContrast=Double.parseDouble(properties.getProperty(prefix+"flatFieldMinimalContrast"));
		if (properties.getProperty(prefix+"flatFieldMinimalAccumulate")!=null)
			this.flatFieldMinimalAccumulate=Double.parseDouble(properties.getProperty(prefix+"flatFieldMinimalAccumulate"));
		if (properties.getProperty(prefix+"flatFieldShrinkForMatching")!=null)
			this.flatFieldShrinkForMatching=Double.parseDouble(properties.getProperty(prefix+"flatFieldShrinkForMatching"));
		if (properties.getProperty(prefix+"flatFieldMaxRelDiff")!=null)
			this.flatFieldMaxRelDiff=Double.parseDouble(properties.getProperty(prefix+"flatFieldMaxRelDiff"));
		if (properties.getProperty(prefix+"flatFieldShrinkMask")!=null)
			this.flatFieldShrinkMask=Integer.parseInt(properties.getProperty(prefix+"flatFieldShrinkMask"));
		if (properties.getProperty(prefix+"flatFieldFadeBorder")!=null)
			this.flatFieldFadeBorder=Double.parseDouble(properties.getProperty(prefix+"flatFieldFadeBorder"));
		if (properties.getProperty(prefix+"flatFieldResetMask")!=null)
			this.flatFieldResetMask=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldResetMask"));
		if (properties.getProperty(prefix+"flatFieldShowSensorMasks")!=null)
			this.flatFieldShowSensorMasks=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldShowSensorMasks"));
		if (properties.getProperty(prefix+"flatFieldShowIndividual")!=null)
			this.flatFieldShowIndividual=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldShowIndividual"));
		if (properties.getProperty(prefix+"flatFieldShowResult")!=null)
			this.flatFieldShowResult=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldShowResult"));
		if (properties.getProperty(prefix+"flatFieldApplyResult")!=null)
			this.flatFieldApplyResult=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldApplyResult"));
		if (properties.getProperty(prefix+"flatFieldUseInterpolate")!=null)
			this.flatFieldUseInterpolate=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldUseInterpolate"));
		if (properties.getProperty(prefix+"flatFieldMaskThresholdOcclusion")!=null)
			this.flatFieldMaskThresholdOcclusion=Double.parseDouble(properties.getProperty(prefix+"flatFieldMaskThresholdOcclusion"));
		if (properties.getProperty(prefix+"flatFieldShrinkOcclusion")!=null)
			this.flatFieldShrinkOcclusion=Integer.parseInt(properties.getProperty(prefix+"flatFieldShrinkOcclusion"));
		if (properties.getProperty(prefix+"flatFieldFadeOcclusion")!=null)
			this.flatFieldFadeOcclusion=Double.parseDouble(properties.getProperty(prefix+"flatFieldFadeOcclusion"));
		if (properties.getProperty(prefix+"flatFieldIgnoreSensorFlatField")!=null)
			this.flatFieldIgnoreSensorFlatField=Boolean.parseBoolean(properties.getProperty(prefix+"flatFieldIgnoreSensorFlatField"));
		if (properties.getProperty(prefix+"repeatFlatFieldSensor")!=null)
			this.repeatFlatFieldSensor=Integer.parseInt(properties.getProperty(prefix+"repeatFlatFieldSensor"));


		if (properties.getProperty(prefix+"specularHighPassSigma")!=null)
			this.specularHighPassSigma=Double.parseDouble(properties.getProperty(prefix+"specularHighPassSigma"));
		if (properties.getProperty(prefix+"specularLowPassSigma")!=null)
			this.specularLowPassSigma=Double.parseDouble(properties.getProperty(prefix+"specularLowPassSigma"));
		if (properties.getProperty(prefix+"specularDiffFromAverageThreshold")!=null)
			this.specularDiffFromAverageThreshold=Double.parseDouble(properties.getProperty(prefix+"specularDiffFromAverageThreshold"));
		if (properties.getProperty(prefix+"specularNumIter")!=null)
			this.specularNumIter=Integer.parseInt(properties.getProperty(prefix+"specularNumIter"));
		if (properties.getProperty(prefix+"specularApplyNewWeights")!=null)
			this.specularApplyNewWeights=Boolean.parseBoolean(properties.getProperty(prefix+"specularApplyNewWeights"));
		if (properties.getProperty(prefix+"specularPositiveDiffOnly")!=null)
			this.specularPositiveDiffOnly=Boolean.parseBoolean(properties.getProperty(prefix+"specularPositiveDiffOnly"));
		if (properties.getProperty(prefix+"specularShowDebug")!=null)
			this.specularShowDebug=Integer.parseInt(properties.getProperty(prefix+"specularShowDebug"));

		if (properties.getProperty(prefix+"is_small")!=null) this.is_small=Boolean.parseBoolean(properties.getProperty(prefix+"is_small"));

		if (properties.getProperty(prefix+"SMALL."+"is_small")!=null) {
			this.refineParametersSmall = new RefineParameters();
			this.refineParametersSmall.getProperties(prefix+"SMALL.", properties);
		}

	}

	public int showDialog(String title, int parMask, int numSeries, double [] averageRGB, boolean include_lwir) { // average RGB?
		int rslt = showDialogThis(title, parMask, numSeries, averageRGB, null);
		if ((rslt < 0) || !include_lwir) {// (this.refineParametersSmall == null) || ) {
			return rslt;
		}
		if (this.refineParametersSmall == null) {
			this.refineParametersSmall = new RefineParameters(); // maybe add small defaults?
		}
		return this.refineParametersSmall.showDialogThis(title+"-small (LWIR) sensors", parMask, numSeries, averageRGB, this);

	}
	public int showDialogThis(String title, int parMask, int numSeries, double [] averageRGB, RefineParameters large) {
		// sensor 0xfff, grid - 0xcc0 // cannot show result (cumulative) grid correction
		GenericDialog gd = new GenericDialog(title);
		if (numSeries>=0) gd.addNumericField("Fitting strategy series number (selects images to process) ", numSeries,0);
		if ((parMask&0x200000)!=0) gd.addNumericField("Repeat target/sensor flat-field calculation", this.repeatFlatFieldSensor,0,3,"times");
		// new parameters for expanding/blurring sensor distortions and photometrics
		if ((parMask&0x1000000)!=0) gd.addMessage("Parameters for extrapolating sensors distortions and vignetting to the frame edges (2020)");
		if ((parMask&0x1000000)!=0) gd.addNumericField("Center fraction (of half-height) used for orthogonal Gaussian blur", this.center_fract, 2,6,"");
		if ((parMask&0x1000000)!=0) gd.addNumericField("Transition fraction (of half-height) between center ortho and peripheral polar extrapolatioon", this.transit_fract, 2,6,"");
		if ((parMask&0x1000000)!=0) gd.addNumericField("Angular Gaussian sigma", this.gaus_ang, 2,6,"radian");
		if ((parMask&0x1000000)!=0) gd.addNumericField("Radial and ortho (central) Gaussian sigma (fraction of half-height)", this.gaus_rad, 2,6,"");
		if ((parMask&0x1000000)!=0) gd.addNumericField("Maximal high-frequency error to remove from the second pass (for distortions)", this.max_diff_err_geom, 2,6,"");
		if ((parMask&0x1000000)!=0) gd.addNumericField("Maximal high-frequency error to remove from the second pass (for vignetting)", this.max_diff_err_photo, 2,6,"");
		if ((parMask&0x1000000)!=0) gd.addMessage("---");


		//sensorExtrapolateDiff
///		if ((parMask&0x80000)!=0) gd.addCheckbox("Extrapolate incremetal (not checked - cumulative) correction",  this.sensorExtrapolateDiff);
///		if ((parMask&0x80000) !=0) gd.addNumericField("Shrink-blur combined sigma", this.sensorShrinkBlurComboSigma, 2,6,"sensor pixels"); // 20
///		if ((parMask&0x80000) !=0) gd.addNumericField("Shrink-blur combined level (-1..+1)", this.sensorShrinkBlurComboLevel, 2,6,""); // 0
///		if ((parMask&0x80000) !=0) gd.addNumericField("Combined alpha extrapolation threshold", this.sensorAlphaThreshold, 2,6,""); // normalize later?
///		if ((parMask&0x80000) !=0) gd.addNumericField("Extrapolation seed step",this.sensorStep, 1,4,"decimated pixels");
///		if ((parMask&0x80000) !=0) gd.addNumericField("Extrapolation gaussian sigma", this.sensorInterpolationSigma, 2,6,"sensor pixels"); // 50
///		if ((parMask&0x80000) !=0) gd.addNumericField("Extrapolation effective radius (doubling sigma in tangential direction)", this.sensorTangentialRadius, 2,6,"fraction of full image radius");
///		if ((parMask&0x80000) !=0) gd.addNumericField("Extrapolation half-square side for polynomial approximation", this.sensorScanDistance, 0,3,"sensor pixels");
///		if ((parMask&0x80000) !=0) gd.addNumericField("Extrapolation half-square side for extrapolation", this.sensorResultDistance, 0,3,"sensor pixels");
///		if ((parMask&0x80000) !=0) gd.addNumericField("Extrapolation polynomial degree", this.sensorInterpolationDegree, 0,1,"");

		if ((parMask&0x100000)!=0) gd.addNumericField("Fitting series number (to select images), negative - use all enabled images (eo -> lwir)", this.flatFieldSerNumber,0);
		if ((parMask&0x100000)!=0) gd.addNumericField("Reference station number (unity target brightness)  (eo -> lwir)", this.flatFieldReferenceStation,0);
		if ((parMask&0x100000)!=0) gd.addNumericField("Shrink sensor mask",    this.flatFieldShrink, 1,6,"sensor pix");
		if ((parMask&0x100000)!=0) gd.addNumericField("Non-vignetted radius", this.flatFieldNonVignettedRadius, 1,6,"sensor pix");
		if ((parMask&0x100000)!=0) gd.addNumericField("Minimal alpha",        100.0*this.flatFieldMinimalAlpha, 3,7,"%");

		if ((parMask&0x100000)!=0) gd.addNumericField("Minimal contrast (occlusion detection)", this.flatFieldMinimalContrast, 3,7,"(0 .. ~0.8");
///		if ((parMask&0x100000)!=0) gd.addNumericField("Minimal alpha for accumulation", 100.0*this.flatFieldMinimalAccumulate, 3,7,"%");

		if ((parMask&0x100000)!=0) gd.addNumericField("Shrink pattern for matching (eo -> lwir)", this.flatFieldShrinkForMatching, 3,7,"grid nodes");
		if ((parMask&0x100000)!=0) gd.addNumericField("Maximal relative difference between nodes (eo -> lwir)", 100.0*this.flatFieldMaxRelDiff, 3,7,"%");
		if ((parMask&0x100000)!=0) gd.addNumericField("Shrink pattern border (eo -> lwir)", this.flatFieldShrinkMask, 0,3,"grid nodes");
		if ((parMask&0x100000)!=0) gd.addNumericField("Fade pattern border (eo -> lwir)", this.flatFieldFadeBorder, 3,7,"grid nodes");
		if ((parMask&0x100000)!=0) gd.addMessage("Update pattern white balance (if the illumination is yellowish, increase red and green here)");

		if ((parMask&0x100000)!=0) gd.addNumericField("Average grid RED   (1.0 for white)",  averageRGB[0], 3,5,"x"); //
		if ((parMask&0x100000)!=0) gd.addNumericField("Average grid GREEN (1.0 for white)",  averageRGB[1], 3,5,"x"); //
		if ((parMask&0x100000)!=0) gd.addNumericField("Average grid BLUE  (1.0 for white)",  averageRGB[2], 3,5,"x"); //
		if ((parMask&0x100000)!=0) gd.addCheckbox("Reset pattern mask  (eo -> lwir)",               this.flatFieldResetMask);
		if ((parMask&0x100000)!=0) gd.addCheckbox("Show non-vignetting sensor masks", this.flatFieldShowSensorMasks);
		if ((parMask&0x100000)!=0) gd.addCheckbox("Show per-sensor patterns",         this.flatFieldShowIndividual);
		if ((parMask&0x100000)!=0) gd.addCheckbox("Show result mask",                 this.flatFieldShowResult);
		if ((parMask&0x100000)!=0) gd.addCheckbox("Apply pattern flat field and mask",this.flatFieldApplyResult);
		if ((parMask&0x100000)!=0) gd.addCheckbox("Use interpolation for sensor correction",this.flatFieldUseInterpolate);
		if ((parMask&0x100000)!=0) gd.addNumericField("Suspect occlusion only if grid is missing in the area where sensor mask is above this threshold",100.0* this.flatFieldMaskThresholdOcclusion, 3,7,"%");
		if ((parMask&0x100000)!=0) gd.addNumericField("Expand suspected occlusion  area", this.flatFieldShrinkOcclusion, 0,3,"grid nodes");
		if ((parMask&0x100000)!=0) gd.addNumericField("Fade grid on image (occlusion handling)", this.flatFieldFadeOcclusion, 3,7,"grid nodes");
		if ((parMask&0x100000)!=0) gd.addCheckbox("Ignore existent sensor flat-field calibration",this.flatFieldIgnoreSensorFlatField);

		if ((parMask&0x400000)!=0) gd.addMessage("Specular reflections removal parameters:");
		if ((parMask&0x400000)!=0) gd.addCheckbox("Apply new (after removal of specular reflections) weights",           this.specularApplyNewWeights);
		if ((parMask&0x400000)!=0) gd.addCheckbox("Process only positive difference from average",                       this.specularPositiveDiffOnly);
		if ((parMask&0x400000)!=0) gd.addNumericField("High-pass sigma for difference from average (to detect specular)",this.specularHighPassSigma, 3,7,"pix");
		if ((parMask&0x400000)!=0) gd.addNumericField("Low-pass sigma for difference from average (to detect specular)",this.specularLowPassSigma, 3,7,"pix");

		if ((parMask&0x400000)!=0) gd.addNumericField("Difference from average threshold",                          100.0*this.specularDiffFromAverageThreshold, 3,7,"%");
		if ((parMask&0x400000)!=0) gd.addNumericField("Number of iterations for calculating average",                    this.specularNumIter, 0);

		if ((parMask&0x400000)!=0) gd.addNumericField("Debug show mode (0 - off, 1 - last iteration only, 2 - all iterations)",this.specularShowDebug, 0);


///		if ((parMask &     1) !=0) gd.addCheckbox    ("Extrapolate correction results", this.extrapolate);
///		if ((parMask &     2) !=0) gd.addNumericField("Threshold alpha (discard pixels with mask below that value)", this.alphaThreshold,3);
///		if ((parMask &0x8000) !=0) gd.addNumericField("Fat zero for color trasfer functions", this.fatZero,3);
///		if ((parMask &     4) !=0) gd.addNumericField("Fitting radius for extrapolation, Gaussian weight function sigma (in non-decimated pixels) ",    this.extrapolationSigma,3);
///		if ((parMask &     8) !=0) gd.addNumericField("Fitting scan half-size of the square, in multiples of Fitting Radius", this.extrapolationKSigma,3);
///		if ((parMask &  0x10) !=0) gd.addCheckbox    ("Apply smoothing to the correction results", this.smoothCorrection);
///		if ((parMask &  0x20) !=0) gd.addNumericField("Smoothing sigma, in non-decimated pixels",  this.smoothSigma,3);
		if ((parMask &  0x40) !=0) gd.addCheckbox    ("Apply correction",                          this.applyCorrection);
		if ((parMask&0x40000) !=0) gd.addCheckbox    ("Apply correction",                          this.targetApplyCorrection);
		if ((parMask &0x4000) !=0) gd.addCheckbox    ("Apply flat-field correction",               this.applyFlatField);
		if ((parMask &  0x80) !=0) gd.addNumericField("Scale correction before applying",          this.correctionScale,3);
		if ((parMask&0x40000) !=0) gd.addNumericField("Scale correction before applying",          this.targetCorrectionScale,3);
		if ((parMask & 0x100) !=0) gd.addCheckbox    ("Show result (cumulative) correction",       this.showCumulativeCorrection);
		if ((parMask & 0x200) !=0) gd.addCheckbox    ("Show additional correction before blurring",this.showUnfilteredCorrection);
		if ((parMask & 0x200) !=0) gd.addCheckbox    ("Show correction extrapolatiuon (polar)",    this.showExtrapolationCorrection);
		if ((parMask & 0x400) !=0) gd.addCheckbox    ("Show this step (additional) correction",      this.showThisCorrection);
		if ((parMask&0x40000) !=0) gd.addCheckbox    ("Show this (additional) target correction",    this.targetShowThisCorrection);
		if ((parMask & 0x800) !=0) gd.addCheckbox    ("Show individual, per-image residuals",        this.showPerImage); // used in 2020
		if ((parMask&0x40000) !=0) gd.addCheckbox    ("Show individual, per-image target residuals", this.targetShowPerImage);
		if ((parMask&0x10000) !=0) gd.addNumericField("Show individual residuals for image number (<0 - all images)", this.showIndividualNumber,0);
		if ((parMask &0x1000) !=0) gd.addCheckbox    ("Correct patetrn grid node locations in 3d (false - in 2d only)",  this.grid3DCorrection);
		if ((parMask &0x1000) !=0) gd.addCheckbox    ("Rotate final 3d pattern correction (?)",   this.rotateCorrection);
		if ((parMask&0x20000) !=0) gd.addNumericField("Maximal Z-axis correction (if more will fall back to 2d correction algorithm)", this.grid3DMaximalZCorr,1,3,"mm");
		if ((parMask&0x20000) !=0) gd.addCheckbox    ("Use Z-variations of the pattern for different stations",   this.useVariations);
		if ((parMask&0x20000) !=0) gd.addNumericField("Penalty for different Z for the same target nodes for different stations", 100.0*this.variationPenalty,3,7,"%");
		if ((parMask&0x20000) !=0) gd.addCheckbox    ("Keep X and Y pattern correction, adjust only Z",this.fixXY);
		if ((parMask&0x20000) !=0) gd.addCheckbox    ("Reset previous Z variations before calculating the new one",   this.resetVariations);
		if ((parMask&0x20000) !=0) gd.addCheckbox    ("Do not fall back to 2-d calculation if 3d fails",   this.noFallBack);
		if (large != null)         gd.addCheckbox    ("Copy all parameters from the main sensors (large), ignore other filedsr",   false);
		//large
///		if ((parMask &0x2000) !=0) gd.addCheckbox    ("Use pattern grid alpha data",  this.usePatternAlpha);
		WindowTools.addScrollBars(gd);
		gd.showDialog();
		if (gd.wasCanceled()) return -1;
		int selectedSeries=0;
		if (numSeries>=0)          selectedSeries=          (int) gd.getNextNumber();
		if ((parMask&0x200000)!=0) this.repeatFlatFieldSensor=  (int) gd.getNextNumber();
		// new in 2020
		if ((parMask&0x1000000)!=0) this.center_fract=                    gd.getNextNumber();
		if ((parMask&0x1000000)!=0) this.transit_fract=                   gd.getNextNumber();
		if ((parMask&0x1000000)!=0) this.gaus_ang=                        gd.getNextNumber();
		if ((parMask&0x1000000)!=0) this.gaus_rad=                        gd.getNextNumber();
		if ((parMask&0x1000000)!=0) this.max_diff_err_geom=               gd.getNextNumber();
		if ((parMask&0x1000000)!=0) this.max_diff_err_photo=              gd.getNextNumber();

///		if ((parMask&0x80000) !=0) this.sensorExtrapolateDiff=          gd.getNextBoolean();
///		if ((parMask&0x80000) !=0) this.sensorShrinkBlurComboSigma=     gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorShrinkBlurComboLevel=     gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorAlphaThreshold=           gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorStep=                     gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorInterpolationSigma=       gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorTangentialRadius=         gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorScanDistance=       (int) gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorResultDistance=     (int) gd.getNextNumber();
///		if ((parMask&0x80000) !=0) this.sensorInterpolationDegree=(int) gd.getNextNumber();

		if ((parMask&0x100000)!=0) this.flatFieldSerNumber=         (int) gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldReferenceStation=  (int) gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldShrink=                  gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldNonVignettedRadius=      gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldMinimalAlpha=       0.01*gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldMinimalContrast=         gd.getNextNumber();
///		if ((parMask&0x100000)!=0) this.flatFieldMinimalAccumulate=  0.01*gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldShrinkForMatching=       gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldMaxRelDiff=         0.01*gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldShrinkMask=        (int) gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldFadeBorder=              gd.getNextNumber();
		if ((parMask&0x100000)!=0) averageRGB[0]=                         gd.getNextNumber();
		if ((parMask&0x100000)!=0) averageRGB[1]=                         gd.getNextNumber();
		if ((parMask&0x100000)!=0) averageRGB[2]=                         gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldResetMask=               gd.getNextBoolean();
		if ((parMask&0x100000)!=0) this.flatFieldShowSensorMasks=         gd.getNextBoolean();
		if ((parMask&0x100000)!=0) this.flatFieldShowIndividual=          gd.getNextBoolean();
		if ((parMask&0x100000)!=0) this.flatFieldShowResult=              gd.getNextBoolean();
		if ((parMask&0x100000)!=0) this.flatFieldApplyResult=             gd.getNextBoolean();
		if ((parMask&0x100000)!=0) this.flatFieldUseInterpolate=          gd.getNextBoolean();
		if ((parMask&0x100000)!=0) this.flatFieldMaskThresholdOcclusion=0.01*gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldShrinkOcclusion=   (int) gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldFadeOcclusion=           gd.getNextNumber();
		if ((parMask&0x100000)!=0) this.flatFieldIgnoreSensorFlatField=   gd.getNextBoolean();
		if ((parMask&0x400000)!=0) this.specularApplyNewWeights=              gd.getNextBoolean();
		if ((parMask&0x400000)!=0) this.specularPositiveDiffOnly=             gd.getNextBoolean();
		if ((parMask&0x400000)!=0) this.specularHighPassSigma=                gd.getNextNumber();
		if ((parMask&0x400000)!=0) this.specularLowPassSigma=                gd.getNextNumber();
		if ((parMask&0x400000)!=0) this.specularDiffFromAverageThreshold=0.01*gd.getNextNumber();;
		if ((parMask&0x400000)!=0) this.specularNumIter=                (int) gd.getNextNumber();
		if ((parMask&0x400000)!=0) this.specularShowDebug=              (int) gd.getNextNumber();


//		if ((parMask &     1) !=0) this.extrapolate=              gd.getNextBoolean();
//		if ((parMask &     2) !=0) this.alphaThreshold=           gd.getNextNumber();
//		if ((parMask &0x8000) !=0) this.fatZero=                  gd.getNextNumber();
//		if ((parMask &     4) !=0) this.extrapolationSigma=       gd.getNextNumber();
//		if ((parMask &     8) !=0) this.extrapolationKSigma=      gd.getNextNumber();
//		if ((parMask &  0x10) !=0) this.smoothCorrection=         gd.getNextBoolean();
//		if ((parMask &  0x20) !=0) this.smoothSigma=              gd.getNextNumber();
		if ((parMask &  0x40) !=0) this.applyCorrection=          gd.getNextBoolean();
		if ((parMask&0x40000) !=0) this.targetApplyCorrection=    gd.getNextBoolean();
		if ((parMask &0x4000) !=0) this.applyFlatField=           gd.getNextBoolean();
		if ((parMask &  0x80) !=0) this.correctionScale=          gd.getNextNumber();
		if ((parMask&0x40000) !=0) this.targetCorrectionScale=    gd.getNextNumber();
		if ((parMask & 0x100) !=0) this.showCumulativeCorrection= gd.getNextBoolean();
		if ((parMask & 0x200) !=0) this.showUnfilteredCorrection= gd.getNextBoolean();
		if ((parMask & 0x200) !=0) this.showExtrapolationCorrection= gd.getNextBoolean();
		if ((parMask & 0x400) !=0) this.showThisCorrection=       gd.getNextBoolean();
		if ((parMask&0x40000) !=0) this.targetShowThisCorrection= gd.getNextBoolean();
		if ((parMask & 0x800) !=0) this.showPerImage=             gd.getNextBoolean();
		if ((parMask&0x40000) !=0) this.targetShowPerImage=       gd.getNextBoolean();
		if ((parMask&0x10000) !=0) this.showIndividualNumber=(int)gd.getNextNumber();
		if ((parMask &0x1000) !=0) this.grid3DCorrection=         gd.getNextBoolean();
		if ((parMask &0x1000) !=0) this.rotateCorrection=         gd.getNextBoolean();
		if ((parMask&0x20000) !=0) this.grid3DMaximalZCorr=       gd.getNextNumber();

		if ((parMask&0x20000) !=0) this.useVariations=            gd.getNextBoolean();
		if ((parMask&0x20000) !=0) this.variationPenalty =   0.01*gd.getNextNumber();
		if ((parMask&0x20000) !=0) this.fixXY=                    gd.getNextBoolean();
		if ((parMask&0x20000) !=0) this.resetVariations=          gd.getNextBoolean();
		if ((parMask&0x20000) !=0) this.noFallBack=               gd.getNextBoolean();
//		if ((parMask &0x2000) !=0) this.usePatternAlpha=          gd.getNextBoolean();
		if ((large != null) && gd.getNextBoolean()) {
			copyFrom(large);
			return this.refineParametersSmall.showDialogThis(title, parMask, numSeries, averageRGB, large);
		}
		return selectedSeries;
	}

	private void copyFrom(RefineParameters other) {
		this.center_fract =                    other.center_fract;
		this.transit_fract =                   other.transit_fract;
		this.gaus_ang =                        other.gaus_ang;
		this.gaus_rad =                        other.gaus_rad;
		this.max_diff_err_geom =               other.max_diff_err_geom;
		this.max_diff_err_photo =              other.max_diff_err_photo;

		this.extrapolate=                      other.extrapolate;
		this.alphaThreshold=                   other.alphaThreshold;
		this.fatZero=                          other.fatZero;        // when extrapolatging color transfer coefficients (flat field) use this for logariphm
		this.extrapolationSigma=               other.extrapolationSigma;
		this.extrapolationKSigma=              other.extrapolationKSigma;
		this.smoothCorrection=                 other.smoothCorrection;
		this.smoothSigma=                      other.smoothSigma;
		this.correctionScale=                  other.correctionScale;
		this.showCumulativeCorrection=         other.showCumulativeCorrection;
		this.showUnfilteredCorrection=         other.showUnfilteredCorrection;
		this.showExtrapolationCorrection=      other.showExtrapolationCorrection;
		this.showThisCorrection=               other.showThisCorrection;
		this.showPerImage=                     other.showPerImage;
		this.showIndividualNumber=             other.showIndividualNumber; // which image to show (-1 - all)
		this.applyCorrection=                  other.applyCorrection;
		this.applyFlatField=                   other.applyFlatField;
		this.grid3DCorrection=                 other.grid3DCorrection;
		this.rotateCorrection=                 other.rotateCorrection; // not clear
		this.grid3DMaximalZCorr=               other.grid3DMaximalZCorr; // Maximal Z-axis correc tion (if more will fall back to 2d correction algorithm)
		this.useVariations=                    other.useVariations;
		this.variationPenalty=                 other.variationPenalty; // "stiffness" of individual (per-station) Z-values of the target pattern
		this.fixXY=                            other.fixXY;
		this.resetVariations=                  other.resetVariations;
		this.noFallBack=                       other.noFallBack; // may have bugs - not tested yet
		this.usePatternAlpha=                  other.usePatternAlpha;
		this.targetShowPerImage=               other.targetShowPerImage;
		this.targetShowThisCorrection=         other.targetShowThisCorrection;
		this.targetApplyCorrection=            other.targetApplyCorrection;
		this.targetCorrectionScale=            other.targetCorrectionScale;
		this.sensorExtrapolateDiff=            other.sensorExtrapolateDiff;
		this.sensorShrinkBlurComboSigma=       other.sensorShrinkBlurComboSigma;
		this.sensorShrinkBlurComboLevel=       other.sensorShrinkBlurComboLevel;
		this.sensorAlphaThreshold=             other.sensorAlphaThreshold;
		this.sensorStep=                       other.sensorStep;
		this.sensorInterpolationSigma=         other.sensorInterpolationSigma;
		this.sensorTangentialRadius=           other.sensorTangentialRadius;
		this.sensorScanDistance=               other.sensorScanDistance;
		this.sensorResultDistance=             other.sensorResultDistance;
		this.sensorInterpolationDegree=        other.sensorInterpolationDegree;
		this.flatFieldSerNumber=               other.flatFieldSerNumber;
		this.flatFieldReferenceStation=        other.flatFieldReferenceStation;
		this.flatFieldShrink=                  other.flatFieldShrink;
		this.flatFieldNonVignettedRadius=      other.flatFieldNonVignettedRadius;
		this.flatFieldMinimalAlpha=            other.flatFieldMinimalAlpha;
		this.flatFieldMinimalContrast=         other.flatFieldMinimalContrast;
		this.flatFieldMinimalAccumulate=       other.flatFieldMinimalAccumulate;
		this.flatFieldShrinkForMatching=       other.flatFieldShrinkForMatching;
		this.flatFieldMaxRelDiff=              other.flatFieldMaxRelDiff;
		this.flatFieldShrinkMask=              other.flatFieldShrinkMask;
		this.flatFieldFadeBorder=              other.flatFieldFadeBorder;
		this.flatFieldResetMask=               other.flatFieldResetMask;
		this.flatFieldShowSensorMasks=         other.flatFieldShowSensorMasks;
		this.flatFieldShowIndividual=          other.flatFieldShowIndividual;
		this.flatFieldShowResult=              other.flatFieldShowResult;
		this.flatFieldApplyResult=             other.flatFieldApplyResult;
		this.flatFieldUseInterpolate=          other.flatFieldUseInterpolate;
		this.flatFieldMaskThresholdOcclusion=  other.flatFieldMaskThresholdOcclusion;
		this.flatFieldShrinkOcclusion=         other.flatFieldShrinkOcclusion;
		this.flatFieldFadeOcclusion=           other.flatFieldFadeOcclusion;
		this.flatFieldIgnoreSensorFlatField=   other.flatFieldIgnoreSensorFlatField;
		this.repeatFlatFieldSensor=            other.repeatFlatFieldSensor;

		this.specularHighPassSigma=            other.specularHighPassSigma;
		this.specularLowPassSigma=             other.specularLowPassSigma;
		this.specularDiffFromAverageThreshold= other.specularDiffFromAverageThreshold;
		this.specularNumIter=                  other.specularNumIter;
		this.specularApplyNewWeights=          other.specularApplyNewWeights;
		this.specularPositiveDiffOnly=         other.specularPositiveDiffOnly;
		this.specularShowDebug=                other.specularShowDebug;
	}


	public boolean isSmall() {
		return is_small;
	}

	public RefineParameters getSmall() {
		return this.refineParametersSmall;
	}
}