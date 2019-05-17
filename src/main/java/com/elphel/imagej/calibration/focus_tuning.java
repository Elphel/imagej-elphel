package com.elphel.imagej.calibration;
/**
** -----------------------------------------------------------------------------**
** focus_tuning.java
**
** Measures focus sharpnes at orthogonal directions, differenct colors,
** Displays results for manual focusing/image plane tilting.
** NOTE: Requires special targets !
**
** Copyright (C) 2010 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**  
**  focus_tuning.java is free software: you can redistribute it and/or modify
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

import ij.*;
import ij.process.*;
import ij.gui.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import ij.plugin.frame.*;
import ij.text.*;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;

import org.w3c.dom.Document;
import org.w3c.dom.Node;


//import java.io.File;
import java.net.*;

import org.xml.sax.SAXException;

import com.elphel.imagej.jp4.JP46_Reader_camera;

import javax.xml.parsers.ParserConfigurationException;

import java.lang.Integer;


public class focus_tuning extends PlugInFrame implements ActionListener {
  /**
	 * 
	 */
	private static final long serialVersionUID = -6722464130506115169L;
JP46_Reader_camera jp4_instance;
//  MTF_Bayer MTF_Bayer_instance;
  Panel panel;
  static Frame instance;
  String crossTalkName="Crosstalk";
  String simulName="Simulation";

 public static int DEBUG_LEVEL = 1;
 public static int MASTER_DEBUG_LEVEL = 1;

 public static int FFTSize=32;
 public static int FFTScanStep=8;
 
 public static int test_x=FFTSize;
 public static int test_y=FFTSize;
 public static double targetFTolerance=0.2; /// tolerance in frequency domain
 public static double targetF = 0.23;         /// target frequency (cpp, )
 public static double targetAngle=5.0;  /// target angle
 public static int     correlationRadius=200; // in original pixels
 public static int     equalizationRadius=500; // in original pixels
 public static int     minTargetSize=   100; // in original pixels, used in minFilter()
 public static int     zonesVert=3;
 public static int     zonesHor=3;
 public static double  targetThreshold=0.04; // 0.01;
// public static int     calib_dx=0; // shift from the center of the pattern to white/black calibration areas
// public static int     calib_dy=90; // shift from the center of the pattern to white/black calibration areas
// public static double  calib_bw2hv=0.85; // relative distance to Black and white from H to V
 public static double  calib_bw2hv=0.77; // 210/270 - relative distance to Black and white from H to V
 public static int     calib_sample=16; // calibration sample size (square, original pixels
 public static int     bayerComponent=-1; //all
 public static int     displayWidth=800;
 public static int     displayHeight=600;
 public static int     displayRotate=0; /* 90- degrees CCV from the image */
 public static double  goalG2R=1.05; /* Level of G/R when R/G indicator changes colors */

 public static int     displayZonesGap=10;
 public static int     displayElementGap=2;
 public static double  displayMin=1.0;
 public static double  displayMax=300.0;
 public static int     displayVDivs=10;
 private static float [] Hamming;
 private static float [] TargetFilter;
 public static float [][] input_bayer; 
 public static double [][][] targets; /// [v][h][y,x,value] (value=0 - invalid)
 public static int [][][] targetLocations=null; /* [number][v,h,white,black][x,y] */
 public static int [][] targetZones=null; /* [number][zone_x,zone_y]*/
 public ImageProcessor ip_display;
 public ImagePlus      imp_display;
 public ImagePlus      imp_camera=null;
 public static String  motorsIP = "192.168.0.236";
 public static double  motorsWaitAfterMove=1.0; // now seconds
 public static int     motors_hysteresis=100;
 public static int     motors_scan_range=5000; // double amplitude. will start current position-motors_scan_range-motors_hysteresis, scan to current position+motors_scan_range and return to current position
 public static int     motors_scan_steps= 100;
 public static int     motors_scan_steps_fast=8;
 public static int     motors_scan_averageNum=5;
 public static double[]scanScales={1.0,1.0,1.0}; // move motors with theses scales during scan
// theoretical values: {-0.808,-0.808,+1.0}
// theoretical values: {+1.0,  -1.0,  -0.212}
 public static int []  motorPositions;
 public static int []  motorTargetPositions= new int[3];
 public static int []  startScan= new int[3];
 public static int []  endScan= new int[3];
 public static double [][][][][] scanResults;
// public static double [][]    motorsRange={{-7000,7000},{-7000,7000},{-7000,7000}};
 public static double [][]    motorsRange={{-15000,15000},{-15000,15000},{-15000,15000}}; // temporarily
 public static boolean        motorScanBothDirections=false;
 public static double [] motors2center={0.25,0.25,0.5}; // how much moves the center target when you move one of the motors

 public static int      motors_calibSteps= 50;     // motor steps / calibration step (all 3 motors the same)
 public static double   stepsPerUM=        89.0;   // number of motor steps per 1 um (all 3 motors the same) theoretical - 68, measured 89. R-G ~=8.15um
 public static int      goalNumber=        12;     // See coments in the function makeOffsetsPerGoal() explaining goals
 public static double   goalGreenWeight=   0.7;    // Goal (for numbers 4,8,12) will be 0.7 from goalGreenWeight:(1-goalGreenWeight) between green and red

 public static double  [][] offsetsPerGoal=null;   // [goalNumer][targetNumber] offsets from selected goal (1==one calibration table line,
                                                   // 1 calibration line =   motors_calibSteps steps, == (motors_calibSteps)/stepsPerUM microns
 public static double   focusOffset =      7.0; // um goal offseet from the best focus at target distance (9um will be for infinity). Will be subtracted from results
                                                // Target was @ 2260mm , infinity - 8.8um, 15m - 7.6
                                                // 10m - 6.95um
                                                // 5m  - 4.92 um.
                                                // Focusing @10 meters would give  5m..infinity (if we have overlap at 5m horizontal?)
 public static int       motors_calibGoalWindow=      5;   // half-window width when looking for maximum in a goal function (first - just maximum, then - mpy by huffaman and find centroid)
 public static double    motors_calibMaxSlopes=    0.85;   // When looking for maximums, use interval araound maximum where points are above (min+motors_calibMaxSlopes*(max-min))
 public static int       motor3rd_signed_hysteresis=100;
 public static int       motor3rd_motorNumber=        3;
 public static int       motor3rd_targetNumber=       3;
 public static int       motor3rd_minStep=          100;
 public static int       motor3rd_maxStep=          500;
 public static int       motor3rd_direction=          1;
 public static int       motor3rd_colorMask=          2; // red
 public static double    motor3rd_V2H=                1.0; // obsolete
 public static double    motor3rd_T2R=                1.5; //tangential resolution is more important than radial (low distortion lens, angular resolution)
 public static int       motor3rd_refine=             6;
 public static int       motor3rd_refineAverage=      4; // number of measurements during refine stage
 public static int       motor3rd_refineExpand=       4; // expand refine range if fails


 public static int       minUpdateStep =            100; // update thisStepsPerUM only if the last step was > this threshold
 public static double    variableStepRatio=         2.0; // 2.0 - maximal thisStepsPerUM change in one step



 public static int       motor3rd_maxIterations=     20;
 public static int       motor3rd_fastIterations=     2;


 public static int       motor3rd_M1range=         4560;
 public static int       motor3rd_M1steps=           10;
 public static int       motor3rd_M2range=         4560;
 public static int       motor3rd_M2steps=           10;
// public static int       motor3rd_initialStep=     1824; //1/2 turn
 public static int       motor3rd_initialStep=      912; //1/4 turn
 public static int       motor3rd_numIterations=     10; // ADJUST ALL ITERATIONS

 public static boolean   motor3rd_useCalibration=    true; /// When moving motors - use calibration data (if available)
 public static double    unreliableFromEnds=        0.0; // when calculating astigmatims use caution when maximums are outside of the scanned range
                                                          // here positive values (in samples) decrease "reliable" range, negative - decrease

 public static  Color [] resolution_colors=  { new Color(1.0f,0.6f,0.6f),
                                               new Color(0.5f,0.0f,0.0f),
                                               new Color(0.6f,1.0f,0.6f),
                                               new Color(0.0f,0.5f,0.0f),
                                               new Color(0.6f,0.6f,1.0f),
                                               new Color(0.0f,0.0f,0.5f)};

 public static int    [] motor3rd_targets_dirs=        {-1,1,1,1,-1};
 public static int    [] motor3rd_targets_ColorMasks=  {2,2,0,2,2};
 public static double [] motor3rd_targets_T2R=         {1.0,1.0,1.05,1.0,1.0};
 public static int       MOTOR_UNDEFINED=0x8000; /// undefined motor position (could not find, out of limits)
 public static int       motors_sideWays=5;
 public static boolean   motors_measureDisplacements=true;
 public boolean run=false;
 public double [][][] targetScanData=null;
 public double [][][][] maxScanTargetOrientColor=null;
 public double badMeasurementError=8.0; //measurement is considered bad if it's difference from the average of others is more than that badMeasurementError * others DMS
 public double conservativeCoefficient=0.7; // when estimating expected measurement RMS, use conservativeCoefficient of the old one and (1-conservativeCoefficient) of the new one
 public static double [] oldErrors=null;

 private ImagePlus imp_src;
 
// Note: For panoramic applications and low geometric distortion lens targets 1,5 needs vertical resolution higher, 2,4 - horizontal

  Plot plotResult;
 

  public focus_tuning() {
    super("Focus Tuning");
    if (IJ.versionLessThan("1.39t")) return;
    if (instance!=null) {
      instance.toFront();
      return;
    }
    instance = this;
    addKeyListener(IJ.getInstance());

    setLayout(new FlowLayout());
    panel = new Panel();
    addButton("Configure");
    addButton("Split Bayer");
//    addButton("Test");
    addButton("Find Targets");
//    addButton("Measure Targets");
    addButton("Acquire and Measure");
//    addButton("Re-acquire and Measure");
    addButton("Configure Motors");
    addButton("Read Motors");
    addButton("Move Motors");
    addButton("Re-center Scan Focus");
    addButton("Scan Focus");
    addButton("Show Plots");
    addButton("Test 3-rd motor");
    addButton("Scan M1&M2, getM3");
    addButton("Adjust all motors");

/*    addButton("Run");
    addButton("Stop"); */
    add(panel);
    pack();
    GUI.center(this);
    setVisible(true);
    Hamming=initHamming(FFTSize);
    TargetFilter=initTargetFilter(FFTSize, targetFTolerance, targetF, targetAngle );
    initDisplay();
    jp4_instance=       new JP46_Reader_camera();
  }
  void addButton(String label) {
    Button b = new Button(label);
    b.addActionListener(this);
    b.addKeyListener(IJ.getInstance());
    panel.add(b);
  }
  public void actionPerformed(ActionEvent e) {
    int i;
    String label = e.getActionCommand();
    if (label==null) return;
    if (label.equals("Configure")) {
      if (showDialog()) {
        Hamming=     initHamming(FFTSize);
        TargetFilter=initTargetFilter(FFTSize, targetFTolerance, targetF, targetAngle );
        initDisplay();
        if (MASTER_DEBUG_LEVEL>2) {
          float [] pixels_f = new float[TargetFilter.length];
          System.arraycopy(TargetFilter, 0, pixels_f, 0, pixels_f.length);
          swapQuads(pixels_f,FFTSize);
          ImageProcessor ip_f= new FloatProcessor(FFTSize,FFTSize);
//          ip_f.setPixels(TargetFilter);
          ip_f.setPixels(pixels_f);
          ImagePlus imp_f=  new ImagePlus("Filter_weights", ip_f);
          imp_f.getProcessor().resetMinAndMax();
          imp_f.show();
        }
      }
      return;
    }
    if (label.equals("Configure Motors")) {
      if (showMotorsDialog()) {
/// Init +enable motors here?
      }
      return;
    }

/*
    if (label.equals("Acquire and Measure")) {
      imp_camera=acquireMeasureAndDisplay(null);
      return;
    } else if (label.equals("Re-acquire and Measure")) {
*/
    if (label.equals("Acquire and Measure")) {
      imp_camera=acquireMeasureAndDisplay(imp_camera);
      return;
    } else if (label.equals("Read Motors")) {
      motorPositions =readElphel10364Motors(motorsIP);

      IJ.showMessage("Motor positions",
                     "motor1="+motorPositions[0]+"\n"+
                     "motor2="+motorPositions[1]+"\n"+
                     "motor3="+motorPositions[2]);


      return;

    } else if (label.equals("Move Motors")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      if (moveMotorsDialog()) {
       // Do something?
      }

      return;

    } else if (label.equals("Re-center Scan Focus")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      motorPositions =readElphel10364Motors(motorsIP);
      startScan[0]=motorPositions[0]- ((int) (motors_scan_range*scanScales[0]/2));
      endScan[0]=  startScan[0]+      ((int) (motors_scan_range*scanScales[0]));
      startScan[1]=motorPositions[1]- ((int) (motors_scan_range*scanScales[1]/2));
      endScan[1]=  startScan[1]+      ((int) (motors_scan_range*scanScales[1]));
      startScan[2]=motorPositions[2]- ((int) (motors_scan_range*scanScales[2]/2));
      endScan[2]=  startScan[2]+      ((int) (motors_scan_range*scanScales[2]));

      motors_calibSteps=motors_scan_range/motors_scan_steps_fast; /// *************** Assuming all 3 motors with scale 1.0 !!!

      scanResults = scanMotors(imp_camera,motorsIP,motors_hysteresis, motorsWaitAfterMove,  startScan, endScan,  motors_scan_steps_fast, motor3rd_refineAverage, motorScanBothDirections);
      targetScanData= null;// invalidate;
      maxScanTargetOrientColor=calcMaxScanTargetOrientColor(scanResults, motors_calibGoalWindow, motors_calibMaxSlopes);
      if (DEBUG_LEVEL>1) {
        if (scanResults != null) showScanTable(scanResults,  startScan, endScan, motors_scan_steps_fast, maxScanTargetOrientColor, "Fast_Scan_Focus_Results" , 4); /// 4 - precision
      }
      int focusTarget=-1;
      double averMax;
      if (maxScanTargetOrientColor[0].length==5) focusTarget=2;
      else if (maxScanTargetOrientColor[0].length==0) focusTarget=0;
      if (focusTarget>=0) {
        averMax=(maxScanTargetOrientColor[0][focusTarget][0][0]+
                 maxScanTargetOrientColor[0][focusTarget][0][1]+
                 maxScanTargetOrientColor[0][focusTarget][1][0]+
                 maxScanTargetOrientColor[0][focusTarget][1][1])/4.0;
        if (averMax<0) averMax=0;
        else if (averMax > motors_scan_steps_fast) averMax = motors_scan_steps_fast;
        for (i=0;i<3;i++) {
          motorPositions[i]+= scanScales[i]* motors_scan_range* (2*averMax - motors_scan_steps_fast)/(2*motors_scan_steps_fast);
        }
        moveElphel10364Motors(motorsIP,motorPositions,1.0,true,"RE-CENTERING", true);
        motorPositions =readElphel10364Motors(motorsIP); /// re-read positions to make sure motors did settle down
        if (DEBUG_LEVEL>0)  IJ.showMessage("Recentered Focus Scan range",
                           "motor1="+motorPositions[0]+"\n"+
                           "motor2="+motorPositions[1]+"\n"+
                           "motor3="+motorPositions[2]);
      } else {
        IJ.showMessage("Error", "re-centering focus scan range works only with a single target or the center one of 5, found "+maxScanTargetOrientColor[0].length+" targets");
      }


      return;
    } else if (label.equals("Scan Focus")) {
      DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
      motorPositions =readElphel10364Motors(motorsIP);
      startScan[0]=motorPositions[0]- ((int) (motors_scan_range*scanScales[0]/2));
      endScan[0]=  startScan[0]+      ((int) (motors_scan_range*scanScales[0]));
      startScan[1]=motorPositions[1]- ((int) (motors_scan_range*scanScales[1]/2));
      endScan[1]=  startScan[1]+      ((int) (motors_scan_range*scanScales[1]));
      startScan[2]=motorPositions[2]- ((int) (motors_scan_range*scanScales[2]/2));
      endScan[2]=  startScan[2]+      ((int) (motors_scan_range*scanScales[2]));

      motors_calibSteps=motors_scan_range/motors_scan_steps; /// *************** Assuming all 3 motors with scale 1.0 !!!

//      scanResults = scanMotors(imp_camera,motorsIP,motors_hysteresis, motorsWaitAfterMove,  startScan, endScan,  motors_scan_steps,targets);
      scanResults = scanMotors(imp_camera,motorsIP,motors_hysteresis, motorsWaitAfterMove,  startScan, endScan,  motors_scan_steps, motor3rd_refineAverage, motorScanBothDirections);

      targetScanData= filterScanTable(scanResults,  motors_scan_averageNum);  // use running average to smouth data
      maxScanTargetOrientColor=calcMaxScanTargetOrientColor(scanResults, motors_calibGoalWindow, motors_calibMaxSlopes);


//      offsetsPerGoal= makeOffsetsPerGoal(scanResults, motor3rd_targets_T2R, motors_calibGoalWindow, motors_calibMaxSlopes);
      offsetsPerGoal= makeOffsetsPerGoal(scanResults, motor3rd_targets_T2R, motors_calibGoalWindow, motors_calibMaxSlopes, goalGreenWeight);

//, int motors_scan_steps, int average, boolean bothDirections
      if (scanResults != null) {
        showScanTable(scanResults,  startScan, endScan, motors_scan_steps, maxScanTargetOrientColor, ((imp_src==null)?"":imp_src.getTitle())+"_scan_focus_results" , 4); /// 4 - precision
        showAstigmatismTable(maxScanTargetOrientColor, //[scan direction][target][v,h][r,g,b] 1.0 - one scan sample
                                   motors_calibSteps/stepsPerUM,  // microns per 1 scan sample
                                   motors_scan_steps,  // scan range (samples, starting from 0)
                                   unreliableFromEnds, // maximums outside of [reliabilityThreshold..range-reliabilityThreshold are
                                                                // considered unreliable
                                   ((imp_src==null)?"":imp_src.getTitle())+"_astigmatism_table",                // window title
                                   2);   // output precision in decimals

      }
      return;

    } else if (label.equals("Show Plots")) {
       if (scanResults == null) {
         IJ.showMessage("Warning", "No results to plot");
       } else {
         plotScanTable(scanResults, motors_calibSteps/stepsPerUM, targetF,  ((imp_src==null)?"":imp_src.getTitle())+" resolution", resolution_colors);
       }
//  public void plotScanTable(double[][][][][] scanData, double micronsPerSample, double sampleFrequency,  String title, Color[] colors) {
       return;
    } else if (label.equals("Test 3-rd motor")) {
       DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
       testThirdMotorsDialog();
       return;
    } else if (label.equals("Scan M1&M2, getM3")) {
       DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
       testThirdByFirstTwoDialog();
       return;
    } else if (label.equals("Adjust all motors")) {
       DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
       testAdjustAllThreeDialog();
       return;
    }



    DEBUG_LEVEL=MASTER_DEBUG_LEVEL;
/* other commands need current image (float) */
    imp_src = WindowManager.getCurrentImage();
    if (imp_src==null) {
      IJ.showStatus("No image");
      IJ.showMessage("Focus Tuning Error","Image required");
      return;
    }
    if (imp_src.getType()!=ImagePlus.GRAY32) {
      IJ.showMessage("Focus Tuning Error","Not a 32-bit grayscale image!");
//      new ImageConverter(imp_src).convertToGray32();
      return;
    }
    imp_src.getProcessor();
    String newTitle= imp_src.getTitle();
    Rectangle r=new Rectangle(imp_src.getWidth(),imp_src.getHeight());

    if (label.equals("Split Bayer")) {
      input_bayer=splitBayer (imp_src);
//      if (DEBUG_LEVEL>5) showBayers(input_bayer, r.width>>1, r.height>>1, newTitle);
      showBayers(input_bayer, r.width>>1, r.height>>1, newTitle);
      return;
    } else if (label.equals("Test")) {
      float[] max2= measuresOne(input_bayer[0], FFTSize, Hamming, TargetFilter,test_x>>1, test_y>>1, r.width>>1, true, newTitle);
      if (DEBUG_LEVEL>5) IJ.showMessage("test_levels","max2[0]="+max2[0]+"\nmax2[1]="+max2[1]);
      return;
    } else if (label.equals("Find Targets")) {
      input_bayer=splitBayer (imp_src);
      if (DEBUG_LEVEL>5) showBayers(input_bayer, r.width>>1, r.height>>1, newTitle);
      float [][][] cvh_pixels=new float[4][][];
          float[][] vh_pixels;
      if ((bayerComponent>=0) && (bayerComponent<4)) {
        IJ.showStatus("Scanning for targets, color="+bayerComponent);
        vh_pixels=scanTargets(input_bayer[bayerComponent], FFTSize, FFTScanStep, Hamming, TargetFilter, r.width>>1, r.height>>1, true, newTitle);
      } else {
        IJ.showStatus("Scanning for targets, color=0");
        cvh_pixels[0]=scanTargets(input_bayer[0], FFTSize, FFTScanStep, Hamming, TargetFilter, r.width>>1, r.height>>1, true, newTitle);
        IJ.showStatus("Scanning for targets, color=1");
        cvh_pixels[1]=scanTargets(input_bayer[1], FFTSize, FFTScanStep, Hamming, TargetFilter, r.width>>1, r.height>>1, true, newTitle);
        IJ.showStatus("Scanning for targets, color=2");
        cvh_pixels[2]=scanTargets(input_bayer[2], FFTSize, FFTScanStep, Hamming, TargetFilter, r.width>>1, r.height>>1, true, newTitle);
        IJ.showStatus("Scanning for targets, color=3");
        cvh_pixels[3]=scanTargets(input_bayer[3], FFTSize, FFTScanStep, Hamming, TargetFilter, r.width>>1, r.height>>1, true, newTitle);
        vh_pixels=useBestComponent(cvh_pixels);
      }

      int mwidth=((r.width>>1)-FFTSize)/FFTScanStep;
      int mheight=((r.height>>1)-FFTSize)/FFTScanStep;

      int eqRadius=equalizationRadius/FFTScanStep/2;
      IJ.showStatus("Equalizing");
      equalize(vh_pixels[0], mwidth, mheight, eqRadius);
      equalize(vh_pixels[1], mwidth, mheight, eqRadius);
      

         // public static int     minTargetSize=   100; // in original pixels, used in minFilter()



      float diff_pixels[]=new float[vh_pixels[0].length];
      for (i=0;i<diff_pixels.length;i++) {
        diff_pixels[i]=vh_pixels[0][i]-vh_pixels[1][i];
      }
      if (DEBUG_LEVEL>2) {
        ImageProcessor ip_v= new FloatProcessor(mwidth,mheight);
        ip_v.setPixels(vh_pixels[0]);
        ImagePlus imp_v=  new ImagePlus(newTitle+"_targets_vertical_byr"+bayerComponent, ip_v);
        imp_v.getProcessor().resetMinAndMax();
        imp_v.show();

        ImageProcessor ip_h= new FloatProcessor(mwidth,mheight);
        ip_h.setPixels(vh_pixels[1]);
        ImagePlus imp_h=  new ImagePlus(newTitle+"_targets_horizontal_byr"+bayerComponent, ip_h);
        imp_h.getProcessor().resetMinAndMax();
        imp_h.show();

        ImageProcessor ip_d= new FloatProcessor(mwidth,mheight);
        ip_d.setPixels(diff_pixels);
        ImagePlus imp_d=  new ImagePlus(newTitle+"_targets_differencel_byr"+bayerComponent, ip_d);
        imp_d.getProcessor().resetMinAndMax();
        imp_d.show();
      }
      int radius=correlationRadius/FFTScanStep/2;
      IJ.showStatus("Correlating Vertical/Horizontal patterns");
      float [] pixels_corr=correlateTargets(vh_pixels, mwidth, mheight, radius);
      if (DEBUG_LEVEL>2) {
        ImageProcessor ip_corr= new FloatProcessor(2*radius+1,2*radius+1);
        ip_corr.setPixels(pixels_corr);
        ImagePlus imp_corr=  new ImagePlus(newTitle+"_targets_correletion_byr"+bayerComponent, ip_corr);
        imp_corr.getProcessor().resetMinAndMax();
        imp_corr.show();
      }
      int xy[]=findMaxPix(pixels_corr, 2*radius+1);
      xy[0]-=radius;
      xy[1]-=radius;
      float [] pixels_product=multiplyTargets(vh_pixels, mwidth, mheight, xy[0], xy[1]);


      if (DEBUG_LEVEL>2) {
        ImageProcessor ip_p= new FloatProcessor(mwidth,mheight);
        ip_p.setPixels(pixels_product);
        ImagePlus imp_p=  new ImagePlus(newTitle+"_targets_product_byr"+bayerComponent, ip_p);
        imp_p.getProcessor().resetMinAndMax();
        imp_p.show();
      }
      int minFilterRadius=(minTargetSize-2*FFTSize)/FFTScanStep/2;
      float [] pixels_product_min=minFilter(pixels_product, mwidth, mheight, minFilterRadius);
      if (DEBUG_LEVEL>1) {
        ImageProcessor ip_pf= new FloatProcessor(mwidth,mheight);
        ip_pf.setPixels(pixels_product_min);
        ImagePlus imp_pf=  new ImagePlus(newTitle+"_targets_product_filt"+bayerComponent, ip_pf);
        imp_pf.getProcessor().resetMinAndMax();
        imp_pf.show();
      }
      targets=locateTargets(pixels_product_min, mwidth, mheight,zonesHor,zonesVert, minFilterRadius ,FFTSize,FFTScanStep);
      if (DEBUG_LEVEL>0) {
        showTargets(targets, targetThreshold, newTitle);
      }
/* make targetLocations*/
      int zv,zh;
      i=0;
      for (zv=0;zv<zonesVert;zv++) for (zh=0;zh<zonesHor;zh++) if (targets[zv][zh][2]>targetThreshold) i++;
      targetLocations=new int [i][4][2];
      targetZones=new int [i][2];

      i=0;
      for (zv=0;zv<zonesVert;zv++) for (zh=0;zh<zonesHor;zh++) if (targets[zv][zh][2]>targetThreshold) {
        targetZones[i][0]=zh;
        targetZones[i][1]=zv;
        targetLocations[i][0][0]=(int)(0.5+targets[zv][zh][1]); // X of vertical lines pattern
        targetLocations[i][0][1]=(int)(0.5+targets[zv][zh][0]); // Y of vertical lines pattern
//      xy[0]*FFTScanStep*2;
//      xy[1]*FFTScanStep*2;

// public static double  calib_bw2hv=0.85; // relative distance to Black and white from H to V
// does not account for H pattern becoming V if turned 90 degrees

        targetLocations[i][1][0]=targetLocations[i][0][0]+xy[0]*FFTScanStep*2; // X of horizontal lines pattern
        targetLocations[i][1][1]=targetLocations[i][0][1]+xy[1]*FFTScanStep*2; // Y of horizontal lines pattern

//        targetLocations[i][2][0]=targetLocations[i][0][0]+calib_dx; // X of white calibration

        targetLocations[i][2][0]=targetLocations[i][0][0]- ((int) (calib_bw2hv*xy[1]*FFTScanStep*2)); // X of white calibration
//        targetLocations[i][2][1]=targetLocations[i][0][1]+calib_dy; // Y of white calibration
        targetLocations[i][2][1]=targetLocations[i][0][1]+ ((int) (calib_bw2hv*xy[0]*FFTScanStep*2)); // Y of white calibration

        targetLocations[i][3][0]=targetLocations[i][2][0]+xy[0]*FFTScanStep*2; // X of black calibration
        targetLocations[i][3][1]=targetLocations[i][2][1]+xy[1]*FFTScanStep*2; // Y of black calibration
        i++;
      }

      if (DEBUG_LEVEL>0) {
         showTargetLocations(targetLocations, newTitle);
      }
      return;
/*    } else if (label.equals("Measure Targets")) {
      measureAndDisplay(imp_src);
      return;
*/
    } else if (label.equals("Run")) {
      run=true;
      while (run) {
        imp_camera=acquireMeasureAndDisplay(imp_camera);
try{
  //do what you want to do before sleeping
  IJ.showStatus("Sleeping");
  Thread.currentThread();
Thread.sleep(5000);//sleep for 1000 ms
  IJ.showStatus("Sleeping Over");
  //do what you want to do after sleeping
}
catch(InterruptedException ie){
//If this thread was interrupted by another thread 
      run=false;
}

      }
      return;
    } else if (label.equals("Stop")) {
      run=false;
      return;
    }

    IJ.showStatus("Focus tuning DONE");
  }
      
  public ImagePlus acquireMeasureAndDisplay(ImagePlus imp_src) {
     ImagePlus imp=jp4_instance.openURL(imp_src);
     if (imp!=null) measureAndDisplay(imp);
     return imp;
  }


  public boolean measureAndDisplay(ImagePlus imp) {
/* Make sure targets are defined (image existence is already tested)*/
    if (targetLocations==null) {
      IJ.showMessage("Error","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return false;
    }
    double [][][] focus_results = measureTargets (targetLocations, imp);
    if (DEBUG_LEVEL>1) showFocusTable(focus_results,  imp.getTitle(), 4 /*precision*/);
    showFocusImage(imp_display, ip_display, focus_results, targetZones, zonesHor, zonesVert, displayRotate, displayZonesGap, displayElementGap , displayMin, displayMax, goalG2R, displayVDivs );
    return true;
  }




  public void processWindowEvent(WindowEvent e) {
    super.processWindowEvent(e);
    if (e.getID()==WindowEvent.WINDOW_CLOSING) {
      instance = null;	
    }
  }
  public boolean showDialog() {
    int i;
    GenericDialog gd = new GenericDialog("Focus tuning parameters");
    gd.addStringField ("Filename prefix:                   ", jp4_instance.getTitle(), 20);
    gd.addNumericField("FFT_Size:",                           FFTSize, 0);
    gd.addNumericField("FFT_Scan_Step:",                      FFTScanStep, 0);
    gd.addNumericField("Target_Frequency Tolerance (%):",     100* targetFTolerance, 1);
    gd.addNumericField("Target_Frequency (cpp):",             targetF, 3);
    gd.addNumericField("Target_Angle   (degr):",              targetAngle, 1);
    gd.addNumericField("Correlation_radius (image pix):",     correlationRadius, 0);
    gd.addNumericField("Equalization_radius (image pix):",    equalizationRadius, 0);
    gd.addNumericField("Minimal_Target_Size (image pix):",    minTargetSize, 0);
    gd.addNumericField("Target_Zones_(vertical):",            zonesVert, 0);
    gd.addNumericField("Target_Zones_(horizontal):",          zonesHor, 0);
    gd.addNumericField("Target_Threshold (%):",               100*targetThreshold, 1);
//    gd.addNumericField("Pattern_to_calibration_(horizontal):",calib_dx, 0);
//    gd.addNumericField("Pattern_to_calibration_(vertical):",  calib_dy, 0);
    gd.addNumericField("Pattern distance to B&W centers relative to H to V centers:",  calib_bw2hv, 4); //0.85; relative distance to Black and white from H to V

    gd.addNumericField("Calibration_sample_size:",            calib_sample, 0);

    gd.addNumericField("Output_display_width  (pixels):",     displayWidth, 0);
    gd.addNumericField("Output_display_height (pixels):",     displayHeight, 0);
    gd.addNumericField("Output_display_rotation (0..4 CCW)):",displayRotate, 0);
    gd.addNumericField("Display_zones_gap (pix):",            displayZonesGap,0);
    gd.addNumericField("Display_elements_gap (pix):",         displayElementGap,0);
    gd.addNumericField("Display_focus_min:",                  displayMin,1);
    gd.addNumericField("Display_focus_max:",                  displayMax,1);
    gd.addNumericField("Level of G/R when color changes:",    goalG2R,3);

    gd.addNumericField("Display_vertical_divisions:",         displayVDivs,0);

    gd.addNumericField("Bayer component:",    bayerComponent, 0);
    gd.addNumericField("test_x:",             test_x, 0);
    gd.addNumericField("test_y:",             test_y, 0);

    gd.addNumericField("Debug Level:",      MASTER_DEBUG_LEVEL, 0);
    gd.showDialog();
    if (gd.wasCanceled()) return false;
    jp4_instance.setTitle(gd.getNextString());

    FFTSize=1;
    for (i=(int) gd.getNextNumber(); i >1; i>>=1) FFTSize <<=1; /* make FFTSize to be power of 2*/
    FFTScanStep=        (int) gd.getNextNumber();

    targetFTolerance=   0.01* gd.getNextNumber();
    targetF=                  gd.getNextNumber();
    targetAngle=              gd.getNextNumber();
    correlationRadius=  (int) gd.getNextNumber();
    equalizationRadius= (int) gd.getNextNumber();
    minTargetSize=      (int) gd.getNextNumber();
    zonesVert=          (int) gd.getNextNumber();
    zonesHor=           (int) gd.getNextNumber();
    targetThreshold=    0.01* gd.getNextNumber();
//    calib_dx=           (int) gd.getNextNumber();
//    calib_dy=           (int) gd.getNextNumber();
    calib_bw2hv=              gd.getNextNumber();
    calib_sample=       (int) gd.getNextNumber();
    displayWidth=       (int) gd.getNextNumber();
    displayHeight=      (int) gd.getNextNumber();
    displayRotate=      (int) gd.getNextNumber();

    displayZonesGap=    (int) gd.getNextNumber();
    displayElementGap=  (int) gd.getNextNumber();
    displayMin=               gd.getNextNumber();
    displayMax=               gd.getNextNumber();
    goalG2R=                  gd.getNextNumber();

    displayVDivs=       (int) gd.getNextNumber();



    bayerComponent=     (int) gd.getNextNumber();
    test_x=             (int) gd.getNextNumber();
    test_y=             (int) gd.getNextNumber();
    MASTER_DEBUG_LEVEL= (int) gd.getNextNumber();
    return true;
  }


  public boolean showMotorsDialog() {
    GenericDialog gd = new GenericDialog("Focus tuning parameters");
    gd.addStringField ("Motor Driver IP Address:                   ", motorsIP, 20);
    gd.addNumericField("Motor wait after move before image (sec):  ", motorsWaitAfterMove, 2);
    gd.addNumericField("Motor hysteresis (steps):                  ", motors_hysteresis, 0);
    gd.addNumericField("Motor focus scan range (form minus to plus)", motors_scan_range, 0);
    gd.addNumericField("Motor focus scan steps                     ", motors_scan_steps, 0);
    gd.addNumericField("Motor focus fast scan steps                ", motors_scan_steps_fast, 0);
    gd.addNumericField("Motor scan average number (running average)", motors_scan_averageNum, 0); //5


    gd.addNumericField("Motor 1 relative scan (-1.0.. 1.0)         ", scanScales[0], 4);
    gd.addNumericField("Motor 2 relative scan (-1.0.. 1.0)         ", scanScales[1], 4);
    gd.addNumericField("Motor 3 relative scan (-1.0.. 1.0)         ", scanScales[2], 4);

    gd.addNumericField("Motor 1 min value (steps)                  ", motorsRange[0][0], 0);
    gd.addNumericField("Motor 1 max value (steps)                  ", motorsRange[0][1], 0);
    gd.addNumericField("Motor 2 min value (steps)                  ", motorsRange[1][0], 0);
    gd.addNumericField("Motor 2 max value (steps)                  ", motorsRange[1][1], 0);
    gd.addNumericField("Motor 3 min value (steps)                  ", motorsRange[2][0], 0);
    gd.addNumericField("Motor 3 max value (steps)                  ", motorsRange[2][1], 0);
    gd.addNumericField("Average measurements at Refine and Scan:   ", motor3rd_refineAverage,     0); // 2
    gd.addCheckbox(    "Scan both directions?:                     ", motorScanBothDirections);

    gd.addNumericField("Number of motors steps (all 3) per 1 micron", stepsPerUM, 1); // number of motor steps per 1 um (all 3 motors the same)
    gd.addNumericField("Adjustment goal number (1:R*R*G*G)         ", goalNumber,0);    //  1 - number of optimizing goal (1 - Rv*Rh*Gv*Gh, 2 - Rv*Rh, 3 - Gv*Gh)
    gd.addNumericField("Goal green weight (0..1.0)                 ", goalGreenWeight, 4);

    gd.addNumericField("Focus offset (+ move focus far, - near)    ", focusOffset,2);    //  9.0 - um goal offseet from the best focus at target distance (9um will be for infinity)
    gd.addNumericField("Half window for maximum goal fn search     ", motors_calibGoalWindow,0); //  5; half-window width when looking for maximum in a goal function (first - just maximum, then - mpy by huffaman and find centroid)
    gd.addNumericField("Max to slopes ratio (for maximum search)   ", motors_calibMaxSlopes,3); //  5; half-window width when looking for maximum in a goal function (first - just maximum, then - mpy by huffaman and find centroid)
    gd.addNumericField("Relative difference of measurement to be bad ", badMeasurementError, 4);
    gd.addNumericField("Conservative coefficient for measuremts RMS ", conservativeCoefficient, 4);
    gd.addNumericField("Unreliable from ends (search for maximums)  ",unreliableFromEnds, 1);


    gd.showDialog();
    if (gd.wasCanceled()) return false;
    motorsIP = gd.getNextString();
    motorsWaitAfterMove=       gd.getNextNumber();
    motors_hysteresis=   (int) gd.getNextNumber();
    motors_scan_range=   (int) gd.getNextNumber();
    motors_scan_steps=   (int) gd.getNextNumber();
    motors_scan_steps_fast=   (int) gd.getNextNumber();
    motors_scan_averageNum=   (int) gd.getNextNumber();
    scanScales[0]=gd.getNextNumber();
    scanScales[1]=gd.getNextNumber();
    scanScales[2]=gd.getNextNumber();
    motorsRange[0][0]= (int) gd.getNextNumber();
    motorsRange[0][1]= (int) gd.getNextNumber();
    motorsRange[1][0]= (int) gd.getNextNumber();
    motorsRange[1][1]= (int) gd.getNextNumber();
    motorsRange[2][0]= (int) gd.getNextNumber();
    motorsRange[2][1]= (int) gd.getNextNumber();
    motor3rd_refineAverage= (int) gd.getNextNumber();
    motorScanBothDirections=gd.getNextBoolean();
    stepsPerUM=              gd.getNextNumber();
    goalNumber=        (int) gd.getNextNumber();
    goalGreenWeight=         gd.getNextNumber();

    focusOffset=             gd.getNextNumber();
    motors_calibGoalWindow=  (int) gd.getNextNumber();
    motors_calibMaxSlopes=   gd.getNextNumber();
    badMeasurementError=     gd.getNextNumber();
    conservativeCoefficient= gd.getNextNumber();
    unreliableFromEnds=      gd.getNextNumber();

    return true;
  }

  public boolean moveMotorsDialog() {
    int i;
    int d;
    double new_motor3rd_T2R;
    double new_motors_calibMaxSlopes;
    double newGoalGreenWeight;
    GenericDialog gd = new GenericDialog("Focus tuning parameters");
/// Read current motor positions, use them as defaults
    motorPositions =readElphel10364Motors(motorsIP);

//    gd.addStringField ("Motor Driver IP Address:                   ", motorsIP, 20);
    gd.addNumericField("Motor 1 (steps):                  ", motorPositions[0], 0);
    gd.addNumericField("Motor 2 (steps):                  ", motorPositions[1], 0);
    gd.addNumericField("Motor 3 (steps):                  ", motorPositions[2], 0);
    if (targetScanData != null) {
      gd.addNumericField("Target Tangential/Radial resolution", motor3rd_T2R,              4); // 1.4
      gd.addCheckbox(    "Measure focus distances?          ", motors_measureDisplacements);
      gd.addNumericField("Try sideways from min difference: ", motors_sideWays, 0);
      gd.addNumericField("Average measurements:   ",          motor3rd_refineAverage,     0); // 2
      gd.addNumericField("Number of motors steps (all 3) per 1 micron", stepsPerUM, 1); // number of motor steps per 1 um (all 3 motors the same)
      gd.addNumericField("Adjustment goal number (1:R*R*G*G)         ", goalNumber,0);    //  1 - number of optimizing goal (1 - Rv*Rh*Gv*Gh, 2 - Rv*Rh, 3 - Gv*Gh)
      gd.addNumericField("Goal green weight (0..1.0)                 ", goalGreenWeight, 4);
      gd.addNumericField("Focus offset (+ move focus far, - near)    ", focusOffset,2);    //  9.0 - um goal offseet from the best focus at target distance (9um will be for infinity)
      gd.addNumericField("Half window for maximum goal fn search     ", motors_calibGoalWindow,0); //  5; half-window width when looking for maximum in a goal function (first - just maximum, then - mpy by huffaman and find centroid)
      gd.addNumericField("Max to slopes ratio (for maximum search)   ", motors_calibMaxSlopes,3); //  5; half-window width when looking for maximum in a goal function (first - just maximum, then - mpy by huffaman and find centroid)
    }
// public static double   goalGreenWeight=   0.7;    // Goal (for numbers 4,8,12) will be 0.7 from goalGreenWeight:(1-goalGreenWeight) between green and red

    gd.showDialog();
    if (gd.wasCanceled()) return false;
    motorTargetPositions[0]= (int) gd.getNextNumber();
    motorTargetPositions[1]= (int) gd.getNextNumber();
    motorTargetPositions[2]= (int) gd.getNextNumber();
    if (targetScanData != null) {
      new_motor3rd_T2R=            gd.getNextNumber();
      motors_measureDisplacements=gd.getNextBoolean();
      motors_sideWays=         (int) gd.getNextNumber();
      motor3rd_refineAverage=  (int) gd.getNextNumber();
      stepsPerUM=                    gd.getNextNumber();
      goalNumber=              (int) gd.getNextNumber();
      newGoalGreenWeight=            gd.getNextNumber();
      focusOffset=                   gd.getNextNumber();
      d=                       (int) gd.getNextNumber();
      new_motors_calibMaxSlopes=     gd.getNextNumber();
      if ((d!=motors_calibGoalWindow) ||
          (new_motor3rd_T2R!=motor3rd_T2R) ||
          (new_motors_calibMaxSlopes!=motors_calibMaxSlopes) ||
          (newGoalGreenWeight !=goalGreenWeight)) { // recalculate if the window had chnaged
         goalGreenWeight=newGoalGreenWeight;
         motor3rd_T2R=           new_motor3rd_T2R;
         motors_calibGoalWindow= d;
         motors_calibMaxSlopes=  new_motors_calibMaxSlopes;
         maxScanTargetOrientColor=calcMaxScanTargetOrientColor(scanResults, motors_calibGoalWindow, motors_calibMaxSlopes); // recalculate maximums 
         showScanTable(scanResults,  startScan, endScan, motors_scan_steps, maxScanTargetOrientColor, "Recalculated Scan_Focus_Results" , 4); /// 4 - precision
// Is it needed here?
         showAstigmatismTable(maxScanTargetOrientColor, //[scan direction][target][v,h][r,g,b] 1.0 - one scan sample
                                   motors_calibSteps/stepsPerUM,  // microns per 1 scan sample
                                   motors_scan_steps,  // scan range (samples, starting from 0)
                                   2.0, // maximums outside of [reliabilityThreshold..range-reliabilityThreshold are
                                                                // considered unreliable
                                   ((imp_src==null)?"":imp_src.getTitle())+"_astigmatism_table",                // window title
                                   2);   // output precision in decimals

         offsetsPerGoal= makeOffsetsPerGoal(scanResults, motor3rd_targets_T2R,motors_calibGoalWindow,motors_calibMaxSlopes,goalGreenWeight);

      }
    }

    motorPositions =moveElphel10364Motors(motorsIP,motorTargetPositions,1.0,true,"manual", true);

    motorPositions =readElphel10364Motors(motorsIP); /// re-read positions to make sure motors did settle down
    if (motors_measureDisplacements && (targetScanData != null)) {
      double[][][] resolutions= measureTargets (targetLocations, imp_camera, motor3rd_refineAverage);
      if (goalNumber <1 ) goalNumber=1;
      if (goalNumber > offsetsPerGoal.length ) goalNumber=offsetsPerGoal.length;

      double [] displacements=findDisplacemntsFromResolutions (targetScanData, resolutions, motors_sideWays, motors_calibSteps,stepsPerUM, offsetsPerGoal[goalNumber-1],focusOffset);
      String srslt= "M1="+motorTargetPositions[0]+"\n"+
                    "M2="+motorTargetPositions[1]+"\n"+
                    "M3="+motorTargetPositions[2]+"\n"+
                    " --- Target focus distances ---\n";
      for (i=0;i < displacements.length; i++)  srslt+=(i+1)+": "+displacements[i]+" um\n";
      if (displacements.length==5) {
        srslt+=      " --- Position/orientation Errors ---\n";
        srslt+="Center: "+(displacements[2]                 )+" um\n";
        srslt+="   1-5: "+(displacements[0]-displacements[4])+" um\n";
        srslt+="   2-4: "+(displacements[1]-displacements[3])+" um\n";
      }
      IJ.showMessage("Distances from focus",srslt);
    }
    return true;
  }
  public boolean testThirdMotorsDialog() {
    int d,i;
    GenericDialog gd = new GenericDialog("Test adjust third motor");
    gd.addNumericField("Signed hysteresis:                ", motor3rd_signed_hysteresis, 0); // 100
    gd.addNumericField("motorNumber (1..3):               ", motor3rd_motorNumber,       0); // 0
    gd.addNumericField("targetNumber (1..5):              ", motor3rd_targetNumber,      0); // 2
    gd.addNumericField("minStep:                          ", motor3rd_minStep,           0); // 30
    gd.addNumericField("maxStep:                          ", motor3rd_maxStep,           0); // 300
    gd.addNumericField("Direction (1,5:-1, 2,4:+1):       ", motor3rd_direction,         0); // 1
    gd.addNumericField("color mask (2-R,4-B,9-G, 0-G-R):  ", motor3rd_colorMask,         0); // 2
    gd.addNumericField("Target V/H (positive Direction):  ", motor3rd_V2H,               4); // 1.0
    gd.addNumericField("Target Tangential/Radial resolution", motor3rd_T2R,              4); // 1.4
    gd.addNumericField("Refine final steps (number):      ", motor3rd_refine,            0); // 6
    gd.addNumericField("Average measurements at Refine:   ", motor3rd_refineAverage,     0); // 2
    gd.addNumericField("Expand refine range if fails:     ", motor3rd_refineExpand,      0); // 4
    gd.addNumericField("Maximal iteractions:              ", motor3rd_maxIterations,     0); // 20
    gd.addNumericField("Motor wait after move before image (sec):  ", motorsWaitAfterMove, 2);
    if (targetScanData != null) {
      gd.addCheckbox(    "Use calibration data?:            ", motor3rd_useCalibration);
      gd.addNumericField("Minimal step to re-calibrate motor", minUpdateStep,            0);    // 100
      gd.addNumericField("Maximal motor scale correction:   ", variableStepRatio,        1); // 1.0
      gd.addNumericField("Adjustment goal number (1:R*R*G*G)         ", goalNumber,0);    //  1 - number of optimizing goal (1 - Rv*Rh*Gv*Gh, 2 - Rv*Rh, 3 - Gv*Gh)
      gd.addNumericField("Focus offset (+ move focus far, - near)    ", focusOffset,2);    //  9.0 - um goal offseet from the best focus at target distance (9um will be for infinity)

    }

// public static int       motor3rd_refineAverage=      2; // number of measurements during refine stage
// public static int       motor3rd_refineExpand=       4; // expand refine range if fails

    gd.showDialog();
    if (gd.wasCanceled()) return false;
    motor3rd_signed_hysteresis= (int) gd.getNextNumber();
    d=(int) gd.getNextNumber(); if ((d>0) && (d<4))  motor3rd_motorNumber=d;
    d=(int) gd.getNextNumber(); if ((d>0) && (d<6))  motor3rd_targetNumber=d;
    motor3rd_minStep=       (int) gd.getNextNumber();
    motor3rd_maxStep=       (int) gd.getNextNumber();
    motor3rd_direction=   (((int) gd.getNextNumber())>=0)?1:-1;
    motor3rd_colorMask=     (int) gd.getNextNumber();
    motor3rd_V2H=                 gd.getNextNumber();
    motor3rd_T2R=                 gd.getNextNumber();
    motor3rd_refine=        (int) gd.getNextNumber();

    motor3rd_refineAverage= (int) gd.getNextNumber();
    motor3rd_refineExpand=  (int) gd.getNextNumber();

    motor3rd_maxIterations= (int) gd.getNextNumber();
    motorsWaitAfterMove=          gd.getNextNumber();
    if (targetScanData != null) {
       motor3rd_useCalibration=         gd.getNextBoolean();
       minUpdateStep =            (int) gd.getNextNumber();
       variableStepRatio=               gd.getNextNumber();
       goalNumber=                (int) gd.getNextNumber();
       focusOffset=                     gd.getNextNumber();
    }
    if (motor3rd_useCalibration) {
       findThirdMotorCalib(  motorsIP,
                             motor3rd_signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                             motorsWaitAfterMove,
                             motor3rd_motorNumber-1,
                             motor3rd_minStep,           // minimal motor rotation - positive
                             motor3rd_targetNumber-1,
                             motor3rd_refineAverage,     // number of measurements to average
                             minUpdateStep,              // ipdate thisStepsPerUM only if the last step was > this threshold
                             variableStepRatio,          // 2.0 - maximal thisStepsPerUM change in one step
                             motor3rd_maxIterations);    // - just for debugging - maximal number of steps to try before giving up
      if (motors_measureDisplacements && (targetScanData != null)) {
        double[][][] resolutions= measureTargets (targetLocations, imp_camera, motor3rd_refineAverage);
        if (goalNumber <1 ) goalNumber=1;
        if (goalNumber > offsetsPerGoal.length ) goalNumber=offsetsPerGoal.length;
        double [] displacements=findDisplacemntsFromResolutions (targetScanData, resolutions, motors_sideWays, motors_calibSteps,stepsPerUM, offsetsPerGoal[goalNumber-1],focusOffset);
        int [] motorPositions =readElphel10364Motors(motorsIP); /// re-read positions to make sure motors did settle down

        String srslt= "M1="+motorPositions[0]+"\n"+
                      "M2="+motorPositions[1]+"\n"+
                      "M3="+motorPositions[2]+"\n"+
                      " --- Target focus distances ---\n";
        for (i=0;i < displacements.length; i++)  srslt+=(i+1)+": "+displacements[i]+" um\n";
        if (displacements.length==5) {
          srslt+=      " --- Position/orientation Errors ---\n";
          srslt+="Center: "+(displacements[2]                 )+" um\n";
          srslt+="   1-5: "+(displacements[0]-displacements[4])+" um\n";
          srslt+="   2-4: "+(displacements[1]-displacements[3])+" um\n";
        }
        IJ.showMessage("Distances from focus",srslt);
      }
    } else {
      findThirdMotor(imp_camera,
                     motorsIP,
                     motor3rd_signed_hysteresis,
                     motorsWaitAfterMove,
                     motor3rd_motorNumber-1,
                     motor3rd_minStep,
                     motor3rd_maxStep,
                     motor3rd_targetNumber-1,
                     motor3rd_colorMask,
                     motor3rd_V2H* motor3rd_direction,
                     motor3rd_refine,
                     motor3rd_refineAverage,
                     motor3rd_refineExpand,
                     motor3rd_maxIterations);
     }
     return true;
  }

  public boolean testThirdByFirstTwoDialog() {
    int i;
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return false;
    }

    GenericDialog gd = new GenericDialog("Test adjust third motor");
    gd.addNumericField("M1 scan range:                     ", motor3rd_M1range,          0);
    gd.addNumericField("M1 scan steps:                     ", motor3rd_M1steps,          0);
    gd.addNumericField("M2 scan range:                     ", motor3rd_M2range,          0);
    gd.addNumericField("M2 scan steps:                     ", motor3rd_M2steps,          0);
    gd.addNumericField("Target Tangential/Radial resolution", motor3rd_T2R,              4); // 1.4

    for (i=0;i<motor3rd_targets_dirs.length; i++) {
      gd.addNumericField("Target "+(i+1)+" Direction (0- disable target): ", motor3rd_targets_dirs[i],       0);
      gd.addNumericField("Target "+(i+1)+" color mask:                    ", motor3rd_targets_ColorMasks[i], 0);
      gd.addNumericField("Target "+(i+1)+" tang/radial resolution:        ", motor3rd_targets_T2R[i],        4);
    }
    gd.addNumericField("Signed hysteresis:                ", motor3rd_signed_hysteresis, 0); // 100
    gd.addNumericField("minStep:                          ", motor3rd_minStep,           0); // 30
    gd.addNumericField("maxStep:                          ", motor3rd_maxStep,           0); // 300
    gd.addNumericField("Refine final steps (number):      ", motor3rd_refine,            0); // 6
    gd.addNumericField("Average measurements at Refine:   ", motor3rd_refineAverage,     0); // 2
    gd.addNumericField("Expand refine range if fails:     ", motor3rd_refineExpand,      0); // 4
    gd.addNumericField("Maximal iterations :              ", motor3rd_maxIterations,     0); // 20
    gd.addNumericField("Motor wait after move before image (sec):  ", motorsWaitAfterMove, 2);
    gd.showDialog();
    if (gd.wasCanceled()) return false;

    motor3rd_M1range= (int) gd.getNextNumber();
    motor3rd_M1steps= (int) gd.getNextNumber();
    motor3rd_M2range= (int) gd.getNextNumber();
    motor3rd_M2steps= (int) gd.getNextNumber();
    motor3rd_T2R=           gd.getNextNumber();
    for (i=0;i<motor3rd_targets_dirs.length; i++) {
      motor3rd_targets_dirs[i]=       (int) gd.getNextNumber();
      motor3rd_targets_ColorMasks[i]= (int) gd.getNextNumber();
//      motor3rd_targets_T2R[i]=              gd.getNextNumber();
    }
    motor3rd_targets_T2R[0]=    motor3rd_T2R; motor3rd_targets_T2R[4]=    motor3rd_T2R;
    motor3rd_targets_T2R[1]=1.0/motor3rd_T2R; motor3rd_targets_T2R[3]=1.0/motor3rd_T2R;
    motor3rd_targets_T2R[2]=1.0;

    motor3rd_signed_hysteresis= (int) gd.getNextNumber();
    motor3rd_minStep=       (int) gd.getNextNumber();
    motor3rd_maxStep=       (int) gd.getNextNumber();
    motor3rd_refine=        (int) gd.getNextNumber();
    motor3rd_refineAverage= (int) gd.getNextNumber();
    motor3rd_refineExpand=  (int) gd.getNextNumber();
    motor3rd_maxIterations= (int) gd.getNextNumber();
    motorsWaitAfterMove=          gd.getNextNumber();
    int [] motor_numbers={0,1,2};
    int [][][] motorsTable; //[M1_index][M2_index][M1,M2,M3[target]]
    double [] V2H=new double[motor3rd_targets_T2R.length];



    for (i=0;i<motor3rd_targets_T2R.length;i++) V2H[i]=motor3rd_targets_T2R[i]*((motor3rd_targets_dirs[i]>0)?1:((motor3rd_targets_dirs[i]<0)?-1:0 ));
//    thirdByFirstTwo(ij.ImagePlus,java.lang.String,int,double,int,  int,int,int,int,int[],int[],double[],int,int,double[][][],int,int)
//                to (ij.ImagePlus,java.lang.String,int,double,int[],int,int,int,int,int[],int[],double[],int,int,<nulltype>,int,int)
    motorsTable=thirdByFirstTwo(
                    imp_camera,
                    motorsIP,
                    motor3rd_signed_hysteresis,
                    motorsWaitAfterMove,
                    motor_numbers,

                    motor3rd_M1range,
                    motor3rd_M1steps,
                    motor3rd_M2range,
                    motor3rd_M2steps,

                    motor3rd_targets_dirs,
                    motor3rd_targets_ColorMasks,
                    V2H,
 
                    motor3rd_minStep,
                    motor3rd_maxStep,
//                    null,
                    motor3rd_refine,
                    motor3rd_refineAverage,
                    motor3rd_refineExpand,
                    motor3rd_maxIterations);

     showThirdByFirstTwoTargets(motorsTable,
                        motor_numbers,
                        motor3rd_targets_dirs,
                        "Motors Table "+motor3rd_M1range+"/"+ motor3rd_M1steps+" "+motor3rd_M2range+"/"+ motor3rd_M2steps);
     showThirdByFirstTwo(motorsTable,
                        motor_numbers,
                        motor3rd_targets_dirs,
                        "Motors Table "+motor3rd_M1range+"/"+ motor3rd_M1steps+" "+motor3rd_M2range+"/"+ motor3rd_M2steps);
     return true;
  }
  public boolean testAdjustAllThreeDialog() {
    int i;
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return false;
    }

    GenericDialog gd = new GenericDialog("Test adjust third motor");
    gd.addNumericField("Initial Steps for M1,M2:           ", motor3rd_initialStep,      0);
    gd.addNumericField("Number of Iterations:              ", motor3rd_numIterations,    0);
    gd.addNumericField("Target Tangential/Radial resolution", motor3rd_T2R,              4); // 1.4

    for (i=0;i<motor3rd_targets_dirs.length; i++) {
      gd.addNumericField("Target "+(i+1)+" Direction (0- disable target): ", motor3rd_targets_dirs[i],       0);
      gd.addNumericField("Target "+(i+1)+" color mask:                    ", motor3rd_targets_ColorMasks[i], 0);
//      gd.addNumericField("Target "+(i+1)+" tang/radial resolution:        ", motor3rd_targets_T2R[i],        4);
    }
    gd.addNumericField("Signed hysteresis:                ", motor3rd_signed_hysteresis, 0); // 100
    gd.addNumericField("minStep:                          ", motor3rd_minStep,           0); // 30
    gd.addNumericField("maxStep:                          ", motor3rd_maxStep,           0); // 300
    gd.addNumericField("Refine final steps (number):      ", motor3rd_refine,            0); // 6
    gd.addNumericField("Average measurements at Refine:   ", motor3rd_refineAverage,     0); // 2
    gd.addNumericField("Expand refine range if fails:     ", motor3rd_refineExpand,      0); // 4
    gd.addNumericField("Maximal iterations:               ", motor3rd_maxIterations,     0); // 20
    gd.addNumericField("Fast iterations (first stage):    ", motor3rd_fastIterations,     0); // 20
    gd.addNumericField("Motor wait after move before image (sec):  ", motorsWaitAfterMove, 2);
    if (targetScanData != null) {
      gd.addCheckbox(    "Use calibration data?:            ", motor3rd_useCalibration);
      gd.addNumericField("Minimal step to re-calibrate motor", minUpdateStep,            0);    // 100
      gd.addNumericField("Maximal motor scale correction:   ", variableStepRatio,        1); // 1.0
      gd.addNumericField("Adjustment goal number (1:R*R*G*G)         ", goalNumber,0);    //  1 - number of optimizing goal (1 - Rv*Rh*Gv*Gh, 2 - Rv*Rh, 3 - Gv*Gh)
      gd.addNumericField("Focus offset (+ move focus far, - near)    ", focusOffset,2);    //  9.0 - um goal offseet from the best focus at target distance (9um will be for infinity)

    }
    gd.showDialog();
    if (gd.wasCanceled()) return false;

    motor3rd_initialStep=   (int) gd.getNextNumber();
    motor3rd_numIterations= (int) gd.getNextNumber();
    motor3rd_T2R=                 gd.getNextNumber();
    for (i=0;i<motor3rd_targets_dirs.length; i++) {
      motor3rd_targets_dirs[i]=       (int) gd.getNextNumber();
      motor3rd_targets_ColorMasks[i]= (int) gd.getNextNumber();
//      motor3rd_targets_T2R[i]=              gd.getNextNumber();
    }
    motor3rd_targets_T2R[0]=    motor3rd_T2R; motor3rd_targets_T2R[4]=    motor3rd_T2R;
    motor3rd_targets_T2R[1]=1.0/motor3rd_T2R; motor3rd_targets_T2R[3]=1.0/motor3rd_T2R;
    motor3rd_targets_T2R[2]=1.0;

    motor3rd_signed_hysteresis= (int) gd.getNextNumber();
    motor3rd_minStep=       (int) gd.getNextNumber();
    motor3rd_maxStep=       (int) gd.getNextNumber();
    motor3rd_refine=        (int) gd.getNextNumber();
    motor3rd_refineAverage= (int) gd.getNextNumber();
    motor3rd_refineExpand=  (int) gd.getNextNumber();
    motor3rd_maxIterations= (int) gd.getNextNumber();
    motor3rd_fastIterations=(int) gd.getNextNumber();
    motorsWaitAfterMove=          gd.getNextNumber();
    if (targetScanData != null) {
       motor3rd_useCalibration=         gd.getNextBoolean();
       minUpdateStep =            (int) gd.getNextNumber();
       variableStepRatio=               gd.getNextNumber();
       goalNumber=                (int) gd.getNextNumber();
       focusOffset=                     gd.getNextNumber();
    }


    int [] motor_numbers={0,1,2};
//    int [][][] motorsTable; //[M1_index][M2_index][M1,M2,M3[target]]
    double [] V2H=new double[motor3rd_targets_T2R.length];

    if (motor3rd_useCalibration) {
      adjustAllThreeCalib( motorsIP,
                      motor3rd_signed_hysteresis,
                      motorsWaitAfterMove,
                      motor_numbers,
                      motor3rd_initialStep,
                      motor3rd_numIterations,
                      motor3rd_minStep,
                      motor3rd_refineAverage,
                      minUpdateStep,              // ipdate thisStepsPerUM only if the last step was > this threshold
                      variableStepRatio,          // 2.0 - maximal thisStepsPerUM change in one step
                      motor3rd_fastIterations,
                      motor3rd_maxIterations);
    } else {
      for (i=0;i<motor3rd_targets_T2R.length;i++) V2H[i]=motor3rd_targets_T2R[i]*((motor3rd_targets_dirs[i]>0)?1:((motor3rd_targets_dirs[i]<0)?-1:0 ));
      adjustAllThree( motorsIP,
                      motor3rd_signed_hysteresis,
                      motorsWaitAfterMove,
                      motor_numbers,

                      motor3rd_initialStep,
                      motor3rd_numIterations,

                      motor3rd_targets_dirs,
                      motor3rd_targets_ColorMasks,
                      V2H,
 
                      motor3rd_minStep,
                      motor3rd_maxStep,
                      motor3rd_refine,
                      motor3rd_refineAverage,
                      motor3rd_refineExpand,
                      motor3rd_maxIterations);
     }
     return true;

  }





  public void initDisplay() {
    if ((imp_display==null) || (imp_display.getWidth()!=displayWidth) || (imp_display.getHeight()!=displayHeight)) {
      if (imp_display!=null)     imp_display.close();
      ip_display=   new ColorProcessor (displayWidth,displayHeight);
      imp_display=  new ImagePlus("Focus Monitor", ip_display);
//      imp_display.getProcessor().resetMinAndMax();
      imp_display.show();

    }
  }

//showFocusImage(focus_results, targetLocations, zonesHor, zonesVert, displayZonesGap, displayElementGap , displayMin, displayMax, displayVDivs );
  public void showFocusImage(ImagePlus imp,  ImageProcessor ip,  double[][][] focusData, int[][]focusZones, int zonesHor, int zonesVert, int rotate, int zonesGap, int elemGap,
                             double min, double max, double g2r, int vdivs ) {
    if (imp==null) return;
    if (focusZones.length==0) return;
//    ImageProcessor ip=imp.getProcessor();
//     ImageProcessor ip=ip_display;
    int iwidth=  imp.getWidth();
    int iheight= imp.getHeight();
    int width;
    int height;
    double [][] bestFocus=new double [2][4];
    int i,x,y,zh,zv;
    int dir,col;
    for (dir=0;dir<2;dir++) for (col=0;col<4;col++) bestFocus[dir][col]=focusData[0][dir][col];
    for (i=1;i<focusZones.length;i++) for (dir=0;dir<2;dir++) for (col=0;col<4;col++) {
      if (bestFocus[dir][col]<focusData[i][dir][col]) bestFocus[dir][col] = focusData[i][dir][col];
    }
    
    if ((rotate &1 )==0) {
      width=  (iwidth- (zonesHor- 1)*zonesGap)/zonesHor;
      height= (iheight-(zonesVert-1)*zonesGap)/zonesVert;
    } else {
      width=  (iwidth- (zonesVert- 1)*zonesGap)/zonesVert;
      height= (iheight-(zonesHor-  1)*zonesGap)/zonesHor;
    }
    for (i=0;i<focusZones.length;i++) {
      zh=   focusZones[i][0];
      zv=   focusZones[i][1];
      switch (rotate) {
//       case 0: zh=            focusZones[i][0]; zv=            focusZones[i][1];  break;
       case 1: zh=            focusZones[i][1]; zv=zonesHor-1- focusZones[i][0];  break;
       case 2: zh=zonesHor-1- focusZones[i][0]; zv=zonesVert-1-focusZones[i][1];  break;
       case 3: zh=zonesVert-1-focusZones[i][1]; zv=            focusZones[i][0];  break;
      }
      x=(width+  zonesGap)*zh;
      y=(height+ zonesGap)*zv;
      displayOneFocus(ip, focusData[i], bestFocus, x, y, width, height, min, max, g2r,  vdivs, elemGap);

    }
    imp.updateAndDraw();
    imp.show();
  }

/// =========================================
  public void displayOneFocus(ImageProcessor ip, double[][] data, double[][] best_data, int x0, int y0, int width, int height, double min, double max,  double g2r, int vdivs, int gap) {
      int [] barColors={0xff0000,0xdd0000,0x88ff00,0x88dd00,0xff88,0xdd88,0xff,0xee};
      double [][] logData=     new double [2][4];
      double [][] logBestData= new double [2][4];
      double [] RminusG =  new double [3];
      float k;
      int dir,col,i,c,x,y,h,h_best;
//      int heye=height/vdivs;
      int heye;
      int rightBarsWidth=width/16-gap;
      int rightBestBarsWidth=(rightBarsWidth+2)/4;
//      Color color=new Color(0xffffff);
      for (dir=0;dir<2;dir++) {
          for (col=0;col<4;col++) {
          logData[dir][col]=    Math.log(Math.max(Math.min(     data[dir][col],max),min)/min)/Math.log(max/min); /* min->0, max ->1.0 */
          logBestData[dir][col]=Math.log(Math.max(Math.min(best_data[dir][col],max),min)/min)/Math.log(max/min); /* min->0, max ->1.0 */
          h=(int) (logData[dir][col]*height);
          h_best=(int) (logBestData[dir][col]*height);
          c=dir+2*(col^1);
          x=x0+(int) (width*(0.5+ 0.0625*c));
//        ip.setColor(color.Color(new Color(barColors[c])));
          ip.setColor(new Color(0));
          ip.setRoi(x,y0,rightBarsWidth,height);
          ip.fill();
          ip.setColor(new Color(barColors[c]));
          ip.setRoi(x,y0+height-h,rightBarsWidth,h);
          ip.fill();
          ip.setRoi(x,y0+height-h_best,rightBestBarsWidth,h_best);
          ip.fill();
        }
        RminusG[dir]= g2r*logData[dir][1]-0.5*(logData[dir][0]+logData[dir][3]);
      }



      RminusG[2]=   0.5*(RminusG[0]+RminusG[1]);
      for (dir=0;dir<3;dir++) {
        for (i=0; i<vdivs; i++) {
          if (i==0) k=(RminusG[dir]>0)?1.0f:0.0f;
          else      k=(float) Math.max(Math.min(0.5+(2.0*vdivs*RminusG[dir])/i,1.0),0.0);
          ip.setColor(new Color(k,1-k,0.0f));
          x=x0+((dir==2)?0:((int) (width*0.125*(dir+2))));
          y=y0+(height*i)/vdivs;
          heye=((height*(i+1))/vdivs)-((height*i)/vdivs);
          ip.setRoi(x,y,((dir==2)?(width/4):(width/8))-gap,heye);
          ip.fill();
        }
      }
  }



/// =========================================

  public void showFocusTable(double[][][] focusData, String title, int precision) {
    String header="#\tVert_Gr\tHor_Gr\tVert_R\tHor_R\tVert_B\tHor_B\tVert_Gb\tHor_Gb";
    StringBuffer sb = new StringBuffer();
    int n;
    for (n=0;n<focusData.length;n++) {
      sb.append((n+1)+
                  "\t"+IJ.d2s(focusData[n][0][0],precision)+ // Vert_Gr
                  "\t"+IJ.d2s(focusData[n][1][0],precision)+ // Hor_Gr
                  "\t"+IJ.d2s(focusData[n][0][1],precision)+ // Vert_R
                  "\t"+IJ.d2s(focusData[n][1][1],precision)+ // Hor_R
                  "\t"+IJ.d2s(focusData[n][0][2],precision)+ // Vert_B
                  "\t"+IJ.d2s(focusData[n][1][2],precision)+ // Hor_B
                  "\t"+IJ.d2s(focusData[n][0][3],precision)+ // Vert_Gb
                  "\t"+IJ.d2s(focusData[n][1][3],precision)+ // Hor_Gb";
                  "\n");
    }
    new TextWindow(title+"_Focus_Table", header, sb.toString(), 600,300);
  }

/// =========================================

  public void showTargetLocations(int[][][] targetLocs, String title) {
    String header="#\tVert_X\tVert_Y\tHor_X\tHor_Y\tWhite_X\tWhite_Y\tBlack_X\tBlack_Y";
    StringBuffer sb = new StringBuffer();
    int n;
    for (n=0;n<targetLocs.length;n++) {
      sb.append((n+1)+
                  "\t"+targetLocs[n][0][0]+
                  "\t"+targetLocs[n][0][1]+
                  "\t"+targetLocs[n][1][0]+
                  "\t"+targetLocs[n][1][1]+
                  "\t"+targetLocs[n][2][0]+
                  "\t"+targetLocs[n][2][1]+
                  "\t"+targetLocs[n][3][0]+
                  "\t"+targetLocs[n][3][1]+"\n");
    }
    new TextWindow(title+"_target_locations", header, sb.toString(), 600,300);
  }


/// =========================================

  public double [][][] locateTargets( float [] pixels, int  width, int height, int zonesHor, int zonesVert, int radius,  int FFTSize, int FFTScanStep) {
    double [][][] result =new double [ zonesVert][ zonesHor][3];
    int zv,zh,x0,y0,x1,y1,x,y, xm, ym, indx;
    double max,p,a,ax,ay; 
    for (zv=0;zv<zonesVert;zv++) for (zh=0;zh<zonesHor;zh++)  {
        y0=(height*zv)/zonesVert;
        y1=(height*(zv+1))/zonesVert;
        x0=(width*zh)/zonesHor;
        x1=(width*(zh+1))/zonesHor;
        ym=y0;
        xm=x0;
        max=pixels[y0*width+x0];
        for (y=y0;y<y1;y++) for (x=x0;x<x1;x++) {
          indx=(y*width+x);
          if (pixels[indx]>max) {
             ym=y;
             xm=x;
             max=pixels[indx];
          }
        }
        y0=ym-radius;
        y1=ym+radius+1;
        x0=xm-radius;
        x1=xm+radius+1;
        if (y0<0) y0=0;
        if (x0<0) x0=0;
        if (y1>height) y1=height;
        if (x1>width)  x1=width;
        a=0.0;
        ax=0.0;
        ay=0.0;
        for (y=y0;y<y1;y++) for (x=x0;x<x1;x++) {
          p=pixels[y*width+x];
          a+=p;
          ax+=p*x;
          ay+=p*y;
        }
        result[zv][zh][0]=(ay/a)*2*FFTScanStep+FFTSize;
        result[zv][zh][1]=(ax/a)*2*FFTScanStep+FFTSize;
        result[zv][zh][2]=max;
    }
    return result;
  }

// imp argument is not used anymore, using global imp_camera
/* measure 'MTF' at target frequencies, for each target, each orientation, each color). Same indexes as in targetLocs */
/**
Tries to re-measure "bad" measurements (i.e. somebody walked through), but false positives  happen when 3 measurements fall very close to each other
so it is difficult to get a replacement image.
TODO: keep running average error, use it fro reference in addition to the other measurements
*/
  public double [][][] measureTargets (int [][][] targetLocs, ImagePlus imp, int N) {
    int T=targetLocs.length;
    double [][][] rslt=new double [targetLocs.length][2][4];
    double [][][][] measAll=new double [N][][][];
    double [][][]   measOne;
    double [][] amplitudes= new double [N][T];
    int t,i,j,n;

    double [] SX_all=     new double [T];
    double [] SX2_all=    new double [T];
    double [] SX_others=  new double [T];
    double [] SX2_others= new double [T];
    double [] rms=        new double [T];
    for (t=0;t<T;t++) {
      SX_all[t]= 0.0;
      SX2_all[t]=0.0;
    }


double ERR_this, ERR_others; 
    boolean allGood=false;
    int retriesLeft;
//    double badMeasurementError2=badMeasurementError*badMeasurementError;
//    System.out.println("measureTargets(), N="+N);
// public double conservativeCoefficient=0.7; // when estimating expected measurement RMS, use conservativeCoefficient of the old one and (1-conservativeCoefficient) of the new one
// public static [] oldErrors;

    for (n=0;n<N;n++) {

      imp_camera=jp4_instance.openURL(imp_camera); // if imp_camera image is closed, causes error
      if (imp_camera==null) {
        IJ.showMessage("Camera Image ERROR","Failed to open camera image");
        return null;
      }

//     System.out.println("measureTargets: n=" + N);
//      measAll[n]=measureTargets (targetLocs, imp);
      measAll[n]=measureTargets (targetLocs, imp_camera);
      for (t=0;t<T;t++) amplitudes[n][t]=0.0;
      for (t=0;t<rslt.length;t++) for (i=0;i<rslt[0].length;i++) for (j=0;j<rslt[0][0].length;j++) amplitudes[n][t]+=measAll[n][t][i][j]*measAll[n][t][i][j];
      for (t=0;t<T;t++)  {
        SX_all[t]+= amplitudes[n][t];
        SX2_all[t]+=amplitudes[n][t]*amplitudes[n][t];
      }
//    System.out.println("--- n="+n+" amplitudes["+n+"]="+amplitudes[n]+" SX_all="+SX_all+", SX2_all=" + SX2_all);
    }
    if (N>2) { // removing bad measurments only works if number of measurements >2
      for(retriesLeft=4;  !allGood && (retriesLeft>0); retriesLeft--) {
        allGood=true;
        for (n=0;n<N;n++) {
          for (t=0;t<T;t++) {
            SX2_others[t]=SX2_all[t]-amplitudes[n][t]*amplitudes[n][t];
            SX_others[t]= SX_all[t]-amplitudes[n][t];
            ERR_others=Math.sqrt(SX2_others[t]/(N-1)  - (SX_others[t]/(N-1))*(SX_others[t]/(N-1)));
// public double conservativeCoefficient=0.7; // when estimating expected measurement RMS, use conservativeCoefficient of the old one and (1-conservativeCoefficient) of the new one
// public static [] oldErrors;
            if (oldErrors!=null) ERR_others=conservativeCoefficient*oldErrors[t]+(1.0-conservativeCoefficient)*ERR_others; // make running average with previous errors squares
//            ERR_this= (SX_all[t]-amplitudes[n][t])/N - amplitudes[n][t];
            ERR_this= SX_others[t]/(N-1) - amplitudes[n][t];
            ERR_this= Math.abs(ERR_this); // distance from the centroid of other measurments
//    System.out.println("--- amplitudes["+n+"]="+amplitudes[n]+" ERR_this="+ERR_this+", ERR_others=" + ERR_others);
//    System.out.println("SX_others="+SX_others+" SX2_others="+SX2_others+" SX_all="+SX_all+" SX2_all="+SX2_all);


            if (ERR_this > badMeasurementError*ERR_others) { // too bad, replace the measuremnet
              allGood=false;
              if (oldErrors!=null) System.out.println("Bad measurement, re-measuring. ERR_this="+ERR_this+", ERR_others=" + ERR_others+" (oldErrors["+t+"]="+oldErrors[t]);
              else                 System.out.println("Bad measurement, re-measuring. ERR_this="+ERR_this+", ERR_others=" + ERR_others);
              for (i=0;i<N;i++) System.out.println(((i==n)?"-*-":"---")+"+ old amplitudes["+i+"]["+t+"]="+amplitudes[i][t]);
              System.out.println("SX_others["+t+"]="+SX_others[t]+" SX2_others["+t+"]="+SX2_others[t]+" SX_all["+t+"]="+SX_all[t]+" SX2_all["+t+"]="+SX2_all[t]);
              imp_camera=jp4_instance.openURL(imp_camera);
              if (imp_camera==null) {
                IJ.showMessage("Camera Image ERROR","Failed to open camera image");
                return null;
              }
              measOne=   measureTargets (targetLocs, imp_camera);




              for (t=0;t<T;t++) for (i=0;i<rslt[0].length;i++) for (j=0;j<rslt[0][0].length;j++) measAll[n][t][i][j]=measOne[t][i][j];
              for (t=0;t<T;t++) {
                amplitudes[n][t]=0.0;
                for (i=0;i<rslt[0].length;i++) for (j=0;j<rslt[0][0].length;j++) amplitudes[n][t]+=measAll[n][t][i][j]*measAll[n][t][i][j];
                SX_all[t]= SX_others[t]+ amplitudes[n][t];
                SX2_all[t]=SX2_others[t]+amplitudes[n][t]*amplitudes[n][t];
                System.out.println(" new amplitudes["+n+"]["+t+"]="+amplitudes[n][t]+" SX_all["+t+"]="+SX_all[t]+" SX2_all["+t+"]="+SX2_all[t]);
//              System.out.println("+++n="+n+" amplitudes["+n+"]="+amplitudes[n]+" SX_all="+SX_all+", SX2_all=" + SX2_all);
              }
              break;
            }
          } // next t
          if (!allGood) break;
        } // next n
      } // good or bad measurements - here we are, let's use whatevcr available
    }
// Update errors
    for (t=0;t<T;t++) rms[t]=Math.sqrt(SX2_all[t]/N  - (SX_all[t]/N)*(SX_all[t]/N));
    if (oldErrors==null) {
      oldErrors=new double[T];
      for (t=0;t<T;t++) oldErrors[t]=rms[t];
    }
    for (t=0;t<T;t++) oldErrors[t]= conservativeCoefficient* oldErrors[t] + (1.0-conservativeCoefficient)*rms[t];

    for (t=0;t<rslt.length;t++) for (i=0;i<rslt[0].length;i++) for (j=0;j<rslt[0][0].length;j++) rslt [t][i][j]=0.0;
    for (n=0;n<N;n++)  for (t=0;t<rslt.length;t++) for (i=0;i<rslt[0].length;i++) for (j=0;j<rslt[0][0].length;j++) rslt[t][i][j]+=measAll[n][t][i][j];
    for (t=0;t<rslt.length;t++) for (i=0;i<rslt[0].length;i++) for (j=0;j<rslt[0][0].length;j++) rslt [t][i][j]/=N;

//  System.out.println("measureTargets("+N+")-->");
    return rslt;
  }

//badMeasurementError

  public double [][][] measureTargets (int [][][] targetLocs, ImagePlus imp) {
    int nTarg, nBayer, i,j,n, indx;
    int  FFTSize2=FFTSize<<1;
    int calib_width=calib_sample&(~3); // just in case
    int calib_height=calib_sample&(~3); // just in case
    double [][][] rslt=new double [targetLocs.length][2][4];
    float [] max2=new float[2];
    double dwb; //w-b calibration
    float [] pixels_black;    // pixels for black calibration (all Bayer components) */
    float [] pixels_white;    // pixels for white calibration (all Bayer components) */
    float [] pixels_vert_all; // pixels for vertical   pattern (all Bayer components), 2*FFTSize x  2*FFTSize */
    float [] pixels_hor_all;  // pixels for horizontal pattern (all Bayer components), 2*FFTSize x  2*FFTSize */
    float [] pixels_vert=new float [FFTSize*FFTSize];     // pixels for vertical   pattern (selected Bayer component), FFTSize x  FFTSize */
    float [] pixels_hor= new float [FFTSize*FFTSize];     // pixels for horizontal pattern (selected Bayer component), FFTSize x  FFTSize */
    ImageProcessor ip=imp.getProcessor();
//    System.out.println("measureTargets(single)");
    
    
    for (nTarg=0;nTarg<rslt.length;nTarg++) {
      ip.setRoi(targetLocs[nTarg][3][0]&(~1)-(calib_width>>1),targetLocs[nTarg][3][1]&(~1)-(calib_height>>1),calib_width,calib_height);
      pixels_black=    (float[])ip.crop().getPixels();
      ip.setRoi(targetLocs[nTarg][2][0]&(~1)-(calib_width>>1),targetLocs[nTarg][2][1]&(~1)-(calib_height>>1),calib_width,calib_height);
      pixels_white=    (float[])ip.crop().getPixels();
      ip.setRoi(targetLocs[nTarg][0][0]&(~1)-FFTSize,targetLocs[nTarg][0][1]&(~1)-FFTSize,FFTSize2,FFTSize2);
      pixels_vert_all=(float[])ip.crop().getPixels();
      ip.setRoi(targetLocs[nTarg][1][0]&(~1)-FFTSize,targetLocs[nTarg][1][1]&(~1)-FFTSize,FFTSize2,FFTSize2);
      pixels_hor_all=(float[])ip.crop().getPixels();
      for (nBayer=0;nBayer<4;nBayer++) {
/* Calibrate White-Black for selected target/color */
        dwb=0.0;
        n=0;
        for (i=(nBayer>>1); i<calib_height;i+=2) for (j=(nBayer & 1); j<calib_width; j+=2) {
          indx=i*calib_width+j;
if (indx>=pixels_white.length) {
    System.out.println("*** BUG in measureTargets(), pixels_white.length="+pixels_white.length+" indx="+indx+" i="+i+" j="+j+" calib_width="+calib_width);
    

        } else {
           dwb+=pixels_white[indx]-pixels_black[indx];
/*
measureTargets(single)
Exception in thread "AWT-EventQueue-0" java.lang.ArrayIndexOutOfBoundsException: 160
        at focus_tuning.measureTargets(focus_tuning.java:1433)
        at focus_tuning.measureTargets(focus_tuning.java:1348)
        at focus_tuning.scanMotors(focus_tuning.java:3602)
        at focus_tuning.actionPerformed(focus_tuning.java:354)

*/
            n++;
          }
        } 
        dwb/=n;
/* Extract bayer component for vertical and horizontal patterns */
        n=0;
        for (i=(nBayer>>1); i<FFTSize2;i+=2) for (j=(nBayer & 1); j<FFTSize2; j+=2) {
          indx=i*FFTSize2+j;
          pixels_vert[n]=  pixels_vert_all[indx]; // out of bound exception 3584
          pixels_hor [n++]=pixels_hor_all[indx];
        }      
/* Calculate vertical pattern MTF at patter frequency (horizontal resolution). Uses global Hamming and TargetFilter arrays */
/* TODO: ? optimize measuresOne() to remove extra processing - copying, both directions */
        max2=measuresOne(pixels_vert, FFTSize, Hamming, TargetFilter, 0, 0, FFTSize, false, "");
        rslt[nTarg][0][nBayer]=max2[0]/dwb; // Is it right 0/1 here (and below)?
        max2=measuresOne(pixels_hor, FFTSize, Hamming, TargetFilter, 0, 0, FFTSize, false, "");
        rslt[nTarg][1][nBayer]=max2[1]/dwb;
      }
    }
    return rslt;
  }


/// =========================================

  public void showTargets(double[][][] targets, double threshold, String title) {
    String header="#\tVert\tHor\tY\tX\tValue";
    StringBuffer sb = new StringBuffer();
    int i,j,n;
    n=0;
    for (i=0;i<targets.length;i++)  for (j=0;j<targets[0].length;j++) if (targets[i][j][2]>threshold) {
      n++;
      sb.append(n+"\t"+i+"\t"+j+"\t"+IJ.d2s(targets[i][j][0],1)+"\t"+IJ.d2s(targets[i][j][1],1)+"\t"+IJ.d2s(targets[i][j][2],4)+"\n");
    }
    new TextWindow(title+"_target_centers", header, sb.toString(), 400,300);
  }


/// =========================================
//public static float [][] input_bayer; 
 /* ignore ROI, use whole image */
  public float[][] splitBayer (ImagePlus imp) {
    ImageProcessor ip=imp.getProcessor();
    Rectangle r=new Rectangle(imp.getWidth(),imp.getHeight());
    float [] pixels;
    pixels=(float[])ip.getPixels();    
    if (DEBUG_LEVEL>10) IJ.showMessage("splitBayer","r.width="+r.width+
                              "\nr.height="+r.height+
                              "\nlength="+pixels.length);
    float [][] bayer_pixels=new float[4][pixels.length>>2];
    int x,y,base,base_b,bv;
    int half_height=r.height>>1;
    int half_width=r.width>>1;
    for (y=0; y<half_height; y++) for (bv=0;bv<2;bv++){
      base=r.width*((y<<1)+bv);
      base_b=half_width*y;
      if (bv==0) for (x=0; x<half_width; x++) {
        bayer_pixels[0][base_b]=  pixels[base++];
        bayer_pixels[1][base_b++]=pixels[base++];
      } else  for (x=0; x<half_width; x++) {
        bayer_pixels[2][base_b]=  pixels[base++];
        bayer_pixels[3][base_b++]=pixels[base++];
      }
    }
    return bayer_pixels;
  }

  public float[] initHamming(int size) {
    float [] hamming =new float [size*size];
    float [] hamming_line=new float [size];
    int i,j;
    for (i=0; i<size; i++) hamming_line[i]=(float) (0.54-0.46*Math.cos((i*2.0*Math.PI)/size));
    for (i=0; i<size; i++) for (j=0; j<size; j++){
       hamming[size*i+j]=hamming_line[i]*hamming_line[j];
    }
    return hamming;
  }


  public float[] initTargetFilter(int size, double tolerance, double frequency, double angle ) {
    float [] filter =new float [size*size];
    double [][] d2 =new double [size][size];
    double  d2a;
    double x,y,di,dd;
    int i,j, ir,jr;
    double k=tolerance*2*size*frequency;
    double max;
    k=1.0/k/k;
    x=2*size*frequency*Math.cos(Math.PI*angle/180);
    y=2*size*frequency*Math.sin(Math.PI*angle/180);
    for (i=0;i<size;i++) {
      di=(i>(size>>1))? (i-size):i;
      for (j=0;j<size;j++){ 
        //  d[i,j]=(j-x)*(j-x)+(di-y)*(di-y) ;
        dd=(j-x)*(j-x)+(di-y)*(di-y) ;
        d2[i][j]=Math.exp(-k*dd);
      }
    }
    /* now add the four of them */
    max=0;
    for (i=0; i<size;i++) for (j=0; j<size; j++) {
      d2a=d2[i][j];
      ir=j;              jr=(size-i)%size;  d2a+=d2[ir][jr];
      ir=(size-i)%size;  jr=(size-j)%size;  d2a+=d2[ir][jr];
      ir=(size-j)%size;  jr=i;              d2a+=d2[ir][jr];
      if (max<d2a) max= d2a;
      filter [i*size+j]=(float) d2a;
    }
    max=1.0/max;
    for(i=0;i<filter.length;i++)  filter[i]*= (float) max;

    return filter;
  }


//        l = ((size-row)%size) * size + (size-c)%size;


// public static double targetFTolerance=0.2; /// tolerance in frequency domain
// public static double targetFX = 0.23;         /// target frequency (cpp, )
// public static double targetAngle=5.0;  /// target angle

// private static float[][] vh_pixels;
//  public int [] findMaxPix(float[] pixels, int width) { // returns x,y
//  public float [] multiplyTargets(float [][] vh_pixels, int width, int height, int dx, int dy) {
  public float [][] useBestComponent(float [][][] cvh_pixels) {
     float[][] best=new float[2][cvh_pixels[0][0].length];
     int i,n,k;
     for (i=0;i<cvh_pixels[0][0].length;i++) for (k=0;k<2;k++)  {
       best[k][i]=cvh_pixels[0][k][i];
       for (n=1;n<4;n++) if (cvh_pixels[n][k][i]>best[k][i]) best[k][i]=cvh_pixels[n][k][i];
     }
    return best;
  }


  public int [] findMaxPix(float[] pixels, int width) { // returns x,y
    int i, imax;
    float max=pixels[0];
    imax=0;
    for (i=1;i<pixels.length;i++) if (pixels[i]>max) {
        max=pixels[i];
        imax=i;
    }
    int [] xy=new int[2];
    xy[1]= imax / width;
    xy[0]= imax % width;
    return xy;
  }

  public void equalize(float [] pixels, int width, int height, int radius) {
    float []max=new float [pixels.length];
    float maxmax=pixels[0];
    float m;
    int i,j,i1,j1,il,ih, jl,jh;
    float threshold=0.0001f;
    for (i=0; i<height; i++) {
      il=i-radius;
      if (il<0) il=0;
      ih=i+radius+1;
      if (ih>height) ih=height;
      for (j=0; j<width; j++) {
        jl=j-radius;
        if (jl<0) jl=0;
        jh=j+radius+1;
        if (jh>width) jh=width;
        m=pixels[il*width+jl];
        for (i1=il;i1<ih;i1++) for (j1=jl;j1<jh;j1++) {
          if (pixels[i1*width+j1]> m) m = pixels[i1*width+j1];
        }
        if (m > maxmax) maxmax=m;        
        max[i*width+j]=m;
      }
    }
    maxmax=1.0f/maxmax;
    for (i=0;i<pixels.length; i++ ) {
      if (max[i]<threshold) pixels[i]=0.0f;
      else pixels[i]/=Math.sqrt(max[i]/maxmax);
    }
  }  

//      minFilter(vh_pixels[0], mwidth, mheight, minFilterRadius);
/* Replaces pixel value with minimal in the square*/
  public float[] minFilter(float [] pixels, int width, int height, int radius) {
    float []min=new float [pixels.length];
    float m;
    int i,j,i1,j1,il,ih, jl,jh;
    for (i=0; i<height; i++) {
      il=i-radius;
      if (il<0) il=0;
      ih=i+radius+1;
      if (ih>height) ih=height;
      for (j=0; j<width; j++) {
        jl=j-radius;
        if (jl<0) jl=0;
        jh=j+radius+1;
        if (jh>width) jh=width;
        m=pixels[il*width+jl];
        for (i1=il;i1<ih;i1++) for (j1=jl;j1<jh;j1++) {
          if (pixels[i1*width+j1] < m) m = pixels[i1*width+j1];
        }
        min[i*width+j]=m;
      }
      
    }
    return min;
  }  




  public float [][] scanTargets(float[] full_pixels, int size, int step, float [] hamming, float [] target_filter, int width, int height, boolean show, String title) {
    int mwidth=(width-size)/step;
    int mheight=(height-size)/step;
    int i,j;
    float [] max2;
    float[][] vh_pixels=new float[2][mwidth*mheight];
    for (i=0; i<mheight; i++) for (j=0; j<mwidth; j++) {
      max2= measuresOne(full_pixels, size, hamming, target_filter, (j*step)+(size>>1), (i*step)+(size>>1), width, false, title);
      vh_pixels[0][i*mwidth+j]=max2[0];
      vh_pixels[1][i*mwidth+j]=max2[1];
    }
    return vh_pixels;
  }
 
  public float [] multiplyTargets(float [][] vh_pixels, int width, int height, int dx, int dy) {
    int i,j,i1,j1;
    float [] product=new float [vh_pixels[0].length];
    for (i=0;i<height;i++) {
      i1=(i+dy+height) %height;
      for (j=0;j<width;j++) {
        j1=(j+dx+width) %width; // does not work with negative
        product[i*width+j]=vh_pixels[0][i*width+j]*vh_pixels[1][i1*width+j1];
      }
    }
    return product;
  }


  public float [] correlateTargets(float [][] vh_pixels, int width, int height, int radius) {
    int w=2*radius+1;
    float [] corr =new float [ w*w];
    double d;
    int i,j, i1, j1, x,y;
    for (y=-radius; y<=radius; y++) for (x=-radius; x<=radius; x++) {
      d=0;
      for (i=0;i<height;i++) {
        i1=(i+y+height) %height;
        for (j=0;j<width;j++) {
          j1=(j+x+width) %width; // does not work with negative
          d+=vh_pixels[0][i*width+j]*vh_pixels[1][i1*width+j1];
        }
      }
      corr[(y+radius)*w+(x+radius)]= (float) d;
    }
    return corr;
  }

/* if (x0==0) full_pixels are considered to be just the square to process, no windowing is performed */
  public float[] measuresOne(float[] full_pixels, int size, float [] hamming, float [] target_filter, int x0, int y0, int width, boolean show, String title) {
    int i,j,m,x,y, base0, base;
    int half_size=size>>1;
    int quater_size=size>>2;
    int threequaters_size=quater_size+half_size;
    int x00=x0-half_size;
    int y00=y0-half_size;
    int row, l, c;
    float p;
    float [] pixels;
    double average=0.0;
    if (x0>0) {
      pixels = new float [size*size];
      base=0;
      for (y=0;y<size;y++) {
        base0=(y00+y)*width+x00;
        for (x=0;x<size; x++ ) {
          average+=full_pixels[base0];
          pixels[base++]=full_pixels[base0++];
        }
      }
    } else {
      for (base0=0; base0<full_pixels.length;base0++) average+=full_pixels[base0];
      pixels=full_pixels;
    }
    average/=pixels.length;
    float faverage=(float) average;

    for (i=0;i<pixels.length; i++ ) {
      pixels[i]=(pixels[i]-faverage)*hamming[i];
    }

    ImageProcessor ip= new FloatProcessor(size,size);
    ip.setPixels(pixels);
    FHT fht=new FHT(ip);
    fht.swapQuadrants();
    if (show && (DEBUG_LEVEL>6)) {
//      fht.resetMinAndMax();
//      ImagePlus imp=  new ImagePlus(title+"_"+x0+"-"+y0, ip);
      ImageProcessor ip1= new FloatProcessor(size,size);
      ip1.setPixels(pixels);
      ImagePlus imp1=  new ImagePlus(title+"_"+x0+"-"+y0, ip1);
      imp1.getProcessor().resetMinAndMax();
      imp1.show();
    }
    fht.transform();
    pixels=(float[])fht.getPixels();
    for (row=0; row<=half_size; row++) {
      base=row*size;
      for (c=0; c<((row==half_size)?(half_size+1):size); c++) {
        i=base+c;
        l = ((size-row)%size) * size + (size-c)%size;
        p=(float) (Math.sqrt(pixels[i]*pixels[i] + pixels[l]*pixels[l]));
        pixels[i]=p;
        pixels[l]=p;
      }
    }
    float [] max2=new float[2];
    max2[0]=0.0f;
    max2[1]=0.0f;

    for (i=0;i<quater_size;i++) for (j=0;j<size;j++) {
      m=((j<quater_size) || (j>=threequaters_size))?1:0;
      base=(((m>0)?quater_size:0)+i)*size+j;
      p=pixels[base]*target_filter[base];
      if (p>max2[m]) max2[m]=p;
    }



    if (show && (DEBUG_LEVEL>3)) {
      ImageStack stack =   fht.getComplexTransform();
      ImagePlus imp2 = new ImagePlus(title+"_fht_stack", stack);
      imp2.getProcessor().resetMinAndMax();
      imp2.show();
    }
    if (show && (DEBUG_LEVEL>4)) {
      ImageProcessor ip3=new FloatProcessor(size,size);
      swapQuads(pixels, size);
      ip3.setPixels(pixels);
      ImagePlus imp3 = new ImagePlus(title+"_abs", ip3);
      imp3.getProcessor().resetMinAndMax();
      imp3.show();
    }
    return max2;
  }

  void swapQuads(float[] pixels, int size) {
    int i, j;
    float p;
    int m=(size>>1) | ((size*size)>>1);
    for (i=0;i<(pixels.length>>1);i++) {
      j= i ^ m;
      p=pixels[i];
      pixels[i]=pixels[j];
      pixels[j]=p;
    }
  }

  public void showBayers(float[][] bayer_pixels, int width, int height, String title) {
    int i;
   if (DEBUG_LEVEL>10) IJ.showMessage("showBayers","width="+width+
                              "\nheight="+height+
                              "\nlength="+bayer_pixels[0].length);

    ImageProcessor[] ip= new ImageProcessor[4];
    ImagePlus[]      imp=new ImagePlus[4];
    for (i=0;i<4;i++) {
      ip[i]=new FloatProcessor(width,height);
      ip[i].setPixels(bayer_pixels[i]);
      ip[i].resetMinAndMax();
      imp[i]=  new ImagePlus(title+"_"+i, ip[i]);
      imp[i].show();
    }

  }

/*
  public void generatePlots(double[][] Vectors, String plotType, String plotName, Color [] colors, boolean asText, String headers, int textWidth, int TextHeight) {
    MTF_Bayer_instance.generatePlots(Vectors, plotType, plotName, colors, asText, headers, textWidth, TextHeight);
  }
 Color [] bayer_colors=  { Color.green, Color.red, Color.blue, Color.cyan, Color.black };
 public static  Color [] resolution_colors=  { Color.Color(1.0,0.3,0.3),
                                               Color.Color(0.7,0.0,0.0),
                                               Color.Color(0.3,1.0,0.3),
                                               Color.Color(0.0,0.7,0.0),
                                               Color.Color(0.3,0.3,1.0),
                                               Color.Color(0.0,0.0,0.7)};
Color(float r, float g, float b)
          Creates an opaque sRGB color with the specified red, green, and blue values in the range (0.0 - 1.0).

    plotResult = new Plot(allTitle, ejeX, ejeY, (double []) null, (double []) null);
    plotResult.setLimits(0,0.5,thisYmin,thisYmax);

    for (n=0; n<finalVectors.length; n++) {
      plotResult.setColor(colors[n]);
      plotResult.addPoints(xValues, finalVectors[n], plotResult.LINE);
    }
    plotResult.draw();
    plotResult.show();

    String ejeX="pixel";
    String ejeY="";
targetF
*/

  public void plotScanTable(double[][][][][] scanData, double micronsPerSample, double sampleFrequency,  String title, Color[] colors) {
    int nSamples=scanData.length;
    int nDirections=scanData[0].length;
    int nTargets=scanData[0][0].length;
    double []   xValues=new double[nSamples];
    double [][][][] yValues=new double [nDirections][3][2][nSamples];
    int sample, target, color, orientation, direction;

    double yMax=0.0;
    String ejeX="distance (um)";
    String ejeY="contrast @"+IJ.d2s(sampleFrequency,2)+"cycles/pix (%)";

    for (sample=0;sample < nSamples; sample++) xValues[sample]= micronsPerSample * (sample - 0.5*scanData.length);
    for (target=0; target < nTargets; target++ ) {
      yMax=0;
      for (sample=0;sample < nSamples; sample++) {
        for (direction=0;direction<nDirections;direction++) {
          for (orientation=0;orientation<2;orientation++) {
            yValues[direction][0][orientation][sample]=     scanData[sample][direction][target][orientation][1]; //red
            yValues[direction][1][orientation][sample]=0.5*(scanData[sample][direction][target][orientation][0]+
                                                            scanData[sample][direction][target][orientation][3]); // green
            yValues[direction][2][orientation][sample]=     scanData[sample][direction][target][orientation][2]; //blue
            for (color=0; color<3; color++) if (yValues[direction][color][orientation][sample] > yMax) yMax=yValues[direction][color][orientation][sample];
          }
        }
      }
      plotResult = new Plot(title+" - target "+(target+1), ejeX, ejeY, (double []) null, (double []) null);
      plotResult.setLimits(xValues[0], xValues[nSamples-1],0,yMax);
      for (direction=0;direction<nDirections;direction++) for (color=0; color<3; color++) for (orientation=0;orientation<2;orientation++) {
        plotResult.setColor(colors[(color<<1)+orientation]);
        plotResult.addPoints(xValues,yValues[direction][color][orientation], Plot.LINE);
      }
      plotResult.draw();
      plotResult.show();
    }
  }

  public void showScanTable(double[][][][][] scanData, int [] startScan, int [] endScan, int motors_scan_steps, double[][][][] maxScanTargetOrientColor, String title, int precision) {
    String header="#\tM1\tM2\tM3"; //FV_R\tFH_R\tFV_G\tFH_G\tFV_B\tFH_B\tRV_R\tRH_R\tRV_G\tRH_G\tRV_B\tRH_B";
    StringBuffer sb = new StringBuffer();
    int [] targetPositions = new int [3];
    int n,i,t;
    if (scanData==null) {
      IJ.showMessage("showScanTable Error","data is null, nothing to display");
      return;
    }
    for (t=1; t<=scanData[0][0].length;t++) {
      header+="\t"+t+"FV_R\t"+t+"FH_R\t"+t+"FV_G\t"+t+"FH_G\t"+t+"FV_B\t"+t+"FH_B";
    }
    if (scanData[0].length>1) for (t=1; t<=scanData[0][0].length;t++) {
      header+="\t"+t+"RV_R\t"+t+"RH_R\t"+t+"RV_G\t"+t+"RH_G\t"+t+"RV_B\t"+t+"RH_B";
    }
//    sb.append(header); // repeat header for easy copying
    for (n=0;n<=motors_scan_steps;n++) {
      for (i=0; i<3;i++)  {
        targetPositions[i]=startScan[i];
        if (startScan[i]!=endScan[i]) targetPositions[i]+=((endScan[i]-startScan[i])*n)/motors_scan_steps;
      }
      sb.append((n+1)+
                  "\t"+targetPositions[0]+
                  "\t"+targetPositions[1]+
                  "\t"+targetPositions[2]);
      for (t=0; t<scanData[0][0].length;t++) { /// forward scan
        sb.append("\t"+IJ.d2s(     scanData[n][0][t][0][1],precision)+                          // vert, R , forward
                  "\t"+IJ.d2s(     scanData[n][0][t][1][1],precision)+                          // hor,  R , forward
                  "\t"+IJ.d2s(0.5*(scanData[n][0][t][0][0]+scanData[n][0][t][0][3]),precision)+ // vert, (Gr+Gb)/2 , forward
                  "\t"+IJ.d2s(0.5*(scanData[n][0][t][1][0]+scanData[n][0][t][1][3]),precision)+ // hor,  (Gr+Gb)/2 , forward
                  "\t"+IJ.d2s(     scanData[n][0][t][0][2],precision)+                          // vert, B , forward
                  "\t"+IJ.d2s(     scanData[n][0][t][1][2],precision));                          // hor,  B , forward
      } // next target
      if (scanData[0].length>1)  for (t=0; t<scanData[0][0].length;t++) { /// reverse scan
        sb.append("\t"+IJ.d2s(     scanData[n][1][t][0][1],precision)+                          // vert, R , reverse
                  "\t"+IJ.d2s(     scanData[n][1][t][1][1],precision)+                          // hor,  R , reverse
                  "\t"+IJ.d2s(0.5*(scanData[n][1][t][0][0]+scanData[n][0][t][0][3]),precision)+ // vert, (Gr+Gb)/2 , reverse
                  "\t"+IJ.d2s(0.5*(scanData[n][1][t][1][0]+scanData[n][0][t][1][3]),precision)+ // hor,  (Gr+Gb)/2 , reverse
                  "\t"+IJ.d2s(     scanData[n][1][t][0][2],precision)+                          // vert, B , reverse
                  "\t"+IJ.d2s(     scanData[n][1][t][1][2],precision));                          // hor,  B , reverse
      } // next target


      sb.append("\n");
    } // next position
    if (maxScanTargetOrientColor!=null) {
        sb.append("Max at\t---\t---\t---");

      for (n=0;n<maxScanTargetOrientColor.length;n++) for (t=0; t<maxScanTargetOrientColor[0].length;t++) {
        sb.append("\t"+IJ.d2s(maxScanTargetOrientColor[n][t][0][0],precision)+  // vert, R , forward
                  "\t"+IJ.d2s(maxScanTargetOrientColor[n][t][1][0],precision)+  // hor,  R , forward
                  "\t"+IJ.d2s(maxScanTargetOrientColor[n][t][0][1],precision)+  // vert, (Gr+Gb)/2 , forward
                  "\t"+IJ.d2s(maxScanTargetOrientColor[n][t][1][1],precision)+  // hor,  (Gr+Gb)/2 , forward
                  "\t"+IJ.d2s(maxScanTargetOrientColor[n][t][0][2],precision)+  // vert, B , forward
                  "\t"+IJ.d2s(maxScanTargetOrientColor[n][t][1][2],precision)); // hor,  B , forward
      }
      sb.append("\n");
    }
    new TextWindow(title, header, sb.toString(), 1200,800);
  }
/* Calculates astigmatism as difference between "best focus" distance for tangential and radial off-center components
     and vertical-horizontal for the center*/
  public void showAstigmatismTable(double[][][][] maxScanTargetOrientColor, //[scan direction][target][v,h][r,g,b] 1.0 - one scan sample
                                   double micronsPerSample,  // microns per 1 scan sample
                                   int range, // scan range (samples, starting from 0)
                                   double reliabilityThreshold, // maximums outside of [reliabilityThreshold..range-reliabilityThreshold are
                                                                // considered unreliable
                                   String title,                // window title
                                   int precision) {             // output precision in decimals
 /* targets 1,3,5:  v=t, h=r; targets 2,4: v=r, h=t  */
    double  [][] astigm =new double [maxScanTargetOrientColor[0].length][3];
    boolean [][] reliable =new boolean [maxScanTargetOrientColor[0].length][3];
    int target, color;
    double  [] maxAstigmColor= new double[3];
    boolean [] reliableColor= new boolean[3];
    double  maxAstigmRedGreen;
    boolean reliableRedGreen;
    for (color=0;color<3;color++) {
      maxAstigmColor[color]=0.0;
      reliableColor[color]=true;
    }
    for (target=0; target <astigm.length; target++) for (color=0;color<astigm[0].length;color++) {
      astigm[target][color]=micronsPerSample*(maxScanTargetOrientColor[0][target][0][color]-maxScanTargetOrientColor[0][target][1][color]);
      if (((target+1)==2) || ((target+1)==4)) astigm[target][color]=-astigm[target][color];
      reliable[target][color] = (maxScanTargetOrientColor[0][target][0][color] >= reliabilityThreshold ) &&
                                (maxScanTargetOrientColor[0][target][0][color] <= range-reliabilityThreshold ) &&
                                (maxScanTargetOrientColor[0][target][1][color] >= reliabilityThreshold ) &&
                                (maxScanTargetOrientColor[0][target][1][color] <= range-reliabilityThreshold );
     if (Math.abs(astigm[target][color])>maxAstigmColor[color]) maxAstigmColor[color]=Math.abs(astigm[target][color]);
     reliableColor[color]=reliableColor[color] && reliable[target][color];
    }
    maxAstigmRedGreen=maxAstigmColor[0];
    if (maxAstigmColor[1]> maxAstigmRedGreen) maxAstigmRedGreen=maxAstigmColor[1];
    reliableRedGreen=reliableColor[0] && reliableColor[1];

    String header="Target\tColor\tAstigmatism(um)\tNotes";
    StringBuffer sb = new StringBuffer();
    String [] colors={"Red","Green","Blue"};
    for (target=0; target <astigm.length; target++) {
      sb.append("Target"+(target+1)+"\t\t\t\n");
      for (color=0;color<astigm[0].length;color++) {
        sb.append("\t"+colors[color]+"\t"+((astigm[target][color]>0)?"+":"")+IJ.d2s(astigm[target][color],precision)+"\t");
        if (!reliable[target][color])sb.append("Maximum out of reliable range");
        sb.append("\n");
      }
    }
    sb.append("Maximal\t\t\t\n");
    for (color=0;color<astigm[0].length;color++) {
      sb.append("\t"+colors[color]+"\t"+IJ.d2s(maxAstigmColor[color],precision)+"\t");
      if (!reliableColor[color])sb.append("Some out of reliable range");
      sb.append("\n");
    }
    sb.append("Result\tRed+Green\t"+IJ.d2s(maxAstigmRedGreen,precision)+"\t");
    if (!reliableRedGreen)sb.append("Some out of reliable range");
    sb.append("\n");
    new TextWindow(title, header, sb.toString(), 800,600);
  }



  public double [][][] filterScanTable(double[][][][][] scanData,
                                 int averageNum) { // use running average to smouth data
    int n,i,j,k,t;
    double [][][] targetScanData=new double[scanData.length][scanData[0][0].length][7]; // r,g,b 
    double [] smouthArr=new double[scanData.length];
    StringBuffer sb = new StringBuffer();
    String header="#";
    for (t=0; t<scanData[0][0].length;t++) header+="\tRv"+(t+1)+"\tRh"+(t+1)+"\tGv"+(t+1)+"\tGh"+(t+1)+"\tBv"+(t+1)+"\tBh"+(t+1)+"\tSum"+(t+1);
    for (n=0;n<scanData.length;n++) for (t=0; t<scanData[0][0].length;t++) { 
        targetScanData[n][t][0]= scanData[n][0][t][0][1]; // R, vertical
        targetScanData[n][t][2]=0.5*(scanData[n][0][t][0][0]+scanData[n][0][t][0][3]); // G, vertical
        targetScanData[n][t][4]= scanData[n][0][t][0][2]; // B, vertical

        targetScanData[n][t][1]= scanData[n][0][t][1][1]; // R, horizontal
        targetScanData[n][t][3]=0.5*(scanData[n][0][t][1][0]+scanData[n][0][t][1][3]); // G, horizontal
        targetScanData[n][t][5]= scanData[n][0][t][1][2]; // B, horizontal
        targetScanData[n][t][6]=0.0;
        for (i=0;i<6;i++) targetScanData[n][t][6]+=targetScanData[n][t][i]*targetScanData[n][t][i];
        targetScanData[n][t][6]=Math.sqrt(targetScanData[n][t][6]);
// normalize componets data
        for (i=0;i<6;i++) targetScanData[n][t][i]/=targetScanData[n][t][6];
    } 
    for (t=0; t<scanData[0][0].length;t++) for (i=0; i<7; i++) {
       for (n=0;n<scanData.length;n++)  {
         smouthArr[n]=0.0;
         for (j=0;j<averageNum;j++) {
           k=n+j;
           if (k<0) k=0;
           else if (k>=scanData.length) k=scanData.length-1;
           smouthArr[n]+=targetScanData[k][t][i];
         }
       }
       for (n=0;n<scanData.length;n++)  targetScanData[n][t][i]=smouthArr[n]/averageNum; // copy back in-place
    }
    if (DEBUG_LEVEL>0 ) {
       for (n=0;n<scanData.length;n++)  {
         sb.append(n);
         for (t=0; t<scanData[0][0].length;t++) for (i=0; i<7; i++) {
           sb.append("\t"+targetScanData[n][t][i]);
         }
         sb.append("\n");
       }
       new TextWindow(((imp_src==null)?"":imp_src.getTitle())+"_signature_table", header, sb.toString(), 1200,800);
       System.out.println("Created Filtered Scan Table");

    }

    return targetScanData;
  }

//  wnd - miniumal half-size of the maximum search area, maxSlopes may only extend it
// finds position of the maximum in a (noisy) arrays. All "real" values are assumed to be >0,  0.0 in array is a special case (masked out), it should not be used in calculations
  double findMaxInWndArray(double[] array, int wnd, double maxSlopes) {
   double max=array[0];
   double min;
   int i;
   double threshold;
   int imx=0;
   int ilow,ihigh;
   double a,b,c;
   double x;
   double SX4=0.0;
   double SX3=0.0;
   double SX2=0.0;
   double SX=0.0;
   double SYX2=0.0;
   double SYX=0.0;
   double SY=0.0;
   double S0=0.0;
   double crslt,rslt;
   StringBuffer sb = new StringBuffer();
   String header="#\tarray";
   for (i=1;i<array.length;i++) {
     if (array[i]>max){
       max=array[i];
       imx=i;
     }
     if (DEBUG_LEVEL>1 )  sb.append(i+"\t"+array[i]+"\n");
   }
   min=max;
   for (i=1;i<array.length;i++)  if ((array[i]!=0.0) && (array[i]<min))  min=array[i]; // 0.0 is a special case - masked out non-existent elements
   threshold=min+ (max-min)*maxSlopes;
   if ((imx==0 ) || (imx==(array.length-1))) return 1.0*imx;

   ilow= imx-wnd;
   ihigh=imx+wnd;
   if (ilow<0) ilow=0;
   if (ihigh>array.length-1) ihigh=array.length-1;

// verify that they are not zero (exactly)
   for (i=imx; (i>ilow) && (array[i-1]>0.0); i--);
   if (i>ilow) ilow=i;
   for (i=imx; (i<ihigh) && (array[i+1]>0.0); i++);
   if (i<ihigh) ihigh=i;

//   while ((wnd>0) && ((array[imx-wnd]==0.0) || (array[imx+wnd] == 0.0))) wnd--; 



   if (imx<wnd) wnd=imx;
   if (imx>(array.length-1-wnd)) wnd=array.length-1-imx;
// verify that they are not zero (exactly)
//   while ((wnd>0) && ((array[imx-wnd]==0.0) || (array[imx+wnd] == 0.0))) wnd--; 


   while ((ilow>1) && (array[ilow-1]>threshold)) ilow--; // extend range downwards until lower than threshold  (or end of array)
   while ((ihigh<(array.length-2)) && (array[ihigh+1]>threshold)) ihigh++;// extend range upwards until lower than threshold (or end of array)
// Now range may be not ceneters ad imx
   for (i=ilow;i <= ihigh ; i++) {
     x= 0.5*(2*i-ilow-ihigh);
     SX4+=x*x*x*x;
     SYX2+=array[i]*x*x;
     SX3+=x*x*x;
     SX2+=x*x;
     SYX+=array[i]*x;
     SY+=array[i];
     SX+=x;
     S0+=1.0;

   }
/**
(1)         a*SX4 +b*SX3 + c*SX2 -SYX2 =0
(2)         a*SX3 +b*SX2 + c*SX  -SYX  =0
(3)         a*SX2 +b*SX  + c*S0  -SY   =0
(1a)        a*SX4*SX  +b*SX3*SX  + c*SX2*SX  -SYX2*SX =0
(2a)        a*SX3*SX2 +b*SX2*SX2 + c*SX*SX2  -SYX*SX2  =0
(1a-2a)     a*(SX4*SX-SX3*SX2) + b* (SX3*SX-SX2*SX2) - (SYX2*SX-SYX*SX2)
(1a2a)      b=  ((SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2)
(2b)        a*SX3*S0 +b*SX2*S0 + c*SX*S0  -SYX*S0  =0
(3b)        a*SX2*SX +b*SX*SX  + c*S0*SX  -SY*SX   =0
(2b-3b)     a*(SX3*S0-SX2*SX) +b*(SX2*S0-SX*SX)  -(SYX*S0-SY*SX)  =0
(2b3b1a2a)  a*(SX3*S0-SX2*SX) +((SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX)  -(SYX*S0-SY*SX)  =0
            a*(SX3*S0-SX2*SX) +(SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) -  a* (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX)  - (SYX*S0-SY*SX)  =0
            a*((SX3*S0-SX2*SX) - (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX)) + (SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) - (SYX*S0-SY*SX)  =0
            a=((SYX*S0-SY*SX)-(SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) )/((SX3*S0-SX2*SX) - (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX))
//(1a2a)   b=  ((SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2)
/**(3)      a*SX2 +b*SX  + c*S0  -SY   =0
            c =(SY -  a*SX2 - b*SX)/S0
*/
   
   if  (SX==0.0) {
     a= (SX2*SY - S0*SYX2)/ (SX2*SX2-S0*SX4);
   } else {
     a=((SYX*S0-SY*SX)-(SYX2*SX-SYX*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX) )/((SX3*S0-SX2*SX) - (SX4*SX-SX3*SX2)/(SX3*SX-SX2*SX2) * (SX2*S0-SX*SX));
   }
   b=  ((SYX2*SX-SYX*SX2) -  a* (SX4*SX-SX3*SX2))/(SX3*SX-SX2*SX2);
   c =(SY -  a*SX2 - b*SX)/S0; // Not needed here
   if (a!=0)  x=-b/(2*a);
   else x=0;
   crslt=0.5*(ilow+ihigh+2*x);
//Try to select a "good" result, use imax as a backup if all wrong
   if      ((ilow>0)                 && (crslt<ilow))              rslt= imx;  // solution outside of the range and the range does not include low end
   else if ((ihigh<(array.length-1)) && (crslt>ihigh))             rslt= imx;  // solution outside of the range and the range does not include high endScan
   else if ((ilow==0.0) && (crslt < (-ihigh)))                     rslt=ilow;  // maximum is too far outside of the range, can not trust
   else if ((ihigh==(array.length-1)) && (crslt > (2*ihigh-ilow))) rslt=ihigh; // maximum is too far outside of the range, can not trust
   else rslt= crslt;
   if (DEBUG_LEVEL>1 ) {
      sb.append("imx\t"+  imx+"\n");
      sb.append("ilow\t"+ ilow+"\n");
      sb.append("ihigh\t"+ihigh+"\n");
      sb.append("max\t"+  max+"\n");
      sb.append("wnd\t"+  wnd+"\n");
      sb.append("a\t"+   a+"\n");
      sb.append("b\t"+   b+"\n");
      sb.append("c\t"+   c+"\n");
      sb.append("SX4\t"+ SX4+"\n");
      sb.append("SX3\t"+ SX3+"\n");
      sb.append("SX2\t"+ SX2+"\n");
      sb.append("SX\t" + SX+ "\n");
      sb.append("SYX2\t"+SYX2+"\n");
      sb.append("SYX\t"+ SYX+"\n");
      sb.append("SY\t"+  SY+"\n");
      sb.append("S0\t"+  S0+"\n");
      sb.append("x\t"+   x+"\n");
      sb.append("crslt\t"+crslt+"\n");
      sb.append("rslt\t"+rslt+"\n");
      new TextWindow("findMaxInWndArray", header, sb.toString(), 350,800);
   }
   return rslt;
  }

double [][][][]calcMaxScanTargetOrientColor(double[][][][][] scanData, int wnd, double calibMaxSlopes) {
   double [][][][]rslt= new double[scanData[0].length][scanData[0][0].length][2][3];
   double [] arrayToMax=new double [scanData.length]; // number of samples
   int scan,t,n,orient;
   for (scan=0; scan<rslt.length; scan++) for (t=0; t<rslt[0].length;t++) for (orient=0;orient<2;orient++) { 
        for (n=0;n<scanData.length;n++) arrayToMax[n]=scanData[n][0][t][orient][1]; // red
        rslt[scan][t][orient][0]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);

        for (n=0;n<scanData.length;n++) arrayToMax[n]=0.5*(scanData[n][0][t][orient][0]+scanData[n][0][t][orient][3]); // (Gr+Gb)/2
        rslt[scan][t][orient][1]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
    
        for (n=0;n<scanData.length;n++) arrayToMax[n]=scanData[n][0][t][orient][2]; // blue
        rslt[scan][t][orient][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   } 
   return rslt;
}
/*
Goals
 1 - 4 - indepent for each target, 5..12 combine opposite targets around the center one (1 with 5, 2 with 4)
 1 - multiply Rv*Rh*Gv*Gh
 2 - Rv*Rh
 3 - Gv*Gh
 4 - between 2 and 3 in the (1-goalGreenWeight):goalGreenWeight
 5 - uses sqrt(Rv^2+Rh^2+Gv^2+Gh^2), targets_T2R to find center target, then for non-center target multiplies opposite arouind the center values
 6 - similar, but only sqrt(Rv^2+Rh^2)
 7 - similar, but only sqrt(Gv^2+Gh^2)
 8 - between 6 and 7 in the (1-goalGreenWeight):goalGreenWeight
 9 - uses sqrt(Rv^2+Rh^2+Gv^2+Gh^2), targets_T2R to find center target, then for non-center target sqrt(target1^2+target2^2) opposite arouind the center values
10 - similar, but only sqrt(Rv^2+Rh^2)
11 - similar, but only sqrt(Gv^2+Gh^2)
12 - between 10 and 11 in the (1-goalGreenWeight):goalGreenWeight



*/

double [][] makeOffsetsPerGoal(double[][][][][] scanData, double[] targets_T2R,int wnd, double calibMaxSlopes, double goalGreenWeight) {

   double [][] offsets = new double [12][scanData[0][0].length]; // number of targets
   double [] arrayToMax=new double [scanData.length]; // number of samples
   int i,i1,t,n;
   int goal;
   double [][][] targetData=new double[scanData.length][scanData[0][0].length][6]; // r,g,b 
   double a,b;
   if (scanData[0][0].length!=5) {
        IJ.showMessage("ERROR","makeOffsetsPerGoal only works with exactly 5 targets, there are "+scanData[0][0].length+" of them");
        return null;

   }
   for (n=0;n<scanData.length;n++) for (t=0; t<scanData[0][0].length;t++) { 
        targetData[n][t][0]= scanData[n][0][t][0][1]; // R, vertical
        targetData[n][t][2]=0.5*(scanData[n][0][t][0][0]+scanData[n][0][t][0][3]); // G, vertical
        targetData[n][t][4]= scanData[n][0][t][0][2]; // B, vertical

        targetData[n][t][1]= scanData[n][0][t][1][1]; // R, horizontal
        targetData[n][t][3]=0.5*(scanData[n][0][t][1][0]+scanData[n][0][t][1][3]); // G, horizontal
        targetData[n][t][5]= scanData[n][0][t][1][2]; // B, horizontal
    } 

// goal # 1 - Rv*Rh*Gv*Gh
   goal=1;
   for (t=0;t<targetData[0].length; t++ ) {
     for (i=0;i<targetData.length;i++) arrayToMax[i]=targetData[i][t][0]*targetData[i][t][1]*targetData[i][t][2]*targetData[i][t][3];
     offsets [goal-1][t]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   }
// goal # 2 - Rv*Rh
   goal=2;
   for (t=0;t<targetData[0].length; t++ ) {
     for (i=0;i<targetData.length;i++) arrayToMax[i]=targetData[i][t][0]*targetData[i][t][1];
     offsets [goal-1][t]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   }
// goal # 3 - Gv*Gh
   goal=3;
   for (t=0;t<targetData[0].length; t++ ) {
     for (i=0;i<targetData.length;i++) arrayToMax[i]=targetData[i][t][2]*targetData[i][t][3];
     offsets[goal-1][t]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   }

// goal # 4 -between #2 and #3
   goal=4;
   for (t=0;t<targetData[0].length; t++ ) {
     offsets[goal-1][t]=offsets[goal-2][t]*goalGreenWeight+offsets[goal-3][t]*(1.0-goalGreenWeight);
   }


// goal # 5..#6 find center on target 2 (center) then multiply 0*4 and 1*3 symmetrically around the center
// # 4 - both g and r, #5 - R only, #6 - G only
// -sqrt( targets_T2R[0]*Rv[0]^2+ (1-targets_T2R[0])*Rh^2+ targets_T2R[0]*Gv[0]^2+ (1-targets_T2R[0])*Gh^2+)
   goal=5;
   t=2;
   for (i=0;i<targetData.length;i++) {
         arrayToMax[i]=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i][t][0]*targetData[i][t][0]+targetData[i][t][2]*targetData[i][t][2])+
                                                                       (targetData[i][t][1]*targetData[i][t][1]+targetData[i][t][3]*targetData[i][t][3]));
   }
   offsets[goal-1][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t  ][0]*targetData[i ][t][0]+targetData[i ][t][2]*targetData[i ][t][2])+
                                                                        (targetData[i ][t  ][1]*targetData[i ][t][1]+targetData[i ][t][3]*targetData[i ][t][3]));
         b=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][0]*targetData[i1][4-t][0]+targetData[i1][4-t][2]*targetData[i1][4-t][2])+
                                                                        (targetData[i1][4-t][1]*targetData[i1][4-t][1]+targetData[i1][4-t][3]*targetData[i1][4-t][3]));
         arrayToMax[i]=a*b;
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t];
   }
   goal=6;
   t=2;
   for (i=0;i<targetData.length;i++) {
         arrayToMax[i]=1.0/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i][t][0]*targetData[i][t][0])+
                                                                       (targetData[i][t][1]*targetData[i][t][1]));
   }
   offsets[goal-1][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][0]*targetData[i ][t][0])+
                                                                        (targetData[i ][t][1]*targetData[i ][t][1]));
         b= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][0]*targetData[i1][4-t][0])+
                                                                        (targetData[i1][4-t][1]*targetData[i1][4-t][1]));
         arrayToMax[i]=a*b;
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t];
   }

   goal=7;
   t=2;
   for (i=0;i<targetData.length;i++) {
         arrayToMax[i]=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i][t][2]*targetData[i][t][2])+
                                                                       (targetData[i][t][3]*targetData[i][t][3]));
   }
   offsets[goal-1][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][2]*targetData[i ][t][2])+
                                                                        (targetData[i ][t][3]*targetData[i ][t][3]));
         b=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][2]*targetData[i1][4-t][2])+
                                                                        (targetData[i1][4-t][3]*targetData[i1][4-t][3]));
         arrayToMax[i]=a*b;
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t];
   }



// goal # 8 -between #6 and #7
   goal=8;
// for center target use 6 & 7
   t=2; // center target
   offsets[goal-1][t]=offsets[goal-2][t]*goalGreenWeight+offsets[goal-3][t]*(1.0-goalGreenWeight); 

   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][0]*targetData[i ][t][0])+ // Rv
                                                                        (targetData[i ][t][1]*targetData[i ][t][1])); //Rh
         b= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][0]*targetData[i1][4-t][0])+
                                                                        (targetData[i1][4-t][1]*targetData[i1][4-t][1]));
         arrayToMax[i]=a*b;
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes); // use for red

     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][2]*targetData[i ][t][2])+ //Gv
                                                                        (targetData[i ][t][3]*targetData[i ][t][3])); //Gh
         b=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][2]*targetData[i1][4-t][2])+
                                                                        (targetData[i1][4-t][3]*targetData[i1][4-t][3]));
         arrayToMax[i]=a*b;
       }
     }
     offsets[goal-1][4-t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes); // temporarily use that for green
     offsets[goal-1][t]= offsets[goal-1][4-t]*goalGreenWeight+  offsets[goal-1][t]*(1.0-goalGreenWeight); // between red and green
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t]; // copy to opposite target
   }








   goal=9;
   t=2;
   for (i=0;i<targetData.length;i++) {
         arrayToMax[i]=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i][t][0]*targetData[i][t][0]+targetData[i][t][2]*targetData[i][t][2])+
                                                                       (targetData[i][t][1]*targetData[i][t][1]+targetData[i][t][3]*targetData[i][t][3]));
   }
   offsets[goal-1][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t  ][0]*targetData[i ][t][0]+targetData[i ][t][2]*targetData[i ][t][2])+
                                                                        (targetData[i ][t  ][1]*targetData[i ][t][1]+targetData[i ][t][3]*targetData[i ][t][3]));
         b=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][0]*targetData[i1][4-t][0]+targetData[i1][4-t][2]*targetData[i1][4-t][2])+
                                                                        (targetData[i1][4-t][1]*targetData[i1][4-t][1]+targetData[i1][4-t][3]*targetData[i1][4-t][3]));
         arrayToMax[i]=Math.sqrt(a*a+b*b);
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t];
   }
   goal=10;
   t=2;
   for (i=0;i<targetData.length;i++) {
         arrayToMax[i]=1.0/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i][t][0]*targetData[i][t][0])+
                                                                       (targetData[i][t][1]*targetData[i][t][1]));
   }
   offsets[goal-1][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][0]*targetData[i ][t][0])+
                                                                        (targetData[i ][t][1]*targetData[i ][t][1]));
         b= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][0]*targetData[i1][4-t][0])+
                                                                        (targetData[i1][4-t][1]*targetData[i1][4-t][1]));
         arrayToMax[i]=Math.sqrt(a*a+b*b);
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t];
   }
   goal=11;
   t=2;
   for (i=0;i<targetData.length;i++) {
         arrayToMax[i]=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i][t][2]*targetData[i][t][2])+
                                                                       (targetData[i][t][3]*targetData[i][t][3]));
   }
   offsets[goal-1][2]=findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][2]*targetData[i ][t][2])+
                                                                        (targetData[i ][t][3]*targetData[i ][t][3]));
         b=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][2]*targetData[i1][4-t][2])+
                                                                        (targetData[i1][4-t][3]*targetData[i1][4-t][3]));
         arrayToMax[i]=Math.sqrt(a*a+b*b);
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes);
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t];
   }

// goal # 12 -between #10 and #11
   goal=12;
// for center target use 10 & 11
   t=2; // center target
   offsets[goal-1][t]=offsets[goal-2][t]*goalGreenWeight+offsets[goal-3][t]*(1.0-goalGreenWeight); 

   for (t=0;t<2;t++) {
     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][0]*targetData[i ][t][0])+ // Rv
                                                                        (targetData[i ][t][1]*targetData[i ][t][1])); //Rh
         b= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][0]*targetData[i1][4-t][0])+
                                                                        (targetData[i1][4-t][1]*targetData[i1][4-t][1]));
         arrayToMax[i]=Math.sqrt(a*a+b*b);
       }
     }
     offsets [goal-1][t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes); // use for red

     for (i=0;i<targetData.length;i++) {
       i1=(int)(2*offsets[goal-1][2])-i;
       if ((i1<0) || (i1>=targetData.length)) {
         arrayToMax[i]=0.0;
       } else {
         a= 0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i ][t][2]*targetData[i ][t][2])+ //Gv
                                                                        (targetData[i ][t][3]*targetData[i ][t][3])); //Gh
         b=0.5/(targets_T2R[t]+1)*Math.sqrt(targets_T2R[2]*(targetData[i1][4-t][2]*targetData[i1][4-t][2])+
                                                                        (targetData[i1][4-t][3]*targetData[i1][4-t][3]));
         arrayToMax[i]=Math.sqrt(a*a+b*b);
       }
     }
     offsets[goal-1][4-t]  =findMaxInWndArray(arrayToMax, wnd, calibMaxSlopes); // temporarily use that for green
     offsets[goal-1][t]= offsets[goal-1][4-t]*goalGreenWeight+  offsets[goal-1][t]*(1.0-goalGreenWeight); // between red and green
     offsets [goal-1][4-t]=2*offsets[goal-1][2] - offsets [goal-1][t]; // copy to opposite target
   }
   return offsets;  
}

 // find index (distance from focus) for using calibration data and measurement of sharpness for different colors.
  double [] findDisplacemntsFromResolutions (double [][][] targetScanData,
                                            double[][][] resolutions, // resolutions for one target [v.h][Gr,R,B,Gb]
                                            int sideWays, // when after finding best correlation index, use two points sideWays above/below to interpolate
                                            int    motors_calibSteps,
                                            double stepsPerUM,
                                            double [] offsetsPerGoal,
                                            double focusOffset
                                            ) {
    double [] rslt=new double[resolutions.length];
    int t;
    for (t=0;t<resolutions.length;t++) rslt[t]=findDisplacemntFromResolutions (targetScanData,
                                         resolutions[t], // resolutions for one target [v.h][Gr,R,B,Gb]
                                         t, // target number
                                         sideWays,
                                         motors_calibSteps,
                                         stepsPerUM,
                                         offsetsPerGoal[t],
                                         focusOffset);
    return rslt;
  }

  double findDisplacemntFromResolutions (double [][][] targetScanData,
                                         double[][] resolutions, // resolutions for one target [v.h][Gr,R,B,Gb]
                                         int target, // target number
                                         int sideWays, // when after finding best correlation index, use two points sideWays above/below to interpolate
                                         int    motors_calibSteps,
                                         double stepsPerUM,
                                         double offsetsPerGoal,
                                         double focusOffset
                                        ){
    int i,n,d,nmin, nlower,nhigher;
    double [] sample=new double[7];
    double diffErr;
    double minErr=0;
    double dindex,rslt;
    if (targetScanData==null) {
        IJ.showMessage("ERROR","Resolutions vs motor positions are not calibrated yet.");
        return 0.0; // or is it ok to fail?
    }
    String header="#\tDiff\tRv\tRh\tGv\tGh\tBv\tBh\tsum";
    StringBuffer sb = new StringBuffer();
    for (d=0;d<2;d++) {
                        sample[d+0]=      resolutions[d][1];
                        sample[d+2]= 0.5*(resolutions[d][0]+resolutions[d][3]);
                        sample[d+4]=      resolutions[d][2];
    }
    sample[6]=0.0;
    for (i=0;i<6;i++) sample[6]+=sample[i]*sample[i];
    sample[6]=Math.sqrt(sample[6]);
    for (i=0;i<6;i++) sample[i]/=sample[6];
    if (DEBUG_LEVEL>2 )  sb.append("sample\t---\t"+sample[0]+"\t"+sample[1]+"\t"+sample[2]+"\t"+sample[3]+"\t"+sample[4]+"\t"+sample[5]+"\t"+sample[6]+"\n");
// find most likely index by minimizing difference between sample resolution vector and element of the array of indexed resolution vectors
    nmin=0;
    for (n=0; n< targetScanData.length; n++) {
      diffErr=0;
      for (i=0;i<6;i++) diffErr+=(sample[i]-targetScanData[n][target][i])*(sample[i]-targetScanData[n][target][i]);
      if ((diffErr<minErr)|| (n==0)) {
        nmin=n;
        minErr=diffErr;
      }
      if (DEBUG_LEVEL>2 )  sb.append(n+"\t"+diffErr+
                                       "\t"+targetScanData[n][target][0]+
                                       "\t"+targetScanData[n][target][1]+
                                       "\t"+targetScanData[n][target][2]+
                                       "\t"+targetScanData[n][target][3]+
                                       "\t"+targetScanData[n][target][4]+
                                       "\t"+targetScanData[n][target][5]+
                                       "\t"+targetScanData[n][target][6]+"\n");

    }
// Stage 2 - refine distance by taking 2 elements in the array above and below the found minimum, then compare using the most distinct components (i.e. in the area where only blue vertical resolution changes - use blue vertical resolution)
    nlower=  nmin-sideWays;
    nhigher= nmin+ sideWays;
    if (nlower<0) nlower=0;
    if (nhigher >=targetScanData.length) nhigher=targetScanData.length-1;
    diffErr=0.0;
    for (i=0;i<6;i++) diffErr+=(targetScanData[nhigher][target][i]-targetScanData[nlower][target][i])*(targetScanData[nhigher][target][i]-targetScanData[nlower][target][i]);
    dindex= 0.0;
    for (i=0;i<6;i++) dindex+=(sample[i]-targetScanData[nlower][target][i])*(targetScanData[nhigher][target][i]-targetScanData[nlower][target][i]);
    dindex= nlower+(nhigher-nlower)*dindex/diffErr;
    rslt= (dindex-offsetsPerGoal)*motors_calibSteps/stepsPerUM - focusOffset;

    if (DEBUG_LEVEL>2 ) {
       sb.append(nmin+"\tmin="+minErr+"\t\t\t\t\t\t\t\t\n");
       sb.append(nlower+"\tnlower\t\t\t\t\t\t\t\t\n");
       sb.append(nhigher+"\tnhigher\t\t\t\t\t\t\t\t\n");
       sb.append(IJ.d2s(dindex,4)+"\tdindex\t\t\t\t\t\t\t\n");
       sb.append(IJ.d2s(offsetsPerGoal,4)+"\toffsetsPerGoal\t\t\t\t\t\t\t\n");
       sb.append(motors_calibSteps+"\tmotors_calibSteps\t\t\t\t\t\t\t\n");
       sb.append(IJ.d2s(stepsPerUM,1)+"\tstepsPerUM\t\t\t\t\t\t\t\n");
       sb.append(IJ.d2s(focusOffset,2)+"\tfocusOffset\t\t\t\t\t\t\t\n");
       sb.append(IJ.d2s(rslt,1)+" um\tresult\t\t\t\t\t\t\t\n");
       new TextWindow("findDisplacemntFromResolutions, target="+(target+1), header, sb.toString(), 1200,800);
    }
//    return dindex;
    return rslt;
  }

double calcGoalDiff(double[][]resolutions, // resolutions for one target [v.h][Gr,R,B,Gb]
                    double goalScale,      // Signed, if positive (tragets 2,4) target ratio V/H resolutions. For center target (if colorMask==0) it is target G/R ratio.
                                           // negative goalScale sets -H/V and reverses sign of the result
                    int colorMask) { //bitmask of colors to compare vertical
    int i;    
    double sumPos=0.0;
    double sumNeg=0.0;
    double diff;

    if (colorMask!=0) {
          for (i=0;i<3;i++)  {
            sumPos+=resolutions[0][i]*((colorMask>>i)&1); // vertical
            sumNeg+=resolutions[1][i]*((colorMask>>i)&1); // horizontal
          }
    } else {
      sumPos=0.5*(resolutions[0][0]+resolutions[1][0]+resolutions[0][3]+resolutions[1][3]); // greens
      sumNeg=resolutions[0][1]+resolutions[1][1]; //red
    }
    if (goalScale>0) {
      diff=(sumPos-goalScale*sumNeg)/(sumPos+goalScale*sumNeg);
    } else {
      diff=(sumNeg+goalScale*sumPos)/(sumNeg-goalScale*sumPos);
    }
    return diff;
}

// Measure four off-center targets, adjusting M3, return a pair of differences M3(T2)-M3(T4) and M3(T1) - M3(T5)

  int []  measureM3Deltas( // Returns null if reached motor limits, otherwise focus results for all colors/targets
                                String        motorsIP,
                                int           signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                                double        motorsWaitAfterMove,
                                int           motorNumber,
                                int           minStep, // minimal motor rotation - positive
                                int           maxStep, // maximal motor rotation (positive)
                                int        [] targets_dirs,
                                int        [] ColorMasks,
                                double     [] V2H,
//                                int           targetNumber,
//                                int           colorMask, // bitmask of colors to average: 2 - red, 4 - blue, 9 - Gr+Gb , 0 - use (G-R) instead - for central target
//                                double        V2H, // target ratio of vertial/horizontal resolutions for positive decrDirection, H/V - for negative. Gree/Red for colorMask==0
                                int           refine, /// number of minimal steps made to refine the final position (will be half each way
                                int           refineAverage,
                                int           refineExpand,
                                int           maxSteps) {// - just for debugging - maximal number of steps to try before giving up
    return measureM3Deltas(motorsIP,signed_hysteresis,motorsWaitAfterMove, motorNumber,minStep, maxStep, targets_dirs, ColorMasks, V2H, refine,refineAverage, refineExpand, maxSteps, true);

  }
  int []  measureM3Deltas( // Returns null if reached motor limits, otherwise focus results for all colors/targets
                                String        motorsIP,
                                int           signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                                double        motorsWaitAfterMove,
                                int           motorNumber,
                                int           minStep, // minimal motor rotation - positive
                                int           maxStep, // maximal motor rotation (positive)
                                int        [] targets_dirs,
                                int        [] ColorMasks,
                                double     [] V2H,
//                                int           targetNumber,
//                                int           colorMask, // bitmask of colors to average: 2 - red, 4 - blue, 9 - Gr+Gb , 0 - use (G-R) instead - for central target
//                                double        V2H, // target ratio of vertial/horizontal resolutions for positive decrDirection, H/V - for negative. Gree/Red for colorMask==0
                                int           refine, /// number of minimal steps made to refine the final position (will be half each way
                                int           refineAverage,
                                int           refineExpand,
                                int           maxSteps,
                                boolean       deltas) {
  int [] m3=new int[5];
  int [] startPosition=readElphel10364Motors(motorsIP);
  int [] currentPosition;
  int [] result=new int[2];
  int    targetNumber;

  for (targetNumber=0;targetNumber<5;targetNumber++) if (!deltas || (targetNumber != 2)){ // skip center when looking for just deltas
     if (      findThirdMotor( // Returns null if reached motor limits, otherwise focus results for all colors/targets
                                imp_camera, // global
                                motorsIP,
                                signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                                motorsWaitAfterMove,
                                motorNumber,
                                minStep, // minimal motor rotation - positive
                                maxStep, // maximal motor rotation (positive)
                                targetNumber,
                                ColorMasks[targetNumber], // bitmask of colors to average: 2 - red, 4 - blue, 9 - Gr+Gb , 0 - use (G-R) instead - for central target
                                V2H[targetNumber], // target ratio of vertial/horizontal resolutions for positive decrDirection, H/V - for negative. Gree/Red for colorMask==0
                                refine, /// number of minimal steps made to refine the final position (will be half each way
                                refineAverage,
                                refineExpand,
                                maxSteps) == null) {
      moveElphel10364Motors(motorsIP, startPosition, 0, true, "return",true); /// restore initial position

      return null;
    } else {
      currentPosition=readElphel10364Motors(motorsIP);
      m3[targetNumber]= currentPosition[motorNumber];
    }
  }
  result[0]=m3[1]-m3[3];
  result[1]=m3[0]-m3[4];
  if (deltas)  return result;
  return       m3; // all 5
}



/// The goal is to adjust motors M1 and M2 so difference betwee value of in-focus M3 for two pairs of the oppposite targets is simultaneously 0
/// Targets are considered to be "in-focus" when resolution for tangential direction is the same as the radial one (they distinctlyt cross)
/// It is possible to match those resolutions with scale, that should work to focus on 2m target for the infinity

  public boolean    adjustAllThree ( String        motorsIP,
                                     int           signed_hysteresis,  // Make it per-motor?
                                     double        motorsWaitAfterMove,
                                     int        [] motor_numbers,
                                     int           initialStep,
                                     int           numIterations,
                                     int        [] targets_dirs,
                                     int        [] ColorMasks,
                                     double     [] V2H,
                                     int           minStep,
                                     int           maxStep, 
                                     int           refine, /// number of minimal steps made to refine the final position (will be half each way
                                     int           refineAverage,
                                     int           refineExpand,
                                     int           maxSteps) { // also applies to initial search
    String header="#\tM1\tM2\tM3\tdelta[0]\tdelta[1]\td1m1\td1m2\td2m1\td2m2\tdm1\tdm2\tmessage";
    StringBuffer sb = new StringBuffer();
//        if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\t---\t---\n");
//    if (DEBUG_LEVEL>1 ) {
//    String header="#\tM1\tM2\tM3\td1m1\td1m2\td2m1\td2m2\tdm1\tdm2\tmessage";
//       TextWindow tw = new TextWindow("motor iterations - TARGET="+(targetNumber+1)+" COLORS="+colorMask+" V2H="+V2H, header, sb.toString(), 1200,800);
//    }


    int [] startPosition;
    int [] targetPosition=new int[3];
    int M1_signed_hysteresis=signed_hysteresis; // split them?
    int M2_signed_hysteresis=signed_hysteresis;
    int M3_signed_hysteresis=signed_hysteresis;
    double d1m1=0.0;
    double d2m1=0.0;
    double d1m2=0.0;
    double d2m2=0.0;
    double dm1=0.0;
    double dm2=0.0;
    boolean go=true;
    boolean abort=false;
    int     numberInside=0; // numer of times new center is inside previous borders (for nuw - second time ends initial iterations)
    int stepNum=0;
    String errorMsg=new String();
    int [] deltas=null;
    int [] m3s;
    int iterNum=0;
//    int [][] fourCorners=new int [4][2];
    int [][] fourCorners=new int [4][];
    int i;
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error (in thirdByFirstTwo)","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return false;
    }
    startPosition=readElphel10364Motors(motorsIP);
    for (i=0;i<3;i++) targetPosition[i]=startPosition[i];
// first stage - trying to get square that includes the zero-zero point  
// public double [][]    motorsRange={{-5000,5000},{-5000,5000},{-5000,5000}};

//    while (go && !abort && (stepNum<maxSteps)) {
    while (go && !abort && (stepNum<numIterations)) {

      if (targetPosition[motor_numbers[0]]>(motorsRange[motor_numbers[0]][1]-initialStep)) {
         errorMsg+="Motor "+(motor_numbers[0]+1)+" exceeds positive limit (arm too close): "+ targetPosition[motor_numbers[0]] +">"+(motorsRange[motor_numbers[0]][1]-initialStep)+"\n";
         abort=true;
      }
      if (targetPosition[motor_numbers[1]]>(motorsRange[motor_numbers[1]][1]-initialStep)) {
         errorMsg+="Motor "+(motor_numbers[1]+1)+" exceeds positive limit (arm too close): "+ targetPosition[motor_numbers[1]] +">"+(motorsRange[motor_numbers[1]][1]-initialStep)+"\n";
         abort=true;
      }
      if (targetPosition[motor_numbers[0]]<(motorsRange[motor_numbers[0]][0]+initialStep)) {
         errorMsg+="Motor "+(motor_numbers[0]+1)+" exceeds negative limit (arm too far): "  + targetPosition[motor_numbers[0]] +"<"+(motorsRange[motor_numbers[0]][0]+initialStep)+"\n";
         abort=true;
      }
      if (targetPosition[motor_numbers[1]]<(motorsRange[motor_numbers[1]][0]+initialStep)) {
         errorMsg+="Motor "+(motor_numbers[1]+1)+" exceeds negative limit (arm too far): "  +targetPosition[motor_numbers[1]]  +"<"+(motorsRange[motor_numbers[1]][0]+initialStep)+"\n";
         abort=true;
      }
      if (abort) break;
      targetPosition[motor_numbers[0]]-=initialStep;
      if (M1_signed_hysteresis>0) { // Make anti-hysteresis move
        targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
        moveElphel10364Motors(motorsIP, targetPosition, 0, true);
        targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
      }
      moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+10, true, "M1:- M2:0");
      fourCorners[0]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
      if (fourCorners[0]==null) abort=true;
      if (!abort) {
        targetPosition[motor_numbers[0]]+=2*initialStep;
        if (M1_signed_hysteresis<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:+ M2:0");
        fourCorners[1]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        if (fourCorners[1]==null) abort=true;
      }
      if (!abort) {
        targetPosition[motor_numbers[0]]-=initialStep;
        if (M1_signed_hysteresis>0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
        }
        targetPosition[motor_numbers[1]]-=initialStep;
        if (M2_signed_hysteresis>0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:0 M2:-");
        fourCorners[2]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        if (fourCorners[2]==null) abort=true;
      }
      if (!abort) {
        targetPosition[motor_numbers[1]]+=2*initialStep;
        if (M2_signed_hysteresis<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:0 M2:+");
        fourCorners[3]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        if (fourCorners[3]==null) abort=true;
      }
// return to the center
      if (!abort) {
        targetPosition[motor_numbers[1]]-=initialStep;
        if (M2_signed_hysteresis>0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:0 M2:0");
        deltas=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
      }
      if (abort) {
         errorMsg+="Motor "+(motor_numbers[2]+1)+" probably hit one of the limits trying to compensate focus position for one of the targets\n";
         break;
      }
      IJ.showMessage("Debug","{{"+fourCorners[0][0]+"},{"+fourCorners[0][1]+"}},\n"+
                              "{"+fourCorners[1][0]+"},{"+fourCorners[1][1]+"}},\n"+
                              "{"+fourCorners[2][0]+"},{"+fourCorners[2][1]+"}},\n"+
                              "{"+fourCorners[3][0]+"},{"+fourCorners[3][1]+"}}");


/// Find double zero
//    double d1m1,d2m1,d1m2,d2m2;
      d1m1=(1.0*(fourCorners[1][0]-fourCorners[0][0]))/(2.0*initialStep);
      d2m1=(1.0*(fourCorners[1][1]-fourCorners[0][1]))/(2.0*initialStep);
      d1m2=(1.0*(fourCorners[3][0]-fourCorners[2][0]))/(2.0*initialStep);
      d2m2=(1.0*(fourCorners[3][1]-fourCorners[2][1]))/(2.0*initialStep);

      dm1 =(d1m2*deltas[1]-d2m2*deltas[0])/(d2m2*d1m1-d1m2*d2m1);
      dm2 =(d1m1*deltas[1]-d2m1*deltas[0])/(d2m1*d1m2-d1m1*d2m2);
      if ((Math.abs(dm1)>initialStep) || (Math.abs(dm1)>initialStep)) { // desired location outside of the square
        if      (dm1 >  initialStep) dm1= initialStep;
        else if (dm1 < -initialStep) dm1= -initialStep;
        if      (dm2 >  initialStep) dm2= initialStep;
        else if (dm2 < -initialStep) dm2= -initialStep;
      } else {
        numberInside++;
        if (numberInside > 1) go = false; // Stage 1 done
      }
// Make step
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+deltas[0]+"\t"+deltas[1]+"\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tStage 1\n");

      targetPosition[motor_numbers[0]]+=dm1;
      if ((M1_signed_hysteresis*dm1)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
      }
      targetPosition[motor_numbers[1]]+=dm2;
      if ((M2_signed_hysteresis*dm2)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
      }
      moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1+"+dm1+" M2+"+dm2+"STEP:"+stepNum);
      stepNum++;
//      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+delta[0]+"\t"+delta[1]+"\t"+"\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tmade step\n");
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t-\t-\t-\t-\t-\t-\t-\t-\t-\tmade step\n");
    }
/// Stage 2 - uae previously found d1m1,d1m2,d2m1,d2m2 to find target relative to current position, exit at certain step
    for (iterNum=0;!abort && (iterNum<numIterations);iterNum++) {
      deltas=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
      dm1 =(d1m2*deltas[1]-d2m2*deltas[0])/(d2m2*d1m1-d1m2*d2m1);
      dm2 =(d1m1*deltas[1]-d2m1*deltas[0])/(d2m1*d1m2-d1m1*d2m2);
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+deltas[0]+"\t"+deltas[1]+"\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tStage 2\n");
// Make step
      targetPosition[motor_numbers[0]]+=dm1;
      if ((M1_signed_hysteresis*dm1)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
      }
      targetPosition[motor_numbers[1]]+=dm2;
      if ((M2_signed_hysteresis*dm2)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
      }
      moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1+"+dm1+" M2+"+dm2+"iterNum:"+iterNum);
//      if (DEBUG_LEVEL>0) sb.append(iterNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tStage 2\n");
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t-\t-\t-\t-\t-\t-\t-\t-\t-\tmade step\n");

    }
    if (DEBUG_LEVEL>0 ) {
       m3s=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps,false);
       sb.append("Final\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\tt1="+m3s[0]+"\tt2="+m3s[1]+"\tt3="+m3s[2]+"\tt4="+m3s[3]+"\tt5="+m3s[4]+"\t-\t-\t-\t-\tmade step\n");
       new TextWindow("Focus adjustmet", header, sb.toString(), 1200,800);
    }
    if (abort) {
      moveElphel10364Motors(motorsIP, startPosition, 0, true, "abort - return",true); /// restore initial position
      IJ.showMessage("Error",errorMsg);
      return false;
    }
    return true;
  }



  public int [][][] thirdByFirstTwo (ImagePlus    imp_src,
                                     String        motorsIP,
                                     int           signed_hysteresis,  // Make it per-motor?
                                     double        motorsWaitAfterMove,
                                     int        [] motor_numbers,
                                     int           M1range,
                                     int           M1steps,
                                     int           M2range,
                                     int           M2steps,
                                     int        [] targets_dirs,
                                     int        [] ColorMasks,
                                     double     [] V2H,
                                     int           minStep,
                                     int           maxStep, 
//                                     double [][][] oldResult,
                                     int           refine, /// number of minimal steps made to refine the final position (will be half each way
                                     int           refineAverage,
                                     int           refineExpand,
                                     int           maxSteps) {
    int [][][] result = new int [M1steps+1][M2steps+1][targets_dirs.length+2];
    int i,M1_index,M2_index;
    int [] currentPosition;
    int [] startPosition;
    int [] targetPosition=new int[3];
    int M1_signed_hysteresis=signed_hysteresis; // split them?
    int M2_signed_hysteresis=signed_hysteresis;
    int M3_signed_hysteresis=signed_hysteresis;
//    int motor_numbers={0,1,2};
    int targetNumber;
    double [][][] resolutions;
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error (in thirdByFirstTwo)","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return null;
    }
    startPosition=readElphel10364Motors(motorsIP);
    for (i=0;i<3;i++) targetPosition[i]=startPosition[i];
    for (M1_index=0; M1_index<=M1steps; M1_index++) {
//      targetPosition[motor_numbers[0]]=startPosition[motor_numbers[0]]+(((M1_signed_hysteresis>0)?1:-1)*(M1range*(2*M1_index-M1steps-1))/(M1steps+1));
      targetPosition[motor_numbers[0]]=startPosition[motor_numbers[0]]+(((M1_signed_hysteresis>0)?1:-1)*(M2range*(2*M1_index-M1steps))/M1steps/2);
      currentPosition=readElphel10364Motors(motorsIP);
      targetPosition[motor_numbers[1]]=currentPosition[motor_numbers[1]]; // so other motors will not move
      targetPosition[motor_numbers[2]]=currentPosition[motor_numbers[2]]; // so other motors will not move
      if (M1_index==0) { // Make anti-hysteresis move
        targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
        moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
        targetPosition[motor_numbers[0]]+=M1_signed_hysteresis; // restore back
      }
      moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here - not really needed
      for (M2_index=0; M2_index<=M2steps; M2_index++) {
        currentPosition=readElphel10364Motors(motorsIP);
        targetPosition[motor_numbers[2]]=currentPosition[motor_numbers[2]]; // so other motor will not move
//        targetPosition[motor_numbers[1]]=startPosition[motor_numbers[1]]+(((M2_signed_hysteresis>0)?1:-1)*(M2range*(2*M2_index-M2steps-1))/(M2steps+1));
        targetPosition[motor_numbers[1]]=startPosition[motor_numbers[1]]+(((M2_signed_hysteresis>0)?1:-1)*(M2range*(2*M2_index-M2steps))/M2steps/2);
        if (M2_index==0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis; // restore back
        }
        currentPosition=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, M1_index+":"+M2_index);
        result[M1_index][M2_index][0]=currentPosition[motor_numbers[0]]; // actual position
        result[M1_index][M2_index][1]=currentPosition[motor_numbers[1]]; // actual position
        for (targetNumber=0; targetNumber<targets_dirs.length; targetNumber++) {
          if (targets_dirs[targetNumber]!=0) {
            resolutions= findThirdMotor(
                   imp_camera,
                   motorsIP,
                   M3_signed_hysteresis,
                   motorsWaitAfterMove,
                   motor_numbers[2],
                   motor3rd_minStep,
                   motor3rd_maxStep,
//                   motor3rd_direction,
                   targetNumber,
                   ColorMasks[targetNumber],
                   V2H[targetNumber]* targets_dirs[targetNumber],
//                   null,
                   refine,
                   refineAverage,
                   refineExpand,
                   maxSteps);
            if (resolutions==null) result[M1_index][M2_index][targetNumber+2]=MOTOR_UNDEFINED;
            else {
              currentPosition=readElphel10364Motors(motorsIP);
              result[M1_index][M2_index][targetNumber+2]=currentPosition[motor_numbers[2]];
            }      
          } else { /// skip this target_filter
            result[M1_index][M2_index][targetNumber+2]=MOTOR_UNDEFINED;
          }
        }
      }
    }
    moveElphel10364Motors(motorsIP, startPosition, 0, true, "return",true); /// restore initial position
    return result;
  }

  public void showThirdByFirstTwoTargets(int [][][] table, int[] motor_numbers, int []  targets_dirs, String title) {
    String header="#\ti\tj\tM"+(motor_numbers[0]+1)+"\tM"+(motor_numbers[1]+1); // +"\tM3";
    StringBuffer sb = new StringBuffer();
    int target,m1,m2,i;
    int lineNumber=0;
    for (target=2;target<table[0][0].length;target++) if (targets_dirs[target-2]!=0) {
      header+="\tM"+(motor_numbers[2]+1)+" T"+(target-1);
    }
    for (m1=0;m1<table.length;m1++) for (m2=0;m2<table[m1].length;m2++) {
      lineNumber++;
      sb.append(lineNumber+
                "\t"+m1+
                "\t"+m2);
      for (i=0;i<table[m1][m2].length; i++)  if ((i<2) || (targets_dirs[i-2]!=0)) sb.append("\t"+((table[m1][m2][i]==MOTOR_UNDEFINED)?"---":table[m1][m2][i]));
      sb.append("\n");
    }
    new TextWindow(title+"(targets)", header, sb.toString(), 1200,800);
  }


  public void showThirdByFirstTwo(int [][][] table, int[] motor_numbers, int []  targets_dirs, String title) {
    String header="t\t---\ti";
    StringBuffer sb = new StringBuffer();
    int target,m1,m2;

    for (m1=0;m1<table.length;m1++)  header+="\t"+m1;
    for (target=0;target<targets_dirs.length;target++) if (targets_dirs[target]!=0) {
      sb.append(header+"\n");
      sb.append((target+1)+"\tj\tM"+(motor_numbers[1]+1)+"\\M"+(motor_numbers[0]+1));
      for (m1=0;m1<table.length;m1++)  sb.append("\t"+table[m1][table[m1].length>>1][0]); // M1 positions for the M2 in the middle
      sb.append("\n");
      for (m2=0;m2<table[0].length;m2++) {
        sb.append((target+1)+"\t"+m2+"\t"+table[table.length>>1][m2][1]);
        for (m1=0;m1<table.length;m1++) {
          sb.append("\t"+((table[m1][m2][target+2]==MOTOR_UNDEFINED)?"---":table[m1][m2][target+2])); ///M3 value
        }
        sb.append("\n");
      }
    }
    new TextWindow(title, header, sb.toString(), 1200,800);
  }




/// The goal is to adjust motors M1 and M2 so difference betwee value of in-focus M3 for two pairs of the oppposite targets is simultaneously 0
/// Targets are considered to be "in-focus" when resolution for tangential direction is the same as the radial one (they distinctlyt cross)
/// It is possible to match those resolutions with scale, that should work to focus on 2m target for the infinity

  public boolean    adjustAllThreeCalib(
                                     String        motorsIP,
                                     int           signed_hysteresis,  // Make it per-motor?
                                     double        motorsWaitAfterMove,
                                     int        [] motor_numbers,
                                     int           initialStep,
                                     int           numIterations,
                                     int           minStep,
                                     int           average, // number of measurements to average
                                     int           minUpdateStep, // ipdate thisStepsPerUM only if the last step was > this threshold
                                     double        variableStepRatio, // 2.0 - maximal thisStepsPerUM change in one step
                                     int           fastSteps, // use reduced number of steps to adjust third motor on the initial stage =2?
                                     int           maxSteps) {// - just for debugging - maximal number of steps to try before giving up
    String header="#\tM1\tM2\tM3\tdelta[0]\tdelta[1]\td1m1\td1m2\td2m1\td2m2\tdm1\tdm2\tmessage";
    StringBuffer sb = new StringBuffer();


    int [] startPosition;
    int [] lastPos=new int[3];
    int [] targetPosition=new int[3];
    int M1_signed_hysteresis=signed_hysteresis; // split them?
    int M2_signed_hysteresis=signed_hysteresis;
    int M3_signed_hysteresis=signed_hysteresis;
    double d1m1=0.0;
    double d2m1=0.0;
    double d1m2=0.0;
    double d2m2=0.0;
    double dm1=0.0;
    double dm2=0.0;
    boolean go=true;
    boolean abort=false;
    int     numberInside=0; // numer of times new center is inside previous borders (for nuw - second time ends initial iterations)
    int stepNum=0;
    String errorMsg=new String();
    double[] deltas=null;
    double [] displacements;
    int iterNum=0;
    double [][] fourCorners=new double [4][];
    int i;
    int fastMinStep = 3*minStep; 
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error (in thirdByFirstTwo)","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return false;
    }
    startPosition=readElphel10364Motors(motorsIP);
    for (i=0;i<3;i++) targetPosition[i]=startPosition[i];
// first stage - trying to get square that includes the zero-zero point  
    while (go && !abort && (stepNum<numIterations)) {

      if (targetPosition[motor_numbers[0]]>(motorsRange[motor_numbers[0]][1]-initialStep)) {
         errorMsg+="Motor "+(motor_numbers[0]+1)+" exceeds positive limit (arm too close): "+ targetPosition[motor_numbers[0]] +">"+(motorsRange[motor_numbers[0]][1]-initialStep)+"\n";
         abort=true;
      }
      if (targetPosition[motor_numbers[1]]>(motorsRange[motor_numbers[1]][1]-initialStep)) {
         errorMsg+="Motor "+(motor_numbers[1]+1)+" exceeds positive limit (arm too close): "+ targetPosition[motor_numbers[1]] +">"+(motorsRange[motor_numbers[1]][1]-initialStep)+"\n";
         abort=true;
      }
      if (targetPosition[motor_numbers[0]]<(motorsRange[motor_numbers[0]][0]+initialStep)) {
         errorMsg+="Motor "+(motor_numbers[0]+1)+" exceeds negative limit (arm too far): "  + targetPosition[motor_numbers[0]] +"<"+(motorsRange[motor_numbers[0]][0]+initialStep)+"\n";
         abort=true;
      }
      if (targetPosition[motor_numbers[1]]<(motorsRange[motor_numbers[1]][0]+initialStep)) {
         errorMsg+="Motor "+(motor_numbers[1]+1)+" exceeds negative limit (arm too far): "  +targetPosition[motor_numbers[1]]  +"<"+(motorsRange[motor_numbers[1]][0]+initialStep)+"\n";
         abort=true;
      }
      if (abort) break;
      targetPosition[motor_numbers[0]]-=initialStep;
      if (M1_signed_hysteresis>0) { // Make anti-hysteresis move
        targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
        moveElphel10364Motors(motorsIP, targetPosition, 0, true);
        targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
      }
      moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+10, true, "M1:- M2:0");
      fourCorners[0]=measureDeltasCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, average, minUpdateStep, variableStepRatio, fastSteps);
      lastPos=readElphel10364Motors(motorsIP);
      targetPosition[motor_numbers[2]]=lastPos[motor_numbers[2]];
      if (fourCorners[0]==null) abort=true;
      if (!abort) {
        targetPosition[motor_numbers[0]]+=2*initialStep;
        if (M1_signed_hysteresis<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:+ M2:0");
//        fourCorners[1]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        fourCorners[1]=measureDeltasCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, average, minUpdateStep, variableStepRatio, fastSteps);
        lastPos=readElphel10364Motors(motorsIP);
        targetPosition[motor_numbers[2]]=lastPos[motor_numbers[2]];
        if (fourCorners[1]==null) abort=true;
      }
      if (!abort) {
        targetPosition[motor_numbers[0]]-=initialStep;
        if (M1_signed_hysteresis>0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
        }
        targetPosition[motor_numbers[1]]-=initialStep;
        if (M2_signed_hysteresis>0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:0 M2:-");
//        fourCorners[2]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        fourCorners[2]=measureDeltasCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, average, minUpdateStep, variableStepRatio, fastSteps);
        lastPos=readElphel10364Motors(motorsIP);
        targetPosition[motor_numbers[2]]=lastPos[motor_numbers[2]];
        if (fourCorners[2]==null) abort=true;
      }
      if (!abort) {
        targetPosition[motor_numbers[1]]+=2*initialStep;
        if (M2_signed_hysteresis<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:0 M2:+");
//        fourCorners[3]=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        fourCorners[3]=measureDeltasCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], fastMinStep, average, minUpdateStep, variableStepRatio, fastSteps);
        lastPos=readElphel10364Motors(motorsIP);
        targetPosition[motor_numbers[2]]=lastPos[motor_numbers[2]];
        if (fourCorners[3]==null) abort=true;
      }
// return to the center
      if (!abort) {
        targetPosition[motor_numbers[1]]-=initialStep;
        if (M2_signed_hysteresis>0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
        }
        moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1:0 M2:0");
//        deltas=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
        deltas=measureDeltasCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, average, minUpdateStep, variableStepRatio, maxSteps);
        lastPos=readElphel10364Motors(motorsIP);
        targetPosition[motor_numbers[2]]=lastPos[motor_numbers[2]];
      }
      if (abort) {
         errorMsg+="Motor "+(motor_numbers[2]+1)+" probably hit one of the limits trying to compensate focus position for one of the targets\n";
         break;
      }
      if (DEBUG_LEVEL>1) IJ.showMessage("Debug","{{"+fourCorners[0][0]+"},{"+fourCorners[0][1]+"}},\n"+
                              "{"+fourCorners[1][0]+"},{"+fourCorners[1][1]+"}},\n"+
                              "{"+fourCorners[2][0]+"},{"+fourCorners[2][1]+"}},\n"+
                              "{"+fourCorners[3][0]+"},{"+fourCorners[3][1]+"}}");


/// Find double zero
//    double d1m1,d2m1,d1m2,d2m2;
      d1m1=(1.0*(fourCorners[1][0]-fourCorners[0][0]))/(2.0*initialStep);
      d2m1=(1.0*(fourCorners[1][1]-fourCorners[0][1]))/(2.0*initialStep);
      d1m2=(1.0*(fourCorners[3][0]-fourCorners[2][0]))/(2.0*initialStep);
      d2m2=(1.0*(fourCorners[3][1]-fourCorners[2][1]))/(2.0*initialStep);

      dm1 =(d1m2*deltas[1]-d2m2*deltas[0])/(d2m2*d1m1-d1m2*d2m1);
      dm2 =(d1m1*deltas[1]-d2m1*deltas[0])/(d2m1*d1m2-d1m1*d2m2);
      if ((Math.abs(dm1)>initialStep) || (Math.abs(dm1)>initialStep)) { // desired location outside of the square
        if      (dm1 >  initialStep) dm1= initialStep;
        else if (dm1 < -initialStep) dm1= -initialStep;
        if      (dm2 >  initialStep) dm2= initialStep;
        else if (dm2 < -initialStep) dm2= -initialStep;
      } else {
        numberInside++;
        if (numberInside > 1) go = false; // Stage 1 done
      }
// Make step
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+deltas[0]+"\t"+deltas[1]+"\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tStage 1\n");

      targetPosition[motor_numbers[0]]+=dm1;
      if ((M1_signed_hysteresis*dm1)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
      }
      targetPosition[motor_numbers[1]]+=dm2;
      if ((M2_signed_hysteresis*dm2)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
      }
      moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1+"+dm1+" M2+"+dm2+"STEP:"+stepNum);
      stepNum++;
//      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+delta[0]+"\t"+delta[1]+"\t"+"\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tmade step\n");
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t-\t-\t-\t-\t-\t-\t-\t-\t-\tmade step\n");
    }
/// Stage 2 - uae previously found d1m1,d1m2,d2m1,d2m2 to find target relative to current position, exit at certain step
    for (iterNum=0;!abort && (iterNum<numIterations);iterNum++) {

//      deltas=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps);
      deltas=measureDeltasCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, average, minUpdateStep, variableStepRatio, maxSteps);
      lastPos=readElphel10364Motors(motorsIP);
      targetPosition[motor_numbers[2]]=lastPos[motor_numbers[2]];
      dm1 =(d1m2*deltas[1]-d2m2*deltas[0])/(d2m2*d1m1-d1m2*d2m1);
      dm2 =(d1m1*deltas[1]-d2m1*deltas[0])/(d2m1*d1m2-d1m1*d2m2);
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+deltas[0]+"\t"+deltas[1]+"\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tStage 2\n");
// Make step
      targetPosition[motor_numbers[0]]+=dm1;
      if ((M1_signed_hysteresis*dm1)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[0]]-=M1_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[0]]+=M1_signed_hysteresis;
      }
      targetPosition[motor_numbers[1]]+=dm2;
      if ((M2_signed_hysteresis*dm2)<0) { // Make anti-hysteresis move
          targetPosition[motor_numbers[1]]-=M2_signed_hysteresis;
          moveElphel10364Motors(motorsIP, targetPosition, 0, true);
          targetPosition[motor_numbers[1]]+=M2_signed_hysteresis;
      }
      moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove+1, true, "M1+"+dm1+" M2+"+dm2+"iterNum:"+iterNum);
//      if (DEBUG_LEVEL>0) sb.append(iterNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t"+d1m1+"\t"+d1m2+"\t"+d2m1+"\t"+d2m2+"\t"+dm1+"\t"+dm2+"\tStage 2\n");
      if (DEBUG_LEVEL>0) sb.append(stepNum+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t-\t-\t-\t-\t-\t-\t-\t-\t-\tmade step\n");

    }
    if (DEBUG_LEVEL>0 ) {

//       m3s=measureM3Deltas(motorsIP,M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, maxStep, targets_dirs, ColorMasks, V2H, refine, refineAverage, refineExpand, maxSteps,false);
       displacements= findThirdMotorCalib(motorsIP, M3_signed_hysteresis, motorsWaitAfterMove, motor_numbers[2], minStep, 2, average, minUpdateStep, variableStepRatio, maxSteps);
       sb.append("Final\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\tt1="+
                 displacements[0]+"\tt2="+displacements[1]+"\tt3="+displacements[2]+"\tt4="+displacements[3]+"\tt5="+displacements[4]+"\t-\t-\t-\t-\tmade step\n");
       new TextWindow("Focus adjustmet", header, sb.toString(), 1200,800);
    }
    if (abort) {
      moveElphel10364Motors(motorsIP, startPosition, 0, true, "abort - return",true); /// restore initial position
      IJ.showMessage("Error",errorMsg);
      return false;
    }
    return true;
  }

// Adjusts motorNumber (should be 2 - last) so the center target is at zero, returns differences
// returns (target2 - target4), (target1 - target5), target3 or NULL if failed
// if maxSteps ==0 - will not move, just measure
  double []  measureDeltasCalib(   String        motorsIP,
                                int           signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                                double        motorsWaitAfterMove,
                                int           motorNumber,
                                int           minStep, // minimal motor rotation - positive
                                int           average, // number of measurements to average
                                int           minUpdateStep, // ipdate thisStepsPerUM only if the last step was > this threshold
                                double        variableStepRatio, // 2.0 - maximal thisStepsPerUM change in one step
                                int           maxSteps) {// - just for debugging - maximal number of steps to try before giving up
    double [] rslt = null;
    double [] displacements= findThirdMotorCalib(
                             motorsIP,
                             signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                             motorsWaitAfterMove,
                             motorNumber,
                             minStep,           // minimal motor rotation - positive
                             2, // center target
                             average,     // number of measurements to average
                             minUpdateStep,              // ipdate thisStepsPerUM only if the last step was > this threshold
                             variableStepRatio,          // 2.0 - maximal thisStepsPerUM change in one step
                             maxSteps);    // - just for debugging - maximal number of steps to try before giving up

    if (displacements!=null) {
      rslt = new double[3];
      rslt[0]=displacements[1]-displacements[3];
      rslt[1]=displacements[0]-displacements[4];
      rslt[2]= displacements[2];
    }
    return rslt;

}





// version that uses motors calibration data
  double [] findThirdMotorCalib( // Returns null if reached motor limits, otherwise focus results for all colors/targets
                                String        motorsIP,
                                int           signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                                double        motorsWaitAfterMove,
                                int           motorNumber,
                                int           minStep, // minimal motor rotation - positive
                                int           targetNumber,
                                int           average, // number of measurements to average
                                int           minUpdateStep, // ipdate thisStepsPerUM only if the last step was > this threshold
                                double        variableStepRatio, // 2.0 - maximal thisStepsPerUM change in one step
                                int           maxSteps) {// - just for debugging - maximal number of steps to try before giving up

// global  motors_sideWays, motors_calibSteps,stepsPerUM, offsetsPerGoal[goalNumber-1],focusOffset,targetLocations, targetScanData, imp_camera

    String header="#\tM1\tM2\tM3\tstepsPerUM\tDisplacement\tT1\tT2\tT3\tT4\tT5";
    StringBuffer sb = new StringBuffer();
    int i;
    int [] currentPosition;
    int [] originalPosition=new int[3];
    int [] lastPos=new int[3];
    int [] targetPosition=new int[3];
    double [] displacements=null;
    double oldDisplacement,newDisplacement;
    double[][][] resolutions=null;
    int numSteps=0;
    int threshold;
    double thisStepsPerUM;
    double newStepsPerUM;
    boolean go=true;
    boolean aborted=false;
    double  step;
    if (offsetsPerGoal==null) {
      IJ.showMessage("Error","Motors are not calibrated yet.\nPlease run \"Scan Focus\" first");
      return null;
    }
    thisStepsPerUM = 0.75*stepsPerUM/motors2center[motorNumber]; // initial value of steps/micron of the current motor. This value will be updated only if the last step was big enough
                                                                 // Use first step smaller to reduce overshoots
    currentPosition=readElphel10364Motors(motorsIP);
    for (i=0; i<3;i++) {
      targetPosition[i]= currentPosition[i];
      originalPosition[i]= currentPosition[i];
    }
    targetPosition[motorNumber]-=signed_hysteresis;
    if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\n");
    moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
    targetPosition[motorNumber]+=signed_hysteresis; // restore back
    lastPos=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove, true, "start position");
    if (DEBUG_LEVEL>1) sb.append("RESTORE"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\n");

    resolutions=     measureTargets (targetLocations, imp_camera, average);
    displacements=findDisplacemntsFromResolutions (targetScanData, resolutions, motors_sideWays, motors_calibSteps,stepsPerUM, offsetsPerGoal[goalNumber-1],focusOffset);
    oldDisplacement=displacements[targetNumber]; // in microns

    if (DEBUG_LEVEL>1) sb.append(numSteps+ "\t"+currentPosition[0]+ "\t"+currentPosition[1]+ "\t"+currentPosition[2]+  "\t"+IJ.d2s(thisStepsPerUM,1)+ "\t"+IJ.d2s(oldDisplacement,2)+
                       "\t"+IJ.d2s(displacements[0],1)+ "\t"+IJ.d2s(displacements[1],1)+ "\t"+IJ.d2s(displacements[2],1)+ "\t"+IJ.d2s(displacements[3],1)+ "\t"+IJ.d2s(displacements[4],1)+ "\n");
    while (go && !aborted && (numSteps<=maxSteps)) {
/// Calculate required move:
      step = -oldDisplacement * thisStepsPerUM; //if displacement was negative, go positive
//      if (Math.abs(step) < (minStep/2)) {
      if (Math.abs(step) < (minStep)) {
        go=false;
        break;
      }

      targetPosition[motorNumber]=currentPosition[motorNumber]+(int) (step+0.5);
      if (((targetPosition[motorNumber]-currentPosition[motorNumber])*signed_hysteresis)<0) {
        targetPosition[motorNumber]-=signed_hysteresis;
        if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\n");
        moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
        targetPosition[motorNumber]+=signed_hysteresis; // restore back
      }

/// Move to the new position
      lastPos=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove, true, numSteps+":");
      threshold= (int) Math.abs(step/2)+2;
      if (Math.abs(lastPos[motorNumber]-targetPosition[motorNumber])>threshold) { /// did not move
        lastPos=readElphel10364Motors(motorsIP, 1); // wait one second more
        if (Math.abs(lastPos[motorNumber]-targetPosition[motorNumber])>threshold) { /// did not move
          IJ.showMessage("Motor ERROR","Did not move: targetPosition[motorNumber]="+targetPosition[motorNumber]+" lastPos[motorNumber]="+lastPos[motorNumber]+" threshold="+threshold);
          aborted=true;
          break;
        }
      }

// Measure at new position
      resolutions=     measureTargets (targetLocations, imp_camera, average);
      displacements=   findDisplacemntsFromResolutions (targetScanData, resolutions, motors_sideWays, motors_calibSteps,stepsPerUM, offsetsPerGoal[goalNumber-1],focusOffset);
      newDisplacement= displacements[targetNumber]; // in microns
/// Update thisStepsPerUM?
      if ((Math.abs(step) > minUpdateStep) &&  ((newDisplacement-oldDisplacement)*step*thisStepsPerUM>0)){ // step big enough and same sign as expected
        newStepsPerUM= step/(newDisplacement-oldDisplacement);
        if (newStepsPerUM/thisStepsPerUM > variableStepRatio) newStepsPerUM= variableStepRatio* thisStepsPerUM;
        else if (newStepsPerUM/thisStepsPerUM < (1.0/ variableStepRatio)) newStepsPerUM= thisStepsPerUM/variableStepRatio;
        thisStepsPerUM=Math.sqrt(thisStepsPerUM*newStepsPerUM);
      }

      currentPosition[motorNumber]=targetPosition[motorNumber];
      oldDisplacement=newDisplacement;
      if (DEBUG_LEVEL>1) sb.append(numSteps+ "\t"+currentPosition[0]+ "\t"+currentPosition[1]+ "\t"+currentPosition[2]+  "\t"+IJ.d2s(thisStepsPerUM,1)+ "\t"+IJ.d2s(oldDisplacement,2)+
                       "\t"+IJ.d2s(displacements[0],1)+ "\t"+IJ.d2s(displacements[1],1)+ "\t"+IJ.d2s(displacements[2],1)+ "\t"+IJ.d2s(displacements[3],1)+ "\t"+IJ.d2s(displacements[4],1)+ "\n");
    }
    if (aborted) {
      moveElphel10364Motors(motorsIP, originalPosition, 0, true, "ABORTED - GOING BACK", true);
      if (DEBUG_LEVEL>1) sb.append(numSteps+ "\t"+originalPosition[0]+ "\t"+originalPosition[1]+ "\t"+originalPosition[2]+  "\t---\t---\t---\t---\t---\t---\t---\n");
    }
    if (DEBUG_LEVEL>1 ) {
       new TextWindow("motor iterations - TARGET="+(targetNumber+1)+" MOTOR="+(motorNumber+1), header, sb.toString(), 1200,800);
    }
    return displacements;
  }


// this (old) version does not use calibrtation data
  double [][][] findThirdMotor( // Returns null if reached motor limits, otherwise focus results for all colors/targets
                                ImagePlus      imp_src,
                                String        motorsIP,
                                int           signed_hysteresis, // if positive move positive directly, negative - with hysteresis. If negatgive - move negative directly, positive with hysteresis.
                                double        motorsWaitAfterMove,
                                int           motorNumber,
                                int           minStep, // minimal motor rotation - positive
                                int           maxStep, // maximal motor rotation (positive)
                                int           targetNumber,
                                int           colorMask, // bitmask of colors to average: 2 - red, 4 - blue, 9 - Gr+Gb , 0 - use (G-R) instead - for central target
                                double        V2H, // target ratio of vertial/horizontal resolutions for positive decrDirection, H/V - for negative. Gree/Red for colorMask==0
//                                double [][][] oldResult, // oldResult - result from the previous run (or null)
                                int           refine, /// number of minimal steps made to refine the final position (will be half each way
                                int           refineAverage,
                                int           refineExpand,
                                int           maxSteps) {// - just for debugging - maximal number of steps to try before giving up
    String header="#\tM1\tM2\tM3\tDiff\tV_Gr\tH_Gr\tV_R\tH_R\tV_B\tH_B\tV_Gb\tH_Gb";
    StringBuffer sb = new StringBuffer();

    ImagePlus imp=imp_src;
    double oldDiff,newDiff;
    double anotherDiff=0.0;
    int i,n;
    int [] currentPosition;
    int [] lastPos=new int[3];
    int [] targetPosition=new int[3];
    int [] anotherPosition=new int[3]; // used in divide
    double [][][] newResult=null;
    double [][][] oldResult=null;
    int numSteps=0;
    int threshold;
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return null;
    }
    if (oldResult==null) {
      imp=jp4_instance.openURL(imp);
      if (imp==null) {
        IJ.showMessage("Camera Image ERROR","Failed to open camera image");
        return null;
      }
    imp_camera=imp;
      oldResult = measureTargets (targetLocations, imp);
    }
    oldDiff= calcGoalDiff(oldResult[targetNumber], V2H, colorMask);

 //OK2moveElphel10364Motors(positions))  
 // move, starting with increasing (limited by maxStep) steps until error changes sign, then divide the remaining interval by 2 until step is lower than minStep. Use signed_hysteresis
//  int[] readElphel10364Motors(String ip)
    currentPosition=readElphel10364Motors(motorsIP);
    if (DEBUG_LEVEL>1) sb.append(numSteps+
           "\t"+currentPosition[0]+ "\t"+currentPosition[1]+ "\t"+currentPosition[2]+ 
           "\t"+IJ.d2s(  oldDiff,4)+
           "\t"+IJ.d2s(oldResult[targetNumber][0][0],4)+"\t"+IJ.d2s(oldResult[targetNumber][1][0],4)+
           "\t"+IJ.d2s(oldResult[targetNumber][0][1],4)+"\t"+IJ.d2s(oldResult[targetNumber][1][1],4)+
           "\t"+IJ.d2s(oldResult[targetNumber][0][2],4)+"\t"+IJ.d2s(oldResult[targetNumber][1][2],4)+
           "\t"+IJ.d2s(oldResult[targetNumber][0][3],4)+"\t"+IJ.d2s(oldResult[targetNumber][1][3],4)+
           "\n");


    for (i=0; i<3;i++) {
      targetPosition[i]= currentPosition[i];
      anotherPosition[i]=currentPosition[i];
    }
    boolean go=true;
    boolean aborted=false;
//    int dir = decrDirection;
    int step=minStep;
    while (go && !aborted) {
      
      if (oldDiff>0) { /// move in specified direction
//        targetPosition[motorNumber]=currentPosition[motorNumber]+dir*step;
        targetPosition[motorNumber]=currentPosition[motorNumber]+step;
      } else { /// move in opposite to specified direction
//        targetPosition[motorNumber]=currentPosition[motorNumber]-dir*step;
        targetPosition[motorNumber]=currentPosition[motorNumber]-step;
      }
      if (!OK2moveElphel10364Motors(targetPosition)) { // motor out of limits
        aborted=true;
        break;
      }
      if (((targetPosition[motorNumber]-currentPosition[motorNumber])*signed_hysteresis)<0) {
        targetPosition[motorNumber]-=signed_hysteresis;
        if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\t---\t---\n");
        moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
        targetPosition[motorNumber]+=signed_hysteresis; // restore back
      }
      threshold= step>>1; if (threshold <10) threshold=10;
      lastPos=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove, true, "stage1, "+numSteps);
      if (Math.abs(lastPos[motorNumber]-targetPosition[motorNumber])>threshold) { /// did not move
        lastPos=readElphel10364Motors(motorsIP, 1); // wait one second more
        if (Math.abs(lastPos[motorNumber]-targetPosition[motorNumber])>threshold) { /// did not move
          IJ.showMessage("Motor ERROR Stage 1","Did not move: targetPosition[motorNumber]="+targetPosition[motorNumber]+" lastPos[motorNumber]="+lastPos[motorNumber]+" (step>>1)="+(step>>1));
          aborted=true;
          break;
        }
      }
      imp=jp4_instance.openURL(imp);
      if (imp==null) {
        IJ.showMessage("Camera Image ERROR","Failed to open camera image");
        aborted=true;
        break;
      }
      imp_camera=imp;

      newResult = measureTargets (targetLocations, imp);
      newDiff= calcGoalDiff(newResult[targetNumber], V2H, colorMask);

      if (DEBUG_LEVEL>1) sb.append(numSteps+
           "\t"+lastPos[0]+ "\t"+lastPos[1]+ "\t"+lastPos[2]+ 
           "\t"+IJ.d2s(  newDiff,4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][0],4)+"\t"+IJ.d2s(newResult[targetNumber][1][0],4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][1],4)+"\t"+IJ.d2s(newResult[targetNumber][1][1],4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][2],4)+"\t"+IJ.d2s(newResult[targetNumber][1][2],4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][3],4)+"\t"+IJ.d2s(newResult[targetNumber][1][3],4)+
           "\n");
      numSteps+=1;
      if (newDiff*oldDiff >0 ) { // same sign of error
         step*=2;
         if (step > maxStep) step=maxStep;
         currentPosition[motorNumber]=targetPosition[motorNumber];
         oldDiff= newDiff;
      } else {
        anotherDiff=oldDiff;
        oldDiff= newDiff;
        anotherPosition[motorNumber]=currentPosition[motorNumber]; // divide between currentPosition and anotherPosition
        currentPosition[motorNumber]= targetPosition[motorNumber];
        go=false;
        break; // go to second part of the algorithm - diminishing steps
      }
      if ((maxSteps>0) && (numSteps>maxSteps)) {
        IJ.showMessage("Error","Number of steps exceeded limit ("+maxSteps+") during stage 1, aborting ");
        aborted=true;
        break;
      }

    }
      if (DEBUG_LEVEL>1) sb.append("STAGE2"+
           "\t"+lastPos[0]+ "\t"+lastPos[1]+ "\t"+lastPos[2]+ 
           "\tanother:"+
           "\t"+anotherPosition[0]+ "\t"+anotherPosition[1]+ "\t"+anotherPosition[2]+ 
           "\tcurrent:"+
           "\t"+currentPosition[0]+ "\t"+currentPosition[1]+ "\t"+currentPosition[2]+ 
           "\t---\n");


    if (aborted) {
       if (DEBUG_LEVEL>1 ) {
         new TextWindow("motor iterations(aborted) - TARGET="+(targetNumber+1)+" COLORS="+colorMask+" V2H="+V2H, header, sb.toString(), 1200,800);
       }
       return null;
    }
    while ((step>minStep) && !aborted) {
/// Divide by 2 method
//      step=Math.abs(anotherPosition[motorNumber]-currentPosition[motorNumber])/2;
//      targetPosition[motorNumber]=(anotherPosition[motorNumber]+currentPosition[motorNumber])/2;
/// Divide by interpolation
      targetPosition[motorNumber]=currentPosition[motorNumber]+
                                  ((int) (((anotherPosition[motorNumber]-currentPosition[motorNumber])*oldDiff)/(oldDiff-anotherDiff)));
      step=Math.abs(targetPosition[motorNumber]-currentPosition[motorNumber]);
      if (((targetPosition[motorNumber]-currentPosition[motorNumber])*signed_hysteresis)<0) {
        targetPosition[motorNumber]-=signed_hysteresis;
        if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\t---\t---\n");
        moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
        targetPosition[motorNumber]+=signed_hysteresis; // restore back
      }
      threshold= step>>1; if (threshold <10) threshold=10;
      lastPos=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove, true, "stage2, "+numSteps);
      if (Math.abs(lastPos[motorNumber]-targetPosition[motorNumber])>threshold) { /// did not move
        lastPos=readElphel10364Motors(motorsIP, 1); // wait one second more
        if (Math.abs(lastPos[motorNumber]-targetPosition[motorNumber])>threshold) { /// did not move
          IJ.showMessage("Motor ERROR Stage 2","Did not move: targetPosition[motorNumber]="+targetPosition[motorNumber]+" lastPos[motorNumber]="+lastPos[motorNumber]+" (step>>1)="+(step>>1));
          aborted=true;
          break;
        }
      }
      imp=jp4_instance.openURL(imp);
      if (imp==null) {
        IJ.showMessage("Camera Image ERROR","Failed to open camera image");
        aborted=true;
        break;
      }
    imp_camera=imp;
      newResult = measureTargets (targetLocations, imp);
      newDiff= calcGoalDiff(newResult[targetNumber], V2H, colorMask);

      if (DEBUG_LEVEL>1) sb.append(numSteps+
           "\t"+lastPos[0]+ "\t"+lastPos[1]+ "\t"+lastPos[2]+ 
           "\t"+IJ.d2s(  newDiff,4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][0],4)+"\t"+IJ.d2s(newResult[targetNumber][1][0],4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][1],4)+"\t"+IJ.d2s(newResult[targetNumber][1][1],4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][2],4)+"\t"+IJ.d2s(newResult[targetNumber][1][2],4)+
           "\t"+IJ.d2s(newResult[targetNumber][0][3],4)+"\t"+IJ.d2s(newResult[targetNumber][1][3],4)+
           "\n");
      if (newDiff*oldDiff <0 ) { // crossed zero, reverse
        anotherPosition[motorNumber]=currentPosition[motorNumber];
        anotherDiff=oldDiff;
      }
      oldDiff= newDiff;
      currentPosition[motorNumber]= targetPosition[motorNumber];
      numSteps+=1;
      if ((maxSteps>0) && (numSteps>maxSteps)) {
        IJ.showMessage("Error","Number of steps exceeded limit ("+maxSteps+") during stage 2, aborting ");
        aborted=true;
        break;
      }
    }
    if (aborted) {
       if (DEBUG_LEVEL>1 ) {
          new TextWindow("motor iterations - TARGET="+(targetNumber+1)+" COLORS="+colorMask+" V2H="+V2H, header, sb.toString(), 1200,800);
       }
       return null; // only if motor got stuck or no jp4
    }
//  int           refineAverage,
//  int           refineExpand,

    if (refine>0) { /// make final refinement, no
// Refine will run up to two passes. First - as specified, making sure correction is not more than 75% of the refine range.
// If that fails - iuncrease refine range and retry
// If it fails again - use original (before refine) value.
      double refineThreshold=0.75;
      int refinePass;
      double x,SX,SX2,SY,SXY;
      int N;
      double a,b;
      int nmeas;
      go=true;
      aborted=false;
      for (refinePass=0; (refinePass<2) && go && !aborted;refinePass++) {
        if (refinePass>0) refine*=refineExpand;
        if (DEBUG_LEVEL>1) sb.append("Refining"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\t---\t---\n");
        anotherPosition[motorNumber]= targetPosition[motorNumber]; // save original position
        SX=0.0;SX2=0.0;SY=0.0;SXY=0.0;
        for (n=0;n<=refine;n++) {
          numSteps++; // no controll for exceedidng the maximum, here number of steps is defined by a parameter
          targetPosition[motorNumber]=anotherPosition[motorNumber]+((signed_hysteresis>=0)?1:-1)*minStep*(2*n-refine)/2;
          if (n==0) {
            targetPosition[motorNumber]-=signed_hysteresis;
            if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\t---\t---\n");
            moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
            targetPosition[motorNumber]+=signed_hysteresis; // restore back
          }
          lastPos=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove, true, "refining"+refinePass+", "+numSteps);

/// Make several resolution measurements and aver4age results
          newDiff=0.0;
          for (nmeas=0; nmeas <refineAverage; nmeas++) {
            imp=jp4_instance.openURL(imp);
            if (imp==null) {
              IJ.showMessage("Camera Image ERROR","Failed to open camera image");
              aborted=true;
              break;
            }
            imp_camera=imp;
            newResult= measureTargets (targetLocations, imp);
            newDiff+= calcGoalDiff(newResult[targetNumber], V2H, colorMask);
          }
          newDiff/=refineAverage;
           if (DEBUG_LEVEL>1) sb.append(numSteps+
             "\t"+lastPos[0]+ "\t"+lastPos[1]+ "\t"+lastPos[2]+ 
             "\t"+IJ.d2s(  newDiff,4)+
             "\t"+IJ.d2s(newResult[targetNumber][0][0],4)+"\t"+IJ.d2s(newResult[targetNumber][1][0],4)+
             "\t"+IJ.d2s(newResult[targetNumber][0][1],4)+"\t"+IJ.d2s(newResult[targetNumber][1][1],4)+
             "\t"+IJ.d2s(newResult[targetNumber][0][2],4)+"\t"+IJ.d2s(newResult[targetNumber][1][2],4)+
             "\t"+IJ.d2s(newResult[targetNumber][0][3],4)+"\t"+IJ.d2s(newResult[targetNumber][1][3],4)+
             "\n");
          x=targetPosition[motorNumber]-anotherPosition[motorNumber];
          SX+=x;
          SX2+=x*x;
          SY+=newDiff;
          SXY+=x*newDiff;
        }
/// Now interpolate diffs
        N=refine+1;
        a=(N*SXY-SY*SX)/(N*SX2-SX*SX);
        b=(SX*SXY-SY*SX2)/(SX*SX-SX2*N);
        x=-b/a;
/// Did that work? Or correction was too hight?
        if (Math.abs(x) < refineThreshold* minStep* refine /2 ) { // OK
          targetPosition[motorNumber]=anotherPosition[motorNumber]+(int) x;
          go=false;
        } else if (refinePass<1) {
            break; //retry
        }else { // failed, last pass
           go=false; // (not needed)
           targetPosition[motorNumber]=anotherPosition[motorNumber];
           if (DEBUG_LEVEL>1) sb.append("Refined FAILED"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+
             "\tN="+  N+
             "\tSX="+ IJ.d2s(SX, 4)+
             "\tSX2="+IJ.d2s(SX2,4)+
             "\tSY="+ IJ.d2s(SY, 4)+
             "\tSXY="+IJ.d2s(SXY,4)+
             "\ta="+  IJ.d2s(a,4)+
             "\tb="+  IJ.d2s(b,4)+
             "\tx="+  IJ.d2s(x,4)+
             "\t---\n");
        }
// Adding hysteresis
        targetPosition[motorNumber]-=signed_hysteresis;
        if (DEBUG_LEVEL>1) sb.append("HYSTERESIS"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+ "\t---\t---\t---\t---\t---\t---\t---\t---\t---\n");
        moveElphel10364Motors(motorsIP, targetPosition, 0, true); /// no sleep here
        targetPosition[motorNumber]+=signed_hysteresis;
        lastPos=moveElphel10364Motors(motorsIP, targetPosition, motorsWaitAfterMove, true, "final: "+numSteps);
        if (DEBUG_LEVEL>1) sb.append("Refined to"+"\t"+targetPosition[0]+ "\t"+targetPosition[1]+ "\t"+targetPosition[2]+
             "\tN="+  N+
             "\tSX="+ IJ.d2s(SX, 4)+
             "\tSX2="+IJ.d2s(SX2,4)+
             "\tSY="+ IJ.d2s(SY, 4)+
             "\tSXY="+IJ.d2s(SXY,4)+
             "\ta="+  IJ.d2s(a,4)+
             "\tb="+  IJ.d2s(b,4)+
             "\tx="+  IJ.d2s(x,4)+
             "\t---\n");

      }
    }
    if (DEBUG_LEVEL>1 ) {
      new TextWindow("motor iterations - TARGET="+(targetNumber+1)+" COLORS="+colorMask+" V2H="+V2H, header, sb.toString(), 1200,800);
    }
    return newResult;
  }







//  double [][][][][] scanMotors(ImagePlus imp_src, String motorsIP, int motors_hysteresis, double motorsWaitAfterMove,int [] startScan, int [] endScan, int motors_scan_steps,  double[][][] targets) {
  double [][][][][] scanMotors(ImagePlus imp_src, String motorsIP, int motors_hysteresis, double motorsWaitAfterMove,int [] startScan, int [] endScan, int motors_scan_steps, int average, boolean bothDirections) {
    int i,n;
    ImagePlus imp=imp_src;
    int [] originalPositions=readElphel10364Motors(motorsIP);
    int [] targetPositions = new int [3];
    if (targetLocations==null) { /// Global array for targets
      IJ.showMessage("Error","Target locations are undefined.\nPlease run \"Find Targets\" first");
      return null;
    }
    double [][][][][] scanResults= new double [motors_scan_steps+1][ bothDirections?2:1][targetLocations.length][2][4]; /// [points]  [direction (for hysteresis) 0 increase, 1 - decrease] [targets] [vert=0/hor=1] [colors - Gr-R-B-Gb]
//    double ar;
/// move to startScan with hysteresis
    for (i=0; i<3;i++)  targetPositions[i]=startScan[i]+((startScan[i]>endScan[i])?1:-1)*motors_scan_steps;
    moveElphel10364Motors(motorsIP, targetPositions, 0, true); /// no sleep here
/// scan forward
/*
     if (jp4_instance==null)  jp4_instance= new JP46_Reader_camera();
     imp=jp4_instance.openURL(imp);
     if (imp==null) {
        IJ.showMessage("Camera Image ERROR","Failed to open camera image");
        moveElphel10364Motors(motorsIP, originalPositions, 0, true, "abort - return",true); /// restore initial position

        return null;
     }
*/
    for (n=0;n<=motors_scan_steps;n++) {
      for (i=0; i<3;i++)  {
        targetPositions[i]=startScan[i];
        if (startScan[i]!=endScan[i]) targetPositions[i]+=((endScan[i]-startScan[i])*n)/motors_scan_steps;
      }
      moveElphel10364Motors(motorsIP, targetPositions, motorsWaitAfterMove, true);
// now acquire and measure image
//      imp=jp4_instance.openURL(imp_src);
/*
      imp=jp4_instance.openURL(imp);
      if (imp==null) {
        IJ.showMessage("Camera Image ERROR","Failed to open camera image");
        moveElphel10364Motors(motorsIP, originalPositions, 0, true, "abort - return",true); /// restore initial position

        return null;
      }
*/
///Make sure targets are defined (image existence is already tested)
      scanResults [n][0] = measureTargets (targetLocations, imp, average);
    }

/// move to endScan with hysteresis
    if (bothDirections) {
      for (i=0; i<3;i++)  targetPositions[i]=endScan[i]+((endScan[i]>startScan[i])?1:-1)*motors_scan_steps;
      moveElphel10364Motors(motorsIP, targetPositions, 0, true); /// no sleep here
/// scan backwards
      for (n=motors_scan_steps;n>=0;n--) {
        for (i=0; i<3;i++)  {
          targetPositions[i]=startScan[i];
          if (startScan[i]!=endScan[i]) targetPositions[i]+=((endScan[i]-startScan[i])*n)/motors_scan_steps;
        }
        moveElphel10364Motors(motorsIP, targetPositions, motorsWaitAfterMove, true);
// now acquire and measure image
//      imp=jp4_instance.openURL(imp_src);
/*
        imp=jp4_instance.openURL(imp);
        if (imp==null) {
          IJ.showMessage("Camera Image ERROR","Failed to open camera image");
          moveElphel10364Motors(motorsIP, originalPositions, 0, true, "abort - return",true); /// restore initial position
          return null;
        }
*/
///Make sure targets are defined (image existence is already tested)
        scanResults [n][1] = measureTargets (targetLocations, imp, average);
      }
    }
    moveElphel10364Motors(motorsIP, originalPositions, 0, true, "return",true); /// restore initial position
    return scanResults;
  }


  int[] moveElphel10364Motors(String ip, int [] positions, double sleep, boolean showStatus, String message, boolean hysteresis)  {
// global motor3rd_signed_hysteresis
      int i;
      int [] curpos=readElphel10364Motors(ip);
      if (hysteresis) {
        for (i=0;i<3;i++) if ((positions[i]-curpos[i])*motor3rd_signed_hysteresis<0) curpos[i]-=motor3rd_signed_hysteresis;
        moveElphel10364Motors(ip, positions, 0, showStatus, "comp. hyster. "+message);
      }
      return moveElphel10364Motors(ip, positions, sleep, showStatus, message);
  }

  int[] moveElphel10364Motors(String ip, int [] positions, double sleep, boolean showStatus, boolean hysteresis)  {
// global motor3rd_signed_hysteresis
      int i;
      int [] curpos=readElphel10364Motors(ip);
      if (hysteresis) {
        for (i=0;i<3;i++) if ((positions[i]-curpos[i])*motor3rd_signed_hysteresis<0) curpos[i]-=motor3rd_signed_hysteresis;
        moveElphel10364Motors(ip, positions, 0, showStatus,"comp. hyster. ");
      }
      return moveElphel10364Motors(ip, positions, sleep, showStatus, "");
  }



  int[] moveElphel10364Motors(String ip, int [] positions, double sleep, boolean showStatus)  {
      return moveElphel10364Motors(ip, positions, sleep, showStatus,"");
  }
  int[] moveElphel10364Motors(String ip, int [] positions, double sleep, boolean showStatus, String message)  {
      IJ.showStatus(message+", Moving m1="+positions[0]+" m2="+positions[1]+" m3="+positions[2]);
      return moveElphel10364Motors(ip, positions, sleep);
  }

  int[] moveElphel10364Motors(String ip, int [] positions, double sleep)  {
      moveElphel10364Motors(ip, positions);
      return readElphel10364Motors( ip, sleep);
  }

  int[] moveElphel10364Motors(String ip, int [] positions, boolean showStatus)  {
      return moveElphel10364Motors( ip, positions, showStatus, "");
  }

  int[] moveElphel10364Motors(String ip, int [] positions, boolean showStatus, String message)  {
      IJ.showStatus(message+", Moving m1="+positions[0]+" m2="+positions[1]+" m3="+positions[2]);
      return moveElphel10364Motors( ip, positions);
  }



  int[] moveElphel10364Motors(String ip, int [] positions)  {
// public double [][]    motorsRange{{-5000,5000},{-5000,5000},{-5000,5000}};
      if (OK2moveElphel10364Motors(positions)) return commandElphel10364Motors("http://"+ip+"/10364.php?m1="+positions[0]+"&m2="+positions[1]+"&m3="+positions[2]);
      return readElphel10364Motors(ip); // just read
  }

  boolean OK2moveElphel10364Motors( int [] positions)  {
// public double [][]    motorsRange{{-5000,5000},{-5000,5000},{-5000,5000}};
       return (positions[0]>=motorsRange[0][0]) && (positions[0]<= motorsRange[0][1]) &&
              (positions[1]>=motorsRange[1][0]) && (positions[1]<= motorsRange[1][1]) &&
              (positions[2]>=motorsRange[2][0]) && (positions[2]<= motorsRange[2][1]);
  }


  int[] readElphel10364Motors(String ip, double sleep)  {
      return commandElphel10364Motors("http://"+ip+"/10364.php?sleep="+sleep);
  }

  int[] readElphel10364Motors(String ip)  {
      return commandElphel10364Motors("http://"+ip+"/10364.php");
  }
  int[] commandElphel10364Motors(String url)  {
      int [] motorPositions= new int [3];
      Document dom=null;
      try {
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        DocumentBuilder db = dbf.newDocumentBuilder();
//        dom = db.parse(url+"/10364.php");
        dom = db.parse(url);
///        dom.getDocumentElement().normalize();
//        IJ.showMessage("Root element :" + dom.getDocumentElement().getNodeName()); 
        if (!dom.getDocumentElement().getNodeName().equals("motors")) {
           IJ.showMessage("Root element: expected 'motors', got'" + dom.getDocumentElement().getNodeName()+"'"); 
           return null;
        }
//        IJ.showMessage("name="+((Node) dom.getDocumentElement().getElementsByTagName("motor1").item(0)).getNodeName());
//        IJ.showMessage("value="+((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor1").item(0)).getChildNodes().item(0))).getNodeValue());
        motorPositions[0]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor1").item(0)).getChildNodes().item(0))).getNodeValue());
        motorPositions[1]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor2").item(0)).getChildNodes().item(0))).getNodeValue());
        motorPositions[2]=Integer.parseInt(((Node) (((Node) dom.getDocumentElement().getElementsByTagName("motor3").item(0)).getChildNodes().item(0))).getNodeValue());
//        IJ.showMessage("motor1="+motorPositions[0]+"\n"+
//                       "motor2="+motorPositions[1]+"\n"+
//                       "motor3="+motorPositions[2]);
      } catch(MalformedURLException e){
	System.out.println("Please check the URL:" + e.toString() );
        return null;
      } catch(IOException  e1){
         IJ.showStatus("");
         String error = e1.getMessage();
         if (error==null || error.equals(""))  error = ""+e1;
         IJ.showMessage("commandElphel10364Motors ERRROR", ""+error);
         return null;
      }catch(ParserConfigurationException pce) {
         pce.printStackTrace();
         return null;
      }catch(SAXException se) {
          se.printStackTrace(); 
          return null;
      }
      return motorPositions;
    }
 	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 * Grabbed from https://github.com/imagej/minimal-ij1-plugin
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Aberration_Calibration.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);
		// start ImageJ
		new ImageJ();
		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}

}

