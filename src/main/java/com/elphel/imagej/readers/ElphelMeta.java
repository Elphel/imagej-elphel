/**
 ** -----------------------------------------------------------------------------**
 ** ElphelTiffReader.java
 **
 ** Parse Elphel MakerNote and other Exif fields for both Tiff and JP4
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ElphelTiffReader.java is free software: you can redistribute it and/or modify
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
package com.elphel.imagej.readers;
import java.io.IOException;
import java.util.Hashtable;

//import ij.IJ;
import loci.formats.FormatException;


public class ElphelMeta {
//	private Hashtable<String, String> property_table = null;

	public static Hashtable<String, String> getMeta (Hashtable<String, String> property_table,
			long[] maker_note,
			double exposure,
			String date_time,
			int bytesPerPixel,
			boolean scale) throws FormatException, IOException {
		if (property_table == null) {
			property_table = new Hashtable<String, String> ();
		}
		if (!Double.isNaN(exposure)) {
			property_table.put("EXPOSURE",  String.format("%f",exposure));
		}
		if (date_time != null) {
			property_table.put("DATE_TIME",  date_time);
		}
		if (maker_note != null) {
			// copied from JP4_Reader_cam
			// Add GPS tags when there will be images to experiment (or while reimplementing JP4 reader)
			double[] gains= new double[4];
			double[] blacks= new double[4];
			double[] blacks256= new double[4];
			double[] gammas= new double[4];
			long []   gamma_scales= new long[4]; /* now not used, was scale _after_ gamma is applied, 0x400(default) corresponds to 1.0 */
			int i;
			double[][] rgammas=new double[4][];
			double min_gain;
			long WOI_LEFT,WOI_WIDTH,WOI_TOP,WOI_HEIGHT,BAYER_MODE,DCM_HOR,DCM_VERT,BIN_HOR,BIN_VERT;
			long COLOR_MODE=0;
			long FLIPH=0;
			long FLIPV=0;
			long  HEIGHT1=0;
			long  HEIGHT2=0;
			long  HEIGHT3=0;
			long  BLANK1=0;
			long  BLANK2=0;
			boolean FLIPH1=false;
			boolean FLIPV1=false;
			boolean FLIPH2=false;
			boolean FLIPV2=false;
			boolean FLIPH3=false;
			boolean FLIPV3=false;
			boolean COMPOSITE=false;
			boolean PORTRAIT=false;
			boolean YTABLEFORC=false;
			long     QUALITY=0;
			long     CQUALITY=0;
			long     CORING_INDEX_Y=0;
			long     CORING_INDEX_C=0;
			double []   satValue={255.0, 255.0, 255.0, 255.0};
			for (i=0;i<4;i++) { /* r,g,gb,b */
				gains[i]= maker_note[i]/65536.0;
				blacks[i]=(maker_note[i+4]>>24)/256.0;
				gammas[i]=((maker_note[i+4]>>16)&0xff)/100.0;
				gamma_scales[i]=maker_note[i+4] & 0xffff;
				property_table.put("gains_"+i,String.format("%f",gains[i]));
				property_table.put("blacks_"+i,String.format("%f",blacks[i]));
				property_table.put("gammas_"+i,String.format("%f",gammas[i]));
				property_table.put("gamma_scales_"+i,String.format("%d",gamma_scales[i]));
			}
			if (maker_note.length>=14) {
				COMPOSITE= ((maker_note[10] & 0xc0000000)!=0);
				if (COMPOSITE) {
					HEIGHT1= maker_note[11] & 0xffff;
					BLANK1= (maker_note[11]>>16) & 0xffff;
					HEIGHT2= maker_note[12] & 0xffff;
					BLANK2= (maker_note[12]>>16) & 0xffff;
					HEIGHT3=( maker_note[9]>>16) - HEIGHT1-BLANK1-HEIGHT2-BLANK2;

					FLIPH1= (((maker_note[10] >> 24) & 1)!=0); // Same value as FLIP_H
					FLIPV1= (((maker_note[10] >> 25) & 1)!=0); // Same value as FLIP_V
					FLIPH2= (((maker_note[10] >> 26) & 1)!=0);
					FLIPV2= (((maker_note[10] >> 27) & 1)!=0);
					FLIPH3= (((maker_note[10] >> 28) & 1)!=0);
					FLIPV3= (((maker_note[10] >> 29) & 1)!=0);
				}
				PORTRAIT=  (((maker_note[13] >>  7) & 1)!=0);
				YTABLEFORC=(((maker_note[13] >> 15) & 1)!=0);
				QUALITY=   (maker_note[13] & 0x7f);
				CQUALITY=  ((maker_note[13]>>8) & 0x7f);
				if (CQUALITY==0) CQUALITY=QUALITY;
				CORING_INDEX_Y= ((maker_note[13]>>16) & 0x7f);
				CORING_INDEX_C= ((maker_note[13]>>24) & 0x7f);
				if (CORING_INDEX_C==0) CORING_INDEX_C=CORING_INDEX_Y;

			}
			if (maker_note.length>=12) {
				WOI_LEFT=   maker_note[8]&0xffff;
				WOI_WIDTH=  maker_note[8]>>16;
				WOI_TOP=    maker_note[9]&0xffff;
				WOI_HEIGHT= maker_note[9]>>16;
				FLIPH=      maker_note[10]      & 1;
				FLIPV=     (maker_note[10]>> 1) & 1;
				BAYER_MODE=(maker_note[10]>> 2) & 3;
				COLOR_MODE=(maker_note[10]>> 4) & 0x0f;
				DCM_HOR=   (maker_note[10]>> 8) & 0x0f;
				DCM_VERT=  (maker_note[10]>>12) & 0x0f;
				BIN_HOR=   (maker_note[10]>>16) & 0x0f;
				BIN_VERT=  (maker_note[10]>>20) & 0x0f;
				property_table.put("WOI_LEFT",  String.format("%d",WOI_LEFT));
				property_table.put("WOI_WIDTH", String.format("%d",WOI_WIDTH));
				property_table.put("WOI_TOP",   String.format("%d",WOI_TOP));
				property_table.put("WOI_HEIGHT",String.format("%d",WOI_HEIGHT));
				property_table.put("FLIPH",     String.format("%d",FLIPH));
				property_table.put("FLIPV",     String.format("%d",FLIPV));
				property_table.put("BAYER_MODE",String.format("%d",BAYER_MODE));
				property_table.put("COLOR_MODE",((COLOR_MODE==2)?"JP46":((COLOR_MODE==5)?"JP4":((COLOR_MODE==0)?"MONO":((COLOR_MODE==15)?"RAW":"OTHER")))));
				property_table.put("DCM_HOR",   String.format("%d",DCM_HOR));
				property_table.put("DCM_VERT",  String.format("%d",DCM_VERT));
				property_table.put("BIN_HOR",   String.format("%d",BIN_HOR));
				property_table.put("BIN_VERT",  String.format("%d",BIN_VERT));

			}
			if (maker_note.length>=14) {
				property_table.put("COMPOSITE",String.format("%d",COMPOSITE?1:0));
				property_table.put("ORIENTATION",(PORTRAIT?"PORTRAIT":"LANDSCAPE" ));
				property_table.put("ORIENTATION",(YTABLEFORC?"1":"0"));
				property_table.put("QUALITY",String.format("%d",QUALITY)); //not full
				property_table.put("CORING_INDEX_Y",String.format("%d",CORING_INDEX_Y));
				property_table.put("CORING_INDEX_C",String.format("%d",CORING_INDEX_C));
			}
			if (maker_note.length>=16) {
				long [] iTemps={
						(maker_note[14]>> 0) & 0xffff,
						(maker_note[14]>>16) & 0xffff,
						(maker_note[15]>> 0) & 0xffff,
						(maker_note[15]>>16) & 0xffff};
				for (i=0;i<iTemps.length;i++) if (iTemps[i]!=0xffff){
					double temperature=(iTemps[i]&0xfff)/16.0;
					property_table.put("TEMPERATURE_"+i,""+temperature);

				}
			}
			if (COMPOSITE) {
				property_table.put("HEIGHT1",String.format("%d",HEIGHT1));
				property_table.put("HEIGHT2",String.format("%d",HEIGHT2));
				property_table.put("HEIGHT3",String.format("%d",HEIGHT3));
				property_table.put("BLANK_ROWS1",String.format("%d",BLANK1));
				property_table.put("BLANK_ROWS2",String.format("%d",BLANK2));
				property_table.put("FLIPH1",FLIPH1?"1":"0");
				property_table.put("FLIPH2",FLIPH2?"1":"0");
				property_table.put("FLIPH3",FLIPH3?"1":"0");
				property_table.put("FLIPV1",FLIPV1?"1":"0");
				property_table.put("FLIPV2",FLIPV2?"1":"0");
				property_table.put("FLIPV3",FLIPV3?"1":"0");
			}
			// If there are FLIPH, FLIPV - swap gains, gammas, blacks accordingly. later the images will be also flipped
			if (FLIPV!=0) {
				swapArrayElements (gains,       1, 3);
				swapArrayElements (gains,       0, 2);
				swapArrayElements (blacks,      1, 3);
				swapArrayElements (blacks,      0, 2);
				swapArrayElements (gammas,      1, 3);
				swapArrayElements (gammas,      0, 2);
				swapArrayElements (gamma_scales,1, 3);
				swapArrayElements (gamma_scales,0, 2);
			}
			if (FLIPH!=0) {
				swapArrayElements (gains,       1, 0);
				swapArrayElements (gains,       3, 2);
				swapArrayElements (blacks,      1, 0);
				swapArrayElements (blacks,      3, 2);
				swapArrayElements (gammas,      1, 0);
				swapArrayElements (gammas,      3, 2);
				swapArrayElements (gamma_scales,1, 0);
				swapArrayElements (gamma_scales,3, 2);
			}
			for (i=0;i<4;i++) rgammas[i]=elphel_gamma_calc (gammas[i], blacks[i], gamma_scales[i]);
			/**adjusting gains to have the result picture in the range 0..256 */
			min_gain=2.0*gains[0];
			for (i=0;i<4;i++) {
				if (bytesPerPixel == 1) {
					if (min_gain > gains[i]*(1.0-blacks[i])) min_gain = gains[i]*(1.0-blacks[i]);
				} else {
					if (min_gain > gains[i]) min_gain = gains[i];

				}
			}
			property_table.put("GAIN",String.format("%f",min_gain)); // common gain

			for (i=0;i<4;i++) gains[i]/=min_gain;
			for (i=0;i<4;i++) blacks256[i]=256.0*blacks[i];


			for (i=0;i<4;i++) {
				if  (maker_note !=null) {
					if (scale) satValue[i]=((rgammas[i][255])-blacks256[i])/gains[i];
					else       satValue[i]=((rgammas[i][255])-blacks256[i]);
				}   else       satValue[i]=255.0;
				property_table.put("saturation_"+i,String.format("%f",satValue[i]));

			}
			// swap satValue to match FLIPH,FLIPV again
			if (FLIPV!=0) {
				swapArrayElements (satValue,       1, 3);
				swapArrayElements (satValue,       0, 2);
			}
			if (FLIPH!=0) {
				swapArrayElements (satValue,       1, 0);
				swapArrayElements (satValue,       3, 2);
			}
			for (i=0;i<4;i++) {
				property_table.put("saturation_"+i,String.format("%f",satValue[i]));
			}
		}
		return property_table;
	}





	//	public Hashtable<String, String> getPropertyTable(){
	//		return property_table;
	//	}

	// -- Helper methods --

	static void swapArrayElements (double[]arr,int i, int j) {
		double tmp=arr[i];
		arr[i]=arr[j];
		arr[j]=tmp;
	}
	static void swapArrayElements (long[]arr,int i, int j) {
		long tmp=arr[i];
		arr[i]=arr[j];
		arr[j]=tmp;
	}
	/* reverses gamma calculations in the camera
    returns double[] table , in the range 0.0..255.996
	 */
	static double [] elphel_gamma_calc (double gamma, double black, long gamma_scale) {
		int i;
		double x, black256 ,k;
		int[] gtable = new int[257];
		double[] rgtable =new double[256];
		int ig;
		black256=black*256.0;
		k=1.0/(256.0-black256);
		if (gamma < 0.13) gamma=0.13;
		if (gamma >10.0)  gamma=10.0;
		for (i=0; i<257; i++) {
			x=k*(i-black256);
			if (x < 0.0 ) x=0.0;
			ig= (int) (0.5+65535.0*Math.pow(x,gamma));
			ig=(ig* (int) gamma_scale)/0x400;
			if (ig > 0xffff) ig=0xffff;
			gtable[i]=ig;
		}
		/* now gtable[] is the same as was used in the camera */
		/* FPGA was using linear interpolation between elements of the gamma table, so now we'll reverse that process */
		// double[] rgtable =new  double[256];
		int indx=0;
		double outValue;
		for (i=0; i<256; i++ ) {
			outValue=128+(i<<8);
			while ((gtable[indx+1]<outValue) && (indx<256)) indx++;
			if (indx>=256) rgtable[i]=65535.0/256;
			else if (gtable[indx+1]==gtable[indx]) rgtable[i]=i;
			else           rgtable[i]=indx+(1.0*(outValue-gtable[indx]))/(gtable[indx+1] - gtable[indx]);
		}
		return rgtable;
	}



}
