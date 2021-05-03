/**
 ** -----------------------------------------------------------------------------**
 ** Boson640Telemetry.java
 **
 ** Parses FLIR Boson640 telemetry data
 **
 ** Copyright (C) 2021 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Boson640Telemetry.java is free software: you can redistribute it and/or modify
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

import java.nio.ByteBuffer;
import java.util.HashMap;

public class Boson640Telemetry {
	private static final int BOSON640_IMAGE_WIDTH =     640;
	private static final int BOSON640_IMAGE_HEIGHT =    512;
	private static final int BOSON640_TELEMETRY_LINES =   1;
	private static final String REVISION =      "REVISION";      private static final int REVISION_OFFSET =       0; //  2
	private static final String CAMSERIAL =     "CAMSERIAL";     private static final int CAMSERIAL_OFFSET =      2; //  4
	private static final String SENSSERIAL =    "SENSSERIAL";    private static final int SENSSERIAL_OFFSET =     6; //  4
	private static final String PART_ASCII =    "PART_ASCII";     private static final int PART_ASCII_OFFSET =     10; // 20

	private static final String SOFT_REV =      "SOFT_REV";      private static final int SOFT_REV_OFFSET =      44; // 12
	private static final String FRAME_RATE =    "FRAME_RATE";    private static final int FRAME_RATE_OFFSET =    56; //  2
	private static final String STATUS =        "STATUS";        private static final int STATUS_OFFSET =        76; // hex,8
	private static final String FRAME =         "FRAME";         private static final int FRAME_OFFSET =         84; // 4
	private static final String FRAME_FFC =     "FRAME_FFC";     private static final int FRAME_FFC_OFFSET =     88; // 4
	private static final String FPA_KELV =      "FPA_KELV";      private static final int FPA_KELV_OFFSET =      94; // 2
	private static final String FFC_KELV =      "FFC_KELV";      private static final int FFC_KELV_OFFSET =      96; // 2
	private static final String PIPELINE =      "PIPELINE";      private static final int PIPELINE_OFFSET =     110; // 4
	private static final String FFC_INTEG =     "FFC_INTEG";     private static final int FFC_INTEG_OFFSET =    114; // 2
	private static final String NUC_CURRENT =   "NUC_CURRENT";   private static final int NUC_CURRENT_OFFSET =  158; // 2
	private static final String NUC_DESIRED =   "NUC_DESIRED";   private static final int NUC_DESIRED_OFFSET =  160; // 2
	private static final String CORE_TEMP =     "CORE_TEMP";     private static final int CORE_TEMP_OFFSET =    162; // 4
	private static final String OVERTEMP =      "OVERTEMP";      private static final int OVERTEMP_OFFSET =     166; // 4
	private static final String ROI_POP_LTH =   "ROI_POP_LTH";   private static final int ROI_POP_LTH_OFFSET =  170; // 4
	private static final String ROI_POP_HTL =   "ROI_POP_HTL";   private static final int ROI_POP_HTL_OFFSET =  174; // 4
	private static final String TGL_PATT =      "TGL_PATT";      private static final int TGL_PATT_OFFSET =     178; // 6
	private static final String ZOOM_FACT =     "ZOOM_FACT";     private static final int ZOOM_FACT_OFFSET =    184; // 4
	private static final String ZOOM_X0 =       "ZOOM_X0";       private static final int ZOOM_X0_OFFSET =      188; // 4
	private static final String ZOOM_Y0 =       "ZOOM_Y0";       private static final int ZOOM_Y0_OFFSET =      192; // 4



	private float [] pixels =  null;
	private byte  [] btm =     null;
	private ByteBuffer bbtm  = null;

	public Boson640Telemetry(ByteBuffer bb, float [] pixels, boolean bottom) {
		this.pixels = pixels;
		int offset = 2*(bottom?(BOSON640_IMAGE_WIDTH * BOSON640_IMAGE_HEIGHT):0);
		if (pixels.length > BOSON640_IMAGE_WIDTH*BOSON640_IMAGE_HEIGHT) {
			this.btm = new byte [BOSON640_IMAGE_WIDTH* BOSON640_TELEMETRY_LINES * 2];
			bb.position(offset);
			byte dbg = bb.get();
			dbg = bb.get();
			dbg = bb.get();
			dbg = bb.get();
			dbg = bb.get();
			bb.position(offset);

//			this.bbtm =
			bb.get(btm, 0,btm.length);
			this.bbtm = ByteBuffer.wrap(btm);
			this.bbtm.order( bb.order()); // or is it already same as in bb?
			
			System.out.println("getShort(94)"+this.bbtm.getShort(94));
			System.out.println("getLong(162)"+this.bbtm.getLong(162));
			System.out.println("getLong(162)"+this.bbtm.getInt(162));
			
			
			
			
			this.pixels = new float [BOSON640_IMAGE_WIDTH * BOSON640_IMAGE_HEIGHT];
			System.arraycopy(
					pixels,
					(bottom? 0: BOSON640_IMAGE_WIDTH * BOSON640_TELEMETRY_LINES),
					this.pixels,
					0,
					this.pixels.length);
		}
	}

	public static int getWidth() {
		return BOSON640_IMAGE_WIDTH;
	}

	public static int getHeight() {
		return BOSON640_IMAGE_HEIGHT;
	}
	
	public static int getTelemetryLines() {
		return BOSON640_TELEMETRY_LINES;
	}
	
	
	public boolean hasTelemetry() {
		return bbtm != null;
	}

	public float [] getPixels() {
		return pixels;
	}

	public HashMap<String, String> parseTelemetry(){
		HashMap<String, String> tm = new HashMap<String, String>();
		byte [] bytes = getBytes(REVISION_OFFSET,2);
		tm.put(REVISION,      "" + bytes[0] + "." + bytes[1]); // "0.2"
		tm.put(CAMSERIAL,     "" + String.format("%d", getU32 (CAMSERIAL_OFFSET)));
		tm.put(SENSSERIAL,    "" + String.format("%d", getU32 (SENSSERIAL_OFFSET)));
		tm.put(PART_ASCII,    "" + getAscii(PART_ASCII_OFFSET,20));
		tm.put(STATUS,        "" + (bbtm.getLong(STATUS_OFFSET))); // getHex (STATUS_OFFSET,8));
		tm.put(SOFT_REV,      "" + getU32 (SOFT_REV_OFFSET)+"."+getU32 (SOFT_REV_OFFSET+4)+"."+getU32 (SOFT_REV_OFFSET+8));
		tm.put(FRAME_RATE,    "" + getU16 (FRAME_RATE_OFFSET));
		tm.put(FRAME,         "" + getU32 (FRAME_OFFSET));
		tm.put(FRAME_FFC,     "" + getU32 (FRAME_FFC_OFFSET));
		tm.put(FPA_KELV,      String.format("%.1f", getKelv10 (FPA_KELV_OFFSET)));
		tm.put(FFC_KELV,      String.format("%.1f", getKelv10 (FFC_KELV_OFFSET)));
		tm.put(PIPELINE,      "" + getU32 (PIPELINE_OFFSET));
		tm.put(FFC_INTEG,     "" + getU16 (FFC_INTEG_OFFSET));
		tm.put(NUC_CURRENT,   "" + getU16 (NUC_CURRENT_OFFSET));
		tm.put(NUC_DESIRED,   "" + getU16 (NUC_DESIRED_OFFSET));
		tm.put(CORE_TEMP,     String.format("%.3f", getCelsius1000 (CORE_TEMP_OFFSET)));
		tm.put(OVERTEMP,      "" + getU32 (OVERTEMP_OFFSET));
		tm.put(ROI_POP_LTH,   "" + getU32 (ROI_POP_LTH_OFFSET));
		tm.put(ROI_POP_HTL,   "" + getU32 (ROI_POP_HTL_OFFSET));
		tm.put(TGL_PATT,      "" + getHex (TGL_PATT_OFFSET,8));
		tm.put(ZOOM_FACT,     "" + getU32 (ZOOM_FACT_OFFSET));
		tm.put(ZOOM_X0,       "" + getU32 (ZOOM_X0_OFFSET));
		tm.put(ZOOM_Y0,       "" + getU32 (ZOOM_Y0_OFFSET));
		return tm;
	}
	// internal methods

	private long getU32(int offset) {
		return ((long) (bbtm.getInt(offset))) & 0xffffffff;
	}

	private long getU16(int offset) {
		return ((long) (bbtm.getShort(offset))) & 0xffff;
	}

	private byte[] getBytes (int offset, int len) { // offset, len in bytes!
		byte [] bytes = new byte [len];
		for (int i = 0; i<len; i++) {
			bytes[i] = bbtm.get(offset+i);
		}
		return bytes;
	}
	private String getHex (int offset, int len) { // offset, len bytes,
		String s = "";
		for (byte b :getBytes(offset, len)) {
			s += String.format("%02x", b);
		}
		return s;
	}

	private String getAscii (int offset, int len) { // offset, len bytes,
		String s = "";
		for (byte b :getBytes(offset, len)) {
			if (b == 0) break;
			s += String.format("%c", b);
		}
		return s;
	}
	
	private double getKelv10(int offset) {
		return 0.1*getU16(offset);
	}
	private double getCelsius1000(int offset) {
		return 0.001*bbtm.getInt(offset);
	}
	

}
