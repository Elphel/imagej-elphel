/**
 ** -----------------------------------------------------------------------------**
 ** Lepton3Telemetry.java
 **
 ** Parses FLIR Lepton3 telemetry data
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  Lepton3Telemetry.java is free software: you can redistribute it and/or modify
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

public class Lepton3Telemetry {
	public static final int LEPTON3_IMAGE_WIDTH =     160;
	public static final int LEPTON3_IMAGE_HEIGHT =    120;
	public static final int LEPTON3_TELEMETRY_LINES =   2;
	public static final String REVISION =      "REVISION";      public static final int REVISION_OFFSET =       0;
	public static final String UPTIME =        "UPTIME";        public static final int UPTIME_OFFSET =         1;
	public static final String STATUS =        "STATUS";        public static final int STATUS_OFFSET =         3;
	public static final String MOD_SER =       "MOD_SER";       public static final int MOD_SER_OFFSET =        5;
	public static final String SOFT_REV =      "SOFT_REV";      public static final int SOFT_REV_OFFSET =      13;
	public static final String FRAME =         "FRAME";         public static final int FRAME_OFFSET =         20;
	public static final String MEAN =          "MEAN";          public static final int MEAN_OFFSET =          22;
	public static final String FPA_TCNT =      "FPA_TCNT";      public static final int FPA_TCNT_OFFSET =      23;
	public static final String FPA_KELV =      "FPA_KELV";      public static final int FPA_KELV_OFFSET =      24;
	public static final String ENCL_TCNT =     "ENCL_TCNT";     public static final int ENCL_TCNT_OFFSET =     25;
	public static final String ENCL_KELV =     "ENCL_KELV";     public static final int ENCL_KELV_OFFSET =     26;
	public static final String FFC_KELV =      "FFC_KELV";      public static final int FFC_KELV_OFFSET =      29;
	public static final String FFC_TIME =      "FFC_TIME";      public static final int FFC_TIME_OFFSET =      30;
	public static final String ENCL_FFC_KELV = "ENCL_FFC_KELV"; public static final int ENCL_FFC_KELV_OFFSET = 32;
	public static final String AGC_ROI_TLBR =  "AGC_ROI_TLBR";  public static final int AGC_ROI_TLBR_OFFSET =  34;
	public static final String AGC_CLIP_HIGH = "AGC_CLIP_HIGH"; public static final int AGC_CLIP_HIGH_OFFSET = 38;
	public static final String AGC_CLIP_LOW =  "AGC_CLIP_LOW";  public static final int AGC_CLIP_LOW_OFFSET =  38;
	public static final String VFORMAT =       "VFORMAT";       public static final int VFORMAT_OFFSET =       72;
	public static final String FFC_LOG2 =      "FFC_LOG2";      public static final int FFC_LOG2_OFFSET =      74;
	public static final String EMISS_8192 =    "EMISS_8192";    public static final int EMISS_8192_OFFSET =    19 + 80; // row B
	public static final String BGND_KELV =     "BGND_KELV";     public static final int BGND_KELV_OFFSET =     20 + 80;
	public static final String ATMOSPH_8192 =  "ATMOSPH_8192";  public static final int ATMOSPH_8192_OFFSET =  21 + 80;
	public static final String ATMOSPH_KELV =  "ATMOSPH_KELV";  public static final int ATMOSPH_KELV_OFFSET =  22 + 80;
	public static final String WND_TRANS_8192 ="WND_TRANS_8192";public static final int WND_TRANS_8192_OFFSET =23 + 80;
	public static final String WND_REFL_8192 = "WND_REFL_8192"; public static final int WND_REFL_8192_OFFSET = 24 + 80;
	public static final String WND_KELV =      "WND_KELV";      public static final int WND_KELV_OFFSET =      25 + 80;
	public static final String WND_REFL_KELV = "WND_REFL_KELV"; public static final int WND_REFL_KELV_OFFSET = 26 + 80;
	public static final String GAIN_MODE =     "GAIN_MODE";     public static final int GAIN_MODE_OFFSET =      5 + 160; // row C
	public static final String GAIN_EFF =      "GAIN_EFF";      public static final int GAIN_EFF_OFFSET =       6 + 160;
	public static final String GAIN_DFLAG =    "GAIN_DFLAG";    public static final int GAIN_DFLAG_OFFSET =     7 + 160;
	public static final String GAIN_THR_HL_C = "GAIN_THR_HL_C"; public static final int GAIN_THR_HL_C_OFFSET =  8 + 160;
	public static final String GAIN_THR_LH_C = "GAIN_THR_LH_C"; public static final int GAIN_THR_LH_C_OFFSET =  9 + 160;
	public static final String GAIN_THR_HL_K = "GAIN_THR_HL_K"; public static final int GAIN_THR_HL_K_OFFSET = 10 + 160;
	public static final String GAIN_THR_LH_K = "GAIN_THR_LH_K"; public static final int GAIN_THR_LH_K_OFFSET = 11 + 160;
	public static final String GAIN_HL_PP =    "GAIN_HL_PP";    public static final int GAIN_HL_PP_OFFSET =    14 + 160;
	public static final String GAIN_LH_PP =    "GAIN_LH_PP";    public static final int GAIN_LH_PP_OFFSET =    15 + 160;
	public static final String GAIN_ROI_TLBR = "GAIN_ROI_TLBR"; public static final int GAIN_ROI_TLBR_OFFSET = 22 + 160;
	public static final String TLIN_EN =       "TLIN_EN";       public static final int TLIN_EN_OFFSET =       48 + 160;
	public static final String TLIN_RESOL =    "TLIN_RESOL";    public static final int TLIN_RESOL_OFFSET =    49 + 160;
	public static final String SPOT_AVG_KELV = "SPOT_AVG_KELV"; public static final int SPOT_AVG_KELV_OFFSET = 50 + 160;
	public static final String SPOT_MAX_KELV = "SPOT_MAX_KELV"; public static final int SPOT_MAX_KELV_OFFSET = 51 + 160;
	public static final String SPOT_MIN_KELV = "SPOT_MIN_KELV"; public static final int SPOT_MIN_KELV_OFFSET = 52 + 160;
	public static final String SPOT_POP_PERC = "SPOT_POP_PERC"; public static final int SPOT_POP_PERC_OFFSET = 53 + 160;
	public static final String SPOT_ROI_TLBR = "SPOT_ROI_TLBR"; public static final int SPOT_ROI_TLBR_OFFSET = 54 + 160;



	private float [] pixels =  null;
	private byte  [] btm =     null;
	private ByteBuffer bbtm  = null;

	public Lepton3Telemetry(ByteBuffer bb, float [] pixels, boolean bottom) {
		this.pixels = pixels;
		int offset = 2*(bottom?(LEPTON3_IMAGE_WIDTH * LEPTON3_IMAGE_HEIGHT):0);
		if (pixels.length > LEPTON3_IMAGE_WIDTH*LEPTON3_IMAGE_HEIGHT) {
			this.btm = new byte [LEPTON3_IMAGE_WIDTH* LEPTON3_TELEMETRY_LINES * 2];
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
			this.pixels = new float [LEPTON3_IMAGE_WIDTH * LEPTON3_IMAGE_HEIGHT];
			System.arraycopy(
					pixels,
					(bottom? 0: LEPTON3_IMAGE_WIDTH * LEPTON3_TELEMETRY_LINES),
					this.pixels,
					0,
					this.pixels.length);
		}
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
		tm.put(REVISION,       "" + bytes[1] + "." + bytes[0]);
		tm.put(UPTIME,         String.format("%.3f", getSeconds (UPTIME_OFFSET)));
		tm.put(STATUS,         "" + getU32     (STATUS_OFFSET));
		tm.put(MOD_SER,        "" + getHex     (MOD_SER_OFFSET,8));
		tm.put(SOFT_REV,       "" + getHex     (SOFT_REV_OFFSET,4));
		tm.put(FPA_TCNT,       "" + getU16     (FPA_TCNT_OFFSET));
		tm.put(FPA_KELV,       String.format("%.2f", getKelv    (FPA_KELV_OFFSET)));
		tm.put(ENCL_TCNT,      "" + getU16     (ENCL_TCNT_OFFSET));
		tm.put(ENCL_KELV,      String.format("%.2f", getKelv    (ENCL_KELV_OFFSET)));
		tm.put(FFC_KELV,       String.format("%.2f", getKelv    (FFC_KELV_OFFSET)));
		tm.put(FFC_TIME,       String.format("%.3f", getSeconds (FFC_TIME_OFFSET)));

		tm.put(ENCL_FFC_KELV,  String.format("%.2f", getKelv    (ENCL_FFC_KELV_OFFSET)));
		tm.put(AGC_ROI_TLBR,
				"" +  getU16 (AGC_ROI_TLBR_OFFSET) +
				" " + getU16 (AGC_ROI_TLBR_OFFSET + 1) +
				" " + getU16 (AGC_ROI_TLBR_OFFSET + 2) +
				" " + getU16 (AGC_ROI_TLBR_OFFSET + 3));
		tm.put(AGC_CLIP_HIGH,  "" + getU16     (AGC_CLIP_HIGH_OFFSET));
		tm.put(AGC_CLIP_LOW,   "" + getU16     (AGC_CLIP_LOW_OFFSET));
		tm.put(VFORMAT,        "" + getU32     (VFORMAT_OFFSET));
		tm.put(FFC_LOG2,       "" + getU16     (FFC_LOG2_OFFSET));

		tm.put(EMISS_8192,     "" + getU16     (EMISS_8192_OFFSET));
		tm.put(BGND_KELV,      String.format("%.2f", getKelv    (BGND_KELV_OFFSET)));
		tm.put(ATMOSPH_8192,   "" + getU16     (ATMOSPH_8192_OFFSET));
		tm.put(ATMOSPH_KELV,   String.format("%.2f", getKelv    (ATMOSPH_KELV_OFFSET)));
		tm.put(WND_TRANS_8192, "" + getU16     (WND_TRANS_8192_OFFSET));
		tm.put(WND_REFL_8192,  "" + getU16     (WND_REFL_8192_OFFSET));
		tm.put(WND_KELV,       String.format("%.2f", getKelv    (WND_KELV_OFFSET)));
		tm.put(WND_REFL_KELV,  String.format("%.2f", getKelv    (WND_REFL_KELV_OFFSET)));
		tm.put(GAIN_MODE,      "" + getU16     (GAIN_MODE_OFFSET));
		tm.put(GAIN_EFF,       "" + getU16     (GAIN_EFF_OFFSET));
		tm.put(GAIN_DFLAG,     "" + getU16     (GAIN_DFLAG_OFFSET));
		tm.put(GAIN_THR_HL_C,  "" + getU16     (GAIN_THR_HL_C_OFFSET));
		tm.put(GAIN_THR_LH_C,  "" + getU16     (GAIN_THR_LH_C_OFFSET));
		tm.put(GAIN_THR_HL_K,  "" + getU16     (GAIN_THR_HL_K_OFFSET));
		tm.put(GAIN_THR_LH_K,  "" + getU16     (GAIN_THR_LH_K_OFFSET));
		tm.put(GAIN_HL_PP,     "" + getU16     (GAIN_HL_PP_OFFSET));
		tm.put(GAIN_LH_PP,     "" + getU16     (GAIN_LH_PP_OFFSET));
		tm.put(GAIN_ROI_TLBR,
				"" +  getU16 (GAIN_ROI_TLBR_OFFSET) +
				" " + getU16 (GAIN_ROI_TLBR_OFFSET + 1) +
				" " + getU16 (GAIN_ROI_TLBR_OFFSET + 2) +
				" " + getU16 (GAIN_ROI_TLBR_OFFSET + 3));
		tm.put(TLIN_EN,        "" + getU16 (TLIN_EN_OFFSET));
		tm.put(TLIN_RESOL,     "" + getU16     (TLIN_RESOL_OFFSET));
		tm.put(SPOT_AVG_KELV,  String.format("%.2f", getKelv    (SPOT_AVG_KELV_OFFSET)));
		tm.put(SPOT_MAX_KELV,  String.format("%.2f", getKelv    (SPOT_MAX_KELV_OFFSET)));
		tm.put(SPOT_MIN_KELV,  String.format("%.2f", getKelv    (SPOT_MIN_KELV_OFFSET)));
		tm.put(SPOT_POP_PERC,  "" + getU16     (SPOT_POP_PERC_OFFSET));
		tm.put(SPOT_ROI_TLBR,
				"" +  getU16 (SPOT_ROI_TLBR_OFFSET) +
				" " + getU16 (SPOT_ROI_TLBR_OFFSET + 1) +
				" " + getU16 (SPOT_ROI_TLBR_OFFSET + 2) +
				" " + getU16 (SPOT_ROI_TLBR_OFFSET + 3));
		return tm;
	}
	// internal methods
	// input byte array depends on tiff endian for short (16-bit) words, but FLIR word sequence is little endian
	private byte[] getBytes (int offset, int len) { // offset in shorts, len in shorts
		byte [] bytes = new byte [2*len];
		for (int i = 0; i<len; i++) {
			int d = (bbtm.getShort(2 * (offset + i))) & 0xffff;
			bytes[2*i+0] = (byte) (d & 0xff);
			bytes[2*i+1] = (byte) ((d >> 8) & 0xff);
		}
		return bytes;
	}

//	private String getString (int offset, int len) { // offset in shorts, len in shorts
//		return new String(getBytes(offset, len));
//	}

	private String getHex (int offset, int len) { // offset in shorts, len in shorts
		String s = "";
		for (byte b :getBytes(offset, len)) {
			s += String.format("%02x", b);
		}
		return s;
	}

	private double getKelv(int offset) {
		return 0.01*((bbtm.getShort(2 * offset)) & 0xffff);
	}
	private int getU16(int offset) {
		return (bbtm.getShort(2 * offset)) & 0xffff;
	}

	private double getSeconds(int offset) {
		return 0.001*(((long) (bbtm.getShort(2 * offset)) & 0xffff) + (((long) (bbtm.getShort(2 * (offset+1))) & 0xffff) << 16));
	}

	private long getU32(int offset) {
		return ((long) (bbtm.getShort(2 * offset)) & 0xffff) + (((long) (bbtm.getShort(2 * (offset+1))) & 0xffff) << 16);
	}


}
