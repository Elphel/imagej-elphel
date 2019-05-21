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
	public static final String REVISION =  "REVISION";   public static final int REVISION_OFFSET =   0;
	public static final String UPTIME =    "UPTIME";     public static final int UPTIME_OFFSET =     1;
	public static final String STATUS =    "STATUS";     public static final int STATUS_OFFSET =     3;
	public static final String MOD_SER =   "MOD_SER";    public static final int MOD_SER_OFFSET =    5;
	public static final String SOFT_REV =  "SOFT_REV";   public static final int SOFT_REV_OFFSET =  13;
	public static final String FRAME =     "FRAME";      public static final int FRAME_OFFSET =     20;
	public static final String MEAN =      "MEAN";       public static final int MEAN_OFFSET =      22;
	public static final String FPA_TCNT =  "FPA_TCNT";   public static final int FPA_TCNT_OFFSET =  23;
	public static final String FPA_KELV =  "FPA_KELV";   public static final int FPA_KELV_OFFSET =  24;
	public static final String ENCL_TCNT = "ENCL_TCNT";  public static final int ENCL_TCNT_OFFSET = 25;
	public static final String ENCL_KELV = "ENCL_KELV";  public static final int ENCL_KELV_OFFSET = 26;
	public static final String FFC_KELV =  "FFC_KELV";   public static final int FFC_KELV_OFFSET =  29;


	private float [] pixels =  null;
	private byte  [] btm =     null;
	private ByteBuffer bbtm  = null;

	public Lepton3Telemetry(ByteBuffer bb, float [] pixels, boolean bottom) {
		this.pixels = pixels;
		int offset = 2*(bottom?(LEPTON3_IMAGE_WIDTH * LEPTON3_IMAGE_HEIGHT):0);
		if (pixels.length > LEPTON3_IMAGE_WIDTH*LEPTON3_IMAGE_HEIGHT) {
			this.btm = new byte [LEPTON3_IMAGE_WIDTH* LEPTON3_TELEMETRY_LINES * 2];
			bb.position(offset);
			this.bbtm = bb.get(btm, 0,btm.length);
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
		tm.put(REVISION,  ""+bytes[1]+"."+bytes[0]);
		tm.put(UPTIME,    ""+getU32(UPTIME_OFFSET));
		tm.put(STATUS,    ""+getU32(STATUS_OFFSET));
		tm.put(MOD_SER,   ""+getHex(MOD_SER_OFFSET,8));
		tm.put(SOFT_REV,  ""+getHex(SOFT_REV_OFFSET,4));
		tm.put(FPA_TCNT,  ""+getU16(FPA_TCNT_OFFSET));
		tm.put(FPA_KELV,  ""+getU16(FPA_KELV_OFFSET));
		tm.put(ENCL_TCNT, ""+getU16(ENCL_TCNT_OFFSET));
		tm.put(ENCL_KELV, ""+getU16(ENCL_KELV_OFFSET));
		tm.put(FFC_KELV,  ""+getU16(FFC_KELV_OFFSET));

/*
	public static final String FPA_TCNT =  "FPA_TCNT";   public static final int FPA_TCNT_OFFSET =  23;
	public static final String FPA_KELV =  "FPA_KELV";   public static final int FPA_KELV_OFFSET =  24;
	public static final String ENCL_TCNT = "ENCL_TCNT";  public static final int ENCL_TCNT_OFFSET = 25;
	public static final String ENCL_KELV = "ENCL_KELV";  public static final int ENCL_KELV_OFFSET = 26;
	public static final String FFC_KELV =  "FFC_KELV";   public static final int FFC_KELV_OFFSET =  29;

 */
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

	private String getString (int offset, int len) { // offset in shorts, len in shorts
		return new String(getBytes(offset, len));
	}
	private String getHex (int offset, int len) { // offset in shorts, len in shorts
		String s = "";
		for (byte b :getBytes(offset, len)) {
			s += String.format("%02x", b);
		}
		return s;
	}

	private int getU16(int offset) {
		return (bbtm.getShort(2 * offset)) & 0xffff;
	}

	private long getU32(int offset) {
		return ((long) (bbtm.getShort(2 * offset)) & 0xffff) + (((long) (bbtm.getShort(2 * (offset+1))) & 0xffff) << 16);
	}


}
