package com.elphel.imagej.gpu;

public class TpTask {
	public int        task; // [0](+1) - generate 4 images, [4..9]+16..+512 - correlation pairs, 2 - generate texture tiles
	public float      target_disparity;
    public int        num_sensors = 4;
	public int        ty;
	public int        tx;
	public float[][]  xy = null;
	public float[][]  xy_aux = null;
	public float [][] disp_dist = null;
//	public float      weight;
	
	public TpTask(
			int num_sensors,
			int tileX,
			int tileY) {
		this.num_sensors = num_sensors;
		this.ty = tileY;
		this.tx = tileX;
	}

	public TpTask(int num_sensors, int tx, int ty, float target_disparity, int task ) {
		this.tx = tx;
		this.ty = ty;
		this.target_disparity = target_disparity;
		this.task = task;
		this.num_sensors = num_sensors; // will not be encoded
		this.disp_dist = new float [num_sensors][4];
	}
	/**
	 * Initialize from the float array (read from the GPU)
	 * @param flt float array containing tasks data
	 * @param indx task number to use
	 */
	public TpTask(float [] flt, int indx, boolean use_aux)
	{
		task =    Float.floatToIntBits(flt[indx++]);
		int txy = Float.floatToIntBits(flt[indx++]);
		ty = txy >> 16;
		tx = txy & 0xffff;
		if (use_aux) {
			xy_aux = new float[num_sensors][2];
    		for (int i = 0; i < num_sensors; i++) {
    			xy_aux[i][0] = flt[indx++];
    			xy_aux[i][1] = flt[indx++];
    		}
		} else {
			xy = new float[num_sensors][2];
    		for (int i = 0; i < num_sensors; i++) {
    			xy[i][0] = flt[indx++];
    			xy[i][1] = flt[indx++];
    		}
		}
		target_disparity = flt[indx++];
		disp_dist = new float [num_sensors][4];
		for (int i = 0; i < num_sensors; i++) {
			for (int j = 0; j < 4; j++) {
				disp_dist[i][j] = flt[indx++];
			}
		}
	}
	public float [][] getDispDist(){
		return disp_dist;
	}
	public void setDispDist(float [][] disp_dist){
		this.disp_dist = disp_dist;
	}
	public void setDispDist(double [][] disp_dist){
		this.disp_dist = new float [disp_dist.length][disp_dist[0].length];
		for (int i = 0; i < disp_dist.length; i++) {
			for (int j = 0; j < disp_dist[0].length; j++) {
				this.disp_dist[i][j] = (float) disp_dist[i][j];		
			}
		}
	}

	public double [][] getDoubleDispDist(){
		if (disp_dist == null) { // can it happen?
			return null;
		}
		double [][] ddisp_dist = new double [disp_dist.length][disp_dist[0].length];
		for (int nsens = 0; nsens < disp_dist.length; nsens++) {
			for (int i = 0; i < disp_dist[nsens].length; i++) {
				ddisp_dist[nsens][i] = disp_dist[nsens][i]; 
			}
		}
		return ddisp_dist;
	}
	@Deprecated
	public float [][] getXY(boolean use_aux){
		return use_aux? xy_aux : xy;
	}
	@Deprecated
	public double [][] getDoubleXY(boolean use_aux){
		float [][] fXY = getXY(use_aux);
		if (fXY == null) {
			return null;
		}
		double [][] dXY = new double [fXY.length][fXY[0].length];
		for (int nsens = 0; nsens < fXY.length; nsens++) {
			for (int i = 0; i < fXY[nsens].length; i++) {
				dXY[nsens][i] = fXY[nsens][i]; 
			}
		}
		return dXY;
	}

	public float [][] getXY(){
		return  xy;
	}
	public double [][] getDoubleXY(){
		float [][] fXY = getXY();
		if (fXY == null) {
			return null;
		}
		double [][] dXY = new double [fXY.length][fXY[0].length];
		for (int nsens = 0; nsens < fXY.length; nsens++) {
			for (int i = 0; i < fXY[nsens].length; i++) {
				dXY[nsens][i] = fXY[nsens][i]; 
			}
		}
		return dXY;
	}
	
	
	
	public int getTileY(){
		return ty;
	}
	public int getTileX(){
		return tx;
	}
	public int getTask() {
		return task;
	}
	public double getTargetDisparity() {
		return target_disparity;
	}

	// convert this class instance to float array to match layout of the C struct
	public float [] asFloatArray(boolean use_aux) {
		float [] flt = new float [GPUTileProcessor.TPTASK_SIZE];
		return asFloatArray(flt, 0, use_aux);
	}
	// convert this class instance to float array to match layout of the C struct,
	// fill existing float array from the specified index
	public float [] asFloatArray(float [] flt, int indx, boolean use_aux) {
		flt[indx++] = Float.intBitsToFloat(task);
		flt[indx++] = Float.intBitsToFloat(tx + (ty << 16));
		float [][] offsets = use_aux? this.xy_aux: this.xy;
		for (int i = 0; i < num_sensors; i++) {
			if (offsets != null) {
				flt[indx++] = offsets[i][0];
				flt[indx++] = offsets[i][1];
			} else {
				indx+= 2;
			}
		}
		flt[indx++] = this.target_disparity;
		/*
		for (int i = 0; i < NUM_CAMS; i++) { // actually disp_dist will be initialized by the GPU
			indx+= 4;
			flt[indx++] = disp_dist[i][0];
			flt[indx++] = disp_dist[i][1];
			flt[indx++] = disp_dist[i][2];
			flt[indx++] = disp_dist[i][3];
		}
		*/
		return flt;
	}
}