package com.elphel.imagej.gpu;

public class TpTask {
	public int        task; // [0](+1) - generate 4 images, [4..9]+16..+512 - correlation pairs, 2 - generate texture tiles
	public float      target_disparity;
    public int        num_sensors = 4;
	public int        ty;
	public int        tx;
	public float []   centerXY = {Float.NaN,Float.NaN};
	public float[][]  xy = null;
	public float[][]  xy_aux = null;
	public float [][] disp_dist = null;
//	public float      weight;
	public static int getSize(int num_sensors) {
		return 5 + 2* num_sensors + 4 * num_sensors;
	}
	public int getSize() {
		return 5 + 2* num_sensors + 4 * num_sensors;
	}
	
	
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
	 * @param num_sensors number of sesnors in an array
	 * @param flt float array containing tasks data
	 * @param indx task number to use
	 * @param use_aux (always false now)
	 */
	public TpTask(int num_sensors, float [] flt, int task_indx, boolean use_aux)
	{
		int indx = task_indx * getSize(num_sensors);
		task =    Float.floatToIntBits(flt[indx++]); // 0
		int txy = Float.floatToIntBits(flt[indx++]); // 1
		ty = txy >> 16;
		tx = txy & 0xffff;
		target_disparity = flt[indx++];              // 2
		centerXY[0] = flt[indx++];                   // 3
		centerXY[1] = flt[indx++];                   // 4
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
	
	public float [] getCenterXY() {
		return centerXY;
	}

	public double [] getDoubleCenterXY() {
		return new double [] {centerXY[0],centerXY[1]};
	}
	
	public void setCenterXY(double [] centerXY) {
		this.centerXY = new float [] {(float) centerXY[0],(float) centerXY[1]};
	}


	// convert this class instance to float array to match layout of the C struct
	public float [] asFloatArray(boolean use_aux) {
		
		float [] flt = new float [getSize()];
		return asFloatArray(flt, 0, use_aux);
	}
	// convert this class instance to float array to match layout of the C struct,
	// fill existing float array from the specified index
	public float [] asFloatArray(float [] flt, int 	task_indx, boolean use_aux) {
		int indx = task_indx * getSize(num_sensors);
		flt[indx++] = Float.intBitsToFloat(task);            // 0
		flt[indx++] = Float.intBitsToFloat(tx + (ty << 16)); // 1
		flt[indx++] = this.target_disparity;                 // 2
		flt[indx++] = centerXY[0];                           // 3
		flt[indx++] = centerXY[1];                           // 4
		
		float [][] offsets = use_aux? this.xy_aux: this.xy;
		for (int i = 0; i < num_sensors; i++) {
			if (offsets != null) {
				flt[indx++] = offsets[i][0];
				flt[indx++] = offsets[i][1];
			} else {
				indx+= 2;
			}
		}
		/*
		for (int i = 0; i < num_sensors; i++) { // actually disp_dist will be initialized by the GPU
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