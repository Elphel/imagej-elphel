package com.elphel.imagej.gpu;

import static jcuda.driver.CUfunction_attribute.CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import static jcuda.driver.JCudaDriver.cuMemAlloc;
import static jcuda.driver.JCudaDriver.cuMemAllocPitch;
import static jcuda.driver.JCudaDriver.cuMemcpy2D;
import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;
import static jcuda.driver.JCudaDriver.cuMemcpyHtoD;
import static jcuda.driver.JCudaDriver.cuModuleGetGlobal;
import static jcuda.driver.JCudaDriver.cuFuncSetAttribute;


import java.awt.Rectangle;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.tileprocessor.CorrVector;
import com.elphel.imagej.tileprocessor.Correlation2d;
import com.elphel.imagej.tileprocessor.DttRad2;
import com.elphel.imagej.tileprocessor.GeometryCorrection;
import com.elphel.imagej.tileprocessor.ImageDtt;
import com.elphel.imagej.tileprocessor.QuadCLT;

import Jama.Matrix;
import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUDA_MEMCPY2D;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUmemorytype;

public class GpuQuad{ // quad camera description
	/**
	 * 
	 */
	private final GPUTileProcessor gpuTileProcessor;
	public           QuadCLT quadCLT;
	public final int img_width;
	public final int img_height;
	public final int kernels_hor;        // int                kernels_hor,
	public final int kernels_vert;       // int                kernels_vert);

	public final int num_cams;
	public final int num_colors; // maybe should always be 3?
	public final int kern_tiles;
	public final int kern_size;

	// public final GPUTileProcessor gPUTileProcessor;
	// CPU arrays of pointers to GPU memory
	// These arrays may go to methods, they are here just to be able to free GPU memory if needed
	private CUdeviceptr [] gpu_kernels_h;
	private CUdeviceptr [] gpu_kernel_offsets_h;
	private CUdeviceptr [] gpu_bayer_h;
	private CUdeviceptr [] gpu_clt_h;
	private CUdeviceptr [] gpu_corr_images_h;
	// GPU pointers to array of GPU pointers
	private CUdeviceptr gpu_kernels;
	private CUdeviceptr gpu_kernel_offsets;
	private CUdeviceptr gpu_bayer;
//	private CUdeviceptr gpu_tasks;
	private CUdeviceptr gpu_ftasks;
	private CUdeviceptr gpu_corrs;
	private CUdeviceptr gpu_corr_weights;
	private CUdeviceptr gpu_corrs_td;
	private CUdeviceptr gpu_corrs_combo;
	private CUdeviceptr gpu_corrs_combo_td;

	private CUdeviceptr gpu_textures;
	private CUdeviceptr gpu_clt;
	private CUdeviceptr gpu_4_images; // may actually be 16
	private CUdeviceptr gpu_corr_indices;
	private CUdeviceptr gpu_corr_combo_indices;
	private CUdeviceptr gpu_num_corr_tiles;
	private CUdeviceptr gpu_texture_indices_ovlp;
	private CUdeviceptr gpu_num_texture_ovlp;
	private CUdeviceptr gpu_texture_indices;
	private CUdeviceptr gpu_texture_indices_len;
	private CUdeviceptr gpu_diff_rgb_combo;
	private CUdeviceptr gpu_color_weights;
	private CUdeviceptr gpu_generate_RBGA_params;
	private CUdeviceptr gpu_woi;
	private CUdeviceptr gpu_num_texture_tiles;
	private CUdeviceptr gpu_textures_rgba;
	private CUdeviceptr gpu_correction_vector;
	private CUdeviceptr gpu_rot_deriv;
	private CUdeviceptr gpu_geometry_correction;
	private CUdeviceptr gpu_rByRDist;
	private CUdeviceptr gpu_active_tiles;
	private CUdeviceptr gpu_num_active_tiles;
	private int mclt_stride;
	private int corr_stride;
	private int corr_stride_td;
	private int corr_stride_combo;
	private int corr_stride_combo_td;
	private int imclt_stride;
	private int texture_stride;
	private int texture_stride_rgba;
	public  int num_task_tiles;
	private int num_corr_tiles;
	private final int num_all_pairs; //  = 6; // number of correlation pairs per tile (should match tasks)
	private int num_corr_combo_tiles;
	private int num_texture_tiles;

	private boolean geometry_correction_set = false;
	private boolean geometry_correction_vector_set = false;
	private int gpu_debug_level = 0; // 1;
	
//    private	boolean []  corr_mask_boolean; // to keep selected correlations between correlation in TD and processing of those correlations
    private	int []  corr_mask_int = new int [4]; // to keep selected correlations between correlation in TD and processing of those correlations
	public boolean [] getCorrMask() {
		boolean [] corr_mask_boolean = new boolean [num_all_pairs];
		for (int i = 0; i < num_all_pairs; i++) {
			corr_mask_boolean[i] = ((corr_mask_int[i>>5] >> (i & 31)) & 1) != 0;
		}
		return corr_mask_boolean;
	}
	public void setCorrMask(int [] mask) {
		corr_mask_int = new int[4];
		for (int i = 0; i < mask.length; i++) {
			corr_mask_int[i] = mask[i];
		}
	}
	public void setCorrMask(boolean [] mask) {
		corr_mask_int = new int [4];
		for (int i = 0; i < mask.length; i++) if (mask[i]){
			corr_mask_int[i >> 5] |= (1 << (i & 31));
		}
	}
	public int getNumUsedPairs() {
		int np = 0;
		for (int i = 0; i < corr_mask_int.length; i++) {
			np += Integer.bitCount(corr_mask_int[i]); 
		}
		return np;
	}
	
	public int [] getUsedPairsList() {
		int [] pairs_list = new int [getNumUsedPairs()];
		boolean [] mask = getCorrMask();
		int indx = 0;
		for (int i =0; i < mask.length; i++) if (mask[i]){
			pairs_list[indx++] = i;
		}
		return pairs_list;
	}
	public int [] getReversedPairsList() {
		int [] pairs_list = getUsedPairsList();  
		int [] reverse_list = new int[num_all_pairs];
		Arrays.fill(reverse_list,-1);
		for (int i = 0; i < pairs_list.length; i++) {
			reverse_list[pairs_list[i]] = i;
		}
		return reverse_list;
	}
	
	
	
	int host_get_textures_shared_size( // in bytes
			//__device__ int get_textures_shared_size( // in bytes
			int                num_cams,     // actual number of cameras
			int                num_colors,   // actual number of colors: 3 for RGB, 1 for LWIR/mono
			int []             offsets){     // in floats
		int DTT_SIZE1 =       GPUTileProcessor.DTT_SIZE + 1;
		int DTT_SIZE2 =       2 * GPUTileProcessor.DTT_SIZE;
		int DTT_SIZE21 =      DTT_SIZE2 + 1;


		int offs = 0;
		if (offsets != null) offsets[0] = offs;
		offs += num_cams * num_colors * 2 * GPUTileProcessor.DTT_SIZE * DTT_SIZE21; //float mclt_tiles         [NUM_CAMS][NUM_COLORS][2*DTT_SIZE][DTT_SIZE21]
		if (offsets != null) offsets[1] = offs;
		offs += num_cams * num_colors * 4 * GPUTileProcessor.DTT_SIZE * DTT_SIZE1;  // float clt_tiles         [NUM_CAMS][NUM_COLORS][4][DTT_SIZE][DTT_SIZE1]
		if (offsets != null) offsets[2] = offs;
		int mclt_tmp_size = num_cams * num_colors * DTT_SIZE2 * DTT_SIZE21;                // [NUM_CAMS][NUM_COLORS][DTT_SIZE2][DTT_SIZE21]
		int rgbaw_size =    (2* (num_colors + 1) + num_cams) * DTT_SIZE2 * DTT_SIZE21;     // [NUM_COLORS + 1 + NUM_CAMS + NUM_COLORS + 1][DTT_SIZE2][DTT_SIZE21]
		offs += (rgbaw_size > mclt_tmp_size) ? rgbaw_size : mclt_tmp_size;
		if (offsets != null) offsets[3] = offs;
		offs += num_cams * 2;                                      // float port_offsets      [NUM_CAMS][2];
		if (offsets != null) offsets[4] = offs;
		offs += num_colors * num_cams;                             // float ports_rgb_shared  [NUM_COLORS][NUM_CAMS];
		if (offsets != null) offsets[5] = offs;
		offs += num_cams;                                          // float max_diff_shared   [NUM_CAMS];
		if (offsets != null) offsets[6] = offs;
		offs += num_cams * GPUTileProcessor.TEXTURE_THREADS_PER_TILE;               // float max_diff_tmp      [NUM_CAMS][TEXTURE_THREADS_PER_TILE]
		if (offsets != null) offsets[7] = offs;
		offs += num_colors * num_cams *  GPUTileProcessor.TEXTURE_THREADS_PER_TILE; //float ports_rgb_tmp     [NUM_COLORS][NUM_CAMS][TEXTURE_THREADS_PER_TILE];
		if (offsets != null) offsets[8] = offs;
		return Sizeof.INT * offs; // shared_floats;
	}
	
	// should only be updated with the same cameras instance
	public void updateQuadCLT(final QuadCLT quadCLT) {
		this.quadCLT =      quadCLT;
		resetGeometryCorrection();
		resetGeometryCorrectionVector();
		this.gpuTileProcessor.bayer_set =        false;
	}
	public void resetBayer() {
		this.gpuTileProcessor.bayer_set =        false;
	}
	public QuadCLT getQuadCLT( ) {
		return this.quadCLT;
	}
	
	public int getTaskSize() {
		return TpTask.getSize(quadCLT.getNumSensors());
	}
	
	public GpuQuad(
			GPUTileProcessor gpuTileProcessor,
			final QuadCLT    quadCLT,
			int              debug_level
			//   			final int img_width,
			//   			final int img_height,
			//   	    	final int kernels_hor,
			//   			final int kernels_vert,
//			final int num_cams,
//			final int num_colors
			) {
		setGpu_debug_level(debug_level);
		this.gpuTileProcessor = gpuTileProcessor;
		this.quadCLT =       quadCLT;
		int [] wh =          quadCLT.getGeometryCorrection().getSensorWH();
		this.img_width =     wh[0]; // img_width;
		this.img_height =    wh[1]; // img_height;
		this.num_cams =      quadCLT.getNumSensors();
		this.num_all_pairs = Correlation2d.getNumPairs(num_cams);
		this.num_colors =    quadCLT.isMonochrome()?1:3; //num_colors; // maybe should always be 3?
		this.kernels_hor =   quadCLT.getCLTKernels()[0][0][0].length; // kernels_hor;
		this.kernels_vert =  quadCLT.getCLTKernels()[0][0].length; // kernels_vert;
		this.kern_tiles =    kernels_hor *  kernels_vert * num_colors;
		this.kern_size =     kern_tiles * 4 * 64;


		// CPU arrays of pointers to GPU memory
		// These arrays may go to methods, they are here just to be able to free GPU memory if needed
		gpu_kernels_h =           new CUdeviceptr[num_cams];
		gpu_kernel_offsets_h =    new CUdeviceptr[num_cams];
		gpu_bayer_h =             new CUdeviceptr[num_cams];
		gpu_clt_h =               new CUdeviceptr[num_cams];
		gpu_corr_images_h=        new CUdeviceptr[num_cams];
		// GPU pointers to array of GPU pointers
		gpu_kernels =             new CUdeviceptr();
		gpu_kernel_offsets =      new CUdeviceptr();
		gpu_bayer =               new CUdeviceptr();
//		gpu_tasks =               new CUdeviceptr(); //  allocate tilesX * tilesY * TPTASK_SIZE * Sizeof.FLOAT
		gpu_ftasks =              new CUdeviceptr(); //  allocate tilesX * tilesY * getTaskSize() * Sizeof.FLOAT
		gpu_corrs =               new CUdeviceptr(); //  allocate tilesX * tilesY * NUM_PAIRS * CORR_SIZE * Sizeof.FLOAT
		gpu_corr_weights =        new CUdeviceptr(); //  allocate tilesX * tilesY * NUM_PAIRS * Sizeof.FLOAT

		gpu_corrs_td =            new CUdeviceptr(); //  allocate tilesX * tilesY * NUM_PAIRS * 4 * DTT_SIZE * DTT_SIZE * Sizeof.FLOAT
		gpu_corrs_combo =         new CUdeviceptr(); //  allocate tilesX * tilesY             * CORR_SIZE * Sizeof.FLOAT
		gpu_corrs_combo_td =      new CUdeviceptr(); //  allocate tilesX * tilesY *             4 * DTT_SIZE * DTT_SIZE * Sizeof.FLOAT

		gpu_textures =            new CUdeviceptr(); //  allocate tilesX * tilesY * ? * 256 * Sizeof.FLOAT
		gpu_clt =                 new CUdeviceptr();
		gpu_4_images =            new CUdeviceptr();
		gpu_corr_indices =        new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
		// May add separate gpu_corr_indices_td here
		gpu_corr_combo_indices =  new CUdeviceptr(); //  allocate tilesX * tilesY * 1 * Sizeof.FLOAT            
		gpu_num_corr_tiles =      new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
		gpu_texture_indices_ovlp =new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
		gpu_num_texture_ovlp =    new CUdeviceptr(); //  8 ints
		gpu_texture_indices =     new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
		gpu_texture_indices_len = new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
		gpu_diff_rgb_combo =      new CUdeviceptr(); //  1 int

		gpu_color_weights =       new CUdeviceptr(); //  allocate 3 * Sizeof.FLOAT
		gpu_generate_RBGA_params =new CUdeviceptr(); //  allocate 5 * Sizeof.FLOAT

		gpu_woi =                 new CUdeviceptr(); //  4 integers (x, y, width, height) Rectangle - in tiles
		gpu_num_texture_tiles =   new CUdeviceptr(); //  8 integers

		gpu_textures_rgba =       new CUdeviceptr(); //  allocate tilesX * tilesY * ? * 256 * Sizeof.FLOAT

		gpu_correction_vector=    new CUdeviceptr();
		gpu_rot_deriv=            new CUdeviceptr(); //  used internally by device, may be read to CPU for testing
		gpu_geometry_correction=  new CUdeviceptr();
		gpu_rByRDist=             new CUdeviceptr(); //  calculated once for the camera distortion model in CPU (move to GPU?)

		gpu_active_tiles =        new CUdeviceptr(); //  TILESX*TILESY*sizeof(int)
		gpu_num_active_tiles =    new CUdeviceptr(); //  1 int

		// Init data arrays for all kernels
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =  img_height / GPUTileProcessor.DTT_SIZE;
		long [] device_stride = new long [1];
		for (int ncam = 0; ncam < num_cams; ncam++) {
			gpu_kernels_h[ncam] =        new CUdeviceptr();
			cuMemAlloc(gpu_kernels_h[ncam],kern_size * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
			gpu_kernel_offsets_h[ncam] = new CUdeviceptr();
			cuMemAlloc(gpu_kernel_offsets_h[ncam],kern_tiles * GPUTileProcessor.CLTEXTRA_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
			gpu_bayer_h[ncam] =          new CUdeviceptr();
			cuMemAllocPitch (
					gpu_bayer_h[ncam],        // CUdeviceptr dptr,
					device_stride,            // long[] pPitch,
					img_width * Sizeof.FLOAT, // long WidthInBytes,
					img_height,               // long Height,
					Sizeof.FLOAT);            // int ElementSizeBytes)
			mclt_stride = (int)(device_stride[0] / Sizeof.FLOAT);
			gpu_corr_images_h[ncam] =  new CUdeviceptr();
			cuMemAllocPitch (
					gpu_corr_images_h[ncam],               // CUdeviceptr dptr,
					device_stride,                         // long[] pPitch,
					(img_width + GPUTileProcessor.DTT_SIZE) * Sizeof.FLOAT, // long WidthInBytes,
					3*(img_height + GPUTileProcessor.DTT_SIZE),// long Height,
					Sizeof.FLOAT);            // int ElementSizeBytes)
			imclt_stride = (int)(device_stride[0] / Sizeof.FLOAT);
			gpu_clt_h[ncam] = new CUdeviceptr();
			cuMemAlloc(gpu_clt_h[ncam],tilesY * tilesX * num_colors * 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
		}
		// now create device arrays pointers
		if (Sizeof.POINTER != Sizeof.LONG) {
			String msg = "Sizeof.POINTER != Sizeof.LONG";
			IJ.showMessage("Error",	msg);
			new IllegalArgumentException (msg);
		}
		cuMemAlloc(gpu_kernels,        num_cams * Sizeof.POINTER);
		cuMemAlloc(gpu_kernel_offsets, num_cams * Sizeof.POINTER);
		cuMemAlloc(gpu_bayer,          num_cams * Sizeof.POINTER);
		cuMemAlloc(gpu_clt,            num_cams * Sizeof.POINTER);
		cuMemAlloc(gpu_4_images,       num_cams * Sizeof.POINTER);

		long [] gpu_kernels_l =        new long [num_cams];
		long [] gpu_kernel_offsets_l = new long [num_cams];
		long [] gpu_bayer_l =          new long [num_cams];
		long [] gpu_clt_l =            new long [num_cams];
		long [] gpu_4_images_l =       new long [num_cams];

		for (int ncam = 0; ncam < num_cams; ncam++) gpu_kernels_l[ncam] =        GPUTileProcessor.getPointerAddress(gpu_kernels_h[ncam]);
		cuMemcpyHtoD(gpu_kernels, Pointer.to(gpu_kernels_l),                     num_cams * Sizeof.POINTER);

		for (int ncam = 0; ncam < num_cams; ncam++) gpu_kernel_offsets_l[ncam] = GPUTileProcessor.getPointerAddress(gpu_kernel_offsets_h[ncam]);
		cuMemcpyHtoD(gpu_kernel_offsets, Pointer.to(gpu_kernel_offsets_l),       num_cams * Sizeof.POINTER);

		for (int ncam = 0; ncam < num_cams; ncam++) gpu_bayer_l[ncam] =          GPUTileProcessor.getPointerAddress(gpu_bayer_h[ncam]);
		cuMemcpyHtoD(gpu_bayer, Pointer.to(gpu_bayer_l),                         num_cams * Sizeof.POINTER);

		for (int ncam = 0; ncam < num_cams; ncam++) gpu_clt_l[ncam] =            GPUTileProcessor.getPointerAddress(gpu_clt_h[ncam]);
		cuMemcpyHtoD(gpu_clt, Pointer.to(gpu_clt_l),                             num_cams * Sizeof.POINTER);

		for (int ncam = 0; ncam < num_cams; ncam++) gpu_4_images_l[ncam] =       GPUTileProcessor.getPointerAddress(gpu_corr_images_h[ncam]);
		cuMemcpyHtoD(gpu_4_images, Pointer.to(gpu_4_images_l),                   num_cams * Sizeof.POINTER);

		// Set GeometryCorrection data
		cuMemAlloc(gpu_geometry_correction,      GeometryCorrection.arrayLength(GPUTileProcessor.MAX_NUM_CAMS) * Sizeof.FLOAT); // always maximal number of cameras (sparse)
		cuMemAlloc(gpu_rByRDist,                 GPUTileProcessor.RBYRDIST_LEN *  Sizeof.FLOAT);
		cuMemAlloc(gpu_rot_deriv,                5*GPUTileProcessor.MAX_NUM_CAMS*3*3 * Sizeof.FLOAT); // always maximal number of cameras (sparse)
		//        	cuMemAlloc(gpu_correction_vector,        CorrVector.LENGTH * Sizeof.FLOAT);
		cuMemAlloc(gpu_correction_vector,        GPUTileProcessor.CORR_VECTOR_MAX_LENGTH * Sizeof.FLOAT); // update CORR_VECTOR_LENGTH to fit 


		// Set task array
//		cuMemAlloc(gpu_tasks,      tilesX * tilesY * GPUTileProcessor.TPTASK_SIZE * Sizeof.FLOAT);
		cuMemAlloc(gpu_ftasks,      tilesX * tilesY * getTaskSize() * Sizeof.FLOAT);
		//=========== Seems that in many places Sizeof.POINTER (==8) is used instead of Sizeof.FLOAT !!! ============
		// Set corrs array
		int num_pairs = Correlation2d.getNumPairs(quadCLT.getNumSensors());
		cuMemAlloc(gpu_corr_indices,         tilesX * tilesY * num_pairs * Sizeof.FLOAT);
		cuMemAlloc(gpu_corr_combo_indices,   tilesX * tilesY *             Sizeof.FLOAT);

		cuMemAlloc(gpu_num_corr_tiles,                    1 * Sizeof.FLOAT);

		//#define TILESYA       ((TILESY +3) & (~3))
		int tilesYa = (tilesY + 3) & ~3;
		cuMemAlloc(gpu_texture_indices,     tilesX * tilesYa * Sizeof.FLOAT); // for non-overlap tiles
		cuMemAlloc(gpu_texture_indices_ovlp,tilesX * tilesYa * Sizeof.FLOAT); // for overlapped tiles

		cuMemAlloc(gpu_diff_rgb_combo, tilesX * tilesYa * num_cams* (num_colors + 1) *  Sizeof.FLOAT);

		cuMemAlloc(gpu_color_weights,             3 * Sizeof.FLOAT);
		cuMemAlloc(gpu_generate_RBGA_params,      5 * Sizeof.FLOAT);


		cuMemAlloc(gpu_woi,                               4 * Sizeof.FLOAT);
		cuMemAlloc(gpu_num_texture_tiles,                 8 * Sizeof.FLOAT);
		cuMemAlloc(gpu_num_texture_ovlp,                  8 * Sizeof.FLOAT);
		cuMemAlloc(gpu_texture_indices_len,               1 * Sizeof.FLOAT);


		cuMemAlloc(gpu_active_tiles,        tilesX * tilesY * Sizeof.FLOAT);
		cuMemAlloc(gpu_num_active_tiles,                  1 * Sizeof.FLOAT);

		cuMemAlloc(gpu_corr_weights, num_pairs* tilesX * tilesY * Sizeof.FLOAT);

		// allocate space for pixel-domain correlations (6 per tile)
		cuMemAllocPitch (
				gpu_corrs,                             // CUdeviceptr dptr,
				device_stride,                         // long[] pPitch,
				GPUTileProcessor.CORR_SIZE * Sizeof.FLOAT,              // long WidthInBytes,
				num_pairs * tilesX * tilesY,           // long Height,
				Sizeof.FLOAT);                         // int ElementSizeBytes)
		corr_stride = (int)(device_stride[0] / Sizeof.FLOAT);

		// allocate space for transform-domain correlations (6 per tile)
		cuMemAllocPitch (
				gpu_corrs_td,                          // CUdeviceptr dptr,
				device_stride,                         // long[] pPitch,
				4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE * Sizeof.FLOAT,// long WidthInBytes,
				num_pairs * tilesX * tilesY,           // long Height,
				Sizeof.FLOAT);                         // int ElementSizeBytes)
		corr_stride_td = (int)(device_stride[0] / Sizeof.FLOAT);

		// allocate space for pixel-domain combined correlations (1 per tile)
		cuMemAllocPitch (
				gpu_corrs_combo,                       // CUdeviceptr dptr,
				device_stride,                         // long[] pPitch,
				GPUTileProcessor.CORR_SIZE * Sizeof.FLOAT,              // long WidthInBytes,
				tilesX * tilesY,                       // long Height,
				Sizeof.FLOAT);                         // int ElementSizeBytes)
		corr_stride_combo = (int)(device_stride[0] / Sizeof.FLOAT);

		// allocate space for transform-domain combined correlations (1 per tile)
		cuMemAllocPitch (
				gpu_corrs_combo_td,                          // CUdeviceptr dptr,
				device_stride,                         // long[] pPitch,
				4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE * Sizeof.FLOAT,// long WidthInBytes,
				tilesX * tilesY,           // long Height,
				Sizeof.FLOAT);                         // int ElementSizeBytes)
		corr_stride_combo_td = (int)(device_stride[0] / Sizeof.FLOAT);

		int max_texture_size = (num_colors + 1 + (num_cams + num_colors + 1)) * (2 * GPUTileProcessor.DTT_SIZE)* (2 * GPUTileProcessor.DTT_SIZE);
		cuMemAllocPitch (
				gpu_textures,                             // CUdeviceptr dptr,
				device_stride,                         // long[] pPitch,
				max_texture_size * Sizeof.FLOAT,              // long WidthInBytes,
				tilesX * tilesY,             // long Height,
				Sizeof.FLOAT);                         // int ElementSizeBytes)
		texture_stride = (int)(device_stride[0] / Sizeof.FLOAT);
		int max_rgba_width  =  (tilesX + 1) * GPUTileProcessor.DTT_SIZE;
		int max_rgba_height =  (tilesY + 1) * GPUTileProcessor.DTT_SIZE;
		int max_rbga_slices =  num_colors + 1;
		cuMemAllocPitch (
				gpu_textures_rgba,                     // CUdeviceptr dptr,
				device_stride,                         // long[] pPitch,
				max_rgba_width * Sizeof.FLOAT,         // long WidthInBytes,
				max_rgba_height * max_rbga_slices,     // long Height,
				Sizeof.FLOAT);                         // int ElementSizeBytes)
		texture_stride_rgba = (int)(device_stride[0] / Sizeof.FLOAT);
	}
	public int getTilesX() {
		return img_width / GPUTileProcessor.DTT_SIZE;
	}

	public int getTilesY() {
		return img_height / GPUTileProcessor.DTT_SIZE;
	}

	public void resetGeometryCorrection() {
		geometry_correction_set = false;
		geometry_correction_vector_set = false;
		if (getGpu_debug_level() > -1) {
			System.out.println("======resetGeometryCorrection()");
		}

	}
	public void resetGeometryCorrectionVector() {
		geometry_correction_vector_set = false;
		if (getGpu_debug_level() > -1) {
			System.out.println("======resetGeometryCorrectionVector()");
		}
	}
	public int getImageWidth()  {return this.img_width;}
	public int getImageHeight() {return this.img_height;}
	public int getDttSize()     {return GPUTileProcessor.DTT_SIZE;}
//	public int getNumCams()     {return GPUTileProcessor.NUM_CAMS;}

	public void setGeometryCorrection() { // will reset geometry_correction_set when running GPU kernel
		//    		if (geometry_correction_set) return;
		setGeometryCorrection(
				quadCLT.getGeometryCorrection(),
				false);
	}

	public void setGeometryCorrectionVector() { // will reset geometry_correction_vector_set when running GPU kernel
		setExtrinsicsVector( // expand to maximal numer of sensors
				quadCLT.getGeometryCorrection().expandSensors(GPUTileProcessor.MAX_NUM_CAMS).getCorrVector());
		if (getGpu_debug_level() > -1) {
			System.out.println("======setGeometryCorrectionVector()");
		}
	}

	// Use NUM_CAMS = 16 gc.expandSensors(NUM_CAMS) here
	public void setGeometryCorrection(GeometryCorrection gc,
			boolean use_java_rByRDist) { // false - use newer GPU execCalcReverseDistortions
//		float [] fgc = gc.toFloatArray();
		float [] fgc = gc.expandSensors(GPUTileProcessor.MAX_NUM_CAMS).toFloatArray();
		
		if (use_java_rByRDist) {
			double [] rByRDist = gc.getRByRDist();
			float [] fFByRDist = new float [rByRDist.length];
			for (int i = 0; i < rByRDist.length; i++) {
				fFByRDist[i] = (float) rByRDist[i];
			}
			cuMemcpyHtoD(gpu_rByRDist,            Pointer.to(fFByRDist), fFByRDist.length * Sizeof.FLOAT);
		}
		cuMemcpyHtoD(gpu_geometry_correction, Pointer.to(fgc),       fgc.length * Sizeof.FLOAT);
///	Already allocated !	
///		cuMemAlloc  (gpu_rot_deriv, 5 * GPUTileProcessor.MAX_NUM_CAMS *3 *3 * Sizeof.FLOAT); // NCAM of 3x3 rotation matrices, plus 4 derivative matrices for each camera
		if (getGpu_debug_level() > -1) {
			System.out.println("======setGeometryCorrection()");
		}
	}
	
	float [] getRbyRDist() { // read back to verify calculation
		float [] frByRDist = new float [GPUTileProcessor.RBYRDIST_LEN];
		cuMemcpyDtoH(Pointer.to(frByRDist), gpu_rByRDist, frByRDist.length * Sizeof.FLOAT);
		return frByRDist;
	}
	
	public double maxRbyRDistErr() {
		double [] rByRDist = quadCLT.getGeometryCorrection().getRByRDist();
		float [] frByRDist = getRbyRDist(); // read from GPU
		double max_err=0;
		int worst_index = -1;
		for (int i = 0; (i < rByRDist.length) && (i < frByRDist.length); i++){
			double e = Math.abs(rByRDist[i] - frByRDist[i]);
			if (e > max_err) {
				max_err = e;
				worst_index = i;
			}
		}
		if (max_err > 0.001) {
			System.out.println("Maximal R-by-RDist error = "+ max_err+", at index = "+worst_index);
		}
		return max_err;
	}
	
	
	/**
	 * Copy extrinsic correction vector to the GPU memory
	 * @param cv correction vector
	 */
	public void setExtrinsicsVector(CorrVector cv) { // expanded to maximal number of sensors
		double [] dcv = cv.toFullRollArray();
		float []  fcv = new float [dcv.length];
		for (int i = 0; i < dcv.length; i++) {
			fcv[i] = (float) dcv[i];
		}
		cuMemcpyHtoD(gpu_correction_vector, Pointer.to(fcv), fcv.length * Sizeof.FLOAT);
		if (getGpu_debug_level() > -1) {
			System.out.println("======setExtrinsicsVector()");
		}
	}

	/**
	 * Copy array of CPU-prepared tasks to the GPU memory
	 * @param tile_tasks array of TpTask prepared by the CPU (before geometry correction is applied)
	 * @param use_aux Use second (aux) camera
	 * @param verify correct data
	 */
	public void setTasks(
			TpTask [] tile_tasks,
			boolean use_aux,  // while is it in class member? - just to be able to free
			boolean verify
			)
	{
		if (verify) checkTasks(tile_tasks);
		num_task_tiles = tile_tasks.length;
		int task_size = getTaskSize();
		float [] ftasks = new float [task_size * num_task_tiles];
		for (int i = 0; i < num_task_tiles; i++) {
			tile_tasks[i].task |= 511; // &= 0x1f; //############## Testing
			
			tile_tasks[i].asFloatArray(ftasks, i, use_aux);
		}
		cuMemcpyHtoD(gpu_ftasks, Pointer.to(ftasks), task_size * num_task_tiles * Sizeof.FLOAT);
		if (getGpu_debug_level() > -1) {
			System.out.println("======setTasks()");
		}
	}

	/**
	 * Update tp_tasks from the GPU after derived coordinates are calculated
	 * @param tile_tasks tile_tasks array of TpTask prepared by the CPU, will be updated with
	 * GPU-calculated geometry data
	 * @param use_aux  Use second (aux) camera (not used now)
	 */
	
	public void updateTasks(
			TpTask [] tile_tasks,
			boolean use_aux    // while is it in class member? - just to be able to free
			)
	{
		num_task_tiles = tile_tasks.length;
		int task_size = getTaskSize();
		float [] ftasks = new float [task_size * num_task_tiles];
		cuMemcpyDtoH(Pointer.to(ftasks), gpu_ftasks, task_size * num_task_tiles * Sizeof.FLOAT);
		for (int i = 0; i < num_task_tiles; i++) {
			tile_tasks[i] = new TpTask(quadCLT.getNumSensors(), ftasks, i, use_aux);
		}
		
		if (getGpu_debug_level() > -1) {
			System.out.println("======updateTasks()");
		}
	}
	
	

	void checkTasks (TpTask [] tile_tasks) {
		float min_disp = tile_tasks[0].target_disparity;
		float max_disp = min_disp; 
		int min_ty =  tile_tasks[0].ty;
		int min_tx =  tile_tasks[0].tx;
		int max_ty = min_ty;
		int max_tx = min_tx;
		int   max_task = 0;
		for (int i = 0; i < tile_tasks.length; i++) {
			if (Float.isNaN(tile_tasks[i].target_disparity)) {
				System.out.println(String.format("Disparity=NaN for tY=%d, tX=%d, task = 0x%x", tile_tasks[i].ty, tile_tasks[i].tx, tile_tasks[i].task));
				System.out.println("Setting disparity and tas to 0");
				tile_tasks[i].target_disparity = 0f;
				tile_tasks[i].task = 0;
			}
			if (tile_tasks[i].task != 0) {
				if (tile_tasks[i].target_disparity < min_disp ) min_disp = tile_tasks[i].target_disparity; 
				if (tile_tasks[i].target_disparity > max_disp ) max_disp = tile_tasks[i].target_disparity; 
				if (tile_tasks[i].task > max_task )             max_task = tile_tasks[i].task; 

				if (tile_tasks[i].ty < min_ty ) min_ty = tile_tasks[i].ty; 
				if (tile_tasks[i].ty > max_ty ) max_ty = tile_tasks[i].ty; 

				if (tile_tasks[i].tx < min_tx ) min_tx = tile_tasks[i].tx; 
				if (tile_tasks[i].tx > max_tx ) max_tx = tile_tasks[i].tx; 
			}
		}

		System.out.println(String.format("---Disparity=[%f..%f] task=[0..0x%x] ty=[%d..%d]  tx=[%d..%d] ---",
				min_disp,max_disp, max_task, min_ty, max_ty, min_tx, max_tx));
	}



	public TpTask [] getTasks (boolean use_aux)
	{
		int task_size = TpTask.getSize(quadCLT.getNumSensors());
		float [] ftasks = new float [task_size * num_task_tiles];
		cuMemcpyDtoH(Pointer.to(ftasks), gpu_ftasks, task_size * num_task_tiles * Sizeof.FLOAT);
		TpTask [] tile_tasks = new TpTask[num_task_tiles];
		for (int i = 0; i < num_task_tiles; i++) {
			tile_tasks[i] = new TpTask(quadCLT.getNumSensors(), ftasks, i, use_aux);
		}
		return tile_tasks;
	}


	public int [] getTextureIndices()
	{
		float [] ftexture_indices_len = new float[1];
		cuMemcpyDtoH(Pointer.to(ftexture_indices_len), gpu_texture_indices_len,  1 * Sizeof.FLOAT);
		int num_tiles =      Float.floatToIntBits(ftexture_indices_len[0]);
		float [] ftexture_indices = new float [num_tiles];
		cuMemcpyDtoH(Pointer.to(ftexture_indices), gpu_texture_indices,  num_tiles * Sizeof.FLOAT);
		int [] texture_indices = new int [num_tiles];
		for (int i = 0; i < num_tiles; i++) {
			texture_indices[i] = Float.floatToIntBits(ftexture_indices[i]);
		}
		num_texture_tiles = num_tiles; //TODO: was it missing?
		return texture_indices;
	}


	public void setConvolutionKernel(
			float [] kernel,  // [tileY][tileX][color][..]
			float [] kernel_offsets,
			int ncam) {
		cuMemcpyHtoD(gpu_kernels_h[ncam],        Pointer.to(kernel),         kern_size * Sizeof.FLOAT);
		cuMemcpyHtoD(gpu_kernel_offsets_h[ncam], Pointer.to(kernel_offsets), kern_tiles * GPUTileProcessor.CLTEXTRA_SIZE * Sizeof.FLOAT);
	}

	/**
	 * Set CLT kernels in GPU memory reading them using quadCLT instance        
	 * @param clt_kernels - [num_cameras==4][num_colors][tilesY][tilesX][quadrants==4][elements==8*8==64]
	 * @param force re-set kernels even if they were already set
	 */
	public void setConvolutionKernels(
			boolean force)
	{
		if (this.gpuTileProcessor.kernels_set && ! force) {
			return;
		}
		setConvolutionKernels(quadCLT.getCLTKernels(), true);
	}

	/**
	 * Set CLT kernels in GPU memory        
	 * @param clt_kernels - [num_cameras==4][num_colors][tilesY][tilesX][quadrants==4][elements==8*8==64]
	 * @param force re-set kernels even if they were already set
	 */
	public void setConvolutionKernels(
			double [][][][][][] clt_kernels,
			boolean force)
	{
		boolean transpose = true;
		if (this.gpuTileProcessor.kernels_set && ! force) {
			return;
		}
		int num_kernels =
				clt_kernels[0][0].length * //tilesY
				clt_kernels[0][0][0].length * //tilesX
				clt_kernels[0].length;  //colors
		int kernel_length =	num_kernels * 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;

		float [] fkernel =  new float [kernel_length];
		float [] foffsets = new float [num_kernels * GPUTileProcessor.CLTEXTRA_SIZE];
		for (int ncam = 0; ncam < clt_kernels.length; ncam++) {
			int indx=0;
			for (int ty = 0; ty <  clt_kernels[ncam][0].length; ty++) {
				for (int tx = 0; tx <  clt_kernels[ncam][0][ty].length; tx++) {
					for (int col = 0; col <  clt_kernels[ncam].length; col++) {
						for (int p = 0; p < 4; p++) {
							double [] pa = clt_kernels[ncam][col][ty][tx][p];
							for (int i0 = 0; i0 < 64; i0++) {
								int i;
								if (transpose) {
									i = ((i0 & 7) << 3) + ((i0 >>3) & 7);
								} else {
									i = i0;
								}
								fkernel[indx++] = (float)pa[i];
							}
						}
					}
				}
			}
			indx = 0;
			for (int ty = 0; ty <  clt_kernels[ncam][0].length; ty++) {
				for (int tx = 0; tx <  clt_kernels[ncam][0][ty].length; tx++) {
					for (int col = 0; col <  clt_kernels[ncam].length; col++) {
						double [] pa = clt_kernels[ncam][col][ty][tx][4];
						for (int i = 0; i < pa.length; i++) {
							foffsets[indx++] = (float)pa[i];
						}
					}
				}
			}

			setConvolutionKernel(
					fkernel,    // float [] kernel,  // [tileY][tileX][color][..]
					foffsets,   // float [] kernel_offsets,
					ncam);      // int ncam)
		}
		this.gpuTileProcessor.kernels_set = true;
		if (getGpu_debug_level() > -1) {
			System.out.println("======setConvolutionKernels()");
		}
	}

	public void setBayerImage(
			float [] bayer_image,
			int ncam) {
		CUDA_MEMCPY2D copyH2D =   new CUDA_MEMCPY2D();
		copyH2D.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyH2D.srcHost =         Pointer.to(bayer_image);
		copyH2D.srcPitch =        img_width*Sizeof.FLOAT; // width_in_bytes;
		copyH2D.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyH2D.dstDevice =       gpu_bayer_h[ncam]; // src_dpointer;
		copyH2D.dstPitch =        mclt_stride *Sizeof.FLOAT; // device_stride[0];
		copyH2D.WidthInBytes =    img_width*Sizeof.FLOAT; // width_in_bytes;
		copyH2D.Height =          img_height; // /4;
		cuMemcpy2D(copyH2D);
	}

	/**
	 * Copy a set of images to the GPU (if they are new)        
	 * @param force set even if there is no new Bayer data available
	 */
	public void setBayerImages(
			boolean                  force) {
		if (!force && this.gpuTileProcessor.bayer_set && !quadCLT.hasNewImageData()) {
			return;
		}
		double [][][]       bayer_data = quadCLT.getImageData(); // resets hasNewImageData()
		setBayerImages(
				bayer_data,
				true);
	}        

	/**
	 * Copy a set of Bayer images to the GPU         
	 * @param bayer_data per camera a set of split-color images [4][3][w*h], may be [4][1][w*h]
	 * @param force copy regardless of if it was copied before
	 */
	// combines 3 bayer channels into one and transfers to GPU memory
	public void setBayerImages(
			double [][][]       bayer_data,
			boolean                  force) {
		if (this.gpuTileProcessor.bayer_set && !force) {
			return;
		}
		float [] fbayer = new float [bayer_data[0][0].length];
		for (int ncam = 0; ncam < bayer_data.length; ncam++) {
			for (int i = 0; i <  bayer_data[ncam][0].length; i++) {
				fbayer[i] = (float) (bayer_data[ncam][0][i]); //  + bayer_data[ncam][1][i] + bayer_data[ncam][2][i]);
				for (int j = 1; j < bayer_data[ncam].length; j++) {
					fbayer[i] += (float) (bayer_data[ncam][j][i]);
				}
			}
			setBayerImage(
					fbayer, // float [] bayer_image,
					ncam); // int ncam)
		}
		this.gpuTileProcessor.bayer_set = true;
		if (getGpu_debug_level() > -1) {
			System.out.println("======setBayerImages()");
		}
	}

	// prepare tasks for full frame, same dispaity.
	// need to run setTasks(TpTask [] tile_tasks, boolean use_aux) to format/transfer to GPU memory
	/**
	 * CPU-side preparation of per-tile tasks , will need GPU kernel run to finalize          
	 * @param woi - rectangle WOI in tiles (x,y,w,h) to process, normally full window (may provide null for full frame, faster)      
	 * @param round_woi - debug feature - use round (elliptic) WOI, ignored for woi == null
	 * @param target_disparity - same disparity for all tiles
	 * @param disparity_corr -   add to all disparities (minus disparity at infinity)
	 * @param out_image from which camera channels to generate image (currently 0/1)
	 * @param corr_mask from which correlation pairs to generate (0x3f - all 6)
	 * @param threadsMax - historic
	 * @param debugLevel -  not yet used
	 * @return per-tile array of TpTask instances (need to be finalized by the GPU kernel) 
	 */
	public TpTask [] setFullFrameImages(
			Rectangle                 woi,
			boolean                   round_woi,
			float                     target_disparity, // apply same disparity to all tiles
			float                     disparity_corr, // add to disparity (at infinity)
			int                       out_image, // from which tiles to generate image (currently 0/1)
			int                       corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
			//    			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel) {
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =  img_height / GPUTileProcessor.DTT_SIZE;
		float [] target_disparities = new float [tilesX * tilesY];
		int [] out_images = new int [tilesX * tilesY];
		int [] corr_masks = new int [tilesX * tilesY];
		if (target_disparity != 0.0) {
			for (int i = 0; i <target_disparities.length; i++ ) target_disparities[i] = target_disparity;
		}
		for (int i = 0; i <out_images.length; i++ ) {
			out_images[i] = out_image; //  0xf;  // all 4 images
			corr_masks[i] = corr_mask; // 0x3f; // all 6 correlations
		}
		return setFullFrameImages(
				woi,                // Rectangle                 woi,
				round_woi,          // boolean                   round_woi,
				target_disparities, // should be tilesX*tilesY long
				disparity_corr,     // float                     disparity_corr, // add to disparity (at infinity)
				out_images,         // int   []                  out_images, // from which tiles to generate image (currently 0/1)
				corr_masks,         // int   []                  corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
				//        			threadsMax,             // maximal number of threads to launch
				debugLevel);
	}

	/**
	 * CPU-side preparation of per-tile tasks, will need GPU kernel run to finalize          
	 * @param woi - rectangle woi in tiles (x,y,w,h) to process, normally full window (may provide null for full frame, faster)      
	 * @param round_woi - debug feature - use round (elliptic) WOI, ignored for woi == null
	 * @param target_disparities 1D array [tilesX*tilesY] of target disparities
	 * @param disparity_corr -   add to all disparities (minus disparity at infinity)
	 * @param out_images per-tile array [tilesX*tilesY] from which tiles to generate image (currently 0/1)
	 * @param corr_mask per-tile array [tilesX*tilesY] from which correlation pairs to generate (0x3f - all 6)
	 * @param threadsMax - historic
	 * @param debugLevel -  not yet used
	 * @return per-tile array of TpTask instances (need to be finalized by the GPU kernel) 
	 */
	public TpTask [] setFullFrameImages(
			Rectangle                 woi, // or null
			boolean                   round_woi,
			float []                  target_disparities, // should be tilesX*tilesY long
			float                     disparity_corr, // add to disparity (at infinity)
			int   []                  out_images, // from which tiles to generate image (currently 0/1)
			int   []                  corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
			final int                 debugLevel)
	{
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =  img_height / GPUTileProcessor.DTT_SIZE;
		boolean [] mask = new boolean[tilesX*tilesY];
		int num_tiles = 0;
		if (woi == null) {
			round_woi = false;
			woi = new Rectangle(0,0,tilesX,tilesY);
			num_tiles = tilesX*tilesY;
			for (int i = 0; i < num_tiles; i++) {
				mask[i] = true;
			}
		} else { // mostly for debug
			if (woi.x < 0) woi.x = 0;
			if (woi.y < 0) woi.y = 0;
			if (woi.x > (tilesX - 2)) woi.x = tilesX - 2;
			if (woi.y > (tilesY - 2)) woi.y = tilesY - 2;
			if ((woi.x + woi.width)  > tilesX) woi.width =  tilesX - woi.x;
			if ((woi.y + woi.height) > tilesY) woi.height = tilesY - woi.y;
			double rx = 0.5*woi.width;
			double ry = 0.5*woi.height;
			double xc = woi.x + rx - 0.5;
			double yc = woi.y + ry - 0.5;
			for (int ty = woi.y; ty < (woi.y +woi.height); ty++) {
				double ry2 = (ty - yc) / ry;
				ry2*=ry2;
				for (int tx = woi.x; tx < (woi.x +woi.width); tx++) {
					double rx2 = (tx - xc) / rx;
					rx2*=rx2;
					if (!round_woi || ((rx2+ry2) < 1.0)) {
						mask[ty * tilesX + tx] = true;
						num_tiles ++;
					}
				}
			}
		}
		TpTask [] tp_tasks = new TpTask[num_tiles];
		int indx = 0;
		for (int ty = 0; ty < tilesY; ty++) {
			for (int tx = 0; tx < tilesX; tx++) {
				int nt = ty * tilesX + tx;
				if (mask[nt]) {
					// Only generate for non-empty tasks, use 1 empty empty as a terminator?
					tp_tasks[indx] = new TpTask(num_cams, tx, ty, target_disparities[indx]+disparity_corr,
							((out_images[nt] & 0x0f) << 0) |
							((corr_mask [nt] & 0x3f) << 4)
							); // task == 1 for now
					indx++;
				}
			}
		}
		return tp_tasks;
	}


	public TpTask [] setFullFrameImages_dbg( // keeping debug features
			boolean                   calc_offsets, // old way, now not needed with GPU calculation
			Rectangle                 woi, // or null
			boolean                   round_woi,
			float []                  target_disparities, // should be tilesX*tilesY long
			int   []                  out_images, // from which tiles to generate image (currently 0/1)
			int   []                  corr_mask,  // which correlation pairs to generate (maybe later - reduce size from 15x15)
			boolean                   use_master,
			boolean                   use_aux,
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)
	{
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =  img_height / GPUTileProcessor.DTT_SIZE;
		if (woi == null) {
			woi = new Rectangle(0,0,tilesX,tilesY);
		}
		if (woi.x < 0) woi.x = 0;
		if (woi.y < 0) woi.y = 0;
		if (woi.x > (tilesX - 2)) woi.x = tilesX - 2;
		if (woi.y > (tilesY - 2)) woi.y = tilesY - 2;

		if ((woi.x + woi.width)  > tilesX) woi.width =  tilesX - woi.x;
		if ((woi.y + woi.height) > tilesY) woi.height = tilesY - woi.y;
		double rx = 0.5*woi.width;
		double ry = 0.5*woi.height;
		double xc = woi.x + rx - 0.5;
		double yc = woi.y + ry - 0.5;
		boolean dbg1 = false; //  true;
		double dbg_frac = 0.0; // 0.25;
		boolean [] mask = new boolean[tilesX*tilesY];
		int num_tiles = 0;
		for (int ty = woi.y; ty < (woi.y +woi.height); ty++) {
			double ry2 = (ty - yc) / ry;
			ry2*=ry2;
			for (int tx = woi.x; tx < (woi.x +woi.width); tx++) {
				double rx2 = (tx - xc) / rx;
				rx2*=rx2;
				if (!round_woi || ((rx2+ry2) < 1.0)) {
					mask[ty * tilesX + tx] = true;
					num_tiles ++;
				}
			}
		}
		if (dbg_frac > 0) {
			Random rnd = new Random(0);
			int num_final = (int) Math.round(num_tiles * (1.0 - dbg_frac));
			while (num_tiles > num_final) {
				int tx = woi.x + rnd.nextInt(woi.width);
				int ty = woi.y + rnd.nextInt(woi.height);
				int indx = ty * tilesX + tx;
				if (mask[indx]) {
					mask[indx] = false;
					num_tiles--;
				}
			}
			// filter out with no neighbors
			for (int indx = 0; indx < mask.length; indx++) if (mask[indx]) {
				int ix = indx % tilesX;
				int iy = indx / tilesX;
				int num_neib = 0;
				if ((ix > 0) && mask[indx-1]) num_neib++;
				if ((ix < (tilesX-1)) && mask[indx+1]) num_neib++;
				if ((iy > 0) && mask[indx-tilesX]) num_neib++;
				if ((iy < (tilesY-1)) && mask[indx+tilesX]) num_neib++;
				if (num_neib == 0) {
					mask[indx] = false;
					num_tiles--;
				}
			}
			//nextInt(int bound)
		}

		if (dbg1) {
			mask[(woi.y+woi.height) * tilesX + (woi.x+woi.width)] = true;
			num_tiles += 1; // 2;

		}

		TpTask [] tp_tasks = new TpTask[num_tiles];

		int indx = 0;
		for (int ty = 0; ty < tilesY; ty++) {
			for (int tx = 0; tx < tilesX; tx++) if (mask[ty * tilesX + tx]) {
				// Only generate for non-empty tasks, use 1 empty empty as a terminator?
				tp_tasks[indx] = new TpTask(num_cams, tx, ty, target_disparities[indx],
						((out_images[indx] & 0x0f) << 0) |
						((corr_mask [indx] & 0x3f) << 4)
						); // task == 1 for now
				indx++;
			}
		}
		if (calc_offsets) {
			getTileSubcamOffsets(
					tp_tasks,                                    // final TpTask[]            tp_tasks,        // will use // modify to have offsets for 8 cameras
					(use_master? geometryCorrection_main: null), // final GeometryCorrection  geometryCorrection_main,
					(use_aux?    geometryCorrection_aux: null),  // final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
					ers_delay,                                   // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
					threadsMax,                                  // final int                 threadsMax,  // maximal number of threads to launch
					debugLevel);                                 // final int                 debugLevel)
		}
		return tp_tasks;
	}

	public TpTask[]  setTpTask(
			//    			final GPUTileProcessor.GpuQuad gpuQuad,
			final double [][]	           disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
			final double                   disparity_corr,
			final boolean []               need_corrs,       // should be initialized to boolean[1] or null
			final int [][]                 tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile
			final int                      corr_mask,        // <0 - use corr mask from the tile tile_op, >=0 - overwrite all with non-zero corr_mask_tp
			final int                      threadsMax)       // maximal number of threads to launch
	{
		final int tilesX = getTilesX();
		final int tilesY = getTilesY();
		final AtomicInteger ai =            new AtomicInteger(0);
		final AtomicBoolean acorrs =        new AtomicBoolean(false);
		final List<TpTask> task_list = new CopyOnWriteArrayList<TpTask>();
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int nTile = ai.getAndIncrement(); nTile < tilesX * tilesY; nTile = ai.getAndIncrement()) {
						int tileY = nTile /tilesX;
						int tileX = nTile % tilesX;
						//						tIndex = tileY * tilesX + tileX;
						if (tile_op[tileY][tileX] == 0) continue; // nothing to do for this tile
						// which images to use
						int img_mask =  ImageDtt.getImgMask(tile_op[tileY][tileX]);
						// which pairs to combine in the combo:  1 - top, 2 bottom, 4 - left, 8 - right
						int corr_mask_tp = ImageDtt.getPairMask(tile_op[tileY][tileX]); // limited to 4 bits only!
						if (corr_mask_tp != 0) {
							if (corr_mask >=0) {
								corr_mask_tp = corr_mask; 	
							}
							if (corr_mask_tp != 0) {
								acorrs.set(true);
							}
						}
						task_list.add(new TpTask(
								num_cams,
								tileX,
								tileY,
								(float) (disparity_array[tileY][tileX] + disparity_corr),
								((img_mask  & 0x0f) << 0) |
								((corr_mask_tp & 0x3f) << 4)
								)); // task == 1 for now
						// mask out pairs that use missing channels
					}

				}
			};
		}
		ImageDtt.startAndJoin(threads);
		if (need_corrs != null) {
			need_corrs[0] = acorrs.get();
		}
		return task_list.toArray(new TpTask[task_list.size()]);		
	}

	/**
	 * Set a single task (partial setup for GPU when geometryCorrection == null)
	 * @param num_cams number of sensors in a system
	 * @param task_code task code (need to be re-defined, currently is mostly zero/non-zero
	 * @param tileX integer tile X coordinate
	 * @param tileY integer tile X coordinate
	 * @param target_disparity target disparity (should include disparity_corr if used)
	 * @param corr_rots Array of rotational matrices, may be null, but it is more efficient to calculate it once. Not used for GPU when geometryCorrection == null 
	 * @param geometryCorrection instance of GeometryCorrection class for the scene. 12/2021 null - rely on GPU for .xy, disp_dist, does not need cor_rots
	 * @return TpTask instance, compatible with both CPU and GPU processing
	 */
	public static TpTask setTask( // sets .xy that can be set inside GPU - fix to remove @Deprecated
			int                 num_cams,
			int                 task_code,
			int                 tileX,
			int                 tileY,
			double              target_disparity, // include disparity_corr
			// use geometryCorrection.getCorrVector().getRotMatrices(); once per whole image set
			Matrix []           corr_rots,        // if null, will be calculated
			GeometryCorrection  geometryCorrection)
	{
		int dbg_tile_x = 32;
		int dbg_tile_y = 88;
		if ((corr_rots == null) && (geometryCorrection != null)) {
			corr_rots = geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		}
		if ((tileX == dbg_tile_x) && (tileY == dbg_tile_y)) {
			System.out.println("Debugging setTask(): tileX="+tileX+", tileY="+tileY);
			System.out.println("Debugging setTask(): tileX="+tileX+", tileY="+tileY);
		}
		TpTask tp_task = new TpTask(num_cams, tileX, tileY);
		tp_task.task = task_code;
		tp_task.target_disparity = (float) target_disparity; // will it be used?
		int transform_size = GPUTileProcessor.DTT_SIZE;
		double	centerX = tileX * transform_size + transform_size/2;
		double	centerY = tileY * transform_size + transform_size/2;
		tp_task.setCenterXY(new double [] {centerX, centerY}); // this pair of coordinates will be used by GPU to set tp_task.xy and task.disp_dist!
		if (geometryCorrection != null) {
			double [][] disp_dist = new double[num_cams][]; // used to correct 3D correlations (not yet used here)
			double [][] centersXY_main = geometryCorrection.getPortsCoordinatesAndDerivatives(
					geometryCorrection, //			GeometryCorrection gc_main,
					false,              // boolean use_rig_offsets,
					corr_rots,          // Matrix []   rots,
					null,               //  Matrix [][] deriv_rots,
					null,               // double [][] pXYderiv,
					disp_dist,          // used to correct 3D correlations
					centerX,
					centerY,
					tp_task.target_disparity); //  + disparity_corr);
			tp_task.setDispDist(disp_dist);
			tp_task.xy = new float [centersXY_main.length][2];
			for (int i = 0; i < centersXY_main.length; i++) {
				tp_task.xy[i][0] = (float) centersXY_main[i][0];
				tp_task.xy[i][1] = (float) centersXY_main[i][1];
			}
		}
		return tp_task;
	}
	
	/**
	 * Set up an array of TpTask. Only partial initialization when GPU is used, provide geometryCorrection == null for that case
	 * @param num_cams number of sensors in a system
	 * @param disparity_array array of target disparities {tiliY][tileX]
	 * @param disparity_corr disparity correction - add to all disparity values
	 * @param tile_op tile operation
	 * @param geometryCorrection geometry correction or null for GPU processing
	 * @param threadsMax maximal number of threads
	 * @return array of prepared TpTask instances
	 */
	public static TpTask[]  setTasks( // sets .xy that can be set inside GPU - fix to remove @Deprecated
			final int                      num_cams,
			final double [][]	           disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
			final double                   disparity_corr,
			final int [][]                 tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile, null - use non-NaN in disparity_array
			final GeometryCorrection       geometryCorrection,
			final int                      threadsMax)       // maximal number of threads to launch
	{
		final int tilesY = disparity_array.length;
		int op = ImageDtt.setImgMask(0, 0xf); // use if tile_op is not provided
		op =     ImageDtt.setPairMask(op,0xf);
		op =     ImageDtt.setForcedDisparity(op,true);
		final int fop = op;
		int tx = -1;
		for (int i = 0; i < disparity_array.length; i++) if (disparity_array[i] != null) {
			tx = disparity_array[i].length;
			break;
		}
		if (tx < 0) {
			return new TpTask[0]; 
		}
		final int tilesX = tx;
		final Matrix [] corr_rots = (geometryCorrection == null) ? null: geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		final TpTask[]  tp_tasks_xy = new TpTask[tilesY * tilesX];
		final AtomicInteger ai =            new AtomicInteger(0);
		final AtomicInteger anum_tiles =    new AtomicInteger(0);
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					for (int iTile = ai.getAndIncrement(); iTile < tp_tasks_xy.length; iTile = ai.getAndIncrement()) {
						int tileY = iTile /tilesX;
						int tileX = iTile % tilesX;
						int op = fop;
						if (tile_op != null) {
							op = tile_op[tileY][tileX];
						}
						if ((op != 0) && !Double.isNaN(disparity_array[tileY][tileX])) {
							tp_tasks_xy[iTile] = setTask(
									num_cams,                                       // int                 num_cams,
//									transform_size,                                 // int                 transform_size,
									op,                                             // int                 task_code,
									tileX,                                          // int                 tileX,
									tileY,                                          // int                 tileY,
									disparity_array[tileY][tileX] + disparity_corr, // double              target_disparity, // include disparity_corr
									corr_rots,                                      // Matrix []           corr_rots,        // if null, will be calculated
									geometryCorrection);                            // GeometryCorrection  geometryCorrection)
							anum_tiles.getAndIncrement();
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		final TpTask[] tp_tasks = new TpTask[anum_tiles.get()];
		int indx = 0;
		for (int i = 0; i < tp_tasks_xy.length; i++) if (tp_tasks_xy[i] != null) {
			tp_tasks[indx++] = tp_tasks_xy[i];
		}
		return tp_tasks;		
	}
	
	
	

	/**
	 * Prepare contents pointers for calculation of the correlation pairs
	 * @param tp_tasks array of tasks that contain masks of the required pairs
	 * @return each element has (tile_number << 8) | (pair_number & 0xff)
	 */
	@Deprecated
	public int [] getCorrTasks(
			TpTask [] tp_tasks) {
		int tilesX = img_width / GPUTileProcessor.DTT_SIZE;
		int num_corr = 0;
		int num_pairs = Correlation2d.getNumPairs(quadCLT.getNumSensors());
		int task_mask = (1 << num_pairs) - 1;
		for (TpTask tt: tp_tasks) {
			int pm = (tt.task >> GPUTileProcessor.TASK_CORR_BITS) & task_mask;
			if (pm != 0) {
				for (int b = 0; b < num_pairs; b++) if ((pm & (1 << b)) != 0) {
					num_corr++;    			}
			}
		}

		int [] iarr = new int[num_corr];
		num_corr = 0;
		for (TpTask tt: tp_tasks) {
			int pm = (tt.task >> GPUTileProcessor.TASK_CORR_BITS) & task_mask;
			if (pm != 0) {
				int tile = (tt.ty * tilesX +tt.tx);
				for (int b = 0; b < num_pairs; b++) if ((pm & (1 << b)) != 0) {
					iarr[num_corr++] = (tile << GPUTileProcessor.CORR_NTILE_SHIFT) | b;
				}
			}
		}
		return iarr;
	}
	/**
	 * Prepare contents pointers for calculation of the texture tiles (RGBA, 16x16)
	 * @param tp_tasks array of tasks that contain masks of the required pairs
	 * @return each element has (tile_number << 8) | (1 << LIST_TEXTURE_BIT)
	 */
	@Deprecated
	public int [] getTextureTasks(
			TpTask [] tp_tasks) {
		int tilesX = img_width / GPUTileProcessor.DTT_SIZE;
		int num_textures = 0;
		for (TpTask tt: tp_tasks) {
			if ((tt.task & GPUTileProcessor.TASK_TEXTURE_BITS) !=0) {
				num_textures++;
			}
		}
		int [] iarr = new int[num_textures];
		num_textures = 0;
		int b = (1 << GPUTileProcessor.LIST_TEXTURE_BIT);
		for (TpTask tt: tp_tasks) {
			if ((tt.task & GPUTileProcessor.TASK_TEXTURE_BITS) !=0) {
				int tile = (tt.ty * tilesX +tt.tx);
				iarr[num_textures++] = (tile << GPUTileProcessor.CORR_NTILE_SHIFT) | b;
			}
		}
		return iarr;
	}

	/**
	 * Calculate rotation matrices and their derivatives
	 * Needed after extrinsics changed
	 * gpu_correction_vector should be set in GPU memory
	 * gpu_rot_deriv will contain results after running
	 */
	public void execRotDerivs() {
		if (geometry_correction_vector_set)	return;
		if (this.gpuTileProcessor.GPU_ROT_DERIV_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_ROT_DERIV_kernel");
			return;
		}
		setGeometryCorrectionVector();
		// kernel parameters: pointer to pointers
		int [] GridFullWarps =    {num_cams, 1, 1}; // round up
		int [] ThreadsFullWarps = {3,        3, 3};
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}), // int                  num_cams,
				Pointer.to(gpu_correction_vector), // struct corr_vector * gpu_correction_vector,
				Pointer.to(gpu_rot_deriv)          // union trot_deriv   * gpu_rot_deriv);
				);

		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_ROT_DERIV_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize(); // remove later
		geometry_correction_vector_set = true;
		if (getGpu_debug_level() > -1) {
			System.out.println("======execRotDerivs()");
		}
	}

	/**
	 *  Calculate reverse radial distortion table from gpu_geometry_correction
	 *  Needed once during initialization of GeometryCorrection
	 *  gpu_geometry_correction should be set in GPU memory
	 *  gpu_rByRDist will contain results after running
	 */
	public void execCalcReverseDistortions() {
		if (geometry_correction_set) return;
		if (this.gpuTileProcessor.GPU_CALC_REVERSE_DISTORTION_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CALC_REVERSE_DISTORTION_kernel");
			return;
		}
		setGeometryCorrection();
		// kernel parameters: pointer to pointers
//		int [] GridFullWarps =    {num_cams, 1, 1}; // round up
		int [] GridFullWarps =    {GPUTileProcessor.MAX_NUM_CAMS, 1, 1};
		int [] ThreadsFullWarps = {3,        3, 3};
		Pointer kernelParameters = Pointer.to(
				Pointer.to(gpu_geometry_correction),     //	struct gc          * gpu_geometry_correction,
				Pointer.to(gpu_rByRDist));                //	float *              gpu_rByRDist)      // length should match RBYRDIST_LEN
		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CALC_REVERSE_DISTORTION_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize(); // remove later
		geometry_correction_set = true;
		if (getGpu_debug_level() > -1) {
			System.out.println("======execCalcReverseDistortions()");
		}
	}

/**
 * Calculate tiles offsets (before each direct conversion run)
 * 
 * @param uniform_grid - true to calculate uniform tile grid, false - use pre-calculated tile centers 
 */
	public void execSetTilesOffsets(boolean uniform_grid) {
		execCalcReverseDistortions(); // will check if it is needed first
		execRotDerivs();              // will check if it is needed first
		if (this.gpuTileProcessor.GPU_CALCULATE_TILES_OFFSETS_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CALCULATE_TILES_OFFSETS_kernel");
			return;
		}
		int [] GridFullWarps =    {1, 1, 1}; 
		int [] ThreadsFullWarps = {1, 1, 1}; // 4,8,1
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { uniform_grid ? 1 : 0 }), // int           uniform_grid, //==0: use provided centers (as for interscene) , !=0 calculate uniform grid
				Pointer.to(new int[] { num_cams}),              // int                  num_cams,
				Pointer.to(gpu_ftasks),                         // float              * gpu_ftasks,
				Pointer.to(new int[] { num_task_tiles }),       // int                  num_tiles,          // number of tiles in task list
				Pointer.to(gpu_geometry_correction),            // struct gc          * gpu_geometry_correction,
				Pointer.to(gpu_correction_vector),              // struct corr_vector * gpu_correction_vector,
				Pointer.to(gpu_rByRDist),                       // float *              gpu_rByRDist)      // length should match RBYRDIST_LEN
				Pointer.to(gpu_rot_deriv));                     // trot_deriv         * gpu_rot_deriv);
		cuCtxSynchronize();
		cuLaunchKernel(this.gpuTileProcessor.GPU_CALCULATE_TILES_OFFSETS_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize(); // remove later
		if (getGpu_debug_level() > -1) {
			System.out.println("======execSetTilesOffsets(), num_task_tiles="+num_task_tiles+", uniform_grid="+uniform_grid);
		}
	}

	/**
	 * Direct CLT conversion and aberration correction 
	 */
	public void execConvertDirect() {
		if (this.gpuTileProcessor.GPU_CONVERT_DIRECT_kernel == null) 
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CONVERT_DIRECT_kernel");
			return;
		}
		setConvolutionKernels(false); // set kernels if they are not set already
		setBayerImages(false); // set Bayer images if this.quadCLT instance has new ones
		// kernel parameters: pointer to pointers
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}),              // int                  num_cams,
				Pointer.to(new int[] { num_colors}),            // int                  num_colors,
				Pointer.to(gpu_kernel_offsets),
				Pointer.to(gpu_kernels),
				Pointer.to(gpu_bayer),
				Pointer.to(gpu_ftasks),
				Pointer.to(gpu_clt),
				Pointer.to(new int[] { mclt_stride }),
				Pointer.to(new int[] { num_task_tiles }),
				// move lpf to 4-image generator kernel - DONE
				Pointer.to(new int[] { 0 }), // lpf_mask
				Pointer.to(new int[] { img_width}),          // int                woi_width,
				Pointer.to(new int[] { img_height}),         // int                woi_height,
				Pointer.to(new int[] { kernels_hor}),        // int                kernels_hor,
				Pointer.to(new int[] { kernels_vert}),       // int                kernels_vert);
				Pointer.to(gpu_active_tiles),
				Pointer.to(gpu_num_active_tiles),
				Pointer.to(new int[] { tilesX })
				);

		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CONVERT_DIRECT_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize(); // remove later
		if (getGpu_debug_level() > -1) {
			System.out.println("======execConvertDirect()");
		}
	}


	/**
	 * Generate corrected image(s) from the CLT representation created by execConvertDirect()
	 * @param is_mono monochrome image (false - RGB)
	 */
	public void execImcltRbgAll(
			boolean is_mono
			) {
		if (this.gpuTileProcessor.GPU_IMCLT_ALL_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_IMCLT_ALL_kernel");
			return;
		}
		int apply_lpf =  1;
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =  img_height / GPUTileProcessor.DTT_SIZE;
		int [] ThreadsFullWarps = {1, 1, 1};
		int [] GridFullWarps =    {1, 1, 1};
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}),                  // int           num_cams,
				Pointer.to(gpu_clt),                                // float      ** gpu_clt, // [num_cams][TILESY][TILESX][num_colors][DTT_SIZE*DTT_SIZE]
				Pointer.to(gpu_4_images),                           // float      ** gpu_corr_images,    // [num_cams][WIDTH, 3 * HEIGHT]
				Pointer.to(new int[] { apply_lpf }),                // int           apply_lpf,
				Pointer.to(new int[] { is_mono ? 1 : num_colors }), // int           colors,
				Pointer.to(new int[] { tilesX }),                   // int           woi_twidth,
				Pointer.to(new int[] { tilesY }),                   // int           woi_theight,
				Pointer.to(new int[] { imclt_stride })              // const size_t  dstride);            // in floats (pixels)
				);
		cuCtxSynchronize();
		
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_IMCLT_ALL_kernel, // failed with LWIR
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                   // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
		if (getGpu_debug_level() > -1) {
			System.out.println("======execImcltRbgAll()");
		}
		
	}

	/**
	 * Generate 2D phase correlations from the CLT representation         
	 * @param scales R,G,B weights
	 * @param fat_zero - absolute fat zero - add to correlations before normalization
	 * @param corr_radius - correlation result size (maximal 7 for 15x15)
	 */
	public void execCorr2D(
			boolean [] pair_select,
			double [] scales,
			double fat_zero,
			int corr_radius) {
		if (this.gpuTileProcessor.GPU_CORRELATE2D_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CORRELATE2D_kernel");
			return;
		}
		double fat_zero2 = fat_zero * fat_zero;
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
//		int num_colors = scales.length;
		int num_colors = this.num_colors;
		if (num_colors > 3) num_colors = 3;
		float fscale0 = (num_colors >1)?((float) scales[0]):1.0f;
		float fscale1 = (num_colors >1)?((float) scales[1]):0.0f;
		float fscale2 = (num_colors >2)?((float) scales[2]):0.0f;
//		int num_pairs = getNumPairs();
		setCorrMask (pair_select);
		
		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}),          // int           num_cams,
				Pointer.to(new int[] { corr_mask_int[0]}),  // int           sel_pairs0, // unused bits should be 0
				Pointer.to(new int[] { corr_mask_int[1]}),  // int           sel_pairs1 // unused bits should be 0
				Pointer.to(new int[] { corr_mask_int[2]}),  // int           sel_pairs2, // unused bits should be 0
				Pointer.to(new int[] { corr_mask_int[3]}),  // int           sel_pairs3, // unused bits should be 0
				Pointer.to(gpu_clt),                        // float          ** gpu_clt,
				Pointer.to(new int[] { num_colors }),       // int               colors,             // number of colors (3/1)
				Pointer.to(new float[] {fscale0  }),        // float             scale0,             // scale for R
				Pointer.to(new float[] {fscale1  }),        // float             scale1,             // scale for B
				Pointer.to(new float[] {fscale2  }),        // float             scale2,             // scale for G
				Pointer.to(new float[] {(float) fat_zero2}),// float             fat_zero,           // here - absolute
				Pointer.to(gpu_ftasks),                     // float           * gpu_tasks,
				Pointer.to(new int[] { num_task_tiles }),   // int               num_tiles           // number of tiles in task
				Pointer.to(new int[] { tilesX }),           // int               tilesx,             // number of tile rows
				Pointer.to(gpu_corr_indices),               // int             * gpu_corr_indices,   // packed tile+pair
				Pointer.to(gpu_num_corr_tiles),             // int             * pnum_corr_tiles,    // pointer to a number of tiles to process
				Pointer.to(new long[] { corr_stride }),     // const size_t      corr_stride,        // in floats
				Pointer.to(new int[] { corr_radius }),      // int               corr_radius,        // radius of the output correlation (7 for 15x15)
				Pointer.to(gpu_corrs)                       // float           * gpu_corrs);         // correlation output data
				);
		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CORRELATE2D_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
	}

	/**
	 * Generate 2D correlations from the CLT representation in transform domain, no normalization         
	 * @param scales R,G,B weights
	 */
	public void execCorr2D_TD(
			double [] scales,
			int    [] sel_pairs_in // [4]
			) {
		if (this.gpuTileProcessor.GPU_CORRELATE2D_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CORRELATE2D_kernel");
			return;
		}
		setCorrMask(sel_pairs_in); // set corr_mask_int
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int num_colors = this.num_colors;
		if (num_colors > 3) num_colors = 3;
		float fscale0 = (num_colors >1)?((float) scales[0]):1.0f;
		float fscale1 = (num_colors >1)?((float) scales[1]):0.0f;
		float fscale2 = (num_colors >2)?((float) scales[2]):0.0f;
		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}),           // int               num_cams,           // number of colors (3/1)
				Pointer.to(new int[] { corr_mask_int[0] }),  // int               sel_pairs0          // unused bits should be 0
				Pointer.to(new int[] { corr_mask_int[1] }),  // int               sel_pairs1          // unused bits should be 0
				Pointer.to(new int[] { corr_mask_int[2] }),  // int               sel_pairs2          // unused bits should be 0
				Pointer.to(new int[] { corr_mask_int[3] }),  // int               sel_pairs3          // unused bits should be 0
				Pointer.to(gpu_clt),                         // float          ** gpu_clt,
				Pointer.to(new int[] { num_colors }),        // int               colors,             // number of colors (3/1)
				Pointer.to(new float[] {fscale0  }),         // float             scale0,             // scale for R
				Pointer.to(new float[] {fscale1  }),         // float             scale1,             // scale for B
				Pointer.to(new float[] {fscale2  }),         // float             scale2,             // scale for G
				Pointer.to(new float[] {(float) 0.0 }),      // float             fat_zero,           // here - absolute does nolt matter if TD (corr_radius==0)
				Pointer.to(gpu_ftasks),                      // float           * gpu_ftasks,
				Pointer.to(new int[] { num_task_tiles }),    // int               num_tiles           // number of tiles in task
				Pointer.to(new int[] { tilesX }),            // int               tilesx,             // number of tile rows
				Pointer.to(gpu_corr_indices),                // int             * gpu_corr_indices,   // packed tile+pair
				Pointer.to(gpu_num_corr_tiles),              // int             * pnum_corr_tiles,    // pointer to a number of tiles to process
				Pointer.to(new long[] { corr_stride_td }),   // const size_t      corr_stride,        // in floats
				Pointer.to(new int[] { 0 }),  // generate TD // int               corr_radius,        // radius of the output correlation (7 for 15x15)
				Pointer.to(gpu_corrs_td)                     // float           * gpu_corrs);         // correlation output data
				);
		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CORRELATE2D_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
		if (getGpu_debug_level() > -1) {
			System.out.println("======execCorr2D_TD()");
		}
		
	}

	/**
	 * Combine intra-scene correlations in transform domain (possible to accumulate more)          
	 * @param init_corr - true: init output to 0 before accumulating, false: add to current value
	 * @param num_pairs_in - typically ==6 - number of pairs per tile (tile task should have same number per each tile).
	 * This number should match correlations in tasks
	 * @param pairs_mask - selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
	 */
	public void execCorr2D_combine(
			boolean init_corr,    // initialize output tiles (false - add to current)
			int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			int     pairs_mask)    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
	{
		execCorr2D_combine (init_corr, num_pairs_in, pairs_mask, false);
	}
	/**
	 * Combine intra-scene correlations in transform domain (possible to accumulate more)          
	 * @param init_corr - true: init output to 0 before accumulating, false: add to current value
	 * @param num_pairs_in - typically ==6 - number of pairs per tile (tile task should have same number per each tile).
	 * This number should match correlations in tasks
	 * @param pairs_mask - selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
	 * @param no_transpose_vertical - do not transpose vertical/diagonal
	 */
	
	@Deprecated // broken
	public void execCorr2D_combine(
			boolean init_corr,    // initialize output tiles (false - add to current)
			int     num_pairs_in, // typically 6 - number of pairs per tile (tile task should have same number per each tile
			int     pairs_mask,    // selected pairs (0x3 - horizontal, 0xc - vertical, 0xf - quad, 0x30 - cross)
			boolean no_transpose_vertical        		
			) {
		if (this.gpuTileProcessor.GPU_CORR2D_COMBINE_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CORR2D_COMBINE_kernel");
			return;
		}
        int num_pairs = getNumUsedPairs();// Number of used pairs		num_pairs = num_pairs_in;
		float [] fnum_corrs = new float[1];
		cuMemcpyDtoH(Pointer.to(fnum_corrs), gpu_num_corr_tiles,  1 * Sizeof.FLOAT);
		num_corr_combo_tiles =      Float.floatToIntBits(fnum_corrs[0])/num_pairs; // number of correlation tiles calculated

		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};
		int init_corr_combo = (init_corr ? 1 : 0) + (no_transpose_vertical ? 2 : 0);
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_corr_combo_tiles }), // num_task_tiles }), // int   num_tiles           // number of tiles in task
				Pointer.to(new int[] { num_pairs }),            // int               num_pairs,          // num pairs per tile (should be the same)
				Pointer.to(new int[] { init_corr_combo }),      // int               init_output,        // 1- reset output tiles to zero before accumulating
				Pointer.to(new int[] { pairs_mask }),           // int               pairs_mask,         // selected pairs
				Pointer.to(gpu_corr_indices),                   // int             * gpu_corr_indices,   // packed tile+pair
				Pointer.to(gpu_corr_combo_indices),             // int             * gpu_combo_indices,  // output if not null: packed tile+pairs_mask
				Pointer.to(new long[] { corr_stride_td }),       // const size_t      corr_stride,        // (in floats) stride for the input TD correlations
				Pointer.to(gpu_corrs_td),                       // float           * gpu_corrs,          // input correlation tiles
				Pointer.to(new long[] { corr_stride_combo_td }), // const size_t      corr_stride_combo,  // (in floats) stride for the output TD
				Pointer.to(gpu_corrs_combo_td));                // float           * gpu_corrs_combo);   // combined correlation output (one per tile)
		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CORR2D_COMBINE_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
	}
	/**
	 * Normalize 2D correlations, transform and unfold 
	 * @param fat_zero - absolute fat zero - add to correlations before normalization
	 * @param corr_radius - correlation result size (maximal 7 for 15x15)
	 */

	public void execCorr2D_normalize(
			boolean combo, // normalize combo correlations (false - per-pair ones) 
			double fat_zero,
//	        int [][][]                num_acc,         // number of accumulated tiles [tilesY][tilesX][pair] (or null)       
  		    float [] fcorr_weights, // null or one per correlation tile (num_corr_tiles) to divide fat zero2
			int corr_radius) {
		if (this.gpuTileProcessor.GPU_CORR2D_NORMALIZE_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_CORR2D_NORMALIZE_kernel");
			return;
		}
		float fat_zero2 = (float) (fat_zero * fat_zero);
		Pointer pFcorr_weights =Pointer.to(new int[]   {0});
		int num_ctiles = combo? num_corr_combo_tiles : num_corr_tiles;
		if (fcorr_weights != null) {
			cuMemcpyHtoD(gpu_corr_weights, Pointer.to(fcorr_weights),  num_ctiles * Sizeof.FLOAT);
			pFcorr_weights = Pointer.to(gpu_corr_weights);
		}
		
		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};
		Pointer kernelParameters;
		if (combo) {
			kernelParameters = Pointer.to(
					Pointer.to(new int[] { num_corr_combo_tiles }), // num_task_tiles }), // int   num_corr_tiles,     // number of correlation tiles to process
					Pointer.to(new long[] { corr_stride_combo_td }),// const size_t      corr_stride_td,     // in floats
					Pointer.to(gpu_corrs_combo_td),                 // float           * gpu_corrs_combo);   // combined correlation output (one per tile)
					pFcorr_weights, 
					Pointer.to(new long[] { corr_stride_combo }),   // const size_t      corr_stride,        // in floats
					Pointer.to(gpu_corrs_combo),                    // float           * gpu_corrs,          // correlation output data (pixel domain)
					Pointer.to(new float[] {fat_zero2}),            // float             fat_zero2,           // here - absolute
					Pointer.to(new int[] { corr_radius }));         // int               corr_radius,        // radius of the output correlation (7 for 15x15)
		} else {
			kernelParameters = Pointer.to(
					Pointer.to(new int[] { num_corr_tiles }),  // num_task_tiles }), // int   num_corr_tiles,     // number of correlation tiles to process
					Pointer.to(new long[] { corr_stride_td }), // const size_t      corr_stride_td,     // in floats
					Pointer.to(gpu_corrs_td), 				   // float           * gpu_corrs_combo);   // combined correlation output (one per tile)
					pFcorr_weights, 
					Pointer.to(new long[] { corr_stride }),    // const size_t      corr_stride,        // in floats
					Pointer.to(gpu_corrs),                     // float           * gpu_corrs,          // correlation output data (pixel domain)
					Pointer.to(new float[] {fat_zero2}),       // float             fat_zero2,           // here - absolute
					Pointer.to(new int[] { corr_radius }));    // int               corr_radius,        // radius of the output correlation (7 for 15x15)
		}

		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CORR2D_NORMALIZE_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
	}

	public void execRBGA(
			double [] color_weights,
			boolean   is_lwir,
			double    min_shot,           // 10.0
			double    scale_shot,         // 3.0
			double    diff_sigma,         // pixel value/pixel change
			double    diff_threshold,     // pixel value/pixel change
			double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean   dust_remove) {
		if (GPUTileProcessor.USE_DS_DP) {
			execRBGA_DP(
					color_weights,  // double [] color_weights,
					is_lwir,        // boolean   is_lwir,
					min_shot,       // double    min_shot,           // 10.0
					scale_shot,     // double    scale_shot,         // 3.0
					diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
					diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove);   // boolean   dust_remove);
		} else {
			execRBGA_noDP(
					color_weights,  // double [] color_weights,
					is_lwir,        // boolean   is_lwir,
					min_shot,       // double    min_shot,           // 10.0
					scale_shot,     // double    scale_shot,         // 3.0
					diff_sigma,     // double    diff_sigma,         // pixel value/pixel change
					diff_threshold, // double    diff_threshold,     // pixel value/pixel change
					min_agree,      // double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove);   // boolean   dust_remove);
		}
	}
	/**
	 * Generate combined (overlapping) texture
	 * @param color_weights - [3] (RGB) or [1] (mono) color weights for matching
	 * @param is_lwir ignore shot noise normalization
	 * @param min_shot - minimal value to use sqrt() for shot noise normalization (default - 10)
	 * @param scale_shot scale shot noise (3.0)
	 * @param diff_sigma pixel value/pixel change (1.5) - normalize channel values difference by the offsets of the cameras 
	 * @param diff_threshold - never used?
	 * @param min_agree  minimal number of channels to agree on a point (real number to work with fuzzy averages
	 * @param dust_remove do not reduce average weight when only one image differs much from the average
	 */
	public void execRBGA_DP(
			double [] color_weights,
			boolean   is_lwir,
			double    min_shot,           // 10.0
			double    scale_shot,         // 3.0
			double    diff_sigma,         // pixel value/pixel change
			double    diff_threshold,     // pixel value/pixel change
			double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean   dust_remove) {
		execCalcReverseDistortions(); // will check if it is needed first
		if (this.gpuTileProcessor.GPU_RBGA_kernel == null) {
			IJ.showMessage("Error", "No GPU kernel: GPU_TEXTURES_kernel");
			return;
		}
//		int num_colors = color_weights.length;
		int num_colors = is_lwir? 1 : color_weights.length;
		if (num_colors > 3) num_colors = 3;
		float [] fcolor_weights = new float[3];
		fcolor_weights[0] = (float) color_weights[0];
		fcolor_weights[1] = (num_colors >1)?((float) color_weights[1]):0.0f;
		fcolor_weights[2] = (num_colors >2)?((float) color_weights[2]):0.0f;
		cuMemcpyHtoD(gpu_color_weights, Pointer.to(fcolor_weights),  fcolor_weights.length * Sizeof.FLOAT);

		float [] generate_RBGA_params = {
				(float)             min_shot,           // 10.0
				(float)             scale_shot,         // 3.0
				(float)             diff_sigma,         // pixel value/pixel change
				(float)             diff_threshold,     // pixel value/pixel change
				(float)             min_agree          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
		};
		cuMemcpyHtoD(gpu_generate_RBGA_params, Pointer.to(generate_RBGA_params),  generate_RBGA_params.length * Sizeof.FLOAT);

		int iis_lwir =      (is_lwir)? 1:0;
		int idust_remove =  (dust_remove)? 1 : 0;

		// uses dynamic parallelization, top kernel is a single-thread one
		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};

	    cuFuncSetAttribute(this.gpuTileProcessor.GPU_RBGA_kernel, CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, 65536);
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}),               // int               num_cams,
				Pointer.to(gpu_ftasks),                          // float            * gpu_ftasks,         // flattened tasks, 27 floats for quad EO, 99 floats for LWIR16
				Pointer.to(new int[] { num_task_tiles }),        // int               num_tiles,          // number of tiles in task list
				// declare arrays in device code?
				Pointer.to(gpu_texture_indices_ovlp),            // int             * gpu_texture_indices_ovlp,// packed tile + bits (now only (1 << 7)
				Pointer.to(gpu_num_texture_ovlp),                // int             * num_texture_tiles,  // number of texture tiles to process (8 elements)
				Pointer.to(gpu_woi),                             // int             * woi,                // x,y,width,height of the woi
				// set smaller for LWIR - it is used to reduce work aread
				Pointer.to(new int[] {img_width / GPUTileProcessor.DTT_SIZE}),     // int                width,  // <= TILESX, use for faster processing of LWIR images (should be actual + 1)
				Pointer.to(new int[] {img_height / GPUTileProcessor.DTT_SIZE}),    // int                height); // <= TILESY, use for faster processing of LWIR images
				// Parameters for the texture generation
				Pointer.to(gpu_clt),                             // float          ** gpu_clt,            // [num_cams] ->[TILESY][TILESX][num_colors][DTT_SIZE*DTT_SIZE]
				Pointer.to(gpu_geometry_correction),             //	struct gc          * gpu_geometry_correction,
				Pointer.to(new int[]   {num_colors}),            // int               colors,             // number of colors (3/1)
				Pointer.to(new int[]   {iis_lwir}),              // int               is_lwir,            // do not perform shot correction
				Pointer.to(gpu_generate_RBGA_params),            // float             generate_RBGA_params[5],
				Pointer.to(gpu_color_weights),                   // float             weights[3],         // scale for R,B,G
				Pointer.to(new int[]   { idust_remove }),        // int               dust_remove,        // Do not reduce average weight when only one image differes much from the average
				Pointer.to(new int[]   {0}),                     // int               keep_weights,       // return channel weights after A in RGBA
				Pointer.to(new int[]   { texture_stride_rgba }), // const size_t      texture_rbga_stride,     // in floats
				Pointer.to(gpu_textures_rgba));                   // float           * gpu_texture_tiles)    // (number of colors +1 + ?)*16*16 rgba texture tiles

		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_RBGA_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
	}
	

	/**
	 * Generate combined (overlapping) texture
	 * Version w/o Dynamic Parallelism to pass shared memory size. May be removed if maximal SM issue will be resolved
	 * @param color_weights - [3] (RGB) or [1] (mono) color weights for matching
	 * @param is_lwir ignore shot noise normalization
	 * @param min_shot - minimal value to use sqrt() for shot noise normalization (default - 10)
	 * @param scale_shot scale shot noise (3.0)
	 * @param diff_sigma pixel value/pixel change (1.5) - normalize channel values difference by the offsets of the cameras 
	 * @param diff_threshold - never used?
	 * @param min_agree  minimal number of channels to agree on a point (real number to work with fuzzy averages
	 * @param dust_remove do not reduce average weight when only one image differs much from the average
	 */
	
	public void execRBGA_noDP(
			double [] color_weights,
			boolean   is_lwir,
			double    min_shot,           // 10.0
			double    scale_shot,         // 3.0
			double    diff_sigma,         // pixel value/pixel change
			double    diff_threshold,     // pixel value/pixel change
			double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean   dust_remove) {
		execCalcReverseDistortions(); // will check if it is needed first
		if (    (this.gpuTileProcessor.GPU_CLEAR_TEXTURE_LIST_kernel == null) &&
				(this.gpuTileProcessor.GPU_MARK_TEXTURE_LIST_kernel == null) &&
				(this.gpuTileProcessor.GPU_MARK_TEXTURE_NEIGHBOR_kernel == null) &&
				(this.gpuTileProcessor.GPU_GEN_TEXTURE_LIST_kernel == null) &&
				(this.gpuTileProcessor.GPU_CLEAR_TEXTURE_RBGA_kernel == null) &&
				(this.gpuTileProcessor.GPU_TEXTURES_ACCUMULATE_kernel == null)) {
			IJ.showMessage("Error", "No GPU kernel(s)");
			return;
		}
		boolean DEBUG8A = false; // true; // false;

		int num_colors = is_lwir? 1 : color_weights.length;
		if (num_colors > 3) num_colors = 3;
		float [] fcolor_weights = new float[3];
		fcolor_weights[0] = (float) color_weights[0];
		fcolor_weights[1] = (num_colors >1)?((float) color_weights[1]):0.0f;
		fcolor_weights[2] = (num_colors >2)?((float) color_weights[2]):0.0f;
		double sw = 0.0; for (int i = 0; i < fcolor_weights.length; i++) sw+= fcolor_weights[i];
		for (int i = 0; i < fcolor_weights.length; i++) fcolor_weights[i]/=sw;
		cuMemcpyHtoD(gpu_color_weights, Pointer.to(fcolor_weights),  fcolor_weights.length * Sizeof.FLOAT);

		float [] generate_RBGA_params = {
				(float)             min_shot,           // 10.0
				(float)             scale_shot,         // 3.0
				(float)             diff_sigma,         // pixel value/pixel change
				(float)             diff_threshold,     // pixel value/pixel change
				(float)             min_agree          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
		};
		cuMemcpyHtoD(gpu_generate_RBGA_params, Pointer.to(generate_RBGA_params),  generate_RBGA_params.length * Sizeof.FLOAT);

		int iis_lwir =      (is_lwir)? 1:0;
		int idust_remove =  (dust_remove)? 1 : 0;
		int width =  img_width / GPUTileProcessor.DTT_SIZE;
		int height = img_height / GPUTileProcessor.DTT_SIZE;

		// Testing parameters
		TpTask [] test_tasks0 = null;
		float [] test_ftasks0 = null;
		int [] cpu_texture_indices_ovlp = null;		
		if (DEBUG8A) {
			test_tasks0 = new TpTask[width * height];
			test_ftasks0 = new float [test_tasks0.length * getTaskSize()];
			cuMemcpyDtoH(Pointer.to(test_ftasks0),            gpu_ftasks,  test_ftasks0.length * Sizeof.FLOAT);
			for (int i = 0; i < test_tasks0.length; i++) {
				test_tasks0[i] = new TpTask(num_cams, test_ftasks0, i, false);	
			}
		}
		
		//		woi.x =      Float.floatToIntBits(fwoi[0]) * GPUTileProcessor.DTT_SIZE;
		//		woi.y =      Float.floatToIntBits(fwoi[1]) * GPUTileProcessor.DTT_SIZE;
		//		woi.width =  Float.floatToIntBits(fwoi[2]) * GPUTileProcessor.DTT_SIZE;
		//		woi.height = Float.floatToIntBits(fwoi[3]) * GPUTileProcessor.DTT_SIZE;

		int blocks_x = (width + ((1 << GPUTileProcessor.THREADS_DYNAMIC_BITS) - 1)) >> GPUTileProcessor.THREADS_DYNAMIC_BITS;
		int []  blocks0 =  {blocks_x, height, 1};
		int []  threads0 = {(1 << GPUTileProcessor.THREADS_DYNAMIC_BITS), 1, 1};
		/*
	    clear_texture_list<<<blocks0,threads0>>>(
	    		gpu_texture_indices,
				width,
				height);
		 */
		Pointer kp_clear_texture_list = Pointer.to(
				Pointer.to(gpu_texture_indices_ovlp), // gpu_texture_indices,
				Pointer.to(new int[] {width}),
				Pointer.to(new int[] {height}));
		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_CLEAR_TEXTURE_LIST_kernel,
				blocks0[0],  blocks0[1],  blocks0[2], // Grid dimension
				threads0[0], threads0[1], threads0[2],// Block dimension
				0, null,                              // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_clear_texture_list, null);         // Kernel- and extra parameters
		cuCtxSynchronize();

		int blocks_t =    (num_task_tiles + ((1 << GPUTileProcessor.THREADS_DYNAMIC_BITS)) -1) >> GPUTileProcessor.THREADS_DYNAMIC_BITS;//
		int []  blocks =  {blocks_t, 1, 1};
		int []  threads = {(1 << GPUTileProcessor.THREADS_DYNAMIC_BITS), 1, 1};
		// mark used tiles in gpu_texture_indices memory
		/*
	    mark_texture_tiles <<<blocks,threads>>>(
	    		num_cams,           // int                num_cams,
				gpu_ftasks,         // float            * gpu_ftasks,         // flattened tasks, 27 floats for quad EO, 99 floats for LWIR16
				num_tiles,          // number of tiles in task list
				width,              // number of tiles in a row
				gpu_texture_indices); // packed tile + bits (now only (1 << 7)
		 */
		Pointer kp_mark_texture_tiles = Pointer.to(
				Pointer.to(new int[] {num_cams}),         // int                num_cams,
				Pointer.to(gpu_ftasks),                   // float            * gpu_ftasks,         // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
				Pointer.to(new int[] { num_task_tiles }), // int                num_tiles,          // number of tiles in task list
				Pointer.to(new int[] {width}),
				Pointer.to(gpu_texture_indices_ovlp));    // gpu_texture_indices,
		cuLaunchKernel(this.gpuTileProcessor.GPU_MARK_TEXTURE_LIST_kernel,
				blocks[0],  blocks[1],  blocks[2],        // Grid dimension
				threads[0], threads[1], threads[2],       // Block dimension
				0, null,                                  // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_mark_texture_tiles, null);             // Kernel- and extra parameters
		cuCtxSynchronize();
		int [] cpu_woi =               new int[4];
		int [] cpu_num_texture_tiles = new int[8];

		//		cuMemcpyDtoH(Pointer.to(cpu_woi), gpu_woi,  cpu_woi.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed
		cpu_woi[0] = width;
		cpu_woi[1] = height;
		cpu_woi[2] = 0;
		cpu_woi[3] = 0;
		cuMemcpyHtoD(gpu_woi, Pointer.to(cpu_woi),  cpu_woi.length * Sizeof.INT);
//		cuMemcpyDtoH(Pointer.to(cpu_woi), gpu_woi,  cpu_woi.length * Sizeof.INT); //just for testing		
		
		/*
    // TODO: create gpu_woi to pass	(copy from woi)
    // set lower 4 bits in each gpu_ftasks task
    mark_texture_neighbor_tiles <<<blocks,threads>>>(
    		num_cams,           // int                num_cams,
			gpu_ftasks,         // float            * gpu_ftasks,         // flattened tasks, 27 floats for quad EO, 99 floats for LWIR16
			num_tiles,           // number of tiles in task list
			width,               // number of tiles in a row
			height,              // number of tiles rows
			gpu_texture_indices, // packed tile + bits (now only (1 << 7)
			gpu_woi);                // min_x, min_y, max_x, max_y

	checkCudaErrors(cudaDeviceSynchronize());
		 */
		Pointer kp_mark_texture_neighbor_tiles = Pointer.to(
				Pointer.to(new int[] {num_cams}),         // int      num_cams,
				Pointer.to(gpu_ftasks),                   // float  * gpu_ftasks,          // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
				Pointer.to(new int[] { num_task_tiles }), // int      num_tiles,           // number of tiles in task list
				Pointer.to(new int[] {width}),            // int      width,               // number of tiles in a row
				Pointer.to(new int[] {height}),           // int      height,              // number of tiles rows
				Pointer.to(gpu_texture_indices_ovlp),     // int    * gpu_texture_indices, // packed tile + bits (now only (1 << 7)
				Pointer.to(gpu_woi));                     // int    * woi,                 // min_x, min_y, max_x, max_y
		cuLaunchKernel(this.gpuTileProcessor.GPU_MARK_TEXTURE_NEIGHBOR_kernel,
				blocks[0],  blocks[1],  blocks[2],        // Grid dimension
				threads[0], threads[1], threads[2],       // Block dimension
				0, null,                                  // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_mark_texture_neighbor_tiles, null);    // Kernel- and extra parameters
		cuCtxSynchronize();
		//		int [] cpu_num_texture_tiles = new int[8];
		Arrays.fill(cpu_num_texture_tiles, 0); 
		cuMemcpyHtoD(gpu_num_texture_tiles, Pointer.to(cpu_num_texture_tiles),  cpu_num_texture_tiles.length * Sizeof.INT);
		/*
    gen_texture_list <<<blocks,threads>>>(
    		num_cams,            // int                num_cams,
			gpu_ftasks,          // float            * gpu_ftasks,         // flattened tasks, 27 floats for quad EO, 99 floats for LWIR16
			num_tiles,           // number of tiles in task list
			width,               // number of tiles in a row
			height,              // int                height,               // number of tiles rows
			gpu_texture_indices, // packed tile + bits (now only (1 << 7)
			gpu_num_texture_tiles,   // number of texture tiles to process
			gpu_woi);                // x,y, here woi[2] = max_X, woi[3] - max-Y
	checkCudaErrors(cudaDeviceSynchronize()); // not needed yet, just for testing
		 */
//		cuMemcpyHtoD(gpu_woi, Pointer.to(cpu_woi),  cpu_woi.length * Sizeof.INT); // just for testing
		
		// Testing parameters
		if (DEBUG8A) {
			TpTask [] test_tasks = new TpTask[width * height];
			float [] test_ftasks = new float [test_tasks.length * getTaskSize()];
			cuMemcpyDtoH(Pointer.to(test_ftasks),            gpu_ftasks,  test_ftasks.length * Sizeof.FLOAT);
			for (int i = 0; i < test_tasks.length; i++) {
				test_tasks[i] = new TpTask(num_cams, test_ftasks, i, false);	
			}

			cuMemcpyDtoH(Pointer.to(cpu_num_texture_tiles), gpu_num_texture_tiles,  cpu_num_texture_tiles.length * Sizeof.INT);		
			cuMemcpyDtoH(Pointer.to(cpu_woi),               gpu_woi,  cpu_woi.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed		
			int tilesX = width;
			int tilesY = height;
			int tilesYa = (tilesY + 3) & ~3;
			cpu_texture_indices_ovlp = new int [tilesX * tilesYa];
			cuMemcpyDtoH(Pointer.to(cpu_texture_indices_ovlp),               gpu_texture_indices_ovlp,  cpu_texture_indices_ovlp.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed		
			cuMemcpyHtoD(gpu_texture_indices_ovlp,                           Pointer.to(cpu_texture_indices_ovlp),  cpu_texture_indices_ovlp.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed		
		}		
		Pointer kp_gen_texture_list = Pointer.to(
				Pointer.to(new int[] {num_cams}),         // int      num_cams,
				Pointer.to(gpu_ftasks),                   // float  * gpu_ftasks,          // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
				Pointer.to(new int[] { num_task_tiles }), // int      num_tiles,           // number of tiles in task list
				Pointer.to(new int[] {width}),            // int      width,               // number of tiles in a row
				Pointer.to(new int[] {height}),           // int      height,              // number of tiles rows
				Pointer.to(gpu_texture_indices_ovlp),     // int    * gpu_texture_indices, // packed tile + bits (now only (1 << 7)
				Pointer.to(gpu_num_texture_tiles),        // int    * gpu_num_texture_tiles,   // number of texture tiles to process
				Pointer.to(gpu_woi));                     // int    * woi,                 // min_x, min_y, max_x, max_y
		cuLaunchKernel(this.gpuTileProcessor.GPU_GEN_TEXTURE_LIST_kernel,
				blocks[0],  blocks[1],  blocks[2],        // Grid dimension
				threads[0], threads[1], threads[2],       // Block dimension
				0, null,                                  // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_gen_texture_list, null);    // Kernel- and extra parameters
		cuCtxSynchronize();
		// copy gpu_woi back to host woi
//		float [] fcpu_woi = new float [4];
//		cuMemcpyDtoH(Pointer.to(fcpu_woi), gpu_woi,  4 * Sizeof.FLOAT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed
		
		cuMemcpyDtoH(Pointer.to(cpu_woi), gpu_woi,  cpu_woi.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed		
		cpu_woi[2] += 1 - cpu_woi[0]; // width  (was min and max)
		cpu_woi[3] += 1 - cpu_woi[1]; // height (was min and max)
		// copy host-modified data back to GPU
		cuMemcpyHtoD(gpu_woi, Pointer.to(cpu_woi),  cpu_woi.length * Sizeof.INT);
		// copy gpu_num_texture_tiles back to host cpu_num_texture_tiles for processing by the CPU code
		cuMemcpyDtoH(Pointer.to(cpu_num_texture_tiles), gpu_num_texture_tiles,  cpu_num_texture_tiles.length * Sizeof.INT);		
		// Testing parameters
		if (DEBUG8A) {
			cuMemcpyDtoH(Pointer.to(cpu_texture_indices_ovlp),               gpu_texture_indices_ovlp,  cpu_texture_indices_ovlp.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed
		}
		// Zero output textures. Trim
		// texture_rbga_stride
		int texture_width =        (cpu_woi[2] + 1) * GPUTileProcessor.DTT_SIZE;
		int texture_tiles_height = (cpu_woi[3] + 1) * GPUTileProcessor.DTT_SIZE;
		int texture_slices =       num_colors + 1;

		int blocks_x2 = ((texture_width + 
				((1 << (GPUTileProcessor.THREADS_DYNAMIC_BITS + GPUTileProcessor.DTT_SIZE_LOG2 )) - 1)) >>
				(GPUTileProcessor.THREADS_DYNAMIC_BITS + GPUTileProcessor.DTT_SIZE_LOG2));
		int [] blocks2  = {blocks_x2, texture_tiles_height * texture_slices, 1}; // each thread - 8 vertical
		int [] threads2 = {(1 << GPUTileProcessor.THREADS_DYNAMIC_BITS), 1, 1};
		/*
		 clear_texture_rbga<<<blocks2,threads2>>>( // illegal value error
				 texture_width,
				 texture_tiles_height * texture_slices, // int               texture_slice_height,
				 texture_rbga_stride,                   // const size_t      texture_rbga_stride,     // in floats 8*stride
				 gpu_texture_tiles) ;                   // float           * gpu_texture_tiles);
		 */

		Pointer kp_clear_texture_rbga = Pointer.to(
				Pointer.to(new int[] {texture_width}),           // int      texture_width,
				Pointer.to(new int[] {texture_tiles_height * texture_slices}),// int               texture_slice_height,
				Pointer.to(new int[]   { texture_stride_rgba }), // const size_t      texture_rbga_stride,     // in floats 8*stride
				Pointer.to(gpu_textures_rgba));                   // float           * gpu_texture_tiles)    // (number of colors +1 + ?)*16*16 rgba texture tiles
		cuLaunchKernel(this.gpuTileProcessor.GPU_CLEAR_TEXTURE_RBGA_kernel,
				blocks2[0],  blocks2[1],  blocks2[2],   // Grid dimension
				threads2[0], threads2[1], threads2[2],  // Block dimension
				0, null,                                // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_clear_texture_rbga, null);           // Kernel- and extra parameters
		cuCtxSynchronize();

		// Testing parameters
		if (DEBUG8A) {
			cuMemcpyDtoH(Pointer.to(cpu_texture_indices_ovlp),               gpu_texture_indices_ovlp,  cpu_texture_indices_ovlp.length * Sizeof.INT); // hope that Float.floatToIntBits(fcorr_indices[i]) is not needed
		}
		
		// Run 8 times - first 4 1-tile offsets  inner tiles (w/o verifying margins), then - 4 times with verification and ignoring 4-pixel
		// oversize (border 16x 16 tiles overhang by 4 pixels)
		int tilesya =  ((height +3) & (~3)); //#define TILES-YA       ((TILES-Y +3) & (~3))

		for (int pass = 0; pass < 8; pass++){
			int num_cams_per_thread = GPUTileProcessor.NUM_THREADS / GPUTileProcessor.TEXTURE_THREADS_PER_TILE; // 4 cameras parallel, then repeat
			int[] threads_texture = {GPUTileProcessor.TEXTURE_THREADS_PER_TILE, num_cams_per_thread, 1}; // TEXTURE_TILES_PER_BLOCK, 1);

			int border_tile =  (pass >> 2);
			int ntt = cpu_num_texture_tiles[((pass & 3) << 1) + border_tile];
			if (ntt > 0) {
				int [] grid_texture = {(ntt + GPUTileProcessor.TEXTURE_TILES_PER_BLOCK-1) / GPUTileProcessor.TEXTURE_TILES_PER_BLOCK,1,1}; // TEXTURE_TILES_PER_BLOCK = 1

				int ti_offset = (pass & 3) * (width * (tilesya >> 2)); //  (TILES-X * (TILES-YA >> 2));  // 1/4
				if (border_tile != 0){
					ti_offset += width * (tilesya >> 2) - ntt; // TILES-X * (TILES-YA >> 2) - ntt;
				}
				int shared_size = host_get_textures_shared_size( // in bytes
						num_cams,     // int                num_cams,     // actual number of cameras
						num_colors,   // int                num_colors,   // actual number of colors: 3 for RGB, 1 for LWIR/mono
						null);           // int *              offsets);     // in floats

				if (DEBUG8A) {
					System.out.println("shared_size = "+shared_size+" ntt="+ntt);
					System.out.println(String.format("generate_RBGA() pass= %d, border_tile= %d, ti_offset= %d, ntt=%d",
							pass, border_tile,ti_offset, ntt));
					//			System.out.println(String.format("generate_RBGA() gpu_texture_indices= %p, gpu_texture_indices + ti_offset= %p\n",
					//					(void *) gpu_texture_indices, (void *) (gpu_texture_indices + ti_offset));
					System.out.println(String.format("\ngenerate_RBGA() grid_texture={%d, %d, %d)\n",
							grid_texture[0], grid_texture[1], grid_texture[2]));
					System.out.println(String.format("\ngenerate_RBGA() threads_texture={%d, %d, %d)\n",
							threads_texture[0], threads_texture[1], threads_texture[2]));

					System.out.println ("pass="+pass);
					for (int i = 0; i < 256; i++){
						int indx0 = ti_offset + i;
						if (indx0 < cpu_texture_indices_ovlp.length) {
							int indx = cpu_texture_indices_ovlp[indx0];
							System.out.println(String.format("%05x | %02d: %4x %3d %3d %2x",indx0, i,indx>>8, (indx>>8) / 80, (indx >> 8) % 80, indx&0xff));
						}
					}
					System.out.println ("\n\n");
				}

				cuFuncSetAttribute(this.gpuTileProcessor.GPU_TEXTURES_ACCUMULATE_kernel, CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, 65536);
				Pointer kp_textures_accumulate = Pointer.to(
						Pointer.to(new int[] {num_cams}),                // int         num_cams,
						Pointer.to(gpu_woi),                             // int       * woi,            // min_x, min_y, max_x, max_y
						Pointer.to(gpu_clt),                             // float    ** gpu_clt,        // [num_cams] ->[TILESY][TILESX][num_colors][DTT_SIZE*DTT_SIZE]
						Pointer.to(new int[] {ntt}),                     // size_t      num_texture_tiles,// number of texture tiles to process
						Pointer.to(new int[] {ti_offset}),                     // size_t      num_texture_tiles,// number of texture tiles to process
						Pointer.to(gpu_texture_indices_ovlp),            //  gpu_texture_indices_offset,// add to gpu_texture_indices
						Pointer.to(gpu_geometry_correction),             // struct gc * gpu_geometry_correction,
						Pointer.to(new int[]   {num_colors}),            // int         colors,         // number of colors (3/1)
						Pointer.to(new int[]   {iis_lwir}),              // int         is_lwir,        // do not perform shot correction
						Pointer.to(new float[] {(float) min_shot}),      // float       min_shot,       // 10.0
						Pointer.to(new float[] {(float) scale_shot}),    // float       scale_shot,     // 3.0
						Pointer.to(new float[] {(float) diff_sigma}),    // float       diff_sigma,     // pixel value/pixel change
						Pointer.to(new float[] {(float) diff_threshold}),// float       diff_threshold, // pixel value/pixel change
						Pointer.to(new float[] {(float) min_agree}),     // float       min_agree,      // minimal number of channels to agree on a point (real number to work with fuzzy averages)
						Pointer.to(gpu_color_weights),                   // float       weights[3],     // scale for R,B,G (or {1.0,0.0,0.0}
						Pointer.to(new int[]   {idust_remove}),          // int         dust_remove,        // Do not reduce average weight when only one image differes much from the average
						Pointer.to(new int[]   {0}),                     // int         keep_weights,       // return channel weights after A in RGBA
						// combining both non-overlap and overlap (each calculated if pointer is not null )
						Pointer.to(new int[]   { texture_stride_rgba }), // const size_t      texture_rbga_stride,     // in floats
						Pointer.to(gpu_textures_rgba),                   // float           * gpu_texture_tiles)    // (number of colors +1 + ?)*16*16 rgba texture tiles
						Pointer.to(new int[]   {0}),                               // size_t      texture_stride,     // in floats (now 256*4 = 1024)
						Pointer.to(new int[]   {0}), // gpu_texture_tiles, // float           * gpu_texture_tiles);  // (number of colors +1 + ?)*16*16 rgba texture tiles
						Pointer.to(new int[]   {0}), // 1, // int               linescan_order,     // if !=0 then output gpu_diff_rgb_combo in linescan order, else  - in gpu_texture_indices order
						Pointer.to(new int[]   {0}), //);//gpu_diff_rgb_combo);             // float           * gpu_diff_rgb_combo) // diff[num_cams], R[num_cams], B[num_cams],G[num_cams]
						Pointer.to(new int[]   {width}));

				cuLaunchKernel(this.gpuTileProcessor.GPU_TEXTURES_ACCUMULATE_kernel, //  jcuda.CudaException: CUDA_ERROR_INVALID_VALUE
						grid_texture[0],    grid_texture[1],    grid_texture[2],   // Grid dimension
						threads_texture[0], threads_texture[1], threads_texture[2],  // Block dimension
						shared_size, null,                      // Shared memory size and stream (shared - only dynamic, static is in code)
						kp_textures_accumulate, null);          // Kernel- and extra parameters

				cuCtxSynchronize();
			}
		}

	}

	
	
	/**
	 * Generate non-overlapping (16x16) texture tiles
	 * @param color_weights - [3] (RGB) or [1] (mono) color weights for matching
	 * @param is_lwir ignore shot noise normalization
	 * @param min_shot - minimal value to use sqrt() for shot noise normalization (default - 10)
	 * @param scale_shot scale shot noise (3.0)
	 * @param diff_sigma pixel value/pixel change (1.5) - normalize channel values difference by the offsets of the cameras 
	 * @param diff_threshold - never used?
	 * @param min_agree  minimal number of channels to agree on a point (real number to work with fuzzy averages
	 * @param dust_remove do not reduce average weight when only one image differs much from the average
	 * @param calc_textures calculate textures (false for macro-only output)
	 * @param calc_extra calculate extra output - low-res for macro
	 */
	public void execTextures(
			double [] color_weights,
			boolean   is_lwir,
			double    min_shot,           // 10.0
			double    scale_shot,         // 3.0
			double    diff_sigma,         // pixel value/pixel change
			double    diff_threshold,     // pixel value/pixel change
			double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
			boolean   calc_textures,
			boolean   calc_extra,
			boolean   linescan_order
			) {
		if (GPUTileProcessor.USE_DS_DP) {
			execTextures_DP(
					color_weights,
					is_lwir,
					min_shot,           // 10.0
					scale_shot,         // 3.0
					diff_sigma,         // pixel value/pixel change
					diff_threshold,     // pixel value/pixel change
					min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove,        // Do not reduce average weight when only one image differs much from the average
					calc_textures,
					calc_extra,
					linescan_order);
		} else {
			execTextures_noDP( // CUDA_ERROR_INVALID_VALUE
					color_weights,
					is_lwir,
					min_shot,           // 10.0
					scale_shot,         // 3.0
					diff_sigma,         // pixel value/pixel change
					diff_threshold,     // pixel value/pixel change
					min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
					dust_remove,        // Do not reduce average weight when only one image differs much from the average
					calc_textures,
					calc_extra,
					linescan_order);
		}
	}
	
	public void execTextures_DP(
			double [] color_weights,
			boolean   is_lwir,
			double    min_shot,           // 10.0
			double    scale_shot,         // 3.0
			double    diff_sigma,         // pixel value/pixel change
			double    diff_threshold,     // pixel value/pixel change
			double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
			boolean   calc_textures,
			boolean   calc_extra,
			boolean   linescan_order)
	{
		execCalcReverseDistortions(); // will check if it is needed first
		if (this.gpuTileProcessor.GPU_TEXTURES_kernel == null)
		{
			IJ.showMessage("Error", "No GPU kernel: GPU_TEXTURES_kernel");
			return;
		}
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int num_colors = color_weights.length;
		if (num_colors > 3) num_colors = 3;
		float [] fcolor_weights = new float[3];
		fcolor_weights[0] = (float) color_weights[0];
		fcolor_weights[1] = (num_colors >1)?((float) color_weights[1]):0.0f;
		fcolor_weights[2] = (num_colors >2)?((float) color_weights[2]):0.0f;
		cuMemcpyHtoD(gpu_color_weights, Pointer.to(fcolor_weights),  fcolor_weights.length * Sizeof.FLOAT);

		float [] generate_RBGA_params = {
				(float)             min_shot,           // 10.0
				(float)             scale_shot,         // 3.0
				(float)             diff_sigma,         // pixel value/pixel change
				(float)             diff_threshold,     // pixel value/pixel change
				(float)             min_agree           // minimal number of channels to agree on a point (real number to work with fuzzy averages)
		};
		cuMemcpyHtoD(gpu_generate_RBGA_params, Pointer.to(generate_RBGA_params),  generate_RBGA_params.length * Sizeof.FLOAT);

		int iis_lwir =      (is_lwir)? 1:0;
		int ilinescan_order = linescan_order? 1 : 0;
		int idust_remove =  (dust_remove)? 1 : 0;

		int [] GridFullWarps =    {1, 1, 1};
		int [] ThreadsFullWarps = {1, 1, 1};

		//        	CUdeviceptr gpu_diff_rgb_combo_local = calc_extra ? gpu_diff_rgb_combo : null;
	    cuFuncSetAttribute(this.gpuTileProcessor.GPU_TEXTURES_kernel, CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, 65536);
		Pointer kernelParameters = Pointer.to(
				Pointer.to(new int[] { num_cams}),               // int                num_cams,
				Pointer.to(gpu_ftasks),                          // float            * gpu_ftasks,
				Pointer.to(new int[] { num_task_tiles }),        // int                num_tiles,          // number of tiles in task list
				// declare arrays in device code?
				Pointer.to(gpu_texture_indices),                 // int              * gpu_texture_indices,// packed tile + bits (now only (1 << 7)
				Pointer.to(gpu_texture_indices_len),
				// Parameters for the texture generation
				Pointer.to(gpu_clt),                             // float          ** gpu_clt,            // [num_cams] ->[TILESY][TILESX][num_colors][DTT_SIZE*DTT_SIZE]
				Pointer.to(gpu_geometry_correction),             //	struct gc          * gpu_geometry_correction,
				Pointer.to(new int[] { num_colors }),
				Pointer.to(new int[] { iis_lwir }),
				Pointer.to(gpu_generate_RBGA_params),            // float             generate_RBGA_params[5],
				Pointer.to(gpu_color_weights),                   // float             weights[3],         // scale for R,B,G
				Pointer.to(new int[] { idust_remove }),
				Pointer.to(new int[] {calc_textures? texture_stride : 0}),
				Pointer.to(gpu_textures),
				Pointer.to(new int[]   {ilinescan_order}),       // 1, // int               linescan_order,     // if !=0 then output gpu_diff_rgb_combo in linescan order, else  - in gpu_texture_indices order
				calc_extra ? Pointer.to(gpu_diff_rgb_combo) : Pointer.to(new int[] { 0 }),
				Pointer.to(new int[] { tilesX }));
		cuCtxSynchronize();
		// Call the kernel function
		cuLaunchKernel(this.gpuTileProcessor.GPU_TEXTURES_kernel,
				GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
				ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
				0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
				kernelParameters, null);   // Kernel- and extra parameters
		cuCtxSynchronize();
	}

	public void execTextures_noDP(
			double [] color_weights,
			boolean   is_lwir,
			double    min_shot,           // 10.0
			double    scale_shot,         // 3.0
			double    diff_sigma,         // pixel value/pixel change
			double    diff_threshold,     // pixel value/pixel change
			double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
			boolean   dust_remove,        // Do not reduce average weight when only one image differs much from the average
			boolean   calc_textures,
			boolean   calc_extra,
			boolean   linescan_order)
	{
		if (num_task_tiles == 0) {
			System.out.println("execTextures_noDP(): num_task_tiles=0!");
			return;
		}
		execCalcReverseDistortions(); // will check if it is needed first
		if (    (this.gpuTileProcessor.GPU_CREATE_NONOVERLAP_LIST_kernel == null) &&
				(this.gpuTileProcessor.GPU_TEXTURES_ACCUMULATE_kernel == null)) {
			IJ.showMessage("Error", "No GPU kernel(s)");
			return;
		}
		int keep_texture_weights = 0; // pass as parameter?
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int num_colors = is_lwir? 1 : color_weights.length;
		if (num_colors > 3) num_colors = 3;
		float [] fcolor_weights = new float[3];
		fcolor_weights[0] = (float) color_weights[0];
		fcolor_weights[1] = (num_colors >1)?((float) color_weights[1]):0.0f;
		fcolor_weights[2] = (num_colors >2)?((float) color_weights[2]):0.0f;
		double sw = 0.0; for (int i = 0; i < fcolor_weights.length; i++) sw+= fcolor_weights[i];
		for (int i = 0; i < fcolor_weights.length; i++) fcolor_weights[i]/=sw;
		int iis_lwir =      (is_lwir)? 1:0;
		int ilinescan_order = linescan_order? 1 : 0;
		int idust_remove =  (dust_remove)? 1 : 0;
		
		cuMemcpyHtoD(gpu_color_weights, Pointer.to(fcolor_weights),  fcolor_weights.length * Sizeof.FLOAT);
		
		/*
		dim3 threads0(CONVERT_DIRECT_INDEXING_THREADS, 1, 1);
		dim3 blocks0 ((tp_task_size + CONVERT_DIRECT_INDEXING_THREADS -1) >> CONVERT_DIRECT_INDEXING_THREADS_LOG2,1, 1);
		int  linescan_order = 1; // output low-res in linescan order, 0 - in gpu_texture_indices order
 		printf("threads0=(%d, %d, %d)\n",threads0.x,threads0.y,threads0.z);
 		printf("blocks0=(%d, %d, %d)\n",blocks0.x,blocks0.y,blocks0.z);
 		int   cpu_pnum_texture_tiles = 0;
 	    int * gpu_pnum_texture_tiles;
 	    checkCudaErrors (cudaMalloc((void **)&gpu_pnum_texture_tiles, sizeof(int)));
#define CONVERT_DIRECT_INDEXING_THREADS_LOG2 5
#define CONVERT_DIRECT_INDEXING_THREADS (1 << CONVERT_DIRECT_INDEXING_THREADS_LOG2) // 32
    int tp_task_size =  TILESX * TILESY; // sizeof(ftask_data)/sizeof(float)/task_size; // number of task tiles
		cpu_pnum_texture_tiles = 0;
		checkCudaErrors(cudaMemcpy(
				gpu_pnum_texture_tiles,
				&cpu_pnum_texture_tiles,
				sizeof(int),
				cudaMemcpyHostToDevice));

	 */
		int CONVERT_DIRECT_INDEXING_THREADS_LOG2 = 5;
		int CONVERT_DIRECT_INDEXING_THREADS = (1 << CONVERT_DIRECT_INDEXING_THREADS_LOG2); // 32

		int [] threads0 = {CONVERT_DIRECT_INDEXING_THREADS, 1, 1};
		int [] blocks0 =  {(num_task_tiles + CONVERT_DIRECT_INDEXING_THREADS -1) >> CONVERT_DIRECT_INDEXING_THREADS_LOG2,1, 1};
		int [] cpu_pnum_texture_tiles = {0};
		cuMemcpyHtoD(gpu_texture_indices_len, Pointer.to(cpu_pnum_texture_tiles),  1 * Sizeof.INT);
		/*
		 create_nonoverlap_list<<<blocks0,threads0>>>(
				 num_cams,                // int                num_cams,
				 gpu_ftasks,              // float            * gpu_ftasks,         // flattened tasks, 27 floats for quad EO, 99 floats for LWIR16
				 tp_task_size,            // int                num_tiles,           // number of tiles in task
				 TILESX,                  // int                width,               // number of tiles in a row
				 gpu_texture_indices,     // int *              nonoverlap_list,     // pointer to the calculated number of non-zero tiles
				 gpu_pnum_texture_tiles); // int *              pnonoverlap_length)  //  indices to gpu_tasks  // should be initialized to zero
		 cudaDeviceSynchronize();
		 */
		Pointer kp_create_nonoverlap_list = Pointer.to(
				Pointer.to(new int[] {num_cams}),         // int     num_cams,           // number of sensors
				Pointer.to(gpu_ftasks),                   // float * gpu_ftasks,         // flattened tasks, 29 floats for quad EO, 101 floats for LWIR16
				Pointer.to(new int[] { num_task_tiles }), // int     num_tiles,          // number of tiles in task list
				Pointer.to(new int[] {tilesX}),           // int     width,              // number of tiles in a row
				Pointer.to(gpu_texture_indices),          // int   * nonoverlap_list,    // pointer to the calculated number of non-zero tiles
				Pointer.to(gpu_texture_indices_len));     // int   * pnonoverlap_length) //  indices to gpu_tasks  // should be initialized to zero
		if (getGpu_debug_level() > -1) {
			System.out.println("======execTextures_noDP(), num_task_tiles="+num_task_tiles);
		}
		cuLaunchKernel(this.gpuTileProcessor.GPU_CREATE_NONOVERLAP_LIST_kernel, //CUDA_ERROR_INVALID_VALUE
				blocks0[0],  blocks0[1],  blocks0[2],        // Grid dimension
				threads0[0], threads0[1], threads0[2],       // Block dimension
				0, null,                                  // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_create_nonoverlap_list, null);             // Kernel- and extra parameters
		cuCtxSynchronize();
		/*
	 checkCudaErrors(cudaMemcpy(
			 &cpu_pnum_texture_tiles,
			 gpu_pnum_texture_tiles,
			 sizeof(int),
			 cudaMemcpyDeviceToHost));
*/
		cuMemcpyDtoH(Pointer.to(cpu_pnum_texture_tiles), gpu_texture_indices_len,  1 * Sizeof.FLOAT);
/*
	 printf("cpu_pnum_texture_tiles = %d\n",  cpu_pnum_texture_tiles);
	 int num_cams_per_thread = NUM_THREADS / TEXTURE_THREADS_PER_TILE; // 4 cameras parallel, then repeat
	 dim3 threads_texture1(TEXTURE_THREADS_PER_TILE, num_cams_per_thread, 1); // TEXTURE_TILES_PER_BLOCK, 1);
	 dim3 grid_texture1((cpu_pnum_texture_tiles + TEXTURE_TILES_PER_BLOCK-1) / TEXTURE_TILES_PER_BLOCK,1,1);
 		printf("threads_texture1=(%d, %d, %d)\n",threads_texture1.x,threads_texture1.y,threads_texture1.z);
 		printf("grid_texture1=(%d, %d, %d)\n",grid_texture1.x,grid_texture1.y,grid_texture1.z);

	 int shared_size = host_get_textures_shared_size( // in bytes
			 num_cams,         // int                num_cams,     // actual number of cameras
			 texture_colors,   // int                num_colors,   // actual number of colors: 3 for RGB, 1 for LWIR/mono
			 0);           // int *              offsets);     // in floats
	printf("\n1. shared_size=%d, num_cams=%d, colors=%d\n",shared_size,num_cams, texture_colors);

	cudaFuncSetAttribute(textures_accumulate, cudaFuncAttributeMaxDynamicSharedMemorySize, 65536); // for CC 7.5
	cudaFuncSetAttribute(textures_accumulate, cudaFuncAttributePreferredSharedMemoryCarveout,cudaSharedmemCarveoutMaxShared);
 */
		int num_cams_per_thread = GPUTileProcessor.NUM_THREADS / GPUTileProcessor.TEXTURE_THREADS_PER_TILE; // 4 cameras parallel, then repeat
		int[] threads_texture = {GPUTileProcessor.TEXTURE_THREADS_PER_TILE, num_cams_per_thread, 1}; // TEXTURE_TILES_PER_BLOCK, 1);
		int[] grid_texture    = {(cpu_pnum_texture_tiles[0] + GPUTileProcessor.TEXTURE_TILES_PER_BLOCK-1) / GPUTileProcessor.TEXTURE_TILES_PER_BLOCK,1,1};
		int shared_size = host_get_textures_shared_size( // in bytes
				num_cams,     // int                num_cams,     // actual number of cameras
				num_colors,   // int                num_colors,   // actual number of colors: 3 for RGB, 1 for LWIR/mono
				null);           // int *              offsets);     // in floats
/*
    				textures_accumulate <<<grid_texture1,threads_texture1,  shared_size>>>( // 65536>>>( //
    						num_cams,                        // 	int               num_cams,           // number of cameras used
							(int *) 0,                       // int             * woi,                // x, y, width,height
							gpu_clt,                         // float          ** gpu_clt,            // [num_cams] ->[TILES-Y][TILES-X][colors][DTT_SIZE*DTT_SIZE]
							cpu_pnum_texture_tiles, // *pnum_texture_tiles,             // size_t            num_texture_tiles,  // number of texture tiles to process
							0,                               //                 gpu_texture_indices_offset,// add to gpu_texture_indices
							gpu_texture_indices,             // int             * gpu_texture_indices,// packed tile + bits (now only (1 << 7)
							gpu_geometry_correction,         // struct gc       * gpu_geometry_correction,
							texture_colors,                  // int               colors,             // number of colors (3/1)
							(texture_colors == 1),           // int               is_lwir,            // do not perform shot correction
							generate_RBGA_params[0], // min_shot,      // float             min_shot,           // 10.0
							generate_RBGA_params[1], // scale_shot,    // float             scale_shot,         // 3.0
							generate_RBGA_params[2], // diff_sigma,    // float             diff_sigma,         // pixel value/pixel change
							generate_RBGA_params[3], // diff_threshold,// float             diff_threshold,     // pixel value/pixel change
							generate_RBGA_params[4], // min_agree,     // float             min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
							
							gpu_color_weights,               // float             weights[3],         // scale for R,B,G
							1,                       // dust_remove,                     // int               dust_remove,        // Do not reduce average weight when only one image differs much from the average
							keep_texture_weights, // 0, // 1                               // int               keep_weights,       // return channel weights after A in RGBA (was removed) (should be 0 if gpu_texture_rbg)?
							// combining both non-overlap and overlap (each calculated if pointer is not null )
							0,                               // size_t      texture_rbg_stride, // in floats
							(float *) 0,                     // float           * gpu_texture_rbg,     // (number of colors +1 + ?)*16*16 rgba texture tiles
							dstride_textures /sizeof(float),          // texture_stride,                  // size_t      texture_stride,     // in floats (now 256*4 = 1024)
							gpu_textures, // (float *) 0,             // gpu_texture_tiles,               //(float *)0);// float           * gpu_texture_tiles);  // (number of colors +1 + ?)*16*16 rgba texture tiles
							linescan_order,          // int               linescan_order,     // if !=0 then output gpu_diff_rgb_combo in linescan order, else  - in gpu_texture_indices order
							gpu_diff_rgb_combo, //);             // float           * gpu_diff_rgb_combo) // diff[num_cams], R[num_cams], B[num_cams],G[num_cams]
							TILESX);
 */
		cuFuncSetAttribute(this.gpuTileProcessor.GPU_TEXTURES_ACCUMULATE_kernel, CU_FUNC_ATTRIBUTE_MAX_DYNAMIC_SHARED_SIZE_BYTES, 65536);
		Pointer kp_textures_accumulate = Pointer.to(
				Pointer.to(new int[] {num_cams}),                // int         num_cams,
				Pointer.to(new int[]   {0}),  // Pointer.to(gpu_woi),  // int       * woi,            // min_x, min_y, max_x, max_y
				Pointer.to(gpu_clt),                             // float    ** gpu_clt,        // [num_cams] ->[TILESY][TILESX][num_colors][DTT_SIZE*DTT_SIZE]
				Pointer.to(new int[] {cpu_pnum_texture_tiles[0]}),                     // size_t      num_texture_tiles,// number of texture tiles to process
				Pointer.to(new int[] {0}),                     // size_t      num_texture_tiles,// number of texture tiles to process
				Pointer.to(gpu_texture_indices),          // int   * nonoverlap_list,    // pointer to the calculated number of non-zero tiles
				Pointer.to(gpu_geometry_correction),             // struct gc * gpu_geometry_correction,
				Pointer.to(new int[]   {num_colors}),            // int         colors,         // number of colors (3/1)
				Pointer.to(new int[]   {iis_lwir}),              // int         is_lwir,        // do not perform shot correction
				Pointer.to(new float[] {(float) min_shot}),      // float       min_shot,       // 10.0
				Pointer.to(new float[] {(float) scale_shot}),    // float       scale_shot,     // 3.0
				Pointer.to(new float[] {(float) diff_sigma}),    // float       diff_sigma,     // pixel value/pixel change
				Pointer.to(new float[] {(float) diff_threshold}),// float       diff_threshold, // pixel value/pixel change
				Pointer.to(new float[] {(float) min_agree}),     // float       min_agree,      // minimal number of channels to agree on a point (real number to work with fuzzy averages)
				Pointer.to(gpu_color_weights),                   // float       weights[3],     // scale for R,B,G (or {1.0,0.0,0.0}
				Pointer.to(new int[]   {idust_remove}),          // int         dust_remove,        // Do not reduce average weight when only one image differes much from the average
				Pointer.to(new int[]   {keep_texture_weights}),  // int         keep_weights,       // return channel weights after A in RGBA
				 // combining both non-overlap and overlap (each calculated if pointer is not null )
				Pointer.to(new int[]   {0}),                     // const size_t      texture_rbga_stride,     // in floats
				Pointer.to(new int[]   {0}),                     // float           * gpu_texture_tiles)    // (number of colors +1 + ?)*16*16 rgba texture tiles
				Pointer.to(new int[] {calc_textures? texture_stride : 0}), // size_t   texture_stride,     // in floats (now 256*4 = 1024)
				Pointer.to(gpu_textures), // gpu_texture_tiles, // float           * gpu_texture_tiles);  // (number of colors +1 + ?)*16*16 rgba texture tiles
				Pointer.to(new int[]   {ilinescan_order}), // 1, // int               linescan_order,     // if !=0 then output gpu_diff_rgb_combo in linescan order, else  - in gpu_texture_indices order
				calc_extra ? Pointer.to(gpu_diff_rgb_combo) : Pointer.to(new int[] { 0 }),
				Pointer.to(new int[] { tilesX }));
		cuLaunchKernel(this.gpuTileProcessor.GPU_TEXTURES_ACCUMULATE_kernel,
				grid_texture[0],    grid_texture[1],    grid_texture[2],   // Grid dimension
				threads_texture[0], threads_texture[1], threads_texture[2],  // Block dimension
				shared_size, null,                      // Shared memory size and stream (shared - only dynamic, static is in code)
				kp_textures_accumulate, null);          // Kernel- and extra parameters
		cuCtxSynchronize();
		
	}
	
	
	public int getNumPairs() {return num_all_pairs;}

	
	// 
	/**
	 * Generating correlation sequence by CPU to correlate all tiles provided in linescan order.
	 * Additionally, if both (num_acc != null) and (pfcorr_weights !=null), pfcorr_weights[0]
	 * will return array of the same length as indices_trim (ntiles * num_pairs) to use as weights
	 * to modulate fat_zero2 (the more tiles accumulated - the lower fat_zero2 is)
	 * @param corr_tiles transform domain correlation tile in linescan order [tilesY][tilesX][used_pairs][4*64]
	 * @param num_acc number or tiles accumulated for the same tileX,tileX,pair or null
	 * @param fcorr_weights input - null or [1][], output: null or [1][ntiles * num_pairs] weights
	 * @return [ntiles * num_pairs] correlation instructions:24 MSb - tile number, 8 low bits - pair numbers
	 */
	public int [] setCorrTilesTd(
			final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
            int [][][]           num_acc,     // number of accumulated tiles [tilesY][tilesX][pair] (or null)
            float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2

	{
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
        int num_pairs = getNumUsedPairs();// Number of used pairs		num_pairs = num_pairs_in;
        // 3-rd index of corr_tiles should be getNumUsedPairs();
		int tilesX = corr_tiles[0].length;
		int tilesY = corr_tiles.length;
		int [] indices = new int [tilesY * tilesX * num_pairs];  // as if all tiles are not null
		float [] fdata = new float [tilesY * tilesX * corr_size_td* num_pairs]; // as if all tiles are not null
		boolean use_fcorr_weights = ((num_acc != null) && (pfcorr_weights != null));
		float [] fcorr_weights = use_fcorr_weights ? (new float[tilesY * tilesX * num_pairs]) : null;
		int ntile = 0;
		for (int ty = 0; ty < corr_tiles.length; ty++) {
			for (int tx = 0; tx < corr_tiles[0].length; tx++) {
				if (corr_tiles[ty][tx]!= null) {
					for (int pair = 0; pair < num_pairs; pair++) {
						int corr_pair = ntile * num_pairs + pair;
						indices[corr_pair]= // ntile * num_pairs + pair] =
								((ty * tilesX + tx) << GPUTileProcessor.CORR_NTILE_SHIFT) +
								(pair & ((1 <<  GPUTileProcessor.CORR_NTILE_SHIFT) -1) );
						if (fcorr_weights != null) {
							fcorr_weights[corr_pair] = num_acc[ty][tx][pair];
						}
//						System.arraycopy(corr_tiles[ty][tx][pair], 0, fdata, (ntile*num_pairs + pair) * corr_size_td, corr_size_td);
						System.arraycopy(corr_tiles[ty][tx][pair], 0, fdata, corr_pair * corr_size_td, corr_size_td);
					}
					ntile++;
				}

			}
		}
		setCorrIndicesTdData(
				ntile * num_pairs,   // int    num_tiles,  // corr_indices, fdata may be longer than needed
				indices, // int [] corr_indices,
				fdata);  // float [] fdata);
		// trim indices
		int [] indices_trim = new int [ntile * num_pairs];
		System.arraycopy(indices, 0, indices_trim, 0 , indices_trim.length);
		if (fcorr_weights != null) {
			pfcorr_weights[0] = new float [indices_trim.length];
			System.arraycopy(fcorr_weights, 0, pfcorr_weights[0], 0 , indices_trim.length);
		}
		return indices_trim;
	}

	/**
	 * Generating correlation sequence by CPU to correlate all tiles provided in linescan order.
	 * This versions follows order of tiles in tp_task array, filling gaps (should be none) with zeros 
	 * Additionally, if both (num_acc != null) and (pfcorr_weights !=null), pfcorr_weights[0]
	 * will return array of the same length as indices_trim (ntiles * num_pairs) to use as weights
	 * to modulate fat_zero2 (the more tiles accumulated - the lower fat_zero2 is)
	 * @param tp_tasks - array of TpTask, should be updated from GPU
	 * @param corr_tiles transform domain correlation tile in linescan order [tilesY][tilesX][used_pairs][4*64]
	 * @param num_acc number or tiles accumulated for the same tileX,tileX,pair or null
	 * @param fcorr_weights input - null or [1][], output: null or [1][ntiles * num_pairs] weights
	 * @return [ntiles * num_pairs] correlation instructions:24 MSb - tile number, 8 low bits - pair numbers
	 */

	public int [] setCorrTilesTd(
			final TpTask []      tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
            // only listed tiles will be processed
			final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
            float [][][]         num_acc,        // number of accumulated tiles [tilesY][tilesX][pair] (or null)
            float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2
	{
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
        int num_pairs = getNumUsedPairs();// Number of used pairs		num_pairs = num_pairs_in;
        // 3-rd index of corr_tiles should be getNumUsedPairs();
		int [] indices = new int [tp_tasks.length * num_pairs];  // as if all tiles are not null
		float [] fdata = new float [indices.length * corr_size_td]; // as if all tiles are not null
		boolean use_fcorr_weights = ((num_acc != null) && (pfcorr_weights != null));
		float [] fcorr_weights = use_fcorr_weights ? (new float[indices.length]) : null;
		if (fcorr_weights != null) {
			pfcorr_weights[0] = fcorr_weights; // new float [indices.length];
		}
		
		int tilesX = corr_tiles[0].length;
		for (int ntile = 0; ntile < tp_tasks.length; ntile++) {
			int ty = tp_tasks[ntile].ty;
			int tx = tp_tasks[ntile].tx;
			for (int pair = 0; pair < num_pairs; pair++) {
				int corr_pair = ntile * num_pairs + pair;
				indices[corr_pair]= // ntile * num_pairs + pair] =
						((ty * tilesX + tx) << GPUTileProcessor.CORR_NTILE_SHIFT) +
						(pair & ((1 <<  GPUTileProcessor.CORR_NTILE_SHIFT) -1) );
				if (fcorr_weights != null) {
					fcorr_weights[corr_pair] = num_acc[ty][tx][pair];
				}
				if (corr_tiles[ty][tx]!= null) {
					System.arraycopy(corr_tiles[ty][tx][pair], 0, fdata, corr_pair * corr_size_td, corr_size_td);
				} else { // fill with zeros
					Arrays.fill(fdata, corr_pair * corr_size_td, (corr_pair + 1) * corr_size_td, 0.0f);
					if (fcorr_weights != null) {
						fcorr_weights[corr_pair] = 1; 
					}
				}
			}
		}
		setCorrIndicesTdData(
				indices.length, // ntile * num_pairs,   // int    num_tiles,  // corr_indices, fdata may be longer than needed
				indices, // int [] corr_indices,
				fdata);  // float [] fdata);
		return indices;
	}
	

	/**
	 * Generating correlation sequence by CPU to correlate all tiles provided in linescan order.
	 * This versions follows order of tiles in tp_task array, filling gaps (should be none) with zeros 
	 * Additionally, if both (num_acc != null) and (pfcorr_weights !=null), pfcorr_weights[0]
	 * will return array of the same length as indices_trim (ntiles * num_pairs) to use as weights
	 * to modulate fat_zero2 (the more tiles accumulated - the lower fat_zero2 is)
	 * @param tp_tasks - array of TpTask, should be updated from GPU
	 * @param corr_tiles transform domain correlation tile in linescan order [tilesY][tilesX][used_pairs][4*64]
	 * @param dcorr_weight number or tiles accumulated for the index of tp_tasks or null (matching weights for CPU)
	 * @param fcorr_weights input - null or [1][], output: null or [1][ntiles * num_pairs] weights
	 * @return [ntiles * num_pairs] correlation instructions:24 MSb - tile number, 8 low bits - pair numbers
	 */

	public int [] setCorrTilesTd(
			final TpTask []      tp_tasks,        // data from the reference frame - will be applied to LMW for the integrated correlations
            // only listed tiles will be processed
			final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
            double []            dcorr_weight, // num_acc,     // [ntile] (or null)
            float [][]           pfcorr_weights) // null or one per correlation tile (num_corr_tiles) to divide fat zero2
	{
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
        int num_pairs = getNumUsedPairs();// Number of used pairs		num_pairs = num_pairs_in;
        // 3-rd index of corr_tiles should be getNumUsedPairs();
		int [] indices = new int [tp_tasks.length * num_pairs];  // as if all tiles are not null
		float [] fdata = new float [indices.length * corr_size_td]; // as if all tiles are not null
		boolean use_fcorr_weights = ((dcorr_weight != null) && (pfcorr_weights != null));
		float [] fcorr_weights = use_fcorr_weights ? (new float[indices.length]) : null;
		if (fcorr_weights != null) {
			pfcorr_weights[0] = fcorr_weights; // new float [indices.length];
		}
		
		int tilesX = corr_tiles[0].length;
		for (int ntile = 0; ntile < tp_tasks.length; ntile++) {
			int ty = tp_tasks[ntile].ty;
			int tx = tp_tasks[ntile].tx;
			for (int pair = 0; pair < num_pairs; pair++) {
				int corr_pair = ntile * num_pairs + pair;
				indices[corr_pair]= // ntile * num_pairs + pair] =
						((ty * tilesX + tx) << GPUTileProcessor.CORR_NTILE_SHIFT) +
						(pair & ((1 <<  GPUTileProcessor.CORR_NTILE_SHIFT) -1) );
				if (fcorr_weights != null) {
					fcorr_weights[corr_pair] = (float) dcorr_weight[ntile];// num_acc[ty][tx][pair];
				}
				if (corr_tiles[ty][tx]!= null) {
					System.arraycopy(corr_tiles[ty][tx][pair], 0, fdata, corr_pair * corr_size_td, corr_size_td);
				} else { // fill with zeros
					Arrays.fill(fdata, corr_pair * corr_size_td, (corr_pair + 1) * corr_size_td, 0.0f);
					if (fcorr_weights != null) {
						fcorr_weights[corr_pair] = 1; 
					}
				}
			}
		}
		setCorrIndicesTdData(
				indices.length, // ntile * num_pairs,   // int    num_tiles,  // corr_indices, fdata may be longer than needed
				indices, // int [] corr_indices,
				fdata);  // float [] fdata);
		return indices;
	}
	
	
	
	
	
	// Generating correlation sequence by CPU to correlate all tiles provided in linescan order
	@Deprecated
	public int [] setCorrTilesTd(
			final float [][][][] corr_tiles, // [tileY][tileX][pair][4*64]
			int [][] pairs_map) // typically {{0,0},{1,1},{2,2},{3,3},{4,4},{5,5} [0] - 3rd index in corr_tiles, [1] - 
	{
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
//		num_pairs = pairs_map.length; // set global num_pairs
        int num_pairs = getNumUsedPairs();// Number of used pairs		num_pairs = num_pairs_in;
		int tilesX = corr_tiles[0].length;
		int tilesY = corr_tiles.length;
		int [] indices = new int [tilesY * tilesX * num_pairs];  // as if all tiles are not null
		float [] fdata = new float [tilesY * tilesX * corr_size_td* num_pairs]; // as if all tiles are not null
		int ntile = 0;
		for (int ty = 0; ty < corr_tiles.length; ty++) {
			for (int tx = 0; tx < corr_tiles[0].length; tx++) {
				if (corr_tiles[ty][tx]!= null) {
					for (int pair = 0; pair < num_pairs; pair++) {
						indices[ntile * num_pairs + pair] = ((ty * tilesX + tx) << GPUTileProcessor.CORR_NTILE_SHIFT) + pairs_map[pair][1];
						System.arraycopy(corr_tiles[ty][tx][pairs_map[pair][0]], 0, fdata, (ntile*num_pairs + pair) * corr_size_td, corr_size_td);
					}
					ntile++;
				}

			}
		}
		setCorrIndicesTdData(
				ntile * num_pairs,   // int    num_tiles,  // corr_indices, fdata may be longer than needed
				indices, // int [] corr_indices,
				fdata);  // float [] fdata);
		// trim indices
		int [] indices_trim = new int [ntile * num_pairs];
		System.arraycopy(indices, 0, indices_trim, 0 , indices_trim.length);        	
		return indices_trim;
	}

	
	
	
	public float [][][][] getCorrTilesTd() // [tileY][tileX][pair][4*64] , read all available pairs
	{
		int tilesX =     img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =     img_height / GPUTileProcessor.DTT_SIZE;
		float [][][][] corr_tiles = new float[tilesY][tilesX][][]; // num_pairs
		return getCorrTilesTd(corr_tiles);
	}

	public float [][][][] getCorrTilesTd( // [tileY][tileX][pair][4*64] , read all available pairs
			float [][][][] corr_tiles)
	{
        int num_pairs = getNumUsedPairs();// Number of used pairs		num_pairs = num_pairs_in;
		final int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		int [] indices = getCorrIndices(); // also sets num_corr_tiles
		float [] fdata = getCorrTdData();
		int num_tiles = num_corr_tiles / num_pairs;
		int width = corr_tiles[0].length;
		for (int nt = 0; nt < num_tiles; nt++ ) {
			int nTile = (indices[nt * num_pairs] >> GPUTileProcessor.CORR_NTILE_SHIFT);
			int ty = nTile / width;
			int tx = nTile % width;

			corr_tiles[ty][tx] = new float [num_pairs][corr_size_td];
			for (int pair = 0; pair < num_pairs; pair++) {
				System.arraycopy(fdata, (nt * num_pairs + pair) * corr_size_td, corr_tiles[ty][tx][pair], 0, corr_size_td);
			}
		}
		return corr_tiles;
	}



	public int [] setCorrTilesComboTd(
			final float [][][] corr_tiles, // [tileY][tileX][4*64]
			int ipair) // just to set in the index low bits
	{
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		int tilesX = corr_tiles[0].length;
		int tilesY = corr_tiles.length;
		int [] indices = new int [tilesY * tilesX];  // as if all tiles are not null
		float [] fdata = new float [tilesY * tilesX * corr_size_td];  // as if all tiles are not null
		int ntile = 0;
		for (int ty = 0; ty < corr_tiles.length; ty++) {
			for (int tx = 0; tx < corr_tiles[0].length; tx++) {
				if (corr_tiles[ty][tx]!= null) {
					indices[ntile] = ((ty * tilesX + tx) << GPUTileProcessor.CORR_NTILE_SHIFT) + ipair;
					System.arraycopy(corr_tiles[ty][tx], 0, fdata, ntile * corr_size_td, corr_size_td);
					ntile++;
				}

			}
		}
		setCorrComboIndicesTdData(
				ntile,   // int    num_tiles,  // corr_indices, fdata may be longer than needed
				indices, // int [] corr_indices,
				fdata);  // float [] fdata);
		int [] indices_trim = new int [ntile];
		System.arraycopy(indices, 0, indices_trim, 0 , indices_trim.length);        	
		return indices_trim;
	}

	public float [][][] getCorrTilesComboTd() // [tileY][tileX][4*64] , read all available pairs
	{
		int tilesX =     img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =     img_height / GPUTileProcessor.DTT_SIZE;
		float [][][] corr_tiles = new float[tilesY][tilesX][]; // num_pairs
		return getCorrTilesComboTd(corr_tiles);
	}

	public float [][][] getCorrTilesComboTd( // [tileY][tileX][4*64] , read all available pairs
			float [][][] corr_tiles // should be initialized as [tilesX][tilesY][]
			)
	{
		final int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		int [] indices = getCorrComboIndices(); // num_corr_combo_tiles should be set earlier (from setCorrTilesComboTd or)
		float [] fdata = getCorrComboTdData();
		int width = corr_tiles[0].length;
		for (int nt = 0; nt < num_corr_combo_tiles; nt++ ) {
			int nTile = (indices[nt] >> GPUTileProcessor.CORR_NTILE_SHIFT);
			int ty = nTile / width;
			int tx = nTile % width;
			corr_tiles[ty][tx] = new float[corr_size_td];
			System.arraycopy(fdata, nt * corr_size_td, corr_tiles[ty][tx], 0, corr_size_td);
		}
		return corr_tiles;
	}



	public void setCorrIndicesTdData(
			int    num_tiles,  // corr_indices, fdata may be longer than needed
			int [] corr_indices,
			float [] fdata)
	{
		num_corr_tiles = num_tiles; // corr_indices.length;
		float [] fcorr_indices = new float [num_corr_tiles];
		for (int i = 0; i < num_corr_tiles; i++) {
			fcorr_indices[i] = Float.intBitsToFloat(corr_indices[i]);
		}
		cuMemcpyHtoD(gpu_corr_indices,   Pointer.to(fcorr_indices),  num_corr_tiles * Sizeof.FLOAT);
		float [] fnum_corr_tiles = {(float) num_corr_tiles};
		cuMemcpyHtoD(gpu_num_corr_tiles, Pointer.to(fnum_corr_tiles), 1 * Sizeof.FLOAT);
		// copy the correlation data
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		CUDA_MEMCPY2D copyH2D =   new CUDA_MEMCPY2D();
		copyH2D.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyH2D.srcHost =         Pointer.to(fdata);
		copyH2D.srcPitch =        corr_size_td*Sizeof.FLOAT; // width_in_bytes;
		copyH2D.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyH2D.dstDevice =       gpu_corrs_td; // src_dpointer;
		copyH2D.dstPitch =        corr_stride_td *Sizeof.FLOAT; // device_stride[0];
		copyH2D.WidthInBytes =    corr_size_td*Sizeof.FLOAT; // width_in_bytes;
		copyH2D.Height =          num_corr_tiles; // /4;
		cuMemcpy2D(copyH2D);
	}

	public int [] getCorrIndices() {
		int [] inum_corrs = new int[1];
		cuMemcpyDtoH(Pointer.to(inum_corrs), gpu_num_corr_tiles,  1 * Sizeof.INT);
		int num_corrs =      inum_corrs[0];
		int [] corr_indices = new int [num_corrs];
		cuMemcpyDtoH(Pointer.to(corr_indices), gpu_corr_indices,  num_corrs * Sizeof.INT);
		num_corr_tiles = num_corrs;
		return corr_indices;
	}
	
	
	public float [] getCorrTdData(){
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		float [] cpu_corrs = new float [ num_corr_tiles * corr_size_td];
		CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
		copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyD2H.srcDevice =       gpu_corrs_td;
		copyD2H.srcPitch =        corr_stride_td * Sizeof.FLOAT;
		copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyD2H.dstHost =         Pointer.to(cpu_corrs);
		copyD2H.dstPitch =        corr_size_td * Sizeof.FLOAT;
		copyD2H.WidthInBytes =    corr_size_td * Sizeof.FLOAT;
		copyD2H.Height =          num_corr_tiles;
		cuMemcpy2D(copyD2H); // run copy
		return cpu_corrs;
	}


	public void setCorrComboIndicesTdData(
			int    num_tiles,  // corr_combo_indices, fdata may be longer than needed
			int [] corr_combo_indices,
			float [] fdata)
	{
		num_corr_combo_tiles = num_tiles; // corr_combo_indices.length;
		float [] fcorr_combo_indices = new float [num_corr_combo_tiles];
		for (int i = 0; i < num_corr_combo_tiles; i++) {
			fcorr_combo_indices[i] = Float.intBitsToFloat(corr_combo_indices[i]);
		}
		cuMemcpyHtoD(gpu_corr_indices,   Pointer.to(fcorr_combo_indices),  num_corr_combo_tiles * Sizeof.FLOAT);
		//        	float [] fnum_corr_tiles = {(float) num_corr_tiles};
		//        	cuMemcpyHtoD(gpu_num_corr_tiles, Pointer.to(fnum_corr_tiles), 1 * Sizeof.FLOAT);
		// copy the correlation data
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		CUDA_MEMCPY2D copyH2D =   new CUDA_MEMCPY2D();
		copyH2D.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyH2D.srcHost =         Pointer.to(fdata);
		copyH2D.srcPitch =        corr_size_td*Sizeof.FLOAT; // width_in_bytes;
		copyH2D.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyH2D.dstDevice =       gpu_corrs_combo_td; // src_dpointer;
		copyH2D.dstPitch =        corr_stride_combo_td *Sizeof.FLOAT; // device_stride[0];
		copyH2D.WidthInBytes =    corr_size_td*Sizeof.FLOAT; // width_in_bytes;
		copyH2D.Height =          num_corr_combo_tiles; // /4;
		cuMemcpy2D(copyH2D);
	}



	public int [] getCorrComboIndices() {
		//        	float [] fnum_corrs = new float[1];
		//        	cuMemcpyDtoH(Pointer.to(fnum_corrs), gpu_num_corr_tiles,  1 * Sizeof.FLOAT);
		//        	int num_corrs =      Float.floatToIntBits(fnum_corrs[0]);
		float [] fcorr_combo_indices = new float [num_corr_combo_tiles];
		cuMemcpyDtoH(Pointer.to(fcorr_combo_indices), gpu_corr_combo_indices,  num_corr_combo_tiles * Sizeof.FLOAT);
		int [] corr_combo_indices = new int [num_corr_combo_tiles];
		for (int i = 0; i < num_corr_combo_tiles; i++) {
			corr_combo_indices[i] = Float.floatToIntBits(fcorr_combo_indices[i]);
		}
		return corr_combo_indices;
	}

	public float [] getCorrComboTdData(){
		int corr_size_td = 4 * GPUTileProcessor.DTT_SIZE * GPUTileProcessor.DTT_SIZE;
		float [] cpu_corrs = new float [ num_corr_combo_tiles * corr_size_td];
		CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
		copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyD2H.srcDevice =       gpu_corrs_combo_td;
		copyD2H.srcPitch =        corr_stride_combo_td * Sizeof.FLOAT;
		copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyD2H.dstHost =         Pointer.to(cpu_corrs);
		copyD2H.dstPitch =        corr_size_td * Sizeof.FLOAT;
		copyD2H.WidthInBytes =    corr_size_td * Sizeof.FLOAT;
		copyD2H.Height =          num_corr_combo_tiles;
		cuMemcpy2D(copyD2H); // run copy
		return cpu_corrs;
	}

	public void eraseGpuCorrs() { // temporary
//		float [] erase_buf = new float [corr_stride * num_corr_tiles];
		float [] erase_buf = new float [corr_stride * 120*80*64];
		Arrays.fill(erase_buf, 0.0f);
		cuMemcpyHtoD(gpu_corrs,   Pointer.to(erase_buf),  erase_buf.length * Sizeof.FLOAT);
	}
	
	
	public float [][] getCorr2D(int corr_rad){
		int corr_size = (2 * corr_rad + 1) * (2 * corr_rad + 1);
		float [] cpu_corrs = new float [ num_corr_tiles * corr_size];
		CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
		copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyD2H.srcDevice =       gpu_corrs;
		copyD2H.srcPitch =        corr_stride * Sizeof.FLOAT;

		copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyD2H.dstHost =         Pointer.to(cpu_corrs);
		copyD2H.dstPitch =        corr_size * Sizeof.FLOAT;

		copyD2H.WidthInBytes =    corr_size * Sizeof.FLOAT;
		copyD2H.Height =          num_corr_tiles;

		cuMemcpy2D(copyD2H); // run copy

		float [][] corrs = new float [num_corr_tiles][ corr_size];
		for (int ncorr = 0; ncorr < num_corr_tiles; ncorr++) {
			System.arraycopy(cpu_corrs, ncorr*corr_size, corrs[ncorr], 0, corr_size);
		}
		return corrs;
	}

	public float [][] getCorr2DCombo(int corr_rad){
		int corr_size = (2 * corr_rad + 1) * (2 * corr_rad + 1);
		float [] cpu_corrs = new float [ num_corr_combo_tiles * corr_size];
		CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
		copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyD2H.srcDevice =       gpu_corrs_combo;
		copyD2H.srcPitch =        corr_stride_combo * Sizeof.FLOAT;

		copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyD2H.dstHost =         Pointer.to(cpu_corrs);
		copyD2H.dstPitch =        corr_size * Sizeof.FLOAT;

		copyD2H.WidthInBytes =    corr_size * Sizeof.FLOAT;
		copyD2H.Height =          num_corr_combo_tiles;

		cuMemcpy2D(copyD2H); // run copy

		float [][] corrs = new float [num_corr_combo_tiles][ corr_size];
		for (int ncorr = 0; ncorr < num_corr_combo_tiles; ncorr++) {
			System.arraycopy(cpu_corrs, ncorr*corr_size, corrs[ncorr], 0, corr_size);
		}
		return corrs;
	}


	//	        
	/**
	 * Read extra data for macro generation: 4 DIFFs, 4 of R,  4 of B, 4 of G
	 * Available after texture        
	 * @return [ num_cams*(num_colors+1)][tilesX*tilesY] array for macro generation
	 */
	public float [][] getExtra(){
		int [] texture_indices = getTextureIndices();
		int num_tile_extra = num_cams*(num_colors+1);
		float [] diff_rgb_combo = new float[texture_indices.length * num_tile_extra];
		cuMemcpyDtoH(Pointer.to(diff_rgb_combo), gpu_diff_rgb_combo,  diff_rgb_combo.length * Sizeof.FLOAT);
		int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		int tilesY =  img_height / GPUTileProcessor.DTT_SIZE;
		float [][] extra = new float[num_tile_extra][tilesX*tilesY];
		for (int i = 0; i < texture_indices.length; i++) {
			if (((texture_indices[i] >> GPUTileProcessor.CORR_TEXTURE_BIT) & 1) != 0) {
				int ntile = (texture_indices[i] >>  GPUTileProcessor.CORR_NTILE_SHIFT);
				for (int l = 0; l < num_tile_extra; l++) {
					extra[l][ntile] = diff_rgb_combo[i * num_tile_extra + l];
				}
			}
		}
		return extra;
	}


	/**
	 * Get woi and RBGA image from the GPU after execRBGA call as 2/4 slices.
	 * device array has 4 pixels margins on each side, skip them here
	 * @param num_colors number of colors (1 or 3)
	 * @param woi should be initialized as Rectangle(). x,y,width, height will be populated (in pixels,)
	 * @return RBGA slices, last (alpha) in 0.0... 1.0 range, colors match input range
	 */
	public float [][] getRBGA(
			int     num_colors,
			Rectangle woi) { // will update to woi
		// first - read woi
		float [] fwoi = new float[4];
		cuMemcpyDtoH(Pointer.to(fwoi), gpu_woi,  4 * Sizeof.FLOAT);
		woi.x =      Float.floatToIntBits(fwoi[0]) * GPUTileProcessor.DTT_SIZE;
		woi.y =      Float.floatToIntBits(fwoi[1]) * GPUTileProcessor.DTT_SIZE;
		woi.width =  Float.floatToIntBits(fwoi[2]) * GPUTileProcessor.DTT_SIZE;
		woi.height = Float.floatToIntBits(fwoi[3]) * GPUTileProcessor.DTT_SIZE;
		float [][] rslt = new float[num_colors + 1][woi.width * woi.height];
		CUDA_MEMCPY2D copy_rbga =  new CUDA_MEMCPY2D();
		copy_rbga.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copy_rbga.srcDevice =       gpu_textures_rgba;
		copy_rbga.srcPitch =        texture_stride_rgba * Sizeof.FLOAT;
		copy_rbga.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copy_rbga.dstPitch =        woi.width * Sizeof.FLOAT;

		copy_rbga.WidthInBytes =    woi.width * Sizeof.FLOAT;
		copy_rbga.Height =          woi.height;
		copy_rbga.srcXInBytes =     4 * Sizeof.FLOAT;

		for (int ncol = 0; ncol<= num_colors; ncol++ ) {
			copy_rbga.dstHost =     Pointer.to(rslt[ncol]);
			copy_rbga.srcY =        4 + (woi.height +GPUTileProcessor.DTT_SIZE) * ncol;
			cuMemcpy2D(copy_rbga); // run copy
		}
		return rslt;
	}

	public float [] getFlatTextures(
			int     num_tiles,
			int     num_colors,
			boolean keep_weights){
		int texture_slices =     (num_colors + 1 + (keep_weights?(num_cams + num_colors + 1):0)); // number of texture slices
		int texture_slice_size = (2 * GPUTileProcessor.DTT_SIZE)* (2 * GPUTileProcessor.DTT_SIZE);        // number of (float) elements in a single slice of a tile
		int texture_tile_size =  texture_slices * texture_slice_size;   // number of (float) elements in a multi-slice tile
		int texture_size =       texture_tile_size * num_texture_tiles; // number of (float) elements in the whole texture
		float [] cpu_textures = new float [texture_size];
		CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
		copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyD2H.srcDevice =       gpu_textures;
		copyD2H.srcPitch =        texture_stride * Sizeof.FLOAT;

		copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyD2H.dstHost =         Pointer.to(cpu_textures);
		copyD2H.dstPitch =        texture_tile_size * Sizeof.FLOAT;

		copyD2H.WidthInBytes =    texture_tile_size * Sizeof.FLOAT;
		copyD2H.Height =          num_tiles; // num_texture_tiles;

		cuMemcpy2D(copyD2H); // run copy
		return cpu_textures;
	}

	public float [][][] getTextures( // todo - get rid of copying by multiple CUDA_MEMCPY2D?
			int     num_tiles,
			int     num_colors,
			boolean keep_weights){

		int texture_slices =     (num_colors + 1 + (keep_weights?(num_cams + num_colors + 1):0));
		int texture_slice_size = (2 * GPUTileProcessor.DTT_SIZE)* (2 * GPUTileProcessor.DTT_SIZE);
		int texture_tile_size =  texture_slices * texture_slice_size;
		//	        	int texture_size =       texture_tile_size * num_texture_tiles;
		float [] cpu_textures = getFlatTextures(
				num_tiles,
				num_colors,
				keep_weights);

		float [][][] textures = new float [num_texture_tiles][texture_slices][texture_tile_size];
		for (int ntile = 0; ntile < num_texture_tiles; ntile++) {
			for (int slice = 0; slice < texture_slices; slice++) {
				System.arraycopy(
						cpu_textures,                                           // src
						ntile * texture_tile_size + slice * texture_slice_size, // src offset
						textures[ntile][slice],                                 // dst
						0,                                                      // dst offset
						texture_slice_size);                                    // length to copy
			}
		}
		return textures;
	}

	public double [][][][] doubleTextures(
			Rectangle    woi,
			int []       indices,
			float [][][] ftextures,
			int          full_width,
			int          num_slices
			){
		int texture_slice_size = (2 * GPUTileProcessor.DTT_SIZE)* (2 * GPUTileProcessor.DTT_SIZE);
		double [][][][] textures = new double [woi.height][woi.width][num_slices][texture_slice_size];
		for (int indx = 0; indx < indices.length; indx++) if ((indices[indx] & (1 << GPUTileProcessor.LIST_TEXTURE_BIT)) != 0){
			int tile = indices[indx] >> GPUTileProcessor.CORR_NTILE_SHIFT;
		int tileX = tile % full_width;
		int tileY = tile / full_width;
		int wtileX = tileX - woi.x;
		int wtileY = tileY - woi.y;
		if ((wtileX < woi.width) && (wtileY < woi.height)) {
			for (int slice = 0; slice < num_slices; slice++) {
				for (int i = 0; i < texture_slice_size; i++) {
					textures[wtileY][wtileX][slice][i] = ftextures[indx][slice][i];
				}
			}
		}
		}
		return textures;
	}

	public double [][][][] doubleTextures( // may be accelerated with multithreading if needed.
			Rectangle    woi, // null or width and height match texture_tiles
			double [][][][] texture_tiles, // null or [tilesY][tilesX]
			int []       indices,
			float []     ftextures,
			int          full_width,
			int          num_slices,
			int          num_src_slices
			){
		if (num_slices > num_src_slices) {
			num_slices =num_src_slices;
		}
		int texture_slice_size = (2 * GPUTileProcessor.DTT_SIZE)* (2 * GPUTileProcessor.DTT_SIZE);
		int texture_tile_size = texture_slice_size * num_src_slices ;
		if ((woi == null) && (texture_tiles==null)) {
			System.out.println("doubleTextures(): woi and texture_tiles can not be simultaneously null");
			return null;
		}
		if (texture_tiles == null) {
			texture_tiles = new double [woi.height][woi.width][][];
		} else {
			woi.height = texture_tiles.length;
			woi.width = texture_tiles[0].length;
		}

		//        	double [][][][] textures = new double [woi.height][woi.width][num_slices][texture_slice_size];
		for (int indx = 0; indx < indices.length; indx++) if ((indices[indx] & (1 << GPUTileProcessor.LIST_TEXTURE_BIT)) != 0){
			int tile = (indices[indx] >> GPUTileProcessor.CORR_NTILE_SHIFT);
			int tileX = tile % full_width;
			int tileY = tile / full_width;
			int wtileX = tileX - woi.x;
			int wtileY = tileY - woi.y;
			texture_tiles[tileY][tileX] = new double [num_slices][texture_slice_size];
			if ((wtileX >=0 ) && (wtileX < woi.width) && (wtileY >= 0) && (wtileY < woi.height)) {
				for (int slice = 0; slice < num_slices; slice++) {
					for (int i = 0; i < texture_slice_size; i++) {
						texture_tiles[wtileY][wtileX][slice][i] = ftextures[indx * texture_tile_size + slice * texture_slice_size + i];
					}
				}
			}
		}
		return texture_tiles;
	}

	public float [][] getRBG (int ncam){
		int height = (img_height + GPUTileProcessor.DTT_SIZE);
		int width =  (img_width + GPUTileProcessor.DTT_SIZE);
		int rslt_img_size =      width * height;
		float [] cpu_corr_image = new float [ num_colors * rslt_img_size];
		int width_in_bytes = width *Sizeof.FLOAT;

		// for copying results to host
		CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
		copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
		copyD2H.srcDevice =       gpu_corr_images_h[ncam]; // ((test & 1) ==0) ? src_dpointer : dst_dpointer; // copy same data
		copyD2H.srcPitch =        imclt_stride*Sizeof.FLOAT;

		copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
		copyD2H.dstHost =         Pointer.to(cpu_corr_image);
		copyD2H.dstPitch =        width_in_bytes;

		copyD2H.WidthInBytes =    width_in_bytes;
//		copyD2H.Height =          3 * height; // /2;
		copyD2H.Height =          num_colors * height; // /2;

		cuMemcpy2D(copyD2H); // run copy

		float [][] fimg = new float [num_colors][ rslt_img_size];
		for (int ncol = 0; ncol < num_colors; ncol++) {
			System.arraycopy(cpu_corr_image, ncol*rslt_img_size, fimg[ncol], 0, rslt_img_size);
		}
		return fimg;
	}

	@Deprecated
	public void  getTileSubcamOffsets(
			final TpTask[]            tp_tasks,        // will use // modify to have offsets for 8 cameras
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets from the main camera
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)
	{
		final int quad_main = (geometryCorrection_main != null)? num_cams:0;
		final int quad_aux =  (geometryCorrection_aux != null)? num_cams:0;


		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(0);
		final Matrix [] corr_rots_main = geometryCorrection_main.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		Matrix [] corr_rots_aux0 = null;
		if (geometryCorrection_aux != null) {
			Matrix rigMatrix = geometryCorrection_aux.getRotMatrix(true);
			corr_rots_aux0 =  geometryCorrection_aux.getCorrVector().getRotMatrices(rigMatrix); // get array of per-sensor rotation matrices
		}
		final Matrix [] corr_rots_aux = corr_rots_aux0;
		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					int tileY,tileX; // , chn;
					//						showDoubleFloatArrays sdfa_instance = new showDoubleFloatArrays(); // just for debugging?
					double centerX; // center of aberration-corrected (common model) tile, X
					double centerY; //
					double disparity_main;
					double disparity_aux = 0.0;

					for (int nTile = ai.getAndIncrement(); nTile < tp_tasks.length; nTile = ai.getAndIncrement()) {

						tileY = tp_tasks[nTile].ty;
						tileX = tp_tasks[nTile].tx;
						if (tp_tasks[nTile].task == 0) {
							continue; // nothing to do for this tile
						}

						centerX = tileX * GPUTileProcessor.DTT_SIZE + GPUTileProcessor.DTT_SIZE/2; //  - shiftX;
						centerY = tileY * GPUTileProcessor.DTT_SIZE + GPUTileProcessor.DTT_SIZE/2; //  - shiftY;
						disparity_main = tp_tasks[nTile].target_disparity;
						if (geometryCorrection_aux != null) {
							disparity_aux =  disparity_main * geometryCorrection_aux.getDisparityRadius()/geometryCorrection_main.getDisparityRadius();
						}

						// TODO: move port coordinates out of color channel loop
						double [][] centersXY_main = null;
						double [][] centersXY_aux =  null;
						double [][] disp_dist_main = new double[quad_main][]; // used to correct 3D correlations
						double [][] disp_dist_aux =  new double[quad_aux][]; // used to correct 3D correlations

						if (geometryCorrection_main != null) {
							centersXY_main = geometryCorrection_main.getPortsCoordinatesAndDerivatives(
									geometryCorrection_main, //			GeometryCorrection gc_main,
									false,          // boolean use_rig_offsets,
									corr_rots_main, // Matrix []   rots,
									null,           //  Matrix [][] deriv_rots,
									null,           // double [][] pXYderiv, // if not null, should be double[8][]
									disp_dist_main,       // used to correct 3D correlations
									centerX,
									centerY,
									disparity_main); //  + disparity_corr);
							tp_tasks[nTile].xy = new float [centersXY_main.length][2];
							for (int i = 0; i < centersXY_main.length; i++) {
								tp_tasks[nTile].xy[i][0] = (float) centersXY_main[i][0];
								tp_tasks[nTile].xy[i][1] = (float) centersXY_main[i][1];
							}
						}
						if (geometryCorrection_aux != null) {
							centersXY_aux =  geometryCorrection_aux.getPortsCoordinatesAndDerivatives(
									geometryCorrection_main, //			GeometryCorrection gc_main,
									true,            // boolean use_rig_offsets,
									corr_rots_aux,   // Matrix []   rots,
									null,            //  Matrix [][] deriv_rots,
									null,            // double [][] pXYderiv, // if not null, should be double[8][]
									disp_dist_aux,   // used to correct 3D correlations
									centerX,
									centerY,
									disparity_aux); //  + disparity_corr);
							tp_tasks[nTile].xy_aux = new float [centersXY_aux.length][2];
							for (int i = 0; i < centersXY_aux.length; i++) {
								tp_tasks[nTile].xy_aux[i][0] = (float) centersXY_aux[i][0];
								tp_tasks[nTile].xy_aux[i][1] = (float) centersXY_aux[i][1];
							}
						}
						// acquisition time of the tiles centers in scanline times
						if (ers_delay != null) {
							for (int i = 0; i < quad_main; i++) ers_delay[0][i][nTile] = centersXY_main[i][1]-geometryCorrection_main.woi_tops[i];
							for (int i = 0; i < quad_aux; i++)  ers_delay[1][i][nTile] = centersXY_aux[i][1]- geometryCorrection_aux.woi_tops[i];
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
	}

	/**
	 * Prepare GPU tasks for interscene accumulation - instead of the uniform grid, pixel coordinates and disparity are provided
	 * by the caller. They are calculated by recalculating from the reference scene after appropriate transformation (shift, rotation
	 * and ERS correction)
	 * @param calcPortsCoordinatesAndDerivatives -  calculate port coordinates when using CPU mode. GPU calculates them 
	 * @param pXpYD Array of per-tile pX, pY and disparity triplets (or nulls for undefined tiles).
	 * @param geometryCorrection GeometryCorrection instance for the camera.
	 * @param disparity_corr Disparity correction at infinity
	 * @param margin Skip tile if at least one channel tile center is closer to the image edge than this margin.
	 * @param valid_tiles Optional (if not null) should be initialized as boolean [tiles] - will contain valid tiles
	 * @param threadsMax Maximal number of threads to run concurrently.
	 * @return Array of TpTask instances (fully prepared) to be fed to the GPU
	 */
	public TpTask[]  setInterTasks(
			final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
			final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
			final GeometryCorrection  geometryCorrection,
			final double              disparity_corr,
			final int                 margin,      // do not use tiles if their centers are closer to the edges
			final boolean []          valid_tiles,            
			final int                 threadsMax)  // maximal number of threads to launch
	{
		return setInterTasks(
				num_cams,           // final int                 num_cams,
				img_width,          // final int                 img_width,
				calcPortsCoordinatesAndDerivatives, // final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
				pXpYD,              // final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
				geometryCorrection, // final GeometryCorrection  geometryCorrection,
				disparity_corr,     // final double              disparity_corr,
				margin,             // final int                 margin,      // do not use tiles if their centers are closer to the edges
				valid_tiles,        // final boolean []          valid_tiles,            
				threadsMax);        // final int                 threadsMax);  // maximal number of threads to launch
	}

	public static TpTask[]  setInterTasks(
			final int                 num_cams,
			final int                 img_width,
			final boolean             calcPortsCoordinatesAndDerivatives, // GPU can calculate them centreXY
			final double [][]         pXpYD, // per-tile array of pX,pY,disparity triplets (or nulls)
			final GeometryCorrection  geometryCorrection,
			final double              disparity_corr,
			final int                 margin,      // do not use tiles if their centers are closer to the edges
			final boolean []          valid_tiles,            
			final int                 threadsMax)  // maximal number of threads to launch
	{
		int num_pairs = Correlation2d.getNumPairs(num_cams);
		//change to fixed 511?
		final int task_code = ((1 << num_pairs)-1) << GPUTileProcessor.TASK_CORR_BITS; //  correlation only
		final double min_px = margin; 
		final double max_px = img_width - 1 - margin;
		final double [] min_py = new double[num_cams] ;
		final double [] max_py = new double[num_cams] ;
		for (int i = 0; i < num_cams; i++) {
			min_py [i] = margin + (calcPortsCoordinatesAndDerivatives? geometryCorrection.getWOITops()[i] : 0);
			// camera_heights array is only set during conditionImageSet(), not called by the intersceneAccumulate()
			// That was correct, as all scenes should be conditioned
//			max_py [i] = geometryCorrection.getWOITops()[i] + geometryCorrection.getCameraHeights()[i] - 1 - margin;
			max_py [i] = geometryCorrection.getSensorWH()[1] - 1 - margin; // same for all channels?
			//.getSensorWH()[0]
		}
		if (valid_tiles!=null) {
			Arrays.fill(valid_tiles, false);
		}
		final int tilesX =  img_width / GPUTileProcessor.DTT_SIZE;
		final int tiles = pXpYD.length;
		final Matrix [] corr_rots = geometryCorrection.getCorrVector().getRotMatrices(); // get array of per-sensor rotation matrices
		final int quad_main = (geometryCorrection != null)? num_cams:0;
		final Thread[] threads = ImageDtt.newThreadArray(threadsMax);
		final AtomicInteger ai = new AtomicInteger(00);
		final AtomicInteger aTiles = new AtomicInteger(0);
		final TpTask[] tp_tasks = new TpTask[tiles]; // aTiles.get()];

		for (int ithread = 0; ithread < threads.length; ithread++) {
			threads[ithread] = new Thread() {
				@Override
				public void run() {
					//    					for (int indx = ai.getAndIncrement(); indx < tp_tasks.length; indx = ai.getAndIncrement()) {
					//   						int nTile = tile_indices[indx];
					for (int nTile = ai.getAndIncrement(); nTile < tiles; nTile = ai.getAndIncrement()) if (pXpYD[nTile] != null) {
						int tileY = nTile / tilesX;
						int tileX = nTile % tilesX;
						TpTask tp_task = new TpTask(num_cams, tileX, tileY);
//						tp_task.ty = tileY;
//						tp_task.tx = tileX;
						tp_task.task = task_code;
						double disparity = pXpYD[nTile][2] + disparity_corr;
						tp_task.target_disparity = (float) disparity; // will it be used?
						double [] centerXY = pXpYD[nTile];
						tp_task.setCenterXY(centerXY); // this pair of coordinates will be used by GPU to set tp_task.xy and task.disp_dist!
						boolean bad_margins = false;
						if (calcPortsCoordinatesAndDerivatives) {
							double [][] disp_dist = new double[quad_main][]; // used to correct 3D correlations (not yet used here)
							double [][] centersXY_main = geometryCorrection.getPortsCoordinatesAndDerivatives(
									geometryCorrection, //			GeometryCorrection gc_main,
									false,          // boolean use_rig_offsets,
									corr_rots, // Matrix []   rots,
									null,           //  Matrix [][] deriv_rots,
									null,           // double [][] pXYderiv, // if not null, should be double[8][]
									disp_dist,       // used to correct 3D correlations
									pXpYD[nTile][0],
									pXpYD[nTile][1],
									disparity); //  + disparity_corr);
							tp_task.setDispDist(disp_dist);
							tp_task.xy = new float [centersXY_main.length][2];
							for (int i = 0; i < centersXY_main.length; i++) {
								if (    (centersXY_main[i][0] < min_px) ||    (centersXY_main[i][0] > max_px) ||
										(centersXY_main[i][1] < min_py[i]) || (centersXY_main[i][1] > max_py[i])) {
									bad_margins = true;
									break;
								}
								tp_task.xy[i][0] = (float) centersXY_main[i][0];
								tp_task.xy[i][1] = (float) centersXY_main[i][1];
							}
						} else { // only check center for margins
							if (    (centerXY[0] < min_px) ||    (centerXY[0] > max_px) ||
									(centerXY[1] < min_py[0]) || (centerXY[1] > max_py[0])) {
								bad_margins = true;
//								break;
							}

						}
						if (bad_margins) {
							continue;
						}
						tp_tasks[aTiles.getAndIncrement()] = tp_task;
						if (valid_tiles!=null) {
							valid_tiles[nTile] = true;
						}
					}
				}
			};
		}
		ImageDtt.startAndJoin(threads);
		final TpTask[] tp_tasks_out = new TpTask[aTiles.get()];
		System.arraycopy(tp_tasks, 0, tp_tasks_out, 0, tp_tasks_out.length);
		return tp_tasks_out;
	}


	
	
	
	



	public void setLpfRbg(
			float [][] lpf_rbg, // 4 64-el. arrays: r,b,g,m
			boolean    debug)
	{

		int l = lpf_rbg[0].length; // 64
		float []   lpf_flat = new float [lpf_rbg.length * l];
		for (int i = 0; i < lpf_rbg.length; i++) {
			for (int j = 0; j < l; j++) {
				lpf_flat[j + i * l] = lpf_rbg[i][j];
			}
		}

		CUdeviceptr constantMemoryPointer = new CUdeviceptr();
		long constantMemorySizeArray[] = { 0 };
		cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray,  this.gpuTileProcessor.module, "lpf_data"); //CUDA_ERROR_ILLEGAL_ADDRESS after re-run
		int constantMemorySize = (int)constantMemorySizeArray[0];
		if (debug) System.out.println("constantMemoryPointer: " + constantMemoryPointer);
		if (debug) System.out.println("constantMemorySize: " + constantMemorySize);
		cuMemcpyHtoD(constantMemoryPointer, Pointer.to(lpf_flat), constantMemorySize);
		if (debug) System.out.println();
		/*
            if (debug) {
        		for (int i = 0; i < lpf_flat.length; i++) {
        			System.out.print(String.format("%8.5f", lpf_flat[i]));
        			if (((i+1) % 16) == 0) {
        				System.out.println();
        			}
        		}
            }
		 */
	}

	public void setLpfCorr(
			String   const_name, // "lpf_corr"
			float [] lpf_flat,
			boolean  debug)

	{
		CUdeviceptr constantMemoryPointer = new CUdeviceptr();
		long constantMemorySizeArray[] = { 0 };
		cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray,  this.gpuTileProcessor.module, const_name);
		int constantMemorySize = (int)constantMemorySizeArray[0];
		if (debug) System.out.println("constantMemoryPointer: " + constantMemoryPointer);
		if (debug) System.out.println("constantMemorySize: " + constantMemorySize);
		cuMemcpyHtoD(constantMemoryPointer, Pointer.to(lpf_flat), constantMemorySize);
		if (debug) System.out.println();
		/*
            if (debug) {
        		for (int i = 0; i < lpf_flat.length; i++) {
        			System.out.print(String.format("%8.5f", lpf_flat[i]));
        			if (((i+1) % 16) == 0) {
        				System.out.println();
        			}
        		}
            }
		 */
	}

	public float [] floatSetCltLpfFd(
			double   sigma) {
		int dct_size = GPUTileProcessor.DTT_SIZE;
		DttRad2 dtt = new DttRad2(dct_size);
		double [] clt_fd = dtt.dttt_iiie(setCltLpf(sigma));
		int l = dct_size*dct_size;
		float []   lpf_flat = new float [l];
		for (int j = 0; j < l; j++) {
			lpf_flat[j] = (float) (clt_fd[j]*2*dct_size);
		}
		return lpf_flat;
	}

	public double [] doubleSetCltLpfFd(
			double   sigma) {
		int dct_size = GPUTileProcessor.DTT_SIZE;
		DttRad2 dtt = new DttRad2(dct_size);
		double [] clt_fd = dtt.dttt_iiie(setCltLpf(sigma));
		int l = dct_size*dct_size;
		double []   lpf_flat = new double [l];
		for (int j = 0; j < l; j++) {
			lpf_flat[j] = (float) (clt_fd[j]*2*dct_size);
		}
		return lpf_flat;
	}

	public double [] setCltLpf(
			double   sigma)
	{
		int dct_size = GPUTileProcessor.DTT_SIZE;
		double [] lpf = new double [dct_size*dct_size];
		int dct_len = dct_size * dct_size;
		if (sigma == 0.0f) {
			lpf[0] = 1.0f;
			for (int i = 1; i < dct_len; i++){
				lpf[i] = 0.0f;
			}
		} else {
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					lpf[i*dct_size+j] = (float) Math.exp(-(i*i+j*j)/(2*sigma));
				}
			}
			// normalize
			double sum = 0;
			for (int i = 0; i < dct_size; i++){
				for (int j = 0; j < dct_size; j++){
					double d = 	lpf[i*dct_size+j];
					d*=Math.cos(Math.PI*i/(2*dct_size))*Math.cos(Math.PI*j/(2*dct_size));
					if (i > 0) d*= 2.0;
					if (j > 0) d*= 2.0;
					sum +=d;
				}
			}
			for (int i = 0; i< dct_len; i++){
				lpf[i] /= sum;
			}
		}
		return lpf;
	}
	public int getGpu_debug_level() {
		return gpu_debug_level;
	}
	public void setGpu_debug_level(int gpu_debug_level) {
		this.gpu_debug_level = gpu_debug_level;
	}

} // end of public class GpuQuad