package com.elphel.imagej.gpu;
/**
** -----------------------------------------------------------------------------**
** GPUTileProcessor.java
**
** GPU acceleration for the Tile Processor
**
**
** Copyright (C) 2018 Elphel, Inc.
**
** -----------------------------------------------------------------------------**
**
**  GPUTileProcessor.java is free software: you can redistribute it and/or modify
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

// Uses code by Marco Hutter - http://www.jcuda.org
import static jcuda.driver.JCudaDriver.cuCtxCreate;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import static jcuda.driver.JCudaDriver.cuMemAlloc;
import static jcuda.driver.JCudaDriver.cuMemAllocPitch;
import static jcuda.driver.JCudaDriver.cuMemcpy2D;
import static jcuda.driver.JCudaDriver.cuMemcpyHtoD;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleGetGlobal;
import static jcuda.driver.JCudaDriver.cuModuleLoadData;
import static jcuda.nvrtc.JNvrtc.nvrtcCompileProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcCreateProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcDestroyProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcGetPTX;
import static jcuda.nvrtc.JNvrtc.nvrtcGetProgramLog;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.tileprocessor.DttRad2;
import com.elphel.imagej.tileprocessor.GeometryCorrection;
import com.elphel.imagej.tileprocessor.ImageDtt;

import Jama.Matrix;
import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUDA_MEMCPY2D;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmemorytype;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import jcuda.nvrtc.JNvrtc;
import jcuda.nvrtc.nvrtcProgram;

public class GPUTileProcessor {
	static String GPU_KERNEL_FILE = "dtt8x8.cuh";
	static String [] GPU_KERNEL_FILES = {"dtt8x8.cuh","TileProcessor.cuh"};
	static String GPU_CONVERT_CORRECT_TILES_NAME = "convert_correct_tiles";
	static String GPU_IMCLT_RBG_NAME = "imclt_rbg";
//  pass some defines to gpu source code with #ifdef JCUDA
	public static int DTT_SIZE =                8;
	static int THREADSX =         DTT_SIZE;
	public static int NUM_CAMS =                4;
	static int NUM_COLORS =              3;
	public static int IMG_WIDTH =            2592;
	public static int IMG_HEIGHT =           1936;
	static int KERNELS_HOR =           164;
	static int KERNELS_VERT =          123;
	static int KERNELS_LSTEP =           4;
	static int THREADS_PER_TILE =        8;
	static int TILES_PER_BLOCK =         4; // 8 - slower
	static int IMCLT_THREADS_PER_TILE = 16;
	static int IMCLT_TILES_PER_BLOCK =   4;

	static int TPTASK_SIZE =   NUM_CAMS * 2 + 2;
	static int CLTEXTRA_SIZE = 8;
	static int KERN_TILES = KERNELS_HOR *  KERNELS_VERT * NUM_COLORS;
	static int KERN_SIZE =  KERN_TILES * 4 * 64;


    int DTTTEST_BLOCK_WIDTH =        32; // may be read from the source code
    int DTTTEST_BLOCK_HEIGHT =       16; // may be read from the source code

    public boolean kernels_set = false;
    public boolean bayer_set =   false;

    private CUfunction GPU_CONVERT_CORRECT_TILES_kernel = null;
    private CUfunction GPU_IMCLT_RBG_kernel =             null;
    // CPU arrays of pointers to GPU memory
    // These arrays may go to method, they are here just to be able to free GPU memory if needed
    private CUdeviceptr [] gpu_kernels_h =        new CUdeviceptr[NUM_CAMS];
    private CUdeviceptr [] gpu_kernel_offsets_h = new CUdeviceptr[NUM_CAMS];
    private CUdeviceptr [] gpu_bayer_h =          new CUdeviceptr[NUM_CAMS];
    private CUdeviceptr [] gpu_clt_h =            new CUdeviceptr[NUM_CAMS];
//    private CUdeviceptr [] gpu_lpf_h =            new CUdeviceptr[NUM_COLORS];
    private CUdeviceptr [] gpu_corr_images_h=     new CUdeviceptr[NUM_CAMS];

    // GPU pointers to array of GPU pointers
    private CUdeviceptr gpu_kernels =             new CUdeviceptr();
    private CUdeviceptr gpu_kernel_offsets =      new CUdeviceptr();
    private CUdeviceptr gpu_bayer =               new CUdeviceptr();

    private CUdeviceptr gpu_tasks =               new CUdeviceptr();

    private CUdeviceptr gpu_clt =                 new CUdeviceptr();
//    private
    CUmodule    module; // to access constants memory
//    private CUdeviceptr gpu_lpf =            new CUdeviceptr();
    private int mclt_stride;
    private int imclt_stride;
    public int num_task_tiles;

    public class TpTask {
    	public int   task;
    	public float target_disparity;

    	public int ty;
    	public int tx;
    	public float[][] xy = null;
    	public float[][] xy_aux = null;
    	public TpTask() {}

    	public TpTask(int tx, int ty, float target_disparity, int task ) {
    		this.tx = tx;
    		this.ty = ty;
    		this.target_disparity = target_disparity;
    		this.task = task;
    	}

    	// convert this class instance to float array to match layout of the C struct
    	public float [] asFloatArray(boolean use_aux) {
    		float [] flt = new float [NUM_CAMS * 2 + 2];
    		return asFloatArray(flt, 0, use_aux);
    	}
    	// convert this class instance to float array to match layout of the C struct,
    	// fill existing float array from the specified index
    	public float [] asFloatArray(float [] flt, int indx, boolean use_aux) {
    		flt[indx++] = Float.intBitsToFloat(task);
    		flt[indx++] = Float.intBitsToFloat(tx + (ty << 16));
    		float [][] offsets = use_aux? this.xy_aux: this.xy;
    		for (int i = 0; i < NUM_CAMS; i++) {
    			flt[indx++] = offsets[i][0];
    			flt[indx++] = offsets[i][1];
    		}
    		return flt;
    	}
    }

    public class CltExtra{
    	public float data_x;   // kernel data is relative to this displacement X (0.5 pixel increments)
    	public float data_y;   // kernel data is relative to this displacement Y (0.5 pixel increments)
    	public float center_x; // actual center X (use to find derivatives)
    	public float center_y; // actual center X (use to find derivatives)
    	public float dxc_dx;   // add this to data_x per each pixel X-shift relative to the kernel center location
    	public float dxc_dy;   // same per each Y-shift pixel
    	public float dyc_dx;
    	public float dyc_dy;
    	public CltExtra() {}
    	public CltExtra(
    	    	float data_x,   // kernel data is relative to this displacement X (0.5 pixel increments)
    	    	float data_y,   // kernel data is relative to this displacement Y (0.5 pixel increments)
    	    	float center_x, // actual center X (use to find derivatives)
    	    	float center_y, // actual center X (use to find derivatives)
    	    	float dxc_dx,   // add this to data_x per each pixel X-shift relative to the kernel center location
    	    	float dxc_dy,   // same per each Y-shift pixel
    	    	float dyc_dx,
    	    	float dyc_dy)
    	{
    		this.data_x = data_x;
    		this.data_y = data_y;
    		this.center_x = center_x;
    		this.center_y = center_y;
    		this.dxc_dx = dxc_dx;
    		this.dxc_dy = dxc_dy;
    		this.dyc_dx = dyc_dx;
    		this.dyc_dy = dyc_dy;
    	}
    	public CltExtra(float [] data, int indx)
    	{
    		this.data_x =   data[indx++];
    		this.data_y =   data[indx++];
    		this.center_x = data[indx++];
    		this.center_y = data[indx++];
    		this.dxc_dx =   data[indx++];
    		this.dxc_dy =   data[indx++];
    		this.dyc_dx =   data[indx++];
    		this.dyc_dy =   data[indx++];
    	}
    	public float [] asFloatArray() {
    		float [] flt = new float [8];
    		return asFloatArray(flt, 0);
    	}
    	public float [] asFloatArray(float [] flt, int indx) {
    		flt[indx++] = this.data_x;
    		flt[indx++] = this.data_y;
    		flt[indx++] = this.center_x;
    		flt[indx++] = this.center_y;
    		flt[indx++] = this.dxc_dx;
    		flt[indx++] = this.dxc_dy;
    		flt[indx++] = this.dyc_dx;
    		flt[indx++] = this.dyc_dy;
    		return flt;
    	}
    };

    private static long getPointerAddress(CUdeviceptr p)

    {
        // WORKAROUND until a method like CUdeviceptr#getAddress exists
        class PointerWithAddress extends Pointer
        {
            PointerWithAddress(Pointer other)
            {
                super(other);
            }
            long getAddress()
            {
                return getNativePointer() + getByteOffset();
            }
        }
        return new PointerWithAddress(p).getAddress();
    }

    public GPUTileProcessor(String cuda_project_directory) throws IOException
    {
    	// From code by Marco Hutter - http://www.jcuda.org
        // Enable exceptions and omit all subsequent error checks
        JCudaDriver.setExceptionsEnabled(true);
        JNvrtc.setExceptionsEnabled(true);

        // Initialize the driver and create a context for the first device.
        cuInit(0);
        CUdevice device = new CUdevice();
        cuDeviceGet(device, 0);
        CUcontext context = new CUcontext();
        cuCtxCreate(context, 0, device);

        // Obtain the CUDA source code from the CUDA file
        // Get absolute path to the file in resource foldder, then read it as a normal file.
        // When using just Eclipse resources - it does not notice that the file
        // was edited (happens frequently during kernel development).
        ClassLoader classLoader = getClass().getClassLoader();
        String kernelSource =
        		"#define JCUDA\n"+
        				"#define DTT_SIZE " +               DTT_SIZE+"\n"+
        				"#define THREADSX " +               THREADSX+"\n"+
        				"#define NUM_CAMS " +               NUM_CAMS+"\n"+
        				"#define NUM_COLORS " +             NUM_COLORS+"\n"+
        				"#define IMG_WIDTH " +              IMG_WIDTH+"\n"+
        				"#define IMG_HEIGHT " +             IMG_HEIGHT+"\n"+
        				"#define KERNELS_HOR " +            KERNELS_HOR+"\n"+
        				"#define KERNELS_VERT " +           KERNELS_VERT+"\n"+
        				"#define KERNELS_LSTEP " +          KERNELS_LSTEP+"\n"+
        				"#define THREADS_PER_TILE " +       THREADS_PER_TILE+"\n"+
        				"#define TILES_PER_BLOCK " +        TILES_PER_BLOCK+"\n"+
        				"#define IMCLT_THREADS_PER_TILE " + IMCLT_THREADS_PER_TILE+"\n"+
        				"#define IMCLT_TILES_PER_BLOCK " +  IMCLT_TILES_PER_BLOCK+"\n";

        for (String src_file:GPU_KERNEL_FILES) {
        	File file = null;
        	if ((cuda_project_directory == null) || (cuda_project_directory == "")) {
        		file = new File(classLoader.getResource(src_file).getFile());
        		System.out.println("Loading resource "+file);
        	} else {
        		File src_dir = new File(cuda_project_directory, "src");
        		file = new File(src_dir.getPath(), src_file);
        		System.out.println("Loading resource "+file);
        	}
        	System.out.println(file.getAbsolutePath());
        	String cuFileName = file.getAbsolutePath(); // /home/eyesis/workspace-python3/nvidia_dct8x8/src/dtt8x8.cuh";// "dtt8x8.cuh";
        	String sourceFile = readFileAsString(cuFileName); // readResourceAsString(cuFileName);
        	if (sourceFile == null) {
        		String msg = "Could not read the kernel source code";
        		IJ.showMessage("Error",	msg);
        		new IllegalArgumentException (msg);
        	}
        	kernelSource += sourceFile;

        }
        // Create the kernel functions (first - just test)
        String [] func_names = {GPU_CONVERT_CORRECT_TILES_NAME, GPU_IMCLT_RBG_NAME};
        CUfunction[] functions = createFunctions(kernelSource, func_names);
        this.GPU_CONVERT_CORRECT_TILES_kernel = functions[0];
        this.GPU_IMCLT_RBG_kernel =             functions[1];
        System.out.println("GPU kernel functions initialized");
        System.out.println("Sizeof.POINTER="+Sizeof.POINTER);
        System.out.println(GPU_IMCLT_RBG_kernel.toString());

        // Init data arrays
        int tilesX =  IMG_WIDTH / DTT_SIZE;
        int tilesY =  IMG_HEIGHT / DTT_SIZE;

        for (int ncam = 0; ncam < NUM_CAMS; ncam++) {
        	gpu_kernels_h[ncam] =        new CUdeviceptr();
        	cuMemAlloc(gpu_kernels_h[ncam],KERN_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
        	gpu_kernel_offsets_h[ncam] = new CUdeviceptr();
        	cuMemAlloc(gpu_kernel_offsets_h[ncam],KERN_TILES * CLTEXTRA_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
        	gpu_bayer_h[ncam] =          new CUdeviceptr();
            long [] device_stride = new long [1];
            cuMemAllocPitch (
            		gpu_bayer_h[ncam],        // CUdeviceptr dptr,
            		device_stride,            // long[] pPitch,
            		IMG_WIDTH * Sizeof.FLOAT, // long WidthInBytes,
            		IMG_HEIGHT,               // long Height,
                    Sizeof.FLOAT);            // int ElementSizeBytes)
            mclt_stride = (int)(device_stride[0] / Sizeof.FLOAT);
            gpu_corr_images_h[ncam] =  new CUdeviceptr();
            cuMemAllocPitch (
            		gpu_corr_images_h[ncam],               // CUdeviceptr dptr,
            		device_stride,                         // long[] pPitch,
            		(IMG_WIDTH + DTT_SIZE) * Sizeof.FLOAT, // long WidthInBytes,
            		3*(IMG_HEIGHT + DTT_SIZE),// long Height,
                    Sizeof.FLOAT);            // int ElementSizeBytes)
            imclt_stride = (int)(device_stride[0] / Sizeof.FLOAT);
            gpu_clt_h[ncam] = new CUdeviceptr();
        	cuMemAlloc(gpu_clt_h[ncam],tilesY * tilesX * NUM_COLORS * 4 * DTT_SIZE * DTT_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)

        }
        // now create device arrays pointers
        if (Sizeof.POINTER != Sizeof.LONG) {
    		String msg = "Sizeof.POINTER != Sizeof.LONG";
    		IJ.showMessage("Error",	msg);
    		new IllegalArgumentException (msg);
        }
    	cuMemAlloc(gpu_kernels,        NUM_CAMS * Sizeof.POINTER);
    	cuMemAlloc(gpu_kernel_offsets, NUM_CAMS * Sizeof.POINTER);
    	cuMemAlloc(gpu_bayer,          NUM_CAMS * Sizeof.POINTER);
    	cuMemAlloc(gpu_clt,            NUM_CAMS * Sizeof.POINTER);
    	long [] gpu_kernels_l =        new long [NUM_CAMS];
    	long [] gpu_kernel_offsets_l = new long [NUM_CAMS];
    	long [] gpu_bayer_l =          new long [NUM_CAMS];
    	long [] gpu_clt_l =            new long [NUM_CAMS];

    	for (int ncam = 0; ncam < NUM_CAMS; ncam++) gpu_kernels_l[ncam] =        getPointerAddress(gpu_kernels_h[ncam]);
        cuMemcpyHtoD(gpu_kernels, Pointer.to(gpu_kernels_l),                     NUM_CAMS * Sizeof.POINTER);

        for (int ncam = 0; ncam < NUM_CAMS; ncam++) gpu_kernel_offsets_l[ncam] = getPointerAddress(gpu_kernel_offsets_h[ncam]);
        cuMemcpyHtoD(gpu_kernel_offsets, Pointer.to(gpu_kernel_offsets_l),       NUM_CAMS * Sizeof.POINTER);

        for (int ncam = 0; ncam < NUM_CAMS; ncam++) gpu_bayer_l[ncam] =          getPointerAddress(gpu_bayer_h[ncam]);
        cuMemcpyHtoD(gpu_bayer, Pointer.to(gpu_bayer_l),                         NUM_CAMS * Sizeof.POINTER);

        for (int ncam = 0; ncam < NUM_CAMS; ncam++) gpu_clt_l[ncam] =            getPointerAddress(gpu_clt_h[ncam]);
        cuMemcpyHtoD(gpu_clt, Pointer.to(gpu_clt_l),                             NUM_CAMS * Sizeof.POINTER);

        // Set task array
    	cuMemAlloc(gpu_tasks,       tilesX * tilesY * TPTASK_SIZE * Sizeof.POINTER);

    }


    public void setTasks(TpTask [] tile_tasks, boolean use_aux)
    {
    	num_task_tiles = tile_tasks.length;
    	float [] ftasks = new float [TPTASK_SIZE * num_task_tiles];
    	for (int i = 0; i < num_task_tiles; i++) {
    		tile_tasks[i].asFloatArray(ftasks, i* TPTASK_SIZE, use_aux);
    	}
        cuMemcpyHtoD(gpu_tasks,        Pointer.to(ftasks),         TPTASK_SIZE * num_task_tiles * Sizeof.FLOAT);
    }

    public void setConvolutionKernel(
    		float [] kernel,  // [tileY][tileX][color][..]
    		float [] kernel_offsets,
    		int ncam) {
        cuMemcpyHtoD(gpu_kernels_h[ncam],        Pointer.to(kernel),         KERN_SIZE * Sizeof.FLOAT);
        cuMemcpyHtoD(gpu_kernel_offsets_h[ncam], Pointer.to(kernel_offsets), KERN_TILES * CLTEXTRA_SIZE * Sizeof.FLOAT);
    }

    public void setConvolutionKernels(
    		double [][][][][][] clt_kernels,
    		boolean force)
    {
    	boolean transpose = true;
    	if (kernels_set && ! force) {
    		return;
    	}
    	int num_kernels =
    			clt_kernels[0][0].length * //tilesY
    			clt_kernels[0][0][0].length * //tilesX
    			clt_kernels[0].length;  //colors
    	int kernel_length =	num_kernels * 4 * DTT_SIZE * DTT_SIZE;

    	float [] fkernel =  new float [kernel_length];
    	float [] foffsets = new float [num_kernels * CLTEXTRA_SIZE];
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
    	kernels_set = true;
    }

    public void setBayerImage(
    		float [] bayer_image,
    		int ncam) {
        CUDA_MEMCPY2D copyH2D =   new CUDA_MEMCPY2D();
        copyH2D.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
        copyH2D.srcHost =         Pointer.to(bayer_image);
        copyH2D.srcPitch =        IMG_WIDTH*Sizeof.FLOAT; // width_in_bytes;
        copyH2D.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
        copyH2D.dstDevice =       gpu_bayer_h[ncam]; // src_dpointer;
        copyH2D.dstPitch =        mclt_stride *Sizeof.FLOAT; // device_stride[0];
        copyH2D.WidthInBytes =    IMG_WIDTH*Sizeof.FLOAT; // width_in_bytes;
        copyH2D.Height =          IMG_HEIGHT; // /4;
        cuMemcpy2D(copyH2D);
    }
    // combines 3 bayer channels into one and transfers to GPU memory
    public void setBayerImages(
			double [][][]       bayer_data,
    		boolean                  force) {
    	if (bayer_set && !force) {
    		return;
    	}
    	float [] fbayer = new float [bayer_data[0][0].length];
		for (int ncam = 0; ncam < bayer_data.length; ncam++) {
			for (int i = 0; i <  bayer_data[ncam][0].length; i++) {
				fbayer[i] = (float) (bayer_data[ncam][0][i] + bayer_data[ncam][1][i] + bayer_data[ncam][2][i]);
			}
		    setBayerImage(
		    		fbayer, // float [] bayer_image,
		    		ncam); // int ncam)
		}
		bayer_set = true;
    }

    // prepare tasks for full frame, same dispaity.
    // need to run setTasks(TpTask [] tile_tasks, boolean use_aux) to format/transfer to GPU memory
    public TpTask [] setFullFrameImages(
    		float                     target_disparity, // apply same disparity to all tiles
    		boolean                   use_master,
    		boolean                   use_aux,
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel) {
        int tilesX =  IMG_WIDTH / DTT_SIZE;
        int tilesY =  IMG_HEIGHT / DTT_SIZE;
    	float [] target_disparities = new float [tilesX * tilesY];
    	if (target_disparity != 0.0) {
    		for (int i = 0; i <target_disparities.length; i++ ) target_disparities[i] = target_disparity;
    	}
    	return setFullFrameImages(
        		target_disparities, // should be tilesX*tilesY long
        		use_master,
        		use_aux,
    			geometryCorrection_main,
    			geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
    			ers_delay,        // if not null - fill with tile center acquisition delay
    			threadsMax,  // maximal number of threads to launch
    			debugLevel);
    }

    public TpTask [] setFullFrameImages(
    		float []                  target_disparities, // should be tilesX*tilesY long
    		boolean                   use_master,
    		boolean                   use_aux,
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)
    {
        int tilesX =  IMG_WIDTH / DTT_SIZE;
        int tilesY =  IMG_HEIGHT / DTT_SIZE;

    	TpTask [] tp_tasks = new TpTask[tilesX*tilesY];
    	int indx = 0;
    	for (int ty = 0; ty < tilesY; ty++) {
        	for (int tx = 0; tx < tilesX; tx++) {
        		tp_tasks[indx] = new TpTask(tx,ty, target_disparities[indx], 1); // task == 1 for now
        		indx++;
        	}
    	}
    	getTileSubcamOffsets(
    			tp_tasks,                                    // final TpTask[]            tp_tasks,        // will use // modify to have offsets for 8 cameras
    			(use_master? geometryCorrection_main: null), // final GeometryCorrection  geometryCorrection_main,
    			(use_aux?    geometryCorrection_aux: null),  // final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
    			ers_delay,                                   // final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
    			threadsMax,                                  // final int                 threadsMax,  // maximal number of threads to launch
    			debugLevel);                                 // final int                 debugLevel)
    	return tp_tasks;
    }

// All data is already copied to GPU memory
    public void execConverCorrectTiles() {
        if (GPU_CONVERT_CORRECT_TILES_kernel == null)
        {
            IJ.showMessage("Error", "No GPU kernel: GPU_CONVERT_CORRECT_TILES_kernel");
            return;
        }
        // kernel parameters: pointer to pointers
        int [] GridFullWarps =    {(num_task_tiles + TILES_PER_BLOCK -1 )/TILES_PER_BLOCK, 1, 1};
        int [] ThreadsFullWarps = {THREADSX, TILES_PER_BLOCK, 1};
        Pointer kernelParameters = Pointer.to(
            Pointer.to(gpu_kernel_offsets),
            Pointer.to(gpu_kernels),
            Pointer.to(gpu_bayer),
            Pointer.to(gpu_tasks),
            Pointer.to(gpu_clt),
            Pointer.to(new int[] { mclt_stride }),
            Pointer.to(new int[] { num_task_tiles }),
            Pointer.to(new int[] { 7 }) // lpf_mask
        );

        cuCtxSynchronize();
        	// Call the kernel function
    	cuLaunchKernel(GPU_CONVERT_CORRECT_TILES_kernel,
    			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
    			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
    			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
    			kernelParameters, null);   // Kernel- and extra parameters
    	cuCtxSynchronize(); // remove later
    }

    public void execImcltRbg() {
    	if (GPU_IMCLT_RBG_kernel == null)
    	{
    		IJ.showMessage("Error", "No GPU kernel: GPU_IMCLT_RBG_kernel");
    		return;
    	}
    	int tilesX =  IMG_WIDTH / DTT_SIZE;
    	int tilesY =  IMG_HEIGHT / DTT_SIZE;
    	int [] ThreadsFullWarps = {IMCLT_THREADS_PER_TILE, IMCLT_TILES_PER_BLOCK, 1};
    	for (int ncam = 0; ncam < NUM_CAMS; ncam++) {
    		for (int color = 0; color < NUM_COLORS; color++) {
    			for (int v_offs = 0; v_offs < 2; v_offs++){
    				for (int h_offs = 0; h_offs < 2; h_offs++){
    					int tilesy_half = (tilesY + (v_offs ^ 1)) >> 1;
	    				int tilesx_half = (tilesX + (h_offs ^ 1)) >> 1;
						int tiles_in_pass = tilesy_half * tilesx_half;
						int [] GridFullWarps =    {(tiles_in_pass + IMCLT_TILES_PER_BLOCK-1) / IMCLT_TILES_PER_BLOCK,1,1};
						Pointer kernelParameters = Pointer.to(
								Pointer.to(gpu_clt_h[ncam]),
								Pointer.to(gpu_corr_images_h[ncam]),
								Pointer.to(new int[] { color }),
								Pointer.to(new int[] { v_offs }),
								Pointer.to(new int[] { h_offs }),
								Pointer.to(new int[] { imclt_stride }) // lpf_mask
								);
						cuCtxSynchronize();
						// Call the kernel function
						cuLaunchKernel(GPU_IMCLT_RBG_kernel,
								GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
								ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
								0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
								kernelParameters, null);   // Kernel- and extra parameters
    				}
    			}
    		}
    	}
    	cuCtxSynchronize();
    }

    public float [][] getRBG (int ncam){
        int height = (IMG_HEIGHT + DTT_SIZE);
        int width =  (IMG_WIDTH + DTT_SIZE);
        int rslt_img_size =      width * height;
        float [] cpu_corr_image = new float [ NUM_COLORS * rslt_img_size];
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
        copyD2H.Height =          3 * height; // /2;

        cuMemcpy2D(copyD2H); // run copy

        float [][] fimg = new float [NUM_COLORS][ rslt_img_size];
        for (int ncol = 0; ncol < NUM_COLORS; ncol++) {
        	System.arraycopy(cpu_corr_image, ncol*rslt_img_size, fimg[ncol], 0, rslt_img_size);
        }
        return fimg;
    }

//    private static CUfunction [] createFunctions(
      private CUfunction [] createFunctions(
            String sourceCode, String [] kernelNames) throws IOException
        {
    	CUfunction [] functions = new CUfunction [kernelNames.length];
        	boolean OK = false;
            // Use the NVRTC to create a program by compiling the source code
            nvrtcProgram program = new nvrtcProgram();
            nvrtcCreateProgram(
                program, sourceCode, null, 0, null, null);
            try {
            	nvrtcCompileProgram(program, 0, null);
            	OK = true;
        	} catch (Exception e) {
        		System.out.println("nvrtcCompileProgram() FAILED");
        	}
            // Compilation log with errors/warnongs
            String programLog[] = new String[1];
            nvrtcGetProgramLog(program, programLog);
            String log = programLog[0].trim();
            if (!log.isEmpty())
            {
                System.err.println("Program compilation log:\n" + log);
            }
            if (!OK) {
            	throw new IOException("Could not compile program");
            }

            // Get the PTX code of the compiled program (not the binary)
            String[] ptx = new String[1];
            nvrtcGetPTX(program, ptx);
            nvrtcDestroyProgram(program);

            // Create a CUDA module from the PTX code
//            CUmodule
            module = new CUmodule();
            cuModuleLoadData(module, ptx[0]);

            for (int i = 0; i < kernelNames.length; i++) {
            // Find the function in the source by name, get its pointer
            	functions[i] = new CUfunction();
            	cuModuleGetFunction(functions[i] , module, kernelNames[i]);
            }

            return functions;
        }

    static String readFileAsString(String path)
    {
    	byte[] encoded;
    	try {
    		encoded = Files.readAllBytes(Paths.get(path));
    	} catch (IOException e) {
    		return null;
    	}
    	return new String(encoded, StandardCharsets.UTF_8);
    }

	public void  getTileSubcamOffsets(
			final TpTask[]            tp_tasks,        // will use // modify to have offsets for 8 cameras
			final GeometryCorrection  geometryCorrection_main,
			final GeometryCorrection  geometryCorrection_aux, // if null, will only calculate offsets fro the main camera
			final double [][][]       ers_delay,        // if not null - fill with tile center acquisition delay
			final int                 threadsMax,  // maximal number of threads to launch
			final int                 debugLevel)
	{
		final int quad_main = (geometryCorrection_main != null)? NUM_CAMS:0;
		final int quad_aux =  (geometryCorrection_aux != null)? NUM_CAMS:0;


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

						centerX = tileX * DTT_SIZE + DTT_SIZE/2; //  - shiftX;
						centerY = tileY * DTT_SIZE + DTT_SIZE/2; //  - shiftY;
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

	public void setLpfRbg(
			float sigma_r,
			float sigma_b,
			float sigma_g)
	{
		int dct_size = DTT_SIZE;
		DttRad2 dtt = new DttRad2(dct_size);

		double [][] lpf_rbg = {
				dtt.dttt_iiie(setCltLpf(sigma_r)),
				dtt.dttt_iiie(setCltLpf(sigma_b)),
				dtt.dttt_iiie(setCltLpf(sigma_g))};
		int l = dct_size*dct_size;
		float []   lpf_flat = new float [3 * l];
		for (int i = 0; i < 3;i++) {
//			System.arraycopy(lpf_rbg[i], 0, lpf_flat, l* i, l);
			for (int j = 0; j < l; j++) {
				lpf_flat[j+i*l] = (float) (lpf_rbg[i][j]*2*dct_size);
			}
		}

		CUdeviceptr constantMemoryPointer = new CUdeviceptr();
		long constantMemorySizeArray[] = { 0 };
		cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray,  module, "lpf_data");
		int constantMemorySize = (int)constantMemorySizeArray[0];
//__constant__ float lpf_data[3][64]={

		System.out.println("constantMemoryPointer: " + constantMemoryPointer);
		System.out.println("constantMemorySize: " + constantMemorySize);
        cuMemcpyHtoD(constantMemoryPointer, Pointer.to(lpf_flat), constantMemorySize);
		System.out.println();

	}

	public double [] setCltLpf(
			double   sigma)
	{
			int dct_size = DTT_SIZE;

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






}
