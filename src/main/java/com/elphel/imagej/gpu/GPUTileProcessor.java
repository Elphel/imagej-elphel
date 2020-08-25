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

import static jcuda.driver.CUdevice_attribute.CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR;
import static jcuda.driver.CUdevice_attribute.CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR;
// Uses code by Marco Hutter - http://www.jcuda.org
import static jcuda.driver.CUjitInputType.CU_JIT_INPUT_LIBRARY;
import static jcuda.driver.CUjitInputType.CU_JIT_INPUT_PTX;
import static jcuda.driver.CUjit_option.CU_JIT_LOG_VERBOSE;
import static jcuda.driver.JCudaDriver.cuCtxCreate;
import static jcuda.driver.JCudaDriver.cuCtxSynchronize;
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuDeviceGetAttribute;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuLaunchKernel;
import static jcuda.driver.JCudaDriver.cuLinkAddData;
import static jcuda.driver.JCudaDriver.cuLinkAddFile;
import static jcuda.driver.JCudaDriver.cuLinkComplete;
import static jcuda.driver.JCudaDriver.cuLinkCreate;
import static jcuda.driver.JCudaDriver.cuLinkDestroy;
import static jcuda.driver.JCudaDriver.cuMemAlloc;
import static jcuda.driver.JCudaDriver.cuMemAllocPitch;
import static jcuda.driver.JCudaDriver.cuMemcpy2D;
import static jcuda.driver.JCudaDriver.cuMemcpyDtoH;
import static jcuda.driver.JCudaDriver.cuMemcpyHtoD;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleGetGlobal;
import static jcuda.driver.JCudaDriver.cuModuleLoadDataEx;
import static jcuda.nvrtc.JNvrtc.nvrtcCompileProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcCreateProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcDestroyProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcGetPTX;
import static jcuda.nvrtc.JNvrtc.nvrtcGetProgramLog;

import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import com.elphel.imagej.tileprocessor.DttRad2;
import com.elphel.imagej.tileprocessor.GeometryCorrection;
import com.elphel.imagej.tileprocessor.ImageDtt;
import com.elphel.imagej.tileprocessor.QuadCLT;

import Jama.Matrix;
import ij.IJ;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.driver.CUDA_MEMCPY2D;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUlinkState;
import jcuda.driver.CUmemorytype;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import jcuda.driver.JITOptions;
import jcuda.nvrtc.JNvrtc;
import jcuda.nvrtc.nvrtcProgram;

public class GPUTileProcessor {
	String LIBRARY_PATH = "/usr/local/cuda/targets/x86_64-linux/lib/libcudadevrt.a"; // linux
	static String GPU_RESOURCE_DIR =              "kernels";
	static String [] GPU_KERNEL_FILES = {"dtt8x8.cuh","TileProcessor.cuh"};
	// "*" - generated defines, first index - separately compiled unit
	static String [][] GPU_SRC_FILES = {{"*","dtt8x8.h","dtt8x8.cu","geometry_correction.h","geometry_correction.cu","TileProcessor.h","TileProcessor.cuh"}};
	static String GPU_CONVERT_DIRECT_NAME =        "convert_direct"; // name in C code
	static String GPU_IMCLT_ALL_NAME =             "imclt_rbg_all";
	static String GPU_CORRELATE2D_NAME =           "correlate2D"; // name in C code
	static String GPU_TEXTURES_NAME =              "textures_nonoverlap"; // name in C code
	static String GPU_RBGA_NAME =                  "generate_RBGA"; // name in C code
	static String GPU_ROT_DERIV =                  "calc_rot_deriv"; // calculate rotation matrices and derivatives
	static String GPU_SET_TILES_OFFSETS =          "get_tiles_offsets"; // calculate pixel offsets and disparity distortions
	static String GPU_CALC_REVERSE_DISTORTION =    "calcReverseDistortionTable"; // calculate reverse radial distortion table from gpu_geometry_correction

//  pass some defines to gpu source code with #ifdef JCUDA
	public static int DTT_SIZE_LOG2 =             3;
	public static int DTT_SIZE =                  (1 << DTT_SIZE_LOG2);
	static int        THREADSX =                  DTT_SIZE;
	public static int NUM_CAMS =                  4;
	public static int NUM_PAIRS =                 6; // top hor, bottom hor, left vert, right vert, main diagonal, other diagonal
	public static int NUM_COLORS =                3;
//	public static int IMG_WIDTH =              2592;
//	public static int IMG_HEIGHT =             1936;
	static int        KERNELS_HOR =             164;
	static int        KERNELS_VERT =            123;
	static int        KERNELS_LSTEP =             4;
	static int        THREADS_PER_TILE =          8;
	static int        TILES_PER_BLOCK =           4; // 8 - slower
	static int        CORR_THREADS_PER_TILE =     8;
	static int        CORR_TILES_PER_BLOCK	=     4;
	static int        TEXTURE_THREADS_PER_TILE =  8; // 16;
	static int        TEXTURE_TILES_PER_BLOCK =   1;
	static int        IMCLT_THREADS_PER_TILE =   16;
	static int        IMCLT_TILES_PER_BLOCK =     4;
	static int        TPTASK_SIZE =                 1+ 1+ NUM_CAMS * 2 + 1 + NUM_CAMS * 4 ; // tp_task structure size in floats
	static int        CLTEXTRA_SIZE =               8;
	static int        CORR_SIZE =                   (2* DTT_SIZE - 1) * (2* DTT_SIZE - 1); // 15x15
	public static int CORR_NTILE_SHIFT =          8;  // also for texture tiles list
	public static int CORR_PAIRS_MASK =        0x3f;  // lower bits used to address correlation pair for the selected tile
	public static int CORR_TEXTURE_BIT =          7;  // bit 7 used to request texture for the tile
	public static int TASK_CORR_BITS =            4;  // start of pair mask
	public static int TASK_TEXTURE_N_BIT =        0; // Texture with North neighbor
	public static int TASK_TEXTURE_E_BIT =        1; // Texture with East  neighbor
	public static int TASK_TEXTURE_S_BIT =        2; // Texture with South neighbor
	public static int TASK_TEXTURE_W_BIT =        3; // Texture with West  neighbor
//	public static int TASK_TEXTURE_BIT =          3;  // bit to request texture calculation int task field of struct tp_task
	public static int LIST_TEXTURE_BIT =          7;  // bit to request texture calculation
	public static int CORR_OUT_RAD =              4;  // output radius of the correlations (implemented)
	public static double FAT_ZERO_WEIGHT =        0.0001; // add to port weights to avoid nan

	public static int THREADS_DYNAMIC_BITS =      5; // treads in block for CDP creation of the texture list

	public static int RBYRDIST_LEN =           5001; //for double, 10001 - float;   // length of rByRDist to allocate shared memory
	public static double RBYRDIST_STEP =          0.0004; //  for double, 0.0002 - for float; // to fit into GPU shared memory (was 0.001);
	public static int TILES_PER_BLOCK_GEOM =     32/NUM_CAMS; // blockDim.x = NUM_CAMS; blockDim.x = TILES_PER_BLOCK_GEOM
	public static int TASK_TEXTURE_BITS = ((1 << TASK_TEXTURE_N_BIT) | (1 << TASK_TEXTURE_E_BIT) | (1 << TASK_TEXTURE_S_BIT) | (1 << TASK_TEXTURE_W_BIT));


    int DTTTEST_BLOCK_WIDTH =        32; // may be read from the source code
    int DTTTEST_BLOCK_HEIGHT =       16; // may be read from the source code

    public boolean kernels_set = false;
    public boolean bayer_set =   false;

    private CUfunction GPU_CONVERT_DIRECT_kernel =          null;
    private CUfunction GPU_IMCLT_ALL_kernel =               null;
    private CUfunction GPU_CORRELATE2D_kernel =             null;
    private CUfunction GPU_TEXTURES_kernel =                null;
    private CUfunction GPU_RBGA_kernel =                    null;
    private CUfunction GPU_ROT_DERIV_kernel =               null;
    private CUfunction GPU_SET_TILES_OFFSETS_kernel =       null;
    private CUfunction GPU_CALC_REVERSE_DISTORTION_kernel = null;

    CUmodule    module; // to access constants memory
    public class TpTask {
    	public int   task; // [0](+1) - generate 4 images, [4..9]+16..+512 - correlation pairs, 2 - generate texture tiles
    	public float target_disparity;

    	public int ty;
    	public int tx;
    	public float[][] xy = null;
    	public float[][] xy_aux = null;
    	public float [][] disp_dist = null;
    	public TpTask() {}

    	public TpTask(int tx, int ty, float target_disparity, int task ) {
    		this.tx = tx;
    		this.ty = ty;
    		this.target_disparity = target_disparity;
    		this.task = task;
    		this.disp_dist = new float [NUM_CAMS][4];
    	}

    	// convert this class instance to float array to match layout of the C struct
    	public float [] asFloatArray(boolean use_aux) {
    		float [] flt = new float [TPTASK_SIZE];
    		return asFloatArray(flt, 0, use_aux);
    	}
    	// convert this class instance to float array to match layout of the C struct,
    	// fill existing float array from the specified index
    	public float [] asFloatArray(float [] flt, int indx, boolean use_aux) {
    		flt[indx++] = Float.intBitsToFloat(task);
    		flt[indx++] = Float.intBitsToFloat(tx + (ty << 16));
    		float [][] offsets = use_aux? this.xy_aux: this.xy;
    		for (int i = 0; i < NUM_CAMS; i++) {
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

    private String getTpDefines() {
        return"#define JCUDA\n"+
        				"#define DTT_SIZE_LOG2 " +            DTT_SIZE_LOG2+"\n"+
        				"#define THREADSX " +                 THREADSX+"\n"+
        				"#define NUM_CAMS " +                 NUM_CAMS+"\n"+
        				"#define NUM_PAIRS " +                NUM_PAIRS+"\n"+
        				"#define NUM_COLORS " +               NUM_COLORS+"\n"+
        				"#define KERNELS_LSTEP " +            KERNELS_LSTEP+"\n"+
        				"#define THREADS_PER_TILE " +         THREADS_PER_TILE+"\n"+
        				"#define TILES_PER_BLOCK " +          TILES_PER_BLOCK+"\n"+
        				"#define CORR_THREADS_PER_TILE " +    CORR_THREADS_PER_TILE+"\n"+
        				"#define CORR_TILES_PER_BLOCK " +     CORR_TILES_PER_BLOCK+"\n"+
        				"#define TEXTURE_THREADS_PER_TILE " + TEXTURE_THREADS_PER_TILE+"\n"+
        				"#define TEXTURE_TILES_PER_BLOCK " +  TEXTURE_TILES_PER_BLOCK+"\n"+
        				"#define IMCLT_THREADS_PER_TILE " +   IMCLT_THREADS_PER_TILE+"\n"+
        				"#define IMCLT_TILES_PER_BLOCK " +    IMCLT_TILES_PER_BLOCK+"\n"+
        				"#define CORR_NTILE_SHIFT " +         CORR_NTILE_SHIFT+"\n"+
        				"#define CORR_PAIRS_MASK " +          CORR_PAIRS_MASK+"\n"+
        				"#define CORR_TEXTURE_BIT " +         CORR_TEXTURE_BIT+"\n"+
        				"#define TASK_CORR_BITS " +           TASK_CORR_BITS+"\n"+
        				"#define TASK_TEXTURE_N_BIT " +       TASK_TEXTURE_N_BIT+"\n"+
        				"#define TASK_TEXTURE_E_BIT " +       TASK_TEXTURE_E_BIT+"\n"+
        				"#define TASK_TEXTURE_S_BIT " +       TASK_TEXTURE_S_BIT+"\n"+
        				"#define TASK_TEXTURE_W_BIT " +       TASK_TEXTURE_W_BIT+"\n"+
        				"#define LIST_TEXTURE_BIT " +         LIST_TEXTURE_BIT+"\n"+
        				"#define CORR_OUT_RAD " +             CORR_OUT_RAD+"\n" +
        				"#define FAT_ZERO_WEIGHT " +          FAT_ZERO_WEIGHT+"\n"+
        				"#define THREADS_DYNAMIC_BITS " +     THREADS_DYNAMIC_BITS+"\n"+
        				"#define RBYRDIST_LEN " +             RBYRDIST_LEN+"\n"+
        				"#define RBYRDIST_STEP " +            RBYRDIST_STEP+"\n"+
        				"#define TILES_PER_BLOCK_GEOM " +     TILES_PER_BLOCK_GEOM+"\n";


    }

    public GPUTileProcessor(
    		String cuda_project_directory) throws IOException
    {

    	// From code by Marco Hutter - http://www.jcuda.org
        // Enable exceptions and omit all subsequent error checks
        JCudaDriver.setExceptionsEnabled(true);
        JNvrtc.setExceptionsEnabled(true);

        // Initialize the driver and create a context for the first device.
        cuInit(0);
        //2020 - making them global
        CUdevice
        device = new CUdevice();
        cuDeviceGet(device, 0);
        CUcontext
        context = new CUcontext();
        cuCtxCreate(context, 0, device);

        int majorArray[] = { 0 };
        int minorArray[] = { 0 };
        cuDeviceGetAttribute(majorArray,  CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device);
        cuDeviceGetAttribute(minorArray,  CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device);
        int major = majorArray[0];
        int minor = minorArray[0];
        int capability = major * 10 + minor;

        // Obtain the CUDA source code from the CUDA file
        // Get absolute path to the file in resource folder, then read it as a normal file.
        // When using just Eclipse resources - it does not notice that the file
        // was edited (happens frequently during kernel development).
        ClassLoader classLoader = getClass().getClassLoader();

        String [] kernelSources = new String[GPU_SRC_FILES.length];

        for (int cunit = 0; cunit < kernelSources.length; cunit++) {
        	kernelSources[cunit] = ""; // use StringBuffer?
            for (String src_file:GPU_SRC_FILES[cunit]) {
            	if (src_file.contentEquals("*")) {
            		kernelSources[cunit] += getTpDefines();
            	}else {
                	File file = null;
                	if ((cuda_project_directory == null) || cuda_project_directory.isEmpty()) {
                		file = new File(classLoader.getResource(GPU_RESOURCE_DIR+"/"+src_file).getFile());
                		System.out.println("Loading resource "+file);
                	} else {
                		File src_dir = new File(cuda_project_directory, "src");
                		file = new File(src_dir.getPath(), src_file);
                		System.out.println("Loading resource "+file);
                	}
//                	System.out.println(file.getAbsolutePath());
                	String cuFileName = file.getAbsolutePath(); // /home/eyesis/workspace-python3/nvidia_dct8x8/src/dtt8x8.cuh";// "dtt8x8.cuh";
                	String sourceFile = readFileAsString(cuFileName); // readResourceAsString(cuFileName);
                	if (sourceFile == null) {
                		String msg = "Could not read the kernel source code from "+cuFileName;
                		IJ.showMessage("Error",	msg);
                		new IllegalArgumentException (msg);
                	}
                	kernelSources[cunit] += sourceFile;
            	}
            }
        }
        // Create the kernel functions (first - just test)
        String [] func_names = {
        		GPU_CONVERT_DIRECT_NAME,
        		GPU_IMCLT_ALL_NAME,
        		GPU_CORRELATE2D_NAME,
        		GPU_TEXTURES_NAME,
        		GPU_RBGA_NAME,
        		GPU_ROT_DERIV,
        		GPU_SET_TILES_OFFSETS,
        		GPU_CALC_REVERSE_DISTORTION
        };
        CUfunction[] functions = createFunctions(kernelSources,
        		                                 func_names,
        		                                 capability); // on my - 75

        GPU_CONVERT_DIRECT_kernel =          functions[0];
        GPU_IMCLT_ALL_kernel =               functions[1];
        GPU_CORRELATE2D_kernel =             functions[2];
        GPU_TEXTURES_kernel=                 functions[3];
        GPU_RBGA_kernel=                     functions[4];
        GPU_ROT_DERIV_kernel =               functions[5];
        GPU_SET_TILES_OFFSETS_kernel =       functions[6];
        GPU_CALC_REVERSE_DISTORTION_kernel = functions[7];

        System.out.println("GPU kernel functions initialized");
        System.out.println(GPU_CONVERT_DIRECT_kernel.toString());
        System.out.println(GPU_IMCLT_ALL_kernel.toString());
        System.out.println(GPU_CORRELATE2D_kernel.toString());
        System.out.println(GPU_TEXTURES_kernel.toString());
        System.out.println(GPU_RBGA_kernel.toString());
        System.out.println(GPU_ROT_DERIV_kernel.toString());
        System.out.println(GPU_SET_TILES_OFFSETS_kernel.toString());
        System.out.println(GPU_CALC_REVERSE_DISTORTION_kernel.toString());
        
        // GPU data structures are now initialized through GpuQuad instances
    }
    

    public static String [] getCorrTitles() {
    	return new String []{"hor-top","hor-bottom","vert-left","vert-right","diag-main","diag-other"};
    }
    public static double [][] getCorr2DView(
    		int tilesX,
    		int tilesY,
    		int [] indices,
    		float [][] corr2d,
    		int [] wh){ // if is [2] - return width, height
    	if ((corr2d == null) || (corr2d.length == 0)) {
    		return new double [NUM_PAIRS][0];
    	}

    	int corr_size = (int)(Math.round(Math.sqrt(corr2d[0].length)));//  make smaller later?
    	int width =  tilesX * (corr_size + 1) + 1;
    	int height = tilesY * (corr_size + 1) + 1;
    	double [][] data = new double [NUM_PAIRS][];
    	data[0] = new double[height*width];
    	for (int ty = 0; ty < tilesY; ty++) {
    		for (int tx = 0; tx < tilesX; tx++) {
    			for (int i = 0; i< corr_size; i++) {
    				for (int j = 0; j < corr_size; j++) {
    					data[0][(ty * (corr_size + 1) + i + 1) * width + (tx * (corr_size + 1) + j + 1)] = Double.NaN;
    				}
    			}
    		}
    	}
		for (int np = 1; np < NUM_PAIRS; np++) {
			data[np] = data[0].clone();
		}
		for (int n = 0; n < indices.length; n++) {
			int nt = indices[n] >> CORR_NTILE_SHIFT;
			int np = indices[n] & CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
			assert np < NUM_PAIRS : "invalid correllation pair";
			int tx = nt % tilesX;
			int ty = nt / tilesX;
			for (int i = 0; i< corr_size; i++) {
				for (int j = 0; j < corr_size; j++) {
					//java.lang.ArrayIndexOutOfBoundsException: 20081634
					int indx1 = (ty * (corr_size + 1) + i + 1) * width + (tx * (corr_size + 1) + j + 1);
					int indx2 = i*corr_size+j;
//    					if ((indx1 > data[0].length) || (indx1 > data[0].length)){
//    						System.out.println("Bugggg!)");
//    					}
					data[np][indx1] = corr2d[n][indx2];
				}
			}
		}
		if (wh != null) {
			wh[0] = width;
			wh[1] = height;
		}
    	return data;
    }



//    private static CUfunction [] createFunctions(
    private CUfunction [] createFunctions(
    		String []  sourceCodeUnits,
    		String []  kernelNames,
    		int        capability
    		) throws IOException
    {
    	CUfunction [] functions = new CUfunction [kernelNames.length];
    	byte[][] ptxDataUnits = new byte [sourceCodeUnits.length][];
    	boolean OK = false;
       	for (int cunit = 0; cunit < ptxDataUnits.length; cunit++) {
       		String sourceCode = sourceCodeUnits[cunit];
    		// Use the NVRTC to create a program by compiling the source code
    		nvrtcProgram program = new nvrtcProgram();
    		nvrtcCreateProgram(	program, sourceCode, null, 0, null, null);
    		String options[] = {"--gpu-architecture=compute_"+capability};

    		try {
    			nvrtcCompileProgram(program, options.length, options);
    			OK = true;
    		} catch (Exception e) {
    			System.out.println("nvrtcCompileProgram() FAILED");
    		}
    		// Compilation log with errors/warnings
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
    		ptxDataUnits[cunit] = ptx[0].getBytes();
    		System.out.println("ptxDataUnits["+cunit+"].length="+ptxDataUnits[cunit].length);
    	}
    	JITOptions jitOptions = new JITOptions();
    	jitOptions.putInt(CU_JIT_LOG_VERBOSE, 1);
    	CUlinkState state = new CUlinkState();
    	cuLinkCreate(jitOptions, state);
    	cuLinkAddFile(state, CU_JIT_INPUT_LIBRARY, LIBRARY_PATH, jitOptions);

       	for (int cunit = 0; cunit < ptxDataUnits.length; cunit++) {
       		cuLinkAddData(state, CU_JIT_INPUT_PTX,     Pointer.to(ptxDataUnits[cunit]), ptxDataUnits[cunit].length, "input"+cunit+".ptx", jitOptions); // CUDA_ERROR_INVALID_PTX
       	}
    	long size[] = { 0 };
    	Pointer image = new Pointer();
    	JCudaDriver.setExceptionsEnabled(false);
    	int cuda_result = cuLinkComplete(state, image, size);
    	System.out.println("cuLinkComplete() -> "+cuda_result);
    	JCudaDriver.setExceptionsEnabled(true);
    	module = new CUmodule();
    	cuModuleLoadDataEx(module, image, 0, new int[0], Pointer.to(new int[0]));
    	cuLinkDestroy(state);

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

   public class GpuQuad{ // quad camera description
        public final QuadCLT quadCLT;
    	public final int img_width;
    	public final int img_height;
    	public final int kernels_hor;        // int                kernels_hor,
		public final int kernels_vert;       // int                kernels_vert);
    	
    	public final int num_cams;
    	public final int num_colors; // maybe should always be 3?
    	public final int kern_tiles;
    	public final int kern_size;
    	
//    	public final GPUTileProcessor gPUTileProcessor;
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
        private CUdeviceptr gpu_tasks;
        private CUdeviceptr gpu_corrs;
        private CUdeviceptr gpu_textures;
        private CUdeviceptr gpu_clt;
        private CUdeviceptr gpu_4_images;
        private CUdeviceptr gpu_corr_indices;
        private CUdeviceptr gpu_num_corr_tiles;
        private CUdeviceptr gpu_texture_indices_ovlp;
        private CUdeviceptr gpu_num_texture_ovlp;
        private CUdeviceptr gpu_texture_indices;
        private CUdeviceptr gpu_texture_indices_len;
        private CUdeviceptr gpu_diff_rgb_combo;
        private CUdeviceptr gpu_color_weights;
        private CUdeviceptr gpu_generate_RBGA_params;
        private CUdeviceptr gpu_woi;
        private CUdeviceptr gpu_textures_rgba;
        private CUdeviceptr gpu_correction_vector;
        private CUdeviceptr gpu_rot_deriv;
        private CUdeviceptr gpu_geometry_correction;
        private CUdeviceptr gpu_rByRDist;
        private CUdeviceptr gpu_active_tiles;
        private CUdeviceptr gpu_num_active_tiles;
        private int mclt_stride;
        private int corr_stride;
        private int imclt_stride;
        private int texture_stride;
        private int texture_stride_rgba;
        private int num_task_tiles;
        private int num_corr_tiles;
        private int num_texture_tiles;
        private boolean geometry_correction_set = false;
        private boolean geometry_correction_vector_set = false;
    	public GpuQuad(
        	final QuadCLT quadCLT,    			
//   			final int img_width,
//   			final int img_height,
//   	    	final int kernels_hor,
//   			final int kernels_vert,
   			final int num_cams,
   			final int num_colors
) {
    		this.quadCLT =      quadCLT;
    		int [] wh =         quadCLT.getGeometryCorrection().getSensorWH();
        	this.img_width =    wh[0]; // img_width;
        	this.img_height =   wh[1]; // img_height;
        	this.num_cams =     num_cams;
        	this.num_colors =   num_colors; // maybe should always be 3?
        	this.kernels_hor =  quadCLT.getCLTKernels()[0][0][0].length; // kernels_hor;
        	this.kernels_vert = quadCLT.getCLTKernels()[0][0].length; // kernels_vert;
        	this.kern_tiles =   kernels_hor *  kernels_vert * num_colors;
        	this.kern_size =    kern_tiles * 4 * 64;
        	
        	
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
            gpu_tasks =               new CUdeviceptr(); //  allocate tilesX * tilesY * TPTASK_SIZE * Sizeof.FLOAT
            gpu_corrs =               new CUdeviceptr(); //  allocate tilesX * tilesY * NUM_PAIRS * CORR_SIZE * Sizeof.FLOAT
            gpu_textures =            new CUdeviceptr(); //  allocate tilesX * tilesY * ? * 256 * Sizeof.FLOAT
            gpu_clt =                 new CUdeviceptr();
            gpu_4_images =            new CUdeviceptr();
            gpu_corr_indices =        new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
            gpu_num_corr_tiles =      new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
            gpu_texture_indices_ovlp =new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
            gpu_num_texture_ovlp =    new CUdeviceptr(); //  8 ints
            gpu_texture_indices =     new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
            gpu_texture_indices_len = new CUdeviceptr(); //  allocate tilesX * tilesY * 6 * Sizeof.FLOAT
            gpu_diff_rgb_combo =      new CUdeviceptr(); //  1 int

            gpu_color_weights =       new CUdeviceptr(); //  allocate 3 * Sizeof.FLOAT
            gpu_generate_RBGA_params =new CUdeviceptr(); //  allocate 5 * Sizeof.FLOAT

            gpu_woi =                 new CUdeviceptr(); //  4 integers (x, y, width, height) Rectangle - in tiles
            gpu_textures_rgba =       new CUdeviceptr(); //  allocate tilesX * tilesY * ? * 256 * Sizeof.FLOAT

            gpu_correction_vector=    new CUdeviceptr();
            gpu_rot_deriv=            new CUdeviceptr(); //  used internally by device, may be read to CPU for testing
            gpu_geometry_correction=  new CUdeviceptr();
            gpu_rByRDist=             new CUdeviceptr(); //  calculated once for the camera distortion model in CPU (move to GPU?)

            gpu_active_tiles =        new CUdeviceptr(); //  TILESX*TILESY*sizeof(int)
            gpu_num_active_tiles =    new CUdeviceptr(); //  1 int
        	
            // Init data arrays for all kernels
            int tilesX =  img_width / DTT_SIZE;
            int tilesY =  img_height / DTT_SIZE;
            long [] device_stride = new long [1];
            for (int ncam = 0; ncam < num_cams; ncam++) {
            	gpu_kernels_h[ncam] =        new CUdeviceptr();
            	cuMemAlloc(gpu_kernels_h[ncam],kern_size * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
            	gpu_kernel_offsets_h[ncam] = new CUdeviceptr();
            	cuMemAlloc(gpu_kernel_offsets_h[ncam],kern_tiles * CLTEXTRA_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
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
                		(img_width + DTT_SIZE) * Sizeof.FLOAT, // long WidthInBytes,
                		3*(img_height + DTT_SIZE),// long Height,
                        Sizeof.FLOAT);            // int ElementSizeBytes)
                imclt_stride = (int)(device_stride[0] / Sizeof.FLOAT);
                gpu_clt_h[ncam] = new CUdeviceptr();
            	cuMemAlloc(gpu_clt_h[ncam],tilesY * tilesX * num_colors * 4 * DTT_SIZE * DTT_SIZE * Sizeof.FLOAT ); //     public static int cuMemAlloc(CUdeviceptr dptr, long bytesize)
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

        	for (int ncam = 0; ncam < num_cams; ncam++) gpu_kernels_l[ncam] =        getPointerAddress(gpu_kernels_h[ncam]);
            cuMemcpyHtoD(gpu_kernels, Pointer.to(gpu_kernels_l),                     num_cams * Sizeof.POINTER);

            for (int ncam = 0; ncam < num_cams; ncam++) gpu_kernel_offsets_l[ncam] = getPointerAddress(gpu_kernel_offsets_h[ncam]);
            cuMemcpyHtoD(gpu_kernel_offsets, Pointer.to(gpu_kernel_offsets_l),       num_cams * Sizeof.POINTER);

            for (int ncam = 0; ncam < num_cams; ncam++) gpu_bayer_l[ncam] =          getPointerAddress(gpu_bayer_h[ncam]);
            cuMemcpyHtoD(gpu_bayer, Pointer.to(gpu_bayer_l),                         num_cams * Sizeof.POINTER);

            for (int ncam = 0; ncam < num_cams; ncam++) gpu_clt_l[ncam] =            getPointerAddress(gpu_clt_h[ncam]);
            cuMemcpyHtoD(gpu_clt, Pointer.to(gpu_clt_l),                             num_cams * Sizeof.POINTER);

            for (int ncam = 0; ncam < num_cams; ncam++) gpu_4_images_l[ncam] =       getPointerAddress(gpu_corr_images_h[ncam]);
            cuMemcpyHtoD(gpu_4_images, Pointer.to(gpu_4_images_l),                   num_cams * Sizeof.POINTER);

            // Set GeometryCorrection data
        	cuMemAlloc(gpu_geometry_correction,      GeometryCorrection.arrayLength(num_cams) * Sizeof.FLOAT);
        	cuMemAlloc(gpu_rByRDist,                 RBYRDIST_LEN *  Sizeof.FLOAT);
        	cuMemAlloc(gpu_rot_deriv,                5*num_cams*3*3 * Sizeof.FLOAT);
        	cuMemAlloc(gpu_correction_vector,        GeometryCorrection.CorrVector.LENGTH * Sizeof.FLOAT);

            // Set task array
        	cuMemAlloc(gpu_tasks,      tilesX * tilesY * TPTASK_SIZE * Sizeof.FLOAT);
    //=========== Seems that in many places Sizeof.POINTER (==8) is used instead of Sizeof.FLOAT !!! ============
        	// Set corrs array
        	cuMemAlloc(gpu_corr_indices,   tilesX * tilesY * NUM_PAIRS * Sizeof.FLOAT);
        	cuMemAlloc(gpu_num_corr_tiles,                    1 * Sizeof.FLOAT);

        	//#define TILESYA       ((TILESY +3) & (~3))
        	int tilesYa = (tilesY + 3) & ~3;
        	cuMemAlloc(gpu_texture_indices,     tilesX * tilesYa * Sizeof.FLOAT); // for non-overlap tiles
        	cuMemAlloc(gpu_texture_indices_ovlp,tilesX * tilesYa * Sizeof.FLOAT); // for overlapped tiles

        	cuMemAlloc(gpu_diff_rgb_combo, tilesX * tilesYa * num_cams* (num_colors + 1) *  Sizeof.FLOAT);

        	cuMemAlloc(gpu_color_weights,             3 * Sizeof.FLOAT);
        	cuMemAlloc(gpu_generate_RBGA_params,      5 * Sizeof.FLOAT);


        	cuMemAlloc(gpu_woi,                               4 * Sizeof.FLOAT);
        	cuMemAlloc(gpu_num_texture_ovlp,                  8 * Sizeof.FLOAT);
        	cuMemAlloc(gpu_texture_indices_len,               1 * Sizeof.FLOAT);


        	cuMemAlloc(gpu_active_tiles,        tilesX * tilesY * Sizeof.FLOAT);
        	cuMemAlloc(gpu_num_active_tiles,                  1 * Sizeof.FLOAT);

            cuMemAllocPitch (
            		gpu_corrs,                             // CUdeviceptr dptr,
            		device_stride,                         // long[] pPitch,
            		CORR_SIZE * Sizeof.FLOAT,              // long WidthInBytes,
            		NUM_PAIRS * tilesX * tilesY,             // long Height,
                    Sizeof.FLOAT);                         // int ElementSizeBytes)
            corr_stride = (int)(device_stride[0] / Sizeof.FLOAT);
            int max_texture_size = (num_colors + 1 + (num_cams + num_colors + 1)) * (2 * DTT_SIZE)* (2 * DTT_SIZE);
            cuMemAllocPitch (
            		gpu_textures,                             // CUdeviceptr dptr,
            		device_stride,                         // long[] pPitch,
            		max_texture_size * Sizeof.FLOAT,              // long WidthInBytes,
            		tilesX * tilesY,             // long Height,
                    Sizeof.FLOAT);                         // int ElementSizeBytes)
            texture_stride = (int)(device_stride[0] / Sizeof.FLOAT);
            int max_rgba_width  =  (tilesX + 1) * DTT_SIZE;
            int max_rgba_height =  (tilesY + 1) * DTT_SIZE;
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
            return img_width / DTT_SIZE;
    	}

    	public int getTilesY() {
            return img_height / DTT_SIZE;
    	}
    	
    	public void resetGeometryCorrection() {
    		geometry_correction_set = false;
    		geometry_correction_vector_set = false;
    	}
    	public void resetGeometryCorrectionVector() {
    		geometry_correction_vector_set = false;
    	}
    	public int getImageWidth()  {return this.img_width;}
    	public int getImageHeight() {return this.img_height;}
        public int getDttSize()     {return DTT_SIZE;}
        public int getNumCams()     {return NUM_CAMS;}
        
    	public void setGeometryCorrection() { // will reset geometry_correction_set when running GPU kernel
//    		if (geometry_correction_set) return;
    		setGeometryCorrection(
    				quadCLT.getGeometryCorrection(),
            		false);
    	}

    	public void setGeometryCorrectionVector() { // will reset geometry_correction_vector_set when running GPU kernel
    		setExtrinsicsVector(
    				quadCLT.getGeometryCorrection().getCorrVector());
    	}
    	
        public void setGeometryCorrection(GeometryCorrection gc,
        		boolean use_java_rByRDist) { // false - use newer GPU execCalcReverseDistortions
        	float [] fgc = gc.toFloatArray();
        	if (use_java_rByRDist) {
        		double [] rByRDist = gc.getRByRDist();
        		float [] fFByRDist = new float [rByRDist.length];
        		for (int i = 0; i < rByRDist.length; i++) {
        			fFByRDist[i] = (float) rByRDist[i];
        		}
        		cuMemcpyHtoD(gpu_rByRDist,            Pointer.to(fFByRDist), fFByRDist.length * Sizeof.FLOAT);
        	}
        	cuMemcpyHtoD(gpu_geometry_correction, Pointer.to(fgc),       fgc.length * Sizeof.FLOAT);
        	cuMemAlloc  (gpu_rot_deriv, 5 * num_cams *3 *3 * Sizeof.FLOAT); // NCAM of 3x3 rotation matrices, plus 4 derivative matrices for each camera
        }
/**
 * Copy extrinsic correction vector to the GPU memory
 * @param cv correction vector
 */
        public void setExtrinsicsVector(GeometryCorrection.CorrVector cv) {
        	double [] dcv = cv.toFullRollArray();
        	float []  fcv = new float [dcv.length];
        	for (int i = 0; i < dcv.length; i++) {
        		fcv[i] = (float) dcv[i];
        	}
        	cuMemcpyHtoD(gpu_correction_vector, Pointer.to(fcv), fcv.length * Sizeof.FLOAT);
        }

/**
 * Copy array of CPU-prepared tasks to the GPU memory
 * @param tile_tasks array of TpTask prepared by the CPU (before geometry correction is applied)
 * @param use_aux Use second (aux) camera
 */
        public void setTasks(TpTask [] tile_tasks, boolean use_aux) // while is it in class member? - just to be able to free
        {
        	num_task_tiles = tile_tasks.length;
        	float [] ftasks = new float [TPTASK_SIZE * num_task_tiles];
        	for (int i = 0; i < num_task_tiles; i++) {
        		tile_tasks[i].asFloatArray(ftasks, i* TPTASK_SIZE, use_aux);
        	}
            cuMemcpyHtoD(gpu_tasks, Pointer.to(ftasks), TPTASK_SIZE * num_task_tiles * Sizeof.FLOAT);
        }

        public void setCorrIndices(int [] corr_indices)
        {
        	num_corr_tiles = corr_indices.length;
        	float [] fcorr_indices = new float [corr_indices.length];
        	for (int i = 0; i < num_corr_tiles; i++) {
        		fcorr_indices[i] = Float.intBitsToFloat(corr_indices[i]);
        	}
            cuMemcpyHtoD(gpu_corr_indices, Pointer.to(fcorr_indices),  num_corr_tiles * Sizeof.FLOAT);
        }

        public void setTextureIndices(int [] texture_indices) // never used
        {
        	num_texture_tiles = texture_indices.length;
        	float [] ftexture_indices = new float [texture_indices.length];
        	for (int i = 0; i < num_texture_tiles; i++) {
        		ftexture_indices[i] = Float.intBitsToFloat(texture_indices[i]);
        	}
            cuMemcpyHtoD(gpu_texture_indices, Pointer.to(ftexture_indices),  num_texture_tiles * Sizeof.FLOAT);
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

    //texture_indices
        public void setConvolutionKernel(
        		float [] kernel,  // [tileY][tileX][color][..]
        		float [] kernel_offsets,
        		int ncam) {
            cuMemcpyHtoD(gpu_kernels_h[ncam],        Pointer.to(kernel),         kern_size * Sizeof.FLOAT);
            cuMemcpyHtoD(gpu_kernel_offsets_h[ncam], Pointer.to(kernel_offsets), kern_tiles * CLTEXTRA_SIZE * Sizeof.FLOAT);
        }

        /**
         * Set CLT kernels in GPU memory reading them using quadCLT instance        
         * @param clt_kernels - [num_cameras==4][num_colors][tilesY][tilesX][quadrants==4][elements==8*8==64]
         * @param force re-set kernels even if they were already set
         */
        public void setConvolutionKernels(
        		boolean force)
        {
        	if (kernels_set && ! force) {
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
        	if (!force && bayer_set && !quadCLT.hasNewImageData()) {
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
        	if (bayer_set && !force) {
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
    		bayer_set = true;
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
            int tilesX =  img_width / DTT_SIZE;
            int tilesY =  img_height / DTT_SIZE;
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
            int tilesX =  img_width / DTT_SIZE;
            int tilesY =  img_height / DTT_SIZE;
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
        				tp_tasks[indx] = new TpTask(tx,ty, target_disparities[indx]+disparity_corr,
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
            int tilesX =  img_width / DTT_SIZE;
            int tilesY =  img_height / DTT_SIZE;
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
            		tp_tasks[indx] = new TpTask(tx,ty, target_disparities[indx],
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
       
    	public GPUTileProcessor.TpTask[]  setTpTask(
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
    		final List<GPUTileProcessor.TpTask> task_list = new CopyOnWriteArrayList<GPUTileProcessor.TpTask>();
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
    						task_list.add(new GPUTileProcessor.TpTask(
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
    		return task_list.toArray(new GPUTileProcessor.TpTask[task_list.size()]);		
    	}

        
        
        /**
         * Prepare contents pointers for calculation of the correlation pairs
         * @param tp_tasks array of tasks that contain masks of the required pairs
         * @return each element has (tile_number << 8) | (pair_number & 0xff)
         */
        public int [] getCorrTasks(
        		TpTask [] tp_tasks) {
        	int tilesX = img_width / DTT_SIZE;
        	int num_corr = 0;
        	int task_mask = (1 << NUM_PAIRS) - 1;
        	for (TpTask tt: tp_tasks) {
        		int pm = (tt.task >> TASK_CORR_BITS) & task_mask;
        		if (pm != 0) {
        			for (int b = 0; b < NUM_PAIRS; b++) if ((pm & (1 << b)) != 0) {
        				num_corr++;    			}
        		}
        	}

        	int [] iarr = new int[num_corr];
        	num_corr = 0;
        	for (TpTask tt: tp_tasks) {
        		int pm = (tt.task >> TASK_CORR_BITS) & task_mask;
        		if (pm != 0) {
        			int tile = (tt.ty * tilesX +tt.tx);
        			for (int b = 0; b < NUM_PAIRS; b++) if ((pm & (1 << b)) != 0) {
        				iarr[num_corr++] = (tile << CORR_NTILE_SHIFT) | b;
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
        public int [] getTextureTasks(
        		TpTask [] tp_tasks) {
        	int tilesX = img_width / DTT_SIZE;
        	int num_textures = 0;
        	for (TpTask tt: tp_tasks) {
        		if ((tt.task & TASK_TEXTURE_BITS) !=0) {
        			num_textures++;
        		}
        	}
        	int [] iarr = new int[num_textures];
        	num_textures = 0;
        	int b = (1 << LIST_TEXTURE_BIT);
        	for (TpTask tt: tp_tasks) {
        		if ((tt.task & TASK_TEXTURE_BITS) !=0) {
        			int tile = (tt.ty * tilesX +tt.tx);
        			iarr[num_textures++] = (tile << CORR_NTILE_SHIFT) | b;
        		}
        	}
        	return iarr;
        }

/**
 * Calculate rotation matrices and their derivatives
 * Needed after extrinsics changed
 */
        public void execRotDerivs() {
            if (geometry_correction_vector_set)	return;
            if (GPU_ROT_DERIV_kernel == null)
            {
                IJ.showMessage("Error", "No GPU kernel: GPU_ROT_DERIV_kernel");
                return;
            }
            setGeometryCorrectionVector();
            // kernel parameters: pointer to pointers
            int [] GridFullWarps =    {num_cams, 1, 1}; // round up
            int [] ThreadsFullWarps = {3,        3, 3};
            Pointer kernelParameters = Pointer.to(
                Pointer.to(gpu_correction_vector),
                Pointer.to(gpu_rot_deriv)
            );

            cuCtxSynchronize();
            	// Call the kernel function
        	cuLaunchKernel(GPU_ROT_DERIV_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize(); // remove later
        	geometry_correction_vector_set = true;
        }
        
/**
 *  Calculate reverse radial distortion table from gpu_geometry_correction
 *  Needed once during initialization of GeometryCorrection       
 */
        public void execCalcReverseDistortions() {
            if (geometry_correction_set) return;
            if (GPU_CALC_REVERSE_DISTORTION_kernel == null)
            {
                IJ.showMessage("Error", "No GPU kernel: GPU_CALC_REVERSE_DISTORTION_kernel");
                return;
            }
            setGeometryCorrection();
            // kernel parameters: pointer to pointers
            int [] GridFullWarps =    {num_cams, 1, 1}; // round up
            int [] ThreadsFullWarps = {3,        3, 3};
            Pointer kernelParameters = Pointer.to(
            		Pointer.to(gpu_geometry_correction),     //	struct gc          * gpu_geometry_correction,
            		Pointer.to(gpu_rByRDist));                //	float *              gpu_rByRDist)      // length should match RBYRDIST_LEN
            cuCtxSynchronize();
            	// Call the kernel function
        	cuLaunchKernel(GPU_CALC_REVERSE_DISTORTION_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize(); // remove later
        	geometry_correction_set = true;
        }
        
/**
 * Calculate tiles offsets (before each direct conversion run)
 */
        public void execSetTilesOffsets() {
        	execCalcReverseDistortions(); // will check if it is needed first
        	execRotDerivs();              // will check if it is needed first
            if (GPU_SET_TILES_OFFSETS_kernel == null)
            {
                IJ.showMessage("Error", "No GPU kernel: GPU_SET_TILES_OFFSETS_kernel");
                return;
            }
            // kernel parameters: pointer to pointers
            int [] GridFullWarps =    {(num_task_tiles + TILES_PER_BLOCK_GEOM - 1)/TILES_PER_BLOCK_GEOM, 1, 1}; // round up
            int [] ThreadsFullWarps = {num_cams, TILES_PER_BLOCK_GEOM, 1}; // 4,8,1
            Pointer kernelParameters = Pointer.to(
            		Pointer.to(gpu_tasks),                   // struct tp_task     * gpu_tasks,
            		Pointer.to(new int[] { num_task_tiles }),// int                  num_tiles,          // number of tiles in task list
            		Pointer.to(gpu_geometry_correction),     //	struct gc          * gpu_geometry_correction,
            		Pointer.to(gpu_correction_vector),       //	struct corr_vector * gpu_correction_vector,
            		Pointer.to(gpu_rByRDist),                //	float *              gpu_rByRDist)      // length should match RBYRDIST_LEN
            		Pointer.to(gpu_rot_deriv));              // trot_deriv         * gpu_rot_deriv);
            cuCtxSynchronize();
        	cuLaunchKernel(GPU_SET_TILES_OFFSETS_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize(); // remove later
        }

/**
 * Direct CLT conversion and aberration correction 
 */
        public void execConvertDirect() {
            if (GPU_CONVERT_DIRECT_kernel == null)
            {
                IJ.showMessage("Error", "No GPU kernel: GPU_CONVERT_DIRECT_kernel");
                return;
            }
            setConvolutionKernels(false); // set kernels if they are not set already
            setBayerImages(false); // set Bayer images if this.quadCLT instance has new ones
            // kernel parameters: pointer to pointers
        	int tilesX =  img_width / DTT_SIZE;
            int [] GridFullWarps =    {1, 1, 1};
            int [] ThreadsFullWarps = {1, 1, 1};
            Pointer kernelParameters = Pointer.to(
            		Pointer.to(gpu_kernel_offsets),
            		Pointer.to(gpu_kernels),
            		Pointer.to(gpu_bayer),
            		Pointer.to(gpu_tasks),
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
        	cuLaunchKernel(GPU_CONVERT_DIRECT_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize(); // remove later
        }


/**
 * Generate corrected image(s) from the CLT representation created by execConvertDirect()
 * @param is_mono monochrome image (false - RGB)
 */
        public void execImcltRbgAll(
        		boolean is_mono
        		) {
        	if (GPU_IMCLT_ALL_kernel == null)
        	{
        		IJ.showMessage("Error", "No GPU kernel: GPU_IMCLT_ALL_kernel");
        		return;
        	}
        	int apply_lpf =  1;
        	int tilesX =  img_width / DTT_SIZE;
        	int tilesY =  img_height / DTT_SIZE;
        	int [] ThreadsFullWarps = {1, 1, 1};
        	int [] GridFullWarps =    {1, 1, 1};
        	Pointer kernelParameters = Pointer.to(
                    Pointer.to(gpu_clt),                                // float  ** gpu_clt, // [num_cams][TILESY][TILESX][num_colors][DTT_SIZE*DTT_SIZE]
                    Pointer.to(gpu_4_images),                           // float           ** gpu_corr_images,    // [num_cams][WIDTH, 3 * HEIGHT]
        			Pointer.to(new int[] { apply_lpf }),                // int                apply_lpf,
        			Pointer.to(new int[] { is_mono ? 1 : num_colors }), // int                colors,
        			Pointer.to(new int[] { tilesX }),                   // int                woi_twidth,
        			Pointer.to(new int[] { tilesY }),                   // int                woi_theight,
        			Pointer.to(new int[] { imclt_stride })              // const size_t       dstride);            // in floats (pixels)
        			);
        	cuCtxSynchronize();
        	// Call the kernel function
        	cuLaunchKernel(GPU_IMCLT_ALL_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                   // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize();
        }

/**
 * Generate 2D phase correlations from the CLT representation         
 * @param scales R,G,B weights
 * @param fat_zero - absolute fat zero - add to correlations before normalization
 * @param corr_radius - correlation result size (maximal 7 for 15x15)
 */
        public void execCorr2D(
        		double [] scales,
        		double fat_zero,
        		int corr_radius) {
        	if (GPU_CORRELATE2D_kernel == null)
        	{
        		IJ.showMessage("Error", "No GPU kernel: GPU_CORRELATE2D_kernel");
        		return;
        	}
            int tilesX =  img_width / DTT_SIZE;
        	int num_colors = scales.length;
        	if (num_colors > 3) num_colors = 3;
        	float fscale0 = (float) scales[0];
        	float fscale1 = (num_colors >1)?((float) scales[1]):0.0f;
        	float fscale2 = (num_colors >2)?((float) scales[2]):0.0f;
    		int [] GridFullWarps =    {1, 1, 1};
        	int [] ThreadsFullWarps = {1, 1, 1};
        	Pointer kernelParameters = Pointer.to(
        			Pointer.to(gpu_clt),                        // float          ** gpu_clt,
        			Pointer.to(new int[] { num_colors }),       // int               colors,             // number of colors (3/1)
        			Pointer.to(new float[] {fscale0  }),        // float             scale0,             // scale for R
        			Pointer.to(new float[] {fscale1  }),        // float             scale1,             // scale for B
        			Pointer.to(new float[] {fscale2  }),        // float             scale2,             // scale for G
        			Pointer.to(new float[] {(float) fat_zero }),// float             fat_zero,           // here - absolute
            		Pointer.to(gpu_tasks),                      // struct tp_task  * gpu_tasks,
            		Pointer.to(new int[] { num_task_tiles }),   // int               num_tiles           // number of tiles in task
            		Pointer.to(new int[] { tilesX }),           // int               tilesx,             // number of tile rows
        			Pointer.to(gpu_corr_indices),               // int             * gpu_corr_indices,   // packed tile+pair
        			Pointer.to(gpu_num_corr_tiles),             // int             * pnum_corr_tiles,    // pointer to a number of tiles to process
        			Pointer.to(new int[] { corr_stride }),      // const size_t      corr_stride,        // in floats
        			Pointer.to(new int[] { corr_radius }),      // int               corr_radius,        // radius of the output correlation (7 for 15x15)
        			Pointer.to(gpu_corrs)                       // float           * gpu_corrs);         // correlation output data
        			);
        	cuCtxSynchronize();
        	// Call the kernel function
        	cuLaunchKernel(GPU_CORRELATE2D_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize();
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
        public void execRBGA(
        		double [] color_weights,
        		boolean   is_lwir,
        		double    min_shot,           // 10.0
        		double    scale_shot,         // 3.0
        		double    diff_sigma,         // pixel value/pixel change
        		double    diff_threshold,     // pixel value/pixel change
        		double    min_agree,          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
        		boolean   dust_remove) {
        	execCalcReverseDistortions(); // will check if it is needed first
        	if (GPU_RBGA_kernel == null) {
        		IJ.showMessage("Error", "No GPU kernel: GPU_TEXTURES_kernel");
        		return;
        	}
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
        			(float)             min_agree          // minimal number of channels to agree on a point (real number to work with fuzzy averages)
        	};
        	cuMemcpyHtoD(gpu_generate_RBGA_params, Pointer.to(generate_RBGA_params),  generate_RBGA_params.length * Sizeof.FLOAT);

        	int iis_lwir =      (is_lwir)? 1:0;
        	int idust_remove =  (dust_remove)? 1 : 0;

        	// uses dynamic parallelization, top kernel is a single-thread one
        	int [] GridFullWarps =    {1, 1, 1};
        	int [] ThreadsFullWarps = {1, 1, 1};

        	Pointer kernelParameters = Pointer.to(
        			Pointer.to(gpu_tasks),                           // struct tp_task   * gpu_tasks,
        			Pointer.to(new int[] { num_task_tiles }),        // int                num_tiles,          // number of tiles in task list
        			// declare arrays in device code?
        			Pointer.to(gpu_texture_indices_ovlp),            // int              * gpu_texture_indices_ovlp,// packed tile + bits (now only (1 << 7)
        			Pointer.to(gpu_num_texture_ovlp),                // int              * num_texture_tiles,  // number of texture tiles to process (8 elements)
        			Pointer.to(gpu_woi),                             // int              * woi,                // x,y,width,height of the woi
        			// set smaller for LWIR - it is used to reduce work aread
        			Pointer.to(new int[] {img_width / DTT_SIZE}),    // int                width,  // <= TILESX, use for faster processing of LWIR images (should be actual + 1)
        			Pointer.to(new int[] {img_height / DTT_SIZE}),   // int                height); // <= TILESY, use for faster processing of LWIR images
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
        	cuLaunchKernel(GPU_RBGA_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize();
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
        		boolean   calc_extra)
        {
        	if (GPU_TEXTURES_kernel == null)
        	{
        		IJ.showMessage("Error", "No GPU kernel: GPU_TEXTURES_kernel");
        		return;
        	}
        	int tilesX =  img_width / DTT_SIZE;
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
        	int idust_remove =  (dust_remove)? 1 : 0;

        	int [] GridFullWarps =    {1, 1, 1};
        	int [] ThreadsFullWarps = {1, 1, 1};
        	
//        	CUdeviceptr gpu_diff_rgb_combo_local = calc_extra ? gpu_diff_rgb_combo : null;
        	Pointer kernelParameters = Pointer.to(
        			Pointer.to(gpu_tasks),                           // struct tp_task   * gpu_tasks,
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
        			calc_extra ? Pointer.to(gpu_diff_rgb_combo) : Pointer.to(new int[] { 0 }),
        			Pointer.to(new int[] { tilesX }));
        	cuCtxSynchronize();
        	// Call the kernel function
        	cuLaunchKernel(GPU_TEXTURES_kernel,
        			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],   // Grid dimension
        			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
        			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
        			kernelParameters, null);   // Kernel- and extra parameters
        	cuCtxSynchronize();
        }

        public int [] getCorrIndices() {
        	float [] fnum_corrs = new float[1];
        	cuMemcpyDtoH(Pointer.to(fnum_corrs), gpu_num_corr_tiles,  1 * Sizeof.FLOAT);
        	int num_corrs =      Float.floatToIntBits(fnum_corrs[0]);
        	float [] fcorr_indices = new float [num_corrs];
        	cuMemcpyDtoH(Pointer.to(fcorr_indices), gpu_corr_indices,  num_corrs * Sizeof.FLOAT);
        	int [] corr_indices = new int [num_corrs];
        	for (int i = 0; i < num_corrs; i++) {
        		corr_indices[i] = Float.floatToIntBits(fcorr_indices[i]);
        	}
        	num_corr_tiles = num_corrs;
        	return corr_indices;

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
        	int tilesX =  img_width / DTT_SIZE;
        	int tilesY =  img_height / DTT_SIZE;
        	float [][] extra = new float[num_tile_extra][tilesX*tilesY];
        	for (int i = 0; i < texture_indices.length; i++) {
        		if (((texture_indices[i] >> CORR_TEXTURE_BIT) & 1) != 0) {
        			int ntile = texture_indices[i] >>  CORR_NTILE_SHIFT;
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
        	woi.x =      Float.floatToIntBits(fwoi[0]) * DTT_SIZE;
        	woi.y =      Float.floatToIntBits(fwoi[1]) * DTT_SIZE;
        	woi.width =  Float.floatToIntBits(fwoi[2]) * DTT_SIZE;
        	woi.height = Float.floatToIntBits(fwoi[3]) * DTT_SIZE;
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
                copy_rbga.srcY =        4 + (woi.height +DTT_SIZE) * ncol;
                cuMemcpy2D(copy_rbga); // run copy
            }
            return rslt;
        }

        public float [] getFlatTextures(
        		int     num_tiles,
        		int     num_colors,
        		boolean keep_weights){
        	int texture_slices =     (num_colors + 1 + (keep_weights?(num_cams + num_colors + 1):0)); // number of texture slices
        	int texture_slice_size = (2 * DTT_SIZE)* (2 * DTT_SIZE);        // number of (float) elements in a single slice of a tile
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
        	int texture_slice_size = (2 * DTT_SIZE)* (2 * DTT_SIZE);
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
        	int texture_slice_size = (2 * DTT_SIZE)* (2 * DTT_SIZE);
        	double [][][][] textures = new double [woi.height][woi.width][num_slices][texture_slice_size];
        	for (int indx = 0; indx < indices.length; indx++) if ((indices[indx] & (1 << LIST_TEXTURE_BIT)) != 0){
        		int tile = indices[indx] >> CORR_NTILE_SHIFT;
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

        public double [][][][] doubleTextures(
        		Rectangle    woi,
        		int []       indices,
        		float []     ftextures,
        		int          full_width,
        		int          num_slices,
        		int          num_src_slices
        		){
        	int texture_slice_size = (2 * DTT_SIZE)* (2 * DTT_SIZE);
        	int texture_tile_size = texture_slice_size * num_src_slices ;
        	double [][][][] textures = new double [woi.height][woi.width][num_slices][texture_slice_size];
        	for (int indx = 0; indx < indices.length; indx++) if ((indices[indx] & (1 << LIST_TEXTURE_BIT)) != 0){
        		int tile = indices[indx] >> CORR_NTILE_SHIFT;
        		int tileX = tile % full_width;
        		int tileY = tile / full_width;
        		int wtileX = tileX - woi.x;
        		int wtileY = tileY - woi.y;
        		if ((wtileX < woi.width) && (wtileY < woi.height)) {
        			for (int slice = 0; slice < num_slices; slice++) {
        				for (int i = 0; i < texture_slice_size; i++) {
        					textures[wtileY][wtileX][slice][i] = ftextures[indx * texture_tile_size + slice * texture_slice_size + i];
        				}
        			}
        		}
        	}
        	return textures;
        }

        public float [][] getRBG (int ncam){
            int height = (img_height + DTT_SIZE);
            int width =  (img_width + DTT_SIZE);
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
            copyD2H.Height =          3 * height; // /2;

            cuMemcpy2D(copyD2H); // run copy

            float [][] fimg = new float [num_colors][ rslt_img_size];
            for (int ncol = 0; ncol < num_colors; ncol++) {
            	System.arraycopy(cpu_corr_image, ncol*rslt_img_size, fimg[ncol], 0, rslt_img_size);
            }
            return fimg;
        }
        
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
    		cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray,  module, "lpf_data");
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
    		cuModuleGetGlobal(constantMemoryPointer, constantMemorySizeArray,  module, const_name);
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
    		int dct_size = DTT_SIZE;
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
    		int dct_size = DTT_SIZE;
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
        
    } // end of public class GpuQuad

} // end of public class GPUTileProcessor
