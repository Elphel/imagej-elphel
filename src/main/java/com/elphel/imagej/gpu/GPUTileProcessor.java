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
import static jcuda.driver.JCudaDriver.cuDeviceGet;
import static jcuda.driver.JCudaDriver.cuDeviceGetAttribute;
import static jcuda.driver.JCudaDriver.cuInit;
import static jcuda.driver.JCudaDriver.cuLinkAddData;
import static jcuda.driver.JCudaDriver.cuLinkAddFile;
import static jcuda.driver.JCudaDriver.cuLinkComplete;
import static jcuda.driver.JCudaDriver.cuLinkCreate;
import static jcuda.driver.JCudaDriver.cuLinkDestroy;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleLoadDataEx;
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

import com.elphel.imagej.tileprocessor.Correlation2d;

import ij.IJ;
import ij.text.TextWindow;
import jcuda.Pointer;
import jcuda.driver.CUcontext;
import jcuda.driver.CUdevice;
import jcuda.driver.CUdeviceptr;
import jcuda.driver.CUfunction;
import jcuda.driver.CUlinkState;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;
import jcuda.driver.JITOptions;
import jcuda.nvrtc.JNvrtc;
import jcuda.nvrtc.nvrtcProgram;

public class GPUTileProcessor {
	public static boolean USE_DS_DP = false; // Use Dynamic Shared memory with Dynamic Parallelism (not implemented)  
	String LIBRARY_PATH = "/usr/local/cuda/targets/x86_64-linux/lib/libcudadevrt.a"; // linux
	// Can be downloaded and twice extracted from
	// https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-cudart-dev-11-2_11.2.152-1_amd64.deb
	// First deb itself, then data.tar.xz, and it will have usr/local/cuda/targets/x86_64-linux/lib/libcudadevrt.a inside
	// Found "cuda-cudart-dev" on https://ubuntu.pkgs.org/
	static String GPU_RESOURCE_DIR =              "kernels";
	static String [] GPU_KERNEL_FILES = {"dtt8x8.cuh","TileProcessor.cuh"};
	// "*" - generated defines, first index - separately compiled unit
	static String [][] GPU_SRC_FILES = {{"*","dtt8x8.h","dtt8x8.cu","geometry_correction.h","geometry_correction.cu","TileProcessor.h","TileProcessor.cuh"}};
	static String GPU_CONVERT_DIRECT_NAME =        "convert_direct";      // name in C code
	static String GPU_IMCLT_ALL_NAME =             "imclt_rbg_all";
	static String GPU_CORRELATE2D_NAME =           "correlate2D";         // name in C code
	static String GPU_CORRELATE2D_INTER_NAME =     "correlate2D_inter";   // name in C code
	static String GPU_CORR2D_COMBINE_NAME =        "corr2D_combine";      // name in C code
	static String GPU_CORR2D_NORMALIZE_NAME =      "corr2D_normalize";    // name in C code
	static String GPU_TEXTURES_NAME =              "textures_nonoverlap"; // name in C code
	static String GPU_RBGA_NAME =                  "generate_RBGA";       // name in C code
	static String GPU_ROT_DERIV =                  "calc_rot_deriv";      // calculate rotation matrices and derivatives
	static String GPU_SET_TILES_OFFSETS =          "get_tiles_offsets";   // calculate pixel offsets and disparity distortions
	static String GPU_CALCULATE_TILES_OFFSETS =    "calculate_tiles_offsets";   // calculate pixel offsets and disparity distortions
	static String GPU_CALC_REVERSE_DISTORTION =    "calcReverseDistortionTable"; // calculate reverse radial distortion table from gpu_geometry_correction
	// Kernels to use w/o Dynamic Parallelism    
	static String GPU_CLEAR_TEXTURE_LIST_NAME =    "clear_texture_list";
	static String GPU_MARK_TEXTURE_LIST_NAME =     "mark_texture_tiles";
	static String GPU_MARK_TEXTURE_NEIGHBOR_NAME = "mark_texture_neighbor_tiles";
	static String GPU_GEN_TEXTURE_LIST_NAME =      "gen_texture_list";
	static String GPU_CLEAR_TEXTURE_RBGA_NAME =    "clear_texture_rbga";
	static String GPU_TEXTURES_ACCUMULATE_NAME =   "textures_accumulate";
	static String GPU_CREATE_NONOVERLAP_LIST_NAME ="create_nonoverlap_list";
	static String GPU_ERASE_CLT_TILES_NAME =       "erase_clt_tiles";

	
	
//  pass some defines to gpu source code with #ifdef JCUDA
	public static int DTT_SIZE_LOG2 =             3;
	public static int DTT_SIZE =                  (1 << DTT_SIZE_LOG2);
	static int        THREADSX =                  DTT_SIZE;
	public static int MAX_NUM_CAMS =             16; // 4; Now - maximal number of sensors
	public static int CORR_VECTOR_MAX_LENGTH =   3 + 4 * MAX_NUM_CAMS; // 67; // 19; // TODO: update to fit for 16-sensor
//	public static int NUM_PAIRS =                 6; // top hor, bottom hor, left vert, right vert, main diagonal, other diagonal
//	public static int NUM_COLORS =                3;
//	public static int IMG_WIDTH =              2592;
//	public static int IMG_HEIGHT =             1936;
	static int        KERNELS_HOR =             164;
	static int        KERNELS_VERT =            123;
///	static int        KERNELS_LSTEP =             3; // 4;// FIXME: Make it dynamic: 3 for LWIR, 4 - for RGB?) 
	static int        THREADS_PER_TILE =          8;
	static int        TILES_PER_BLOCK =           4; // 8 - slower
	static int        CORR_THREADS_PER_TILE =     8;
	static int        CORR_TILES_PER_BLOCK	=     4;
	static int        CORR_TILES_PER_BLOCK_NORMALIZE = 4; // maybe change to 8?	
	static int        CORR_TILES_PER_BLOCK_COMBINE = 4; // increase to 16?
	static int        NUM_THREADS              = 32;

	static int        TEXTURE_THREADS_PER_TILE =  8; // 16;
	static int        TEXTURE_TILES_PER_BLOCK =   1;
	static int        IMCLT_THREADS_PER_TILE =   16;
	static int        IMCLT_TILES_PER_BLOCK =     4;
//	static int        TPTASK_SIZE =                 1+ 1+ NUM_CAMS * 2 + 1 + NUM_CAMS * 4 ; // tp_task structure size in floats
	static int        CLTEXTRA_SIZE =               8;
	static int        CORR_SIZE =                   (2* DTT_SIZE - 1) * (2* DTT_SIZE - 1); // 15x15
	public static int CORR_NTILE_SHIFT =          8;  // also for texture tiles list - not anymore 11/18/2022
	// FIXME: CORR_PAIRS_MASK will not work !!!	
	public static int CORR_PAIRS_MASK =        0x3f;  // lower bits used to address correlation pair for the selected tile
//	public static int CORR_TEXTURE_BIT =          7;  // bit 7 used to request texture for the tile GET RID !!!
	public static int TASK_CORR_BITS =            4;  // start of pair mask
	public static int TASK_TEXTURE_N_BIT =        0; // Texture with North neighbor
	public static int TASK_TEXTURE_E_BIT =        1; // Texture with East  neighbor
	public static int TASK_TEXTURE_S_BIT =        2; // Texture with South neighbor
	public static int TASK_TEXTURE_W_BIT =        3; // Texture with West  neighbor
//	public static int TASK_TEXTURE_BIT =          3;  // bit to request texture calculation int task field of struct tp_task
	
	public static int LIST_TEXTURE_BIT =          8; // 7;  // bit to request texture calculation
	public static int TEXT_NTILE_SHIFT =          9; // 8;  // split from CORR_NTILE_SHIFT
//	public static int CORR_OUT_RAD =              4;  // output radius of the correlations (implemented)
	public static double FAT_ZERO_WEIGHT =        0.0001; // add to port weights to avoid nan

	public static int THREADS_DYNAMIC_BITS =      5; // treads in block for CDP creation of the texture list

	public static int RBYRDIST_LEN =           5001; //for double, 10001 - float;   // length of rByRDist to allocate shared memory
	public static double RBYRDIST_STEP =          0.0004; //  for double, 0.0002 - for float; // to fit into GPU shared memory (was 0.001);
	public static int TILES_PER_BLOCK_GEOM =     32/MAX_NUM_CAMS; // blockDim.x = NUM_CAMS; blockDim.x = TILES_PER_BLOCK_GEOM
	public static int TASK_TEXTURE_BITS = ((1 << TASK_TEXTURE_N_BIT) | (1 << TASK_TEXTURE_E_BIT) | (1 << TASK_TEXTURE_S_BIT) | (1 << TASK_TEXTURE_W_BIT));


    int DTTTEST_BLOCK_WIDTH =        32; // may be read from the source code
    int DTTTEST_BLOCK_HEIGHT =       16; // may be read from the source code

    public boolean kernels_set = false;
    public boolean bayer_set =   false;

    CUfunction GPU_CONVERT_DIRECT_kernel =          null; // "convert_direct"
    CUfunction GPU_IMCLT_ALL_kernel =               null; // "imclt_rbg_all"
    CUfunction GPU_CORRELATE2D_kernel =             null; // "correlate2D"
    CUfunction GPU_CORRELATE2D_INTER_kernel =       null; // "correlate2D_inter"
    CUfunction GPU_CORR2D_COMBINE_kernel =          null; // "corr2D_combine"
    CUfunction GPU_CORR2D_NORMALIZE_kernel =        null; // "corr2D_normalize";
    CUfunction GPU_TEXTURES_kernel =                null; // "textures_nonoverlap"
    CUfunction GPU_RBGA_kernel =                    null; // "generate_RBGA"
    CUfunction GPU_ROT_DERIV_kernel =               null; // "calc_rot_deriv"
    CUfunction GPU_CALCULATE_TILES_OFFSETS_kernel = null; // "calculate_tiles_offsets"
    CUfunction GPU_CALC_REVERSE_DISTORTION_kernel = null; // "calcReverseDistortionTable"
// Kernels to use w/o Dynamic Parallelism    
    CUfunction GPU_CLEAR_TEXTURE_LIST_kernel =      null; // "clear_texture_list"
    CUfunction GPU_MARK_TEXTURE_LIST_kernel =       null; // "mark_texture_tiles"
    CUfunction GPU_MARK_TEXTURE_NEIGHBOR_kernel =   null; // "mark_texture_neighbor_tiles"
    CUfunction GPU_GEN_TEXTURE_LIST_kernel =        null; // "gen_texture_list"
    CUfunction GPU_CLEAR_TEXTURE_RBGA_kernel =      null; // "clear_texture_rbga"
    CUfunction GPU_TEXTURES_ACCUMULATE_kernel =     null; // "textures_accumulate"
    CUfunction GPU_CREATE_NONOVERLAP_LIST_kernel =  null; // "create_nonoverlap_list"
    CUfunction GPU_ERASE_CLT_TILES_kernel =         null; // "erase_clt_tiles"
    
    

    CUmodule    module; // to access constants memory
    //    private
    static long getPointerAddress(CUdeviceptr p)

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
        				"#define DTT_SIZE_LOG2 " +                  DTT_SIZE_LOG2+"\n"+
        				"#define THREADSX " +                       THREADSX+"\n"+
        				"#define NUM_CAMS " +                       MAX_NUM_CAMS+"\n"+
//        				"#define NUM_PAIRS " +                      NUM_PAIRS+"\n"+
//        				"#define NUM_COLORS " +                     NUM_COLORS+"\n"+
///        				"#define KERNELS_LSTEP " +                  KERNELS_LSTEP+"\n"+
        				"#define THREADS_PER_TILE " +               THREADS_PER_TILE+"\n"+
        				"#define TILES_PER_BLOCK " +                TILES_PER_BLOCK+"\n"+
        				"#define CORR_THREADS_PER_TILE " +          CORR_THREADS_PER_TILE+"\n"+
        				"#define CORR_TILES_PER_BLOCK " +           CORR_TILES_PER_BLOCK+"\n"+
        				"#define CORR_TILES_PER_BLOCK_NORMALIZE " + CORR_TILES_PER_BLOCK_NORMALIZE+"\n"+
        				"#define CORR_TILES_PER_BLOCK_COMBINE " +   CORR_TILES_PER_BLOCK_COMBINE+"\n"+
        				"#define NUM_THREADS " +                    NUM_THREADS+"\n"+
        				"#define TEXTURE_THREADS_PER_TILE " +       TEXTURE_THREADS_PER_TILE+"\n"+
        				"#define TEXTURE_TILES_PER_BLOCK " +        TEXTURE_TILES_PER_BLOCK+"\n"+
        				"#define IMCLT_THREADS_PER_TILE " +         IMCLT_THREADS_PER_TILE+"\n"+
        				"#define IMCLT_TILES_PER_BLOCK " +          IMCLT_TILES_PER_BLOCK+"\n"+
        				"#define CORR_NTILE_SHIFT " +               CORR_NTILE_SHIFT+"\n"+
//        				"#define CORR_PAIRS_MASK " +                CORR_PAIRS_MASK+"\n"+
//        				"#define CORR_TEXTURE_BIT " +               CORR_TEXTURE_BIT+"\n"+
        				"#define TASK_CORR_BITS " +                 TASK_CORR_BITS+"\n"+
        				"#define TASK_TEXTURE_N_BIT " +             TASK_TEXTURE_N_BIT+"\n"+
        				"#define TASK_TEXTURE_E_BIT " +             TASK_TEXTURE_E_BIT+"\n"+
        				"#define TASK_TEXTURE_S_BIT " +             TASK_TEXTURE_S_BIT+"\n"+
        				"#define TASK_TEXTURE_W_BIT " +             TASK_TEXTURE_W_BIT+"\n"+
        				"#define LIST_TEXTURE_BIT " +               LIST_TEXTURE_BIT+"\n"+
        				"#define TEXT_NTILE_SHIFT " +               TEXT_NTILE_SHIFT+"\n"+
//        				"#define CORR_OUT_RAD " +                   CORR_OUT_RAD+"\n" +
        				"#define FAT_ZERO_WEIGHT " +                FAT_ZERO_WEIGHT+"\n"+
        				"#define THREADS_DYNAMIC_BITS " +           THREADS_DYNAMIC_BITS+"\n"+
        				"#define RBYRDIST_LEN " +                   RBYRDIST_LEN+"\n"+
        				"#define RBYRDIST_STEP " +                  RBYRDIST_STEP+"\n"+
        				"#define TILES_PER_BLOCK_GEOM " +           TILES_PER_BLOCK_GEOM+"\n";
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
        boolean show_source = false; // true;
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
            if (show_source) {
//            	String lines[] = kernelSources[cunit].split("\\r?\\n");
//            	String body = "";
//            	for (int l = 0; l < lines.length; l++) {
//            		body += (l+1)+"\t"+lines[l]+"\n";
//            	}
//            	new TextWindow("GPU_Source_Code", "#\tline", body,400,800);
            	new TextWindow("GPU_Source_Code", "", kernelSources[cunit],400,800);
            }
        }

        // Create the kernel functions (first - just test)
        String [] func_names = {
        		GPU_CONVERT_DIRECT_NAME,
        		GPU_IMCLT_ALL_NAME,
        		GPU_CORRELATE2D_NAME,
        		GPU_CORRELATE2D_INTER_NAME,
        		GPU_CORR2D_COMBINE_NAME,
        		GPU_CORR2D_NORMALIZE_NAME,
        		GPU_TEXTURES_NAME,
        		GPU_RBGA_NAME,
        		GPU_ROT_DERIV,
        		GPU_CALCULATE_TILES_OFFSETS,
        		GPU_CALC_REVERSE_DISTORTION,
        		// Kernels to use w/o Dynamic Parallelism    
        		GPU_CLEAR_TEXTURE_LIST_NAME,
        		GPU_MARK_TEXTURE_LIST_NAME,
        		GPU_MARK_TEXTURE_NEIGHBOR_NAME,
        		GPU_GEN_TEXTURE_LIST_NAME,
        		GPU_CLEAR_TEXTURE_RBGA_NAME,
        		GPU_TEXTURES_ACCUMULATE_NAME,
        		GPU_CREATE_NONOVERLAP_LIST_NAME,
        		GPU_ERASE_CLT_TILES_NAME
        };
        CUfunction[] functions = createFunctions(kernelSources,
        		                                 func_names,
        		                                 capability); // on my - 75

        GPU_CONVERT_DIRECT_kernel =          functions[0];
        GPU_IMCLT_ALL_kernel =               functions[1];
        GPU_CORRELATE2D_kernel =             functions[2];
        GPU_CORRELATE2D_INTER_kernel =       functions[3];
        GPU_CORR2D_COMBINE_kernel =          functions[4];
        GPU_CORR2D_NORMALIZE_kernel =        functions[5];
        GPU_TEXTURES_kernel=                 functions[6];
        GPU_RBGA_kernel=                     functions[7];
        GPU_ROT_DERIV_kernel =               functions[8];
        GPU_CALCULATE_TILES_OFFSETS_kernel = functions[9];
        GPU_CALC_REVERSE_DISTORTION_kernel = functions[10];
     // Kernels to use w/o Dynamic Parallelism    
        GPU_CLEAR_TEXTURE_LIST_kernel =      functions[11];
        GPU_MARK_TEXTURE_LIST_kernel =       functions[12];
        GPU_MARK_TEXTURE_NEIGHBOR_kernel =   functions[13];
        GPU_GEN_TEXTURE_LIST_kernel =        functions[14];
        GPU_CLEAR_TEXTURE_RBGA_kernel =      functions[15];
        GPU_TEXTURES_ACCUMULATE_kernel =     functions[16];
        GPU_CREATE_NONOVERLAP_LIST_kernel =  functions[17];
        GPU_ERASE_CLT_TILES_kernel =         functions[18];
        

        System.out.println("GPU kernel functions initialized");
        System.out.println(GPU_CONVERT_DIRECT_kernel.toString());
        System.out.println(GPU_IMCLT_ALL_kernel.toString());
        System.out.println(GPU_CORRELATE2D_kernel.toString());
        System.out.println(GPU_CORRELATE2D_INTER_kernel.toString());
        System.out.println(GPU_CORR2D_COMBINE_kernel.toString());
        System.out.println(GPU_CORR2D_NORMALIZE_kernel.toString());
        System.out.println(GPU_TEXTURES_kernel.toString());
        System.out.println(GPU_RBGA_kernel.toString());
        System.out.println(GPU_ROT_DERIV_kernel.toString());
        System.out.println(GPU_CALCULATE_TILES_OFFSETS_kernel.toString());
        System.out.println(GPU_CALC_REVERSE_DISTORTION_kernel.toString());
        // Kernels to use w/o Dynamic Parallelism    
        System.out.println(GPU_CLEAR_TEXTURE_LIST_kernel.toString());
        System.out.println(GPU_MARK_TEXTURE_LIST_kernel.toString());
        System.out.println(GPU_MARK_TEXTURE_NEIGHBOR_kernel.toString());
        System.out.println(GPU_GEN_TEXTURE_LIST_kernel.toString());
        System.out.println(GPU_CLEAR_TEXTURE_RBGA_kernel.toString());
        System.out.println(GPU_TEXTURES_ACCUMULATE_kernel.toString());
        System.out.println(GPU_CREATE_NONOVERLAP_LIST_kernel.toString());
        System.out.println(GPU_ERASE_CLT_TILES_kernel.toString());
        // GPU data structures are now initialized through GpuQuad instances
    }
    

    public static String [] getCorrTitles() {
    	return new String []{"hor-top","hor-bottom","vert-left","vert-right","diag-main","diag-other","quad","cross"};
    }
    public static double [][] getCorr2DView(
    		int num_sensors,
    		int tilesX,
    		int tilesY,
    		int [] indices,
    		float [][] corr2d,
    		int [] wh){ // if is [2] - return width, height
    	int max_num_pairs = Correlation2d.getNumPairs(num_sensors);
    	if ((corr2d == null) || (corr2d.length == 0)) {
    		return new double [max_num_pairs][0];
    	}
    	int num_pairs = -1; // corr2d.length;
		for (int n = 0; n < indices.length; n++) {
			int np = indices[n] & CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
			if (np > num_pairs) num_pairs = np;
		}
		num_pairs++;
		if (num_pairs < 1) {
    		return new double [max_num_pairs][0];
		}
    	boolean [] bpairs = new boolean[num_pairs];
		for (int n = 0; n < indices.length; n++) {
			bpairs[indices[n] & CORR_PAIRS_MASK] = true;
		}    	
    	int first_pair = -1;
    	for (int i = 0; (i < bpairs.length) && (first_pair < 0); i++) {
    		if (bpairs[i]) first_pair = i; 
    	}
    	
    	int corr_size = (int)(Math.round(Math.sqrt(corr2d[0].length)));//  make smaller later?
    	
    	int width =  tilesX * (corr_size + 1) + 1;
    	int height = tilesY * (corr_size + 1) + 1;
    	double [][] data = new double [num_pairs][];
    	data[first_pair] = new double[height*width];
    	for (int ty = 0; ty < tilesY; ty++) {
    		for (int tx = 0; tx < tilesX; tx++) {
    			for (int i = 0; i< corr_size; i++) {
    				for (int j = 0; j < corr_size; j++) {
    					data[first_pair][(ty * (corr_size + 1) + i + 1) * width + (tx * (corr_size + 1) + j + 1)] = Double.NaN;
    				}
    			}
    		}
    	}
		for (int np = first_pair+1; np < num_pairs; np++) {
			if (bpairs[np]) {
				data[np] = data[first_pair].clone();
			}
		}
		for (int n = 0; n < indices.length; n++) {
			int nt = indices[n] >> CORR_NTILE_SHIFT;
			int np = indices[n] & CORR_PAIRS_MASK; // ((1 << CORR_NTILE_SHIFT) - 1); // np should
			assert np < num_pairs : "invalid correllation pair";
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
    		System.out.println("Looking for GPU kernel ["+i+"]: "+kernelNames[i]);
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

} // end of public class GPUTileProcessor
