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
import static jcuda.driver.JCudaDriver.cuMemAllocPitch;
import static jcuda.driver.JCudaDriver.cuMemFree;
import static jcuda.driver.JCudaDriver.cuMemcpy2D;
import static jcuda.driver.JCudaDriver.cuModuleGetFunction;
import static jcuda.driver.JCudaDriver.cuModuleLoadData;
import static jcuda.nvrtc.JNvrtc.nvrtcCompileProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcCreateProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcDestroyProgram;
import static jcuda.nvrtc.JNvrtc.nvrtcGetPTX;
import static jcuda.nvrtc.JNvrtc.nvrtcGetProgramLog;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

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
	static String  GPU_DTT24_NAME = "GPU_DTT24_DRV";       // this.kernelFunction = createFunction(sourceCode, "GPU_DTT24_DRV"); // "invert");
    int DTTTEST_BLOCK_WIDTH =        32; // may be read from the source code
    int DTTTEST_BLOCK_HEIGHT =       16; // may be read from the source code
    int DTT_SIZE =                    8; // may be read from the source code

    private CUfunction GPU_DTT24_kernel = null;

    public GPUTileProcessor() throws IOException
    {
    	int su = setup();
    	if (su < 0) {
    		new IllegalArgumentException ("setup() returned "+su);
    	}
    }

    // run kernel with dttx
    public void exec_dtt24(float src_pixels[],float dst_pixels[], int width, int dtt_mode)
    {
        if (GPU_DTT24_kernel == null)
        {
            IJ.showMessage("Error", "No GPU kernel");
            return;
        }
        int height = src_pixels.length / width;
        long width_in_bytes = width * Sizeof.FLOAT;
        CUdeviceptr src_dpointer = new CUdeviceptr();
        CUdeviceptr dst_dpointer = new CUdeviceptr();
        long [] device_stride = new long [1];
        cuMemAllocPitch (
        		src_dpointer,   // CUdeviceptr dptr,
        		device_stride,  // long[] pPitch,
        		width_in_bytes, // long WidthInBytes,
        		height,         // long Height,
                Sizeof.FLOAT);  // int ElementSizeBytes)
        int pitchInElements = (int)(device_stride[0] / Sizeof.FLOAT);
        cuMemAllocPitch (
        		dst_dpointer,   // CUdeviceptr dptr,
        		device_stride,  // long[] pPitch,
        		width_in_bytes, // long WidthInBytes,
        		height,         // long Height,
                Sizeof.FLOAT);  // int ElementSizeBytes)
        CUDA_MEMCPY2D copyH2D =   new CUDA_MEMCPY2D();
        copyH2D.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
        copyH2D.srcHost =         Pointer.to(src_pixels);
        copyH2D.srcPitch =        width_in_bytes;

        copyH2D.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
        copyH2D.dstDevice =       src_dpointer;
        copyH2D.dstPitch =        device_stride[0];

        copyH2D.WidthInBytes =    width_in_bytes;
        copyH2D.Height =          height; // /4;

// for copying results back to host
        CUDA_MEMCPY2D copyD2H =   new CUDA_MEMCPY2D();
        copyD2H.srcMemoryType =   CUmemorytype.CU_MEMORYTYPE_DEVICE;
        copyD2H.srcDevice =       dst_dpointer; // ((test & 1) ==0) ? src_dpointer : dst_dpointer; // copy same data
        copyD2H.srcPitch =        device_stride[0];

        copyD2H.dstMemoryType =   CUmemorytype.CU_MEMORYTYPE_HOST;
        copyD2H.dstHost =         Pointer.to(dst_pixels);
        copyD2H.dstPitch =        width_in_bytes;

        copyD2H.WidthInBytes =    width_in_bytes;
        copyD2H.Height =          height; // /2;

        // Set up the kernel parameters: A pointer to an array
        // of pointers which point to the actual values.
        Pointer kernelParameters = Pointer.to(
            Pointer.to(dst_dpointer),
            Pointer.to(src_dpointer),
            Pointer.to(new int[] { pitchInElements }),
            Pointer.to(new int[] { dtt_mode })
        );

        int [] GridFullWarps =    {width / DTTTEST_BLOCK_WIDTH, height / DTTTEST_BLOCK_HEIGHT, 1};
        int [] ThreadsFullWarps = {DTT_SIZE, DTTTEST_BLOCK_WIDTH/DTT_SIZE, DTTTEST_BLOCK_HEIGHT/DTT_SIZE};

        // Actual work starts here:
        cuMemcpy2D(copyH2D);
        cuCtxSynchronize();
        	// Call the kernel function
    	cuLaunchKernel(GPU_DTT24_kernel,
    			GridFullWarps[0],    GridFullWarps[1],   GridFullWarps[2],  // Grid dimension
    			ThreadsFullWarps[0], ThreadsFullWarps[1],ThreadsFullWarps[2],// Block dimension
    			0, null,                 // Shared memory size and stream (shared - only dynamic, static is in code)
    			kernelParameters, null);   // Kernel- and extra parameters

        // Copy the data from the device back to the host
        cuMemcpy2D(copyD2H);
        // clean up
        cuMemFree(src_dpointer);
        cuMemFree(dst_dpointer);
    }

    public int setup() throws IOException // String arg, ImagePlus imagePlus)
    {

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
        String cuFileName = "/home/eyesis/workspace-python3/nvidia_dct8x8/src/dtt8x8.cuh";// "dtt8x8.cuh";
        String sourceCode = readFileAsString(cuFileName); // readResourceAsString(cuFileName);
        if (sourceCode == null)
        {
            IJ.showMessage("Error",
                "Could not read the kernel source code");
            return -1;
        }

        // Create the kernel function
        this.GPU_DTT24_kernel = createFunction(sourceCode, GPU_DTT24_NAME); // "invert");

    	return 0;
    }

    /**
     * Create the CUDA function object for the kernel function with the
     * given name that is contained in the given source code
     *
     * @param sourceCode The source code
     * @param kernelName The kernel function name
     * @return
     * @throws IOException
     */
    private static CUfunction createFunction(
        String sourceCode, String kernelName) throws IOException
    {
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
        CUmodule module = new CUmodule();
        cuModuleLoadData(module, ptx[0]);

        // Find the function in the source by name, get its pointer
        CUfunction function = new CUfunction();
        cuModuleGetFunction(function, module, kernelName);

        return function;
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

}
