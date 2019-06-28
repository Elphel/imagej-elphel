package com.elphel.imagej.common;
/**
 * The code below is extracted form ImageJ plugin GaussianBlur.java, stripped of ImageProcessor and used (double) instead of (float) arrays.
 * The following are notes from the original file:
 *
 *
 *  This plug-in filter uses convolution with a Gaussian function for smoothing.
 * 'Radius' means the radius of decay to exp(-0.5) ~ 61%, i.e. the standard
 * deviation sigma of the Gaussian (this is the same as in Photoshop, but
 * different from the previous ImageJ function 'Gaussian Blur', where a value
 * 2.5 times as much has to be entered.
 * - Like all convolution operations in ImageJ, it assumes that out-of-image
 * pixels have a value equal to the nearest edge pixel. This gives higher
 * weight to edge pixels than pixels inside the image, and higher weight
 * to corner pixels than non-corner pixels at the edge. Thus, when smoothing
 * with very high blur radius, the output will be dominated by the edge
 * pixels and especially the corner pixels (in the extreme case, with
 * a blur radius of e.g. 1e20, the image will be raplaced by the average
 * of the four corner pixels).
 * - For increased speed, except for small blur radii, the lines (rows or
 * columns of the image) are downscaled before convolution and upscaled
 * to their original length thereafter.
 *
 * Version 03-Jun-2007 M. Schmid with preview, progressBar stack-aware,
 * snapshot via snapshot flag; restricted range for resetOutOfRoi
 *
 */


public class DoubleGaussianBlur {
	   /* the standard deviation of the Gaussian*/
//    private static double sigma = 2.0;
    /* whether sigma is given in units corresponding to the pixel scale (not pixels)*/
//    private static boolean sigmaScaled = false;
    /* The flags specifying the capabilities and needs */
//    private int flags = DOES_ALL|SUPPORTS_MASKING|PARALLELIZE_STACKS|KEEP_PREVIEW;
//    private ImagePlus imp;              // The ImagePlus of the setup call, needed to get the spatial calibration
//    private boolean hasScale = false;   // whether the image has an x&y scale
    private int nPasses = 1;            // The number of passes (filter directions * color channels * stack slices)
//	private int nChannels = 1;        // The number of color channels
    private int pass;                   // Current pass

    /* Default constructor */
    public DoubleGaussianBlur() {
    }


	  public void blurDouble(double[] pixels,
			                 int width,
			                 int height,
			                 double sigmaX,
			                 double sigmaY,
			                 double accuracy) {
	        if (sigmaX > 0)
	            blur1Direction(pixels, width,height, sigmaX, accuracy, true);
	        if (Thread.currentThread().isInterrupted()) return; // interruption for new parameters during preview?
	        if (sigmaY > 0)
	            blur1Direction(pixels, width,height, sigmaY, accuracy, false);
	        return;
	    }

	    /* Blur an image in one direction (x or y) by a Gaussian.
	     * @param ip        The Image with the original data where also the result will be stored
	     * @param sigma     Standard deviation of the Gaussian
	     * @param accuracy  Accuracy of kernel, should not be > 0.02
	     * @param xDirection True for blurring in x direction, false for y direction
	     * @param extraLines Number of lines (parallel to the blurring direction)
	     *                  below and above the roi bounds that should be processed.
	     */
	    public void blur1Direction(double [] pixels,
	    		                   int        width,
	    		                   int       height,
	    		                   double     sigma,
	    		                   double   accuracy,
	                               boolean xDirection
//	                               int extraLines
	                               ) {
	        final int UPSCALE_K_RADIUS = 2;             //number of pixels to add for upscaling
	        final double MIN_DOWNSCALED_SIGMA = 4.;     //minimum standard deviation in the downscaled image
//	        float[] pixels = (float[])ip.getPixels();
//	        int width = ip.getWidth();
//	        int height = ip.getHeight();
//	        Rectangle roi = ip.getRoi();
	        int length = xDirection ? width : height;   //number of points per line (line can be a row or column)
	        int pointInc = xDirection ? 1 : width;      //increment of the pixels array index to the next point in a line
	        int lineInc = xDirection ? width : 1;       //increment of the pixels array index to the next line
//	        int lineFrom = (xDirection ? roi.y : roi.x) - extraLines;  //the first line to process
//	        if (lineFrom < 0) lineFrom = 0;
	        int lineFrom = 0;  //the first line to process
//	        int lineTo = (xDirection ? roi.y+roi.height : roi.x+roi.width) + extraLines; //the last line+1 to process
//	        if (lineTo > (xDirection ? height:width)) lineTo = (xDirection ? height:width);
	        int lineTo = (xDirection ? height:width);
//	        int writeFrom = xDirection? roi.x : roi.y;  //first point of a line that needs to be written
//	        int writeTo = xDirection ? roi.x+roi.width : roi.y+roi.height;
	        int writeFrom = 0;  //first point of a line that needs to be written
	        int writeTo = xDirection ? width : height;
//	        int inc = Math.max((lineTo-lineFrom)/(100/(nPasses>0?nPasses:1)+1),20);
	        pass++;
	        if (pass>nPasses) pass =1;
//	        Thread thread = Thread.currentThread();     // needed to check for interrupted state
	        if (sigma > 2*MIN_DOWNSCALED_SIGMA + 0.5) {
	             /* large radius (sigma): scale down, then convolve, then scale up */
	            int reduceBy = (int)Math.floor(sigma/MIN_DOWNSCALED_SIGMA); //downscale by this factor
	            if (reduceBy > length) reduceBy = length;
	            /* Downscale gives std deviation sigma = 1/sqrt(3); upscale gives sigma = 1/2. (in downscaled pixels) */
	            /* All sigma^2 values add to full sigma^2  */
	            double sigmaGauss = Math.sqrt(sigma*sigma/(reduceBy*reduceBy) - 1./3. - 1./4.);
	            int maxLength = (length+reduceBy-1)/reduceBy + 2*(UPSCALE_K_RADIUS + 1); //downscaled line can't be longer
	            double[][] gaussKernel = makeGaussianKernel(sigmaGauss, accuracy, maxLength);
	            int kRadius = gaussKernel[0].length*reduceBy; //Gaussian kernel radius after upscaling
	            int readFrom = (writeFrom-kRadius < 0) ? 0 : writeFrom-kRadius; //not including broadening by downscale&upscale
	            int readTo = (writeTo+kRadius > length) ? length : writeTo+kRadius;
	            int newLength = (readTo-readFrom+reduceBy-1)/reduceBy + 2*(UPSCALE_K_RADIUS + 1);
	            int unscaled0 = readFrom - (UPSCALE_K_RADIUS + 1)*reduceBy; //input point corresponding to cache index 0
	            //IJ.log("reduce="+reduceBy+", newLength="+newLength+", unscaled0="+unscaled0+", sigmaG="+(float)sigmaGauss+", kRadius="+gaussKernel[0].length);
	            double[] downscaleKernel = makeDownscaleKernel(reduceBy);
	            double[] upscaleKernel = makeUpscaleKernel(reduceBy);
	            double[] cache1 = new double[newLength];  //holds data after downscaling
	            double[] cache2 = new double[newLength];  //holds data after convolution
	            int pixel0 = lineFrom*lineInc;
	            for (int line=lineFrom; line<lineTo; line++, pixel0+=lineInc) {
	                downscaleLine(pixels, cache1, downscaleKernel, reduceBy, pixel0, unscaled0, length, pointInc, newLength);
	                convolveLine(cache1, cache2, gaussKernel, 0, newLength, 1, newLength-1, 0, 1);
	                upscaleLine(cache2, pixels, upscaleKernel, reduceBy, pixel0, unscaled0, writeFrom, writeTo, pointInc);
	            }
	        } else {
	            /* small radius: normal convolution */
	            double[][] gaussKernel = makeGaussianKernel(sigma, accuracy, length);
	            int kRadius = gaussKernel[0].length;
	            double[] cache = new double[length];          //input for convolution, hopefully in CPU cache
	            int readFrom = (writeFrom-kRadius < 0) ? 0 : writeFrom-kRadius;
	            int readTo = (writeTo+kRadius > length) ? length : writeTo+kRadius;
	            int pixel0 = lineFrom*lineInc;
	            for (int line=lineFrom; line<lineTo; line++, pixel0+=lineInc) {
	                int p = pixel0 + readFrom*pointInc;
	                for (int i=readFrom; i<readTo; i++ ,p+=pointInc)
	                    cache[i] = pixels[p];
	                convolveLine(cache, pixels, gaussKernel, readFrom, readTo, writeFrom, writeTo, pixel0, pointInc);
	            }
	        }
	        return;
	    }
	    /* Create a 1-dimensional normalized Gaussian kernel with standard deviation sigma
	     *  and the running sum over the kernel
	     *  Note: this is one side of the kernel only, not the full kernel as used by the
	     *  Convolver class of ImageJ.
	     *  To avoid a step due to the cutoff at a finite value, the near-edge values are
	     *  replaced by a 2nd-order polynomial with its minimum=0 at the first out-of-kernel
	     *  pixel. Thus, the kernel function has a smooth 1st derivative in spite of finite
	     *  length.
	     *
	     * @param sigma     Standard deviation, i.e. radius of decay to 1/sqrt(e), in pixels.
	     * @param accuracy  Relative accuracy; for best results below 0.01 when processing
	     *                  8-bit images. For short or float images, values of 1e-3 to 1e-4
	     *                  are better (but increase the kernel size and thereby the
	     *                  processing time). Edge smoothing will fail with very poor
	     *                  accuracy (above approx. 0.02)
	     * @param maxRadius Maximum radius of the kernel: Limits kernel size in case of
	     *                  large sigma, should be set to image width or height. For small
	     *                  values of maxRadius, the kernel returned may have a larger
	     *                  radius, however.
	     * @return          A 2*n array. Array[0][n] is the kernel, decaying towards zero,
	     *                  which would be reached at kernel.length (unless kernel size is
	     *                  limited by maxRadius). Array[1][n] holds the sum over all kernel
	     *                  values > n, including non-calculated values in case the kernel
	     *                  size is limited by <code>maxRadius</code>.
	     */
	    public double[][] makeGaussianKernel(double sigma, double accuracy, int maxRadius) {
	        int kRadius = (int)Math.ceil(sigma*Math.sqrt(-2*Math.log(accuracy)))+1;
	        if (maxRadius < 50) maxRadius = 50;         // too small maxRadius would result in inaccurate sum.
	        if (kRadius > maxRadius) kRadius = maxRadius;
	        double[][] kernel = new double[2][kRadius];
	        for (int i=0; i<kRadius; i++)               // Gaussian function
	            kernel[0][i] = (Math.exp(-0.5*i*i/sigma/sigma));
	        if (kRadius < maxRadius && kRadius > 3) {   // edge correction
	            double sqrtSlope = Double.MAX_VALUE;
	            int r = kRadius;
	            while (r > kRadius/2) {
	                r--;
	                double a = Math.sqrt(kernel[0][r])/(kRadius-r);
	                if (a < sqrtSlope)
	                    sqrtSlope = a;
	                else
	                    break;
	            }
	            for (int r1 = r+2; r1 < kRadius; r1++)
	                kernel[0][r1] = (kRadius-r1)*(kRadius-r1)*sqrtSlope*sqrtSlope;
	        }
	        double sum;                                 // sum over all kernel elements for normalization
	        if (kRadius < maxRadius) {
	            sum = kernel[0][0];
	            for (int i=1; i<kRadius; i++)
	                sum += 2*kernel[0][i];
	        } else
	            sum = sigma * Math.sqrt(2*Math.PI);

	        double rsum = 0.5 + 0.5*kernel[0][0]/sum;
	        for (int i=0; i<kRadius; i++) {
	            double v = (kernel[0][i]/sum);
	            kernel[0][i] = v;
	            rsum -= v;
	            kernel[1][i] = rsum;
	            //IJ.log("k["+i+"]="+(float)v+" sum="+(float)rsum);
	        }
	        return kernel;
	    }
	    /* Scale a line (row or column of a FloatProcessor or part thereof)
	     * down by a factor <code>reduceBy</code> and write the result into
	     * <code>cache</code>.
	     * Input line pixel # <code>unscaled0</code> will correspond to output
	     * line pixel # 0. <code>unscaled0</code> may be negative. Out-of-line
	     * pixels of the input are replaced by the edge pixels.
	     */
	    void downscaleLine(double[] pixels, double[] cache, double[] kernel,
	            int reduceBy, int pixel0, int unscaled0, int length, int pointInc, int newLength) {
	    	if (pixel0 > pixels.length) {
	    		System.out.println("++++++ Error in DoubleGaussianBlur, pixel0="+pixel0+", pixels.length="+(pixels.length));
	    		return;
	    	}
	        double first = pixels[pixel0];
	        double last = pixels[pixel0 + pointInc*(length-1)];
	        int xin = unscaled0 - reduceBy/2;
	        int p = pixel0 + pointInc*xin;
	        for (int xout=0; xout<newLength; xout++) {
	            double v = 0;
	            for (int x=0; x<reduceBy; x++, xin++, p+=pointInc) {
	                v += kernel[x] * ((xin-reduceBy < 0) ? first : ((xin-reduceBy >= length) ? last : pixels[p-pointInc*reduceBy]));
	                v += kernel[x+reduceBy] * ((xin < 0) ? first : ((xin >= length) ? last : pixels[p]));
	                v += kernel[x+2*reduceBy] * ((xin+reduceBy < 0) ? first : ((xin+reduceBy >= length) ? last : pixels[p+pointInc*reduceBy]));
	            }
	            cache[xout] = v;
	        }
	    }

	    /* Create a kernel for downscaling. The kernel function preserves
	     * norm and 1st moment (i.e., position) and has fixed 2nd moment,
	     * (in contrast to linear interpolation).
	     * In scaled space, the length of the kernel runs from -1.5 to +1.5,
	     * and the standard deviation is 1/2.
	     * Array index corresponding to the kernel center is
	     * unitLength*3/2
	     */
	    double[] makeDownscaleKernel (int unitLength) {
	        int mid = unitLength*3/2;
	        double[] kernel = new double[3*unitLength];
	        for (int i=0; i<=unitLength/2; i++) {
	            double x = i/(double)unitLength;
	            double v = (0.75-x*x)/unitLength;
	            kernel[mid-i] = v;
	            kernel[mid+i] = v;
	        }
	        for (int i=unitLength/2+1; i<(unitLength*3+1)/2; i++) {
	            double x = i/(double)unitLength;
	            double v = (0.125 + 0.5*(x-1)*(x-2))/unitLength;
	            kernel[mid-i] = v;
	            kernel[mid+i] = v;
	        }
	        return kernel;
	    }

	    /* Scale a line up by factor <code>reduceBy</code> and write as a row
	     * or column (or part thereof) to the pixels array of a FloatProcessor.
	     */
	    void upscaleLine (double[] cache, double[] pixels, double[] kernel,
	            int reduceBy, int pixel0, int unscaled0, int writeFrom, int writeTo, int pointInc) {
	        int p = pixel0 + pointInc*writeFrom;
	        for (int xout = writeFrom; xout < writeTo; xout++, p+=pointInc) {
	            int xin = (xout-unscaled0+reduceBy-1)/reduceBy; //the corresponding point in the cache (if exact) or the one above
	            int x = reduceBy - 1 - (xout-unscaled0+reduceBy-1)%reduceBy;
	            pixels[p] = cache[xin-2]*kernel[x]
	                    + cache[xin-1]*kernel[x+reduceBy]
	                    + cache[xin]*kernel[x+2*reduceBy]
	                    + cache[xin+1]*kernel[x+3*reduceBy];
	        }
	    }

	    /* Create a kernel for upscaling. The kernel function is a convolution
	     *  of four unit squares, i.e., four uniform kernels with value +1
	     *  from -0.5 to +0.5 (in downscaled coordinates). The second derivative
	     *  of this kernel is smooth, the third is not. Its standard deviation
	     *  is 1/sqrt(3) in downscaled cordinates.
	     *  The kernel runs from [-2 to +2[, corresponding to array index
	     *  0 ... 4*unitLength (whereby the last point is not in the array any more).
	     */
	    double[] makeUpscaleKernel (int unitLength) {
	        double[] kernel = new double[4*unitLength];
	        int mid = 2*unitLength;
	        kernel[0] = 0;
	        for (int i=0; i<unitLength; i++) {
	            double x = i/(double)unitLength;
	            double v = ((2./3. -x*x*(1-0.5*x)));
	            kernel[mid+i] = v;
	            kernel[mid-i] = v;
	        }
	        for (int i=unitLength; i<2*unitLength; i++) {
	            double x = i/(double)unitLength;
	            double v = (2.-x)*(2.-x)*(2.-x)/6.;
	            kernel[mid+i] = v;
	            kernel[mid-i] = v;
	        }
	        return kernel;
	    }

	    /* Convolve a line with a symmetric kernel and write to a separate array,
	     * possibly the pixels array of a FloatProcessor (as a row or column or part thereof)
	     *
	     * @param input     Input array containing the line
	     * @param pixels    Float array for output, can be the pixels of a FloatProcessor
	     * @param kernel    "One-sided" kernel array, kernel[0][n] must contain the kernel
	     *                  itself, kernel[1][n] must contain the running sum over all
	     *                  kernel elements from kernel[0][n+1] to the periphery.
	     *                  The kernel must be normalized, i.e. sum(kernel[0][n]) = 1
	     *                  where n runs from the kernel periphery (last element) to 0 and
	     *                  back. Normalization should include all kernel points, also these
	     *                  not calculated because they are not needed.
	     * @param readFrom  First array element of the line that must be read.
	     *                  <code>writeFrom-kernel.length</code> or 0.
	     * @param readTo    Last array element+1 of the line that must be read.
	     *                  <code>writeTo+kernel.length</code> or <code>input.length</code>
	     * @param writeFrom Index of the first point in the line that should be written
	     * @param writeTo   Index+1 of the last point in the line that should be written
	     * @param point0    Array index of first element of the 'line' in pixels (i.e., lineNumber * lineInc)
	     * @param pointInc  Increment of the pixels array index to the next point (for an ImageProcessor,
	     *                  it should be <code>1</code> for a row, <code>width</code> for a column)
	     */
	    public void convolveLine(double[] input, double[] pixels, double[][] kernel, int readFrom,
	            int readTo, int writeFrom, int writeTo, int point0, int pointInc) {
	        int length = input.length;
	        double first = input[0];                 //out-of-edge pixels are replaced by nearest edge pixels
	        double last = input[length-1];
	        double[] kern = kernel[0];               //the kernel itself
	        double kern0 = kern[0];
	        double[] kernSum = kernel[1];            //the running sum over the kernel
	        int kRadius = kern.length;
	        int firstPart = kRadius < length ? kRadius : length;
	        int p = point0 + writeFrom*pointInc;
	        int i = writeFrom;
	        for (; i<firstPart; i++,p+=pointInc) {  //while the sum would include pixels < 0
	            double result = input[i]*kern0;
	            result += kernSum[i]*first;
	            if (i+kRadius>length) result += kernSum[length-i-1]*last;
	            for (int k=1; k<kRadius; k++) {
	                double v = 0;
	                if (i-k >= 0) v += input[i-k];
	                if (i+k<length) v+= input[i+k];
	                result += kern[k] * v;
	            }
	            pixels[p] = result;
	        }
	        int iEndInside = length-kRadius<writeTo ? length-kRadius : writeTo;
	        for (;i<iEndInside;i++,p+=pointInc) {   //while only pixels within the line are be addressed (the easy case)
	            double result = input[i]*kern0;
	            for (int k=1; k<kRadius; k++)
	                result += kern[k] * (input[i-k] + input[i+k]);
	            pixels[p] = result;
	        }
	        for (; i<writeTo; i++,p+=pointInc) {    //while the sum would include pixels >= length
	            double result = input[i]*kern0;
	            if (i<kRadius) result += kernSum[i]*first;
	            if (i+kRadius>=length) result += kernSum[length-i-1]*last;
	            for (int k=1; k<kRadius; k++) {
	                double v = 0;
	                if (i-k >= 0) v += input[i-k];
	                if (i+k<length) v+= input[i+k];
	                result += kern[k] * v;
	            }
	            pixels[p] = result;
	        }
	    }



}
