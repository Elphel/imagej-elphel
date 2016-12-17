import ij.IJ;
import ij.ImagePlus;

import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.SwingUtilities;


public class DenseCorrespondence {
	public double [][] disparityScales=       null; // for each image - a pair of {scaleX, scaleY} or null if undefined (interSensor has the same)
	public ImagePlus impDisparity=null;
	public int corrFFTSize; // to properties    - 32 *4
	public int overlapStep; // to properties    - 16
	public int paddedSize;  // same as in zTile - 20 *4
	public int subpixel;    // same as in cyclopeanTile - subdivide mask visibility pixels
	private double disparityPerEntry; //disparity increment (in pix) per array element
	private int tilesX;
	private int tilesY;
	public String title;
//	private float [][] pixels=null; 
	private float [] centerPixels=null; // disparity arrays combined for the center virtual image
//	private float [][] syntheticPixels=null; // disparity arrays for individual images restored from the centerPixels
	private BitSet innerMask=null; // will be provided to zMap instances to quickly find that there are no inner (sans padding) pixels
	private int  []   borderMask=null; // will be provided to zMap +1:top+2:bottom+8:left+16:right (to prevent roll over when iterating around
//	private double centerPixelsFatZero=0.0; 
	//	private int [][] imagePairs=null;
	private int [][] imagePairIndices=null;
//	private double [] doubleTileWindow=null;
	private CyclopeanTile [][] cyclopeanMap=null;
	private Rectangle zMapWOI=null; // full: 0,0,tilesX,tilesY
	public  Photometric photometric=null;
	public int getPadding() {return (this.paddedSize-this.overlapStep)/2;}



	public void disparitySweepTile ( // never used so far
			int tileX,
			int tileY,
			CyclopeanTile [][] allCyclopeanMap,
			//				int nImg,
			//				int [] sImgSet, // list of second images
			int [][]imgPairs,
			double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
			int imageFullWidth,
			double blurVarianceSigma,
			//				int refineTilePeriod,
			double subTilePhaseCoeff,
			double subTileHighPassSigma,
			double subTileLowPassSigma,
			//				double refineCorrMaxDistance,
			//				double refineCorrThreshold,
			int refineSubPixel,
			double zMapMinForeground,
			int zMapVarMask,
			double [] zMapVarThresholds,
			double [] zMapVarWeights,
			int       auxVarMode,
			int     normalize,
			int zMapCorrMask,
			double [] zMapCorrThresholds,
			double [] zMapCorrThresholdsRel,
			double [] zMapCorrWeights,
			double [] window,
			DoubleFHT doubleFHT,
			double [] subTileWindow,
			DoubleFHT subTileFHT,
			// new arguments
			int combineMode, // different image pairs - 0 
			double disparityMax,
			double disparityMin,
			double minAbsolute, // or NaN - will use enabled/disabled state of the tile
			double minRelative,
			boolean filterByForeground, // apply known certain masks
			double filterByForegroundMargin,
			boolean filterByDisabled,
			double disparityTolearnce,
			double maskBlurSigma,
			double corrHighPassSigma, // subtract blurred version to minimize correlation caused by masks
			//tone-matching statistical parameters
			double varianceBlurScale,
			double kLocal,
			boolean refineTileDisparity,
			int matchStatMode,
			int threadsMax,
			boolean showProgress,
			int debugLevel){
		// TODO: also calculate "unlikely" - high autocorrelation, not occluded, low inter-correlation
		double [] normVarWeights= normalizeWeights(zMapVarMask,zMapVarWeights);
		double [] normCorrWeights=normalizeWeights(zMapCorrMask,zMapCorrWeights);
		//          int auxChannelNumber=3;
		CyclopeanTile cyclopeanTile=allCyclopeanMap[tileY][tileX];
		int tileOverlap=cyclopeanTile.getOverlap();
		int paddedSize=cyclopeanTile.getPaddedSize();
		int paddedLength=cyclopeanTile.getPaddedLength();
		int size=this.overlapStep*2;
		int length=size*size;
		double[] zeros=new double [length];
		for (int i=0;i<length;i++) zeros[i]=0.0;
		int margin=size/4-tileOverlap;
		//			double [][][][] slices=new double [sImgSet.length][][][];
		double [][][][] slices=new double [imgPairs.length][][][];
		int [] dirs1={1,size+1,size,size-1,-1,-size-1,-size,-size+1,0};
		//			double [][][][] variance=new double [sImgSet.length][][][]; // [sIndex][pair 0/1][chn][pixel
		//			int length=4*this.overlapStep*this.overlapStep; // initialize it to a correct value right here?
		int [] borderMask=cyclopeanTile.getBorderMask();
		// Iterate through pairsd to this image - oter pairs will be processed separately and the result (likelyhood of belonguing to
		// the particular plane can be evaluated  (they use the same "radar" data
		cyclopeanTile.setMinCorrelations(minAbsolute, minRelative, true); // Set here or before? include FG
		if (debugLevel>3) {
			System.out.println ("zTile.setMinCorrelations("+minAbsolute+","+ minRelative+")");
			boolean [] enabled=cyclopeanTile.enabledPlane;
			for (int plane=0;plane<cyclopeanTile.getNumberOfPlanes();plane++) {
				System.out.println(plane+": disparity="+cyclopeanTile.getPlaneDisparity(plane)+" strength="+cyclopeanTile.getPlaneStrength(plane)+" enabled="+
						enabled[plane]);

			}
			System.out.println("cyclopeanTile.minAbsolute="+cyclopeanTile.minAbsolute+" cyclopeanTile.minRelative="+cyclopeanTile.minRelative);
		}
		if (Double.isNaN(window[0])) {
			int index=0;
			int quarterSize=size/4; //8
			int halfSize=size/2; // 16
			int size34=3*size/4; // 24
			double[] window1d=doubleFHT.getHamming1d(halfSize); // Hamming
			for (int iy=0;iy<size;iy++) {
				double wy=(iy<quarterSize)?window1d[iy]:((iy>size34)?window1d[iy-halfSize]:1.0);
				for (int ix=0;ix<size;ix++) {
					double wx=(ix<quarterSize)?window1d[ix]:((ix>size34)?window1d[ix-halfSize]:1.0);
					window[index++]=wx*wy;
				}
			}
			if (debugLevel>2){ // one per thread
				(new showDoubleFloatArrays()).showArrays(
						window,
						size,
						size,
						"Window_X"+tileX+"-Y"+tileY
				);
			}
		}
		// create image list;
		int maxImage=0;
		for (int i=0;i<imgPairs.length;i++) for (int n=0;n<2;n++) if (imgPairs[i][n]>maxImage) maxImage=imgPairs[i][n];
		boolean [] imgList= new boolean [maxImage+1]; // which images are used
		for (int i=0;i<imgList.length;i++) imgList[i]=false;
		for (int i=0;i<imgPairs.length;i++) for (int n=0;n<2;n++) imgList[imgPairs[i][n]]=true;
		int plane=cyclopeanTile.getForegroundIndex(); //foregroundIndex should be updated to point next to the  last processed foreground (may be after the last)
//		for (;(plane<cyclopeanTile.getNumberOfPlanes()) && (cyclopeanTile.getPlaneDisparity(plane)<disparityMin);plane++);{
		for (;(plane<cyclopeanTile.getNumberOfPlanes()) && (cyclopeanTile.getPlaneDisparity(plane)>=disparityMin);plane++);{
//			if ((plane>=cyclopeanTile.getNumberOfPlanes()) || (cyclopeanTile.getPlaneDisparity(plane)>disparityMax)){
//				return;// nothing left in this tile
//			}
			double disparity=cyclopeanTile.getPlaneDisparity(plane);
			double [][][] planeStrength=new double [cyclopeanTile.getNumberOfPlanes()][imgPairs.length][];
			//disparityMin..disparityMax should contain not more than 1 plane (verify? enforce?)
			for (int planeOther=plane; planeOther<cyclopeanTile.getNumberOfPlanes() ;planeOther++)  planeStrength[plane]=null;
			double [][] visibility=new double[imgList.length][];
			for (int nImg=0;nImg<imgList.length;nImg++){
				if (imgList[nImg]){
//					double dX=disparity*this.disparityScales[nImg][0];
//					double dY=disparity*this.disparityScales[nImg][1];
					//					double dXYlen=Math.sqrt(dX*dX+dY*dY);
					//					double [] udXY={dX/dXYlen,dY/dXYlen};
					// gets just cumulative data from already processed (closer) layers. Will update visibility later
					visibility[nImg]=	getImageTransparency( // later make option with full color calculation?
							nImg,//ZTile [][] thisZMap,
							allCyclopeanMap,
							(tileX+0.5)*this.overlapStep, //double pX, // result tile center in image pixels
							(tileY+0.5)*this.overlapStep, //double pY, // result tile center in image pixels
//							tileX,
//							tileY,
							size, // resultSize
							
							disparity, //double disparity,
///							filterByForegroundMargin,
///							filterByDisabled?disparityTolearnce:Double.NaN,
							debugLevel);
					if (debugLevel>3){ // +2 for selected tile
						System.out.println("Iterating tileX="+tileX+", tileY="+tileY+" disparity="+disparity+" nImg="+nImg);
					} else visibility[nImg]=null; 
				}
			}
			if (debugLevel>5){
				(new showDoubleFloatArrays()).showArrays(
						visibility, // TODO: add slice names
						paddedSize,
						paddedSize,
						true,
						"PEM_X"+tileX+"-Y"+tileY
				);
			}
			if (refineTileDisparity) {
// read image tiles (upscaled as visibility)s, 				
//				double newDisparity=
				
				// see if it did change - recalculate visibility	
			}

// read 			
// update foregroundIndex, calculate and update cumulative transparency for each image	
// cumnulative transparency - mark, but update from the result tiles, to mamke it thread safe			

		}
	} // end of private void disparitySweepTile	
	
	
	/**
	 * Calculate disparity correction by processing selected pair of images (normally - all pairs)
	 * @param visibility - square visibility array for each image twice the tile period (2*16=32), up-scaled by this.subpixel (normally 4)
	 * @param disparity - original disparity value
	 * @param xc - selection center (in pixels) on virtual cyclopean image 
	 * @param yc
	 * @param refineSubPixel subdivide disparity arrays (and final correlation array) from image pixels (normally 8)
	 * @param maskBlurSigma  blur visibility mask for correlation (~=corrHighPassSigma), default 2.0 pixels (original image)
	 * @param corrHighPassSigma high pass image data before masking for correlation, default 2.0 pixels (original image)
	 * @param zMapCorrMask  bitmap of used channels (+1 - Y, +2 - Cb, +4 - Cr, +8 - Aux) 
	 * @param zMapCorrWeights array of correlation weights - {Y,Cb,Cr,Aux} - will be normalized
	 * @param refinePhaseCoeff  mix phase/normal correlation (0.0 - normal, 1.0 - phase)
	 * @param refineHighPassSigma correlation high pass filter - counts in frequency domain, will not be scaled with this.subpixel 
	 * @param refineLowPassSigma correlation low pass filter - fraction of full frequency range will not be down-scaled with this.subpixel
	 * @param refineCorrMaxDistance maximum disparity correction distance (in image pixels) - 1.5 pix
	 * @param refineCorrThreshold relative correlation strength to enable disparity correction 0.01
	 * @param refineCorrMinPixels minimal number of visible (in both images) pixels to enable disparity correction 
	 * @param window double array of twice tile period squared (normally 32*32=1024), Double.isNaN(window[0]) triggers initialization
	 * @param scaledWindow double array of twice scaled tile period squared (normally 128*128*32=16384), Double.isNaN(window[0]) triggers initialization
	 * @param doubleFHT DoubleFHT instance to be reused for the thread
	 * @param imgPairs pairs of images to use - normally {{0,1},{0,2},{1,2}}
	 * @param imageData original image data - for each image, each channel (0 - alpha, ...4 - Aux) pixels in scanline order: [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
	 * @param imageFullWidth original image width
	 * @param debugLevel debug level
	 * @return corrected disparity value (or same as input disparity if correction is impossible
	 */
	
	public double refineDisparity(
			double [][] visibility,
			double disparity,
			double xc, // "cyclopean" center of the tile
			double yc,
			int refineSubPixel, // 8 (normally > than this.subPixel)
			double maskBlurSigma,
			double corrHighPassSigma,
			int zMapCorrMask, 
//			double [] zMapCorrThresholds,
//			double [] zMapCorrThresholdsRel,
			double [] zMapCorrWeights,
			double refinePhaseCoeff,
			double refineHighPassSigma,
			double refineLowPassSigma,
			double refineCorrMaxDistance,
			double refineCorrThreshold,
			double refineCorrMinPixels,// minimal non-occluded overlapping area to use for disparity correction
			double [] window, // flat top,  now 32x32
			double [] scaledWindow, // 
			DoubleFHT doubleFHT,
			int [][]imgPairs,
			double [][][] imageData, // [img]{Alpha, Y,Cb,Cr, Ext}. Alpha may be null - will handle later?
			int imageFullWidth,
			int debugLevel){
		double [] normCorrWeights=normalizeWeights(zMapCorrMask,zMapCorrWeights);
		boolean [] imgList= getImageList(imgPairs); // which images are used
		int [] channelIndex=getChannelList(zMapCorrMask<<1); // skip alpha, start with Y==1
		int [] channelIndexWeights=getChannelList(zMapCorrMask); // Y==0
		double [][][] filteredImageData = new double [imgList.length][][]; // does not include unused channels
		double [][][] imageTiles=new double [imgList.length][][];
		int size=2*this.overlapStep;
		int scaledSize=this.subpixel*size;
		double scaledMinPixels=refineCorrMinPixels*this.subpixel*this.subpixel;
		double scaledMaskBlurSigma=     maskBlurSigma*this.subpixel;
		double scaledCorrHighPassSigma= corrHighPassSigma*this.subpixel;
		double scaledMaxDistance=refineCorrMaxDistance*this.subpixel;
		double scaledMaxDist2=scaledMaxDistance*scaledMaxDistance;
		if (refineSubPixel<this.subpixel) refineSubPixel=this.subpixel;
		int refineUpsampleExtra=refineSubPixel/this.subpixel;
		int halfDispRange=(int) Math.ceil(2*refineCorrMaxDistance*refineSubPixel);
		for (int nImg=0;nImg<imgList.length;nImg++){
			if (imgList[nImg]){
				filteredImageData[nImg]=new double [channelIndex.length][];
				for (int iChn=0;iChn<channelIndex.length;iChn++){
					filteredImageData[nImg][iChn]=imageData[nImg][channelIndex[iChn]];
				}
			}	else{
				filteredImageData[nImg]=null;
			}
		}
		for (int nImg=0;nImg<imgList.length;nImg++){
			if (imgList[nImg]){
				double xcImg= xc + disparity*this.disparityScales[nImg][0]*this.subpixel;
				double ycImg= yc + disparity*this.disparityScales[nImg][1]*this.subpixel; 
				imageTiles[nImg]=getImageTile(
						xcImg,
						ycImg,
						size, // power of 2 - full size (was half)
						this.subpixel, // power of 2
						filteredImageData[nImg], // only needed channels
						imageFullWidth, 
						window, // should be size*size long;
						scaledWindow, // should be size*size*subPixel*subPixel long;
						doubleFHT, // to reuse tables
						debugLevel);
			} else imageTiles[nImg]=null;
		}
		if (debugLevel>3){
			String [] channelNames={"Alpha","Y","Cb","Cr","Aux"};
			String [] debugTitles=new String [imgList.length*channelIndex.length];
			double [][] debugData=new double [imgList.length*channelIndex.length][];
			for (int iChn=0;iChn<channelIndex.length;iChn++) for (int nImg=0;nImg<imgList.length;nImg++) {
				debugTitles[iChn*imgList.length+nImg]=channelNames[channelIndex[iChn]]+"-"+nImg;
				debugData[iChn*imgList.length+nImg]=imageTiles[nImg][iChn];
			}
			(new showDoubleFloatArrays()).showArrays(
					debugData,
					scaledSize,
					scaledSize,
					true,
					"SIMG-"+IJ.d2s(disparity,1)+"_x"+IJ.d2s(xc,1)+"_y"+IJ.d2s(yc,1),
					debugTitles);
		}
		double [][]  visibilityMask=new double [visibility.length][];
		for (int nImg=0;nImg<imgList.length;nImg++) if (imgList[nImg]){
			// Create a blurred version of visibility mask for each image
			visibilityMask[nImg]=visibility[nImg];
			if (maskBlurSigma>0.0){
				visibilityMask[nImg]=visibility[nImg].clone();
				(new DoubleGaussianBlur()).blurDouble(
						visibilityMask[nImg],
						scaledSize,
						scaledSize,
						scaledMaskBlurSigma,
						scaledMaskBlurSigma,
						0.01);
			}
			// Optionally (normally needed) high-pass each channel of each image, multiply but flat-top window
			for (int iChn=0;iChn<imageTiles[nImg].length;iChn++){
				if (scaledCorrHighPassSigma>0.0){
					double [] loPass=imageTiles[nImg][iChn].clone();
					(new DoubleGaussianBlur()).blurDouble(
							loPass,
							scaledSize,
							scaledSize,
							scaledCorrHighPassSigma,
							scaledCorrHighPassSigma,
							0.01);
//					for (int i=0;i<loPass.length;i++) imageTiles[nImg][iChn][i]=(imageTiles[nImg][iChn][i]-loPass[i])*visibilityMask[nImg][i];
					for (int i=0;i<loPass.length;i++) imageTiles[nImg][iChn][i]=(imageTiles[nImg][iChn][i]-loPass[i])*scaledWindow[i];
				} else {
//					normalizeAndWindow (imageTiles[nImg][iChn], null, true); // only remove DC
					normalizeAndWindow (imageTiles[nImg][iChn], scaledWindow, true); // remove DC, multiply by scaled flat-top window
				}
			}
		}
		if (debugLevel>3){ // debug show visibility masks and high-passed channels for each image 
			String [] channelNames={"Alpha","Y","Cb","Cr","Aux"};
			String [] debugTitles=new String [imgList.length*(channelIndex.length+1)];
			double [][] debugData=new double [imgList.length*(channelIndex.length+1)][];
			for (int iChn=0;iChn<channelIndex.length;iChn++) for (int nImg=0;nImg<imgList.length;nImg++) {
				debugTitles[iChn*imgList.length+nImg]=channelNames[channelIndex[iChn]]+"-"+nImg;
				debugData[iChn*imgList.length+nImg]=imageTiles[nImg][iChn];
			}
			for (int nImg=0;nImg<imgList.length;nImg++) {
				debugTitles[channelIndex.length*imgList.length+nImg]="v-mask-"+nImg;
				debugData[channelIndex.length*imgList.length+nImg]=visibilityMask[nImg];
			}
			(new showDoubleFloatArrays()).showArrays(
					debugData,
					scaledSize,
					scaledSize,
					true,
					"HPM-"+IJ.d2s(disparity,1)+"_x"+IJ.d2s(xc,1)+"_y"+IJ.d2s(yc,1),
					debugTitles);
		}
		double [][]centerCorrs=new double [imgPairs.length][3];
		double [] maskTotal=new double [imgPairs.length];
		double [] disparityArray=null;
		double [] minCorr=new double [imgPairs.length];
		// calculate center correlations for each enabled channel, use blurred visibility masks
		for (int nPair=0;nPair<imgPairs.length;nPair++){
			int nImg=imgPairs[nPair][0];
			int sImg=imgPairs[nPair][1];
			centerCorrs[nPair][0]=0.0;
			centerCorrs[nPair][1]=0.0;
			centerCorrs[nPair][2]=0.0;
			maskTotal[nPair]=0.0;
			double [] pairWeight=new double [scaledWindow.length];
			for (int i=0;i<pairWeight.length;i++) {
				//					double w= scaledWindow[i]*scaledWindow[i]*visibilityMask[nImg][i]*visibilityMask[sImg][i];
				double w= visibilityMask[nImg][i]*visibilityMask[sImg][i];
				pairWeight[i]=w;
				maskTotal[nPair]+=scaledWindow[i]*scaledWindow[i]*w;
			}
			for (int iChn=0;iChn<imageTiles[nImg].length;iChn++){
				double [] centerCorrsChn ={0.0,0.0,0.0};
				for (int i=0;i<pairWeight.length;i++) {
					centerCorrsChn[0]+=pairWeight[i]*imageTiles[nImg][iChn][i]*imageTiles[nImg][iChn][i];
					centerCorrsChn[1]+=pairWeight[i]*imageTiles[sImg][iChn][i]*imageTiles[sImg][iChn][i];
					centerCorrsChn[2]+=pairWeight[i]*imageTiles[nImg][iChn][i]*imageTiles[sImg][iChn][i];
				}
				for (int n=0;n<3;n++) {
					centerCorrs[nPair][n]+=normCorrWeights[channelIndexWeights[iChn]]*centerCorrsChn[n];
				}

			}
			for (int n=0;n<3;n++) {
				centerCorrs[nPair][n]=Math.sqrt(centerCorrs[nPair][n]/maskTotal[nPair]);
			}
			// for each image pair - see if correlation is sufficient, then - correlate and find disparity correction				
			//				disparityArray[nImg][sImg]=null;
			if (maskTotal[nPair]>=scaledMinPixels){
				minCorr[nPair]=1.0;
				for (int n=0;n<centerCorrs[nPair].length;n++) {
					centerCorrs[nPair][n]=Math.sqrt(centerCorrs[nPair][n]/maskTotal[nPair]);
					minCorr[nPair]*=centerCorrs[nPair][n];
				}
				minCorr[nPair]=Math.pow(minCorr[nPair],1.0/centerCorrs[nPair].length); // Use other metrics?
			} else {
				minCorr[nPair]=0.0;
			}
			if (debugLevel>3){ // +2 for selected tile
				System.out.println(
						" maskTotal["+nPair+"]="+maskTotal[nPair]+
						" centerCorrs["+nPair+"][0]="+centerCorrs[nPair][0]+
						" centerCorrs["+nPair+"][1]="+centerCorrs[nPair][1]+
						" centerCorrs["+nPair+"][2]="+centerCorrs[nPair][2]+
						" minCorr["+nPair+"]="+minCorr[nPair]+
						" refineCorrThreshold="+refineCorrThreshold);
			}
			// TODO: now look if the correlation is strong enough to perform correction (minimal of 3?)
			// scaledMinPixels
			// decide if correction is possible?		
			if (minCorr[nPair]>=refineCorrThreshold){
				double [] combinedChnCorr=new double [scaledWindow.length];
				for (int i=0;i<combinedChnCorr.length;i++) combinedChnCorr[i]=0.0;
				double[][]  corr=new double [imageTiles[nImg].length][];
				for (int iChn=0;iChn<imageTiles[nImg].length;iChn++){
					corr[iChn]=doubleFHT.correlate (
							imageTiles[sImg][iChn],
							imageTiles[nImg][iChn],
							refineHighPassSigma, // not scaled - they are in frequency domain
							refineLowPassSigma/this.subpixel,  // scaled to keep the same
							refinePhaseCoeff); // second will be modified
					for (int i=0;i<combinedChnCorr.length;i++) combinedChnCorr[i]+=normCorrWeights[channelIndexWeights[iChn]]*corr[iChn][i];
				}
				if (debugLevel>3){ // for each channel - show per-channel and combine correlation shows imageStack for each image pair
					String [] channelNames={"Alpha","Y","Cb","Cr","Aux"};
					String [] debugTitles=new String [channelIndex.length+1];
					double [][] debugData=new double [channelIndex.length+1][];
					debugTitles[0]="combo";
					debugData[0]=combinedChnCorr;
					for (int iChn=0;iChn<channelIndex.length;iChn++){
						debugTitles[iChn+1]=channelNames[channelIndex[iChn]];
						debugData[iChn+1]=corr[iChn];
					}
					(new showDoubleFloatArrays()).showArrays(
							debugData,
							scaledSize,
							scaledSize,
							true,
							"CORR_"+nImg+"-"+sImg+"_d"+IJ.d2s(disparity,1)+"_x"+IJ.d2s(xc,1)+"_y"+IJ.d2s(yc,1),
							debugTitles);
				}
				// See if maximum is close enough, then subpixel and calculate partial disparityArray - to be combined into a single one later
				double max=0;
				int iMax=0;
				for (int i=0;i<combinedChnCorr.length;i++) if (combinedChnCorr[i]>max){
					max=combinedChnCorr[i];
					iMax=i;
				}
				int ixc=iMax%scaledSize-scaledSize/2;
				int iyc=iMax/scaledSize-scaledSize/2;
				if (debugLevel>3) System.out.println("refineDisparity(): max="+max+" iMax="+iMax+" refineCorrMaxDistance="+refineCorrMaxDistance+" scaledMaxDist2="+scaledMaxDist2+
						" r2="+(ixc*ixc+iyc*iyc)+" r="+Math.sqrt(ixc*ixc+iyc*iyc));
				if ((ixc*ixc+iyc*iyc)<=scaledMaxDist2){ // maximum close enough
					double dX=disparity*(this.disparityScales[sImg][0]-this.disparityScales[nImg][0]);
					double dY=disparity*(this.disparityScales[sImg][1]-this.disparityScales[nImg][1]);
					double dXYlen=Math.sqrt(dX*dX+dY*dY);
					double [] udXY={dX/dXYlen,dY/dXYlen};
					double [] upsampled=doubleFHT.upsample(combinedChnCorr,refineUpsampleExtra);
					int interpolatedSize=scaledSize*refineUpsampleExtra;
					int interpolatedCenter=(interpolatedSize+1)*interpolatedSize/2;
					if (debugLevel>3) { // show upsampled combined for all channels) correlation for enabled pairs
						(new showDoubleFloatArrays()).showArrays(
								upsampled,
								interpolatedSize,
								interpolatedSize,
								"US_N"+nImg+"-"+sImg+"_xc"+IJ.d2s(xc,2)+"_yc"+IJ.d2s(yc,2));
					}
					if (disparityArray==null) { // initilize disparity array (linear, to find maximal value) at the first use
						disparityArray=new double [2*halfDispRange+1];
						for (int i=0;i<disparityArray.length;i++) disparityArray[i]=0.0;
					}
					// build 1-d disparity section by bi-linear interpolation of the 2d correlation results (accumulate for all enabled
					// image pairs
					for (int i=0;i<2*halfDispRange+1;i++){
						double deltaX=udXY[0]*(i-halfDispRange);
						double deltaY=udXY[1]*(i-halfDispRange);
						int iX=(int) Math.floor(deltaX); // zero in the center
						int iY=(int) Math.floor(deltaY);
						deltaX-=iX;
						deltaY-=iY;
						int index00= interpolatedCenter+iY*interpolatedSize+iX;
						int index01=index00+interpolatedSize;
						int index10=index00+1;
						int index11=index01+1;
						disparityArray[i]+= // ACCUMULATE bi-linear interpolated data 
							(upsampled[index00]*(1.0-deltaX)+upsampled[index10]*deltaX)*(1.0-deltaY)+
							(upsampled[index01]*(1.0-deltaX)+upsampled[index11]*deltaX)*     deltaY;
						if (debugLevel>5){
							System.out.println("disparityArray["+i+"]="+disparityArray[i]+
									" deltaX="+deltaX+" deltaY="+deltaY+" iX="+iX+" iY="+iY+" index00="+index00+" index01="+index01+
									" index10="+index10+" index11="+index11);
						}
					}
				}
			} // if (minCorr[nPair]>=refineCorrThreshold){
		} // for (int nPair=0;nPair<imgPairs.length;nPair++)
		
		// find best fitting disparity as a maximum on the disparity array (if it exists)
		double disparCorr=0.0;
		if (disparityArray!=null){
			int iMax=0;
			for (int i=1;i<disparityArray.length;i++){
				if (disparityArray[i]>disparityArray[iMax]) iMax=i;
			}
			disparCorr=((double) (iMax-halfDispRange))/refineSubPixel;
			if (debugLevel>5){
				for (int i=0;i<disparityArray.length;i++){
					System.out.println("combined_disparityArray["+i+"]="+disparityArray[i]);
				}
				System.out.println("\nDisparity correction="+disparCorr+" (iMax="+iMax+"), limit ="+refineCorrMaxDistance +" - VERIFY SIGN IS CORRECT!");
			}
			if (disparCorr<=refineCorrMaxDistance){
				if (debugLevel>3) System.out.println("Old disparity="+disparity+" new disparity="+(disparity+disparCorr));

			} else {
				if (debugLevel>3) System.out.println("Old disparity="+disparity+" unchanged, attempted to correct by "+disparCorr);
				disparCorr=0;
			}
		}
		//just add correction to the original disparity value
		
		return disparity+disparCorr;
	}

	
	
	boolean [] getImageList(int [][]imgPairs){
		int maxImage=0;
		for (int i=0;i<imgPairs.length;i++) for (int n=0;n<2;n++) if (imgPairs[i][n]>maxImage) maxImage=imgPairs[i][n];
		boolean [] imgList= new boolean [maxImage+1]; // which images are used
		for (int i=0;i<imgList.length;i++) imgList[i]=false;
		for (int i=0;i<imgPairs.length;i++) for (int n=0;n<2;n++) imgList[imgPairs[i][n]]=true;
		return imgList;
	}
	int [] getChannelList(int mask){
		int numChannels=0;
		for (int d=mask;d!=0;d>>=1) if ((d&1)!=0)numChannels++;
		int [] channelIndex=new int [numChannels];

		numChannels=0;
		//		int chn=0;
		int d=mask;
		for (int i=0;i<numChannels;i++) {
			if ((d&1)!=0) channelIndex[numChannels++]=i;
			d>>=1;
		}
		return channelIndex;	
	}

	/**
	 * Create 1d disparity array for subcamera from Cyclopean data (will need low-pass filter to merge close maximums)
	 * @param nImg subcamera number
	 * @param cyclopeanMap 2-d array of cyclopean disparity data (uses only disparity/strengths and enabled - no visibility (supposed not to be
	 *                     available yet - easy to modify if needed)
	 * @param pX           result tile center X in original pixels
	 * @param pY           result tile center Y in original pixels
	 * @param window       flat top window encoded in line-scan order. Normally 32 if subpixel==1, has to be the same resolution as visibilityMask (subpixel) 
	 * @param pwr          multiply strength by strength of the cyclopean tile to this power (strength depends on the features and the area )
	 * @param visibilityMask square (encoded in line-scan order) for the visibility of the tile of interest (masked by higher disparity leyaers)
	 * @param subpixel     use resolution higher than pixels (not needed?), should match both window and cisibilityMask
	 * @param disparity    disparity of the layer of interest 
	 * @param disparityStep subdivide disparity pixels by this number when encoding result array
	 * @param debugLevel   debug level (not yet used)
	 * @return             1d array of strengths of disparity values [0] - @infinity, [1] - @disparityStep, ... [disparity/disparityStep] - @disparity
	 *                     Will need low-pass filter to merge close disparities from neighbor cyclopean tiles 
	 */
	
	public double [] getImageDisparityFromCyclopean(
			int nImg,
			CyclopeanTile [][] cyclopeanMap,
			double pX, // center of the result X (in image pixels)
			double pY, // center of the result Y (in image pixels) 
			double [] window, // flat-top tile mask to apply to cyclopean tiles
			double pwr, // multiply strength of the cyclopean tile by strength to this power (strength depends on features and area)
			double [] visibilityMask, // masked visibility at disparity (should have the same size as window)?
			int subpixel, //>1 if window,  wisibilityMask have higher resolution than original pixels 
			double disparity,
			double disparityStep,
			int debugLevel
			){
		if (visibilityMask==null) visibilityMask=window;
		double dSizeVisibility=Math.sqrt(visibilityMask.length);
		double dSizeWindow=Math.sqrt(window.length);
		double dSizeSum=dSizeVisibility+dSizeWindow;
		int sizeVisibility=(int) dSizeVisibility;
		int sizeWindow=(int) dSizeWindow;
		double [] disparityArray=new double [((int) Math.ceil(disparity/disparityStep))+1];
		for (int i=0;i<disparityArray.length;i++) disparityArray[i]=0.0;
		Point2D.Double disparity2DSubpixCenter=new Point2D.Double(
				pX*subpixel,
				pY*subpixel); // in subpixels, center
		Point2D.Double infinity2DSubpixCenter= new Point2D.Double(
				pX*subpixel+this.disparityScales[nImg][0]*disparity*subpixel,
				pY*subpixel+this.disparityScales[nImg][1]*disparity*subpixel); // in subpixels!

		//rectangle including tile TL corners at disparity that can influence result 
		Rectangle2D tilesDisparity2DSubpixCenter=new Rectangle2D.Double(
				disparity2DSubpixCenter.getX(),
				disparity2DSubpixCenter.getY(),
				dSizeSum,
				dSizeSum);
		//rectangle including tile centers at infinity that can influence result
		Rectangle2D tilesInfinity2DSubpixCenter=new Rectangle2D.Double(
				infinity2DSubpixCenter.getX(),
				infinity2DSubpixCenter.getY(),
				dSizeSum,
				dSizeSum);
		Rectangle2D tiles2DSubpixCenter=new Rectangle2D.Double();
		// This rectangle includes all tile centers that may have layers influence the result  
		Rectangle2D.union(
				tilesDisparity2DSubpixCenter,
				tilesInfinity2DSubpixCenter,
				tiles2DSubpixCenter);
		double periodSubpixel=this.overlapStep*subpixel;
		double pXInfinity=pX+this.disparityScales[nImg][0]*disparity; //TODO: check sign!
		double pYInfinity=pX+this.disparityScales[nImg][1]*disparity; //TODO: check sign!
		//rectangle including tile centers at disparity that can influence result 
		Rectangle2D tilesCenterDisparity2D=new Rectangle2D.Double(pXInfinity-(0.5*dSizeSum/subpixel), pYInfinity-(0.5*dSizeSum/subpixel), dSizeSum, dSizeSum);
		//rectangle including tile centers at infinity that can influence result 
		Rectangle2D tilesCenterInfinity2D=new Rectangle2D.Double(pX-(0.5*dSizeSum/subpixel), pY-(0.5*dSizeSum/subpixel), dSizeSum, dSizeSum);
		Rectangle2D tilesCenter2D=new Rectangle2D.Double();
		// This rectangle includes all tile centers that may have layers influence the result  
		Rectangle2D.union(tilesCenterDisparity2D,tilesCenterInfinity2D,tilesCenter2D);
		// convert to tiles, intersect with overall bounds
		Point tilesTL=new Point(
				(int) Math.floor((tilesCenter2D.getMinX()-0.5*dSizeWindow)/periodSubpixel ),
				(int) Math.floor((tilesCenter2D.getMinY()-0.5*dSizeWindow)/periodSubpixel));
		Point tilesBR=new Point(
				(int) Math.ceil((tilesCenter2D.getMaxX()+0.5*dSizeWindow)/periodSubpixel ),
				(int) Math.ceil((tilesCenter2D.getMaxY()+0.5*dSizeWindow)/periodSubpixel));
		Rectangle rectTiles=new Rectangle(
				tilesTL.x,
				tilesTL.y,
				tilesBR.x-tilesTL.x,
				tilesBR.y-tilesTL.y);
		// limit by cyclopeanMap bounds
		rectTiles=rectTiles.intersection(new Rectangle(0,0,cyclopeanMap[0].length,cyclopeanMap.length));
		Rectangle rectResultTile=new Rectangle(0,0,sizeVisibility,sizeVisibility);
		double sumWindow=Double.NaN;
		for (int tY=0;tY<rectTiles.y+rectTiles.height;tY++){
			for (int tX=0;tX<rectTiles.x+rectTiles.width;tX++){
				for (int nPlane=0;
				(nPlane<cyclopeanMap[tY][tX].getNumberOfPlanes()) && (cyclopeanMap[tY][tX].getPlaneDisparity(nPlane)<disparity);
				nPlane++) if (cyclopeanMap[tY][tX].getPlaneEnabled(nPlane)){
					double planeDisparity=cyclopeanMap[tY][tX].getPlaneDisparity(nPlane);
					double planeStrength=Math.pow(cyclopeanMap[tY][tX].getPlaneStrength(nPlane),pwr);
					Point2D.Double tileCenter=new Point2D.Double(
							(tX+0.5)*periodSubpixel-this.disparityScales[nImg][0]*planeDisparity*subpixel, //TODO: check sign!
							(tY+0.5)*periodSubpixel-this.disparityScales[nImg][1]*planeDisparity*subpixel); //TODO: check sign!
					if (tilesCenter2D.contains(tileCenter)){
						int disparityIndex=(int)Math.round(planeDisparity*disparityStep);
						// calculate overlapping masks: - find shifty, shift x (destination to source),  intersection rectangle, scan and accumulate strength
						// calculate subpixel coordinate of the tile (tX,tY) top-left corner in result coordinates (relative to result top-left corner)
						int dx=(sizeVisibility-sizeWindow)/2+(int) Math.round(subpixel*this.disparityScales[nImg][0]*(disparity-planeDisparity));
						int dy=(sizeVisibility-sizeWindow)/2+(int) Math.round(subpixel*this.disparityScales[nImg][1]*(disparity-planeDisparity));
						// wrong - round each separately?
						Rectangle rectThisTile=new Rectangle(dx,dy,sizeWindow,sizeWindow);
						Rectangle rectIntersect=rectResultTile.intersection(rectThisTile);
						if (!rectIntersect.isEmpty()){
							if (Double.isNaN(sumWindow)){
								sumWindow=0.0;
								for (int i=0;i<window.length;i++) sumWindow+=window[i];
							}
							double d=0;
							for (int y=rectIntersect.y;y<(rectIntersect.y+rectIntersect.height);y++){
								int rIndex=y*sizeVisibility+rectIntersect.x;
								int sIndex=(y-dy)*sizeWindow+(rectIntersect.x-dx);
								for (int i=0;i<rectIntersect.width;i++){
									d+=visibilityMask[rIndex++]*window[sIndex++];
								}
								d*=planeStrength/sumWindow;
							}
							disparityArray[disparityIndex]+=d;
						}
					}
				}
			}
		}
		return disparityArray;
	}

	
	/**
	 * Calculate transparency of the closer layers (in the specified disparity range) for the square tile around specified point
	 * using the "cyclopean" tile array
	 * @param nImg           number of image (sub-camera)
	 * @param cyclopeanMap   cyclopean map tiles
	 * @param pX             result tile center, in original image pixels, X
	 * @param pY             result tile center, in original image pixels, Y
	 * @param resultSize     side of the result square (in original image pixels, will be multiplied by this.subpixel (full size - 32, paddedSize - 20)
	 * @param disparity      disparity for which to map transparency
//	 * @param minDisparity   farthest layers to process
//	 * @param maxDisparity   closest layers to process
//	 * @param mergeTolerance maximal distance between layers to consider them the same
	 * @param debugLevel     debug level
	 * @return               square visibility tile in scanline order
	 */
	
	public double [] getImageTransparency(
			int nImg, // --> this.disparityScales
			CyclopeanTile [][] cyclopeanMap,
//			int tX, // current tile
//			int tY,
			double pX, // center of the result X (in image pixels)
			double pY, // center of the result Y (in image pixels) 
			int resultSize, // will be multiplied by this.subpixel, full size - 32, paddedSize - 20
			double disparity, 
//			double minDisparity, // not used?
//			double maxDisparity, // not used?
//			double mergeTolerance,
			int debugLevel
	){
		// assuming to be exactly scaledTilePeriod*scaledTilePeriod  public float [][]imageTransparency; // [image number][pixel number] - cumulative per-image transparency

		int scaledResultSize=this.subpixel*resultSize;
//		int scaledPaddedSize=this.subpixel*this.paddedSize;
		int scaledTilePeriod=this.subpixel*this.overlapStep;
		double [] spxTL={ // relative to  imageTransparency[nImg] top left corner
				(pX+this.disparityScales[nImg][0]*disparity)*this.subpixel-(scaledResultSize)/2,
				(pY+this.disparityScales[nImg][1]*disparity)*this.subpixel-(scaledResultSize)/2};
/* Version using tile tX,tY numbers
  		double [] spxTL={ // relative to  imageTransparency[nImg] top left corner
				tX*scaledTilePeriod+this.disparityScales[nImg][0]*this.subpixel-(scaledResultSize-scaledTilePeriod)/2,
				tY*scaledTilePeriod+this.disparityScales[nImg][1]*this.subpixel-(scaledResultSize-scaledTilePeriod)/2};
*/				
		int [] ispTL={(int) Math.floor(spxTL[0]),(int) Math.floor(spxTL[1])};
		double [] dpxTL={spxTL[0]-ispTL[0],spxTL[1]-ispTL[1]};
		double [] result = new double [resultSize*resultSize];
		double [][]k={
				{(1-dpxTL[0])*(1-dpxTL[1]), (  dpxTL[0])*(1-dpxTL[1])},
				{(1-dpxTL[0])*(  dpxTL[1]), (  dpxTL[0])*(  dpxTL[1])}};
		double [] ones=null;
		int [][] whereY=new int [resultSize+1][2];
		int tileY00=0,tileX00=0;
		for (int y=0;y<resultSize+1;y++){
			int tileY=(y+ispTL[1])/scaledTilePeriod-((y<0)?1:0);
			if (y==0) tileY00=tileY;	
			whereY[y][0]=(y+ispTL[1])-tileY*scaledTilePeriod;
			whereY[y][1]=tileY-tileY00;
		}
		int [][] whereX=new int [resultSize+1][2];
		for (int x=0;x<resultSize+1;x++){
			int tileX=(x+ispTL[0])/scaledTilePeriod-((x<0)?1:0);
			if (x==0) tileX00=tileX;	
			whereX[x][0]=(x+ispTL[1])-tileX*scaledTilePeriod;
			whereX[x][1]=tileX-tileX00;
		}
		double [][][] transparencyTiles=new double[whereY[whereY.length-1][2]][whereX[whereX.length-1][1]][scaledTilePeriod*scaledTilePeriod];
		for (int y=0;y<transparencyTiles.length;y++) for (int x=0;x<transparencyTiles[y].length;x++){
			int tileY=tileY00+y;
			int tileX=tileX00+x;
			transparencyTiles[y][x]=null;
			if ((tileY>=0) && (tileY<cyclopeanMap.length) && (tileX>=0) && (tileX<cyclopeanMap[tileY].length) &&
					(cyclopeanMap[tileY][tileX]!=null)){
				transparencyTiles[y][x]=cyclopeanMap[tileY][tileX].getDoubleTransparency(nImg); // get cumulative transparency stored for each image
			}
			if (transparencyTiles[y][x]==null) {
				if (ones==null){
					ones = new double [scaledTilePeriod*scaledTilePeriod];
					for (int i=0;i<ones.length;i++) ones[i]=1.0;
				}
				transparencyTiles[y][x]=ones;
			}
		}

		for (int y=0;y<resultSize;y++){
			for (int x=0;x<resultSize;y++){
				int index=y*resultSize+x;
				result[index]=                k[0][0]*transparencyTiles[whereY[  y][1]][whereX[  x][1]][whereY[  y][0]*scaledTilePeriod+whereX[x  ][0]];
				if (k[0][1]>0) result[index]+=k[0][1]*transparencyTiles[whereY[  y][1]][whereX[x+1][1]][whereY[  y][0]*scaledTilePeriod+whereX[x+1][0]];
				if (k[1][0]>0) result[index]+=k[1][0]*transparencyTiles[whereY[y+1][1]][whereX[  x][1]][whereY[y+1][0]*scaledTilePeriod+whereX[x  ][0]];
				if (k[1][1]>0) result[index]+=k[1][1]*transparencyTiles[whereY[y+1][1]][whereX[x+1][1]][whereY[y+1][0]*scaledTilePeriod+whereX[x+1][0]];
			}
		}
		return result;
	}

	/**
	 * Create a stack of square tiles from the single image, sub-pixel shifted
	 * @param xc selection center X 
	 * @param yc selection center Y 
	 * @param size result tile size before upsampling (actual size will be size*subPixel
	 * @param subPixel subdivide from original pixels
	 * @param imageData stack of full size image per-channel data
	 * @param imageWidth source image width
	 * @param window double array of size*size, window[0] should be Double.NaN to be initialized
	 * @param iWindowUpSample inverse values of the scaled window - should be double[size*subPixel*size*subPixel] 
	 * @param doubleFHT Double FHT class instance to be reused by the same thread multiple times
	 * @param debugLevel debug level
	 * @return stack of double[size*subPixel*size*subPixel], outer pixels are unreliable  
	 */
	public double [][] getImageTile(
			double xc,
			double yc,
			int size, // power of 2 - full size (was half)
			int subPixel, // power of 2
			double [][] imageData, // only needed channels
			int imageWidth,
			double [] window, // should be size*size long;
			double [] windowUpSample, // should be size*size*subPixel*subPixel long;
			DoubleFHT doubleFHT, // to reuse tables
			int debugLevel
	){
		if (Double.isNaN(windowUpSample[0])) {
			int upsapledSize=size*subPixel;
			double[] iWindowUpSample1d=doubleFHT.getHamming1d(upsapledSize); // Hamming
			int index=0;
			int quarterSize=size/4;
			int size34=3*size/4;
			for (int iy=0;iy<size;iy++) {
				double wy=(iy<quarterSize)?iWindowUpSample1d[iy]:((iy>size34)?iWindowUpSample1d[iy-size]:1.0);
				for (int ix=0;ix<size;ix++) {
					double wx=(ix<quarterSize)?iWindowUpSample1d[ix]:((ix>size34)?iWindowUpSample1d[ix-size]:1.0);
					windowUpSample[index++]=wx*wy;
				}
			}
		}
		if (Double.isNaN(window[0])) {
			double[] window1d=doubleFHT.getHamming1d(size); // Hamming
			int index=0;
			int quarterSize=size/4;
			int size34=3*size/4;
			for (int iy=0;iy<size;iy++) {
				double wy=(iy<quarterSize)?window1d[iy]:((iy>size34)?window1d[iy-size]:1.0);
				for (int ix=0;ix<size;ix++) {
					double wx=(ix<quarterSize)?window1d[ix]:((ix>size34)?window1d[ix-size]:1.0);
					window[index++]=wx*wy;
				}
			}
		}
		if (debugLevel>4){
			System.out.print("getImageTile()");

			if (debugLevel>3){
				(new showDoubleFloatArrays()).showArrays(
						window,
						size,
						size,
						"window");
			}
		}
		double [][] result=new double [imageData.length][];
		double [] sliceTile=new double [size*size];
		int ixc=(int) Math.round(xc);
		int iyc=(int) Math.round(yc);
		double dxc=xc-ixc;
		double dyc=yc-iyc;
		boolean shift=(dxc!=0.0) || (dyc!=0.0);
		boolean scale=(subPixel!=1);
		for (int chn =0; chn<imageData.length;chn++) { 
			sliceTile= getSelection(
					imageData[chn], // source image/channel slice
					sliceTile,
					size, //int width,
					size, //int height,
					imageWidth,
					ixc , //xc,
					iyc); //yc);
			if (debugLevel>4) {
				(new showDoubleFloatArrays()).showArrays(
						sliceTile,
						size,
						size,
				"slice-getImageTile");
			}
			// normalize and save DC
			double dc=0.0;
			if (shift || scale){
				dc=normalizeAndWindowGetDC (sliceTile, window); //windowInterpolation
				if (debugLevel>4) {
					(new showDoubleFloatArrays()).showArrays(
							sliceTile,
							size,
							size,
							"slice-normalized-dc"+dc);
				}

				result[chn]=doubleFHT.shift(sliceTile, subpixel, -dxc, -dyc);
				if (debugLevel>4) {
					(new showDoubleFloatArrays()).showArrays(
							sliceTile,
							size,
							size,
					"slice-shifted");
				}
				for (int i=0;i<result[chn].length;i++) result[chn][i]=result[chn][i]/windowUpSample[i]+dc;
			} else {
				if (debugLevel>3) System.out.println("No shift or scale is needed");
				result[chn]=sliceTile.clone();
			}
			if (debugLevel>3) System.out.println("getImageTile() dc="+dc);
		}
		return result;
	}
	
	
   	public double [] getSelection(
			double [] imageSlice, // one image/channel slice
			double [] selection, // null or array to reuse
			int width,
			int height,
			int fullWidth,
			int xc,
			int yc){
		int length=width*height;
		int fullHeight=imageSlice.length/fullWidth;
		if (selection ==null) selection = new double[length];
			int y0=yc-height/2;
			int x0=xc-width/2;
    		for (int iy=0;iy<height;iy++) {
    			int srcY=iy+y0;
    			boolean oob=(srcY<0) || (srcY>=fullHeight);
    			for (int ix=0;ix<width;ix++){
    				int oIndex=iy*width+ix;
					int srcX=x0+ix;
					if (oob ||(srcX<0) || (srcX>=fullWidth)) {
						if (oIndex>=selection.length) System.out.println("\ngetSelection(imageSlice["+imageSlice.length+"],selection["+selection.length+"]"+
								","+fullWidth+","+width+","+height+","+xc+","+yc+") ix="+ix+" iy="+iy+" oIndex="+oIndex);
						selection[oIndex]=0.0;
					} else {
    						selection[oIndex]=imageSlice[srcY*fullWidth+srcX];
					}
    			}
    		}
		return selection;
	}
	public double[] normalizeAndWindow (double [] pixels, double [] windowFunction, boolean removeDC) {
		int j;
		if (pixels==null) return null;
		double s=0.0,s0=0.0;
		if (removeDC) {
			for (j=0;j<pixels.length;j++){
				s+=pixels[j]*windowFunction[j];
				s0+=windowFunction[j];
				
			}
			s/=s0;
		}
		if (windowFunction!=null) for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
		else for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s);
		return pixels;
	}

	public double normalizeAndWindowGetDC (double [] pixels, double [] windowFunction) {
		int j;
		if (pixels==null) return 0;
		double s=0.0,s0=0.0;
		for (j=0;j<pixels.length;j++){
			s+=pixels[j]*windowFunction[j];
			s0+=windowFunction[j];
		}
		s/=s0;
		for (j=0;j<pixels.length;j++) pixels[j]=(pixels[j]-s)*windowFunction[j];
		return s;
	}




	public double [] normalizeWeights(int mask,double [] weights){
		double [] maskedWeights=weights.clone();
		for (int i=0;i<maskedWeights.length;i++) if ((mask & (1<<i))==0) maskedWeights[i]=0.0;
		return normalizeWeights(maskedWeights);
	}
	public double [] normalizeWeights(double [] weights){
		double [] normalizedWeights=new double [weights.length];
		double sum=0.0;
		for (int i=0;i<weights.length;i++) sum+=weights[i];
		for (int i=0;i<weights.length;i++) normalizedWeights[i]=(sum==0.0)?0.0:(weights[i]/sum);
		return normalizedWeights;
	}

	public int    getDisparityPoints() {return this.impDisparity.getWidth()/this.tilesX;}
	public double [] getTileDisparity(
			int itx,
			int ity,
			float [] data){ //npair <0 - use center data
		int width=this.impDisparity.getWidth();
		int disparityPoints=getDisparityPoints();
		double [] result=new double [disparityPoints];
		int index=ity*width+itx*disparityPoints;
		for (int i=0;i<disparityPoints;i++) result[i]=data[index++]; //oob -146583
		return result;
	}

	public void setTileDisparity(
			double []disparityArray,
			int itx,
			int ity,
			float [] data){ //npair <0 - use center data
		int width=this.impDisparity.getWidth();
		int disparityPoints=getDisparityPoints();
		int index=ity*width+itx*disparityPoints;
		for (int i=0;i<disparityPoints;i++) {
			if ((index>=data.length) || (i>=disparityArray.length)){
				System.out.println("setTileDisparity(disparityArray,"+itx+","+ity+", data): index="+index+
						" data.length="+data.length+" i="+i+" disparityArray.length="+disparityArray.length+" disparityPoints="+disparityPoints+" width="+width);
			}
			data[index++]=(float) disparityArray[i]; //oob 20375037
		}
	}
   	
	public void setupCyclopeanTile(
			Rectangle woi, // in tiles - may be
//			final int nImg,
			final int maxNumber,
			final double minFirst,
			final double minAbsolute,
			final double minRelative,
			final double mergeMax,
			final int overlap,
			final double zMapMinForeground,
			final int threadsMax,
			final boolean showProgress,
			final int debugLevel){
		if (this.centerPixels==null) {
			String msg="Centered disparity data is not defined";
			IJ.showMessage("Error",msg);
			throw new IllegalArgumentException (msg);
		}
//		if (this.syntheticPixels==null) moveDisparityFromCenter(threadsMax,showProgress,debugLevel);
		final float [] thisCenterPixels= this.centerPixels;
		if (woi==null) woi=new Rectangle (0, 0, this.tilesX,this.tilesY);
		this.zMapWOI=new Rectangle();
		this.zMapWOI.x=(woi.x>=0)?woi.x:0;
		this.zMapWOI.y=(woi.y>=0)?woi.y:0;
		this.zMapWOI.width= (((woi.x+woi.width) <=this.tilesX)?(woi.x+woi.width): this.tilesX)-this.zMapWOI.x;
		this.zMapWOI.height=(((woi.y+woi.height)<=this.tilesY)?(woi.y+woi.height):this.tilesY)-this.zMapWOI.y;
		final Rectangle zMapWOI=this.zMapWOI;
		if (this.cyclopeanMap==null) initCyclopeanMap();
		final CyclopeanTile [][] thisCyclopeanMap=this.cyclopeanMap;
		final int tilesX=this.zMapWOI.width;
		final int tiles=this.zMapWOI.width*this.zMapWOI.height;
		if (debugLevel>2) System.out.println("setupZMap() woi.x="+woi.x+" woi.y="+woi.y+
				" woi.width="+woi.width+" woi.height="+woi.height);
		if (debugLevel>2) System.out.println("setupZMap() zMapWOI.x="+zMapWOI.x+" zMapWOI.y="+zMapWOI.y+
				" zMapWOI.width="+zMapWOI.width+" zMapWOI.height="+zMapWOI.height);
		this.paddedSize=this.overlapStep+2*overlap;
		this.innerMask=(BitSet) new BitSet(this.paddedSize*this.paddedSize);
		for (int i=0;i<this.overlapStep;i++) for (int j=0;j<this.overlapStep;j++){
			this.innerMask.set((overlap+i)*this.paddedSize+overlap+j);
		}
//	private int  []   borderMask=null; // will be provided to zMap +1:top+2:bottom+8:left+16:right (to prevent roll over when iterating around
		int tileSize=2*this.overlapStep;
		this.borderMask=new int [tileSize*tileSize];
		for (int i=0;i<tileSize;i++) for (int j=0;j<tileSize;j++){
			this.borderMask[i*tileSize+j]=((i==0)?5:0)+((i==(tileSize-1))?2:0)+((j==0)?0x50:0)+((j==(tileSize-1))?0x20:0);
		}
		

		final Thread[] threads = newThreadArray(threadsMax);
   		final AtomicInteger tileIndexAtomic     = new AtomicInteger(0);
   		if (showProgress) IJ.showProgress(0.0);
   		if (showProgress) IJ.showStatus("Setting up cyclopeanMap ...");
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				public void run() {
   					for (int tile=tileIndexAtomic.getAndIncrement(); tile<tiles;tile=tileIndexAtomic.getAndIncrement()){
   						setupCyclopeanTile(
   								zMapWOI.x+tile%tilesX, //int tileX,
   								zMapWOI.y+tile/tilesX, //int tileY,
   								centerPixels,
   								thisCyclopeanMap,
   			        			maxNumber,
   			        			minFirst,
   			        			minAbsolute,
   			        			minRelative,
   			        			mergeMax,
   			        			overlap,
   			        			zMapMinForeground,
   			        			debugLevel);
   						if (showProgress){
   							final int finalTile=tile;
   							SwingUtilities.invokeLater(new Runnable() {
   								public void run() {
   									IJ.showProgress(finalTile,tiles);
   								}
   							});
   						}
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
   		IJ.showProgress(1.0);
	}


	private void setupCyclopeanTile(
			int tileX,
			int tileY,
			float [] thisDisparityTiles,
			CyclopeanTile [][] thisCyclopeanMap,
			int maxNumber,
			double minFirst,
			double minAbsolute,
			double minRelative,
			double mergeMax,
			int overlap,
			double zMapMinForeground,
			int debugLevel
	){
		if (debugLevel>2) System.out.println("setupCyclopeanTile("+tileX+","+tileY+",...)");
		boolean debugThisTile=debugLevel>3;
		thisCyclopeanMap[tileY][tileX]=new CyclopeanTile();
		CyclopeanTile cyclopeanTile=thisCyclopeanMap[tileY][tileX];
		double [] correlationSection = getTileDisparity(tileX, tileY, thisDisparityTiles);
		//					double threshold=minAbsolute;
		boolean [] isMax =new boolean[correlationSection.length];
		int iMax=-1;
		int numMax=0;
		for (int i=0;i<correlationSection.length;i++){
			if ((correlationSection[i]>=minAbsolute) &&
					((i==0)||(correlationSection[i]>=correlationSection[i-1])) &&
					((i==(correlationSection.length-1)) || (correlationSection[i]>=correlationSection[i+1]))){
				isMax[i]=true; // local max above absolute threshold
				if ((numMax==0) || (correlationSection[i]>correlationSection[iMax])) iMax=i;
				numMax=1;
			} else {
				isMax[i]=false;
			}
		}
		//						double mergeMax=1.5; //merge maximums for interpolation if closer than mergeMax;
		int iMergeMax=(int) Math.round(mergeMax/this.disparityPerEntry);
		if (debugThisTile){
			System.out.println("mergeMax="+mergeMax+" iMergeMax="+iMergeMax);
		}
		if ((numMax>0) && (correlationSection[iMax]<minFirst)) numMax=0; // no maximums if the largest is below thershold
		// for now - now maximums above threshold - use infinity (disparity 0)
		if (numMax==0) {
			isMax[0]=true; // infinity
			numMax=1;
			iMax=0;
		} else {
			numMax=0;
//			double threshold=correlationSection[iMax]*minRelative; // no relative filtering here - will be done later
			for (int i=0;i<correlationSection.length;i++) if (isMax[i]){
//				if (correlationSection[i]>=threshold) numMax++;
				numMax++;
				//else isMax[i]=false;
			}
		}
		if (numMax>maxNumber) numMax=maxNumber; // limit to specified number of correlation maximums to subpixel
		int [] maxIndices=new int [numMax];
		maxIndices[0]=iMax;
		isMax[iMax]=false;
		int maxNum;
		for (maxNum=1;maxNum<maxIndices.length;maxNum++){
			// merge previous one if possible
			int nearDown=-1, nearUp=-1;
			for (int i=0;i<=iMergeMax;i++){
				if ((maxIndices[maxNum-1]-i)<0) break;
				if (isMax[maxIndices[maxNum-1]-i]){
					nearDown = maxIndices[maxNum-1]-i;
					break;
				}
			}
			for (int i=0;i<=iMergeMax;i++){
				if ((maxIndices[maxNum-1]+i)>=isMax.length) break;
				if (isMax[maxIndices[maxNum-1]+i]){
					nearUp = maxIndices[maxNum-1]+i;
					break;
				}
			}
			int n=1+((nearDown>=0)?1:0)+((nearUp>=0)?1:0);
			if (n>1){
				double s=maxIndices[maxNum-1]+((nearDown>=0)?nearDown:0)+((nearUp>=0)?nearUp:0);
				int iMerged=(int) Math.round(s/n);

				if (debugThisTile){
					System.out.println("Merging close maximums: "+maxIndices[maxNum-1]+" with "+
							((nearDown>=0)?nearDown:"")+" "+((nearUp>=0)?nearUp:"")+" to "+iMerged);
				}
				maxIndices[maxNum-1]=iMerged;
				isMax[iMerged]=false;
				if (nearDown>=0) isMax[nearDown]=false;
				if (nearUp>=0)   isMax[nearUp]=false;
			}
			iMax=-1;
			for (int i=0;i<correlationSection.length;i++) {
				if (isMax[i]&&((iMax<0) || (correlationSection[i]>=correlationSection[iMax]))) iMax=i;
			}
			if (iMax<0) break; // no more maximums
			isMax[iMax]=false;
			maxIndices[maxNum]=iMax;
		}
		if (debugThisTile){
			System.out.println("List maximums on correlation section");
			for (int n=0;n<maxNum;n++) {
				System.out.println(n+" "+maxIndices[n]+"("+(this.disparityPerEntry*maxIndices[n])+" pix) "+correlationSection[maxIndices[n]]);
			}
		}
		//TODO: Always add infinity?
		//TODO: use quadratic interpolation?
		// reorder maximums by decreasing disparity (they are now ordered by strength
		for (boolean ordered=false;!ordered;){
			ordered=true;
			for (int i=0;i<(maxIndices.length-1);i++) if (maxIndices[i]<maxIndices[i+1]){
				ordered=false;
				int tmp=maxIndices[i];
				maxIndices[i]=maxIndices[i+1];
				maxIndices[i+1]=tmp;
			}
		}
		// is there infinity?
		boolean noInfinity=(maxIndices.length==0) || (maxIndices[maxIndices.length-1]!=0); 
		cyclopeanTile.numPlanes=noInfinity?(maxIndices.length+1):maxIndices.length;
		cyclopeanTile.maximums=new double [cyclopeanTile.numPlanes][2];
		for (int i=0;i<maxIndices.length;i++){
			cyclopeanTile.maximums[i][0]= (float) (maxIndices[i]*this.disparityPerEntry); // disparity in pixels
			cyclopeanTile.maximums[i][1]= (float) correlationSection[maxIndices[i]];      // strength
		}
		cyclopeanTile.maximums[cyclopeanTile.numPlanes-1][0]= 0.0f; // disparity in pixels
		cyclopeanTile.maximums[cyclopeanTile.numPlanes-1][1]= (float) correlationSection[0];      // strength
		//overlap
		cyclopeanTile.overlap=overlap;
		cyclopeanTile.size=this.overlapStep;
		cyclopeanTile.subpixel=this.subpixel;
		cyclopeanTile.reset(zMapMinForeground, minAbsolute, minRelative);
		cyclopeanTile.setInnerMask(this.innerMask);
		cyclopeanTile.setBorderMask(this.borderMask);
	}

	private int imagesToNPair(int firstImage, int secondImageIndex){
		for (int i=0;i<this.imagePairIndices.length;i++) if ((this.imagePairIndices[i][0]==firstImage)&& (this.imagePairIndices[i][1]==secondImageIndex)) return i;
		return -1;
	}
	public void initCyclopeanMap(){
		if (this.cyclopeanMap==null){
			this.cyclopeanMap=new CyclopeanTile[this.tilesY][this.tilesX];
			for (int tileY=0;tileY<this.tilesY;tileY++) for (int tileX=0;tileX<this.tilesX;tileX++){
				this.cyclopeanMap[tileY][tileX]=null;
			}
		}
	}
	public Rectangle pixelsToTilesWOI(
			Rectangle pixelWOI){
		Rectangle woi=new Rectangle (
				pixelWOI.x/this.overlapStep,
				pixelWOI.y/this.overlapStep,
				(pixelWOI.x+pixelWOI.width -1)/this.overlapStep - pixelWOI.x/this.overlapStep +1, 
				(pixelWOI.y+pixelWOI.height-1)/this.overlapStep - pixelWOI.y/this.overlapStep +1);
		return woi;
	}
	
	
	public class CyclopeanTile{
		public int numPlanes; //=0;
		public int size;      // 16
		public int overlap; //=0;
		public int subpixel;
		public double foregroundThreshold; // current foreground threshold (may reconsider if occlusion happens)
		public double minAbsolute=0;
		public double minRelative=0;
		public int foregroundIndex; // current foreground index;
		public double [][] maximums; //={}; // {position, strength} in the order of decreasing disparity // now should always include infinity
		public float [] zmap; //=null; // square, centered in the center of the tile - may have margins - final z-value for each pixel - will not work with transparency
		public boolean [] enabledPlane;
//		public BitSet [] enabledPixels; //=null;
//		public BitSet [] certainPixels; //=null;
		public float [][] likely; //=null; 
		public float [][] unlikely; //=null; 
		public float [][] auxData=null;
		private BitSet innerMask=null;
		private int []borderMask=null;

// TODO: replace - use "global" (cyclopean) layrered opaqueness and cumulative per-image ones (from the largest disparity to the current one)
// Maybe - later add cumulative pixel color also
// another idea for thin wire / tree branches. Guess color by increasing contrast and/or using the same object over different background, then shrink
// the width (sub-pixel) to have the same average color as on the images.
		
		public float [][] globPlaneOpaque; // for each plane, each pixel inside window 1.0 for opaque, 0.0 - for transparent, no overlap, no - with overlap (used before merge
		public double [][] globPlaneEdgeDisparity; // for each plane - top,left,right,bottom - NaN - use center, otherwise - linear interpolate from the center
		
		// after the pass - see if any is NaN, while neighbor exists - then fill in both this and neighbor.
		// Races will not harm, as the result will be the same 
		// what if two tiles are worked on simultaneously?
		// first - correct disparity, then add globPlaneDisparity[plane]
		
		public AtomicBoolean  planeModified; // set during first pass - do not use it yet when calculating visibility (to prevent races), reset at second pass
		public int lastModifiedPlane; // valid with planeModified set
		// no need for atomic? just "-1" for lastModifiedPlane? Just set it before modifying disparity (it only matters if disparity correction made it farther,
		// so it can be counted again? Or is it OK?
		
		
		public AtomicBoolean [] planeSet; //set after the disparity is corrected (second pass, resets planeModified) and globPlaneOpaque is calculated (maybe does not need to be atomic?)
		
		public AtomicBoolean [][] planeEdgeSet; // [plane][edge], where edge is 0/1 either right or bottom
		// or one of 4 - E,SE,S,SW (wit diagonals?
		
		
		// after the whole pass, try to merge free-hanging sides - with the larger disparity - within tolerance, with smaller disparities (possible after
		// disparity correction - any?
		// synchronize to prevent races when two neigb. tiles are corrected simultaneously?
/*
 *  No need to synchronize - on the		
 */
		
		
		
		public float [][]imageTransparency; // [image number][pixel number] - cumulative per-image transparency, mapped to zero disparity, no overlap 

		public double [] getDoubleTransparency(int nImg){
			if ((this.imageTransparency==null) || (this.imageTransparency[nImg]==null)) return null;
			double [] result = new double [this.imageTransparency.length];
			for (int i=0;i<result.length;i++) result[i]=this.imageTransparency[nImg][i];
			return result;
		}
		public int getNumberOfPlanes(){return this.maximums.length;}
		public int getForegroundPlane(){return this.foregroundIndex;}
		public double getPlaneDisparity (int plane){return this.maximums[plane][0];}
		public void setPlaneDisparity (double disparity, int plane){this.maximums[plane][0]=disparity;}
		public double getPlaneStrength (int plane){return this.maximums[plane][1];}
		public boolean getPlaneEnabled (int plane){return this.enabledPlane[plane];}
		public int getSize(){return this.size;}
		public int getPaddedSize(){return this.size+2*this.overlap;}
		public int getOverlap(){return this.overlap;}
		public int getPaddedLength(){return (this.size+2*this.overlap)*(this.size+2*this.overlap);}
		public void setInnerMask(BitSet mask){
			this.innerMask=mask;
		}
		public void setBorderMask(int [] mask){
			this.borderMask=mask;
		}
		public int [] getBorderMask(){return this.borderMask;}
		public int getForegroundIndex() {return this.foregroundIndex;}
		public void setForegroundPlane(int plane){
			this.foregroundIndex=plane;
			/*
			if (this.maximums[plane][1]<this.foregroundThreshold) this.foregroundThreshold=this.maximums[plane][1];
			this.enabledPlane[plane]=true;
			*/
		}
		public double [] getDisparityListDescending(
				double minDisparity,
				double maxDisparity){
			int numResult=0;
			int plane;
			for (plane=0;plane<this.numPlanes;plane++){
				if (this.maximums[plane][0]>maxDisparity) continue;
				if (this.maximums[plane][0]<minDisparity) break;
				if (this.enabledPlane[plane]) numResult++;
			}
			double [] result=new double[numResult];
			numResult=0;
			for (plane=0;plane<this.numPlanes;plane++){
				if (this.maximums[plane][0]>maxDisparity) continue;
				if (this.maximums[plane][0]<minDisparity) break;
				if (this.enabledPlane[plane]) result[numResult++]=this.maximums[plane][0];
			}
			return result;
		}
		
		
		/**
		 * Filter pixels by occlusion of the foreground and disabled in the current plane
		 * @param disparityFG - threshold to consider pixel being in front
		 * @param disparity - target disparity
		 * @param disparityTolerance - consider "this" as being withing disparityTolerance of disparity. NaN - do not filter by this
		 * @return
		 */
		/*
		public double [] getVisibility(
				double disparityFG,
				double disparity,
				double disparityTolerance){
    		return getEnabledNonOccluded(
    				disparityFG,
    				disparity,
    				disparityTolerance,
    				0);
		}
		public double [] getVisibility(
				double disparityFG,
				double disparity,
				double disparityTolerance,
				int debugLevel){
			if (debugLevel>3) System.out.println(" ---- getEnabledNonOccluded("+disparityFG+","+
				disparity+","+
				disparityTolerance+","+
				debugLevel+")");
			int paddedSize=this.size+2*this.overlap;
			int paddedLength=paddedSize*paddedSize;
			BitSet nonOccludedBits=new BitSet(paddedLength);
			nonOccludedBits.set(0,paddedLength);
			for (int plane=0;(plane<getNumberOfPlanes()) && (getPlaneDisparity(plane)>disparityFG);plane++) if (this.enabledPlane[plane]) {
				if (this.certainPixels[plane]!=null) nonOccludedBits.andNot(this.certainPixels[plane]);
    			if (debugLevel>3) System.out.println(" ------ plane="+plane+" cumul. nonOccludedBits.cardinality()="+nonOccludedBits.cardinality());
			}
			if (!Double.isNaN(disparityTolerance)){
				BitSet enabledBits=new BitSet(paddedLength);
				for (int plane=0;plane<getNumberOfPlanes();plane++)
					if (this.enabledPlane[plane] && (Math.abs(getPlaneDisparity(plane)-disparity)<=disparityTolerance)) {
						if (this.enabledPixels[plane]!=null) enabledBits.or(this.enabledPixels[plane]);
						else enabledBits.set(0,paddedLength); // all enabled
            			if (debugLevel>3) System.out.println(" ------ plane="+plane+" cumul. enabledBits.cardinality()="+enabledBits.cardinality());
					}
				nonOccludedBits.and(enabledBits);
    			if (debugLevel>3) System.out.println(" ------ result nonOccludedBitscardinality()="+nonOccludedBits.cardinality());
			}
			boolean [] nonOccluded=new boolean[paddedLength];
			for (int i=0;i<paddedLength;i++) nonOccluded[i]=nonOccludedBits.get(i);
			return nonOccluded;
		}
*/
		
    	public void reset(
    			double minForeground,
    	    	double minAbsolute,
    	    	double minRelative){
    		this.numPlanes=this.maximums.length;
    		this.globPlaneOpaque=new float[this.numPlanes][];
    		this.imageTransparency=null;
    		this.enabledPlane=new boolean[this.numPlanes];
    		for (int i=0;i<this.numPlanes;i++){
    			this.globPlaneOpaque[i]= null;
    			this.enabledPlane[i]=true;
    		}
    		this.zmap=null;
    		setMinCorrelations(minAbsolute, minRelative,false);
    		setForeGround(minForeground);
    		this.likely=new float[this.numPlanes][]; // maybe use later
    		for (int i=0;i<this.likely.length;i++)this.likely[i]=null;
    		this.unlikely=new float[this.numPlanes][]; // maybe use later
    		for (int i=0;i<this.unlikely.length;i++)this.unlikely[i]=null;
    	}
    	public void setMinCorrelations(
    	    	double minAbsolute,
    	    	double minRelative,
    	    	boolean keepFG){
    		if (!Double.isNaN(minAbsolute))   this.minAbsolute=minAbsolute;
    		if (!Double.isNaN(minRelative))   this.minRelative=minRelative;
    		double aMax=0;
    		for (int i=0;i<this.numPlanes;i++) if (this.maximums[i][1]>aMax) aMax=this.maximums[i][1];
    		aMax=Math.max(aMax*this.minRelative,this.minAbsolute);
    		for (int i=0;i<this.numPlanes;i++)	this.enabledPlane[i]=this.maximums[i][1]>=aMax;
    		if (keepFG && (this.foregroundIndex<this.enabledPlane.length) && (this.foregroundIndex>=0)) {
    			this.enabledPlane[this.foregroundIndex]=true;
    		}
    	}
    	public void setForeGround(
    			double minForeground){
    		if (!Double.isNaN(minForeground)) this.foregroundThreshold=minForeground;
    		for (this.foregroundIndex=0;
  		  (this.foregroundIndex<this.maximums.length) && (!this.enabledPlane[this.foregroundIndex] || (this.maximums[this.foregroundIndex][1]<this.foregroundThreshold))
  		  ;this.foregroundIndex++);
    	}

		
		
		

	}


	public class Photometric{
		public String [] channelNames={"Y","Cb","Cr","Aux"};
		public double  [][] valueLimits=null; // lo/high limits for each channel
		public int subdivAverage=256; 
		public int subdivHalfDifference=128; // 257
		//		public int subdivVariance = 256;
		public double smoothVarianceSigma=10.0;
		public double scaleVariance=3.0; // sigma along difference (vertical axis) will be scaleVariance*variance for this average value
		public int numImages=0;
		public double [][][] standardVariance=null;// for each image each channel - average variance (among 9 neighbors) for different values
		public double [][] averageVariance=null;
		public double [][][][] matchingQuality=null; // for each image, each second image index, each channel [subdivAverage*subdivDifference] <=1.0 values
		// horizontal average, vertical - difference "2-1" (for "1-2" flip vertical)
		public int histogramSize=1000;
		public double getAverageVariance(int nImg,int chn){
			return this.averageVariance[nImg][chn];
		}
		public double [] initStaging(){
			double [] staging=new double[this.subdivAverage*(2*this.subdivHalfDifference+1)];
			for (int i=0;i<staging.length;i++) staging[i]=0.0;
			return staging;
		}
		public void addStagingSample(
				double []staging,
				int chn,
				double weight,
				double v1,
				double v2){
			double min=this.valueLimits[chn][0];
			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
			double av=0.5*(v1+v2);
			int countAvg=(int) Math.round((av-min)/step);
			if (countAvg<0) countAvg=0;
			if (countAvg>=this.subdivAverage) countAvg=this.subdivAverage-1;
			int countDiff= (int)Math.round(0.5*(v2-v1))+subdivHalfDifference;
			if (countDiff<0) countDiff=0;
			if (countDiff>2*this.subdivHalfDifference) countDiff=2*this.subdivHalfDifference;
			staging[countDiff*this.subdivAverage+countAvg]+=weight;
		}
		public void addStagingSample(
				double []staging,
				int chn,
				double weight,
				double v1,
				double v2,
				// reduce weight depending on difference, scale to measured variance
				double scaleVariance,
				double kLocal, // 0 - use global varaiance, 1.0 - use local 
				int nImg1, // first image number
				int nImg2, // second image number
				double var1, // first image variance
				double var2, // second image variance
				boolean debug){
			double min=this.valueLimits[chn][0];
			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
			double av=0.5*(v1+v2);
			int countAvg=(int) Math.round((av-min)/step);
			if (countAvg<0) countAvg=0;
			if (countAvg>=this.subdivAverage) countAvg=this.subdivAverage-1;
			double diff=0.5*(v2-v1);
			int countDiff= (int)Math.round(diff/step)+subdivHalfDifference;
			if (countDiff<0) countDiff=0;
			if (countDiff>2*this.subdivHalfDifference) countDiff=2*this.subdivHalfDifference;
			if (debug){
				System.out.println("addStagingSample(...,"+chn+","+weight+","+v1+","+v2+","+scaleVariance+","+
						kLocal+","+nImg1+","+nImg2+","+var1+","+var2);
				System.out.println(" +++ min="+min+" step="+step+" av="+av+" diff="+diff+" countAvg="+countAvg+" countDiff="+countDiff);
			}
			if (scaleVariance>0.0){
				if (kLocal<1.0){ // mix local variance with average over all image
					int count1=(int) Math.round((v1-min)/step);
					if (count1<0) count1=0;
					if (count1>=this.subdivAverage) count1=this.subdivAverage-1;
					int count2=(int) Math.round((v2-min)/step);
					if (count2<0) count2=0;
					if (count2>=this.subdivAverage) count2=this.subdivAverage-1;
					if (kLocal<0) kLocal=0;
					var1=kLocal*var1+(1.0-kLocal)*this.standardVariance[nImg1][chn][count1];
					var2=kLocal*var2+(1.0-kLocal)*this.standardVariance[nImg2][chn][count2];
				}
				double sigma=scaleVariance*Math.sqrt(var1*var1);
				weight*=Math.exp(-diff*diff/(2.0*sigma*sigma))/sigma/Math.sqrt(2*Math.PI);
				if (debug)	System.out.println(" +++ sigma="+sigma+" weight="+weight);
			}
			staging[countDiff*this.subdivAverage+countAvg]+=weight;
			if (debug)	System.out.println("     staging["+(countDiff*this.subdivAverage+countAvg)+"]="+staging[countDiff*this.subdivAverage+countAvg]);

		}
		public double getStrength(
				double []staging,
				int chn,
				double v1,
				double v2){
			double min=this.valueLimits[chn][0];
			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
			double av=0.5*(v1+v2);
			int countAvg=(int) Math.floor((av-min)/step);
			int countAvg1=(int) Math.ceil((av-min)/step);
			double dx=(av-min)-step*countAvg;
			if (countAvg<0) {
				countAvg=0;
				countAvg1=0;
				dx=0.0;
			}
			if (countAvg1>=this.subdivAverage){
				countAvg=this.subdivAverage-1;
				countAvg1=this.subdivAverage-1;
				dx=0.0;
			}
			double diff=0.5*(v2-v1);
			int countDiff= (int)Math.floor(diff/step)+this.subdivHalfDifference;
			int countDiff1= (int)Math.ceil(diff/step)+this.subdivHalfDifference;
			double dy=diff-step*(countDiff-this.subdivHalfDifference);
			if (countDiff<0) {
				countDiff=0;
				countDiff1=0;
				dy=0.0;
			}
			if (countDiff1>2*this.subdivHalfDifference){
				countDiff=2*this.subdivHalfDifference;
				countDiff1=2*this.subdivHalfDifference;
				dy=0.0;
			}
			if ((countDiff1*this.subdivAverage+countAvg1)>=staging.length){
				System.out.println("BUG: getStrength() countAvg="+countAvg+" countAvg1="+countAvg1+
						" countDiff="+countDiff+" countDiff1="+countDiff1+" staging.length="+staging.length);
			}
			return
			staging[countDiff *this.subdivAverage+countAvg ]*(1-dy)*(1-dx)+
			staging[countDiff *this.subdivAverage+countAvg1]*(1-dy)*(  dx)+
			staging[countDiff1*this.subdivAverage+countAvg ]*(  dy)*(1-dx)+
			staging[countDiff1*this.subdivAverage+countAvg1]*(  dy)*(  dx);
		}

		public void blurStaging(
				double []staging,
				int chn,
				double sigma){
			if (sigma<=0) return;

			double min=this.valueLimits[chn][0];
			double step=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);

			(new DoubleGaussianBlur()).blurDouble(
					staging,
					this.subdivAverage,
					2*this.subdivHalfDifference+1,
					sigma/step,
					sigma/step,
					0.01);
		}
		public  Photometric(
				double [][][] images,
				int imageWidth,
				int margins,
				double ignoreFraction,
				int subdivAverage, 
				int subdivHalfDifference,
				//        		int subdivVariance,
				double smoothVarianceSigma,
				double scaleVariance,
				int debugLevel
		){
			if (margins<1) margins=1;
			this.subdivAverage=subdivAverage; 
			this.subdivHalfDifference=subdivHalfDifference;
			//    		this.subdivVariance = subdivVariance;
			this.smoothVarianceSigma=smoothVarianceSigma;
			this.scaleVariance=scaleVariance; // sigma along difference (vertical axis) will be scaleVariance*variance for this average value
			this.numImages=images.length;
			this.valueLimits=new double[channelNames.length][];
			int [] dirs1={1,imageWidth+1,imageWidth,imageWidth-1,-1,-imageWidth-1,-imageWidth,-imageWidth+1,0};
			this.standardVariance=new double[this.numImages][this.channelNames.length][];
			this.averageVariance=new double[this.numImages][this.channelNames.length];
			for (int i=0;i<this.standardVariance.length;i++) for (int j=0;j<this.standardVariance[0].length;j++){
				this.standardVariance[i][j]=null;
				this.averageVariance[i][j]=0.0;
			}

			for (int chn=0;chn<this.valueLimits.length;chn++){
				this.valueLimits[chn]=null;
				double min=Double.NaN, max=Double.NaN;
				for (int nImg=0;nImg<this.numImages;nImg++) if ((images[nImg]!=null) && (images[nImg][chn]!=null)){
					double [] image=images[nImg][chn];
					int imageHeight=image.length/imageWidth;
					if (Double.isNaN(min)){
						min=image[margins*imageWidth+margins];
						max=min;
					}
					if (debugLevel>2){
						System.out.println("nImg="+nImg+" chn="+chn+
								" min0="+IJ.d2s(min,3)+" max0="+IJ.d2s(max,3)); 
					}
					for (int y=margins;y<(imageHeight-margins);y++) for (int x=margins;x<(imageWidth-margins);x++){
						if (image[y*imageWidth+x]<min) {
							min=image[y*imageWidth+x];
							if (debugLevel>2)	System.out.println("y="+y+" x="+x+" min="+IJ.d2s(min,3)); 

						}
						if (image[y*imageWidth+x]>max) {
							max=image[y*imageWidth+x];
							if (debugLevel>2)	System.out.println("y="+y+" x="+x+" max="+IJ.d2s(max,3)); 
						}
					}
					if (debugLevel>2){
						System.out.println("nImg="+nImg+" chn="+chn+
								" min="+IJ.d2s(min,3)+" max="+IJ.d2s(max,3)); 
					}

				}
				int [] histogram=new int [this.histogramSize];
				for (int i=0;i<histogram.length;i++) histogram[i]=0;
				double step=(max-min)/(this.histogramSize-0.0001);
				if (debugLevel>2){
					System.out.println(	" min="+IJ.d2s(min,3)+" max="+IJ.d2s(max,3)+" step="+IJ.d2s(step,6)); 
				}

				for (int nImg=0;nImg<this.numImages;nImg++) if ((images[nImg]!=null) && (images[nImg][chn]!=null)){
					double [] image=images[nImg][chn];
					int imageHeight=image.length/imageWidth;
					for (int y=margins;y<(imageHeight-margins);y++) for (int x=margins;x<(imageWidth-margins);x++){
						int index=(int) Math.floor((image[y*imageWidth+x]-min)/step);
						if (index<0) index=0;
						if (index>=this.histogramSize) index=this.histogramSize-1;
						histogram[index]++; //java.lang.ArrayIndexOutOfBoundsException: 1005
					}
				}
				int totalNum=0;
				for (int i=0;i<histogram.length;i++)totalNum+= histogram[i];
				int ignoreNumPix= (int)(Math.floor(totalNum*ignoreFraction));
				if (debugLevel>2){
					System.out.println(	" totalNum="+totalNum+" ignoreNumPix="+ignoreNumPix); 
				}

				this.valueLimits[chn]=new double[2];
				int num=0;
				for (int i=0;i<histogram.length;i++){
					num+=histogram[i];
					if (debugLevel>2)	System.out.println("---  "+num+" min+"+i+"*step="+IJ.d2s(min+i*step,3)); 
					if (num>=ignoreNumPix){
						this.valueLimits[chn][0]=min+i*step;
						if (debugLevel>2){
							System.out.println("i="+i+" min+i*step="+IJ.d2s(min+i*step,3)); 
						}
						break;
					}
				}
				num=0;
				for (int i=histogram.length-1;i>=0;i--){
					num+=histogram[i];
					if (debugLevel>2)	System.out.println("---  "+num+" min+"+i+"*step="+IJ.d2s(min+i*step,3)); 
					if (num>=ignoreNumPix){
						this.valueLimits[chn][1]=min+i*step;
						if (debugLevel>2){
							System.out.println("i="+i+" min+i*step="+IJ.d2s(min+i*step,3)); 
						}
						break;
					}
				}
				if (debugLevel>1){
					System.out.println("Channel '"+this.channelNames[chn]+
							"' min="+IJ.d2s(this.valueLimits[chn][0],3)+" max="+IJ.d2s(this.valueLimits[chn][1],3)); 
				}
				// calculatye variance for each emage/channel
				for (int nImg=0;nImg<this.numImages;nImg++) if ((images[nImg]!=null) && (images[nImg][chn]!=null)){
					double [] image=images[nImg][chn];
					int imageHeight=image.length/imageWidth;
					this.standardVariance[nImg][chn]=new double[this.subdivAverage];
					int [] samples=new int [this.subdivAverage];
					int totalSamples=0;
					for (int i=0;i<this.subdivAverage;i++){
						this.standardVariance[nImg][chn][i]=0.0;
						samples[i]=0;
					}
					double stepVar=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/(this.subdivAverage-1);
					double minVar=this.valueLimits[chn][0];
					for (int y=margins;y<(imageHeight-margins);y++) for (int x=margins;x<(imageWidth-margins);x++){
						int index=y*imageWidth+x;
						double  S1=0.0;
						double  S2=0.0;
						for (int d=0;d<dirs1.length;d++){
							int indexOther=index+dirs1[d];
							S1+=image[indexOther];
							S2+=image[indexOther]*image[indexOther];
						}
						double v2=(S2-(S1*S1)/dirs1.length)/dirs1.length; //dirs1.length;
						double avg=S1/dirs1.length;
						int count=(int) Math.round((avg-minVar)/stepVar);
						if (count<0) count=0;
						if (count>=this.subdivAverage) count=this.subdivAverage-1;
						this.standardVariance[nImg][chn][count]+=v2; // squared
						samples[count]++;
						this.averageVariance[nImg][chn]+=v2;
						totalSamples++;
					}
					for (int i=0;i<this.subdivAverage;i++){
						if (samples[i]>0) this.standardVariance[nImg][chn][i]=Math.sqrt(this.standardVariance[nImg][chn][i]/samples[i]);
						else this.standardVariance[nImg][chn][i]=0.0;
						//            			if ((debugLevel>1) && (nImg==0) && (chn==0)){
						if ((debugLevel>2) && (nImg==0)){
							System.out.println("samples["+i+"]="+samples[i]+
									" this.standardVariance["+nImg+"]["+chn+"]["+i+"]="+IJ.d2s(this.standardVariance[nImg][chn][i],3)); 
						}

					}
					this.averageVariance[nImg][chn]=Math.sqrt(this.averageVariance[nImg][chn]/totalSamples);
					for (int i=0;i<this.subdivAverage;i++) if (samples[i]==0){
						int j=i+1;
						for (;(j<this.subdivAverage) && (samples[j]==0);j++);
						if (j==this.subdivAverage) j--;
						int i0=i-1;
						if (i0<0) i0=0;
						if (samples[i0]==0) this.standardVariance[nImg][chn][i0]=this.standardVariance[nImg][chn][j];
						if (samples[j]==0)  this.standardVariance[nImg][chn][j]=this.standardVariance[nImg][chn][i0];
						if (i0<j){
							double a=(this.standardVariance[nImg][chn][j]-this.standardVariance[nImg][chn][i0])/(j-i0);
							for (int k=0;k<(j-i0);k++){
								this.standardVariance[nImg][chn][i0+k]=this.standardVariance[nImg][chn][i0]+a*k;
							}
							i=j+1;
						}
					}
					// fill gaps (if any)
					if (this.smoothVarianceSigma>0.0){
						(new DoubleGaussianBlur()).blur1Direction(
								this.standardVariance[nImg][chn], //double [] pixels,
								this.subdivAverage, //int        width,
								1, //int       height,
								this.smoothVarianceSigma, //double     sigma,
								0.01, //double   accuracy,
								true //boolean xDirection: true - horizontal
						);
					}
				}            			
			}
			// now cteate matchingQuality arrays for each image pair/channel
			//    		this.matchingQuality=new double [this.numImages][this.numImages-1][this.valueLimits.length][this.subdivAverage*this.subdivDifference];
			this.matchingQuality=new double [this.numImages][this.numImages-1][][]; //[this.subdivAverage*this.subdivDifference];
			int subdivDifference=2*this.subdivHalfDifference+1;
			int matchingQualityLength=this.subdivAverage*subdivDifference;
			for (int nImg=0;nImg<this.numImages;nImg++) for (int sIndex=0;sIndex<(this.numImages-1);sIndex++){
				int sImg=(sIndex>=nImg)?(sIndex+1):sIndex;
				if ((images[nImg]!=null) && (images[sImg]!=null)){
					this.matchingQuality[nImg][sIndex]=new double [this.valueLimits.length][];
					for (int chn=0;chn<this.valueLimits.length;chn++){
						if ((images[nImg][chn]!=null) && (images[sImg][chn]!=null)){
							this.matchingQuality[nImg][sIndex][chn]=new double [matchingQualityLength];
							double diffStep=(this.valueLimits[chn][1]-this.valueLimits[chn][0])/this.subdivHalfDifference;
							double k=0.5*diffStep*diffStep/(this.scaleVariance*this.scaleVariance);
							for (int averageIndex=0;averageIndex<this.subdivAverage;averageIndex++){
								double variance2=this.standardVariance[nImg][chn][averageIndex]*this.standardVariance[sImg][chn][averageIndex];
								//		public double scaleVariance=3.0; // sigma along difference (vertical axis) will be scaleVariance*variance for this average value
								double a=k/variance2;
								for (int i=0;i<=this.subdivHalfDifference;i++){
									double d=Math.exp(-a*i*i);
									this.matchingQuality[nImg][sIndex][chn][(this.subdivHalfDifference-i)*this.subdivAverage+averageIndex]=d;
									this.matchingQuality[nImg][sIndex][chn][(this.subdivHalfDifference+i)*this.subdivAverage+averageIndex]=d;
								}
							}
						} else {
							this.matchingQuality[nImg][sIndex][chn]=null;
						}
					}
				} else {
					this.matchingQuality[nImg][sIndex]=null;
				}
			}
		}
		public void showmatchingQuality(){
			int subdivDifference=2*this.subdivHalfDifference+1;
			int matchingQualityLength=this.subdivAverage*subdivDifference;
			double [] zero = new double [matchingQualityLength];
			int numPairs=this.numImages*(this.numImages-1);
			double [][] debugData=new double [numPairs*this.channelNames.length][];
			String [] titles=new String [numPairs*this.channelNames.length];
			int index=0;
			for (int nImg=0;nImg<this.numImages;nImg++) for (int sIndex=0;sIndex<(this.numImages-1);sIndex++){
				int sImg=(sIndex>=nImg)?(sIndex+1):sIndex;
				for (int chn=0;chn<this.channelNames.length;chn++){
					titles[index]=this.channelNames[chn]+"-"+nImg+"-"+sImg;
					debugData[index++]=(this.matchingQuality[nImg][sIndex][chn]!=null)?this.matchingQuality[nImg][sIndex][chn]:zero;
				}

			}
			(new showDoubleFloatArrays()).showArrays(
					debugData,
					this.subdivAverage,
					subdivDifference,
					true,
					"MQ-"+IJ.d2s(this.scaleVariance,2)+"_"+IJ.d2s(this.smoothVarianceSigma,2),
					titles);
		}
	} // end of private class Photometric
	/* Create a Thread[] array as large as the number of processors available.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private Thread[] newThreadArray(int maxCPUs) {
		int n_cpus = Runtime.getRuntime().availableProcessors();
		if (n_cpus>maxCPUs)n_cpus=maxCPUs;
		return new Thread[n_cpus];
	}
/* Start all given threads and wait on each of them until all are done.
	 * From Stephan Preibisch's Multithreading.java class. See:
	 * http://repo.or.cz/w/trakem2.git?a=blob;f=mpi/fruitfly/general/MultiThreading.java;hb=HEAD
	 */
	private static void startAndJoin(Thread[] threads)
	{
		for (int ithread = 0; ithread < threads.length; ++ithread)
		{
			threads[ithread].setPriority(Thread.NORM_PRIORITY);
			threads[ithread].start();
		}

		try
		{   
			for (int ithread = 0; ithread < threads.length; ++ithread)
				threads[ithread].join();
		} catch (InterruptedException ie)
		{
			throw new RuntimeException(ie);
		}
	}

}
