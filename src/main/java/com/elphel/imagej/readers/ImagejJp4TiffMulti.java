/**
 ** -----------------------------------------------------------------------------**
 ** ImagejJp4TiffMulti.java
 **
 ** Uses loci.format compatible readers for Elphel 8/16 bpp monochrome Tiff and
 ** JP4 files to read/parse camera files in a single url read operation by
 ** buffering camera data with Location.mapFile()
 **
 **
 ** Copyright (C) 2019 Elphel, Inc.
 **
 ** -----------------------------------------------------------------------------**
 **
 **  ImagejJp4TiffMulti.java is free software: you can redistribute it and/or modify
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

import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ij.ImagePlus;
import loci.formats.ClassList;
import loci.formats.FormatException;

public class ImagejJp4TiffMulti {
	private static final int MAX_THREADS =   100;
	private static final Logger LOGGER =     LoggerFactory.getLogger(ClassList.class);
	private final ImagejJp4Tiff [] imagejJp4Tiff = new ImagejJp4Tiff [MAX_THREADS];


	public ImagePlus [] getMultiImages(
			final String [] urls,
			final ImagePlus [] imps,
//			final boolean telemetry,
			final boolean scale,
			final String std) throws IOException, FormatException  // std - include non-elphel properties with prefix std
	{
//		final ImagePlus [] imps = new ImagePlus [urls.length];
   		final Thread[] threads = newThreadArray(MAX_THREADS);
   		final AtomicInteger indxAtomic = new AtomicInteger(0);
   		final AtomicInteger threadAtomic = new AtomicInteger(0);
   		for (int ithread = 0; ithread < threads.length; ithread++) {
   			threads[ithread] = new Thread() {
   				@Override
				public void run() {
   					int threadIndx = threadAtomic.getAndIncrement();
   					if (imagejJp4Tiff[threadIndx] == null ) {
   						imagejJp4Tiff[threadIndx] = new ImagejJp4Tiff();
   					}
   					for (int indx = indxAtomic.getAndIncrement(); indx < urls.length; indx = indxAtomic.getAndIncrement())
   					{
   						try {
							imps[indx] = imagejJp4Tiff[threadIndx].readTiffJp4(urls[indx], scale, std);
						} catch (IOException e) {
							LOGGER.error("getMultiImages IOException " + urls[indx]);
						} catch (FormatException e) {
							LOGGER.error("getMultiImages FormatException " + urls[indx]);
						}
   					}
   				}
   			};
   		}
   		startAndJoin(threads);
		return imps;
	}




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
