package com.elphel.imagej.gpu;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.Channels;
import java.nio.channels.WritableByteChannel;

import com.elphel.imagej.cameras.CLTParameters;
import com.elphel.imagej.tileprocessor.GeometryCorrection;
import com.elphel.imagej.tileprocessor.QuadCLT;

import ij.ImagePlus;

public class ExportForGPUDevelopment {
	public static void processCLTQuadCorrPairForGPU(
			String                                         save_prefix,
			QuadCLT                                        quadCLT,
			ImagePlus []                                   imp_quad,
			CLTParameters                                  clt_parameters){
		int tilesX = quadCLT.tp.getTilesX();
		int tilesY = quadCLT.tp.getTilesY();
		double z_correction =  clt_parameters.z_correction;
		if (clt_parameters.z_corr_map.containsKey(quadCLT.getImageName())){ // not used in lwir
			z_correction +=clt_parameters.z_corr_map.get(quadCLT.getImageName());
		}
		double disparity_corr = (z_correction == 0) ? 0.0 : quadCLT.geometryCorrection.getDisparityFromZ(1.0/z_correction);
		int [][]    tile_op = quadCLT.tp.setSameTileOp(clt_parameters,  clt_parameters.tile_task_op, -1); // debugLevel); 		
		double [][] disparity_array = quadCLT.tp.setSameDisparity(clt_parameters.disparity); // 0.0); // [tp.tilesY][tp.tilesX] - individual per-tile expected disparity
		
		//		  int mcorr_sel = Correlation2d.corrSelEncode(clt_parameters.img_dtt, getNumSensors());
		TpTask[] tp_tasks = GpuQuad.setTasks(
				quadCLT.getNumSensors(),       // final int                      num_cams,
				clt_parameters.transform_size, // final int                      transform_size,
				disparity_array,               // final double [][]	           disparity_array,  // [tilesY][tilesX] - individual per-tile expected disparity
				disparity_corr,                // final double                   disparity_corr,
				tile_op,                       // final int [][]                 tile_op,          // [tilesY][tilesX] - what to do - 0 - nothing for this tile
				quadCLT.geometryCorrection,    // final GeometryCorrection       geometryCorrection,
				100); // threadsMax);                   // final int                      threadsMax)       // maximal number of threads to launch
		double [][][]       port_xy = new double [tilesX*tilesY][][];
		for (int nTile = 0; nTile < port_xy.length; nTile++) {
			port_xy[nTile] = tp_tasks[nTile].getDoubleXY();
		}
		if ((save_prefix != null) && (save_prefix != "")) {

			String kernel_dir = save_prefix+"clt/";
			File kdir = new File(kernel_dir);
			kdir.mkdir();
			//		boolean [][] what_to_save = {{false,false,true}, {false,false,true}};
			boolean [][] what_to_save = {{true,true,true}, {true,true,true}};
			try {
				saveFloatKernels(
						kernel_dir + (quadCLT.isAux()?"aux":"main"), // String file_prefix,
						(what_to_save[0][0]?quadCLT.getCLTKernels(): null), // double [][][][][][] clt_kernels, // null
						(what_to_save[0][1]?quadCLT.image_data:      null),
						(what_to_save[0][2]?port_xy:                 null), // double [][][]       port_xy,
						true);
			} catch (IOException e) {
				System.out.println("Failed to save flattened kernels tp "+kernel_dir);
				// TODO Auto-generated catch block
				e.printStackTrace();
			} // boolean transpose);

			// make it same length of 16 sensors (for fixed-size struct gc in GPU kernel code 
			GeometryCorrection ext_gc = quadCLT.getGeometryCorrection().expandSensors(GPUTileProcessor.NUM_CAMS) ;
			
			try {
//				quadCLT.getGeometryCorrection().saveFloatsGPU(kernel_dir + (quadCLT.isAux()?"aux":"main"));
				ext_gc.saveFloatsGPU(kernel_dir + (quadCLT.isAux()?"aux":"main"));
			} catch (IOException e) {
				System.out.println("Failed to save geometry correction data (float) to "+kernel_dir);
				e.printStackTrace();
			}
			ext_gc.getCorrVector().getRotMatricesDbg();  
			ext_gc.getCorrVector().getRotDeriveMatricesDbg();
//			quadCLT.getGeometryCorrection().getCorrVector().getRotMatricesDbg();       // ?
//			quadCLT.getGeometryCorrection().getCorrVector().getRotDeriveMatricesDbg(); // ?
		}
		
		
		
	}
	public static void saveFloatKernels(String file_prefix,
			double [][][][][][] clt_kernels,
			double [][][]       image_data,
			double [][][]       port_xy,
			boolean transpose) throws IOException {
		if (clt_kernels != null) {
			for (int chn = 0; chn < clt_kernels.length; chn++) {
				String kern_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel";
				String offs_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel_offsets";
				FileOutputStream fos = new FileOutputStream(kern_path);
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				int float_buffer_size = clt_kernels[chn].length * clt_kernels[chn][0].length* clt_kernels[chn][0][0].length * 4 * 64;
				ByteBuffer bb = ByteBuffer.allocate(float_buffer_size * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							for (int p = 0; p < 4; p++) {
								double [] pa = clt_kernels[chn][col][ty][tx][p];
								for (int i0 = 0; i0 < 64; i0++) {
									int i;
									if (transpose) {
										i = ((i0 & 7) << 3) + ((i0 >>3) & 7);
									} else {
										i = i0;
									}
									bb.putFloat((float)pa[i]);
								}
							}
						}
					}
				}
				bb.flip();
				channel.write(bb);
				dos.close();

				fos = new FileOutputStream(offs_path);
				dos = new DataOutputStream(fos);
				channel = Channels.newChannel(dos);
				float_buffer_size = clt_kernels[chn][0].length * clt_kernels[chn][0].length* clt_kernels[chn][0][0].length * 4 * clt_kernels[chn][0][0][0][4].length;
				bb = ByteBuffer.allocate(float_buffer_size * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							double [] pa = clt_kernels[chn][col][ty][tx][4];
							for (int i = 0; i < pa.length; i++) {
								bb.putFloat((float)pa[i]);
							}
						}
					}
				}
				bb.flip();
				channel.write(bb);
				dos.close();
			}
		}

		if (image_data != null) {
			for (int chn = 0; chn < image_data.length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".bayer";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				ByteBuffer bb = ByteBuffer.allocate(image_data[chn][0].length * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int i = 0; i <  image_data[chn][0].length; i++) {
					double d = 0;
					for (int c = 0; c < image_data[chn].length; c++) {
						d += image_data[chn][c][i];
					}
					bb.putFloat((float) d);
				}
				bb.flip();
				channel.write(bb);
				dos.close();
			}
		}
		if (port_xy != null) {
			for (int chn = 0; chn < port_xy[0].length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".portsxy";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				WritableByteChannel channel = Channels.newChannel(dos);
				ByteBuffer bb = ByteBuffer.allocate(port_xy.length * 2 * 4);
				bb.order(ByteOrder.LITTLE_ENDIAN);
				bb.clear();
				for (int i = 0; i <  port_xy.length; i++) {
					bb.putFloat((float) (port_xy[i][chn][0])); // x-offset
					bb.putFloat((float) (port_xy[i][chn][1])); // y-offset
				}
				bb.flip();
				channel.write(bb);
				dos.close();
			}
		}
	}

	public static void saveFloatKernelsBigEndian(String file_prefix, // never used
			double [][][][][][] clt_kernels,
			double [][][]       image_data,
			double [][][]       port_xy,
			boolean transpose) throws IOException {
		if (clt_kernels != null) {
			for (int chn = 0; chn < clt_kernels.length; chn++) {
				String kern_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel";
				String offs_path = file_prefix+"_chn"+chn+(transpose?"_transposed":"")+".kernel_offsets";
				FileOutputStream fos = new FileOutputStream(kern_path);
				DataOutputStream dos = new DataOutputStream(fos);
				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							for (int p = 0; p < 4; p++) {
								double [] pa = clt_kernels[chn][col][ty][tx][p];
								for (int i0 = 0; i0 < 64; i0++) {
									int i;
									if (transpose) {
										i = ((i0 & 7) << 3) + ((i0 >>3) & 7);
									} else {
										i = i0;
									}
									dos.writeFloat((float)pa[i]);
								}
							}
						}
					}
				}
				dos.close();
				fos = new FileOutputStream(offs_path);
				dos = new DataOutputStream(fos);

				for (int ty = 0; ty <  clt_kernels[chn][0].length; ty++) {
					for (int tx = 0; tx <  clt_kernels[chn][0][ty].length; tx++) {
						for (int col = 0; col <  clt_kernels[chn].length; col++) {
							double [] pa = clt_kernels[chn][col][ty][tx][4];
							for (int i = 0; i < pa.length; i++) {
								dos.writeFloat((float)pa[i]);
							}
						}
					}
				}

				dos.close();
			}
		}

		if (image_data != null) {
			for (int chn = 0; chn < image_data.length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".bayer";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				for (int i = 0; i <  image_data[chn][0].length; i++) {
					dos.writeFloat((float) (image_data[chn][0][i] + image_data[chn][1][i] + image_data[chn][2][i]));
				}
				dos.close();
			}
		}
		if (port_xy != null) {
			for (int chn = 0; chn < port_xy[0].length; chn++) {
				String img_path =  file_prefix+"_chn"+chn+".portsxy";
				FileOutputStream fos = new FileOutputStream(img_path);
				DataOutputStream dos = new DataOutputStream(fos);
				for (int i = 0; i <  port_xy.length; i++) {
					dos.writeFloat((float) (port_xy[i][chn][0])); // x-offset
					dos.writeFloat((float) (port_xy[i][chn][1])); // y-offset
				}
				dos.close();
			}

		}

	}
	
}
