/**
 * Copyright (C) 2018 Elphel, Inc.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

import org.tensorflow.Graph;
import org.tensorflow.Session;
import org.tensorflow.Tensor;
import org.tensorflow.Tensors;
import org.tensorflow.TensorFlow;
import org.tensorflow.SavedModelBundle;
import org.tensorflow.OperationBuilder;
import org.tensorflow.Shape;
import org.tensorflow.Output;
import org.tensorflow.Operation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**

Note 0: articles & examples:

- https://divis.io/2017/11/enterprise-tensorflow-1/
- https://divis.io/2018/01/enterprise-tensorflow-code-examples/
- https://divis.io/2018/01/enterprise-tensorflow-2-saving-a-trained-model/
- https://divis.io/2018/01/enterprise-tensorflow-3-loading-a-savedmodel-in-java/
- https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/
- https://www.programcreek.com/java-api-examples/?api=org.tensorflow.SavedModelBundle
- https://github.com/imagej/imagej-tensorflow
- simple: http://www.riptutorial.com/tensorflow/example/32154/load-and-use-the-model-in-java-

Note 1: How to feed:
 a. https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/
 b. https://github.com/tensorflow/tensorflow/blob/master/tensorflow/java/src/main/java/org/tensorflow/examples/LabelImage.java



Note 2: https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/

...
Two types of objects will need manual closing for proper resource handling:
Sessions and Tensors
...



*/

public class TensorflowExamplePlugin
{

    public final static String EXPORTDIR = "/home/oleg/GIT/python3-imagej-tiff/data_sets/tf_data_5x5_main_13_heur/exportdir";
    // tf.saved_model.tag_constants.SERVING = "serve"
    public final static String SERVING = "serve";
    

    public static void run()
    {
        System.out.println("TensorflowExamplePlugin run");
        try {
            main();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /** 
     * From https://github.com/DIVSIO/tensorflow_java_cli_example/blob/master/src/main/java/divisio/example/tensorflow/cli/RunRegression.java
     */
    
	/**
	 * wraps a single float in a tensor
	 * @param f the float to wrap
	 * @return a tensor containing the float
	 */
	private static Tensor<Float> toTensor(final float f, final Collection<Tensor<?>> tensorsToClose) {
		final Tensor<Float> t = Tensors.create(f);
		if (tensorsToClose != null) {
			tensorsToClose.add(t);
		}
		return t;
	}

	private static Tensor<Float> toTensor2DFloat(final float [][] f, final Collection<Tensor<?>> tensorsToClose) {
		final Tensor<Float> t = Tensors.create(f);
		if (tensorsToClose != null) {
			tensorsToClose.add(t);
		}
		return t;
	}
	
	private static Tensor<Integer> toTensor1DInt(final int [] f, final Collection<Tensor<?>> tensorsToClose) {
		final Tensor<Integer> t = Tensors.create(f);
		if (tensorsToClose != null) {
			tensorsToClose.add(t);
		}
		return t;
	}
	
	private static void closeTensors(final Collection<Tensor<?>> ts) {		
		for (final Tensor<?> t : ts) {
			try {
				t.close();
			} catch (final Exception e) {
				System.err.println("Error closing Tensor.");
				e.printStackTrace();
			}
		}
		ts.clear();
	}
	
	
	
    public static void main() throws Exception {

        final Graph smpb;

        // init for variable?
        float [][] rv_stage1_out = new float[78408][32];
        
        // from infer_qcds_01.py
        float [][] img_corr2d = new float[78408][324];
        float [][] img_target = new float[78408][  1];
        int     [] img_ntile  = new   int[78408];
        
        // init ntile for testing?
        for(int i=0;i<img_ntile.length;i++){
        	img_ntile[i] = i;
        }
        
        final SavedModelBundle bundle = SavedModelBundle.load(EXPORTDIR,SERVING); 
        
        final List<Tensor<?>> tensorsToClose = new ArrayList<Tensor<?>>(5);
        
        System.out.println("OK");
        
        try {
        	
        	System.out.println("S0:");
        	// read Variable info test
        	Operation opr = bundle.graph().operation("rv_stage1_out");
        	System.out.println(opr.toString());
        	
        	System.out.println("S1:");

        	//opr = bundle.graph().operation("rv_stageY_out");
        	//System.out.println(opr.toString());
        	
        	// init variable via constant?
        	//Tensor<Float> tsr = toTensor2DFloat(rv_stage1_out, tensorsToClose);
        	/*
        	Output builder_init = bundle.graph()
        			                    .opBuilder("Const", "rv_stage1_out_init")
        			                    .setAttr  ("dtype", tsr.dataType())
        			                    .setAttr  ("value", tsr)
        			                    .build()
        			                    .output(0);
        	*/
        	//System.out.println(builder_init);
        	
        	// variable
        	//OperationBuilder builder2 = bundle.graph().opBuilder("Variable", "rv_stage1_out_extra_variable");
        											  //.addInput(builder_init);
        	
        	//builder2.
        	//bundle.graph().opBuilder("Assign", "Assign/" + builder2.op().name()).addInput(variable).addInput(value).build().output(0);
        	
        	//Tensor<Float> tensorVal = tsr;
        	//Output oValue = bundle.graph().opBuilder("Const", "rv_stage1_out_2").setAttr("dtype", tensorVal.dataType()).setAttr("value", tensorVal).build().output(0);
        	//System.out.println(oValue);
        	//Output oValue = bundle.graph().opBuilder("Variable", "rv_stage1_out").setAttr("value", tensorVal).build().output(0);
        	//bundle.graph().opBuilder("Assign", "Assign/rv_stage1_out").setAttr("value", tsr).build();
        	
        	System.out.println("Stage 0.1");
        	//bundle.session().runner().fetch("rv_stageY_out").run();
        	System.out.println("Stage 0.2");
        	bundle.session().runner().fetch("rv_stage1_out").run();
        	System.out.println("Stage 1");
        	
        	// stage 1        	
        	bundle.session().runner()
	                  .feed("ph_corr2d",toTensor2DFloat(img_corr2d, tensorsToClose))
	                  .feed("ph_target_disparity",toTensor2DFloat(img_target, tensorsToClose))
	                  .feed("ph_ntile",toTensor1DInt(img_ntile, tensorsToClose))
	                  .fetch("Disparity_net/stage1done:0")
	                  .run()
	                  .get(0);
        	
        	System.out.println("Stage 1 DONE");
        	
        	System.out.println("Stage 2");
        	
        	// stage 2
        	final Tensor<?> result = bundle.session().runner()
	                  .feed("ph_ntile_out",toTensor1DInt(img_ntile, tensorsToClose))
	                  .fetch("Disparity_net/stage2_out_sparse:0")
	                  .run()
	                  .get(0);
        	
        	System.out.println("Stage 2 DONE: "+result.shape());
        	
        	tensorsToClose.add(result);
        	
        	System.out.println("Copy result to variable");
        	float [][] resultValues = (float[][]) result.copyTo(new float[78408][1]);
        	
        	System.out.println("DONE");
        	
        //} catch (final IllegalStateException ise) {
        //	System.out.println("Very Bad Error (VBE): "+ise);
        //	closeTensors(tensorsToClose);
        } catch (final NumberFormatException nfe) {
			//just skip unparsable lines ?!
		} finally {
			closeTensors(tensorsToClose);
		}
        
        //try (){
            
            //smpb = b.graph();
           
            //Session sess = b.session();
            //System.out.println(b.metaGraphDef());
            
            //final List<String> labels = tensorFlowService.loadLabels(source,
    		//		MODEL_NAME, "imagenet_comp_graph_label_strings.txt");
            //System.out.println("Loaded graph and " + labels.size() + " labels");
            
            //output = sess.runner().feed(o, t).fetch().run().get(0).copyTo()
            
            /*
            try (
            	final Session s = new Session(g);
        			@SuppressWarnings("unchecked")
        			final Tensor<Float> result = (Tensor<Float>) s.runner().feed("input", image)//
        				.fetch("output").run().get(0)
        		){
            	...
            }
            */
            
        //}

        try (Graph g = new Graph()) {
            final String value = "Hello from " + TensorFlow.version();

            // Construct the computation graph with a single operation, a constant
            // named "MyConst" with a value "value".
            try (Tensor t = Tensor.create(value.getBytes("UTF-8"))) {
              // The Java API doesn't yet include convenience functions for adding operations.
              g.opBuilder("Const", "MyConst").setAttr("dtype", t.dataType()).setAttr("value", t).build();
            }

            // Execute the "MyConst" operation in a Session.
            try (Session s = new Session(g);
                // Generally, there may be multiple output tensors, all of them must be closed to prevent resource leaks.
                Tensor output = s.runner().fetch("MyConst").run().get(0)) {
                    System.out.println(new String(output.bytesValue(), "UTF-8"));
                    s.close();
                    output.close();
            }
        }
    }

}


















