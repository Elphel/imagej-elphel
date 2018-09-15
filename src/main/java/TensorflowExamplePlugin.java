/**
 * Copyright (C) 2018 Elphel, Inc.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

import org.tensorflow.Graph;
import org.tensorflow.Session;
import org.tensorflow.Tensor;
import org.tensorflow.TensorFlow;
import org.tensorflow.SavedModelBundle;

/**

Note 0: articles & examples:

- https://divis.io/2017/11/enterprise-tensorflow-1/
- https://divis.io/2018/01/enterprise-tensorflow-code-examples/
- https://divis.io/2018/01/enterprise-tensorflow-2-saving-a-trained-model/
- https://divis.io/2018/01/enterprise-tensorflow-3-loading-a-savedmodel-in-java/
- https://divis.io/2018/01/enterprise-tensorflow-4-executing-a-tensorflow-session-in-java/
- https://www.programcreek.com/java-api-examples/?api=org.tensorflow.SavedModelBundle
- https://github.com/imagej/imagej-tensorflow

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
    public final static String PB_TAG = "model_pb";

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

    public static void main() throws Exception {

        final Graph smpb;

        float [][] img_corr2d = new float[78408][324];
        float [][] img_target = new float[78408][  1];
        int     [] img_ntile  = new   int[78408];

        // init ntile
        for(int i=0;i<img_ntile.length;i++){
        	img_ntile[i] = i;
        }

        try (SavedModelBundle b = SavedModelBundle.load(EXPORTDIR,PB_TAG)){
            System.out.println("OK");
            smpb = b.graph();
           
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
            
        }

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


















