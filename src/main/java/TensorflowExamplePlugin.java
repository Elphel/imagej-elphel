/**
 * Copyright (C) 2018 Elphel, Inc.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

import org.tensorflow.Graph;
import org.tensorflow.Session;
import org.tensorflow.Tensor;
import org.tensorflow.TensorFlow;

public class TensorflowExamplePlugin
{	
    public void run() 
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
          }
        }
      }
    
	
}
