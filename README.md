This is a set of Elphel plugins for ImageJ 1.x as a Maven project with specified
dependencies.
If cloned and imported as existent Maven project into Eclipse IDE, it allows to launch
ImageJ with plugins in both run and debug modes.

## Java TensorFlow CUDA testing

- As of 2018/09/11 TF for GPU on Maven supports CUDA 9.0 (vs latest 9.2).
Workaround:

```
 ~$ sudo dpkg -i cuda-repo-ubuntu1604-9-0-local_9.0.176-1_amd64.deb
 ~$ sudo apt install cuda-9.0
 # install or link back to cuda-9.2:
 ~$ sudo rm /usr/local/cuda; sudo ln -sf /usr/local/cuda-9.2 /usr/local/cuda
 Then in Eclipse's Run configurations... add an environment variable:
 LD_LIBRARY_PATH = /usr/local/cuda-9.0/lib64
```
