This is a basic C++ program that uses "old-style" OpenGL to render a fluid. It uses Eigen (http://eigen.tuxfamily.org/) to represent vectors, and GLFW (http://www.glfw.org/) for window management and input handling.

Compiling the code is as simple as running the command

	g++ main.cpp -std=c++11 `pkg-config --cflags --libs eigen3 glfw3 gl glu`

Maybe this will work for you! I don't know the exact procedure for compiling in Visual Studio on Windows or in XCode on Mac (or even on the IIT labs' Linux machines), but it will basically involve

  i. providing the include paths for Eigen and GLFW, and

  ii. linking with the libraries for GLFW (usually libglfw), GLU (libGLU), and OpenGL (libGL).

On Mac, instead of adding the OpenGL and GLU libraries individually, you will have to add the OpenGL framework instead (either in XCode, or by adding '-framework OpenGL' to the g++ command-line arguments).
