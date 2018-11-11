all:
	g++ main.cpp -std=c++11 `pkg-config --cflags --libs eigen3 glfw3 gl glu`
	./a.out

clean:
	rm a.out

render:
	g++ render.cpp -o render.out -std=c++11 `pkg-config --cflags --libs eigen3 glfw3 gl glu`
