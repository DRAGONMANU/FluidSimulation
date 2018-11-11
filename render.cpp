#include "camera.hpp"
#include "draw.hpp"
#include "gui.hpp"
#include "lighting.hpp"
#include "text.hpp"
#include "render.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cmath>
#include <unordered_map>

using namespace std;

Window window;
Camera camera;
Lighting lighting;
Text text;


float dt = 0.001;
float t = 0;
bool paused = false;

/////////////////////////////////////////////////

ifstream in("drop.txt");

void drawStuff() 
{
    
    setColor(vec3(0.2f,0.5f,1));
    int n=0;
    while(n<1024)
    {
        n++;
        render obj;
        in>>obj;
        drawSphere(vec3(obj.x,obj.y,obj.z),0.6);
    }
}

void drawWorld() 
{
    camera.apply(window);
    lighting.apply();
    clear(vec3(0.9,0.9,0.9));
    drawStuff();
    setColor(vec3(0,0,0));
    text.draw("WASD and LShift/LCtrl to move camera", -0.9, 0.90);
    text.draw("Mouse to rotate view", -0.9, 0.85);
    text.draw("Space to play/pause animation", -0.9, 0.80);
    text.draw(to_string(t), -0.9, 0.70);
}


void keyPressed(int key) 
{
    // See http://www.glfw.org/docs/latest/group__keys.html for key codes
    if (key == GLFW_KEY_SPACE)
        paused = !paused;
    if (key == GLFW_KEY_Q)
        exit(0);
}

int main(int argc, char **argv) 
{
    window.create("Animation", 1024, 768);
    window.onKeyPress(keyPressed);
    camera.lookAt(vec3(0,20,30), vec3(0,0,0));
    lighting.createDefault();
    text.initialize();

    while (!window.shouldClose() && t < 5) 
    {
        camera.processInput(window);
        if (!paused)
            t+=dt;
        window.prepareDisplay();
        drawWorld();
        window.updateDisplay();
        window.waitForNextFrame(dt);
    }
}
