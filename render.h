#ifndef Render_H
#define Render_H
#include "common.hpp"
#include <iostream>

using namespace std;
class render
{
public:
    int type;
    float x;
    float y;
    float z;
    render() {}
};

    istream& operator>>(istream& is,render& r);
    ostream& operator<<(ostream& is,const render& r);
#endif