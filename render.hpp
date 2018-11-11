#ifndef RENDER_HPP
#define RENDER_HPP

#include "common.hpp"
#include "render.h"
#include <iostream>


using namespace std;

istream& operator>>(istream& is,render& r)
{
    is>>r.type;
    is>>r.x;
    is>>r.y;
    is>>r.z;

    return is;
}

ostream& operator<<(ostream& os,const render& r)
{
    os<<r.type<<" "<<r.x<<" "<<r.y<<" "<<r.z;
    return os;
}


#endif