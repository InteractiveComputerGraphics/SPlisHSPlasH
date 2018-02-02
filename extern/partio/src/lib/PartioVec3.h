/*
PARTIO SOFTWARE
Copyright 2010 Disney Enterprises, Inc. All rights reserved

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.

* The names "Disney", "Walt Disney Pictures", "Walt Disney Animation
Studios" or the names of its contributors may NOT be used to
endorse or promote products derived from this software without
specific prior written permission from Walt Disney Pictures.

Disclaimer: THIS SOFTWARE IS PROVIDED BY WALT DISNEY PICTURES AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, NONINFRINGEMENT AND TITLE ARE DISCLAIMED.
IN NO EVENT SHALL WALT DISNEY PICTURES, THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND BASED ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
*/

#pragma once

#include <algorithm>
#include <iostream>
#include <cmath>

namespace Partio{

class Vec3
{
public:
    float x,y,z;

    inline Vec3()
        :x(0),y(0),z(0)
    {}

    inline Vec3(const float x,const float y,const float z)
        :x(x),y(y),z(z)
    {}

    inline Vec3(const float v[3])
        :x(v[0]),y(v[1]),z(v[2])
    {}

    inline float length()
    {return std::sqrt(x*x+y*y+z*z);}

    inline float normalize()
    {float l=length();x/=l;y/=l;z/=l;return l;}

    inline Vec3 normalized() const
    {Vec3 foo(x,y,z);foo.normalize();return foo;}

    inline Vec3 operator*(const float a) const
    {return Vec3(a*x,a*y,a*z);}

    inline Vec3 operator-(const Vec3& v) const
    {return Vec3(x-v.x,y-v.y,z-v.z);}

    inline Vec3 operator+(const Vec3& v) const
    {return Vec3(x+v.x,y+v.y,z+v.z);}

    inline Vec3 operator+=(const Vec3& v)
    {x+=v.x;y+=v.y;z+=v.z;return *this;}

    inline Vec3 cross(const Vec3& v) const
    {Vec3 ret(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);return ret;}

    inline Vec3 min(const Vec3& v) const
    {
        return Vec3(std::min(x,v.x),std::min(y,v.y),std::min(z,v.z));
    }

    inline Vec3 max(const Vec3& v) const
    {
        return Vec3(std::max(x,v.x),std::max(y,v.y),std::max(z,v.z));
    }
};


std::ostream& operator<<(std::ostream& stream,const Vec3& v)
{return stream<<v.x<<" "<<v.y<<" "<<v.z;}

Vec3 operator*(const float a,const Vec3& v)
{
    return Vec3(a*v.x,a*v.y,a*v.z);
}

}
