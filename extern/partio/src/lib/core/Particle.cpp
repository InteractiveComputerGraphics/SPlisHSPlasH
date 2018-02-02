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
#ifdef PARTIO_WIN32
#    define NOMINMAX
#endif
#define __STDC_LIMIT_MACROS
#include "../PartioVec3.h"
#include "ParticleSimple.h"
#include "ParticleSimpleInterleave.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cassert>
#include <algorithm>
namespace Partio{

std::string
TypeName(ParticleAttributeType attrType)
{
    switch(attrType){
        case NONE: return "NONE";
        case VECTOR: return "VECTOR";
        case FLOAT: return "FLOAT";
        case INT: return "INT";
        case INDEXEDSTR: return "INDEXEDSTR";
        default: return 0;
    }
}

ParticlesDataMutable*
create()
{
   return new ParticlesSimple;
}

ParticlesDataMutable*
createInterleave()
{
    return new ParticlesSimpleInterleave;
}


ParticlesDataMutable*
cloneSchema(const ParticlesData& other)
{
    ParticlesDataMutable* p = create();

    FixedAttribute detail;
    for(int i=0;i<other.numFixedAttributes();++i) {
        other.fixedAttributeInfo(i,detail);
        p->addFixedAttribute(detail.name.c_str(), detail.type, detail.count);
    }

    ParticleAttribute attr;
    for(int j=0;j<other.numAttributes();++j) {
        other.attributeInfo(j,attr);
        p->addAttribute(attr.name.c_str(), attr.type, attr.count);
    }

    return p;
}

ParticlesDataMutable*
clone(const ParticlesData& other, bool particles)
{
    ParticlesDataMutable* p = create();
    // Fixed attributes
    FixedAttribute srcFixedAttr, dstFixedAttr;
    for (int i(0), iend(other.numFixedAttributes()); i < iend; ++i) {
        other.fixedAttributeInfo(i, srcFixedAttr);

        dstFixedAttr = p->addFixedAttribute(
            srcFixedAttr.name.c_str(), srcFixedAttr.type, srcFixedAttr.count);
        assert(srcFixedAttr.type == dstFixedAttr.type);
        assert(srcFixedAttr.count == dstFixedAttr.count);
        // Register indexed strings
        if (srcFixedAttr.type == Partio::INDEXEDSTR) {
            const std::vector<std::string>& values = other.fixedIndexedStrs(srcFixedAttr);
            for (int j = 0, jend = values.size(); j < jend; ++j) {
                p->registerFixedIndexedStr(dstFixedAttr, values[j].c_str());
            }
        }
        // Copy fixed data
        const void* src = other.fixedData<void>(srcFixedAttr);
        void* dst = p->fixedDataWrite<void>(dstFixedAttr);
        size_t size = Partio::TypeSize(dstFixedAttr.type) * dstFixedAttr.count;
        std::memcpy(dst, src, size);
    }
    if (!particles) {
        return p;
    }
    // Particle data
    Partio::ParticleAttribute srcAttr, dstAttr;
    const int numAttributes = other.numAttributes();
    const size_t numParticles = other.numParticles();
    std::vector<Partio::ParticleAttribute> dstAttrs;

    p->addParticles(numParticles);

    // We can't assume that the particle backend stores data contiguously, so
    // we copy one particle at a time.  A bulk memcpy would be faster.
    for (int i = 0; i < numAttributes; ++i) {
        other.attributeInfo(i, srcAttr);
        // Register indexed strings
        if (srcAttr.type == Partio::INDEXEDSTR) {
            const std::vector<std::string>& values = other.indexedStrs(srcAttr);
            for (int m = 0, mend = values.size(); m < mend; ++m) {
                p->registerIndexedStr(dstAttr, values[m].c_str());
            }
        }
        size_t size = Partio::TypeSize(srcAttr.type) * srcAttr.count;
        dstAttr = p->addAttribute(srcAttr.name.c_str(), srcAttr.type, srcAttr.count);

        for (Partio::ParticleIndex j = 0; j < numParticles; ++j) {
            const void *src = other.data<void>(srcAttr, j);
            void *dst = p->dataWrite<void>(dstAttr, j);
            std::memcpy(dst, src, size);
        }
    }

    return p;
}


template<ParticleAttributeType ETYPE> void
printAttr(const ParticlesData* p,const ParticleAttribute& attr,const int particleIndex)
{
    typedef typename ETYPE_TO_TYPE<ETYPE>::TYPE TYPE;
    const TYPE* data=p->data<TYPE>(attr,particleIndex);
    for(int k=0;k<attr.count;k++) std::cout<<" "<<data[k];
}

void
print(const ParticlesData* particles)
{
    std::cout<<"Particle count "<<particles->numParticles()<<std::endl;
    std::cout<<"Attribute count "<<particles->numAttributes()<<std::endl;

    std::vector<ParticleAttribute> attrs;
    for(int i=0;i<particles->numAttributes();i++){
        ParticleAttribute attr;
        particles->attributeInfo(i,attr);
        attrs.push_back(attr);
        std::cout<<"attribute "<<attr.name<<" "<<int(attr.type)<<" "<<attr.count<<std::endl;
    }

    int numToPrint=std::min(10,particles->numParticles());
    std::cout<<"num to print "<<numToPrint<<std::endl;

    ParticlesData::const_iterator it=particles->begin(),end=particles->end();
    std::vector<ParticleAccessor> accessors;
    for(size_t k=0;k<attrs.size();k++) accessors.push_back(ParticleAccessor(attrs[k]));
    for(size_t k=0;k<attrs.size();k++) it.addAccessor(accessors[k]);

    for(int i=0;i<numToPrint && it != end;i++){
        std::cout<<i<<": ";
        for(unsigned int k=0;k<attrs.size();k++){
            switch(attrs[k].type){
            case NONE:break;
            case FLOAT:
            case VECTOR:
                for(int c=0;c<attrs[k].count;c++) std::cout<<accessors[k].raw<float>(it)[c];
                break;
            case INT:
                for(int c=0;c<attrs[k].count;c++) std::cout<<accessors[k].raw<int>(it)[c];
                break;
            case INDEXEDSTR:
                for(int c=0;c<attrs[k].count;c++) std::cout<<accessors[k].raw<int>(it)[c];
                break;
            }
        }
        std::cout<<std::endl;
    }
}

double hash(int n, double* args)
{
    // combine args into a single seed
    uint32_t seed = 0;
    for (int i = 0; i < n; i++) {
        // make irrational to generate fraction and combine xor into 32 bits
        int exp=0;
        double frac = std::frexp(args[i] * double(M_E*M_PI), &exp);
        uint32_t s = (uint32_t) (frac * UINT32_MAX) ^ (uint32_t) exp;

        // blend with seed (constants from Numerical Recipes, attrib. from Knuth)
        static const uint32_t M = 1664525, C = 1013904223;
        seed = seed * M + s + C;
    }

    // tempering (from Matsumoto)
    seed ^= (seed >> 11);
    seed ^= (seed << 7) & 0x9d2c5680UL;
    seed ^= (seed << 15) & 0xefc60000UL;
    seed ^= (seed >> 18);

    // permute
    static unsigned char p[256] = {
        148,201,203,34,85,225,163,200,174,137,51,24,19,252,107,173,
        110,251,149,69,180,152,141,132,22,20,147,219,37,46,154,114,
        59,49,155,161,239,77,47,10,70,227,53,235,30,188,143,73,
        88,193,214,194,18,120,176,36,212,84,211,142,167,57,153,71,
        159,151,126,115,229,124,172,101,79,183,32,38,68,11,67,109,
        221,3,4,61,122,94,72,117,12,240,199,76,118,5,48,197,
        128,62,119,89,14,45,226,195,80,50,40,192,60,65,166,106,
        90,215,213,232,250,207,104,52,182,29,157,103,242,97,111,17,
        8,175,254,108,208,224,191,112,105,187,43,56,185,243,196,156,
        246,249,184,7,135,6,158,82,130,234,206,255,160,236,171,230,
        42,98,54,74,209,205,33,177,15,138,178,44,116,96,140,253,
        233,125,21,133,136,86,245,58,23,1,75,165,92,217,39,0,
        218,91,179,55,238,170,134,83,25,189,216,100,129,150,241,210,
        123,99,2,164,16,220,121,139,168,64,190,9,31,228,95,247,
        244,81,102,145,204,146,26,87,113,198,181,127,237,169,28,93,
        27,41,231,248,78,162,13,186,63,66,131,202,35,144,222,223};
    union {
        uint32_t i;
        unsigned char c[4];
    } u1, u2;
    u1.i = seed;
    u2.c[3] = p[u1.c[0]];
    u2.c[2] = p[(u1.c[1]+u2.c[3])&0xff];
    u2.c[1] = p[(u1.c[2]+u2.c[2])&0xff];
    u2.c[0] = p[(u1.c[3]+u2.c[1])&0xff];

    // scale to [0.0 .. 1.0]
    return u2.i * (1.0/UINT32_MAX);
}

struct IdAndIndex {
    IdAndIndex(int id, int index) : _id(id), _index(index) {}
    int _id, _index;
    bool operator<(const IdAndIndex &other) const {
        return _id < other._id;
    }
};

template<class T> T smoothstep(T t){return (3.-2.*t)*t*t;}

void addClusterAttribute(ParticlesDataMutable* cluster, ParticleAttribute& clusterAttribute, const ParticlesDataMutable* particle, const int index, const ParticleAttribute& attribute, const int neighborIndex, const std::vector<std::pair<ParticleIndex,float> >& indexAndInterp)
{
    switch(attribute.type){
    case Partio::VECTOR:
    case Partio::FLOAT:
        {
            float* data = particle->dataWrite<float>(attribute,index);
            float* neighborData = particle->dataWrite<float>(attribute,neighborIndex);
            for (size_t i=0; i<indexAndInterp.size(); i++) {
                float* clusterData = cluster->dataWrite<float>(clusterAttribute,indexAndInterp[i].first);
                for (int j=0; j<attribute.count; j++) {
                    clusterData[j] = indexAndInterp[i].second ? data[j]+indexAndInterp[i].second*(neighborData[j]-data[j]) : data[j];
                }
            }
            break;
        }
    case Partio::INT:
    case Partio::INDEXEDSTR:
        {
            int* data = particle->dataWrite<int>(attribute,index);
            for (size_t i=0; i<indexAndInterp.size(); i++) {
                int* clusterData = cluster->dataWrite<int>(clusterAttribute,indexAndInterp[i].first);
                for (int j=0; j<attribute.count; j++) {
                    clusterData[j] = data[j];
                }
            }
            break;
        }
    }
}

ParticlesDataMutable*
computeClustering(ParticlesDataMutable* particles, const int numNeighbors,const double radiusSearch,const double radiusInside,const int connections,const double density)
{
    ParticleAttribute posAttr;
    bool hasPosAttr = false;
    ParticleAttribute idAttr;
    bool hasIdAttr = false;
    ParticlesDataMutable* cluster = create();
    std::vector<ParticleAttribute> attributes;
    std::vector<ParticleAttribute> clusterAttributes;
    for (int i=0; i<particles->numAttributes(); i++) {
        ParticleAttribute attr;
        if (particles->attributeInfo(i,attr)) {
            if (attr.type == Partio::NONE) continue;
            attributes.push_back(attr);
            clusterAttributes.push_back(cluster->addAttribute(attributes[i].name.c_str(),attributes[i].type,attributes[i].count));
            if (attr.type == Partio::INDEXEDSTR) {
                const std::vector<std::string>& strings = particles->indexedStrs(attr);
                for (size_t j=0; j<strings.size(); j++) {
                    cluster->registerIndexedStr(clusterAttributes.back(),strings[j].c_str());
                }
            }
            if (attr.name == "position") {
                posAttr = attr;
                hasPosAttr = true;
            } else if (attr.name == "id") {
                idAttr = attr;
                hasIdAttr = true;
            }
        }
    }
    if (!hasPosAttr || posAttr.type != VECTOR && posAttr.type != FLOAT || posAttr.count !=3) {
        cluster->release();
        return 0;
    }
    if (!hasIdAttr ||idAttr.type != INT || idAttr.count != 1) {
        cluster->release();
        return 0;
    }
    particles->sort();
    ParticleAttribute clusterIdAttr = cluster->addAttribute("clusterId", Partio::INT, 1);
    for (int index=0; index<particles->numParticles(); index++) {
        const float* center=particles->data<float>(posAttr,index);
        Vec3 position(center[0], center[1], center[2]);
        int id = particles->data<int>(idAttr,index)[0];
        double radius = std::max(0., radiusInside);
        radius = std::min(100., radiusInside);
        double innerRadius = .01 * radius * radiusSearch;
        double invRadius = 1 / (radiusSearch - innerRadius);

        std::vector<Partio::ParticleIndex> points;
        std::vector<float> distSq;
        particles->findNPoints(center, numNeighbors, radiusSearch, points, distSq);

        std::vector<IdAndIndex> idAndIndex;
        idAndIndex.reserve(points.size());
        idAndIndex.push_back(IdAndIndex(id, 0));
        for (unsigned int i = 0; i < points.size(); i++) {
            const int pointid = particles->data<int>(idAttr,points[i])[0];
            if (pointid != id) idAndIndex.push_back(IdAndIndex(pointid, points[i]));
        }
        std::sort(++idAndIndex.begin(), idAndIndex.end());

        std::vector<std::pair<ParticleIndex,float> > originalPoint;
        originalPoint.push_back(std::make_pair(cluster->addParticle(),0));
        int clusterId = 0;
        cluster->dataWrite<int>(clusterIdAttr,originalPoint.back().first)[0] = clusterId++;
        for (size_t j = 0; j < attributes.size(); j++) {
            addClusterAttribute(cluster, clusterAttributes[j], particles, index, attributes[j], 0, originalPoint);
        }

        double hashArgs[3];
        hashArgs[0] = id;

        int foundConnections = std::min((int) idAndIndex.size(), connections+1);
        for (int i = 1; i < foundConnections; i++) {
            const float* neighbor=particles->data<float>(posAttr,idAndIndex[i]._index);
            Vec3 neighborPosition(neighbor[0], neighbor[1], neighbor[2]);
            Vec3 dir = neighborPosition - position;

            // calculate number of instances based on density
            int numInstances = 0;
            double len = dir.length();
            if (len < innerRadius) {
                numInstances = density * len;
            } else {
                numInstances = density * len * smoothstep(1.-(len-innerRadius)*invRadius);
            }
            std::vector<std::pair<ParticleIndex,float> > indexAndInterp;
            for (int j = 0; j < numInstances; j++) {
                hashArgs[1] = idAndIndex[i]._id;
                hashArgs[2] = j;
                indexAndInterp.push_back(std::make_pair(cluster->addParticle(),hash(3,hashArgs)));
                cluster->dataWrite<int>(clusterIdAttr,indexAndInterp.back().first)[0] = clusterId++;
            }
            for (size_t j = 0; j < attributes.size(); j++) {
                addClusterAttribute(cluster, clusterAttributes[j], particles, index, attributes[j], idAndIndex[i]._index, indexAndInterp);
            }
        }
    }
    return cluster;;
}

}
