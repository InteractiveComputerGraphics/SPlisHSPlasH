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

#include "../Partio.h"
#include "PartioEndian.h"
#include "../core/ParticleHeaders.h"
#include "ZIP.h"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include <string.h>
#include <algorithm>

namespace Partio
{

using namespace std;

void writeHoudiniStr(ostream& ostream,const string& s)
{
    write<BIGEND>(ostream,(short)s.size());
    ostream.write(s.c_str(),s.size());
}

template<class T>
struct Helper
{
    T addAttribute(ParticlesDataMutable* simple, const char* name, ParticleAttributeType type, int size);
    int registerIndexedStr(const T& attribute,const char* str);
};
template<>
struct Helper<ParticleAttribute>
{
    ParticleAttribute addAttribute(ParticlesDataMutable* simple, const char* name, ParticleAttributeType type, int size) {return simple->addAttribute(name,type,size);}
    int registerIndexedStr(ParticlesDataMutable* simple, const ParticleAttribute& attribute,const char* str) {return simple->registerIndexedStr(attribute,str);}
};
template<>
struct Helper<FixedAttribute>
{
    FixedAttribute addAttribute(ParticlesDataMutable* simple, const char* name, ParticleAttributeType type, int size) {return simple->addFixedAttribute(name,type,size);}
    int registerIndexedStr(ParticlesDataMutable* simple, const FixedAttribute& attribute,const char* str) {return simple->registerFixedIndexedStr(attribute,str);}
};

class DummyAttribute {};
template<>
struct Helper<DummyAttribute>
{
    DummyAttribute addAttribute(ParticlesDataMutable* simple, const char* name, ParticleAttributeType type, int size) { return DummyAttribute(); }
    int registerIndexedStr(ParticlesDataMutable* simple, const DummyAttribute& attribute,const char* str) { return 0; }
};
struct DummyAccessor{
    template<class T>
    DummyAccessor(const T& /*attr*/){}
};

template<class TAttribute, class TAccessor>
bool getAttributes(int& particleSize, vector<int>& attrOffsets, vector<TAttribute>& attrHandles, vector<TAccessor>& accessors, int nAttrib, istream* input, ParticlesDataMutable* simple, bool headersOnly, std::ostream* errorStream)
{
    Helper<TAttribute> helper;
    for(int i=0;i<nAttrib;i++){
        unsigned short nameLength;
        read<BIGEND>(*input,nameLength);
        char* name=new char[nameLength+1];
        input->read(name,nameLength);name[nameLength]=0;
        unsigned short size;
        int houdiniType;
        read<BIGEND>(*input,size,houdiniType);
        if(houdiniType==0 || houdiniType==1 || houdiniType==5){
            // read default values. don't do anything with them
            for(int i=0;i<size;i++) {
                int defaultValue;
                input->read((char*)&defaultValue,sizeof(int));
            }
            ParticleAttributeType type=NONE;
            if(houdiniType==0) type=FLOAT;
            else if(houdiniType==1) type=INT;
            else if(houdiniType==5) type=VECTOR;
            attrHandles.push_back(helper.addAttribute(simple,name,type,size));
            accessors.push_back(TAccessor(attrHandles.back()));
            attrOffsets.push_back(particleSize);
            particleSize+=size;
        }else if(houdiniType==4){
            TAttribute attribute=helper.addAttribute(simple,name,INDEXEDSTR,size);
            attrHandles.push_back(attribute);
            accessors.push_back(TAccessor(attrHandles.back()));
            attrOffsets.push_back(particleSize);
            int numIndices=0;
            read<BIGEND>(*input,numIndices);
            for(int ii=0;ii<numIndices;ii++){
                unsigned short indexNameLength;
                read<BIGEND>(*input,indexNameLength);
                char* indexName=new char[indexNameLength+1];;
                input->read(indexName,indexNameLength);
                indexName[indexNameLength]=0;
                if (!headersOnly) {
                    int id=helper.registerIndexedStr(simple,attribute,indexName);
                    if(id != ii){
                        if(errorStream) *errorStream <<"Partio: error on read, expected registerIndexStr to return index "<<ii<<" but got "<<id<<" for string "<<indexName<<std::endl;
                    }
                }
                delete [] indexName;
            }
            particleSize+=size;
        }else if(houdiniType==2){
            if(errorStream) *errorStream <<"Partio: found attr of type 'string', aborting"<<endl;
            delete [] name;
            simple->release();
            return 0;
        }else{
            if(errorStream) *errorStream <<"Partio: unknown attribute "<<houdiniType<<" type... aborting"<<endl;
            delete [] name;
            simple->release();
            return 0;
        }
        delete[] name;
    }

    return true;
}

// read buffer, seekg doesn't work with gzip
void skip(istream *input, size_t numChars)
{
    static const size_t bufferSize = 4096;
    static char buffer[bufferSize];
    while (numChars>0) {
        int toRead=std::min(numChars,bufferSize);
        input->read(buffer,toRead);
        numChars-=toRead;
    }
}

// ignore primitive attributes, only know about Particle Systems currently
bool skipPrimitives(int nPoints, int nPrims, int nPrimAttrib, istream* input,std::ostream* errorStream)
{
    int particleSize=0;
    vector<int> primAttrOffsets; // offsets in # of 32 bit offsets
    vector<DummyAttribute> primAttrHandles;
    vector<DummyAccessor> primAccessors;
    getAttributes(particleSize, primAttrOffsets, primAttrHandles, primAccessors, nPrimAttrib, input, 0, true, errorStream);

    for(int i=0;i<nPrims;i++) {
        int primType;
        read<BIGEND>(*input,primType);
        if(primType==0x00008000) {
            int size;
            read<BIGEND>(*input,size);
            if(nPoints>=(int)1<<16)
                skip(input,size*sizeof(int));
            else
                skip(input,size*sizeof(unsigned short));
            skip(input,particleSize*sizeof(int));
        } else {
            if(errorStream) *errorStream << "Partio: Unrecognized Primitive Type: 0x" << std::hex << primType << " - Cannot process detail attributes" << std::endl;
            return false;
        }
    }
    return true;
}

ParticlesDataMutable* readBGEO(const char* filename,const bool headersOnly,std::ostream* errorStream)
{
    unique_ptr<istream> input(Gzip_In(filename,ios::in|ios::binary));
    if(!*input){
        if(errorStream) *errorStream<<"Partio: Unable to open file "<<filename<<endl;
        return 0;
    }

    // header values
    char magic[5];
    magic[4]=0;
    char versionChar;
    int version;
    int nPoints;
    int nPrims;
    int nPointGroups;
    int nPrimGroups;
    int nPointAttrib;
    int nVertexAttrib;
    int nPrimAttrib;
    int nAttrib;
    read<BIGEND>(*input,magic[0],magic[1],magic[2],magic[3]);
    read<BIGEND>(*input,versionChar,version,nPoints,nPrims,nPointGroups);
    read<BIGEND>(*input,nPrimGroups,nPointAttrib,nVertexAttrib,nPrimAttrib,nAttrib);


    // Check header magic and version
    const char bgeo_magic[5]={'B','g','e','o',0};
    if(strcmp(magic,bgeo_magic)){
        const char new_bgeo_magic[5]={0x7f,0x4e,0x53,0x4a,0};
        if(!strcmp(magic,new_bgeo_magic)){
            if(errorStream) *errorStream<<"Partio: Attempting to read new BGEO format, we only support old BGEO format. Try writing .bhclassic from Houdini."<<std::endl;
        }else{
            if(errorStream) *errorStream<<"Partio: Magic number '"<<magic<<" of '"<<filename<<"' doesn't match bgeo magic '"<<bgeo_magic<<endl;
        }
        return 0;
    }
    if(version!=5){
        if(errorStream) *errorStream<<"Partio: BGEO must be version 5"<<endl;
        return 0;
    }

    // Allocate a simple particle with the appropriate number of points
    ParticlesDataMutable* simple=0;
    if(headersOnly) simple=new ParticleHeaders;
    else simple=create();

    simple->addParticles(nPoints);


    // Read attribute definitions
    int particleSize=4; // Size in # of 32 bit primitives 
    vector<int> attrOffsets; // offsets in # of 32 bit offsets
    vector<ParticleAttribute> attrHandles;
    vector<ParticleAccessor> accessors;
    attrOffsets.push_back(0); // pull values from byte offset
    attrHandles.push_back(simple->addAttribute("position",VECTOR,3)); // we always have one
    accessors.push_back(ParticleAccessor(attrHandles[0]));
    getAttributes(particleSize, attrOffsets, attrHandles, accessors, nPointAttrib, input.get(), simple, headersOnly, errorStream);

    if(headersOnly) {
        skip(input.get(),nPoints*particleSize*sizeof(int));
    } else {
        // Read the points
        int *buffer=new int[particleSize];

        // make iterator and register accessors
        ParticlesDataMutable::iterator iterator=simple->begin();
        for(size_t i=0;i<accessors.size();i++) iterator.addAccessor(accessors[i]);

        for(ParticlesDataMutable::iterator end=simple->end();iterator!=end;++iterator){
            input->read((char*)buffer,particleSize*sizeof(int));
            for(unsigned int attrIndex=0;attrIndex<attrHandles.size();attrIndex++){
                ParticleAttribute& handle=attrHandles[attrIndex];
                ParticleAccessor& accessor=accessors[attrIndex];
                // TODO: this violates strict aliasing, we could just go to char* and make
                // a different endian swapper
                int* data=accessor.raw<int>(iterator);
                for(int k=0;k<handle.count;k++){
                    BIGEND::swap(buffer[attrOffsets[attrIndex]+k]);
                    data[k]=buffer[attrOffsets[attrIndex]+k];
                }
            }
        }
        delete [] buffer;
    }

    if (!skipPrimitives(nPoints, nPrims, nPrimAttrib, input.get(),errorStream)) return simple;

    particleSize=0;
    vector<int> fixedAttrOffsets; // offsets in # of 32 bit offsets
    vector<FixedAttribute> fixedAttrHandles;
    vector<DummyAccessor> fixedAccessors;
    getAttributes(particleSize, fixedAttrOffsets, fixedAttrHandles, fixedAccessors, nAttrib, input.get(), simple, headersOnly, errorStream);

    if (headersOnly) return simple;

    // Read the points
    int *fixedBuffer=new int[particleSize];
    input->read((char*)fixedBuffer,particleSize*sizeof(int));
    for(unsigned int attrIndex=0;attrIndex<fixedAttrHandles.size();attrIndex++){
        FixedAttribute& handle=fixedAttrHandles[attrIndex];
        // TODO: this violates strict aliasing, we could just go to char* and make
        // a different endian swapper
        for(int k=0;k<handle.count;k++){
            BIGEND::swap(fixedBuffer[fixedAttrOffsets[attrIndex]+k]);
            simple->fixedDataWrite<int>(fixedAttrHandles[attrIndex])[k]=fixedBuffer[fixedAttrOffsets[attrIndex]+k];
        }
    }
    delete [] fixedBuffer;

    // return the populated simpleParticle
    return simple;
}

bool writeBGEO(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream)
{
    unique_ptr<ostream> output(
        compressed ? 
        Gzip_Out(filename,ios::out|ios::binary)
        :new ofstream(filename,ios::out|ios::binary));

    if(!*output){
        if(errorStream) *errorStream <<"Partio Unable to open file "<<filename<<endl;
        return false;
    }

    int magic=((((('B'<<8)|'g')<<8)|'e')<<8)|'o';
    char versionChar='V';
    int version=5;
    int nPoints=p.numParticles();
    int nPrims=0;
    int nPointGroups=0;
    int nPrimGroups=0;
    int nPointAttrib=p.numAttributes()-1;
    int nVertexAttrib=0;
    int nPrimAttrib=0;
    int nAttrib=p.numFixedAttributes();

    write<BIGEND>(*output,magic,versionChar,version,nPoints,nPrims,nPointGroups);
    write<BIGEND>(*output,nPrimGroups,nPointAttrib,nVertexAttrib,nPrimAttrib,nAttrib);

    vector<ParticleAttribute> handles;
    vector<ParticleAccessor> accessors;
    vector<int> attrOffsets;
    bool foundPosition=false;
    int particleSize=4;
    for(int i=0;i<p.numAttributes();i++){
        ParticleAttribute attr;
        p.attributeInfo(i,attr);
        if(attr.name=="position"){
            attrOffsets.push_back(0);
            foundPosition=true;
        }else{
            writeHoudiniStr(*output,attr.name);
            if(attr.type==INDEXEDSTR){
                int houdiniType=4;
                unsigned short size=attr.count;
                const std::vector<std::string>& indexTable=p.indexedStrs(attr);
                int numIndexes=indexTable.size();
                write<BIGEND>(*output,size,houdiniType,numIndexes);
                for(int i=0;i<numIndexes;i++)
                    writeHoudiniStr(*output,indexTable[i]);
            }else{
                int houdiniType=0;
                switch(attr.type){
                    case FLOAT: houdiniType=0;break;
                    case INT: houdiniType=1;break;
                    case VECTOR: houdiniType=5;break;
                    case INDEXEDSTR:
                    case NONE: assert(false);houdiniType=0;break;
                }
                unsigned short size=attr.count;
                write<BIGEND>(*output,size,houdiniType);
                for(int i=0;i<attr.count;i++){
                    int defaultValue=0;
                    write<BIGEND>(*output,defaultValue);
                }
            }
            attrOffsets.push_back(particleSize);
            particleSize+=attr.count;
        }
        handles.push_back(attr);
        accessors.push_back(ParticleAccessor(handles.back()));
    }
    if(!foundPosition){
        if(errorStream) *errorStream <<"Partio: didn't find attr 'position' while trying to write GEO"<<endl;
        return false;
    }

    ParticlesData::const_iterator iterator=p.begin();
    for(size_t i=0;i<accessors.size();i++) iterator.addAccessor(accessors[i]);

    int *buffer=new int[particleSize];
    for(ParticlesData::const_iterator end=p.end();iterator!=end;++iterator){
        for(unsigned int attrIndex=0;attrIndex<handles.size();attrIndex++){
            ParticleAttribute& handle=handles[attrIndex];
            ParticleAccessor& accessor=accessors[attrIndex];
            // TODO: this violates strict aliasing, we could just go to char* and make
            // a different endian swapper
            const int* data=accessor.raw<int>(iterator);
            for(int k=0;k<handle.count;k++){
                buffer[attrOffsets[attrIndex]+k]=data[k];
                BIGEND::swap(buffer[attrOffsets[attrIndex]+k]);
            }
        }
        // set homogeneous coordinate
        float *w=(float*)&buffer[3];
        *w=1.;
        BIGEND::swap(*w);
        output->write((char*)buffer,particleSize*sizeof(int));
    }
    delete [] buffer;

    vector<FixedAttribute> fixedHandles;
    vector<int> fixedAttrOffsets;
    particleSize=0;
    for(int i=0;i<p.numFixedAttributes();i++){
        FixedAttribute attr;
        p.fixedAttributeInfo(i,attr);

        writeHoudiniStr(*output,attr.name);
        if(attr.type==INDEXEDSTR){
            int houdiniType=4;
            unsigned short size=attr.count;
            const std::vector<std::string>& indexTable=p.fixedIndexedStrs(attr);
            int numIndexes=indexTable.size();
            write<BIGEND>(*output,size,houdiniType,numIndexes);
            for(int ii=0;ii<numIndexes;ii++) {
                writeHoudiniStr(*output,indexTable[ii]);
            }
        }else{
            int houdiniType=0;
            switch(attr.type){
            case FLOAT: houdiniType=0;break;
            case INT: houdiniType=1;break;
            case VECTOR: houdiniType=5;break;
            case INDEXEDSTR:
            case NONE: assert(false);houdiniType=0;break;
            }
            unsigned short size=attr.count;
            write<BIGEND>(*output,size,houdiniType);
            for(int i=0;i<attr.count;i++){
                int defaultValue=0;
                write<BIGEND>(*output,defaultValue);
            }
        }
        fixedAttrOffsets.push_back(particleSize);
        particleSize+=attr.count;

        fixedHandles.push_back(attr);
    }

    int *fixedBuffer=new int[particleSize];

    for(unsigned int attrIndex=0;attrIndex<fixedHandles.size();attrIndex++){
        FixedAttribute& handle=fixedHandles[attrIndex];
        // TODO: this violates strict aliasing, we could just go to char* and make
        // a different endian swapper

        for(int k=0;k<handle.count;k++){

            fixedBuffer[fixedAttrOffsets[attrIndex]+k]=p.fixedData<int>(fixedHandles[attrIndex])[k];
            BIGEND::swap(fixedBuffer[fixedAttrOffsets[attrIndex]+k]);
        }
    }
    output->write((char*)fixedBuffer,particleSize*sizeof(int));

    delete [] fixedBuffer;

    // Write extra
    write<BIGEND>(*output,(char)0x00);
    write<BIGEND>(*output,(char)0xff);

    // success
    return true;
}

}
