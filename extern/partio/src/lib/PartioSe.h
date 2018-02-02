/**
PARTIO SOFTWARE
Copyright 202 Disney Enterprises, Inc. All rights reserved

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
#include <Partio.h>
#include <SeExpression.h>
#include <map>

namespace Partio{

template<class T>
class AttribVar:public SeExprVarRef
{
    Partio::ParticlesDataMutable* parts;
    Partio::ParticleAttribute attr;
    int& currentIndex;
    int clampedCount;
public:
    AttribVar(Partio::ParticlesDataMutable* parts,
        Partio::ParticleAttribute attr,int& currentIndex)
        :parts(parts),attr(attr),currentIndex(currentIndex),clampedCount(std::min(attr.count,3))
    {}

    bool isVec(){return attr.count!=1;}
    void eval(const SeExprVarNode* node,SeVec3d& result){
        const T* ptr=parts->data<T>(attr,currentIndex);
        //std::cerr<<"in eval for "<<attr.name<<" count is "<<clampedCount<<" cur "<<currentIndex<<std::endl;
        for(int k=0;k<clampedCount;k++){
            result[k]=ptr[k];
        }
        // set any remaining fields (i.e. if clampedCount is 2)
        for(int k=clampedCount;k<3;k++){
            result[k]=0;
        }
    }
};

struct SimpleVar:public SeExprVarRef{
    double val;
    SimpleVar():val(0){}
    bool isVec(){return false;}
    void eval(const SeExprVarNode* node,SeVec3d& result){
        result[0]=result[1]=result[2]=val;
    }
};

/// Class that maps back to the partio data
template<class T> class VarToPartio;

/// NOTE: This class is experimental and may be deleted/modified in future versions
class PartioSe:public SeExpression{
    bool isPaired;
    int currentIndex;
    Partio::ParticleAttribute pairH1,pairH2;
    int pairIndex1,pairIndex2;
    Partio::ParticlesDataMutable* parts;
    Partio::ParticlesDataMutable* partsPairing;
    typedef std::map<std::string,AttribVar<int>*> IntVarMap;
    mutable IntVarMap intVars;
    typedef std::map<std::string,AttribVar<float>*> FloatVarMap;
    mutable FloatVarMap floatVars;

    typedef std::vector<VarToPartio<int>*> IntVarToPartio;
    typedef std::vector<VarToPartio<float>*> FloatVarToPartio;
    IntVarToPartio intVarToPartio;
    FloatVarToPartio floatVarToPartio;

    mutable SimpleVar indexVar,countVar,timeVar;

public:
    typedef  SeExpression::LocalVarTable::const_iterator LocalVarTableIterator;
    
    PartioSe(Partio::ParticlesDataMutable* parts,const char* expr);
    PartioSe(Partio::ParticlesDataMutable* partsPairing,Partio::ParticlesDataMutable* parts,const char* expr);
    void addSet(const char* prefix,Partio::ParticlesDataMutable* parts,int& setIndex);
    void addExport(const std::string& name,LocalVarTableIterator it,Partio::ParticlesDataMutable* parts,int& setIndex);
    virtual ~PartioSe();
    bool runAll();
    bool runRandom();
    void run(int i);
    bool runRange(int istart,int iend);
    void setTime(float val);
    SeExprVarRef*  resolveVar(const std::string& s) const;
private:
    PartioSe(const PartioSe&);
    PartioSe& operator=(const PartioSe&);
};
}
