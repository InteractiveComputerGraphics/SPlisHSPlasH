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
#include "PartioSe.h"
#include <cstdlib>

namespace Partio{

template<class T> class VarToPartio
{
    Partio::ParticlesDataMutable* parts;
    const SeExprLocalVarRef* local;
    Partio::ParticleAttribute attr;
    int& currentIndex;
    int clampedCount;

public:
    VarToPartio(Partio::ParticlesDataMutable* parts,const SeExprLocalVarRef* local,
        Partio::ParticleAttribute attr,int& currentIndex)
        :parts(parts),local(local),attr(attr),currentIndex(currentIndex),
        clampedCount(std::min(3,attr.count))
        
    {}

    void mapBack(){
        T* ptr=parts->dataWrite<T>(attr,currentIndex);
        for(int k=0;k<clampedCount;k++){
            ptr[k]=local->val[k];
        }
        if(clampedCount == 1) for(int k=clampedCount;k<attr.count;k++)
            ptr[k]=ptr[clampedCount-1];
        else
            for(int k=clampedCount;k<attr.count;k++) ptr[k]=0;
    }
};

void PartioSe::addSet(const char* suffix,Partio::ParticlesDataMutable* parts,int& setIndex)
{
    int attrn=parts->numAttributes();
    for(int i=0;i<attrn;i++){
        Partio::ParticleAttribute attr;
        parts->attributeInfo(i,attr);
        if(attr.type == Partio::FLOAT || attr.type == Partio::VECTOR){
            //std::cerr<<"mapping float var "<<attr.name<<std::endl;
            floatVars[attr.name+suffix]=new AttribVar<float>(parts,attr,setIndex);
        }else if(attr.type == Partio::INT || attr.type == Partio::INDEXEDSTR){
            //std::cerr<<"mapping int var "<<attr.name<<std::endl;
            intVars[attr.name+suffix]=new AttribVar<int>(parts,attr,setIndex);
        }else{
            std::cerr<<"failed to map variable "<<attr.name<<std::endl;
        }
    }
}

PartioSe::PartioSe(Partio::ParticlesDataMutable* parts,Partio::ParticlesDataMutable* partsPairing,const char* expression)
:SeExpression(expression),parts(parts),partsPairing(partsPairing),isPaired(true)
{
    parts->attributeInfo("p1",pairH1);
    parts->attributeInfo("p2",pairH2);
    // add variable references
    addSet("_1",partsPairing,pairIndex1);
    addSet("_2",partsPairing,pairIndex2);
    addSet("",parts,currentIndex);
    // checkout what we build
    isValid();
    // go through local variables and look for exports
    const SeExpression::LocalVarTable& vars=getLocalVars();
    typedef  SeExpression::LocalVarTable::const_iterator LocalVarTableIterator;
    for(LocalVarTableIterator it=vars.begin(),itend=vars.end();it != itend;++it){
        //std::cerr<<"assignment of "<<it->first<<std::endl;
        size_t idx=it->first.rfind("_");
        std::string realName=it->first;
        int whichNum;
        if(idx != std::string::npos){
            realName=it->first.substr(0,idx);
            whichNum=atoi(it->first.substr(idx+1).c_str());
            //std::cerr<<"adding as stripped pairing "<<realName<<" which "<<whichNum<<std::endl;
            if(whichNum != 1 && whichNum != 2){
                std::cerr<<"Expect _1 and _2 on the end of particle variables got "<<whichNum<<std::endl;;
            }else{
                int *idxPtr=whichNum == 1 ? &pairIndex1 : &pairIndex2;
                addExport(realName,it,partsPairing,*idxPtr);
            }
        }else{
            //std::cerr<<"adding as normal pairing "<<realName<<std::endl;
            addExport(realName,it,parts,currentIndex);
        }
    }
}

PartioSe::PartioSe(Partio::ParticlesDataMutable* parts,const char* expression)
:SeExpression(expression),parts(parts),partsPairing(0),isPaired(false)
{
    addSet("",parts,currentIndex);
    isValid();
    const SeExpression::LocalVarTable& vars=getLocalVars();
    typedef  SeExpression::LocalVarTable::const_iterator LocalVarTableIterator;
    for(LocalVarTableIterator it=vars.begin(),itend=vars.end();it != itend;++it){
        //std::cerr<<"assignment of "<<it->first<<std::endl;
        addExport(it->first,it,parts,currentIndex);
    }
        
}

void PartioSe::addExport(const std::string& name,LocalVarTableIterator it,Partio::ParticlesDataMutable* parts,int& setIndex){
        Partio::ParticleAttribute attr;
        bool isParticleAttr=parts->attributeInfo(name.c_str(),attr);
        if(isParticleAttr){
            if(attr.type==Partio::FLOAT || attr.type==Partio::VECTOR){
                floatVarToPartio.push_back(new VarToPartio<float>(parts,&it->second,attr,setIndex));
            }else if(attr.type==Partio::INT || attr.type==Partio::INDEXEDSTR){
                intVarToPartio.push_back(new VarToPartio<int>(parts,&it->second,attr,setIndex));
            }else{
                std::cerr<<"unknown particle attribute type "<<attr.name<<std::endl;
            }
        }else{
            //std::cerr<<"not found "<<name<<std::endl;

        }

}

bool PartioSe::runAll()
{
    return runRange(0,parts->numParticles());
}

void PartioSe::setTime(float val)
{
    timeVar.val=val;
}

void PartioSe::run(int i){
        currentIndex=i;
        indexVar.val=currentIndex;
        if(isPaired){
            pairIndex1=parts->data<int>(pairH1,currentIndex)[0];
            pairIndex2=parts->data<int>(pairH2,currentIndex)[0];
            if(pairIndex1 < 0 || pairIndex1 >= partsPairing->numParticles()){
                // TODO: return false
                std::cerr<<"Index "<<currentIndex<<" Referenced particle index "<<pairIndex1<<" is out of bounds"<<std::endl;
            }
            if(pairIndex2 < 0 || pairIndex2 >= partsPairing->numParticles()){
                std::cerr<<"Index "<<currentIndex<<" Referenced particle index "<<pairIndex1<<" is out of bounds"<<std::endl;
            }
        }
        SeVec3d value=evaluate();
        // map data into the particles again
        for(IntVarToPartio::iterator it=intVarToPartio.begin(),itend=intVarToPartio.end();
            it!=itend;++it)
            (*it)->mapBack();
        for(FloatVarToPartio::iterator it=floatVarToPartio.begin(),itend=floatVarToPartio.end();
            it!=itend;++it)
            (*it)->mapBack();
}

bool PartioSe::runRange(int istart,int iend)
{
    countVar.val=iend-istart+1;
    if(!isValid()){
        std::cerr<<"Not running expression because it is invalid"<<std::endl;
        std::cerr<<parseError()<<std::endl;
        return false;
    }
    if(istart < 0  || iend > parts->numParticles() || iend < istart){
        std::cerr<<"Invalid range ["<<istart<<","<<iend<<") specified. Max valid range is [0,"<<parts->numParticles()<<")"<<std::endl;
        return false;
    }
    for(int i=istart;i<iend;i++){
        run(i);
    }
    return true;
}

bool PartioSe::runRandom(){
    countVar.val=parts->numParticles();
    if(!isValid()){
        std::cerr<<"Not running expression because it is invalid"<<std::endl;
        std::cerr<<parseError()<<std::endl;
        return false;
    }
    std::vector<int> order(parts->numParticles());
    for(size_t i=0;i<order.size();i++) order[i]=i;
    for(size_t i=0;i<order.size();i++){
        size_t other=float(rand())/RAND_MAX*order.size();
        if(other>=order.size()) other=order.size()-1;
        std::swap(order[i],order[other]);
    }
    for(size_t i=0;i<order.size();i++){
        run(order[i]);
    }
    return true;
}

SeExprVarRef*  PartioSe::resolveVar(const std::string& s) const{
    {
        IntVarMap::iterator it=intVars.find(s);
        if(it != intVars.end()) return it->second;
    }
    {
        FloatVarMap::iterator it=floatVars.find(s);
        if(it != floatVars.end()) return it->second;
    }
    if(s=="PT") return &indexVar;
    if(s=="NPT") return &countVar;
    if(s=="time") return &timeVar;
    return 0;
}


PartioSe::~PartioSe()
{
    for(IntVarToPartio::iterator it=intVarToPartio.begin(),itend=intVarToPartio.end();
        it!=itend;++it)
        delete *it;
    for(FloatVarToPartio::iterator it=floatVarToPartio.begin(),itend=floatVarToPartio.end();
        it!=itend;++it)
        delete *it;
    for(IntVarMap::iterator it=intVars.begin(),itend=intVars.end();it!=itend;++it)
        delete it->second;
    for(FloatVarMap::iterator it=floatVars.begin(),itend=floatVars.end();it!=itend;++it)
        delete it->second;
}

} // end namespace Partio
