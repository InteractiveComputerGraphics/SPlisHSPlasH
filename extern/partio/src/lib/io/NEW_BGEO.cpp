/*
Copyright Disney Enterprises, Inc. All rights reserved.

This license governs use of the accompanying software. If you use the software, you
accept this license. If you do not accept the license, do not use the software.

1. Definitions
The terms "reproduce," "reproduction," "derivative works," and "distribution" have
the same meaning here as under U.S. copyright law. A "contribution" is the original
software, or any additions or changes to the software. A "contributor" is any person
that distributes its contribution under this license. "Licensed patents" are a
contributor's patent claims that read directly on its contribution.

2. Grant of Rights
(A) Copyright Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free copyright license to reproduce its contribution, prepare
derivative works of its contribution, and distribute its contribution or any derivative
works that you create.
(B) Patent Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free license under its licensed patents to make, have made,
use, sell, offer for sale, import, and/or otherwise dispose of its contribution in the
software or derivative works of the contribution in the software.

3. Conditions and Limitations
(A) No Trademark License- This license does not grant you rights to use any
contributors' name, logo, or trademarks.
(B) If you bring a patent claim against any contributor over patents that you claim
are infringed by the software, your patent license from such contributor to the
software ends automatically.
(C) If you distribute any portion of the software, you must retain all copyright,
patent, trademark, and attribution notices that are present in the software.
(D) If you distribute any portion of the software in source code form, you may do
so only under this license by including a complete copy of this license with your
distribution. If you distribute any portion of the software in compiled or object code
form, you may only do so under a license that complies with this license.
(E) The software is licensed "as-is." You bear the risk of using it. The contributors
give no express warranties, guarantees or conditions. You may have additional
consumer rights under your local laws which this license cannot change.
To the extent permitted under your local laws, the contributors exclude the
implied warranties of merchantability, fitness for a particular purpose and non-
infringement.
*/
#include "PartioBinaryJson.h"
#include <fstream>
#include <cstdlib>
#include <cstring>
#include "../Partio.h"

ENTER_PARTIO_NAMESPACE

struct BGEOAttributeType:public JSONParser<BGEOAttributeType>
{
    std::string attrName;

    BGEOAttributeType(JSONParserState& state)
        :JSONParser<BGEOAttributeType>(state)
    {}

    bool arrayBegin(const char* key){
        throw JSONParseError("unexpected array begin",state);
    }

    void string(const char* key,const std::string& s){
        //std::cerr<<"key "<<keyString(key)<<" string "<<s<<std::endl;
        if(strcmp(key,"name")==0){
            attrName=s;
        }
    }
    template<class T> void numberData(const char* key,T x){
        throw JSONParseError("unexpected number data",state);
    }
    void boolData(const char* key,bool x){
        throw JSONParseError("unexpected number data",state);
    }
    bool mapBegin(const char* key){
        //std::cerr<<"map "<<key<<std::endl;
        return false;}
    void mapEnd(){
        //std::cerr<<"mapend"<<std::endl;
    }
    void arrayEnd(){
    }

    template<class T> bool uniformArray(const char* currKey,int length){
        return false;
    }


};
 
struct BGEOAttributeValues:public JSONParser<BGEOAttributeValues>
{
    const char* attrName;
    std::string storage;
    int size;
    Partio::ParticlesDataMutable* particles;
    std::vector<int> packing;

    BGEOAttributeValues(JSONParserState& state,Partio::ParticlesDataMutable* particles,const char* attrName)
        :JSONParser<BGEOAttributeValues>(state),attrName(attrName),storage(""),size(0),particles(particles)
    {
        std::cerr<<"BGEOAttributeValues"<<std::endl;
    }


    void string(const char* key,const std::string& s){
        if(strcmp(key,"storage")==0){
            storage=s;
        }
        std::cerr<<"key "<<keyString(key)<<" string "<<s<<std::endl;
    }
    template<class T> bool uniformArray(const char* key,int length){
        std::cerr<<"uniform array starting with "<<keyString(key)<<" lenght "<<length<<std::endl;
        //for(int64 i=0;i<length;i++){
        //    T val=((JSONParser<BGEOAttributeValues>*)this)->numberData<T>();
        //    std::cerr<<"  key "<<keyString(key)<<" "<<val<<std::endl;
        //}
        //return true;
        if(strcmp(key,"rawpagedata")==0){
            // TODO: check for size and storage
            int nparts=particles->numParticles();
            std::cerr<<"attr "<<attrName<<" reading raw page data with packing ";
            for(size_t k=0;k<packing.size();k++) std::cerr<<" "<<packing[k];
            std::cerr<<" storage "<<storage<<" length "<<length<<std::endl;
            if(strcmp(attrName,"P")==0) attrName="position";


            if(storage.substr(0,2)=="fp"){
                Partio::ParticleAttribute attrH=particles->addAttribute(attrName,Partio::FLOAT,size);
                int num=length/size;
                int offset=0;
                if(packing.size()==0) packing.push_back(size); // if no packing then full
                for(size_t pack=0;pack<packing.size();pack++){
                    int packsize=packing[pack];
                    for(int i=0;i<num;i++){
                        float* partioVal=particles->dataWrite<float>(attrH,i);
                        for(int k=0;k<packsize;k++){
                            T val=((JSONParser<BGEOAttributeValues>*)this)->numberData<T>();
                            partioVal[offset+k]=val;
                        }
                    }
                    offset+=packing[pack];
                }
            }else if(storage.substr(0,3)=="int" || storage.substr(0,4)=="uint"){
                Partio::ParticleAttribute attrH=particles->addAttribute(attrName,Partio::INT,size);
                int num=length/size;
                int offset=0;
                if(packing.size()==0) packing.push_back(size); // if no packing then full
                for(size_t pack=0;pack<packing.size();pack++){
                    int packsize=packing[pack];
                    for(int i=0;i<num;i++){
                        float* partioVal=particles->dataWrite<float>(attrH,i);
                        for(int k=0;k<packsize;k++){
                            T val=((JSONParser<BGEOAttributeValues>*)this)->numberData<T>();
                            partioVal[offset+k]=val;
                        }
                    }
                    offset+=packing[pack];
                }
            }else{
                throw JSONParseError("Unexpected data storage "+storage,state);
            }
            return true;
        }else if(strcmp(key,"packing")==0){
            for(int i=0;i<length;i++){
                T val=((JSONParser<BGEOAttributeValues>*)this)->numberData<T>();
                packing.push_back(val);
            }
            return true;
        }
        return false;
    }
    template<class T> void numberData(const char* key,T x){
        std::cerr<<"key "<<keyString(key)<<" "<<int64(x)<<std::endl;
        if(strcmp(key,"size")==0){size=x;}
    }
    void boolData(const char* key,bool x){
        std::cerr<<"key "<<keyString(key)<<" "<<x<<std::endl;
    }
    bool mapBegin(const char* key){
        throw JSONParseError("Got unknown map begin with key "+keyString(key),state);
        return false;}
    void mapEnd(){
        std::cerr<<"mapend"<<std::endl;
    }
    bool arrayBegin(const char* key){
        throw JSONParseError("Got unknown key in values block in attribute"+keyString(key),state);
        return false;
        //std::cerr<<"arraybegin "<<key<<std::endl;
    }
    void arrayEnd(){}
};



struct BGEOAttributeValue:public JSONParser<BGEOAttributeValue>
{
    Partio::ParticlesDataMutable* particles;
    const char* attrName;
    BGEOAttributeValue(JSONParserState& state,Partio::ParticlesDataMutable* particles,const char* attrName)
        :JSONParser<BGEOAttributeValue>(state),particles(particles),attrName(attrName)
    {
        std::cerr<<"BGEOAttributeValue attr="<<attrName<<std::endl;
    }


    void string(const char* key,const std::string& s){
        std::cerr<<"key "<<keyString(key)<<" string "<<s<<std::endl;
    }
    template<class T> void numberData(const char* key,T x){
        std::cerr<<"key "<<keyString(key)<<" "<<int64(x)<<std::endl;
    }
    void boolData(const char* key,bool x){
        std::cerr<<"key "<<keyString(key)<<" "<<x<<std::endl;
    }
    bool mapBegin(const char* key){
        std::cerr<<"map "<<key<<std::endl;
        return false;}
    void mapEnd(){
        std::cerr<<"mapend"<<std::endl;
    }
    bool arrayBegin(const char* key){
        if(strcmp(key,"defaults")==0){
            NULLParser(state).arrayMap();
            return true;
        }else if(strcmp(key,"values")==0){
            //NULLParser(state).arrayMap();
            BGEOAttributeValues(state,particles,attrName).arrayMap();
            return true;
        }
        throw JSONParseError("Got unknown key in values block in attribute"+keyString(key),state);
        return false;
        //std::cerr<<"arraybegin "<<key<<std::endl;
    }
    void arrayEnd(){
        std::cerr<<"arrayend"<<std::endl;
    }
    template<class T> bool uniformArray(const char* currKey,int length){
        return false;
    }

};

struct BGEOAttributes:public JSONParser<BGEOAttributes>
{
    Partio::ParticlesDataMutable* particles;
    enum ATTRTYPE{NONE,POINT};
    ATTRTYPE attrtype;
    enum ATTRSTATE{TOP,ATTRS,MID,MID2,MID3};
    ATTRSTATE attrstate;
    std::string attrName;
    
    BGEOAttributes(JSONParserState& state,Partio::ParticlesDataMutable* particles)
        :JSONParser<BGEOAttributes>(state),particles(particles),attrtype(NONE),attrstate(TOP)
    {}

    template<class T> bool uniformArray(const char* currKey,int length){
        return false;
    }


    bool arrayBegin(const char* key){
        if(attrstate==TOP){
            if(attrtype==NONE){
                if(strcmp(key,"pointattributes")==0){
                    std::cerr<<"parsing "<<key<<std::endl;
                    attrtype=POINT;
                    attrstate=ATTRS;
                    return false;
                }else{
                    std::cerr<<"ignoring "<<key<<std::endl;
                    NULLParser(state).array();
                    return true;
                }
            }
        }else if(attrstate==ATTRS){
            attrstate=MID;
            return false;
        }else if(attrstate==MID){
            BGEOAttributeType attrType(state);
            attrType.arrayMap();
            attrName=attrType.attrName;
            attrstate=MID2;
            return true;
        }else if(attrstate==MID2){
            BGEOAttributeValue(state,particles,attrName.c_str()).arrayMap();
            attrstate=MID3;
            return true;
        }
        std::stringstream ss;
        ss<<"unknown key "<<keyString(key)<<" at type "<<attrtype<<" attrstate "<<attrstate;
        throw JSONParseError(ss.str(),state);
    }
    void arrayEnd()
    {
        if(attrstate==MID){
            throw JSONParseError("Trying to exit array in MID state",state);
        }else if(attrstate==MID2){
            throw JSONParseError("Trying to exit array in MID2 state",state);
        }else if(attrstate==MID3){
            attrstate=ATTRS;
        }else if(attrstate==ATTRS){
            attrtype=NONE;
            attrstate=TOP;
        }
    }


    void string(const char* key,const std::string& s){

        throw JSONParseError("unexpected string data",state);
    }
    template<class T> void numberData(const char* key,T x){
        throw JSONParseError("unexpected number data",state);
    }
    void boolData(const char* key,bool x){
        throw JSONParseError("unexpected number data",state);
    }
    bool mapBegin(const char* key){return false;}
    void mapEnd(){}


};

struct BGEOMainParser:public JSONParser<BGEOMainParser>
{
    int64 pointcount,vertexcount,primitivecount;
    Partio::ParticlesDataMutable* particles;

    BGEOMainParser(JSONParserState& state, ParticlesDataMutable* particles)
        :JSONParser<BGEOMainParser>(state),particles(particles)
    {}

    ~BGEOMainParser(){
        std::cerr<<"pt count "<<pointcount<<" vert "<<vertexcount<<" prim count "<<primitivecount<<std::endl;
    }

    void string(const char* key,const std::string& s){
        //if(strcmp(key,"fileversion")==0) std::cerr<<"saw file version "<<s<<std::endl;
    }
    template<class T> void numberData(const char* key,T x){
        if(strcmp(key,"primitivecount")==0) primitivecount=x;
        else if(strcmp(key,"vertexcount")==0) vertexcount=x;
        else if(strcmp(key,"pointcount")==0){
            pointcount=x;
            particles->addParticles(pointcount);
        }else{
            std::cerr<<"WARNING unknown number data ignored in key "<<keyString(key)<<std::endl;
        }
    }
    void boolData(const char* key,bool x){
        std::cerr<<"WARNING ignoring boolData with key="<<keyString(key)<<" value="<<x<<std::endl;
    }
    bool mapBegin(const char* key){
        if(strcmp(key,"info")==0){
            NULLParser parser(state);
            parser.map();
            return true;
        }else{
            std::cerr<<"WARNING unknown key "<<keyString(key)<<" with type array"<<std::endl;
        }
        std::cerr<<"WARNING found mapbegin with key "<<key<<std::endl;
        return false;
    }
    void mapEnd(){throw JSONParseError("logic error",state);}
    bool arrayBegin(const char* key){
        if(strcmp(key,"topology")==0){
            NULLParser parser(state);
            parser.arrayMap();
            return true;
        }else if(strcmp(key,"attributes")==0){
            BGEOAttributes(state,particles).arrayMap();
            return true;
        }else if(strcmp(key,"primitives")==0){
            NULLParser parser(state);
            parser.array();
            return true;
        }else if(strcmp(key,"pointgroups")==0){
            NULLParser parser(state);
            parser.array();
            return true;
        }
        std::cerr<<"WARNING found arraybegin with key "<<key<<std::endl;
        return false;
    }
    void arrayEnd(){}

    template<class T> bool uniformArray(const char* currKey,int length){
        return false;
    }

};

/// This parser ignores everything it sees
struct BGEOParser:public JSONParser<BGEOParser>
{
    Partio::ParticlesDataMutable* particles;

    BGEOParser(JSONParserState& state,Partio::ParticlesDataMutable* particles)
        :JSONParser<BGEOParser>(state),particles(particles)
    {}

    void string(const char* key,const std::string& s){throw JSONParseError("invalid string at top level",state);}
    template<class T> void numberData(const char* key,T x){throw JSONParseError("invalid numberData at top level",state);}
    void boolData(const char* key,bool x){throw JSONParseError("invalid boolData  at top level",state);}
    bool mapBegin(const char* key){throw JSONParseError("invalid mapBegin at top level",state);return false;}
    void mapEnd(){throw JSONParseError("invalid mapEnd at top level",state);}
    bool arrayBegin(const char* key){
        std::cerr<<"starting array"<<std::endl;
        BGEOMainParser main(state,particles);
        main.arrayMap();
        return true;
    }
    void arrayEnd(){throw JSONParseError("invalid arrayEnd at top level",state);}
    template<class T> bool uniformArray(const char* currKey,int length){
        return false;
    }
};

ParticlesDataMutable*  testRead(const char* filename)
{
    std::ifstream fp(filename,std::ios::in);
    if(!fp){
        std::cerr<<"Error opening file"<<std::endl;
        exit(1);
    }
    Partio::JSONParserState state(fp);
    //Partio::GEOParser parser(state);
    Partio::ParticlesDataMutable* particles=Partio::create();
    Partio::BGEOParser parser(state,particles);
    //Partio::PrintParser parser(state);
    //PrintDelegate delegate;
    //Partio::JSONParser<PrintDelegate> json(delegate,fp);
    parser.parse();
    std::cerr<<"particle count "<<particles->numParticles()<<std::endl;
    std::cerr<<"attr count "<<particles->numAttributes()<<std::endl;
    return particles;
}

EXIT_PARTIO_NAMESPACE
