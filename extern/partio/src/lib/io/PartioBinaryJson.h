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
#include <map>
#include <iostream>
#include "PartioEndian.h"
//#include <src/lib/io/PartioEndian.h>
#include <cstdio>
#include <typeinfo>
#include <sstream>
#include <stdexcept>
ENTER_PARTIO_NAMESPACE

/// TODO: figure out what this isfor
template<class TDERIVED>
struct JSONDelegate{
    
};

/// format a string as hex for debugging
template<class T>
std::string hex(const T val)
{
    std::stringstream ss;
    ss<<"0x"<<std::hex<<val;
    return ss.str();
}

/// state of the parser (where in the file it is)
struct JSONParserState
{
    std::map<long long,std::string> tokens; /// code to string
    std::istream& fp; /// file pointer

    JSONParserState(std::istream& fp)
        :fp(fp)
    {}

    /// where in the file am I?
    size_t getOffset() const{
        return fp.tellg();
    }
private:
    /* this class should nto be copyable */
    JSONParserState(const JSONParserState&);
    const JSONParserState& operator=(const JSONParserState&);
};

/// Parse error exception
class JSONParseError : public std::runtime_error
{
    std::string offsetStr(const JSONParserState& state){
        size_t offset=state.getOffset();
        char buf[1024];
        sprintf(buf,"offset %d ",offset);
        return buf;
    }
public:
    JSONParseError(const std::string& s,const JSONParserState& state): runtime_error(hex(state.getOffset()) +": "+s) {}
};


/// JSON parser. Calls through static polymorphism into the delegate
/// this is similar to many XML parsers structure (except we avoid virtuals)
/// it is duck typed... see printparserdelegate
template<class DELEGATE>
struct JSONParser
{

    typedef unsigned char uint8;
    typedef char int8;
    typedef unsigned short uint16;
    typedef short int16;
    typedef unsigned int uint32;
    typedef int int32;
    typedef long long int64;
    typedef float real32;
    typedef double real64;

    static const uint32 BINARY_MAGIC            = 0x624a534e;
    static const uint32 BINARY_MAGIC_SWAP       = 0x4e534a62;


    DELEGATE& Derived(){
        return static_cast<DELEGATE&>(*this);
    }

    enum TOKEN{
        JID_NULL                = 0x00,
        JID_MAP_BEGIN           = 0x7b,
        JID_MAP_END             = 0x7d,
        JID_ARRAY_BEGIN         = 0x5b,
        JID_ARRAY_END           = 0x5d,
        JID_BOOL                = 0x10,
        JID_INT8                = 0x11,
        JID_INT16               = 0x12,
        JID_INT32               = 0x13,
        JID_INT64               = 0x14,
        JID_REAL16              = 0x18, // TODO: this one
        JID_REAL32              = 0x19,
        JID_REAL64              = 0x1a,
        JID_UINT8               = 0x21,
        JID_UINT16              = 0x22,
        JID_STRING              = 0x27,
        JID_FALSE               = 0x30,
        JID_TRUE                = 0x31,
        JID_TOKENDEF            = 0x2b,
        JID_TOKENREF            = 0x26,
        JID_TOKENUNDEF          = 0x2d, // TODO: handle this
        JID_UNIFORM_ARRAY       = 0x40,
        JID_KEY_SEPARATOR       = 0x3a, // TODO: these don't seem to occur maybe in AScii
        JID_VALUE_SEPARATOR     = 0x2c, // TODO: these also don't seem to occur maybe in ascii
        JID_MAGIC               = 0x7f
    };

    TOKEN currToken;
    //DELEGATE delegate;
    JSONParserState& state;
    const char* currKey;
    JSONParser(JSONParserState& state)
        :state(state),currKey(0)
    {
        //nextToken();
    }

    std::string readString(){
        int64 strLength=readLength();
//        std::cerr<<"str length is "<<strLength<<std::endl;
        std::string s;
        char* buf=new char[strLength+1];
        buf[strLength]=0;
        state.fp.read(buf,strLength);
        std::string sAsStr=buf;
        return sAsStr;

    }

    void defineToken(){
        int64 id=readLength();
        //std::cerr<<"id is "<<id<<std::endl;
        std::string s=readString();
        //std::cerr<<"string is \""<<s<<"\""<<std::endl;
        state.tokens[id]=s;
    }

    void undefineToken(){
        // TODO: do
    }

    int64 tokenRefNoNotify(){
        int64 id=readLength();
        return id;
    }

    void tokenRef(){
        int64 id=readLength();
        // TODO: need check if token exists
        //std::cerr<<"token ref "<<tokens[id]<<std::endl;
        Derived().string(currKey,state.tokens[id]);
    }

    bool nextToken()
    {
        // TODO: bette buffering
        unsigned char c;
        state.fp>>c;
        currToken=TOKEN(c);
        while(currToken==JID_TOKENDEF || currToken==JID_TOKENUNDEF){
//            std::cerr<<"iteration of def/undef"<<std::endl;
            if(currToken==JID_TOKENDEF) defineToken();
            else undefineToken();
            state.fp>>c;
            currToken=TOKEN(c);
        }
        //std::cerr<<"NEXT TOKEN "<<currToken<<" "<<hex(currToken)<<std::endl;
//        std::cerr<<"doen with tokendef"<<std::endl;
        return true;
    }

    void parse()
    {
        // first check for magic
        //unsigned char c;
        //Partio::read<BIGEND>(state.fp,c);
        nextToken();
        printf("c 0x%02x\n",currToken);
        if(currToken==JID_MAGIC){
            uint32 magic;
            Partio::read<BIGEND>(state.fp,magic);
            printf("magic 0x%08x\n",magic);
            if(magic == BINARY_MAGIC_SWAP){
                std::cerr<<"swap which is big endian"<<std::endl;
            }else if(magic==BINARY_MAGIC){
                std::cerr<<"little endian not supported"<<std::endl;
            }else{
                std::cerr<<"bad magic # "<<magic<<std::endl;
            }
            nextToken();
        }else{
            std::cerr<<"no magic so ascii"<<std::endl;
        }

        value();
    }

    int64 readLength(){
        uint8 lengthCode;
        uint32 length32;
        int64 length64;
        read<LITEND>(state.fp,lengthCode);
        //std::cerr<<"length code "<<int(lengthCode)<<std::endl;
        if(lengthCode<0xf1) return lengthCode;
        else if(lengthCode==0xf2){
            uint16 length16;
            read<LITEND>(state.fp,length16);
            //std::cerr<<"   16bit length "<<length16<<std::endl;
            return length16;
        }else if(lengthCode==0xf4){
            uint32 length32;
            read<LITEND>(state.fp,length32);
            return length32;
        }else if(lengthCode==0xf8){
            int64 length64;
            read<LITEND>(state.fp,length64);
            return length64;
        }

        // TODO: error check
        throw JSONParseError(std::string("unexpected length code ")+hex(lengthCode),state);
    }

    /// parse an array map  of data i.e.   [<key>,<value>,<key>,value>,...]
    /// (expects currToken to have been JID_ARRAY_BEGIN)
    void arrayMap(){
        while(1){
            nextToken();
            if(currToken==JID_ARRAY_END) break;
            requireKey();

            nextToken();
            TOKEN val=value();
        }
        Derived().arrayEnd();
    }

    /// parse an array of data (expects currToken to have been JID_ARRAY_BEGIN)
    void array(){
        while(1){
            nextToken();
            if(currToken==JID_ARRAY_END) break;
            value();
        }
        Derived().arrayEnd();
    }

    /// check to make sure currToken
    void requireKey()
    {
        if(currToken!= JID_TOKENREF){
            throw JSONParseError("needed a token ref got token "+hex(currToken),state);
        }
        int64 token=tokenRefNoNotify();
        currKey=state.tokens[token].c_str();
    }

    void map(){
        for(;;){
            nextToken();
            if(currToken == JID_MAP_END) break;
            requireKey();
            nextToken();
            TOKEN finalValue=value();
            currKey=0;
        }
        Derived().mapEnd();
    }

    template<class T>
    T numberData(){
        T val;
        read<LITEND>(state.fp,val);
        return val;
    }

    template<class T>
    void numberDataNotify(){
        T val;
        read<LITEND>(state.fp,val);
        Derived().numberData<T>(currKey,val);
    }

    template<class T> void uniformArray(){
        int64 length=readLength();
        // give delegate chance to do something
        if(!Derived().uniformArray<T>(currKey,length)){
            uniformArrayDefaultImpl<T>(currKey,length);
        }
    }
    template<class T> bool uniformArrayDefaultImpl(const char* currKey,int length){
        int count=length; // /sizeof(T);
        for(int64 i=0;i<count;i++) numberDataNotify<T>();
        return true;
    }

    void uniformArray(){
        nextToken();
        switch(currToken){
            case JID_REAL32: uniformArray<real32>();break;
            case JID_REAL64: uniformArray<real64>();break;
            case JID_INT8: uniformArray<int8>();break;
            case JID_UINT8: uniformArray<uint8>();break;
            case JID_INT16: uniformArray<int16>();break;
            case JID_UINT16: uniformArray<uint16>();break;
            case JID_INT32: uniformArray<int32>();break;
            case JID_INT64: uniformArray<int64>();break;
            case JID_TOKENREF: break; // TODO: implement
            case JID_STRING: break; // TODO: implement
            default: throw JSONParseError("UNKNOWN DATA TYPE "+hex(currToken),state);
        }
    }

    void boolData(){
        uint8 val;
        read<LITEND>(state.fp,val);;
        bool result=val != 0;
        Derived().boolData(currKey,result);
    }

    TOKEN value(){
        TOKEN currTokenCopy=currToken;
        switch(currTokenCopy){
            case JID_ARRAY_BEGIN: 
                if(!Derived().arrayBegin(currKey)) array();              
                return JID_ARRAY_BEGIN;break;
            case JID_MAP_BEGIN: 
                if(!Derived().mapBegin(currKey)) map();              
               return JID_MAP_BEGIN;break;
            case JID_TOKENREF: tokenRef();return JID_TOKENREF;break;
            case JID_BOOL: boolData();return JID_BOOL;break;
            case JID_TRUE: 
                Derived().boolData(currKey,true);
                return JID_BOOL;break;
            case JID_FALSE: 
                Derived().boolData(currKey,false);
                return JID_BOOL;break;
            case JID_REAL32: numberDataNotify<real32>();return JID_REAL32;break;
            case JID_REAL64: numberDataNotify<real64>();return JID_REAL64;break;
            case JID_INT32: numberDataNotify<int32>();return JID_INT32;break;
            case JID_INT16: numberDataNotify<int16>();return JID_INT16;break;
            case JID_INT8: numberDataNotify<int8>();return JID_INT8;break;
            case JID_UINT16: numberDataNotify<uint16>();return JID_UINT16;break;
            case JID_UINT8: numberDataNotify<uint8>();return JID_UINT8;break;
            case JID_UNIFORM_ARRAY: uniformArray();return JID_UNIFORM_ARRAY;break;
            case JID_STRING: Derived().string(currKey,readString());break;
                
            default:
                throw JSONParseError("Unknown token 0x"+hex(currToken),state);
        }
        return JID_NULL;
    }

};


/// This parser ignores everything it sees 
/// it serves as an example of how a delegate might be written
struct NULLParser:public JSONParser<NULLParser>
{
    NULLParser(JSONParserState& state)
        :JSONParser<NULLParser>(state)
    {}

    /// handle string data
    void string(const char* key,const std::string& s){}
    /// handle a numeric data of a specific type
    template<class T> void numberData(const char* key,T x){}
    /// handle bool data
    void boolData(const char* key,bool x){}
    /// handle map data (Returns true ifthis parser has handled reading through another delegate)
    bool mapBegin(const char* key){return false;}
    /// see the end of a map
    void mapEnd(){}
    /// begin of array return true if you have already read the array data to through another delegate
    bool arrayBegin(const char* key){return false;}
    /// end of arrays
    void arrayEnd(){}
    /// got a uniform array of data
    template<class T> bool uniformArray(const char* currKey,int length){
        return false;
    }
};

namespace{
std::string keyString(const char* key) {
    return std::string("key=")+(key?key:"null");
}
}

/// This parser prints everything it sees
struct PrintParser:public JSONParser<PrintParser>
{
    int _indent;

    PrintParser(JSONParserState& state)
        :JSONParser<PrintParser>(state),_indent(0)
    {}


    void indent(){
        for(int i=0;i<_indent*4;i++) std::cerr<<" ";
    }

    void string(const char* key,const std::string& s){
        indent();
        std::cerr<<keyString(key)<<" string "<<s<<std::endl;
    }
    template<class T> void numberData(const char* key,T x){
        indent();
        std::cerr<<keyString(key)<<" num "<<x<<std::endl;
    }
    void boolData(const char* key,bool x){
        indent();
        std::cerr<<keyString(key)<<" num "<<x<<std::endl;
    }
    bool mapBegin(const char* key){
        indent();
        std::cerr<<keyString(key)<<" {"<<std::endl;;
        _indent++;
        return false;
    }
    void mapEnd(){
        _indent--;
        indent();
        std::cerr<<"]"<<std::endl;
    }
    bool arrayBegin(const char* key){
        indent();
        std::cerr<<keyString(key)<<" ["<<std::endl;;
        _indent++;
        return false;
    }
    void arrayEnd(){
        _indent--;
        indent();
        std::cerr<<"]"<<std::endl;
    }
    template<class T> bool uniformArray(const char* key,int length){
        std::cerr<<"uniform array begin key="<<keyString(key)<<" length "<<length<<std::endl;
        return false;
    }
};

EXIT_PARTIO_NAMESPACE
