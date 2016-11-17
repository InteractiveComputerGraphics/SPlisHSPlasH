/*
PARTIO SOFTWARE
Copyright 2013 Disney Enterprises, Inc. All rights reserved

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
#include "../core/ParticleHeaders.h"
#include "ZIP.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <memory>

///////////////////////////////////////////
/// PCD File format  example
///<pointclouds.org> PCL point cloud library
/*
# .PCD v.7 - Point Cloud Data file format
VERSION .7
FIELDS x y z rgb
SIZE 4 4 4 4
TYPE F F F F
COUNT 1 1 1 1
WIDTH 213
HEIGHT 1
VIEWPOINT 0 0 0 1 0 0 0
POINTS 213
DATA ascii
0.93773 0.33763 0 4.2108e+06
0.90805 0.35641 0 4.2108e+06
0.81915 0.32 0 4.2108e+06.........

*/

ENTER_PARTIO_NAMESPACE

using namespace std;

// TODO: convert this to use iterators like the rest of the readers/writers

ParticlesDataMutable* readPCD(const char* filename,const bool headersOnly)
{
	//cout <<  "readPCD" << endl;
    auto_ptr<istream> input(Gzip_In(filename,ios::in|ios::binary));
    if(!*input){
        cerr<<"Partio: Can't open particle data file: "<<filename<<endl;
        return 0;
    }

    ParticlesDataMutable* simple=0;
    if(headersOnly) simple=new ParticleHeaders;
    else simple=create();

    // read NPoints and NPointAttrib
    string word;

    while(input->good())
	{
        *input>>word;
        if(word=="VERSION") break;;
    }
    // minversion is  0.7
	float versionNum = 0.0;
	if(input->good())
	{
        *input>>versionNum;
	}
	if (versionNum < 0.1)
	{
		simple->release();return 0;
	}
	//cout << "PCD version: " << versionNum << endl;
	if(input->good())
	{
        *input>>word;
        if(word!="FIELDS")
		{
			cout << "no FIELDS" << endl;
			simple->release();return 0;
		}
    }

	vector<string> attrNames;
	vector<int> attrSizes;
	vector<string> attrTypes;
	vector<string> attrCounts;
    vector<ParticleAttribute> attrs;
	unsigned int width = 0;
	unsigned int points = 0;
	bool binaryData = false;

	// we have to make some assumptions here
	while(input->good())
	{
        *input>>word;
        if(word=="SIZE") break;;
		attrNames.push_back(word);
    }
    while(input->good())
	{
        *input>>word;
        if(word=="TYPE") break;;
		if (word == "1") attrSizes.push_back(1);
		else if (word == "2") attrSizes.push_back(2);
		else if (word == "4") attrSizes.push_back(4);
		else if (word == "8") attrSizes.push_back(8);
		else
		{
			cout << "missing TYPE" << endl;
			simple->release();return 0;
		}
    }
    while(input->good())
	{
        *input>>word;
        if(word=="COUNT") break;;
		attrTypes.push_back(word);
    }
    //cout << "counts " ;
    while(input->good())
	{
        *input>>word;
        if(word=="WIDTH") break;;
		//cout << word << " ";
		attrCounts.push_back(word);
    }
    //cout << endl;
    if(input->good())
	{
		*input>>width;
	}

	while(input->good())
	{
        *input>>word;
        if(word=="POINTS") break;;
    }
    if(input->good())
	{
		*input>>points;
	}
	if(input->good())
	{
		*input>>word;
		if (word == "DATA")
		{
			if(input->good())
			{
				*input>>word;
				if(word != "ascii")
				{
					binaryData = true;
				}
			}
		}
	}

	simple->addParticles(points);
	if(headersOnly) {return simple;}// escape before we try to touch data

	// WE CONTINUE ON.... TO DATA!!!!!

	int hasPosition = 0;
	vector<string>::iterator namesIt;
	int counter = 0;
	for (namesIt = attrNames.begin(); namesIt != attrNames.end(); namesIt++)
	{

		if (*namesIt == "x" || *namesIt == "y" || *namesIt == "z")
		{
			hasPosition ++;
			if (hasPosition == 3)
			{
				attrs.push_back(simple->addAttribute("position",Partio::VECTOR,3));
			}
		}
		else
		{
			int attrCount = 1;
			if (attrCounts[counter].c_str() ==  string("1")) { attrCount = 1; }
			else if (attrCounts[counter].c_str() == string("3")) { attrCount = 3; }

			//cout << attrNames[counter].c_str() << " " << attrCount << " " << endl;

			if (attrCount == 1)
			{
				//cout << "float/int " <<  attrNames[counter].c_str() << endl;
				if (attrTypes[counter] == "I")
				{
					attrs.push_back(simple->addAttribute(attrNames[counter].c_str(),Partio::INT,1));
				}
				else
				{
					attrs.push_back(simple->addAttribute(attrNames[counter].c_str(),Partio::FLOAT,1));
				}
			}
			else if (attrCount == 3)
			{
				//cout << "vector " <<  attrNames[counter].c_str() << endl;
				attrs.push_back(simple->addAttribute(attrNames[counter].c_str(),Partio::VECTOR,3));
			}
			else
			{
				//cout << "non standard " <<  attrNames[counter].c_str() << endl;
				if (attrTypes[counter] == "I")
				{
					attrs.push_back(simple->addAttribute(attrNames[counter].c_str(),Partio::INT,attrCount));
				}
				else
				{
					attrs.push_back(simple->addAttribute(attrNames[counter].c_str(),Partio::FLOAT,attrCount));
				}
			}
		}
		counter ++;
	}


    // Read actual particle data
    if(!input->good())
	{
		simple->release();return 0;
	}
    for(unsigned int particleIndex=0;input->good() && particleIndex<points; particleIndex++)
	{
        for(unsigned int attrIndex=0;attrIndex<attrs.size();attrIndex++)
		{
            if(attrs[attrIndex].type==Partio::INT)
			{
                int* data=simple->dataWrite<int>(attrs[attrIndex],particleIndex);
                for(int count=0;count<attrs[attrIndex].count;count++)
				{
                    int ival;
                    *input>>ival;
					if (ival)
					{
						data[count]=ival;
					}
					else
					{
						data[count] = 0;
					}
                }
            }
            else if(attrs[attrIndex].type==Partio::FLOAT || attrs[attrIndex].type==Partio::VECTOR)
			{
                float* data=simple->dataWrite<float>(attrs[attrIndex],particleIndex);
                for(int count=0;count<attrs[attrIndex].count;count++)
				{
                    float fval;
                    *input>>fval;
					if (fval)
					{
						data[count]=fval;
					}
					else
					{
						data[count] = 0.0;
					}
                }
            }
        }
    }

    return simple;
}

bool writePCD(const char* filename,const ParticlesData& p,const bool compressed)
{
	//cout << "write PCD" << endl;
	bool comp = false ;
    auto_ptr<ostream> output(
        comp ?
        Gzip_Out(filename,ios::out|ios::binary)
        :new ofstream(filename,ios::out|ios::binary));

    *output<<"# .PCD v.7 - Point Cloud Data file format"<<endl;
	*output<<"VERSION .7"<<endl;
	*output<<"FIELDS ";

    vector<ParticleAttribute> attrs;
    for (int aIndex=0;aIndex<p.numAttributes();aIndex++)
	{
        attrs.push_back(ParticleAttribute());
        p.attributeInfo(aIndex,attrs[aIndex]);
		if(attrs[aIndex].name == "position")
		{
			*output<< " x y z";
		}
		else if (attrs[aIndex].name == "normal")
		{
			*output<< " normal_x normal_y normal_z";
		}
		else
		{
        *output<<" "<<attrs[aIndex].name;
		}
    }
    *output<<endl;

	*output<<"SIZE ";
	for (int aIndex=0;aIndex<p.numAttributes();aIndex++)
	{

		if(attrs[aIndex].name == "position" || attrs[aIndex].name == "normal")
		{
			*output<< " 4 4 4";
		}
		else
		{
        	switch(attrs[aIndex].type)
			{
            case FLOAT: *output<<" 4";break;
            case VECTOR: *output<<" 4";break;
            case INDEXEDSTR:
            case INT: *output<<" 4";break;
            case NONE: assert(false); break; // TODO: more graceful
			}
		}

	}
	*output<<endl;

    // TODO: assert right count
    *output<<"TYPE ";
    for (int aIndex=0;aIndex<p.numAttributes();aIndex++)
	{
		if(attrs[aIndex].name == "position" || attrs[aIndex].name == "normal")
		{
			*output<< " F F F";
		}
		else
		{
			switch(attrs[aIndex].type)
			{
            case FLOAT: *output<<" F";break;
            case VECTOR: *output<<" F";break;
            case INDEXEDSTR:
            case INT: *output<<" I";break;
            case NONE: assert(false); break; // TODO: more graceful
			}
		}
    }
    *output<<endl;
	*output<<"COUNT ";
	for (int aIndex=0;aIndex<p.numAttributes();aIndex++)
	{
		if(attrs[aIndex].name == "position" || attrs[aIndex].name == "normal")
		{
			*output<< " 1 1 1";
		}
		else
		{
			switch(attrs[aIndex].type)
			{
            case FLOAT: *output<<" 1";break;
            case VECTOR: *output<<" 3";break;
            case INDEXEDSTR:
            case INT: *output<<" 1";break;
            case NONE: assert(false); break; // TODO: more graceful
			}
		}
    }
    *output<<endl;

	*output<<"WIDTH ";
	*output<<p.numParticles();
    *output<<endl;

	*output<<"POINTS ";
	*output<<p.numParticles();
    *output<<endl;

	*output<<"DATA ascii"<<endl;

    for(int particleIndex=0;particleIndex<p.numParticles();particleIndex++)
	{
        for(unsigned int attrIndex=0;attrIndex<attrs.size();attrIndex++)
		{
            if(attrs[attrIndex].type==Partio::INT || attrs[attrIndex].type==Partio::INDEXEDSTR)
			{
                const int* data=p.data<int>(attrs[attrIndex],particleIndex);
                for(int count=0;count<attrs[attrIndex].count;count++)
                    *output<<data[count]<<" ";
            }
            else if(attrs[attrIndex].type==Partio::FLOAT || attrs[attrIndex].type==Partio::VECTOR)
			{
                const float* data=p.data<float>(attrs[attrIndex],particleIndex);
                for(int count=0;count<attrs[attrIndex].count;count++)
				{
                    *output<<data[count]<<" ";
				}
            }
        }
        *output<<endl;
    }
    return true;

}

EXIT_PARTIO_NAMESPACE
