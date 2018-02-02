/*
PARTIO SOFTWARE
Copyright 2011 Disney Enterprises, Inc. All rights reserved

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
#ifndef _READERS_h_
#define _READERS_h_

namespace Partio{
ParticlesDataMutable* readBGEO(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readGEO(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPDB(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPDB32(const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPDB64(const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPDA(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readMC(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPTC(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPDC(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPRT(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readBIN(	const char* filename,const bool headersOnly,std::ostream* errorStream);
ParticlesDataMutable* readPTS(  const char* filename,const bool headersOnly,std::ostream* errorStream);

bool writeBGEO(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writeGEO(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePDB(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePDB32(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePDB64(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePDA(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePTC(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writeRIB(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePDC(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writePRT(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
bool writeBIN(const char* filename,const ParticlesData& p,const bool compressed,std::ostream* errorStream);
}

#endif
