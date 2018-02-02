Partio - A library for particle IO and manipulation
===================================================

This is the initial source code release of partio a tool we used for particle
reading/writing.  It started out as an abstraction for the commonalities in
particle models (i.e. accessing many attributes associated with an index or
entity).

Super impatient building guide
==============================

    # Install Location ~ adjust accordingly
    prefix=$HOME/local
    # Partio depends on SeExpr v2.11
    git clone https://github.com/wdas/seexpr.git
    pushd seexpr
      git checkout v2.11
      make -j prefix=$prefix install
    popd
    git clone https://github.com/wdas/partio.git
    cd partio
    make -j prefix=$prefix RP_SeExpr=$prefix install

Getting Started
===============

CMake is used to build the project, but we provide a top-level Makefile
for convenience that takes care of all the steps.

See the Makefile for the user-tweakable variables and corresponding
cmake options.

The typical usage for an installation into `/usr/local`
with a temporary staging directory of `/tmp/stage` is:

    make DESTDIR=/tmp/stage RP_SeExpr=/usr/local prefix=/usr/local install

SeExpr
======
Partio depends on SeExpr v2 (current `v2.11`).  SeExpr v2 and v3 can co-exist,
so you can have v2 for building partio installed alongside SeExpr v3+ (master).

SeExpr is not strictly required. SeExpr support is disabled when `RP_SeExpr`
is not specified to the Makefile.  In CMake, the corresponding variables are
`-DSEEXPR_BASE=<path>` for the SeExpr location, and `-DPARTIO_SE_ENABLED=1`
to enable the PartioSe class.

Source code overview
====================

src/
   lib/      Library code (public API in root)
   lib/core  Core library (KDtree traversal, data representations)
   lib/io    Input/Output (Different file formats)
   py/       SWIG based python bindings
   doc/      Doxygen documentation and (the start of) a manual
   tests/    Start of regression tests (I need more)
   tools/    Useful tools 
             partconv <input format> <output format>
             partinfo <particle file>
             partview <particle file>

Class Model
-----------

The goal of the library is to abstract the particle interface from the data
representation. That is why Partio represents particles using three classes that
inherit and provide more functionality

ParticlesInfo - Information about # of particles and attributes
ParticlesData - Read only access to all particle data
ParticlesDataMutable - Read/write access to all particle data

The functions used to get particle access are these:

readHeaders() 
   returns ParticlesInfo
   reads only the minimum data necessary to get number of particles and
   attributes

readCached()
   returns ParticlesData
   For multiple users in different threads using the same particle file
   ParticlesData

create() and read()
   returns ParticlesDataMutable   
   allows read/write access

Behind the scenes you could implement these classes however you like.  Headers
only representation is called core/ParticleHeader.{h,cpp}.  Simple
non-interleaved attributes is core/ParticleSimple.{h,cpp}.

Attribute Data Model
--------------------

All particles have the same data attributes.  They have the model that they are
of three basic types with a count of how many scalar values they have.

VECTOR[3]    
FLOAT[d]
INT[d]

VECTOR[3] and FLOAT[3] have the same data representations.
VECTOR[4] is invalid however FLOAT[4] is valid as is FLOAT[1...infinity]

This seems to encompass the most common file formats for particles


Iterating
---------

There are multiple ways to access data in the API. Here are
some tips

- Use SIMD functions when possible prefer dataAsFloat(),data(arrayOfIndices) as
  opposed to data(int singleIndex) which accesses multiple pieces of data at
  once

- Cache ParticleAttributes for quick access instead of calling attributeInfo()
  over a loop of particles

- Use iterators to do linear operations over all particles They are much more
  optimized than both data() and the dataAsFloat or


Backends
--------

Behind the scenes there are SimpleParticles, ParticleHeaders, and
SimpleParticlesInterleaved.  In the future I would like to write a disk-based
cached back end that can dynamically only load the data that is necessary.
create(), read() and readCached could be augmented to create different
structures in these cases.

Readers/Writers
---------------

New readers and writers can be added in the io/ directory. You simply need to
implement the interface ParticlesInfo, ParticlesData and ParticlesDataMutable
(or as many as you need). Editing the io/readers.h to add prototypes and
io/ParticleIO.cpp to add file extension bindings should be easy.


- Andrew Selle, Walt Disney Animation Studios
