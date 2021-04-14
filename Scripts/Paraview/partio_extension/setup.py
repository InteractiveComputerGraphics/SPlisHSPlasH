from setuptools import setup

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11 import get_cmake_dir

import sys
import platform

__version__ = "0.0.1"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)


defs = []
cxx_args = []

defs.append(('PARTIO_USE_ZLIB', None))

if platform.system() == 'Windows':
    defs.append(('PARTIO_WIN32', None))
    defs.append(('_USE_MATH_DEFINES', None))
elif platform.system() == 'Linux':
    cxx_args = ["-fPIC", "-w"]

ext_modules = [
    Pybind11Extension("partio",
                      ["partio_bindings.cpp"] +
                      [
                          '../../../extern/partio/src/lib/core/Particle.cpp',
                          '../../../extern/partio/src/lib/core/ParticleCaching.cpp',
                          '../../../extern/partio/src/lib/core/ParticleHeaders.cpp',
                          '../../../extern/partio/src/lib/core/ParticleSimple.cpp',
                          '../../../extern/partio/src/lib/core/ParticleSimpleInterleave.cpp',
                          '../../../extern/partio/src/lib/io/BGEO.cpp',
                          '../../../extern/partio/src/lib/io/BIN.cpp',
                          '../../../extern/partio/src/lib/io/GEO.cpp',
                          '../../../extern/partio/src/lib/io/MC.cpp',
                          '../../../extern/partio/src/lib/io/ParticleIO.cpp',
                          '../../../extern/partio/src/lib/io/PDA.cpp',
                          '../../../extern/partio/src/lib/io/PDB.cpp',
                          '../../../extern/partio/src/lib/io/PDC.cpp',
                          '../../../extern/partio/src/lib/io/PRT.cpp',
                          '../../../extern/partio/src/lib/io/PTC.cpp',
                          '../../../extern/partio/src/lib/io/PTS.cpp',
                          '../../../extern/partio/src/lib/io/RIB.cpp',
                          '../../../extern/partio/src/lib/io/ZIP.cpp',
                          '../../../extern/zlib/src/adler32.c',
                          '../../../extern/zlib/src/compress.c',
                          '../../../extern/zlib/src/crc32.c',
                          '../../../extern/zlib/src/deflate.c',
                          '../../../extern/zlib/src/gzio.c',
                          '../../../extern/zlib/src/infback.c',
                          '../../../extern/zlib/src/inffast.c',
                          '../../../extern/zlib/src/inflate.c',
                          '../../../extern/zlib/src/inftrees.c',
                          '../../../extern/zlib/src/trees.c',
                          '../../../extern/zlib/src/uncompr.c',
                          '../../../extern/zlib/src/zutil.c'
                      ],
                      include_dirs=['../../../extern/partio/src/lib', '../../../extern/zlib/src'],
                      # Example: passing in the version to the compiled code
                      define_macros=[('VERSION_INFO', __version__)] + defs,
                      cxx_std=14
                      ),
]

setup(
    name="partio",
    version=__version__,
    author="Stefan Jeske",
    author_email="jeske@cs.rwth-aachen.de",
    description="Python Bindings for Partio using pybind11",
    long_description="",
    ext_modules=ext_modules,
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)