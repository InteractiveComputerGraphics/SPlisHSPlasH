import os
import re
import sys
import platform
import subprocess
import multiprocessing as mp
import argparse

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

# Extract cmake arguments
parser = argparse.ArgumentParser()
parser.add_argument("-D", action='append', dest='cmake',
                    help="CMake Options")
parser.add_argument("--manylinux-build", action='store_true', dest='manylinux_build')
args, other_args = parser.parse_known_args(sys.argv)
cmake_clargs = args.cmake
sys.argv = other_args

# Project binding name
name = "pySPlisHSPlasH"
internal_name = "pysplishsplash"


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        bin_dir_windows = os.path.join(os.path.abspath(self.build_temp), "bin")
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        # Add cmake command line arguments
        if cmake_clargs is not None:
            cmake_args += ['-D{}'.format(arg) for arg in cmake_clargs]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                           '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + bin_dir_windows,
                           '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), bin_dir_windows)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j{}'.format(mp.cpu_count())]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        # Add position independent code flags if using gcc on linux probably
        if platform.system() == "Linux":
            cmake_args += ['-DCMAKE_CXX_FLAGS=-fPIC', '-DCMAKE_C_FLAGS=-fPIC']

            # Using relative rpath messes up repairing the wheel file. The relative rpath is only necessary when
            # building locally from source
            if not args.manylinux_build:
                cmake_args += ['-DCMAKE_INSTALL_RPATH={}'.format("$ORIGIN"),
                               '-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON',
                               '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=OFF']

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        if not args.manylinux_build:
            subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
            subprocess.check_call(['cmake', '--build', '.', '--target', internal_name] + build_args, cwd=self.build_temp)
        else:
            subprocess.check_call(['cmake3', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
            subprocess.check_call(['cmake3', '--build', '.', '--target', internal_name] + build_args, cwd=self.build_temp)

        # Copy dlls to ext directory so they are installed alongside the bindings
        if platform.system() == "Windows":
            subprocess.check_call(['cmake', '-E', 'copy', os.path.join(bin_dir_windows, "glew.dll"), extdir])
            subprocess.check_call(['cmake', '-E', 'copy', os.path.join(bin_dir_windows, "freeglut.dll"), extdir])


# List the files that should be installed alongside the package
models = [os.path.join('data/models/', file) for file in os.listdir('data/models/')]
Scenes = [os.path.join('data/Scenes/', file) for file in os.listdir('data/Scenes/')]
shaders = [os.path.join('data/shaders/', file) for file in os.listdir('data/shaders/')]
emitter_boundary = [os.path.join('data/emitter_boundary/', file) for file in os.listdir('data/emitter_boundary/')]

# Install paths depending on system
models_dest = 'data/models' if platform.system() == "Windows" else "bin/data/models"
scenes_dest = 'data/Scenes' if platform.system() == "Windows" else "bin/data/Scenes"
shaders_dest = 'resources/shaders' if platform.system() == "Windows" else 'bin/resources/shaders'
emitter_boundary_dest = 'resources/emitter_boundary' if platform.system() == "Windows" else 'bin/resources/emitter_boundary'

# Get Readme text for long description
cur_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(cur_dir, "README.md"), 'r') as f:
    long_description = f.read()

setup(
    name=name,
    version='2.7.2',
    author='Interactive Computer Graphics',
    author_email='',
    description='SPlisHSPlasH Project Python Bindings',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/InteractiveComputerGraphics/SPlisHSPlasH',
    license="MIT",
    keywords="sph fluids sph-fluids smoothed-particle-hydrodynamics fluid-simulation fluid-dynamics multiphase-flow viscous-fluids deformable-solids simulation",
    ext_modules=[CMakeExtension(name)],
    cmdclass=dict(build_ext=CMakeBuild),
    packages=find_packages(),
    entry_points={'console_scripts': 'splash = pySPlisHSPlasH.scripts.simulator:main'},
    data_files=[(models_dest, [m for m in models if os.path.isfile(m)]),
                (scenes_dest, [s for s in Scenes if os.path.isfile(s)]),
                (shaders_dest, [s for s in shaders if os.path.isfile(s)]),
				(emitter_boundary_dest, [s for s in emitter_boundary if os.path.isfile(s)])],
    zip_safe=False,
)