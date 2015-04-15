#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, Extension, find_packages
from Cython.Distutils import build_ext
import numpy

setup(
    name="poisson",
    version="0.1",
    description="POISSON alignment utility for nanopore sequencing data",
    author="Tamas Szalay",
    packages=find_packages(),
    install_requires=["Cython>0.18"],
    entry_points={
        'console_scripts': [
        'poisson = poisson.cmdline:main',
    ],},
    ext_modules=[Extension(name="poisson.poisscpp",
                           sources=["poisson/_poisscpp.pyx",
                                    "cpp/MakeMutations.cpp",
                                    "cpp/FindMutations.cpp",
                                    "cpp/Alignment.cpp",
                                    "cpp/swlib.cpp",
                                    "cpp/EventUtil.cpp",
                                    "cpp/Viterbi.cpp"],
                           language="c++",
                           include_dirs=[numpy.get_include(),'.'],
                           extra_compile_args=["-std=c++0x","-O3"],
                           extra_link_args=["-g"])],
    cmdclass = {'build_ext': build_ext}
)
