#from distutils.core import setup
#from distutils.extension import Extension
from setuptools import setup, Extension, find_packages
from Cython.Distutils import build_ext
import numpy

setup(
    name="poreseq",
    version="0.1",
    description="PoreSeq alignment utility for nanopore sequencing data",
    author="Tamas Szalay",
    packages=find_packages(),
    install_requires=["Cython>0.18","Biopython>1.6", "pysam>0.8.2", "h5py>2.0"],
    scripts=[
        'scripts/poreseq_align',
        'scripts/poreseq_assemble',
    ],
    data_files=[
        ('resources', ['poreseq.spec', 'scripts/convertFastaAndQualToFastq.jar']),
    ],
    zip_safe=False, # so that the resources can be fetched by the scripts
    entry_points={
        'console_scripts': [
            'poreseq = poreseq.cmdline:main',
        ]
    },
    ext_modules=[Extension(name="poreseq.poreseqcpp",
                           sources=["poreseq/_poreseqcpp.pyx",
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
