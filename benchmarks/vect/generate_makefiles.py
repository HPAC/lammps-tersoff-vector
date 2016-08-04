# Build all the binaries for our vectorization experiment, and copy them to bin/
# Input Command Line Argument: LAMMPS root directory
# Of course, the compilation environment has to be set up correctly.
# As a word of warning: This compilation process might take some time.
# Only invoke from its directory, as this program will look there for template files.

import sys
import os

archs = ["-mmic", "-xAVX", "-xCORE-AVX2", "-xSSE4.2"]

my_dir = os.getcwd()
template = open('template-Makefile.intel_eval', "r").read()
os.chdir(sys.argv[1])
os.chdir("src")
lmp_dir = os.getcwd()

def gen2(opt, name):
    with open("MAKE/MINE/Makefile.intel_eval", 'w') as f:
        f.write('CCCFG=%s\n' % opt)
        if opt.startswith("-mmic"):
            f.write("MKLKIND=mic\n")
        else:
            f.write("MKLKIND=intel64\n")
        f.write(template)
    os.system("touch pair_tersoff_intel.cpp")
    os.system("make intel_eval -j16")
    os.system("mv lmp_intel_eval "+my_dir)
    os.chdir(my_dir)
    os.system("mv lmp_intel_eval bin/lmp_intel_eval_"+name)
    os.chdir(lmp_dir)

def gen(opt, name):
    gen2(opt + " -DLMP_INTEL_TERSOFF_PACK_I=true ", name+"i")
    gen2(opt + " -DLMP_INTEL_TERSOFF_PACK_I=false", name+"j")

for opt in archs:
    os.system("make clean-intel_eval")
    # generate scalar code
    gen(opt + " -DLMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE", opt[2]+"s")
    # generate vectorized code
    gen(opt, opt[2]+"v")
