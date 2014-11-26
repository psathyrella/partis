from os import getenv, environ
from glob import glob

# http://www.scons.org/wiki/PhonyTargets
def PhonyTarget(target, action):
    phony = Environment(ENV = environ,
                        BUILDERS = { 'phony' : Builder(action = action) })
    AlwaysBuild(phony.phony(target = target, source = 'SConstruct'))


VariantDir('_build/yaml-cpp', 'yaml-cpp')
VariantDir('_build', 'src')

SConscript('_build/yaml-cpp/SConscript', duplicate=0)
SConscript('_build/SConscript', duplicate=0)

SConscript('test/SConscript', duplicate=0)

# What just `scons` will build.
Default('hample')

# `scons test`
Alias('test', 'test/_results/ALL.passed')

# Phony targets
astyle = """
astyle  -A2 \
        --pad-oper \
        --unpad-paren \
        --keep-one-line-blocks \
        --keep-one-line-statements \
        --suffix=none \
        --formatted \
        --lineend=linux \
        --indent=spaces=2 \
"""
ham_files = glob(getenv('PWD')+'/src/*.cc')+glob(getenv('PWD')+'/include/*.h')

# `scons style`
PhonyTarget('style', astyle + " ".join(ham_files))

