from os import getenv, environ
import glob
import os

env = Environment(ENV=os.environ)

VariantDir('_build', 'src')

SConscript('_build/SConscript', duplicate=0)

SConscript('test/SConscript', duplicate=0)

# What just `scons` will build.
Default('hample')

# `scons test`
Alias('test', 'test/_results/ALL.passed')
