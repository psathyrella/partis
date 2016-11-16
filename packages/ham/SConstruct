from os import getenv, environ
import glob
import os

env = Environment(ENV=os.environ)
yaml_cpp = env.Command(['_build/libyaml-cpp.a'],
                       glob.glob('yaml-cpp/src/{*.cpp,*.h}'),
                       ['mkdir -p _build/yaml-cpp; cd _build/yaml-cpp; cmake ../../yaml-cpp; make yaml-cpp; ln -sf yaml-cpp/libyaml-cpp.a ../'])

VariantDir('_build', 'src')

SConscript('_build/SConscript', duplicate=0)

SConscript('test/SConscript', duplicate=0)

# What just `scons` will build.
Default('hample')

# `scons test`
Alias('test', 'test/_results/ALL.passed')
