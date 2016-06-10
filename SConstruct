env = Environment(
    CPPPATH = ['src'],
    CCFLAGS='-std=c++11 -g',
#    LIBS='pll',
    LINKFLAGS='-g')

env.Program(
    target='test',
    source=[Glob('src/*.cpp')])
