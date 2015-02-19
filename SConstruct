SConscript('test/SConscript', duplicate=0)

# `scons test`
Alias('test', 'test/_results/ALL.passed')

Alias('validate', '_output/validation/valid.out')
Command('_output/validation/valid.out', './bin/run-driver.py', './bin/run-driver.py --label validation --plotdir _output/validation/plots --datafname test/A-every-100-subset-0.tsv.bz2 && touch $TARGET')
