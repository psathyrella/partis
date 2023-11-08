fixers="fissix.fixes.fix_except \
fissix.fixes.fix_ne \
fissix.fixes.fix_numliterals \
fissix.fixes.fix_repr \
fissix.fixes.fix_standarderror \
fissix.fixes.fix_tuple_params \
libmodernize.fixes.fix_basestring \
libmodernize.fixes.fix_import \
libmodernize.fixes.fix_imports_six \
libmodernize.fixes.fix_raise \
libmodernize.fixes.fix_basestring \
libmodernize.fixes.fix_open \
fissix.fixes.fix_has_key \
libmodernize.fixes.fix_dict_six \
"
fixers=libmodernize.fixes.fix_dict_six
# python/alleleclusterer.py python/allelefinder.py python/alleleremover.py python/annotationclustering.py python/bar.py python/baseutils.py python/baz.py python/cached_uncertainties.py python/clusterpath.py python/coar.py python/corrcounter.py python/event.py python/foo.py python/fraction_uncertainty.py python/gex.py python/glomerator.py python/glutils.py python/hist.py python/hmmwriter.py python/hutils.py python/indelutils.py python/__init__.py python/lbplotting.py python/mds.py python/mutefreqer.py python/paircluster.py python/parametercounter.py python/paramutils.py python/partitiondriver.py python/partitionplotter.py python/performanceplotter.py python/plotconfig.py python/plotting.py python/processargs.py python/prutils.py python/recombinator.py python/scanplot.py python/seqfileopener.py python/treegenerator.py python/treeutils.py python/utils.py python/viterbicluster.py python/waterer.py
paths="python/*.py test bin/*.py datascripts/*.py projects/*.py"
for fixer in $fixers; do
    modernize -wnf $fixer $paths
done
