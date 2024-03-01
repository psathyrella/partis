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
libmodernize.fixes.fix_map \
libmodernize.fixes.fix_filter \
libmodernize.fixes.fix_zip \
libmodernize.fixes.fix_xrange_six \
libmodernize.fixes.fix_next \
libmodernize.fixes.fix_print \
libmodernize.fixes.fix_input_six \
libmodernize.fixes.fix_file \
fissix.fixes.fix_apply \
fissix.fixes.fix_reduce \
fissix.fixes.fix_exec \
fissix.fixes.fix_paren \
fissix.fixes.fix_funcattrs
"
# fixers=libmodernize.fixes.fix_dict_six
paths="/loc/dralph/work/linearham/scripts/*.py" #"datascripts/meta/*/*.py" #"python/*.py test bin/*.py datascripts/*.py projects/*.py"
for fixer in $fixers; do
    modernize -wnf $fixer $paths
done
