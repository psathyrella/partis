import contextlib
import os.path
import shutil

from pkg_resources import resource_stream, resource_exists

from .. import util

_PKG = __package__ + '.data'

class InvalidLocusException(ValueError):
    pass

def _check_locus(locus):
    loci = frozenset(['igh', 'igk', 'igl'])
    llocus = locus.lower()

    if llocus not in loci:
        raise InvalidLocusException(locus)

    return llocus

def _get_file_name(locus, segment, collection=None, extension='.fasta'):
    if collection:
        return '{0}{1}-{2}{3}'.format(_check_locus(locus), segment, collection,
                                      extension)
    else:
        return '{0}{1}{2}'.format(_check_locus(locus), segment, extension)

def _handle(locus, segment, collection, extension, vdj_dir=None):
    file_name = _get_file_name(locus, segment, collection=collection,
                               extension=extension)
    if vdj_dir != None:
        return open(vdj_dir + '/' + file_name)
    else:
        return resource_stream(_PKG, file_name)

def fasta_handle(locus, segment, collection=None, vdj_dir=None):
    return _handle(locus, segment, collection, extension='.fasta', vdj_dir=vdj_dir)

def gff3_handle(locus, segment, collection=None):
    return _handle(locus, segment, collection, extension='.gff3')

def has_gff3(locus, segment, collection=None):
    file_name= _get_file_name(locus, segment, collection=collection,
                              extension='.gff3')
    return resource_exists(_PKG, file_name)

@contextlib.contextmanager
def temp_fasta(locus, segment, collection=None, vdj_dir=None):
    """
    Copy an IGH segment collection to a temporary file
    """
    base = os.path.splitext(_get_file_name(locus, segment, collection))[0]
    with fasta_handle(locus, segment, collection=collection, vdj_dir=vdj_dir) as ifp, \
            util.tempdir(prefix=base) as td:
        with open(td(base + '.fasta'), 'w') as ofp:
            shutil.copyfileobj(ifp, ofp)
        yield ofp.name
