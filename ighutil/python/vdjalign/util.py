import bz2
import gzip
import collections
import contextlib
import functools
import os.path
import shutil
import sys
import tempfile

# From https://github.com/lh3/readfq
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


@contextlib.contextmanager
def tempdir(**kwargs):
    td = tempfile.mkdtemp(**kwargs)
    try:
        yield functools.partial(os.path.join, td)
    finally:
        shutil.rmtree(td)


@contextlib.contextmanager
def tmpfifo(**kwargs):
    """
    Context manager - yields a temporary filesystem FIFO object, with optional
    named argument ``name`` from ``**kwargs``. Remaining arguments passed to
    :func:`tempdir`.

    On exiting the context manager, the temporary directory containing the FIFO
    is deleted.
    """
    fifo_name = kwargs.pop('name', 'fifo')
    with tempdir(**kwargs) as j:
        f = j(fifo_name)
        os.mkfifo(f)
        yield f


def opener(mode, *args, **kwargs):
    """
    Open a file, with optional compression based on extension
    """
    exts = {'.bz2': bz2.BZ2File,
            '.gz': gzip.open}
    def open_file(path):
        if path == '-':
            if mode.startswith('r'):
                return sys.stdin
            else:
                return sys.stdout

        open_fn = exts.get(os.path.splitext(path)[1], open)
        return open_fn(path, mode, *args, **kwargs)
    return open_file

# From http://wiki.python.org/moin/PythonDecoratorLibrary
class memoized(object):
   """Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).
   """
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)


@contextlib.contextmanager
def maybe_with(fn, arg=None, **kwargs):
    """
    if ``arg`` is truthy, yields the managed object resulting from applying context
    manager ``fn`` to ``arg`` with extra arguments ``**kwargs``.
    Otherwise, yields None.
    """
    if not arg:
        yield arg
    else:
        with fn(arg, **kwargs) as result:
            yield result

@contextlib.contextmanager
def with_if(condition, fn, *args, **kwargs):
    """
    If ``condition`` yields the managed object resulting from applying context
    manager ``fn`` with arguments ``*args`` and ``**kwargs``, otherwise yields
    None.
    """
    if condition:
        with fn(*args, **kwargs) as result:
            yield result
    else:
        yield

