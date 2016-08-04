import collections
import csv
import urllib

from .util import memoized

# Bits for handling attribute encoding
_ENCODE_CHARS = ',=;\t'
_ENCODING_MAP = {k: urllib.quote(k) for k in _ENCODE_CHARS}

def attribute_encode(s):
    return ''.join(_ENCODING_MAP.get(i, i) for i in s)

def parse_attributes(attr_string):
    kvs = (i.split('=', 1) for i in attr_string.split(';'))
    result = collections.OrderedDict()
    for k, v in kvs:
        result[urllib.unquote(k)] = urllib.unquote(v)

    return result

def encode_attributes(attr_dict):
    result = []
    for k, v in attr_dict.iteritems():
        result.append(''.join((attribute_encode(k), '=', attribute_encode(v))))
    return ';'.join(result)

# GFF3 records
_COL_CLASSES = {'start': int, 'end': int}

_GFF3Base = collections.namedtuple('GFF3Record',
        ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand',
         'phase', 'attributes'])

class GFF3Record(_GFF3Base):
    """
    Primitive GFF3 representation
    """
    @property
    def attr(self):
        return self.attribute_dict()

    @memoized
    def attribute_dict(self):
        return parse_attributes(self.attributes)

    def update_attributes(self, **kwargs):
        """Update attributes"""
        attr = self.attribute_dict()
        attr.update(**kwargs)
        return self._replace(attributes=encode_attributes(attr))

    @classmethod
    def _make(cls, row):
        """
        Add some types
        """
        if len(row) != len(cls._fields):
            # Let namedtuple handle
            return super(GFF3Record, cls)._make(row)

        res = []
        for field, v in zip(cls._fields, row):
            if field in _COL_CLASSES:
                v = _COL_CLASSES[field](v)
            res.append(v)

        return super(GFF3Record, cls)._make(res)

    @property
    def start0(self):
        return self.start - 1

def parse(fp):
    lines = (line for line in fp if not line.startswith('#'))
    reader = csv.reader(lines, delimiter='\t')
    for row in reader:
        yield GFF3Record._make(row)
