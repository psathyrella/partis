from .codons import translate

def ref_map(list ap, int qstart=0):
    """Build a map from reference base index to query base index"""
    cdef dict result = {}
    for i, j in ap:
        if j is not None:
            result[j] = i + qstart if i is not None else None
    return result

def fst(*args):
    for arg in args:
        if arg is not None:
            return arg

def intersect_feature(dict rmap, feature):
    """Update a GFF feature to refer to query coordinates"""
    cdef int start0 = feature.start0, end0 = feature.end - 1
    cdef int rmin = min(rmap), rmax = max(rmap)
    cdef int any_overlap = rmin <= end0 and rmax >= start0
    cdef int complete_overlap = rmin <= start0 and rmax >= end0

    if not any_overlap:
        return None

    cdef int sstart0 = fst(rmap.get(rmin if rmin > start0 else start0), -1)
    cdef int send0 = fst(rmap.get(rmax if rmax < end0 else end0), -1)

    if sstart0 < 0 or send0 < 0: # gap in query at start or end
        return None

    return feature.update_attributes(complete_overlap=str(complete_overlap))\
            ._replace(start=sstart0 + 1, end=send0 + 1)

def intersect_features(dict rmap, list features):
    cdef dict r = {}
    for feature in features:
        n = feature.attribute_dict()['Name']
        r[n] = intersect_feature(rmap, feature)
    return r

def classify_record(v, d, j, dict v_annot, dict j_annot):
    vm = ref_map(v.aligned_pairs, v.qstart)
    jm = ref_map(j.aligned_pairs, j.qstart)

    # Translate features
    v_trans = intersect_features(vm, v_annot.values())
    j_trans = intersect_features(jm, j_annot.values())

    result = {}

    # CDR3 length
    cys = v_trans['Cys']
    tryp = j_trans['J_Tryp']

    # Find maximum aligned base in v at or before the cysteine
    vmax_r, vmax_q = max(((r, q) for r, q in vm.iteritems()
                          if r is not None and q is not None and r < v_annot['Cys'].start0))

    result['cdr3_start'] = v_annot['Cys'].start0 - (vmax_r - vmax_q)
    result['cdr3_end'] = j_annot['J_Tryp'].end - j.pos + j.qstart
    result['frame'] = result['cdr3_start'] % 3
    result['amino_acid'] = translate(v.seq[result['frame']:])
    result['cdr3_length'] = result['cdr3_end'] - result['cdr3_start']
    result['cdr3_aa'] = translate(v.seq[result['cdr3_start']:result['cdr3_end']])
    result['cdr3_coverage'] = (cys and tryp and
        cys.attr['complete_overlap'] == '1' and
        tryp.attr['complete_overlap'] == '1')

    features = ((n, feat) for i in [v_trans, j_trans]
                for n, feat in i.iteritems())

    for (gene, alignment) in [('v', v), ('d', d), ('j', j)]:
        result['{0}_nm'.format(gene)] = alignment.opt('NM')
        result['{0}_qstart'.format(gene)] = alignment.qstart
        result['{0}_qend'.format(gene)] = alignment.qend

    for n, feature in features:
        result['{0}_begin'.format(n)] = feature.start0 if feature else None
        result['{0}_end'.format(n)] = feature.end if feature else None
        result['{0}_complete'.format(n)] = feature.attr['complete_overlap'] == '1' if feature else None

    return result
