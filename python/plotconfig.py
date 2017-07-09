import utils

rstrings = ['', 'cdr3_'] + [r + '_' for r in utils.regions]
gene_usage_columns = [r + '_gene' for r in utils.regions]
int_columns = [c for c in utils.index_columns if c not in gene_usage_columns]

xtitles = {
    'v_gene' : '',
    'd_gene' : '',
    'j_gene' : '',
    'hamming_to_true_naive' : 'hamming',
    'v_hamming_to_true_naive' : 'hamming',
    'd_hamming_to_true_naive' : 'hamming',
    'j_hamming_to_true_naive' : 'hamming',
    'cdr3_hamming_to_true_naive' : 'hamming',
    'seq_content' : '',
    'dj_insertion_content' : '',
    'vd_insertion_content' : '',
    'mute_freqs' : 'mutation freq',
    'v_mute_freqs' : 'mutation freq',
    'd_mute_freqs' : 'mutation freq',
    'j_mute_freqs' : 'mutation freq',
    'muted_bases' : 'total mutations',
    'v_muted_bases' : 'V mutations',
    'd_muted_bases' : 'D mutations',
    'j_muted_bases' : 'J mutations',
    'cdr3_muted_bases' : 'CDR3 mutations',
    'all-mean-freq' : 'mutation freq',
    'v-mean-freq' : 'mutation freq',
    'd-mean-freq' : 'mutation freq',
    'j-mean-freq' : 'mutation freq',
    'v_fraction_correct_vs_mute_freq' : 'mut freq',
    'd_fraction_correct_vs_mute_freq' : 'mut freq',
    'j_fraction_correct_vs_mute_freq' : 'mut freq'
}
for rstr in [r + '_' for r in utils.regions] + ['', ]:
    xtitles[rstr + 'hamming_to_true_naive'] = 'hamming distance'

plot_titles = {
    'v_gene' : 'V gene',
    'd_gene' : 'D gene',
    'j_gene' : 'J gene',
    'hamming_to_true_naive' : 'total distance to true naive',
    'v_hamming_to_true_naive' : 'V distance to true naive',
    'd_hamming_to_true_naive' : 'D distance to true naive',
    'j_hamming_to_true_naive' : 'J distance to true naive',
    'cdr3_hamming_to_true_naive' : 'distance to true naive within CDR3',
    'seq_content' : 'sequence base content',
    'dj_insertion_content' : 'DJ insert base content',
    'vd_insertion_content' : 'VD insert base content',
    'mute_freqs' : 'sequence mutation freq',
    'v_mute_freqs' : 'V mutation freq',
    'd_mute_freqs' : 'D mutation freq',
    'j_mute_freqs' : 'J mutation freq',
    'muted_bases' : 'total mutations',
    'v_muted_bases' : 'V mutations',
    'd_muted_bases' : 'D mutations',
    'j_muted_bases' : 'J mutations',
    'cdr3_muted_bases' : 'CDR3 mutations',
    'all-mean-freq' : 'sequence mutation freq',
    'v-mean-freq' : 'V mutation freq',
    'd-mean-freq' : 'D mutation freq',
    'j-mean-freq' : 'J mutation freq',
    'v_fraction_correct_vs_mute_freq' : 'V gene',
    'd_fraction_correct_vs_mute_freq' : 'D gene',
    'j_fraction_correct_vs_mute_freq' : 'J gene',
    'v_allele_fraction_correct_vs_per_gene_support' : 'V gene',
    'd_allele_fraction_correct_vs_per_gene_support' : 'D gene',
    'j_allele_fraction_correct_vs_per_gene_support' : 'J gene',
    'shm_indel_length' : 'total/net SHM indel length',
}
for region in utils.regions:
    for end in ['5', '3']:
        xtitles[region + '_' + end + 'p_del'] = 'bases'
        plot_titles[region + '_' + end + 'p_del'] = region.upper() + ' ' + end + '\' deletion'
for boundary in utils.all_boundaries:
    xtitles[boundary + '_insertion'] = 'bases'
    plot_titles[boundary + '_insertion'] = boundary.upper() + ' N length'


true_vs_inferred_hard_bounds = {
    # gene call
    'XXX_allele_fraction_correct_vs_per_gene_support' : (0., 1.),
    'XXX_fraction_correct_vs_mute_freq' :  (0, 1./3),
    # mutation
    'v_hamming_to_true_naive' : (-0.5, 4.5),
    'd_hamming_to_true_naive' : (-0.5, 5.5),
    'j_hamming_to_true_naive' : (-0.5, 5.5),
    'cdr3_hamming_to_true_naive' : (-0.5, 10.5),
    'hamming_to_true_naive' : (-0.5, 10.5),
    'v_muted_bases' : (-2.5, 2.5),
    'd_muted_bases' : (-4.5, 4.5),
    'j_muted_bases' : (-5.5, 5.5),
    'cdr3_muted_bases' : (-4.5, 4.5),
    'mute_freqs' : (-.015, .015),  # NOTE make sure you know where the decimal place is here!
    'muted_bases' : (-5.5, 5.5),
    # boundaries
    'v_3p_del' : (-3.5, 3.5),
    'd_5p_del' : (-8.5, 8.5),
    'd_3p_del' : (-8.5, 8.5),
    'j_5p_del' : (-10.5, 15.5),
    'vd_insertion' : (-8.5, 8.5),
    'dj_insertion' : (-10.5, 15.5)
}
for name, bounds in true_vs_inferred_hard_bounds.items():
    if name.find('XXX_') == 0:
        for rstr in rstrings:  # adds some we don't need, but that's ok
            true_vs_inferred_hard_bounds[name.replace('XXX_', rstr)] = bounds
        del true_vs_inferred_hard_bounds[name]

default_hard_bounds = {
    # 'hamming_to_true_naive' : (-0.5, 19.5),
    'hamming_to_true_naive' : (-0.5, 10),
    'v_hamming_to_true_naive' : (-0.5, 8.5),
    'd_hamming_to_true_naive' : (-0.5, 10.5),
    'j_hamming_to_true_naive' : (-0.5, 12.5),
    'd_3p_del' : (-1, 15),
    'd_5p_del' : (-1, 18),
    'dj_insertion' : (-1, 13),
    'jf_insertion' : (-1, 13),
    'fv_insertion' : (-1, 13),
    'j_5p_del' : (-1, 15),
    'all-mean-freq' : (0.0, 0.4),  # NOTE make sure you know where the decimal place is here!
    'v-mean-freq' : (0.0, 0.4),  # NOTE make sure you know where the decimal place is here!
    'd-mean-freq' : (0.0, 0.5),  # NOTE make sure you know where the decimal place is here!
    'j-mean-freq' : (0.0, 0.4),  # NOTE make sure you know where the decimal place is here!
    'v_3p_del' : (-1, 6),
    'vd_insertion' : (-1, 15),
    'all_insertions' : (-0.5, 20),
    'IGHJ6*02' : (-0.5, 39.5),
    'IGHJ3*02' : (-0.5, 26),
    'IGHJ1*01' : (-0.5, 28),
    'IGHJ2*01' : (-0.5, 28),
    'IGHJ3*01' : (-0.5, 25),
    'IGHJ4*01' : (-0.5, 23),
    'IGHJ4*02' : (-0.5, 24),
    'IGHJ4*03' : (-0.5, 24),
    'IGHJ5*01' : (-0.5, 27),
    'IGHJ5*02' : (-0.5, 28)
    # 'IGHV3-23D*01' : (260, 296),
    # 'IGHV3-33*06' : (260, 296)
}
