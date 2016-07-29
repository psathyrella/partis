import utils

xtitles = {
    'v_gene' : '',
    'd_gene' : '',
    'j_gene' : '',
    'hamming_to_true_naive' : 'hamming',
    'v_hamming_to_true_naive' : 'hamming',
    'd_hamming_to_true_naive' : 'hamming',
    'j_hamming_to_true_naive' : 'hamming',
    'v_hamming_to_true_naive_normed' : '% hamming',
    'd_hamming_to_true_naive_normed' : '% hamming',
    'j_hamming_to_true_naive_normed' : '% hamming',
    'd_3p_del' : 'bases',
    'd_5p_del' : 'bases',
    'dj_insertion' : 'bases',
    'dj_insertion_content' : '',
    'j_5p_del' : 'bases',
    'mute_freqs' : 'mutation freq',
    'v_mute_freqs' : 'mutation freq',
    'd_mute_freqs' : 'mutation freq',
    'j_mute_freqs' : 'mutation freq',
    'v_3p_del' : 'bases',
    'vd_insertion' : 'bases',
    'vd_insertion_content' : '',
    'all-mean-freq' : 'mutation freq',
    'v-mean-freq' : 'mutation freq',
    'd-mean-freq' : 'mutation freq',
    'j-mean-freq' : 'mutation freq',
    
}

plot_titles = {
    'v_gene' : 'V gene',
    'd_gene' : 'D gene',
    'j_gene' : 'J gene',
    'hamming_to_true_naive' : 'Distance to true naive',
    'v_hamming_to_true_naive' : 'V Distance to true naive',
    'd_hamming_to_true_naive' : 'D Distance to true naive',
    'j_hamming_to_true_naive' : 'J Distance to true naive',
    'v_hamming_to_true_naive_normed' : 'V Distance to true naive',
    'd_hamming_to_true_naive_normed' : 'D Distance to true naive',
    'j_hamming_to_true_naive_normed' : 'J Distance to true naive',
    'dj_insertion_content' : 'DJ insert base content',
    'mute_freqs' : 'sequence mutation freq',
    'v_mute_freqs' : 'V mutation freq',
    'd_mute_freqs' : 'D mutation freq',
    'j_mute_freqs' : 'J mutation freq',
    'vd_insertion_content' : 'VD insert base content',
    'all-mean-freq' : 'sequence mutation freq',
    'v-mean-freq' : 'V mutation freq',
    'd-mean-freq' : 'D mutation freq',
    'j-mean-freq' : 'J mutation freq',
    'v_gene_fraction_vs_mute_freq' : 'V gene',
    'd_gene_fraction_vs_mute_freq' : 'D gene',
    'j_gene_fraction_vs_mute_freq' : 'J gene'
}
for region in utils.regions:
    for end in ['5', '3']:
        plot_titles[region + '_' + end + 'p_del'] = region.upper() + ' ' + end + '\' deletion'
for boundary in utils.boundaries + utils.effective_boundaries:
    plot_titles[boundary + '_insertion'] = boundary.upper() + ' N length'


true_vs_inferred_hard_bounds = {
    'hamming_to_true_naive' : (-0.5, 19.5),
    'v_hamming_to_true_naive' : (-0.5, 8.5),
    'd_hamming_to_true_naive' : (-0.5, 10.5),
    'j_hamming_to_true_naive' : (-0.5, 12.5),
    'v_hamming_to_true_naive_normed' : (-0.5, 8.5),
    'd_hamming_to_true_naive_normed' : (-0.5, 50.5),
    'j_hamming_to_true_naive_normed' : (-0.5, 20),
    'd_3p_del' : (-8.5, 8.5),
    'd_5p_del' : (-8.5, 8.5),
    'dj_insertion' : (-10.5, 15.5),
    'j_5p_del' : (-10.5, 15.5),
    'mute_freqs' : (-.05, .05),  # NOTE make sure you know where the decimal place is here!
    'v_3p_del' : (-3.5, 3.5),
    'vd_insertion' : (-8.5, 8.5),
    'v_gene_fraction_vs_mute_freq' : (0, 1./3),
    'd_gene_fraction_vs_mute_freq' : (0, 1./3),
    'j_gene_fraction_vs_mute_freq' : (0, 1./3)
}

default_hard_bounds = {
    # 'hamming_to_true_naive' : (-0.5, 19.5),
    'hamming_to_true_naive' : (-0.5, 10),
    'v_hamming_to_true_naive' : (-0.5, 8.5),
    'd_hamming_to_true_naive' : (-0.5, 10.5),
    'j_hamming_to_true_naive' : (-0.5, 12.5),
    'v_hamming_to_true_naive_normed' : (-0.5, 8.5),
    'd_hamming_to_true_naive_normed' : (-0.5, 50.5),
    'j_hamming_to_true_naive_normed' : (-0.5, 20),
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
