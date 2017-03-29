import utils
from subprocess import check_output

basedatadir = '/fh/fast/matsen_e/data'
baseprocdatadir = '/fh/fast/matsen_e/processed-data/partis/clustering-paper'

datasets = ['vollmers', 'adaptive', 'stern']

humans = {}
humans['vollmers'] = ['021-018', '021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-081', '021-084']
humans['adaptive'] = ['A', 'B', 'C']
humans['stern'] = ['SRR1383326', 'SRR1383447', 'SRR1383448', 'SRR1383449', 'SRR1383450', 'SRR1383451', 'SRR1383452', 'SRR1383453', 'SRR1383454', 'SRR1383455', 'SRR1383456', 'SRR1383457', 'SRR1383458', 'SRR1383459', 'SRR1383460', 'SRR1383461', 'SRR1383462', 'SRR1383463', 'SRR1383464', 'SRR1383465', 'SRR1383466', 'SRR1383467', 'SRR1383468', 'SRR1383469', 'SRR1383470', 'SRR1383471', 'SRR1383472', 'SRR1383473', 'SRR1383474', 'SRR1383475', 'SRR1383476', 'SRR1383477']

dataset_dirs = {}
dataset_dirs['vollmers'] = 'vollmers'
dataset_dirs['adaptive'] = 'adaptive-billion-read'
dataset_dirs['stern'] = '2016-01-06-stern-bcell/proc_data'
colors = {
    'A': '595', 'B': '807', 'C': '834',
    # '021-018': '1', '021-019': '1', '021-044': '1', '021-048': '1', '021-050': '1', '021-055': '1', '021-057': '1', '021-059': '1',
    # '021-060': '1', '021-061': '1', '021-063': '1', '021-068': '1', '021-071': '1', '021-081': '1', '021-084': '1'
    '021-018': '8', '021-019': '9', '021-044': '28', '021-048': '30', '021-050': '33', '021-055': '34', '021-057': '36', '021-059': '38',
    '021-060': '39', '021-061': '41', '021-063': '42', '021-068': '45', '021-071': '46', '021-081': '47', '021-084': '49'
}

all_subdirs = [ '.', ] \
              + [ e + '_del' for e in utils.real_erosions ] \
              + [ i + '_insertion' for i in utils.boundaries] \
              + [ 'mute-freqs', ] \
              + [ 'mute-freqs/' + r for r in utils.regions ]

def get_datafname(human, dataset=None):
    if dataset is None:
        dataset = get_dataset(human)
    basepath = basedatadir + '/' + dataset_dirs[dataset] + '/' + human
    if dataset == 'adaptive':
        return basepath + '/shuffled.csv'
    elif dataset == 'vollmers':
        return basepath + '/' + human + '_Lineages.fasta'
    elif dataset == 'stern':
        return basepath + '_collapse-unique_atleast-2.fastq'
    else:
        assert False

def get_outdir(human, dataset=None):
    if dataset is None:
        dataset = get_dataset(human)
    return baseprocdatadir + '/' + dataset_dirs[dataset] + '/' + human

def get_nseqs(human, dataset=None):
    fname = get_datafname(human, dataset)
    n_lines = int(check_output(['wc', '-l', fname]).split()[0])
    suffix = fname.split('.')[-1]
    if suffix == 'fasta' or suffix == 'fastq':
        return n_lines / 2
    elif suffix == 'csv':
        return n_lines - 1
    else:
        raise Exception('bad suffix %s in %s' % (suffix, fname))

def get_dataset(human):
    for dset in datasets:
        if human in humans[dset]:
            return dset
