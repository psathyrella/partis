import random
import re

import pytest

from partis import disjointgrouper
from partis.disjointgrouper import _split_procs


def test_one_big_group_gets_most_threads():
    # big group gets every thread the small ones do not need to end before it, the rest one thread each
    threads, n_jobs = _split_procs([2100000] + [50000] * 19, 16)
    assert threads[0] == 12
    assert set(threads[1:]) == {1}
    assert n_jobs == 5
    assert sum(threads[:n_jobs]) == 16


def test_small_groups_wait_so_the_largest_get_two_threads():
    # 15 groups on 16 procs: one thread each would leave the largest running long after the small ones end
    sizes = [451, 425, 411, 406, 394, 279, 271, 248, 206, 189, 155, 105, 98, 93, 67]
    threads, n_jobs = _split_procs(sizes, 16)
    assert threads == [2] * 5 + [1] * 10
    assert n_jobs == 11


def test_equal_small_groups_get_one_thread_each():
    threads, n_jobs = _split_procs([1000] * 40, 16)
    assert set(threads) == {1}
    assert n_jobs == 16


def test_single_group_gets_whole_budget():
    assert _split_procs([500000], 16) == ([16], 1)


def test_invariants_on_random_sizes():
    rng = random.Random(12345)
    for _ in range(200):
        budget = rng.randint(1, 64)
        sizes = sorted([rng.randint(1, 5000000) for _ in range(rng.randint(1, 60))], reverse=True)
        threads, n_jobs = _split_procs(sizes, budget)
        assert len(threads) == len(sizes)
        assert min(threads) >= 1
        assert 1 <= n_jobs <= len(sizes)
        assert sum(threads[:n_jobs]) <= budget


def test_round2_threads_follow_the_squared_centroid_counts(tmp_path, monkeypatch):
    # on the counts this plan would be 4, 2, 1, 1 threads with 4 jobs at once
    n_centroids = {30 : 60, 33 : 30, 36 : 20, 39 : 5}
    specs, r1_by_c3len = {}, {}
    for c3len, n in n_centroids.items():
        uids = ['c%d-%d' % (c3len, i) for i in range(n)]
        specs[c3len] = {'seqfos' : [{'name' : u, 'naive_seq' : 'ACGT'} for u in uids], 'naive_fasta' : None}
        r1_by_c3len[c3len] = [(u, [u]) for u in uids]
    ran = {}
    monkeypatch.setattr(disjointgrouper.utils, 'run_cmds', lambda cmdfos, n_max_procs: ran.update(cmdfos=cmdfos, n_max_procs=n_max_procs))
    comps, _ = disjointgrouper._run_round2_tcm(specs, r1_by_c3len, 3., 0.07, str(tmp_path), 8)
    threads = {int(re.search(r'cdr3-(\d+)/', c['workdir']).group(1)) : int(re.search(r'--threads (\d+)', c['cmd_str']).group(1)) for c in ran['cmdfos']}
    assert threads == {30 : 6, 33 : 1, 36 : 1, 39 : 1}
    assert ran['n_max_procs'] == 3
    assert {c3len : len(cs) for c3len, cs in comps.items()} == n_centroids  # no pairs file, so every centroid is its own component


@pytest.mark.parametrize('sizes, budget, expected', [
    ([100], 0, ([1], 1)),
    ([100, 50], 1, ([1, 1], 1)),
    ([], 16, ([], 1)),
    ([0, 0, 0], 16, ([1, 1, 1], 3)),
    ([0], 16, ([1], 1)),
])
def test_edge_cases(sizes, budget, expected):
    assert _split_procs(sizes, budget) == expected
