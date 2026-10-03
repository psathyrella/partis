import random

import pytest

from partis.disjointgrouper import _split_procs


def test_one_big_group_gets_most_threads():
    threads, n_jobs = _split_procs([2100000] + [50000] * 19, 16)
    assert threads[0] == 11
    assert n_jobs == 6
    assert sum(threads[:n_jobs]) == 16


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
        assert n_jobs >= 1
        assert sum(threads[:n_jobs]) <= budget


@pytest.mark.parametrize('sizes, budget, expected', [
    ([100], 0, ([1], 1)),
    ([100, 50], 1, ([1, 1], 1)),
    ([], 16, ([], 1)),
    ([0, 0, 0], 16, ([1, 1, 1], 3)),
    ([0], 16, ([1], 1)),
])
def test_edge_cases(sizes, budget, expected):
    assert _split_procs(sizes, budget) == expected
