import numpy as np

from matchmaker.utils import (
    check_slash,
    check_strictly_positive,
    check_underscore,
    difference,
    eval_flux_given_spectral_index,
    intersection,
    load,
    save,
    share_items,
    union,
)


def test_check_underscore():
    assert check_underscore('prefix') == 'prefix_'
    assert check_underscore('prefix_') == 'prefix_'
    assert check_underscore('') == ''


def test_check_slash():
    assert check_slash('folder') == 'folder/'
    assert check_slash('folder/') == 'folder/'
    assert check_slash('') == ''


def test_set_helpers():
    a, b = [1, 2, 3], [2, 3, 4]
    assert share_items(a, b)
    assert not share_items([1], [2])
    assert intersection(a, b).tolist() == [2, 3]
    assert difference(a, b).tolist() == [1]
    assert union(a, b).tolist() == [1, 2, 3, 4]


def test_check_strictly_positive():
    assert check_strictly_positive(5) == 5
    assert check_strictly_positive(0) == 10**-4
    assert check_strictly_positive(-3) == 10**-4


def test_save_load_roundtrip(tmp_path):
    folder = f'{tmp_path}/'
    data = {'a': np.arange(3), 'b': 'text'}
    save('myvar', data, folder=folder)
    loaded = load('myvar', folder=folder)
    assert loaded['b'] == 'text'
    assert loaded['a'].tolist() == [0, 1, 2]


def test_load_missing_returns_none(tmp_path):
    assert load('nope', folder=f'{tmp_path}/') is None


def test_eval_flux_given_spectral_index():
    # alpha=0: flux unchanged; alpha=-1: flux halves when frequency doubles
    assert eval_flux_given_spectral_index(1.0, 2.0, 10.0, 0.0) == 10.0
    assert eval_flux_given_spectral_index(1.0, 2.0, 10.0, -1.0) == 5.0
