"""
Tests salt intrusion lenght computation with gap detection and gap threshold.

Tuomas Karna 2016-03-03
"""
from crane.data.generateSaltIntrusion import compute_intrusion_length
import numpy as np
import pytest

def assert_sil(x, salt, salt_threshold, min_gap_length, correct_sil):
    for i in range(len(x)):
        print('{:4d} {:12.3f} {:12.3f}'.format(i, x[i], salt[i]))

    x_sil = compute_intrusion_length(x, salt, salt_threshold, min_gap_length=min_gap_length, verbose=True)
    assert x_sil == correct_sil, 'wrong sil: {:} != {:}'.format(x_sil, correct_sil)

MAX_SALT = 20.0

@pytest.fixture(scope='function')
def x_arr(request):
    n = 101
    x = np.linspace(0, 100e3, n)
    return x


@pytest.fixture(scope='function')
def salt_arr(x_arr):
    s = np.zeros_like(x_arr)
    return s
    

def test_no_gaps_case(x_arr, salt_arr):
    salt_arr[:30] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 50.0e3, 29500.0)

def test_interpolation(x_arr, salt_arr):
    salt_arr[:30] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/4, 50.0e3, 29750.0)

def test_salt_patch_in_end_ingored(x_arr, salt_arr):
    salt_arr[:30] = MAX_SALT
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 50.0e3, 29500.0)

def test_salt_patch_in_end_long_gap(x_arr, salt_arr):
    salt_arr[:30] = MAX_SALT
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 500.0e3, 100.0e3)

def test_salt_gap_at_mouth(x_arr, salt_arr):
    salt_arr[10:30] = MAX_SALT
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 8.0e3, 0.0)

def test_salt_gap_at_mouth_ignored(x_arr, salt_arr):
    salt_arr[10:30] = MAX_SALT
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 11.0e3, 29500.0)

def test_many_gaps(x_arr, salt_arr):
    salt_arr[:5] = MAX_SALT
    salt_arr[10::12] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 4.0e3, 4500.0)

def test_many_gaps_bigger_gap(x_arr, salt_arr):
    salt_arr[:5] = MAX_SALT
    salt_arr[10::12] = MAX_SALT
    salt_arr[-20:] = 0.0
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 6.0e3, 10500.0)

def test_many_gaps_bigger_gap(x_arr, salt_arr):
    salt_arr[:5] = MAX_SALT
    salt_arr[10::12] = MAX_SALT
    salt_arr[-20:] = 0.0
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 12.0e3, 70500.0)

def test_many_gaps_bigger_gap2(x_arr, salt_arr):
    salt_arr[:5] = MAX_SALT
    salt_arr[10::12] = MAX_SALT
    salt_arr[-20:] = 0.0
    salt_arr[-2:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 29.0e3, 100.0e3)

def test_zero_gap(x_arr, salt_arr):
    salt_arr[:5] = MAX_SALT
    salt_arr[10::12] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 0.0, 4500.0)

def test_full_salt(x_arr, salt_arr):
    salt_arr[:] = MAX_SALT
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 5.0e3, 100.0e3)

def test_no_salt(x_arr, salt_arr):
    assert_sil(x_arr, salt_arr, MAX_SALT/2, 5.0e3, 0.0)
