import json
import pytest
import numpy as np

from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe4_pssm_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def msa_sequences(json_data):
    return json_data['msa_sequences']


@pytest.fixture(scope="module")
def invalid_msa(json_data):
    return json_data['invalid_msa']


@pytest.fixture(scope="module")
def primary_seq(json_data):
    return json_data['primary_seq']


@pytest.fixture(scope="module")
def msa_size(json_data):
    return tuple(json_data['msa_size'])


@pytest.fixture(scope="module")
def sequence_weights(json_data):
    return np.array(json_data['sequence_weights'])


@pytest.fixture(scope="module")
def num_observations(json_data):
    return np.float64(json_data['num_observations'])


@pytest.fixture(scope="module")
def bg_matrix(json_data):
    return np.array(json_data['bg_matrix'])


@pytest.fixture(scope="module")
def pssm_01(json_data):
    return np.array(json_data['pssm_01'])


@pytest.fixture(scope="module")
def pssm_02(json_data):
    return np.array(json_data['pssm_02'])


@pytest.fixture(scope="module")
def pssm_03(json_data):
    return np.array(json_data['pssm_03'])


@pytest.fixture(scope="module")
def pssm_04(json_data):
    return np.array(json_data['pssm_04'])


@pytest.fixture(scope="module")
def pssm_05(json_data):
    return np.array(json_data['pssm_05'])


@pytest.fixture(scope="module")
def pssm_06(json_data):
    return np.array(json_data['pssm_06'])


@pytest.fixture(scope="module")
def pssm_07(json_data):
    return np.array(json_data['pssm_07'])


@pytest.fixture(scope="module")
def pssm_08(json_data):
    return np.array(json_data['pssm_08'])


@pytest.fixture(scope="module")
def pssm_09(json_data):
    return np.array(json_data['pssm_09'])


@pytest.fixture(scope="module")
def pssm_10(json_data):
    return np.array(json_data['pssm_10'])


@pytest.fixture(scope="module")
def msa(msa_sequences):
    try:
        from exe4_pssm import MSA
        msa = MSA(msa_sequences)
        return msa
    except Exception:
        raise AssertionError('Error while creating MSA. MSA initialization failed.') from None


def check_error(sequences):
    try:
        from exe4_pssm import MSA
        msa = MSA(sequences)
        return False
    except TypeError:
        return True
    except Exception:
        raise AssertionError('Error while creating MSA. MSA initialization failed.') from None


@pytest.fixture(scope="module")
def pssm_raise_error(msa_sequences, invalid_msa):
    passed = all([check_error(inv_msa) for inv_msa in invalid_msa])
    return passed, 'Failed to raise TypeError for invalid MSA.'


def test_pssm_raise_error(pssm_raise_error):
    passed, assertion_msg, *_ = pssm_raise_error
    assert passed, f'Failed test pssm_raise_error(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_size(msa, msa_size):
    try:
        size = msa.get_size()
    except Exception:
        raise AssertionError('Error while testing MSA.get_size().') from None

    passed = (size == msa_size)
    return passed, 'Incorrect MSA size.'


def test_pssm_get_size(pssm_get_size):
    passed, asssertion_msg, *_ = pssm_get_size
    assert passed, f'Failed test pssm_get_size(). {asssertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_primary_sequence(msa, primary_seq):
    try:
        seq = msa.get_primary_sequence()
    except Exception:
        raise AssertionError('Error while testing MSA.get_primary_sequence().') from None

    passed = (seq == primary_seq)
    return passed, 'Incorrect primary sequence.'


def test_pssm_get_primary_sequence(pssm_get_primary_sequence):
    passed, assertion_msg, *_ = pssm_get_primary_sequence
    assert passed, f'Failed test pssm_get_primary_sequence(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_sequence_weights(msa, sequence_weights):
    try:
        weights = msa.get_sequence_weights()
    except Exception:
        raise AssertionError('Error while testing MSA.get_sequence_weights().') from None

    try:
        passed = (isinstance(weights, np.ndarray) and (weights.dtype == np.float64))
        assert passed, 'Return value not a numpy.ndarray of dtype numpy.float64.'

        passed = weights.shape == sequence_weights.shape
        assert passed, 'Incorrect dimensions of returned array.'

        passed = all(np.isclose(weights, sequence_weights, atol=1e-8, rtol=0))
        return passed, 'Incorrect sequence weights.'
    except AssertionError as msg:
        return False, msg


def test_pssm_get_sequence_weights(pssm_get_sequence_weights):
    passed, assertion_msg, *_ = pssm_get_sequence_weights
    assert passed, f'Failed test pssm_get_sequence_weights(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_number_of_observations(msa, num_observations):
    try:
        num_obs = msa.get_number_of_observations()
    except Exception:
        raise AssertionError('Error while testing MSA.get_number_of_observations().') from None

    try:
        passed = isinstance(num_obs, np.float64)
        assert passed, 'Return value not a numpy.float64.'

        passed = np.isclose(num_obs, num_observations, atol=1e-8, rtol=0)
        return passed, f'Incorrect number of independent observations.'
    except AssertionError as msg:
        return False, msg


def test_pssm_get_number_of_observations(pssm_get_number_of_observations):
    passed, assertion_msg, *_ = pssm_get_number_of_observations
    assert passed, f'Failed test pssm_get_number_of_observations(). {assertion_msg}'


def check_pssm(pssm, case, case_text):
    try:
        passed = (isinstance(pssm, np.ndarray) and (pssm.dtype == np.int64))
        assert passed, 'Return value not a numpy.ndarray of dtype numpy.int64.'

        passed = np.array_equal(pssm, case)
        return passed, f'Incorrect PSSM ({case_text}).'
    except AssertionError as msg:
        return False, msg


@pytest.fixture(scope="module")
def pssm_get_pssm_basic(msa, pssm_01):
    try:
        pssm = msa.get_pssm()
        return check_pssm(pssm, pssm_01, 'basic')
    except Exception:
        assert False, 'Error while testing MSA.get_pssm().'


def test_pssm_get_pssm_basic(pssm_get_pssm_basic):
    passed, assertion_msg, *_ = pssm_get_pssm_basic
    assert passed, f'Failed test pssm_get_pssm_basic(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_bg_matrix(msa, pssm_02, bg_matrix):
    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix)
        return check_pssm(pssm, pssm_02, 'bg_matrix')
    except Exception:
        raise AssertionError('Error while testing MSA.get_pssm(bg_matrix).') from None


def test_pssm_get_pssm_with_bg_matrix(pssm_get_pssm_with_bg_matrix):
    passed, assertion_msg, *_ = pssm_get_pssm_with_bg_matrix
    assert passed, f'Failed test pssm_get_pssm_with_bg_matrix(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_redistribute_gaps(msa, pssm_03):
    try:
        pssm = msa.get_pssm(redistribute_gaps=True)
        return check_pssm(pssm, pssm_03, 'redistribute_gaps')
    except Exception:
        raise AssertionError('Error while testing MSA.get_pssm(redistribute_gaps).') from None


def test_pssm_get_pssm_with_redistribute_gaps(pssm_get_pssm_with_redistribute_gaps):
    passed, assertion_msg, *_ = pssm_get_pssm_with_redistribute_gaps
    assert passed, f'Failed test pssm_get_pssm_with_redistribute_gaps(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_sequence_weights(msa, pssm_04):
    try:
        pssm = msa.get_pssm(use_sequence_weights=True)
        return check_pssm(pssm, pssm_04, 'use_sequence_weights')
    except Exception:
        raise AssertionError('Error while testing get_pssm(use_sequence_weights).') from None


def test_pssm_get_pssm_with_sequence_weights(pssm_get_pssm_with_sequence_weights):
    passed, assertion_msg, *_ = pssm_get_pssm_with_sequence_weights
    assert passed, f'Failed test pssm_get_pssm_with_sequence_weights(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_bg_matrix_and_redistribute_gaps(msa, pssm_05, bg_matrix):
    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix, redistribute_gaps=True)
        return check_pssm(pssm, pssm_05, 'bg_matrix, redistribute_gaps')
    except Exception:
        raise AssertionError('Error while testing MSA.get_pssm'
                             '(bg_matrix, redistribute_gaps).') from None


def test_pssm_get_pssm_with_bg_matrix_and_redistribute_gaps(pssm_get_pssm_with_bg_matrix_and_redistribute_gaps):
    passed, assertion_msg, *_ = pssm_get_pssm_with_bg_matrix_and_redistribute_gaps
    assert passed, f'Failed test pssm_get_pssm_with_bg_matrix_and_redistribute_gaps(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_bg_matrix_and_sequence_weights(msa, pssm_06, bg_matrix):
    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix, use_sequence_weights=True)
        return check_pssm(pssm, pssm_06, 'bg_matrix, use_sequence_weights')
    except Exception:
        raise AssertionError('Error while testing get_pssm'
                             '(bg_matrix, use_sequence_weights).') from None


def test_pssm_get_pssm_with_bg_matrix_and_sequence_weights(pssm_get_pssm_with_bg_matrix_and_sequence_weights):
    passed, assertion_msg, *_ = pssm_get_pssm_with_bg_matrix_and_sequence_weights
    assert passed, f'Failed test pssm_get_pssm_with_bg_matrix_and_sequence_weights(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_pseudocounts(msa, pssm_07, bg_matrix):
    try:
        pssm = msa.get_pssm(add_pseudocounts=True)
        return check_pssm(pssm, pssm_07, 'add_pseudocounts')
    except Exception:
        raise AssertionError('Error while testing get_pssm(add_pseudocounts).') from None


def test_pssm_get_pssm_with_pseudocounts(pssm_get_pssm_with_pseudocounts):
    passed, assertion_msg, *_ = pssm_get_pssm_with_pseudocounts
    assert passed, f'Failed test pssm_get_pssm_with_pseudocounts(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_bg_matrix_and_redistribute_gaps_and_pseudocounts(msa, pssm_08, bg_matrix):
    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix,
                            redistribute_gaps=True,
                            add_pseudocounts=True)
        return check_pssm(pssm, pssm_08, 'bg_matrix, redistribute_gaps, add_pseudocounts')
    except Exception:
        raise AssertionError('Error while testing get_pssm'
                             '(bg_matrix, redistribute_gaps, add_pseudocounts).') from None


def test_pssm_get_pssm_with_bg_matrix_and_redistribute_gaps_and_pseudocounts(
        pssm_get_pssm_with_bg_matrix_and_redistribute_gaps_and_pseudocounts):
    passed, assertion_msg, *_ = pssm_get_pssm_with_bg_matrix_and_redistribute_gaps_and_pseudocounts
    assert passed, f'Failed test pssm_get_pssm_with_bg_matrix_and_redistribute_gaps_and_pseudocounts(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_pseudocounts(msa, pssm_09, bg_matrix):
    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix,
                            use_sequence_weights=True,
                            add_pseudocounts=True)
        return check_pssm(pssm, pssm_09, 'bg_matrix, use_sequence_weights, add_pseudocounts')
    except Exception:
        raise AssertionError('Error while testing get_pssm'
                             '(bg_matrix, use_sequence_weights, add_pseudocounts).') from None


def test_pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_pseudocounts(
        pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_pseudocounts):
    passed, assertion_msg, *_ = pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_pseudocounts
    assert passed, f'Failed test pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_pseudocounts(). {assertion_msg}'


@pytest.fixture(scope="module")
def pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_redistribute_gaps_and_pseudocounts(msa, pssm_10, bg_matrix):
    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix,
                            redistribute_gaps=True,
                            use_sequence_weights=True,
                            add_pseudocounts=True)
        return check_pssm(pssm, pssm_10, 'bg_matrix, redistribute_gaps, use_sequence_weights, add_pseudocounts')
    except Exception:
        raise AssertionError('Error while testing get_pssm'
                             '(bg_matrix, redistribute_gaps, '
                             'use_sequence_weights, add_pseudocounts).') from None


def test_pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_redistribute_gaps_and_pseudocounts(
        pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_redistribute_gaps_and_pseudocounts):
    passed, assertion_msg, *_ = pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_redistribute_gaps_and_pseudocounts
    assert passed, f'Failed test pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_redistribute_gaps_and_pseudocounts(). {assertion_msg}'
