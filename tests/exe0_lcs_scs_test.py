import json
import pytest

from math import ceil
from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe0_lcs_scs_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def sequences(json_data):
    return json_data['sequence_1'], json_data['sequence_2']


@pytest.fixture(scope="module")
def lcs_length(json_data):
    return json_data['lcs']['lcs_length']


@pytest.fixture(scope="module")
def lcs_list(json_data):
    return json_data['lcs']['lcs_list']


@pytest.fixture(scope="module")
def lcs_count(json_data):
    return json_data['lcs']['lcs_count']


@pytest.fixture(scope="module")
def scs_length(json_data):
    return json_data['scs']['scs_length']


@pytest.fixture(scope="module")
def scs_list(json_data):
    return json_data['scs']['scs_list']


@pytest.fixture(scope="module")
def scs_count(json_data):
    return json_data['scs']['scs_count']


@pytest.fixture(scope="module")
def get_lcs_length(sequences, lcs_length):
    try:
        from exe0_lcs_scs import longest_common_subseq_length
        student_lcs_length = longest_common_subseq_length(sequences[0], sequences[1])
    except Exception:
        raise AssertionError('Error in longest_common_subseq_length().') from None

    try:
        passed = isinstance(student_lcs_length, int)
        assert passed, 'Return type is not an int.'

        passed = student_lcs_length == lcs_length
        return passed, 'Length of longest common subsequence incorrect.'
    except AssertionError as msg:
        return False, msg


def test_get_lcs_length(get_lcs_length):
    passed, assertion_msg, *_ = get_lcs_length
    assert passed, f'Failed test get_lcs_length(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_lcs_count(sequences, lcs_count):
    try:
        from exe0_lcs_scs import longest_common_subseq_count
        student_lcs_count = longest_common_subseq_count(sequences[0], sequences[1])
    except Exception:
        raise AssertionError('Error in longest_common_subseq_count().') from None

    try:
        passed = isinstance(student_lcs_count, int)
        assert passed, 'Return type is not an int.'

        passed = student_lcs_count == lcs_count
        return passed, 'Number of longest common subsequences incorrect.'
    except AssertionError as msg:
        return False, msg


def test_get_lcs_count(get_lcs_count):
    passed, assertion_msg, *_ = get_lcs_count
    assert passed, f'Failed test get_lcs_count(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_lcs_list(sequences, lcs_list):
    try:
        from exe0_lcs_scs import longest_common_subseq
        student_lcs = longest_common_subseq(sequences[0], sequences[1])
    except Exception:
        raise AssertionError('Error in longest_common_subseq().') from None

    try:
        passed = isinstance(student_lcs, list) or isinstance(student_lcs, set)
        assert passed, 'Return type is not a list or set.'

        passed = (len(student_lcs) > 0)
        assert passed, 'List is empty.'

        passed = all([isinstance(seq, str) for seq in student_lcs])
        assert passed, 'List does not contain only strings.'

        passed = len(set(map(len, student_lcs))) == 1
        assert passed, 'Not all sequences have the same length.'

        passed = (len(student_lcs) == len(set(student_lcs)))
        assert passed, 'List contains duplicates.'

        valid_lcs = len((set(lcs_list) & set(student_lcs)))
        total_lcs = len((set(lcs_list) | set(student_lcs)))

        passed = (valid_lcs > 0)
        assert passed, 'List or set does not contain any valid longest common subsequences.'

        return passed, '', valid_lcs, total_lcs
    except AssertionError as msg:
        return False, msg, 0, 999999


def test_get_lcs_list(get_lcs_list):
    passed, assertion_msg, *_ = get_lcs_list
    assert passed, f'Failed test get_lcs_list(). {assertion_msg}'


def test_lcs_50(get_lcs_list):
    passed, _, valid_lcs, total_lcs, *_ = get_lcs_list
    passed = (valid_lcs >= ceil(0.50 * total_lcs))
    assert passed, 'Less than 50% valid longest common subsequences.'


@pytest.fixture(scope="module")
def test_lcs_100(get_lcs_list):
    _, _, valid_lcs, total_lcs, *_ = get_lcs_list
    passed = (valid_lcs == total_lcs)
    return passed, 'Less than 100% valid longest common subsequence.'


def test_lcs_100_x1(test_lcs_100):
    passed, assertion_msg, *_ = test_lcs_100
    assert passed, f'Failed test test_lcs_100(). {assertion_msg}'


def test_lcs_100_x2(test_lcs_100):
    passed, assertion_msg, *_ = test_lcs_100
    assert passed, f'Failed test test_lcs_100(). {assertion_msg}'

########################################################################################################################


@pytest.fixture(scope="module")
def get_scs_length(sequences, scs_length):
    try:
        from exe0_lcs_scs import shortest_common_supersequence_length
        student_scs_length = shortest_common_supersequence_length(sequences[0], sequences[1])
    except Exception:
        raise AssertionError('Error in shortest_common_supersequence_length().') from None

    try:
        passed = isinstance(student_scs_length, int)
        assert passed, 'Return type is not an int.'

        passed = student_scs_length == scs_length
        return passed, 'Length of shortest common supersequences incorrect.'
    except AssertionError as msg:
        return False, msg


def test_get_scs_length(get_scs_length):
    passed, assertion_msg, *_ = get_scs_length
    assert passed, f'Failed test get_scs_length(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_scs_count(sequences, scs_count):
    try:
        from exe0_lcs_scs import shortest_common_supersequence_count
        student_scs_count = shortest_common_supersequence_count(sequences[0], sequences[1])
    except Exception:
        raise AssertionError('Error in shortest_common_supersequence_count().') from None

    try:
        passed = isinstance(student_scs_count, int)
        assert passed, 'Return type is not an int.'

        passed = student_scs_count == scs_count
        return passed, 'Number of shortest common supersequences incorrect.'
    except AssertionError as msg:
        return False, msg


def test_get_scs_count(get_scs_count):
    passed, assertion_msg, *_ = get_scs_count
    assert passed, f'Failed test get_scs_count(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_scs_list(sequences, scs_list):
    try:
        from exe0_lcs_scs import shortest_common_supersequence
        student_scs = shortest_common_supersequence(sequences[0], sequences[1])
    except Exception:
        raise AssertionError('Error in shortest_common_supersequence().') from None

    try:
        passed = isinstance(student_scs, list) or isinstance(student_scs, set)
        assert passed, 'Return type is not a list or set.'

        passed = (len(student_scs) > 0)
        assert passed, 'List is empty.'

        passed = all([isinstance(seq, str) for seq in student_scs])
        assert passed, 'List does not contain only strings.'

        passed = len(set(map(len, student_scs))) == 1
        assert passed, 'Not all sequences have the same length.'

        passed = (len(student_scs) == len(set(student_scs)))
        assert passed, 'List contains duplicates.'

        valid_scs = len((set(scs_list) & set(student_scs)))
        total_scs = len((set(scs_list) | set(student_scs)))

        passed = (valid_scs > 0)
        assert passed, 'List or set does not contain any valid shortest common supersequence.'

        return passed, '', valid_scs, total_scs
    except AssertionError as msg:
        return False, msg, 0, 999999


def test_get_scs_list(get_scs_list):
    passed, assertion_msg, *_ = get_scs_list
    assert passed, f'Failed test get_scs_list(). {assertion_msg}'


def test_scs_50(get_scs_list):
    passed, _, valid_scs, total_scs = get_scs_list
    passed = (valid_scs >= ceil(0.50 * total_scs))
    assert passed, 'Failed test_scs_50(). Less than 50% valid shortest common supersequence.'


@pytest.fixture(scope="module")
def test_scs_100(get_scs_list):
    _, _, valid_scs, total_scs, *_ = get_scs_list
    passed = (valid_scs == total_scs)
    return passed, 'Less than 100% valid shortest common supersequence.'


def test_scs_100_x1(test_scs_100):
    passed, assertion_msg, *_ = test_scs_100
    assert passed, f'Failed test_scs_100(). {assertion_msg}'

    
def test_scs_100_x2(test_scs_100):
    passed, assertion_msg, *_ = test_scs_100
    assert passed, f'Failed test_scs_100(). {assertion_msg}'
