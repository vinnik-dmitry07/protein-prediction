import json
import pytest

from math import ceil
from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe1_orffinder_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def test_cases(json_data):
    cases = [case for case in json_data['sequences'].values()]
    return cases


@pytest.fixture(scope="module")
def sequences(test_cases):
    return [case.get('sequence', None) for case in test_cases]


@pytest.fixture(scope="module")
def invalid_genome(json_data):
    return json_data['invalid_genome']


@pytest.fixture(scope="module")
def min_num_aa(json_data):
    return json_data['min_num_aa']


@pytest.fixture(scope="module")
def orf_list(test_cases):
    cases = [case.get('orf_list', None) for case in test_cases]
    orfs = [tuple(e) for case in cases for e in case]
    return set(orfs)


@pytest.fixture(scope="module")
def orffinder_get_orfs(sequences, min_num_aa):
    orfs_total = []
    for sequence in sequences:
        try:
            from exe1_orffinder import get_orfs
            orfs_case = get_orfs(sequence, min_num_aa)
        except Exception:
            raise AssertionError('Error in get_orfs().') from None

        try:
            passed = isinstance(orfs_case, list)
            assert passed, 'Return type is not a list.'

            passed = all([isinstance(t, tuple) for t in orfs_case])
            assert passed, 'List does not contain only tuples.'

            passed = all([(len(t) == 4) for t in orfs_case])
            assert passed, 'Not all tuples contain exactly 4 elements.'

            for p1, p2, seq, flag in orfs_case:
                passed = all([isinstance(p1, int),
                              isinstance(p2, int),
                              isinstance(seq, str),
                              isinstance(flag, bool)])
                assert passed, 'Tuple malformed. Expected: (int, int, str, bool).'

            orfs_total.extend(orfs_case)
        except AssertionError as msg:
            return False, msg, None
    return True, '', orfs_total


def test_orffinder_get_orfs(orffinder_get_orfs):
    passed, assertion_msg, *_ = orffinder_get_orfs
    assert passed, f'Failed test orffinder_get_orfs(). {assertion_msg}'


def test_orffinder_raise_error(orffinder_get_orfs, invalid_genome, min_num_aa):
    assert orffinder_get_orfs[0], 'Failed test orffinder_get_orfs().'
    try:
        from exe1_orffinder import get_orfs
        get_orfs(invalid_genome, min_num_aa)
    except TypeError:
        pass
    except Exception:
        raise AssertionError('Failed to raise TypeError for invalid genome.') from None


@pytest.fixture(scope="module")
def orffinder_valid_orfs(orf_list, orffinder_get_orfs):
    try:
        passed, assertion_msg, student_orfs = orffinder_get_orfs
        assert passed, f'Failed test orffinder_get_orfs(). {assertion_msg}'

        num_student = len(student_orfs)
        student_set = set(student_orfs)

        passed = (num_student > 0)
        assert passed, 'List is empty.'

        passed = (num_student == len(student_set))
        assert passed, 'List contains duplicates.'

        valid_orfs = len((orf_list & student_set))
        total_orfs = len((orf_list | student_set))

        passed = (valid_orfs > 0)
        assert passed, 'List does not contain any valid orfs.'

        return passed, '', valid_orfs, total_orfs
    except AssertionError as msg:
        return False, msg, 0, 999999


def test_orffinder_valid_orfs(orffinder_valid_orfs):
    passed, assertion_msg, *_ = orffinder_valid_orfs
    assert passed, f'Failed test_orffinder_valid_orfs(). {assertion_msg}'


# The following tests will always fail if the ORF list contains duplicates!
def test_orffinder_valid_orfs_25(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs >= ceil(0.25 * total_orfs))
    assert passed, 'Failed test_orffinder_valid_orfs_25(). Less than 25% valid ORFs.'


def test_orffinder_valid_orfs_50(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs >= ceil(0.50 * total_orfs))
    assert passed, 'Failed test_orffinder_valid_orfs_50(). Less than 50% valid ORFs.'


def test_orffinder_valid_orfs_75(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs >= ceil(0.75 * total_orfs))
    assert passed, 'Failed test_orffinder_valid_orfs_75(). Less than 75% valid ORFs.'


def test_orffinder_valid_orfs_95(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs >= ceil(0.95 * total_orfs))
    assert passed, 'Failed test_orffinder_valid_orfs_95(). Less than 95% valid ORFs.'


def test_orffinder_valid_orfs_100_x1(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs == total_orfs)
    assert passed, 'Failed test_orffinder_valid_orfs_100(). Less than 100% valid ORFs.'


def test_orffinder_valid_orfs_100_x2(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs == total_orfs)
    assert passed, 'Failed test_orffinder_valid_orfs_100(). Less than 100% valid ORFs.'


def test_orffinder_valid_orfs_100_x3(orffinder_valid_orfs):
    _, _, valid_orfs, total_orfs, *_ = orffinder_valid_orfs
    passed = (valid_orfs == total_orfs)
    assert passed, 'Failed test_orffinder_valid_orfs_100(). Less than 100% valid ORFs.'
