import json
import pytest
import numpy as np

from pathlib import Path
from .matrices import MATRICES
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe3_global_alignment_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def small(json_data):
    return json_data['small']


@pytest.fixture(scope="module")
def large(json_data):
    return json_data['large']


@pytest.fixture(scope="module")
def matching(json_data):
    return json_data['matching']


@pytest.fixture(scope="module")
def mismatching(json_data):
    return json_data['mismatching']


@pytest.fixture(scope="module")
def indel(json_data):
    return json_data['indel']


@pytest.fixture(scope="module")
def all_changes(json_data):
    return json_data['all_changes']


def create_GA(string_1, string_2, gap_penalty, matrix):
    try:
        from exe3_global_alignment import GlobalAlignment
        ga = GlobalAlignment(string_1, string_2, gap_penalty, matrix)
        return ga is not None, 'GlobalAlignment initialization failed.', ga
    except Exception:
        return False, 'Error while creating GlobalAlignment.', None


@pytest.fixture(scope="module")
def small_ga(small):
    passed, assertion_msg, small_ga = create_GA(*small['strings'],
                                                small['gap_penalty'],
                                                MATRICES[small['matrix']])
    assert passed, assertion_msg
    return small_ga


@pytest.fixture(scope="module")
def get_best_score_on_small_strings(small_ga, small):
    try:
        score = small_ga.get_best_score()
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_best_score().') from None

    try:
        passed = (small['best_score'] == score)
        return passed, 'Incorrect best score (small strings).'
    except AssertionError as msg:
        return False, msg


def test_get_best_score_on_small_strings(get_best_score_on_small_strings):
    passed, assertion_msg, *_ = get_best_score_on_small_strings
    assert passed, f'Failed test get_best_score_on_small_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def large_ga(large):
    passed, assertion_msg, large_ga = create_GA(*large['strings'],
                                                large['gap_penalty'],
                                                MATRICES[large['matrix']])
    assert passed, assertion_msg
    return large_ga


@pytest.fixture(scope="module")
def get_best_score_on_large_strings(large_ga, large):
    try:
        score = large_ga.get_best_score()
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_best_score().') from None

    try:
        passed = (large['best_score'] == score)
        return passed, 'Incorrect best score (large strings).'
    except AssertionError as msg:
        return False, msg


def test_get_best_score_on_large_strings(get_best_score_on_large_strings):
    passed, assertion_msg, *_ = get_best_score_on_large_strings
    assert passed, f'Failed test get_best_score_on_large_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_number_of_alignments_on_small_strings(small_ga, small):
    try:
        num_alignments = small_ga.get_number_of_alignments()
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_number_of_alignments().') from None

    try:
        passed = isinstance(num_alignments, int)
        assert passed, 'Return type is not int.'

        passed = (small['number_of_alignments'] == num_alignments)
        return passed, 'Incorrect number of alignments (small strings).'
    except AssertionError as msg:
        return False, msg


def test_get_number_of_alignments_on_small_strings(get_number_of_alignments_on_small_strings):
    passed, assertion_msg, *_ = get_number_of_alignments_on_small_strings
    assert passed, f'Failed test get_number_of_alignments_on_small_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_number_of_alignments_on_large_strings(large_ga, large):
    try:
        num_alignments = large_ga.get_number_of_alignments()
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_number_of_alignments().') from None

    try:
        passed = isinstance(num_alignments, int)
        assert passed, 'Return type is not int.'

        passed = (large['number_of_alignments'] == num_alignments)
        return passed, 'Incorrect number of alignments (large strings).'
    except AssertionError as msg:
        return False, msg


def test_get_number_of_alignments_on_large_strings(get_number_of_alignments_on_large_strings):
    passed, assertion_msg, *_ = get_number_of_alignments_on_large_strings
    assert passed, f'Failed test get_number_of_alignments_on_large_strings(). {assertion_msg}'


def check_alignments(alignments, case, case_text):
    try:
        passed = isinstance(alignments, list)
        assert passed, 'Return type is not a list.'

        passed = all([isinstance(t, tuple) for t in alignments])
        assert passed, 'List does not contain only tuples.'

        passed = all([(len(t) == 2) for t in alignments])
        assert passed, 'Not all tuples contain exactly 2 elements.'

        for s1, s2 in alignments:
            passed = all([isinstance(s1, str),
                          isinstance(s2, str)])
            assert passed, 'Tuple malformed. Expected: (str, str).'

        passed = ({tuple(x) for x in case['alignments']} == set(alignments))
        return passed, f'Incorrect alignments ({case_text} strings).'
    except AssertionError as msg:
        return False, msg


@pytest.fixture(scope="module")
def get_alignments_on_small_strings(small_ga, small):
    try:
        alignments = small_ga.get_alignments()
        return check_alignments(alignments, small, 'small')
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_alignments().') from None


def test_get_alignments_on_small_strings(get_alignments_on_small_strings):
    passed, assertion_msg, *_ = get_alignments_on_small_strings
    assert passed, f'Failed test get_alignments_on_small_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_alignments_on_large_strings(large_ga, large):
    try:
        alignments = large_ga.get_alignments()
        return check_alignments(alignments, large, 'large')
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_alignments().') from None


def test_get_alignments_on_large_strings(get_alignments_on_large_strings):
    passed, assertion_msg, *_ = get_alignments_on_large_strings
    assert passed, f'Failed test get_alignments_on_large_strings(). {assertion_msg}'


def check_score_matrix(score_matrix, case, case_text):
    passed = np.array_equal(np.array(score_matrix),
                            np.array(case['score_matrix']))
    return passed, f'Incorrect score matrix ({case_text} strings).'


@pytest.fixture(scope="module")
def score_matrix_on_matching_strings(matching):
    try:
        _, _, ga = create_GA(*matching['strings'],
                             matching['gap_penalty'],
                             MATRICES[matching['matrix']])
        score_matrix = ga.get_score_matrix()
        return check_score_matrix(score_matrix, matching, 'matching')
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_score_matrix().') from None


def test_score_matrix_on_matching_strings(score_matrix_on_matching_strings):
    passed, assertion_msg, *_ = score_matrix_on_matching_strings
    assert passed, f'Failed test score_matrix_on_matching_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def score_matrix_on_mismatching_strings(mismatching):
    try:
        _, _, ga = create_GA(*mismatching['strings'],
                             mismatching['gap_penalty'],
                             MATRICES[mismatching['matrix']])
        score_matrix = ga.get_score_matrix()
        return check_score_matrix(score_matrix, mismatching, 'mismatching')
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_score_matrix().') from None


def test_score_matrix_on_mismatching_strings(score_matrix_on_mismatching_strings):
    passed, assertion_msg, *_ = score_matrix_on_mismatching_strings
    assert passed, f'Failed test score_matrix_on_mismatching_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def score_matrix_on_indel_strings(indel):
    try:
        _, _, ga = create_GA(*indel['strings'],
                             indel['gap_penalty'],
                             MATRICES[indel['matrix']])
        score_matrix = ga.get_score_matrix()
        return check_score_matrix(score_matrix, indel, 'indel')
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_score_matrix().') from None


def test_score_matrix_on_indel_strings(score_matrix_on_indel_strings):
    passed, assertion_msg, *_ = score_matrix_on_indel_strings
    assert passed, f'Failed test score_matrix_on_indel_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def score_matrix_on_all_changes_strings(all_changes):
    try:
        _, _, ga = create_GA(*all_changes['strings'],
                             all_changes['gap_penalty'],
                             MATRICES[all_changes['matrix']])
        score_matrix = ga.get_score_matrix()
        return check_score_matrix(score_matrix, all_changes, 'all_changes')
    except Exception:
        raise AssertionError('Error in GlobalAlignment.get_score_matrix().') from None


def test_score_matrix_on_all_changes_strings(score_matrix_on_all_changes_strings):
    passed, assertion_msg, *_ = score_matrix_on_all_changes_strings
    assert passed, f'Failed test score_matrix_on_all_changes_strings(). {assertion_msg}'
