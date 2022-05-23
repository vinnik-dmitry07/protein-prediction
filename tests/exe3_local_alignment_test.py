import json
import pytest
import numpy as np

from pathlib import Path
from .matrices import MATRICES
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe3_local_alignment_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def null(json_data):
    return json_data['null']


@pytest.fixture(scope="module")
def short(json_data):
    return json_data['short']


def create_LA(string_1, string_2, gap_penalty, matrix):
    try:
        from exe3_local_alignment import LocalAlignment
        la = LocalAlignment(string_1, string_2, gap_penalty, matrix)
        return la is not None, 'LocalAlignment initialization failed.', la
    except Exception:
        return False, 'Error while creating LocalAlignment.', None


@pytest.fixture(scope="module")
def null_la(null):
    passed, assertion_msg, null_la = create_LA(*null['strings'],
                                               null['gap_penalty'],
                                               MATRICES[null['matrix']])
    assert passed, assertion_msg
    return null_la


@pytest.fixture(scope="module")
def short_la(short):
    passed, assertion_msg, short_la = create_LA(*short['strings'],
                                                short['gap_penalty'],
                                                MATRICES[short['matrix']])
    assert passed, assertion_msg
    return short_la


@pytest.fixture(scope="module")
def has_alignment(null_la, short_la):
    try:
        null_has_alignment = null_la.has_alignment()
        short_has_alignment = short_la.has_alignment()
    except Exception:
        raise AssertionError('Error in LocalAlignment.has_alignment().') from None

    try:
        passed = (not null_has_alignment)
        assert passed, 'Incorrect value for has_alignment(); null strings.'

        passed = (short_has_alignment)
        assert passed, 'Incorrect value for has_alignment(); small strings.'

        return True, ''
    except AssertionError as msg:
        return False, msg


def test_has_alignment(has_alignment):
    passed, assertion_msg, *_ = has_alignment
    assert passed, f'Failed test has_alignment(). {assertion_msg}'


def check_alignment(alignment, case, case_text):
    try:
        passed = isinstance(alignment, tuple)
        assert passed, 'Return type is not a tuple.'

        passed = len(alignment) == 2
        assert passed, f'Tuple contains {len(alignment)} elements. Expected: 2'

        passed = all([isinstance(s, str) for s in alignment])
        assert passed, 'Tuple does not contain only strings.'

        passed = (tuple(case['alignment']) == alignment)
        return passed, f'Incorrect alignment ({case_text} strings).'
    except AssertionError as msg:
        return False, msg


@pytest.fixture(scope="module")
def get_alignment_on_small_strings(short_la, short):
    try:
        alignment = short_la.get_alignment()
        return check_alignment(alignment, short, 'short')
    except Exception:
        raise AssertionError('Error in LocalAlignment.get_alignment().') from None


def test_get_alignment_on_small_strings(get_alignment_on_small_strings):
    passed, assertion_msg, *_ = get_alignment_on_small_strings
    assert passed, f'Failed test get_alignment_on_small_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def get_alignment_on_null_strings(null_la, null):
    try:
        alignment = null_la.get_alignment()
        return check_alignment(alignment, null, 'null')
    except Exception:
        raise AssertionError('Error in LocalAlignment.get_alignment().') from None


def test_get_alignment_on_null_strings(get_alignment_on_null_strings):
    passed, assertion_msg, *_ = get_alignment_on_null_strings
    assert passed, f'Failed test get_alignment_on_null_strings(). {assertion_msg}'


@pytest.fixture(scope="module")
def is_residue_aligned_on_first_string(short_la, short):
    [s_1, p_1, res_1], [s_2, p_2, res_2] = short["residue_aligned_on_first"]

    try:
        is_aligned_1 = short_la.is_residue_aligned(s_1, p_1)
        is_aligned_2 = short_la.is_residue_aligned(s_2, p_2)
    except Exception:
        raise AssertionError('Error in LocalAlignment.is_residue_aligned().') from None

    passed = ((res_1 == is_aligned_1) and (res_2 == is_aligned_2))
    return passed, 'Incorrect value for is_residue_aligned(); first string.'


def test_is_residue_aligned_on_first_string(is_residue_aligned_on_first_string):
    passed, assertion_msg, *_ = is_residue_aligned_on_first_string
    assert passed, f'Failed test is_residue_aligned_on_first_string(). {assertion_msg}'


@pytest.fixture(scope="module")
def is_residue_aligned_on_second_string(short_la, short):
    [s_1, p_1, res_1], [s_2, p_2, res_2] = short["residue_aligned_on_second"]

    try:
        is_aligned_1 = short_la.is_residue_aligned(s_1, p_1)
        is_aligned_2 = short_la.is_residue_aligned(s_2, p_2)
    except Exception:
        raise AssertionError('Error in LocalAlignment.is_residue_aligned().') from None

    passed = ((res_1 == is_aligned_1) and (res_2 == is_aligned_2))
    return passed, 'Incorrect value for is_residue_aligned(); second string.'


def test_is_residue_aligned_on_second_string(is_residue_aligned_on_second_string):
    passed, assertion_msg, *_ = is_residue_aligned_on_second_string
    assert passed, f'Failed test is_residue_aligned_on_second_string(). {assertion_msg}'
