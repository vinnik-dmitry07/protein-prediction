import json
import pytest

from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe0_motifs_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def test_cases(json_data):
    cases = [json_data[i] for i in json_data]
    return cases


@pytest.fixture(scope="module")
def test_subsearch(test_cases):
    for cases in test_cases:
        try:
            from exe0_motifs import find_similar_motifs
            student_motifs = find_similar_motifs(*cases['input'])
        except Exception:
            raise AssertionError('Error in exe0_motifs.find_similar_motifs().') from None

        try:
            passed = isinstance(student_motifs, list)

            assert passed, 'Return type is not a list.'

            passed = (len(student_motifs) == len(cases['output']))

            assert passed, 'Number of motifs not correct.'

            passed = (sorted(student_motifs) == sorted(cases['output']))

            assert passed, 'Motifs incorrect.'
        except AssertionError as msg:
            return False, msg

    return True, ''


def test_subsearch_x1(test_subsearch):
    passed, assertion_msg, *_ = test_subsearch
    assert passed, f'Failed test_subsearch(). {assertion_msg}'


def test_subsearch_x2(test_subsearch):
    passed, assertion_msg, *_ = test_subsearch
    assert passed, f'Failed test_subsearch(). {assertion_msg}'


def test_subsearch_x3(test_subsearch):
    passed, assertion_msg, *_ = test_subsearch
    assert passed, f'Failed test_subsearch(). {assertion_msg}'


def test_subsearch_x4(test_subsearch):
    passed, assertion_msg, *_ = test_subsearch
    assert passed, f'Failed test_subsearch(). {assertion_msg}'


def test_subsearch_x5(test_subsearch):
    passed, assertion_msg, *_ = test_subsearch
    assert passed, f'Failed test_subsearch(). {assertion_msg}'
