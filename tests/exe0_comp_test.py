import json
import pytest

from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe0_comp_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def test_cases(json_data):
    cases = [json_data[i] for i in json_data]
    return cases


@pytest.fixture(scope="module")
def test_complementary(test_cases):
    passed = True

    for cases in test_cases:
        try:
            from exe0_comp import complementary
            student_complement = complementary(cases['input']).upper()
        except Exception:
            raise AssertionError('Error in exe0_comp.complementary().') from None

        try:
            passed = isinstance(student_complement, str)
            assert passed, 'Return type is not a string.'

            passed = (student_complement == cases['output'])
            assert passed, 'Complement is wrong.'
        except AssertionError as msg:
            return False, msg

    return passed, ''


def test_complementary_x1(test_complementary):
    passed, assertion_msg, *_ = test_complementary
    assert passed, f'Failed test_complementary(). {assertion_msg}'
        

def test_complementary_x2(test_complementary):
    passed, assertion_msg, *_ = test_complementary
    assert passed, f'Failed test_complementary(). {assertion_msg}'


def test_complementary_x3(test_complementary):
    passed, assertion_msg, *_ = test_complementary
    assert passed, f'Failed test_complementary(). {assertion_msg}'
