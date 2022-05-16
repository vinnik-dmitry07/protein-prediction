import json
import pytest

from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def relative_path():
    return Path(__file__).parent


@pytest.fixture(scope="module")
def json_data(relative_path):
    test_json = 'exe2_swissprot_test.json'

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def filename(relative_path, json_data):
    return Path(relative_path, json_data['filename'])


@pytest.fixture(scope="module")
def identification(json_data):
    return json_data['identification']


@pytest.fixture(scope="module")
def sequence_info(json_data):
    return json_data['seq_info']


@pytest.fixture(scope="module")
def organism(json_data):
    return json_data['organism']


@pytest.fixture(scope="module")
def localization(json_data):
    return set(json_data['localization'])


@pytest.fixture(scope="module")
def pdb_support(json_data):
    return set(json_data['pdb_support'])


@pytest.fixture(scope="module")
def sp_parser(filename):
    try:
        from exe2_swissprot import SwissProt_Parser
        sp_parser = SwissProt_Parser(filename)
        return sp_parser
    except Exception:
        raise AssertionError('Error while creating SwissProt_Parser. Failed to initialize SwissProt_Parser.') from None


@pytest.fixture(scope="module")
def sp_get_sp_identification(sp_parser, identification):
    try:
        student_id = sp_parser.get_sp_identification()
    except Exception:
        raise AssertionError('Error in SwissProt_Parser.get_sp_identification().') from None

    try:
        passed = isinstance(student_id, tuple)
        assert passed, 'Return value is not a tuple.'

        passed = all([isinstance(content, str) for content in student_id])
        assert passed, 'Not all entries are strings.'

        length_ident = len(identification)
        passed = (length_ident == len(student_id))
        assert passed, f'Expected {length_ident} many items but {len(student_id)} received.'

        passed = (identification == list(student_id))
        return passed, 'Incorrect SwissProt identification.'
    except AssertionError as msg:
        return False, msg


def test_sp_get_sp_identification(sp_get_sp_identification):
    passed, assertion_msg, *_ = sp_get_sp_identification
    assert passed, f'Failed test sp_get_sp_identification(). {assertion_msg}'


@pytest.fixture(scope="module")
def sp_get_sp_sequence_info(sp_parser, sequence_info):
    try:
        student_info = sp_parser.get_sp_sequence_info()
    except Exception:
        raise AssertionError('Error in SwissProt_Parser.get_sp_sequence_info().') from None

    try:
        passed = isinstance(student_info, tuple)
        assert passed, 'Return value is not a tuple.'

        length_info = len(sequence_info)
        passed = (length_info == len(student_info))
        assert passed, f'Expected {length_info} many items but {len(student_info)} received.'

        passed = all([isinstance(content, instance) for content, instance in zip(student_info, [str, int, int])])
        assert passed, 'Incorrect data types of the entries.'

        passed = (sequence_info == list(student_info))
        return passed, 'Incorrect sequence information.'
    except AssertionError as msg:
        return False, msg


def test_sp_get_sp_sequence_info(sp_get_sp_sequence_info):
    passed, assertion_msg, *_ = sp_get_sp_sequence_info
    assert passed, f'Failed test sp_get_sp_sequence_info(). {assertion_msg}'


@pytest.fixture(scope="module")
def sp_get_organism(sp_parser, organism):
    try:
        org = sp_parser.get_organism()
    except Exception:
        raise AssertionError('Error in SwissProt_Parser.get_organism().') from None

    hint = 'Check the field [organism].'
    passed = (organism == org)
    return passed, f'Incorrect organism. {hint}'


def test_sp_get_organism(sp_get_organism):
    passed, assertion_msg, *_ = sp_get_organism
    assert passed, f'Failed test sp_get_organism(). {assertion_msg}'


@pytest.fixture(scope="module")
def sp_get_localization(sp_parser, localization):
    try:
        loc = sp_parser.get_localization()
    except Exception:
        raise AssertionError('Error in SwissProt_Parser.get_localization().') from None

    try:
        hint = 'Check the field [comment_subcellularlocation_location].'
        passed = (len(localization) == len(loc))
        assert passed, f'Incorrect number of localizations. {hint}'

        passed = (localization == set(loc))
        return passed, f'Incorrect localization(s). {hint}'
    except AssertionError as msg:
        return False, msg


def test_sp_get_localization(sp_get_localization):
    passed, assertion_msg, *_ = sp_get_localization
    assert passed, f'Failed test sp_get_localization(). {assertion_msg}'


@pytest.fixture(scope="module")
def sp_get_pdb_support(sp_parser, pdb_support):
    try:
        pdb = sp_parser.get_pdb_support()
    except Exception:
        raise AssertionError('Error in SwissProt_Parser.get_pdb_support().') from None

    try:
        hint = 'Check the field [dbxrefs].'
        passed = (len(pdb_support) == len(pdb))
        assert passed, f'Incorrect number of PDB IDs. {hint}'

        passed = (pdb_support == set(pdb))
        return passed, f'Incorrect PDB ID(s). {hint}'
    except AssertionError as msg:
        return False, msg


def test_sp_get_pdb_support(sp_get_pdb_support):
    passed, assertion_msg, *_ = sp_get_pdb_support
    assert passed, f'Failed test sp_get_pdb_support(). {assertion_msg}'
