import json
import pytest
import numpy as np

from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def relative_path():
    return Path(__file__).parent


@pytest.fixture(scope="module")
def json_data(relative_path):
    test_json = 'exe2_pdb_test.json'

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def filename(relative_path, json_data):
    return Path(relative_path, json_data['filename'])


@pytest.fixture(scope="module")
def pdb_parser(filename):
    try:
        from exe2_pdb import PDB_Parser
        pdb_parser = PDB_Parser(filename)
        return pdb_parser
    except Exception:
        raise AssertionError('Error while creating PDB_Parser. Failed to initialize PDB_Parser.') from None


@pytest.fixture(scope="module")
def pdb_get_number_of_chains(pdb_parser, json_data):
    try:
        num_chains = pdb_parser.get_number_of_chains()
    except Exception:
        raise AssertionError('Error in PDB_Parser.get_number_of_chains().') from None

    try:
        passed = isinstance(num_chains, int)
        assert passed, 'Return value is not an integer.'

        hint = 'General structure of Bio.PDB: Structure -> Model -> Chains -> Residues -> Atoms.'
        passed = (json_data['number_of_chains'] == num_chains)
        return passed, f'Incorrect number of chains. {hint}'
    except AssertionError as msg:
        return False, msg


def test_pdb_get_number_of_chains(pdb_get_number_of_chains):
    passed, assertion_msg, *_ = pdb_get_number_of_chains
    assert passed, f'Failed test pdb_get_number_of_chains(). {assertion_msg}'


@pytest.fixture(scope="module")
def pdb_get_sequence(pdb_parser, json_data):
    try:
        seq = pdb_parser.get_sequence(json_data['sequence_chain'])
    except Exception:
        raise AssertionError('Error in PDB_Parser.get_sequence().') from None

    try:
        passed = isinstance(seq, str)
        assert passed, 'Return value is not a string.'

        passed = (json_data['sequence'] == seq)
        return passed, 'Incorrect sequence.'
    except AssertionError as msg:
        return False, msg


def test_pdb_get_sequence(pdb_get_sequence):
    passed, assertion_msg, *_ = pdb_get_sequence
    assert passed, f'Failed test pdb_get_sequence(). {assertion_msg}'


@pytest.fixture(scope="module")
def pdb_get_number_of_water_molecules(pdb_parser, json_data):
    try:
        chain = json_data['water_molecules_chain']
        num_water = pdb_parser.get_number_of_water_molecules(chain)
    except Exception:
        raise AssertionError('Error in PDB_Parser.get_number_of_water_molecules().') from None

    try:
        passed = isinstance(num_water, int)
        assert passed, 'Return value is not an integer.'

        hint = 'Water molecules are called HOH and are part of Bio.PDBs chains.'
        passed = (json_data['number_of_water_molecules'] == num_water)
        return passed, f'Incorrect number of water molecules. {hint}'
    except AssertionError as msg:
        return False, msg


def test_pdb_get_number_of_water_molecules(pdb_get_number_of_water_molecules):
    passed, assertion_msg, *_ = pdb_get_number_of_water_molecules
    assert passed, f'Failed test pdb_get_number_of_water_molecules(). {assertion_msg}'


@pytest.fixture(scope="module")
def pdb_get_ca_distance_same_chain(pdb_parser, json_data):
    try:
        ca_dis = pdb_parser.get_ca_distance(
            json_data["ca_distance_same_chains"]["chain_1"],
            json_data["ca_distance_same_chains"]["id_1"],
            json_data["ca_distance_same_chains"]["chain_2"],
            json_data["ca_distance_same_chains"]["id_2"]
        )
    except Exception:
        raise AssertionError('Error in PDB_Parser.get_ca_distance().') from None

    try:
        passed = isinstance(ca_dis, int)
        assert passed, 'Return value is not an integer.'

        passed = (json_data["ca_distance_same_chains"]["result"] == ca_dis)
        return passed, 'Incorrect C-alpha distance (same chain).'
    except AssertionError as msg:
        return False, msg


def test_pdb_get_ca_distance_same_chain(pdb_get_ca_distance_same_chain):
    passed, assertion_msg, *_ = pdb_get_ca_distance_same_chain
    assert passed, f'Failed test pdb_get_ca_distance_same_chain(). {assertion_msg}'


@pytest.fixture(scope="module")
def pdb_get_ca_distance_different_chains(pdb_parser, json_data):
    try:
        ca_dis = pdb_parser.get_ca_distance(
            json_data["ca_distance_diff_chains"]["chain_1"],
            json_data["ca_distance_diff_chains"]["id_1"],
            json_data["ca_distance_diff_chains"]["chain_2"],
            json_data["ca_distance_diff_chains"]["id_2"]
        )
    except Exception:
        raise AssertionError('Error in PDB_Parser.get_ca_distance().') from None

    try:
        passed = isinstance(ca_dis, int)
        assert passed, 'Return value is not an integer.'

        passed = (json_data["ca_distance_diff_chains"]["result"] == ca_dis)
        return passed, 'Incorrect C-alpha distance (different chains).'
    except AssertionError as msg:
        return False, msg


def test_pdb_get_ca_distance_different_chains(pdb_get_ca_distance_different_chains):
    passed, assertion_msg, *_ = pdb_get_ca_distance_different_chains
    assert passed, f'Failed test pdb_get_ca_distance_different_chains(). {assertion_msg}'


@pytest.fixture(scope="module")
def test_pdb_get_bfactors(pdb_parser, json_data, relative_path):
    bfactors = json_data['bfactors']
    for chain, file in bfactors.items():
        try:
            bf = pdb_parser.get_bfactors(chain)
        except Exception:
            raise AssertionError('Error in PDB_Parser.get_bfactors().') from None

        try:
            passed = (type(bf) == np.ndarray and bf.dtype == np.int64)
            assert passed, 'B-Factors are not a numpy.ndarray of dtype numpy.int64.'

            bfactor = np.load(Path(relative_path, file))
            passed = (bfactor.shape == bf.shape)
            assert passed, 'The shape of the numpy array is incorrect.'

            passed = np.all((bfactor == bf) | (np.isnan(bfactor) & np.isnan(bf)))

            assert passed, 'Incorrect B-Factors.'
        except AssertionError as msg:
            return False, msg
    return True, ''


def test_pdb_get_bfactors_x1(test_pdb_get_bfactors):
    passed, assertion_msg, *_ = test_pdb_get_bfactors
    assert passed, f'Failed test test_pdb_get_bfactors(). {assertion_msg}'


def test_pdb_get_bfactors_x2(test_pdb_get_bfactors):
    passed, assertion_msg, *_ = test_pdb_get_bfactors
    assert passed, f'Failed test test_pdb_get_bfactors(). {assertion_msg}'


def test_pdb_get_bfactors_x3(test_pdb_get_bfactors):
    passed, assertion_msg, *_ = test_pdb_get_bfactors
    assert passed, f'Failed test test_pdb_get_bfactors(). {assertion_msg}'


def test_pdb_get_bfactors_x4(test_pdb_get_bfactors):
    passed, assertion_msg, *_ = test_pdb_get_bfactors
    assert passed, f'Failed test test_pdb_get_bfactors(). {assertion_msg}'


def test_pdb_get_bfactors_x5(test_pdb_get_bfactors):
    passed, assertion_msg, *_ = test_pdb_get_bfactors
    assert passed, f'Failed test test_pdb_get_bfactors(). {assertion_msg}'


@pytest.fixture(scope="module")
def test_pdb_get_contact_map(pdb_parser, json_data, relative_path):
    contacts = json_data['contacts']
    for chain, file in contacts.items():
        try:
            c_map = pdb_parser.get_contact_map(chain)
        except Exception:
            raise AssertionError('Error in PDB_Parser.get_contact_map().') from None

        try:
            passed = (type(c_map) == np.ndarray and c_map.dtype == np.int64)
            assert passed, 'Contact map is not a numpy.ndarray of dtype numpy.int64.'

            contact_map = np.load(Path(relative_path, file))
            passed = (contact_map.shape == c_map.shape)
            assert passed, f'The shape of the numpy array is incorrect.'

            passed = np.all((contact_map == c_map)
                            | (np.isnan(contact_map) & np.isnan(c_map)))
            assert passed, 'Incorrect contact map.'
        except AssertionError as msg:
            return False, msg
    return True, ''


def test_pdb_get_contact_map_x1(test_pdb_get_contact_map):
    passed, assertion_msg, *_ = test_pdb_get_contact_map
    assert passed, f'Failed test test_pdb_get_contact_map(). {assertion_msg}'


def test_pdb_get_contact_map_x2(test_pdb_get_contact_map):
    passed, assertion_msg, *_ = test_pdb_get_contact_map
    assert passed, f'Failed test test_pdb_get_contact_map(). {assertion_msg}'


def test_pdb_get_contact_map_x3(test_pdb_get_contact_map):
    passed, assertion_msg, *_ = test_pdb_get_contact_map
    assert passed, f'Failed test test_pdb_get_contact_map(). {assertion_msg}'


def test_pdb_get_contact_map_x4(test_pdb_get_contact_map):
    passed, assertion_msg, *_ = test_pdb_get_contact_map
    assert passed, f'Failed test test_pdb_get_contact_map(). {assertion_msg}'


def test_pdb_get_contact_map_x5(test_pdb_get_contact_map):
    passed, assertion_msg, *_ = test_pdb_get_contact_map
    assert passed, f'Failed test test_pdb_get_contact_map(). {assertion_msg}'
