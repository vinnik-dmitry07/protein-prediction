import json
import pytest

from math import isclose
from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe1_genome_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def test_cases(json_data):
    cases = [json_data[i] for i in json_data]
    return cases


@pytest.fixture(scope="module")
def sequences(test_cases):
    return [case.get('sequence', None) for case in test_cases]


@pytest.fixture(scope="module")
def at_content(test_cases):
    return [case.get('at_content', None) for case in test_cases]


@pytest.fixture(scope="module")
def codon_dist(test_cases):
    return [case.get('codon_dist', None) for case in test_cases]


@pytest.fixture(scope="module")
def amino_acid_dist(test_cases):
    return [case.get('amino_acid_dist', None) for case in test_cases]


@pytest.fixture(scope="module")
def genome(sequences):
    try:
        from exe1_genome import Genome
        genomes = [Genome(sequence) for sequence in sequences]
        return genomes
    except Exception:
        raise AssertionError('Error while creating Genome. Genome initialization failed.') from None


@pytest.fixture(scope="module")
def genome_get_at_content(genome, at_content):
    try:
        ats = [gen_case.get_at_content() for gen_case in genome]
    except Exception:
        raise AssertionError('Error in Genome.get_at_content().') from None

    try:
        passed = all([isclose(at_case, at_solution, rel_tol=0, abs_tol=1e-5)
                      for at_case, at_solution in zip(ats, at_content)])
        return passed, 'Incorrect AT content.'
    except Exception:
        raise AssertionError('Incorrect AT content.') from None


def test_genome_get_at_content(genome_get_at_content):
    passed, assertion_msg, *_ = genome_get_at_content
    assert passed, f'Failed test genome_get_at_content(). {assertion_msg}'


def compare_dicts(obj_a, obj_b):
    if not isinstance(obj_a, dict) or not isinstance(obj_b, dict):
        try:
            return isclose(obj_a, obj_b, rel_tol=0, abs_tol=1e-5)
        except Exception:
            return False

    keys_a = set(obj_a.keys())
    keys_b = set(obj_b.keys())

    if keys_a != keys_b:
        return False

    return all([compare_dicts(obj_a[key], obj_b[key]) for key in keys_a])


@pytest.fixture(scope="module")
def test_genome_get_codon_dist(genome, codon_dist):
    try:
        cds = [gen_case.get_codon_dist() for gen_case in genome]
    except Exception:
        raise AssertionError('Error in Genome.get_codon_dist().') from None

    passed = all([compare_dicts(codon_solution, cd_case)
                  for codon_solution, cd_case in zip(codon_dist, cds)])
    return passed, 'Incorrect codon distribution.'


def test_genome_get_codon_dist_x1(test_genome_get_codon_dist):
    passed, assertion_msg, *_ = test_genome_get_codon_dist
    assert passed, f'Failed test_genome_get_codon_dist(). {assertion_msg}'


def test_genome_get_codon_dist_x2(test_genome_get_codon_dist):
    passed, assertion_msg, *_ = test_genome_get_codon_dist
    assert passed, f'Failed test_genome_get_codon_dist(). {assertion_msg}'


@pytest.fixture(scope="module")
def test_genome_get_amino_acid_dist(genome, amino_acid_dist):
    try:
        aads = [gen_case.get_amino_acid_dist() for gen_case in genome]
    except Exception:
        raise AssertionError('Error in Genome.get_amino_acid_dist().') from None

    passed = all([compare_dicts(aad_solution, aad_case)
                  for aad_solution, aad_case in zip(amino_acid_dist, aads)])
    return passed, 'Incorrect amino acid distribution.'


def test_genome_get_amino_acid_dist_x1(test_genome_get_amino_acid_dist):
    passed, assertion_msg, *_ = test_genome_get_amino_acid_dist
    assert passed, f'Failed test_genome_get_amino_acid_dist(). {assertion_msg}'


def test_genome_get_amino_acid_dist_x2(test_genome_get_amino_acid_dist):
    passed, assertion_msg, *_ = test_genome_get_amino_acid_dist
    assert passed, f'Failed test_genome_get_amino_acid_dist(). {assertion_msg}'
