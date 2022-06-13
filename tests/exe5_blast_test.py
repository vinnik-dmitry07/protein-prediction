import json
import pytest
import numpy as np

from time import time as ttime
from pathlib import Path
from tests import *


@pytest.fixture(scope="module")
def json_data():
    test_json = 'exe5_blast_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def db_sequences(json_data):
    return json_data['db_sequences']


@pytest.fixture(scope="module")
def db_stats(json_data):
    return tuple(json_data['db_stats'])


@pytest.fixture(scope="module")
def db_seqs_for_word(json_data):
    return json_data['db_seqs_for_word']


@pytest.fixture(scope="module")
def sub_matrix(json_data):
    return np.array(json_data['sub_matrix'], dtype=np.int64)


@pytest.fixture(scope="module")
def query_seq(json_data):
    return json_data['query_seq']


@pytest.fixture(scope="module")
def query_word(json_data):
    return json_data['query_word']


@pytest.fixture(scope="module")
def query_pssm(json_data):
    return np.array(json_data['query_pssm'], dtype=np.int64)


@pytest.fixture(scope="module")
def blast_words(json_data):
    return json_data['blast_words']


@pytest.fixture(scope="module")
def blast_words_pssm(json_data):
    return json_data['blast_words_pssm']


def table_list_tuple(data):
    for key, value in data.items():
        data[key] = [tuple(x) for x in value]
    return data


@pytest.fixture(scope="module")
def blast_hsp_one_hit(json_data):
    return table_list_tuple(json_data['blast_hsp_one_hit'])


@pytest.fixture(scope="module")
def blast_hsp_two_hit(json_data):
    return table_list_tuple(json_data['blast_hsp_two_hit'])


@pytest.fixture(scope="module")
def blast_hsp_one_hit_pssm(json_data):
    return table_list_tuple(json_data['blast_hsp_one_hit_pssm'])


@pytest.fixture(scope="module")
def blast_hsp_two_hit_pssm(json_data):
    return table_list_tuple(json_data['blast_hsp_two_hit_pssm'])


@pytest.mark.timeout(30)
@pytest.fixture(scope="module")
def blast_db_get_db_stats(db_sequences, db_stats):
    try:
        from exe5_blast import BlastDB
        blast_db = BlastDB()
    except Exception:
        raise AssertionError('Error while creating BlastDB. BlastDB initialization failed.') from None

    try:
        for s in db_sequences:
            blast_db.add_sequence(s)
    except Exception:
        raise AssertionError('Error in BlastDB.add_sequence().') from None

    try:
        stats = blast_db.get_db_stats()
    except Exception:
        raise AssertionError('Error in BlastDB.get_db_stats().') from None

    passed = (db_stats == stats)
    return passed, 'Incorrect BlastDB statistics.', blast_db


def test_blast_db_get_db_stats(blast_db_get_db_stats):
    passed, assertion_msg, *_ = blast_db_get_db_stats
    assert passed, f'Failed test blast_db_get_db_stats(). {assertion_msg}'


@pytest.mark.timeout(5)
@pytest.fixture(scope="module")
def blast_db_get_sequences(blast_db_get_db_stats, query_word, db_seqs_for_word):
    _, _, blast_db, *_ = blast_db_get_db_stats
    try:
        seqs_for_word = blast_db.get_sequences(query_word)
    except Exception:
        raise AssertionError('Error in BlastDB.get_sequences().') from None

    try:
        passed_1 = (len(db_seqs_for_word) == len(seqs_for_word))
        passed_2 = (set(db_seqs_for_word) == set(seqs_for_word))
    except Exception:
        raise AssertionError('Error while comparing BlastDB.get_sequences() output.') from None

    passed = (passed_1 and passed_2)
    return passed, 'Incorrect sequences returned.'


def test_blast_db_get_sequences(blast_db_get_sequences):
    passed, assertion_msg, *_ = blast_db_get_sequences
    assert passed, f'Failed test blast_db_get_sequences(). {assertion_msg}'


@pytest.mark.timeout(5)
@pytest.fixture(scope="module")
def blast_get_words(sub_matrix, query_seq, blast_words):
    try:
        from exe5_blast import Blast
        blast = Blast(sub_matrix)
    except Exception:
        raise AssertionError('Error while creating Blast. Blast initialization failed.') from None

    try:
        words = blast.get_words(sequence=query_seq, T=13)
    except Exception:
        raise AssertionError('Error in Blast.get_words(sequence).') from None

    try:
        passed_1 = (len(blast_words) == len(words))
        passed_2 = (set(blast_words) == set(words))
    except Exception:
        raise AssertionError('Error while comparing Blast.get_words(sequence) output.') from None

    passed = (passed_1 and passed_2)
    return passed, 'Incorrect words returned for sequence.', blast


def test_blast_get_words(blast_get_words):
    passed, assertion_msg, *_ = blast_get_words
    assert passed, f'Failed test blast_get_words(). {assertion_msg}'


@pytest.mark.timeout(5)
@pytest.fixture(scope="module")
def blast_get_words_with_pssm(blast_get_words, query_pssm, blast_words_pssm):
    _, _, blast, *_ = blast_get_words
    try:
        words = blast.get_words(pssm=query_pssm, T=11)
    except Exception:
        raise AssertionError('Error in Blast.get_words(pssm).') from None

    try:
        passed_1 = (len(blast_words_pssm) == len(words))
        passed_2 = (set(blast_words_pssm) == set(words))
    except Exception:
        raise AssertionError('Error while comparing Blast.get_words(pssm) output.') from None

    passed = (passed_1 and passed_2)
    return passed, 'Incorrect words returned for PSSM.'


def test_blast_get_words_with_pssm(blast_get_words_with_pssm):
    passed, assertion_msg, *_ = blast_get_words_with_pssm
    assert passed, f'Failed test blast_get_words_with_pssm(). {assertion_msg}'


@pytest.fixture(scope="module")
def blast_blastdb(blast_get_words, blast_db_get_db_stats):
    _, _, blast, *_ = blast_get_words
    _, _, blast_db, *_ = blast_db_get_db_stats
    return blast, blast_db


def compare_blast_results(blast_results, results, hint):
    try:
        try:
            passed_1 = (len(blast_results) == len(results))
            passed_2 = (set(blast_results) == set(results))
        except Exception:
            return False, f'Error while comparing Blast results ({hint}).', True

        passed = (passed_1 and passed_2)
        assert passed, f'Incorrect target sequences returned ({hint}).'

        for target, hsp_list in results.items():
            blast_hsp_list = blast_results[target]

            try:
                passed_1 = (len(blast_hsp_list) == len(hsp_list))
                passed_2 = (set(blast_hsp_list) == set(hsp_list))
            except Exception:
                return False, f'Error while comparing Blast HSP list ({hint}).', True

            passed = (passed_1 and passed_2)
            assert passed, f'Incorrect HSPs returned ({hint}).'

        return True, '', False
    except AssertionError as msg:
        return False, msg, False


@pytest.mark.timeout(60)
@pytest.fixture(scope="module")
def blast_search_one_hit(blast_blastdb, query_seq, blast_hsp_one_hit):
    blast, blast_db = blast_blastdb
    run_time_t1 = ttime()
    try:
        results = blast.search_one_hit(blast_db,
                                       query=query_seq,
                                       T=13,
                                       X=5,
                                       S=30)
    except Exception:
        raise AssertionError('Error in Blast.search_one_hit(sequence)') from None
    run_time_t2 = ttime()

    passed, assertion_msg, error, *_ = compare_blast_results(blast_hsp_one_hit, results, 'one-hit, sequence')
    assert not error, assertion_msg
    return passed, assertion_msg, int(run_time_t2 - run_time_t1)


def test_blast_search_one_hit(blast_search_one_hit):
    passed, assertion_msg, *_ = blast_search_one_hit
    assert passed, f'Failed test blast_search_one_hit(). {assertion_msg}'


def test_blast_search_one_hit_30s(blast_search_one_hit):
    passed, assertion_msg, runtime, *_ = blast_search_one_hit
    assert passed, f'Failed test blast_search_one_hit(). {assertion_msg}'
    assert runtime >= 0, 'Invalid runtime.'
    assert runtime <= 30, f'Too slow (one-hit, sequence): {runtime}s'


def test_blast_search_one_hit_15s(blast_search_one_hit):
    passed, assertion_msg, runtime, *_ = blast_search_one_hit
    assert passed, f'Failed test blast_search_one_hit(). {assertion_msg}'
    assert runtime >= 0, 'Invalid runtime.'
    assert runtime <= 15, f'Too slow (one-hit, sequence): {runtime}s'


@pytest.mark.timeout(60)
@pytest.fixture(scope="module")
def blast_search_two_hit(blast_blastdb, query_seq, blast_hsp_two_hit):
    blast, blast_db = blast_blastdb
    run_time_t1 = ttime()
    try:
        results = blast.search_two_hit(blast_db,
                                       query=query_seq,
                                       T=11,
                                       X=5,
                                       S=30,
                                       A=40)
    except Exception:
        raise AssertionError('Error in Blast.search_two_hit(sequence)') from None
    run_time_t2 = ttime()

    passed, assertion_msg, error, *_ = compare_blast_results(blast_hsp_two_hit, results, 'two-hit, sequence')
    assert not error, assertion_msg
    return passed, assertion_msg, int(run_time_t2 - run_time_t1)


def test_blast_search_two_hit(blast_search_two_hit):
    passed, assertion_msg, *_ = blast_search_two_hit
    assert passed, f'Failed test blast_search_two_hit(). {assertion_msg}'


def test_blast_search_two_hit_20s(blast_search_two_hit):
    passed, assertion_msg, runtime, *_ = blast_search_two_hit
    assert passed, f'Failed test blast_search_two_hit(). {assertion_msg}'
    assert runtime >= 0, 'Invalid runtime.'
    assert runtime <= 20, f'Too slow (two-hit, sequence): {runtime}s'


def test_blast_search_two_hit_10s(blast_search_two_hit):
    passed, assertion_msg, runtime, *_ = blast_search_two_hit
    assert passed, f'Failed test blast_search_two_hit(). {assertion_msg}'
    assert runtime >= 0, 'Invalid runtime.'
    assert runtime <= 10, f'Too slow (two-hit, sequence): {runtime}s'


@pytest.mark.timeout(60)
@pytest.fixture(scope="module")
def blast_search_one_hit_with_pssm(blast_blastdb, query_pssm, blast_hsp_one_hit_pssm):
    blast, blast_db = blast_blastdb
    run_time_t1 = ttime()
    try:
        results = blast.search_one_hit(blast_db,
                                       pssm=query_pssm,
                                       T=13,
                                       X=5,
                                       S=30)
    except Exception:
        raise AssertionError('Error in Blast.search_one_hit(pssm)') from None
    run_time_t2 = ttime()

    passed, assertion_msg, error, *_ = compare_blast_results(blast_hsp_one_hit_pssm, results, 'one-hit, pssm')
    assert not error, assertion_msg
    return passed, assertion_msg, int(run_time_t2 - run_time_t1)


def test_blast_search_one_hit_with_pssm(blast_search_one_hit_with_pssm):
    passed, assertion_msg, *_ = blast_search_one_hit_with_pssm
    assert passed, f'Failed test blast_search_one_hit_with_pssm(). {assertion_msg}'


def test_blast_search_one_hit_with_pssm_30s(blast_search_one_hit_with_pssm):
    passed, assertion_msg, runtime, *_ = blast_search_one_hit_with_pssm
    assert passed, f'Failed test blast_search_one_hit_with_pssm(). {assertion_msg}'
    assert runtime >= 0, 'Invalid runtime.'
    assert runtime <= 30, f'Too slow (one-hit, pssm): {runtime}s'


@pytest.mark.timeout(60)
@pytest.fixture(scope="module")
def blast_search_two_hit_with_pssm(blast_blastdb, query_pssm, blast_hsp_two_hit_pssm):
    blast, blast_db = blast_blastdb
    run_time_t1 = ttime()
    try:
        results = blast.search_two_hit(blast_db,
                                       pssm=query_pssm,
                                       T=11,
                                       X=5,
                                       S=30,
                                       A=40)
    except Exception:
        raise AssertionError('Error in Blast.search_two_hit(pssm)') from None
    run_time_t2 = ttime()

    passed, assertion_msg, error, *_ = compare_blast_results(blast_hsp_two_hit_pssm, results, 'two-hit, pssm')
    assert not error, assertion_msg
    return passed, assertion_msg, int(run_time_t2 - run_time_t1)


def test_blast_search_two_hit_with_pssm(blast_search_two_hit_with_pssm):
    passed, assertion_msg, *_ = blast_search_two_hit_with_pssm
    assert passed, f'Failed test blast_search_two_hit_with_pssm(). {assertion_msg}'


def test_blast_search_two_hit_with_pssm_20s(blast_search_two_hit_with_pssm):
    passed, assertion_msg, runtime, *_ = blast_search_two_hit_with_pssm
    assert passed, f'Failed test blast_search_two_hit_with_pssm(). {assertion_msg}'
    assert runtime >= 0, 'Invalid runtime.'
    assert runtime <= 20, f'Too slow (two-hit, pssm): {runtime}s'


def test_blast_search_all_15s(blast_search_one_hit, blast_search_two_hit, blast_search_one_hit_with_pssm,
                              blast_search_two_hit_with_pssm):
    passed, _, run_time_one_hit, *_ = blast_search_one_hit
    run_time_one_hit = 9999 if not passed else run_time_one_hit
    passed, _, run_time_two_hit, *_ = blast_search_two_hit
    run_time_two_hit = 9999 if not passed else run_time_two_hit
    passed, _, run_time_one_hit_pssm, *_ = blast_search_one_hit_with_pssm
    run_time_one_hit_pssm = 9999 if not passed else run_time_one_hit_pssm
    passed, _, run_time_two_hit_pssm, *_ = blast_search_two_hit_with_pssm
    run_time_two_hit_pssm = 9999 if not passed else run_time_two_hit_pssm

    run_time_max = max(run_time_one_hit,
                       run_time_two_hit,
                       run_time_one_hit_pssm,
                       run_time_two_hit_pssm)

    assert run_time_max >= 0, 'Invalid runtime.'
    assert run_time_max <= 15, f'Too slow (max runtime): {run_time_max}s'
