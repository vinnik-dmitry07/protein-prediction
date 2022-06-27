# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 16:11:18 2018

@author: Michael
"""

import pytest
import json

import os
import math
import numpy as np

from tests import *


############ HELPER FUNCTIONS ##################
@pytest.fixture(scope="module")
def relative_path():
    return os.path.dirname(__file__)


@pytest.fixture(scope="module")
def json_data(relative_path):
    with open(os.path.join(relative_path, 'exe6_test.json')) as json_file:
        json_data = json.load(json_file)
    return json_data


@pytest.fixture(scope="module")
def bias(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["bias"]))


@pytest.fixture(scope="module")
def hinge_loss(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["hinge"]))


@pytest.fixture(scope="module")
def delta_hinge(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["delta_hinge"]))


@pytest.fixture(scope="module")
def l2_loss(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["l2_loss"]))


@pytest.fixture(scope="module")
def delta_l2(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["delta_l2"]))


@pytest.fixture(scope="module")
def sigmoid(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["sigmoid"]))


@pytest.fixture(scope="module")
def perceptron(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["perceptron"]))


@pytest.fixture(scope="module")
def perceptron_bias(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["perceptron_bias"]))


@pytest.fixture(scope="module")
def multiperceptron_bias(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["multiperceptron_bias"]))


@pytest.fixture(scope="module")
def multiperceptron_bias_nonlin(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["multiperceptron_bias_nonlin"]))


############ INIT STUDENT PERCEPTRON ######################
@pytest.fixture(scope="module")
def student_perceptron(relative_path, json_data):
    # learning rate, number of epochs and random seed for single perceptrons
    LEARNING_RATE = json_data['parameters']['learning_rate']
    NEPOCHS = json_data['parameters']['nepochs']
    SEED = json_data['parameters']['seed']

    try:
        from exe6_perceptron import Perceptron
        perc = Perceptron(LEARNING_RATE, NEPOCHS, SEED)
        return perc
    except:
        raise AssertionError('Something went wrong with initializing your Perceptron!') from None


############ TESTS ##################
@pytest.fixture(scope="module")
def adding_bias_term(bias, student_perceptron):
    # create 1D test array and 2D test array
    test_array_1D = bias['bias_1d_in']
    test_array_2D = bias['bias_2d_in']

    try:
        student_answer_1D = student_perceptron._add_bias(test_array_1D)
        student_answer_2D = student_perceptron._add_bias(test_array_2D)
    except Exception:
        raise AssertionError('Error in bias adding procedure.') from None

    try:
        assert student_answer_1D is not None, "No array was returned by the _add_bias function."
        assert student_answer_2D is not None, "No array was returned by the _add_bias function."
    except AssertionError as msg:
        return False, msg

    try:
        correct_answer_1D = np.all(np.isclose(bias['bias_1d_out'], student_answer_1D))
        correct_answer_2D = np.all(np.isclose(bias['bias_2d_out'], student_answer_2D))
    except Exception:
        raise AssertionError('Something with your bias return value went wrong.') from None

    try:
        assert correct_answer_1D, ("When adding a bias term to a 1D array you have to " +
                                   "add a 1 to the end of the array.")
        assert correct_answer_2D, ("When adding a bias term to a 2D array you have to " +
                                   "add a 1 to the end of each row in the array.")
        return True, ''
    except AssertionError as msg:
        return False, msg


def test_adding_bias_term(adding_bias_term):
    passed, assertion_msg, *_ = adding_bias_term
    assert passed, f'Failed test adding_bias_term(). {assertion_msg}'


@pytest.fixture(scope="module")
def _hinge_loss(hinge_loss, student_perceptron):
    for _, data in hinge_loss.items():
        y, y_pred = data[0], data[1]

        try:
            student_answer_hinge_loss = student_perceptron._hinge_loss(y, y_pred)
        except Exception:
            raise AssertionError('Your hinge loss function produces an error.') from None

        try:
            assert student_answer_hinge_loss is not None, (
                "You returned None instead of the hinge loss.")
        except AssertionError as msg:
            return False, msg

        try:
            correct_answer = math.isclose(data[2], student_answer_hinge_loss)
        except Exception:
            raise AssertionError('Your hinge loss does not seem to be well formated.') from None

        try:
            assert correct_answer, ("Your definition of hinge loss is not correct. " +
                                    "Please keep in mind that this is the loss and not the derivative! " +
                                    "Also keep in mind that groundtruth labels and predicted labels are " +
                                    "within [-1, 1], not within [0, 1].")
        except AssertionError as msg:
            return False, msg
    return True, ''


def test_hinge_loss(_hinge_loss):
    passed, assertion_msg, *_ = _hinge_loss
    assert passed, f'Failed test hinge_loss(). {assertion_msg}'


@pytest.fixture(scope="module")
def _delta_hinge(delta_hinge, student_perceptron):
    for _, data in delta_hinge.items():
        y, y_pred = data[0], data[1]

        try:
            student_answer_delta_hinge = student_perceptron._delta_hinge(y, y_pred)
        except Exception:
            raise AssertionError('Your delta hinge function seems to produce an error.') from None

        try:
            assert student_answer_delta_hinge is not None, "You returned None instead of the derivative hinge loss."
        except AssertionError as msg:
            return False, msg

        try:
            correct_answer = math.isclose(data[2], student_answer_delta_hinge)
        except Exception:
            raise AssertionError('Your delta hinge does not seem to be well formated.') from None

        try:
            assert correct_answer, "Your definition of the derivative of the hinge loss is not correct."
        except AssertionError as msg:
            return False, msg
    return True, ''


def test_delta_hinge(_delta_hinge):
    passed, assertion_msg, *_ = _delta_hinge
    assert passed, f'Failed test delta_hinge(). {assertion_msg}'


@pytest.fixture(scope="module")
def _l2_loss(l2_loss, student_perceptron):
    for _, data in l2_loss.items():
        y = data[0]
        y_pred = data[1]

        try:
            student_answer_l2_loss = student_perceptron._l2_loss(y, y_pred)
        except Exception:
            raise AssertionError('Your l2 loss seems to produce an error.') from None

        try:
            assert student_answer_l2_loss is not None, "You returned None instead of the l2 loss."
        except AssertionError as msg:
            return False, msg

        try:
            correct_answer = math.isclose(data[2], student_answer_l2_loss)
        except Exception:
            raise AssertionError('Your l2 loss does not seem to be well formated.') from None

        try:
            assert correct_answer, ("Your definition of l2 loss is not correct. " +
                                    "Please keep in mind that this is the loss and not the derivative.")
        except AssertionError as msg:
            return False, msg
    return True, ''


def test_l2_loss(_l2_loss):
    passed, assertion_msg, *_ = _l2_loss
    assert passed, f'Failed test l2_loss(). {assertion_msg}'


@pytest.fixture(scope="module")
def _delta_l2(delta_l2, student_perceptron):
    for _, data in delta_l2.items():
        y = data[0]
        y_pred = data[1]
        try:
            student_answer_delta_l2 = student_perceptron._delta_l2(y, y_pred)
        except Exception:
            raise AssertionError('Your delta l2 function seems to produce an error.') from None

        try:
            assert student_answer_delta_l2 is not None, (
                "You returned None instead of the correct derivative of the l2 loss.")
        except AssertionError as msg:
            return False, msg

        try:
            correct_answer = math.isclose(data[2], student_answer_delta_l2)
        except Exception:
            raise AssertionError('Your delta l2 does not seem to be well formated.') from None
        try:
            assert correct_answer, "Your definition of the derivative of the l2 loss is not correct."
        except AssertionError as msg:
            return False, msg
    return True, ''


def test_delta_l2(_delta_l2):
    passed, assertion_msg, *_ = _delta_l2
    assert passed, f'Failed test delta_l2(). {assertion_msg}'


@pytest.fixture(scope="module")
def _sigmoid(sigmoid, student_perceptron):
    try:
        student_answer_sigmoid = student_perceptron._sigmoid(sigmoid[0, :])
    except Exception:
        raise AssertionError('Your sigmoid function seems to produce an error.') from None

    try:
        assert student_answer_sigmoid is not None, (
            "You returned None instead of the correct sigmoid.")
    except AssertionError as msg:
        return False, msg
    try:
        correct_answer = np.all(np.isclose(sigmoid[1, :], student_answer_sigmoid))
    except Exception:
        raise AssertionError('Your sigmoid does not seem to be well formated.') from None
    return correct_answer, "Your definition of sigmoid is not correct."


def test_sigmoid(_sigmoid):
    passed, assertion_msg, *_ = _sigmoid
    assert passed, f'Failed test sigmoid(). {assertion_msg}'


##############################################################################
################################ 2 POINTS ####################################

@pytest.fixture(scope="module")
def _single_perceptron(perceptron, student_perceptron):
    try:
        student = student_perceptron.single_perceptron()
    except Exception:
        raise AssertionError('Your delta single_perceptron function seems to produce an error.') from None

    try:
        assert student is not None, (
            "You did not return any weights for the single_perceptron.")
    except AssertionError as msg:
        return False, msg

    try:
        correct_answer = np.all(np.isclose(perceptron, student))
    except Exception:
        raise AssertionError('Your single_perceptron does not seem to be well formated.') from None

    return correct_answer, ("Your single_perceptron did not return the correct weights " +
                            "based on the given learning rate and number of epochs." +
                            "Double check whether you are using the correct loss " +
                            " function ( here: hinge loss) and the correct targets (OR gate).")


def test_single_perceptron_x1(_single_perceptron):
    passed, assertion_msg, *_ = _single_perceptron
    assert passed, f'Failed test single_perceptron(). {assertion_msg}'


def test_single_perceptron_x2(_single_perceptron):
    passed, assertion_msg, *_ = _single_perceptron
    assert passed, f'Failed test single_perceptron(). {assertion_msg}'


############################################################################## 
################################ 2 POINTS ####################################

@pytest.fixture(scope="module")
def _single_perceptron_with_bias(perceptron_bias, student_perceptron):
    try:
        student = student_perceptron.single_perceptron_with_bias()
    except Exception:
        raise AssertionError('Your single_perceptron_with_bias function seems to produce an error.') from None

    try:
        assert student is not None, (
                "You did not return any weights for the single_perceptron with a bias term." +
                "Probably, you did not implement the _add_bias function, yet.")
    except AssertionError as msg:
        return False, msg

    try:
        correct_answer = np.all(np.isclose(perceptron_bias, student))
    except Exception:
        raise AssertionError('Your single_perceptron_with_bias does not seem to be well formated.') from None

    return correct_answer, ("Your single_perceptron_with_bias did not return the " +
                            "correct weights based on the given learning rate " +
                            "and number of epochs. Double check whether you " +
                            "are using the correct loss function ( here: hinge loss) " +
                            "and the correct targets (OR gate).")


def test_single_perceptron_with_bias_x1(_single_perceptron_with_bias):
    passed, assertion_msg, *_ = _single_perceptron_with_bias
    assert passed, f'Failed test single_perceptron_with_bias(). {assertion_msg}'


def test_single_perceptron_with_bias_x2(_single_perceptron_with_bias):
    passed, assertion_msg, *_ = _single_perceptron_with_bias
    assert passed, f'Failed test single_perceptron_with_bias(). {assertion_msg}'


##################### START TESTING MULTILAYER PERCEPTRONS ###################
################################ 5 POINTS ####################################

@pytest.fixture(scope="module")
def _multi_perceptron_with_bias(multiperceptron_bias, student_perceptron):
    try:
        student = student_perceptron.multi_perceptron_with_bias()
    except Exception:
        raise AssertionError('Your multi_perceptron_with_bias function seems to produce an error.') from None

    try:
        assert student is not None, (
                "You did not return any weights for the multi_perceptron with a bias term. " +
                "Probably, you did not implement the _add_bias function, yet.")
    except AssertionError as msg:
        return False, msg

    try:
        correct_answer = np.all(np.isclose(multiperceptron_bias, student))
    except Exception:
        raise AssertionError('Your multi_perceptron_with_bias does not seem to be well formated.') from None

    return correct_answer, ("Your multi_perceptron did not return the correct weights " +
                            "based on the given learning rate and number of epochs. " +
                            "Double check whether you are using the correct loss " +
                            " function ( here: l2 loss) and the correct targets (XOR gate).")


def test_multi_perceptron_with_bias_x1(_multi_perceptron_with_bias):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias
    assert passed, f'Failed test multi_perceptron_with_bias(). {assertion_msg}'


def test_multi_perceptron_with_bias_x2(_multi_perceptron_with_bias):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias
    assert passed, f'Failed test multi_perceptron_with_bias(). {assertion_msg}'


def test_multi_perceptron_with_bias_x3(_multi_perceptron_with_bias):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias
    assert passed, f'Failed test multi_perceptron_with_bias(). {assertion_msg}'


def test_multi_perceptron_with_bias_x4(_multi_perceptron_with_bias):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias
    assert passed, f'Failed test multi_perceptron_with_bias(). {assertion_msg}'


def test_multi_perceptron_with_bias_x5(_multi_perceptron_with_bias):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias
    assert passed, f'Failed test multi_perceptron_with_bias(). {assertion_msg}'


############################################################################## 
################################ 5 POINTS ####################################

@pytest.fixture(scope="module")
def _multi_perceptron_with_bias_and_nonlinearity(multiperceptron_bias_nonlin, student_perceptron):
    try:
        student = student_perceptron.multi_perceptron_with_bias_and_nonlinearity()
    except Exception:
        raise AssertionError('Your multi_perceptron_with_bias_and_nonlinearity '
                             'function seems to produce an error.') from None

    try:
        assert student is not None, (
                "You did not return any weights for the multi_perceptron with a bias term. " +
                "Probably, you did not implement the _add_bias function, yet.")
    except AssertionError as msg:
        return False, msg

    try:
        correct_answer = np.all(np.isclose(multiperceptron_bias_nonlin, student))
    except:
        raise AssertionError('Your multi_perceptron_with_bias_and_nonlinearity '
                             'does not seem to be well formated.') from None

    return correct_answer, ("Your multi_perceptron did not return the correct weights " +
                            "based on the given learning rate and number of epochs. " +
                            "Double check whether you are using the correct loss " +
                            " function ( here: l2 loss) and the correct targets (XOR gate).")


def test_multi_perceptron_with_bias_and_nonlinearity_x1(_multi_perceptron_with_bias_and_nonlinearity):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias_and_nonlinearity
    assert passed, f'Failed test multi_perceptron_with_bias_and_nonlinearity(). {assertion_msg}'


def test_multi_perceptron_with_bias_and_nonlinearity_x2(_multi_perceptron_with_bias_and_nonlinearity):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias_and_nonlinearity
    assert passed, f'Failed test multi_perceptron_with_bias_and_nonlinearity(). {assertion_msg}'


def test_multi_perceptron_with_bias_and_nonlinearity_x3(_multi_perceptron_with_bias_and_nonlinearity):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias_and_nonlinearity
    assert passed, f'Failed test multi_perceptron_with_bias_and_nonlinearity(). {assertion_msg}'


def test_multi_perceptron_with_bias_and_nonlinearity_x4(_multi_perceptron_with_bias_and_nonlinearity):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias_and_nonlinearity
    assert passed, f'Failed test multi_perceptron_with_bias_and_nonlinearity(). {assertion_msg}'


def test_multi_perceptron_with_bias_and_nonlinearity_x5(_multi_perceptron_with_bias_and_nonlinearity):
    passed, assertion_msg, *_ = _multi_perceptron_with_bias_and_nonlinearity
    assert passed, f'Failed test multi_perceptron_with_bias_and_nonlinearity(). {assertion_msg}'
