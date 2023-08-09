# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 19:55:16 2018

@author: Michael
"""

import numpy as np

class Perceptron:
    '''
        This exercise focuses on the very basic concepts of feed-forward neural 
        networks. This class of networks summarizes the traditional multi-layer 
        neural networks (MLPs or fully-connected network ) as well as 
        state-of-the-art convolutional neural networks (CNNs).
        During this exercise you will implement the most basic building block
        of these networks: the perceptron.
        After completing this exercise you will understand why you need to add
        a bias-term to your network, why you need more than one hidden layer to
        find an approximation/solution for the XOR gate and why you need a
        non-linearity to connect these layers.
    '''

    def __init__( self, learning_rate=0.1, nEpochs=1000, seed=42 ):
        '''
            The Perceptron class is a wrapper for different neural network
            architectures. Summarizing the different architectures in one class
            allows you to easily experiement with different loss-functions or 
            activation functions.
            During the initialization we will only define basic variables which
            are shared by all network architectures, such as number of training
            epochs, learning rate and inputs as well as targets.

            Loss function: Function which measures the difference between the
            label (or value) predicted by the network and the true class label 
            ( or true value ).

            Number of training epochs: During every epoch all training samples
            are fed once in the network and their predicted outputs are compared 
            to the true class label (supervised learning). The weights in the 
            network are adjusted according to this difference using the 
            backpropagation algorithm. Here the weights are updated after each 
            sample (in contrast to batch learning).

            Learning rate: This variable allows you to adjust the magnitude of 
            every learning step ( == update of weights). If you set this value 
            close to 1 you speed up the learning process at the risk of missing 
            the actual minimum of your cost function based on the given inputs.
            If you set the value too small (e.g. 1e-10) your network will take
            very long to find a minimum.
        '''
        # inputs; here, boolean variables are used as an input to a logic operation  
        self.IN   = np.array([[0,0], [0,1], [1,0], [1,1]])
        # targets for single_layer perceptron; a simple OR function encoded within [-1, 1]
        self.OR   = np.array([ -1,     1,     1,     1 ])
        # targets for the multilayer perceptron; a XOR function encoded within [0,1]
        # The difference in the encoding is due to sigmoid non-linearity which outputs within [0,1]
        self.XOR  = np.array([ 0,     1,     1,    0 ])

        self.nEpochs       = nEpochs # number of training epochs (desc. see above)
        self.learning_rate = learning_rate # learning reate of the network (desc. see above)
        self.seed          = seed # random seed for random weight init.
        self.__reset_rnd_state()

    def __reset_rnd_state( self ):
        '''
            Helper function to ensure that random seed is set correctly before
            random weight initialization.
        '''
        self.rnd = np.random.RandomState( self.seed ) # seed to control randomness

    def _add_bias( self, arr ):
        '''
            This function is used to add a bias term to a given numpy array arr.
            A bias term is a constant factor which is introduced in each layer
            by adding an additional node which does not receive any inputs from
            previously layers. Usually, a 1 is used as a bias term. 
            The introduction of this bias term allows the network to
            find decision borders which do not cross the origin.
            Your function should handle 1D and 2D arrays as described below.
            Input:
                arr: Numpy array (1D or 2D).
            Return:
                padded version of the input numpy array. In case of a 1D numpy
                array ( e.g. [0, 0, 0] ), add a single 1 to the end of the
                array ( e.g. [0, 0, 0, 1] ).
                In case of a 2D numpy array, add a single 1 to every row.
                E.g. [[ 0, 0 ], [0, 0]] should become [[0, 0, 1], [0, 0, 1]].
        '''
        if len(arr.shape) == 1:
            arr = np.append(arr, 1)
        elif len(arr.shape) == 2:
            arr = np.column_stack((arr, [1] * arr.shape[0]))
        else:
            raise ValueError(arr.shape)
        return arr


    def _hinge_loss( self, y, y_pred ):
        '''
            This function returns the 'Hinge loss' or SVM loss.
            Please note, that this loss requires the class label to be either
            positive (e.g. +1) or negative (e.g. -1). If you encode one class
            with 0, this loss does not return any loss for this class.

            Input:
                y: the true class label of a sample
                y_pred: the prediction for this sample

            Output:
                Hinge loss based on the predicted class (y_pred) and true class (y)
        '''
        hinge_loss = np.maximum(1 - y * y_pred, 0)
        return hinge_loss

    def _delta_hinge( self, y, y_pred ):
        '''
            This function should return the derivative of the hinge loss based
            on the true class (y) and the prediction (y_pred).
            Please note again, that the two classes have to be encoded as 
            positive (e.g. +1) and negative values (e.g. -1) in order to produce
            correct gradients for both classes.

            Input:
                y: the true class label of a sample
                y_pred: the prediction for this sample

            Output:
                Derivative of the hinge loss function based on the 
                predicted class (y_pred) and true class (y)
        '''
        delta_hinge = np.zeros_like(y)
        delta_hinge[y_pred * y < 1] = -y
        return delta_hinge  # calculate derivative of hinge_loss based on y & y_pred


    def _l2_loss( self, y, y_pred ):
        '''
            This function returns the L2 loss function based on the true
            class (or value) and the prediction (or value).
            Please understand that the constant factor of 0.5 is multiplied
            to allow straightforward derivation later on.

            Input:
                y: the true class label of a sample
                y_pred: the predicted class label for this sample

            Output:
                L2 loss based on true class label (y) and prediction (y_pred)

        '''
        l2_loss = np.mean(np.sum(np.square(y - y_pred))) / 2
        return l2_loss # calculate l2 loss based on y & y_pred

    def _delta_l2( self, y, y_pred ):
        '''
            This function should return the derivative of L2 loss function based
            on the true class label (y) and the predicted class label (y_pred).
            Please understand that previously a constant factor of 0.5 was multiplied
            previously to the L2 norm in order to allow straightforward derivation.

            Input:
                y: the true class label of a sample
                y_pred: the prediction for this sample

            Output:
                Derivative of L2 loss based on true class label (y) and 
                the prediction (y_pred)
        '''
        delta_l2 = -np.mean(y - y_pred)  # calculate l2 loss based on y & y_pred
        return delta_l2

    def _sigmoid( self, x):
        '''
            Uses a sigmoid non-linearity to tranform a given input x.
            A non-linearity is necessary to stack several layers of neurons.
            If no non-linearity is used, the benefit of stacking several layers
            vanishes as multiple successive linear transformations can always
            be replaced by a single linear transformation.
            Only the introduction of a non-linearity allows the network to
            actually benefit from the increased complexity.

            Input:
                x: input to sigmoid transformation
            Output:
                Sigmoid transformation of input x
        '''
        sigmoid = 1 / (1 + np.exp(-x))
        return sigmoid # calculate sigmoid based on input x

    def single_perceptron( self ):
        '''
            In this function you have to implement the most basic unit of every
            feed-forward neural network: A perceptron.
            This perceptron will be trained on two boolean inputs [0,1] in order
            to approximate a OR gate:

                        IN_1    IN_2    OUT
                         0       0       -1
                         0       1        1
                         1       0        1
                         1       1        1

            Train the network for n epochs (defined during initialization; see
            self.nEpochs). Process every boolean combination seperately and 
            update the weights after each sample (see 'Online Machine Learning').
            Do not include a bias term. Please understand that this very easy
            model will not be able to find the optimal solution for the given
            problem. This is intended. Only after you've implemented a bias term
            ( see next task ), your model will be able to find an optimal solution
            for the OR gate.

            Inputs:
                None ( only the object itself aka self.)

            Outputs:
                return your weights w1 after training
        '''
        d_in = 2 # number of inputs
        self.__reset_rnd_state() # reset random state to make weight init. deterministic
        w1 = self.rnd.rand( d_in ) # randomly inititalize weights for given number of inputs

        for epoch in range(self.nEpochs): # for each training epoch
            for index, x in enumerate( self.IN ): # for each 'sample' in the input data

                y      =  self.OR[index]  # get target label for input; use OR gate states

                # forward pass: Make predictions for input
                y_pred = w1 @ x # get predictions

                # backward pass: Update weight(s) based on discrepancy between
                # predicted and true target; use _hinge_loss to measure loss 
                # and use _delta_hinge to update weights
                loss = self._hinge_loss( y, y_pred ) # error
                dw   = self._delta_hinge( y, y_pred  ) # gradient of error

                # update weights here based on learning rate, dw and 
                # inputs from previous layer
                w1   -= self.learning_rate * dw * x

        # return weights after training
        return w1

    def single_perceptron_with_bias( self ):
        '''
            Now we will extend the function 'single_perceptron' by adding a
            constant factor ( bias term ) to the input. This allows the network
            to find decision borders which do not cross the origin.
            This step is equivalent to adding an intersection with the y-axis
            to a linear equation of this form: 
                without bias-term:                  with bias term:
                    y(x) = a*x                          y(x) = a*x + b
            Inputs:
                None ( only the object itself aka self.)

            Outputs:
                return your weights w1
        '''
        d_in = 3 # number of inputs. 2 boolean inputs + 1 bias term
        self.__reset_rnd_state() # reset random state to make weight init. deterministic
        w1 = self.rnd.rand( d_in ) # randomly initialize weights for every input and bias
        X  = self._add_bias( self.IN ) # add bias ( 1 ) to every input

        if X is None:
            return X

        for epoch in range(self.nEpochs): # for every epoch
            for index, x in enumerate( X ): # for every 'sample' in the input data

                y      = self.OR[index]  # get target label for input; use OR gate states

                # forward pass: Make predictions for input
                y_pred = w1 @ x  # get predictions

                # backward pass: Update weight(s) based on discrepancy between
                # predicted and true target
                loss = self._hinge_loss( y, y_pred )   # error
                dw   = self._delta_hinge( y, y_pred  ) # gradient of error  

                # update weights here based on learning rate, dw and 
                # inputs from previous layer
                w1   -= self.learning_rate * dw * x

        # return weights after training
        return w1

    def multi_perceptron_with_bias( self ):
        '''
            After completing the previous tasks, you've successfully approximated
            an OR gate using a perceptron. This was possible as the data was
            seperable by a single decision border.
            Now we want to look at a slightly more complicated case: the XOR gate.

                           IN_1    IN_2    OUT
                            0       0        0
                            0       1        1
                            1       0        1
                            1       1        0

            Please understand that this problem is not solvable by drawing a
            single decision border. 
            Still, a neural network can solve this problem by adding an additional 
            hidden layer to the network. Again, we will use a constant bias
            term (bias==1) for every layer in this network. This bias term does
            not receive any inputs from previous layers.
            Please understand that this problem will not be solved by this simple
            architecture. As outlined above, multiple successive linear layers 
            can be replaced by a single linear layer, meaning that we did not 
            win anything from adding complexity to the architecture.
            However, only a slight modification to the architecture will be 
            required to solve the problem (see next task).
            The solution for the problem will require the targets to be within
            [0,1] instead of [-1,1]; therefore, the class labels are now within
            [0,1].
            Also, the loss function changes now from Hinge loss to L2 loss.
            This requires a smaller learning rate of e.g 0.1 or 0.01.
        '''
        d_in  = 3 # numer of input nodes (again: 2 inputs + 1 bias )
        d_out = 2 # number of nodes in the hidden layer

        # number of nodes in the hidden layer + 1 bias term.
        # This number is not equal to d_out as bias terms do not reveive any
        # input from previous layers.
        d_out_bias = d_out + 1 

        self.__reset_rnd_state() # reset random state to make weight init. deterministic

        w1 = self.rnd.rand( d_in,      d_out ) # input to hidden weights
        w2 = self.rnd.rand( d_out_bias       ) # hidden to output weights
        X  = self._add_bias( self.IN ) # add bias to inputs

        if X is None:
            return X

        for epoch in range(self.nEpochs): # for every epoch
            for index, x in enumerate( X ): # for every sample in the dataset

                # Forward pass
                y      = self.XOR[index] # get target label for input; now XOR
                h      = x[np.newaxis] @ w1  # forward pass; inputs to hidden
                h_bias = self._add_bias( h ) # add bias to intermediate layer
                y_pred = h_bias @ w2 # hidden to output

                # backward pass: Update weight(s) based on discrepancy between
                # predicted and true target. Watch out: Backward pass differs
                # between output to hidden update and hidden to hidden updates.
                loss        = self._l2_loss( y, y_pred  )
                grad_y_pred = self._delta_l2( y, y_pred )
                dw2 = grad_y_pred
                dw1 = dw2 * w2

                w1 -= self.learning_rate * x[:, np.newaxis] @ dw1[:-1][np.newaxis]
                w2 -= self.learning_rate * dw2 * h_bias.ravel()

        # return the weights in the hidden layer after training
        return w2

    def multi_perceptron_with_bias_and_nonlinearity( self ):
        '''
            As already mentioned, stacking several linear layers does not 
            increase the predictive power of the network as many linear layers
            can simply be replaced a single linear layer. However, by adding a
            non-linearity (such as sigmoid or rectified linear unit) in the 
            hidden layer allows the network to find non-linear solutions to a problem.
            In the following task you will extend your previous multi-layer 
            architecture ( 'multi_perceptron_with_bias' ) by adding a
            non-linearity after the hidden layer. This allows you to solve non-linear 
            problems.
        '''
        d_in  = 3 # numer of input nodes (again: 2 inputs + 1 bias )
        d_out = 2 # number of nodes in the hidden layer

        # number of nodes in the hidden layer + 1 bias term.
        # This number is not equal to d_out as bias terms do not reveive any
        # input from previous layers.
        d_out_bias = d_out + 1 

        self.__reset_rnd_state() # reset random state to make weight init. deterministic
        w1 = self.rnd.rand( d_in,      d_out ) # input to hidden weights
        w2 = self.rnd.rand( d_out_bias       ) # hidden to output weights
        X  = self._add_bias( self.IN ) # add bias to inputs

        if X is None:
            return X

        for epoch in range(self.nEpochs): # for every epoch 
            for index, x in enumerate( X ): # for every sample in the data set

                # Forward pass
                y      = self.XOR[index] # get target label for input; now XOR

                h      = x[np.newaxis] @ w1  # forward pass; inputs to hidden.
                h_non_lin = self._sigmoid(h) # use non-linear transformation (here: sigmoid)
                h_bias = self._add_bias( h_non_lin ) # add bias to intermediate layer
                y_pred = h_bias @ w2 # hidden to output

                # backward pass: Update weight(s) based on discrepancy between
                # predicted and true target. Watch out: Backward pass differs
                # between output to hidden update and hidden to hidden updates.
                loss        = self._l2_loss( y, y_pred  )
                grad_y_pred = self._delta_l2( y, y_pred )
                dw2 = grad_y_pred
                dw1 = dw2 * w2

                sigma = h_non_lin * (1 - h_non_lin)

                w1 -= self.learning_rate * x[:, np.newaxis] @ (sigma * dw1[:-1][np.newaxis])
                w2 -= self.learning_rate * dw2 * h_bias.ravel()

        # return the weights in the hidden layer after training
        return w2


def main():
    perceptron = Perceptron()
    perceptron.multi_perceptron_with_bias_and_nonlinearity()


if __name__ == '__main__':
    main()
