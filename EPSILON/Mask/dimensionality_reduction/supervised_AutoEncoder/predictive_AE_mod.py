import tensorflow as tf
import numpy as np
from tensorflow.keras.layers import Input, Dense, Dropout, Embedding, Flatten, LayerNormalization
from tensorflow.keras.callbacks import EarlyStopping, History, ModelCheckpoint
from tensorflow.keras.utils import plot_model
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.initializers import VarianceScaling
from tensorflow.keras.optimizers import Adam, AdamW, Nadam



class PredictiveAE(tf.keras.Model):
    
    def __init__(self, original_input_size, hidden_layer_size, bottleneck_size, activation_fun = 'relu', dropout_rate = 0.0):
        super().__init__()
    

        # the encoder
        self.encoder = tf.keras.Sequential([
            tf.keras.layers.Dense(hidden_layer_size, activation=activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(bottleneck_size, activation=activation_fun)
        ])


        # the decoder
        self.decoder = tf.keras.Sequential([
            tf.keras.layers.Dense(hidden_layer_size, activation=activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(original_input_size)
        ])

        # use the encoded input to make a prediction
        self.predictor = tf.keras.Sequential([
            tf.keras.layers.Dense(16, activation=activation_fun),
            tf.keras.layers.Dense(2, activation="sigmoid")  # two output nodes for the two parameters of a zero-inflated binomial model
        ])


    def call(self, x, training = True):
        
        encoded_input = self.encoder(x, training = training)
        reconstructed_input = self.decoder(encoded_input, training = training)
        ZIB_estimates= self.predictor(encoded_input, training = training)

        return reconstructed_input, ZIB_estimates    

'''
            ============ how to compile the model =========================

            model.compile(
            optimizer=tf.keras.optimizers.Adam(1e-3),

            loss=[
                tf.keras.losses.MeanSquaredError(),   # -------> reconstruction loss
                tf.keras.losses.MeanSquaredError()   # ------------> prediction loss
            ],

            loss_weights=[1.0, 1.0]      # a hyperparameter
                )
'''


'''
          =================== how to run the model =====================

          model.fit(
            x=train_x,
            y=[train_x, train_y],   # ----> note that one of the targets is the input itself
            batch_size=32,
            epochs=50
        )
'''








class ZeroInfBinomialLoss(tf.keras.losses.Loss):
    '''
    https://doi.org/10.1111/j.0006-341X.2000.01030.x
    '''

    def __init__(self, pseudocount=1e-8, normalize_by_coverage=False, name="zero_inflated_biomial_loss"):
        super().__init__(name=name)
        self.eps = pseudocount
        self.normalize = normalize_by_coverage

    def call(self, y_true, y_pred):

        '''
        n --> the coverage (from the data)
        k --> the alt read count (from the data)
        pi --> the probability of k=0
        p --> the binomial success probability
        '''

        y_true = tf.cast(y_true, tf.float64)
        y_pred = tf.cast(y_pred, tf.float64)

        # make it robust for weird samples
        n = tf.maximum(y_true[:, 0], 0.0)
        k = tf.maximum(y_true[:, 1], 0.0)
        k = tf.minimum(k, n)

        #y_pred = tf.clip_by_value(y_pred, self.eps, 1.0 - self.eps)

        pi = y_pred[:,0]
        p = y_pred[:,1]
        pi = tf.clip_by_value(pi, self.eps, 1.0 - self.eps)
        p  = tf.clip_by_value(p, self.eps, 1.0 - self.eps)
        

        # the log-likelihood if k=0
        log_lik_zero = tf.math.log(pi + (1.0 - pi) * tf.math.exp(n*tf.math.log(1.0-p)))

        # the log-likelihood if k>0
        standard_binomial_likelihood = k*tf.math.log(p) + (n-k)*tf.math.log(1.0-p)
        log_lik_nonzero = (tf.math.log(1.0-pi)) + standard_binomial_likelihood

        # if k=0 --> log_lik_zero,    else (i.e. k>0) ----> log_lik_nonzero
        log_likelihood = tf.where(k == 0, log_lik_zero, log_lik_nonzero)

        if self.normalize:
            log_likelihood = log_likelihood / (n + self.eps)

        log_likelihood = tf.clip_by_value(log_likelihood, -1e6, 1e6) # clip extreme values to stabilize the loss


        #tf.print(" mean pi per batch:", tf.reduce_mean(pi), " mean p per batch:", tf.reduce_mean(p), ' mean error per batch: ', tf.reduce_mean((1 - pi) * p), '\n') # ----> for debugging !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        '''
        THEN DURING INFERENCE, IN ORDER TO GET ONE ESTIMATE OF THE ERROR RATE (BINOMIAL SUCCESS PROBABILITY PARAMETER ESTIMATE) YOU NEED TO DO THE FOLLOWING:

        y_pred = model(inputs, training=False)

        pi = y_pred[:, 0]
        p  = y_pred[:, 1]

        error_rate_estimate = (1 - pi) * p
        '''


        return -tf.reduce_mean(log_likelihood)