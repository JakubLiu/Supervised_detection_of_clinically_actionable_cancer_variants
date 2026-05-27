import tensorflow as tf
import numpy as np
from tensorflow.keras.layers import Input, Dense, Dropout, Embedding, Flatten, LayerNormalization
from tensorflow.keras.callbacks import EarlyStopping, History, ModelCheckpoint
from tensorflow.keras.utils import plot_model
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.initializers import VarianceScaling
from tensorflow.keras.optimizers import Adam, AdamW, Nadam



class StandardAE(tf.keras.Model):
    
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


    def call(self, x, training = True):
        
        encoded_input = self.encoder(x, training = training)
        reconstructed_input = self.decoder(encoded_input, training = training)

        return reconstructed_input
    


'''
            ============ how to compile the model =========================

            model.compile(
            optimizer=tf.keras.optimizers.Adam(1e-3),
            loss=tf.keras.losses.MeanSquaredError()
        )
'''


'''
          =================== how to run the model =====================

          model.fit(
            x=train_x,
            y=train_x,
            batch_size=32,
            epochs=50
        )
'''