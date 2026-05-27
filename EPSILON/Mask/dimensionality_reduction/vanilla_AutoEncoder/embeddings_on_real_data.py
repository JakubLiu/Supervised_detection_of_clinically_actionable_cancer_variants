import numpy as np
import pandas as pd

import os
#os.environ["TF_XLA_FLAGS"] = "--tf_xla_enable_xla_devices=false"
#os.environ["TF_ENABLE_ONEDNN_OPTS"] = "0"

import tensorflow as tf
#tf.config.optimizer.set_jit(False)
#tf.config.run_functions_eagerly(True)


import random

from standard_AE_mod import *

X_orig = np.loadtxt('X.csv', dtype = np.float32, skiprows = 1, delimiter = ',')

model = StandardAE(
    original_input_size = X_orig.shape[1],
    hidden_layer_size = 128,
    bottleneck_size = 62
)


model.compile(
            optimizer=tf.keras.optimizers.Adam(1e-3),
            loss=tf.keras.losses.MeanSquaredError()
        )


early_stopping = EarlyStopping(
                        monitor='loss',  
                        patience=5,
                        min_delta=1e-3,           # minimum change to count as improvement
                        restore_best_weights=True,
                        verbose=1                 # prints when stopping
                    )


model.fit(
            x=X_orig,
            y=X_orig,
            batch_size=32,
            epochs=1000,
            callbacks=[early_stopping]
        )


# use the trained model to create the embedding and save the embedding
embedding = model.encoder(X_orig, training = False).numpy()
np.savetxt('standard_AE_embedding.txt', embedding, fmt = '%.8f')

print('embedding saved')