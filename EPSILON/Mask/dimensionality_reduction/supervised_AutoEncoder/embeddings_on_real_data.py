import numpy as np
import pandas as pd

import os

import tensorflow as tf
#tf.config.optimizer.set_jit(False)
#tf.config.run_functions_eagerly(True)


import random

from predictive_AE_mod import *


"""
                                        Notes about the model

        - it has two losses
            - one is the reconstruction loss (MSE)  ---> pushes the model to reconstruct the input well 
            - the other is the prediction loss (custom zero-inflated binomial loss (same as in Epsilon DL with a negative control cohort)) ---> pushes the model to predicts the error rate well*
                    * more specifically to predict the parameters of a zero-inflated binomial model well
        - the ZIB loss is usually 10x-20x smaller in scale than the MSE loss
        - therefore I upweight it 10x to balance out the contribution of it in the total loss
             - because I upweight it 10x then it is often still smaller than the MSE loss
             - but this is okay since I want the AE to create embeddings that are good for the prediction
        - I also do not do any train-test split, because I do not care about the generalization, because I just want to extract good embeddings for this specific dataset
        - offcourse I do not want the model to create embeddings that are based on learning noise, thatshwy I regularize the model and use early stopping**
                    ** these things do not guarantee that the model will not learn noise, but at least they do no harm

"""



X_orig = np.loadtxt('X.csv', dtype = np.float32, skiprows = 1, delimiter = ',')
alt_counts = np.loadtxt('alt_counts.txt', dtype = np.float32, skiprows = 1)
ref_counts = np.loadtxt('ref_counts.txt', dtype = np.float32, skiprows = 1)
n = (alt_counts + ref_counts).reshape(-1,1)
k = alt_counts.reshape(-1,1)
y = np.concatenate([n,k], axis = 1) # the custom zero-inflated binomial loss function needs y to be in this form y[:,0] = coverage, y[:,1] = alt counts

print(X_orig.shape)
print(n.shape)
print(k.shape)
print(y.shape)

model = PredictiveAE(
    original_input_size = X_orig.shape[1],
    hidden_layer_size = 128,
    bottleneck_size = 62,
    dropout_rate = 0.3
)


model.compile(
            optimizer=tf.keras.optimizers.Adam(1e-3),

            loss=[
                tf.keras.losses.MeanSquaredError(),   # -------> reconstruction loss
                ZeroInfBinomialLoss(pseudocount=1e-5, normalize_by_coverage=False)   # ------------> prediction loss
            ],

            loss_weights=[1.0, 10.0]      # contribution of each loss
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
            y=[X_orig, y],   # ----> note that one of the targets is the input itself
            batch_size=32,
            epochs=100,  # probably never used in practice due to early stopping
            callbacks = [early_stopping]
        )


# use the trained model to create the embedding and save the embedding
embedding = model.encoder(X_orig, training = False).numpy()
np.savetxt('predictive_AE_embedding.txt', embedding, fmt = '%.8f')

print('embedding saved')