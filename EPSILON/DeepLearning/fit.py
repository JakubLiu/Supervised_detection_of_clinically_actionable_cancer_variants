import numpy as np
import pandas as pd
import tensorflow as tf
import random
from mod import *  # for Epsilon



# ===================================== load datasets =====================================================
training_data = np.load('training_data.npz')

x_instructive_genomic_train = training_data["x_instructive_genomic_train"]
x_instructive_sample_train = training_data["x_instructive_sample_train"]
x_free_sample_padded_train = training_data["x_free_sample_padded_train"]
x_free_genomic_train = training_data["x_free_genomic_train"]
x_free_sample_mask_train = training_data["x_free_sample_mask_train"]
y_true_train = training_data["y_true_train"]


validation_data = np.load('validation_data.npz')

x_instructive_genomic_val = validation_data["x_instructive_genomic_val"]
x_instructive_sample_val = validation_data["x_instructive_sample_val"]
x_free_sample_padded_val = validation_data["x_free_sample_padded_val"]
x_free_genomic_val = validation_data["x_free_genomic_val"]
x_free_sample_mask_val = validation_data["x_free_sample_mask_val"]
y_true_val = validation_data["y_true_val"]





# ============================ instantiate and compile the model ======================================
window_size = x_free_genomic_train.shape[1]
num_genomic_features = x_instructive_genomic_train.shape[1]
num_features_per_sample = x_instructive_sample_train.shape[1]
num_unique_bases = len(['A', 'T', 'G', 'C', 'N'])
padded_size = 10000
features_per_read = x_free_sample_padded_train.shape[2]


epsilon_model = Epsilon(
    intructive_genomic_input_dim=num_genomic_features,
    instructive_genomic_hidden_dim1=64,
    instructive_genomic_hidden_dim2=32,
    instructive_genomic_output_dim=16,

    instructive_sample_input_dim=num_features_per_sample,
    instructive_sample_hidden_dim1=64,
    instructive_sample_hidden_dim2=32,
    instructive_sample_output_dim=16,

    free_sample_hidden_dim1_psi=64,
    free_sample_hidden_dim2_psi=32,
    free_sample_output_dim_psi=16,
    free_sample_hidden_dim1_phi=64,
    free_sample_hidden_dim2_phi=32,
    free_sample_output_dim_phi=16,

    genomic_context_size=window_size,
    num_unique_bases=num_unique_bases,
    base_embedding_dim=128,
    dense_dim=32,
    num_heads=4,

    final_proj_dim1=32,
    final_proj_dim2=16
)


instructive_genomic_input = tf.keras.Input(shape=(num_genomic_features,), name="instructive_genomic_input")
instructive_sample_input = tf.keras.Input(shape=(num_features_per_sample,), name="instructive_sample_input")
free_sample_padded_input = tf.keras.Input(shape=(padded_size, features_per_read), name="free_sample_padded_input")
free_genomic_input = tf.keras.Input(shape=(window_size,), name="free_genomic_input")
free_sample_padding_mask = tf.keras.Input(shape=(padded_size,), name="free_sample_padding_mask")

output = epsilon_model(
    instructive_genomic_input,
    instructive_sample_input,
    free_sample_padded_input,
    free_genomic_input,
    free_sample_padding_mask
)

functional_model = tf.keras.Model(
    inputs={
        'instructive_genomic_input': instructive_genomic_input,
        'instructive_sample_input': instructive_sample_input,
        'free_sample_padded_input': free_sample_padded_input,
        'free_genomic_input': free_genomic_input,
        'free_sample_padding_mask': free_sample_padding_mask
    },
    outputs=output
)

functional_model.compile(
    optimizer=tf.keras.optimizers.Adam(1e-3),
    loss=BinomialNegLogLoss(pseudocount=1e-8, normalize_by_coverage=False),
    metrics=['mae']
)



# ====================================================== batch the data ==========================================================================

batch_size = 8

training_dataset = tf.data.Dataset.from_tensor_slices((
    {
        'instructive_genomic_input': x_instructive_genomic_train,
        'instructive_sample_input': x_instructive_sample_train,
        'free_sample_padded_input': x_free_sample_padded_train,
        'free_genomic_input': x_free_genomic_train,
        'free_sample_padding_mask': x_free_sample_mask_train
    },
    y_true_train
))

training_dataset = training_dataset.shuffle(buffer_size=y_true_train.shape[0]).batch(batch_size)


validation_dataset = tf.data.Dataset.from_tensor_slices((
    {
        'instructive_genomic_input': x_instructive_genomic_val,
        'instructive_sample_input': x_instructive_sample_val,
        'free_sample_padded_input': x_free_sample_padded_val,
        'free_genomic_input': x_free_genomic_val,
        'free_sample_padding_mask': x_free_sample_mask_val
    },
    y_true_val
))

validation_dataset = validation_dataset.shuffle(buffer_size=y_true_val.shape[0]).batch(batch_size)



# ======================================================== train the model ===========================================================

history = functional_model.fit(
    x = training_dataset,
    validation_data=validation_dataset,
    epochs=200)

train_loss = history.history['loss']
val_loss = history.history['val_loss']

np.save("train_loss.npy", train_loss)
np.save("val_loss.npy", val_loss)

print("all done.")
