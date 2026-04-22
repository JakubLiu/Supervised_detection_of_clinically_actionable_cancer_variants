import numpy as np
import pandas as pd
import tensorflow as tf
import random
from mod import *  # for Epsilon
#tf.keras.mixed_precision.set_global_policy('mixed_float16')  # for mixed precision


class PrintPredictionsCallback(tf.keras.callbacks.Callback):
    def __init__(self, dataset):
        super().__init__()
        self.dataset = dataset

    def on_epoch_end(self, epoch, logs=None):
        print(f"\nEpoch {epoch + 1} predictions:")

        for x_batch, y_true in self.dataset:
            y_pred = self.model.predict(x_batch, verbose=0)

            print("True:", (y_true[:,1]/y_true[:,0]).numpy())
            print("Pred:", y_pred)



# ===================================== load datasets =====================================================
training_data = np.load('training_data.npz')

#train_subset_size = 500
#train_subset_size = training_data["x_instructive_genomic_train"].shape[0]

alt =  training_data["y_true_train"][:,1]
cov =  training_data["y_true_train"][:,0]
pileup_mask = (alt <= cov * 0.05) & (cov > 0) # keep only rows where the error rate is at most 5% and there is some coverage
all_zero_padding_mask = training_data["x_free_sample_mask_train"].sum(axis=1) > 0 # keep only these samples that have a nonzero padding mask*
mask = pileup_mask & all_zero_padding_mask # merge both masks



x_instructive_genomic_train = training_data["x_instructive_genomic_train"][mask]
x_instructive_sample_train = training_data["x_instructive_sample_train"][mask]
x_free_sample_padded_train = training_data["x_free_sample_padded_train"][mask]
x_free_genomic_train = training_data["x_free_genomic_train"][mask]
x_free_sample_mask_train = training_data["x_free_sample_mask_train"][mask]
y_true_train = training_data["y_true_train"][mask]


validation_data = np.load('validation_data.npz')


alt =  validation_data["y_true_val"][:,1]
cov =  validation_data["y_true_val"][:,0]
pileup_mask = (alt <= cov * 0.05) & (cov > 0) # keep only rows where the error rate is at most 5% and there is some coverage
all_zero_padding_mask = validation_data["x_free_sample_mask_val"].sum(axis=1) > 0 # keep only these samples that have a nonzero padding mask*
mask = pileup_mask & all_zero_padding_mask  # merge both masks



x_instructive_genomic_val = validation_data["x_instructive_genomic_val"][mask]
x_instructive_sample_val = validation_data["x_instructive_sample_val"][mask]
x_free_sample_padded_val = validation_data["x_free_sample_padded_val"][mask]
x_free_genomic_val = validation_data["x_free_genomic_val"][mask]
x_free_sample_mask_val = validation_data["x_free_sample_mask_val"][mask]
y_true_val = validation_data["y_true_val"][mask]



'''
* I should go into the actual data and see why some samples fail these filters, but since the training set
is large eitherway then I dont give a fuck and just remove these samples. These samples need to be removed
because the loss is not computed for them.
'''


# ============================ instantiate and compile the model ======================================
window_size = x_free_genomic_train.shape[1]
num_genomic_features = x_instructive_genomic_train.shape[1]
num_features_per_sample = x_instructive_sample_train.shape[1]
num_unique_bases = len(['A', 'T', 'G', 'C', 'N'])
padded_size = 10000
features_per_read = x_free_sample_padded_train.shape[2]


epsilon_model = Epsilon(
    intructive_genomic_input_dim=num_genomic_features,
    instructive_genomic_hidden_dim1=256,
    instructive_genomic_hidden_dim2=128,
    instructive_genomic_output_dim=64,

    instructive_sample_input_dim=num_features_per_sample,
    instructive_sample_hidden_dim1=256,
    instructive_sample_hidden_dim2=128,
    instructive_sample_output_dim=64,

    free_sample_hidden_dim1_psi=256,
    free_sample_hidden_dim2_psi=128,
    free_sample_output_dim_psi=64,
    free_sample_hidden_dim1_phi=256,
    free_sample_hidden_dim2_phi=128,
    free_sample_output_dim_phi=64,

    genomic_context_size=window_size,
    num_unique_bases=num_unique_bases,
    base_embedding_dim=256,
    dense_dim=128,
    num_heads=10,

    final_proj_dim1=256,
    final_proj_dim2=128
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
    optimizer=tf.keras.optimizers.AdamW(3e-5, clipnorm = 1), 
    loss=ZeroInfBinomialLoss(pseudocount=1e-5, normalize_by_coverage=True),  # the zero-inflated Binomial loss
    metrics=['mae']
)



# ====================================================== batch the data ==========================================================================

batch_size = 32

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

training_dataset = training_dataset.shuffle(buffer_size=y_true_train.shape[0]).batch(batch_size).prefetch(tf.data.AUTOTUNE)



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



validation_dataset = validation_dataset.shuffle(buffer_size=y_true_val.shape[0]).batch(batch_size).prefetch(tf.data.AUTOTUNE)



# ======================================================== train the model ===========================================================

history = functional_model.fit(
    x = training_dataset,
    validation_data=validation_dataset,
    epochs=15#,
    #callbacks = [PrintPredictionsCallback(validation_dataset)]
    )

train_loss = history.history['loss']
val_loss = history.history['val_loss']

np.save("train_loss.gpu.200ep.npy", train_loss)
np.save("val_loss.gpu.200ep.npy", val_loss)

functional_model.save("model1.keras")


print("all done.")


