import numpy as np
import keras
import tensorflow as tf
from poligon import *  # for Epsilon

print('Im here')

# load the data
validation_data = np.load('../fit_to_actual_data/validation_data.npz')
alt =  validation_data["y_true_val"][:,1]
cov =  validation_data["y_true_val"][:,0]
pileup_mask = (alt <= cov * 0.05) & (cov > 0) # keep only rows where the error rate is at most 5% and there is some coverage
all_zero_padding_mask = validation_data["x_free_sample_mask_val"].sum(axis=1) > 0 # keep only these samples that have a nonzero padding mask*
mask = pileup_mask & all_zero_padding_mask
x_instructive_genomic_val = validation_data["x_instructive_genomic_val"][mask]
x_instructive_sample_val = validation_data["x_instructive_sample_val"][mask]
x_free_sample_padded_val = validation_data["x_free_sample_padded_val"][mask]
x_free_genomic_val = validation_data["x_free_genomic_val"][mask]
x_free_sample_mask_val = validation_data["x_free_sample_mask_val"][mask]
y_true_val = validation_data["y_true_val"][mask]



print('Im here')

X_val_negative_control = tf.data.Dataset.from_tensor_slices((
    {
        'instructive_genomic_input': x_instructive_genomic_val,
        'instructive_sample_input': x_instructive_sample_val,
        'free_sample_padded_input': x_free_sample_padded_val,
        'free_genomic_input': x_free_genomic_val,
        'free_sample_padding_mask': x_free_sample_mask_val
    }
)).batch(32)


print('Im here')



# save the model and build it from config
model = tf.keras.models.load_model(
    "../saved_models/current_model.keras",
    custom_objects={
        "Epsilon": Epsilon,
        "BinomialNegLogLoss": BinomialNegLogLoss
    },
    compile=False,
    safe_mode=False
)


predictions = model.predict(X_val_negative_control)

np.savetxt('predicted_background_error_rate_negative_controls.txt', predictions)

print('done.')


