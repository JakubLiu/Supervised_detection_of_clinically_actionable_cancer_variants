import numpy as np
import tensorflow as tf
from tensorflow.keras import layers
import os
import glob
import argparse
import datetime
from keras.models import load_model
import scipy
from model_definitions import EpsilonTO



print(" ")
print(" ")

print(r"""
      ███████╗██████╗ ███████╗██╗██╗      ██████╗ ███╗   ██╗
      ██╔════╝██╔══██╗██╔════╝██║██║     ██╔═══██╗████╗  ██║
      █████╗  ██████╔╝███████╗██║██║     ██║   ██║██╔██╗ ██║
      ██╔══╝  ██╔═══╝ ╚════██║██║██║     ██║   ██║██║╚██╗██║
      ███████╗██║     ███████║██║███████╗╚██████╔╝██║ ╚████║
      ╚══════╝╚═╝     ╚══════╝╚═╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝
""")

print(" ")
print(datetime.datetime.now())
print(" ")


# ===================================================================================================================================

def get_npz_files(data_dir):

    files = glob.glob(os.path.join(data_dir, "rank*", "*.npz"))

    np.random.shuffle(files)

    split = int(0.7 * len(files))

    train_files = files[:split]
    val_files = files[split:]

    return train_files, val_files



# reads the numpy objects and stores them as numpy arrays
def npz_generator(npz_files):

    for npz_file in npz_files:

        bundle = np.load(npz_file, allow_pickle=True)

        reads = bundle["reads"].astype(np.float32)

        if reads.shape == (0,):
            continue

        context = bundle["context"].astype(np.float32)
        error_rate = np.float32(bundle["error_rate"])

        yield reads, context, error_rate


# stores the numpy arrays as tensorflow datasets
def create_dataset(npz_files):

    return tf.data.Dataset.from_generator(

        lambda: npz_generator(npz_files),

        output_signature=(
            tf.TensorSpec(shape=(None,None,9), dtype=tf.float32),
            tf.TensorSpec(shape=(None,4), dtype=tf.float32),
            tf.TensorSpec(shape=(), dtype=tf.float32)
        )
    )


# ===================================================================================================================================
# the same parameters as during training need to be used
datdir = 'data_call'
files, val_files = get_npz_files(datdir)
dataset = create_dataset(files)


model = EpsilonTO(
    _projection_dim=256,
    _num_attention_heads_set=8,
    _genomic_context_size=1001,
    _base_embedding_dim=128,
    _dense_dim=256,
    _num_attention_heads_context=8
)

reads, context, error_rate = next(iter(dataset))

_ = model(
    pileup_input=reads,
    genomic_context_input=context
)

model.load_weights("best_model.weights.h5")

pvalues_raw = []
posteriors = []
prior = 0.3

for reads, context, error_rate in dataset:

    predictions = model(
        pileup_input=reads,
        genomic_context_input=context,
        training=False
    )



    coverage_obs = reads.shape[0]                            # observed coverage
    error_rate_obs = error_rate.numpy()                      # observed error rate
    alt_count_obs = int(np.floor(coverage_obs * error_rate_obs))  # observed alternative read count
    error_rate_pred = predictions.numpy()[0][0]              # predicted error rate



    # Binomial test
    pval = float(scipy.stats.binomtest(k = alt_count_obs, n = coverage_obs, p = error_rate_pred, alternative = 'greater').pvalue)

    # Bayesian posterior
    P_H0 = scipy.stats.binom.logpmf(
            k=alt_count_obs,
            n=coverage_obs,
            p=error_rate_obs
    )

    P_H1 = scipy.stats.binom.logpmf(
            k=alt_count_obs,
            n=coverage_obs,
            p=error_rate_pred
    )

    BF = np.exp(P_H0 - P_H1)
    posterior = float(BF * prior / (BF * prior + (1 - prior)))

    if alt_count_obs == 0.0:
            pvalues_raw.append(1.0)
            posteriors.append(0.0)
    else:
            pvalues_raw.append(pval)
            posteriors.append(posterior)


    print("-----------------------------------")
    print("Coverage:", coverage_obs)
    print("Observed error rate:", error_rate_obs)
    print("Observed alternative count:", alt_count_obs)
    print("Predicted error rate:", error_rate_pred)
    print("Binomial p-value:", pval)
    print("Bayesian posterior:", posterior)
    print("-----------------------------------")




print(pvalues_raw)
print(posteriors)

