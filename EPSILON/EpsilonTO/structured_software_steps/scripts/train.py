import numpy as np
import os
import glob
import tensorflow as tf
from tensorflow.keras import layers
import argparse
import datetime
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











parser = argparse.ArgumentParser()

parser.add_argument(
            "--data_dir",
            required=True,
            type = str,
            help="Path to the directory containing all the .npz files"
        )

parser.add_argument(
            "--max_epochs",
            required=False,
            type = int,
            default=100,
            help="Maximum number of epochs to train the model for, incase the early stopping criterion is not invoked"
        )


parser.add_argument(
            "--zero_keep_prop_multiplier",
            required=False,
            type = float,
            default = 10.0,
            help="The value by which to multiply the probability of keeping a zero error rate site"
        )


parser.add_argument(
    "--optimizer",
    required=False,
    type=str,
    default="adam",
    help="Possible choices: ['adam', 'adamw', 'sgd', 'rmsprop', 'adagrad', 'adadelta', 'adamax', 'nadam', 'ftrl', 'lion']"
)

parser.add_argument(
            "--learning_rate",
            required=False,
            type = float,
            default = 1e-4,
            help="The learning rate."
        )


parser.add_argument(
            "--patience",
            required=False,
            type = int,
            default = 5,
            help="Patience for the early stopping"
        )

parser.add_argument(
            "--max_error_rate",
            required=False,
            type = float,
            default = 0.05,
            help="Loci that have an error rate above this threshold ill not be considered."
        )

parser.add_argument(
            "--output_weights",
            required=False,
            type = str,
            default = 'fitted.weights.h5'
        )



args = parser.parse_args()
# ============================================================== functions =====================================================================================

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



@tf.function(reduce_retracing=True)
def train_step(reads, context, y_true):

    with tf.GradientTape() as tape:

        y_pred = model(reads, context)

        loss = loss_function(
            tf.reshape(y_true, (1,1)),
            y_pred
        )

    gradients = tape.gradient(
        loss,
        model.trainable_variables
    )

    optimizer.apply_gradients(
        zip(gradients, model.trainable_variables)
    )


    tf.print('train loss:', loss, 'y_pred:', y_pred, 'y_true:', y_true)

    return loss



@tf.function(reduce_retracing=True)
def val_step(reads, context, y_true):

    y_pred = model(reads, context)

    loss = loss_function(
        tf.reshape(y_true, (1,1)),
        y_pred
    )

    tf.print('val loss:', loss)

    return loss




def estimate_zero_inflation(npz_files):

    print('####################################################################')
    print('Estimating zero inflation...')

    total = 0
    zero_count = 0
    nonzero_count = 0

    for reads, context, y in npz_generator(npz_files):

        total += 1

        if y == 0.0:
            zero_count += 1
        else:
            nonzero_count += 1


    zero_fraction = zero_count / total

    print("Total loci:", total)
    print("Zero error loci:", zero_count)
    print("Non-zero loci:", nonzero_count)
    print("Zero fraction:", zero_fraction)
    print('####################################################################')

    return zero_fraction


def subsample_zero_errors(reads, context, y, zero_keep_probability):

    def subsample():
        return tf.random.uniform([]) < zero_keep_probability

    # Keep all positive error sites
    keep_positive = (y > 0.0) & (y < args.max_error_rate)

    # Randomly keep zero error sites
    keep_zero = subsample()

    return tf.logical_or(keep_positive, keep_zero)


# ============================================================================================================================================================


# create the dataset
train_files, val_files = get_npz_files(args.data_dir)

train_dataset = create_dataset(train_files)
val_dataset = create_dataset(val_files)

zero_fraction = estimate_zero_inflation(train_files)

zero_keep_probability = min(1.0, (1-zero_fraction) / zero_fraction)

print("Keeping zero loci with probability:", zero_keep_probability*args.zero_keep_prop_multiplier)

train_dataset = train_dataset.filter(
    lambda reads, context, y:
        subsample_zero_errors(reads, context, y, zero_keep_probability=zero_keep_probability*args.zero_keep_prop_multiplier)
)


train_dataset = train_dataset.shuffle(buffer_size=1000,reshuffle_each_iteration=True)
train_dataset = train_dataset.prefetch(tf.data.AUTOTUNE)
val_dataset = val_dataset.prefetch(tf.data.AUTOTUNE)

# define the model
model = EpsilonTO(
    _projection_dim=256,
    _num_attention_heads_set=8,
    _genomic_context_size=1001,
    _base_embedding_dim=128,
    _dense_dim=256,
    _num_attention_heads_context=8
)

for reads, context, y in train_dataset.take(1):
    _ = model(reads, context)




if args.optimizer == "adam":
    optimizer = tf.keras.optimizers.Adam(learning_rate=args.learning_rate)

elif args.optimizer == "adamw":
    optimizer = tf.keras.optimizers.AdamW(learning_rate=args.learning_rate)

elif args.optimizer == "sgd":
    optimizer = tf.keras.optimizers.SGD(learning_rate=args.learning_rate)

elif args.optimizer == "rmsprop":
    optimizer = tf.keras.optimizers.RMSprop(learning_rate=args.learning_rate)

elif args.optimizer == "adagrad":
    optimizer = tf.keras.optimizers.Adagrad(learning_rate=args.learning_rate)

elif args.optimizer == "adadelta":
    optimizer = tf.keras.optimizers.Adadelta(learning_rate=args.learning_rate)

elif args.optimizer == "adamax":
    optimizer = tf.keras.optimizers.Adamax(learning_rate=args.learning_rate)

elif args.optimizer == "nadam":
    optimizer = tf.keras.optimizers.Nadam(learning_rate=args.learning_rate)

elif args.optimizer == "ftrl":
    optimizer = tf.keras.optimizers.Ftrl(learning_rate=args.learning_rate)

elif args.optimizer == "lion":
    optimizer = tf.keras.optimizers.Lion(learning_rate=args.learning_rate)

else:
    raise ValueError(f"Unknown optimizer: {args.optimizer}")



loss_function = tf.keras.losses.MeanSquaredError()


epochs = args.max_epochs

best_val_loss = np.inf
patience = args.patience
patience_counter = 0


for epoch in range(epochs):

    print("\nEpoch", epoch+1)


    train_losses = []

    for reads, context, y in train_dataset:

        loss = train_step(
            reads,
            context,
            y
        )

        train_losses.append(float(loss))


    mean_train_loss = np.mean(train_losses)



    # --------------------
    # VALIDATION
    # --------------------

    val_losses = []

    for reads, context, y in val_dataset:

        loss = val_step(
            reads,
            context,
            y
        )

        val_losses.append(float(loss))


    mean_val_loss = np.mean(val_losses)



    print(
        "train loss:",
        mean_train_loss,
        "val loss:",
        mean_val_loss
    )


    if mean_val_loss < best_val_loss:

        print("Validation improved")

        best_val_loss = mean_val_loss
        patience_counter = 0

        model.save_weights(args.output_weights)


    else:

        patience_counter += 1

        print(
            "No improvement:",
            patience_counter,
            "/",
            patience
        )


        if patience_counter >= patience:

            print("Early stopping")
            break

train_loss = np.array(train_losses, dtype = np.float32)
validation_loss = np.array(val_losses, dtype = np.float32)
np.savetxt('train_loss.txt', train_loss)
np.savetxt('validation_loss.txt', validation_loss)
