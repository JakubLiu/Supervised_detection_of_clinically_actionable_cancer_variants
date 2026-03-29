import tensorflow as tf
import numpy as np
from tensorflow.keras.layers import Input, Dense, Dropout, Embedding, Flatten, LayerNormalization
from tensorflow.keras.callbacks import EarlyStopping, History, ModelCheckpoint
from tensorflow.keras.utils import plot_model
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.initializers import VarianceScaling
from tensorflow.keras.optimizers import Adam, AdamW, Nadam


# ====================================== defining the individual components of the full model =============================================================

# the handcrafted genomic features sub-model________________________________________________________________
class InstructiveGenomic(tf.keras.Model):
    
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, output_dim, activation_fun = 'relu', dropout_rate = 0.0):
        super().__init__()
    

        self.MLP = tf.keras.Sequential([
            tf.keras.layers.Dense(hidden_dim1, activation = activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(hidden_dim2, activation = activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(output_dim, activation = activation_fun)
        ])

    def call(self, x, training = True):
        # shape of x: (batch_size, input_dim (number of genomic features))
        y = self.MLP(x, training = training)
        return y

# how to use as part of the full model
'''
x_dummy = tf.random.normal((32,7))

model_dummy = InstructiveGenomic(input_dim = 7,
                                hidden_dim1 = 64,
                                hidden_dim2 = 32,
                                output_dim = 16,
                                activation_fun = 'relu',
                                dropout_rate = 0.0)

model_dummy.compile(
    optimizer = 'adam',  # just an example
    loss = 'mse'  # just an example
)

y = model_dummy(x_dummy, training = True)   # during inference add training = False
'''









# the handcrated sample-specific features sub-model_____________________________________________________________________
class InstructiveSample(tf.keras.Model):
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, output_dim, activation_fun = 'relu', dropout_rate = 0.0):
        super().__init__()
    

        self.MLP = tf.keras.Sequential([
            tf.keras.layers.Dense(hidden_dim1, activation = activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(hidden_dim2, activation = activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(output_dim, activation = activation_fun)
        ])

    def call(self, x, training = True):
        # shape of x: (batch_size, input_dim (number of sample-specific features))
        y = self.MLP(x, training = training)
        return y



# how to use as part of the full model
'''
x_dummy = tf.random.normal((32,7))

model_dummy = InstructiveSample(input_dim = 7,
                                hidden_dim1 = 64,
                                hidden_dim2 = 32,
                                output_dim = 16,
                                activation_fun = 'relu',
                                dropout_rate = 0.0)

model_dummy.compile(
    optimizer = 'adam',  # just an example
    loss = 'mse'  # just an example
)

y = model_dummy(x_dummy, training = True)   # during inference add training = False
'''








# the free sample-specific feature sub-model (following the DeepSets idea: https://scibits.blog/posts/deepsets/)_________________________________________
class FreeSample(tf.keras.Model):
    def __init__(self,
                input_dim, hidden_dim1_psi, hidden_dim2_psi, output_dim_psi, dropout_rate_psi = 0.0, activation_fun_psi = 'relu',  # parameters of the psi MLP's
                aggregation_function = 'mean',                                                                                     # the aggregation function
                hidden_dim1_phi, hidden_dim2_phi, output_dim_phi, dropout_rate_phi = 0.0, activation_fun_phi = 'relu'):            # parameters of the phi MLP
                super().__init__()

                self.psi = tf.keras.Sequential[(
                    tf.keras.layers.Dense(hidden_dim1_psi, activation = activation_fun_psi),
                    tf.keras.layers.Dropout(dropout_rate_psi),
                    tf.keras.layers.Dense(hidden_dim2_psi, activation = activation_fun_psi),
                    tf.keras.layers.Dropout(dropout_rate_psi),
                    tf.keras.layers.Dense(output_dim_psi, activation = activation_fun_psi)
                )]

                if aggregation_function == 'mean':
                    self.aggregation_function = tf.reduce_mean
                elif aggregation_function == 'max':
                    self.aggregation_function = tf.reduce_max
                elif aggregation_function == 'sum':
                    self.aggregation_function = tf.reduce_sum
                else:
                    raise ValueError(f"Unknown aggregation function: {aggregation_function}")

                
                self.phi = tf.keras.Sequential[(
                    tf.keras.layers.Dense(hidden_dim1_phi, activation = activation_fun_phi),
                    tf.keras.layers.Dropout(dropout_rate_phi),
                    tf.keras.layers.Dense(hidden_dim2_phi, activation = activation_fun_phi),
                    tf.keras.layers.Dropout(dropout_rate_phi),
                    tf.keras.layers.Dense(output_dim_phi, activation = activation_fun_phi)
                )]

    def call(self, x, training = True):
        # x: (batch_size, set_size (number of reads), dimensions (number of features per read))
        psi_embedding = self.psi(x, training = training)

        # psi_embedding: (batch_size, set_size, output_dim_psi)
        psi_embedding_agg = self.aggregation_function(psi_embedding, axis = 1)

        # psi_embedding_agg: (batch_size, output_dim_psi)
        y = self.phi(psi_embedding_agg, training = training)

        return y


# how to use as part of the full model
'''
x_dummy = tf.random.normal((32,200,7))

model_dummy = FreeSample(
    input_dim = 7,
    hidden_dim1_psi = 128,
    hidden_dim2_psi = 64,
    output_dim_psi = 32,
    activation_psi = 'relu',
    dropout_rate_psi = 0.0,
    aggregation_function = 'mean',
    hidden_dim1_phi = 32,
    hidden_dim2_phi = 16,
    output_dim_phi = 8,
    activation_fun_phi = 'relu',
    dropout_rate_phi = 0.0
)


model_dummy.compile(
    optimizer = 'adam',  # just an example
    loss = 'mse'  # just an example
)

y = model_dummy(x_dummy, training = True)   # during inference add training = False
'''








# the free genomic context sub-model______________________________________________________________________________________________


# the learned positional encoding class

class GetSemanticAndLearnedPositionalEncodings(layers.Layer):

    '''
    genomic_context_size --> just the genomic window size around the variant of interest
    num_unique_bases ---> 5 (A,T,G,C,N)
    base_embedding_dim --> the size of the vector to which each individual base will be projected to
                           (In this specific case this is also the size of the vector of the learned positonal encoding)
    '''
    
    def __init__(self, genomic_context_size, num_unique_bases, base_embedding_dim, **kwargs):
        super()/__init__(**kwargs)

        # embedding for the semantic meaning
        self.semantic_embedding = layers.Embedding(
                                                    input_dim = num_unique_bases,
                                                    output_dim = base_embedding_dim
        )

        # learned positional encoding
        self.learned_positional_encoding = layers.Embedding(
                                                    input_dim = genomic_context_size,
                                                    output_dim = base_embedding_dim  # the dimensionality of the positional and semantic encodings are the same
        )
        
        self.genomic_context_size = genomic_context_size
        self.num_unique_bases = num_unique_bases
        self.base_embedding_dim = base_embedding_dim

    
    def call(self, x):
        length = tf.shape(x)[-1]
        positions = tf.range(start = 0, limit = length, delta = 1)  # a placeholder array for the positional embeddings
        base_semantic_embeddings = self.semantic_embedding(x)
        positional_encodings = self.learned_positional_encoding(x)
        return base_semantic_embeddings + positional_encodings
        '''
        Each base will therefore be represented by two (equally sized) vectors, one for the semantic meaning and one for the
        learned positional encoding.
        '''


class TransformerEncoder(layers.Layer):

    '''
    base_embedding_dim --> the size of the vector to which each individual base will be projected to
                           (In this specific case this is also the size of the vector of the learned positonal encoding)
    dense_dim --> the output dimension of the dense projection
    num_heads --> the number of the self attention heads
    '''

    def __init__(self, base_embedding_dim, dense_dim, num_heads, **kwargs):
        super().__init__(**kwargs)

        self.base_embedding_dim = base_embedding_dim
        self.dense_dim = dense_dim
        self.num_heads = num_heads

        self.attention = layers.MultiHeadAttention(
                                                    num_heads = num_heads,
                                                    key_dim = base_embedding_dim
        )

        self.dense_projection = tf.keras.Sequenctial([
                                                    tf.keras.layers.Dense(dense_dim, activation = 'relu'),
                                                    tf.keras.layers.Dense(base_embedding_dim)  # back to the original embedding dim
        ])

        self.layer_norm1 = tf.keras.layers.LayerNormalization()
        self.layer_norm2 = tf.keras.layers.LayerNormalization()


    def call(self,x):

        attention_output = self.attention(x,x)  # 'x,x' because its self attention
        normalization_output1 = self.layer_norm1(x,attention_output)  # the '+' is a residual connection
        projection_output = self.dense_projection(normalization_output1)
        multi_head_self_attention_block_output = self.layer_norm2(projection_output + normalization_output1)  # residual connection
        return multi_head_self_attention_block_output


# how to use as part of the full model
'''
num_unique_bases = 5  # (A,T,G,C,N)  # !!!!! CHECK WHAT IS THE MEANING OD THIS AFTER ONE HOT ENCODING
genomic_context_windowsize = 21
embed_dim = 256
num_heads = 4
dense_dim = 32
dropout_rate = 0.0

inputs = keras.Input(shape=(None,), dtype="int64")  # I DONT UNDERSTAND THIS !!!
x = GetSemanticAndLearnedPositionalEncodings(genomic_context_windowsize, num_unique_bases, embed_dim)(inputs)
x = TransformerEncoder(embed_dim, dense_dim, num_heads)
x = layers.GlobalMaxPooling1D()(x)  # *
context_aware_genomic_context = layers.Dropout(dropout_rate)(x)  # this is the output of the free genomic context sub-model

# * before this line, each token/base is represented by a vector, each base must be represented by a scalar, so that
# the entire sequence will be represented by a vector.
'''




# the final projection of the concatenated output_____________________________________________________________________________
class FinalProjection(tf.keras.Model):
    def __init__(self, input_dim, hidden_dim1, hidden_dim2, hidden_activation_fun = 'relu', dropout_rate = 0.0):
        super().__init__()
    

        self.MLP = tf.keras.Sequential([
            tf.keras.layers.Dense(hidden_dim1, activation = hidden_activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(hidden_dim2, activation = hidden_activation_fun),
            tf.keras.layers.Dropout(dropout_rate),
            tf.keras.layers.Dense(1, activation = 'sigmoid')  # the point estimator of the Binomial success probability parameter
        ])

    def call(self, x, training = True):
        # shape of x: (batch_size, input_dim (dimensionality of the concatenated embeddings))
        y = self.MLP(x, training = training)
        return y



# ===============================================================================================================================================


# ===================================================== defining the full model =================================================================


class Epsilon(tf.keras.Model):
    def __init__(self,
                 
                 # inputs for the instructive genomic submodel
                 intructive_genomic_input_dim,
                 instructive_genomic_hidden_dim1,
                 instructive_genomic_hidden_dim2,
                 instructive_genomic_output_dim,
                 instructive_genomic_activation_fun = 'relu',
                 intructive_genomic_dropout_rate = 0.0,


                 # inputs for the instructive sample-specific submodel
                 instructive_sample_input_dim,
                 instructive_sample_hidden_dim1,
                 instructive_sample_hidden_dim2,
                 instructive_sample_output_dim,
                 instructive_sample_activation_fun = 'relu',
                 instructive_sample_dropout_rate = 0.0,


                 # inputs for the free sample-specific submodel
                 free_sample_input_dim,
                 free_sample_hidden_dim1_psi,
                 free_sample_hidden_dim2_psi,
                 free_sample_output_dim_psi,
                 free_sample_dropout_rate_psi = 0.0,
                 free_sample_activation_fun_psi = 'relu',
                 free_sample_aggregation_function = 'mean',
                 free_sample_hidden_dim1_phi,
                 free_sample_hidden_dim2_phi,
                 free_sample_output_dim_phi,
                 free_sample_droupout_rate_phi = 0.0,
                 free_sample_activation_fun_phi = 'relu',


                 # inputs for the free genomic submodel
                 genomic_context_size,
                 num_unique_bases,
                 base_embedding_dim,
                 dense_dim,
                 num_heads
                 
                 
                 
                 final_proj_dim1,
                 final_proj_dim2,
                 final_proj_dropout_rate = 0.0,
                 final_proje_hidden_activation_fun = 'relu'
                 ):
                
                super().__init__()


                self.InstructiveGenomic = InstructiveGenomic(
                     input_dim = intructive_genomic_input_dim,
                     hidden_dim1 = instructive_genomic_hidden_dim2,
                     hidden_dim2 = instructive_genomic_hidden_dim2,
                     output_dim = instructive_genomic_output_dim,
                     activation_fun = instructive_genomic_activation_fun,
                     dropout_rate = intructive_genomic_dropout_rate
                )


                self.InstructiveSample = InstructiveSample(
                     input_dim = intructive_sample_input_dim,
                     hidden_dim1 = instructive_sample_hidden_dim2,
                     hidden_dim2 = instructive_sample_hidden_dim2,
                     output_dim = instructive_sample_output_dim,
                     activation_fun = instructive_sample_activation_fun,
                     dropout_rate = intructive_sample_dropout_rate
                )


                self.FreeSample = FreeSample(
                     input_dim = free_sample_input_dim,
                     hidden_dim1_psi = free_sample_hidden_dim1_psi,
                     hidden_dim2_psi = free_sample_hidden_dim2_psi,
                     output_dim_psi = free_sample_output_dim_psi,
                     dropout_rate_psi = free_sample_dropout_rate_psi,
                     activation_fun_psi = free_sample_activation_fun_psi,
                     aggregation_function = free_sample_aggregation_function,                                                                                     # the aggregation function
                     hidden_dim1_phi = free_sample_hidden_dim1_phi,
                     hidden_dim2_phi = free_sample_hidden_dim2_phi,
                     output_dim_phi = free_sample_output_dim_phi,
                     dropout_rate_phi = free_sample_dropout_rate_phi,
                     activation_fun_phi = free_sample_activation_fun_phi
                )


                self.FinalProjection = FinalProjection(
                     input_dim = (instructive_genomic_output_dim + instructive_sample_output_dim + free_sample_input_dim + base_embedding_dim),
                     hidden_dim1 = final_proj_dim1,
                     hidden_dim2 = final_proj_dim2,
                     dropout_rate = final_proj_dropout_rate,
                     hidden_activation_function = final_proje_hidden_activation_fun)        

    def call(self,
            instructive_genomic_input, # input for the instructive genomic submodel
            instructive_sample_input,   # input for the instructive sample-specific submodel
            free_sample_input,   # input for the free sample-specific submodel
            free_genomic_input,   # input for the free genomic submodel
            training=True):
         

         intructive_genomic_output = self.InstructiveGenomic(x = instructive_genomic_input, training = training)
         instructive_sample_output = self.InstructiveSample(x = instructive_sample_input, training = training)
         free_sample_output = self.FreeSample(x = free_sample_input, training = training)

         embeddings_for_transformer = GetSemanticAndLearnedPositionalEncodings(
              genomic_context_size=self.genomic_context_size,
              num_unique_bases=self.num_unique_bases,
              base_embedding_dim=self.base_embedding_dim)(free_genomic_input)
         
         free_genomic_context_output = TransformerEncoder(
              base_embedding_dim=self.base_embedding_dim,
              dense_dim = self.dense_dim,
              num_heads = self.num_heads)(embeddings_for_transformer)
         

         free_genomic_context_output = tf.reduce_max(free_genomic_context_output, axis=1)

         concatenated = tf.concat([
              intructive_genomic_output,
              instructive_sample_output,
              free_sample_output,
              free_genomic_context_output], axis=-1)
         

         binomial_estimate = self.FinalProjection(x = concatenated, training = training)

         return binomial_estimate
    


         


         


