import tensorflow as tf
import numpy as np
from tensorflow.keras.layers import Input, Dense, Dropout, Embedding, Flatten, LayerNormalization
from tensorflow.keras.callbacks import EarlyStopping, History, ModelCheckpoint
from tensorflow.keras.utils import plot_model
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.initializers import VarianceScaling
from tensorflow.keras.optimizers import Adam, AdamW, Nadam






class GetSemanticAndLearnedPositionalEncodings(tf.keras.layers.Layer):

    '''
    genomic_context_size --> just the genomic window size around the variant of interest
    num_unique_bases ---> 5 (A,T,G,C,N)
    base_embedding_dim --> the size of the vector to which each individual base will be projected to
                           (In this specific case this is also the size of the vector of the learned positonal encoding)
    '''
    
    def __init__(self, genomic_context_size, num_unique_bases, base_embedding_dim, **kwargs):
        super().__init__(**kwargs)

        # embedding for the semantic meaning
        self.semantic_embedding = tf.keras.layers.Embedding(
                                                    input_dim = num_unique_bases,
                                                    output_dim = base_embedding_dim
        )

        # learned positional encoding
        self.learned_positional_encoding = tf.keras.layers.Embedding(
                                                    input_dim = genomic_context_size,
                                                    output_dim = base_embedding_dim  # the dimensionality of the positional and semantic encodings are the same
        )
        
        self.genomic_context_size = genomic_context_size
        self.num_unique_bases = num_unique_bases
        self.base_embedding_dim = base_embedding_dim

    
    def call(self, x):
        positions = tf.range(start = 0, limit = tf.shape(x)[-1], delta = 1)  # a placeholder array for the positional embeddings
        positions = tf.expand_dims(positions, axis=0)
        base_semantic_embeddings = self.semantic_embedding(x)
        positional_encodings = self.learned_positional_encoding(positions)
        return base_semantic_embeddings + positional_encodings
        '''
        Each base will therefore be represented by two (equally sized) vectors, one for the semantic meaning and one for the
        learned positional encoding.
        '''


class TransformerEncoder(tf.keras.layers.Layer):

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

        self.attention = tf.keras.layers.MultiHeadAttention(
                                                    num_heads = num_heads,
                                                    key_dim = base_embedding_dim
        )

        self.dense_projection = tf.keras.Sequential([
                                                    tf.keras.layers.Dense(dense_dim, activation = 'relu'),
                                                    tf.keras.layers.Dense(base_embedding_dim)  # back to the original embedding dim
        ])

        self.layer_norm1 = tf.keras.layers.LayerNormalization()
        self.layer_norm2 = tf.keras.layers.LayerNormalization()


    def call(self,x):

        attention_output = self.attention(x,x)  # 'x,x' because its self attention
        normalization_output1 = self.layer_norm1(x + attention_output)  # the '+' is a residual connection
        projection_output = self.dense_projection(normalization_output1)
        multi_head_self_attention_block_output = self.layer_norm2(projection_output + normalization_output1)  # residual connection
        return multi_head_self_attention_block_output






class FreeSample(tf.keras.Model):
    """
    ===========================================  Notes  ============================================================

    - the coverage (number of reads) is different for each sample and this is a bit of a problem for the model,
      becasue it needs a fixed sized input
    - so I do a padding up to 10K reads
    - I embedd the padded input
    - then during the aggregation I mask out the padded embeddings
    - then I proceed as normal
    """
    def __init__(self,
                hidden_dim1_psi, hidden_dim2_psi, output_dim_psi,  # parameters of the psi MLP's
                                                                                                     # the aggregation function
                hidden_dim1_phi, hidden_dim2_phi, output_dim_phi,activation_fun_psi = 'relu',
                aggregation_function = 'mean', dropout_rate_psi = 0.0,
                dropout_rate_phi = 0.0, activation_fun_phi = 'relu',
                ):

                super().__init__()

                self.psi = tf.keras.Sequential([
                    tf.keras.layers.Dense(hidden_dim1_psi, activation = activation_fun_psi),
                    tf.keras.layers.Dropout(dropout_rate_psi),
                    tf.keras.layers.Dense(hidden_dim2_psi, activation = activation_fun_psi),
                    tf.keras.layers.Dropout(dropout_rate_psi),
                    tf.keras.layers.Dense(output_dim_psi, activation = activation_fun_psi)
                ])

                self.aggregation_function = aggregation_function

                
                self.phi = tf.keras.Sequential([
                    tf.keras.layers.Dense(hidden_dim1_phi, activation = activation_fun_phi),
                    tf.keras.layers.Dropout(dropout_rate_phi),
                    tf.keras.layers.Dense(hidden_dim2_phi, activation = activation_fun_phi),
                    tf.keras.layers.Dropout(dropout_rate_phi),
                    tf.keras.layers.Dense(output_dim_phi, activation = activation_fun_phi)
                ])

    def call(self, x_padded, padding_mask, training = True):
        # x: (batch_size, set_size (number of reads), dimensions (number of features per read))
        psi_embedding = self.psi(x_padded, training = training)

        # expand the padding mask to match the dimensions and cast it to a dtype that is compatible with psi_embedding
        padding_mask = tf.expand_dims(padding_mask, axis = -1)
        padding_mask = tf.cast(padding_mask, psi_embedding.dtype)

        # multiply the psi_embedding by the padding mask
        psi_embedding_masked = psi_embedding * padding_mask

        # psi_embedding: (batch_size, set_size, output_dim_psi)
        # mask out the embeddings that come from the padded regions
        if self.aggregation_function == 'mean':
            sum_over_reads = tf.reduce_sum(psi_embedding_masked, axis=1)
            num_real_reads = tf.reduce_sum(padding_mask, axis=1) + 1e-8      # 1e-8 is a pseudcount to not get division by zero error
            psi_embedding_masked_agg = sum_over_reads / num_real_reads

        elif self.aggregation_function == 'max':
            neg_inf = tf.constant(-1e9, dtype=psi_embedding.dtype)
            psi_embedding_masked = tf.where(padding_mask > 0, psi_embedding_masked, neg_inf)
            psi_embedding_masked_agg = tf.reduce_max(psi_embedding_masked, axis=1)

        elif self.aggregation_function == 'sum':
            psi_embedding_masked_agg = tf.reduce_sum(psi_embedding_masked, axis=1)

        else:
             print('Wrong value for the aggregation function parameter.')


        # psi_embedding_agg: (batch_size, output_dim_psi)
        y = self.phi(psi_embedding_masked_agg, training = training)

        return y





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
    



class FinalProjection(tf.keras.Model):
    def __init__(self, hidden_dim1, hidden_dim2, hidden_activation_fun = 'relu', dropout_rate = 0.0):
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


class Epsilon(tf.keras.Model):
    def __init__(self,
                 
                 # inputs for the instructive genomic submodel
                 intructive_genomic_input_dim,
                 instructive_genomic_hidden_dim1,
                 instructive_genomic_hidden_dim2,
                 instructive_genomic_output_dim,


                 # inputs for the instructive sample-specific submodel
                 instructive_sample_input_dim,
                 instructive_sample_hidden_dim1,
                 instructive_sample_hidden_dim2,
                 instructive_sample_output_dim,


                 # inputs for the free sample-specific submodel
                 free_sample_hidden_dim1_psi,
                 free_sample_hidden_dim2_psi,
                 free_sample_output_dim_psi,
                 free_sample_hidden_dim1_phi,
                 free_sample_hidden_dim2_phi,
                 free_sample_output_dim_phi,


                 # inputs for the free genomic submodel
                 genomic_context_size,
                 num_unique_bases,
                 base_embedding_dim,
                 dense_dim,
                 num_heads,
                 
                 
                 # inputs for the final projection
                 final_proj_dim1,
                 final_proj_dim2,


                 # the parameters with default values must be last
                 instructive_genomic_activation_fun = 'relu',
                 intructive_genomic_dropout_rate = 0.0,
                 instructive_sample_activation_fun = 'relu',
                 instructive_sample_dropout_rate = 0.0,
                 free_sample_dropout_rate_psi = 0.0,
                 free_sample_activation_fun_psi = 'relu',
                 free_sample_aggregation_function = 'mean',
                 free_sample_dropout_rate_phi = 0.0,
                 free_sample_activation_fun_phi = 'relu',
                 final_proj_dropout_rate = 0.0,
                 final_projection_hidden_activation_fun = 'relu'
                 ):
                
                super().__init__()


                self.InstructiveGenomic = InstructiveGenomic(
                     input_dim = intructive_genomic_input_dim,
                     hidden_dim1 = instructive_genomic_hidden_dim1,
                     hidden_dim2 = instructive_genomic_hidden_dim2,
                     output_dim = instructive_genomic_output_dim,
                     activation_fun = instructive_genomic_activation_fun,
                     dropout_rate = intructive_genomic_dropout_rate
                )


                self.InstructiveSample = InstructiveSample(
                     input_dim = instructive_sample_input_dim,
                     hidden_dim1 = instructive_sample_hidden_dim1,
                     hidden_dim2 = instructive_sample_hidden_dim2,
                     output_dim = instructive_sample_output_dim,
                     activation_fun = instructive_sample_activation_fun,
                     dropout_rate = instructive_sample_dropout_rate
                )


                self.FreeSample = FreeSample(
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


                self.embedding_layer = GetSemanticAndLearnedPositionalEncodings(
                    genomic_context_size=genomic_context_size,
                    num_unique_bases=num_unique_bases,
                    base_embedding_dim=base_embedding_dim
                )

                self.transformer_encoder = TransformerEncoder(
                    base_embedding_dim=base_embedding_dim,
                    dense_dim=dense_dim,
                    num_heads=num_heads
                )


                self.FinalProjection = FinalProjection(
                     hidden_dim1 = final_proj_dim1,
                     hidden_dim2 = final_proj_dim2,
                     dropout_rate = final_proj_dropout_rate,
                     hidden_activation_fun = final_projection_hidden_activation_fun)


    def call(self,
            instructive_genomic_input, # input for the instructive genomic submodel
            instructive_sample_input,   # input for the instructive sample-specific submodel
            free_sample_padded_input,   # input for the free sample-specific submodel
            free_genomic_input,   # input for the free genomic submodel
            free_sample_padding_mask,
            training=True):
         

         intructive_genomic_output = self.InstructiveGenomic(x = instructive_genomic_input, training = training)
         instructive_sample_output = self.InstructiveSample(x = instructive_sample_input, training = training)
         free_sample_output = self.FreeSample(x_padded = free_sample_padded_input, padding_mask = free_sample_padding_mask, training = training)         
         embeddings_for_transformer = self.embedding_layer(free_genomic_input)
         free_genomic_context_output = self.transformer_encoder(embeddings_for_transformer)
         free_genomic_context_output = tf.reduce_max(free_genomic_context_output, axis=1)

         concatenated = tf.concat([
              intructive_genomic_output,
              instructive_sample_output,
              free_sample_output,
              free_genomic_context_output], axis=-1)
         

         binomial_estimate = self.FinalProjection(x = concatenated, training = training)

         return binomial_estimate
    




# ========================================================= the loss function =======================================================================


class BinomialNegLogLoss(tf.keras.losses.Loss):
    def __init__(self, pseudocount=1e-8, normalize_by_coverage=False, name="binomial_nll"):
        super().__init__(name=name)
        self.eps = pseudocount
        self.normalize = normalize_by_coverage

    def call(self, y_true, y_pred):
        '''
        y_pred ---> the error rate predicted by Epsilon
        n ---> the coverage (from data)
        k ---> the number of alternative reads (from data)

        normalize_by_coverage = False
                --> samples with higher coverage contribute more to the loss
        
        normalize_by_coverage = True
                --> the coverage has no effect on the contribution to the loss
        '''

        # make it robust for weird samples
        n = tf.maximum(y_true[:, 0], 0.0)
        k = tf.maximum(y_true[:, 1], 0.0)
        k = tf.minimum(k, n)

        y_pred = tf.squeeze(y_pred, axis=-1)
        y_pred = tf.clip_by_value(y_pred, self.eps, 1.0 - self.eps)

        log_likelihood = k * tf.math.log(y_pred) + (n - k) * tf.math.log(1 - y_pred)

        if self.normalize:
            log_likelihood = log_likelihood / (n + self.eps)

        log_likelihood = tf.clip_by_value(log_likelihood, -1e6, 1e6) # clip extreme values to stabilize the loss

        return -tf.reduce_mean(log_likelihood)

        log_likelihood = k * tf.math.log(y_pred) + (n - k) * tf.math.log(1 - y_pred)

        if self.normalize:
            log_likelihood = log_likelihood / (n + self.eps)

        return -tf.reduce_mean(log_likelihood)
