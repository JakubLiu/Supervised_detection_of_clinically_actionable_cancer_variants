import tensorflow as tf
from tensorflow.keras import layers



'''
Takes as input a single read that is represented by a (d,900) 2D tensor, where d is the variable read length.
It outputs a 1D tensor of fixed size (projection_dim,).
These 1D embeddings then go into the Set Transformer.
'''
class CNN_read_embedder(tf.keras.Model):
    def __init__(self, projection_dim):
        super(CNN_read_embedder, self).__init__()

        self.conv32 = layers.Conv1D(filters=32, kernel_size=3, padding="same", activation="relu") # padding 'same' padds with zeros, this is important, because I dont want to loose the ends of the reads
        self.conv64 = layers.Conv1D(filters=64, kernel_size=3, padding="same", activation="relu")
        self.conv128 = layers.Conv1D(filters=128, kernel_size=3, padding="same", activation="relu")
        self.conv256 = layers.Conv1D(filters=256, kernel_size=3, padding="same", activation="relu")
        self.local_mean_pool = layers.AveragePooling1D(pool_size = 4, strides = 4)
        self.global_mean_pool = layers.GlobalAveragePooling1D()
        self.projection = layers.Dense(projection_dim, activation='relu')

    def call(self, x):
        #x = tf.expand_dims(x, axis=0)  # manually 'add' the batch dimension, because the tensorflow CNN expects it (although in EpsilonTO one_locus = one_batch)
        x = self.conv32(x)
        x = self.conv64(x)
        x = self.conv128(x)
        x = self.local_mean_pool(x)
        x = self.conv256(x)
        x = self.global_mean_pool(x)
        x = self.projection(x)
        #x = tf.squeeze(x, axis=0)   # explicitly squeeze to one dimension

        return x



'''
Takes as input a pileup (set) of one dimensional read embeddings. So, the set has a shape of (coverage, projection_dim).
Where projection_dim is the size of the 1D embedding returned by CNN_read_embedder.
It does all the standard Transformer stuff but with no positional encoding, because I dont give a fuck about the order of the reads
in the pileup, because it carries no biological meaning.
It outputs a 1-dimensional representation of the read pileup.
'''
class ReadSetTransformer(tf.keras.Model):

    def __init__(self, read_embedding_dim, num_attention_heads):
        super().__init__()

        self.MHA = layers.MultiHeadAttention(num_heads=num_attention_heads, key_dim=read_embedding_dim//num_attention_heads)
        self.projection = tf.keras.Sequential([layers.Dense(2*read_embedding_dim, activation='relu'), layers.Dense(read_embedding_dim)])
        self.layer_normalization1 = layers.LayerNormalization()
        self.layer_normalization2 = layers.LayerNormalization()

    def call(self,x):

        MHA_output = self.MHA(x,x) # multi head attention
        resid_norm_output = self.layer_normalization1(x + MHA_output) # first residual connection + layer normalization
        projection_output = self.projection(resid_norm_output)   # MLP embedding
        final_output = self.layer_normalization2(resid_norm_output + projection_output)  # second residual connection + layer normalization
        return final_output


'''
This takes as input teh 2-dimensional genomic context tensor (the first dimension is the size of the genomic window and the 2nd dimension is the one-hot encoding of the bases (4)).
It returns a 1-dimensional tensor that encodes both semantic and positional meaning.
'''
class ContextEmbeddingsAndPositionalEncodings(layers.Layer):

    def __init__(self, genomic_context_size, base_embedding_dim, **kwargs):
        super().__init__(**kwargs)

        self.semantic_embedding = layers.Dense(base_embedding_dim)  # project the one-hot encoded bases into an embedding of higher dimensionality
        self.learned_positional_encoding = layers.Embedding(input_dim=genomic_context_size, output_dim=base_embedding_dim)
        self.genomic_context_size = genomic_context_size

    def call(self,x):

        positions = tf.range(start = 0, limit = tf.shape(x)[1], delta = 1)  # a placeholder array for the positional embeddings tf.shape(x) = (batch_size, genomic_context_size, one-hot-encoding_dim)
        positions = tf.expand_dims(positions, axis=0)
        base_embeddings = self.semantic_embedding(x)
        positional_encodings = self.learned_positional_encoding(positions)
        output_embedding = base_embeddings + positional_encodings

        return output_embedding


'''
This takes as input the embedding of the genomic context that encapsulates positional meaning.
Then it does the stnadard Transformer stuff.
It returns a 1-dimensional tensor.
'''
class ContextTransformer(layers.Layer):

    def __init__(self, base_embedding_dim, dense_dim, num_attention_heads, **kwargs):
        super().__init__(**kwargs)

        self.base_embedding_dim = base_embedding_dim
        self.dense_dim = dense_dim
        self.num_heads = num_attention_heads

        self.MHA = layers.MultiHeadAttention(num_heads=num_attention_heads, key_dim=base_embedding_dim//num_attention_heads)
        self.projection = tf.keras.Sequential([layers.Dense(2*base_embedding_dim, activation='relu'), layers.Dense(base_embedding_dim)])
        self.layer_norm1 = tf.keras.layers.LayerNormalization()
        self.layer_norm2 = tf.keras.layers.LayerNormalization()

    def call(self,x):

        MHA_output = self.MHA(x,x)
        resid_norm_output1 = self.layer_norm1(x + MHA_output)
        projection_output = self.projection(resid_norm_output1)
        output_embedding = self.layer_norm2(resid_norm_output1 + projection_output)

        return output_embedding


'''
Takes as input the concatenated tensor that is the concatenation of the context and pileup tensors.
It returns a scalar error rate estimatie [0,1].
'''
class FinalProjection(tf.keras.Model):

    def __init__(self):
        super().__init__()

        self.MLP = tf.keras.Sequential([
            layers.Dense(128, activation='relu'),
            layers.Dropout(0.2),
            layers.Dense(64, activation='relu'),
            layers.Dropout(0.2),
            layers.Dense(32, activation='relu'),
            layers.Dropout(0.2),
            layers.Dense(1, activation='sigmoid')  # the error rate
        ])

    def call(self,x):
        error_rate = self.MLP(x)
        return error_rate



class EpsilonTO(tf.keras.Model):

    def __init__(self,
                 _projection_dim,  # the dimension of the read tensor (from the CNN)
                 _num_attention_heads_set,  # for the Set Transformer
                 _genomic_context_size,  # for the positional and semantic embeddings
                 _base_embedding_dim,    # for the positional and semantic embeddings
                 _dense_dim,          # for the context Transformer
                 _num_attention_heads_context,   # for the context Transformer
                 **kwargs
                 ):

        super().__init__(**kwargs)

        self.config = {"_projection_dim": _projection_dim, "_num_attention_heads_set": _num_attention_heads_set, "_genomic_context_size": _genomic_context_size, "_base_embedding_dim": _base_embedding_dim, "_dense_dim": _dense_dim, "_num_attention_heads_context": _num_attention_heads_context}
        self.cnn_read_embedder = CNN_read_embedder(projection_dim=_projection_dim)
        self.read_set_transformer = ReadSetTransformer(read_embedding_dim=_projection_dim, num_attention_heads=_num_attention_heads_set)
        self.pileup_pool = layers.GlobalAveragePooling1D()
        self.context_transformer_encodings = ContextEmbeddingsAndPositionalEncodings(genomic_context_size=_genomic_context_size, base_embedding_dim=_base_embedding_dim)
        self.context_transformer = ContextTransformer(base_embedding_dim=_base_embedding_dim, dense_dim=_dense_dim, num_attention_heads=_num_attention_heads_context)
        self.context_pool = layers.GlobalAveragePooling1D()
        self.final_projection = FinalProjection()


    def call(self, pileup_input, genomic_context_input):

        #single_read_embeddings = tf.map_fn(
        #    self.cnn_read_embedder,
        #    pileup_input,
        #    fn_output_signature=tf.float32
        #)

        #if len(pileup_input.shape) == 4:
        #    pileup_input = tf.squeeze(pileup_input, axis=0)

        single_read_embeddings = self.cnn_read_embedder(pileup_input)

        #single_read_embeddings = tf.expand_dims(
        #    single_read_embeddings,
        #    axis=0
        #)


        pileup_set_embedding = self.read_set_transformer(
           tf.expand_dims(single_read_embeddings, axis=0)
        )

        #pileup_set_embedding = self.read_set_transformer(
        #    single_read_embeddings
        #)


        pooled_pileup_set_embedding = self.pileup_pool(pileup_set_embedding)


        genomic_context_input = tf.expand_dims(
            genomic_context_input,
            axis=0
        )


        context_encodings = self.context_transformer_encodings(genomic_context_input)
        context_embedding = self.context_transformer(context_encodings)
        pooled_context_embedding = self.context_pool(context_embedding)
        concatenation = tf.concat([pooled_pileup_set_embedding, pooled_context_embedding], axis = -1)
        y = self.final_projection(concatenation)

        return y  # the error rate estimate

    def get_config(self):
        config = super().get_config()
        config.update(self.config)
        return config
    
    @classmethod
    def from_config(cls, config):
        return cls(**config)


