from keras import Sequential, Input
from keras.layers import Dense, Dropout, Input
from keras.optimizers import Adam

def build_ann(input_dim):
    """
    Build a deeper and more robust ANN model.
    """
    model = Sequential([
        Input(shape=(input_dim,)),
        
        Dense(64, activation='relu'), # first layer more neurons
        Dropout(0.3),
        Dense(32, activation='relu'), 
        Dropout(0.2),
        Dense(16, activation='relu'), # extra layer
        Dropout(0.1),
        Dense(1, activation='sigmoid')
    ])
    
    optimizer = Adam(learning_rate=0.0001)
    
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])
    return model