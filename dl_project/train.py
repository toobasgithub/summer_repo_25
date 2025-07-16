import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import matplotlib.pyplot as plt
from sklearn.utils import class_weight
from model import build_ann
from keras import callbacks

def load_data(csv_file):
    df = pd.read_csv(csv_file, low_memory=False) 
    df['label'] = df['label'].map({'deleterious': 1, 'neutral': 0})
    X = df.drop(['chrom', 'pos', 'label'], axis=1).values
    y = df['label'].values
    return X, y

def plot_confusion(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(4,4))
    plt.imshow(cm, cmap='Blues')
    plt.title('Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.xticks([0,1], ['neutral', 'deleterious'])
    plt.yticks([0,1], ['neutral', 'deleterious'])
    for i in range(2):
        for j in range(2):
            plt.text(j, i, cm[i, j], ha='center', va='center', color='red')
    plt.tight_layout()
    plt.savefig('confusion_matrix.png')
    plt.close()

def plot_history(history):
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.title('Model Accuracy')
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Validation'], loc='upper left')
    plt.subplot(1, 2, 2)
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('Model Loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Validation'], loc='upper left')
    plt.tight_layout()
    plt.savefig('training_history.png')
    plt.close()

def run_training_pipeline(csv_file):
    X, y = load_data(csv_file)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    print("Calculating class weights...")
    class_weights = class_weight.compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
    class_weights = dict(enumerate(class_weights))
    print(f"Weights calculated for classes: {class_weights}")

    model = build_ann(X.shape[1])

    # 1. Early Stopping
    early_stopping = callbacks.EarlyStopping(monitor='val_loss', patience=3, restore_best_weights=True)
    # 2. Reduce Learning Rate
    reduce_lr = callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=2, min_lr=0.00001)
    
    print("Starting Advanced Model Training")
    
    history = model.fit(
        X_train, y_train, 
        epochs=10, 
        batch_size=32, 
        validation_split=0.1, 
        class_weight=class_weights, 
        callbacks=[early_stopping, reduce_lr], 
        verbose=1
    )
    
    print("Model Training Complete.")
    y_pred_prob = model.predict(X_test).flatten()
    y_pred = (y_pred_prob > 0.5).astype(int)
    
    acc = accuracy_score(y_test, y_pred)
    prec = precision_score(y_test, y_pred)
    rec = recall_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    
    print("Final Evaluation Metrics")
    print(f"Accuracy: {acc:.4f}")
    print(f"Precision: {prec:.4f}")
    print(f"Recall: {rec:.4f}")
    print(f"F1-Score: {f1:.4f}")
    
    plot_confusion(y_test, y_pred)
    print("Confusion matrix saved as 'confusion_matrix.png'")
    plot_history(history) 
    print("Training history saved as 'training_history.png'")