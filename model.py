from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def train_model(feature_df):
    feature_df.fillna(0, inplace=True)
    feature_df['label'] = feature_df['expression'].apply(lambda x: 'High' if x >= 1 else 'Low')

    X = feature_df.drop(columns=['gene', 'expression', 'label'])
    y = feature_df['label']
    genes = feature_df['gene']

    # Split data
    X_train, X_test, y_train, y_test, genes_train, genes_test = train_test_split(
        X, y, genes, test_size=0.2, random_state=42
    )

    # Train classification model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    # Predict
    y_pred = model.predict(X_test)

    # Evaluate
    acc = accuracy_score(y_test, y_pred)
    print(f"Test Accuracy: {acc:.4f}")
    print("\nClassification Report:\n", classification_report(y_test, y_pred))

    # Confusion Matrix Plot
    cm = confusion_matrix(y_test, y_pred, labels=['High', 'Low'])

    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['High', 'Low'],
                yticklabels=['High', 'Low'])

    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.savefig("confusion_matrix.png")
    plt.show()

    # Combine predictions with genes
    test_features = X_test.copy()
    test_features["gene"] = genes_test.values
    test_features["actual_label"] = y_test.values
    test_features["predicted_label"] = y_pred

    # Reorder columns
    cols = ["gene", "actual_label", "predicted_label"] + [
        col for col in test_features.columns if col not in ["gene", "actual_label", "predicted_label"]
    ]
    test_features = test_features[cols]

    # Save predictions
    test_features.to_csv("predicted_gene_activity.csv", index=False)
    print("Predicted labels saved to predicted_gene_activity.csv")

    # Show sample
    #print(test_features.head())

    return model
