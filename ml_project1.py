import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt

# --- Step 1: Extract features from SDF files ---
def extract_features_from_sdf(folder_path, label):
    data = []
    for file in os.listdir(folder_path):
        if file.endswith('.sdf'):
            sdf = Chem.SDMolSupplier(os.path.join(folder_path, file))
            for mol in sdf:
                if mol is None:
                    continue
                # Descriptors
                desc = {
                    'MolWt': Descriptors.MolWt(mol),
                    'LogP': Descriptors.MolLogP(mol),
                    'NumHDonors': Descriptors.NumHDonors(mol),
                    'NumHAcceptors': Descriptors.NumHAcceptors(mol),
                    'TPSA': Descriptors.TPSA(mol),
                    'NumRotatableBonds': Descriptors.NumRotatableBonds(mol)
                }
                # Fingerprint
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fp_bits = np.array(list(fp.ToBitString())).astype(int)
                for i, bit in enumerate(fp_bits):
                    desc[f'FP_{i}'] = bit
                desc['Label'] = label
                data.append(desc)
    return pd.DataFrame(data)

# --- Step 2: Load, label, and combine data ---
drugs_df = extract_features_from_sdf('drugs_sdf', 1)
non_drugs_df = extract_features_from_sdf('nondrugs_sdf', 0)
df = pd.concat([drugs_df, non_drugs_df], ignore_index=True).fillna(0)
df.to_csv('ligand_dataset_with_fingerprints.csv', index=False)

print("Sample of the dataset:")
print(df.head())

# --- Step 3: Train-test split & model training ---
X = df.drop(columns=['Label'])
y = df['Label']
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)
model = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)
y_pred = model.predict(X_test)

# --- Step 4: Model evaluation ---
print("\nModel Evaluation")
print(f"Accuracy: {accuracy_score(y_test, y_pred) * 100:.2f}%")
print("\nClassification Report:\n", classification_report(y_test, y_pred, target_names=['Non-Drug', 'Drug-Like']))

# confusion matrix
cm = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(6, 5))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
            xticklabels=['Predicted Non-Drug', 'Predicted Drug-Like'],
            yticklabels=['Actual Non-Drug', 'Actual Drug-Like'])
plt.title('Confusion Matrix')
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.tight_layout()
plt.show()

# --- Step 5: Feature importance ---
importances = model.feature_importances_
top_features = pd.DataFrame({
    'Feature': X.columns,
    'Importance': importances
}).sort_values(by='Importance', ascending=False).head(20)

plt.figure(figsize=(8, 6))
sns.barplot(x='Importance', y='Feature', data=top_features)
plt.title('Top 20 Feature Importances')
plt.tight_layout()
plt.show()

# --- Step 6: Cross-validation ---
cv_scores = cross_val_score(model, X, y, cv=5)
print("\nCross-Validation Scores:", cv_scores)
print(f"Average CV Accuracy: {np.mean(cv_scores) * 100:.2f}%")