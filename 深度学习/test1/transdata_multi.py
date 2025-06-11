import os
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
from functools import partial
from itertools import combinations
from tqdm import tqdm
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
from concurrent.futures import ThreadPoolExecutor
import optuna
import random
import logging
import shap
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


# Configure logging
logging.basicConfig(
    filename="run_log.txt",
    level=logging.INFO,
    format="%(asctime)s - %(message)s"
)

# 正确初始化后，日志记录
logging.info("Logging system initialized.")

def log_message(message):
    logging.info(message)
    tqdm.write(message)  # Show in tqdm progress bar if visible
    
# Set random seed for reproducibility
def set_random_seed(seed=42):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# Load datasets and merge them
def load_and_merge_datasets(folder_path):
    merged_data, feature_sources, target = [], [], None
    for file_name in filter(lambda f: f.endswith(".csv"), os.listdir(folder_path)):
        file_path = os.path.join(folder_path, file_name)
        data = pd.read_csv(file_path)
        if target is None:
            target = data.iloc[:, 0].values
        else:
            assert np.array_equal(target, data.iloc[:, 0].values), "Inconsistent targets between files."
        feature_sources.extend([(col, file_name) for col in data.columns[1:]])
        merged_data.append(data.iloc[:, 1:].values)
    log_message(f"Loaded and merged datasets from {folder_path}.")
    return np.hstack(merged_data), target.astype(np.int64), feature_sources

# Standardize features and split data
def preprocess_data(X, y, test_size=0.2, random_state=42):
    X = StandardScaler().fit_transform(X)
    return train_test_split(X, y, test_size=test_size, random_state=random_state)

# PyTorch Dataset
class GeneExpressionDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.tensor(X, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.long)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

# Transformer-based Classifier
class TransformerClassifier(nn.Module):
    def __init__(self, input_size, num_classes, num_heads=4, num_layers=2, hidden_size=64, dropout_rate=0.1):
        super().__init__()
        self.input_proj = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Dropout(dropout_rate)
        )
        self.position_embeddings = nn.Parameter(torch.randn(1, 1, hidden_size))
        self.encoder_layers = nn.ModuleList([
            nn.ModuleDict({
                'self_attention': nn.MultiheadAttention(hidden_size, num_heads, dropout=dropout_rate, batch_first=True),
                'norm1': nn.LayerNorm(hidden_size),
                'ffn': nn.Sequential(
                    nn.Linear(hidden_size, hidden_size * 4),
                    nn.GELU(),
                    nn.Dropout(dropout_rate),
                    nn.Linear(hidden_size * 4, hidden_size),
                ),
                'norm2': nn.LayerNorm(hidden_size)
            }) for _ in range(num_layers)
        ])
        self.classifier = nn.Sequential(
            nn.Linear(hidden_size, num_classes)
        )

    def forward(self, x):
        x = self.input_proj(x).unsqueeze(1) + self.position_embeddings
        for layer in self.encoder_layers:
            attn_out, _ = layer['self_attention'](x, x, x)
            x = layer['norm1'](x + attn_out)
            ff_out = layer['ffn'](x)
            x = layer['norm2'](x + ff_out)
        return self.classifier(x.squeeze(1))

# Train and evaluate the model
def train_and_evaluate(model, train_loader, val_loader, device, max_epochs=3, lr=0.001):
    model = model.to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.AdamW(model.parameters(), lr=lr)
    model.train()
    for _ in range(max_epochs):
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            loss = criterion(model(X_batch), y_batch)
            loss.backward()
            optimizer.step()
    return evaluate_model(model, val_loader, device)

def evaluate_model(model, loader, device):
    model.eval()
    y_true, y_pred, y_prob = [], [], []
    with torch.no_grad():
        for X_batch, y_batch in loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            outputs = model(X_batch)
            probs = torch.softmax(outputs, dim=1)
            y_true.extend(y_batch.cpu().tolist())
            y_pred.extend(torch.argmax(probs, dim=1).cpu().tolist())
            y_prob.extend(probs[:, 1].cpu().tolist())
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Precision": precision_score(y_true, y_pred, zero_division=0),
        "Recall": recall_score(y_true, y_pred, zero_division=0),
        "F1 Score": f1_score(y_true, y_pred, zero_division=0),
        "ROC AUC": roc_auc_score(y_true, y_prob) if len(set(y_true)) > 1 else 0.0
    }

# Compute variable importance
def compute_variable_importance(model, feature_sources, selected_indices):
    """
    Computes variable importance for each selected feature.
    """
    weights = model.input_proj[0].weight.detach().cpu().numpy()
    importances = np.mean(np.abs(weights), axis=0)
    return pd.DataFrame([
        {"source": feature_sources[idx][1], "variable": feature_sources[idx][0], "importance": importance}
        for idx, importance in zip(selected_indices, importances)
    ])

# Save variable importance to CSV
def save_variable_importance(importance_df, file_name="variable_importance.csv"):
    """
    Save variable importance to a CSV file.
    """
    # Aggregate importance by summing up values for the same source and variable
    aggregated_df = importance_df.groupby(['source', 'variable'], as_index=False)['importance'].sum()

    # Save the aggregated data to a CSV file
    aggregated_df.to_csv(file_name, index=False)

    print(f"Variable importance saved to '{file_name}'.")
    logging.info(f"Aggregated variable importance saved to '{file_name}'.")

# SHAP computation and visualization
def calculate_shap_importance(model, X, feature_sources, file_name="shap_importance.csv"):
    """
    Calculate and save SHAP values for feature importance.
    Generate a combined SHAP plot with a bee swarm and a bar chart, saved as a high-resolution PDF.
    """
    log_message("Calculating SHAP values for feature importance...")
    model.eval()

    # Wrapper for SHAP compatibility with PyTorch
    def model_forward(inputs):
        with torch.no_grad():
            inputs = torch.tensor(inputs, dtype=torch.float32).to(next(model.parameters()).device)
            return model(inputs).cpu().numpy()

    # Initialize SHAP explainer
    explainer = shap.Explainer(model_forward, X, feature_names=[f[0] for f in feature_sources])
    shap_values = explainer(X)

    # Check SHAP values shape
    log_message(f"SHAP values shape: {shap_values.values.shape}")

    # Aggregate SHAP values across classes for multi-class models
    if len(shap_values.values.shape) == 3:  # Multi-class case: (n_samples, n_features, n_classes)
        shap_values_mean = np.mean(np.abs(shap_values.values), axis=2)
    elif len(shap_values.values.shape) == 2:  # Single-class case: (n_samples, n_features)
        shap_values_mean = shap_values.values
    else:
        raise ValueError(f"Unexpected SHAP values shape: {shap_values.values.shape}")

    # Aggregate across samples
    shap_values_mean_aggregated = np.mean(shap_values_mean, axis=0)

    # Create DataFrame for variable importance
    shap_importance_df = pd.DataFrame([
        {"source": feature_sources[i][1], "variable": feature_sources[i][0], "importance": shap_values_mean_aggregated[i]}
        for i in range(len(feature_sources))
    ])

    # Aggregate importance by summing up values for the same source and variable
    aggregated_importance = (
        shap_importance_df
        .groupby(["source", "variable"], as_index=False)["importance"]
        .sum()
        .sort_values(by="importance", ascending=False)
    )

    # Save the aggregated data to a CSV file
    aggregated_importance.to_csv(file_name, index=False)
    log_message(f"SHAP importance saved to '{file_name}'.")

    # Beautify SHAP Combined Plot (Bee Swarm + Bar Chart) and save as PDF
    log_message("Generating combined SHAP plot...")

    # Convert SHAP values to numpy array
    shap_values_numpy = shap_values_mean

    # Beautify SHAP Combined Plot (Bee Swarm + Bar Chart) and save as PDF
    with PdfPages("shap_combined_plot_optimized.pdf") as pdf:
        # Create the main figure
        fig, ax1 = plt.subplots(figsize=(10, 8), dpi=300)

        # Draw the bee swarm plot
        shap.summary_plot(
            shap_values_numpy, X, feature_names=[f[0] for f in feature_sources],
            plot_type="dot", show=False, color_bar=True
        )

        # Adjust position for color bar and other elements
        plt.gca().set_position([0.2, 0.2, 0.65, 0.65])   # Leave space for the bar chart and color bar

        # Add second axis on top of the bee swarm plot
        ax2 = ax1.twiny()

        # Draw the bar chart on the second axis
        shap.summary_plot(
            shap_values_numpy, X, feature_names=[f[0] for f in feature_sources],
            plot_type="bar", show=False
        )

        # Align bar chart position with the bee swarm plot
        plt.gca().set_position([0.2, 0.2, 0.65, 0.65]) 

        # Add a horizontal line on the bar chart for emphasis
        ax2.axhline(y=10, color='gray', linestyle='--', linewidth=1)  # Adjust `y` value as necessary

        # Adjust transparency for the bar chart
        bars = ax2.patches
        for bar in bars:
            bar.set_alpha(0.4)  # Set bar transparency to make it visually balanced

        # Customize axis labels and font sizes
        ax1.set_xlabel('Shapley Value Contribution (Bee Swarm)', fontsize=12)
        ax2.set_xlabel('Mean Shapley Value (Feature Importance)', fontsize=12)
        ax1.set_ylabel('Features', fontsize=12)
        ax1.tick_params(labelsize=10)
        ax2.tick_params(labelsize=10)

        # Move top axis ticks and labels to the top
        ax2.xaxis.set_label_position('top')
        ax2.xaxis.tick_top()

        # Adjust layout for better visualization
        plt.tight_layout()

        # Save the plot to the PDF
        pdf.savefig(fig, bbox_inches='tight')

        # Close the plot to free memory
        plt.close()

        log_message("Combined SHAP plot saved to 'shap_combined_plot.pdf'.")

    return aggregated_importance


# Optuna objective function for hyperparameter optimization
def objective(trial, X_train, y_train, X_val, y_val, device):
    try:
        '''
        # Define a broader range of hyperparameters
        num_heads = trial.suggest_categorical("num_heads", [2, 4, 8, 16])  # Include larger values for more complexity
        num_layers = trial.suggest_int("num_layers", 2, 6)  # Allow up to 6 transformer layers
        hidden_size = trial.suggest_categorical("hidden_size", [64, 128, 256, 512])  # Larger hidden dimensions
        dropout_rate = trial.suggest_float("dropout_rate", 0.0, 0.5, step=0.05)  # Wider range for dropout
        lr = trial.suggest_float("lr", 1e-5, 1e-2)  # Log scale for learning rates
        batch_size = trial.suggest_categorical("batch_size", [8, 16, 32, 64])  # Include smaller/larger batch sizes
        '''
       
        # Define a narrower range of hyperparameters for testing
        num_heads = trial.suggest_categorical("num_heads", [2, 4])  # Narrow down to simpler options
        num_layers = trial.suggest_int("num_layers", 2, 3)  # Allow only 2 to 3 transformer layers
        hidden_size = trial.suggest_categorical("hidden_size", [64, 128])  # Reduce hidden dimension choices
        dropout_rate = trial.suggest_float("dropout_rate", 0.1, 0.3, step=0.05)  # Narrow down dropout range
        lr = trial.suggest_float("lr", 1e-4, 1e-3)  # Smaller range for learning rates
        batch_size = trial.suggest_categorical("batch_size", [16, 32])  # Smaller batch sizes for faster testing

        # Define the model with the current set of hyperparameters
        model = TransformerClassifier(
            input_size=X_train.shape[1],
            num_classes=len(set(y_train)),
            num_heads=num_heads,
            num_layers=num_layers,
            hidden_size=hidden_size,
            dropout_rate=dropout_rate,
        ).to(device)

        # Create data loaders
        train_loader = DataLoader(GeneExpressionDataset(X_train, y_train), batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(GeneExpressionDataset(X_val, y_val), batch_size=batch_size, shuffle=False)

        # Train and evaluate the model
        metrics = train_and_evaluate(model, train_loader, val_loader, device, max_epochs=10, lr=lr)
        log_message(f"Trial {trial.number} finished with ROC AUC: {metrics['ROC AUC']} "
                    f"and parameters: {trial.params}")
        return metrics["ROC AUC"]

    except Exception as e:
        log_message(f"Trial {trial.number} failed with error: {e}")
        return 0.0


def perform_optuna_search(X_train, y_train, X_val, y_val, device, n_trials=50):
    """
    Perform Optuna hyperparameter optimization.
    """
    study = optuna.create_study(direction="maximize")
    study.optimize(
        partial(objective, X_train=X_train, y_train=y_train, X_val=X_val, y_val=y_val, device=device),
        n_trials=n_trials,
    )
    return study.best_params, study.best_value


# Evaluate all source combinations
def evaluate_source_combinations(X, y, feature_sources, device, max_workers=10):
    """
    Evaluate all source combinations and save results to CSV files.
    """
    combinations_results = []
    all_importances = []
    source_combinations = generate_source_combinations(feature_sources)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(evaluate_source_combination, combo, X, y, feature_sources, device) for combo in source_combinations]
        for future in tqdm(futures, desc="Evaluating Source Combinations"):
            result = future.result()
            # Append combination metrics only (exclude Variable Importance)
            combinations_results.append({
                "Combination": result["Combination"],
                "Accuracy": result["Accuracy"],
                "Precision": result["Precision"],
                "Recall": result["Recall"],
                "F1 Score": result["F1 Score"],
                "ROC AUC": result["ROC AUC"],
                "Best Params": result["Best Params"]
            })
            # Collect variable importance data if available
            if result.get("Variable Importance") is not None:
                all_importances.append(result["Variable Importance"])

    # Save combination results to CSV
    results_df = pd.DataFrame(combinations_results)
    results_df.to_csv("source_combinations_results.csv", index=False)
    log_message("Source combination results saved to 'source_combinations_results.csv'.")

   # Save and log top 10 variable importance to CSV
    if all_importances:
        variable_importance_df = pd.concat(all_importances, ignore_index=True)
        aggregated_importance = (
            variable_importance_df
            .groupby(["source", "variable"], as_index=False)["importance"]
            .mean()
            .sort_values(by="importance", ascending=False)
        )
        aggregated_importance.to_csv("variable_importance.csv", index=False)
        log_message("Aggregated variable importance saved to 'variable_importance.csv'.")

        # Print and log the top 10 variables
        top_10_variables = aggregated_importance.head(10)
        top_10_text = "\nTop 10 Variables by Importance:\n" + top_10_variables.to_string(index=False)
        log_message(top_10_text)
        print(top_10_text)

    return results_df


def evaluate_source_combination(source_combination, X, y, feature_sources, device):
    """
    Evaluate a single source combination using Optuna for hyperparameter tuning.
    """
    selected_indices = [idx for idx, (_, source) in enumerate(feature_sources) if source in source_combination]
    if not selected_indices:
        return {"Combination": source_combination, "Accuracy": 0.0, "Precision": 0.0, 
                "Recall": 0.0, "F1 Score": 0.0, "ROC AUC": 0.0, "Variable Importance": None}

    X_combination = X[:, selected_indices]
    X_train, X_val, y_train, y_val = preprocess_data(X_combination, y)

    # Perform hyperparameter search
    best_params, best_score = perform_optuna_search(X_train, y_train, X_val, y_val, device)

    # Train and evaluate the model with the best parameters
    model = TransformerClassifier(
        input_size=X_combination.shape[1],
        num_classes=len(set(y)),
        num_heads=best_params["num_heads"],
        num_layers=best_params["num_layers"],
        hidden_size=best_params["hidden_size"],
        dropout_rate=best_params["dropout_rate"],
    ).to(device)

    train_loader = DataLoader(GeneExpressionDataset(X_train, y_train), batch_size=16, shuffle=True)
    val_loader = DataLoader(GeneExpressionDataset(X_val, y_val), batch_size=16, shuffle=False)

    metrics = train_and_evaluate(model, train_loader, val_loader, device)

    # Compute variable importance
    importance_df = compute_variable_importance(model, feature_sources, selected_indices)

    return {
        "Combination": source_combination,
        "Accuracy": metrics["Accuracy"],
        "Precision": metrics["Precision"],
        "Recall": metrics["Recall"],
        "F1 Score": metrics["F1 Score"],
        "ROC AUC": metrics["ROC AUC"],
        "Best Params": best_params,
        "Variable Importance": importance_df
    }

def generate_source_combinations(feature_sources):
    unique_sources = list({source for _, source in feature_sources})
    return [list(combo) for r in range(1, len(unique_sources) + 1) for combo in combinations(unique_sources, r)]

# Main function
def main():
    # Set random seed for reproducibility
    set_random_seed(42)
    
    # Check if CUDA is available and set the device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # Specify the folder path for loading datasets
    folder_path = "/home/xxu/AD/PD/GEOdata/module/"  # Update this path as necessary
    
    # Load and merge datasets
    X, y, feature_sources = load_and_merge_datasets(folder_path)
    
    # Initialize a list to collect all output messages for final printing
    final_outputs = []
    
    # Evaluate all source combinations
    log_message("Starting evaluation of source combinations...")
    results_df = evaluate_source_combinations(X, y, feature_sources, device, max_workers=20)
    
    # Ensure there are results
    if results_df.empty:
        log_message("No valid combinations found. Exiting.")
        final_outputs.append("No valid combinations found. Exiting.")
        print("\n".join(final_outputs))
        return
    
    # Log and print the top 10 combinations by ROC AUC
    top_combinations = results_df.nlargest(10, "ROC AUC")[["Combination", "ROC AUC"]]
    top_combinations_text = "\nTop 10 Combinations by ROC AUC:\n" + top_combinations.to_string(index=False)
    log_message(top_combinations_text)
    final_outputs.append(top_combinations_text)
    
    # Use the best combination (highest ROC AUC) to retrain the model for SHAP computation
    log_message("Retraining model with the best combination for SHAP computation...")
    best_combination = top_combinations.iloc[0]["Combination"]  # Get the combination with the highest ROC AUC
    selected_indices = [
        idx for idx, (_, source) in enumerate(feature_sources) if source in best_combination
    ]
    selected_feature_sources = [feature_sources[idx] for idx in selected_indices]
    
    # Filter data for the best combination
    X_best_combination = X[:, selected_indices]
    
    # Split the data
    X_train, X_val, y_train, y_val = preprocess_data(X_best_combination, y)
    
    # Perform hyperparameter search to find the best parameters
    best_params, _ = perform_optuna_search(X_train, y_train, X_val, y_val, device)
    
    # Retrain the model with the best parameters
    model = TransformerClassifier(
        input_size=X_best_combination.shape[1],
        num_classes=len(set(y)),
        num_heads=best_params["num_heads"],
        num_layers=best_params["num_layers"],
        hidden_size=best_params["hidden_size"],
        dropout_rate=best_params["dropout_rate"]
    ).to(device)
    
    # Create data loaders
    train_loader = DataLoader(GeneExpressionDataset(X_train, y_train), batch_size=best_params["batch_size"], shuffle=True)
    val_loader = DataLoader(GeneExpressionDataset(X_val, y_val), batch_size=best_params["batch_size"], shuffle=False)
    
    # Train the model with the best parameters
    train_and_evaluate(model, train_loader, val_loader, device, max_epochs=10, lr=best_params["lr"])
    
    # Calculate SHAP feature importance
    shap_importance = calculate_shap_importance(model, X_best_combination, selected_feature_sources)
    
    # Print top 10 variables by importance
    top_variables = shap_importance.nlargest(10, "importance")[["variable", "importance"]]
    top_variables_text = "\nTop 10 Variables by Importance:\n" + top_variables.to_string(index=False)
    log_message(top_variables_text)
    final_outputs.append(top_variables_text)
    
    # At the end of the program, print all collected outputs
    print("\n".join(final_outputs))


if __name__ == "__main__":
    main()

