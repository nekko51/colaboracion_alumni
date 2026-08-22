#!/usr/bin/env python3
"""
logreg_meta.py
--------------
Logistic regression meta-classifier that stacks the decision function outputs
of 10 SVM kernels into a single binary classifier for antibody humanness:
  +1 (Human) vs -1 (Murine).

Model:  z(x) = b + sum_i  w_i * hat{f}_i(x)
        P(Human | x) = sigmoid(z(x)) = 1 / (1 + exp(-z(x)))

where hat{f}_i(x) is the z-score normalized i-th SVM decision function output.

Training: L2-regularized binary cross-entropy minimized via gradient descent.
          Training set = validation_vectors_10d.txt  (667 samples)
          Test set     = test_vectors_10d.txt         (668 samples)

Outputs (in results/logreg/):
  logreg_model.txt           - saved model weights and normalization stats
  logreg_val_scores.txt      - per-sample logit, probability, tag (validation)
  logreg_test_scores.txt     - per-sample logit, probability, tag (test)
  logreg_loss_curve.txt      - training loss per epoch
  logreg_weights.txt         - feature weights table (for pgfplots)
  logreg_summary.txt         - human-readable metrics summary
"""

import math
import os
from pathlib import Path


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

KERNEL_NAMES = [
    "Cauchy",
    "Gaussian (RBF)",
    "Linear",
    "Sigmoid",
    "Tanimoto",
    "Piecewise C=0.80 W=0.20",
    "Piecewise C=0.80 W=0.08",
    "Piecewise C=0.82 W=0.08",
    "Piecewise C=0.84 W=0.08",
    "Piecewise C=0.78 W=0.08",
]

DEFAULT_PARAMS = {
    "lambda": 0.01,
    "learning_rate": 0.5,
    "epochs": 2000,
    "threshold": 0.5,
}


# ---------------------------------------------------------------------------
# Pure-Python logistic regression (no external dependencies)
# ---------------------------------------------------------------------------

def sigmoid(z: float) -> float:
    """Numerically stable sigmoid."""
    if z >= 0.0:
        return 1.0 / (1.0 + math.exp(-z))
    else:
        e = math.exp(z)
        return e / (1.0 + e)


def _dot(w, x):
    return sum(wi * xi for wi, xi in zip(w, x))


def logreg_fit(X, y_binary, lam=0.01, lr=0.5, epochs=2000):
    """
    Fits logistic regression via gradient descent on z-score normalized features.

    Parameters
    ----------
    X        : list of lists, shape (n_samples, n_features) - raw (unnormalized) features
    y_binary : list of int, values in {0, 1}
    lam      : L2 regularization coefficient
    lr       : gradient descent learning rate
    epochs   : number of gradient descent steps

    Returns
    -------
    w            : weight vector (in normalized feature space)
    b            : bias scalar
    mu           : per-feature means (for normalization at inference)
    sigma        : per-feature standard deviations
    loss_history : list of per-epoch loss values
    """
    n = len(X)
    d = len(X[0])

    # Compute normalization statistics
    mu = [sum(X[i][j] for i in range(n)) / n for j in range(d)]
    sigma = [
        math.sqrt(sum((X[i][j] - mu[j]) ** 2 for i in range(n)) / n) + 1e-10
        for j in range(d)
    ]

    # Normalize features
    X_norm = [[(X[i][j] - mu[j]) / sigma[j] for j in range(d)] for i in range(n)]

    # Initialize weights
    w = [0.0] * d
    b = 0.0

    loss_history = []

    for _ in range(epochs):
        grad_w = [0.0] * d
        grad_b = 0.0
        loss = 0.0

        for i in range(n):
            z = b + _dot(w, X_norm[i])
            p = sigmoid(z)
            err = p - y_binary[i]

            p_clip = max(1e-12, min(1.0 - 1e-12, p))
            loss -= y_binary[i] * math.log(p_clip) + (1 - y_binary[i]) * math.log(1 - p_clip)

            for j in range(d):
                grad_w[j] += err * X_norm[i][j]
            grad_b += err

        loss = loss / n + 0.5 * lam * sum(wi ** 2 for wi in w)

        for j in range(d):
            grad_w[j] = grad_w[j] / n + lam * w[j]
            w[j] -= lr * grad_w[j]
        grad_b /= n
        b -= lr * grad_b

        loss_history.append(loss)

    return w, b, mu, sigma, loss_history


def logreg_predict(x_raw, w, b, mu, sigma, threshold=0.5):
    """Predict label (+1 or -1) for a single raw feature vector."""
    x_norm = [(x_raw[j] - mu[j]) / sigma[j] for j in range(len(x_raw))]
    z = b + _dot(w, x_norm)
    p = sigmoid(z)
    return 1 if p >= threshold else -1


def logreg_predict_prob(x_raw, w, b, mu, sigma):
    """Returns P(Human = +1) for a single raw feature vector."""
    x_norm = [(x_raw[j] - mu[j]) / sigma[j] for j in range(len(x_raw))]
    z = b + _dot(w, x_norm)
    return sigmoid(z)


def logreg_logit(x_raw, w, b, mu, sigma):
    """Returns the linear logit z(x) for a single raw feature vector."""
    x_norm = [(x_raw[j] - mu[j]) / sigma[j] for j in range(len(x_raw))]
    return b + _dot(w, x_norm)


# ---------------------------------------------------------------------------
# Data I/O
# ---------------------------------------------------------------------------

def load_dataset(filepath, n_features=10):
    """
    Loads a labeled 10D feature dataset.
    Lines starting with '#' are skipped.
    Returns (X, y_pm1, y_binary) as lists.
    """
    X, y_pm1, y_binary = [], [], []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < n_features + 1:
                parts = line.split()
            if len(parts) < n_features + 1:
                continue
            label = int(parts[0])
            # Normalize label: 0 -> -1
            label_pm1 = 1 if label > 0 else -1
            label_bin = 1 if label > 0 else 0
            feats = [float(v) for v in parts[1: n_features + 1]]
            X.append(feats)
            y_pm1.append(label_pm1)
            y_binary.append(label_bin)
    return X, y_pm1, y_binary


def save_scores(filepath, tags_pm1, logits, probs, y_binary, preds_pm1):
    """Saves per-sample scores to a tab-delimited file."""
    with open(filepath, "w") as f:
        f.write("tag\tlogit\tprob\tcorrect\n")
        for tag, logit, prob, yb, pred in zip(tags_pm1, logits, probs, y_binary, preds_pm1):
            correct = 1 if pred == tag else 0
            f.write(f"{tag}\t{logit:.15e}\t{prob:.15e}\t{correct}\n")


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def confusion_matrix(y_true_pm1, y_pred_pm1, target_class=1):
    tp = sum(1 for t, p in zip(y_true_pm1, y_pred_pm1) if t == target_class and p == target_class)
    tn = sum(1 for t, p in zip(y_true_pm1, y_pred_pm1) if t != target_class and p != target_class)
    fp = sum(1 for t, p in zip(y_true_pm1, y_pred_pm1) if t != target_class and p == target_class)
    fn = sum(1 for t, p in zip(y_true_pm1, y_pred_pm1) if t == target_class and p != target_class)
    return tp, tn, fp, fn


def f1(tp, fp, fn):
    denom = 2 * tp + fp + fn
    return 2 * tp / denom if denom > 0 else 0.0


def print_metrics(name, y_true_pm1, y_pred_pm1):
    """Prints a human-readable classification report."""
    tp_h, tn_h, fp_h, fn_h = confusion_matrix(y_true_pm1, y_pred_pm1, target_class=1)
    tp_m, tn_m, fp_m, fn_m = confusion_matrix(y_true_pm1, y_pred_pm1, target_class=-1)
    total = len(y_true_pm1)
    correct = tp_h + tn_h
    acc = correct / total if total > 0 else 0.0
    f1_h = f1(tp_h, fp_h, fn_h)
    f1_m = f1(tp_m, fp_m, fn_m)

    print("=" * 70)
    print(f" Logistic Regression Evaluation: {name}")
    print("=" * 70)
    print(f" Total: {total}  (Human [+1]: {tp_h + fn_h}, Murine [-1]: {tp_m + fn_m})")
    print(f" Accuracy:  {acc:.6f}  ({acc * 100:.2f}%)")
    print(f" F1 Human:  {f1_h:.6f}")
    print(f" F1 Murine: {f1_m:.6f}")
    print("-" * 70)
    print(" Confusion Matrix (Human = Positive [+1]):")
    print("                  Pred Human   Pred Murine")
    print(f"   True Human:    {tp_h:6d} (TP)    {fn_h:6d} (FN)")
    print(f"   True Murine:   {fp_h:6d} (FP)    {tn_h:6d} (TN)")
    print("=" * 70)
    print()
    return acc, f1_h, f1_m, tp_h, tn_h, fp_h, fn_h


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run(
    val_path="results/cart/validation_vectors_10d.txt",
    test_path="results/cart/test_vectors_10d.txt",
    output_dir="results/logreg",
    params=None,
    kernel_names=None,
):
    """
    Full logistic regression training and evaluation pipeline.

    Parameters
    ----------
    val_path    : path to 10D validation/training vectors
    test_path   : path to 10D test vectors
    output_dir  : directory for output files
    params      : dict with keys lambda, learning_rate, epochs, threshold
    kernel_names: list of 10 feature/kernel names
    """
    if params is None:
        params = DEFAULT_PARAMS
    if kernel_names is None:
        kernel_names = KERNEL_NAMES

    os.makedirs(output_dir, exist_ok=True)

    print("=" * 70)
    print(" LOGISTIC REGRESSION META-CLASSIFIER (10D SVM Kernel Stacking)")
    print("=" * 70)
    print(f"  lambda={params['lambda']}, lr={params['learning_rate']}, epochs={params['epochs']}")
    print()

    # --- Load data ---
    print(f"Loading training set: {val_path}")
    X_val, y_val_pm1, y_val_bin = load_dataset(val_path)
    print(f"  {len(X_val)} samples loaded.")

    print(f"Loading test set:     {test_path}")
    X_test, y_test_pm1, y_test_bin = load_dataset(test_path)
    print(f"  {len(X_test)} samples loaded.\n")

    # --- Train ---
    print("Fitting logistic regression...")
    w, b, mu, sigma, loss_hist = logreg_fit(
        X_val, y_val_bin,
        lam=params["lambda"],
        lr=params["learning_rate"],
        epochs=params["epochs"],
    )
    print("Training complete.\n")

    # --- Predict ---
    threshold = params["threshold"]

    def _score_dataset(X, y_pm1):
        logits = [logreg_logit(x, w, b, mu, sigma) for x in X]
        probs  = [sigmoid(z) for z in logits]
        preds  = [1 if p >= threshold else -1 for p in probs]
        return logits, probs, preds

    logits_val, probs_val, preds_val = _score_dataset(X_val, y_val_pm1)
    logits_test, probs_test, preds_test = _score_dataset(X_test, y_test_pm1)

    # --- Evaluate ---
    val_acc, val_f1h, val_f1m, val_tp, val_tn, val_fp, val_fn = print_metrics(
        "Validation Set (Training)", y_val_pm1, preds_val
    )
    test_acc, test_f1h, test_f1m, test_tp, test_tn, test_fp, test_fn = print_metrics(
        "Test Set (Generalization)", y_test_pm1, preds_test
    )

    # --- Save loss curve ---
    loss_path = os.path.join(output_dir, "logreg_loss_curve.txt")
    with open(loss_path, "w") as f:
        f.write("epoch\tloss\n")
        for ep, loss in enumerate(loss_hist):
            f.write(f"{ep + 1}\t{loss:.15e}\n")
    print(f"Saved: {loss_path}")

    # --- Save scores ---
    val_scores_path = os.path.join(output_dir, "logreg_val_scores.txt")
    save_scores(val_scores_path, y_val_pm1, logits_val, probs_val, y_val_bin, preds_val)
    print(f"Saved: {val_scores_path}")

    test_scores_path = os.path.join(output_dir, "logreg_test_scores.txt")
    save_scores(test_scores_path, y_test_pm1, logits_test, probs_test, y_test_bin, preds_test)
    print(f"Saved: {test_scores_path}")

    # --- Save weights table ---
    w_abs = [abs(wi) for wi in w]
    w_total = sum(w_abs)
    weights_path = os.path.join(output_dir, "logreg_weights.txt")
    with open(weights_path, "w") as f:
        f.write("feature\tweight\tweight_abs\tweight_norm\n")
        for j, name in enumerate(kernel_names):
            f.write(f"{name}\t{w[j]:.10f}\t{w_abs[j]:.10f}\t{w_abs[j] / w_total:.10f}\n")
    print(f"Saved: {weights_path}")

    # --- Save normalization stats ---
    norm_path = os.path.join(output_dir, "logreg_normalization.txt")
    with open(norm_path, "w") as f:
        f.write("feature\tmean\tstd\n")
        for j, name in enumerate(kernel_names):
            f.write(f"{name}\t{mu[j]:.10f}\t{sigma[j]:.10f}\n")
    print(f"Saved: {norm_path}")

    # --- Save model ---
    model_path = os.path.join(output_dir, "logreg_model.txt")
    with open(model_path, "w") as f:
        f.write("LOGREG_MODEL_V1\n")
        f.write(f"n_features {len(w)}\n")
        f.write(f"params {params['lambda']:.15e} {params['learning_rate']:.15e} "
                f"{params['epochs']} {params['threshold']:.15e}\n")
        f.write(f"bias {b:.15e}\n")
        f.write("weights " + " ".join(f"{wi:.15e}" for wi in w) + "\n")
        f.write("feature_mean " + " ".join(f"{m:.15e}" for m in mu) + "\n")
        f.write("feature_std " + " ".join(f"{s:.15e}" for s in sigma) + "\n")
    print(f"Saved: {model_path}")

    # --- Save human-readable summary ---
    summary_path = os.path.join(output_dir, "logreg_summary.txt")
    with open(summary_path, "w") as f:
        f.write("LOGISTIC REGRESSION META-CLASSIFIER SUMMARY\n")
        f.write("=" * 60 + "\n")
        f.write(f"Hyperparameters: lambda={params['lambda']}, lr={params['learning_rate']}, "
                f"epochs={params['epochs']}, threshold={params['threshold']}\n\n")
        f.write(f"Training (Validation) set: {len(X_val)} samples\n")
        f.write(f"  Accuracy:  {val_acc:.6f}\n")
        f.write(f"  F1 Human:  {val_f1h:.6f}\n")
        f.write(f"  F1 Murine: {val_f1m:.6f}\n")
        f.write(f"  TP={val_tp}, TN={val_tn}, FP={val_fp}, FN={val_fn}\n\n")
        f.write(f"Test set: {len(X_test)} samples\n")
        f.write(f"  Accuracy:  {test_acc:.6f}\n")
        f.write(f"  F1 Human:  {test_f1h:.6f}\n")
        f.write(f"  F1 Murine: {test_f1m:.6f}\n")
        f.write(f"  TP={test_tp}, TN={test_tn}, FP={test_fp}, FN={test_fn}\n\n")
        f.write("Weights (in normalized space):\n")
        for j, name in enumerate(kernel_names):
            f.write(f"  w[{j:2d}] {name:<25s} = {w[j]:+.6f}  ({w_abs[j]/w_total*100:.2f}% of |w|)\n")
        f.write(f"  bias                           = {b:+.6f}\n")
    print(f"Saved: {summary_path}")
    print()

    return {
        "w": w, "b": b, "mu": mu, "sigma": sigma,
        "val_acc": val_acc, "val_f1h": val_f1h, "val_f1m": val_f1m,
        "test_acc": test_acc, "test_f1h": test_f1h, "test_f1m": test_f1m,
        "loss_history": loss_hist,
    }


if __name__ == "__main__":
    run()
