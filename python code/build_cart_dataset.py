#!/usr/bin/env python3
"""
build_cart_dataset.py
--------------------
Extracts 10-dimensional SVM decision function vectors and class labels (+1 for Human, -1 for Murine)
from the SVM test results directory for both the validation set (667 samples) and test set (668 samples).
Outputs the formatted dataset files ready for the CART Decision Tree algorithm in c-code/cart.c.
"""

import os
from pathlib import Path


DEFAULT_10_KERNELS = [
    'kernel-c',                 # 0: Cauchy
    'kernel-g',                 # 1: Gaussian / RBF
    'kernel-l',                 # 2: Linear
    'kernel-s',                 # 3: Sigmoid
    'kernel-t',                 # 4: Tanimoto
    'kernel-w_C-0.80_W-0.20',   # 5: Piecewise C=0.80 W=0.20
    'kernel-w_C-0.80_W-0.08',   # 6: Piecewise C=0.80 W=0.08
    'kernel-w_C-0.82_W-0.08',   # 7: Piecewise C=0.82 W=0.08
    'kernel-w_C-0.84_W-0.08',   # 8: Piecewise C=0.84 W=0.08
    'kernel-w_C-0.78_W-0.08',   # 9: Piecewise C=0.78 W=0.08
]


def find_svm_runs(tests_dir, split_len, kernel_list=None):
    """
    Finds one valid SVM run file for each specified kernel matching the requested sample count (split_len).
    """
    if kernel_list is None:
        kernel_list = DEFAULT_10_KERNELS

    runs_by_kernel = {}

    for k_name in kernel_list:
        k_path = os.path.join(tests_dir, k_name)
        if not os.path.isdir(k_path):
            print(f"Warning: Kernel directory not found: {k_path}")
            continue

        # Look through date/time subdirectories for the LAST matching test-results.txt
        subdirs = sorted(os.listdir(k_path))
        matched_file = None
        matched_subdir = None
        for subdir in subdirs:
            candidate = os.path.join(k_path, subdir, 'test-results.txt')
            if os.path.isfile(candidate):
                with open(candidate, 'r') as f:
                    lines = [line.strip() for line in f if line.strip()]
                if len(lines) == split_len:
                    matched_file = candidate
                    matched_subdir = subdir

        if matched_file:
            runs_by_kernel[k_name] = matched_file
            print(f"  [{k_name}] Selected latest run: {matched_subdir} ({split_len} samples)")
        else:
            print(f"Warning: No run with {split_len} samples found for {k_name}")

    return runs_by_kernel


def extract_vectors(runs_by_kernel, split_len, kernel_list=None):
    """
    Extracts 10-dimensional decision function vectors and tags from the list of matched runs.
    Returns:
        tags: list of int (+1 or -1)
        vectors: list of lists of float (10 dimensions per sample)
        kernel_names: list of kernel names corresponding to each feature dimension
    """
    if kernel_list is None:
        kernel_list = DEFAULT_10_KERNELS

    valid_kernels = [k for k in kernel_list if k in runs_by_kernel]
    if len(valid_kernels) != len(kernel_list):
        print(f"Notice: Found {len(valid_kernels)} out of {len(kernel_list)} requested kernels.")

    tags = None
    feature_columns = []

    for k_name in valid_kernels:
        file_path = runs_by_kernel[k_name]
        with open(file_path, 'r') as f:
            lines = [line.strip().split('\t') for line in f if line.strip()]

        current_tags = [int(parts[0]) for parts in lines]
        # Column 1 is the raw decision function output f(x)
        current_features = [float(parts[1]) if len(parts) > 1 else float(parts[0]) for parts in lines]

        if tags is None:
            tags = current_tags
        elif tags != current_tags:
            print(f"Warning: Tag mismatch in {file_path}")

        feature_columns.append(current_features)

    # Transpose feature_columns from (n_features, n_samples) to (n_samples, n_features)
    n_samples = len(tags)
    n_features = len(feature_columns)
    vectors = []
    for i in range(n_samples):
        sample_vec = [feature_columns[j][i] for j in range(n_features)]
        vectors.append(sample_vec)

    return tags, vectors, valid_kernels


def save_dataset_file(output_path, tags, vectors, kernel_names):
    """
    Saves the extracted dataset in tab-delimited format:
    # tag <feat_0> <feat_1> ... <feat_9>
    <tag>\t<f0>\t<f1>\t...\t<f9>
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        # Header comment with feature names
        header = "# tag\t" + "\t".join(kernel_names) + "\n"
        f.write(header)
        for tag, vec in zip(tags, vectors):
            vec_str = "\t".join(f"{val:.15e}" for val in vec)
            f.write(f"{tag}\t{vec_str}\n")
    print(f"Saved {len(tags)} samples to: {output_path}")


def build_all_cart_datasets(base_dir=None, kernel_list=None):
    """
    Builds both validation (667 samples) and test (668 samples) 10D vector datasets.
    """
    if base_dir is None:
        base_dir = Path(__file__).resolve().parent.parent
    else:
        base_dir = Path(base_dir)

    if kernel_list is None:
        kernel_list = DEFAULT_10_KERNELS

    tests_dir = os.path.join(base_dir, 'results', 'svm', 'decision-function-tests')
    output_dir = os.path.join(base_dir, 'results', 'cart')

    print("=" * 70)
    print("Extracting 10-Dimensional SVM Vectors for CART Decision Tree")
    print("=" * 70)
    print(f"Source tests directory: {tests_dir}")
    print(f"Output directory:       {output_dir}")
    print(f"Selected kernels ({len(kernel_list)}):")
    for idx, k in enumerate(kernel_list):
        print(f"  Feature {idx:2d}: {k}")
    print("-" * 70)

    # 1. Validation Set (667 samples) - Used as Training set for CART
    val_runs = find_svm_runs(tests_dir, split_len=667, kernel_list=kernel_list)
    val_tags, val_vectors, val_kernels = extract_vectors(val_runs, split_len=667, kernel_list=kernel_list)
    val_output_path = os.path.join(output_dir, 'validation_vectors_10d.txt')
    save_dataset_file(val_output_path, val_tags, val_vectors, val_kernels)

    # 2. Test Set (668 samples) - Used for Final Testing of CART
    test_runs = find_svm_runs(tests_dir, split_len=668, kernel_list=kernel_list)
    test_tags, test_vectors, test_kernels = extract_vectors(test_runs, split_len=668, kernel_list=kernel_list)
    test_output_path = os.path.join(output_dir, 'test_vectors_10d.txt')
    save_dataset_file(test_output_path, test_tags, test_vectors, test_kernels)

    print("-" * 70)
    print("Summary:")
    print(f"  Validation (Train) samples: {len(val_tags)} ({val_tags.count(1)} Human, {val_tags.count(-1)} Murine)")
    print(f"  Test samples:               {len(test_tags)} ({test_tags.count(1)} Human, {test_tags.count(-1)} Murine)")
    print("=" * 70)
    return val_output_path, test_output_path


if __name__ == '__main__':
    build_all_cart_datasets()
