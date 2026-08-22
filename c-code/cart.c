#include "head.h"

/**
 * ============================================================================
 * CART (Classification and Regression Trees) Decision Tree Algorithm
 * ============================================================================
 * 
 * Implements a binary classification decision tree for antibody humanness
 * classification based on 10-dimensional SVM decision score vectors.
 * 
 * Features:
 *   - 10-dimensional feature vector input (representing outputs of 10 SVM models)
 *   - Binary target classification: +1 (Human) vs -1 (Murine)
 *   - Gini impurity based optimal split selection
 *   - Hyperparameters: max_depth, min_samples_split, min_samples_leaf, min_impurity_decrease
 *   - Model training, evaluation (ConfusionMatrix, F1-scores, Accuracy), and tree visualization
 *   - Dataset loading/saving and tree serialization/deserialization
 * ============================================================================
 */

#define CART_DEFAULT_DIM 10
#define CART_MAX_FEATURES 64

/* Struct for a single sample with a feature vector and class label */
typedef struct {
    double features[CART_MAX_FEATURES];
    int n_features;
    int label; // +1 (Human) or -1 (Murine)
} CartSample;

/* Struct for a dataset containing multiple samples */
typedef struct {
    CartSample *samples;
    int n_samples;
    int n_features;
} CartDataset;

/* Struct for CART hyperparameters */
typedef struct {
    int max_depth;                 /* Maximum tree depth (-1 for unlimited) */
    int min_samples_split;         /* Minimum samples required to split an internal node */
    int min_samples_leaf;          /* Minimum samples required in a leaf node */
    double min_impurity_decrease;  /* Minimum impurity decrease required for a split */
} CartParams;

/* Struct for a node in the Decision Tree */
typedef struct CartNode {
    bool is_leaf;
    int feature_index;             /* Index of feature used for splitting (0 .. n_features-1) */
    double threshold;              /* Split threshold: x[feature] <= threshold -> left, else right */
    int prediction;                /* Class prediction: +1 or -1 */
    double probability_pos;        /* Probability of class +1 at this node */
    double impurity;               /* Gini impurity at this node */
    int n_samples;                 /* Total samples at this node */
    int n_pos;                     /* Number of +1 samples */
    int n_neg;                     /* Number of -1 samples */
    struct CartNode *left;         /* Left child (x[feature] <= threshold) */
    struct CartNode *right;        /* Right child (x[feature] > threshold) */
} CartNode;

/* Struct for a trained Decision Tree model */
typedef struct {
    CartNode *root;
    CartParams params;
    int n_features;
    int total_nodes;
    int max_depth_reached;
} CartTree;

/* Struct used internally for sorting feature values */
typedef struct {
    double value;
    int label;
    int original_index;
} CartSortItem;

/* Forward declarations */
CartParams cart_default_params(void);
double cart_gini(int n_pos, int n_neg);
CartNode* cart_create_leaf(int n_pos, int n_neg);
void cart_free_node(CartNode *node);
void cart_free_tree(CartTree *tree);
CartTree* cart_fit(const CartDataset *dataset, const CartParams *params);
int cart_predict(const CartTree *tree, const double *vector, int n_features);
int cart_predict_10d(const CartTree *tree, const double vector10d[10]);
double cart_predict_probability(const CartTree *tree, const double *vector, int n_features);
ConfusionMatrix cart_evaluate(const CartTree *tree, const CartDataset *dataset, int target_class);
void cart_print_metrics(const CartTree *tree, const CartDataset *dataset, const char *dataset_name);
void cart_print_tree_recursive(const CartNode *node, int depth, const char **feature_names, FILE *out);
void cart_print_tree(const CartTree *tree, const char **feature_names);
int cart_load_dataset(const char *filepath, CartDataset *out_dataset, int n_features);
void cart_free_dataset(CartDataset *dataset);
int cart_save_tree(const CartTree *tree, const char *filepath);
CartTree* cart_load_tree(const char *filepath);
void cart_train_and_test(const char *validation_file, const char *test_file, const CartParams *params);


/**
 * @brief Returns default hyperparameters for CART decision tree.
 */
CartParams cart_default_params(void) {
    CartParams params;
    params.max_depth = 5;
    params.min_samples_split = 2;
    params.min_samples_leaf = 1;
    params.min_impurity_decrease = 1e-7;
    return params;
}

/**
 * @brief Computes Gini impurity for a binary classification node.
 * Gini = 1 - p(+1)^2 - p(-1)^2 = 2 * p(+1) * p(-1)
 */
double cart_gini(int n_pos, int n_neg) {
    int total = n_pos + n_neg;
    if (total <= 0) return 0.0;
    double p_pos = (double)n_pos / (double)total;
    double p_neg = (double)n_neg / (double)total;
    return 1.0 - (p_pos * p_pos + p_neg * p_neg);
}

/**
 * @brief Creates a leaf node with class prediction and statistics.
 */
CartNode* cart_create_leaf(int n_pos, int n_neg) {
    CartNode *node = (CartNode*)malloc(sizeof(CartNode));
    if (!node) {
        fprintf(stderr, "Error: Memory allocation failed for CartNode.\n");
        return NULL;
    }
    node->is_leaf = true;
    node->feature_index = -1;
    node->threshold = 0.0;
    node->n_pos = n_pos;
    node->n_neg = n_neg;
    node->n_samples = n_pos + n_neg;
    node->impurity = cart_gini(n_pos, n_neg);
    node->prediction = (n_pos >= n_neg) ? 1 : -1;
    node->probability_pos = (node->n_samples > 0) ? ((double)n_pos / (double)node->n_samples) : 0.5;
    node->left = NULL;
    node->right = NULL;
    return node;
}

/**
 * @brief Comparator for qsort to sort CartSortItem ascending by value.
 */
static int compare_sort_items(const void *a, const void *b) {
    const CartSortItem *item_a = (const CartSortItem*)a;
    const CartSortItem *item_b = (const CartSortItem*)b;
    if (item_a->value < item_b->value) return -1;
    if (item_a->value > item_b->value) return 1;
    return 0;
}

/**
 * @brief Recursive tree builder using the CART algorithm.
 */
static CartNode* cart_build_recursive(const CartDataset *dataset, const int *indices, int n_samples,
                                     int depth, const CartParams *params, int *total_nodes, int *max_depth) {
    if (n_samples <= 0) return NULL;

    if (depth > *max_depth) {
        *max_depth = depth;
    }
    (*total_nodes)++;

    // Count class distribution
    int n_pos = 0;
    int n_neg = 0;
    for (int i = 0; i < n_samples; i++) {
        int idx = indices[i];
        if (dataset->samples[idx].label > 0) {
            n_pos++;
        } else {
            n_neg++;
        }
    }

    double current_impurity = cart_gini(n_pos, n_neg);

    // Check stopping conditions
    bool should_stop = false;
    if (n_pos == 0 || n_neg == 0) {
        should_stop = true; // Pure node
    } else if (params->max_depth >= 0 && depth >= params->max_depth) {
        should_stop = true; // Reached maximum depth
    } else if (n_samples < params->min_samples_split) {
        should_stop = true; // Not enough samples to split
    }

    if (should_stop) {
        return cart_create_leaf(n_pos, n_neg);
    }

    // Search for best feature and threshold split
    int best_feature = -1;
    double best_threshold = 0.0;
    double best_impurity_decrease = -1.0;

    CartSortItem *sort_items = (CartSortItem*)malloc(n_samples * sizeof(CartSortItem));
    if (!sort_items) {
        fprintf(stderr, "Error: Failed to allocate memory for sorting split items.\n");
        return cart_create_leaf(n_pos, n_neg);
    }

    for (int f = 0; f < dataset->n_features; f++) {
        for (int i = 0; i < n_samples; i++) {
            int idx = indices[i];
            sort_items[i].value = dataset->samples[idx].features[f];
            sort_items[i].label = dataset->samples[idx].label;
            sort_items[i].original_index = idx;
        }

        qsort(sort_items, n_samples, sizeof(CartSortItem), compare_sort_items);

        int left_pos = 0;
        int left_neg = 0;
        int right_pos = n_pos;
        int right_neg = n_neg;

        for (int i = 0; i < n_samples - 1; i++) {
            if (sort_items[i].label > 0) {
                left_pos++;
                right_pos--;
            } else {
                left_neg++;
                right_neg--;
            }

            // Only split between distinct feature values
            if (sort_items[i].value < sort_items[i + 1].value) {
                int left_total = left_pos + left_neg;
                int right_total = right_pos + right_neg;

                if (left_total < params->min_samples_leaf || right_total < params->min_samples_leaf) {
                    continue;
                }

                double left_impurity = cart_gini(left_pos, left_neg);
                double right_impurity = cart_gini(right_pos, right_neg);
                double split_impurity = ((double)left_total / (double)n_samples) * left_impurity +
                                       ((double)right_total / (double)n_samples) * right_impurity;
                double impurity_decrease = current_impurity - split_impurity;

                if (impurity_decrease > best_impurity_decrease) {
                    best_impurity_decrease = impurity_decrease;
                    best_feature = f;
                    best_threshold = (sort_items[i].value + sort_items[i + 1].value) / 2.0;
                }
            }
        }
    }

    free(sort_items);

    // If no valid split found or improvement is below threshold, make a leaf
    if (best_feature == -1 || best_impurity_decrease < params->min_impurity_decrease) {
        return cart_create_leaf(n_pos, n_neg);
    }

    // Partition indices into left and right subsets based on best split
    int *left_indices = (int*)malloc(n_samples * sizeof(int));
    int *right_indices = (int*)malloc(n_samples * sizeof(int));
    if (!left_indices || !right_indices) {
        if (left_indices) free(left_indices);
        if (right_indices) free(right_indices);
        return cart_create_leaf(n_pos, n_neg);
    }

    int left_count = 0;
    int right_count = 0;
    for (int i = 0; i < n_samples; i++) {
        int idx = indices[i];
        if (dataset->samples[idx].features[best_feature] <= best_threshold) {
            left_indices[left_count++] = idx;
        } else {
            right_indices[right_count++] = idx;
        }
    }

    // Allocate internal decision node
    CartNode *node = (CartNode*)malloc(sizeof(CartNode));
    if (!node) {
        free(left_indices);
        free(right_indices);
        return cart_create_leaf(n_pos, n_neg);
    }

    node->is_leaf = false;
    node->feature_index = best_feature;
    node->threshold = best_threshold;
    node->n_samples = n_samples;
    node->n_pos = n_pos;
    node->n_neg = n_neg;
    node->impurity = current_impurity;
    node->prediction = (n_pos >= n_neg) ? 1 : -1;
    node->probability_pos = (double)n_pos / (double)n_samples;

    // Recursively build children
    node->left = cart_build_recursive(dataset, left_indices, left_count, depth + 1, params, total_nodes, max_depth);
    node->right = cart_build_recursive(dataset, right_indices, right_count, depth + 1, params, total_nodes, max_depth);

    free(left_indices);
    free(right_indices);

    return node;
}

/**
 * @brief Fits a CART Decision Tree to the given dataset.
 * 
 * @param dataset The training dataset containing feature vectors and labels.
 * @param params Hyperparameters for the tree (or NULL for default parameters).
 * @return CartTree* Pointer to the trained tree, or NULL on error.
 */
CartTree* cart_fit(const CartDataset *dataset, const CartParams *params) {
    if (!dataset || dataset->n_samples <= 0 || dataset->n_features <= 0) {
        fprintf(stderr, "Error: Invalid dataset provided to cart_fit.\n");
        return NULL;
    }

    CartTree *tree = (CartTree*)malloc(sizeof(CartTree));
    if (!tree) {
        fprintf(stderr, "Error: Failed to allocate CartTree struct.\n");
        return NULL;
    }

    if (params) {
        tree->params = *params;
    } else {
        tree->params = cart_default_params();
    }

    tree->n_features = dataset->n_features;
    tree->total_nodes = 0;
    tree->max_depth_reached = 0;

    int *all_indices = (int*)malloc(dataset->n_samples * sizeof(int));
    if (!all_indices) {
        free(tree);
        return NULL;
    }
    for (int i = 0; i < dataset->n_samples; i++) {
        all_indices[i] = i;
    }

    tree->root = cart_build_recursive(dataset, all_indices, dataset->n_samples, 0,
                                      &tree->params, &tree->total_nodes, &tree->max_depth_reached);

    free(all_indices);
    return tree;
}

/**
 * @brief Predicts class label (+1 or -1) for an input feature vector.
 */
int cart_predict(const CartTree *tree, const double *vector, int n_features) {
    if (!tree || !tree->root || !vector) return 0;
    if (n_features < tree->n_features) {
        fprintf(stderr, "Warning: Vector features (%d) less than tree features (%d).\n",
                n_features, tree->n_features);
    }

    const CartNode *curr = tree->root;
    while (!curr->is_leaf) {
        if (vector[curr->feature_index] <= curr->threshold) {
            curr = curr->left;
        } else {
            curr = curr->right;
        }
        if (!curr) return 0;
    }
    return curr->prediction;
}

/**
 * @brief Predicts class label (+1 for Human, -1 for Murine) for a 10-dimensional vector.
 */
int cart_predict_10d(const CartTree *tree, const double vector10d[10]) {
    return cart_predict(tree, vector10d, 10);
}

/**
 * @brief Convenience classifier taking 10 explicit double values.
 */
int cart_classify_10d(const CartTree *tree, double v0, double v1, double v2, double v3, double v4,
                      double v5, double v6, double v7, double v8, double v9) {
    double vec[10] = {v0, v1, v2, v3, v4, v5, v6, v7, v8, v9};
    return cart_predict_10d(tree, vec);
}

/**
 * @brief Predicts the class probability P(class = +1) for an input feature vector.
 */
double cart_predict_probability(const CartTree *tree, const double *vector, int n_features) {
    if (!tree || !tree->root || !vector) return 0.5;

    const CartNode *curr = tree->root;
    while (!curr->is_leaf) {
        if (vector[curr->feature_index] <= curr->threshold) {
            curr = curr->left;
        } else {
            curr = curr->right;
        }
        if (!curr) return 0.5;
    }
    return curr->probability_pos;
}

/**
 * @brief Predicts labels for all samples in a dataset.
 */
void cart_predict_dataset(const CartTree *tree, const CartDataset *dataset, int *out_predictions) {
    if (!tree || !dataset || !out_predictions) return;
    for (int i = 0; i < dataset->n_samples; i++) {
        out_predictions[i] = cart_predict(tree, dataset->samples[i].features, dataset->n_features);
    }
}

/**
 * @brief Evaluates tree predictions on a dataset and computes a ConfusionMatrix.
 */
ConfusionMatrix cart_evaluate(const CartTree *tree, const CartDataset *dataset, int target_class) {
    ConfusionMatrix cm = {0, 0, 0, 0};
    if (!tree || !dataset) return cm;

    for (int i = 0; i < dataset->n_samples; i++) {
        int true_tag = dataset->samples[i].label;
        int pred_tag = cart_predict(tree, dataset->samples[i].features, dataset->n_features);

        if (true_tag == target_class) {
            if (pred_tag == target_class) cm.TP++;
            else cm.FN++;
        } else {
            if (pred_tag == target_class) cm.FP++;
            else cm.TN++;
        }
    }
    return cm;
}

/**
 * @brief Prints detailed classification metrics (Accuracy, F1-Human, F1-Murine, Confusion Matrix).
 */
void cart_print_metrics(const CartTree *tree, const CartDataset *dataset, const char *dataset_name) {
    if (!tree || !dataset) return;

    ConfusionMatrix cm_h = cart_evaluate(tree, dataset, 1);
    ConfusionMatrix cm_m = cart_evaluate(tree, dataset, -1);

    int total = dataset->n_samples;
    int correct = cm_h.TP + cm_h.TN;
    double accuracy = (total > 0) ? ((double)correct / (double)total) : 0.0;

    double f1_human = (2.0 * cm_h.TP + cm_h.FP + cm_h.FN > 0) ?
                      (2.0 * cm_h.TP / (2.0 * cm_h.TP + cm_h.FP + cm_h.FN)) : 0.0;

    double f1_murine = (2.0 * cm_m.TP + cm_m.FP + cm_m.FN > 0) ?
                       (2.0 * cm_m.TP / (2.0 * cm_m.TP + cm_m.FP + cm_m.FN)) : 0.0;

    printf("======================================================================\n");
    printf(" CART Decision Tree Evaluation: %s\n", dataset_name ? dataset_name : "Dataset");
    printf("======================================================================\n");
    printf(" Total Samples: %d (Human [+1]: %d, Murine [-1]: %d)\n",
           total, cm_h.TP + cm_h.FN, cm_h.FP + cm_h.TN);
    printf(" Correctly Classified: %d (%.2f%%)\n", correct, accuracy * 100.0);
    printf(" Accuracy:  %.6f (%.2f%%)\n", accuracy, accuracy * 100.0);
    printf(" F1 Human:  %.6f\n", f1_human);
    printf(" F1 Murine: %.6f\n", f1_murine);
    printf("----------------------------------------------------------------------\n");
    printf(" Confusion Matrix (Human = Positive [+1]):\n");
    printf("                  Predicted Human (+1)   Predicted Murine (-1)\n");
    printf("   Actual Human:       %6d (TP)                %6d (FN)\n", cm_h.TP, cm_h.FN);
    printf("   Actual Murine:      %6d (FP)                %6d (TN)\n", cm_h.FP, cm_h.TN);
    printf("======================================================================\n\n");
}

/**
 * @brief Helper for recursive tree printing.
 */
void cart_print_tree_recursive(const CartNode *node, int depth, const char **feature_names, FILE *out) {
    if (!node || !out) return;

    for (int i = 0; i < depth; i++) {
        fprintf(out, "|   ");
    }

    if (node->is_leaf) {
        const char *class_str = (node->prediction > 0) ? "Human (+1)" : "Murine (-1)";
        fprintf(out, "|--- Leaf: %s [P(+1)=%.3f, N=%d, pos=%d, neg=%d, gini=%.4f]\n",
                class_str, node->probability_pos, node->n_samples, node->n_pos, node->n_neg, node->impurity);
    } else {
        char feat_str[64];
        if (feature_names && feature_names[node->feature_index]) {
            snprintf(feat_str, sizeof(feat_str), "%s (f%d)", feature_names[node->feature_index], node->feature_index);
        } else {
            snprintf(feat_str, sizeof(feat_str), "feature_%d", node->feature_index);
        }

        fprintf(out, "|--- [%s] <= %.6lf [N=%d, pos=%d, neg=%d, gini=%.4f]\n",
                feat_str, node->threshold, node->n_samples, node->n_pos, node->n_neg, node->impurity);

        cart_print_tree_recursive(node->left, depth + 1, feature_names, out);

        for (int i = 0; i < depth; i++) {
            fprintf(out, "|   ");
        }
        fprintf(out, "|--- [%s] > %.6lf\n", feat_str, node->threshold);
        cart_print_tree_recursive(node->right, depth + 1, feature_names, out);
    }
}

/**
 * @brief Prints ASCII representation of the decision tree to stdout.
 */
void cart_print_tree(const CartTree *tree, const char **feature_names) {
    if (!tree || !tree->root) {
        printf("Empty or uninitialized CART tree.\n");
        return;
    }

    printf("======================================================================\n");
    printf(" CART Decision Tree Structure (Total Nodes: %d, Max Depth: %d)\n",
           tree->total_nodes, tree->max_depth_reached);
    printf("======================================================================\n");
    cart_print_tree_recursive(tree->root, 0, feature_names, stdout);
    printf("======================================================================\n\n");
}

/**
 * @brief Loads a dataset from a text file.
 * Format of each line: <tag>\t<f0>\t<f1>\t...\t<f9>
 * Lines starting with '#' are treated as comments and skipped.
 */
int cart_load_dataset(const char *filepath, CartDataset *out_dataset, int n_features) {
    if (!filepath || !out_dataset) return 0;
    if (n_features <= 0 || n_features > CART_MAX_FEATURES) {
        n_features = CART_DEFAULT_DIM;
    }

    FILE *f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "Error: Could not open dataset file: %s\n", filepath);
        return 0;
    }

    int capacity = 512;
    int count = 0;
    CartSample *samples = (CartSample*)malloc(capacity * sizeof(CartSample));
    if (!samples) {
        fclose(f);
        return 0;
    }

    char line_buf[4096];
    while (fgets(line_buf, sizeof(line_buf), f)) {
        char *ptr = line_buf;
        while (*ptr == ' ' || *ptr == '\t') ptr++;
        if (*ptr == '#' || *ptr == '\n' || *ptr == '\r' || *ptr == '\0') {
            continue; // Skip comments and blank lines
        }

        int label_val;
        int chars_read = 0;
        if (sscanf(ptr, "%d%n", &label_val, &chars_read) != 1) {
            continue;
        }
        ptr += chars_read;

        // Standardize labels: 0 -> -1, 1 -> +1, -1 -> -1
        int final_label = (label_val == 0) ? -1 : ((label_val > 0) ? 1 : -1);

        CartSample current_sample;
        current_sample.label = final_label;
        current_sample.n_features = n_features;

        bool parse_success = true;
        for (int j = 0; j < n_features; j++) {
            double feat_val;
            if (sscanf(ptr, "%lf%n", &feat_val, &chars_read) != 1) {
                parse_success = false;
                break;
            }
            current_sample.features[j] = feat_val;
            ptr += chars_read;
        }

        if (parse_success) {
            if (count >= capacity) {
                capacity *= 2;
                CartSample *new_samples = (CartSample*)realloc(samples, capacity * sizeof(CartSample));
                if (!new_samples) {
                    fprintf(stderr, "Error: Memory realloc failed for dataset samples.\n");
                    break;
                }
                samples = new_samples;
            }
            samples[count++] = current_sample;
        }
    }

    fclose(f);

    out_dataset->samples = samples;
    out_dataset->n_samples = count;
    out_dataset->n_features = n_features;

    return count;
}

/**
 * @brief Frees all allocated memory in a CartDataset.
 */
void cart_free_dataset(CartDataset *dataset) {
    if (dataset) {
        if (dataset->samples) {
            free(dataset->samples);
            dataset->samples = NULL;
        }
        dataset->n_samples = 0;
        dataset->n_features = 0;
    }
}

/**
 * @brief Recursively frees a node and all its descendants.
 */
void cart_free_node(CartNode *node) {
    if (!node) return;
    if (node->left) cart_free_node(node->left);
    if (node->right) cart_free_node(node->right);
    free(node);
}

/**
 * @brief Frees the entire CART decision tree.
 */
void cart_free_tree(CartTree *tree) {
    if (!tree) return;
    if (tree->root) {
        cart_free_node(tree->root);
        tree->root = NULL;
    }
    free(tree);
}

/**
 * @brief Helper to serialize a node recursively.
 */
static void cart_save_node_recursive(const CartNode *node, FILE *f) {
    if (!node) {
        fprintf(f, "NULL\n");
        return;
    }
    if (node->is_leaf) {
        fprintf(f, "LEAF %d %.15e %d %d %d %.15e\n",
                node->prediction, node->probability_pos, node->n_samples,
                node->n_pos, node->n_neg, node->impurity);
    } else {
        fprintf(f, "NODE %d %.15e %d %.15e %d %d %d %.15e\n",
                node->feature_index, node->threshold, node->prediction,
                node->probability_pos, node->n_samples, node->n_pos,
                node->n_neg, node->impurity);
        cart_save_node_recursive(node->left, f);
        cart_save_node_recursive(node->right, f);
    }
}

/**
 * @brief Saves a trained CART tree to a text file.
 */
int cart_save_tree(const CartTree *tree, const char *filepath) {
    if (!tree || !filepath) return -1;
    FILE *f = fopen(filepath, "w");
    if (!f) {
        fprintf(stderr, "Error: Cannot open %s for writing model.\n", filepath);
        return -1;
    }

    fprintf(f, "CART_MODEL_V1\n");
    fprintf(f, "PARAMS %d %d %d %.15e %d %d %d\n",
            tree->params.max_depth, tree->params.min_samples_split,
            tree->params.min_samples_leaf, tree->params.min_impurity_decrease,
            tree->n_features, tree->total_nodes, tree->max_depth_reached);

    cart_save_node_recursive(tree->root, f);
    fclose(f);
    return 0;
}

/**
 * @brief Helper to deserialize a node recursively.
 */
static CartNode* cart_load_node_recursive(FILE *f) {
    char type[32];
    if (fscanf(f, "%31s", type) != 1) return NULL;

    if (strcmp(type, "NULL") == 0) {
        return NULL;
    }

    CartNode *node = (CartNode*)malloc(sizeof(CartNode));
    if (!node) return NULL;

    if (strcmp(type, "LEAF") == 0) {
        node->is_leaf = true;
        node->feature_index = -1;
        node->threshold = 0.0;
        if (fscanf(f, "%d %lf %d %d %d %lf",
                   &node->prediction, &node->probability_pos, &node->n_samples,
                   &node->n_pos, &node->n_neg, &node->impurity) != 6) {
            free(node);
            return NULL;
        }
        node->left = NULL;
        node->right = NULL;
        return node;
    } else if (strcmp(type, "NODE") == 0) {
        node->is_leaf = false;
        if (fscanf(f, "%d %lf %d %lf %d %d %d %lf",
                   &node->feature_index, &node->threshold, &node->prediction,
                   &node->probability_pos, &node->n_samples, &node->n_pos,
                   &node->n_neg, &node->impurity) != 8) {
            free(node);
            return NULL;
        }
        node->left = cart_load_node_recursive(f);
        node->right = cart_load_node_recursive(f);
        return node;
    }

    free(node);
    return NULL;
}

/**
 * @brief Loads a saved CART tree model from a text file.
 */
CartTree* cart_load_tree(const char *filepath) {
    if (!filepath) return NULL;
    FILE *f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open model file: %s\n", filepath);
        return NULL;
    }

    char header[64];
    if (fscanf(f, "%63s", header) != 1 || strcmp(header, "CART_MODEL_V1") != 0) {
        fprintf(stderr, "Error: Invalid CART model format in %s\n", filepath);
        fclose(f);
        return NULL;
    }

    CartTree *tree = (CartTree*)malloc(sizeof(CartTree));
    if (!tree) {
        fclose(f);
        return NULL;
    }

    char tag[32];
    if (fscanf(f, "%31s %d %d %d %lf %d %d %d",
               tag, &tree->params.max_depth, &tree->params.min_samples_split,
               &tree->params.min_samples_leaf, &tree->params.min_impurity_decrease,
               &tree->n_features, &tree->total_nodes, &tree->max_depth_reached) != 8) {
        free(tree);
        fclose(f);
        return NULL;
    }

    tree->root = cart_load_node_recursive(f);
    fclose(f);
    return tree;
}

/**
 * @brief End-to-end training and testing workflow for the 10D CART Decision Tree.
 * 
 * @param validation_file Path to validation dataset file (used as training set for CART).
 * @param test_file Path to test dataset file (used for final evaluation).
 * @param params Hyperparameters (or NULL for default).
 */
void cart_train_and_test(const char *validation_file, const char *test_file, const CartParams *params) {
    const char *default_val = "results/cart/validation_vectors_10d.txt";
    const char *default_test = "results/cart/test_vectors_10d.txt";

    const char *val_path = validation_file ? validation_file : default_val;
    const char *test_path = test_file ? test_file : default_test;

    printf("\n======================================================================\n");
    printf("        CART DECISION TREE FOR ANTIBODY CLASSIFICATION (10D)          \n");
    printf("======================================================================\n\n");

    CartDataset train_set;
    printf("Loading training (validation) set from: %s\n", val_path);
    int n_train = cart_load_dataset(val_path, &train_set, CART_DEFAULT_DIM);
    if (n_train <= 0) {
        fprintf(stderr, "Error: Failed to load training dataset.\n");
        return;
    }
    printf("Successfully loaded %d training samples (features: %d).\n\n",
           train_set.n_samples, train_set.n_features);

    CartDataset test_set;
    printf("Loading test set from: %s\n", test_path);
    int n_test = cart_load_dataset(test_path, &test_set, CART_DEFAULT_DIM);
    if (n_test <= 0) {
        fprintf(stderr, "Warning: Failed to load test dataset from %s.\n", test_path);
    } else {
        printf("Successfully loaded %d test samples.\n\n", test_set.n_samples);
    }

    // Default feature names corresponding to the 10 SVM kernel models
    const char *kernel_names[10] = {
        "Cauchy",
        "Gaussian (RBF)",
        "Linear",
        "Rational Quadratic",
        "Sigmoid",
        "Tanimoto",
        "Piecewise (C=0.40)",
        "Piecewise (C=0.60)",
        "Piecewise (C=0.80)",
        "Piecewise (C=0.90)"
    };

    printf("Fitting CART Decision Tree...\n");
    CartTree *tree = cart_fit(&train_set, params);
    if (!tree) {
        fprintf(stderr, "Error: Tree fitting failed.\n");
        cart_free_dataset(&train_set);
        if (n_test > 0) cart_free_dataset(&test_set);
        return;
    }
    printf("Tree fitting complete! (Total nodes: %d, Max depth: %d)\n\n",
           tree->total_nodes, tree->max_depth_reached);

    // Print the learned decision tree structure
    cart_print_tree(tree, kernel_names);

    // Evaluate on training (validation) set
    cart_print_metrics(tree, &train_set, "Validation Set (Training)");

    // Evaluate on test set
    if (n_test > 0) {
        cart_print_metrics(tree, &test_set, "Test Set (Generalization)");
    }

    // Save model to results/cart/cart_tree_model.txt
    mkdir_p("results/cart");
    const char *model_save_path = "results/cart/cart_tree_model.txt";
    if (cart_save_tree(tree, model_save_path) == 0) {
        printf("Model successfully saved to: %s\n\n", model_save_path);
    }

    // Clean up
    cart_free_tree(tree);
    cart_free_dataset(&train_set);
    if (n_test > 0) cart_free_dataset(&test_set);
}
