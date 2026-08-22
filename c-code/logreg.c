#include "head.h"

/**
 * Logistic Regression Meta-Classifier
 *
 * Stacks the decision function outputs of multiple SVM kernels into a single binary classifier for antibody humanness: +1 (Human) vs -1 (Murine).
 *
 * The model learns a weight w_i for each of the 10 SVM outputs f_i(x), producing a linear combination:
 *
 *     z(x) = b + sum_i  w_i * f_i(x)
 *     P(Human | x) = sigmoid(z(x)) = 1 / (1 + exp(-z(x)))
 *
 * Training minimizes L2-regularized binary cross-entropy via gradient descent. Feature normalization (zero-mean, unit-variance) is applied using training set statistics and stored in the model for consistent inference.
 *
 * Files:
 *   results/logreg/logreg_model.txt     - saved model
 *   results/logreg/logreg_val_scores.txt - validation set scores
 *   results/logreg/logreg_test_scores.txt - test set scores
 *   results/logreg/logreg_loss_curve.txt  - training loss per epoch
 */


#define LOGREG_MAX_DIM 64
#define LOGREG_DEFAULT_DIM 10




/* --- Utility --- */

static double logreg_sigmoid(double z) {
    /* Numerically stable sigmoid */
    if (z >= 0.0) {
        double e = exp(-z);
        return 1.0 / (1.0 + e);
    } else {
        double e = exp(z);
        return e / (1.0 + e);
    }
}

/**
 * @brief Returns default training hyperparameters.
 */
LogregParams logreg_default_params(void) {
    LogregParams p;
    p.lambda = 0.01;
    p.learning_rate = 0.5;
    p.epochs = 2000;
    p.threshold = 0.5;
    return p;
}


/* --- Training --- */

/**
 * @brief Fits a logistic regression model via L2-regularized gradient descent.
 *
 * Computes feature normalization statistics from the training set, normalizes
 * all features, and runs gradient descent to minimize binary cross-entropy.
 *
 * Decision convention: label > 0 is mapped to y=1, else y=0 internally.
 * Labels are expected to be +1 (Human) or -1 (Murine).
 *
 * @param dataset       Training dataset.
 * @param params        Hyperparameters (NULL uses defaults).
 * @param out_loss_history  If non-NULL, caller-allocated array of size params->epochs
 *                          that receives per-epoch loss values.
 * @param out_n_history     Set to the number of epochs actually run.
 * @return Pointer to trained LogregModel (caller must free with logreg_free_model).
 */
LogregModel* logreg_fit(const LogregDataset *dataset, const LogregParams *params, double *out_loss_history, int *out_n_history) {
    if (!dataset || dataset->n_samples <= 0) {
        fprintf(stderr, "Error: Invalid dataset in logreg_fit.\n");
        return NULL;
    }

    LogregModel *model = (LogregModel*)calloc(1, sizeof(LogregModel));
    if (!model) {
        fprintf(stderr, "Error: Failed to allocate LogregModel.\n");
        return NULL;
    }

    model->n_features = dataset->n_features;
    model->params = params ? *params : logreg_default_params();
    const LogregParams *hp = &model->params;

    int n = dataset->n_samples;
    int d = dataset->n_features;

    /* 1. Compute normalization statistics from training data */
    for (int j = 0; j < d; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += dataset->samples[i].features[j];
        model->feature_mean[j] = sum / n;

        double sq = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = dataset->samples[i].features[j] - model->feature_mean[j];
            sq += diff * diff;
        }
        model->feature_std[j] = sqrt(sq / n) + 1e-10;
    }

    /* 2. Allocate normalized feature matrix and binary labels */
    double *X = (double*)malloc((size_t)n * (size_t)d * sizeof(double));
    double *y = (double*)malloc((size_t)n * sizeof(double));
    if (!X || !y) {
        if (X) free(X);
        if (y) free(y);
        free(model);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < d; j++) {
            X[i * d + j] = (dataset->samples[i].features[j] - model->feature_mean[j]) / model->feature_std[j];
        }
        y[i] = (dataset->samples[i].label > 0) ? 1.0 : 0.0;
    }

    /* 3. Initialize weights to zero */
    double *w = (double*)calloc(d, sizeof(double));
    double bias = 0.0;
    if (!w) {
        free(X); free(y); free(model);
        return NULL;
    }

    /* 4. Gradient descent */
    double *grad_w = (double*)malloc(d * sizeof(double));
    if (!grad_w) {
        free(X); free(y); free(w); free(model);
        return NULL;
    }

    for (int ep = 0; ep < hp->epochs; ep++) {
        /* Forward pass: compute predictions and loss */
        double loss = 0.0;
        for (int j = 0; j < d; j++) grad_w[j] = 0.0;
        double grad_b = 0.0;

        for (int i = 0; i < n; i++) {
            double z = bias;
            for (int j = 0; j < d; j++) z += w[j] * X[i * d + j];

            double p = logreg_sigmoid(z);
            double err = p - y[i];

            /* Loss: cross-entropy component */
            double p_clip = p < 1e-12 ? 1e-12 : (p > 1.0 - 1e-12 ? 1.0 - 1e-12 : p);
            loss -= (y[i] * log(p_clip) + (1.0 - y[i]) * log(1.0 - p_clip));

            /* Accumulate gradients */
            for (int j = 0; j < d; j++) {
                grad_w[j] += err * X[i * d + j];
            }
            grad_b += err;
        }

        /* Normalize by sample count, add L2 regularization to loss */
        loss /= n;
        for (int j = 0; j < d; j++) {
            loss += 0.5 * hp->lambda * w[j] * w[j];
        }

        /* Gradient step */
        for (int j = 0; j < d; j++) {
            grad_w[j] = grad_w[j] / n + hp->lambda * w[j];
            w[j] -= hp->learning_rate * grad_w[j];
        }
        grad_b /= n;
        bias -= hp->learning_rate * grad_b;

        if (out_loss_history) out_loss_history[ep] = loss;
    }

    if (out_n_history) *out_n_history = hp->epochs;

    /* 5. Store final weights in model */
    for (int j = 0; j < d; j++) model->weights[j] = w[j];
    model->bias = bias;

    free(X); free(y); free(w); free(grad_w);
    return model;
}


/* --- Inference --- */

/**
 * @brief Normalizes a raw feature vector and computes the linear logit z.
 */
static double logreg_logit(const LogregModel *model, const double *raw_features) {
    double z = model->bias;
    for (int j = 0; j < model->n_features; j++) {
        double x_norm = (raw_features[j] - model->feature_mean[j]) / model->feature_std[j];
        z += model->weights[j] * x_norm;
    }
    return z;
}

/**
 * @brief Returns +1 (Human) or -1 (Murine) for a raw 10D feature vector.
 */
int logreg_predict(const LogregModel *model, const double *raw_features) {
    if (!model || !raw_features) return 0;
    double p = logreg_sigmoid(logreg_logit(model, raw_features));
    return (p >= model->params.threshold) ? 1 : -1;
}

/**
 * @brief Returns P(Human=+1) ∈ [0,1] for a raw 10D feature vector.
 */
double logreg_predict_prob(const LogregModel *model, const double *raw_features) {
    if (!model || !raw_features) return 0.5;
    return logreg_sigmoid(logreg_logit(model, raw_features));
}


/* --- Evaluation --- */

/**
 * @brief Computes a ConfusionMatrix for a given target class (+1 or -1).
 */
ConfusionMatrix logreg_evaluate(const LogregModel *model, const LogregDataset *dataset, int target_class) {
    ConfusionMatrix cm = {0, 0, 0, 0};
    for (int i = 0; i < dataset->n_samples; i++) {
        int true_label = dataset->samples[i].label;
        int pred_label = logreg_predict(model, dataset->samples[i].features);

        if (true_label == target_class) {
            if (pred_label == target_class) cm.TP++;
            else cm.FN++;
        } else {
            if (pred_label == target_class) cm.FP++;
            else cm.TN++;
        }
    }
    return cm;
}

/**
 * @brief Prints full evaluation metrics: confusion matrix, F1 scores, accuracy.
 */
void logreg_print_metrics(const LogregModel *model, const LogregDataset *dataset, const char *name) {
    ConfusionMatrix cm_h = logreg_evaluate(model, dataset, 1);
    ConfusionMatrix cm_m = logreg_evaluate(model, dataset, -1);

    int total = dataset->n_samples;
    int correct = cm_h.TP + cm_h.TN;
    double acc = (total > 0) ? (double)correct / total : 0.0;
    double f1_h = (2.0*cm_h.TP + cm_h.FP + cm_h.FN > 0) ? 2.0*cm_h.TP / (2.0*cm_h.TP + cm_h.FP + cm_h.FN) : 0.0;
    double f1_m = (2.0*cm_m.TP + cm_m.FP + cm_m.FN > 0) ? 2.0*cm_m.TP / (2.0*cm_m.TP + cm_m.FP + cm_m.FN) : 0.0;

    printf("======================================================================\n");
    printf(" Logistic Regression Evaluation: %s\n", name ? name : "Dataset");
    printf("======================================================================\n");
    printf(" Total Samples: %d (Human [+1]: %d, Murine [-1]: %d)\n",
           total, cm_h.TP + cm_h.FN, cm_h.FP + cm_h.TN);
    printf(" Accuracy:  %.6f (%.2f%%)\n", acc, acc * 100.0);
    printf(" F1 Human:  %.6f\n", f1_h);
    printf(" F1 Murine: %.6f\n", f1_m);
    printf("----------------------------------------------------------------------\n");
    printf(" Confusion Matrix (Human = Positive [+1]):\n");
    printf("                  Predicted Human (+1)   Predicted Murine (-1)\n");
    printf("   Actual Human:       %6d (TP)                %6d (FN)\n", cm_h.TP, cm_h.FN);
    printf("   Actual Murine:      %6d (FP)                %6d (TN)\n", cm_h.FP, cm_h.TN);
    printf("======================================================================\n\n");
}

/**
 * @brief Prints the learned model weights and feature names.
 */
void logreg_print_weights(const LogregModel *model, const char **feature_names) {
    printf("----------------------------------------------------------------------\n");
    printf(" Logistic Regression Model Weights (in normalized feature space):\n");
    printf("----------------------------------------------------------------------\n");
    double total = 0.0;
    for (int j = 0; j < model->n_features; j++) total += fabs(model->weights[j]);
    for (int j = 0; j < model->n_features; j++) {
        const char *name = (feature_names && feature_names[j]) ? feature_names[j] : "feature";
        printf("   w[%2d] %-22s = %+.6f  (%.2f%% of |w|)\n",
               j, name, model->weights[j], fabs(model->weights[j]) / total * 100.0);
    }
    printf("   bias                       = %+.6f\n", model->bias);
    printf("----------------------------------------------------------------------\n\n");
}


/* --- I/O --- */

/**
 * @brief Loads a 10D vector dataset from a text file.
 * Format: lines of "<tag>\t<f0>\t<f1>\t...<f9>"
 * Lines starting with '#' are skipped.
 */
int logreg_load_dataset(const char *filepath, LogregDataset *out, int n_features) {
    if (!filepath || !out) return 0;
    if (n_features <= 0 || n_features > LOGREG_MAX_DIM) n_features = LOGREG_DEFAULT_DIM;

    FILE *f = fopen(filepath, "r");
    if (!f) {
        fprintf(stderr, "Error: Cannot open dataset: %s\n", filepath);
        return 0;
    }

    int capacity = 512, count = 0;
    LogregSample *samples = (LogregSample*)malloc(capacity * sizeof(LogregSample));
    if (!samples) { fclose(f); return 0; }

    char buf[4096];
    while (fgets(buf, sizeof(buf), f)) {
        char *ptr = buf;
        while (*ptr == ' ' || *ptr == '\t') ptr++;
        if (*ptr == '#' || *ptr == '\n' || *ptr == '\r' || *ptr == '\0') continue;

        int label, nr;
        if (sscanf(ptr, "%d%n", &label, &nr) != 1) continue;
        ptr += nr;

        /* Normalize label: 0 -> -1, 1 -> +1, -1 -> -1 */
        int final_label = (label == 0) ? -1 : ((label > 0) ? 1 : -1);

        LogregSample s;
        s.label = final_label;
        s.n_features = n_features;
        bool ok = true;
        for (int j = 0; j < n_features; j++) {
            if (sscanf(ptr, "%lf%n", &s.features[j], &nr) != 1) { ok = false; break; }
            ptr += nr;
        }
        if (!ok) continue;

        if (count >= capacity) {
            capacity *= 2;
            LogregSample *tmp = (LogregSample*)realloc(samples, capacity * sizeof(LogregSample));
            if (!tmp) break;
            samples = tmp;
        }
        samples[count++] = s;
    }

    fclose(f);
    out->samples = samples;
    out->n_samples = count;
    out->n_features = n_features;
    return count;
}

/**
 * @brief Frees memory used by a LogregDataset.
 */
void logreg_free_dataset(LogregDataset *dataset) {
    if (dataset && dataset->samples) {
        free(dataset->samples);
        dataset->samples = NULL;
        dataset->n_samples = 0;
        dataset->n_features = 0;
    }
}

/**
 * @brief Saves a LogregModel to a text file.
 */
int logreg_save_model(const LogregModel *model, const char *filepath) {
    if (!model || !filepath) return -1;
    FILE *f = fopen(filepath, "w");
    if (!f) { fprintf(stderr, "Error: Cannot open model file %s for writing.\n", filepath); return -1; }

    fprintf(f, "LOGREG_MODEL_V1\n");
    fprintf(f, "n_features %d\n", model->n_features);
    fprintf(f, "params %.15e %.15e %d %.15e\n",
            model->params.lambda, model->params.learning_rate,
            model->params.epochs, model->params.threshold);
    fprintf(f, "bias %.15e\n", model->bias);
    fprintf(f, "weights");
    for (int j = 0; j < model->n_features; j++) fprintf(f, " %.15e", model->weights[j]);
    fprintf(f, "\n");
    fprintf(f, "feature_mean");
    for (int j = 0; j < model->n_features; j++) fprintf(f, " %.15e", model->feature_mean[j]);
    fprintf(f, "\n");
    fprintf(f, "feature_std");
    for (int j = 0; j < model->n_features; j++) fprintf(f, " %.15e", model->feature_std[j]);
    fprintf(f, "\n");

    fclose(f);
    return 0;
}

/**
 * @brief Loads a LogregModel previously saved with logreg_save_model.
 */
LogregModel* logreg_load_model(const char *filepath) {
    if (!filepath) return NULL;
    FILE *f = fopen(filepath, "r");
    if (!f) { fprintf(stderr, "Error: Cannot open model file %s.\n", filepath); return NULL; }

    char header[64];
    if (fscanf(f, "%63s", header) != 1 || strcmp(header, "LOGREG_MODEL_V1") != 0) {
        fprintf(stderr, "Error: Invalid model file format.\n"); fclose(f); return NULL;
    }

    LogregModel *model = (LogregModel*)calloc(1, sizeof(LogregModel));
    if (!model) { fclose(f); return NULL; }

    char tag[64];
    fscanf(f, "%63s %d", tag, &model->n_features);
    fscanf(f, "%63s %lf %lf %d %lf", tag,
           &model->params.lambda, &model->params.learning_rate,
           &model->params.epochs, &model->params.threshold);
    fscanf(f, "%63s %lf", tag, &model->bias);
    fscanf(f, "%63s", tag);
    for (int j = 0; j < model->n_features; j++) fscanf(f, "%lf", &model->weights[j]);
    fscanf(f, "%63s", tag);
    for (int j = 0; j < model->n_features; j++) fscanf(f, "%lf", &model->feature_mean[j]);
    fscanf(f, "%63s", tag);
    for (int j = 0; j < model->n_features; j++) fscanf(f, "%lf", &model->feature_std[j]);

    fclose(f);
    return model;
}

/**
 * @brief Frees a LogregModel.
 */
void logreg_free_model(LogregModel *model) {
    if (model) free(model);
}

/**
 * @brief Writes per-sample prediction scores to a tab-delimited file.
 * Columns: tag, logit, prob, correct
 */
void logreg_save_scores(const LogregModel *model, const LogregDataset *dataset, const char *filepath) {
    FILE *f = fopen(filepath, "w");
    if (!f) { fprintf(stderr, "Error: Cannot open %s for writing.\n", filepath); return; }

    fprintf(f, "tag\tlogit\tprob\tcorrect\n");
    for (int i = 0; i < dataset->n_samples; i++) {
        int true_label = dataset->samples[i].label;
        double logit = logreg_logit(model, dataset->samples[i].features);
        double prob = logreg_sigmoid(logit);
        int pred = (prob >= model->params.threshold) ? 1 : -1;
        fprintf(f, "%d\t%.15e\t%.15e\t%d\n",
                true_label, logit, prob, (pred == true_label) ? 1 : 0);
    }
    fclose(f);
}

/**
 * @brief Writes per-epoch loss history to a tab-delimited file.
 * Columns: epoch, loss
 */
void logreg_save_loss_curve(const double *loss_history, int n, const char *filepath) {
    FILE *f = fopen(filepath, "w");
    if (!f) { fprintf(stderr, "Error: Cannot open %s for writing.\n", filepath); return; }

    fprintf(f, "epoch\tloss\n");
    for (int i = 0; i < n; i++) {
        fprintf(f, "%d\t%.15e\n", i + 1, loss_history[i]);
    }
    fclose(f);
}


/* --- Full Pipeline --- */

/**
 * @brief Full training-and-evaluation pipeline for the 10D logistic regression classifier.
 *
 * Loads the 10D SVM vector datasets, trains logistic regression, prints metrics,
 * and saves the model and score files to results/logreg/.
 *
 * @param val_file  Path to validation set (used as training data), or NULL for default.
 * @param test_file Path to test set, or NULL for default.
 * @param params    Hyperparameters, or NULL for defaults.
 */
void logreg_train_and_test(const char *val_file, const char *test_file, const LogregParams *params) {
    const char *default_val  = "results/cart/validation_vectors_10d.txt";
    const char *default_test = "results/cart/test_vectors_10d.txt";

    const char *val_path  = val_file  ? val_file  : default_val;
    const char *test_path = test_file ? test_file : default_test;

    const char *kernel_names[10] = {
        "Cauchy",
        "Gaussian (RBF)",
        "Linear",
        "Sigmoid",
        "Tanimoto",
        "Piecewise (C=0.80, W=0.20)",
        "Piecewise (C=0.80, W=0.08)",
        "Piecewise (C=0.82, W=0.08)",
        "Piecewise (C=0.84, W=0.08)",
        "Piecewise (C=0.78, W=0.08)"
    };

    printf("\n======================================================================\n");
    printf("      LOGISTIC REGRESSION META-CLASSIFIER FOR ANTIBODY CLASS          \n");
    printf("======================================================================\n\n");

    /* Load datasets */
    LogregDataset train_set;
    printf("Loading training set from: %s\n", val_path);
    int n_train = logreg_load_dataset(val_path, &train_set, LOGREG_DEFAULT_DIM);
    if (n_train <= 0) { fprintf(stderr, "Error: Failed to load training dataset.\n"); return; }
    printf("Loaded %d training samples.\n\n", n_train);

    LogregDataset test_set;
    printf("Loading test set from: %s\n", test_path);
    int n_test = logreg_load_dataset(test_path, &test_set, LOGREG_DEFAULT_DIM);
    if (n_test > 0) printf("Loaded %d test samples.\n\n", n_test);
    else fprintf(stderr, "Warning: Could not load test dataset.\n\n");

    /* Set up params */
    LogregParams hp = params ? *params : logreg_default_params();

    /* Allocate loss history buffer */
    double *loss_history = (double*)malloc(hp.epochs * sizeof(double));
    int n_epochs_run = 0;

    /* Train */
    printf("Fitting Logistic Regression (lambda=%.4f, lr=%.4f, epochs=%d)...\n",
           hp.lambda, hp.learning_rate, hp.epochs);
    LogregModel *model = logreg_fit(&train_set, &hp, loss_history, &n_epochs_run);
    if (!model) {
        fprintf(stderr, "Error: Logistic regression fitting failed.\n");
        if (loss_history) free(loss_history);
        logreg_free_dataset(&train_set);
        if (n_test > 0) logreg_free_dataset(&test_set);
        return;
    }
    printf("Training complete!\n\n");

    /* Print weights */
    logreg_print_weights(model, kernel_names);

    /* Evaluate */
    logreg_print_metrics(model, &train_set, "Validation Set (Training)");
    if (n_test > 0) logreg_print_metrics(model, &test_set, "Test Set (Generalization)");

    /* Save outputs */
    mkdir_p("results/logreg");
    logreg_save_model(model, "results/logreg/logreg_model.txt");
    printf("Model saved to results/logreg/logreg_model.txt\n");

    logreg_save_scores(model, &train_set, "results/logreg/logreg_val_scores.txt");
    printf("Validation scores saved to results/logreg/logreg_val_scores.txt\n");

    if (n_test > 0) {
        logreg_save_scores(model, &test_set, "results/logreg/logreg_test_scores.txt");
        printf("Test scores saved to results/logreg/logreg_test_scores.txt\n");
    }

    if (loss_history) {
        logreg_save_loss_curve(loss_history, n_epochs_run, "results/logreg/logreg_loss_curve.txt");
        printf("Loss curve saved to results/logreg/logreg_loss_curve.txt\n");
        free(loss_history);
    }

    printf("\n");

    /* Clean up */
    logreg_free_model(model);
    logreg_free_dataset(&train_set);
    if (n_test > 0) logreg_free_dataset(&test_set);
}
