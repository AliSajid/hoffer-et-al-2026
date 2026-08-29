suppressPackageStartupMessages({
    library(e1071)
    library(MASS)
    library(KRSA)
    library(mixOmics)
})

datafile <- file.path("raw", "_ExportsNorm_Default_Mean_260825124141.txt")

raw_data <- krsa_read(datafile)

clean_data <- raw_data |>
    dplyr::select(SampleName, Group = `Comment 3`, Peptide, Signal) |>
    pivot_wider(names_from = Peptide, values_from = Signal)

data_matrix <- clean_data |>
    column_to_rownames("SampleName") |>
    dplyr::select(-Group) |>
    as.matrix()

label_vector <- clean_data[["Group"]]


##############################################
# PLS-DA Analysis
##############################################

pls_fit <- plsda(
    X = data_matrix,
    Y = label_vector
)


##############################################
# PCA Analysis
##############################################

pca_result <- prcomp(
    data_matrix,
    scale = TRUE
)

pca_data <- data.frame(
    Sample = rownames(pca_result[["x"]]),
    Group = as.factor(clean_data[["Group"]]),
    X = pca_result[["x"]][, 1L],
    Y = pca_result[["x"]][, 2L]
)

pca_variance <- pca_result$sdev^2L

pca_variance_percentage <- round(
    pca_variance / sum(pca_variance) * 100L,
    1L
)


##############################################
# Common plotting objects
##############################################

pca_basis <- ggplot(
    pca_data,
    aes(
        x = X,
        y = Y,
        color = Group,
        label = Sample
    )
)

pca_scales <- list(
    scale_color_brewer(
        limits = c("T", "TCL"),
        labels = c("Tumor", "Tumor Cell Line"),
        palette = "Dark2",
        name = ""
    ),
    labs(
        x = str_glue(
            "PC1 - ({pca_variance_percentage[1]} %)"
        ),
        y = str_glue(
            "PC2 - ({pca_variance_percentage[2]} %)"
        )
    ),
    theme_minimal()
)


##############################################
# Figure 1: Baseline PCA
##############################################

pca_plot <- pca_basis +
    geom_point(size = 5L) +
    pca_scales

pca_plot


############################################################
# LDA
############################################################

lda_fit <- lda(
    Group ~ X + Y,
    data = pca_data
)

lda_fit


############################################################
# Calculate LDA classification boundary
############################################################

groups <- levels(pca_data$Group)

group_1 <- pca_data |>
    filter(Group == groups[1L]) |>
    dplyr::select(X, Y)

group_2 <- pca_data |>
    filter(Group == groups[2L]) |>
    dplyr::select(X, Y)

mu_1 <- colMeans(group_1)
mu_2 <- colMeans(group_2)

n_1 <- nrow(group_1)
n_2 <- nrow(group_2)

############################################################
# Pooled within-group covariance
############################################################

S_1 <- cov(group_1)
S_2 <- cov(group_2)

pooled_covariance <- (
    (n_1 - 1L) * S_1 +
        (n_2 - 1L) * S_2
) / (n_1 + n_2 - 2L)

pooled_covariance_inverse <- solve(
    pooled_covariance
)

############################################################
# LDA coefficients
############################################################

lda_w <- pooled_covariance_inverse %*% (
    mu_1 - mu_2
)

lda_priors <- lda_fit$prior

lda_b <- (
    -0.5 *
        (
            t(mu_1) %*%
                pooled_covariance_inverse %*%
                mu_1 -
                t(mu_2) %*%
                pooled_covariance_inverse %*%
                mu_2
        )
) +
    log(lda_priors[1L] / lda_priors[2L])

lda_slope <- -lda_w[1L] / lda_w[2L]

lda_intercept <- -lda_b / lda_w[2L]


############################################################
# Figure 2: PCA + LDA boundary
############################################################

lda_plot <- pca_basis +
    geom_point(size = 5L) +
    geom_abline(
        slope = lda_slope,
        intercept = lda_intercept[1L, 1L],
        linewidth = 1L,
        linetype = "dashed",
        inherit.aes = FALSE
    ) +
    pca_scales

lda_plot


############################################################
# LDA Bootstrap
############################################################

set.seed(12345L)

n_boot <- 2000L

bootstrap_lda <- map_dfr(
    seq_len(n_boot),
    function(i) {
        # Stratified bootstrap:
        # resample each group independently
        bootstrap_data <- pca_data |>
            group_by(Group) |>
            slice_sample(
                prop = 1L,
                replace = TRUE
            ) |>
            ungroup()

        group_1_boot <- bootstrap_data |>
            filter(Group == groups[1L]) |>
            dplyr::select(X, Y)

        group_2_boot <- bootstrap_data |>
            filter(Group == groups[2L]) |>
            dplyr::select(X, Y)

        n_1_boot <- nrow(group_1_boot)
        n_2_boot <- nrow(group_2_boot)

        mu_1_boot <- colMeans(group_1_boot)
        mu_2_boot <- colMeans(group_2_boot)

        S_1_boot <- cov(group_1_boot)
        S_2_boot <- cov(group_2_boot)

        pooled_covariance_boot <- (
            (n_1_boot - 1L) * S_1_boot +
                (n_2_boot - 1L) * S_2_boot
        ) / (
            n_1_boot + n_2_boot - 2L
        )

        pooled_inverse_boot <- solve(
            pooled_covariance_boot
        )

        w_boot <- pooled_inverse_boot %*% (
            mu_1_boot - mu_2_boot
        )

        b_boot <- (
            -0.5 *
                (
                    t(mu_1_boot) %*%
                        pooled_inverse_boot %*%
                        mu_1_boot -
                        t(mu_2_boot) %*%
                        pooled_inverse_boot %*%
                        mu_2_boot
                )
        ) +
            log(
                (n_1_boot / (n_1_boot + n_2_boot)) /
                    (n_2_boot / (n_1_boot + n_2_boot))
            )

        slope_boot <- -w_boot[1L] / w_boot[2L]

        intercept_boot <- -b_boot / w_boot[2L]

        tibble(
            iteration = i,
            slope = slope_boot,
            intercept = intercept_boot
        )
    }
)


############################################################
# Bootstrap LDA confidence band
############################################################

x_grid <- seq(
    min(pca_data$X),
    max(pca_data$X),
    length.out = 200L
)

lda_bootstrap_lines <- bootstrap_lda |>
    mutate(
        line_id = row_number()
    ) |>
    select(line_id, slope, intercept) |>
    crossing(X = x_grid) |>
    mutate(
        Y = slope * X + intercept
    )

lda_bootstrap_band <- lda_bootstrap_lines |>
    group_by(X, slope, intercept) |>
    summarise(
        Y_lower = quantile(
            Y,
            probs = 0.025,
            na.rm = TRUE
        ),
        Y_upper = quantile(
            Y,
            probs = 0.975,
            na.rm = TRUE
        ),
        .groups = "drop"
    )


############################################################
# Figure 3: PCA + bootstrapped LDA
############################################################

lda_bootstrap_plot <- pca_basis +
    geom_ribbon(
        data = lda_bootstrap_band,
        aes(
            x = X,
            ymin = Y_lower,
            ymax = Y_upper
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +
    geom_abline(
        slope = lda_slope,
        intercept = lda_intercept,
        linewidth = 1L,
        linetype = "dashed",
        inherit.aes = FALSE
    ) +
    geom_point(size = 5L) +
    pca_scales

lda_bootstrap_plot


############################################################
# SVM
############################################################

svm_fit <- svm(
    Group ~ X + Y,
    data = pca_data,
    kernel = "linear",
    cost = 1L
)

svm_fit


############################################################
# Extract SVM hyperplane
############################################################

svm_coefficients <- coef(svm_fit)

svm_w <- svm_coefficients[1L:2L]

svm_b <- svm_coefficients[3L]

svm_slope <- svm_w[1L] / svm_w[2L]

svm_intercept <- svm_b / svm_w[2L]


############################################################
# Figure 4: PCA + SVM boundary
############################################################

svm_plot <- pca_basis +

    geom_point(size = 5L) +

    geom_abline(
        slope = svm_slope,
        intercept = svm_intercept,
        linewidth = 1L,
        linetype = "dashed",
        inherit.aes = FALSE
    ) +

    pca_scales

svm_plot


############################################################
# SVM Bootstrap
############################################################

bootstrap_svm <- map_dfr(
    seq_len(n_boot),
    function(i) {
        # Stratified bootstrap
        bootstrap_data <- pca_data |>
            group_by(Group) |>
            slice_sample(
                prop = 1,
                replace = TRUE
            ) |>
            ungroup()

        fit <- tryCatch(
            svm(
                Group ~ X + Y,
                data = bootstrap_data,
                kernel = "linear",
                cost = 1L
            ),
            error = function(e) NULL
        )

        if (is.null(fit)) {
            return(
                tibble(
                    iteration = i,
                    slope = NA_real_,
                    intercept = NA_real_
                )
            )
        }

        coefficients <- coef(fit)

        w <- coefficients[1L:2L]

        b <- coefficients[3L]

        tibble(
            iteration = i,
            slope = -w[1L] / w[2L],
            intercept = -b / w[2L]
        )
    }
)


############################################################
# Remove failed bootstrap fits
############################################################

bootstrap_svm <- bootstrap_svm |>
    filter(
        is.finite(slope),
        is.finite(intercept)
    )


############################################################
# Bootstrap SVM confidence band
############################################################

svm_bootstrap_lines <- bootstrap_svm |>
    crossing(X = x_grid) |>
    mutate(
        Y = slope * X + intercept
    )

svm_bootstrap_band <- svm_bootstrap_lines |>
    group_by(X) |>
    summarise(
        Y_lower = quantile(
            Y,
            probs = 0.025,
            na.rm = TRUE
        ),
        Y_upper = quantile(
            Y,
            probs = 0.975,
            na.rm = TRUE
        ),
        .groups = "drop"
    )


############################################################
# Figure 5: PCA + bootstrapped SVM
############################################################

svm_bootstrap_plot <- pca_basis +

    geom_ribbon(
        data = svm_bootstrap_band,
        aes(
            x = X,
            ymin = Y_lower,
            ymax = Y_upper
        ),
        inherit.aes = FALSE,
        alpha = 0.2
    ) +

    geom_abline(
        slope = svm_slope,
        intercept = svm_intercept,
        linewidth = 1L,
        linetype = "dashed",
        inherit.aes = FALSE
    ) +

    geom_point(size = 5L) +

    pca_scales

svm_bootstrap_plot
