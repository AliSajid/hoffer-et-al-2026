suppressPackageStartupMessages({
    library(tidyverse)
    library(e1071)
    library(KRSA)
    library(mixOmics)
})

datafile <- file.path("raw", "_ExportsNorm_Default_Mean_260825124141.txt")

raw_data <- krsa_read(datafile)

clean_data <- raw_data |>
    select(SampleName, Group = `Comment 3`, Peptide, Signal) |>
    pivot_wider(names_from = Peptide, values_from = Signal)

data_matrix <- clean_data |>
    column_to_rownames("SampleName") |>
    select(-Group) |>
    as.matrix()

label_vector <- clean_data[["Group"]]

##############################################
# PLS-DA Analysis
##############################################

pls_fit <- plsda(X = data_matrix, Y = label_vector)

##############################################
# PCA Analysis
##############################################

pca_result <- prcomp(data_matrix, scale = TRUE)

pca_data <- data.frame(
    Sample = rownames(pca_result[["x"]]),
    Group = as.factor(clean_data[["Group"]]),
    X = pca_result[["x"]][, 1L],
    Y = pca_result[["x"]][, 2L]
)

pca_variance <- pca_result$sdev^2L
pca_variance_percentage <- round(pca_variance / sum(pca_variance) * 100L, 1L)

pca_boundary <- svm(
    Group ~ X + Y,
    data = pca_data,
    kernel = "linear",
    cost = 1L
)

pca_boundry_coefs <- coef(pca_boundary)

pca_basis <- ggplot(pca_data, aes(x = X, y = Y, color = Group, label = Sample))

pca_plot <- pca_basis +
    geom_point(size = 5L) +
    geom_abline(
        intercept = pca_boundry_coefs[1L] / -pca_boundry_coefs[3],
        slope = pca_boundry_coefs[2] / -pca_boundry_coefs[3]
    ) +
    scale_color_brewer(
        limits = c("T", "TCL"),
        labels = c("Tumor", "Tumor Cell Line"),
        palette = "Dark2", name = ""
    ) +
    labs(
        x = str_glue("PC1 - ({pca_variance_percentage[1]} %)"),
        y = str_glue("PC2 - ({pca_variance_percentage[2]} %)")
    ) +
    theme_minimal()
pca_plot
