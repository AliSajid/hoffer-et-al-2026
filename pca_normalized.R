suppressPackageStartupMessages({
    library(e1071)
    library(KRSA)
    library(tidyverse)
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
    geom_abline(
        intercept = pca_boundry_coefs[1L] / -pca_boundry_coefs[3L],
        slope = pca_boundry_coefs[2L] / -pca_boundry_coefs[3L],
        linetype = 2L,
        linewidth = 2L,
        color = "grey80"
    ) +
    geom_point(size = 3L) +
    scale_color_manual(
        values = c("coral3", "blue3"),
        limits = c("T", "TCL"),
        labels = c("Tumor", "Tumor Cell Line"),
        name = ""
    ) +
    labs(
        x = str_glue("Principal Component 1 - ({pca_variance_percentage[1]} %)"),
        y = str_glue("Principal Component 2 - ({pca_variance_percentage[2]} %)")
    ) +
    guides(color = guide_legend(position = "bottom")) +
    theme_bw() +
    theme(text = element_text(size = 20))
pca_plot

ggsave("Tumor_v_CellLine_PCA_plot_MeanScaled.jpg",
    width = 7L, height = 7L,
    units = "in",
    path = "figures",
    plot = pca_plot
)
ggsave("Tumor_v_CellLine_PCA_plot_MeanScaled.svg",
    width = 7L, height = 7L,
    units = "in",
    path = "figures",
    plot = pca_plot
)
