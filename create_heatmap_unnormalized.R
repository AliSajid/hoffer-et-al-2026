# Tumor vs CellLine heatmap and volcano plot
#
# Heatmap:
#   - ComplexHeatmap
#   - Rows split by PClass
#   - PClass shown as a row annotation
#   - Descriptive PClass labels
#   - Peptides ordered by p-value within each PClass

suppressPackageStartupMessages({
    library(tidyverse)
    library(KRSA)
    library(broom)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})


# -------------------------------------------------------------------------
# Data import
# -------------------------------------------------------------------------

datafile <- file.path(
    "raw",
    "_Export_Default_Mean_260827153505.txt"
)

raw_data <- krsa_read(datafile)


# -------------------------------------------------------------------------
# Clean data
# -------------------------------------------------------------------------

clean_data <- select(
    raw_data,
    SampleName,
    Peptide,
    Signal,
    Group = `Comment 3`
)

# -------------------------------------------------------------------------
# Statistical testing
# -------------------------------------------------------------------------

tested_data <- clean_data |>
    select(-SampleName) |>
    nest(.by = Peptide) |>
    mutate(
        test = map(data, \(x) {
            tidy(t.test(Signal ~ Group, data = x))
        })
    ) |>
    unnest_wider(test) |>
    mutate(
        adjusted.p.value = p.adjust(p.value)
    ) |>
    select(-data)


# -------------------------------------------------------------------------
# Prepare figure data
# -------------------------------------------------------------------------

figure_data <- tested_data |>
    select(
        Peptide,
        CellLine = estimate1,
        Tumor = estimate2,
        p.value,
        adjusted.p.value
    ) |>
    mutate(
        FC = Tumor / CellLine,
        LFC = log2(FC),
        significant = p.value < 0.05,
        PClass = case_when(
            p.value < 0.001 ~ "P0.001",
            p.value < 0.01 ~ "P0.01",
            p.value < 0.05 ~ "P0.05",
            .default = "NS"
        )
    ) |>
    filter(is.finite(LFC)) |>
    mutate(
        PClass = factor(
            PClass,
            levels = c(
                "P0.001",
                "P0.01",
                "P0.05",
                "NS"
            )
        )
    ) |>
    group_by(PClass) |>
    arrange(
        p.value,
        desc(abs(LFC)),
        .by_group = TRUE
    ) |>
    ungroup()


# -------------------------------------------------------------------------
# Heatmap matrix
# -------------------------------------------------------------------------

heatmap_matrix <- figure_data |>
    select(
        Peptide,
        CellLine,
        Tumor
    ) |>
    column_to_rownames("Peptide") |>
    as.matrix()


# -------------------------------------------------------------------------
# Row metadata
# -------------------------------------------------------------------------

row_pclass <- figure_data |>
    select(
        Peptide,
        PClass
    ) |>
    column_to_rownames("Peptide")

row_pclass <- row_pclass[
    rownames(heatmap_matrix), ,
    drop = FALSE
]


# -------------------------------------------------------------------------
# P-value labels
# -------------------------------------------------------------------------

pclass_labels <- c(
    P0.001 = "p < 0.001",
    P0.01  = "p < 0.01",
    P0.05  = "p < 0.05",
    NS     = "p ≥ 0.05"
)


# -------------------------------------------------------------------------
# P-value colors
# -------------------------------------------------------------------------

pclass_col <- c(
    P0.001 = "#7F0000",
    P0.01  = "#D7301F",
    P0.05  = "#FC8D59",
    NS     = "grey70"
)


# -------------------------------------------------------------------------
# Heatmap colors
# -------------------------------------------------------------------------

signal_range <- range(
    heatmap_matrix,
    na.rm = TRUE
)

signal_col <- colorRamp2(
    c(
        signal_range[1L],
        mean(signal_range),
        signal_range[2L]
    ),
    c(
        "#2166AC",
        "white",
        "#B2182B"
    )
)


# -------------------------------------------------------------------------
# P-value row annotation
# -------------------------------------------------------------------------

row_annotation <- rowAnnotation(
    PClass = anno_block(
        gp = gpar(
            fill = pclass_col,
            col = NA
        )
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    annotation_name_gp = gpar(
        fontsize = 10L,
        fontface = "bold"
    )
)
# -------------------------------------------------------------------------
# Row split
#
# Use the descriptive labels for the actual slice titles while retaining
# the original PClass values for the annotation.
# -------------------------------------------------------------------------

row_split <- factor(
    row_pclass$PClass,
    levels = c(
        "P0.001",
        "P0.01",
        "P0.05",
        "NS"
    ),
    labels = c(
        "p < 0.001",
        "p < 0.01",
        "p < 0.05",
        "p ≥ 0.05"
    )
)
# -------------------------------------------------------------------------
# ComplexHeatmap
# -------------------------------------------------------------------------

heatmap <- Heatmap(
    heatmap_matrix,
    name = "Signal",
    col = signal_col,
    row_split = row_split,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    row_gap = unit(1.5, "mm"),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8L),
    column_names_side = "bottom",
    column_names_gp = gpar(
        fontsize = 10L,
        fontface = "bold"
    ),
    column_title_rot = 90L,
    row_title_rot = 0L,
    row_title_gp = gpar(
        fontsize = 11L,
        fontface = "bold"
    ),
    rect_gp = gpar(
        col = "white",
        lwd = 1L
    )
)

# -------------------------------------------------------------------------
# Combine annotation and heatmap
# -------------------------------------------------------------------------

ht <- row_annotation + heatmap


# -------------------------------------------------------------------------
# Display heatmap
# -------------------------------------------------------------------------

draw(
    ht,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
)


# -------------------------------------------------------------------------
# Save heatmap
# -------------------------------------------------------------------------

pdf(
    file = file.path(
        "figures",
        "Tumor-v-CellLine-heatmap-MeanUnscaled.pdf"
    ),
    width = 7L,
    height = 10L
)

draw(
    ht,
    merge_legends = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "left"
)

dev.off()


# -------------------------------------------------------------------------
# Volcano plot
# -------------------------------------------------------------------------

volcano_base <- figure_data |>
    select(
        Peptide,
        LFC,
        significant,
        p.value
    ) |>
    mutate(
        Loading = -log10(p.value)
    ) |>
    ggplot(
        aes(
            x = LFC,
            y = Loading,
            color = significant
        )
    )


volcano <- volcano_base +
    geom_vline(
        xintercept = 0L,
        color = "grey50"
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        color = "grey40",
        linetype = 2L
    ) +
    geom_point(size = 3L) +
    scale_color_manual(
        breaks = c(TRUE, FALSE),
        values = c(
            "red",
            "black"
        )
    ) +
    scale_x_continuous(
        limits = c(-2L, 2L)
    ) +
    theme_bw() +
    labs(
        x = "Relative Phosphorylation",
        y = expression(-log[10L](p.value))
    ) +
    guides(
        color = "none"
    )


# -------------------------------------------------------------------------
# Display volcano plot
# -------------------------------------------------------------------------

volcano


# -------------------------------------------------------------------------
# Save volcano plot
# -------------------------------------------------------------------------

ggsave(
    "Tumor-v-CellLine-volcanoplot-MeanUnscaled.jpg",
    path = "figures",
    plot = volcano,
    width = 7L,
    height = 7L,
    units = "in"
)

ggsave(
    "Tumor-v-CellLine-volcanoplot-MeanUnscaled.svg",
    path = "figures",
    plot = volcano,
    width = 7L,
    height = 7L,
    units = "in"
)
