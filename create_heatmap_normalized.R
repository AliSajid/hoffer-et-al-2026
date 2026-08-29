# Create a heatmap, ordered by t tests.

suppressPackageStartupMessages({
    library(tidyverse)
    library(KRSA)
    library(broom)
})

datafile <- file.path("raw", "_ExportsNorm_Default_Mean_260825124141.txt")

raw_data <- krsa_read(datafile)

clean_data <- raw_data |>
    select(SampleName, Peptide, Signal) |>
    mutate(
        Group = case_when(
            str_detect(SampleName, "HSC4") ~ "CellLine",
            str_detect(SampleName, "SAS") ~ "CellLine",
            str_detect(SampleName, "UT") ~ "CellLine",
            .default = "Tumor"
        )
    )

tested_data <- clean_data |>
    select(-SampleName) |>
    nest(.by = c(Peptide)) |>
    mutate(test = map(data, \(x) {
        tidy(t.test(Signal ~ Group, data = x))
    })) |>
    unnest_wider(test) |>
    mutate(adjusted.p.value = p.adjust(p.value)) |>
    select(-data)

figure_data <- tested_data |>
    select(Peptide, CellLine = estimate1, Tumor = estimate2, p.value, adjusted.p.value) |>
    mutate(
        FC = Tumor / CellLine, LFC = log2(FC), significant = p.value < 0.05,
        PClass = case_when(
            p.value < 0.001 ~ "P0.001",
            p.value < 0.01 ~ "P0.01",
            p.value < 0.05 ~ "P0.05",
            .default = "NS"
        )
    ) |>
    group_by(PClass) |>
    arrange(desc(p.value), desc(abs(LFC))) |>
    ungroup() |>
    select(-PClass) |>
    filter(!is.infinite(LFC))

heatmap_base <- figure_data |>
    select(Peptide, CellLine, Tumor) |>
    pivot_longer(c(CellLine, Tumor), names_to = "Group", values_to = "Signal") |>
    ggplot(aes(x = Group, y = Peptide, fill = Signal))

heatmap <- heatmap_base +
    geom_tile(width = 1L, height = 1L) +
    scale_fill_distiller(palette = "RdBu") +
    scale_y_discrete(limits = figure_data[["Peptide"]]) +
    theme_minimal() +
    theme(plot.background = element_blank(), panel.grid = element_blank())

heatmap
ggsave(
    "Tumor-v-CellLine-heatmap-MeanScaled.jpg",
    path = "figures", plot = heatmap,
    width = 7L, height = 20L, units = "in"
)

volcano_base <- figure_data |>
    select(Peptide, LFC, significant, p.value) |>
    mutate(Loading = -log10(p.value)) |>
    ggplot(aes(x = LFC, y = Loading, color = significant))

volcano <- volcano_base +
    geom_vline(xintercept = 0L, color = "grey50") +
    geom_hline(yintercept = -log10(0.05), color = "grey40", linetype = 2L) +
    geom_point() +
    scale_color_manual(breaks = c(TRUE, FALSE), values = c("red", "black")) +
    scale_x_continuous(limits = c(-2L, 2L)) +
    theme_bw() +
    labs(x = "Relative Phosphorylation", y = expression(-log[10L](p.value))) +
    guides(color = "none")

volcano
ggsave("Tumor-v-CellLine-volcanoplot-MeanScaled.jpg",
    path = "figures", plot = volcano,
    width = 5L, height = 5L, units = "in"
)
