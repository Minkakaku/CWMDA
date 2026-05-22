# go_data:orderd GO/KEGG Matrix
set <- list(1)
# set: make one X value
p <- ggplot(go_data[1:10, ], aes(set, reorder(Description, -pvalue))) +
    geom_point(aes(size = count[1:10], color = -1 * log2(pvalue))) +
    # low = "#AB05F9", high = "#F97805"
    scale_colour_gradient(low = "#059CFE", high = "#FE05EB") +
    # ,breaks = c(200,230,260)
    labs(color = "-log2(p)", size = "Gene number", x = "", y = "") +
    theme(
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(
            face = "bold", size = 9,
            colour = "black", lineheight = 0.7
        )
    ) +
    theme(
        legend.position = c(1, 0),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.justification = c(0.1, 1),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8, face = "bold"),
        plot.margin = margin(0.3, 0.3, 1.2, 0.3, "cm")
    ) +
    guides(
        color = guide_colorbar(
            barheight = 0.5, barwidth = 5,
            title.position = "top", order = 1
        ),
        size = guide_legend(
            override.aes = list(shape = 21), title.position = "top",
            label.position = "bottom", order = 0
        )
    ) +
    scale_size(range = c(3, 5)) + # , breaks = c(5, 10, 20)) #
    scale_y_discrete(
        position = "right",
        labels = function(df) str_wrap(df, width = 40)
    ) +
    coord_fixed(ratio = 1, xlim = c(0.1, 1.4))
ggsave(p,
    filename = paste0(data.name, ".pdf"), dpi = 600, width = 12,
    height = 12, units = c("cm")
)
