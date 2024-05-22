##### --Libraries --#####
library(ggplot2)
library(devtools)
install_github("galelab/AnchorNCounterNorm")
library(AnchorNCounterNorm)
library(reshape2)
library(stringr)

##### -- Functions --#####
generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

##### -- Count -- #####
resultsqc <- "1.initial_qc"
generate_folder(resultsqc)

data <- load_nCounter_files(
    pathtoRCC = "./rccfiles/",
    meta.data = "./Targets-file.csv",
    save.fig = TRUE, output_dir = resultsqc
)

data$p2 + scale_x_discrete(labels = data$meta.data$Project.name)
ggsave(file.path(resultsqc, "barplot_log2rawcountsProject.name.png"), width = 4.5, height = 3, bg = "white", dpi = 300)


HKstats <- hk_gene_stats(data, min_num_hk_genes = 10, group.by = "day", output_dir = resultsqc)

##### -- Normalization -- #####

normdata <- ratio_normalization(data, hkgenes = c("HPRT1", "PGK1", "ABCF1", "SDHA"), output_dir = resultsqc)
saveRDS(normdata, "normdata.RDS")
saveRDS(data, "data.RDS")

norm_matrix_melt <- reshape2::melt(normdata$log_counts_ratio)
facetvariable <- c()
for (s in norm_matrix_melt$variable) {
    facetvariable <- c(facetvariable, data$meta.data[data$meta.data$sample==s, "Project.name"])
}
norm_matrix_melt$facetvariable <- facetvariable
pl <- ggplot(norm_matrix_melt, aes(x = variable, y = value)) +
    geom_boxplot(fill = "#ff615d") + facet_wrap(~facetvariable, scales="free_x", ncol=2) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        axis.title.x = element_blank()
    )
ggsave(file.path(resultsqc, "boxplotnormlog2Project.name.png"), width = 6.5, height = 3, bg = "white", dpi = 300)
ggsave(file.path(resultsqc, "boxplotnormlog2Project.name.pdf"), width = 6.5, height = 3, bg = "white", dpi = 300)

facetvariable <- c()
for (s in norm_matrix_melt$variable) {
    facetvariable <- c(facetvariable, data$meta.data[data$meta.data$sample == s, "Treatment"])
}
norm_matrix_melt$facetvariable <- facetvariable
pl <- ggplot(norm_matrix_melt, aes(x = variable, y = value)) +
    geom_boxplot(fill = "#ff615d") +
    facet_wrap(~facetvariable, scales = "free_x", ncol = 2) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        axis.title.x = element_blank()
    )
ggsave(file.path(resultsqc, "boxplotnormlog2Treatment.png"), width = 6.5, height = 6, bg = "white", dpi = 300)
ggsave(file.path(resultsqc, "boxplotnormlog2Treatment.pdf"), width = 6.5, height = 6, bg = "white", dpi = 300)

facetvariable <- c()
for (s in norm_matrix_melt$variable) {
    facetvariable <- c(facetvariable, data$meta.data[data$meta.data$sample == s, "Sex"])
}
norm_matrix_melt$facetvariable <- facetvariable
pl <- ggplot(norm_matrix_melt, aes(x = variable, y = value)) +
    geom_boxplot(fill = "#ff615d") +
    facet_wrap(~facetvariable, scales = "free", ncol = 2) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        axis.title.x = element_blank()
    )
ggsave(file.path(resultsqc, "boxplotnormlog2Sex.png"), width = 6.5, height = 6, bg = "white", dpi = 300)
ggsave(file.path(resultsqc, "boxplotnormlog2Sex.pdf"), width = 6.5, height = 6, bg = "white", dpi = 300)

facetvariable <- c()
for (s in norm_matrix_melt$variable) {
    facetvariable <- c(facetvariable,
         paste(data$meta.data[data$meta.data$sample == s, "Treatment"], 
            data$meta.data[data$meta.data$sample == s, "Time.Course"], sep="_"))
}
norm_matrix_melt$facetvariable <- facetvariable
pl <- ggplot(norm_matrix_melt, aes(x = variable, y = value)) +
    geom_boxplot(fill = "#ff615d") +
    facet_wrap(~facetvariable, scales = "free", ncol = 4) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3),
        axis.title.x = element_blank()
    )
ggsave(file.path(resultsqc, "boxplotnormlog2TreatmentTime.png"), width = 8.5, height = 10, bg = "white", dpi = 300)
ggsave(file.path(resultsqc, "boxplotnormlog2TreatmentTime.pdf"), width = 8.5, height = 10, bg = "white", dpi = 300)

data$meta.data$TreatmentTime <- paste(data$meta.data$Treatment, data$meta.data$Time.Course, sep="_")

##### -- Dimensionality Reduction -- #####

resultsdim <- "1.resultsdim"
generate_folder(resultsdim)


dimred <- dim_reduction(normdata$log_counts_ratio, data$meta.data,
    target_columns = c(5, 6), output_dir = resultsdim, file.name = "Treatment-Time", 
    ordertargetcolumn2 = c("W-3", "D-14", "D-13", "D-12",  "D2", "D4", "D7", "D10", "D14", "D21"), 
    colorgradient = TRUE
)
target <- data$meta.data
target$Time.Course <- str_replace_all(target$Time.Course, "D-12", "D-14")
target$Time.Course <- str_replace_all(target$Time.Course, "D-13", "D-14")

dimred <- dim_reduction(normdata$log_counts_ratio, target,
    target_columns = c(3, 6), output_dir = resultsdim, file.name = "Sex-Time", 
    ordertargetcolumn2 = c("W-3", "D-14", "D2", "D4", "D7", "D10", "D14", "D21"), 
    colorgradient = TRUE
)


dimred <- dim_reduction(normdata$log_counts_ratio, target,
    target_columns = c(3, 4), output_dir = resultsdim, file.name = "Sex-Animal", 
)
# dimred <- dim_reduction(normdata$log_counts_ratio, data$meta.data,
#     target_columns = c(11, 3), output_dir = resultsdim, file.name = "Sex-TreatmentTime"
# )

# dimred <- dim_reduction(normdata$log_counts_ratio, data$meta.data,
#     target_columns = c(10, 3), output_dir = resultsdim, file.name = "Sex-Group"
# )
# targetsub <- data$meta.data[data$meta.data$Treatment == "preZIKV", ]
# subnorm <- normdata$log_counts_ratio[, targetsub$sample]

# dimred <- dim_reduction(subnorm, targetsub,
#     target_columns = c(3, 6), output_dir = resultsdim, file.name = "preZIKA-sex_time"
# )

# dimred <- dim_reduction(subnorm, targetsub,
#     target_columns = c(4,6), output_dir = resultsdim, file.name = "preZIKA-AnimalID_time"
# )

# #remove animal A11230
# targetsub <- targetsub[targetsub$Animal.ID != "A11230", ]
# subnorm <- subnorm[, targetsub$sample]

# dimred <- dim_reduction(subnorm, targetsub,
#     target_columns = c(4, 6), output_dir = resultsdim, file.name = "preZIKA-AnimalID_timeremovA11230"
# )

##### -- DE -- #####
resultsde <- "1.DE_groups1&2_noadj_0.01"
generate_folder(resultsde)
DEgroup1and2nonadj0.01 <- run_DE_analysis(normdata$log_counts_ratio, data$meta.data,
    compare.column = "grouptime", pval.cutoff = 0.01, adj.pval.method = "none",
    contrastslist = c(
        "Group 1_D_14 - Group 1_W_3", "Group 1_D2 - Group 1_W_3", 
        "Group 1_D4 - Group 1_W_3", "Group 1_D7 - Group 1_W_3", "Group 1_D10 - Group 1_W_3",
        "Group 1_D14 - Group 1_W_3", "Group 1_D21 - Group 1_W_3",
        "Group 2_D2 - Group 2_D_14", "Group 2_D4 - Group 2_D_14",
        "Group 2_D7 - Group 2_D_14", "Group 2_D10 - Group 2_D_14",
        "Group 2_D14 - Group 2_D_14", "Group 2_D21 - Group 2_D_14"
    ), DE.test = "ttest", output_dir = resultsde, pval.cuttof.for.naming.vp=0.01
) # no sig genes


resultsde <- "1.DE_groups1vs2_tps_nonadj_0.01"
generate_folder(resultsde)
DEgroup1and2 <- run_DE_analysis(normdata$log_counts_ratio, data$meta.data,
    compare.column = "grouptime", pval.cutoff = 0.01, adj.pval.method = "none",
    contrastslist = c(
        "Group 1_D_14 - Group 2_D_14", "Group 1_D2 - Group 2_D2", 
        "Group 1_D4 -Group 2_D4", "Group 1_D7 - Group 2_D7", "Group 1_D10 - Group 2_D10",
        "Group 1_D14 - Group 2_D14", "Group 1_D21 - Group 2_D21"
    ), DE.test = "ttest", output_dir = resultsde, pval.cuttof.for.naming.vp=0.01
) 

#####SIG PLOTS FOR NON ADJUSTED PVAL 0.01 CUTOFF

results_folder <- "1.SigLinePlots_nonAdjPval"
generate_folder(results_folder)
normsig <- normdata$counts_ratio[unique(DEgroup1and2nonadj0.01$sigresults$gene), ]
normsig$genes <- rownames(normsig)
normsigmelt <- data.table::melt(normsig)
time <- c()
group <- c()
for (v in normsigmelt$variable) {
    time <- c(time, data$meta.data[data$meta.data$sample==v, "Time.Course"])
    group <- c(group, data$meta.data[data$meta.data$sample == v, "Group"])
}
time <- str_replace_all(time, "D\\-13", "D\\-14")
time <- str_replace_all(time, "D\\-12", "D\\-14")
group <- str_replace_all(group, "Group 1", "SIV+ & ZIKV+")
group <- str_replace_all(group, "Group 2", "ZIKV+")

normsigmelt$time <- time
normsigmelt$group <- group
normsigmelt$time <- factor(normsigmelt$time, 
    levels=c("W-3", "D-14", "D2", "D4", "D7", "D-12", "D-13", "D10", "D14", "D21"))

normsigmelt$genes <- factor(normsigmelt$genes, levels=c(
    "ISG15", "MX2","IFIT1", "IFI44", "MX1", "IFIT2", "IL1B", "RSAD2", 
    "OAS1","DDX58","DHX58","IFIH1", "TNF", "TRIM25", "ISG20", "GBP2",
    "IFNL3", "IRF7", "NLRP3","IRF2","IFNA4", "CXCL10",  "IFITM1"))

for (g in unique(normsigmelt$genes)) {
    df <- normsigmelt[normsigmelt$genes == g, ]
    ggplot(normsigmelt, aes(x=time, y=value, fill=group, color=group)) +
        theme_minimal(base_size = 8) +
        stat_summary(aes(y = value, group = group), fun.data = "mean_se", geom = "ribbon", alpha=0.4, colour = NA) +
        stat_summary(aes(y = value, group = group), fun = "mean", geom = "line", linewidth=1) +
 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), 
        axis.title.x=element_blank(), legend.key.size = unit(0.3, "cm"), legend.text = element_text(size=6),
        legend.title =element_blank()) +labs(y="norm ratios")
    ggsave(file.path(results_folder, paste0(g, ".png")), width=2.5, height=1, dpi=300,bg="white")
}

ggplot(normsigmelt, aes(x = time, y = value, fill = group, color = group)) +
    theme_minimal(base_size = 12) + facet_wrap(~genes, scales="free_y", ncol=5) +
    stat_summary(aes(y = value, group = group), fun.data = "mean_se", geom = "ribbon", alpha = 0.4, colour = NA) +
    stat_summary(aes(y = value, group = group), fun = "mean", geom = "line", linewidth = 1) +
    scale_color_manual(values = c("SIV+ & ZIKV+" = "red", "ZIKV+" = "black")) +
    scale_fill_manual(values = c("SIV+ & ZIKV+" = "red", "ZIKV+" = "black")) +
    labs(y="Normalized Counts") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.title.x = element_blank(), legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 6),
        legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal"
    )
ggsave(file.path(results_folder, "all.png"), width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(results_folder, "all.pdf"), width = 8, height = 8, dpi = 300, bg = "white")

normsigmeltlog <- normsigmelt
normsigmeltlog$value <- log2(normsigmeltlog$value)
ggplot(normsigmeltlog, aes(x = time, y = value, fill = group, color = group)) +
    theme_minimal(base_size = 12) + facet_wrap(~genes, scales="free_y", ncol=5) +
    stat_summary(aes(y = value, group = group), fun.data = "mean_se", geom = "ribbon", alpha = 0.4, colour = NA) +
    stat_summary(aes(y = value, group = group), fun = "mean", geom = "line", linewidth = 1) +
    scale_color_manual(values = c("SIV+ & ZIKV+" = "red", "ZIKV+" = "black")) +
    scale_fill_manual(values = c("SIV+ & ZIKV+" = "red", "ZIKV+" = "black")) +
    labs(y="log2(Normalized Counts)") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.title.x = element_blank(), legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 6),
        legend.title = element_blank(), legend.position = "top", legend.direction = "horizontal"
    )
ggsave(file.path(results_folder, "alllog.png"), width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(results_folder, "alllog.pdf"), width = 8, height = 8, dpi = 300, bg = "white")


