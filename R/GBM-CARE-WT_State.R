source("R/GBM-CARE-WT-State_utils.R")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Generate NULL distribution for cell cycle and state classification
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

mdata <- as_tibble(meta_data) %>%
    filter(CellID %in% MP_scores$CellID)
dim(mdata)

mdata <- mdata %>%
    left_join(
        MP_scores %>%
            dplyr::select(CellID, starts_with("MP_"), -MP_16_GlioNeural),
        by = "CellID"
    )

sigs <- MP_list_named[-length(MP_list_named)]

permuted_data <- generate_null_dist(
    umi_data_list = umi_data_all,
    md = mdata,
    sigs = sigs,
    n_iter = 20,
    n_cells = 5000,
    verbose = T
)

saveRDS(permuted_data, paste0(DATA_ROOT, "state_class_permuted_data.RDS"))

state_programs <- sigs

permuted_data_merged <- permuted_data

vars <- names(state_programs)

scores_nd <- lapply(colnames(permuted_data_merged), function(mp) {
    x <- permuted_data_merged %>% pull(mp)

    fit <- MASS::fitdistr(x, "normal")
    class(fit)

    para <- fit$estimate

    tibble(MP = mp, Mean = para[1], SD = para[2])
})
scores_nd <- do.call(rbind, scores_nd)

mean_vec <- setNames(scores_nd$Mean, scores_nd$MP)
sd_vec <- setNames(scores_nd$SD, scores_nd$MP)

norm_fit <- lapply(colnames(permuted_data_merged), function(mp) {
    tibble(
        MP = mp,
        Sig = rnorm(n = 100000, mean = mean_vec[mp], sd = sd_vec[mp])
    )
})
norm_fit <- do.call(rbind, norm_fit)

dm <- rbind(
    melt(data = permuted_data_merged, measure.vars = vars) %>%
        mutate(DataType = "Permuted"),
    melt(data = mdata %>% select(starts_with("MP")), measure.vars = vars) %>%
        mutate(DataType = "Actual"),
    norm_fit %>%
        rename(variable = MP, value = Sig) %>%
        mutate(DataType = "Classifier")
) %>%
    mutate(DataType = factor(DataType, c("Permuted", "Actual", "Classifier")))

dm_stats <- dm %>%
    group_by(variable) %>%
    filter(DataType == "Classifier") %>%
    summarise(Q95 = quantile(value, .95), Q99 = quantile(value, .99))

ggplot(
    data = dm,
    aes(
        x = value,
        y = after_stat(ncount),
        color = DataType,
        linetype = DataType
    )
) +
    facet_wrap(~variable, scales = "free_x", nrow = 3) +
    geom_freqpoly(
        bins = 100,
        size = 1,
        show.legend = c(color = T, linetype = F)
    ) +
    scale_color_manual(
        name = "",
        values = c(
            "Permuted" = "dodgerblue",
            "Actual" = "red",
            "Classifier" = "black"
        )
    ) +
    scale_linetype_manual(
        values = c(
            "Permuted" = "solid",
            "Actual" = "solid",
            "Classifier" = "dashed"
        )
    ) +
    scale_fill_discrete(name = "Data distribution") +
    geom_vline(
        data = dm_stats,
        mapping = aes(xintercept = Q99),
        linetype = "dashed",
        size = 1
    ) +
    xlab("Program score") +
    ylab("Count (scaled to 1)") +
    theme_gbm_pvsr(panel.spacing = unit(1, "lines")) +
    theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Cell cycle classification
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

cc_data <- mdata %>% select(CellID, Sample, Patient, Timepoint, MP_3_CC)

mn <- mean_vec["MP_3_CC"]
sd <- sd_vec["MP_3_CC"]

cc_data$Z_score <- sapply(1:nrow(cc_data), function(i) {
    (cc_data$MP_3_CC[i] - mn) / sd
})

cc_data$p.val <- sapply(1:nrow(cc_data), function(i) {
    pnorm(cc_data$MP_3_CC[i], mean = mn, sd = sd, lower.tail = F)
})

cc_data$p.adj <- p.adjust(cc_data$p.val, "fdr")

cc_data$p.sig <- cc_data$p.adj < .05

cc_data$isCC <- cc_data$p.sig == T

table(cc_data$Timepoint, cc_data$isCC)
table(cc_data$Timepoint, cc_data$isCC) /
    rowSums(table(cc_data$Timepoint, cc_data$isCC))

cc_data <- cc_data %>%
    mutate(
        CCstate = case_when(
            p.adj < .05 ~ "Cycling",
            p.val < .05 ~ "Borderline",
            TRUE ~ "Non-cycling"
        )
    ) %>%
    mutate(CCstate = factor(CCstate, c("Non-cycling", "Borderline", "Cycling")))

table(cc_data$Timepoint, cc_data$CCstate)
table(cc_data$Timepoint, cc_data$CCstate) /
    rowSums(table(cc_data$Timepoint, cc_data$CCstate))

ggplot(cc_data, aes(x = Timepoint, y = MP_3_CC, color = CCstate)) +
    ggdist::stat_halfeye(
        adjust = .5,
        width = .75,
        justification = -.2,
        .width = 0,
        point_colour = NA
    ) +
    geom_boxplot(width = .2, outlier.color = NA) +
    scale_color_manual(
        name = "",
        values = c(
            "Non-cycling" = "dodgerblue",
            "Borderline" = "black",
            "Cycling" = "red"
        )
    ) +
    coord_cartesian(xlim = c(1.2, NA)) +
    xlab("") +
    ylab("Cell cycle score") +
    scale_y_continuous(
        breaks = seq(-.5, 1.5, .25),
        limits = c(-.5, 1.5),
        oob = squish
    ) +
    theme_gbm_pvsr() +
    theme(panel.grid.major = element_line())

cc_data$log10Padj <- -log10(cc_data$p.adj)

ggscatter(
    data = cc_data,
    x = "MP_3_CC",
    y = "log10Padj",
    color = "CCstate",
    facet.by = "Timepoint"
) +
    theme_gbm_pvsr() +
    theme(panel.grid.major = element_line())

mdata$isCC <- cc_data$isCC
table(mdata$isCC)

mdata$CCstate <- cc_data$CCstate
table(mdata$CCstate)

mdata$isCCplot <- ifelse(mdata$isCC == TRUE, "Cycling", "Non-cycling")
mdata$isCCplot <- as.factor(mdata$isCCplot)

which(is.na(mdata$isCC))
which(is.na(mdata$isCCplot))
which(is.na(mdata$CCstate))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Cell state classification
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

state_vars <- names(state_programs)[grep("MP_", names(state_programs))]
state_vars <- state_vars[state_vars %ni% c("MP_3_CC")]

state_data <- melt(
    data = mdata %>%
        select(CellID, state_vars),
    id.vars = "CellID",
    variable.name = "Program",
    value.name = "Score",
    measure.vars = state_vars
)
state_data <- as_tibble(state_data)
state_data$Program <- as.character(state_data$Program)
table(state_data$Program)

state_data <- state_data %>%
    mutate(Score_z = (Score - mean_vec[Program]) / sd_vec[Program])

state_data <- state_data %>%
    mutate(p.val = pnorm(Score_z, lower.tail = F))

state_data <- state_data %>%
    group_by(CellID) %>%
    mutate(p.adj = p.adjust(p.val, "holm"))

state_data$p.sig <- state_data$p.adj < .05

state_stats <- state_data %>%
    group_by(Program) %>%
    filter(p.sig == T) %>%
    summarise(n = n()) %>%
    left_join(
        state_data %>%
            group_by(Program) %>%
            summarise(N = n()),
        by = "Program"
    ) %>%
    mutate(Freq = n / N)

state_data_classify <- state_data %>%
    group_by(CellID) %>%
    filter(p.sig == T) %>%
    arrange(desc(Score)) %>%
    filter(!duplicated(CellID)) %>%
    ungroup()

state_vec <- setNames(rep("Unresolved", nrow(mdata)), mdata$CellID)
state_vec[state_data_classify$CellID] <- state_data_classify$Program
table(state_vec)
table(state_vec) / length(state_vec)

mdata$State <- state_vec[mdata$CellID]
table(mdata$State)

table(mdata$Timepoint, mdata$State)
table(mdata$Timepoint, mdata$State) /
    rowSums(table(mdata$Timepoint, mdata$State))

mdata %>%
    mutate(State = factor(State, c(names(MP_list_named), "Unresolved"))) %>%
    ggplot(aes(x = State, y = nFeature_RNA)) +
    ggdist::stat_halfeye(
        adjust = .5,
        width = .75,
        justification = -.2,
        .width = 0,
        point_colour = NA
    ) +
    geom_boxplot(width = .2, outlier.color = NA) +
    coord_cartesian(xlim = c(1.2, NA)) +
    xlab("") +
    ylab("Complexity") +
    geom_hline(
        yintercept = median(mdata$nFeature_RNA),
        linetype = "dashed",
        size = 1,
        color = "black"
    ) +
    scale_y_continuous(breaks = seq(0, 7000, 1000)) +
    theme_gbm_pvsr(
        axis.text.x = element_text(size = 16, angle = 90),
        panel.border = element_blank()
    ) +
    theme(panel.grid.major = element_line())

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Hybrids classification
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

mp_class <- c(
    "MP_1_RP" = "LQ",
    "MP_2_OPC" = "OPC",
    "MP_4_AC" = "AC",
    "MP_5_Hypoxia" = "Hypoxia",
    "MP_6_MES" = "MES",
    "MP_7_NPC" = "NPC",
    "MP_8_GPC" = "GPC",
    "MP_9_ExN" = "Neuron",
    "MP_10_Stress1" = "Stress",
    "MP_11_MIC" = "Doublet",
    "MP_12_LQ" = "LQ",
    "MP_13_Cilia" = "Cilia",
    "MP_14_NRGN" = "Neuron",
    "MP_15_Stress2" = "Stress"
)
table(mp_class)

state_data <- state_data %>%
    left_join(mdata %>% dplyr::select(CellID, Sample), by = "CellID")

####################################################################################################################################
# First pass - find the two strongest signals
####################################################################################################################################

state_data_hybrid <- state_data %>%
    filter(p.sig == T) %>%
    group_by(CellID) %>%
    summarise(N = n())

state_data_hybrid <- state_data_hybrid %>%
    filter(N > 1)
table(state_data_hybrid$N)

state_data_hybrid$Sample <- sapply(1:nrow(state_data_hybrid), function(i) {
    get_sname(state_data_hybrid$CellID[i], "-")
})

sdata_arranged <- state_data %>%
    filter(p.sig == T) %>%
    group_by(CellID) %>%
    arrange(CellID, desc(Score))

max_prog_class <- sdata_arranged %>%
    group_by(CellID) %>%
    summarise(
        H1 = first(Program),
        H2 = nth(Program, 2),
        S1 = first(Score),
        S2 = nth(Score, 2),
        Sdiff = abs(S1 - S2),
        Class1 = mp_class[H1],
        Class2 = mp_class[H2]
    )
length(which(is.na(max_prog_class$Class1)))

max_prog_class <- max_prog_class %>%
    filter(CellID %in% state_data_hybrid$CellID)
dim(max_prog_class)
length(which(is.na(max_prog_class$Class1)))

table(max_prog_class$Class1)
table(max_prog_class$Class2)
table(max_prog_class$Class1, max_prog_class$Class2)

# If the strongest signal is LQ then try to reclassify using the 2nd strongest signal
# If one of the two strongest signals is doublet then this is the final classification
max_prog_class <- max_prog_class %>%
    ungroup() %>%
    mutate(
        FinalClass = case_when(
            is.na(Class2) ~ Class1,
            Class1 == Class2 ~ Class1,
            Class1 == "LQ" ~ Class2,
            Class1 == "Doublet" ~ Class1,
            Class2 == "Doublet" ~ Class2,
            Class1 != Class2 ~ Class1
        )
    )
length(which(is.na(max_prog_class$H1)))
length(which(is.na(max_prog_class$FinalClass)))

table(max_prog_class$FinalClass)
table(max_prog_class$FinalClass) / nrow(max_prog_class)

class_vec <- setNames(mp_class[mdata$State], mdata$CellID)
class_vec[is.na(class_vec)] <- "Unresolved"
table(class_vec)
length(which(is.na(class_vec)))

class_vec[max_prog_class$CellID] <- max_prog_class$FinalClass
length(which(is.na(class_vec)))
table(class_vec)
table(class_vec) / sum(table(class_vec))
length(class_vec)

####################################################################################################################################
# Remove LQ and residual doublets
####################################################################################################################################

class_vec[class_vec == "LQ" & mdata$isCC] <- "Unresolved"
table(class_vec)
length(class_vec)

class_vec <- class_vec[class_vec %ni% c("LQ", "Doublet")]
length(class_vec)
which(is.na(class_vec))
table(class_vec)

####################################################################################################################################
# Second pass - find the two strongest signals without the doublet and LQ scores
####################################################################################################################################

state_data$Class <- mp_class[state_data$Program]

state_data_hybrid <- state_data %>%
    filter(p.sig == T, Class %ni% c("LQ", "Doublet")) %>%
    group_by(CellID) %>%
    summarise(N = n())
state_data_hybrid <- state_data_hybrid %>% filter(N > 1)
state_data_hybrid <- state_data_hybrid %>% filter(CellID %in% names(class_vec))
dim(state_data_hybrid)
table(state_data_hybrid$N)

state_data_hybrid$Sample <- sapply(1:nrow(state_data_hybrid), function(i) {
    get_sname(state_data_hybrid$CellID[i], "-")
})

sdata_arranged <- state_data %>%
    ungroup() %>%
    filter(
        p.sig == T,
        Class %ni% c("LQ", "Doublet"),
        CellID %in% state_data_hybrid$CellID
    ) %>%
    group_by(CellID) %>%
    arrange(CellID, desc(Score))
dim(sdata_arranged)

max_prog_class <- sdata_arranged %>%
    group_by(CellID) %>%
    summarise(
        H1 = first(Program),
        H2 = nth(Program, 2),
        S1 = first(Score),
        S2 = nth(Score, 2),
        Sdiff = S1 - S2,
        Class1 = mp_class[H1],
        Class2 = mp_class[H2]
    )

max_prog_class <- max_prog_class %>%
    filter(CellID %in% state_data_hybrid$CellID)
dim(max_prog_class)

table(max_prog_class$Class1, max_prog_class$Class2)

permuted_data_merged_sampled <- permuted_data_merged %>%
    select(-MP_3_CC) %>%
    sample_n(10000)

permuted_diffs <- apply(permuted_data_merged_sampled, 1, function(x) {
    x <- sort(x, decreasing = T)
    x[1] - x[2]
})

dm <- rbind(
    max_prog_class %>%
        select(Sdiff) %>%
        mutate(DataType = "Observed"),
    tibble(Sdiff = permuted_diffs, DataType = "Expected")
)

ggplot(
    data = dm,
    aes(
        x = Sdiff,
        y = after_stat(ncount),
        color = DataType,
        linetype = DataType
    )
) +
    geom_freqpoly(
        bins = 100,
        size = 1,
        show.legend = c(color = T, linetype = F)
    ) +
    scale_color_manual(
        name = "",
        values = c("Observed" = "red", "Expected" = "black")
    ) +
    scale_linetype_manual(
        values = c("Observed" = "solid", "Expected" = "dashed")
    ) +
    scale_fill_discrete(name = "Data distribution") +
    geom_vline(
        xintercept = quantile(permuted_diffs, .95),
        linetype = "dashed",
        size = 1,
        color = "black"
    ) +
    xlab("Distribution of score differences") +
    ylab("Count (scaled to 1)") +
    theme_gbm_pvsr(panel.border = element_blank()) +
    theme(panel.grid.major = element_line())

max_prog_class$SdiffWeight <- max_prog_class$Sdiff / abs(max_prog_class$S1)

permuted_diffs_weight <- apply(permuted_data_merged_sampled, 1, function(x) {
    x <- sort(x, decreasing = T)
    (x[1] - x[2]) / abs(x[1])
})

dm <- rbind(
    max_prog_class %>%
        select(SdiffWeight) %>%
        mutate(DataType = "Observed"),
    tibble(SdiffWeight = permuted_diffs_weight, DataType = "Expected")
)

ggplot(
    data = dm,
    aes(
        x = SdiffWeight,
        y = after_stat(ncount),
        color = DataType,
        linetype = DataType
    )
) +
    geom_freqpoly(
        bins = 100,
        size = 1,
        show.legend = c(color = T, linetype = F)
    ) +
    scale_color_manual(
        name = "",
        values = c("Observed" = "red", "Expected" = "black")
    ) +
    scale_linetype_manual(
        values = c("Observed" = "solid", "Expected" = "dashed")
    ) +
    scale_fill_discrete(name = "Data distribution") +
    geom_vline(
        xintercept = quantile(permuted_diffs_weight, .95),
        linetype = "dashed",
        size = 1,
        color = "black"
    ) +
    xlab("Distribution of score differences") +
    ylab("Count (scaled to 1)") +
    theme_gbm_pvsr(panel.border = element_blank()) +
    theme(panel.grid.major = element_line())

max_prog_class <- max_prog_class %>%
    ungroup() %>%
    mutate(
        isHybrid = !is.na(Class2) &
            Class1 != Class2 &
            Sdiff < quantile(permuted_diffs, .95) &
            SdiffWeight < .5
    )
table(max_prog_class$isHybrid)
table(max_prog_class$isHybrid) / nrow(max_prog_class)

max_prog_class %>%
    filter(isHybrid == T) %>%
    dplyr::select(Class1, Class2) %>%
    table()
(max_prog_class %>%
    filter(isHybrid == T) %>%
    dplyr::select(Class1, Class2) %>%
    table()) /
    nrow(max_prog_class %>% filter(isHybrid == T)) *
    100

####################################################################################################################################
# Update the mdata object
####################################################################################################################################

mdata$State_full <- mdata$State
dim(mdata)

mdata <- mdata %>% filter(CellID %in% names(class_vec))
dim(mdata)

mdata$State <- class_vec[mdata$CellID]
table(mdata$State)
length(which(is.na(mdata$State)))

hybrids <- max_prog_class %>% filter(isHybrid == T)
hybrids$edge <- sapply(1:nrow(hybrids), function(i) {
    d <- hybrids[i, c("Class1", "Class2")]
    d <- sort(d)
    paste0(d[, 1], "-", d[, 2])
})

hybrid_vec <- setNames(rep("Singular", nrow(mdata)), mdata$CellID)
hybrid_vec[hybrids$CellID] <- hybrids$edge
table(hybrid_vec)

mdata$HybridID <- hybrid_vec[mdata$CellID]
table(mdata$HybridID) %>% sort()

table(mdata$State, mdata$Timepoint)
table(mdata$Timepoint, mdata$State) /
    rowSums(table(mdata$Timepoint, mdata$State)) *
    100
table(mdata$isCC, mdata$Timepoint)
table(mdata$Timepoint, mdata$isCC) /
    rowSums(table(mdata$Timepoint, mdata$isCC)) *
    100
table(mdata$Timepoint)
