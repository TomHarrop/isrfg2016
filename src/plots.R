library(ggplot2)
library(data.table)
library(DESeq2)

#################
# plot defaults #
#################

set.seed(1)
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
cscale <- c(alpha("#999999", 0.75), "#E41A1C")

#############
# load data #
#############

# lmd tpm
lmd.tpm.wide <- data.table(readRDS("data/lmd/tpm.Rds"), keep.rownames = TRUE)
setnames(lmd.tpm.wide, "rn", "gene")
lmd.exp.wide <- data.table(readRDS("data/lmd/expGenTT.Rds"))
setnames(lmd.exp.wide, "id", "gene")

# lmd lrt
lmd.lrt <- readRDS("data/lmd/ddsLrt.Rds")
lmd.lrt.res <- data.table(data.frame(DESeq2::results(lmd.lrt)),
                      keep.rownames = TRUE, key = "rn")
setnames(lmd.lrt.res, "rn", "gene")
lmd.lrt.p <- unique(lmd.lrt.res[, .(gene, padj)])

# 5acc tpm
fa.tpm <- readRDS("data/fa/tpm_with_calls.Rds")

# 5acc wald results
fa.stage <- readRDS("data/fa/stage_species_results_table.Rds")

# qc
lib.stats <- readRDS("data/lmd/libStats.Rds")

# cluster results
expressionMatrix <- readRDS('data/lmd/expressionMatrix.Rds')
c1 <- readRDS('data/lmd/c1.Rds')

# 20 accessions phenotyping
branching.data <- as.data.table(read.table("data/OsOgObOrPTRAPdata.txt",
                                           sep = "\t",
                                           header = TRUE,
                                           stringsAsFactors = FALSE,
                                           dec = ","))

# crowell gwas genes
crowell.gwas <- readRDS("data/gwas_crowell_genes.Rds")

#########
# munge #
#########

# prepare lmd tpm data
stage.levels <- c(
  "n1" = "RM",
  "n2" = "PBM",
  "n3" = "ePBM & AM",
  "n4" = "SM"
)
lmd.tpm <- melt(lmd.tpm.wide, id.vars = "gene", variable.name = "library",
                value.name = "TPM")
lmd.exp <- melt(lmd.exp.wide, id.vars = "gene", variable.name = "library",
                value.name = "expressed")
setkey(lmd.tpm, gene, library)
setkey(lmd.exp, gene, library)
lmd.tpm.exp <- lmd.exp[lmd.tpm]
lmd.tpm.exp[, stage := substr(library, 1, 2)]
lmd.tpm.exp[, replicate := as.numeric(substr(library, 4, 4))]
lmd.tpm.exp[, stage := factor(plyr::revalue(stage, stage.levels),
                              levels = stage.levels)]
lmd.tpm.exp[, library := NULL]

# prepare 5acc column plot data
fa.to <- c("O. rufipogon", "O. sativa japonica", "O. sativa indica",
           "O. barthii", "O. glaberrima")

fa.tpm <- fa.tpm[, .(
  gene, call, tpm, stage,
  accession = factor(
    plyr::mapvalues(species, from = c("R", "J", "I", "B", "G"), to = fa.to),
    levels = fa.to)
)]
setkey(fa.tpm)
fa.tpm <- unique(fa.tpm)

fa.stage <- fa.stage[, .(
  gene, baseMean, log2FoldChange, lfcSE, padj, 
  accession = factor(
    plyr::mapvalues(accession,
                    from = c("rufipogon", "japonica", "indica",
                             "barthii", "glaberrima"),
                    to = fa.to),
    levels = fa.to
  )
)]
setkey(fa.stage)
fa.stage <- unique(fa.stage)

setkey(fa.stage, gene, accession)
setkey(fa.tpm, gene, accession)
fa.cpdata <- fa.stage[fa.tpm]

############
# QC PLOTS #
############

mv.order <- c("Reads (M)", "Reads in genes (M)", "Detected genes",
              "Dissected area (mm²)", "Yield (ng)", "RIN")
s.order <- structure(c("RM", "PBM", "ePBM & AM", "SM"),
                     .Names = c("RM", "PBM", "ePBM/SBM", "SM"))

qc.pd <- melt(
  lib.stats[RIN >= 7],
  id.vars = c("Sample", "Replicate"),
  measure.vars = mv.order
)
qc.pd[, variable := factor(variable, levels = mv.order)]
qc.pd[, Sample := factor(
  plyr::revalue(Sample, s.order),
  levels = s.order
)]

qc.means <- qc.pd[, .(mean = rutils::GeometricMean(value),
                      lt = "Mean"), by = variable]

qc.plot <- ggplot(qc.pd, aes(
  x = Sample, group = Replicate, y = value, colour = Sample)) +
  theme_slide + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(0.75))
  ) +
  xlab(NULL) + ylab(NULL) +
  geom_hline(data = qc.means,
             aes(yintercept = mean, linetype = lt),
             size = 0.5, colour = "#FF7F00BF") +
  scale_linetype_manual(
    name = NULL,
    values = c(2)
  ) +
  geom_point(position = position_dodge(width = .5)) +
  facet_wrap("variable", scales = "free_y", switch = "y", as.table = FALSE) +
  scale_colour_brewer(palette = "Set1", guide = FALSE)

#############
# G1L PLOTS #
#############

# get G1L genes
g1ls <-c('LOC_Os07g04670', 'LOC_Os02g07030', 'LOC_Os06g46030',
         'LOC_Os02g41460', 'LOC_Os04g43580', 'LOC_Os10g33780',
         'LOC_Os02g56610', 'LOC_Os01g61310', 'LOC_Os05g39500',
         'LOC_Os05g28040')
g1l.order <- c("G1", "G1L1", "G1L2", "G1L3", "G1L4", "G1L5", "G1L6", "G1L9",
               "LOC_Os01g61310", "LOC_Os05g39500")

# make lmd plot
setkey(lmd.tpm.exp, gene)
g1l.lmd <- lmd.tpm.exp[g1ls]
g1l.lmd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
g1l.lmd[is.na(symbol), symbol := gene]
g1l.lmd[, symbol := factor(symbol, levels = g1l.order)]

g1l.lmd.plot <- ggplot(g1l.lmd, aes(x = stage, y = TPM, colour = expressed)) +
  theme_slide + 
  theme(
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab(NULL) +
  scale_color_manual(values = cscale, guide = FALSE) +
  geom_point(size = 2, position = position_jitter(width = 0.25)) +
  facet_wrap("symbol", nrow = 2)

# make svp plot
svps <- c("LOC_Os02g52340", "LOC_Os03g08754", "LOC_Os06g11330")
setkey(lmd.tpm.exp, gene)
svp.lmd <- lmd.tpm.exp[svps]
svp.lmd[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
svp.lmd[is.na(symbol), symbol := gene]
svp.lmd[, symbol := factor(symbol, levels = c(
  "MADS22", "MADS47", "MADS55"
))]

svp.lmd.plot <- ggplot(svp.lmd, aes(x = stage, y = TPM, colour = expressed)) +
  theme_slide + 
  theme(
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab(NULL) +
  scale_color_manual(values = cscale, guide = FALSE) +
  geom_point(size = 2, position = position_jitter(width = 0.25)) +
  facet_wrap("symbol", nrow = 1)


# make fa plot
setkey(fa.cpdata, gene)
g1l.fa <- fa.cpdata[g1ls]
g1l.fa[, symbol := oryzr::LocToGeneName(gene)$symbols, by = gene]
g1l.fa[is.na(symbol), symbol := gene]
g1l.fa[, symbol := factor(symbol, levels = g1l.order)]

my.pal <- c("#fcbba1", "#ef3b2c", "#a50f15",
            "#9ecae1", "#08519c")

my.labels <- plyr::revalue(fa.to, c(
  "O. rufipogon" = bquote(italic("O. rufipogon")),
  "O. sativa japonica" = bquote(italic("O. sativa")~"japonica"),
  "O. sativa indica" = bquote(italic("O. sativa")~"indica"),
  "O. barthii" = bquote(italic("O. barthii")),
  "O. glaberrima" = bquote(italic("O. glaberrima"))))

g1l.fa.plot <- ggplot(g1l.fa, aes(x = stage, y = tpm, colour = accession)) +
  theme_slide + 
  theme(
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.text = element_text(face = "italic"),
    legend.margin = unit(0, "mm")
  ) +
  xlab(NULL) + ylab("TPM") +
  geom_point(size = 2, position = position_jitterdodge(jitter.width = 0.5)) +
  scale_color_manual(values = my.pal,
                     breaks = fa.to,
                     labels = my.labels,
                     guide = guide_legend(title = NULL)) +
  facet_wrap("symbol", nrow = 2)

############
# CLUSTERS #
############

memCutoff = 0.5

# get clusters and membership
cluster <- data.table(id = names(c1$cluster), Cluster = c1$cluster,
                      Membership = apply(c1$membership, 1, max), key = "id")

# get standardised VST counts
exprs <- data.table(Biobase::exprs(expressionMatrix), keep.rownames = TRUE)
exprs[,id := rn][, rn := NULL]
setkey(exprs, "id")

plotData.wide <- exprs[cluster[Membership > memCutoff]]
plotData.wide[, number := length(id), by = Cluster]
plotData.wide[, label := factor(paste0("Cluster ", Cluster,
                                       "\n(", number, " genes)"))]

# relevel clusters
centres.wide <- data.table(c1$centers)

# re-order the cluster for the plot
centres.wide[, Cluster := paste("Cluster", 1:nrow(centres.wide))]

# find the changes between RM and PBM and PBM and SM for each cluster
centres.wide[, c("n1n2", "n2n4") :=
               list(PBM - RM,
                    SM - PBM)]
# divide these changes into categories
centres.wide[, c("dn1n2", "dn2n4") :=
               list(c('dec', 'small', 'inc')[cut(
                 n1n2, breaks = c(-Inf, -0.5, 0.5, Inf))],
                 c('dec', 'small', 'inc')[cut(
                   n2n4, breaks = c(-Inf, -1, 1, Inf))])]               

# first, show gradual increase / decrease
centres.wide[dn1n2 == dn2n4, cOrder := c(1,2)[order(RM)]]

# next, big changes in n1n2 but small in n2n4
centres.wide[dn2n4 == 'small', cOrder := c(3,4)[order(RM)]]

# small changes in n1n2, then large change
centres.wide[dn1n2 == 'small', cOrder := c(5,6)[order(RM)]]

# complex patterns 
centres.wide[!dn1n2 == dn2n4 & !dn1n2 == "small" & !dn2n4 == "small",
             cOrder := c(7,8)[order(SM)]]

# order any leftovers on RM
if (any(is.na(centres.wide[, cOrder]))) {
  orderMax <- max(centres.wide[,cOrder], na.rm = TRUE)
  centres.wide[is.na(cOrder), cOrder := c((orderMax + 1):nrow(centres.wide))]
}

# relevel the clusters by cOrder
plotData.wide[, label := factor(
  label, levels = levels(label)[order(centres.wide[,cOrder])])]

# add label to centres.wide
setkey(centres.wide, 'cOrder')
centres.wide[, label := plotData.wide[, levels(label)]]
centres.wide[, label := factor(label, levels = label)]

# make long
plotData <- reshape2::melt(plotData.wide,
                           id.vars = c("id", "Cluster", "Membership", "label", "number"),
                           variable.name = "Stage",
                           value.name = "Scaled read count")

# fix stage label
plotData[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM & AM")]

# add centres to plot
centres <- reshape2::melt(centres.wide, id.vars = 'label',
                          measure.vars = c("RM", "PBM", "ePBM.SBM", "SM"),
                          variable.name = 'Stage',
                          value.name = "Scaled read count")
centres[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM & AM")]

# main cluster figure
mfuzzClusters <- ggplot(
  plotData, aes(x = Stage, y = `Scaled read count`,
                colour = Membership, group = id)) +
  theme_slide +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab(NULL) +
  scale_colour_gradientn(colours = heatscale, limits = c(0, 1),
                         breaks = seq(0, 1, 0.2), guide = FALSE) +
  geom_line(alpha = 0.8) +
  geom_line(data = centres, mapping = aes(group = 1),
            colour = "#999999") +
  facet_wrap("label", ncol = 4)

##################
# VENN DIAGRAMME #
##################

det.lmd <- lmd.tpm.exp[, .(detected = sum(expressed) > 1),
                       by = . (gene, stage)]

cols <- RColorBrewer::brewer.pal(4, "Set1")[c(1, 4, 2, 3)]
vd <- VennDiagram::venn.diagram(
  list(
    "RM" = det.lmd[detected == TRUE & stage == "RM", unique(gene)],
    "SM" = det.lmd[detected == TRUE & stage == "SM", unique(gene)],
    "PBM" = det.lmd[detected == TRUE & stage == "PBM", unique(gene)],
    "ePBM & AM" = det.lmd[detected == TRUE & stage == "ePBM & AM", unique(gene)]
  ),
  filename = NULL,
  fill = cols,
  lty = "solid",
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  fontfamily = 'Lato',
  cat.fontfamily = 'Lato',
  alpha = 0.5,
  margin = 0.01
)

######################
# PANICLE PHENOTYPES #
######################

# set legend labels
branching.data[, Species := plyr::revalue(
  Origin, c(Ob = "O. barthii", Os = "O. sativa",
            Og = "O. glaberrima", Or = "O. rufipogon"))]

# add number of accessions per species
branching.data[, n_acc := length(unique(Bar_Code)), by = "Species"]
setkey(branching.data, "Species", "Bar_Code", "Sowing_nb", "Repet_nb",
       "Plant_nb", "Panicle_nb")

# select data to plot
branching.pd.wide <- unique(branching.data[, .(
  Species, Bar_Code, Sowing_nb, Repet_nb, Plant_nb, Panicle_nb, Sp_nb,
  SA_nb, TA_nb, n_acc)])

# melt for pa vs sa
branching.pd <- melt(
  branching.pd.wide,
  id.vars = c("Species", "Bar_Code", "Sowing_nb", "Plant_nb", "Panicle_nb",
              "Repet_nb", "n_acc", "Sp_nb"),
  measure.vars = c("SA_nb", "TA_nb"))

# fix labels
var.levels <- c("SA_nb" = "Primary branches", "TA_nb" = "Secondary branches")
branching.pd[, guide_label := paste0(Species, " (", n_acc, ")")]
branching.pd[, variable := factor(
  plyr::revalue(variable, var.levels), levels = var.levels)]
setnames(branching.pd, "Sp_nb", "Spikelets")

# plot
branching <- ggplot(branching.pd,
                    aes(x = value, y = Spikelets, colour = guide_label)) + 
  theme_slide +
  theme(
    strip.background = element_blank(),
    legend.text = element_text(face = "italic"),
    axis.title.y = element_text(size = rel(0.8))
  ) +
  xlab(NULL) +
  scale_color_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", colour = alpha("black", 0.5), fill = cscale[1]) +
  facet_grid( . ~ variable, switch = "x")

################
# GWAS Z-SCORE #
################

# get lmd expression data
lmd.vst.object <- DESeq2::varianceStabilizingTransformation(lmd.lrt, blind = FALSE)
lmd.vst.wide <- data.table(SummarizedExperiment::assay(lmd.vst.object), keep.rownames = TRUE)
setnames(lmd.vst.wide, "rn", "gene")
lmd.vst <- melt(lmd.vst.wide, id.vars = "gene", variable.name = "lib", value.name = "vst")
lmd.vst[, stage := substr(lib, 1, 2)]
lmd.vst[, stage := factor(plyr::revalue(stage, stage.levels),
                              levels = stage.levels)]
lmd.vst.mean <- lmd.vst[, .(vst.mean = rutils::GeometricMean(vst)),
                        by = .(gene, stage)]

# scale lmd values
lmd.mean.wide <- dcast(lmd.vst.mean, stage ~ gene, value.var = "vst.mean")
lmd.mean.mat <- as.matrix(data.frame(lmd.mean.wide, row.names = "stage"))
lmd.vst.scaled.wide <- data.table(t(scale(lmd.mean.mat, center = TRUE)),
                                  keep.rownames = TRUE)
setnames(lmd.vst.scaled.wide, "rn", "gene")
lmd.vst.scaled <- melt(lmd.vst.scaled.wide, id.vars = "gene",
                       variable.name = "stage",
                       value.name = "vst.scaled")

# get expressed genes in qtls
lmd.det.all <- det.lmd[, .(call = any(detected)), by = gene]

# number significant
setkey(lmd.det.all, gene)
lmd.p.call <- lmd.det.all[lmd.lrt.p]
lmd.n.sig <- lmd.p.call[call == TRUE & padj < 0.05, length(unique(gene))]

# merge sig and exp calls with crowell qtls
setkey(crowell.gwas, gene)
crowell.gwas.called <- lmd.p.call[crowell.gwas]
setkey(crowell.gwas.called, gene, Bin_id)
crowell.gwas.called <- unique(crowell.gwas.called)

# number of genes and number expressed per qtl
crowell.hyp <- crowell.gwas.called[call == TRUE, .(
  qtl.k.exp = length(unique(gene, na.rm = TRUE)),
  qtl.q = sum(padj < 0.05, na.rm = TRUE)
), by = Bin_id]

# number of genes expressed in the dataset
n.exp.lmd <- lmd.det.all[, sum(call == TRUE, na.rm = TRUE)]
crowell.hyp[, lmd.n := n.exp.lmd]
crowell.hyp[, lmd.m := lmd.n.sig]

# perform the hypergeometric test
crowell.hyp[, hyper.p := phyper(
  q = qtl.q - 1,
  m = lmd.m,
  n = lmd.n - lmd.m,
  k = qtl.k.exp,
  lower.tail = FALSE
)]

# merge back to gwas data
setkey(crowell.gwas.called,  Bin_id)
setkey(crowell.hyp, Bin_id)
crowell.with.hyper <- crowell.hyp[crowell.gwas.called]

# find an interesting qtl
# chr3: 13143000-13346000 ?
setkey(crowell.with.hyper, hyper.p, Bin_id)
unique(crowell.with.hyper[!is.na(hyper.p) & hyper.p < 0.1])

gwas.total.genes <- crowell.with.hyper[Bin_id == "chr3: 13143000-13346000",
                                       length(unique(gene))]

crowell.with.hyper[Bin_id == "chr3: 13143000-13346000" & padj < 0.05]

# QTL title
gwas.plot.title <- crowell.with.hyper[Bin_id == "chr3: 13143000-13346000", paste0(
  unique(Bin_id), " | ",
  paste(unique(TRAIT), collapse = " "), " | ",
  paste(unique(SUB_POP), collapse = " ")
)]

# extract genes
gwas.plot.genes <- crowell.with.hyper[Bin_id == "chr3: 13143000-13346000",
                   data.table(
                     oryzr::LocToGeneName(unique(gene)), keep.rownames = TRUE)]
setnames(gwas.plot.genes, "rn", "gene")
gwas.plot.genes <- gwas.plot.genes[gene %in% lmd.det.all[call == TRUE, gene]]
setkey(gwas.plot.genes, gene)
setkey(lmd.vst.scaled, gene)

# merge vst values
gwas.pd <- lmd.vst.scaled[gwas.plot.genes]
gwas.pd[is.na(symbols), symbols := gene]
gwas.pd[, stage := factor(stage, levels = rev(s.order))]

# order genes
gwas.pd.wide <- dcast(gwas.pd, stage ~ symbols, value.var = "vst.scaled")
gwas.mat <- as.matrix(data.frame(gwas.pd.wide, row.names = "stage"))
hc <- hclust(dist(t(gwas.mat), method = "minkowski"), method = "ward.D2")
gene.order <- hc$labels[hc$order]
gwas.pd[, symbols := factor(symbols, levels = gene.order)]

# make plot
gwas.lmd.plot <- ggplot(gwas.pd, aes(x = symbols, y = stage, fill = vst.scaled)) +
  theme_slide +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")
  ) +
  ggtitle(gwas.plot.title) +
  xlab(NULL) + ylab(NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_raster() +
  scale_fill_gradientn(colours = heatscale,
                       guide = guide_legend(title = "Scaled\nread count"))


