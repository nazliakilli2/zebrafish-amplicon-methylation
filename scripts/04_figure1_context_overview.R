#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = "."), "config.R"))
CpG_FILE <- CPG_FILTERED
CHH_FILE <- CHH_FILTERED
CHG_FILE <- CHG_FILTERED
PLOT_DIR <- FIGURE_DIR
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

# ======================================
# SHARED SETTINGS
# ======================================
GENE_ORDER <- c('hsp70.1','hsp70.2','hsp70.3','actb2','b2m','cyp1a')
FUNC_COLOURS <- c("Stress response"="#1f78b4","Cytoskeletal ref."="#888888",
                  "Immune function"="#e31a1c","Xenobiotic metab."="#33a02c")
CTX_COLOURS  <- c(CpG="#2E86C1", CHG="#1E8449", CHH="#8E44AD")

GROUP_ORDER <- c('Control 7d',
                 'Fipronil Low 7d','Fipronil High 7d',
                 'Cyfluthrin Low 7d','Cyfluthrin High 7d')

ROW_ORDER <- c(
  "CpG - Control 7d",
  "CpG - Fip Low 7d","CpG - Fip High 7d",
  "CpG - Cyt Low 7d","CpG - Cyt High 7d",
  "CHG - Control 7d",
  "CHG - Fip Low 7d","CHG - Fip High 7d",
  "CHG - Cyt Low 7d","CHG - Cyt High 7d",
  "CHH - Control 7d",
  "CHH - Fip Low 7d","CHH - Fip High 7d",
  "CHH - Cyt Low 7d","CHH - Cyt High 7d"
)

theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor=element_blank(),
          strip.background=element_blank(),
          axis.title=element_text(size=9),
          axis.text=element_text(size=8),
          plot.title=element_text(size=10, face="bold"),
          legend.text=element_text(size=8),
          legend.title=element_text(size=9, face="bold"),
          legend.position="bottom")
}

# ======================================
# AMPLICON METADATA
# ======================================
amplicons <- data.table(
  gene      = factor(GENE_ORDER, levels=GENE_ORDER),
  chr       = c("3","3","3","3","4","18"),
  start     = c(26113353,26120181,26126142,40528408,12794094,5594770),
  end       = c(26113669,26120519,26126458,40528826,12794494,5595161),
  size_bp   = c(316,338,316,418,400,391),
  function_ = c("Stress response","Stress response","Stress response",
                "Cytoskeletal ref.","Immune function","Xenobiotic metab."),
  cpg_sites = c(15,38,15,18,28,49),
  chg_sites = c(12,35,12,27,26,47),
  chh_sites = c(33,75,32,40,87,78),
  mean_cov  = c(125,4478,90,886,11024,25632),
  n_samples = c(19,19,19,16,17,18)
)

# ======================================
# PANEL A: Improved amplicon annotation table
# Visual coloured table — no gene body strip
# ======================================
message(">>> Panel A: Improved annotation table...")

ann_tbl <- amplicons[, .(gene, chr, start, end, size_bp, function_,
                          cpg_sites, n_samples)]
ann_tbl[, coord := paste0("chr",chr,":",
                           formatC(start, big.mark=",", format="d"),
                           "–",
                           formatC(end,   big.mark=",", format="d"))]
ann_tbl[, n_row := seq_len(nrow(ann_tbl))]

# Build as a ggplot tile table with coloured function sidebar
pA <- ggplot(ann_tbl, aes(y = reorder(gene, -n_row))) +
  # Function colour sidebar
  geom_tile(aes(x = 0, fill = function_),
            width = 0.18, height = 0.75, alpha = 0.85) +
  # Gene name
  geom_text(aes(x = 0.18, label = gene),
            hjust = 0, size = 3.5, fontface = "bold") +
  # Coordinates
  geom_text(aes(x = 1.35, label = coord),
            hjust = 0.5, size = 2.8, color = "grey30") +
  # Size
  geom_text(aes(x = 2.5, label = paste0(size_bp, " bp")),
            hjust = 0.5, size = 2.8, color = "grey30") +
  # CpG sites
  geom_text(aes(x = 3.1, label = paste0(cpg_sites, " CpGs")),
            hjust = 0.5, size = 2.8, color = "#2166AC", fontface = "bold") +
  # n samples
  geom_text(aes(x = 3.7, label = paste0("n=", n_samples)),
            hjust = 0.5, size = 2.8, color = "grey40") +
  # Function label
  geom_text(aes(x = 4.4, label = function_, color = function_),
            hjust = 0, size = 2.8, fontface = "bold") +
  # Column headers
  annotate("text", x = 0,    y = 6.65, label = "Gene",
           hjust = 0,   size = 3, fontface = "bold", color = "grey20") +
  annotate("text", x = 1.35, y = 6.65, label = "Coordinates (GRCz11)",
           hjust = 0.5, size = 3, fontface = "bold", color = "grey20") +
  annotate("text", x = 2.5,  y = 6.65, label = "Size",
           hjust = 0.5, size = 3, fontface = "bold", color = "grey20") +
  annotate("text", x = 3.1,  y = 6.65, label = "CpG sites",
           hjust = 0.5, size = 3, fontface = "bold", color = "grey20") +
  annotate("text", x = 3.7,  y = 6.65, label = "Samples",
           hjust = 0.5, size = 3, fontface = "bold", color = "grey20") +
  annotate("text", x = 4.4,  y = 6.65, label = "Function",
           hjust = 0,   size = 3, fontface = "bold", color = "grey20") +
  # Divider line under headers
  annotate("segment", x = -0.1, xend = 5.8, y = 6.55, yend = 6.55,
           color = "grey60", linewidth = 0.5) +
  scale_fill_manual(values  = FUNC_COLOURS, guide = "none") +
  scale_color_manual(values = FUNC_COLOURS, guide = "none") +
  scale_x_continuous(limits = c(-0.1, 5.9)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.8))) +
  theme_void() +
  theme(plot.title = element_text(size=11, face="bold",
                                  margin=margin(b=6))) +
  labs(title = "a   Amplicon panel")

# ======================================
# PANEL C: Cytosine site density
# ======================================
message(">>> Panel C...")
site_counts <- melt(amplicons[,.(gene,CpG=cpg_sites,CHG=chg_sites,CHH=chh_sites)],
                   id.vars="gene", variable.name="context", value.name="n_sites")
site_counts[, context := factor(context, levels=c("CpG","CHG","CHH"))]
site_counts[, gene    := factor(gene, levels=GENE_ORDER)]

pC <- ggplot(site_counts, aes(x=gene,y=n_sites,fill=context)) +
  geom_col(position="dodge", width=0.7, alpha=0.85) +
  geom_text(aes(label=n_sites), position=position_dodge(width=0.7),
            vjust=-0.4, size=2.3) +
  scale_fill_manual(values=CTX_COLOURS, name="Context") +
  scale_y_continuous(expand=expansion(mult=c(0,0.18))) +
  theme_pub() +
  theme(panel.grid.major.x=element_blank()) +
  labs(title="C   Cytosine sites per amplicon",
       x="", y="N sites (coverage >= 5x)")

# ======================================
# PANEL D: Coverage
# ======================================
message(">>> Panel D...")
pD <- ggplot(amplicons, aes(x=gene,y=mean_cov,fill=function_)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_text(aes(label=ifelse(mean_cov>=1000,
                              paste0(round(mean_cov/1000,1),"k"),
                              as.character(mean_cov))),
            vjust=-0.4, size=2.5) +
  geom_hline(yintercept=10, linetype="dashed", linewidth=0.5) +
  annotate("text",x=0.6,y=18,label="10x min",size=2.3,color="grey30",hjust=0) +
  scale_fill_manual(values=FUNC_COLOURS, name="Function",
                    guide=guide_legend(nrow=2)) +
  scale_y_log10(labels=scales::comma,
                expand=expansion(mult=c(0,0.15))) +
  theme_pub() +
  theme(panel.grid.major.x=element_blank()) +
  labs(title="D   Mean coverage per amplicon (log scale)",
       x="", y="Mean read depth (log10)")

# ======================================
# PANEL E: Methylation heatmap (pheatmap -> grob)
# ======================================
message(">>> Panel E...")

parse_group <- function(s) {
  if (grepl("-K-",s)) {
    if (grepl("-24-",s))                   return(NA_character_)
    if (grepl("-96-",s)|grepl("96$",s))    return(NA_character_)
    if (grepl("-7-",s) |grepl("7$",s))     return("Control 7d")
    return("Control 7d")
  }
  if (grepl("-170-",s)  &(grepl("-7-",s)|grepl("7$",s)))                          return("Fipronil Low 7d")
  if (grepl("-1700-",s) &(grepl("-7-",s)|grepl("7$",s)))                          return("Fipronil High 7d")
  if (grepl("-18-",s)&!grepl("-180-|-180$",s)&(grepl("-96-",s)|grepl("96$",s)))   return(NA_character_)
  if (grepl("-18-",s)&!grepl("-180-|-180$",s)&(grepl("-7-",s) |grepl("7$",s)))    return("Cyfluthrin Low 7d")
  if (grepl("-180-|-180$",s)&(grepl("-96-",s)|grepl("96$",s)))                    return(NA_character_)
  if (grepl("-180-|-180$",s)&(grepl("-7-",s) |grepl("7$",s)))                     return("Cyfluthrin High 7d")
  return(NA_character_)
}

load_ctx <- function(path, ctx) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)[coverage>=5]
  dt[, meth_percent := 100*meth/coverage]
  meta <- data.table(sample=unique(dt$sample))
  meta[, group := sapply(sample, parse_group)]
  dt <- merge(dt, meta[!is.na(group)], by="sample")
  means <- dt[,.(mean_meth=mean(meth_percent,na.rm=TRUE)),
              by=.(amplicon,group)]
  means[, row_label := paste(ctx, group, sep=" - ")]
  means
}

all_data <- rbindlist(Filter(Negate(is.null),list(
  load_ctx(CpG_FILE,"CpG"),
  load_ctx(CHH_FILE,"CHH"),
  load_ctx(CHG_FILE,"CHG")
)), use.names=TRUE, fill=TRUE)

mat_long <- all_data[row_label %in% ROW_ORDER & amplicon %in% GENE_ORDER]
mat_wide <- dcast(mat_long, row_label ~ amplicon, value.var="mean_meth")
mat_wide <- mat_wide[match(ROW_ORDER[ROW_ORDER %in% mat_wide$row_label],
                           mat_wide$row_label)]
mat <- as.matrix(mat_wide[,-1])
rownames(mat) <- mat_wide$row_label
mat <- mat[, GENE_ORDER]

# Build heatmap as ggplot2 tiles — fully compatible with patchwork
heat_long <- as.data.table(as.data.frame(as.table(mat)))
setnames(heat_long, c("row_label","gene","value"))
heat_long[, value := as.numeric(value)]
heat_long[, gene  := factor(gene, levels=GENE_ORDER)]

# Short row labels
heat_long[, row_short := gsub("Cyfluthrin","Cyt",
                         gsub("Fipronil","Fip",
                         gsub("Control","Ctrl",
                         gsub(" - ","\n", row_label))))]
heat_long[, row_short := factor(row_short, levels=rev(unique(
  gsub("Cyfluthrin","Cyt",
  gsub("Fipronil","Fip",
  gsub("Control","Ctrl",
  gsub(" - ","\n", ROW_ORDER[ROW_ORDER %in% unique(heat_long$row_label)])))))))]

# Context and compound annotation columns
heat_long[, ctx      := sub("\n.*","", row_short)]
heat_long[, compound := ifelse(grepl("Ctrl",row_short),"Control",
                        ifelse(grepl("Fip", row_short),"Fipronil","Cyfluthrin"))]

# Gap indicator rows (white tiles) between context blocks
gap_rows <- c(paste0("CpG\n",c("Ctrl 7d","Fip High 7d","Cyt High 7d")),
              paste0("CHG\n",c("Ctrl 7d","Fip High 7d","Cyt High 7d")))

# Context strip data
ctx_strip <- unique(heat_long[,.(row_short, ctx, compound)])

CTX_COL <- c(CpG="#2E86C1", CHG="#1E8449", CHH="#8E44AD")
CMP_COL <- c(Control="#888888", Fipronil="#1f78b4", Cyfluthrin="#e31a1c")

pE <- ggplot(heat_long, aes(x=gene, y=row_short, fill=value)) +
  geom_tile(color="white", linewidth=0.4) +
  geom_text(aes(label=round(value,0)), size=2.8, color="black") +
  # Context sidebar
  geom_tile(data=ctx_strip,
            aes(x=0.05, y=row_short, fill=NULL),
            width=0.35, height=0.92,
            fill=CTX_COL[ctx_strip$ctx]) +
  # Compound sidebar
  geom_tile(data=ctx_strip,
            aes(x=-0.4, y=row_short, fill=NULL),
            width=0.35, height=0.92,
            fill=CMP_COL[ctx_strip$compound]) +
  scale_fill_gradientn(
    colors=c("#2166AC","#F7F7F7","#D6604D"),
    limits=c(0,100), name="Mean meth (%)") +
  scale_x_discrete(expand=expansion(add=c(0.8,0.5))) +
  theme_pub() +
  theme(axis.text.x=element_text(angle=35,hjust=1,size=7),
        axis.text.y=element_text(size=6),
        panel.border=element_blank(),
        panel.grid=element_blank(),
        legend.position="right") +
  labs(title="b   Mean methylation (%) per group and context",
       x="", y="")



# ======================================
# PANEL F: Co-methylation (7d samples only)
# ======================================
message(">>> Panel F...")

raw_7d <- fread(CpG_FILE)[coverage>=5]
raw_7d[, meth_percent := 100*meth/coverage]
meta7 <- data.table(sample=unique(raw_7d$sample))
meta7[, group := sapply(sample, function(s) {
  g <- parse_group(s)
  if (is.na(g)) return(NA_character_)
  if (grepl("7d",g)) return(g)
  return(NA_character_)
})]
meta7 <- meta7[!is.na(group)]
raw_7d <- merge(raw_7d, meta7, by="sample")

sm7 <- raw_7d[,.(mean_meth=mean(meth_percent,na.rm=TRUE)),
              by=.(sample,amplicon)]
wide7 <- dcast(sm7, sample~amplicon, value.var="mean_meth")
mat7  <- as.matrix(wide7[,-1])
rownames(mat7) <- wide7$sample
mat7 <- mat7[, GENE_ORDER[GENE_ORDER %in% colnames(mat7)], drop=FALSE]

n_g <- ncol(mat7)
genes_present <- colnames(mat7)
cor_mat <- matrix(NA,n_g,n_g,dimnames=list(genes_present,genes_present))
for (i in seq_len(n_g)) for (j in seq_len(n_g)) {
  ok <- complete.cases(mat7[,i],mat7[,j])
  if (sum(ok)>=3) cor_mat[i,j] <- cor(mat7[ok,i],mat7[ok,j])
}

cor_dt <- as.data.table(as.data.frame(as.table(cor_mat)))
setnames(cor_dt,c("gene1","gene2","r"))
cor_dt[, gene1 := factor(gene1, levels=genes_present)]
cor_dt[, gene2 := factor(gene2, levels=rev(genes_present))]

pF <- ggplot(cor_dt, aes(x=gene1,y=gene2,fill=r)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=round(r,2)),size=2.5,
            color=ifelse(abs(cor_dt$r)>0.5,"white","grey30")) +
  scale_fill_gradientn(colors=c("#D9705A","#F7F7F7","#1D6FA4"),
                       limits=c(-1,1), name="Pearson r",
                       breaks=c(-1,-0.5,0,0.5,1)) +
  theme_pub() +
  theme(axis.text.x=element_text(angle=35,hjust=1),
        panel.border=element_blank(),
        panel.grid=element_blank()) +
  labs(title="c   Co-methylation (7d samples)",
       subtitle="Pearson r across 7d samples only",
       x="",y="")

# ======================================
# COMBINE
# ======================================
message(">>> Combining all 6 panels...")

# ======================================
# PANEL E (new): hsp70 cluster schematic
# Shows relative positions of hsp70.1/2/3 and actb2 on chr3
# ======================================
message(">>> Panel E: hsp70 cluster schematic...")

cluster <- data.table(
  gene      = c("hsp70.1","hsp70.2","hsp70.3","actb2"),
  start     = c(26113353, 26120181, 26126142, 40528408),
  end       = c(26113669, 26120519, 26126458, 40528826),
  g_start   = c(26112780, 26119603, 26123736, 40526107),
  g_end     = c(26115559, 26122353, 26127091, 40529804),
  strand    = c(-1, -1, 1, -1),
  function_ = c("Stress response","Stress response",
                "Stress response","Cytoskeletal ref.")
)

# Use genomic coordinates directly on x axis (kb scale)
cluster[, g_start_kb := g_start / 1000]
cluster[, g_end_kb   := g_end   / 1000]
cluster[, amp_s_kb   := start   / 1000]
cluster[, amp_e_kb   := end     / 1000]
cluster[, mid_kb     := (g_start_kb + g_end_kb) / 2]

# Arrow direction by strand
cluster[, arr_x    := ifelse(strand==1, g_end_kb,   g_start_kb)]
cluster[, arr_xend := ifelse(strand==1, g_end_kb+0.05, g_start_kb-0.05)]

pE_new <- ggplot(cluster) +
  # Chromosome backbone
  geom_segment(aes(x=min(g_start_kb)-0.5, xend=max(g_end_kb)+0.5,
                   y=0, yend=0),
               color="grey70", linewidth=3, inherit.aes=FALSE) +
  # Gene body
  geom_rect(aes(xmin=g_start_kb, xmax=g_end_kb,
                ymin=-0.15, ymax=0.15, fill=function_),
            alpha=0.4, color="grey50", linewidth=0.3) +
  # Amplicon
  geom_rect(aes(xmin=amp_s_kb, xmax=amp_e_kb,
                ymin=-0.22, ymax=0.22, fill=function_),
            alpha=0.9, color="grey30", linewidth=0.4) +
  # Strand arrows
  geom_segment(aes(x=arr_x, xend=arr_xend, y=0.28, yend=0.28),
               arrow=arrow(length=unit(0.12,"cm"), type="closed"),
               color="grey40", linewidth=0.5) +
  # Gene labels
  geom_text(aes(x=mid_kb, y=-0.38, label=gene, color=function_),
            size=2.8, fontface="bold") +
  # Distance annotation between hsp70.1 and hsp70.3
  annotate("segment", x=26113.669, xend=26126.142, y=0.55, yend=0.55,
           color="grey50", linewidth=0.4,
           arrow=arrow(ends="both", length=unit(0.1,"cm"), type="open")) +
  annotate("text", x=(26113.669+26126.142)/2, y=0.65,
           label="~12.5 kb", size=2.3, color="grey40") +
  scale_fill_manual(values=FUNC_COLOURS, name="Function") +
  scale_color_manual(values=FUNC_COLOURS, guide="none") +
  scale_x_continuous(labels=function(x) paste0(x," kb")) +
  scale_y_continuous(limits=c(-0.55, 0.8)) +
  theme_pub() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position="bottom") +
  labs(title="E   hsp70 cluster and actb2 on chromosome 3",
       subtitle="Gene bodies (faded) with amplicons (solid). Arrows show transcription direction.",
       x="Genomic position (GRCz11, chr3)")

# ======================================
# PANEL F (new): Sample group overview
# ======================================
message(">>> Panel F: Sample group overview...")


sample_design_7d <- data.table(
  compound  = c("Control","Control",
                "Fipronil","Fipronil",
                "Fipronil","Fipronil",
                "Cyfluthrin","Cyfluthrin",
                "Cyfluthrin","Cyfluthrin"),
  dose      = c("--","--",
                "Low\n170 µg/L","Low\n170 µg/L",
                "High\n1700 µg/L","High\n1700 µg/L",
                "Low\n18 µg/L","Low\n18 µg/L",
                "High\n180 µg/L","High\n180 µg/L"),
  group     = c("Control 7d","Control 7d",
                "Fipronil Low 7d","Fipronil Low 7d",
                "Fipronil High 7d","Fipronil High 7d",
                "Cyfluthrin Low 7d","Cyfluthrin Low 7d",
                "Cyfluthrin High 7d","Cyfluthrin High 7d")
)

GROUP_COLOURS_7D <- c(
  "Control 7d"         = "#444444",
  "Fipronil Low 7d"    = "#a6cee3",
  "Fipronil High 7d"   = "#1f78b4",
  "Cyfluthrin Low 7d"  = "#fdbf6f",
  "Cyfluthrin High 7d" = "#ff7f00"
)

grp_counts <- sample_design_7d[, .(n = .N, dose = dose[1]),
                                  by = .(group, compound)]
grp_counts[, compound := factor(compound,
             levels = c("Control","Fipronil","Cyfluthrin"))]
grp_counts[, group   := factor(group, levels = names(GROUP_COLOURS_7D))]

pF_new <- ggplot(grp_counts, aes(x = "7d", y = compound, fill = group)) +
  geom_tile(color = "white", linewidth = 1.5,
            width = 0.88, height = 0.75) +
  geom_text(aes(label = paste0(dose, "\nn=", n)),
            size = 3.8, lineheight = 1.4,
            color = "white", fontface = "bold") +
  scale_fill_manual(values = GROUP_COLOURS_7D, guide = "none") +
  scale_x_discrete(position = "bottom") +
  theme_pub() +
  theme(panel.grid   = element_blank(),
        panel.border = element_blank(),
        axis.ticks   = element_blank(),
        axis.text.x  = element_text(size = 12, face = "bold"),
        axis.text.y  = element_text(size = 12, face = "bold"),
        plot.title   = element_text(size = 11, face = "bold")) +
  labs(title    = "d   Experimental design (7-day exposures)",
       subtitle = "",
       x = "Exposure duration", y = "")


# ======================================
# COMBINE — 2x2 layout, 4 panels
# A=annotation table, B=heatmap, C=co-methylation, D=design
# ======================================
message(">>> Combining panels...")

pB_final <- pE    # methylation heatmap
pC_final <- pF    # co-methylation
pD_final <- pF_new  # experimental design

layout <- "
AABB
AABB
CCDD
CCDD
"

combined <- pA + pB_final + pC_final + pD_final +
  plot_layout(design=layout) +
  plot_annotation(
    title    = "Figure 1 — Target loci and methylation landscape",
    subtitle = "Six amplicons in Danio rerio GRCz11. Differential methylation analysis in subsequent figures.",
    theme    = theme(plot.title   =element_text(size=13,face="bold"),
                     plot.subtitle=element_text(size=9,color="grey40"))
  )

pdf(file.path(PLOT_DIR,"Figure1_complete.pdf"), width=20, height=14)
print(combined)
dev.off()

# Also save PNG at 350 DPI
png(file.path(PLOT_DIR,"Figure1_complete.png"),
    width=20, height=14, units="in", res=350)
print(combined)
dev.off()

message("Saved: Figure1_complete.pdf and Figure1_complete.png")
message("Done.")
