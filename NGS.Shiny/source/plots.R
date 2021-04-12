# ma plot
# input: a data table (with column mu, log_folchange, q_value, clr and pch), q_value_column, fold_change_column, q_value_thresh, fold_change_thresh, plot title
# output: a MA plot, ggplot object

ma <- function(dt.tb, q_value_thresh, fold_change_thresh, tt){
  
  
  p <- ggplot(dt.tb,
              aes(x = mu,
                  y = log_foldchange,
                  colour = clr,
                  shape = pch)) +
    scale_shape_manual(name = "Legend:",
                       labels = c("No significance",
                                  paste("p-Value < ", q_value_thresh),
                                  paste("p-Value < ", q_value_thresh, " & abs(log2) >= ", fold_change_thresh)),
                       values = c(46, 3, 4)) +
    scale_color_manual(name = "Legend:",
                       values = c("black",
                                  "purple",
                                  "red"),
                       labels = c("No significance",
                                  paste("p-Value < ", q_value_thresh),
                                  paste("p-Value < ", q_value_thresh, " & abs(log2) >= ", fold_change_thresh))) +
    geom_hline(yintercept = c(-fold_change_thresh, fold_change_thresh),
               lty = 2) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.position = "top") +
    geom_point() +
    labs(x = "log2Mean", 
         y = "log2FoldChange",
         title = tt)
  
  return(p)
  
}

# Change-in-gene-expression plot, a two-circle heatmap 
# input: a sorted up.dn_dn.up_table, like the one produced by cige function
# output: a two cricle heatmap of the normalized logfoldcgange
two.cl.heat <- function(up.dn_dn.up_table){
  up.dn_dn.up_table_long <- melt(data = up.dn_dn.up_table,
                                 id.vars = 1,
                                 measure.vars = 2:3,
                                 variable.name = "Comparison",
                                 value.name = "Gene Expression Diff")
  up.dn_dn.up_table_long$Comparison <- factor(up.dn_dn.up_table_long$Comparison,
                                              # inner circle trt3-trt2 and outer circle trt2-trt1
                                              levels = c(colnames(up.dn_dn.up_table)[3],
                                                         colnames(up.dn_dn.up_table)[2]))
  lvls <- up.dn_dn.up_table_long[up.dn_dn.up_table_long$Comparison == colnames(up.dn_dn.up_table)[2], ]
  up.dn_dn.up_table_long$gene <- factor(up.dn_dn.up_table_long$gene,
                                        levels = lvls$gene[order(lvls$`Gene Expression Diff`)])
  
  two_circle_heatmap <- ggplot(data = up.dn_dn.up_table_long) +
    coord_polar("y",
                start = 0,
                direction = -1) +
    geom_tile(aes(x = as.numeric(Comparison),
                  y = gene,
                  fill = `Gene Expression Diff`),
              color = "white") +
    geom_text(data = up.dn_dn.up_table_long[Comparison == colnames(up.dn_dn.up_table)[2], ],
              aes(x = rep(1.75,
                          nlevels(gene)),
                  y = gene,
                  label = unique(gene),
                  angle = 90 + seq(from = 0,
                                   to = 360,
                                   length.out = nlevels(gene))[as.numeric(gene)]),
              hjust = 0) +
    scale_fill_gradient2(low = "red",
                         high = "green",
                         mid = "grey",
                         midpoint = 0,
                         name = "Gene Expr Diff") +
    scale_x_continuous(limits = c(0,
                                  max(as.numeric(up.dn_dn.up_table_long$Comparison)) + 0.5),
                       expand = c(0, 0)) +
    scale_y_discrete("",
                     expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return(two_circle_heatmap)
  
}


# two-column heatmep, sort the input dt decreasing based on diff1, then plot heatmap
# input: a dataframe or tibble, names of gene col, diff1 and diff2 col in the dataframe or tibble, title of the plot, x_axis_titles
# output: a two column heatmap, first column is diff1, second column is diff2
two_column_heatmap <- function(df,
                               gene_col,
                               diff_1,
                               diff_2,
                               title,
                               x_text_1,
                               x_text_2){
  gene_col = enquo(gene_col)
  diff_1 = enquo(diff_1)
  diff_2 = enquo(diff_2)
  
  
  df <- df %>% 
    as.tibble() %>% 
    mutate(gene_reorder = fct_reorder(!!gene_col, !!diff_1)) %>% 
    select(gene_reorder, !!diff_1, !!diff_2) %>% 
    rename(gene = gene_reorder)
  
  df <- df %>% 
    gather(key = "contrast", value = "diff", -gene)
  
 
  
  ggplot(data = df,
       aes(x = contrast, 
           y = gene )) +
  geom_tile(aes(fill = diff), color = "white") +
  scale_fill_gradient2(low = "red", high = "green", mid = "grey", midpoint = 0) +
  labs(title = title) +
  scale_x_discrete(labels = c(x_text_1, x_text_2)) +  
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
        )
}




# starburst plot
# input: a dataframe or tibble, gene col, RNA exp diff col, DNA methyl ratio diff col, 
#        absolute threshhold of RNA diff, absolute threshhold of DNA diff %， string of plot title
starburst <- function(df,
                      gene,
                      rna_exp_diff, 
                      dna_methyl_diff, 
                      region,
                      rna_thresh,
                      dna_thresh,
                      title){
  gene = enquo(gene)
  rna_exp_diff = enquo(rna_exp_diff)
  dna_methyl_diff = enquo(dna_methyl_diff)
  region = enquo(region)
  
  df <- df %>%
    as_tibble() %>%
    mutate(DNA_methyl_diff_pct = !!dna_methyl_diff * 100)
  
  ggplot(data = df,
         aes(text = !!gene,
             x = DNA_methyl_diff_pct,
             y = !!rna_exp_diff,
             fill = !!region)) +
    geom_point(alpha = 0.7, size = 2, shape = 21) +
    geom_hline(yintercept = c(-rna_thresh, rna_thresh), linetype = "dashed") +
    geom_vline(xintercept = c(-dna_thresh, dna_thresh), linetype = "dashed") +
    scale_x_continuous("DNA Methylation Difference (%)") +
    scale_y_continuous("RNA Expression Difference (log2)") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}



