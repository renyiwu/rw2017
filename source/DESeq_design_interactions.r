# 2
library("DESeq2")
dds <- makeExampleDESeqDataSet(n=100,m=8) # 100 genes and 12 columns. 
dds$genotype <- factor(rep(rep(c("I","II"),each=2),2))
dds$timepoint <- factor(rep(1:2, 4))
dds$group <- factor(paste0(dds$condition, dds$genotype, dds$timepoint))
colData(dds)
# DataFrame with 8 rows and 4 columns
# condition genotype    group timepoint
# <factor> <factor> <factor>  <factor>
# sample1         A        I      AI1         1
# sample2         A        I      AI2         2
# sample3         A       II     AII1         1
# sample4         A       II     AII2         2
# sample5         B        I      BI1         1
# sample6         B        I      BI2         2
# sample7         B       II     BII1         1
# sample8         B       II     BII2         2
#
# Design 0
design(dds) <- ~ condition
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"        "condition_B_vs_A"
results(dds)
# gene1      48.13121    -0.04799584 0.3856385 -0.1244581 0.90095257 0.9584602


# Design 1
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds) 

resultsNames(dds)
# [1] "Intercept"             "genotype_II_vs_I"      "condition_B_vs_A"     
# [4] "genotypeII.conditionB"

results(dds, contrast=c("condition","B","A")) # condition effect For genotype I only
# gene1   157.346759     -0.8257938 0.5081115 -1.6252218  0.1041153 0.7436808

# the condition effect for genotype II
results(dds, list(c("condition_B_vs_A", "genotypeII.conditionB"))) # list is opposite 
# to contrast ie, addition vs substraction.
# gene1   157.346759     -0.2051861 0.5058975  -0.4055883  0.6850452 0.9521637


results(dds, contrast=c("genotype","II","I"))
# gene1   157.346759     0.03299175 0.5056155 0.06525068 0.9479744 0.9934227

results(dds, name = "genotypeII.conditionB")
# gene1   157.346759      0.6206077 0.7170141  0.8655447  0.3867399 0.9074027


# Design 2
design(dds) <- ~ condition + genotype
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"        "condition_B_vs_A" "genotype_II_vs_I"

results(dds, contrast=c("condition","B","A")) # B vs A reagrdless of genotype. calculation on data from both genotypes.
# gene1   157.346759     -0.5140979 0.3615334 -1.4219927 0.1550284 0.7281889
# Note the results are different from in design 0.

results(dds, contrast=c("genotype","II","I")) # II vs I regardless of condition. Calculation on data from both conditions.
# gene1   157.346759      0.3416206 0.3615315  0.9449262 0.3446965   0.96208


# Design 3
design(dds) <- ~ group
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"       "group_AII_vs_AI" "group_BI_vs_AI"  "group_BII_vs_AI"

results(dds, contrast=c("group","AII","AI"))
# gene1   157.346759     0.03299162 0.5056156 0.06525040 0.9479746 0.9934225
results(dds, name = "group_AII_vs_AI")
# gene1   157.346759     0.03299162 0.5056156 0.06525040 0.9479746 0.9934225

results(dds, name = "group_BI_vs_AI")
# gene1   157.346759     -0.8257940 0.5081116 -1.6252218  0.1041153 0.7436808

results(dds, name = "group_BII_vs_AI")
# gene1   157.346759     -0.1721943 0.5057238 -0.3404908 0.7334870 0.9687264

results(dds, contrast = c("group", "BII", "AII"))
# gene1   157.346759     -0.2051859 0.5058976  -0.4055878  0.6850455 0.9521637
# This is the same as in the previous section. note the padj value is also the same.

# Design 4
design(dds) <- ~ condition*genotype # equals to "condition + genotype + condition:genotype"
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"             "condition_B_vs_A"      "genotype_II_vs_I"     
# [4] "conditionB.genotypeII"

### 

# Design 5
design(dds) <- ~ condition + condition:genotype # 
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"             "condition_B_vs_A"      "conditionA.genotypeII"
# [4] "conditionB.genotypeII"

results(dds, name = "condition_B_vs_A")[1,]
# gene1  17.36441      0.6167375 0.8695201 0.7092849 0.4781477 0.9696017

results(dds, contrast = c("condition", "B", "A"))[1,]
# gene1  17.36441      0.6167375 0.8695201 0.7092849 0.4781477 0.9696017
# the same as above

results(dds, name = "conditionA.genotypeII")[1,]
# gene1  17.36441    -0.04546674 0.8881536 -0.05119243 0.9591722 0.9888516

results(dds, name = "conditionB.genotypeII")[1,]
# gene1  17.36441      0.3708084 0.8405654 0.4411416 0.6591105 0.9791759

# Design 6
design(dds) <- ~ condition + condition:genotype + condition:timepoint# 
colData(dds)
dds <- DESeq(dds) 
resultsNames(dds)
# [1] "Intercept"             "condition_B_vs_A"      "conditionA.genotypeII"
# [4] "conditionB.genotypeII" "conditionA.timepoint2" "conditionB.timepoint2"

results(dds, name = "condition_B_vs_A")[1,]
# gene1  17.36441      0.8609872  1.248727 0.6894918 0.4905139 0.9907193

results(dds, name = "conditionA.genotypeII")[1,]
# gene1  17.36441    0.003773875  1.036011 0.003642699 0.9970936 0.9970936

results(dds, name = "conditionB.genotypeII")[1,]
# gene1  17.36441      0.3586995   0.99382 0.3609301 0.7181517 0.9997351

results(dds, name = "conditionA.timepoint2")[1,]
# gene1  17.36441      0.4819115  1.036205 0.4650737 0.6418787 0.9390541

results(dds, name = "conditionB.timepoint2")[1,]
# gene1  17.36441      0.0891689 0.9937896 0.08972613 0.9285048 0.9971533

results(dds, tidy = T, contrast = list("conditionB.timepoint2","conditionA.timepoint2" ))[1,]
# gene1  17.36441     -0.3927426  1.435736 -0.273548  0.784432 0.9884171
