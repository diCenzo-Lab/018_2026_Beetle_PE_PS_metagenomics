PE_PS_metagenomics_statistics
================
2026-01-15

# Statistical analyses

This document incudes the code for running all the statistical analyses
of the mealworm and superworm gut microbiome data reported in this
study. It includes the alpha-diversity, beta-diversity, and differential
abundance analyses.

# Load required libraries

``` r
library(lme4)
library(phyloseq)
library(car)
library(emmeans)
library(ggplot2)
library(patchwork)
library(vegan)
library(dplyr)
library(tidyr)
library(kableExtra)
library(ANCOMBC)
```

# Analyses between insect species

## Prepare a phyloseq object for the full dataset across both insects

Use the MAG abundance data and the GTDB-Tk output data to generate a
phyloseq object.

``` r
# Load the OTU table
otumat <- as.matrix(read.table("combined_MAG_abundances_depth_absolute.txt",
                        header = TRUE, row.names = 1))

# Load the taxonomy table
taxmat <- as.matrix(read.table("combined_MAG_taxonomy.txt",
                        header = TRUE, row.names = 1))

# Load the sample metadata
metadata <- read.table("combined_sample_metadata.txt",
                       header = TRUE, row.names = 1)

# Combine the data into a phyloseq object
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, META)

# Convert MAG abundances to percentages
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
```

## Calculate alpha-diversity and perform statistical analysis

Calculate Shannon and Simpson’s diversity metrics and test if the values
are statistically different between insect species. Insect batcth is
included as a random effect, but replicate was excluded as it had 0
variance and led to non-convergence of the model.

``` r
# Calculate alpha-diversity metrics
richness <- estimate_richness(ps.rel, measures = c("Shannon", "Simpson"))

# Add columns with relevant metadata
richness$insect <- metadata$Insect
richness$replicate <- metadata$Replicate
richness$batch <- metadata$Insect_batch
richness$depth <- metadata$Read_count
richness$diet <- metadata$Treatment

# Create factors for treatment, ordering the treatments as desired
richness$treatment <- factor(richness$insect, 
                                 levels = c("Mealworm", "Superworm"))

# Perform statistical analysis on the Shannon data
model.shannon <- lmer(log(Shannon) ~ insect + depth + insect:depth + (1 | batch),
                      data=richness)
summary(model.shannon)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(Shannon) ~ insect + depth + insect:depth + (1 | batch)
    ##    Data: richness
    ## 
    ## REML criterion at convergence: 10
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8875 -0.5979  0.2070  0.5490  2.1578 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  batch    (Intercept) 0.03011  0.1735  
    ##  Residual             0.03636  0.1907  
    ## Number of obs: 46, groups:  batch, 6
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           0.332097   0.210999   1.574
    ## insectSuperworm       0.228621   0.281683   0.812
    ## depth                 0.007972   0.003582   2.225
    ## insectSuperworm:depth 0.002289   0.005199   0.440
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) insctS depth 
    ## insctSprwrm -0.749              
    ## depth       -0.859  0.643       
    ## insctSprwr:  0.592 -0.836 -0.689

``` r
model.shannon.aov <- Anova(model.shannon, type = 3)
model.shannon.aov
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(Shannon)
    ##               Chisq Df Pr(>Chisq)  
    ## (Intercept)  2.4772  1    0.11551  
    ## insect       0.6587  1    0.41701  
    ## depth        4.9527  1    0.02605 *
    ## insect:depth 0.1938  1    0.65979  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
model.shannon.bt <- update(ref_grid(model.shannon), tran = "log")
model.shannon.em <- emmeans(model.shannon.bt, pairwise ~ insect,
                            infer = TRUE, type = "response")
model.shannon.em$emmeans
```

    ##  insect    response    SE   df lower.CL upper.CL null t.ratio p.value
    ##  Mealworm      2.00 0.220 4.18     1.48     2.70    1   6.296  0.0028
    ##  Superworm     2.79 0.303 4.05     2.06     3.76    1   9.401  0.0007
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
model.shannon.em$contrasts
```

    ##  contrast             ratio    SE   df lower.CL upper.CL null t.ratio p.value
    ##  Mealworm / Superworm 0.717 0.111 4.12    0.469      1.1    1  -2.145  0.0966
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
# Perform statistical analysis on the Simpson data
model.simpson <- lmer(log(Simpson) ~ insect + depth + insect:depth + (1 | batch),
                      data=richness)
summary(model.simpson)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(Simpson) ~ insect + depth + insect:depth + (1 | batch)
    ##    Data: richness
    ## 
    ## REML criterion at convergence: -14.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3506 -0.4401  0.2107  0.4997  2.4932 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  batch    (Intercept) 0.02032  0.1426  
    ##  Residual             0.01978  0.1406  
    ## Number of obs: 46, groups:  batch, 6
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)           -0.580869   0.160168  -3.627
    ## insectSuperworm        0.086153   0.214348   0.402
    ## depth                  0.005750   0.002651   2.169
    ## insectSuperworm:depth  0.001498   0.003843   0.390
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) insctS depth 
    ## insctSprwrm -0.747              
    ## depth       -0.837  0.626       
    ## insctSprwr:  0.577 -0.812 -0.690

``` r
model.simpson.aov <- Anova(model.simpson, type = 3)
model.simpson.aov
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(Simpson)
    ##                Chisq Df Pr(>Chisq)    
    ## (Intercept)  13.1524  1  0.0002871 ***
    ## insect        0.1615  1  0.6877372    
    ## depth         4.7054  1  0.0300688 *  
    ## insect:depth  0.1520  1  0.6966738    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
model.simpson.bt <- update(ref_grid(model.simpson), tran = "log")
model.simpson.em <- emmeans(model.simpson.bt, pairwise ~ insect,
                            infer = TRUE, type = "response")
model.simpson.em$emmeans
```

    ##  insect    response     SE   df lower.CL upper.CL null t.ratio p.value
    ##  Mealworm     0.725 0.0644 4.15    0.569    0.925    1  -3.614  0.0211
    ##  Superworm    0.846 0.0746 4.05    0.663    1.079    1  -1.897  0.1299
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
model.simpson.em$contrasts
```

    ##  contrast             ratio    SE  df lower.CL upper.CL null t.ratio p.value
    ##  Mealworm / Superworm 0.857 0.107 4.1    0.608     1.21    1  -1.229  0.2848
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

## Plot the alpha-diversity metrics

Generate plots showing the Shannon and Simpson’s diversity metrics

``` r
# Collect the emmeans for Shannon data
emms.shannon <- as.data.frame(model.shannon.em$emmeans)
emms.shannon$metric <- 'Shannon'

# Collect the emmeans for teh Simpson's data
emms.simpson <- as.data.frame(model.simpson.em$emmeans)
emms.simpson$metric <- 'Simpson'

# Plot the Shannon data
p.shannon <- ggplot() +
  geom_bar(data = emms.shannon,
           aes(x = insect, y = response, fill = insect),
           stat = "identity",
           colour = "black") +
  geom_errorbar(data = emms.shannon,
                aes(x = insect, ymin = response-SE, ymax = response+SE),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_jitter(data = richness,
              aes(x = insect, y = Shannon),
              width = 0.15, size = 1.25) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  theme(text = element_text(size = 17), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  scale_y_continuous(lim=c(0,3.5), breaks=c(0.0,0.7,1.4,2.1,2.8,3.5)) + 
  xlab("") + ylab("")

# Plot the Simpson's data
p.simpson <- ggplot() +
  geom_bar(data = emms.simpson,
           aes(x = insect, y = response, fill = insect),
           stat = "identity",
           colour = "black") +
  geom_errorbar(data = emms.simpson,
                aes(x = insect, ymin = response-SE, ymax = response+SE),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_jitter(data = richness,
              aes(x = insect, y = Simpson),
              width = 0.15, size = 1.25) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  theme(text = element_text(size = 17), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black")) +
  scale_y_continuous(lim=c(0,1.0), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0)) + 
  xlab("") + ylab("")

# Combine the plots as two panels in one figure
p.shannon + p.simpson
```

![](Statistics_files/figure-gfm/Insect%20alpha%20diversity%20plots-1.png)<!-- -->

``` r
# Export a PDF of the figure
pdf('Insect.alpha.pdf', width=10, height=5)
p.shannon + p.simpson
dev.off()
```

## Calculate beta-diversity and make a PCoA plot at the species level to compare insect species

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.rel.sp <- tax_glom(ps.rel, "Species", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.sp, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.sp))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Insect, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm  Pr(>F)  
    ## Groups     1 0.03284 0.032844 3.6265  10000 0.06619 .
    ## Residuals 44 0.39849 0.009057                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           Mealworm Superworm
    ## Mealworm              0.0651
    ## Superworm 0.063417

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Insect + Read_count + Insect:Read_count, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
```

    ##                   Df   SumOfSqs         R2         F     Pr(>F)
    ## Insect             1  5.8557853 0.41851350 33.720298 0.00009999
    ## Read_count         1  0.4995291 0.03570139  2.876518 0.02229777
    ## Insect:Read_count  1  0.3429345 0.02450956  1.974774 0.07369263
    ## Residual          42  7.2936184 0.52127556        NA         NA
    ## Total             45 13.9918672 1.00000000        NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.sp, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.sp, pso, type="samples", color="Insect", axes = c(1,2))
p.beta.sp <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 4) +
  geom_point(size = 4, shape=21, stroke=0.5, color="black") +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x))
```

## Calculate beta-diversity and make a PCoA plot at the genus level to compare insect species

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.rel.genus <- tax_glom(ps.rel, "Genus", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.genus, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.genus))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Insect, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
    ## Groups     1 0.04750 0.047504 4.6425  10000 0.0348 *
    ## Residuals 44 0.45023 0.010233                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           Mealworm Superworm
    ## Mealworm              0.0348
    ## Superworm 0.036704

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Insect + Read_count + Insect:Read_count, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
```

    ##                   Df  SumOfSqs         R2         F     Pr(>F)
    ## Insect             1 2.5258023 0.31422831 23.336719 0.00009999
    ## Read_count         1 0.8437962 0.10497443  7.796111 0.00049995
    ## Insect:Read_count  1 0.1227287 0.01526835  1.133931 0.31616838
    ## Residual          42 4.5457845 0.56552890        NA         NA
    ## Total             45 8.0381117 1.00000000        NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.genus, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.genus, pso, type="samples", color="Insect", axes = c(1,2))
p.beta.genus <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 4) +
  geom_point(size = 4, shape=21, stroke=0.5, color="black") +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.60,0.60), breaks=c(-0.60,-0.30,0.00,0.30,0.60),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.60,0.60), breaks=c(-0.60,-0.30,0.00,0.30,0.60),
                     labels = function(x) sprintf("%.2f", x))
```

## Calculate beta-diversity and make a PCoA plot at the family level to compare insect species

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.rel.fam <- tax_glom(ps.rel, "Family", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.fam, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.fam))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Insect, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.02850 0.028505 2.1438  10000 0.1509
    ## Residuals 44 0.58504 0.013296                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           Mealworm Superworm
    ## Mealworm              0.1517
    ## Superworm  0.15026

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Insect + Read_count + Insect:Read_count, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
```

    ##                   Df  SumOfSqs          R2          F     Pr(>F)
    ## Insect             1 0.8331004 0.162824877 10.5936034 0.00009999
    ## Read_count         1 0.9703743 0.189654314 12.3391623 0.00009999
    ## Insect:Read_count  1 0.0101107 0.001976081  0.1285665 0.95450455
    ## Residual          42 3.3029569 0.645544729         NA         NA
    ## Total             45 5.1165423 1.000000000         NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.fam, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.fam, pso, type="samples", color="Insect", axes = c(1,2))
p.beta.fam <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 4) +
  geom_point(size = 4, shape=21, stroke=0.5, color="black") +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.40,0.40), breaks=c(-0.40,-0.20,0.00,0.20,0.40), 
                     labels = function(x) sprintf("%.2f", x))
```

## Calculate beta-diversity and make a PCoA plot at the order level to compare insect species

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.rel.order <- tax_glom(ps.rel, "Order", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.order, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.order))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Insect, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00777 0.0077727 0.5735  10000 0.4524
    ## Residuals 44 0.59638 0.0135540                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           Mealworm Superworm
    ## Mealworm              0.4479
    ## Superworm  0.45292

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Insect + Read_count + Insect:Read_count, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
```

    ##                   Df   SumOfSqs          R2          F     Pr(>F)
    ## Insect             1 0.68120240 0.153020087  9.8142709 0.00019998
    ## Read_count         1 0.84305863 0.189378230 12.1461782 0.00019998
    ## Insect:Read_count  1 0.01226405 0.002754903  0.1766916 0.90740926
    ## Residual          42 2.91519374 0.654846780         NA         NA
    ## Total             45 4.45171883 1.000000000         NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.order, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.order, pso, type="samples", color="Insect", axes = c(1,2))
p.beta.order <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 4) +
  geom_point(size = 4, shape=21, stroke=0.5, color="black") +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.40,0.40), breaks=c(-0.40,-0.20,0.00,0.20,0.40), 
                     labels = function(x) sprintf("%.2f", x))
```

## Calculate beta-diversity and make a PCoA plot at the class level to compare insect species

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.rel.class <- tax_glom(ps.rel, "Class", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.class, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.class))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Insect, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00178 0.0017799 0.1108  10000  0.742
    ## Residuals 44 0.70688 0.0160654                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           Mealworm Superworm
    ## Mealworm              0.7419
    ## Superworm  0.74083

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Insect + Read_count + Insect:Read_count, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
```

    ##                   Df    SumOfSqs           R2         F     Pr(>F)
    ## Insect             1  0.21276008  0.085056140  5.098741 0.02319768
    ## Read_count         1  0.54904315  0.219493675 13.157680 0.00039996
    ## Insect:Read_count  1 -0.01296997 -0.005185068 -0.310822 1.00000000
    ## Residual          42  1.75257436  0.700635253        NA         NA
    ## Total             45  2.50140761  1.000000000        NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.class, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.class, pso, type="samples", color="Insect", axes = c(1,2))
p.beta.class <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 4) +
  geom_point(size = 4, shape=21, stroke=0.5, color="black") +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.60,0.60), breaks=c(-0.60,-0.30,0.00,0.30,0.60),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.20,0.20), breaks=c(-0.20,-0.10,0.00,0.10,0.20), 
                     labels = function(x) sprintf("%.2f", x))
```

## Calculate beta-diversity and make a PCoA plot at the phylum level to compare insect species

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.rel.phylum <- tax_glom(ps.rel, "Phylum", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.phylum, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.phylum))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Insect, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     1 0.00605 0.006055 0.3651  10000 0.5386
    ## Residuals 44 0.72968 0.016584                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           Mealworm Superworm
    ## Mealworm              0.5409
    ## Superworm  0.54878

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Insect + Read_count + Insect:Read_count, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
```

    ##                   Df     SumOfSqs           R2          F     Pr(>F)
    ## Insect             1  0.135136801  0.057386547  3.3612645 0.06649335
    ## Read_count         1  0.539493065  0.229098542 13.4188382 0.00039996
    ## Insect:Read_count  1 -0.008352841 -0.003547077 -0.2077606 0.99990001
    ## Residual          42  1.688574560  0.717061988         NA         NA
    ## Total             45  2.354851585  1.000000000         NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.phylum, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.phylum, pso, type="samples", color="Insect", axes = c(1,2))
p.beta.phylum <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 4) +
  geom_point(size = 4, shape=21, stroke=0.5, color="black") +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.60,0.60), breaks=c(-0.60,-0.30,0.00,0.30,0.60),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.20,0.20), breaks=c(-0.20,-0.10,0.00,0.10,0.20), 
                     labels = function(x) sprintf("%.2f", x))
```

## Make a figure of all the PCoA plots at different taxonomic levels

``` r
p.beta.sp + p.beta.genus + p.beta.fam + p.beta.order + p.beta.class + p.beta.phylum
```

![](Statistics_files/figure-gfm/Insect%20beta%20diversity%20plot-1.png)<!-- -->

``` r
pdf('Insect.beta.pdf', width=9, height=5)
p.beta.sp + p.beta.genus + p.beta.fam + p.beta.order + p.beta.class + p.beta.phylum
dev.off()
```

## Compare taxon abundance across insects

Summarize the most abundant taxa at the phylum level in the mealworms
and superworms across diets.

``` r
# Agglomerate data at phylum level 
ps.rel.phylum <- tax_glom(ps.rel, "Phylum", NArm = FALSE)

# Create a dataframe from phyloseq object
ps.rel.phylum.df <- psmelt(ps.rel.phylum)

# Get summary statistics
df <- ps.rel.phylum.df %>%
  group_by(Insect, Phylum) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            stdev = sd(Abundance, na.rm = TRUE),
            counts = n(),
            se = stdev / sqrt(counts))

# Create a table rounded to 1 decimal place
df.w <- df %>%
  mutate(stat = paste0(round(100 * mean, 2), " +/- ", round(100 * stdev, 2))) %>%
  select(Insect, Phylum, stat) %>%
  pivot_wider(
    names_from = Insect,
    values_from = stat
  )
kable(df.w, align = "l",
      caption = "Abundance of phyla across insect species") %>%
  kable_classic(full_width = FALSE)
```

<table class=" lightable-classic" style="color: black; font-family: &quot;Arial Narrow&quot;, &quot;Source Sans Pro&quot;, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Abundance of phyla across insect species
</caption>

<thead>

<tr>

<th style="text-align:left;">

Phylum
</th>

<th style="text-align:left;">

Mealworm
</th>

<th style="text-align:left;">

Superworm
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Actinomycetota
</td>

<td style="text-align:left;">

1.86 +/- 3.4
</td>

<td style="text-align:left;">

0.95 +/- 0.79
</td>

</tr>

<tr>

<td style="text-align:left;">

Bacillota
</td>

<td style="text-align:left;">

63.46 +/- 23.89
</td>

<td style="text-align:left;">

65.23 +/- 20.26
</td>

</tr>

<tr>

<td style="text-align:left;">

Bacteroidota
</td>

<td style="text-align:left;">

0.11 +/- 0.53
</td>

<td style="text-align:left;">

5.32 +/- 2.85
</td>

</tr>

<tr>

<td style="text-align:left;">

Fusobacteriota
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.29 +/- 0.25
</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudomonadota
</td>

<td style="text-align:left;">

34.57 +/- 22.61
</td>

<td style="text-align:left;">

28.2 +/- 18.14
</td>

</tr>

</tbody>

</table>

Summarize the taxa at the genus level in the mealworms and superworms
for the 40 most abundant genera

``` r
# Agglomerate data at genus level 
ps.rel.genus <- tax_glom(ps.rel, "Genus", NArm = FALSE)

# Identify the most abundant orders
top.taxa = names(sort(taxa_sums(ps.rel.genus), TRUE)[1:40])
top.taxa.table = cbind(tax_table(ps.rel.genus), topTaxa = NA)
top.taxa.table[top.taxa, "topTaxa"] <- as(tax_table(ps.rel.genus)[top.taxa, "Genus"], 
                                          "character")
tax_table(ps.rel.genus) <- tax_table(top.taxa.table)

# Create a dataframe from phyloseq object
ps.rel.genus.df <- psmelt(ps.rel.genus)
ps.rel.genus.df$Genus <- as.character(ps.rel.genus.df$Genus) #convert to character
ps.rel.genus.df.noNA <- na.omit(ps.rel.genus.df)
ps.rel.genus.df <- ps.rel.genus.df.noNA

# Get summary statistics
df <- ps.rel.genus.df %>%
  group_by(Insect, Genus) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            stdev = sd(Abundance, na.rm = TRUE),
            counts = n(),
            se = stdev / sqrt(counts))

# Create a table rounded to 1 decimal place
df.w <- df %>%
  mutate(stat = paste0(round(100 * mean, 2), " +/- ", round(100 * stdev, 2))) %>%
  select(Insect, Genus, stat) %>%
  pivot_wider(
    names_from = Insect,
    values_from = stat
  )
kable(df.w, align = "l",
      caption = "Abundance of genera across insect species") %>%
  kable_classic(full_width = FALSE)
```

<table class=" lightable-classic" style="color: black; font-family: &quot;Arial Narrow&quot;, &quot;Source Sans Pro&quot;, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Abundance of genera across insect species
</caption>

<thead>

<tr>

<th style="text-align:left;">

Genus
</th>

<th style="text-align:left;">

Mealworm
</th>

<th style="text-align:left;">

Superworm
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Apibacter
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.64 +/- 0.52
</td>

</tr>

<tr>

<td style="text-align:left;">

CALYQQ01
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

1.88 +/- 1.12
</td>

</tr>

<tr>

<td style="text-align:left;">

Chryseobacterium
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.58 +/- 0.67
</td>

</tr>

<tr>

<td style="text-align:left;">

Citrobacter
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

2.83 +/- 2.34
</td>

</tr>

<tr>

<td style="text-align:left;">

Citrobacter_B
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

2.16 +/- 3.54
</td>

</tr>

<tr>

<td style="text-align:left;">

Corynebacterium
</td>

<td style="text-align:left;">

1.7 +/- 3.07
</td>

<td style="text-align:left;">

0.95 +/- 0.79
</td>

</tr>

<tr>

<td style="text-align:left;">

DAIDKD01
</td>

<td style="text-align:left;">

0.36 +/- 1.31
</td>

<td style="text-align:left;">

0 +/- 0
</td>

</tr>

<tr>

<td style="text-align:left;">

Dysgonomonas
</td>

<td style="text-align:left;">

0.11 +/- 0.53
</td>

<td style="text-align:left;">

3.68 +/- 2.89
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterobacter
</td>

<td style="text-align:left;">

5.84 +/- 6.51
</td>

<td style="text-align:left;">

0.79 +/- 0.71
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterococcus
</td>

<td style="text-align:left;">

2.83 +/- 2
</td>

<td style="text-align:left;">

0.85 +/- 0.58
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterococcus_A
</td>

<td style="text-align:left;">

0.37 +/- 0.45
</td>

<td style="text-align:left;">

3.56 +/- 2.1
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterococcus_B
</td>

<td style="text-align:left;">

2.41 +/- 2.02
</td>

<td style="text-align:left;">

0.44 +/- 0.59
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterococcus_C
</td>

<td style="text-align:left;">

0.02 +/- 0.08
</td>

<td style="text-align:left;">

2.8 +/- 2.16
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterococcus_D
</td>

<td style="text-align:left;">

1.63 +/- 1.53
</td>

<td style="text-align:left;">

2.4 +/- 2.13
</td>

</tr>

<tr>

<td style="text-align:left;">

Enterococcus_J
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.54 +/- 0.56
</td>

</tr>

<tr>

<td style="text-align:left;">

Entomomonas
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

3.62 +/- 4.95
</td>

</tr>

<tr>

<td style="text-align:left;">

Hafnia
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.24 +/- 0.86
</td>

</tr>

<tr>

<td style="text-align:left;">

Intestinirhabdus
</td>

<td style="text-align:left;">

7.25 +/- 7.34
</td>

<td style="text-align:left;">

1.79 +/- 2.65
</td>

</tr>

<tr>

<td style="text-align:left;">

JAGNPU01
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

1 +/- 0.64
</td>

</tr>

<tr>

<td style="text-align:left;">

JALAGN01
</td>

<td style="text-align:left;">

30.18 +/- 25.4
</td>

<td style="text-align:left;">

42.17 +/- 21.28
</td>

</tr>

<tr>

<td style="text-align:left;">

Jejubacter
</td>

<td style="text-align:left;">

8.92 +/- 9.56
</td>

<td style="text-align:left;">

0 +/- 0
</td>

</tr>

<tr>

<td style="text-align:left;">

Klebsiella
</td>

<td style="text-align:left;">

0.83 +/- 0.95
</td>

<td style="text-align:left;">

5.69 +/- 6.32
</td>

</tr>

<tr>

<td style="text-align:left;">

Kluyvera
</td>

<td style="text-align:left;">

0.92 +/- 1.85
</td>

<td style="text-align:left;">

5.31 +/- 4.41
</td>

</tr>

<tr>

<td style="text-align:left;">

Lactococcus
</td>

<td style="text-align:left;">

13.73 +/- 8.8
</td>

<td style="text-align:left;">

7.75 +/- 10.54
</td>

</tr>

<tr>

<td style="text-align:left;">

Latilactobacillus
</td>

<td style="text-align:left;">

6.56 +/- 7.35
</td>

<td style="text-align:left;">

0 +/- 0
</td>

</tr>

<tr>

<td style="text-align:left;">

Listeria
</td>

<td style="text-align:left;">

0.59 +/- 1.41
</td>

<td style="text-align:left;">

0 +/- 0
</td>

</tr>

<tr>

<td style="text-align:left;">

Merdenecus
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.65 +/- 0.44
</td>

</tr>

<tr>

<td style="text-align:left;">

Mixta
</td>

<td style="text-align:left;">

2.72 +/- 3.17
</td>

<td style="text-align:left;">

0.4 +/- 0.35
</td>

</tr>

<tr>

<td style="text-align:left;">

Morganella
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

1.17 +/- 3.16
</td>

</tr>

<tr>

<td style="text-align:left;">

Pediococcus
</td>

<td style="text-align:left;">

3.3 +/- 4.94
</td>

<td style="text-align:left;">

1.3 +/- 2.47
</td>

</tr>

<tr>

<td style="text-align:left;">

Pseudomonas
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.3 +/- 0.57
</td>

</tr>

<tr>

<td style="text-align:left;">

Romboutsia_A
</td>

<td style="text-align:left;">

0.01 +/- 0.05
</td>

<td style="text-align:left;">

1 +/- 1.09
</td>

</tr>

<tr>

<td style="text-align:left;">

Scandinavium
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.35 +/- 0.41
</td>

</tr>

<tr>

<td style="text-align:left;">

Sebaldella
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.29 +/- 0.25
</td>

</tr>

<tr>

<td style="text-align:left;">

Sphingobacterium
</td>

<td style="text-align:left;">

0 +/- 0
</td>

<td style="text-align:left;">

0.27 +/- 0.31
</td>

</tr>

<tr>

<td style="text-align:left;">

Staphylococcus
</td>

<td style="text-align:left;">

1 +/- 2.89
</td>

<td style="text-align:left;">

0.2 +/- 0.84
</td>

</tr>

<tr>

<td style="text-align:left;">

Tenebrionibacter
</td>

<td style="text-align:left;">

8.1 +/- 10.61
</td>

<td style="text-align:left;">

0 +/- 0
</td>

</tr>

<tr>

<td style="text-align:left;">

Weissella
</td>

<td style="text-align:left;">

0.28 +/- 0.76
</td>

<td style="text-align:left;">

0.06 +/- 0.13
</td>

</tr>

</tbody>

</table>

## Clean the workspace to switch to mealworm data

``` r
rm(list = ls())
```

# Analyse of the mealworm data across treatments

## Prepare a phyloseq object for the mealworm data

Use the MAG abundance data and the GTDB-Tk output data to generate a
phyloseq object.

``` r
# Load the OTU table
otumat <- as.matrix(read.table("mealworm_MAG_abundances_depth_absolute.txt",
                        header = TRUE, row.names = 1))

# Load the taxonomy table
taxmat <- as.matrix(read.table("mealworm_MAG_taxonomy.txt",
                        header = TRUE, row.names = 1))

# Load the sample metadata
metadata <- read.table("mealworm_sample_metadata.txt",
                       header = TRUE, row.names = 1)

# Combine the data into a phyloseq object
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, META)

# Convert MAG abundances to percentages
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
```

## Calculate alpha-diversity for the mealworm data and perform statistical analysis

Calculate Shannon and Simpson’s diversity metrics and test if the values
are statistically different between samples. Includes insect batch as a
random effect. Replicate was not included as a random effect as it had
no variance and resulted in non-convergence of the model.

``` r
# Calculate alpha-diversity metrics
richness <- estimate_richness(ps.rel, measures = c("Shannon", "Simpson"))

# Add columns with relevant metadata
richness$treatment <- metadata$Treatment
richness$replicate <- metadata$Replicate
richness$batch <- metadata$Insect_batch
richness$depth <- metadata$Read_count

# Create factors for treatment, ordering the treatments as desired
richness$treatment <- factor(richness$treatment, 
                                 levels = c("Standard", "Polyethylene", "Polystyrene", "Starvation"))

# Perform statistical analysis on the Shannon data
model.shannon <- lmer(log(Shannon) ~ treatment + depth + (1 | batch), data=richness)
summary(model.shannon)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(Shannon) ~ treatment + depth + (1 | batch)
    ##    Data: richness
    ## 
    ## REML criterion at convergence: 5.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.1384 -0.6332 -0.1285  0.4187  1.7872 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  batch    (Intercept) 0.05638  0.2374  
    ##  Residual             0.02399  0.1549  
    ## Number of obs: 22, groups:  batch, 3
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)            0.231927   0.225012   1.031
    ## treatmentPolyethylene  0.162824   0.101788   1.600
    ## treatmentPolystyrene  -0.029917   0.101405  -0.295
    ## treatmentStarvation    0.277123   0.100739   2.751
    ## depth                  0.007724   0.002989   2.584
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) trtmntPlyt trtmntPlys trtmnS
    ## trtmntPlyth -0.375                             
    ## trtmntPlyst -0.357  0.612                      
    ## trtmntStrvt -0.306  0.605      0.606           
    ## depth       -0.713  0.152      0.125      0.051

``` r
model.shannon.aov <- Anova(model.shannon, type = 2)
model.shannon.aov
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: log(Shannon)
    ##             Chisq Df Pr(>Chisq)   
    ## treatment 14.5482  3   0.002246 **
    ## depth      6.6767  1   0.009768 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
model.shannon.bt <- update(ref_grid(model.shannon), tran = "log")
model.shannon.em <- emmeans(model.shannon.bt, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.shannon.em$emmeans
```

    ##  treatment    response    SE   df lower.CL upper.CL null t.ratio p.value
    ##  Standard         1.86 0.294 3.08     1.13     3.06    1   3.921  0.0281
    ##  Polyethylene     2.19 0.331 2.58     1.29     3.71    1   5.185  0.0200
    ##  Polystyrene      1.81 0.273 2.58     1.06     3.06    1   3.911  0.0389
    ##  Starvation       2.45 0.371 2.58     1.45     4.16    1   5.943  0.0144
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
model.shannon.em$contrasts
```

    ##  contrast                   ratio     SE   df lower.CL upper.CL null t.ratio p.value
    ##  Standard / Polyethylene    0.850 0.0867 15.1    0.633    1.140    1  -1.596  0.4098
    ##  Standard / Polystyrene     1.030 0.1050 15.1    0.769    1.381    1   0.294  0.9907
    ##  Standard / Starvation      0.758 0.0764 15.0    0.567    1.014    1  -2.748  0.0642
    ##  Polyethylene / Polystyrene 1.213 0.1080 15.0    0.937    1.569    1   2.154  0.1811
    ##  Polyethylene / Starvation  0.892 0.0803 15.0    0.688    1.156    1  -1.269  0.5950
    ##  Polystyrene / Starvation   0.736 0.0660 15.0    0.568    0.953    1  -3.421  0.0178
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: tukey method for comparing a family of 4 estimates 
    ## Intervals are back-transformed from the log scale 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log scale

``` r
# Perform statistical analysis on the Simpson data
model.simpson <- lmer(log(Simpson) ~ treatment  + depth + (1 | batch), data=richness)
summary(model.simpson)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(Simpson) ~ treatment + depth + (1 | batch)
    ##    Data: richness
    ## 
    ## REML criterion at convergence: 3.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.32466 -0.76106  0.04457  0.48536  1.66878 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  batch    (Intercept) 0.04136  0.2034  
    ##  Residual             0.02173  0.1474  
    ## Number of obs: 22, groups:  batch, 3
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)           -0.557145   0.206213  -2.702
    ## treatmentPolyethylene  0.025380   0.096858   0.262
    ## treatmentPolystyrene  -0.125898   0.096497  -1.305
    ## treatmentStarvation    0.126512   0.095870   1.320
    ## depth                  0.005135   0.002840   1.808
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) trtmntPlyt trtmntPlys trtmnS
    ## trtmntPlyth -0.389                             
    ## trtmntPlyst -0.370  0.612                      
    ## trtmntStrvt -0.317  0.605      0.606           
    ## depth       -0.738  0.151      0.124      0.050

``` r
model.simpson.aov <- Anova(model.simpson, type = 2)
model.simpson.aov
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: log(Simpson)
    ##            Chisq Df Pr(>Chisq)  
    ## treatment 8.8791  3    0.03094 *
    ## depth     3.2692  1    0.07059 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
model.simpson.bt <- update(ref_grid(model.simpson), tran = "log")
model.simpson.em <- emmeans(model.simpson.bt, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.simpson.em$emmeans
```

    ##  treatment    response     SE   df lower.CL upper.CL null t.ratio p.value
    ##  Standard        0.742 0.1030 3.33    0.487     1.13    1  -2.142  0.1125
    ##  Polyethylene    0.761 0.1000 2.72    0.487     1.19    1  -2.071  0.1396
    ##  Polystyrene     0.654 0.0863 2.71    0.418     1.02    1  -3.218  0.0561
    ##  Starvation      0.842 0.1110 2.71    0.539     1.32    1  -1.305  0.2915
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
model.simpson.em$contrasts
```

    ##  contrast                   ratio     SE   df lower.CL upper.CL null t.ratio p.value
    ##  Standard / Polyethylene    0.975 0.0947 15.1    0.737    1.290    1  -0.261  0.9935
    ##  Standard / Polystyrene     1.134 0.1100 15.1    0.858    1.499    1   1.302  0.5758
    ##  Standard / Starvation      0.881 0.0846 15.0    0.668    1.162    1  -1.318  0.5663
    ##  Polyethylene / Polystyrene 1.163 0.0991 15.0    0.910    1.487    1   1.776  0.3218
    ##  Polyethylene / Starvation  0.904 0.0775 15.0    0.706    1.157    1  -1.180  0.6480
    ##  Polystyrene / Starvation   0.777 0.0664 15.0    0.607    0.994    1  -2.955  0.0436
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: tukey method for comparing a family of 4 estimates 
    ## Intervals are back-transformed from the log scale 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log scale

## Plot the alpha-diversity metrics for the mealworm samples

Generate plots showing the Shannon and Simpson’s diversity metrics

``` r
# Collect the emmeans for Shannon data
emms.shannon <- as.data.frame(model.shannon.em$emmeans)
emms.shannon$metric <- 'Shannon'

# Collect the emmeans for teh Simpson's data
emms.simpson <- as.data.frame(model.simpson.em$emmeans)
emms.simpson$metric <- 'Simpson'

# Plot the Shannon data
p.shannon <- ggplot() +
  geom_bar(data = emms.shannon,
           aes(x = treatment, y = response, fill = treatment),
           stat = "identity",
           colour = "black") +
  geom_errorbar(data = emms.shannon,
                aes(x = treatment, ymin = response-SE, ymax = response+SE),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_jitter(data = richness,
              aes(x = treatment, y = Shannon),
              width = 0.15, size = 1.25) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  theme(text = element_text(size = 17), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color="black"),
        axis.text.y = element_text(color="black")) +
  scale_y_continuous(lim=c(0,3.5), breaks=c(0.0,0.7,1.4,2.1,2.8,3.5)) + 
  xlab("") + ylab("")

# Plot the Simpson's data
p.simpson <- ggplot() +
  geom_bar(data = emms.simpson,
           aes(x = treatment, y = response, fill = treatment),
           stat = "identity",
           colour = "black") +
  geom_errorbar(data = emms.simpson,
                aes(x = treatment, ymin = response-SE, ymax = response+SE),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_jitter(data = richness,
              aes(x = treatment, y = Simpson),
              width = 0.15, size = 1.25) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  theme(text = element_text(size = 17), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color="black"),
        axis.text.y = element_text(color="black")) +
  scale_y_continuous(lim=c(0,1), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0)) + 
  xlab("") + ylab("")

# Combine the plots as two panels in one figure
p.shannon + p.simpson
```

![](Statistics_files/figure-gfm/Mealworm%20alpha%20diversity%20plots-1.png)<!-- -->

``` r
# Export a PDF of the figure
pdf('Mealworm.alpha.pdf', width=8.5, height=5)
p.shannon + p.simpson
dev.off()
```

## Calculate and plot beta-diversity for the mealworm data and perform statistical analysis

``` r
# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Treatment, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm   Pr(>F)   
    ## Groups     3 0.13078 0.043593 5.6005  10000 0.005499 **
    ## Residuals 18 0.14011 0.007784                          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Polyethylene Polystyrene Standard Starvation
    ## Polyethylene                 0.543846 0.045895     0.0494
    ## Polystyrene      0.502756             0.026397     0.2488
    ## Standard         0.059182    0.036783              0.0077
    ## Starvation       0.069541    0.242129 0.011910

The dispersion is significantly different between samples, which means
the results of PERMANOVA are not reliable. Nevertheless, run the
PERMANOVA (accounting for treatment and batch) and pair it with a PCoA
to visualize the data.

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Treatment + Insect_batch + Read_count, data = meta, permutations = 10000, by = "margin")
as.data.frame(adonis.bray)
```

    ##              Df  SumOfSqs         R2        F     Pr(>F)
    ## Treatment     3 0.7929616 0.21907371 2.446673 0.00709929
    ## Insect_batch  2 0.6927609 0.19139098 3.206257 0.00319968
    ## Read_count    1 0.2626746 0.07256985 2.431437 0.03799620
    ## Residual     15 1.6204898 0.44769723       NA         NA
    ## Total        21 3.6196110 1.00000000       NA         NA

``` r
# Create a PCoA plot
sample_data(ps.rel)$Treatment <- factor(sample_data(ps.rel)$Treatment, 
                                        levels = c("Standard", "Polyethylene", "Polystyrene", "Starvation"))
pso <- ordinate(ps.rel, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel, pso, type="samples", 
                            color="Treatment", shape="Insect_batch",
                            axes = c(1,2))
p.beta.pcoa <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.40,0.40), breaks=c(-0.40,-0.20,0.00,0.20,0.40),
                     labels = function(x) sprintf("%.2f", x))
p.beta.pcoa
```

![](Statistics_files/figure-gfm/Mealworm%20beta%20diversity%20figure-1.png)<!-- -->

``` r
# Save the figure
pdf('Mealworm.beta.pcoa.pdf', width=7.35, height=5)
p.beta.pcoa
dev.off()
```

Now run a capscale analysis since it is less sensitive to differences in
dispersion and thus should be more appropriate than PERMANOVA in this
case.

``` r
# Set the formula for the model
class(fo <- dist_bc ~ meta$Treatment + meta$Insect_batch + meta$Read_count)
```

    ## [1] "formula"

``` r
# Run the capscale analysis followed by an ANOVA
pso_cap <- ordinate(ps.rel, "CAP", "bray", formula=fo)
anova(pso_cap, by = "margin", type = 2)
```

    ## Permutation test for capscale under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: capscale(formula = OTU ~ meta$Treatment + meta$Insect_batch + meta$Read_count, data = data, distance = distance)
    ##                   Df SumOfSqs      F Pr(>F)   
    ## meta$Treatment     3  0.80758 2.2989  0.005 **
    ## meta$Insect_batch  2  0.69749 2.9783  0.002 **
    ## meta$Read_count    1  0.26533 2.2659  0.048 * 
    ## Residual          15  1.75643                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the relative eigenvalues
all_eig <- eigenvals(pso_cap)
rel_eig <- all_eig / sum(all_eig)
var_explained <- round(rel_eig * 100, 1)

# Plot the capscale results
ordplot5 <- plot_ordination(ps.rel, pso_cap, type="samples", color="Treatment",
                            shape="Insect_batch", axes = c(1,2))
p.beta.cap <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" CAP 1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" CAP 2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-1.50,1.50), breaks=c(-1.50,-0.75,0.00,0.75,1.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-1.50,1.50), breaks=c(-1.50,-0.75,0.00,0.75,1.50),
                     labels = function(x) sprintf("%.2f", x))
p.beta.cap
```

![](Statistics_files/figure-gfm/Mealworm%20beta%20diversity%20capscale%20figure-1.png)<!-- -->

``` r
# Save the figure
pdf('Mealworm.beta.cap.pdf', width=7.35, height=5)
p.beta.cap
dev.off()
```

To mimic what will be done for the superworm data, see what the PCoA and
CAP plots look like if summarized at the genus level

``` r
# Summarize at the genus level, removing taxa not classified at the genus level
ps.rel.genus <- tax_glom(ps.rel, "Genus", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.genus, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.genus))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Treatment, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
    ## Groups     3 0.16242 0.054141 6.1243  10000 0.0031 **
    ## Residuals 18 0.15912 0.008840                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Standard Polyethylene Polystyrene Starvation
    ## Standard                  0.030797    0.023598     0.0049
    ## Polyethylene 0.041019                 0.936506     0.0805
    ## Polystyrene  0.036024     0.917809                 0.0957
    ## Starvation   0.010228     0.094030    0.103670

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Treatment + Insect_batch + Read_count, data = meta, permutations = 10000, by = "margin")
as.data.frame(adonis.bray)
```

    ##              Df  SumOfSqs         R2        F     Pr(>F)
    ## Treatment     3 0.5956565 0.18796774 2.044537 0.03639636
    ## Insect_batch  2 0.6612374 0.20866271 3.404456 0.00469953
    ## Read_count    1 0.2664178 0.08407186 2.743365 0.03119688
    ## Residual     15 1.4567027 0.45968292       NA         NA
    ## Total        21 3.1689293 1.00000000       NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.genus, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.genus, pso, type="samples", 
                            color="Treatment", shape="Insect_batch", 
                            axes = c(1,2))
p.beta.pcoa.genus <- ordplot5 + theme_bw() + 
  theme(text = element_text(size = 17), aspect.ratio = 1) +
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.40,0.40), breaks=c(-0.40,-0.20,0.00,0.20,0.40),
                     labels = function(x) sprintf("%.2f", x))

# Set the formula for the model for the capscale analysis
class(fo <- dist_bc ~ meta$Treatment + meta$Insect_batch + meta$Read_count)
```

    ## [1] "formula"

``` r
# Run the capscale analysis followed by an ANOVA
pso_cap <- ordinate(ps.rel.genus, "CAP", "bray", formula=fo)
anova(pso_cap, by = "margin", type = 2)
```

    ## Permutation test for capscale under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: capscale(formula = OTU ~ meta$Treatment + meta$Insect_batch + meta$Read_count, data = data, distance = distance)
    ##                   Df SumOfSqs      F Pr(>F)   
    ## meta$Treatment     3  0.62002 1.8980  0.031 * 
    ## meta$Insect_batch  2  0.66978 3.0755  0.008 **
    ## meta$Read_count    1  0.27678 2.5418  0.038 * 
    ## Residual          15  1.63336                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the relative eigenvalues
all_eig <- eigenvals(pso_cap)
rel_eig <- all_eig / sum(all_eig)
var_explained <- round(rel_eig * 100, 1)

# Plot the capscale results
ordplot5 <- plot_ordination(ps.rel.genus, pso_cap, type="samples", color="Treatment",
                            shape="Insect_batch", axes = c(1,2))
p.beta.cap.genus <- ordplot5 + theme_bw() + 
  theme(text = element_text(size = 17), aspect.ratio = 1) +
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" CAP 1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" CAP 2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-1.50,1.50), breaks=c(-1.50,-0.75,0.00,0.75,1.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-1.50,1.50), breaks=c(-1.50,-0.75,0.00,0.75,1.50),
                     labels = function(x) sprintf("%.2f", x))
```

``` r
# Display both graphs
p.beta.pcoa.genus + p.beta.cap.genus
```

![](Statistics_files/figure-gfm/Mealworm%20beta%20diversity%20genus%20figure-1.png)<!-- -->

``` r
# Save the figure but just the capscale
pdf('Mealworm.beta.genus.pdf', width=7.35, height=5)
p.beta.cap.genus
dev.off()
```

## Differential abundance analysis of taxa in the mealworm data

Test whether any of the MAGs are enriched in the PE or PS diets relative
to both the starvation and standard diets.

``` r
# Set the factors to get standard as the reference
sample_data(ps)$Treatment <- factor(sample_data(ps)$Treatment, 
                                                 levels = c("Standard", "Polyethylene",
                                                            "Polystyrene","Starvation"))

# Run ANCOMBC2 with standard as the reference
ancombc.standard <- ancombc2(data = ps,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.standard.res <- ancombc.standard$res

# Set the factors to get standard as the reference
sample_data(ps)$Treatment <- factor(sample_data(ps)$Treatment, 
                                                 levels = c("Starvation", "Polyethylene",
                                                            "Polystyrene","Standard"))

# Run ANCOMBC2 with standard as the reference
ancombc.starvation <- ancombc2(data = ps,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.starvation.res <- ancombc.starvation$res

# Get differentialy abundant taxa for PE
taxa.diff_abn.PE = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.PE <- rbind(taxa.diff_abn.PE, res.asv$taxon[n])
  }
}

# Print the results for PE
if(nrow(taxa.diff_abn.PE) == 0) {
  print("No MAGs met the criteria for differential abundance for PE")
} else {
  print(taxa.diff_abn.PE)
}
```

    ## [1] "No MAGs met the criteria for differential abundance for PE"

``` r
# Get differentialy abundant taxa for PS
taxa.diff_abn.PS = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.PS <- rbind(taxa.diff_abn.PS, res.asv$taxon[n])
  }
}

# Print the results for PS
if(nrow(taxa.diff_abn.PS) == 0) {
  print("No MAGs met the criteria for differential abundance for PS")
} else {
  print(taxa.diff_abn.PS)
}
```

    ## [1] "No MAGs met the criteria for differential abundance for PS"

Summarize the data at the species level and then test whether any
species are enriched in the PE or PS diets relative to both the
starvation and standard diets.

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.species <- tax_glom(ps, "Species", NArm = TRUE)

# Set the factors to get standard as the reference
sample_data(ps.species)$Treatment <- factor(sample_data(ps.species)$Treatment, 
                                                 levels = c("Standard", "Polyethylene",
                                                            "Polystyrene","Starvation"))

# Run ANCOMBC2 with standard as the reference
ancombc.standard <- ancombc2(data = ps.species,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.standard.res <- ancombc.standard$res

# Set the factors to get standard as the reference
sample_data(ps.species)$Treatment <- factor(sample_data(ps.species)$Treatment, 
                                                 levels = c("Starvation", "Polyethylene",
                                                            "Polystyrene","Standard"))

# Run ANCOMBC2 with standard as the reference
ancombc.starvation <- ancombc2(data = ps.species,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.starvation.res <- ancombc.starvation$res

# Get differentialy abundant taxa for PE
taxa.diff_abn.species.PE = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.species.PE <- rbind(taxa.diff_abn.species.PE, res.asv$taxon[n])
  }
}

# Print the results for PE
if(nrow(taxa.diff_abn.species.PE) == 0) {
  print("No species met the criteria for differential abundance for PE")
} else {
  print(taxa.diff_abn.species.PE)
}
```

    ## [1] "No species met the criteria for differential abundance for PE"

``` r
# Get differentialy abundant taxa for PS
taxa.diff_abn.species.PS = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.species.PS <- rbind(taxa.diff_abn.species.PS, res.asv$taxon[n])
  }
}

# Print the results for PS
if(nrow(taxa.diff_abn.species.PS) == 0) {
  print("No species met the criteria for differential abundance for PS")
} else {
  print(taxa.diff_abn.species.PS)
}
```

    ## [1] "No species met the criteria for differential abundance for PS"

## Clean the workspace to switch to superworm data

``` r
rm(list = ls())
```

# Analyse of the superworm data across treatments

## Prepare a phyloseq object for the superworm data

Use the MAG abundance data and the GTDB-Tk output data to generate a
phyloseq object.

``` r
# Load the OTU table
otumat <- as.matrix(read.table("superworm_MAG_abundances_depth_absolute.txt",
                        header = TRUE, row.names = 1))

# Load the taxonomy table
taxmat <- as.matrix(read.table("superworm_MAG_taxonomy.txt",
                        header = TRUE, row.names = 1))

# Load the sample metadata
metadata <- read.table("superworm_sample_metadata.txt",
                       header = TRUE, row.names = 1)

# Combine the data into a phyloseq object
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(metadata)
ps <- phyloseq(OTU, TAX, META)

# Convert MAG abundances to percentages
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
```

## Calculate alpha-diversity for the superworm data and perform statistical analysis

Calculate Shannon and Simpson’s diversity metrics and test if the values
are statistically different between samples. Includes insect batch as a
random effect. Replicate was not included as a random effect as it had
no variance and resulted in non-convergence of the model.

``` r
# Calculate alpha-diversity metrics
richness <- estimate_richness(ps.rel, measures = c("Shannon", "Simpson"))

# Add columns with relevant metadata
richness$treatment <- metadata$Treatment
richness$replicate <- metadata$Replicate
richness$batch <- metadata$Insect_batch
richness$depth <- metadata$Read_count

# Create factors for treatment, ordering the treatments as desired
richness$treatment <- factor(richness$treatment, 
                                 levels = c("Standard", "Polyethylene", "Polystyrene", "Starvation"))

# Perform statistical analysis on the Shannon data
model.shannon <- lmer(log(Shannon) ~ treatment + depth + (1 | batch), data=richness)
summary(model.shannon)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(Shannon) ~ treatment + depth + (1 | batch)
    ##    Data: richness
    ## 
    ## REML criterion at convergence: 4.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5147 -0.4980  0.1155  0.5857  1.5856 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  batch    (Intercept) 0.008477 0.09207 
    ##  Residual             0.029114 0.17063 
    ## Number of obs: 24, groups:  batch, 3
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)            0.40452    0.16318   2.479
    ## treatmentPolyethylene  0.20378    0.09866   2.066
    ## treatmentPolystyrene   0.21379    0.09853   2.170
    ## treatmentStarvation    0.11262    0.09862   1.142
    ## depth                  0.01085    0.00335   3.237
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) trtmntPlyt trtmntPlys trtmnS
    ## trtmntPlyth -0.347                             
    ## trtmntPlyst -0.287  0.498                      
    ## trtmntStrvt -0.341  0.501      0.499           
    ## depth       -0.844  0.054     -0.018      0.046

``` r
model.shannon.aov <- Anova(model.shannon, type = 2)
model.shannon.aov
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: log(Shannon)
    ##             Chisq Df Pr(>Chisq)   
    ## treatment  6.1075  3   0.106495   
    ## depth     10.4799  1   0.001207 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
model.shannon.bt <- update(ref_grid(model.shannon), tran = "log")
model.shannon.em <- emmeans(model.shannon.bt, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.shannon.em$emmeans
```

    ##  treatment    response    SE   df lower.CL upper.CL null t.ratio p.value
    ##  Standard         2.32 0.204 6.41     1.88     2.87    1   9.625 <0.0001
    ##  Polyethylene     2.85 0.250 6.42     2.31     3.52    1  11.944 <0.0001
    ##  Polystyrene      2.88 0.253 6.43     2.33     3.56    1  12.055 <0.0001
    ##  Starvation       2.60 0.228 6.42     2.11     3.21    1  10.908 <0.0001
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
model.shannon.em$contrasts
```

    ##  contrast                   ratio     SE df lower.CL upper.CL null t.ratio p.value
    ##  Standard / Polyethylene    0.816 0.0805 17    0.616     1.08    1  -2.065  0.2041
    ##  Standard / Polystyrene     0.808 0.0796 17    0.610     1.07    1  -2.170  0.1716
    ##  Standard / Starvation      0.893 0.0881 17    0.675     1.18    1  -1.142  0.6697
    ##  Polyethylene / Polystyrene 0.990 0.0978 17    0.748     1.31    1  -0.101  0.9996
    ##  Polyethylene / Starvation  1.095 0.1080 17    0.828     1.45    1   0.925  0.7919
    ##  Polystyrene / Starvation   1.106 0.1090 17    0.836     1.46    1   1.025  0.7376
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: tukey method for comparing a family of 4 estimates 
    ## Intervals are back-transformed from the log scale 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log scale

``` r
# Perform statistical analysis on the Simpson data
model.simpson <- lmer(log(Simpson) ~ treatment + depth + (1 | batch), data=richness)
summary(model.simpson)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(Simpson) ~ treatment + depth + (1 | batch)
    ##    Data: richness
    ## 
    ## REML criterion at convergence: -17
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.84445 -0.34298  0.01028  0.63900  1.23334 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance  Std.Dev.
    ##  batch    (Intercept) 0.0008618 0.02936 
    ##  Residual             0.0102772 0.10138 
    ## Number of obs: 24, groups:  batch, 3
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)           -0.565861   0.091982  -6.152
    ## treatmentPolyethylene  0.091003   0.058613   1.553
    ## treatmentPolystyrene   0.106227   0.058539   1.815
    ## treatmentStarvation    0.002045   0.058591   0.035
    ## depth                  0.007774   0.001956   3.975
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) trtmntPlyt trtmntPlys trtmnS
    ## trtmntPlyth -0.364                             
    ## trtmntPlyst -0.303  0.498                      
    ## trtmntStrvt -0.358  0.501      0.499           
    ## depth       -0.874  0.053     -0.018      0.046

``` r
model.simpson.aov <- Anova(model.simpson, type = 2)
model.simpson.aov
```

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: log(Simpson)
    ##             Chisq Df Pr(>Chisq)    
    ## treatment  5.6281  3     0.1312    
    ## depth     15.7990  1  7.044e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
model.simpson.bt <- update(ref_grid(model.simpson), tran = "log")
model.simpson.em <- emmeans(model.simpson.bt, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.simpson.em$emmeans
```

    ##  treatment    response     SE   df lower.CL upper.CL null t.ratio p.value
    ##  Standard        0.778 0.0348 10.9    0.705    0.859    1  -5.612  0.0002
    ##  Polyethylene    0.852 0.0381 10.9    0.772    0.940    1  -3.576  0.0044
    ##  Polystyrene     0.865 0.0387 10.9    0.784    0.955    1  -3.235  0.0080
    ##  Starvation      0.780 0.0349 10.9    0.706    0.860    1  -5.565  0.0002
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## Tests are performed on the log scale

``` r
model.simpson.em$contrasts
```

    ##  contrast                   ratio     SE df lower.CL upper.CL null t.ratio p.value
    ##  Standard / Polyethylene    0.913 0.0535 17    0.773     1.08    1  -1.552  0.4302
    ##  Standard / Polystyrene     0.899 0.0526 17    0.761     1.06    1  -1.815  0.3007
    ##  Standard / Starvation      0.998 0.0585 17    0.845     1.18    1  -0.035  1.0000
    ##  Polyethylene / Polystyrene 0.985 0.0578 17    0.834     1.16    1  -0.259  0.9937
    ##  Polyethylene / Starvation  1.093 0.0640 17    0.926     1.29    1   1.520  0.4480
    ##  Polystyrene / Starvation   1.110 0.0651 17    0.939     1.31    1   1.776  0.3181
    ## 
    ## Degrees-of-freedom method: kenward-roger 
    ## Confidence level used: 0.95 
    ## Conf-level adjustment: tukey method for comparing a family of 4 estimates 
    ## Intervals are back-transformed from the log scale 
    ## P value adjustment: tukey method for comparing a family of 4 estimates 
    ## Tests are performed on the log scale

## Plot the alpha-diversity metrics for the superworm samples

Generate plots showing the Shannon and Simpson’s diversity metrics

``` r
# Collect the emmeans for Shannon data
emms.shannon <- as.data.frame(model.shannon.em$emmeans)
emms.shannon$metric <- 'Shannon'

# Collect the emmeans for teh Simpson's data
emms.simpson <- as.data.frame(model.simpson.em$emmeans)
emms.simpson$metric <- 'Simpson'

# Plot the Shannon data
p.shannon <- ggplot() +
  geom_bar(data = emms.shannon,
           aes(x = treatment, y = response, fill = treatment),
           stat = "identity",
           colour = "black") +
  geom_errorbar(data = emms.shannon,
                aes(x = treatment, ymin = response-SE, ymax = response+SE),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_jitter(data = richness,
              aes(x = treatment, y = Shannon),
              width = 0.15, size = 1.25) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  theme(text = element_text(size = 17), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color="black"),
        axis.text.y = element_text(color="black")) +
  scale_y_continuous(lim=c(0,3.5), breaks=c(0.0,0.7,1.4,2.1,2.8,3.5)) + 
  xlab("") + ylab("")

# Plot the Simpson's data
p.simpson <- ggplot() +
  geom_bar(data = emms.simpson,
           aes(x = treatment, y = response, fill = treatment),
           stat = "identity",
           colour = "black") +
  geom_errorbar(data = emms.simpson,
                aes(x = treatment, ymin = response-SE, ymax = response+SE),
                stat = "identity", width = 0.5, size = 0.5) +
  geom_jitter(data = richness,
              aes(x = treatment, y = Simpson),
              width = 0.15, size = 1.25) +
  facet_wrap(~metric, scales="free_y") +
  theme_bw() +
  theme(text = element_text(size = 17), legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face="bold"),
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, color="black"),
        axis.text.y = element_text(color="black")) +
  scale_y_continuous(lim=c(0,1), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0)) + 
  xlab("") + ylab("")

# Combine the plots as two panels in one figure
p.shannon + p.simpson
```

![](Statistics_files/figure-gfm/Superworm%20alpha%20diversity%20plots-1.png)<!-- -->

``` r
# Export a PDF of the figure
pdf('Superworm.alpha.pdf', width=8.5, height=5)
p.shannon + p.simpson
dev.off()
```

## Calculate and plot beta-diversity for the superworm data and perform statistical analysis

``` r
# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Treatment, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
    ## Groups     3 0.010225 0.0034085 0.383  10000  0.767
    ## Residuals 20 0.177967 0.0088983                    
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Polyethylene Polystyrene Standard Starvation
    ## Polyethylene                  0.54045  0.50445     0.7909
    ## Polystyrene       0.53713              0.82882     0.4820
    ## Standard          0.49579     0.82759              0.4432
    ## Starvation        0.78884     0.47717  0.44243

Unlike the mealworm data, the dispersion is not significantly different
between samples, which means the results of PERMANOVA should be
reliable.

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Treatment + Insect_batch + Read_count, data = meta, permutations = 10000, by = "margin")
as.data.frame(adonis.bray)
```

    ##              Df  SumOfSqs        R2        F     Pr(>F)
    ## Treatment     3 0.4129333 0.1314637 1.909898 0.02239776
    ## Insect_batch  2 0.9924674 0.3159674 6.885537 0.00009999
    ## Read_count    1 0.3070895 0.0977667 4.261049 0.00089991
    ## Residual     17 1.2251728 0.3900528       NA         NA
    ## Total        23 3.1410438 1.0000000       NA         NA

``` r
# Create a PCoA plot
sample_data(ps.rel)$Treatment <- factor(sample_data(ps.rel)$Treatment, 
                                        levels = c("Standard", "Polyethylene", "Polystyrene", "Starvation"))
pso <- ordinate(ps.rel, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel, pso, type="samples", 
                            color="Treatment", shape="Insect_batch",
                            axes = c(1,2))
p.beta.pcoa <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.40,0.40), breaks=c(-0.40,-0.20,0.00,0.20,0.40),
                     labels = function(x) sprintf("%.2f", x))
p.beta.pcoa
```

![](Statistics_files/figure-gfm/Superworm%20beta%20diversity%20figure-1.png)<!-- -->

``` r
# Save the figure
pdf('Superworm.beta.pcoa.pdf', width=7.3, height=5)
p.beta.pcoa
dev.off()
```

Now run a capscale analysis since it is less sensitive to differences in
dispersion and thus should be more appropriate than PERMANOVA in this
case.

``` r
# Set the formula for the model
class(fo <- dist_bc ~ meta$Treatment + meta$Insect_batch + meta$Read_count)
```

    ## [1] "formula"

``` r
# Run the capscale analysis followed by an ANOVA
pso_cap <- ordinate(ps.rel, "CAP", "bray", formula=fo)
anova(pso_cap, by = "margin", type = 2)
```

    ## Permutation test for capscale under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: capscale(formula = OTU ~ meta$Treatment + meta$Insect_batch + meta$Read_count, data = data, distance = distance)
    ##                   Df SumOfSqs      F Pr(>F)    
    ## meta$Treatment     3  0.42078 1.8368  0.031 *  
    ## meta$Insect_batch  2  0.99767 6.5327  0.001 ***
    ## meta$Read_count    1  0.30824 4.0366  0.002 ** 
    ## Residual          17  1.29813                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the relative eigenvalues
all_eig <- eigenvals(pso_cap)
rel_eig <- all_eig / sum(all_eig)
var_explained <- round(rel_eig * 100, 1)

# Plot the capscale results
ordplot5 <- plot_ordination(ps.rel, pso_cap, type="samples", color="Treatment",
                            shape="Insect_batch", axes = c(1,2))
p.beta.cap <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" CAP 1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" CAP 2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-1.00,1.40), breaks=c(-1.00,-0.40,0.20,0.80,1.40),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-1.4,1.0), breaks=c(-1.40,-0.80,-0.20,0.40,1.00),
                     labels = function(x) sprintf("%.2f", x))
p.beta.cap
```

![](Statistics_files/figure-gfm/Superworm%20beta%20diversity%20capscale%20figure-1.png)<!-- -->

``` r
# Save the figure
pdf('Superworm.beta.cap.pdf', width=7.34, height=5)
p.beta.cap
dev.off()
```

Given the strong impact of batch, see what the PCoA and CAP plots look
like if summarized at the genus level

``` r
# Summarize at the genus level, removing taxa not classified at the genus level
ps.rel.genus <- tax_glom(ps.rel, "Genus", NArm = TRUE)

# Calculate the bray curtis distances
dist_bc <- phyloseq::distance(ps.rel.genus, method = "bray")

# Extract sample data
meta <- data.frame(sample_data(ps.rel.genus))

# Perform a permutation test to test if the data meet the assumptions of PERMANOVA
beta.bray <- betadisper(dist_bc, meta$Treatment, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 10000
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     3 0.020585 0.0068617 0.7232  10000 0.5568
    ## Residuals 20 0.189765 0.0094883                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##              Standard Polyethylene Polystyrene Starvation
    ## Standard                   0.40756     0.52055     0.2329
    ## Polyethylene  0.40113                  0.94771     0.4506
    ## Polystyrene   0.51580      0.94916                 0.5196
    ## Starvation    0.23234      0.45331     0.52514

``` r
# Run the PERMANOVA test
adonis.bray <- adonis2(dist_bc ~ Treatment + Insect_batch + Read_count, data = meta, permutations = 10000, by = "margin")
as.data.frame(adonis.bray)
```

    ##              Df  SumOfSqs        R2        F     Pr(>F)
    ## Treatment     3 0.2770033 0.1182068 1.753170 0.07629237
    ## Insect_batch  2 0.6483302 0.2766646 6.154979 0.00009999
    ## Read_count    1 0.2881039 0.1229437 5.470279 0.00139986
    ## Residual     17 0.8953412 0.3820726       NA         NA
    ## Total        23 2.3433801 1.0000000       NA         NA

``` r
# Create a PCoA plot
pso <- ordinate(ps.rel.genus, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel.genus, pso, type="samples", 
                            color="Treatment", shape="Insect_batch", 
                            axes = c(1,2))
p.beta.pcoa.genus <- ordplot5 + theme_bw() + 
  theme(text = element_text(size = 17), aspect.ratio = 1) +
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" PC1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" PC2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-0.50,0.50), breaks=c(-0.50,-0.25,0.00,0.25,0.50),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-0.2,0.4), breaks=c(-0.20,-0.05,0.10,0.25,0.40),
                     labels = function(x) sprintf("%.2f", x))

# Set the formula for the model for the capscale analysis
class(fo <- dist_bc ~ meta$Treatment + meta$Insect_batch + meta$Read_count)
```

    ## [1] "formula"

``` r
# Run the capscale analysis followed by an ANOVA
pso_cap <- ordinate(ps.rel.genus, "CAP", "bray", formula=fo)
anova(pso_cap, by = "margin", type = 2)
```

    ## Permutation test for capscale under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: capscale(formula = OTU ~ meta$Treatment + meta$Insect_batch + meta$Read_count, data = data, distance = distance)
    ##                   Df SumOfSqs      F Pr(>F)    
    ## meta$Treatment     3  0.28249 1.6562  0.092 .  
    ## meta$Insect_batch  2  0.65447 5.7557  0.001 ***
    ## meta$Read_count    1  0.29009 5.1023  0.004 ** 
    ## Residual          17  0.96652                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the relative eigenvalues
all_eig <- eigenvals(pso_cap)
rel_eig <- all_eig / sum(all_eig)
var_explained <- round(rel_eig * 100, 1)

# Plot the capscale results
ordplot5 <- plot_ordination(ps.rel.genus, pso_cap, type="samples", color="Treatment",
                            shape="Insect_batch", axes = c(1,2))
p.beta.cap.genus <- ordplot5 + theme_bw() + 
  theme(text = element_text(size = 17), aspect.ratio = 1) +
  geom_point(size = 7) +
  scale_shape_manual(values=c(16,15,17)) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  labs(x = paste0(" CAP 1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" CAP 2 [", round(var_explained[2], 1), "%]")) +
  scale_x_continuous(lim=c(-1.40,1.40), breaks=c(-1.40,-0.70,0.00,0.70,1.40),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_y_continuous(lim=c(-1.6,0.8), breaks=c(-1.60,-1.00,-0.40,0.20,0.80),
                     labels = function(x) sprintf("%.2f", x))
```

``` r
# Display both graphs
p.beta.pcoa.genus + p.beta.cap.genus
```

![](Statistics_files/figure-gfm/Superworm%20beta%20diversity%20genus%20figure-1.png)<!-- -->

``` r
# Save the figure but just the capscale
pdf('Superworm.beta.genus.pdf', width=7.35, height=5)
p.beta.cap.genus
dev.off()
```

## Differential abundance analysis of taxa in the superworm data

``` r
# Set the factors to get standard as the reference
sample_data(ps)$Treatment <- factor(sample_data(ps)$Treatment, 
                                                 levels = c("Standard", "Polyethylene",
                                                            "Polystyrene","Starvation"))

# Run ANCOMBC2 with standard as teh reference
ancombc.standard <- ancombc2(data = ps,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.standard.res <- ancombc.standard$res

# Set the factors to get standard as the reference
sample_data(ps)$Treatment <- factor(sample_data(ps)$Treatment, 
                                                 levels = c("Starvation", "Polyethylene",
                                                            "Polystyrene","Standard"))

# Run ANCOMBC2 with standard as the reference
ancombc.starvation <- ancombc2(data = ps,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch/Replicate)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.starvation.res <- ancombc.starvation$res

# Get differentialy abundant taxa for PE
taxa.diff_abn.PE = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.PE <- rbind(taxa.diff_abn.PE, res.asv$taxon[n])
  }
}

# Print the results for PE
if(nrow(taxa.diff_abn.PE) == 0) {
  print("No MAGs met the criteria for differential abundance for PE")
} else {
  print(taxa.diff_abn.PE)
}
```

    ## [1] "No MAGs met the criteria for differential abundance for PE"

``` r
# Get differentialy abundant taxa for PS
taxa.diff_abn.PS = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.PS <- rbind(taxa.diff_abn.PS, res.asv$taxon[n])
  }
}

# Print the results for PS
if(nrow(taxa.diff_abn.PS) == 0) {
  print("No MAGs met the criteria for differential abundance for PS")
} else {
  print(taxa.diff_abn.PS)
}
```

    ## [1] "No MAGs met the criteria for differential abundance for PS"

Summarize the data at the species level and then test whether any
species are enriched in the PE or PS diets relative to both the
starvation and standard diets.

``` r
# Summarize at the species level, removing taxa not classified at the species level
ps.species <- tax_glom(ps.rel, "Species", NArm = TRUE)

# Set the factors to get standard as the reference
sample_data(ps.species)$Treatment <- factor(sample_data(ps.species)$Treatment, 
                                                 levels = c("Standard", "Polyethylene",
                                                            "Polystyrene","Starvation"))

# Run ANCOMBC2 with standard as teh reference
ancombc.standard <- ancombc2(data = ps.species,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.standard.res <- ancombc.standard$res

# Set the factors to get standard as the reference
sample_data(ps.species)$Treatment <- factor(sample_data(ps.species)$Treatment, 
                                                 levels = c("Starvation", "Polyethylene",
                                                            "Polystyrene","Standard"))

# Run ANCOMBC2 with standard as the reference
ancombc.starvation <- ancombc2(data = ps.species,
                             fix_formula = "Treatment", 
                             rand_formula = "(1 | Insect_batch)",
                             p_adj_method = "BH",
                             alpha = 0.025,
                             global = FALSE)
ancombc.starvation.res <- ancombc.starvation$res

# Get differentialy abundant taxa for PE
taxa.diff_abn.species.PE = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.species.PE <- rbind(taxa.diff_abn.species.PE, res.asv$taxon[n])
  }
}

# Print the results for PE
if(nrow(taxa.diff_abn.species.PE) == 0) {
  print("No species met the criteria for differential abundance for PE")
} else {
  print(taxa.diff_abn.species.PE)
}
```

    ## [1] "No species met the criteria for differential abundance for PE"

``` r
# Get differentialy abundant taxa for PS
taxa.diff_abn.species.PS = data.frame()
for(n in 1:nrow(ancombc.standard.res)) { 
  if(ancombc.standard.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.standard.res$passed_ss_TreatmentPolyethylene[n] == TRUE &
     ancombc.starvation.res$diff_TreatmentPolyethylene[n] == TRUE & 
     ancombc.starvation.res$passed_ss_TreatmentPolyethylene[n] == TRUE) { 
    taxa.diff_abn.species.PS <- rbind(taxa.diff_abn.species.PS, res.asv$taxon[n])
  }
}

# Print the results for PS
if(nrow(taxa.diff_abn.species.PS) == 0) {
  print("No species met the criteria for differential abundance for PS")
} else {
  print(taxa.diff_abn.species.PS)
}
```

    ## [1] "No species met the criteria for differential abundance for PS"

## Clean the workspace to switch to the antiSMASH analysis

``` r
rm(list = ls())
```

# Analysis of the antiSMASH output

## Create a table as a binary presence/absence matrix to summarize the antismash output

``` r
# Import the data
data <- read.table("antismash.summary.tsv")

# Remove the contig column
data$V2 <- NULL

# Get just the unique rows to not count twice BGCs split across contigs
data.unique <- unique(data)

# Now merge the class and product columns and remove unwanted columns
data.unique$V5 <- paste0(data.unique$V3, '__', data.unique$V4)
data.unique$V3 <- NULL
data.unique$V4 <- NULL

# Convert into a binary matrix
data.matrix <- as.data.frame.matrix(table(data.unique))

# Save the matrix
write.csv(data.matrix, file="antiSMASH.matrix.csv")
```

## Clean the workspace to switch to CFU analysis

``` r
rm(list = ls())
```

# Analysis of CFU counts per insect

``` r
# Load the data
data <- read.table("cfu.counts.tsv", header = TRUE)
```

First run the t-test on just the standard diet samples

``` r
# Paired t-test for the standard diet
t.test(data$Ratio[1:3], data$Ratio[4:6], paired=TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data$Ratio[1:3] and data$Ratio[4:6]
    ## t = 2.2984, df = 2, p-value = 0.1483
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -415665.7 1368999.0
    ## sample estimates:
    ## mean difference 
    ##        476666.7

Now run the t-test on just the PS diet samples

``` r
# Paired t-test for the plastic diet
t.test(data$Ratio[7:9], data$Ratio[10:12], paired=TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data$Ratio[7:9] and data$Ratio[10:12]
    ## t = -1.1433, df = 2, p-value = 0.3713
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -667443.1  387194.5
    ## sample estimates:
    ## mean difference 
    ##       -140124.3

There is a clear outlier in the PS diet samples, so remove it and rerun
the t-test in case it has an impact

``` r
# Remove outlier
data.cleaned <- data
data.cleaned <- data.cleaned[-c(8), ]
data.cleaned <- data.cleaned[-c(10), ]

# Paired t-test for the plastic diet
t.test(data.cleaned$Ratio[7:8], data.cleaned$Ratio[9:10], paired=TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data.cleaned$Ratio[7:8] and data.cleaned$Ratio[9:10]
    ## t = -0.32782, df = 1, p-value = 0.7983
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -1615885  1534603
    ## sample estimates:
    ## mean difference 
    ##       -40641.03

Given there is no effect per diet, check if combining the data gives
more power to detect an effect

``` r
# Rearrange the data to make this easier to run
data.rearranged <- data %>% arrange(desc(Spin_condition))

# Paired t-test overall
t.test(data.rearranged$Ratio[1:6], data.rearranged$Ratio[7:12], paired=TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data.rearranged$Ratio[1:6] and data.rearranged$Ratio[7:12]
    ## t = 0.96151, df = 5, p-value = 0.3805
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -281601.3  618143.6
    ## sample estimates:
    ## mean difference 
    ##        168271.2

Still not statistically significant; confirm this is true with the
outlier removed

``` r
# Remove outlier
data.rearranged.cleaned <- data.rearranged
data.rearranged.cleaned <- data.rearranged.cleaned[-c(5), ]
data.rearranged.cleaned <- data.rearranged.cleaned[-c(10), ]

# Paired t-test overall
t.test(data.rearranged.cleaned$Ratio[1:5], data.rearranged.cleaned$Ratio[6:10], paired=TRUE)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  data.rearranged.cleaned$Ratio[1:5] and data.rearranged.cleaned$Ratio[6:10]
    ## t = 1.5446, df = 4, p-value = 0.1973
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -215115.7  754602.9
    ## sample estimates:
    ## mean difference 
    ##        269743.6

Based on these results, while there is a massive difference in bacterial
versus fungal CFU, there is no statistically significant effect of the
centrifugation step.
