# Load libraries
library(lme4)
library(phyloseq)
library(dplyr)
library(car)
library(emmeans)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(multcompView)
library(patchwork)
library(vegan)

# Load the OTU table
otumat <- as.matrix(read.table("mealworm_MAG_abundances_depth.txt",
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

# Calculate alpha-diversity metrics
richness <- estimate_richness(ps, measures = c("Shannon", "InvSimpson", "Simpson"))
richness$treatment <- metadata$Treatment
richness$replicate <- metadata$Replicate
richness$batch <- metadata$Insect_batch
richness$treatment <- factor(richness$treatment, 
                                 levels = c("Standard", "Polyethylene", "Polystyrene", "Starvation"))

# Do the statistics on the Shannon
model.shannon <- lmer(log(Shannon) ~ treatment + (1 | batch/replicate), data=richness)
summary(model.shannon)
model.shannon.aov <- Anova(model.shannon, type = 2)
model.shannon.em <- emmeans(model.shannon, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.shannon.em$emmeans
model.shannon.em$contrasts

# Do the statistics on the Simpson
model.simpson <- lmer(log(Simpson) ~ treatment + (1 | batch/replicate), data=richness)
summary(model.simpson)
model.simpson.aov <- Anova(model.simpson, type = 2)
model.simpson.em <- emmeans(model.simpson, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.simpson.em$emmeans
model.simpson.em$contrasts

# Do the statistics on the Inverse Simpson
model.invsimpson <- lmer(log(InvSimpson) ~ treatment + (1 | batch/replicate), data=richness)
summary(model.invsimpson)
model.invsimpson.aov <- Anova(model.invsimpson, type = 2)
model.invsimpson.em <- emmeans(model.invsimpson, pairwise ~ treatment,
                            infer = TRUE, type = "response")
model.invsimpson.em$emmeans
model.invsimpson.em$contrasts

# Collect the emmeans for all metrics
emms.shannon <- as.data.frame(model.shannon.em$emmeans)
emms.shannon$metric <- 'Shannon'
emms.simpson <- as.data.frame(model.simpson.em$emmeans)
emms.simpson$metric <- 'Simpson'
emms <- rbind(emms.shannon, emms.simpson)

# Plot the alpha diversity
p1 <- ggplot() +
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
  scale_y_continuous(lim=c(0,3), breaks=c(0.0,0.6,1.2,1.8,2.4,3.0)) + 
  xlab("") + ylab("")
p2 <- ggplot() +
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
p1 + p2
pdf('Mealworm.alpha.pdf', width=10, height=5)
p1 + p2
dev.off()











# Calculate the beta-diversity
dist_bc <- phyloseq::distance(ps.rel, method = "bray")
meta <- data.frame(sample_data(ps.rel))
beta.bray <- betadisper(dist_bc, meta$Treatment, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
adonis.bray <- adonis2(dist_bc ~ Treatment * Insect_batch, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)

# Beta diversity plot - PCoA
class(fo <- dist_bc ~ meta$Treatment * meta$Insect_batch/meta$Replicate)
pso <- ordinate(ps.rel, "PCoA", "bray", formula=fo)

pso <- ordinate(ps.rel, "PCoA", "bray")
var_explained <- pso$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel, pso, type="samples", color="Treatment", axes = c(1,2))
plot_beta_diversity <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 5) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = paste0(" Axis 1 [", round(var_explained[1], 1), "%]"),
       y = paste0(" Axis 2 [", round(var_explained[2], 1), "%]"))
plot_beta_diversity
pdf('Sed16S.beta.pcoa.pdf', width=6.75, height=5)
plot_beta_diversity
dev.off()

cap_result <- capscale(dist_bc ~ Treatment + Insect_batch, data = meta)
# Beta diversity plot

ps.rel.genus <- tax_glom(ps.rel, "Species", NArm = FALSE)


class(fo <- dist_bc ~ meta$Treatment * meta$Insect_batch)
pso_cap <- ordinate(ps.rel, "CAP", "bray", formula=fo)
var_explained <- pso_cap$values$Relative_eig * 100
ordplot5 <- plot_ordination(ps.rel, pso_cap, type="samples", color="Treatment",
                            shape="Insect_batch", axes = c(1,2))
plot_beta_diversity3 <- ordplot5 +  theme_bw() + theme(text = element_text(size = 17), aspect.ratio = 1) + 
  geom_point(size = 5) +
  theme_bw() +
  theme(text = element_text(size = 17),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_beta_diversity3










# Stats for beta diversity
dist_bc <- phyloseq::distance(sediment.physeqR.relative, method = "bray")
meta <- data.frame(sample_data(sediment.physeqR.relative))
beta.bray <- betadisper(dist_bc, meta$dose, type = "centroid", bias.adjust = TRUE)
permdisp <- permutest(beta.bray, pairwise = TRUE, permutations = 10000)
permdisp
adonis.bray <- adonis2(dist_bc ~ dose, data = meta, permutations = 10000, by = "terms")
as.data.frame(adonis.bray)
otu_table <- data.frame(otu_table(sediment.physeqR.relative))
otu_table.t <- t(otu_table)
nmds <- metaMDS(otu_table.t, k = 2, trymax = 100)
plot(nmds)
fit <- envfit(nmds, meta$dose, permutations = 10000,
              na.rm = TRUE)
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$Study_ID <- rownames(nmds_scores)
fit$vectors$r
fit$vectors$pvals
write.csv(as.data.frame(adonis.bray), file='Sed16S.adonis.csv')
write.csv(as.data.frame(fit$vectors$pvals), file='Sed16S.nmds.csv')
write.csv(as.data.frame(permdisp$tab), file='Sed16S.permdisp.csv')
cap_result <- capscale(dist_bc ~ dose, data = meta)
(total_constrained <- cap_result$CCA$tot.chi)
anova_terms <- anova(cap_result, by = "term")
sum_sqs <- anova_terms$"SumOfSqs"
total_variation <- sum(sum_sqs)
prop_explained <- sum_sqs / total_variation
explained_df <- data.frame(
  Term = rownames(anova_terms),
  SumOfSqs = sum_sqs,
  Proportion = prop_explained
)
print(explained_df)
write.csv(anova_terms, file='Sed16S.cap.anova.csv')
write.csv(explained_df, file='Sed16S.cap.explained.csv')




