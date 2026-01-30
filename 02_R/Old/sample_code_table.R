

# glom data at genus level 
ps.rel.genus <- tax_glom(ps.rel, "Order", NArm = FALSE)

# Make bubble plot at the phylum level
top.taxa = names(sort(taxa_sums(ps.rel.genus), TRUE)[1:20])
top.taxa.table = cbind(tax_table(ps.rel.genus), topTaxa = NA)
top.taxa.table[top.taxa, "topTaxa"] <- as(tax_table(ps.rel.genus)[top.taxa, "Order"], 
                                          "character")
tax_table(ps.rel.genus) <- tax_table(top.taxa.table)
ps.rel.genus.glom <- psmelt(ps.rel.genus) # create dataframe from phyloseq object
ps.rel.genus.glom$Phylum <- as.character(ps.rel.genus.glom$Order) #convert to character
ps.rel.genus.glom.noNA <- na.omit(ps.rel.genus.glom)
ps.rel.genus.glom <- ps.rel.genus.glom.noNA


plot_bubble = ggplot(ps.rel.genus.glom, aes(x = Sample, y = topTaxa)) + 
  geom_point(aes(size = Abundance*100, fill=Insect), alpha = 1, shape = 21) + # abundance *100 to convert to percentage, make alpha = 0.75 if there is overlap
  scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks = c(0.5,1,10,30,50)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Sample")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(colour = "black", size = 12, face = "italic"), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 11, face = "bold"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey95"),
        strip.text.x = element_text(face="bold"),
        strip.background = element_rect(color="black", size=1)) +
  guides(fill = "none") +
  facet_wrap(~Insect, scales="free_x", nrow=1) +
  scale_y_discrete(limits=rev) # makes y axis in alphabetical from top down 
plot_bubble
pdf('Sed16S.phylum.bubble.pdf', width=10, height=4)
plot_bubble
dev.off()

df <- ps.rel.genus.glom.noNA %>%
  group_by(Insect, Phylum) %>%
  summarize(mean = mean(Abundance, na.rm = TRUE),
            median = median(Abundance, na.rm = TRUE),
            stdev = sd(Abundance, na.rm = TRUE),
            counts = n(),
            se = stdev / sqrt(counts))


kable(df, digits = 3, align = "c",
      caption = "test") %>%
  kable_classic(full_width = FALSE)

df.w <- df %>%
  mutate(stat = paste0(round(100 * mean, 1), " +/- ", round(100 * stdev, 1))) %>%
  select(Insect, Phylum, stat) %>%
  pivot_wider(
    names_from = Insect,
    values_from = stat
  )
kable(df.w, digits = 3, align = "c",
      caption = "test") %>%
  kable_classic(full_width = FALSE)
