# code to make prevalence/ plots from INRC bioreactor 16S data

# load packages
library(tidyverse)
library(reshape2)
library(vegan)
library(phyloseq)
library(data.table)
library(viridis)
library(RColorBrewer)
theme_set(theme_light())

# load data
setwd("~/Box Sync/INRCAnalysis/")
ps <- readRDS("16_18inrc.rds")
sample_data(ps)$year <- as.factor(sample_data(ps)$year)
taxa_names(ps) <- paste("OTU-", 1:ntaxa(ps), sep="")
otu <- otu_table(ps)

hist(sample_sums(ps), main="Histogram: Read Counts", xlab="Total Reads", 
     +      border="blue", col="green", las=1, breaks=12)
rarecurve(t(otu_table(ps)))
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)

chips_pcoa <- ordinate(
  physeq = ps.rarefied, 
  method = "PCoA", 
  distance = "bray"
)

ps_scale = ps.rarefied
# Calculate bray curtis distance matrix
chips_bray <- phyloseq::distance(ps_scale, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps))

# Adonis test
adonis(chips_bray ~ year + hrt + location, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(chips_bray, sampledf$year)
permutest(beta)

# Plot 
plot_ordination(
  physeq = ps_scale,
  ordination = chips_pcoa,
  color = "year",
  title = "PCoA of Denitrifying Bioreactor Bacterial Communities"
) + 
  scale_color_manual(name = "Sampling Year",
                     values = c("#C8102E", "#8da0cb")
  ) +
  geom_point(aes(color = year), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() +
  stat_ellipse(geom = "polygon", type="t", level = 0.95, alpha=0.03, show.legend = FALSE)

# Plot 
plot_ordination(
  physeq = ps_scale,
  ordination = chips_pcoa,
  color = "hrt",
  title = "PCoA of Denitrifying Bioreactor Bacterial Communities, Rarefied Dataset"
) + 
  scale_color_manual(name = "Hydraulic Retention Time, h",
                     values = c("#C8102E", "#8da0cb", "#E09F3E")
  ) +
  scale_shape_manual(name = "Sampling Year", values = c(16, 17)) +
  geom_point(aes(color = hrt, shape = year), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() +
  stat_ellipse(geom = "polygon", type="t", level = 0.95, alpha=0.03, show.legend = FALSE) 

ps2 = transform_sample_counts(ps, function(x) x / sum(x) )
ps = ps2
# need to make Schulyer's code run
sep = ''

# load functions from PhyloSmith
source("ineternal.R")
# i did rev(graph_colors) in the source code
source("graphs.R")
source("data_wrangling.R")
source("data_analytics.R")

# should either be 'year' or 'hrt' depending on which analysis you want
myvar = 'year'

# create list with phylosmith functions
keepTaxa1 = common_taxa(ps, treatment = myvar, subset = NULL, n = 'all')

# make prevalence dataframe
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps)) %>% 
  rownames_to_column(var="OTU")
prevdf$OTU <- as.factor(prevdf$OTU)
# make list output from common_taxa a dataframe
keepTaxa2 <- data.frame(OTU = keepTaxa1)
prevdf <- mutate(prevdf, core = ifelse(prevdf$OTU %in% keepTaxa2$OTU, TRUE, FALSE))

# now make dataframe with differences in years; could be a function.....
abundance_to_prevalence <- function(dataset, myvariable, myfactor) {
  subset <- subset_samples(dataset, dataset$myvariable == myfactor) 
  df = apply(X = otu_table(subset),
                   MARGIN = ifelse(taxa_are_rows(subset), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
                   FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  df = data.frame(Prevalence = df,
                        TotalAbundance = taxa_sums(subset),
                        tax_table(subset))
  df <- df %>% rownames_to_column(var="OTU")
  df <- mutate(df, myvariable = myfactor)
  
  return(df)
}
  
# actually do the thing i wanted this function to do
samp2016 <- subset_samples(ps, year == 2016) #sort based on year?
prevdf16 = apply(X = otu_table(samp2016),
                 MARGIN = ifelse(taxa_are_rows(samp2016), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf16 = data.frame(Prevalence = prevdf16,
                      TotalAbundance = taxa_sums(samp2016),
                      tax_table(samp2016))
df16 <- prevdf16 %>% rownames_to_column(var="OTU")
df16 <- mutate(df16, year = 2016)


samp2018 <- subset_samples(ps, year == 2018) #sort based on year?
prevdf18 = apply(X = otu_table(samp2018),
                 MARGIN = ifelse(taxa_are_rows(samp2018), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf18 = data.frame(Prevalence = prevdf18,
                      TotalAbundance = taxa_sums(samp2018),
                      tax_table(samp2016))
df18 <- prevdf18 %>% rownames_to_column(var="OTU")
df18 <- mutate(df18, year = 2018)

total <- rbind(df16, df18)
total$OTU <- as.factor(total$OTU)
total$year <- as.factor(total$year)
total <- mutate(total, core = ifelse(total$OTU %in% keepTaxa2$OTU, TRUE, FALSE))

ggplot(total, aes(x=log(TotalAbundance), y=Prevalence / nsamples(samp2016), color=core)) +
  # Include a guess for parameter
  geom_point(size = 1.5, alpha = 0.25) +
  xlab("ln (Relative Abundance)") + 
  ylab("Occupancy (n = 18)") +
  labs(color = "Genera") +
  facet_wrap(~year) + 
  theme_bw() +
  scale_color_manual(labels = c("Other", "Shared Core"), values = c("red", "black")) +
  theme(text = element_text(size=14))+theme(legend.position = "none")
  
  

# should either be 'year' or 'hrt' depending on which analysis you want
myvar = 'hrt'

# create list with phylosmith functions
keepTaxa1 = common_taxa(ps, treatment = myvar, subset = NULL, n = 'all')

# make prevalence dataframe
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps)) %>% 
  rownames_to_column(var="OTU")
prevdf$OTU <- as.factor(prevdf$OTU)
# make list output from common_taxa a dataframe
keepTaxa2 <- data.frame(OTU = keepTaxa1)
prevdf <- mutate(prevdf, core = ifelse(prevdf$OTU %in% keepTaxa2$OTU, TRUE, FALSE))

# now make dataframe with differences in years; could be a function.....
abundance_to_prevalence <- function(dataset, myvariable, myfactor) {
  subset <- subset_samples(dataset, dataset$myvariable == myfactor) 
  df = apply(X = otu_table(subset),
             MARGIN = ifelse(taxa_are_rows(subset), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
             FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  df = data.frame(Prevalence = df,
                  TotalAbundance = taxa_sums(subset),
                  tax_table(subset))
  df <- df %>% rownames_to_column(var="OTU")
  df <- mutate(df, myvariable = myfactor)
  
  return(df)
}

# actually do the thing i wanted this function to do
samp2016 <- subset_samples(ps, hrt == 2) #sort based on year?
prevdf16 = apply(X = otu_table(samp2016),
                 MARGIN = ifelse(taxa_are_rows(samp2016), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf16 = data.frame(Prevalence = prevdf16,
                      TotalAbundance = taxa_sums(samp2016),
                      tax_table(samp2016))
df16 <- prevdf16 %>% rownames_to_column(var="OTU")
df16 <- mutate(df16, hrt = 2)


samp2018 <- subset_samples(ps, hrt == 8) #sort based on year?
prevdf18 = apply(X = otu_table(samp2018),
                 MARGIN = ifelse(taxa_are_rows(samp2018), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf18 = data.frame(Prevalence = prevdf18,
                      TotalAbundance = taxa_sums(samp2018),
                      tax_table(samp2016))
df18 <- prevdf18 %>% rownames_to_column(var="OTU")
df18 <- mutate(df18, hrt = 8)

sampx <- subset_samples(ps, hrt == 16) #sort based on year?
prevdfx = apply(X = otu_table(sampx),
                 MARGIN = ifelse(taxa_are_rows(sampx), yes = 1, no = 2), #applies the sum function to the columns and counts which samples have abundances of that taxa greater than zero
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdfx = data.frame(Prevalence = prevdfx,
                      TotalAbundance = taxa_sums(sampx),
                      tax_table(sampx))
dfx <- prevdfx %>% rownames_to_column(var="OTU")
dfx <- mutate(dfx, hrt = 16)

total <- rbind(df16, df18, dfx)
total$OTU <- as.factor(total$OTU)
total$hrt <- as.factor(total$hrt)
total <- mutate(total, core = ifelse(total$OTU %in% keepTaxa2$OTU, TRUE, FALSE))
total = total[total$TotalAbundance != 0,]
ggplot(total, aes(x=log(TotalAbundance), y=Prevalence / nsamples(samp2016), color=core)) +
  # Include a guess for parameter
  geom_point(size = 1.5, alpha = 0.25) +
  xlab("ln(Relative Abundance)") + 
  ylab("Occupancy (n = 12)") +
  labs(color = "Genera") +
  facet_wrap(~hrt) + 
  theme_bw() +
  scale_color_manual(labels = c("Other", "Shared Core"), values = c("red", "black")) +
  theme(text = element_text(size=14))+theme(legend.position = "none")

