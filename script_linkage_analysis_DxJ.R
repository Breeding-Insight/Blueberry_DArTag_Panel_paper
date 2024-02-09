setwd("~/")

#install.packages("devtools")
#devtools::install_github("mmollina/mappoly2", dependencies=TRUE)

#loading library
library(mappoly2)

#loading data
data <- read_geno_csv(file.in = "draperjewel_dose_mp2.csv",
                      ploidy.p1 = 4,
                      ploidy.p2 = 4,
                      name.p1 = "Draper",
                      name.p2 = "Jewel")


plot(data)

#Quality control
p.val<-0.05/3502
data <- filter_data(data, mrk.thresh = 0.05, chisq.pval.thresh = p.val, ind.thresh = 0.07)
plot(data)
data <- filter_individuals(data)
plot(data, type = "density")
data

#calculating pairwise rf
data.all <- pairwise_rf(data, mrk.scope = "all", ncpus = 15)
plot(data.all)

#filter for rf values
data.all <- rf_filter(data.all, 
                      thresh.LOD.ph = 5, 
                      thresh.LOD.rf = 5, 
                      thresh.rf = 0.15, 
                      probs = c(0.05, 0.95))

data.all
plot(data.all)

#assigning linkage groups
g <- group(x = data.all, expected.groups = 12, comp.mat = TRUE, inter = FALSE)
g
plot(g)

#making mappoly2 sequence file
s <- make_sequence(g, 
                   lg = list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                   ch = list(12, 4, 6, 11, 20, 2, 22, 7, 1, 13, 21, 17))

#sorting based on genomic order
s <- order_sequence(s, type = "genome")
print(s, type = "genome")
plot_rf_matrix(s, type = "genome", fact = 2)
s <- rf_filter(s, type = "genome", probs = c(0.05, 0.95), diag.markers = 50)
mappoly2:::plot_rf_matrix(s, type = "genome", fact = 2)

s <- pairwise_phasing(s, 
                      type = "genome",
                      thresh.LOD.ph = 3, 
                      thresh.LOD.rf = 3, 
                      thresh.rf = 0.5, 
                      max.search.expansion.p1 = 5, 
                      max.search.expansion.p2 = 5)
print(s, type = "genome")

#mapping using p1
s <- mapping(s, type = "genome", parent = "p1", ncpus = 15)
print(s, type = "genome")

#mapping using p2
s <- mapping(s, type = "genome", parent = "p2", ncpus = 15)
print(s, type = "genome")

#merging two maps
s <- merge_single_parent_maps(s, type = "genome", ncpus = 15, error = 0.05)
plot_map_list(s, type = "genome", parent = "p1p2")

#dropping weird marker
s<-drop_marker( s,
                "VaccDscaff17_037962776",
                lg = 12,
                type = "genome",
                parent = "p1p2",
                reestimate.map.and.haplo = TRUE,
                verbose = TRUE,
                tol = 0.001
)

plot_map_list(s, type = "genome", parent = "p1p2")

#calculating haplotype probability
s <- calc_haplotypes(s, type = "genome", ncpus = 8)

#augmenting more markers (double simplex)
s <- augment_phased_map(s, type = "genome", ncpus = 8)

#plotting final map
plot_map_list(s, type = "genome", parent = "p1p2",col=mp_pal(12),horiz = T)
map_summary(s, type="genome")
plot_genome_vs_map(s, type = "genome", parent="p1p2", same.ch.lg=TRUE)

#recalculating haplotype probability
s <- calc_haplotypes(s, type = "genome", ncpus = 8)

#plotting haplotypes
plot_haplotypes(s, ind = 1, type="genome", parent="p1p2")
plot_haplotypes(s, ind = 2, type="genome", parent="p1p2")

#exporting map position
t<- plot_map_list(s, type = "genome", parent = "p1p2",col=mp_pal(12),horiz = T)
write.csv(t, "map_position.csv")

#exporting dosages
l1<-as.data.frame(s$maps$lg1$genome$p1p2$rf.phase)
l2<-as.data.frame(s$maps$lg2$genome$p1p2$rf.phase)
l3<-as.data.frame(s$maps$lg3$genome$p1p2$rf.phase)
l4<-as.data.frame(s$maps$lg4$genome$p1p2$rf.phase)
l5<-as.data.frame(s$maps$lg5$genome$p1p2$rf.phase)
l6<-as.data.frame(s$maps$lg6$genome$p1p2$rf.phase)
l7<-as.data.frame(s$maps$lg7$genome$p1p2$rf.phase)
l8<-as.data.frame(s$maps$lg8$genome$p1p2$rf.phase)
l9<-as.data.frame(s$maps$lg9$genome$p1p2$rf.phase)
l10<-as.data.frame(s$maps$lg10$genome$p1p2$rf.phase)
l11<-as.data.frame(s$maps$lg11$genome$p1p2$rf.phase)
l12<-as.data.frame(s$maps$lg12$genome$p1p2$rf.phase)

# Extract the first 9 columns from each dataframe
l1_subset <- l1[, 1:8]
l2_subset <- l2[, 1:8]
l3_subset <- l3[, 1:8]
l4_subset <- l4[, 1:8]
l5_subset <- l5[, 1:8]
l6_subset <- l6[, 1:8]
l7_subset <- l7[, 1:8]
l8_subset <- l8[, 1:8]
l9_subset <- l9[, 1:8]
l10_subset <- l10[, 1:8]
l11_subset <- l11[, 1:8]
l12_subset <- l12[, 1:8]

# Combine the first 9 columns vertically
combined_df <- rbind(l1_subset, l2_subset, l3_subset, l4_subset, l5_subset, l6_subset, l7_subset, l8_subset, l9_subset, l10_subset, l11_subset, l12_subset)

# Export the combined dataframe to a CSV file
write.csv(combined_df, "combined_data_subset.csv", row.names = TRUE)

save.image("mp2.Rdata")



