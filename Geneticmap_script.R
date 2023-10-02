setwd("~/")
library("mappoly")

#working with all dosage calls
dat<- read_geno_csv(file.in  = "draperjewel_all.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)


#working with individual chromosome
#LG1
dat<- read_geno_csv(file.in  = "draperjewel1.csv", ploidy = 4) #loading data
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE) #filtering based on segregation
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores) #calculating recombination fraction
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init) #extracting genomic order
s.o <- make_seq_mappoly(id) 
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3) #plotting recombination fraction based on genomic order

#plotting genetic map
map1<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise,  extend.tail=10)
plot(map1) #genetic map
plot(map1,mrk.names = TRUE)
plot_genome_vs_map(map1) #genomic vs genetic distance map
map1.err<-est_full_hmm_with_global_error(map1,error=0.05) #applying global error to reestimate genetic distance
plot(map1.err) 
plot(map1.err,mrk.names = TRUE)
plot_genome_vs_map(map1.err)

#LG2
dat<- read_geno_csv(file.in  = "draperjewel2.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map2<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise, thres.twopt = 3, extend.tail=15)
plot(map2)
plot(map2,mrk.names = TRUE)
plot_genome_vs_map(map2)
map2.err<-est_full_hmm_with_global_error(map2,error=0.05)
plot(map2.err)
plot(map2.err,mrk.names = TRUE)
plot_genome_vs_map(map2.err)

#LG3
dat<- read_geno_csv(file.in  = "draperjewel4.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map4<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise,  extend.tail=10)
plot(map4)
plot(map4,mrk.names = TRUE)
plot_genome_vs_map(map4)
map4.err<-est_full_hmm_with_global_error(map4,error=0.05)
plot(map4.err)
plot(map4.err,mrk.names = TRUE)
plot_genome_vs_map(map4.err)

#LG4
dat<- read_geno_csv(file.in  = "draperjewel6.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map6<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise, thres.twopt = 4, extend.tail=10)
plot(map6)
plot(map6,mrk.names = TRUE)
plot_genome_vs_map(map6)
map6.err<-est_full_hmm_with_global_error(map6,error=0.05)
plot(map6.err)
plot(map6.err,mrk.names = TRUE)
plot_genome_vs_map(map6.err)

#LG5
dat<- read_geno_csv(file.in  = "draperjewel7.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map7<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise, thres.twopt = 3, extend.tail=20)
plot(map7)
plot(map7,mrk.names = TRUE)
plot_genome_vs_map(map7)
map7.err<-est_full_hmm_with_global_error(map7,error=0.05)
plot(map7.err)
plot(map7.err,mrk.names = TRUE)
plot_genome_vs_map(map7.err)


#LG6
dat<- read_geno_csv(file.in  = "draperjewel11.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map11<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise,  extend.tail=15)
plot(map11)
plot(map11,mrk.names = TRUE)
plot_genome_vs_map(map11)
map1.err<-est_full_hmm_with_global_error(map11,error=0.05)
plot(map11.err)
plot(map11.err,mrk.names = TRUE)
plot_genome_vs_map(map11.err)


#LG7
dat<- read_geno_csv(file.in  = "draperjewel12.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map12<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise, extend.tail=10)
plot(map12)
plot(map12,mrk.names = TRUE)
plot_genome_vs_map(map12)
map12.err<-est_full_hmm_with_global_error(map12,error=0.05)
plot(map12.err)
plot(map12.err,mrk.names = TRUE)
plot_genome_vs_map(map12.err)


#LG8
dat<- read_geno_csv(file.in  = "draperjewel13.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map13<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise, thres.twopt=3, extend.tail=10)
plot(map13)
plot(map13,mrk.names = TRUE)
plot_genome_vs_map(map13)
map13.err<-est_full_hmm_with_global_error(map13,error=0.05)
plot(map13.err)
plot(map13.err,mrk.names = TRUE)
plot_genome_vs_map(map13.err)


#LG9
dat<- read_geno_csv(file.in  = "draperjewel17.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map17<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise,  extend.tail=10)
plot(map17)
plot(map17,mrk.names = TRUE)
plot_genome_vs_map(map17)
map17.err<-est_full_hmm_with_global_error(map17,error=0.05)
plot(map17.err)
plot(map17.err,mrk.names = TRUE)
plot_genome_vs_map(map17.err)


#LG10
dat<- read_geno_csv(file.in  = "draperjewel20.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map20<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise,  extend.tail=10)
plot(map20)
plot(map20,mrk.names = TRUE)
plot_genome_vs_map(map20)
map20.err<-est_full_hmm_with_global_error(map20,error=0.05)
plot(map20.err)
plot(map20.err,mrk.names = TRUE)
plot_genome_vs_map(map20.err)


#LG11
dat<- read_geno_csv(file.in  = "draperjewel21.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map21<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise,  extend.tail=10)
plot(map21)
plot(map21,mrk.names = TRUE)
plot_genome_vs_map(map21)
map21.err<-est_full_hmm_with_global_error(map21,error=0.05)
plot(map21.err)
plot(map21.err,mrk.names = TRUE)
plot_genome_vs_map(map21.err)

#LG12
dat<- read_geno_csv(file.in  = "draperjewel22.csv", ploidy = 4)
print(dat, detailed = TRUE)
pval.bonf <- 0.05/dat$n.mrk
mrks.chi.filt <- filter_segregation(dat, chisq.pval.thres =  pval.bonf, inter = TRUE)
seq.init <- make_seq_mappoly(mrks.chi.filt)
plot(seq.init)
n.cores = parallel::detectCores() - 1
all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = n.cores)
all.rf.pairwise
mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)
id<-get_genomic_order(seq.init)
s.o <- make_seq_mappoly(id)
plot(mat, ord = s.o$seq.mrk.names, fact = 3)
s1<- make_seq_mappoly(seq.init)
plot(mat, ord=s1$seq.mrk.names, fact = 3)
map22<-est_rf_hmm_sequential(s1, twopt = all.rf.pairwise, thres.twopt = 3, extend.tail=10)
plot(map22)
plot(map22,mrk.names = TRUE)
plot_genome_vs_map(map22)
map22.err<-est_full_hmm_with_global_error(map22,error=0.05)
plot(map22.err)
plot(map22.err,mrk.names = TRUE)
plot_genome_vs_map(map22.err)

#constructing combined map
MAP <-vector("list", 12)
MAP[[1]] <- map1.err
MAP[[2]] <- map2.err
MAP[[3]] <- map4.err
MAP[[4]] <- map6.err
MAP[[5]] <- map7.err
MAP[[6]] <- map11.err
MAP[[7]] <- map12.err
MAP[[8]] <- map13.err
MAP[[9]] <- map17.err
MAP[[10]] <- map20.err
MAP[[11]] <- map21.err
MAP[[12]] <- map22.err

#plotting all genetic map
plot_genome_vs_map(MAP, same.ch.lg = T)
plot_map_list(MAP, col=rainbow(12), horiz = FALSE)
knitr::kable(summary_maps(MAP))
export_map_list(MAP, file = "output_file.csv")

#F1 individual haplotype reconstruction
genoprob.err <- vector("list", 12)
for(i in 1:12)
  genoprob.err[[i]] <- calc_genoprob_error(input.map = MAP[[i]], error = 0.05)

homoprobs = calc_homologprob(genoprob.err)
plot(homoprobs, lg = "all", ind = 1) #plotting haplotype 
plot(homoprobs, lg = "all", ind = 175)
