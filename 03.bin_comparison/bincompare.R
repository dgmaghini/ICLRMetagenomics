library(tidyverse)

df <- read.table("/Users/dylanmaghini/scg4/projects/benchmarking/AssemblyMethods/03.bin_comparison/data_tables/Ndb.csv", sep=",", header=TRUE)

# filter out self-comparisons
df <- df %>% filter(ani != 1)

# make a unified list of all the comparisons and remove duplicates (since it compares bin A to bin B, and bin B to bin A)
df <- mutate(df, comp=ifelse(querry > reference, paste(querry, reference), paste(reference, querry)))
df <- df[!duplicated(df[c('comp')]),]

# get the donor for the query and donor for the reference
df <- mutate(df, donor_querry=sub("_.*", "" , sub("nanopore_", "", sub("infinity_", "", querry))))
df <- mutate(df, donor_ref=sub("_.*", "" , sub("nanopore_", "", sub("infinity_", "", reference))))

# get the method for the query and the method for the reference
df <- mutate(df, querry_method=sub("_.*", "", querry))
df <- mutate(df, ref_method=sub("_.*", "", reference))

df <- df %>% filter(donor_querry == donor_ref) %>% filter(querry_method != ref_method) %>% filter(ani > 0.99)

comp_table <- df %>% count("donor_ref")

write.csv(comp_table, "/Users/dylanmaghini/scg4/projects/benchmarking/AssemblyMethods/03.bin_comparison/matching_bins.csv", quote=FALSE, row.names = FALSE)

paths <- df %>% select(querry, reference) %>% 
  mutate(querry_path=paste("/labs/asbhatt/jackshan/projects/Illumina_ONP_Comparision/Dereplication/genomes/", querry, sep="")) %>%
  mutate(ref_path=paste("/labs/asbhatt/jackshan/projects/Illumina_ONP_Comparision/Dereplication/genomes/", reference, sep=""))

names(paths) <- c("Query", "Reference", "Query_Path", "Reference_Path")

write.table(paths, "/Users/dylanmaghini/scg4/projects/benchmarking/AssemblyMethods/03.bin_comparison/matching_bins_paths.tsv", quote=FALSE, row.names = FALSE, sep="\t")


##### COMPARE BINNING STATISTICS BETWEEN BIN MATCHES #####
# plan
# read in the binning statistics tables, label the bins with the appropriate donor name/method so they match the above table
# merge the tables using the above 'comp' value
# calculate fold change/difference in statistics
# make plots

