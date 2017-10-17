library(dplyr)
library(ggplot2)
library(reshape2)
library(plyr)

#set working directory
setwd("/Users/falk/BCB546X/BCB546X-Fall2017/R_Assignment")

####################### Part 1 ##################################

#read in fang_genotypes text file and set it as genotypes
genotypes <- read.table("fang_et_al_genotypes.txt", header = T, stringsAsFactors = T)

#take a look at the data
head(genotypes)
str(genotypes)
nrow(genotypes)
ncol(genotypes)
geno_table <- tbl_df(genotypes)
head(geno_table)
glimpse(geno_table)
data.frame(head(geno_table))

#read in snp_position text file and set it as snps
snps <- read.delim("snp_position.txt", header = T, stringsAsFactors = F)

#take a look at the data
head(snps)
str(genotypes)
nrow(snps)
ncol(snps)
snps_table <- tbl_df(snps)
head(snps_table)
glimpse(snps_table)

#create a list of the SNP Names for later use
SNP_names <- colnames(genotypes)[-c(1:3)]

####################################################################
############################# MAIZE ################################

#use the filter verb to screen for just the maize SNPs
maize <- filter(genotypes, Group %in% c("ZMMIL","ZMMLR","ZMMMR"))
#lets take a look at our new df
nrow(maize)
head(maize)

#transpose maize df so that it reads horizontally
maize_t <- t(maize)

#use merge function to combine both files together using the SNP_ID column
#in the snps_table and the dataframe's actual row names in the maize_t df.
maize_m <- merge(snps_table, maize_t, by.x = "SNP_ID", by.y = "row.names")

#cut out unwanted columns
maize_c <- maize_m[,-c(2,5:15)]

#change the ? to - using the apply function
df <- apply(maize_c,2, function(x) gsub("\\?/\\?", "-/-", x))
#apply outputs a matrix, so change it to a dataframe
genotypes2 <- as.data.frame(df)

#reorder dataframe by increasing position number
maize_increasing_snps <- maize_c[order(as.numeric(as.character(maize_c$Position))),]

#reorder dataframe by decreasing position number
maize_decreasing_snps <- genotypes2[order(-as.numeric(as.character(genotypes2$Position))),]

#instead of manually creating each dataframe, we can use a for loop
for (i in 1:10) {
  x <- maize_increasing_snps[maize_increasing_snps$Chromosome == i,]
#we can create a new CSV file for each chromosome
  write.csv(x, sprintf("maize_increasing_snps_chromosome_%d.csv", i), row.names = F)
}

#instead of manually creating each dataframe, we can use a for loop
for (i in 1:10) {
  x <- maize_decreasing_snps[maize_decreasing_snps$Chromosome == i,]
  #we can create a new CSV file for each chromosome
  write.csv(x, sprintf("maize_decreasing_snps_chromosome_%d.csv", i), row.names = F)
}


####################################################################
########################### TEOSINTE ###############################

#use the filter verb to screen for just the teosinte SNPs
teosinte <- filter(genotypes, Group %in% c("ZMPBA","ZMPIL","ZMPJA"))
#lets take a look at our new df
head(teosinte)
nrow(teosinte)
ncol(teosinte)

#transpose teosinte df so that it reads horizontally
teosinte_t <- t(teosinte)

#use merge function to combine both files together using the SNP_ID column
#in the snps_table and the dataframe's actual row names in the teosinte_t df.
teosinte_m <- merge(snps_table, teosinte_t, by.x = "SNP_ID", by.y = "row.names")

#cut out unwanted columns
teosinte_c <- teosinte_m[,-c(2,5:15)]

#change the ? to - using the apply function
df <- apply(teosinte_c,2, function(x) gsub("\\?/\\?", "-/-", x))
#apply outputs a matrix, so change it to a dataframe
genotypes2 <- as.data.frame(df)

#reorder dataframe by increasing position number
teosinte_increasing_snps <- teosinte_c[order(as.numeric(as.character(teosinte_c$Position))),]

#reorder dataframe by decreasing position number
teosinte_decreasing_snps <- genotypes2[order(-as.numeric(as.character(genotypes2$Position))),]

#instead of manually creating each dataframe, we can use a for loop
for (i in 1:10) {
  x <- teosinte_increasing_snps[teosinte_increasing_snps$Chromosome == i,]
  #we can create a new CSV file for each chromosome
  write.csv(x, sprintf("teosinte_increasing_snps_chromosome_%d.csv", i), row.names = F)
}

#instead of manually creating each dataframe, we can use a for loop
for (i in 1:10) {
  x <- teosinte_decreasing_snps[teosinte_decreasing_snps$Chromosome == i,]
  #we can create a new CSV file for each chromosome
  write.csv(x, sprintf("teosinte_decreasing_snps_chromosome_%d.csv", i), row.names = F)
}

####################### Part 2 ##################################

#transpose genotypes df so that it reads horizontally
genotypes_t <- t(genotypes)

#use merge function to combine both files together using the SNP_ID column
#in the snps_table and the dataframe's actual row names in the maize_t df.
genotypes_m <- merge(snps, genotypes_t, by.x = "SNP_ID", by.y = "row.names", all.y = TRUE)

#Pull the name of all the SNPs from the original file
SNP_names <- colnames(genotypes)[-c(1:3)]
#reshape the data using the melt command
alleles <- melt(genotypes, measure.vars = SNP_names)

#look at the just the OTU and Groups and SNPs
genotypes_n <- genotypes_m[,-c(1:15)]
#Make the OTU Groups the header
colnames(genotypes_n) <- as.character(unlist(genotypes_n[14,]))
#Remove the OTU and Groups rows
genotypes_n <- genotypes_n[-c(14,16),]
View(genotypes_n)

#look at the just the other stuff
genotypes_o <- genotypes_m[,c(1:15)]
#remove the rows of OTU and Groups
genotypes_o <- genotypes_o[-c(14,16),]

#bind both dataframes together in new format
genotypes_p <- cbind(genotypes_o, genotypes_n)
genotypes.p <- melt(genotypes_p)

#reorder dataframe by increasing position number
genotypes_s <- genotypes_m[order(as.numeric(as.character(genotypes_m$Chromosome))),]

genotypes_s <- tbl_df(genotypes_s)
#plot the SNPs per chromosome
ggplot(genotypes_s, aes(x=Chromosome)) + geom_bar()

#######################################################################################

#create a summary table of SNPs per chromosome
x <- genotypes_s %>% group_by(Chromosome) %>% summarise(SNP_count = n())
#turn columns into factors so that we can plot them
x$Chromosome <- as.factor(x$Chromosome)
x$SNP_count <- as.factor(x$SNP_count)
#plot the SNPs per chromosome
plot(x$Chromosome, x$SNP_count, xlab = "Chromosome", ylab = "SNP Count")
#in general, there is a decreasing number of SNPs from the first to tenth chromosome

######################################################################################
#Missing data and amount of heterozygosity

#Create a new column to indicate whether a particular site is homozygous
alleles <- apply(alleles,2, function(x) gsub("C/C", TRUE, x))
alleles <- apply(alleles,2, function(x) gsub("T/T", TRUE, x))
alleles <- apply(alleles,2, function(x) gsub("A/A", TRUE, x))
alleles <- apply(alleles,2, function(x) gsub("G/G", TRUE, x))
alleles <- apply(alleles,2, function(x) gsub("G/T", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("A/C", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("A/G", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("A/T", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("C/G", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("C/T", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("C/A", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("G/A", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("G/C", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("T/A", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("T/C", FALSE, x))
alleles <- apply(alleles,2, function(x) gsub("T/G", FALSE, x))  
#Recode the missing data as NA
alleles <- apply(alleles,2, function(x) gsub("\\?/\\?", NA, x))  

#Name columns in the data set
colnames(alleles_df)[1:5] <- c("Sample_ID","JG_OTU","Group","SNP_ID","Homozygous")
alleles_df <- as.data.frame(alleles)

#Sort your dataframe using Group and Species_ID values
alleles_s <-arrange(alleles_df, Sample_ID, Group)
#Count the number of Homozygous (TRUE) and Heterozygous (FALSE)
alleles_tf <- alleles_s %>% 
  group_by(Homozygous) %>% 
  tally(sort=TRUE)

#Use ggplot to plot the number of Homozygous, Non-homozygous and NA samples
ggplot(data = alleles_tf, aes(x = Homozygous, y = as.numeric(n))) + geom_bar(stat = "identity", position = "stack")

#create dataframe that compares Group with Homozygousity
group_homo <- alleles_s %>%
  group_by(Group) %>%
  select(Homozygous) %>%
  data.frame()

#Visualize one other feature of the dataset.
#Histogram comparing homozygous loci and species grouping
ggplot(data = group_homo, aes(x = Group, y = Homozygous))+ geom_bar(stat = "identity", position = "stack")

