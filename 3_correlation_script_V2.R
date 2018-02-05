#give as arguments the phenotype file and then the folder with the counts
#Define working directory
setwd("./"); 


#the first argumente is the file with the phenotypes
#the second is the folder with the assoc results
args <- commandArgs(TRUE)

#read the phenotype input file
#The input file must have no header row and column.
#the columns xorespond to: bigsid sucretion_levels  phenotypic_group
data <- read.table(args[1], sep="\t");


#save the values into vectors
bigs_ids <- data[,c("V1")]
sucretion_levels <- data[,c("V2")]
phenotypes <- data[,c("V3")]

#get the file names with the numbers of words per gene per isolate
files <- list.files(path=args[2], pattern="*.txt", full.names=T, recursive=FALSE)

#for every associated gene
for (i in 1:length(files)){
  #print(files[i])
  
  #get the gene name
  #gene_name <- strsplit(files[i],"/")
  #print(gene_name[1]) 
  # stop("Message")
  
  
  #read the number of words per gene file
  #words_count <- read.table(files[i], sep="\t");
  #print(files[i]);
  words_count <- try(read.table(files[i], sep="\t"))
  if (inherits(words_count, 'try-error')){  print(files[i]); next;  } 
  #print(words_count); stop("test");
  
  #save the word countes per isolate
  number_of_words <- words_count[,c("V2")]
  
  #print(number_of_words)
  #print(sucretion_levels)
  
  #concatenate the table 
  concatanated_table = cbind(sucretion_levels, number_of_words)
  concatanated_table <- as.data.frame(concatanated_table)
  
  #uncomment if you want to do the percentile analisys (or exclude some isolates!!! Give tehm the value 0.5)
  concatanated_table <- concatanated_table[concatanated_table$sucretion_levels != 0.5,]
  #print(concatanated_table)
  #quit()
  
  #initialise 2by2 table
  a <- 0
  b <- 0
  c <- 0
  d <- 0
  
  #calculate the 2 by 2 matrix 
  for (j in 1:nrow(concatanated_table)) {
	#print(concatanated_table[i, "sucretion_levels"])
	if(concatanated_table[j, "sucretion_levels"] == 1 & concatanated_table[j, "number_of_words"] == 1) {a <- a + 1; next;}
	else if(concatanated_table[j, "sucretion_levels"] == 1 & concatanated_table[j, "number_of_words"] == 0) {b <- b + 1; next;}
	else if(concatanated_table[j, "sucretion_levels"] == 0 & concatanated_table[j, "number_of_words"] == 1) {c <- c + 1; next;}
	else if(concatanated_table[j, "sucretion_levels"] == 0 & concatanated_table[j, "number_of_words"] == 0) {d <- d + 1; next;}
	else{print("Something is wrong!!!!"); quit()}
  }
  
   #form the 2 by 2 table
   two_by_two_table <- t(matrix(c(a, b, c, d), nrow = 2, ncol = 2, dimnames = list(c("Present", "Absent"), c("High", "Low"))))
   
   #perform fisher test
   fisher_test = fisher.test(two_by_two_table, alternative = "two.sided")
   #and store the results
   fisher_pvalue <- fisher_test$p.value
   fisher_log_odds_ratio <- unname(fisher_test$estimate)
   #quit();

  #run the linear regretion
  model = lm(concatanated_table$sucretion_levels ~ concatanated_table$number_of_words)
  #print(summary(model))
  
  
  #calculate the correlation coeficient
  #print(cor(sucretion_levels,number_of_words))
  
  #get the regretion stats  
  adj_r_squared <- summary(model) $adj.r.squared
  r_squared <- summary(model) $r.squared
  pvalue <- anova(model)$'Pr(>F)'
  
  #cat(files[i],"\t",r_squared, "\t", pvalue, "\n")
  
  
  if(fisher_pvalue <= 0.05){
    cat(files[i],"\t",r_squared, "\t", pvalue[1], "\t", fisher_pvalue, "\t", fisher_log_odds_ratio, "\n")
	
    #write.table(concatanated_table, file = "foo.csv")
    #pdf("plots.pdf")
    #plot the correlation
    #plot(number_of_words, sucretion_levels)
    #title(main=files[i], col.main="red", font.main=4)
    #abline(model, col="red")
    #dev.off()
    #quit()	
  }
  
  
 
}

