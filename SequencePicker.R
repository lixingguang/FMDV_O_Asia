library(ape)
library(plyr)

args = commandArgs(trailingOnly=TRUE)


seqs <- read.FASTA(args[1])
names <- labels(seqs)
mat.names <- as.matrix(names)
split.names <- cbind(names,t(as.data.frame(strsplit(mat.names, "|", fixed=TRUE), stringsAsFactors=FALSE)))
data <- as.data.frame(split.names, stringsAsFactors=FALSE)
rownames(data) <- NULL

headings <- c("Full.label", "C1", "C2", "C3", "Country", "Region", "Region.group", "Year", "Accession")

colnames(data) <- headings

data$Year <- strtoi(data$Year)
countries <- unique(data$Country)
year.bins <- seq(1995, 2013, by=2)
year.bins <- sample(year.bins)

sample.counts.list <- list()
sampled.years.list <- list()

replicates <- 10
max.sequences.per.bin <- 10

start <- TRUE

selected.names <- vector()

for(country in countries){

	sampled.years <- 0
	
	print(paste("Country: ",country,sep=""))
	
	country.selections <- vector()
	empty <- TRUE
	first.run <- TRUE
	
	country.data <- data[which(data$Country==country),]
	
	total <- nrow(country.data)
	print(paste("Limit: ",total,sep=""))

	
	repeat{
		if(empty){
			print(paste("Starting loop, ",max.sequences.per.bin," to be added", sep=""))
		} else {
			first.run <- FALSE
			print(paste("Restarting loop, ",max.sequences.per.bin-nrow(country.selections)," to be added",sep=""))
		}
							
		year.bins <- sample(year.bins)
		for (bin.start in year.bins){

			print(paste("Years: ",bin.start,"-",bin.start+1,sep=""))
			
			eligible.sequences <- country.data[which(country.data$Year >= bin.start & country.data$Year<= bin.start +1), ]
			
			if(nrow(eligible.sequences)>=1){
				if(first.run){
					sampled.years <- sampled.years + 1
				}
				
				print(paste("Found ",nrow(eligible.sequences)," sequences",sep=""))
				
				sampled.row <- sample(nrow(eligible.sequences), 1)

				print(paste("Added ",eligible.sequences[sampled.row,]$Full.label,". ",sep=""))
				
			
				country.data <- country.data[country.data[,1]!=eligible.sequences[sampled.row,1], ]

			
				country.selections <- rbind(country.selections, eligible.sequences[sampled.row,])
				print(paste("Have now added ",nrow(country.selections)," sequences", sep=""))
				
				empty <- FALSE
				
						
			} else {
				print("Found no sequences")
				
			}		
	
			if(!empty){
				if(nrow(country.selections)>=total){
					break
				}
				
			 	if(nrow(country.selections)>=max.sequences.per.bin){
			 		break
				}

			}
		}
		if(!empty){
			if(nrow(country.selections)>=total){
				print("Stopping because no sequences remain to be added")
				sample.counts.list[[country]] <- total
				sampled.years.list[[country]] <- sampled.years
				break
			}
				
			 if(nrow(country.selections)>=max.sequences.per.bin){
			 	print("Stopping because the limit has been reached")
			 	sample.counts.list[[country]] <- max.sequences.per.bin
			 	sampled.years.list[[country]] <- sampled.years
			 	break
			}
		}
	
	
	}
	
if(!start){
	selected.names <- rbind(selected.names, country.selections)		
} else {
	selected.names <- country.selections
	start <- FALSE
}	
cat("\n")	
}

selected.seqs <- seqs[selected.names$Full.label]


write.dna(selected.seqs, args[2], format='fasta')




