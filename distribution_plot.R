#!/usr/bin/env Rscript

# A script allowing to create a Distribution plot according to Legitimate and Illegitimate files, output from plsit_junctions.py

##############################IMPORTS##############################
write("\nLibraries loading...",stderr());
library(getopt, warn.conflicts = FALSE)
suppressWarnings(suppressMessages(library(chromPlot, warn.conflicts = FALSE)))
library(regioneR, warn.conflicts = FALSE)
##############################OPTIONS##############################

spec = matrix(c(
'help' , 'h', 0, "logical", "display help",
'file_metadata' , 'm', 1, "character", "metadata file",
'genome' , 'g', 1, "character", "only filter librairies results with this genome",
'output_mark' , 'o', 1, "character", "mark added to the output file",
'dir_post' , 'p', 1, "character", "postprocess directory",
'file_cytoband' , 'y', 1, "character", "the cytoband file corresponding to your genome",
'file_gap' , 'w', 1, "character", "the gap file corresponding to your genome",
'input_mark' , 'i', 2, "character", "marks from input file"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

##############################FUNCTIONS##############################

# SUBSTRING FROM RIGHT
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

##############################CHECK UP/SET UP##############################

if ( is.null(opt$file_metadata ) ) { write("Error : -m|--file_metadata option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$genome ) ) { write("Error : -g|--genome option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$output_mark ) ) { write("Error : -o|--output_mark option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$dir_post ) ) { write("Error : -p|--dir_post option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$file_cytoband ) ) { write("Error : -y|--file_cytoband option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$file_gap ) ) { write("Error : -w|--file_gap option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$input_mark)) { write("Warning : You will process the raw file !",stderr()); opt$input_mark = "" }

# CHECK METADATA FILE
if (file.exists(opt$file_metadata)){
	df_metadata = read.csv(file=opt$file_metadata, sep="\t", header=T)
	if (is.data.frame(df_metadata) && nrow(df_metadata)==0){
		write("Error : The metadata file is empty !",stderr())
		write("\n",stderr())
		cat(getopt(spec, usage=TRUE))
		q(status=1)
	}
	
} else{
	write("Error : The metadata file can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# FILTER METADATA FILE IF GENOME INPUT
metadata <- df_metadata[df_metadata$Assembly == opt$genome,,drop=FALSE]
if (nrow(metadata)==0){
	write("Error : This assembly does not exist in the metadata file !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# CHECK POSTPROCESS DIRECTORY
if (file.exists(opt$dir_post)){
	if (substrRight(opt$dir_post, 1) != "/"){
		opt$dir_post <- paste0(opt$dir_post, '/')
	}
} else{
	write("Error : The postprocess directory can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# CHECK CYTOBAND FILE
if (!(file.exists(opt$file_cytoband))){
	write("Error : The cytoband file can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# CHECK GAP FILE
if (!(file.exists(opt$file_gap))){
	write("Error : The gap file can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# CHECK INPUT MARK HISTORY
if (opt$input_mark == ""){
	write("Warning : You will process the raw file !",stderr())
}

# CHECK OUTPUT MARK
if (opt$output_mark==""){
	write("Error : You have to set an output mark !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# SELECT INPUT FILES
if (is.null(opt$input_mark )){
	file_input_extension <- ".csv"
} else{
	file_input_extension <- paste0("_",paste0(gsub(",", "_", opt$input_mark),".csv"))
}

# TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
check_input_mark <- "False"

for(i in 1:nrow(metadata[,1,drop=FALSE])){
	if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Distribution",file_input_extension), collapse=''))){
		check_input_mark = "True"
	}
}
if (check_input_mark == "False"){
	write("Error : Your input mark can not localize a good input file !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# SELECT OUTPUT FILES
if (is.null(file_input_extension )){
	file_output_extension <- paste(c("_Distribution_",opt$output_mark,".pdf"), collapse='')
} else{
	file_output_extension <- paste(c(substr(file_input_extension, 1, nchar(file_input_extension)-4),"_Distribution_",opt$output_mark,".pdf"), collapse='')
}

##############################PRINTS##############################

write('\n-----------------------------------------',stderr())
write(paste0('Metadata file : ',opt$file_metadata),stderr())
write(paste0('Genome : ',opt$genome),stderr())
write(paste0('Reference file : ',opt$file_reference),stderr())
write(paste0('Postprocess directory : ',opt$dir_post),stderr())
write(paste0('Cytoband file : ',opt$file_cytoband),stderr())
write(paste0('Gap file : ',opt$file_gap),stderr())
write(paste0('Input file extension : ',file_input_extension),stderr())
write(paste0('Output file extension : ',file_output_extension),stderr())
write('-----------------------------------------\n',stderr())

##############################PROGRAM##############################

for(library in 1:nrow(metadata[,1,drop=FALSE])){
# TESTING ONE LIBRARY
#for(library in 1:1){
	# CHECK POSTPROCESS DIRECTORY EXISTS
	if (!file.exists(paste0(opt$dir_post,as.vector(metadata$Library[library])))){
		write(paste(c("Warning : ",opt$dir_post," does not contains {",as.vector(metadata$Library[library]),"}"), collapse=''),stderr())
		write(paste(c("Warning :  {",library,"} will not be filtered"),collapse=''),stderr())
	} else{
		# CHECK INPUT FILE EXISTS
		if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]),"_Distribution",file_input_extension), collapse=''))){
			write(as.vector(metadata$Library[library]),stderr())
			
			# IMPORT CSV
			data <- read.table(paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]),"_Distribution",file_input_extension), collapse=''), header = TRUE, sep="\t", stringsAsFactors = FALSE)

			# GENOME DISPLAY
			cytoband = read.table(opt$file_cytoband, header = TRUE, sep="\t", stringsAsFactors = FALSE)
			gap = read.table(opt$file_gap, header = TRUE, sep="\t", stringsAsFactors = FALSE)

			# PREPARE CYTOBAND
			colnames(cytoband)[1] <- "Chrom"
			colnames(cytoband)[4] <- "Name"
			colnames(cytoband)[5] <- "gieStain"

			# PREPARE GAP
			colnames(gap)[1] <- "Chrom"
			colnames(gap)[4] <- "Name"

			# PREPARE DATA FORWARD/REVERSE
			data_forward <- c()
			data_reverse <- c()
			strand_array <- c()

			for (i in 1:length(data[,1])){
				if (strtoi(data[i,"Strand"]) == 1){
					data_forward <- rbind(data_forward, data[i,])
				}
				else{
					data_reverse <- rbind(data_reverse, data[i,])
				}
				strand_array <- append(strand_array, strtoi(data[i,"Strand"]))
			}

			colnames(data_forward)[1] <- "Chrom"
			colnames(data_reverse)[1] <- "Chrom"

			colnames(data)[1] <- "Chrom"

			#strand_array <- as.vector(data[,"Strand"])
			# OUTPUT LINK PDF
			
			# CREATE PDF NAME
			pdf(file=paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]), file_output_extension), collapse=''), height=8, width=8);

			if (length(data_forward[,1]) >= length(data_reverse[,1])){
				chromPlot(gaps=gap, bands=cytoband, annot1=data_forward, annot2=data_reverse, chrSide=c(1,-1,1,1,1,1,1,1), figCols=11, cex=0.67, colAnnot1="green", colAnnot2="red", title="Distribution plot of double stranded break areas")
			} else{
				chromPlot(gaps=gap, bands=cytoband, annot1=data_reverse, annot2=data_forward, chrSide=c(-1,1,1,1,1,1,1,1), figCols=11, cex=0.67, colAnnot1="red", colAnnot2="green", title="Distribution plot of double stranded break areas")
			}
			# CLOSE THE GRAPHIC DEVICE AND CLEAR MEMORY
			dev.off()
		}
	}
}
