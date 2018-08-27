#!/usr/bin/env Rscript

### A script allowing to create a Circos plot according to 2 files (Histo and link) output from generate_visualization_files.py
##########################################################################IMPORTS##########################################################################

library(getopt, quietly = TRUE)
library(RCircos, quietly = TRUE)

###################################################################################OPTIONS###################################################################################

spec = matrix(c(
'help' , 'h', 0, "logical", "display help",
'metadata' , 'f', 1, "character", "metadata file",
'pos_directory' , 'k', 1, "character", "postprocess directory",
'outputMark' , 'm', 1, "character", "mark added to the output file",
'genome' , 'g', 1, "character", "only filter librairies results with this genome",
'cytoband' , 'c', 1, "character", "the cytoband file corresponding to your genome",
'inputMark' , 'n', 2, "character", "marks from input file"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if ( !is.null(opt$help) ) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

##################################################################################FUNCTIONS##################################################################################

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

###################################################################################CHECK UP/SET UP###################################################################################

if ( is.null(opt$metadata ) ) { write("Error : -f|--metadata option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$pos_directory ) ) { write("Error : -k|--pos_directory option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$outputMark ) ) { write("Error : -m|--outputMark option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$genome ) ) { write("Error : -g|--genome option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$cytoband ) ) { write("Error : -c|--cytoband option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$inputMark ) ) { write("Warning : You will process the raw file !",stderr()); opt$inputMark = "" }

###CHECK METADATA FILE
if (file.exists(opt$metadata)){
	df_metadata = read.csv(file=opt$metadata, sep="\t", header=T)
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
###CHECK POSTPROCESS DIRECTORY
if (file.exists(opt$pos_directory)){
	if (substrRight(opt$pos_directory, 1) != "/"){
		opt$pos_directory <- paste0(opt$pos_directory, '/')
	}
} else{
	write("Error : The post process directory can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
###FILTER METADATA FILE IF GENOME INPUT
metadata <- df_metadata[df_metadata$Assembly == opt$genome,]
if (nrow(metadata)==0){
	write("Error : This assembly does not exist in the metadata file !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

###CHECK CYTOBAND FILE
if (!(file.exists(opt$cytoband))){
	write("Error : The cytoband file can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

###CHECK INPUT MARKS HISTORY
if (opt$inputMark == ""){
	write("Warning : You will process the raw file !\n")
}

###CHECK MARK
if (opt$outputMark==""){
	print("Error : You have to set a mark !\n")
	usage()
	sys.exit(2)
}

###SELECT INPUT FILES
if (is.null(opt$inputMark )){
	file_input_extension <- ".csv"
} else{
	file_input_extension <- paste0("_",paste0(gsub(",", "_", opt$inputMark),".csv"))
}

###TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
check_inputMark <- "False"

for(i in 1:nrow(metadata[1])){	
	if (file.exists(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension), collapse='')) & file.exists(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension), collapse=''))){
		check_inputMark = "True"
	}
}
if (check_inputMark == "False"){
	write("Error : Your input marks can not localize a good input file !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
###SELECT OUTPUT FILES
if (is.null(file_input_extension )){
	file_output_extension <- paste(c("_",opt$outputMark,".pdf"), collapse='')
} else{
	file_output_extension <- paste(c(substr(file_input_extension, 1, nchar(file_input_extension)-4),"_",opt$outputMark,".pdf"), collapse='')
}

###################################################################################PRINTS###################################################################################

write('\n-----------------------------------------',stderr())
write(paste0('Metadata File tlx : ',opt$metadata),stderr())
write(paste0('Postprocess Directory : ',opt$pos_directory),stderr())
write(paste0('Genome : ',opt$genome),stderr())
write(paste0('Cytoband : ',opt$cytoband),stderr())
write(paste0('Input file extension: ',file_input_extension),stderr())
write(paste0('Output file extension : ',file_output_extension),stderr())
write('-----------------------------------------\n',stderr())

###################################################################################PROGRAM###################################################################################

for(i in 1:nrow(metadata[1])){
###CHECK DIRECTORY EXISTS
	if (!file.exists(paste0(opt$pos_directory,as.vector(metadata$Library[i])))){
		write(paste(c("Warning : ",opt$pos_directory," does not contains {",as.vector(metadata$Library[i]),"}"), collapse=''),stderr())
		write(paste(c("Warning :  {",i,"} will not be filtered"),collapse=''),stderr())
	} else{
		###CHECK INPUT FILE EXISTS
		if (file.exists(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension), collapse='')) & file.exists(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension), collapse=''))){
			write(as.vector(metadata$Library[i]),stderr())
			#Import CSV
			histo <- read.table(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension), collapse=''), header = TRUE)
			link <- read.table(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension), collapse=''), header = TRUE)
			cytoband <- read.table(opt$cytoband, header = TRUE)
			colnames(cytoband) <- c("Chromosome","Start","End","Band","Stain")

			RCircos.Set.Core.Components(cyto.info=cytoband, chr.exclude=NULL, tracks.inside=3, tracks.outside=0);

			write(paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),substr(file_input_extension,1,nchar(file_input_extension)-4),"_",opt$outputMark,".pdf"), collapse=''),stderr())
			pdf(file=paste(c(opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),substr(file_input_extension,1,nchar(file_input_extension)-4),"_",opt$outputMark,".pdf"), collapse=''), height=8, width=8);

			#Call an area for the plot
			RCircos.Set.Plot.Area();

			#Write a title
			title(paste(c("Circos plot library ",as.vector(metadata$Library[i])), collapse=' '));

			#   Draw chromosome ideogram
			RCircos.Chromosome.Ideogram.Plot();

			#################################

			#   Histogram plot
			RCircos.Histogram.Plot(hist.data=histo, data.col=4, track.num=1, side="in", min.value=0, max.value=100);

			#   Link lines. Link data has only paired chromosome locations in each row and link lines are always drawn inside of chromosome ideogram.
			RCircos.Link.Plot(link.data=link, track.num=2, by.chromosome=TRUE);

			#   Add ribbon link to the center of plot area (link lines).
			RCircos.Ribbon.Plot(ribbon.data=link, track.num=3, by.chromosome=TRUE);

			#################################

			#   Close the graphic device and clear memory
			dev.off();
		} else{
			write(paste(c("Warning : ",opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension," does not exist or ",opt$pos_directory,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension," does not exist"), collapse=''),stderr())
			write(paste(c("Warning :  {",opt$pos_directory,as.vector(metadata$Library[i]),"} will not be filtered"), collapse=''),stderr())
		}
	}
}
