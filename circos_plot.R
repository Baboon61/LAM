#!/usr/bin/env Rscript

### A script allowing to create a Circos plot according to 3 files (histo, link and cytoband file) output from generate_visualization_files.py

##############################IMPORTS##############################

library(getopt, quietly = TRUE)
library(RCircos, quietly = TRUE)

##############################OPTIONS##############################

spec = matrix(c(
'help' , 'h', 0, "logical", "display help",
'file_metadata' , 'm', 1, "character", "metadata file",
'genome' , 'g', 1, "character", "only filter librairies results with this genome",
'output_mark' , 'o', 1, "character", "mark added to the output file",
'dir_post' , 'p', 1, "character", "postprocess directory",
'file_cytoband' , 'y', 1, "character", "the cytoband file corresponding to your genome",
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
if ( is.null(opt$input_mark ) ) { write("Warning : You will process the raw file !",stderr()); opt$input_mark = "" }


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
metadata <- df_metadata[df_metadata$Assembly == opt$genome,]
if (nrow(metadata)==0){
	write("Error : This assembly does not exist in the metadata file !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# CHECK OUTPUT MARK
if (opt$output_mark==""){
	print("Error : You have to set an output mark !\n")
	usage()
	sys.exit(2)
}

# CHECK POSTPROCESS DIRECTORY
if (file.exists(opt$dir_post)){
	if (substrRight(opt$dir_post, 1) != "/"){
		opt$dir_post <- paste0(opt$dir_post, '/')
	}
} else{
	write("Error : The post process directory can not be found !",stderr())
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

# CHECK INPUT MARK HISTORY
if (opt$input_mark == ""){
	write("Warning : You will process the raw file !\n")
}

# SELECT INPUT FILES
if (is.null(opt$input_mark )){
	file_input_extension <- ".csv"
} else{
	file_input_extension <- paste0("_",paste0(gsub(",", "_", opt$input_mark),".csv"))
}

# TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
check_input_mark <- "False"

for(i in 1:nrow(metadata[1])){	
	if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension), collapse='')) & file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension), collapse=''))){
		check_input_mark = "True"
	}
}
if (check_input_mark == "False"){
	write("Error : Your input marks can not localize a good input file !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}
# SELECT OUTPUT FILES
if (is.null(file_input_extension )){
	file_output_extension <- paste(c("_",opt$output_mark,".pdf"), collapse='')
} else{
	file_output_extension <- paste(c(substr(file_input_extension, 1, nchar(file_input_extension)-4),"_",opt$output_mark,".pdf"), collapse='')
}

##############################PRINTS##############################

write('\n-----------------------------------------',stderr())
write(paste0('Metadata file : ',opt$file_metadata),stderr())
write(paste0('Genome : ',opt$genome),stderr())
write(paste0('Postprocess directory : ',opt$dir_post),stderr())
write(paste0('Cytoband file : ',opt$file_cytoband),stderr())
write(paste0('Input file extension: ',file_input_extension),stderr())
write(paste0('Output file extension : ',file_output_extension),stderr())
write('-----------------------------------------\n',stderr())

##############################PROGRAM##############################

for(i in 1:nrow(metadata[1])){
	# CHECK POSTPROCESS DIRECTORY EXISTS
	if (!file.exists(paste0(opt$dir_post,as.vector(metadata$Library[i])))){
		write(paste(c("Warning : ",opt$dir_post," does not contains {",as.vector(metadata$Library[i]),"}"), collapse=''),stderr())
		write(paste(c("Warning :  {",i,"} will not be filtered"),collapse=''),stderr())
	} else{
		# CHECK INPUT FILE EXISTS
		if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension), collapse='')) & file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension), collapse=''))){
			write(as.vector(metadata$Library[i]),stderr())

			# IMPORT CSV
			histo <- read.table(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension), collapse=''), header = TRUE)
			link <- read.table(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension), collapse=''), header = TRUE)
			cytoband <- read.table(opt$file_cytoband, header = TRUE)
			colnames(cytoband) <- c("Chromosome","Start","End","Band","Stain")

			# INITIALIZE CIRCOS PARAMETERS
			RCircos.Set.Core.Components(cyto.info=cytoband, chr.exclude=NULL, tracks.inside=3, tracks.outside=0);

			# SET UP OUTPUT FILE
			write(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),file_output_extension), collapse=''),stderr())
			pdf(file=paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),file_output_extension), collapse=''), height=8, width=8);

			# CALL AN AREA FOR CIRCOS PLOT
			RCircos.Set.Plot.Area();

			# WRITE A TITLE
			title(paste(c("Circos plot library ",as.vector(metadata$Library[i])), collapse=' '));

			# DRAW CHROMOSOME IDEOGRAM
			RCircos.Chromosome.Ideogram.Plot();

			# HISTOGRAM PLOT
			RCircos.Histogram.Plot(hist.data=histo, data.col=4, track.num=1, side="in", min.value=0, max.value=100);

			# LINK LINES. LINK DATA HAS ONLY PAIRED CHROMOSOME LOCATIONS IN EACH ROW AND LINK LINES ARE ALWAYS DRAWN INSIDE OF CHROMOSOME IDEOGRAM
			RCircos.Link.Plot(link.data=link, track.num=2, by.chromosome=TRUE);

			# ADD RIBBON LINK TO THE CENTER OF THE PLOT AREA
			RCircos.Ribbon.Plot(ribbon.data=link, track.num=3, by.chromosome=TRUE);

			#  CLOSE THE GRAPHIC DEVICE AND CLEAR MEMORY
			dev.off();

		} else{
			write(paste(c("Warning : ",opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Histo",file_input_extension," does not exist or ",opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Link",file_input_extension," does not exist"), collapse=''),stderr())
			write(paste(c("Warning :  {",opt$dir_post,as.vector(metadata$Library[i]),"} will not be filtered"), collapse=''),stderr())
		}
	}
}
