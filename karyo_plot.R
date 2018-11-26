#!/usr/bin/env Rscript

# A script allowing to create a Karyoplot plot according to Karyoplot_link file, output from generate_visualization_files.py

##############################IMPORTS##############################
write("\nLibraries loading...",stderr());
library(getopt, warn.conflicts = FALSE)
library(RCircos, warn.conflicts = FALSE)
suppressWarnings(suppressMessages(library(karyoploteR)))
library(Cairo, warn.conflicts = FALSE)
library(biomaRt, warn.conflicts = FALSE)
library(regioneR, warn.conflicts = FALSE)
suppressWarnings(suppressMessages(library(Biostrings)))
library(stringr, warn.conflicts = FALSE)

###Infos
#https://raw.githubusercontent.com/bernatgel/karyoploter_examples/master/Examples/Tutorial/CustomGenomes/mycytobands.txt
#http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

##############################OPTIONS##############################

spec = matrix(c(
'help' , 'h', 0, "logical", "display help",
'file_metadata' , 'm', 1, "character", "metadata file",
'genome' , 'g', 1, "character", "only filter librairies results with this genome",
'output_mark' , 'o', 1, "character", "mark added to the output file",
'file_reference' , 'r', 1, "character", "the reference file used to call junctions",
'dir_post' , 'p', 1, "character", "postprocess directory",
'file_cytoband' , 'y', 1, "character", "the cytoband file corresponding to your genome",
'file_chrom_length', 'n', 1, "character", "the chromosome length file corresponding to your genome",
'visualization' , 'v', 1, "character", "visualization type for your karyo plot (full_genome, selected_chromosomes or zoom_in)",
'chr_bait' , 'a', 1, "character", "chromosomes bait to display, separated by comma, write 'all' to display all chromosomes",
'chr_prey' , 'b', 1, "character", "chromosomes prey to display, separated by comma, write 'all' to display all chromosomes",
'input_mark' , 'i', 2, "character", "marks from input file",
'UCSC', 'z', 0, "logical", "add UCSC gene loci",
'file_construction', 'c', 2, "character", "construction fasta file for modified genome (in $BOWTIE2_INDEXES)",
'file_locus', 'l', 2, "character", "file to set up some locus labels",
'greek', 'k', 0, "logical", "(only if -l option), transform locus name to greek name (alpha, Alpha, beta...)",
'unlink' , 'u', 0, "logical", "unlink the bait1/prey1, bait2/prey2 process and allows to search for bait1/prey2 junctions",
'file_rename', 'x', 2, "character", "file to rename chromosome to a better display"
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
# CREATE LINKS AND REGIONS ON KARYOPLOT
LinksAndRegions <- function(kp, df, value, arch_circ, r0_Links, r1_Links, r0_Regions, r1_Regions){
    total_bin <- GRanges()
    for (i in 1:length(value_color)){
        if (i==1){
            sub_data <- df[df$value<=value_color[i], ,drop=FALSE]
        }
        else{
            sub_data <- df[df$value>value_color[i-1] & df$value<=value_color[i], ,drop=FALSE]
        }
        if (nrow(sub_data)!=0){
            starts <- toGRanges(data.frame(sub_data$B_Rname,sub_data$B_Rstart,sub_data$B_Rend))
            ends <- toGRanges(data.frame(sub_data$Rname,sub_data$Rstart,sub_data$Rend))
            kpPlotLinks(kp, data=starts, data2=ends, col=names(value_color[i]), border=names(value_color[i]), arch.height=arch_circ, r0=r0_Links, r1=r1_Links)
        }
    }
    for (i in 1:length(value_color)){
        if (i==1){
            sub_data <- df[df$value<=value_color[i], ,drop=FALSE]
        }
        else{
            sub_data <- df[df$value>value_color[i-1] & df$value<=value_color[i], ,drop=FALSE]
        }
        if (nrow(sub_data)!=0){
            starts <- toGRanges(data.frame(sub_data$B_Rname,sub_data$B_Rstart,sub_data$B_Rend))
            ends <- toGRanges(data.frame(sub_data$Rname,sub_data$Rstart,sub_data$Rend))
            suppressWarnings(total_bin <- append(total_bin, ends))
            kpPlotRegions(kp, ends, r0=r0_Regions, r1=r1_Regions, col="#ff8d92", avoid.overlapping=FALSE)
            kpPlotRegions(kp, starts, r0=r0_Regions, r1=r1_Regions, col="#57d53b", avoid.overlapping=FALSE)
        }
    }
    return(total_bin)
}

# CREATE A GENE LIST USING ENSEMBL
RetrieveGenes <- function(total_bin){
    ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
    genes <- GRanges()
    geneList <- c()
    GmList <- c()
    BMresult <- data.frame()
    for (i in 1:length(total_bin)){
        BMresult <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'), filters=c('chromosome_name', 'start', 'end'), values=list(substr(seqnames(total_bin)[i], 4, 10),start(total_bin)[i],end(total_bin)[i]), mart = ensembl)
        geneList <- append(geneList, BMresult[['external_gene_name']])
        suppressWarnings(genes <- append(genes, toGRanges(BMresult)))
    }
    if (length(geneList) == 0){
        return(genes)
    }
    else{
        for (i in 1:length(geneList)){
            if (grepl('Gm*',geneList[i])){
                GmList <- append(GmList,geneList[i])
            }
        }
        genes <- genes[!(elementMetadata(genes)[, "external_gene_name"] %in% GmList)]
        return(genes)
    }
}

# REVERSE GENOME MODIFY POSITIONS
Modify_genome_position <- function(df_construction, chrom, start, end){
	aux_start <- 0
	aux_end <- 0
	aux_df_construction <- df_construction[df_construction$chr == chrom, ]
	if (nrow(aux_df_construction) != 0){
		for (i in 1:nrow(aux_df_construction)){
			if (!(as.integer(aux_df_construction$start[i]) > as.integer(end))){
				if ((as.integer(start) > as.integer(aux_df_construction$start[i])) && (as.integer(end) < as.integer(aux_df_construction$end[i]))){
					return(c(-1,-1))
				}
				aux_start = as.integer(aux_start) + as.integer(aux_df_construction$length[i])
				aux_end = as.integer(aux_end) + as.integer(aux_df_construction$length[i])
			}
		}
	}
	return(c(as.integer(start)-as.integer(aux_start), as.integer(end)-as.integer(aux_end)))
}

# MODIFY GENES POSITIONS (WITH A RANGE)
Modify_genes_position <- function(df_construction, chrom, start, end){
	aux_start <- 0
	aux_end <- 0
	aux_df_construction <- df_construction[df_construction$chr == chrom, ]
	if (nrow(aux_df_construction) != 0){
		for (i in 1:nrow(aux_df_construction)){
			if (!(as.integer(end) < as.integer(aux_df_construction$start[i]))){
				if (as.integer(start) <= as.integer(aux_df_construction$start[i])){
					if (as.integer(end) < as.integer(aux_df_construction$end[i])){
						aux_end = as.integer(aux_end) + (as.integer(end)-as.integer(aux_df_construction$start[i]))
					}else{
						aux_end = as.integer(aux_end) + as.integer(aux_df_construction$length[i])
					}
				}else if ((as.integer(start) > as.integer(aux_df_construction$start[i])) && (as.integer(start) < as.integer(aux_df_construction$end[i]))){
					if (as.integer(end) >= as.integer(aux_df_construction$start[i])){
						aux_start = as.integer(aux_start) + (as.integer(aux_df_construction$end[i]) - as.integer(start))
						aux_end = as.integer(aux_end) + as.integer(aux_df_construction$length[i])
					}
				}else{
					aux_start = as.integer(aux_start) + as.integer(aux_df_construction$length[i])
					aux_end = as.integer(aux_end) + as.integer(aux_df_construction$length[i])
				}
			}
		}
	}
	return(c(as.integer(start)-as.integer(aux_start), as.integer(end)-as.integer(aux_end)))
}


##############################CHECK UP/SET UP##############################

if ( is.null(opt$file_metadata ) ) { write("Error : -m|--file_metadata option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$genome ) ) { write("Error : -g|--genome option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$output_mark ) ) { write("Error : -o|--output_mark option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$file_reference ) ) { write("Error : -r|--file_reference option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$dir_post ) ) { write("Error : -p|--dir_post option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$file_cytoband ) ) { write("Error : -y|--file_cytoband option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$file_chrom_length ) ) { write("Error : -n|--file_chrom_length option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$visualization ) ) { write("Error : -v|--visualization option can not be null (full_genome, selected_chromosomes or zoom_in)",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$chr_bait ) ) { write("Error : -a|--chr_bait option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$chr_prey ) ) { write("Error : -b|--chr_prey option can not be null",stderr()); write("\n",stderr()); cat(getopt(spec, usage=TRUE)); q(status=1) }
if ( is.null(opt$input_mark)) { write("Warning : You will process the raw file !",stderr()); opt$input_mark = "" }
if ( is.null(opt$UCSC ) ) { write("Warning : Do not add UCSC gene loci !",stderr()); opt$UCSC=FALSE }
if ( is.null(opt$file_construction ) ) { write("Warning : You will not add a construction file !",stderr()); opt$file_construction = NULL }
if ( is.null(opt$file_locus ) ) { write("Warning : You will not add personal loci to the karyo plot !",stderr()); opt$file_locus = NULL; opt$greek=FALSE }
if ( is.null(opt$greek ) ) { write("Warning : Do not change loci name with greek letters !",stderr()); opt$greek=FALSE }
if ( is.null(opt$unlink ) ) { write("Warning : Bait and prey are linked !",stderr()); unlink=FALSE }
if ( is.null(opt$file_rename ) ) { write("Warning : You will not rename some chromosome names !",stderr()); opt$file_rename = NULL }


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

# CHECK REFERENCE FILE
if (!(file.exists(opt$file_reference))){
	write("Error : The reference file can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
} else{
	chrom_fasta_list <- c('all')
	chrom_fasta_length <- c()
	chrom_fasta = readDNAStringSet(opt$file_reference)
	for (i in 1:length(chrom_fasta)){
		chrom_fasta_list[i+1] <- strsplit(names(chrom_fasta), " ")[[i]][1]
		chrom_fasta_length[[strsplit(names(chrom_fasta), " ")[[i]][1]]] <- width(chrom_fasta[i])
	}
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

# CHECK CHROMOSOME LENGTH FILE
if (!(file.exists(opt$file_chrom_length))){
	write("Error : The chromosome length file can not be found !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

# CHECK VISUALIZATION
if (opt$visualization==""){
	write("Error : You have to set a visualization (full_genome, selected_chromosomes or zoom_in) !",stderr())
	write("\n",stderr())
	cat(getopt(spec, usage=TRUE))
	q(status=1)
} else{
	if ((opt$visualization!="full_genome") & (opt$visualization!="selected_chromosomes") & (opt$visualization!="zoom_in")){
		write("Error : You have to set a visualization (full_genome, selected_chromosomes or zoom_in) !",stderr())
		write("\n",stderr())
		cat(getopt(spec, usage=TRUE))
		q(status=1)
	}
}

# CHECK CHROMOSOME BAITS OPTION
if (opt$chr_bait != ""){
	chromosomesBait = strsplit(opt$chr_bait, ",")
	matrix_Bait=matrix(,nrow = length(chromosomesBait[[1]]), ncol = 3);
	if (opt$chr_bait == 'all'){
		matrix_Bait[1,] <- c('all', '', '')
	} else{
		for (i in 1:length(chromosomesBait[[1]])){
			name <- strsplit(chromosomesBait[[1]][i], ":")[[1]][1]
			start <- ""
			end <- ""
			if (!('all' == name)){
				if ((length(strsplit(chromosomesBait[[1]][i], ":")[[1]]) > 1) & (length(strsplit(chromosomesBait[[1]][i], ":")[[1]]) <= 3)){
					if (length(strsplit(chromosomesBait[[1]][i], ":")[[1]]) == 2){
						write(paste(c("Error : The chromosome bait ", name," start or end argument missing !"), collapse=''),stderr())
						write("\n",stderr())
						cat(getopt(spec, usage=TRUE))
						q(status=1)
					}
					start <- strsplit(chromosomesBait[[1]][i], ":")[[1]][2]
					end <- strsplit(chromosomesBait[[1]][i], ":")[[1]][3]
					if (strtoi(end) <= strtoi(start)){
						write(paste(c("Error : The chromosome bait ", name," end can not be equal or inferior to start !"), collapse=''),stderr())
						write("\n",stderr())
						cat(getopt(spec, usage=TRUE))
						q(status=1)
					}
				} else if (length(strsplit(chromosomesBait[[1]][i], ":")[[1]]) > 3){
						write(paste(c("Error : The chromosome bait ", name," can not have more than start and stop values !"), collapse=''),stderr())
						write("\n",stderr())
						cat(getopt(spec, usage=TRUE))
						q(status=1)
				}
				if (!(name %in% chrom_fasta_list)){
					write(paste(c("Error : The chromosome bait ", name," can not be found in reference file !"), collapse=''),stderr())
					write("\n",stderr())
					cat(getopt(spec, usage=TRUE))
					q(status=1)
				}
				matrix_Bait[i,] <- c(name, start, end)
			} else{
				matrix_Bait[i,] <- c('all', '', '')
			}
		}
	}
}
#print(matrix_Bait)

# CHECK CHROMOSOME PREYS OPTION
if (opt$chr_prey != ""){
	chromosomesPrey = strsplit(opt$chr_prey, ",")
	matrix_Prey=matrix(,nrow = length(chromosomesPrey[[1]]), ncol = 3);
	if (opt$chr_prey == 'all'){
		matrix_Prey[1,] <- c('all', '', '')
	} else{
		for (i in 1:length(chromosomesPrey[[1]])){
			name <- strsplit(chromosomesPrey[[1]][i], ":")[[1]][1]
			start <- ""
			end <- ""
			if (!('all' == name)){
				if ((length(strsplit(chromosomesPrey[[1]][i], ":")[[1]]) > 1) & (length(strsplit(chromosomesPrey[[1]][i], ":")[[1]]) <= 3)){
					if (length(strsplit(chromosomesPrey[[1]][i], ":")[[1]]) == 2){
						write(paste(c("Error : The chromosome prey ", name," start or end argument missing !"), collapse=''),stderr())
						write("\n",stderr())
						cat(getopt(spec, usage=TRUE))
						q(status=1)
					}
					start <- strsplit(chromosomesPrey[[1]][i], ":")[[1]][2]
					end <- strsplit(chromosomesPrey[[1]][i], ":")[[1]][3]
					if (strtoi(end) <= strtoi(start)){
						write(paste(c("Error : The chromosome prey ", name," end can not be equal or inferior to start !"), collapse=''),stderr())
						write("\n",stderr())
						cat(getopt(spec, usage=TRUE))
						q(status=1)
					}
				} else if (length(strsplit(chromosomesPrey[[1]][i], ":")[[1]]) > 3){
						write(paste(c("Error : The chromosome prey ", name," can not have more than start and stop values !"), collapse=''),stderr())
						write("\n",stderr())
						cat(getopt(spec, usage=TRUE))
						q(status=1)
				}
				if (!(name %in% chrom_fasta_list)){
					write(paste(c("Error : The chromosome prey ", name," can not be found in reference file !"), collapse=''),stderr())
					write("\n",stderr())
					cat(getopt(spec, usage=TRUE))
					q(status=1)
				}
				matrix_Prey[i,] <- c(name, start, end)
			} else{
			matrix_Prey[i,] <- c('all', '', '')
			}
		}
	}
}
#print(matrix_Prey)

# CREATE CONSTRUCTION OBJECT IF UCSC GENES SELECTED AND GENOME IS MODIFIED (NO NEED IF NORMAL GENOME)
if (opt$UCSC) {
	if (!is.null(opt$file_construction)) {
		if (!(file.exists(opt$file_construction))){
			write("Error : The construction file can not be found !",stderr())
			write("\n",stderr())
			cat(getopt(spec, usage=TRUE))
			q(status=1)
		} else{
			# CHECK CONSTRUCTION FILE
			df_construction=setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("chr", "modif_type", "length", "start", "end", "seq"))
			removed_list <- c()
			added_dict <- c()
			#names(added_dict) <- c()
			construction_fasta = readDNAStringSet(opt$file_construction)
			# FOR EACH ELEMENT OF CONSTRUCTION FILE
			for (i in 1:length(construction_fasta)){
				modif_type <- ""
				chrom <- ""
				start <- 0
				end <- 0
				check_good <- TRUE
				if (length(strsplit(names(construction_fasta)[i],":")[[1]]) == 2){
					check_good <- FALSE
					if (strsplit(names(construction_fasta)[i],":")[[1]][1] == "add"){
						if (!strsplit(names(construction_fasta)[i],":")[[1]][2] %in% names(added_dict)){
							added_dict[[strsplit(names(construction_fasta)[i],":")[[1]][2]]] <- as.character(construction_fasta[i])
						}
					
					} else if(strsplit(names(construction_fasta)[i],":")[[1]][1] == "remove"){
						removed_list <- c(removed_list, strsplit(names(construction_fasta)[i],":")[[1]][2])
					}
				}else if (length(strsplit(names(construction_fasta)[i],":")[[1]]) == 4){
					modif_type=strsplit(names(construction_fasta)[i],":")[[1]][1]
					chrom=strsplit(names(construction_fasta)[i],":")[[1]][2]
					if (suppressWarnings(!is.na(as.numeric(strsplit(names(construction_fasta)[i],":")[[1]][3])))){
						start=strsplit(names(construction_fasta)[i],":")[[1]][3]
					}else{
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} start need to be integer, please use this pattern modif_type:chrmX:start:end !"),stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"will not be used"),stderr())
					}
					if (suppressWarnings(!is.na(as.numeric(strsplit(names(construction_fasta)[i],":")[[1]][4])))){
						end=strsplit(names(construction_fasta)[i],":")[[1]][4]
					}else{
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} end need to be integer, please use this pattern modif_type:chrmX:start:end !"),stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"will not be used"),stderr())
					}
					if (modif_type != "insertion" && modif_type != "deletion"){
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} does not have the good format, please use modif_type:chrmX:start:end !"),stderr())
						write("Warning : To delete or insert use modif_type:chrmX:start:end !",stderr())
						write("Warning : To add or remove use modif_type:chrmX !",stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
					}
					if (modif_type == "insertion" && as.character(construction_fasta[i])==""){
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} does not have a sequence to insert"),stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
					}
					if (!chrom %in% chrom_fasta_list){
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} chrom is not in the reference genome file, please use this pattern modif_type:chrmX:start:end !"),stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
					}
					if (start < 0 && end < 1){
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} start < 0 or end position < 1 is not allowed, please use this pattern modif_type:chrmX:start:end !"),stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
					}
					if (end < start){
						check_good=FALSE
						write(paste0("Warning : {",names(construction_fasta)[i],"} end < start is not allowed, please use this pattern modif_type:chrmX:start:end !"),stderr())
						write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
					}
				}else{
					check_good=FALSE
					write(paste0("Warning : {",names(construction_fasta)[i],"} does not have the good format"),stderr())
					write("Warning : To delete or insert use modif_type:chrmX:start:end !",stderr())
					write("Warning : To add or remove use modif_type:chrmX !",stderr())
					write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
				}
				if (check_good){
					# FOR INSERTION DO THE FOLLOWING
					if (modif_type == "insertion"){
						if (as.integer(start) == as.integer(end)){
							write(paste0("Warning : {",names(construction_fasta)[i],"} can not insert between ",str(start)," and ",str(end)," !"),stderr())
							write(paste0("Warning : {",names(construction_fasta)[i],"} insertion between ",str(start)," and ",str(end+1)," !"),stderr())
							end = end + 1
						}
						if (as.integer(start) > as.integer(chrom_fasta_length[[chrom]])){
							write(paste0("Warning : {",names(construction_fasta)[i],"} start > chromosome length !"),stderr())
							write(paste0("Warning : {",names(construction_fasta)[i],"} insertion at the end !"),stderr())
							start=chrom_fasta_length[[chrom]]
							end=chrom_fasta_length[[chrom]]+1
						}
						if (as.integer(end) > as.integer(chrom_fasta_length[[chrom]])+1){
							write(paste0("Warning : {",names(construction_fasta)[i],"} end > chromosome length +1 !"),stderr())
							write(paste0("Warning : {",names(construction_fasta)[i],"} insertion at the end !"),stderr())
							start=chrom_fasta_length[[chrom]]
							end=chrom_fasta_length[[chrom]]+1
						}
						df_construction[nrow(df_construction)+1,] <- c(chrom, modif_type, nchar(as.character(construction_fasta)[i])-(as.integer(end)-as.integer(start)-1), as.integer(start), as.integer(end), as.character(construction_fasta)[i])
					}
					# FOR DELETION DO THE FOLLOWING
					else if (modif_type == "deletion"){
						if (as.integer(start) == as.integer(end) || as.integer(start) == as.integer(end)-1){
							write(paste0("Warning : {",names(construction_fasta)[i],"} can not delete between ",start," and ",end,", nothing to delete !"),stderr())
						}else if (as.integer(start) >= as.integer(chrom_fasta_length[[chrom]]) && as.integer(end) > as.integer(chrom_fasta_length[[chrom]])){
							write(paste0("Warning : {",names(construction_fasta)[i],"} can not delete between ",start," and ",end,", nothing to delete !"),stderr())
						}else{
							if (as.integer(start) < as.integer(chrom_fasta_length[[chrom]]) && as.integer(end) > as.integer(chrom_fasta_length[[chrom]])){
								write(paste0("Warning : {",names(construction_fasta)[i],"} end > chromosome length +1 !"),stderr())
								write(paste0("Warning : {",names(construction_fasta)[i],"} deletion from ",start," to the end !"),stderr())
								end=chrom_fasta_length[[chrom]]+1
							}
							df_construction[nrow(df_construction)+1,] <- c(chrom, modif_type, -(as.integer(end)-as.integer(start)-1), as.integer(start), as.integer(end), "")
						}
					}else{
						write(paste0("Warning : {",names(construction_fasta)[i],"} does not have the good format"),stderr())
						write("Warning : To delete or insert use modif_type:chrmX:start:end !")
						write("Warning : To add or remove use modif_type:chrmX !")
						write(paste0("Warning : {",names(construction_fasta)[i],"} will not be used"),stderr())
						sys.exit()
					}
				}	
			}
			# SORT DATAFRAME
			df_construction <- df_construction[with(df_construction, order(chr, start, end)), ]
		}
	}
}

# CHECK LOCUS FILE
if (!is.null(opt$file_locus)){
	if (!(file.exists(opt$file_locus))){
		write("Error : The locus file can not be found !",stderr())
		write("\n",stderr())
		cat(getopt(spec, usage=TRUE))
		q(status=1)
	}
}

# CHECK UNLINK BAIT/PREY
if (is.null(opt$unlink)){
	if ((chromosomesBait[[1]][1] == 'all') | (chromosomesPrey[[1]][1] == 'all')){
		write("Warning : Default unlink if -b all or -p all !",stderr())
		unlink=TRUE
	} else{
		if (length(chromosomesBait[[1]]) != length(chromosomesPrey[[1]])){
			write("Error : The number of element in chromosome bait and chromosome prey have to be the same if you link them !",stderr())
			write("\n",stderr())
			cat(getopt(spec, usage=TRUE))
			q(status=1)
		}
		unlink=FALSE
	}
	
} else{
	unlink=TRUE
}

# CHECK RENAME FILE
if (!is.null(opt$file_rename)){
	if (!(file.exists(opt$file_rename))){
		write("Error : The rename file can not be found !",stderr())
		write("\n",stderr())
		cat(getopt(spec, usage=TRUE))
		q(status=1)
	} else{
		df_rename = read.csv(file=opt$file_rename, sep="\t", header=F, col.names=c('current', 'new'), stringsAsFactors = FALSE)
		tmp_df_rename_chr = c()
		tmp_df_rename_chrpart = c()
		if (!(is.data.frame(df_rename) && nrow(df_rename)==0)){
			for (i in 1:nrow(df_rename)){
				if (grepl(":", df_rename[i,1])){
					tmp_df_rename_chrpart <- rbind(tmp_df_rename_chrpart, df_rename[i,])
				} else{
					tmp_df_rename_chr <- rbind(tmp_df_rename_chr, df_rename[i,])
				}
			}
		}
	}
} else {
	if (opt$visualization=="zoom_in"){
		write("Error : The rename file is requiered for the zoom_in visualization !",stderr())
		write("\n",stderr())
		cat(getopt(spec, usage=TRUE))
		q(status=1)
	}
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
	if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[i]),"/",as.vector(metadata$Library[i]),"_Karyoplot_link",file_input_extension), collapse=''))){
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
	file_output_extension_link <- paste(c("_Link_",opt$output_mark,".pdf"), collapse='')
	file_output_extension_histo <- paste(c("_Histo_",opt$output_mark,".pdf"), collapse='')
} else{
	file_output_extension_link <- paste(c(substr(file_input_extension, 1, nchar(file_input_extension)-4),"_Link_",opt$output_mark,".pdf"), collapse='')
	file_output_extension_histo <- paste(c(substr(file_input_extension, 1, nchar(file_input_extension)-4),"_Histo_",opt$output_mark,".pdf"), collapse='')
}

##############################PRINTS##############################

write('\n-----------------------------------------',stderr())
write(paste0('Metadata file : ',opt$file_metadata),stderr())
write(paste0('Genome : ',opt$genome),stderr())
write(paste0('Reference file : ',opt$file_reference),stderr())
write(paste0('Postprocess directory : ',opt$dir_post),stderr())
write(paste0('Cytoband file : ',opt$file_cytoband),stderr())
write(paste0('Chromosomes length file : ',opt$file_chrom_length),stderr())
write(paste0('Visualization : ',opt$visualization),stderr())
write(paste0('Chromosomes bait : ',opt$chr_bait),stderr())
write(paste0('Chromosomes prey : ',opt$chr_prey),stderr())
if (opt$UCSC) {
	write('UCSC genes : Yes',stderr())
	if (!is.null(opt$file_construction)) {
		write(paste0('Construction file : ',opt$file_construction),stderr())
	}else{
		write(paste0('Construction file : No'),stderr())
	}	
} else{
	write('UCSC genes : No',stderr())
}
if (!is.null(opt$file_locus)) {
	write(paste0('Locus file : ',opt$file_locus),stderr())
}
if (!is.null(opt$file_locus)) {
	if (!is.null(opt$greek)) {
		write('Greek locus : Yes',stderr())
	} else{
		write('Greek locus : No',stderr())
	}
}
if (unlink) {
	write('Unlink : Yes',stderr())
} else{
	write('Unlink : No',stderr())
}
if (!is.null(opt$file_rename)) {
	write(paste0('Rename file : ',opt$file_rename),stderr())
}
write(paste0('Input file extension: ',file_input_extension),stderr())
write(paste0('Output file extension for link : ',file_output_extension_link),stderr())
write(paste0('Output file extension for histo : ',file_output_extension_histo),stderr())
write('-----------------------------------------\n',stderr())

##############################PROGRAM##############################

# INITIALIZE COLOR VALUES
value_color <- c(0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100)
names(value_color) <- c("#e5e5ff", "#ccccff", "#b2b2ff", "#9999ff", "#7f7fff", "#6666ff", "#4c4cff", "#3232ff", "#1919ff", "#0000ff", "#0000cc", "#000099", "#000066")
custom.frequency <- NULL

for(library in 1:nrow(metadata[,1,drop=FALSE])){
# TESTING ONE LIBRARY
#for(library in 1:1){
	# CHECK POSTPROCESS DIRECTORY EXISTS
	if (!file.exists(paste0(opt$dir_post,as.vector(metadata$Library[library])))){
		write(paste(c("Warning : ",opt$dir_post," does not contains {",as.vector(metadata$Library[library]),"}"), collapse=''),stderr())
		write(paste(c("Warning :  {",library,"} will not be filtered"),collapse=''),stderr())
	} else{
		# CHECK INPUT FILE EXISTS
		if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]),"_Karyoplot_link",file_input_extension), collapse=''))){
			write(as.vector(metadata$Library[library]),stderr())
			# IMPORT LOCUS IF ANY
			if (!is.null(opt$file_locus)){
				if (!(file.exists(opt$file_locus))){
					write("Error : The locus file can not be found !",stderr())
					write("\n",stderr())
					cat(getopt(spec, usage=TRUE))
					q(status=1)
				} else{
					df_locus = read.csv(file=opt$file_locus, sep="\t", header=F, col.names=c('chr', 'start', 'end', 'strand', 'name'), stringsAsFactors = FALSE)
		
				}
			}
			# IMPORT CSV
			karyo_data <- read.table(paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]),"_Karyoplot_link",file_input_extension), collapse=''), header = TRUE, sep="\t", stringsAsFactors = FALSE)

			#print("BEFORE FILTER")
			#print(karyo_data)

			# FILTER BY BAIT/PREY
			
			# IF NOT LINKED
			if (unlink){
				tmp_karyo_data <- c()
				duprows <- c()
				# TEST ALL POSSIBILITIES BAIT/PREY
				for (i in 1:length(matrix_Bait[,1])){
					for (j in 1:length(matrix_Prey[,1])){
						#print(matrix_Bait[i,1])
						#print(matrix_Prey[j,1])
						# IF ALL IN BAIT
						if ('all' == matrix_Bait[i,1]){
							#print("Bait all")
							if ('all' != matrix_Prey[i,1]){
								if (matrix_Prey[j,2] == ""){
									duprows <- rownames(karyo_data[(karyo_data$Rname==matrix_Prey[j,1]),]) %in% rownames(tmp_karyo_data)
									tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$Rname==matrix_Prey[j,1]),][!duprows,])
								} else{
									duprows <- rownames(karyo_data[(karyo_data$Rname==matrix_Prey[j,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[j,2])) & (karyo_data$Junction <= matrix_Prey[j,3]),]) %in% rownames(tmp_karyo_data)
									tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$Rname==matrix_Prey[j,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[j,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[j,3])),][!duprows,])
								}
							}else{
								tmp_karyo_data <- karyo_data
							}
						}
						# IF ALL IN PREY
						else if ('all' == matrix_Prey[j,1]){
							#print("Prey all")
							if ('all' != matrix_Bait[i,1]){
								if (matrix_Bait[i,2] == ""){
									duprows <- rownames(karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]),]) %in% rownames(tmp_karyo_data)
									tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]),][!duprows,])
								} else{
									duprows <- rownames(karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3])),]) %in% rownames(tmp_karyo_data)
									tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3])),][!duprows,])
								}
							}else{
								tmp_karyo_data <- karyo_data
							}
						}
						# NO ALL
						else{
							#print("No all")
							if (matrix_Bait[i,2] == "" & matrix_Prey[j,2] == ""){
								#print("Both no position")
								duprows <- rownames(karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$Rname==matrix_Prey[j,1])),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$Rname==matrix_Prey[j,1])),][!duprows,])
							} else if (matrix_Bait[i,2] == "" & matrix_Prey[j,2] != ""){
								#print("Bait no position")
								duprows <- rownames(karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & ((karyo_data$Rname==matrix_Prey[j,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[j,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[j,3])))),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & ((karyo_data$Rname==matrix_Prey[j,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[j,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[j,3])))),][!duprows,])
							} else if (matrix_Bait[i,2] != "" & matrix_Prey[j,2] == ""){
								#print("Prey no position")
								duprows <- rownames(karyo_data[(((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & (karyo_data$Rname==matrix_Prey[j,1])),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & (karyo_data$Rname==matrix_Prey[j,1])),][!duprows,])
							} else{
								#print("Both position")
								duprows <- rownames(karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & ((karyo_data$Rname==matrix_Prey[j,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[j,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[j,3]))),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & ((karyo_data$Rname==matrix_Prey[j,1]) & (karyo_data$Junction >= matrix_Prey[j,2]) & (karyo_data$Junction <= matrix_Prey[j,3])),][!duprows,])
							}
						}
					}
				}
				# SAVE NEW DATA
				karyo_data <- tmp_karyo_data

			# IF LINKED
			} else{
				tmp_karyo_data <- c()
				duprows <- c()
				# TEST BAIT1/PREY1, BAIT2/PREY2 ETC...
				for (i in 1:length(matrix_Bait[,1])){
					#print(matrix_Bait[i,1])
					#print(matrix_Prey[i,1])
					# IF ALL IN BAIT
					if ('all' == matrix_Bait[i,1]){
						if ('all' != matrix_Prey[i,1]){
							#print("Bait all")
							if (matrix_Prey[i,2] == ""){
								duprows <- rownames(karyo_data[(karyo_data$Rname==matrix_Prey[i,1]),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$Rname==matrix_Prey[i,1]),][!duprows,])
							} else{
								duprows <- rownames(karyo_data[(karyo_data$Rname==matrix_Prey[i,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[i,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[i,3])),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$Rname==matrix_Prey[i,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[i,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[i,3])),][!duprows,])
							}
						}else{
							tmp_karyo_data <- karyo_data
						}
					}
					# IF ALL IN PREY
					else if ('all' == matrix_Prey[i,1]){
						#print("Prey all")
						if (matrix_Bait[i,2] == ""){
							if (matrix_Bait[i,2] == ""){
								duprows <- rownames(karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]),][!duprows,])
							} else{
								duprows <- rownames(karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3])),]) %in% rownames(tmp_karyo_data)
								tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(karyo_data$B_Rname==matrix_Bait[i,1]) & (strtoi(karyo_data$B_Rstart) >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3])),][!duprows,])
							}
						}else{
							tmp_karyo_data <- karyo_data
						}
					# NO ALL
					} else{
						#print("No all")
						if (matrix_Bait[i,2] == "" & matrix_Prey[i,2] == ""){
							#print("Both no position")
							duprows <- rownames(karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$Rname==matrix_Prey[i,1])),]) %in% rownames(tmp_karyo_data)
							tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$Rname==matrix_Prey[i,1])),][!duprows,])
						} else if (matrix_Bait[i,2] == "" & matrix_Prey[i,2] != ""){
							#print("Bait no position")
							duprows <- rownames(karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & ((karyo_data$Rname==matrix_Prey[i,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[i,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[i,3])))),]) %in% rownames(tmp_karyo_data)
							tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & ((karyo_data$Rname==matrix_Prey[i,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[i,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[i,3])))),][!duprows,])
						} else if (matrix_Bait[i,2] != "" & matrix_Prey[i,2] == ""){
							#print("Prey no position")
							duprows <- rownames(karyo_data[(((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & (karyo_data$Rname==matrix_Prey[i,1])),]) %in% rownames(tmp_karyo_data)
							tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[(((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & (karyo_data$Rname==matrix_Prey[i,1])),][!duprows,])
						} else{
							#print("Both position")
							duprows <- rownames(karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & ((karyo_data$Rname==matrix_Prey[i,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[i,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[i,3]))),]) %in% rownames(tmp_karyo_data)
							tmp_karyo_data <- rbind(tmp_karyo_data, karyo_data[((karyo_data$B_Rname==matrix_Bait[i,1]) & (karyo_data$B_Rstart >= strtoi(matrix_Bait[i,2])) & (karyo_data$B_Rend <= strtoi(matrix_Bait[i,3]))) & ((karyo_data$Rname==matrix_Prey[i,1]) & (karyo_data$Junction >= strtoi(matrix_Prey[i,2])) & (karyo_data$Junction <= strtoi(matrix_Prey[i,3]))),][!duprows,])
						}
					}
				}
				karyo_data <- tmp_karyo_data
			}
			#print("after bait/prey")
			#print(karyo_data)
			write(paste(c("Number of megajunction : ",nrow(karyo_data)),collapse=''),stderr())
			
			###FOR UCSC GENES LATER (IF ANY)
			genes_karyo_data <- karyo_data

			# GENOME DISPLAY
			chromosomesLength = read.table(opt$file_chrom_length, header = TRUE, sep="\t", stringsAsFactors = FALSE)
			cytoband = read.table(opt$file_cytoband, header = TRUE, sep="\t", stringsAsFactors = FALSE)

			# FULL GENOME DISPLAY
			if (opt$visualization == "full_genome"){
				custom.genome <- toGRanges(chromosomesLength)
				custom.cytobands <- toGRanges(cytoband)
			# SELECTED CHROMOSOME DISPLAY
			} else if (opt$visualization == "selected_chromosomes"){
				tmp_custom_genome <- c()
				duprows <- c()
				# SUB GENOME WITH SELECTED BAIT CHROMOSOME(S)
				if (!(is.data.frame(matrix_Bait) && nrow(matrix_Bait)==0)){
					for (i in 1:nrow(matrix_Bait[,1,drop=FALSE])){
						if ('all' == matrix_Bait[i,1,drop=FALSE]){
							duprows <- rownames(chromosomesLength) %in% rownames(tmp_custom_genome)
							tmp_custom_genome <- rbind(tmp_custom_genome, chromosomesLength[!duprows,])
							break
						} else{
							duprows <- rownames(chromosomesLength[(chromosomesLength$chr==matrix_Bait[i,1]),]) %in% rownames(tmp_custom_genome)
							tmp_custom_genome <- rbind(tmp_custom_genome, chromosomesLength[(chromosomesLength$chr==matrix_Bait[i,1]),][!duprows,])
						}
					}
				}
				# SUB GENOME WITH SELECTED PREY CHROMOSOME(S)
				if (!(is.data.frame(matrix_Prey) && nrow(matrix_Prey)==0)){
					for (i in 1:nrow(matrix_Prey[,1,drop=FALSE])){
						if ('all' == matrix_Prey[i,1,drop=FALSE]){
							duprows <- rownames(chromosomesLength) %in% rownames(tmp_custom_genome)
							tmp_custom_genome <- rbind(tmp_custom_genome, chromosomesLength[!duprows,])
							break
						} else{
							duprows <- rownames(chromosomesLength[(chromosomesLength$chr==matrix_Prey[i,1]),]) %in% rownames(tmp_custom_genome)
							tmp_custom_genome <- rbind(tmp_custom_genome, chromosomesLength[(chromosomesLength$chr==matrix_Prey[i,1]),][!duprows,])
						}
					}
				}
				tmp_custom_genome <- tmp_custom_genome[order(as.numeric(row.names(tmp_custom_genome))),]
				custom.genome <- toGRanges(tmp_custom_genome)
				custom.cytobands <- toGRanges(cytoband)
			}
			# SELECTED ZOOM IN DISPLAY
			else if (opt$visualization == "zoom_in"){
				tmp_custom_genome <- c()
				# IF ALL IN BAIT, DELETE IT AND REPLACE BY ALL THE CHROMS
				if ('all' %in% matrix_Bait[,1]){
					matrix_Bait <- matrix_Bait[!matrix_Bait[,1] == "all", ]
					for (i in 2:length(chrom_fasta_list)){
						matrix_Bait <- rbind(matrix_Bait, c(chrom_fasta_list[i],"",""))
					}
				}
				# IF ALL IN PREY, DELETE IT AND REPLACE BY ALL THE CHROMS
				if ('all' %in% matrix_Prey[,1]){
					matrix_Prey <- matrix_Prey[!matrix_Prey[,1] == "all", ]
					for (i in 2:length(chrom_fasta_list)){
						matrix_Prey <- rbind(matrix_Prey, c(chrom_fasta_list[i],"",""))
					}
				}
				# TEST EVERY CHROMOSOME
				for (i in 1:length(chrom_fasta_list)){
					chr <- ""
					start <- NULL
					end <- NULL
					# SUB BAIT AND PREY WITH CURRENT CHROMOSOME
					tmp_bait <- matrix_Bait[(matrix_Bait[,1] == chrom_fasta_list[i]), ,drop=FALSE]
					tmp_prey <- matrix_Prey[(matrix_Prey[,1] == chrom_fasta_list[i]), ,drop=FALSE]
					# SUB GENOME WITH SELECTED BAIT CURRENT CHROMOSOME
					if (nrow(tmp_bait) > 0){
						if (!(is.data.frame(tmp_bait) && nrow(tmp_bait)==0)){
							for (j in 1:nrow(tmp_bait)){
								if (is.null(tmp_custom_genome)){
									tmp_custom_genome <- rbind(tmp_custom_genome, c(chrom_fasta_list[i], tmp_bait[j,2], tmp_bait[j,3]))
								} else{								
									if (nrow(tmp_custom_genome[((tmp_custom_genome[,1] == chrom_fasta_list[i]) & (tmp_custom_genome[,2] == tmp_bait[j,2]) & (tmp_custom_genome[,3] == tmp_bait[j,3])), ,drop=FALSE]) == 0){
										tmp_custom_genome <- rbind(tmp_custom_genome, c(chrom_fasta_list[i], tmp_bait[j,2], tmp_bait[j,3]))
									}
								}
							}
						}
					}
					# SUB GENOME WITH SELECTED PREY CURRENT CHROMOSOME
					if (nrow(tmp_prey) > 0){
						if (!(is.data.frame(tmp_prey) && nrow(tmp_prey)==0)){
							for (j in 1:nrow(tmp_prey)){
								if (is.null(tmp_custom_genome)){
									tmp_custom_genome <- rbind(tmp_custom_genome, c(chrom_fasta_list[i], tmp_prey[j,2], tmp_prey[j,3]))
								} else{								
									if (nrow(tmp_custom_genome[((tmp_custom_genome[,1] == chrom_fasta_list[i]) & (tmp_custom_genome[,2] == tmp_prey[j,2]) & (tmp_custom_genome[,3] == tmp_prey[j,3])), ,drop=FALSE]) == 0){
										tmp_custom_genome <- rbind(tmp_custom_genome, c(chrom_fasta_list[i], tmp_prey[j,2], tmp_prey[j,3]))
									}
								}
							}
						}
					}
				}
				if (!(is.data.frame(tmp_custom_genome) && nrow(tmp_custom_genome)==0)){
					for (i in 1:nrow(tmp_custom_genome)){
						if ((tmp_custom_genome[i,2] == "") & (tmp_custom_genome[i,3] == "")){
							if (!(is.data.frame(chromosomesLength) && nrow(chromosomesLength)==0)){
								for (j in 1:nrow(chromosomesLength)){
									if (tmp_custom_genome[i,1] == chromosomesLength[j,1]){
										tmp_custom_genome[i,2] <- chromosomesLength[j,2]
										tmp_custom_genome[i,3] <- chromosomesLength[j,3]
									}
								}
							}
						}
					}
				}
				custom.genome <- toGRanges(data.frame(tmp_custom_genome[,1], as.integer(unlist(tmp_custom_genome[,2])), as.integer(unlist(tmp_custom_genome[,3]))))
				custom.cytobands <- toGRanges(cytoband)
				if (opt$visualization == "zoom_in"){ 
					if (file.exists(paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]),"_Karyoplot_freq",file_input_extension), collapse=''))){
						custom.frequency <- toGRanges(read.table(paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]),"_Karyoplot_freq",file_input_extension), collapse=''), header = TRUE, sep="\t", stringsAsFactors = FALSE))
					} else{
						custom.frequency <- NULL
					}
				} else {
					custom.frequency <- NULL
				}
			} else {
				write("Error : You have to set a visualization (full_genome, selected_chromosomes or zoom_in) !",stderr())
				write("\n",stderr())
				cat(getopt(spec, usage=TRUE))
				q(status=1)
			}

			# RENAME CUSTOM.GENOME, CUSTOM.CYTOBANDS, KARYO DATA, LOCUS DATA, KARYO FREQ AND CHROMOSOME NAMES
			if (!is.null(opt$file_rename)) {
				# TO CHANGE FULL CHR (ex:chr9)
				if (!is.null(tmp_df_rename_chr)){
					if (!(is.data.frame(tmp_df_rename_chr) && nrow(tmp_df_rename_chr)==0)){
						for (i in 1:nrow(tmp_df_rename_chr)){
							# CHECK FULL CHROM IN FASTA LIST
							#print("CHECK FULL CHROM IN FASTA LIST")
							if (tmp_df_rename_chr$current[i] %in% chrom_fasta_list){
								# RENAME CUSTOM.GENOME
								suppressWarnings(seqlevels(custom.genome) <- sub(tmp_df_rename_chr$current[i],tmp_df_rename_chr$new[i],seqlevels(custom.genome)))
								# RENAME CUSTOM.CYTOBANDS
								suppressWarnings(seqlevels(custom.cytobands) <- sub(tmp_df_rename_chr$current[i],tmp_df_rename_chr$new[i],seqlevels(custom.cytobands)))
								# RENAME CUSTOM.FREQUENCY
								if (!is.null(custom.frequency)){
									suppressWarnings(seqlevels(custom.frequency) <- sub(tmp_df_rename_chr$current[i],tmp_df_rename_chr$new[i],seqlevels(custom.frequency)))
								}
								# RENAME KARYO DATA
								if (!(is.data.frame(karyo_data) && nrow(karyo_data)==0)){
									for (j in 1:nrow(karyo_data)){
										if (tmp_df_rename_chr$current[i] == karyo_data$Rname[j]){
											karyo_data$Rname[j] <- tmp_df_rename_chr$new[i]
										}
										if (tmp_df_rename_chr$current[i] == karyo_data$B_Rname[j]){
											karyo_data$B_Rname[j] <- tmp_df_rename_chr$new[i]
										}
									}
								}
								# RENAME LOCUS
								if (!is.null(opt$file_locus)) {
									if (!(is.data.frame(df_locus) && nrow(df_locus)==0)){
										for (j in 1:nrow(df_locus)){
											if (tmp_df_rename_chr$current[i] == df_locus$chr[j]){
												df_locus$chr[j] <- tmp_df_rename_chr$new[i]
											}
										}
									}
								}
							}
							# WARN IF THE RENAME LINE IS NOT USE
							else{
								write(paste(c("Warning : The chromosome ", tmp_df_rename_chr$current[i], " in the rename file does not exist in your reference genome !"), collapse=""),stderr())
								write(paste(c("Warning : {", tmp_df_rename_chr$current[i], "} will not be use for the rename!"), collapse=""),stderr())
							}
						}
					}
				}
				# GOT TO CHANGE ZOOMED CHR (ex:chr6:150:200)
				if (!is.null(tmp_df_rename_chrpart)){
					if (!(is.data.frame(tmp_df_rename_chrpart) && nrow(tmp_df_rename_chrpart)==0)){
						for (i in 1:nrow(tmp_df_rename_chrpart)){
							new_chr_name <- ""
							start <- 0
							end <- 0
							# CHECK IF CHRPART IS NOT MISSPELLED
							if (length(strsplit(tmp_df_rename_chrpart$current[i],":")[[1]]) != 3){
								write(paste(c("Error : The chromosome ", tmp_df_rename_chrpart$current[i], " is misspelled !"), collapse=""),stderr())
								write("\n",stderr())
								cat(getopt(spec, usage=TRUE))
								q(status=1)
							}
							# IF CHR IS IN REFERENCE FILE
							if ((strsplit(tmp_df_rename_chrpart$current[i],":")[[1]][1] %in% chrom_fasta_list) & (tmp_df_rename_chrpart$current[i] %in% strsplit(opt$chr_bait,",")[[1]] | tmp_df_rename_chrpart$current[i] %in% strsplit(opt$chr_prey,",")[[1]])){
								#print("IF CHR IS IN REFERENCE FILE")
								# CHANGE CHR NAME IF ALREADY CHANGE WITH FULL CHR
								if (strsplit(tmp_df_rename_chrpart$current[i],":")[[1]][1] %in% tmp_df_rename_chr$current){
									new_chr_name <- tmp_df_rename_chr[tmp_df_rename_chr$current == strsplit(tmp_df_rename_chrpart$current[i],":")[[1]][1],]$new
								} else{
									new_chr_name <- strsplit(tmp_df_rename_chrpart$current[i],":")[[1]][1]
								}
								start <- strsplit(tmp_df_rename_chrpart$current[i],":")[[1]][2]
								end <- strsplit(tmp_df_rename_chrpart$current[i],":")[[1]][3]

								# RENAME CUSTOM.GENOME
								#print("--------------------------------------------------------------custom.genome")
								#print(custom.genome)
								for(j in 1:length(custom.genome)) {
									if ((as.character(seqnames(custom.genome[j])) == new_chr_name) & (start(custom.genome[j]) == strtoi(start)) & (end(custom.genome[j]) == strtoi(end))){
										suppressWarnings(seqlevels(custom.genome[j])<- sub(new_chr_name,tmp_df_rename_chrpart$new[i],seqlevels(custom.genome[j])))
										if (!new_chr_name %in% seqnames(custom.genome)){
											custom.genome <- dropSeqlevels(custom.genome, new_chr_name)
										}
									}
								}
								#print("OUT")
								#print(custom.genome)
								# RENAME CUSTOM.CYTOBANDS
								#print("--------------------------------------------------------------custom.cytobands")
								copy_cytoband <- custom.cytobands[(seqnames(custom.cytobands) == new_chr_name) & (end(custom.cytobands) >= strtoi(start)) & (start(custom.cytobands) < strtoi(end))]
								
								start(copy_cytoband[1]) <- strtoi(start)
								end(copy_cytoband[length(copy_cytoband)]) <- strtoi(end)
								seqlevels(copy_cytoband)<- sub(new_chr_name,tmp_df_rename_chrpart$new[i],seqlevels(copy_cytoband)) 
								suppressWarnings(custom.cytobands <- append(custom.cytobands, copy_cytoband))
								#print(custom.cytobands)

								# RENAME CUSTOM.FREQUENCY
								#print("--------------------------------------------------------------custom.frequency")
								if (!is.null(custom.frequency)){
									copy_frequency <- custom.frequency[(seqnames(custom.frequency) == new_chr_name) & (end(custom.frequency) >= strtoi(start)) & (start(custom.frequency) < strtoi(end))]				
									start(copy_frequency[1]) <- strtoi(start)
									end(copy_frequency[length(copy_frequency)]) <- strtoi(end)
									seqlevels(copy_frequency)<- sub(new_chr_name,tmp_df_rename_chrpart$new[i],seqlevels(copy_frequency)) 
									suppressWarnings(custom.frequency <- append(custom.frequency, copy_frequency))
								}
								# RENAME KARYO DATA
								#print("--------------------------------------------------------------karyo_data")
								if (!(is.data.frame(karyo_data) && nrow(karyo_data)==0)){
									for (j in 1:nrow(karyo_data)){
										if ((new_chr_name == karyo_data$Rname[j]) & (strtoi(start) <= strtoi(karyo_data$Junction[j])) & (strtoi(end) >= strtoi(karyo_data$Junction[j]))){
											karyo_data$Rname[j] <- tmp_df_rename_chrpart$new[i]
										}
										if ((new_chr_name == karyo_data$B_Rname[j]) & (strtoi(start) <= strtoi(karyo_data$B_Rstart[j])) & (strtoi(end) >= strtoi(karyo_data$B_Rend[j]))){
											karyo_data$B_Rname[j] <- tmp_df_rename_chrpart$new[i]
										}
									}
								}

								# RENAME LOCUS
								#print("--------------------------------------------------------------df_locus")
								# CHOICE TO MAKE HERE ! I COULD RATHER DUPLICATE LOCUS WITH NEW NAME OR REPLACE THEM, I REPLACE THEM FOR NOW
								
								if (!is.null(opt$file_locus)) {
									if (!(is.data.frame(df_locus) && nrow(df_locus)==0)){
										for (j in 1:nrow(df_locus)){
											if ((new_chr_name == df_locus$chr[j]) & (strtoi(start) <= strtoi(df_locus$start[j])) & (strtoi(end) >= strtoi(df_locus$end[j]))){
												df_locus$chr[j] <- tmp_df_rename_chrpart$new[i]
											}
										}
									}
								}
							}
							# WARN IF THE RENAME LINE IS NOT USE
							else{
								write(paste(c("Warning : The chromosome ", tmp_df_rename_chrpart$current[i], " in the rename file does not exist in your reference genome !"), collapse=""),stderr())
								write(paste(c("Warning : {", tmp_df_rename_chrpart$current[i], "} will not be use for the rename !"), collapse=""),stderr())
							}
						}
					}
				}		
			}

			# ADD LOCUS
			if (!is.null(opt$file_locus)) {
				locus_regions <- toGRanges(data.frame(df_locus))
			}

			#print("DATA BEFORE")
			#print(karyo_data)
			#print("UNSCALED GENOMES")
			#print(custom.genome)

			# APPLY SCALE MODIFICATION TO ALL VALUES (GENOME, CYTOBANDS, LOCUS, KARYO DATA) IF MORE THAN 1 CHR TO DISPLAY
			#print("SCALE MODIF")
			tmp_granges_cytobands <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
			if (!is.null(custom.frequency)){
				tmp_granges_frequency <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
			}
			if (!is.null(opt$file_locus)) {
				tmp_locus <- df_locus[FALSE,]
			}
			
			if (opt$visualization == "zoom_in" & length(custom.genome) > 1){
				#print("HARD CHANGE...")
				max_width = max(width(custom.genome))
				# SCALE CUSTOM.GENOME
				for (i in 1:length(custom.genome)){
					#print(custom.genome[i])
					#print("WHAT WE WANT TO MODIFY")
					#print(custom.genome[i])
					# SAVE OLD VALUES
					old_start <- start(custom.genome[i])
					old_end <- end(custom.genome[i])
					# SET UP SCALE
					width_current <- width(custom.genome[i])
					scale_current <- trunc(max_width/width_current)
					#print("SCALE")
					#print(scale_current)

					# MODIFY START AND END OF GENOME
					start(custom.genome[i]) <- 1
					end(custom.genome[i]) <- width_current*scale_current
					# SCALE CUSTOM.CYTOBANDS
					tmp_cytobands <- custom.cytobands[(seqnames(custom.cytobands) == as.character(seqnames(custom.genome[i]))) & (end(custom.cytobands) >= strtoi(old_start)) & (start(custom.cytobands) < strtoi(old_end))]
					# AJUST START(1) AND END(-1) WITH THE GENOME START AND END
					#print("UNSCALED CYTOBAND")
					end(tmp_cytobands[length(tmp_cytobands)]) <- old_end
					#print(tmp_cytobands)
					#print("CYTOBANDS")
					# CHANGE CYTOBANDS
					if (length(tmp_cytobands) == 1){
						start(tmp_cytobands[1]) <- 1
						end(tmp_cytobands[1]) <- width_current*scale_current
					} else{
						for (j in 1:length(tmp_cytobands)){
							if (j == 1){
								start(tmp_cytobands[j]) <- strtoi(((start(tmp_cytobands[j]) - strtoi(old_start)) * scale_current)+1)
								end(tmp_cytobands[j]) <- strtoi(((end(tmp_cytobands[j]) - strtoi(old_start))+1) * scale_current)
							}
							else{
								if (strtoi(((start(tmp_cytobands[j]) - strtoi(old_start))+1) * scale_current) > end(tmp_cytobands[j])){
									end(tmp_cytobands[j]) <- strtoi(((end(tmp_cytobands[j]) - strtoi(old_start))+1) * scale_current)
									start(tmp_cytobands[j]) <- strtoi(((start(tmp_cytobands[j]) - strtoi(old_start))+1) * scale_current)
								} else{
									start(tmp_cytobands[j]) <- strtoi(((start(tmp_cytobands[j]) - strtoi(old_start))+1) * scale_current)
									end(tmp_cytobands[j]) <- strtoi(((end(tmp_cytobands[j]) - strtoi(old_start))+1) * scale_current)
								}
							}	
						}
					}
					tmp_granges_cytobands <- append(tmp_granges_cytobands, tmp_cytobands)
					if (!is.null(custom.frequency)){
						#print("FREQUENCY")
						# SCALE CUSTOM.FREQUENCY
						tmp_frequency <- custom.frequency[(seqnames(custom.frequency) == as.character(seqnames(custom.genome[i]))) & (end(custom.frequency) >= strtoi(old_start)) & (start(custom.frequency) < strtoi(old_end))]
						#print(tmp_frequency)
						# AJUST START(1) AND END(-1) WITH THE GENOME START AND END
						#print("UNSCALED FREQUENCY")
						start(tmp_frequency[1]) <- old_start
						end(tmp_frequency[length(tmp_frequency)]) <- old_end
						#print(tmp_frequency)
						# CHANGE FREQUENCY
						start1 <- strtoi(((start(tmp_frequency[1]) - strtoi(old_start)) * scale_current)+1)
						start(tmp_frequency) <- strtoi(((start(tmp_frequency) - strtoi(old_start))+1) * scale_current)
						end(tmp_frequency) <- strtoi(((end(tmp_frequency) - strtoi(old_start))+1) * scale_current)
						start(tmp_frequency[1]) <- start1
						#print("OUT FREQUENCY")
						tmp_granges_frequency <- append(tmp_granges_frequency, tmp_frequency)
					}
					#print("DATA")
					# SCALE KARYO DATA
					if (!(is.data.frame(karyo_data) && nrow(karyo_data)==0)){
						for (j in 1:nrow(karyo_data)){
							if (karyo_data$Rname[j] == as.character(seqnames(custom.genome[i]))){
								karyo_data$Junction[j] <- strtoi(((strtoi(karyo_data$Junction[j]) - strtoi(old_start)) * scale_current)+1)
								karyo_data$Rstart[j] <- strtoi(((strtoi(karyo_data$Rstart[j]) - strtoi(old_start)) * scale_current)+1)
								karyo_data$Rend[j] <- strtoi(((strtoi(karyo_data$Rend[j]) - strtoi(old_start)) * scale_current)+1)
							}
							if (karyo_data$B_Rname[j] == as.character(seqnames(custom.genome[i]))){
								karyo_data$B_Rstart[j] <- strtoi(((strtoi(karyo_data$B_Rstart[j]) - strtoi(old_start)) * scale_current)+1)
								karyo_data$B_Rend[j] <- strtoi(((strtoi(karyo_data$B_Rend[j]) - strtoi(old_start)) * scale_current)+1)
							}
						}
					}
					#print("LOCUS")
					# SCALE LOCUS
					if (!is.null(opt$file_locus)) {
						if (!(is.data.frame(df_locus) && nrow(df_locus)==0)){
							for (j in 1:nrow(df_locus)){
								if (df_locus$chr[j] == as.character(seqnames(custom.genome[i]))){
									df_locus$start[j] <- strtoi(((strtoi(df_locus$start[j]) - strtoi(old_start)) * scale_current)+1)
									df_locus$end[j] <- strtoi((strtoi(df_locus$end[j]) - strtoi(old_start)+1) * scale_current)
									if (!(is.na(df_locus$start[j])) | !(is.na(df_locus$end[j]))){
										tmp_locus <- rbind(tmp_locus, df_locus[j,])
									}
								}
							}
						}
					}

				}
				custom.cytobands <- tmp_granges_cytobands
				if (!is.null(custom.frequency)){
					custom.frequency <- tmp_granges_frequency
					tmp_granges_frequency <- NULL
				}
				if (!is.null(opt$file_locus)) {
					df_locus <- tmp_locus
				}
			}
			#print("END SCALE MODIF")
			#print("SCALED GENOMES")
			#print(custom.genome)
			#print("SCALED CYTOBANDS")
			#print(custom.cytobands)
			#print("SCALED KARYO DATA")
			#print(karyo_data)
			#print("---------------------------------------------------------")
			#print("OLD LOCUS")
			#print(save_locus)
			#print("SCALED LOCUS")
			#print(df_locus)

			# MODIFY CUSTOME.FREQUENCY TO FAKE READ NUMBER TO PLOT A DENSITY
			if (!is.null(custom.frequency)){

				custom.frequency.plus <- custom.frequency[custom.frequency$Strand == '+']
				tmp_granges_frequency <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL, Data=NULL, Strand=NULL))
				for (i in 1:length(custom.frequency.plus)){
					if (custom.frequency.plus[i]$Data != 0){
						factor <- round(custom.frequency.plus[i]$Data*1000)
						tmp_granges_frequency <- append(tmp_granges_frequency, custom.frequency.plus[rep(i, factor), ])
					}
				}
				custom.frequency.plus <- tmp_granges_frequency

				custom.frequency.minus <- custom.frequency[custom.frequency$Strand == '-']
				tmp_granges_frequency <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL, Data=NULL, Strand=NULL))
				for (i in 1:length(custom.frequency.minus)){
					if (custom.frequency.minus[i]$Data != 0){
						factor <- round(custom.frequency.minus[i]$Data*1000)
						tmp_granges_frequency <- append(tmp_granges_frequency, custom.frequency.minus[rep(i, factor), ])
					}
				}
				custom.frequency.minus <- tmp_granges_frequency
			}

			# DISPLAY LOCUS FILE INFORMATION
			if (!is.null(opt$file_locus)) {
				if (!(is.data.frame(df_locus) && nrow(df_locus)==0)){
					# CHANGE NAMES TO GREEK NAMES IF GREEK OPTION TRUE
					if (opt$greek){
						df_locus$name <- gsub("\'", "\\\\'", df_locus$name)
						base_name <- ""
						base_name <- df_locus$name
 						df_locus$name <- gsub("alpha", "'~alpha~'", df_locus$name)
						df_locus$name <- gsub("Alpha", "'~Alpha~'", df_locus$name)
						df_locus$name <- gsub("beta", "'~beta~'", df_locus$name)
						df_locus$name <- gsub("Beta", "'~Beta~'", df_locus$name)
						df_locus$name <- gsub("gamma", "'~gamma~'", df_locus$name)
						df_locus$name <- gsub("Gamma", "'~Gamma~'", df_locus$name)
						df_locus$name <- gsub("delta", "'~delta~'", df_locus$name)
						df_locus$name <- gsub("Delta", "'~Delta~'", df_locus$name)
						df_locus$name <- gsub("epsilon", "'~epsilon~'", df_locus$name)
						df_locus$name <- gsub("Epsilon", "'~Epsilon~'", df_locus$name)
						df_locus$name <- gsub("zeta", "'~zeta~'", df_locus$name)
						df_locus$name <- gsub("Zeta", "'~Zeta~'", df_locus$name)
						df_locus$name <- gsub("theta", "'~theta~'", df_locus$name)
						df_locus$name <- gsub("Theta", "'~Theta~'", df_locus$name)
						df_locus$name <- gsub("iota", "'~iota~'", df_locus$name)
						df_locus$name <- gsub("Iota", "'~Iota~'", df_locus$name)
						df_locus$name <- gsub("kappa", "'~kappa~'", df_locus$name)
						df_locus$name <- gsub("Kappa", "'~Kappa~'", df_locus$name)
						df_locus$name <- gsub("lambda", "'~lambda~'", df_locus$name)
						df_locus$name <- gsub("Lambda", "'~Lambda~'", df_locus$name)
						df_locus$name <- gsub("mu", "'~mu~'", df_locus$name)
						df_locus$name <- gsub("Mu", "'~Mu~'", df_locus$name)
						df_locus$name <- gsub("nu", "'~nu~'", df_locus$name)
						df_locus$name <- gsub("Nu", "'~Nu~'", df_locus$name)
						df_locus$name <- gsub("xi", "'~xi~'", df_locus$name)
						df_locus$name <- gsub("Xi", "'~Xi~'", df_locus$name)
						df_locus$name <- gsub("omicron", "'~omicron~'", df_locus$name)
						df_locus$name <- gsub("Omicron", "'~Omicron~'", df_locus$name)
						df_locus$name <- gsub("pi", "'~pi~'", df_locus$name)
						df_locus$name <- gsub("Pi", "'~Pi~'", df_locus$name)
						df_locus$name <- gsub("rho", "'~rho~'", df_locus$name)
						df_locus$name <- gsub("Rho", "'~Rho~'", df_locus$name)
						df_locus$name <- gsub("sigma", "'~sigma~'", df_locus$name)
						df_locus$name <- gsub("Sigma", "'~Sigma~'", df_locus$name)
						df_locus$name <- gsub("tau", "'~tau~'", df_locus$name)
						df_locus$name <- gsub("Tau", "'~Tau~'", df_locus$name)
						df_locus$name <- gsub("upsilon", "'~upsilon~'", df_locus$name)
						df_locus$name <- gsub("Upsilon", "'~Upsilon~'", df_locus$name)
						df_locus$name <- gsub("phi", "'~phi~'", df_locus$name)
						df_locus$name <- gsub("Phi", "'~Phi~'", df_locus$name)
						df_locus$name <- gsub("chi", "'~chi~'", df_locus$name)
						df_locus$name <- gsub("Chi", "'~Chi~'", df_locus$name)
						df_locus$name <- gsub("omega", "'~omega~'", df_locus$name)
						df_locus$name <- gsub("Omega", "'~Omega~'", df_locus$name)

						# PROBLEMS WITH ETA AND PSI (BECAUSE ETA IN BETA/ZETA AND PSI IN EPSILON)
						# LOWER CASE ETA
						locate_eta <- str_locate_all(df_locus$name, "eta")
						char_added=0
						for (i in 1:length(locate_eta)){
							if (nrow(locate_eta[[i]]) != 0){
								for (j in 1:nrow(locate_eta[[i]])){
									if (!(substr(df_locus$name[i],locate_eta[[i]][j,1]+char_added-1,locate_eta[[i]][j,2]+char_added+3) == "beta" | substr(df_locus$name[i],locate_eta[[i]][j,1]+char_added-1,locate_eta[[i]][j,2]+char_added+3) == "zeta" | substr(df_locus$name[i],locate_eta[[i]][j,1]+char_added-1,locate_eta[[i]][j,2]+char_added+3) == "theta")){
										df_locus$name[i] <- paste0(substr(df_locus$name[i], 0, locate_eta[[i]][j,1]+char_added-1), "'~psi~'", substr(df_locus$name[i], locate_eta[[i]][j,2]+char_added+1, nchar(df_locus$name[i])))
										char_added=char_added+4
									}
								}
							}
						}
						# UPPER CASE ETA
						df_locus$name <- gsub("Eta", "'~Eta~'", df_locus$name)

						# LOWER CASE PSI
						locate_psi <- str_locate_all(df_locus$name, "psi")
						char_added=0
						for (i in 1:length(locate_psi)){
							if (nrow(locate_psi[[i]]) != 0){
								for (j in 1:nrow(locate_psi[[i]])){
									if (!(substr(df_locus$name[i],locate_psi[[i]][j,1]+char_added-1,locate_psi[[i]][j,2]+char_added+3) == "epsilon" | substr(df_locus$name[i],locate_psi[[i]][j,1]+char_added-1,locate_psi[[i]][j,2]+char_added+3) == "upsilon")){
										df_locus$name[i] <- paste0(substr(df_locus$name[i], 0, locate_psi[[i]][j,1]+char_added-1), "'~psi~'", substr(df_locus$name[i], locate_psi[[i]][j,2]+char_added+1, nchar(df_locus$name[i])))
										char_added=char_added+4
									}
								}
							}
						}
						# UPPER CASE PSI
						df_locus$name <- gsub("Psi", "'~Psi~'", df_locus$name)

						# TRANSFORM IN EXPRESSION
						names_markers <- c()
						for (i in 1:nrow(df_locus)){
							if (substr(df_locus$name[i],1,1) == "'"){
								df_locus$name[i] <- paste0("expression(",substr(df_locus$name[i],3,nchar(df_locus$name[i])))
							} else{
								df_locus$name[i] <- paste0("expression('",df_locus$name[i])
							}
							if (str_sub(df_locus$name[i],-1,-1) == "'"){
								df_locus$name[i] <- paste0(substr(df_locus$name[i],1,nchar(df_locus$name[i])-2),")")
							} else{
								df_locus$name[i] <- paste0(df_locus$name[i],"')")
							}
							names_markers <- c(names_markers,eval(parse(text=df_locus$name[i])))
						}
					} else{
						#df_locus$name <- gsub("\'", "\'", df_locus$name)
						names_markers <- df_locus$name
					}
				}
			}

			pp <- getDefaultPlotParams(plot.type=1)

			# CHANGE THE IDEOGRAM HEIGHT PARAM TO CREATE THICKER IDEOGRAMS 
			pp$ideogramheight <- 0

			# SET UP KARYOPLOT
			#print("custom.genome")
			#print(custom.genome)

			# OUTPUT LINK PDF
			
			# CREATE PDF NAME
			pdf(file=paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]), file_output_extension_link), collapse=''), height=8, width=8);
			
			#kp_link <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, main = expression("Karyoplot of double stranded break areas"), cex=0.6, plot.params=pp)
			kp_link <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, main = expression("Karyoplot of double stranded break areas"), cex=0.6, ideogram.plotter=NULL)

			# ADD TICK MARKERS IF NOT ZOOM IN AND SOLO CUSTOM GENOME
			if (!(opt$visualization == "zoom_in" & length(custom.genome) > 1)){
				metric <- 10^(strtoi(nchar(toString(max(width(custom.genome)))))-2)
				kpAddBaseNumbers(kp_link, tick.dist = metric, tick.len = 15, tick.col="red", cex=0.4, minor.tick.dist = metric/10, minor.tick.len = 5, minor.tick.col = "gray")
			}
			if (!is.null(opt$file_locus)) {
				kpPlotRegions(kp_link, toGRanges(df_locus), col="#7c7c7c", border="#000000", r0=0, r1=0.05)
			}
			#print("DATA")
			#print(karyo_data)

			# ADD LINKS AND REGIONS
			total_bin <- LinksAndRegions(kp_link, karyo_data, value_color, 6, 0.11, 0.26, 0.06, 0.11)
			# UCSC GENES
			if (opt$UCSC) {
				# SET UP DATAFRAME
				df_genes <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("name", "start", "end"))
				if (length(total_bin)!=0){
					# FOR EACH JUNCTION
					for (i in 1:length(total_bin)){
						line_gene = c()
						check_found_rename=FALSE
						if (!is.null(opt$file_rename)) {
							# IF CHROM EXIST IN RENAME
							for (j in 1:length(df_rename$new)){
								if (as.character(seqnames(total_bin[i])) == df_rename$new[j]){
									check_found_rename = TRUE
									# REVERSE CHROMOSOME NAME
									line_gene <- append(line_gene, strsplit(df_rename$current[j],":")[[1]][1])
								}
							}
						}
						# IF NOT, THIS IS THE REAL CHROMOSOME NAME
						if (!check_found_rename){
							line_gene <- append(line_gene, as.character(seqnames(total_bin[i])))
						}

						# TAKE POSITIONS TO FIND ASSOCIATED GENE
						line_gene <- append(line_gene, strtoi(start(total_bin[i])))
						line_gene <- append(line_gene, strtoi(end(total_bin[i])))

						# CREATE A DATAFRAME
						df_genes[nrow(df_genes)+1,] <- line_gene
					}
					# MAKE COLUMNS WITH NUMBER INTEGER
					df_genes$start <- as.integer(df_genes$start)
					df_genes$end <- as.integer(df_genes$end)

					# IF THERE WAS A CONSTRUCTION, IT IS NEEDED TO REVERSE JUNCTION POSITIONS
					if (!is.null(opt$file_construction)) {
						add_to_start=0
						add_to_end=0
						current_length=0
						# MODIFY THE CONSTRUCTION DATAFRAME TO FIT THE ACTUAL DATA
						df_construction_genes <- df_construction
						for (i in 1:nrow(df_construction)){
							current_length <- as.integer(df_construction$length[i])
							df_construction$start[i] <- as.integer(df_construction$start[i])+as.integer(add_to_start)
							add_to_end = as.integer(add_to_end) + as.integer(current_length)
							df_construction$end[i] <- as.integer(df_construction$end[i])+as.integer(add_to_end)
							add_to_start = as.integer(add_to_start) + as.integer(current_length)
						}
						#print("MODIFIED CONSTRUCTION")
						#print(df_construction)
						# REVERSE JUNCTION POSITION
						i <- 1
						while (i<=nrow(genes_karyo_data)){
							result_Modify_genome_position <- Modify_genome_position(df_construction, genes_karyo_data$Rname[i],genes_karyo_data$Junction[i], genes_karyo_data$Junction[i])
							genes_karyo_data$Rstart[i] <- result_Modify_genome_position[1]
							genes_karyo_data$Rend[i] <- result_Modify_genome_position[2]
							if (genes_karyo_data$Rstart[i] == -1 && genes_karyo_data$Rend[i] == -1){
								genes_karyo_data <- genes_karyo_data[-c(i), ]
							}else{
								i = i+1
							}
						}
					}

					# REMOVE UNUSED COLUMNS
					genes_karyo_data$Junction <- NULL
					genes_karyo_data$B_Rname <- NULL
					genes_karyo_data$B_Rstart <- NULL
					genes_karyo_data$B_Rend <- NULL
					genes_karyo_data$value <- NULL
					genes_karyo_data$Rstart <- as.numeric(as.character(genes_karyo_data$Rstart))
					genes_karyo_data$Rend <- as.numeric(as.character(genes_karyo_data$Rend))

					# RETRIEVE GENES
					if (nrow(genes_karyo_data) > 0){
						genes <- RetrieveGenes(toGRanges(genes_karyo_data))
					}else{
						genes=c()
					}
					# IF GENES RESULT
					if (length(genes) != 0){
						seqlevelsStyle(genes) <- "UCSC"
						# MODIFY GENE POSITIONS IF CONSTRUCTION FILE TO FIT THE KARYOPLOT
						if (!is.null(opt$file_construction)) {
							for (i in 1:length(genes)){
								result_Modify_genes_position <- Modify_genes_position(df_construction_genes ,as.character(seqnames(genes[i])), start(genes[i]), end(genes[i]))
								start(genes[i]) <- result_Modify_genes_position[1]
								end(genes[i]) <- result_Modify_genes_position[2]
							}
						}
						# PRINT GENES
						suppressWarnings(kpPlotMarkers(kp_link, data=genes, text.orientation = "vertical", label.color="black", labels=genes$external_gene_name, r0=0, r1=-0.05, cex=0.3, adjust.label.position = TRUE))
					}
				}
			}
			# PRINT LOCI
			if (!is.null(opt$file_locus)) {
				#kpPlotMarkers(kp_link, data=toGRanges(df_locus), text.orientation = "horizontal", label.color="black", labels=names_markers, label.margin=0, r0=0.15, r1=0.1, cex=0.4, marker.parts=c(0,0,0), adjust.label.position = TRUE)
				kpPlotMarkers(kp_link, data=toGRanges(df_locus), text.orientation = "horizontal", label.color="black", labels=names_markers, r0=0, r1=-0.05, marker.parts=c(0,00,0), cex=0.4, adjust.label.position = FALSE)
			}

			# CLOSE THE GRAPHIC DEVICE AND CLEAR MEMORY
			dev.off()

			if (!is.null(custom.frequency)){
			
				# OUTPUT HISTO PDF

				# CREATE PDF NAME
				pdf(file=paste(c(opt$dir_post,as.vector(metadata$Library[library]),"/",as.vector(metadata$Library[library]), file_output_extension_histo), collapse=''), height=8, width=8);

				seqlevels(custom.genome) <- c(as.character(seqnames(custom.genome)))

				kp_histo <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, main = expression("Karyoplot of double stranded break areas"), cex=0.6, ideogram.plotter=NULL)

				if (length(custom.frequency.minus) > 0){
					kpPlotCoverage(kp_histo, data=custom.frequency.minus, r0=0.50, r1=0.05, col="#eb0c0c",ymax=15000, clipping=FALSE, show.0.cov=FALSE)
				}
				if (length(custom.frequency.plus) > 0){
					kpPlotCoverage(kp_histo, data=custom.frequency.plus, r0=0.505, r1=0.955, col="#20eb0d",ymax=15000, clipping=FALSE, show.0.cov=FALSE)
				}

				# ADD TICK MARKERS IF NOT ZOOM IN AND SOLO CUSTOM GENOME
				if (!(opt$visualization == "zoom_in" & length(custom.genome) > 1)){
					metric <- 10^(strtoi(nchar(toString(max(width(custom.genome)))))-2)
					kpAddBaseNumbers(kp_histo, tick.dist = metric, tick.len = 15, tick.col="red", cex=0.4, minor.tick.dist = metric/10, minor.tick.len = 5, minor.tick.col = "gray")
				}
				if (!is.null(opt$file_locus)) {
					kpPlotRegions(kp_histo, toGRanges(df_locus), col="#7c7c7c", border="#000000", r0=0, r1=0.05)
				}
				# PRINT LOCI
				if (!is.null(opt$file_locus)) {
					kpPlotMarkers(kp_histo, data=toGRanges(df_locus), text.orientation = "horizontal", label.color="black", labels=names_markers, r0=0, r1=-0.05, marker.parts=c(0,00,0), cex=0.4, adjust.label.position = FALSE)
				}
				# CLOSE THE GRAPHIC DEVICE AND CLEAR MEMORY
				dev.off()
			}
		}
	}
}
