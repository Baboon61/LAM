### A script allowing to split junctions to legitimate, illegitimate or trash junctions. Also, able to categorize legitimate junctions to different input locis, group illegitimate junctions to pool in order to find repeat locis and gene location close to this pool.

##########################################################################IMPORTS##########################################################################

import os
import sys
import getopt
import os.path
import shutil
import subprocess
import time
import numpy as np
import pandas as pd
from Bio import BiopythonWarning 
from Bio import SeqIO 
import pyensembl as pe
import sys
from datetime import datetime
import csv
import collections

###Experimentaly we want to find chr6/chr12 or chr12/chr12 junctions depends on the reference genome file. To split junctions, a file has to be create with these fields : bait chromosome, prey chromosome, start and end position on prey chromosome to delimit legitimate junction. Once a junction is called legitimite, it will be send to the closest locus set up in the input locus file. While illegitimate junctions will be pool according to the size pool filter, then every pool will be match against repeatMasker database to find repeat event but also against Ensembl to find attach gene to pools. In the same times validation on junctions are made to see if the junction is a viable junction (min_gap option). It is also possible for duplicated parts to send junctions from donor to the acceptor position of this junction (-t). To translate a junction the distance between start locus to junction on donor chromosome will be add to the start of the acceptor chromosome to find a new junction. As some locus between donor and acceptor are not the same length, if a junction fall after the end of a locus, the junction will be add at the end of it.

#data=pe.EnsemblRelease(release=85, species=pe.Species(latin_name='homo_sapiens', synonyms=['human'], reference_assemblies={'GRCh38':(76, 85), 'GRCh37': (55, 75), 'NCBI36': (54, 54)}), server='ftp://ftp.ensembl.org')

def usage():
	print('Usage:\n')
	print('\tpython '+sys.argv[0]+' -f <metadata.txt> -d <results directory> -k <postprocess directory> -l <locus file> -s <legitimate file> -m <mark> -g <genome type> -v <species name> -a <release number> [-n <marks input> -z <pool size> -m <min gap> -r <repeatMasker file> -t <duplicate locus file> -c <contruction fasta file> -j <reference genome fasta file>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-f or --metadata : metadata file')
	print('\t\t-d or --res_directory : results directory')
	print('\t\t-k or --pos_directory : postprocess directory')
	print('\t\t-l or --file_locus : locus bed file (locus, start, end, strand, flag)')
	print('\t\t-s or --file_legitimate : legitimate positions (bait, prey, start, end, flag)')
	print('\t\t-m or --outputMark : mark added to the output file')
	print('\t\t-g or --genome : only filter librairies results with this genome')
	print('\t\t-v or --species : species in latin known in Ensembl database (homo_sapiens, mus_musculus...)')
	print('\t\t-a or --release : release number of your species known in Ensembl database (87,86...)')
	print('\t\t-z or --size_pool : the number of bases between two junctions to pool them for illegitimate junctions (Default : 100)')
	print('\t\t-p or --min_gap : minimum of bases between bait and prey to call junction (Default : 5)')
	print('\t\t-n or --inputMark : marks from input file')
	print('\t\t-r or --file_repeat : repeatMasker library (.csv)')	
	print('\t\t-t or --translate : duplicate locus file for modified genome, this option will translate junctions from donor site to acceptor site')
	print('\t\t-c or --construction : construction fasta file for modified genome (in $BOWTIE2_INDEXES), this option will ajust genes position to modified genome')
	print('\t\t-j or --fasta : reference genome fasta file (mandatory with -c option)')

def main(argv):

	file_metadata = ""
	path_res_directory = ""
	path_pos_directory = ""
	file_locus = ""
	file_legitimate = ""
	mark=""
	genome = ""
	species = ""
	release = ""
	size_pool = 100
	min_gap = 5
	input_marks = ""
	file_input_extension = ""
	file_output_extension = ""
	file_repeat = ""
	file_translate = ""
	file_construction = ""
	file_fasta = ""
	chrom_dict={}
	df_construction=pd.DataFrame(columns=["chr", "modif_type", "length", "start", "end", "seq"])
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'f:d:k:l:s:m:g:v:a:z:p:n:r:t:c:j:', ['metadata=', 'res_directory=', 'pos_directory=', 'file_locus=', 'file_legitimate=', 'outputMark=', 'genome=', 'species=', 'release=', 'size_pool=', 'min_gap=', 'inputMark=', 'file_repeat=', 'translate=', 'construction=', 'fasta=', 'help'])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

###################################################################################OPTIONS###################################################################################

	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit(2)
		elif opt in ('-f', '--metadata'):
			file_metadata = arg
		elif opt in ('-d', '--res_directory'):
			path_res_directory = arg
		elif opt in ('-k', '--pos_directory'):
			path_pos_directory = arg
		elif opt in ('-l', '--file_locus'):
			file_locus = arg
		elif opt in ('-s', '--file_legitimate'):
			file_legitimate = arg
		elif opt in ('-m', '--outputMark'):
			mark = arg
		elif opt in ('-g', '--genome'):
			genome = arg
		elif opt in ('-v', '--species'):
			species = arg
		elif opt in ('-a', '--release'):
			release = arg
		elif opt in ('-z', '--size_pool'):
			size_pool = arg
		elif opt in ('-p', '--min_gap'):
			min_gap = arg
		elif opt in ('-n', '--inputMark'):
			input_marks = arg
		elif opt in ('-r', '--file_repeat'):
			file_repeat = arg
		elif opt in ('-t', '--translate'):
			file_translate = arg
		elif opt in ('-c', '--construction'):
			file_construction = arg
		elif opt in ('-j', '--fasta'):
			file_fasta = arg
		else:
			print("Error : Bad option -> "+opt)
			usage()
			sys.exit(2)

###################################################################################CHECK UP/SET UP###################################################################################

###CHECK METADATA FILE
	if file_metadata=="" or not os.path.exists(file_metadata):
		print("Error : You have to set a metadata file !\n")
		usage()
		sys.exit(2)
	else:
		###READ METADATA FILE
		metadata = pd.read_table(file_metadata, sep='\t')

###CHECK RESULTS DIRECTORY
	if not os.path.exists(path_res_directory):
		print("Error : You have to set a results directory !\n")
		usage()
		sys.exit(2)
	else:
		if path_res_directory[-1] != "/":
			path_res_directory+="/"

###CHECK POSTPROCESS DIRECTORY
	if not os.path.exists(path_pos_directory):
		print("Error : You have to set a postprocess directory !\n")
		usage()
		sys.exit(2)
	else:
		if path_pos_directory[-1] != "/":
			path_pos_directory+="/"

###CHECK LOCUS FILE
	if file_locus=="" or not os.path.exists(file_locus):
		print("Error : You have to set a locus file !\n")
		usage()
		sys.exit(2)
	else:
		###READ LOCUS FILE
		df_locus = pd.read_table(file_locus, sep='\t', header=None)
		for index, row in df_locus.iterrows():
			try :
				if int(row[1]) < 1:
					print("Error : line."+str(index+1)+", col.2 of your locus file !\n")
					print("Error : Start position has to be positive integer !\n")
					usage()
					sys.exit(2)
			except :
				print("Error : line."+str(index+1)+", col.2 of your locus file !\n")
				print("Error : Unknown start position : "+str(row[1])+" !\n")
				usage()
				sys.exit(2)
			try : 
				if int(row[2]) < 1:
					print("Error : line."+str(index+1)+", col.3 of your locus file !\n")
					print("Error : End position has to be positive integer !\n")
					usage()
					sys.exit(2)
			except :
				print("Error : line."+str(index+1)+", col.3 of your locus file !\n")
				print("Error : Unknown end position : "+str(row[2])+" !\n")
				usage()
				sys.exit(2)
			if int(row[2]) < int(row[1]):
				print("Error : line."+str(index+1)+" of your locus file !\n")
				print("Error : End position is smaller than start position !\n")
				usage()
				sys.exit(2)
			if row[3] != '+' and row[3] != '-':
				print("Error : line."+str(index+1)+", col.4 of your locus file !\n")
				print("Error : Unknown strand (+,-) : "+str(row[3])+" !\n")
				usage()
				sys.exit(2)

###CHECK LEGITIME FILE
	if file_legitimate=="" or not os.path.exists(file_legitimate):
		print("Error : You have to set a legitimate file !\n")
		usage()
		sys.exit(2)
	else:
		###READ LEGITIME FILE
		df_legitimate_locus = pd.read_table(file_legitimate, sep='\t', header=None)
		for index, row in df_legitimate_locus.iterrows():
			if row[0][0:3] != 'chr':
				print("Error : line."+str(index+1)+", col.1 of your legitimate locus file !\n")
				print("Error : Unknown chromosome : "+row[0]+" !\n")
				usage()
				sys.exit(2)
			if row[1][0:3] != 'chr':
				print("Error : line."+str(index+1)+", col.2 of your legitimate locus file !\n")
				print("Error : Unknown chromosome : "+row[1]+" !\n")
				usage()
				sys.exit(2)
			try :
				if int(row[2]) < 1:
					print("Error : line."+str(index+1)+", col.3 of your legitimate locus file !\n")
					print("Error : Start position has to be positive integer !\n")
					usage()
					sys.exit(2)
			except :
				print("Error : line."+str(index+1)+", col.3 of your legitimate locus file !\n")
				print("Error : Unknown start position : "+str(row[2])+" !\n")
				usage()
				sys.exit(2)
			
			try :
				if int(row[3]) < 1:
					print("Error : line."+str(index+1)+", col.4 of your legitimate locus file !\n")
					print("Error : End position has to be positive integer !\n")
					usage()
					sys.exit(2)
			except :
				print("Error : line."+str(index+1)+", col.4 of your legitimate locus file !\n")
				print("Error : Unknown end position : "+str(row[3])+" !\n")
				usage()
				sys.exit(2)

			if int(row[3]) < int(row[2]):
				print("Error : line."+str(index+1)+" of your legitimate locus file !\n")
				print("Error : End position is smaller than start position !\n")
				usage()
				sys.exit(2)

###FILTER METADATA FILE
	if genome != "":
		metadata = metadata.loc[metadata['Assembly'] == genome]
		if metadata.empty:
			print("Error : This assembly does not exist in the metadata file !\n")
			usage()
			sys.exit(2)
	else:
		print("Error : You have to set a genome name !\n")
		usage()
		sys.exit(2)

###CHECK SIZE POOL
	try:
		size_pool = int(size_pool)
		if size_pool < 1:
			print("Error : Pool size option needs to be more than 0 !\n")
			usage()
			sys.exit(2)
	except:
		print("Error : You have to set an integer to pool size option !\n")
		usage()
		sys.exit(2)

###CHECK MINIMUM GAP
	try:
		min_gap = int(min_gap)
		if int(min_gap) < 1:
			print("Error : Minimum gap option needs to be more than 0 !\n")
			usage()
			sys.exit(2)
	except:
		print("Error : You have to set an integer to minimum gap option !\n")
		usage()
		sys.exit(2)

###CHECK SPECIES AND RELEASE
	if species != "" and release != "" :
		try:
			data_ensembl = pe.EnsemblRelease(release, species=species) # will take a while to get the genome data installed in your system
			#data_ensembl=pe.EnsemblRelease(release=85, species=pe.Species(latin_name='homo_sapiens', synonyms=['human'], reference_assemblies={'GRCh38':(76, 85), 'GRCh37': (55, 75), 'NCBI36': (54, 54)}), server='ftp://ftp.ensembl.org')
		except:
			print("Error : Species name and/or release number not correct !\n")
			usage()
			sys.exit(2)
	else:
		if species == "":
			print("Error : You have to set a species name !\n")
			usage()
			sys.exit(2)
		if release == "":
			print("Error : You have to set a release number according to your species name !\n")
			usage()
			sys.exit(2)
	
###CHECK INPUT MARKS HISTORY
	if input_marks == "":
		print("Warning : You will process the raw file !\n")

###CHECK MARK
	if mark=="" :
		print("Error : You have to set a mark !\n")
		usage()
		sys.exit(2)

###SELECT INPUT FILES
	if input_marks == "":
		file_input_extension = ".tlx"
	else:
		file_input_extension = "_"+"_".join(input_marks.split(","))+".tlx"

###TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
	check_inputMark = False
	for i in metadata['Library'].tolist():	
		if os.path.exists(path_res_directory+i+"/"+i+file_input_extension):
			check_inputMark = True
	if not check_inputMark:
		print("Error : Your input marks can not localize a good input file !\n")
		usage()
		sys.exit(2)

###SELECT OUTPUT FILES
	if file_input_extension != "":
		file_output_extension = file_input_extension[:-4]+"_"+mark+".tlx"
	else:
		file_output_extension = "_"+mark+".tlx"

###CHECK REPEAT FILE
	if file_repeat != "":
		if not os.path.exists(file_repeat):
			print("Error : You activate -r option but the repeatMasker file is missing !\n")
			usage()
			sys.exit(2)
		else:
			f_repeat = open(file_repeat, 'rt')
			try:
				df_repeat = pd.read_csv(f_repeat, sep='\t', header=0, index_col=None)
			finally:
				f_repeat.close()

###CHECK TRANSLATE FILE
	if file_translate != "":
		if not os.path.exists(file_translate):
			print("Error : You activate -t option but the duplicate construction file is missing !\n")
			usage()
			sys.exit(2)
		else:
			###READ CONSTRUCTION FILE
			duplicate_locus = pd.read_table(file_translate, sep='\t', header=None)
			for index, row in duplicate_locus.iterrows():
				###DONOR CHECK
				if row[0][0:3] != 'chr':
					print("Error : line."+str(index+1)+", col.1 of your construction file !\n")
					print("Error : Unknown chromosome : "+row[0]+" !\n")
					usage()
					sys.exit(2)
				try :
					if int(row[1]) < 1:
						print("Error : line."+str(index+1)+", col.2 of your construction file !\n")
						print("Error : Start position has to be positive integer !\n")
						usage()
						sys.exit(2)
				except :
					print("Error : line."+str(index+1)+", col.2 of your construction file !\n")
					print("Error : Unknown start position : "+str(row[1])+" !\n")
					usage()
					sys.exit(2)
				try :
					if int(row[2]) < 1:
						print("Error : line."+str(index+1)+", col.3 of your construction file !\n")
						print("Error : End position has to be positive integer !\n")
						usage()
						sys.exit(2)
				except :
					print("Error : line."+str(index+1)+", col.3 of your construction file !\n")
					print("Error : Unknown end position : "+str(row[2])+" !\n")
					usage()
					sys.exit(2)
				if int(row[2]) < int(row[1]):
					print("Error : line."+str(index+1)+" of your construction file !\n")
					print("Error : End position is smaller than start position !\n")
					usage()
					sys.exit(2)
				if row[3] != '+' and row[3] != '-':
					print("Error : line."+str(index+1)+", col.4 of your construction file !\n")
					print("Error : Unknown strand (+,-) : "+str(row[3])+" !\n")
					usage()
					sys.exit(2)
				###ACCEPTOR CHECK
				if row[4][0:3] != 'chr':
					print("Error : line."+str(index+1)+", col.5 of your construction file !\n")
					print("Error : Unknown chromosome : "+row[4]+" !\n")
					usage()
					sys.exit(2)
				try :
					if int(row[5]) < 1:
						print("Error : line."+str(index+1)+", col.6 of your construction file !\n")
						print("Error : Start position has to be positive integer !\n")
						usage()
						sys.exit(2)
				except :
					print("Error : line."+str(index+1)+", col.6 of your construction file !\n")
					print("Error : Unknown start position : "+str(row[5])+" !\n")
					usage()
					sys.exit(2)
				try :
					if int(row[6]) < 1:
						print("Error : line."+str(index+1)+", col.7 of your construction file !\n")
						print("Error : End position has to be positive integer !\n")
						usage()
						sys.exit(2)
				except :
					print("Error : line."+str(index+1)+", col.7 of your construction file !\n")
					print("Error : Unknown end position : "+str(row[6])+" !\n")
					usage()
					sys.exit(2)
				if int(row[6]) < int(row[5]):
					print("Error : line."+str(index+1)+" of your construction file !\n")
					print("Error : End position is smaller than start position !\n")
					usage()
					sys.exit(2)
				if row[7] != '+' and row[7] != '-':
					print("Error : line."+str(index+1)+", col.8 of your construction file !\n")
					print("Error : Unknown strand (+,-) : "+str(row[7])+" !\n")
					usage()
					sys.exit(2)

	###IF CONSTRUCTION OPTION IS SET
	if file_construction != "":
		if not os.path.exists(file_construction):
			print("Error : You activate -c option but the construction file is missing !\n")
			usage()
			sys.exit(2)
		elif (file_fasta[-3:]!=".fa" and file_fasta[-4:]!=".fna" and file_fasta[-6:]!=".fasta") or not os.path.exists(file_fasta):
			print("Error : You activate -c option but the fasta file is missing (-j option) !\n")
			usage()
			sys.exit(2)
		else:
			###DICTIONNARY ALL CHROMOSOMES AND LENGTHS
			chrom_dict = {
				record.id: len(record.seq)
				for record in SeqIO.parse(file_fasta, 'fasta')
				if record.id not in chrom_dict
			}

			###CHECK CONSTRUCTION FILE
			removed_list=[]
			added_dict={}
			
			###FILTER CHANGES_FILE, KEEP GOOD DELETION AND INSERTION. DISPATCH ADD AND REMOVE IN LIST
			for record in SeqIO.parse(file_construction, 'fasta'):
				modif_type=""
				chrom=""
				start=0
				end=0
				check_good=True
				if len(record.id.split(":")) == 2:
					check_good=False
					if record.id.split(":")[0] == "add":
						if record.id.split(":")[0] not in added_dict:
							added_dict[record.id.split(":")[1]]=record.seq
					elif record.id.split(":")[0] == "remove":
						removed_list.append(record.id.split(":")[1])
			
				elif len(record.id.split(":")) == 4:
					modif_type=record.id.split(":")[0]
					chrom=record.id.split(":")[1]
					try:
						start=int(record.id.split(":")[2])
					except ValueError:
						check_good=False
						print("Warning : {"+record.id+"} start need to be integer, please use this pattern modif_type:chrmX:start:end !")
						print("Warning : {"+record.id+"} will not be used")
					try:
						end=int(record.id.split(":")[3])
					except ValueError:
						check_good=False
						print("Warning : {"+record.id+"} end need to be integer, please use this pattern modif_type:chrmX:start:end !")
						print("Warning : {"+record.id+"} will not be used")
			
					if modif_type != "insertion" and modif_type != "deletion":
						check_good=False
						print("Warning : {"+record.id+"} does not have the good format, please use modif_type:chrmX:start:end !")
						print("Warning : To delete or insert use modif_type:chrmX:start:end !")
						print("Warning : To add or remove use modif_type:chrmX !")
						print("Warning : {"+record.id+"} will not be used")

					if modif_type == "insertion" and record.seq=="":
						check_good=False
						print("Warning : {"+record.id+"} does not have a sequence to insert")
						print("Warning : {"+record.id+"} will not be used")

					if chrom not in chrom_dict:
						check_good=False
						print("Warning : {"+record.id+"} chrom is not in the reference genome file, please use this pattern modif_type:chrmX:start:end !")
						print("Warning : {"+record.id+"} will not be used")
				
					if start < 0 or end < 1:
						check_good=False
						print("Warning : {"+record.id+"} start < 0 or end position < 1 is not allowed, please use this pattern modif_type:chrmX:start:end !")
						print("Warning : {"+record.id+"} will not be used")
					if end < start:
						check_good=False
						print("Warning : {"+record.id+"} end < start is not allowed, please use this pattern modif_type:chrmX:start:end !")
						print("Warning : {"+record.id+"} will not be used")
		
				else:
					check_good=False
					print("Warning : {"+record.id+"} does not have the good format")
					print("Warning : To delete or insert use modif_type:chrmX:start:end !")
					print("Warning : To add or remove use modif_type:chrmX !")
					print("Warning : {"+record.id+"} will not be used")

				if check_good:
					###FOR INSERTION DO THE FOLLOWING			
					if modif_type == "insertion":
						if start == end:
							print("Warning : {"+record.id+"} can not insert between "+str(start)+" and "+str(end)+" !")
							print("Warning : {"+record.id+"} insertion between "+str(start)+" and "+str(end+1)+" !")
							end+=1
						if start > chrom_dict[chrom]:
								print("Warning : {"+record.id+"} start > chromosome length !")
								print("Warning : {"+record.id+"} insertion at the end !")
								start=chrom_dict[chrom]
								end=chrom_dict[chrom]+1
						if end > chrom_dict[chrom]+1:
								print("Warning : {"+record.id+"} end > chromosome length +1 !")
								print("Warning : {"+record.id+"} insertion at the end !")
								start=chrom_dict[chrom]
								end=chrom_dict[chrom]+1

						df_construction = df_construction.append(pd.Series([chrom, modif_type, len(record.seq)-(end-start-1), start, end, str(record.seq)], index=["chr", "modif_type", "length", "start", "end", "seq"]), ignore_index=True)
					###FOR DELETION DO THE FOLLOWING
					elif modif_type == "deletion":
						if start == end or start == end-1:
							print("Warning : {"+record.id+"} can not delete between "+str(start)+" and "+str(end)+", nothing to delete !")
					
						elif start >= chrom_dict[chrom] and end > chrom_dict[chrom]:
							print("Warning : {"+record.id+"} can not delete between "+str(start)+" and "+str(end)+", nothing to delete !")
						else:
							if start < chrom_dict[chrom] and end > chrom_dict[chrom]:
								print("Warning : {"+record.id+"} end > chromosome length +1 !")
								print("Warning : {"+record.id+"} deletion from "+str(start)+" to the end !")
								end=chrom_dict[chrom]+1

							df_construction = df_construction.append(pd.Series([chrom, modif_type, -(end-start-1), start, end, ""], index=["chr", "modif_type", "length", "start", "end", "seq"]), ignore_index=True)
					else:
						print("Warning : {"+record.id+"} does not have the good format")
						print("Warning : To delete or insert use modif_type:chrmX:start:end !")
						print("Warning : To add or remove use modif_type:chrmX !")
						print("Warning : {"+record.id+"} will not be used")
						sys.exit()

			###SORT DATAFRAME
			df_construction = df_construction.sort_values(by=['chr', 'start', 'end'], ascending=[1, 1, 1])

###################################################################################PRINTS###################################################################################

	print('\n-----------------------------------------')
	print('Metadata File tlx : '+file_metadata)
	print('Results Directory : '+path_res_directory)
	print('Postprocess Directory : '+path_pos_directory)
	print('Locus File bed : '+file_locus)
	print('Legitimate locus File : '+file_legitimate)
	print('Genome : '+genome)
	print('Species : '+species)
	print('Size pool : '+str(size_pool))
	print('Minimum gap : '+str(min_gap))
	print('Input file extension: '+file_input_extension)
	print('Output file extension : '+file_output_extension)
	if not file_repeat == "":
		print('RepeatMasker file : '+file_repeat)
	else:
		print('No repeatMasker information')
	if not file_translate == "":
		print('Duplicate locus file : '+file_translate)
		print('\tWarning : Some junctions will be change')
	else:
		print('No junction translation file needed')
	if not file_construction == "":
		print('Construction file : '+file_construction)
	else:
		print('No construction file needed')
	if not file_fasta == "":
		print('Reference genome file : '+file_fasta)
	else:
		print('No reference genome file needed')
	print('-----------------------------------------\n')

###################################################################################PROGRAM###################################################################################

	##WRITE DATAFRAME ON ONE LINE
	pd.set_option('expand_frame_repr', False)

	df = pd.DataFrame()
	dict_position=[]
	table_aux=[]
	check_legitimate=False
	check_trash=False
	max_type=""
	count_legitimate = 0
	count_trash = 0
	count_illegitimate = 0

	###LOOP OVER EACH LIBRARIES
	for library in metadata['Library'].tolist():
		###CHECK DIRECTORY EXISTS
		if not os.path.exists(path_res_directory+library):
			print("Warning : "+path_res_directory+" does not contains {"+library+"}")
			print("Warning :  {"+library+"} will not be filtered")
		else:
			count_junction = 0
			###CHECK INPUT FILE EXISTS
			if os.path.exists(path_res_directory+library+"/"+library+file_input_extension):
				print(library)
				if not os.path.exists(path_pos_directory+library):
					os.system("mkdir "+path_pos_directory+library)
				###WRITE HEADER LEGITIME
				try:
					with open(path_pos_directory+library+"/"+library+"_Legitimate"+file_output_extension, 'wb') as csvfile:
						spamwriter = csv.writer(csvfile, delimiter='\t')
						spamwriter.writerow(['Qname', 'JuncID', 'Rname', 'Junction', 'Strand', 'Rstart', 'Rend', 'B_Rname', 'B_Rstart', 'B_Rend', 'B_Strand', 'Locus Type', 'Position Type', 'Position Delta', 'Seq', 'LenSeq'])
				finally:
					csvfile.close()

				###WRITE HEADER TRASH
				try:
					with open(path_pos_directory+library+"/"+library+"_Trash"+file_output_extension, 'wb') as csvfile:
						spamwriter = csv.writer(csvfile, delimiter='\t')
						spamwriter.writerow(['Qname', 'JuncID', 'Rname', 'Junction', 'Strand', 'Rstart', 'Rend', 'B_Rname', 'B_Rstart', 'B_Rend', 'B_Strand', 'Locus Type', 'Position Type', 'Position Delta', 'Seq', 'LenSeq'])
				finally:
					csvfile.close()

				df = pd.read_csv(path_res_directory+library+"/"+library+file_input_extension, sep='\t', header=0, index_col=None)
				
				###CREATE ILLEGITIME TABLE
				df_illegitimates = pd.DataFrame(columns=list(df.columns.values))

				df = df.sort_values('Qname', ascending=False)

				###FOR EACH JUNCTION OF TLX FILE
				for index, row in df.iterrows():
					count_junction+=1
					#print("--------------------------------------------------------------------")
					#print(row)
					current_junction = row
					Flag = ''
					###FOCUS ON PRIMARY ASSEMBLY BECAUSE VISUALIZATION ALLOWS ONLY PRIMARY ASSEMBLY (REMOVE THIS LINE IF NECESSARY)
					if len(row["Rname"]) <= 5 and len(row["B_Rname"]) <= 5 and row["Rname"] != "chrM" and row["B_Rname"] != "chrM":
						###DODGE GENOMIC BACKGROUD B6/SV129dodge genomic background B6/SV129
						if min_gap != 0 :
							if row["Rname"] == row["B_Rname"] :
								if row["Rend"] <= row["B_Rstart"]:
									if row["B_Rstart"]-row["Rend"] <= min_gap:
										check_legitimate=False
										check_trash=True
										max_type = "None"
										max_position = "None"
										best_locus = "B6/SV129_construct"
										#print("Remove B6/SV129_construct")
								elif row["B_Rend"] <= row["Rstart"]:
									if row["Rstart"]-row["B_Rend"] <= min_gap:
										check_legitimate=False
										check_trash=True
										max_type = "None"
										max_position = "None"
										best_locus = "B6/SV129_construct"
										#print("Remove B6/SV129_construct")
						if not check_trash:
							###SEND IT TO LEGITIME OR ILLEGITIME TABLE
							for index_locus, row_locus in df_legitimate_locus.iterrows():
								###READ EACH LINE OF LOCUS LEGITIME
								if row["B_Rname"] == row_locus[0] and row["Rname"] == row_locus[1] and int(row["Junction"]) >= int(row_locus[2]) and int(row["Junction"]) <= int(row_locus[3]):
									###LEGITIME JUNCTION
									check_legitimate=True
									basic_locus_name = row_locus[4]
									if file_translate != "" :
										###TRANSLATE JUNCTION USING DUPLICATE AREA
										for index_dup_locus, row_dup_locus in duplicate_locus.iterrows():
											if row_dup_locus[0] == current_junction["Rname"] and int(current_junction["Junction"]) >= int(row_dup_locus[1]) and int(current_junction["Junction"]) <= int(row_dup_locus[2]):
												distance_from_start = int(current_junction["Junction"]) - int(row_dup_locus[1])
												distance_end_start = int(current_junction["Rend"]) - int(current_junction["Rstart"])

												current_junction["Junction"] = int(row_dup_locus[5]) + distance_from_start
												current_junction["Rname"] = row_dup_locus[4]
												if int(current_junction["Strand"]) == 1:
													current_junction["Rstart"] = current_junction["Junction"]
													current_junction["Rend"] = current_junction["Junction"] + distance_end_start

												else:
													current_junction["Rend"] = current_junction["Junction"]
													current_junction["Rstart"] = current_junction["Junction"] - distance_end_start
												if current_junction["Rend"] > int(row_dup_locus[6]):
													current_junction["Rend"] = int(row_dup_locus[6])
													if int(current_junction["Strand"]) == -1:
														current_junction["Junction"] = current_junction["Rend"]
										#print(current_junction)
								
									###FIND THE BEST LOCUS
									max_type, max_position, best_locus = getLocus(current_junction, basic_locus_name, df_locus)
									#print(str(max_type), str(max_position), best_locus)
										

						#if not check_legitimate:
						#	print("-------------------------------ILLEGITIME")
					else:
						#print("-------------------------------TRASH")
						check_trash=True
						check_legitimate=False
						max_type = "None"
						max_position = "None"
						best_locus = "Not Primary"

					###IF NOT LEGITIME NEITHER TRASH
					if not check_legitimate and not check_trash:
						count_illegitimate += 1
						df_illegitimates = df_illegitimates.append(row, ignore_index=True)

					###Write for locus file
					if check_legitimate and not check_trash:
						count_legitimate+=1
						try:
							with open(path_pos_directory+library+"/"+library+"_Legitimate"+file_output_extension, 'a') as csvfile:
								spamwriter = csv.writer(csvfile, delimiter='\t')
								table=np.array(current_junction, dtype=pd.Series)[0:11].tolist()
								table.append(best_locus)
								table.append(max_type)
								table.append(max_position)
								table.append(np.array(current_junction, dtype=pd.Series)[18])
								table.append(int(len(np.array(current_junction, dtype=pd.Series)[18])))
					
								spamwriter.writerow(table)
						finally:
							csvfile.close()

					###Write for trash file
					if check_trash:
						count_trash+=1
						try:
							with open(path_pos_directory+library+"/"+library+"_Trash"+file_output_extension, 'a') as csvfile:
								spamwriter = csv.writer(csvfile, delimiter='\t')
								table=np.array(current_junction, dtype=pd.Series)[0:11].tolist()
								table.append(best_locus)
								table.append(max_type)
								table.append(max_position)
								table.append(np.array(current_junction, dtype=pd.Series)[18])
								table.append(int(len(np.array(current_junction, dtype=pd.Series)[18])))
								spamwriter.writerow(table)
						finally:
							csvfile.close()

					best_locus=""
					max_type=""
					max_position=0
					check_legitimate=False
					check_trash=False

				###Fill Illegtime file with repeatMasker informations
				res=""
				repeat_array=[]
				for index, row in df_illegitimates.iterrows():
					###CHANGE JUNCTION POSITION IF MODIFIED GENOME
					if not df_construction.empty:
						junc_start, junc_end = modified_genome_position(df_construction, row['Rname'], row['Junction'], row['Junction'])
					else:
						junc_start = row['Junction']
						junc_end = row['Junction']
					if junc_start != -1 and junc_end != -1:
						res = df_repeat.loc[(df_repeat['begin'].astype(int) <= junc_start) & (df_repeat['end'].astype(int) >= junc_start) & (df_repeat['sequence'] == row['Rname'])]
					if len(res) != 0:
						repeat_array.append(res["class/family"].values[0])
					else:
						repeat_array.append("None")
				df_illegitimates.insert(11, "RepeatEvent", repeat_array)

				df_illegitimates.to_csv(path_pos_directory+library+"/"+library+"_Illegitimate"+file_output_extension, sep='\t', float_format='%.0f', index=False)

				###Illegitimates management
				#print(df_illegitimates)

				dict_position = listChromIllegitimate(df_illegitimates)
				#print(dict_position)
				binIllegitimate(dict_position, path_pos_directory+library+"/"+library+"_Percent_illegitimate"+file_output_extension, len(df_illegitimates), size_pool, data_ensembl, df_construction)

			else:
				print("Warning : "+path_res_directory+library+"/"+library+file_input_extension+" does not exist")
				print("Warning :  {"+path_res_directory+library+"/"+library+file_input_extension+"} will not be filtered")

		print("Number of junctions : "+str(count_junction))
		print("Number of legitimate junctions : "+str(count_legitimate))
		print("Number of illegitimate junctions : "+str(count_illegitimate))
		print("Number of trash junctions : "+str(count_trash))
		#sys.exit()

	#print("DONE")

##################################################################################FUNCTIONS##################################################################################

###Define the best locus or the closest for a junction
def getLocus(line, basic_locus_name, df_locus):

	dic_bed={}
	df_accurate_chrom = df_locus.loc[df_locus[0] == line['Rname']]
	for index, row in df_accurate_chrom.iterrows():
		dic_bed[row[4]]=str(row[1])+"-"+str(row[2])

	best_locus=basic_locus_name
	max_type=""
	max_position=0
	break_bait_position = line["Junction"]
	#print(break_bait_position)
	for key,value in dic_bed.items():
		#print("******************")
		
		locusS=int(value.split("-")[0])
		locusE=int(value.split("-")[1])
		#print(key)
		#print(locusS)
		#print(locusE)
		#print(key,locusS, locusE)

		if break_bait_position<locusS :
			if max_type=="" or (locusS-break_bait_position)<max_position :
				#print("5'Outside")
				max_type="5'Outside"
				max_position=int(locusS-break_bait_position)
				#print(max_position)
				best_locus=key
			#else:
			#	print("5'Outside not better")
		elif break_bait_position>locusE :
			if max_type=="" or (break_bait_position-locusE)<max_position :
				#print("3'Outside")
				max_type="3'Outside"
				max_position=int(break_bait_position-locusE)
				#print(max_position)
				best_locus=key
			#else:
			#	print("3'Outside not better")

		elif break_bait_position>=locusS and break_bait_position<=locusE:
			###Manage inclusive locus
			#print("Inside")
			max_type="Inside"
			max_position=0
			best_locus=key
			break
		#else:
		#	print("Impossible")
	#print(max_type, max_position, best_locus)
	return max_type, max_position, best_locus

#####################ILLEGITIMES

def listChromIllegitimate(df):
	dict_position={}
	for index, row in df.iterrows():
		if row["Rname"] not in dict_position:
			dict_position[str(row["Rname"])] = {}
		dict_position[row["Rname"]][int(row["Junction"])] = 0
	return dict_position

#ENABLE ITERATION ON DATAFRAME
def iterate(iterable):
    iterator = iter(iterable)
    item = iterator.next()

    for next_item in iterator:
        yield item, next_item
        item = next_item

    yield item, None

def modified_genome_position(df_construction, chrom, start, end):
	aux_start=0
	aux_end=0
	for index, row in df_construction.loc[df_construction['chr'] == chrom].iterrows():
		if not int(row['start']) > int(end):
			if int(start) > int(row['start']) and int(end) < int(row['end']):
				return -1,-1
			aux_start = int(aux_start) + int(row['length'])
			aux_end = int(aux_end) + int(row['length'])
	return int(start)-int(aux_start),int(end)-int(aux_end)	

def binIllegitimate(dict_position, output, nb_illegitimates, size_pool, data_ensembl, df_construction):
	
	table_aux=[]

	###SORT DICTIONNARY
	dict_position = collections.OrderedDict(sorted(dict_position.items()))

	###SORT SUBDICTIONNARY
	for key, value in dict_position.iteritems():
		dict_position[key] = collections.OrderedDict(sorted(dict_position[key].items()))

	if len(dict_position)!=0:
		for key, value in dict_position.iteritems():
			dic_size_pool={}
			nb_transloc=1
			min_pos=0
			j=0
			for key_sub, next_key_sub in iterate(dict_position[key].iteritems()):
				if j < len(dict_position[key]):
					if j == len(dict_position[key])-1:
						#print("LAST ONE")
						if min_pos==0:
							#print("MIN_POS -> 0")
							min_pos=key_sub[0]
							#print("MIN_POS -> "+str(min_pos))
						dic_size_pool[str(min_pos)+"-"+str(key_sub[0])]=nb_transloc
						#print("DIC SIZE POOL")
						#print(dic_size_pool)
						min_pos=key_sub[0]
						#print("NEW MIN_POS -> "+str(min_pos))
						nb_transloc=1
					else:
						#print("NOT THE LAST ONE")
						if int(key_sub[0])+size_pool<int(next_key_sub[0]):
							#print("CAN'T GRAB FURTHER")
							if min_pos==0:
								#print("MIN_POS -> 0")
								min_pos=key_sub[0]
								#print("MIN_POS -> "+str(min_pos))
							dic_size_pool[str(min_pos)+"-"+str(key_sub[0])]=nb_transloc
							#print("DIC SIZE POOL")
							#print(dic_size_pool)
							min_pos=next_key_sub[0]
							#print("NEW MIN_POS -> "+str(min_pos))
							nb_transloc=1
						else:
							#print("+1")
							nb_transloc+=1

					j+=1
			if len(dic_size_pool) != 0:
				with open(output, 'a') as csvfile:
					spamwriter = csv.writer(csvfile, delimiter='\t')
					#print(data_ensembl.genes(contig=i[0][3:]))
					for key_print, value_print in dic_size_pool.iteritems():
						###TOO HARD TO MODIFY ENSEMBL COORDINATES TO FIT MY MODIFIED GENOME, SO I WILL MODIFY POSITION TO FIT ENSEMBL COORDINATES
						if not df_construction.empty:
							start, end = modified_genome_position(df_construction, key, key_print.split("-")[0], key_print.split("-")[1])
						else:
							start = key_print.split("-")[0]
							end = key_print.split("-")[1]

						gene_names = data_ensembl.gene_names_at_locus(contig=key[3:], position=int(start), end=int(end))
						genes= ",".join(gene_names)
						spamwriter.writerow([key,key_print.split("-")[0],key_print.split("-")[1],round(float(int(value_print)*100/float(nb_illegitimates)),4), genes])
			else:
				print("Nothing")

			#sys.exit(2)	

	else:
		print("Nothing")

if __name__ =='__main__':
	main(sys.argv[1:])
