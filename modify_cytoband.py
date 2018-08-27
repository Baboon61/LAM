### A script allowing to modify a cytoband file according the modification apply to a genome. Modification information saved in an other file

##########################################################################IMPORTS##########################################################################

import os
import sys
import getopt
import os.path
import numpy as np
import pandas as pd
from Bio import BiopythonWarning 
from Bio import SeqIO 
import pyensembl as pe
import csv
import collections

def usage():
	print('Usage:\n')
	print('\tpython '+sys.argv[0]+' -f <cytoband.txt> -g <genome type> -j <reference genome fasta file> [-c <contruction fasta file>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-f or --cytoband : cytoband file')
	print('\t\t-g or --genome : will add the genome to the cytoband file')
	print('\t\t-c or --construction : construction fasta file for modified genome (in $BOWTIE2_INDEXES)')
	print('\t\t-j or --fasta : reference genome fasta file')

def main(argv):

	file_cytoband = ""
	genome = ""
	file_construction = ""
	file_fasta = ""
	chrom_dict={}
	df_construction=pd.DataFrame(columns=["chr", "modif_type", "length", "start", "end", "seq"])
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'f:g:c:j:', ['cytoband=', 'genome=', 'construction=', 'fasta=', 'help'])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

###################################################################################OPTIONS###################################################################################

	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit(2)
		elif opt in ('-f', '--cytoband'):
			file_cytoband = arg
		elif opt in ('-g', '--genome'):
			genome = arg
		elif opt in ('-c', '--construction'):
			file_construction = arg
		elif opt in ('-j', '--fasta'):
			file_fasta = arg
		else:
			print("Error : Bad option -> "+opt)
			usage()
			sys.exit(2)

###################################################################################CHECK UP/SET UP###################################################################################

###CHECK CYTOBAND FILE
	if file_cytoband=="" or not os.path.exists(file_cytoband):
		print("Error : You have to set a cytoband file !\n")
		usage()
		sys.exit(2)
	else:
		###READ CYTOBAND FILE
		df_cytoband = pd.read_table(file_cytoband, sep='\t', names = ["Chromosome", "Start", "End", "Band", "Stain"])
		df_cytoband = df_cytoband.sort_values(['Chromosome', "Start"], ascending=[True, True])

###CHECK GENOME NAME
	if genome =="":
		print("Error : You have to set the name of the genome you are working on !\n")
		usage()
		sys.exit(2)

###CHECK REFERENCE GENOME FILE
	if (file_fasta[-3:]!=".fa" and file_fasta[-4:]!=".fna" and file_fasta[-6:]!=".fasta") or not os.path.exists(file_fasta):
		print("Error : You have to set a reference fasta file !\n")
		usage()
		sys.exit(2)
	else:
		###DICTIONNARY ALL CHROMOSOMES AND LENGTHS
		chrom_dict = {
			record.id: len(record.seq)
			for record in SeqIO.parse(file_fasta, 'fasta')
			if record.id not in chrom_dict
		}


###IF CONSTRUCTION OPTION IS SET
	if file_construction != "":
		if not os.path.exists(file_construction):
			print("Error : You activate -c option but the construction file is missing !\n")
			usage()
			sys.exit(2)
		else:
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
						print("{"+record.id+"} start need to be integer, please use this pattern modif_type:chrmX:start:end !")
						print("{"+record.id+"} will not be used")
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
	print('Cytoband File : '+file_cytoband)
	print('Genome : '+genome)
	if not file_construction == "":
		print('Construction file : '+file_construction)
	else:
		print('No construction file needed')
	if not file_fasta == "":
		print('Reference genome file : '+file_fasta)
	else:
		print('No reference genome file needed')
	print('Output file : '+genome+"_"+file_cytoband)
	print('-----------------------------------------\n')

###################################################################################PROGRAM###################################################################################

	##WRITE DATAFRAME ON ONE LINE
	pd.set_option('expand_frame_repr', False)

	pd.options.mode.chained_assignment = None  # default='warn'

	###CREATE AN EMPTY CYTOBAND MODIFIED DATAFRAME
	df_new_cytoband = pd.DataFrame(columns=list(df_cytoband))

	###RECOVER LIST OF MODIFIED CHROMOSOMES
	list_chr_construction = df_construction.chr.unique()

	modified_chromosome={}

	if file_construction!= "":
		
		####################################################
		###MANAGE ADD AND REMOVE CASES
		for key, value in added_dict.iteritems():
			df_new_cytoband = df_new_cytoband.append(pd.Series([key, 0, len(value), "Nan", "Nan"], index=list(df_new_cytoband.columns.values)), ignore_index=True)
			#print("ADD : "+str(key))
	
		for i in removed_list:
			df_cytoband = df_cytoband[df_cytoband.Chromosome != i]
			#print("REMOVE : "+str(i))

		###FOR EACH CHROMOSOME
		last_index_df_cyto = len(df_cytoband.index)+1

		for i in df_cytoband.Chromosome.unique():
			###IF IT IS IN FASTA FILE
			if i in list(chrom_dict.keys()) and i in list_chr_construction:
				###GET THE CONSTRUCTION FOR THIS CHROMOSOME
				df_construction_by_chr = df_construction.loc[df_construction['chr'] == i]
				###GET THE CYTOBAND FOR THIS CHROMOSOME
				df_cytoband_by_chr = df_cytoband.loc[df_cytoband['Chromosome'] == i]
				###CALL MODIFY FUNCTION
				#print("ANALYZE CHR : "+str(i))
				#print(df_cytoband_by_chr)
				modified_chromosome[i], length_df_cytoband = modify_cytoband_by_chr(df_cytoband_by_chr, df_construction_by_chr, last_index_df_cyto)

		#print(modified_chromosome)

		###FOR EACH CYTOBAND
		for index, row in df_cytoband.iterrows():
			###APPEND IF IT IS NOT MODIFIED
			if row["Chromosome"] not in list_chr_construction and row["Chromosome"] in list(chrom_dict.keys()):
				df_new_cytoband = df_new_cytoband.append(row)

		###FOR EACH MODIFIED CHROMOSOME
		for i in list_chr_construction:
			###APPEND EACH LINE MODIFIED
			for index, row in modified_chromosome[i].iterrows():
				df_new_cytoband = df_new_cytoband.append(row)
	else:
		for index, row in df_cytoband.iterrows():
			if row['Chromosome'] in list(chrom_dict.keys()):
				df_new_cytoband = df_new_cytoband.append(row)			

	###SORT ALL
	df_new_cytoband_int=pd.DataFrame(columns=list(df_new_cytoband.columns.values))
	df_new_cytoband_str=pd.DataFrame(columns=list(df_new_cytoband.columns.values))
	df_after_sort=pd.DataFrame(columns=list(df_new_cytoband.columns.values))
	df_new_cytoband['Chromosome'].replace(regex=True,inplace=True,to_replace=r'chr',value=r'')
	for index, row in df_new_cytoband.iterrows():
		try:
			int(row['Chromosome'])
			df_new_cytoband_int = df_new_cytoband_int.append(pd.Series(row, index=list(df_new_cytoband.columns.values)))
		except:
			df_new_cytoband_str = df_new_cytoband_str.append(pd.Series(row, index=list(df_new_cytoband.columns.values)))

	df_new_cytoband_int['Chromosome'] = df_new_cytoband_int['Chromosome'].astype(int)
	df_new_cytoband_int = df_new_cytoband_int.sort_values(by=['Chromosome', 'Start', 'End'], ascending=[1, 1, 1])
	df_new_cytoband_str = df_new_cytoband_str.sort_values(by=['Chromosome', 'Start', 'End'], ascending=[1, 1, 1])
	df_after_sort = pd.concat([df_new_cytoband_int,df_new_cytoband_str])
	df_after_sort['Chromosome'] = 'chr' + df_after_sort['Chromosome'].astype(str)
	df_after_sort.to_csv(genome+"_"+file_cytoband, sep='\t', float_format='%.0f', index=False, header=True)


##################################################################################FUNCTIONS##################################################################################

def modify_cytoband_by_chr(df_cyto, df_cons, last_index_df_cyto):
	#print(df_cyto)
	add_length=0

	for index, row in df_cons.iterrows():
		df_cons.at[index, 'start'] = int(row["start"])+add_length
		df_cons.at[index, 'end'] = int(row["end"])+add_length
		add_length=add_length+int(row["length"])

	#print(df_cons)
	for index_cons, row_cons in df_cons.iterrows():
		#check_modif=False
		to_modified_cyto = df_cyto[(df_cyto["End"] > row_cons["start"]) & (df_cyto["Start"] < row_cons["end"])]
		if row_cons["modif_type"] == "insertion":
			#print("INSERTION")
			if len(list(to_modified_cyto.index.values)) == 1:
				#print("SOLO")
				End = df_cyto.at[list(to_modified_cyto.index.values)[0],"End"]
				df_cyto.at[list(to_modified_cyto.index.values)[0],"End"] = int(row_cons["start"])
				###I removed the length of the construct because I add it later
				df_cyto.loc[last_index_df_cyto] = pd.Series([row_cons["chr"], (int(row_cons["start"])-int(row_cons["length"])), int(row_cons["end"]), "Nan", "Nan"], index=list(df_cyto.columns.values))
				last_index_df_cyto += 1
				df_cyto.loc[last_index_df_cyto] = pd.Series([df_cyto.at[list(to_modified_cyto.index.values)[0],"Chromosome"], int(row_cons["end"]), int(End), df_cyto.at[list(to_modified_cyto.index.values)[0],"Band"], df_cyto.at[list(to_modified_cyto.index.values)[0],"Stain"]], index=list(df_cyto.columns.values))
				last_index_df_cyto += 1

			elif len(list(to_modified_cyto.index.values)) > 1:
				#print("MULTIPLE")
				End = df_cyto.at[list(to_modified_cyto.index.values)[-1],"End"]
				#print(list(to_modified_cyto.index.values))
				df_cyto.at[list(to_modified_cyto.index.values)[0],"End"] = int(row_cons["start"])
				df_cyto.drop(list(to_modified_cyto.index.values)[1:-1], inplace = True)
				df_cyto.loc[last_index_df_cyto] = pd.Series([row_cons["chr"], (int(row_cons["start"])-int(row_cons["length"])), int(row_cons["end"]), "Nan", "Nan"], index=list(df_cyto.columns.values))
				last_index_df_cyto += 1
				df_cyto.at[list(to_modified_cyto.index.values)[-1],"Start"] = int(row_cons["end"])+int(row_cons["length"])
				df_cyto.at[list(to_modified_cyto.index.values)[-1],"End"] += int(row_cons["length"])
			else:
				print("Error in cytoband management for insertion")
				sys.exit()
				
			
		elif row_cons["modif_type"] == "deletion":
			#print("DELETION")
			if len(list(to_modified_cyto.index.values)) == 1:
				#print("SOLO")
				df_cyto.at[list(to_modified_cyto.index.values)[0],"End"] = int(df_cyto.at[list(to_modified_cyto.index.values)[0],"End"]) + int(row_cons["length"])

			elif len(list(to_modified_cyto.index.values)) > 1:
				#print("MULTIPLE")
				df_cyto.at[list(to_modified_cyto.index.values)[0],"End"] = int(row_cons["start"])
				df_cyto.drop(list(to_modified_cyto.index.values)[1:-1], inplace = True)
				df_cyto.at[list(to_modified_cyto.index.values)[-1],"Start"] = int(row_cons["start"])
				df_cyto.at[list(to_modified_cyto.index.values)[-1],"End"] += int(row_cons["length"])
			else:
				print("Error in cytoband management for deletion")
				sys.exit()
		else:
			print("Error in cytoband management not insertion or deletion")
			sys.exit()

		to_modified_rest = df_cyto[(df_cyto["End"] >= row_cons["start"])]
		for i in list(to_modified_rest.index.values):
			if i not in list(to_modified_cyto.index.values):
				df_cyto.at[i,"Start"] = int(df_cyto.at[i,"Start"])+int(row_cons["length"])
				df_cyto.at[i,"End"] = int(df_cyto.at[i,"End"])+int(row_cons["length"])
		
			
		df_cyto = df_cyto.sort_values(by=['Chromosome', 'Start', 'End'], ascending=[1, 1, 1])

		#print(df_cyto)

	return df_cyto, last_index_df_cyto

	
	
if __name__ =='__main__':
	main(sys.argv[1:])
