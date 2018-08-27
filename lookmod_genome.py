### A script allowing to look at the reference genome and to modify it

##########################################################################IMPORTS##########################################################################

import os
import sys
import getopt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pandas as pd

#########################################################################VISUALIZE#########################################################################

##########
###LIST###
##########
###If you want to list the chromosomes in your file, with their length
###python lookmod_genome.py -m list -f file.fasta

##########
###LOOK###
##########
###If you want to look at a spot with 10 bases before and after (between 12 and 22) where you want to insert a sequence you can write
###python lookmod_genome.py -m look -f file.fasta -p 'chr1:12:22' -l 10
###You will get the following area 12-9, so 3, to 22+9, so 31
###from3to12**from22to31

#########################################################################INTERACT#########################################################################

###For the next mode, you will need a fasta file with the modifications you want
###Then, launch this command
###python lookmod_genome.py -m modify -f file.fasta -c changes_file.fasta


############
###MODIFY###
############
###All informations will be contained in the header of fasta sequences
###Exemples:
###------DELETION------
###Delete sequence in chr4 from 40 to 50, will delete base from 41 to 49
###>deletion:chr4:40:50
###INFO : No sequence needed for deletion
###INFO : The result of >deletion:chr4:40:40 and >deletion:chr4:40:41 will delete nothing
###------INSERTION------
###Insert sequence in chr1 at position 25
###>insertion:chr1:25:26
###GCTAGCTAGC
###Insert sequence in chr2 from 5 to 15 (so deletion from 6 to 14 then insertion after base 5)
###>insertion:chr1:5:15
###GTTGCATGCTATC
###INFO : >insertion:chr1:15:15 will insert position 25
###------ADD------
###If you want to add a new chromosome, you will need to create a header like this, don't forget your sequence below
###>add:your_chromosome_name
###GTCGATCGTCATGGTT
###------REMOVE------
###If you want to remove a new chromosome, you will need to create a header like this
###>remove:your_chromosome_name
###INFO : No sequence needed for deletion

###INFO : You can combine DELETION, INSERTION, ADD and REMOVE in the same fasta file

###All positions will be calculated on reference genome, this way a modification will not impact the genome position for the next modification
###Use ":" to separate your entities



def usage():
	print('Usage:')
	print('\tpython '+sys.argv[0]+' -m <mode : list, look, modify> -f <genome fasta file> -c <changes file> -p <position> [-i <info file>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-m or --mode : list, look, modify')
	print('\t\t-f or --fasta : genome fasta file')
	print('\t\t-c or --change : changes file (if mode is modify)')
	print('\t\t-p or --position : position chrmX:start:end (if mode is look)')
	print('\t\t-l or --length : position (default : 10) (if mode is look)')
	print('\t\t-i or --info : output in file the terminal ouput')

def main(argv):

	mode = ""
	file_fasta = ""
	changes_file = ""
	position = ""
	full_seq=""
	length=10
	chrom_dict={}
	info_file = ""
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'm:f:c:p:l:i:', ['mode=', 'fasta=', 'change=', 'position=', 'length=', 'info=', 'help'])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

###################################################################################OPTIONS###################################################################################

	if not opts :
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit(2)
		elif opt in ('-m', '--mode'):
			mode = arg
		elif opt in ('-f', '--fasta'):
			file_fasta = arg
		elif opt in ('-c', '--change'):
			changes_file = arg
		elif opt in ('-p', '--position'):
			position = arg
		elif opt in ('-l', '--length'):
			length = int(arg)
		elif opt in ('-i', '--info'):
			info_file = arg
		else:
			print("Error : Bad option -> "+opt)
			usage()
			sys.exit(2)

###OPEN INFORMATION FILE IF SELECTED AS OPTION
	if info_file != "":
		if not os.path.exists(info_file):
			info_handle = open(info_file, 'a')
		else:
			print("Error : This info file already exist !\n")
			usage()
			sys.exit(2)
###CHOOSE A MODE	
	if mode!="list" and mode!="modify" and mode!="look":
		print("Error : You have to choose a mode : list, look or modify !\n")
		usage()
		sys.exit(2)

###CHECK THE REFERENCE FASTA FILE
	if (file_fasta[-3:]!=".fa" and file_fasta[-4:]!=".fna" and file_fasta[-6:]!=".fasta") or not os.path.exists(file_fasta):
		print("Error : The reference fasta file is missing, not .fa, .fna or .fasta !\n")
		usage()
		sys.exit(2)
###DICTIONNARY ALL CHROMOSOMES AND LENGTHS
	else:
		chrom_dict = {
			record.id: len(record.seq)
			for record in SeqIO.parse(file_fasta, 'fasta')
			if record.id not in chrom_dict
		}

###################################################################################CHECK###################################################################################

###CHECK FOR MODIFY MODE
	if mode=="modify":
		if (changes_file[-3:]!=".fa" and changes_file[-4:]!=".fna" and changes_file[-6:]!=".fasta") or not os.path.exists(changes_file):
			print("Error : The changes file is missing, not .fa, .fna or .fasta !\n")
			usage()
			sys.exit(2)
###CHECK FOR LOOK MODE
	elif mode=="look":
		if length < 1 or length > 2000:
			print("Error : length has to be between 1 and 100 !\n")
			usage()
			sys.exit(2)
		if position!="":
			if len(position.split(":")) == 3 :
				chrom=position.split(":")[0]
				try:
					start=int(position.split(":")[1])
				except ValueError:
					print("Error : start need to be integer, please use this pattern 'chrX:start:end' !\n")
					usage()
					sys.exit(2)
				try:
					end=int(position.split(":")[2])
				except ValueError:
					print("Error : end need to be integer, please use this pattern 'chrX:start:end' !\n")
					usage()
					sys.exit(2)
				
				if chrom not in chrom_dict:
					print("Error : "+chrom+" is not in the reference genome file, please use this pattern 'chrX:start:end' !\n")
					usage()
					sys.exit(2)
				if start < 1 or end < 1:
					print("Error : start or end position < 1 is not allowed for look mode, please use this pattern 'chrX:start:end' !\n")
					usage()
					sys.exit(2)
				if end < start:
					print("Error : end < start is not allowed, please use this pattern 'chrX:start:end' !\n")
					usage()
					sys.exit(2)
			else:
				print("Error : The position entry is not correct, please use this pattern 'chrX:start:end' !\n")
				usage()
				sys.exit(2)
		else:
			print("Error : The position entry is empty, please use this pattern 'chrX:start:end' !\n")
			usage()
			sys.exit(2)

	output_rename = file_fasta[:-3]+"_modified"+".fa"

###################################################################################PRINTS###################################################################################

	print('\n-----------------------------------------')
	print('Mode : '+mode)
	print('fasta file : '+file_fasta)
	if mode=="modify":
		print('change file : '+changes_file)
	if mode=="look":
		print('position : '+position)
		print('length : '+str(length))
	if mode=="modify":
		print('Output file rename: '+output_rename)
	if info_file != "":
		print('Output Information File : '+info_file)
	print('-----------------------------------------\n')

###################################################################################PROGRAM###################################################################################

### LIST MODE
	if mode=="list":
		print("----------------------")
		print("------   LIST   ------")
		print("----------------------")
		for key, value in chrom_dict.iteritems():
			print(key+"\t->\t"+str(value))
			if info_file != "":
				info_handle.write(key+"\t"+str(value)+"\n")

### LOOK MODE
	elif mode=="look":
		for record in SeqIO.parse(file_fasta, 'fasta'):
			if chrom == record.id:
				if start > len(record.seq):
					print("Error : start is bigger than the reference sequence length, please use this pattern 'chrX:start:end' !\n")
					usage()
					sys.exit(2)
				if end > len(record.seq):
					print("Warning : end is bigger than the reference sequence length, only start will be shown' !\n")
					if start<length:
						full_seq+=record.seq[0:start]+"**"
					else:
						full_seq+=record.seq[start-length:start]+"**"

				else:
					if start<length:
						full_seq+=record.seq[0:start]+"**"
					else:
						full_seq+=record.seq[start-length:start]+"**"

					if len(record.seq)-end<length:
						full_seq+=record.seq[end-1:len(record.seq)]
					else:
					
						full_seq+=record.seq[end-1:end-1+length]
		print("----------------------")
		print("------   LOOK   ------")
		print("----------------------")
		print("-> "+full_seq)
		if info_file != "":
			info_handle.write(str(full_seq))

### MODIFY MODE
	elif mode=="modify":
		removed_list=[]
		added_dict={}
		df_file=pd.DataFrame(columns=["chr", "modif_type", "length", "start", "end", "seq"])

		if info_file != "":
			orig_stdout = sys.stdout
			sys.stdout = info_handle


	###FILTER CHANGES_FILE, KEEP GOOD DELETION AND INSERTION. DISPATCH ADD AND REMOVE IN LIST
		for record in SeqIO.parse(changes_file, 'fasta'):
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

					df_file = df_file.append(pd.Series([chrom, modif_type, len(record.seq)-(end-start-1), start, end, str(record.seq)], index=["chr", "modif_type", "length", "start", "end", "seq"]), ignore_index=True)
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

						df_file = df_file.append(pd.Series([chrom, modif_type, -(end-start-1), start, end, ""], index=["chr", "modif_type", "length", "start", "end", "seq"]), ignore_index=True)
				else:
					print("Warning : {"+record.id+"} does not have the good format")
					print("Warning : To delete or insert use modif_type:chrmX:start:end !")
					print("Warning : To add or remove use modif_type:chrmX !")
					print("Warning : {"+record.id+"} will not be used")
					sys.exit()

	###SORT DATAFRAME
		df_file = df_file.sort_values(by=['chr', 'start', 'end'], ascending=[1, 1, 1])

		pd.options.mode.chained_assignment = None

		#print(df_file)

		records=[]
	###MODIFY EACH RECORD WITH THE DATAFRAME, DO NOT OUTPUT THE REMOVED CHROMOSOME
		for record in SeqIO.parse(file_fasta, 'fasta'):
			###if the record is not in the deleted list
			#print("Normal")
			#print(record.seq)
			if record.id not in removed_list:
				df_temp = df_file.loc[df_file['chr'] == record.id]
				for index, row in df_temp.iterrows():
					if row['modif_type']=='deletion':
						record.seq=record.seq[:row['start']]+record.seq[row['end']-1:]
					elif row['modif_type']=='insertion':
						record.seq=record.seq[:row['start']]+row['seq']+record.seq[row['end']-1:]
					else:
						print("Warning : {"+record.id+"} does not have the good format")
						print("Warning : To delete or insert use modif_type:chrmX:start:end !")
						print("Warning : To add or remove use modif_type:chrmX !")
						print("Warning : {"+record.id+"} will not be used")
						sys.exit()
					df_temp.start = df_temp.start + row['length']
					df_temp.end = df_temp.end + row['length']
				records.append(record)
			

		for key, value in added_dict.iteritems():
			records.append(SeqRecord(Seq(str(value),IUPAC.unambiguous_dna), id = key, description = key, name = key))
		SeqIO.write(records, output_rename, "fasta")

		sys.stdout = orig_stdout

		print("------------------------")
		print("------   MODIFY   ------")
		print("------------------------")
		print("-> Done")
		print(output_rename)
		if info_file != "":
			info_handle.write("-> Done\n")
			info_handle.write(output_rename)
	else:
		print("Error : You have to choose a mode : list, look or modify !")
		usage()
		sys.exit(2)


if __name__ =='__main__':
	main(sys.argv[1:])
