### A script allowing to filter out chromosomes and contigs outside the "conventional" chromosomes (chr1-19, X, Y for mouse)

#########################################################################IMPORTS#########################################################################

import os
import sys
import getopt
import pandas as pd
from Bio import SeqIO

###To keep the most junctions as possible only "conventional chromosomes" (readable with Circos Plot) will be kept. According to the GRC documentation (https://www.ncbi.nlm.nih.gov/grc/help/definitions), chromosomes and contigs outside the conventional assemblies will be deleted from the reference genome.

#Reads will not have the possibility to be align on fixed or patched chromosomes. This action will "force" reads to be align on the conventional chromosomes (to be correct : will not have the possibility to be align on alternative chromosomes). So in case a read is aligned on a secondary assembly but also have a good match on conventional chromosome, this solution will save a read that we would have lost downstream in the analysis (CircosPlot, Karyoplot, locus detection... all focus on conventional chromosomes)

#This is a manipulation to beware of. This works in that particular case because we want to know were is the read on our genome, not the intrinsic read composition. Typicaly, in SNPs discovery the step is very dangerous and could lead to false positive discovery.

def usage():
	print('Usage:')
	print('\tpython '+sys.argv[0]+' -r <reference fasta file> -a <autosome number> [-s -i <info file>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-r or --file_reference : reference fasta file')
	print('\t\t-a or --autosomes : number of automome in the reference genome')
	print('\t\t-s or --sexual : remove sexual chromosomes')
	print('\t\t-i or --info : output in file the terminal ouput')

def main(argv):

	file_reference = ""
	number_autosomes = 0
	sexual = False
	conventional_chrom_list = []
	records = []
	output_file = ""
	info_file = ""
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'r:a:si:', ['file_reference=', 'autosomes=', 'sexual', 'info=', 'help'])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

#########################################################################OPTIONS#########################################################################
	
	if not opts :
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit(2)
		elif opt in ('-r', '--file_reference'):
			file_reference = arg
		elif opt in ('-a', '--autosomes'):
			number_autosomes = arg
		elif opt in ('-s', '--sexual'):
			sexual = True
		elif opt in ('-i', '--info'):
			info_file = arg
		else:
			print("Error : Bad option")
			usage()
			sys.exit(2)

#########################################################################CHECK UP/SET UP#########################################################################

	#OPEN INFORMATION FILE IF SELECTED AS OPTION
	if info_file != "":
		if not os.path.exists(info_file):
			info_handle = open(info_file, 'a')
		else:
			print("Error : This info file already exist !\n")
			usage()
			sys.exit(2)

	#CHECK REFERENCE FASTA FILE
	if file_reference=="" or not os.path.exists(file_reference):
		print("Error : You have to set a reference fasta file !\n")
		usage()
		sys.exit(2)

	#CHECK AUTOSOMES NUMBER
	try: 
		number_autosomes = int(number_autosomes)
	except ValueError:
		print("Error : You have set a non integer autosome number !\n")
		usage()
		sys.exit(2)

	if number_autosomes < 1:
		print("Error : You have set an autosome number > 1 !\n")
		usage()
		sys.exit(2)

	output_file = '.'.join(file_reference.split(".")[:-1])+"_filtered."+file_reference.split(".")[-1]

#########################################################################PRINTS#########################################################################

	print('\n-----------------------------------------')
	print('Reference fasta file : '+file_reference)
	print('Autosomes : '+ str(number_autosomes))
	if sexual:
		print('Remove sexual chromosomes !')
	else:
		print('Keep sexual chromosomes')
	print('Output modified file : '+output_file)
	if info_file != "":
		print('Output information file : '+info_file)
	else:
		print('No output information file')
	print('-----------------------------------------\n')

#########################################################################PROGRAM#########################################################################

	#CREATE LIST OF CONVENTIONAL CHROMOSOMES
	conventional_chrom_list = [
		"chr"+str(i+1)
		for i in range(number_autosomes)
	]
	if not sexual:
		conventional_chrom_list.append("chrX")
		conventional_chrom_list.append("chrY")

	if info_file != "":
		for i in conventional_chrom_list:
			info_handle.write(i+"\n")

	#FILTER OUT CHROMOSOMES NOT IN CONVENTIONAL CHROMOSOME LIST
	records = [
		record
		for record in SeqIO.parse(file_reference, 'fasta')
		if record.id in conventional_chrom_list
	]
	SeqIO.write(records, output_file, "fasta")

if __name__ =='__main__':
	main(sys.argv[1:])
