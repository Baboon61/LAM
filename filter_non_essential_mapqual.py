### A script allowing to filter out junctions flagged mapqual when there are outside the experimentaly duplicated parts

##########################################################################IMPORTS##########################################################################
import os
import sys
import getopt
import numpy as np
import pandas as pd


import sys
from datetime import datetime
import csv

###Experimentaly we duplicate part of the genome (chr12 to chr6). This lead to trigger the mapqual filter of the translocPipeline for this regions. We want to keep these junctions. So we disabled the mapqual filter but then we need to remove the "real" mapqual junctions.

def usage():
	print('Usage:\n')
	print('\tpython '+sys.argv[0]+' -f <metadata.txt> -d <results directory> -m <mark> -g <genome type> -c <duplicate file> [-n <marks input>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-f or --metadata : metadata file')
	print('\t\t-d or --directory : results directory')
	print('\t\t-m or --outputMark : mark added to the output file')
	print('\t\t-g or --genome : only filter librairies results with this genome')
	print('\t\t-n or --inputMark : marks from input file')
	print('\t\t-c or --construction : duplicate position from the experimental construction (chr,start,end)')

def main(argv):

	file_metadata = ""
	directory_path = ""
	mark=""
	genome=""
	input_marks = ""
	file_construction = ""
	file_input_extension = ""
	file_output_extension = ""

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'f:d:m:g:n:c:', ['metadata=', 'directory=', 'outputMark=', 'genome=', 'inputMark', 'construction=', 'help'])
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
		elif opt in ('-d', '--directory'):
			directory_path = arg
		elif opt in ('-m', '--outputMark'):
			mark = arg
		elif opt in ('-g', '--genome'):
			genome = arg
		elif opt in ('-n', '--inputMarK'):
			input_marks = arg
		elif opt in ('-c', '--construction'):
			file_construction = arg
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
	if not os.path.exists(directory_path):
		print("Error : You have to set a results directory !\n")
		usage()
		sys.exit(2)
	else:
		if directory_path[-1] != "/":
			directory_path+="/"

###FILTER METADATA FILE IF GENOME INPUT
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

###CHECK INPUT MARKS HISTORY
	if input_marks == "":
		print("Warning : You will filter the raw file !\n")

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
		if os.path.exists(directory_path+i+"/"+i+file_input_extension):
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

###CHECK CONSTRUCTION FILE
	if file_construction=="" or not os.path.exists(file_construction):
		print("Error : You have to set a construction file !\n")
		usage()
		sys.exit(2)
	else:
		###READ CONSTRUCTION FILE
		construction_locus = pd.read_table(file_construction, sep='\t', header=None)
		for index, row in construction_locus.iterrows():
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

###################################################################################PRINTS###################################################################################

	print('\n-----------------------------------------')
	print('Metadata File tlx : '+file_metadata)
	print('Results Directory : '+directory_path)
	print('Genome : '+genome)
	print('Construction File : '+file_construction)
	print('Input file extension: '+file_input_extension)
	print('Output file extension : '+file_output_extension)
	print('-----------------------------------------\n')

###################################################################################PROGRAMS###################################################################################

	###LOOP OVER EACH LIBRARIES
	for i in metadata['Library'].tolist():
		###CHECK DIRECTORY EXISTS
		if not os.path.exists(directory_path+i):
			print("Warning : "+directory_path+" does not contains {"+i+"}")
			print("Warning :  {"+i+"} will not be filtered")
		else:
			###CHECK INPUT FILE EXISTS
			if os.path.exists(directory_path+i+"/"+i+file_input_extension):
				file_output = open(directory_path+i+"/"+i+file_output_extension, 'a')
				dataframe_input = pd.read_table(directory_path+i+"/"+i+file_input_extension, sep='\t')
				file_output.write("\t".join(list(dataframe_input))+"\n")
				for index, row in dataframe_input.iterrows():
					if row.mapqual == 1:
						if construction_locus[(construction_locus[0] == row['Rname']) & (construction_locus[1] < row['Junction']) & (construction_locus[2] > row['Junction'])].empty or construction_locus[(construction_locus[4] == row['Rname']) & (construction_locus[5] < row['Junction']) & (construction_locus[6] > row['Junction'])].empty:
							print(row['Qname'])							
							pd.DataFrame(row).T.to_csv(file_output, sep='\t', header=None, index=False)
					else:
						pd.DataFrame(row).T.to_csv(file_output, sep='\t', header=None, index=False)
						
			else:
				print("Warning : "+directory_path+i+"/"+i+file_input_extension+" does not exist")
				print("Warning :  {"+directory_path+i+"/"+i+file_input_extension+"} will not be filtered")

if __name__ =='__main__':
	main(sys.argv[1:])
