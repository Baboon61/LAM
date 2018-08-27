### A script allowing to run FastQC software on all the librairies of metadata.txt

##########################################################################IMPORTS##########################################################################

import os
import sys
import getopt
import pandas as pd

def usage():
	print('Usage:')
	print('\tpython '+sys.argv[0]+' -f <metadata.txt> -d <results directory> -o <output directory>')
	print('\t\t-h or --help : display this help')
	print('\t\t-f or --metadata : metadata txt file')
	print('\t\t-d or --directory : preprocess directory')
	print('\t\t-o or --output : output directory')

def main(argv):

	command_line = ""
	metadata_file = ""
	directory_path = ""
	output_path = "/home/lam/Documents/LAM/Data/temp/FastQC/"
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'f:d:o:', ['metadata=', 'directory=', 'output=', 'help'])
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
		elif opt in ('-f', '--metadata'):
			metadata_file = arg
		elif opt in ('-d', '--directory'):
			directory_path = arg
		elif opt in ('-o', '--output'):
			output_path = arg
		else:
			print("Error : Bad option")
			usage()
			sys.exit(2)
###CHECK METADATA FILE
	if metadata_file=="" or not os.path.exists(metadata_file):
		print("Error : You have to set a metadata file !\n")
		usage()
		sys.exit(2)

###CHECK RESULTS DIRECTORY
	if not os.path.exists(directory_path):
		print("Error : You have to set a preprocess directory !\n")
		usage()
		sys.exit(2)
	else:
		if directory_path[-1] != "/":
			directory_path+="/"

###CHECK OUTPUT DIRECTORY
	if not os.path.exists(output_path):
		os.system('mkdir '+output_path)
	if output_path[-1] != "/":
		output_path+="/"

###################################################################################PRINTS###################################################################################

	print('\n-----------------------------------------')
	print('Metadata File tlx : '+metadata_file)
	print('Preprocess Directory : '+directory_path)
	print('Output Directory : '+output_path)
	print('-----------------------------------------\n')

###################################################################################PROGRAM###################################################################################

###READ METADATA FILE
	metadata = pd.read_table(metadata_file, sep='\t')

###LOOP OVER EACH LIBRARIES
	for i in metadata['Library'].tolist():
		if not os.path.exists(directory_path+i+"_R1.fq.gz"):
			print("Warning : "+directory_path+" does not contains {"+i+"_R1.fq.gz}")
			print("Warning :  {"+i+"} will not be filtered")
		elif not os.path.exists(directory_path+i+"_R2.fq.gz"):
			print("Warning : "+directory_path+" does not contains {"+i+"_R2.fq.gz}")
			print("Warning :  {"+i+"} will not be filtered")
		else:
###CREATE THE COMMAND LINE
			command_line = "fastqc "+directory_path+i+"_R1.fq.gz "+directory_path+i+"_R2.fq.gz "
			print(command_line)
			os.system(command_line)
			command_line = "mv "+directory_path+i+"_R1_fastqc.* "+output_path
			print(command_line)
			os.system(command_line)
			command_line = "mv "+directory_path+i+"_R2_fastqc.* "+output_path
			print(command_line)
			os.system(command_line)

if __name__ =='__main__':
	main(sys.argv[1:])
