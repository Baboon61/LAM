### A script allowing to run CSReport software on all the librairies of metadata.txt

##########################################################################IMPORTS##########################################################################

import os
import sys
import shutil
import getopt
import pandas as pd


def usage():
	print('Usage:')
	print('\tpython '+sys.argv[0]+' -d <CSReport directory> -f <metadata.txt> -k <postprocess directory> -t <annotation file> -a <from annotation> -b <to annotation> -g <genome type> -j <reference fasta file> [-u <unlink> -c <cluster size> -r -n <marks input>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-d or --CSReport : CSReport directory')
	print('\t\t-f or --metadata : metadata txt file')
	print('\t\t-k or --pos_directory : postprocess directory')
	print('\t\t-t or --file_annot : annotation bed file (locus, start, end)')
	print('\t\t-a or --from_annot : search junction from this locus ("all" to get all)')
	print('\t\t-b or --to_annot : search junction to this annotation ("all" to get all)')
	print('\t\t-u or --unlink : unlink the from_annot1/to_annot1, from_annot2/to_annot2 process and allows to search for from_annot1/to_annot2 junctions')
	print('\t\t-g or --genome : only filter librairies results with this genome')
	print('\t\t-j or --fasta_reference : reference fasta file')
	print('\t\t-c or --cluster_size : count number of clustered junctions (Default : 2)')
	print('\t\t-r or --reverse_CSR : do normal + reverse CSReport (change from_annot and to_annot)')
	print('\t\t-n or --inputMark : marks from input file')

def main(argv):

	path_CSReport = ""
	metadata_file = ""
	path_pos_directory = ""
	file_annot = ""
	from_annot = ""
	to_annot = ""
	unlink = False
	genome = ""
	file_reference = ""
	cluster_size = 2
	input_marks = ""
	reverse=False
	file_input_extension = ""

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'd:f:k:t:a:b:ug:j:c:rn:', ['CSReport=', 'metadata=', 'pos_directory=', 'file_annot=', 'from_annot=', 'to_annot=', 'unlink', 'genome=', 'fasta_reference=', 'cluster_size=', 'reverse=', 'inputMark=', 'help'])
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
		elif opt in ('-d', '--CSReport'):
			path_CSReport = arg
		elif opt in ('-f', '--metadata'):
			metadata_file = arg
		elif opt in ('-k', '--pos_directory'):
			path_pos_directory = arg
		elif opt in ('-t', '--file_annot'):
			file_annot = arg
		elif opt in ('-a', '--from_annot'):
			from_annot = arg
		elif opt in ('-b', '--to_annot'):
			to_annot = arg
		elif opt in ('-u', '--unlink'):
			unlink = True
		elif opt in ('-g', '--genome'):
			genome = arg
		elif opt in ('-j', '--fasta_reference'):
			file_reference = arg
		elif opt in ('-c', '--cluster_size'):
			cluster_size = arg
		elif opt in ('-r', '--reverse'):
			reverse = True
		elif opt in ('-n', '--inputMark'):
			input_marks = arg
		else:
			print("Error : Bad option")
			usage()
			sys.exit(2)

###CHECK CSREPORT DIRECTORY
	if not os.path.exists(path_CSReport):
		print("Error : You have to set a CSReport directory !\n")
		usage()
		sys.exit(2)
	else:
		if path_CSReport[-1] != "/":
			path_CSReport+="/"
	sys.path.append(path_CSReport)

###CHECK METADATA FILE
	if metadata_file=="" or not os.path.exists(metadata_file):
		print("Error : You have to set a metadata file !\n")
		usage()
		sys.exit(2)
	else:
		###READ METADATA FILE
		metadata = pd.read_table(metadata_file, sep='\t')

###CHECK POSTPROCESS DIRECTORY
	if not os.path.exists(path_pos_directory):
		print("Error : You have to set a postprocess directory !\n")
		usage()
		sys.exit(2)
	else:
		if path_pos_directory[-1] != "/":
			path_pos_directory+="/"

###CHECK ANNOTATION FILE
	if file_annot=="" or not os.path.exists(file_annot):
		print("Error : You have to set an annotation file !\n")
		usage()
		sys.exit(2)
	else:
		###READ ANNOTATION FILE
		df_annot = pd.read_table(file_annot, sep='\t', header=None)
		for index, row in df_annot.iterrows():
			try :
				if int(row[1]) < 1:
					print("Error : line."+str(index+1)+", col.2 of your annotation file !\n")
					print("Error : Start position has to be positive integer !\n")
					usage()
					sys.exit(2)
			except :
				print("Error : line."+str(index+1)+", col.2 of your annotation file !\n")
				print("Error : Unknown start position : "+str(row[1])+" !\n")
				usage()
				sys.exit(2)
			try : 
				if int(row[2]) < 1:
					print("Error : line."+str(index+1)+", col.3 of your annotation file !\n")
					print("Error : End position has to be positive integer !\n")
					usage()
					sys.exit(2)
			except :
				print("Error : line."+str(index+1)+", col.3 of your annotation file !\n")
				print("Error : Unknown end position : "+str(row[2])+" !\n")
				usage()
				sys.exit(2)
			if int(row[2]) < int(row[1]):
				print("Error : line."+str(index+1)+" of your annotation file !\n")
				print("Error : End position is smaller than start position !\n")
				usage()
				sys.exit(2)

###CHECK FROM ANNOTATION OPTION
	if from_annot != "":
		if 'all' in from_annot:
			list_from_annot = df_annot[0].tolist()
		else:
			list_from_annot = from_annot.split(",")
			for i in list_from_annot:
				if i not in df_annot[0].tolist():
					print("Error : {"+i+"} does not exist in the annotation file !\n")
					usage()
					sys.exit(2)

###CHECK TO ANNOTATION OPTION
	if to_annot != "":
		if 'all' in to_annot:
			list_to_annot = df_annot[0].tolist()
		else:
			list_to_annot = to_annot.split(",")
			for i in list_to_annot:
				if i not in df_annot[0].tolist():
					print("Error : {"+i+"} does not exist in the annotation file !\n")
					usage()
					sys.exit(2)

###CHECK UNLINK
	if not unlink:
		if 'all' in from_annot or 'all' in to_annot:
			print("Warning : Default unlink if -a all or -b all !")
			unlink=True			
		else:
			if len(list_from_annot) != len(list_to_annot):
				print("Error : The number of element in from_annotation and to_annotation have to be the same if you link them !")
				usage()
				sys.exit(2)
			unlink=False

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

###CHECK REFERENCE FILE
	if (file_reference[-3:]!=".fa" and file_reference[-4:]!=".fna" and file_reference[-6:]!=".fasta") or not os.path.exists(file_reference):
		print("Error : The reference file is not fasta or does not exist !\n")
		usage()
		sys.exit(2)

###CHECK CLUSTER SIZE
	try:
		cluster_size = int(cluster_size)
		if cluster_size < 1:
			print("Error : Cluster size option needs to be more than 0 !\n")
			usage()
			sys.exit(2)
	except:
		print("Error : You have to set an integer to cluster size option !\n")
		usage()
		sys.exit(2)

###CHECK INPUT MARKS HISTORY
	if input_marks == "":
		print("Warning : You will process the raw file !\n")

###SELECT INPUT FILES
	if input_marks == "":
		file_input_extension = ".tlx"
	else:
		file_input_extension = "_"+"_".join(input_marks.split(","))+".fasta"

###################################################################################PRINTS###################################################################################

	print('\n-----------------------------------------')
	print('Metadata File tlx : '+metadata_file)
	print('Postprocess Directory : '+path_pos_directory)
	print('Annotation File : '+file_annot)
	print('From Annotation : '+','.join(list_from_annot))
	print('To Annotation : '+','.join(list_to_annot))
	if unlink:
		print('Unlink : Yes')
	else:
		print('Unlink : No')
	print('Genome : '+genome)
	print('Reference fasta file : '+file_reference)
	print('Cluster size : '+str(cluster_size))
	if reverse:
		print('Reverse CSReport : Yes')
	else:
		print('Reverse CSReport : No')
	print('Input file extension: '+file_input_extension)
	print('-----------------------------------------\n')

###################################################################################PROGRAM###################################################################################

###IMPORT CSReport
	from CSReport import CSReport_run,CSReport_summary,histogramStructures,sketchBP,evaluateDiversity,summarizeMotifs

###LOOP OVER EACH LIBRARIES
	for library in metadata['Library'].tolist():
		if not os.path.exists(path_pos_directory+library):
			print("Warning : "+path_pos_directory+" does not contains {"+library+"}")
			print("Warning :  {"+library+"} will not be filtered")
		else:
###CREATE THE COMMAND LINE
			os.mkdir(path_pos_directory+library+"/CSReport"+file_input_extension[0:-6])
			for j in list_from_annot:
				for k in list_to_annot:
					print(library,j,k)
					path_to_fasta = path_pos_directory+library+"/"
					fasta = library+"_"+k.replace("/","-").replace(" ","+").replace("'","!")+file_input_extension
					if os.stat(path_to_fasta+fasta).st_size != 0:
						if not os.path.exists(path_to_fasta+fasta):
							print("Warning : "+path_to_fasta+fasta+" does not exist")
							print("Warning :  {"+path_to_fasta+fasta+"} will not be reported")
						else:
							print("NORMAL")
							CSReport_run(Seq=path_to_fasta+fasta[0:-6],Ref=file_reference[0:-6],Region1=j,Region2=k,Form='fasta')
							if os.path.exists(path_to_fasta+fasta[0:-6]+"/"+fasta[0:-6]+"_readsJunction.csv"):
								df_junction = pd.read_table(path_to_fasta+fasta[0:-6]+"/"+fasta[0:-6]+"_readsJunction.csv", sep="\t", header=0)
								if len(df_junction.index) > 0:
									current_path = os.getcwd()
									os.chdir(path_to_fasta)

									#CSReport_summary(fasta[0:-6],cluster=cluster_size)
									#Structures=histogramStructures(Seq=fasta[0:-6],bw=False)
									#sketchBP(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=j,Region2=k,density2D=False)

									#evaluateDiversity(Seq=fasta[0:-6])
									#summarizeMotifs(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=j,Region2=k)

									os.chdir(current_path)
									shutil.move(path_to_fasta+fasta[0:-6], path_to_fasta+"CSReport"+file_input_extension[0:-6]+"/"+j.replace("/","-").replace(" ","+").replace("'","!")+"_"+k.replace("/","-").replace(" ","+").replace("'","!")+"_"+fasta[0:-6])
								else:
									shutil.rmtree(path_to_fasta+fasta[0:-6])
					

							
							if reverse and j!=k:
								print("REVERSE")
								CSReport_run(Seq=path_to_fasta+fasta[0:-6],Ref=file_reference[0:-6],Region1=k,Region2=j,Form='fasta')
								if os.path.exists(path_to_fasta+fasta[0:-6]+"/"+fasta[0:-6]+"_readsJunction.csv"):
									df_junction = pd.read_table(path_to_fasta+fasta[0:-6]+"/"+fasta[0:-6]+"_readsJunction.csv", sep="\t", header=0)
									if len(df_junction.index) > 0:
										current_path = os.getcwd()
										os.chdir(path_to_fasta)

										#CSReport_summary(fasta[0:-6],cluster=cluster_size)
										#Structures=histogramStructures(Seq=fasta[0:-6],bw=False)
										#sketchBP(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=k,Region2=j,density2D=False)

										#evaluateDiversity(Seq=fasta[0:-6])
										#summarizeMotifs(Seq=fasta[0:-6],Ref=file_reference[0:-6],Region1=k,Region2=j)

										os.chdir(current_path)
										shutil.move(path_to_fasta+fasta[0:-6], path_to_fasta+"CSReport"+file_input_extension[0:-6]+"/"+k.replace("/","-").replace(" ","+").replace("'","!")+"_"+j.replace("/","-").replace(" ","+").replace("'","!")+"_"+fasta[0:-6])
									else:
										shutil.rmtree(path_to_fasta+fasta[0:-6])

			###MOVE ALL FASTA FILES
			print("mv "+path_to_fasta+"*"+file_input_extension+" "+path_to_fasta+"CSReport"+file_input_extension[0:-6]+"/")
			os.system("mv "+path_to_fasta+"*"+file_input_extension+" "+path_to_fasta+"CSReport"+file_input_extension[0:-6]+"/")							

if __name__ =='__main__':
	main(sys.argv[1:])
