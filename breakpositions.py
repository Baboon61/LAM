###A script that manage to create an horizontal bar visualisation of double strand break in LAM-HTGTS method

##########################################################################IMPORTS##########################################################################

import os
import sys
import getopt
import os.path
import re
import csv
import pandas as pd
import numpy as np
import math
import operator
import collections
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.ticker import MultipleLocator

###USED
#(0.0, 0.0, 0.0, 1.0) = black = intra
#(1.0, 0.0, 0.0, 1.0)) = red = inter

def usage():
	print('Usage:\n')
	print('\tpython '+sys.argv[0]+' -f <metadata.txt> -k <pos_directory> -m <outputMark> -g <genome type> [-m [0,1] -n <marks input>]')
	print('\t\t-h or --help : display this help')
	print('\t\t-f or --metadata : metadata txt file')
	print('\t\t-k or --pos_directory : postprocess directory')
	print('\t\t-m or --outputMark : mark added to the output file')
	print('\t\t-g or --genome : only filter librairies results with this genome')
	print('\t\t-d [0,1] or --method [0,1] : 0 = output picture using counts, 1 = output picture using log2(counts)')
	print('\t\t-n or --inputMark : marks from input file')

###################################################################################FUNCTIONS###################################################################################

###Allow to divide all value with max_distance_count_all then remove max_distance_count and max_distance from dictionnary
def manageArray(dictionnary, max_distance_count_all, method):
	#print(dictionnary)
	del dictionnary['max_distance']
	if dictionnary['max_distance_count'] !=0:
		for key, value in dictionnary.items():
			if method == 0 :
				###Counts
				dictionnary[key]=float(float(dictionnary[key])/float(max_distance_count_all))
			elif method == 1:
				###Log2
				dictionnary[key]=float(float(np.log2(dictionnary[key]+1))/float(np.log2(max_distance_count_all+1)))
			else :
				print("Wrong method use 0 or 1")
				sys.exit(2)
	del dictionnary['max_distance_count']
	return dictionnary

###############################################################################################################################################################################

def main(argv):

	file_metadata = ""
	path_pos_directory = ""
	mark = ""
	genome = ""
	method=-1
	file_input_extension = ""
	file_output_extension = ""
	title_label = ""

	good_directory_list=[]
	check_result_tlx=False
	check_result_filter_stats_txt=False
	total_junction=0

	#pd.options.mode.chained_assignment = None
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'f:k:m:g:d:n:', ['metadata=', 'pos_directory=', 'outputMark=', 'genome=', 'method=', 'inputMark=', 'help'])
	except getopt.GetoptError:
		usage()
		sys.exit(2)

###################################################################################OPTIONS###################################################################################

	print("\n")
	for opt, arg in opts:
		if opt in ('-h', '--help'):
			usage()
			sys.exit(2)
		elif opt in ('-f', '--metadata'):
			file_metadata = arg
		elif opt in ('-k', '--pos_directory'):
			path_pos_directory = arg
		elif opt in ('-m', '--outputMark'):
			mark = arg
		elif opt in ('-g', '--genome'):
			genome = arg
		elif opt in ('-d', '--method'):
			method = arg
		elif opt in ('-n', '--inputMark'):
			input_marks = arg
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

	###CHECK POSTPROCESS DIRECTORY
	if not os.path.exists(path_pos_directory):
		print("Error : You have to set a postprocess directory !\n")
		usage()
		sys.exit(2)
	else:
		if path_pos_directory[-1] != "/":
			path_pos_directory+="/"

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

	###CHECK METHOD
	try:
		method = int(method)	
		if method == 0 :
			title_label="Number of double stranded breaks at position divided by the maximum of breaks of all position"
		elif method == 1:
			title_label="Log2 of number of double stranded breaks at position divided by log2 of the maximum of breaks of all position"
		else:
			print("Error : Method option needs to be 0 or 1 !\n")
			usage()
			sys.exit(2)
	except:
		print("Error : You have to set an integer (0 or 1) to method option !\n")
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
	for library in metadata['Library'].tolist():	
		if os.path.exists(path_pos_directory+library+"/"+library+"_Legitimate"+file_input_extension) or os.path.exists(path_pos_directory+library+"/"+library+"_Illegitimate"+file_input_extension):
			check_inputMark = True
	if not check_inputMark:
		print("Error : Your input marks can not localize a good legitimate or illegitimate file !\n")
		usage()
		sys.exit(2)

	###SELECT OUTPUT FILES
	if file_input_extension != "":
		file_output_extension = file_input_extension[:-4]+"_"+mark+".png"
	else:
		file_output_extension = "_"+mark+".png"

###################################################################################PRINTS###################################################################################

	print('\n-----------------------------------------')
	print('Metadata File tlx : '+file_metadata)
	print('Postprocess Directory : '+path_pos_directory)
	print('Genome : '+genome)
	print('Method : '+str(method))
	print('Input file extension: '+file_input_extension)
	print('Output file extension : '+file_output_extension)
	print('-----------------------------------------\n')

###################################################################################PROGRAMS###################################################################################

	fig = plt.figure()
	ax = fig.add_subplot(111)
	d = {}
	yticks = []
	yticklabels = []
	track_number=0
	yrange = []
	coff_leg_ille_graph=0

	###LOOP OVER EACH LIBRARIES
	for library in metadata['Library'].tolist():
		#print(library)
		###CHECK DIRECTORY EXISTS
		if not os.path.exists(path_pos_directory+library):
			print("Warning : "+path_pos_directory+" does not contains {"+library+"}")
			print("Warning :  {"+library+"} will not be filtered")
		else:
			###CHECK INPUT FILE EXISTS
			if os.path.exists(path_pos_directory+library+"/"+library+"_Legitimate"+file_input_extension):
				df_legitimate = pd.read_csv(path_pos_directory+library+"/"+library+"_Legitimate"+file_input_extension, sep='\t', header=0, index_col=None)
				df_legitimate = df_legitimate.drop(columns=df_legitimate.columns[11:])

			else:
				print("Error : The Legitimate file for "+library+" is missing !\n")
				usage()
				sys.exit(2)

			if os.path.exists(path_pos_directory+library+"/"+library+"_Illegitimate"+file_input_extension):
				df_illegitimate = pd.read_csv(path_pos_directory+library+"/"+library+"_Illegitimate"+file_input_extension, sep='\t', header=0, index_col=None)
				df_illegitimate = df_illegitimate.drop(columns=df_illegitimate.columns[11:])

			else:
				print("Error : The Illegitimate file for "+library+" is missing !\n")
				usage()
				sys.exit(2)

	metadata = metadata.sort_values(['Library'], ascending=[True])
	global_list=path_pos_directory + metadata['Library']
	#for i in global_list:
	#	print(i)
	
	nb_files=len(global_list)
	###Dictionnary with illegitimates and legitimates junctions
	distance_dict={"illegitimates":{}, "legitimates":{}}

	###Open one by one tlx files
	for result_file in global_list:
		label=result_file.split("/")[-1]

		###yrange used at the end to space horizontal bar

		yrange.append(((track_number*1.8)+(track_number*1.8)+1.8+coff_leg_ille_graph,1))
		yrange.append(((track_number*1.8)+(track_number*1.8)+3.6+coff_leg_ille_graph,1))

		coff_leg_ille_graph+=1

		###Dictionnary label inside illegitimates and legitimates dictionnaries
		distance_dict["illegitimates"][label]={}
		distance_dict["legitimates"][label]={}

		#with open(result_file+"/"+label+"_result.tlx", 'r') as f_tlx:
		with open(result_file+"/"+label+"_Legitimate"+file_input_extension, 'r') as f_tlx:
			f_tlx.readline()
			for line in f_tlx:
				distance=int(int(line.split("\t")[9])-int(line.split("\t")[8]))
				if distance not in distance_dict["legitimates"][label] :
					distance_dict["legitimates"][label][distance]=1
				else:
					distance_dict["legitimates"][label][distance]+=1
				total_junction+=1

		with open(result_file+"/"+label+"_Illegitimate"+file_input_extension, 'r') as f_tlx:
			f_tlx.readline()
			for line in f_tlx:
				distance=int(int(line.split("\t")[9])-int(line.split("\t")[8]))
				if distance not in distance_dict["illegitimates"][label] :
					distance_dict["illegitimates"][label][distance]=1
				else:
					distance_dict["illegitimates"][label][distance]+=1
				total_junction+=1

		### ADD max_distance and max_distance_count to both
		max_distance=0
		max_distance_count=0
		if len(distance_dict["illegitimates"][label]) > 0:
			#max_distance
			max_distance=max(distance_dict["illegitimates"][label].iteritems(), key=operator.itemgetter(0))[0]
			#max_distance_count
			max_distance_count=max(distance_dict["illegitimates"][label].iteritems(), key=operator.itemgetter(1))[1]
		distance_dict["illegitimates"][label]["max_distance"]=max_distance
		distance_dict["illegitimates"][label]["max_distance_count"]=max_distance_count

		max_distance=0
		max_distance_count=0
		if len(distance_dict["legitimates"][label]) > 0:
			#max_distance
			max_distance=max(distance_dict["legitimates"][label].iteritems(), key=operator.itemgetter(0))[0]
			#max_distance_count
			max_distance_count=max(distance_dict["legitimates"][label].iteritems(), key=operator.itemgetter(1))[1]
		distance_dict["legitimates"][label]["max_distance"]=max_distance
		distance_dict["legitimates"][label]["max_distance_count"]=max_distance_count
		track_number+=1

	#print(distance_dict)

	max_distance_all=0
	max_distance_count_all=0
	###Define max_distance_all and max_distance_count_all
	for legitimate, label_dict in distance_dict.items():
		for key,value in label_dict.items():
			if max_distance_all < int(label_dict[key]['max_distance']):
				max_distance_all = int(label_dict[key]['max_distance'])
			if max_distance_count_all < int(label_dict[key]['max_distance_count']):
				max_distance_count_all = int(label_dict[key]['max_distance_count'])

	#print(max_distance_all)
	#print(max_distance_count_all)

	yrange_count=0

	for bad_label in global_list:
		label=bad_label.split("/")[-1]

		distance_dict["illegitimates"][label] = manageArray(distance_dict["illegitimates"][label], max_distance_count_all, method)
		distance_dict["legitimates"][label] = manageArray(distance_dict["legitimates"][label], max_distance_count_all, method)

		"""distance_dict["illegitimates"][label] = manageArray(distance_dict["illegitimates"][label], total_junction, method)
		distance_dict["legitimates"][label] = manageArray(distance_dict["legitimates"][label], total_junction, method)"""

		distance_dict["illegitimates"][label] = collections.OrderedDict(sorted(distance_dict["illegitimates"][label].items()))
		distance_dict["legitimates"][label] = collections.OrderedDict(sorted(distance_dict["legitimates"][label].items()))
		
		###Concatenate legitimate and illegitimate distance and count(use for colors)
		xranges_leg=[]
		colors_leg=[]
		xranges_illeg=[]
		colors_illeg=[]

		#print(distance_dict["illegitimates"][label])
		for key,value in distance_dict["legitimates"][label].items():
			xranges_leg.append((key,1))
			colors_leg.append((0.0, 0.0, 0.0, value))
		
		for key,value in distance_dict["illegitimates"][label].items():
			xranges_illeg.append((key,1))
			colors_illeg.append((1.0, 0.0, 0.0, value))
		

		#############################################WORKS#############################################
		#xranges=[(29, 1),(30, 1),(31, 1)]
		#yrange=(1.8, 1)
		#colors=[(0.0, 0.0, 0.0, 0.012),(0.0, 0.0, 0.0, 0.008),(0.0, 0.0, 0.0, 0.012)]
		###############################################################################################

		#print("xranges")
		#print(xranges)
		#print(type(xranges))
		#print("yrange")
		#print(yrange[yrange_count])
		#print(type(yrange[yrange_count]))
		#print("colors")
		#print(colors)
		#print(type(colors))
		
		coll = BrokenBarHCollection(xranges_leg, yrange[yrange_count], facecolors=colors_leg, edgecolors=colors_leg)
		ax.add_collection(coll)
		center = yrange[yrange_count][0] + yrange[yrange_count][1]/2.0
		yticks.append(center)
		yticklabels.append(label+"_legi")
 		d[label+"legi"] = xranges_leg

		yrange_count+=1

		coll = BrokenBarHCollection(xranges_illeg, yrange[yrange_count], facecolors=colors_illeg, edgecolors=colors_illeg)
		ax.add_collection(coll)
		center = yrange[yrange_count][0] + yrange[yrange_count][1]/2.0
		yticks.append(center)
		yticklabels.append(label+"_illegi")
 		d[label+"illegi"] = xranges_illeg

		yrange_count+=1

	ax.axis('tight')
	#ax.set_xlim([0, max_distance_all+10])
	ax.set_yticks(yticks)
	ax.set_yticklabels(yticklabels)

	ax.set_xticks(range(0,max_distance_all+10,5))
	ax.get_xaxis().get_major_formatter().set_scientific(False)
	plt.xlabel("Position from bait start")
	plt.axes().xaxis.set_minor_locator(MultipleLocator(1))
	plt.title(title_label, fontdict=None, loc='center')
	fig = plt.gcf()
	fig.set_size_inches(22, 10)

	fig.savefig(path_pos_directory+genome+file_output_extension, format='png', bbox_inches='tight')
	#plt.show()


if __name__ =='__main__':
    main(sys.argv[1:])
