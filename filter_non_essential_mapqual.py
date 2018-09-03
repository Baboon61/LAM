# A script allowing to filter out junctions flagged mapqual when there are
# outside the experimental duplicated parts

##############################IMPORTS##############################
import os
import sys
import getopt
import numpy as np
import pandas as pd
import csv

# Experimentaly, we duplicate part of the genome (chr12 to chr6). This
# lead to trigger the mapqual filter of the translocPipeline for this
# regions. We want to keep these junctions. So we disabled the mapqual
# filter but then we need to remove the "real" mapqual junctions.


def usage():
    print('Usage:\n')
    print('\tpython ' +
          sys.argv[0] + ' -m <metadata file> -g <genome type> -o <output mark> -c <construction fasta file> -t <results directory> [-i <input mark>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-g or --genome : only filter librairies results with this genome')
    print('\t\t-o or --output_mark : mark added to the output file name')
    print('\t\t-c or --file_construction : duplicate position from the experimental construction (chr,start,end)')
    print('\t\t-t or --dir_results : results directory')
    print('\t\t-i or --input_mark : marks from input file')


def main(argv):

    file_metadata = ""
    genome = ""
    output_mark = ""
    dir_results = ""
    input_mark = ""
    file_construction = ""
    file_input_extension = ""
    file_output_extension = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:g:o:c:t:i:', [
                                   'file_metadata=', 'genome=', 'output_mark=', 'file_construction=', 'dir_results=', 'input_mark', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

##############################OPTIONS##############################

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-m', '--file_metadata'):
            file_metadata = arg
        elif opt in ('-g', '--genome'):
            genome = arg
        elif opt in ('-o', '--output_mark'):
            output_mark = arg
        elif opt in ('-c', '--file_construction'):
            file_construction = arg
        elif opt in ('-t', '--dir_results'):
            dir_results = arg
        elif opt in ('-i', '--input_mark'):
            input_mark = arg

        else:
            print("Error : Bad option -> " + opt)
            usage()
            sys.exit(2)

##############################CHECK UP/SET UP##############################

    # CHECK METADATA FILE
    if file_metadata == "" or not os.path.exists(file_metadata):
        print("Error : You have to set a metadata file !\n")
        usage()
        sys.exit(2)
    else:
        # READ METADATA FILE
        metadata = pd.read_table(file_metadata, sep='\t')

     # FILTER METADATA FILE IF GENOME INPUT
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

    # CHECK RESULTS DIRECTORY
    if not os.path.exists(dir_results):
        print("Error : You have to set a results directory !\n")
        usage()
        sys.exit(2)
    else:
        if dir_results[-1] != "/":
            dir_results += "/"

    # CHECK INPUT MARKS HISTORY
    if input_mark == "":
        print("Warning : You will filter the raw file !\n")

    # CHECK OUTPUT MARK
    if output_mark == "":
        print("Error : You have to set an output mark !\n")
        usage()
        sys.exit(2)

    # SELECT INPUT FILES

    if input_mark == "":
        file_input_extension = ".tlx"
    else:
        file_input_extension = "_" + "_".join(input_mark.split(",")) + ".tlx"

    # TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
    check_input_mark = False
    for i in metadata['Library'].tolist():
        if os.path.exists(dir_results + i + "/" + i + file_input_extension):
            check_input_mark = True
    if not check_input_mark:
        print("Error : Your input marks can not localize a good input file !\n")
        usage()
        sys.exit(2)

    # SELECT OUTPUT FILES
    if file_input_extension != "":
        file_output_extension = file_input_extension[
            :-4] + "_" + output_mark + ".tlx"
    else:
        file_output_extension = "_" + output_mark + ".tlx"

    # CHECK CONSTRUCTION FILE
    if file_construction == "" or not os.path.exists(file_construction):
        print("Error : You have to set a construction file !\n")
        usage()
        sys.exit(2)
    else:
        # READ CONSTRUCTION FILE
        construction_locus = pd.read_table(
            file_construction, sep='\t', header=None)
        for index, row in construction_locus.iterrows():
            # DONOR CHECK
            if row[0][0:3] != 'chr':
                print("Error : line." + str(index + 1) +
                      ", col.1 of your construction file !\n")
                print("Error : Unknown chromosome : " + row[0] + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[1]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.2 of your construction file !\n")
                    print("Error : Start position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.2 of your construction file !\n")
                print("Error : Unknown start position : " +
                      str(row[1]) + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[2]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.3 of your construction file !\n")
                    print("Error : End position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.3 of your construction file !\n")
                print("Error : Unknown end position : " + str(row[2]) + " !\n")
                usage()
                sys.exit(2)
            if int(row[2]) < int(row[1]):
                print("Error : line." + str(index + 1) +
                      " of your construction file !\n")
                print("Error : End position is smaller than start position !\n")
                usage()
                sys.exit(2)
            if row[3] != '+' and row[3] != '-':
                print("Error : line." + str(index + 1) +
                      ", col.4 of your construction file !\n")
                print("Error : Unknown strand (+,-) : " + str(row[3]) + " !\n")
                usage()
                sys.exit(2)
            # ACCEPTOR CHECK
            if row[4][0:3] != 'chr':
                print("Error : line." + str(index + 1) +
                      ", col.5 of your construction file !\n")
                print("Error : Unknown chromosome : " + row[4] + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[5]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.6 of your construction file !\n")
                    print("Error : Start position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.6 of your construction file !\n")
                print("Error : Unknown start position : " +
                      str(row[5]) + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[6]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.7 of your construction file !\n")
                    print("Error : End position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.7 of your construction file !\n")
                print("Error : Unknown end position : " + str(row[6]) + " !\n")
                usage()
                sys.exit(2)
            if int(row[6]) < int(row[5]):
                print("Error : line." + str(index + 1) +
                      " of your construction file !\n")
                print("Error : End position is smaller than start position !\n")
                usage()
                sys.exit(2)
            if row[7] != '+' and row[7] != '-':
                print("Error : line." + str(index + 1) +
                      ", col.8 of your construction file !\n")
                print("Error : Unknown strand (+,-) : " + str(row[7]) + " !\n")
                usage()
                sys.exit(2)

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Genome : ' + genome)
    print('Construction file : ' + file_construction)
    print('Results directory : ' + dir_results)
    print('Input file extension: ' + file_input_extension)
    print('Output file extension : ' + file_output_extension)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # LOOP OVER EACH LIBRARIES
    for library in metadata['Library'].tolist():
        print(library)
        # CHECK RESULTS DIRECTORY EXISTS
        if not os.path.exists(dir_results + library):
            print("Warning : " + dir_results +
                  " does not contains {" + library + "}")
            print("Warning :  {" + library + "} will not be filtered")
        else:
            # CHECK INPUT FILE EXISTS
            if os.path.exists(dir_results + library + "/" + library + file_input_extension):
                file_output = open(dir_results + library + "/" +
                                   library + file_output_extension, 'a')
                dataframe_input = pd.read_table(
                    dir_results + library + "/" + library + file_input_extension, sep='\t')
                file_output.write("\t".join(list(dataframe_input)) + "\n")
                for index, row in dataframe_input.iterrows():
                    # IF MAPQUAL
                    if row.mapqual == 1:
                        # CHECK IF IN CONSTRUCTION
                        if construction_locus[(construction_locus[0] == row['Rname']) & (construction_locus[1] < row['Junction']) & (construction_locus[2] > row['Junction'])].empty or construction_locus[(construction_locus[4] == row['Rname']) & (construction_locus[5] < row['Junction']) & (construction_locus[6] > row['Junction'])].empty:
                            pd.DataFrame(row).T.to_csv(
                                file_output, sep='\t', header=None, index=False)
                    else:
                        pd.DataFrame(row).T.to_csv(
                            file_output, sep='\t', header=None, index=False)

            else:
                print("Warning : " + dir_results + library + "/" + library +
                      file_input_extension + " does not exist")
                print("Warning :  {" + dir_results + library + "/" + library +
                      file_input_extension + "} will not be filtered")

if __name__ == '__main__':
    main(sys.argv[1:])
