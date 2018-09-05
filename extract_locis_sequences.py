# A script allowing to extract sequences from Legitimate files according
# to their locus

##############################IMPORTS##############################

import os
import sys
import getopt
import os.path
import numpy as np
import pandas as pd


def usage():
    print('Usage:\n')
    print('\tpython ' +
          sys.argv[0] + ' -m <metadata file> -g <genome type> -p <postprocess directory> -l <locus file> [-i <mark input>]')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-g or --genome : only filter librairies results with this genome')
    print('\t\t-p or --dir_post : postprocess directory')
    print('\t\t-l or --file_locus : locus bed file (locus, start, end, strand, flag)')
    print('\t\t-i or --input_mark : mark from input file')


def main(argv):

    file_metadata = ""
    genome = ""
    dir_post = ""
    file_locus = ""
    input_mark = ""
    file_input_extension = ""
    file_output_extension = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:g:p:l:i:', [
                                   'file_metadata=', 'genome=', 'dir_post=', 'file_locus=', 'input_mark=', 'help'])
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
        elif opt in ('-p', '--dir_post'):
            dir_post = arg
        elif opt in ('-l', '--file_locus'):
            file_locus = arg
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

    # CHECK POSTPROCESS DIRECTORY
    if not os.path.exists(dir_post):
        print("Error : You have to set a postprocess directory !\n")
        usage()
        sys.exit(2)
    else:
        if dir_post[-1] != "/":
            dir_post += "/"

    # CHECK LOCUS FILE
    if file_locus == "" or not os.path.exists(file_locus):
        print("Error : You have to set a locus file !\n")
        usage()
        sys.exit(2)
    else:
        # READ LOCUS FILE
        df_locus = pd.read_table(file_locus, sep='\t', header=None)
        for index, row in df_locus.iterrows():
            try:
                if int(row[1]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.2 of your locus file !\n")
                    print("Error : Start position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.2 of your locus file !\n")
                print("Error : Unknown start position : " +
                      str(row[1]) + " !\n")
                usage()
                sys.exit(2)
            try:
                if int(row[2]) < 1:
                    print("Error : line." + str(index + 1) +
                          ", col.3 of your locus file !\n")
                    print("Error : End position has to be positive integer !\n")
                    usage()
                    sys.exit(2)
            except:
                print("Error : line." + str(index + 1) +
                      ", col.3 of your locus file !\n")
                print("Error : Unknown end position : " + str(row[2]) + " !\n")
                usage()
                sys.exit(2)
            if int(row[2]) < int(row[1]):
                print("Error : line." + str(index + 1) +
                      " of your locus file !\n")
                print("Error : End position is smaller than start position !\n")
                usage()
                sys.exit(2)
            if row[3] != '+' and row[3] != '-':
                print("Error : line." + str(index + 1) +
                      ", col.4 of your locus file !\n")
                print("Error : Unknown strand (+,-) : " + str(row[3]) + " !\n")
                usage()
                sys.exit(2)

    # CHECK INPUT MARKS HISTORY
    if input_mark == "":
        print("Warning : You will process the raw file !\n")

    # SELECT INPUT FILES
    if input_mark == "":
        file_input_extension = ".tlx"
    else:
        file_input_extension = "_" + "_".join(input_mark.split(",")) + ".tlx"

    # TEST IF INPUT EXIST IN AT LEAST ON LIBRARY
    check_input_mark = False
    for library in metadata['Library'].tolist():
        if os.path.exists(dir_post + library + "/" + library + "_Legitimate" + file_input_extension):
            check_input_mark = True
    if not check_input_mark:
        print("Error : Your input marks can not localize a good input file !\n")
        usage()
        sys.exit(2)

    # SELECT OUTPUT FILES
    if file_input_extension != "":
        file_output_extension = file_input_extension[:-4] + ".fasta"
    else:
        file_output_extension = ".fasta"

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Genome : ' + genome)
    print('Postprocess directory : ' + dir_post)
    print('Locus file : ' + file_locus)
    print('Input file extension: ' + file_input_extension)
    print('Output file extension : ' + file_output_extension)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # WRITE DATAFRAME ON ONE LINE
    pd.set_option('expand_frame_repr', False)

    # LOOP OVER EACH LIBRARIES
    for library in metadata['Library'].tolist():
        print(library)
        # CHECK DIRECTORY EXISTS
        if not os.path.exists(dir_post + library):
            print("Warning : " + dir_post +
                  " does not contains {" + library + "}")
            print("Warning :  {" + library + "} will not be filtered")
        else:
            # CHECK INPUT FILE EXISTS
            if os.path.exists(dir_post + library + "/" + library + "_Legitimate" + file_input_extension):
                df_legitimate = pd.read_csv(dir_post + library + "/" + library +
                                            "_Legitimate" + file_input_extension, sep='\t', header=0, index_col=None)
                df_legitimate = df_legitimate[
                    ['Qname', 'Locus Type', 'Position Type', 'Seq']]
                df_legitimate = df_legitimate[
                    df_legitimate['Position Type'] == "Inside"]

            else:
                print("Error : The Legitimate file for " +
                      library + " is missing !\n")
                usage()
                sys.exit(2)

            for index_locus, row_locus in df_locus.iterrows():
                print("-->" + row_locus[4] + " : " + str(
                    len(df_legitimate[df_legitimate['Locus Type'] == row_locus[4]])))
                # Open file
                f = open(dir_post + library + "/" + library + "_" + row_locus[4].replace(
                    "/", "-").replace(" ", "+").replace("'", "!") + file_output_extension, 'a')
                # WRITE SEQUENCES
                for index, row in df_legitimate[df_legitimate['Locus Type'] == row_locus[4]].iterrows():
                    f.write(">" + row['Qname'] + "\n" + row['Seq'] + "\n")

if __name__ == '__main__':
    main(sys.argv[1:])
