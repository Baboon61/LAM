# A script allowing to run FastQC software on all the librairies of
# metadata.txt

##############################IMPORTS##############################

import os
import sys
import getopt
import pandas as pd


def usage():
    print('Usage:')
    print('\tpython ' +
          sys.argv[0] + ' -m <metadata file> -d <preprocess directory> -o <output directory>')
    print('\t\t-h or --help : display this help')
    print('\t\t-m or --file_metadata : metadata file')
    print('\t\t-p or --dir_prep : preprocess directory')
    print('\t\t-o or --dir_output : output directory')


def main(argv):

    command_line = ""
    file_metadata = ""
    dir_prep = ""
    dir_output = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'm:p:o:', [
                                   'file_metadata=', 'dir_prep=', 'dir_output=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

##############################OPTIONS##############################

    if not opts:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-m', '--file_metadata'):
            file_metadata = arg
        elif opt in ('-p', '--dir_prep'):
            dir_prep = arg
        elif opt in ('-o', '--dir_output'):
            dir_output = arg
        else:
            print("Error : Bad option")
            usage()
            sys.exit(2)

##############################CHECK UP/SET UP##############################

    # CHECK METADATA FILE
    if file_metadata == "" or not os.path.exists(file_metadata):
        print("Error : You have to set a metadata file !\n")
        usage()
        sys.exit(2)

    # CHECK PREPROCESS DIRECTORY
    if not os.path.exists(dir_prep):
        print("Error : You have to set a preprocess directory !\n")
        usage()
        sys.exit(2)
    else:
        if dir_prep[-1] != "/":
            dir_prep += "/"

    # CHECK OUTPUT DIRECTORY
    if dir_output != "":
        if not os.path.exists(dir_output):
            os.system('mkdir ' + dir_output)
        if dir_output[-1] != "/":
            dir_output += "/"
    else:
        print("Error : You have to set an output directory !\n")
        usage()
        sys.exit(2)

##############################PRINTS##############################

    print('\n-----------------------------------------')
    print('Metadata file : ' + file_metadata)
    print('Preprocess directory : ' + dir_prep)
    print('Output directory : ' + dir_output)
    print('-----------------------------------------\n')

##############################PROGRAM##############################

    # READ METADATA FILE
    metadata = pd.read_table(file_metadata, sep='\t')

    # LOOP OVER EACH LIBRARIES
    for library in metadata['Library'].tolist():
        if not os.path.exists(dir_prep + library + "_R1.fq.gz"):
            print("Warning : " + dir_prep +
                  " does not contains {" + library + "_R1.fq.gz}")
            print("Warning :  {" + library + "} will not be filtered")
        elif not os.path.exists(dir_prep + library + "_R2.fq.gz"):
            print("Warning : " + dir_prep +
                  " does not contains {" + library + "_R2.fq.gz}")
            print("Warning :  {" + library + "} will not be filtered")
        else:
            # CREATE THE COMMAND LINE
            command_line = "fastqc " + dir_prep + library + \
                "_R1.fq.gz " + dir_prep + library + "_R2.fq.gz "
            print(command_line)
            os.system(command_line)
            command_line = "mv " + dir_prep + library + "_R1_fastqc.* " + dir_output
            print(command_line)
            os.system(command_line)
            command_line = "mv " + dir_prep + library + "_R2_fastqc.* " + dir_output
            print(command_line)
            os.system(command_line)

if __name__ == '__main__':
    main(sys.argv[1:])
