#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math, platform, re, decimal

Description = """
Pipeline for testing the suffixient index locate
"""
base_path = os.path.dirname(os.path.abspath(__file__))

## BINARY FILES DECLARATION
build_sA_exe    = base_path + "/build/index/build_store_index"
locate_sA_exe    = base_path + "/build/index/locate"
locate_ri_exe    = base_path + "/competitors/r-index/build/ri-locate"
build_ri_exe     = base_path + "/competitors/r-index/build/ri-build"

# TIMEOUR DECLARATION
timeout = 86400 # 24h

# SET SOFTWARE TO GET TIME/SPACE BASED ON OS
if platform.system() == "Darwin":
    time_space_bin = "gtime"
else:
    time_space_bin = "\\usr\time\bin"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_file', nargs='?', help='the input file', type=str)
    parser.add_argument('--locate_sA',  help='test locate one occurrence queries (def. False)',action='store_true')
    parser.add_argument('--locate_ri',  help='test locate one occurrence queries (def. False)',action='store_true')
    parser.add_argument('--logs_dir_name', help='Define the directory name containing the logs (no default)', type=str, required=True)
    args = parser.parse_args()

    # Specify the log directory path
    log_dir_path = os.path.join(base_path, args.logs_dir_name)
    # Check if the directory exists
    if not os.path.exists(log_dir_path):
        os.makedirs(log_dir_path) 
        print("Logs directory created successfully!")
    # Specify the logfile and result files paths
    log_file_path = os.path.join(log_dir_path, args.input_file.split('/')[-1] + ".log")
    log_res_path  = os.path.join(log_dir_path, args.input_file.split('/')[-1] + ".res")
    log_csv_path  = os.path.join(log_dir_path, args.input_file.split('/')[-1] + ".csv")

    # Open the file containing the table storing the results
    with open(log_csv_path,"w+") as res:

        # CSV file header
        header = "index,oracle,dataset,dataset_size,pattern_length,no_patterns,time,time_per_patt,time_per_char,Memory peak\n"
        res.write(header)
        # Compute the dataset size
        command = "wc -c {file}".format(file = args.input_file)
        dataset_size = str(subprocess.check_output(command.split()).decode()).split(' ')[0]
        print("Dataset size (in bytes) =",dataset_size)

        pattern_lengths = [10,100,1000]

        for pattern_length in pattern_lengths:
            args.pattern_file = args.input_file + "_pat" + str(pattern_length) + ".fasta"
            print(args.pattern_file)

            command = "wc -l {file}".format(file = args.pattern_file)
            no_patterns = str(subprocess.check_output(command.split()).decode()).split(' ')[0]
            no_patterns = int(no_patterns)
            if(no_patterns%2 != 0):
                no_patterns += 1
            no_patterns = int(no_patterns/2)
            print("No patterns =",no_patterns)

            if args.locate_sA:
                #for ot in [["plain",""],["bitpacked",""],["lz77",".lz77"],["Hk",".hkfv"]]:
                for ot in [["rlz",".rlz"]]:
                    for it in [["prefix-array",".pai"],["baseline",".bai"],["elias-fano",".efi"]]:

                        if not os.path.exists(args.input_file + it[1]) and not os.path.exists(args.input_file + ot[1]):
                            print("Index file not detected! Computing the",it[0],"index for",args.input_file)
                            command = "{exe} -i {input_path} -t {it} -o {ot}".format(exe = build_sA_exe, input_path =args.input_file, it = it[0], ot = ot[0])
                            print(command)
                            subprocess.check_output(command.split())
                        if not os.path.exists(args.input_file + it[1]) and os.path.exists(args.input_file + ot[1]):
                            print("Index file not detected! Computing the",it[0],"index for",args.input_file,"However the oracle data structure was detected")
                            command = "{exe} -i {input_path} -t {it} -o plain".format(exe = build_sA_exe, input_path =args.input_file, it = it[0])
                            print(command)
                            subprocess.check_output(command.split())
                        if os.path.exists(args.input_file + it[1]) and not os.path.exists(args.input_file + ot[1]):
                            print("Index file not detected! Computing the",it[0],"index for",args.input_file,"However the oracle data structure was detected")
                            command = "{exe} -i {input_path} -t prefix-array -o {ot}".format(exe = build_sA_exe, input_path =args.input_file, ot = ot[0])
                            print(command)
                            subprocess.check_output(command.split())


                        command = "{exe} -i {base_path} -t {index_variant} -o {oracle_type} -p {pattern_file}".format(
                                   exe = locate_sA_exe, base_path = args.input_file, index_variant = it[0], oracle_type = ot[0], pattern_file = args.pattern_file)
                        print("#####",command)
                        manage_file_cache(args.input_file)
                        output_str = str(subprocess.check_output(command.split())).split(' ')
                        print(output_str)

                        peak_memory = output_str[19]
                        tot_time = output_str[28]
                        time_per_pattern = output_str[44]
                        time_per_character = output_str[50]
                        #print("-->",peak_memory,tot_time,time_per_pattern,time_per_character)

                        csv_line = it[0] + "," + ot[0] + "," + args.input_file + "," + str(dataset_size) + "," + str(pattern_length) + \
                                   "," + str(no_patterns) + "," + str(decimal.Decimal(tot_time)) + "," + str(time_per_pattern) + \
                                   "," + time_per_character + "," + peak_memory + "\n"
                        res.write(csv_line) 
                        res.flush()

            if args.locate_ri:

                if not os.path.exists(args.input_file + ".ri"):
                    print("Index file not detected! Computing the r-index for",args.input_file)
                    command = "{exe} {input_path}".format(exe = build_ri_exe, input_path = args.input_file)
                    subprocess.check_output(command.split())

                command = "{exe} {index_path} {pattern_path}".format(
                            exe = locate_ri_exe, index_path = args.input_file + ".ri", pattern_path = args.pattern_file)
                manage_file_cache(args.input_file)

                print(command)
                output_str = str(subprocess.check_output(command.split())).split(' ')

                peak_memory = output_str[12]
                tot_time = output_str[21]
                time_per_pattern = output_str[37]
                time_per_character = output_str[43]

                csv_line = "r-index," + args.input_file + "," + str(dataset_size) + "," + str(pattern_length) + \
                           "," + str(no_patterns) + "," + str(decimal.Decimal(tot_time)) + "," + str(time_per_pattern) + \
                           "," + time_per_character + "," + peak_memory + "\n"
                res.write(csv_line) 
                res.flush()

def manage_file_cache(filename):
    '''
    Remove cached data for the input file.

    Args:
        filename (str): path to the input file.

    Returns:
        None

    IMPORTANT: This function may not work on certain MacOS
    '''
    try:
        # Open the file in read mode
        with open(filename, "r") as fd:
            # Get the length of the file
            length = os.path.getsize(filename)
            
            # Flush any buffered data to disk
            os.fdatasync(fd.fileno())
            #os.fsync(fd.fileno())
            
            # Advise the OS about dropping and loading cache
            os.posix_fadvise(fd.fileno(), 0, length, os.POSIX_FADV_DONTNEED)
            #fcntl.fadvise(fd.fileno(), 0, length, os.POSIX_FADV_DONTNEED)
            os.posix_fadvise(fd.fileno(), 0, length, os.POSIX_FADV_WILLNEED)
            #fcntl.fadvise(fd.fileno(), 0, length, os.POSIX_FADV_DONTNEED)
            
            # Advise the OS about reading the file sequentially
            # os.posix_fadvise(fd.fileno(), 0, length, os.POSIX_FADV_SEQUENTIAL)

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    main()