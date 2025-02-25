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
    #parser.add_argument('--linear',  help='test linear time algorithm (def. False)',action='store_true')
    #parser.add_argument('--PFP',  help='test PFP algorithm (def. False)',action='store_true')
    #parser.add_argument('--pattern_file', help='File path containing the patterns', type=str, required=True)
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
        header = "index,dataset,dataset_size(bytes),pattern_length,no_patterns,time(sec),time_per_patt(microsec),time_per_char(nanosec),space(bytes)\n"
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
            '''
            with open(args.pattern_file,'r') as patt:
                line = patt.readline()
                line = patt.readline()
                pattern_length = int(len(line))-1
            print("Pattern length =",pattern_length)
            '''

            if args.locate_sA:
                #for ot in [["plain",""],["bitpacked",""],["lz77",".lz77"],["Hk",".hkfv"]]:
                for ot in [["bitpacked",""]]:
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

                        peak_memory = output_str[17]
                        tot_time = output_str[26]
                        time_per_pattern = output_str[42]
                        time_per_character = output_str[48]

                        csv_line = it[0]+"+"+ot[0] + "," + args.input_file + "," + str(dataset_size) + "," + str(pattern_length) + \
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



            '''
            if args.one_pass:
                command = "{exe} -p -o {output}".format(exe = one_pass_exe, output = args.input_file+".suffixient")
                print("#######",command,"<",args.input_file)
                # Start writing the csv line
                csv_preamble = "one-pass," + args.input_file.split('/')[-1] + "," + str(dataset_size) + ","
                # Run the testing command
                run_test(command,args.input_file,log_file_path,log_csv_path,res,log_res_path,csv_preamble,True)
                # Remove output files
                os.remove(args.input_file+".suffixient")
            if args.linear:
                command = "{exe} -p -o {output}".format(exe = linear_time_exe, output = args.input_file+".suffixient")
                print("#######",command,"<",args.input_file)
                # Start writing the csv line
                csv_preamble = "linear," + args.input_file.split('/')[-1] + "," + str(dataset_size) + ","
                # Run the testing command
                run_test(command,args.input_file,log_file_path,log_csv_path,res,log_res_path,csv_preamble,True)
                # Remove output files
                os.remove(args.input_file+".suffixient")
            '''

def run_test(command, input_file, log_file, csv_file, resfile, res_file, csv_preamble, stdIN = False):
    '''
    Start the testing procedure for a specified command and store the results.
    This function executes a command with a given input file, captures its output, 
    and logs the results in a csv table.

    Args:
        command (str): the command line to execute.
        input_file (str): the path to the input file.
        log_file (str): the path to the log file where the execution logs are recorded.
        csv_file (str): the path to the CSV file where all collected results are stored.
        resfile (pointer): a file pointer to the currently opened CSV result file.
        res_file (str): the path to the file that stores time and space results for the current input.
        csv_preamble (str): a string containing the preamble for the current entry in the CSV file.
        stdIN (bool): A boolean indicating whether the current executable expects input to be streamed 
                      from standard input (True) or read from a file (False).

    Returns:
        None
    '''
    # Execution elapsed time
    elapsed_time = 0.0
    # Clear and load cache associated to a file
    manage_file_cache(input_file)

    # Open the log file in append mode
    with open(log_file,"a") as logfile:
        # Write the command to get the execution time/space with /usr/bin/time
        command = "{bin} -v -o {info} ".format(bin=time_space_bin,info=res_file) + command
        # Start command execution
        print("Command:", command, flush=True)
        start = time.time()
        if stdIN:
            ret = execute_command_stdin(command,input_file,logfile,log_file,timeout,None)
        else:
            ret = execute_command(command,logfile,log_file,timeout,None)
        # Check if the execution ended correctly
        if not ret[0]:
            print("Some error occurred during the execution...Check the log file:",log_file)
            exit(1)
        # Print the elapsed time
        elapsed_time = time.time()-start
        print("Elapsed time: {0:.4f}".format(elapsed_time), flush=True)
        logfile.write("Elapsed time: {0:.4f}".format(elapsed_time))
        # End command execution

        # Extract time/space from resfile
        WCT,rss,cpu,cpu_eff = get_time_space(res_file);
        # Complete the CSV line and store it in the table storing the results
        csv_preamble += str(ret[1]) + "," + str(WCT) + "," + str(rss/10**3) + "," + str(cpu_eff) + "\n"
        resfile.write(csv_preamble) 
        resfile.flush()

    return 

# execute command: return True is everything OK, False otherwise, then return the execution stdout
def execute_command(command,logfile,logfile_path,timeout=None,env=None):
    '''
    Execute command and returns information given through the standard output.

    Args:
        command (str): command to execute.
        logfile (pointer): pointer to a currently opened logfile.
        logfile_path (str): path to a currently opened logfile.
        timeout (int): time in seconds before stopping the execution for timeout.
        env (environ): custom execution environment.

    Returns:
        bool: boolen equals to False if some error happened, True otherwise.
        int: dummy int.
    '''
    try:
        p = subprocess.Popen(command.split(),stdout=logfile,stderr=logfile,env=env)
        if timeout != None:
            try:
                p.communicate(timeout=timeout)
            except subprocess.TimeoutExpired as e:
                print("Timeout executing command line:")
                print("\t"+ command)
                # get all descendant pids
                pids_o = subprocess.check_output("pstree -p {pid} | grep -o '([0-9]\\+)' | grep -o '[0-9]\\+'".format(pid=p.pid),shell=True)
                pids = re.findall('\\d+', pids_o.decode("utf-8"))
                print("Killing PIDs: " + "".join(str(pid) + " " for pid in pids))
                print("Check log file: " + logfile)
                for pid in pids:
                    os.kill(int(pid),signal.SIGKILL)
                    os.kill(int(pid),signal.SIGTERM)
                return False,-1

    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile)
        return False,-1

    return True,0

def execute_command_stdin(command,inputfile,logfile,logfile_path,timeout=None,env=None):
    '''
    Execute command and returns information given through the standard output.
    This function executes binaries which expect the input coming as a stream from
    the standard input.

    Args:
        command (str): command to execute.
        inputfile (pointer): path to the input file.
        logfile (str): pointer to a currently opened logfile.
        logfile_path (str): path to a currently opened logfile.
        timeout (int): time in seconds before stopping the execution for timeout.
        env (environ): custom execution environment.

    Returns:
        bool: boolen equals to False if some error happened, True otherwise.
        int: Suffixient set size.
    '''
    try:
        with open(inputfile, 'rb', 0) as a:
            p = subprocess.Popen(command.split(),stdin=a,stdout=subprocess.PIPE,stderr=logfile,env=env)

        if timeout != None:
            try:
                stdout = p.communicate(timeout=timeout)[0].decode()
                logfile.write(stdout)
                suffixient_set_size = stdout.split(' ')[-1].split('\n')[0]
            except subprocess.TimeoutExpired as e:
                print("Timeout executing command line:")
                print("\t"+ command)
                # get all descendant pids
                pids_o = subprocess.check_output("pstree -p {pid} | grep -o '([0-9]\\+)' | grep -o '[0-9]\\+'".format(pid=p.pid),shell=True)
                pids = re.findall('\\d+', pids_o.decode("utf-8"))
                print("Killing PIDs: " + "".join(str(pid) + " " for pid in pids))
                print("Check log file: " + logfile)
                for pid in pids:
                    os.kill(int(pid),signal.SIGKILL)
                    os.kill(int(pid),signal.SIGTERM)
                return False,-1

    except subprocess.CalledProcessError as e:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile)
        return False,-1
    
    return True,suffixient_set_size

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

def get_time_space(resfile):
    '''
    Extract execution time/space information from gtime or \\usr\\bin\\time output.

    Args:
        resfile (str): path to the file containing the gtime or \\usr\\bin\\time output.

    Returns:
        int: Wall clock time in seconds.
        int: Resident set size in bytes.
        int: CPU time in seconds.
        str: CPU efficiency percentage.
    '''
    # Get resident set size from res file
    rss_o = subprocess.check_output("awk -F ':' '{{print $2}}' {log} | sed -n '10p'".format(log=resfile),shell=True)
    rss = int(re.search(r'\d+',rss_o.decode("utf-8")).group())
    # Get CPU time
    cpu_o = subprocess.check_output("awk -F ':' '{{print $2}}' {log} | sed -n '2p'".format(log=resfile),shell=True)
    cpu = cpu_o.decode("utf-8")
    cpu = cpu[1:-1]
    # Get CPU efficiency
    cpu_o = subprocess.check_output("awk -F ':' '{{print $2}}' {log} | sed -n '4p'".format(log=resfile),shell=True)
    cpu_eff = cpu_o.decode("utf-8")
    cpu_eff = cpu_eff[1:-1]
    # Get Wall clock time
    wall_clock_o = subprocess.check_output("grep \"Elapsed\" {log}".format(log=resfile), shell=True)
    wall_clock_time = wall_clock_o.decode("utf-8").split('(h:mm:ss or m:ss):')[-1].split('\n')[0]
    wall_clock_time = time_to_seconds(wall_clock_time)

    return (wall_clock_time,rss,cpu,cpu_eff)

def time_to_seconds(time_str):
    '''
    Convert a time string in ormat h:mm:ss or m:ss format to seconds.

    Args:
        time_str (str): the time string to convert.

    Returns:
        int: the total number of seconds.
    '''
    # Split the time string by colon
    parts = time_str.split(':')
    
    # Initialize seconds variable
    time_in_seconds = 0
    
    # Process the parts depending on their length
    if len(parts) == 3:  # h:mm:ss format
        hours = int(parts[0])
        minutes = int(parts[1])
        seconds = float(parts[2])  # float captures decimal seconds
        time_in_seconds = hours * 3600 + minutes * 60 + seconds
    elif len(parts) == 2:  # m:ss format
        minutes = int(parts[0])
        seconds = float(parts[1])  # float captures decimal seconds
        time_in_seconds = minutes * 60 + seconds
    else:
        raise ValueError("Invalid time format. h:mm:ss or m:ss format expected.")
    
    return time_in_seconds

if __name__ == '__main__':
    main()