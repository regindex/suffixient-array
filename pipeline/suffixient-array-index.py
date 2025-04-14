#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, threading

Description = """
Tool to build the suffixient-array index of a text.
The input file cannot contain the characters 0, 1 or 2 which are
used internally by our construction algorithms.
If opt-sA is selected, the input must contain DNA characters (A,C,G,T) only.
"""

dirname          =  os.path.dirname(os.path.abspath(__file__))

bigbwt_dirname   =  os.path.join(dirname, "_deps/bigbwt-build")
suff_src_dirname =  os.path.join(dirname, "suff-set-src")
sA_src_dirname   =  os.path.join(dirname, "sA-index-src")

step1_exe   =  os.path.join(suff_src_dirname, "one-pass-build-index")
step2_exe   =  os.path.join(sA_src_dirname,   "build_store_sA_index")
locate_exe  =  os.path.join(sA_src_dirname,   "locate")
mems_exe    =  os.path.join(sA_src_dirname,   "mems")

mapping   = { "sA": ".sA", "opt-sA": ".opt_sA", "PA": ".pa" }
mapping_2 = { "sA": "suffixient-array", "opt-sA": "elias-fano-opt", "PA": "prefix-array" }

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input files basepath', type=str)
  parser.add_argument('-v', '--index-variant', 
                      help='suffixient-array index variant: sA(def.)|opt-sA|PA',
                      default="sA", type=str)
  parser.add_argument('-o', '--oracle-variant', 
                      help='random access text oracle variant: lz77(def.)|rlz|bitpacked-text',
                      default="lz77", type=str)
  parser.add_argument('-b', '--build-index', 
                      help='build suffixient-array index', action='store_true')
  parser.add_argument('-e', '--extra-ef-space', 
                      help='opt-sA: maximum allowed extra space (percentage) on top of the suffixient-array: Def. 30',
                      default="30", type=str)
  parser.add_argument('-l', '--locate-one-occ', 
                      help='run locate one occurrence queries', action='store_true')
  parser.add_argument('-m', '--find-mems', 
                      help='run find MEMs queries', action='store_true')
  parser.add_argument('-p', '--pattern-file', 
                      help='file containing the patterns to locate (fasta format required)', type=str)
  #parser.add_argument('-o', '--output', 
  #                    help='output files basepath (def. input)', default="", type=str)
  args = parser.parse_args()

  logfile_name = args.input + ".suffixient-array.log"
  # get main directory
  args.bigbwt_dir = os.path.split(sys.argv[0])[0]
  print("Sending logging messages to file:", logfile_name)

  start = time.time()

  with open(logfile_name,"a") as logfile:

    if args.build_index:

      # check for the input files
      if not os.path.exists(args.input):
        print("Error: the input file is missing... exiting.")
        exit(1)

      print("Computing a suffixient-array index of:",args.input)

      command = "{exe} -t {itype} -o {ofile}".format(
                exe   = step1_exe,
                itype = args.index_variant,
                ofile = args.input)
      print("Running step 1 (suffixient set computation) on:",args.input)
      print("### command:",command)
      if(execute_command_stdin(args.input,command,logfile,logfile_name)!=True):
        return

      command = "{exe} -i {ifile} -t {itype} -o {oracle} -l {maxes}".format(
                exe    = step2_exe,
                ifile  = args.input,
                itype  = args.index_variant,
                oracle = args.oracle_variant,
                maxes  = args.extra_ef_space)
      print("Running step 2 (suffixient-array index computation) on:",args.input)
      print("### command:",command)
      if(execute_command(command,logfile,logfile_name)!=True):
        return

      print("The resulting suffixient-array index was sent to: "+args.input+mapping[args.index_variant])

      # ---- delete temp files
      delete_temp_files(args,logfile,logfile_name)

    if args.locate_one_occ:

      # check for the input files
      if( not os.path.exists(args.input+mapping[args.index_variant]) or
          not os.path.exists(args.input+"."+args.oracle_variant) or
          not os.path.exists(args.input) or
          not os.path.exists(args.pattern_file)):
        print("Error: some of the input files are missing... exiting.")
        exit(1)

      print("Locating one occurrence of the patterns in:",args.pattern_file)
      command = "{exe} -i {ifile} -t {itype} -o {oracle} -p {pattf}".format(
                exe    = locate_exe,
                ifile  = args.input,
                itype  = mapping_2[args.index_variant],
                oracle = args.oracle_variant,
                pattf  = args.pattern_file)
      print("### command:",command)
      if(execute_command(command,logfile,logfile_name)!=True):
        return

      print("The results were sent to: "+args.pattern_file+".occs")

    if args.find_mems:

      # check for the input files
      if( not os.path.exists(args.input+mapping[args.index_variant]) or
          not os.path.exists(args.input+"."+args.oracle_variant) or
          not os.path.exists(args.input) or
          not os.path.exists(args.pattern_file)):
        print("Error: some of the input files are missing... exiting.")
        exit(1)

      print("Finding MEMs of the patterns in:",args.pattern_file)
      command = "{exe} -i {ifile} -t {itype} -o {oracle} -p {pattf}".format(
                exe    = mems_exe,
                ifile  = args.input,
                itype  = mapping_2[args.index_variant],
                oracle = args.oracle_variant,
                pattf  = args.pattern_file)
      print("### command:",command)
      if(execute_command(command,logfile,logfile_name)!=True):
        return

      print("The results were sent to: "+args.pattern_file+".mems")

    elapsed_time = time.time() - start
    print("### Elapsed time: {time} seconds".format(time=elapsed_time))

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name,env=None):
  try:
    #subprocess.run(command.split(),stdout=logfile,stderr=logfile,check=True,env=env)
    subprocess.check_call(command.split(),stdout=logfile,stderr=logfile,env=env)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True

def execute_command_stdin(input_file,command,logfile,logfile_name,env=None):
  try:
    with open(input_file,"r") as i:
      subprocess.check_call(command.split(),stdin=i,stdout=logfile,stderr=logfile,env=env)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True

def delete_temp_files(args,logfile,logfile_name):

    print("==== Deleting temporary files.") # no need to show the command

    if args.index_variant == "PA":
      command = "rm -f {file}.pa".format(file=args.input)
      if(execute_command(command,logfile,logfile_name)!=True):
        return

    if args.index_variant == "sA" or args.index_variant == "opt-sA":
      command = "rm -f {file}.suff {file}.mult {file}.lcs".format(file=args.input)
      if(execute_command(command,logfile,logfile_name)!=True):
        return

    if args.oracle_variant == "lz77":
      command = "rm -f {file}.alph {file}.alph.sa {file}.alph.start " \
          "{file}.alph.start {file}.alph.len {file}.alph.char".format(file=args.input)
      if(execute_command(command,logfile,logfile_name)!=True):
        return


##########################
if __name__ == '__main__':
    main()