#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, threading

Description = """
Tool to build the suffixient of a text.
The input file cannot contain the characters 0, 1 or 2 which are
used internally by our construction algorithms.
"""

dirname         =  os.path.dirname(os.path.abspath(__file__))

bigbwt_dirname  =  os.path.join(dirname, "_deps/bigbwt-build")
src_dirname   =  os.path.join(dirname, "suff-set-src")

parse_exe       =  os.path.join(bigbwt_dirname, "pscan.x")
parsebwt_exe    =  os.path.join(bigbwt_dirname, "bwtparse")
parsebwt_exe64  =  os.path.join(bigbwt_dirname, "bwtparse64")
pfbwt_exe       =  os.path.join(bigbwt_dirname, "pfbwt.x")
pfbwtNT_exe     =  os.path.join(bigbwt_dirname, "pfbwtNT.x")
pfbwt_exe64     =  os.path.join(bigbwt_dirname, "pfbwt64.x")
pfbwtNT_exe64   =  os.path.join(bigbwt_dirname, "pfbwtNT64.x")

pfp_exe         =  os.path.join(src_dirname,  "pfp_suffixient")
one_pass_exe    =  os.path.join(src_dirname,  "one-pass")
linear_exe      =  os.path.join(src_dirname,  "linear-time")
fm_exe          =  os.path.join(src_dirname,  "fm")

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('-a', '--algorithm', 
                      help='suffixient set construction algorithm (one-pass,linear,fm,PFP)',
                      default="one-pass", type=str)
  parser.add_argument('-o', '--output', 
                      help='output files basepath (def. input)', default="", type=str)

  parser.add_argument('-w', '--wsize',  
                      help='PFP: sliding window size (def. 10)',     default=10,  type=int)
  parser.add_argument('-p', '--mod',    
                      help='PFP: hash modulus (def. 100)',           default=100, type=int)
  parser.add_argument('-t', '--threads',
                      help='PFP: number of threads (def. 1)', default=1,   type=int)
  parser.add_argument('-i', '--no-invert',
                      help='PFP: do not invert the text before running the algorithm', action='store_true')
  args = parser.parse_args()

  # define output basepath
  if args.output == "": args.output = args.input

  logfile_name = args.input + ".suffixient.log"
  # get main directory
  args.bigbwt_dir = os.path.split(sys.argv[0])[0]
  print("Sending logging messages to file:", logfile_name)

  with open(logfile_name,"a") as logfile:

    print("Computing a smallest suffixient set of:",args.input)
    start = time.time()

    if args.algorithm == "one-pass":

      command = "{exe} -o {ofile}".format(
                exe = one_pass_exe,
                ofile=args.output+".suff")
      execute_command_stdin(args.input,command,logfile,logfile_name)

    elif args.algorithm == "linear":

      command = "{exe} -o {ofile}".format(
                exe = linear_exe,
                ofile=args.output+".suff")
      start = time.time()
      execute_command_stdin(args.input,command,logfile,logfile_name)

    elif args.algorithm == "fm":

      command = "{exe} -o {ofile}".format(
                exe = fm_exe,
                ofile=args.output+".suff")
      start = time.time()
      execute_command_stdin(args.input,command,logfile,logfile_name)

    elif args.algorithm == "PFP":

      if not args.no_invert:
        text = ""
        with open(args.input,"r") as file:
          text = file.read()
        text = reversed(text)
        args.input += ".inv"
        with open(args.input,"w+") as file:
          file.write("".join(text))

      # ---------- parsing of the input file
      command = "{exe} {file} -w {wsize} -p {modulus} -t {th}".format(
              exe = os.path.join(args.bigbwt_dir,parse_exe),
              wsize=args.wsize, modulus = args.mod, th=args.threads, file=args.input)

      command += " -s"
      print("==== Parsing. Command:", command)
      if(execute_command(command,logfile,logfile_name)!=True):
        return
      print("Elapsed time: {0:.4f}".format(time.time()-start))

      # ---- delete intermediate files
      delete_temp_files(args,logfile,logfile_name)

      # ---- run suffixient construction
      command = "{exe} -i {file} -w {wsize} -n {size}".format(
              exe = os.path.join(args.bigbwt_dir,pfp_exe),
              file = args.input, wsize = args.wsize, size =  os.path.getsize(args.input)+1)
      command += " -o {out_file}".format(out_file=args.output+".suff")
      #if args.o != "":
      #  command += " -o {out_file}".format(out_file=args.o)
      #if args.c:
      #  command += " -p"
      #if args.r:
      #  command += " -r"
      print("==== Compute suffixient set. Command:", command)
      if(execute_command(command,logfile,logfile_name)!=True):
        return
      #subprocess.run(command.split())

  elapsed_time = time.time() - start

  print("The resulting suffixient set was sent to: "+args.output+".suff")
  print("### Elapsed time: {time} seconds".format(time=elapsed_time))

######### AUXILIARY FUNCTIONS #########

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name,env=None):
  try:
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

# delete intermediate files
def delete_temp_files(args,logfile,logfile_name):
    #if args.k==False:
    if True:
      print("==== Deleting temporary files.") # no need to show the command
      command = "rm -f {file}.parse_old {file}.last {file}.bwlast {file}.ilist".format(file=args.input)
      #command = "rm -f {file}.parse {file}.parse_old {file}.last {file}.bwlast {file}.dict {file}.ilist {file}.occ".format(file=args.input)
      if(execute_command(command,logfile,logfile_name)!=True):
        return
      for i in range(args.threads):
        command = "rm -f {file}.{i}.parse_old {file}.{i}.last".format(file=args.input, i=i)
        if(execute_command(command,logfile,logfile_name)!=True):
          return
      
      command = "rm -f {file}.sai {file}.bwsai".format(file=args.input);
      if(execute_command(command,logfile,logfile_name)!=True):
        return
      for i in range(args.threads):
        command = "rm -f {file}.{i}.sai".format(file=args.input, i=i)
        if(execute_command(command,logfile,logfile_name)!=True):
          return

##########################
if __name__ == '__main__':
    main()