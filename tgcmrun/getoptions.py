import argparse,sys,os
from classes import Run,Job
from utils import getenv

#-----------------------------------------------------------------------
def get_args():
  description='''\
      Make selected benchmark or test runs of the TIEGCM/TIMEGCM models.

      For each run, a fortran namelist input file and a csh job script are generated, 
      and the job script is optionally executed to submit the run. The user selects 
      from a list of available runs, and provides the necessary information through 
      command-line options, interactive prompts, or environment variables.

      This Python application can be run interactively by typing "tgcmrun", and 
      optionally providing command-line arguments (-arg=value pairs). The user will 
      be prompted for any required arguments that are not on the command-line. 

      The "required" arguments are: -run_name or -run_number, -model_name, and -model_res. 
      If at least these three arguments are provided, AND the four below environment 
      variables are set, then everything else will default, and the user will only be 
      prompted for submitting the job. This allows multiple executions from a shell 
      script, see below example (the -submit=yes argument means submit the job without 
      prompting the user):

      tgcmrun -run_name=decsol_smin -model_name=tiegcm -model_res=5.0 -submit=yes
      tgcmrun -run_name=decsol_smax -model_name=tiegcm -model_res=5.0 -submit=yes
      tgcmrun -run_name=decsol_smin -model_name=tiegcm -model_res=2.5 -submit=yes
      tgcmrun -run_name=decsol_smax -model_name=tiegcm -model_res=2.5 -submit=yes

      For the above shell script to work, the user should have these four environment variables set:
        TGCMTEMP:     Large temporary directory where the model can be built, executed, and output stored.
        TGCMDATA:     Directory containing data files required by the model (netcdf data and history files)
        TIEGCM_ROOT:  Path to the tiegcm model root directory containing source code, scripts, tgcmrun, etc.
        TIMEGCM_ROOT: Path to the timegcm model root directory containing source code, scripts, tgcmrun, etc.

      A word of caution: this code was developed with Python 2.7.7 on the NCAR machine yellowstone,
      and may be sensitive to version change (As of Nov 2015, it has not been tested with Python 3.x)
'''
  parser = argparse.ArgumentParser(
    description=description,
    formatter_class=argparse.RawTextHelpFormatter)

  help_run_name   = "Short name of run (from run list)"
  help_run_number = "Number of run (from run list)"
  help_model_name = "Model name (either 'tiegcm' or 'timegcm')"
  help_model_res  = "Model resolution (either 5.0 or 2.5 degrees)"
  help_model_root = "Model root directory (default: env vars TIEGCM_ROOT or TIMEGCM_ROOT)"
  help_machine    = "Machine or platform (either 'ys' (yellowstone) or 'linux' (generic Linux))"
  help_execdir    = "Directory where model will be built and executed (default: env var TGCMTEMP)"
  help_tgcmdata   = "Path to data files needed by the model (default: env var TGCMDATA)"
  help_nprocs     = "Number of processors (total MPI tasks)"
  help_project    = "Authorized NCAR project number, e.g.: #BSUB -P P28100036 (ys only)"
  help_queue      = "LSF queue name, e.g.: #BSUB -q regular (ys only)"
  help_wc         = "Wallclock limit e.g.: '01:30' is 1 hour, 30 minutes (ys only)"
  help_step       = "Model timestep (seconds) default=60 for res5.0, default=30 for res2.5"
  help_submit     = "Submit job without prompting user? (yes/no))"
  help_execute    = "Execution flag for job script ('yes'/'no') (default: 'yes')"
  help_stdout_dir = "Directory in which stdout files and run scripts are to be saved (default: cwd)"
  help_hist_dir   = "Directory in which output history files are to be saved (default: cwd)"
  help_compiler   = "Compiler to be used on Linux desktop platform (not supercomputer)\n(valid values: 'intel', 'pgi', 'gfort') (default 'intel')"

  arg_list = [
    ['run_name'   , help_run_name],
    ['run_number' , help_run_number],
    ['model_name' , help_model_name],
    ['model_res'  , help_model_res],
    ['model_root' , help_model_root],
    ['machine'    , help_machine],
    ['compiler'   , help_compiler],
    ['execdir'    , help_execdir],
    ['execute'    , help_execute],
    ['tgcmdata'   , help_tgcmdata],
    ['nprocs'     , help_nprocs],
    ['project'    , help_project],
    ['queue'      , help_queue],
    ['wc'         , help_wc],
    ['step'       , help_step],
    ['submit'     , help_submit],
    ['stdout_dir' , help_stdout_dir],
    ['hist_dir'   , help_hist_dir]]

  for arg in arg_list:
    parser.add_argument('-'+arg[0], help=arg[1])

  return parser.parse_args()

#-----------------------------------------------------------------------
def get_options(arg,run,job,option):
#
# arg is the needed argument (string).
# 'run' and 'job' are partially defined objects (if promting, order of 
# obtaining options is important because of dependencies).
#
# If true (non-None), "option" is the argument value from the command line.
# If option is false (not provided on command line), then the user is
#   prompted, and/or env vars are used to get a value for the arg.
#
# If run_name is not on command line, get it by printing a list, and 
# prompting the user for the run number:
# 
  if arg == 'run_number':
    if run.name:            # use this if it was set on the command line
        return run.get_number(run.name)
    if option:              # check for command-line argument
      if not run.validate_number(option):
        print '>>> Bad run number found on command-line: ',option,' (Must be an integer >= 0 and <= ',run.nruns-1,')'
        sys.exit()
      else:
        return option
    else:                    # prompt user for run number
      run.print_runs()
      answer = ''
      while not run.validate_number(answer):
#
# Note: type returned by raw_input is always a string.
# run.number has to be an integer (it will used to index into lists, etc)
#
        answer = raw_input("\nEnter number of desired run (0-"+str(run.nruns-1)+") ('q' to quit, 'p' to print list, default=0): ")
        if answer == 'q': 
          sys.exit()        # quit the app
        elif answer == '':
          return 0 # accept default
        elif answer == 'p': 
          run.print_runs()
          answer = ''       # go back up to prompt ("while not" loop) w/o printing list
        else:
          if run.validate_number(answer): # if this succeeds, it returns integer (not a string)
            return answer
          else:
            print '>>> Bad run number: ',answer,' (Must be an integer >= 0 and <= ',run.nruns-1,')'
            answer = ''     # go back up to prompt ("while not" loop) w/o printing list

  elif arg == 'run_name':
    if option:
      if run.validate_name(option):
        run.name = option
        return run.name
      else:
        print '>>> Invalid run name: ',option,' found on command line.'
        sys.exit()
#
# Get run.model_name:
#
  elif arg == 'model_name':
    if option:
      if option != 'tiegcm' and option != 'timegcm':
        print ">>> Unknown model name on the command line: '",option,"' (must be 'tiegcm' or 'timegcm')"
        sys.exit()
      run.model_name = option
    else:
      job.user = os.getlogin() 
      answer = ''
      while answer != 'tiegcm' and answer != 'timegcm':
        answer = raw_input("Run "+run.name+": Enter model name ('tiegcm' or 'timegcm', default=tiegcm): ")
        if answer == 'q':
          sys.exit()
        elif answer == '':
          answer = 'tiegcm'
          continue
        elif answer != 'tiegcm' and answer != 'timegcm':
          print ">>> Unknown model name ",answer," (must be 'tiegcm' or 'timegcm')"
          get_options(arg,run,job,option)
      run.model_name = answer
    return run.model_name
#
# Get model root directory (both model_root and model_root/scripts have to exist
# as directories.
#   
#   1. If model_root set on command-line, check that it exists as a directory, and
#      check that model_root/scripts exists as a directory. If these fail, exit the program.
#   2. If not on command-line, check for env var TIEGCM_ROOT or TIMEGCM_ROOT.
#      If they are set and meet the existence checks of part 1, return the value.
#   3. If not on command-line and env vars are not set, prompt the user.
#
# Check for command-line option:
#
  elif arg == 'model_root':
    if option:
      if not os.path.isdir(option):
        print '>>> Cannot find model_root directory provided on the command line: ',option
        sys.exit()
      else:
        scripts_dir = option+'/scripts'
        if not os.path.isdir(scripts_dir):
          print '>>> Model root ',option,' from command line exists, but scripts directory ',scripts_dir,' does not exist.'
          sys.exit()
      return option
#
# Check for env var:
#
    else:
      if job.model_name=='tiegcm':
        envvar = 'TIEGCM_ROOT'
      else:
        envvar = 'TIMEGCM_ROOT'
      envvalue = getenv(envvar)
      if envvalue != '':
        rootdir = envvalue   
#       print 'Got env var ',envvar,' = ',rootdir,' (needed for ',job.model_name,' model_root)'
#
# model_root not on command-line, and env var not set, so prompt user:
#
      else:
        print 'Environment variable ',envvar,' is not set.'
        answer = raw_input('Run '+run.name+': Enter model_root directory for '+job.model_name+' (q to quit): ')
        if answer == 'q':
          sys.exit()
        if answer == '':
          get_options(arg,run,job,option) # try again
        rootdir = answer
        print 'Got response to prompt: rootdir=',rootdir
#
# Check for directory existence of either env var value, or prompt response:
#
      if not os.path.isdir(rootdir):
        print '>>> Cannot find model root directory ',rootdir
        get_options(arg,run,job,option) # try again
      
      scripts_dir = rootdir+'/scripts'
      if not os.path.isdir(scripts_dir):
        print '>>> Model root ',rootdir,' exists, but scripts directory ',scripts_dir,' does not exist.'
        get_options(arg,run,job,option) # try again
      return rootdir
#
# Get desired model resolution (5.0 or 2.5):
#
  elif arg == 'model_res':
    if option:
      if option != '5.0' and option != '2.5':
        print ">>> Unknown model resolution on the command line: '",option,"' (must be 5.0 or 2.5)"
        sys.exit()
      run.model_res = option
    else:
      answer = ''
      while answer != '5.0' and answer != '2.5':
        answer = raw_input("Run "+run.name+": Enter model resolution (5.0 or 2.5, default=5.0): ")
        if answer == 'q':
          sys.exit()
        elif answer == '':
          model_res = '5.0' 
          answer = model_res
        elif answer != '5.0' and answer != '2.5':
          print ">>> Unknown model resolution (must be 5.0 or 2.5)"
        run.model_res = answer
    return run.model_res
#
# Machine/platform:
#
  elif arg == 'machine':
    if option:
      if option != 'ys' and option != 'linux':
        print '>>> Unrecognized machine type found on command line: ',option
        print "    Machine must be either 'ys' (yellowstone) or 'linux'"
        sys.exit()
      job.machine = option
    else:
      job.machine = ''
      for line in os.popen('uname -a'): uname = line
      loc = uname.find('yslogin')
      if loc >= 0:
        job.machine = 'ys'
      else:
        loc = uname.find('Linux')
        if loc >= 0:
          job.machine = 'linux'
      if job.machine == '':
        print ">>> Could not determine machine (must be either 'linux' or 'ys')"
        sys.exit()
    return job.machine
#
# Compiler ('linux' platform only):
# ToDo: if compiler option is specified and execdir exists, force gmake clean in execdir.
#
  elif arg == 'compiler':
    if option:
      if job.machine != 'linux':
        print ">>> The '-compiler' option is valid only for the linux desktop platform"
        print "    Machine = ",job.machine
	sys.exit()
      if option != 'intel' and option != 'pgi' and option != 'gfort':
        print ">>> Invalid value for 'compiler' option: ",option
        print ">>> Value for 'compiler' option must be 'intel', 'pgi', or 'gfort'"
	sys.exit()
      return option
    else:                     # default is intel
      return 'intel'
#
# Get execution directory:
#
  elif arg == 'execdir':
    if option: # job script will validate
      return option
    else:
      tgcmtemp = getenv('TGCMTEMP') # User's large scratch temp directory
      if tgcmtemp != '':            # env var is set
        execdir = tgcmtemp+'/'+job.model_name+'_'+'res'+job.model_res+'_'+run.name+'/run'
        print 'Set ',execdir,' as execdir for ',run.name
        return execdir
      else: # prompt for execdir
        print 'Env var ',TGCMTEMP,' is not set.'
        answer = raw_input('Enter execution directory (execdir) for run '+run.name+': ')
        if not os.path.isdir:
          os.path.makedirs(answer)
          print "Made execdir directory ",answer," for run ",run.name
        else:
          print 'Set ',answer,' as execdir for ',run.name,' Note: the job script will make this directory for the run.'
        return answer
#
# If execdir exists, prompt for whether to empty the directory prior to the run.
#
#   if os.path.isdir(job.execdir):
#     print "Use existing execdir ",job.execdir,"? (yes/no)" 
#     answer = raw_input("  If no, all files in the directory will be removed prior to running the run) (default=yes): ")
#     if answer == 'q':
#       sys.exit()
#     elif answer == 'no':
#       for file in os.listdir(job.execdir):
#         os.remove(job.execdir+'/'+file)
#       print "Removed files in execdir ",job.execdir
#     else:
#       print "Files NOT removed from execdir ",job.execdir
#   else:
#     os.makedirs(job.execdir) # this makes parents and leaf as necessary
#     print "Made directory ",job.execdir," for run ",run_name
#
# Execution flag TRUE/FALSE (command-line only):
#
  elif arg == 'exec':
    if option:
      if option != 'yes' and option != 'no':
        print ">>> Bad value for ",arg,": '",option,"' (must be 'yes' or 'no' )"
        sys.exit()
      return option
    else:
      return 'yes'
#
# Get env var TGCMDATA:
#
  elif arg == 'tgcmdata':
    if option:
      if not os.path.isdir(option):
        print '>>> Could not find tgcmdata directory found on command line: ',option
        sys.exit()
      tgcmdata = option
    else:
      default_tgcmdata = ''
      if job.machine == 'ys':
        default_tgcmdata = '/glade/p/hao/tgcm/data'
      tgcmdata = getenv('TGCMDATA',default=default_tgcmdata)
      if tgcmdata:
        if not os.path.isdir(tgcmdata):
          print '>>> Could not find TGCMDATA directory ',tgcmdata
          print '    Am setting tgcmdata to execdir ',job.execdir
          tgcmdata = job.execdir
      else:
        print 'Env var TGCMDATA is not set.'
        answer = raw_input('Enter directory containing data for model '+job.model_name+' (default='+execdir+' or q to quit): ')
        if answer == 'q':
          sys.exit()
        if answer == '':
          print "NOTE: Setting tgcmdata to execdir (",tgcmdata,")"
          tgcmdata = job.execdir
        else:
          get_options(arg,run,job,option)
    return tgcmdata
#
# Number of processors (mpi tasks) (command-line only):
#
  elif arg == 'nprocs':
    if option:
      try: 
        job.nprocs = int(option)
      except:
        print '\n>>> Ooops, bad value for -nprocs option: ',option,' (must be an integer)' 
        sys.exit()
      return job.nprocs
#
# Project number (ys only) (command-line only):
#
  elif arg == 'project':
    if option:
      job.project = option
      return job.project
#
# LSF queue, e.g.: #BSUB -q regular (ys only) (command-line only):
#
  elif arg == 'queue':
    if option:
      job.queue = option
      return job.queue
#
# Wallclock limit (ys only):
#
  elif arg == 'wc':
    if option:
      job.wallclock = option
      return job.wallclock
#   else:
#     return job.calc_wallclock
#
# Execute flag (TRUE/FALSE):
#
  elif arg == 'execute':
    if option:
      job.execute = option
      return job.execute
#
# Submit flag (yes/no):
#
  elif arg == 'submit':
    if option:
      if option != 'yes' and option != 'no':
        print "\n>>> Oooops, bad value for -submit: ",option," (must be 'yes' or 'no')"
        sys.exit()
      job.submitflag = option
      return job.submitflag
#
# Stdout directory:
#
  elif arg == 'stdout_dir':
    if option:
      if not os.path.isdir(option):
        os.makedirs(option) # this makes parents and leaf as necessary
        print 'Made stdout directory ',option
      else: # do not clean stdout dir 
        print 'Will use existing stdout directory ',option
      job.stdout_dir = option
      return job.stdout_dir
#
# Directory in which to save history files:
# (absolute path or relative to execdir)
#
  elif arg == 'hist_dir':
    if option:
      histdir = option
      if option[0]=='/':  # absolute path
        histdir = option
      else:               # relative to execdir
        histdir = job.execdir+'/'+option

      if not os.path.isdir(histdir):
        os.makedirs(histdir) # this makes parents and leaf as necessary
        print 'Made history directory ',histdir
      else: # do not clean hist dir 
        print 'Will use existing history directory ',histdir

      if option[0] != '/':    # relative to execdir
        histdir = option      # this goes in the namelist read file

      histdir = histdir + '/' # file name will be appended
    else:
      histdir = './'
    return histdir
#
# Model timestep:
#
  elif arg == 'step':
    if option:
      return option
    else:
      if job.model_res == '5.0': step = '60'
      if job.model_res == '2.5': step = '30'
      return step
#
# Source history file:
#
# elif arg == 'source':
#   if option:
#     print 'source file from command-line: ',option
#     if not os.path.isfile(option):
#       print '>>> Cannot find source history file "',option,'" (this was read from the command-line)'
#       sys.exit()
#     else:
#       print 'Source history file (SOURCE) from command line: ',option
#       print 'NOTE: You may need to change SOURCE_START in the namelist input file(s).'
#       job.source = option
#       return job.source
#
# Unknown option:
#
  else:
    print '>>> Unknown arg: ',arg
    return None

#-----------------------------------------------------------------------
