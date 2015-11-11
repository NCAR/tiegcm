import re,os
#-----------------------------------------------------------------------
def replace_string(file,oldstr,newstr):
#
# Replace all occurrences of "oldstr" in oldfile to "newstr" in newfile.
# Note oldfile and newfile may be the same, i.e., overwrite oldfile. 
#
# Read contents first, then close input file.
#
  file_read = open(file,'r')
  contents = file_read.read()
  file_read.close()
# print contents
#
# Open output file and change chars per line.
#
  fw = open(file,'w')
  for line in contents:
    newline = line.replace(oldstr,newstr)
    fw.write(newline)
  fw.close()
#-----------------------------------------------------------------------
def remove_comments(file):
#
# Strip '!' comments from input file.
#
# Read input file into list of lines:
#
  f = open(file,'r')
  lines = f.readlines()  
  lines = [line.rstrip() for line in lines]
  f.close()
#
# Reopen input file for writing, and loop over lines,
# removing comments.
#
  f = open(file,'w')
  nlines = 0
  for line in lines:
    nlines = nlines+1
    loc = 0
    loc = line.find('!')
    if loc == 0:           # If comment is in first column (loc==0), skip the line
      continue
    elif loc > 0:          # Remove line from comment char to end
      line = line[:loc-1]
    f.write(line+'\n')
  f.close()
#-----------------------------------------------------------------------
def getenv(var,default=''):
#
# Return setting of environment variable var. If it is not set and
# there is a default, set the env var to the default and return that
# value. If the env var is not set and there is no default, return ''.
#
  try:
    value = os.environ[var]
    return value
  except KeyError:
    if default:
      os.environ[var] = default
      print 'Note: Set env var ',var,' to default: ',default
      return os.environ[var]
    else:
      return ''
#-----------------------------------------------------------------------
