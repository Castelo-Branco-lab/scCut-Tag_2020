#!/usr/bin/env python3

#import pandas as pd
import sys
import pandas as pd

def convert_to_integer(l):
  try:
    return([int(x) for x in l])
  except ValueError:
    return(False)

# Parse arguments

def main():
  # Error return:
  err_msg = '*** Usage is merge_bedpe.py file1 column file2 column file3 column ... ***\n'
  
  to_merge = sys.argv[1:]
  to_merge =  [(to_merge[x],to_merge[x+1]) for x in range(0,len(to_merge),2)]
  
  if len(to_merge) %2 != 0:
    sys.stderr.write(err_msg)
    sys.exit(1)

  
  columns = convert_to_integer([x[1] for x in to_merge])
  if not columns:
    sys.stderr.write("Column numbers can't be converted to integers\n" + err_msg)
    sys.exit(1)
    
  files   = [x[0] for x in to_merge]
  
  
  dfs = []
  for x,f in enumerate(files):
    df                     = pd.read_table(f, header=None)
    colnames               = [x for x in df.columns]
    colnames[columns[x]-1] = 'merge_column'
    df.columns             = colnames
    dfs.append(df)






main()