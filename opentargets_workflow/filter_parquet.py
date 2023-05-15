import pandas as pd
import pyarrow.parquet as pq
from sys import argv

snplist = argv[3]
with open(snplist, 'r') as f:
  lines = f.readlines()

lines = [l.split(':', 1)[1] for l in lines]

lines = [l.split('_', 1)[0] for l in lines]

parquet_file = argv[1]
pqfile = pq.read_table(parquet_file)

df = pqfile.to_pandas()

# Iterate through the lines in the text file and search for them in the parquet file

results = []

for line in lines:
  result = df[df['position'] == int(line.strip())]
  results.append(result)

output_file = argv[2] + '.tsv'

# Concatenate the resulting DataFrames and save to a new file
final_result = pd.concat(results)
final_result.to_csv(output_file, index=False, sep="\t")
