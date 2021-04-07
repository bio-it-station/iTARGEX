import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t')
output_file = sys.argv[1].replace(".tsv", ".csv")
df.to_csv(output_file, index=False)
