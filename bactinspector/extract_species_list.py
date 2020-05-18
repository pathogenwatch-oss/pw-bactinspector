import sys
import pandas as pd

species_df = pd.read_parquet('data/all_complete_refseq.k21s1000.species.pqt')
print(species_df.head(2))
