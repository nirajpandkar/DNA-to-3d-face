import glob
import io
import os

import pandas as pd

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def create_id(row):
    return str(row.CHROM) + "_" + str(row.POS)

def generate(path):
    # Read in the raw vcf into pandas
    df = read_vcf(path)
    
    # Drop unnecessary columns
    df = df.drop(["ID", "FILTER", "INFO", "FORMAT", "QUAL"], axis=1)
    
    # Create a unique ID
    df["ID"] = df.apply(create_id, axis=1)
    
    # Reposition ID column
    ID = df.pop("ID")
    df.insert(0, ID.name, ID)
    df.head()
    
    # Drop unnecessary columns
    df = df.drop(["CHROM", "POS"], axis=1)
    
    # Clean the sample calls
    for col in df.columns:
        if col in ['ID', "REF", "ALT"]:
            pass
        else:
            df[col] = df[col].apply(lambda x: x.split(":")[0])
    return df


if __name__ == '__main__':
    final_df = pd.DataFrame()

    vcf_input_directory = "raw_vcf_files"
    for filepath in glob.glob(os.path.join(vcf_input_directory, "*.gz")):
        print(f"Processing {filepath}")
        try:
            df = generate(filepath)
        except:
            print(f"Something went wrong with {filepath}")
        final_df = final_df.append(df) 

    final_df.to_csv("masterfile.csv", index=False)
