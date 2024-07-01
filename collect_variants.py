import os
import argparse
import pandas as pd

def main():
    # parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_directory", 
        type=str, 
        help="Pass path to a directory containing parsed pcgr outputs"
        )
    parser.add_argument(
        "output_file", 
        type=str, 
        help="Pass path to a output xlsx file"
        )
    args = parser.parse_args()

    # collecting all xlsx files in directory and subdirectories
    xlsxfiles = [os.path.join(root, name)
             for root, dirs, files in os.walk(args.input_directory)
             for name in files
             if name.endswith((".xlsx"))]
    
    # Combining dataframes of all samples into one
    all_vars = pd.DataFrame()
    for file in xlsxfiles:
        df = pd.read_excel(file)
        sample = file.split("/")[-1].replace('.xlsx', '')
        df.insert(0, 'sample', sample)
        all_vars = pd.concat([all_vars, df], ignore_index=True)

    all_vars['TIER'] = pd.Categorical(all_vars['TIER'], ["TIER 1", "TIER 2", "TIER 3", "TIER 4t", "NONCODING"])
    all_vars = all_vars.sort_values('TIER') # Sort dataframe by TIER

    all_vars.to_excel(args.output_file, index=False) 

    

if __name__ == '__main__':
    main()