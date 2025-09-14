#!/bin/bash

# Make sure required Python packages are installed
pip install pandas openpyxl

# Create and run a Python script to combine all Excel files
python3 - << 'EOF'
import pandas as pd
import os
import glob

# Output file name
master_file = "master_gpcrdb_mor_lig_contacts.xlsx"

# Initialize an empty list to store all dataframes
all_dfs = []

# Get all Excel files in the current directory
excel_files = glob.glob("Interaction_data_*.xlsx")

print(f"Found {len(excel_files)} Excel files to process")

# Process each file
for file in excel_files:
    try:
        print(f"Reading {file}...")
        # Read the Excel file
        df = pd.read_excel(file)
        
        # Check if PDB and organism columns exist, if not, extract them from filename
        if 'PDB' not in df.columns or 'organism' not in df.columns:
            # Extract PDB and organism from filename
            parts = os.path.basename(file).split('_')
            if len(parts) >= 3:
                pdb_code = parts[2][:4]  # Extract the 4-character PDB code
                organism = "human" if "human" in file else "mouse"
                
                # Add the columns if they don't exist
                if 'PDB' not in df.columns:
                    df['PDB'] = pdb_code
                if 'organism' not in df.columns:
                    df['organism'] = organism
        
        # Add the dataframe to our list
        all_dfs.append(df)
        print(f"  ✓ Added {len(df)} rows from {file}")
        
    except Exception as e:
        print(f"  ✗ Error processing {file}: {e}")

# Combine all dataframes
if all_dfs:
    print("Combining all data...")
    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    # Save the combined data to the master file
    print(f"Saving to {master_file}...")
    combined_df.to_excel(master_file, index=False)
    
    print(f"✓ Master file created successfully with {len(combined_df)} total rows")
else:
    print("✗ No data to combine")

EOF