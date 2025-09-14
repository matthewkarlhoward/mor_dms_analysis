#!/bin/bash

# Make sure required Python packages are installed
pip install pandas openpyxl

# Create and run a Python script to process all Excel files
python3 - << 'EOF'
import pandas as pd
import os
import re

# Get all Excel files in the current directory
excel_files = [f for f in os.listdir('.') if f.endswith('.xlsx')]

# Process each file
for file in excel_files:
    # Extract PDB and organism from filename
    match = re.search(r'Interaction_data_([A-Z0-9]{4})_([a-z]+)\.xlsx', file)
    if match:
        pdb_code = match.group(1)
        organism = match.group(2)
        
        print(f"Processing {file} - PDB: {pdb_code}, Organism: {organism}")
        
        # Read the Excel file
        try:
            df = pd.read_excel(file)
            
            # Add the new columns
            df['PDB'] = pdb_code
            df['organism'] = organism
            
            # Save the modified file
            df.to_excel(file, index=False)
            print(f"  ✓ Successfully updated {file}")
            
        except Exception as e:
            print(f"  ✗ Error processing {file}: {e}")
    else:
        print(f"  ✗ Skipping {file} - doesn't match expected pattern")

print("Processing complete!")
EOF