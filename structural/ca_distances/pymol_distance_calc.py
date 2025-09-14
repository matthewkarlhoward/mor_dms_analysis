import csv
from pymol import cmd

def calculate_ca_distances(structure1, structure2, output_file="ca_distances.csv"):
    """
    Calculate distances between alpha carbons of corresponding residues
    in two protein structures and save to CSV.
    
    Parameters:
    structure1 (str): Name of first structure in PyMOL
    structure2 (str): Name of second structure in PyMOL  
    output_file (str): Output CSV filename
    """
    
    # Get alpha carbon coordinates for both structures
    ca_coords1 = {}
    ca_coords2 = {}
    
    # Extract CA coordinates from structure 1
    cmd.iterate_state(1, f"{structure1} and name CA", 
                     "ca_coords1[resi] = (x, y, z)", 
                     space=locals())
    
    # Extract CA coordinates from structure 2
    cmd.iterate_state(1, f"{structure2} and name CA", 
                     "ca_coords2[resi] = (x, y, z)", 
                     space=locals())
    
    # Find common residues and convert keys to integers
    resi_keys1 = {int(k): v for k, v in ca_coords1.items()}
    resi_keys2 = {int(k): v for k, v in ca_coords2.items()}
    common_residues = set(resi_keys1.keys()) & set(resi_keys2.keys())
    
    # Calculate distances
    distances = []
    for resi in sorted(common_residues):
        coord1 = resi_keys1[resi]
        coord2 = resi_keys2[resi]
        
        # Calculate Euclidean distance
        distance = ((coord1[0] - coord2[0])**2 + 
                   (coord1[1] - coord2[1])**2 + 
                   (coord1[2] - coord2[2])**2)**0.5
        
        distances.append([resi, distance])
    
    # Write to CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Residue_Number', 'CA_Distance'])
        writer.writerows(distances)
    
    print(f"Distances calculated for {len(distances)} residues")
    print(f"Results saved to {output_file}")
    
    return distances

# Example usage for your structures (comment out to avoid auto-execution):
# distances = calculate_ca_distances("5c1m", "4dkl", "5c1m_4dkl_ca_distances.csv")

# Alternative one-liner approach if you prefer:
def quick_ca_distances(struct1, struct2, output="distances.csv"):
    """Simplified version for quick calculations"""
    distances = []
    
    # Get residue numbers that exist in both structures
    residues1 = set()
    residues2 = set()
    
    cmd.iterate(f"{struct1} and name CA", "residues1.add(resi)", space=locals())
    cmd.iterate(f"{struct2} and name CA", "residues2.add(resi)", space=locals())
    
    # Convert to integers
    residues1 = {int(r) for r in residues1}
    residues2 = {int(r) for r in residues2}
    
    common_residues = residues1 & residues2
    
    for resi in sorted(common_residues):
        # Use PyMOL's built-in distance calculation
        dist_obj = f"dist_{resi}"
        cmd.distance(dist_obj, 
                    f"{struct1} and resi {resi} and name CA",
                    f"{struct2} and resi {resi} and name CA")
        
        # Get the distance value
        dist_value = cmd.get_distance(f"{struct1} and resi {resi} and name CA",
                                     f"{struct2} and resi {resi} and name CA")
        
        distances.append([resi, dist_value])
        
        # Clean up the distance object
        cmd.delete(dist_obj)
    
    # Write to CSV
    with open(output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Residue_Number', 'CA_Distance'])
        writer.writerows(distances)
    
    print(f"Quick calculation complete: {len(distances)} residues written to {output}")